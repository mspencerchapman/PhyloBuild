# Getting the BAMs. Two options, set by `bam_staging` in the config:
#   none    use the BAM paths in the sample sheet as they are
#   sanger  stage each BAM out of iRODS with stageBam.pl first (needs a `project` column in the sample sheet;
#           the in silico normal uses reference.normal_project)
#           Samples in a test project (e.g. nst_links/test/9999) need irods_db = test in the sample sheet; the
#           in silico normal uses reference.normal_irods_db (default live).
# Every later step asks bam_path() for a sample's BAM, so they work the same either way.

STAGING = config.get("bam_staging", "none")
STAGE_ROOT = config.get("staging_dir", f"{OUT}/staged_bams")
NORMAL = config["reference"]["normal_name"]

if STAGING not in ("none", "sanger"):
    raise ValueError(f"bam_staging must be 'none' or 'sanger', not '{STAGING}'")

#stageBam.pl writes its sample list to a fixed file name in its output directory, so concurrent calls sharing one
#directory overwrite each other's lists. Each sample therefore gets its own stageBam.pl output directory.
def staged_bam(project, sample):
    return f"{STAGE_ROOT}/{sample}/{project}/{sample}/mapped_sample/{sample}.sample.dupmarked.bam"

def bam_path(donor, sample):
    if STAGING == "sanger":
        project = sample_row(donor, sample).get("project")
        if not project:
            raise ValueError(f"bam_staging is 'sanger' but sample {sample} has no project in the sample sheet")
        return staged_bam(project, sample)
    return sample_row(donor, sample)["bam"]

def normal_bam_path():
    if STAGING == "sanger":
        return staged_bam(config["reference"]["normal_project"], NORMAL)
    return config["reference"]["normal_bam"]

#The BAM with its index (.bai) and stats (.bas) files. Rules list all the files they read, so that with
#keep_staged_bams: false Snakemake deletes a staged BAM only once the last job that reads it has finished.
def bam_files(donor, sample, exts=(".bai", ".bas")):
    bam = bam_path(donor, sample)
    return [bam] + [bam + e for e in exts]

def normal_bam_files():
    bam = normal_bam_path()
    return [bam, bam + ".bai", bam + ".bas"]

def irods_db(project, sample):
    """iRODS database (stageBam.pl -l) holding a sample's BAM: the sample sheet's irods_db column, default live."""
    if sample == NORMAL and str(project) == str(config["reference"].get("normal_project")):
        db = config["reference"].get("normal_irods_db", "live")
    else:
        rows = [r for donor_rows in SAMPLES.values() for r in donor_rows.values()
                if r["sample"] == sample and str(r.get("project")) == str(project)]
        db = (rows[0].get("irods_db") or "live") if rows else "live"
    if db not in ("live", "test"):
        raise ValueError(f"irods_db must be live or test, not '{db}' ({project}/{sample})")
    return db

def maybe_temp(path):
    return path if config.get("keep_staged_bams", True) else temp(path)

rule stage_bam:
    """Stage one sample's mapped BAM (with .bai/.bas) from iRODS.
    stageBam.pl only submits its own LSF transfer job and exits, so this waits for that job to end."""
    output:
        bam = maybe_temp(f"{STAGE_ROOT}/{{sample}}/{{project}}/{{sample}}/mapped_sample/{{sample}}.sample.dupmarked.bam"),
        bai = maybe_temp(f"{STAGE_ROOT}/{{sample}}/{{project}}/{{sample}}/mapped_sample/{{sample}}.sample.dupmarked.bam.bai"),
        bas = maybe_temp(f"{STAGE_ROOT}/{{sample}}/{{project}}/{{sample}}/mapped_sample/{{sample}}.sample.dupmarked.bam.bas"),
    params:
        root = f"{STAGE_ROOT}/{{sample}}",
        sample_dir = f"{STAGE_ROOT}/{{sample}}/{{project}}/{{sample}}/mapped_sample",
        lsf_logs = f"{STAGE_ROOT}/logs/lsf",
        db = lambda w: irods_db(w.project, w.sample),
    wildcard_constraints:
        project = r"\d+",
        sample = r"[^/]+",
    envmodules: *modules("staging")
    resources:
        mem_mb = 1000,
        runtime = 720,
        irods_transfers = 1,     # capped in total by the profile's `resources: irods_transfers`
    retries: 2
    log: f"{STAGE_ROOT}/logs/{{project}}_{{sample}}.log"
    shell:
        r"""
        # Start clean: stageBam.pl keeps its own record of what it has downloaded (tmpExportData), so after a failed
        # attempt (when Snakemake deletes the outputs) it would otherwise skip the files as "already downloaded"
        rm -rf {params.root}
        mkdir -p {params.lsf_logs}
        stageBam.pl -p {wildcards.project} -s {wildcards.sample} -t m -l {params.db} -o {params.root} -fo \
            -G {resources.lsf_group} -lo {params.lsf_logs} > {log} 2>&1
        # Wait for the transfer job(s) stageBam.pl submitted, polling bjobs (bwait reports these array jobs as ended
        # while they are still running). The rule's runtime limit stops this if a transfer hangs.
        # exportData.pl sometimes hangs, idle, once it has started its last file (the .met.gz, which comes after the
        # BAM, .bai and .bas). So if its log has been quiet for 15 minutes after "(N of N)" and the BAM is complete
        # (ends with the BGZF end-of-file block), the transfer job is killed and the files are used.
        bgzf_eof=1f8b08040000000000ff0600424302001b0003000000000000000000
        for id in $(grep -oE 'Job <[0-9]+> is submitted' {log} | grep -oE '[0-9]+'); do
            echo "Waiting for transfer job $id" >> {log}
            stat=""
            until [ "$stat" = "DONE" ] || [ "$stat" = "EXIT" ]; do
                sleep 30
                stat=$(bjobs -noheader -o stat $id 2>/dev/null | tail -n1 || true)
                xlog=$(ls {params.lsf_logs}/export.$id.* 2>/dev/null | head -n1 || true)
                if [ "$stat" = "RUN" ] && [ -n "$xlog" ] && [ -n "$(find $xlog -mmin +15)" ] \
                    && grep -qE '\(([0-9]+) of \1\)$' $xlog; then
                    bam=$(find {params.root} -path '*/tmpExportData' -prune -o \
                        -name {wildcards.sample}.sample.dupmarked.bam -print | head -n1)
                    if [ -n "$bam" ] && [ -s $bam.bai ] && [ -s $bam.bas ] \
                        && [ "$(tail -c 28 $bam | od -An -tx1 | tr -d ' \n')" = "$bgzf_eof" ]; then
                        echo "Transfer job $id idle after starting its last file, and the BAM is complete: killing it" >> {log}
                        bkill $id >> {log} 2>&1 || true
                        stat=DONE
                    fi
                fi
            done
            echo "Transfer job $id: $stat" >> {log}
            if [ "$stat" != "DONE" ]; then exit 1; fi
        done
        # The test database writes to <output>/mapped_sample/ rather than <output>/<project>/<sample>/mapped_sample/
        if [ ! -e {output.bam} ] && [ -d {params.root}/mapped_sample ]; then
            mkdir -p {params.sample_dir}
            mv {params.root}/mapped_sample/* {params.sample_dir}/
        fi
        if [ ! -s {output.bam} ] && ls {params.sample_dir}/*.cram > /dev/null 2>&1; then
            echo "{wildcards.sample} was staged as a CRAM, not a BAM. CRAM input hasn't been tested with the hairpin filter yet." >> {log}
            exit 1
        fi
        for f in {output.bam} {output.bai} {output.bas}; do
            if [ ! -s $f ]; then echo "Missing after staging: $f" >> {log}; exit 1; fi
        done
        """
