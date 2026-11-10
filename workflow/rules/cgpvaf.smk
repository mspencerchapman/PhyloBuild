# Read counts at every site in every sample, with cgpVAF.
# createVafCmd.pl only symlinked each sample's files into one folder and wrote out cgpVaf.pl commands; this does
# the same directly from the sample sheet, so there is no interactive step, no database lookup, no split config
# file and no editing of run_bsub.sh. Failed chunks are retried automatically with more memory.

VCF_EXTENSION = {"snv": ".caveman_c.annot.vcf.gz", "indel": ".pindel.annot.vcf.gz"}
CGPVAF_TYPE = {"snv": "snp", "indel": "indel"}

def cgpvaf_dir(w):
    return f"{OUT}/{w.donor}/cgpvaf/{w.vartype}"

rule cgpvaf_stage:
    """Folder of symlinks in the layout cgpVaf.pl expects: <sample>.bam(.bai/.bas) and <sample><vcf extension>."""
    input:
        bams = lambda w: [bam_path(w.donor, s) for s in samples_of(w.donor)],
        vcfs = lambda w: [sample_row(w.donor, s)["caveman_vcf" if w.vartype == "snv" else "pindel_vcf"] for s in samples_of(w.donor)],
        normal = normal_bam_path(),
        index_files = lambda w: [f for s in samples_of(w.donor) for f in bam_files(w.donor, s)[1:]] + normal_bam_files()[1:],
    output:
        f"{OUT}/{{donor}}/cgpvaf/{{vartype}}/input/staged.txt",
    run:
        import os
        stage = os.path.dirname(output[0])
        ext = VCF_EXTENSION[wildcards.vartype]
        def link(src, dest):
            if not os.path.exists(src):
                raise FileNotFoundError(f"{src} not found (needed by cgpVAF)")
            dest = os.path.join(stage, dest)
            if os.path.lexists(dest):
                os.remove(dest)
            #Link to the real file, not through other links: the cgpVAF container only mounts /nfs and /lustre, so
            #a link via /software/SGP/nst_links (BAMs used in place) would be broken inside it
            os.symlink(os.path.realpath(src), dest)
        entries = list(zip(samples_of(wildcards.donor), input.bams))
        entries.append((NORMAL, input.normal))
        for name, bam in entries:
            link(bam, f"{name}.bam")
            link(bam + ".bai", f"{name}.bam.bai")
            link(bam + ".bas", f"{name}.bam.bas")
        for s, vcf in zip(samples_of(wildcards.donor), input.vcfs):
            link(vcf, s + ext)
            if os.path.exists(vcf + ".tbi"):
                link(vcf + ".tbi", s + ext + ".tbi")
        with open(output[0], "w") as fh:
            fh.write("\n".join(n for n, _ in entries) + "\n")

#Files cgpVaf.pl reads through the symlinks in input/: the batch's BAMs and the normal's. Listed as inputs so that
#temporary staged BAMs are kept until the last cgpVAF job of their batch has finished.
def cgpvaf_batch_bams(w):
    return [f for s in cgpvaf_batches(w.donor)[int(w.batch)] for f in bam_files(w.donor, s)] + normal_bam_files()

def cgpvaf_args(w):
    return " ".join([
        "-d", f"{cgpvaf_dir(w)}/input",
        "-o", f"{cgpvaf_dir(w)}/batch{w.batch}",
        "-a", CGPVAF_TYPE[w.vartype],
        "-e", VCF_EXTENSION[w.vartype],
        "-g", config["reference"]["genome"],
        "-hdr", config["reference"]["high_depth_regions"],
        "-b", f"{OUT}/{w.donor}/calls/{w.donor}_{w.vartype}_sites.bed",
        "-bo", "1",
        "-mq", str(donor_cfg(w.donor)["cgpvaf_mapq"]),
        "-nn", config["reference"]["normal_name"],
        "-tn", " ".join(cgpvaf_batches(w.donor)[int(w.batch)]),
    ])

rule cgpvaf_chromosome:
    input:
        staged = f"{OUT}/{{donor}}/cgpvaf/{{vartype}}/input/staged.txt",
        bed = f"{OUT}/{{donor}}/calls/{{donor}}_{{vartype}}_sites.bed",
        bams = cgpvaf_batch_bams,
    output:
        touch(f"{OUT}/{{donor}}/cgpvaf/{{vartype}}/batch{{batch}}/done/{{chrom}}.done"),
    params:
        args = cgpvaf_args,
        cgpvaf = config["tools"]["cgpvaf"],
    envmodules: *modules("cgpvaf")
    resources:
        mem_mb = lambda w, attempt: 4000 * attempt,
        runtime = lambda w, attempt: 720 * attempt,
        lsf_queue = lambda w, attempt: "normal" if attempt == 1 else "long",  #a retry gets more time than normal allows
    log: f"{OUT}/{{donor}}/logs/cgpvaf/{{vartype}}_batch{{batch}}_{{chrom}}.log"
    shell:
        "{params.cgpvaf} {params.args} -chr {wildcards.chrom} > {log} 2>&1"

rule cgpvaf_concat:
    input:
        done = lambda w: expand(f"{OUT}/{{donor}}/cgpvaf/{{vartype}}/batch{{batch}}/done/{{chrom}}.done", donor=w.donor, vartype=w.vartype, batch=w.batch, chrom=CHROMOSOMES),
        bams = cgpvaf_batch_bams,
    output:
        touch(f"{OUT}/{{donor}}/cgpvaf/{{vartype}}/batch{{batch}}/done/concat.done"),
    params:
        args = cgpvaf_args,
        cgpvaf = config["tools"]["cgpvaf"],
    envmodules: *modules("cgpvaf")
    resources:
        mem_mb = lambda w, attempt: 1000 * attempt,
    log: f"{OUT}/{{donor}}/logs/cgpvaf/{{vartype}}_batch{{batch}}_concat.log"
    shell:
        "{params.cgpvaf} {params.args} -ct 1 > {log} 2>&1"

rule cgpvaf_merge:
    """One matrix of MTR/DEP per sample for the donor, found by column name (replaces the cut/paste merge)."""
    input:
        lambda w: expand(f"{OUT}/{{donor}}/cgpvaf/{{vartype}}/batch{{batch}}/done/concat.done", donor=w.donor, vartype=w.vartype, batch=cgpvaf_batches(w.donor)),
    output:
        f"{OUT}/{{donor}}/cgpvaf/{{donor}}_merged_{{vartype}}.tsv",
    params:
        batch_dirs = lambda w: [f"{OUT}/{w.donor}/cgpvaf/{w.vartype}/batch{b}" for b in cgpvaf_batches(w.donor)],
        expected_samples = lambda w: samples_of(w.donor),
    envmodules: *modules("r")
    resources:
        mem_mb = lambda w, attempt: mem_for_samples(w.donor, 16000, 60, attempt),
    log: f"{OUT}/{{donor}}/logs/cgpvaf/{{vartype}}_merge.log"
    script: "../scripts/merge_cgpvaf.R"
