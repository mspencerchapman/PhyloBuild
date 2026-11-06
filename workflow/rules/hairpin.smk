# LCM hairpin artefact filter (the Sanger LCM filtering container), one job per sample.
# Ported from runlcmfilter.sh: the same five steps, with the paths taken from the config.

rule hairpin_filter:
    input:
        vcf = lambda w: sample_row(w.donor, w.sample)["caveman_vcf"],
        bam = lambda w: bam_path(w.donor, w.sample),
        bai = lambda w: bam_path(w.donor, w.sample) + ".bai",
    output:
        passed = f"{OUT}/{{donor}}/hairpin/{{sample}}_passed.vcf",
        pre = temp(f"{OUT}/{{donor}}/hairpin/{{sample}}.pre.vcf"),
        annovar = temp(f"{OUT}/{{donor}}/hairpin/{{sample}}.pre.annovar.txt"),
        annot = temp(f"{OUT}/{{donor}}/hairpin/{{sample}}.pre.annovar.annot.txt"),
        full = temp(f"{OUT}/{{donor}}/hairpin/{{sample}}.pre.annovar.annot.full.txt"),
    params:
        outdir = f"{OUT}/{{donor}}/hairpin",
        sif = config["tools"]["hairpin_sif"],
        singularity = config["tools"]["singularity"],
        bind = config["tools"]["singularity_bind"],
        genome = config["reference"]["genome"],
        snps = config["reference"]["snp_database"],
        asmd = lambda w: donor_cfg(w.donor)["asmd_threshold"],
        frag = lambda w: donor_cfg(w.donor)["hairpin_fragment_threshold"],
    envmodules: *modules("hairpin")
    threads: 1
    resources:
        mem_mb = 16000,
        runtime = 720,
    log: f"{OUT}/{{donor}}/logs/hairpin/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        run="{params.singularity} run --bind {params.bind}"
        echo "Running filters with fragment threshold {params.frag} and ASMD threshold {params.asmd}" > {log}
        $run --app preselect {params.sif} -a {params.asmd} -v {input.vcf} > {output.pre} 2>> {log}
        $run --app imitateANNOVAR {params.sif} -v {output.pre} > {output.annovar} 2>> {log}
        $run --app annotateBAMStatistics {params.sif} -a {output.annovar} -b {input.bam} -t {threads} > {output.annot} 2>> {log}
        $run --app additionalBAMStatistics {params.sif} -a {output.annot} -b {input.bam} -t {threads} -r {params.genome} -s {params.snps} > {output.full} 2>> {log}
        $run --app filtering {params.sif} -a {output.full} -v {output.pre} -o {params.outdir} -p {wildcards.sample} -f {params.frag} >> {log} 2>&1
        """

rule hairpin_pass:
    """PASS lines from the hairpin output (replaces filter_for_bed_file.pl)."""
    input:
        f"{OUT}/{{donor}}/hairpin/{{sample}}_passed.vcf",
    output:
        f"{OUT}/{{donor}}/hairpin/{{sample}}.hairpin_pass.vcf",
    shell:
        r"""awk -F'\t' '!/^#/ && $7 == "PASS"' {input} > {output}"""

rule hairpin_sites:
    """SNV sites passing the hairpin filter in at least one sample."""
    input:
        lambda w: expand(f"{OUT}/{{donor}}/hairpin/{{sample}}.hairpin_pass.vcf", donor=w.donor, sample=samples_of(w.donor)),
    output:
        f"{OUT}/{{donor}}/hairpin/{{donor}}_hairpin_pass_sites.bed",
    resources:
        mem_mb = 4000,
    shell:
        "cat {input} | cut -f 1,2,4,5 | sort -u > {output}"
