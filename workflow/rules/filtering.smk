# Sex check, filter parameters (with/without upfront colony filtering), hairpin reduction, depth window and
# sensitivity. Replaces Mutation_filtering_get_parameters_v2.R, Reducing_mutset_from_MSfilters.R and
# Sensitivity_analysis_from_SNPs.R.

MERGED = {
    "snv": f"{OUT}/{{donor}}/cgpvaf/{{donor}}_merged_snv.tsv",
    "indel": f"{OUT}/{{donor}}/cgpvaf/{{donor}}_merged_indel.tsv",
}
MATS = f"{OUT}/{{donor}}/mats_and_params"

rule sex_check:
    input:
        snv = MERGED["snv"],
        indel = MERGED["indel"],
    output:
        summary = f"{OUT}/{{donor}}/qc/sex_check.tsv",
        per_colony = f"{OUT}/{{donor}}/qc/sex_check_per_colony.tsv",
        plot = f"{OUT}/{{donor}}/qc/sex_check_per_colony.pdf",
        sex = f"{OUT}/{{donor}}/qc/sex.txt",
    params:
        normal = config["reference"]["normal_name"],
        build = config["genome_build"],
        sex = lambda w: donor_cfg(w.donor)["sex"],
        mode = lambda w: donor_cfg(w.donor)["sex_check"],
    envmodules: *modules("r")
    resources:
        mem_mb = 16000,
    log: f"{OUT}/{{donor}}/logs/sex_check.log"
    script: "../scripts/sex_check.R"

rule filter_params:
    input:
        snv = MERGED["snv"],
        indel = MERGED["indel"],
        sex = f"{OUT}/{{donor}}/qc/sex.txt",
    output:
        mats = f"{MATS}/mats_and_params_{{donor}}_{{cf}}.rds",
        vaf_plots = f"{MATS}/colony_vaf_plots_{{donor}}_{{cf}}.pdf",
        peak_vaf = f"{MATS}/colony_peak_vaf_{{donor}}_{{cf}}.tsv",
    params:
        normal = config["reference"]["normal_name"],
        colony_filter = lambda w: w.cf,
        exclude_samples = lambda w: donor_cfg(w.donor)["exclude_samples"],
        legacy_rdata = lambda w: donor_cfg(w.donor)["legacy_rdata"],
    envmodules: *modules("r")
    resources:
        mem_mb = lambda w, attempt: mem_for_samples(w.donor, 8000, 100, attempt),
        runtime = 720,
    log: f"{OUT}/{{donor}}/logs/filter_params_{{cf}}.log"
    script: "../scripts/filter_params.R"

rule reduce:
    """SNVs that passed the hairpin filter (all indels kept), then drop near-certain germline variants."""
    input:
        mats = f"{MATS}/mats_and_params_{{donor}}_{{cf}}.rds",
        hairpin_sites = f"{OUT}/{{donor}}/hairpin/{{donor}}_hairpin_pass_sites.bed",
    output:
        postMS = f"{MATS}/mats_and_params_{{donor}}_{{cf}}_postMS.rds",
        reduced = f"{MATS}/mats_and_params_{{donor}}_{{cf}}_postMS_reduced.rds",
    params:
        legacy_rdata = lambda w: donor_cfg(w.donor)["legacy_rdata"],
    envmodules: *modules("r")
    resources:
        mem_mb = lambda w, attempt: 16000 * attempt,
    log: f"{OUT}/{{donor}}/logs/reduce_{{cf}}.log"
    script: "../scripts/reduce.R"

rule depth_cutoffs:
    """Mean-depth window for the tree step, from the full (no colony filter) mutation set."""
    input:
        mats = f"{MATS}/mats_and_params_{{donor}}_none.rds",
    output:
        tsv = f"{OUT}/{{donor}}/qc/depth_cutoffs.tsv",
        plot = f"{OUT}/{{donor}}/qc/depth_cutoffs.pdf",
    params:
        setting = lambda w: donor_cfg(w.donor)["tree"]["depth_cutoffs"],
        k = lambda w: donor_cfg(w.donor)["tree"]["depth_mad_k"],
    envmodules: *modules("r")
    resources:
        mem_mb = 8000,
    log: f"{OUT}/{{donor}}/logs/depth_cutoffs.log"
    script: "../scripts/depth_cutoffs.R"

rule sensitivity:
    """Per-sample sensitivity: the proportion of germline SNVs/indels each sample's calls picked up."""
    input:
        mats = f"{MATS}/mats_and_params_{{donor}}_none_postMS.rds",
        snv_calls = lambda w: expand(f"{OUT}/{{donor}}/hairpin/{{sample}}.hairpin_pass.vcf", donor=w.donor, sample=samples_of(w.donor)),
        indel_calls = lambda w: expand(f"{OUT}/{{donor}}/calls/indel/{{sample}}.pass.vcf", donor=w.donor, sample=samples_of(w.donor)),
    output:
        f"{OUT}/{{donor}}/qc/sensitivity_{{donor}}.txt",
    params:
        samples = lambda w: samples_of(w.donor),
        depth_cutoff = lambda w: donor_cfg(w.donor)["sensitivity_min_total_depth"],
    envmodules: *modules("r")
    resources:
        mem_mb = 16000,
    log: f"{OUT}/{{donor}}/logs/sensitivity.log"
    script: "../scripts/sensitivity.R"
