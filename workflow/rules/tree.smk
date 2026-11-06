# Tree building: one job per donor x upfront colony filter x filtering type.
# Replaces filtering_from_table_mix_remove_GRCh38.R.

TREE_DIR = f"{OUT}/{{donor}}/trees/{{donor}}_{{cf}}"

rule build_tree:
    input:
        mats = f"{OUT}/{{donor}}/mats_and_params/mats_and_params_{{donor}}_{{cf}}_postMS_reduced.rds",
        sensitivity = f"{OUT}/{{donor}}/qc/sensitivity_{{donor}}.txt",
        depth_cutoffs = f"{OUT}/{{donor}}/qc/depth_cutoffs.tsv",
    output:
        tree = f"{TREE_DIR}/tree_files/tree_{{donor}}_{{cf}}_{{filtering}}_post_mix_post_dup.tree",
        muts = f"{TREE_DIR}/annotated_muts/annotated_mut_set_{{donor}}_{{cf}}_{{filtering}}_post_mix_post_dup.rds",
        plots = f"{TREE_DIR}/plots/{{donor}}_{{cf}}_{{filtering}}_final_summary_trees.pdf",
    params:
        outdir = TREE_DIR,
        run_id = lambda w: f"{w.donor}_{w.cf}",
        filtering = lambda w: w.filtering,
        tree = lambda w: donor_cfg(w.donor)["tree"],
        group_pattern = lambda w: donor_cfg(w.donor)["tree"]["group_pattern"] or w.donor,
        age = lambda w: donor_cfg(w.donor)["age"],
        legacy_rdata = lambda w: donor_cfg(w.donor)["legacy_rdata"],
        build = config["genome_build"],
        mpboot = config["tools"]["mpboot"],
        vagrent = config["tools"]["vagrent"],
    envmodules: *modules("tree")
    resources:
        mem_mb = lambda w, attempt: mem_for_samples(w.donor, 16000, 200, attempt),
        runtime = lambda w, attempt: 2880 * attempt,
    log: f"{OUT}/{{donor}}/logs/build_tree_{{cf}}_{{filtering}}.log"
    script: "../scripts/build_tree.R"
