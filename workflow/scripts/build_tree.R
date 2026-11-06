# Tree building for one donor / upfront colony filter / filtering type (rule build_tree)
log = file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
library(phylobuild)

p = snakemake@params
t = p$tree
x = read_mats_and_params(snakemake@input$mats)
cut = read.delim(snakemake@input$depth_cutoffs)
depth_cutoffs = setNames(cut$value, cut$name)[c("auto_low", "auto_high", "xy_low", "xy_high")]
empty_to_null = function(v) if (length(v) == 0) NULL else unlist(v)

annotate = NULL
if (isTRUE(p$vagrent$enabled)) {
  annotate = vagrent_annotator(work_prefix = file.path(p$outdir, "mutation_vcfs", paste0("mutations_all_", p$run_id, "_", p$filtering)),
                               cache = p$vagrent$cache, species = p$vagrent$species, assembly = p$vagrent$assembly, command = p$vagrent$command)
  dir.create(file.path(p$outdir, "mutation_vcfs"), recursive = TRUE, showWarnings = FALSE)
}

out = build_tree(COMB_mats = x$COMB_mats,
                 filter_params = x$filter_params,
                 run_id = p$run_id,
                 filtering_type = p$filtering,
                 output_dir = p$outdir,
                 sensitivity_path = snakemake@input$sensitivity,
                 donor_id = p$group_pattern,
                 recipient_id = t$recipient_pattern,
                 build = p$build,
                 min_sample_mean_cov = t$min_sample_mean_cov,
                 exclude_muts = empty_to_null(t$exclude_muts),
                 keep_muts = empty_to_null(t$keep_muts),
                 polytomous = t$polytomous,
                 bb_cutoff = t$bb_cutoff,
                 germline_log10_cutoff = t$germline_log10_cutoff,
                 depth_cutoffs = depth_cutoffs,
                 filter_y = t$filter_y,
                 keep_ancestral = t$keep_ancestral,
                 remove_nonclonal = t$remove_nonclonal,
                 germline_addback = t$germline_addback,
                 snvs_only = t$snvs_only,
                 CN_table = if (is.null(t$cn_table)) NULL else read.csv(t$cn_table, stringsAsFactors = FALSE),
                 genes_to_annotate = if (is.null(t$genes_to_annotate)) NULL else readLines(t$genes_to_annotate),
                 donor_age = if (is.null(p$age)) NA else p$age,
                 annotate = annotate,
                 mpboot = p$mpboot,
                 mpboot_bootstraps = t$mpboot_bootstraps,
                 seed = t$seed,
                 treefit_pval_cutoff = t$poor_fit_pval,
                 remove_poor_fit = isTRUE(t$remove_poor_fit),
                 addback_require_filters = isTRUE(t$addback_require_filters),
                 clonality_exclude_cn_xy = isTRUE(t$clonality_exclude_cn_xy),
                 duplicate_min_private = t$duplicate_min_private,
                 duplicate_prop = t$duplicate_prop,
                 legacy_rdata = isTRUE(p$legacy_rdata))
message("Final tree: ", out$tree)
