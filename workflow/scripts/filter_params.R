# Filter parameters for every mutation, with or without upfront colony filtering (rule filter_params)
log = file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
library(phylobuild)

p = snakemake@params
cm = import_cgpvaf_SNV_and_INDEL(snakemake@input$snv, snakemake@input$indel, ref_ID = p$normal)
sex = readLines(snakemake@input$sex)
cutoff = if (p$colony_filter == "none") NULL else as.numeric(sub("^cf", "", p$colony_filter))

res = compute_filter_params(cm, sex = sex, exclude_samples = unlist(p$exclude_samples), mixed_colony_vaf_cutoff = cutoff,
                            vaf_plot_file = snakemake@output$vaf_plots)
save_object(list(COMB_mats = res$COMB_mats, filter_params = res$filter_params), snakemake@output$mats,
            legacy_name = if (isTRUE(p$legacy_rdata)) TRUE)

if (is.null(cutoff)) {
  pdf(snakemake@output$vaf_plots); plot.new(); text(0.5, 0.5, "No upfront colony filtering for this run"); dev.off()
  peak = data.frame(sample = colnames(res$COMB_mats$NV), peak_vaf = NA, removed = FALSE)
} else {
  peak = data.frame(sample = names(res$peak_vaf), peak_vaf = unname(res$peak_vaf), removed = names(res$peak_vaf) %in% res$mixed_samples)
}
write.table(peak, snakemake@output$peak_vaf, sep = "\t", quote = FALSE, row.names = FALSE)
