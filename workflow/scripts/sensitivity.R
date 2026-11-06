# Per-sample sensitivity from the proportion of germline variants called (rule sensitivity)
log = file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
library(phylobuild)

x = read_mats_and_params(snakemake@input$mats)
samples = snakemake@params$samples
snv_files = setNames(unlist(snakemake@input$snv_calls), samples)
indel_files = setNames(unlist(snakemake@input$indel_calls), samples)
keep = samples %in% colnames(x$COMB_mats$NV) #samples excluded upstream are left out
res = germline_sensitivity(x$COMB_mats, x$filter_params, snv_files[keep], indel_files[keep], depth_cutoff = snakemake@params$depth_cutoff)
print(summary(res))
write.table(res, snakemake@output[[1]], row.names = FALSE) #space separated, as before
