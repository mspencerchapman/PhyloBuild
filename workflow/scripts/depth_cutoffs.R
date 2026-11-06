# Mean-depth window for the tree step: "auto" from the data, or fixed values from the config (rule depth_cutoffs)
log = file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
library(phylobuild)

x = read_mats_and_params(snakemake@input$mats)
setting = snakemake@params$setting
auto = depth_cutoffs_auto(x$filter_params, x$COMB_mats, k = snakemake@params$k)
if (identical(setting, "auto")) {
  cutoffs = auto
} else {
  cutoffs = c(unlist(setting)[c("auto_low", "auto_high", "xy_low", "xy_high")], median_depth = auto[["median_depth"]])
}
message("Median autosomal depth ", round(auto[["median_depth"]], 1), " (MAD ", round(auto[["mad_depth"]], 2), "); using cutoffs: ", paste(names(cutoffs), cutoffs, sep = "=", collapse = ", "))
write.table(data.frame(name = names(cutoffs), value = unname(cutoffs)), snakemake@output$tsv, sep = "\t", quote = FALSE, row.names = FALSE)
plot_depth_cutoffs(x$filter_params, x$COMB_mats, cutoffs, snakemake@output$plot)
