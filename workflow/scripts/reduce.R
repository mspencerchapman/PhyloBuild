# Keep SNVs that passed the hairpin filter; then drop near-certain germline variants (rule reduce)
log = file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
library(phylobuild)

x = read_mats_and_params(snakemake@input$mats)
hairpin_pass = read_mut_refs(snakemake@input$hairpin_sites)
post = reduce_to_hairpin_pass(x$COMB_mats, x$filter_params, hairpin_pass)
legacy = if (isTRUE(snakemake@params$legacy_rdata)) TRUE
save_object(post, snakemake@output$postMS, legacy_name = legacy)
reduced = reduce_mutation_set(post$COMB_mats, post$filter_params)
save_object(reduced, snakemake@output$reduced, legacy_name = legacy)
message(nrow(x$filter_params), " mutations -> ", nrow(post$filter_params), " after hairpin filter -> ", nrow(reduced$filter_params), " after removing near-certain germline")
