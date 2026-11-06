# Check the donor's sex against the data (called by Snakemake rule sex_check)
log = file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
library(phylobuild)

cm = import_cgpvaf_SNV_and_INDEL(snakemake@input$snv, snakemake@input$indel, ref_ID = snakemake@params$normal)
res = check_donor_sex(cm, configured_sex = snakemake@params$sex, build = snakemake@params$build,
                      strict = snakemake@params$mode == "strict")
write.table(cbind(res$summary, configured = res$configured, inferred = res$inferred, used = res$sex, ok = res$ok),
            snakemake@output$summary, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res$per_colony, snakemake@output$per_colony, sep = "\t", quote = FALSE, row.names = FALSE)
plot_sex_check(res, snakemake@output$plot)
writeLines(res$sex, snakemake@output$sex)
