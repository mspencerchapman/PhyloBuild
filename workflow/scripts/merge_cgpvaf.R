# Merge the per-batch cgpVAF files for one donor (called by Snakemake rule cgpvaf_merge)
log = file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
library(phylobuild)

files = unlist(lapply(snakemake@params$batch_dirs, function(d) list.files(d, pattern = "_vaf\\.tsv$", full.names = TRUE)))
message("Merging ", length(files), " cgpVAF files:\n", paste(files, collapse = "\n"))
if (length(files) != length(snakemake@params$batch_dirs)) stop("Expected one *_vaf.tsv per batch folder")

merged = merge_cgpvaf_files(files)
found = sub("_MTR$", "", grep("_MTR$", colnames(merged), value = TRUE))
missing = setdiff(snakemake@params$expected_samples, found)
if (length(missing) > 0) stop("Samples missing from the cgpVAF output: ", paste(missing, collapse = ", "))
write_merged_cgpvaf(merged, snakemake@output[[1]])
message("Wrote ", nrow(merged), " sites x ", length(found), " samples (including the in silico normal)")
