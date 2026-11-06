# Reading and merging cgpVAF output.
# Columns are found by name (<sample>_MTR / <sample>_DEP), so this works for any number of
# samples per batch and any cgpVAF column layout. It replaces the old `cut -f 24,26,39,...` merge.

#Read one cgpVAF *_vaf.tsv file: skips the "##" header lines and keeps Chrom/Pos/Ref/Alt + MTR/DEP columns
read_cgpvaf_tsv = function(file) {
  con = file(file, "r")
  n_skip = 0
  while (length(line <- readLines(con, n = 1)) > 0 && startsWith(line, "##")) n_skip = n_skip + 1
  close(con)
  dt = data.table::fread(file, skip = n_skip, header = TRUE, sep = "\t", data.table = FALSE, check.names = FALSE)
  keep = c("Chrom", "Pos", "Ref", "Alt", grep("_(MTR|DEP)$", colnames(dt), value = TRUE))
  missing = setdiff(c("Chrom", "Pos", "Ref", "Alt"), colnames(dt))
  if (length(missing) > 0) stop(file, " is missing column(s): ", paste(missing, collapse = ", "))
  dt[, keep]
}

#Merge the per-batch cgpVAF files for one donor into a single matrix.
#Output layout matches the legacy merged_SNVs_*.tsv: Chrom Pos Ref Alt, then <sample>_MTR <sample>_DEP pairs
#(the in silico normal first, once).
merge_cgpvaf_files = function(files) {
  if (length(files) == 0) stop("No cgpVAF files supplied to merge")
  tables = lapply(files, read_cgpvaf_tsv)
  key = function(t) paste(t$Chrom, t$Pos, t$Ref, t$Alt, sep = "-")
  ref_key = key(tables[[1]])
  if (anyDuplicated(ref_key)) stop("Duplicated variants in ", files[1])
  merged = tables[[1]]
  for (i in seq_along(tables)[-1]) {
    t = tables[[i]]
    k = key(t)
    if (!setequal(k, ref_key) || length(k) != length(ref_key)) {
      stop("Variant sites differ between ", files[1], " and ", files[i], " - were they run with the same bed file?")
    }
    t = t[match(ref_key, k), , drop = FALSE]
    new_cols = setdiff(grep("_(MTR|DEP)$", colnames(t), value = TRUE), colnames(merged)) #drops the repeated normal columns
    merged = cbind(merged, t[, new_cols, drop = FALSE])
  }
  merged
}

write_merged_cgpvaf = function(merged, file) {
  data.table::fwrite(merged, file, sep = "\t", quote = FALSE)
}

#Import a merged cgpVAF matrix into the COMB_mats structure (mat, NV, NR), dropping the in silico normal
import_cgpvaf_output = function(cgpvaf_output_file, ref_ID) {
  mat = read.delim(cgpvaf_output_file, stringsAsFactors = FALSE, check.names = FALSE)
  normal_cols = grepl(paste0("^", ref_ID, "_"), colnames(mat))
  if (!any(normal_cols)) warning("No columns for the in silico normal '", ref_ID, "' found in ", cgpvaf_output_file)
  mat = mat[, !normal_cols]
  NV = mat[, grepl("_MTR$", colnames(mat)), drop = FALSE]
  NR = mat[, grepl("_DEP$", colnames(mat)), drop = FALSE]
  colnames(NV) = gsub("_MTR$", "", colnames(NV))
  colnames(NR) = gsub("_DEP$", "", colnames(NR))
  if (!identical(colnames(NV), colnames(NR))) stop("MTR and DEP columns don't match in ", cgpvaf_output_file)
  mat$mut_ref = paste(mat$Chrom, mat$Pos, mat$Ref, mat$Alt, sep = "-")
  mat = mat[, c("Chrom", "Pos", "Ref", "Alt", "mut_ref")]
  rownames(NV) = rownames(NR) = mat$mut_ref
  list(mat = mat, NV = NV, NR = NR)
}

#Import SNV and (optionally) indel matrices and combine them, keeping only samples present in both
import_cgpvaf_SNV_and_INDEL = function(SNV_output_file, INDEL_output_file = NULL, ref_ID) {
  SNV_mats = import_cgpvaf_output(SNV_output_file, ref_ID = ref_ID)
  SNV_mats$mat$Mut_type = "SNV"
  if (is.null(INDEL_output_file)) return(SNV_mats)

  INDEL_mats = import_cgpvaf_output(INDEL_output_file, ref_ID = ref_ID)
  INDEL_mats$mat$Mut_type = "INDEL"

  samples_in_both = intersect(colnames(SNV_mats$NV), colnames(INDEL_mats$NV))
  only_one = setdiff(union(colnames(SNV_mats$NV), colnames(INDEL_mats$NV)), samples_in_both)
  if (length(only_one) > 0) message("Dropping samples not in both SNV and indel matrices: ", paste(only_one, collapse = ", "))
  message(length(samples_in_both), " samples in both SNV and INDEL cgpVAF output matrices, and will be combined")
  for (m in c("NV", "NR")) {
    SNV_mats[[m]] = SNV_mats[[m]][, samples_in_both, drop = FALSE]
    INDEL_mats[[m]] = INDEL_mats[[m]][, samples_in_both, drop = FALSE]
  }
  mapply(SNV_mats, INDEL_mats, FUN = rbind, SIMPLIFY = FALSE)
}
