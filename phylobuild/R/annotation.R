# Mutation annotation with VAGrENT (Sanger). Optional: build_tree() works without it.

VAGRENT_HEADER = c(
  "##fileformat=VCFv4.1",
  '##INFO=<ID=VD,Number=1,Type=String,Description="Vagrent Default Annotation">',
  '##INFO=<ID=VW,Number=1,Type=String,Description="Vagrent Most Deleterious Annotation">',
  '##INFO=<ID=VT,Number=1,Type=String,Description="Variant type based on the Vagrent Default Annotation">',
  '##INFO=<ID=VC,Number=.,Type=String,Description="Variant consequence based on the Vagrent Default Annotation">'
)

#Split the VAGrENT default annotation (VD) into columns
split_vagrent_output = function(vd, col_IDs = c("Gene", "Transcript", "RNA", "CDS", "Protein", "Type", "SO_codes")) {
  parts = stringr::str_split(vd, pattern = "\\|", simplify = TRUE)
  parts = parts[, seq_len(min(ncol(parts), length(col_IDs))), drop = FALSE]
  out = as.data.frame(matrix("", nrow = length(vd), ncol = length(col_IDs), dimnames = list(NULL, col_IDs)), stringsAsFactors = FALSE)
  out[, seq_len(ncol(parts))] = parts
  out
}

#Returns a function(mat) that annotates mutations with VAGrENT, for use as build_tree(annotate = ...).
#  work_prefix: path prefix for the VCF files written along the way
vagrent_annotator = function(work_prefix, cache, species = "Human", assembly = "GRCh38_full_analysis_set_plus_decoy_hla", command = "AnnotateVcf.pl") {
  #VAGrENT runs in a container that mounts /lustre but not /nfs: give it the cache's real path, not a symlink
  cache = normalizePath(cache, mustWork = FALSE)
  function(mat) {
    vcf_in = paste0(work_prefix, "_header.vcf")
    vcf_out = paste0(vcf_in, ".annot")
    writeLines(VAGRENT_HEADER, vcf_in)
    suppressWarnings(write.table(create_vcf_files(mat), file = vcf_in, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
    args = c("-i", vcf_in, "-o", vcf_out, "-sp", species, "-as", assembly, "-c", cache)
    message("Running: ", command, " ", paste(args, collapse = " "))
    status = system2(command, args)
    if (status != 0 || !file.exists(vcf_out)) stop("VAGrENT annotation failed (exit status ", status, ")")

    annot = data.table::fread(vcf_out, skip = "#CHROM", data.table = FALSE, colClasses = "character")
    #INFO is "VT=...;VD=...;VC=...;VW=..." - parse by key rather than position
    vd = sub("^VD=", "", sapply(strsplit(annot$INFO, ";"), function(x) {v = grep("^VD=", x, value = TRUE); if (length(v)) v[1] else "VD="}))
    key_out = paste(annot[["#CHROM"]], annot$POS, annot$REF, annot$ALT, sep = "-")
    out = split_vagrent_output(vd)[match(mat$mut_ref, key_out), ]
    if (anyNA(match(mat$mut_ref, key_out))) warning(sum(is.na(match(mat$mut_ref, key_out))), " mutations missing from the VAGrENT output")
    rownames(out) = NULL
    out
  }
}
