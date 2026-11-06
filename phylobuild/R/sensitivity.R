# Hairpin-filter reduction and germline-based sensitivity (ported from Reducing_mutset_from_MSfilters.R and
# Sensitivity_analysis_from_SNPs.R).

#Read Chrom/Pos/Ref/Alt from a VCF body or bed-like file and return "Chrom-Pos-Ref-Alt" mut_refs
read_mut_refs = function(file, skip = 0) {
  if (!file.exists(file)) stop("Call file not found: ", file)
  if (file.info(file)$size == 0) return(character(0))
  calls = data.table::fread(file, skip = skip, header = FALSE, sep = "\t", colClasses = "character", data.table = FALSE)
  pos = gsub("^\\s+|\\s+$", "", calls[, 2])
  if (ncol(calls) >= 5) paste(calls[, 1], pos, calls[, 4], calls[, 5], sep = "-") else paste(calls[, 1], pos, calls[, 3], calls[, 4], sep = "-")
}

#Keep only SNVs that passed the hairpin (LCM artefact) filter; all indels are kept.
#Removes the need to re-run cgpVAF on the reduced set.
reduce_to_hairpin_pass = function(COMB_mats, filter_params, hairpin_pass_mut_refs) {
  in_matrix = mean(hairpin_pass_mut_refs %in% COMB_mats$mat$mut_ref)
  message(round(100 * in_matrix, 1), "% of hairpin-pass SNVs are in the cgpVAF matrix")
  if (in_matrix < 0.9) warning("Fewer than 90% of hairpin-pass SNVs are in the matrix - check the chromosome naming matches")
  keep = COMB_mats$mat$mut_ref %in% hairpin_pass_mut_refs | COMB_mats$mat$Mut_type == "INDEL"
  list(COMB_mats = list_subset(COMB_mats, keep), filter_params = filter_params[keep, ])
}

#Germline variants: high germline p-value, low overdispersion, well covered
define_germline = function(COMB_mats, filter_params, mut_type, binom_cutoff = 6, rho_cutoff = 0.05, depth_cutoff = 500) {
  COMB_mats$mat$mut_ref[COMB_mats$mat$Mut_type == mut_type & log10(filter_params$germline_pval) > -binom_cutoff &
                          filter_params$bb_rhoval < rho_cutoff & rowSums(COMB_mats$NR) > depth_cutoff]
}

#Per-sample sensitivity: the proportion of germline variants that were called in that sample.
#  snv_files / indel_files: named vectors of per-sample call files (names = sample IDs)
germline_sensitivity = function(COMB_mats, filter_params, snv_files, indel_files = NULL, binom_cutoff = 6, rho_cutoff = 0.05, depth_cutoff = 500) {
  germline_SNVs = define_germline(COMB_mats, filter_params, "SNV", binom_cutoff, rho_cutoff, depth_cutoff)
  germline_INDELs = define_germline(COMB_mats, filter_params, "INDEL", binom_cutoff, rho_cutoff, depth_cutoff)
  message(length(germline_SNVs), " germline SNVs and ", length(germline_INDELs), " germline INDELs used to estimate sensitivity")
  if (length(germline_SNVs) == 0) stop("No germline SNVs found - is depth_cutoff (", depth_cutoff, ") too high for this number of samples?")

  pick_up = function(files, germline) sapply(files, function(f) {
    if (length(germline) == 0) return(NA_real_)
    mean(germline %in% read_mut_refs(f))
  })
  snv = pick_up(snv_files, germline_SNVs)
  if (is.null(indel_files)) return(data.frame(Sample = names(snv_files), SNV_sensitivity = unname(snv)))

  indel = pick_up(indel_files, germline_INDELs)
  both = intersect(names(snv_files), names(indel_files))
  data.frame(Sample = both, SNV_sensitivity = unname(snv[both]), INDEL_sensitivity = unname(indel[both]))
}
