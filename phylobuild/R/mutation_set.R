# Applying the filter cut-offs, genotyping each sample, and building the dummy DNA strings for MPBoot
# (ported from foetal.filters.parallel.R).

#Apply a pair of thresholds (autosomes, X/Y) to a value: haploid rows in males (Y, non-PAR X) use the second threshold.
#Returns 1 (pass) / 0 (fail), as the original assess_* functions did.
pass_by_chrom_class = function(value, COMB_mats, threshold_auto, threshold_xy) {
  xy = haploid_rows(COMB_mats)
  ifelse(value >= ifelse(xy, threshold_xy, threshold_auto), 1, 0)
}

assess_mean_depth = function(filter_params, COMB_mats, AUTO_low_depth_cutoff, AUTO_high_depth_cutoff, XY_low_depth_cutoff, XY_high_depth_cutoff) {
  xy = haploid_rows(COMB_mats)
  low = ifelse(xy, XY_low_depth_cutoff, AUTO_low_depth_cutoff)
  high = ifelse(xy, XY_high_depth_cutoff, AUTO_high_depth_cutoff)
  ifelse(filter_params$mean_depth >= low & filter_params$mean_depth <= high, 1, 0)
}

get_filtered_mut_set = function(input_set_ID,
                                COMB_mats,
                                filter_params,
                                gender,

                                ##PARAMETERS FOR FILTERING THE MUTATION SET. These arguments set the cut-offs for selecting the set of "true somatic mutations"
                                #NB. If parameter is set to NA, filter will not be applied.
                                retain_muts = NA, #vector of "mut_refs" (in format Chrom-Pos-Ref-Alt) to retain even if fail the filtering
                                exclude_muts = NA, #vector of "mut_refs" (in format Chrom-Pos-Ref-Alt) to exclude even if pass the filtering
                                germline_pval = -10,
                                rho = 0.1, #Beta-binomial filter, rho value cut-off (i.e. rho must be > cut-off)
                                mean_depth = NA, #Numeric vector of length 4 in the order (1) AUTO min depth cutoff, (2) AUTO max depth cutoff, (3) XY min depth cutoff, (4) XY max depth cutoff
                                pval_dp2 = NA, #Numeric vector of length 1.
                                pval_dp3 = 0.01, #Numeric vector of length 1.
                                min_depth = c(6, 4), #Numeric vector of length 2 for (1) AUTO and (2) XY muts. Minimum depth that at least one positive sample must have.
                                min_pval_for_true_somatic = 0.1, #Numeric vector of length 1. At least one sample must have a p-value of >this cutoff for coming from a true somatic distribution.
                                min_vaf = c(0.2, 0.8), #Numeric vector of length 2 for (1) AUTO and (2) XY muts. At least one sample must meet the minimum vaf for the mutation.

                                #These arguments decide the thresholds for genotyping each sample. This may be less stringent.
                                min_variant_reads_SHARED = 2, #should be equal to, or more relaxed than the above. At least one parameter must be used.
                                min_pval_for_true_somatic_SHARED = 0.05,
                                min_vaf_SHARED = c(0.2, 0.7) #Numeric vector of length 2 for AUTO and XY muts
) {
  COMB_mats$gender = gender
  #List (1) the names of each of the initial filters that can be applied,
  #(2) the names of their input parameters,
  #(3) the variable name that must be present in the filter_params dataframe if filter is being used.
  filter_name = c("Germline filter", "Beta-binomial filter", "Mean depth filter", "P-value within positive samples [defining positive as >=2 variant reads] filter", "P-value within positive samples [defining positive as >=3 variant reads] filter", "Minimum depth filter", "Minimum p-value for true somatic mutation filter", "Minimum VAF filter")
  set_parameter = list(germline_pval, rho, mean_depth, pval_dp2, pval_dp3, min_depth, min_pval_for_true_somatic, min_vaf)
  var_name = c("germline_pval", "bb_rhoval", "mean_depth", "pval_within_pos_dp2", "pval_within_pos_dp3", "max_depth_in_pos_samples", "max_pval_in_pos_sample", "max_mut_vaf")

  #Functions that decide if each mutation is a pass (1) or fail (0) for each filter
  test_function = list(function(x) {ifelse(log10(filter_params$germline_pval) < x, 1, 0)},
                       function(x) {ifelse(filter_params$bb_rhoval > x, 1, 0)},
                       function(x) {assess_mean_depth(filter_params, COMB_mats, x[1], x[2], x[3], x[4])},
                       function(x) {ifelse(filter_params$pval_within_pos_dp2 > x, 1, 0)},
                       function(x) {ifelse(filter_params$pval_within_pos_dp3 > x, 1, 0)},
                       function(x) {pass_by_chrom_class(filter_params$max_depth_in_pos_samples, COMB_mats, x[1], x[2])},
                       function(x) {ifelse(filter_params$max_pval_in_pos_sample > x, 1, 0)},
                       function(x) {pass_by_chrom_class(filter_params$max_mut_vaf, COMB_mats, x[1], x[2])}
  )

  #CHECK THE INPUT DATA is all present for the set filters
  if (any(!c("NR", "NV") %in% names(COMB_mats)) || !isTRUE(all.equal(dim(COMB_mats$NR), dim(COMB_mats$NV)))) {
    stop("COMB_mats object must contain matched NR and NV matrices of equal dimensions")
  }
  if (nrow(filter_params) != nrow(COMB_mats$NV)) {
    stop("The filter_params data frame and the COMB_mats$NV and NR matrices must contain the same number of mutations")
  }
  for (i in 1:length(filter_name)) {
    if (!is.na(set_parameter[[i]][1]) & !any(names(filter_params) == var_name[i])) {
      stop(paste(filter_name[i], "needs a variable named", var_name[i], "in the filter_params data frame. Update the filter_params data frame or set the appropriate argument to NA"))
    }
  }
  if ((!is.na(min_pval_for_true_somatic_SHARED[1]) | !is.na(min_pval_for_true_somatic[1])) & !any(names(COMB_mats) == "PVal")) {
    stop("If using the 'min_pval_for_true_somatic' or 'min_pval_for_true_somatic_SHARED' parameters, the COMB_mats list must contain a matrix called PVal")
  }

  #Apply each of the filter functions; for those with NA parameters, a vector of 1's is returned (i.e. all mutations "pass" the filter)
  out = mapply(FUN = function(param, test_function) {
    if (is.na(param[1])) rep(1, nrow(filter_params)) else test_function(param)
  }, param = set_parameter, test_function = test_function, SIMPLIFY = FALSE)

  filter_pass = Reduce(cbind, out); rownames(filter_pass) = rownames(filter_params); colnames(filter_pass) = var_name
  select_muts = (apply(filter_pass, 1, function(x) all(x == 1)) & !rownames(filter_pass) %in% exclude_muts) | rownames(filter_pass) %in% retain_muts #Mutations must pass ALL the applied filters.
  filter_code = apply(filter_pass, 1, paste, collapse = "-") #Quick record of which filter is removing most mutations.
  COMB_mats.tree.build = list_subset(COMB_mats, select_vector = select_muts) #Subset matrices to include only the PASS mutations

  #Set the rownames to the mut_ref, and colnames to the sample names (without MTR or DEP)
  rownames(COMB_mats.tree.build$NV) = rownames(COMB_mats.tree.build$NR) = rownames(COMB_mats.tree.build$PVal) = COMB_mats.tree.build$mat$mut_ref
  colnames = gsub(pattern = "_MTR", replacement = "", x = colnames(COMB_mats.tree.build$NV))
  colnames(COMB_mats.tree.build$NV) = colnames(COMB_mats.tree.build$NR) = colnames(COMB_mats.tree.build$PVal) = colnames

  #BUILD THE GENOTYPE MATRIX
  #"Definite positive" samples get genotype 1: each criterion is a matrix of 1/0 (or 1 if not used)
  if (!is.na(min_variant_reads_SHARED[1])) {min_variant_reads_mat = COMB_mats.tree.build$NV >= min_variant_reads_SHARED} else {min_variant_reads_mat = 1}
  if (!is.na(min_pval_for_true_somatic_SHARED[1])) {min_pval_for_true_somatic_mat = COMB_mats.tree.build$PVal > min_pval_for_true_somatic_SHARED} else {min_pval_for_true_somatic_mat = 1}
  if (!is.na(min_vaf_SHARED[1]) & gender == "female") {
    depth_no_zero = COMB_mats.tree.build$NR
    depth_no_zero[depth_no_zero == 0] = 1
    min_vaf_mat = COMB_mats.tree.build$NV / depth_no_zero > min_vaf_SHARED[1]
  } else if (!is.na(min_vaf_SHARED[1]) & gender == "male") {
    min_vaf_mat = matrix(0, ncol = ncol(COMB_mats.tree.build$NV), nrow = nrow(COMB_mats.tree.build$NV))
    xy_muts = haploid_rows(COMB_mats.tree.build)
    depth_no_zero = COMB_mats.tree.build$NR
    depth_no_zero[depth_no_zero == 0] = 1
    min_vaf_mat[xy_muts, ] = as.matrix(COMB_mats.tree.build$NV[xy_muts, ] / depth_no_zero[xy_muts, ] > min_vaf_SHARED[2])
    min_vaf_mat[!xy_muts, ] = as.matrix(COMB_mats.tree.build$NV[!xy_muts, ] / depth_no_zero[!xy_muts, ] > min_vaf_SHARED[1])
  } else {min_vaf_mat = 1}

  COMB_mats.tree.build$Genotype_bin = min_variant_reads_mat * min_pval_for_true_somatic_mat * min_vaf_mat

  #"Not sure" samples get genotype 0.5. These thresholds were set empirically.
  COMB_mats.tree.build$Genotype_bin[COMB_mats.tree.build$NV > 0 & COMB_mats.tree.build$PVal > 0.01 & COMB_mats.tree.build$Genotype_bin != 1] = 0.5 #Any mutant reads: "?" as long as p-value > 0.01
  COMB_mats.tree.build$Genotype_bin[COMB_mats.tree.build$NV >= 3 & COMB_mats.tree.build$PVal > 0.001 & COMB_mats.tree.build$Genotype_bin != 1] = 0.5 #Many mutant reads: "?" even if p-value is low (may be biased sequencing)
  COMB_mats.tree.build$Genotype_bin[(COMB_mats.tree.build$NV == 0) & (COMB_mats.tree.build$PVal > 0.05)] = 0.5 #Inadequate depth to exclude the mutation, even with no variant reads

  Genotype_shared_bin = COMB_mats.tree.build$Genotype_bin[rowSums(COMB_mats.tree.build$Genotype_bin == 1) > 1, , drop = FALSE]

  params = list(input_set_ID = input_set_ID, gender = gender, retain_muts = retain_muts, exclude_muts = exclude_muts,
                germline_pval = germline_pval, rho = rho, mean_depth = mean_depth, pval_dp2 = pval_dp2, pval_dp3 = pval_dp3,
                min_depth = min_depth, min_pval_for_true_somatic = min_pval_for_true_somatic, min_vaf = min_vaf,
                min_variant_reads_SHARED = min_variant_reads_SHARED, min_pval_for_true_somatic_SHARED = min_pval_for_true_somatic_SHARED,
                min_vaf_SHARED = min_vaf_SHARED)

  summary = data.frame(total_mutations = sum(select_muts),
                       total_SNVs = sum(COMB_mats.tree.build$mat$Mut_type == "SNV"),
                       total_INDELs = sum(COMB_mats.tree.build$mat$Mut_type == "INDEL"),
                       shared_mutations = nrow(Genotype_shared_bin),
                       shared_SNVs = sum(COMB_mats.tree.build$mat$mut_ref %in% rownames(Genotype_shared_bin) & COMB_mats.tree.build$mat$Mut_type == "SNV"),
                       shared_INDELs = sum(COMB_mats.tree.build$mat$mut_ref %in% rownames(Genotype_shared_bin) & COMB_mats.tree.build$mat$Mut_type == "INDEL"))
  dna_strings = dna_strings_from_genotype(Genotype_shared_bin)

  message(summary$total_mutations, " total mutations (", summary$total_SNVs, " SNVs, ", summary$total_INDELs, " INDELs); ",
          summary$shared_mutations, " shared (", summary$shared_SNVs, " SNVs, ", summary$shared_INDELs, " INDELs)")

  list(COMB_mats.tree.build = COMB_mats.tree.build, Genotype_shared_bin = Genotype_shared_bin, filter_code = filter_code, params = params, summary = summary, dna_strings = dna_strings)
}

#Dummy DNA strings from the shared genotype matrix: W = wild type, V = variant, ? = uncertain
dna_strings_from_genotype = function(genotype_mat) {
  Muts = rownames(genotype_mat)
  Ref = rep("W", length(Muts))
  Alt = rep("V", length(Muts))

  dna_strings = list()
  dna_strings[1] = paste(Ref, sep = "", collapse = "")
  for (k in 1:ncol(genotype_mat)) {
    Mutations = Ref
    Mutations[genotype_mat[, k] == 0.5] = "?"
    Mutations[genotype_mat[, k] == 1] = Alt[genotype_mat[, k] == 1]
    dna_strings[k + 1] = paste(Mutations, sep = "", collapse = "")
  }
  names(dna_strings) = c("Ancestral", colnames(genotype_mat))
  return(dna_strings)
}

#Minimal VCF body (no header) from the "mat" object
create_vcf_files = function(mat, select_vector = NULL) {
  if (is.null(select_vector)) {vcf_file = mat[, c("Chrom", "Pos", "Ref", "Alt")]} else {vcf_file = mat[select_vector, c("Chrom", "Pos", "Ref", "Alt")]}
  names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")
  vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."
  vcf_file = vcf_file[, c(1, 2, 8, 3, 4, 7, 6, 5)]
  return(vcf_file)
}
