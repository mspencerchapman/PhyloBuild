# Per-mutation filtering parameters (ported from foetal.filters.parallel.R and
# Mutation_filtering_get_parameters_v2.R).
#
# The statistics are the same as the originals. The loops of binom.test() calls are replaced with
# pbinom(): for alternative = "less", binom.test(x, n, p)$p.value is exactly pbinom(x, n, p).
# Sex is passed in explicitly (from the donor config, checked by the sex check step) instead of being
# guessed from the number of Y chromosome calls.

#Function to subset the list of matrices
list_subset = function(list, select_vector) {
  for (i in 1:length(list)) {
    if (!is.null(dim(list[[i]]))) {
      list[[i]] = list[[i]][select_vector, , drop = FALSE]
    }
  }
  return(list)
}

#Expected VAF of a heterozygous germline/clonal somatic variant for each row: 0.95 at haploid sites in males
#(Y and non-PAR X), else 0.5. (The original treated the X PARs as haploid too.)
expected_vaf = function(COMB_mats) {
  ifelse(haploid_rows(COMB_mats), 0.95, 0.5)
}

#Minimum variant reads for a sample to count as positive: 2 at haploid sites in males, otherwise 3
min_variant_reads = function(COMB_mats, auto = 3, xy = 2) {
  ifelse(haploid_rows(COMB_mats), xy, auto)
}

#P-value matrix: probability of the observed variant reads (or fewer) if the sample carried the mutation clonally
pval_matrix = function(COMB_mats) {
  NV = as.matrix(COMB_mats$NV); NR = as.matrix(COMB_mats$NR)
  NR[NR == 0] = 1
  p = matrix(expected_vaf(COMB_mats), nrow = nrow(NV), ncol = ncol(NV))
  pval_mat = matrix(pbinom(NV, NR, p), nrow = nrow(NV), ncol = ncol(NV), dimnames = dimnames(NV))
  return(pval_mat)
}

#Exact binomial test of the pooled read counts against the germline expectation. Returns a p-value per mutation.
germline.binomial.filter = function(COMB_mats) {
  NV_vec = rowSums(COMB_mats$NV)
  NR_vec = rowSums(COMB_mats$NR)
  pval = rep(1, length(NV_vec))
  covered = NR_vec > 0
  pval[covered] = pbinom(NV_vec[covered], NR_vec[covered], expected_vaf(COMB_mats)[covered])
  return(pval)
}

#Beta-binomial overdispersion. Returns the maximum likelihood rho per mutation, on the same grid as the original.
#Works best on sets of clonal samples (ideally >10 or so).
beta.binom.filter = function(COMB_mats, chunk_size = 5000) {
  NV = as.matrix(COMB_mats$NV); NR = as.matrix(COMB_mats$NR)
  NR[NR == 0] = 1
  rhovec = 10^seq(-6, -0.05, by = 0.05) # rho will be bounded within 1e-6 and 0.89
  rho_est = rep(NA_real_, nrow(NV))
  for (start in seq(1, nrow(NV), by = chunk_size)) {
    rows = start:min(start + chunk_size - 1, nrow(NV))
    nv = NV[rows, , drop = FALSE]; nr = NR[rows, , drop = FALSE]
    mu = matrix(rowSums(nv) / rowSums(nr), nrow = length(rows), ncol = ncol(nv))
    ll = sapply(rhovec, function(rhoj) rowSums(matrix(VGAM::dbetabinom(x = nv, size = nr, rho = rhoj, prob = mu, log = TRUE), nrow = length(rows))))
    ll = matrix(ll, nrow = length(rows))
    rho_est[rows] = rhovec[apply(ll, 1, which.max)]
  }
  return(rho_est)
}

#Pooled binomial test across "positive" samples (>= define_pos variant reads). Filters mutations shared by
#several samples whose combined VAF is too low for a true clonal somatic mutation.
low_vaf_in_pos_samples = function(COMB_mats, define_pos) {
  NV = as.matrix(COMB_mats$NV); NR = as.matrix(COMB_mats$NR)
  pos = NV >= define_pos
  any_pos = rowSums(pos) > 0
  pval = rep(0, nrow(NV))
  pval[any_pos] = pbinom(rowSums(NV * pos)[any_pos], rowSums(NR * pos)[any_pos], expected_vaf(COMB_mats)[any_pos])
  return(pval)
}
low_vaf_in_pos_samples_dp2 = function(COMB_mats) low_vaf_in_pos_samples(COMB_mats, define_pos = 2)
low_vaf_in_pos_samples_dp3 = function(COMB_mats) low_vaf_in_pos_samples(COMB_mats, define_pos = 3)

get_mean_depth = function(COMB_mats) rowMeans(COMB_mats$NR)

#Row-wise maximum of `values` over positive samples; -Inf where no sample is positive (as max() of nothing)
max_in_positive_samples = function(values, COMB_mats) {
  NV = as.matrix(COMB_mats$NV)
  pos = NV >= min_variant_reads(COMB_mats)
  values = as.matrix(values)
  values[!pos] = -Inf
  apply(values, 1, max)
}
get_max_depth_in_pos = function(COMB_mats) max_in_positive_samples(COMB_mats$NR, COMB_mats)
get_max_pval_in_pos = function(COMB_mats) max_in_positive_samples(COMB_mats$PVal, COMB_mats)

get_max_vaf = function(COMB_mats) {
  NR = as.matrix(COMB_mats$NR)
  NR[NR == 0] = 1
  apply(as.matrix(COMB_mats$NV) / NR, 1, max)
}

#Remove samples (low coverage or listed), then mutations no longer positive in any remaining sample.
#Returns COMB_mats, or list(COMB_mats, filter_params) if filter_params is given.
remove_low_coverage_samples = function(COMB_mats,
                                       filter_params = NULL,
                                       min_sample_mean_cov,
                                       other_samples_to_remove = NULL,
                                       min_variant_reads_auto = 3,
                                       min_variant_reads_xy = 2) {
  mean_sample_cov = colMeans(COMB_mats$NR)
  low_cov_names = gsub(pattern = "_DEP", replacement = "", x = names(mean_sample_cov)[mean_sample_cov < min_sample_mean_cov])
  remove = gsub(pattern = "_MTR", replacement = "", x = colnames(COMB_mats$NV)) %in% c(low_cov_names, other_samples_to_remove)
  unknown = setdiff(other_samples_to_remove, gsub("_MTR", "", colnames(COMB_mats$NV)))
  if (length(unknown) > 0) warning("Samples to exclude not found in the matrix: ", paste(unknown, collapse = ", "))
  if (any(remove)) {
    message("Removing ", sum(remove), " samples: ", paste(colnames(COMB_mats$NV)[remove], collapse = ", "))
    COMB_mats$NV = COMB_mats$NV[, !remove, drop = FALSE]
    COMB_mats$NR = COMB_mats$NR[, !remove, drop = FALSE]
    if (!is.null(COMB_mats$PVal)) COMB_mats$PVal = COMB_mats$PVal[, !remove, drop = FALSE]

    null_remove = rowSums(COMB_mats$NV >= min_variant_reads_auto | (COMB_mats$NV >= min_variant_reads_xy & haploid_rows(COMB_mats))) == 0
    message(sum(null_remove), " mutations removed as no positives in any remaining samples.")
    COMB_mats = list_subset(COMB_mats, select_vector = !null_remove)
    if (!is.null(filter_params)) filter_params = filter_params[!null_remove, ]
  } else {
    message("No samples removed")
  }
  if (is.null(filter_params)) COMB_mats else list(COMB_mats = COMB_mats, filter_params = filter_params)
}

#Peak of the VAF density of a sample's likely-somatic mutations - screens for non-clonal (mixed) colonies
check_peak_vaf = function(sample, COMB_mats, filter_params, rho_cutoff = 0.3) {
  colnames(COMB_mats$NV) = gsub(pattern = "_MTR", replacement = "", x = colnames(COMB_mats$NV))
  dens = density((COMB_mats$NV / COMB_mats$NR)[COMB_mats$NV[, sample] >= 2 & !is_xy(COMB_mats$mat$Chrom) & log10(filter_params$germline_pval) < (-10) & filter_params$bb_rhoval > rho_cutoff, sample])
  return(dens$x[which.max(dens$y)])
}

vaf_density_plot = function(sample, COMB_mats, filter_params, rho_cutoff = 0.3) {
  colnames(COMB_mats$NV) = colnames(COMB_mats$NR) = gsub(pattern = "_MTR", replacement = "", x = colnames(COMB_mats$NV))
  sample_mean_depth = mean(COMB_mats$NR[, sample])
  dens = density((COMB_mats$NV / COMB_mats$NR)[COMB_mats$NV[, sample] >= 2 & !is_xy(COMB_mats$mat$Chrom) & log10(filter_params$germline_pval) < (-10) & filter_params$bb_rhoval > rho_cutoff, sample])
  plot(dens, main = sample, xlim = c(-0.05, 1.05))
  abline(v = dens$x[which.max(dens$y)])
  text(0.7, max(dens$y) - 0.2, paste("Peak VAF dens=", round(dens$x[which.max(dens$y)], digits = 2), "\nMean coverage=", round(sample_mean_depth, digits = 2)), col = "red", cex = 0.7)
}

#Remove mutations with no sample meeting the variant read threshold
remove_null_mutations = function(COMB_mats, extra = NULL) {
  null_remove = rowSums(as.matrix(COMB_mats$NV) >= min_variant_reads(COMB_mats)) == 0
  message(sum(null_remove), " mutations removed as no sample meets minimum threshold of number of variant reads")
  out = list(COMB_mats = list_subset(COMB_mats, select_vector = !null_remove))
  if (!is.null(extra)) out$extra = extra[!null_remove, , drop = FALSE]
  out
}

#Main function: compute the filter_params table for a set of cgpVAF matrices.
#  sex: "male" or "female" (from the donor config)
#  mixed_colony_vaf_cutoff: if not NULL, colonies whose peak VAF is below this are removed first ("upfront
#    colony filtering", the old -m/-v options)
#  vaf_plot_file: with upfront colony filtering, a PDF of VAF density plots (removed colonies first, then kept)
#Returns list(COMB_mats, filter_params, peak_vaf, mixed_samples)
compute_filter_params = function(COMB_mats,
                                 sex,
                                 min_sample_mean_cov = 0,
                                 exclude_samples = NULL,
                                 mixed_colony_vaf_cutoff = NULL,
                                 vaf_plot_file = NULL) {
  if (!sex %in% c("male", "female")) stop("sex must be 'male' or 'female', not '", sex, "'")
  COMB_mats$gender = sex

  if (min_sample_mean_cov > 0 || length(exclude_samples) > 0) {
    COMB_mats = remove_low_coverage_samples(COMB_mats = COMB_mats,
                                            min_sample_mean_cov = min_sample_mean_cov,
                                            other_samples_to_remove = exclude_samples,
                                            min_variant_reads_auto = 3,
                                            min_variant_reads_xy = 2)
  }
  COMB_mats$nsamp = ncol(COMB_mats$NV)
  COMB_mats = remove_null_mutations(COMB_mats)$COMB_mats

  peak_vaf = NULL
  mixed_samples = character(0)
  if (!is.null(mixed_colony_vaf_cutoff)) {
    #Germline and beta-binomial filters first, to recognise and remove mixed colonies.
    #As in the original script, these two are NOT recalculated after the mixed colonies are removed.
    basic_params = data.frame(germline_pval = germline.binomial.filter(COMB_mats),
                              bb_rhoval = beta.binom.filter(COMB_mats))
    message("Testing for mixed colonies: removing samples with a peak VAF < ", mixed_colony_vaf_cutoff)
    peak_vaf = sapply(colnames(COMB_mats$NV), check_peak_vaf, COMB_mats = COMB_mats, filter_params = basic_params)
    mixed_samples = names(peak_vaf)[peak_vaf < mixed_colony_vaf_cutoff]
    for (s in mixed_samples) message("Removing sample ", s, " as peak VAF is ", round(peak_vaf[s], 2), ", suggesting a non-clonal sample")
    if (!is.null(vaf_plot_file)) {
      pdf(vaf_plot_file)
      for (s in c(mixed_samples, setdiff(names(peak_vaf), mixed_samples))) vaf_density_plot(s, COMB_mats, filter_params = basic_params)
      dev.off()
    }

    #Only drop columns when there is something to drop (the original `[, -which(...)]` removed every sample when none failed)
    if (length(mixed_samples) > 0) {
      keep = !colnames(COMB_mats$NV) %in% mixed_samples
      COMB_mats$NV = COMB_mats$NV[, keep, drop = FALSE]
      COMB_mats$NR = COMB_mats$NR[, keep, drop = FALSE]
    }
    nulls = remove_null_mutations(COMB_mats, extra = basic_params)
    COMB_mats = nulls$COMB_mats
    basic_params = nulls$extra
  } else {
    basic_params = data.frame(germline_pval = germline.binomial.filter(COMB_mats),
                              bb_rhoval = beta.binom.filter(COMB_mats))
  }

  COMB_mats$PVal = pval_matrix(COMB_mats)
  filter_params = data.frame(basic_params,
                             mean_depth = get_mean_depth(COMB_mats),
                             pval_within_pos_dp2 = low_vaf_in_pos_samples_dp2(COMB_mats),
                             pval_within_pos_dp3 = low_vaf_in_pos_samples_dp3(COMB_mats),
                             max_depth_in_pos_samples = get_max_depth_in_pos(COMB_mats))
  rownames(filter_params) = COMB_mats$mat$mut_ref
  filter_params$max_pval_in_pos_sample = get_max_pval_in_pos(COMB_mats)
  filter_params$max_mut_vaf = get_max_vaf(COMB_mats)

  list(COMB_mats = COMB_mats, filter_params = filter_params, peak_vaf = peak_vaf, mixed_samples = mixed_samples)
}

#The "reduced" set drops mutations that are near-certainly germline (fail both the binomial and beta-binomial
#tests at very stringent cut-offs), to make file sizes manageable
reduce_mutation_set = function(COMB_mats, filter_params) {
  reduced = filter_params$germline_pval < 0.01 | filter_params$bb_rhoval > 0.01
  list(COMB_mats = list_subset(COMB_mats, reduced), filter_params = filter_params[reduced, ])
}
