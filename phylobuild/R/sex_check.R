# Checking the donor's sex from the data.
#
# The old rule (Y calls / X calls > 1%) labelled 7 of 9 female donors as male (checked 2026-10-30).
# These measures separated the sexes cleanly in all 14 donors tested:
#   x_auto_depth  median depth at non-PAR X sites / autosomal sites   male 0.49-0.56, female 0.92-0.98
#   x_het_frac    fraction of germline-like non-PAR X variants that are heterozygous   male <=0.31, female >=0.83
#   y_auto_depth  median depth at Y sites / autosomal sites   male >=0.33, female <=0.04
# None of them use the germline p-values, which are themselves calculated using the assumed sex.

sex_check_metrics = function(COMB_mats, build = guess_build(COMB_mats$mat$Chrom)) {
  NV = as.matrix(COMB_mats$NV); NR = as.matrix(COMB_mats$NR)
  chrom = COMB_mats$mat$Chrom
  auto = is_autosome(chrom)
  x_nonpar = is_x(chrom) & !in_par(chrom, COMB_mats$mat$Pos, build)
  y = is_y(chrom)
  mean_depth = rowMeans(NR)

  #Germline-like: seen in >80% of adequately covered colonies
  covered = NR >= 5
  frac_pos = rowSums(NV >= 1 & covered) / pmax(rowSums(covered), 1)
  germ_x = frac_pos > 0.8 & rowSums(covered) >= 5 & x_nonpar
  pooled_vaf = rowSums(NV) / pmax(rowSums(NR), 1)

  list(
    summary = data.frame(
      n_x_sites = sum(x_nonpar), n_y_sites = sum(y), n_germline_x = sum(germ_x),
      x_auto_depth = median(mean_depth[x_nonpar]) / median(mean_depth[auto]),
      x_het_frac = if (sum(germ_x) > 0) mean(pooled_vaf[germ_x] > 0.3 & pooled_vaf[germ_x] < 0.7) else NA_real_,
      y_auto_depth = if (sum(y) > 0) median(mean_depth[y]) / median(mean_depth[auto]) else 0
    ),
    #Per colony: flags X loss, or colonies of the other sex (maternal cells, transplant recipient, contamination)
    per_colony = data.frame(
      sample = colnames(NR),
      mean_depth = colMeans(NR[auto, , drop = FALSE]),
      x_auto_depth = apply(NR, 2, function(d) median(d[x_nonpar]) / median(d[auto]))
    )
  )
}

#Call sex from each measure; "unclear" if a measure falls in the gap between the sexes
infer_sex = function(summary) {
  call = function(value, female_if, male_if) ifelse(is.na(value), "unclear", ifelse(female_if(value), "female", ifelse(male_if(value), "male", "unclear")))
  c(x_auto_depth = call(summary$x_auto_depth, function(v) v > 0.8, function(v) v < 0.65),
    x_het_frac = call(summary$x_het_frac, function(v) v > 0.6, function(v) v < 0.4),
    y_auto_depth = call(summary$y_auto_depth, function(v) v < 0.15, function(v) v > 0.25))
}

#Compare the configured sex with the data and decide which sex the pipeline uses (result$sex).
#  configured_sex: "male", "female", or "unknown" to use the inferred sex
#  strict: stop on a mismatch (otherwise warn and use the configured sex)
check_donor_sex = function(COMB_mats, configured_sex, build = guess_build(COMB_mats$mat$Chrom), strict = TRUE) {
  if (!configured_sex %in% c("male", "female", "unknown")) stop("sex must be 'male', 'female' or 'unknown', not '", configured_sex, "'")
  m = sex_check_metrics(COMB_mats, build)
  calls = infer_sex(m$summary)
  decided = unique(calls[calls != "unclear"])
  #Every measure that gives a call must agree, and the depth ratio (the most robust) must give one
  inferred = if (length(decided) == 1 && calls["x_auto_depth"] != "unclear") decided else "unclear"

  #Colonies that look like the other sex (or have lost an X)
  r = m$per_colony$x_auto_depth
  m$per_colony$flag = if (inferred == "female") ifelse(r < 0.75, "low X: possible X loss or male cells", "")
                      else if (inferred == "male") ifelse(r > 0.75, "high X: possible female cells", "")
                      else ""

  msg = sprintf("Configured sex: %s. Inferred from data: %s (X:autosome depth %.2f -> %s; X het fraction %.2f -> %s; Y:autosome depth %.2f -> %s)",
                configured_sex, inferred,
                m$summary$x_auto_depth, calls["x_auto_depth"], m$summary$x_het_frac, calls["x_het_frac"], m$summary$y_auto_depth, calls["y_auto_depth"])
  if (configured_sex == "unknown") {
    if (inferred == "unclear") stop("Sex is 'unknown' in the config and can't be inferred reliably from the data. ", msg,
                                    ". Set the donor's sex in the config. (A donor with colonies of both sexes, e.g. a sex-mismatched transplant, needs the per-colony plot.)")
    ok = TRUE
    sex = inferred
  } else {
    ok = configured_sex == inferred
    if (!ok && strict) stop("Sex check failed. ", msg, ". Fix the donor's 'sex' in the config, set it to 'unknown' to use the inferred sex, or set sex_check: warn to continue anyway.")
    sex = configured_sex
  }
  if (!ok) warning(msg) else message(msg, ". Using: ", sex)
  c(m, list(calls = calls, inferred = inferred, configured = configured_sex, sex = sex, ok = ok, message = msg))
}

plot_sex_check = function(result, file) {
  pdf(file, width = 10, height = 5)
  pc = result$per_colony
  plot(seq_len(nrow(pc)), pc$x_auto_depth, ylim = c(0, 1.3), pch = 19, col = ifelse(pc$flag == "", "grey40", "red"),
       xlab = "Colony", ylab = "X : autosome depth", main = paste("Per-colony X dosage -", result$message), cex.main = 0.6, xaxt = "n")
  abline(h = c(0.5, 1), lty = 2, col = "grey")
  flagged = which(pc$flag != "")
  if (length(flagged) > 0) text(flagged, pc$x_auto_depth[flagged], pc$sample[flagged], pos = 3, cex = 0.6, col = "red")
  dev.off()
}
