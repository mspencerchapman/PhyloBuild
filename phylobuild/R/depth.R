# Mean-depth window for keeping a site (the tree step's "mean_depth" filter).
# The filter removes outlying sites: very low mean depth across colonies means poor mapping, very high mean depth
# means mismapping (repeats, collapsed duplications). "auto" treats this as outlier detection: the window is
# median +/- k * MAD of the sites' mean depth, so it follows both the donor's coverage and how variable the
# per-site mean is (which depends on the number of colonies). k = 4 gives 6.8-23.5 for a 47-colony test donor
# (~15x), close to the fixed 6-23 used before.
# X/Y in males get their own window from the X/Y sites (falling back to half the autosomal values if there are
# too few X/Y sites to estimate the spread).

depth_window = function(d, k) {
  c(low = max(0, median(d) - k * mad(d)), high = median(d) + k * mad(d))
}

depth_cutoffs_auto = function(filter_params, COMB_mats, k = 4, min_xy_sites = 200) {
  xy = haploid_rows(COMB_mats)
  auto_depth = filter_params$mean_depth[is_autosome(COMB_mats$mat$Chrom)]
  auto = depth_window(auto_depth, k)
  if (!any(xy)) {
    xy_win = auto
  } else if (sum(xy) >= min_xy_sites) {
    xy_win = depth_window(filter_params$mean_depth[xy], k)
  } else {
    xy_win = auto / 2
  }
  c(auto_low = unname(auto["low"]), auto_high = unname(auto["high"]),
    xy_low = unname(xy_win["low"]), xy_high = unname(xy_win["high"]),
    median_depth = median(auto_depth), mad_depth = mad(auto_depth))
}

plot_depth_cutoffs = function(filter_params, COMB_mats, cutoffs, file) {
  xy = haploid_rows(COMB_mats)
  pdf(file, width = 10, height = 5)
  par(mfrow = c(1, 2))
  for (set in c("autosomes", "X/Y")) {
    d = filter_params$mean_depth[if (set == "autosomes") !xy else xy]
    lims = if (set == "autosomes") cutoffs[c("auto_low", "auto_high")] else cutoffs[c("xy_low", "xy_high")]
    if (length(d) == 0) {plot.new(); text(0.5, 0.5, paste("No", set, "sites treated separately")); next}
    hist(d[d < 3 * cutoffs["median_depth"]], breaks = 100, main = paste("Mean depth across colonies:", set),
         xlab = "Mean depth at site", sub = sprintf("window %.1f-%.1f keeps %.1f%% of sites", lims[1], lims[2], 100 * mean(d >= lims[1] & d <= lims[2])))
    abline(v = lims, col = "red", lwd = 2)
  }
  dev.off()
}
