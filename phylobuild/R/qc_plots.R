# QC plots for the tree-building step.

#Heatmap of mutations that fit the tree poorly: one row per mutation, one column per sample (in tree tip order,
#under the tree). Colour = VAF (grey = no reads). A dot marks the samples that the branch treemut assigned the
#mutation to says should carry it, so mismatches between the read counts and the assignment are visible.
plot_poor_fit_heatmap = function(tree, COMB_mats.tree.build, mut_refs, file, max_muts = 300) {
  pdf(file, width = 14, height = 10)
  on.exit(dev.off())
  if (length(mut_refs) == 0) {
    plot.new(); text(0.5, 0.5, "No mutations fit the tree poorly")
    return(invisible())
  }
  details = COMB_mats.tree.build$mat
  rows = match(mut_refs, details$mut_ref)
  rows = rows[order(details$node[rows], details$Chrom[rows], details$Pos[rows])] #group by assigned branch
  shown = utils::head(rows, max_muts)

  plot_tree = if ("Ancestral" %in% tree$tip.label) drop.tip(tree, "Ancestral") else tree
  n = length(plot_tree$tip.label)
  layout(matrix(1:2), heights = c(1, 2))
  par(mar = c(0, 12, 3, 1))
  plot(plot_tree, direction = "downwards", show.tip.label = FALSE,
       main = sprintf("%d mutations fit the tree poorly%s", length(rows), if (length(rows) > max_muts) sprintf(" (first %d shown)", max_muts) else ""))
  xx = get("last_plot.phylo", envir = .PlotPhyloEnv)$xx[seq_len(n)]
  samples = plot_tree$tip.label[order(xx)]
  usr_x = par("usr")[1:2] #the heatmap uses the tree's x coordinates, so columns sit under their tips
  xs = sort(xx)

  NV = as.matrix(COMB_mats.tree.build$NV[shown, samples, drop = FALSE])
  NR = as.matrix(COMB_mats.tree.build$NR[shown, samples, drop = FALSE])
  vaf = ifelse(NR > 0, NV / NR, NA)
  #Samples expected to carry each mutation, from its assigned branch (in the tree it was assigned on)
  expected = t(sapply(details$node[shown], function(node) samples %in% getTips(tree, node)))
  if (length(shown) == 1) expected = matrix(expected, nrow = 1)

  par(mar = c(6, 12, 0, 1))
  m = length(shown)
  plot(NA, xlim = usr_x, ylim = c(0.5, m + 0.5), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  pal = grDevices::colorRampPalette(c("white", "red3"))(50)
  cell = cbind(col = rep(seq_len(n), each = m), row = rep(seq_len(m), n))
  v = vaf[cbind(cell[, "row"], cell[, "col"])]
  fill = ifelse(is.na(v), "grey85", pal[pmax(1, ceiling(v * 50))])
  half = if (n > 1) min(diff(xs)) / 2 else 0.5
  rect(xs[cell[, "col"]] - half, cell[, "row"] - 0.5, xs[cell[, "col"]] + half, cell[, "row"] + 0.5, col = fill, border = NA)
  dots = which(expected, arr.ind = TRUE)
  if (nrow(dots) > 0) points(xs[dots[, 2]], dots[, 1], pch = 16, cex = 0.3)
  axis(1, at = xs, labels = samples, las = 2, cex.axis = 0.5, tick = FALSE)
  labels = paste0(substr(details$mut_ref[shown], 1, 28), ifelse(nchar(details$mut_ref[shown]) > 28, "...", ""), " (", details$Mut_type[shown], ")")
  if (m <= 120) axis(2, at = seq_len(m), labels = labels, las = 1, cex.axis = 0.45, tick = FALSE)
  legend("topright", inset = c(0, -0.04), xpd = TRUE, horiz = TRUE, bty = "n", cex = 0.7,
         legend = c("VAF 0", "VAF 1", "no reads", "expected from assigned branch"),
         pch = c(22, 22, 22, 16), pt.bg = c("white", "red3", "grey85", NA), col = "black", pt.cex = c(1.5, 1.5, 1.5, 0.6))
  invisible()
}
