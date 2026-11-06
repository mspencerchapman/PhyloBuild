# Running MPBoot and assigning mutations to tree branches with treemut
# (ported from foetal.filters.parallel.R and filtering_from_table_mix_remove_GRCh38.R).

#Run MPBoot on a fasta file of dummy DNA strings and return the tree.
#A fixed seed makes the tree reproducible (MPBoot otherwise uses a random seed).
run_mpboot = function(fasta, mpboot = "mpboot", bootstraps = 1000, seed = NULL) {
  #remove outputs of any earlier run on this file, so a failed run can't leave an old tree behind
  unlink(Sys.glob(paste0(fasta, ".*")))
  args = c("-s", fasta, "-bb", bootstraps)
  if (!is.null(seed)) args = c(args, "-seed", seed)
  message("Running: ", mpboot, " ", paste(args, collapse = " "))
  status = system2(mpboot, args)
  treefile = paste0(fasta, ".treefile")
  if (status != 0 || !file.exists(treefile)) stop("MPBoot failed (exit status ", status, ") - see the log above")
  read.tree(treefile)
}

#Add the ancestral tip at the end (keeps tip numbering of the existing samples)
add_ancestral_outgroup = function(tree, outgroup_name = "Ancestral") {
  tmp = tree$edge
  N = length(tree$tip.label)
  newroot = N + 2
  renamedroot = N + 3
  ancestral_tip = N + 1
  tmp = ifelse(tmp > N, tmp + 2, tmp)

  tree$edge = rbind(c(newroot, renamedroot), tmp, c(newroot, ancestral_tip))
  tree$edge.length = c(0, tree$edge.length, 0)

  tree$tip.label = c(tree$tip.label, outgroup_name)
  tree$Nnode = tree$Nnode + 1
  mode(tree$Nnode) = "integer"
  mode(tree$edge) = "integer"
  tree
}

get_ancestral_nodes = function(node, edge, exclude_root = TRUE) {
  idx = which(edge[, 2] == node)
  parents = node
  while (length(idx) > 0) {
    if (length(idx) > 1) stop("multiple parents!")
    parent = edge[idx, 1]
    parents = c(parents, parent)
    idx = which(edge[, 2] == parent)
  }
  if (exclude_root) parents[-length(parents)] else parents
}

#Assign mutations to the branches of the tree described by df (from reconstruct_genotype_summary).
#Uses treemut's assign_to_df (nangalialab/treemut >= 1.1), the function the old scripts called as
#assign_to_tree(mtr, dep, df, ...) from a copied treemut.R. On a 47-colony test donor it puts 99.98% of mutations
#on the same branch as that copy; the 2 that differ were near 50:50 between two branches. treemut's exported
#assign_to_tree(tree, ...) is not used because it ignores the per-sample error rates, which the pipeline needs
#for the Ancestral pseudo-sample.
assign_mutations = function(mtr, depth, df, error_rate, maxits = 5) {
  ns = asNamespace("treemut")
  if (!exists("assign_to_df", envir = ns, inherits = FALSE)) {
    stop("This version of treemut has no assign_to_df(). Install the current version with remotes::install_github('nangalialab/treemut')")
  }
  f = get("assign_to_df", envir = ns)
  args = list(mtr, depth, df, error_rate = error_rate, maxits = maxits)
  if ("bverbose" %in% names(formals(f))) args$bverbose = FALSE
  do.call(f, args)
}

#Read count matrices in the column order treemut expects, with a dummy "Ancestral" sample if needed
treemut_inputs = function(COMB_mats.tree.build, df, keep_ancestral) {
  mtr = COMB_mats.tree.build$NV
  depth = COMB_mats.tree.build$NR
  if (keep_ancestral) {mtr$Ancestral = 0; depth$Ancestral = 10}
  list(mtr = as.matrix(mtr)[, df$samples, drop = FALSE], depth = as.matrix(depth)[, df$samples, drop = FALSE])
}

assign_mutations_to_branches = function(tree, filtered_muts, keep_ancestral = TRUE, create_multi_tree = TRUE, p.error.value = 0.01, treefit_pval_cutoff = 1e-3) {
  tree = drop.tip(tree, "Ancestral")
  tree = multi2di(tree)
  if (keep_ancestral) {
    message("Assigning mutations with an ancestral branch")
    tree = add_ancestral_outgroup(tree) #Re-add the ancestral outgroup after making tree dichotomous - avoids the random way that baseline polytomy is resolved
  } else {
    message("Assigning mutations without an ancestral branch")
  }
  tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for treemut to work

  df = reconstruct_genotype_summary(tree)
  m = treemut_inputs(filtered_muts$COMB_mats.tree.build, df, keep_ancestral)
  p.error = sapply(df$samples, function(x) ifelse(x == "Ancestral", 1e-6, p.error.value))
  res = assign_mutations(m$mtr, m$depth, df, error_rate = p.error)

  if (create_multi_tree) {
    message("Converting to a multi-furcating tree structure")
    tree$edge.length = res$df$df$edge_length
    #Maintain the dichotomy with the ancestral branch
    if (keep_ancestral) {
      ROOT = tree$edge[1, 1]
      not_anc = tree$edge[, 1] == ROOT & tree$edge[, 2] != which(tree$tip.label == "Ancestral")
      tree$edge.length[not_anc] = ifelse(tree$edge.length[not_anc] == 0, 1, tree$edge.length[not_anc])
    }
    tree = di2multi(tree)
    df = reconstruct_genotype_summary(tree)
    m = treemut_inputs(filtered_muts$COMB_mats.tree.build, df, keep_ancestral)
    res = assign_mutations(m$mtr, m$depth, df, error_rate = p.error)
  }

  tree$edge.length = res$df$df$edge_length
  res$tree = tree

  poor_fit = res$summary$pval < treefit_pval_cutoff
  message(sum(poor_fit), " mutations do not have read counts that fit any tree branch well")
  return(res)
}

#Store the branch assignment of each mutation in the filtered_muts info table
store_assignment = function(filtered_muts, res) {
  filtered_muts$COMB_mats.tree.build$mat$node = res$tree$edge[res$summary$edge_ml, 2]
  filtered_muts$COMB_mats.tree.build$mat$pval = res$summary$pval
  filtered_muts
}

#Edge lengths = number of mutations assigned to each branch
edge_lengths_from_assignment = function(tree, details) {
  tree$edge.length = sapply(tree$edge[, 2], function(node) sum(details$node == node))
  tree
}

#Look for mutations removed by the germline filter that are convincingly absent from a small clade near the
#root - these are early somatic mutations present in most, but not all, samples.
check_for_false_germline_calls = function(tree,
                                          COMB_mats,
                                          filter_params,
                                          max_clade_prop = 0.1, #only test root clades containing fewer than this proportion of samples
                                          SNVs_only = TRUE, #Only re-add SNVs (indels are likely to be high frequency artefacts)
                                          CN_table = NULL,
                                          germline_log10_cutoff = -10,
                                          candidate_filter = NULL) { #optional function(mut_refs, clade_samples) -> logical keep
  get_root_clades = function(tree) {
    tree = di2multi(tree)
    ROOT = tree$edge[1, 1]
    clades = tree$edge[tree$edge[, 1] == ROOT, 2]
    lapply(clades, function(node) getTips(tree, node))
  }

  tree.noancestral = drop.tip(tree, "Ancestral")
  root_clades = get_root_clades(tree.noancestral)
  nsamp = length(tree.noancestral$tip.label)
  root_clades = root_clades[sapply(root_clades, length) < max_clade_prop * nsamp]

  if (length(root_clades) == 0) {
    message("There are no clades from the root that include <", max_clade_prop * 100, "% of samples")
    return(NULL)
  }
  rownames(COMB_mats$NR) = rownames(COMB_mats$NV) = rownames(COMB_mats$PVal) = COMB_mats$mat$mut_ref
  colnames(COMB_mats$NR) = colnames(COMB_mats$PVal) = colnames(COMB_mats$NV) = gsub("_MTR", "", colnames(COMB_mats$NV))

  res = lapply(root_clades, function(outlier_sample_group) {
    message("Testing outlier group: ", paste(outlier_sample_group, collapse = " "))
    germline_filtered = rownames(filter_params)[log10(filter_params$germline_pval) > germline_log10_cutoff]

    #Exclude those that are subject to a deletion/ LOH in any of the outlier sample group
    if (!is.null(CN_table)) {
      CN_restricted = CN_table[CN_table$Sample %in% outlier_sample_group & CN_table$Type %in% c("LOH", "DEL"), ]
      for (k in seq_len(nrow(CN_restricted))) {
        muts_df = split_mut_ref(germline_filtered)
        exclude = strip_chr(muts_df$Chrom) == strip_chr(CN_restricted$Chrom[k]) & muts_df$Pos > CN_restricted$Pos_min[k] & muts_df$Pos < CN_restricted$Pos_max[k]
        germline_filtered = germline_filtered[!exclude]
      }
    }

    NR_outlier = rowSums(COMB_mats$NR[germline_filtered, outlier_sample_group, drop = FALSE])
    NV_outlier = rowSums(COMB_mats$NV[germline_filtered, outlier_sample_group, drop = FALSE])
    outlier_pvals = mapply(FUN = function(NV, NR) {if (NR == 0) 1 else binom.test(NV, NR, alternative = "less")$p.value}, NV = NV_outlier, NR = NR_outlier)
    outlier_pval.adj = p.adjust(outlier_pvals, method = "BH")

    mut_refs = germline_filtered[outlier_pval.adj < 0.05 & NV_outlier == 0]
    if (SNVs_only && length(mut_refs) > 0) mut_refs = mut_refs[is.snv(mut_refs)]
    if (!is.null(candidate_filter) && length(mut_refs) > 0) {
      keep = candidate_filter(mut_refs, outlier_sample_group)
      message(sum(!keep), " of ", length(mut_refs), " candidates fail the depth/overdispersion/clonal VAF checks and are not added back")
      mut_refs = mut_refs[keep]
    }
    if (length(mut_refs) > 0) {
      message(length(mut_refs), " germline-filtered mutations are convincingly absent in this group")
      return(mut_refs)
    }
    message("No mutations called as germline are robustly absent in this group (this relies on adequate coverage)")
    NULL
  })
  return(res)
}
