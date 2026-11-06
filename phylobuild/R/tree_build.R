# The tree-building step (ported from filtering_from_table_mix_remove_GRCh38.R).
#
# 1. Filter the mutations (p-value or VAF based) and genotype each sample
# 2. Build a tree with MPBoot and assign mutations to branches with treemut
# 3. Remove poorly fitting indels; add back early mutations wrongly removed by the germline filter; rebuild
# 4. Annotate mutations (VAGrENT, optional)
# 5. Remove non-clonal ("mixed") samples by testing their VAFs against the tree; rebuild
# 6. Correct branch lengths for per-sample sensitivity
# 7. Remove duplicate samples, then plot
#
# Output files keep the names used by the original script (objects as .rds), under output_dir/{filtered_muts,mpboot_files,
# tree_files,annotated_muts,mutation_vcfs,plots}.

#Thresholds for the two filtering strategies
filtering_thresholds = function(filtering_type) {
  switch(filtering_type,
         pval = list(min_pval_for_true_somatic = 0.1, min_pval_for_true_somatic_SHARED = 0.05, min_vaf = NA, min_vaf_SHARED = NA),
         vaf = list(min_pval_for_true_somatic = NA, min_pval_for_true_somatic_SHARED = NA, min_vaf = c(0.2, 0.4), min_vaf_SHARED = c(0.15, 0.30)),
         stop("Unrecognised filtering type '", filtering_type, "' - use 'pval' or 'vaf'"))
}

#Genotypes for the shared mutations, adjusted as requested, and the DNA strings MPBoot uses.
#Applied the same way every time the tree is (re)built.
prepare_tree_input = function(Genotype_shared_bin, filter_y = FALSE, CN_table = NULL, snvs_only = FALSE, build = "hg38") {
  muts = split_mut_ref(rownames(Genotype_shared_bin))
  if (filter_y) {
    message("Removing variants on the Y chromosome and X chromosome PARs from the dummy dna strings used for tree building")
    keep = !(is_y(muts$Chrom) | in_par(muts$Chrom, muts$Pos, build))
    Genotype_shared_bin = Genotype_shared_bin[keep, , drop = FALSE]
    muts = muts[keep, ]
  }
  #Samples with an LOH or deletion over a region can't be called negative there, so set to "uncertain".
  #Rows for samples no longer in the set (low coverage, excluded, non-clonal) are skipped.
  if (!is.null(CN_table)) CN_table = CN_table[CN_table$Type %in% c("LOH", "DEL") & CN_table$Sample %in% colnames(Genotype_shared_bin), , drop = FALSE]
  for (k in seq_len(if (is.null(CN_table)) 0 else nrow(CN_table))) {
    rows = strip_chr(muts$Chrom) == strip_chr(CN_table$Chrom[k]) & muts$Pos > CN_table$Pos_min[k] & muts$Pos < CN_table$Pos_max[k]
    if (sum(rows) == 0) warning("CN table row for ", CN_table$Sample[k], " (", CN_table$Chrom[k], ":", CN_table$Pos_min[k], "-", CN_table$Pos_max[k], ") matches no shared mutations - check the chromosome names and coordinates")
    message("Altering genotypes at ", sum(rows), " locations in ", CN_table$Sample[k])
    Genotype_shared_bin[rows, CN_table$Sample[k]] = pmax(0.5, Genotype_shared_bin[rows, CN_table$Sample[k]])
  }
  list(Genotype_shared_bin = Genotype_shared_bin,
       dna_strings = dna_strings_from_genotype(if (snvs_only) Genotype_shared_bin[is.snv(rownames(Genotype_shared_bin)), , drop = FALSE] else Genotype_shared_bin))
}

#Test each sample's branch VAFs for consistency with a clonal sample. Returns one row per sample with result PASS/FAIL.
#  exclude_rows: optional function(sample) giving a logical over the rows of details to leave out of that sample's
#  test (e.g. copy-number altered regions, where a clonal sample's VAF isn't 0.5)
test_sample_clonality = function(tree, details, matrices, keep_ancestral, min.mean.vaf = 0.425, exclude_rows = NULL) {
  sample_test = lapply(setdiff(tree$tip.label, "Ancestral"), function(samples, min.depth = 1) {
    excluded = if (is.null(exclude_rows)) rep(FALSE, nrow(details)) else exclude_rows(samples)
    node_test = sapply(tree$edge[, 2], function(node) {
      info = get_edge_info(tree, details, node)
      idx = info$idx.in.details[!excluded[info$idx.in.details]]
      if (length(idx) == 0) return(NA)
      df = data.frame(mtr = sum(matrices$mtr[idx, samples], na.rm = TRUE),
                      dep = sum(matrices$dep[idx, samples], na.rm = TRUE), stringsAsFactors = FALSE)
      df = df[which(df$dep >= min.depth), ]
      df$vaf = df$mtr / df$dep
      df = df[which(!is.na(df$vaf)), ]
      MTR = sum(df$mtr)
      DEP = sum(df$dep)
      if (DEP == 0) return(NA)
      z1 = binom.test(MTR, DEP, alternative = "less", p = min.mean.vaf)
      z2 = binom.test(MTR, DEP, alternative = "greater", p = 0.05)
      comb.p.value = max(z1$p.value, z2$p.value)
      if (comb.p.value < 0.05) {
        if (comb.p.value < 0.05 / dim(tree$edge)[1]) 2 else 1
      } else if (z1$p.value > 0.05 & z2$p.value < 0.05) {
        3
      } else {
        0
      }
    })
    names(node_test) = tree$edge[, 2]
    node_test = node_test[!is.na(node_test)]

    #The 'consistent branch' test - checks that all the "POSITIVE" branches are those expected to be positive
    consistent_branches = as.numeric(names(node_test[node_test == 3]))
    consistent_branch_lengths = sapply(consistent_branches, function(node) sum(details$node == node))
    consistent_branch_test = all(consistent_branches[consistent_branch_lengths > 3] %in% get_ancestral_nodes(which(tree$tip.label == samples), tree$edge))

    #The terminal branch is allowed a few subclonal (potentially in vitro) mutations, so is ignored. If there is an
    #ancestral branch, ignore the root branches too: mutations there are likely to be artefactual.
    if (keep_ancestral) {
      root_branches = tree$edge[which(tree$edge[, 1] == tree$edge[1, 1]), 2]
      branches_to_ignore = c(root_branches, which(tree$tip.label == samples))
    } else {
      branches_to_ignore = which(tree$tip.label == samples)
    }
    data.frame(sample = samples,
               high_prob = sum(node_test == 2 & !as.numeric(names(node_test)) %in% branches_to_ignore),
               low_prob = sum(node_test == 1 & !as.numeric(names(node_test)) %in% branches_to_ignore),
               consistent_pos_branches = consistent_branch_test)
  })
  sample_test_df = Reduce(rbind, sample_test)
  sample_test_df$result = ifelse(sample_test_df$high_prob > 0 | sample_test_df$low_prob > 2 | !sample_test_df$consistent_pos_branches, "FAIL", "PASS")
  sample_test_df
}

#Remove samples and any mutations then no longer positive in any sample
drop_samples_from_mut_set = function(COMB_mats.tree.build, drop) {
  for (m in c("Genotype_bin", "NV", "NR", "PVal")) {
    COMB_mats.tree.build[[m]] = COMB_mats.tree.build[[m]][, !colnames(COMB_mats.tree.build[[m]]) %in% drop, drop = FALSE]
  }
  still_positive = rowSums(COMB_mats.tree.build$Genotype_bin == 1) > 0
  list_subset(COMB_mats.tree.build, select_vector = still_positive)
}

#Mean mutation burden of the real samples (the Ancestral tip, at height 0, is left out)
mean_sample_burden = function(tree) {
  heights = nodeHeights(tree)[match(seq_along(tree$tip.label), tree$edge[, 2]), 2]
  mean(heights[tree$tip.label != "Ancestral"])
}

#Find sets of samples that are near-identical (fewer private mutations than the threshold) and pick which to drop.
#Threshold: max(min_private, prop x mean burden)
find_duplicate_samples = function(tree, details, min_private = 30, prop = 0.03) {
  mean_mut_burden = mean_sample_burden(tree)
  private_mut_threshold = max(min_private, prop * mean_mut_burden)
  message("Threshold for determining duplicate samples is ", private_mut_threshold, " mutations.")
  internal = tree$edge[, 2][!tree$edge[, 2] %in% 1:length(tree$tip.label)]
  pseudo_terminal_nodes = sapply(internal, function(node) {
    node_height = nodeheight(tree = tree, node = node)
    samples = get_edge_info(tree, details, node)$samples
    sample_heights = nodeHeights(tree)[tree$edge[, 2] %in% which(tree$tip.label %in% samples), 2]
    if (all((sample_heights - node_height) < private_mut_threshold)) node else NA
  })
  pseudo_terminal_nodes = pseudo_terminal_nodes[!is.na(pseudo_terminal_nodes)]
  duplicate_sets = lapply(pseudo_terminal_nodes, function(node) getTips(tree, node = node))
  for (x in duplicate_sets) message(paste(x, collapse = " & "), " are recognised as duplicate samples")

  #Keep the sample with the shortest branch from each set
  drop = unlist(lapply(duplicate_sets, function(samples) {
    included = samples[samples %in% tree$tip.label]
    if (length(included) < 2) return(NULL)
    heights = sapply(included, function(s) nodeheight(tree = tree, node = which(tree$tip.label == s)))
    retain = included[heights == min(heights)][1]
    included[included != retain]
  }))
  unique(drop)
}

#Mark coding changes, and coding changes in genes of interest, for plotting
annotate_coding_changes = function(mat, genes_to_annotate = NULL) {
  coding_types = c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                   "protein_coding:exon:CDS:insertion:frameshift_variant",
                   "protein_coding:exon:CDS:deletion:frameshift_variant",
                   "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                   "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                   "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss")
  is_coding = mat$Type %in% coding_types | grepl("splice_site_variant", mat$Type)
  mat$coding_change = ifelse(is_coding, "Coding change", "no")
  mat$coding_change_chip = ifelse(mat$Gene %in% genes_to_annotate & is_coding, "Coding change mutation in driver", "no")
  mat$ChromPos = paste(mat$Chrom, mat$Pos, sep = "-")
  mat$variant_ID = paste(mat$Gene, mat$Protein, sep = " ")
  mat
}

#Read the sensitivity table and add the dummy "Ancestral" sample; zero sensitivities become 0.001
read_sensitivity = function(path) {
  df = read.delim(path, header = TRUE, stringsAsFactors = FALSE, sep = " ")
  df = rbind(df, data.frame(Sample = "Ancestral", SNV_sensitivity = 1, INDEL_sensitivity = 1))
  df$SNV_sensitivity[df$SNV_sensitivity == 0] = 0.001
  df$INDEL_sensitivity[df$INDEL_sensitivity == 0] = 0.001
  df
}

build_tree = function(COMB_mats,
                      filter_params,
                      run_id,
                      filtering_type = "pval",
                      output_dir,
                      sensitivity_path,
                      donor_id,
                      recipient_id = "none",
                      build = "hg38",
                      min_sample_mean_cov = 0,
                      exclude_samples = NULL,
                      exclude_muts = NULL,
                      keep_muts = NULL,
                      polytomous = TRUE,
                      bb_cutoff = 0.1,
                      germline_log10_cutoff = -10,
                      depth_cutoffs = c(auto_low = 6, auto_high = 23, xy_low = 3, xy_high = 10),
                      filter_y = FALSE,
                      keep_ancestral = TRUE,
                      remove_nonclonal = TRUE,
                      germline_addback = TRUE,
                      snvs_only = FALSE,
                      CN_table = NULL,
                      genes_to_annotate = NULL,
                      donor_age = NA,
                      annotate = NULL, #function(mat) returning Gene/Transcript/RNA/CDS/Protein/Type/SO_codes columns, or NULL
                      mpboot = "mpboot",
                      mpboot_bootstraps = 1000,
                      seed = NULL,
                      treefit_pval_cutoff = 1e-3, #treemut p-value below which a mutation "fits the tree poorly"
                      remove_poor_fit = FALSE, #drop poorly fitting mutations (all types) before rebuilding
                      addback_require_filters = FALSE, #added-back germline-filtered mutations must pass depth/rho checks and look clonal elsewhere
                      clonality_exclude_cn_xy = FALSE, #clonality test ignores CN-table regions and male haploid sites
                      duplicate_min_private = 30, #duplicate samples: fewer private mutations than max(this, ...)
                      duplicate_prop = 0.03, #... and this proportion of the mean burden
                      legacy_rdata = FALSE) { #also write the old save() files alongside the .rds
  save_muts = function(filtered_muts, path) save_object(filtered_muts, path, legacy_name = if (legacy_rdata) "filtered_muts")
  thresholds = filtering_thresholds(filtering_type)
  #Sample IDs match from the start of the name, so e.g. PD1234 doesn't also match PD12345 colonies
  donor_pattern = paste0("^", donor_id)
  recipient_pattern = paste0("^", recipient_id)
  id = paste0(run_id, "_", filtering_type)
  for (folder in c("filtered_muts", "mpboot_files", "tree_files", "annotated_muts", "mutation_vcfs", "plots")) {
    dir.create(file.path(output_dir, folder), recursive = TRUE, showWarnings = FALSE)
  }
  path = function(folder, name) file.path(output_dir, folder, name)
  filtered_muts_file = path("filtered_muts", paste0("filtered_muts_", id, ".rds"))
  fasta = path("mpboot_files", paste0("Filter_", id, ".fa"))
  fasta_post_mix = path("mpboot_files", paste0("Filter_", id, "_post_mix.fa"))
  out = list(tree = path("tree_files", paste0("tree_", id, "_post_mix_post_dup.tree")),
             annotated_muts = path("annotated_muts", paste0("annotated_mut_set_", id, "_post_mix_post_dup.rds")),
             final_plots = path("plots", paste0(id, "_final_summary_trees.pdf")))

  #Remove low coverage / excluded samples and their private mutations
  if (min_sample_mean_cov > 0 || length(exclude_samples) > 0) {
    removed = remove_low_coverage_samples(COMB_mats = COMB_mats, filter_params = filter_params,
                                          min_sample_mean_cov = min_sample_mean_cov, other_samples_to_remove = exclude_samples,
                                          min_variant_reads_auto = 3, min_variant_reads_xy = 2)
    COMB_mats = removed$COMB_mats
    filter_params = removed$filter_params
  }

  filter_and_genotype = function(retain_muts, exclude_muts) {
    fm = get_filtered_mut_set(input_set_ID = run_id, COMB_mats = COMB_mats, filter_params = filter_params, gender = COMB_mats$gender,
                              retain_muts = retain_muts, exclude_muts = exclude_muts,
                              germline_pval = germline_log10_cutoff, rho = bb_cutoff,
                              mean_depth = unname(depth_cutoffs[c("auto_low", "auto_high", "xy_low", "xy_high")]),
                              pval_dp2 = NA, pval_dp3 = 0.001, min_depth = c(6, 4),
                              min_pval_for_true_somatic = thresholds$min_pval_for_true_somatic, min_vaf = thresholds$min_vaf,
                              min_variant_reads_SHARED = 2,
                              min_pval_for_true_somatic_SHARED = thresholds$min_pval_for_true_somatic_SHARED,
                              min_vaf_SHARED = thresholds$min_vaf_SHARED)
    prep = prepare_tree_input(fm$Genotype_shared_bin, filter_y = filter_y, CN_table = CN_table, snvs_only = snvs_only, build = build)
    fm$Genotype_shared_bin = prep$Genotype_shared_bin
    fm$dna_strings = prep$dna_strings
    save_muts(fm, filtered_muts_file)
    fm
  }
  build_and_assign = function(filtered_muts, fasta) {
    write.fasta(filtered_muts$dna_strings, names = names(filtered_muts$dna_strings), fasta)
    tree = run_mpboot(fasta, mpboot = mpboot, bootstraps = mpboot_bootstraps, seed = seed)
    assign_mutations_to_branches(tree, filtered_muts, keep_ancestral = keep_ancestral, create_multi_tree = polytomous, treefit_pval_cutoff = treefit_pval_cutoff)
  }

  # ---- 1-2. Filter, build the initial tree and assign mutations ----
  filtered_muts = filter_and_genotype(keep_muts, exclude_muts)
  res = build_and_assign(filtered_muts, fasta)
  filtered_muts = store_assignment(filtered_muts, res)
  tree = res$tree

  # ---- 3a. Remove poorly fitting indels ----
  #Indels that fit the tree poorly AND have borderline rho values, or that are assigned to the ancestral branch
  details = filtered_muts$COMB_mats.tree.build$mat
  low_pval_indels = intersect(details$mut_ref[details$pval < 1e-5 & details$Mut_type == "INDEL"], rownames(filter_params)[filter_params$bb_rhoval < 0.15])
  remove_muts = low_pval_indels
  if ("Ancestral" %in% tree$tip.label) {
    ROOT = tree$edge[1, 1]
    ancestral_branch = tree$edge[tree$edge[, 1] == ROOT & tree$edge[, 2] != which(tree$tip.label == "Ancestral"), 2]
    remove_muts = unique(c(details$mut_ref[details$node == ancestral_branch & details$Mut_type == "INDEL"], low_pval_indels))
  }
  #All mutations that fit the tree poorly: shown as a heatmap whether or not they are removed
  poor_fit = details$mut_ref[details$pval < treefit_pval_cutoff]
  plot_poor_fit_heatmap(tree, filtered_muts$COMB_mats.tree.build, poor_fit, file = path("plots", paste0(id, "_poor_fit_mutations.pdf")))
  write.table(details[details$mut_ref %in% poor_fit, c("mut_ref", "Mut_type", "node", "pval")], path("plots", paste0(id, "_poor_fit_mutations.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  if (remove_poor_fit && length(poor_fit) > 0) {
    message("Removing ", length(poor_fit), " mutations that fit the tree poorly (treemut p < ", treefit_pval_cutoff, ")")
    remove_muts = unique(c(remove_muts, poor_fit))
  }
  if (length(remove_muts) > 0) {
    exclude_muts = c(exclude_muts, remove_muts)
    message("Removing ", length(remove_muts), " mutations: poorly fitting indels with borderline rho values < 0.15, ancestral-branch indels",
            if (remove_poor_fit) ", and poorly fitting mutations")
    pdf(file = path("plots", paste0(id, "_muts_removed_as_artefact.pdf")), width = 15, height = 7)
    for (mut in remove_muts) {
      t = plot_tree(tree, cex.label = 0)
      add_annotation(tree = t, details = details, matrices = list(NV = filtered_muts$COMB_mats.tree.build$NV, NR = filtered_muts$COMB_mats.tree.build$NR),
                     annot_function = plot_MAV_mut, mut1 = mut)
    }
    dev.off()
  }

  # ---- 3b. Add back early mutations wrongly removed by the germline filter ----
  false_germline_muts = NULL
  if (germline_addback) {
    candidate_filter = NULL
    if (addback_require_filters) {
      candidate_filter = function(mut_refs, clade) {
        fp = filter_params[mut_refs, ]
        rows = match(mut_refs, COMB_mats$mat$mut_ref)
        hap = haploid_rows(COMB_mats)[rows]
        depth_ok = fp$mean_depth >= ifelse(hap, depth_cutoffs[["xy_low"]], depth_cutoffs[["auto_low"]]) &
                   fp$mean_depth <= ifelse(hap, depth_cutoffs[["xy_high"]], depth_cutoffs[["auto_high"]])
        others = setdiff(colnames(COMB_mats$NV), clade)
        nv = rowSums(COMB_mats$NV[rows, others, drop = FALSE]); nr = rowSums(COMB_mats$NR[rows, others, drop = FALSE])
        clonal_elsewhere = pbinom(nv, pmax(nr, 1), expected_vaf(COMB_mats)[rows]) > 0.001 #VAF outside the clade not too low for a clonal mutation
        depth_ok & fp$bb_rhoval > bb_cutoff & clonal_elsewhere
      }
    }
    false_germline_muts = unlist(check_for_false_germline_calls(tree, COMB_mats, filter_params, max_clade_prop = 0.1, CN_table = CN_table,
                                                                germline_log10_cutoff = germline_log10_cutoff, candidate_filter = candidate_filter))
  }
  if (!is.null(false_germline_muts) || length(remove_muts) > 0) {
    message("Saving the initial tree, then rebuilding with false germline calls added and artefacts removed")
    write.tree(tree, file = path("tree_files", paste0("tree_", id, "_v1.tree")))
    #mutations excluded by the user (or removed above) are not added back
    filtered_muts = filter_and_genotype(setdiff(c(keep_muts, false_germline_muts), exclude_muts), exclude_muts)
    res = build_and_assign(filtered_muts, fasta)
    filtered_muts = store_assignment(filtered_muts, res)
    tree = edge_lengths_from_assignment(res$tree, filtered_muts$COMB_mats.tree.build$mat)
  }
  write.tree(tree, file = path("tree_files", paste0("tree_", id, ".tree")))

  # ---- 4. Annotation ----
  annot_cols = c("Gene", "Transcript", "RNA", "CDS", "Protein", "Type", "SO_codes")
  if (is.null(annotate)) {
    filtered_muts$COMB_mats.tree.build$mat[annot_cols] = NA_character_
  } else {
    filtered_muts$COMB_mats.tree.build$mat = cbind(filtered_muts$COMB_mats.tree.build$mat, annotate(filtered_muts$COMB_mats.tree.build$mat)[annot_cols])
  }
  save_muts(filtered_muts, path("annotated_muts", paste0("annotated_mut_set_", id, ".rds")))

  # ---- 5. Remove non-clonal samples ----
  details = filtered_muts$COMB_mats.tree.build$mat
  matrices = list(mtr = filtered_muts$COMB_mats.tree.build$NV, dep = filtered_muts$COMB_mats.tree.build$NR)
  if (remove_nonclonal) {
    exclude_rows = NULL
    if (clonality_exclude_cn_xy) {
      hap = haploid_rows(filtered_muts$COMB_mats.tree.build)
      exclude_rows = function(sample) {
        out = hap
        cn = if (is.null(CN_table)) NULL else CN_table[CN_table$Sample == sample, , drop = FALSE]
        for (k in seq_len(if (is.null(cn)) 0 else nrow(cn))) {
          out = out | (strip_chr(details$Chrom) == strip_chr(cn$Chrom[k]) & details$Pos > cn$Pos_min[k] & details$Pos < cn$Pos_max[k])
        }
        out
      }
    }
    sample_test_df = test_sample_clonality(tree, details, matrices, keep_ancestral, exclude_rows = exclude_rows)
    write.table(sample_test_df, path("plots", paste0(id, "_clonality_test.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    mixed_samples = as.character(sample_test_df$sample[sample_test_df$result == "FAIL"])
    if (length(mixed_samples) == 0) {
      message("No samples failed the test for clonality")
    } else {
      message("Removing samples that failed the test for clonality: ", paste(mixed_samples, collapse = ", "))
      pdf(path("plots", paste0(id, "_mixed_vaf_trees.pdf")), width = 30, height = 10)
      for (sample in mixed_samples) {
        t = plot_tree(tree, cex.label = 0)
        plot_tree_vaf(t, details = details, matrices = matrices, samples = sample)
      }
      dev.off()

      filtered_muts$COMB_mats.tree.build = drop_samples_from_mut_set(filtered_muts$COMB_mats.tree.build, mixed_samples)
      Genotype_shared_bin = filtered_muts$COMB_mats.tree.build$Genotype_bin[rowSums(filtered_muts$COMB_mats.tree.build$Genotype_bin == 1) > 1, , drop = FALSE]
      prep = prepare_tree_input(Genotype_shared_bin, filter_y = filter_y, CN_table = CN_table, snvs_only = snvs_only, build = build)
      filtered_muts$Genotype_shared_bin = prep$Genotype_shared_bin
      filtered_muts$dna_strings = prep$dna_strings

      message("Rerunning tree-building now that mixed samples have been removed")
      res = build_and_assign(filtered_muts, fasta_post_mix)
      filtered_muts = store_assignment(filtered_muts, res)
      tree = edge_lengths_from_assignment(res$tree, filtered_muts$COMB_mats.tree.build$mat)
    }
    write.tree(tree, file = path("tree_files", paste0("tree_", id, "_post_mix.tree")))
    save_muts(filtered_muts, path("annotated_muts", paste0("annotated_mut_set_", id, "_post_mix.rds")))
  }

  # ---- 6. Correct the edge lengths for sample sensitivity ----
  #Correction is done BEFORE removing duplicates; duplicates are then dropped without changing branch lengths
  details = filtered_muts$COMB_mats.tree.build$mat
  sensitivity_df = read_sensitivity(sensitivity_path)
  tips = setdiff(tree$tip.label, "Ancestral")
  missing = setdiff(tips, sensitivity_df$Sample[!is.na(sensitivity_df$SNV_sensitivity)])
  if (length(missing) > 0) stop("No SNV sensitivity for: ", paste(missing, collapse = ", "), " - check ", sensitivity_path)
  if (anyNA(sensitivity_df$INDEL_sensitivity[sensitivity_df$Sample %in% tips])) {
    warning("Some samples have no indel sensitivity (too few germline indels?) - indel counts are left uncorrected")
    sensitivity_df$INDEL_sensitivity[is.na(sensitivity_df$INDEL_sensitivity)] = 1
  }
  tree_c = get_corrected_tree(tree = tree, details = details, include_indels = TRUE, sensitivity_df = sensitivity_df)
  tree_SNV_c = get_corrected_tree(tree = tree, details = details, include_indels = FALSE, sensitivity_df = sensitivity_df)
  tree_INDEL_c = get_corrected_tree(tree = tree, details = details, include_SNVs = FALSE, sensitivity_df = sensitivity_df)

  groups = function(t) highlight_groups(t, group1 = grep(donor_pattern, t$tip.label, value = TRUE), group2 = grep(recipient_pattern, t$tip.label, value = TRUE))
  pdf(path("plots", paste0(id, "_summary_trees_with_duplicates.pdf")), width = 30, height = 10)
  par(mfrow = c(2, 2))
  plot(tree, show.tip.label = FALSE, main = "Uncorrected tree: SNVs and INDELs", direction = "downwards", edge.color = groups(tree))
  plot(tree_c, show.tip.label = FALSE, main = "Corrected tree: SNVs and INDELs", direction = "downwards", edge.color = groups(tree))
  plot(tree_SNV_c, show.tip.label = FALSE, main = "Corrected tree: SNVs only", direction = "downwards", edge.color = groups(tree))
  plot(tree_INDEL_c, show.tip.label = FALSE, main = "Corrected tree: INDELs only", direction = "downwards", edge.color = groups(tree))
  par(mfrow = c(1, 1))
  dev.off()

  # ---- 7. Remove duplicate samples ----
  drop_samples = find_duplicate_samples(tree, details, min_private = duplicate_min_private, prop = duplicate_prop)
  if (length(drop_samples) == 0) {
    message("No duplicate samples identified")
  } else {
    message("Dropping as one of a duplicate set: ", paste(drop_samples, collapse = ", "))
    filtered_muts$COMB_mats.tree.build = drop_samples_from_mut_set(filtered_muts$COMB_mats.tree.build, drop_samples)
    tree_SNV_c = drop.tip(tree_SNV_c, drop_samples)

    #Re-assign mutations to the corrected tree without changing its branch lengths
    df = reconstruct_genotype_summary(tree_SNV_c)
    m = treemut_inputs(filtered_muts$COMB_mats.tree.build, df, keep_ancestral)
    p.error = sapply(df$samples, function(x) ifelse(x == "Ancestral", 1e-6, 0.01))
    res = assign_mutations(m$mtr, m$depth, df, error_rate = p.error)
    filtered_muts$COMB_mats.tree.build$mat$pval = res$summary$pval
    filtered_muts$COMB_mats.tree.build$mat$node = tree_SNV_c$edge[res$summary$edge_ml, 2]
  }
  write.tree(tree_SNV_c, file = out$tree)

  filtered_muts$COMB_mats.tree.build$mat = annotate_coding_changes(filtered_muts$COMB_mats.tree.build$mat, genes_to_annotate)
  save_muts(filtered_muts, out$annotated_muts)

  plot_final_trees(tree_SNV_c, filtered_muts, donor_id = donor_pattern, recipient_id = recipient_pattern, donor_age = donor_age, file = out$final_plots)
  out
}

#Final summary plots: corrected and ultrametric trees, with coding changes in genes of interest labelled.
#If there are recipient samples (transplant studies), donor and recipient trees are also drawn separately.
plot_final_trees = function(tree_SNV_c, filtered_muts, donor_id, recipient_id, donor_age, file) {
  details = filtered_muts$COMB_mats.tree.build$mat
  matrices = list(mtr = filtered_muts$COMB_mats.tree.build$NV, dep = filtered_muts$COMB_mats.tree.build$NR)
  if (is.na(donor_age)) donor_age = mean_sample_burden(tree_SNV_c) #If no age is provided, scale to the mean mutation burden

  ultra = function(t) {
    u = make.ultrametric.tree(t)
    u$edge.length = u$edge.length * donor_age
    u$edge.length[u$edge[, 2] == which(u$tip.label == "Ancestral")] = 0
    u
  }
  D_or_R = ifelse(grepl(donor_id, tree_SNV_c$tip.label), "D", ifelse(grepl("Ancestral", tree_SNV_c$tip.label), "Ancestral", "R"))
  tree_SNV_c$coords = NULL #reset the coords so that plot_tree will recalculate them
  trees = list(tree_SNV_c, ultra(tree_SNV_c))
  if (!all(D_or_R %in% c("D", "Ancestral"))) {
    tree_R = drop.tip(tree_SNV_c, tree_SNV_c$tip.label[D_or_R == "D"])
    tree_D = drop.tip(tree_SNV_c, tree_SNV_c$tip.label[D_or_R == "R"])
    trees = c(trees, list(tree_R, tree_D, ultra(tree_R), ultra(tree_D)))
  }

  sharing = function(t) add_annotation(t, details, matrices, annot_function = plot_sharing_info, donor_ID = donor_id, recip_ID = recipient_id,
                                       sharing_cols = c("black", "#11a0aa80", "#c8256580"))
  driver_labels = function(t) plot_tree_labels(t, details = details, type = "line", query.field = "coding_change_chip",
                                               data.frame(value = "Coding change mutation in driver", col = "red", pch = 17, stringsAsFactors = FALSE),
                                               label.field = "variant_ID", cex.label = 0.7, lty = 2)

  pdf(file = file, width = 15, height = 7)
  for (i in seq_along(trees)) {
    t = plot_tree(trees[[i]], cex.label = 0)
    sharing(t)
    if (i <= 2) driver_labels(t)
  }
  #All coding changes shared by >1 sample, on the ultrametric tree
  t = plot_tree(trees[[2]], cex.label = 0)
  details$coding_change_shared = ifelse(sapply(details$node, function(node) length(getTips(tree_SNV_c, node)) > 1) & details$coding_change == "Coding change", "yes", "no")
  sharing(t)
  plot_tree_labels(t, details = details, type = "label", query.field = "coding_change_shared",
                   data.frame(value = "yes", col = "red", lwd = 1, pch = 1, stringsAsFactors = FALSE),
                   label.field = "variant_ID", cex.label = 0.5)
  dev.off()
}
