test_that("vectorised binomial tests equal the original binom.test loops", {
  cm = simulate_COMB_mats("male")
  cm$gender = "male"
  pv = pval_matrix(cm)
  NR = as.matrix(cm$NR); NR[NR == 0] = 1
  for (i in c(1, 500, nrow(NR))) for (j in c(1, 5)) {
    p = if (is_xy(cm$mat$Chrom[i])) 0.95 else 0.5
    expect_equal(pv[i, j], binom.test(cm$NV[i, j], NR[i, j], p = p, alternative = "less")$p.value)
  }
  g = germline.binomial.filter(cm)
  i = 900 #an X chromosome row
  expect_equal(g[i], binom.test(sum(cm$NV[i, ]), sum(cm$NR[i, ]), p = 0.95, alternative = "less")$p.value)
})

test_that("beta-binomial rho matches the original grid search", {
  cm = simulate_COMB_mats("female")
  rho = beta.binom.filter(cm)
  original = function(NV_vec, NR_vec) {
    rhovec = 10^seq(-6, -0.05, by = 0.05)
    mu = sum(NV_vec) / sum(NR_vec)
    ll = sapply(rhovec, function(rhoj) sum(VGAM::dbetabinom(x = NV_vec, size = NR_vec, rho = rhoj, prob = mu, log = TRUE)))
    rhovec[ll == max(ll)][1]
  }
  NR = as.matrix(cm$NR); NR[NR == 0] = 1
  for (i in c(1, 100, 400)) expect_equal(rho[i], original(as.numeric(cm$NV[i, ]), as.numeric(NR[i, ])))
})

test_that("upfront colony filtering keeps every sample when none fail (old -which() bug)", {
  cm = simulate_COMB_mats("female")
  res = compute_filter_params(cm, sex = "female", mixed_colony_vaf_cutoff = 0.1)
  expect_length(res$mixed_samples, 0)
  expect_equal(ncol(res$COMB_mats$NV), ncol(cm$NV))
  expect_equal(nrow(res$filter_params), nrow(res$COMB_mats$mat))
})

test_that("compute_filter_params refuses an unset sex", {
  expect_error(compute_filter_params(simulate_COMB_mats(), sex = "unknown"), "male' or 'female")
})

test_that("automatic depth window scales with coverage and widens with fewer colonies", {
  win = function(depth, n) {
    cm = simulate_COMB_mats("female", n_samples = n, depth = depth)
    cm$gender = "female"
    depth_cutoffs_auto(data.frame(mean_depth = get_mean_depth(cm)), cm)
  }
  w15 = win(15, 40); w30 = win(30, 40); w15_few = win(15, 8)
  expect_gt(w30[["auto_high"]], 1.7 * w15[["auto_high"]])
  expect_gt(w15_few[["auto_high"]] - w15_few[["auto_low"]], w15[["auto_high"]] - w15[["auto_low"]])
  expect_equal(w15[["xy_high"]], w15[["auto_high"]]) #female: X/Y treated as autosomes
})
