test_that("sex is inferred correctly for simulated female and male donors", {
  f = check_donor_sex(simulate_COMB_mats("female"), "female", build = "hg38")
  expect_true(f$ok); expect_equal(f$sex, "female")
  m = check_donor_sex(simulate_COMB_mats("male"), "male", build = "hg38")
  expect_true(m$ok); expect_equal(m$sex, "male")
})

test_that("a female configured as male stops the run, or warns if not strict", {
  cm = simulate_COMB_mats("female")
  expect_error(check_donor_sex(cm, "male", build = "hg38"), "Sex check failed")
  expect_warning(r <- check_donor_sex(cm, "male", build = "hg38", strict = FALSE))
  expect_equal(r$sex, "male")
})

test_that("sex 'unknown' uses the inferred sex", {
  expect_equal(check_donor_sex(simulate_COMB_mats("male"), "unknown", build = "hg38")$sex, "male")
  expect_equal(check_donor_sex(simulate_COMB_mats("female"), "unknown", build = "hg38")$sex, "female")
})
