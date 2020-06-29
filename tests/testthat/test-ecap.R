test_that("ecap model generation and prediction works", {
  p_obs <- runif(1000, 0, 1)
  p_new <- runif(1000, 0, 1)
  win_var <- rbinom(length(p_obs), 1, p_obs)
  ecap_fit <- try(ecap(unadjusted_prob = p_obs, win_var = win_var, win_id = 1, bias_indicator = FALSE))
  ecap_fit_biased <- try(ecap(unadjusted_prob = p_obs, win_var = win_var, win_id = 1, bias_indicator = FALSE))
  ecap_new <- try(predict(object=ecap_fit, new_unadjusted=p_new))

  expect_equivalent(length(p_obs), 1000)
  expect_equivalent(class(ecap_fit), "ecap")
  expect_equivalent(class(ecap_fit_biased),"ecap")
  expect_equivalent(class(ecap_new),"numeric")
})
