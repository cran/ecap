## ----eval=FALSE---------------------------------------------------------------
#  library(devtools)
#  install_github("bradleyrava/ecap")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(ecap)

## -----------------------------------------------------------------------------
x <- ecap::elections_2018
p_obs <- x$Democrat_WinProbability[x$version=="classic"]
win_var <- x$Democrat_Won[x$version=="classic"]

## -----------------------------------------------------------------------------
set.seed(1)
train_rows <- sample(1:length(p_obs), length(p_obs)/2)
test_rows <- (1:length(p_obs))[-train_rows]

p_obs_train <- p_obs[train_rows]
win_var_train <- win_var[train_rows]

p_obs_test <- p_obs[test_rows]

## -----------------------------------------------------------------------------
ecap_fit <- ecap(unadjusted_prob = p_obs_train, win_var = win_var_train, win_id = 1, bias_indicator = T)

## -----------------------------------------------------------------------------
print(ecap_fit)

summary(ecap_fit)

## -----------------------------------------------------------------------------
plot(ecap_fit)

## -----------------------------------------------------------------------------
prob_ecap <- predict(object=ecap_fit, new_unadjusted=p_obs_test)
head(cbind.data.frame(ecap_p=prob_ecap, unadjusted_p=p_obs_test))

## -----------------------------------------------------------------------------
x <- ecap::elections_2018
y <- ecap::espn

