#' Estimate parameters for the ECAP adjustment
#'
#' This function estimates 3 parameters that are needed to preform the ECAP adjustment. As such, it is meant to be used
#' together with the predict.ecap() function included in the "ecap" package. The parameters estimated are the level of corruption
#' and bias in the given unadjusted probability estimates, along with a tuning parameter needed to adjust the level of
#' smoothness in the ECAP estimation.
#'
#' @param unadjusted_prob Numeric vector of probability estimates that you want estimate the ECAP parameters from.
#' @param win_var A binary vector of wins and losses that correspond to the probabilities in the unadjusted_prob vector
#' @param win_id A value that denotes a "win" (or if the event occurred) in the win_var vector.
#' @param bias_indicator Set this equal to F if you don't want to consider bias in your estimation. Set it equal to T if you do.
#' @param lambda_grid This is already predefined. However, you can adjust the grid of tuning parameters lambda that ECAP searches over if needed.
#' @param gamma_grid This is already predefined. However, you can adjust the grid of gamma that ECAP searches over if needed.
#' @param theta_grid This is already predefined. However, you can adjust the grid of theta that ECAP searches over if needed.
#' @return An ecap object that can be used to adjust new probability estimates. It contains all of the tuning parameters needed to calibrate
#' ECAP as well as diagnostic information on the estimate of g. The probabilities used to calibrate ECAP have also been ECAP corrected and
#' are given as part of the output.
#' @importFrom splines splineDesign
#' @importFrom quadprog solve.QP
#' @importFrom stats quantile
#' @author Bradley Rava, Peter Radchenko and Gareth M. James.
#' @references http://faculty.marshall.usc.edu/gareth-james/Research/Probs.pdf
#' @examples
#' \donttest{
#' set.seed(1)
#' p_obs <- runif(1000, 0, 1)
#' win_var <- rbinom(length(p_obs), 1, p_obs)
#' ecap_fit <- ecap(unadjusted_prob = p_obs, win_var = win_var, win_id = 1, bias_indicator = FALSE)
#' }
#' @export
ecap <- function(unadjusted_prob, win_var, win_id, bias_indicator = F, lambda_grid=10^seq(-6, 0, by=0.5),
                 gamma_grid=seq(0.001, 0.05, by=0.001), theta_grid=seq(-4, 2, 0.1)) {
  ## Warning messages
  if (min(theta_grid) < -4 | max(theta_grid) > 2) {
    return(warning("The bias (theta) must be greater than -4 or less than 2. Please adjust theta grid accordingly."))
  }
  if (min(gamma_grid) < 0) {
    return(warning("The corruption level (gamma) must be greater than 0. Please adjust gamma grid accordingly."))
  }
  if (min(lambda_grid) < 0) {
    return(warning("The smoothness parameter (lambda) must be greater than 0. Please adjust lambda grid accordingly."))
  }
  if (length(unadjusted_prob) != length(win_var)) {
    return(warning("The vector of wins and losses (win_var) must be the same length as the vector of unadjusted probability estimates."))
  }
  if (min(unadjusted_prob) < 0 | max(unadjusted_prob) > 1) {
    return(warning("Error. At least one of the unadjusted probabilities given is not a valid probability."))
  }

  ## Store the data
  probs <- cbind.data.frame(p.tilde=unadjusted_prob, greater_half = ifelse(unadjusted_prob > 0.5, T, F), win_var = win_var)

  ## Make sure we return the probabilities in the original order
  probs$orig_order <- 1:nrow(probs)

  ## Convert all probabilities to between 0 and 1/2
  p_flip <- sapply(probs$p.tilde, function(p) {
    if (p > 0.5) {
      p <- 1-p
    }
    return(p)
  })
  probs$p_flip <- p_flip
  probs <- probs[order(probs$p_flip),]

  ## For MLE Later -- obtain the win variable
  win_index <- which(probs$win_var == win_id)
  lose_index <- (1:length(win_var))[-win_index]

  ## Generate basis function / omega matrix from p.tilde
  # Get Knot Locations
  probs_flip <- probs$p_flip
  knot.range <- range(probs$p_flip)
  quantiles <- seq(knot.range[1]+0.0001, knot.range[2]-0.0001, length = 50)

  # Generate the basis matrix and its correspoding 1st and 2nd deriv's
  basis_0 <- myns(probs_flip, knots = quantiles, intercept = T, Boundary.knots = knot.range)
  basis_1 <- myns(probs_flip, knots = quantiles, deriv = 1, intercept = T, Boundary.knots = knot.range)
  basis_2 <- myns(probs_flip, knots = quantiles, deriv = 2, intercept = T, Boundary.knots = knot.range)
  basis_sum <- t(basis_0)%*%basis_0
  sum_b_d1 <- t(basis_1)%*%rep(1,nrow(basis_1))

  # We also want to calculate Omega on a fine grid of points
  fine_grid <- seq(0, 0.5, by=0.001)
  basis_fine_grid <- myns(fine_grid, knots = quantiles, intercept = T, Boundary.knots = knot.range)
  basis_fine_grid_d2 <- myns(fine_grid, knots = quantiles, deriv = 2, intercept = T, Boundary.knots = knot.range)
  omega <- (1/nrow(basis_fine_grid)) * (t(basis_fine_grid_d2) %*% basis_fine_grid_d2)

  ## Grid for the optimization algorithm
  pt <- c(seq(10^(-12), 0.5, by = 0.001), 0.5)
  basis_0.grid <- myns(pt, knots = quantiles, intercept = T, Boundary.knots = knot.range)
  basis_1.grid <- myns(pt, knots = quantiles, deriv = 1, intercept = T, Boundary.knots = knot.range)
  basis_sum.grid <- t(basis_0.grid)%*%basis_0.grid

  #### Risk function for lambda and grid for gamma####
  ## CV SET UP to get min value of lambda from risk function
  # Randomize the rows
  rows <- 1:nrow(basis_0)
  rows_rand <- rows[sample(rows)]

  ## Declate the number of groups that we want
  n.group <- 10
  ## Return a list with 10 approx equal vectors of rows.
  partitions <- split(rows_rand,
                      cut(rows_rand,quantile(rows_rand,(0:n.group)/n.group),
                          include.lowest=TRUE, labels=FALSE))

  for(jjj in 1:n.group) {
    partitions[[jjj]] <- rows_rand[partitions[[jjj]]]
  }

  ### Here we are going to pick the best value of lambda through cross validation
  risk_cvsplit <- lapply(1:n.group, function(j) {
    return(risk_cvsplit_fcn(j, partitions, basis_0, basis_1, probs$p_flip, lambda_grid, pt, omega, basis_0.grid, basis_1.grid))
  })

  r_cv_split_matrix <- do.call(cbind, risk_cvsplit)
  r_cv_split_vec <- apply(r_cv_split_matrix, 1, mean)
  names(r_cv_split_vec) <- lambda_grid

  ## Get the value of lambda that corresponds to the smallest risk
  lambda_opt <- lambda_grid[which.min(r_cv_split_vec)]


  ## 2D grid search for gamma and theta (1D if user specifies there is no bias)
  if (bias_indicator == F) {
    theta_grid <- 0
  }
  ## Iterate over every value of lambda and theta
  gamma_theta_matrix <- matrix(NA, ncol=length(theta_grid), nrow=length(gamma_grid))
  colnames(gamma_theta_matrix) <- theta_grid
  rownames(gamma_theta_matrix) <- gamma_grid

  for (i in 1:length(gamma_grid)) {
    g <- gamma_grid[i]
    for (j in 1:length(theta_grid)) {
      t <- theta_grid[j]
      score <- tweed.adj.fcn(lambda_opt, g, t, probs$p.tilde, probs$p_flip, pt,
                             probs, omega, basis_0, basis_1, basis_sum, basis_0.grid,
                             basis_1.grid, win_index, lose_index)
      gamma_theta_matrix[i,j] <- score
    }
  }

  ## Get the value of gamma and theta that maximized the MLE
  max_pair <- which(gamma_theta_matrix == max(gamma_theta_matrix), arr.ind = TRUE)
  gamma_opt <- gamma_grid[max_pair[1]]
  theta_opt <- theta_grid[max_pair[2]]

  ## ECAP adjust the training unadjusted probabilities
  ecap_adj_p <- tweedie_est(lambda_opt, gamma_opt, theta_opt, unadjusted_prob, p_flip, pt,
                            omega, basis_0, basis_1, basis_sum, basis_0.grid, basis_1.grid)

  eta_hat <- eta_min_fcn(lambda_opt, gamma_opt, p_flip, pt, omega, basis_0, basis_1, basis_sum, basis_0.grid, basis_1.grid)
  g_hat <- basis_0%*%eta_hat
  g_hat_d1 <- basis_1%*%eta_hat

  ## Object we return to the user
  return_data <- list(lambda=lambda_opt,
                      gamma=gamma_opt,
                      theta=theta_opt,
                      g_hat=as.vector(g_hat),
                      g_hat_d1=as.vector(g_hat_d1),
                      unadjusted_prob=unadjusted_prob,
                      unadjusted_flip=p_flip,
                      ecap_training_probabilities=ecap_adj_p,
                      lambda_grid=lambda_grid,
                      gamma_grid=gamma_grid,
                      theta_grid=theta_grid,
                      gamma_theta_matrix=gamma_theta_matrix,
                      lambda_cv_error=r_cv_split_vec)
  class(return_data) <- "ecap"

  return(return_data)
}


#' Printing ECAP Object
#'
#' Prints summary information about the ECAP object
#'
#' @param x An object of class ecap.
#' @param digits The number of significant digits that should be displayed.
#' @param ... Additional arguments
#' @author Bradley Rava, Peter Radchenko and Gareth M. James.
#' @references http://faculty.marshall.usc.edu/gareth-james/Research/Probs.pdf
#' @examples
#' \donttest{
#' set.seed(1)
#' p_obs <- runif(1000, 0, 1)
#' win_var <- rbinom(length(p_obs), 1, p_obs)
#' ecap_fit <- ecap(unadjusted_prob = p_obs, win_var = win_var, win_id = 1, bias_indicator = TRUE)
#' print(ecap_fit)
#' }
#' @export
print.ecap <- function(x, digits, ...) {
  if (missing(digits)) {
    digits <- 4
  }
  cat("\nECAP Adjustment Formula:\n",
      "E(p|p_tilde) + Var(p|p_tilde) / E(p|p_tilde) \n \n",
      "The optimal parameters picked for the ECAP adjustment are:",
      paste("lambda = ", round(x$lambda,digits), ", ", sep = ""),
      paste("gamma =", round(x$gamma,digits)), "and",
      paste("theta = ", round(x$theta,digits)))

  cat("\n")
}


#' Summary of ECAP Object
#'
#' Prints summary information about the calibration of an ECAP object.
#' *** Denotes that one of the parameter estimates has hit the end of the given grid of tuning parameters.
#' The grid can be adjusted in the ecap function.
#'
#' @param object An object of class ecap.
#' @param digits The number of significant digits that should be displayed.
#' @param ... Additional arguments
#' @author Bradley Rava, Peter Radchenko and Gareth M. James.
#' @references http://faculty.marshall.usc.edu/gareth-james/Research/Probs.pdf
#' @examples
#' \donttest{
#' set.seed(1)
#' p_obs <- runif(1000, 0, 1)
#' win_var <- rbinom(length(p_obs), 1, p_obs)
#' ecap_fit <- ecap(unadjusted_prob = p_obs, win_var = win_var, win_id = 1, bias_indicator = FALSE)
#' summary(ecap_fit)
#' }
#' @export
summary.ecap <- function(object, digits, ...) {
  if (missing(digits)) {
    digits <- 4
  }
  ## Head Text:
  title_ecap <- "\n ECAP Adjustment Summary: \n\n"
  ## Show given estimated values
  optimal_params <- paste("The optimal parameters picked for the ECAP adjustment are:",
                          paste("lambda = ", round(object$lambda,digits) , ", ", sep = ""),
                          paste("gamma =", round(object$gamma, digits)),
                          "and",
                          paste("theta = ", round(object$theta,digits)))

  ## Message for start of grid summary
  grid_summary_text <- "\n\n The optimal parameters were picked from these grids (*** estimate hit end of grid)"

  ## Summary of the grids used to pick the optimal parameters
  lambda_summary <- paste("\n Lambda Grid - Length", length(object$lambda_grid), "From", round(range(object$lambda_grid)[1],digits), "to",
                          round(range(object$lambda_grid)[2],digits))
  gamma_summary <- paste("\n Gamma Grid - Length", length(object$gamma_grid),  "From", round(range(object$gamma_grid)[1],digits), "to",
                         round(range(object$gamma_grid)[2],digits))
  theta_summary <- paste("\n Theta Grid - Length", length(object$theta_grid),  "From", round(range(object$theta_grid)[1],digits), "to",
                         round(range(object$theta_grid)[2],digits))

  ## Warnings in case a parameter hits the end of the grid search
  if (min(object$lambda_grid) == object$lambda | max(object$lambda_grid) == object$lambda) {
    lambda_summary <- NULL
    lambda_summary <- paste("\n Lambda Grid - Length", length(object$lambda_grid), "From", round(range(object$lambda_grid)[1],digits), "to",
                            round(range(object$lambda_grid)[2],digits), "***")
  }
  if (min(object$gamma_grid) == object$gamma | max(object$gamma_grid) == object$gamma) {
    gamma_summary <- NULL
    gamma_summary <- paste("\n Gamma Grid - Length", length(object$gamma_grid),  "From", round(range(object$gamma_grid)[1],digits), "to",
                           round(range(object$gamma_grid)[2],digits), "***")
  }
  if (min(object$theta_grid) == object$theta | max(object$theta_grid) == object$theta) {
    theta_summary <- NULL
    theta_summary <- paste("\n Theta Grid - Length", length(object$theta_grid),  "From", round(range(object$theta_grid)[1],digits), "to",
                           round(range(object$theta_grid)[2], digits), "***")
  }

  ## Summary of given probability estimates
  unadjusted_summary <- paste("\n \n The unadjusted probabilities range from", round(min(object$unadjusted_prob),digits), "to",
                              round(max(object$unadjusted_prob),digits), "with an average of", round(mean(object$unadjusted_prob),digits),
                              " - Length", length(object$unadjusted_prob))
  ## Summary of the ecap training probability estimates
  ecap_training_summary <- paste("\n The ECAP training probabilities range from", round(min(object$ecap_training_probabilities),digits), "to",
                                 round(max(object$ecap_training_probabilities),digits), "with an average of", round(mean(object$ecap_training_probabilities),digits),
                                 " - Length", length(object$ecap_training_probabilities))


  return(cat(title_ecap,
             optimal_params,
             grid_summary_text,
             lambda_summary,
             gamma_summary,
             theta_summary,
             unadjusted_summary,
             ecap_training_summary))
}

#' Plotting ECAP Object
#'
#' Plots diagnostic information of an ECAP object. Two plots are produced. The first plot displays the estimate of the
#' function g that the ecap procedure produced. The second compares the unadjusted probability estimates to the ECAP adjusted
#' probability estimates that were used to train the model.
#'
#' @param x An object of class ecap.
#' @param ... Additional arguments
#' @author Bradley Rava, Peter Radchenko and Gareth M. James.
#' @references http://faculty.marshall.usc.edu/gareth-james/Research/Probs.pdf
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme_minimal
#' @examples
#' \donttest{
#' set.seed(1)
#' p_obs <- runif(1000, 0, 1)
#' win_var <- rbinom(length(p_obs), 1, p_obs)
#' ecap_fit <- ecap(unadjusted_prob = p_obs, win_var = win_var, win_id = 1, bias_indicator = FALSE)
#' plot(ecap_fit)
#' }
#' @export
plot.ecap <- function(x, ...) {
  ## Empty containers
  p_tilde <- NULL
  g_hat <- NULL
  p <- NULL
  adjustment_type <- NULL

  ## Plot of estimate of g()
  plot_data <- cbind.data.frame(p_tilde=sort(x$unadjusted_flip), g_hat=x$g_hat)
  g_hat_plot <- ggplot2::ggplot(plot_data, aes(x=p_tilde, y=g_hat)) +
    ggplot2::geom_line() + ggplot2::xlab("Unadjusted probabilities flipped between 0 and 0.5") +
    ggplot2::ylab("Estimate of g") + ggplot2::ggtitle("Estimate of g on the unadjusted probabilities") +
    ggplot2::theme_minimal()

  ## Plot of the unadjusted and training ecap probabilities
  plot_data2 <- cbind.data.frame(p=c(x$unadjusted_prob, x$ecap_training_probabilities),
                                 adjustment_type=factor(rep(c("Unadjusted", "ECAP Training"), each=length(x$unadjusted_prob))))
  probs_plot <- ggplot2::ggplot(plot_data2, aes(x=p, color=adjustment_type)) +
    ggplot2::geom_line(stat="density") + ggplot2::theme_minimal() + ggplot2::scale_color_manual(values=c("#009E73","#D55E00")) +
    ggplot2::ggtitle("Density of unadjusted probabilities against the training ECAP adjusted") + xlab("Probability")

  ## Return object
  plot_return <- list(g_hat_plot, probs_plot)

  return(plot_return)
}

#' Implementing the ECAP procedure
#'
#' Takes in an ECAP object and a new set of probability estimates that the user wishes to adjust. The model uses the
#' calibration from the ecap object to ECAP adjust the new probability estimates given to the function predict.
#' @return A vector of ECAP adjusted probability estimates.
#' @param object An object of class ecap.
#' @param ... Additional arguments
#' @param new_unadjusted A numerical vector of unadjusted probabilities that you want to ECAP adjust.
#' @author Bradley Rava, Peter Radchenko and Gareth M. James.
#' @references http://faculty.marshall.usc.edu/gareth-james/Research/Probs.pdf
#' @importFrom utils tail
#' @examples
#' \donttest{
#' set.seed(1)
#' p_obs <- runif(1000, 0, 1)
#' win_var <- rbinom(length(p_obs), 1, p_obs)
#' ecap_fit <- ecap(unadjusted_prob = p_obs, win_var = win_var, win_id = 1, bias_indicator = FALSE)
#'
#' p_new <- runif(1000, 0, 1)
#' ecap_new <- predict(object=ecap_fit, new_unadjusted=p_new)
#' }
#' @export
predict.ecap <- function(object, new_unadjusted, ...) {
  if (min(new_unadjusted) < 0 | max(new_unadjusted) > 1) {
    return(warning("Error. At least one of the unadjusted probabilities given is not a valid probability."))
  }

  new_flip <- sapply(new_unadjusted, function(p) {
    if (p > 0.5) {
      p <- 1-p
    }
    return(p)
  })

  ## Combine new probs with the old ones
  p_old_new <- c(object$unadjusted_prob, new_unadjusted)
  p_old_new_flip <- c(object$unadjusted_flip, new_flip)
  probs_flip <- sort(p_old_new_flip)

  ## Generate basis function / omega matrix from p.tilde
  # Get Knot Locations
  knot.range <- range(probs_flip)
  quantiles <- seq(knot.range[1]+0.0001, knot.range[2]-0.0001, length = 50)

  # Generate the basis matrix and its correspoding 1st and 2nd deriv's
  basis_0 <- myns(probs_flip, knots = quantiles, intercept = T, Boundary.knots = knot.range)
  basis_1 <- myns(probs_flip, knots = quantiles, deriv = 1, intercept = T, Boundary.knots = knot.range)
  basis_2 <- myns(probs_flip, knots = quantiles, deriv = 2, intercept = T, Boundary.knots = knot.range)
  basis_sum <- t(basis_0)%*%basis_0
  sum_b_d1 <- t(basis_1)%*%rep(1,nrow(basis_1))

  # We also want to calculate Omega on a fine grid of points
  fine_grid <- seq(0, 0.5, by=0.001)
  basis_fine_grid <- myns(fine_grid, knots = quantiles, intercept = T, Boundary.knots = knot.range)
  basis_fine_grid_d2 <- myns(fine_grid, knots = quantiles, deriv = 2, intercept = T, Boundary.knots = knot.range)
  omega <- (1/nrow(basis_fine_grid)) * (t(basis_fine_grid_d2) %*% basis_fine_grid_d2)

  ## Grid for the optimization algorithm
  pt <- c(seq(10^(-12), 0.5, by = 0.001), 0.5)
  basis_0.grid <- myns(pt, knots = quantiles, intercept = T, Boundary.knots = knot.range)
  basis_1.grid <- myns(pt, knots = quantiles, deriv = 1, intercept = T, Boundary.knots = knot.range)
  basis_sum.grid <- t(basis_0.grid)%*%basis_0.grid

  ecap_old_new <- tweedie_est(object$lambda, object$gamma, object$theta, p_old_new, p_old_new_flip, pt,
                          omega, basis_0, basis_1, basis_sum, basis_0.grid, basis_1.grid)

  ecap_new <- tail(ecap_old_new, length(new_unadjusted))

  ## Return valid probabilities
  ecap_new[ecap_new < 0] <- 0
  ecap_new[ecap_new > 1] <- 1

  return(ecap_new)
}







