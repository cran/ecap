## Functions for implimenting the ECAP Adjustment

myns <- function (x, df = NULL, knots = NULL, intercept = FALSE,
                  Boundary.knots = range(x), deriv=0)
{
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax))
    x <- x[!nax]
  if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x >
                                                   Boundary.knots[2L])
  }
  else outside <- FALSE
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - 1L - intercept
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d",
                       1L + intercept), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(0, 1, length.out = nIknots + 2L)[-c(1L,
                                                           nIknots + 2L)]
      quantile(x[!outside], knots)
    }
  }
  else nIknots <- length(knots)
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
  if (any(outside)) {
    basis <- array(0, c(length(x), nIknots + 4L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, x[ol] - k.pivot)
      tt <- splines::splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0,
                                                        1), derivs=deriv)
      basis[ol, ] <- xl %*% tt
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, x[or] - k.pivot)
      tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0,
                                                        1), derivs=deriv)
      basis[or, ] <- xr %*% tt
    }
    if (any(inside <- !outside))
      basis[inside, ] <- splineDesign(Aknots, x[inside],
                                      4, derivs=deriv)
  }
  else basis <- splineDesign(Aknots, x, 4, derivs=deriv)
  const <- splineDesign(Aknots, Boundary.knots, 4, c(2, 2))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L),
                                                     drop = FALSE])
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3L, knots = if (is.null(knots)) numeric() else knots,
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("ns", "basis", "matrix")
  basis
}


risk_cvsplit_fcn <- function(i, partitions, basis_0, basis_1, probs_flip, lambda_grid, pt, omega, basis_0.grid, basis_1.grid) {
  test.rows <- partitions[[i]]
  train.rows <- assign_train(i, partitions)
  b_g <- basis_0[train.rows,]
  b_g_d1 <- basis_1[train.rows,]
  p_train <- probs_flip[train.rows]

  # I can do the below without thinking about lambda because it doesn't depend on it (only on test and train).
  # Calculate every term in the sum
  basis_sum_train <- t(b_g)%*%b_g

  # calculate the column wise sum of our the first derivative of our basis matrix
  sum_b_d1 <- t(b_g_d1)%*%rep(1,nrow(b_g_d1))

  # Now, we need to do this inversion for every value of lambda
  eta_g <- matrix(nrow = ncol(b_g), ncol = length(lambda_grid))
  for (ii in 1:length(lambda_grid)) {
    l <- lambda_grid[ii]
    ## Note the argument of gamma here is redundant without running the full QP
    eta_g[,ii] <- eta_min_fcn(l, 0.1, p_train, pt, omega, b_g, b_g_d1, basis_sum_train, basis_0.grid, basis_1.grid)
  }

  # Now, let's consider the test data
  b_g_test <- basis_0[test.rows,]
  b_g_d1_test <- basis_1[test.rows,]
  p_test <- probs_flip[test.rows]
  n <- nrow(b_g_test)

  # We just need to calculate the basis sum part again.
  basis_sum <- t(b_g_test)%*%b_g_test

  # Now, we can calculate the risk where eta was calculatd on our test data and we take the
  # basis matrix to calculate the score function from our test data
  one_vector <- as.vector(t(rep(1, nrow(b_g_test))))

  # Returns a vector of the risk for every value of lambda for a single test and train dataset
  risk_hat <- apply(eta_g, 2, function(eta) {
    g_hat <- b_g_test%*%eta
    g_hat_d1 <- b_g_d1_test%*%eta
    return((1/n)*sum(g_hat^2)+(2/n)*sum(g_hat*(1-2*p_test)+
                                          p_test*(1-p_test)*g_hat_d1))
  })
  # We now want to loop over this for all permutations of our cross validation
  return(risk_hat)
}


tweed.adj.fcn <- function(lambda, gamma, theta, p.tilde, p_flip, pt, probs, omega, basis_0, basis_1, basis_sum, basis_0.grid, basis_1.grid,
                          win_index, lose_index) {
  eta_hat <- eta_min_fcn(lambda, gamma, p_flip, pt, omega, basis_0, basis_1, basis_sum, basis_0.grid, basis_1.grid)
  g_hat <- basis_0%*%eta_hat
  g_hat_d1 <- basis_1%*%eta_hat

  mu_hat <-  p_flip + gamma*(g_hat + 1 - 2*p_flip)
  sigma2_hat <- gamma*p_flip*(1-p_flip) + gamma^2*p_flip*(1-p_flip)*(g_hat_d1-2)

  exp_p_hat <- mu_hat + 0.5*theta*(-mu_hat-6*mu_hat*sigma2_hat-2*mu_hat^3+3*sigma2_hat+3*mu_hat^2)
  var_p_hat <- (1-0.5*theta)^2*sigma2_hat + theta*sigma2_hat*(9*mu_hat^4*theta-18*mu_hat^3*theta+
                                                                9*mu_hat^2*theta-(1-theta/2)*(3*mu_hat^2-3*mu_hat))

  p.hat <- sapply((exp_p_hat + var_p_hat/exp_p_hat), function(uu) {min(uu,0.5)})
  ifelse(sum(p.hat<0)>0, p.hat <- p_flip, p.hat <- p.hat)

  # Flip back
  greater <- which(p.tilde > 0.5)
  p.hat[greater] <- 1-p.hat[greater]

  Q.gamma <- mle_binomial(p.hat, win_index, lose_index)

  return(Q.gamma)
}

eta_min_fcn <- function(lambda, gamma, p.tilde, pt, omega, basis_0, basis_1, basis_sum, basis_0.grid, basis_1.grid) {
  ## Define vars used below
  n <- nrow(basis_1)
  n_grid <- nrow(basis_1.grid)
  p_tilde <- p.tilde
  end_row <- which(pt == 0.5)

  ## Set up into the correct form
  Dmat <- 2*((1/n) * basis_sum + lambda * omega)
  dvec_terms <-  (1-2*p_tilde)*basis_0 + (p_tilde*(1-p_tilde))*basis_1
  dvec <- -(2/n)*colSums(dvec_terms)

  ## Constraint vectors
  b.vec1 <- 0
  bvec <- b.vec1

  Amat.part1 <- basis_0.grid[end_row,]
  Amat <- as.matrix(Amat.part1)

  return(quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq=1)$solution)
}

mle_binomial <- function(est, win_index, lose_index) {
  # Calculate wins with true p
  log.term <- sum(log(est[win_index]))
  minus.log.term <- sum(log(1-est[lose_index]))
  sum <- log.term+minus.log.term
  return(sum)
}

pFlip <- function(p) {
  if (p > 0.5) {
    p <- 1-p
  }
  return(p)
}

tweedie_est <- function(lambda, gamma, theta, p.tilde, p_flip, pt, omega, basis_0, basis_1, basis_sum,
                        basis_0.grid, basis_1.grid) {
  eta_hat <- eta_min_fcn(lambda, gamma, p_flip, pt, omega, basis_0, basis_1, basis_sum, basis_0.grid, basis_1.grid)
  g_hat <- basis_0%*%eta_hat
  g_hat_d1 <- basis_1%*%eta_hat

  mu_hat <-  p_flip + gamma*(g_hat + 1 - 2*p_flip)
  sigma2_hat <- gamma*p_flip*(1-p_flip) + gamma^2*p_flip*(1-p_flip)*(g_hat_d1-2)

  exp_p_hat <- mu_hat + 0.5*theta*(-mu_hat-6*mu_hat*sigma2_hat-2*mu_hat^3+3*sigma2_hat+3*mu_hat^2)
  var_p_hat <- (1-0.5*theta)^2*sigma2_hat + theta*sigma2_hat*(9*mu_hat^4*theta-18*mu_hat^3*theta+
                                                                9*mu_hat^2*theta-(1-theta/2)*(3*mu_hat^2-3*mu_hat))

  p.hat <- sapply((exp_p_hat + var_p_hat/exp_p_hat), function(uu) {min(uu, 0.5)})

  # Flip back
  greater <- which(p.tilde > 0.5)
  p.hat[greater] <- 1-p.hat[greater]

  return(p.hat)
}

assign_train <- function(test_index, partitions) {
  partitions_train <- partitions
  partitions_train[test_index] <- NA
  partitions_train <- unlist(partitions_train)
  partitions_train <- partitions_train[!is.na(partitions_train)]
  return(partitions_train)
}



