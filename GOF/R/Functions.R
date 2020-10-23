### Test statistics for Cramer von Mises- and Anderson Darling test for Goodness of fit

#' Cramer-von Mises test statistic
#'
#' Computes the Cramér-von Mises test statistic \deqn{\sum_{i=1}^{n}\left(Z_{i}-\frac{2i-1}{2n}\right)^2 + \frac{1}{12n},} for values \eqn{Z_{1},\ldots ,Z_{n}} where \eqn{Z_{i} = F(x_i)}, \eqn{F} is the assumed CDF and \eqn{x} is the ordered sample.
#'
#' @param x Vector of ordered real numbers
#' @return Value of test statistic
#' @examples
#' W2(c(0.4,0.5,0.6))
W2 <- function(x){
  n <- length(x)
  result <- sum((x - (seq(1, 2 * n - 1, 2) / (2 * n))) ** 2)
  return(result + 1 / (12 * n))
}


#' Anderson-Darling test statistic
#'
#' Computes the Anderson-Darling test statistic \deqn{A^2 = - n - \frac{1}{n}\sum_{i=1}^{n}(2i-1)(\log Z_{i}+\log (1-Z_{n+1-i})),} for values \eqn{Z_{1},\ldots ,Z_{n}} where \eqn{Z_{i} = F(x_i)}, \eqn{F} is the assumed CDF and \eqn{x} is the ordered sample.
#'
#' @param x Vector of ordered real numbers
#' @return Value of test statistic
#' @examples
#' A2(c(0.4, 0.5, 0.6))
A2 <- function(x){
  n <- length(x)
  return(- n - (1 / n) * sum(seq(1, 2 * n - 1, 2) * (log(x) + log(1 - rev(x)))))
}


### Create a 'print' method for class "gof_test"
print.gof_test <- function(x){
  cat('\nResults of Goodness-of-Fit Test\n-------------------------------\n\nTest Method:                     ',x$Method,'\n\nNull distribution:               ',x$Distribution,'\n\nSample Size:                     ',x$n,'\n\nEstimated Parameters:            ',x$Estimated,'\n\nParameters (and Estimations):    ',x$Parameters,'\n\nTest Statistic:                  ',x$Statistic_Name,' = ',x$Statistic,'\n\n\nIntervals for H0:\n',x$Critical,'', sep='')
}

### Distribution functions for the supported distributions

#' Cumulative distribution function of a Gumbel distribution
#'
#' Computes the CDF of a Gumbel distribution given as \deqn{F(x) = \exp\left(-\exp\left(-\frac{x-\alpha}{\beta}\right)\right).}
#'
#' @param x Vector of real numbers
#' @param location location parameter \eqn{\alpha} of the distribution
#' @param scale scale parameter \eqn{\beta > 0} of the distribution
#' @return Vector of values of the distribution function
#' @examples
#' F_extreme(0.5, 1, 1)
F_extreme <- function(x, location, scale){
  return(exp(-exp(-(x - location) / scale)))
}

#' Cumulative distribution function of a Generalized Pareto distribution
#'
#' Computes the CDF of a generalized Pareto distribution given as \deqn{F(x) = 1 - \left(1 - \frac{kx}{a}\right)^{\frac{1}{k}}.}
#'
#' @param x Vector of real numbers
#' @param scale scale parameter \eqn{a} of the distribution
#' @param shape shape parameter \eqn{k} of the distribution
#' @return Vector of values of the distribution function at \eqn{x}
#' @examples
#' F_extreme(c(0.5, 0.6), 1, 1)
# Generalized Pareto
F_gpd <- function(x, scale, shape){
  return(1 - (1 - shape * x / scale) ** (1 / shape))
}





### Estimation Functions

# Extreme value / Gumbel

#' Maximum Likelihood Estimator for a Gumbel distribution
#'
#' Computes Maximum Likelihood Estimators for any unknown parameters of a Gumbel distribution with CDF \deqn{F(x) = \exp\left(-\exp\left(-\frac{x-\alpha}{\beta}\right)\right),} for a sample \eqn{x}.
#'
#' @param x Vector of real numbers
#' @param param Vector of location \eqn{\alpha} and scale \eqn{\beta > 0} of the distribution. For the unknown parameter(s), put NA.
#' @return Vector of (estimated) parameters
#' @examples
#' est_extreme(c(0.5, 0.6, 1.1), c(NA, NA))
est_extreme <- function(x, param = c(NA,NA)){
  n <- length(x)

  # location unknown, scale known, Case 1
  if (is.na(param[1]) == TRUE & is.na(param[2]) != TRUE){
    location <- -param[2] * log(mean(exp(- x / param[2])))
    return(c(location, param[2]))
  }

  # location known, scale unknown, Case 2
  if (is.na(param[1]) != TRUE & is.na(param[2]) == TRUE){
    L <- function(scale){
      return(- sum((x - param[1]) / scale) - n * log(scale) - sum(exp(- (x - param[1]) / scale)))
    }
    scale_est <- optimize(L, c(0.00001, 100), maximum = TRUE)$maximum
    return(c(param[1], scale_est))

  }

  # both unknown, Case 3
  if (is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
    L <- function(params){
      location <- params[1]
      scale <- params[2]
      return(-(- sum((x - location) / scale) - n * log(scale) - sum(exp(- (x - location) / scale))))
    }
    location_est <- optim(c(1,1), L)$par[1]
    scale_est <- optim(c(1,1), L)$par[2]

    return(c(location_est, scale_est))
  }
}


# Pareto

#' Maximum Likelihood Estimator for a Pareto distribution
#'
#' Computes Maximum Likelihood Estimators for any unknown parameters of a Pareto distribution with CDF \deqn{F(x) = 1 - \left(\frac{x_m}{x}\right)^\alpha,} for (\eqn{x > x_m}), for a sample \eqn{x}.
#'
#' @param x Vector of real numbers
#' @param param Vector of scale \eqn{x_m > 0} and shape \eqn{\alpha > 0} of the distribution. For the unknown parameter(s), put NA.
#' @return Vector of (estimated) parameters
#' @examples
#' est_epareto(c(0.5, 0.6, 1.1), c(NA, NA))
est_pareto <- function(x, param = c(NA,NA)){
  n <- length(x)

  # scale unknown, shape known, Case 1
  if (is.na(param[1]) == TRUE & is.na(param[2]) != TRUE){
    scale <- min(x, na.rm = TRUE)
    return(c(scale, param[2]))
  }

  # scale known, shape unknown, Case 2
  if (is.na(param[1]) != TRUE & is.na(param[2]) == TRUE){
    shape <- n / sum(log(x / param[1]))
    return(c(param[1], shape))
  }

  # both unknown, Case 3
  if (is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
    scale <- min(x, na.rm = TRUE)
    shape <- n / sum(log(x / scale))
    return(c(scale, shape))
  }
}


# Generalized Pareto

#' Maximum Likelihood Estimator for a Generalized Pareto distribution
#'
#' Computes Maximum Likelihood Estimators for any unknown parameters of a Generalized Pareto distribution with CDF \deqn{F(x) = 1 - \left(1 - \frac{kx}{a}\right)^{\frac{1}{k}},} for a sample \eqn{x}.
#'
#' @param x Vector of real numbers
#' @param param Vector of scale \eqn{a > 0} and shape \eqn{k < 1} of the distribution. For the unknown parameter(s), put NA.
#' @return Vector of (estimated) parameters
#' @examples
#' est_gpd(c(0.5, 0.6, 1.1), c(NA, NA))
est_gpd <- function(x, param=c(NA,NA)){
  # Parameters given as scale > 0, shape < 1
  n <- length(x)

  # scale unknown, shape known, Case 1
  if (is.na(param[1]) == TRUE & is.na(param[2]) != TRUE){
    L <- function(a){
      if (param[2] == 0){
        return (- n * log(a) - sum(x / a))
      } else{
        return(- n * log(a) - (1 - 1 / param[2]) * sum(log(1 - (param[2] * x) / a)))
      }
    }
    scale <- optimize(L, interval = c(0.000001, 100), maximum = TRUE)$maximum
    return(c(scale, param[2]))
  }

  # scale known, shape unknown, Case 2
  if (is.na(param[1]) != TRUE & is.na(param[2]) == TRUE){
    L <- function(k){
      - n * log(param[1]) - (1 - 1 / k) * sum(log(1 - k * x / param[1]))
    }
    shape <- optimize(L, interval = c(- 100, 0.999999), maximum = TRUE)$maximum
    return(c(param[1], shape))
  }

  # both unknown, Case 3
  if (is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
    L <- function(theta){
      return(- n - sum(log(1 - theta * x)) - n * log(- (1 / (n * theta)) * sum(log(1 - theta * x))))
    }
    theta_max <- optimize(L, interval = c(-100, 1 / max(x) - 0.00001), maximum = TRUE)$maximum
    shape <- - (1 / n) * sum(log(1 - theta_max * x))
    scale <- shape / theta_max
    return(c(scale, shape))
  }
}

# Gamma
# Should not appear in documentation, implementation is inconsistent with others thus confusing

#' Maximum Likelihood Estimator for a Gamma distribution
#'
#' Computes Maximum Likelihood Estimators for any unknown parameters of a Gamma distribution with CDF \deqn{F(x) = \frac{1}{\beta\Gamma(m)}\left(\frac{x}{\beta}\right)^{m-1}}, for a sample \eqn{x}.
#'
#' @param x Vector of real numbers
#' @param param Vector of scale \eqn{\beta > 0} and shape \eqn{m > 0} of the distribution. For the unknown parameter(s), put NA.
#' @return Vector of (estimated) parameters
#' @examples
#' est_gamma(c(0.5, 0.6, 1.1), c(NA, NA))
#' #' @keywords internal
est_gamma <- function(x, param = c(NA,NA)){
  # Parameters given as scale > 0, shape > 0
  n <- length(x)

  # scale unknown, shape known, Case 1
  if (is.na(param[1]) == TRUE & is.na(param[2]) != TRUE){
    return(sum(x) / (param[2] * n))
  }

  # scale known, shape unknown, Case 2
  if (is.na(param[1]) != TRUE & is.na(param[2]) == TRUE){
    return(NA)
  }

  # both unknown, Case 3
  if (is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
    return(NA)
  }
}

### Normal Distribution
# Parameters are given as mu, sigmasq

#' Goodness-of-fit Test for the normal distribution
#'
#' Computes either Cramér-von Mises or Anderson-Darling test for Goodness-of-Fit of the normal distribution for a sample \eqn{x}. Unknown parameters are first estimated. The procedure accounts for shifts of critical values due to parameter estimation.
#'
#' @param x Vector of real numbers
#' @param statistic Either 'cvm' or 'ad'. Specifies which of the two tests to compute
#' @param param vector of mean \eqn{\mu} and variance \eqn{\sigma^2} of the distribution. For the unknown parameter(s), put NA.
#' @return Object of class gof_test
#' @examples
#' sim <- rnorm(100, 0, 1)
#' gof_normal(sim, 'cvm', c(NA, NA))
gof_normal <- function(x, statistic, param = c(NA, NA)){
  x <- sort(x)
  n <- length(x)

  # Cramer von Mises
  if(statistic == "cvm"){

    # mu unknown, sigma known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
      mu_est <- mean(x)
      w = (x - mu_est) / sqrt(param[2])
      Z <- pnorm(w, 0, 1)
      stat <- W2(Z)
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "normal", n = n, Estimated = "mu", Parameters = c(as.character(mu_est),", ",as.character(param[2])), Statistic = stat, Statistic_Name = "W^2", Critical = ".25: [0,.094] ; .15: [0,.117] ; .10: [0,.134] ; .05: [0,.165] ; .025: [0,.197] ; .01: [0,.238] ; .005: [0,.270] ; .0.0025: [0,.302]")
      class(out) <- "gof_test"
      return(out)
    }

    # mu known, sigma unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      s1sq <- sum((x - param[1])**2) / n
      w <- (x - param[1]) / sqrt(s1sq)
      Z <- pnorm(w, 0, 1)
      stat <- W2(Z)
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "normal", n = n, Estimated = "sigma^2", Parameters = c(as.character(param[1]),", ",as.character(s1sq)), Statistic = stat, Statistic_Name = "W^2", Critical = ".25: [0,.190] ; .15: [0,.263] ; .10: [0,.327] ; .05: [0,.442] ; .025: [0,.562] ; .01: [0,.725] ; .005: [0,.851] ; .0.0025: [0,.978]")
      class(out) <- "gof_test"
      return(out)
    }

    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      mu_est <- mean(x)
      ssq <- var(x)
      w <- (x - mu_est) / sqrt(ssq)
      Z <- pnorm(w, 0, 1)
      stat <- W2(Z) * (1 + 0.5 / n)
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "normal", n = n, Estimated = "mu, sigma^2", Parameters = c(as.character(mu_est),", ",as.character(ssq)), Statistic = stat, Statistic_Name = "W^2", Critical = ".25: [0,.074] ; .15: [0,.091] ; .10: [0,.104] ; .05: [0,.126] ; .025: [0,.148] ; .01: [0,.179] ; .005: [0,.201]")
      class(out) <- "gof_test"
      return(out)
    }
  }



  # Anderson Darling
  if(statistic == "ad"){

    # mu unknown, sigma known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
      mu_est <- mean(x)
      w = (x - mu_est) / sqrt(param[2])
      Z <- pnorm(w, 0, 1)
      stat <- A2(Z)
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "normal", n = n, Estimated = "mu", Parameters = c(as.character(mu_est),", ",as.character(param[2])), Statistic = stat, Statistic_Name = "A^2", Critical = ".25: [0,.644] ; .15: [0,.782] ; .10: [0,.894] ; .05: [0,1.087] ; .025: [0,1.285] ; .01: [0,1.551] ; .005: [0,1.756] ; .0.0025: [0,1.964]")
      class(out) <- "gof_test"
      return(out)
    }

    # mu known, sigma unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      s1sq <- sum((x - param[1])**2) / n
      w <- (x - param[1]) / sqrt(s1sq)
      Z <- pnorm(w, 0, 1)
      stat <- A2(Z)
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "normal", n = n, Estimated = "sigma^2", Parameters = c(as.character(param[1]),", ",as.character(s1sq)), Statistic = stat, Statistic_Name = "A^2", Critical = ".25: [0,1.072] ; .15: [0,1.430] ; .10: [0,1.743] ; .05: [0,2.308] ; .025: [0,2.898] ; .01: [0,3.702] ; .005: [0,4.324] ; .0.0025: [0,4.954]")
      class(out) <- "gof_test"
      return(out)
    }

    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      mu_est <- mean(x)
      ssq <- var(x)
      w <- (x - mu_est) / sqrt(ssq)
      Z <- pnorm(w, 0, 1)
      stat <- A2(Z) * (1 + 0.75 / n + 2.25 / (n**2))
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "normal", n = n, Estimated = "mu, sigma^2", Parameters = c(as.character(mu_est),", ",as.character(ssq)), Statistic = stat, Statistic_Name = "A^2", Critical = ".25: [0,.470] ; .15: [0,.561] ; .10: [0,.631] ; .05: [0,.752] ; .025: [0,.873] ; .01: [0,1.035] ; .005: [0,1.159]")
      class(out) <- "gof_test"
      return(out)
    }
  }
}





### Exponential distribution

#' Goodness-of-fit Test for the exponential distribution
#'
#' Computes either Cramér-von Mises or Anderson-Darling test for Goodness-of-Fit of the exponential distribution for a sample \eqn{x}. Unknown parameters are first estimated. The CDF is given by \deqn{F(x) = 1 - \exp(-\lambda x).}
#'
#' @param x Vector of real numbers
#' @param statistic Either 'cvm' or 'ad'. Specifies which of the two tests to compute
#' @return Object of class gof_test
#' @examples
#' gof_exp(c(0.02, 0.03, 1.1), 'cvm')
gof_exp <- function(x, statistic) {
  x <- sort(x)
  n <- length(x)

  # Cramer von Mises
  if(statistic == "cvm"){

    # Lambda not known, Case 1
    lambda_est <- 1 / mean(x)
    Z <- pexp(x, rate = lambda_est)
    stat <- W2(Z) * (1 + 0.16 / n)
    out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "exponential", n = n, Estimated = "lambda", Parameters = as.character(lambda_est), Statistic = stat, Statistic_Name = "W^2", Critical = ".25: [0,.116] ; .15: [0,.148] ; .10: [0,.175] ; .05: [0,.222] ; .025: [0,.271] ; .01: [0,0.338] ; .005: [0,0.390]")
    class(out) <- "gof_test"
    return(out)
  }

  # Anderson Darling
  if(statistic == "ad"){

    # Lambda not known, Case 1
    lambda_est <- 1 / mean(x)
    Z <- pexp(x, rate = lambda_est)
    stat <- A2(Z) * (1 + 0.6 / n)
    out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "exponential", n = n, Estimated = "lambda", Parameters = as.character(lambda_est), Statistic = stat, Statistic_Name = "A^2", Critical = ".25: [0,.736] ; .15: [0,.916] ; .10: [0,1.062] ; .05: [0,1.321] ; .025: [0,1.591] ; .01: [0,1.959] ; .005: [0,2.244]")
    class(out) <- "gof_test"
    return(out)
  }
}





### Uniform Distribution
# No parameters accepted, no statistic accepted

#' Goodness-of-fit Test for the uniform distribution
#'
#' Computes the Cramér-von Mises test for Goodness-of-Fit of the uniform distribution for a sample \eqn{x}.
#' @param x Vector of real numbers
#' @return Object of class gof_test
#' @examples
#' sim <- runif(100, 0, 1)
#' gof_uniform(sim)
gof_uniform <- function(x){
  x <- sort(x)
  n <- length(x)

  # Only Cramer von Mises is available
  # Only acceptable case is case of no known parameters
  a <- min(x)
  b <- max(x)
  Z <- punif(x, min = a, max = b)
  stat <- W2(Z)
  out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "uniform", n = n, Estimated = "xmin, xmax", Parameters = c(as.character(a),", ",as.character(b)), Statistic = stat, Statistic_Name = "W^2", Critical = "https://www.tandfonline.com/doi/pdf/10.1080/03610919608813360?casa_token=L9tnjIeXo-oAAAAA:tbZlrkJvhnm4iZ9TQ-Q1wgpreds5kHiWYGUDPMC7eRKYoJynt2gJAoYDFHbC1n9e5vCoPgKKcYl8")
  class(out) <- "gof_test"
  return(out)
}





### Extreme-value / Gumbel Distribution
# Parameters are given as location, scale (alpha, beta)

#' Goodness-of-fit Test for the Gumbel/Extreme value distribution
#'
#' Computes either Cramér-von Mises or Anderson-Darling test for Goodness-of-Fit of the Gumbel/Extreme value distribution for a sample \eqn{x}. Unknown parameters are first estimated. The CDF is given by \deqn{F(x) = \exp\left(-\exp\left(-\frac{x-\alpha}{\beta}\right)\right).}
#'
#' @param x Vector of real numbers
#' @param statistic Either 'cvm' or 'ad'. Specifies which of the two tests to compute
#' @param param vector of location \eqn{\alpha} and scale \eqn{\beta > 0} of the distribution. For the unknown parameter(s), put NA.
#' @return Object of class gof_test
#' @examples
#' gof_extremel(c(1, 2, 3), 'ad', c(NA, NA))
gof_extreme <- function(x, statistic, param = c(NA, NA)){
  x <- sort(x)
  n <- length(x)

  # Cramer von Mises
  if(statistic == "cvm"){

    # location unknown, scale known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
      location <- est_extreme(x, c(NA, param[2]))[1]
      Z <- F_extreme(x, location, param[2])
      stat <- W2(Z) * (1 + 0.16 / n)
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Extreme Value", n = n, Estimated = "location", Parameters = c(as.character(location),", ",as.character(param[2])), Statistic = stat, Statistic_Name = "W^2", Critical = ".25: [0,.116] ; .10: [0,.175] ; .05: [0,.222] ; .025: [0,.271] ; .01: [0,.338]")
      class(out) <- "gof_test"
      return(out)
    }

    # location known, scale unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      scale <- est_extreme(x, c(param[1],NA))[2]
      Z <- F_extreme(x, param[1], scale)
      stat <- W2(Z)
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Extreme Value", n = n, Estimated = "scale", Parameters = c(as.character(param[1]),", ",as.character(scale)), Statistic = stat, Statistic_Name = "W^2", Critical = ".25: [0,.186]; .10: [0,.320] ; .05: [0,.431] ; .025: [0,.547] ; .01: [0,.705]")
      class(out) <- "gof_test"
      return(out)
    }

    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      location <- est_extreme(x)[1]
      scale <- est_extreme(x)[2]
      Z <- F_extreme(x, location, scale)
      stat <- W2(Z) * (1 + 0.2 / sqrt(n))
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Extreme Value", n = n, Estimated = "location, scale", Parameters = c(as.character(location),", ",as.character(scale)), Statistic = stat, Statistic_Name = "W^2", Critical = ".25: [0,.073]; .10: [0,.102] ; .05: [0,.124] ; .025: [0,.146] ; .01: [0,.175]")
      class(out) <- "gof_test"
      return(out)
    }
  }



  # Anderson Darling
  if(statistic == "ad"){

    # location unknown, scale known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
      location <- est_extreme(x, c(NA, param[2]))[1]
      Z <- F_extreme(x, location, param[2])
      stat <- A2(Z) * (1 + 0.3 / n)
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Extreme Value", n = n, Estimated = "location", Parameters = c(as.character(location),", ",as.character(param[2])), Statistic = stat, Statistic_Name = "A^2", Critical = ".25: [0,.736]; .10: [0,1.062] ; .05: [0,1.321] ; .025: [0,1.591] ; .01: [0,1.959]")
      class(out) <- "gof_test"
      return(out)
    }

    # location known, scale unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      scale <- est_extreme(x, c(param[1],NA))[2]
      Z <- F_extreme(x, param[1], scale)
      stat <- A2(Z)
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Extreme Value", n = n, Estimated = "scale", Parameters = c(as.character(param[1]),", ",as.character(scale)), Statistic = stat, Statistic_Name = "A^2", Critical = ".25: [0,.1.06]; .10: [0,.1.725] ; .05: [0,2.277] ; .025: [0,2.854] ; .01: [0,3.640]")
      class(out) <- "gof_test"
      return(out)
    }

    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      location <- est_extreme(x, c(NA, NA))[1]
      scale <- est_extreme(x, c(NA, NA))[2]
      Z <- F_extreme(x, location, scale)
      stat <- A2(Z) * (1 + 0.2 / sqrt(n))
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Extreme Value", n = n, Estimated = "location, scale", Parameters = c(as.character(location),", ",as.character(scale)), Statistic = stat, Statistic_Name = "A^2", Critical = ".25: [0,.474]; .10: [0,.637] ; .05: [0,.757] ; .025: [0,.877] ; .01: [0,1.038]")
      class(out) <- "gof_test"
      return(out)
    }
  }
}





### Weibull
# Parameters are given as scale, shape (beta, k)
# Just perform a transformation to get an extreme value distribution, then call
# the function gof_extreme

#' Goodness-of-fit Test for the Weibull distribution
#'
#' Computes either Cramér-von Mises or Anderson-Darling test for Goodness-of-Fit of the Weibull distribution for a sample \eqn{x}. Unknown parameters are first estimated. The CDF is given by \deqn{F(x) = 1 - \exp\left(\left(\frac{x}{\beta}\right)^{m}\right).}
#'
#' @param x Vector of real numbers
#' @param statistic Either 'cvm' or 'ad'. Specifies which of the two tests to compute
#' @param param vector of scale \eqn{\beta > 0} and shape \eqn{m > 1} of the distribution. For the unknown parameter(s), put NA.
#' @return Object of class gof_test
#' @examples
#' gof_weibulll(c(1,2,3), 'cvm', c(NA, NA))
gof_weibull <- function(x, statistic, param = c(NA, NA)){
  x <- sort(x)
  n <- length(x)
  y <- sort(-log(x))


  # Cramer von Mises
  if(statistic == "cvm"){

    # scale unknown, shape known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
      extreme_value_trans <- gof_extreme(y, "cvm", c(NA, 1 / param[2]))
      scale <- exp(as.numeric(extreme_value_trans$Parameters[1]))
      result <- extreme_value_trans$Statistic
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "scale", Parameters = c(as.character(scale),", ",as.character(param[2])), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.116] ; .10: [0,.175] ; .05: [0,.222] ; .025: [0,.271] ; .01: [0,.338]")
      class(out) <- "gof_test"
      return(out)
    }

    # scale known, shape unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      extreme_value_trans <- gof_extreme(y, "cvm", c(- log(1 / param[1]), NA))
      shape <- 1 / as.numeric(extreme_value_trans$Parameters[3])
      result <- extreme_value_trans$Statistic
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "shape", Parameters = c(as.character(param[1]),", ",as.character(shape)), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.186]; .10: [0,.320] ; .05: [0,.431] ; .025: [0,.547] ; .01: [0,.705]")
      class(out) <- "gof_test"
      return(out)
    }

    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      extreme_value_trans <- gof_extreme(y, "cvm", c(NA,NA))
      scale <- exp(as.numeric(extreme_value_trans$Parameters[1]))
      shape <- 1 / as.numeric(extreme_value_trans$Parameters[3])
      result <- extreme_value_trans$Statistic
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "scale, shape", Parameters = c(as.character(scale),", ",as.character(shape)), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.073]; .10: [0,.102] ; .05: [0,.124] ; .025: [0,.146] ; .01: [0,.175]")
      class(out) <- "gof_test"
      return(out)
    }
  }



  # Anderson Darling
  if(statistic == "ad"){

    # scale unknown, shape known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
      extreme_value_trans <- gof_extreme(y, "ad", c(NA, 1 / param[2]))
      scale <- exp(as.numeric(extreme_value_trans$Parameters[1]))
      result <- extreme_value_trans$Statistic
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "scale", Parameters = c(as.character(scale),", ",as.character(param[2])), Statistic = result, Statistic_Name = "A^2", Critical = ".25: [0,.736]; .10: [0,1.062] ; .05: [0,1.321] ; .025: [0,1.591] ; .01: [0,1.959]")
      class(out) <- "gof_test"
      return(out)
    }

    # scale known, shape unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      extreme_value_trans <- gof_extreme(y, "ad", c(- log(1 / param[1]), NA))
      shape <- 1 / as.numeric(extreme_value_trans$Parameters[3])
      result <- extreme_value_trans$Statistic
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "shape", Parameters = c(as.character(param[1]),", ",as.character(shape)), Statistic = result, Statistic_Name = "A^2", Critical = ".25: [0,.1.06]; .10: [0,.1.725] ; .05: [0,2.277] ; .025: [0,2.854] ; .01: [0,3.640]")
      class(out) <- "gof_test"
      return(out)
    }

    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      extreme_value_trans <- gof_extreme(y, "ad", c(NA,NA))
      scale <- exp(as.numeric(extreme_value_trans$Parameters[1]))
      shape <- 1 / as.numeric(extreme_value_trans$Parameters[3])
      result <- extreme_value_trans$Statistic
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "scale, shape", Parameters = c(as.character(scale),", ",as.character(shape)), Statistic = result, Statistic_Name = "A^2", Critical = ".25: [0,.474]; .10: [0,.637] ; .05: [0,.757] ; .025: [0,.877] ; .01: [0,1.038]")
      class(out) <- "gof_test"
      return(out)
    }
  }
}





### Pareto
# Parameters are given as scale, shape (xmin, alpha)
# Just perform a transformation to get an exponential distribution, then call
# the function exponential_gof
# Can only implement Case 2 atm. Case 1 and 3 should not give accurate critical values
# as the transformation depends on estimated parameter

#' Goodness-of-fit Test for the Pareto distribution
#'
#' Computes either Cramér-von Mises or Anderson-Darling test for Goodness-of-Fit of the Pareto distribution for a sample \eqn{x}. The CDF is given as \deqn{F(x) = 1 - \left(\frac{x_m}{x}\right)^\alpha,} (\eqn{x > x_m}). The only case that is implemented is the one of known scale \eqn{x_m} and unknown shape \eqn{\alpha}.
#'
#' @param x Vector of real numbers
#' @param statistic Either 'cvm' or 'ad'. Specifies which of the two tests to compute
#' @param param scale parameter \eqn{x_m} of the distribution
#' @return Object of class gof_test
#' @examples
#' sim <- unif(100, 1, 1)
#' gof_pareto(sim, 'cvm', 1)
gof_pareto <- function(x, param){
  x <- sort(x)
  n <- length(x)

  # Cramer von Mises
  if(statistic == "cvm"){

    # scale known, shape unknown, Case 2
    y <- log(x / param)
    exp_trans <- gof_exp(y, "cvm")
    shape <- exp_trans$Parameters
    result <- exp_trans$Statistic
    out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Pareto", n = n, Estimated = "shape", Parameters = c(as.character(param),", ",as.character(shape)), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.116] ; .15: [0,.148] ; .10: [0,.175] ; .05: [0,.222] ; .025: [0,.271] ; .01: [0,0.338] ; .005: [0,0.390]")
    class(out) <- "gof_test"
    return(out)
  }


  # Anderson Darling
  if(statistic == "ad"){

    # scale known, shape unknown, Case 2
    y <- log(x / param)
    exp_trans <- gof_exp(y, "ad")
    shape <- exp_trans$Parameters
    result <- exp_trans$Statistic
    out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Pareto", n = n, Estimated = "shape", Parameters = c(as.character(param),", ",as.character(shape)), Statistic = result, Statistic_Name = "A^2", Critical = ".25: [0,.736] ; .15: [0,.916] ; .10: [0,1.062] ; .05: [0,1.321] ; .025: [0,1.591] ; .01: [0,1.959] ; .005: [0,2.244]")
    class(out) <- "gof_test"
    return(out)
  }
}





### Generalized Pareto distribution
# Parameters are given as scale, shape (a,k)
# Case 2 only works for negative shapes which is undesirable.

#' Goodness-of-fit Test for the Generalized Pareto distribution
#'
#' Computes either Cramér-von Mises or Anderson-Darling test for Goodness-of-Fit of the GPD for a sample \eqn{x}. Unknown parameters are first estimated. The CDF is given as \deqn{F(x) = 1 - \left(1 - \frac{kx}{a}\right)^{\frac{1}{k}}.}
#'
#' @param x Vector of real numbers
#' @param statistic Either 'cvm' or 'ad'. Specifies which of the two tests to compute
#' @param param vector of scale \eqn{a > 0} and shape \eqn{k < 1} of the distribution. For the unknown parameter(s), put NA.
#' @return Object of class gof_test
#' @examples
#' gof_gpd(c(1,2,3), 'ad', c(NA, NA))
gof_gpd <- function(x, statistic, param = c(NA, NA)){
  x <- sort(x)
  n <- length(x)

  # Cramer von Mises
  if(statistic == "cvm"){

    # scale unknown, shape known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
      scale_est <- est_gpd(x, param = c(NA, param[2]))[1]
      Z <- F_gpd(x, scale_est, param[2])
      stat <- W2(Z)
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Generalized Pareto", n = n, Estimated = "scale", Parameters = c(as.character(scale_est),", ",as.character(param[2])), Statistic = stat, Statistic_Name = "W^2", Critical = "https://www.tandfonline.com/doi/pdf/10.1198/00401700152672573?casa_token=FKsCzUZmFBYAAAAA:jAf6TceZfULvHuXR3_fu80ykBxgiNhzxTumiQVVFL_iUbUt0Uw_SVawlOFml_cQcp0EXg6IJgkJH")
      class(out) <- "gof_test"
      return(out)
    }

    # scale known, shape unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      shape_est <- est_gpd(x, param = c(param[1], NA))[2]
      Z <- F_gpd(x, param[1], shape_est)
      stat <- W2(Z)
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Generalized Pareto", n = n, Estimated = "shape", Parameters = c(as.character(param[1]),", ",as.character(shape_est)), Statistic = stat, Statistic_Name = "W^2", Critical = "https://www.tandfonline.com/doi/pdf/10.1198/00401700152672573?casa_token=FKsCzUZmFBYAAAAA:jAf6TceZfULvHuXR3_fu80ykBxgiNhzxTumiQVVFL_iUbUt0Uw_SVawlOFml_cQcp0EXg6IJgkJH")
      class(out) <- "gof_test"
      return(out)
    }

    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      scale_est <- est_gpd(x, param = c(NA, NA))[1]
      shape_est <- est_gpd(x, param = c(NA, NA))[2]
      Z <- F_gpd(x, scale_est, shape_est)
      stat <- W2(Z)
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Generalized Pareto", n = n, Estimated = "scale, shape", Parameters = c(as.character(scale_est),", ",as.character(shape_est)), Statistic = stat, Statistic_Name = "W^2", Critical = "https://www.tandfonline.com/doi/pdf/10.1198/00401700152672573?casa_token=FKsCzUZmFBYAAAAA:jAf6TceZfULvHuXR3_fu80ykBxgiNhzxTumiQVVFL_iUbUt0Uw_SVawlOFml_cQcp0EXg6IJgkJH")
      class(out) <- "gof_test"
      return(out)
    }
  }



  # Anderson Darling
  if(statistic == "ad"){


    # scale unknown, shape known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
      scale_est <- est_gpd(x, param = c(NA, param[2]))[1]
      Z <- F_gpd(x, scale_est, param[2])
      stat <- A2(Z)
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Generalized Pareto", n = n, Estimated = "scale", Parameters = c(as.character(scale_est),", ",as.character(param[2])), Statistic = stat, Statistic_Name = "A^2", Critical = "https://www.tandfonline.com/doi/pdf/10.1198/00401700152672573?casa_token=FKsCzUZmFBYAAAAA:jAf6TceZfULvHuXR3_fu80ykBxgiNhzxTumiQVVFL_iUbUt0Uw_SVawlOFml_cQcp0EXg6IJgkJH")
      class(out) <- "gof_test"
      return(out)
    }

    # scale known, shape unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      shape_est <- est_gpd(x, param = c(param[1], NA))[2]
      Z <- F_gpd(x, param[1], shape_est)
      stat <- A2(Z)
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Generalized Pareto", n = n, Estimated = "shape", Parameters = c(as.character(param[1]),", ",as.character(shape_est)), Statistic = stat, Statistic_Name = "A^2", Critical = "https://www.tandfonline.com/doi/pdf/10.1198/00401700152672573?casa_token=FKsCzUZmFBYAAAAA:jAf6TceZfULvHuXR3_fu80ykBxgiNhzxTumiQVVFL_iUbUt0Uw_SVawlOFml_cQcp0EXg6IJgkJH")
      class(out) <- "gof_test"
      return(out)
    }

    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      scale_est <- est_gpd(x, param = c(NA, NA))[1]
      shape_est <- est_gpd(x, param = c(NA, NA))[2]
      Z <- F_gpd(x, scale_est, shape_est)
      stat <- A2(Z)
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Generalized Pareto", n = n, Estimated = "scale, shape", Parameters = c(as.character(scale_est),", ",as.character(shape_est)), Statistic = stat, Statistic_Name = "A^2", Critical = "https://www.tandfonline.com/doi/pdf/10.1198/00401700152672573?casa_token=FKsCzUZmFBYAAAAA:jAf6TceZfULvHuXR3_fu80ykBxgiNhzxTumiQVVFL_iUbUt0Uw_SVawlOFml_cQcp0EXg6IJgkJH")
      class(out) <- "gof_test"
      return(out)
    }
  }
}





### Gamma
# Parameters are given as scale, shape (beta, m)

#' Goodness-of-fit Test for the Gamma distribution
#'
#' Computes either Cramér-von Mises or Anderson-Darling test for Goodness-of-Fit of the Gamma distribution for a sample \eqn{x}. Unknown parameters are first estimated. The CDF is given as \deqn{F(x) = \frac{1}{\beta\Gamma(m)}\left(\frac{x}{\beta}\right)^{m-1}.}
#'
#' @param x Vector of real numbers
#' @param statistic Either 'cvm' or 'ad'. Specifies which of the two tests to compute
#' @param param vector of scale \eqn{\beta > 0} and shape \eqn{m > 0} of the distribution. For the unknown parameter(s), put NA.
#' @return Object of class gof_test
#' @examples
#' gof_gamma(c(1,2,3), 'ad', c(NA, NA))
gof_gamma <- function(x, statistic, param = c(NA, NA)){
  x <- sort(x)
  n <- length(x)

  # Cramer von Mises
  if(statistic == "cvm"){

    # scale unknown, shape known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
      scale_est <- est_gamma(x, param = c(NA, param[2]))
      Z <- pgamma(x, scale = scale_est, shape = param[2])
      W <- W2(Z)
      result <- W
      if (param[2] == 1){
        result <- W * (1 + 0.16 / n)
      }
      if (param[2] >= 2){
        result <- (1.8 * n * W - 0.14) / (1.8 * n - 1)
      }
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Gamma", n = n, Estimated = "scale", Parameters = c(as.character(scale_est),", ",as.character(param[2])), Statistic = result, Statistic_Name = "W^2", Critical = "check documentation for critical values")
      class(out) <- "gof_test"
      return(out)
    }

    # scale known, shape unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      L <- function(m){
        (m - 1) * sum(log(x)) - sum(x / param[1]) - n * m * log(param[1]) - n * log(gamma(m))
      }
      shape_est <- optimize(L, c(0.00001, 100), maximum = TRUE)$maximum
      Z <- pgamma(x, scale = param[1], shape = shape_est)
      result <- W2(Z)
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Gamma", n = n, Estimated = "shape", Parameters = c(as.character(param[1]),", ",as.character(shape_est)), Statistic = result, Statistic_Name = "W^2", Critical = "check documentation for critical values")
      class(out) <- "gof_test"
      return(out)
    }

    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      L <- function(m){
        (m - 1) * sum(log(x)) - n * m - n * m * log(sum(x / (m * n))) - n * log(gamma(m))
      }
      shape_est <- optimize(L, c(0.00001, 100), maximum = TRUE)$maximum
      scale_est <- mean(x) / shape_est
      Z <- pgamma(x, scale = scale_est, shape = shape_est)
      result <- W2(Z)
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Gamma", n = n, Estimated = "scale, shape", Parameters = c(as.character(scale_est),", ",as.character(shape_est)), Statistic = result, Statistic_Name = "W^2", Critical = "check documentation for critical values")
      class(out) <- "gof_test"
      return(out)
    }
  }



  # Anderson Darling
  if(statistic == "ad"){

    # scale unknown, shape known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
      scale_est <- est_gamma(x, param = c(NA, param[2]))
      Z <- pgamma(x, scale = scale_est, shape = param[2])
      A <- A2(Z)
      result <- A
      if (param[2] == 1){
        result <- A * (1 + 0.6 / n)
      }
      if (param[2] >= 2){
        result <- A + (1 / n) * (0.2 + (0.3 / param[2]))
      }
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Gamma", n = n, Estimated = "scale", Parameters = c(as.character(scale_est),", ",as.character(param[2])), Statistic = result, Statistic_Name = "A^2", Critical = "check documentation for critical values")
      class(out) <- "gof_test"
      return(out)
    }

    # scale known, shape unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      L <- function(m){
        (m - 1) * sum(log(x)) - sum(x / param[1]) - n * m * log(param[1]) - n * log(gamma(m))
      }
      shape_est <- optimize(L, c(0.00001, 100), maximum = TRUE)$maximum
      Z <- pgamma(x, scale = param[1], shape = shape_est)
      result <- A2(Z)
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Gamma", n = n, Estimated = "shape", Parameters = c(as.character(param[1]),", ",as.character(shape_est)), Statistic = result, Statistic_Name = "A^2", Critical = "check documentation for critical values")
      class(out) <- "gof_test"
      return(out)
    }

    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      L <- function(m){
        (m - 1) * sum(log(x)) - n * m - n * m * log(sum(x / (m * n))) - n * log(gamma(m))
      }
      shape_est <- optimize(L, c(0.00001, 100), maximum = TRUE)$maximum
      scale_est <- mean(x) / shape_est
      Z <- pgamma(x, scale = scale_est, shape = shape_est)
      result <- A2(Z)
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Gamma", n = n, Estimated = "scale, shape", Parameters = c(as.character(scale_est),", ",as.character(shape_est)), Statistic = result, Statistic_Name = "A^2", Critical = "check documentation for critical values")
      class(out) <- "gof_test"
      return(out)
    }
  }
}

#' A_GOF package
#'
#' @description When dealing with a problem regarding Goodness-of-Fit, one must first check if the distribution under the null hypothesis
#' is fully specified. If all parameters are known prior and need not be estimated, the standard EDF tests may be used. This package does not cover
#' this case, so the reader is referred to the package 'goftest' (available on CRAN). If (some) parameters must be estimated, check if the distribution
#' in question is listed in the table of contents. If so, the corresponding function can be used to estimate the unspecified parameters and subsequently
#' calculate the desired test statistic. The null hypothesis of a sample belonging to the hypothesized distribution is rejected at some level if the test statistic does not lie in the interval given in the test
#' output for that particular level. Some estimation functions and functions for calculating CVM and AD statistics are also documented and may be used, but are
#' rather supplementary to the main test functions. Mathematical derivations of the methods used in this package are mostly given in D'Agostino, Stephens (1986).
#' @author Alexander Lappe
#' @references D'Agostino, Stephens (1986) Goodness-of-Fit Techniques
#' @references Choulakian, Stephens (2012) Goodness-of-Fit Tests for the Generalized Pareto Distribution
#' @docType package
#' @name A GOF Package
NULL


