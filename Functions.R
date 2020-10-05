library(tidyverse)
library(EnvStats)
library(rootSolve)


# Test statistics for Cramer von Mises- and Anderson Darling test for Goodness of fit
W2 <- function(x){
  n <- length(x)
  result <- sum((x - (seq(1, 2 * n - 1, 2) / (2 * n))) ** 2)
  return(result + 1 / (12 * n))
}

A2 <- function(x){
  n <- length(x)
  return(- n - (1 / n) * sum(seq(1, 2 * n - 1, 2) * (log(x) + log(1 - rev(x)))))
}


# Create a 'print' method for class "gof_test"
print.gof_test <- function(x){
  cat('Results of Goodness-of-Fit Test\n-------------------------------\n\nTest Method:                     ',x$Method,'\n\nNull distribution:               ',x$Distribution,'\n\nSample Size:                     ',x$n,'\n\nEstimated Parameters:            ',x$Estimated,'\n\nParameters (and Estimations):    ',x$Parameters,'\n\nTest Statistic:                  ',x$Statistic_Name,' = ',x$Statistic,'\n\n\nIntervals for H0:\n',x$Critical,'', sep='')
}



### Distribution functions for the supported distributions
# Extreme value
F_extreme <- function(x, location, scale){
  return(exp(-exp(-(x - location) / scale)))
}


### Estimation Functions

# Extreme value
est_extreme <- function(x, param = c(NA,NA)){
  n <- length(x)
  
  # location unknown, scale known, Case 1
  if (is.na(param[1]) == TRUE & is.na(param[2]) != TRUE){
    location <- -param[2] * log(mean(exp(- x / param[2])))
    return(c(location, param[2]))
  }
  
  # location known, scale unknown, Case 2
  if (is.na(param[1]) != TRUE & is.na(param[2]) == TRUE){
    y <- x - param[1]
    estimator <- function(scale){
      return((sum(y) - sum(y * exp(- y / scale))) / n - scale)
    }
    scale <- uniroot(estimator, c(0,100))$root
    return(c(param[1], scale))
  }
  
  # both unknown, Case 3
  if (is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
    estimator <- function(scale){
      return(mean(x) - (sum(x * exp(- x / scale)) / sum(exp(- x / scale))) - scale)
    }
    scale <- uniroot(estimator, c(0.1,100))$root
    location <- - scale * log(1 / n * sum(exp(- x / scale)))
    return(c(location, scale))
  }
}


# Pareto
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

# Placeholder for generalized pareto from SEAL script
#
#
#
#
#
#
#

### Normal Distribution
# Parameters are given as mu, sigmasq
normal_gof <- function(x, statistic, param = c(NA, NA)){
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
exp_gof <- function(x, statistic) {
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
    out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "exponential", n = n, Estimated = "lambda", Parameters = as.character(lambda_est), Statistic = stat, Statistic_Name = "A^2", Critical = ".10: [.208,1.321] ; .05: [.178,1.591]")
    class(out) <- "gof_test"
    return(out)
  }
}


x)[2]
      Z <- F_extreme(x, location, scale)
      stat <- W2(Z) * (1 + 0.2 / sqrt(n))
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Extreme Value", n = n, Estimated = "location, scale", Parameters = c(as.character(location),", ",as.character(scale)), Statistic = stat, Statistic_Name = "W^2", Critical = ".25: [0,.736] ; .15: [0,.916] ; .10: [0,1.062] ; .05: [0,1.321] ; .025: [0,1.591] ; .01: [0,1.959] ; .005: [0,2.244]")
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


### Extreme-value Distribution
# Parameters are given as location, scale (alpha, beta)
extreme_gof <- function(x, statistic, param = c(NA, NA)){
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
# the function extreme_gof
weibull_gof <- function(x, statistic, param = c(NA, NA)){
  x <- sort(x)
  n <- length(x)
  y <- sort(-log(x))
  
  
  # Cramer von Mises
  if(statistic == "cvm"){
    
    # scale unknown, shape known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
    extreme_value_trans <- extreme_gof(y, "cvm", c(NA, 1 / param[2]))
    scale <- exp(as.numeric(extreme_value_trans$Parameters[1]))
    result <- extreme_value_trans$Statistic
    out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "scale", Parameters = c(as.character(scale),", ",as.character(param[2])), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.116] ; .10: [0,.175] ; .05: [0,.222] ; .025: [0,.271] ; .01: [0,.338]")
      class(out) <- "gof_test"
      return(out)
    }
    
    # scale known, shape unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
    extreme_value_trans <- extreme_gof(y, "cvm", c(- log(1 / param[1]), NA))
    shape <- 1 / as.numeric(extreme_value_trans$Parameters[3])
    result <- extreme_value_trans$Statistic
    out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "shape", Parameters = c(as.character(param[1]),", ",as.character(shape)), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.186]; .10: [0,.320] ; .05: [0,.431] ; .025: [0,.547] ; .01: [0,.705]")
      class(out) <- "gof_test"
      return(out)
    }
    
    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      extreme_value_trans <- extreme_gof(y, "cvm", c(NA,NA))
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
    extreme_value_trans <- extreme_gof(y, "ad", c(NA, 1 / param[2]))
    scale <- exp(as.numeric(extreme_value_trans$Parameters[1]))
    result <- extreme_value_trans$Statistic
    out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "scale", Parameters = c(as.character(scale),", ",as.character(param[2])), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.736]; .10: [0,1.062] ; .05: [0,1.321] ; .025: [0,1.591] ; .01: [0,1.959]")
      class(out) <- "gof_test"
      return(out)
    }
    
    # scale known, shape unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
    extreme_value_trans <- extreme_gof(y, "ad", c(- log(1 / param[1]), NA))
    shape <- 1 / as.numeric(extreme_value_trans$Parameters[3])
    result <- extreme_value_trans$Statistic
    out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "shape", Parameters = c(as.character(param[1]),", ",as.character(shape)), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.1.06]; .10: [0,.1.725] ; .05: [0,2.277] ; .025: [0,2.854] ; .01: [0,3.640]")
      class(out) <- "gof_test"
      return(out)
    }
    
    # both unknown, Case 3 
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      extreme_value_trans <- extreme_gof(y, "ad", c(NA,NA))
      scale <- exp(as.numeric(extreme_value_trans$Parameters[1]))
      shape <- 1 / as.numeric(extreme_value_trans$Parameters[3])
      result <- extreme_value_trans$Statistic
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "scale, shape", Parameters = c(as.character(scale),", ",as.character(shape)), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.474]; .10: [0,.637] ; .05: [0,.757] ; .025: [0,.877] ; .01: [0,1.038]")
      class(out) <- "gof_test"
      return(out)
    } 
  }
}


### Pareto
# Parameters are given as scale, shape (xmin, alpha)
# Just perform a transformation to get an exponential value distribution, then call     
# the function exponential_gof
# Can only implement Case 2 atm. Case 1 and 3 should not give accurate critical values
# as the transformation depends on estimated parameter
pareto_gof <- function(x, statistic, param = c(NA, NA)){
  x <- sort(x)
  n <- length(x)
  
  # Cramer von Mises
  if(statistic == "cvm"){
    
    # scale unknown, shape known, Case 1
    if(is.na(param[1]) == TRUE & is.na(param[2]) == FALSE){
    extreme_value_trans <- extreme_gof(y, "cvm", c(NA, 1 / param[2]))
    scale <- exp(as.numeric(extreme_value_trans$Parameters[1]))
    result <- extreme_value_trans$Statistic
    out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "scale", Parameters = c(as.character(scale),", ",as.character(param[2])), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.116] ; .10: [0,.175] ; .05: [0,.222] ; .025: [0,.271] ; .01: [0,.338]")
      class(out) <- "gof_test"
      return(out)
    }
    
    # scale known, shape unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      y <- log(x / param[1])
      exp_trans <- exp_gof(y, "cvm")
      shape <- exp_trans$Parameters
      result <- exp_trans$Statistic
      out <- list(Method = "Cramér-von Mises Test for Goodness-of-Fit", Distribution = "Pareto", n = n, Estimated = "shape", Parameters = c(as.character(param[1]),", ",as.character(shape)), Statistic = result, Statistic_Name = "W^2", Critical = ".10: [.0276,.222] ; .05: [.0233,.251]")
      class(out) <- "gof_test"
      return(out)
    }
    
    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      extreme_value_trans <- extreme_gof(y, "cvm", c(NA,NA))
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
    extreme_value_trans <- extreme_gof(y, "ad", c(NA, 1 / param[2]))
    scale <- exp(as.numeric(extreme_value_trans$Parameters[1]))
    result <- extreme_value_trans$Statistic
    out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "scale", Parameters = c(as.character(scale),", ",as.character(param[2])), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.736]; .10: [0,1.062] ; .05: [0,1.321] ; .025: [0,1.591] ; .01: [0,1.959]")
      class(out) <- "gof_test"
      return(out)
    }
    
    # scale known, shape unknown, Case 2
    if(is.na(param[1]) == FALSE & is.na(param[2]) == TRUE){
      y <- log(x / param[1])
      exp_trans <- exp_gof(y, "ad")
      shape <- exp_trans$Parameters
      result <- exp_trans$Statistic
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Pareto", n = n, Estimated = "shape", Parameters = c(as.character(param[1]),", ",as.character(shape)), Statistic = result, Statistic_Name = "W^2", Critical = ".10: [.0276,.222] ; .05: [.0233,.251]")
      class(out) <- "gof_test"
      return(out)
    }
    
    # both unknown, Case 3
    if(is.na(param[1]) == TRUE & is.na(param[2]) == TRUE){
      extreme_value_trans <- extreme_gof(y, "ad", c(NA,NA))
      scale <- exp(as.numeric(extreme_value_trans$Parameters[1]))
      shape <- 1 / as.numeric(extreme_value_trans$Parameters[3])
      result <- extreme_value_trans$Statistic
      out <- list(Method = "Anderson-Darling Test for Goodness-of-Fit", Distribution = "Weibull", n = n, Estimated = "scale, shape", Parameters = c(as.character(scale),", ",as.character(shape)), Statistic = result, Statistic_Name = "W^2", Critical = ".25: [0,.474]; .10: [0,.637] ; .05: [0,.757] ; .025: [0,.877] ; .01: [0,1.038]")
      class(out) <- "gof_test"
      return(out)
    } 
  }
}
