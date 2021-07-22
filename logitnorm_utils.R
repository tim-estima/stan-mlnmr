#########################################################################################
# Utility functions for the logit-Normal distribution.
#
#   pars_logitnorm : derive mu, sigma parameters for logit-Normal from sample mean, sd
#   qlogitnorm     : quantile function, taking sample mean and sd
#   dlogitnorm     : density funciton, taking sample mean and sd
#   plogitnorm     : distribution function, taking sample mean and sd
#
#########################################################################################


.lndiff <- function(est, smean, ssd){
  x <- logitnorm::momentsLogitnorm(est[1], est[2])
  (x[1] - smean)^2 + (sqrt(x[2]) - ssd)^2
}

.lnopt <- function(sample_mean, sample_sd) {
  opt <- optim(c(sample_mean, sample_sd), .lndiff, smean = sample_mean, ssd = sample_sd)
  
  if (opt$convergence != 0) {
    warning("Optimisation did not converge, NAs produced.")
    c("mu" = NA, "sigma" = NA)
  } else {
    c("mu" = opt$par[1], "sigma" = opt$par[2])
  }
}


# Estimate mu and sigma parameters of logit-normal from sample mean and sd
pars_logitnorm <- function(sample_mean, sample_sd) {
  if (length(sample_mean) != length(sample_sd)) stop("Parameters not same length.")
  if (any(sample_mean > 1 | sample_mean < 0)) stop("Sample mean not in [0,1]. Have you rescaled?")
  
  as.data.frame(do.call(rbind, mapply(.lnopt, sample_mean, sample_sd, SIMPLIFY = FALSE)))

}


# Wrapper for logitnorm::qlogitnorm which takes mean and sd rather than mu and sigma as parameters
qlogitnorm <- function(p, mean, sd){
  pars <- pars_logitnorm(mean, sd)
  logitnorm::qlogitnorm(p, pars$mu, pars$sigma)
}

# Wrapper for logitnorm::dlogitnorm which takes mean and sd rather than mu and sigma as parameters
dlogitnorm <- function(q, mean, sd) {
  pars <- pars_logitnorm(mean, sd)
  logitnorm::dlogitnorm(q, pars$mu, pars$sigma)
}

# Wrapper for logitnorm::plogitnorm which takes mean and sd rather than mu and sigma as parameters
plogitnorm <- function(q, mean, sd) {
  pars <- pars_logitnorm(mean, sd)
  logitnorm::plogitnorm(q, pars$mu, pars$sigma)
}

