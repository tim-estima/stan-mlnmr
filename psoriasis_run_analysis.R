###############################################################################
# Plaque psoriasis analysis using ML-NMR 
# Author: David Phillippo
###############################################################################

# Install required packages on first use
# pkgs <- c("tidyverse", "rstan", "shinystan", "broom", "randtoolbox", "copula",
#           "xtable", "boot", "parallel", "logitnorm")
# install.packages(pkgs)

# For information on setting up RStan, see 
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

# Load packages
library(tidyverse)
library(rstan)
library(shinystan)
library(broom)
library(randtoolbox)
library(copula)
library(xtable)
library(boot) # Bootstrap for MAIC variance estimates
library(parallel)

getwd()

# Set up Stan options
Sys.setenv(USE_CXX14 = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Check for output subdirectories
check_subdir <- function(paths) {
  for (i in seq_along(paths)) {
    if (!dir.exists(paths[i])) {
      dir.create(paths[i])
      message("Created subdirectory ", paths[i])
    }
  }
}
check_subdir(c("figures", "tables"))

#----------------------
#---- Read in data ----
#----------------------

# Treatment codings
make_trtn <- function(trtc) {
  recode(trtc,
         PBO = 1,
         IXE_Q2W = 2,
         IXE_Q4W = 3,
         ETN = 4,
         SEC_150 = 5,
         SEC_300 = 6,
         UST = 7)
}

# Individual Patient Data
uncover_dat <- 
  read_csv("/dummy_ipd.csv") %>% 
  # Add in a numeric treatment code
  mutate(trtn = make_trtn(trtc)) %>% 
  # Add in numeric study IDs
  mutate(studyn = recode(study,
                         `UNCOVER-1` = 1L,
                         `UNCOVER-2` = 2L,
                         `UNCOVER-3` = 3L)) %>% 
  select(study, studyn, trtc, trtn, everything())
  
ipd_studyn <- unique(uncover_dat$studyn)

# Aggregate data
fixture_dat <- read_csv("fixture_agd.csv") %>%
  mutate(studyn = 4L, trtn = make_trtn(trtc)) %>%
  select(study, studyn, trtc, trtn, everything())

# Set some parameters---number of IPD and AgD studies, number of treatments.
ns_ipd <- length(unique(uncover_dat$studyn))
ns_agd <- length(unique(fixture_dat$studyn))
ntrt <- length(unique(c(uncover_dat$trtn, fixture_dat$trtn)))

#--------------------------------------------------
#---- ML-NMR with generalised numerical method ----
#--------------------------------------------------

# Models with the AgD model explicitly integrated can become complex and/or intractable.
# We have proposed a generalised numerical approach, where the aggregation integral is performed using numerical integration.
# We use Quasi Monte Carlo integration, based on a Sobol' sequence, for reduced integration error over standard MC.
#
# First, we must generate a grid of integration points from the joint covariate distribution in each arm of the AgD studies.
n_X <- 5  # number of covariates
n_int <- 10000  # number of sample points for numerical integration

sobol_points <- sobol(n_int, n_X)

# Since we know nothing of the correlation structure of the covariates in the AgD trials, we will assume that this is the same as that observed in the IPD. We impose the correlation structure of the UNCOVER trials on the generated Sobol' points using a simple Gaussian copula.
uncover_X_names <- c("durnpso", "prevsys_01", "bsa", "weight", "psa_01")

# Unfortunately, weight has some missing values. Since they make up such a tiny proportion of the IPD, we will simply remove these incomplete cases for this analysis.
uncover_dat <- uncover_dat %>% mutate(is.complete = select(., !! uncover_X_names) %>% complete.cases())
uncover_dat %>% 
  group_by(study) %>% 
  summarise(sum(!is.complete), mean(!is.complete))

# Now compute the correlation matrix, and apply to the integration points using a Gaussian copula
uncover_cor <- cor(uncover_dat %>% filter(is.complete) %>% select(!! uncover_X_names), 
                   method = "spearman")

uncover_copula <- normalCopula(P2p(uncover_cor), dim = n_X, dispstr = "un")

sobol_points_cor <- as_tibble(cCopula(sobol_points, uncover_copula, inverse = TRUE)) %>% setNames(uncover_X_names)

# Use inverse CDF method to generate points with AgD covariate distributions from the Sobol' points in $[0,1]$.
source("logitnorm_utils.R")

fixture_xpoints <- fixture_dat %>%
  rowwise() %>%
  do(durnpso = qgamma(sobol_points_cor$durnpso, 
                      (.$durnpso_mean / .$durnpso_sd)^2 , 
                      .$durnpso_mean / .$durnpso_sd^2),
     prevsys = qbinom(sobol_points_cor$prevsys_01, 1, .$prevsys / 100),
     bsa = qlogitnorm(sobol_points_cor$bsa,
                      .$bsa_mean / 100, .$bsa_sd / 100) * 100,
     weight = qgamma(sobol_points_cor$weight, 
                     (.$weight_mean / .$weight_sd)^2 , 
                     .$weight_mean / .$weight_sd^2),
     psa = qbinom(sobol_points_cor$psa_01, 1, .$psa / 100))  
fixture_xpoints$ag_id[1] <- 1   
fixture_xpoints$ag_id[2] <- 2    
fixture_xpoints$ag_id[3] <- 3
fixture_xpoints$ag_id[4] <- 4


    fixture_xpoints <- fixture_xpoints %>%
     tidyr::unnest(cols = c(durnpso, prevsys, bsa, weight, psa)) %>%
     dplyr::select(ag_id, everything(), -matches("^ag_id[0-9]+$"))

# Quick plot of continuous covariates to make sure these look ok
fixture_xpoints %>%
  gather(key = .var, value = x, durnpso, bsa, weight) %>%
  ggplot(aes(x = x, colour = factor(ag_id))) + 
  geom_density() +
  facet_wrap(~.var, scales = "free")

# And binary covariate proportions
fixture_xpoints %>% 
  group_by(ag_id) %>%
  select(ag_id, prevsys, psa) %>%
  summarise_all(mean)

fixture_dat %>% select(prevsys, psa)

# The Stan model takes the QR decomposition as its input, since the QR function in Stan is very memory intensive. We need to compute the thin QR decomposition of the design matrix for the centred covariates here in R, then pass the result to Stan.
# 
# First, we need to calculate a global mean value for each continuous covariate, to use for centring.

# Function to get global mean of a covariate from IPD and AgD datasets
global_mean <- function(v,    # Name of covariate, will be appended with _mean for AgD name
                        ss,   # Sample size variable in AgD
                        ipd,  # IPD dataset
                        agd,  # AgD dataset
                        na.rm = FALSE) { # Remove missing values?
  
  v <- enquo(v)
  ss <- enquo(ss)
  
  v_ipd <- pull(ipd, !! v)
  
  if (na.rm & any(is.na(v_ipd))) {
    num_na <- sum(is.na(v_ipd))
    v_ipd <- v_ipd[!is.na(v_ipd)]
    message("Removed ", num_na, " missing values in IPD.")
  } else if (any(is.na(v_ipd))) {
    warning("IPD has missing values, NA returned.", call. = FALSE)
    return(NA)
  }
  
  v_agd <- pull(agd, !! paste0(quo_name(enquo(v)), "_mean"))
  ss_agd <- pull(agd, !! ss)
  
  drop((sum(v_ipd) + ss_agd %*% v_agd) / (length(v_ipd) + sum(ss_agd)))
}


gmean_durnpso <- global_mean(durnpso, sample_size_w0,
                             filter(uncover_dat, is.complete), fixture_dat)
gmean_bsa <- global_mean(bsa, sample_size_w0,
                         filter(uncover_dat, is.complete), fixture_dat)
gmean_weight <- global_mean(weight, sample_size_w0,
                            filter(uncover_dat, is.complete), fixture_dat)


# Now we'll create a big data frame with the AgD integration points tacked on to the bottom of the IPD data set.
stan_xdat <- uncover_dat %>% filter(is.complete) %>%
  select(studyn, trtn, prevsys = prevsys_01, durnpso, bsa, weight, psa = psa_01) %>%
  bind_rows(
    fixture_dat %>% transmute(ag_id = 1:n(), studyn, trtn) %>%
      full_join(fixture_xpoints, by = "ag_id")
  ) %>%
  mutate(study = as.factor(studyn), 
         trt = as.factor(trtn), 
         trtclass = recode_factor(trtn,
                                  "1" = 1,  # Placebo
                                  "2" = 2, "3" = 2, "5" = 2, "6" = 2,  # IL-* blockers
                                  "4" = 3),  # TNFa blocker
         prevsys = prevsys,
         durnpso = (durnpso - gmean_durnpso),
         bsa = (bsa - gmean_bsa),
         weight = (weight - gmean_weight),
         psa = psa)


# Now we can simply use the `model.matrix` function to create the model matrix with a formula. To be compatible with the Stan code, the resulting model matrix columns should be in the order: study baselines, treatment parameters, PVs, EM interactions.
X_all <- model.matrix(~ -1 + study + trt + 
                        prevsys + durnpso + bsa + weight + psa +
                        (prevsys + durnpso + bsa + weight + psa):trtclass, 
                      data = stan_xdat)

# Then get the thin QR decomposition
X_all_qr <- qr(X_all)
X_all_Q <- qr.Q(X_all_qr) * sqrt(nrow(X_all) - 1)
X_all_R <- qr.R(X_all_qr)[, sort.list(X_all_qr$pivot)] / sqrt(nrow(X_all) - 1)
X_all_R_inv <- solve(X_all_R)

# Save a lookup table of the parameters for better formatting later
X_parnames <- data_frame(model_term = colnames(X_all),
                         stan_term = c(paste0("beta0[", 1:(ns_ipd + ns_agd), "]"),
                                       paste0("gamma[", 1:(ntrt - 1), "]"),
                                       paste0("beta1[", 1:n_X, "]"),
                                       paste0("beta2[", 1:(2 * n_X), "]")))

# Construct the data list for Stan
pasi75_standat <- list(
  # Constants
  ns_ipd = ns_ipd,
  ns_agd = ns_agd,
  ni_ipd = nrow(filter(uncover_dat, is.complete)),
  ni_agd = nrow(fixture_dat),
  nt = ntrt,
  nint = n_int,
  nPV = n_X,
  nEM = 2 * n_X,
  int_thin = 100,
  # IPD
  y = filter(uncover_dat, is.complete) %>% pull(pasi75_w12_nri_01),
  # AgD
  ag_n = fixture_dat$pasi75_n,
  ag_r = fixture_dat$pasi75_r,
  # QR decomposition
  Q = X_all_Q,
  R_inv = X_all_R_inv)

# Now run Stan
pasi75_stan <- stan("ML-NMR_binomial_probit_twoparbin.stan",
  data = pasi75_standat,
  pars = c("beta0", "beta1", "beta2", "gamma", "nprime", "pprime",
           "p_bar_cum", "p2_bar_cum",
           "log_lik", "resdev", "r_hat", "lp__"),
  iter = 100,
  chains = 1,
  init_r = 0.5)
