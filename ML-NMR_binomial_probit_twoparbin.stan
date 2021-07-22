data {
  // -- Constants --
  int<lower=1> ns_ipd; // number of IPD studies
  int<lower=1> ns_agd; // number of AgD studies
  int<lower=1> ni_ipd; // total number of IPD individuals
  int<lower=2> ni_agd; // total number of AgD data points
  int<lower=1> nt; // number of treatments
  int<lower=1> nint; // number of samples for numerical integration
  int<lower=0> nPV; // number of prognostic variables
  int<lower=0> nEM; // number of effect modifier *interactions* (NOT number of EMs!)
  int<lower=1> int_thin; // thinning factor for saved p_ii integration points

  // -- IPD --
  int<lower=0, upper=1> y[ni_ipd]; // binary outcome

  // -- AgD --
  int<lower=0> ag_n[ni_agd]; // outcome denominator
  int<lower=0> ag_r[ni_agd]; // outcome numerator (tried setting upper=ag_n, no such luck...)

  // The following are only needed if no PVs are included in the model (improves sampling
  // efficiency by not doing numerical integration on AgD reference trt 1 arms)

  // int<lower=1> ag_trt[ni_agd]; // treatment indicator
  // int<lower=2> ag_study[ni_agd]; // study indicator

  // -- Thin QR decomposition --
  matrix[ni_ipd + nint * ni_agd, ns_ipd + ns_agd + nPV + nEM + (nt - 1)] Q;
  matrix[ns_ipd + ns_agd + nPV + nEM + (nt - 1), ns_ipd + ns_agd + nPV + nEM + (nt - 1)] R_inv;
}
transformed data {
  // Total number of parameters and data points
  int totnpar = ns_ipd + ns_agd + nPV + nEM + (nt - 1);
  int totni = ni_ipd + nint * ni_agd;

  // Split Q matrix into IPD and AgD rows
  matrix[ni_ipd, totnpar] Q_ipd = Q[1:ni_ipd];
  matrix[nint * ni_agd, totnpar] Q_agd = Q[(ni_ipd + 1):totni];

  // nint/int_thin for numerical integration checks
  // This will give a warning about integer division, which cannot be avoided
  int n_int_thin = nint / int_thin;
}
parameters {
  // Parameters on QR scale
  vector[totnpar] beta_tilde;
}
transformed parameters {
  // -- Likelihood parameters needed later for log lik calculation --
  vector[ni_ipd] eta; // IPD linear predictor
  vector<lower=0, upper=1>[ni_ipd] theta; // IPD predicted probability
  vector[ni_agd] nprime; // AgD adjusted binomial denominator
  vector<lower=0, upper=1>[ni_agd] pprime; // AgD adjusted binomial probability


  // -- Back-transformed parameters --
  vector[totnpar] allbeta = R_inv * beta_tilde;
  // Study baselines
  vector[ns_ipd + ns_agd] beta0 = allbeta[1:(ns_ipd + ns_agd)];
  // Treatment effects
  vector[nt - 1] gamma = allbeta[(ns_ipd + ns_agd +1):(ns_ipd + ns_agd + nt - 1)];
  // Prognostic variables
  vector[nPV] beta1 = allbeta[(ns_ipd + ns_agd + nt):(ns_ipd + ns_agd + nt - 1 + nPV)];
  // EM interactions
  vector[nEM] beta2 = allbeta[(ns_ipd + ns_agd + nt + nPV):totnpar];

  // -- AgD integration --
  vector[nint * ni_agd] p_ii = Phi(Q_agd * beta_tilde);
  vector[ni_agd] p_bar;
  vector[ni_agd] p2_bar;

  // -- IPD model --
  // We define the IPD and AgD models here in the transformed parameters block,
  // as the linear predictors are required to calculate the log likelihood
  // later on. This is slightly more inefficient than defining the models
  // locally in the model block.
  eta = Q_ipd * beta_tilde;
  theta = Phi(eta);

  // -- AgD model --
  // Using the one-parameter Binomial approximation to the Poisson Binomial.
  for (i in 1:ni_agd) {
    // Uncomment if no PVs are included in the model, don't do numerical integration on reference arms

    // if (nPV == 0 && ag_trt[i] == 1) {
    //   p_bar[i] = inv_logit(beta0[ag_study[i]]);
    //   p2_bar[i] = inv_logit(beta0[ag_study[i]])^2;
    //   nprime[i] = ag_n[i];
    //   pprime[i] = p_bar[i];
    // } else {

      p_bar[i] = mean(p_ii[(1 + (i-1)*nint):(i*nint)]);
      p2_bar[i] = dot_self(p_ii[(1 + (i-1)*nint):(i*nint)]) / nint;

      // Calculate adjusted n and p
      nprime[i] = ag_n[i] * p_bar[i]^2 / p2_bar[i];
      pprime[i] = p2_bar[i] / p_bar[i];

    // }

    // Reject sample if nprime less than number of observed events (shouldn't be necessary...)
    if (nprime[i] < ag_r[i]) reject("nprime = ", nprime[i], " less than ag_r = ", ag_r[i]);
  }
}
model {
  // -- Priors --
  // These prior statements will cause Stan to raise warnings regarding transformed
  // parameters possibly needing Jacobian adjustments. These should be ignored, as
  // the transformation is entirely linear.
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);
  gamma ~ normal(0, 100);

  // -- IPD likelihood --
  y ~ bernoulli(theta);

  // -- AgD likelihood --
  // We have to hand code the log likelihood contribution for the adjusted
  // binomial here, as N is not necessarily an integer (which Stan doesn't
  // like). The following is exactly equivalent to:
  // ag_r ~ binomial(nprime, pprime);
  for (i in 1:ni_agd)
    target += lchoose(nprime[i], ag_r[i]) + lmultiply(ag_r[i], pprime[i]) + (nprime[i] - ag_r[i]) * log1m(pprime[i]);
}
generated quantities {
  // -- Log likelihood and residual deviance calculation --
  vector[ni_ipd + ni_agd] log_lik;
  vector[ni_ipd + ni_agd] resdev;

  // -- Estimate integration error --
  vector[ni_agd * n_int_thin] p_bar_cum;
  vector[ni_agd * n_int_thin] p2_bar_cum;

  // -- Predicted probabilities and numbers of events --
  vector[ni_ipd + ni_agd] p_hat;
  vector[ni_ipd + ni_agd] r_hat;

  for (i in 1:ni_ipd) {
	  p_hat[i] = theta[i];
	  r_hat[i] = theta[i];
    log_lik[i] = bernoulli_lpmf(y[i] | theta[i]);
    resdev[i] = -2 * log_lik[i];
  }

  for (i in 1:ni_agd) {
    log_lik[ni_ipd + i] = lchoose(nprime[i], ag_r[i]) + lmultiply(ag_r[i], pprime[i]) + (nprime[i] - ag_r[i]) * log1m(pprime[i]);

    r_hat[ni_ipd + i] = nprime[i] * pprime[i];
    p_hat[ni_ipd + i] = r_hat[i] / ag_n[i];

    // Approximate residual deviance for AgD, letting nprime be fixed
    resdev[ni_ipd + i] = 2 * (lmultiply(ag_r[i], ag_r[i] / r_hat[ni_ipd + i]) + lmultiply(ag_n[i] - ag_r[i], (ag_n[i] - ag_r[i]) / (ag_n[i] - r_hat[ni_ipd + i])));

	for (j in 1:n_int_thin) {
      p_bar_cum[(i-1)*n_int_thin + j] = mean(p_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]);
      p2_bar_cum[(i-1)*n_int_thin + j] = (dot_self(p_ii[(1 + (i-1)*nint):((i-1)*nint + j*int_thin)]) / (j*int_thin));
    }
  }

}
