# estimate abundance from depletion electrofishing
# includes random effect for spatial replication of sites
# and fixed effects for year (or other strata)

data {
  int<lower=1> I;   # number of unique site-years
  int<lower=1> J;   # number of rows in data matrix (~ I * 3)
  int<lower=1> Y;   # number of years
  int<lower=1> S;   # number of sites
  int<lower=1> Pass[J];  # electrofishing pass (usually 1, 2 or 3)
  int<lower=0> r[J];     # number of captures during this pass
  real<lower=0> area[I]; # habitat area (or length) for computing population density
  int<lower=1> siyr[J];  # unique site-year identifier for each electrofishing pass.
                         # must be: min(siyr[])=1, max(siyr[]=I), with integers mapping to rows in area[I]
  int<lower=1> si[J];    # unique site identifier for each electrofishing pass.
  int<lower=1> yr[J];    # unique site identifier for each electrofishing pass.
}

transformed data {
  real<lower=0> larea[I];
  int<lower=1> nPass[I];
  int<lower=0> R[I];     # total number of captured fish per site-year, all passes
  for(i in 1:I) {
    R[i] <- 0;
    larea[i] <- log(area[i]);
  }
  for( j in 1:J) {
    nPass[ siyr[j]] <- max(nPass[siyr[j]], Pass[j]); # compute number of passes at each site-year
    R[ siyr[j] ] <- R[siyr[j]] + r[j];
  }
}

parameters {
  real<lower=0.01, upper=0.99> p[I];  // probability of capture
  real lrho[I]; // log of rho
}
transformed parameters {
  real<lower=0> rho[I];                 // expected density
  real<lower=0> lambda[I];        // expected number of fish in site-year
  for(i in 1:I) {
    lambda[i] <- area[i]*rho[i];  # mean density multiplied by area
    rho[i] <- exp(lrho[i]);
  }
}

model {
  real q[I];  // probability of not capture
  real plam[J]; // expected number of captures each pass

  for(i in 1:I) {
    q[i] <- 1 - p[i];   # prob of no capture on each pass
    lrho[i] ~ normal(alpha[yr ])
  }

  # estimate lambda, capture prob 
  for(j in 1:J) {
    plam[j] <- p[siyr[j] ] * pow(q[siyr[j] ], (Pass[j]-1) ) * lambda[siyr[j]];
    r[j] ~ poisson(plam[j]);
  }
}

generated quantities {
  int<lower=0> U[I];  // estimated number of uncaptured fish
  int<lower=0> N[I];  // estimated fish population size
  real<lower=0> Dn[I]; // estimated fish density (expected density + Poisson error)

  for(i in 1:I) {
      U[i] <- poisson_rng(lambda[i] * pow((1-p[i]), nPass[i]));  // uncaptured fish
      N[i] <- R[i] + U[i];      // pop estimate = total captures + estimate of uncaptured
      Dn[i] <- N[i] / area[i];
  }
}

  
  
  
  
  
  
  
  
  
  
  