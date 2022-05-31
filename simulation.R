# Clara Tapia Pérez
# Master's degree thesis 
# 202206
##----------------------------------------------------------------------------##


##----------------------------------------------------------------------------##
## Libraries
##----------------------------------------------------------------------------##

library(Rcpp);
library(Bessel);
library(CircularDDM);
library(survival);
library(fitdistrplus);
library(ggplot2);
library(actuar);
library(GofCens);
library(actuar);
library(GoFKernel);
library(CompQuadForm);
library(survsim);


##----------------------------------------------------------------------------##
## Function definitions
##----------------------------------------------------------------------------##

# The function cvmcens from the file cvmcens_imhof.R must be loaded


##----------------------------------------------------------------------------##
## Simulation Study
##----------------------------------------------------------------------------##

m <- 20; # Value m for the estimator of C:
n <- 500; # Sample size for each simulation 
N <- 1000; # Number of simulations for each scenario 
alpha <- 0.05; # Significance level 


# (1) Survival times with Weibull distribution of scale 1 and shape 1
# (exponential distribution of rate 1)
# ---------------------------------------------------------------------------- #

# (1.1) 50 % of censored data 
# ---------------------------------------------------------------------------- #
# (1.1.1) Censored times with Weibull distribution: Koziol and Green model
# ---------------------------------------------------------------------------- #

set.seed(1234);

p_111_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution using CvM 
p_111_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution using CvM
p2_111_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_111_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_111 <- simple.surv.sim(n, Inf, dist.ev = "weibull", 1, 0,
                              dist.cens = "weibull", 1, 0);
  
  times <- data_111$stop;         # Observed survival times
  delta <- data_111$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 1));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 1));
  
  # Testing the Weibull data against the log-normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-values with Imhof's method:
  p_111_w[k] = CvM_w$`p-value`;
  p2_111_w[k] = KS_w$Test[1];
  p_111_l[k] = CvM_l$`p-value`;
  p2_111_l[k] = KS_l$Test[1];
  
};

# Results: 
result_111_w <- as.numeric(p_111_w > alpha);
prop_111_w <- sum(result_111_w)/N;
result2_111_w <- as.numeric(p2_111_w > alpha);
prop2_111_w <- sum(result2_111_w)/N;
result_111_l <- as.numeric(p_111_l > alpha);
prop_111_l <- sum(result_111_l)/N;
result2_111_l <- as.numeric(p2_111_l > alpha);
prop2_111_l <- sum(result2_111_l)/N;


# (1.1.2) Censored times with Log-normal distribution
# ---------------------------------------------------------------------------- #

set.seed(1234);

p_112_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_112_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution
p2_112_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_112_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_112 <- simple.surv.sim(n, Inf,  dist.ev = "weibull", 1, 0,
                              dist.cens = "lnorm", 0.5, -0.3986672);
  
  times <- data_112$stop;         # Observed survival times
  delta <- data_112$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 1));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 1));
  
  # Testing the Weibull data against the log-normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-value with Imhof's method:
  p_112_w[k] = CvM_w$`p-value`;
  p2_112_w[k] = KS_w$Test[1];
  p_112_l[k] = CvM_l$`p-value`;
  p2_112_l[k] = KS_l$Test[1];
};

# Results:
result_112_w <- as.numeric(p_112_w > alpha);
prop_112_w <- sum(result_112_w)/N;
result2_112_w <- as.numeric(p2_112_w > alpha);
prop2_112_w <- sum(result2_112_w)/N;
result_112_l <- as.numeric(p_112_l > alpha);
prop_112_l <- sum(result_112_l)/N;
result2_112_l <- as.numeric(p2_112_l > alpha);
prop2_112_l <- sum(result2_112_l)/N;


# (1.2) 33 % of censored data 
# ---------------------------------------------------------------------------- #
# (1.2.1) Censored times with Weibull distribution: Koziol and Green model

set.seed(1234);

p_121_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_121_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution 
p2_121_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_121_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_121 <- simple.surv.sim(n, Inf, dist.ev = "weibull", 1, 0,
                              dist.cens = "weibull", 1, -log(0.492537));
  
  times <- data_121$stop;         # Observed survival times
  delta <- data_121$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 1));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 1));
  
  # Testing the Weibull data against the log-normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-values with Imhof's method:
  p_121_w[k] = CvM_w$`p-value`;
  p2_121_w[k] = KS_w$Test[1];
  p_121_l[k] = CvM_l$`p-value`;
  p2_121_l[k] = KS_l$Test[1];
};

# Results:
result_121_w <- as.numeric(p_121_w > alpha);
prop_121_w <- sum(result_121_w)/N;
result2_121_w <- as.numeric(p2_121_w > alpha);
prop2_121_w <- sum(result2_121_w)/N;
result_121_l <- as.numeric(p_121_l > alpha);
prop_121_l <- sum(result_121_l)/N;
result2_121_l <- as.numeric(p2_121_l > alpha);
prop2_121_l <- sum(result2_121_l)/N;


# (1.2.2) Censored times with Log-normal distribution -------------------------- #

set.seed(1234);

p_122_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_122_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution
p2_122_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_122_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_122 <- simple.surv.sim(n, Inf,  dist.ev = "weibull", 1, 0,
                              dist.cens = "lnorm", 0.5, 0.1224371);
  
  times <- data_122$stop;         # Observed survival times
  delta <- data_122$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 1));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 1));
  
  # Testing the Weibull data against the log-normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal"); 
  
  # Computing the p-value with Imhof's method:
  p_122_w[k] = CvM_w$`p-value`;
  p2_122_w[k] = KS_w$Test[1];
  p_122_l[k] = CvM_l$`p-value`;
  p2_122_l[k] = KS_l$Test[1];
};

# Results:
result_122_w <- as.numeric(p_122_w > alpha);
prop_122_w <- sum(result_122_w)/N;
result2_122_w <- as.numeric(p2_122_w > alpha);
prop2_122_w <- sum(result2_122_w)/N;
result_122_l <- as.numeric(p_122_l > alpha);
prop_122_l <- sum(result_122_l)/N;
result2_122_l <- as.numeric(p2_122_l > alpha);
prop2_122_l <- sum(result2_122_l)/N;


# (2) Survival times with Weibull distribution of scale 1 and shape 2
# ---------------------------------------------------------------------------- #

# (2.1) 50 % of censored data 
# ---------------------------------------------------------------------------- #
# (2.1.1) Censored times with Weibull distribution: Koziol and Green model
# ---------------------------------------------------------------------------- #

set.seed(1234);

p_211_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_211_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution 
p2_211_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_211_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_211 <- simple.surv.sim(n, Inf, dist.ev = "weibull", 2, 0,
                              dist.cens = "weibull", 2, 0);
  
  times <- data_211$stop;         # Observed survival times
  delta <- data_211$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 2));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 2));
  
  # Testing the Weibull data against the log-normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-values with Imhof's method:
  p_211_w[k] = CvM_w$`p-value`;
  p2_211_w[k] = KS_w$Test[1];
  p_211_l[k] = CvM_l$`p-value`;
  p2_211_l[k] = KS_l$Test[1];
};

# Results:
result_211_w <- as.numeric(p_211_w > alpha);
prop_211_w <- sum(result_211_w)/N;
result2_211_w <- as.numeric(p2_211_w > alpha);
prop2_211_w <- sum(result2_211_w)/N;
result_211_l <- as.numeric(p_211_l > alpha);
prop_211_l <- sum(result_211_l)/N;
result2_211_l <- as.numeric(p2_211_l > alpha);
prop2_211_l <- sum(result2_211_l)/N;


# (2.1.2) Censored times with Log-normal distribution
# ---------------------------------------------------------------------------- #

set.seed(1234);

p_212_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_212_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution
p2_212_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_212_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_212 <- simple.surv.sim(n, Inf,  dist.ev = "weibull", 2, 0,
                              dist.cens = "lnorm", 0.5, -0.2260421);
  
  times <- data_212$stop;         # Observed survival times
  delta <- data_212$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 2));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 2)); 
  
  # Testing the Weibull data against the log-normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-value with Imhof's method:
  p_212_w[k] = CvM_w$`p-value`;
  p2_212_w[k] = KS_w$Test[1];
  p_212_l[k] = CvM_l$`p-value`;
  p2_212_l[k] = KS_l$Test[1];
};

# Results:
result_212_w <- as.numeric(p_212_w > alpha);
prop_212_w <- sum(result_212_w)/N;
result2_212_w <- as.numeric(p2_212_w > alpha);
prop2_212_w <- sum(result2_212_w)/N;
result_212_l <- as.numeric(p_212_l > alpha);
prop_212_l <- sum(result_212_l)/N;
result2_212_l <- as.numeric(p2_212_l > alpha);
prop2_212_l <- sum(result2_212_l)/N;


# (2.2) 33 % of censored data 
# ---------------------------------------------------------------------------- #
# (2.2.1) Censored times with Weibull distribution: Koziol and Green model

set.seed(1234);

p_221_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_221_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution 
p2_221_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_221_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_221 <- simple.surv.sim(n, Inf, dist.ev = "weibull", 2, 0,
                              dist.cens = "weibull", 1, -log(0.4783016));
  
  times <- data_221$stop;         # Observed survival times
  delta <- data_221$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 2));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 2));
  
  # Testing the Weibull data against the log-normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-values with Imhof's method:
  p_221_w[k] = CvM_w$`p-value`;
  p2_221_w[k] = KS_w$Test[1];
  p_221_l[k] = CvM_l$`p-value`;
  p2_221_l[k] = KS_l$Test[1];
};

result_221_w <- as.numeric(p_221_w > alpha);
prop_221_w <- sum(result_221_w)/N;
result2_221_w <- as.numeric(p2_221_w > alpha);
prop2_221_w <- sum(result2_221_w)/N;
result_221_l <- as.numeric(p_221_l > alpha);
prop_221_l <- sum(result_221_l)/N;
result2_221_l <- as.numeric(p2_221_l > alpha);
prop2_221_l <- sum(result2_221_l)/N;


# (2.2.2) Censored times with Log-normal distribution ------------------------ #

set.seed(1234);

p_222_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_222_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution
p2_222_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_222_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_222 <- simple.surv.sim(n, Inf,  dist.ev = "weibull", 2, 0,
                              dist.cens = "lnorm", 0.5, 0.1015375);
  
  times <- data_222$stop;         # Observed survival times
  delta <- data_222$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 2));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 2));  
  
  # Testing the Weibull data against the log-normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-value with Imhof's method:
  p_222_w[k] = CvM_w$`p-value`;
  p2_222_w[k] = KS_w$Test[1];
  p_222_l[k] = CvM_l$`p-value`;
  p2_222_l[k] = KS_l$Test[1];
};

# Results:
result_222_w <- as.numeric(p_222_w > alpha);
prop_222_w <- sum(result_222_w)/N;
result2_222_w <- as.numeric(p2_222_w > alpha);
prop2_222_w <- sum(result2_222_w)/N;
result_222_l <- as.numeric(p_222_l > alpha);
prop_222_l <- sum(result_222_l)/N;
result2_222_l <- as.numeric(p2_222_l > alpha);
prop2_222_l <- sum(result2_222_l)/N;



# (3) Survival times with Weibull distribution of scale 1 and shape 3.5
# ---------------------------------------------------------------------------- #

# (3.1) 50 % of censored data 
# ---------------------------------------------------------------------------- #
# (3.1.1) Censored times with Weibull distribution: Koziol and Green model
# ---------------------------------------------------------------------------- #

set.seed(1234);

p_311_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_311_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution 
p2_311_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_311_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_311 <- simple.surv.sim(n, Inf, dist.ev = "weibull", 3.5, 0,
                              dist.cens = "weibull", 3.5, 0);
  
  times <- data_311$stop;         # Observed survival times
  delta <- data_311$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 3.5));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 3.5)); 
  
  # Testing the Weibull data against the log normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-values with Imhof's method:
  p_311_w[k] = CvM_w$`p-value`;
  p2_311_w[k] = KS_w$Test[1];
  p_311_l[k] = CvM_l$`p-value`;
  p2_311_l[k] = KS_l$Test[1];
};

# Results:
result_311_w <- as.numeric(p_311_w > alpha);
prop_311_w <- sum(result_311_w)/N;
result2_311_w <- as.numeric(p2_311_w > alpha);
prop2_311_w <- sum(result2_311_w)/N;
result_311_l <- as.numeric(p_311_l > alpha);
prop_311_l <- sum(result_311_l)/N;
result2_311_l <- as.numeric(p2_311_l > alpha);
prop2_311_l <- sum(result2_311_l)/N;


# (3.1.2) Censored times with Log-normal distribution
# ---------------------------------------------------------------------------- #

set.seed(1234);

p_312_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_312_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution
p2_312_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_312_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_312 <- simple.surv.sim(n, Inf,  dist.ev = "weibull", 3.5, 0,
                              dist.cens = "lnorm", 0.5, -0.144319);
  
  times <- data_312$stop;         # Observed survival times
  delta <- data_312$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 3.5));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 3.5)); 
  
  # Testing the Weibull data against the log normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-value with Imhof's method:
  p_312_w[k] = CvM_w$`p-value`;
  p2_312_w[k] = KS_w$Test[1];
  p_312_l[k] = CvM_l$`p-value`;
  p2_312_l[k] = KS_l$Test[1];
  
};

# Results: 
result_312_w <- as.numeric(p_312_w > alpha);
prop_312_w <- sum(result_312_w)/N;
result2_312_w <- as.numeric(p2_312_w > alpha);
prop2_312_w <- sum(result2_312_w)/N;
result_312_l <- as.numeric(p_312_l > alpha);
prop_312_l <- sum(result_312_l)/N;
result2_312_l <- as.numeric(p2_312_l > alpha);
prop2_312_l <- sum(result2_312_l)/N;


# (3.2) 33 % of censored data 
# ---------------------------------------------------------------------------- #
# (3.2.1) Censored times with Weibull distribution: Koziol and Green model

set.seed(1234);

p_321_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_321_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution 
p2_321_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_321_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_321 <- simple.surv.sim(n, Inf, dist.ev = "weibull", 3.5, 0,
                              dist.cens = "weibull", 1, -log(0.454399));
  
  times <- data_321$stop;         # Observed survival times
  delta <- data_321$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 3.5));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 3.5)); 
  
  # Testing the Weibull data against the log-normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-values with Imhof's method:
  p_321_w[k] = CvM_w$`p-value`;
  p2_321_w[k] = KS_w$Test[1];
  p_321_l[k] = CvM_l$`p-value`;
  p2_321_l[k] = KS_l$Test[1];
};

# Results: 
result_321_w <- as.numeric(p_321_w > alpha);
prop_321_w <- sum(result_321_w)/N;
result2_321_w <- as.numeric(p2_321_w > alpha);
prop2_321_w <- sum(result2_321_w)/N;
result_321_l <- as.numeric(p_321_l > alpha);
prop_321_l <- sum(result_321_l)/N;
result2_321_l <- as.numeric(p2_321_l > alpha);
prop2_321_l <- sum(result2_321_l)/N;


# (3.2.2) Censored times with Log-normal distribution ------------------------ #

set.seed(1234);

p_322_w <- vector("numeric", N);  # array to store the p-values comparing against
                                  # the weibull distribution
p_322_l <- vector("numeric", N);  # array to store the p-values comparing against 
                                  # the log-normal distribution
p2_322_w <- vector("numeric", N); # array to store the p-values comparing 
                                  # against the weibull distribution using KS
p2_322_l <- vector("numeric", N); # array to store the p-values comparing  
                                  # against the log-normal distribution  using KS

for(k in 1:N){
  # Sample generation:
  data_322 <- simple.surv.sim(n, Inf,  dist.ev = "weibull", 3.5, 0,
                              dist.cens = "lnorm", 0.5, 0.1190697);
  
  times <- data_322$stop;         # Observed survival times
  delta <- data_322$status;       # Status of each observation 
                                  # (0 cens 1 not cens)
  
  # Testing the Weibull data against the Weibull distribution
  # specifying the parameters:
  CvM_w <- cvmcens(times, delta, distr = "weibull", 
                   params = list(scale = 1, 
                                 shape = 3.5));
  KS_w <- KScens(times, delta, distr = "weibull", 
                 params = list(scale = 1, 
                               shape = 3.5));
  
  # Testing the Weibull data against the log normal distribution
  # not specifying the parameters:
  CvM_l <- cvmcens(times, delta, distr = "lognormal");
  KS_l <- KScens(times, delta, distr = "lognormal");
  
  # Computing the p-value with Imhof's method:
  p_322_w[k] = CvM_w$`p-value`;
  p2_322_w[k] = KS_w$Test[1];
  p_322_l[k] = CvM_l$`p-value`;
  p2_322_l[k] = KS_l$Test[1];
};

# Results: 
result_322_w <- as.numeric(p_322_w > alpha);
prop_322_w <- sum(result_322_w)/N;
result2_322_w <- as.numeric(p2_322_w > alpha);
prop2_322_w <- sum(result2_322_w)/N;
result_322_l <- as.numeric(p_322_l > alpha);
prop_322_l <- sum(result_322_l)/N;
result2_322_l <- as.numeric(p2_322_l > alpha);
prop2_322_l <- sum(result2_322_l)/N;