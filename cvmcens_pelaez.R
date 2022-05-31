# Clara Tapia Pérez
# Master's degree thesis 
# 202206
##----------------------------------------------------------------------------##

cvmcens1 <- function(times, cens = rep(1, length(times)), m = 20, distr, 
                    betaLimits = c(0, 1), igumb = c(10, 10), 
                    params = list(shape = NULL, shape2 = NULL, 
                                  location = NULL, scale = NULL)){
  
# Libraries: 
require(CircularDDM);
require(GofCens);
  

## Characteristic function of the estimated distribution of the CvM statistic
##----------------------------------------------------------------------------##
  
char_psi <- function(t, lambda_j, lambda, m, gamma) {
  res <- 1 + 0i;
  for(j in 1:m){
    res <- res * ((1 - 2i * lambda_j[j] * t)^(-1/2));
  }
  res <- res * ((1 - 2i * lambda * t)^(-gamma/2));
  return(res)
};
  
##----------------------------------------------------------------------------##
  
  
## Estimating the distribution function applying Gil-Pelaez's theorem 
## numerically integrating with the R function integrate
##----------------------------------------------------------------------------##

F_x <- function(x, lambda_j, lambda, gamma, m){
  integral_kernel <- function(t, x, lambda_j, lambda, gamma, m){
    res <- Im(exp(-t * x * 1i) * char_psi(t, lambda_j, lambda, m, gamma));
    res <- res/t;
    return(res)
  }
  res <- integrate(integral_kernel, 0, Inf, x = x, lambda_j = lambda_j, 
                   lambda =  lambda, gamma = gamma, m = m)$value;
  res <- 1/2 - (1/pi) * res;
  return(res)
};
  
##----------------------------------------------------------------------------##
  

# Sample size
n <- length(times);   

# Proportion of censored observations 
prop_cens <- 1 - sum(cens)/n;
  
# Estimation of beta 
beta <- prop_cens/(1 - prop_cens); 
  
# Beta must be strictly lower than 2
if(beta >= 2) 
  stop('The proportion of censroded observations must be under 66%')

# Beta must be strictly greater than 0
if(beta <= 0) 
  stop('The data is not censored')

# Order of the Bessel function
order <- (1 + beta)/(2 - beta);    
  
# Computation of the m eigenvalues 
lambda_j <- (1/(besselzero(order, m, 1)^2)) * ((2/(2 - beta))^2);
  
# Computation of lambda and gamma 
mean_psi <- 1/(3 * (2 - beta));   
var_psi <- (2/((1 + beta)^2)) * ((beta - 1)/((2 - beta) * (5 - beta)) + 1/9); 
lambda <- (var_psi - 2 * sum(lambda_j^2))/(2 * (mean_psi - sum(lambda_j)));
gamma <- (mean_psi - sum(lambda_j))/lambda; 
  
# Applying the gofcens function
CvM_list <- gofcens(times, cens, distr, 
                    betaLimits = betaLimits, igumb = igumb, 
                    params = params);
  
# Computation of the Cramér-von Mises test statistic
CvM <- CvM_list$`Test statistics`[[2]];
  
# Computation of the p-value
p_value <- 1 - F_x(CvM, lambda_j, lambda, gamma, m);
  
# Computation of the MLE of the parameters of F_0 in case they are not 
# specified
parameters <- CvM_list$Parameters;
  
# Results
output <- list('CvM' = CvM, 'p-value' = p_value, Distribution = distr, 
                 Parameters = parameters);
  
return(output)
}