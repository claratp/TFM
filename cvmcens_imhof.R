# Clara Tapia Pérez
# Master's degree thesis 
# 202206
##----------------------------------------------------------------------------##

cvmcens <- function(times, cens = rep(1, length(times)), m = 20, distr, 
                    betaLimits = c(0, 1), igumb = c(10, 10), 
                    params = list(shape = NULL, shape2 = NULL, 
                                  location = NULL, scale = NULL)){
  
# Libraries: 
require(CircularDDM);
require(CompQuadForm);
require(GofCens);
  
  
# Function to estimate the the survival function using numerical Fourier 
# inversion following Imhof (1961) using the function imhom from the R 
# package CompQuadForm
##----------------------------------------------------------------------------##
  
S_t <- function(t, lambda_j, lambda, gamma, m){
  res <- imhof(t, c(lambda_j, lambda), h = c(rep(1, length(lambda_j)), gamma))
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
p_value <- S_t(CvM, lambda_j, lambda, gamma, m)$Qq;

# Computation of the MLE of the parameters of F_0 in case they are not 
# specified
parameters <- CvM_list$Parameters;

# Results
output <- list('CvM' = CvM, 'p-value' = p_value, Distribution = distr, 
               Parameters = parameters);

return(output)
}