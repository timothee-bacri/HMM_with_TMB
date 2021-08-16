#include <TMB.hpp> //import the TMB template

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y); // Data vector y passed from R
  DATA_VECTOR(x); // Data vector x passed from R
  
  PARAMETER(a);         // Parameter a passed from R 
  PARAMETER(b);         // Parameter b passed from R 
  PARAMETER(tsigma);    // Parameter sigma (transformed, on log-scale)
                        // passed from R
  
  // Transform tsigma back to natural scale
  Type sigma = exp(tsigma);
  
  // Declare negative log-likelihood
  Type nll = - sum(dnorm(y,
                         a + b * x,
                         sigma,
                         true));
  
  // Necessary for inference on sigma, not only tsigma
  ADREPORT(sigma);
  
  return nll;
}
