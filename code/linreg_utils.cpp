#include <TMB.hpp> //import the TMB template
#include <code/utils_linreg.cpp>

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
  
  /* This is a useless example to show how to manipulate matrices
   * in C++
   * This creates a matrix of 2 rows and 3 columns.
   * The Eigen library is used to manipulate vectors, arrays, matrices...
   */
  matrix<Type> mat_example(2, 3);
  mat_example << 1, 2, 3,
                 4, 5, 6;
  matrix<Type> mat_example2 = function_example(mat_example);
  
  // This lets us retrieve any variables in a nice format
  REPORT(mat_example);
  REPORT(mat_example2);
  
  return nll;
}
