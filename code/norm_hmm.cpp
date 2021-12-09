#include <TMB.hpp>
#include "../functions/norm_utils.cpp"

// Likelihood for a normal hidden markov model. 
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // Data
  DATA_VECTOR(x);          // timeseries vector
  DATA_INTEGER(m);         // Number of states m
  
  // Parameters
  PARAMETER_VECTOR(tmu);            // conditional means
  PARAMETER_VECTOR(tsigma);         // conditional log_sd's
  PARAMETER_VECTOR(tgamma);         // m(m-1) working parameters of TPM
  // PARAMETER_VECTOR(tdelta);      // m-1 working parameters of initial distribution
  
 
  // Transform working parameters to natural parameters:
  vector<Type> mu = tmu;
  vector<Type> sigma = tsigma.exp();
  matrix<Type> gamma = gamma_w2n(m, tgamma);
  // vector<Type> delta = delta_w2n(m, tdelta);
  
  // Construct stationary distribution 
  vector<Type> delta = stat_dist(m, gamma);
  
  // Get number of timesteps (n)
  int n = x.size();
  
  // Evaluate conditional distribution: Put conditional probabilities
  // of observed x in n times m matrix (one column for each state)
  matrix<Type> emission_probs(n, m);
  matrix<Type> row1vec(1, m);
  row1vec.setOnes();
  for (int i = 0; i < n; i++) {
    if (x[i] != x[i]) { // f != f returns true if and only if f is NaN. 
      // Replace missing values (NA in R, NaN in C++) with 1
      emission_probs.row(i) = row1vec;
    }
    else {
      emission_probs.row(i) = dnorm(x(i), mu, sigma, false);
    }
  }
  
  // Corresponds to the book page 333
  matrix<Type> foo, P;
  Type mllk, sumfoo, lscale;

  foo = (delta * vector<Type>(emission_probs.row(0))).matrix();
  sumfoo = foo.sum();
  lscale = log(sumfoo);
  foo.transposeInPlace();
  foo /= sumfoo;
  for (int i = 2; i <= n; i++) {
    P = emission_probs.row(i - 1);
    foo = ((foo * gamma).array() * P.array()).matrix();
    sumfoo = foo.sum();
    lscale += log(sumfoo);
    foo /= sumfoo;
  }
  mllk = -lscale;
  
  // Use adreport on variables for which we want standard errors
  ADREPORT(mu);
  ADREPORT(sigma);
  ADREPORT(gamma);
  ADREPORT(delta);
  
  // Variables we need for local decoding and in a convenient format
  REPORT(mu);
  REPORT(sigma);
  REPORT(gamma);
  REPORT(delta);
  REPORT(n);
  // REPORT(emission_probs);
  REPORT(mllk);
  
  return mllk;
}
