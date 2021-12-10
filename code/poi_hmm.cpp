#include <TMB.hpp>
#include "../functions/utils.cpp"

// Likelihood for a poisson hidden markov model. 
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // Data
  DATA_VECTOR(x);          // timeseries vector
  DATA_INTEGER(m);         // Number of states m
  
  // Parameters
  PARAMETER_VECTOR(tlambda);     // conditional log_lambdas's
  PARAMETER_VECTOR(tgamma);      // m(m-1) working parameters of TPM
  
  // Uncomment only when using a non-stationary distribution
  //PARAMETER_VECTOR(tdelta);    // transformed stationary distribution,
  
  // Transform working parameters to natural parameters:
  vector<Type> lambda = tlambda.exp();
  matrix<Type> gamma = gamma_w2n(m, tgamma);
  
  // Construct stationary distribution
  vector<Type> delta = stat_dist(m, gamma);
  // If using a non-stationary distribution, use this instead
  //vector<Type> delta = delta_w2n(m, tdelta);
  
  // Get number of timesteps (n)
  int n = x.size();
  
  // Evaluate conditional distribution: Put conditional
  // probabilities of observed x in n times m matrix
  // (one column for each state, one row for each datapoint):
  matrix<Type> emission_probs(n, m);
  matrix<Type> row1vec(1, m);
  row1vec.setOnes();
  for (int i = 0; i < n; i++) {
    if (x[i] != x[i]) { // f != f returns true if and only if f is NaN. 
      // Replace missing values (NA in R, NaN in C++) with 1
      emission_probs.row(i) = row1vec;
    }
    else {
      emission_probs.row(i) = dpois(x[i], lambda, false);
    }
  }
  
  // Corresponds to Zucchini's book page 333
  matrix<Type> foo, P;
  Type mllk, sumfoo, lscale;
  
  // // No need for an "else" statement because we return
  // // the likelihood directly if m is 1, thus ending the function
  // if (m == 1) {
  //   mllk = - emission_probs.col(0).array().log().sum();
  //   
  //   // Use adreport on variables we are interested in:
  //   ADREPORT(lambda);
  //   ADREPORT(gamma);
  //   ADREPORT(delta);
  //   
  //   // Things we need for local decoding
  //   REPORT(lambda);
  //   REPORT(gamma);
  //   REPORT(delta);
  //   
  //   return mllk;
  // }
  
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
  ADREPORT(lambda);
  ADREPORT(gamma);
  ADREPORT(delta);
  
  // Variables we need for local decoding and in a convenient format
  REPORT(lambda);
  REPORT(gamma);
  REPORT(delta);
  REPORT(n);
  REPORT(emission_probs);
  REPORT(mllk);
  
  return mllk;
}
