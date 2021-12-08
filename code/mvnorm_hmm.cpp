#include <TMB.hpp>
#include "../functions/mvnorm_utils.cpp"


// Likelihood for a poisson hidden markov model. 
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_MATRIX(x);         // timeseries array (n rows * p cols)
  DATA_INTEGER(m);        // Number of states m
  
  // Parameters
  PARAMETER_MATRIX(tmu);         // mp conditional mu's (matrix: p rows, m columns)
  PARAMETER_MATRIX(tsigma);      // mp(p+1)/2 working parameters of covariance matrices (matrix: p(p+1)/2 rows, m columns)
  PARAMETER_VECTOR(tgamma);     // m(m-1) working parameters of TPM (vector: m*m-m columns)
  
  // Uncomment only when using a non-stationary distribution
  //PARAMETER_VECTOR(tdelta);    // m-1 working parameters (vector)
  
  // Load namespace which contains the multivariate distributions
  using namespace density;
  
  // Get number of covariates (p)
  int p = x.cols();
  
  // Transform working parameters to natural parameters:
  matrix<Type> mu = tmu;
  matrix<Type> gamma = gamma_w2n(m, tgamma);
  array<Type> sigma = sigma_w2n(m, p, tsigma); // Construct m matrices of size p x p (array: p x p x m)
  
  // Construct stationary distribution
  vector<Type> delta = stat_dist(m, gamma);
  // If using a non-stationary distribution, use this instead
  //vector<Type> delta = delta_w2n(m, tdelta);
  
  // Get number of timesteps (n)
  int n = x.rows();
  
  // Evaluate conditional distribution: Put conditional
  // probabilities of observed x in n times m matrix
  // (one column for each state, one row for each datapoint):
  matrix<Type> emission_probs(n, m);
  // emission_probs.setOnes();
  // matrix<Type> row1vec(1, m * p);
  // row1vec.setOnes();
  matrix<Type> sigma_m(p, p);
  vector<Type> residual_vec(p);
  Type nll = 0;
  
  // Evaluate and store the conditional distributions column-wise
  bool NA_appears = false;
  for (int m_idx = 0; m_idx < m; m_idx++) {
    MVNORM_t<Type> neg_log_dmvnorm(sigma.col(m_idx).matrix());

    for (int i = 0; i < n; i++) {
      // Replace missing values (NA in R, NaN in C++) with 1
      NA_appears = false;
      for (int p_idx = 0; p_idx < p; p_idx++) {
        if (x(i, p_idx) != x(i, p_idx)) { // f != f returns true if and only if f is NaN.
          NA_appears = true;
        }
      }

      if (NA_appears) {
        emission_probs(i, m_idx) = 1;
      } else {
        // Fill the emission probability matrix
        residual_vec = vector<Type>(x.row(i));
        residual_vec -= vector<Type>(mu.col(m_idx));
        nll = neg_log_dmvnorm(residual_vec);
        // MVNORM_t returns the negative log-likelihood, we only want the likelihood
        emission_probs(i, m_idx) = exp(-nll);
      }
    }
  }
  
  // Corresponds to Zucchini's book page 333
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
