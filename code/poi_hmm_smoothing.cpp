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
  
  // Corresponds to (Zucchini et al., 2016, p 333)
  // temp is used in the computation of log-backward probabilities
  matrix<Type> foo, P, temp;
  Type mllk, sumfoo, lscale;
  
  // Log-forward and log-backward probabilities
  matrix<Type> lalpha(m, n);
  matrix<Type> lbeta(m, n);
  
  //// Log-forward probabilities (scaling used) and nll computation
  foo.setZero();
  // foo is a matrix where only the first column is relevant (vectors are
  // treated as columns)
  // Specifying "row(0)" or "col(0)" is optional, but it helps to remember
  // if it is a row or a column
  foo = (delta * vector<Type>(emission_probs.row(0))).matrix();
  sumfoo = foo.col(0).sum();
  lscale = log(sumfoo);
  // foo = foo.transpose(); is not reliable because of
  // aliasing, c.f https://eigen.tuxfamily.org/dox/group__TopicAliasing.html
  foo.transposeInPlace();
  // foo is now a matrix where only the first row is relevant.
  foo.row(0) /= sumfoo;
  // We can alternatively use foo.row(0).array().log()
  lalpha.col(0) = vector<Type>(foo.row(0)).log() + lscale;
  
  for (int i = 2; i <= n; i++) {
    P = emission_probs.row(i - 1);
    foo = ((foo * gamma).array() * P.array()).matrix();
    sumfoo = foo.row(0).sum();
    lscale += log(sumfoo);
    foo.row(0) /= sumfoo;
    lalpha.col(i - 1) = vector<Type>(foo.row(0)).log() + lscale;
  }
  mllk = -lscale;
  
  
  //// Log-backward probabilities (scaling used)
  lbeta.setZero();
  // "Type(m)" is used because of a similar issue at 
  // https://kaskr.github.io/adcomp/_book/Errors.html#missing-casts-for-vectorized-functions
  foo.row(0).setConstant(1 / Type(m));
  lscale = log(m);
  for (int i = n - 1; i >= 1; i--) {
    P = emission_probs.row(i);
    temp = (P.array() * foo.array()).matrix(); // temp is a row matrix
    // temp must be a column for the matrix product below
    temp.transposeInPlace();
    foo = gamma * temp;
    lbeta.col(i - 1) = foo.col(0).array().log() + lscale;
    // "foo" is a column, but must be turned into a row in order to do
    // the element-wise product above ("P.array() * foo.array()")
    foo.transposeInPlace();
    sumfoo = foo.row(0).sum();
    foo.row(0) /= sumfoo;
    lscale += log(sumfoo);
  }
  
  //// Local decoding (part 1)
  matrix<Type> smoothing_probs(m, n);
  smoothing_probs.setZero();
  Type llk = - mllk;
  for (int i = 0; i <= n - 1; i++) {
    if (x[i] != x[i]) {
      // Missing data will get a smoothing probability of NA
      smoothing_probs.col(i).setConstant(x[i]);
    } else {
      smoothing_probs.col(i) = ((lalpha.col(i) + lbeta.col(i)).array() - llk).exp();
    }
  }
  
  // Local decoding (part 2)
  vector<Type> ldecode(n);
  int col_idx;
  for (int i = 0; i < n; i++) {
    // https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
    // we don't save the result because we are not interested in the max value
    // only the index
    smoothing_probs.col(i).maxCoeff(&col_idx);
    // Columns start from 1 in R, but from 0 in C++ so we adjust to be similar
    // to R results.
    ldecode(i) = col_idx + 1;
  }
  
  // Use adreport on variables for which we want standard errors
  // Be careful, matrices are read column-wise and displayed (like 
  // the rest) as part of one vector
  ADREPORT(lambda);
  ADREPORT(gamma);
  ADREPORT(delta);
  ADREPORT(mllk);
  ADREPORT(smoothing_probs);

  // Variables we need for local decoding and in a convenient format
  REPORT(mllk);
  REPORT(lambda);
  REPORT(gamma);
  REPORT(delta);
  REPORT(n);
  REPORT(ldecode);
  
  return mllk;
}
