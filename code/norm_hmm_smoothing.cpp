#include <TMB.hpp>
#include "../functions/norm_utils.cpp"

// Likelihood for a poisson hidden markov model. 
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // Data
  DATA_VECTOR(x);          // timeseries vector
  DATA_INTEGER(m);         // Number of states m
  // With a large data set, there are too many smoothing probabilities to return
  // So we truncate them.
  // They are returned as a block of m rows and (size of data set) columns
  // The Eigen library lets us crop that block and return a sub-matrix
  // DATA_INTEGER(start_row); // The crop will start with this row index
  // DATA_INTEGER(start_col); // The crop will start with this column index
  // DATA_INTEGER(nb_rows); // The crop will take this many rows
  // DATA_INTEGER(nb_cols); // The crop will take this many columns

  // Parameters
  PARAMETER_VECTOR(tmu);            // conditional means
  PARAMETER_VECTOR(tsigma);         // conditional log_sd's
  PARAMETER_VECTOR(tgamma);      // m(m-1) working parameters of TPM
  
  // Uncomment only when using a non-stationary distribution
  //PARAMETER_VECTOR(tdelta);    // transformed stationary distribution,
  
  
  // Transform working parameters to natural parameters:
  vector<Type> mu = tmu;
  vector<Type> sigma = tsigma.exp();
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
      emission_probs.row(i) = dnorm(x(i), mu, sigma, false);
    }
  }
  
  // Corresponds to (Zucchini et al., 2016, p 333)
  matrix<Type> foo, P, temp; // temp is used in the computation of log-backward probabilities
  Type mllk, sumfoo, lscale;

  
  // Log-forward and log-backward probabilities
  matrix<Type> lalpha(m, n);
  matrix<Type> lbeta(m, n);
  
  //// Log-forward probabilities (scaling used) and negative log-likelihood computation /////////////////
  foo.setZero();
  // foo is a matrix where only the first column is relevant (vectors are treated as columns)
  foo = (delta * vector<Type>(emission_probs.row(0))).matrix(); // Specifying "row(0)" or "col(0)" is optional, but it helps to remember if it is a row or a column
  sumfoo = foo.col(0).sum();
  lscale = log(sumfoo);
  foo.transposeInPlace(); // foo = foo.transpose(); is not reliable because of aliasing, c.f https://eigen.tuxfamily.org/dox/group__TopicAliasing.html
  // foo is now a matrix where only the first row is relevant.
  foo.row(0) /= sumfoo;
  lalpha.col(0) = vector<Type>(foo.row(0)).log() + lscale; // We can alternatively use foo.row(0).array().log()
  
  for (int i = 2; i <= n; i++) {
    P = emission_probs.row(i - 1);
    foo = ((foo * gamma).array() * P.array()).matrix();
    sumfoo = foo.row(0).sum();
    lscale += log(sumfoo);
    foo.row(0) /= sumfoo;
    lalpha.col(i - 1) = vector<Type>(foo.row(0)).log() + lscale;
  }
  mllk = -lscale;
  
  
  //// Log-backward probabilities (scaling used) ////////////
  lbeta.setZero();
  // "Type(m)" -> Similar issue at https://kaskr.github.io/adcomp/_book/Errors.html#missing-casts-for-vectorized-functions
  foo.row(0).setConstant(1 / Type(m));
  lscale = log(m);
  for (int i = n - 1; i >= 1; i--) {
    P = emission_probs.row(i);
    temp = (P.array() * foo.array()).matrix(); // temp is a row matrix
    temp.transposeInPlace(); // temp must be a column for the matrix product below
    foo = gamma * temp;
    lbeta.col(i - 1) = foo.col(0).array().log() + lscale;
    foo.transposeInPlace(); // "foo" is a column, but must be turned into a row in order to do the element-wise product above ("P.array() * foo.array()")
    sumfoo = foo.row(0).sum();
    foo.row(0) /= sumfoo;
    lscale += log(sumfoo);
  }
  
  //// Local decoding ////////////////
  matrix<Type> stateprobs(m, n);
  stateprobs.setZero();
  Type llk = - mllk;
  for (int i = 0; i <= n - 1; i++) {
    stateprobs.col(i) = ((lalpha.col(i) + lbeta.col(i)).array() - llk).exp();
  }
  
  vector<Type> ldecode(n);
  int col_idx;
  for (int i = 0; i < n; i++) {
    // https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
    // we don't save the result because we are not interested in the max value
    // only the index
    stateprobs.col(i).maxCoeff(&col_idx);
    // Columns start from 1 in R, but from 0 in C++ so we adjust to be similar
    // to R results.
    ldecode(i) = col_idx + 1;
  }
  
  // If there is a large amount of values in stateprobs and ldecode, retrieving submatrices is faster
  // https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
  // block starting at (start_row, start_col) taking nb_rows rows and nb_cols cols
  // matrix<Type> truncated_stateprobs = stateprobs.block(start_row, start_col, nb_rows, nb_cols);
  
  // Use adreport on variables for which we want standard errors
  // Be careful, matrices are read column-wise and displayed (like the rest) as part of one vector
  ADREPORT(mu);
  ADREPORT(sigma);
  ADREPORT(gamma);
  ADREPORT(delta);
  ADREPORT(stateprobs);
  // ADREPORT(truncated_stateprobs);
  ADREPORT(mllk);

  // Variables we need for local decoding and in a convenient format
  REPORT(mu);
  REPORT(sigma);
  REPORT(gamma);
  REPORT(delta);
  REPORT(n);
  REPORT(emission_probs);
  // REPORT(lalpha);
  // REPORT(lbeta);
  REPORT(stateprobs);
  // REPORT(truncated_stateprobs);
  REPORT(ldecode);
  REPORT(mllk);
  
  return mllk;
}
