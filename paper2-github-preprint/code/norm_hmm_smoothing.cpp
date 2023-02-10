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
  // DATA_INTEGER does not allow for NA, but DATA_SCALAR does, so we use it
  DATA_SCALAR(start_row); // The crop will start with this row index
  DATA_SCALAR(start_col); // The crop will start with this column index
  DATA_SCALAR(nb_rows); // The crop will take this many rows
  DATA_SCALAR(nb_cols); // The crop will take this many columns
  // TODO:
  // Use DATA_UPDATE to adjust the submatrix dimension without reestimating
  // c.f. https://kaskr.github.io/adcomp/group__macros.html#ga3a3716228d46fbe32e23856b2f53be8d

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
  
  // Are any of these smoothing matrix parameters NA? If yes, set them to the maximum
  int int_start_row, int_start_col, int_nb_rows, int_nb_cols;
  if (start_row != start_row or start_col != start_col or nb_rows != nb_rows or nb_cols != nb_cols) {
    int_start_row = 0;
    int_start_col = 0;
    int_nb_rows = m;
    int_nb_cols = n;
  } else {
    int_start_row = CppAD::Integer(start_row);
    int_start_col = CppAD::Integer(start_col);
    int_nb_rows = CppAD::Integer(nb_rows);
    int_nb_cols = CppAD::Integer(nb_cols);
  }
  
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
  // temp is used in the computation of log-backward probabilities
  matrix<Type> foo, P, temp;
  Type mllk, sumfoo, lscale;
  
  // Log-forward and log-backward probabilities
  matrix<Type> lalpha(m, n);
  matrix<Type> lbeta(m, n);
  
  //// Log-forward probabilities (scaling used) and nll computation
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
  Type llk = - mllk;
  for (int i = 0; i <= n - 1; i++) {
    if (x[i] != x[i]) {
      // Missing data will get a smoothing probability of NA
      smoothing_probs.col(i).setConstant(x[i]);
    } else {
      smoothing_probs.col(i) = ((lalpha.col(i) + lbeta.col(i)).array() - llk).exp();
    }
  }
  // If there is a large amount of values in smoothing_probs due to e.g. a large
  // amount of data, retrieving submatrices is faster
  // https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
  // block starting at (start_row, start_col) taking nb_rows rows and nb_cols cols
  // If any submatrix parameters are NA, return the entire matrix
  // f != f returns true if and only if f is NaN (NA in R, NaN in C++).
  matrix<Type> truncated_smoothing_probs;
  truncated_smoothing_probs = smoothing_probs.block(int_start_row,
                                                    int_start_col,
                                                    int_nb_rows,
                                                    int_nb_cols);
  
  // Local decoding (part 2)
  vector<Type> ldecode(int_nb_cols);
  int col_idx;
  for (int i = int_start_row; i < int_start_row + int_nb_cols; i++) {
    // If the data is missing, set the state to -1
    if (x[i] != x[i]) {
      ldecode(i) = x[i];
    } else {
      // https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
      // we don't save the result because we are not interested in the max value
      // only the index
      truncated_smoothing_probs.col(i).maxCoeff(&col_idx);
      // Columns start from 1 in R, but from 0 in C++ so we adjust to be similar
      // to R results.
      ldecode(i) = col_idx + 1;
    }
  }
  
  // Use adreport on variables for which we want standard errors
  // Be careful, matrices are read column-wise and displayed (like 
  // the rest) as part of one vector
  ADREPORT(mu);
  ADREPORT(sigma);
  ADREPORT(gamma);
  ADREPORT(delta);
  ADREPORT(mllk);
  ADREPORT(truncated_smoothing_probs);

  // Variables we need for local decoding and in a convenient format
  REPORT(mu);
  REPORT(sigma);
  REPORT(gamma);
  REPORT(delta);
  REPORT(n);
  // REPORT(emission_probs);
  // REPORT(lalpha);
  // REPORT(lbeta);
  // REPORT(smoothing_probs);
  // REPORT(truncated_smoothing_probs);
  REPORT(ldecode);
  
  return mllk;
}
