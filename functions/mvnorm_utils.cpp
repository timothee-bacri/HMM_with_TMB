// Function transforming working parameters in initial distribution
// to natural parameters
template<class Type>
vector<Type> delta_w2n(int m, vector<Type> tdelta) {

  vector<Type> delta(m);
  vector<Type> foo(m);
  
  if (m == 1)
    return Type(1);
  
  // set first element to one.
  // Fill in the last m - 1 elements with working parameters
  // and take exponential
  foo << Type(1), tdelta.exp();

  // normalize
  delta = foo / foo.sum();

  return delta;
}

// Function transforming the working parameters in TPM to
// natural parameters (w2n)
template<class Type>
matrix<Type> gamma_w2n(int m, vector<Type> tgamma) {

  // Construct m x m identity matrix
  matrix<Type> gamma(m, m);
  gamma.setIdentity();
  
  if (m == 1)
    return gamma;

  // Fill offdiagonal elements with working parameters column-wise:
  int idx = 0;
  for (int i = 0; i < m; i++){
    for (int j = 0; j < m; j++){
      if (j != i){
        // Fill gamma according to mapping and take exponential
        gamma(j, i) = tgamma.exp()(idx);
        idx++;
      }
    }
  }
  // Normalize each row:
  vector<Type> cs = gamma.rowwise().sum();
  for (int i = 0; i < m; i++) gamma.row(i) /= cs[i];

  return gamma;
}

// Function transforming the working parameters in the covariance matrices to
// natural parameters (w2n) (from row-wise upper triangle as a vector to
// the original covariance matrix)
template<class Type>
array<Type> sigma_w2n(int m, int p, matrix<Type> tsigma) {
  
  // Construct m matrices of size p x p (p x p x m array)
  array<Type> sigma_array(p, p, m);
  
  matrix<Type> temporary_matrix(p, p);
  matrix<Type> sigma_matrix(p, p);
  
  // Fill upper triangular elements with working parameters column-wise
  int idx;
  for (int m_idx = 0; m_idx < m; m_idx++) {
    idx = 0;
    for (int j = 0; j < p; j++) { // tsigma is filled column-wise from sigma
      for (int i = 0; i <= j; i++) {
        // Fill sigma_array according to mapping, undo the log transformation
        if (i == j) {
          sigma_array(i, j, m_idx) = exp(tsigma(idx, m_idx));
          idx++;
        } else if (i < j) {
          // Fill sigma_array according to mapping
          sigma_array(i, j, m_idx) = tsigma(idx, m_idx);
          idx++;
        } else {
          sigma_array(i, j, m_idx) = 0;
        }
      }
    }

    // Undo the Cholesky transformation
    sigma_matrix = sigma_array.col(m_idx).matrix(); // col() selects the last index, row() doesn't exist
    temporary_matrix = sigma_matrix.transpose() * sigma_matrix;
    sigma_array.col(m_idx) = temporary_matrix.array();
  }
  return sigma_array;
}

// Function computing the stationary distribution of a Markov chain
template<class Type>
vector<Type> stat_dist(int m, matrix<Type> gamma) {
  
  // Construct stationary distribution
  matrix<Type> I(m, m);
  matrix<Type> U(m, m);
  matrix<Type> row1vec(1, m);
  U = U.setOnes();
  I = I.setIdentity();
  row1vec.setOnes();
  matrix<Type> A =  I - gamma + U;
  matrix<Type> Ainv = A.inverse();
  matrix<Type> deltamat = row1vec * Ainv;
  vector<Type> delta = deltamat.row(0);
  
  return delta;
}
