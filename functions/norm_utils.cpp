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
