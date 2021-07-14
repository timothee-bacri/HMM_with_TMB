template<class Type>
matrix<Type> function_example(matrix<Type> mat_example) {
  
  // This function doesn't do anything meaningful
  matrix<Type> mat(2, 3);
  mat.setOnes();
  mat.row(1) << 5, 5, 5;
  mat(0, 2) = mat.row(1).sum();
  
  return mat;
}

template<class Type>
Type logistic(Type tsigma) {
  
  return 1 / (1 + exp (-tsigma));
}
