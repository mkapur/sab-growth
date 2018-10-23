#include <TMB.hpp>
template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(Length_cm);
  DATA_VECTOR(Age);
  DATA_VECTOR(DES);
  DATA_INTEGER(nStrata);
  
  // things to estimate (per stratum)
  PARAMETER_VECTOR(t0); // t0 can be negative
  PARAMETER_VECTOR(log_k);
  PARAMETER_VECTOR(log_Linf);
  PARAMETER(sigma0);
  PARAMETER(sigma1);
  
  // exponentiate params
  vector<Type> k(nStrata);
  k = exp(log_k);  
  
  vector<Type> Linf(nStrata);
  Linf = exp(log_Linf);

  Type Sigma = 0;
  Type yfit = 0;
  vector<Type> ypreds(Age.rows());
  Type obj_fun = 0.0;
  Scalar idx;
  for(int i = 0; i < Age.rows(); i++){
    idx = DES(i);
    ypreds(i) = Linf(1)*(1-exp(-k(1)*(Age(i) - t0(1))));
    // Sigma = sigma0*pow(ypreds(i),sigma1); // Francis 1988
    // obj_fun -= log(dnorm(Length_cm(i),  ypreds(i), Sigma, true));
  }
  
        
  REPORT(ypreds);
  REPORT(idx);
  ADREPORT(Linf);
  ADREPORT(k);
  ADREPORT(t0);
  
  return(obj_fun);
}