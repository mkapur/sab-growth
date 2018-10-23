#include <TMB.hpp>
template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(Length_cm);
  DATA_VECTOR(Age);
  DATA_IVECTOR(DES); // make sure this knows it is an IVECTOR
  DATA_INTEGER(nStrata);
  
  // things to estimate (per stratum)
  PARAMETER_VECTOR(t0); // t0 can be negative
  PARAMETER_VECTOR(log_k);
  PARAMETER_VECTOR(log_Linf);
  // PARAMETER(sigma0);
  // PARAMETER(sigma1);
  PARAMETER(Sigma);
  // Type Sigma = 0.0;
  
  // exponentiate params
  vector<Type> k(nStrata);
  k = exp(log_k);  
  
  vector<Type> Linf(nStrata);
  Linf = exp(log_Linf);


  vector<Type> ypreds(Age.rows());
  Type obj_fun = 0.0;
  Type yfit = 0.0;
  for(int i = 0; i < Age.rows(); i++){
    yfit = Linf(DES(i))*(1-exp(-k(DES(i))*(Age(i) - t0(DES(i)))));
    // Sigma = sigma0*pow(yfit,sigma1); // Francis 1988
    obj_fun -= dnorm(Length_cm(i),  yfit, Sigma, true);
    ypreds(i) = yfit;
  }
  
        
  REPORT(ypreds);
  REPORT(Sigma);
  ADREPORT(Linf);
  ADREPORT(k);
  ADREPORT(t0);
  
  return(obj_fun);
}