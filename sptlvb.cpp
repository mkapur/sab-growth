#include <TMB.hpp>
template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(Length_cm);
  DATA_VECTOR(Age);
  // DATA_VECTOR(Sel); // externally generated for BC only
  DATA_IVECTOR(DES); // make sure this knows it is an IVECTOR
  // DATA_IVECTOR(REG); // for conditional truncation in BC region
  DATA_INTEGER(nStrata);
  
  // things to estimate (per stratum)
  PARAMETER_VECTOR(t0); // t0 can be negative
  PARAMETER_VECTOR(log_k);
  PARAMETER_VECTOR(log_Linf);
  // PARAMETER(log_sigma0);
  // Type sigma0 = exp(log_sigma0);
  // PARAMETER(sigma1);
  PARAMETER(log_Sigma);
  Type Sigma = exp(log_Sigma);
  
  // exponentiate params
  vector<Type> k(nStrata);
  k = exp(log_k);  
  
  vector<Type> Linf(nStrata);
  Linf = exp(log_Linf);
  

  vector<Type> ypreds(Age.rows());
  Type obj_fun = 0.0;
  Type yfit = 0.0;
  Type aic;
  
  for(int i = 0; i < Age.rows(); i++){
    yfit = Linf(DES(i))*(1-exp(-k(DES(i))*(Age(i) - t0(DES(i)))));
    ypreds(i) = yfit;
    obj_fun -= dnorm(Length_cm(i),yfit,Sigma,true);
  }
  aic = 2*(2*nStrata*3+1) - 2*obj_fun;
  
  REPORT(ypreds);
  ADREPORT(t0);
  ADREPORT(log_k);
  ADREPORT(log_Linf);
  ADREPORT(Sigma);
  REPORT(aic);
  return(obj_fun);
}

    
    
    
    