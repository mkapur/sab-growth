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
  
  vector<Type> L1(nStrata);
  vector<Type> L2(nStrata);
  

  for(int i = 0; i < Age.rows(); i++){
    yfit = Linf(DES(i))*(1-exp(-k(DES(i))*(Age(i) - t0(DES(i)))));
    ypreds(i) = yfit;
    obj_fun -= dnorm(Length_cm(i),yfit,Sigma,true);
    L1(DES(i)) =  Linf(DES(i))*(1-exp(-k(DES(i))*(0 - t0(DES(i)))));
    L2(DES(i))  =  Linf(DES(i))*(1-exp(-k(DES(i))*(15 - t0(DES(i)))));
  }
  

  REPORT(ypreds);
  ADREPORT(t0);
  ADREPORT(k);
  ADREPORT(Linf);
  ADREPORT(L1);
  ADREPORT(L2);
  ADREPORT(Sigma);
  return(obj_fun);
}

    
    
    
    