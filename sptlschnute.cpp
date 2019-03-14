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
  // PARAMETER_VECTOR(log_Linf);
  PARAMETER_VECTOR(log_Ltwo);
  PARAMETER_VECTOR(log_Lone);
  PARAMETER(log_Sigma);
  Type Sigma = exp(log_Sigma);
  
  // exponentiate params
  vector<Type> k(nStrata);
  k = exp(log_k);  
  
  vector<Type> Linf(nStrata);
  // Linf = exp(log_Linf);
  
  vector<Type> Ltwo(nStrata);
  Ltwo = exp(log_Ltwo);
  
  vector<Type> Lone(nStrata);
  Lone = exp(log_Lone);
  
  vector<Type> ypreds(Age.rows());
  Type obj_fun = 0.0;
  Type yfit = 0.0;
  
  for(int i = 0; i < Age.rows(); i++){
    Linf(DES(i)) = (Ltwo(DES(i)) - Lone(DES(i))*exp(-k(DES(i))*(15)))/(1-exp(-k(DES(i))*(15)));
    // Linf(DES(i)) = (Ltwo(DES(i)) - 50*exp(-k(DES(i))*(15)))/(1-exp(-k(DES(i))*(15)));
    yfit = Linf(DES(i))*(1-exp(-k(DES(i))*(Age(i) - t0(DES(i)))));
    ypreds(i) = yfit;
    obj_fun -= dnorm(Length_cm(i),yfit,Sigma,true);
  }

  REPORT(ypreds);
  REPORT(Linf);
  ADREPORT(t0);
  ADREPORT(log_k);
  // ADREPORT(log_Linf);
  ADREPORT(log_Ltwo);
  ADREPORT(log_Lone);
  ADREPORT(Sigma);
  return(obj_fun);
  }




