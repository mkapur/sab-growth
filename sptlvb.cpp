#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // load data; each column is a strata
  DATA_MATRIX(Length_cm);
  DATA_MATRIX(Age);
  // DATA_VECTOR(st);
  DATA_INTEGER(nStrata);
  // int n = Age.size();
  int n = 500; // for now, need to change this to ,size() or something
  
  // things to estimate (one for each strata)
  PARAMETER(dummy);
  PARAMETER_VECTOR(log_t0); 
  PARAMETER_VECTOR(log_k);
  PARAMETER_VECTOR(log_Linf);
  PARAMETER(log_Sigma);
  
  // exponentiate params
  vector<Type> t0(nStrata);
  t0 = exp(log_t0);  
  
  vector<Type> k(nStrata);
  k = exp(log_k);  
  
  vector<Type> Linf(nStrata);
  Linf = exp(log_Linf);
  
  Type Sigma = exp(log_Sigma); // assuming universal sigma for now, later vectorze

  matrix<Type> yfit(n,nStrata); 
  
  Type obj_fun = 0.0;
  
  for(int j = 0; j < nStrata-1; j++){     // iterate strata
    for(int i = 0; i < n-1; i++){ // iterate rows
      yfit(i,j) = Linf(j)*(1-exp(-k(j)*(Age(i,j) - t0(j))));
    }
    obj_fun -= dnorm(Length_cm(j), yfit(j), Sigma, true); // calc LL per strata once populated
  }
  
  REPORT(yfit);
 
  return(obj_fun);
}