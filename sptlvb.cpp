#include <TMB.hpp>
template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // load data; each column is a strata
  
  // DATA_MATRIX(Length_cm); 
  // DATA_MATRIX(Age);
  DATA_ARRAY(Length_cm);
  DATA_ARRAY(Age);
  DATA_VECTOR(minSel); // obtained via 95%CI from raw data
  DATA_VECTOR(maxSel);
  DATA_INTEGER(nStrata);
  DATA_INTEGER(selType);
  DATA_VECTOR(nmat);
  
  // things to estimate (per stratum)
  PARAMETER_VECTOR(t0); // t0 can be negative 
  PARAMETER_VECTOR(log_k);
  PARAMETER_VECTOR(log_Linf);
  PARAMETER(log_Sigma);
  
  // exponentiate params
  vector<Type> k(nStrata);
  k = exp(log_k);  
  
  vector<Type> Linf(nStrata);
  Linf = exp(log_Linf);
  

  
  
  vector<Type> log_minSel(nStrata);
  log_minSel = exp(minSel);
  
  vector<Type> log_maxSel(nStrata);
  log_maxSel = exp(maxSel);
  
  // vector<Type> n(nStrata);
  Type Sigma = exp(log_Sigma); // assuming universal sigma for now, later vectorize
  Type yfit = 0;

  matrix<Type> ypreds(8481,nStrata); // used for plotting later
  
  Type obj_fun = 0.0;
  Type unif_obj_fun = 0.0; // used as numerator in all simulations
  Type trunc_fun = 0.0; // denominator
  
  for(int s = 1; s < 2; s++) // iterate sexes
  for(int j = 0; j < nStrata; j++){     // iterate strata columns
    for(int i = 0; i < nmat(j); i++){ // iterate rows, unique to strata
      yfit = Linf(j)*(1-exp(-k(j)*(Age(i,j,s) - t0(j)))); 
      unif_obj_fun = dnorm(Length_cm(i,j,s), yfit, Sigma);
      switch( selType ){
      case 1 : // uniform selectivity(no correction)
        obj_fun -= log(unif_obj_fun);
        ypreds(i,j) = yfit; 
        break;
      case 2 : // minimum truncation
        trunc_fun = 1-pnorm(minSel(j), yfit, Sigma); 
        obj_fun -= unif_obj_fun/log(trunc_fun);
        ypreds(i,j) = yfit;
        break;
      case 3 : // maximum truncation
        trunc_fun = unif_obj_fun/pnorm(maxSel(j), yfit, Sigma);
        obj_fun -= log(trunc_fun);
        ypreds(i,j) = yfit; // for plotting
        break;
      case 4 : // dome shaped selectivity
        trunc_fun = pnorm(maxSel(j), yfit, Sigma) - pnorm(minSel(j), yfit, Sigma) ;
        obj_fun = unif_obj_fun/log(trunc_fun);
        ypreds(i,j) = yfit; // for plotting
        break;
      }
    }
  }
  
  REPORT(ypreds);
  REPORT(trunc_fun)
  ADREPORT(Linf);
  ADREPORT(k);
  ADREPORT(t0);
  
  return(obj_fun);
}