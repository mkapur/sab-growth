#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // load data; each column is a strata
  DATA_MATRIX(Length_cm); //  may have truncated data
  DATA_MATRIX(Age);
  
  // DATA_MATRIX(len_unif); // full dataset
  // DATA_MATRIX(age_unif); // full dataset
  
  DATA_VECTOR(minSel); // obtained via 95%CI from raw data
  DATA_VECTOR(maxSel);
  DATA_INTEGER(nStrata);
  DATA_INTEGER(selType);

  int n = 500; // for now, need to change this to ,size() or something
  
  // things to estimate (one for each strata)
  PARAMETER_VECTOR(t0); // t0 can be negative 
  PARAMETER_VECTOR(log_k);
  PARAMETER_VECTOR(log_Linf);
  PARAMETER(log_Sigma);
  
  // exponentiate params
  vector<Type> k(nStrata);
  k = exp(log_k);  
  
  vector<Type> Linf(nStrata);
  Linf = exp(log_Linf);
  
  Type Sigma = exp(log_Sigma); // assuming universal sigma for now, later vectorize
  Type yfit = 0;
  matrix<Type> ypreds(n,nStrata); // only used for plotting later
  
  Type obj_fun = 0.0;
  Type unif_obj_fun = 0.0; // used as numerator in all simulations
  Type trunc_fun = 0.0; // denominator
  
  for(int j = 0; j < nStrata; j++){     // iterate strata columns
    for(int i = 0; i < n; i++){ // iterate rows
      yfit = Linf(j)*(1-exp(-k(j)*(Age(i,j) - t0(j)))); //use full data set
      
      unif_obj_fun -= dnorm(Length_cm(i,j), yfit, Sigma, true); // uniform selectivity
      // switch( selType ){
      // case 1 : // uniform selectivity(no correction)
        obj_fun = unif_obj_fun;
        ypreds(i,j) = yfit; // for plotting
      // case 2 : // minimum truncation
      //   yfit = Linf(j)*(1-exp(-k(j)*(Age(i,j) - t0(j)))); // return single value
      //   // trunc_fun = 1-pnorm(minSel(j), yfit, Sigma); //not logged yet
      //   // obj_fun -= unif_obj_fun/log(trunc_fun);
      //   obj_fun -= unif_obj_fun/log(3);
      //   ypreds(i,j) = yfit; // for plotting
      //   break;
      // case 3 : // maximum truncation
      //   yfit = Linf(j)*(1-exp(-k(j)*(Age(i,j) - t0(j)))); // return single value
      //   trunc_fun = pnorm(maxSel(i), yfit, Sigma, true);
      //   obj_fun -= unif_obj_fun/trunc_fun;
      //   break;
      // case 4 : // dome shaped selectivity
      //   trunc_fun = unif_obj_fun -  1-pnorm(minSel(i), yfit, Sigma, true) - pnorm(maxSel(i), yfit, Sigma, true);;
      // 
      //   obj_fun = unif_obj_fun/trunc_fun;
        // break;
      // }
    }
  }
  
  REPORT(ypreds);
  // REPORT(Length_cm);
  ADREPORT(Linf);
  ADREPORT(k);
  ADREPORT(t0);
  
  return(obj_fun);
}