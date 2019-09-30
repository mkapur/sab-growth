#include <TMB.hpp>
template <class Type> Type square(Type x){return x*x;}
template <class Type> Type selex(Type x){return 1/(1+exp( 52.976 - x)) ;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_VECTOR(Length_cm); // observed lengths
  DATA_VECTOR(Age); // observed ages
  DATA_VECTOR(Sel); // selectivity at observed length
  DATA_IVECTOR(selType); // turns selectivity on/off for DFO data
  DATA_IVECTOR(DES); // to organize data
  DATA_INTEGER(nStrata);
  DATA_INTEGER(a2); // fixed, for reporing L1/L2
  
  // things to estimate (per stratum)
  PARAMETER_VECTOR(t0); // t0 can be negative
  PARAMETER_VECTOR(log_k);
  PARAMETER_VECTOR(log_Linf);
  PARAMETER(log_Sigma);
  Type Sigma = exp(log_Sigma); // use exponentiated here
  
  // exponentiate params
  vector<Type> k(nStrata);
  k = exp(log_k);  
  
  vector<Type> Linf(nStrata);
  Linf = exp(log_Linf);
  
  
  vector<Type> ypreds(Age.rows());
  Type obj_fun = 0.0; // gets minimized
  Type unif_obj_fun = 0.0; // used as numerator in all simulations
  Type trunc_fun = 0.0; // denominator
  Type selPred = 0.0; // fill with selectivity of predicted lengths
  Type yfit = 0.0;
  
  vector<Type> L1(nStrata);
  vector<Type> L2(nStrata);
  
  
  for(int i = 0; i < Age.rows(); i++){ // loop rows
    yfit = Linf(DES(i))*(1-exp(-k(DES(i))*(Age(i) - t0(DES(i))))); 
    
    if(selType(i) == 2){
      selPred = selex(yfit); // compute based on function
    }
    if(selType(i) != 2){
      selPred = 1.0; // coerce to 1.0 for non DFO regions
    }
    unif_obj_fun = Sel(i)*dnorm(Length_cm(i),yfit,Sigma,true); // traditional dnorm; Sel(i) is 1.0 for non-DFO regions 
    trunc_fun = selPred*pnorm(Length_cm(i),yfit,Sigma); // penalize by likelihood of seeing it

    obj_fun -= unif_obj_fun - trunc_fun; // logged already so subtract
    ypreds(i) = yfit; // store estimated length 

    // for reporting endpoints
    L1(DES(i)) =  Linf(DES(i))*(1-exp(-k(DES(i))*(0.5 - t0(DES(i))))); // change 0, 0.5
    L2(DES(i))  =  Linf(DES(i))*(1-exp(-k(DES(i))*(a2 - t0(DES(i)))));
  } // end rows
  
  REPORT(ypreds);
  ADREPORT(t0);
  ADREPORT(k);
  ADREPORT(Linf);
  ADREPORT(L1);
  ADREPORT(L2);
  ADREPORT(Sigma);
  return(obj_fun);
}




