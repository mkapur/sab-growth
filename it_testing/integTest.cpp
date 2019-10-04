#include <TMB.hpp>


// template<class Type>
// struct univariate {
//   Type Sigma;               // Parameter in integrand
//   univariate(Type Sigma_)   // Constructor of integrand
//     : Sigma (Sigma_) {}       // Initializer list
//   Type operator()(Type x){  // Evaluate integrand
//     return selPred*exp(-(yfit- x)/(2.0*pow(Sigma*Age(i),2.0)));
//   }
// };

template <class Type> Type square(Type x){return x*x;}
template <class Type> Type selex(Type x){return 1/(1+exp( 52.976 - x)) ;}

template<class Type>
// first tell it what a 'univariate' is (via struct)
struct univariate {
  Type Sigma;               // Parameter in integrand
  univariate(Type Sigma_)   // Constructor of integrand
    : Sigma (Sigma_) {}       // Initializer list
  
  // then build this operator which will accept x from a to b
  Type operator()(Type x, Type yfit, Type Agei){  // Evaluate integrand
    return selex(yfit)*exp(-(yfit- x)/(2.0*pow(Sigma*Agei,2.0)));
  }
};


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SCALAR(a); // lower bound for x (cursive L)
  DATA_SCALAR(b); // upper bound for x (cursive L)
  
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
  Type numerator = 0.0; // used as numerator in all simulations
  Type denominator = 0.0; // denominator
  // Type selPred = 0.0; // fill with selectivity of predicted lengths
  Type yfit = 0.0;
  
  vector<Type> L1(nStrata);
  vector<Type> L2(nStrata);
  
  // first loop over all available length and generate lookup table
  // of integrals at each possible length
  // this will be given by 
  

  
  
  for(int i = 0; i < Age.rows(); i++){ // loop rows
    yfit = Linf(DES(i))*(1-exp(-k(DES(i))*(Age(i) - t0(DES(i))))); 
    Type Agei = Age(i);
    // traditional dnorm; Sel(i) is 1.0 for non-DFO regions 
    numerator = Sel(i)*dnorm(Length_cm(i),yfit,Sigma,true); 
    univariate<Type> f(Sigma, yfit, Agei);
    
    if(selType(i) == 2){
      // selPred = selex(yfit); // compute for obs based on function

      denominator = romberg::integrate(f, a, b,7,2, yfit, Agei); // return integral
    }
    if(selType(i) != 2){
      // selPred = 1.0; // coerce to 1.0 for non DFO regions
      denominator = 1.0;
    }
    
 
    
    // trunc_fun = selPred*pnorm(Length_cm(i),yfit,Sigma); // penalize by likelihood of seeing it
    
    obj_fun -= numerator/denominator; // logged already so subtract
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




