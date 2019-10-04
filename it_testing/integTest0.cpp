// Similar to example 'adaptive_integration' using CppAD Romberg integration. REGISTER_ATOMIC is used to reduce tape size.
#include <TMB.hpp>

// x = Linf, n = K, mu = Age, t = t0
// This is so it knows where to deal with things created within myDenom
template <class Type> Type selex(Type x){return 1/(1+exp( 52.976 - x)) ;}

// this is ONLY get get y fits and evaluate integral for selex type != 1. otherwise coerced.
template<class Type>
struct univariate {
  Type log_Linf, n, mu, t0, Sigma; // variables in integrand,  
  Type operator() (Type u) {
    Type yfit = 0;
    Type selPred = 0;
    Type ans = 0;
    yfit = log_Linf*(1-exp(-n*(mu - t0))); // generate predicted length
    selPred = selex(yfit);
    ans = selPred * exp(-(yfit - u)/(2*pow(Sigma*mu,2.0))); // u is what is plugged a and b
    
    // yfit = Linf(DES(i))*(1-exp(-k(DES(i))*(Age(i) - t0(DES(i))))); 
    return ans;
  }
};

// This is a full function which returns an integral given inputs and does some manipulation to it
template<class Type>
vector<Type> myDenom(vector<Type> input) {    // changed from GaussBinomial to myDenom (just the name)
  Type log_Linf  = input[0], n  = input[1];         // Data & Parameters (order mixed)
  Type mu = input[2], t0 = input[3], Sigma = input[4];        
  univariate<Type> f = {log_Linf, n, mu, t0, Sigma};  // references the function made above, and relevant pars
  Type a = 0, b = 200; // limits for integration
  vector<Type> res(1);
  res[0] = romberg::integrate(f, a, b);
  // res[0] /= exp(lgamma(n+1) - (lgamma(x+1) + lgamma(n-x+1)));
  return res;
}

// now estimate parameters
REGISTER_ATOMIC(myDenom)
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    
    // UNIVERSAL DATA
    DATA_VECTOR(Length_cm); // observed lengths
    DATA_MATRIX(A);
    PARAMETER_VECTOR(b);
    vector<Type> mu = A * b;
    
    
    PARAMETER_VECTOR(log_Linf); // log_Linf
    PARAMETER_VECTOR(t0); // t0 can be negative
    PARAMETER_VECTOR(log_k);
    PARAMETER(logSigma);
    
    // exponentiate params
    vector<Type> k(log_k.size); // nstrata
    k = exp(log_k);  
    
    vector<Type> Linf(log_Linf.size); // nstrata
    Linf = exp(log_Linf);
    
    Type Sigma = exp(logSigma);
    Type yfit = 0;
    Type ans2 = 0;
    // Type tiny = 0.0; // Set to 1e-12 for robustness
    // Type totLike = 0.0; // for quotient
    for(int i=0; i < log_Linf.size(); i++) {
      // yfit = Linf(i)*(1-exp(-k(i)*(mu(i) - t0(i)))); // generate predicted length
      
      vector<Type> input(5); // ensure this matches the number of inputs
      input << Linf(i), k(i), mu(i), t0(i), Sigma; // this is how you pass several things
      
      
      ans2 -= log( myDenom(input)[0]  ); // calc deno
    }
    return yfit;
  }