// https://kaskr.github.io/adcomp/register_atomic_8cpp-example.html#a3
// Similar to example 'adaptive_integration' using CppAD Romberg integration. REGISTER_ATOMIC is used to reduce tape size.
#include <TMB.hpp>
template <class Type> Type selex(Type x){return 1/(1+exp( 52.976 - x)) ;}


template<class Type>
  struct univariate {
    Type Linf, k, mu, t0, Sigma; // variables in integrand,  
    Type operator() (Type u) {
      Type yfit = 0;
      Type selPred = 0;
      Type ans = 0;
      yfit = Linf*(1-exp(-k*(mu - t0))); // generate predicted length
      selPred = selex(yfit);
      ans = selPred * exp(-(yfit - u)/(2*pow(Sigma*mu,2.0))); // u is what is plugged a and b
      
      // yfit = Linf(DES(i))*(1-exp(-k(DES(i))*(Age(i) - t0(DES(i))))); 
      return ans;
    }
  };
/*
  This function evaluates the marginal density of x where
u     ~ Normal( mu, sd^2 )
x | u ~ Binom ( n , plogis(u) )
*/
  template<class Type>
  vector<Type> mydenom(vector<Type> input) {
    Type log_Linf  = input[0], n  = input[1];         // Data
    Type mu = input[2], sd = input[3];         // Parameters
    univariate<Type> f = {log_Linf, n, mu, sd};
    Type a = -5, b = 5;
    vector<Type> res(1);
    res[0] = romberg::integrate(f, a, b);
    // res[0] /= exp(lgamma(n+1) - (lgamma(x+1) + lgamma(n-x+1)));
    return res;
  }
REGISTER_ATOMIC(mydenom)
template<class Type>
  Type objective_function<Type>::operator() ()
{
  

  DATA_VECTOR(t0);
  DATA_MATRIX(A);
  PARAMETER_VECTOR(b);
  vector<Type> mu = A * b;
  
  PARAMETER(logsd);
  Type sd = exp(logsd);
  
  PARAMETER_VECTOR(log_Linf);
  vector<Type> Linf(log_Linf.size()); //nstrata
  Linf = exp(log_Linf);
  
  PARAMETER_VECTOR(log_k);
  vector<Type> k(log_k.size());
  k = exp(log_k);  
  
  Type ans = 0;
  Type tiny = 0.0; // Set to 1e-12 for robustness
  for(int i=0; i < log_Linf.size(); i++) {
    vector<Type> input(5);
    input << Linf(i), k(i), mu(i), t0(i), sd;
    ans -= log( mydenom(input)[0] + tiny );
  }
  return ans;
  }