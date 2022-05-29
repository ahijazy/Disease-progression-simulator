// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

// M(t) = E(exp(t * X)) = int exp(t * x) * f(x) dx, f(x) is the p.d.f.
class Bintegrand: public Func
{
private:
  const double k;
  const double t0;  
  const double a;
  const double b;

  
  
public:
  Bintegrand(double k_,double t0_,double a_, double b_) : k(k_),t0(t0_), a(a_), b(b_){}
  
  double operator()(const double& x) const
  {
    double inter_screen=2  ;
    double 	tk= t0+inter_screen*(k);
    return  0.15* R::dexp(x,a,0)* (R::pexp(tk-x,b,1,0));
  }
};

// [[Rcpp::export]]
double integrate_w1_q(double lower,double upper,double k,double t0 ,double a, double b)
{
  Bintegrand f(k,t0,a, b);
  double err_est;
  int err_code;
  return integrate(f, lower, upper, err_est, err_code);
}
