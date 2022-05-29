// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

// M(t) = E(exp(t * X)) = int exp(t * x) * f(x) dx, f(x) is the p.d.f.
class Mintegrand: public Func
{
private:
  const double z;
  const double k;
  const double t0;  
  const double a;
  const double b;
  
  
public:
  Mintegrand(double z_,double k_,double t0_,double a_, double b_) : z(z_),k(k_),t0(t0_), a(a_), b(b_){}
  
  double operator()(const double& x) const
  {
    return  0.15* R::dexp(x,a,0)* (R::dexp(z-x,b,0));
  }
};

// [[Rcpp::export]]
double dens_w1_q(double z,double lower,double upper,double k,double t0 ,double a, double b)
{
  Mintegrand f(z,k,t0,a, b );
  double err_est;
  int err_code;
  return integrate(f, lower, upper, err_est, err_code);
}
