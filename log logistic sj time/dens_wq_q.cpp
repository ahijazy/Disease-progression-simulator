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
  const double c;
  const double d;
  
  
public:
  Mintegrand(double z_,double k_,double t0_,double a_, double b_, double c_,double d_) : z(z_),k(k_),t0(t0_), a(a_), b(b_), c(c_), d(d_){}
  
  double operator()(const double& x) const
  {
  	double Q;
  	if(z>x)
  	{
  		Q= ((d/c)* pow(((z-x)/c),d-1))/ pow(1 +pow((z-x)/c,d),2);
	  } else{
	  	Q=0;
	  }
  	
    return  0.15* R::dlnorm(x,a,b,0)* Q;
  }
};

// [[Rcpp::export]]
double dens_wq_q(double z,double lower,double upper,double k,double t0 ,double a, double b, double c,double d)
{
  Mintegrand f(z,k,t0,a, b, c, d);
  double err_est;
  int err_code;
  return integrate(f, lower, upper, err_est, err_code);
}
