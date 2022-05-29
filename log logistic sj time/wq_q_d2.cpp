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
  const double c;
  const double d;
  
  
public:
  Bintegrand(double k_,double t0_,double a_, double b_, double c_,double d_) : k(k_),t0(t0_), a(a_), b(b_), c(c_), d(d_){}
  double operator()(const double& x) const
  {
    double inter_screen=2  ;
    double 	tk= t0+inter_screen*(k);
    double 	tk_1=t0+inter_screen*(k-1);
    double 	Q_k=1-1/((1+ pow((tk-x)/c,-d )))  ;	
	double 	Q_k_1=1-1/((1+ pow((tk_1-x)/c,-d )))  ;	

      return  0.15* R::dlnorm(x,a,b,0)* (Q_k_1-Q_k);
  }
};

// [[Rcpp::export]]
double integrate_wq_q(double lower,double upper,double k,double t0 ,double a, double b, double c,double d)
{
  Bintegrand f(k,t0,a, b, c, d);
  double err_est;
  int err_code;
  return integrate(f, lower, upper, err_est, err_code);
}





