// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

// M(t) = E(exp(t * X)) = int exp(t * x) * f(x) dx, f(x) is the p.d.f.
class Mintegrand: public Func
{
private:
const double k;
const double t0;  
const double a;
const double b;
const double c;
const double d;


public:
Mintegrand(double k_,double t0_,double a_, double b_, double c_,double d_) : k(k_),t0(t0_), a(a_), b(b_), c(c_), d(d_){}

double operator()(const double& x) const
{

 double inter_screen=1  ;
 double tk_1= t0+inter_screen*(k);
 double Q=1-1/((1+ pow((tk_1-x)/c,-d )))  ;	
return  0.15* R::dlnorm(x,a,b,0)* (1-Q);
}
};

// [[Rcpp::export]]
double integrate_w1_q(double lower,double upper,double k,double t0 ,double a, double b, double c,double d)
{
Mintegrand f(k,t0,a, b, c, d);
double err_est;
int err_code;
return integrate(f, lower, upper, err_est, err_code);
}





