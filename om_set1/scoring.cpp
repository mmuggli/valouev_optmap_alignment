using namespace std;

#include "scoring.h"
/*
int right_int(double x){
  return static_cast<int>(floor(x)+1);
}
double double_mmax(double x, double y){
  if (x>y) return x;
  else return y;
}

double int_pow(double x, int p){
  double res = 1;
  if (p==0) return res;
  int n;
  if(p>0) n = p;
  if(p<0) n = -p;

  int i;
  for(i=1; i<=p; i++){
    res*= x;
  }
  if (p>0) return res;
  return (1/res);
}

double poiss_pr(double zeta, int k){
  double log_pr;
  log_pr = -zeta + k*log(zeta) - logGamma(k+1);
  if(!(log_pr<=0)){
    cout<<"assertion in poiss_pr"<<endl;
    cout<<"zeta="<<zeta<<" k="<<k<<" log_pr:"<<log_pr<<endl;
  }
  assert(log_pr<=0);
  return exp(log_pr);
}
double log_poiss_pr(double zeta, int k){
  double log_pr;
  log_pr = -zeta + k*log(zeta) - logGamma(k+1);
  if(!(log_pr<=0.01)){
    cout<<"assertion in log_poiss_pr:";
    cout<<"zeta:"<<zeta<<" k:"<<k;
    cout<<" log_pr:"<<log_pr<<endl;
  }
  assert(zeta > 0);
  assert(log_pr<=0.01);
  return (log_pr);
}
double log_binomial_pr(int x, int n, double p){
  assert(p>=0 && p<=1);
  assert(x>=0 && n>=0);
  assert(x<=n);

  double res=logGamma(n)-logGamma(x)-logGamma(n-x) + x*log(p) + (n-x)*log(1-p);
  assert(res <= 0.01);
  return res;
}
double binomial_pr(int x, int n, double p){
  assert(p>=0 && p<=1);
  assert(x>=0 && n>=0);
  assert(x<=n);
  return exp(log_binomial_pr(x,n,p));
}
double BesselI(int n, double z){
  //calculates modified Bessel function of the first kind 
  //through series expansion
  assert (n>=0);
  int upper_lim = int_min(n+3,10);
  //sets the limit on the number of terms in expansion
  double res = 0;

  int k;
  for(k=0; k<upper_lim; k++){
    res+= pow(0.25*z*z,k)/(Gamma(k+1)*Gamma(n+k+1));
  }
  
  res *= int_pow(0.5*z,n);
  return res;
}
double _BesselK0(double x){
  assert(x>0);
  if(x<=0.2) return (-log(x));
  if(x>0.2 && x<=0.6) return (-2.44*x + 2.24);
  else return (exp(-x)/sqrt(2*x/pi));
}

double BesselK0(double x){
  assert(x>0);
  double res;
  if(x<=2) 
    res = -log(x/2)*BesselI(0,x) - 0.57721566 
      + 0.42278420*int_pow(x/2,2) + 0.23069756*int_pow(x/2,4)
      + 0.03488590*int_pow(x/2,6) + 0.00262698*int_pow(x/2,8);
  else{
    res = 1.25331414 - 0.07832358*(2/x)
      + 0.02189568*int_pow(2/x,2) - 0.01062446*int_pow(2/x,3)
      + 0.00587872*int_pow(2/x,4) - 0.00251540*int_pow(2/x,5)
      + 0.00053208*int_pow(2/x,6);
    res = res*exp(-x)/sqrt(x);
  }
  return res;
}
double BesselK1(double x){
  assert (x>0);
  double res;
  if(x<=2){
    res = x*log(x/2)*BesselI(1,x) + 1 + 0.15443144*int_pow(x/2,2)
      - 0.67278579*int_pow(x/2,4) - 0.18156897*int_pow(x/2,6)
      - 0.01919402*int_pow(x/2,8) - 0.00110404*int_pow(x/2,10);
    res = res/x;
  }
  else{
    res=1.25331414 + 0.23498619*(2/x)
      - 0.03655620*int_pow(2/x,2) + 0.01504268*int_pow(2/x,3)
      - 0.00780353*int_pow(2/x,4) + 0.00325614*int_pow(2/x,5)
      - 0.00068245*int_pow(2/x,6);
    res = res*exp(-x)/sqrt(x);
  }

  return res;
}
double BesselK(int n, double z){
  //calculates the values of the Bessel
  //function of the second type for the
  //integer order
  assert(z>0);

  int m;
  if(n>=0) m = n;
  else m = -n;

  if(m==0) return BesselK0(z);
  if(m==1) return BesselK1(z);

  double res = 0;
  double Knu = BesselK1(z);
  double Knumo = BesselK0(z);
  int i;
  
  double mult = 2.0;

  for(i=2; i<=m; i++){
    mult = i-1;
    res = Knumo + Knu*mult*2.0/z;
    Knumo = Knu;
    Knu = res;
  }
  return res;
}
double BesselK_semi_int(double nu, double z){
  //this function calculates the values of 
  //the Bessel function of the second type,
  //when nu = n-1/2 , where m is integer.

  assert (z>=0);

  int i=0;
  double res;
  res = 0;
  
  int n;
  if(nu>=0) n = (int)(nu+0.5);
  else n = (int)(-nu+0.5);

  for(i=0; i<n; i++){
    res += pow(2.0*z,-i)*Gamma(n+i)/(Gamma(i+1)*Gamma(n-i));
  }

  res *= sqrt(pi/(2*z))*exp(-z);
  return res;
}
double logBesselK_semi_int(double nu, double z){
  //this function calculates the values of 
  //log Bessel function of the second type,
  //when nu = n-1/2 , where m is integer.

  assert (z>=0);

  int i=0;
  double sum;
  sum = 0;// sqrt(pi/(2*z))*exp(-z);

  
  int n;
  if(nu>=0) n = (int)(nu+0.5);
  else n = (int)(-nu+0.5);

  for(i=0; i<n; i++){
    sum += pow(2.0*z,-i)*Gamma(n+i)/(Gamma(i+1)*Gamma(n-i));
  }

  double res;
  res = log(sum) + 0.5*log(pi/(2*z)) - z;

  return res;
}
double B0(double x){
  //calculates BesselK_0(z)*exp(z)*sqrt(z) to
  //avoid underflow
  assert(x>2);
  double res;

  res = 1.25331414 - 0.07832358*(2/x)
    + 0.02189568*int_pow(2/x,2) - 0.01062446*int_pow(2/x,3)
    + 0.00587872*int_pow(2/x,4) - 0.00251540*int_pow(2/x,5)
    + 0.00053208*int_pow(2/x,6);
  return res;
}
double B1(double x){
  assert (x>2);
  double res;

  res=1.25331414 + 0.23498619*(2/x)
    - 0.03655620*int_pow(2/x,2) + 0.01504268*int_pow(2/x,3)
    - 0.00780353*int_pow(2/x,4) + 0.00325614*int_pow(2/x,5)
    - 0.00068245*int_pow(2/x,6);

  return res;
}
double BK(int n, double z){
  //calculates the values of the 
  //scaled Bessel BesselK_n(z)*exp(z)*sqrt(z)
  //function of the second type for the
  //integer order
  assert(z>2);

  int m;
  if(n>=0) m = n;
  else m = -n;

  if(m==0) return B0(z);
  if(m==1) return B1(z);

  double res = 0;
  double Knu = B1(z);
  double Knumo = B0(z);
  int i;
  
  double mult = 2.0;

  for(i=2; i<=m; i++){
    mult = i-1;
    res = Knumo + Knu*mult*2.0/z;
    Knumo = Knu;
    Knu = res;
  }
  return res;
}
double new_BesselK(int n, double z){
  if (z<=2) return BesselK(n,z);
  else return BK(n,z)*exp(-z)/sqrt(z);
}
double logBesselK(int n, double z){
  assert(z>0);
  if (z>=2)return log(BK(n,z))-z-0.5*log(z);
  else return log(BesselK(n,z));
}

double norm_dens(double x, double mu, double sigma){
  double f;
  f = (1/(sqrt(2*pi)*sigma))*exp(-sqr(x-mu)/(2*sqr(sigma)));
  return f;
}
double log_norm_dens(double x, double mu, double sigma){
  return -sqr(x-mu)/(2*sqr(sigma)) - 0.5*log(2*pi*sigma*sigma);
}
double gamma_dens(double x, double beta, int k){
  double log_gamma_dens = -logGamma(k) - k*log(beta) + (k-1)*x - x/beta;
  return exp(log_gamma_dens);
}

double log_exp_dens(double x, double lambda){
  assert(x>=0);
  assert(lambda > 0);
  return -x/lambda - log(lambda);
}
*/
/*
double scoring_params::total_ref_matching_score1(double x, double y, int k_x, int k_y){
  //if((x-y)*(x-y) > fit_score_mult_thresh*y)
  //return minf;
  //else
    return total_ref_matching_score_high1(x,y,k_x,k_y);
}
*/
double scoring_params::ref_total_score_high(double x, double y, int k_x, int k_y){
  double neg_score;

  neg_score = c1 - (k_y-1)*c3 - (k_x-1)*log_zeta
    - k_x*log_theta + (k_x-1)*log(x) - (k_x-1.5)*log(y)
   + sqr(x-y)/(c2*y) - x/theta + zeta*y;

  return(-neg_score);
}

double scoring_params::optimized_opt_size_score(double x1, double x2, int m1, int m2, double size_std_mult){
  assert(dig_p > 0 && dig_p < 1);
  //double sigma = distr_sigma;
  //double _lambda = distr_lambda;
  //double tau = 1/(zeta + dig_p/distr_lambda);
  //double phi = distr_lambda/sqr(dig_p);  
  //double theta = distr_sigma/(sqrt(2.0/tau+1/sqr(distr_sigma))-1/distr_sigma);

  assert (x1> 0 && x2>0);

  //double mult = 5.0;

  if(fabs(x1-x2) > size_std_mult*distr_sigma*sqrt(double_max(x1,x2))) return -1000;

  //double dm1, dm2;
  //dm1 = (double) m1;
  //dm2 = (double) m2;

  double up_thresh = 4.0;
  if(double_max(x1,x2)<=up_thresh)
    return opt_size_score_low(x1,x2,m1,m2);
  else
    return opt_size_score_high(x1,x2,m1,m2);
}

double scoring_params::opt_size_score(double x1, double x2, int m1, int m2){
  assert(dig_p > 0 && dig_p < 1);
  assert (x1> 0 && x2>0);

  double up_thresh = 4.0;
  if(double_max(x1,x2)<=up_thresh)
    return opt_size_score_low(x1,x2,m1,m2);
  else
    return opt_size_score_high(x1,x2,m1,m2);
}


double scoring_params::optimized_opt_size_score_low(double x1, double x2, int m1, int m2, double size_std_mult){
  assert(dig_p > 0 && dig_p < 1);
  //double sigma = distr_sigma;
  //double _lambda = distr_lambda;
  //double tau = 1/(zeta + dig_p/distr_lambda);
  //double phi = distr_lambda/sqr(dig_p);  
  //double theta = distr_sigma/(sqrt(2.0/tau+1/sqr(distr_sigma))-1/distr_sigma);

  assert (x1> 0 && x2>0);

  //double mult = 5.0;

  if(fabs(x1-x2) > size_std_mult*distr_sigma*sqrt(double_max(x1,x2))) return -1000;

  double dm1, dm2;
  dm1 = (double) m1;
  dm2 = (double) m2;
  
  double score;

  score = logGamma(m1) + logGamma(m2) + (m1+m2)*log(theta)
    - log(pi*phi*sqr(distr_sigma))
    + logBesselK(0,2*sqrt((1/phi+1/sqr(distr_sigma))
			  *((x1*x1+x2*x2)/(2*sqr(distr_sigma)))))
    + (x1+x2)*(1/theta + 1/sqr(distr_sigma))
    - (m1-1)*log(x1) - (m2-1)*log(x2)+(m2-1)*log(x2)+log(0.25+0.375*double_max(x1,x2));
  return score;
 
}

double scoring_params::optimized_opt_size_score_high(double x1, double x2, int m1, int m2, double size_std_mult){
  assert(dig_p > 0 && dig_p < 1);
  //double sigma = distr_sigma;
  //double _lambda = distr_lambda;
  //double tau = 1/(zeta + dig_p/distr_lambda);
  //double phi = distr_lambda/sqr(dig_p);  
  //double theta = distr_sigma/(sqrt(2.0/tau+1/sqr(distr_sigma))-1/distr_sigma);

  assert (x1> 0 && x2>0);

  //double mult = 5.0;

  if(fabs(x1-x2) > size_std_mult*distr_sigma*sqrt(double_max(x1,x2))) return -1000;

  double dm1, dm2;
  dm1 = (double) m1;
  dm2 = (double) m2;
  
  double score;

  score = logGamma(m1) + logGamma(m2) + (m1+m2)*log(theta)
    - log(pi*phi*sqr(distr_sigma))
    + logBesselK(0,2*sqrt((1/phi+1/sqr(distr_sigma))
			  *((x1*x1+x2*x2)/(2*sqr(distr_sigma)))))
    + (x1+x2)*(1/theta + 1/sqr(distr_sigma))
    - (m1-1)*log(x1) - (m2-1)*log(x2);
  return score;
 
}

/*
double scoring_params::alt_optimized_opt_size_score(double x1, double x2, int m1, int m2, 
						    double size_std_mult){
  if(fabs(x1-x2) > size_std_mult*distr_sigma*sqrt(double_max(x1,x2))) return -1000;
  
  double s1 = ref_size_score_high(x1,x2,m1);
  double s2 = ref_size_score_high(x2,x1,m2);
  
  return 0.5*(s1+s2);
}
*/

double scoring_params::opt_size_score_high(double x1, double x2, int m1, int m2){
  assert(dig_p > 0 && dig_p < 1);
  //double sigma = distr_sigma;
  //double _lambda = distr_lambda;
  //double tau = 1/(zeta + dig_p/distr_lambda);
  //double phi = distr_lambda/sqr(dig_p);  
  //double theta = distr_sigma/(sqrt(2.0/tau+1/sqr(distr_sigma))-1/distr_sigma);

  assert (x1> 0 && x2>0);

  //double mult = 5.0;

  //if(fabs(x1-x2) > mult*distr_sigma*sqrt(double_max(x1,x2))) return -1000;

  double dm1, dm2;
  dm1 = (double) m1;
  dm2 = (double) m2;
  
  double score;

  score = logGamma(m1) + logGamma(m2) + (m1+m2)*log(theta)
    - log(pi*phi*sqr(distr_sigma))
    + logBesselK(0,2*sqrt((1/phi+1/sqr(distr_sigma))
			  *((x1*x1+x2*x2)/(2*sqr(distr_sigma)))))
    + (x1+x2)*(1/theta + 1/sqr(distr_sigma))
    - (m1-1)*log(x1) - (m2-1)*log(x2);
  return score;
 
}

double scoring_params::opt_size_score_low(double x1, double x2, int m1, int m2){
  assert(dig_p > 0 && dig_p < 1);
  //double sigma = distr_sigma;
  //double _lambda = distr_lambda;
  //double tau = 1/(zeta + dig_p/distr_lambda);
  //double phi = distr_lambda/sqr(dig_p);  
  //double theta = distr_sigma/(sqrt(2.0/tau+1/sqr(distr_sigma))-1/distr_sigma);

  assert (x1> 0 && x2>0);

  double mult = 5.0;

  //if(fabs(x1-x2) > mult*distr_sigma*sqrt(double_max(x1,x2))) return -1000;

  double dm1, dm2;
  dm1 = (double) m1;
  dm2 = (double) m2;
  
  double score;

  score = logGamma(m1) + logGamma(m2) + (m1+m2)*log(theta)
    - log(pi*phi*sqr(distr_sigma))
    + logBesselK(0,2*sqrt((1/phi+1/sqr(distr_sigma))
			  *((x1*x1+x2*x2)/(2*sqr(distr_sigma)))))
    + (x1+x2)*(1/theta + 1/sqr(distr_sigma))
    - (m1-1)*log(x1) - (m2-1)*log(x2)+log(0.25+0.375*double_max(x1,x2));
  return score;
 
}


double scoring_params::ref_size_score_high(double x, double y, int k_x) const{
  double score;
  score = 0 - c5 - (k_x-1)*log(x) + k_x*log_theta 
    + logGamma(k_x) - sqr(x-y)/(c2*y) + x/theta;
  return score;
}

/*
double scoring_params::matching_opt_score(double x1, double x2){
  assert(x1>=0 && x2>=0);

  double sigma = distr_sigma;//sqrt(0.306);
  double _lambda = distr_lambda;//8;

  if(x1==0 && x2==0) return 0;

  double num;
  double den;

  double lr;

  num = pi*exp(-(1/sigma)*sqrt(2/_lambda+1/sqr(sigma))*(x1+x2));
  den = _lambda*(2/_lambda+1/sqr(sigma))*
    BesselK0(sqrt(2.0)*sqrt(1/_lambda+1/sqr(sigma))*sqrt((x1*x1+x2*x2)/sqr(sigma)));

  lr = num/den;
  return -log(lr);
}
*/

/*
double scoring_params::matching_ref_score(double x, double y){
  double lr;

  double inf=1000;
  double _lambda = distr_lambda;//8;
  double sigma = distr_sigma;//sqrt(0.306);

  if (y==0 && x==0) return 0;
  if (y==0 && x>0) return (-inf);
  
  lr = sqrt(2*pi*y)/(_lambda*sqrt(2/_lambda+1/(sigma*sigma))) *
    exp( x*(1/(sigma*sigma)-(1/sigma)*sqrt(2/_lambda + 1/(sigma*sigma))) +
	 sqr(x-y)/(2*sigma*sigma*y));

  return (-log(lr));
}
*/
/*
double scoring_params::matching_ref_score_mult_low(double x, double y, int k){
  
  double m = (double) k;
  double loglr = (m-1)*log(x) - m*log(theta) - logGamma(k)
    + log(sqrt(2*pi)*eta)
    - x/theta + sqr(x-y)/(2*sqr(eta));
  
  return -loglr;
}
*/
/*
double scoring_params::matching_ref_score_mult_high(double x, double y, int k){

  assert(x>0 && y>0);
  double m = (double) k;
  double inf = 10000;
  
  double loglr = (m-1)*log(x) - m*log(theta) - logGamma(k)
    + log(sqrt(2*pi*y)*distr_sigma)
    - x/theta + sqr(x-y)/(2*sqr(distr_sigma)*y);
  
  if (loglr<-inf) return (-inf);
  else return -loglr;
}
*/
/*
double scoring_params::ref_total_score(double x, double y, int k){
  if(y>lo_fr_thresh) return matching_ref_score_mult_high(x,y,k);
  else return matching_ref_score_mult_low(x,y,k);
}
*/
/*
double scoring_params::total_end_matching_ref_score(double x, double y,
						    int k_x, int k_y){
  return end_matching_size_ref_score(x,y,k_x) + site_ref_match_score(k_x, k_y, y);
}
*/

double scoring_params::ref_site_score(int k_x, int k_y, double y){
  double log_alt = (k_y - 1)*log(1-dig_p) + 
    log_poiss_pr(zeta*y, k_x-1);
  double log_nul = log(1/(((double)delta)));

  return log_alt - log_nul;
}

double scoring_params::opt_site_score(int m1, int m2){
  assert(m1>=1 && m2>=1);
  return nu - lambda*(m1+m2-2);

}
/*
double scoring_params::site_opt_match_score(int m1, int m2, int type){
  assert(m1>=1 && m2>=1);
  //return nu - lambda*(m1+m2-2);
  //if(m1==1 && m2==1) return nu;
  //else return -lambda*(m1+m2-2);
  //return nu-lambda*(m1+m2-2);
  if(type == 1)
    return opt_site_match_vec[m1-1][m2-1];
  if(type == 2) return nu-lambda*(m1+m2-2);

  assert(false);
}
*/
double scoring_params::optimized_opt_total_score(double x1, double x2, int m1, int m2, double size_std_mult){
  assert(m1>0 && m2>0 && x1>0 && x2>0);
  int gap1 = m1-1;
  int gap2 = m2-1;

  double score = optimized_opt_size_score(x1,x2,m1,m2,size_std_mult)
    + opt_site_score(m1,m2);
  return score;
}

double scoring_params::opt_total_score(double x1, double x2, int m1, int m2){
  assert(m1>0 && m2>0 && x1>0 && x2>0);
  int gap1 = m1-1;
  int gap2 = m2-1;

  double score = opt_size_score(x1,x2,m1,m2)
    + opt_site_score(m1,m2);
  return score;
}
/*
double scoring_params::total_opt_matching_score_mult
(double x1, double x2, int m1, int m2, int type){
  assert(m1>0 && m2>0 && x1>0 && x2>0);
  assert(!opt_site_match_vec.empty());
  int gap1 = m1-1;
  int gap2 = m2-1;

  //double score = matching_opt_score_mult(x1,x2,m1,m2)
  //  + site_opt_match_score(m1,m2,type);
  //return score;
}
*/
/*
double scoring_params::old_total_ref_matching_score_mult
(double x, double y, int k_x, int k_y){
  return nu - mu*fabs(x-y) - lambda*(k_x+k_y-2);
}
*/
scoring_params::scoring_params(){//empty constructor
}

scoring_params::scoring_params(double _mu, double _lambda, double _nu, 
			       int _delta, double _distr_lambda, 
			       double _distr_sigma, double _zeta, 
			       double _dig_p, double _eta,
			       double _lo_fr_thresh){
  mu = _mu;
  lambda = _lambda;
  nu = _nu;
  delta = _delta;
  distr_lambda = _distr_lambda;
  distr_sigma = _distr_sigma;
  zeta = _zeta;
  dig_p = _dig_p;
  eta = _eta;
  lo_fr_thresh = _lo_fr_thresh;
  if(_dig_p < 0.5){
    cout<<"unreasonably low digestion rate"<<endl;
    assert(false);
  }

  tau = 1/(zeta + dig_p/distr_lambda);
  theta = distr_sigma/(sqrt(2.0/tau+1/sqr(distr_sigma))-1/distr_sigma);
  phi = distr_lambda/sqr(dig_p);

  c1 = log(sqrt(2*pi)*distr_sigma/delta);
  c2 = 2*sqr(distr_sigma);
  c3 = log(1-dig_p);
  c4 = c3 + log(zeta);
  c5 = log(sqrt(2*pi)*distr_sigma);

  log_zeta = log(zeta);
  log_theta = log(theta);

  minf = -1000;
  fit_score_std_mult_thresh = 16.0;
  fit_score_mult_thresh = distr_sigma*distr_sigma*fit_score_std_mult_thresh;
  
}

scoring_params & scoring_params::operator=(const scoring_params &sp){
  mu = sp.mu;
  lambda = sp.lambda;
  nu = sp.nu;
  delta = sp.delta;
  distr_lambda = sp.distr_lambda;
  distr_sigma = sp.distr_sigma;
  zeta = sp.zeta;
  dig_p = sp.dig_p;
  eta = sp.eta;
  lo_fr_thresh = sp.lo_fr_thresh;
  //opt_site_match_vec = sp.opt_site_match_vec;

  theta = sp.theta;
  phi = sp.phi;
  tau = sp.tau;

  c1 = sp.c1;
  c2 = sp.c2;
  c3 = sp.c3;
  c4 = sp.c4;
  c5 = sp.c5;

  log_dig_p = sp.log_dig_p;
  log_theta = sp.log_theta;
  log_zeta = sp.log_zeta;

  minf = sp.minf;
  fit_score_std_mult_thresh = sp.fit_score_std_mult_thresh;
  fit_score_mult_thresh = sp.fit_score_mult_thresh;

  return *this;
}

double scoring_params::log_size_dens(double x, double y){
  return  - sqr(x-y)/(2*sqr(distr_sigma)*y) - log(sqrt(2*pi*y)*distr_sigma);
}
