#include "msfl.h"

#define pi 3.14156

int int_max(int x, int y){
  if(x>y) return x;
  else return y;
}
int int_min(int x, int y){
  if(x<y) return x;
  else return y;
}
double double_max(double x, double y){
  if(x>y) return x;
  else return y;
}
double double_min(double x, double y){
  if(x<y) return x;
  else return y;
}

double square(double x){
  return x*x;
}
double sqr(double x){
  return x*x;
}

double Z(double x){
  //normal probability density function
  //Z(x) = \frac{1}{2\pi} exp[-\frac{x^2}{2}]
  return (1.0/sqrt(2.0*Pi))*exp(-square(x)/2.0);
}

double P(double x){
  //area under the normal density curve to the left of x
  //$P(x) = \frac{1}{2\pi}\int\limits_{-\infty}^x e^{-t^2/2}dt$
  //error < 7.5*10^{-8}
  double p =0.2316419;
  double b1 = 0.319381530;
  double b2 = -0.356563782;
  double b3 = 1.781477937;
  double b4 = -1.821255978;
  double b5 = 1.330274429;

  double t = 1/(1+p);

  double res = 1-Z(x)*(b1*t+b2*t*t+b3*t*t*t+b4*t*t*t*t+b5*t*t*t*t*t);
  return res;
}

double Q(double x){
  //area under the normal density curve to the right of x
  //$P(x) = \frac{1}{2\pi}\int\limits_{0}^\infty e^{-t^2/2}dt$
  //error < 7.5*10^{-8}

  return 1-P(x);
}

double I(double x, double a, double b){
  //incomplete Beta function
  //$I_x(a,b) = \frac{1}{B(a,b)}\int\limits_0^x t^{\a-1}(1-t)^{b-1}dt$
  //error < 5*10^{-3}
}

double A(double t, int nu){
  assert(nu>=1);
  double res;
  double theta = atan(t/sqrt((double)nu));

  if(nu==1){
    res = (2.0/Pi)*theta;
    return res;
  }
  if(nu%2==1){
    double sum=0;
    double cur_term = 1;
    for(int i=1; i<=(nu-1)/2; i++){
      if(i==1){
	cur_term = cos(theta);
	sum += cur_term;
      }
      if(i==2){
	cur_term = cos(theta)*cos(theta)*cos(theta)*2.0/3.0;
	sum += cur_term;
      }
      if(i>2){
	cur_term *= cos(theta)*cos(theta)*((i-1)*2)/(2*(i-1)+1);
	sum += cur_term;
      }
    }
    res = (2/Pi)*(theta+sin(theta)*sum);
    return res;
  }
  else{
    double sum=0;
    double cur_term = 1;
    for(int i=1; i<=(nu-2)/2+1; i++){
      if(i==1){
	cur_term = 1;
	sum += cur_term;
      }
      else{
	cur_term *= cos(theta)*cos(theta)*((i-1)*2-1)/(2*(i-1));
	sum += cur_term;
      }
    }
    res = sin(theta)*sum;
    return res;
  }  
}

double t_prob(double t, int nu){
  //the area under the Student t - distr with nu degrees of freedom
  //between points t and -t
  assert(nu>0);
  return A(t,nu);
}

double Gamma(int n){
  //calculates factorial
  assert(n>=1);
  double max_prec_n = 11;
  if(n<=max_prec_n){
    switch(n){
    case 1: return 1;
    case 2: return 1;
    case 3: return 2;
    case 4: return 6;
    case 5: return 24;
    case 6: return 120;
    case 7: return 720;
    case 8: return 5040;
    case 9: return 40320;
    case 10: return 362880;
    case 11: return 3628800;
    }
  }
  assert(n<=11);
}

double logGamma(int n){
  //calculates log factorial
  assert(n>=1);
 
  int max_prec_n = 11;

  if(n<=max_prec_n){
    switch(n){
    case 0: assert(false);
    case 1: return 0;
    case 2: return 0;
    case 3: return 0.693147181;
    case 4: return 1.791759469;
    case 5: return 3.17805383;
    case 6: return 4.787491743;
    case 7: return 6.579251212;
    case 8: return 8.525161361;
    case 9: return 10.6046029;
    case 10: return 12.80182748;
    case 11: return 15.10441257;
    }
  }
  else{ 
    double sum = 0;
    int i;
    for(i=2; i<n; i++){
      double di;
      di = (double)i;
      
      sum += log(di);
    }
    return sum;
  }
}

double log_factorial(int n){
  return logGamma(n+1);
}
double log_poiss_pr(int x, double lambda){
  assert(lambda > 0);
  assert(x >= 0);

  double log_pr;
  log_pr = -lambda + x*log(lambda)-logGamma(x+1);

  assert(log_pr <= 0);

  return log_pr;
}
double poiss_pr(int x, double lambda){
  return exp(log_poiss_pr(x, lambda));
}

double poiss_p_value(int x, double lambda){
  double max_iterations = 15;
  int it = 0;
  double epsilon = 0.000001;
  //precision
  double cur_diff = 1;
  double p_value = 0;
  
  int cur_x = x;

  double cur_pr = 0;
  double last_pr = 10;
  while (cur_diff > epsilon){
    //while(it<= max_iterations){
    it++;
    cur_pr = poiss_pr(cur_x, lambda);
   
    
    p_value += cur_pr;
    cur_diff = fabs(cur_pr-last_pr);
    last_pr = cur_pr;
    cur_x++;
  }

  //cout<<"iterations: "<<it<<endl;
  return p_value;
}

double log_n_choose_k(int n, int k){
  if(k==0) return 0;
  else{
    assert(k<=n);
    assert(k>=1);

    double log_n_fact = 0;
    double log_k_fact = 0;
    double log_nmink_fact = 0;
    
    for(int i=1; i<=n; i++){
      log_n_fact += log((double)i);
    }
    for(int i=1; i<=k; i++){
      log_k_fact += log((double)i);
    }
    for(int i=1; i<=n-k; i++){
      log_nmink_fact += log((double)i);
    }
    return log_n_fact - log_k_fact - log_nmink_fact;
  }
}

double right_bin_p_value(int x, double p, int n){
  assert(x>=0 && n>=0 && p>=0);
  assert(x<=n);
  assert(p<=1);

  double sum=0;
  for(int i=x; i<=n; i++){
    double cur_pr;
    double cur_log_pr = log_n_choose_k(n,i)+i*log(p)+(n-i)*log(1-p);

    cur_pr = exp(cur_log_pr);
    sum += cur_pr;
    //cout<<"cur_pr:"<<cur_pr<<" sum:"<<sum<<endl;
  }
  assert(sum <= 1.1);
  return sum;
}
double left_bin_p_value(int x, double p, int n){
  return right_bin_p_value(n-x, 1-p, n);
}


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

double average(const vector<double>& v){
  double res = 0;
  int vec_size = v.size();
  for(vector<double>::const_iterator it = v.begin(); it!= v.end(); it++){
    res += (*it)/vec_size;
  }
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
