// scoring.h: interface for the scoring class.
//
//////////////////////////////////////////////////////////////////////
#define _SCORING_H

//#define pi 3.14156

#ifndef __IOSTREAM__
#include <iostream>
#endif

#ifndef __CASSERT__
#include <cassert>
#endif

#ifndef __CSTDIO__
#include <cstdio>
#endif

#ifndef __CMATH__
#include <cmath>
#endif

#ifndef __FSTREAM__
#include <fstream>
#endif

#ifndef __STRING__
#include <string>
#endif

#ifndef __SGI_STL_VECTOR
#include <vector>
#endif

#ifndef __SGI_STL_ALGORITHM
#include <algorithm>
#endif


class scoring_params{
 public:

  double mu; //for old scoring
  double lambda; //penalty for the missing restriction site
  double nu;  //reward for each matched stretch of a map
  int delta; //how much back in the scoring matrix you can trace

  double distr_lambda; //parameters for the
  double distr_sigma;  //length distribution

  double zeta; //average number of random breaks per Kb
  double dig_p; //digestion rate

  double lo_fr_thresh; //the size thresh below which the
  //low size fr model is employed (refer to the alignment paper)
  double eta; //the std for fragments below lo_fr_thresh

 private:
  double theta;
  double phi;
  double tau;


  //below are the terms to aid fast score calculations
  double c1; // c1 = log(\sqrt{2\pi}*sigma/delta)
  double c2;//2*\sigma^2
  double c3; //log(1-dig_p);
  double c4; //c_4=c_3+log\zeta;
  double c5; //c5 = log(\sqrt{2\pi}sigma)
  
  double log_dig_p;
  double log_zeta;
  double log_theta;

  double minf;

 public:
  double fit_score_std_mult_thresh;
  double fit_score_mult_thresh;

 public:
  scoring_params();

  scoring_params(double _mu, double _lambda, double _nu, int _delta, 
		 double _distr_lambda, double _distr_sigma, 
		 double _zeta, double _dig_p, double _eta, 
		 double _lo_fr_thresh);
  
  //reference match scores: i.e. optical-reference map matching score
  //here and below:
  //x - optical map region size
  //y - reference map region size
  //k_x - optical region fragments
  //k_y - reference region fragments

  double ref_total_score_high(double x, double y, int k_x, int k_y);
  //total score for fragments >= lo_fr_thresh: includes size+site scores

  double ref_site_score(int k_x, int k_y, double y);
  //site score

  double ref_size_score_high(double x, double y, int k_x) const;

  //optical match score: i.e. optical-optical map score
  //here and below:
  //x1 - optical map #1 region size
  //x2 - optical map #2 region size
  //m1 - optical map #1 region fragments
  //m2 - optical map #2 region fragments

  double opt_total_score(double x1, double x2, int m1, int m2);
  //total score: total_score = size_score + site_score
  double opt_size_score_high(double x1, double x2, int m1, int m2);
  //this is the sizing score for comparing two optical map fragments

  double optimized_opt_total_score(double x1, double x2, int m1, 
				   int m2,double size_std_mult);

  double opt_size_score(double x1, double x2, int m1, int m2);
  double opt_size_score_low(double x1, double x2, int m1, int m2);
  //size score
  double optimized_opt_size_score(double x1, double x2, int m1, int m2, 
				  double size_std_mult);
  double optimized_opt_size_score_low(double x1, double x2, int m1, int m2, 
				  double size_std_mult);
  double optimized_opt_size_score_high(double x1, double x2, int m1, int m2,
				       double size_std_mult);

  double alt_optimized_opt_size_score(double x1, double x2, int m1, int m2, 
				      double size_std_mult);
  //this is an alternativer reference match based score for faster calc.

  double opt_site_score(int m1, int m2);
  //site score
  
  //void init(); //initializes parameters/tables

  //double matching_ref_score(double lm, double lr); 
  //the scoring for matching against the reference map
  //double matching_ref_score_mult(double x, double y, int k); 
  //the scoring for multiple matching against the reference map
  
  //double total_ref_matching_score_mult(double x, double y, int k_x, int k_y);
  //double old_total_ref_matching_score_mult(double x, double y, int k_x, int k_y);
  //double total_ref_matching_score_high1(double x, double y, int k_x, int k_y);
  //double total_ref_matching_score_low1(double x, double y, int k_x, int k_y);
  //double total_ref_matching_score1(double x, double y, int k_x, int k_y);
  //fast score to replace the old one

  //double end_matching_size_ref_score(double x, double y, int k);
  //allows to match the end positions
  //double total_end_matching_ref_score(double x, double y, int k_x, int k_y);

  //double matching_ref_score_mult_low(double x, double y, int k);
  //double matching_ref_score_mult_high(double x, double y, int k);
  //the last two functions give 2 different solutions for 2 error
  //models: one for short fragments and the other for long ones

  
  //double site_opt_match_score(int m1, int m2);
  //double site_opt_match_score(int m1, int m2, int type);
  //double matching_opt_score(double x1, double x2);
  //the scoring for matching between two random fragments
  //double matching_opt_score_mult(double x1, double x2, int m1, int m2);
  //the scoring for multiple matching between two random maps
  //double total_opt_matching_score_mult(double x1, double x2, int m1, int m2);
  //double total_opt_matching_score_mult(double x1, double x2, int m1, int m2, int type);
  //double site_ref_match_score(int map_fr, int ref_fr, double ref_y);


  //service methods
 double log_size_dens(double x, double y);
  //inline double sqr(double x) { return x * x; }

  scoring_params & operator=(const scoring_params& sp);
};
