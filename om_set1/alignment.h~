// alignment.h: interface for the alignment class.
//
//////////////////////////////////////////////////////////////////////

#define _ALIGNMENT_H

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

#ifndef _M_READ_H
#include "m_read.h"
#endif

#ifndef _SCORING_H
#include "scoring.h"
#endif

//#ifndef _SIM_TABLE_H
//#include "sim_table.h"
//#endif

class rm_alignment{
 public:
  om_read ref_map;
  om_read target_map;

 public:
  
  vector<vector< double > > S; //total score alignment matrix
  vector<vector< double > > T; //t-score matrix
  
  vector<vector< int > > fromi;
  vector<vector< int > > fromj; //two traceback matrices

  vector<double> s_scores;
  vector<double> t_scores;

  double Smax; //maximum total alignment score
  double Tmax; //t-score of the optimal alignment

  scoring_params score_pars;

  vector< int > ref_restr_al_sites; 
  //contains aligned restr. sites for ref. map
  vector< int > tar_restr_al_sites; 
  //contains aligned restr. sites for tar. map

  //vector< double > declumped_best_scores;

 public:
  rm_alignment(om_read& rm, om_read& tm, scoring_params& sp);

  double al_ref_size();
  double al_tar_size();

  void overlap_alignment();
  void optimized_overlap_alignment();

  void fit_alignment(); //does a fit into the reference map
  void optimized_fit_alignment();
  void optimized_local_ref_alignment();
  void optimized_local_alignment();

  void fast_fit_alignment();
  void localized_fit_alignment(int ref_left_ind, int ref_right_ind,
			       int tar_left_ind, int tar_right_ind);
  void fast_localized_fit_alignment(int ref_left_ind, int ref_right_ind,
				    int tar_left_ind, int tar_right_ind);

  void localized_gap_alignment(double gap_open, double gap_ext_kb,
			       int ref_left_ind, int ref_right_ind, 
			       int tar_left_ind, int tar_right_ind);

  void gap_alignment(double gap_open, double gap_ext_kb);
  
  void fit_t_score();
  void overlap_t_score();
  
  double t_score_drop();

  void output_kmers(ostream& out_str);
  void output_alignment(ostream& out_str);

  int ovlp_size();

  double fit_p_value();
  double ovlp_p_value();

  double ref_size(); //size of the fit region of the alignment

  void traceback(int max_i, int max_j);
};
