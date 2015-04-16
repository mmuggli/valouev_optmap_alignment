
using namespace std;

#define _GAPPEDALIGNMENT_H

#ifndef __SGI_STL_VECTOR
#include <vector>
#endif

#ifndef __SGI_STL_ALGORITHM
#include <algorithm>
#endif

#ifndef _ALIGNMENT_H
#include "alignment.h"
#endif

class alignment_seed {
 public:
  int bt;
  int br;
  int et;
  int er;
  int maxi;
  int maxj;
  double maxscore;
  int chainseed;
  int chainseedendtar;
  int chainseedendref;
  double chainseed_gap;
  double chainseed_Sscore;
  double chainseed_Tscore;

 public:
  alignment_seed() { bt = br = et = er = chainseed = maxi = maxj = -1; maxscore = 0; }
  alignment_seed & operator=(const alignment_seed &a) { bt = a.bt; br = a.br; et = a.et; er = a.er; return *this; }
    
  void setbegin(int tar_ind, int ref_ind) { bt = tar_ind;  br = ref_ind; }
  void setend(int tar_ind, int ref_ind) { et = tar_ind;  er = ref_ind; }
  void setchainedseed(int seed, int endtar, int endref, double gapscore, double Sscore, double Tscore) { 
    chainseed = seed; 
    chainseedendtar = endtar;
    chainseedendref = endref;
    chainseed_gap = gapscore;
    chainseed_Sscore = Sscore; 
    chainseed_Tscore = Tscore; 
  }
  int beg_tar_ind() { return bt; }
  int beg_ref_ind() { return br; }
  int end_tar_ind() { return et; }
  int end_ref_ind() { return er; }
  int chainedseed() { return chainseed; }
  int chainedseedendtar() { return chainseedendtar; }
  int chainedseedendref() { return chainseedendref; }
  double chainedseed_gapscore() { return chainseed_gap; }
  double chainedseed_Sscore() { return chainseed_Sscore; }
  double chainedseed_Tscore() { return chainseed_Tscore; }  
};

class gapped_alignment: public rm_alignment {
  vector<double> q;
  vector<double> r;
  int m, n, delta;
  vector< vector<int> > SeedInds;
  vector<alignment_seed> seed_list;
  vector< vector<double> > Mt;

  double align_cutoff_tscore;
  double min_seed_tscore;
  double min_seed_kb;
  double gap_dist_penalty;
  double gap_diag_penalty;
  double max_gap_kb;

 public:

  double best_seed_ind;
  bool best_seed_used;
  int max_i;
  int max_j;
  int best_chain_seed;
  vector< vector<bool> > GapEnd;

  gapped_alignment(om_read& rm, om_read& tm, scoring_params& sp) : rm_alignment(rm, tm, sp) {}
  void setparams(double align_cutoff_tscore, double min_seed_tscore, double min_seed_kb, 
		 double gap_dist_penalty, double gap_diag_penalty, double max_gap_kb);
  void align();
  void localized_align(int ref_left_ind, int ref_right_ind, int tar_left_ind, int tar_right_ind);

  void traceback(int max_i, int max_j);
  void output_alignment(ostream& out_str);

  int al_ref_beg_frag() { return ref_restr_al_sites[ref_restr_al_sites.size()-1]; }
  int al_ref_end_frag() { return ref_restr_al_sites[0]-1; }
  int al_tar_beg_frag() { return tar_restr_al_sites[tar_restr_al_sites.size()-1]; }
  int al_tar_end_frag() { return tar_restr_al_sites[0]-1; }
  double gapped_score(int tar_gap, int ref_gap);

 private:
  void init(int ref_left_ind, int ref_right_ind, int tar_left_ind, int tar_right_ind);
  void DP_main();
  void DP_chain_seeds();
  bool is_seed_long(int seed);

};

class seedcomp : public binary_function<int,int,bool> {
  const vector<alignment_seed>	&seeds;
  const bool by_ref;

 public:
  seedcomp( const vector<alignment_seed> &s, bool r ) : seeds(s), by_ref(r) {}

  bool operator()(int a, int b) const {
    if(by_ref) {
      return (seeds.at(a).br) < (seeds.at(b).br);
    } else {
      return (sqrt(seeds.at(a).bt^2 + seeds.at(a).br^2) < sqrt(seeds.at(b).bt^2 + seeds.at(b).br^2));
    }
  }
};
