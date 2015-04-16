
using namespace std;

#include "alignment.cpp"
#include "gappedalignment.h"

void gapped_alignment::setparams(double align_cutoff_tscore, 
				 double min_seed_tscore, double min_seed_kb,
				 double gap_dist_penalty, double gap_diag_penalty, double max_gap_kb) {

  this->align_cutoff_tscore = align_cutoff_tscore;
  this->min_seed_tscore =  min_seed_tscore;
  this->min_seed_kb =  min_seed_kb;
  this->gap_dist_penalty =  gap_dist_penalty;
  this->gap_diag_penalty =  gap_diag_penalty;
  this->max_gap_kb = max_gap_kb;
}

void gapped_alignment::align() {

  localized_align(0, ref_map.map_read.size()-1, 0, target_map.map_read.size()-1);

}

void gapped_alignment::localized_align(int ref_left_ind, int ref_right_ind, 
				      int tar_left_ind, int tar_right_ind) {

  assert(ref_left_ind >= 0);
  assert(ref_right_ind <= ref_map.map_read.size());
  assert(ref_left_ind < ref_right_ind);

  assert(tar_left_ind >= 0);
  assert(tar_right_ind <= target_map.map_read.size());
  assert(tar_left_ind < tar_right_ind);

  init(ref_left_ind, ref_right_ind, tar_left_ind, tar_right_ind);
  DP_main();
  DP_chain_seeds();
  traceback(max_i, max_j);

}

void gapped_alignment::init(int ref_left_ind, 
			    int ref_right_ind, 
			    int tar_left_ind, 
			    int tar_right_ind){

  int i, j;
  double counter;

  m = tar_right_ind - tar_left_ind;
  n = ref_right_ind - ref_left_ind;
  delta = score_pars.delta;

  Smax = 0;
  Tmax = 0;
  max_i = 0;
  max_j = 0;
  best_seed_ind = -1;
  best_seed_used = false;
 
  S.clear();
  fromi.clear();
  fromj.clear();
  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();
  SeedInds.clear();
  seed_list.clear();

  r.push_back(0);
  counter = 0;
  for(i=ref_left_ind; i<ref_right_ind; i++){
    counter += ref_map.map_read[i];
    r.push_back(counter);
  }

  q.push_back(0);
  counter=0;
  for(j=tar_left_ind; j<tar_right_ind; j++){
    counter += target_map.map_read[j];
    q.push_back(counter);
  }

  for(i=0; i<=m; i++){
    vector<double> zs;
    vector<double> zt;
    vector<double> zm;
    vector<int> zseed;
    vector<int> mminone;
    vector<bool> zgap;

    for(j=0; j<=n; j++){
      zs.push_back(0);
      zt.push_back(0);
      zm.push_back(0);
      zseed.push_back(-1);
      mminone.push_back(-1);
      zgap.push_back(false);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    S.push_back(zs);
    T.push_back(zt);
    Mt.push_back(zm);
    SeedInds.push_back(zseed);
    GapEnd.push_back(zgap);
  }

    for(j=0; j<=n; j++){
      S[0][j] = 0;
      T[0][j] = 0;
      Mt[0][j] = 0;
    }
}

void gapped_alignment::DP_main() {

  int i, j, g, h;
  int k, b;

  double min_score_thresh = 0;//-10;
  double min_t_score_thresh = 0;//-5;

  for(i=1; i<=m; i++){
    for(j=1; j<=n; j++){
      S[i][j] = 0;
      T[i][j] = 0;

      for(g=int_max(0, i-delta); g<i; g++){
	for(h=int_max(0, j-delta); h<j; h++){
	  
	  int map_gap = i-g;
	  int ref_gap = j-h;
	  double tar_size = q[i]-q[g];    
	  double ref_size = r[j]-r[h];
	  
	  double s = S[g][h] + score_pars.ref_size_score_high
	    (tar_size, ref_size, int_max(map_gap, map_gap/*ref_gap*/));	  
	  double t = T[g][h] + score_pars.opt_site_score(map_gap, ref_gap);
	  //double t = T[g][h] + score_pars.ref_site_score(map_gap, ref_gap, ref_size);

	  //Set maximum score as next match in alignment
	  if(s >= S[i][j] && 
	     s >= min_score_thresh &&
	     t >= min_t_score_thresh &&
	     t - Mt[g][h] >= align_cutoff_tscore) {
	    fromi[i][j] = g;
	    fromj[i][j] = h;
	    S[i][j] = s;
	    T[i][j] = t;
	    Mt[i][j] = double_max(t, Mt[g][h]);

	    if(i==1430 && j==49){
	      cout<<"debug: t="<<t<<endl;
	    }
	  }
	}
      }

      int fi = fromi[i][j];
      int fj = fromj[i][j];

      //Current Alignment is good enough to be a seed
      if(T[i][j] >= min_seed_tscore ||
	 (fi >= 0 && fj >= 0 && SeedInds[fi][fj] >= 0)) {
	int seed_ind = SeedInds[fi][fj];

	//New Seed
	if(seed_ind < 0) {
	  alignment_seed newseed;
	  seed_list.push_back(newseed);
	  seed_ind = seed_list.size()-1;

	  bool end_found = false;
	  int curposi = i;
	  int curposj = j;

	  while(end_found == false){
	    int newi;
	    int newj;

	    SeedInds[curposi][curposj] = seed_ind;
      
	    newi = fromi[curposi][curposj];
	    newj = fromj[curposi][curposj];
      
	    if(newi == -1 || newj == -1){
	      end_found = true;
	    } else {
	      curposi = newi;
	      curposj = newj;
	    }
	  }

	  seed_list[seed_ind].setbegin(curposi, curposj);
	}

	//New or continued seed
	SeedInds[i][j] = seed_ind;

	if(T[i][j] > seed_list[seed_ind].maxscore) {
	  seed_list[seed_ind].setend(i, j);
	  seed_list[seed_ind].maxscore = T[i][j];
	  seed_list[seed_ind].maxi = i;
	  seed_list[seed_ind].maxj = j;
	}
	if(T[i][j] > Tmax) 
	  best_seed_ind = seed_ind;
      }

      //Overall aligned maximum is the end of the best local alignment
      if(T[i][j] > Tmax) {
	Tmax = T[i][j];
	Smax = S[i][j];
	max_i = i;
	max_j = j;
      }

    }
  }
  cout<<"Seeds: "<<seed_list.size()<<endl;
}

void gapped_alignment::DP_chain_seeds() {

  int i, j;

  //sorted lists of seeds by beg_ref_ind or beg_tar_ind
  vector<int> ref_seeds;
  vector<int> origin_seeds;
  vector<int> seed_ref_map;
  vector<int> seed_origin_map;

  for(i=0; i<seed_list.size(); i++) {
    cout<<"Debug: seed="<<i<<" t="<<T[seed_list[i].et][seed_list[i].er]<<endl;

    ref_seeds.push_back(i);
    origin_seeds.push_back(i);
    seed_ref_map.push_back(-1);
    seed_origin_map.push_back(-1);
  }

  cout<<"Debug: best seed="<<best_seed_ind<<endl;

  sort(ref_seeds.begin(), ref_seeds.end(), seedcomp(seed_list, true));
  sort(origin_seeds.begin(), origin_seeds.end(), seedcomp(seed_list, false));

  for(i=0; i<seed_list.size(); i++) {
    seed_origin_map[origin_seeds[i]] = i;
    seed_ref_map[ref_seeds[i]] = i;
  }


  //score offsets for chained seeds
  vector< vector<double> > S_offset;
  vector< vector<double> > T_offset;
  for(i=0; i<=m; i++){
    vector<double> so;
    vector<double> to;    
    for(j=0; j<=n; j++) {
      so.push_back(0);
      to.push_back(0);
    }
    S_offset.push_back(so);
    T_offset.push_back(to);
  }

  int thisseed;
  for(i=1; i<m; i++) {
    for(j=1; j<=n; j++){
      thisseed = SeedInds[i][j];

      if(thisseed >= 0 && is_seed_long(thisseed)) {
	//Find gap only for beginning of seed alignment
	if(fromi[i][j] < 0 || fromj[i][j] < 0) {
	  double max_gs = 0;
	  double max_gsS = 0;
	  double max_gsT = 0;
	  int max_gap = -1;
	  int max_gap_endtar = -1;
	  int max_gap_endref = -1;

	  int end_tar_pos_ind = seed_list[thisseed].beg_tar_ind();
	  int end_ref_seed_ind = seed_ref_map[thisseed];
	  int end_origin_seed_ind = seed_origin_map[thisseed];

	  //Set minimum ref ind as next closest seed to the origin, so even if this is too far away it will be checked
	  int min_ref_ind;
	  if(end_origin_seed_ind > 0) {
	    min_ref_ind = seed_ref_map[origin_seeds[end_origin_seed_ind - 1]];
	  } else {
	    min_ref_ind = end_ref_seed_ind;
	  }

	  //find thisseed in the ref ordered list, and go to the furthest one with the same beginning ref ind
	  while(end_ref_seed_ind+1 < ref_seeds.size() && 
		seed_list[ref_seeds[end_ref_seed_ind+1]].beg_ref_ind() == seed_list[ref_seeds[end_ref_seed_ind]].beg_ref_ind()) {
	    end_ref_seed_ind++;
	  }

	  //start with one to the left of thisseed in ref order
	  int cur_ref_seed_ind = end_ref_seed_ind - 1;
	  while(cur_ref_seed_ind >= 0 && 
		(cur_ref_seed_ind >= min_ref_ind ||
		 seed_list[thisseed].beg_ref_ind() <= seed_list[ref_seeds[cur_ref_seed_ind]].end_ref_ind() ||
		 (r[seed_list[thisseed].beg_ref_ind()] - r[seed_list[ref_seeds[cur_ref_seed_ind]].end_ref_ind()]) <= max_gap_kb)
		) {

	    int gapseed = ref_seeds[cur_ref_seed_ind];

	    //if also less than thisseed in tar order, and is long enough then use gapseed
	    if(cur_ref_seed_ind != end_ref_seed_ind &&
	       is_seed_long(gapseed) &&
	       seed_list[gapseed].beg_tar_ind() <= end_tar_pos_ind &&
	       (cur_ref_seed_ind >= min_ref_ind ||
		seed_list[thisseed].beg_tar_ind() <= seed_list[gapseed].end_tar_ind() ||
		(q[seed_list[thisseed].beg_tar_ind()] - q[seed_list[gapseed].end_tar_ind()]) <= max_gap_kb)
	       ) {

	      int gbt = seed_list[gapseed].beg_tar_ind();
	      int gbr = seed_list[gapseed].beg_ref_ind();
	
	      if(gbt < i && gbr < j) {
		//find end point of gapseed
		int get = seed_list[gapseed].end_tar_ind();
		int ger = seed_list[gapseed].end_ref_ind();
		while(get > i || ger > j) {
		  int newget = fromi[get][ger];
		  int newger = fromj[get][ger];
		  get = newget;
		  ger = newger;
		  assert(get >= 0 && ger >= 0);
		}

		//Calc gapped score and find maximum
		double gs = gapped_score(i - get, j - ger);
		double gsS = S[get][ger];
		double gsT = T[get][ger];
		if(gsT + gs > max_gsT + max_gs) {
		  max_gs = gs;
		  max_gsS = gsS;
		  max_gsT = gsT;
		  max_gap = gapseed;
		  max_gap_endtar = get;
		  max_gap_endref = ger;
		}
	      }
	    }

	    cur_ref_seed_ind--;
	  }

	  fromi[i][j] = max_gap_endtar;
	  fromj[i][j] = max_gap_endref;

	  if(max_gap_endtar >= 0 && max_gap_endref >= 0) {
	    GapEnd[i][j] = true;
	    S_offset[i][j] = max_gsS;
	    T_offset[i][j] = max_gsT + max_gs;
	  }

	} else {
	  //Continue with the same seed alignment
	  S_offset[i][j] = S_offset[fromi[i][j]][fromj[i][j]];
	  T_offset[i][j] = T_offset[fromi[i][j]][fromj[i][j]];

	}
	
	//Update score for either continued alignment or gap  
	S[i][j] = S[i][j] + S_offset[i][j];
	T[i][j] = T[i][j] + T_offset[i][j];


	//Overall aligned maximum is the end of the best chained alignment
	if(T[i][j] > Tmax) {
	  Tmax = T[i][j];
	  Smax = S[i][j];
	  max_i = i;
	  max_j = j;
	}
      }
    }
  }
}

void gapped_alignment::traceback(int max_i, int max_j) {
  bool end_found = false;
  int curposi = max_i;
  int curposj = max_j;

  //curposi = 1430+1;
  //curposj = 49+1;

  tar_restr_al_sites.clear();
  ref_restr_al_sites.clear();

  while(end_found == false){
    tar_restr_al_sites.push_back(curposi);
    ref_restr_al_sites.push_back(curposj);

    if(SeedInds[curposi][curposj] == best_seed_ind)
      best_seed_used = true;
      
    int newi;
    int newj;
    
    newi = fromi[curposi][curposj];
    newj = fromj[curposi][curposj];
      
    curposi = newi;
    curposj = newj;
      
    if(curposi == -1 || curposj == -1){
      end_found = true;
    }

  }

  //Set Smax and Tmax for this alignment
  Smax = S[max_i][max_j];
  Tmax = T[max_i][max_j];

}  

double gapped_alignment::gapped_score(int tar_gap, int ref_gap) {
  return (gap_dist_penalty * (tar_gap + ref_gap)) + (gap_diag_penalty * abs(tar_gap - ref_gap)) ;
}

bool gapped_alignment::is_seed_long(int seed) {
  return ((r[seed_list[seed].end_ref_ind()] - r[seed_list[seed].beg_ref_ind()] >= min_seed_kb) &&
	  (q[seed_list[seed].end_tar_ind()] - q[seed_list[seed].beg_tar_ind()] >= min_seed_kb));
}

void gapped_alignment::output_alignment(ostream& out_str){
  int gap_counter = 0;
  bool s_score_output;
  
  if(s_scores.size() == tar_restr_al_sites.size())
    s_score_output = true;
  else
    s_score_output = false;

  int i;
  for(i=tar_restr_al_sites.size()-2; i>=0; i--){
    int igap = ref_restr_al_sites[i] - ref_restr_al_sites[i+1];
    int jgap = tar_restr_al_sites[i] - tar_restr_al_sites[i+1];
    gap_counter++;
  }

  out_str<<"gapped alignment:"<<endl;
  out_str<<target_map.read_name<<" -> "<<ref_map.read_name<<endl;
  
  int gapcount = 0;

  for(i=ref_restr_al_sites.size()-1; i>=0; i--){
    if(GapEnd[tar_restr_al_sites[i]][ref_restr_al_sites[i]])
      cout<<"----------------------------------------------------------------------"<<endl;

    if(i < ref_restr_al_sites.size()-1){
      int igap = tar_restr_al_sites[i] - tar_restr_al_sites[i+1];
      int jgap = ref_restr_al_sites[i] - ref_restr_al_sites[i+1];
      	
      int count = 0;
      out_str<<"[ ";
      while(igap>=1){
	count ++;
	igap--;
	out_str<<tar_restr_al_sites[i+1]+count-1<<":";
	out_str<<target_map.map_read[tar_restr_al_sites[i+1]+count-1];
	if(igap>=1) out_str<<", ";
      }
      
      out_str<<" ]->[ ";
      count = 0;
      while(jgap>=1){
	count++;
	jgap--;
	out_str<<ref_restr_al_sites[i+1]+count-1<<":";
	out_str<<ref_map.map_read[ref_restr_al_sites[i+1]+count-1];
	if(jgap >= 1) out_str<<", ";
      }
      out_str<<" ]";
      //out_str<<"  s: "<<S[tar_restr_al_sites[i]][ref_restr_al_sites[i]];
      out_str<<"  t: "<<T[tar_restr_al_sites[i]][ref_restr_al_sites[i]];
      //out_str<<"  m: "<<Mt[tar_restr_al_sites[i]][ref_restr_al_sites[i]];
      
      out_str<<endl;		
    }
    if(GapEnd[tar_restr_al_sites[i]][ref_restr_al_sites[i]]) {
      cout<<"----------------------------------------------------------------------"<<endl;
      gapcount++;
    }
  }

  cout<<"s-score: "<<Smax<<endl;
  cout<<"t-score: "<<Tmax<<endl;
  cout<<"gaps: "<<gapcount<<endl;
  if(best_seed_used)
    out_str<<"includes best scoring local alignment"<<endl;
  else 
    out_str<<"does *NOT* include best scoring local alignment"<<endl;

}
