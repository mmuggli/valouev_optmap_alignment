using namespace std;

#include "alignment.h"

double d_max(double x, double y){
  if(x<y) return y;
  else return x;
}

void rm_alignment::gap_alignment(double gap_open, double gap_ext_kb){
  
  localized_gap_alignment(gap_open, gap_ext_kb, 
			  0, ref_map.map_read.size()-1, 
			  0, target_map.map_read.size()-1);
}

double rm_alignment::t_score_drop(){
  double max_t_score = 0;
  double max_drop = 0;

  for(int i=0; i<t_scores.size(); i++){
    double cur_t_score = t_scores[i];
    if(cur_t_score > max_t_score) max_t_score = cur_t_score;
    double cur_drop = cur_t_score - max_t_score;
    if(cur_drop < max_drop) max_drop = cur_drop;
  }
  return max_drop;
}

double rm_alignment::ref_size(){
  int al_size = ref_restr_al_sites.size();

  if(al_size > 0){
    int al_start = ref_restr_al_sites[al_size-1];
    int al_end = ref_restr_al_sites[0];
    assert(al_start < al_end);

    double res = 0;
    for(int i=al_start; i<al_end; i++){
      res+=ref_map.map_read[i];
    }
    return res;
  }
  else return 0;
}

double rm_alignment::al_ref_size(){
  int al_size = ref_restr_al_sites.size();
  assert(al_size > 1);
  int ref_al_begin = ref_restr_al_sites[al_size-1];
  int ref_al_end = ref_restr_al_sites[0];

  if(ref_al_begin > ref_al_end){
    cerr<<ref_al_begin<<" "<<ref_al_end<<endl;
  }

  assert(ref_al_begin <= ref_al_end);
  assert(ref_al_begin >= 0);
  assert(ref_al_end <= ref_map.map_read.size());

  double res = 0;
  for(int i=ref_al_begin; i<ref_al_end; i++){
    res += ref_map.map_read[i];
  }
  assert(res >= 0);
  return res;
}

void rm_alignment::
localized_gap_alignment(double gap_open, double gap_ext_kb,
			int ref_left_ind, int ref_right_ind, 
			int tar_left_ind, int tar_right_ind){
  Smax = -10000; 
  //set a bad score in case
  //no alignment found
  //equivalent to minus infinity
  
  double score_low_thresh = -5;

  int delta = score_pars.delta;

  int m = target_map.map_read.size();
  int n = ref_map.map_read.size();
  int i, j, g, h, l;

  S.clear();

  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  vector<double> F;
  vector<int> fromF;
  double E;
  int fromE;

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;
 
  F.push_back(0);
  fromF.push_back(-1);

  //fill out the vector of site positions for reference
  for(i=0; i<n; i++){
    counter += ref_map.map_read[i];
    r.push_back(counter);
    F.push_back(0);
    fromF.push_back(-1);
  }
  
  //fill out the vector of site positions for optical map
  q.push_back(0);
  counter=0;
  for(j=0; j<m; j++){
    counter += target_map.map_read[j];
    q.push_back(counter);
  }

  vector<vector< bool > > mminf_matrix; //true if score is mminus inf
  
  for(i=0; i<=m; i++){
    vector<double> z;
    vector<int> mminone;
    vector<bool> bmminf;
    for(j=0; j<=n; j++){      
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    S.push_back(z);
    mminf_matrix.push_back(bmminf);
  }
  

  //initialize S and inf matrix
  for(j=0; j<=n; j++){
    S[0][j] = 0;
    mminf_matrix[0][j] = false; //has been initialized
  }
  
  //end of matrix initialization

  int nextbesti;
  int nextbestj;

  //begin DP computation
  for(i=1; i<=m; i++){
    E = 0;
    fromE = -1;

    for(j=1; j<=n; j++){
      double y; //holds the optimal score
      
      bool y_minf = mminf_matrix[i][j];
      for(g=int_max(0, i-delta); g<i; g++){
	for(h=int_max(0, j-delta); h<j; h++){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    
	    double tar_size = q[i]-q[g];
	    double ref_size = r[j]-r[h];

	    double cur_score = score_pars.ref_total_score_high
	      (tar_size,ref_size,map_gap,ref_gap);
	    //double cur_score = score_pars.total_ref_matching_score1
	    //  (tar_size, ref_size, map_gap, ref_gap);
	      
	    double s = S[g][h] + cur_score;
	    
	    if(s>score_low_thresh){
	      if(y_minf == true){
		y = s;
		fromi[i][j] = g;
		fromj[i][j] = h;	    
		y_minf = false; //y has been set
	      }
	      else{
		if(y < s){
		  fromi[i][j] = g;
		  fromj[i][j] = h;
		  y = s;	    //better choice of y
		}
	      }
	    }
	  }  
	}
      }  

      //Gap (Deletion in Reference)
      if(S[i][j-1] - gap_open > E - (gap_ext_kb * (r[j] - r[j-1])) ) {
	E = S[i][j-1] - gap_open;
	fromE = j-1;
      } else {
	E = E - (gap_ext_kb * (r[j] - r[j-1]));
      }
      
      if(E > y) {
	y = E;
	fromi[i][j] = i;
	fromj[i][j] = fromE;
	y_minf = false;
      }

      //Gap (Deletion in Target)
      if(S[i-1][j] - gap_open > F[j] - (gap_ext_kb * (q[i] - q[i-1])) ) {
	F[j] = S[i-1][j] - gap_open;
	fromF[j] = i-1;
      } else {
	F[j] = F[j] - (gap_ext_kb * (q[i] - q[i-1]));
      }

      if(F[j] > y) {
	y = F[j];
	fromi[i][j] = fromF[j];
	fromj[i][j] = j;
	y_minf = false;
      }

      //Assign best score y from match with possible missing/false or indel gap
      S[i][j] = y;
      mminf_matrix[i][j] = y_minf;

    }
  }
  //end of DP computation

  //search for the optimal alignment
  bool opt_score_minf = true;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(j=0; j<=n; j++){
    if(mminf_matrix[m][j] == false){ //score has been set
      if(opt_score_minf == true){ //opt score not initialized
	Smmax = S[m][j];
	immax = m;
	jmmax = j;
	opt_score_minf = false;
      }
      else{ //opt score already initialized
	if(S[m][j] > Smmax){
	  Smmax = S[m][j];
	  immax = m;
	  jmmax = j;
	}
      }
    }
  }  
  //end of optimal alignment search

  if(immax !=-1 && jmmax != -1){
    traceback(immax, jmmax);
  }
}

void rm_alignment::traceback(int max_i, int max_j) {
  bool end_found = false;
  int curposi = max_i;
  int curposj = max_j;

  tar_restr_al_sites.clear();
  ref_restr_al_sites.clear();

  while(end_found == false){
    tar_restr_al_sites.push_back(curposi);
    ref_restr_al_sites.push_back(curposj);
      
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
  //  fit_t_score();
  Tmax = 0;
}

double rm_alignment::ovlp_p_value(){
  double map1_size = 0;
  double map2_size = 0;

  for(int i=ref_restr_al_sites[ref_restr_al_sites.size()-1];
      i < ref_restr_al_sites[0]; i++){
    double cur_fr = ref_map.map_read[i];
    map1_size += cur_fr;
  }

  for(int i=tar_restr_al_sites[tar_restr_al_sites.size()-1]; 
      i< tar_restr_al_sites[0]; i++){
    double cur_fr = target_map.map_read[i];
    map2_size += cur_fr;
  }

  double ref_size_est = 0.5*(map1_size + map2_size);

  int map1_cuts = 
    ref_restr_al_sites[0]-ref_restr_al_sites[ref_restr_al_sites.size()-1]-1;
  int map2_cuts = 
    tar_restr_al_sites[0]-tar_restr_al_sites[tar_restr_al_sites.size()-1]-1;
  
  int matches = ref_restr_al_sites.size()-2;

  if(matches <= 3){
    cout<<"overlap too short"<<endl;
    return 1;
  }

  int map1_nonmatching = map1_cuts - matches;
  int map2_nonmatching = map2_cuts - matches;

  cout<<char(9)<<"ref_size: "<<ref_size_est;
  cout<<" bad_m1: "<<map1_nonmatching;
  cout<<" ("<<map1_cuts<<")";
  cout<<" bad_m2: "<<map2_nonmatching;
  cout<<" ("<<map2_cuts<<")";
  cout<<endl;

  assert(map1_nonmatching >= 0);
  assert(map2_nonmatching >= 0);
  
  double phi = ((1-score_pars.dig_p)/score_pars.distr_lambda)*ref_size_est +
    score_pars.zeta*ref_size_est;
  
  //double map1_p_value = 
  //poiss_p_value(map1_nonmatching, phi);
  //double map2_p_value = 
  //poiss_p_value(map2_nonmatching, phi);

  double map1_p_value = left_bin_p_value(map1_nonmatching, 0.33,map1_cuts);
  double map2_p_value = left_bin_p_value(map2_nonmatching, 0.33,map2_cuts);
  
  cout<<"m1_p_v: "<<map1_p_value;
  cout<<" m2_p_v: "<<map2_p_value<<endl;

  double match_p_value = map1_p_value*map2_p_value;

  return match_p_value;
				       
}

double rm_alignment::fit_p_value(){
  double min_fr_size = 2.0;
  
  int internal_ref_restr_sites;
  //number of internal (not including the left and right ones)
  //sites of the matching reference region not including
  //sites next to short fragments (<min_fr_size)

  int internal_opt_map_sites;
  //number of internal (not including the left and right ones)
  //sites of the matching optical map region

  int matching_sites; 
  //number of internal matching sites (not including the left and right ones)
  //between optical and reference map
  int miss_cuts = 0;
  int false_cuts = 0;

  int corrected_ref_fr = 0;
  

  //counts the number of matching ref fragments more than min_fr_size
  //not to penalize for short fragments missing in optical maps
  
  double ref_size = 0;
  for(int i=ref_restr_al_sites[ref_restr_al_sites.size()-1];
      i < ref_restr_al_sites[0]; i++){
    double cur_ref_fr = ref_map.map_read[i];
    if(cur_ref_fr > min_fr_size){
      corrected_ref_fr++;
      ref_size += cur_ref_fr;
    }
  }
  internal_ref_restr_sites = corrected_ref_fr - 1;
  
  internal_opt_map_sites =
    tar_restr_al_sites[0] -
    tar_restr_al_sites[tar_restr_al_sites.size()-1] - 1;
  
  matching_sites = ref_restr_al_sites.size()-2;

  miss_cuts = internal_ref_restr_sites - matching_sites;
  false_cuts = internal_opt_map_sites - matching_sites;
  
  cout<<"ref_region: "<<ref_size;
  cout<<" ref_sites: "<<internal_ref_restr_sites;
  cout<<" mc:"<<miss_cuts;
  cout<<" fc:"<<false_cuts;
 
  cout<<endl;

  assert(miss_cuts >= 0);
  assert(false_cuts >= 0);

  double fc_p_value = poiss_p_value(false_cuts,
				    score_pars.zeta*ref_size);

  double mc_p_value = right_bin_p_value
    (miss_cuts, 1-score_pars.dig_p,
     internal_ref_restr_sites);

  cout<<"fc_p_v: "<<fc_p_value<<endl;
  cout<<"mc_p_v: "<<mc_p_value<<endl;

  double p_value=fc_p_value*mc_p_value;

  return p_value;
}

int rm_alignment::ovlp_size(){
  //assert(ref_restr_al_sites.size() >= 2);

  if(ref_restr_al_sites.size() < 2) return 0;
  else{
    int size = ref_restr_al_sites.size();
    return int_max(abs(ref_restr_al_sites[0]-ref_restr_al_sites[size-1]),
		   abs(tar_restr_al_sites[0]-tar_restr_al_sites[size-1]));
  }
}

/*
double rm_alignment::ov_p_value(sim_table& st){
  //this function calculated the p-value corresponding to current overlap
  //before calling this function make sure you compute the fit
  //by calling overlap_alignment
  //also, make sure that you have a precomputed score from simulation
  //that are loaded by calling load_fit_sim_results.
  //the table is computed by additional module in
  //ov_sim.cc

  //assert(!S.empty());

  //cout<<"collection size:"<<st.fit_scores_table[0].size();
  
  assert(!S.empty());
  assert(!st.fit_scores_table.empty());
  assert(!st.fit_lengths.empty());

  int i,j; //for iterative purposes
  double p_value;

  int tar_size = target_map.map_read.size();
  int ref_size = ref_map.map_read.size();

  int cur_clone_size = (tar_size+ref_size)/2; //target_map.map_read.size();
  cur_clone_size = tar_restr_al_sites.size();

  assert(cur_clone_size > 0);
  double score = Smax;

  int index_min = -1; 
  int index_max = -1;
  //indexes of two closest lines from which p-value
  //will be calculated
  
  if(cur_clone_size <= st.fit_lengths[0]){
    index_min = 0;
    index_max = 0;
  }
  if(cur_clone_size >= st.fit_lengths[st.fit_lengths.size()-1]){
    index_min = st.fit_lengths.size()-1;
    index_max = st.fit_lengths.size()-1;
  }
  if(cur_clone_size > st.fit_lengths[0] &&
     cur_clone_size < st.fit_lengths[st.fit_lengths.size()-1]){

    bool found_min = false;
    bool found_max = false;
    //cout<<"size:"<<fit_lengths.size()<<endl;
    for(i=0; i<st.fit_lengths.size(); i++){
      //cout<<i<<endl;
      int cur_len = st.fit_lengths[i];
      if(cur_len <= cur_clone_size){// && found_min == false){
	//found_min = true;
	index_min = i;
      }
      if(cur_len > cur_clone_size && found_max == false){
	found_max = true;
	index_max = i;
      }
    }
  }
  assert(index_min >=0);
  assert(index_max >=0);

  //cout<<"num_cl "<<st.fit_lengths[0]<<endl;
  //cout<<"clone_l:"<<cur_clone_size;
  //cout<<" imin:"<<index_min<<" "<<st.fit_lengths[index_min];
  //cout<<" imax:"<<index_max<<" "<<st.fit_lengths[index_max]<<endl;

  double pmin;
  double pmax;

  int size;
  int counter;

  size = 0;
  counter = 0;
  size = st.fit_scores_table[index_min].size();
  //cout<<"compare to "<<st.fit_scores_table[0].size()<<endl;
  //cout<<size<<endl;
  for(i=0; i<size; i++){
    //cout<<"comparing "<<score<<" to ";
    //cout<<st.fit_scores_table[index_min][i]<<endl;
    if(score < st.fit_scores_table[index_min][i]) counter++;
  }
  //cout<<"counter1:"<<counter;
  
  pmin = (double) counter/size;

  size = 0;
  counter = 0;
  size = st.fit_scores_table[index_max].size();
  for(i=0; i<size; i++){
    if(score < st.fit_scores_table[index_max][i]) counter++;
  }
  //cout<<" counter2:"<<counter<<endl;
  
  pmax = (double) counter/size;

  p_value = (pmin+pmax)/2;
  return p_value;  
}
*/
/*
double rm_alignment::fit_p_value(sim_table& st){
  //this function calculated the p-value corresponding to current fit
  //before calling this function make sure you compute the fit
  //by calling fit_alignment
  //also, make sure that you have a precomputed score from simulation
  //that are loaded by calling load_fit_sim_results.
  //the table is computed by additional module in
  //fit_sim.cc

  //assert(!S.empty());
  assert(!T.empty());
  assert(!st.fit_scores_table.empty());
  assert(!st.fit_lengths.empty());

  int i,j; //for iterative purposes
  double p_value;

  int cur_clone_size = target_map.map_read.size();
  assert(cur_clone_size > 0);
  double score = Tmax;

  int index_min = -1; 
  int index_max = -1;
  //indexes of two closest lines from which p-value
  //will be calculated
  
  if(cur_clone_size <= st.fit_lengths[0]){
    index_min = 0;
    index_max = 0;
  }
  if(cur_clone_size >= st.fit_lengths[st.fit_lengths.size()-1]){
    index_min = st.fit_lengths.size()-1;
    index_max = st.fit_lengths.size()-1;
  }
  if(cur_clone_size > st.fit_lengths[0] &&
     cur_clone_size < st.fit_lengths[st.fit_lengths.size()-1]){

    bool found_min = false;
    bool found_max = false;
    //cout<<"size:"<<fit_lengths.size()<<endl;
    for(i=0; i<st.fit_lengths.size(); i++){
      //cout<<i<<endl;
      int cur_len = st.fit_lengths[i];
      if(cur_len <= cur_clone_size){// && found_min == false){
	//found_min = true;
	index_min = i;
      }
      if(cur_len > cur_clone_size && found_max == false){
	found_max = true;
	index_max = i;
      }
    }
  }
  assert(index_min >=0);
  assert(index_max >=0);

  //cout<<"clone_l:"<<cur_clone_size;
  //cout<<" imin:"<<fit_lengths[index_min];
  //cout<<" imax:"<<fit_lengths[index_max]<<endl;

  double pmin;
  double pmax;

  int size;
  int counter;

  size = 0;
  counter = 0;
  size = st.fit_scores_table[index_min].size();
  for(i=0; i<size; i++){
    if(score < st.fit_scores_table[index_min][i]) counter++;
  }
  
  pmin = (double) counter/size;

  size = 0;
  counter = 0;
  size = st.fit_scores_table[index_max].size();
  for(i=0; i<size; i++){
    if(score < st.fit_scores_table[index_max][i]) counter++;
  }
  
  pmax = (double) counter/size;

  p_value = (pmin+pmax)/2;
  return p_value;
}
*/

void rm_alignment::output_kmers(ostream& out_str){

  int mers = 3; //shortest kmer
  double ref_size_thresh = 2;
  //out_str<<"kmers shared by "<<target_map.read_name;//<<":";
  //out_str<<" and "<<ref_map.read_name<<endl;

  vector<double> ref_tuple;
  vector<double> tar_tuple;

  for(int i=ref_restr_al_sites.size()-1; i>=0; i--){
    if(i < ref_restr_al_sites.size()-1){//print the gap
      int igap = tar_restr_al_sites[i] - tar_restr_al_sites[i+1];
      int jgap = ref_restr_al_sites[i] - ref_restr_al_sites[i+1];
      
      //cmers output
      int _igap = igap;
      int _jgap = jgap;
      //end cmers output

      {



	//cmer output part
       
	if(_igap == 1 && _jgap == 1){ //keep adding
	  double cur_tar = target_map.map_read[tar_restr_al_sites[i+1]];
	  double cur_ref = ref_map.map_read[ref_restr_al_sites[i+1]];
	  
	  //cout<<"adding tar:"<<cur_tar<<" ref:"<<cur_ref<<endl;
	  ref_tuple.push_back(cur_ref);
	  tar_tuple.push_back(cur_tar);
	}
	else{ //flush the data if > # mers, then print 
	  //cout<<"flushing"<<endl;
	  assert(ref_tuple.size() == tar_tuple.size());
	  int tuple_size = ref_tuple.size();
	  if(tuple_size >= mers){
	    for(int m=0; m<tuple_size; m++){
	      out_str<<tar_tuple[m]<<" "<<ref_tuple[m]<<endl;
	    }
	  }
	  //cout<<"flushing tuple of size "<<tuple_size<<endl;
	  ref_tuple.clear();
	  tar_tuple.clear();
	  
	}
	
	//end of cmer output part
      }
      
    }
    else{      
    }
  }
  if(ref_tuple.size()>= mers){
    assert(ref_tuple.size() == tar_tuple.size());
    int tuple_size = ref_tuple.size();
    if(tuple_size >= mers){
      for(int m=0; m<tuple_size; m++){
	double cur_ref_fr = ref_tuple[m];
	if(cur_ref_fr >= ref_size_thresh){
	  out_str<<tar_tuple[m]<<" "<<ref_tuple[m]<<endl;
	}
      }
    }
  }
  ref_tuple.clear();
  tar_tuple.clear();
  
}

void rm_alignment::output_alignment(ostream& out_str){
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


  out_str<<"alignment for "<<target_map.read_name;
  out_str<<" and "<<ref_map.read_name<<endl;

  for(i=ref_restr_al_sites.size()-1; i>=0; i--){
    if(i < ref_restr_al_sites.size()-1){//print the gap
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
      
      out_str<<" ";
      out_str<<"]->[ ";
      count = 0;
      while(jgap>=1){
	count++;
	jgap--;
	out_str<<ref_restr_al_sites[i+1]+count-1<<":";
	out_str<<ref_map.map_read[ref_restr_al_sites[i+1]+count-1];
	if(jgap >= 1) out_str<<", ";
      }
      out_str<<" ]";
      //out_str<<" t: "<<S[tar_restr_al_sites[i]][ref_restr_al_sites[i]];
      //  out_str<<" s: "<<T[tar_restr_al_sites[i]][ref_restr_al_sites[i]];
      if(s_score_output)
	out_str<<" s: "<<s_scores[i];
      
      out_str<<endl;		
    }
  }

  out_str<<"s-score:"<<Smax<<endl;
  out_str<<"t-score:"<<Tmax<<endl;
  out_str<<endl;
}

void rm_alignment::fit_t_score(){
  t_scores.clear();
  t_scores.push_back(0);

  double total_t_score = 0;
  for(int i=ref_restr_al_sites.size()-1; i>0; i--){

    int cur_ref_gap = ref_restr_al_sites[i-1]-ref_restr_al_sites[i];
    int cur_tar_gap = tar_restr_al_sites[i-1]-tar_restr_al_sites[i];

    double cur_ref_size = 0;

    for(int k=0; k<cur_ref_gap; k++){
      cur_ref_size += ref_map.map_read[ref_restr_al_sites[i]+k];
    }
    

    double cur_t_score = score_pars.ref_site_score
      (cur_tar_gap, cur_ref_gap, cur_ref_size);
    //double cur_t_score = score_pars.site_ref_match_score
    //  (cur_tar_gap, cur_ref_gap, cur_ref_size);

    total_t_score += cur_t_score;
    t_scores.push_back(total_t_score);
  }

  Tmax = total_t_score;
}

void rm_alignment::overlap_t_score(){
  //assert(ref_restr_al_sites.size() == tar_restr_al_sites.size());
  t_scores.clear();
  t_scores.push_back(0);
  double t_score = 0;

  for(int i=ref_restr_al_sites.size()-1; i>0; i--){

    int cur_ref_gap = ref_restr_al_sites[i-1]-ref_restr_al_sites[i];
    int cur_tar_gap = tar_restr_al_sites[i-1]-tar_restr_al_sites[i];

    double cur_t_score = 
      score_pars.opt_site_score(cur_tar_gap, cur_ref_gap);

    //double cur_t_score = 
    //  score_pars.site_opt_match_score(cur_tar_gap, cur_ref_gap);
    t_score += cur_t_score;
    t_scores.push_back(t_score);
  }

  Tmax = t_score;
}

void rm_alignment::fit_alignment(){
  Smax = -10000; 
  //set a bad score in case
  //no alignment found
  //equivalent to minus infinity

  double score_low_thresh = -5;

  int delta = score_pars.delta;

  int end_delta = 2;

  int m = target_map.map_read.size();
  int n = ref_map.map_read.size();
  int i, j, g, h, l;

  S.clear();

  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;
 
  //fill out the vector of site positions for reference
  for(i=0; i<n; i++){
    counter += ref_map.map_read[i];
    r.push_back(counter);
  }
  
  //fill out the vector of site positions for optical map
  q.push_back(0);
  counter=0;
  for(j=0; j<m; j++){
    counter += target_map.map_read[j];
    q.push_back(counter);
  }

  vector<vector< bool > > mminf_matrix; //true if score is mminus inf
  
  for(i=0; i<=m; i++){
    vector<double> z;
    vector<int> mminone;
    vector<bool> bmminf;
    for(j=0; j<=n; j++){
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    S.push_back(z);
    mminf_matrix.push_back(bmminf);
  }
  

  //initialize T and inf matrix
  for(j=0; j<=n; j++){
    S[0][j] = 0;
    mminf_matrix[0][j] = false; //has been initialized
  }
  for(j=0; j<=m; j++){
    S[j][0] = 0;
    mminf_matrix[j][0] = false;
  }
  
  //end of matrix initialization

  //begin DP computation
  for(i=1; i<=m; i++){
    for(j=1; j<=n; j++){
      double y; //holds the optimal score
      
      bool y_minf = mminf_matrix[i][j];
      for(g=int_max(0, i-delta); g<i; g++){
	for(h=int_max(0, j-delta); h<j; h++){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    
	    double tar_size = q[i]-q[g];
	    double ref_size = r[j]-r[h];

	    double cur_score = 	
	      score_pars.ref_total_score_high
	      (tar_size, ref_size, map_gap, ref_gap);
	      //score_pars.total_ref_matching_score_mult
	      //(tar_size, ref_size, map_gap, ref_gap);
	      //score_pars.total_ref_matching_score_high1
	      //(tar_size,ref_size, map_gap, ref_gap);
	      //score_pars.total_ref_matching_score1
	      //(tar_size,ref_size, map_gap, ref_gap);
	      
	    double s = S[g][h] + cur_score;
	    
	    if(s>score_low_thresh){
	    // {
	      if(y_minf == true){
		y = s;
		fromi[i][j] = g;
		fromj[i][j] = h;	    
		y_minf = false; //y has been set
	      }
	      else{
		if(y < s){
		  fromi[i][j] = g;
		  fromj[i][j] = h;
		  y = s;	    //better choice of y
		}
	      }
	    }
	  }	  
	}
      }      

      S[i][j] = y;
      mminf_matrix[i][j] = y_minf;
    }
  }
  //end of DP computation
  
  //search for the optimal alignment
  bool opt_score_minf = true;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(j=0; j<=n; j++){
    for(i=int_max(0,m-end_delta); i<=m; i++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][j] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }  

  for(i=0; i<=m; i++){
    for(j=int_max(0,n-end_delta); j<=n; j++){
      if(mminf_matrix[i][j] == false){
	if(opt_score_minf == true){
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{
	  if(S[i][j] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  //end of optimal alignment search
  //cout<<"ref_s: "<<n<<" map_s: "<<m;
  //cout<<" end found: "<<jmmax<<" "<<immax<<endl;

  if(immax !=-1 && jmmax != -1){
    Smax = Smmax; //set the optimal score

    //trace back procedure to recover the full alignment
    bool end_found = false;
    
    int curposi = immax;
    int curposj = jmmax;
    
    while(end_found == false){
      tar_restr_al_sites.push_back(curposi);
      ref_restr_al_sites.push_back(curposj);
      
      int newi;
      int newj;
      
      newi = fromi[curposi][curposj];
      newj = fromj[curposi][curposj];
      
      curposi = newi;
      curposj = newj;

      //cout<<"curposi:"<<curposi<<" curposj:"<<curposj<<endl;
      
      if(curposi == -1 || curposj == -1){
	end_found = true;
      }
    }
  }
}




void rm_alignment::
localized_fit_alignment(int ref_left_ind, int ref_right_ind,
			int tar_left_ind, int tar_right_ind){

  assert(ref_left_ind >= 0);
  assert(ref_right_ind <= ref_map.map_read.size());
  assert(ref_left_ind < ref_right_ind);

  assert(tar_left_ind >= 0);
  assert(tar_right_ind <= target_map.map_read.size());
  assert(tar_left_ind < tar_right_ind);

  Smax = -10000; 
  //set a bad score in case
  //no alignment found
  //equivalent to minus infinity

  double score_low_thresh = -1500;
  score_pars.fit_score_mult_thresh = 
    sqr(score_pars.distr_sigma)*15;

  int delta = score_pars.delta;

  int m = tar_right_ind - tar_left_ind;
  int n = ref_right_ind - ref_left_ind;
  int i, j, g, h, l;

  S.clear();

  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
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

  vector<vector< bool > > mminf_matrix; //true if score is mminus inf
  
  for(i=0; i<=m; i++){
    vector<double> z;
    vector<int> mminone;
    vector<bool> bmminf;
    for(j=0; j<=n; j++){
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    S.push_back(z);
    mminf_matrix.push_back(bmminf);
  }
  

  //initialize T and inf matrix
  for(j=0; j<=n; j++){
    S[0][j] = 0;
    mminf_matrix[0][j] = false; //has been initialized
  }
  
  //end of matrix initialization

  //begin DP computation
  for(i=1; i<=m; i++){
    for(j=1; j<=n; j++){
      double y; //holds the optimal score
      
      bool y_minf = mminf_matrix[i][j];
      for(g=int_max(0, i-delta); g<i; g++){
	for(h=int_max(0, j-delta); h<j; h++){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    
	    double tar_size = q[i]-q[g];
	    double ref_size = r[j]-r[h];

	    double cur_score = 	
	      score_pars.ref_total_score_high
	      (tar_size, ref_size, map_gap, ref_gap);
	      //score_pars.total_ref_matching_score_mult
	      //(tar_size, ref_size, map_gap, ref_gap);
	      //score_pars.total_ref_matching_score_high1
	      //(tar_size,ref_size, map_gap, ref_gap);
	      //score_pars.total_ref_matching_score1
	      //(tar_size,ref_size, map_gap, ref_gap);
	      
	    double s = S[g][h] + cur_score;
	    
	    if(s>score_low_thresh){
	    // {
	      if(y_minf == true){
		y = s;
		fromi[i][j] = g;
		fromj[i][j] = h;	    
		y_minf = false; //y has been set
	      }
	      else{
		if(y < s){
		  fromi[i][j] = g;
		  fromj[i][j] = h;
		  y = s;	    //better choice of y
		}
	      }
	    }
	  }	  
	}
      }      

      S[i][j] = y;
      mminf_matrix[i][j] = y_minf;
    }
  }
  //end of DP computation
  

  //search for the optimal alignment
  bool opt_score_minf = true;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(j=0; j<=n; j++){
    if(mminf_matrix[m][j] == false){ //score has been set
      if(opt_score_minf == true){ //opt score not initialized
	Smmax = S[m][j];
	immax = m;
	jmmax = j;
	opt_score_minf = false;
      }
      else{ //opt score already initialized
	if(S[m][j] > Smmax){
	  Smmax = S[m][j];
	  immax = m;
	  jmmax = j;
	}
      }
    }
  }  
  //end of optimal alignment search

  if(immax !=-1 && jmmax != -1){
    Smax = Smmax; //set the optimal score

    //trace back procedure to recover the full alignment

    s_scores.clear();

    bool end_found = false;
    
    int curposi = immax;
    int curposj = jmmax;
    
    while(end_found == false){

      tar_restr_al_sites.push_back(curposi+tar_left_ind);
      ref_restr_al_sites.push_back(curposj+ref_left_ind);

      s_scores.push_back(S[curposi][curposj]);
      
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
  }
}

void rm_alignment::
fast_localized_fit_alignment(int ref_left_ind, int ref_right_ind,
			     int tar_left_ind, int tar_right_ind){
  Smax = -10000; //set a bad score in case
  //no alignment found

  assert(ref_left_ind >= 0);
  assert(ref_right_ind <= ref_map.map_read.size());
  assert(ref_left_ind < ref_right_ind);

  assert(tar_left_ind >= 0);
  assert(tar_right_ind <= target_map.map_read.size());
  assert(tar_left_ind < tar_right_ind);

  double score_low_thresh = -5;
  double t_score_low_thresh = -5;
  double t_score_mult = 0.5;
  double t_offset = -5;

  double minf = -10000;

  int delta1 = 3;
  int delta2 = 5;
  
  int m = tar_right_ind - tar_left_ind;
  int n = ref_right_ind - ref_left_ind;

  int i, j, g, h, l;

  vector<int> matrix_mask;
  vector<int> indexj;
  vector<int> indexj_helper;
  vector<bool> indexj_mask;

  S.clear();
  T.clear();

  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
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

  vector<vector< bool > > mminf_matrix; //true if score is mminus inf
  
  for(i=0; i<=m; i++){
    vector<double> z;
    vector<int> mminone;
    vector<bool> bmminf;
    for(j=0; j<=n; j++){
      //z.push_back(0);
      z.push_back(minf);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    S.push_back(z);
    T.push_back(z);
    mminf_matrix.push_back(bmminf);
  }
  

  //initialize S, T and inf matrix
  for(j=0; j<=n; j++){
    matrix_mask.push_back(delta1+1);
    indexj.push_back(j);
    indexj_mask.push_back(true);
    S[0][j] = 0;
    T[0][j] = 0;
    mminf_matrix[0][j] = false; //has been initialized
  }
  
  //end of matrix initialization

  //begin DP computation
  for(i=1; i<=m; i++){
   
    indexj_helper.clear();

    for(int a=0; a<indexj.size(); a++){
      j = indexj[a];
      matrix_mask[j]--;
      if(matrix_mask[j]>=0){ 
	indexj_helper.push_back(j);
      }      
      else{
	indexj_mask[j] = false; //kicked_out
      }

      double y; //holds the optimal score
      
      bool y_minf = mminf_matrix[i][j];
      bool score_to_set = false;
     
      for(h=j-1; h>=int_max(0,j-delta2); h--){
	for(g=i-1; g>=int_max(0,i-delta1); g--){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    double tar_size = q[i]-q[g];
	    double ref_size = r[j]-r[h];

	    double cur_score = 	
	      score_pars.ref_total_score_high
	      (tar_size, ref_size, map_gap, ref_gap);
	      //score_pars.total_ref_matching_score_mult
	      //(tar_size, ref_size, map_gap, ref_gap);
	      //score_pars.total_ref_matching_score_high1
	      //(tar_size,ref_size, map_gap, ref_gap);
	      //score_pars.total_ref_matching_score1
	      //(tar_size,ref_size, map_gap, ref_gap);
	    
	      
	    double s = S[g][h] + cur_score;
	    
	    if(s>score_low_thresh){
	      score_to_set = true;

	      if(y_minf == true){
		y = s;
		fromi[i][j] = g;
		fromj[i][j] = h;	    
		y_minf = false; //y has been set
	      }
	      else{
		if(y < s){
		  fromi[i][j] = g;
		  fromj[i][j] = h;
		  y = s;	    //better choice of y
		}
	      }
	    }
	  }	  
	}
      }
      
      if(score_to_set){
	int ilink = fromi[i][j];
	int jlink = fromj[i][j];
       
	double ref_size = r[j] - r[jlink];

	double loc_t_score = score_pars.ref_site_score
	  (i-ilink, j-jlink,ref_size);
	//double loc_t_score = score_pars.site_ref_match_score
	//  (i-ilink,j-jlink,ref_size);
	double cur_t_score = T[ilink][jlink]+loc_t_score;
       
	if(cur_t_score > t_offset + t_score_mult*i){
     
	  for(int b=j+1; b<=int_min(j+delta2,n); b++){
	 
	    if(indexj_mask[b] == false){ //cannot be pushed above
	      matrix_mask[b] = delta1;
	      indexj_helper.push_back(b);
	      indexj_mask[b] = true;
	    }
	    else{
	      matrix_mask[b] = delta1+1;
	    }
	  }
	  T[i][j] = cur_t_score;
	  S[i][j] = y;
	  mminf_matrix[i][j] = y_minf;
	}
      }
    }
    indexj = indexj_helper;
  }
  //end of DP computation
  

  //search for the optimal alignment
  bool opt_score_minf = true;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(j=0; j<=n; j++){
    if(mminf_matrix[m][j] == false){ //score has been set
      if(opt_score_minf == true){ //opt score not initialized
	Smmax = S[m][j];
	immax = m;
	jmmax = j;
	opt_score_minf = false;
      }
      else{ //opt score already initialized
	if(S[m][j] > Smmax){
	  Smmax = S[m][j];
	  immax = m;
	  jmmax = j;
	}
      }
    }
  }  
  //end of optimal alignment search

  if(immax !=-1 && jmmax != -1){
    Smax = Smmax; //set the optimal score

    //trace back procedure to recover the full alignment
    bool end_found = false;
    
    int curposi = immax;
    int curposj = jmmax;
    
    while(end_found == false){
      tar_restr_al_sites.push_back(curposi+tar_left_ind);
      ref_restr_al_sites.push_back(curposj+ref_left_ind);
      
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
  }
}

void rm_alignment::fast_fit_alignment(){
  double lambda = score_pars.lambda;
  int delta = score_pars.delta;

  double std_mult = 3.0;

  int m = target_map.map_read.size();
  int n = ref_map.map_read.size();
  int i, j, g, h, l;

  Smax = 0;

  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;
  
  for(i=0; i<n; i++){
    counter += ref_map.map_read[i];
    r.push_back(counter);
  }
  q.push_back(0);
  counter=0;
  for(j=0; j<m; j++){
    counter += target_map.map_read[j];
    q.push_back(counter);
  }

  bool first_al = true;
  double optimal_score;
  bool total_score_set = false;
  for(i=0; i<=n; i++){
  
    bool alignment_acceptable = true;

    vector<int> tar_al_inds;
    vector<int> ref_al_inds;

    bool alignment_end = false;
    
    int last_tar_al_site = 0;
    int last_ref_al_site = i;
    
    tar_al_inds.push_back(last_tar_al_site);
    ref_al_inds.push_back(last_ref_al_site);

    double last_al_score = 0;

    while(!alignment_end){
      double cur_ref_size = 0;
      double cur_max_score;
      bool score_set = false;

      int opt_tar_ind;
      int opt_ref_ind;

      for(g=last_ref_al_site+1; g<=int_min(n, last_ref_al_site+delta); g++){
	double cur_ref_fr = r[g]-r[g-1];
	cur_ref_size += cur_ref_fr;

	bool no_further_match = false;
	int for_tar_ind = 1;
	double cur_tar_size = 0;

	while(!no_further_match){
	  double cur_tar_fr = q[last_tar_al_site + for_tar_ind] - 
	    q[last_tar_al_site+for_tar_ind-1];
	  cur_tar_size += cur_tar_fr;
	 
	  //cout<<"["<<last_tar_al_site+for_tar_ind<<": "<<cur_tar_fr<<"] -> [";
	  //cout<<g<<": "<<cur_ref_fr<<"]"<<endl;
	  
	  //cout<<"total: tar_fr: "<<for_tar_ind;
	  //cout<<" tar: "<<cur_tar_size;
	  //cout<<" -> ref_fr: "<<g-last_ref_al_site;
	  //cout<<" ref: "<<cur_ref_size;
	 
	

	  if(cur_tar_size < cur_ref_size 
	     + sqrt(cur_ref_size)*score_pars.distr_sigma*std_mult){
	    double cur_score = score_pars.ref_total_score_high
	      (cur_tar_size, cur_ref_size, for_tar_ind, g-last_ref_al_site);
	    //double cur_score = score_pars.total_ref_matching_score1
	    //  (cur_tar_size, cur_ref_size ,for_tar_ind, g-last_ref_al_site);
	    //cout<<" score:"<<cur_score;
	    bool update_score = false;
	    if(!score_set) update_score = true;
	    if(cur_score > cur_max_score) update_score = true;
	    if(update_score){
	      score_set = true;
	      //cout<<endl<<"setting score to "<<cur_score<<endl;
	      cur_max_score = cur_score;
	      opt_tar_ind = last_tar_al_site + for_tar_ind;
	      opt_ref_ind = g;
	    }
	    for_tar_ind++;
	    if(for_tar_ind + last_tar_al_site > m){
	      no_further_match = true;
	      //alignment_end = true;
	      //if(score_set) total_score_set = true;
	    }
	  }
	  else{
	    no_further_match = true;
	  }
	  //cout<<endl<<endl;
	  
	}
      }

      if(!score_set){ //couldn't find reasonable alignment
	alignment_end = true;
      }
      else{
	if(opt_tar_ind == m)
	  alignment_end = true;
      }


      last_tar_al_site = opt_tar_ind;
      last_ref_al_site = opt_ref_ind;
      last_al_score += cur_max_score;

      //cout<<"max_score:"<<cur_max_score;
      //cout<<" total: "<<last_al_score<<endl;
      
      if(score_set){
	//cout<<"index store: ("<<last_tar_al_site;
	//cout<<", "<<last_ref_al_site<<")"<<endl;
	tar_al_inds.push_back(last_tar_al_site);
	ref_al_inds.push_back(last_ref_al_site);
      }
      //cout<<"******************"<<endl;
    }
    

    //close the gap between the last site and
    //the end of the alignment to have a "fit"
    if(last_tar_al_site <m && last_ref_al_site >= n)
      alignment_acceptable = false;

    if(last_tar_al_site < m && alignment_acceptable){
      int map_gap = m - last_tar_al_site;
      int gap_multiplier = 10;
      assert(map_gap >= 0);
      
      double tar_gap_size = 0;
      for(int k=last_tar_al_site+1; k<=m; k++){
	double cur_tar_fr = q[k]-q[k-1];
	tar_gap_size += cur_tar_fr;
      }
      
      double best_cur_score;
      int best_ref_ind;
      double ref_gap_size = 0;
      for(int k= last_ref_al_site+1; 
	  k<=int_min(n,last_ref_al_site + map_gap*gap_multiplier); k++){
	double cur_ref_fr = r[k]-r[k-1];
	int cur_ref_gap = k-last_ref_al_site;
	ref_gap_size += cur_ref_fr;
	assert(cur_ref_fr>0);
	
	double cur_score = score_pars.ref_total_score_high
	  (ref_gap_size, tar_gap_size, map_gap, cur_ref_gap);
	//double cur_score = score_pars.total_ref_matching_score_high1
	//  (ref_gap_size, tar_gap_size, map_gap, cur_ref_gap);
	if(k==last_ref_al_site+1){
	  best_cur_score = cur_score;
	  best_ref_ind = k;
	}
	
	if(cur_score > best_cur_score){
	  best_cur_score = cur_score;
	  best_ref_ind = k;
	}
      }
      //cout<<"found best match to site: "<<best_ref_ind<<endl;
      tar_al_inds.push_back(m);
      ref_al_inds.push_back(best_ref_ind);

      last_al_score += best_cur_score;
    }
    //end of gap close

    if(total_score_set = true){
      bool update_opt_al = false;
      if(first_al){
	first_al = false;
	optimal_score = last_al_score;
	update_opt_al = true;
      }
      else{
	if(last_al_score > optimal_score){
	  optimal_score = last_al_score;
	  update_opt_al = true;
	}
      }

      if(update_opt_al && alignment_acceptable){
	ref_restr_al_sites.clear();
	tar_restr_al_sites.clear();
	for(int k=tar_al_inds.size()-1; k>=0; k--){
	  //for(int k = 0; k<tar_al_inds.size(); k++){
	  int cur_tar_site = tar_al_inds[k];
	  int cur_ref_site = ref_al_inds[k];
	  ref_restr_al_sites.push_back(cur_ref_site);
	  tar_restr_al_sites.push_back(cur_tar_site);
	}
	Smax = optimal_score;
      }
    }

    //cout<<i<<": total score: "<<last_al_score<<endl;
  }
  //cout<<"opt score set: "<<total_score_set;
  //cout<<" opt_score: "<<optimal_score<<endl;
}


/*
void rm_alignment::open_end_fit_alignment(){
  double lambda = score_pars.lambda;
  int delta = score_pars.delta;

  int m = target_map.map_read.size();
  int n = ref_map.map_read.size();
  int i, j, g, h, l;
 
  vector<vector< double > > X;

  T.clear();
  
  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;
  
  for(i=0; i<n; i++){
    counter += ref_map.map_read[i];
    r.push_back(counter);
  }
  q.push_back(0);
  counter=0;
  for(j=0; j<m; j++){
    counter += target_map.map_read[j];
    q.push_back(counter);
  }

  vector<vector< bool > > mminf_matrix; //true if score is mminus inf
  
  for(i=0; i<=m; i++){
    vector<double> z;
    vector<int> mminone;
    vector<bool> bmminf;
    for(j=0; j<=n; j++){
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    T.push_back(z);
    //S.push_back(z);
    mminf_matrix.push_back(bmminf);
  }
  

  //initialize T and inf matrix
  for(j=0; j<=n; j++){
    T[0][j] = 0;
    mminf_matrix[0][j] = false; //has been initialized
  }
  for(i=0; i<=m; i++){
    T[i][0] = -i*lambda;
    mminf_matrix[i][0] = false; //has been initialized
  }
  
  //end of matrix initialization

  //begin DP computation
  for(i=1; i<=m; i++){
    for(j=1; j<=n; j++){
      double y; //holds the optimal score

      bool y_minf = mminf_matrix[i][j];
      for(g=int_mmax(0, i-delta); g<i; g++){
	for(h=int_mmax(0, j-delta); h<j; h++){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    assert(ref_gap > 0);
	    assert(map_gap > 0);
	    double tar_size = q[i]-q[g];
	    double ref_size = r[j]-r[h];

	    assert(ref_size > 0);
	    assert(tar_size > 0);
	    	    
	    double cur_score;
	    
	    if(g==0 || i==m){
	      cur_score = score_pars.total_end_matching_ref_score
		(tar_size, ref_size, map_gap, ref_gap);
	    }
	    else{
	      cur_score = score_pars.total_ref_matching_score_mult
		(tar_size, ref_size, map_gap, ref_gap);
	    }
	    
	    double s = T[g][h] + cur_score;

	    double cur_t_match_score = 
	      score_pars.site_ref_match_score(map_gap,ref_gap,ref_size);

	    
	    if(y_minf == true){
	      y = s;
	      //z = site_match_score;
	      fromi[i][j] = g;
	      fromj[i][j] = h;	    
	      y_minf = false; //y has been set
	    }
	    else{
	      if(y < s){
		fromi[i][j] = g;
		fromj[i][j] = h;
		y = s;	    //better choice of y
		//z = site_match_score;
	      }
	    }
	  }	  
	}
      }

      if(i>1 || j>1 && delta >0 ) assert(y_minf == false);

      S[i][j] = y;
      
      mminf_matrix[i][j] = y_minf;
      //message end

    }
  }
  for(i=0; i<=m; i++){
    S[i][n] -= score_pars.lambda*(m-i);
  }
  //end of DP computation
  


  //search for the optimal alignment
  bool opt_score_minf = true;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(j=0; j<=n; j++){
    if(mminf_matrix[m][j] == false){ //score has been set
      if(opt_score_minf == true){ //opt score not initialized
	Smmax = S[m][j];
	immax = m;
	jmmax = j;
	opt_score_minf = false;
      }
      else{ //opt score already initialized
	if(S[m][j] > Smmax){
	  Smmax = S[m][j];
	  immax = m;
	  jmmax = j;
	}
      }
    }
  }
  for(i=0; i<=m; i++){
    if(mminf_matrix[i][n] == false){ //score has been set
      if(opt_score_minf == true){ //opt score not initialized
	Smmax = S[i][n];
	immax = i;
	jmmax = n;
	opt_score_minf = false;
      }
      else{ //opt score already initialized
	if(S[i][n] > Smmax){
	  Smmax = S[i][n];
	  immax = i;
	  jmmax = n;
	}
      }
    }
  }
  assert(opt_score_minf == false);
  //end of optimal alignment search

  Smax = Smmax; 

  //trace back procedure to recover the full alignment
  bool end_found = false;
  
  int curposi = immax;
  int curposj = jmmax;

  while(end_found == false){
    if(curposi == m || curposi == 0){
    }
    else{
      tar_restr_al_sites.push_back(curposi);
      ref_restr_al_sites.push_back(curposj);
    }

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
}
*/
/*
void rm_alignment::declumping_fit_alignment(int max_als, bool new_scoring){
  //int max_als = 10;
  assert(max_als > 0);

  double lambda = score_pars.lambda;
  int delta = score_pars.delta;

  int m = target_map.map_read.size();
  int n = ref_map.map_read.size();
  int i, j, g, h, l;
 
  vector<vector< bool > > M;
  //masking matrix for declumping

  declumped_best_scores.clear();

  T.clear();
  //S.clear();
  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;
  
  for(i=0; i<n; i++){
    counter += ref_map.map_read[i];
    r.push_back(counter);
  }
  q.push_back(0);
  counter=0;
  for(j=0; j<m; j++){
    counter += target_map.map_read[j];
    q.push_back(counter);
  }

  vector<vector< bool > > mminf_matrix; //true if score is mminus inf
  
  for(i=0; i<=m; i++){
    vector<double> z;
    vector<int> mminone;
    vector<bool> bmminf;
  
    for(j=0; j<=n; j++){
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    T.push_back(z);
    //S.push_back(z);
    M.push_back(bmminf);
    mminf_matrix.push_back(bmminf);
  }
  
  for(int u=0; u<max_als; u++){ //declumping cycle
    //initialize T and inf matrix

    for(i=0; i<=m; i++){
      for(j=0; j<=n; j++){
	mminf_matrix[i][j] = true;
      }
    }
    for(j=0; j<=n; j++){
      T[0][j] = 0;
      mminf_matrix[0][j] = false; //has been initialized
    }


    //end of matrix initialization
    
    //begin DP computation
    for(i=1; i<=m; i++){
      for(j=1; j<=n; j++){
	double y; //holds the optimal score
	double z; //holds additional score
       
	if(M[i][j] == false){//marked for no use, set minf
	  assert(mminf_matrix[i][j] == true);
	}
	else{
	  bool y_minf = mminf_matrix[i][j];
	  for(g=int_mmax(0, i-delta); g<i; g++){
	    for(h=int_mmax(0, j-delta); h<j; h++){
	      
	      if(mminf_matrix[g][h] == false 
		 && M[g][h] == true){ 
		//score has been set before and not masked
		int map_gap  = i-g;
		int ref_gap = j-h;
		assert(ref_gap > 0);
		assert(map_gap > 0);
		double tar_size = q[i]-q[g];
		double ref_size = r[j]-r[h];
		
		assert(ref_size > 0);
		assert(tar_size > 0);
		
		double cur_score;
		if(new_scoring == true){
		  cur_score = score_pars.total_ref_matching_score_mult
		  (tar_size, ref_size, map_gap, ref_gap);
		}
		else{
		  cur_score = score_pars.old_total_ref_matching_score_mult
		  (tar_size, ref_size, map_gap, ref_gap);		  
		}
		
		double s = T[g][h] + cur_score;
		
		double cur_t_match_score = 0;
		//score_pars.site_ref_match_score(map_gap,ref_gap,ref_size);
		//double site_match_score = S[g][h]+cur_t_match_score;
		
		if(y_minf == true){
		  y = s;
		  //z = site_match_score;
		  fromi[i][j] = g;
		  fromj[i][j] = h;	    
		  y_minf = false; //y has been set
		}
		else{
		  if(y < s){
		    fromi[i][j] = g;
		    fromj[i][j] = h;
		    y = s;	    //better choice of y
		    //z = site_match_score;
		  }
		}
	      }	  
	    }
	  }
		
	  T[i][j] = y;
	 
	  mminf_matrix[i][j] = y_minf;
	  //message end
	}	
      }
    }
    for(i=0; i<=m; i++){
      T[i][n] -= score_pars.lambda*(m-i);
    }
    //end of DP computation
    
    
    //search for the optimal alignment
    bool opt_score_minf = true;
    double Tmmax;
    int immax = -1;
    int jmmax = -1;
    
    for(j=0; j<=n; j++){
      if(mminf_matrix[m][j] == false){ //score has been set
	assert(M[m][j] == true);
	if(opt_score_minf == true){ //opt score not initialized
	  Tmmax = T[m][j];
	  immax = m;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(T[m][j] > Tmmax){
	    Tmmax = T[m][j];
	    immax = m;
	    jmmax = j;
	  }
	}
      }
    }

    if(opt_score_minf == true){
      cout<<"no more alignments"<<endl;
    }
    assert(opt_score_minf == false);
    //end of optimal alignment search
    
    cout<<"alignment pass:"<<u<<" Tmax:"<<Tmmax;
    cout<<" eding @ "<<immax<<","<<jmmax<<endl;
    Tmax = Tmmax; //to be removed?
    //Smax = S[immax][jmmax];

    declumped_best_scores.push_back(Tmax);
    
    //cout<<"The highest scoring alignment for ";
    //cout<<target_map.read_name<<":"<<endl;
    //cout<<Smmax<<" at ("<<immax<<", "<<jmmax<<")"<<endl;
    
    //trace back procedure to recover the full alignment
    bool end_found = false;
    
    int curposi = immax;
    int curposj = jmmax;

    tar_restr_al_sites.clear();
    ref_restr_al_sites.clear();

    while(end_found == false){
      tar_restr_al_sites.push_back(curposi);
      ref_restr_al_sites.push_back(curposj);
      
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
    assert(tar_restr_al_sites.size() >= 2);
    for(int k=tar_restr_al_sites.size()-1; k>=0; k--){
      int curi = tar_restr_al_sites[k];
      int curj = ref_restr_al_sites[k];
      M[curi][curj] = false;
      //cout<<"masking ("<<curi<<","<<curj<<")"<<endl;
    }
    if(u==0){
      tmp1 = ref_restr_al_sites[ref_restr_al_sites.size()-1];
      tmp2 = ref_restr_al_sites[0];
    }
  }
}
*/


void rm_alignment::overlap_alignment(){
  //double lambda = score_pars.lambda;
  int delta = score_pars.delta;

  int m = target_map.map_read.size();
  int n = ref_map.map_read.size();
  int i, j, g, h, l;
 
  double local_score_thresh = -100;

  //T.clear();
  S.clear();
  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  s_scores.clear();

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;
  
  for(i=0; i<n; i++){
    counter += ref_map.map_read[i];
    r.push_back(counter);
  }
  q.push_back(0);
  counter=0;
  for(j=0; j<m; j++){
    counter += target_map.map_read[j];
    q.push_back(counter);
  }

  vector<vector< bool > > mminf_matrix; //true if score is mminus inf
  
  for(i=0; i<=m; i++){
    vector<double> z;
    vector<int> mminone;
    vector<bool> bmminf;
    for(j=0; j<=n; j++){
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    // T.push_back(z);
    S.push_back(z);
    mminf_matrix.push_back(bmminf);
  }
  

  //initialize S and inf matrix
  for(j=0; j<=n; j++){
    S[0][j] = 0;
    mminf_matrix[0][j] = false; //has been initialized
  }
  for(i=0; i<=m; i++){
    S[i][0] = 0;
    mminf_matrix[i][0] = false; //has been initialized
  }
  
  //end of matrix initialization

  //begin DP computation
  for(i=1; i<=m; i++){
    for(j=1; j<=n; j++){
      double y;
      //double t_score;

      bool y_minf = mminf_matrix[i][j];
      for(g=int_max(0, i-delta); g<i; g++){
	for(h=int_max(0, j-delta); h<j; h++){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    assert(ref_gap > 0);
	    assert(map_gap > 0);	  

	    double s = S[g][h] +
	      score_pars.opt_size_score
	      (q[i]-q[g],r[j]-r[h],map_gap,ref_gap);
	    //double s = S[g][h] + 
	    //  score_pars.total_opt_matching_score_mult
	    //  (q[i]-q[g], r[j]-r[h], map_gap, ref_gap);

	    //double cur_score = score_pars.total_opt_matching_score_mult
	    //(q[i]-q[g], r[j]-r[h], map_gap, ref_gap);
	    //cout<<cur_score<<endl;
	    
	    
	    //double cur_t = S[g][h] + 
	    //score_pars.site_opt_match_score(map_gap, ref_gap);
	    //if(s > local_score_thresh){
	    {  
	      if(y_minf == true){
		y = s;
		//t_score = cur_t;
		fromi[i][j] = g;
		fromj[i][j] = h;	    
		y_minf = false; //y has been set
	      }
	      else{
		if(y < s){
		  fromi[i][j] = g;
		  fromj[i][j] = h;
		  y = s;	    //better choice of y
		  //t_score = cur_t;
		}
	      }
	    }
	  }	  
	}
      }

      if(y_minf == false){
	S[i][j] = y;
	//T[i][j] = score_pars.site_opt_match_score(i-fromi[i][j], j-fromj[i][j]);
	mminf_matrix[i][j] = y_minf;	
      }
     

    }
  }
  //end of DP computation
  

  //search for the optimal alignment
  bool opt_score_minf = true;
  //double Tmmax;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(j=0; j<=n; j++){
    if(mminf_matrix[m][j] == false){ //score has been set
      if(opt_score_minf == true){ //opt score not initialized
	Smmax = S[m][j];
	immax = m;
	jmmax = j;
	opt_score_minf = false;
      }
      else{ //opt score already initialized
	if(S[m][j] > Smmax){
	  Smmax = S[m][j];
	  immax = m;
	  jmmax = j;
	}
      }
    }
  }
  for(int i=0; i<=m; i++){
    if(mminf_matrix[i][n] == false){ //score has been set
      if(opt_score_minf == true){ //opt score not initialized
	Smmax = S[i][n];
	immax = i;
	jmmax = n;
	opt_score_minf = false;
      }
      else{ //opt score already initialized
	if(S[i][n] > Smmax){
	  Smmax = S[i][n];
	  immax = i;
	  jmmax = n;
	}
      }
    }
  }
  assert(opt_score_minf == false);
  //end of optimal alignment search

  Smax = Smmax; 

  //trace back procedure to recover the full alignment
  bool end_found = false;
  
  int curposi = immax;
  int curposj = jmmax;

  while(end_found == false){
    tar_restr_al_sites.push_back(curposi);
    ref_restr_al_sites.push_back(curposj);

    s_scores.push_back(S[curposi][curposj]);

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
}

void rm_alignment::optimized_overlap_alignment(){

  int end_delta = 0;
  //double lambda = score_pars.lambda;
  int delta = score_pars.delta;

  int m = target_map.map_read.size();
  int n = ref_map.map_read.size();
  int i, j, g, h, l;
 
  double local_score_thresh = -10;
  double local_t_score_thresh = -5;//-5;

  T.clear();
  S.clear();
  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  s_scores.clear();

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;
  
  for(i=0; i<n; i++){
    counter += ref_map.map_read[i];
    r.push_back(counter);
  }
  q.push_back(0);
  counter=0;
  for(j=0; j<m; j++){
    counter += target_map.map_read[j];
    q.push_back(counter);
  }

  vector<vector< bool > > mminf_matrix; //true if score is mminus inf
  
  for(i=0; i<=m; i++){
    vector<double> z;
    vector<int> mminone;
    vector<bool> bmminf;
    for(j=0; j<=n; j++){
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    T.push_back(z);
    S.push_back(z);
    mminf_matrix.push_back(bmminf);
  }
  

  //initialize S and inf matrix
  for(i=0; i<=int_min(end_delta,m); i++){
    for(j=0; j<=n; j++){
      T[0][j] = 0;
      S[0][j] = 0;
      mminf_matrix[0][j] = false; //has been initialized
    }
  }
  for(i=0; i<=m; i++){
    for(j=0; j<=int_min(end_delta,n); j++){
      T[i][j] = 0;
      S[i][j] = 0;
      mminf_matrix[i][j] = false; //has been initialized
    }
  }

    
  //cout<<"starting"<<endl;
  
  //end of matrix initialization

  //begin DP computation
  for(i=1; i<=m; i++){
    for(j=1; j<=n; j++){
      double y;
      //double t_score;
 
      bool y_minf = mminf_matrix[i][j];
      for(g=int_max(0, i-delta); g<i; g++){
	for(h=int_max(0, j-delta); h<j; h++){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    assert(ref_gap > 0);
	    assert(map_gap > 0);	  

	    double s = S[g][h] +
	      score_pars.opt_size_score
	      (q[i]-q[g],r[j]-r[h],map_gap,ref_gap);

	    {  
	      if(y_minf == true){
		y = s;
		fromi[i][j] = g;
		fromj[i][j] = h;	    
		y_minf = false; //y has been set
	      }
	      else{
		if(y < s){
		  fromi[i][j] = g;
		  fromj[i][j] = h;
		  y = s;	    //better choice of y
		  //t_score = cur_t;
		}
	      }
	    }
	  }	  
	}
      }

      if(y_minf == false){
	if(fromi[i][j] != -1 && fromj[i][j] != -1){
	  double cur_t_score = T[fromi[i][j]][fromj[i][j]] +
	    score_pars.opt_site_score(i-fromi[i][j], j-fromj[i][j]);
	  
	  //double cur_t_score = T[fromi[i][j]][fromj[i][j]] +
	  //  score_pars.site_opt_match_score(i-fromi[i][j], j-fromj[i][j]);
	  if(cur_t_score >= local_t_score_thresh &&
	     y >= local_score_thresh){
	    S[i][j] = y;
	    T[i][j] = cur_t_score;
	    
	    mminf_matrix[i][j] = false;	
	  }
	}
      }     
    }
  }
  //end of DP computation
  

  //search for the optimal alignment
  bool opt_score_minf = true;
  //double Tmmax;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(j=0; j<=n; j++){
    for(i=int_max(0, m-end_delta); i<=m; i++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][j] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  for(int i=0; i<=m; i++){
    for(int j=int_max(0, n-end_delta); j<=n; j++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][n] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  assert(opt_score_minf == false);
  //end of optimal alignment search

  Smax = Smmax; 
  Tmax = T[immax][jmmax];
  //trace back procedure to recover the full alignment
  bool end_found = false;
  
  int curposi = immax;
  int curposj = jmmax;

  while(end_found == false){
    tar_restr_al_sites.push_back(curposi);
    ref_restr_al_sites.push_back(curposj);

    s_scores.push_back(S[curposi][curposj]);

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
}

void rm_alignment::optimized_fit_alignment(){

  int end_delta = 0;
  //double lambda = score_pars.lambda;
  int delta = score_pars.delta;

  int m = target_map.map_read.size();
  int n = ref_map.map_read.size();
  int i, j, g, h, l;
 
  double local_score_thresh = -10;
  double local_t_score_thresh = -5;

  T.clear();
  S.clear();
  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  s_scores.clear();

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;
  
  for(i=0; i<n; i++){
    counter += ref_map.map_read[i];
    r.push_back(counter);
  }
  q.push_back(0);
  counter=0;
  for(j=0; j<m; j++){
    counter += target_map.map_read[j];
    q.push_back(counter);
  }

  vector<vector< bool > > mminf_matrix; //true if score is mminus inf
  
  for(i=0; i<=m; i++){
    vector<double> z;
    vector<int> mminone;
    vector<bool> bmminf;
    for(j=0; j<=n; j++){
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    T.push_back(z);
    S.push_back(z);
    mminf_matrix.push_back(bmminf);
  }
  

  //initialize S and inf matrix
  for(i=0; i<=int_min(end_delta,m); i++){
    for(j=0; j<=n; j++){
      T[0][j] = 0;
      S[0][j] = 0;
      mminf_matrix[0][j] = false; //has been initialized
    }
  }
  for(i=0; i<=m; i++){
    for(j=0; j<=int_min(end_delta,n); j++){
      T[i][j] = 0;
      S[i][j] = 0;
      mminf_matrix[i][j] = false; //has been initialized
    }
  }

    
  //cout<<"starting"<<endl;
  
  //end of matrix initialization

  //begin DP computation
  for(i=1; i<=m; i++){
    for(j=1; j<=n; j++){
      double y;
      //double t_score;
 
      bool y_minf = mminf_matrix[i][j];
      for(g=int_max(0, i-delta); g<i; g++){
	for(h=int_max(0, j-delta); h<j; h++){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    assert(ref_gap > 0);
	    assert(map_gap > 0);	  

	    double s = S[g][h] +
	      score_pars.ref_size_score_high
	      (q[i]-q[g],r[j]-r[h],map_gap);
	      //score_pars.opt_size_score
	      //(q[i]-q[g],r[j]-r[h],map_gap,ref_gap);

	    {  
	      if(y_minf == true){
		y = s;
		fromi[i][j] = g;
		fromj[i][j] = h;	    
		y_minf = false; //y has been set
	      }
	      else{
		if(y < s){
		  fromi[i][j] = g;
		  fromj[i][j] = h;
		  y = s;	    //better choice of y
		  //t_score = cur_t;
		}
	      }
	    }
	  }	  
	}
      }

      if(y_minf == false){
	if(fromi[i][j] != -1 && fromj[i][j] != -1){
	  //double cur_t_score = T[fromi[i][j]][fromj[i][j]] +
	  //  score_pars.ref_site_score
	  //  (i-fromi[i][j], j-fromj[i][j], 
	  //  r[j]-r[fromj[i][j]]);
	  double cur_t_score = T[fromi[i][j]][fromj[i][j]] +
	    score_pars.opt_site_score(i-fromi[i][j], j-fromj[i][j]);
	  
	  if(cur_t_score >= local_t_score_thresh &&
	     y >= local_score_thresh){
	    S[i][j] = y;
	    T[i][j] = cur_t_score;
	    
	    mminf_matrix[i][j] = false;	
	  }
	}
      }     
    }
  }
  //end of DP computation
  

  //search for the optimal alignment
  bool opt_score_minf = true;
  //double Tmmax;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(j=0; j<=n; j++){
    for(i=int_max(0, m-end_delta); i<=m; i++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][j] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  for(int i=0; i<=m; i++){
    for(int j=int_max(0, n-end_delta); j<=n; j++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][n] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  assert(opt_score_minf == false);
  //end of optimal alignment search

  Smax = Smmax; 
  Tmax = T[immax][jmmax];
  //trace back procedure to recover the full alignment
  bool end_found = false;
  
  int curposi = immax;
  int curposj = jmmax;

  while(end_found == false){
    tar_restr_al_sites.push_back(curposi);
    ref_restr_al_sites.push_back(curposj);

    s_scores.push_back(S[curposi][curposj]);

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
}

void rm_alignment::optimized_local_ref_alignment(){

  int end_delta = 0;
  //double lambda = score_pars.lambda;
  int delta = score_pars.delta;

  int m = target_map.map_read.size();
  int n = ref_map.map_read.size();
  int i, j, g, h, l;
 
  double local_score_thresh = 0;//-10;
  double local_t_score_thresh = 0;//-5;

  double local_t_score_drop_thresh = -5;

  T.clear();
  S.clear();

  vector<vector<double> > Max_T;

  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  s_scores.clear();

  vector<double> q;
  vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;
  
  for(i=0; i<n; i++){
    counter += ref_map.map_read[i];
    r.push_back(counter);
  }
  q.push_back(0);
  counter=0;
  for(j=0; j<m; j++){
    counter += target_map.map_read[j];
    q.push_back(counter);
  }

  vector<vector< bool > > mminf_matrix; //true if score is mminus inf

  for(i=0; i<=m; i++){
    vector<double> z;
    vector<int> mminone;
    vector<bool> bmminf;
    for(j=0; j<=n; j++){
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(false);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    T.push_back(z);
    S.push_back(z);
    Max_T.push_back(z);
    mminf_matrix.push_back(bmminf);
  }
  

  //initialize S and inf matrix
  /*
  for(i=0; i<=int_min(end_delta,m); i++){
    for(j=0; j<=n; j++){
      T[0][j] = 0;
      S[0][j] = 0;
      mminf_matrix[0][j] = false; //has been initialized
    }
  }
  for(i=0; i<=m; i++){
    for(j=0; j<=int_min(end_delta,n); j++){
      T[i][j] = 0;
      S[i][j] = 0;
      mminf_matrix[i][j] = false; //has been initialized
    }
  }
  */
  //end of matrix initialization

  //begin DP computation
  for(i=1; i<=m; i++){
    for(j=1; j<=n; j++){
      double y;
      //double t_score;
      int besti=-1;
      int bestj=-1;
      bool y_minf = true;
      for(g=int_max(0, i-delta); g<i; g++){
	for(h=int_max(0, j-delta); h<j; h++){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    assert(ref_gap > 0);
	    assert(map_gap > 0);	  

	    double s = S[g][h] +
	      score_pars.ref_size_score_high
	      (q[i]-q[g],r[j]-r[h],map_gap);

	    {  
	      if(y_minf == true){
		y = s;
		besti = g;
		bestj = h;	    
		y_minf = false; //y has been set
	      }
	      else{
		if(y < s){
		  besti = g;
		  bestj = h;
		  y = s;	    //better choice of y
		}
	      }
	    }
	  }	  
	}
      }

      if(y_minf == false){
	if(besti != -1 && bestj != -1){
	  double cur_t_score = T[besti][bestj] +
	    score_pars.opt_site_score(i-besti, j-bestj);
	  double cur_t_drop = cur_t_score - Max_T[besti][bestj];

	  if(cur_t_score >= local_t_score_thresh &&
	     y >= local_score_thresh && 
	     cur_t_drop >= local_t_score_drop_thresh){
	    S[i][j] = y;
	    T[i][j] = cur_t_score;	    
	    mminf_matrix[i][j] = false;
	    fromi[i][j] = besti;
	    fromj[i][j] = bestj;
	    if(cur_t_score >  Max_T[besti][bestj]) Max_T[i][j] = cur_t_score;
	    else Max_T[i][j] =  Max_T[besti][bestj];
	  }
	}
      }
    }
  }
  //end of DP computation
  

  //search for the optimal alignment
  //bool opt_score_minf = true;
  //double Tmmax;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(i=0; i<=m; i++){
    for(j=0; j<=n; j++){
      if(S[i][j] > Smmax){
	immax = i;
	jmmax = j;
	Smmax = S[i][j];
      }
    }
  }

  /*
  for(j=0; j<=n; j++){
    for(i=int_max(0, m-end_delta); i<=m; i++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][j] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  for(int i=0; i<=m; i++){
    for(int j=int_max(0, n-end_delta); j<=n; j++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][n] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  assert(opt_score_minf == false);
  */
  //end of optimal alignment search

  Smax = Smmax; 
  Tmax = T[immax][jmmax];
  //trace back procedure to recover the full alignment
  bool end_found = false;
  
  int curposi = immax;
  int curposj = jmmax;

  while(end_found == false){
    tar_restr_al_sites.push_back(curposi);
    ref_restr_al_sites.push_back(curposj);

    s_scores.push_back(S[curposi][curposj]);

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
}


rm_alignment::rm_alignment(om_read& rm, om_read& tm, scoring_params& sp){
  ref_map.map_read = rm.map_read;
  ref_map.read_name = rm.read_name;
  target_map.map_read = tm.map_read;
  target_map.read_name = tm.read_name;
  
  score_pars = sp;
}
