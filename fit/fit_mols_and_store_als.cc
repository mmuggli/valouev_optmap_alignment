 #include <iostream>
#include <fstream>
#include <string>

#include <vector>
#include <algorithm>
#include <math.h>

using namespace std;

#include "./../om_set1/msfl.cpp"
#include "./../om_set1/m_read.cpp"
#include "./../om_set1/scoring.cpp"
#include "./../om_set1/alignment.cpp"

int main(int argc, char *argv[])
{
  if(argc < 3){
    cerr<<"Usage:"<<endl;
    cerr<<argv[0]<<" <your reference map file> <your optical map file> ";
    return 0;
  }


  //init files
  const char* reference_map_file = argv[1];
  const char* optmap_file = argv[2];
  
  ifstream ifs(optmap_file);
  assert(ifs.good());

  ifstream ref_if(reference_map_file);
  assert(ref_if.good());

  //om_read_collection maps(ifs);
  om_read_collection ref_maps(ref_if);
  //ref_if.close();


  cout<<"reference maps: "<<ref_maps.collection.size()<<endl;
  //cout<<"optical maps: "<<maps.collection.size()<<endl;

  ref_maps.erase_short_fragments(2.0); //erase short fragments if necessary

  scoring_params sp(0.2,2,1,5,17.43,0.579, 0.005, 0.8, 3, 1); //scoring parameters
  //modify here your scoring function

  bool _end_of_map_file = false;
  int counter = 0;
  int good_counter = 0;
  while(!_end_of_map_file){
    string header;
    string frags;
    string dummy;

    bool end1 = getline(ifs, header).eof();
    bool end2 = getline(ifs, frags).eof();
    getline(ifs, dummy);
    
    _end_of_map_file = end1 || end2;

    if(!_end_of_map_file){
      om_read cur_for_opt_map(header, frags, 0);
      om_read cur_rev_opt_map = cur_for_opt_map.reverse();

      //cur_for_opt_map.trimm_start(1);
      //cur_for_opt_map.trimm_end(1);
      //cur_rev_opt_map.trimm_start(1);
      //cur_rev_opt_map.trimm_end(1);

      vector<rm_alignment> for_contig_fits;
      vector<rm_alignment> rev_contig_fits;

      for(int i=0; i<ref_maps.collection.size(); i++){
	rm_alignment cur_for_al(ref_maps.collection[i], cur_for_opt_map, sp);
	rm_alignment cur_rev_al(ref_maps.collection[i], cur_rev_opt_map, sp);

	cur_for_al.fit_alignment();
	cur_rev_al.fit_alignment();

	for_contig_fits.push_back(cur_for_al);
	rev_contig_fits.push_back(cur_rev_al);
      }

      if(!for_contig_fits.empty()){
	double max_score;
	bool for_orient;
	int best_al_ind;

	for(int i=0; i<for_contig_fits.size(); i++){
	  double cur_for_score = for_contig_fits[i].Smax;
	  double cur_rev_score = rev_contig_fits[i].Smax;

	  if(i==0){
	    best_al_ind = i;
	    if(cur_for_score > cur_rev_score){
	      max_score = cur_for_score;
	      for_orient = true;
	    }
	    else{
	      max_score = cur_rev_score;
	      for_orient = false;
	    }
	  }
	  else{
	    if(cur_for_score > max_score){
	      max_score = cur_for_score;
	      for_orient = true;
	      best_al_ind = i;
	    }
	    if(cur_rev_score > max_score){
	      max_score = cur_rev_score;
	      for_orient = false;
	      best_al_ind = i;
	    }
	  }
	}

	//output best alignment
	cout<<"cur_map:"<<counter<<" good:"<<good_counter;
	cout<<" max_score:"<<max_score;
	cout<<endl;

	double min_score_thresh = 10;
	double t_mult = 0;//.1;
	if(max_score > min_score_thresh)
	{
	  double t_score;
	  int fr_num = cur_for_opt_map.map_read.size();	
	    
	  bool store = false;
	  vector<int> al_ref_sites;
	  vector<int> al_map_sites;

	  if(for_orient == true){
	    //store = true;

	    for_contig_fits[best_al_ind].fit_t_score();
	    t_score = for_contig_fits[best_al_ind].Tmax;
	    if(t_score >= ((double)fr_num)*t_mult){
	      store = true;
	      al_ref_sites = for_contig_fits[best_al_ind].ref_restr_al_sites;
	      al_map_sites = for_contig_fits[best_al_ind].tar_restr_al_sites;
	      cout<<"forward:"<<endl;
	      for_contig_fits[best_al_ind].output_alignment(cout);
	      good_counter++;
	    }
	  }
	  else{ 	  
	    rev_contig_fits[best_al_ind].fit_t_score();
	    t_score = rev_contig_fits[best_al_ind].Tmax;
	    if(t_score >= ((double)fr_num)*t_mult){
	      store = true;
	      al_ref_sites = rev_contig_fits[best_al_ind].ref_restr_al_sites;
	      al_map_sites = rev_contig_fits[best_al_ind].tar_restr_al_sites;
	      cout<<"reverse:"<<endl;
	      rev_contig_fits[best_al_ind].output_alignment(cout);
	      good_counter++;
	    }
	  }

	  if(store == true){
	    //store names of maps
	    //out_als<<for_contig_fits[best_al_ind].ref_map.read_name;
	    //out_als<<char(9);
	    //out_als<<cur_for_opt_map.read_name;	    	    
	    //out_als<<char(9)<<1<<char(9)<<for_orient;
	    //out_als<<endl;

	    assert(al_ref_sites.size() == al_map_sites.size());
	    for(int m=al_ref_sites.size()-1; m>=0; m--){
	      if(m == al_ref_sites.size()-1){
	      }
	      else{
		//out_als<<char(9);
	      }
	      //out_als<<al_ref_sites[m]<<char(9)<<al_map_sites[m];
	    }
	    //out_als<<endl<<endl;
	  }
	}
      }
    }
    counter++;
  }

  ifs.close();
  ref_if.close();

  return 0;
}
