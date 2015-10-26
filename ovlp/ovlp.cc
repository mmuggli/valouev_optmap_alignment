#include <iostream>
#include <fstream>
#include <string>

#include <vector>
#include <algorithm>
#include <math.h>
#include <time.h>

using namespace std;

#include "./../om_set1/msfl.cpp"
#include "./../om_set1/m_read.cpp"
#include "./../om_set1/scoring.cpp"
#include "./../om_set1/alignment.cpp"


int main(int argc, char *argv[])
{
  if(argc < 4){
    cerr<<"Usage:"<<endl;
    cerr<<"ovlp <Map_File> <Ouput> <Range num1> <Range num2>"<<endl;
  }
  else{
    const char* detailed_ovlps_filename = "./detailed.ovlps";
    
    remove(detailed_ovlps_filename);
    ofstream det_str(detailed_ovlps_filename);
    assert(det_str.good());

    //ifstream ifs("./../datasets/pestis/rough_map_2");
    //ifstream ifs("./../datasets/pestis/rough_maps.v2");
    //ifstream ifs("./../datasets/pestis/YPestis.kim.rfmap");
    //ifstream ifs("./../datasets/pestis/rm2.cons");
    //ifstream ifs("./../datasets/Ypestis_XhoI.maps.ex"); //file with the om data
    //ifstream ifs("./../datasets/pestis/Ypestis_XhoI.maps_original");
    //ifstream ifs("./../datasets/susan/chr19/mole_candidate_cut5_chr19.maps");
    //ifstream ifs("./../datasets/pestis/Ypestisoptical.maps");
    //ifstream ifs("./../datasets/Ypestis_XhoIcontigged.maps.short");
    //ifstream ifs("./../datasets/Ypestis_XhoIcontigged.maps");
    //ifstream ifs("./datasets/Ypestis.perm");
    //ifstream ifs("./datasets/Ypestis_XhoIcontigged.maps.77.nonflip");
    //ifstream ifs("./datasets/ch13.om.sim");
    //ifstream ifs("./../datasets/pfluor/Pfluorescensoptical.cut");
    ifstream ifs;
    ifs.open(argv[1]);
    if(!ifs.good()){
      cerr<<"wrong maps file name: "<<argv[1]<<endl;
      assert(false);
    }
    //assert(ifs.good());
    
    om_read_collection maps(ifs);

    ifs.close();

    //ifstream ref_if("./../datasets/pestis/YPestis.kim.rfmap.ext");
    //ifstream ref_if("./../datasets/susan/chr19/chr19_silico.maps");
    //ifstream ref_if("./datasets/YPestis.kim.rfmap");
    //ifstream ref_if("./../datasets/YPestis.kim.rfmap.ex");
    //ifstream ref_if("./datasets/ch13.ref");
    //assert(ref_if.good());

    remove(argv[2]);
ofstream times;
times.open("valuev_row_times.csv");
    ofstream ovlps;
    ovlps.open(argv[2]);
    

    int start_num = atoi(argv[3]);
    int end_num = atoi(argv[4]);

    start_num = 0;
    end_num = maps.collection.size();

    start_num = 0;
    end_num = maps.collection.size();

    assert(start_num <= end_num);
    assert(end_num <= maps.collection.size());
    assert(start_num >= 0);

    scoring_params sp(.2,1.2,.9,3,17.43,0.58, 0.0015, 0.8, 1, 3);
    //sp.init();

    int stored = 0;

    for(int i=start_num; i<end_num; i++){
      clock_t start_time = clock();
      //for(int j=0; j<i/*maps.collection.size()*/; j++){
      for(int j=i; j<end_num/*maps.collection.size()*/; j++){
	if(i!=j && maps.collection[i].map_read.size() > 5
	   && maps.collection[j].map_read.size() > 5){

        if (j % 1000 == 0 ) {
            cout<<endl;
	  
            cout<<"stored: "<<stored<<endl;

            cout<<"start:"<<start_num<<" end:"<<end_num;
            cout<<" cur: "<<i<<"->"<<j<<endl;
        }
	  om_read tar_map = maps.collection[i];
	  om_read for_map = maps.collection[j];
	  om_read rev_map = for_map.reverse();

	  rm_alignment for_alignment(tar_map, for_map, sp);
	  rm_alignment rev_alignment(tar_map, rev_map, sp);

	  //for_alignment.localized_overlap_alignment
	  //  (0,tar_map.map_read.size(),0,for_map.map_read.size());
	  //rev_alignment.localized_overlap_alignment
	  //  (0,tar_map.map_read.size(),0,rev_map.map_read.size());

	  for_alignment.optimized_overlap_alignment();
	  rev_alignment.optimized_overlap_alignment();

	  //for_alignment.overlap_alignment();
	  //rev_alignment.overlap_alignment();

	  for_alignment.overlap_t_score();
	  rev_alignment.overlap_t_score();
	      
	  double for_score = for_alignment.Smax;
	  double rev_score = rev_alignment.Smax;

	  double for_t_score = for_alignment.Tmax;
	  double rev_t_score = rev_alignment.Tmax;

	  //double for_p_value = for_alignment.ovlp_p_value();
	  //double rev_p_value = rev_alignment.ovlp_p_value();

	  double for_ovlp_size = for_alignment.ovlp_size();
	  double rev_ovlp_size = rev_alignment.ovlp_size();

	  // cout<<"fs: "<<for_score<<" ft: "<<for_t_score;
	  // //cout<<" p_v: "<<for_p_value;
	  // cout<<endl;
	  // cout<<"rs: "<<rev_score<<" rt: "<<rev_t_score;
	  // //cout<<" p_v: "<<rev_p_value;
	  // cout<<endl;

	  //rev_alignment.output_alignment(cout);

	  double score_thresh = 25;// 21; // originally 25
	  double t_score_thresh = 8;//7; // originally 8
	  double t_mult = 0;

	  if(for_score > rev_score){
	    //for_alignment.output_alignment(cout);
	  }
	  else{
	    //rev_alignment.output_alignment(cout);
	  }

	  if(for_score > rev_score &&
	     for_t_score > t_score_thresh &&
	     for_score > score_thresh){
	    //&& for_t_score > t_mult*for_ovlp_size){
	    stored++;
	    int ovlp_start1 = 
	      for_alignment.ref_restr_al_sites
	      [for_alignment.ref_restr_al_sites.size()-1];
	    int ovlp_end1 = for_alignment.ref_restr_al_sites[0];

	    int ovlp_start2 = 
	      for_alignment.tar_restr_al_sites
	      [for_alignment.tar_restr_al_sites.size()-1];
	    int ovlp_end2 = for_alignment.tar_restr_al_sites[0];

	    ovlps<<tar_map.read_name.c_str();	    
	    ovlps<<" "<<for_map.read_name.c_str();
	    ovlps<<" "<<tar_map.map_read.size();
	    ovlps<<" "<<for_map.map_read.size();
	    ovlps<<" 1 1 "<<for_score;
	    ovlps<<" "<<for_t_score<<endl;//" ";
	    for(int k=for_alignment.ref_restr_al_sites.size()-1; k>=0; k--){
	      if(k!=for_alignment.ref_restr_al_sites.size()-1)
		ovlps<<" ";
	      ovlps<<for_alignment.ref_restr_al_sites[k];
	      ovlps<<" ";
	      ovlps<<for_alignment.tar_restr_al_sites[k];
	    }
	    //ovlps<<tar_map.map_read.size()<<" ";
	    //ovlps<<for_map.map_read.size()<<" ";
	    //ovlps<<ovlp_start1<<" "<<ovlp_end1<<" ";
	    //ovlps<<ovlp_start2<<" "<<ovlp_end2<<endl;
	    ovlps<<endl<<endl;

	    //if(for_score>score_thresh)
	    //if(for_t_score > t_score_thresh)
	    for_alignment.output_alignment(cout);
	    for_alignment.output_alignment(det_str);
	  }
	  if(for_score <= rev_score && 
	     rev_t_score > t_score_thresh &&
	     rev_score > score_thresh ){
	    //&& rev_t_score > t_mult*rev_ovlp_size){
	    stored++;
	    int rev_map_size = rev_map.map_read.size();
	    int ovlp_start1 = 
	      rev_alignment.ref_restr_al_sites
	      [rev_alignment.ref_restr_al_sites.size()-1];
	    int ovlp_end1 = rev_alignment.ref_restr_al_sites[0];
	    
	    int ovlp_start2 = 
	      rev_alignment.tar_restr_al_sites
	      [rev_alignment.tar_restr_al_sites.size()-1];
	    int ovlp_end2 = rev_alignment.tar_restr_al_sites[0];

	    assert(ovlp_start2 >=0 && ovlp_end2 >= 0);


	    ovlps<<tar_map.read_name.c_str();	    
	    ovlps<<" "<<rev_map.read_name.c_str();
	    ovlps<<" "<<tar_map.map_read.size();
	    ovlps<<" "<<rev_map.map_read.size();
	    ovlps<<" 1 0 "<<rev_score;
	    ovlps<<" "<<rev_t_score<<endl;//" ";
	    for(int k=rev_alignment.ref_restr_al_sites.size()-1; k>=0; k--){
	      if(k!=rev_alignment.ref_restr_al_sites.size()-1)
		ovlps<<" ";
	      ovlps<<rev_alignment.ref_restr_al_sites[k];
	      ovlps<<" ";
	      ovlps<<rev_alignment.tar_restr_al_sites[k];
	    }
	    //ovlps<<tar_map.map_read.size()<<" ";
	    //ovlps<<(rev_map.map_read.size())<<" ";
	    //ovlps<<ovlp_start1<<" "<<ovlp_end1<<" ";
	    //ovlps<<ovlp_start2<<" "<<ovlp_end2<<endl;
	    ovlps<<endl<<endl;
	    	    
	    //if(rev_score > score_thresh)
	    //if(rev_t_score > t_score_thresh)
	    rev_alignment.output_alignment(cout);
	    rev_alignment.output_alignment(det_str);
	  }
	}
      }
      clock_t end_time = clock();
      times << i << ", " << end_time - start_time  << std::endl;
    }
    ovlps.close();
    det_str.close();
    times.close();
  }

  return 0;  
}
  
