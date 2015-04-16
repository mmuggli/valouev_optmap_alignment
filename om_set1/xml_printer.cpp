using namespace std;

#include "xml_printer.h"

xml_printer::xml_printer(const char* of_str_name){
  remove(of_str_name);
  of_str.open(of_str_name);
  assert(of_str.good());
}
void xml_printer::print_start(){
  assert(of_str.good());
  of_str<<"<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>"<<endl;
  of_str<<"<aligned_maps_document version=\"0.5\">"<<endl;
  of_str<<endl;
}
void xml_printer::print_finish(){
  assert(of_str.good());
  of_str<<"</aligned_maps_document>"<<endl;
  of_str.close();
}
void xml_printer::print_header(){

  assert(of_str.good());
  of_str<<"<experiment>"<<endl;
  of_str<<" <creator>usc_aligner</creator>"<<endl;
  of_str<<" <uuid>00000000-0000-0000-0000-000000000000</uuid>"<<endl;
  of_str<<" <reference_mapset_file>ref.maps</reference_mapset_file>"<<endl;
  of_str<<" <optical_mapset_file>opmaps.maps</optical_mapset_file>"<<endl;
  of_str<<" <alignment_type>fit</alignment_type>"<<endl;
  of_str<<" <allow_danglers>false</allow_danglers>"<<endl;
  of_str<<" <genspect_file></genspect_file>"<<endl;
  of_str<<" <log_file></log_file>"<<endl;
  of_str<<" <verbose>true</verbose>"<<endl;
  of_str<<" <number_alignments></number_alignments>"<<endl;
  of_str<<" <minimum_score></minimum_score>"<<endl;
  of_str<<" <p_value_cutoff_probability></p_value_cutoff_probability>"<<endl;
  of_str<<" <p_value_lambda></p_value_lambda>"<<endl;
  of_str<<" <p_value_K></p_value_K>"<<endl;
  of_str<<" <p_value_parameter_file>NONE</p_value_parameter_file>"<<endl;
  of_str<<" <enzyme_cut_rate></enzyme_cut_rate>"<<endl;
  of_str<<" <small_fragment_size_error></small_fragment_size_error>"<<endl;
  of_str<<" <small_fragment_size></small_fragment_size>"<<endl;
  of_str<<" <medium_fragment_size_error></medium_fragment_size_error>"<<endl;
  of_str<<" <medium_fragment_size></medium_fragment_size>"<<endl;
  of_str<<" <large_fragment_size_error></large_fragment_size_error>"<<endl;
  of_str<<" <missing_cut_penalty></missing_cut_penalty>"<<endl;
  of_str<<" <false_cut_penalty></false_cut_penalty>"<<endl;
  of_str<<" <small_fragment_median></small_fragment_median>"<<endl;
  of_str<<" <maximum_missing_cuts></maximum_missing_cuts>"<<endl;
  of_str<<" <maximum_false_cuts></maximum_false_cuts>"<<endl;
  of_str<<" <create_date>Fri Jan 23 15:01:48 CST 2004</create_date>"<<endl;
  of_str<<"</experiment>"<<endl;
  of_str<<endl;
}

void xml_printer::print_consensus(om_read& ref_map){

  assert(of_str.good());
  
  of_str<<"<restriction_map>"<<endl;
  of_str<<" <type>consensus</type>"<<endl;
  //of_str<<" <name>"<<ref_map.read_name.c_str()<<"</name>"<<endl;
  of_str<<" <name>reference</name>"<<endl;
  of_str<<" <circular>false</circular>"<<endl;
  of_str<<" <orientation>N</orientation>"<<endl;
  of_str<<" <num_frags>"<<ref_map.map_read.size()<<"</num_frags>"<<endl;
  of_str<<" <enzymes>SwaI</enzymes>"<<endl;

  of_str<<" <map_block>";
  for(int i=0; i<ref_map.map_read.size(); i++){
    of_str<<ref_map.map_read[i];
    if(i<ref_map.map_read.size()-1) of_str<<" ";
  }
  of_str<<"</map_block>"<<endl;
  
  of_str<<"</restriction_map>"<<endl;
  of_str<<endl;
}

void xml_printer::print_alignment(rm_alignment& al){
  of_str<<"<restriction_map>"<<endl;
  of_str<<" <type>opmap</type>"<<endl;
  of_str<<" <name>"<<al.target_map.read_name.c_str()<<"</name>"<<endl;
  of_str<<" <circular>false</circular>"<<endl;
  of_str<<" <orientation>N</orientation>"<<endl;
  of_str<<" <num_frags>"<<al.target_map.map_read.size()<<"</num_frags>"<<endl;
  of_str<<" <enzymes>SwaI</enzymes>"<<endl;
  of_str<<" <map_block>";
  
  for(int i=0; i<al.target_map.map_read.size(); i++){
    of_str<<al.target_map.map_read[i];
    if(i<al.target_map.map_read.size()-1)
      of_str<<" ";
  }
  of_str<<"</map_block>"<<endl;

  of_str<<"</restriction_map>"<<endl;
  of_str<<endl;


  assert(of_str.good());
  of_str<<"<map_alignment>"<<endl;
  of_str<<" <uuid>36a0b669-7eb6-4b26-91d8-9e6be5703715</uuid>"<<endl;
  of_str<<" <reference_map>"<<endl;
  of_str<<"  <name>reference</name>"<<endl;
  of_str<<" </reference_map>"<<endl;
  of_str<<" <aligned_map>"<<endl;
  of_str<<"  <name>"<<al.target_map.read_name.c_str();
  of_str<<"</name>"<<endl;
  of_str<<"  <orientation>N</orientation>"<<endl;
  of_str<<" </aligned_map>"<<endl;
  of_str<<" <usc_aligner_score>"<<al.Tmax<<"</usc_aligner_score>"<<endl;
  of_str<<" <count>"<<al.ref_restr_al_sites.size();
  of_str<<"</count>"<<endl;

  for(int i=al.tar_restr_al_sites.size()-1; i>=1; i--){

    int cur_ref_site = al.ref_restr_al_sites[i];
    int cur_tar_site = al.tar_restr_al_sites[i];
    int next_ref_site = al.ref_restr_al_sites[i-1];
    int next_tar_site = al.tar_restr_al_sites[i-1];

    if(cur_ref_site == next_ref_site - 1){
      for(int j=cur_tar_site; j<next_tar_site; j++){
	of_str<<"<f><i>"<<j<<"</i>";
	of_str<<"<l>"<<cur_ref_site<<"</l>";
	of_str<<"<r>"<<cur_ref_site<<"</r></f>";
	of_str<<endl;
      }
    }
    else{
      if(cur_tar_site == next_tar_site - 1){
	of_str<<"<f><i>"<<cur_tar_site<<"</i>";
	of_str<<"<l>"<<cur_ref_site<<"</l>";
	of_str<<"<r>"<<next_ref_site-1<<"</r></f>";
	of_str<<endl;
      }
    }
    

  }

  of_str<<"</map_alignment>"<<endl;
  of_str<<endl;
}

xml_printer::~xml_printer(){
}
