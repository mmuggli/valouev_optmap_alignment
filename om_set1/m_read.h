// om_read.h: interface for the om_read class.
//
//////////////////////////////////////////////////////////////////////
#define _M_READ_H

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

#ifndef _GLIBCXX_SSTREAM
#include <sstream>
#endif

class om_read{

 public:
  string read_name;
  string Enz_name;
  string Enz_acr;
  vector<double> map_read; //stores the array of rf lengths

  //int index;
  
  om_read() {};
  om_read(string r_name, string r, int ind = 0); //reads the om_read from the lines from the file
  om_read(vector<double> &read, int id, string name, string Ename , string Eacr);
  om_read & operator=(const om_read &read);
  void print();
  void save(const char* f_name);
  void save(ofstream& of_str);
  om_read reverse(); //returns inverse read
  void revert();
  void trimm_start(int num);
  void trimm_end(int num);
  void erase_short_fragments(double thresh);
  double av_size();
  double total_size();//returns total size in Kb
 private:
  bool num_symbol(char c){
    if(c=='0' || c=='1' || c=='2' || c=='3' || c=='4' || c=='5' ||
       c=='6' || c=='7' || c=='8' || c=='9' || c=='.') return true;
    else return false;
  }
};

class om_read_collection{
 public:
  vector<om_read> collection;
  
  om_read_collection(ifstream& fstr);
  om_read_collection();
  void load(ifstream& fstr);
  void output_lengths();
  void report_identical_reads();
  
  void trimm(int start, int end);
  double av_size();
  double est_ref_av_size(double dig_p, double sigma, double zeta);
  inline double sqr(double x){return x*x;};
  void erase_short_fragments(double thresh);
};

//#endif // !defined(AFX_OM_READ_H__15BB8AD0_CBF5_4260_96BB_11384BCB7F1B__INCLUDED_)
