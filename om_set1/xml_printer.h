#define _XML_PRINTER_H

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

class xml_printer{
  ofstream of_str;
 public:
  xml_printer(const char* of_str_name);
  ~xml_printer();
  void print_start(); //prints pre-header
  void print_header();
  void print_consensus(om_read & ref_map);
  void print_alignment(rm_alignment& al);
  void print_finish();
};
