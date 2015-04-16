
using namespace std;

#include "m_read.h"

double om_read::av_size(){
  int size = map_read.size();
  double res = 0;

  for(int i=0; i<size; i++){
    double cur_fr = map_read[i];
    res += cur_fr/size;
  }
  return res;
}
double om_read::total_size(){
  double res = 0;
  for(int i=0; i<map_read.size(); i++){
    res += map_read[i];
  }
  return res;
}
void om_read::erase_short_fragments(double thresh){
  assert(thresh >= 0);
  vector<double> new_map;
  
  assert(map_read.size() > 0);

  int i;
  for(i=0; i<map_read.size(); i++){
    double cur_fr = map_read[i];
    if(cur_fr > thresh) new_map.push_back(cur_fr);
  }
  assert(!new_map.empty());
  map_read = new_map;
}
void om_read::revert(){
  vector<double> new_read;
  for(int i=map_read.size()-1; i>=0 ; i--){
    double cur_fr = map_read[i];
    new_read.push_back(cur_fr);
  }
  assert(new_read.size()==map_read.size());
  map_read = new_read;
}
void om_read::trimm_start(int num){
  assert(num>=0);
  if( num > 0){
    int trim;
    int size = map_read.size();
    
    if(num > size) trim = size;
    else trim = num;
    
    vector<double>::iterator it = map_read.begin();
    
    map_read.erase(it, it+num);
  }
}

void om_read::trimm_end(int num){
  assert(num>=0);
  if(num>0){
    int trim;
    int size = map_read.size();
    
    if(num > size) trim = size;
    else trim = num;
    
    vector<double>::iterator it = map_read.end();
    map_read.erase(it-trim, it);
  }
}

void om_read::save(ofstream& of_str){
  assert(of_str.good());
  of_str<<read_name<<endl;
  int counter = 0;
  //fstr<<Enz_name<<" "<<counter<<endl;
  of_str<<char(9)<<Enz_name<<char(9)<<Enz_acr<<char(9);

  int i;//indexer
  for(i=0; i<map_read.size();i++){
    double cur_length = map_read[i];

    char number[20];
    sprintf(number,"%4.2f",(float)cur_length);
    of_str<<number<<char(9);
    //of_str.precision(3);
    //of_str<<cur_length<<char(9);
  }
  of_str<<endl;
  of_str<<endl;
}
/*
void om_read::save(const char* f_name){
  remove(f_name);
  ofstream fstr(f_name);
  fstr<<read_name<<endl;
  int counter = 0;
  fstr<<Enz_name<<" "<<counter<<endl;
  for(int i=0; i<map_read.size(); i++){
    double cur_length = map_read[i];
    //int int_length = (int) 1000*cur_length;
    //counter += int_length;
    fstr<<Enz_name<<" "<<counter<<endl;
  }
  fstr<<"/"<<"/";
  fstr.close();
}
*/

void om_read::print(){
  cout<<"printing the rf_read:"<<endl;
  cout<<read_name.c_str()<<endl;
  cout<<"size:"<<map_read.size()<<endl;
  for(int i=0; i<map_read.size(); i++){
    cout<<map_read[i]<<" ";
  }
  
  cout<<endl;
  cout<<"total size:"<<total_size()<<endl;
  cout<<"end of rf_read printout"<<endl;
}

om_read::om_read(string r_name, string read, int ind){
  read_name = r_name;
  //index = ind;
  
  map_read.clear();

  istringstream cur_stream;//(read);
  cur_stream.str(read);
  
  cur_stream>>Enz_name;
  cur_stream>>Enz_acr;

  int last_pos = cur_stream.tellg();
  while(cur_stream.good()){
    double cur_fr;
    cur_stream>>cur_fr;

    int cur_pos = cur_stream.tellg();
    //if(cur_stream.good())
    if(cur_pos > last_pos){
      map_read.push_back(cur_fr);      
    }
    last_pos = cur_pos;
  }

  /*
  //cout<<"outputting string by letters"<<endl;
  string::iterator it_read = read.begin();
  int length = read.length();
  
  if((int)(*it_read)==9){
    //cout<<"tab in the beginning of the read name"<<endl;
    //cout<<r_name<<endl;
  }
  assert((int)(*it_read)==9);
  it_read++;//skip the tab;


  //read the enzyme name
  string e_name;
  while((int)(*it_read) != 9){
    assert(it_read < read.end()); //assert not the end of the line
    e_name = e_name + (*it_read);
    it_read++;
  }
  //cout<<"the enzyme name is "<<e_name<<endl;
  Enz_name = e_name;
  //cout<<Enz_name<<endl;

  //skip the tab
  assert((int)(*it_read) == 9);
  it_read++;

  //read the enzyme acronym
  Enz_acr = (*it_read);

  //skip the tab
  it_read++;
  assert((int)(*it_read) == 9);
  it_read++;

  //read the optical map
  while(it_read < read.end()){
    string cur_num;

    //while((int)(*it_read) != 9){
    while(num_symbol(*it_read)){
      cur_num = cur_num + (*it_read);
      it_read++;
    }
    //cout<<"current number is: "<<cur_num.c_str()<<endl;
    float c_num;
    const char* read_num = cur_num.c_str();
    //cout<<read_num<<endl;
    char r_num[10];

    for(int i=0; i<10; i++){
      //if(is_number(read_num[i]))
	r_num[i]=read_num[i];
	//else r_num[i] = ' ';
    }
    sscanf ((r_num),"%f",&c_num);
    //cout<<"current number is: "<<(c_num)<<endl;
    assert(c_num > 0);

    double double_c_num = (double)c_num;
    map_read.push_back(double_c_num);

    //cout<<"current number is: "<<(double_c_num)<<endl;
    it_read++;
  }
  */
}

om_read::om_read(vector<double> &read, int id, string name, string Ename, string Eacr)
{	
  read_name = name;
  Enz_name = Ename;
  Enz_acr = Eacr;
  //index = id;
  map_read = read;
}

om_read & om_read::operator=(const om_read &read)
{	
  read_name = read.read_name;
  Enz_name = read.Enz_name;
  Enz_acr = read.Enz_acr;
  //index = read.index;
  map_read = read.map_read;
  
  return *this;
}

om_read om_read::reverse(){
  int i;
  vector<double> inv_map;
  for(i=map_read.size()-1; i>=0; i--){
    inv_map.push_back(map_read[i]);
  }

  // append '.r' to the end of the original name. //
  //string new_name = read_name + ".r";
  string new_name = read_name;

  //  om_read new_read(inv_map, index, read_name, Enz_name, Enz_acr);
  int cur_id; //to remove
  
  om_read new_read(inv_map, cur_id, new_name, Enz_name, Enz_acr);
  return new_read;
}

void om_read_collection::trimm(int start, int end){
  assert(start>=0);
  assert(end>=0);
  for(int i=0; i<collection.size(); i++){
    collection[i].trimm_start(start);
    collection[i].trimm_end(end);
  }
}
void om_read_collection::report_identical_reads(){
  int count=0;
  for(int i=0; i<collection.size(); i++){
    for(int j=i+1; j<collection.size(); j++){
      if(collection[i].read_name ==  collection[j].read_name){
	cout<<"found identical reads with name: "
	    <<collection[i].read_name<<endl;
	collection[i].print();
	collection[j].print();
	count++;
      }
    }
  }
  cout<<"found "<<count<<" identical reads"<<endl;
}

om_read_collection::om_read_collection(){
}//default constructor

void om_read_collection::output_lengths(){
  vector<double> lengths;
  int i;
  for(i=0; i<collection.size(); i++){
    for(int j=0; j<collection[i].map_read.size(); j++){
      
      lengths.push_back(collection[i].map_read[j]);
    }
  }
  sort(lengths.begin(), lengths.end());
  
  ofstream ofstr("./datasets/lengths");
  
  for(i=0; i<lengths.size(); i++){
    ofstr<<lengths[i]<<endl;
    //cout<<lengths[i]<<endl;
  }
  ofstr.close();
}

om_read_collection::om_read_collection(ifstream& fstr){
  load(fstr);
}

void om_read_collection::load(ifstream& fstr){
  string line1;
  string line2;
  string line3;
  
  int read_indexer = 0;
  while(!getline(fstr, line1).eof() &&
	!getline(fstr, line2).eof() ){
    //!getline(fstr, line3).eof() ){
    //getline(fstr, line2);
    getline(fstr, line3);//empty line

    //cout<<line1.c_str()<<endl;
    //cout<<line2.c_str()<<endl;
    //cout<<line3.c_str()<<endl;
    if (read_indexer % 100 == 0) cout<<read_indexer<<endl;
    om_read cur_read(line1, line2, read_indexer);

    //bool first_occurance = true; 
    //make sure was not loaded before
    //*
    //for(int i=0; i<collection.size(); i++){
    //  if(cur_read.read_name == collection[i].read_name){
    //	first_occurance = false;
    // }
    //}
    //*/
    //if(first_occurance == true)
    // {
    collection.push_back(cur_read);
    //      }
    
    //cout<<"current read name: "<<cur_read.read_name;
    //cout<<" length: "<<cur_read.map_read.size()<<endl;

    read_indexer++;
  }
  cout<<"the number of reads in the file: "<<collection.size()<<endl;
}

double om_read_collection::av_size(){
  assert(!collection.empty());
  int i;
  double total_frs = 0;
  for(i=0; i<collection.size(); i++){
    total_frs += collection[i].map_read.size();
  }
  double total_av_size = 0;
  for(i=0; i<collection.size(); i++){
    double fr_num = collection[i].map_read.size();
    double cur_av_fr_size = collection[i].av_size();
    total_av_size += cur_av_fr_size * fr_num / total_frs; 
  }

  assert(total_av_size > 0);
  return total_av_size;
}

double om_read_collection::est_ref_av_size(double dig_p, double sigma,
					   double zeta){
  double cur_lambda;
  double theta = this->av_size();
  assert(theta > 0);
  cur_lambda = 2*dig_p/(sqr((1/theta+1/sqr(sigma))*sigma) 
			- 1/sqr(sigma) - 2*zeta);
  return cur_lambda;
}

void om_read_collection::erase_short_fragments(double thresh){
  for(int i=0; i<collection.size(); i++){
    collection[i].erase_short_fragments(thresh);
  }
}
