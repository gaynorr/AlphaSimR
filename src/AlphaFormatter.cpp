/*
 * Adapted from AlphaFormatter.exe
 * 
 * Originally based on msformatter which can be found at 
 * https://github.com/gchen98/macs.
 */

#include <Rcpp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>

using namespace std;

const string MUTATIONSITE="SITE:";
const string SNPBEGIN="BEGIN_SELECTED_SITES";
const string TEMPFILE_PREFIX="haplo.";
const string SNPEND="END_SELECTED_SITES";
const string TOTALSAMPLES="TOTAL_SAMPLES:";
const string TOTALSITES="TOTAL_SITES:";
const string FIELD_DELIMITER="\t";

// [[Rcpp::export]]
int AlphaFormatter(){
  int totalpos;
  int totalpersons;
  ostringstream oss;
  srand(time(NULL));
  double r = rand();
  oss<<TEMPFILE_PREFIX<<r;
  string tempfile=oss.str();
  ofstream haploDataFile(tempfile.data()); 
  double * positions = NULL;
  bool ** genotypes = NULL;
  ofstream PhysicalMapInput;
  ofstream SegSites;
  ofstream MacsHaplotypes;
  ifstream inFile; 
  // nHaplos should be first argument
  int nHaplos = 0;    
  
  // Files to output data to
  PhysicalMapInput.open("PhysicalMapInput.txt");
  MacsHaplotypes.open("MacsHaplotypes.txt");  
  SegSites.open("SegSites.txt");
  inFile.open("output.txt");
  
  if (!inFile) {
    cout<<"ERROR- alphaformatter- output.txt not found.\n";
    return -1;
  }
  
  if (!MacsHaplotypes) {
    cout<<"ERROR- alphaformatter- MacsHaplotypes.txt cannot be created.\n";
    return -1;
  }
  
  if (!SegSites) {
    cout<<"ERROR- alphaformatter- SegSites.txt cannot be created.\n";
    return -1;
  }
  
  if (!PhysicalMapInput) {
    cout<<"ERROR- alphaformatter- PhysicalMapInput.txt cannot be created.\n";
    return -1;
  }
  
  
  string line,constant;
  while(getline(inFile,line)){
    istringstream tokens(line);
    tokens>>constant;
    if(constant.compare(MUTATIONSITE)==0){
      string index,position,mutationTime,mutations;
      tokens>>index>>position>>mutationTime>>mutations;
      haploDataFile<<index<<FIELD_DELIMITER<<position<<FIELD_DELIMITER<<mutations<<endl;
    }else if (constant.compare(TOTALSAMPLES)==0){
      haploDataFile.close();
      tokens>>totalpersons;
      // checks for atoi errors or if specified nHaplos is to large
      if (nHaplos == 0 || (nHaplos > totalpersons)) {
        nHaplos = totalpersons;
      }
    }else if (constant.compare(TOTALSITES)==0){
      tokens>>totalpos;
    }else if (line.compare(SNPBEGIN)==0){          
      SegSites<<totalpos<<endl; 
      try {         
        positions = new double[totalpos];
        genotypes = new bool * [nHaplos];
        for (int i=0;i<totalpersons;++i){
          genotypes[i] = new bool[totalpos];
        }
        
      }
      catch (const std::bad_alloc&) {
        cerr<< "Bad Alloc";
        return -1;
      }
      
      ifstream haploDataFile(tempfile.data());
      // read in the positions
      getline(inFile,line);
      istringstream input(line);
      int cur_snp_index = -1;
      double cur_snp_pos = 0.0;
      string cur_haplo;
      for (int i=0;i<totalpos;++i){
        int selectedsnp;
        input>>selectedsnp;
        string haploline;
        do{
          getline(haploDataFile,haploline);
          istringstream haploInput(haploline);
          haploInput>>cur_snp_index>>cur_snp_pos>>cur_haplo;
        }while(selectedsnp!=cur_snp_index);
        positions[i] = cur_snp_pos;
        for (int j=0;j<nHaplos;++j){
          char cur_char = cur_haplo[j];             
          switch(cur_char){
          case '0':
            genotypes[j][i] = 0;
            break;
          case '1':
            genotypes[j][i] = 1;
            break;
          default:
            cerr<<"Error -Unknown characters being given to formatter";
          return -1;
          break;
          }
        }
      }
      haploDataFile.close();
      remove(tempfile.data());
    }
    else if (line.compare(SNPEND)==0){
      for (int i=0;i<totalpos;++i){
        PhysicalMapInput<<" "<< setprecision(14)<<positions[i];
      }
      PhysicalMapInput<<endl;
      
      delete [] positions;
      for (int j=0;j<nHaplos;++j){
        MacsHaplotypes<<" ";
        for (int i=0;i<totalpos;++i){
          MacsHaplotypes<<genotypes[j][i]<<" ";
        }
        MacsHaplotypes<<endl;
      }
      for (int i=0;i<nHaplos;++i){
        genotypes[i] = new bool[totalpos];
      }
      delete [] genotypes;
      haploDataFile.open(tempfile.data());
    }
  }
  
  
  SegSites.close();
  haploDataFile.close();
  MacsHaplotypes.close();
  remove(tempfile.data());
  return 0;
}
