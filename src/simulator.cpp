#include <RcppArmadillo.h>
#include <bitset>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <tuple>
#include <math.h>
#include <algorithm> 
#include <cctype>
#include <locale>

#include <boost/algorithm/string/split.hpp> // Include for boost::split
#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include "simulator.h"
#include <boost/algorithm/string/split.hpp> // Include for boost::split
#include "misc.h"

const double Node::MAX_HEIGHT=1e50;

AlphaSimRReturn::AlphaSimRReturn(){}

RandNumGenerator::RandNumGenerator(unsigned long iRandomSeed){
  boost::mt19937 mt(static_cast<unsigned long int>(iRandomSeed));
  unif = new boost::uniform_01<boost::mt19937>(mt);
}

RandNumGenerator::~RandNumGenerator(){
  delete unif;
}

Configuration::Configuration(){
  bVariableRecomb = false;
  bSNPAscertainment = false;
  bFlipAlleles = false;
  bNewickFormat = false;
  bMigrationChangeEventDefined = false;
  dTheta = 0; // mutation rate parameter theta
  dRecombRateRAcrossSites = 0.0; // recombination rate parameter r
  dSeqLength = 0.;
  dGeneConvRatio = 0.;
  iGeneConvTract = 1;
  iTotalPops = 1; // total populations declared
  dBasesToTrack = 1;
  iRandomSeed = time(NULL);
  iIterations = 1;
  pAlleleFreqBinPtrSet = NULL;
  pEventList = NULL;
}

Configuration::~Configuration(){
  if (pEventList) {
    delete pEventList;
  }
  if (bSNPAscertainment){
    AlleleFreqBinPtrSet::iterator it;
    for (it=pAlleleFreqBinPtrSet->begin();it!=pAlleleFreqBinPtrSet->end();++it){
      delete(*it);
    }
    delete pAlleleFreqBinPtrSet;
  }
  if (bVariableRecomb){
    HotSpotBinPtrList::iterator it;
    for(it=pHotSpotBinPtrList->begin();it!=pHotSpotBinPtrList->end();++it){
      delete(*it);
    }
    if (pHotSpotBinPtrList) {
      delete pHotSpotBinPtrList;
    }
  }
}

void Simulator::readInputParameters(CommandArguments arguments){
  unsigned int popId;
  
  int iSampleSize;
  double dDefaultPopSize,dDefaultGrowthAlpha,dDefaultMigrationRate,popSize;
  bool bAcceptFullMigrMatrix;
  
  unsigned int iTotalArgs = arguments.size();
  
  pConfig=new Configuration();
  
  
  dDefaultPopSize = 1.0;
  dDefaultGrowthAlpha =0.0;
  
  
  pConfig->iTotalPops = 1;
  
  EventPtrList * pEventList = new EventPtrList;
  
  dDefaultMigrationRate = 0.0;
  for (unsigned int i=0;i<pConfig->iTotalPops;++i){
    vector<double> newRow;
    for (unsigned int j=0;j<pConfig->iTotalPops;++j)
      newRow.push_back(dDefaultMigrationRate);
    pConfig->dMigrationMatrix.push_back(newRow);
  }
  
  iSampleSize = atoi(arguments[0][0].data());
  Population newPop;
  newPop.setChrSampled(iSampleSize);
  newPop.setPopSize(dDefaultPopSize);
  newPop.setGrowthAlpha(dDefaultGrowthAlpha);
  
  pConfig->pPopList.push_back(newPop);
  pConfig->iSampleSize = iSampleSize;
  
  pConfig->dSeqLength = atof(arguments[0][1].data());
  
  set<float> eventTimes;
  for (unsigned int iCurrentArg = 1;iCurrentArg<iTotalArgs;++iCurrentArg){
    try{
      double dTime;
      char chType;
      EventPtr wrapper;
      short unsigned int iNoMigrPops;
      short unsigned int iMigrPops;
      short unsigned int iTotalCells;
      short unsigned int iSubOption;
      int iRunningSample=0;
      string command;
      const char * filename;
      ifstream inFile;
      bool  flipAllele;
      switch (arguments[iCurrentArg][0][1] ){
      case 'T' :
        pConfig->bNewickFormat = true;
        // example:
        // (2:1.766,(4:0.505,(3:0.222,(1:0.163,5:0.163):0.059):0.283):1.261);
        break;
      case 'h' :
        if (arguments[iCurrentArg].size()!=2) {
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
            ", you must enter a single integer for retaining the number of previous trees\n";
          Rcpp::stop("Argument error");
          return;
        }
        pConfig->dBasesToTrack = atof(arguments[iCurrentArg][1].data());
        break;
      case 's' :
        if (arguments[iCurrentArg].size()<2) {
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
            ", you must enter a single integer for the random seed\n";
          Rcpp::stop("Argument error");
        }
        pConfig->iRandomSeed = atoi(arguments[iCurrentArg][1].data());
        break;
      case 't' :  // set mutation parameter
        if (arguments[iCurrentArg].size()!=2) {
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
            ", you must enter a single float value for the mutation parameter\n";
          Rcpp::stop("Argument error");
        }
        pConfig->dTheta = pConfig->dSeqLength * atof(arguments[iCurrentArg][1].data());
        break;
      case 'F':
        if (arguments[iCurrentArg].size()!=3){
          Rcpp::Rcerr<<"For the SNP ascertainment feature you must enter the filename of the SNP "<<
            "frequency list and whether to flip the alleles."<<endl;
          Rcpp::stop("Argument error");
        }
        filename = arguments[iCurrentArg][1].data();
        flipAllele = pConfig->bFlipAlleles = atoi(arguments[iCurrentArg][2].data());
        inFile.open(filename);
        if (!inFile)
          throw "SNP freq input file not found\n";
        if (inFile.is_open()) {
          pConfig->pAlleleFreqBinPtrSet = new AlleleFreqBinPtrSet;
          string line;
          int total=0;
          double lastStart=0.;
          double cumFreq=0.;
          double maxFreq = flipAllele?0.5:1.;
          while(getline(inFile,line)){
            istringstream inStr(line);
            double start = lastStart,end,freq;
            inStr>>end>>freq;
            cumFreq+=freq;
            if (end>maxFreq) end = maxFreq;
            if (start>=end) throw "The freq range entered is incorrect.";
            AlleleFreqBinPtr bin = AlleleFreqBinPtr(new AlleleFreqBin(start,end,freq));
            pConfig->pAlleleFreqBinPtrSet->insert(bin);
            lastStart = end;
            ++total;
          }
          inFile.close();
          pConfig->bSNPAscertainment = true;
          if (cumFreq>1.0) throw "The total frequency entered exceeds one";
          if (cumFreq<1.0){
            if (lastStart==maxFreq){
            }else{
              AlleleFreqBinPtr bin = AlleleFreqBinPtr(new AlleleFreqBin(lastStart,maxFreq,1-cumFreq));
              pConfig->pAlleleFreqBinPtrSet->insert(bin);
            }
          }
        }
        break;
      case 'r' :
        if (arguments[iCurrentArg].size()!=2) {
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
            ", you must enter a single float value for the recombination parameter\n";
          Rcpp::stop("Argument error");
        }
        pConfig->dRecombRateRAcrossSites = pConfig->dSeqLength * atof(arguments[iCurrentArg][1].data());
        break;
      case 'c' :
        if (arguments[iCurrentArg].size()!=3) {
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
            ", you must enter the conversion to xover ratio followed by the mean tract length in bp.\n";
          Rcpp::stop("Argument error");
        }
        pConfig->dGeneConvRatio = atof(arguments[iCurrentArg][1].data());
        pConfig->iGeneConvTract = atoi(arguments[iCurrentArg][2].data());
        if (pConfig->dGeneConvRatio<0){
          Rcpp::Rcerr<<"The gene conversion parameters must be positive\n";
          Rcpp::stop("Argument error");
        }
        break;
      case 'R':
        if (arguments[iCurrentArg].size()!=2){
          Rcpp::Rcerr<<"For the hotspot feature you must enter the filename of the hotspot list."<<endl;
          Rcpp::stop("Argument error");
        }
        filename = arguments[iCurrentArg][1].data();
        inFile.open(filename);
        if (!inFile)
          throw "HotSpot input file not found\n";
        if (inFile.is_open()) {
          pConfig->pHotSpotBinPtrList = new HotSpotBinPtrList;
          string line;
          int total=0;
          while(getline(inFile,line)){
            istringstream inStr(line);
            double start,end,ratio;
            inStr>>start>>end>>ratio;
            HotSpotBinPtr bin(new HotSpotBin(start,end,ratio));
            pConfig->pHotSpotBinPtrList->push_back(bin);
            ++total;
          }
          inFile.close();
          pConfig->bVariableRecomb = true;
        }
        break;
      case 'i' :
        if (arguments[iCurrentArg].size()!=2) {
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
            ", you must enter a single int value for the number of iterations\n";
          Rcpp::stop("Argument error");
        }
        pConfig->iIterations = atoi(arguments[iCurrentArg][1].data());
        break;
      case 'I' :
        if (arguments[iCurrentArg].size()<2) {
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
            ", the first parameter needs to the number of population islands\n";
          Rcpp::stop("Argument error");
        }
        pConfig->iTotalPops = atoi( arguments[iCurrentArg][1].data());
        iNoMigrPops=2+pConfig->iTotalPops;
        iMigrPops=3+pConfig->iTotalPops;
        
        if (arguments[iCurrentArg].size()==iNoMigrPops||
            arguments[iCurrentArg].size()==iMigrPops) {
          pConfig->pPopList.clear();
          for(unsigned int i=0; i<pConfig->iTotalPops; ++i) {
            Population newPop;
            const char * arg = arguments[iCurrentArg][2+i].data();
            newPop.setChrSampled(atoi(arg));
            iRunningSample+=newPop.getChrSampled();
            //Rcpp::Rcerr<<"INPUT: Setting chr sampled for pop "<<(i+1)<<" to "<<newPop.getChrSampled()<<endl;
            newPop.setPopSize(dDefaultPopSize) ;
            newPop.setGrowthAlpha(dDefaultGrowthAlpha);
            newPop.setLastTime(0);
            pConfig->pPopList.push_back(newPop);
          }
          if (arguments[iCurrentArg].size()==iMigrPops){
            pConfig->dGlobalMigration = atof(arguments[iCurrentArg][iMigrPops-1].data());
          }else{
            pConfig->dGlobalMigration = dDefaultMigrationRate;
          }
        }else{
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
            ", the number of island sample sizes entered does not match the first parameter\n";
          Rcpp::stop("Argument error");
        }
        if (iRunningSample!=iSampleSize){
          throw "The number of chromosomes entered in the -I option doesn't match the total sample size";
        }
        
        // Allocate migration rate matrix
        pConfig->dMigrationMatrix.clear();
        for (int i=0;i<pConfig->iTotalPops;++i){
          vector<double> newRow;
          for (int j=0;j<pConfig->iTotalPops;++j)
            newRow.push_back(pConfig->dGlobalMigration/
              (pConfig->iTotalPops-1));
          pConfig->dMigrationMatrix.push_back(newRow);
        }
        break;
      case 'm' :
        if( pConfig->iTotalPops < 2 ) {
          Rcpp::Rcerr<<"You must use -I option first (i.e. specify more than one population)."<<endl;
          Rcpp::stop("Argument error");
        }
        if (arguments[iCurrentArg][0][2]=='a') {
          iTotalCells = pConfig->iTotalPops * pConfig->iTotalPops + 1;
          if (arguments[iCurrentArg].size()!=iTotalCells){
            Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
              ", the number of matrix cells does not match the total populations squared\n";
            Rcpp::stop("Argument error");
          }
          iSubOption = 0;
          for(int pop1 = 0; pop1 <pConfig->iTotalPops; ++pop1){
            for(int pop2 = 0; pop2 <pConfig->iTotalPops; ++pop2){
              pConfig->dMigrationMatrix[pop1][pop2]=
                atof( arguments[iCurrentArg][++iSubOption].data() ) ;
            }
          }
          // set the diagonals as the sum of the off-diagonals
          for(int pop1 = 0; pop1 < pConfig->iTotalPops; ++pop1) {
            pConfig->dMigrationMatrix[pop1][pop1] = 0.0 ;
            for(int pop2 = 0; pop2 < pConfig->iTotalPops; ++pop2){
              if( pop1 != pop2 )
                pConfig->dMigrationMatrix[pop1][pop1] +=
                  pConfig->dMigrationMatrix[pop1][pop2];
            }
          }
        } else {
          //                    // lets the user enter the entire migration by specified element
          if (arguments[iCurrentArg].size()!=4){
            Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
              ", you need the source pop, dest pop, and the migration rate.\n";
            Rcpp::stop("Argument error");
          }else{
            int i = atoi( arguments[iCurrentArg][1].data() ) -1;
            int j = atoi( arguments[iCurrentArg][2].data() ) -1;
            double mij = atof( arguments[iCurrentArg][3].data() );
            pConfig->dMigrationMatrix[i][i] += mij -
              pConfig->dMigrationMatrix[i][j];
            pConfig->dMigrationMatrix[i][j] = mij;
          }
        }
        break;
      case 'n' :
        //                    // specify population size for each population
        if (arguments[iCurrentArg].size()!=3){
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
            ", you need to specify the pop ID and the population size.\n";
          Rcpp::stop("Argument error");
        }else{
          popId = atoi( arguments[iCurrentArg][1].data() ) -1;
          popSize = atof( arguments[iCurrentArg][2].data() );
          if (popId>pConfig->iTotalPops){
            Rcpp::Rcerr<<"Invalid pop ID"<<endl;
            Rcpp::stop("Argument error");
          }
          pConfig->pPopList[popId].setPopSize(popSize) ;
        }
        break;
      case 'g' :
        //                    // specify growth rates
        if (arguments[iCurrentArg].size()!=3){
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
            ", you need to specify the pop ID and the population growth rate.\n";
          Rcpp::stop("Argument error");
        }else{
          popId = atoi( arguments[iCurrentArg][1].data() ) -1;
          dDefaultGrowthAlpha = atof( arguments[iCurrentArg][2].data() );
          if (popId>pConfig->iTotalPops){
            Rcpp::Rcerr<<"Invalid pop ID"<<endl;
            Rcpp::stop("Argument error");
          }
          pConfig->pPopList[popId].setGrowthAlpha(dDefaultGrowthAlpha);
        }
        break;
      case 'G' :
        //                    // specify growth rates across all populations
        if (arguments[iCurrentArg].size()!=2){
          Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
            ", you need to specify a single growth rate for all populations.\n";
          Rcpp::stop("Argument error");
        }else{
          float g = atof(arguments[iCurrentArg][1].data());
          if (g<0) throw "Global growth rate must be positive";
          dDefaultGrowthAlpha = atof( arguments[iCurrentArg][1].data() );
          //                    Rcpp::Rcerr<<"INPUT: Growth rate for all pop "<<dDefaultGrowthAlpha<<endl;
          for(int i=0; i<pConfig->iTotalPops; ++i){
            pConfig->pPopList[i].setGrowthAlpha(dDefaultGrowthAlpha);
          }
        }
        break;
      case 'e' :
        // these are events.  Be sure the times are unique
        chType = arguments[iCurrentArg][0][2];
        if ((arguments[iCurrentArg][0][3])=='a') bAcceptFullMigrMatrix = true;
        else bAcceptFullMigrMatrix = false;
        if (arguments[iCurrentArg].size()<2){
          Rcpp::Rcerr<<"For event flags, you need to specify at least a time after "<<
            arguments[iCurrentArg][0]<<endl;
          Rcpp::stop("Argument error");
        }
        dTime = atof(arguments[iCurrentArg][1].data());
        if (eventTimes.find(dTime)==eventTimes.end()){
          eventTimes.insert(dTime);
        }else{
          Rcpp::Rcerr<<"Error, this event is redundant with a previous time.  Please increment it slightly from "<<dTime<<" to prevent unpredictable results\n";
          throw "Invalid input";
        }
        int iPop1,iPop2;
        double dProportion;
        switch(chType){
        case 'N': // global population size
          if (arguments[iCurrentArg].size()!=3){
            Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
              ", you need to specify a single pop size for all populations.\n";
            Rcpp::stop("Argument error");
          }else{
            //int iType = Event::GLOBAL_POPSIZE;
            wrapper = EventPtr(new GenericEvent(
              Event::GLOBAL_POPSIZE,dTime,
              atof(arguments[iCurrentArg][2].data())));
            //Rcpp::Rcerr<<"Global pop size is "<<
            //  atof( arguments[iCurrentArg][2].data() )<<endl;
          }
          break;
        case 'G': // global growth rate
          if (arguments[iCurrentArg].size()!=3){
            Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
              ", you need to specify a single growth rate for all populations.\n";
            Rcpp::stop("Argument error");
          }else{
            float g = atof(arguments[iCurrentArg][2].data());
            //if (g<0) throw "Global event growth rate must be positive";
            //int iType = Event::GLOBAL_POPGROWTH;
            wrapper = EventPtr(new GenericEvent(
              Event::GLOBAL_POPGROWTH,dTime,
              g));
            //Rcpp::Rcerr<<"Global growth rate is "<<
            //  g<<endl;
          }
          break;
        case 'M': // global migration rate
          pConfig->bMigrationChangeEventDefined = true;
          if (arguments[iCurrentArg].size()!=3){
            Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
              ", you need to specify a single migration rate for all populations.\n";
            Rcpp::stop("Argument error");
          }else{
            //iType = Event::GLOBAL_MIGRATIONRATE;
            wrapper = EventPtr(new GenericEvent(
              Event::GLOBAL_MIGRATIONRATE,dTime,
              atof(arguments[iCurrentArg][2].data())));
            //Rcpp::Rcerr<<"Global migration rate is "<<
            //  atof( arguments[iCurrentArg][2].data() )<<endl;
          }
          break;
        case 'n' :  // subpopulation size
          if (arguments[iCurrentArg].size()!=4){
            Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
              ", you need to specify pop id followed by the new size.\n";
            Rcpp::stop("Argument error");
          }else{
            //iType = Event::POPSIZE;
            wrapper = EventPtr(new PopSizeChangeEvent(
              Event::POPSIZE,dTime,atoi( arguments[iCurrentArg][2].data() ) -1,
              atof( arguments[iCurrentArg][3].data() )));
            //Rcpp::Rcerr<<"For population "<<arguments[iCurrentArg][2]<<
            //  ", pop size is now "<<atof( arguments[iCurrentArg][3].data() )<<endl;
          }
          break;
        case 'g' :  // subpopulation growth
          if (arguments[iCurrentArg].size()!=4){
            Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
              ", you need to specify pop id followed by the new growth rate.\n";
            Rcpp::stop("Argument error");;
          }else{
            //iType = Event::GROWTH;
            wrapper = EventPtr(new PopSizeChangeEvent(
              Event::GROWTH,dTime,atoi( arguments[iCurrentArg][2].data() ) -1,
              atof( arguments[iCurrentArg][3].data() )));
            //Rcpp::Rcerr<<"For population "<<arguments[iCurrentArg][2]<<
            //  ", pop growth rate is now "<<atof( arguments[iCurrentArg][3].data() )<<endl;
          }
          break;
        case 's' :  // split
          pConfig->bMigrationChangeEventDefined = true;
          if (arguments[iCurrentArg].size()!=4){
            Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
              ", you need to specify pop id followed by the proportion of the split.\n";
            Rcpp::stop("Argument error");
          }else{
            //iType = Event::POPSPLIT;
            iPop1 = atoi( arguments[iCurrentArg][2].data() )-1;
            dProportion = atof( arguments[iCurrentArg][3].data() );
            wrapper = EventPtr(new PopSizeChangeEvent(
              Event::POPSPLIT,dTime,iPop1,dProportion));
            //Rcpp::Rcerr<<"Population "<<arguments[iCurrentArg][2]<<
            //  " splits at proportion "<<dProportion<<endl;
          }
          break;
        case 'j':   // move lineages from pop1 to pop2
          
          if (arguments[iCurrentArg].size()!=4){
            Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
              ", you need to specify source pop id followed by the destination pop id.\n";
            Rcpp::stop("Argument error");
          }else{
            //iType = Event::POPJOIN;
            iPop1 = atoi( arguments[iCurrentArg][2].data() ) -1;
            iPop2 = atoi( arguments[iCurrentArg][3].data() ) -1;
            if (iPop1>=pConfig->iTotalPops||
                iPop2>=pConfig->iTotalPops){
              Rcpp::Rcerr<<"WARNING: The pop IDs used in pop join is greater than the number specified in -I.  You must have a split event before this join event.\n";
            }
            wrapper = EventPtr(new PopJoinEvent(
              Event::POPJOIN,dTime,iPop1,iPop2));
            //Rcpp::Rcerr<<"Population "<<
            //  arguments[iCurrentArg][2].data()<<
            //    " will merge with "<<
            //      arguments[iCurrentArg][3].data()<<endl;
          }
          break;
        case 'm':
          pConfig->bMigrationChangeEventDefined = true;
          if (bAcceptFullMigrMatrix){ // the -ema iTotalPops
            //                            //<matrix element list>
            if (arguments[iCurrentArg].size()<3){
              Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
                ", you need to at least specify the total number of populations.\n";
              Rcpp::stop("Argument error");
            }
            int iTotalPops = atoi(arguments[iCurrentArg][2].data());
            iTotalCells = iTotalPops * iTotalPops + 3;
            if (arguments[iCurrentArg].size()!=iTotalCells){
              Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
                ", the number of cells do not match the number of pops specified squared.\n";
              Rcpp::stop("Argument error");
            }
            //iType = Event::MIGRATION_MATRIX_RATE;
            iSubOption = 2;
            MatrixDouble dMigrationMatrix;
            for (int i=0;i<iTotalPops;++i){
              vector<double> newRow;
              for(int j=0; j<iTotalPops; ++j) {
                if (i==j){
                  newRow.push_back(0.0);
                  ++iSubOption;
                }
                else{
                  newRow.push_back(atof(arguments[iCurrentArg][++iSubOption].data()));
                }
              }
              dMigrationMatrix.push_back(newRow);
            }
            for(int i=0; i< iTotalPops; ++i) {
              for(int j=0; j<iTotalPops; ++j) {
                if (i!=j) dMigrationMatrix[i][i] +=
                  dMigrationMatrix[i][j];
              }
            }
            wrapper = EventPtr(new
                                 MigrationRateMatrixEvent(
                                   Event::MIGRATION_MATRIX_RATE,dTime,
                                   dMigrationMatrix));
            //Rcpp::Rcerr<<"Full migration matrix provided by the user\n";
            
          }else{
            // the -em t i j x option specify just
            //part of the migration matrix
            if (arguments[iCurrentArg].size()!=5){
              Rcpp::Rcerr<<"For flag "<<arguments[iCurrentArg][0]<<
                ", you must specify the source pop, dest pop, and the migration rate.\n";
              Rcpp::stop("Argument error");
            }else{
              //iType = Event::MIGRATION_RATE;
              //Rcpp::Rcerr<<"Mig rate of source pop "<<arguments[iCurrentArg][2] <<" to dest pop "<<arguments[iCurrentArg][3]<<" set to "<<arguments[iCurrentArg][4]<<".\n";
              wrapper = EventPtr(new MigrationRateEvent(
                Event::MIGRATION_RATE,dTime,
                atoi( arguments[iCurrentArg][2].data() ) -1,
                atoi( arguments[iCurrentArg][3].data() ) -1,atof( arguments[iCurrentArg][4].data() ) ));
            }
          }
          break;
        default:
          Rcpp::Rcerr<<"Invalid suboption, you entered"<<chType<<endl;
        break;
        }
        pEventList->push_back(wrapper);
        break;
      default:
        Rcpp::Rcerr<<"Invalid option, you entered "<<arguments[iCurrentArg][0][1]<<endl;
      }
    }catch(const out_of_range & e){
      Rcpp::Rcerr<<"There were too many arguments.\n";
    }
  }
  
  // Final sanity checks for the program before we begin:
  
  if (pConfig->iGeneConvTract>pConfig->dBasesToTrack) {
    Rcpp::Rcerr<<"Warning: the gene conversion tract (-c 2nd parameter) cannot be "<<
      "longer than the length of sequence (-h parameter) to retain. ";
    pConfig->dBasesToTrack=2.0*pConfig->iGeneConvTract;
    //Rcpp::Rcerr<<"The -h parameter is now revised to the recommend value of 2*tractlen = "
    //    <<pConfig->dBasesToTrack<<endl;
  }
  
  pEventList->sort(byEventTime());
  pConfig->pEventList = pEventList;
}


vector<AlphaSimRReturn> Simulator::beginSimulationMemory() {
  
  vector<AlphaSimRReturn> toRet;
  
  try {
    RandNumGenerator *rg = new RandNumGenerator(pConfig->iRandomSeed);
    for (unsigned int i = 0; i < pConfig->iIterations; ++i) {
      GraphBuilder graphBuilder = GraphBuilder(pConfig, rg);
      graphBuilder.build();
      vector<AlphaSimRReturn> tmp = graphBuilder.getMutations();
      graphBuilder.printHaplotypes();
      toRet.insert(toRet.end(), tmp.begin(), tmp.end());
    }
    delete rg;
  } catch (const char *message) {
    Rcpp::Rcerr << "Simulator caught exception with message:" << endl << message << endl;
  }
  return  toRet;
}


void Simulator::beginSimulation() {
  try {
    RandNumGenerator *rg = new RandNumGenerator(pConfig->iRandomSeed);
    for (unsigned int i = 0; i < pConfig->iIterations; ++i) {
      GraphBuilder graphBuilder = GraphBuilder(pConfig, rg);
      graphBuilder.build();
      graphBuilder.printHaplotypes();
      
    }
    delete rg;
  } catch (const char *message) {
    Rcpp::Rcerr << "Simulator caught exception with message:" << endl << message << endl;
  }
}

Simulator::Simulator() {
}


Simulator::~Simulator() {
  delete pConfig;
}

// AlphaSimR specific functions

vector<AlphaSimRReturn> runFromAlphaSimR(string in) {
  vector<std::string> words;
  Simulator simulator;
  
  if (in == ""){
    Rcpp::stop("Not enough args for macs call");
  }
  if (in.empty()) {
    Rcpp::stop("Not enough args for macs call");
  }
  boost::split(words, in, boost::is_any_of(", "), boost::token_compress_on);
  CommandArguments arguments;
  vector<string> subOption;
  // sample size
  subOption.emplace_back(words[0]);
  // seq length
  subOption.emplace_back(words[1]);
  arguments.push_back(subOption);
  subOption.clear();
  for (unsigned int i=2;i<words.size();++i){
    subOption.emplace_back(words[i]);
    if (i==words.size()-1 || (words[i+1][0]=='-' && words[i+1][1]>=65)){
      arguments.push_back(subOption);
      subOption.clear();
    }
  }
  if (arguments.size() == 0) {
    Rcpp::stop("Not enough args for macs call");
  }
  
  simulator.readInputParameters(arguments);
  vector<AlphaSimRReturn> test = simulator.beginSimulationMemory();
  return test;
}

// [[Rcpp::export]]
Rcpp::List MaCS(Rcpp::String args, arma::uvec maxSites, bool inbred, 
                arma::uword ploidy, int nThreads, Rcpp::StringVector seed){
  //Check input
  string t = args;
  if (t == "") {
    Rcpp::stop("error passing argument string");
  }
  
  // Output objects
  arma::uword nChr = maxSites.n_elem;
  arma::field<arma::Cube<unsigned char> > geno(nChr);
  arma::field<arma::vec > genMap(nChr);
  
  //Loop through chromosomes
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads)
#endif
  for(arma::uword chr=0; chr<nChr; chr++){
    // Run MaCS and check for valid output
    vector<AlphaSimRReturn> macsOutput;
    macsOutput = runFromAlphaSimR(args+seed(chr));
    
    arma::uword nSites, nBins, nHap, nInd;
    nSites = macsOutput.size();
    nHap = macsOutput[0].haplotypes.size();
    if(inbred){
      nInd = nHap;
    }else{
      nInd = nHap/ploidy;
    }
    arma::uvec selSites;
    if(maxSites(chr)>0){
      if(nSites<maxSites(chr)){
        maxSites(chr) = nSites;
      }
      selSites = sampleInt(maxSites(chr),nSites);
      nSites = maxSites(chr);
    }else{
      selSites.set_size(nSites);
      for(arma::uword i=0; i<nSites; ++i)
        selSites(i) = i;
    }
    
    // Fill genMap
    genMap(chr).set_size(nSites);
    for(arma::uword site=0; site<nSites; ++site){
      genMap(chr).at(site) = macsOutput[selSites(site)].length;
    }
    
    // Fill Geno
    nBins = nSites/8;
    if((nSites%8) > 0){
      ++nBins;
    }
    
    geno(chr).set_size(nBins,ploidy,nInd);
    arma::uword grp=0;
    arma::uword locus;
    std::bitset<8> workBits;
    for(arma::uword i=0; i<nHap; ++i){
      if(inbred){
        for(grp=0; grp<ploidy; ++grp){
          locus = 0;
          // Fill in bins known to be complete
          if(nBins > 1){
            for(arma::uword j=0; j<(nBins-1); ++j){
              for(arma::uword k=0; k<8; ++k){
                workBits[k] = macsOutput[selSites(locus)].haplotypes[i];
                ++locus;
              }
              geno(chr).slice(i).col(grp).row(j) = 
                toByte(workBits);
            }
          }
          // Fill in potentially incomplete bins
          for(arma::uword k=0; k<8; ++k){
            if(locus<nSites){
              workBits[k] = macsOutput[selSites(locus)].haplotypes[i];
              ++locus;
            }else{
              workBits[k] = 0;
            }
          }
          geno(chr).slice(i).col(grp).row(nBins-1) = 
            toByte(workBits);
        }
      }else{
        locus = 0;
        // Fill in bins known to be complete
        if(nBins > 1){
          for(arma::uword j=0; j<(nBins-1); ++j){
            for(arma::uword k=0; k<8; ++k){
              workBits[k] = macsOutput[selSites(locus)].haplotypes[i];
              ++locus;
            }
            geno(chr).slice(i/ploidy).col(grp).row(j) = 
              toByte(workBits);
          }
        }
        // Fill in potentially incomplete bins
        for(arma::uword k=0; k<8; ++k){
          if(locus<nSites){
            workBits[k] = macsOutput[selSites(locus)].haplotypes[i];
            ++locus;
          }else{
            workBits[k] = 0;
          }
        }
        geno(chr).slice(i/ploidy).col(grp).row(nBins-1) = 
          toByte(workBits);
        ++grp;
        grp = grp%ploidy; 
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("geno")=geno,
                            Rcpp::Named("genMap")=genMap);
}
