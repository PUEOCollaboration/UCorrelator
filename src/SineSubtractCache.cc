#include "pueo/SineSubtractCache.h"
#include <iostream>
#include "pueo/Dataset.h"
#include "TFile.h"
#include "TTree.h"
#include "pueo/FilteredEvent.h"
#include "pueo/FilterStrategy.h"
#include "pueo/RawHeader.h"
#include "TDirectory.h"


const TString ssrTreeName = "sineSubResultTree";

TString pueo::UCorrelator::SineSubtractCache::branchName(pol::pol_t pol, Int_t ant){
  return TString::Format("ssr_%d_%d", (int)pol, ant);
}

TString pueo::UCorrelator::SineSubtractCache::fileName(const char* specDir, UInt_t hash, Int_t run){
  int av = version::get();
  return TString::Format("%s/sineSubResults_%u_pueo%d_run%d.root", specDir, hash, av, run);
}


void pueo::UCorrelator::SineSubtractCache::makeCache(int run, SineSubtractFilter* ssf){

  pueo::UCorrelator::SineSubtractFilter::setUseCache(false);
  FFTtools::SineSubtractResult* results[pol::kNotAPol][k::NUM_HORNS] = {{NULL}};
  
  const char* ssDesc = ssf->description();
  if(!ssDesc){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", can't cache or use cached results without SineSubtract description!" << std::endl;
  }  
  else{
    UInt_t hash = TString(ssDesc).Hash();
    
    const char* specDir = getenv("UCORRELATOR_SPECAVG_DIR");
    if(!specDir){
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", can't cache or use cached results without UCORRELATOR_SPECAVG_DIR environment variable!" << std::endl;
    }
    else{

      std::cout << "Info in " << __PRETTY_FUNCTION__ << ", will generate cached sine subtraction results!" << std::endl;
      Dataset d(run);

      FilterStrategy fs;
      fs.addOperation(ssf);
  
      TFile* fOut = new TFile(fileName(specDir, hash, run), "recreate");
      TTree* tOut = new TTree(ssrTreeName, ssrTreeName);
      UInt_t eventNumber = 0;
      tOut->Branch("eventNumber", &eventNumber);
      for(int polInd=0; polInd < pol::kNotAPol; polInd++){
        pol::pol_t pol = (pol::pol_t) polInd;
        for(int ant=0; ant < k::NUM_HORNS; ant++){
          tOut->Branch(branchName(pol, ant), &results[pol][ant]);
        }
      }


      const int n = d.N();
      const double deltaPrint = double(n)/1000;
      double nextPrint = 0;
      for(int entry=0; entry < d.N(); entry++){
        d.getEntry(entry);

        FilteredEvent fEv(d.useful(), &fs, d.gps(), d.header());

        for(int polInd=0; polInd < pol::kNotAPol; polInd++){
          pol::pol_t pol = (pol::pol_t) polInd;          
          for(int ant=0; ant < k::NUM_HORNS; ant++){
            const FFTtools::SineSubtract* ss = ssf->sinsub(pol, ant);
            FFTtools::SineSubtractResult* ssr = const_cast<FFTtools::SineSubtractResult*>(ss->getResult());
            results[pol][ant] = ssr;
          }
        }

        eventNumber = d.header()->eventNumber;
        tOut->Fill();
        if(entry >= nextPrint){
          const int nm = 50;
          int m = nm*nextPrint/n;
          fprintf(stderr, "\r%4.2f %% complete", 100*nextPrint/n);
          std::cerr << "[";
          for(int i=0; i < m; i++) std::cerr << "=";
          for(int i=m; i < nm; i++) std::cerr << " ";
          std::cerr << "]";
          nextPrint += deltaPrint;
          if(nextPrint >= n){std::cerr << std::endl;}
        }
        // if(entry > 100) break;
      }
      tOut->BuildIndex("eventNumber"); // does this get saved?
      fOut->Write();
      fOut->Close();

    }
  }
}









/** 
 * Constructor
 * 
 * @param ssDesc description string, used for identifying file of cached results
 */
pueo::UCorrelator::SineSubtractCache::SineSubtractCache(const char* ssDesc)
  : fDebug(false), fFile(NULL), fTree(NULL), fDescHash(0), fSpecDir(), fCurrentRun(-1), fLastAttemptedRun(-1), fLastEventNumber(0)
{


  for(int pol=0; pol < pol::kNotAPol; pol++){
    for(int ant=0; ant < k::NUM_HORNS; ant++){
      results[pol][ant] = NULL;
    }
  }
  
  if(!ssDesc){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", can't cache or use cached results without SineSubtract description!" << std::endl;
  }  
  else{
    fDescHash = TString(ssDesc).Hash();
    
    const char* specDir = getenv("UCORRELATOR_SPECAVG_DIR");
    if(!specDir){
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", can't cache or use cached results without UCORRELATOR_SPECAVG_DIR environment variable!" << std::endl;
    }
    else{
      fSpecDir = specDir;
    }
  }
}



pueo::UCorrelator::SineSubtractCache::~SineSubtractCache(){

  if(fFile){
    fFile->Close();
  }
}



const FFTtools::SineSubtractResult* pueo::UCorrelator::SineSubtractCache::getResult(UInt_t eventNumber, pol::pol_t pol, Int_t ant){

  // hard to check whether anita version is correct...
  // this should happen

  // std::cerr << eventNumber << "\t" << pol << "\t" << ant << std::endl;

  if(eventNumber != fLastEventNumber){
    int run = Dataset::getRunContainingEventNumber(eventNumber);
    if(run!=fLastAttemptedRun){//fCurrentRun){
      loadRun(run);
    }
    if(fTree){
      Int_t entry = fTree->GetEntryNumberWithIndex(eventNumber);
      if(entry >= 0){
        fTree->GetEntry(entry);

        if(eventNumber != fLastEventNumber){
          std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", loaded sine subtract cache: run = " << fCurrentRun << ", eventNumber requested = "
                    << eventNumber << ", but eventNumber read = " << fLastEventNumber << std::endl;
        }
      }
      else{
        static int numEntryWarnings = 0;
        const  int maxEntryWarnings = 20;
        if(numEntryWarnings < maxEntryWarnings || fDebug){
          std::cerr << "Warning  " << (numEntryWarnings+1) << " of " << maxEntryWarnings << " in " << __PRETTY_FUNCTION__ << ", can't find entry "
                    << entry << " in " << fTree->GetName() << " in file " << fFile->GetName()
                    << " for eventNumber " << eventNumber << std::endl;
          numEntryWarnings++;
        }
        return NULL;
      }
    }
    else{
      static int numTreeWarnings = 0;
      const  int maxTreeWarnings = 20;
      if(numTreeWarnings < maxTreeWarnings || fDebug){
	std::cerr << "Warning  " << (numTreeWarnings+1) << " of " << maxTreeWarnings << " in " << __PRETTY_FUNCTION__ << ", can't find  "
		  << ssrTreeName << " in file " << fFile->GetName() << std::endl;
	numTreeWarnings++;
      }
    }
  }
  
  return results[pol][ant];
}



void pueo::UCorrelator::SineSubtractCache::loadRun(Int_t run){

  if(fSpecDir.Length() > 0 && fDescHash > 0){
    if(fFile){
      fFile->Close();
      fFile = NULL;
      fTree = NULL;
    }

    const TString theRootPwd = gDirectory->GetPath();
    TString fName = fileName(fSpecDir, fDescHash, run);
    fFile = TFile::Open(fName, "read");
    if(fFile){
      fTree = (TTree*) fFile->Get(ssrTreeName);

      fTree->SetBranchAddress("eventNumber", &fLastEventNumber);
      for(int polInd=0; polInd < pol::kNotAPol; polInd++){
        pol::pol_t pol = (pol::pol_t) polInd;
        for(int ant=0; ant < k::NUM_HORNS; ant++){
          fTree->SetBranchAddress(branchName(pol, ant), &results[pol][ant]);
        }
      }
      fTree->BuildIndex("eventNumber");
      fTree->GetEntry(0);
      fCurrentRun = run;

      if(fDebug){
        std::cerr << "Loaded first entry in run " << fCurrentRun << ", which has eventNumber " << fLastEventNumber << std::endl;
      }
      gDirectory->cd(theRootPwd); 
    }
    else{
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't open file " << fName << std::endl;
    }
    fLastAttemptedRun = run;
  }
}
