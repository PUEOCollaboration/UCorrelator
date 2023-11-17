#include "pueo/Baseline.h" 
#include "TGraph.h" 
#include <cstdio>
#include "pueo/FilteredEvent.h" 
#include "TSystem.h" 
#include "pueo/FilterStrategy.h" 
#include "pueo/Dataset.h"
#include "pueo/UsefulEvent.h"
#include "FFTtools.h"
#include "pueo/RawHeader.h"
#include "pueo/UCFlags.h"
#include "TFile.h" 
#include "FFTWComplex.h" 

#define NSAMPLES 4096  
#define NOMINAL_DT 1./3


static int makeBaselines(int run, TGraph ** hpol, TGraph ** vpol, int N = 5000) 
{
  char *  datadir = getenv("PUEO_ROOT_DATA"); 
  if (!datadir) 
  {
    fprintf(stderr,"The environmental variable PUEO_ROOT_DATA must be defined! Aborting."); 
    exit(1); 
  }


  pueo::Dataset d(run); 
  pueo::FilterStrategy empty; 

  int nevents = 0; 
  int i = 10; 


  while (nevents < N && i < d.N()) 
  {
    d.getEntry(i++); 
 
    // skip non-RF triggers and a small fraction of times? is this really right? This is what Abby has though
    if (d.header()-> trigType != 1 && d.header()->triggerTimeNs <= 1e6)  
       continue; 

    //skip saturated events 
    
    pueo::FilteredEvent fae(d.useful(), &empty, d.gps(), d.header()); 

    if (fae.checkSaturation()) 
      continue; 

    for (int ant = 0; ant < pueo::k::NUM_HORNS; ant++) 
    {
       if (nevents == 0)
       {
           hpol[ant] = new TGraph(NSAMPLES/2+1); 
           vpol[ant] = new TGraph(NSAMPLES/2+1); 
       }
       
       TGraph * gh = d.useful()->makeGraph(ant, pueo::pol::kHorizontal); 
       TGraph * gv = d.useful()->makeGraph(ant, pueo::pol::kVertical); 

       TGraph * igh = FFTtools::getInterpolatedGraph(gh, NOMINAL_DT); 
       TGraph * igv = FFTtools::getInterpolatedGraph(gv, NOMINAL_DT); 

       igh->Set(NSAMPLES); 
       igv->Set(NSAMPLES); 

       FFTWComplex * fft_h = FFTtools::doFFT(igh->GetN(), igh->GetY()); 
       FFTWComplex * fft_v = FFTtools::doFFT(igv->GetN(), igv->GetY()); 

       for (int j = 0; j < NSAMPLES/2+1; j++) 
       {
         if (nevents == 0) 
         {
           hpol[ant]->GetX()[j] = j  / (NOMINAL_DT * NSAMPLES);  //GHz
           vpol[ant]->GetX()[j] = j  / (NOMINAL_DT * NSAMPLES);  //GHz
         }

         hpol[ant]->GetY()[j] += fft_h[j].getAbs() / (NSAMPLES/2.+1)  * (j == 0 || j == NSAMPLES/2 ? 1 : 2); 
         vpol[ant]->GetY()[j] += fft_v[j].getAbs() / (NSAMPLES/2.+1)  * (j == 0 || j == NSAMPLES/2 ? 1 : 2); 

       }

       delete [] fft_h; 
       delete [] fft_v; 
       delete gh; 
       delete gv; 
       delete igh; 
       delete igv; 

    }

    nevents++; 
  }

  for (int ant = 0; ant < pueo::k::NUM_HORNS; ant++)
  {
    for (int i = 0; i < NSAMPLES / 2 + 1; i++) 
    {
      hpol[ant]->GetY()[i] /= nevents; 
      vpol[ant]->GetY()[i] /= nevents; 
    }
  }


  return nevents; 
}


void pueo::UCorrelator::Baseline::saveToDir(const char * dir) 
{

  gSystem->mkdir(dir,true); 
  TFile f(TString::Format("%s/%d_%d.root", dir, run,navg),"RECREATE"); 
  for (int ant = 0; ant < k::NUM_HORNS; ant++) 
  {
    hpol[ant]->Write(TString::Format("h%d",ant)); 
    vpol[ant]->Write(TString::Format("v%d",ant)); 
  }
}

pueo::UCorrelator::Baseline::Baseline(int run, int navg, const char * persistdir) 
  : navg(navg) , run(run)
{
  bool foundit = false;
  if (persistdir)  
  {
    TFile f(TString::Format("%s/%d_%d.root", persistdir, run,navg)); 
    if (f.IsOpen())
    {
      for (int ant = 0; ant < k::NUM_HORNS; ant++) 
      {
        TGraph * found_hpol = (TGraph*) f.Get(TString::Format("h%d",ant)); 
        hpol[ant] = new TGraph(*found_hpol);  

        TGraph * found_vpol = (TGraph*) f.Get(TString::Format("v%d",ant)); 
        vpol[ant] = new TGraph(*found_vpol);  
      }
      foundit = true;
    }
  }

  if (!foundit) 
  {
    printf("Didn't find baselines... creating!\n"); 
    if (!persistdir) 
    {
      printf("Define UCORRELATOR_BASELINE_DIR to persist somewhere.\n"); 
    }
    makeBaselines(run, hpol, vpol, navg); 
  }


  hpol_avg = new TGraph(NSAMPLES/2+1); 
  vpol_avg = new TGraph(NSAMPLES/2+1); 

  memcpy(hpol_avg->GetX(), hpol[0]->GetX(), sizeof(double) * hpol_avg->GetN());
  memcpy(vpol_avg->GetX(), vpol[0]->GetX(), sizeof(double) * vpol_avg->GetN());

  for (int ant = 0; ant < k::NUM_HORNS; ant++) 
  {
    for (int i =0; i < NSAMPLES/2 +1; i++) 
    {
      hpol_avg->GetY()[i] += hpol[ant]->GetY()[i] / k::NUM_HORNS; 
      vpol_avg->GetY()[i] += vpol[ant]->GetY()[i] / k::NUM_HORNS; 
    }
  }

  if (persistdir) saveToDir(persistdir); 
}
            
