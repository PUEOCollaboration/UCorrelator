#include "pueo/TimeDependentAverage.h" 
#include "pueo/Version.h" 
#include "TFile.h" 
#include "pueo/UCImageTools.h" 
#include "TSystem.h"
#include "pueo/AnalysisWaveform.h" 
#include "pueo/FilterStrategy.h"
#include "TH2.h" 
#include "pueo/RawHeader.h" 
#include "pueo/Dataset.h" 
#include "TCut.h" 
#include "pueo/FilteredEvent.h" 
#include "TMutex.h" 
#include "pueo/BasicFilters.h"
#include <math.h>// for isnan



#define DEBUG_SPEC_AVG 


pueo::UCorrelator::TimeDependentAverage::~TimeDependentAverage()
{

  for (int i = 0; i < k::NUM_ANTS; i++) 
  {
    for (int j = 0; j < 2; j++) 
    {
      delete avgs[i][j]; 
      delete avgs_minbias[i][j]; 
      delete rms[i][j]; 
      if (peakiness[i][j]) delete peakiness[i][j]; 
      if (peakiness_minbias[i][j]) delete peakiness_minbias[i][j]; 

    }
  }

  delete nblasts; 
  delete norms; 
  delete norms_minbias; 
}

int pueo::UCorrelator::TimeDependentAverage::computeAverage(double max_r, int min_norm, double max_power) 
{
  Dataset d(run); 




  //figure out how long it is. 

  d.getEntry(0); 
  double startTime = d.header()->triggerTime; 
  d.getEntry(d.N()-1);
  double endTime = d.header()->triggerTime+1; 

  int nbins = (endTime-startTime) / nsecs; 

  TString name; 
  TString title; 
  name.Form("nblasts%d_%d",run,nsecs); 
  title.Form("N Blast Candidates run %d, nsec=%d", run,nsecs); 
  nblasts = new TH1I(name, title, nbins, startTime,endTime) ; 
  nblasts->SetDirectory(0); 
  nblasts->GetYaxis()->SetTitle("NBlasts"); 
  nblasts->GetXaxis()->SetTitle("Time"); 
  nblasts->GetXaxis()->SetTimeDisplay(true); 
  nblasts->GetXaxis()->SetTimeOffset(0); 

  name.Form("norm_%d_%d",run,nsecs); 
  title.Form("RF Normalization run %d, nsec=%d", run,nsecs); 
  
  norms = new TH1I(name, title,nbins,startTime,endTime); 
  norms->SetDirectory(0); 
  norms->GetYaxis()->SetTitle("NEvents"); 
  norms->GetXaxis()->SetTitle("Time"); 
  norms->GetXaxis()->SetTimeDisplay(true); 
  norms->GetXaxis()->SetTimeOffset(0); 

  name.Form("norm_minbias_%d_%d",run,nsecs); 
  title.Form("MinBias Normalization run %d, nsec=%d", run,nsecs); 
  
  norms_minbias = new TH1I(name, title,nbins,startTime,endTime); 
  norms_minbias->SetDirectory(0); 
  norms_minbias->GetYaxis()->SetTitle("NEvents"); 
  norms_minbias->GetXaxis()->SetTitle("Time"); 
  norms_minbias->GetXaxis()->SetTimeDisplay(true); 
  norms_minbias->GetXaxis()->SetTimeOffset(0); 


//  printf("%d\n", nbins); 

  for (int ant = 0; ant < k::NUM_HORNS; ant++)
  {
    for (int pol = 0; pol < 2; pol++) 
    {

      name.Form("specavg_%d_%d_%d_%d",run,nsecs,ant,pol); 
      title.Form("RF Spectrum Average ant=%d, pol=%d, run %d, nsec=%d", ant,pol,run,nsecs); 

        
      avgs[ant][pol] = new TH2F(name, title, 
                        131,0,1.31,//TODO don't hardcode here
                        nbins, startTime,endTime) ; 
      avgs[ant][pol]->SetDirectory(0); 
      avgs[ant][pol]->GetYaxis()->SetTitle("Time"); 
      avgs[ant][pol]->GetYaxis()->SetTimeDisplay(true); 
      avgs[ant][pol]->GetYaxis()->SetTimeOffset(0); 
      avgs[ant][pol]->GetXaxis()->SetTitle("Frequency"); 


      name.Form("specavg_minbias_%d_%d_%d_%d",run,nsecs,ant,pol); 
      title.Form("MinBias Spectrum Average ant=%d, pol=%d, run %d, nsec=%d", ant,pol,run,nsecs); 

      avgs_minbias[ant][pol] = new TH2F(name, title, 
                        131,0,1.31,//TODO don't hardcode here
                        nbins, startTime,endTime) ; 
      avgs_minbias[ant][pol]->SetDirectory(0); 
      avgs_minbias[ant][pol]->GetYaxis()->SetTitle("Time"); 
      avgs_minbias[ant][pol]->GetYaxis()->SetTimeDisplay(true); 
      avgs_minbias[ant][pol]->GetYaxis()->SetTimeOffset(0); 
      avgs_minbias[ant][pol]->GetXaxis()->SetTitle("Frequency"); 

      name.Form("rms_%d_%d_%d_%d",run,nsecs,ant,pol); 
      title.Form("MinBias RMS ant=%d, pol=%d, run %d, nsec=%d", ant,pol,run,nsecs); 
      rms[ant][pol] = new TH1D(name, title, nbins, startTime,endTime) ; 
      rms[ant][pol]->SetDirectory(0); 
      rms[ant][pol]->GetYaxis()->SetTitle("RMS (mVish)"); 
      rms[ant][pol]->GetXaxis()->SetTitle("Time"); 
      rms[ant][pol]->GetXaxis()->SetTimeDisplay(true); 
      rms[ant][pol]->GetXaxis()->SetTimeOffset(0); 

 


    }
  }



                       
  FilterStrategy str; 


  int N = d.N(); 

   N = d.N(); //For all the events.
  for (int i = 0; i < N; i++)
  {
    d.getEntry(i); 
    FilteredEvent ev(d.useful(), &str, d.gps(), d.header()); 
    
    double t= d.header()->triggerTime; 

    //cut out possible blasts. This is perhaps a bit aggressive, but it's worth it 
    double max_ratio_hpol, max_ratio_vpol; 
    ev.getMinMaxRatio(pol::kHorizontal, &max_ratio_hpol,0,0,0);
    if (max_ratio_hpol > max_r || ev.getAveragePower(pol::kHorizontal)  > max_power) 
    {
      nblasts->Fill(t); 
      continue; 
    }
    ev.getMinMaxRatio(pol::kVertical, &max_ratio_vpol,0,0,0); 
    if (max_ratio_vpol > max_r || ev.getAveragePower(pol::kVertical) > max_power)
    {
      nblasts->Fill(t); 
      continue; 
    }

    bool isRF = d.header()->trigType & 1; 

    (isRF ? norms : norms_minbias)->Fill(t); 
    int tbin = norms->FindFixBin(t); 

#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for 
#endif
    for (int j = 0; j < k::NUM_ANTS * 2; j++) 
    {
      int ant = j /2; 
      int pol = j %2; 
      const AnalysisWaveform * wf = ev.getFilteredGraph(ant,pol::pol_t(pol)); 


      if (!isRF) 
      {
        rms[ant][pol]->SetBinContent(tbin, rms[ant][pol]->GetBinContent(tbin) + wf->even()->GetRMS(2)); 
      }

      const TGraphAligned * g = wf->power(); 
      for (int k = 0; k<avgs[ant][pol]->GetNbinsX() ; k++)
      {
        TH2F * h = isRF ? avgs[ant][pol] : avgs_minbias[ant][pol]; 
        h->SetBinContent(k,tbin,h->GetBinContent(k,tbin) + g->GetY()[k]); 
      }
    }

    if (i % 50 == 0)
    {
      printf("\r%d/%d",i,N);
    }
    else {
      printf(".");
      fflush(stdout);
    }
    
  }
  printf("..Done\n"); 


  for (int isRF = 0; isRF <=1; isRF++)
  {
    for (int j = 0; j < k::NUM_ANTS * 2; j++) 
    {
      int ant = j /2; 
      int pol = j %2; 

      TH2F * h = isRF ? avgs[ant][pol] : avgs_minbias[ant][pol]; 
      TH1I * hnorm = isRF ? norms : norms_minbias; 


      for (int ybin = 1; ybin <= h->GetNbinsY(); ybin++)
      {
        if (hnorm->GetBinContent(ybin) < min_norm) 
        {
          continue; 
        }

        for (int k = 0; k<h->GetNbinsX() ; k++)
        {
            h->SetBinContent(k,ybin, h->GetBinContent(k,ybin) / hnorm->GetBinContent(ybin)); 
        }

        if (!isRF) 
        {
            rms[ant][pol]->SetBinContent(ybin, rms[ant][pol]->GetBinContent(ybin) / hnorm->GetBinContent(ybin)); 
        }

      }

      //now go through and fill in the rows with no content with an average of ones that do 
      for (int ybin = 1; ybin <= h->GetNbinsY(); ybin++)
      {

        if (hnorm->GetBinContent(ybin) < min_norm)
        {
          int last_bin = ybin-1; 
          while(hnorm->GetBinContent(last_bin) < min_norm && last_bin > 0) last_bin--; 

          int next_bin = ybin+1; 
          while (!hnorm->GetBinContent(next_bin) && next_bin <= h->GetNbinsY()) next_bin++; 

          //lost cause, no good bins
          if (last_bin == 0 && next_bin > h->GetNbinsY()) break; 

          if (last_bin == 0) 
          {
            for (int k = 0; k<h->GetNbinsX() ; k++)
            {
                h->SetBinContent(k,ybin, h->GetBinContent(k,next_bin)); 
                if (!isRF) rms[ant][pol]->SetBinContent(ybin, rms[ant][pol]->GetBinContent(next_bin)); 
            }
          }
          else if (next_bin > h->GetNbinsY()) 
          {
            for (int k = 0; k<h->GetNbinsX() ; k++)
            {
                h->SetBinContent(k,ybin, h->GetBinContent(k,last_bin)); 
                if (!isRF) rms[ant][pol]->SetBinContent(ybin, rms[ant][pol]->GetBinContent(last_bin)); 
            }
          }
          else
          {
            double last_frac = double(ybin - last_bin) / (next_bin-last_bin); 
            
            for (int k = 0; k<h->GetNbinsX() ; k++)
            {
                h->SetBinContent(k,ybin, last_frac * h->GetBinContent(k,last_bin)  + (1-last_frac) * h->GetBinContent(k,next_bin)); 
                if (!isRF) rms[ant][pol]->SetBinContent(ybin, last_frac * rms[ant][pol]->GetBinContent(last_bin) + (1-last_frac) * rms[ant][pol]->GetBinContent(next_bin)); 
            }
          }
        }
      }
    }

  }

  return 0; 
}



void pueo::UCorrelator::TimeDependentAverage::saveToDir(const char * dir) 
{

  gSystem->mkdir(dir,true); 

  TString fstr; 
  fstr.Form("%s/%d_%d.root", dir, run,nsecs); 
  TFile f(fstr.Data(),"RECREATE"); 
  f.cd(); 

  //make sure averages and RMS are loaded
  //getSpectrogram(pol::kHorizontal,0); 
  //getRMS(pol::kHorizontal,0); 
  
  f.mkdir("specavg"); 
  f.cd("specavg"); 

  for (int ant = 0; ant < k::NUM_ANTS; ant++) 
  {
     avgs[ant][0]->Write(TString::Format("h%d",ant)); 
     avgs[ant][1]->Write(TString::Format("v%d",ant)); 
  }

  f.cd(); 

  f.mkdir("specavg_minbias"); 
  f.cd("specavg_minbias"); 
  for (int ant = 0; ant < k::NUM_ANTS; ant++) 
  {
     avgs_minbias[ant][0]->Write(TString::Format("h%d",ant)); 
     avgs_minbias[ant][1]->Write(TString::Format("v%d",ant)); 
  }
  f.cd();

  f.mkdir("rms"); 
  f.cd("rms"); 
  for (int ant = 0; ant < k::NUM_ANTS; ant++) 
  {
     rms[ant][0]->Write(TString::Format("h%d",ant)); 
     rms[ant][1]->Write(TString::Format("v%d",ant)); 
  }
  f.cd();

  nblasts->Write("nblasts"); 
  norms->Write("norms"); 
  norms_minbias->Write("norms_minbias"); 


}

const TH2F * pueo::UCorrelator::TimeDependentAverage::getSpectrogram(pol::pol_t pol, int ant, bool minbias) const
{

  
  __sync_synchronize(); //memory barrier

  if (!avgs_loaded)
  {
    m.Lock(); 
    if (!avgs_loaded) 
    {
      TFile * f = TFile::Open(fname); 
      f->cd("specavg"); 
      for (int ant = 0; ant < k::NUM_ANTS; ant++) 
      {
        TH2F * found_hpol = (TH2F*) gDirectory->Get(TString::Format("h%d",ant)); 
        avgs[ant][0] = new TH2F(*found_hpol);  
        avgs[ant][0]->SetDirectory(0); 

        TH2F * found_vpol = (TH2F*) gDirectory->Get(TString::Format("v%d",ant)); 
        avgs[ant][1] = new TH2F(*found_vpol);  
        avgs[ant][1]->SetDirectory(0); 
      }

      f->cd("specavg_minbias"); 

      for (int ant = 0; ant < k::NUM_ANTS; ant++) 
      {

        TH2F * found_hpol = (TH2F*) gDirectory->Get(TString::Format("h%d",ant)); 
        avgs_minbias[ant][0] = new TH2F(*found_hpol);  
        avgs_minbias[ant][0]->SetDirectory(0); 

        TH2F * found_vpol = (TH2F*) gDirectory->Get(TString::Format("v%d",ant)); 
        avgs_minbias[ant][1] = new TH2F(*found_vpol);  
        avgs_minbias[ant][1]->SetDirectory(0); 
      }
      delete f ; 
      __sync_synchronize(); //memory barrier
      avgs_loaded = true; 
    }
    m.UnLock(); 
  }


  return minbias ? avgs_minbias[ant][pol] : avgs[ant][pol]; 
}


const TH1D * pueo::UCorrelator::TimeDependentAverage::getRMS(pol::pol_t pol, int ant) const
{
  __sync_synchronize(); //memory barrier
  if (!rms_loaded)
  {
    m.Lock(); 

    if (!rms_loaded)
    {
      TFile *f = TFile::Open(fname); 
      f->cd("rms"); 

      for (int ant = 0; ant < k::NUM_ANTS; ant++) 
      {

        TH1D * found_hpol = (TH1D*) gDirectory->Get(TString::Format("h%d",ant)); 
        rms[ant][0] = new TH1D(*found_hpol);  
        rms[ant][0]->SetDirectory(0); 

        TH1D * found_vpol = (TH1D*) gDirectory->Get(TString::Format("v%d",ant)); 
        rms[ant][1] = new TH1D(*found_vpol);  
        rms[ant][1]->SetDirectory(0); 
      }
      __sync_synchronize(); //memory barrier
      rms_loaded = true; 
      delete f; 
    }
    m.UnLock(); 
  }

  return rms[ant][pol]; 

}

pueo::UCorrelator::TimeDependentAverage::TimeDependentAverage(int run, int nsecs, const char * persistdir,
      double max_bottom_top_ratio, int min_norm, double max_power) 
  :  nsecs(nsecs) , run(run)
{
  bool foundit = false;
  rms_loaded = false;
  peakiness_loaded = false;
  avgs_loaded = false; 

  const char * check_dir = persistdir ?: getenv("UCORRELATOR_TIMEAVG_DIR") ?: getenv("UCORRELATOR_SPECAVG_DIR") ?: 0; 


  if (check_dir)
  {
    fname.Form("%s/%d_%d.root", check_dir, run,nsecs); 
    TFile * f = TFile::Open(fname); 
    if (f && f->IsOpen())
    {
      f->cd(); 

      if (gDirectory->Get("norms") && gDirectory->Get("norms_minbias") && gDirectory->Get("nblasts"))
      {
        TH1I * found= (TH1I*) gDirectory->Get("norms"); 
        norms = new TH1I(*found); 
        norms->SetDirectory(0); 

        found= (TH1I*) gDirectory->Get("norms_minbias"); 
        norms_minbias = new TH1I(*found); 
        norms_minbias->SetDirectory(0); 

        found= (TH1I*) gDirectory->Get("nblasts"); 
        nblasts = new TH1I(*found); 
        nblasts->SetDirectory(0); 

        foundit = true;
      }
      else
      {
        fprintf(stderr,"Could not find one or more of norms, norms_minbias or nblasts inside %s/%d_%d.root. Are your spectrum averages out of date? ", check_dir,run,nsecs); 

      }
    }
    else
    {
      printf("Could not open %s/%d_%d.root", check_dir, run, nsecs); 
    }
    memset(rms,0,sizeof(rms)); 

    memset(avgs,0,sizeof(avgs)); 
    memset(avgs_minbias,0,sizeof(avgs_minbias)); 

    if (f) delete f; 
  }

  if (!foundit) 
  {
    printf("Didn't find spectrum_averages... creating!\n"); 
    if (!check_dir) 
    {
      printf("Define UCORRELATOR_TIMEAVG_DIR to persist somewhere.\n"); 
    }
    computeAverage(  max_bottom_top_ratio, min_norm, max_power); 
    if (check_dir) saveToDir(check_dir); 
  }

  memset(peakiness,0,sizeof(peakiness)); 
  memset(peakiness_minbias,0,sizeof(peakiness_minbias)); 

}
            


TH1* pueo::UCorrelator::TimeDependentAverage::getSpectrumAverage(pol::pol_t pol, int ant, double t, bool db, bool minbias)  const
{

  const TH2 * h =getSpectrogram(pol,ant,minbias); 
  //figure out which bin we are in 
   int bin =  h->GetYaxis()->FindFixBin(t); 

   if (bin == 0 || bin == h->GetNbinsY()+1) return 0 ; 

   TString title; 
   title.Form("Ant%d %c-Pol Average at t=%g", ant, pol == pol::kHorizontal ? 'H' : 'V', t); 
   TString name;
   name.Form("specavg_%d_%d_%g",ant,pol,t); 
   TH1D * answer =  new TH1D(name,title, 131,0,1.31); //todo don't hardcode 

   for (int i = 1; i <= answer->GetNbinsX(); i++)
   {
     answer->SetBinContent(i, h->GetBinContent(i,bin)); 
   }

   answer->GetXaxis()->SetTitle("Frequency"); 
   answer->GetYaxis()->SetTitle(db ? "Power (dBish)" :"Power (linear)"); 


   if (db) 
   {
     for (int i = 1; i <= answer->GetNbinsX(); i++) 
     {
       answer->SetBinContent(i, 10 * TMath::Log10(answer->GetBinContent(i))); 
     }
   }


   return answer;
  

}

TH1 *pueo::UCorrelator::TimeDependentAverage::getSpectrumPercentile(pol::pol_t pol, int ant, double pct , bool db, bool minbias ) const
{

  TH1 * answer = pueo::UCorrelator::image::getPctileProjection( getSpectrogram(pol,ant,minbias), 1, pct, true, getNorms(minbias)); 

  answer->GetXaxis()->SetTitle("Frequency"); 
  answer->GetYaxis()->SetTitle(db ? "Pctile Power (dBish)" :"Pctile Power (linear)"); 

   answer->SetTitle(TString::Format("Ant%d %c-Pol %gth Pctile", ant, pol == pol::kHorizontal ? 'H' : 'V', pct*100)); 

   if (db) 
   {
     for (int i = 1; i <= answer->GetNbinsX(); i++) 
     {
       answer->SetBinContent(i, 10 * TMath::Log10(answer->GetBinContent(i))); 
     }
   }


   return answer;
 
}


double pueo::UCorrelator::TimeDependentAverage::getStartTime() const 
{

  return norms ? norms->GetXaxis()->GetXmin() : -1; 
}

double pueo::UCorrelator::TimeDependentAverage::getEndTime() const 
{
  return norms ? norms->GetXaxis()->GetXmax() : -1;   
}


// to take advantage of magic statics
struct ThermalAverageWrapper
{
  pueo::UCorrelator::TimeDependentAverage * avg; 
  ThermalAverageWrapper()
  {
     if (pueo::version::get() > 1) 
     {
       fprintf(stderr,"warning: using default terminated thermal spectrum for P1 for peakiness\n"); 
     }

     TString dir; 
     dir.Form("%s/share/UCorrelator/terminated_noise/", getenv("PUEO_UTIL_INSTALL_DIR")); 

     avg = new pueo::UCorrelator::TimeDependentAverage(11382,60, dir.Data()); 
  }

  static ThermalAverageWrapper & get() 
  {
    static ThermalAverageWrapper t;
    return t; 
  }

};

const pueo::UCorrelator::TimeDependentAverage* pueo::UCorrelator::TimeDependentAverage::defaultThermal() 
{
  return ThermalAverageWrapper::get().avg;
}



void pueo::UCorrelator::TimeDependentAverage::computePeakiness(const TimeDependentAverage * thermalSpec, double fractionForNormalization) const
{

  if (!thermalSpec) thermalSpec = defaultThermal(); 


//#ifdef UCORRELATOR_OPENMP
//#pragma omp parallel for
//#endif 
  for (int ant  = 0; ant < k::NUM_ANTS; ant++)
  {
//    printf("%d\n",ant);
    for (int ipol = 0; ipol < 2; ipol++) 
    {

      for (int minbias = 0; minbias < 2; minbias++)
      {
//      printf("%d %d\n",ant,ipol);

        //get median spectrum 
        TH1 * median = getSpectrumPercentile(pol::pol_t(ipol), ant, 0.5,false,minbias); 
        TH1 * thermal = thermalSpec->getSpectrumPercentile(pol::pol_t(ipol),ant,0.5,false,minbias);

        // now we have to figure out scale. To do this, we compute the mean of the middle frac of points. 
        int index_spec[median->GetNbinsX()]; 
        int index_therm[thermal->GetNbinsX()]; 

        TMath::Sort(thermal->GetNbinsX(),((TH1F*) thermal)->GetArray()+1, index_therm); 
        TMath::Sort(median->GetNbinsX(), ((TH1F*)median)->GetArray()+1, index_spec); 

        double sum_spec = 0; 
        double sum_therm =0;  

        for (int i = 0; i <= thermal->GetNbinsX()*fractionForNormalization; i++) 
        {
          sum_therm+=thermal->GetBinContent(index_therm[ (int) (i + (thermal->GetNbinsX() * (0.5-fractionForNormalization)))  ]); 
        }

        for (int i = 0; i <= median->GetNbinsX()*fractionForNormalization; i++) 
        {
          sum_spec+=median->GetBinContent(index_spec[ (int) (i + median->GetNbinsX() * (0.5-fractionForNormalization))]); 
        }

        sum_therm /= thermal->GetNbinsX(); 
        sum_spec /= median->GetNbinsX(); 
        double ratio = sum_spec / sum_therm; 

        TString name; 
        name.Form("peakiness_%d_%d_%s\n", ant,ipol, minbias ? "minbias" : "rf"); 
        TString title; 
        title.Form("%s peakiness ant=%d pol=%d\n", minbias ? "minbias" : "RF" , ant,ipol); 

        const TH2 * avg = getSpectrogram(pol::pol_t(ipol),ant,minbias); 
        TH2D * peaky = new TH2D(name,title,
                                         avg->GetNbinsX(), avg->GetXaxis()->GetXmin(), avg->GetXaxis()->GetXmax(), 
                                         avg->GetNbinsY(), avg->GetYaxis()->GetXmin(), avg->GetYaxis()->GetXmax());  


        peaky->SetDirectory(0); 

        for (int ii = 1; ii < peaky->GetNbinsX(); ii++)
        {
          for (int jj = 1; jj < peaky->GetNbinsY(); jj++)
          {
            peaky->SetBinContent(ii,jj, avg->GetBinContent(ii,jj)  / ( ratio * thermal->GetBinContent(ii))); 
            if (std::isnan(peaky->GetBinContent(ii,jj)))
                peaky->SetBinContent(ii,jj,0); 
          }
        }
         (minbias ? peakiness_minbias[ant][ipol] : peakiness[ant][ipol]) = peaky; 
        delete median; 
        delete thermal; 
      }

    }
  }

}

pueo::UCorrelator:: TimeDependentAverageLoader::TimeDependentAverageLoader(const char * the_dir, int secs) 
  : dir(the_dir), nsecs(secs)
{
  tavg = 0; 
}


static TMutex mut; 
const pueo::UCorrelator::TimeDependentAverage* pueo::UCorrelator::TimeDependentAverageLoader::avg(double t) const
{
  TLockGuard l(&mut); 

//  printf("%g\n",t); 

  if (tavg && t >= tavg->getStartTime()-5 && t <= tavg->getEndTime()+5 ) 
  {
    return tavg; 
  }

  int run = Dataset::getRunAtTime(t); 
//  printf("loading average from run %d\n",run); 
  
  if (tavg) delete tavg; 
  tavg = new TimeDependentAverage(run,nsecs, dir); 
  
  return tavg; 

}

double pueo::UCorrelator::TimeDependentAverage::getBlastFraction(double t) const
{
  //Estimate by estimating nblasts and n, taking their ratio
  int bin = nblasts->GetXaxis()->FindFixBin(t); 

  //the other closest bin
  int other_bin = t < nblasts->GetXaxis()->GetBinCenter(bin) ? bin-1 : bin+1; 


  //the fraction of blasts in the closest bin
  double blast_frac_bin = nblasts->GetBinContent(bin) / (nblasts->GetBinContent(bin) + norms->GetBinContent(bin)); 

  //f is the weight of the main bin. 
  double f = norms->GetBinContent(other_bin) == 0 ? 1 :  1-fabs(t-nblasts->GetXaxis()->GetBinCenter(bin)) / nsecs; 
  if (f < 1 )
  {
    double blast_frac_other_bin = nblasts->GetBinContent(other_bin) / (nblasts->GetBinContent(other_bin) + norms->GetBinContent(other_bin)); 
    return f * blast_frac_bin + (1-f) * blast_frac_other_bin; 
  }

  return blast_frac_bin; 
}


TMutex loader_map_lock; 
std::map<int,pueo::UCorrelator::TimeDependentAverageLoader*> loader_map; 

static pueo::UCorrelator::TimeDependentAverageLoader * getLoader(int nsecs) 
{
  TLockGuard l(&loader_map_lock); 

  if (loader_map.count(nsecs))
  {
    return loader_map[nsecs]; 
  }

  pueo::UCorrelator::TimeDependentAverageLoader *  ldr = new pueo::UCorrelator::TimeDependentAverageLoader(0,nsecs); 
  loader_map[nsecs] = ldr; 
  return ldr; 
}

double pueo::UCorrelator::TimeDependentAverage::getRMS(pol::pol_t pol, int ant, double t) const
{
  return ((TH1*) getRMS(pol,ant))->Interpolate(t); 
}



double pueo::UCorrelator::TimeDependentAverageLoader::getRMS(double t, int ipol, int ant, int nsecs) 
{
  return getLoader(nsecs)->avg(t)->getRMS(pol::pol_t(ipol),ant,t); 
}

double pueo::UCorrelator::TimeDependentAverageLoader::getPayloadBlastFraction(double t,  int nsecs) 
{
  return getLoader(nsecs)->avg(t)->getBlastFraction(t); 
}


const TH2D * pueo::UCorrelator::TimeDependentAverage::getPeakiness(pol::pol_t pol, int ant, bool minbias) const
{

  //make sure we have averages loaded 
  getSpectrogram(pol::kHorizontal,0); 

  __sync_synchronize(); //memory barrier
  if (!peakiness_loaded)
  {
    m.Lock();
    if (!peakiness_loaded)
    {
      computePeakiness(); 
     __sync_synchronize(); //memory barrier
      peakiness_loaded = true; 
    }
    m.UnLock(); 
  }

  return minbias ? peakiness_minbias[ant][pol] : peakiness[ant][pol]; 
}
