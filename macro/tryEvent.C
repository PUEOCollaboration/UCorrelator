#include "FFTtools.h" 

pueo::EventSummary * tryEvent(int run = 2, int event = -1 , bool nofilter = true, pueo::Dataset::DataDirectory datadir  = pueo::Dataset::PUEO_MC_DATA) 
{

//  AnalysisWaveform::enableDebug(true); 

  FFTtools::loadWisdom("wisdom.dat"); 

  gStyle->SetOptStat(0); 

  pueo::Dataset d(run, false, datadir); 
  if (event > 0)
  {
    d.getEvent(event); 
  }
  else
  {
    d.getEntry(-event); 
  }


  TFile out("test.root","RECREATE"); 
  pueo::FilterStrategy strategy(&out); 
  double fmins[]={0.23, 0.43}; 
  double fmaxs[]={0.29, 0.49}; 
  

  printf("strategy applied!\n"); 

//  AnalysisWaveform::enableDebug(true); 
  printf("creating event\n"); 
  pueo::FilteredEvent * fae = new pueo::FilteredEvent(d.useful(), &strategy, d.gps(), d.header(),true); 



  printf("processed strategy!\n"); 

  pueo::EventSummary * sum = new pueo::EventSummary; 
  pueo::UCorrelator::AnalysisConfig cfg; 
  cfg.use_forced_trigger_rms = false; 
  pueo::UCorrelator::Analyzer al(&cfg, true); ;
  al.analyze(fae, sum); 


  FFTtools::saveWisdom("wisdom.dat"); 
  return sum; 

}
