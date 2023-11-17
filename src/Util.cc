#include "pueo/UCUtil.h" 
#include "pueo/UsefulAttitude.h" 
#include "pueo/RawHeader.h" 
#include "pueo/AnalysisConfig.h" 

/* nothing to see here, please move along */ 


static const pueo::UCorrelator::AnalysisConfig defaultConfig;

static const double  C_in_m_ns = 0.299792;


double pueo::UCorrelator::getWAISDt(const UsefulAttitude * pat, const RawHeader * hdr, pol::pol_t pol, const AnalysisConfig * cfg, double * distance)
{

  if (!cfg) cfg = &defaultConfig; 

  double time_to_wais = ((UsefulAttitude*) pat)->getWaisDivideTriggerTimeNs(); 
  if (distance) *distance = time_to_wais * C_in_m_ns; 
  unsigned trig_time = hdr->triggerTimeNs; 
  return trig_time + (pol == pol::kHorizontal ? cfg->wais_hpol : cfg->wais_vpol).GPS_offset- time_to_wais; 
}


bool pueo::UCorrelator::isWAISVPol(const UsefulAttitude * pat, const RawHeader * hdr, const AnalysisConfig * cfg )  
{

  if (!cfg) cfg = &defaultConfig; 

  double distance; 
  double dt = getWAISDt(pat,hdr, pol::kVertical, cfg,&distance); 
  if (distance > cfg->wais_vpol.max_distance) return false; 
  if (fabs(dt) > cfg->wais_vpol.max_dt) return false; 

  return true; 
}

bool pueo::UCorrelator::isWAISHPol(const UsefulAttitude * pat, const RawHeader * hdr, const AnalysisConfig * cfg )  
{
  if (!cfg) cfg = &defaultConfig; 

  double distance; 
  double dt = getWAISDt(pat,hdr, pol::kHorizontal, cfg,&distance); 

//  printf("dt: %f; distance: %f\n", dt, distance); 
  if (distance > cfg->wais_hpol.max_distance) return false; 
  if (fabs(dt) > cfg->wais_hpol.max_dt) return false; 

  return true; 
}




bool pueo::UCorrelator::isLDB(const RawHeader * hdr, const AnalysisConfig * cfg )  
{
  if (!cfg) cfg = &defaultConfig; 

  if (hdr->run > cfg->ldb_max_run) return false; 
  if (! cfg->ldb_hist()) return false; 

  int binx = cfg->ldb_hist()->GetXaxis()->FindFixBin(hdr->triggerTime); 
  int biny = cfg->ldb_hist()->GetYaxis()->FindFixBin(hdr->triggerTimeNs); 
  return cfg->ldb_hist()->GetBinContent(binx,biny); 
}





