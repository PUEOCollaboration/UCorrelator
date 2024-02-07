#include "pueo/AnalysisConfig.h" 
#include "pueo/SystemResponse.h" 
#include "TFile.h" 
#include "TH2.h" 
#include "pueo/Version.h"

static const char * peakfinders[] = {"Abby","Bicubic","Gaussian","QuadraticFit9","QuadraticFit16","QuadraticFit25", "Histogram" }; 
static const char * responses[] = {"None","SingleBRotter","IndividualBRotter","HarmSignalOnly", "TUFFs", "A4noNotches", "A4ImpulseTUFFs", "A4OldImpulseTUFFs"};
//BinnedAnalysis Changes - JCF 9/27/2021
// added sammy  -----------------------------------------------
static const char * normalizations[] = {"None", "Magnitude", "Overlap magnitude", "Overlap-weighted magnitude"}; 
//-------------------------------------------------------------
//end BinnedAnalysis changes.

#ifdef ENABLE_LIBCONFIG
#include "libconfig.h++"
#define LOOKUP(X) cfg.lookupValue(#X,X) 

static void lookupEnum(libconfig::Config * cfg, const char * key, int * val, int N, const char ** allowed) 
{
  if (!cfg->exists(key)) return; 

  libconfig::Setting & set = cfg->lookup(key); 

  for (int i = 0; i < N; i++)
  {
    if (!strcasecmp(set, allowed[i]))
    {
      *val = i; 
      return; 
    }
  }

  fprintf(stderr,"Config Parse Error on line %d: %s must be one of: \n\t", set.getSourceLine(), key); 


  for (int i = 0; i < N; i++)
  {
    fprintf(stderr,"%s ", allowed[i]); 
  }

  fprintf(stderr, "\n But %s was given.", (const char * ) set); 

}


void pueo::UCorrelator::AnalysisConfig::loadFromFile(const char * config_file) 
{

  libconfig::Config cfg; 
  cfg.setAutoConvert(true); 
  cfg.readFile(config_file); 


  LOOKUP(correlator_nphi); 
  LOOKUP(correlator_ntheta); 
  LOOKUP(correlator_theta_lowest); 
  LOOKUP(correlator_theta_highest); 
  LOOKUP(enable_group_delay); 
  LOOKUP(zoomed_nphi); 
  LOOKUP(zoomed_ntheta); 
  LOOKUP(zoomed_dphi); 
  LOOKUP(zoomed_dtheta); 
  LOOKUP(zoomed_nant); 
  LOOKUP(combine_nantennas); 
  LOOKUP(combine_npad); 
  LOOKUP(saturation_threshold); 
  LOOKUP(peak_isolation_requirement); 
  LOOKUP(nmaxima); 
  LOOKUP(use_bin_center); 
  LOOKUP(bw_ndb);
  LOOKUP(spectral_fit_start);
  LOOKUP(spectral_fit_stop); 
  LOOKUP(noise_estimate_t0); 
  LOOKUP(noise_estimate_t1); 
  LOOKUP(response_npad); 
  LOOKUP(scale_by_cos_theta); 
  LOOKUP(use_offline_mask); 
  LOOKUP(max_mean_power_filtered); 
  LOOKUP(max_median_power_filtered); 
  LOOKUP(max_bottom_to_top_ratio); 
  LOOKUP(baseline_weight); 
  LOOKUP(max_peak_trigger_angle); 
  LOOKUP(min_peak_distance_from_unmasked); 
  LOOKUP(fill_blast_fraction); 
  LOOKUP(fill_channel_info); 
  LOOKUP(compute_shape_parameters); 
  LOOKUP(trace_to_continent); 
  LOOKUP(max_theta_adjustment); 
  LOOKUP(set_bottom_first);
  LOOKUP(delay_to_center);
  LOOKUP(only_use_usable);
  LOOKUP(use_antenna_level_snr);
  LOOKUP(stokes_fracI); 
  LOOKUP(use_forced_trigger_rms); 
  LOOKUP(use_coherent_spectra); 
  LOOKUP(combine_t0); 
  LOOKUP(combine_t1); 
  LOOKUP(use_hilbert_for_antenna_average); 
  LOOKUP(r_time_shift_correction); 
  LOOKUP(simulation_time_shift_correction); 
  LOOKUP(correlator_gain_correction); 
  const char * pols[] = {"horizontal", "vertical" }; 
  lookupEnum(&cfg, "start_pol", (int*) &start_pol, 2,pols); 
  lookupEnum(&cfg, "end_pol", (int*) &end_pol, 2,pols); 
  lookupEnum(&cfg, "fine_peak_finding_option", (int*) &fine_peak_finding_option, sizeof(peakfinders)/sizeof(char *), peakfinders); 
  lookupEnum(&cfg, "response_option", (int*) &response_option, sizeof(responses)/sizeof(char *), responses); 
  //BinnedAnalysis change - JCF 9/27/2021
  // added sammy
  lookupEnum(&cfg, "normalization_option", normalization_option, sizeof(normalizations)/sizeof(char *), normalizations); 
  //end BinnedAnalysis change
}

#else
void pueo::UCorrelator::AnalysisConfig::loadFromFile(const char * config_file) 
{

  fprintf(stderr, "Not compiled with support for reading config files. You need libconfig for that.\n"); 
  (void) config_file; 
}
#endif


const int wais_hpol_time_offset[5] = {0,0,0,93,11000}; 
const int wais_vpol_time_offset[5] = {0,0,0,-99757,1000}; 
const int siple_hpol_time_offset = -41;
const int siple_vpol_time_offset = +328;

pueo::UCorrelator::AnalysisConfig::AnalysisConfig(const char * config) 
  : 
    wais_hpol(wais_hpol_time_offset[version::get()], 800e3, 1e3), 
    wais_vpol(wais_vpol_time_offset[version::get()], 800e3, 1e3), 
    siple_hpol(siple_hpol_time_offset, 800e3, 1e3), 
    siple_vpol(siple_vpol_time_offset, 800e3, 1e3)
{
  correlator_nphi = 180; 
  correlator_ntheta = 100; 
  correlator_theta_lowest = 60; 
  correlator_theta_highest = 40; 
  correlation_gain_correction = 0; 
  enable_group_delay = true; 
  use_offline_mask = true; 
  zoomed_nphi = 40; 
  zoomed_ntheta = 40; 
  zoomed_dphi = 0.5; 
  zoomed_dtheta = 0.5; 
  zoomed_nant = 12; 


  combine_nantennas = 24; 
  combine_npad = 3; 

  saturation_threshold = 1500; 

  start_pol = pol::kHorizontal; 
  end_pol = pol::kVertical; 

  peak_isolation_requirement = 10; 

  fine_peak_finding_option = FinePeakFindingQuadraticFit25; 
  nmaxima = 2; 
  use_bin_center = false; 
  scale_by_cos_theta = false; 
  max_peak_trigger_angle = 0; 
  min_peak_distance_from_unmasked = -1; 

  bw_ndb = 6; //let the bandwidth search stop immediately
  spectral_fit_start = 0.22; 
  spectral_fit_stop =0.55; 

  noise_estimate_t0 = 70; 
  noise_estimate_t1 = 100; 

  response_option = ResponseNone; 
  response_string = 0; 
  response_npad = 50; 

  max_mean_power_filtered = 1e6; 
  max_median_power_filtered = 1e6; 
  max_bottom_to_top_ratio = 5; 
  max_theta_adjustment = 3; 

  baseline_weight = 0; 

  if (config) loadFromFile(config);

  deconvolution_method = &pueo::kDefaultDeconvolution; 
  
  the_ldb_hist = 0; 

  ldb_max_run = 160; 

  fill_blast_fraction = false; 
  fill_channel_info = true; 
  set_bottom_first = true;
  stokes_fracI = 0.2; 

  compute_shape_parameters = true;
  trace_to_continent = true;

  delay_to_center = true;
  r_time_shift_correction = true; 
  simulation_time_shift_correction = false; 
  cross_correlate_hv = 0; 

  use_forced_trigger_rms = true; 
  
  only_use_usable = false;
  use_antenna_level_snr = 0;

  use_coherent_spectra = false; 
  combine_t0 = -25; 
  combine_t1 = 425; 
  use_hilbert_for_antenna_average = true; 

}

static int nag = 0; 
TH2 * pueo::UCorrelator::AnalysisConfig::ldb_hist() const 
{
  if (the_ldb_hist) return the_ldb_hist; 

  TString fname; 
  fname.Form("%s/share/UCorrelator/ldbSelection.root", getenv("PUEO_UTIL_INSTALL_DIR")); 
  TFile f(fname); 

  if (!f.IsOpen())
  {
    if (nag++ < 3) 
      fprintf(stderr,"WARNING: data/ldbSelection.root doesn't exist. This is probably because you haven't generated it.\n. If you don't care about LDB pulsers, ignore this. Otherwise, you can generate one using the makeLDBSelection.C macro.\n"); 

    the_ldb_hist = 0; 
  }
  else
  {
    the_ldb_hist = (TH2*) f.Get("ldbHist")->Clone("LDBHistogram"); 
    the_ldb_hist->SetDirectory(0); 
  }
  return the_ldb_hist; 
}


const char * pueo::UCorrelator::AnalysisConfig::getPeakFindingString(FinePeakFindingOption_t opt) 
{
  return peakfinders[opt]; 
}

const char * pueo::UCorrelator::AnalysisConfig::getResponseString(ResponseOption_t opt) 
{

  return opt == ResponseCustomString ? 0 :  responses[opt]; 
}

//BinnedAnalysis addition - JCF 9/27/2021
// added sammy ------------------------------------------------------------------
const char * pueo::UCorrelator::AnalysisConfig::getNormalizationString(NormalizationOption_t opt) 
{
  return normalizations[opt]; 
}
// ------------------------------------------------------------------------------
//End BinnedAnalysis addition.

pueo::UCorrelator::AnalysisConfig::~AnalysisConfig()
{
  if (the_ldb_hist) delete the_ldb_hist; 
}

