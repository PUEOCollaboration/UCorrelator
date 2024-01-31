#include "pueo/Analyzer.h" 
#include "pueo/AnalysisConfig.h" 
#include "pueo/Polarimetry.h" 
#include "TPaveText.h"
#include "pueo/AnalysisWaveform.h"
#include "pueo/Version.h" 
//#include "pueo/FlightInfo.h" 
#include "TCanvas.h"
#include "pueo/ImpulsivityMeasure.h" 
#include "pueo/BandwidthMeasure.h" 
#include "TStyle.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "pueo/FilteredEvent.h"
#include "pueo/UsefulEvent.h"
#include "pueo/UCUtil.h"
#include "DigitalFilter.h" 
#include "FFTtools.h" 
#include "pueo/UsefulAttitude.h"
#include "pueo/RawHeader.h" 
#include "pueo/PeakFinder.h"
#include "pueo/UCFlags.h"
#include "pueo/ShapeParameters.h" 
#include "pueo/SpectrumParameters.h" 
#include "TF1.h" 
#include "TGraphErrors.h"
#include "pueo/UCImageTools.h"
#include <stdint.h>
#include "pueo/TimeDependentAverage.h" 

#ifdef UCORRELATOR_OPENMP
#include <omp.h>
#include "TThread.h" 

#define SECTIONS _Pragma("omp parallel sections")
#define SECTION _Pragma("omp section") 

#else 

#define SECTIONS if(true) 
#define SECTION if(true) 

#endif 

static pueo::UCorrelator::AnalysisConfig defaultConfig; 
static int instance_counter = 0; 

#ifndef DEG2RAD
#define DEG2RAD (TMath::Pi()/ 180) 
#endif

#ifndef RAD2DEG
#define RAD2DEG (180 / TMath::Pi()) 
#endif






  pueo::UCorrelator::Analyzer::Analyzer(const AnalysisConfig * conf, bool interactive_mode) 
: cfg(conf ? conf: &defaultConfig),
  corr(cfg->correlator_nphi,0,360,  cfg->correlator_ntheta, -cfg->correlator_theta_lowest, cfg->correlator_theta_highest, cfg->use_bin_center, cfg->scale_by_cos_theta, cfg->baseline_weight, cfg->correlation_gain_correction ) , 
  responses(cfg->response_option == AnalysisConfig::ResponseCustomString ? cfg->response_string : AnalysisConfig::getResponseString(cfg->response_option), cfg->response_npad, cfg->deconvolution_method), 
  wfcomb(cfg->combine_nantennas, cfg->combine_npad, true, cfg->response_option!=AnalysisConfig::ResponseNone, &responses), 
  wfcomb_xpol(cfg->combine_nantennas, cfg->combine_npad, true, cfg->response_option!=AnalysisConfig::ResponseNone, &responses), 
  wfcomb_filtered(cfg->combine_nantennas, cfg->combine_npad, false, cfg->response_option!=AnalysisConfig::ResponseNone, &responses), 
  wfcomb_xpol_filtered(cfg->combine_nantennas, cfg->combine_npad, false, cfg->response_option!=AnalysisConfig::ResponseNone, &responses), 
  interactive(interactive_mode)
{
#ifdef UCORRELATOR_OPENMP
  TThread::Initialize(); 
  ROOT::EnableThreadSafety(); 
#endif
  zoomed = new TH2D(TString::Format("zoomed_%d", instance_counter), "Zoomed!", cfg->zoomed_nphi, 0 ,1, cfg->zoomed_ntheta, 0, 1);
  zoomed->SetDirectory(0); 

  avg_spectra[0] = 0; 
  avg_spectra[1] = 0; 

  disallowedAnts[0] = 0;
  disallowedAnts[1] = 0;

  phiRange[0] = 0.;
  phiRange[1] = 0.;
  thetaRange[0] = 0.;
  thetaRange[1] = 0.;
  exclude = false;	
  dPhi = 0;
  dTheta = 0;
  trackSource = false;
  trackSun = false;
  sourceLon = 0; sourceLat = 0; sourceAlt = 0;

  corr.setGroupDelayFlag(cfg->enable_group_delay); 
  if (cfg->correlation_gain_correction) corr.setMaxAntennaMaxPhiDistance(3*cfg->correlation_gain_correction); 

  //BinnedAnalysis addition - JCF 9/27/2021
  corr.setNormalization(cfg->normalization_option);
  //End BinnedAnalysis addition.

  wfcomb.setGroupDelayFlag(cfg->enable_group_delay); 
  wfcomb_xpol.setGroupDelayFlag(cfg->enable_group_delay); 
  wfcomb_filtered.setGroupDelayFlag(cfg->enable_group_delay); 
  wfcomb_xpol_filtered.setGroupDelayFlag(cfg->enable_group_delay); 

  wfcomb.setBottomFirst(cfg->set_bottom_first);
  wfcomb_xpol.setBottomFirst(cfg->set_bottom_first);
  wfcomb_filtered.setBottomFirst(cfg->set_bottom_first);
  wfcomb_xpol_filtered.setBottomFirst(cfg->set_bottom_first);

  wfcomb.setDelayToCenter(cfg->delay_to_center);
  wfcomb_xpol.setDelayToCenter(cfg->delay_to_center);
  wfcomb_filtered.setDelayToCenter(cfg->delay_to_center);
  wfcomb_xpol_filtered.setDelayToCenter(cfg->delay_to_center);

  wfcomb.setRTimeShiftFlag(cfg->r_time_shift_correction);
  wfcomb_xpol.setRTimeShiftFlag(cfg->r_time_shift_correction);
  wfcomb_filtered.setRTimeShiftFlag(cfg->r_time_shift_correction);
  wfcomb_xpol_filtered.setRTimeShiftFlag(cfg->r_time_shift_correction);
  wfcomb.setSimulationTimeShiftFlag(cfg->simulation_time_shift_correction);
  wfcomb_xpol.setSimulationTimeShiftFlag(cfg->simulation_time_shift_correction);
  wfcomb_filtered.setSimulationTimeShiftFlag(cfg->simulation_time_shift_correction);
  wfcomb_xpol_filtered.setSimulationTimeShiftFlag(cfg->simulation_time_shift_correction);



  instance_counter++; 
  power_filter = new FFTtools::GaussianFilter(2,3) ; //TODO make this configurable

  if (interactive) 
  {
    correlation_maps[0] = 0;
    correlation_maps[1] = 0;

    for (int i = 0; i < cfg->nmaxima; i++)
    {
      zoomed_correlation_maps[0].push_back(0); 
      zoomed_correlation_maps[1].push_back(0); 

      for (int ipol =0; ipol <2;ipol++)
      {
        for (int ifilt = 0; ifilt <2; ifilt++) 
        {
          coherent[ipol][ifilt].push_back(new AnalysisWaveform); 
          deconvolved[ipol][ifilt].push_back(new AnalysisWaveform); 
          coherent_xpol[ipol][ifilt].push_back(new AnalysisWaveform); 
          deconvolved_xpol[ipol][ifilt].push_back(new AnalysisWaveform); 
          coherent_power[ipol][ifilt].push_back(new TGraphAligned); 
          deconvolved_power[ipol][ifilt].push_back(new TGraphAligned); 
          coherent_power_xpol[ipol][ifilt].push_back(new TGraphAligned); 
          deconvolved_power_xpol[ipol][ifilt].push_back(new TGraphAligned); 
        }
      }
    }
  }
}



static double computeCombinedRMS(double t, pueo::pol::pol_t pol, const pueo::UCorrelator::WaveformCombiner * wfcomb, int use_antenna_level_snr) 
{
  double rms = 0 ;

  for (int ant =0; ant < wfcomb->getNAntennas(); ant++) 
  {
    double ant_rms = pueo::UCorrelator::TimeDependentAverageLoader::getRMS(t, pol, wfcomb->getUsedAntennas()[ant]); 
    if(use_antenna_level_snr) rms += ant_rms; 
    else rms += ant_rms * ant_rms; 
  }

  if(use_antenna_level_snr) rms = rms / wfcomb->getNAntennas(); 
  else rms = sqrt(rms) / wfcomb->getNAntennas(); 

  return rms; 
}



void pueo::UCorrelator::Analyzer::analyze(const FilteredEvent * event, EventSummary * summary, const TruthEvent * truth) 
{

  const RawHeader * hdr = event->getHeader(); 
  //this is for time dependent responses (right now only TUFFs)	
  responses.checkTime(hdr->payloadTime);
  //we need a UsefulAdu5Pat for this event
  UsefulAttitude * pat =  (UsefulAttitude*) event->getGPS();  //unconstifying it .. hopefully that won't cause problems

  /* Initialize the summary */ 
  summary = new (summary) EventSummary(hdr, (UsefulAttitude*) event->getGPS(),truth); 

  //just right away store where the payload is
  summary->location.update(pat);

  //check if were blocking out/looking at a source
  if(trackSource)
  {
    pat->getThetaAndPhiWave(sourceLon, sourceLat, sourceAlt, thetaRange[0], phiRange[0]);
    phiRange[0] *= 180./TMath::Pi();
    thetaRange[0] *= -180./TMath::Pi();
    phiRange[1] = phiRange[0] + dPhi;
    thetaRange[1] = thetaRange[0] + dTheta;
    phiRange[0] -= dPhi;
    thetaRange[0] -= dTheta;
  }
  //sun is done separately because its easier this way
  if(trackSun)
  {
    pat->getSunPosition(summary->sun.phi, summary->sun.theta);
    phiRange[0] = summary->sun.phi - dPhi;
    phiRange[1] = summary->sun.phi + dPhi;
    thetaRange[0] = -summary->sun.theta - dTheta;
    thetaRange[1] = -summary->sun.theta + dTheta;
  }

  //check for saturation
  std::bitset<k::NUM_ANTS> saturated[2] = {0,0}; 
  event->checkSaturation( &saturated[pol::kHorizontal], 
      &saturated[pol::kVertical], 
      cfg->saturation_threshold); 

  //also disable missing sectors 
  pueo::UCorrelator::flags::checkEmpty(event->getUsefulEvent(), &saturated[pol::kHorizontal], &saturated[pol::kVertical]); 

  // loop over wanted polarizations 
  for (int pol = cfg->start_pol; pol <= cfg->end_pol; pol++) 
  {

    UShort_t triggeredPhi = event->getHeader()->phiTrigMask[pol];
    UShort_t triggeredPhiXpol =event->getHeader()->phiTrigMask[!pol];


    UShort_t maskedL2 = 0;   
    UShort_t maskedPhi = 0 ; 
    UShort_t maskedL2Xpol = 0;   
    UShort_t maskedPhiXpol = 0 ; 

    maskedPhi =  event->getHeader()->phiTrigMask[pol];
    maskedPhiXpol = event->getHeader()->phiTrigMask[!pol];


    //alright, we need a little song and dance here to combine Phi masks and L2 masks 
    //
    //  Someone should check my logic here. 
    //
    // An L2 mask effectively masks both phi sector N and N+1 (since it'll prevent L3 triggers in both). 
    // An L3 mask has contributions from both phi sector N and N-1. 
    //
    // I'm going to take the aggressive approach where an L3 mask means we mark a pointing hypothesis as masked
    // if it falls in phi sector N or N-1. 


#ifdef __cpp_static_assert
    static_assert(sizeof(maskedPhi == k::NUM_PHI),"masked phi must be same size as num phi "); 
#endif

    maskedPhi |= ( (maskedPhi >> 1) | ( maskedPhi << (k::NUM_PHI-1))) ; 
    //or with l2 mask
    maskedPhi |= maskedL2; 


#ifdef __cpp_static_assert
    static_assert(sizeof(maskedPhiXpol ==k::NUM_PHI),"masked phi xpol must be same size as num phi "); 
#endif

    //ditto for xpol 
    maskedPhiXpol |=  (maskedPhiXpol >> 1) | ((maskedPhiXpol << (k::NUM_PHI-1))); 
    maskedPhiXpol |= maskedL2Xpol; 

    //TODO: check this 
    TVector2 triggerAngle(0,0); 

    int ntriggered = __builtin_popcount(triggeredPhi); 
    int ntriggered_xpol = __builtin_popcount(triggeredPhiXpol); 

    UInt_t which_trigger = ntriggered ? triggeredPhi : triggeredPhiXpol; 
    int which_ntriggered = ntriggered ?: ntriggered_xpol; 

    const double phi_sector_width = 360. / k::NUM_PHI; 

    for (int i = 0; i < k::NUM_PHI; i++) 
    {
      if (which_trigger & (1 << i))
      {
        //TODO: this 45 is hardcoded here. Should come from GeomTool or something... 
        double ang = (i * phi_sector_width - 45) * TMath::Pi()/180;
        triggerAngle += TVector2(cos(ang), sin(ang)) / which_ntriggered; 
      }
    }

    double avgHwAngle = triggerAngle.Phi() * RAD2DEG; 

    // tell the correlator not to use saturated events or disallowed antennas and make the correlation map
    saturated[pol] |= disallowedAnts[pol];
//    if(cfg->only_use_usable) saturated[pol] |= ~FlightInfo::getUsableAntennas(hdr, event->getUsefulEvent(), pol::pol_t(pol));


    //if we are only considering antennas that are unmasked, figure out which ones are 
    std::bitset<k::NUM_ANTS> maskedAnts= 0; 
    if (cfg->min_peak_distance_from_unmasked >=0) 
    {
      for (uint64_t iphi = 0; iphi < 24; iphi++)
      {
        bool unmasked = !(maskedPhi & (1ul << iphi)) ; 

        for (int neighboring = 0; neighboring < cfg->min_peak_distance_from_unmasked; neighboring++)
        {
          if (unmasked) break; 
          unmasked = unmasked || !(maskedPhi & (1ul << ( (iphi + neighboring) % k::NUM_PHI))); 
          unmasked = unmasked || !(maskedPhi & (1ul << ( (iphi + neighboring + k::NUM_PHI - 1) % k::NUM_PHI))); 
        }

        if (!unmasked) 
        {
          maskedAnts |= 1ul << iphi; 
          maskedAnts |= 1ul << (iphi+k::NUM_PHI); 
          maskedAnts |= 1ul << (iphi+2*k::NUM_PHI); 
          maskedAnts |= 1ul << (iphi+3*k::NUM_PHI); 
        }
      }
    }

    corr.setDisallowedAntennas(saturated[pol] | disallowedAnts[pol] | maskedAnts); 
    //BinnedAnalysis addition - JCF 9/29/2021
    corr.disallowAntennas(saturated[pol]);
    //End BinnedAnalysis addition
    corr.compute(event, pol::pol_t(pol)); 

    //compute RMS of correlation map 
    //    maprms = corr.getHist()->GetRMS(3); //This doesn't work!  Probably because ROOT is dumb


    if (cfg->max_peak_trigger_angle)
    {
      phiRange[0] = FFTtools::wrap(avgHwAngle - cfg->max_peak_trigger_angle,360); 
      phiRange[1] = FFTtools::wrap(avgHwAngle + cfg->max_peak_trigger_angle,360); 
    }

    // Find the isolated peaks in the image 
    peakfinder::RoughMaximum maxima[cfg->nmaxima]; 
    int npeaks = pueo::UCorrelator::peakfinder::findIsolatedMaxima((const TH2D*) corr.getHist(),
        cfg->peak_isolation_requirement,
        cfg->nmaxima, maxima, phiRange[0], phiRange[1], 
        thetaRange[0], thetaRange[1], exclude, cfg->use_bin_center); 
    //    printf("npeaks: %d\n", npeaks); 
    summary->nPeaks[pol] = npeaks; 

    rough_peaks[pol].clear(); 


    //get the average spectra 
    if (!avg_spectra[pol]) 
    {
      avg_spectra[pol] = new TGraph; 
      avg_spectra[pol]->GetXaxis()->SetTitle("Frequency (GHz)"); 
      avg_spectra[pol]->GetYaxis()->SetTitle("Power (dBish)"); 
      avg_spectra[pol]->SetTitle(TString::Format("Average spectra for %s", pol == pol::kHorizontal ? "HPol" : "VPol")); 
    }

    event->getMedianSpectrum(avg_spectra[pol], pol::pol_t(pol),0.5); 

    //optionally fill the channel info 
    if (cfg->fill_channel_info) 
    {
      fillChannelInfo(event, summary);
    }



    int rms_lims[4]; 
    rms_lims[0] = (cfg->correlator_nphi/6 + maxima[0].xbin) % cfg->correlator_nphi; 
    rms_lims[1] = (-cfg->correlator_nphi/6 + maxima[0].xbin) % cfg->correlator_nphi; 
    rms_lims[2] = 1; 
    rms_lims[3] = cfg->correlator_ntheta; 
    maprms = pueo::UCorrelator::getZRMS(corr.getHist(), rms_lims);

    // Loop over found peaks 
    for (int i = 0; i < npeaks; i++) 
    {
      // zoom in on the values 
      //      printf("rough phi:%f, rough theta: %f\n", maxima[i].x, -maxima[i].y); 


      fillPointingInfo(maxima[i].x, maxima[i].y, &summary->peak[pol][i], pat, avgHwAngle, triggeredPhi, maskedPhi, triggeredPhiXpol, maskedPhiXpol); 

      if (interactive) 
      {
        if (zoomed_correlation_maps[pol][i]) delete zoomed_correlation_maps[pol][i]; 
        zoomed_correlation_maps[pol][i] = new gui::Map(*zoomed, event, &wfcomb, &wfcomb_filtered,pol::pol_t(pol), summary); 
        zoomed_correlation_maps[pol][i]->SetName(TString::Format("zoomed_%d_%d", pol,i)); 
      }


      //fill in separation 
      summary->peak[pol][i].phi_separation = 1000; 
      for (int j = 0; j < i; j++)
      {
        summary->peak[pol][i].phi_separation = TMath::Min(summary->peak[pol][i].phi_separation, fabs(FFTtools::wrap(summary->peak[pol][i].phi - summary->peak[pol][j].phi, 360, 0))); 
      }

      //      printf("phi:%f, theta:%f\n", summary->peak[pol][i].phi, summary->peak[pol][i].theta); 


    }

    for (int i = 0; i < npeaks; i++) 
    {
      rough_peaks[pol].push_back(std::pair<double,double>(maxima[i].x, maxima[i].y)); 
      //now make the combined waveforms 
      if(cfg->use_antenna_level_snr > 0)
      {
        wfcomb.setCheckVpp(cfg->use_antenna_level_snr);
        wfcomb_xpol.setCheckVpp(cfg->use_antenna_level_snr);
        wfcomb_filtered.setCheckVpp(cfg->use_antenna_level_snr);
        wfcomb_xpol_filtered.setCheckVpp(cfg->use_antenna_level_snr);
      }


      SECTIONS
      {
        SECTION
          wfcomb.combine(summary->peak[pol][i].phi, summary->peak[pol][i].theta, event, (pol::pol_t) pol, saturated[pol], 
              cfg->combine_t0, cfg->combine_t1, &summary->peak[pol][i].antennaPeakAverage, cfg->use_hilbert_for_antenna_average); 
        SECTION
          wfcomb_xpol.combine(summary->peak[pol][i].phi, summary->peak[pol][i].theta, event, (pol::pol_t) (1-pol), saturated[pol], cfg->combine_t0, cfg->combine_t1); 
        SECTION
          wfcomb_filtered.combine(summary->peak[pol][i].phi, summary->peak[pol][i].theta, event, (pol::pol_t) pol, saturated[pol], cfg->combine_t0, cfg->combine_t1); 
        SECTION
          wfcomb_xpol_filtered.combine(summary->peak[pol][i].phi, summary->peak[pol][i].theta, event, (pol::pol_t) (1-pol), saturated[pol], cfg->combine_t0, cfg->combine_t1); 
      }

      double rms =  cfg->use_forced_trigger_rms ? computeCombinedRMS(event->getHeader()->triggerTime, (pol::pol_t) pol, &wfcomb, cfg->use_antenna_level_snr) : 0 ; 

      SECTIONS 
      {
        SECTION
          fillWaveformInfo(wfcomb.getCoherent(), wfcomb_xpol.getCoherent(), wfcomb.getCoherentAvgSpectrum(), &summary->coherent[pol][i], (pol::pol_t) pol, rms, wfcomb.getMaxAntennaVpp()); 
        SECTION
          fillWaveformInfo(wfcomb.getDeconvolved(), wfcomb_xpol.getDeconvolved(), wfcomb.getDeconvolvedAvgSpectrum(), &summary->deconvolved[pol][i],  (pol::pol_t)pol, 0); 
        SECTION
          fillWaveformInfo(wfcomb_filtered.getCoherent(), wfcomb_xpol_filtered.getCoherent(), wfcomb_filtered.getCoherentAvgSpectrum(), &summary->coherent_filtered[pol][i], (pol::pol_t) pol, rms, wfcomb_filtered.getMaxAntennaVpp()); 
        SECTION
          fillWaveformInfo(wfcomb_filtered.getDeconvolved(), wfcomb_xpol_filtered.getDeconvolved(), wfcomb_filtered.getDeconvolvedAvgSpectrum(), &summary->deconvolved_filtered[pol][i],  (pol::pol_t)pol, 0); 
      }


      if (interactive) //copy everything
      {
        coherent[pol][0][i]->~AnalysisWaveform(); 
        coherent[pol][0][i] = new (coherent[pol][0][i]) AnalysisWaveform(*wfcomb.getCoherent()); 

        coherent[pol][1][i]->~AnalysisWaveform(); 
        coherent[pol][1][i] = new (coherent[pol][1][i]) AnalysisWaveform(*wfcomb_filtered.getCoherent()); 

        coherent_power[pol][0][i]->~TGraphAligned(); 
        coherent_power[pol][0][i] = new (coherent_power[pol][0][i]) TGraphAligned(*wfcomb.getCoherentAvgSpectrum()); 
        coherent_power[pol][0][i]->dBize(); 

        coherent_power[pol][1][i]->~TGraphAligned(); 
        coherent_power[pol][1][i] = new (coherent_power[pol][1][i]) TGraphAligned(*wfcomb_filtered.getCoherentAvgSpectrum()); 
        coherent_power[pol][1][i]->dBize(); 


        if (wfcomb.getDeconvolved())
        {
          deconvolved[pol][0][i]->~AnalysisWaveform(); 
          deconvolved[pol][0][i] = new (deconvolved[pol][0][i]) AnalysisWaveform(*wfcomb.getDeconvolved()); 
          deconvolved[pol][0][i]->updateEven()->SetLineColor(2); 
          deconvolved[pol][0][i]->updateEven()->SetMarkerColor(2); 

          deconvolved[pol][1][i]->~AnalysisWaveform(); 
          deconvolved[pol][1][i] = new (deconvolved[pol][1][i]) AnalysisWaveform(*wfcomb_filtered.getDeconvolved()); 
          deconvolved[pol][1][i]->updateEven()->SetLineColor(2); 
          deconvolved[pol][1][i]->updateEven()->SetMarkerColor(2); 

          deconvolved_power[pol][0][i]->~TGraphAligned(); 
          deconvolved_power[pol][0][i] = new (deconvolved_power[pol][0][i]) TGraphAligned(*wfcomb.getDeconvolvedAvgSpectrum()); 
          deconvolved_power[pol][0][i]->dBize(); 
          deconvolved_power[pol][0][i]->SetLineColor(2); 

          deconvolved_power[pol][1][i]->~TGraphAligned(); 
          deconvolved_power[pol][1][i] = new (deconvolved_power[pol][1][i]) TGraphAligned(*wfcomb_filtered.getDeconvolvedAvgSpectrum()); 
          deconvolved_power[pol][1][i]->dBize(); 
          deconvolved_power[pol][1][i]->SetLineColor(2); 

          interactive_deconvolved = true; 
        }
        else
        {
          interactive_deconvolved = false; 
        }


        coherent_xpol[pol][0][i]->~AnalysisWaveform(); 
        coherent_xpol[pol][0][i] = new (coherent_xpol[pol][0][i]) AnalysisWaveform(*wfcomb_xpol.getCoherent()); 
        coherent_xpol[pol][0][i]->updateEven()->SetLineColor(11); 
        coherent_xpol[pol][0][i]->updateEven()->SetLineStyle(3); 

        coherent_xpol[pol][1][i]->~AnalysisWaveform(); 
        coherent_xpol[pol][1][i] = new (coherent_xpol[pol][1][i]) AnalysisWaveform(*wfcomb_xpol_filtered.getCoherent()); 
        coherent_xpol[pol][1][i]->updateEven()->SetLineColor(11); 
        coherent_xpol[pol][1][i]->updateEven()->SetLineStyle(3); 


        coherent_power_xpol[pol][0][i]->~TGraphAligned(); 
        coherent_power_xpol[pol][0][i] = new (coherent_power_xpol[pol][0][i]) TGraphAligned(*wfcomb_xpol.getCoherentAvgSpectrum()); 
        coherent_power_xpol[pol][0][i]->dBize(); 
        coherent_power_xpol[pol][0][i]->SetLineStyle(3); 
        coherent_power_xpol[pol][0][i]->SetLineColor(11); 

        coherent_power_xpol[pol][1][i]->~TGraphAligned(); 
        coherent_power_xpol[pol][1][i] = new (coherent_power_xpol[pol][1][i]) TGraphAligned(*wfcomb_xpol_filtered.getCoherentAvgSpectrum()); 
        coherent_power_xpol[pol][1][i]->dBize(); 
        coherent_power_xpol[pol][1][i]->SetLineStyle(3); 
        coherent_power_xpol[pol][1][i]->SetLineColor(11); 


        if (wfcomb_xpol.getDeconvolved())
        {
          deconvolved_xpol[pol][0][i]->~AnalysisWaveform(); 
          deconvolved_xpol[pol][0][i] = new (deconvolved_xpol[pol][0][i]) AnalysisWaveform(*wfcomb_xpol.getDeconvolved()); 
          deconvolved_xpol[pol][0][i]->updateEven()->SetLineColor(45); 
          deconvolved_xpol[pol][0][i]->updateEven()->SetMarkerColor(45); 
          deconvolved_xpol[pol][0][i]->updateEven()->SetLineStyle(3); 

          deconvolved_power_xpol[pol][0][i]->~TGraphAligned(); 
          deconvolved_power_xpol[pol][0][i] = new (deconvolved_power_xpol[pol][0][i]) TGraphAligned(*wfcomb_xpol.getDeconvolvedAvgSpectrum()); 
          deconvolved_power_xpol[pol][0][i]->dBize(); 
          deconvolved_power_xpol[pol][0][i]->SetLineStyle(3); 
          deconvolved_power_xpol[pol][0][i]->SetLineColor(46); 

          deconvolved_xpol[pol][1][i]->~AnalysisWaveform(); 
          deconvolved_xpol[pol][1][i] = new (deconvolved_xpol[pol][1][i]) AnalysisWaveform(*wfcomb_xpol_filtered.getDeconvolved()); 
          deconvolved_xpol[pol][1][i]->updateEven()->SetLineColor(45); 
          deconvolved_xpol[pol][1][i]->updateEven()->SetMarkerColor(45); 
          deconvolved_xpol[pol][1][i]->updateEven()->SetLineStyle(3); 

          deconvolved_power_xpol[pol][1][i]->~TGraphAligned(); 
          deconvolved_power_xpol[pol][1][i] = new (deconvolved_power_xpol[pol][1][i]) TGraphAligned(*wfcomb_xpol_filtered.getDeconvolvedAvgSpectrum()); 
          deconvolved_power_xpol[pol][1][i]->dBize(); 
          deconvolved_power_xpol[pol][1][i]->SetLineStyle(3); 
          deconvolved_power_xpol[pol][1][i]->SetLineColor(46); 



          interactive_xpol_deconvolved = true; 
        }
        else
        {
          interactive_xpol_deconvolved = false; 
        }
      }
    }
    if (interactive) 
    {
      if (correlation_maps[pol]) delete correlation_maps[pol];
      correlation_maps[pol] = new gui::Map(*corr.getHist(), event, &wfcomb, &wfcomb_filtered,pol::pol_t(pol), summary ); 
    }
  }

  fillFlags(event, summary, pat); 


  if (truth)
  { 
    //    SECTIONS
    {
      //      SECTION
      if(cfg->use_antenna_level_snr > 0)
      {
        wfcomb.setCheckVpp(cfg->use_antenna_level_snr);
        wfcomb_xpol.setCheckVpp(cfg->use_antenna_level_snr);
      }

      wfcomb.combine(summary->mc.phi, summary->mc.theta, event, pol::kHorizontal, 0, cfg->combine_t0, cfg->combine_t1); 
      //      SECTION 
      wfcomb_xpol.combine(summary->mc.phi, summary->mc.theta, event, pol::kVertical, 0, cfg->combine_t0, cfg->combine_t1); 
    }

    double hpol_rms =  cfg->use_forced_trigger_rms ? computeCombinedRMS(event->getHeader()->triggerTime, pol::kHorizontal, &wfcomb, cfg->use_antenna_level_snr) : 0 ; 
    double vpol_rms =  cfg->use_forced_trigger_rms ? computeCombinedRMS(event->getHeader()->triggerTime, pol::kVertical, &wfcomb_xpol, cfg->use_antenna_level_snr) : 0 ; 

    //    SECTIONS
    {
      //      SECTION
      fillWaveformInfo(wfcomb.getCoherent(), wfcomb_xpol.getCoherent(), wfcomb.getCoherentAvgSpectrum(), &(summary->mc.wf[pol::kHorizontal]), pol::kHorizontal, hpol_rms, wfcomb.getMaxAntennaVpp()); 
      //      SECTION
      fillWaveformInfo(wfcomb_xpol.getCoherent(), wfcomb.getCoherent(), wfcomb_xpol.getCoherentAvgSpectrum(), &(summary->mc.wf[pol::kVertical]), pol::kVertical, vpol_rms, wfcomb_xpol.getMaxAntennaVpp()); 
    }
  }

  if (interactive) last = *summary; 

}

static bool outside(const TH2 * h, double x, double y) 
{

  return x > h->GetXaxis()->GetXmax() || 
    x < h->GetXaxis()->GetXmin() || 
    y < h->GetYaxis()->GetXmin() ||
    y > h->GetYaxis()->GetXmax(); 

}

void pueo::UCorrelator::Analyzer::fillPointingInfo(double rough_phi, double rough_theta, EventSummary::PointingHypothesis * point, 
    UsefulAttitude * pat, double hwPeakAngle, UShort_t triggered_sectors, UShort_t masked_sectors, UShort_t triggered_sectors_xpol, UShort_t masked_sectors_xpol)
{
  corr.computeZoomed(rough_phi, rough_theta, cfg->zoomed_nphi, cfg->zoomed_dphi,  cfg->zoomed_ntheta, cfg->zoomed_dtheta, cfg->zoomed_nant, zoomed); 

  //get pointer to the pointing hypothesis we are about to fill 

  // This will fill in phi, theta, value, var_theta, var_phi and covar 

  peakfinder::FineMaximum max; 
  switch (cfg->fine_peak_finding_option)
  {
    case AnalysisConfig::FinePeakFindingAbby: 
      pueo::UCorrelator::peakfinder::doInterpolationPeakFindingAbby(zoomed, &max); 
      break; 
    case AnalysisConfig::FinePeakFindingBicubic: 
      pueo::UCorrelator::peakfinder::doInterpolationPeakFindingBicubic(zoomed, &max); 
      break; 
    case AnalysisConfig::FinePeakFindingHistogram: 
      pueo::UCorrelator::peakfinder::doPeakFindingHistogram(zoomed, &max); 
      break; 
    case AnalysisConfig::FinePeakFindingQuadraticFit16: 
      pueo::UCorrelator::peakfinder::doPeakFindingQuadratic16(zoomed, &max); 
      break; 
    case AnalysisConfig::FinePeakFindingQuadraticFit25: 
      pueo::UCorrelator::peakfinder::doPeakFindingQuadratic25(zoomed, &max); 
      break; 
    case AnalysisConfig::FinePeakFindingGaussianFit: 
      pueo::UCorrelator::peakfinder::doPeakFindingGaussian(zoomed, &max); 
      break; 
    case AnalysisConfig::FinePeakFindingQuadraticFit9: 
    default: 
      pueo::UCorrelator::peakfinder::doPeakFindingQuadratic9(zoomed, &max); 
      break; 
  }; 


  //Check to make sure that fine max isn't OUTSIDE of zoomed window
  // If it is, revert to very stupid method  of just using histogram 

  if (outside(zoomed, max.x, max.y))
  {
    pueo::UCorrelator::peakfinder::doPeakFindingHistogram(zoomed, &max); 
  }




  max.copyToPointingHypothesis(point); 

  //snr is ratio of point value to map rms
  point->snr = point->value / maprms; 
  point->mapRMS = maprms;
  point->dphi_rough = FFTtools::wrap(point->phi - rough_phi, 360,0); 
  point->dtheta_rough = FFTtools::wrap(point->theta - (-rough_theta), 360,0); //sign reversal. doh. 

  point->hwAngle = FFTtools::wrap(point->phi - hwPeakAngle,360,0); 

  //TODO: I don't believe this really yet
  int sector = 2+fmod(point->phi + 11.25,360) / 22.5; 

  point->masked = masked_sectors & ( 1 << sector); 
  point->triggered = triggered_sectors & ( 1 << sector); 
  point->masked_xpol = masked_sectors_xpol & ( 1 << sector); 
  point->triggered_xpol = triggered_sectors_xpol & ( 1 << sector); 

  if(cfg->trace_to_continent)
  {
    //Compute intersection with continent, or set values to -9999 if no intersection
    if (!pat->traceBackToContinent3(point->phi * DEG2RAD, point->theta * DEG2RAD, 
          &point->longitude, &point->latitude, &point->altitude, &point->theta_adjustment_needed)
        || point->theta_adjustment_needed > cfg->max_theta_adjustment * DEG2RAD) 
    {
      point->latitude = -9999; 
      point->longitude = -9999;  
      point->altitude = -9999; 
      point->distanceToSource = -9999; 
      point->theta_adjustment_needed = -9999; 
    }
    else
    {
      point->distanceToSource=pat->getDistanceFromSource(point->latitude, point->longitude, point->altitude); 
      point->theta_adjustment_needed *= RAD2DEG; 
    }
  }
  else
  {
    point->latitude = 0; 
    point->longitude = 0;  
    point->altitude = 0; 
    point->distanceToSource = 0; 
    point->theta_adjustment_needed = 0; 
  }
}


void pueo::UCorrelator::Analyzer::fillWaveformInfo(const AnalysisWaveform * wf, const AnalysisWaveform * xpol_wf, const TGraph* pwr, EventSummary::WaveformInfo * info, pol::pol_t pol, double rms, double vpp)
{
  if (!wf || wf->Neven() == 0)
  {
    if (wf && !wf->Neven()) 
      fprintf(stderr,"wf passed to fillWaveformInfo has no points\n");  

    *info = {}; 
    return;
  }
  const TGraphAligned * even = wf->even(); 
  const TGraphAligned * xpol_even= xpol_wf->even(); 
  const TGraphAligned* hilbert = wf->hilbertEnvelope(); //BinnedAnalysis addition - JCF 9/27/2021
  int peakBin;

  int peakHilbertBin; 
  info->peakVal = FFTtools::getPeakVal((TGraph*) even,&peakBin); 
  info->xPolPeakVal = FFTtools::getPeakVal( xpol_even); 
  info->peakHilbert = FFTtools::getPeakVal((TGraph*) wf->hilbertEnvelope(),&peakHilbertBin); 
  double minHilbert = *std::min_element(wf->hilbertEnvelope()->GetY(), wf->hilbertEnvelope()->GetY() + wf->Neven()); 

  info->xPolPeakHilbert = FFTtools::getPeakVal((TGraph*) xpol_wf->hilbertEnvelope()); 
  info->numAntennasInCoherent = cfg->combine_nantennas; 

  info->totalPower = even->getSumV2(); 
  info->totalPowerXpol = xpol_even->getSumV2(); 
  info->peakTime = even->GetX()[peakHilbertBin]; 

  double hilbertRange = info->peakHilbert - minHilbert; 

  if(cfg->compute_shape_parameters)
  {
    info->riseTime_10_90 = shape::getRiseTime((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.1*hilbertRange, minHilbert + 0.9*hilbertRange,peakHilbertBin); 
    info->riseTime_10_50 = shape::getRiseTime((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.1*hilbertRange, minHilbert + 0.5*hilbertRange,peakHilbertBin); 
    info->fallTime_90_10 = shape::getFallTime((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.1*hilbertRange, minHilbert + 0.9*hilbertRange,peakHilbertBin); 
    info->fallTime_50_10 = shape::getFallTime((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.1*hilbertRange, minHilbert + 0.5*hilbertRange,peakHilbertBin); 

    int ifirst, ilast; 
    info->width_50_50 = shape::getWidth((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.5*hilbertRange, &ifirst, &ilast,peakHilbertBin); 
    info->power_50_50 = even->getSumV2(ifirst, ilast); 
    even->getMoments(sizeof(info->peakMoments)/sizeof(double), info->peakTime, info->peakMoments); 
    info->width_10_10 = shape::getWidth((TGraph*) wf->hilbertEnvelope(), minHilbert + 0.1*hilbertRange, &ifirst, &ilast,peakHilbertBin); 
    info->power_10_10 = even->getSumV2(ifirst, ilast); 
  }
  if(!cfg->compute_shape_parameters)
  {
    info->riseTime_10_90 = 0; 
    info->riseTime_10_50 = 0; 
    info->fallTime_90_10 = 0; 
    info->fallTime_50_10 = 0; 

    info->width_50_50 = 0; 
    info->power_50_50 = 0; 
    info->width_10_10 = 0; 
    info->power_10_10 = 0; 
  }

  //BinnedAnalysis additions - JCF 9/27/2021
  // fix fail on unequal lengths - sammy
  int nMin = std::min(even->GetN(), xpol_even->GetN());
  if (pol == pol::kHorizontal)
  {
  	FFTtools::stokesParameters(nMin,
				   even->GetY(),
				   wf->hilbertTransform()->even()->GetY(),
				   xpol_even->GetY(),
				   xpol_wf->hilbertTransform()->even()->GetY(),
				   &(info->I), &(info->Q), &(info->U), &(info->V));
  }
  else
  {
	FFTtools::stokesParameters(nMin,
				   xpol_even->GetY(),
				   xpol_wf->hilbertTransform()->even()->GetY(),
				   even->GetY(),
				   wf->hilbertTransform()->even()->GetY(),
				   &(info->I), &(info->Q), &(info->U), &(info->V));
  }

  double dt = wf->deltaT();
  double t0 = even->GetX()[0];

  int i0 = TMath::Max(0.,floor((cfg->noise_estimate_t0 - t0)/dt));
  int i1 = TMath::Min(even->GetN()-1.,ceil((cfg->noise_estimate_t1 - t0)/dt));
  int n = i1 - i0 + 1;

  // changed by sammy to use SNR similar to B. Dailey, A. Vieregg, for like-to-like comparison
  double rmsSam;
  if (n<0) {
	printf("   -----Warning! - finagling index limits for rms calculation \n");
	printf("         even->GetN()=%i \n", even->GetN());
	printf("         %d-%d -> %d \n", i0, i1, n);
	rmsSam = TMath::RMS(-n, even->GetY() + i0 + n);
	printf("    rms value is %f \n", rmsSam);
  } else {
	rmsSam = TMath::RMS(n, even->GetY() + i0);
  }
  double thatSnr = info->peakVal / rmsSam; 

  // new code sammy
  // get rms of a 1/5 sampling of the waveform, use the end 1/5 since our pulses are early in the sampling window
  int k0 = (int)(0.75 * even->GetN());
  int k1 = (int)(0.95 * even->GetN());
  double thisNoise = TMath::RMS(k1-k0, &(even->GetY()[k0]));
  double thisSig = 0.5 * even->pk2pk(0, 0, NULL, NULL);
  double thisSnr = thisSig / thisNoise;

  // third way: use (noiseInterval) seconds starting with the first nonzero bucket
  // fourth way: integrate the hilbert envelope acroos the peak for S; for N, do like third way except use hilbert env instead of waveform
  float noiseInterval = 10.;
  int j0 = 0;
  for (; j0<even->GetN() && even->GetY()[j0]==0; ++j0) {}
  int j1 = j0 + noiseInterval/dt;
  float otherSnr = 1.0;
  float snr3 = 1.0;
  if (j1 >= even->GetN()) {
	printf("warning: could not obtain noise sample interval \n");
  } else {
	rmsSam = TMath::RMS(j1-j0+1, even->GetY() + j0);
	otherSnr = thisSig / rmsSam;
	double noise3 = hilbert->Integral(j0, j1);
	noise3 /= (j1-j0+1);
	float hilbPeakTime = 0;
	hilbPeakTime = info->peakTime;

	int hilbPeakIndex = (info->peakTime - hilbert->GetX()[0])/dt;
	int j3 = hilbPeakIndex - 7.5/dt;
	int j4 = hilbPeakIndex + 12.5/dt;
	double sig3 = hilbert->Integral(j3, j4);
	sig3 /= (j4-j3+1);
	snr3 = sig3/noise3;
	int wk;
  }

  info->snrStafford = thatSnr; //LE 5-17 should I be messing with this?
  info->snrFifthRMS = thisSnr;
  info->snrTenRMS = otherSnr;
  info->snrHilbert = snr3;
  //End BinnedAnalysis additions.

  polarimetry::StokesAnalysis stokes( pol == pol::kHorizontal ? wf : xpol_wf,  pol == pol::kHorizontal ? xpol_wf: wf, cfg->cross_correlate_hv); 
  info->I = stokes.getAvgI(); 
  info->Q = stokes.getAvgQ(); 
  info->U = stokes.getAvgU(); 
  info->V = stokes.getAvgV(); 
  info->NPointsMaxStokes = stokes.computeWindowedAverage(cfg->stokes_fracI, &info->max_dI, &info->max_dQ, &info->max_dU, &info->max_dV, &info->polErr); 

  TGraph distance_cdf; 
  info->impulsivityMeasure = impulsivity::impulsivityMeasure(wf, &distance_cdf); 
  info->bandwidthMeasure = bandwidth::lowness(wf); 

  //fill in narrowest widths
  for (int iw = 0; iw < EventSummary::numFracPowerWindows; iw++)
  {
    int half_width =  TMath::BinarySearch(distance_cdf.GetN(), distance_cdf.GetY(), 0.1 * (iw+1)); 
    info->fracPowerWindowBegins[iw] = peakHilbertBin < half_width ? wf->even()->GetX()[0] : wf->even()->GetX()[peakHilbertBin - half_width]; 
    info->fracPowerWindowEnds[iw] = peakHilbertBin + half_width >= distance_cdf.GetN() ?  wf->even()->GetX()[distance_cdf.GetN()-1] : wf->even()->GetX()[half_width + peakHilbertBin];
  }

  if (!rms) 
  {
    double dt = wf->deltaT(); 
    double t0 = even->GetX()[0]; 

    int i0 = TMath::Max(0.,floor((cfg->noise_estimate_t0 - t0)/dt)); 
    int i1 = TMath::Min(even->GetN()-1.,ceil((cfg->noise_estimate_t1 - t0)/dt)); 
    int n = i1 - i0 + 1; 
    //  printf("%d-%d -> %d \n", i0, i1, n); 
    if (n < 0 || n > even->GetN()) n = even->GetN();
    rms = TMath::RMS(n, even->GetY() + i0); 
  }

  info->snr = info->peakHilbert / rms;
  if(vpp > 0 && cfg->use_antenna_level_snr) info->snr = vpp/(2*rms);

  if(cfg->use_coherent_spectra) pwr = wf->powerdB(); 

  TGraphAligned power(pwr->GetN(),pwr->GetX(),pwr->GetY()); 

  if (!cfg->use_coherent_spectra) power.dBize(); 

  if (power_filter)
  {
    power_filter->filterGraph(&power); 
  }

  spectrum::fillSpectrumParameters(&power, avg_spectra[pol], info, cfg); 
}

/** 
 * Fill Peng's ChannelInfo object, it might come in handy...
 */
void pueo::UCorrelator::Analyzer::fillChannelInfo(const FilteredEvent* event, EventSummary* summary){
  for(int polInd=0; polInd < pol::kNotAPol; polInd++){
    pol::pol_t pol = (pol::pol_t) polInd;
#ifdef UCORRELATOR_OPENMP
#pragma omp parallel for
#endif
    for (int ant=0; ant<k::NUM_HORNS; ant++)
    {
      const AnalysisWaveform* wf = event->getFilteredGraph(ant, pol);
      const TGraphAligned* hilbertEnvelope= wf->hilbertEnvelope();
      const TGraphAligned* gr = wf->even();
      double rms =  cfg->use_forced_trigger_rms ? pueo::UCorrelator::TimeDependentAverageLoader::getRMS( event->getHeader()->triggerTime, pol, ant) : gr->GetRMS(2);

      summary->channels[polInd][ant].rms = rms; 
      summary->channels[polInd][ant].avgPower = gr->getSumV2() / gr->GetN();
      summary->channels[polInd][ant].peakHilbert = hilbertEnvelope->peakVal();
      summary->channels[polInd][ant].snr = FFTtools::getPeakVal(gr) / rms; 
    }
  }
}


pueo::UCorrelator::Analyzer::~Analyzer()
{

  delete zoomed; 
  if (interactive)
  {
    delete correlation_maps[0];
    delete correlation_maps[1]; 


    for (int pol = 0; pol < 2; pol++)
    {
      for (int i = 0; i < cfg->nmaxima; i++)
      {
        delete zoomed_correlation_maps[pol][i];

        for (int ifilt = 0; ifilt < 2; ifilt++) 
        {

          delete coherent[pol][ifilt][i];
          delete deconvolved[pol][ifilt][i];

          delete coherent_xpol[pol][ifilt][i];
          delete deconvolved_xpol[pol][ifilt][i];

          delete coherent_power[pol][ifilt][i];
          delete deconvolved_power[pol][ifilt][i];

          delete coherent_power_xpol[pol][ifilt][i];
          delete deconvolved_power_xpol[pol][ifilt][i];
        }

      }
    }

    clearInteractiveMemory(1); 
  }

  if (power_filter)
    delete power_filter; 
}

void pueo::UCorrelator::Analyzer::clearInteractiveMemory(double frac) const
{

  for (unsigned i = (1-frac) * delete_list.size(); i < delete_list.size(); i++) 
  {
    delete delete_list[i]; 
  }

  delete_list.clear(); 
}


/* Nevermind this... wanted to zoom in on analyzer canvas on click, but too much work :( 
   static void setOnClickHandler(TPad * pad) 
   {

   }
   */





void pueo::UCorrelator::Analyzer::drawSummary(TPad * ch, TPad * cv, int draw_filtered) const
{
  TPad * pads[2] = {ch,cv}; 

  clearInteractiveMemory(); 

  //  gStyle->SetOptStat(0); 
  for (int ipol = cfg->start_pol; ipol <= cfg->end_pol; ipol++)
  {
    if (!pads[ipol])
    {
      pads[ipol] = new TCanvas(ipol == 0 ? "analyzer_ch" : "analyzer_cv", ipol == 0 ? "hpol" : "vpol",1920,500); 
    }

    pads[ipol]->Clear(); 
    pads[ipol]->Divide(2,1); 

    pads[ipol]->cd(1)->Divide(1,2); 

    pads[ipol]->cd(1)->cd(1); 
    correlation_maps[ipol]->SetTitle(ipol == 0 ? "HPol map" : "VPol map" ); 
    correlation_maps[ipol]->addRough(rough_peaks[ipol]); 
    correlation_maps[ipol]->Draw("colz"); 

    for (int i = 0; i < last.nPeaks[ipol]; i++) 
    {
      pads[ipol]->cd(1)->cd(2); 
      pueo::UCorrelator::gui::SummaryText * pt  = new gui::SummaryText(i, pol::pol_t(ipol), this, draw_filtered); 
      delete_list.push_back(pt); 
      pt->Draw(); 
    }


    pads[ipol]->cd(2)->Divide(last.nPeaks[ipol], interactive_deconvolved ? 5 : 3); 

    for (int i = 0; i < last.nPeaks[ipol]; i++) 
    {
      pads[ipol]->cd(2)->cd(i+1); 

      zoomed_correlation_maps[ipol][i]->SetTitle(TString::Format("Zoomed peak %d", i+1)); 
      zoomed_correlation_maps[ipol][i]->addFine(last.peak[ipol][i]); 
      zoomed_correlation_maps[ipol][i]->Draw("colz"); 

      pads[ipol]->cd(2)->cd(i+last.nPeaks[ipol]+1); 

      ((TGraph*) coherent[ipol][draw_filtered][i]->even())->SetTitle(TString::Format ( "Coherent (+ xpol) %d", i+1)); 
      coherent[ipol][draw_filtered][i]->drawEven("al"); 


      coherent_xpol[ipol][draw_filtered][i]->drawEven("lsame"); 


      pads[ipol]->cd(2)->cd(i+2*last.nPeaks[ipol]+1); 


      if (cfg->use_coherent_spectra) 
      {
        ((TGraph*) coherent[ipol][draw_filtered][i]->powerdB())->SetTitle(TString::Format ( "Power Coherent (+ xpol) %d", i+1)); 
        ((TGraph*) coherent[ipol][draw_filtered][i]->powerdB())->GetXaxis()->SetRangeUser(0,1.3); 
        coherent[ipol][draw_filtered][i]->drawPowerdB("al"); 
        ((TGraph*)coherent_xpol[ipol][draw_filtered][i]->powerdB())->SetLineColor(15); 
        ((TGraph*)coherent_xpol[ipol][draw_filtered][i]->powerdB())->SetName("powerdBXpol"); 
        coherent_xpol[ipol][draw_filtered][i]->drawPowerdB("lsame"); 
      }
      else
      {
        (((TGraph*)coherent_power[ipol][draw_filtered][i]))->SetTitle(TString::Format ( "Power Coherent (+ xpol) %d", i+1)); 
        ((TGraph*)coherent_power[ipol][draw_filtered][i])->Draw("al"); 
        ((TGraph*)coherent_power[ipol][draw_filtered][i])->GetXaxis()->SetRangeUser(0,1.3); 
        coherent_power_xpol[ipol][draw_filtered][i]->Draw("lsame"); 
        ((TGraph*)avg_spectra[ipol])->SetLineColor(2); 
        ((TGraph*)avg_spectra[ipol])->Draw("lsame"); 
      }



      /*
         TF1 * spectral_slope = new TF1(TString::Format("__slope_%d", i), "pol1",0.2,0.7); 
         spectral_slope->SetParameter(0, last.coherent[ipol][i].spectrumIntercept) ;
         spectral_slope->SetParameter(1, last.coherent[ipol][i].spectrumSlope) ;


         TGraphErrors *gbw = new TGraphErrors(EventSummary::peaksPerSpectrum); 
         gbw->SetTitle("Bandwidth Peaks"); 
         for (int bwpeak = 0; bwpeak < EventSummary::peaksPerSpectrum; bwpeak++) 
         {
         double bwf = last.coherent[ipol][i].peakFrequency[bwpeak]; 
         gbw->SetPoint(bwpeak, bwf, avg_spectra[ipol]->Eval(bwf)+ last.coherent[ipol][i].peakPower[bwpeak]); 
         gbw->SetPointError(bwpeak  , last.coherent[ipol][i].bandwidth[bwpeak]/2,0);
         }
         gbw->SetMarkerColor(4); 
         gbw->SetMarkerStyle(20); 
         gbw->Draw("psame"); 


         delete_list.push_back(spectral_slope); 
         delete_list.push_back(gbw);

*/  

      if (interactive_deconvolved)
      {
        pads[ipol]->cd(2)->cd(i+3*last.nPeaks[ipol]+1); 
        ((TGraph*) deconvolved[ipol][draw_filtered][i]->even())->SetTitle(TString::Format ( "Deconvolved (+ xpol) %d", i+1)); 
        deconvolved[ipol][draw_filtered][i]->drawEven("alp"); 
        if (interactive_xpol_deconvolved)
        {
          deconvolved_xpol[ipol][draw_filtered][i]->drawEven("lsame"); 
        }

        pads[ipol]->cd(2)->cd(i+4*last.nPeaks[ipol]+1); 

        if (cfg->use_coherent_spectra) 
        {

          ((TGraph*) deconvolved[ipol][draw_filtered][i]->powerdB())->SetTitle(TString::Format ( "Power Deconvolved (+ xpol) %d", i+1)); 
          ((TGraph*) deconvolved[ipol][draw_filtered][i]->powerdB())->GetXaxis()->SetRangeUser(0,1.3); 
          (deconvolved[ipol][draw_filtered][i])->drawPowerdB();; 
          if (interactive_xpol_deconvolved)
          {
            ((TGraph*) deconvolved_xpol[ipol][draw_filtered][i]->powerdB())->SetLineColor(15); 
            ((TGraph*) deconvolved_xpol[ipol][draw_filtered][i]->powerdB())->SetName("powerdBXpol"); 
            (deconvolved_xpol[ipol][draw_filtered][i])->drawPowerdB("lsame"); 
          }

        }
        else
        {
          (((TGraph*)deconvolved_power[ipol][draw_filtered][i]))->SetTitle(TString::Format ( "Power Deconvolved (+ xpol) %d", i+1)); 
          ((TGraph*)deconvolved_power[ipol][draw_filtered][i])->Draw();; 
          if (interactive_xpol_deconvolved)
          {
            ((TGraph*)deconvolved_power_xpol[ipol][draw_filtered][i])->Draw("lsame"); 
          }
        }


      }

    }
  }
}

void pueo::UCorrelator::Analyzer::fillFlags(const FilteredEvent * fae, EventSummary* summary, UsefulAttitude * pat) 
{
  EventSummary::EventFlags * flags = &summary->flags;
  flags->nadirFlag = true; // we should get rid of htis I guess? 

  if (cfg->fill_blast_fraction) 
    flags->blastFraction = TimeDependentAverageLoader::getPayloadBlastFraction(fae->getHeader()->triggerTime); 

  flags->meanPower[0] = fae->getAveragePower(); 
  flags->medianPower[0] = fae->getMedianPower(); 
  flags->meanPowerFiltered[0] = fae->getAveragePower(pol::kNotAPol, ring::kNotARing, true); 
  flags->medianPowerFiltered[0] = fae->getMedianPower(pol::kNotAPol, ring::kNotARing, true); 

  for (int ring = ring::kTopRing; ring <ring::kNotARing; ring++)
  {
    flags->meanPower[1+ring] = fae->getAveragePower(pol::kNotAPol, ring::ring_t(ring)); 
    flags->medianPower[1+ring] = fae->getMedianPower(pol::kNotAPol, ring::ring_t(ring)); 
    flags->meanPowerFiltered[1+ring] = fae->getAveragePower(pol::kNotAPol, ring::ring_t(ring), true); 
    flags->medianPowerFiltered[1+ring] = fae->getMedianPower(pol::kNotAPol, ring::ring_t(ring), true); 
  }


  flags->nSectorsWhereBottomExceedsTop = 0;
  for (int pol = pol::kHorizontal; pol <= pol::kVertical; pol++)
  {
    int n_above_this_pol; 
    fae->getMinMaxRatio(pol::pol_t(pol), &flags->maxBottomToTopRatio[pol], &flags->minBottomToTopRatio[pol], &flags->maxBottomToTopRatioSector[pol], &flags->minBottomToTopRatioSector[pol], ring::kBottomRing, ring::kTopRing,1, &n_above_this_pol); 
    flags->nSectorsWhereBottomExceedsTop += n_above_this_pol; 
  }


  if ( isLDB(fae->getHeader(), cfg))
  {
    flags->pulser = EventSummary::EventFlags::LDB; 
  }
  else if ( isWAISHPol(pat, fae->getHeader(), cfg) )
  {
    flags->pulser = EventSummary::EventFlags::WAIS; 
  }
  else if( isWAISVPol (pat, fae->getHeader(), cfg)){
    flags->pulser = EventSummary::EventFlags::WAIS_V;
  }
  else
  {
    flags->pulser = EventSummary::EventFlags::NONE; 
  }

  // more than 80 percent filterd out 
  flags->strongCWFlag = flags->medianPowerFiltered[0] / flags->medianPower[0] < 0.2; 

  if(version::get() == 4){
    //Blast event multi variant selector
    // 1)more than 30 channels that waveform rms > 100mV  2)p2p bottom to top ratio 3)meanPower bottom to top ratio
    flags->isPayloadBlast =
      summary->countChannelAboveThreshold(100)>30 ||
      flags->maxBottomToTopRatio[0] <0.9 || flags->maxBottomToTopRatio[0] >2.6 ||
      flags->maxBottomToTopRatio[1] <0.9 || flags->maxBottomToTopRatio[1] >2.6 ||
      flags->meanPower[3]/flags->meanPower[1] <0.2 || flags->meanPower[3]/flags->meanPower[1] >1.7; 
  }else{
    flags->isPayloadBlast =  
      (cfg->max_mean_power_filtered && flags->meanPowerFiltered[0] > cfg->max_mean_power_filtered) ||
      (cfg->max_median_power_filtered && flags->medianPowerFiltered[0] > cfg->max_median_power_filtered) ||
      (cfg->max_bottom_to_top_ratio && flags->maxBottomToTopRatio[0] > cfg->max_bottom_to_top_ratio) || 
      (cfg->max_bottom_to_top_ratio && flags->maxBottomToTopRatio[1] > cfg->max_bottom_to_top_ratio); 
  }

  flags->isVarner = false; 
  flags->isVarner2 = false; 

  flags->isGood = !flags->isVarner && !flags->isVarner2 && !flags->strongCWFlag; 

  flags->isStepFunction = 0; 

  //might as well use the power flags that Ben Strutt added
  const double bandsLowGHz[EventSummary::numBlastPowerBands] = {0.15-1e-10, 1.1-1e-10, 0};
  const double bandsHighGHz[EventSummary::numBlastPowerBands] = {0.26-1e-10, 999, 999};

  //reset values
  for(int band = 0; band < EventSummary::numBlastPowerBands; band++)
  {
    flags->middleOrBottomPower[band] = 0;
    flags->middleOrBottomAnt[band] = -1;
    flags->middleOrBottomPol[band] = 2;
    flags->topPower[band] = 0;
  }

  //FIXME study this code, and figure out if we want something like it in PUEO 
  /*
  for(Int_t polInd = 0; polInd < pol::kNotAPol; polInd++)
  {
    pol::pol_t pol = (pol::pol_t) polInd;
    for(UInt_t phi = 0; phi < NUM_PHI; phi++)
    {
      for(UInt_t ring = ring::kMiddleRing; ring < ring::kNotARing; ring++)
      {
        Int_t ant = ring*NUM_PHI + phi;

        const AnalysisWaveform* wf = fae->getRawGraph(ant,pol);
        const TGraphAligned* grPow = wf->power();
        const double df_GHz = grPow->GetX()[1] - grPow->GetX()[0];

        Double_t powThisAnt[EventSummary::numBlastPowerBands] = {0};
        for(int i = 0; i < grPow->GetN(); i++)
        {
          const double f_GHz = grPow->GetX()[i];
          for(int band = 0; band < EventSummary::numBlastPowerBands; band++)
          {
            if(f_GHz >= bandsLowGHz[band] && f_GHz < bandsHighGHz[band])
            {
              powThisAnt[band] += grPow->GetY()[i]*df_GHz;
            }
          }
        }

        for(int band = 0; band < EventSummary::numBlastPowerBands; band++)
        {
          if(powThisAnt[band] > flags->middleOrBottomPower[band])
          {
            flags->middleOrBottomPower[band] = powThisAnt[band];
            flags->middleOrBottomAnt[band] = ant;
            flags->middleOrBottomPol[band] = polInd;
          }
        }
      }
    }
  }

  for(int band = 0; band < EventSummary::numBlastPowerBands; band++)
  {
    int ant = flags->middleOrBottomAnt[band] % NUM_PHI;
    pol::pol_t pol = (pol::pol_t) flags->middleOrBottomPol[band];
    const AnalysisWaveform* wf = fae->getRawGraph(ant, pol);
    const TGraphAligned* grPow = wf->power();
    const double df_GHz = grPow->GetX()[1] - grPow->GetX()[0];

    for(int i = 0; i < grPow->GetN(); i++)
    {
      const double f_GHz = grPow->GetX()[i];
      for(int band = 0; band < EventSummary::numBlastPowerBands; band++)
      {
        if(f_GHz >= bandsLowGHz[band] && f_GHz < bandsHighGHz[band])
        {
          flags->topPower[band] += grPow->GetY()[i]*df_GHz;
        }
      }
    }
  }

  */
}

void pueo::UCorrelator::Analyzer::setExtraFilters(FilterStrategy* extra)
{
  wfcomb_filtered.setExtraFilters(extra);
  wfcomb_xpol_filtered.setExtraFilters(extra);
}

void pueo::UCorrelator::Analyzer::setExtraFiltersDeconvolved(FilterStrategy* extra)
{
  wfcomb_filtered.setExtraFiltersDeconvolved(extra);
  wfcomb_xpol_filtered.setExtraFiltersDeconvolved(extra);
}
