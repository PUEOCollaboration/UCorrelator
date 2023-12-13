#include "pueo/WaveformCombiner.h" 
#include "pueo/AntennaPositions.h" 
#include <vector> 
#include "pueo/BasicFilters.h" 
#include "pueo/FilteredEvent.h"
#include "pueo/Version.h" 
#include "FFTtools.h"
#include "pueo/ResponseManager.h" 
#include "pueo/SystemResponse.h"
#include "pueo/DeltaT.h"
#include <dirent.h>
#include <sys/types.h>
#include <stdio.h>

  pueo::UCorrelator::WaveformCombiner::WaveformCombiner(int nantennas, int npad, bool useUnfiltered, bool deconvolve, const ResponseManager * response)
: coherent(2048), deconvolved(2048)
{
  setNAntennas(nantennas); 
  setNPad(npad); 
  setDeconvolve(deconvolve); 
  setUseUnfiltered(useUnfiltered); 
  setGroupDelayFlag(true); 
  setRTimeShiftFlag(true); 
  use_raw= useUnfiltered; 
  setResponseManager(response); 
  setBottomFirst(false);
  setDelayToCenter(false);
  extra_filters = 0;
  extra_filters_deconvolved = 0;
  enable_r_time_shift = true;
  enable_simulation_time_shift = false;
  max_vpp = 0;
  check_vpp = 0;
}





pueo::UCorrelator::WaveformCombiner::~WaveformCombiner()
{
}

const pueo::AnalysisWaveform * pueo::UCorrelator::WaveformCombiner::getDeconvolved() const 
{

  if (!do_deconvolution)
  {
    //    fprintf(stderr,"WARNING! Deconvolution has not been enabled!\n"); 
    return 0; 
  }

  return &deconvolved; 
}




static void scaleGraph(TGraph * g, double C)
{

  int max = g->GetN(); 

  for (int i = 0; i < max; i++) 
  {
    g->GetY()[i] *= C; 
  }
}

static void addToGraph(TGraph * g, const TGraph * h)
{
  int max = TMath::Min(g->GetN(), h->GetN()); 
  for (int i = 0; i < max; i++)
  {
    g->GetY()[i] += h->GetY()[i]; 
  }
}


const pueo::AnalysisWaveform * pueo::UCorrelator::WaveformCombiner::wf (const FilteredEvent * event, int ant, pol::pol_t pol) 
{
  return use_raw ? event->getRawGraph(ant, pol) : event->getFilteredGraph(ant,pol); 


}



void pueo::UCorrelator::WaveformCombiner::combine(double phi, double theta, const FilteredEvent * event, pol::pol_t pol, std::bitset<k::NUM_ANTS> disallowed, double t0, double t1, double * avg_of_peaks, bool use_hilbert)
{

  std::vector<AnalysisWaveform> padded(nant);
  const int numDeconv = do_deconvolution ? nant : 0;
  std::vector<AnalysisWaveform> deconv(numDeconv);

  const pueo::UCorrelator::AntennaPositions * antpos = pueo::UCorrelator::AntennaPositions::instance(); 
  nant = antpos->getClosestAntennas(phi, nant, &antennas[0], disallowed); 
  double delays[nant]; 

  if (bottom_first) {
    if (antennas[0] < 32) {
      int toMove = antennas[0];
      if (antennas[1] > 32) {
        antennas[0] = antennas[1];
        antennas[1] = toMove; }
      else {
        antennas[0] = antennas[2];
        antennas[2] = toMove; }
    }
  }

  double sum_peak = 0; 
  for (int i = 0; i < nant; i++) 
  {
    //ensure transform already calculated so we don't have to repeat when deconvolving
    (void) wf(event,antennas[i],pol)->freq(); 
    int ipol = (pol == pol::kVertical) ? 1 : 0;

    padded[i].~AnalysisWaveform(); 
    new (&padded[i]) AnalysisWaveform(*wf(event,antennas[i],pol));
    if(extra_filters)
    {
      for(int j = 0; j < (int) extra_filters->nOperations(); j++)
      {
        FilterOperation * filterOp = (FilterOperation*) extra_filters->getOperation(j);
        filterOp->processOne(&padded[i], event->getHeader(), antennas[i], ipol);
        //delete filterOp;
      }
    }

    if (i == 0)
    {
      const TGraphAligned *p = wf(event,antennas[0],pol)->power(); 
      coherent_avg_spectrum.adopt(p); 
    }
    else
    {
      addToGraph(&coherent_avg_spectrum,wf(event,antennas[0],pol)->power()); 
    }

    

    if (avg_of_peaks) 
    {
      sum_peak += (use_hilbert ? padded[i].hilbertEnvelope() : padded[i].even())->peakVal(); 
    } 


    padded[i].padFreq(npad);

    if (do_deconvolution)
    {
      deconv[i].~AnalysisWaveform(); 
      new (&deconv[i]) AnalysisWaveform(*wf(event,antennas[i],pol));
      if(extra_filters_deconvolved)
      {
        for(int j = 0; j < (int) extra_filters_deconvolved->nOperations(); j++)
        {
          FilterOperation * filterOp = (FilterOperation*) extra_filters_deconvolved->getOperation(j);
          filterOp->processOne(&deconv[i], event->getHeader(), antennas[i], ipol);
          //delete filterOp;
        }
      }
      responses->response(pol,antennas[i])->deconvolveInPlace(&deconv[i], responses->getDeconvolutionMethod(), theta ); //TODO add angle  
      if (i == 0)
      {
        const TGraphAligned *p =deconv[i].power(); 
        deconvolved_avg_spectrum.adopt(p); 
      }
      else
      {
        addToGraph(&deconvolved_avg_spectrum, deconv[i].power()); 
      }
      deconv[i].padFreq(npad); 
    }

    if (delay_to_center) delays[i] = getDeltaTtoCenter(antennas[i], phi, theta, pol, enable_group_delay, enable_r_time_shift, enable_simulation_time_shift);
    else delays[i] = i == 0 ? 0 : getDeltaT(antennas[i], antennas[0], phi, theta, pol, enable_group_delay); 

  }

  max_vpp = combineWaveforms(nant, &padded[0], delays, 0, &coherent, t0, t1, check_vpp); 
  scaleGraph(&coherent_avg_spectrum, 1./nant); 
  if (do_deconvolution)
  {
    combineWaveforms(nant, &deconv[0], delays, 0, &deconvolved, t0, t1, 0); 
    scaleGraph(&deconvolved_avg_spectrum, 1./nant); 
  }

  if (avg_of_peaks) 
  {
    *avg_of_peaks = sum_peak / nant; 
  }
}


double pueo::UCorrelator::WaveformCombiner::combineWaveforms(int nwf, const AnalysisWaveform * wf, const double * delays, const double * scales, AnalysisWaveform * out, double min, double max, int checkvpp)
{
  // we want to make the waveform as big as the entire span of all the waveforms to combine 

  double dt = wf[0].deltaT(); 

  int N = ceil((max-min)/dt); 
  double maxvpp = 0;

  if (!out) out = new AnalysisWaveform(N); 
  TGraph * gupdate = out->updateEven(); 
  gupdate->Set(N); 

  for (int i = 0; i < N; i++) 
  {
    gupdate->GetX()[i] = min + i * dt; 

    double val = 0; 
    for (int j = 0; j < nwf; j++) 
    {
      double scale = scales ? scales[j] : 1; 
      double tempvpp = 0;
      val += wf[j].evalEven(gupdate->GetX()[i] + delays[j])/nwf*scale; 
      if(i == 0)
      {
        if(checkvpp == 1)
        {
          tempvpp = wf[j].even()->pk2pk();
          if(maxvpp < tempvpp) maxvpp = tempvpp;
        }
        if(checkvpp == 2)
        {
          tempvpp = wf[j].even()->pk2pk()/double(nwf); 
          maxvpp += tempvpp;
        }
      }
    }

    //    printf("%f %f\n", gupdate->GetX()[i], val); 
    gupdate->GetY()[i] = val; 
  }

  return maxvpp; 
}

void pueo::UCorrelator::WaveformCombiner::setExtraFilters(FilterStrategy* extra)
{
  if(extra_filters) delete extra_filters;
  extra_filters = extra;
}

void pueo::UCorrelator::WaveformCombiner::setExtraFiltersDeconvolved(FilterStrategy* extra)
{
  if(extra_filters_deconvolved) delete extra_filters_deconvolved;
  extra_filters_deconvolved = extra;
}

