#ifndef PUEO__UCORRELATOR_SPECTRUM_PARAMETERS_H
#define PUEO__UCORRELATOR_SPECTRUM_PARAMETERS_H

class TGraph; 

#include "pueo/EventSummary.h" 

namespace pueo {
namespace UCorrelator
{
  class AnalysisConfig; 
  namespace spectrum
  {
    void fillSpectrumParameters( const TGraph * spectrum, const TGraph * baseline, EventSummary::WaveformInfo * winfo, const AnalysisConfig * config); 
  }
}
}



#endif
