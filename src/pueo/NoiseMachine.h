#ifndef PUEO_NOISE_MACHINE
#define PUEO_NOISE_MACHINE

#include "pueo/Conventions.h"
#include "pueo/NoiseSummary.h"
#include "pueo/EventSummary.h"
#include "pueo/FilteredEvent.h"
#include "pueo/Analyzer.h"
#include "pueo/UCorrelatorGUI.h"
#include "pueo/UCImageTools.h"

namespace pueo
{

  namespace UCorrelator
  {
    /*=====================
      A class to process and save information about the thermal environment*/
    class NoiseMachine
    {
     public:

      /* Constructor */
      NoiseMachine(const int length = 60);

      const int fifoLength; //one minute of noise averaging

      bool quiet;// = true; //whether to complain about things that I randomly decided upon

      //do you want to save the interferometric maps?  They are very large.  Also multiple ways to save them
      bool fillMap ;    //save the min bias maps as TH2D (~18kB per event)
      bool fillAvgMap; //save it as the average of fifoLength min bias maps (Won't work if fillMap==true)
      bool fillArray;  //save it as a double array (25% smaller, ~14kB per event)
      
      /* Reset */
      void zeroInternals();

      /* makes an averaged TProfile2D out of the histograms in the mapFifo */
      TProfile2D *getAvgMapNoiseProfile(pol::pol_t pol);


      /* updates all the fifos with the current event */
      void updateMachine(Analyzer *analyzer,FilteredEvent *filtered);


      /* Moves things into the summary */
      void fillNoiseSummary(NoiseSummary *noiseSummary); //crabcrabcrab

      /* Fill the mapHistoryVal value in EventSummary (eventSummary should be mostly filled already) */
      void fillEventSummary(EventSummary *eventSummary);
      //crab


      /* check if it is still all zeros or what */
      bool isJustInitialized() { return fJustInitialized; };


     private:

      //Makes sure the fifos start at zero
      bool fJustInitialized;

      /** Saves the last updated heading value */
      double lastHeading;
      
      //induvidual update functions which need to all be called at once so the fifos incriment properly
      /* for building up an enormous memory block of histograms from a bunch of events, then making an average 
                  also does the double array   */
      void updateAvgMapFifo(Analyzer *analyzer, FilteredEvent *filtered);
      /* for calculating rms of waveform from a minute average before event capture */
      void updateAvgRMSFifo(FilteredEvent *filtered);

      /* fills up the history for a single source, called by fillEventSummary */
      void setSourceMapHistoryVal(EventSummary::SourceHypothesis& source);

      //internals for time domain waveform rms fifo
      double *rmsFifo; //where the info is saved
      int rmsFifoPos; //where in the fifo the most recent write was
      bool rmsFifoFillFlag ; //whether you've completely filled the fifo once
      int rmsFifoIndex(int ant,int poli,int pos); //so I know where to write in the 1d array
      double rmsAvg[k::NUM_HORNS][k::NUM_POLS]; //a running count of the average per channel


      //internals for interferometric map fifo (probably enormous in memory so maybe make a flag)
      TH2 **mapFifo[k::NUM_POLS]; //where the info is saved
      int mapFifoPos;  //where in the fifo the most recent write was
      bool mapFifoFillFlag ; //whether you've completely filled the fifo once.
      int mapFifoIndex(); 


      //maybe will be more compressed
      static const int nPhi = 180; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic
      static const int nTheta = 100; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic
      //where the quickly calculated double array is stored.  should be maxDirections*nPhi*nTheta*NUM_POLS long
      double *rollingMapAvg;
      int rollingMapIndex(int poli,int iPhi,int iTheta); 

      ClassDefNV(NoiseMachine, 0); 

    };
/*------------------*/
  }
}


#endif
