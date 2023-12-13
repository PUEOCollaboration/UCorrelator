#ifndef PUEO_UCORRELATOR_CORRELATOR_H
#define PUEO_UCORRELATOR_CORRELATOR_H


#include "pueo/Conventions.h"
#include "TH2.h"
#include "pueo/AnalysisConfig.h"
#include <bitset>

namespace pueo 
{
  class FilteredEvent; 
  class AnalysisWaveform; 


  namespace UCorrelator
  {

    class TrigCache; 

    //This is just to hide the OpenMP stuff from cling
    class CorrelatorLocks; 

    /** This creates the inteferometric map for an ANITA event */ 
    class Correlator
    {
      public:
        /** Create a correlator with the following options for the rough map */ 
        Correlator(int nphi, double phimin, double phimax, int ntheta, double theta_lowest, double theta_highest, bool use_bin_center = false, bool scale_by_cos_theta = false, double baseline_weight = 0, double gain_sigma = 0); 

        /** Compute the rough correlation map for the event and pol */ 
        void compute(const FilteredEvent * event, pol::pol_t pol); 

        /** Get the rough correlation map */ 
        const TH2D * getHist() const { return hist; } 
   
        /** Get the rough correlation map normalization */ 
        const TH2D * getNorm() const { return norm; } 
        

        /** Compute a zoomed in map around phi and theta. nphi,dphi,ntheta,dtheta. If nant is non-zero, only the nearest nant antennas are used. You can use useme to avoid allocating a new TH2.  */ 
        TH2D* computeZoomed(double phi, double theta, int nphi, double dphi,  int ntheta, double dtheta, int nant = 0, TH2D * useme = 0); 

        /** Disable the antennas given by the bitmap */ 
        void setDisallowedAntennas(std::bitset<k::NUM_ANTS> disallowed) { disallowed_antennas = disallowed; } 

        /** Enable only the antennas given by the bitmap */ 
        void setAllowedAntennas(std::bitset<k::NUM_ANTS> allowed) { disallowed_antennas = ~allowed; } 

        /** An antenna only contributes to an angle if it's within max_phi  of it */
        void setMaxAntennaMaxPhiDistance(double max_ant_phi) { max_phi = max_ant_phi;  max_phi2 = max_phi * max_phi; } 

        /** Use the group delay in computing the delay */ 
        void setGroupDelayFlag(bool flag) { groupDelayFlag = flag; } 


        /** Get the correlation between two antennas */
        const AnalysisWaveform * getCorrelationGraph(int ant1, int ant2) { return getCorrelation(ant1,ant2); }

        /** Set the supersampling factor */ 
        void setPadFactor(int pad) { pad_factor = pad; } 

        //Added for BinnedAnalysis - JCF 9/27/2021
        /** Disable the antennas given by the bitmap, but leave all others alone */    // sammy 2016/10/27
        void disallowAntennas(std::bitset<k::NUM_ANTS> disallowed) { disallowed_antennas |= disallowed; }

        /** Enable the antennas given by the bitmap, but leave all others alone */    // sammy 2016/10/27
        void allowAntennas(std::bitset<k::NUM_ANTS> allowed) { disallowed_antennas &= ~allowed; }

        /** Set the normalization      added sammy */
        void setNormalization(AnalysisConfig::NormalizationOption_t n) { normalization_option = n; }
        const char * getNormalizationString() {return AnalysisConfig::getNormalizationString(normalization_option);}
        //----------------------------------------------------------------------

        //End BinnedAnalysis Additions

        /** Debugging method to dump out some info to a file */ 
        void dumpDeltaTs(const char * file) const; 
        virtual ~Correlator(); 

      private: 
        AnalysisWaveform* padded_waveforms[k::NUM_ANTS]; 
        AnalysisWaveform* correlations[k::NUM_ANTS][k::NUM_ANTS]; 

        TH2D *hist; 
        TH2D *norm; 

        std::vector<TH2D*> hists; 
        std::vector<TH2D*> norms; 

        TrigCache * trigcache[k::NUM_PUEO+1]; 
        double rms[k::NUM_ANTS]; 

        double max_phi, max_phi2;
        std::bitset<k::NUM_ANTS> disallowed_antennas;
        int pad_factor;
        const FilteredEvent * ev; 
        pol::pol_t pol; 
        bool groupDelayFlag; 
        bool use_bin_center; 
        bool scale_cos_theta; 
        double baselineWeight;
        double gainSigma; 

        AnalysisWaveform * getCorrelation(int ant1, int ant2); 
        void doAntennas(int ant1, int ant2, TH2D ** hist, TH2D ** norm, const UCorrelator::TrigCache * tc, const double * center_point  = 0); 
        void reset(); 

        CorrelatorLocks * locks; 

        //Added for BinnedAnalysis - JCF 9/27/2021
        // added sammy -----------------------------------------------
        AnalysisConfig::NormalizationOption_t normalization_option;
        //------------------------------------------------------------
        //End BinnedAnalysis addition.
    }; 



  }
}



#endif
