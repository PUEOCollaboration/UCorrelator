#ifndef _PUEO_UCORRELATOR_ANTENNA_POSITIONS_H
#define _PUEO_UCORRELATOR_ANTENNA_POSITIONS_H

#include <bitset>
#include "pueo/Conventions.h" 
#include "pueo/Version.h" 

namespace pueo 
{
  namespace UCorrelator
  {
    /** This class keeps the positions of the antennas, and has some related methods. It is a singleton. */ 
    class AntennaPositions
    {

      AntennaPositions(int v);

      public: 

        static const AntennaPositions * instance (int version = 0); 
        
        /** Retrieve an instance */ 
       /** Find closest N antennas to phi. Results put into closest, which should have sufficient room. Disallowed is a bitmap of antenna numbers that should be excluded. Returns number found (could be less than number requested if too many disallowed)*/
        int getClosestAntennas(double phi, int N, int * closest, std::bitset<k::NUM_HORNS> disallowed = 0, pol::pol_t pol = pol::kHorizontal) const; 

        /** antenna phi positions (degrees)*/
        double phiAnt[k::NUM_POLS][k::NUM_HORNS]; 

        /** antenna r positions (m) */
        double rAnt[k::NUM_POLS][k::NUM_HORNS]; 

        /** antenna z positions (m) */
        double zAnt[k::NUM_POLS][k::NUM_HORNS]; 

        double distance(int ant1, int ant2, pol::pol_t pol = pol::kHorizontal) const; 
        int  v; 

    }; 
  }
}


#endif
