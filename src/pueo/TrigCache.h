#ifndef PUEO__UCORRELATOR_TRIGCACHE_H
#define PUEO__UCORRELATOR_TRIGCACHE_H

#include "pueo/AntennaPositions.h"
#include <cmath>

#ifndef DEG2RAD
#define DEG2RAD  (M_PI / 180.)
#endif

namespace pueo {
namespace UCorrelator
{

/** This class just defines a trig function cache for the correlator nothing interesting! 
 *
 *   Angles all in degrees. 
 **/ 
  struct TrigCache
  {
    TrigCache(int n_phi, double dphi, double phi_start, int ntheta, double dtheta, double theta_start, const AntennaPositions * aps, bool use_bin_center= false, int nant2use = 0, const int * ants = 0) 
      : nphi(n_phi) ,ap(aps)
    {

      phi = new double[nphi]; 
      int num_ants = nant2use == 0 ? k::NUM_HORNS : nant2use; 
      cos_phi = new double[nphi * k::NUM_HORNS * 2];  //allocate enough space even if unused... the non-wanted ants will be full of junk
      theta = new double[ntheta]; 
      tan_theta = new double[ntheta]; 
      cos_theta = new double[ntheta]; 

      for (int i = 0; i < nphi; i++) 
      {
        double current_phi = phi_start + dphi * i; 
        if (use_bin_center) current_phi += 0.5 * dphi; 

        phi[i] = current_phi; 
        for (int j = 0; j < num_ants; j++)
        {
          int ant =  nant2use ?  ants[j] : j; 
          cos_phi[2 * (i + ant * nphi)] = cos( DEG2RAD *  ( current_phi - ap->phiAnt[0][ant] )); 
          cos_phi[2 * (i + ant * nphi) + 1] = cos( DEG2RAD *  ( current_phi - ap->phiAnt[1][ant] )); 
        }
      }

      for (int i = 0; i < ntheta; i++) 
      {
        double current_theta = theta_start + dtheta * i; 
        if (use_bin_center) current_theta += 0.5 * dphi; 
        theta[i] = -current_theta; //lol
        cos_theta[i] = cos(theta[i] * DEG2RAD);  //naturally, we reverse the sign... 
        tan_theta[i] = tan(theta[i] * DEG2RAD);  //naturally, we reverse the sign... 
      }


    }; 


    ~TrigCache()
    {

      delete [] phi; 
      delete [] cos_phi; 
      delete [] theta; 
      delete [] cos_theta; 
      delete [] tan_theta; 

    }


    int nphi; 

    double * phi; //degrees!
    double * cos_phi; 
    //lookup since depends on antenna / polarization
    // but caerful, if you restricted antennas and ask for one you didn't want, you get junk!
    double cosPhi(int phibin, int ant, int pol) { return cos_phi[ 2 * (phibin + ant * nphi) + pol]; }
    double * theta; //degrees! 
    double * cos_theta; 
    double * tan_theta;

    const AntennaPositions * ap; 

  }; 



} 
}


#endif
