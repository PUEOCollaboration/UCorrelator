#include "pueo/AntennaPositions.h" 
#include "pueo/GeomTool.h"
#include "assert.h"
#include "FFTtools.h"
#include <map>
#include "TMutex.h"


#ifndef RAD2DEG
#define RAD2DEG (180/M_PI)
#endif




static const pueo::UCorrelator::AntennaPositions *instances[pueo::k::NUM_PUEO+1] = {0,0}; 


pueo::UCorrelator::AntennaPositions::AntennaPositions(int version)
{
  v = version; 
  printf("AntennaPositions(%d)\n",version); 

#ifdef MULTIVERSION_PUEO_ENABLED 
        const GeomTool * geom = &GeomTool::Instance(version); 
#else
        int old_ver = version::get(); 
        version::set(version); 
        const GeomTool * geom = &GeomTool::Instance(); 
        version::set(old_ver); 
#endif



  for (int pol = pol::kHorizontal; pol <= pol::kVertical; pol++)
  {
    for (int ant  = 0; ant < k::NUM_HORNS; ant++) 
    {
//      printf("%d %d\n",pol,ant); 
//
      phiAnt[pol][ant] = geom->getAntPhiPositionRelToAftFore(ant,(pol::pol_t)pol) * RAD2DEG; 
      rAnt[pol][ant] = geom->getAntR(ant,(pol::pol_t) pol); 
      zAnt[pol][ant] = geom->getAntZ(ant,(pol::pol_t) pol); 
    }
  }

}

int pueo::UCorrelator::AntennaPositions::getClosestAntennas(double phi, int N, int * closest, std::bitset<k::NUM_HORNS> disallowed , pol::pol_t pol) const
{

  assert(N < k::NUM_HORNS); 

  int pol_ind = pol == pol::kHorizontal ? 0 : 1; 
  std::multimap<double,int> dphis; 
  for (int i = 0; i < k::NUM_HORNS; i++) 
  {

    double dphi = disallowed.test(i) ? 360 : fabs(FFTtools::wrap(phi-phiAnt[pol_ind][i], 360, 0)); 
    dphis.insert(std::pair<double,int>(dphi,i)); 
  }

  int Nused = 0; 

//  printf("Closest antennas to %f: ",phi); 
  for (std::multimap<double,int>::const_iterator it = dphis.begin(); it !=dphis.end(); it++) 
  {
//    printf("  %d %f\n", (*it).second, (*it).first); 
    
    if ((*it).second >=360) break; 
    closest[Nused++] = (*it).second; 
    if (Nused == N) break; 
  }
//  printf("\n"); 
  return Nused; 
}



double pueo::UCorrelator::AntennaPositions::distance(int i, int j, pol::pol_t pol) const
{

#ifdef MULTIVERSION_PUEO_ENABLED 
  const GeomTool * geom = &GeomTool::Instance(v); 
#else
  int old_ver = version::get(); 
  version::set(v); 
  const GeomTool * geom = &GeomTool::Instance(); 
  version::set(old_ver); 
#endif

  double x0,y0,z0; 
  double x1,y1,z1; 
  geom->getAntXYZ(i,x0,y0,z0, pol); 
  geom->getAntXYZ(j,x1,y1,z1, pol); 


  double dx = x0-x1; 
  double dy = y0-y1; 
  double dz = z0-z1; 
  return sqrt( dx*dx + dy*dy + dz*dz); 

}

static TMutex instance_lock; 

const pueo::UCorrelator::AntennaPositions * pueo::UCorrelator::AntennaPositions::instance(int v)
{ 
  if (!v) v = version::get(); 

  const AntennaPositions * tmp = instances[v]; 
  __asm__ __volatile__ ("" ::: "memory"); //memory fence! 
  if (!tmp) 
  {
    instance_lock.Lock(); 
    tmp = instances[v]; 
    if (!tmp) 
    {
      tmp = new AntennaPositions(v); 
      __asm__ __volatile__ ("" ::: "memory");
      instances[v] = tmp; 
    }
    instance_lock.UnLock(); 
  }

  return instances[v];


}
     
