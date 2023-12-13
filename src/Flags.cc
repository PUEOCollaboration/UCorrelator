#include "pueo/UCFlags.h" 
#include "pueo/UsefulEvent.h"
#include "pueo/GeomTool.h" 




int pueo::UCorrelator::flags::checkEmpty(const UsefulEvent *ev, std::bitset<k::NUM_ANTS> * save_h, std::bitset<k::NUM_ANTS> * save_v)
{
  

  int nmissing = 0; 
  for (int i = 0; i < k::NUM_HORNS; i++) 
  {
      int hindex = GeomTool::getChanIndexFromAntPol(i,pol::kHorizontal); 
      int vindex = GeomTool::getChanIndexFromAntPol(i,pol::kVertical); 

      const double *yh = &ev->volts[hindex][0]; 
      const double *yv = &ev->volts[vindex][0]; 

      bool hbad = true; 
      for (int j = 0; j < ev->volts[hindex].size(); j++)
      {
        if (yh[j]) 
        {
          hbad = false; 
          break; 
        }
      }

      bool vbad = true; 
      for (int j = 0; j < ev->volts[vindex].size(); j++)
      {
        if (yv[j])
        {
          vbad = false; 
          break; 
        }
      }


      if (hbad && save_h) save_h->set(i);
      if (vbad && save_v) save_v->set(i); 


      if (hbad) nmissing++; 
      if (vbad) nmissing++; 

  }


  return nmissing; 
}
