#ifndef _PUEO_UCORRELATOR_FLAGS_H
#define _PUEO_UCORRELATOR_FLAGS_H

/** \file This file is full of a bunch of methods that check various things about events */ 

// Where a bunch is apparently defined as one 

#include "TObject.h" // for RTypes
#include <bitset> 
#include "pueo/Conventions.h" 

namespace pueo 
{
class UsefulEvent; 

namespace UCorrelator
{
  namespace flags 
  {
    /** Checks for missing antennas and marks the bitmasks */ 
    int checkEmpty(const UsefulEvent *ev,  std::bitset<k::NUM_HORNS> *hempty, std::bitset<k::NUM_HORNS> *vempty); 


  }
}
}




#endif
