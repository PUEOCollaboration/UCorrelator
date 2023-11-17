#include "pueo/ShapeParameters.h" 
#include "FFTtools.h" 
#include "TGraph.h" 


double pueo::UCorrelator::shape::getRiseTime(const TGraph * g, double min, double max, int peak) 
{
  int ifirst = -1; 
  int ilast = -1; 
  if (peak < 0 || peak >= g->GetN()) 
  {
    FFTtools::getPeakVal(g,&peak); 
  }
  

  for (int i = peak; i >=0; i--) 
  {

    if (g->GetY()[i] <= min)
    {
      ifirst = i; 
      break; 
    }
    if (g->GetY()[i] <= max && ilast < 0)
    {
      ilast = i; 
    }
  }

  if (ifirst < 0 || ilast < 0) return -1; 

  return g->GetX()[ilast] - g->GetX()[ifirst]; 
}

double pueo::UCorrelator::shape::getFallTime(const TGraph * g, double min, double max, int peak) 
{
  int ifirst = -1; 
  int ilast = -1; 
 if (peak < 0 || peak >= g->GetN()) 
  {
    FFTtools::getPeakVal(g,&peak); 
  }
  

  for (int i = peak; i < g->GetN(); i++) 
  {

    if (g->GetY()[i] <= min)
    {
      ilast = i; 
      break; 
    }
    if (g->GetY()[i] <= max && ifirst < 0)
    {
      ifirst = i; 
    }
  }

  if (ifirst < 0 || ilast < 0) return -1; 

  return g->GetX()[ilast] - g->GetX()[ifirst]; 
}


double pueo::UCorrelator::shape::getWidth(const TGraph * g, double val, int * start, int * end, int peak)
{

  int ifirst = -1; 
  int ilast = -1; 

  if (peak < 0 || peak >= g->GetN()) 
  {
    FFTtools::getPeakVal(g,&peak); 
  }
  

  for (int i = peak; i < g->GetN(); i++) 
  {
    if (g->GetY()[i] <=val) 
    {
      ilast = i; 
      break;

    }
  }

  for (int i = peak; i >= 0; i--)
  {
    if (g->GetY()[i] <=val) 
    {
      ifirst = i; 
      break; 
    }
  }

  if (ifirst < 0) 
    ifirst = 0; 

  if (ilast < 0)
    ilast = g->GetN()-1; 

  double t0 = g->GetX()[ifirst]; 
  double t1 = g->GetX()[ilast]; 

  if (start) *start = ifirst; 
  if (end) *end = ilast; 

  return t1-t0; 

}
