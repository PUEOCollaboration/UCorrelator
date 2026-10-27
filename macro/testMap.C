#if !defined(__CLING__)
// For some fucked up reason, CLING doesn't like this header,
// so we'll only include it if we want to compile via g++ and hide it if running macro through CLING.
#include "pueo/Conventions.h"
#endif
#include "pueo/EventSummary.h" 
#include "pueo/Dataset.h" 
#include "pueo/ProbabilityMap.h"
#include "pueo/PointingResolutionModel.h"

#include "TF1.h"
#include "TFile.h"
#include "TTree.h"


void testMap() 
{
  int run=802;
  int event_number=197894;
  pueo::pol::pol_t pol = pueo::pol::kVertical;
  int peak = 0;

  TFile myFile(Form("/work/rootified/eventSummaries/%d.root", run));
  TTree * myTree = myFile.Get<TTree>("eventTree");
  myTree->BuildIndex("eventNumber");

  pueo::EventSummary * mySummary = nullptr;
  myTree->SetBranchAddress("eventSummary", &mySummary);
  myTree->GetEntryWithIndex(event_number); 

  pueo::nav::Attitude location;
  location = mySummary->location.att();

  TFile f("debug_data.root","RECREATE"); 

  StereographicGrid  * g = new StereographicGrid(1024,1024); 

  TF1 f_dtheta("ftheta", "[0] / x^[1]", 1, 50);
  TF1 f_dphi("fphi", "[0] / x^[1]", 1, 50);
  f_dtheta.SetParameter(0, 0.3936); 
  f_dtheta.SetParameter(1, 0.2102); 
  f_dphi.SetParameter(0, 1.065); 
  f_dphi.SetParameter(1, 0.2479); 
  pueo::UCorrelator::PointingResolutionParSNRModel * m1 = new pueo::UCorrelator::PointingResolutionParSNRModel(f_dtheta, f_dphi, true,true);
  pueo::UCorrelator::PointingResolutionModelPlusHeadingError*  m = new pueo::UCorrelator::PointingResolutionModelPlusHeadingError(20, m1); 

  Refraction::SphRay * ref = new Refraction::SphRay; 

  pueo::UCorrelator::ProbabilityMap::Params p; 
  p.refract = ref; 
  p.seg = g; 
  p.point = m; 
  p.collision_detection = false; 
  p.backwards_params.enhance_threshold = 1e-2; 
//  p.backwards_params.num_samples_per_bin=16;
  p.verbosity = 4;

  printf("THETA: %g\n", mySummary->peak[pol][peak].theta); 
  pueo::UCorrelator::ProbabilityMap *map = new pueo::UCorrelator::ProbabilityMap(&p); 
  map->add(mySummary, &location, pol, peak,1,&f); 
  map->segmentationScheme()->Draw("colz", map->getProbSums(true)); 
  //
  // f.cd(); 
  // map->Write("map"); 
  //
  // return map ; 

}

int main(){
  testMap();
  return 0;
}
