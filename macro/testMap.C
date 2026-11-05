// For some fucked up reason, CLING doesn't like this header, 
// so we'll only include it if we are compiling this test, and hide it if running as a macro.
#if !defined(__CLING__)
#include "pueo/Conventions.h"
#else
R__LOAD_LIBRARY(libpueoEvent);
#endif
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "AntarcticaGeometry.h"           // AntarcticaRoot
#include "pueo/EventSummary.h"            // pueoAnalysisFramework
#include "pueo/ProbabilityMap.h"          // UCorrelator
#include "pueo/PointingResolutionModel.h"


namespace // Just some globals cuz I am lazy and don't wanna pass a bunch of shit around
{
const auto WHICH_POL  = pueo::pol::kVertical; // Build probability distribution based on the VPol skymap.
const int  WHICH_PEAK = 0 % 3; // Center the probability distribution around this peak of the skymap.

// const int RUN = 802, EVENT_NUMBER = 197894;
// const TString PLOT_NAME("maybe_TD_lmao.png");

const int RUN = 784, EVENT_NUMBER = 115752;
const TString PLOT_NAME("maybe_LDB_lmao.png");
}

namespace{  // forward declarations of unimportant functions, implemented later
pueo::EventSummary retrieve_one_event_summary();
pueo::UCorrelator::ProbabilityMap::Params specify_probability_map_settings();
}


void testMap() // This is the meat of this script; it makes and plots probability map
{
  TFile debug_file("debug_data.root","RECREATE"); 
  pueo::UCorrelator::ProbabilityMap::Params settings = specify_probability_map_settings();
  pueo::UCorrelator::ProbabilityMap probability_map(&settings); 

  // Why do we need a separate payload_attitude object if the info is already in the EventSummary?
  // IDK but this design is annoying af lol
  pueo::EventSummary mySummary = retrieve_one_event_summary();
  pueo::nav::Attitude payload_attitude = mySummary.location.att();

  // Use the info in the summary of this event (eg skymap peak location, snr, etc)
  // to build the probability map
  probability_map.add(&mySummary, &payload_attitude, WHICH_POL, WHICH_PEAK, 1, &debug_file);

  TCanvas c1("","", 1920, 1080);
  probability_map.segmentationScheme()->Draw("mapcolz", probability_map.getProbSums(true)); 
  c1.SaveAs(PLOT_NAME);

  debug_file.cd(); 
  probability_map.Write("map"); 
}


namespace
{
pueo::EventSummary retrieve_one_event_summary()
{
  // Nothing fancy here
  TFile myFile(Form("/work/rootified/eventSummaries/%d.root", RUN));
  TTree * myTree = myFile.Get<TTree>("eventTree");
  pueo::EventSummary * mySummary = nullptr;

  // Maybe we should just build the damn index when we make these summary files
  myTree->BuildIndex("eventNumber");
  myTree->SetBranchAddress("eventSummary", &mySummary);
  myTree->GetEntryWithIndex(EVENT_NUMBER); 

  // Making a stupid copy: not necessary but makes memory management easy :)
  pueo::EventSummary mySummaryCopy = *mySummary;

  // Manually adjust theta to be positive. This code is obviously very fagile lol.
  if (mySummaryCopy.peak[WHICH_POL][WHICH_PEAK].theta < 0) mySummaryCopy.peak[WHICH_POL][WHICH_PEAK].theta *= -1;
  printf("Source from payload's perspective (down is positive theta): (phi, theta) = (%g, %g)\n",
    mySummaryCopy.peak[WHICH_POL][WHICH_PEAK].phi, mySummaryCopy.peak[WHICH_POL][WHICH_PEAK].theta);

  return mySummaryCopy;
}


pueo::UCorrelator::ProbabilityMap::Params specify_probability_map_settings()
{
  pueo::UCorrelator::ProbabilityMap::Params probability_map_settings;
  probability_map_settings.collision_detection = false;
  probability_map_settings.backwards_params.enhance_threshold = 1e-2;
  probability_map_settings.backwards_params.num_samples_per_bin = 16;
  probability_map_settings.verbosity = 4;

  // The static keyword makes objects live forever.
  // Needed because the probability map class is designed to hold pointers to these external objects,
  // and these objects get destroyed without the static keyword upon function termination.
  // Alternatively we can have these objects in the global scope and drop the static keyword.
  static Refraction::SphRay ref;
  probability_map_settings.refract = &ref;

  // continent binning (ie "segmentation")
  int num_x_bins = 1024, num_y_bins = 1024;
  static StereographicGrid antarctica_map_grid(num_x_bins, num_y_bins);
  probability_map_settings.seg = &antarctica_map_grid;

  // Pointing resolution is a function of SNR, these two functions are just examples.
  int _snr_min = 1, _snr_max = 50;
  TF1 dtheta_vs_snr("ftheta", "[0] / x^[1]", _snr_min, _snr_max); // f(x) = A / x^B
  TF1 dphi_vs_snr  ("fphi",   "[0] / x^[1]", _snr_min, _snr_max);
  dtheta_vs_snr.SetParameter(0, 0.3936); // set parameter A
  dtheta_vs_snr.SetParameter(1, 0.2102); // set parameter B
  dphi_vs_snr.SetParameter(0, 1.065);
  dphi_vs_snr.SetParameter(1, 0.2479);
  static pueo::UCorrelator::PointingResolutionParSNRModel m1(dtheta_vs_snr, dphi_vs_snr, true, true);
  static pueo::UCorrelator::PointingResolutionModelPlusHeadingError m(20, &m1); 
  probability_map_settings.point = &m; 

  return probability_map_settings;
}
}


int main() // main function needed in case we want to compile the macro
{
  testMap();
  return 0;
}
