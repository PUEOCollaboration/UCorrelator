// @file testRaytracing.C
// @brief testMap.C currently doesn't work, so I'm testing AntarcticaRoot's ray tracer.
// Currently raytracing doesn't work either, it seems...

#include "TFile.h"
#include "TTree.h"

#if !defined(__CLING__)
#include "AntarcticaGeometry.h"
#else
#include "Rtypes.h"
R__LOAD_LIBRARY(libAntarcticaRoot);
#endif

#if !defined(__CLING__)
#include "pueo/UsefulAttitude.h"
#include "pueo/EventSummary.h"
#else
R__LOAD_LIBRARY(libpueoAnalysisFramework);
#endif


void testRaytracing(){
  int run = 802;
  int event_number = 197894;                      // This is a known Taylor Dome VPol pulser event.
  int peak = 0;                                   // Use the first peak of the skymap.
  pueo::pol::pol_t pol = pueo::pol::kVertical;

  TFile myFile(Form("/work/rootified/eventSummaries/%d.root", run));
  TTree * myTree = myFile.Get<TTree>("eventTree");
  myTree->BuildIndex("eventNumber");

  pueo::EventSummary * mySummary = nullptr;
  myTree->SetBranchAddress("eventSummary", &mySummary);
  myTree->GetEntryWithIndex(event_number); 

  // sign convention: DOWN IS POSITIVE (don't ask), so we have to change pueoTALON's sign convention
  double elevation_deg = mySummary->peak[pol][peak].theta * -1;
  double azimuth_deg = mySummary->peak[pol][peak].phi;
  printf("\033[1;32mRun %d Event %d\n\033[0m", run, event_number);
  printf(
    "Skymap Peak (positive elevation is downwards):\n"
    "\033[1;32m(phi, theta) = (%g, %g)\033[0m\n",
    azimuth_deg, elevation_deg
  ); 

  pueo::nav::Attitude _navatt = mySummary->location.att();
  pueo::UsefulAttitude payload_attitude(&_navatt);

  using Raytracer = PayloadParameters;
  Raytracer myRaytraceResult;

  int _status = Raytracer::findSourceOnContinent(
    elevation_deg, azimuth_deg,          // source direction from payload's perspective
    payload_attitude.asAntarcticCoord(), // payload lat lon alt
    &myRaytraceResult,
    payload_attitude.asPayloadAttitude() // payload heading etc
  );

  AntarcticCoord MaybeTaylorDome = myRaytraceResult.source.as(AntarcticCoord::WGS84);

  printf("\033[1;33mRaytracing result using (phi, theta) from above:\n");
  printf("lat: %g, lon: %g, alt: %g\n\033[0m", MaybeTaylorDome.x, MaybeTaylorDome.y, MaybeTaylorDome.z);
  printf("Actual Taylor Dome Location (copied from geographyUtil/locations.py):\n");
  printf("\033[1;32mlat: -77.78698333, lon: 158.6655333, alt: 2367\033[0m\n");

}
