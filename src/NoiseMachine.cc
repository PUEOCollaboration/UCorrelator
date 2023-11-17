#include "pueo/NoiseMachine.h"


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
pueo::UCorrelator::NoiseMachine::NoiseMachine(const int length)
    : fifoLength(length), quiet(true)
{
  
  if (fifoLength < 2) {
    std::cout << "Warning from constructor of pueo::UCorrelator::NoiseMachine(int length): ";
    std::cout << "You picked an average less than 2, this is a bad idea, and values will likely be wrong" << std::endl;
  }

  rmsFifo = (double*)malloc(k::NUM_HORNS*k::NUM_POLS*fifoLength*sizeof(double));
  rollingMapAvg = (double*)malloc(k::NUM_POLS*nPhi*nTheta*sizeof(double));


  for (int poli=0; poli<k::NUM_POLS; poli++) {
    mapFifo[poli] = (TH2**)malloc(fifoLength*sizeof(TH2*));
    for (int fifoPos=0; fifoPos<fifoLength; fifoPos++) {
      mapFifo[poli][fifoPos] = NULL;
    }
  }
  
  zeroInternals();
}

//Resetting to initial state
void pueo::UCorrelator::NoiseMachine::zeroInternals() {

  fJustInitialized = true;

  //reset save flags
  fillArray = false;
  fillMap = false;
  fillAvgMap = false;

  //heading
  lastHeading = -999;

  //reset rms fifo
  rmsFifoPos = 0;
  rmsFifoFillFlag = false;
  memset(rmsFifo,0,k::NUM_HORNS*k::NUM_POLS*fifoLength*sizeof(double)); 
  memset(rmsAvg,0,k::NUM_HORNS*k::NUM_POLS*sizeof(double));

  //reset map histogram fifo
  for (int poli=0; poli<k::NUM_POLS; poli++) {
    for (int loc=0; loc<fifoLength; loc++) {
      if (mapFifo[poli][loc] != NULL) {
	delete mapFifo[poli][loc];
	mapFifo[poli][loc] = NULL;
      }
    }
  }
  mapFifoPos = 0;
  mapFifoFillFlag = false;

  //reset map double array fifo
  memset(rollingMapAvg,0,k::NUM_POLS*nPhi*nTheta*sizeof(double));

}




//so that I always am refering to the same index
int pueo::UCorrelator::NoiseMachine::rmsFifoIndex(int ant, int poli, int pos) {
  int index = ant*(k::NUM_POLS*fifoLength) + poli*(fifoLength) + pos;
  int maximumIndex = k::NUM_HORNS*k::NUM_POLS*fifoLength;
  if (index >= maximumIndex) {
    std::cout << "rmsFifoIndex: " << index << " > " << maximumIndex << std::endl;
    return (maximumIndex - 1); }
  
  return index;
}




//so that I always am refering to the same index
int pueo::UCorrelator::NoiseMachine::rollingMapIndex(int poli, int iPhi, int iTheta) {
  int index = poli*(nPhi*nTheta) + iPhi*(nTheta) + iTheta;
  int maximumIndex = nPhi*nTheta*k::NUM_POLS;
  if (index >= maximumIndex) {
    std::cout << "rollingMapIndex(): uhh you asked for an index that is bigger than the array ";
    std::cout << "(" << index << " > " << maximumIndex << ")" << std::endl;
    std::cout << "poli:" << poli << " iPhi:"  << iPhi << " iTheta:" << iTheta << std::endl;
    return (maximumIndex - 1);
  }
  return index;
}



/*===============
  Fill everything */

void pueo::UCorrelator::NoiseMachine::updateMachine(UCorrelator::Analyzer *analyzer,FilteredEvent *filtered) {
  updateAvgRMSFifo(filtered);
  updateAvgMapFifo(analyzer,filtered);

  fJustInitialized = false;
  return;
}









/*=======================
  Updates Waveform RMS stuff */
void pueo::UCorrelator::NoiseMachine::updateAvgRMSFifo(FilteredEvent *filtered) {
  
  if (!fJustInitialized) rmsFifoPos++;
  if (rmsFifoPos >= fifoLength) {
    rmsFifoPos = 0;
    rmsFifoFillFlag = true;
  }


  for (int iant=0; iant<k::NUM_HORNS; iant++) {
      for (int poli=0; (pol::pol_t)poli != pol::kNotAPol; poli++) {

        pol::pol_t pol = (pol::pol_t)poli;

        //subtract fifo position that is expiring (if there is one)
        if (rmsFifoFillFlag) {
          double valueSub = rmsFifo[rmsFifoIndex(iant,poli,rmsFifoPos)];
          rmsAvg[iant][poli] -= valueSub;
        }	

        //get the rms
        const TGraphAligned *currWave = filtered->getFilteredGraph(iant,pol)->even();
        double value = pow(currWave->GetRMS(2),2)/fifoLength;

        //add the new fifo position to the running average
        rmsAvg[iant][poli] += value;

        //save it for the fifo for when you gotta subtract it later
        rmsFifo[rmsFifoIndex(iant,poli,rmsFifoPos)] = value;
    }
  }

  return;
}

/*---------------------*/
	

/*======================
  Get the correlation maps from the analyzer and update everything internally 
  updates:  mapFifo - the collection of recent map histograms
            rollingMapAvg - the double array that is in essence a histogram I can also subtract from

*/
void pueo::UCorrelator::NoiseMachine::updateAvgMapFifo(UCorrelator::Analyzer *analyzer, FilteredEvent *filtered) {
  std::stringstream name;
  
  //I gotta adjust it for heading too, so pull it out of filtered (do it like this for later when it gets integrated)
  double heading = filtered->getGPS()->heading;
  lastHeading = heading;

  //update position unless you if you haven't even written once though, then stay at zero
  if (!fJustInitialized) mapFifoPos++;

  //wrap if you're at the end of the fifo
  if (mapFifoPos >= fifoLength) {
    mapFifoPos = 0;
    mapFifoFillFlag = true;
    //    std::cout << "Buffer has been filled!" << std::endl;
  }

  //do this for all polarization maps
  for (int poli=0; poli<k::NUM_POLS; poli++) {

    
    //if you've filled up the whole fifo and are now scanning
    if (mapFifoFillFlag) {

    //subtract fifo position that is expiring from the averages(if there is one)
      for (Int_t iPhi = 0; iPhi < nPhi; iPhi++) {
	for (Int_t iTheta = 0; iTheta < nTheta; iTheta++) {     
	  double valueSub = mapFifo[poli][mapFifoPos]->GetBinContent(iPhi+1,iTheta+1);
	  rollingMapAvg[rollingMapIndex(poli,iPhi,iTheta)] -= valueSub;
	  if (!quiet) {
	    if (poli==0 && iPhi==1 && iTheta==61 ) {
	      std::cout << "subtracting " << valueSub << std::endl; } }
	}
      }

      //and delete the old histogram if it is there
      if (mapFifo[poli][mapFifoPos] != NULL) {
	delete mapFifo[poli][mapFifoPos];
	mapFifo[poli][mapFifoPos] = NULL;
      }
    
    } //if (mapFifoFillFlag)


    //syntactically weirdly copy the new map out of UCorrelator    
    const UCorrelator::gui::Map *currMap = analyzer->getCorrelationMap((pol::pol_t)poli);
    
    //rotate it and put the rotated map into the fifo for later
    mapFifo[poli][mapFifoPos] = UCorrelator::rotateHistogram(currMap,heading);
    TH2 *tempHist = mapFifo[poli][mapFifoPos]; //because this is easier to read
    name.str("");
    name << "mapFifo[" << poli << "][" << mapFifoPos << "]";
    tempHist->SetName(name.str().c_str());
       

    //and add the new addition to the fifo
    for (Int_t iPhi = 0; iPhi < nPhi; iPhi++) {
      for (Int_t iTheta = 0; iTheta < nTheta; iTheta++) {     
	double valueAdd = tempHist->GetBinContent(iPhi+1,iTheta+1);
	rollingMapAvg[rollingMapIndex(poli,iPhi,iTheta)] += valueAdd;

	if (!quiet) {
	  if (poli==0 && iPhi==1 && iTheta==61 ) {
	    std::cout << "adding " << valueAdd;
	    std::cout << " => " << rollingMapAvg[rollingMapIndex(poli,iPhi,iTheta)] << std::endl; } }
  

      }//end theta
    }//end phi

  }//end of pol
  
  
  return;
  
}

/*------------*/

TProfile2D* pueo::UCorrelator::NoiseMachine::getAvgMapNoiseProfile(pol::pol_t pol) {

  int nBinX = mapFifo[(int)pol][0]->GetNbinsX();
  int nBinY = mapFifo[(int)pol][0]->GetNbinsY();
  double xMin  = mapFifo[(int)pol][0]->GetXaxis()->GetBinLowEdge(1);
  double yMin  = mapFifo[(int)pol][0]->GetYaxis()->GetBinLowEdge(1);
  double xMax  = mapFifo[(int)pol][0]->GetXaxis()->GetBinUpEdge(nBinX);
  double yMax  = mapFifo[(int)pol][0]->GetYaxis()->GetBinUpEdge(nBinY);

  std::stringstream name;
  name.str("");
  name << "mapAvg_" << pol;
  
  TProfile2D *outProfile  = new TProfile2D(name.str().c_str(),"Average Interferometric Map", nBinX, xMin, xMax, nBinY, yMin, yMax);

  int lengthToDo = fifoLength;
  if (mapFifoFillFlag == false) lengthToDo = mapFifoPos;

  //stolen from Rene
  for (Int_t ix = 1; ix <= nBinX; ix++) {
    Double_t x = outProfile->GetXaxis()->GetBinCenter(ix);
    for (Int_t iy = 1; iy <= nBinY; iy++) {     
      Double_t y = outProfile->GetYaxis()->GetBinCenter(iy);
      for (Int_t i=0; i<lengthToDo; i++) {
	outProfile->Fill(x,y,mapFifo[(int)pol][i]->GetBinContent(ix,iy),1);
      }
    }
  }

  return outProfile;
}




/*===========================
  Copy the relevant things into the noise summary */
void pueo::UCorrelator::NoiseMachine::fillNoiseSummary(NoiseSummary *noiseSummary) {
  if (!rmsFifoFillFlag && !quiet) {
    std::cout << "WARNING in pueo::UCorrelator::NoiseMachine::fillNoiseSummary(): Fifo hasn't been filled entirely yet, ";
    std::cout << " only gotten " << rmsFifoPos << " out of " << fifoLength << " values" << std::endl;
  }
  
  //update flags and length
  noiseSummary->fifoLength = fifoLength;
  noiseSummary->mapFifoFillFlag = mapFifoFillFlag;
  noiseSummary->rmsFifoFillFlag = rmsFifoFillFlag;

  //update RMS stuff (always do this, it compresses well and is important)
  for (int ant=0; ant<k::NUM_HORNS; ant++) {
      for (int poli=0; poli<k::NUM_POLS; poli++) {
        noiseSummary->avgRMSNoise[ant][poli] = rmsAvg[ant][poli];
    }
  }
  

  //fill a map, but only if one of the flags is on, fillMap is used if both are selected
  if (fillAvgMap || fillMap) {
    for (int poli=0; (pol::pol_t)poli != pol::kNotAPol; poli++) {
      pol::pol_t pol = (pol::pol_t)poli;
      if (noiseSummary->avgMapProf[poli] != NULL) {
	delete noiseSummary->avgMapProf[poli];
	noiseSummary->avgMapProf[poli] = NULL;
      }
      if (!quiet) std::cout << "filling map: " << poli <<  std::endl;
      if (fillMap) noiseSummary->avgMapProf[poli] = (TH2D*)mapFifo[poli][mapFifoPos]->Clone(); //non-average is default
      else noiseSummary->avgMapProf[poli] = getAvgMapNoiseProfile(pol);
    }
  }
  
  return;
  
}
 
/*-----------------------------*/


/*=========================
  Added a thing to EventSummary that this machine needs to fill 
*/
 
 
void pueo::UCorrelator::NoiseMachine::fillEventSummary(EventSummary *eventSummary) {
   
  

  for (int poli=0; poli<k::NUM_POLS; poli++) {
    
    for (int dir=0; dir<EventSummary::maxDirectionsPerPol; dir++) {
      double peakPhi = eventSummary->peak[poli][dir].phi;
      while (peakPhi < 0)    peakPhi += 360; //because the sun is translated stupidly
      while (peakPhi >= 360) peakPhi -= 360;
      double peakTheta = eventSummary->peak[poli][dir].theta;
      
      if (mapFifo[0][0] == NULL) {//can't do it if you haven't filled anything yet
	eventSummary->peak[poli][dir].mapHistoryVal = -999;
      }
      else {
	int iPhi = mapFifo[poli][0]->GetXaxis()->FindBin(peakPhi) - 1;
	int iTheta = mapFifo[poli][0]->GetYaxis()->FindBin(-peakTheta) - 1;
	if (iTheta >= nTheta) iTheta = nTheta -1; //Hmm why does Theta overflow occasionally?  This will fix it with edge effects
	double avgNoise = rollingMapAvg[rollingMapIndex(poli,iPhi,iTheta)];
	eventSummary->peak[poli][dir].mapHistoryVal = avgNoise;
      }
    }

  }

  //also fill up the history for the sources
  setSourceMapHistoryVal(eventSummary->sun);
  setSourceMapHistoryVal(eventSummary->wais);
  setSourceMapHistoryVal(eventSummary->ldb);

  return;
}



void pueo::UCorrelator::NoiseMachine::setSourceMapHistoryVal(EventSummary::SourceHypothesis& source) {

  for (int poli=0; poli<k::NUM_POLS; poli++) {
    double sourcePhi = source.phi - lastHeading;
    while (sourcePhi < 0)   sourcePhi += 360; //because the sun is translated stupidly
    while (sourcePhi >= 360) sourcePhi -= 360;

    double sourceTheta = source.theta;  
    if (mapFifo[0][0] == NULL) {//can't do it if you haven't filled anything yet
      source.mapHistoryVal[poli] = -999;
    }
    else {
      int iPhi = mapFifo[poli][0]->GetXaxis()->FindBin(sourcePhi) - 1;
      int iTheta = mapFifo[poli][0]->GetYaxis()->FindBin(-sourceTheta) - 1; //needs *-1 for sign reasons
      double avgNoise = rollingMapAvg[rollingMapIndex(poli,iPhi,iTheta)];
      source.mapHistoryVal[poli] = avgNoise;
      
      source.mapValue[poli] = mapFifo[poli][mapFifoPos]->GetBinContent(iPhi+1,iTheta+1);
    }

  }
    return;
}
