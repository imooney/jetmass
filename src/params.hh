// Parameters for efficiency analysis
// Veronica Verkest    5-20-18


//#include "TStopwatch.h"
#include "ktTrackEff.hh"
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include <utility>      // std::pair
#include <string>
#include <iostream>
#include <sstream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TRandom.h>

// TStarJetPico
#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

// Define a namespace for the variables

#ifndef PARAMS_HH
#define PARAMS_HH


#define __ERR(message) {std::cerr << "[" << __FILE__ << "::" << __func__ << "()] -- ERR: " << message << std::endl;}
#define __OUT(message) {std::cout << "[" << __FILE__ << "::" << __func__ << "()] -- OUT: " << message << std::endl;}

namespace Analysis {

  const std::string outFileName = "out/analysis.root";
  const double nEvents = -1;       //  NUMBER OF EVENTS  (-1 runs all)

  //consts                                                                                                                                                                        
  const double Pi = 3.141592653;
  const double PionMass = 0.13957018;

  const double R = 0.4;                   //jet resolution parameter                                                                                                              

  //quality cuts                                                  
  const std::string sim_triggerString = "All";
  const std::string dat_triggerString = "ppHT";
  const int refMultCut = 0;
  const double sim_evEtMin = -1;
  const double dat_evEtMin = 5.4;
  const double sim_maxEtTow = 9999;
  const double dat_maxEtTow = 100;
  const std::string sim_badTowers = "src/dummy_tower_list.txt";
  //truth 
  const double truth_absMaxVz = 1000;           //|Vz|<=30 cm
  const double truth_vZDiff = 1000;
  const double truth_evPtMax = 1000;
  const double truth_evEtMax = 1000;
  const double truth_DCA = 100;
  const double truth_NFitPts = -1;
  const double truth_FitOverMaxPts = -1;
  // const std::string truth_badTowers = "src/dummy_tower_list.txt";
  //detector
  const double det_absMaxVz = 30.0;           //|Vz|<=30 cm
  const double det_vZDiff = 31.0;             //max diff between selected TPC vertex and most probable VPD vertex (in ppRun6, VPD vz = 0, so vZDiff should be > absMaxVz)
  const double det_evPtMax = 30.0;
  const double det_evEtMax = 30.0;
  const double det_DCA = 3.0;
  const double det_NFitPts = 20;
  const double det_FitOverMaxPts = 0.52;
  const std::string det_badTowers = "src/y7_AuAu_HT_hot_list.txt";

  //particle cuts                                                                                                                                                                 
  const double max_track_rap = 1.0;
  const double partMinPt = 0.2;           //particle pT >= 0.2 GeV                                                                                                                

  //jet cuts                                                                                                                                                                      
  const double jet_ptmin = 2.0;           //jet pT >= 2.0 GeV                                                                                                                    
  const double jet_ptmax = 1000.0;        //DEBUG                                                                                                                                
  const double max_rap = max_track_rap-R; //|eta_jet| < 1-R
                                                                                                                                                                                 
  //ghosts                                                                                                                                                                        
  const int ghost_repeat = 1;
  const double ghost_area = 0.01;
  const double ghost_maxrap = max_rap + 2.0 * R;

}

#endif
