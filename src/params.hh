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
  const int refMultCut = 0;
  const double truth_evEtMin = -1;
  const double det_evEtMin = -1;
  const double sim_maxEtTow = 9999;
  const double dat_maxEtTow = 100;
  const std::string sim_badTowers = "src/dummy_tower_list.txt";
  const std::string sim_bad_run_list = "dummy_badrun.list";
  const std::string dat_bad_run_list = "pp200Y12_badrun.list";
  
  //truth 
  const std::string truth_triggerString = "All";
  const double truth_absMaxVz = 1000;           //|Vz|<=30 cm
  const double truth_vZDiff = 1000;
  const double truth_evPtMax = 1000;
  const double truth_evEtMax = 1000;
  const double truth_DCA = 100;
  const double truth_NFitPts = -1;
  const double truth_FitOverMaxPts = -1;
  // const std::string truth_badTowers = "src/dummy_tower_list.txt";
 
  //detector
  
  const std::string det_triggerString = "ppJP2";
  const double det_absMaxVz = 30.0;           //|Vz|<=30 cm
  const double det_vZDiff = 31.0;             //max diff between selected TPC vertex and most probable VPD vertex (in ppRun6, VPD vz = 0, so vZDiff should be > absMaxVz)
  const double det_evPtMax = 30.0;
  const double det_evEtMax = 30.0;
  const double det_DCA = 1.0;
  const double det_NFitPts = 20;
  const double det_FitOverMaxPts = 0.52;
  const std::string det_badTowers = "Combined_pp200Y12_badtower.list";//"src/y7_AuAu_HT_hot_list.txt";
  //particle cuts                                                                                                                                                                 
  const double max_track_rap = 1.0;
  const double partMinPt = 0.2;           //30.0 GeV >= particle pT >= 0.2 GeV 
  const double partMaxPt = 30.0;         

  //jet cuts                                                                                                                                                                      
  const double jet_ptmin = 5.0;           //jet pT >= 2.0 GeV                 
  const double jet_ptmax = 1000.0;        //DEBUG
  const double max_rap = max_track_rap-R; //|eta_jet| < 1-R
  const double NEF_max = 0.9;             //neutral energy fraction of jet must be < 90%                                                                                                                                                                                      
  //ghosts                                                                                                                                                                        
  const int ghost_repeat = 1;
  const double ghost_area = 0.01;
  const double ghost_maxrap = max_rap + 2.0 * R;

  //softdrop params
  
  const double z_cut = 0.10;
  const double beta = 0.0;
  const double R0 = 1.0;

}

#endif
