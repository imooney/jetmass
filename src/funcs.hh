//  functions.hh
//  Veronica Verkest May 13, 2018
//  Adapted by Isaac Mooney June, 2018

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"

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

//#include "ktTrackEff.hh"
#include <string>
#include <utility>      // std::pair
#include <iostream>
#include <sstream>

#ifndef funcs_hh
#define funcs_hh

namespace Analysis {

  // IO/OS MANIP Functions
                                                                                                                                                                                
  // Helper to build the TChain, used to decide which input format                                                                                                                                        
  bool HasEnding (std::string const &full_string, std::string const &ending);
  
  void FillTrees ( std::vector<fastjet::PseudoJet> jets, TTree* Tree, double &jPt, double &jEta, double &jPhi, double &jM, double &jE, int &jncons, double &wt, double weight); 

  void AnalysisSummary( int events, int pJets, int eJets, int gJets, int pgMatchedJets, int epMatchedJets, int egMatchedJets, std::string outName );
  
  std::vector<fastjet::PseudoJet> GatherParticles ( TStarJetVectorContainer<TStarJetVector> * container , double etaCutVal, double partMinPtVal, std::vector<fastjet::PseudoJet> & rawParticles );

  double LookupXsec(TString);
  
  void InitReader( TStarJetPicoReader & reader, TChain* chain, int nEvents, const std::string, const double, const double, const double, const double, const double, const double, const double, const double, const double, const std::string );

}

#endif
