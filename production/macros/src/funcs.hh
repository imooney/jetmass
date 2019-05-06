//  functions.hh
//  Veronica Verkest May 13, 2018
//  Adapted by Isaac Mooney June, 2018

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "fastjet/contrib/SoftDrop.hh"

#include "RooUnfold.h"

#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"

#include <string>
#include <utility>      // std::pair
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <memory>
#include <set>
#include <fstream>

#ifndef funcs_hh
#define funcs_hh

const double R = 0.4;

namespace Analysis {
  std::vector<int> MatchJets(const std::vector<fastjet::PseudoJet>, const std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet> &, std::vector<fastjet::PseudoJet> &);

  std::vector<int> FakesandMisses(const std::vector<fastjet::PseudoJet>, const std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet> &);
  
  void ConstructResponses(std::vector<RooUnfoldResponse*>, const std::vector<fastjet::PseudoJet>, const std::vector<fastjet::PseudoJet>, const std::vector<fastjet::PseudoJet>, const std::vector<fastjet::PseudoJet>, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, const double);

  void HistFromTree(TFile *onfile, TH1D*, TH1D*,TH1D*,TH1D*,TH1D*,TH1D* /*TFile *offfile, std::vector<RooUnfoldResponse*> res_vec*/);
}
#endif
