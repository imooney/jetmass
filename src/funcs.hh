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
#include <unordered_map>
#include <memory>
#include <set>
#include <fstream>

#ifndef funcs_hh
#define funcs_hh

using std::unordered_map; using std::make_shared; using std::shared_ptr; 

namespace Analysis {

  // IO/OS MANIP Functions
                                                                                                                                                                                
  // Helper to build the TChain, used to decide which input format                                                                                                                                        
  bool HasEnding (std::string const &, std::string const &);

  template <typename T>
  std::set<T> ParseCSV(std::string);

  template<typename T>
  bool CanCast(std::string);

  template<typename T>
  T CastTo(std::string);  
  
  void FillMatchedTree (const fastjet::PseudoJet, const fastjet::PseudoJet, TTree *, double &, double &, double &, double &, double &, double &, double &, double &, const double);
  
  void FillTrees ( std::vector<fastjet::PseudoJet>, TTree*, double &, double &, double &, double &, double &, int &, double &, const double); 

  void AnalysisSummary( int, int, int, int, int, int, int, std::string);
  
  void GatherParticles ( TStarJetVectorContainer<TStarJetVector> *, TStarJetVector*, std::vector<fastjet::PseudoJet> &, const bool, const bool);

  double LookupRun6Xsec(TString);

  double LookupRun12Xsec( TString );
  
  bool GetTriggers(std::vector<fastjet::PseudoJet> &, const std::vector<fastjet::PseudoJet>);

  bool GetTriggerJet(std::vector<fastjet::PseudoJet> &, const std::vector<fastjet::PseudoJet>);

  void MatchJets(const std::vector<fastjet::PseudoJet>, const fastjet::PseudoJet, int &);

  void ConstructResponses(RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, const fastjet::PseudoJet, const fastjet::PseudoJet, const double);

  void ConstructResponsesSD(RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, const fastjet::PseudoJet, const fastjet::PseudoJet, const double);

  void Misses(RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, const fastjet::PseudoJet, const double);

  void MissesSD(RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, const fastjet::PseudoJet, const double);

  void Fakes(RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, const fastjet::PseudoJet, const double);

  void FakesSD(RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, RooUnfoldResponse &, const fastjet::PseudoJet, const double);  
  
  void InitReader( TStarJetPicoReader &, TChain*, int, const std::string, const double, const double, const double, const double, const double, const double, const double, const double, const double, const std::string, const std::string);

  //HISTOGRAMS
  template<class Key, class H, class hash=std::hash<Key>>
    class Collection {
    public:
      Collection() : collection_() { };
      ~Collection() { };
  
      H* get(Key key) {
	if (keyExists(key))
	  return collection_[key].get();
	return nullptr;
      }
  
      template <typename... Args>
      void add(Key key, Args... args) {
	collection_[key] = make_shared<H>(key.c_str(), args...);
	collection_[key]->SetDirectory(0);
      }
  
      template <typename ...Args>
      bool fill(Key key, Args... args) {
	if (!keyExists(key))
	  return false;
	collection_[key]->Fill(args...);
	return true;
      }
  
      bool write(Key key) {
	if (!keyExists(key))
	  return false;
	collection_[key]->Write();
	return true;
      }

      void clear() {
	collection_.clear();
      }
  
    private:
  
      unordered_map<Key, shared_ptr<H>, hash> collection_ = {};
  
      bool keyExists(std::string key) {
	for (auto& h : collection_) {
	  if (h.first == key)
	    return true;
	}
	return false;
      }
  
    };

  void FillHistsHelper(Collection<std::string, TH1D> &, Collection<std::string, TH2D> &, Collection<std::string, TH3D> &, const std::string, const std::string, const std::string, const fastjet::PseudoJet, const double);

  void FillHists(Collection<std::string, TH1D> &, Collection<std::string, TH2D> &, Collection<std::string, TH3D> &, const std::string, const std::string, const std::vector<fastjet::PseudoJet>, const double);

  
  void FillSDHistsHelper(Collection<std::string, TH1D> &, Collection<std::string, TH2D> &, Collection<std::string, TH3D> &, const std::string, const std::string, const std::string, const fastjet::PseudoJet, const double);

  void FillSDHists(Collection<std::string, TH1D> &, Collection<std::string, TH2D> &, Collection<std::string, TH3D> &, const std::string, const std::string, const std::vector<fastjet::PseudoJet>, const double);
  
}

#endif
