#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TFile.h>

#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TChain.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <chrono>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"

#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

#include "fastjet/contrib/SoftDrop.hh"

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

#include "funcs.hh"
#include "params.hh"

using namespace std;
using namespace fastjet;
using namespace Analysis;
typedef fastjet::contrib::SoftDrop SD;

// -------------------------
// Command line arguments: ( Defaults
// Defined for debugging in main )
// [0]: output directory
// [1]: name for the histogram file
// [2]: input data: can be a single .root or a .txt or .list of root files - should always be last argument

// DEF MAIN()
int main ( int argc, const char** argv) {
    
  //Start a timer
  TStopwatch TimeKeeper;
  TimeKeeper.Start( );
    
  //starting timing for function duration
  typedef std::chrono::high_resolution_clock clock;

  // Read in command line arguments
  // ------------------------------
  // Defaults
  std::string     executable    = "./bin/ppjetmass";     // placeholder
  std::string        outputDir         = "out/";                                        // directory where everything will be saved
  std::string     outFileName        = "test.root";                        // histograms will be saved here
  std::string         chainList            = "list.txt";            // input file: can be .root, .txt, .list
  std::string  chainName     = "JetTree";                                // Tree name in input file
  bool full = 1;
  // Now check to see if we were given modifying arguments
  switch ( argc ) {
  case 1: // Default case
    __OUT("Using Default Settings");
      break;
  case 6: { // Custom case
    __OUT("Using Custom Settings");
      std::vector<std::string> arguments( argv+1, argv+argc );
            
    // Set non-default values
    // ----------------------
            
    // output and file names
    outputDir         = arguments[0];
    outFileName       = arguments[1];
    if (arguments[2] == "ch") {full = 0;} else {full = 1;}
    chainList         = arguments[4];
            
    break;
  }
  default: { // Error: invalid custom settings
    __ERR("Invalid number of command line arguments");
      return -1;
    break;
  }
  }

  // Build our input now
  // --------------------
  TChain* chain = new TChain( chainName.c_str() );
  
  // Check to see if the input is a .root file or a .txt
  bool inputIsRoot = Analysis::HasEnding( chainList.c_str(), ".root" );
  bool inputIsTxt  = Analysis::HasEnding( chainList.c_str(), ".txt"  );
  bool inputIsList = Analysis::HasEnding( chainList.c_str(), ".list" );
  
  // If its a recognized file type, build the chain
  // If its not recognized, exit
  if ( inputIsRoot ) { chain->Add( chainList.c_str() ); }
  else if ( inputIsTxt )  { chain = TStarJetPicoUtils::BuildChainFromFileList( chainList.c_str() ); }
  else if ( inputIsList)  { chain = TStarJetPicoUtils::BuildChainFromFileList( chainList.c_str() ); }
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }
	
  // Build the event structure w/ cuts
  // ---------------------------------
  TStarJetPicoReader reader;
  InitReader(reader, chain, nEvents, pAu_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, dat_maxEtTow, 0.9999, false, /*"dummy_badtows.list"*/pAu_badTowers, pAu_bad_run_list);

  // Data classes
  // ------------
  TStarJetVectorContainer<TStarJetVector>* container;
  TStarJetVector* sv; // TLorentzVector* would be sufficient
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  
  // Histograms
  // ----------
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  // Trees
  // -----
  int dummy_int;
  double dummy_double;

  double n_jets;
  vector<double> Pt; vector<double> Eta; vector<double> Phi; vector<double> M; vector<double> E;
  vector<double> ch_e_frac;
  vector<double> zg; vector<double> rg; vector<double> mg; vector<double> ptg;
  vector<double> mcd;
  vector<double> tau0; vector<double> tau05; vector<double> tau_05; vector<double> tau_1;
  vector<double> tau0_g; vector<double> tau05_g; vector<double> tau_05_g; vector<double> tau_1_g;
  
  TTree *eventTree = new TTree("event","event");
  eventTree->Branch("n_jets", &n_jets);
  eventTree->Branch("Pt", &Pt); eventTree->Branch("Eta",&Eta); eventTree->Branch("Phi",&Phi); eventTree->Branch("M",&M); eventTree->Branch("E",&E);
  eventTree->Branch("ch_e_frac",&ch_e_frac);
  eventTree->Branch("zg", &zg); eventTree->Branch("rg", &rg); eventTree->Branch("mg", &mg); eventTree->Branch("ptg",&ptg);
  eventTree->Branch("mcd",&mcd);
  eventTree->Branch("tau0",&tau0); eventTree->Branch("tau05",&tau05); eventTree->Branch("tau_05",&tau_05); eventTree->Branch("tau_1",&tau_1);
  eventTree->Branch("tau0_g",&tau0_g); eventTree->Branch("tau05_g",&tau05_g); eventTree->Branch("tau_05_g",&tau_05_g); eventTree->Branch("tau_1_g",&tau_1_g);
  
  double evt_vtx_BBCMB; double bbc_coinc_BBCMB;
  double n_trks_BBCMB; double n_tows_BBCMB;
  vector<double> trackPt_BBCMB; vector<double> trackEta_BBCMB; vector<double> trackPhi_BBCMB; vector<double> trackDCA_BBCMB;
  vector<double> towerEt_BBCMB; vector<double> towerEta_BBCMB; vector<double> towerPhi_BBCMB; vector<double> towerId_BBCMB;

  double evt_vtx_VPDMB; double bbc_coinc_VPDMB;
  double n_trks_VPDMB; double n_tows_VPDMB;
  vector<double> trackPt_VPDMB; vector<double> trackEta_VPDMB; vector<double> trackPhi_VPDMB; vector<double> trackDCA_VPDMB;
  vector<double> towerEt_VPDMB; vector<double> towerEta_VPDMB; vector<double> towerPhi_VPDMB; vector<double> towerId_VPDMB;


  // ; towerEt vs towerID; average Et per event in towers vs luminosity;   ;trackDCA event_vertex  2D of z_vertex and BBC_coincidence.

  TTree *QABBCMB = new TTree("QABBCMB","QABBCMB");
  QABBCMB->Branch("n_trks", &n_trks_BBCMB);
  QABBCMB->Branch("n_tows", &n_tows_BBCMB);
  QABBCMB->Branch("bbc_coinc", &bbc_coinc_BBCMB);
  QABBCMB->Branch("evt_vtx", &evt_vtx_BBCMB); //vpdvz
  QABBCMB->Branch("trackPt", &trackPt_BBCMB);
  QABBCMB->Branch("trackEta", &trackEta_BBCMB);
  QABBCMB->Branch("trackPhi", &trackPhi_BBCMB);
  QABBCMB->Branch("trackDCA", &trackDCA_BBCMB);
  QABBCMB->Branch("towerEta", &towerEta_BBCMB);
  QABBCMB->Branch("towerPhi", &towerPhi_BBCMB);
  QABBCMB->Branch("towerEt", &towerEt_BBCMB);
  QABBCMB->Branch("towerId", &towerId_BBCMB);

  TTree *QAVPDMB = new TTree("QAVPDMB","QAVPDMB");
  QAVPDMB->Branch("n_trks", &n_trks_VPDMB);
  QAVPDMB->Branch("n_tows", &n_tows_VPDMB);
  QAVPDMB->Branch("bbc_coinc", &bbc_coinc_VPDMB);
  QAVPDMB->Branch("evt_vtx", &evt_vtx_VPDMB); //vpdvz
  QAVPDMB->Branch("trackPt", &trackPt_VPDMB);
  QAVPDMB->Branch("trackEta", &trackEta_VPDMB);
  QAVPDMB->Branch("trackPhi", &trackPhi_VPDMB);
  QAVPDMB->Branch("trackDCA", &trackDCA_VPDMB);
  QAVPDMB->Branch("towerEta", &towerEta_VPDMB);
  QAVPDMB->Branch("towerPhi", &towerPhi_VPDMB);
  QAVPDMB->Branch("towerEt", &towerEt_VPDMB);
  QAVPDMB->Branch("towerId", &towerId_VPDMB);


  // Helpers
  // -------
  vector<PseudoJet> particles;
  
  int nJets = 0;

  // Constituent selectors
  // ---------------------
  Selector select_track_rap = fastjet::SelectorAbsRapMax(max_track_rap);
  Selector select_lopt      = fastjet::SelectorPtMin( partMinPt );
  Selector select_loptmax   = fastjet::SelectorPtMax( partMaxPt );
  Selector slo = select_track_rap * select_lopt * select_loptmax;

  // Jet candidate selectors
  // -----------------------
  Selector select_jet_rap     = fastjet::SelectorAbsRapMax(max_rap);
  Selector select_jet_pt_min  = fastjet::SelectorPtMin( det_jet_ptmin );
  Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
  Selector select_jet_m_min = fastjet::SelectorMassMin( mass_min );
  Selector sjet = select_jet_rap && select_jet_pt_min && select_jet_pt_max && select_jet_m_min;
  
  // Choose a jet and area definition
  // --------------------------------
  JetDefinition jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, R);
  
  // create an area definition for the clustering
  //----------------------------------------------------------
  // ghosts should go up to the acceptance of the detector or
  // (with infinite acceptance) at least 2R beyond the region
  // where you plan to investigate jets.
  GhostedAreaSpec area_spec = fastjet::GhostedAreaSpec( ghost_maxrap, ghost_repeat, ghost_area );
  AreaDefinition  area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);

  //Creating SoftDrop grooming object
  contrib::SoftDrop sd(Beta,z_cut,R0);
  cout << "SoftDrop groomer is: " << sd.description() << endl;
  
  cout << "Performing analysis." << endl;
  // Cycle through events
  // --------------------  
  //int nEvents = -1;
  int nEventsUsed = 0;

  //for later use looking up PDG masses using particle PID                                                                                                                 
  TDatabasePDG *pdg = new TDatabasePDG();
  
  //initialize the reader
  //reader.Init( nEvents ); //runs through all events with -1

  try{
    while ( reader.NextEvent() ) {
      
      //clearing vectors
      Pt.clear(); Eta.clear(); Phi.clear(); M.clear(); E.clear();
      zg.clear(); rg.clear(); mg.clear(); ptg.clear();
      ch_e_frac.clear(); mcd.clear();
      tau0.clear(); tau05.clear(); tau_05.clear(); tau_1.clear();
      tau0_g.clear(); tau05_g.clear(); tau_05_g.clear(); tau_1_g.clear();
      //initializing variables to -9999
      n_jets = -9999;

      n_trks_BBCMB = -9999; n_tows_BBCMB = -9999; bbc_coinc_BBCMB = -9999; evt_vtx_BBCMB = -9999;
      trackPt_BBCMB.clear(); trackEta_BBCMB.clear(); trackPhi_BBCMB.clear(); trackDCA_BBCMB.clear();
      towerEt_BBCMB.clear(); towerEta_BBCMB.clear(); towerPhi_BBCMB.clear(); towerId_BBCMB.clear();
      
      n_trks_VPDMB = -9999; n_tows_VPDMB = -9999; bbc_coinc_VPDMB = -9999; evt_vtx_VPDMB = -9999;
      trackPt_VPDMB.clear(); trackEta_VPDMB.clear(); trackPhi_VPDMB.clear(); trackDCA_VPDMB.clear();
      towerEt_VPDMB.clear(); towerEta_VPDMB.clear(); towerPhi_VPDMB.clear(); towerId_VPDMB.clear();

      reader.PrintStatus(10);
      //get the event header
      
      header = reader.GetEvent()->GetHeader();
      event = reader.GetEvent();
      
      particles.clear();
      
      // Get the output container from the reader
      // ----------------------------------------
      container = reader.GetOutputContainer();
      
      
      //for pAu, need to ask if the event has the trigger. The trigger IDs are:
      //HT2*BBCMB : 500205, 500215
      //JP2 : 500401, 500411
      //BBCMB : 500008, 500018
      //VPDMB :  500904
      
      //JET INFO IS USELESS UNTIL YOU PICK ONE TRIGGER AND "continue" AFTER IT. RIGHT NOW ALL TRIGGERS ARE TAKEN DURING QA
      //ONLY USE THIS FOR pAu, trigger ids can change by run!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //      if (!(header->HasTriggerId(500401) || header->HasTriggerId(500411))) {//!JP2
	//continue;
      // }
      if (header->HasTriggerId(500008) || header->HasTriggerId(500018)) {//!BBCMB
	n_trks_BBCMB = header->GetNOfTowers(); n_tows_BBCMB = header->GetNOfPrimaryTracks(); bbc_coinc_BBCMB = header->GetBbcCoincidenceRate(); evt_vtx_BBCMB = header->GetPrimaryVertexZ();
	
	for (int i = 0; i < header->GetNOfPrimaryTracks(); ++ i) {
	  trackDCA_BBCMB.push_back(event->GetPrimaryTrack(i)->GetDCA());
	  trackPt_BBCMB.push_back(event->GetPrimaryTrack(i)->GetPt());
	  trackEta_BBCMB.push_back(event->GetPrimaryTrack(i)->GetEta());
	  trackPhi_BBCMB.push_back(event->GetPrimaryTrack(i)->GetPhi());
	}
	for (int i = 0; i < header->GetNOfTowers(); ++ i) {
	  towerEta_BBCMB.push_back(event->GetTower(i)->GetEta());
	  towerPhi_BBCMB.push_back(event->GetTower(i)->GetPhi());
	  towerEt_BBCMB.push_back(event->GetTower(i)->GetEt());
	  towerId_BBCMB.push_back(event->GetTower(i)->GetId());

	}
	QABBCMB->Fill();
	//continue;
      }
      if (header->HasTriggerId(500904)) {//!VPDMB
	n_trks_VPDMB = header->GetNOfTowers(); n_tows_VPDMB = header->GetNOfPrimaryTracks(); bbc_coinc_VPDMB = header->GetBbcCoincidenceRate(); evt_vtx_VPDMB = header->GetPrimaryVertexZ();
	
	for (int i = 0; i < header->GetNOfPrimaryTracks(); ++ i) {
	  trackDCA_VPDMB.push_back(event->GetPrimaryTrack(i)->GetDCA());
	  trackPt_VPDMB.push_back(event->GetPrimaryTrack(i)->GetPt());
	  trackEta_VPDMB.push_back(event->GetPrimaryTrack(i)->GetEta());
	  trackPhi_VPDMB.push_back(event->GetPrimaryTrack(i)->GetPhi());
	}
	for (int i = 0; i < header->GetNOfTowers(); ++ i) {
	  towerEta_VPDMB.push_back(event->GetTower(i)->GetEta());
	  towerPhi_VPDMB.push_back(event->GetTower(i)->GetPhi());
	  towerEt_VPDMB.push_back(event->GetTower(i)->GetEt());
	  towerId_VPDMB.push_back(event->GetTower(i)->GetId());

	}
	QAVPDMB->Fill();

//continue;
      }

      // Transform TStarJetVectors into (FastJet) PseudoJets
      // ----------------------------------------------------------
      GatherParticles(container, sv, particles, full,0, pdg); //ch+ne

      // Analysis
      // --------
      // Apply selector to the full particle set
      vector<PseudoJet> pLo = slo( particles );
      
      // find corresponding jets with soft constituents
      // ----------------------------------------------
      ClusterSequence/*Area*/ csaLo ( pLo, jet_def/*, area_def */); // WITHOUT background subtraction
      
      /*
      // Background initialization
      // -------------------------
      // Background selector - Exclude two hardest jets for background extermination
      Selector selector_bkgd = fastjet::SelectorAbsRapMax( max_rap ) * (!fastjet::SelectorNHardest(2));
      // Area - same as for jets
      AreaDefinition area_def_bkgd ( area_def );
      // Jet definition - use kT instead of anti-kT algorithm here
      JetDefinition jet_def_bkgd (fastjet::kt_algorithm, R );
      // Energy density estimate from median ( pt_i / area_i )
      JetMedianBackgroundEstimator bkgd_estimator (selector_bkgd, jet_def_bkgd, area_def_bkgd);
      bkgd_estimator.set_particles( pLo );
      // Subtract A*rho from the original pT & m
      Subtractor bkgd_subtractor (&bkgd_estimator);
      bkgd_subtractor.set_use_rho_m();*/
      
      vector<PseudoJet> LoInitial = fastjet::sorted_by_pt(sjet(/*bkgd_subtractor(*/csaLo.inclusive_jets())/*)*/);
      vector<PseudoJet> LoResult;

      //Implementing a neutral energy fraction cut of 90% on inclusive jets
      ApplyNEFSelection(LoInitial, LoResult);
      
      vector<PseudoJet> GroomedJets;
      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
      for (int i = 0; i < LoResult.size(); ++ i) {
	GroomedJets.push_back(sd(LoResult[i]));
      }
      
      if (LoResult.size() != 0) {
	//SoftDrop is a groomer not a tagger, so if we have at least one ungroomed jet, we should also have a SoftDrop'd jet.
	nJets += LoResult.size();
	n_jets = LoResult.size();
	for (int i = 0; i < n_jets; ++ i) {
	  Pt.push_back(LoResult[i].pt()); Eta.push_back(LoResult[i].eta()); Phi.push_back(LoResult[i].phi());
	  M.push_back(LoResult[i].m()); E.push_back(LoResult[i].e());
	  zg.push_back(GroomedJets[i].structure_of<SD>().symmetry()); rg.push_back(GroomedJets[i].structure_of<SD>().delta_R());
	  mg.push_back(GroomedJets[i].m()); ptg.push_back(GroomedJets[i].pt());
	  
	  double m2 = (LoResult[i].m())*(LoResult[i].m()); double gm2 = (GroomedJets[i].m())*(GroomedJets[i].m());
	  double m_cd = (double) sqrt(m2 - gm2); if ((m2 - gm2) < 10e-10) {m_cd = 0;}
	  mcd.push_back(m_cd);
	  
	  double ch_e = 0; double tot_e = 0;
	  double sum_tau0 = 0; double sum_tau05 = 0; double sum_tau_05 = 0; double sum_tau_1 = 0;
	  double sum_tau0_g = 0; double sum_tau05_g = 0; double sum_tau_05_g = 0; double sum_tau_1_g = 0;
	  vector<PseudoJet> cons = LoResult[i].constituents();
	  vector<PseudoJet> cons_g = GroomedJets[i].constituents();
	  //using pT, not energy, so some names are misnomers
	  for (int j = 0; j < cons.size(); ++ j) {
	    if (cons[j].user_index() != 0) {ch_e += cons[j].pt();}
	    tot_e += cons[j].pt();
	    //angularity:
	    sum_tau0 += (cons[j].pt()*pow(LoResult[i].delta_R(cons[j]), 2 - 0));
	    sum_tau05 += (cons[j].pt()*pow(LoResult[i].delta_R(cons[j]), 2 - 0.5));
	    sum_tau_05 += (cons[j].pt()*pow(LoResult[i].delta_R(cons[j]), 2 + 0.5));
	    sum_tau_1 += (cons[j].pt()*pow(LoResult[i].delta_R(cons[j]), 2 + 1));
	  }
	  for (int j = 0; j < cons_g.size(); ++j) { 
	    sum_tau0_g += (cons_g[j].pt()*pow(GroomedJets[i].delta_R(cons_g[j]), 2 - 0));
	    sum_tau05_g += (cons_g[j].pt()*pow(GroomedJets[i].delta_R(cons_g[j]), 2 - 0.5));
	    sum_tau_05_g += (cons_g[j].pt()*pow(GroomedJets[i].delta_R(cons_g[j]), 2 + 0.5));
	    sum_tau_1_g += (cons_g[j].pt()*pow(GroomedJets[i].delta_R(cons_g[j]), 2 + 1));
	  }
	  tau0.push_back(sum_tau0 / (double) LoResult[i].pt());
	  tau05.push_back(sum_tau05 / (double) LoResult[i].pt());
	  tau_05.push_back(sum_tau_05 / (double) LoResult[i].pt());
	  tau_1.push_back(sum_tau_1 / (double) LoResult[i].pt());
	  tau0_g.push_back(sum_tau0_g / (double) GroomedJets[i].pt());
	  tau05_g.push_back(sum_tau05_g / (double) GroomedJets[i].pt());
	  tau_05_g.push_back(sum_tau_05_g / (double) GroomedJets[i].pt());
	  tau_1_g.push_back(sum_tau_1_g / (double) GroomedJets[i].pt());

	  ch_e_frac.push_back(ch_e/(double)tot_e);
	}
      }
      
      // And we're done!
      // -----------------------------
      if (LoResult.size() != 0) {
	nEventsUsed++;
	eventTree->Fill();
      }
    } // Event loop
  }catch ( std::exception& e) {
    std::cerr << "Caught " << e.what() << std::endl;
    return -1;
  }

  // Output
  // ------                                                                                                                                                                                               
  TFile* fout = new TFile((outputDir + outFileName).c_str(), "RECREATE");

  // Close up shop
  // -------------
  eventTree->Write();
  QABBCMB->Write(); QAVPDMB->Write();
  //fout->Write();
  fout->Close();

  cout << "In " << nEventsUsed << " events, found " << endl
       << nJets << " jets above 5 GeV, with constituents above 0.2 GeV." << endl;  

  cout << "Wrote to " << fout->GetName() << endl;
  cout << "Bye :)" << endl;
    
  return 0;
}
  
