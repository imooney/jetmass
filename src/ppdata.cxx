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
  InitReader(reader, chain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, dat_maxEtTow, 0.9999, false, /*"dummy_badtows.list"*/det_badTowers, dat_bad_run_list);

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
  
  double evt_vtx_JP2; double bbc_coinc_JP2; double runID_JP2;
  double n_trks_JP2; double n_tows_JP2; double n_vertices_JP2; double n_globals_JP2;
  vector<double> trackPt_JP2; vector<double> trackEta_JP2; vector<double> trackPhi_JP2; vector<double> trackDCA_JP2; vector<double> trackNhits_JP2; vector<double> trackNhitsposs_JP2;
  vector<double> towerEt_JP2; vector<double> towerEta_JP2; vector<double> towerPhi_JP2; vector<double> towerId_JP2;


  // ; towerEt vs towerID; average Et per event in towers vs luminosity;   ;trackDCA event_vertex  2D of z_vertex and bbc_coincidence.

  TTree *QAJP2 = new TTree("QAJP2","QAJP2");
  QAJP2->Branch("n_globals", &n_globals_JP2);
  QAJP2->Branch("n_trks", &n_trks_JP2);
  QAJP2->Branch("n_tows", &n_tows_JP2);
  QAJP2->Branch("n_vertices", &n_vertices_JP2);
  QAJP2->Branch("bbc_coinc", &bbc_coinc_JP2);
  QAJP2->Branch("evt_vtx", &evt_vtx_JP2); //vpdvz
  QAJP2->Branch("runID", &runID_JP2);
  QAJP2->Branch("trackPt", &trackPt_JP2);
  QAJP2->Branch("trackEta", &trackEta_JP2);
  QAJP2->Branch("trackPhi", &trackPhi_JP2);
  QAJP2->Branch("trackDCA", &trackDCA_JP2);
  QAJP2->Branch("trackNhits", &trackNhits_JP2);
  QAJP2->Branch("trackNhitsposs", &trackNhitsposs_JP2);
  QAJP2->Branch("towerEta", &towerEta_JP2);
  QAJP2->Branch("towerPhi", &towerPhi_JP2);
  QAJP2->Branch("towerEt", &towerEt_JP2);
  QAJP2->Branch("towerId", &towerId_JP2);


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

      n_globals_JP2 = -9999; n_trks_JP2 = -9999; n_tows_JP2 = -9999; bbc_coinc_JP2 = -9999; evt_vtx_JP2 = -9999; runID_JP2 = -9999; n_vertices_JP2 = -9999;
      trackPt_JP2.clear(); trackEta_JP2.clear(); trackPhi_JP2.clear(); trackDCA_JP2.clear(); trackNhits_JP2.clear(); trackNhitsposs_JP2.clear();
      towerEt_JP2.clear(); towerEta_JP2.clear(); towerPhi_JP2.clear(); towerId_JP2.clear();

      reader.PrintStatus(10);
      //get the event header
      
      header = reader.GetEvent()->GetHeader();
      event = reader.GetEvent();
      
      particles.clear();
      
      // Get the output container from the reader
      // ----------------------------------------
      container = reader.GetOutputContainer();

      //what are the trigger IDs for pp Y12?
      //JP2: 370621
      
      // if it HAS the trigger, do the QA for the triggered events. This should ALWAYS be satisfied for pp, since I'm selecting the JP2 trigger!
      if (header->HasTriggerId(370621)) {//!JP2 - only for ppY12! change if you change datasets! 
	/*n_tows_JP2 = header->GetNOfTowers(); n_trks_JP2 = header->GetNOfPrimaryTracks();*/ bbc_coinc_JP2 = header->GetBbcCoincidenceRate(); evt_vtx_JP2 = header->GetPrimaryVertexZ();
	n_vertices_JP2 = header->GetNumberOfVertices(); n_globals_JP2 = header->GetNGlobalTracks();
	
	for (int i = 0; i < header->GetNOfPrimaryTracks(); ++ i) {
          trackNhits_JP2.push_back(event->GetPrimaryTrack(i)->GetNOfFittedHits());
          trackNhitsposs_JP2.push_back(event->GetPrimaryTrack(i)->GetNOfPossHits());
        }

	TList *trks = reader.GetListOfSelectedTracks();
        TIter nxt_trk(trks);
	int n_trks = 0;
	
        while (TStarJetPicoPrimaryTrack* trk = (TStarJetPicoPrimaryTrack*) nxt_trk()) {
          if (fabs(trk->GetEta()) > 1.0 || trk->GetPt() < 0.2)
            continue;

          trackDCA_JP2.push_back(trk->GetDCA());
          trackPt_JP2.push_back(trk->GetPt());
          trackEta_JP2.push_back(trk->GetEta());
          trackPhi_JP2.push_back(trk->GetPhi());
	  //  trackNhits_JP2.push_back(trk->GetNOfFittedHits());
          //trackNhitsposs_JP2.push_back(trk->GetNOfPossHits());
          ++n_trks;
        }
	n_trks_JP2 = n_trks;
	
        TList *tows = reader.GetListOfSelectedTowers();
        TIter nxt_tow(tows);
	int n_tows = 0;
	
        while (TStarJetPicoTower* tow = (TStarJetPicoTower*) nxt_tow()) {
          if (fabs(tow->GetEta()) > 1.0 || tow->GetEt() < 0.2)
            continue;

          towerEta_JP2.push_back(tow->GetEta());
          towerPhi_JP2.push_back(tow->GetPhi());
          towerEt_JP2.push_back(tow->GetEt());
          towerId_JP2.push_back(tow->GetId());
	  ++n_tows;
        }
	n_tows_JP2 = n_tows;

	QAJP2->Fill();
	//continue;
      }
      else {
	cerr << "You should never see this message! The JP2 trigger is selected, so should never have a case in which an event doesn't have that trigger ID!" << endl; exit(1);
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
  QAJP2->Write();
  //fout->Write();
  fout->Close();

  cout << "In " << nEventsUsed << " events, found " << endl
       << nJets << " jets above 5 GeV, with constituents above 0.2 GeV." << endl;  

  cout << "Wrote to " << fout->GetName() << endl;
  cout << "Bye :)" << endl;
    
  return 0;
}
  
