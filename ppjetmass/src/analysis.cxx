#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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
    
  // Read in command line arguments
  // ------------------------------
  // Defaults
  std::string     executable    = "./bin/ppjetmass";     // placeholder
  std::string        outputDir         = "out/";                                        // directory where everything will be saved
  std::string     outFileName        = "test.root";                        // histograms will be saved here
  std::string         chainList            = "list.txt";            // input file: can be .root, .txt, .list
  std::string  chainName     = "JetTree";                                // Tree name in input file
  std::string badTowerList  = "src/y7_AuAu_HT_hot_list.txt";
    
  // Now check to see if we were given modifying arguments
  switch ( argc ) {
  case 1: // Default case
    __OUT("Using Default Settings");
      break;
  case 4: { // Custom case
    __OUT("Using Custom Settings");
      std::vector<std::string> arguments( argv+1, argv+argc );
            
    // Set non-default values
    // ----------------------
            
    // output and file names
    outputDir         = arguments[0];
    outFileName       = arguments[1];
    chainList         = arguments[2];
            
    break;
  }
  default: { // Error: invalid custom settings
    __ERR("Invalid number of command line arguments");
      return -1;
    break;
  }
  }

  /*
int main ( int argc, const char** argv ) {
  
	// set of default arguments that will work
        const char* defaults[5] = {"SimpleJetFinder", "1", "out/tmp.root", "list.txt" };
	if (argc == 1) {
		argv = defaults;
		argc = 4;
	}
 
	std::vector<std::string> arguments(argv+1, argv + argc);
 
	// load variables
 
	bool    chainFromList = atof( arguments.at(0).c_str() );
	TString outFileName   = arguments.at(1);
	TString chainList     = arguments.at(2);
	TString chainName     = "JetTree";
	TString badTowerList  = "src/y7_AuAu_HT_hot_list.txt";
  
	// Input
	// -----
	// create the data handling classes
	TChain* chain         = new TChain( chainName );
	if ( chainFromList == false ) {
		for ( int i = 2; i < arguments.size(); ++i ) {
			chain->Add( arguments.at(i).data() );
		}
	}
	else if ( chainFromList == true ) {
		chain = TStarJetPicoUtils::BuildChainFromFileList( chainList );
	}
  */

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
	reader.SetInputChain (chain);
	reader.SetApplyFractionHadronicCorrection( true );
	reader.SetFractionHadronicCorrection( 0.9999 );
	reader.SetRejectTowerElectrons( kFALSE );

	// Event and track selection
	// -------------------------
	//ref mult cut and trigger selection
	int refMultCut = 0;
	TString triggerString = "All";
 
	TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
	evCuts->SetTriggerSelection( triggerString ); //All, MB, HT, pp, ppHT, ppJP
	// Additional cuts
	evCuts->SetVertexZCut ( absMaxVz );
	evCuts->SetRefMultCut ( refMultCut );
	evCuts->SetVertexZDiffCut( 1000 );
	evCuts->SetMaxEventPtCut( 30 );
	evCuts->SetMaxEventEtCut( 30 );
 
	// Tracks cuts
	TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
	trackCuts->SetDCACut( DCA );
	trackCuts->SetMinNFitPointsCut( NFitPts ); //AjParameters::NMinFit
	trackCuts->SetFitOverMaxPointsCut( FitOverMaxPts ); //AjParameters::FitOverMaxPointsCut
	trackCuts->SetMaxPtCut ( MaxPt );
 
	std::cout << "Using these track cuts:" << std::endl;
	std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
	std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
	std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
	std::cout << " maxpt : " << trackCuts->GetMaxPtCut (  ) << std::endl;
 
	// Towers
	TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
	towerCuts->SetMaxEtCut( 100 );
	towerCuts->AddBadTowers( badTowerList );
 
	std::cout << "Using these tower cuts:" << std::endl;
	std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
	std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
 
	// V0s: Turn off
	reader.SetProcessV0s(false);
	
	// Data classes
	// ------------
	TStarJetVectorContainer<TStarJetVector>* container;
	TStarJetVector* sv; // TLorentzVector* would be sufficient
	TStarJetPicoEventHeader* header;

  // Histograms
  // ----------
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  TH3 * PtEtaPhi_incl = new TH3D("PtEtaPhi_incl","",80,0,80,5,-1,1,10,0,2*Pi);
  TH3 * PtEtaPhi_lead = new TH3D("PtEtaPhi_lead","",80,0,80,5,-1,1,10,0,2*Pi);

  TH3 * cons_PtEtaPhi_incl = new TH3D("cons_PtEtaPhi_incl","",80,0,80,5,-1,1,10,0,2*Pi);
  TH3 * cons_PtEtaPhi_lead = new TH3D("cons_PtEtaPhi_lead","",80,0,80,5,-1,1,10,0,2*Pi);
  
  TH1 * m_incl = new TH1D("m_inclusive","",80,0,80);
  TH1 * m_lead = new TH1D("m_leading","",80,0,80);

  TH2 * m_v_pt_lead = new TH2D("m_v_pt_lead","",20,0,20,80,0,80);
  TH2 * m_v_pt_incl = new TH2D("m_v_pt_incl","",20,0,20,80,0,80);

  // Trees
  // -----
  double lead_Pt, lead_Eta, lead_Phi, lead_M, lead_E;
  double cons_lead_Pt, cons_lead_Eta, cons_lead_Phi, cons_lead_M, cons_lead_E;
  double incl_Pt, incl_Eta, incl_Phi, incl_M, incl_E;
  double cons_incl_Pt, cons_incl_Eta, cons_incl_Phi, cons_incl_M, cons_incl_E;



  TTree *leadTree = new TTree("lead","lead");
  TTree *cons_leadTree = new TTree("cons_lead","cons_lead");
  TTree *inclTree = new TTree("incl","incl");
  TTree *cons_inclTree = new TTree("cons_incl","cons_incl");

  TBranch *leadPt;  TBranch *leadEta;  TBranch *leadPhi; TBranch *leadM; TBranch *leadE;
  leadPt = leadTree->Branch("lead_Pt", &lead_Pt);  leadEta = leadTree->Branch("lead_Eta", &lead_Eta);  leadPhi = leadTree->Branch("lead_Phi", &lead_Phi);
  leadM = leadTree->Branch("lead_M", &lead_M); leadE = leadTree->Branch("lead_E", &lead_E);

  TBranch *cons_leadPt;  TBranch *cons_leadEta;  TBranch *cons_leadPhi; TBranch *cons_leadM; TBranch *cons_leadE;
  cons_leadPt = cons_leadTree->Branch("cons_lead_Pt", &cons_lead_Pt); cons_leadEta = cons_leadTree->Branch("cons_lead_Eta", &cons_lead_Eta); cons_leadPhi = cons_leadTree->Branch("cons_lead_Phi", &cons_lead_Phi);
  cons_leadM = cons_leadTree->Branch("cons_lead_M", &cons_lead_M); cons_leadE = cons_leadTree->Branch("cons_lead_E", &cons_lead_E);

  TBranch *inclPt;  TBranch *inclEta;  TBranch *inclPhi; TBranch *inclM; TBranch *inclE;
  inclPt = inclTree->Branch("incl_Pt", &incl_Pt);  inclEta = inclTree->Branch("incl_Eta", &incl_Eta);  inclPhi = inclTree->Branch("incl_Phi", &incl_Phi);
  inclM = inclTree->Branch("incl_M", &incl_M);  inclE = inclTree->Branch("incl_E", &incl_E);

  TBranch *cons_inclPt;  TBranch *cons_inclEta;  TBranch *cons_inclPhi; TBranch *cons_inclM; TBranch *cons_inclE;
  cons_inclPt = cons_inclTree->Branch("cons_incl_Pt", &cons_incl_Pt);  cons_inclEta = cons_inclTree->Branch("cons_incl_Eta", &cons_incl_Eta);  cons_inclPhi = cons_inclTree->Branch("cons_incl_Phi", &cons_incl_Phi);
  cons_inclM = cons_inclTree->Branch("cons_incl_M", &cons_incl_M); cons_inclE = cons_inclTree->Branch("cons_incl_E", &cons_incl_E);
  
  // Helpers
  // -------
  vector<PseudoJet> particles;

  int nJets = 0;

  // Constituent selectors
  // ---------------------
  Selector select_track_rap = fastjet::SelectorAbsRapMax(max_track_rap);
  Selector select_lopt      = fastjet::SelectorPtMin( partMinPt );
  
  Selector slo = select_track_rap * select_lopt;

  // Jet candidate selectors
  // -----------------------
  Selector select_jet_rap     = fastjet::SelectorAbsRapMax(max_rap);
  Selector select_jet_pt_min  = fastjet::SelectorPtMin( jet_ptmin );
  Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
  Selector sjet = select_jet_rap && select_jet_pt_min && select_jet_pt_max;
  
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

  cout << "Performing analysis." << endl;
  // Cycle through events
  // --------------------  
  int nEvents = 10000;
  int nEventsUsed = 0;
	
  //initialize the reader
  reader.Init( nEvents ); //runs through all events with -1
  
  try{
    while ( reader.NextEvent() ) {
      reader.PrintStatus(10);
      //get the event header
      header = reader.GetEvent()->GetHeader();
      particles.clear();
      
      // Get the output container from the reader
      // ----------------------------------------
      container = reader.GetOutputContainer();
      

    	// Transform TStarJetVectors into (FastJet) PseudoJets
    	// ----------------------------------------------------------
    	for ( int i=0; i < container->GetEntries() ; ++i ){
            sv = container->Get(i);
	    fastjet::PseudoJet current = fastjet::PseudoJet( *sv );
	    current.set_user_index(sv->GetCharge());
	    if (sv->GetCharge() != 0) { //for charged particles, we assign the pion mass
	      current.reset_PtYPhiM(sqrt(current.perp2()),current.rap(),current.phi(), PionMass);
	    }	    
            particles.push_back( current );
    	}
    
    	// Analysis
    	// --------
    	// Apply selector to the full particle set
    	vector<PseudoJet> pLo = slo( particles );        

	    // find corresponding jets with soft constituents
  	  // ----------------------------------------------
    	ClusterSequenceArea csaLo ( pLo, jet_def, area_def ); // WITH background subtraction

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
	bkgd_subtractor.set_use_rho_m();
    	vector<PseudoJet> LoResult = fastjet::sorted_by_pt(sjet(bkgd_subtractor(csaLo.inclusive_jets())));
	
	if (LoResult.size() != 0) {
	  ++ nJets;
	  //leading
	  m_v_pt_lead->Fill(LoResult[0].m(), LoResult[0].pt());
	  PtEtaPhi_lead->Fill(LoResult[0].pt(),LoResult[0].eta(),LoResult[0].phi());
	  vector<PseudoJet> LoLead; LoLead.push_back(LoResult[0]); 
	  FillTrees(LoLead, leadTree, lead_Pt, lead_Eta, lead_Phi, lead_E, lead_M);
	  FillTrees(LoLead[0].constituents(), cons_leadTree, cons_lead_Pt, cons_lead_Eta, cons_lead_Phi, cons_lead_E, cons_lead_M);
	  for (int cons = 0; cons < LoResult[0].constituents().size(); ++ cons) {
	    if (LoResult[0].constituents()[cons].pt() < partMinPt) {continue;} //ignores contributions from ghosts
	    cons_PtEtaPhi_lead->Fill(LoResult[0].constituents()[cons].pt(), LoResult[0].constituents()[cons].eta(), LoResult[0].constituents()[cons].phi());
	  }
	  m_lead->Fill(LoResult[0].m());
	  
	  //inclusive
	  FillTrees(LoResult, inclTree, incl_Pt, incl_Eta, incl_Phi, incl_E, incl_M);
	  for (int j = 0; j < LoResult.size(); ++ j) {
	    m_v_pt_incl->Fill(LoResult[j].m(), LoResult[j].pt());
	    PtEtaPhi_incl->Fill(LoResult[j].pt(),LoResult[j].eta(),LoResult[j].phi()); 
	     FillTrees(LoResult[j].constituents(), cons_inclTree, cons_incl_Pt, cons_incl_Eta, cons_incl_Phi, cons_incl_E, cons_incl_M);
	    for (int cons = 0; cons < LoResult[j].constituents().size(); ++ cons) {
	      if (LoResult[j].constituents()[cons].pt() < partMinPt) {continue;} //ignores contributions from ghosts
	      cons_PtEtaPhi_incl->Fill(LoResult[j].constituents()[cons].pt(), LoResult[j].constituents()[cons].eta(), LoResult[j].constituents()[cons].phi());
	    }
	    m_incl->Fill(LoResult[j].m());
	  }
	}
	  // And we're done! Fill histograms
	  // -----------------------------
    
	nEventsUsed++;
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
 
  leadTree->Write("lead"); cons_leadTree->Write("cons_lead"); 
  inclTree->Write("incl"); cons_inclTree->Write("cons_incl");

  PtEtaPhi_incl->Write(); PtEtaPhi_lead->Write();
  cons_PtEtaPhi_incl->Write(); cons_PtEtaPhi_lead->Write();
  m_incl->Write(); m_lead->Write(); m_v_pt_incl->Write(); m_v_pt_lead->Write();

  //fout->Write();
  fout->Close();

  cout << "In " << nEventsUsed << " events, found " << endl
       << nJets << " leading jets with constituents above 0.2 GeV." << endl;  

  cout << "Wrote to " << fout->GetName() << endl;
  cout << "Bye :)" << endl;
    
  return 0;
}
  
