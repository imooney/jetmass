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
  InitReader(reader, chain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, dat_maxEtTow, /*"dummy_badtows.list"*/det_badTowers, dat_bad_run_list);

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

  Collection<string, TH1D> hists1D; Collection<string, TH2D> hists2D; Collection<string, TH3D> hists3D;
  Collection<string, THnSparseD> histsnD;

  vector<string> flag_i = {"ch", "full"};
  vector<string> flag_j = {"lead", "sublead", "trig", "rec", "incl"};   
  vector<string> flag_k = {"jet", "cons", "sd"};
  for (int i = 0; i < flag_i.size(); ++ i) {
    for (int j = 0; j < flag_j.size(); ++ j) {
      for (int k = 0; k < flag_k.size(); ++ k) {
	hists1D.add(("m_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]).c_str(),"",20,0,10); //mass
	hists2D.add(("m_v_pt_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]).c_str(),";M;p_{T}",20,0,10,11,5,60); //mass vs. pT
	hists2D.add(("m_v_pt_rebin_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]).c_str(),";M;p_{T}",20,0,10,9,15,60);
	hists3D.add(("PtEtaPhi_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]).c_str(),"",11,5,60,30,-0.6,0.6,50,0,2*Pi); //pT vs. eta vs. phi
      }
    }
  }
  const unsigned nDim = 5;
  int bins[nDim] = {20, 20, 20, 11, 11};
  double min[nDim] = {0,0,0,5,5};
  double max[nDim] = {1,10,1,60,60};
  THnSparse * SDnD = new THnSparseD("zg_mg_thetag_ptg_pt_full_incl_sd", "", nDim, bins, min, max);
  SDnD->Sumw2();
  
  TH1D * dPhi_trig_rec = new TH1D("dPhi_trig_rec",";#Delta #phi;arb.", 28, -Pi - 0.4, Pi + 0.4); //defined as trigger - recoil
  //TEST
  TH3D *PtEtaPhi_tracks = new TH3D("PtEtaPhi_tracks",";p^{track}_{T} [GeV/c]; #eta; #phi",80,0,80,43,-1,1,50,0,2*Pi);
  TH2D *nCons_v_pt = new TH2D("nCons_v_pt",";num. constituents; p_{T}^{jet} [GeV/c]",20,0,20,80,0,80);
  TH2D *ch_frac_v_pt = new TH2D("ch_frac_v_pt",";charged frac.;p_{T}^{jet} [GeV/c]",10,0,1,11,5,60);
  TH1D *towers_check = new TH1D("towers_check","",200,0,200);
  TH2D *tow_id_v_e = new TH2D("tow_id_v_e",";tower ID;tower E_{T} [GeV]",4800,1,4801,140,0,140);
  TH1D *tow_freq = new TH1D("tow_freq",";tower ID;tower frequency",4800,1,4801);
  TH2D *tow_eta_phi_e_w = new TH2D("tow_eta_phi_e_w",";tower #eta;tower #phi", 40,-1,1,126,-Pi,Pi);
  TH2D *tow_eta_phi_e_wo = new TH2D("tow_eta_phi_e_wo",";tower #eta;tower #phi",40,-1,1,126,-Pi,Pi);
  
  // Trees
  // -----
  double lead_Pt, lead_Eta, lead_Phi, lead_M, lead_E;
  double sublead_Pt, sublead_Eta, sublead_Phi, sublead_M, sublead_E;
  double cons_lead_Pt, cons_lead_Eta, cons_lead_Phi, cons_lead_M, cons_lead_E;
  double incl_Pt, incl_Eta, incl_Phi, incl_M, incl_E;
  double cons_incl_Pt, cons_incl_Eta, cons_incl_Phi, cons_incl_M, cons_incl_E;
  int nCons_incl, nCons_lead, nCons_sublead, dummy_int;
  double dummy_double;

  TTree *leadTree = new TTree("lead","lead");
  TTree *subleadTree = new TTree("sublead","sublead");
  TTree *cons_leadTree = new TTree("cons_lead","cons_lead");
  TTree *inclTree = new TTree("incl","incl");
  TTree *cons_inclTree = new TTree("cons_incl","cons_incl");

  //TBranch *leadPt;  TBranch *leadEta;  TBranch *leadPhi; TBranch *leadM; TBranch *leadE; TBranch *NCons_lead;
  leadTree->Branch("Pt", &lead_Pt); leadTree->Branch("Eta", &lead_Eta);  leadTree->Branch("Phi", &lead_Phi);
  leadTree->Branch("M", &lead_M); leadTree->Branch("E", &lead_E); leadTree->Branch("nCons", &nCons_lead);

  //TBranch *cons_leadPt;  TBranch *cons_leadEta;  TBranch *cons_leadPhi; TBranch *cons_leadM; TBranch *cons_leadE;
  cons_leadTree->Branch("Pt", &cons_lead_Pt); cons_leadTree->Branch("Eta", &cons_lead_Eta); cons_leadTree->Branch("Phi", &cons_lead_Phi);
  cons_leadTree->Branch("M", &cons_lead_M); cons_leadTree->Branch("E", &cons_lead_E);

  //TBranch *subleadPt;  TBranch *subleadEta;  TBranch *subleadPhi; TBranch *subleadM; TBranch *subleadE; TBranch *NCons_sublead;
  subleadTree->Branch("Pt", &sublead_Pt); subleadTree->Branch("Eta", &sublead_Eta); subleadTree->Branch("Phi", &sublead_Phi);
  subleadTree->Branch("M", &sublead_M); subleadTree->Branch("E", &sublead_E); subleadTree->Branch("nCons", &nCons_sublead);

  //TBranch *inclPt;  TBranch *inclEta;  TBranch *inclPhi; TBranch *inclM; TBranch *inclE; TBranch *NCons_incl;
  inclTree->Branch("Pt", &incl_Pt); inclTree->Branch("Eta", &incl_Eta); inclTree->Branch("Phi", &incl_Phi);
  inclTree->Branch("M", &incl_M); inclTree->Branch("E", &incl_E); inclTree->Branch("nCons", &nCons_incl);

  //TBranch *cons_inclPt;  TBranch *cons_inclEta;  TBranch *cons_inclPhi; TBranch *cons_inclM; TBranch *cons_inclE;
  cons_inclTree->Branch("Pt", &cons_incl_Pt); cons_inclTree->Branch("Eta", &cons_incl_Eta); cons_inclTree->Branch("Phi", &cons_incl_Phi);
  cons_inclTree->Branch("M", &cons_incl_M); cons_inclTree->Branch("E", &cons_incl_E);
  
  // Helpers
  // -------
  vector<PseudoJet> particles, ch_particles;

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

  //Creating SoftDrop grooming object
  contrib::SoftDrop sd(beta,z_cut,R0);
  cout << "SoftDrop groomer is: " << sd.description() << endl;
  
  cout << "Performing analysis." << endl;
  // Cycle through events
  // --------------------  
  //int nEvents = -1;
  int nEventsUsed = 0;
	
  //initialize the reader
  //  reader.Init( nEvents ); //runs through all events with -1
 
  
  try{
    while ( reader.NextEvent() ) {
      //auto start_timetotal = clock::now();

      //TEMP!!!!!!!!!!!!!!!!!
      //METHOD OF GETTING RAW TOWER INFO FROM EVENT
      TStarJetPicoEvent * current_event = reader.GetEvent();
      //      cout << "event " << reader.GetNOfCurrentEvent() << " has these tows" << endl;
      
      for (int i = 0; i < current_event->GetTowers()->GetEntries(); ++ i) {
	if (current_event->GetTower(i)->GetEnergy() > 0) {
	  tow_freq->Fill(current_event->GetTower(i)->GetId());
	}
	tow_id_v_e->Fill(current_event->GetTower(i)->GetId(), current_event->GetTower(i)->GetEnergy());
	//cout << current_event->GetTower(i)->GetId() << " " << current_event->GetTower(i)->GetEnergy() << endl;
	tow_eta_phi_e_wo->Fill(current_event->GetTower(i)->GetEta(), current_event->GetTower(i)->GetPhi(), current_event->GetTower(i)->GetEt());
      }
      
      //METHOD OF GETTING TOWER INFO FROM LIST OF SELECTED TOWERS (BUT NOT ALL SELECTIONS APPLIED YET)
      
      TList* select_tows = reader.GetListOfSelectedTowers();
      TIter next(select_tows);
      TStarJetPicoTower* object = 0;
      while ((object = (TStarJetPicoTower*) next())){
	//tow_freq->Fill(object->GetId());
	//tow_id_v_e->Fill(object->GetId(), object->GetEnergy());
	tow_eta_phi_e_w->Fill(object->GetEta(), object->GetPhi(), object->GetEt());
      }
      
      //METHOD OF GETTING TOWER INFO FROM CONTAINER AFTER ALL SELECTIONS ARE APPLIED
      /*
      container = reader.GetOutputContainer();
      for (int i = 0; i < container->GetEntries(); ++ i) {
	if (container->Get(i)->GetCharge() == 0) {
	  tow_freq->Fill(container->Get(i)->GetTowerID());
	  tow_id_v_e->Fill(container->Get(i)->GetTowerID(),container->Get(i)->Et());
	}
      }
      */
      
      
      reader.PrintStatus(10);
      //get the event header
      
      header = reader.GetEvent()->GetHeader();
      particles.clear(); ch_particles.clear();
      
      // Get the output container from the reader
      // ----------------------------------------
      container = reader.GetOutputContainer();
      
      // Transform TStarJetVectors into (FastJet) PseudoJets
      // ----------------------------------------------------------
      GatherParticles(container, sv, particles, 1,0); //ch+ne
      GatherParticles(container, sv, ch_particles, 0,0); //ch
      
      for(int i = 0; i < particles.size(); ++ i) {
	PtEtaPhi_tracks->Fill(particles[i].pt(), particles[i].eta(), particles[i].phi());
      }
             
      //auto start_timefunc = clock::now();
      
      // Analysis
      // --------
      // Apply selector to the full particle set
      vector<PseudoJet> pLo = slo( particles ); vector<PseudoJet> ch_pLo = slo( ch_particles );
      
      // find corresponding jets with soft constituents
      // ----------------------------------------------
      ClusterSequence/*Area*/ csaLo ( pLo, jet_def/*, area_def */); // WITHOUT background subtraction
      ClusterSequence/*Area*/ ch_csaLo ( ch_pLo, jet_def/*, area_def */); //WITHOUT background subtraction
      
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
      for (int i = 0; i < LoInitial.size(); ++ i) {
        double towersum = 0; double ptsum = 0;
        for (int j = 0; j < LoInitial[i].constituents().size(); ++ j) {
          if (LoInitial[i].constituents()[j].user_index() == 0) {
            towersum += LoInitial[i].constituents()[j].pt();
          }
          ptsum += LoInitial[i].constituents()[j].pt();
        }
        if (towersum / (double) ptsum < NEF_max) {
          LoResult.push_back(LoInitial[i]);
        }
      } 
      vector<PseudoJet> GroomedJets;
      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
      for (int i = 0; i < LoResult.size(); ++ i) {
	GroomedJets.push_back(sd(LoResult[i]));
      }

      vector<PseudoJet> ch_LoResult = fastjet::sorted_by_pt(sjet(ch_csaLo.inclusive_jets()));
      vector<PseudoJet> ch_GroomedJets;
      for (int i = 0; i < ch_LoResult.size(); ++ i) {
	ch_GroomedJets.push_back(sd(ch_LoResult[i]));
      }      
      //      auto execution_func = std::chrono::duration_cast<std::chrono::microseconds>(clock::now() - start_timefunc).count();
      //cout << "Function time = " << execution_func << " microseconds" << endl;
            
      if (LoResult.size() != 0) {
	//SoftDrop is a groomer not a tagger, so if we have at least one ungroomed jet, we should also have a SoftDrop'd jet.
	++ nJets;

	//TEMP!!!!
	for (int i = 0; i < LoResult.size(); ++ i) {
	  if (LoResult[i].pt() > 30) {
	    for (int j = 0; j < LoResult[i].constituents().size(); ++ j) {
	      towers_check->Fill(LoResult[i].constituents()[j].Et());
	    }  
	  }
	}
	FillHists(hists1D, hists2D, hists3D, "full", "", LoResult, 1.0); 
	FillSDHists(hists1D, hists2D, hists3D, "full", "", GroomedJets, 1.0);

	//leading
	vector<PseudoJet> LoLead; LoLead.push_back(LoResult[0]); //I do this so I can use the same function for jets & constituents (takes a vector of pseudojets)
	FillTrees(LoLead, leadTree, lead_Pt, lead_Eta, lead_Phi, lead_M, lead_E, nCons_lead, dummy_double, dummy_double);
	FillTrees(LoLead[0].constituents(), cons_leadTree, cons_lead_Pt, cons_lead_Eta, cons_lead_Phi, cons_lead_M, cons_lead_E, dummy_int, dummy_double, dummy_double);
	//subleading
	if (LoResult.size() > 1) {
	  vector<PseudoJet> LoSublead; LoSublead.push_back(LoResult[1]);
	  FillTrees(LoSublead, subleadTree, sublead_Pt, sublead_Eta, sublead_Phi, sublead_M, sublead_E, nCons_sublead, dummy_double, dummy_double);
	}
	//inclusive
	FillTrees(LoResult, inclTree, incl_Pt, incl_Eta, incl_Phi, incl_M, incl_E, nCons_incl, dummy_double, dummy_double);
	for (int j = 0; j < LoResult.size(); ++ j) {
	  double val_list[nDim] = {GroomedJets[j].structure_of<contrib::SoftDrop>().symmetry(),GroomedJets[j].m(),GroomedJets[j].structure_of<contrib::SoftDrop>().delta_R(),GroomedJets[j].pt(), LoResult[j].pt()}; //Groomed jet not guaranteed to be highest pT even though ungroomed one is, but for inclusive this doesn't matter
	  SDnD->Fill(val_list);
	  FillTrees(LoResult[j].constituents(), cons_inclTree, cons_incl_Pt, cons_incl_Eta, cons_incl_Phi, cons_incl_M, cons_incl_E, dummy_int, dummy_double, dummy_double);
	  //TEMP
	  nCons_v_pt->Fill(LoResult[j].constituents().size(), LoResult[j].pt());
	  int numch = 0;
	  for (int k = 0; k < LoResult[j].constituents().size(); ++ k) {
	    if (LoResult[j].constituents()[k].user_index() != 0) {
	      numch ++;
	    }
	  }
	  ch_frac_v_pt->Fill(numch/(double) LoResult[j].constituents().size(),LoResult[j].pt());
	}
      
	
	//TEMP
	//trigger & recoil                                                                                             
	std::vector<fastjet::PseudoJet> candidates;
	bool which_one = GetTriggerJet(candidates, LoResult);
	if (candidates.size() == 1 && LoResult.size() > 1) { //potential trigger     
	  if (fabs(fabs(LoResult[which_one].delta_phi_to(LoResult[(which_one + 1) % 2])) - Pi) < R) { //found a recoil
	    dPhi_trig_rec->Fill(LoResult[which_one].phi() - LoResult[(which_one+1)%2].phi());
	  }
	}
	else if (candidates.size() == 2) {
	  if (fabs(fabs(candidates[0].delta_phi_to(candidates[1])) - Pi) < R) { //trigger & recoil found!              
	    dPhi_trig_rec->Fill(candidates[0].phi() - candidates[1].phi());   
	  }
	}
	
	}
      
      if (ch_LoResult.size() != 0) {
	FillHists(hists1D, hists2D, hists3D, "ch", "", ch_LoResult, 1.0);
	FillSDHists(hists1D, hists2D, hists3D, "ch","",ch_GroomedJets, 1.0);
      }
      
      // And we're done!
      // -----------------------------
      
      nEventsUsed++;
      
	
      //      auto execution_total = std::chrono::duration_cast<std::chrono::microseconds>(clock::now() - start_timetotal).count();
      //      cout << "Total time = " << execution_total << " microseconds" << endl;
      
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
  subleadTree->Write("sublead");
  inclTree->Write("incl"); cons_inclTree->Write("cons_incl");
    
  for (int i = 0; i < flag_i.size(); ++ i) {
    for (int j = 0; j < flag_j.size(); ++ j) {
      for (int k = 0; k < flag_k.size(); ++ k) {
	hists1D.write(("m_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]).c_str());
	hists2D.write(("m_v_pt_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]).c_str());
	hists2D.write(("m_v_pt_rebin_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]).c_str());
	hists3D.write(("PtEtaPhi_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]).c_str());
      }
    }
  }
  SDnD->Write();

  dPhi_trig_rec->Write(); PtEtaPhi_tracks->Write(); nCons_v_pt->Write(); ch_frac_v_pt->Write();
  towers_check->Write(); tow_id_v_e->Write(); tow_freq->Write(); tow_eta_phi_e_w->Write(); tow_eta_phi_e_wo->Write();

  
  //fout->Write();
  fout->Close();

  cout << "In " << nEventsUsed << " events, found " << endl
       << nJets << " leading jets with constituents above 0.2 GeV." << endl;  

  cout << "Wrote to " << fout->GetName() << endl;
  cout << "Bye :)" << endl;
    
  return 0;
}
  
