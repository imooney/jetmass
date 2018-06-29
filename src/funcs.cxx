//  functions.cxx
//  Veronica Verkest May 13, 2018
//  Adapted by Isaac Mooney June, 2018

#include "params.hh"
#include "funcs.hh"
//#include "ktTrackEff.hh"

namespace Analysis {

  // -------------------------                                                                                                                                                                            
  // IO/OS Manip functionality                                                                                                                                                                            
  // -------------------------                                                                                                                                                                             
  // Used to understand which format of input file is being used                                                                                                                                          
  // ( .root file, .txt, .list, etc )                                                                                                                                                                     
  // ---------------------------------------------------------------------                                                                                                                                 
  bool HasEnding (std::string const &full_string, std::string const &ending) {
    if (full_string.length() >= ending.length()) {
      return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending) );
    } else {
      return false;
    }
  }

  void FillTrees ( std::vector<fastjet::PseudoJet> jets, TTree* Tree, double &jPt, double &jEta, double &jPhi, double &jM, double &jE, int &jncons, double &wt, double weight) {
    for ( int j = 0; j< jets.size(); ++ j) {   // FILL JET INFO
      int nGhosts = 0;
      if (jets[j].pt() < 0.2) continue;
      jPt = jets[j].pt();    jEta = jets[j].eta();    jPhi = jets[j].phi();
      jE = jets[j].e();    jM = jets[j].m(); wt = weight;
      std::vector<fastjet::PseudoJet> Cons = jets[j].constituents(); //jncons = Cons.size();
      //cons without ghosts:
      for (int c = 0; c < Cons.size(); ++ c) {
	if (Cons[c].pt() < 0.2) ++nGhosts;
      }
      jncons = Cons.size() - nGhosts;
      Tree->Fill();
    }
  }

  void AnalysisSummary( int events, int pJets, int eJets, int gJets, int pgMatchedJets, int epMatchedJets, int egMatchedJets, std::string outName ) {
    std::cout << std::endl << std::endl << " Of " << events << " events: "<< std::endl;
    std::cout << pJets << " jets have been found for the Pythia6 data" << std::endl;
    std::cout << eJets << " jets have been found for the Pythia6 + Efficiency data" << std::endl;
    std::cout << gJets << " jets have been found for the Pythia6 + GEANT data" << std::endl << std::endl;
    std::cout << pgMatchedJets << " GEANT jets have been matched to Pythia6 jets" << std::endl;
    std::cout << epMatchedJets << " Pythia6+Efficiency jets have been matched to Pythia6 jets" << std::endl;
    std::cout << egMatchedJets << " GEANT jets have been matched to Pythia6+Efficiency jets" << std::endl;
    std::cout <<std::endl << "Writing to:  " << outName << std::endl << std::endl;
  }

  void GatherParticles ( TStarJetVectorContainer<TStarJetVector> * container , TStarJetVector *sv, std::vector<fastjet::PseudoJet> & Particles, const bool full){
    for ( int i=0; i < container->GetEntries() ; ++i ) {
      sv = container->Get(i);
      fastjet::PseudoJet current = fastjet::PseudoJet( *sv );
      current.set_user_index( sv->GetCharge() );
      
      if (sv->GetCharge() != 0) {
	current.reset_PtYPhiM(sqrt(current.perp2()),current.rap(),current.phi(), PionMass); //assigning pion mass to charged particles
      }
      if ((sv->GetCharge() == 0) && (full == 0)) { continue; } // if we don't want full jets, skip neutrals

      Particles.push_back(current);
    }
    return;
  }

  double LookupXsec(TString currentfile ) {

    static const Double_t Xsec[12] = {
      1.0,        // Placeholder for 2-3
      1.30E+09, // 3-4
      3.15E+08, // 4-5
      1.37E+08, // 5-7
      2.30E+07, // 7-9
      5.53E+06, // 9-11
      2.22E+06, // 11-15
      3.90E+05, // 15-25
      1.02E+04, // 25-35
      5.01E+02, // 35-45
      2.86E+01, // 45-55
      1.46E+00 // 55-65
    };

    static const Double_t Nmc[12] = {
      1, // 2-3
      672518, // 3-4
      672447, // 4-5
      393498, // 5-7
      417659, // 7-9
      412652, // 9-11
      419030, // 11-15
      396744, // 15-25
      399919, // 25-35
      119995, // 35-45
      117999, // 45-55
      119999 // 55-65
    };

    Double_t w[12];
    for ( int i=0; i<12 ; ++i ){
      w[i] = Xsec[i] / Nmc[i];
      // w[i] = Nmc[i] / Xsec[i];
    }

    if ( currentfile.Contains("picoDst_3_4") ) return w[1];
    if ( currentfile.Contains("picoDst_4_5") ) return w[2];
    if ( currentfile.Contains("picoDst_5_7") ) return w[3];
    if ( currentfile.Contains("picoDst_7_9") ) return w[4];
    if ( currentfile.Contains("picoDst_9_11") ) return w[5];
    if ( currentfile.Contains("picoDst_11_15") ) return w[6];
    if ( currentfile.Contains("picoDst_15_25") ) return w[7];
    if ( currentfile.Contains("picoDst_25_35") ) return w[8];
    if ( currentfile.Contains("picoDst_35_45") ) return w[9];
    if ( currentfile.Contains("picoDst_45_55") ) return w[10];
    if ( currentfile.Contains("picoDst_55_65") ) return w[11];
    return 1;
  }

  //finds potential trigger jets out of the two highest pT jets
  bool GetTriggerJet(std::vector<fastjet::PseudoJet> & triggers, const std::vector<fastjet::PseudoJet> jets) {
    triggers.clear();
    bool placeholder = 0; //to keep track of which jet was the trigger in the case of only one trigger 
    for (int i = 0; i < jets.size(); ++ i) {
      if (i == 2) {return placeholder;} //only want to look at the two highest pT jets
      for (int j = 0; j < jets[i].constituents().size(); ++ j) {
	if (jets[i].constituents()[j].pt() > dat_evEtMin) {//has a trigger
	  placeholder = i;
	  triggers.push_back(jets[i]);
	  break;
	}
      }
      if (triggers.size() == 2) {return placeholder;} //we only need at most two objects in the triggers vector - one to be trigger, one to be recoil
    }
    return placeholder;
  }

  
  //  INITIATE READER
  void InitReader( TStarJetPicoReader & reader, TChain* chain, int nEvents, const std::string trig, const double vZ, const double vZDiff, const double Pt, const double Et, const double Etmin,  const double DCA, const double NFit, const double NFitRatio, const double maxEtTow, const std::string badTows) {
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader.SetApplyFractionHadronicCorrection( true );
    reader.SetFractionHadronicCorrection( 0.9999 );
    reader.SetRejectTowerElectrons( kFALSE );
    
    // Event and track selection
    // -------------------------
    
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetTriggerSelection( trig.c_str() ); //All, MB, HT, pp, ppHT, ppJP
    evCuts->SetVertexZCut ( vZ );
    evCuts->SetVertexZDiffCut( vZDiff );
    evCuts->SetRefMultCut( refMultCut );
    evCuts->SetMaxEventPtCut( Pt );
    evCuts->SetMaxEventEtCut( Et );
    evCuts->SetMinEventEtCut( Etmin );

    // Tracks cuts
    TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
    trackCuts->SetDCACut( DCA );
    trackCuts->SetMinNFitPointsCut( NFit );
    trackCuts->SetFitOverMaxPointsCut( NFitRatio );    
    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
    // Towers
    TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
    towerCuts->SetMaxEtCut( maxEtTow );
    towerCuts->AddBadTowers( badTows );

    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader.SetProcessV0s(false);
    
    // Initialize the reader
    reader.Init( nEvents ); //runs through all events with -1
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FILL HISTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  void FillHistsHelper(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> &c2D, Collection<std::string, TH3D> & c3D, const std::string flag1, const std::string flag2, const std::string flag3, const fastjet::PseudoJet jet, const double weight) {
    c1D.fill(("m_" + flag1 + "_" + flag3 + "_" + "jet" + flag2).c_str(), jet.m(), weight);
    c2D.fill(("m_v_pt_" + flag1 + "_" + flag3 + "_" + "jet" + flag2).c_str(), jet.m(), jet.pt(), weight);
    c3D.fill(("PtEtaPhi_" + flag1 + "_" + flag3 + "_" + "jet" + flag2).c_str(), jet.pt(), jet.eta(), jet.phi(), weight);
    for (int cons = 0; cons < jet.constituents().size(); ++ cons) {
      if (jet.constituents()[cons].pt() < partMinPt) {continue;} //ignores contributions from ghosts                       
      c3D.fill(("PtEtaPhi_" + flag1 + "_" + flag3 + "_" + "cons" + flag2).c_str(), jet.constituents()[cons].pt(), jet.constituents()[cons].eta(), jet.constituents()[cons].phi(), weight);
    }
    return;
  }
  
  void FillHists(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> & c2D, Collection<std::string, TH3D> & c3D, const std::string flag1, const std::string flag2, const std::vector<fastjet::PseudoJet> jets, const double weight) {
    //leading
    FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "lead", jets[0], weight);
    //subleading
    if (jets.size() > 1) {
      FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "sublead", jets[1], weight);
    }
    //inclusive
    for (int i = 0; i < jets.size(); ++ i) {
      FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "incl", jets[i], weight);
    }
    //trigger & recoil
    std::vector<fastjet::PseudoJet> candidates;
    bool which_one = GetTriggerJet(candidates, jets);
    if (candidates.size() == 0 || jets.size() < 2) { // means there isn't a trigger or there isn't a recoil
      return;
    }
    if (candidates.size() == 1 && jets.size() > 1) { //potential trigger
      if (fabs(fabs(jets[which_one].delta_phi_to(jets[(which_one + 1) % 2])) - Pi) < R) { //found a recoil
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "trig", jets[which_one], weight); //filling hists for trigger
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "rec", jets[(which_one + 1) % 2], weight); //filling hists for recoil
	return;
      }
    }
    if (candidates.size() == 2) {
      if (fabs(fabs(candidates[0].delta_phi_to(candidates[1])) - Pi) < R) { //trigger & recoil found!
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "trig", candidates[0], weight); //filling hists for trigger
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "rec", candidates[1], weight); //filling hists for recoil
	return;
      }
    }
    return;
  }
}
