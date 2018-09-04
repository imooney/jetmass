//  Isaac Mooney 8/28/2018 - for jet mass analysis

#include "params.hh"
#include "funcs.hh"
#include "RooUnfold.h"

#include "fastjet/contrib/SoftDrop.hh"

using namespace fastjet;
using namespace std;
using namespace Analysis;
typedef fastjet::contrib::SoftDrop SD;

// -------------------------
// Command line arguments: ( Defaults 
// Defined for debugging in main )           
// [0]: output directory
// [1]: name for the histogram file                                   
// [2]: input data: can be a single .root or a .txt or .list of root files - should always be last argument  

int main (int argc, const char ** argv) {
  
  //  TStarJetPicoDefinitions::SetDebugLevel(0);
  TH1::SetDefaultSumw2( );  // Histograms will calculate gaussian errors
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );

  // Read in command line arguments       
  // ------------------------------                                                                                                                                                                       
  // Defaults 
  std::string        outputDir        = "out/";                        // directory where everything will be saved
  std::string     outFileName        = "test.root"; 
  std::string         chainList            = "simlist.txt";            // input file: can be .root, .txt, .list
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
    //    cout << outputDir << " " << outFileName << " " << full << " " << arguments[3] << " " << chainList << endl;
    break;
  }
  default: { // Error: invalid custom settings                                                                                                                                                             
    __ERR("Invalid number of command line arguments");
    return -1;
    break;
  }
  }



  //  Initialize readers and provide chains
  TStarJetPicoReader P6Reader;      TChain* P6Chain = new TChain( "JetTreeMc" );    //  PURE PYTHIA DATA  (particle)
  TStarJetPicoReader GEANTReader;   TChain* GEANTChain = new TChain( "JetTree" );     //  CORRESPONDING GEANT DATA  (detector)
  
  // Check to see if the input is a .root file or a .txt                                                                                                                                                   
  bool inputIsRoot = Analysis::HasEnding( chainList.c_str(), ".root" );
  bool inputIsTxt  = Analysis::HasEnding( chainList.c_str(), ".txt"  );
  bool inputIsList = Analysis::HasEnding( chainList.c_str(), ".list" );

  // If its a recognized file type, build the chain                                                                                                                                                       
  // If its not recognized, exit                                                                                                                                                                           
  if ( inputIsRoot ) { P6Chain->Add( chainList.c_str()); GEANTChain->Add( chainList.c_str());}
  else if ( inputIsTxt )  { P6Chain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str()); GEANTChain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str());}
  else if ( inputIsList)  { P6Chain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str()); GEANTChain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str());}
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }

  //initialize the readers!
  InitReader(P6Reader, P6Chain, nEvents, truth_triggerString, truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, sim_badTowers, sim_bad_run_list);
  InitReader(GEANTReader, GEANTChain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, sim_badTowers, sim_bad_run_list);
  
  TStarJetPicoEventHeader* p_header;  TStarJetPicoEvent* p_event;  TStarJetPicoEventHeader* g_header;  TStarJetPicoEvent* g_event;
  TStarJetVectorContainer<TStarJetVector> * p_container;         TStarJetVector* p_sv;
  TStarJetVectorContainer<TStarJetVector> * g_container;        TStarJetVector* g_sv;

  //note: there is only one match per event, so none of these vectors should have more than one entry. It is only done for later convenience.
  //  vector<double> deltaPt; vector<double> deltaM; vector<double> deltaZg; vector<double> deltaRg;
  // vector<double> ratioPt; vector<double> ratioM; vector<double> ratioZg; vector<double> ratioRg;
  vector<double> pyPt; vector<double> pyM; vector<double> pyZg; vector<double> pyRg;
  vector<double> gePt; vector<double> geM; vector<double> geZg; vector<double> geRg;
  vector<double> pyPtg; vector<double> pyMg; vector<double> gePtg; vector<double> geMg;
  double mc_weight; int p_EventID;
  
  TTree *eventTree = new TTree("event", "event");
  //eventTree->Branch("deltaPt", &deltaPt); eventTree->Branch("deltaM", &deltaM); eventTree->Branch("deltaZg", &deltaZg); eventTree->Branch("deltaRg", &deltaRg);
  //eventTree->Branch("ratioPt", &ratioPt); eventTree->Branch("ratioM", &ratioM); eventTree->Branch("ratioZg", &ratioZg); eventTree->Branch("ratioRg", &ratioRg);
  eventTree->Branch("pyPt", &pyPt); eventTree->Branch("pyM", &pyM); eventTree->Branch("pyZg", &pyZg); eventTree->Branch("pyRg", &pyRg);
  eventTree->Branch("gePt", &gePt); eventTree->Branch("geM", &geM); eventTree->Branch("geZg", &geZg); eventTree->Branch("geRg", &geRg);
  eventTree->Branch("pyPtg", &pyPtg); eventTree->Branch("pyMg", &pyMg); eventTree->Branch("gePtg", &gePtg); eventTree->Branch("geMg", &geMg);
  eventTree->Branch("weight", &mc_weight); eventTree->Branch("EventID", &p_EventID);

  //temp hists
    TH1D *pt_gen_odd = new TH1D("pt_gen_odd","",15,5,80); TH1D *pt_det_odd = new TH1D("pt_det_odd","",9,15,60);
  
  //responses for MC Closure test
  RooUnfoldResponse pt_odd(15,5,80,15,5,80,"pt_odd","");
  RooUnfoldResponse pt_even(15,5,80,15,5,80,"pt_even","");
  RooUnfoldResponse m_odd(20,0,10,20,0,10,"m_odd","");
  RooUnfoldResponse m_even(20,0,10,20,0,10,"m_even","");
  
  //Creating SoftDrop grooming object                                                                                                                                                        
  contrib::SoftDrop sd(beta,z_cut,R0);
   
  //SELECTORS
  // Constituent selectors                                                                                                                                                     
  // ---------------------                                                                                                                                                        
  Selector select_track_rap = fastjet::SelectorAbsRapMax(max_track_rap);
  Selector select_lopt      = fastjet::SelectorPtMin( partMinPt );
  Selector select_loptmax   = fastjet::SelectorPtMax( partMaxPt );
  Selector spart = select_track_rap * select_lopt * select_loptmax;
  
  // Jet candidate selectors                                                                                                                                                   
  // -----------------------                                                                                                                                                      
  Selector select_jet_rap     = fastjet::SelectorAbsRapMax(max_rap);
  Selector select_det_jet_pt_min  = fastjet::SelectorPtMin( det_jet_ptmin );
  Selector select_gen_jet_pt_min = fastjet::SelectorPtMin( jet_ptmin );
  Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
  Selector sjet_gen = select_jet_rap && select_gen_jet_pt_min && select_jet_pt_max;
  Selector sjet_det = select_jet_rap && select_det_jet_pt_min && select_jet_pt_max;
  
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  TString geantFilename, pythiaFilename;
  
  // Particle containers & counters
  vector<PseudoJet> p_Particles, g_Particles, p_JetsInitial, g_JetsInitial, dummy;
  int nEvents = 0;   int p_NJets = 0;  int g_NJets = 0;  /*int p_EventID; see above for declaration*/   int g_EventID;
  //1=inclusive, 2=lead
  int counter_debug1 = 0, counter_debug2 = 0;
  double p_wt = -1, g_wt = -1; 
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVEN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( GEANTReader.NextEvent() ) {      //    GEANTReader    P6Reader
    //initialize values to -9999
    pyPt.clear(); pyM.clear(); pyZg.clear(); pyRg.clear();
    gePt.clear(); geM.clear(); geZg.clear(); geRg.clear();
    pyPtg.clear(); pyMg.clear(); gePtg.clear(); geMg.clear();
    mc_weight = -9999;
    
    g_EventID = GEANTReader.GetNOfCurrentEvent();
    
    if ( P6Reader.ReadEvent( g_EventID ) != 1 ) continue;   //  ENSURES BOTH DATASETS HAVE AN EVENT
    
    p_Particles.clear(); g_Particles.clear();
    p_JetsInitial.clear(); g_JetsInitial.clear(); //clear all containers
    dummy.clear();
    
    nEvents++;  P6Reader.PrintStatus(10);  GEANTReader.PrintStatus(10);     // Print out reader status every 10 seconds
    
    p_event = P6Reader.GetEvent();       p_header = p_event->GetHeader();           // Get the PYTHIA header and event
    g_event = GEANTReader.GetEvent();    g_header = g_event->GetHeader();           // Get GEANT event header and event
    
    p_EventID = P6Reader.GetNOfCurrentEvent();
    if ( p_EventID != g_EventID ) { cout << endl << "ERROR: READING DIFFERENT EVENTS " <<endl; }
    
    //BEGIN EVENS
    if (p_EventID % 2 == 0) { //even events will be used for the response
      p_container = P6Reader.GetOutputContainer();      // Pythia container
      g_container = GEANTReader.GetOutputContainer();      // GEANT container
      
      pythiaFilename =  P6Reader.GetInputChain()->GetCurrentFile()->GetName();	
      geantFilename =  GEANTReader.GetInputChain()->GetCurrentFile()->GetName();	
      
      //TEMP! CHANGE BACK IF DOESN'T WORK!
      if (((string) pythiaFilename).find("2_3_") != std::string::npos) {continue; }
      
      if (pythiaFilename != geantFilename) {std::cerr << "NOT WHAT I EXPECTED" << std::endl; exit(1);}
      
      p_wt = LookupRun12Xsec( pythiaFilename );
      g_wt = LookupRun12Xsec( geantFilename );
      if (p_wt != g_wt) {std::cerr << "WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; exit(1);}
      mc_weight = p_wt;
      
      //  GATHER PARTICLES
      GatherParticles ( p_container, p_sv, p_Particles, full,1);    //  Pythia particles. full = 0 signifies charged-only, 1 signifies ch+ne
      GatherParticles ( g_container, g_sv, g_Particles, full,0);    //  GEANT particles
      
      vector<PseudoJet> p_cut_Particles = spart(p_Particles); vector<PseudoJet> g_cut_Particles = spart(g_Particles); //applying constituent cuts
      
      ClusterSequence p_Cluster(p_cut_Particles, jet_def); ClusterSequence g_Cluster(g_cut_Particles, jet_def);           //  CLUSTER BOTH
      p_JetsInitial = sorted_by_pt(sjet_gen(p_Cluster.inclusive_jets())); g_JetsInitial = sorted_by_pt(sjet_det(g_Cluster.inclusive_jets()));    // EXTRACT JETS
      vector<PseudoJet> p_Jets; vector<PseudoJet> g_Jets;
      
      //Implementing a neutral energy fraction cut of 90% on inclusive jets
      ApplyNEFSelection(p_JetsInitial, p_Jets); ApplyNEFSelection(g_JetsInitial, g_Jets);
      
      vector<PseudoJet> p_GroomedJets; vector<PseudoJet> g_GroomedJets;
      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)                                                                              
      for (int i = 0; i < p_Jets.size(); ++ i) {
	p_GroomedJets.push_back(sd(p_Jets[i]));
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	g_GroomedJets.push_back(sd(g_Jets[i]));
      }
      /*      
      //TEST!
      bool bad_event = 0;
      std::string tail = ((string) pythiaFilename).substr(((string) pythiaFilename).size() - 10);
      std::string upstring = tail.substr(0,2);
      std::string upstring_copy = upstring;
      if (upstring.find("_") != std::string::npos || upstring.find("-") != std::string::npos) { if (upstring.substr(1,1) != "_") {upstring = upstring.substr(1,1);} else {upstring = \
	    upstring.substr(0,1);}}
      int upbin = std::stoi(upstring);
      
      for (int i = 0; i < p_Jets.size(); ++ i) {
	if ((p_Jets[i].pt() > 2*upbin) && upstring_copy != "-1") {
	  std::cout << "from " << pythiaFilename << " removing pythia event " << p_EventID << " with weight " << mc_weight << " and jet with pt, eta, phi, and m: " << p_Jets[i].pt() << " " << p_Jets[i].eta() << " " << p_Jets[i].phi() << " " << p_Jets[i].m() << std::endl;
	  bad_event = 1;
	}
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	if ((g_Jets[i].pt() > 2*upbin) && upstring_copy != "-1") {
	  std::cout << "from " << pythiaFilename << " removing geant event " << g_EventID << " with weight " << mc_weight << " and jet with pt, eta, phi, and m: " << g_Jets[i].pt() << " " << g_Jets[i].eta() << " " << g_Jets[i].phi() << " " << g_Jets[i].m() << std::endl;
	  bad_event = 1;
	}
      }
      
      if (bad_event == 1) {counter_debug1 ++; continue; }
      if (bad_event == 1) {std::cout << "should never see this message" << std::endl;}
      */
      //constructing even population sample responses for the MC closure test
      if (p_Jets.size() != 0) {
	int position = -1; 
	MatchJets(g_Jets, p_Jets[0], position);
	if (position == -1) { //didn't find a match
          pt_even.Miss(p_Jets[0].pt(), mc_weight);
          m_even.Miss(p_Jets[0].m(), mc_weight);
	}
	else { //found a match
          pt_even.Fill(g_Jets[position].pt(), p_Jets[0].pt(), mc_weight);
          m_even.Fill(g_Jets[position].m(), p_Jets[0].m(), mc_weight);
	}
      }
      
      //fake rate  
      if (g_Jets.size() != 0) {
	int position = -1;
	MatchJets(p_Jets, g_Jets[0], position);
	if (position == -1) { //didn't find a match
          pt_even.Fake(g_Jets[0].pt(), mc_weight);
          m_even.Fake(g_Jets[0].m(), mc_weight);
	}
      }
    } //END EVENS
    
      //BEGIN ODDS
    if (p_EventID % 2 != 0) { //odd events will be used for the 'data'
      p_container = P6Reader.GetOutputContainer();      // Pythia container
      g_container = GEANTReader.GetOutputContainer();      // GEANT container
      
      pythiaFilename =  P6Reader.GetInputChain()->GetCurrentFile()->GetName();
      geantFilename =  GEANTReader.GetInputChain()->GetCurrentFile()->GetName();
      
      //TEMP! CHANGE BACK IF DOESN'T WORK!
      if (((string) pythiaFilename).find("2_3_") != std::string::npos) {continue; }
      
      if (pythiaFilename != geantFilename) {std::cerr << "NOT WHAT I EXPECTED" << std::endl; exit(1);}
      
      p_wt = LookupRun12Xsec( pythiaFilename );
      g_wt = LookupRun12Xsec( geantFilename );
      if (p_wt != g_wt) {std::cerr << "WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; exit(1);}
      mc_weight = p_wt;
      
      //  GATHER PARTICLES
      GatherParticles ( p_container, p_sv, p_Particles, full,1);    //  Pythia particles. full = 0 signifies charged-only, 1 signifies ch+ne
      GatherParticles ( g_container, g_sv, g_Particles, full,0);    //  GEANT particles
      
      vector<PseudoJet> p_cut_Particles = spart(p_Particles); vector<PseudoJet> g_cut_Particles = spart(g_Particles); //applying constituent cuts
      
      ClusterSequence p_Cluster(p_cut_Particles, jet_def); ClusterSequence g_Cluster(g_cut_Particles, jet_def);           //  CLUSTER BOTH
      p_JetsInitial = sorted_by_pt(sjet_gen(p_Cluster.inclusive_jets())); g_JetsInitial = sorted_by_pt(sjet_det(g_Cluster.inclusive_jets()));    // EXTRACT JETS
      vector<PseudoJet> p_Jets; vector<PseudoJet> g_Jets;
      
      //Implementing a neutral energy fraction cut of 90% on inclusive jets
      ApplyNEFSelection(p_JetsInitial, p_Jets); ApplyNEFSelection(g_JetsInitial, g_Jets);
      
      vector<PseudoJet> p_GroomedJets; vector<PseudoJet> g_GroomedJets;
      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
      for (int i = 0; i < p_Jets.size(); ++ i) {
	p_GroomedJets.push_back(sd(p_Jets[i]));
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	g_GroomedJets.push_back(sd(g_Jets[i]));
      }
      /*
      //TEST!
      bool bad_event = 0;
      std::string tail = ((string) pythiaFilename).substr(((string) pythiaFilename).size() - 10);
      std::string upstring = tail.substr(0,2);
      std::string upstring_copy = upstring;
      if (upstring.find("_") != std::string::npos || upstring.find("-") != std::string::npos) { if (upstring.substr(1,1) != "_") {upstring = upstring.substr(1,1);} else {upstring = \
	    upstring.substr(0,1);}}
      int upbin = std::stoi(upstring);
      
      for (int i = 0; i < p_Jets.size(); ++ i) {
	if ((p_Jets[i].pt() > 2*upbin) && upstring_copy != "-1") {
	  std::cout << "from " << pythiaFilename << " removing pythia event " << p_EventID << " with weight " << mc_weight << " and jet with pt, eta, phi, and m: " << p_Jets[i].pt() << " " << p_Jets[i].eta() << " " << p_Jets[i].phi() << " " << p_Jets[i].m() << std::endl;
	  bad_event = 1;
	}
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	if ((g_Jets[i].pt() > 2*upbin) && upstring_copy != "-1") {
	  std::cout << "from " << pythiaFilename << " removing geant event " << g_EventID << " with weight " << mc_weight << " and jet with pt, eta, phi, and m: " << g_Jets[i].pt() << " " << g_Jets[i].eta() << " " << g_Jets[i].phi() << " " << g_Jets[i].m() << std::endl;
	  bad_event = 1;
	}
      }
      
      if (bad_event == 1) {counter_debug1 ++; continue; }
      if (bad_event == 1) {std::cout << "should never see this message" << std::endl;}
      */
      
      //unmatched filling
      for (int i = 0; i < p_Jets.size(); ++ i) {
	pt_gen_odd->Fill(p_Jets[i].pt(), mc_weight);
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	pt_det_odd->Fill(g_Jets[i].pt(), mc_weight);
      }
      
      //constructing odd population sample responses for the MC Closure test
      if (p_Jets.size() != 0) {
	int position = -1; 
	MatchJets(g_Jets, p_Jets[0], position);
	if (position == -1) { //didn't find a match
          pt_odd.Miss(p_Jets[0].pt(), mc_weight);
          m_odd.Miss(p_Jets[0].m(), mc_weight);
	}
	else { //found a match
          pt_odd.Fill(g_Jets[position].pt(), p_Jets[0].pt(), mc_weight);
          m_odd.Fill(g_Jets[position].m(), p_Jets[0].m(), mc_weight);
	}
      }
      
      //fake rate  
      if (g_Jets.size() != 0) {
	int position = -1;
	MatchJets(p_Jets, g_Jets[0], position);
	if (position == -1) { //didn't find a match
          pt_odd.Fake(g_Jets[0].pt(), mc_weight);
          m_odd.Fake(g_Jets[0].m(), mc_weight);
	}
      }
      
      p_NJets += p_Jets.size(); g_NJets += g_Jets.size();               //  Save jet info and add jets to total
    } //END ODDS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
  }
  
  //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  
  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
  
  std::cout << std::endl << std::endl << "Of " << nEvents << " events" << std::endl;
  std::cout << p_NJets << " gen jets have been found" << std::endl;
  std::cout << g_NJets << " det jets have been found" << std::endl << std::endl;
  std::cout <<std::endl << "Writing to:  " << fout->GetName() << std::endl << std::endl;
  std::cout << "Discarded " << counter_debug1 << " events on grounds of the found jets being too much higher than the pT-hat range" << std::endl;

  eventTree->Write("event");

  pt_odd.Write(); pt_even.Write(); m_odd.Write(); m_even.Write();

  pt_gen_odd->Write(); pt_det_odd->Write();
  fout->Close();

  return 0;
}
