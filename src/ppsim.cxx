//  Veronica Verkest        May 13, 2018
//  Compare:   p6  VS  p6+efficiency  VS  p6+GEANT
//  Functions in src/functions.cxx
//  Parameters in src/parameters.cxx
//  Adapted by Isaac Mooney June, 2018 for jet mass analysis

#include "params.hh"
#include "funcs.hh"

using namespace fastjet;
using namespace std;
using namespace Analysis;

// -------------------------                                                                                                                                                                              
// Command line arguments: ( Defaults                                                                                                                                                                     
// Defined for debugging in main )                                                                                                                                                                        
// [0]: output directory                                                                                                                                                                                  
// [1]: name for the histogram file                                                                                                                                                                       
// [2]: input data: can be a single .root or a .txt or .list of root files - should always be last argument  

int main (int argc, const char ** argv) {
    
  TH1::SetDefaultSumw2( );  // Histograms will calculate gaussian errors
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );

  // Read in command line arguments                                                                                                                                                                       
  // ------------------------------                                                                                                                                                                       
  // Defaults 
  std::string        outputDir         = "out/";                                        // directory where everything will be saved                                                                       
  std::string     outFileName        = "test.root"; 
  std::string         chainList            = "list.txt";            // input file: can be .root, .txt, .list
  
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
  InitReader(P6Reader, P6Chain, nEvents, sim_triggerString, truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, sim_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, sim_badTowers);
  InitReader(GEANTReader, GEANTChain, nEvents, sim_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, sim_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, sim_badTowers);

  TStarJetPicoEventHeader* p_header;  TStarJetPicoEvent* p_event;  TStarJetPicoEventHeader* g_header;  TStarJetPicoEvent* g_event;
  TStarJetVectorContainer<TStarJetVector> * p_container;         TStarJetVector* p_sv;
  TStarJetVectorContainer<TStarJetVector> * g_container;        TStarJetVector* g_sv;

  // Histograms 

  Collection<string, TH1D> hists1D; Collection<string, TH2D> hists2D; Collection<string, TH3D> hists3D;

  vector<string> flag_i = {"ch", "full"};
  vector<string> flag_j = {"lead", "sublead", "trig", "rec", "incl"};
  vector<string> flag_k = {"jet", "cons"};
  vector<string> flag_l = {"py", "ge"};
  for (int i = 0; i < flag_i.size(); ++ i) {
    for (int j = 0; j < flag_j.size(); ++ j) {
      for (int k = 0; k < flag_k.size(); ++ k) {
	for (int l = 0; l < flag_l.size(); ++ l) {
	  hists1D.add(("m_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]+"_"+flag_l[l]).c_str(),"",80,0,80); //mass
	  hists2D.add(("m_v_pt_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]+"_"+flag_l[l]).c_str(),"",20,0,20,80,0,80); //mass vs. pT
	  hists3D.add(("PtEtaPhi_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]+"_"+flag_l[l]).c_str(),"",80,0,80,5,-1,1,10,0,2*Pi);
	}
      }
    }
  }

  double jetPt_incl, jetEta_incl, jetPhi_incl, jetM_incl, jetE_incl, consPt_incl, consEta_incl, consPhi_incl, consM_incl, consE_incl, pythia_wt, geant_wt;
  int nCons_incl, cons_dummy;
  double jetPt_lead, jetEta_lead, jetPhi_lead, jetM_lead, jetE_lead, consPt_lead, consEta_lead, consPhi_lead, consM_lead, consE_lead;
  double jetPt_sublead, jetEta_sublead, jetPhi_sublead, jetM_sublead, jetE_sublead;
  int nCons_lead, nCons_sublead;

  TTree *p_inclTree = new TTree("py_inclTree","py_inclTree");
  TTree *g_inclTree = new TTree("ge_inclTree","ge_inclTree");
  TTree *p_leadTree = new TTree("py_leadTree","py_leadTree");
  TTree *g_leadTree = new TTree("ge_leadTree","ge_leadTree");
  TTree *p_subleadTree = new TTree("py_subleadTree","py_subleadTree");
  TTree *g_subleadTree = new TTree("ge_subleadTree","ge_subleadTree");
  TTree *p_cons_inclTree = new TTree("py_cons_inclTree","py_cons_inclTree");
  TTree *g_cons_inclTree = new TTree("ge_cons_inclTree","ge_cons_inclTree");
  TTree *p_cons_leadTree = new TTree("py_cons_leadTree","py_cons_leadTree");
  TTree *g_cons_leadTree = new TTree("ge_cons_leadTree","ge_cons_leadTree");
  
  p_leadTree->Branch("Pt", &jetPt_lead); p_leadTree->Branch("Eta", &jetEta_lead); p_leadTree->Branch("Phi", &jetPhi_lead);
  p_leadTree->Branch("M", &jetM_lead); p_leadTree->Branch("E", &jetE_lead); p_leadTree->Branch("nCons", &nCons_lead);
  p_leadTree->Branch("weight", &pythia_wt);

  g_leadTree->Branch("Pt", &jetPt_lead); g_leadTree->Branch("Eta", &jetEta_lead); g_leadTree->Branch("Phi", &jetPhi_lead);
  g_leadTree->Branch("M", &jetM_lead); g_leadTree->Branch("E", &jetE_lead); g_leadTree->Branch("nCons", &nCons_lead);
  g_leadTree->Branch("weight", &geant_wt);

  p_subleadTree->Branch("Pt", &jetPt_sublead); p_subleadTree->Branch("Eta", &jetEta_sublead); p_subleadTree->Branch("Phi", &jetPhi_sublead);
  p_subleadTree->Branch("M", &jetM_sublead); p_subleadTree->Branch("E", &jetE_sublead); p_subleadTree->Branch("nCons", &nCons_sublead);
  p_subleadTree->Branch("weight", &pythia_wt);

  g_subleadTree->Branch("Pt", &jetPt_sublead); g_subleadTree->Branch("Eta", &jetEta_sublead); g_subleadTree->Branch("Phi", &jetPhi_sublead);
  g_subleadTree->Branch("M", &jetM_sublead); g_subleadTree->Branch("E", &jetE_sublead); g_subleadTree->Branch("nCons", &nCons_sublead);
  g_subleadTree->Branch("weight", &geant_wt);

  p_cons_leadTree->Branch("Pt", &consPt_lead); p_cons_leadTree->Branch("Eta", &consEta_lead); p_cons_leadTree->Branch("Phi", &consPhi_lead);
  p_cons_leadTree->Branch("M", &consM_lead); p_cons_leadTree->Branch("E", &consE_lead);
  p_cons_leadTree->Branch("weight", &pythia_wt);

  g_cons_leadTree->Branch("Pt", &consPt_lead); g_cons_leadTree->Branch("Eta", &consEta_lead); g_cons_leadTree->Branch("Phi", &consPhi_lead);
  g_cons_leadTree->Branch("M", &consM_lead); g_cons_leadTree->Branch("E", &consE_lead);
  g_cons_leadTree->Branch("weight", &geant_wt);
 
  p_inclTree->Branch("Pt", &jetPt_incl); p_inclTree->Branch("Eta", &jetEta_incl); p_inclTree->Branch("Phi", &jetPhi_incl);
  p_inclTree->Branch("M", &jetM_incl); p_inclTree->Branch("E", &jetE_incl); p_inclTree->Branch("nCons", &nCons_incl);
  p_inclTree->Branch("weight", &pythia_wt);

  g_inclTree->Branch("Pt", &jetPt_incl); g_inclTree->Branch("Eta", &jetEta_incl); g_inclTree->Branch("Phi", &jetPhi_incl);
  g_inclTree->Branch("M", &jetM_incl); g_inclTree->Branch("E", &jetE_incl); g_inclTree->Branch("nCons", &nCons_incl);
  g_inclTree->Branch("weight", &geant_wt);

  p_cons_inclTree->Branch("Pt", &consPt_incl); p_cons_inclTree->Branch("Eta", &consEta_incl); p_cons_inclTree->Branch("Phi", &consPhi_incl);
  p_cons_inclTree->Branch("M", &consM_incl); p_cons_inclTree->Branch("E", &consE_incl);
  p_cons_inclTree->Branch("weight", &pythia_wt);

  g_cons_inclTree->Branch("Pt", &consPt_incl); g_cons_inclTree->Branch("Eta", &consEta_incl); g_cons_inclTree->Branch("Phi", &consPhi_incl);
  g_cons_inclTree->Branch("M", &consM_incl); g_cons_inclTree->Branch("E", &consE_incl);
  g_cons_inclTree->Branch("weight", &geant_wt);
  
  //SELECTORS
  // Constituent selectors                                                                                                                                                     
  // ---------------------                                                                                                                                                        
  Selector select_track_rap = fastjet::SelectorAbsRapMax(max_track_rap);
  Selector select_lopt      = fastjet::SelectorPtMin( partMinPt );
  Selector spart = select_track_rap * select_lopt;
  
  // Jet candidate selectors                                                                                                                                                   
  // -----------------------                                                                                                                                                      
  Selector select_jet_rap     = fastjet::SelectorAbsRapMax(max_rap);
  Selector select_jet_pt_min  = fastjet::SelectorPtMin( jet_ptmin );
  Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
  Selector sjet = select_jet_rap && select_jet_pt_min && select_jet_pt_max;
  
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  double p_wt, g_wt;
  TString geantFilename, pythiaFilename;
  
  // Particle containers & counters
  vector<PseudoJet> p_Particles, g_Particles, ch_p_Particles, ch_g_Particles, p_Jets, g_Jets, ch_p_Jets, ch_g_Jets;//, p_Cons, g_Cons;
  int nEvents = 0;   int p_NJets = 0;  int g_NJets = 0;  int p_EventID;   int g_EventID;
  //1=inclusive, 2=lead
  int counter_debug1 = 0, counter_debug2 = 0;
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( GEANTReader.NextEvent() ) {      //    GEANTReader    P6Reader
    
    g_EventID = GEANTReader.GetNOfCurrentEvent();

    if ( P6Reader.ReadEvent( g_EventID ) != 1 ) continue;   //  ENSURES BOTH DATASETS HAVE AN EVENT

    p_Particles.clear(); g_Particles.clear();
    ch_p_Particles.clear(); ch_g_Particles.clear();
    p_Jets.clear(); g_Jets.clear(); 
    ch_p_Jets.clear(); ch_g_Jets.clear();
    //p_Cons.clear(); g_Cons.clear();   //  clear all containers

    nEvents++;  P6Reader.PrintStatus(10);  GEANTReader.PrintStatus(10);     // Print out reader status every 10 seconds

    p_event = P6Reader.GetEvent();       p_header = p_event->GetHeader();           // Get the PYTHIA header and event
    g_event = GEANTReader.GetEvent();    g_header = g_event->GetHeader();           // Get GEANT event header and event
    
    p_EventID = P6Reader.GetNOfCurrentEvent();
    if ( p_EventID != g_EventID ) { cout << endl << "ERROR: READING DIFFERENT EVENTS " <<endl; }
    
    p_container = P6Reader.GetOutputContainer();      // Pythia container
    g_container = GEANTReader.GetOutputContainer();      // GEANT container

    pythiaFilename =  P6Reader.GetInputChain()->GetCurrentFile()->GetName();	
    geantFilename =  GEANTReader.GetInputChain()->GetCurrentFile()->GetName();	
    
    p_wt = LookupXsec( pythiaFilename );
    g_wt = LookupXsec( geantFilename );

    //  GATHER PARTICLES
    GatherParticles ( p_container, p_sv, p_Particles, 1);    //  Pythia particles
    GatherParticles ( p_container, p_sv, ch_p_Particles, 0); //0 signifies charged-only, 1 signifies ch+ne
    GatherParticles ( g_container, g_sv, g_Particles, 1);    //  GEANT particles
    GatherParticles ( g_container, g_sv, ch_g_Particles, 0);
    
    vector<PseudoJet> p_cut_Particles = spart(p_Particles); vector<PseudoJet> g_cut_Particles = spart(g_Particles); //applying constituent cuts
    vector<PseudoJet> ch_p_cut_Particles = spart(ch_p_Particles); vector<PseudoJet> ch_g_cut_Particles = spart(ch_g_Particles);
    
    ClusterSequence p_Cluster(p_cut_Particles, jet_def); ClusterSequence g_Cluster(g_cut_Particles, jet_def);           //  CLUSTER BOTH
    ClusterSequence ch_p_Cluster(ch_p_cut_Particles, jet_def); ClusterSequence ch_g_Cluster(ch_g_cut_Particles, jet_def);
    p_Jets = sorted_by_pt(sjet(p_Cluster.inclusive_jets())); g_Jets = sorted_by_pt(sjet(g_Cluster.inclusive_jets()));    // EXTRACT JETS
    ch_p_Jets = sorted_by_pt(sjet(p_Cluster.inclusive_jets())); ch_g_Jets = sorted_by_pt(sjet(ch_g_Cluster.inclusive_jets()));

    if (p_Jets.size() != 0) {
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TREES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pythia~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      //leading
      vector<PseudoJet> pLead; pLead.push_back(p_Jets[0]);
      FillTrees(pLead, p_leadTree, jetPt_lead, jetEta_lead, jetPhi_lead, jetM_lead, jetE_lead, nCons_lead, pythia_wt, p_wt);
      //subleading
      if (p_Jets.size() > 1) {
	vector<PseudoJet> pSublead; pSublead.push_back(p_Jets[1]);
	FillTrees(pSublead, p_subleadTree, jetPt_sublead, jetEta_sublead, jetPhi_sublead, jetM_sublead, jetE_sublead, nCons_sublead, pythia_wt, p_wt);
      }
      //constituents
      FillTrees(pLead[0].constituents(), p_cons_leadTree, consPt_lead, consEta_lead, consPhi_lead, consM_lead, consE_lead, cons_dummy, pythia_wt, p_wt);
      //inclusive
      FillTrees(p_Jets, p_inclTree, jetPt_incl, jetEta_incl, jetPhi_incl, jetM_incl, jetE_incl, nCons_incl, pythia_wt, p_wt);
      //constituents
      for(int j = 0; j < p_Jets.size(); ++ j) {
	FillTrees(p_Jets[j].constituents(), p_cons_inclTree, consPt_incl, consEta_incl, consPhi_incl, consM_incl, consE_incl, cons_dummy, pythia_wt, p_wt); 
      }
    }
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    if (g_Jets.size() != 0) { 
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~geant~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      //leading
      vector<PseudoJet> gLead; gLead.push_back(g_Jets[0]);
      FillTrees(gLead, g_leadTree, jetPt_lead, jetEta_lead, jetPhi_lead, jetM_lead, jetE_lead, nCons_lead, geant_wt, g_wt);
      //subleading
       if (g_Jets.size() > 1) {
	vector<PseudoJet> gSublead; gSublead.push_back(g_Jets[1]);
	FillTrees(gSublead, g_subleadTree, jetPt_sublead, jetEta_sublead, jetPhi_sublead, jetM_sublead, jetE_sublead, nCons_sublead, geant_wt, g_wt);
      }
      //constituents
      FillTrees(gLead[0].constituents(), g_cons_leadTree, consPt_lead, consEta_lead, consPhi_lead, consM_lead, consE_lead, cons_dummy, geant_wt, g_wt);  
      //inclusive
      FillTrees(g_Jets, g_inclTree, jetPt_incl, jetEta_incl, jetPhi_incl, jetM_incl, jetE_incl, nCons_incl, geant_wt, g_wt);
      //constituents
      for(int j = 0; j < g_Jets.size(); ++ j) {
	FillTrees(g_Jets[j].constituents(), g_cons_inclTree, consPt_incl, consEta_incl, consPhi_incl, consM_incl, consE_incl, cons_dummy, geant_wt, g_wt);
      }
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HISTOS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    if(p_Jets.size() != 0) {
      FillHists(hists1D, hists2D, hists3D, "full", "_py", p_Jets, p_wt);
    }
    if(g_Jets.size() != 0) {
      FillHists(hists1D, hists2D, hists3D, "full", "_ge", g_Jets, g_wt);
    }
    if(ch_p_Jets.size() != 0) {
      FillHists(hists1D, hists2D, hists3D, "ch", "_py", ch_p_Jets, p_wt);
    }
    if(ch_g_Jets.size() != 0) {
      FillHists(hists1D, hists2D, hists3D, "ch", "_ge", ch_g_Jets, g_wt);
    }
    p_NJets += p_Jets.size(); g_NJets += g_Jets.size();               //  Save jet info and add jets to total
  }

  //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");

  std::cout << std::endl << std::endl << "Of " << nEvents << " events: "<< std::endl;
  std::cout << p_NJets << " jets have been found for the Pythia6 data" << std::endl;
  std::cout << g_NJets << " jets have been found for the Pythia6 + GEANT data" << std::endl << std::endl;
  std::cout <<std::endl << "Writing to:  " << fout->GetName() << std::endl << std::endl;

  p_leadTree->Write("py_leadTree");
  g_leadTree->Write("ge_leadTree");
  p_subleadTree->Write("py_subleadTree");
  g_subleadTree->Write("ge_subleadTree");
  p_cons_leadTree->Write("py_cons_leadTree");
  g_cons_leadTree->Write("ge_cons_leadTree");
  p_inclTree->Write("py_inclTree");
  g_inclTree->Write("ge_inclTree");
  p_cons_inclTree->Write("py_cons_inclTree");
  g_cons_inclTree->Write("ge_cons_inclTree");

  
  //hist1->Write(); hist2->Write(); etc, goes here
  for (int i = 0; i < flag_i.size(); ++ i) {
    for (int j = 0; j < flag_j.size(); ++ j) {
      for (int k = 0; k < flag_k.size(); ++ k) {
	for (int l = 0; l < flag_l.size(); ++ l) {
	  hists1D.write(("m_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]+"_"+flag_l[l]).c_str());
	  hists2D.write(("m_v_pt_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]+"_"+flag_l[l]).c_str());
	  hists3D.write(("PtEtaPhi_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]+"_"+flag_l[l]).c_str());
	}
      }
    }
  }

  fout->Close();

  return 0;
}
