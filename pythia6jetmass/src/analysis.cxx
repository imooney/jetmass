//  Veronica Verkest        May 13, 2018
//  Compare:   p6  VS  p6+efficiency  VS  p6+GEANT
//  Functions in src/functions.cxx
//  Parameters in src/parameters.cxx
//  Adapted by Isaac Mooney June, 2018 for jet mass analysis

#include "parameters.hh"
#include "functions.hh"

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


  /*P6Chain->Add("AddedGeantPythia/picoDst_11_15_0.root");*/ InitReaderPythia(P6Reader, P6Chain, numEvents);
  /*GEANTChain->Add("AddedGeantPythia/picoDst_11_15_0.root");*/ InitReaderGeant(GEANTReader, GEANTChain, numEvents);
  
  //P6Chain->Add( "AddedGeantPythia/picoDst*" );          InitReaderPythia( P6Reader, P6Chain, numEvents );
  //GEANTChain->Add( "AddedGeantPythia/picoDst*" );   InitReaderGeant( GEANTReader, GEANTChain, numEvents );
  
  TStarJetPicoEventHeader* p_header;  TStarJetPicoEvent* p_event;  TStarJetPicoEventHeader* g_header;  TStarJetPicoEvent* g_event;
  TStarJetVectorContainer<TStarJetVector> * p_container;         TStarJetVector* p_sv;
  TStarJetVectorContainer<TStarJetVector> * g_container;        TStarJetVector* g_sv;

  // Histograms 

  TH3 * p_PtEtaPhi_incl = new TH3D("p_PtEtaPhi_incl","",80,0,80,5,-1,1,10,0,2*Pi);
  TH3 * p_PtEtaPhi_lead = new TH3D("p_PtEtaPhi_lead","",80,0,80,5,-1,1,10,0,2*Pi);

  TH3 * p_cons_PtEtaPhi_incl = new TH3D("p_cons_PtEtaPhi_incl","",80,0,80,5,-1,1,10,0,2*Pi);
  TH3 * p_cons_PtEtaPhi_lead = new TH3D("p_cons_PtEtaPhi_lead","",80,0,80,5,-1,1,10,0,2*Pi);

  TH1 * p_m_incl = new TH1D("p_m_incl","",80,0,80);
  TH1 * p_m_lead = new TH1D("p_m_lead","",80,0,80);

  TH3 * g_PtEtaPhi_incl = new TH3D("g_PtEtaPhi_incl","",80,0,80,5,-1,1,10,0,2*Pi);
  TH3 * g_PtEtaPhi_lead = new TH3D("g_PtEtaPhi_lead","",80,0,80,5,-1,1,10,0,2*Pi);
  
  TH3 * g_cons_PtEtaPhi_incl = new TH3D("g_cons_PtEtaPhi_incl","",80,0,80,5,-1,1,10,0,2*Pi);
  TH3 * g_cons_PtEtaPhi_lead = new TH3D("g_cons_PtEtaPhi_lead","",80,0,80,5,-1,1,10,0,2*Pi);

  TH1 * g_m_incl = new TH1D("g_m_incl","",80,0,80);
  TH1 * g_m_lead = new TH1D("g_m_lead","",80,0,80);

  TH2 * p_m_v_pt_lead = new TH2D("p_m_v_pt_lead","",20,0,20,80,0,80);
  TH2 * p_m_v_pt_incl = new TH2D("p_m_v_pt_incl","",20,0,20,80,0,80);

  TH2 * g_m_v_pt_lead = new TH2D("g_m_v_pt_lead","",20,0,20,80,0,80);
  TH2 * g_m_v_pt_incl = new TH2D("g_m_v_pt_incl","",20,0,20,80,0,80);

  double jetPt_incl, jetEta_incl, jetPhi_incl, jetM_incl, jetE_incl, consPt_incl, consEta_incl, consPhi_incl, consM_incl, consE_incl, pythia_wt, geant_wt;
  int nCons_incl, cons_dummy;
  double jetPt_lead, jetEta_lead, jetPhi_lead, jetM_lead, jetE_lead, consPt_lead, consEta_lead, consPhi_lead, consM_lead, consE_lead;
  int nCons_lead;

  TTree *p_inclTree = new TTree("py_inclTree","py_inclTree");
  TTree *g_inclTree = new TTree("ge_inclTree","ge_inclTree");
   TTree *p_leadTree = new TTree("py_leadTree","py_leadTree");
  TTree *g_leadTree = new TTree("ge_leadTree","ge_leadTree");
   TTree *p_cons_inclTree = new TTree("py_cons_inclTree","py_cons_inclTree");
  TTree *g_cons_inclTree = new TTree("ge_cons_inclTree","ge_cons_inclTree");
   TTree *p_cons_leadTree = new TTree("py_cons_leadTree","py_cons_leadTree");
  TTree *g_cons_leadTree = new TTree("ge_cons_leadTree","ge_cons_leadTree");
  
  TBranch *pythiaPt_lead;  TBranch *pythiaEta_lead;  TBranch *pythiaPhi_lead; TBranch * pythiaM_lead; TBranch *pythiaE_lead; TBranch *pythiaNCons_lead; TBranch *pythiaWeight;
  pythiaPt_lead = p_leadTree->Branch("jetPt", &jetPt_lead);  pythiaEta_lead = p_leadTree->Branch("jetEta", &jetEta_lead);  pythiaPhi_lead = p_leadTree->Branch("jetPhi", &jetPhi_lead);
  pythiaM_lead = p_leadTree->Branch("jetM", &jetM_lead); pythiaE_lead = p_leadTree->Branch("jetE", &jetE_lead); pythiaNCons_lead = p_leadTree->Branch("nCons", &nCons_lead);
  pythiaWeight = p_leadTree->Branch("pythia_wt", &pythia_wt);

  TBranch *geantPt_lead;  TBranch *geantEta_lead;  TBranch *geantPhi_lead; TBranch *geantM_lead; TBranch *geantE_lead; TBranch *geantNCons_lead; TBranch *geantWeight;
  geantPt_lead = g_leadTree->Branch("jetPt", &jetPt_lead);  geantEta_lead = g_leadTree->Branch("jetEta", &jetEta_lead);  geantPhi_lead = g_leadTree->Branch("jetPhi", &jetPhi_lead);
  geantM_lead = g_leadTree->Branch("jetM", &jetM_lead); geantE_lead = g_leadTree->Branch("jetE", &jetE_lead); geantNCons_lead = g_leadTree->Branch("nCons", &nCons_lead);
  geantWeight = g_leadTree->Branch("geant_wt", &geant_wt);

  TBranch *c_pythiaPt_lead;  TBranch *c_pythiaEta_lead;  TBranch *c_pythiaPhi_lead; TBranch * c_pythiaM_lead; TBranch *c_pythiaE_lead;
  c_pythiaPt_lead = p_cons_leadTree->Branch("consPt", &consPt_lead);  c_pythiaEta_lead = p_cons_leadTree->Branch("consEta", &consEta_lead);  c_pythiaPhi_lead = p_cons_leadTree->Branch("consPhi", &consPhi_lead);
  c_pythiaM_lead = p_cons_leadTree->Branch("consM", &consM_lead); c_pythiaE_lead = p_cons_leadTree->Branch("consE", &consE_lead);
  pythiaWeight = p_cons_leadTree->Branch("pythia_wt", &pythia_wt);

  TBranch *c_geantPt_lead;  TBranch *c_geantEta_lead;  TBranch *c_geantPhi_lead; TBranch *c_geantM_lead; TBranch *c_geantE_lead;
  c_geantPt_lead = g_cons_leadTree->Branch("consPt", &consPt_lead);  c_geantEta_lead = g_cons_leadTree->Branch("consEta", &consEta_lead);  c_geantPhi_lead = g_cons_leadTree->Branch("consPhi", &consPhi_lead);
  c_geantM_lead = g_cons_leadTree->Branch("consM", &consM_lead); c_geantE_lead = g_cons_leadTree->Branch("consE", &consE_lead);
  geantWeight = g_cons_leadTree->Branch("geant_wt", &geant_wt);
 
   TBranch *pythiaPt_incl;  TBranch *pythiaEta_incl;  TBranch *pythiaPhi_incl; TBranch * pythiaM_incl; TBranch *pythiaE_incl; TBranch *pythiaNCons_incl;
  pythiaPt_incl = p_inclTree->Branch("jetPt", &jetPt_incl);  pythiaEta_incl = p_inclTree->Branch("jetEta", &jetEta_incl);  pythiaPhi_incl = p_inclTree->Branch("jetPhi", &jetPhi_incl);
  pythiaM_incl = p_inclTree->Branch("jetM", &jetM_incl); pythiaE_incl = p_inclTree->Branch("jetE", &jetE_incl); pythiaNCons_lead = p_inclTree->Branch("nCons", &nCons_incl);
  pythiaWeight = p_inclTree->Branch("pythia_wt", &pythia_wt);

  TBranch *geantPt_incl;  TBranch *geantEta_incl;  TBranch *geantPhi_incl; TBranch *geantM_incl; TBranch *geantE_incl; TBranch *geantNCons_incl;
  geantPt_incl = g_inclTree->Branch("jetPt", &jetPt_incl);  geantEta_incl = g_inclTree->Branch("jetEta", &jetEta_incl);  geantPhi_incl = g_inclTree->Branch("jetPhi", &jetPhi_incl);
  geantM_incl = g_inclTree->Branch("jetM", &jetM_incl); geantE_incl = g_inclTree->Branch("jetE", &jetE_incl); geantNCons_incl = g_inclTree->Branch("nCons", &nCons_incl);
  geantWeight = g_inclTree->Branch("geant_wt", &geant_wt);

  TBranch *c_pythiaPt_incl;  TBranch *c_pythiaEta_incl;  TBranch *c_pythiaPhi_incl; TBranch * c_pythiaM_incl; TBranch *c_pythiaE_incl;
  c_pythiaPt_incl = p_cons_inclTree->Branch("consPt", &consPt_incl);  c_pythiaEta_incl = p_cons_inclTree->Branch("consEta", &consEta_incl);  c_pythiaPhi_incl = p_cons_inclTree->Branch("consPhi", &consPhi_incl);
  c_pythiaM_incl = p_cons_inclTree->Branch("consM", &consM_incl); c_pythiaE_incl = p_cons_inclTree->Branch("consE", &consE_incl);
  pythiaWeight = p_cons_inclTree->Branch("pythia_wt", &pythia_wt);

  TBranch *c_geantPt_incl;  TBranch *c_geantEta_incl;  TBranch *c_geantPhi_incl; TBranch *c_geantM_incl; TBranch *c_geantE_incl;
  c_geantPt_incl = g_cons_inclTree->Branch("consPt", &consPt_incl);  c_geantEta_incl = g_cons_inclTree->Branch("consEta", &consEta_incl);  c_geantPhi_incl = g_cons_inclTree->Branch("consPhi", &consPhi_incl);
  c_geantM_incl = g_cons_inclTree->Branch("consM", &consM_incl); c_geantE_incl = g_cons_inclTree->Branch("consE", &consE_incl);
  geantWeight = g_cons_inclTree->Branch("geant_wt", &geant_wt);
  

  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  double p_wt, g_wt;
  TString geantFilename, pythiaFilename;
  
  // Particle containers & counters
  vector<PseudoJet> p_Particles, g_Particles, p_Jets, g_Jets, p_Cons, g_Cons;
  int nEvents = 0;   int p_NJets = 0;  int g_NJets = 0;  int p_EventID;   int g_EventID;
  //1=inclusive, 2=lead
  int counter_debug1 = 0, counter_debug2 = 0;
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( GEANTReader.NextEvent() ) {      //    GEANTReader    P6Reader
    
    g_EventID = GEANTReader.GetNOfCurrentEvent();

    if ( P6Reader.ReadEvent( g_EventID ) != 1 ) continue;   //  ENSURES BOTH DATASETS HAVE AN EVENT

    p_Particles.clear(), g_Particles.clear(),
    p_Jets.clear(), g_Jets.clear(), 
    p_Cons.clear(),  g_Cons.clear();   //  clear all containers

    nEvents++;  P6Reader.PrintStatus(10);  GEANTReader.PrintStatus(10);     // Print out reader status every 10 seconds

    p_event = P6Reader.GetEvent();       p_header = p_event->GetHeader();           // Get the PYTHIA header and event
    g_event = GEANTReader.GetEvent();    g_header = g_event->GetHeader();           // Get GEANT event header and event
    
    p_EventID = P6Reader.GetNOfCurrentEvent();
    if ( p_EventID != g_EventID ) { cout << endl << "ERROR: READING DIFFERENT EVENTS " <<endl; }
    
    p_container = P6Reader.GetOutputContainer();      // Pythia container
    g_container = GEANTReader.GetOutputContainer();      // GEANT container

    //  APPLY VERTEX-Z CUT
    if ( Vz_candidate( g_header, absMaxVz ) == false ) { continue; }
    // if ( Vz_candidate( p_header, absMaxVz ) == false || Vz_candidate( g_header, absMaxVz ) == false ) { continue; }

    pythiaFilename =  P6Reader.GetInputChain()->GetCurrentFile()->GetName();	
    geantFilename =  GEANTReader.GetInputChain()->GetCurrentFile()->GetName();	
    
    p_wt = LookupXsec( pythiaFilename );
    g_wt = LookupXsec( geantFilename );

    //  GATHER PARTICLES
    GatherParticles ( p_container, etaCut, partMinPt, p_Particles);    //  Pythia particles
    GatherParticles ( g_container, etaCut, partMinPt, g_Particles);    //  GEANT particles

    //  CREATE JET SELECTOR
    Selector etaSelector = SelectorAbsEtaMax( 1.0-R );    Selector ptSelector = SelectorPtMin(jetMinPt);    Selector etaPtSelector = etaSelector && ptSelector;
    
    ClusterSequence p_Cluster(p_Particles, jet_def); ClusterSequence g_Cluster(g_Particles, jet_def);           //  CLUSTER BOTH
    p_Jets = sorted_by_pt(etaPtSelector(p_Cluster.inclusive_jets())); g_Jets = sorted_by_pt(etaPtSelector(g_Cluster.inclusive_jets()));    // EXTRACT JETS

    if (p_Jets.size() != 0) {
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TREES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pythia~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      //leading
      vector<PseudoJet> pLead; pLead.push_back(p_Jets[0]);
      FillTrees(pLead, p_leadTree, jetPt_lead, jetEta_lead, jetPhi_lead, jetM_lead, jetE_lead, nCons_lead, pythia_wt, p_wt);
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
    //pythia
    for (int i = 0; i < p_Jets.size(); ++ i) {
      if (i == 0) { //lead
	p_m_lead->Fill(p_Jets[0].m(), p_wt);
	p_m_v_pt_lead->Fill(p_Jets[0].m(), p_Jets[0].pt(), p_wt);
	p_PtEtaPhi_lead->Fill(p_Jets[0].pt(), p_Jets[0].eta(), p_Jets[0].phi(), p_wt); //jets
	vector<PseudoJet> cons_lead = p_Jets[0].constituents();
	for (int j = 0; j < cons_lead.size(); ++ j) {
	  p_cons_PtEtaPhi_lead->Fill(cons_lead[j].pt(), cons_lead[j].eta(), cons_lead[j].phi(), p_wt); //constituents
	}
      }
      //inclusive
      p_m_incl->Fill(p_Jets[i].m(), p_wt);
      p_m_v_pt_incl->Fill(p_Jets[i].m(), p_Jets[i].pt(), p_wt);
      p_PtEtaPhi_incl->Fill(p_Jets[i].pt(),p_Jets[i].eta(),p_Jets[i].phi(),p_wt); //jets
      for (int j = 0; j < p_Jets[i].constituents().size(); ++ j) {
	p_cons_PtEtaPhi_incl->Fill(p_Jets[i].constituents()[j].pt(), p_Jets[i].constituents()[j].eta(), p_Jets[i].constituents()[j].phi(), p_wt);
      }
    }

    //geant
    for (int i = 0; i < g_Jets.size(); ++ i) {
      if (i == 0) { //lead
	g_m_lead->Fill(g_Jets[0].m(), g_wt);
	g_m_v_pt_lead->Fill(g_Jets[0].m(), g_Jets[0].pt(), g_wt);
	g_PtEtaPhi_lead->Fill(g_Jets[0].pt(), g_Jets[0].eta(), g_Jets[0].phi(), g_wt); //jets
	vector<PseudoJet> cons_lead = g_Jets[0].constituents();
	for (int j = 0; j < cons_lead.size(); ++ j) {
	  g_cons_PtEtaPhi_lead->Fill(cons_lead[j].pt(), cons_lead[j].eta(), cons_lead[j].phi(), g_wt); //constituents
	}
      }
      //inclusive
      g_m_incl->Fill(g_Jets[i].m(), g_wt);
      g_m_v_pt_incl->Fill(g_Jets[i].m(), g_Jets[i].pt(), g_wt);
      g_PtEtaPhi_incl->Fill(g_Jets[i].pt(),g_Jets[i].eta(),g_Jets[i].phi(),g_wt); //jets
      for (int j = 0; j < g_Jets[i].constituents().size(); ++ j) {
	g_cons_PtEtaPhi_incl->Fill(g_Jets[i].constituents()[j].pt(), g_Jets[i].constituents()[j].eta(), g_Jets[i].constituents()[j].phi(), g_wt);
      }
    }
    p_NJets += p_Jets.size(); g_NJets += g_Jets.size();               //  Save jet info and add jets to total
  }

  //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");

  std::cout << std::endl << std::endl << "Of " << nEvents << " events: "<< std::endl;
  std::cout << p_NJets << " jets have been found for the Pythia6 data" << std::endl;
  std::cout << g_NJets << " jets have been found for the Pythia6 + GEANT data" << std::endl << std::endl;
  std::cout <<std::endl << "Writing to:  " << fout->GetName() << std::endl << std::endl;

  std::cout << "TESTING: number of leading jets = " << counter_debug2 << " and number of inclusive jets = " << counter_debug1 << std::endl;
  p_leadTree->Write("py_leadTree");
  g_leadTree->Write("ge_leadTree");
  p_cons_leadTree->Write("py_cons_leadTree");
  g_cons_leadTree->Write("ge_cons_leadTree");
  p_inclTree->Write("py_inclTree");
  g_inclTree->Write("ge_inclTree");
  p_cons_inclTree->Write("py_cons_inclTree");
  g_cons_inclTree->Write("ge_cons_inclTree");

  
  //once there are hists: hist1->Write(); hist2->Write(); etc, goes here
  p_PtEtaPhi_lead->Write(); p_cons_PtEtaPhi_lead->Write(); p_PtEtaPhi_incl->Write(); p_cons_PtEtaPhi_incl->Write();
  g_PtEtaPhi_lead->Write(); g_cons_PtEtaPhi_lead->Write(); g_PtEtaPhi_incl->Write(); g_cons_PtEtaPhi_incl->Write(); 
  p_m_incl->Write(); p_m_lead->Write(); g_m_incl->Write(); g_m_lead->Write();
  p_m_v_pt_incl->Write(); p_m_v_pt_lead->Write(); g_m_v_pt_incl->Write(); g_m_v_pt_lead->Write();
  

  fout->Close();

  return 0;
}
