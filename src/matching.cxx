//  Veronica Verkest        May 13, 2018
//  Compare:   p6  VS  p6+efficiency  VS  p6+GEANT
//  Functions in src/functions.cxx
//  Parameters in src/parameters.cxx
//  Adapted by Isaac Mooney June, 2018 for jet mass analysis

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
  InitReader(P6Reader, P6Chain, nEvents, truth_triggerString, truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, 0.9999, false, sim_badTowers, sim_bad_run_list);
  InitReader(GEANTReader, GEANTChain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, 0.9999, false, det_badTowers, dat_bad_run_list/*sim_badTowers, sim_bad_run_list*/);
  
  TStarJetPicoEventHeader* p_header;  TStarJetPicoEvent* p_event;  TStarJetPicoEventHeader* g_header;  TStarJetPicoEvent* g_event;
  TStarJetVectorContainer<TStarJetVector> * p_container;         TStarJetVector* p_sv;
  TStarJetVectorContainer<TStarJetVector> * g_container;        TStarJetVector* g_sv;
  
  //detector response fits from file
  TFile * funcs_in = new TFile("~/jetmass/macros/funcs/funcs.root","READ");
  
  TF1 * scale_lo_0 = (TF1*) funcs_in->Get("scale_lo_0");
  TF1 * scale_lo_1 = (TF1*) funcs_in->Get("scale_lo_1");
  TF1 * scale_lo_2 = (TF1*) funcs_in->Get("scale_lo_2");
  TF1 * scale_lo_3 = (TF1*) funcs_in->Get("scale_lo_3");
  TF1 * scale_lo_4 = (TF1*) funcs_in->Get("scale_lo_4");
  TF1 * scale_mid_0 = (TF1*) funcs_in->Get("scale_mid_0");
  TF1 * scale_mid_1 = (TF1*) funcs_in->Get("scale_mid_1");
  TF1 * scale_mid_2 = (TF1*) funcs_in->Get("scale_mid_2");
  TF1 * scale_mid_3 = (TF1*) funcs_in->Get("scale_mid_3");
  TF1 * scale_mid_4 = (TF1*) funcs_in->Get("scale_mid_4");
  TF1 * scale_hi_0 = (TF1*) funcs_in->Get("scale_hi_0");
  TF1 * scale_hi_1 = (TF1*) funcs_in->Get("scale_hi_1");
  TF1 * scale_hi_2 = (TF1*) funcs_in->Get("scale_hi_2");
  TF1 * scale_hi_3 = (TF1*) funcs_in->Get("scale_hi_3");
  TF1 * scale_hi_4 = (TF1*) funcs_in->Get("scale_hi_4");
  
  TF1 * groom_scale_lo_0 = (TF1*) funcs_in->Get("groom_scale_lo_0");
  TF1 * groom_scale_lo_1 = (TF1*) funcs_in->Get("groom_scale_lo_1");
  TF1 * groom_scale_lo_2 = (TF1*) funcs_in->Get("groom_scale_lo_2");
  TF1 * groom_scale_lo_3 = (TF1*) funcs_in->Get("groom_scale_lo_3");
  TF1 * groom_scale_lo_4 = (TF1*) funcs_in->Get("groom_scale_lo_4");
  TF1 * groom_scale_mid_0 = (TF1*) funcs_in->Get("groom_scale_mid_0");
  TF1 * groom_scale_mid_1 = (TF1*) funcs_in->Get("groom_scale_mid_1");
  TF1 * groom_scale_mid_2 = (TF1*) funcs_in->Get("groom_scale_mid_2");
  TF1 * groom_scale_mid_3 = (TF1*) funcs_in->Get("groom_scale_mid_3");
  TF1 * groom_scale_mid_4 = (TF1*) funcs_in->Get("groom_scale_mid_4");
  TF1 * groom_scale_hi_0 = (TF1*) funcs_in->Get("groom_scale_hi_0");
  TF1 * groom_scale_hi_1 = (TF1*) funcs_in->Get("groom_scale_hi_1");
  TF1 * groom_scale_hi_2 = (TF1*) funcs_in->Get("groom_scale_hi_2");
  TF1 * groom_scale_hi_3 = (TF1*) funcs_in->Get("groom_scale_hi_3");
  TF1 * groom_scale_hi_4 = (TF1*) funcs_in->Get("groom_scale_hi_4");
  
  
  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
  fout->cd();  
  //  vector<double> deltaPt; vector<double> deltaM; vector<double> deltaZg; vector<double> deltaRg;
  // vector<double> ratioPt; vector<double> ratioM; vector<double> ratioZg; vector<double> ratioRg;
  vector<double> pyPt; vector<double> pyM; vector<double> pyZg; vector<double> pyRg;
  vector<double> gePt; vector<double> geM; vector<double> geZg; vector<double> geRg;
  vector<double> pyPtg; vector<double> pyMg; vector<double> gePtg; vector<double> geMg;
  vector<double> pyEta; vector<double> geEta;
  vector<double> gePt_corr; vector<double> geM_corr; vector<double> geZg_corr; vector<double> geRg_corr;
  vector<double> gePtg_corr; vector<double> geMg_corr; vector<double> geEta_corr;
  double mc_weight; int p_EventID;

  //vector of vector of the doubles which will fill the branches to make it easier to pass all of them to the function
  vector<vector<double> > py_arr = {pyPt, pyM, pyZg, pyRg, pyPtg, pyMg, pyEta};
  vector<vector<double> > ge_arr = {gePt, geM, geZg, geRg, gePtg, geMg, geEta};
  vector<vector<double> > ge_corr = {gePt_corr, geM_corr, geZg_corr, geRg_corr, gePtg_corr, geMg_corr, geEta_corr};

  TTree *eventTree = new TTree("event", "event");
  eventTree->Branch("pyPt", &py_arr[0]); eventTree->Branch("pyM", &py_arr[1]); eventTree->Branch("pyZg", &py_arr[2]); eventTree->Branch("pyRg", &py_arr[3]);
  eventTree->Branch("gePt", &ge_arr[0]); eventTree->Branch("geM", &ge_arr[1]); eventTree->Branch("geZg", &ge_arr[2]); eventTree->Branch("geRg", &ge_arr[3]);
  eventTree->Branch("pyPtg", &py_arr[4]); eventTree->Branch("pyMg", &py_arr[5]); eventTree->Branch("gePtg", &ge_arr[4]); eventTree->Branch("geMg", &ge_arr[5]);
  eventTree->Branch("pyEta", &py_arr[6]); eventTree->Branch("geEta",&ge_arr[6]);
  eventTree->Branch("weight", &mc_weight); eventTree->Branch("EventID", &p_EventID);
  
  TTree *correctedTree = new TTree("corrected", "corrected");
  correctedTree->Branch("pyPt", &py_arr[0]); correctedTree->Branch("pyM", &py_arr[1]); correctedTree->Branch("pyZg", &py_arr[2]); correctedTree->Branch("pyRg", &py_arr[3]);
  correctedTree->Branch("gePt", &ge_corr[0]); correctedTree->Branch("geM", &ge_corr[1]); correctedTree->Branch("geZg", &ge_corr[2]); correctedTree->Branch("geRg", &ge_corr[3]);
  correctedTree->Branch("pyPtg", &py_arr[4]); correctedTree->Branch("pyMg", &py_arr[5]); correctedTree->Branch("gePtg", &ge_corr[4]); correctedTree->Branch("geMg", &ge_corr[5]);
  correctedTree->Branch("pyEta", &py_arr[6]); correctedTree->Branch("geEta",&ge_corr[6]);
  correctedTree->Branch("weight", &mc_weight); correctedTree->Branch("EventID", &p_EventID);

  //TEST HISTS
  TH1D* misses = new TH1D("misses","",15,5,80); TH1D* misses_cts = new TH1D("misses_cts","",15,5,80);
  
  //Hists for use in responses
  TH2D *pyMvPt = new TH2D("pyMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt = new TH2D("geMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_counts = new TH2D("pyMvPt_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_counts = new TH2D("geMvPt_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyZgvPt = new TH2D("pyZgvPt", ";z_{g};p_{T} [GeV/c]",20,0,1,15,5,80);
  TH2D *geZgvPt = new TH2D("geZgvPt", ";z_{g};p_{T} [GeV/c]",20,0,1,9,15,60);
  TH2D *pyRgvPt = new TH2D("pyRgvPt", ";R_{g};p_{T} [GeV/c]",20,0,1,15,5,80);
  TH2D *geRgvPt = new TH2D("geRgvPt", ";R_{g};p_{T} [GeV/c]",20,0,1,9,15,60);
  TH2D *pyPtgvPt = new TH2D("pyPtgvPt", ";p_{T,g} [GeV/c];p_{T} [GeV/c]",15,5,80,15,5,80);
  TH2D *gePtgvPt = new TH2D("gePtgvPt", ";p_{T,g} [GeV/c];p_{T} [GeV/c]",9,15,60,9,15,60);
  TH2D *pyMgvPt = new TH2D("pyMgvPt", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt = new TH2D("geMgvPt", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_counts = new TH2D("pyMgvPt_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_counts = new TH2D("geMgvPt_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  
  // Responses
  RooUnfoldResponse *pt_response = new RooUnfoldResponse(15,5,80,15,5,80,"pt_response","");
  RooUnfoldResponse *m_response = new RooUnfoldResponse(14,0,14,14,0,14,"m_response","");
  RooUnfoldResponse *zg_response = new RooUnfoldResponse(20,0,1,20,0,1, "zg_response","");
  RooUnfoldResponse *rg_response = new RooUnfoldResponse(20,0,1,20,0,1, "rg_response","");
  RooUnfoldResponse *ptg_response = new RooUnfoldResponse(9,15,60,15,5,80, "ptg_response","");
  RooUnfoldResponse *mg_response = new RooUnfoldResponse(14,0,14,14,0,14, "mg_response","");
  
  RooUnfoldResponse *pt_res_coarse = new RooUnfoldResponse(9,15,60,15,5,80,"pt_res_coarse","");
  
  RooUnfoldResponse *pt_m_response = new RooUnfoldResponse(geMvPt, pyMvPt, "pt_m_response");
  RooUnfoldResponse *pt_m_response_counts = new RooUnfoldResponse(geMvPt_counts, pyMvPt_counts, "pt_m_response_counts");
  RooUnfoldResponse *pt_zg_response = new RooUnfoldResponse(geZgvPt, pyZgvPt, "pt_zg_response");
  RooUnfoldResponse *pt_rg_response = new RooUnfoldResponse(geRgvPt, pyRgvPt, "pt_rg_response");
  RooUnfoldResponse *pt_ptg_response = new RooUnfoldResponse(gePtgvPt, pyPtgvPt, "pt_ptg_response");
  RooUnfoldResponse *pt_mg_response = new RooUnfoldResponse(geMgvPt, pyMgvPt, "pt_mg_response");
  RooUnfoldResponse *pt_mg_response_counts = new RooUnfoldResponse(geMgvPt_counts, pyMgvPt_counts, "pt_mg_response_counts");


  RooUnfoldResponse *m_response1520 = new RooUnfoldResponse(14,0,14,14,0,14,"m_response1520","");
  RooUnfoldResponse *m_response2025 = new RooUnfoldResponse(14,0,14,14,0,14,"m_response2025","");
  RooUnfoldResponse *m_response2530 = new RooUnfoldResponse(14,0,14,14,0,14,"m_response2530","");
  RooUnfoldResponse *m_response3040 = new RooUnfoldResponse(14,0,14,14,0,14,"m_response3040","");
  RooUnfoldResponse *m_response4060 = new RooUnfoldResponse(14,0,14,14,0,14,"m_response4060","");
  
  vector<RooUnfoldResponse*> m_response_ptbinned = {m_response1520, m_response2025, m_response2530, m_response3040, m_response4060};
  
  //vector of responses to make it easy to call ConstructResponses and fill all responses at the same time
  std::vector<RooUnfoldResponse*> res = {pt_res_coarse, m_response, pt_m_response, zg_response, rg_response, ptg_response, mg_response, pt_zg_response, pt_rg_response, pt_ptg_response, pt_mg_response, pt_m_response_counts, pt_mg_response_counts, pt_response};
  
  //Creating SoftDrop grooming object         
  contrib::SoftDrop sd(Beta,z_cut,R0);
   
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
  Selector select_det_jet_m_min = fastjet::SelectorMassMin( mass_min );
  Selector select_gen_jet_m_min = fastjet::SelectorMassMin( 0.0 );
  
  Selector sjet_gen = select_jet_rap && select_gen_jet_pt_min && select_jet_pt_max && select_gen_jet_m_min;
  Selector sjet_det = select_jet_rap && select_det_jet_pt_min && select_jet_pt_max && select_det_jet_m_min;
  
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  TString geantFilename, pythiaFilename;

  //for later use looking up PDG masses using particle PID                                                                            
  TDatabasePDG *pdg = new TDatabasePDG();
  
  // Particle containers & counters
  vector<PseudoJet> p_Particles, g_Particles, p_JetsInitial, g_JetsInitial, dummy;
  int n_accepted = 0;   int p_NJets = 0;  int g_NJets = 0;  /*int p_EventID; see above for declaration*/   int g_EventID;
  //1=inclusive, 2=lead
  int counter_debug1 = 0, counter_debug2 = 0;
  double p_wt = -1, g_wt = -1; 
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  //  while ( GEANTReader.NextEvent() ) {      //    GEANTReader    P6Reader
  for (int event = 0; event < P6Chain->GetEntries(); ++ event) {
    P6Reader.ReadEvent(event);
    GEANTReader.ReadEvent(event);
    
    //initialize values to -9999
    pyPt.clear(); pyM.clear(); pyZg.clear(); pyRg.clear();
    gePt.clear(); geM.clear(); geZg.clear(); geRg.clear();
    pyPtg.clear(); pyMg.clear(); gePtg.clear(); geMg.clear();
    pyEta.clear(); geEta.clear();
    gePt_corr.clear(); geM_corr.clear(); geZg_corr.clear(); geRg_corr.clear();
    gePtg_corr.clear(); geMg_corr.clear(); geEta_corr.clear();
    for (int i = 0; i < py_arr.size(); ++ i) {py_arr[i].clear(); ge_arr[i].clear(); ge_corr[i].clear();}
    mc_weight = -9999;
    
    g_EventID = GEANTReader.GetNOfCurrentEvent();
    p_EventID = P6Reader.GetNOfCurrentEvent();
    //cout <<"b4 potential skip" << endl;
    if ( GEANTReader.ReadEvent(p_EventID) != 1 ) {/*cout << "no corresponding geant event...skipping."<<endl;*/ continue;}   //  ENSURES BOTH DATASETS HAVE AN EVENT
    //cout <<"after potential skip" << endl;
    if ( p_EventID != g_EventID ) { cout << endl << "ERROR: READING DIFFERENT EVENTS " <<endl; }
    //if (p_EventID == 6868 || p_EventID == 5290 || g_EventID == 8999) { continue; }
    //could in reality use p_EventID both times, but this reminds me that the bad event in Pythia is 6868 and the bad event in Geant is 8999.

    p_Particles.clear(); g_Particles.clear();
    p_JetsInitial.clear(); g_JetsInitial.clear(); //clear all containers
    dummy.clear();
    
    n_accepted++;  P6Reader.PrintStatus(10);  GEANTReader.PrintStatus(10);     // Print out reader status every 10 seconds
    
    p_event = P6Reader.GetEvent();       p_header = p_event->GetHeader();           // Get the PYTHIA header and event
    g_event = GEANTReader.GetEvent();    g_header = g_event->GetHeader();           // Get GEANT event header and event
    
    p_container = P6Reader.GetOutputContainer();      // Pythia container
    g_container = GEANTReader.GetOutputContainer();      // GEANT container
    
    pythiaFilename =  P6Reader.GetInputChain()->GetCurrentFile()->GetName();	
    geantFilename =  GEANTReader.GetInputChain()->GetCurrentFile()->GetName();	

    //TEMP! CHANGE BACK IF DOESN'T FIX THINGS
    //if (((string) pythiaFilename).find("2_3_") != std::string::npos) {continue;}
 
    if (pythiaFilename != geantFilename) {std::cerr << "NOT WHAT I EXPECTED" << std::endl; exit(1);}
    
    p_wt = LookupRun12Xsec( pythiaFilename );
    g_wt = LookupRun12Xsec( geantFilename );
    if (p_wt != g_wt) {std::cerr << "WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; exit(1);}
    mc_weight = p_wt;

    //  GATHER PARTICLES
    GatherParticles ( p_container, p_sv, p_Particles, full,1, pdg);    //  Pythia particles. full = 0 signifies charged-only, 1 signifies ch+ne
    GatherParticles ( g_container, g_sv, g_Particles, full,0, pdg);    //  GEANT particles
    
    vector<PseudoJet> p_cut_Particles = spart(p_Particles); vector<PseudoJet> g_cut_Particles = spart(g_Particles); //applying constituent cuts
    
    ClusterSequence p_Cluster(p_cut_Particles, jet_def); ClusterSequence g_Cluster(g_cut_Particles, jet_def);           //  CLUSTER BOTH
    p_JetsInitial = sorted_by_pt(sjet_gen(p_Cluster.inclusive_jets())); g_JetsInitial = sorted_by_pt(sjet_det(g_Cluster.inclusive_jets()));    // EXTRACT JETS
    vector<PseudoJet> p_Jets; vector<PseudoJet> g_Jets;
    
    //Implementing a neutral energy fraction cut of 90% on inclusive det-level jets
    /*ApplyNEFSelection(p_JetsInitial, p_Jets);*/ p_Jets = p_JetsInitial; ApplyNEFSelection(g_JetsInitial, g_Jets);
    
    vector<PseudoJet> p_GroomedJets; vector<PseudoJet> g_GroomedJets;
    //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)                                                                              
    for (int i = 0; i < p_Jets.size(); ++ i) {
      p_GroomedJets.push_back(sd(p_Jets[i]));
    }
    for (int i = 0; i < g_Jets.size(); ++ i) {
      g_GroomedJets.push_back(sd(g_Jets[i]));
    }

    if (DiscardEvent(pythiaFilename, p_Jets, g_Jets)) { counter_debug1 ++; continue; }

    ConstructResponses(res, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, ge_arr, py_arr, mc_weight, m_response_ptbinned, misses, misses_cts);
    
    ge_corr = ge_arr;
    
    if (ge_arr[1].size() != 0) { //found at least one match. Should be equivalent to asking if pyPt.size() != 0 (& <=> to asking about any other observable)
      
      eventTree->Fill();
      
      //do geant correction
      for (int i = 0; i < ge_arr[1].size(); ++ i) {
	if (ge_arr[6][i] < -0.2) {
	  if (ge_arr[0][i] > 15 && ge_arr[0][i] < 20) {
	    ge_corr[1][i] /= (double) scale_lo_0->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_lo_0->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 20 && ge_arr[0][i] < 25) {
	    ge_corr[1][i] /= (double) scale_lo_1->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_lo_1->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 25 && ge_arr[0][i] < 30) {
	    ge_corr[1][i] /= (double) scale_lo_2->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_lo_2->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 30 && ge_arr[0][i] < 40) {
	    ge_corr[1][i] /= (double) scale_lo_3->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_lo_3->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 40 && ge_arr[0][i] < 60) {
	    ge_corr[1][i] /= (double) scale_lo_4->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_lo_4->Eval(ge_arr[5][i]);
	  }
	}
	if (ge_arr[6][i] > -0.2 && ge_arr[6][i] < 0.2) {
	  if (ge_arr[0][i] > 15 && ge_arr[0][i] < 20) {
	    ge_corr[1][i] /= (double) scale_mid_0->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_mid_0->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 20 && ge_arr[0][i] < 25) {
	    ge_corr[1][i] /= (double) scale_mid_1->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_mid_1->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 25 && ge_arr[0][i] < 30) {
	    ge_corr[1][i] /= (double) scale_mid_2->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_mid_2->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 30 && ge_arr[0][i] < 40) {
	    ge_corr[1][i] /= (double) scale_mid_3->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_mid_3->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 40 && ge_arr[0][i] < 60) {
	    ge_corr[1][i] /= (double) scale_mid_4->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_mid_4->Eval(ge_arr[5][i]);
	  }
	}
	if (ge_arr[6][i] > 0.2) {
	  if (ge_arr[0][i] > 15 && ge_arr[0][i] < 20) {
	    ge_corr[1][i] /= (double) scale_hi_0->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_hi_0->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 20 && ge_arr[0][i] < 25) {
	    ge_corr[1][i] /= (double) scale_hi_1->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_hi_1->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 25 && ge_arr[0][i] < 30) {
	    ge_corr[1][i] /= (double) scale_hi_2->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_hi_2->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 30 && ge_arr[0][i] < 40) {
	    ge_corr[1][i] /= (double) scale_hi_3->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_hi_3->Eval(ge_arr[5][i]);
	  }
	  if (ge_arr[0][i] > 40 && ge_arr[0][i] < 60) {
	    ge_corr[1][i] /= (double) scale_hi_4->Eval(ge_arr[1][i]);
	    ge_corr[5][i] /= (double) groom_scale_hi_4->Eval(ge_arr[5][i]);
	  }
	}
      }
      correctedTree->Fill();
    }
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    p_NJets += p_Jets.size(); g_NJets += g_Jets.size();               //  Save jet info and add jets to total
    
  }
  
  //DROPPING STATISTICS-LIMITED BINS FROM THE RESPONSE!!!
  //  DropLowStatsBins(pt_m_response, pt_m_response_counts);
  // DropLowStatsBins(pt_mg_response, pt_mg_response_counts);
  
  //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  std::cout << std::endl << std::endl << "Of " << n_accepted << " events" << std::endl;
  std::cout << p_NJets << " gen jets have been found" << std::endl;
  std::cout << g_NJets << " det jets have been found" << std::endl << std::endl;
  std::cout <<std::endl << "Writing to:  " << fout->GetName() << std::endl << std::endl;
  std::cout << "Discarded " << counter_debug1 << " events on grounds of the found jets being too much higher than the pT-hat range" << std::endl;

  eventTree->Write("event");
  correctedTree->Write("corrected");

  misses->Write(); misses_cts->Write();

  pt_res_coarse->Write(); pt_m_response->Write();
  pt_response->Write(); m_response->Write(); zg_response->Write(); rg_response->Write();
  ptg_response->Write(); mg_response->Write();
  pt_zg_response->Write(); pt_rg_response->Write(); pt_ptg_response->Write(); pt_mg_response->Write();
  m_response1520->Write(); m_response2025->Write(); m_response2530->Write(); m_response3040->Write(); m_response4060->Write();  

  pt_m_response_counts->Write(); pt_mg_response_counts->Write();
  
  fout->Close();

  return 0;
}
