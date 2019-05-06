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
  InitReader(P6Reader, P6Chain, nEvents, truth_triggerString, truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, 0.9999, false, sim_badTowers, sim_bad_run_list);
  InitReader(GEANTReader, GEANTChain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, 0.9999, false, sim_badTowers, sim_bad_run_list);
  
  TStarJetPicoEventHeader* p_header;  TStarJetPicoEvent* p_event;  TStarJetPicoEventHeader* g_header;  TStarJetPicoEvent* g_event;
  TStarJetVectorContainer<TStarJetVector> * p_container;         TStarJetVector* p_sv;
  TStarJetVectorContainer<TStarJetVector> * g_container;        TStarJetVector* g_sv;

  //  vector<double> deltaPt; vector<double> deltaM; vector<double> deltaZg; vector<double> deltaRg;
  // vector<double> ratioPt; vector<double> ratioM; vector<double> ratioZg; vector<double> ratioRg;
  vector<double> pyPt; vector<double> pyM; vector<double> pyZg; vector<double> pyRg;
  vector<double> gePt; vector<double> geM; vector<double> geZg; vector<double> geRg;
  vector<double> pyPtg; vector<double> pyMg; vector<double> gePtg; vector<double> geMg;
  vector<double> pyEta; vector<double> geEta;
  
  vector<vector<double> > py_arr = {pyPt, pyM, pyZg, pyRg, pyPtg, pyMg, pyEta};
  vector<vector<double> > ge_arr = {gePt, geM, geZg, geRg, gePtg, geMg, geEta};
  
  double mc_weight; int p_EventID;
  vector<double> pyDummy; vector<double> geDummy;
  
  TTree *eventTree = new TTree("event", "event");
  //eventTree->Branch("deltaPt", &deltaPt); eventTree->Branch("deltaM", &deltaM); eventTree->Branch("deltaZg", &deltaZg); eventTree->Branch("deltaRg", &deltaRg);
  //eventTree->Branch("ratioPt", &ratioPt); eventTree->Branch("ratioM", &ratioM); eventTree->Branch("ratioZg", &ratioZg); eventTree->Branch("ratioRg", &ratioRg);
  eventTree->Branch("pyPt", &pyPt); eventTree->Branch("pyM", &pyM); eventTree->Branch("pyZg", &pyZg); eventTree->Branch("pyRg", &pyRg);
  eventTree->Branch("gePt", &gePt); eventTree->Branch("geM", &geM); eventTree->Branch("geZg", &geZg); eventTree->Branch("geRg", &geRg);
  eventTree->Branch("pyPtg", &pyPtg); eventTree->Branch("pyMg", &pyMg); eventTree->Branch("gePtg", &gePtg); eventTree->Branch("geMg", &geMg);
  eventTree->Branch("pyEta",&pyEta); eventTree->Branch("geEta",&geEta);
  eventTree->Branch("weight", &mc_weight); eventTree->Branch("EventID", &p_EventID);
  
  double nAccepted_ge_even = 0; double nAccepted_py_even = 0; double nJets_ge_even = 0; double nJets_py_even = 0;
  double nEntries_even = 0; double nFakes_even = 0; double nMisses_even = 0; double nMatches_even = 0;
  double nAccepted_ge_odd = 0; double nAccepted_py_odd = 0; double nJets_ge_odd = 0; double nJets_py_odd = 0;
  double nEntries_odd = 0; double nFakes_odd = 0; double nMisses_odd = 0; double nMatches_odd = 0;
  
  TTree *checks = new TTree("checks","checks");
  checks->Branch("nAccepted_ge_even",&nAccepted_ge_even); checks->Branch("nAccepted_py_even",&nAccepted_py_even);
  checks->Branch("nJets_ge_even",&nJets_ge_even); checks->Branch("nJets_py_even",&nJets_py_even); checks->Branch("nEntries_even",&nEntries_even);
  checks->Branch("nFakes_even",&nFakes_even); checks->Branch("nMisses_even",&nMisses_even); checks->Branch("nMatches_even",&nMatches_even);
  checks->Branch("nAccepted_ge_odd",&nAccepted_ge_odd); checks->Branch("nAccepted_py_odd",&nAccepted_py_odd);
  checks->Branch("nJets_ge_odd",&nJets_ge_odd); checks->Branch("nJets_py_odd",&nJets_py_odd); checks->Branch("nEntries_odd",&nEntries_odd);
  checks->Branch("nFakes_odd",&nFakes_odd); checks->Branch("nMisses_odd",&nMisses_odd); checks->Branch("nMatches_odd",&nMatches_odd);
  
  //hists
  TH1D *pt_gen_odd = new TH1D("pt_gen_odd","",15,5,80); TH1D *pt_det_odd = new TH1D("pt_det_odd","",9,15,60);
  TH1D *pt_gen_even = new TH1D("pt_gen_even","",15,5,80); TH1D *pt_det_even = new TH1D("pt_det_even","",9,15,60);
  TH1D *pt_gen_odd_counts = new TH1D("pt_gen_odd_counts","",15,5,80); TH1D *pt_det_odd_counts = new TH1D("pt_det_odd_counts","",9,15,60);
  TH1D *pt_gen_even_counts = new TH1D("pt_gen_even_counts","",15,5,80); TH1D *pt_det_even_counts = new TH1D("pt_det_even_counts","",9,15,60);
  TH1D *m_gen = new TH1D("m_gen","",14,0,14); TH1D *m_det = new TH1D("m_det","",14,0,14);
  TH1D *zg_gen = new TH1D("zg_gen","",20,0,1); TH1D *zg_det = new TH1D("zg_det","",20,0,1);
  TH1D *rg_gen = new TH1D("rg_gen","",20,0,1); TH1D *rg_det = new TH1D("rg_det","",20,0,1);
  TH1D *ptg_gen = new TH1D("ptg_gen","",15,5,80); TH1D *ptg_det = new TH1D("ptg_det","",9,15,60);
  TH1D *mg_gen = new TH1D("mg_gen","",14,0,14); TH1D *mg_det = new TH1D("mg_det","",14,0,14);
  TH1D *hdummy1D = new TH1D("hdummy1D","",1,0,1);
  
  TH2D *pt_m_gen_even = new TH2D("pt_m_gen_even","",14,0,14,15,5,80); TH2D *pt_m_det_even = new TH2D("pt_m_det_even","",14,0,14,9,15,60);
  TH2D *pt_m_gen_odd = new TH2D("pt_m_gen_odd","",14,0,14,15,5,80); TH2D *pt_m_det_odd = new TH2D("pt_m_det_odd","",14,0,14,9,15,60);
  TH2D *pt_mg_gen_even = new TH2D("pt_mg_gen_even","",14,0,14,15,5,80); TH2D *pt_mg_det_even = new TH2D("pt_mg_det_even","",14,0,14,9,15,60);
  TH2D *pt_mg_gen_odd = new TH2D("pt_mg_gen_odd","",14,0,14,15,5,80); TH2D *pt_mg_det_odd = new TH2D("pt_mg_det_odd","",14,0,14,9,15,60);
  TH2D *pt_m_gen_even_counts = new TH2D("pt_m_gen_even_counts","",14,0,14,15,5,80); TH2D *pt_m_det_even_counts = new TH2D("pt_m_det_even_counts","",14,0,14,9,15,60);
  TH2D *pt_m_gen_odd_counts = new TH2D("pt_m_gen_odd_counts","",14,0,14,15,5,80); TH2D *pt_m_det_odd_counts = new TH2D("pt_m_det_odd_counts","",14,0,14,9,15,60);
  TH2D *pt_mg_gen_even_counts = new TH2D("pt_mg_gen_even_counts","",14,0,14,15,5,80); TH2D *pt_mg_det_even_counts = new TH2D("pt_mg_det_even_counts","",14,0,14,9,15,60);
  TH2D *pt_mg_gen_odd_counts = new TH2D("pt_mg_gen_odd_counts","",14,0,14,15,5,80); TH2D *pt_mg_det_odd_counts = new TH2D("pt_mg_det_odd_counts","",14,0,14,9,15,60);
  TH2D *pt_zg_gen = new TH2D("pt_zg_gen","",20,0,1,15,5,80); TH2D *pt_zg_det = new TH2D("pt_zg_det","",20,0,1,9,15,60);
  TH2D *pt_rg_gen = new TH2D("pt_rg_gen","",20,0,1,15,5,80); TH2D *pt_rg_det = new TH2D("pt_rg_det","",20,0,1,9,15,60);
  TH2D *pt_ptg_gen = new TH2D("pt_ptg_gen","",15,5,80,15,5,80); TH2D *pt_ptg_det = new TH2D("pt_ptg_det","",9,15,60,9,15,60);

  TH2D *dummygen = new TH2D("dummygen","",1,0,1,1,0,1);
  TH2D *dummydet = new TH2D("dummydet","",1,0,1,1,0,1);
  
  //hists for use in responses
  TH2D *pyMvPt = new TH2D("pyMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt = new TH2D("geMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);  
  TH2D *pyMvPt_odd = new TH2D("pyMvPt_odd",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_odd = new TH2D("geMvPt_odd",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyZgvPt = new TH2D("pyZgvPt", ";z_{g};p_{T} [GeV/c]",20,0,1,15,5,80);
  TH2D *geZgvPt = new TH2D("geZgvPt", ";z_{g};p_{T} [GeV/c]",20,0,1,9,15,60);
  TH2D *pyRgvPt = new TH2D("pyRgvPt", ";R_{g};p_{T} [GeV/c]",20,0,1,15,5,80);
  TH2D *geRgvPt = new TH2D("geRgvPt", ";R_{g};p_{T} [GeV/c]",20,0,1,9,15,60);
  TH2D *pyPtgvPt = new TH2D("pyPtgvPt", ";p_{T,g} [GeV/c];p_{T} [GeV/c]",15,5,80,15,5,80);
  TH2D *gePtgvPt = new TH2D("gePtgvPt", ";p_{T,g} [GeV/c];p_{T} [GeV/c]",15,5,80,9,15,60);
  TH2D *pyMgvPt = new TH2D("pyMgvPt", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt = new TH2D("geMgvPt", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_odd = new TH2D("pyMgvPt_odd", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_odd = new TH2D("geMgvPt_odd", ";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
 
  TH2D *pyMvPt_counts = new TH2D("pyMvPt_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_counts = new TH2D("geMvPt_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);  
  TH2D *pyMvPt_odd_counts = new TH2D("pyMvPt_odd_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_odd_counts = new TH2D("geMvPt_odd_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
 
  TH2D *pyMgvPt_counts = new TH2D("pyMgvPt_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_counts = new TH2D("geMgvPt_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_odd_counts = new TH2D("pyMgvPt_odd_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_odd_counts = new TH2D("geMgvPt_odd_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
 
 
  //responses for MC Closure test
  RooUnfoldResponse * pt_odd = new RooUnfoldResponse(9,15,60,15,5,80,"pt_odd","");
  RooUnfoldResponse * pt_even = new RooUnfoldResponse(9,15,60,15,5,80,"pt_even",""); 
  RooUnfoldResponse * m_res = new RooUnfoldResponse(14,0,14,14,0,14,"m_res",""); 
  RooUnfoldResponse * zg_res = new RooUnfoldResponse(20,0,1,20,0,1,"zg_res",""); 
  RooUnfoldResponse * rg_res = new RooUnfoldResponse(20,0,1,20,0,1,"rg_res",""); 
  RooUnfoldResponse * ptg_res = new RooUnfoldResponse(9,15,60,15,5,80,"ptg_res",""); 
  RooUnfoldResponse * mg_res = new RooUnfoldResponse(14,0,14,14,0,14,"mg_res",""); 

  RooUnfoldResponse * pt_m_response = new RooUnfoldResponse(geMvPt, pyMvPt, "pt_m_response");
  RooUnfoldResponse * pt_m_response_odd = new RooUnfoldResponse(geMvPt_odd, pyMvPt_odd, "pt_m_response_odd");
  RooUnfoldResponse * pt_zg_response = new RooUnfoldResponse(geZgvPt, pyZgvPt, "pt_zg_response");
  RooUnfoldResponse * pt_rg_response = new RooUnfoldResponse(geRgvPt, pyRgvPt, "pt_rg_response");
  RooUnfoldResponse * pt_ptg_response = new RooUnfoldResponse(gePtgvPt, pyPtgvPt, "pt_ptg_response");
  RooUnfoldResponse * pt_mg_response = new RooUnfoldResponse(geMgvPt, pyMgvPt, "pt_mg_response");
  RooUnfoldResponse * pt_mg_response_odd = new RooUnfoldResponse(geMgvPt_odd, pyMgvPt_odd, "pt_mg_response_odd");

  RooUnfoldResponse * pt_m_response_counts = new RooUnfoldResponse(geMvPt_counts, pyMvPt_counts, "pt_m_response_counts");
  RooUnfoldResponse * pt_m_response_odd_counts = new RooUnfoldResponse(geMvPt_odd_counts, pyMvPt_odd_counts, "pt_m_response_odd_counts");
  RooUnfoldResponse * pt_mg_response_counts = new RooUnfoldResponse(geMgvPt_counts, pyMgvPt_counts, "pt_mg_response_counts");
  RooUnfoldResponse * pt_mg_response_odd_counts = new RooUnfoldResponse(geMgvPt_odd_counts, pyMgvPt_odd_counts, "pt_mg_response_odd_counts");
  
  
  RooUnfoldResponse * dummy_res = new RooUnfoldResponse(1,0,1,1,0,1,"dummy_res","");

  RooUnfoldResponse * dummy_res2D = new RooUnfoldResponse(dummydet,dummygen,"dummy_res2D","");
  
  std::vector<RooUnfoldResponse*> dummy_resvec = {dummy_res,dummy_res,dummy_res,dummy_res,dummy_res};

  //vector of responses to make it easy to call ConstructResponses and fill all responses at the same time
  std::vector<RooUnfoldResponse*> res = {pt_even, m_res, pt_m_response, zg_res, rg_res, ptg_res, mg_res, pt_zg_response, pt_rg_response, pt_ptg_response, pt_mg_response, pt_m_response_counts, pt_mg_response_counts, dummy_res};
  std::vector<RooUnfoldResponse*> res_odd = {pt_odd, dummy_res, pt_m_response_odd, dummy_res,dummy_res,dummy_res,dummy_res,dummy_res2D,dummy_res2D,dummy_res2D, pt_mg_response_odd, pt_m_response_odd_counts, pt_mg_response_odd_counts, dummy_res};
  
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
  int nEvents = 0;   int p_NJets = 0;  int g_NJets = 0;  /*int p_EventID; see above for declaration*/   int g_EventID;
  //1=inclusive, 2=lead
  int counter_debug1 = 0, counter_debug2 = 0;
  double p_wt = -1, g_wt = -1; 
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVEN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  for (int event = 0; event < P6Chain->GetEntries(); ++ event) {      //    GEANTReader    P6Reader
    P6Reader.ReadEvent(event);
    GEANTReader.ReadEvent(event);
    //    if (GEANTReader.ReadEvent(event) == 0) { cout << "No Geant event! Should have only misses following this" << endl;}
	
    //initialize values to -9999
    pyPt.clear(); pyM.clear(); pyZg.clear(); pyRg.clear();
    gePt.clear(); geM.clear(); geZg.clear(); geRg.clear();
    pyPtg.clear(); pyMg.clear(); gePtg.clear(); geMg.clear();
    pyDummy.clear(); geDummy.clear();
    for (int i = 0; i < py_arr.size(); ++ i) {py_arr[i].clear(); ge_arr[i].clear();}
    mc_weight = -9999;
    
    g_EventID = GEANTReader.GetNOfCurrentEvent();
    p_EventID = P6Reader.GetNOfCurrentEvent();
    
    //    if (p_EventID == 5290 || p_EventID == 6868 || g_EventID == 8999) { continue; }
    //could in reality use p_EventID both times, but this reminds me that the bad event in Pythia is 6868 and the bad event in Geant is 8999.        
    
    //if ( p_EventID != g_EventID ) { cerr << endl << "ERROR: READING DIFFERENT EVENTS " << endl; exit(1);}
    if ( GEANTReader.ReadEvent( p_EventID ) != 1 ) continue;   //  ENSURES BOTH DATASETS HAVE AN EVENT

    p_Particles.clear(); g_Particles.clear();
    p_JetsInitial.clear(); g_JetsInitial.clear(); //clear all containers
    dummy.clear();
    
    nEvents++;  P6Reader.PrintStatus(10); // GEANTReader.PrintStatus(10);     // Print out reader status every 10 seconds
    
    p_event = P6Reader.GetEvent();       p_header = p_event->GetHeader();           // Get the PYTHIA header and event
    g_event = GEANTReader.GetEvent();    g_header = g_event->GetHeader();           // Get GEANT event header and event
       
    //BEGIN EVENS
    if (p_EventID % 2 == 0) { //even events will be used for the response
      if (GEANTReader.ReadEvent(event) != 0) { nAccepted_ge_even ++; }
      if (P6Reader.ReadEvent(event) != 0) {nAccepted_py_even ++; }
      
      p_container = P6Reader.GetOutputContainer();      // Pythia container
      g_container = GEANTReader.GetOutputContainer();      // GEANT container
      
      pythiaFilename =  P6Reader.GetInputChain()->GetCurrentFile()->GetName();	
      geantFilename =  GEANTReader.GetInputChain()->GetCurrentFile()->GetName();	
      
      //TEMP! CHANGE BACK IF DOESN'T WORK!
      //if (((string) pythiaFilename).find("2_3_") != std::string::npos) {continue; }
      
      //  if (pythiaFilename != geantFilename) {std::cerr << "NOT WHAT I EXPECTED" << std::endl; exit(1);}
      
      p_wt = LookupRun12Xsec( pythiaFilename );
      g_wt = LookupRun12Xsec( geantFilename );
      //if (p_wt != g_wt) {std::cerr << "WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; exit(1);}
      mc_weight = p_wt;
      
      //  GATHER PARTICLES
      GatherParticles ( p_container, p_sv, p_Particles, full,1,pdg);    //  Pythia particles. full = 0 signifies charged-only, 1 signifies ch+ne
      GatherParticles ( g_container, g_sv, g_Particles, full,0,pdg);    //  GEANT particles
      
      vector<PseudoJet> p_cut_Particles = spart(p_Particles); vector<PseudoJet> g_cut_Particles = spart(g_Particles); //applying constituent cuts
      
      ClusterSequence p_Cluster(p_cut_Particles, jet_def); ClusterSequence g_Cluster(g_cut_Particles, jet_def);           //  CLUSTER BOTH
      p_JetsInitial = sorted_by_pt(sjet_gen(p_Cluster.inclusive_jets())); g_JetsInitial = sorted_by_pt(sjet_det(g_Cluster.inclusive_jets()));    // EXTRACT JETS
      vector<PseudoJet> p_Jets; vector<PseudoJet> g_Jets;
      
      //Implementing a neutral energy fraction cut of 90% on inclusive det-level jets
      /*      ApplyNEFSelection(p_JetsInitial, p_Jets);*/ p_Jets = p_JetsInitial; ApplyNEFSelection(g_JetsInitial, g_Jets);
      
      vector<PseudoJet> p_GroomedJets; vector<PseudoJet> g_GroomedJets;

      if (DiscardEvent(pythiaFilename, p_Jets, g_Jets)) { counter_debug1 ++; continue; }

      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)                                    
                                          
      for (int i = 0; i < p_Jets.size(); ++ i) {
	p_GroomedJets.push_back(sd(p_Jets[i]));
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	g_GroomedJets.push_back(sd(g_Jets[i]));
      }

      nJets_py_even += p_Jets.size(); nJets_ge_even += g_Jets.size();
      
      //unmatched filling
      for (int i = 0; i < p_Jets.size(); ++ i) {
	pt_gen_even->Fill(p_Jets[i].pt(), mc_weight);
	pt_m_gen_even->Fill(p_Jets[i].m(), p_Jets[i].pt(), mc_weight);
	pt_gen_even_counts->Fill(p_Jets[i].pt());
	pt_m_gen_even_counts->Fill(p_Jets[i].m(), p_Jets[i].pt());
	pt_mg_gen_even->Fill(p_GroomedJets[i].m(), p_Jets[i].pt(), mc_weight);
	pt_mg_gen_even_counts->Fill(p_GroomedJets[i].m(), p_Jets[i].pt());

	//cout << "PT_GEN_EVEN " << p_Jets[i].pt() << endl;
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	pt_det_even->Fill(g_Jets[i].pt(), mc_weight);
	pt_m_det_even->Fill(g_Jets[i].m(), g_Jets[i].pt(), mc_weight);
	pt_det_even_counts->Fill(g_Jets[i].pt());
	pt_m_det_even_counts->Fill(g_Jets[i].m(), g_Jets[i].pt());
	pt_mg_det_even->Fill(g_GroomedJets[i].m(), g_Jets[i].pt(), mc_weight);
	pt_mg_det_even_counts->Fill(g_GroomedJets[i].m(), g_Jets[i].pt());

	//pt_det_even_pt15->Fill(g_Jets[i].pt(), mc_weight);
	//cout << "PT_DET_EVEN " << g_Jets[i].pt() << endl;
      }
      //      cout << "EVEN RESPONSE: " << endl;
      //constructing even population sample responses for the MC closure test
      ConstructResponses(res/*pt_even*/, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, ge_arr, py_arr, mc_weight/*, nEntries_even, nFakes_even, nMisses_even, nMatches_even*/, dummy_resvec, hdummy1D, hdummy1D);

    } //END EVENS
    
      //BEGIN ODDS
    else if (p_EventID % 2 != 0) { //odd events will be used for the 'data'
      if (GEANTReader.ReadEvent(event) != 0) { nAccepted_ge_odd ++; }
      if (P6Reader.ReadEvent(event) != 0) {nAccepted_py_odd ++; }
      
      p_container = P6Reader.GetOutputContainer();      // Pythia container
      g_container = GEANTReader.GetOutputContainer();      // GEANT container
      
      pythiaFilename =  P6Reader.GetInputChain()->GetCurrentFile()->GetName();
      geantFilename =  GEANTReader.GetInputChain()->GetCurrentFile()->GetName();
      
      //TEMP! CHANGE BACK IF DOESN'T WORK!
      //if (((string) pythiaFilename).find("2_3_") != std::string::npos) {continue; }
      
      //if (pythiaFilename != geantFilename) {std::cerr << "NOT WHAT I EXPECTED" << std::endl; exit(1);}
      
      p_wt = LookupRun12Xsec( pythiaFilename );
      g_wt = LookupRun12Xsec( geantFilename );
      // if (p_wt != g_wt) {std::cerr << "WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; exit(1);}
      mc_weight = p_wt;
      
      //  GATHER PARTICLES
      GatherParticles ( p_container, p_sv, p_Particles, full,1,pdg);    //  Pythia particles. full = 0 signifies charged-only, 1 signifies ch+ne
      GatherParticles ( g_container, g_sv, g_Particles, full,0,pdg);    //  GEANT particles
      
      vector<PseudoJet> p_cut_Particles = spart(p_Particles); vector<PseudoJet> g_cut_Particles = spart(g_Particles); //applying constituent cuts
      
      ClusterSequence p_Cluster(p_cut_Particles, jet_def); ClusterSequence g_Cluster(g_cut_Particles, jet_def);           //  CLUSTER BOTH
      p_JetsInitial = sorted_by_pt(sjet_gen(p_Cluster.inclusive_jets())); g_JetsInitial = sorted_by_pt(sjet_det(g_Cluster.inclusive_jets()));    // EXTRACT JETS
      vector<PseudoJet> p_Jets; vector<PseudoJet> g_Jets;
      
      //Implementing a neutral energy fraction cut of 90% on inclusive det-level jets
      /*      ApplyNEFSelection(p_JetsInitial, p_Jets);*/ p_Jets = p_JetsInitial; ApplyNEFSelection(g_JetsInitial, g_Jets);
      
      vector<PseudoJet> p_GroomedJets; vector<PseudoJet> g_GroomedJets;

      if (DiscardEvent(pythiaFilename, p_Jets, g_Jets)) { counter_debug1 ++; continue; }

      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
      for (int i = 0; i < p_Jets.size(); ++ i) {
	p_GroomedJets.push_back(sd(p_Jets[i]));
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	g_GroomedJets.push_back(sd(g_Jets[i]));
      }
      
      nJets_py_odd += p_Jets.size(); nJets_ge_odd += g_Jets.size();
      
      //unmatched filling
      for (int i = 0; i < p_Jets.size(); ++ i) {
	pt_gen_odd->Fill(p_Jets[i].pt(), mc_weight);
	pt_gen_odd_counts->Fill(p_Jets[i].pt());
	m_gen->Fill(p_Jets[i].pt(), mc_weight);
	zg_gen->Fill(p_GroomedJets[i].structure_of<SD>().symmetry(), mc_weight);
	rg_gen->Fill(p_GroomedJets[i].structure_of<SD>().delta_R(), mc_weight);
	ptg_gen->Fill(p_GroomedJets[i].pt(), mc_weight);
	mg_gen->Fill(p_GroomedJets[i].m(), mc_weight);
	
	pt_m_gen_odd->Fill(p_Jets[i].m(), p_Jets[i].pt(), mc_weight);
	pt_m_gen_odd_counts->Fill(p_Jets[i].m(), p_Jets[i].pt());
	pt_zg_gen->Fill(p_GroomedJets[i].structure_of<SD>().symmetry(), p_Jets[i].pt(), mc_weight);
	pt_rg_gen->Fill(p_GroomedJets[i].structure_of<SD>().delta_R(), p_Jets[i].pt(), mc_weight);
	pt_ptg_gen->Fill(p_GroomedJets[i].pt(), p_Jets[i].pt(), mc_weight);
	pt_mg_gen_odd->Fill(p_GroomedJets[i].m(), p_Jets[i].pt(), mc_weight);
	pt_mg_gen_odd_counts->Fill(p_GroomedJets[i].m(), p_Jets[i].pt());
	
	//cout << "PT_GEN_ODD " << p_Jets[i].pt() << endl;
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	pt_det_odd->Fill(g_Jets[i].pt(), mc_weight);
	pt_det_odd_counts->Fill(g_Jets[i].pt());
	m_det->Fill(g_Jets[i].pt(), mc_weight);
	zg_det->Fill(g_GroomedJets[i].structure_of<SD>().symmetry(), mc_weight);
	rg_det->Fill(g_GroomedJets[i].structure_of<SD>().delta_R(), mc_weight);
	ptg_det->Fill(g_GroomedJets[i].pt(), mc_weight);
	mg_det->Fill(g_GroomedJets[i].m(), mc_weight);
	
	pt_m_det_odd->Fill(g_Jets[i].m(), g_Jets[i].pt(), mc_weight);
	pt_m_det_odd_counts->Fill(g_Jets[i].m(), g_Jets[i].pt());
	pt_zg_det->Fill(g_GroomedJets[i].structure_of<SD>().symmetry(), g_Jets[i].pt(), mc_weight);
	pt_rg_det->Fill(g_GroomedJets[i].structure_of<SD>().delta_R(), g_Jets[i].pt(), mc_weight);
	pt_ptg_det->Fill(g_GroomedJets[i].pt(), g_Jets[i].pt(), mc_weight);
	pt_mg_det_odd->Fill(g_GroomedJets[i].m(), g_Jets[i].pt(), mc_weight);
	pt_mg_det_odd_counts->Fill(g_GroomedJets[i].m(), g_Jets[i].pt());

	
	//pt_det_odd_pt15->Fill(g_Jets[i].pt(), mc_weight);
	//cout << "PT_DET_ODD " << g_Jets[i].pt() << endl;
      }
      //      cout << "ODD RESPONSE: " << endl;
      //constructing odd population sample responses for the MC Closure test
      //ConstructResponses(/*pt_odd*/, g_Jets, p_Jets, geDummy, pyDummy, mc_weight, nEntries_odd, nFakes_odd, nMisses_odd, nMatches_odd);
      ConstructResponses(res_odd, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, ge_arr, py_arr, mc_weight, dummy_resvec, hdummy1D, hdummy1D);
      
       p_NJets += p_Jets.size(); g_NJets += g_Jets.size();               //  Save jet info and add jets to total 
    } //END ODDS
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  }
  //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  
  //nEntries_even = pt_even->Hresponse()->GetEntries(); nEntries_odd = pt_odd->Hresponse()->GetEntries();
  checks->Fill();

  
  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
  
  std::cout << std::endl << std::endl << "Of " << nEvents << " events" << std::endl;
  std::cout << p_NJets << " gen jets have been found in sample of odds" << std::endl;
  std::cout << g_NJets << " det jets have been found in sample of odds" << std::endl << std::endl;
  std::cout <<std::endl << "Writing to:  " << fout->GetName() << std::endl << std::endl;
  std::cout << "Discarded " << counter_debug1 << " events on grounds of the found jets being too much higher than the pT-hat range" << std::endl;

  eventTree->Write("event"); checks->Write("checks");

  pt_odd->Write(); pt_even->Write(); m_res->Write(); zg_res->Write(); rg_res->Write(); ptg_res->Write(); mg_res->Write();
  
  pt_m_response->Write(); pt_m_response_odd->Write(); pt_zg_response->Write(); pt_rg_response->Write();
  pt_ptg_response->Write(); pt_mg_response->Write(); pt_mg_response_odd->Write();
  
  pt_gen_odd->Write(); pt_det_odd->Write(); pt_gen_even->Write(); pt_det_even->Write();
  m_gen->Write(); m_det->Write(); zg_gen->Write(); zg_det->Write();
  rg_gen->Write(); rg_det->Write(); ptg_gen->Write(); ptg_det->Write();
  mg_gen->Write(); mg_det->Write();
  
  pt_m_gen_odd->Write(); pt_m_gen_even->Write(); pt_m_det_odd->Write(); pt_m_det_even->Write();
  pt_zg_gen->Write(); pt_zg_det->Write();
  pt_rg_gen->Write(); pt_rg_det->Write();
  pt_ptg_gen->Write(); pt_ptg_det->Write();
  pt_mg_gen_odd->Write(); pt_mg_gen_even->Write(); pt_mg_det_odd->Write(); pt_mg_det_even->Write();

  pt_m_response_counts->Write(); pt_mg_response_counts->Write();
  pt_m_response_odd_counts->Write(); pt_mg_response_odd_counts->Write();

  pt_gen_odd_counts->Write(); pt_det_odd_counts->Write();
  pt_gen_even_counts->Write(); pt_det_even_counts->Write();
  pt_m_gen_odd_counts->Write(); pt_m_det_odd_counts->Write();
  pt_m_gen_even_counts->Write(); pt_m_det_even_counts->Write();

  pt_mg_gen_odd_counts->Write(); pt_mg_gen_even_counts->Write(); pt_mg_det_odd_counts->Write(); pt_mg_det_even_counts->Write();
  
  fout->Close();

  return 0;
}
