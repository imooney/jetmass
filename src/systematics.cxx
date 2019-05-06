//evaluates systematic uncertainties by varying values in the Geant.
//Isaac Mooney, 11/11/2018

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

  TStarJetPicoEventHeader* p_header;  TStarJetPicoEvent* p_event;  TStarJetPicoEventHeader* g_header;  TStarJetPicoEvent* g_event;
  TStarJetVectorContainer<TStarJetVector> * p_container;         TStarJetVector* p_sv;
  TStarJetVectorContainer<TStarJetVector> * g_container;        TStarJetVector* g_sv;
  
  //TFile *p8in = new TFile("~/jetmass/production/macros/hists/hists.root","READ");
  //p8in->cd();
  //TH2D *p8MvPt = (TH2D*) p8in->Get("m_v_pt_un");
  TFile *p6in = new TFile("~/jetmass/macros/hists/hists_w_o_bin_drop_R04.root","READ"); //CHANGE THIS IF YOU CHANGE THE RADIUS!!!
  //p6in->cd();
  //  TH2D *p6MvPt = (TH2D*) p6in->Get("m_v_pt_p");
  TH2D *pt_res_py2D = (TH2D*) p6in->Get("deltaPtvPyPt");
  TH2D *pt_res_ge2D = (TH2D*) p6in->Get("deltaPtvGePt");
  //TProfile *pt_res_py = (TProfile*) pt_res_py2D->ProfileX("pt_res_py",1,220);
  //TProfile *pt_res_ge = (TProfile*) pt_res_ge2D->ProfileX("pt_res_ge",1,220);

  TH2D *p8ratiop6 = (TH2D*) p6in->Get("p8MvPyPt_clone"); 
  TH2D *p8ratiop6_g = (TH2D*) p6in->Get("p8MgvPyPt_clone"); 
  TH1D *p8ratiop6_2030 = (TH1D*) p6in->Get("p8Mproj");
  TH1D *p8ratiop6_g_2030 = (TH1D*) p6in->Get("p8Mgproj"); 
  
  pt_res_py2D->SetDirectory(0);
  pt_res_ge2D->SetDirectory(0);
  p8ratiop6->SetDirectory(0); p8ratiop6_g->SetDirectory(0);
  p8ratiop6_2030->SetDirectory(0); p8ratiop6_g_2030->SetDirectory(0);
  p6in->Close();
  
  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
  fout->cd();  
  //  vector<double> deltaPt; vector<double> deltaM; vector<double> deltaZg; vector<double> deltaRg;
  // vector<double> ratioPt; vector<double> ratioM; vector<double> ratioZg; vector<double> ratioRg;
  vector<double> pyPt; vector<double> pyM;
  vector<double> gePt; vector<double> geM;
  
  double mc_weight; int p_EventID;

  //vector of vector of the doubles which will fill the branches to make it easier to pass all of them to the function
  vector<vector<double> > py_arr = {pyPt, pyM};
  vector<vector<double> > ge_arr = {gePt, geM};

  TTree *eventTree = new TTree("event", "event");
  eventTree->Branch("pyPt", &py_arr[0]); eventTree->Branch("pyM", &py_arr[1]);
  eventTree->Branch("gePt", &ge_arr[0]); eventTree->Branch("geM", &ge_arr[1]);
  eventTree->Branch("weight", &mc_weight); eventTree->Branch("EventID", &p_EventID);
  
  //TESTS!
  vector<double> Ntows_0; vector<double> Ntracks_0; double mc_weight_0; vector<double> pt_0;
  vector<double> Ntows_50; vector<double> Ntracks_50; double mc_weight_50; vector<double> pt_50;
  vector<double> Ntows_nom; vector<double> Ntracks_nom; double mc_weight_nom; vector<double> pt_nom; vector<double> M_nom; vector<double> NEF_nom;
  vector<double> Ntows_TU; vector<double> Ntracks_TU; vector<double> M_TU; vector<double> pt_TU; vector<double> NEF_TU; double mc_weight_TU;
  TTree *HC0Tree = new TTree("HC0","HC0");
  HC0Tree->Branch("Ntows",&Ntows_0); HC0Tree->Branch("Ntracks",&Ntracks_0);
  HC0Tree->Branch("mc_weight",&mc_weight_0); HC0Tree->Branch("pt",&pt_0);
  TTree *HC50Tree = new TTree("HC50","HC50");
  HC50Tree->Branch("Ntows",&Ntows_50); HC50Tree->Branch("Ntracks",&Ntracks_50);
  HC50Tree->Branch("mc_weight",&mc_weight_50); HC50Tree->Branch("pt",&pt_50);
  TTree *TUTree = new TTree("TU","TU");
  TUTree->Branch("Ntows",&Ntows_TU); TUTree->Branch("Ntracks",&Ntracks_TU);
  TUTree->Branch("M",&M_TU); TUTree->Branch("pt",&pt_TU); TUTree->Branch("NEF",&NEF_TU);
  TUTree->Branch("mc_weight",&mc_weight_TU);
  TTree *nomTree = new TTree("nom","nom");
  nomTree->Branch("Ntows",&Ntows_nom); nomTree->Branch("Ntracks",&Ntracks_nom);
  nomTree->Branch("mc_weight",&mc_weight_nom); nomTree->Branch("pt",&pt_nom);
  nomTree->Branch("M",&M_nom); nomTree->Branch("NEF",&NEF_nom);

  TH2D *dummy2D = new TH2D("dummy2D","",1,0,1,1,0,1);
  TH1D *dummy1D = new TH1D("dummy1D","",1,0,1);
  
  //Hists for use in responses
  TH2D *pyMvPt = new TH2D("pyMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt = new TH2D("geMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_TS = new TH2D("pyMvPt_TS",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_TS = new TH2D("geMvPt_TS",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_TU = new TH2D("pyMvPt_TU",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_TU = new TH2D("geMvPt_TU",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_HC50 = new TH2D("pyMvPt_HC50",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_HC50 = new TH2D("geMvPt_HC50",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_HC0 = new TH2D("pyMvPt_HC0",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_HC0 = new TH2D("geMvPt_HC0",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_DS = new TH2D("pyMvPt_DS",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_DS = new TH2D("geMvPt_DS",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_GS = new TH2D("pyMvPt_GS",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_GS = new TH2D("geMvPt_GS",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_MS = new TH2D("pyMvPt_MS",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_MS = new TH2D("geMvPt_MS",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  
  //Hists for use in responses
  TH2D *pyMgvPt = new TH2D("pyMgvPt",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt = new TH2D("geMgvPt",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_TS = new TH2D("pyMgvPt_TS",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_TS = new TH2D("geMgvPt_TS",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_TU = new TH2D("pyMgvPt_TU",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_TU = new TH2D("geMgvPt_TU",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_HC50 = new TH2D("pyMgvPt_HC50",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_HC50 = new TH2D("geMgvPt_HC50",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_HC0 = new TH2D("pyMgvPt_HC0",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_HC0 = new TH2D("geMgvPt_HC0",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_DS = new TH2D("pyMgvPt_DS",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_DS = new TH2D("geMgvPt_DS",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_GS = new TH2D("pyMgvPt_GS",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_GS = new TH2D("geMgvPt_GS",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_MS = new TH2D("pyMgvPt_MS",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_MS = new TH2D("geMgvPt_MS",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  
  //Hists for use in responses
  TH2D *pyMvPt_counts = new TH2D("pyMvPt_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_counts = new TH2D("geMvPt_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_TS_counts = new TH2D("pyMvPt_TS_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_TS_counts = new TH2D("geMvPt_TS_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_TU_counts = new TH2D("pyMvPt_TU_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_TU_counts = new TH2D("geMvPt_TU_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_HC50_counts = new TH2D("pyMvPt_HC50_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_HC50_counts = new TH2D("geMvPt_HC50_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_HC0_counts = new TH2D("pyMvPt_HC0_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_HC0_counts = new TH2D("geMvPt_HC0_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_DS_counts = new TH2D("pyMvPt_DS_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_DS_counts = new TH2D("geMvPt_DS_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_GS_counts = new TH2D("pyMvPt_GS_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_GS_counts = new TH2D("geMvPt_GS_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMvPt_MS_counts = new TH2D("pyMvPt_MS_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMvPt_MS_counts = new TH2D("geMvPt_MS_counts",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  
  //Hists for use in responses
  TH2D *pyMgvPt_counts = new TH2D("pyMgvPt_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_counts = new TH2D("geMgvPt_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_TS_counts = new TH2D("pyMgvPt_TS_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_TS_counts = new TH2D("geMgvPt_TS_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_TU_counts = new TH2D("pyMgvPt_TU_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_TU_counts = new TH2D("geMgvPt_TU_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_HC50_counts = new TH2D("pyMgvPt_HC50_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_HC50_counts = new TH2D("geMgvPt_HC50_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_HC0_counts = new TH2D("pyMgvPt_HC0_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_HC0_counts = new TH2D("geMgvPt_HC0_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_DS_counts = new TH2D("pyMgvPt_DS_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_DS_counts = new TH2D("geMgvPt_DS_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_GS_counts = new TH2D("pyMgvPt_GS_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_GS_counts = new TH2D("geMgvPt_GS_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  TH2D *pyMgvPt_MS_counts = new TH2D("pyMgvPt_MS_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D *geMgvPt_MS_counts = new TH2D("geMgvPt_MS_counts",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,9,15,60);
  
  
  // Responses  
  RooUnfoldResponse *pt_m_res_nom = new RooUnfoldResponse(geMvPt, pyMvPt, "pt_m_res_nom"); //nominal
  RooUnfoldResponse *pt_m_res_TS = new RooUnfoldResponse(geMvPt_TS, pyMvPt_TS, "pt_m_res_TS"); //tower scale
  RooUnfoldResponse *pt_m_res_TU = new RooUnfoldResponse(geMvPt_TU, pyMvPt_TU, "pt_m_res_TU"); //tracking uncertainty
  RooUnfoldResponse *pt_m_res_HC50 = new RooUnfoldResponse(geMvPt_HC50, pyMvPt_HC50, "pt_m_res_HC50"); //hadronic correction
  RooUnfoldResponse *pt_m_res_HC0 = new RooUnfoldResponse(geMvPt_HC0, pyMvPt_HC0, "pt_m_res_HC0"); //hadronic correction
  RooUnfoldResponse *pt_m_res_DS = new RooUnfoldResponse(geMvPt_DS, pyMvPt_DS, "pt_m_res_DS"); //smear detector spectrum
  RooUnfoldResponse *pt_m_res_GS = new RooUnfoldResponse(geMvPt_GS, pyMvPt_GS, "pt_m_res_GS"); //shift generator spectrum
  RooUnfoldResponse *pt_m_res_MS = new RooUnfoldResponse(geMvPt_MS, pyMvPt_MS, "pt_m_res_MS"); //shift generator mass spectrum
  
  RooUnfoldResponse *pt_mg_res_nom = new RooUnfoldResponse(geMgvPt, pyMgvPt, "pt_mg_res_nom"); //nominal
  RooUnfoldResponse *pt_mg_res_TS = new RooUnfoldResponse(geMgvPt_TS, pyMgvPt_TS, "pt_mg_res_TS"); //tower scale
  RooUnfoldResponse *pt_mg_res_TU = new RooUnfoldResponse(geMgvPt_TU, pyMgvPt_TU, "pt_mg_res_TU"); //tracking uncertainty
  RooUnfoldResponse *pt_mg_res_HC50 = new RooUnfoldResponse(geMgvPt_HC50, pyMgvPt_HC50, "pt_mg_res_HC50"); //hadronic correction
  RooUnfoldResponse *pt_mg_res_HC0 = new RooUnfoldResponse(geMgvPt_HC0, pyMgvPt_HC0, "pt_mg_res_HC0"); //hadronic correction
  RooUnfoldResponse *pt_mg_res_DS = new RooUnfoldResponse(geMgvPt_DS, pyMgvPt_DS, "pt_mg_res_DS"); //smear detector spectrum
  RooUnfoldResponse *pt_mg_res_GS = new RooUnfoldResponse(geMgvPt_GS, pyMgvPt_GS, "pt_mg_res_GS"); //shift generator spectrum
  RooUnfoldResponse *pt_mg_res_MS = new RooUnfoldResponse(geMgvPt_MS, pyMgvPt_MS, "pt_mg_res_MS"); //shift generator mass spectrum
  
  RooUnfoldResponse *pt_m_res_nom_counts = new RooUnfoldResponse(geMvPt_counts, pyMvPt_counts, "pt_m_res_nom_counts"); //nominal
  RooUnfoldResponse *pt_m_res_TS_counts = new RooUnfoldResponse(geMvPt_TS_counts, pyMvPt_TS_counts, "pt_m_res_TS_counts"); //tower scale
  RooUnfoldResponse *pt_m_res_TU_counts = new RooUnfoldResponse(geMvPt_TU_counts, pyMvPt_TU_counts, "pt_m_res_TU_counts"); //tracking uncertainty
  RooUnfoldResponse *pt_m_res_HC50_counts = new RooUnfoldResponse(geMvPt_HC50_counts, pyMvPt_HC50_counts, "pt_m_res_HC50_counts"); //hadronic correction
  RooUnfoldResponse *pt_m_res_HC0_counts = new RooUnfoldResponse(geMvPt_HC0_counts, pyMvPt_HC0_counts, "pt_m_res_HC0_counts"); //hadronic correction
  RooUnfoldResponse *pt_m_res_DS_counts = new RooUnfoldResponse(geMvPt_DS_counts, pyMvPt_DS_counts, "pt_m_res_DS_counts"); //smear detector spectrum
  RooUnfoldResponse *pt_m_res_GS_counts = new RooUnfoldResponse(geMvPt_GS_counts, pyMvPt_GS_counts, "pt_m_res_GS_counts"); //shift generator spectrum
  RooUnfoldResponse *pt_m_res_MS_counts = new RooUnfoldResponse(geMvPt_MS_counts, pyMvPt_MS_counts, "pt_m_res_MS_counts"); //shift generator mass spectrum
  
  RooUnfoldResponse *pt_mg_res_nom_counts = new RooUnfoldResponse(geMgvPt_counts, pyMgvPt_counts, "pt_mg_res_nom_counts"); //nominal
  RooUnfoldResponse *pt_mg_res_TS_counts = new RooUnfoldResponse(geMgvPt_TS_counts, pyMgvPt_TS_counts, "pt_mg_res_TS_counts"); //tower scale
  RooUnfoldResponse *pt_mg_res_TU_counts = new RooUnfoldResponse(geMgvPt_TU_counts, pyMgvPt_TU_counts, "pt_mg_res_TU_counts"); //tracking uncertainty
  RooUnfoldResponse *pt_mg_res_HC50_counts = new RooUnfoldResponse(geMgvPt_HC50_counts, pyMgvPt_HC50_counts, "pt_mg_res_HC50_counts"); //hadronic correction
  RooUnfoldResponse *pt_mg_res_HC0_counts = new RooUnfoldResponse(geMgvPt_HC0_counts, pyMgvPt_HC0_counts, "pt_mg_res_HC0_counts"); //hadronic correction
  RooUnfoldResponse *pt_mg_res_DS_counts = new RooUnfoldResponse(geMgvPt_DS_counts, pyMgvPt_DS_counts, "pt_mg_res_DS_counts"); //smear detector spectrum
  RooUnfoldResponse *pt_mg_res_GS_counts = new RooUnfoldResponse(geMgvPt_GS_counts, pyMgvPt_GS_counts, "pt_mg_res_GS_counts"); //shift generator spectrum
  RooUnfoldResponse *pt_mg_res_MS_counts = new RooUnfoldResponse(geMgvPt_MS_counts, pyMgvPt_MS_counts, "pt_mg_res_MS_counts"); //shift generator mass spectrum
  
  std::vector<RooUnfoldResponse*> res_nom = {pt_m_res_nom,pt_mg_res_nom,pt_m_res_nom_counts,pt_mg_res_nom_counts};
  std::vector<RooUnfoldResponse*> res_TS = {pt_m_res_TS,pt_mg_res_TS,pt_m_res_TS_counts,pt_mg_res_TS_counts};
  std::vector<RooUnfoldResponse*> res_TU = {pt_m_res_TU,pt_mg_res_TU,pt_m_res_TU_counts,pt_mg_res_TU_counts};
  std::vector<RooUnfoldResponse*> res_HC50 = {pt_m_res_HC50,pt_mg_res_HC50,pt_m_res_HC50_counts,pt_mg_res_HC50_counts};
  std::vector<RooUnfoldResponse*> res_HC0 = {pt_m_res_HC0,pt_mg_res_HC0,pt_m_res_HC0_counts,pt_mg_res_HC0_counts};
  std::vector<RooUnfoldResponse*> res_DS = {pt_m_res_DS,pt_mg_res_DS,pt_m_res_DS_counts,pt_mg_res_DS_counts};
  std::vector<RooUnfoldResponse*> res_GS = {pt_m_res_GS,pt_mg_res_GS,pt_m_res_GS_counts,pt_mg_res_GS_counts};
  std::vector<RooUnfoldResponse*> res_MS = {pt_m_res_MS,pt_mg_res_MS,pt_m_res_MS_counts,pt_mg_res_MS_counts};

  
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
  
  double hc = 0.9999; bool mip_correction = false;
  const int nSources = 8; //includes the nominal settings as a "systematic". 
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  //  while ( GEANTReader.NextEvent() ) {      //    GEANTReader    P6Reader
  for (int iSyst = 0; iSyst < nSources; ++ iSyst) {
    if (iSyst == 0) {cout << endl << "RUNNING WITH NOMINAL SETTINGS!" << endl << endl;}
    if (iSyst == 1) {cout << endl << "RUNNING WITH INCREASED TOWER SCALE!" << endl << endl;}
    if (iSyst == 2) {cout << endl << "RUNNING WITH DECREASED TRACKING EFFICIENCY!" << endl << endl;}
    if (iSyst == 3) {cout << endl << "RUNNING WITH 50% HADRONIC CORRECTION!" << endl << endl;}
    if (iSyst == 4) {cout << endl << "RUNNING WITH 0% HADRONIC CORRECTION!" << endl << endl;}
    if (iSyst == 5) {cout << endl << "RUNNING WITH SMEARED DETECTOR SPECTRUM!" << endl << endl;}
    if (iSyst == 6) {cout << endl << "RUNNING WITH SMEARED GENERATOR SPECTRUM!" << endl << endl;}
    if (iSyst == 7) {cout << endl << "RUNNING WITH SMEARED GENERATOR ~MASS~ SPECTRUM!" << endl << endl;}
    
    p_NJets = 0; g_NJets = 0; n_accepted = 0;
    //set parameters back to their nominal values after the previous iteration changed them.
    hc = 0.9999; mip_correction = false;
    
    //change the nominal values
    if (iSyst == 3) { //this means the systematic we're examining is the hadronic correction. Set it to 50%
      hc = 0.5;
    }
    if (iSyst == 4) { //Set it to 0%
      hc = 0.0; mip_correction = true;
    }
    
    //initialize the readers!
    InitReader(P6Reader, P6Chain, nEvents, truth_triggerString, truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, hc, mip_correction, sim_badTowers, sim_bad_run_list);
    InitReader(GEANTReader, GEANTChain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, hc, mip_correction, sim_badTowers, sim_bad_run_list);
    
    for (int event = 0; event < P6Chain->GetEntries(); ++ event) {
      P6Reader.ReadEvent(event);
      GEANTReader.ReadEvent(event);
      
      //initialize values to -9999
      Ntows_0.clear(); Ntracks_0.clear(); pt_0.clear();
      Ntows_50.clear(); Ntracks_50.clear(); pt_50.clear();
      Ntows_nom.clear(); Ntracks_nom.clear(); pt_nom.clear(); M_nom.clear(); NEF_nom.clear();
      Ntows_TU.clear(); Ntracks_TU.clear(); M_TU.clear(); pt_TU.clear(); NEF_TU.clear();
      pyPt.clear(); pyM.clear();
      gePt.clear(); geM.clear();
      for (int i = 0; i < py_arr.size(); ++ i) {py_arr[i].clear(); ge_arr[i].clear();}
      mc_weight = -9999; mc_weight_0 = -9999; mc_weight_50 = -9999; mc_weight_nom = -9999; mc_weight_TU = -9999;
      
      g_EventID = GEANTReader.GetNOfCurrentEvent();
      p_EventID = P6Reader.GetNOfCurrentEvent();

      if ( GEANTReader.ReadEvent(p_EventID) != 1 ) {continue;} //!ENSURES BOTH DATASETS HAVE AN EVENT
      
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
      
      p_wt = LookupRun12Xsec( pythiaFilename );
      g_wt = LookupRun12Xsec( geantFilename );
      // if (p_wt != g_wt) {std::cerr << "WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; exit(1);}
      mc_weight = p_wt;

      if (iSyst == 1) {//varying the gain of the towers
	for (int i = 0; i < g_container->GetEntries(); ++ i) {
	  g_sv = g_container->Get(i);
	  if (!(g_sv->IsCharged())) {
	    //	    cout << "Pre-change: " << g_sv->E() << " " << g_sv->Eta() << " " << g_sv->Phi() << " " << g_sv->M() << endl;
	    double Enew = 1.038*g_sv->E();
	    g_sv->SetE(Enew);
	    //	    g_sv->SetPtEtaPhiM(sqrt(Etnew*Etnew - g_sv->M()*g_sv->M()), g_sv->Eta(), g_sv->Phi(), g_sv->M());
	    //cout << "Post-change: " << g_sv->E() << " " << g_sv->Eta() << " " << g_sv->Phi() << " " << g_sv->M() << endl;
	    //cout << "Percent diff: " << ((g_sv->E()/(double)(Enew/(double)1.038)) - 1)*100 << "%" << endl;
	  }
	}
      }     

      //  GATHER PARTICLES
      GatherParticles ( p_container, p_sv, p_Particles, full,1, pdg);    //  Pythia particles. full = 0 signifies charged-only, 1 signifies ch+ne
      GatherParticles ( g_container, g_sv, g_Particles, full,0, pdg);    //  GEANT particles
            
      if (iSyst == 2) {//varying the tracking efficiency randomly by 4%
	double effic_num;
	for (int i = 0; i < g_Particles.size(); ++ i) {
	  if (g_Particles[i].user_index() != 0) {
	    effic_num = gRandom->Uniform(0.0, 1.0);
	    if (effic_num > 0.96) {
	      g_Particles.erase(g_Particles.begin() + i);
	      i --; //need to account for the shrinking of the list.
	    }
	  }
	}
      }
      
      vector<PseudoJet> p_cut_Particles = spart(p_Particles); vector<PseudoJet> g_cut_Particles = spart(g_Particles); //applying constituent cuts
      
      ClusterSequence p_Cluster(p_cut_Particles, jet_def); ClusterSequence g_Cluster(g_cut_Particles, jet_def);           //  CLUSTER BOTH
      p_JetsInitial = sorted_by_pt(sjet_gen(p_Cluster.inclusive_jets())); g_JetsInitial = sorted_by_pt(sjet_det(g_Cluster.inclusive_jets()));    // EXTRACT JETS
      vector<PseudoJet> p_Jets; vector<PseudoJet> g_Jets; vector<double> prior_adjust;
      
      //Implementing a neutral energy fraction cut of 90% on inclusive det-level jets
      p_Jets = p_JetsInitial; ApplyNEFSelection(g_JetsInitial, g_Jets);

      if (DiscardEvent(pythiaFilename, p_Jets, g_Jets)) { counter_debug1 ++; continue; }
      /*    
      if (iSyst == 5) { //detector smearing
	double ptsmear = 0;
	for (int i = 0; i < g_Jets.size(); ++ i) {
	  //cout << "PRE-CALL" << endl;
	  double res_for_this_jet = pt_res_ge->GetBinContent(pt_res_ge2D->GetXaxis()->FindBin(g_Jets[i].pt()));
	  //cout << "POST-CALL" << endl;
	  //	  cout << "det res: " << res_for_this_jet << " for jet with pT " << g_Jets[i].pt() << endl;
	  ptsmear = gRandom->Gaus(0,fabs(res_for_this_jet*g_Jets[i].pt()));
	  prior_adjust.push_back(ptsmear);
	  //	  cout << "POST-2ND-CALL" <<endl;
	  //cout << "smearing with " << ptsmear << endl;
	  //g_Jets[i].reset_momentum_PtYPhiM(g_Jets[i].pt() + ptsmear, g_Jets[i].rap(), g_Jets[i].phi(), g_Jets[i].m());
	  //cout << "results in jet with " << g_Jets[i].pt() << endl;
	}
	cout << prior_adjust.size() << " = " << g_Jets.size() << "?" << endl;
      }

      if (iSyst == 6) {
	double ptsmear = 0;
	for (int i = 0; i < p_Jets.size(); ++ i) {
	  double res_for_this_jet = pt_res_py->GetBinContent(pt_res_py->GetXaxis()->FindBin(p_Jets[i].pt()));
	  //cout << "gen res: " << res_for_this_jet << " for jet with pT " << p_Jets[i].pt() << endl;
	  ptsmear = fabs(gRandom->Gaus(0,fabs(res_for_this_jet*p_Jets[i].pt())));
	  prior_adjust.push_back(ptsmear);
	  //cout << "smearing with " << ptsmear << endl;
	  //p_Jets[i].reset_momentum_PtYPhiM(p_Jets[i].pt() - ptsmear, p_Jets[i].rap(), p_Jets[i].phi(), p_Jets[i].m());
	  //cout << "results in jet with " << p_Jets[i].pt() << endl;
	}
	cout << prior_adjust.size() << " = " << p_Jets.size() << "?" << endl;
      }
      
      if (iSyst == 7) {
	double msmear = 0;
	for (int i = 0; i < p_Jets.size(); ++ i) {
	  msmear = p8ratiop6->GetBinContent(p8ratiop6->GetXaxis()->FindBin(p_Jets[i].m()),p8ratiop6->GetYaxis()->FindBin(p_Jets[i].pt()));
	  prior_adjust.push_back(msmear);
	  //  cout << "smearing mass with " << msmear << endl;
	  //p_Jets[i].reset_momentum_PtYPhiM(p_Jets[i].pt(), p_Jets[i].rap(), p_Jets[i].phi(), msmear*p_Jets[i].m());
	}
	cout << prior_adjust.size() << " = " << p_Jets.size() << "?" << endl;
      }
*/
      vector<PseudoJet> p_GroomedJets; vector<PseudoJet> g_GroomedJets;
      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)                                                                              
      for (int i = 0; i < p_Jets.size(); ++ i) {
	p_GroomedJets.push_back(sd(p_Jets[i]));
      }
      for (int i = 0; i < g_Jets.size(); ++ i) {
	g_GroomedJets.push_back(sd(g_Jets[i]));
      }
      

      //QA HISTOGRAMS & TREES~~~
      if (iSyst == 0) {mc_weight_nom = mc_weight;}
      if (iSyst == 2) {mc_weight_TU = mc_weight;}
      if (iSyst == 3) {mc_weight_50 = mc_weight;}
      if (iSyst == 4) {mc_weight_0 = mc_weight;}
      for (int i = 0; i < g_Jets.size(); ++ i) {
	double nTows = 0; double nTracks = 0; double towsum = 0; double ptsum = 0;
	for (int j = 0; j < g_Jets[i].constituents().size(); ++ j) {
	  if (g_Jets[i].constituents()[j].user_index() == 0) { //have a tower
	    ++ nTows; 
	    towsum += g_Jets[i].constituents()[j].pt();
	  }
	  else if (g_Jets[i].constituents()[j].user_index() != 0) { //have a track
	    ++ nTracks;
	  }
	  ptsum += g_Jets[i].constituents()[j].pt();
	}
	if (iSyst == 0) {
	  Ntows_nom.push_back(nTows); Ntracks_nom.push_back(nTracks); pt_nom.push_back(g_Jets[i].pt());
	  M_nom.push_back(g_Jets[i].m()); NEF_nom.push_back(towsum / (double) ptsum);
	}
	if (iSyst == 2) {
	  Ntows_TU.push_back(nTows); Ntracks_TU.push_back(nTracks); M_TU.push_back(g_Jets[i].m()); pt_TU.push_back(g_Jets[i].pt());
	  NEF_TU.push_back(towsum / (double) ptsum);
	}
	if (iSyst == 3) {
	  Ntows_50.push_back(nTows); Ntracks_50.push_back(nTracks); pt_50.push_back(g_Jets[i].pt()); 
	}
	if (iSyst == 4) {
	  Ntows_0.push_back(nTows); Ntracks_0.push_back(nTracks); pt_0.push_back(g_Jets[i].pt()); 
	}
	
      }
      if (g_Jets.size() != 0) {
	if (iSyst == 0) {nomTree->Fill();}
	if (iSyst == 2) {TUTree->Fill();}
	if (iSyst == 3) {HC50Tree->Fill();}
	if (iSyst == 4) {HC0Tree->Fill();}
      }
      //~~~

      
      //ConstructResponses(res, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, ge_arr, py_arr, mc_weight, m_response_ptbinned);
      if (iSyst == 0) {ConstructSystematicsResponses(res_nom, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, mc_weight, dummy2D, dummy2D, dummy1D, dummy1D, iSyst);}
      if (iSyst == 1) {ConstructSystematicsResponses(res_TS, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, mc_weight, dummy2D, dummy2D, dummy1D, dummy1D, iSyst);}
      if (iSyst == 2) {ConstructSystematicsResponses(res_TU, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, mc_weight, dummy2D, dummy2D, dummy1D, dummy1D, iSyst);}
      if (iSyst == 3) {ConstructSystematicsResponses(res_HC50, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, mc_weight, dummy2D, dummy2D, dummy1D, dummy1D, iSyst);}
      if (iSyst == 4) {ConstructSystematicsResponses(res_HC0, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, mc_weight, dummy2D, dummy2D, dummy1D, dummy1D, iSyst);}
      if (iSyst == 5) {ConstructSystematicsResponses(res_DS, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, mc_weight, pt_res_ge2D, pt_res_ge2D, dummy1D, dummy1D, iSyst);}
      if (iSyst == 6) {ConstructSystematicsResponses(res_GS, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, mc_weight, pt_res_py2D, pt_res_py2D, dummy1D, dummy1D, iSyst);}
      if (iSyst == 7) {ConstructSystematicsResponses(res_MS, g_Jets, p_Jets, g_GroomedJets, p_GroomedJets, mc_weight, dummy2D, dummy2D, p8ratiop6_2030, p8ratiop6_g_2030, iSyst);}
      
      if (ge_arr[1].size() != 0) { //found at least one match. Should be equivalent to asking if pyPt.size() != 0 (& <=> to asking about any other observable)
	//  eventTree->Fill();
      }
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      p_NJets += p_Jets.size(); g_NJets += g_Jets.size();               //  Save jet info and add jets to total
          
    }
    //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    
    std::cout << std::endl << std::endl << "Of " << n_accepted << " events" << std::endl;
    std::cout << p_NJets << " gen jets have been found" << std::endl;
    std::cout << g_NJets << " det jets have been found" << std::endl << std::endl;
    std::cout <<std::endl << "Writing to:  " << fout->GetName() << std::endl << std::endl;
    fout->cd();
    std::cout << "Discarded " << counter_debug1 << " events on grounds of the found jets being too much higher than the pT-hat range" << std::endl;
    
    //eventTree->Write("event");
    
    //First drop low statistics bins from the responses, then write to file.    
    if (iSyst == 0) {//DropLowStatsBins(pt_m_res_nom, pt_m_res_nom_counts); DropLowStatsBins(pt_mg_res_nom, pt_mg_res_nom_counts);
      pt_m_res_nom->Write(); pt_mg_res_nom->Write(); nomTree->Write();
      pt_m_res_nom_counts->Write(); pt_mg_res_nom_counts->Write();
    }
    if (iSyst == 1) {//DropLowStatsBins(pt_m_res_TS, pt_m_res_TS_counts); DropLowStatsBins(pt_mg_res_TS, pt_mg_res_TS_counts);
      pt_m_res_TS->Write(); pt_mg_res_TS->Write();
      pt_m_res_TS_counts->Write(); pt_mg_res_TS_counts->Write();
    }
    if (iSyst == 2) {//DropLowStatsBins(pt_m_res_TU, pt_m_res_TU_counts); DropLowStatsBins(pt_mg_res_TU, pt_mg_res_TU_counts);
      pt_m_res_TU->Write(); pt_mg_res_TU->Write(); TUTree->Write();
      pt_m_res_TU_counts->Write(); pt_mg_res_TU_counts->Write();
    }
    if (iSyst == 3) {//DropLowStatsBins(pt_m_res_HC50, pt_m_res_HC50_counts); DropLowStatsBins(pt_mg_res_HC50, pt_mg_res_HC50_counts);
      pt_m_res_HC50->Write(); pt_mg_res_HC50->Write(); HC50Tree->Write();
      pt_m_res_HC50_counts->Write(); pt_mg_res_HC50_counts->Write();
    }
    if (iSyst == 4) {//DropLowStatsBins(pt_m_res_HC0, pt_m_res_HC0_counts); DropLowStatsBins(pt_mg_res_HC0, pt_mg_res_HC0_counts);
      pt_m_res_HC0->Write(); pt_mg_res_HC0->Write(); HC0Tree->Write();
      pt_m_res_HC0_counts->Write(); pt_mg_res_HC0_counts->Write();

    }
    if (iSyst == 5) {//DropLowStatsBins(pt_m_res_DS, pt_m_res_DS_counts); DropLowStatsBins(pt_mg_res_DS, pt_mg_res_DS_counts);
      pt_m_res_DS->Write(); pt_mg_res_DS->Write();
      pt_m_res_DS_counts->Write(); pt_mg_res_DS_counts->Write();
    }
    if (iSyst == 6) {//DropLowStatsBins(pt_m_res_GS, pt_m_res_GS_counts); DropLowStatsBins(pt_mg_res_GS, pt_mg_res_GS_counts);
      pt_m_res_GS->Write(); pt_mg_res_GS->Write();
      pt_m_res_GS_counts->Write(); pt_mg_res_GS_counts->Write();
    }
    if (iSyst == 7) {
      pt_m_res_MS->Write(); pt_mg_res_MS->Write();
      pt_m_res_MS_counts->Write(); pt_mg_res_MS_counts->Write();
      //p8ratiop6->Write();
    }
  }
  //fout->Flush();
  fout->Close();
  
  return 0;
}
