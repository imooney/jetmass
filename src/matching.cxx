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
  InitReader(P6Reader, P6Chain, nEvents, truth_triggerString, truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, sim_badTowers, sim_bad_run_list);
  InitReader(GEANTReader, GEANTChain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, sim_badTowers, sim_bad_run_list);
  
  TStarJetPicoEventHeader* p_header;  TStarJetPicoEvent* p_event;  TStarJetPicoEventHeader* g_header;  TStarJetPicoEvent* g_event;
  TStarJetVectorContainer<TStarJetVector> * p_container;         TStarJetVector* p_sv;
  TStarJetVectorContainer<TStarJetVector> * g_container;        TStarJetVector* g_sv;

  const int nBins = 7;
  double edges[nBins + 1] = {5,10,15,20,25,30,40,60};
  
  //hists
  TH2D *deltaMvPyPt = new TH2D("deltaMvPyPt",";Gen. p^{jet}_{T} [GeV/c];#Delta M_{jet} (Det - Gen) / M^{gen}_{jet}",11,5,60,220,-1,1);
  TH2D *ratioMvPyPt = new TH2D("ratioMvPyPt",";M^{det}_{jet} / M^{gen}_{jet};Gen. p^{jet}_{T} [GeV/c]",51,0,2,nBins, edges);
  TH2D *deltaPtvPyPt = new TH2D("deltaPtvPyPt",";Gen. p^{jet}_{T} [GeV/c];#Delta p_{T}^{jet} (Det - Gen) / p_{T}^{gen-jet}",11,5,60,220,-1,1);
  TH2D *ratioPtvPyPt = new TH2D("ratioPtvPyPt",";p_{T}^{det-jet} / p_{T}^{gen-jet};Gen. p^{jet}_{T} [GeV/c]",51,0,2,nBins, edges);
  TH2D *deltaZgvPyPt = new TH2D("deltaZgvPyPt",";Groomed gen. p^{jet}_{T} [GeV/c];#Delta z_{g} (Det - Gen)",11,5,60,220,-1,1);
  TH2D *ratioZgvPyPt = new TH2D("ratioZgvPyPt",";z_{g}^{det} / z_{g}^{gen};Groomed gen. p^{jet}_{T} [GeV/c]",51,0,2,nBins, edges); 
  TH2D *deltaRgvPyPt = new TH2D("deltaRgvPyPt",";Groomed gen. p^{jet}_{T} [GeV/c];#Delta R_{g} (Det - Gen)",11,5,60,220,-1,1);
  TH2D *ratioRgvPyPt = new TH2D("ratioRgvPyPt",";R_{g}^{det} / R_{g}^{gen};Groomed gen. p^{jet}_{T} [GeV/c]",51,0,2,nBins, edges); 
  TH2D *pyMvPt = new TH2D("pyMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",20,0,10,15,5,80);
  TH2D *geMvPt = new TH2D("geMvPt",";M [GeV/c^{2}];p_{T} [GeV/c]",20,0,10,9,15,60);
  TH2D *pyZgvPtg = new TH2D("pyZgvPtg", ";z_{g};p_{T} [GeV/c]",20,0,1,15,5,80);
  TH2D *geZgvPtg = new TH2D("geZgvPtg", ";z_{g};p_{T} [GeV/c]",20,0,1,9,15,60);
  TH2D *pyRgvPtg = new TH2D("pyRgvPtg", ";R_{g};p_{T} [GeV/c]",20,0,1,15,5,80);
  TH2D *geRgvPtg = new TH2D("geRgvPtg", ";R_{g};p_{T} [GeV/c]",20,0,1,9,15,60);

  //note: there is only one match per event, so none of these vectors should have more than one entry. It is only done for later convenience.
  vector<double> deltaPt; vector<double> deltaM; vector<double> deltaZg; vector<double> deltaRg;
  vector<double> ratioPt; vector<double> ratioM; vector<double> ratioZg; vector<double> ratioRg;
  vector<double> pyPt; vector<double> pyM; vector<double> pyZg; vector<double> pyRg;
  vector<double> gePt; vector<double> geM; vector<double> geZg; vector<double> geRg;
  vector<double> pyPtg; vector<double> gePtg;
  double mc_weight;
  
  TTree *eventTree = new TTree("event", "event");
  eventTree->Branch("deltaPt", &deltaPt); eventTree->Branch("deltaM", &deltaM); eventTree->Branch("deltaZg", &deltaZg); eventTree->Branch("deltaRg", &deltaRg);
  eventTree->Branch("ratioPt", &ratioPt); eventTree->Branch("ratioM", &ratioM); eventTree->Branch("ratioZg", &ratioZg); eventTree->Branch("ratioRg", &ratioRg);
  eventTree->Branch("pyPt", &pyPt); eventTree->Branch("pyM", &pyM); eventTree->Branch("pyZg", &pyZg); eventTree->Branch("pyRg", &pyRg);
  eventTree->Branch("gePt", &gePt); eventTree->Branch("geM", &geM); eventTree->Branch("geZg", &geZg); eventTree->Branch("geRg", &geRg);
  eventTree->Branch("pyPtg", &pyPtg); eventTree->Branch("gePtg", &gePtg);
  eventTree->Branch("weight", &mc_weight);
  
  // Responses
  RooUnfoldResponse pt_response(60,0,60,80,0,80,"pt_response","");
  RooUnfoldResponse m_response(20,0,10,20,0,10,"m_response","");
  RooUnfoldResponse zg_response(20,0,1,20,0,1, "zg_response","");
  RooUnfoldResponse rg_response(20,0,1,20,0,1, "rg_response","");
  
  RooUnfoldResponse pt_res_coarse(9,15,60,15,5,80,"pt_res_coarse","");
  
  RooUnfoldResponse pt_m_response(geMvPt, pyMvPt, "pt_m_response");
  RooUnfoldResponse ptg_zg_response(geZgvPtg, pyZgvPtg, "ptg_zg_response");
  RooUnfoldResponse ptg_rg_response(geRgvPtg, pyRgvPtg, "ptg_rg_response");

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
  Selector select_jet_pt_min  = fastjet::SelectorPtMin( jet_ptmin );
  Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
  Selector sjet = select_jet_rap && select_jet_pt_min && select_jet_pt_max;
  
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  TString geantFilename, pythiaFilename;
  
  // Particle containers & counters
  vector<PseudoJet> p_Particles, g_Particles, p_JetsInitial, g_JetsInitial;
  int nEvents = 0;   int p_NJets = 0;  int g_NJets = 0;  int p_EventID;   int g_EventID;
  //1=inclusive, 2=lead
  int counter_debug1 = 0, counter_debug2 = 0;
  double p_wt = -1, g_wt = -1;
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( GEANTReader.NextEvent() ) {      //    GEANTReader    P6Reader
    //initialize values to -9999
    deltaPt.clear(); deltaM.clear(); deltaZg.clear(); deltaRg.clear();
    ratioPt.clear(); ratioM.clear(); ratioZg.clear(); ratioRg.clear();
    pyPt.clear(); pyM.clear(); pyZg.clear(); pyRg.clear();
    gePt.clear(); geM.clear(); geZg.clear(); geRg.clear();
    pyPtg.clear(); gePtg.clear();
    mc_weight = -9999;
    
    g_EventID = GEANTReader.GetNOfCurrentEvent();
    
    if ( P6Reader.ReadEvent( g_EventID ) != 1 ) continue;   //  ENSURES BOTH DATASETS HAVE AN EVENT
    
    p_Particles.clear(); g_Particles.clear();
    p_JetsInitial.clear(); g_JetsInitial.clear(); //clear all containers
    
    nEvents++;  P6Reader.PrintStatus(10);  GEANTReader.PrintStatus(10);     // Print out reader status every 10 seconds
    
    p_event = P6Reader.GetEvent();       p_header = p_event->GetHeader();           // Get the PYTHIA header and event
    g_event = GEANTReader.GetEvent();    g_header = g_event->GetHeader();           // Get GEANT event header and event
    
    p_EventID = P6Reader.GetNOfCurrentEvent();
    if ( p_EventID != g_EventID ) { cout << endl << "ERROR: READING DIFFERENT EVENTS " <<endl; }
    
    p_container = P6Reader.GetOutputContainer();      // Pythia container
    g_container = GEANTReader.GetOutputContainer();      // GEANT container
    
    pythiaFilename =  P6Reader.GetInputChain()->GetCurrentFile()->GetName();	
    geantFilename =  GEANTReader.GetInputChain()->GetCurrentFile()->GetName();	
    
    p_wt = LookupRun12Xsec( pythiaFilename );
    g_wt = LookupRun12Xsec( geantFilename );
    if (p_wt != g_wt) {std::cerr << "WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; exit(1);}
    mc_weight = p_wt;

    //  GATHER PARTICLES
    GatherParticles ( p_container, p_sv, p_Particles, full,1);    //  Pythia particles. full = 0 signifies charged-only, 1 signifies ch+ne
    GatherParticles ( g_container, g_sv, g_Particles, full,0);    //  GEANT particles
    
    vector<PseudoJet> p_cut_Particles = spart(p_Particles); vector<PseudoJet> g_cut_Particles = spart(g_Particles); //applying constituent cuts
    
    ClusterSequence p_Cluster(p_cut_Particles, jet_def); ClusterSequence g_Cluster(g_cut_Particles, jet_def);           //  CLUSTER BOTH
    p_JetsInitial = sorted_by_pt(sjet(p_Cluster.inclusive_jets())); g_JetsInitial = sorted_by_pt(sjet(g_Cluster.inclusive_jets()));    // EXTRACT JETS
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
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //matching leading Pythia jet to Pythia+Geant jet
    if (p_Jets.size() != 0) {
      int position = -1;
      MatchJets(g_Jets, p_Jets[0], position); //MatchJets returns with position = -1 if there is no geant jet to match in this event
      if (position == -1) {//didn't find a match
	Misses(pt_res_coarse, pt_response, m_response, pt_m_response, p_Jets[0], mc_weight);
	MissesSD(zg_response, rg_response, ptg_zg_response, ptg_rg_response, p_GroomedJets[0], mc_weight);
      }
      else { //found a match
	deltaMvPyPt->Fill(p_Jets[0].pt(),(g_Jets[position].m() - p_Jets[0].m()) / (double) p_Jets[0].m(), mc_weight);
	ratioMvPyPt->Fill(g_Jets[position].m() / (double) p_Jets[0].m(), p_Jets[0].pt(), mc_weight);
	deltaPtvPyPt->Fill(p_Jets[0].pt(),(g_Jets[position].pt() - p_Jets[0].pt()) /(double) p_Jets[0].pt(), mc_weight);
	ratioPtvPyPt->Fill(g_Jets[position].pt() / (double) p_Jets[0].pt(), p_Jets[0].pt(), mc_weight);
	//	FillMatchedTree(p_Jets[0], g_Jets[position], matchedTree, pyPt, gePt, pyM, geM, deltaPt, deltaM, Mratio, wt, mc_weight);
	ConstructResponses(pt_res_coarse, pt_response, m_response, pt_m_response, g_Jets[position], p_Jets[0], mc_weight);
	pyMvPt->Fill(p_Jets[0].m(), p_Jets[0].pt());
	geMvPt->Fill(g_Jets[position].m(), g_Jets[position].pt());
	
	deltaZgvPyPt->Fill(p_GroomedJets[0].pt(), g_GroomedJets[position].structure_of<SD>().symmetry() - p_GroomedJets[0].structure_of<SD>().symmetry(), mc_weight); //NOTE THAT WE'RE FILLING WITH THE ~GROOMED~ JET PT!
	ratioZgvPyPt->Fill(g_GroomedJets[position].structure_of<SD>().symmetry() / (double) p_GroomedJets[0].structure_of<SD>().symmetry(), p_GroomedJets[0].pt(), mc_weight); //NOTE THAT WE'RE FILLING WITH THE ~GROOMED~ JET PT!
	deltaRgvPyPt->Fill(p_GroomedJets[0].pt(), g_GroomedJets[position].structure_of<SD>().delta_R() - p_GroomedJets[0].structure_of<SD>().delta_R(), mc_weight); //NOTE THAT WE'RE FILLING WITH THE ~GROOMED~ JET PT!
	ratioRgvPyPt->Fill(g_GroomedJets[position].structure_of<SD>().delta_R() / (double) p_GroomedJets[0].structure_of<SD>().delta_R(), p_GroomedJets[0].pt(), mc_weight); //NOTE THAT WE'RE FILLING WITH THE ~GROOMED~ JET PT!
	
	ConstructResponsesSD(zg_response, rg_response, ptg_zg_response, ptg_rg_response, g_GroomedJets[position], p_GroomedJets[0], mc_weight);
	
	
	deltaPt.push_back((g_Jets[position].pt() - p_Jets[0].pt()) / (double) p_Jets[0].pt());
	deltaM.push_back((g_Jets[position].m() - p_Jets[0].m()) /(double) p_Jets[0].m());
	deltaZg.push_back(g_GroomedJets[position].structure_of<SD>().symmetry() - p_GroomedJets[0].structure_of<SD>().symmetry());
	deltaRg.push_back(g_GroomedJets[position].structure_of<SD>().delta_R() - p_GroomedJets[0].structure_of<SD>().delta_R());
	ratioPt.push_back(g_Jets[position].pt() / (double) p_Jets[0].pt());
	ratioM.push_back(g_Jets[position].m() / (double) p_Jets[0].m());
	ratioZg.push_back(g_GroomedJets[position].structure_of<SD>().symmetry() / (double) p_GroomedJets[0].structure_of<SD>().symmetry());
	ratioRg.push_back(g_GroomedJets[position].structure_of<SD>().delta_R() / (double) p_GroomedJets[0].structure_of<SD>().delta_R());
	pyPt.push_back(p_Jets[0].pt()); pyM.push_back(p_Jets[0].m()); pyZg.push_back(p_GroomedJets[0].structure_of<SD>().symmetry()); pyRg.push_back(p_GroomedJets[0].structure_of<SD>().delta_R());
	gePt.push_back(g_Jets[position].pt()); geM.push_back(g_Jets[position].m()); geZg.push_back(g_GroomedJets[position].structure_of<SD>().symmetry()); geRg.push_back(g_GroomedJets[position].structure_of<SD>().delta_R());
	pyPtg.push_back(p_GroomedJets[0].pt()); gePtg.push_back(g_GroomedJets[position].pt());
	
	eventTree->Fill();
	
      }
    }
    
    //NEED TO ALSO LOOP OVER THE GEANT JETS TO GET THE FAKE RATE
    if (g_Jets.size() != 0) {
      int position = -1;
      MatchJets(p_Jets, g_Jets[0], position);
      if (position == -1) {
	Fakes(pt_res_coarse, pt_response, m_response, pt_m_response, g_Jets[0], mc_weight);
	FakesSD(zg_response, rg_response, ptg_zg_response, ptg_rg_response, g_GroomedJets[0], mc_weight);
      }
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    p_NJets += p_Jets.size(); g_NJets += g_Jets.size();               //  Save jet info and add jets to total
    
  }
  

  //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");

  std::cout << std::endl << std::endl << "Of " << nEvents << " events" << std::endl;
  std::cout << p_NJets << " gen jets have been found" << std::endl;
  std::cout << g_NJets << " reco jets have been found" << std::endl << std::endl;
  std::cout <<std::endl << "Writing to:  " << fout->GetName() << std::endl << std::endl;

  //  matchedTree->Write("py_ge_matchedTree");
  eventTree->Write("event");
  
  pt_res_coarse.Write(); pt_m_response.Write();
  pt_response.Write(); m_response.Write(); zg_response.Write(); rg_response.Write();
  ptg_zg_response.Write(); ptg_rg_response.Write();
  //deltaMvPyPt->Write(); ratioMvPyPt->Write(); deltaPtvPyPt->Write(); ratioPtvPyPt->Write();
  //deltaZgvPyPt->Write(); ratioZgvPyPt->Write(); deltaRgvPyPt->Write(); ratioRgvPyPt->Write();
  //pyMvPt->Write(); geMvPt->Write();
  
  fout->Close();

  return 0;
}
