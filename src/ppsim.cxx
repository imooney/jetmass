//  Veronica Verkest        May 13, 2018
//  Compare:   p6  VS  p6+efficiency  VS  p6+GEANT
//  Functions in src/functions.cxx
//  Parameters in src/parameters.cxx
//  Adapted by Isaac Mooney June, 2018 for jet mass analysis

#include "params.hh"
#include "funcs.hh"
#include "TStarJetPicoDefinitions.h"

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
  InitReader(P6Reader, P6Chain, nEvents, truth_triggerString, truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, sim_badTowers, sim_bad_run_list);
  InitReader(GEANTReader, GEANTChain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, sim_badTowers, sim_bad_run_list);
  
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
    if (i == 0 && full == 1) {continue;} if (i == 1 && full == 0) {continue;}
    for (int j = 0; j < flag_j.size(); ++ j) {
      for (int k = 0; k < flag_k.size(); ++ k) {
	for (int l = 0; l < flag_l.size(); ++ l) {
	  hists1D.add(("m_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]+"_"+flag_l[l]).c_str(),"",20,0,10); //mass
	  hists2D.add(("m_v_pt_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]+"_"+flag_l[l]).c_str(),"",20,0,10,11,5,60); //mass vs. pT
	  hists3D.add(("PtEtaPhi_"+flag_i[i]+"_"+flag_j[j]+"_"+flag_k[k]+"_"+flag_l[l]).c_str(),"",11,5,60,30,-0.6,0.6,50,0,2*Pi);
	}
      }
    }
  }

  const unsigned nDim = 5;
  int bins[nDim] = {20, 20, 20, 11, 11};
  double min[nDim] = {0,0,0,5,5};
  double max[nDim] = {1,10,1,60,60};
  THnSparse * SDnD_py = new THnSparseD("zg_mg_thetag_ptg_pt_full_incl_sd_py", "", nDim, bins, min, max);
  THnSparse * SDnD_ge = new THnSparseD("zg_mg_thetag_ptg_pt_full_incl_sd_ge","", nDim, bins, min, max);
  SDnD_py->Sumw2(); SDnD_ge->Sumw2();
  
  //matching histograms
  TH2D * pt_response = new TH2D("pt_response", ";Py+Ge p^{jet}_{T} [GeV/c];Py p^{jet}_{T} [GeV/c]", 60,0,60,80,0,80);
  
  //TESTS!
  TH1D * py_dPhi_trig_rec = new TH1D("py_dPhi_trig_rec",";#Delta #phi;arb.", 28, -Pi - 0.4, Pi + 0.4); //defined as trigger - recoil
  TH1D * ge_dPhi_trig_rec = new TH1D("ge_dPhi_trig_rec",";#Delta #phi;arb.", 28, -Pi - 0.4, Pi + 0.4); //defined as trigger - recoil
  TH3D * py_PtEtaPhi_tracks = new TH3D("py_PtEtaPhi_tracks",";p^{track}_{T} [GeV/c]; #eta; #phi", 80, 0, 80, 30, -0.6,0.6,50,0,2*Pi);
  TH3D * ge_PtEtaPhi_tracks = new TH3D("ge_PtEtaPhi_tracks",";p^{track}_{T} [GeV/c]; #eta; #phi", 80, 0, 80, 30, -0.6,0.6,50,0,2*Pi);
  TH2D * ch_frac_v_pt_py = new TH2D("ch_frac_v_pt_py",";charged fraction; p_{T}^{jet} [GeV/c]",10,0,1,11,5,60);
  TH2D * ch_frac_v_pt_ge = new TH2D("ch_frac_v_pt_ge",";charged fraction; p_{T}^{jet} [GeV/c]",10,0,1,11,5,60);
  TH2D *tow_id_v_e = new TH2D("tow_id_v_e",";tower ID; tower E_{T} [GeV]",4800,1,4801,140,0,140);


  double jetPt_incl, jetEta_incl, jetPhi_incl, jetM_incl, jetE_incl, consPt_incl, consEta_incl, consPhi_incl, consM_incl, consE_incl, pythia_wt, geant_wt;
  int nCons_incl, cons_dummy;
  double jetPt_lead, jetEta_lead, jetPhi_lead, jetM_lead, jetE_lead, consPt_lead, consEta_lead, consPhi_lead, consM_lead, consE_lead;
  double jetPt_sublead, jetEta_sublead, jetPhi_sublead, jetM_sublead, jetE_sublead;
  int nCons_lead, nCons_sublead;
  double deltaM, deltaPt, geM, pyM, gePt, pyPt;

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
  TTree *matchedTree = new TTree("py_ge_matchedTree","py_ge_matchedTree");
				 
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

  matchedTree->Branch("DeltaM", &deltaM); matchedTree->Branch("DeltaPt", &deltaPt); matchedTree->Branch("GePt", &gePt); matchedTree->Branch("PyPt",&pyPt);
  matchedTree->Branch("GeM", &geM); matchedTree->Branch("PyM", &pyM); matchedTree->Branch("py_weight", &pythia_wt); matchedTree->Branch("ge_weight", &geant_wt);
 
  ch_matchedTree->Branch("DeltaM", &deltaM); ch_matchedTree->Branch("DeltaPt", &deltaPt); ch_matchedTree->Branch("GePt", &gePt); ch_matchedTree->Branch("PyPt",&pyPt);
  ch_matchedTree->Branch("GeM", &geM); ch_matchedTree->Branch("PyM", &pyM); ch_matchedTree->Branch("py_weight", &pythia_wt); ch_matchedTree->Branch("ge_weight", &geant_wt);


  
  //Creating SoftDrop grooming object                                                                                                                        
  contrib::SoftDrop sd(beta,z_cut,R0);
  cout << "SoftDrop groomer is: " << sd.description() << endl;
  
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
  double p_wt, g_wt;
  TString geantFilename, pythiaFilename;
  
  // Particle containers & counters
  vector<PseudoJet> p_Particles, g_Particles, ch_p_Particles, ch_g_Particles, p_Jets, g_Jets, ch_p_Jets, ch_g_Jets, p_GroomedJets, g_GroomedJets, ch_p_GroomedJets, ch_g_GroomedJets;//, p_Cons, g_Cons;
  int nEventsPythia = 0, nEventsGeant = 0;   int p_NJets = 0;  int g_NJets = 0;  int p_EventID;   int g_EventID;
  //1=inclusive, 2=lead
  int counter_debug1 = 0, counter_debug2 = 0;
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  //TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  while ( GEANTReader.NextEvent() ) {
    g_EventID = GEANTReader.GetNOfCurrentEvent(); p_EventID = P6Reader.GetNOfCurrentEvent();
    
    g_event = GEANTReader.GetEvent();    g_header = g_event->GetHeader();
    p_event = P6Reader.GetEvent(); p_header = p_event->GetHeader();

    //TEMP!!!!!!!!!!!!!!!!!                                                                                                                                           
    TStarJetPicoEvent * current_event = GEANTReader.GetEvent();
    //      cout << "event " << reader.GetNOfCurrentEvent() << " has these tows" << endl;                                                                             
    for (int i = 0; i < current_event->GetTowers()->GetEntries(); ++ i) {
      tow_id_v_e->Fill(current_event->GetTower(i)->GetId(), current_event->GetTower(i)->GetEnergy());
      //cout << current_event->GetTower(i)->GetId() << " " << current_event->GetTower(i)->GetEnergy() << endl;                                                        
    }    


    while (p_header->GetEventId() != g_header->GetEventId()) {
      //DO PYTHIA STUFF!!!
      
       p_Particles.clear(); ch_p_Particles.clear();
       p_Jets.clear(); ch_p_Jets.clear();
       p_GroomedJets.clear(); ch_p_GroomedJets.clear();
       
       nEventsPythia ++;  P6Reader.PrintStatus(10);     // Print out reader status every 10 seconds
       p_event = P6Reader.GetEvent();       p_header = p_event->GetHeader();           // Get the PYTHIA header and event   
       p_EventID = P6Reader.GetNOfCurrentEvent();
       p_container = P6Reader.GetOutputContainer();      // Pythia container
       pythiaFilename =  P6Reader.GetInputChain()->GetCurrentFile()->GetName();	
       p_wt = LookupRun12Xsec( pythiaFilename );
       
       //  GATHER PARTICLES
       GatherParticles ( p_container, p_sv, p_Particles, 1, 1);    //  Pythia particles
       GatherParticles ( p_container, p_sv, ch_p_Particles, 0, 1); // first bool flag: 0 signifies charged-only, 1 signifies ch+ne. second bool flag: 0 = geant, 1 = pythia
       
       //TEST
       for (int i = 0; i < p_Particles.size(); ++ i) {
	 py_PtEtaPhi_tracks->Fill(p_Particles[i].pt(), p_Particles[i].eta(), p_Particles[i].phi(), p_wt);
       }

       vector<PseudoJet> p_cut_Particles = spart(p_Particles); //applying constituent cuts
       vector<PseudoJet> ch_p_cut_Particles = spart(ch_p_Particles);
       
       ClusterSequence p_Cluster(p_cut_Particles, jet_def);           //  CLUSTER
       ClusterSequence ch_p_Cluster(ch_p_cut_Particles, jet_def);
       p_Jets = sorted_by_pt(sjet(p_Cluster.inclusive_jets()));    // EXTRACT JETS
       ch_p_Jets = sorted_by_pt(sjet(p_Cluster.inclusive_jets()));
       //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
       for (int i = 0; i < p_Jets.size(); ++ i) {
	 p_GroomedJets.push_back(sd(p_Jets[i]));
       }
       for (int i = 0; i < ch_p_Jets.size(); ++ i) {
	 ch_p_GroomedJets.push_back(sd(ch_p_Jets[i]));
       }
       
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
	   //TEMP                                                                                                                                                      
	   int numch_py = 0;
	   for (int k = 0; k < p_Jets[j].constituents().size(); ++ k) {
	     if (p_Jets[j].constituents()[k].user_index() != 0) {
	       numch_py ++;
	     }
	   }
	   ch_frac_v_pt_py->Fill(numch_py/(double) p_Jets[j].constituents().size(),p_Jets[j].pt(), p_wt);
	 }
       }
       
       //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
       //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HISTOS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
       
       if(p_Jets.size() != 0) {
	 FillHists(hists1D, hists2D, hists3D, "full", "_py", p_Jets, p_wt);
	 for (int j = 0; j < p_Jets.size(); ++ j) {
	   double val_list[nDim] = {p_GroomedJets[j].structure_of<contrib::SoftDrop>().symmetry(),p_GroomedJets[j].m(),p_GroomedJets[j].structure_of<contrib::SoftDrop>().delta_R(),p_GroomedJets[j].pt(), p_Jets[j].pt()}; //Groomed jet not guaranteed to be highest pT even though ungroomed one is
	   SDnD_py->Fill(val_list); 
	 }
       }
       if(ch_p_Jets.size() != 0) {
	 FillHists(hists1D, hists2D, hists3D, "ch", "_py", ch_p_Jets, p_wt);
       }
       p_NJets += p_Jets.size();               //  Save jet info and add jets to total

       P6Reader.NextEvent();        
    }
    //DO BOTH PYTHIA & GEANT STUFF
    g_EventID = GEANTReader.GetNOfCurrentEvent();
     
    if ( P6Reader.ReadEvent( g_EventID ) != 1 ) continue;   //  ENSURES BOTH DATASETS HAVE AN EVENT

    p_Particles.clear(); g_Particles.clear();
    ch_p_Particles.clear(); ch_g_Particles.clear();
    p_Jets.clear(); g_Jets.clear(); 
    ch_p_Jets.clear(); ch_g_Jets.clear();
    p_GroomedJets.clear(); g_GroomedJets.clear();
    ch_p_GroomedJets.clear(); ch_g_GroomedJets.clear();

    nEventsGeant ++; nEventsPythia ++; P6Reader.PrintStatus(10);  GEANTReader.PrintStatus(10);     // Print out reader status every 10 seconds

    p_event = P6Reader.GetEvent();       p_header = p_event->GetHeader();           // Get the PYTHIA header and event
    g_event = GEANTReader.GetEvent();    g_header = g_event->GetHeader();           // Get GEANT event header and event
    
    p_EventID = P6Reader.GetNOfCurrentEvent();
    if ( p_EventID != g_EventID ) { std::cerr << endl << "ERROR: READING DIFFERENT EVENTS " << std::endl; }
    
    p_container = P6Reader.GetOutputContainer();      // Pythia container
    g_container = GEANTReader.GetOutputContainer();      // GEANT container

    pythiaFilename =  P6Reader.GetInputChain()->GetCurrentFile()->GetName();	
    geantFilename =  GEANTReader.GetInputChain()->GetCurrentFile()->GetName();	
    
    p_wt = LookupRun12Xsec( pythiaFilename );
    g_wt = LookupRun12Xsec( geantFilename );
    if (p_wt != g_wt) {std::cout<< "WRONG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;}

    //  GATHER PARTICLES
    GatherParticles ( p_container, p_sv, p_Particles, 1, 1);    //  Pythia particles: second bool flag: pythia = 1, geant = 0
    GatherParticles ( p_container, p_sv, ch_p_Particles, 0, 1); //0 signifies charged-only, 1 signifies ch+ne
    GatherParticles ( g_container, g_sv, g_Particles, 1, 0);    //  GEANT particles
    GatherParticles ( g_container, g_sv, ch_g_Particles, 0, 0);

    //TEST
    for (int i = 0; i < p_Particles.size(); ++ i) {
      py_PtEtaPhi_tracks->Fill(p_Particles[i].pt(), p_Particles[i].eta(), p_Particles[i].phi(), p_wt);
    }
    //TEST
    for (int i = 0; i < g_Particles.size(); ++ i) {
      ge_PtEtaPhi_tracks->Fill(g_Particles[i].pt(), g_Particles[i].eta(), g_Particles[i].phi(), g_wt);
    }
    

    
    vector<PseudoJet> p_cut_Particles = spart(p_Particles); vector<PseudoJet> g_cut_Particles = spart(g_Particles); //applying constituent cuts
    vector<PseudoJet> ch_p_cut_Particles = spart(ch_p_Particles); vector<PseudoJet> ch_g_cut_Particles = spart(ch_g_Particles);
    
    ClusterSequence p_Cluster(p_cut_Particles, jet_def); ClusterSequence g_Cluster(g_cut_Particles, jet_def);           //  CLUSTER BOTH
    ClusterSequence ch_p_Cluster(ch_p_cut_Particles, jet_def); ClusterSequence ch_g_Cluster(ch_g_cut_Particles, jet_def);
    p_Jets = sorted_by_pt(sjet(p_Cluster.inclusive_jets())); g_Jets = sorted_by_pt(sjet(g_Cluster.inclusive_jets()));    // EXTRACT JETS
    ch_p_Jets = sorted_by_pt(sjet(p_Cluster.inclusive_jets())); ch_g_Jets = sorted_by_pt(sjet(ch_g_Cluster.inclusive_jets()));

    //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
    for (int i = 0; i < p_Jets.size(); ++ i) {
      p_GroomedJets.push_back(sd(p_Jets[i]));
    }
    for (int i = 0; i < ch_p_Jets.size(); ++ i) {
      ch_p_GroomedJets.push_back(sd(ch_p_Jets[i]));
    }
    for (int i = 0; i < g_Jets.size(); ++ i) {
      g_GroomedJets.push_back(sd(g_Jets[i]));
    }
    for (int i = 0; i < ch_g_Jets.size(); ++ i) {
      ch_g_GroomedJets.push_back(sd(ch_g_Jets[i]));
    } 
    
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
	//TEMP                                                              
	int numch_py = 0;
	for (int k = 0; k < p_Jets[j].constituents().size(); ++ k) {
	  if (p_Jets[j].constituents()[k].user_index() != 0) {
	    numch_py ++;
	  }
	}
	ch_frac_v_pt_py->Fill(numch_py/(double)p_Jets[j].constituents().size(),p_Jets[j].pt(), p_wt);
	
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
      
	//TEMP                                                              
	int numch_ge = 0;
	for (int k = 0; k < g_Jets[j].constituents().size(); ++ k) {
	  if (g_Jets[j].constituents()[k].user_index() != 0) {
	    numch_ge ++;
	  }
	}
	ch_frac_v_pt_ge->Fill(numch_ge/(double) g_Jets[j].constituents().size(),g_Jets[j].pt(), g_wt);
	
      }
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HISTOS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    if(p_Jets.size() != 0) {
      FillHists(hists1D, hists2D, hists3D, "full", "_py", p_Jets, p_wt);
      for (int j = 0; j < p_Jets.size(); ++ j) {
	double val_list[nDim] = {p_GroomedJets[j].structure_of<contrib::SoftDrop>().symmetry(),p_GroomedJets[j].m(),p_GroomedJets[j].structure_of<contrib::SoftDrop>().delta_R(),p_GroomedJets[j].pt(), p_Jets[j].pt()}; //Groomed jet not guaranteed to be highest pT even though ungroomed one is
	SDnD_py->Fill(val_list);
      }
    }
    if(g_Jets.size() != 0) {
      FillHists(hists1D, hists2D, hists3D, "full", "_ge", g_Jets, g_wt);
      for (int j = 0; j < g_Jets.size(); ++ j) {
	double val_list[nDim] = {g_GroomedJets[j].structure_of<contrib::SoftDrop>().symmetry(),g_GroomedJets[j].m(),g_GroomedJets[j].structure_of<contrib::SoftDrop>().delta_R(),g_GroomedJets[j].pt(), g_Jets[j].pt()}; //Groomed jet not guaranteed to be highest pT even though ungroomed one is
	SDnD_ge->Fill(val_list);
      }
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

  std::cout << std::endl << std::endl << "Of " << nEventsPythia << " Pythia events, and " << nEventsGeant << " Pythia+Geant events" << std::endl;
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
  ch_matchedTree->Write("ch_py_ge_matchedTree");
  matchedTree->Write("py_ge_matchedTree");
  
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

  SDnD_py->Write(); SDnD_ge->Write();

  py_dPhi_trig_rec->Write(); ge_dPhi_trig_rec->Write();
  py_PtEtaPhi_tracks->Write(); ge_PtEtaPhi_tracks->Write();
  ch_frac_v_pt_py->Write(); ch_frac_v_pt_ge->Write();
  tow_id_v_e->Write();
  
  pt_response->Write();
  fout->Close();

  return 0;
}
