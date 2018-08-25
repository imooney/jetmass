//Isaac Mooney June, 2018 for jet mass analysis in simulation

#include "params.hh"
#include "funcs.hh"
#include "TStarJetPicoDefinitions.h"

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
    bool full = 1;  //full = 1 => ch+ne
    bool ge_or_py = 1; //ge_or_py = 1 => ge. ge_or_py = 0 => py
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
            if (arguments[3] == "ge") {ge_or_py = 1;} else if (arguments[3] == "py") {ge_or_py = 0;} else {cerr << "Not a valid flag!" << endl; exit(1);}
            chainList         = arguments[4];
            
            cout << outputDir << " " << outFileName << " " << full << " " << " " << ge_or_py << " " << chainList << endl;
            break;
        }
        default: { // Error: invalid custom settings
            __ERR("Invalid number of command line arguments");
            return -1;
            break;
        }
    }
    
    //  Initialize readers and provide chains
    TStarJetPicoReader Reader; TChain* Chain;
    if (ge_or_py == 0) {//pythia
        Chain = new TChain( "JetTreeMc" );    //  PURE PYTHIA DATA  (particle)
    }
    else {
        Chain = new TChain( "JetTree" );     //  CORRESPONDING GEANT DATA  (detector)
    }
    
    // Check to see if the input is a .root file or a .txt
    bool inputIsRoot = Analysis::HasEnding( chainList.c_str(), ".root" );
    bool inputIsTxt  = Analysis::HasEnding( chainList.c_str(), ".txt"  );
    bool inputIsList = Analysis::HasEnding( chainList.c_str(), ".list" );
    
    // If its a recognized file type, build the chain
    // If its not recognized, exit
    if ( inputIsRoot ) { Chain->Add( chainList.c_str()); }
    else if ( inputIsTxt )  { Chain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str());}
    else if ( inputIsList)  { Chain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str());}
    else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }
    
    //initialize the readers!
    if (ge_or_py == 0) {//pythia
        InitReader(Reader, Chain, nEvents, truth_triggerString, truth_absMaxVz, truth_vZDiff, truth_evPtMax, truth_evEtMax, truth_evEtMin, truth_DCA, truth_NFitPts, truth_FitOverMaxPts, sim_maxEtTow, sim_badTowers, sim_bad_run_list);
    }
    if (ge_or_py == 1) {//geant
        InitReader(Reader, Chain, nEvents, det_triggerString, det_absMaxVz, det_vZDiff, det_evPtMax, det_evEtMax, det_evEtMin, det_DCA, det_NFitPts, det_FitOverMaxPts, sim_maxEtTow, sim_badTowers, sim_bad_run_list);
    }
    TStarJetPicoEventHeader* header;  TStarJetPicoEvent* event;
    TStarJetVectorContainer<TStarJetVector> * container;        TStarJetVector* sv;
                                                                
    double n_jets, wt;
    vector<double> Pt; vector<double> Eta; vector<double> Phi; vector<double> M; vector<double> E;
    vector<double> ch_e_frac;
    vector<double> zg; vector<double> rg; vector<double> mg; vector<double> ptg;
    
    TTree *eventTree = new TTree("event","event");
    eventTree->Branch("n_jets", &n_jets);
    eventTree->Branch("Pt", &Pt); eventTree->Branch("Eta",&Eta); eventTree->Branch("Phi",&Phi); eventTree->Branch("M",&M); eventTree->Branch("E",&E);
    eventTree->Branch("ch_e_frac", &ch_e_frac);
    eventTree->Branch("zg", &zg); eventTree->Branch("rg", &rg); eventTree->Branch("mg", &mg); eventTree->Branch("ptg",&ptg);
    eventTree->Branch("weight", &wt);
    
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
    
    Selector select_jet_pt_min;
    if (ge_or_py == 0) {//pythia
      select_jet_pt_min  = fastjet::SelectorPtMin( jet_ptmin );
    }
    else { select_jet_pt_min = fastjet::SelectorPtMin( det_jet_ptmin ); }
    Selector select_jet_pt_max  = fastjet::SelectorPtMax( jet_ptmax );
    Selector sjet = select_jet_rap && select_jet_pt_min && select_jet_pt_max;
    
    JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
    TString Filename;
    
    // Particle containers & counters
    vector<PseudoJet> Particles, Jets, GroomedJets;
    int nEvents = 0;   int NJets = 0;  int EventID;
    //1=inclusive, 2=lead
    int counter_debug1 = 0, counter_debug2 = 0;

    // LOOP!
    
    while ( Reader.NextEvent() ) {
      
      //clearing vectors                                                                                                                                                                              
      Pt.clear(); Eta.clear(); Phi.clear(); M.clear(); E.clear();
      zg.clear(); rg.clear(); mg.clear(); ptg.clear();
      ch_e_frac.clear();
      //initializing variables to -9999                                                                                                                                                               
      n_jets = -9999; wt = -9999;
      
      EventID = Reader.GetNOfCurrentEvent();
      
      event = Reader.GetEvent();    header = event->GetHeader();

      EventID = Reader.GetNOfCurrentEvent();
      
      Particles.clear();
      Jets.clear();
      GroomedJets.clear();
      
      nEvents ++; Reader.PrintStatus(10);     // Print out reader status every 10 seconds
      
      event = Reader.GetEvent();    header = event->GetHeader();           // Get  event header and event
      
      container = Reader.GetOutputContainer();      //  container
      
      Filename =  Reader.GetInputChain()->GetCurrentFile()->GetName();
      
      wt = LookupRun12Xsec( Filename );

      //  GATHER PARTICLES
      GatherParticles ( container, sv, Particles, 1, ge_or_py);    // first bool flag: 0 signifies charged-only, 1 = ch+ne;  particles: second bool flag: pythia = 1,  = 0
      
      vector<PseudoJet> cut_Particles = spart(Particles); //applying constituent cuts
      
      ClusterSequence Cluster(cut_Particles, jet_def);           //  CLUSTER BOTH
      vector<PseudoJet> JetsInitial;
      
      if (ge_or_py == 0) {//pythia
	JetsInitial = sorted_by_pt(sjet(Cluster.inclusive_jets()));    // EXTRACT JETS
      }
      else {
	JetsInitial = sorted_by_pt(sjet(Cluster.inclusive_jets())); //apply jet cuts in Geant (not Pythia!)
      }
      
      vector<PseudoJet> Jets;
      
      //Implementing a neutral energy fraction cut of 90% on inclusive jets                                                                                                                           
      if (ge_or_py == 1) { //geant
	ApplyNEFSelection(JetsInitial, Jets);
      }
      else {Jets = JetsInitial;}
      
      //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
      for (int i = 0; i < Jets.size(); ++ i) {
	GroomedJets.push_back(sd(Jets[i]));
      }
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TREES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      
      NJets += Jets.size();               //  Save jet info and add jets to total
      n_jets = Jets.size();
      
      for (int i = 0; i < n_jets; ++ i) {
	if (ge_or_py == 0) { //pythia
	  //	  if (Jets[i].pt() > 24 && Jets[i].pt() < 26) {cout << Filename << " " << wt << " " << Jets[i].pt() << " " << Jets[i].eta() << " " << Jets[i].phi() << " " << Jets[i].m() << endl;}
	  std::string tail = ((string) Filename).substr(((string) Filename).size() - 10);
	  std::string upstring = tail.substr(0,2);
	  if (upstring.find("_") != std::string::npos || upstring.find("-") != std::string::npos) { if (upstring.substr(1,1) != "_") {upstring = upstring.substr(1,1);} else {upstring = upstring.substr(0,1);}}
	  double upbin = std::stoi(upstring);
	  if (Jets[i].pt() > 2*upbin) { std::cout << Filename << " " << EventID << " " << wt << " " << Jets[i].pt() << " " << Jets[i].eta() << " " << Jets[i].phi() << " " << Jets[i].m() << std::endl; }
	}
	Pt.push_back(Jets[i].pt()); Eta.push_back(Jets[i].eta()); Phi.push_back(Jets[i].phi());
	M.push_back(Jets[i].m()); E.push_back(Jets[i].e());
	zg.push_back(GroomedJets[i].structure_of<SD>().symmetry()); rg.push_back(GroomedJets[i].structure_of<SD>().delta_R());
	mg.push_back(GroomedJets[i].m()); ptg.push_back(GroomedJets[i].pt());
	double ch_e = 0; double tot_e = 0;
	vector<PseudoJet> cons = Jets[i].constituents();
	for (int j = 0; j < cons.size(); ++ j) {
	  if (cons[j].user_index() != 0) {ch_e += cons[j].e();}
	  tot_e += cons[j].e();
	}
	ch_e_frac.push_back(ch_e/(double)tot_e);
      }
      if (Jets.size() != 0) {
	eventTree->Fill();
      }
    }
    
    //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
    
    TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
    
    std::cout << std::endl << std::endl << "Of " << nEvents << " events" << std::endl;
    std::cout << NJets << " jets have been found" << std::endl;
    std::cout <<std::endl << "Writing to:  " << fout->GetName() << std::endl << std::endl;
   
    eventTree->Write();
    
    fout->Close();
    
    return 0;
    
}
