//Isaac Mooney June, 2018 for jet mass analysis in simulation

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
    
    // Histograms
    
    Collection<string, TH1D> hists1D; Collection<string, TH2D> hists2D; Collection<string, TH3D> hists3D;
    
    vector<string> flag_j = {"lead", "sublead", "trig", "rec", "incl"};
    vector<string> flag_k = {"jet", "cons"};
    for (int j = 0; j < flag_j.size(); ++ j) {
        for (int k = 0; k < flag_k.size(); ++ k) {
                hists1D.add(("m_"+flag_j[j]+"_"+flag_k[k]).c_str(),"",20,0,10); //mass
                hists2D.add(("m_v_pt_"+flag_j[j]+"_"+flag_k[k]).c_str(),"",20,0,10,11,5,60); //mass vs. pT
                hists3D.add(("PtEtaPhi_"+flag_j[j]+"_"+flag_k[k]).c_str(),"",11,5,60,30,-0.6,0.6,50,0,2*Pi);
            }
        }
    
    const unsigned nDim = 5;
    int bins[nDim] = {20, 20, 20, 11, 11};
    double min[nDim] = {0,0,0,5,5};
    double max[nDim] = {1,10,1,60,60};
    THnSparse * SDnD = new THnSparseD("zg_mg_thetag_ptg_pt_incl_sd", "", nDim, bins, min, max);
    SDnD->Sumw2();
    
    //TESTS!
    TH1D * dPhi_trig_rec = new TH1D("dPhi_trig_rec",";#Delta #phi;arb.", 28, -Pi - 0.4, Pi + 0.4); //defined as trigger - recoil
    TH3D * PtEtaPhi_tracks = new TH3D("PtEtaPhi_tracks",";p^{track}_{T} [GeV/c]; #eta; #phi", 80, 0, 80, 30, -0.6,0.6,50,0,2*Pi);
    TH2D * ch_frac_v_pt = new TH2D("ch_frac_v_pt_ge",";charged fraction; p_{T}^{jet} [GeV/c]",10,0,1,11,5,60);
    TH2D * tow_id_v_e = new TH2D("tow_id_v_e",";tower ID; tower E_{T} [GeV]",4800,1,4801,140,0,140);
    
    double jetPt_incl, jetEta_incl, jetPhi_incl, jetM_incl, jetE_incl, consPt_incl, consEta_incl, consPhi_incl, consM_incl, consE_incl, wt;
    int nCons_incl, cons_dummy;
    double jetPt_lead, jetEta_lead, jetPhi_lead, jetM_lead, jetE_lead, consPt_lead, consEta_lead, consPhi_lead, consM_lead, consE_lead;
    double jetPt_sublead, jetEta_sublead, jetPhi_sublead, jetM_sublead, jetE_sublead;
    int nCons_lead, nCons_sublead;
    
    TTree *inclTree = new TTree("inclTree","inclTree");
    TTree *leadTree = new TTree("leadTree","leadTree");
    TTree *subleadTree = new TTree("subleadTree","subleadTree");
    TTree *cons_inclTree = new TTree("cons_inclTree","cons_inclTree");
    TTree *cons_leadTree = new TTree("cons_leadTree","cons_leadTree");
    
    leadTree->Branch("Pt", &jetPt_lead); leadTree->Branch("Eta", &jetEta_lead); leadTree->Branch("Phi", &jetPhi_lead);
    leadTree->Branch("M", &jetM_lead); leadTree->Branch("E", &jetE_lead); leadTree->Branch("nCons", &nCons_lead);
    leadTree->Branch("weight", &wt);
    
    subleadTree->Branch("Pt", &jetPt_sublead); subleadTree->Branch("Eta", &jetEta_sublead); subleadTree->Branch("Phi", &jetPhi_sublead);
    subleadTree->Branch("M", &jetM_sublead); subleadTree->Branch("E", &jetE_sublead); subleadTree->Branch("nCons", &nCons_sublead);
    subleadTree->Branch("weight", &wt);
    
    cons_leadTree->Branch("Pt", &consPt_lead); cons_leadTree->Branch("Eta", &consEta_lead); cons_leadTree->Branch("Phi", &consPhi_lead);
    cons_leadTree->Branch("M", &consM_lead); cons_leadTree->Branch("E", &consE_lead);
    cons_leadTree->Branch("weight", &wt);
    
    inclTree->Branch("Pt", &jetPt_incl); inclTree->Branch("Eta", &jetEta_incl); inclTree->Branch("Phi", &jetPhi_incl);
    inclTree->Branch("M", &jetM_incl); inclTree->Branch("E", &jetE_incl); inclTree->Branch("nCons", &nCons_incl);
    inclTree->Branch("weight", &wt);
    
    cons_inclTree->Branch("Pt", &consPt_incl); cons_inclTree->Branch("Eta", &consEta_incl); cons_inclTree->Branch("Phi", &consPhi_incl);
    cons_inclTree->Branch("M", &consM_incl); cons_inclTree->Branch("E", &consE_incl);
    cons_inclTree->Branch("weight", &wt);
    
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
    //double wt;
    TString Filename;
    
    // Particle containers & counters
    vector<PseudoJet> Particles, Jets, GroomedJets;
    int nEvents = 0;   int NJets = 0;  int EventID;
    //1=inclusive, 2=lead
    int counter_debug1 = 0, counter_debug2 = 0;

// LOOP!

while ( Reader.NextEvent() ) {
    EventID = Reader.GetNOfCurrentEvent();
    
    event = Reader.GetEvent();    header = event->GetHeader();
    
    //TEMP!!!!!!!!!!!!!!!!!
    /*
    TStarJetPicoEvent * current_event = Reader.GetEvent();
    //      cout << "event " << reader.GetNOfCurrentEvent() << " has these tows" << endl;
    for (int i = 0; i < current_event->GetTowers()->GetEntries(); ++ i) {
        tow_id_v_e->Fill(current_event->GetTower(i)->GetId(), current_event->GetTower(i)->GetEnergy());
        //cout << current_event->GetTower(i)->GetId() << " " << current_event->GetTower(i)->GetEnergy() << endl;
    }
     */
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
    
    //TEST
    for (int i = 0; i < Particles.size(); ++ i) {
      PtEtaPhi_tracks->Fill(Particles[i].pt(), Particles[i].eta(), Particles[i].phi(), wt);
    }
    
    vector<PseudoJet> cut_Particles = spart(Particles); //applying constituent cuts
    
    ClusterSequence Cluster(cut_Particles, jet_def);           //  CLUSTER BOTH
    vector<PseudoJet> JetsInitial;
    
    if (ge_or_py == 0) {//pythia
      JetsInitial = sorted_by_pt(Cluster.inclusive_jets());    // EXTRACT JETS
    }
    else {
      JetsInitial = sorted_by_pt(sjet(Cluster.inclusive_jets())); //apply jet cuts in Geant (not Pythia!)
    }
    
    vector<PseudoJet> Jets;

    //Implementing a neutral energy fraction cut of 90% on inclusive jets                                                                                                                           
    ApplyNEFSelection(JetsInitial, Jets);
    /*
    for (int i = 0; i < JetsInitial.size(); ++ i) {
      double towersum = 0; double ptsum = 0;
      for (int j = 0; j < JetsInitial[i].constituents().size(); ++ j) {
	if (JetsInitial[i].constituents()[j].user_index() == 0) {
	  towersum += JetsInitial[i].constituents()[j].pt();
	}
	ptsum += JetsInitial[i].constituents()[j].pt();
      }
      if (towersum / (double) ptsum < NEF_max) {
	Jets.push_back(JetsInitial[i]);
      }
    }
    */
    
    
    //loop over the jets which passed cuts, groom them, and add to a vector (sorted by pt of the original jet)
    for (int i = 0; i < Jets.size(); ++ i) {
      GroomedJets.push_back(sd(Jets[i]));
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    if (Jets.size() != 0) {
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      //leading
      vector<PseudoJet> Lead; Lead.push_back(Jets[0]);
      FillTrees(Lead, leadTree, jetPt_lead, jetEta_lead, jetPhi_lead, jetM_lead, jetE_lead, nCons_lead, wt, wt);
      //subleading
      if (Jets.size() > 1) {
	vector<PseudoJet> Sublead; Sublead.push_back(Jets[1]);
	FillTrees(Sublead, subleadTree, jetPt_sublead, jetEta_sublead, jetPhi_sublead, jetM_sublead, jetE_sublead, nCons_sublead, wt, wt);
      }
      //constituents
      FillTrees(Lead[0].constituents(), cons_leadTree, consPt_lead, consEta_lead, consPhi_lead, consM_lead, consE_lead, cons_dummy, wt, wt);
      //inclusive
      FillTrees(Jets, inclTree, jetPt_incl, jetEta_incl, jetPhi_incl, jetM_incl, jetE_incl, nCons_incl, wt, wt);
      //constituents
      for(int j = 0; j < Jets.size(); ++ j) {
	FillTrees(Jets[j].constituents(), cons_inclTree, consPt_incl, consEta_incl, consPhi_incl, consM_incl, consE_incl, cons_dummy, wt, wt);
	//TEMP
	int numch = 0;
	for (int k = 0; k < Jets[j].constituents().size(); ++ k) {
	  if (Jets[j].constituents()[k].user_index() != 0) {
	    numch ++;
	  }
	}
	ch_frac_v_pt->Fill(numch/(double) Jets[j].constituents().size(),Jets[j].pt(), wt);
        
      }
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HISTOS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    if(Jets.size() != 0) {
      FillHists(hists1D, hists2D, hists3D, Jets, wt);
      for (int j = 0; j < Jets.size(); ++ j) {
	double val_list[nDim] = {GroomedJets[j].structure_of<contrib::SoftDrop>().symmetry(),GroomedJets[j].m(),GroomedJets[j].structure_of<contrib::SoftDrop>().delta_R(),GroomedJets[j].pt(), Jets[j].pt()}; //Groomed jet not guaranteed to be highest pT even though ungroomed one is
	SDnD->Fill(val_list);
      }
    }
    NJets += Jets.size();               //  Save jet info and add jets to total
 }
 
//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
 
 
 TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");
 
 std::cout << std::endl << std::endl << "Of " << nEvents << " events" << std::endl;
 std::cout << NJets << " jets have been found" << std::endl;
 std::cout <<std::endl << "Writing to:  " << fout->GetName() << std::endl << std::endl;
 
 leadTree->Write("leadTree");
 subleadTree->Write("subleadTree");
 cons_leadTree->Write("cons_leadTree");
 inclTree->Write("inclTree");
 cons_inclTree->Write("cons_inclTree");
 
 //hist1->Write(); hist2->Write(); etc, goes here
 for (int j = 0; j < flag_j.size(); ++ j) {
   for (int k = 0; k < flag_k.size(); ++ k) {
     hists1D.write(("m_"+flag_j[j]+"_"+flag_k[k]).c_str());
     hists2D.write(("m_v_pt_"+flag_j[j]+"_"+flag_k[k]).c_str());
     hists3D.write(("PtEtaPhi_"+flag_j[j]+"_"+flag_k[k]).c_str());
   }
 }
 
 SDnD->Write();
 
 dPhi_trig_rec->Write();
 PtEtaPhi_tracks->Write();
 ch_frac_v_pt->Write();
 tow_id_v_e->Write();
 
 fout->Close();
 
 return 0;
 
}
