//! macro to read in PY8 or HW7 files and create the recursive histograms from the trees that can be read by the plotter. 


#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TProfile.h>
#include <TLine.h>

#include <TROOT.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TChain.h>
#include <TBranch.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>

//#include "TStarJetVector.h"
//#include "TStarJetVectorJet.h"

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <exception>

const int npthatbins = 11;
int pthatbins_HW7[npthatbins+1]        = {3, 4, 5, 7,  9, 11, 15, 25, 35, 45, 55, 65};
int pthatbins_HW7_ranges[npthatbins+1] = {4, 5, 6, 8, 10, 12, 17, 27, 38, 48, 58, 80};
int pthatbins_PY8[npthatbins+1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 80};


// double pthatweight_HW7[npthatbins+1] = {3.67764e+07, 4.39797e+06, 707461, 14372.6, 0, 77.2581, 3.95503, 0.0375627, 0.00117562, 4.87868e-05, 1.9473e-06};
double pthatweight_HW7[npthatbins+1] = {3.67764e+07, 4.94731e+06, 866978, 26833.9, 0, 119.408, 4.35634, 0.0444399, 0.00156817, 7.82549e-05, 4.3366e-06}; 
// double pthatweight_PY8[npthatbins+1] = {118402, 184.206, 3.35249, 0.235818, 0.029338, 0.00475577, 0.000888841, 0.000177773, 3.66509e-05, 9.05136e-06, 9.21373e-07};
double pthatweight_PY8[npthatbins+1] = {1240.43, 9.18156, 0.615729, 0.0779501, 0.0132134, 0.00261512, 0.000563591, 0.000126293, 2.8601e-05, 7.7384e-06, 9.04096e-07};//{446225, 184.244, 3.35316, 0.23586, 0.0293411, 0.00475607, 0.000888898, 0.000177778, 3.6652e-05, 9.05141e-06, 9.21373e-07};


void combineMCpthat(std::string mctype = "pythia8", int radius = 4)
{

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  bool printDebug = false;
  bool doHC = true;
  bool doMC = true;
  
  // vector<int> pthatbins = 0;

  int *pthatbins=0;
  if(mctype == "pythia8")
    pthatbins = pthatbins_PY8;
  if(mctype == "herwig7")
    pthatbins = pthatbins_HW7;


  double *pthatweight = 0;
  if(mctype == "pythia8")
    pthatweight = pthatweight_PY8;
  if(mctype == "herwig7")
    pthatweight = pthatweight_HW7;
  
  
  // if(mctype == "herwig7"){
  //   pthatbins->push_back(3);
  //   pthatbins->push_back(4);
  //   pthatbins->push_back(5);
  //   pthatbins->push_back(7);
  //   pthatbins->push_back(9);
  //   pthatbins->push_back(11);
  //   pthatbins->push_back(15);
  //   pthatbins->push_back(25);
  //   pthatbins->push_back(35);
  //   pthatbins->push_back(45);
  //   pthatbins->push_back(55);
  //   pthatbins->push_back(65);
  // }
  // if(mctype == "pythia8"){
  //   pthatbins->push_back(5);
  //   pthatbins->push_back(10);
  //   pthatbins->push_back(15);
  //   pthatbins->push_back(20);
  //   pthatbins->push_back(25);
  //   pthatbins->push_back(30);
  //   pthatbins->push_back(35);
  //   pthatbins->push_back(40);
  //   pthatbins->push_back(45);
  //   pthatbins->push_back(50);
  //   pthatbins->push_back(60);
  //   pthatbins->push_back(80);
  // }

  // std::cout<<"pthatbin = "<<pthatbins->size()<<std::endl;

  //! load in the files.
  //! loop over the pthat bins and fill in the histograms 
  //! bin boundaries
  // Set up response matrix
  // ----------------------
  int nZgBinsTrue = 10;
  float zgminTrue = 0.00;
  float zgmaxTrue = 0.5;
  // int nZgBinsTrue = 20;
  // float zgminTrue = 0.05;
  // float zgmaxTrue = 0.55;
  int nZgBinsMeas = nZgBinsTrue;
  float zgminMeas = zgminTrue;
  float zgmaxMeas = zgmaxTrue;

  int nRgBinsTrue = 20;
  float rgminTrue = 0.0;
  float rgmaxTrue = 1.0;
  int nRgBinsMeas = 20;
  float rgminMeas = 0.0;
  float rgmaxMeas = 1.0;

  int nMBinsTrue = 40;
  float MminTrue = 0.0;
  float MmaxTrue = 16;
  int nMBinsMeas = 40;
  float MminMeas = 0.0;
  float MmaxMeas = 16;

  bool addpTCut = true;
  
  int nPtBins = 80;
  float ptmin=0;
  float ptmax=80;
  int nPtBinsTrue =  nPtBins;
  float ptminTrue =  ptmin;
  float ptmaxTrue =  ptmax;
  int nPtBinsMeas =  60;
  float ptminMeas =  0;
  float ptmaxMeas =  60;

  const int nSplits = 4;
  
  if(addpTCut){
    nPtBins=15;
    ptmin=5;
    ptmax=80;
    nPtBinsTrue =  nPtBins;
    ptminTrue =  ptmin;
    ptmaxTrue =  ptmax;
    nPtBinsMeas =  9;
    ptminMeas =  15;
    ptmaxMeas =  60;
  }


  //! Declare ptbins for publication 
  double pubptbins[] = {15.0, 20.0, 25.0, 30.0, 40.0, 60.0};
  const int npubptbins = sizeof(pubptbins)/sizeof(double)-1;

  //! declare the output file and histograms
  TFile * fout = new TFile(Form("Results/%s_recSD_histograms_R0%d_decays_off.root", mctype.c_str(), radius),"RECREATE");
  fout->cd();

  double groompTbins[] = {15.0, 20.0, 25.0, 30.0, 40.0, 60.0};
  const int ngroompTbins = sizeof(groompTbins)/sizeof(double)-1;
  
  double pt_testbins[] = {0.2,0.7,1.2,1.7,2.2,3.2,5.2,8.2,12.2,17.2,23.2,30.2,38.2,50.2,80.2};
  const int npt_testbins = sizeof(pt_testbins)/sizeof(double)-1;
  
  TH2D * hrecpTvtrigpT = new TH2D("hrecpTvtrigpT",";p^{rec. jet}_{T} [GeV/c];p^{trig.}_{T} [GeV/c]",npt_testbins,pt_testbins,21,9,30);
  
  
  TH2D * hGenZg_Split_InitpTbins[ngroompTbins];
  TH2D * hGenRg_Split_InitpTbins[ngroompTbins];

  for(int i = 0; i<ngroompTbins; ++i){

    hGenZg_Split_InitpTbins[i] = new TH2D(Form("hGenZg_Split_InitpTbins_%d", i), "", 8, 0, 8, nZgBinsTrue, zgminTrue, zgmaxTrue);
    hGenRg_Split_InitpTbins[i] = new TH2D(Form("hGenRg_Split_InitpTbins_%d", i), "", 8, 0, 8, nRgBinsTrue, rgminTrue, rgmaxTrue);

  }

  
  TH2D * hGen_Zg_vs_Split[npubptbins];
  TH2D * hGen_Rg_vs_Split[npubptbins];
  TH1D * hGen_Zg[npubptbins];
  TH1D * hGen_Rg[npubptbins];
  TH1D * hGen_M[npubptbins];
  TH1D * hGenJetpT = new TH1D("hGenJetpT","", 100, 0, 100);
  TH1D * hGenJeteta = new TH1D("hGenJeteta","", 50, -1, 1);
  TH1D * hGenJetphi = new TH1D("hGenJetphi","", 50, 0, 6.5);

  
  TH1D * hHCGen_Zg[npubptbins];
  TH1D * hHCGen_Rg[npubptbins];
  TH1D * hHCGen_M[npubptbins];
  TH1D * hHCGenJetpT = new TH1D("hHCGenJetpT","", 100, 0, 100);
  TH1D * hHCGenJetpT_wWeight = new TH1D("hHCGenJetpT_wWeight","", 100, 0, 100);
  TH1D * hHCGenJeteta = new TH1D("hHCGenJeteta","", 50, -1, 1);
  TH1D * hHCGenJetphi = new TH1D("hHCGenJetphi","", 50, 0, 6.5);

  TH1D * hMCGen_Zg[npubptbins];
  TH1D * hMCGen_Rg[npubptbins];
  TH1D * hMCGen_M[npubptbins];
  TH1D * hMCGenJetpT = new TH1D("hMCGenJetpT","", 100, 0, 100);
  TH1D * hMCGenJetpT_wWeight = new TH1D("hMCGenJetpT_wWeight","", 100, 0, 100);
  TH1D * hMCGenJeteta = new TH1D("hMCGenJeteta","", 50, -1, 1);
  TH1D * hMCGenJetphi = new TH1D("hMCGenJetphi","", 50, 0, 6.5);

  
  TH1D * hGenJetpthat = new TH1D("hGenJetpthat","", 100, 0, 100);
  TH1D * hGenJetpT_wWeight = new TH1D("hGenJetpT_wWeight","", 100, 0, 100);
  TH1D * hGenJetpthat_wWeight = new TH1D("hGenJetpthat_wWeight","", 100, 0, 100);

  for(int i = 0; i<npubptbins; ++i){
    hGen_Zg_vs_Split[i] = new TH2D(Form("hGen_Zg_vs_Split_%d", i), "", 8, 0, 8, nZgBinsTrue, zgminTrue, zgmaxTrue);
    hGen_Rg_vs_Split[i] = new TH2D(Form("hGen_Rg_vs_Split_%d", i), "", 8, 0, 8, nRgBinsTrue, rgminTrue, rgmaxTrue);
    hGen_Zg[i] = new TH1D(Form("hGen_Zg_ptbin_%d", i), "", nZgBinsMeas, zgminMeas, zgmaxMeas);
    hGen_Rg[i] = new TH1D(Form("hGen_Rg_ptbin_%d", i), "", nRgBinsMeas, rgminMeas, rgmaxMeas);
    hGen_M[i] = new TH1D(Form("hGen_M_ptbin_%d", i), "", nMBinsMeas, MminMeas, MmaxMeas);

    if(doHC){      
      hHCGen_Zg[i] = new TH1D(Form("hHCGen_Zg_ptbin_%d", i), "", nZgBinsMeas, zgminMeas, zgmaxMeas);
      hHCGen_Rg[i] = new TH1D(Form("hHCGen_Rg_ptbin_%d", i), "", nRgBinsMeas, rgminMeas, rgmaxMeas);
      hHCGen_M[i] = new TH1D(Form("hHCGen_M_ptbin_%d", i), "", nMBinsMeas, MminMeas, MmaxMeas);
    }

    if(doMC){      
      hMCGen_Zg[i] = new TH1D(Form("hMCGen_Zg_ptbin_%d", i), "", nZgBinsMeas, zgminMeas, zgmaxMeas);
      hMCGen_Rg[i] = new TH1D(Form("hMCGen_Rg_ptbin_%d", i), "", nRgBinsMeas, rgminMeas, rgmaxMeas);
      hMCGen_M[i] = new TH1D(Form("hMCGen_M_ptbin_%d", i), "", nMBinsMeas, MminMeas, MmaxMeas);
    }

    
  }


  double xsec_avg[npthatbins] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  long nevts_pthatbin[npthatbins] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  
  //! File name
  // std::string fileloc = "/wsu/home/gp/gp15/gp1574/WORK/MC/ANALYSIS/MCrootfilesr_recSD_norecGroom";
  std::string fileloc = "~/jetmass/production";
  
  int tot_trigs = 0;
  
  for(unsigned inpt = 0; inpt < npthatbins; ++inpt){

    if(printDebug)
      std::cout<<"inpt = "<<inpt<<", pthatbins = ("<<pthatbins[inpt]<<", "<<pthatbins[inpt+1]<<")"<<std::endl;
    
    std::string pthatbinname = Form("star_subjetvars_%s_%dpthatbin%d_R0%d_decays_off.root",
				    mctype.c_str(),
				    pthatbins[inpt],
				    pthatbins[inpt+1],
				    radius
				    );

    // if(printDebug)
    std::cout<<"running on "<<pthatbinname<<std::endl;

    TFile* fin = TFile::Open(Form("%s/%s", fileloc.c_str(), pthatbinname.c_str()));

    if(printDebug)
      fin->ls();
    
    TTree * trigrecTree = (TTree*)fin->Get("TrigRecTree");
    TTree * jetTree = (TTree*)fin->Get("ResultTree");
    TTree * HCjetTree;
    if(doHC)
      HCjetTree = (TTree*)fin->Get("HCResultTree");
    TTree * MCjetTree;
    if(doMC)
      MCjetTree = (TTree*)fin->Get("MCResultTree");
    //! Set branch address for the tree 

    double cent;
    double mcweight;
    double pthat;

    double HCcent;
    double HCmcweight;
    double HCpthat;

    double MCcent;
    double MCmcweight;
    double MCpthat;

    double mcweight_test;
    double pthat_test;
    
    vector<double> *recpT = 0;
    vector<double> *trigpT = 0;
    trigrecTree->SetBranchAddress("recpT",&recpT);
    trigrecTree->SetBranchAddress("trigpT",&trigpT);
    trigrecTree->SetBranchAddress("pthat",&pthat_test);
    trigrecTree->SetBranchAddress("mcweight",&mcweight_test);
    
    vector<double> *jetpT = 0;
    vector<double> *jetM = 0;
    vector<double> *sdjetM = 0;
    vector<double> *jeteta = 0;
    vector<double> *jetphi = 0;
    vector<double> *zg = 0;
    vector<double> *rg = 0;
    vector<vector<double> > *rec_zg = 0;
    vector<vector<double> > *rec_rg = 0;
    vector<vector<double> > *sdjetpT = 0;
    jetTree->SetBranchAddress("cent", &cent);
    jetTree->SetBranchAddress("mcweight", &mcweight);
    jetTree->SetBranchAddress("pthat", &pthat);
    jetTree->SetBranchAddress("jetpT", &jetpT);
    jetTree->SetBranchAddress("jetM", &jetM);
    jetTree->SetBranchAddress("rec_sdpt", &sdjetpT);
    jetTree->SetBranchAddress("sdjetM", &sdjetM);
    jetTree->SetBranchAddress("jeteta", &jeteta);
    jetTree->SetBranchAddress("jetphi", &jetphi);
    jetTree->SetBranchAddress("zg", &zg);
    jetTree->SetBranchAddress("rg", &rg);
    jetTree->SetBranchAddress("rec_zg", &rec_zg);
    jetTree->SetBranchAddress("rec_rg", &rec_rg);
    
    vector<double> *HCjetpT = 0;
    vector<double> *HCjetM = 0;
    vector<double> *HCsdjetM = 0;
    vector<double> *HCsdjetpT = 0;
    vector<double> *HCjeteta = 0;
    vector<double> *HCjetphi = 0;
    vector<double> *HCzg = 0;
    vector<double> *HCrg = 0;    
    if(doHC){
      HCjetTree->SetBranchAddress("cent", &HCcent);
      HCjetTree->SetBranchAddress("mcweight", &HCmcweight);
      HCjetTree->SetBranchAddress("pthat", &HCpthat);
      HCjetTree->SetBranchAddress("jetpT", &HCjetpT);
      HCjetTree->SetBranchAddress("jetM", &HCjetM);
      HCjetTree->SetBranchAddress("sdjetpT", &HCsdjetpT);
      HCjetTree->SetBranchAddress("sdjetM", &HCsdjetM);
      HCjetTree->SetBranchAddress("jeteta", &HCjeteta);
      HCjetTree->SetBranchAddress("jetphi", &HCjetphi);
      HCjetTree->SetBranchAddress("zg", &HCzg);
      HCjetTree->SetBranchAddress("rg", &HCrg);
    }

    vector<double> *MCjetpT = 0;
    vector<double> *MCjetM = 0;
    vector<double> *MCsdjetM = 0;
    vector<double> *MCsdjetpT = 0;
    vector<double> *MCjeteta = 0;
    vector<double> *MCjetphi = 0;
    vector<double> *MCzg = 0;
    vector<double> *MCrg = 0;    
    if(doMC){
      MCjetTree->SetBranchAddress("cent", &MCcent);
      MCjetTree->SetBranchAddress("mcweight", &MCmcweight);
      MCjetTree->SetBranchAddress("pthat", &MCpthat);
      MCjetTree->SetBranchAddress("jetpT", &MCjetpT);
      MCjetTree->SetBranchAddress("jetM", &MCjetM);
      MCjetTree->SetBranchAddress("sdjetpT", &MCsdjetpT);
      MCjetTree->SetBranchAddress("sdjetM", &MCsdjetM);
      MCjetTree->SetBranchAddress("jeteta", &MCjeteta);
      MCjetTree->SetBranchAddress("jetphi", &MCjetphi);
      MCjetTree->SetBranchAddress("zg", &MCzg);
      MCjetTree->SetBranchAddress("rg", &MCrg);
    }
    
    double n_trig = trigrecTree->GetEntries();

    cout <<"DEBUG BEGIN" << endl;
    cout << trigrecTree->GetEntries() << endl;
    for(long Nev = 0; Nev<trigrecTree->GetEntries(); ++Nev) {
      trigrecTree->GetEntry(Nev);
      if(Nev%10000 == 0)
	std::cout<<"Event Number "<<Nev<<std::endl;
      
      if(printDebug){
	//std::cout<<"There are "<<recpT->size()<<" recoil jets in this event"<<std::endl;
      }

      if(mctype == "herwig7"){
	if(pthat < pthatbins[inpt])
	  continue;
      }else if(mctype == "pythia8"){
	if(pthat_test >= pthatbins[inpt+1] || pthat_test < pthatbins[inpt])
	  continue;
      }

      // double xsecweight = pthatweight[inpt];
      //double xsecweight = mcweight_test;

      if(printDebug)
	std::cout<<"Running through the event now"<<std::endl;

      //      xsec_avg[inpt] += mcweight_test;
      //nevts_pthatbin[inpt] += 1;

      for (int j = 0; j < recpT->size(); ++j) {
	if (recpT->at(j) > 0) {
	  hrecpTvtrigpT->Fill(recpT->at(j), trigpT->at(0), mcweight_test);
	}
      }
      
    }
    cout << "NUMBER OF TRIGGERS: " << n_trig << endl;
    tot_trigs += n_trig;
    
    //! loop over the events and fill the histograms
    for(long Nev = 0; Nev<jetTree->GetEntries(); ++Nev){
    
      jetTree->GetEntry(Nev);
      if(Nev%10000 == 0)
	std::cout<<"Event Number "<<Nev<<std::endl;
      
      if(printDebug){
	std::cout<<"There are "<<jetpT->size()<<" Jets in this event"<<std::endl;
      }

      if(mctype == "herwig7"){
	if(pthat < pthatbins[inpt])
	  continue;
      }else if(mctype == "pythia8"){
	if(pthat >= pthatbins[inpt+1] || pthat < pthatbins[inpt])
	  continue;
      }

      // double xsecweight = pthatweight[inpt];
      double xsecweight = mcweight;

      if(printDebug)
	std::cout<<"Running through the event now"<<std::endl;
      
      hGenJetpthat->Fill(pthat, 1);
      hGenJetpthat_wWeight->Fill(pthat, xsecweight);
      
      if(jetpT->size()==0) continue;

      xsec_avg[inpt] += mcweight;
      nevts_pthatbin[inpt] += 1;
      
      //! loop over the jets in that event 
      for(unsigned ij = 0; ij<jetpT->size(); ++ij){
      
	if(printDebug){
	  std::cout<<"    Jet pt = "<<jetpT->at(ij)<<std::endl;
	}

	if(printDebug)
	  std::cout<<"Inside the jet loop"<<std::endl;


	if(mctype == "herwig7")	  
	  if(jetpT->at(ij) > 3*pthatbins[inpt+1])
	    continue;
	  
	if(mctype == "pythia8")
	  if(jetpT->at(ij) > 3*pthatbins[inpt+1])
	    continue;

	if(printDebug)
	  cout <<"running through the jet now - 1"<<endl;

	hGenJetpT->Fill(jetpT->at(ij), 1);
	hGenJetpT_wWeight->Fill(jetpT->at(ij), xsecweight);
	if(printDebug)
	  cout <<"running through the jet now - 2"<<endl;

	hGenJeteta->Fill(jeteta->at(ij), xsecweight);
	hGenJetphi->Fill(jetphi->at(ij), xsecweight);
	if(printDebug)
	  cout <<"running through the jet now - 3"<<endl;
	
	int ptbin=-1;
	for(int i = 0; i<npubptbins; ++i){
	  if(jetpT->at(ij) > pubptbins[i])
	    ptbin = i;
	}
	if(ptbin == -1)
	  continue;

	if(printDebug){
	  std::cout<<    "pt bin = "<<ptbin<<std::endl;
	}

	if(zg->at(ij) > 0){
	  hGen_Zg[ptbin]->Fill(zg->at(ij), xsecweight);
	  hGen_Rg[ptbin]->Fill(rg->at(ij), xsecweight);
	}
	hGen_M[ptbin]->Fill(jetM->at(ij), xsecweight);

	vector<double> reczg = rec_zg->at(ij);
	vector<double> recrg = rec_rg->at(ij);
	vector<double> recsdpt = sdjetpT->at(ij);

	if(printDebug){
	  std::cout<<"    There are a total of "<<reczg.size()<<" zg splits in this jet"<<std::endl;
	  std::cout<<"    There are a total of "<<recrg.size()<<" rg splits in this jet"<<std::endl;
	}      
      
	for(unsigned in = 0; in<reczg.size(); ++in){
	  if(printDebug){
	    std::cout<<"           split no.  "<<in<<" zg = "<<reczg.at(in)<<std::endl;
	    std::cout<<"           split no.  "<<in<<" rg = "<<recrg.at(in)<<std::endl;
	    std::cout<<"           split no.  "<<in<<" pt = "<<recsdpt.at(in)<<std::endl;
	  }

	  hGen_Zg_vs_Split[ptbin]->Fill(in+1, reczg.at(in), mcweight);
	  hGen_Rg_vs_Split[ptbin]->Fill(in+1, recrg.at(in), mcweight);
	  
	  int gptbin=-1;
	  for(int ig = 0; ig<ngroompTbins; ++ig){
	    if(recsdpt.at(in) > groompTbins[ig])
	      gptbin = ig;
	  }
	  if(gptbin == -1)
	    continue;
	  
	  hGenZg_Split_InitpTbins[gptbin]->Fill(in+1, reczg.at(in), mcweight);
	  hGenRg_Split_InitpTbins[gptbin]->Fill(in+1, recrg.at(in), mcweight);
	  
	}
	
      }

    }
    
    if(doHC){
      
      cout<<"total number of events in HT tree = "<<HCjetTree->GetEntries()<<endl;
      
      //! loop over the events and fill the histograms
      for(long Nev = 0; Nev<HCjetTree->GetEntries(); ++Nev){
	
	HCjetTree->GetEntry(Nev);
	if(Nev%10000 == 0)
	  std::cout<<"HT Event Number "<<Nev<<std::endl;
      
	if(printDebug)
	  std::cout<<"There are "<<HCjetpT->size()<<" high tower Jets in this event"<<std::endl;
	
	if(mctype == "herwig7"){
	  if(HCpthat < pthatbins[inpt])
	    continue;
	}else if(mctype == "pythia8"){
	  if(HCpthat >= pthatbins[inpt+1] || HCpthat < pthatbins[inpt])
	    continue;
	}

	double xsecweight = HCmcweight;
	
	//! loop over the jets in that event 
	for(unsigned ij = 0; ij<HCjetpT->size(); ++ij){
	  
	  if(printDebug){
	    std::cout<<"    hard core Jet pt = "<<HCjetpT->at(ij)<<std::endl;
	  }

	  if(printDebug)
	    std::cout<<"Inside the hard core jet loop"<<std::endl;


	  if(mctype == "herwig7")	  
	    if(HCjetpT->at(ij) > 3*pthatbins[inpt+1])
	      continue;
	  
	  if(mctype == "pythia8")
	    if(HCjetpT->at(ij) > 3*pthatbins[inpt+1])
	      continue;

	  if(printDebug)
	    cout <<"running through the HC jet now - 1"<<endl;

	  hHCGenJetpT->Fill(HCjetpT->at(ij), 1);
	  hHCGenJetpT_wWeight->Fill(HCjetpT->at(ij), xsecweight);
	  if(printDebug)
	    cout <<"running through the HC jet now - 2"<<endl;

	  hHCGenJeteta->Fill(HCjeteta->at(ij), xsecweight);
	  hHCGenJetphi->Fill(HCjetphi->at(ij), xsecweight);
	  if(printDebug)
	    cout <<"running through the HC jet now - 3"<<endl;
	
	  int ptbin=-1;
	  for(int i = 0; i<npubptbins; ++i){
	    if(HCjetpT->at(ij) > pubptbins[i])
	      ptbin = i;
	  }
	  if(ptbin == -1)
	    continue;

	  if(printDebug){
	    std::cout<<    "     HCzg->size() =  "<<HCzg->size()<<std::endl;
	    std::cout<<    "     zg  = "<<HCzg->at(ij)<<std::endl;
	  }

	  if(printDebug){
	    std::cout<<    "pt bin = "<<ptbin<<std::endl;
	  }

	  if(HCzg->at(ij) > 0){
	    hHCGen_Zg[ptbin]->Fill(HCzg->at(ij), xsecweight);
	    hHCGen_Rg[ptbin]->Fill(HCrg->at(ij), xsecweight);
	    if(printDebug){
	      std::cout<<    "     rg  = "<<HCrg->at(ij)<<std::endl;
	    }
	  }
	  hHCGen_M[ptbin]->Fill(HCjetM->at(ij), xsecweight);
	  if(printDebug){
	    std::cout<<    "     mj  = "<<HCjetM->at(ij)<<std::endl;
	  }

	}

      }

    }
    //! doHC

    if(doMC){
      
      cout<<"total number of events in Matched tree = "<<MCjetTree->GetEntries()<<endl;
      
      //! loop over the events and fill the histograms
      for(long Nev = 0; Nev<MCjetTree->GetEntries(); ++Nev){
	
	MCjetTree->GetEntry(Nev);
	if(Nev%10000 == 0)
	  std::cout<<"Matched Event Number "<<Nev<<std::endl;
      
	if(printDebug)
	  std::cout<<"There are "<<MCjetpT->size()<<" matched Jets in this event"<<std::endl;
	
	if(mctype == "herwig7"){
	  if(MCpthat < pthatbins[inpt])
	    continue;
	}else if(mctype == "pythia8"){
	  if(MCpthat >= pthatbins[inpt+1] || MCpthat < pthatbins[inpt])
	    continue;
	}

	double xsecweight = MCmcweight;
	
	//! loop over the jets in that event 
	for(unsigned ij = 0; ij<MCjetpT->size(); ++ij){
	  
	  if(printDebug){
	    std::cout<<"    matched Jet pt = "<<MCjetpT->at(ij)<<std::endl;
	  }

	  if(printDebug)
	    std::cout<<"Inside the matched jet loop"<<std::endl;


	  if(mctype == "herwig7")	  
	    if(MCjetpT->at(ij) > 3*pthatbins[inpt+1])
	      continue;
	  
	  if(mctype == "pythia8")
	    if(MCjetpT->at(ij) > 3*pthatbins[inpt+1])
	      continue;

	  if(printDebug)
	    cout <<"running through the MC jet now - 1"<<endl;

	  hMCGenJetpT->Fill(MCjetpT->at(ij), 1);
	  hMCGenJetpT_wWeight->Fill(MCjetpT->at(ij), xsecweight);
	  if(printDebug)
	    cout <<"running through the MC jet now - 2"<<endl;

	  hMCGenJeteta->Fill(MCjeteta->at(ij), xsecweight);
	  hMCGenJetphi->Fill(MCjetphi->at(ij), xsecweight);
	  if(printDebug)
	    cout <<"running through the MC jet now - 3"<<endl;
	
	  int ptbin=-1;
	  for(int i = 0; i<npubptbins; ++i){
	    if(MCjetpT->at(ij) > pubptbins[i])
	      ptbin = i;
	  }
	  if(ptbin == -1)
	    continue;

	  if(printDebug){
	    std::cout<<    "pt bin = "<<ptbin<<std::endl;
	  }

	  if(MCzg->at(ij) > 0){
	    hMCGen_Zg[ptbin]->Fill(MCzg->at(ij), xsecweight);
	    hMCGen_Rg[ptbin]->Fill(MCrg->at(ij), xsecweight);
	  }
	  hMCGen_M[ptbin]->Fill(MCjetM->at(ij), xsecweight);

	}
      }
    }//! doMC
  }//! pthat bin loop
  cout << "total triggers! " << tot_trigs << endl;
  
  fout->Write();
  fout->Close();
  
  
  for(unsigned inpt = 0; inpt < npthatbins; ++inpt){

    // xsec_avg[inpt];
    // nevts_pthatbin[inpt];

    std::cout<<"pthatbin = ["<<pthatbins[inpt]<<", "<<pthatbins[inpt+1]<<"]"<<std::endl;
    if(nevts_pthatbin[inpt]>0)
      xsec_avg[inpt] = xsec_avg[inpt]/nevts_pthatbin[inpt];
    std::cout<<"     average x-sec = "<<xsec_avg[inpt]<<", with a total of "<<nevts_pthatbin[inpt]<<" events"<<std::endl;
    if(nevts_pthatbin[inpt]>0)
      xsec_avg[inpt] = xsec_avg[inpt]/nevts_pthatbin[inpt];

  }

  for(unsigned inpt = 0; inpt < npthatbins; ++inpt){
    std::cout<<xsec_avg[inpt]<<", ";
  }
  std::cout<<std::endl;
  
}
