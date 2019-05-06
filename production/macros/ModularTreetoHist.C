#include "RooUnfold.h"
#include <string>
#include <iostream>
#include "math.h"

using namespace std;

void HistsFromTreePL(TFile *file, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D, std::vector<TH3D*> hists3D) {
  vector<double> *Pt = 0; vector<double> *M = 0; vector<double> *Mg = 0; vector<double> *Zg = 0;
  double weight = 1; double n_jets = 0;
  
  TTree *t = (TTree*) file->Get("virtPartonTree");
  
  t->SetBranchAddress("virtPLpt",&Pt);
  t->SetBranchAddress("virtPLm",&M);
  t->SetBranchAddress("virtPLmg",&Mg);
  t->SetBranchAddress("virtPLzg",&Zg);
  
  cout << t->GetEntries() << endl;
  for (int i = 0; i < t->GetEntries(); ++ i) {
    if (i % 100000 == 0) { cout << "still chuggin. event " << i << endl;}
    t->GetEntry(i);
    for (int j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size   
      hists2D[0]->Fill(M->at(j), Pt->at(j), weight);
      if (Zg->at(j) >= 0.1) {
        hists2D[1]->Fill(Mg->at(j), Pt->at(j), weight);
      }
    }
  }
  t->ResetBranchAddresses();
  return;
}

void HistsFromTree(TFile *file, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D, std::vector<TH3D*> hists3D) {
  vector<double> *Pt = 0; vector<double> *M = 0;
  vector<double> *Eta = 0; vector<double> *Phi = 0;
  vector<double> *Mg = 0; vector<double> *Rg = 0; vector<double> *Zg = 0;
  double weight = 1; double n_jets = 0;
  
  TTree *t = (TTree*) file->Get("ResultTree");
  
  t->SetBranchAddress("jetpT",&Pt); t->SetBranchAddress("jetM",&M);
  t->SetBranchAddress("jeteta",&Eta); t->SetBranchAddress("jetphi",&Phi);
  t->SetBranchAddress("sdjetM",&Mg); t->SetBranchAddress("zg",&Zg); t->SetBranchAddress("rg",&Rg);
  t->SetBranchAddress("mcweight",&weight);

  cout << t->GetEntries() << endl;
  for (int i = 0; i < t->GetEntries(); ++ i) {
    if (i % 100000 == 0) { cout << "still chuggin. event " << i << endl;}
    t->GetEntry(i);
    for (int j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size      
      hists2D[0]->Fill(M->at(j), Pt->at(j), weight);
      hists2D[2]->Fill(Eta->at(j), Pt->at(j), weight);
      hists2D[3]->Fill(Phi->at(j), Pt->at(j), weight);
      if (Zg->at(j) >= 0.1) {
	hists2D[1]->Fill(Mg->at(j), Pt->at(j), weight);
	hists2D[4]->Fill(Zg->at(j), Pt->at(j), weight);
	hists2D[5]->Fill(Rg->at(j), Pt->at(j), weight);
	hists3D[0]->Fill(M->at(j), Mg->at(j), Pt->at(j), weight);
	hists3D[1]->Fill(M->at(j), Zg->at(j), Pt->at(j), weight);
	hists3D[2]->Fill(M->at(j), Rg->at(j), Pt->at(j), weight);
	hists3D[3]->Fill(Mg->at(j), Zg->at(j), Pt->at(j), weight);
	hists3D[4]->Fill(Mg->at(j), Rg->at(j), Pt->at(j), weight);
	hists3D[5]->Fill(Zg->at(j), Rg->at(j), Pt->at(j), weight);
      }
    }
  }
  t->ResetBranchAddresses();
  return;
}


void HistsFromTreeP6off(TFile *file, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D, std::vector<TH3D*> hists3D) {
  vector<double> *Pt = 0; vector<double> *M = 0;
  vector<double> *Eta = 0; vector<double> *Phi = 0;
  vector<double> *Mg = 0; vector<double> *Rg = 0; vector<double> *Zg = 0;
  double weight = 1; double n_jets = 0;
  
  TTree *t = (TTree*) file->Get("event");
  
  t->SetBranchAddress("Pt",&Pt); t->SetBranchAddress("M",&M);
  t->SetBranchAddress("Eta",&Eta); t->SetBranchAddress("Phi",&Phi);
  t->SetBranchAddress("mg",&Mg); t->SetBranchAddress("zg",&Zg); t->SetBranchAddress("rg",&Rg);
  t->SetBranchAddress("weight",&weight);

  cout << t->GetEntries() << endl;
  for (int i = 0; i < t->GetEntries(); ++ i) {
    if (i % 100000 == 0) { cout << "still chuggin. event " << i << endl;}
    t->GetEntry(i);
    for (int j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size      
      hists2D[0]->Fill(M->at(j), Pt->at(j), weight);
      hists2D[2]->Fill(Eta->at(j), Pt->at(j), weight);
      hists2D[3]->Fill(Phi->at(j), Pt->at(j), weight);
      if (Zg->at(j) >= 0.1) {
	hists2D[1]->Fill(Mg->at(j), Pt->at(j), weight);
	hists2D[4]->Fill(Zg->at(j), Pt->at(j), weight);
	hists2D[5]->Fill(Rg->at(j), Pt->at(j), weight);
	hists3D[0]->Fill(M->at(j), Mg->at(j), Pt->at(j), weight);
	hists3D[1]->Fill(M->at(j), Zg->at(j), Pt->at(j), weight);
	hists3D[2]->Fill(M->at(j), Rg->at(j), Pt->at(j), weight);
	hists3D[3]->Fill(Mg->at(j), Zg->at(j), Pt->at(j), weight);
	hists3D[4]->Fill(Mg->at(j), Rg->at(j), Pt->at(j), weight);
	hists3D[5]->Fill(Zg->at(j), Rg->at(j), Pt->at(j), weight);
      }
    }
  }
  t->ResetBranchAddresses();
  return;
}


void ModularTreetoHist() {
  cout << "a" << endl;
  string loc = "~/jetmass/production/Results/";
  string locp6off = "~/jetmass/out/sim/py/";
  string p8onin = "pythia8_decays_on.root";
  string p8offin = "pythia8_decays_off_and_noHad_R02_isaac2.root";//"pythia8_decays_off_R04.root"; //!
  string p6onin = "py6_decayed_jewel_pthatbin580_R04.root";
  string p6offin = "full_w_o_bin_drop_R02.root"; //!
  string h7onin = "herwig7_decays_on.root";
  string h7offin = "herwig7_decays_off_R02.root"; //!
  
  string out = "~/jetmass/production/macros/hists/";
  string outname = "hists_allsim_lowzgremoved_R02.root"; //!

  TFile *fp8on = new TFile((loc+p8onin).c_str(),"READ");
  TFile *fp8off = new TFile((loc+p8offin).c_str(),"READ");
  TFile *fp6on = new TFile((loc+p6onin).c_str(),"READ");
  TFile *fp6off = new TFile((locp6off+p6offin).c_str(),"READ");
  TFile *fh7on = new TFile((loc+h7onin).c_str(),"READ");
  TFile *fh7off = new TFile((loc+h7offin).c_str(),"READ");

  //TH1Ds!
  TH1D* hdummy = new TH1D("hdummy","",1,0,1); TH3D* hdummy3D = new TH3D("hdummy3D","",1,0,1,1,0,1,1,0,1);
  
  //TH2Ds!

  TH2D* mvpt_p8on = new TH2D("mvpt_p8on",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mvpt_p8off = new TH2D("mvpt_p8off",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mvpt_p6on = new TH2D("mvpt_p6on",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mvpt_p6off = new TH2D("mvpt_p6off",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mvpt_h7on = new TH2D("mvpt_h7on",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mvpt_h7off = new TH2D("mvpt_h7off",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mvpt_p8offPL = new TH2D("mvpt_p8offPL",";M [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  
  TH2D* mgvpt_p8on = new TH2D("mgvpt_p8on",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mgvpt_p8off = new TH2D("mgvpt_p8off",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mgvpt_p6on = new TH2D("mgvpt_p6on",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mgvpt_p6off = new TH2D("mgvpt_p6off",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mgvpt_h7on = new TH2D("mgvpt_h7on",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mgvpt_h7off = new TH2D("mgvpt_h7off",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  TH2D* mgvpt_p8offPL = new TH2D("mgvpt_p8offPL",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]",14,0,14,15,5,80);
  
  TH2D* etavpt_p8on = new TH2D("etavpt_p8on",";#eta;p_{T} [GeV/c]",50,-1,1,15,5,80);
  TH2D* etavpt_p8off = new TH2D("etavpt_p8off",";#eta;p_{T} [GeV/c]",50,-1,1,15,5,80);
  TH2D* etavpt_p6on = new TH2D("etavpt_p6on",";#eta;p_{T} [GeV/c]",50,-1,1,15,5,80);
  TH2D* etavpt_p6off = new TH2D("etavpt_p6off",";#eta;p_{T} [GeV/c]",50,-1,1,15,5,80);
  TH2D* etavpt_h7on = new TH2D("etavpt_h7on",";#eta;p_{T} [GeV/c]",50,-1,1,15,5,80);
  TH2D* etavpt_h7off = new TH2D("etavpt_h7off",";#eta;p_{T} [GeV/c]",50,-1,1,15,5,80);

  TH2D* phivpt_p8on = new TH2D("phivpt_p8on",";#phi;p_{T} [GeV/c]",50,0,2*M_PI,15,5,80);
  TH2D* phivpt_p8off = new TH2D("phivpt_p8off",";#phi;p_{T} [GeV/c]",50,0,2*M_PI,15,5,80);
  TH2D* phivpt_p6on = new TH2D("phivpt_p6on",";#phi;p_{T} [GeV/c]",50,0,2*M_PI,15,5,80);
  TH2D* phivpt_p6off = new TH2D("phivpt_p6off",";#phi;p_{T} [GeV/c]",50,0,2*M_PI,15,5,80);
  TH2D* phivpt_h7on = new TH2D("phivpt_h7on",";#phi;p_{T} [GeV/c]",50,0,2*M_PI,15,5,80);
  TH2D* phivpt_h7off = new TH2D("phivpt_h7off",";#phi;p_{T} [GeV/c]",50,0,2*M_PI,15,5,80);
  
  TH2D* zgvpt_p8on = new TH2D("zgvpt_p8on",";z_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);
  TH2D* zgvpt_p8off = new TH2D("zgvpt_p8off",";z_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);
  TH2D* zgvpt_p6on = new TH2D("zgvpt_p6on",";z_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);
  TH2D* zgvpt_p6off = new TH2D("zgvpt_p6off",";z_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);
  TH2D* zgvpt_h7on = new TH2D("zgvpt_h7on",";z_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);
  TH2D* zgvpt_h7off = new TH2D("zgvpt_h7off",";z_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);

  TH2D* rgvpt_p8on = new TH2D("rgvpt_p8on",";R_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);
  TH2D* rgvpt_p8off = new TH2D("rgvpt_p8off",";R_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);
  TH2D* rgvpt_p6on = new TH2D("rgvpt_p6on",";R_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);
  TH2D* rgvpt_p6off = new TH2D("rgvpt_p6off",";R_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);
  TH2D* rgvpt_h7on = new TH2D("rgvpt_h7on",";R_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);
  TH2D* rgvpt_h7off = new TH2D("rgvpt_h7off",";R_{g} [GeV/c^{2}];p_{T} [GeV/c]",10,0,0.5,15,5,80);

  //TH3Ds!

  TH3D* mvmgvpt_p8on = new TH3D("mvmgvpt_p8on",";M [GeV/c^{2}];M_{g} [GeV/c^{2}]; p_{T} [GeV/c]",14,0,14,14,0,14,15,5,80);
  TH3D* mvmgvpt_p8off = new TH3D("mvmgvpt_p8off",";M [GeV/c^{2}];M_{g} [GeV/c^{2}]; p_{T} [GeV/c]",14,0,14,14,0,14,15,5,80);
  TH3D* mvmgvpt_p6on = new TH3D("mvmgvpt_p6on",";M [GeV/c^{2}];M_{g} [GeV/c^{2}]; p_{T} [GeV/c]",14,0,14,14,0,14,15,5,80);
  TH3D* mvmgvpt_p6off = new TH3D("mvmgvpt_p6off",";M [GeV/c^{2}];M_{g} [GeV/c^{2}]; p_{T} [GeV/c]",14,0,14,14,0,14,15,5,80);
  TH3D* mvmgvpt_h7on = new TH3D("mvmgvpt_h7on",";M [GeV/c^{2}];M_{g} [GeV/c^{2}]; p_{T} [GeV/c]",14,0,14,14,0,14,15,5,80);
  TH3D* mvmgvpt_h7off = new TH3D("mvmgvpt_h7off",";M [GeV/c^{2}];M_{g} [GeV/c^{2}]; p_{T} [GeV/c]",14,0,14,14,0,14,15,5,80);
  
  TH3D* mvzgvpt_p8on = new TH3D("mvzgvpt_p8on",";M [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mvzgvpt_p8off = new TH3D("mvzgvpt_p8off",";M [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mvzgvpt_p6on = new TH3D("mvzgvpt_p6on",";M [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mvzgvpt_p6off = new TH3D("mvzgvpt_p6off",";M [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mvzgvpt_h7on = new TH3D("mvzgvpt_h7on",";M [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mvzgvpt_h7off = new TH3D("mvzgvpt_h7off",";M [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  
  TH3D* mvrgvpt_p8on = new TH3D("mvrgvpt_p8on",";M [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mvrgvpt_p8off = new TH3D("mvrgvpt_p8off",";M [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mvrgvpt_p6on = new TH3D("mvrgvpt_p6on",";M [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mvrgvpt_p6off = new TH3D("mvrgvpt_p6off",";M [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mvrgvpt_h7on = new TH3D("mvrgvpt_h7on",";M [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mvrgvpt_h7off = new TH3D("mvrgvpt_h7off",";M [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  
  TH3D* mgvzgvpt_p8on = new TH3D("mgvzgvpt_p8on",";M_{g} [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mgvzgvpt_p8off = new TH3D("mgvzgvpt_p8off",";M_{g} [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mgvzgvpt_p6on = new TH3D("mgvzgvpt_p6on",";M_{g} [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mgvzgvpt_p6off = new TH3D("mgvzgvpt_p6off",";M_{g} [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mgvzgvpt_h7on = new TH3D("mgvzgvpt_h7on",";M_{g} [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mgvzgvpt_h7off = new TH3D("mgvzgvpt_h7off",";M_{g} [GeV/c^{2}];z_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  
  TH3D* mgvrgvpt_p8on = new TH3D("mgvrgvpt_p8on",";M_{g} [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mgvrgvpt_p8off = new TH3D("mgvrgvpt_p8off",";M_{g} [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mgvrgvpt_p6on = new TH3D("mgvrgvpt_p6on",";M_{g} [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mgvrgvpt_p6off = new TH3D("mgvrgvpt_p6off",";M_{g} [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mgvrgvpt_h7on = new TH3D("mgvrgvpt_h7on",";M_{g} [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  TH3D* mgvrgvpt_h7off = new TH3D("mgvrgvpt_h7off",";M_{g} [GeV/c^{2}];R_{g}; p_{T} [GeV/c]",14,0,14,10,0,0.5,15,5,80);
  
  TH3D* zgvrgvpt_p8on = new TH3D("zgvrgvpt_p8on",";z_{g};R_{g}; p_{T} [GeV/c]",10,0,0.5,10,0,0.5,15,5,80);
  TH3D* zgvrgvpt_p8off = new TH3D("zgvrgvpt_p8off",";z_{g};R_{g}; p_{T} [GeV/c]",10,0,0.5,10,0,0.5,15,5,80);
  TH3D* zgvrgvpt_p6on = new TH3D("zgvrgvpt_p6on",";z_{g};R_{g}; p_{T} [GeV/c]",10,0,0.5,10,0,0.5,15,5,80);
  TH3D* zgvrgvpt_p6off = new TH3D("zgvrgvpt_p6off",";z_{g};R_{g}; p_{T} [GeV/c]",10,0,0.5,10,0,0.5,15,5,80);
  TH3D* zgvrgvpt_h7on = new TH3D("zgvrgvpt_h7on",";z_{g};R_{g}; p_{T} [GeV/c]",10,0,0.5,10,0,0.5,15,5,80);
  TH3D* zgvrgvpt_h7off = new TH3D("zgvrgvpt_h7off",";z_{g};R_{g}; p_{T} [GeV/c]",10,0,0.5,10,0,0.5,15,5,80);
  
  //vectors!

  vector<TH1D*> hists1D_dummy = {hdummy}; vector<TH3D*> hists3D_dummy = {hdummy3D};
  
  vector<TH2D*> hists2D_p8on = {mvpt_p8on,mgvpt_p8on,etavpt_p8on,phivpt_p8on,zgvpt_p8on,rgvpt_p8on};
  vector<TH2D*> hists2D_p8off = {mvpt_p8off,mgvpt_p8off,etavpt_p8off,phivpt_p8off,zgvpt_p8off,rgvpt_p8off};
  vector<TH2D*> hists2D_p6on = {mvpt_p6on,mgvpt_p6on,etavpt_p6on,phivpt_p6on,zgvpt_p6on,rgvpt_p6on};
  vector<TH2D*> hists2D_p6off = {mvpt_p6off,mgvpt_p6off,etavpt_p6off,phivpt_p6off,zgvpt_p6off,rgvpt_p6off};
  vector<TH2D*> hists2D_h7on = {mvpt_h7on,mgvpt_h7on,etavpt_h7on,phivpt_h7on,zgvpt_h7on,rgvpt_h7on};
  vector<TH2D*> hists2D_h7off = {mvpt_h7off,mgvpt_h7off,etavpt_h7off,phivpt_h7off,zgvpt_h7off,rgvpt_h7off};
  vector<TH2D*> hists2D_p8offPL = {mvpt_p8offPL,mgvpt_p8offPL};

  vector<TH3D*> hists3D_p8on = {mvmgvpt_p8on,mvzgvpt_p8on,mvrgvpt_p8on,mgvzgvpt_p8on,mgvrgvpt_p8on,zgvrgvpt_p8on};
  vector<TH3D*> hists3D_p8off = {mvmgvpt_p8off,mvzgvpt_p8off,mvrgvpt_p8off,mgvzgvpt_p8off,mgvrgvpt_p8off,zgvrgvpt_p8off};
  vector<TH3D*> hists3D_p6on = {mvmgvpt_p6on,mvzgvpt_p6on,mvrgvpt_p6on,mgvzgvpt_p6on,mgvrgvpt_p6on,zgvrgvpt_p6on};
  vector<TH3D*> hists3D_p6off = {mvmgvpt_p6off,mvzgvpt_p6off,mvrgvpt_p6off,mgvzgvpt_p6off,mgvrgvpt_p6off,zgvrgvpt_p6off};
  vector<TH3D*> hists3D_h7on = {mvmgvpt_h7on,mvzgvpt_h7on,mvrgvpt_h7on,mgvzgvpt_h7on,mgvrgvpt_h7on,zgvrgvpt_h7on};
  vector<TH3D*> hists3D_h7off = {mvmgvpt_h7off,mvzgvpt_h7off,mvrgvpt_h7off,mgvzgvpt_h7off,mgvrgvpt_h7off,zgvrgvpt_h7off};
  
  cout << "b" << endl;
  // HistsFromTree(fp8on,hists1D_dummy, hists2D_p8on, hists3D_p8on);
  cout << "c" << endl;
  HistsFromTree(fp8off,hists1D_dummy, hists2D_p8off, hists3D_p8off);
  cout << "d" << endl;
  //HistsFromTree(fp6on,hists1D_dummy, hists2D_p6on, hists3D_p6on);
  cout << "e" << endl;
  HistsFromTreeP6off(fp6off,hists1D_dummy, hists2D_p6off, hists3D_p6off);
  cout << "f" << endl;
  //HistsFromTree(fh7on,hists1D_dummy, hists2D_h7on, hists3D_h7on);
  cout << "g" << endl;
  HistsFromTree(fh7off,hists1D_dummy, hists2D_h7off, hists3D_h7off);

  HistsFromTreePL(fp8off,hists1D_dummy,hists2D_p8offPL, hists3D_dummy);
  
  cout << "h" << endl;
  TFile *fout = new TFile((out+outname).c_str(),"RECREATE");

  for (int i = 0; i < hists2D_p8on.size(); ++ i) {
    //hists2D_p8on[i]->Write();
    hists2D_p8off[i]->Write();
    //hists2D_p6on[i]->Write();
    hists2D_p6off[i]->Write();
    //hists2D_h7on[i]->Write();
    hists2D_h7off[i]->Write();
  }

  for (int i = 0; i < hists2D_p8offPL.size(); ++ i) {
    hists2D_p8offPL[i]->Write();
  }

  for (int i = 0; i < hists3D_p8on.size(); ++ i) {
    //hists3D_p8on[i]->Write();
    hists3D_p8off[i]->Write();
    // hists3D_p6on[i]->Write();
    hists3D_p6off[i]->Write();
    //hists3D_h7on[i]->Write();
    hists3D_h7off[i]->Write();
  }

  cout << "i" << endl;
  fout->Close();
  
  return;
}
