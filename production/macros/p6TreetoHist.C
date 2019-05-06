#include "RooUnfold.h"
#include <string>
#include <iostream>
#include "math.h"

using namespace std;

void HistsFromTreeP6dec(TFile *file, const char * treestring, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D) {
  vector<double> *Pt = 0; vector<double> *M = 0;
  vector<double> *Eta = 0; vector<double> *Phi = 0;
  double weight = 1; double n_jets = 0;
  
  TTree *t = (TTree*) file->Get(treestring);
  
  t->SetBranchAddress("jetpT",&Pt); t->SetBranchAddress("jetM",&M);
  t->SetBranchAddress("jeteta",&Eta); t->SetBranchAddress("jetphi",&Phi);
  t->SetBranchAddress("mcweight",&weight);
  
  //cout << t->GetEntries() << endl;
  for (int i = 0; i < t->GetEntries(); ++ i) {
    //if (i % 1000000 == 0) { cout << "still chuggin. " << i << endl;}
    t->GetEntry(i);
    for (int j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size      
      hists1D[0]->Fill(M->at(j), weight);
      hists1D[1]->Fill(Pt->at(j), weight);
      hists1D[2]->Fill(Pt->at(j), weight);
      hists1D[3]->Fill(Eta->at(j), weight);
      hists1D[4]->Fill(Phi->at(j), weight);
      //      if (Eta->at(j) < 0.00001 && Eta->at(j) > -0.00001) {
      hists2D[0]->Fill(M->at(j), Pt->at(j), weight);
      hists2D[1]->Fill(Eta->at(j), Pt->at(j), weight);
      hists2D[2]->Fill(Phi->at(j), Pt->at(j), weight);
	//}
    }
  }
  t->ResetBranchAddresses();
  return;
}

//consdist consgirth jetmult jetgirth
void HistsFromTreeP6P6decCompare(TFile *file, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D) {
  vector<double> *Pt = 0; vector<double> *jetMult = 0;
  vector<double> *jetGirth = 0; vector<vector<double> > *consDist = 0;
  vector<vector<double> > *consGirth = 0;
  vector<vector<double> > *conspT = 0;
  
  double weight = 1; double n_jets = 0;

  TTree *t = (TTree*) file->Get("ResultTree");
  
  t->SetBranchAddress("jetpT",&Pt); t->SetBranchAddress("conspT",&conspT);
  t->SetBranchAddress("consGirth",&consGirth); t->SetBranchAddress("consDist",&consDist);
  t->SetBranchAddress("jetGirth",&jetGirth); t->SetBranchAddress("jetMult",&jetMult);
  t->SetBranchAddress("mcweight",&weight);

  cout << t->GetEntries() << endl;
  for (int i = 0; i < t->GetEntries(); ++ i) {
    if (i % 1000000 == 0) { cout << "still chuggin. " << i << endl;}
    t->GetEntry(i);
    for (int j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size                                           
      for (int k = 0; k < consGirth->at(j).size(); ++ k) {
        hists1D[0]->Fill(consDist->at(j).at(k),weight);
        hists1D[1]->Fill(consGirth->at(j).at(k),weight);
	hists1D[2]->Fill(conspT->at(j).at(k),weight);

        hists2D[0]->Fill(consDist->at(j).at(k), Pt->at(j), weight);
        hists2D[1]->Fill(consGirth->at(j).at(k), Pt->at(j), weight);
	hists2D[2]->Fill(conspT->at(j).at(k), Pt->at(j), weight);
      }
      hists1D[3]->Fill(jetMult->at(j),weight);
      hists1D[4]->Fill(jetGirth->at(j),weight);

      hists2D[3]->Fill(jetMult->at(j), Pt->at(j), weight);
      hists2D[4]->Fill(jetGirth->at(j), Pt->at(j), weight);
    }
  }
  t->ResetBranchAddresses();
  return;
}

void p6TreetoHist() {
  string dir = "~/jetmass/";
  string p6in = "production/";
  string file_dec = "py6_decayed_jewel_pthatbin1080_R04.root";
  string out = "~/jetmass/production/macros/p6_dec_hists/";
  string filetype = ".pdf";

  TFile* p6_decayed = new TFile( (dir + p6in + file_dec).c_str(), "READ");

  //making variable bin size histogram
  const int nBins_pt = 8;
  double edges[nBins_pt + 1] = {5,10,15,20,25,30,40,60,100};
  
  TH1D* m_un = new TH1D("m_un",";M [GeV/c^{2}]; arb.",10,0,10);
  TH1D* pt_coarse_un = new TH1D("pt_coarse_un",";p_{T} [GeV/c];arb.", nBins_pt, edges);
  TH1D* pt_fine_un = new TH1D("pt_fine_un",";p_{T} [GeV/c];arb.", 20, 0, 80);
  TH1D* m_effic_un = new TH1D("m_effic_un",";M [GeV/c^{2}]; arb.",10,0,10);
  TH1D* pt_coarse_effic_un = new TH1D("pt_coarse_effic_un",";p_{T} [GeV/c];arb.", nBins_pt, edges);
  TH1D* pt_fine_effic_un = new TH1D("pt_fine_effic_un",";p_{T} [GeV/c];arb.", 20, 0, 80);
  
  TH2D* m_v_pt_un = new TH2D("m_v_pt_un",";M [GeV/c^{2}];p_{T} [GeV/c]", 10, 0, 10, 15,5,80);
  TH2D* m_v_pt_effic_un = new TH2D("m_v_pt_effic_un",";M [GeV/c^{2}];p_{T} [GeV/c]", 10, 0, 10, 15,5,80);
  
  TH1D* m_dec = new TH1D("m_dec",";M [GeV/c^{2}]; arb.",10,0,10);
  TH1D* pt_coarse_dec = new TH1D("pt_coarse_dec",";p_{T} [GeV/c];arb.", nBins_pt, edges);
  TH1D* pt_fine_dec = new TH1D("pt_fine_dec",";p_{T} [GeV/c];arb.", 20, 0, 80);
  TH1D* m_effic_dec = new TH1D("m_effic_dec",";M [GeV/c^{2}]; arb.",10,0,10);
  TH1D* pt_coarse_effic_dec = new TH1D("pt_coarse_effic_dec",";p_{T} [GeV/c];arb.", nBins_pt, edges);
  TH1D* pt_fine_effic_dec = new TH1D("pt_fine_effic_dec",";p_{T} [GeV/c];arb.", 20, 0, 80);
  
  TH2D* m_v_pt_dec = new TH2D("m_v_pt_dec",";M [GeV/c^{2}];p_{T} [GeV/c]", 10, 0, 10, 15,5,80);
  TH2D* m_v_pt_effic_dec = new TH2D("m_v_pt_effic_dec",";M [GeV/c^{2}];p_{T} [GeV/c]", 10, 0, 10, 15,5,80);

  TH1D* eta_dec = new TH1D("eta_dec","",50,-1,1);
  TH1D* eta_un = new TH1D("eta_un","",50,-1,1);
  TH1D* eta_effic_dec = new TH1D("eta_effic_dec","",50,-1,1);
  TH1D* eta_effic_un = new TH1D("eta_effic_un","",50,-1,1);

  TH2D* eta_dec_v_pt = new TH2D("eta_dec_v_pt","",50,-1,1,15,5,80);
  TH2D* eta_un_v_pt = new TH2D("eta_un_v_pt","",50,-1,1,15,5,80);
  
  TH1D* phi_dec = new TH1D("phi_dec","",50,0,2*M_PI);
  TH1D* phi_un = new TH1D("phi_un","",50,0,2*M_PI);
  TH1D* phi_effic_dec = new TH1D("phi_effic_dec","",25,0,2*M_PI);
  TH1D* phi_effic_un = new TH1D("phi_effic_un","",25,0,2*M_PI);

  TH2D* phi_dec_v_pt = new TH2D("phi_dec_v_pt","",24,0,2*M_PI,15,5,80);
  TH2D* phi_un_v_pt = new TH2D("phi_un_v_pt","",24,0,2*M_PI,15,5,80);

  TH1D* conspT = new TH1D("conspTp6dec","",30,0.2,30.2);
  TH1D* consDist = new TH1D("consDistp6dec","",20,0,0.4);
  TH1D* consGirth = new TH1D("consGirthp6dec","",20,0,0.12);
  TH1D* jetMult = new TH1D("jetMultp6dec","",20,0,20);
  TH1D* jetGirth = new TH1D("jetGirthp6dec","",20,0,0.3);

  TH2D* conspT_v_pt = new TH2D("conspT_v_ptp6dec","",30,0.2,30.2,15,5,80);
  TH2D* consDist_v_pt = new TH2D("consDist_v_ptp6dec","",20,0,0.4,15,5,80);
  TH2D* consGirth_v_pt = new TH2D("consGirth_v_ptp6dec","",20,0,0.12,15,5,80);
  TH2D* jetMult_v_pt = new TH2D("jetMult_v_ptp6dec","",20,0,20,15,5,80);
  TH2D* jetGirth_v_pt = new TH2D("jetGirth_v_ptp6dec","",20,0,0.3,15,5,80);

  
  vector<TH1D*> un_hists1D = {m_un, pt_coarse_un, pt_fine_un,eta_un,phi_un}; vector<TH2D*> un_hists2D = {m_v_pt_un,eta_un_v_pt,phi_un_v_pt};  
  vector<TH1D*> dec_hists1D = {m_dec, pt_coarse_dec, pt_fine_dec,eta_dec,phi_dec}; vector<TH2D*> dec_hists2D = {m_v_pt_dec,eta_dec_v_pt,phi_dec_v_pt};
  vector<TH1D*> un_effic_hists1D = {m_effic_un, pt_coarse_effic_un, pt_fine_effic_un,eta_effic_un}; vector<TH2D*> un_effic_hists2D = {m_v_pt_effic_un};  
  vector<TH1D*> dec_effic_hists1D = {m_effic_dec, pt_coarse_effic_dec, pt_fine_effic_dec,eta_effic_dec}; vector<TH2D*> dec_effic_hists2D = {m_v_pt_effic_dec};
  
  vector<TH1D*> p6_1D = {consDist, consGirth, conspT, jetMult, jetGirth};
  vector<TH2D*> p6_2D = {consDist_v_pt, consGirth_v_pt, conspT_v_pt, jetMult_v_pt, jetGirth_v_pt};
  
  HistsFromTreeP6P6decCompare(p6_decayed, p6_1D, p6_2D);  
  
  HistsFromTreeP6dec(p6_decayed, "ResultTree", dec_hists1D, dec_hists2D);

  TFile *fout = new TFile( ( out + "hists.root" ).c_str() ,"RECREATE");
  
  for (int i = 0; i < un_hists1D.size(); ++ i) {
    un_hists1D[i]->Write(); dec_hists1D[i]->Write();
    //un_effic_hists1D[i]->Write(); dec_effic_hists1D[i]->Write();
  }
  for (int i = 0; i < un_hists2D.size(); ++ i) {
    un_hists2D[i]->Write(); dec_hists2D[i]->Write();
    //un_effic_hists2D[i]->Write(); dec_effic_hists2D[i]->Write();
  }
  
  for (int i = 0; i < p6_1D.size(); ++ i) {
    p6_1D[i]->Write(); p6_2D[i]->Write();
  }
    
  fout->Close();
  
  return;
}
