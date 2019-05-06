#include "RooUnfold.h"
#include <string>
#include <iostream>
#include "math.h"

using namespace std;

void HistsFromTreeHer(TFile *file, const char * treestring, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D) {
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
void HistsFromTreeP6HerCompare(TFile *file, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D) {
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

void HerwigTreetoHist() {
  string dir = "~/jetmass/";
  string herin = "production/Results/";
  string file_undec = "herwig7_decays_off.root";
  string out = "~/jetmass/production/macros/herwig_hists/";
  string filetype = ".pdf";

  TFile* her = new TFile( (dir + herin + file_dec).c_str(), "READ");

  //making variable bin size histogram
  const int nBins_pt = 8;
  double edges[nBins_pt + 1] = {5,10,15,20,25,30,40,60,100};
  
  TH1D* m_her = new TH1D("m_her",";M [GeV/c^{2}]; arb.",10,0,10);
  TH1D* pt_coarse_her = new TH1D("pt_coarse_her",";p_{T} [GeV/c];arb.", nBins_pt, edges);
  TH1D* pt_fine_her = new TH1D("pt_fine_her",";p_{T} [GeV/c];arb.", 20, 0, 80);
  
  TH2D* m_v_pt_her = new TH2D("m_v_pt_her",";M [GeV/c^{2}];p_{T} [GeV/c]", 10, 0, 10, 15,5,80);

  TH1D* eta_her = new TH1D("eta_her","",25,-1,1);
  TH2D* eta_her_v_pt = new TH2D("eta_her_v_pt","",25,-1,1,15,5,80);
  
  TH1D* phi_her = new TH1D("phi_her","",24,0,2*M_PI);
  TH2D* phi_her_v_pt = new TH2D("phi_her_v_pt","",24,0,2*M_PI,15,5,80);

  TH1D* conspT = new TH1D("conspTher","",30,0.2,30.2);
  TH1D* consDist = new TH1D("consDisther","",20,0,0.4);
  TH1D* consGirth = new TH1D("consGirthher","",20,0,0.12);
  TH1D* jetMult = new TH1D("jetMulther","",25,0,25);
  TH1D* jetGirth = new TH1D("jetGirthher","",20,0,0.3);

  TH2D* conspT_v_pt = new TH2D("conspT_v_pther","",30,0.2,30.2,15,5,80);
  TH2D* consDist_v_pt = new TH2D("consDist_v_pther","",20,0,0.4,15,5,80);
  TH2D* consGirth_v_pt = new TH2D("consGirth_v_pther","",20,0,0.12,15,5,80);
  TH2D* jetMult_v_pt = new TH2D("jetMult_v_pther","",25,0,25,15,5,80);
  TH2D* jetGirth_v_pt = new TH2D("jetGirth_v_pther","",20,0,0.3,15,5,80);

  
  vector<TH1D*> her_hists1D = {m_her, pt_coarse_her, pt_fine_her,eta_her,phi_her}; vector<TH2D*> her_hists2D = {m_v_pt_her,eta_her_v_pt,phi_her_v_pt};  
  
  vector<TH1D*> her_1D = {consDist, consGirth, conspT, jetMult, jetGirth};
  vector<TH2D*> her_2D = {consDist_v_pt, consGirth_v_pt, conspT_v_pt, jetMult_v_pt, jetGirth_v_pt};
  
  //  HistsFromTreeP6HerCompare(her, her_1D, her_2D);  
  
  HistsFromTreeHer(her, "ResultTree", her_hists1D, her_hists2D);

  TFile *fout = new TFile( ( out + "hists.root" ).c_str() ,"RECREATE");
  
  for (int i = 0; i < her_hists1D.size(); ++ i) {
    her_hists1D[i]->Write();
  }
  for (int i = 0; i < her_hists2D.size(); ++ i) {
    her_hists2D[i]->Write();
  }
  
  for (int i = 0; i < her_1D.size(); ++ i) {
    her_1D[i]->Write(); her_2D[i]->Write();
  }
    
  fout->Close();
  
  return;
}
