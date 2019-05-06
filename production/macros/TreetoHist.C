#include "RooUnfold.h"
#include <string>
#include <iostream>
#include "math.h"

using namespace std;

void HistsFromTreeP8(TFile *file, const char * treestring, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D, std::vector<TH3D*> hists3D) {
  vector<double> *Pt = 0; vector<double> *M = 0;
  vector<double> *Eta = 0; vector<double> *Phi = 0;
  vector<double> *Mg = 0;
  vector<double> *qvg = 0;
  /*
  vector<double> *Tau1 = 0;
  vector<double> *Tau05 = 0;
  vector<double> *Tau0 = 0;
  vector<double> *Tau_05 = 0;
  vector<double> *Tau_1 = 0;
  */
  double weight = 1; double n_jets = 0;
  
  TTree *t = (TTree*) file->Get(treestring);
  
  t->SetBranchAddress("jetpT",&Pt); t->SetBranchAddress("jetM",&M);
  t->SetBranchAddress("jeteta",&Eta); t->SetBranchAddress("jetphi",&Phi);
  t->SetBranchAddress("sdjetM",&Mg);
  t->SetBranchAddress("qvg",&qvg);
  t->SetBranchAddress("mcweight",&weight);
  /*
  t->SetBranchAddress("jetTau1",&Tau1);
  t->SetBranchAddress("jetTau05",&Tau05);
  t->SetBranchAddress("jetTau0",&Tau0);
  t->SetBranchAddress("jetTau_05",&Tau_05);
  t->SetBranchAddress("jetTau_1",&Tau_1);
  */
  cout << t->GetEntries() << endl;
  for (int i = 0; i < t->GetEntries(); ++ i) {
    if (i % 100000 == 0) { cout << "still chuggin. " << i << endl;}
    t->GetEntry(i);
    for (int j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size      
      hists1D[0]->Fill(M->at(j), weight);
      hists1D[1]->Fill(Pt->at(j), weight);
      hists1D[2]->Fill(Pt->at(j), weight);
      hists1D[3]->Fill(Eta->at(j), weight);
      hists1D[4]->Fill(Phi->at(j), weight);
      hists2D[0]->Fill(M->at(j), Pt->at(j), weight);
      hists2D[1]->Fill(Mg->at(j), Pt->at(j), weight);
      hists2D[2]->Fill(Eta->at(j), Pt->at(j), weight);
      hists2D[3]->Fill(Phi->at(j), Pt->at(j), weight);
      if (qvg->at(j) == 0) {//q jet!
	hists2D[4]->Fill(M->at(j),Pt->at(j),weight);
      } 
      else if (qvg->at(j) == 1) {//g jet!
	hists2D[5]->Fill(M->at(j),Pt->at(j),weight);
      }
      else {//neither!
	hists2D[6]->Fill(M->at(j),Pt->at(j),weight);
      }

      //CHANGE BACK LATER! JUST NEED TO MAKE UNDECAYED AND DECAYED ON EQUAL FOOTING UNTIL I CAN RUN RIVET AGAIN ONCE IT IS FIXED
      hists3D[0]->Fill(log10(pow(M->at(j),2)*pow(Pt->at(j),-2)), log10(1/*Tau1->at(j)*/), Pt->at(j), weight);
      hists3D[1]->Fill(log10(pow(M->at(j),2)*pow(Pt->at(j),-2)), log10(1/*Tau05->at(j)*/), Pt->at(j), weight);
      hists3D[2]->Fill(log10(pow(M->at(j),2)*pow(Pt->at(j),-2)), log10(1/*Tau0->at(j)*/), Pt->at(j), weight);
      hists3D[3]->Fill(log10(pow(M->at(j),2)*pow(Pt->at(j),-2)), log10(1/*Tau_05->at(j)*/), Pt->at(j), weight);
      hists3D[4]->Fill(log10(pow(M->at(j),2)*pow(Pt->at(j),-2)), log10(1/*Tau_1->at(j)*/), Pt->at(j), weight);		       
    }
  }
  t->ResetBranchAddresses();
  return;
}

//consdist consgirth jetmult jetgirth
void HistsFromTreeP6P8Compare(TFile *file, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D) {
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

void TreetoHist() {
  string dir = "~/jetmass/";
  string p8in = "production/Results/";
  string file_dec = "pythia8_decays_on.root";
  string file_undec = "pythia8_decays_off_R06.root";
  string out = "~/jetmass/production/macros/hists/";
  string filetype = ".pdf";

  TFile* p8_decayed = new TFile( (dir + p8in + file_dec).c_str(), "READ");
  TFile* p8_undecayed = new TFile( (dir + p8in + file_undec).c_str(), "READ");

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

  //quark v gluon mass hists
  TH2D* m_v_pt_un_q = new TH2D("m_v_pt_un_q","",10,0,10,15,5,80);
  TH2D* m_v_pt_un_g = new TH2D("m_v_pt_un_g","",10,0,10,15,5,80);
  TH2D* m_v_pt_un_neither = new TH2D("m_v_pt_un_neither","",10,0,10,15,5,80);
  
  TH2D* m_v_pt_dec_q = new TH2D("m_v_pt_dec_q","",10,0,10,15,5,80);
  TH2D* m_v_pt_dec_g = new TH2D("m_v_pt_dec_g","",10,0,10,15,5,80);
  TH2D* m_v_pt_dec_neither = new TH2D("m_v_pt_dec_neither","",10,0,10,15,5,80);
  
  TH2D* mg_v_pt_un = new TH2D("mg_v_pt_un",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]", 10, 0, 10, 15,5,80);
  TH2D* mg_v_pt_effic_un = new TH2D("mg_v_pt_effic_un",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]", 10, 0, 10, 15,5,80);
  TH2D* mg_v_pt_dec = new TH2D("mg_v_pt_dec",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]", 10, 0, 10, 15,5,80);
  TH2D* mg_v_pt_effic_dec = new TH2D("mg_v_pt_effic_dec",";M_{g} [GeV/c^{2}];p_{T} [GeV/c]", 10, 0, 10, 15,5,80);


  TH1D* eta_dec = new TH1D("eta_dec","",25,-1,1);
  TH1D* eta_un = new TH1D("eta_un","",25,-1,1);
  TH1D* eta_effic_dec = new TH1D("eta_effic_dec","",25,-1,1);
  TH1D* eta_effic_un = new TH1D("eta_effic_un","",25,-1,1);

  TH2D* eta_dec_v_pt = new TH2D("eta_dec_v_pt","",25,-1,1,15,5,80);
  TH2D* eta_un_v_pt = new TH2D("eta_un_v_pt","",25,-1,1,15,5,80);
  
  TH1D* phi_dec = new TH1D("phi_dec","",24,0,2*M_PI);
  TH1D* phi_un = new TH1D("phi_un","",24,0,2*M_PI);
  TH1D* phi_effic_dec = new TH1D("phi_effic_dec","",24,0,2*M_PI);
  TH1D* phi_effic_un = new TH1D("phi_effic_un","",24,0,2*M_PI);

  TH2D* phi_dec_v_pt = new TH2D("phi_dec_v_pt","",24,0,2*M_PI,15,5,80);
  TH2D* phi_un_v_pt = new TH2D("phi_un_v_pt","",24,0,2*M_PI,15,5,80);

  TH1D* conspT = new TH1D("conspTp8","",30,0.2,30.2);
  TH1D* consDist = new TH1D("consDistp8","",20,0,0.4);
  TH1D* consGirth = new TH1D("consGirthp8","",20,0,0.12);
  TH1D* jetMult = new TH1D("jetMultp8","",25,0,25);
  TH1D* jetGirth = new TH1D("jetGirthp8","",20,0,0.3);

  TH2D* conspT_v_pt = new TH2D("conspT_v_ptp8","",30,0.2,30.2,15,5,80);
  TH2D* consDist_v_pt = new TH2D("consDist_v_ptp8","",20,0,0.4,15,5,80);
  TH2D* consGirth_v_pt = new TH2D("consGirth_v_ptp8","",20,0,0.12,15,5,80);
  TH2D* jetMult_v_pt = new TH2D("jetMult_v_ptp8","",25,0,25,15,5,80);
  TH2D* jetGirth_v_pt = new TH2D("jetGirth_v_ptp8","",20,0,0.3,15,5,80);

  TH1D* conspT_dec = new TH1D("conspTp8dec","",30,0.2,30.2);
  TH1D* consDist_dec = new TH1D("consDistp8dec","",20,0,0.4);
  TH1D* consGirth_dec = new TH1D("consGirthp8dec","",20,0,0.12);
  TH1D* jetMult_dec = new TH1D("jetMultp8dec","",25,0,25);
  TH1D* jetGirth_dec = new TH1D("jetGirthp8dec","",20,0,0.3);

  TH2D* conspT_v_pt_dec = new TH2D("conspT_v_ptp8dec","",30,0.2,30.2,15,5,80);
  TH2D* consDist_v_pt_dec = new TH2D("consDist_v_ptp8dec","",20,0,0.4,15,5,80);
  TH2D* consGirth_v_pt_dec = new TH2D("consGirth_v_ptp8dec","",20,0,0.12,15,5,80);
  TH2D* jetMult_v_pt_dec = new TH2D("jetMult_v_ptp8dec","",25,0,25,15,5,80);
  TH2D* jetGirth_v_pt_dec = new TH2D("jetGirth_v_ptp8dec","",20,0,0.3,15,5,80);

  TH3D* jetM_v_a1_v_pt_dec = new TH3D("jetM_v_a1_v_pt_dec","",30,-8,0,30,-8,0,15,5,80);
  TH3D* jetM_v_a05_v_pt_dec = new TH3D("jetM_v_a05_v_pt_dec","",30,-8,0,30,-8,0,15,5,80);
  TH3D* jetM_v_a0_v_pt_dec = new TH3D("jetM_v_a0_v_pt_dec","",30,-8,0,30,-8,0,15,5,80);
  TH3D* jetM_v_a_05_v_pt_dec = new TH3D("jetM_v_a_05_v_pt_dec","",30,-8,0,30,-8,0,15,5,80);
  TH3D* jetM_v_a_1_v_pt_dec = new TH3D("jetM_v_a_1_v_pt_dec","",30,-8,0,30,-8,0,15,5,80);
  
  TH3D* hdummy3D = new TH3D("hdummy3D","",1,0,1,1,0,1,1,0,1);
  
  vector<TH1D*> un_hists1D = {m_un, pt_coarse_un, pt_fine_un,eta_un,phi_un}; vector<TH2D*> un_hists2D = {m_v_pt_un, mg_v_pt_un, eta_un_v_pt,phi_un_v_pt, m_v_pt_un_q, m_v_pt_un_g, m_v_pt_un_neither};  
  vector<TH3D*> un_hists3D = {hdummy3D,hdummy3D,hdummy3D,hdummy3D,hdummy3D};
  vector<TH1D*> dec_hists1D = {m_dec, pt_coarse_dec, pt_fine_dec,eta_dec,phi_dec}; vector<TH2D*> dec_hists2D = {m_v_pt_dec, mg_v_pt_dec, eta_dec_v_pt,phi_dec_v_pt, m_v_pt_dec_q, m_v_pt_dec_g, m_v_pt_dec_neither};
  vector<TH3D*> dec_hists3D = {jetM_v_a1_v_pt_dec, jetM_v_a05_v_pt_dec, jetM_v_a0_v_pt_dec, jetM_v_a_05_v_pt_dec, jetM_v_a_1_v_pt_dec};
  vector<TH1D*> un_effic_hists1D = {m_effic_un, pt_coarse_effic_un, pt_fine_effic_un,eta_effic_un}; vector<TH2D*> un_effic_hists2D = {m_v_pt_effic_un, mg_v_pt_effic_un};  
  vector<TH1D*> dec_effic_hists1D = {m_effic_dec, pt_coarse_effic_dec, pt_fine_effic_dec,eta_effic_dec}; vector<TH2D*> dec_effic_hists2D = {m_v_pt_effic_dec, mg_v_pt_effic_dec};
  
  vector<TH1D*> p8_1D = {consDist, consGirth, conspT, jetMult, jetGirth};
  vector<TH2D*> p8_2D = {consDist_v_pt, consGirth_v_pt, conspT_v_pt, jetMult_v_pt, jetGirth_v_pt};
  vector<TH1D*> p8_dec_1D = {consDist_dec, consGirth_dec, conspT_dec, jetMult_dec, jetGirth_dec};
  vector<TH2D*> p8_dec_2D = {consDist_v_pt_dec, consGirth_v_pt_dec, conspT_v_pt_dec, jetMult_v_pt_dec, jetGirth_v_pt_dec};
  
  //  HistsFromTreeP6P8Compare(p8_undecayed, p8_1D, p8_2D);  
  //HistsFromTreeP6P8Compare(p8_decayed, p8_dec_1D, p8_dec_2D);  
  
  HistsFromTreeP8(p8_decayed, "ResultTree", dec_hists1D, dec_hists2D, dec_hists3D);
  HistsFromTreeP8(p8_undecayed, "ResultTree", un_hists1D, un_hists2D, dec_hists3D);

  TFile *fout = new TFile( ( out + "hists_R06.root" ).c_str() ,"RECREATE");
  
  for (int i = 0; i < un_hists1D.size(); ++ i) {
    un_hists1D[i]->Write(); //dec_hists1D[i]->Write();
    //un_effic_hists1D[i]->Write(); dec_effic_hists1D[i]->Write();
  }
  for (int i = 0; i < un_hists2D.size(); ++ i) {
    un_hists2D[i]->Write(); //dec_hists2D[i]->Write();
    //un_effic_hists2D[i]->Write(); dec_effic_hists2D[i]->Write();
  }

  for (int i = 0; i < dec_hists3D.size(); ++ i) {
    un_hists3D[i]->Write(); //dec_hists3D[i]->Write();
  }
  
  for (int i = 0; i < p8_1D.size(); ++ i) {
    p8_1D[i]->Write(); p8_2D[i]->Write(); //p8_dec_1D[i]->Write(); p8_dec_2D[i]->Write();
  }
    
  fout->Close();
  
  return;
}
