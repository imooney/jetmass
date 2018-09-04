#include "RooUnfold.h"
#include <string>
#include <iostream>
#include "Plots.h"
#include "math.h"

using namespace std;

void HistsFromTree(TFile* file, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D, const bool dataflag) {
  vector<double> *Pt = 0; vector<double> *Eta = 0; vector<double> *Phi = 0;
  vector<double> *M = 0; vector<double> *zg = 0; vector<double> *rg = 0; vector<double> *mg = 0;
  vector<double> *ptg = 0; vector<double> *E = 0; vector<double> *ch_e_frac = 0; vector<double> *mcd = 0;
  double weight = 1; double n_jets = 0;

  TTree* t = (TTree*) file->Get("event");
  
  t->SetBranchAddress("n_jets",&n_jets);
  t->SetBranchAddress("Pt",&Pt); t->SetBranchAddress("Eta", &Eta); t->SetBranchAddress("Phi", &Phi); t->SetBranchAddress("E", &E);
  t->SetBranchAddress("M",&M); t->SetBranchAddress("zg", &zg); t->SetBranchAddress("rg", &rg);
  t->SetBranchAddress("mg",&mg); t->SetBranchAddress("ptg",&ptg); t->SetBranchAddress("ch_e_frac", &ch_e_frac);
  t->SetBranchAddress("mcd",&mcd);
  if (dataflag == 0) {t->SetBranchAddress("weight", &weight);}
  
  for (int i = 0; i < t->GetEntries(); ++ i) {
    t->GetEntry(i);
    //Filling 1Ds
    for (int j = 0; j < Pt->size(); ++ j) { //all vectors of doubles in the branches should have the same size
      //inclusive plots (1D)
      hists1D[0]->Fill(Pt->at(j), weight);
      hists1D[1]->Fill(Pt->at(j), weight);
      hists1D[2]->Fill(Eta->at(j), weight);
      hists1D[3]->Fill(Phi->at(j), weight);
      hists1D[4]->Fill(E->at(j), weight);
      hists1D[5]->Fill(M->at(j), weight);
      hists1D[6]->Fill(ch_e_frac->at(j), weight);
      if (zg->at(j) >= 0.1) {
	hists1D[7]->Fill(zg->at(j), weight);
	hists1D[8]->Fill(rg->at(j), weight);
	hists1D[9]->Fill(mg->at(j), weight);
	hists1D[10]->Fill(ptg->at(j), weight);
	hists1D[11]->Fill(ptg->at(j) / (double) Pt->at(j));
      }
      hists1D[12]->Fill(mcd->at(j), weight);
      
      //inclusive > 15 GeV (1D)
      if (Pt->at(j) >= 15) {
	hists1D[13]->Fill(Pt->at(j), weight);
	hists1D[14]->Fill(Eta->at(j), weight);
	hists1D[15]->Fill(Phi->at(j), weight);
	hists1D[16]->Fill(E->at(j), weight);
	hists1D[17]->Fill(M->at(j), weight);
	hists1D[18]->Fill(ch_e_frac->at(j), weight);
	if (zg->at(j) >= 0.1) {
	  hists1D[19]->Fill(zg->at(j), weight);
	  hists1D[20]->Fill(rg->at(j), weight);
	  hists1D[21]->Fill(mg->at(j), weight);
	  hists1D[22]->Fill(ptg->at(j), weight);
	  hists1D[23]->Fill(ptg->at(j) / (double) Pt->at(j));
	}
	hists1D[24]->Fill(mcd->at(j), weight);
      }
      //Filling 2Ds!
      hists2D[0]->Fill(M->at(j), Pt->at(j), weight);
      hists2D[1]->Fill(ch_e_frac->at(j), Pt->at(j), weight);
      if (zg->at(j) >= 0.1) {
	hists2D[2]->Fill(zg->at(j), Pt->at(j), weight);
	hists2D[3]->Fill(rg->at(j), Pt->at(j), weight);
	hists2D[4]->Fill(mg->at(j), Pt->at(j), weight);
	hists2D[5]->Fill(ptg->at(j), Pt->at(j), weight);
	hists2D[6]->Fill(ptg->at(j) / (double) Pt->at(j), Pt->at(j), weight);
      }
      hists2D[7]->Fill(mcd->at(j), Pt->at(j), weight);
    }
  }
  t->ResetBranchAddresses();
  return;
}

void HistsFromTreeMatched(TFile* file, std::vector<TH1D*> hists1D, std::vector<TH2D*> hists2D) {
  vector<double> *deltaPt = 0; vector<double> *deltaM = 0; vector<double> *deltaZg = 0; vector<double> *deltaRg = 0;
  vector<double> *ratioPt = 0; vector<double> *ratioM = 0; vector<double> *ratioZg = 0; vector<double> *ratioRg = 0;
  vector<double> *pyPt = 0; vector<double> *pyM = 0;
  vector<double> *pyZg = 0; vector<double> *pyRg = 0; vector<double> *pyPtg = 0; vector<double> *pyMg = 0;
  vector<double> *gePt = 0; vector<double> *geM = 0;
  vector<double> *geZg = 0; vector<double> *geRg = 0; vector<double> *gePtg = 0; vector<double> *geMg = 0;
  double weight = 1;

  TTree* t = (TTree*) file->Get("event");
  
  t->SetBranchAddress("pyPt",&pyPt); t->SetBranchAddress("pyM",&pyM); t->SetBranchAddress("pyZg", &pyZg); t->SetBranchAddress("pyRg", &pyRg);
  t->SetBranchAddress("pyPtg",&pyPtg); t->SetBranchAddress("pyMg",&pyMg);
  t->SetBranchAddress("gePt",&gePt); t->SetBranchAddress("geM",&geM); t->SetBranchAddress("geZg", &geZg); t->SetBranchAddress("geRg", &geRg);
  t->SetBranchAddress("gePtg",&gePtg); t->SetBranchAddress("geMg",&geMg);
  t->SetBranchAddress("weight", &weight);

  double count_sd_loss = 0; double count_sd = 0;
  
  for (int i = 0; i < t->GetEntries(); ++ i) {
    t->GetEntry(i);
    for(int j = 0; j < pyM->size(); ++ j) {
      //DELTAS
      hists1D[0]->Fill((gePt->at(j) - pyPt->at(j))/(double)pyPt->at(j), weight);
      hists1D[1]->Fill((geM->at(j) - pyM->at(j))/(double)pyM->at(j), weight);
      if (pyZg->at(j) >= 0.1 && geZg->at(j) < 0.1) {count_sd_loss ++;} count_sd ++;
      if (pyZg->at(j) >= 0.1 && geZg->at(j) >= 0.1) { //only fill softdrop'd deltas if we started & finished with a groomed jet with zg > 0.1
	hists1D[2]->Fill(geZg->at(j) - pyZg->at(j)/(double)pyZg->at(j), weight);
	hists1D[3]->Fill(geRg->at(j) - pyRg->at(j)/(double)pyRg->at(j), weight);
	hists1D[4]->Fill((gePtg->at(j) - pyPtg->at(j))/(double)pyPtg->at(j), weight);
	hists1D[5]->Fill((geMg->at(j) - pyMg->at(j))/(double)pyMg->at(j), weight);
      }
      //RATIOS
      hists1D[6]->Fill(gePt->at(j) / (double) pyPt->at(j), weight);
      hists1D[7]->Fill(geM->at(j) / (double) pyM->at(j), weight);
      if (pyZg->at(j) >= 0.1 && geZg->at(j) >= 0.1) { //only fill softdrop'd ratios if we started & finished with a groomed jet with zg > 0.1
	hists1D[8]->Fill(geZg->at(j) / (double) pyZg->at(j), weight);
	hists1D[9]->Fill(geRg->at(j) / (double) pyRg->at(j), weight);
	hists1D[10]->Fill(gePtg->at(j) / (double) pyPtg->at(j), weight);
	hists1D[11]->Fill(geMg->at(j) / (double) pyMg->at(j), weight);
      } 
      //PYTHIA
      //if (pyPt->at(j) > 15) { //inclusive pT plots now must be of jets with pT > 15 GeV
	hists1D[12]->Fill(pyPt->at(j), weight);
	hists1D[13]->Fill(pyM->at(j), weight);
	if (pyZg->at(j) >= 0.1) { //only fill softdrop'd pythia jet info if the zg > 0.1
	  hists1D[14]->Fill(pyZg->at(j), weight);
	  hists1D[15]->Fill(pyRg->at(j), weight);
	  hists1D[16]->Fill(pyPtg->at(j), weight);
	  hists1D[17]->Fill(pyMg->at(j), weight);
	}
	//}
      //GEANT
      //if (gePt->at(j) > 15) { //inclusive pT plots now must be of jets with pT > 15 GeV
	hists1D[18]->Fill(gePt->at(j), weight);
	hists1D[19]->Fill(geM->at(j), weight);
	if (geZg->at(j) >= 0.1) { //only fill softdrop'd geant jet info if the zg > 0.1
	  hists1D[20]->Fill(geZg->at(j), weight);
	  hists1D[21]->Fill(geRg->at(j), weight);
	  hists1D[22]->Fill(gePtg->at(j), weight);
	  hists1D[23]->Fill(geMg->at(j), weight);
	}
	//}
      hists2D[0]->Fill(pyPt->at(j),(gePt->at(j) - pyPt->at(j)) / (double) pyPt->at(j), weight);
      hists2D[1]->Fill(gePt->at(j) / (double) pyPt->at(j), pyPt->at(j), weight);
      hists2D[2]->Fill(pyPt->at(j),(geM->at(j) - pyM->at(j)) / (double) pyM->at(j), weight);
      hists2D[3]->Fill(geM->at(j) / (double) pyM->at(j), pyPt->at(j), weight);
      if (pyZg->at(j) >= 0.1) {
	hists2D[4]->Fill(pyPt->at(j),(geZg->at(j) - pyZg->at(j))/(double) pyZg->at(j), weight);
	hists2D[5]->Fill(geZg->at(j) / (double) pyZg->at(j), pyPt->at(j), weight);
	hists2D[6]->Fill(pyPt->at(j),(geRg->at(j) - pyRg->at(j))/(double) pyRg->at(j), weight);
	hists2D[7]->Fill(geRg->at(j) / (double) pyRg->at(j), pyPt->at(j), weight);
	hists2D[8]->Fill(pyPt->at(j),(gePtg->at(j) - pyPtg->at(j)) / (double) pyPtg->at(j), weight);
	hists2D[9]->Fill(gePtg->at(j) / (double) pyPtg->at(j), pyPt->at(j), weight);
	hists2D[10]->Fill(pyPt->at(j),(geMg->at(j) - pyMg->at(j)) / (double) pyMg->at(j), weight);
	hists2D[11]->Fill(geMg->at(j) / (double) pyMg->at(j), pyPt->at(j), weight);		
      }
    }
  }
  double percent_sd_loss = 100 * (double) count_sd_loss / (double) count_sd;
  
  std::cout << "Lost " << percent_sd_loss << "% of SD'd jets due to detector" << std::endl;
  
  t->ResetBranchAddresses();
  return;
}

void HistsForMCClosure(TFile* file,std::vector<TH1D*> vec) {
  vector<double> *Pt = 0; vector<double> *M = 0; vector<double> *zg = 0; vector<double> *rg = 0; vector<double> *mg = 0;
  vector<double> *ptg = 0;
  double weight; int EventID;

  TTree* t = (TTree*) file->Get("event");

  t->SetBranchAddress("Pt",&Pt); t->SetBranchAddress("M",&M); t->SetBranchAddress("zg", &zg); t->SetBranchAddress("rg", &rg);
  t->SetBranchAddress("mg",&mg); t->SetBranchAddress("ptg",&ptg);
  t->SetBranchAddress("weight", &weight); t->SetBranchAddress("EventID",&EventID);

  for (int i = 0; i < t->GetEntries(); ++ i) {
    t->GetEntry(i);
    if ((EventID % 2) == 0) { //it's an even event
      for (int j = 0; j < Pt->size(); ++ j) {
	vec[0]->Fill(Pt->at(j), weight);
	vec[2]->Fill(M->at(j), weight);
	vec[4]->Fill(zg->at(j), weight);
	vec[6]->Fill(rg->at(j), weight);
        vec[8]->Fill(ptg->at(j), weight);
	vec[10]->Fill(mg->at(j), weight);
      }
    }
    else {
      for (int j = 0; j < Pt->size(); ++ j) {
	vec[1]->Fill(Pt->at(j), weight);
	vec[3]->Fill(M->at(j), weight);
	vec[5]->Fill(zg->at(j), weight);
	vec[7]->Fill(rg->at(j), weight);
	vec[9]->Fill(ptg->at(j), weight);
	vec[11]->Fill(mg->at(j), weight);
      }
    }
    
  }
  
  t->ResetBranchAddresses();
  return;
}

/*
void HistsForMCClosure(TFile* file,std::vector<TH1D*> py, std::vector<TH1D*> ge) {
  vector<double> *pyPt = 0; vector<double> *pyM = 0;
  vector<double> *pyZg = 0; vector<double> *pyRg = 0; vector<double> *pyPtg = 0; vector<double> *pyMg = 0;
  vector<double> *gePt = 0; vector<double> *geM = 0;
  vector<double> *geZg = 0; vector<double> *geRg = 0; vector<double> *gePtg = 0; vector<double> *geMg = 0;
  double weight; int EventID;

  TTree* t = (TTree*) file->Get("event");

  t->SetBranchAddress("pyPt",&pyPt); t->SetBranchAddress("pyM",&pyM); t->SetBranchAddress("pyZg", &pyZg); t->SetBranchAddress("pyRg", &pyRg);
  t->SetBranchAddress("pyPtg",&pyPtg); t->SetBranchAddress("pyMg",&pyMg);
  t->SetBranchAddress("gePt",&gePt); t->SetBranchAddress("geM",&geM); t->SetBranchAddress("geZg", &geZg); t->SetBranchAddress("geRg", &geRg);
  t->SetBranchAddress("gePtg",&gePtg); t->SetBranchAddress("geMg",&geMg);
  t->SetBranchAddress("weight", &weight); t->SetBranchAddress("EventID", &EventID);

  for (int i = 0; i < t->GetEntries(); ++ i) {
    t->GetEntry(i);
    if ((EventID % 2) == 0) { //it's an even event
      for (int j = 0; j < pyPt->size(); ++ j) {
	py[0]->Fill(pyPt->at(j), weight);
	py[2]->Fill(pyM->at(j), weight);
	py[4]->Fill(pyZg->at(j), weight);
	py[6]->Fill(pyRg->at(j), weight);
        py[8]->Fill(pyPtg->at(j), weight);
	py[10]->Fill(pyMg->at(j), weight);
	ge[0]->Fill(gePt->at(j), weight);
	ge[2]->Fill(geM->at(j), weight);
	ge[4]->Fill(geZg->at(j), weight);
	ge[6]->Fill(geRg->at(j), weight);
        ge[8]->Fill(gePtg->at(j), weight);
	ge[10]->Fill(geMg->at(j), weight);
	
      }
    }
    else {
      for (int j = 0; j < pyPt->size(); ++ j) {
	py[1]->Fill(pyPt->at(j), weight);
	py[3]->Fill(pyM->at(j), weight);
	py[5]->Fill(pyZg->at(j), weight);
	py[7]->Fill(pyRg->at(j), weight);
	py[9]->Fill(pyPtg->at(j), weight);
	py[11]->Fill(pyMg->at(j), weight);
	ge[1]->Fill(gePt->at(j), weight);
	ge[3]->Fill(geM->at(j), weight);
	ge[5]->Fill(geZg->at(j), weight);
	ge[7]->Fill(geRg->at(j), weight);
	ge[9]->Fill(gePtg->at(j), weight);
	ge[11]->Fill(geMg->at(j), weight);
      }
    }
    
  }
  
  t->ResetBranchAddresses();
  return;
}
*/
void TreetoHist() {
  string dir = "~/jetmass/";
  string matchin = "out/matching/";
  string pyin = "out/sim/py/";
  string gein = "out/sim/ge/";
  string datain = "out/data/";
  string file = "full.root";
  string out = "~/jetmass/macros/hists/";
  string filetype = ".pdf";

  TFile* matchFile = new TFile( (dir + matchin + file).c_str(), "READ");
  TFile* pyFile = new TFile( (dir + pyin + file).c_str(), "READ");
  TFile* geFile = new TFile( (dir + gein +file).c_str(), "READ");
  TFile* dataFile = new TFile( (dir + datain + file).c_str(), "READ");
  
  //making variable bin size histogram
  const int nBins_pt = 8;
  double edges[nBins_pt + 1] = {5,10,15,20,25,30,40,60,100};
  
  //1D inclusives (including below 15 GeV pT jets)
  TH1D* pt_d_var_bin = new TH1D("pt_d_var_bin","",nBins_pt,edges);
  TH1D* pt_g_var_bin = new TH1D("pt_g_var_bin","",nBins_pt,edges);
  TH1D* pt_p_var_bin = new TH1D("pt_p_var_bin","",nBins_pt,edges);
  TH1D* pt_d = new TH1D("pt_d","",9,15,60);
  TH1D* pt_g = new TH1D("pt_g","",9,15,60);
  TH1D* pt_p = new TH1D("pt_p","",15,5,80);
  TH1D* eta_d = new TH1D("eta_d","",50,-1,1);
  TH1D* eta_g = new TH1D("eta_g","",50,-1,1);
  TH1D* eta_p = new TH1D("eta_p","",50,-1,1);
  TH1D* phi_d = new TH1D("phi_d","",50,0,2*M_PI); 
  TH1D* phi_g = new TH1D("phi_g","",50,0,2*M_PI);
  TH1D* phi_p = new TH1D("phi_p","",50,0,2*M_PI);
  TH1D* e_d = new TH1D("e_d","",100,0,100);
  TH1D* e_g = new TH1D("e_g","",100,0,100);
  TH1D* e_p = new TH1D("e_p","",100,0,100);
  TH1D* m_d = new TH1D("m_d","",20,0,10);
  TH1D* m_g = new TH1D("m_g","",20,0,10);
  TH1D* m_p = new TH1D("m_p","",20,0,10);
  TH1D* ch_e_frac_d = new TH1D("ch_e_frac_d","",10,0,1);
  TH1D* ch_e_frac_g = new TH1D("ch_e_frac_g","",10,0,1);
  TH1D* ch_e_frac_p = new TH1D("ch_e_frac_p","",10,0,1);
  TH1D* zg_d = new TH1D("zg_d","",20,0,1);
  TH1D* zg_g = new TH1D("zg_g","",20,0,1);
  TH1D* zg_p = new TH1D("zg_p","",20,0,1);
  TH1D* rg_d = new TH1D("rg_d","",20,0,1);
  TH1D* rg_g = new TH1D("rg_g","",20,0,1);
  TH1D* rg_p = new TH1D("rg_p","",20,0,1);
  TH1D* mg_d = new TH1D("mg_d","",20,0,10);
  TH1D* mg_g = new TH1D("mg_g","",20,0,10);
  TH1D* mg_p = new TH1D("mg_p","",20,0,10);
  TH1D* ptg_d = new TH1D("ptg_d","",9,15,60);
  TH1D* ptg_g = new TH1D("ptg_g","",9,15,60);
  TH1D* ptg_p = new TH1D("ptg_p","",15,5,80);
  TH1D* ratio_ptg_pt_d = new TH1D("ratio_ptg_pt_d","",21,0,2);
  TH1D* ratio_ptg_pt_g = new TH1D("ratio_ptg_pt_g","",21,0,2);
  TH1D* ratio_ptg_pt_p = new TH1D("ratio_ptg_pt_p","",21,0,2);
  TH1D* mcd_d = new TH1D("mcd_d","",20,0,5);
  TH1D* mcd_g = new TH1D("mcd_g","",20,0,5);
  TH1D* mcd_p = new TH1D("mcd_p","",20,0,5);
  
  //1D inclusives (for jets above pT = 15 GeV)
  TH1D* pt_d_pt15 = new TH1D("pt_d_pt15","",15,5,80);
  TH1D* pt_g_pt15 = new TH1D("pt_g_pt15","",15,5,80);
  TH1D* pt_p_pt15 = new TH1D("pt_p_pt15","",15,5,80);
  TH1D* eta_d_pt15 = new TH1D("eta_d_pt15","",50,-1,1);
  TH1D* eta_g_pt15 = new TH1D("eta_g_pt15","",50,-1,1);
  TH1D* eta_p_pt15 = new TH1D("eta_p_pt15","",50,-1,1);
  TH1D* phi_d_pt15 = new TH1D("phi_d_pt15","",50,0,2*M_PI); 
  TH1D* phi_g_pt15 = new TH1D("phi_g_pt15","",50,0,2*M_PI);
  TH1D* phi_p_pt15 = new TH1D("phi_p_pt15","",50,0,2*M_PI);
  TH1D* e_d_pt15 = new TH1D("e_d_pt15","",100,0,100);
  TH1D* e_g_pt15 = new TH1D("e_g_pt15","",100,0,100);
  TH1D* e_p_pt15 = new TH1D("e_p_pt15","",100,0,100);
  TH1D* m_d_pt15 = new TH1D("m_d_pt15","",20,0,10);
  TH1D* m_g_pt15 = new TH1D("m_g_pt15","",20,0,10);
  TH1D* m_p_pt15 = new TH1D("m_p_pt15","",20,0,10);
  TH1D* ch_e_frac_d_pt15 = new TH1D("ch_e_frac_d_pt15","",10,0,1);
  TH1D* ch_e_frac_g_pt15 = new TH1D("ch_e_frac_g_pt15","",10,0,1);
  TH1D* ch_e_frac_p_pt15 = new TH1D("ch_e_frac_p_pt15","",10,0,1);
  TH1D* zg_d_pt15 = new TH1D("zg_d_pt15","",20,0,1);
  TH1D* zg_g_pt15 = new TH1D("zg_g_pt15","",20,0,1);
  TH1D* zg_p_pt15 = new TH1D("zg_p_pt15","",20,0,1);
  TH1D* rg_d_pt15 = new TH1D("rg_d_pt15","",20,0,1);
  TH1D* rg_g_pt15 = new TH1D("rg_g_pt15","",20,0,1);
  TH1D* rg_p_pt15 = new TH1D("rg_p_pt15","",20,0,1);
  TH1D* mg_d_pt15 = new TH1D("mg_d_pt15","",20,0,10);
  TH1D* mg_g_pt15 = new TH1D("mg_g_pt15","",20,0,10);
  TH1D* mg_p_pt15 = new TH1D("mg_p_pt15","",20,0,10);
  TH1D* ptg_d_pt15 = new TH1D("ptg_d_pt15","",15,5,80);
  TH1D* ptg_g_pt15 = new TH1D("ptg_g_pt15","",15,5,80);
  TH1D* ptg_p_pt15 = new TH1D("ptg_p_pt15","",15,5,80);
  TH1D* ratio_ptg_pt_d_pt15 = new TH1D("ratio_ptg_pt_d_pt15","",20,0,20);
  TH1D* ratio_ptg_pt_g_pt15 = new TH1D("ratio_ptg_pt_g_pt15","",20,0,20);
  TH1D* ratio_ptg_pt_p_pt15 = new TH1D("ratio_ptg_pt_p_pt15","",20,0,20);
  TH1D* mcd_d_pt15 = new TH1D("mcd_d_pt15","",20,0,5);
  TH1D* mcd_g_pt15 = new TH1D("mcd_g_pt15","",20,0,5);
  TH1D* mcd_p_pt15 = new TH1D("mcd_p_pt15","",20,0,5);
  
  //2D hists mostly for slices in pT
  TH2D* m_v_pt_d = new TH2D("m_v_pt_d","",20,0,10,15,5,80);
  TH2D* m_v_pt_g = new TH2D("m_v_pt_g","",20,0,10,15,5,80);
  TH2D* m_v_pt_p = new TH2D("m_v_pt_p","",20,0,10,15,5,80);
  TH2D* ch_frac_v_pt_d = new TH2D("ch_frac_v_pt_d","",10,0,1,15,5,80);
  TH2D* ch_frac_v_pt_g = new TH2D("ch_frac_v_pt_g","",10,0,1,15,5,80);
  TH2D* ch_frac_v_pt_p = new TH2D("ch_frac_v_pt_p","",10,0,1,15,5,80);
  TH2D* zg_v_pt_d = new TH2D("zg_v_pt_d","",20,0,1,15,5,80);
  TH2D* zg_v_pt_g = new TH2D("zg_v_pt_g","",20,0,1,15,5,80);
  TH2D* zg_v_pt_p = new TH2D("zg_v_pt_p","",20,0,1,15,5,80);
  TH2D* rg_v_pt_d = new TH2D("rg_v_pt_d","",20,0,1,15,5,80);
  TH2D* rg_v_pt_g = new TH2D("rg_v_pt_g","",20,0,1,15,5,80);
  TH2D* rg_v_pt_p = new TH2D("rg_v_pt_p","",20,0,1,15,5,80);
  TH2D* mg_v_pt_d = new TH2D("mg_v_pt_d","",20,0,10,15,5,80);
  TH2D* mg_v_pt_g = new TH2D("mg_v_pt_g","",20,0,10,15,5,80);
  TH2D* mg_v_pt_p = new TH2D("mg_v_pt_p","",20,0,10,15,5,80);
  TH2D* ptg_v_pt_d = new TH2D("ptg_v_pt_d","",9,15,60,15,5,80);
  TH2D* ptg_v_pt_g = new TH2D("ptg_v_pt_g","",9,15,60,15,5,80);
  TH2D* ptg_v_pt_p = new TH2D("ptg_v_pt_p","",15,5,80,15,5,80);
  TH2D* ratio_ptgpt_v_pt_d = new TH2D("ratio_ptgpt_v_pt_d","",21,0,2,15,5,80);
  TH2D* ratio_ptgpt_v_pt_g = new TH2D("ratio_ptgpt_v_pt_g","",21,0,2,15,5,80);
  TH2D* ratio_ptgpt_v_pt_p = new TH2D("ratio_ptgpt_v_pt_p","",21,0,2,15,5,80);
  TH2D* mcd_v_pt_d = new TH2D("mcd_v_pt_d","",20,0,5,15,5,80);
  TH2D* mcd_v_pt_g = new TH2D("mcd_v_pt_g","",20,0,5,15,5,80);
  TH2D* mcd_v_pt_p = new TH2D("mcd_v_pt_p","",20,0,5,15,5,80);
  
  //matched hists (1D)
  TH1D* deltaPt = new TH1D("deltaPt","",220,-1,1);
  TH1D* deltaM = new TH1D("deltaM","",220,-1,1);
  TH1D* deltaZg = new TH1D("deltaZg","",220,-1,1);
  TH1D* deltaRg = new TH1D("deltaRg","",220,-1,1);
  TH1D* deltaPtg = new TH1D("deltaPtg","",220,-1,1);
  TH1D* deltaMg = new TH1D("deltaMg","",220,-1,1);
  TH1D* ratioPt = new TH1D("ratioPt","",51,0,2);
  TH1D* ratioM = new TH1D("ratioM","",51,0,2);
  TH1D* ratioZg = new TH1D("ratioZg","",51,0,2);
  TH1D* ratioRg = new TH1D("ratioRg","",51,0,2);
  TH1D* ratioPtg = new TH1D("ratioPtg","",51,0,2);
  TH1D* ratioMg = new TH1D("ratioMg","",51,0,2);
  TH1D* pyPt = new TH1D("pyPt","",15,5,80);
  TH1D* pyM = new TH1D("pyM","",20,0,10);
  TH1D* pyZg = new TH1D("pyZg","",20,0,1);
  TH1D* pyRg = new TH1D("pyRg","",20,0,1);
  TH1D* pyPtg = new TH1D("pyPtg","",15,5,80);
  TH1D* pyMg = new TH1D("pyMg","",20,0,10);
  TH1D* gePt = new TH1D("gePt","",9,15,60);
  TH1D* geM = new TH1D("geM","",20,0,10);
  TH1D* geZg = new TH1D("geZg","",20,0,1);
  TH1D* geRg = new TH1D("geRg","",20,0,1);
  TH1D* gePtg = new TH1D("gePtg","",9,15,60);
  TH1D* geMg = new TH1D("geMg","",20,0,10);
  
  const int nBins_m = 7;
  double edges_m[nBins_m + 1] = {5,10,15,20,25,30,40,60};

  //matched hists (2D)
  TH2D *deltaPtvPyPt = new TH2D("deltaPtvPyPt",";Gen. p^{jet}_{T} [GeV/c];#Delta p_{T}^{jet} (Det - Gen) / p_{T}^{gen-jet}",11,5,60,220,-1,1);
  TH2D *ratioPtvPyPt = new TH2D("ratioPtvPyPt",";p_{T}^{det-jet} / p_{T}^{gen-jet};Gen. p^{jet}_{T} [GeV/c]",51,0,2,nBins_m, edges_m);
  TH2D *deltaMvPyPt = new TH2D("deltaMvPyPt",";Gen. p^{jet}_{T} [GeV/c];#Delta M_{jet} (Det - Gen) / M^{gen}_{jet}",11,5,60,220,-1,1);
  TH2D *ratioMvPyPt = new TH2D("ratioMvPyPt",";M^{det}_{jet} / M^{gen}_{jet};Gen. p^{jet}_{T} [GeV/c]",51,0,2,nBins_m, edges_m);
  TH2D *deltaZgvPyPt = new TH2D("deltaZgvPyPt",";Gen. p^{jet}_{T} [GeV/c];#Delta z_{g} (Det - Gen)",11,5,60,220,-1,1);
  TH2D *ratioZgvPyPt = new TH2D("ratioZgvPyPt",";z_{g}^{det} / z_{g}^{gen};Gen. p^{jet}_{T} [GeV/c]",51,0,2,nBins_m, edges_m);
  TH2D *deltaRgvPyPt = new TH2D("deltaRgvPyPt",";Gen. p^{jet}_{T} [GeV/c];#Delta R_{g} (Det - Gen)",11,5,60,220,-1,1);
  TH2D *ratioRgvPyPt = new TH2D("ratioRgvPyPt",";R_{g}^{det} / R_{g}^{gen};Gen. p^{jet}_{T} [GeV/c]",51,0,2,nBins_m, edges_m);
  TH2D *deltaPtgvPyPt = new TH2D("deltaPtgvPyPt",";Gen. p^{jet}_{T} [GeV/c];#Delta p_{T,g} (Det - Gen) / p_{T,g}^{gen-jet}",11,5,60,220,-1,1);
  TH2D *ratioPtgvPyPt = new TH2D("ratioPtgvPyPt",";p_{T,g}^{det} / p_{T,g}^{gen};Gen. p^{jet}_{T} [GeV/c]",51,0,2,nBins_m, edges_m);
  TH2D *deltaMgvPyPt = new TH2D("deltaMgvPyPt",";Gen. p^{jet}_{T} [GeV/c];#Delta M_{g} (Det - Gen) / M_{g}^{gen-jet}",11,5,60,220,-1,1);
  TH2D *ratioMgvPyPt = new TH2D("ratioMgvPyPt",";M_{g}^{det} / M_{g}^{gen};Gen. p^{jet}_{T} [GeV/c]",51,0,2,nBins_m, edges_m);
  /*
  //hists to be split into evens and odds for MC Closure tests
  TH1D* pyPtEven = new TH1D("pyPtEven","",15,5,80); TH1D* pyPtOdd = new TH1D("pyPtOdd","",15,5,80);
  TH1D* pyMEven = new TH1D("pyMEven","",20,0,10); TH1D* pyMOdd = new TH1D("pyMOdd","",20,0,10);
  TH1D* pyZgEven = new TH1D("pyZgEven","",20,0,1); TH1D* pyZgOdd = new TH1D("pyZgOdd","",20,0,1);
  TH1D* pyRgEven = new TH1D("pyRgEven","",20,0,1); TH1D* pyRgOdd = new TH1D("pyRgOdd","",20,0,1);
  TH1D* pyPtgEven = new TH1D("pyPtgEven","",15,5,80); TH1D* pyPtgOdd = new TH1D("pyPtgOdd","",15,5,80);
  TH1D* pyMgEven = new TH1D("pyMgEven","",20,0,10); TH1D* pyMgOdd = new TH1D("pyMgOdd","",20,0,10);
  TH1D* gePtEven = new TH1D("gePtEven","",9,15,60); TH1D* gePtOdd = new TH1D("gePtOdd","",9,15,60);
  TH1D* geMEven = new TH1D("geMEven","",20,0,10); TH1D* geMOdd = new TH1D("geMOdd","",20,0,10);
  TH1D* geZgEven = new TH1D("geZgEven","",20,0,1); TH1D* geZgOdd = new TH1D("geZgOdd","",20,0,1);
  TH1D* geRgEven = new TH1D("geRgEven","",20,0,1); TH1D* geRgOdd = new TH1D("geRgOdd","",20,0,1);
  TH1D* gePtgEven = new TH1D("gePtgEven","",9,15,60); TH1D* gePtgOdd = new TH1D("gePtgOdd","",9,15,60);
  TH1D* geMgEven = new TH1D("geMgEven","",20,0,10); TH1D* geMgOdd = new TH1D("geMgOdd","",20,0,10);
  */
  vector<TH1D*> d_hists1D = {pt_d_var_bin,pt_d,eta_d,phi_d,e_d,m_d,ch_e_frac_d,zg_d,rg_d,mg_d,ptg_d,ratio_ptg_pt_d,mcd_d,pt_d_pt15,eta_d_pt15,phi_d_pt15,e_d_pt15,m_d_pt15,ch_e_frac_d_pt15,zg_d_pt15,rg_d_pt15,mg_d_pt15,ptg_d_pt15,ratio_ptg_pt_d_pt15,mcd_d_pt15};
  vector<TH2D*> d_hists2D = {m_v_pt_d,ch_frac_v_pt_d,zg_v_pt_d,rg_v_pt_d,mg_v_pt_d,ptg_v_pt_d,ratio_ptgpt_v_pt_d,mcd_v_pt_d};

  vector<TH1D*> g_hists1D = {pt_g_var_bin,pt_g,eta_g,phi_g,e_g,m_g,ch_e_frac_g,zg_g,rg_g,mg_g,ptg_g,ratio_ptg_pt_g,mcd_g,pt_g_pt15,eta_g_pt15,phi_g_pt15,e_g_pt15,m_g_pt15,ch_e_frac_g_pt15,zg_g_pt15,rg_g_pt15,mg_g_pt15,ptg_g_pt15,ratio_ptg_pt_g_pt15,mcd_g_pt15};
  vector<TH2D*> g_hists2D = {m_v_pt_g,ch_frac_v_pt_g,zg_v_pt_g,rg_v_pt_g,mg_v_pt_g,ptg_v_pt_g,ratio_ptgpt_v_pt_g,mcd_v_pt_g};
  
  vector<TH1D*> p_hists1D = {pt_p_var_bin,pt_p,eta_p,phi_p,e_p,m_p,ch_e_frac_p,zg_p,rg_p,mg_p,ptg_p,ratio_ptg_pt_p,mcd_p,pt_p_pt15,eta_p_pt15,phi_p_pt15,e_p_pt15,m_p_pt15,ch_e_frac_p_pt15,zg_p_pt15,rg_p_pt15,mg_p_pt15,ptg_p_pt15,ratio_ptg_pt_p_pt15,mcd_p_pt15};
  vector<TH2D*> p_hists2D = {m_v_pt_p,ch_frac_v_pt_p,zg_v_pt_p,rg_v_pt_p,mg_v_pt_p,ptg_v_pt_p,ratio_ptgpt_v_pt_p,mcd_v_pt_p};

  vector<TH1D*> m_hists1D = {deltaPt,deltaM,deltaZg,deltaRg,deltaPtg,deltaMg,ratioPt,ratioM,ratioZg,ratioRg,ratioPtg,ratioMg,pyPt,pyM,pyZg,pyRg,pyPtg,pyMg,gePt,geM,geZg,geRg,gePtg,geMg};
  vector<TH2D*> m_hists2D = {deltaPtvPyPt,ratioPtvPyPt,deltaMvPyPt,ratioMvPyPt,deltaZgvPyPt,ratioZgvPyPt,deltaRgvPyPt,ratioRgvPyPt,deltaPtgvPyPt,ratioPtgvPyPt,deltaMgvPyPt,ratioMgvPyPt};
  /*
  vector<TH1D*> c_hists1D_py = {pyPtEven,pyPtOdd,pyMEven,pyMOdd,pyZgEven,pyZgOdd,pyRgEven,pyRgOdd,pyPtgEven,pyPtgOdd,pyMgEven,pyMgOdd};
  vector<TH1D*> c_hists1D_ge = {gePtEven,gePtOdd,geMEven,geMOdd,geZgEven,geZgOdd,geRgEven,geRgOdd,gePtgEven,gePtgOdd,geMgEven,geMgOdd};
  */
  HistsFromTree(dataFile, d_hists1D, d_hists2D, 1);  
  HistsFromTree(geFile, g_hists1D, g_hists2D, 0);
  HistsFromTree(pyFile, p_hists1D, p_hists2D, 0);
  HistsFromTreeMatched(matchFile, m_hists1D, m_hists2D);
  /*
  HistsForMCClosure(pyFile,c_hists1D_py);
  HistsForMCClosure(geFile,c_hists1D_ge);
  */
  //HistsForMCClosure(matchFile, c_hists1D_py, c_hists1D_ge); //this function now uses the ~matched~ jets to create spectra
  
  TFile *fout = new TFile( ( out + "hists.root" ).c_str() ,"RECREATE");

  for (int i = 0; i < d_hists1D.size(); ++ i) {
    // TCanvas * c = new TCanvas(((string) "c_" + (string) d_hists1D[i]->GetName()).c_str(),"",800,800);
    d_hists1D[i]->Write(); g_hists1D[i]->Write(); p_hists1D[i]->Write();
  }
  for (int i = 0; i < d_hists2D.size(); ++ i) {
    //TCanvas * c = new TCanvas(((string) "c_" + (string) d_hists2D[i]->GetName()).c_str(),"",800,800);
    d_hists2D[i]->Write();
  }
  for (int i = 0; i < g_hists2D.size(); ++ i) {
    //TCanvas * c = new TCanvas(((string) "c_" + (string) g_hists2D[i]->GetName()).c_str(),"",800,800);
    g_hists2D[i]->Write();
  }
  for (int i = 0; i < p_hists2D.size(); ++ i) {
    //TCanvas * c = new TCanvas(((string) "c_" + (string) p_hists2D[i]->GetName()).c_str(),"",800,800);
    p_hists2D[i]->Write();
  }
  for (int i = 0; i < m_hists1D.size(); ++ i) {
    m_hists1D[i]->Write();
  }
  for (int i = 0; i < m_hists2D.size(); ++ i) {
    m_hists2D[i]->Write();
  }
  /*
  for (int i = 0; i < c_hists1D_py.size(); ++ i) {
    c_hists1D_py[i]->Write();
    c_hists1D_ge[i]->Write(); //should be the same size vectors of histograms (py & ge)
  }
  */
  fout->Close();
  
  return;
}
