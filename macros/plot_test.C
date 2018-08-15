#include <iostream>
#include <string>
#include "Plots.h"
#include <vector>

using namespace std;

void plot_test() {
  string dir = "~/Desktop/jetmass_local/";
  string datain = "out/data/";
  string file = "tree_test.root";

  TFile* dataFile = new TFile( (dir + datain + file).c_str(), "READ");
 
  TCanvas *test = MakeCanvas("test", "0",800,800);
  
  TTree *d = (TTree*) dataFile->Get("event");
  
  vector<double> *pT = 0;
  d->SetBranchAddress("Pt", &pT);
  
  TH1D* d_pT = new TH1D("d_pT","",100,0,40);
  
  for (int i = 0; i < d->GetEntries(); ++ i) {
    d->GetEntry(i);
    for (int j = 0; j < pT->size(); ++ j) { 
      d_pT->Fill(pT->at(j));
    }
  }
  
  d_pT->Draw();

  return;
  
}
