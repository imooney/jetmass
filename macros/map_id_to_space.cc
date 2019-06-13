#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "bemc_helper.h" //thanks Nick!

#include <iostream>

using namespace std;

int main () {

  TFile *f = new TFile("~/jetmass/macros/QA_pp_pAu_w_MB_and_vpdvzcut_3.root","READ");
  
  TGraph* JP2 = (TGraph*) f->Get("htowerId_meanEt_pAuJP2");
  TGraph* BBCMB = (TGraph*) f->Get("htowerId_meanEt_pAuBBCMB");

  //  JP2->SetDirectory(0); BBCMB->SetDirectory(0); //to avoid clashing with ownership of fout later

  TH2D* spaceJP2 = new TH2D("spaceJP2",";#eta;#phi",40,-1,1,120,-M_PI,M_PI);
  TH2D* spaceBBCMB = new TH2D("spaceBBCMB",";#eta;#phi",40,-1,1,120,-M_PI,M_PI);
  
  jetreader::BemcHelper * lookup = new jetreader::BemcHelper(); 
  
  double xjp2 = 0, yjp2 = 0; double xmb = 0, ymb = 0;
  if (JP2->GetN() != BBCMB->GetN()) {cerr << "there should be the same number of towers in each case\n"; exit(1);}
  for (int i = 0; i < JP2->GetN(); ++ i) {
    if (JP2->GetPoint(i,xjp2,yjp2) == -1 || BBCMB->GetPoint(i,xmb,ymb) == -1) {cout << "invalid request for point " << i << '\n';}
    if (yjp2 < 0.4 && yjp2 > 0.0) {
      //cout << "tower: " << (int) xjp2 << " has average Et: " << yjp2 << " for JP2\n";
      spaceJP2->Fill(lookup->towerEta((int) xjp2), lookup->towerPhi((int) xjp2));
    }
    if (ymb < 0.38 && ymb > 0.0) {
      //cout << "tower: " << (int) xmb << " has average Et: " << ymb << " for BBCMB\n";
      spaceBBCMB->Fill(lookup->towerEta((int) xmb), lookup->towerPhi((int) xmb));
    }
  }

  TFile *fout = new TFile("~/jetmass/macros/lowE_towers.root","RECREATE");
  spaceJP2->SetDirectory(gDirectory);
  spaceBBCMB->SetDirectory(gDirectory);

  cout << "TEST ENDED AFTER SCANNING OVER " << BBCMB->GetN() << " TOWERS\n";
  
  //spaceJP2->Write();
  fout->Write();

  return 0;
}
