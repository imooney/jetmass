#include "RooUnfoldResponse.h"
#include "Plots.h"

using namespace std;

void systematics () {
  gStyle->SetPalette(kPastel);
  /*
  TFile *p6in = new TFile("~/jetmass/macros/hists/hists_w_o_bin_drop.root","READ");
  TH2D *pt_res_py2D = (TH2D*) p6in->Get("deltaPtvPyPt");
  TH2D *pt_res_ge2D = (TH2D*) p6in->Get("deltaPtvGePt");
  TProfile *pt_res_py = (TProfile*) pt_res_py2D->ProfileX("pt_res_py",1,220);
  TProfile *pt_res_ge = (TProfile*) pt_res_ge2D->ProfileX("pt_res_ge",1,220);

  TH1D *p8ratiop6 = (TH1D*) p6in->Get("p8M_clone");
  
  pt_res_py->SetDirectory(0);
  pt_res_ge->SetDirectory(0);
  p8ratiop6->SetDirectory(0);
  p6in->Close();
  
  TCanvas *cmsmear = new TCanvas ("cmsmear","cmsmear",1400,1000);
  cmsmear->cd();
  Prettify1D (p8ratiop6, kBlue, kOpenCircle, 2, kBlue, "m_{j}^{gen.} [GeV/c^{2}]","PY8 / PY6", 0,10,-1,-1);
  p8ratiop6->Draw();

  TCanvas *cptsmear = new TCanvas ("cptsmear","cptsmear",1400,1000);
  cptsmear->cd();
  PrettifyTProfile (pt_res_py, kBlue, kOpenCircle, 2, kBlue,"p_{T,j} [GeV/c]", "<(p_{T}^{det}/p_{T}^{gen.}) - 1>",5,80,-1,-1);
  PrettifyTProfile (pt_res_ge, kRed, kOpenCircle, 2, kRed,"p_{T,j} [GeV/c]", "<(p_{T}^{det}/p_{T}^{gen.}) - 1>",5,80,-1,-1);
  TLegend *tres = new TLegend(0.5,0.6,0.8,0.8); tres->SetBorderSize(0);
  tres->AddEntry(pt_res_py,"Gen. level","p");
  tres->AddEntry(pt_res_ge,"Det. level","p");
  pt_res_py->Draw("same"); pt_res_ge->Draw("same"); tres->Draw("same");
  */

  TFile *fres = new TFile("~/jetmass/out/systematics/full_R04.root"/*ch.root"*/,"READ");
  RooUnfoldResponse *rnom = (RooUnfoldResponse*) fres->Get("pt_m_res_nom"); //change back1
  RooUnfoldResponse *rTS = (RooUnfoldResponse*) fres->Get("pt_m_res_TS"); //!
  RooUnfoldResponse *rTU = (RooUnfoldResponse*) fres->Get("pt_m_res_TU"); //!
  RooUnfoldResponse *rHC50 = (RooUnfoldResponse*) fres->Get("pt_m_res_HC50"); //!
  RooUnfoldResponse *rHC0 = (RooUnfoldResponse*) fres->Get("pt_m_res_HC0"); //!
  RooUnfoldResponse *rDS = (RooUnfoldResponse*) fres->Get("pt_m_res_DS"); //!
  RooUnfoldResponse *rGS = (RooUnfoldResponse*) fres->Get("pt_m_res_GS"); //!
  RooUnfoldResponse *rMS = (RooUnfoldResponse*) fres->Get("pt_m_res_MS"); //!
  
  //smearing the priors:
  //cout << endl << endl << endl << endl << "DETECTOR SMEARING" << endl << endl << endl << endl;
  //smear(rDS, pt_res_ge,"d");
  //cout << endl << endl << endl << endl << "PARTICLE SMEARING" << endl << endl << endl << endl;
  //smear(rGS, pt_res_py,"p");
  //cout << endl << endl << endl << endl << "MASS SMEARING" << endl << endl << endl << endl << endl;
  //smear(rMS, p8ratiop6);
  

  TFile *fdat = new TFile("~/jetmass/macros/hists/hists_R04.root"/*ch_hists.root"*/,"READ");
  TFile *herFile = new TFile("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved_R04.root","READ"); //m_v_pt_her
  TFile *p8File = new TFile("~/jetmass/production/macros/hists/hists_allsim_lowzgremoved_R04.root","READ"); //m_v_pt_un
  //TFile *matchFile = new TFile("~/jetmass/out/matching/full.root","READ");
  //TFile *closureFile = new TFile("~/jetmass/out/closure/full.root","READ");
  
  TH2D* m_pt_dat = (TH2D*) fdat->Get("m_v_pt_d"); TH2D* m_pt_ge = (TH2D*) fdat->Get("m_v_pt_g");
  TH2D* m_pt_py = (TH2D*) fdat->Get("m_v_pt_p"); //!
  TH2D* m_pt_her = (TH2D*) herFile->Get("mvpt_h7off"); TH2D* m_pt_p8 = (TH2D*) p8File->Get("mvpt_p8off");
  TH2D* m_pt_p8PL = (TH2D*) p8File->Get("mvpt_p8offPL"); //<--- DON'T CHANGE WHEN GOING TO GROOMED MASS!!!
  
  //TH2D* m_pt_dat_copy = (TH2D*) m_pt_dat->Clone("m_v_pt_d_copy");
  
  RooUnfoldBayes *unfold_nom = new RooUnfoldBayes(rnom, m_pt_dat, 4, false, "unfold_nom","");
  RooUnfoldBayes *unfold_IP2 = new RooUnfoldBayes(rnom, m_pt_dat, 2, false, "unfold_IP2","");
  RooUnfoldBayes *unfold_IP6 = new RooUnfoldBayes(rnom, m_pt_dat, 6, false, "unfold_IP6","");
  RooUnfoldBayes *unfold_TS = new RooUnfoldBayes(rTS, m_pt_dat, 4, false, "unfold_TS","");
  RooUnfoldBayes *unfold_TU = new RooUnfoldBayes(rTU, m_pt_dat, 4, false, "unfold_TU","");
  RooUnfoldBayes *unfold_HC50 = new RooUnfoldBayes(rHC50, m_pt_dat, 4, false, "unfold_HC50","");
  RooUnfoldBayes *unfold_HC0 = new RooUnfoldBayes(rHC0, m_pt_dat, 4, false, "unfold_HC0","");
  RooUnfoldBayes *unfold_DS = new RooUnfoldBayes(rDS, m_pt_dat, 4, false, "unfold_DS","");
  RooUnfoldBayes *unfold_GS = new RooUnfoldBayes(rGS, m_pt_dat, 4, false, "unfold_GS","");
  RooUnfoldBayes *unfold_MS = new RooUnfoldBayes(rMS, m_pt_dat, 4, false, "unfold_MS","");
  
  TH2D *reco_nom = (TH2D*) unfold_nom->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_IP2 = (TH2D*) unfold_IP2->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_IP6 = (TH2D*) unfold_IP6->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_TS = (TH2D*) unfold_TS->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_TU = (TH2D*) unfold_TU->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_HC50 = (TH2D*) unfold_HC50->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_HC0 = (TH2D*) unfold_HC0->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_DS = (TH2D*) unfold_DS->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_GS = (TH2D*) unfold_GS->Hreco((RooUnfold::ErrorTreatment) 3);
  TH2D *reco_MS = (TH2D*) unfold_MS->Hreco((RooUnfold::ErrorTreatment) 3);

  
  const int nBins = 3;
  TAxis* reco_axis = reco_nom->GetYaxis(); TAxis* det_axis = m_pt_dat->GetYaxis();  
  double ranges[nBins + 1] = {(double) reco_axis->FindBin(20), (double) reco_axis->FindBin(25), (double) reco_axis->FindBin(30), (double) reco_axis->FindBin(40)};//3,4,5,6,8,12};
  double ranges_d[nBins + 1] = {(double) det_axis->FindBin(20), (double) det_axis->FindBin(25), (double) det_axis->FindBin(30), (double) det_axis->FindBin(40)};//1,2,3,4,6,10};
  string pts[nBins + 1] = {"20","25","30","40"};
  
  vector<TH1D*> reco_noms = Projection2D(reco_nom,nBins,ranges,"x");
  vector<TH1D*> reco_IP2s = Projection2D(reco_IP2,nBins,ranges,"x");
  vector<TH1D*> reco_IP6s = Projection2D(reco_IP6,nBins,ranges,"x");
  vector<TH1D*> reco_TSs = Projection2D(reco_TS,nBins,ranges,"x");
  vector<TH1D*> reco_TUs = Projection2D(reco_TU,nBins,ranges,"x");
  vector<TH1D*> reco_HC50s = Projection2D(reco_HC50,nBins,ranges,"x");
  vector<TH1D*> reco_HC0s = Projection2D(reco_HC0,nBins,ranges,"x");
  vector<TH1D*> reco_DSs = Projection2D(reco_DS,nBins,ranges,"x");
  vector<TH1D*> reco_GSs = Projection2D(reco_GS,nBins,ranges,"x");
  vector<TH1D*> reco_MSs = Projection2D(reco_MS,nBins,ranges,"x");

  vector<TH1D*> reco_noms_copy;
  vector<TH1D*> reco_IP2s_copy;
  vector<TH1D*> reco_IP6s_copy;
  vector<TH1D*> reco_TSs_copy;
  vector<TH1D*> reco_TUs_copy;
  vector<TH1D*> reco_HC50s_copy;
  vector<TH1D*> reco_HC0s_copy;
  vector<TH1D*> reco_DSs_copy;
  vector<TH1D*> reco_GSs_copy;
  vector<TH1D*> reco_MSs_copy;

  
  for (int i = 0; i < nBins; ++ i) {
    reco_noms_copy.push_back((TH1D*) reco_noms[i]->Clone(("nom_"+to_string(i)).c_str()));
    reco_IP2s_copy.push_back((TH1D*) reco_IP2s[i]->Clone(("IP2_"+to_string(i)).c_str()));
    reco_IP6s_copy.push_back((TH1D*) reco_IP6s[i]->Clone(("IP6_"+to_string(i)).c_str()));
    reco_TSs_copy.push_back((TH1D*) reco_TSs[i]->Clone(("TS_"+to_string(i)).c_str()));
    reco_TUs_copy.push_back((TH1D*) reco_TUs[i]->Clone(("TU_"+to_string(i)).c_str()));
    reco_HC50s_copy.push_back((TH1D*) reco_HC50s[i]->Clone(("HC50_"+to_string(i)).c_str()));
    reco_HC0s_copy.push_back((TH1D*) reco_HC0s[i]->Clone(("HC0_"+to_string(i)).c_str()));
    reco_DSs_copy.push_back((TH1D*) reco_DSs[i]->Clone(("DS_"+to_string(i)).c_str()));
    reco_GSs_copy.push_back((TH1D*) reco_GSs[i]->Clone(("GS_"+to_string(i)).c_str()));
    reco_MSs_copy.push_back((TH1D*) reco_MSs[i]->Clone(("MS_"+to_string(i)).c_str()));

  }

  for (int i = 0; i < nBins; ++ i) {
    reco_noms[i]->Scale(1/(double)reco_noms[i]->Integral());
    reco_IP2s[i]->Scale(1/(double)reco_IP2s[i]->Integral());
    reco_IP6s[i]->Scale(1/(double)reco_IP6s[i]->Integral());
    reco_TSs[i]->Scale(1/(double)reco_TSs[i]->Integral());
    reco_TUs[i]->Scale(1/(double)reco_TUs[i]->Integral());
    reco_HC50s[i]->Scale(1/(double)reco_HC50s[i]->Integral());
    reco_HC0s[i]->Scale(1/(double)reco_HC0s[i]->Integral());
    reco_DSs[i]->Scale(1/(double)reco_DSs[i]->Integral());
    reco_GSs[i]->Scale(1/(double)reco_GSs[i]->Integral());
    reco_MSs[i]->Scale(1/(double)reco_MSs[i]->Integral());
    
    reco_IP2s[i]->Divide(reco_noms[i]);
    reco_IP6s[i]->Divide(reco_noms[i]);
    reco_TSs[i]->Divide(reco_noms[i]);
    reco_TUs[i]->Divide(reco_noms[i]);
    reco_HC50s[i]->Divide(reco_noms[i]);
    reco_HC0s[i]->Divide(reco_noms[i]);
    reco_DSs[i]->Divide(reco_noms[i]);
    reco_GSs[i]->Divide(reco_noms[i]);
    reco_MSs[i]->Divide(reco_noms[i]);
  }
  
  for (int i = 0; i < nBins; ++ i) {
    Double_t stats[5] = {0,0,0,0,0};
    reco_IP2s[i]->PutStats(stats);
    reco_IP6s[i]->PutStats(stats);
    reco_TSs[i]->PutStats(stats);
    reco_TUs[i]->PutStats(stats);
    reco_HC50s[i]->PutStats(stats);
    reco_HC0s[i]->PutStats(stats);
    reco_DSs[i]->PutStats(stats);
    reco_GSs[i]->PutStats(stats);
    reco_MSs[i]->PutStats(stats);

    reco_IP2s[i]->Sumw2(0);
    reco_IP6s[i]->Sumw2(0);
    reco_TSs[i]->Sumw2(0);
    reco_TUs[i]->Sumw2(0);
    reco_HC50s[i]->Sumw2(0);
    reco_HC0s[i]->Sumw2(0);
    reco_DSs[i]->Sumw2(0);
    reco_GSs[i]->Sumw2(0);
    reco_MSs[i]->Sumw2(0);

  
    for (int j = 1; j <= 14; ++ j) {
      //turning the ratio into a percentage.
      reco_IP2s[i]->SetBinContent(j,fabs(reco_IP2s[i]->GetBinContent(j) - 1));
      reco_IP6s[i]->SetBinContent(j,fabs(reco_IP6s[i]->GetBinContent(j) - 1));
      reco_TSs[i]->SetBinContent(j,fabs(reco_TSs[i]->GetBinContent(j) - 1));
      reco_TUs[i]->SetBinContent(j,fabs(reco_TUs[i]->GetBinContent(j) - 1));
      reco_HC50s[i]->SetBinContent(j,fabs(reco_HC50s[i]->GetBinContent(j) - 1));
      reco_HC0s[i]->SetBinContent(j,fabs(reco_HC0s[i]->GetBinContent(j) - 1));
      reco_DSs[i]->SetBinContent(j,fabs(reco_DSs[i]->GetBinContent(j) - 1));
      reco_GSs[i]->SetBinContent(j,fabs(reco_GSs[i]->GetBinContent(j) - 1));
      reco_MSs[i]->SetBinContent(j,fabs(reco_MSs[i]->GetBinContent(j) - 1));
    }
    
  }
  
  for (int i = 0; i < nBins; ++ i) {
    Prettify1DwLineStyle(reco_IP2s[i],2, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1); 
    reco_IP2s[i]->SetFillColor(2); reco_IP2s[i]->SetFillStyle(3305);
    Prettify1DwLineStyle(reco_IP6s[i],3, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_IP6s[i]->SetFillColor(3); reco_IP6s[i]->SetFillStyle(3395);
    Prettify1DwLineStyle(reco_TSs[i],4, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_TSs[i]->SetFillColor(4); reco_TSs[i]->SetFillStyle(3490);
    Prettify1DwLineStyle(reco_TUs[i],5, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_TUs[i]->SetFillColor(5); reco_TUs[i]->SetFillStyle(3436);
    Prettify1DwLineStyle(reco_HC50s[i],6, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_HC50s[i]->SetFillColor(6); reco_HC50s[i]->SetFillStyle(3335);
    Prettify1DwLineStyle(reco_HC0s[i],7, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1); 
    reco_HC0s[i]->SetFillColor(7); reco_HC0s[i]->SetFillStyle(3353);
    Prettify1DwLineStyle(reco_DSs[i],8, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1); 
    reco_DSs[i]->SetFillColor(8); reco_DSs[i]->SetFillStyle(3944);
    Prettify1DwLineStyle(reco_GSs[i],9, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_GSs[i]->SetFillColor(9); reco_GSs[i]->SetFillStyle(3544);
    Prettify1DwLineStyle(reco_MSs[i],11, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    reco_MSs[i]->SetFillColor(11); reco_MSs[i]->SetFillStyle(3690);

  }
  
  TLegend *l1 = new TLegend(0.1,0.5,0.28,0.8); l1->SetBorderSize(0);
  l1->AddEntry(reco_IP2s[0],"IP2","f");
  l1->AddEntry(reco_IP6s[0],"IP6","f");
  l1->AddEntry(reco_TSs[0],"TS","f");
  l1->AddEntry(reco_TUs[0],"TU","f");

  TLegend *l2 = new TLegend(0.1,0.5,0.28,0.8); l2->SetBorderSize(0);
  l2->AddEntry(reco_HC50s[0],"HC50","f");
  //  l2->AddEntry(reco_HC0s[0],"HC0","f");
  l2->AddEntry(reco_DSs[0],"DS","f");
  l2->AddEntry(reco_GSs[0],"GS","f");
  l2->AddEntry(reco_MSs[0],"MS","f");
 

  TLatex *slice = new TLatex();
  
  TCanvas *csys = MakeCanvas ("csys","0",1200,500);
  DivideCanvas(csys,"0", 3,1);
  
  TH1D* hdummy = new TH1D("hdummy",";M [GeV/c^{2}];relative uncertainty",14,0,14);
  hdummy->GetYaxis()->SetRangeUser(0,1);
    
  TLatex *p = new TLatex();
  p->SetTextAlign(11);
  p->SetTextSize(0.07);
  
  /*csys->cd(1); hdummy->Draw();*/ TLatex *t = new TLatex();
  for (int i = 0; i < nBins; ++ i) {
    csys->cd(i+1); 
    //reco_HC0s[i]->Draw("lf2same   ");
    reco_HC50s[i]->Draw("lf2same   "); reco_TSs[i]->Draw("lf2same   "); reco_TUs[i]->Draw("lf2same   ");
    reco_DSs[i]->Draw("lf2same   "); reco_GSs[i]->Draw("lf2same   "); reco_IP2s[i]->Draw("lf2same   "); reco_IP6s[i]->Draw("lf2same   "); reco_MSs[i]->Draw("lf2same   ");
    if (i == 0) {
      p->DrawLatexNDC(0.2,0.65, "pp 200 GeV run12 JP2");
      p->DrawLatexNDC(0.2,0.55, "anti-k_{T}, R = 0.6");
      p->DrawLatexNDC(0.2,0.45, "Ch+Ne jets, |#eta| < 0.4");
    }
    if (i==1) {l1->Draw("same");} if (i==2) {l2->Draw("same");} slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }
    
  //  csys->SaveAs("~/jetmass/plots/systematics/systematics_R06.pdf");


  

  //taking maximum envelopes!
  TH2D* env_HC = new TH2D("env_HC","",14,0,14,15,5,80);
  vector<TH1D*> env_HCs = Projection2D(env_HC,nBins,ranges,"x");
  TH2D* env_un = new TH2D("env_un","",14,0,14,15,5,80);
  vector<TH1D*> env_uns = Projection2D(env_un,nBins,ranges,"x");
  TH2D* net = new TH2D("net","",14,0,14,15,5,80);
  vector<TH1D*> nets = Projection2D(net,nBins,ranges,"x");
  TH2D* stat = new TH2D("stat","",14,0,14,15,5,80);
  vector<TH1D*> stats = Projection2D(stat,nBins,ranges,"x");
  
  const int nBins_m = 14;

  vector<vector< double> > syst_errs2D;
  
  for (int i = 0; i < nBins; ++ i) {
    vector<double> syst_errs1D; 
    reco_noms_copy[i]->Scale(1/(double)reco_noms_copy[i]->Integral());
    
    for (int j = 1; j < nBins_m + 1; ++ j) {
      //hadronic correction envelope - using an ordered set here to automatically get the largest value. 
      double hcs [/*2*/1] = {/*reco_HC0s[i]->GetBinContent(j),*/ reco_HC50s[i]->GetBinContent(j)};
      set<double> hc_sort (hcs, hcs+1);
      set<double>::iterator hc = hc_sort.end(); hc --;
      double hc_envelope = *hc;
      env_HCs[i]->SetBinContent(j, hc_envelope);
      //unfolding envelope
      double uns [5] = {reco_DSs[i]->GetBinContent(j), reco_GSs[i]->GetBinContent(j), reco_MSs[i]->GetBinContent(j), reco_IP2s[i]->GetBinContent(j), reco_IP6s[i]->GetBinContent(j)};
      set<double> un_sort (uns, uns+5);
      set<double>::iterator un = un_sort.end(); un --;
      double un_envelope = *un;
      env_uns[i]->SetBinContent(j, un_envelope);
      //total uncertainty = TU + TS + un envelope + hc envelope
      double square = pow(hc_envelope,2) + pow(un_envelope,2) + pow(reco_TUs[i]->GetBinContent(j),2) + pow(reco_TSs[i]->GetBinContent(j),2);
      nets[i]->SetBinContent(j,sqrt(square));
      stats[i]->SetBinContent(j,reco_noms_copy[i]->GetBinError(j));
      syst_errs1D.push_back(nets[i]->GetBinContent(j));
    }
    
    syst_errs2D.push_back(syst_errs1D);
  }

  for (int i = 0; i < nBins; ++ i) {
    Prettify1DwLineStyle(env_HCs[i],7, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    env_HCs[i]->SetFillColor(7); env_HCs[i]->SetFillStyle(3353);
    Prettify1DwLineStyle(env_uns[i],2, kSolid, 2,"M [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    env_uns[i]->SetFillColor(2); env_uns[i]->SetFillStyle(3305);
    Prettify1DwLineStyle(nets[i], kBlack, kSolid, 2, "M [GeV/c^{2}]","relative uncertainty",0,14,0,1);
    Prettify1DwLineStyle(stats[i], kBlack, kDashed, 2, "M [GeV/c^{2}]","relative uncertainty",0,14,0,1);
  }
  
  TLegend *tenvs = new TLegend(0.1,0.4,0.4,0.75); tenvs->SetBorderSize(0);
  tenvs->AddEntry(env_HCs[0],"Hadronic correction","f");
  tenvs->AddEntry(reco_TSs[0],"Tower scale","f");
  tenvs->AddEntry(reco_TUs[0],"Tracking","f");
  tenvs->AddEntry(env_uns[0],"Unfolding","f");
  tenvs->AddEntry(nets[0],"Total systematic uncertainty","l");
  //tenvs->AddEntry(stats[0],"Total statistical uncertainty","l");
  
  TCanvas *cenvs = new TCanvas("cenvs","cenvs",1200,500);
  DivideCanvas(cenvs,"0",3,1);
  
  TH1D* hdummyenvs = new TH1D("hdummyenvs",";;relative uncertainty",1,0,14);
  hdummyenvs->GetYaxis()->SetRangeUser(0,1);

  for (int i = 0; i < nBins; ++ i) {
    cenvs->cd(i+1); 
    env_HCs[i]->Draw("lf2same   "); reco_TSs[i]->Draw("lf2same   "); reco_TUs[i]->Draw("lf2same   "); env_uns[i]->Draw("lf2same   "); nets[i]->Draw("lf2same"); /*stats[i]->Draw("lf2same");*/ slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    if (i == 0) {
      p->DrawLatexNDC(0.2,0.65, "pp 200 GeV run12 JP2");
      p->DrawLatexNDC(0.2,0.55, "anti-k_{T}, R = 0.6");
      p->DrawLatexNDC(0.2,0.45, "Ch+Ne jets, |#eta| < 0.4");
    }
    if (i == 1) {tenvs->Draw("same");}
  }
  
  for (int i = 1; i <= nets[1]->GetNbinsX(); ++ i) {
    if (i == 2 || i == 5 || i == 8) {
      cout << 100*env_HCs[1]->GetBinContent(i) << "% " << 100*reco_TSs[1]->GetBinContent(i) << "% " << 100*reco_TUs[1]->GetBinContent(i) << "% " << 100*env_uns[1]->GetBinContent(i) << "% " << 100*nets[1]->GetBinContent(i) << endl;
    }
  }
  
  // cenvs->SaveAs("~/jetmass/plots/systematics/systematic_envelopes_R06.pdf");
  
  //unfolded result with systematic errors!
  vector<TH1D*> w_systs;
  for (int i = 0; i < nBins; ++ i) { 
    w_systs.push_back((TH1D*) reco_noms_copy[i]->Clone(("w_systs_"+to_string(i)).c_str()));
    for (int j = 1; j < nBins_m + 1; ++ j) {
      w_systs[i]->SetBinError(j, syst_errs2D[i][j-1]*w_systs[i]->GetBinContent(j));
    }
  }

  vector<TH1D*> dats = Projection2D(m_pt_dat,nBins,ranges_d,"x");
  vector<TH1D*> ges = Projection2D(m_pt_ge,nBins,ranges_d,"x");
  vector<TH1D*> pys = Projection2D(m_pt_py,nBins,ranges,"x");
  vector<TH1D*> p8s = Projection2D(m_pt_p8,nBins,ranges,"x");
  vector<TH1D*> hers = Projection2D(m_pt_her,nBins,ranges,"x");
  vector<TH1D*> p8PLs = Projection2D(m_pt_p8PL,nBins,ranges,"x");
  
  for (int i = 0; i < nBins; ++ i) {
    //temp
    //CONFUSED: seems like "setlimits" changes the actual values of the axis. Mean of the data seems to shift somehow. Don't use until you understand better.
    /*
    w_systs[i]->GetXaxis()->SetLimits(0,14);
    pys[i]->GetXaxis()->SetLimits(0,14);
    hers[i]->GetXaxis()->SetLimits(0,14);
    p8s[i]->GetXaxis()->SetLimits(0,14);
    reco_noms_copy[i]->GetXaxis()->SetLimits(0,14);
    */
    Prettify1D(reco_noms_copy[i],kRed,kFullStar,4,kRed,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.5);
    Prettify1DwLineStyle(pys[i], kBlue, kSolid,2, "M [GeV/c^{2}]", "1/N dN/dM",0,14,0,0.5);
    //Prettify1D(pys[i], kBlue, kFullCircle,3, kBlue, "M [GeV/c^{2}]", "1/N dN/dM",0,14,0,0.5);
    Prettify1D(ges[i], kRed, kOpenCircle, 3/*3*/, kRed, "M [GeV/c^{2}]", "1/N dN/dM",0,14,0,0.5);
    Prettify1D(dats[i], kBlack, kOpenStar, 4/*4*/, kBlack, "M [GeV/c^{2}]", "1/N dN/dM",0,14,0,0.5);
    Prettify1D(w_systs[i],kRed,kFullStar,0,kRed,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.5);
    w_systs[i]->SetFillColor(kRed - 10); w_systs[i]->SetFillStyle(/*3352*/1001);
    Prettify1DwLineStyle(p8s[i],kBlack,kSolid,2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.5);
    Prettify1DwLineStyle(hers[i],kMagenta,kSolid,2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.5);
    Prettify1DwLineStyle(p8PLs[i],kBlack,kDashed,2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.5);
  }
   
  //scaling errors
  /*
  for (int i = 0; i < nBins; ++ i) {
    cout << "i = " << i << endl << endl;
    for (int j = 1; j <= reco_noms_copy[i]->GetNbinsX(); ++ j) {
      double scaling = -1;
      if (i == 0) {scaling = 1.122;} //these numbers are calculated using the bin content of the ratio of gen. matched spectrum to gen. inclusive (unmatched).
      if (i == 1) {scaling = 1.082;}
      if (i == 2) {scaling = 1.062;}
      double binerror = reco_noms_copy[i]->GetBinError(j);
      reco_noms_copy[i]->SetBinError(j,(double) binerror*scaling);
      cout << "bin j error: " << binerror << " " << reco_noms_copy[i]->GetBinError(j) << endl;
      cout << "corresponding raw data error: " << dats[i]->GetBinError(j) << endl;
    }
  }
  */
  TCanvas *cws = new TCanvas("cws","cws",1200,500);
  DivideCanvas(cws,"0",3,1);

  TLegend *twsysts = new TLegend(0.45,0.55,0.7,0.7); twsysts->SetBorderSize(0);
  //  twsysts->AddEntry(hers[0],"Herwig7","l");
  //twsysts->AddEntry(p8s[0],"Pythia8","l");
  //twsysts->AddEntry(p8PLs[0],"Pythia8 Parton Jets","l");
  TLegend *twsysts2 = new TLegend(0.5,0.55,0.75,0.7); twsysts2->SetBorderSize(0);
  twsysts2->AddEntry(pys[0],"Pythia6","l");
  twsysts2->AddEntry(ges[0],"Pythia6+Geant","p");
  twsysts2->AddEntry(dats[0],"Raw data","p");
  TH1D* for_legend = (TH1D*) w_systs[0]->Clone("for_legend"); for_legend->SetMarkerSize(2);
  //twsysts2->AddEntry(for_legend,"Unfolded data","pf");
  //twsysts->AddEntry(stats[0],"Unfold stat. err.", "l");
  
  TLatex *tpost = new TLatex(); tpost->SetTextColor(kRed);
  
  TLatex *ttitle = new TLatex(); ttitle->SetTextAlign(11); ttitle->SetTextSize(0.05);
  t->SetTextAlign(11);
  t->SetTextSize(0.07);
  
  TH1D* hdummycws = new TH1D("hdummycws",";;1/N dN/dM",1,0,14);
  hdummycws->GetYaxis()->SetRangeUser(0,0.4);
  hdummycws->GetYaxis()->SetTitleSize(0.06);

  //cws->cd(1); hdummycws->Draw(); t = PanelTitle();
  for (int i = 0; i < nBins; ++ i) {

    cws->cd(i+1); /*w_systs[i]->Draw("E3same9");*/ pys[i]->Draw("Csame9"); /*stats[i]->Draw("lf2same");*/ /*hers[i]->Draw("Csame9"); p8s[i]->Draw("Csame9"); p8PLs[i]->Draw("Csame9");*/ ges[i]->Draw("same"); dats[i]->Draw("same"); /*reco_noms_copy[i]->Draw("same9");*/ slice->DrawLatexNDC(0.5,0.77,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
    //w_systs[i]->GetXaxis()->SetTitleSize(0.08); w_systs[i]->GetYaxis()->SetTitleSize(0.08);
    if (i == 0) {ttitle->DrawLatex(4.2,0.46, "pp 200 GeV run12 JP2");ttitle->DrawLatex(4.2,0.43, "anti-k_{T}, R = 0.6");ttitle->DrawLatex(4.2,0.4, "Ch+Ne jets, |#eta| < 0.4");}
    if (i == 1) { /*tpost->DrawLatexNDC(0.15,0.87,"PREVIEW - Work in Progress");*/twsysts->Draw("same9");}
    if (i == 2) { twsysts2->Draw("same9");}
    
  }

  //  cws->SaveAs("~/jetmass/plots/systematics/unfolded_w_systs_etc_and_PL_R06.pdf");
 
  
  //SPECTRA!
  /*
  for (int i = 0; i < nBins; ++ i) {
    reco_noms_copy[i]->Sumw2(0);
    reco_IP2s_copy[i]->Sumw2(0);
    reco_IP6s_copy[i]->Sumw2(0);
    reco_TSs_copy[i]->Sumw2(0);
    reco_TUs_copy[i]->Sumw2(0);
    reco_HC50s_copy[i]->Sumw2(0);
    reco_HC0s_copy[i]->Sumw2(0);
    reco_DSs_copy[i]->Sumw2(0);
    reco_GSs_copy[i]->Sumw2(0);
    reco_MSs_copy[i]->Sumw2(0);


    cout << endl << "i: " << i << endl;
    cout << reco_noms_copy[i]->Integral() << endl;
    cout << reco_IP2s_copy[i]->Integral() << endl;
    cout << reco_IP6s_copy[i]->Integral() << endl;
    cout << reco_TSs_copy[i]->Integral() << endl;
    cout << reco_TUs_copy[i]->Integral() << endl;
    cout << reco_HC50s_copy[i]->Integral() << endl;
    cout << reco_HC0s_copy[i]->Integral() << endl;
    cout << reco_DSs_copy[i]->Integral() << endl;
    cout << reco_GSs_copy[i]->Integral() << endl;
    cout << reco_MSs_copy[i]->Integral() << endl;


    Prettify1DwLineStyle(reco_noms_copy[i],1, kSolid, 2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.6); 
    //reco_noms_copy[i]->SetFillColor(1); reco_noms_copy[i]->SetFillStyle(3244);
    Prettify1DwLineStyle(reco_IP2s_copy[i],2, kSolid, 2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.6); 
    reco_IP2s_copy[i]->SetFillColor(2); reco_IP2s_copy[i]->SetFillStyle(3305);
    Prettify1DwLineStyle(reco_IP6s_copy[i],3, kSolid, 2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.6);
    reco_IP6s_copy[i]->SetFillColor(3); reco_IP6s_copy[i]->SetFillStyle(3395);
    Prettify1DwLineStyle(reco_TSs_copy[i],4, kSolid, 2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.6);
    reco_TSs_copy[i]->SetFillColor(4); reco_TSs_copy[i]->SetFillStyle(3490);
    Prettify1DwLineStyle(reco_TUs_copy[i],5, kSolid, 2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.6);
    reco_TUs_copy[i]->SetFillColor(5); reco_TUs_copy[i]->SetFillStyle(3436);
    Prettify1DwLineStyle(reco_HC50s_copy[i],6, kSolid, 2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.6);
    reco_HC50s_copy[i]->SetFillColor(6); reco_HC50s_copy[i]->SetFillStyle(3335);
    Prettify1DwLineStyle(reco_HC0s_copy[i],7, kSolid, 2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.6); 
    reco_HC0s_copy[i]->SetFillColor(7); reco_HC0s_copy[i]->SetFillStyle(3353);
    Prettify1DwLineStyle(reco_DSs_copy[i],8, kSolid, 2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.6); 
    reco_DSs_copy[i]->SetFillColor(8); reco_DSs_copy[i]->SetFillStyle(3944);
    Prettify1DwLineStyle(reco_GSs_copy[i],9, kSolid, 2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.6);
    reco_GSs_copy[i]->SetFillColor(9); reco_GSs_copy[i]->SetFillStyle(3544);
    Prettify1DwLineStyle(reco_MSs_copy[i],10, kSolid, 2,"M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.6);
    reco_MSs_copy[i]->SetFillColor(10); reco_MSs_copy[i]->SetFillStyle(3690);
  }

  TLegend *specl1 = new TLegend(0.7,0.5,0.88,0.8); specl1->SetBorderSize(0);
  specl1->AddEntry(reco_IP2s_copy[0],"IP2","f");
  specl1->AddEntry(reco_IP6s_copy[0],"IP6","f");
  specl1->AddEntry(reco_TSs_copy[0],"TS","f");
  specl1->AddEntry(reco_TUs_copy[0],"TU","f");

  TLegend *specl2 = new TLegend(0.7,0.5,0.88,0.8); specl2->SetBorderSize(0);
  specl2->AddEntry(reco_HC50s_copy[0],"HC50","f");
  specl2->AddEntry(reco_HC0s_copy[0],"HC0","f");
  specl2->AddEntry(reco_DSs_copy[0],"DS","f");
  specl2->AddEntry(reco_GSs_copy[0],"GS","f");
  specl2->AddEntry(reco_MSs_copy[0],"MS","f");


  TLegend *specl3 = new TLegend(0.15,0.65,0.35,0.75); specl3->SetBorderSize(0);
  specl3->AddEntry(reco_noms_copy[0],"nominal","l");
  
  TLatex *slice_spec = new TLatex();
  
  TCanvas *cspec = MakeCanvas ("cspec","0",1400,1000);
  DivideCanvas(cspec,"0", 3,2);
  
  TH1D* hdummy_spec = new TH1D("hdummy_spec",";M [GeV/c^{2}];1/N dN/dM",10,0,10);
  hdummy_spec->GetYaxis()->SetRangeUser(0,0.6);
  
  cspec->cd(1); hdummy_spec->Draw(); TLatex *tspec = new TLatex(); tspec = PanelTitle();
  for (int i = 0; i < nBins; ++ i) {
    cspec->cd(i+2); 
    reco_HC0s_copy[i]->Draw("lf2same   ");
    reco_HC50s_copy[i]->Draw("lf2same   "); reco_TSs_copy[i]->Draw("lf2same   "); reco_TUs_copy[i]->Draw("lf2same   ");
    reco_DSs_copy[i]->Draw("lf2same   "); reco_GSs_copy[i]->Draw("lf2same   "); reco_IP2s_copy[i]->Draw("lf2same   "); reco_IP6s_copy[i]->Draw("lf2same   "); 
    reco_noms_copy[i]->Draw("lsame"); reco_MSs_copy[i]->Draw("lf2same   ");
    if (i == 0) {specl1->Draw("same");} if (i==1) {specl2->Draw("same");} if (i==2) {specl3->Draw("same");} slice_spec->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str());
  }

  //cspec->SaveAs("~/jetmass/plots/systematics/spectra.pdf");
  */
  
  /*
  //QA!!!
  
  TFile *f = new TFile("~/jetmass/out/systematics/full_w_o_bin_drop.root","READ");
  
  vector<double> *pt_nom = 0; vector<double> *M_nom = 0;
  vector<double> *NEF_nom = 0; vector<double> *Ntracks_nom = 0; vector<double> *Ntows_nom = 0;
  vector<double> *pt_TU = 0; vector<double> *M_TU = 0;
  vector<double> *NEF_TU = 0; vector<double> *Ntracks_TU = 0; vector<double> *Ntows_TU = 0;
  double wt_nom = -9999; double wt_TU = -9999;

  TH2D* Ntracks_nom_v_pt = new TH2D("Ntracks_nom_v_pt","",30,0,30,9,15,60);
  TH2D* Ntows_nom_v_pt = new TH2D("Ntows_nom_v_pt","",30,0,30,9,15,60);
  TH2D* M_nom_v_pt = new TH2D("m_nom_v_pt","",10,0,14,9,15,60);
  TH2D* NEF_nom_v_pt = new TH2D("NEF_nom_v_pt","",10,0,1,9,15,60);
  TH2D* Ntracks_TU_v_pt = new TH2D("Ntracks_TU_v_pt","",30,0,30,9,15,60);
  TH2D* Ntows_TU_v_pt = new TH2D("Ntows_TU_v_pt","",30,0,30,9,15,60);
  TH2D* M_TU_v_pt = new TH2D("m_TU_v_pt","",10,0,14,9,15,60);
  TH2D* NEF_TU_v_pt = new TH2D("NEF_TU_v_pt","",10,0,1,9,15,60);

  TTree *nom = (TTree*) f->Get("nom");
  TTree *TU = (TTree*) f->Get("TU");

  nom->SetBranchAddress("pt",&pt_nom); nom->SetBranchAddress("M",&M_nom);
  nom->SetBranchAddress("NEF",&NEF_nom); nom->SetBranchAddress("Ntracks",&Ntracks_nom); nom->SetBranchAddress("Ntows",&Ntows_nom);
  nom->SetBranchAddress("mc_weight",&wt_nom);
  TU->SetBranchAddress("pt",&pt_TU); TU->SetBranchAddress("M",&M_TU);
  TU->SetBranchAddress("NEF",&NEF_TU); TU->SetBranchAddress("Ntracks",&Ntracks_TU); TU->SetBranchAddress("Ntows",&Ntows_TU);
  TU->SetBranchAddress("mc_weight",&wt_TU);

  for (int i = 0; i < nom->GetEntries(); ++ i) {
    nom->GetEntry(i);
    for (int j = 0; j < pt_nom->size(); ++ j) { //all vectors of doubles in the branches should have the same size
      Ntracks_nom_v_pt->Fill(Ntracks_nom->at(j),pt_nom->at(j),wt_nom);
      Ntows_nom_v_pt->Fill(Ntows_nom->at(j),pt_nom->at(j),wt_nom);
      M_nom_v_pt->Fill(M_nom->at(j),pt_nom->at(j),wt_nom);
      NEF_nom_v_pt->Fill(NEF_nom->at(j),pt_nom->at(j),wt_nom);
    }
  }
  for (int i = 0; i < TU->GetEntries(); ++ i) {
    TU->GetEntry(i);
    for (int j = 0; j < pt_TU->size(); ++ j) { //all vectors of doubles in the branches should have the same size
      Ntracks_TU_v_pt->Fill(Ntracks_TU->at(j),pt_TU->at(j),wt_TU);
      Ntows_TU_v_pt->Fill(Ntows_TU->at(j),pt_TU->at(j),wt_TU);
      M_TU_v_pt->Fill(M_TU->at(j),pt_TU->at(j),wt_TU);
      NEF_TU_v_pt->Fill(NEF_TU->at(j),pt_TU->at(j),wt_TU);
    }
  }

  TAxis* nom_axis = Ntracks_nom_v_pt->GetYaxis();
  double ranges_d[nBins + 1] = {(double) nom_axis->FindBin(15), (double) nom_axis->FindBin(20), (double) nom_axis->FindBin(25), (double) nom_axis->FindBin(30), (double) nom_axis->FindBin(40), (double) nom_axis->FindBin(60)};//1,2,3,4,6,10};                                                                     

  vector<TH1D*> Ntracks_noms = Projection2D(Ntracks_nom_v_pt,nBins,ranges_d,"x");
  vector<TH1D*> Ntows_noms = Projection2D(Ntows_nom_v_pt,nBins,ranges_d,"x");
  vector<TH1D*> M_noms = Projection2D(M_nom_v_pt,nBins,ranges_d,"x");
  vector<TH1D*> NEF_noms = Projection2D(NEF_nom_v_pt,nBins,ranges_d,"x");
  vector<TH1D*> Ntracks_TUs = Projection2D(Ntracks_TU_v_pt,nBins,ranges_d,"x");
  vector<TH1D*> Ntows_TUs = Projection2D(Ntows_TU_v_pt,nBins,ranges_d,"x");
  vector<TH1D*> M_TUs = Projection2D(M_TU_v_pt,nBins,ranges_d,"x");
  vector<TH1D*> NEF_TUs = Projection2D(NEF_TU_v_pt,nBins,ranges_d,"x");
  
  for (int i = 0; i < nBins; ++ i) {
    cout << "INTEGRALS FOR BIN " << i << ": nom - " << Ntracks_noms[i]->Integral() << " tu - " << Ntracks_TUs[i]->Integral() << endl;
    Prettify1D(Ntracks_noms[i],kBlue, kOpenCircle, 2, kBlue, "N_{tracks}","1/N dN/dN_{tracks}",0,15,0,0.3);
    Prettify1D(Ntows_noms[i],kBlue, kOpenCircle, 2, kBlue, "N_{towers}","1/N dN/dN_{tows}",0,15,0,0.2);
    Prettify1D(M_noms[i],kBlue, kOpenCircle, 2, kBlue, "M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.5);
    Prettify1D(NEF_noms[i],kBlue, kOpenCircle, 2, kBlue, "NEF","1/N dN/dNEF",0,1,0,3);
    Prettify1D(Ntracks_TUs[i],kRed, kOpenCircle, 2, kRed, "N_{tracks}","1/N dN/dN_{tracks}",0,15,0,0.3);
    Prettify1D(Ntows_TUs[i],kRed, kOpenCircle, 2, kRed, "N_{towers}","1/N dN/dN_{tows}",0,15,0,0.2);
    Prettify1D(M_TUs[i],kRed, kOpenCircle, 2, kRed, "M [GeV/c^{2}]","1/N dN/dM",0,14,0,0.5);
    Prettify1D(NEF_TUs[i],kRed, kOpenCircle, 2, kRed, "NEF","1/N dN/dNEF",0,1,0,3);
  }
  

  TCanvas *ctr = new TCanvas ("ctr","ctr",1400,1000);
  DivideCanvas(ctr,"0",3,2);
  TCanvas *ctow = new TCanvas ("ctow","ctow",1400,1000);
  DivideCanvas(ctow,"0",3,2);
  TCanvas *cM = new TCanvas ("cM","cM",1400,1000);
  DivideCanvas(cM,"0",3,2);
  TCanvas *cNEF = new TCanvas ("cNEF","cNEF",1400,1000);
  DivideCanvas(cNEF,"0",3,2);

  TLatex * p = new TLatex ();
  
  TLegend *tleg = new TLegend(0.7,0.6,0.85,0.8); tleg->SetBorderSize(0);
  tleg->AddEntry(Ntracks_noms[0],"nominal","p");
  tleg->AddEntry(Ntracks_TUs[0],"TU","p");
  
  TH1D* hdummytr = new TH1D("hdummytr",";;1/N dN/dN_{tracks}",1,0,15);
  TH1D* hdummytow = new TH1D("hdummytow",";;1/N dN/dN_{tows}",1,0,30);
  TH1D* hdummyM = new TH1D("hdummyM",";;1/N dN/dM",1,0,10);
  TH1D* hdummyNEF = new TH1D("hdummyNEF",";;1/N dN/dNEF",1,0,1);
  
  hdummytr->GetYaxis()->SetRangeUser(0,0.3);
  hdummytow->GetYaxis()->SetRangeUser(0,0.2);
  hdummyM->GetYaxis()->SetRangeUser(0,0.5);
  hdummyNEF->GetYaxis()->SetRangeUser(0,3);

  ctr->cd(1); hdummytr->Draw(); p = PanelTitle(); 
  for (int i = 0; i < nBins; ++ i) {
    ctr->cd(i+2); Ntracks_noms[i]->Draw("same"); Ntracks_TUs[i]->Draw("same"); slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if (i == 0) { tleg->Draw("same"); }
  }
  
  ctow->cd(1); hdummytow->Draw(); p = PanelTitle(); 
  for (int i = 0; i < nBins; ++ i) {
    ctow->cd(i+2); Ntows_noms[i]->Draw("same"); Ntows_TUs[i]->Draw("same"); slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if (i == 0) { tleg->Draw("same"); }
  }

  cM->cd(1); hdummyM->Draw(); p = PanelTitle(); 
  for (int i = 0; i < nBins; ++ i) {
    cM->cd(i+2); M_noms[i]->Draw("same"); M_TUs[i]->Draw("same"); slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if (i == 0) { tleg->Draw("same"); }
  }

  cNEF->cd(1); hdummyNEF->Draw(); p = PanelTitle(); 
  for (int i = 0; i < nBins; ++ i) {
    cNEF->cd(i+2); NEF_noms[i]->Draw("same"); NEF_TUs[i]->Draw("same"); slice->DrawLatexNDC(0.3,0.8,(pts[i]+" < p_{T} < "+pts[i+1]+" GeV/c").c_str()); if (i == 0) { tleg->Draw("same"); }
  }


  
  
  ctr->SaveAs("~/jetmass/plots/systematics/tracking_uncertainty_QA_tracks.pdf");
  ctow->SaveAs("~/jetmass/plots/systematics/tracking_uncertainty_QA_tows.pdf");
  cM->SaveAs("~/jetmass/plots/systematics/tracking_uncertainty_QA_mass.pdf");
  cNEF->SaveAs("~/jetmass/plots/systematics/tracking_uncertainty_QA_NEF.pdf");
  */
  
  return;
}
