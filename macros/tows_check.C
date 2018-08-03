
void tows_check () {
  TFile* f = new TFile("../out/data/full.root", "READ");
  
  TCanvas * c1d = new TCanvas ("c1d","c1d",800,800);
  TCanvas * c2d = new TCanvas ("c2d","c2d",800,800); c2d->SetLogz();

  TH1D* tow_freq = (TH1D*)f->Get("tow_freq");
  TH2D* tow_id_v_e = (TH2D*)f->Get("tow_id_v_e");

  double run_mean = 0, run_var = 0;
  
  for (int i = 1; i < tow_freq->GetNbinsX(); ++ i) {
    run_mean += tow_freq->GetBinContent(i);
  }
  run_mean /= (double) tow_freq->GetNbinsX();

  cout << setprecision(8) << "mean! " << run_mean << endl;

  for (int i = 1; i < tow_freq->GetNbinsX(); ++ i) {
    run_var += (tow_freq->GetBinContent(i) - run_mean)*(tow_freq->GetBinContent(i) - run_mean);
  }
  run_var /= (double) tow_freq->GetNbinsX();
  run_var = sqrt(run_var);
  
  cout << "sigma! " << run_var << endl;
  cout << "3sigma! " << 3*run_var << endl;
  cout << "mean - 3sigma! " << run_mean - 3*run_var << endl;
  cout << "mean + 3sigma! " << run_mean + 3*run_var << endl;

  const double low = run_mean - 3*run_var, high = run_mean + 3*run_var;
  
  TF1 * lo = new TF1("lo","[0]",0,4801);
  TF1 * hi = new TF1("hi","[0]",0,4801);
  
  lo->SetParameter(0, low);
  hi->SetParameter(0, high);

  cout << "bad tows are..." << endl;
  for (int i = 1; i < tow_freq->GetNbinsX(); ++ i) {
    if (tow_freq->GetBinContent(i) < low || tow_freq->GetBinContent(i) > high) {
      cout << i << ", ";
    }
  }
  cout << endl;
  cout << "dead tows are..." << endl;
  for (int i = 1; i < tow_freq->GetNbinsX(); ++ i) {
    if (tow_freq->GetBinContent(i) == 0) {
      cout << i << ", ";
    }
  }
  cout << endl;
  cout << "mostly dead tows are..." << endl;
  for (int i = 1; i < tow_freq->GetNbinsX(); ++ i) {
    if (tow_freq->GetBinContent(i) < low) {
      cout << i << ", ";
    }
  }
  cout << endl;
  cout << "hot tows are..." << endl;
  for (int i = 1; i < tow_freq->GetNbinsX(); ++ i) {
    if (tow_freq->GetBinContent(i) > high) {
      cout << i << ", ";
    }
  }
  cout << endl;

  c1d->cd(); tow_freq->Draw(); lo->Draw("same"); hi->Draw("same");
  c2d->cd(); tow_id_v_e->Draw("colz");

  return;
}
