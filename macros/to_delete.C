



void to_delete () {
  string dir = "~/jetmass/";
  string closein = "out/closure/";
  string in = "macros/hists/";
  string file = "full.root";

  TFile* closureFile = new TFile( (dir + closein + file).c_str(), "READ");

  double nAccepted_py_even = 0; double nAccepted_ge_even = 0; double nAccepted_py_odd = 0; double nAccepted_ge_odd = 0;
  double nJets_py_even = 0; double nJets_ge_even = 0; double nJets_py_odd = 0; double nJets_ge_odd = 0;
  double nEntries_even = 0; double nEntries_odd = 0; double nFakes_even = 0; double nFakes_odd = 0;
  double nMisses_even = 0; double nMisses_odd = 0; double nMatches_even = 0; double nMatches_odd = 0;

  double t_nAccepted_py_even = 0; double t_nAccepted_ge_even = 0; double t_nAccepted_py_odd = 0; double t_nAccepted_ge_odd = 0;
  double t_nJets_py_even = 0; double t_nJets_ge_even = 0; double t_nJets_py_odd = 0; double t_nJets_ge_odd = 0;
  double t_nEntries_even = 0; double t_nEntries_odd = 0; double t_nFakes_even = 0; double t_nFakes_odd = 0;
  double t_nMisses_even = 0; double t_nMisses_odd = 0; double t_nMatches_even = 0; double t_nMatches_odd = 0;
  
  
  TTree* t = (TTree*) closureFile->Get("checks");

  t->SetBranchAddress("nAccepted_py_even",&nAccepted_py_even); t->SetBranchAddress("nAccepted_ge_even",&nAccepted_ge_even);
  t->SetBranchAddress("nAccepted_py_odd",&nAccepted_py_odd); t->SetBranchAddress("nAccepted_ge_odd",&nAccepted_ge_odd);
  t->SetBranchAddress("nJets_py_even",&nJets_py_even); t->SetBranchAddress("nJets_ge_even", &nJets_ge_even);
  t->SetBranchAddress("nJets_py_odd", &nJets_py_odd); t->SetBranchAddress("nJets_ge_odd", &nJets_ge_odd);
  t->SetBranchAddress("nEntries_even",&nEntries_even); t->SetBranchAddress("nEntries_odd", &nEntries_odd);
  t->SetBranchAddress("nFakes_even", &nFakes_even); t->SetBranchAddress("nFakes_odd",&nFakes_odd);
  t->SetBranchAddress("nMisses_even",&nMisses_even); t->SetBranchAddress("nMisses_odd", &nMisses_odd);
  t->SetBranchAddress("nMatches_even",&nMatches_even); t->SetBranchAddress("nMatches_odd",&nMatches_odd);

  for (int i = 0; i < t->GetEntries(); ++ i) {
    t->GetEntry(i);
    //Filling 1Ds
    t_nAccepted_py_even += nAccepted_py_even; t_nAccepted_ge_even += nAccepted_ge_even;
    t_nAccepted_py_odd += nAccepted_py_odd; t_nAccepted_ge_odd += nAccepted_ge_odd;
    t_nJets_py_even += nJets_py_even; t_nJets_ge_even += nJets_ge_even;
    t_nJets_py_odd += nJets_py_odd; t_nJets_ge_odd += nJets_ge_odd;
    t_nEntries_even += nEntries_even; t_nEntries_odd += nEntries_odd;
    t_nFakes_even += nFakes_even; t_nFakes_odd += nFakes_odd;
    t_nMisses_even += nMisses_even; t_nMisses_odd += nMisses_odd;
    t_nMatches_even += nMatches_even; t_nMatches_odd += nMatches_odd;
  }
  
  std::cout << "CLOSURE INFORMATION: " << std::endl;
  std::cout << "Number of accepted events (py even, ge even, py odd, ge odd, py total, ge total): " << std::endl;
  std::cout << t_nAccepted_py_even << " " << t_nAccepted_ge_even << " " << t_nAccepted_py_odd << " " << t_nAccepted_ge_odd << " " << (t_nAccepted_py_even + t_nAccepted_py_odd) << " " << (t_nAccepted_ge_even + t_nAccepted_ge_odd) << std::endl;
  std::cout << "Number of accepted jets (py even, ge even, py odd, ge odd, py total, ge total): " << std::endl;
  std::cout << t_nJets_py_even << " " << t_nJets_ge_even << " " << t_nJets_py_odd << " " << t_nJets_ge_odd << " " << (t_nJets_py_even + t_nJets_py_odd) << " " << (t_nJets_ge_even + t_nJets_ge_odd) << std::endl;
  std::cout << "Number of entries in the response (even, odd, total): " << std::endl;
  std::cout << t_nEntries_even << " " << t_nEntries_odd << " " << (t_nEntries_even + t_nEntries_odd) << std::endl;
  std::cout << "Number of fakes in the response (even, odd, total): " << std::endl;
  std::cout << t_nFakes_even << " " << t_nFakes_odd << " " << (t_nFakes_even + t_nFakes_odd) << std::endl;
  std::cout << "Number of misses in the response (even, odd, total): " << std::endl;
  std::cout << t_nMisses_even << " " << t_nMisses_odd << " " << (t_nMisses_even + t_nMisses_odd) << std::endl;
  std::cout << "Number of matches in the response (even, odd, total): " << std::endl;
  std::cout << t_nMatches_even << " " << t_nMatches_odd << " " << (t_nMatches_even + t_nMatches_odd) << std::endl;
  std::cout << "Ratio of fakes, misses, and matches to the number of entries in the response (even, odd, total)" << std::endl;
  std::cout << (t_nFakes_even / (double) t_nEntries_even) << " " << (t_nMisses_even / (double) t_nEntries_even) << " " << (t_nMatches_even / (double) t_nEntries_even) << std::endl;
  std::cout << (t_nFakes_odd / (double) t_nEntries_odd) << " " << (t_nMisses_odd / (double) t_nEntries_odd) << " " << (t_nMatches_odd / (double) t_nEntries_odd) << std::endl;
  std::cout << ((t_nFakes_even+t_nFakes_odd) / (double) (t_nEntries_even+t_nEntries_odd)) << " " << ((t_nMisses_even+t_nMisses_odd) / (double) (t_nEntries_even+t_nEntries_odd)) << " " << ((t_nMatches_even+t_nMatches_odd) / (double) (t_nEntries_even+t_nEntries_odd)) << std::endl;
  
}
