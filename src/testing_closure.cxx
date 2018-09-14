//  Isaac Mooney 9/6/2018 - tests closure in general

#include "params.hh"
#include "funcs.hh"
#include "RooUnfold.h"
#include "TRandom.h"

using namespace std;
using namespace Analysis;



int main () {
  //first tested simple case in which all jets are matched. (No misses, no fakes). Exactly produces (within errors) ratio of 1.
  //next, testing case in which 5% of "jets" are missed and 5% are fakes.
  

  std::string        outputDir        = "out/closure/";                        // directory where everything will be saved
  std::string     outFileName        = "test.root";
  
  //  TStarJetPicoDefinitions::SetDebugLevel(0);
  TH1::SetDefaultSumw2( );  // Histograms will calculate gaussian errors
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );

  //temp hists
  TH1D *gen_odd = new TH1D("gen_odd","",100,0,20); TH1D *det_odd = new TH1D("det_odd","",100,0,20);
  TH1D *gen_even = new TH1D("gen_even","",100,0,20); TH1D *det_even = new TH1D("det_even","",100,0,20);
  
  //responses for MC Closure test
  RooUnfoldResponse odd(100,0,20,100,0,20,"odd","");
  RooUnfoldResponse even(100,0,20,100,0,20,"even","");
 
  TRandom *r = new TRandom();
  
  for(int nTimes = 0; nTimes < 1e5; nTimes ++) {
    //EVENS
    if (nTimes % 2 == 0) {
      if (nTimes % 10 != 0) { //not a miss or fake.
	double truth = r->Exp(2);
	double measured = r->Exp(1);
	//filling
	gen_even->Fill(truth);
	det_even->Fill(measured);
	
	even.Fill(measured, truth);
      }
      else {
	if (nTimes % 20 != 0) { //let's say this is a miss
	  double truth = r->Exp(2);
	  gen_even->Fill(truth);
	  
	  even.Miss(truth);
	}
	else { //so this will be a fake
	  double measured = r->Exp(1);
	  det_even->Fill(measured);
	  
	  even.Fake(measured);
	}
      }
    } //END EVENS
  
    //BEGIN ODDS
    else if (nTimes % 2 != 0) {
      if ((nTimes - 1) % 10 != 0) { //not a miss or fake
	double truth = r->Exp(2);
	double measured = r->Exp(1);
	//filling
	gen_odd->Fill(truth);
	det_odd->Fill(measured);
	
	odd.Fill(measured, truth);
      }
      else {
	if ((nTimes - 1) % 20 != 0) { //let's say this is a miss
	  double truth = r->Exp(2);
	  gen_odd->Fill(truth);
	  
	  odd.Miss(truth);
	}
	else { //so this will be a fake
	  double measured = r->Exp(1);
	  det_odd->Fill(measured);
	  
	  odd.Fake(measured);
	}
      }
    } //END ODDS
  }
  
  //~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END EVENT LOOP! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  
  TFile *fout = new TFile( ( outputDir + outFileName ).c_str() ,"RECREATE");

  odd.Write(); even.Write();

  gen_odd->Write(); det_odd->Write(); gen_even->Write(); det_even->Write();
  fout->Close();
  
  std::cout << "Wrote to " << fout->GetName() << std::endl;

  return 0;
}
