#include "funcs.hh"
#include "params.hh"

namespace Analysis {

  // -------------------------
  // IO/OS Manip functionality
  // -------------------------

  
  // Used to understand which format of input file is being used
  // ( .root file, .txt, .list, etc )
  // ---------------------------------------------------------------------
  bool HasEnding (std::string const &full_string, std::string const &ending) {
    if (full_string.length() >= ending.length()) {
      return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending) );
    } else {
      return false;
    }
  }

  void FillTrees ( std::vector<fastjet::PseudoJet> jets, TTree* Tree, double &jPt, double &jEta, double &jPhi, double &jE, double &jM) {
    for ( int j = 0; j< jets.size(); ++ j) {   // FILL JET INFO
      if (jets[j].pt() < 0.2) continue;
      jPt = jets[j].pt();    jEta = jets[j].eta();    jPhi = jets[j].phi();
      jE = jets[j].e();    jM = jets[j].m();
      Tree->Fill();
    }
  }



}
