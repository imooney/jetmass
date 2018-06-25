

#ifndef PARAMS_HH
#define PARAMS_HH

#define __ERR(message) {std::cerr << "[" << __FILE__ << "::" << __func__ << "()] -- ERR: " << message << std::endl;}
#define __OUT(message) {std::cout << "[" << __FILE__ << "::" << __func__ << "()] -- OUT: " << message << std::endl;}

namespace Analysis {

  //consts
  const double Pi = 3.141592653;
  const double PionMass = 0.13957018;

  const double R = 0.4;                   //jet resolution parameter
  
  //quality cuts
  const double absMaxVz = 30.0;           //|Vz|<=30 cm                                                                                                             
  const double DCA = 3.0;
  const double NFitPts = 20;                                                                                                                            
  const double FitOverMaxPts = 0.52;                                                                                                           
  const double MaxPt = 1000;

  
  //particle cuts
  const double max_track_rap = 1.0;
  const double partMinPt = 0.2;           //particle pT >= 0.2 GeV                                                                                              
 
  //jet cuts
  const double jet_ptmin = 2.0;           //jet pT >= 2.0 GeV
  const double jet_ptmax = 1000.0;        //DEBUG
  const double max_rap = max_track_rap-R; //|eta_jet| < 1-R                                                                                                                                                                                 
  //ghosts
  const int ghost_repeat = 1;
  const double ghost_area = 0.01;
  const double ghost_maxrap = max_rap + 2.0 * R;
   
}

#endif
