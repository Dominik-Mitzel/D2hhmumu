//class for the implemetation of 2D dm-mD0 fit for D2hhmumu decays
#include <iostream>
#include "sWeights.h"
#include <cmath>
#include <iostream>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include <sstream>
#include <RooDataSet.h>
#include "RooGaussModel.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddModel.h"
#include "RooPolynomial.h"
#include "RooTruthModel.h"
#include "RooFitResult.h"
#include "RooDecay.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooDstD0BG.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooCategory.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooHist.h"
#include "RooStats/SPlot.h"
#include "RooTreeDataStore.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooGlobalFunc.h"
#include "Tools.h"
#include "RooProdPdf.h"
#include "RooJohnsonSU.h"
#include "RooThreshold.h"
#include "D2hhmumuModel.h"

class D2hhmumuFitter {

 public:
 D2hhmumuFitter() ;
  
  ~D2hhmumuFitter();

  //variables used in the fitter

  RooRealVar D0_M;
  RooRealVar deltaM;
  
  //signal
  RooRealVar deltaM_xi;
  RooRealVar deltaM_lambda;
  RooRealVar deltaM_gamma;
  RooRealVar deltaM_delta;

  RooRealVar D0_M_xi;
  RooRealVar D0_M_lambda;
  RooRealVar D0_M_gamma;
  RooRealVar D0_M_delta;

  //purely combinatorial background                                                                                                                                             
  RooRealVar deltaM_threshold;
  RooRealVar deltaM_alpha;

  RooRealVar D0_M_chebyA;
  RooRealVar D0_M_chebyB;
  RooRealVar D0_M_chebyC;
 
  RooRealVar mean1;
  RooRealVar sigma1;
   
  /*
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1800., 1940.,"MeV");
  RooRealVar deltaM("deltaM","#delta m", 139.8,150,"MeV");
  
  //signal
  RooRealVar deltaM_xi("deltaM_xi","deltaM_xi",1.45437e+02,144,146);
  RooRealVar deltaM_lambda("deltaM_lambda","deltaM_lambda",5.35039e-01,0.1,2);
  RooRealVar deltaM_gamma("deltaM_gamma","deltaM_gamma",1.00164e-01,-2,2);
  RooRealVar deltaM_delta("deltaM_delta","deltaM_delta",1.07829e+00,0.,10);

  RooRealVar D0_M_xi("D0_M_xi","D0_M_xi",1.86560e+03,1860,1870);
  RooRealVar D0_M_lambda("D0_M_lambda","D0_M_lambda",8.54839e+00,0.1,20);
  RooRealVar D0_M_gamma("D0_M_gamma","D0_M_gamma",-5.42758e-02,-2,2);
  RooRealVar D0_M_delta("D0_M_delta","D0_M_delta",6.09288e-01,0.,10);

  //purely combinatorial background                                                                                                                                             
  RooRealVar deltaM_threshold("deltaM_threshold","deltaM_threshold",139.57018);
  RooRealVar deltaM_alpha("deltaM_alpha","deltaM_alpha",3.9761e+00,0,10.);

  RooRealVar D0_M_chebyA("D0_M_chebyA","D0_M_chebyA",-3.5906e-02,-1,1);
  RooRealVar D0_M_chebyB("D0_M_chebyB","D0_M_chebyB",-1.7004e-02,-1,1);
  RooRealVar D0_M_chebyC("D0_M_chebyC","D0_M_chebyC",-1.7882e-02,-1,1);
 
  //peaking , in mD0 just a Gaussian... to be improved!                                                                                                                         
  RooRealVar mean1("mu1", "mean1", 1.8350e3,1835.,1845.);
  RooRealVar sigma1("sigma_{1}", "sigma1", 1.4187e+01,3.,25.);
  */
  
  void setKpimumuStartParameters(); //to be implemented, especially when there will be other channles..
  
  void setStyle();
  void fit_MC();
  void fit_Data();

};
 

  
  

