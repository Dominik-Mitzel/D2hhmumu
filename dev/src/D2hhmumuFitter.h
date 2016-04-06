#ifndef D2HHMUMUFITTER_H
#define D2HHMUMUFITTER_H

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
  //
  //************************

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


  RooRealVar deltaM_xi_bkg;
  RooRealVar deltaM_lambda_bkg;
  RooRealVar deltaM_gamma_bkg;
  RooRealVar deltaM_delta_bkg;

  RooRealVar D0_M_xi_bkg;
  RooRealVar D0_M_lambda_bkg;
  RooRealVar D0_M_gamma_bkg;
  RooRealVar D0_M_delta_bkg;
 
  RooRealVar mean1;
  RooRealVar sigma1;

  
  //functions
  //                                                                                                                                                                            
  //************************                                                                                                                                                    
  

  D2hhmumuModel*  initializeModel();
  void setStyle();

  TString pathToSignalData;
  TString pathToNormData;
  TString pathToSignalMC;
  TString pathToInvData;
  TString pathToSidebandData;

  void setPathToSignalMC(TString path);
  void setPathToSignalData(TString path);
  void setPathToNormData(TString path);
  void setPathToInvData(TString path);
  void setPathToSidebandData(TString path);


  void setKpimumuStartParameters(); //to be implemented, especially when there will be other channles..    
  void setKKmumuStartParameters(); //to be implemented, especially when there will be other channles..                                                                        
  void setpipimumuStartParameters(); //to be implemented, especially when there will be other channles..                                                                      

  double getMisIDbkgExp(TString cut,TString namePlot);
  double getCombBkg(TString cut,TString namePlot);
  void fit_MC(bool fixShape);
  void fit_PIDinverted_Data(bool fixShape);
  void fit_Data(TString cut);
  void toyStudy();
  
  
};
 
#endif


  
  

