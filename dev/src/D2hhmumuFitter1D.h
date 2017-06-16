#ifndef D2HHMUMUFITTER1D_H
#define D2HHMUMUFITTER1D_H

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
#include "D2hhmumuModel1D.h"
#include "RooWorkspace.h"
#include <fstream>


class D2hhmumuFitter1D {

 public:
  //D2hhmumuFitter1D() ;
  D2hhmumuFitter1D(TString results="") ;

  ~D2hhmumuFitter1D();

  TRandom3* generator;
  double Nsig_offset;
  double BFsig_offset;

  RooRealVar Nsig_blinding; 
  RooRealVar BFsig_blinding;

  //variables used in the fitter
  //
  //************************

  
  //signal

  RooRealVar D0_M_xi;
  RooRealVar D0_M_lambda;
  RooRealVar D0_M_gamma;
  RooRealVar D0_M_delta;

  RooRealVar D0_M_xi_norm;
  RooRealVar D0_M_lambda_norm;
  RooRealVar D0_M_gamma_norm;
  RooRealVar D0_M_delta_norm;

  RooRealVar ResolutionScale;
  RooRealVar globalShift;

  //purely combinatorial background                                                                                                                                             

  RooRealVar D0_M_chebyA;
  RooRealVar D0_M_chebyB;
  RooRealVar D0_M_chebyC;

  RooRealVar D0_M_xi_bkg;
  RooRealVar D0_M_lambda_bkg;
  RooRealVar D0_M_gamma_bkg;
  RooRealVar D0_M_delta_bkg;
  RooRealVar D0_M_ExpoLambda;

  RooRealVar D0_M_xi_bkg_norm;
  RooRealVar D0_M_lambda_bkg_norm;
  RooRealVar D0_M_gamma_bkg_norm;
  RooRealVar D0_M_delta_bkg_norm;

  RooRealVar D0_M_ExpoLambda_norm;

  RooRealVar EffRatio;
  RooRealVar nNorm;
  RooRealVar BFsig;
  RooRealVar BFnorm;

  RooRealVar D0_M_mean;
  RooRealVar D0_M_sigma;
  RooRealVar D0_M_alphaR;
  RooRealVar D0_M_alphaL;
  RooRealVar D0_M_nL;
  RooRealVar D0_M_nR;

  RooRealVar D0_M_mean_bkg;
  RooRealVar D0_M_sigma_bkg;
  RooRealVar D0_M_alphaR_bkg;
  RooRealVar D0_M_alphaL_bkg;
  RooRealVar D0_M_nL_bkg;
  RooRealVar D0_M_nR_bkg;


  // binning scheme

  std::vector<double> rangesKpi_low = {675};
  std::vector<double> rangesKpi_high = {875};
  std::vector<double> rangespipi_low = {200,525,565,950,1100};
  std::vector<double> rangespipi_high = {525,565,950,1100,1600};
  std::vector<double> rangesKK_low = {200,525,565};
  std::vector<double> rangesKK_high = {525,565,900};//!!//THIS IS THE NEW BINNING                                                                                                             
  double binsKpi[2]={675,875};
  double binsKK[4]={200,525,565,900};
  double binspipi[6]={200,525,565,950,1100,1600};

  
  //functions
  //                                                                                                                                                                            
  //************************                                                                                                                                        
            
  std::ofstream myFile;
  

  RooWorkspace  initializeModel(D2hhmumuModel1D* myModel, RooRealVar D0_M);
  RooWorkspace  initializeNormalizationModel(D2hhmumuModel1D* myModel, RooRealVar D0_M);
  void setStyle();

  TString pathToSignalData;
  TString pathToNormData;
  TString pathToSignalMC;
  TString pathToNormMC;
  TString pathToInvData;
  TString pathToSidebandData;
  TString pathToKpipipiData;
  TString pathToKpipipiHistoData;
  TString pathToHHpipiData;
  
  TString q2RangeNormalizationMode;

  void setPathToHHpipiData(TString path);
  void setPathToSignalMC(TString path);
  void setPathToNormMC(TString path);
  void setPathToSignalData(TString path);
  void setPathToNormData(TString path);
  void setPathToInvData(TString path);
  void setPathToKpipipiData(TString path);
  void setPathToSidebandData(TString path);
  void setPathToKpipipiHistoData(TString path);

  void setKpimumuStartParameters(); //to be implemented, especially when there will be other channles..    
  void setKKmumuStartParameters(); //to be implemented, especially when there will be other channles..                                                    
  void setpipimumuStartParameters(); //to be implemented, especially when there will be other channles..                                              
  //used for selection optimisation
  double getMisIDbkgExp(TString cut,TString namePlot);
  double getCombBkg(TString cut,TString namePlot);
  double getCombBkgFromDeltaM(TString cut,TString namePlot);

  double get_Significance(TString kind, TString cut,TString q2Range,TString namePlot,TString xLabel,TString legend,bool bkgShapeFromInvertedBDTCut);

  //actually the fits
  void fit_MC(TString cut="",bool fixShape = true,TString namePlot="",TString xLabel="m(hh#mu#mu)",TString legend="");
  void fit_normalization_MC(TString cut, bool fixShape, TString namePlot);
  void fit_PIDinverted_Data(bool fixShape,TString namePlot);
  void fit_Data(TString kind="D2KKmumu", TString cut="", TString q2Range="",TString namePlot="",TString xLabel="m(hh#mu#mu)",TString legend="",bool bkgShapeFromInvertedBDTCut=false);
  void fit_unblinded_Data(TString kind="D2KKmumu", TString cut="", TString q2Range="",TString namePlot="",TString xLabel="m(hh#mu#mu)",TString legend="",bool bkgShapeFromInvertedBDTCut=false);
  double fit_resonant_Data(TString kind="D2KKmumu", TString cut="",TString q2Range="",TString namePlot="",TString xLabel="m(hh#mu#mu)",TString legend="",bool fixBkgShape=false);
  void fit_Kpipipi_misID(TString cut="",bool fixShape=false,TString namePlot="", bool fixCombBkgShape=true);
  void fit_HHpipi_misID(TString cut="",bool fixShape = true,TString namePlot="",TString xLabel="m_{hhhh}(hh#mu#mu)",TString legend="",bool fixCombBkgShape=true);
  double fit_HHpipi(TString cut,TString namePlot);
  void fit_Kpipipi_misID_fromHistogramm(TString cut,bool fixShape,TString namePlot);
  double fit_normalization_Data(TString cut,TString namePlot);  
  //std::pair<std::pair<double,double>,std::pair<double,double>>getBkgFromBlindedFit(TString cut,TString q2Range, TString namePlot);
  std::pair<double,double> getBkgFromBlindedFit(TString cut,TString q2Range, TString namePlot);
  //void fillWorkspace(RooWorkspace &ws,RooRealVar D0_M,RooRealVar deltaM);
  double fit_invertedBDT_data(TString cut,TString q2Range,TString namePlot,TString xLabel,TString legend);
  void GausExpModel(int nsig ,int nbkg ) ;
  void fillModelConfig(TString kind, TString dataCut,TString nomalizationCut,Int_t q2Bin, TString fileName);
  void fillModelConfigFake2011(TString kind, TString dataCut,TString nomalizationCut,TString misIDCut, TString name);
  void makeToyStudy(TString kind, TString dataCut,TString nomalizationCut,TString misIDCut,TString targetFile,double nSig_exp, double nCombBkg_exp, double nMisID_exp);  
  void makeSimpleToyStudy(TString kind, TString dataCut,TString q2Range, TString misIDCut,TString targetFile,double nSig_exp, double nCombBkg_exp, double nMisID_exp);
  void addNormalizationSWeights(TString dataCut, TString misIDCut,TString fIn,TString fOut);
  void addNormalizationSWeightsHadronicChannel(TString dataCut,TString fIn, TString fOut, TString plot);
  void addSignalSWeights(TString kind, TString cut,TString q2Range,TString namePlot,TString xLabel,TString legend, bool bkgShapeFromInvertedBDTCut, TString fOut);
  void defineBinning();  
};
 
#endif


  
  

