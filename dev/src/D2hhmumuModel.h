#ifndef DST2D0PIMODEL_H
#define DST2D0PIMODEL_H


#include <map>
#include <string>
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

class D2hhmumuModel
{
 public:
  // Constructors and Destructor
  D2hhmumuModel();
  ~D2hhmumuModel();

  //PDF in mD0 
  //signal PDF 
  RooAbsPdf* signal_D0_M( RooRealVar D0_M, RooRealVar D0_M_xi, RooRealVar D0_M_lambda, RooRealVar D0_M_gamma, RooRealVar D0_M_delta);
  //combinatoric Background
  RooAbsPdf* CombinatoricBackground_D0_M( RooRealVar D0_M, RooRealVar D0_M_chebyA, RooRealVar D0_M_chebyB ,RooRealVar D0_M_chebyC);
  //random Pions
  RooAbsPdf* randomPion_D0_M(RooRealVar D0_M, RooRealVar D0_M_xi, RooRealVar D0_M_lambda, RooRealVar D0_M_gamma, RooRealVar D0_M_delta);
  //misidentified D2hhhh background
  RooAbsPdf* D2hhhh_D0_M(RooRealVar D0_M,RooRealVar D0_M_mean,RooRealVar D0_M_sigma);
  //total BKG
  RooAbsPdf* totalBkg_D0_M( RooAbsPdf* );
  RooAbsPdf* model_D0_M();

  //PDF in deltaM  
  //signal                                                                                                                                                  
  RooAbsPdf* signal_deltaM( RooRealVar deltaM, RooRealVar deltaM_xi, RooRealVar deltaM_lambda, RooRealVar deltaM_gamma, RooRealVar deltaM_delta);
  //combinatoric Background                                                                                                                                     
  RooAbsPdf* CombinatoricBackground_deltaM();
  //random Pions                                                                                                                                                
  RooAbsPdf* randomPion_deltaM();
  //misidentified D2hhhh background                                                                                                                             
  RooAbsPdf* D2hhhh_deltaM();
  //total BKG                                                                                                                                                   
  RooAbsPdf* totalBkg_deltaM();

  //2D PDFs
 //signal                                                                                                                                                  
  RooAbsPdf* signal( RooRealVar deltaM, RooRealVar deltaM_xi, RooRealVar deltaM_lambda, RooRealVar deltaM_gamma, RooRealVar deltaM_delta);
  //combinatoric Background                                                                                                                                     
  RooAbsPdf* CombinatoricBackground();
  //random Pions                                                                                                                                                
  RooAbsPdf* randomPion();
  //misidentified D2hhhh background                                                                                                                             
  RooAbsPdf* D2hhhh();
  //total BKG                                                                                                                                                   
  RooAbsPdf* totalBkg();

  


  // Methods
  RooAbsPdf * CombinatoricBackground(RooRealVar m, RooRealVar dm, 
				     RooRealVar mSlope, RooRealVar dmThreshold, RooRealVar dmAlpha);
  RooAbsPdf * D0Peaking(RooRealVar m, RooRealVar dm, 
			RooRealVar mMean, RooRealVar mWidth, 
			RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mDeltaJSU, RooRealVar mGammaJSU, 
			RooRealVar dmThreshold, RooRealVar dmAlpha, RooRealVar FractionJSU);
  RooAbsPdf * Ds2KKPiPiPiBkg(RooRealVar m, RooRealVar dm, 
			     RooRealVar mMean, RooRealVar mWidth, RooRealVar meanCorr, RooRealVar widthCorr,
			     RooRealVar dmThreshold, RooRealVar dmAlpha);
  RooAbsPdf * D02KKPiPiPiBkg(RooRealVar m, RooRealVar dm, 
			     RooRealVar mMean, RooRealVar mWidth,
			     RooRealVar dmMean, RooRealVar dmWidth);
  //     RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar deltaJSU, RooRealVar gammaJSU);
  RooAbsPdf * MisreconstructedCFD0(RooRealVar m, RooRealVar dm, 
				   RooRealVar mMean, RooRealVar mWidth, RooRealVar dmMean, RooRealVar dmWidth, RooRealVar corr,
				   RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, 
				   RooRealVar deltaJSU, RooRealVar gammaJSU, RooRealVar dmAlpha, RooRealVar dmThreshold, RooRealVar FractionJSU );
  //  RooRealVar mMean, RooRealVar mWidth, RooRealVar dmThreshold, RooRealVar dmAlpha);
  RooAbsPdf * DmPeaking(RooRealVar m, RooRealVar dm, 
			RooRealVar mSlope, RooRealVar dmMean, RooRealVar dmWidth, RooRealVar dmThreshold, RooRealVar dmAlpha );//,
  //RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar deltaJSU, RooRealVar gammaJSU, RooRealVar fJSU);
  RooAbsPdf * Signal(RooRealVar m, RooRealVar dm, 
		     RooRealVar mMean, RooRealVar mWidth, RooRealVar dmMean, RooRealVar dmWidth, RooRealVar corr,
		     RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mDeltaJSU, RooRealVar mGammaJSU, 
		     RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar dmDeltaJSU, RooRealVar dmGammaJSU, 
		     RooRealVar dmAlpha, RooRealVar dmThreshold, RooRealVar FractionJSU);
  RooAbsPdf * Model(const char* components);
  RooWorkspace GetWorkspace() {return m_ws;}
 private:
  // Variables
  RooWorkspace m_ws;
};

#endif
