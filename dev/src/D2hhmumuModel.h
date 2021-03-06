#ifndef D2HHMUMUMODEL_H
#define D2HHMUMUMODEL_H


#include <map>
#include <string>
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"

class D2hhmumuModel
{
 public:
  // Constructors and Destructor
  D2hhmumuModel();
  ~D2hhmumuModel();

  //PDF in mD0 
  //signal PDF 
 
  RooAbsPdf * Signal(RooRealVar m, RooRealVar dm,
                     RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mDeltaJSU, RooRealVar mGammaJSU,
                     RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar dmDeltaJSU, RooRealVar dmGammaJSU,
		     RooRealVar ResolutionScale, RooRealVar globalShift
                    );


  /*
  RooAbsPdf * CombinatoricBackground(RooRealVar m, RooRealVar dm,
				     RooRealVar mChebyA,RooRealVar mChebyB ,RooRealVar mChebyC ,
				     RooRealVar dmThreshold, RooRealVar dmAlpha
				     );
  */

  RooAbsPdf * CombinatoricBackground(RooRealVar m, RooRealVar dm,                                                                                                                                                 RooRealVar mChebyA,RooRealVar mChebyB ,                                                                                           
                                     RooRealVar dmThreshold, RooRealVar dmAlpha                                                                                                            
                                     );   


  RooAbsPdf* RandomPionBackground(RooRealVar m, RooRealVar dm,
				  RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mDeltaJSU, RooRealVar mGammaJSU,
				  RooRealVar dmThreshold, RooRealVar dmAlpha
				  );

  /*
  RooAbsPdf* D2hhhhBackground(RooRealVar m, RooRealVar dm,
			      RooRealVar mMeanGauss, RooRealVar mSigmaGauss,
			      RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar dmDeltaJSU, RooRealVar dmGammaJSU
			      );
			      
  */

  RooAbsPdf* D2hhhhBackground(RooRealVar m, RooRealVar dm,
			      RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mDeltaJSU, RooRealVar mGammaJSU,
			      RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar dmDeltaJSU, RooRealVar dmGammaJSU
			      );

			      
  RooAbsPdf* D2hhhhRandomPionBackground(RooRealVar m, RooRealVar dm,
                              RooRealVar mMeanGauss, RooRealVar mSigmaGauss,
  	 		      RooRealVar dmThreshold, RooRealVar dmAlpha
                              );
  
   RooAbsPdf* Model(std::string components);
  RooWorkspace GetWorkspace() {return m_ws;}

 private:
  void initializeDefaultYields();
  RooWorkspace m_ws;
};

#endif
