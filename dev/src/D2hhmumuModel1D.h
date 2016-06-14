#ifndef D2HHMUMUMODEL1D_H
#define D2HHMUMUMODEL1D_H


#include <map>
#include <string>
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooAddPdf.h"

class D2hhmumuModel1D
{
 public:
  // Constructors and Destructor
  D2hhmumuModel1D();
  ~D2hhmumuModel1D();

  //PDF in mD0 
  //signal PDF 
 
  RooAbsPdf * Signal(RooRealVar m, 
                     RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mDeltaJSU, RooRealVar mGammaJSU,
		     RooRealVar ResolutionScale, RooRealVar globalShift
                    );

    //RooAbsPdf * Signal(RooRealVar m,
    //		     RooRealVar mean, RooRealVar sigma, RooRealVar alphaR, RooRealVar alphaL, RooRealVar nR, RooRealVar nL
    //		   );

  RooAbsPdf * CombinatoricBackground(RooRealVar m,                                                                                                                                                           RooRealVar mChebyA,RooRealVar mChebyB );   

  //RooAbsPdf* D2hhhhBackground(RooRealVar m,
  //					       RooRealVar mean, RooRealVar sigma, RooRealVar alphaR, RooRealVar alphaL, RooRealVar nR, RooRealVar nL
  //			      );

  //RooAbsPdf *D2hhhhBackground(RooRealVar m,
  //		     RooRealVar mean, RooRealVar sigma,  RooRealVar alphaL,  RooRealVar nL
  //		   );

   RooAbsPdf* D2hhhhBackground(RooRealVar m,
              RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mDeltaJSU, RooRealVar mGammaJSU
              );

			        
  RooAbsPdf* Model(std::string components);
  RooWorkspace GetWorkspace() {return m_ws;}

 private:
  void initializeDefaultYields();
  RooWorkspace m_ws;
};

#endif
