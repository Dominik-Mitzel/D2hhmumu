#ifndef DST2D0PIMODEL_H
#define DST2D0PIMODEL_H

/*
  Author: Maurizio Martinelli
  Institute: Nikhef
  Email: maurizio.martinelli@nikhef.nl

  Name: Dst2D0Pimodel
  Description: a class that creates the model to be fit to the 2-dim
     data distribution of the events.

  Date: 13/6/2013
 */

// STL
#include <map>
#include <string>
// ROOT
// RooFit
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

class Dst2D0piModel
{
 public:
  // Constructors and Destructor
  Dst2D0piModel();
  ~Dst2D0piModel();
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
  //			     RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar deltaJSU, RooRealVar gammaJSU);
  RooAbsPdf * MisreconstructedCFD0(RooRealVar m, RooRealVar dm, 
				   RooRealVar mMean, RooRealVar mWidth, RooRealVar dmMean, RooRealVar dmWidth, RooRealVar corr,
				   RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, 
				   RooRealVar deltaJSU, RooRealVar gammaJSU, RooRealVar dmAlpha, RooRealVar dmThreshold, RooRealVar FractionJSU );
  //				  RooRealVar mMean, RooRealVar mWidth, RooRealVar dmThreshold, RooRealVar dmAlpha);
  RooAbsPdf * DmPeaking(RooRealVar m, RooRealVar dm, 
			RooRealVar mSlope, RooRealVar dmMean, RooRealVar dmWidth, RooRealVar dmThreshold, RooRealVar dmAlpha );//,
    //			RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar deltaJSU, RooRealVar gammaJSU, RooRealVar fJSU);
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
