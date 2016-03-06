#include "Dst2D0piModel.h"

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
#include <vector>
// ROOT
// RooFit
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
// USR
#include "RooDstD0BG.h"
#include "RooDmPeakingBkg.h"
#include "RooGauss2D.h"
#include "RooJohnsonSU.h"
#include "RooThreshold.h"
#include "Utilities.h"

// Constructors and Destructor
Dst2D0piModel::Dst2D0piModel()
{
}

Dst2D0piModel::~Dst2D0piModel()
{
}

// Methods
RooAbsPdf *
Dst2D0piModel::CombinatoricBackground(RooRealVar m, RooRealVar dm, 
				      RooRealVar mSlope, RooRealVar dmThreshold, RooRealVar dmAlpha)
{
  if (m_ws.pdf("CombinatoricBkg") == 0) {
    //RooChebychev CombinatoricBkgM("CombinatoricBkgM", "Combinatoric Background (M)", m, RooArgList(mSlope));
    RooExponential CombinatoricBkgM("CombinatoricBkgM", "Combinatoric Background (M)", m, mSlope);
    RooThreshold CombinatoricBkgDm("CombinatoricBkgDm", "Combinatoric Background (#Deltam)", dm, dmThreshold, dmAlpha);
    RooProdPdf CombinatoricBkg("CombinatoricBkg", "Combinatoric Background", RooArgList(CombinatoricBkgM, CombinatoricBkgDm));
    m_ws.import(CombinatoricBkg);
    RooRealVar nCombinatoricBkg("nCombinatoricBkg", "Combinatoric Background Events", 1000, 0., 10000000.);
    m_ws.import(nCombinatoricBkg);
  }
  return m_ws.pdf("CombinatoricBkg");
}

RooAbsPdf *
Dst2D0piModel::Signal(RooRealVar m, RooRealVar dm, 
		      RooRealVar mMean, RooRealVar mWidth, RooRealVar dmMean, RooRealVar dmWidth, RooRealVar corr,
		      RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mNuJSU, RooRealVar mTauJSU, 
		      RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar dmNuJSU, RooRealVar dmTauJSU, 
		      RooRealVar FractionJSU )
//		      RooRealVar dmAlpha, RooRealVar dmThreshold, RooRealVar FractionJSU)
{
  if (m_ws.pdf("Signal") == 0) {
    RooGauss2D SignalCore("SignalCore", "Signal Core", m, dm, mMean, mWidth, dmMean, dmWidth, corr);
    //    RooGaussian SignalD0_JSU("SignalD0_JSU", "Signal D^{0} JSU", m, mMeanJSU, mWidthJSU);
    RooJohnsonSU SignalD0_JSU("SignalD0_JSU", "Signal D^{0} JSU", m, mMeanJSU, mWidthJSU,mNuJSU,mTauJSU);
    RooJohnsonSU SignalDm_JSU("SignalDm_JSU", "Signal #Deltam JSU", dm, dmMeanJSU, dmWidthJSU, dmNuJSU, dmTauJSU);//, dmThreshold, dmAlpha);
    RooProdPdf SignalJSU("SignalJSU","Signal Johnson SU Component",RooArgList(SignalD0_JSU,SignalDm_JSU));
    RooAddPdf Signal("Signal","Signal",RooArgList(SignalJSU,SignalCore),FractionJSU);
    m_ws.import(Signal);
    RooRealVar nSignal("nSignal", "Signal Events", 1000, 0., 1000000.);
    m_ws.import(nSignal);
  }
  return m_ws.pdf("Signal");
}

RooAbsPdf *
Dst2D0piModel::D0Peaking(RooRealVar m, RooRealVar dm, 
			 RooRealVar mMean, RooRealVar mWidth, 
			 RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mNuJSU, RooRealVar mTauJSU, 
			 RooRealVar dmThreshold, RooRealVar dmAlpha, RooRealVar FractionJSU)
{
  if (m_ws.pdf("D0Peaking") == 0) {
    RooGaussian UntaggedD0_Gau("UntaggedD0_Gau", "Untagged D^{0} gaussian", m, mMean, mWidth);    
    if (m_ws.pdf("SignalD0_JSU") == 0) {
      RooJohnsonSU SignalD0_JSU("SignalD0_JSU", "Signal D^{0} JSU", m, mMeanJSU, mWidthJSU,mNuJSU,mTauJSU);
      m_ws.import(SignalD0_JSU);
    }
    RooAddPdf UntaggedD0("UntaggedD0", "Untagged D^{0}", RooArgList(UntaggedD0_Gau, *m_ws.pdf("SignalD0_JSU")), FractionJSU);
    if (m_ws.pdf("CombinatoricBkgDm") == 0) {
      //      RooDstD0BG CombinatoricBkgDm("CombinatoricBkgDm", "Combinatoric Background (#Deltam)", dm, dmThreshold, dmAlpha, dmBeta, dmGamma);
      RooThreshold CombinatoricBkgDm("CombinatoricBkgDm", "Combinatoric Background (#Deltam)", dm, dmThreshold, dmAlpha);
      m_ws.import(CombinatoricBkgDm);
    }
    RooProdPdf D0Peaking("D0Peaking","D^{0} Peaking",RooArgList(*m_ws.pdf("SignalD0_JSU"),*m_ws.pdf("CombinatoricBkgDm")));
    m_ws.import(D0Peaking);
    RooRealVar nD0Peaking("nD0Peaking", "D^{0} Peaking Events", 100, 0., 1000000.);
    m_ws.import(nD0Peaking);
  }
  return m_ws.pdf("D0Peaking");
}

RooAbsPdf *
Dst2D0piModel::Ds2KKPiPiPiBkg(RooRealVar m, RooRealVar dm, 
			      RooRealVar mMean, RooRealVar mWidth, RooRealVar meanCorr, RooRealVar widthCorr,
			      RooRealVar dmThreshold, RooRealVar dmAlpha)
{
  // This distribution has a combinatorial shape in dm and a Gaussian with mean and width correlated to dm
  if (m_ws.pdf("Ds2KKPiPiPiBkg") == 0) {
    // definition
    if (m_ws.pdf("CombinatoricBkgDm") == 0) {
      RooThreshold CombinatoricBkgDm("CombinatoricBkgDm", "Combinatoric Background (#Deltam)", dm, dmThreshold, dmAlpha);
      m_ws.import(CombinatoricBkgDm);
    }
    RooFormulaVar DsMean("DsMean","DsMean","@0+@1*(@2-@3)",RooArgList(mMean,meanCorr,dm,dmThreshold));
    RooFormulaVar DsWidth("DsWidth","Dswidth","@0+@1*(@2-@3)",RooArgList(mWidth,widthCorr,dm,dmThreshold));
    RooGaussian Ds2KKPiPiPiBkgM("Ds2KKPiPiPiBkgM","D_{s}^{+}#rigtharrowK^{+}K^{-}#pi^{+}pi^{+}pi^{-} (M)", m, DsMean, DsWidth);
    RooProdPdf Ds2KKPiPiPiBkg("Ds2KKPiPiPiBkg","D_{s}^{+}#rigtharrowK^{+}K^{-}#pi^{+}pi^{+}pi^{-} bkg", RooArgList(Ds2KKPiPiPiBkgM, *m_ws.pdf("CombinatoricBkgDm")));
    m_ws.import(Ds2KKPiPiPiBkg);
    RooRealVar nDs2KKPiPiPiBkg("nDs2KKPiPiPiBkg", "D_{s}^{+}#rigtharrowK^{+}K^{-}#pi^{+}pi^{+}pi^{-} Events", 10, 0., 1000000.);
    m_ws.import(nDs2KKPiPiPiBkg);
  }
  return m_ws.pdf("Ds2KKPiPiPiBkg");
}

RooAbsPdf *
Dst2D0piModel::D02KKPiPiPiBkg(RooRealVar m, RooRealVar dm, 
			      RooRealVar mMean, RooRealVar mWidth,
			      //			      RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar deltaJSU, RooRealVar gammaJSU
			      RooRealVar dmMean, RooRealVar dmWidth)
{
  // This distribution has a signal-like peaking shape in dm and a Gaussian for m
  if (m_ws.pdf("D02KKPiPiPiBkg") == 0) {
    // definition
    /*
    if (m_ws.pdf("SignalDm_JSU") == 0) {
      RooJohnsonSU SignalDm_JSU("SignalDm_JSU", "Signal #Deltam JSU", dm, dmMeanJSU, dmWidthJSU, deltaJSU, gammaJSU);
      m_ws.import(SignalDm_JSU);
    }
    */
    RooGaussian D02KKPiPiPiBkgM("D02KKPiPiPiBkgM","D^{0}#rigtharrowK^{+}K^{-}#pi^{+}pi^{+}pi^{0} (M)", m, mMean, mWidth);
    RooGaussian D02KKPiPiPiBkgDm("D02KKPiPiPiBkgDm","D^{0}#rigtharrowK^{+}K^{-}#pi^{+}pi^{+}pi^{0} (Dm)", dm, dmMean, dmWidth);
    RooProdPdf D02KKPiPiPiBkg("D02KKPiPiPiBkg","D^{0}#rigtharrowK^{+}K^{-}#pi^{+}pi^{+}pi^{0} bkg", RooArgList(D02KKPiPiPiBkgM, D02KKPiPiPiBkgDm) );//*m_ws.pdf("SignalDm_JSU")));
    m_ws.import(D02KKPiPiPiBkg);
    RooRealVar nD02KKPiPiPiBkg("nD02KKPiPiPiBkg", "D^{0}#rigtharrowK^{+}K^{-}#pi^{+}pi^{+}pi^{0} Events", 10, 0., 1000000.);
    m_ws.import(nD02KKPiPiPiBkg);
  }
  return m_ws.pdf("D02KKPiPiPiBkg");
}

RooAbsPdf *
Dst2D0piModel::MisreconstructedCFD0(RooRealVar m, RooRealVar dm, 
				    RooRealVar mMean, RooRealVar mWidth, RooRealVar dmMean, RooRealVar dmWidth, RooRealVar corr,
				    RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, 
				    RooRealVar deltaJSU, RooRealVar gammaJSU, RooRealVar dmAlpha, RooRealVar dmThreshold, RooRealVar FractionJSU)
//				   RooRealVar mMean, RooRealVar mWidth, RooRealVar dmThreshold, RooRealVar dmAlpha)
{
  if (m_ws.pdf("MisreconstructedCFD0") == 0) {
    /* old
    RooGaussian MisreconstructedD0("MisreconstructedD0", "Misreconstructed D^{0}", m, mMean, mWidth);
    if (m_ws.pdf("CombinatoricBkgDm") == 0) {
      RooThreshold CombinatoricBkgDm("CombinatoricBkgDm", "Combinatoric Background (#Deltam)", dm, dmThreshold, dmAlpha);
      m_ws.import(CombinatoricBkgDm);
    }
    RooProdPdf MisreconstrucedCFD0("MisreconstrucedCFD0","Misreconstructed CF D^{0}",RooArgList(MisreconstructedD0,*m_ws.pdf("CombinatoricBkgDm")));
    */
    // new with signal-like model
    RooGauss2D MisreconstructedCFD0("MisreconstructedCFD0", "Misreconstructed D^{0} Core", m, dm, mMean, mWidth, dmMean, dmWidth, corr);
    /*
    RooGaussian MisreconstructedCFD0_D0_JSU("MisreconstructedCFD0_D0_JSU", "Misreconstructed D^{0} JSU", m, mMeanJSU, mWidthJSU);
    RooJohnsonSU MisreconstructedCFD0_Dm_JSU("MisreconstructedCFD0_Dm_JSU", "Misreconstructed D^{0} #Deltam JSU", dm, dmMeanJSU, dmWidthJSU, deltaJSU, gammaJSU);//, dmThreshold, dmAlpha);
    RooProdPdf MisreconstructedCFD0_JSU("MisreconstructedCFD0_JSU","Misreconstructed D^{0} Johnson SU Component",RooArgList(MisreconstructedCFD0_D0_JSU,MisreconstructedCFD0_Dm_JSU));
    RooAddPdf MisreconstructedCFD0("MisreconstructedCFD0","Misreconstructed D^{0}",RooArgList(MisreconstructedCFD0_JSU,MisreconstructedCFD0Core),FractionJSU);
    */
    m_ws.import(MisreconstructedCFD0);
    RooRealVar nMisreconstructedCFD0("nMisreconstructedCFD0", "Misreconstructed CF D^{0} Events", 10, 0., 100000.);
    m_ws.import(nMisreconstructedCFD0);
  }
  return m_ws.pdf("MisreconstructedCFD0");
}

RooAbsPdf *
Dst2D0piModel::DmPeaking(RooRealVar m, RooRealVar dm, 
			 RooRealVar mSlope, RooRealVar dmThreshold, RooRealVar dmAlpha, RooRealVar dmBeta, RooRealVar dmGamma )
			 //			 RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar deltaJSU, RooRealVar gammaJSU, RooRealVar fJSU)
{
  if (m_ws.pdf("DmPeaking") == 0) {
    if (m_ws.pdf("CombinatoricBkgM") == 0) {
      RooChebychev CombinatoricBkgM("CombinatoricBkgM", "Combinatoric Background (M)", m, RooArgList(mSlope));
      m_ws.import(CombinatoricBkgM);
    }
    /*
    RooGaussian DmMissingD0Gaussian("DmMissingD0Gaussian", "#Deltam Missing D^{0} Gaussian", dm, dmMean, dmWidth);
    RooThreshold DmMissingD0Threshold("DmMissingD0Threshold", "#Deltam Missing D^{0} Threshold", dm, dmThreshold, dmAlpha);
    RooProdPdf DmMissingD0("DmMissingD0", "#Deltam Missing D^{0}", RooArgList(DmMissingD0Gaussian,DmMissingD0Threshold) );
    */
    //    RooDmPeakingBkg DmMissingD0("DmMissingD0", "#Deltam Missing D^{0}",dm,dmThreshold,dmAlpha,dmMean,dmWidth);
    RooDstD0BG DmMissingD0("DmMissingD0", "#Deltam Missing D^{0}",dm,dmThreshold,dmAlpha,dmBeta,dmGamma);
    // Add a peaking component like signal - CF senza una particella, che viene poi recuperata dal fondo.
    //    RooJohnsonSU DmPeaking_JSU("DmPeaking_JSU", "#Deltam Peaking JSU", dm, dmMeanJSU, dmWidthJSU, deltaJSU, gammaJSU);
    //    RooAddPdf DmPeaking_dm("DmPeaking_dm", "#Deltam Peaking Sum", RooArgList(DmPeaking_JSU, DmMissingD0), fJSU);
    RooProdPdf DmPeaking("DmPeaking","#Deltam Peaking",RooArgList(*m_ws.pdf("CombinatoricBkgM"),DmMissingD0));
    m_ws.import(DmPeaking);
    RooRealVar nDmPeaking("nDmPeaking", "#Deltam Peaking Events", 100, 0., 1000000.);
    m_ws.import(nDmPeaking);
  }
  return m_ws.pdf("DmPeaking");
}

RooAbsPdf *
Dst2D0piModel::Model(const char* components)
{
  Utilities ut;
  std::vector<std::string> comps = ut.splitString(components, ",");
  if (m_ws.pdf("Dst2D0piModel") == 0) {
    RooArgList pdfList, coefList;
    for (std::vector<std::string>::iterator it = comps.begin();
	 it != comps.end(); ++it) {
      if (m_ws.pdf(it->c_str()) == 0) {
	std::cout << "Error(Dst2D0piModel::Model): no component found with name " << *it << std::endl;
	return 0;
      }
      pdfList.add(*m_ws.pdf(it->c_str()));
      coefList.add(*m_ws.var((std::string("n")+*it).c_str()) );
    }
    pdfList.Print();
    coefList.Print();
    RooAddPdf Model("Dst2D0piModel", "D*^{+} #rightarrow D^{0}#pi^{+} Model", pdfList, coefList);
    m_ws.import(Model);
  }
  return m_ws.pdf("Dst2D0piModel");
}
