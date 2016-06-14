#include "D2hhmumuModel.h"           
#include <vector>
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooJohnsonSU.h"
#include "RooThreshold.h"
#include "RooExtendPdf.h"
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "RooUnblindPrecision.h"
// Constructors and Destructor                                                                                                                       
D2hhmumuModel::D2hhmumuModel()
{
}

D2hhmumuModel::~D2hhmumuModel()
{
}

RooAbsPdf* D2hhmumuModel::Signal(RooRealVar m, RooRealVar dm,
                      RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mNuJSU, RooRealVar mTauJSU,
				 RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar dmNuJSU, RooRealVar dmTauJSU,
				 RooRealVar ResolutionScale, RooRealVar globalShift)
{
  if (m_ws.pdf("Signal") == 0) {
    RooFormulaVar mMean("mMean","@0*@1",RooArgList(mMeanJSU,globalShift));
    RooFormulaVar mWidth("mWidth","@0*@1",RooArgList(mWidthJSU,ResolutionScale));
    RooFormulaVar dmMean("dmMean","@0*@1",RooArgList(dmMeanJSU,globalShift));
    RooFormulaVar dmWidth("dmWidth","@0*@1",RooArgList(dmWidthJSU,ResolutionScale));
    RooJohnsonSU SignalD0_JSU("SignalD0_JSU", "Signal D^{0} JSU", m, mMean, mWidth,mNuJSU,mTauJSU);
    RooJohnsonSU SignalDm_JSU("SignalDm_JSU", "Signal #Deltam JSU", dm, dmMean, dmWidth, dmNuJSU, dmTauJSU);         
    RooProdPdf Signal("Signal","Signal PDF",RooArgList(SignalD0_JSU,SignalDm_JSU));
    m_ws.import(Signal);
  }
  return m_ws.pdf("Signal");
}



RooAbsPdf* D2hhmumuModel::CombinatoricBackground(RooRealVar m, RooRealVar dm,
				      RooRealVar mChebyA,RooRealVar mChebyB ,
				      RooRealVar dmThreshold, RooRealVar dmAlpha)
{
  if (m_ws.pdf("CombinatoricBkg") == 0) {
    RooChebychev CombinatoricBkgM("CombinatoricBkgM", "Combinatoric Background (M)", m, RooArgList(mChebyA,mChebyB));
    RooThreshold CombinatoricBkgDm("CombinatoricBkgDm", "Combinatoric Background (#Deltam)", dm, dmThreshold, dmAlpha);
    RooProdPdf CombinatoricBkg("CombinatoricBkg", "Combinatoric Background", RooArgList(CombinatoricBkgM, CombinatoricBkgDm));
    m_ws.import(CombinatoricBkg);
  }
  return m_ws.pdf("CombinatoricBkg");
}



RooAbsPdf* D2hhmumuModel::RandomPionBackground(RooRealVar m, RooRealVar dm,
				RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mDeltaJSU, RooRealVar mGammaJSU,
				RooRealVar dmThreshold, RooRealVar dmAlpha
				)
{
  if (m_ws.pdf("RandomPionBkg") == 0) {
    RooJohnsonSU RandomPionM_JSU("RandomPion_JSU", "Random Pion D^{0} JSU", m, mMeanJSU, mWidthJSU,mDeltaJSU,mGammaJSU);
    RooThreshold RandomPionBkgDm("RandomPionBkgDm", "Random Pion Background (#Deltam)", dm, dmThreshold, dmAlpha);
    RooProdPdf RandomPionBkg("RandomPionBkg", "RandomPion Background", RooArgList(RandomPionM_JSU,RandomPionBkgDm));
    m_ws.import(RandomPionBkg);
  }
  return m_ws.pdf("RandomPionBkg");
}


RooAbsPdf* D2hhmumuModel::D2hhhhBackground(RooRealVar m, RooRealVar dm,
					   RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mDeltaJSU, RooRealVar mGammaJSU,
					   RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar dmDeltaJSU, RooRealVar dmGammaJSU
					   )
{

  if (m_ws.pdf("D2hhhhBkg") == 0) {
    RooJohnsonSU D2hhhhBackgroundM("D2hhhhBackgroundM", "D^{0} misidentified D2hhhhh", m, mMeanJSU, mWidthJSU,mDeltaJSU,mGammaJSU);
    RooJohnsonSU D2hhhhBackgroundDm_JSU("D2hhhhBackgroundDm_JSU", "misidentified D2hhhhh (#Deltam)", dm, dmMeanJSU, dmWidthJSU,dmDeltaJSU,dmGammaJSU);
    RooProdPdf D2hhhhBkg("D2hhhhBkg", "misidentified D2hhhhh background", RooArgList(D2hhhhBackgroundM,D2hhhhBackgroundDm_JSU));
    m_ws.import(D2hhhhBkg);
  }
  return m_ws.pdf("D2hhhhBkg");
}

void D2hhmumuModel::initializeDefaultYields() {

  //if no special yields are desired, a default set can be initialized here. the function Model return then a fully working model PDF. If blinded variables are to be used,
  //they can be given to the workspace later in the Fitter. See D2hhmumuFitter.cc for details. DONT use D2hhmumuModel in this case
  if(m_ws.var("nD2hhhhBkg") == 0){
    RooRealVar nD2hhhhBkg("nD2hhhhBkg", "misidentified D2hhmumu background events", 500, 0., 150000.);
    m_ws.import(nD2hhhhBkg);
  }   
  if(m_ws.var("nRandomPionBkg") == 0){
    RooRealVar nRandomPionBkg("nRandomPionBkg", "Random Pion Background Events", 200, 0., 10000000.);
    m_ws.import(nRandomPionBkg);
  }
  if(m_ws.var("nCombinatoricBkg") == 0){
    RooRealVar nCombinatoricBkg("nCombinatoricBkg", "Combinatoric Background Events", 1000, 0., 10000000.);
    m_ws.import(nCombinatoricBkg);  
  }
  if(m_ws.var("nSignal") == 0){
   RooRealVar nSignal("nSignal", "Signal Events", 3000, 0., 1000000.);
   m_ws.import(nSignal);
  }
   
}


/*
RooAbsPdf* D2hhmumuModel::D2hhhhBackground(RooRealVar m, RooRealVar dm,
			    RooRealVar mMeanGauss, RooRealVar mSigmaGauss,
			    RooRealVar dmMeanJSU, RooRealVar dmWidthJSU, RooRealVar dmDeltaJSU, RooRealVar dmGammaJSU
			    )
{

  if (m_ws.pdf("D2hhhhBkg") == 0) {
    RooGaussian D2hhhhBackgroundM("D2hhhhBackgroundM", "misidentified D2hhhhh D^{0} background",m ,mMeanGauss ,mSigmaGauss );
    RooJohnsonSU D2hhhhBackgroundDm_JSU("D2hhhhBackgroundDm_JSU", "misidentified D2hhhhh (#Deltam)", dm, dmMeanJSU, dmWidthJSU,dmDeltaJSU,dmGammaJSU);
    RooProdPdf D2hhhhBkg("D2hhhhBkg", "misidentified D2hhhhh background", RooArgList(D2hhhhBackgroundM,D2hhhhBackgroundDm_JSU));
    m_ws.import(D2hhhhBkg);
    RooRealVar nD2hhhhBkg("nD2hhhhBkg", "misidentified D2hhmumu background events", 500, 0., 100000.);
    m_ws.import(nD2hhhhBkg);
  }
  return m_ws.pdf("D2hhhhBkg");
}
*/


RooAbsPdf* D2hhmumuModel::D2hhhhRandomPionBackground(RooRealVar m, RooRealVar dm,
				      RooRealVar mMeanGauss, RooRealVar mSigmaGauss,
				      RooRealVar dmThreshold, RooRealVar dmAlpha
				      )
{
  if (m_ws.pdf("D2hhhhRandomPionBkg") == 0) {
    RooGaussian D2hhhhRandomPionBackgroundM("D2hhhhRandomPionBackgroundM", "misidentified D2hhhhh D^{0} background associated with random pion",m ,mMeanGauss ,mSigmaGauss );
    RooThreshold D2hhmumuRandomPionBkgDm("D2hhmumuRandomPionBkgDm", "misidentified D2hhhhh D^{0} with Random Pion Background (#Deltam)", dm, dmThreshold, dmAlpha);
    RooProdPdf D2hhhhRandomPionBkg("D2hhhhRandomPionBkg", "misidentified D2hhhhh background", RooArgList(D2hhhhRandomPionBackgroundM,D2hhmumuRandomPionBkgDm));
    m_ws.import(D2hhhhRandomPionBkg);
    RooRealVar nD2hhhhRandomPionBkg("nD2hhhhRandomPionBkg", "misidentified D2hhmumu background events with random pion", 100, 0., 150000.);
    m_ws.import(nD2hhhhRandomPionBkg);
  }
}


RooAbsPdf* D2hhmumuModel::Model(std::string components)
{
 
  initializeDefaultYields(); // here the default yields are initialized. So do not use this function if you want to give special yields, e.g. for blinded variables

  std::vector<std::string> comps;
  std::istringstream iss(components);
  copy(std::istream_iterator<std::string>(iss),
       std::istream_iterator<std::string>(),
       back_inserter(comps));


  if (m_ws.pdf("D2hhmumuModel") == 0) {
    RooArgList pdfList, coefList;
    for (std::vector<std::string>::iterator it = comps.begin();
         it != comps.end(); ++it) {
      if (m_ws.pdf(it->c_str()) == 0) {
	std::cout << "Error(D2hhmumuModel::Model): no component found with name " << *it << std::endl;
        return 0;
      }
      pdfList.add(*m_ws.pdf(it->c_str()));
      coefList.add(*m_ws.var((std::string("n")+*it).c_str()) );
     }

    pdfList.Print();
    coefList.Print();
    RooAddPdf Model("D2hhmumuModel", "D*^{+} #rightarrow D^{0}(#rightarrow hh #mu #mu)#pi^{+} Model", pdfList, coefList);
    m_ws.import(Model);
  } 
  return m_ws.pdf("D2hhmumuModel");
}
