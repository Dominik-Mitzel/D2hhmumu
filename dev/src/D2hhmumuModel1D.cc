#include "D2hhmumuModel1D.h"           
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
D2hhmumuModel1D::D2hhmumuModel1D()
{
}

D2hhmumuModel1D::~D2hhmumuModel1D()
{
}

RooAbsPdf* D2hhmumuModel1D::Signal(RooRealVar m, 
				   RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mNuJSU, RooRealVar mTauJSU,
				   RooRealVar ResolutionScale, RooRealVar globalShift
				   )
{
  RooFormulaVar mean("mean","@0*@1",RooArgList(mMeanJSU,globalShift));
  RooFormulaVar width("width","@0*@1",RooArgList(mWidthJSU,ResolutionScale));
  if (m_ws.pdf("Signal") == 0) {
    //RooJohnsonSU Signal("Signal", "Signal D^{0} JSU", m, mMeanJSU,mWidthJSU,mNuJSU,mTauJSU);
    RooJohnsonSU Signal("Signal", "Signal D^{0} JSU", m, mean,width,mNuJSU,mTauJSU);
    m_ws.import(Signal);
  }
  return m_ws.pdf("Signal");
}

/* //alternative signal PDF 
RooAbsPdf* D2hhmumuModel1D::Signal(RooRealVar m,
				   RooRealVar mean, RooRealVar sigma, RooRealVar alphaR, RooRealVar alphaL, RooRealVar nR, RooRealVar nL  )
{
  if (m_ws.pdf("Signal") == 0) {
    RooCBShape CBleft("CBleft", "Signal D^{0} CBleft", m, mean, sigma,alphaL,nL);
    RooCBShape CBright("CBright", "Signal D^{0} CBright", m, mean, sigma,alphaR,nR);
    RooRealVar fraction("fraction","fraction",0.5);
    RooAddPdf  Signal ("Signal", "D*^{+} #rightarrow D^{0}(#rightarrow hh #mu #mu)#pi^{+} Signal", RooArgSet(CBleft,CBright),fraction);
    m_ws.import(Signal);
  }
  return m_ws.pdf("Signal");
}
*/


RooAbsPdf* D2hhmumuModel1D::CombinatoricBackground(RooRealVar m, 
				      RooRealVar mChebyA,RooRealVar mChebyB 
						 )			
{
  if (m_ws.pdf("CombinatoricBkg") == 0) {
    RooChebychev CombinatoricBkg("CombinatoricBkg", "Combinatoric Background (M)", m, RooArgList(mChebyA,mChebyB));
    m_ws.import(CombinatoricBkg);
  }
  return m_ws.pdf("CombinatoricBkg");
}


RooAbsPdf* D2hhmumuModel1D::CombinatoricExponentialBackground(RooRealVar m, 
				      RooRealVar mExpoLambda 
						 )			
{
  if (m_ws.pdf("CombinatoricExpoBkg") == 0) {
    RooExponential CombinatoricExpoBkg("CombinatoricExpoBkg", "Combinatoric Background (M) Exponential", m, mExpoLambda);
    m_ws.import(CombinatoricExpoBkg);
  }
  return m_ws.pdf("CombinatoricExpoBkg");
}




RooAbsPdf* D2hhmumuModel1D::D2hhhhBackground(RooRealVar m, 
					   RooRealVar mMeanJSU, RooRealVar mWidthJSU, RooRealVar mDeltaJSU, RooRealVar mGammaJSU
								   )
{

  if (m_ws.pdf("D2hhhhBkg") == 0) {
    RooJohnsonSU D2hhhhBkg("D2hhhhBkg", "D^{0} misidentified D2hhhhh", m, mMeanJSU, mWidthJSU,mDeltaJSU,mGammaJSU);
    m_ws.import(D2hhhhBkg);
  }
  return m_ws.pdf("D2hhhhBkg");
}

 //alternative D2hhhhBkg models
RooAbsPdf* D2hhmumuModel1D::D2hhhhDoubleCBBackground(RooRealVar m,
                                   RooRealVar mean, RooRealVar sigma, RooRealVar alphaR, RooRealVar alphaL, RooRealVar nR, RooRealVar nL
                                   )
{
  if (m_ws.pdf("D2hhhhDoubleCBBkg") == 0) {
    RooCBShape CBleftBkg("CBleftBkg", "D2hhhhbkg D^{0} CBleft", m, mean, sigma,alphaL,nL);
    RooCBShape CBrightBkg("CBrightBkg", "D2hhhhbkg D^{0} CBright", m, mean, sigma,alphaR,nR);
    RooRealVar fractionBkg("fraction","fraction",0.5);
    RooAddPdf D2hhhhDoubleCBBkg("D2hhhhDoubleCBBkg", "D*^{+} #rightarrow D^{0}(#rightarrow hh #mu #mu)#pi^{+} Signal", RooArgSet(CBleftBkg,CBrightBkg),fractionBkg);
    m_ws.import(D2hhhhDoubleCBBkg);
  }
  return m_ws.pdf("D2hhhhDoubleCBBkg");
}
 

RooAbsPdf* D2hhmumuModel1D::D2hhhhSingleCBBackground(RooRealVar m,
				   RooRealVar mean, RooRealVar sigma, RooRealVar alphaL, RooRealVar nL
				   )
{
  if (m_ws.pdf("D2hhhhSingleCBBkg") == 0) {
    RooCBShape D2hhhhSingleCBBkg("D2hhhhSingleCBBkg", "Signal D^{0} CBleft", m, mean, sigma,alphaL,nL);
    m_ws.import(D2hhhhSingleCBBkg);
  }
  return m_ws.pdf("D2hhhhSingleCBBkg");
}



void D2hhmumuModel1D::initializeDefaultYields() {

  //if no special yields are desired, a default set can be initialized here. the function Model return then a fully working model PDF. If blinded variables are to be used,
  //they can be given to the workspace later in the Fitter. See D2hhmumuFitter.cc for details. DONT use D2hhmumuModel in this case
  if(m_ws.var("nD2hhhhBkg") == 0){
    RooRealVar nD2hhhhBkg("nD2hhhhBkg", "misidentified D2hhmumu background events", 500, 0., 2500000.);
    m_ws.import(nD2hhhhBkg);
  }   
  if(m_ws.var("nD2hhhhDoubleCBBkg") == 0){
    RooRealVar nD2hhhhDoubleCBBkg("nD2hhhhDoubleCBBkg", "misidentified D2hhmumu background events", 500, 0., 2500000.);
    m_ws.import(nD2hhhhDoubleCBBkg);
  }   
  if(m_ws.var("nD2hhhhSingleCBBkg") == 0){
    RooRealVar nD2hhhhSingleCBBkg("nD2hhhhSingleCBBkg", "misidentified D2hhmumu background events", 500, 0., 2500000.);
    m_ws.import(nD2hhhhSingleCBBkg);
  }   
  if(m_ws.var("nCombinatoricBkg") == 0){
    RooRealVar nCombinatoricBkg("nCombinatoricBkg", "Combinatoric Background Events", 1000, 0., 10000000.);
    m_ws.import(nCombinatoricBkg);  
  }
  if(m_ws.var("nCombinatoricExpoBkg") == 0){
    RooRealVar nCombinatoricExpoBkg("nCombinatoricExpoBkg", "Combinatoric Background Events", 1000, 0., 10000000.);
    m_ws.import(nCombinatoricExpoBkg);  
  }
  if(m_ws.var("nSignal") == 0){
   RooRealVar nSignal("nSignal", "Signal Events", 3000, 0., 1000000.);
   m_ws.import(nSignal);
  }
   
}




RooAbsPdf* D2hhmumuModel1D::Model(std::string components)
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
