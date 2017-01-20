#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "Tools.h"


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

#include "TMVA/MethodBase.h"
#include "TMVA/MethodCategory.h"


//training, testing and application for D2KKmumu
void Classification_D2KKmumu(int part);
void D2KKmumuCrosstraining();

//training, testing and application for D2pipimumu
void Classification_D2pipimumu(int part);
void D2pipimumuCrosstraining();

//BDT is applied to data, MC, data sideband and normalization channel
void Application_D2KKmumu(TString treeName, TString fileIn, TString fileOut, int part,bool isMC, bool isNormalizationMode, bool skipPID);
void D2KKmumuCrossapplication();

//BDT is applied to data, MC, data sideband and normalization channel
void Application_D2pipimumu(TString treeName, TString fileIn, TString fileOut, int part,bool isMC, bool isNormalizationMode, bool skipPID);
void D2pipimumuCrossapplication();

void CrossapplicationForEfficiencyStudies();
void CrossapplicationForEfficiencyStudiesNoTruthmatching();
void CrossapplicationForEfficiencyStudiesNoQ2Splitting();
