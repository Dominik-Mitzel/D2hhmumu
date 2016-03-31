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

//BDT is applied to data, MC, data sideband and normalization channel
void Application_D2KKmumu(TString treeName, TString fileIn, TString fileOut, int part);
void D2KKmumuCrossapplication();

