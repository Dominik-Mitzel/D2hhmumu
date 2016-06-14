#include "D2hhmumuFitter.h"
#include "D2KKmumuReader.h"
#include "D2pipimumuReader.h"
#include "D2KpimumuReader.h"
#include "D2KKpipiReader.h"
#include "D2KpipipiReader.h"

/////////////////////////////////////////////////////////
//
// in optimizeSelection, a lot of functions used in the 
// selection are defined:
//  - the D2hh(`)ll/hh() functions merge the raw tuples, add relevant new branches and apply a very 
//    lose preselection (trigger selection, cut on ghost probability). Moreover, the data set is split
//    by even and odd event number which is needed fot the MVA selection. Also a tuple containing of sideband data
//    is added to all the tuples.
//  - optimizeSelection performs a scan of the FOM is the 2D space PID/BDT to decide on an optimal cut value
//    Different versions are available: 1D fitter, 2D fitter, binned in dimuon mass
//  - draw_BDT_crosschecks is used to check possible correlation of the BDT with masses. The BDT is designed to have
//    flat efficincy in mD, deltaM, mHH and mMuMu
//  - check_peakingBackground is reconstrcting the data under several mass hypotheses to check for possible sources of peaking bkg
////////////////////////////////////////////////////////



void optimizeSelection(TString cutNtracks);
void optimizeSelection2D(TString cutNtracks);
void optimizeSelectionInBins(TString cutNtracks,TString q2Range);

double getMCSignalEfficiency(TString cut);
double EffD2KKpipiToEffD2Kpipipi(TString cut);

void D2pipimumuMC();
void D2pipimumuData();
void D2KpimumuMC();
void D2KpimumuData();
void D2KKmumuMC(); 
void D2KKmumuData(); 
void D2KKpipiMC();
void D2KpipipiMC();
void D2KKpipiData();
void D2KpipipiData();

void draw_BDT_crosschecks_forVariable(TString variable, int xLow, int xHigh, double BDTcut);
void draw_BDT_crosschecks();

void check_peakingBackground(TString kind,bool PIDCut);
void createGeneratorLevelMCTuple(TString kind);
