#include "optimizeSelection.h"
#include "D2hhmumuFitter.h"
#include "TEventList.h"
#include "TPaletteAxis.h"
#include "D2KKmumuReader.h"
#include "D2pipimumuReader.h"
#include "D2KpimumuReader.h"
#include "D2KKpipiReader.h"
#include "D2KpipipiReader.h"
#include "TProfile.h"
#include "D2hhmumuFitter1D.h"
#include "TGenPhaseSpace.h"

///////////////////////////////////////////////////////
//
//Class to calcululate efficiency for signal, normalization
//mode and misID D2hhhh decays.
//Create an object of EfficiencyCalculator("D2KKmumu") ( EfficiencyCalculator("D2pipimumu")) and 
// - D2KKmumu(D2pipimumu) is set as signal
// - D2Kpimumu is normlization 
// - D2KKpipi (D2pipipi) is set as signal misID
// - D2Kpipipi is set as normalization misiD 
///////////////////////////////////////////////////////


class EfficiencyCalculator {

 public:

  EfficiencyCalculator(TString kind);
  ~EfficiencyCalculator();

  double genLevelCutEfficiency_KKmumu;
  double genLevelCutEfficiency_Kpimumu;
  double genLevelCutEfficiency_pipimumu;

  
  TString m_kind;
  TString pathToFiles;
  
  TChain* tree_recoSignal;
  TChain* tree_genSignal;

  TChain* tree_recoNorm;
  TChain* tree_genNorm;

  TChain* tree_recoSigMisID;
  TChain* tree_genSigMisID;
 
  TChain* tree_recoNormMisID;
  TChain* tree_genNormMisID;
  
  // selection Cut, splitting of dataset (as used in signal optimisation) and q2 range gan be specified
  // for normalization channel, the dimuon mass range is restricted to 675-875 MeV by default, however, can be specified to what ever is needed
  double getMCSignalEfficiency(TString selectionCut, TString splitting, TString q2Range);  
  double getMCNormalizationEfficiency(TString selectionCut, TString splitting, TString q2Range); 
  double getMCRelativeSigToNormEfficiency(TString selectionCut, TString splitting, TString q2RangeSig,TString q2RangeNorm);
  
  double getMCSignalMisIDEfficiency(TString selectionCut, TString splitting, TString q2Range);
  double getMCNormalizationMisIDEfficiency(TString selectionCut, TString splitting, TString q2Range); 
  double getMCRelativeSigToNormMisIDEfficiency(TString selectionCut, TString splitting, TString q2RangeSig,TString q2RangeNorm);
  double getMisIDFractionQ2Range(TString q2Range);

  double getMCSignalEfficiencyError(TString selectionCut, TString splitting, TString q2Range);
  double getMCNormalizationEfficiencyError(TString selectionCut, TString splitting, TString q2Range);
  double getMCRelativeSigToNormEfficiencyError(TString selectionCut, TString splitting, TString q2RangeSig,TString q2RangeNorm);

};
  


  
  
