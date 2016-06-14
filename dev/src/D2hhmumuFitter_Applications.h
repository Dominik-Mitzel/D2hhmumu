// Fill workspaces do the fits 
#include "D2hhmumuFitter1D.h"
#include "D2hhmumuFitter.h"

class D2hhmumuFitter_Applications {

 public:
  D2hhmumuFitter_Applications(TString m_kind, TString m_year);
  ~D2hhmumuFitter_Applications();

  TString year;
  TString kind;
  TString targetFolder;
  TString q2RangeNormalizationMode;

  TString pathToSignalData;
  TString pathToNormData;
  TString pathToInvData;
  TString pathToSidebandData;
  TString pathToKpipipiData;
  TString pathToKKpipiData;
  TString pathToSignalMC;
  TString pathToKpimumuMC;

  std::vector<TString> q2Ranges;

  void setDefaultPathToData(TString m_kind);
  void setQ2Ranges(TString m_kind);
  void saveModelConfig(TString dataCut,TString nomalizationCut,TString misIDCut,TString q2Range);
  void saveAllModelConfig(TString dataCut,TString nomalizationCut,TString misIDCut);

  void compare_1D_and_2D_fit(TString dataCut,TString nomalizationCut,TString misIDCut);
  void compare_misID_shapes();
  void compare_misID_shapes_2D();
  void ExtractExpectedLimit(TString fIn);
  void ExtractAllExpectedLimit(TString dataCut);
  void runFull1DFits(TString dataCut,TString normalizationCut,TString misIDCut,TString q2Range);
  void runAllFull1DFits(TString dataCut,TString normalizationCut,TString misIDCut);  
  void performAllToyStudies();
  };
  
