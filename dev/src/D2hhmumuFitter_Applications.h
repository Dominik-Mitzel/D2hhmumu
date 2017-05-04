// Fill workspaces do the fits 
#include "D2hhmumuFitter1D.h"
#include "D2hhmumuFitter.h"

class D2hhmumuFitter_Applications {

 public:
  D2hhmumuFitter_Applications(TString m_kind, TString m_year);
  ~D2hhmumuFitter_Applications();

  ///////////////////////////////////////////////////
  //
  // class to manage everything related to the fitter
  //
  ///////////////////////////////////////////////////


  TString year;
  TString kind;
  TString targetFolder;
  TString q2RangeNormalizationMode;

  TString pathToSignalData;
  TString pathToNormData;
  TString pathToNormMC;
  TString pathToInvData;
  TString pathToSidebandData;
  TString pathToKpipipiData;
  TString pathToKKpipiData;
  TString pathToSignalMC;
  TString pathToKpimumuMC;

  std::vector<TString> q2Ranges;

  void setDefaultPathToData(TString m_kind);
  void setQ2Ranges(TString m_kind);
  std::vector<double> rangespipi_low = {200,525,565,950,1100};
  std::vector<double> rangespipi_high = {525,565,950,1100,1600};
  std::vector<double> rangesKK_low = {200,525,565};
  std::vector<double> rangesKK_high = {525,565,900};//!!//THIS IS THE NEW BINNING 

  void saveModelConfig(TString dataCut,TString misIDCut);
  void compare_1D_and_2D_fit(TString dataCut,TString nomalizationCut,TString misIDCut);
  void ExtractExpectedLimit();
  void runFull1DFits(TString dataCut,TString misIDCut);
  void runFullResonant1DFits(TString dataCut,TString misIDCut);
  void constrainCombBkgShapes(TString dataCut,TString misIDCut);
  void studyNormalizationFits(TString dataCut,TString misIDCut, bool doNSharedCut); 

  void drawMisIDShapes(TString cut);
  void drawMCSignalShapes(TString cut);

  //for systematics studies
  void studyResolutionScale(TString cut);
  void performAllToyStudies();
  void studyCombBkgShape();    
  void compare_misID_shapes(TString dataCut);
  void draw_misID_shapes_singlePlot() ;
  
  //obsolet
  void studyTriggerEfficiencyNormalizationMode();
  void studyTISEfficiencyNormalizationModeBinned();

};
  
