#include "D2hhmumuFitter.h"
#include "D2KKmumuReader.h"
#include "D2pipimumuReader.h"
#include "D2KpimumuReader.h"
#include "D2KKpipiReader.h"
#include "D2KpipipiReader.h"
#include <utility> 

#include "D2hhmumuFitter1D.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooStats/ModelConfig.h"
#include "RooUnblindPrecision.h"
#include "RooUnblindOffset.h"
#include "RooCategory.h"
#include "TBox.h"
#include "TLatex.h"
#include "RooRandom.h"
#include "RooMinuit.h"
#include "RooGenericPdf.h"


using namespace RooFit ;

//#include "dcastyle.cc"

//help functions
void createMCTuplesForEffStudies();//created a tuples for each q2 bin and channel, applied to the flagged MC sample
void createMCTuplesForEffStudiesNoTruthmatching(); //also no q2 splitting is done
void createDataTuplesForEffStudies();
void copyMCTuplesWithCut(TString fileIn,TString target, TString nameTree,TString cut) ;
//the following one is applied to the flagged MC sample, needed to get the reco and stripping efficiency
void createTuplesForRecoEfficiency();

void MC_PID_efficiency(double ghostProbCut,double hadronPID,double muonPID);
//MC BDT efficiency
void MC_BDT_efficiency(double ghostProbCut,double hadronPID,double muonPID, double BDT, bool applyPID, bool applyTrigger);

void MC_total_efficiency(double ghostProbCut,double hadronPID,double muonPID,double BDT);


//PID efficiency with PID calib
std::pair<double,double> evaluatePIDCalibEfficiency(TString fIn);
std::pair<double,double> getMCEfficiency(TString fIn,TString nameTrree, TString cut_sel,TString cut_norm );
std::pair<double,double> getMCEfficiencyFrom2Files(TString fIn1,TString fIn2,TString nameTrree, TString cut_sel,TString cut_norm );
std::pair<double,double> fitMC(TString file,TString treeName, TString cut, TString namePlot);

// muon PID efficiency studies
void splitMCtuplesInPt(bool cutNshared);
void splitMCtuplesInPtOfBothMuons(bool cutNshared); //pt>800 on both 
void splitMCtuplesPtCutOnSingleMuon(bool cutNshared,int index); //pt>800 cut only on one of the muons (chose by index 0 or 1)

void evaluateMCSingleMuonEfficiency(int muonIndex);
void evaluatePIDCalibMuonEfficiencyForPTBins(TString polarity);
void evaluateMCDoubleMuonEfficiency();
void evaluatePIDCalibSingleMuonEfficiency(TString polarity,int muonIndex,TString binningScheme);
void evaluatePIDCalibDoubleMuonEfficiency(TString polarity);
void drawLowPtCorrectedEfficiency(TString binningScheme);
void drawBothPolaritiesPIDCalibMuonEfficiencyForPTBins();
void drawBothPolaritiesPIDCalibDoubleMuonEfficiency();
void drawBothPolaritiesPIDCalibSingleMuonEfficiency(int muonIndex,TString binningScheme);
void computeGlobalLowPtEfficiency();
void plotPIDCalibMuonSampleAndSignalMC();
void checkMuonPIDFactorization();
void evaluatePIDEfficiencySystematicUncertainty();

//hadron PID
void evaluatePIDCalibHadronEfficiency(TString polarity,TString binningScheme);
void plotPIDCalibKaonSampleAndSignalMC();
void plotPIDCalibPionSampleAndSignalMC();
void drawBothPolaritiesPIDCalibHadronEfficiency(TString binningScheme);

void drawPtSpectra(int muonIndex) ;
void drawTotalPIDEfficiency(TString binningScheme,double dEffMu_sys,double dEffHad_sys);
void performFullPIDEfficiencyStudy();
void compareHadonAndMuonEfficienciesWithDifferentBinning();

//generator level cuts Efficiency
void drawGenLevelCutEfficiencies();

//stripping and reco
void addDiMuonMassBranch(TChain *chain, TString fOut, bool isReco);
void drawRecoAndStrippingEfficiency(int cat);
void drawRecoAndStrippingEfficiencyWithFit(bool useFit,TString additionalCut);
void compareDifferentRecoAndStrippingEfficiencies();

//multiple candidates
void writeMultipleCanidatesToFile(bool withGhosts);
void chose_multiple_events(const char *input,TString channel,bool withGhosts);
void addMultipleCandidateBranch(TChain *chain, TString fOut,const char *f_chosenCand);
void createTuplesWithSelectedMultipleCand();
void writeAllTrueCanidatesToFile();


///TAG and Probe method to get L0 trigger efficiecny, hadron PID and BDT cut hardcoded
void calibrateTagAndProbeL0EfficiencyWithMC();
void mergeDs2PhiPiTuples();
void calibrateTagAndProbeL0EfficiencyWithDsPhiPi();
void applyTagAndProbeL0Efficiency(bool useData); //if use data is set to false a MC closure test will be performed
//void applyTagAndProbeL0EfficiencyToyStudy(bool useData, int nToys);

//MC bases HLT level efficiency (HLT1 and 2)
void MC_HltTrigger_efficiency(double ghostProbCut,double hadronPID,double muonPID,double BDT,bool applyPID, bool applyBDT, TString triggerLevel);
void MC_Combined_Hlt_BDT_efficiency(double ghostProbCut,double hadronPID,double muonPID,double BDT,bool applyPID,bool withGhosts);



//Old version where the efficienc was tried to be computed with TISTOS methon. No satisfying results
void MC_L0Trigger_efficiency(double ghostProbCut,double hadronPID,double muonPID, double BDT, bool applyPID,bool applyTrigger, TString cut_Trigger_TIS);


void drawTotalEfficiency();


