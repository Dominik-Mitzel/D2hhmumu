#include <TChain.h>
#include <iostream>
#include <TString.h>
#include "LumiTupleReader.C"
#include "D2KKmumuReader.h"
#include "D2pipimumuReader.h"
#include "D2KpimumuReader.h"
#include "D2KKpipiReader.h"
#include "D2KpipipiReader.h"
#include "D2pipipipiReader.h"
#include "D2hhmumuFitter.h"
#include "sWeights.h"
#include "TMVA_applications.h"
#include "optimizeSelection.h"
#include "D2hhmumuFitter1D.h"
#include "D2hhmumuFitter_Applications.h"
#include "EfficiencyCalculator.h"
#include "EfficiencyStudies.h"

using namespace std;

using namespace RooFit ;
using namespace RooStats;


void checkLuminosity(){

  TChain* LumiTree_D2hhmumu = new TChain("GetIntegratedLuminosity/LumiTuple");

  for (int i=0; i<1500; ++i) {

    LumiTree_D2hhmumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/2012/magUp/2012Data_D2hhmumu_%i.root",i));
    LumiTree_D2hhmumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/2012/magDw/2012Data_D2hhmumu_%i.root",i));
 
  }

  LumiTupleReader* KK_TupleReader = new LumiTupleReader(LumiTree_D2hhmumu);
  double totalLumiKK = KK_TupleReader->getTotalLumi();


  std::cout<<"Total Luminosity KK: "<<totalLumiKK<<std::endl;

}


int main() {


  time_t startTime = time(0);


  ////////////////////////////////////////
  //                                    //
  //  create preselected samples        //
  //  split by run numbers              //
  //  create TMVA training samples      //
  //  defined in optimizeSelection.cc/h //
  ////////////////////////////////////////

  
  //signal and norm data
  //D2pipimumuMC();      
  //D2pipimumuData();                                                                                                                

  /*
  D2KpimumuData();
  D2KpimumuMC();        

  D2KKmumuMC();                                                                                                                     
  D2KKmumuData();

  
  //misID
  D2pipipipiData();
  D2pipipipiMC();

  D2KKpipiMC();                                                                                   
  D2KpipipiMC();                         

  */
  //D2pipipipiRandomizedData();
  //D2KpipipiRandomizedData();

  //D2KKpipiData();                                                                                                               
  //D2KpipipiData();                                                                                                               
  

  //create the tuples binned in defined D_DiMuon mass binning scheme
  //createTuplesForRecoEfficiency();
  

  //////////////////////////////////////
  //                                  //
  //  TMVA training and application   //
  //  training and application defined//
  //  in TMV_Aapplications.cc/.h      //
  //////////////////////////////////////


  //TMVA
  //D2KKmumuCrosstraining();
  //D2pipimumuCrosstraining();  
  
  //D2KKmumuCrossapplication();
  //D2pipimumuCrossapplication();
  
  //apply also the BDT for the efficiency study tuples
  //CrossapplicationForEfficiencyStudies();
  //createMCTuplesForEffStudiesNoTruthmatching();
  //CrossapplicationForEfficiencyStudiesNoTruthmatching();


  
   ///////////////////////////////
   //
   //    Efficiency Studies
   //
   //////////////////////////////

  //plotPIDCalibMuonSampleAndSignalMC();
  //plotPIDCalibKaonSampleAndSignalMC();
  //plotPIDCalibPionSampleAndSignalMC();
  
  
  // GEN LEVEL EFFICIENCIES
  //Draw the generator level efficiency. Taken from privatey produced gauss sample, where a TTree was saved before and after the cut tool
  //Results saved in "/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/GenLevelCutEfficiency/generatorLevelCutEfficiency.root"

  //drawGenLevelCutEfficiencies();
  

  //STRIPPING AND RECO AND PRESELECTION EFFICIENCY
  //Different ways to get the combined reco and stripping eff. are tried. Especially different treatments of multiple candidates in MC

  //here the different method are listed. the first argument is set true if a fit is ti be used to get the signal yield. second argument is cut (on top of the muon Nshared and delta M cut )
  //results can be found in ../img/EfficiencyStudies/strippingEfficiency/RecoAndStrippingEfficiency

  //drawRecoAndStrippingEfficiencyWithFit(true,"(Dst_BKGCAT<11||Dst_BKGCAT==60)");
  //drawRecoAndStrippingEfficiencyWithFit(true,"isSelectedMultipleCandidate");
  //drawRecoAndStrippingEfficiencyWithFit(false,"Dst_BKGCAT<11"); 
  //drawRecoAndStrippingEfficiencyWithFit(false,"Dst_BKGCAT<11&&isSelectedMultipleCandidate"); 
  //drawRecoAndStrippingEfficiencyWithFit(true,"isMatchedCandidate");
  //drawRecoAndStrippingEfficiencyWithFit(true,"isMatchedCandidate&&isSelectedMultipleCandidate");
  //compareDifferentRecoAndStrippingEfficiencies();

  ////////////////////////////////////////////////
  //multiple Candidate MC
  //
  //
  ////////////////////////////////////////////////

  //writeMultipleCanidatesToFile(true);
  //writeMultipleCanidatesToFile(false);

  //chose_multiple_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_200.0_1600.0_D2pipimumuBDT.root","D2pipimumu",true);
  //chose_multiple_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_200.0_1600.0_D2KKmumuBDT.root","D2KKmumu", true);
  //chose_multiple_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2pipimumuBDT.root","D2Kpimumu_D2pipimumuBDT",true);
  //chose_multiple_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2KKmumuBDT.root","D2Kpimumu_D2KKmumuBDT",true);

  //chose_multiple_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_200.0_1600.0_D2pipimumuBDT.root","D2pipimumu",false);
  //chose_multiple_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_200.0_1600.0_D2KKmumuBDT.root","D2KKmumu", false);
  //chose_multiple_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2KKmumuBDT.root","D2Kpimumu_D2KKmumuBDT",false);
  //chose_multiple_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2pipimumuBDT.root","D2Kpimumu_D2pipimumuBDT",false);

  //createTuplesWithSelectedMultipleCand();
  
  //writeAllTrueCanidatesToFile(); //write all candidates of the true MC tuple to be able to perform a matching later with the reconstructed candidates
  
  //////////////////////////////////////////
  //
  // PID EFFICIENCIES
  //
  //////////////////////////////////////////

  //this is done using the PID CALIB package 
  //a description of the procedure can be found in EfficienyStudy.cc performFullPIDEfficiencyStudy()
  //performFullPIDEfficiencyStudy();
  //evaluatePIDEfficiencySystematicUncertainty();

  //checkMuonPIDFactorization();

  //////////////////////////////////////////
  //
  //HLT efficiency
  //
  //////////////////////////////////////////

  //MC_HltTrigger_efficiency(0.5,0.2,0.5,0.,true,true,"Hlt1");
  //MC_HltTrigger_efficiency(0.5,0.2,0.5,0.,true,true,"Hlt2");
  //MC_HltTrigger_efficiency(0.5,0.2,0.5,0.,true,false,"FullHlt");
  
  //BDT EFFICIENCY
  //MC_BDT_efficiency(0.5,0.2,0.5,0.4,true,true);

  ////////////////////////////////////////////
  //
  //LEVEL 0 WITH TAGPROBE
  //
  /////////////////////////////////////////////

  //mergeDs2PhiPiTuples();
  //calibrateTagAndProbeL0EfficiencyWithMC();
  //calibrateTagAndProbeL0EfficiencyWithDsPhiPi();
  //applyTagAndProbeL0Efficiency(false);
  //applyTagAndProbeL0Efficiency(true);
  //applyTagAndProbeL0EfficiencyToyStudy(true,100); not yet implemented


  //MC_Combined_Hlt_BDT_efficiency(0.5,0.2,0.5,0.4,true,false);
  //MC_Combined_Hlt_BDT_efficiency(0.5,0.2,0.5,0.4,true,true);
  
 
  //createDataTuplesForEffStudies();

  //drawTotalEfficiency();


  //norm mode and D->hhpipi sWeights
  
  /*
  D2hhmumuFitter1D myFitter1D_1;
  myFitter1D_1.addNormalizationSWeightsHadronicChannel("eventNumber%5==0","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/temp/D2KKpipi_D2KKmumuBDT.root","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/temp/D2KKpipi_D2KKmumuBDT_sWeights.root","drawingMacros/sWeights_D2KKpipi.eps");

  D2hhmumuFitter1D myFitter1D_2;
  myFitter1D_2.addNormalizationSWeightsHadronicChannel("eventNumber%500==0","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/temp/D2Kpipipi_D2KKmumuBDT.root","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/temp/D2Kpipipi_D2KKmumuBDT_sWeights.root","drawingMacros/sWeights_D2KKpipi.eps");
  
  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2KKmumuBDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.addNormalizationSWeights("mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&Dst_DTF_D0_M<1950&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875","mu0_ProbNNmu>0.5","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT_sWeights.root");
  */
 



  //MC_L0Trigger_efficiency(0.5,0.2,0.5,0.,true,false,"D_L0Global_TIS");
  //MC_L0Trigger_efficiency(0.5,0.2,0.5,0.,false,false,"D_L0MuonDecision_TIS");
  //MC_L0Trigger_efficiency(0.5,0.2,0.5,0.,true,false,"D_L0MuonDecision_TIS");
  //MC_L0Trigger_efficiency(0.5,0.2,0.5,0.,true,false,"D_L0HadronDecision_TIS");
  //MC_L0Trigger_efficiency(0.5,0.2,0.5,0.,true,false,"(D_L0MuonDecision_TIS||D_L0DiMuonDecision_TIS||D_L0HadronDecision_TIS)");
  //MC_L0Trigger_efficiency(0.5,0.2,0.5,0.,true,true,"(D_L0MuonDecision_TIS||D_L0DiMuonDecision_TIS||D_L0HadronDecision_TIS)");


  //MC_L0Trigger_efficiency(0.5,0.2,0.5,0.,true,true,"(D_L0MuonDecision_TIS&&D_L0DiMuonDecision_TIS&&D_L0HadronDecision_TIS)");

  //drawPtSpectra(0);
  //drawPtSpectra(1);
  //createMCTuplesForEffStudies();  
  //CrossapplicationForEfficiencyStudies(); //apply BDT to the tuples splitted in q2 ranges
  //plotPIDCalibSampleAndSignalMC("");
  //splitMCtuplesInPt(true);//binning is defined in EfficiencyStudies.cc, optionally apply nShared==o cut to muons
  //splitMCtuplesInPtOfBothMuons(true);//binning is defined in EfficiencyStudies.cc, optionally apply nShared==o cut to muons
  //splitMCtuplesPtCutOnSingleMuon(true,1);

  //drawBothPolaritiesPIDCalibDoubleMuonEfficiency();
  //evaluatePIDCalibSingleMuonEfficiency("magDw",0);
  //evaluatePIDCalibSingleMuonEfficiency("magUp",0);
  //evaluatePIDCalibSingleMuonEfficiency("magDw",1);
  //evaluatePIDCalibDoubleMuonEfficiency("magDw");
  //evaluatePIDCalibSingleMuonEfficiency("magUp",1);
  //evaluatePIDCalibMuonEfficiencyForPTBins("magDw");
  //evaluatePIDCalibMuonEfficiencyForPTBins("magUp");
  //evaluatePIDCalibMuonEfficiencyForPTBins("magUp");
  //drawBothPolaritiesPIDCalibMuonEfficiencyForPTBins();
  //evaluateMCSingleMuonEfficiency(0);
  //evaluateMCSingleMuonEfficiency(1);
  //evaluateMCDoubleMuonEfficiency();
  //drawBothPolaritiesPIDCalibSingleMuonEfficiency(1);
  //drawBothPolaritiesPIDCalibSingleMuonEfficiency(0);
  //drawLowPtCorrectedEfficiency();
  //checkMuonPIDFactorization();
  //drawBothPolaritiesPIDCalibHadronEfficiency();
  //drawBothPolaritiesPIDCalibMuonEfficiencyForPTBins();
  //evaluateMCSingleMuonEfficiency();
  //evaluateMCDoubleMuonEfficiency();
  //evaluatePIDCalibHadronEfficiency("magDw");
  //evaluatePIDCalibHadronEfficiency("magUp");
  //plotPIDCalibKaonSampleAndSignalMC();
  //plotPIDCalibPionSampleAndSignalMC();
  //getMCEfficiency("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/PTBins/MC_D2KKmumu_0.0_525.0_800.0_1200.0_D2KKmumuBDT_magUp.root","BDT_Tree","mu0_ProbNNmu>0.5","D_DiMuon_Mass>0");
  //getMCEfficiencyError("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/PTBins/MC_D2KKmumu_0.0_525.0_800.0_1200.0_D2KKmumuBDT_magUp.root","BDT_Tree","mu0_ProbNNmu>0.5","D_DiMuon_Mass>0");

  
  //createGeneratorLevelMCTuple("D2KKmumu");
  //createGeneratorLevelMCTuple("D2Kpimumu");
  //createGeneratorLevelMCTuple("D2KKpipi");
  //createGeneratorLevelMCTuple("D2Kpipipi");
  //createGeneratorLevelMCTuple("D2pipimumu");
  //createGeneratorLevelMCTuple("D2pipipipi");  

  //EfficiencyCalculator myEfficiencies("D2KKmumu");
  //myEfficiencies.getMCSignalEfficiency("BDT>0.&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","","D_DiMuon_Mass<565");
  //myEfficiencies.getMCSignalEfficiencyError("BDT>0.&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","","D_DiMuon_Mass<565");

  //myEfficiencies.getMCNormalizationEfficiency("BDT>0.&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","","D_DiMuon_Mass>675&&D_DiMuon_Mass<875 ");
  //myEfficiencies.getMCNormalizationEfficiencyError("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","","D_DiMuon_Mass>675&&D_DiMuon_Mass<875 ");
  //myEfficiencies.getMCRelativeSigToNormEfficiency("BDT>0.&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","","D_DiMuon_Mass<565","D_DiMuon_Mass>675&&D_DiMuon_Mass<875 ");
  //myEfficiencies.getMCRelativeSigToNormEfficiencyError("BDT>0.&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","","D_DiMuon_Mass<565","D_DiMuon_Mass>675&&D_DiMuon_Mass<875 ");



  ////////////////////////////////////
  // data MC comparison with        //
  //   sweights                     //                                                                                                                                
  //                                //
  ////////////////////////////////////

  //sweights part                     
  //calculateSweights("../Data_MC_Comparison/D2Kpimumu_BDT_selected.root");
  //data_MC_comparison();
  // checkLuminosity();      


  ////////////////////////////////////////
  //                                    //
  // Fitter                             //
  // Code for fit model: D2hhmumuModel  //
  // fitting is done with D2hhmumuFitter//
  ////////////////////////////////////////

 
  
  /*         
  D2hhmumuFitter1D* myFitter1D = new D2hhmumuFitter1D();
  myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");
  myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root");
  myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_new.root");  
  myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT.root");  
  myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root");  
  myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root");
  
  myFitter1D->fit_MC("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,"test.eps","m(KK#mu#mu)","950MeV<m(#mu#mu)<1100MeV");
  //myFitter1D->fit_normalization_MC("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5",true,"test1.eps");
  
  //myFitter1D->fit_Kpipipi_misID("BDT>0.4&&nTracks%15==0",true,"test2.eps");
  myFitter1D->fit_HHpipi_misID("BDT>0.4&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,"test3.eps","m_{#pi#pi#mu#mu}(#pi#pi#pi#pi)","950MeV<m(#pi#pi)<1100MeV");
  
  //myFitter1D->fit_normalization_Data("BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","test4.eps");
  //myFitter1D->fit_resonant_Data("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100","test5.eps","m(#pi#pi#mu#mu)","950MeV<m(#mu#mu)<1100MeV");
  myFitter1D->fit_Data("D2pipimumu","BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass>950&D_DiMuon_Mass<1100","test5.eps");
  
  */
  
  D2hhmumuFitter_Applications myApplication("D2KKmumu","2012");
  //myApplication.drawMisIDShapes("BDT>0.4&&mu1_ProbNNmu>0.5");
  myApplication.drawMCSignalShapes("BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5");
  //myApplication.studyResolutionScale("BDT>0.4&&mu0_ProbNNmu<0.5&&mu1_ProbNNmu<0.5");
  //myApplication.runFull1DFits("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication.saveModelConfig("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5"); 
  //myApplication.ExtractExpectedLimit();

  D2hhmumuFitter_Applications myApplication2("D2pipimumu","2012");
  //myApplication2.drawMisIDShapes("BDT>0.4&&mu1_ProbNNmu>0.5");
  myApplication2.drawMCSignalShapes("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5");
  //myApplication2.studyResolutionScale("BDT>0.4&&mu0_ProbNNmu<0.5&&mu1_ProbNNmu<0.5");
  //myApplication2.performAllToyStudies();
  //myApplication2.runFull1DFits("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication2.compare_misID_shapes("BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5");
  //myApplication2.saveModelConfig("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication2.ExtractExpectedLimit();

  /*

  // myFitter1D->ResolutionScale.setVal(1.1007);
   myFitter1D->ResolutionScale.setVal(1.0918);
   myFitter1D->ResolutionScale.setConstant();
  //myFitter1D->globalShift.setVal(1.0004);
  myFitter1D->globalShift.setVal(1.0003);
  myFitter1D->globalShift.setConstant();
  myFitter1D->D0_M_chebyB.setVal(0); 
  myFitter1D->D0_M_chebyB.setConstant();
  myFitter1D->D0_M_chebyA.setVal(-4.8960e-01);
  myFitter1D->D0_M_chebyA.setConstant();
  myFitter1D->fit_normalization_Data("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","test4.eps");
  */
  //myFitter1D->fit_HHpipi("BDT>0.6&&D_DiMuon_Mass>565&&eventNumber<5e5","test2.eps");
  //myFitter1D->fit_HHpipi_misID("BDT>0.6&&D_DiMuon_Mass>565&&mu1_ProbNNmu>0.5",true,"test2.eps");
  //myFitter1D->fit_Data("BDT>0.6&&D_DiMuon_Mass>565&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu<0.5","test3.eps");
    
  //myApplication.runFull1DFits("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.6&&mu0_ProbNNmu>0.5","D_DiMuon_Mass>0&&D_DiMuon_Mass<2000");
  // myApplication.runAllFull1DFits("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.6&&mu1_ProbNNmu>0.5");
  //myApplication.runAllFull1DFits("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&mu0_MuonNShared==0&&mu1_MuonNShared==0","BDT>0.6&&mu1_ProbNNmu>0.5");
  //myApplication.compare_misID_shapes();
  //myApplication.compare_misID_shapes_2D();
  //myApplication.saveAllModelConfig("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.6&&mu1_ProbNNmu>0.5");
  //myApplication.compare_1D_and_2D_fit("BDT>0.7&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.7&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.7&&mu0_ProbNNmu>0.5");
  //myApplication.ExtractAllExpectedLimit("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5");
 
  /******/ 
  

  //StandardHypoTestInvDemo();  
  

 
  //////////////////////////////////////
  //                                  //
  //  Selection optimisation          //
  //  optimizes BDT and PID sim.      //
  //  defined in selectionOptimization//
  //                                  //
  //////////////////////////////////////

  /*
  //pipimumu bins
  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass<525");
  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass>565&&D_DiMuon_Mass<950");
  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100");
  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass>1100");
  newCutOptimization("D2pipimumu","Polarity==1","D_DiMuon_Mass<525");
  newCutOptimization("D2pipimumu","Polarity==1","D_DiMuon_Mass>565&&D_DiMuon_Mass<950");
  newCutOptimization("D2pipimumu","Polarity==1","D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
  newCutOptimization("D2pipimumu","Polarity==1","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100");
  newCutOptimization("D2pipimumu","Polarity==1","D_DiMuon_Mass>1100");

  //KKmumu bins
  newCutOptimization("D2KKmumu","Polarity==-1","D_DiMuon_Mass<525");
  newCutOptimization("D2KKmumu","Polarity==-1","D_DiMuon_Mass>565&&D_DiMuon_Mass<950");
  newCutOptimization("D2KKmumu","Polarity==-1","D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
  newCutOptimization("D2KKmumu","Polarity==1","D_DiMuon_Mass<525");
  newCutOptimization("D2KKmumu","Polarity==1","D_DiMuon_Mass>565&&D_DiMuon_Mass<950");
  newCutOptimization("D2KKmumu","Polarity==1","D_DiMuon_Mass>525&&D_DiMuon_Mass<565");

  //full di muon mass range
  newCutOptimization("D2KKmumu","Polarity==1","D_DiMuon_Mass>0");
  newCutOptimization("D2KKmumu","Polarity==-1","D_DiMuon_Mass>0");
  newCutOptimization("D2pipimumu","Polarity==1","D_DiMuon_Mass>0");
  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass>0");

  */


  /*
  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root");
  myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT.root");
  myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root");
  myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root");

  TString dataCut = "deltaM>144.5&&deltaM<146.5&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&BDT>0.4";
  TString misIDCut = "deltaM>144.5&&deltaM<146.5&&mu0_ProbNNmu>0.5&&BDT>0.4";
  TString q2Cut = "D_DiMuon_Mass<1100&&D_DiMuon_Mass>950";
  myFitter1D.fit_MC(dataCut+"&&"+q2Cut,true,"test.eps");
  myFitter1D.fit_normalization_MC(dataCut,true,"test2.eps");
  myFitter1D.fit_Kpipipi_misID(misIDCut,true,"test.eps");
  myFitter1D.fit_HHpipi_misID(misIDCut+"&&"+q2Cut,true,"test.eps"); ////q2 cut missing? check!                                            
  myFitter1D.fit_normalization_Data(dataCut,"test2.eps");
  myFitter1D.fit_resonant_Data(dataCut+"&&"+q2Cut,"test3.eps");
  */


  //define_binning("D2KKmumu");

  // Draw_MC_L0Muon_Efficiencies_2D();
  //optimizeSelection2D("nTracksOdd");
  //optimizeSelection("nTracksEven"); 
  //optimizeSelection2D("nTracksEven"); 
  //draw_BDT_crosschecks();
  //check_peakingBackground("D2KKmumu",true); //last argument decides if PID cuts are aplied or not
  //check_peakingBackground("D2KKmumu",false);

  //studyKKMCEfficiency(0.,0.5,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/EfficiencyStudies/offlineSelected.eps");
  //studypipiMCEfficiency(0.,0.5,"/work/mitzel/D2hhmumu/dev/D2pipimumu/img/EfficiencyStudies/offlineSelected.eps");
  
  //studyKKMCEfficiency(-10,-10,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/EfficiencyStudies/stripTrigSelected.eps");
  
  cout << "==============================================" << endl;
  cout << " Done " 
       << " \n Time since start " << (time(0) - startTime)/60.0
       << " min." << endl;
  cout << "==============================================" << endl;

   return 0;

}
