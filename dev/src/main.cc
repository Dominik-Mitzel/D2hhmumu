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


  //ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(1);
  //ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-12);

  // checkLuminosity();       

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
  //D2KpimumuData();
  //D2KpimumuMC();        
  //D2KKmumuMC();                                                                                                                     
  //D2KKmumuData();
  
  //misID
  //D2pipipipiData();
  //D2KKpipiData();                                                                                                               
  //D2KpipipiData();                                                                                                               
  //D2pipipipiMC();
  //D2KKpipiMC();                                                                                   
  //D2KpipipiMC();                         
  //randomize the order of the pions for the double misID K3Pi and 4pi
  //these are the one that should be used in the final analysis for K3pi and 4pi
  //D2pipipipiRandomizedData();
  //D2KpipipiRandomizedData();

  //create the generator level MC tuples for all channels.
  //createGeneratorLevelMCTuple("D2KKmumu");
  //createGeneratorLevelMCTuple("D2Kpimumu");
  //createGeneratorLevelMCTuple("D2KKpipi");
  //createGeneratorLevelMCTuple("D2Kpipipi");
  //createGeneratorLevelMCTuple("D2pipimumu");
  //createGeneratorLevelMCTuple("D2pipipipi");  
  
  //create the tuples binned in defined D_DiMuon mass binning scheme, also needed to efficiency evaluation
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
  
  //Application_D2KKmumu("DecayTree","flaggedMC/MC12_DstD2KKmumu_RerunCalo.root","flaggedMC/MC12_DstD2KKmumu_RerunCalo_BDT.root",1,true,false,false);
  //Application_D2pipimumu("DecayTree","flaggedMC/MC12_DstD2pipimumu_RerunCalo.root","flaggedMC/MC12_DstD2pipimumu_RerunCalo_BDT.root",1,true,false,false);

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

  //here the different method are listed. the first argument is set true if a fit is to be used to get the signal yield. second argument is cut (on top of the muon Nshared and delta M cut )
  //results can be found in ../img/EfficiencyStudies/strippingEfficiency/RecoAndStrippingEfficiency

  //drawRecoAndStrippingEfficiencyWithFit(true,"(Dst_BKGCAT<11||Dst_BKGCAT==60)");
  //drawRecoAndStrippingEfficiencyWithFit(false,"Dst_BKGCAT<11"); 
  //drawRecoAndStrippingEfficiencyWithFit_alternativeFitModel(true,"(Dst_BKGCAT<11||Dst_BKGCAT==60)");                                                                                   
  //drawRecoAndStrippingEfficiencyWithFit(true,"(Dst_BKGCAT<11||Dst_BKGCAT==60||Dst_BKGCAT==50)");
  //drawRecoAndStrippingEfficiencyWithFit_alternativeFitModel(true,"(Dst_BKGCAT<11||Dst_BKGCAT==60||Dst_BKGCAT==50)");                       
  
  //drawRecoAndStrippingEfficiencyWithFit(true,"(Dst_BKGCAT<11||Dst_BKGCAT==60||Dst_BKGCAT==50)");

  /*
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","PHSP","");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","phi","");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","rho","");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","omega","");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","all","");

  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","","PHSP");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","","phi");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","","rho");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","","omega");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","","all");
  */
  /*
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","omega","phi");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","phi","phi");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","rho","phi");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","PHSP","phi");

  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","omega","rho");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","phi","rho");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","rho","rho");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","PHSP","rho");

  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","omega","PHSP");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","phi","PHSP");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","rho","PHSP");
  drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","PHSP","PHSP");
  */
  
  // drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","all","PHSP");
  //drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","all","phi");
  //drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","all","all");
  //drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","all","rho");

  //drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","rho","all");
  //drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","omega","all");
  //drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","phi","all");
  //drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","PHSP","all");
  //drawRecoAndStrippingEfficiencyWithFit_PhaseSpaceModel(false,"Dst_BKGCAT<11","all","all");

  //compareDifferentModelsForRecoEfficiency();



  //compareDifferentRecoAndStrippingEfficiencies();

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

  //calibratePIDEfficiencyWithDsPhiPi();
  //evaluateDsPhiPiMuonEfficiencyForLowPtBins("magDw");
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
  
  //this is the one which is supposed to be used. Takes the filtere MC samples and optionally allows for cat60 or not
  //MC_Combined_Hlt_BDT_efficiency(0.5,0.2,0.5,0.4,true,false);
  //MC_Combined_Hlt_BDT_efficiency(0.5,0.2,0.5,0.4,true,true);

  
  //check for model dependency of MC efficiency
  //MC_Combined_efficiency_onlyPhaseSpaceModel(.5,.2,.5,.4,true,false);
  //MC_Combined_efficiency_onlyPhaseSpaceModel(.5,.2,.5,.4,false,true);
  //MC_Combined_efficiency_onlyPhaseSpaceModel(.5,.2,.5,.4,false,false);
  
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
  
  //applyTagAndProbeL0EfficiencyToyStudy(true,500); 
  //applyTagAndProbeL0EfficiencyToyStudy(false,500); 
  


  //createDataTuplesForEffStudies();
  // Evaluate TOTAL relative efficiency
  //drawTotalEfficiency();
  drawTotalEfficiencyForCorrelationStudies();
  //drawBothPolaritiesPIDCalibHadronEfficiency("default");

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

  //add a branch with a flag to get rid of multiple candidates
  //createTuplesWithSelectedMultipleCand();
  
  //writeAllTrueCanidatesToFile(); //write all candidates of the true MC tuple to be able to perform a matching later with the reconstructed candidates. This thing is not uesd at the moment
  
  //write_ghosts_withMultipleCand_toFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_200.0_1600.0_D2pipimumuBDT.root","D2pipimumu");
  //write_ghosts_withMultipleCand_toFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_200.0_1600.0_D2KKmumuBDT.root","D2KKmumu");
 

  ////////////////////////////////////////////////
  //
  //multiple Candidate Data
  //
  //                                                                                                                                                             
  ////////////////////////////////////////////////  

  //writeMultipleDataCanidatesToFile();
  //chose_multiple_Data_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root","D2pipimumu");
  //chose_multiple_Data_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT.root","D2KKmumu");
  //chose_multiple_Data_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root","D2Kpimumu_D2KKmumuBDT");
  //chose_multiple_Data_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root","D2Kpimumu_D2pipimumuBDT"); 
  //createDataTuplesWithSelectedMultipleCand();

  /////////////////////////////////////////////////

  //check nature of multiple candidate on a track by track basis
  //writeFullKinematicMultipleDataCanidatesToFile();
  //chose_FullKinematic_multiple_Data_events("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root","D2pipimumu");


  ////////////////////////////////////
  //
  //norm mode and D->hhpipi sWeights to check MC data differences
  // 
  ///////////////////////////////////

  /* 
  D2hhmumuFitter1D myFitter1D_1;
  myFitter1D_1.addNormalizationSWeightsHadronicChannel("eventNumber%5==0","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/temp/D2KKpipi_D2KKmumuBDT.root","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/temp/D2KKpipi_D2KKmumuBDT_sWeights.root","drawingMacros/sWeights_D2KKpipi.eps");

  D2hhmumuFitter1D myFitter1D_2;
  myFitter1D_2.addNormalizationSWeightsHadronicChannel("eventNumber%500==0","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/temp/D2Kpipipi_D2KKmumuBDT.root","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/temp/D2Kpipipi_D2KKmumuBDT_sWeights.root","drawingMacros/sWeights_D2KKpipi.eps");
  */

  /*
  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2pipimumuBDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root");
  myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");                                                          
  myFitter1D.addNormalizationSWeights("mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&Dst_DTF_D0_M<1940&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875","mu0_ProbNNmu>0.5","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root","/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT_sWeights.root");
  */

  //////////////////////////////////////////////////


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
  myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_Randomized.root");  
  //myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT.root");  
  myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_Randomized.root");  
  myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT_noMultCand.root");  
  myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT_noMultCand.root");
  //myFitter1D->fillModelConfig("D2pipimumu","BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5",0,"test.root");

  myFitter1D->fit_MC("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>525&&D_DiMuon_Mass<565",true,"test.eps","m(KK#mu#mu)","950MeV<m(#mu#mu)<1100MeV");
  myFitter1D->fit_normalization_MC("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5",true,"test1.eps");
  myFitter1D->fit_Kpipipi_misID("BDT>0.4&&mu1_ProbNNmu>0.5",true,"test2.eps",true);
  myFitter1D->fit_HHpipi_misID("BDT>0.4&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>525&&D_DiMuon_Mass<565",true,"test3.eps","m_{#pi#pi#mu#mu}(#pi#pi#pi#pi)","950MeV<m(#pi#pi)<1100MeV",true);
  myFitter1D->fit_normalization_Data("BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","test4.eps");
  //myFitter1D->fit_Data("D2pipimumu","BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","D_DiMuon_Mass>525&&D_DiMuon_Mass<565","test5","m(#pi^{+}#pi^{-}#mu^{+}#mu^{-})","950<m(#mu^{+}#mu^{-})<1100 MeV/c^{2}",true);  
  //myFitter1D->get_Significance("D2pipimumu","BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","D_DiMuon_Mass>525&&D_DiMuon_Mass<565","test5","m(#pi^{+}#pi^{-}#mu^{+}#mu^{-})","950<m(#mu^{+}#m^{-})<1100 MeV/c^{2}",true);  
  myFitter1D->fit_unblinded_Data("D2pipimumu","BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","D_DiMuon_Mass>525&&D_DiMuon_Mass<565","test5","m(K^{+}K^{-}#mu^{+}#mu^{-})","525<m(#mu^{+}#mu^{-})<565 MeV/c^{2}",false);  

  //myFitter1D->fillModelConfig("D2pipimumu","BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5",0,"testConfig.root");
  */
  /*            
  D2hhmumuFitter_Applications myApplication("D2KKmumu","2012");
  //myApplication.drawMisIDShapes("BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication.drawMCSignalShapes("BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5");
  myApplication.runFullResonant1DFits("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication.runFull1DFits("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication.addAllSignalSweights("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication.saveModelConfig("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5"); 
  //myApplication.ExtractExpectedLimit();
  //myApplication.performAllToyStudies();
 
  D2hhmumuFitter_Applications myApplication2("D2pipimumu","2012");
  //myApplication2.draw_misID_shapes_singlePlot();
  //myApplication2.drawMisIDShapes("BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication2.drawMCSignalShapes("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5");
  //myApplication2.performAllToyStudies();
  //myApplication2.runFull1DFits("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication2.addAllSignalSweights("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5");
  myApplication2.runFullResonant1DFits("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication2.saveModelConfig("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.4&&mu1_ProbNNmu>0.5");
  //myApplication2.ExtractExpectedLimit();
  //myApplication2.studyCombBkgShape();
  //myApplication2.compare_misID_shapes("BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5");
  //myApplication2.studyResolutionScale("");
  */
 
  /*
  /******/ 

 
  //////////////////////////////////////
  //                                  //
  //  Selection optimisation          //
  //  optimizes BDT and PID sim.      //
  //  defined in selectionOptimization//
  //                                  //
  //////////////////////////////////////

  /*  

  //pipimumu bins

  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass>565&&D_DiMuon_Mass<950");
  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100");
  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass>1100");
  newCutOptimization("D2pipimumu","Polarity==-1","D_DiMuon_Mass<525");

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
  
//draw_BDT_crosschecks();

 

  
  cout << "==============================================" << endl;
  cout << " Done " 
       << " \n Time since start " << (time(0) - startTime)/60.0
       << " min." << endl;
  cout << "==============================================" << endl;

   return 0;

}
