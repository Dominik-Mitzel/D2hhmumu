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


  //D2pipimumuMC();      
  //D2pipimumuData();                                                                                                                
  //D2KpimumuMC();        
  //D2KpimumuData();
  //D2KKmumuMC();                                                                                                                     
  //D2KKmumuData();
  //D2pipipipiData();

  //D2KKpipiMC();                                                                                   
  //D2KpipipiMC();                         

  //D2KKpipiData();                                                                                                               
  //D2KpipipiData();                                                                                                               
  
  //createGeneratorLevelMCTuple("D2KKmumu");
  //createGeneratorLevelMCTuple("D2Kpimumu");
  //createGeneratorLevelMCTuple("D2KKpipi");
  //createGeneratorLevelMCTuple("D2Kpipipi");
  
  //EfficiencyCalculator myEfficiencies("D2KKmumu");
  //myEfficiencies.getMCSignalMisIDEfficiency("BDT>0.5","nTracks%2!=0","D_DiMuon_Mass>565");
  //myEfficiencies.getMCRelativeSigToNormMisIDEfficiency("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","","D_DiMuon_Mass>565","D_DiMuon_Mass>675&&D_DiMuon_Mass<875");

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


  //D2hhmumuFitter_Applications myApplication("D2pipimumu","2012");
  //myApplication.runAllResonantFull1DFits("BDT>0&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","BDT>0.9&&mu1_ProbNNmu>0.5");
  
  D2hhmumuFitter_Applications myApplication("D2KKmumu","2012");
  //myApplication.runAllFull1DFits("BDT>0&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","BDT>0.95&&mu1_ProbNNmu>0.5");
  myApplication.performAllToyStudies();
  
  //D2hhmumuFitter_Applications myApplication("D2KKmumu","2012");
  //myApplication.runAllResolutionScaleStudies("BDT>0.6","BDT>0.6");
  //myApplication.studyResolutionScale("BDT>0.6","BDT>0.6","D_DiMuon_Mass>0","D2KKpipi");
  //myApplication.studyResolutionScale("BDT>0.6","BDT>0.6","D_DiMuon_Mass>0","D2Kpipipi");  
  //myApplication.studyResolutionScale("BDT>0.6","BDT>0.6","D_DiMuon_Mass>525&&D_DiMuon_Mass<565","D2KKpipi");
  //myApplication.studyResolutionScale("BDT>0.6","BDT>0.6","D_DiMuon_Mass>565","D2Kpipipi");
  
  //myApplication.studyNormalizationFits("BDT>0.6&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","BDT>0.6&&mu1_ProbNNmu>0.5",true); 
  //myApplication.studyNormalizationFits("BDT>0.6&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","BDT>0.6&&mu1_ProbNNmu>0.5",false); 
  /*
   D2hhmumuFitter1D* myFitter1D = new D2hhmumuFitter1D();
   myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
   myFitter1D->setPathToKKpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_D2KKmumuBDT_even1.root");  
   myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_D2KKmumuBDT_even1.root");  
   myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");  
  
   //myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmummu_BDT.root");
  myFitter1D->fit_MC("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5",true,"test.eps");
  myFitter1D->fit_Kpipipi_misID("BDT>0.6&&mu1_ProbNNmu>0.5",true,"test2.eps");

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

  /*
  //Fit part 2D
  D2hhmumuFitter myFitter;
  myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter.fit_MC("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",true,"test1_2D.eps");
  myFitter.fit_Kpipipi_misID("BDT>0.5&&mu1_ProbNNmu>0.5",true,"test2_2D.eps");
  //myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  //myFitter.getCombBkg("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","temp");
  //myFitter.fit_PIDinverted_Data(true);
  //myFitter.getCombBkg("BDT>0.6 && mu0_ProbNNmu>0.4 && mu1_ProbNNmu>0.4","");
  //myFitter.fit_Kpipipi_misID("BDT>0.7 && mu1_ProbNNmu>0.5",true);
  //myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
    myFitter.fit_normalization_Data("BDT>0.5&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","test3_2D.eps");
  //myFitter.fit_Data("BDT>0.7 && mu0_ProbNNmu>0.5 && mu1_ProbNNmu>0.5 "); 
  */

  /****** 
         
  //Fit part 1D
  D2hhmumuFitter1D myFitter1D;

  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root");
  myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root");
  //myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root");
  myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");

  myFitter1D.fit_MC("BDT>0.&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5",true,"test.eps");
  myFitter1D.fit_normalization_MC("BDT>0&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5",true,"test2.eps");
  myFitter1D.fit_Kpipipi_misID("BDT>0.95&&mu1_ProbNNmu>0.5",true,"test3.eps");
  myFitter1D.fit_HHpipi_misID("BDT>0.95&&mu1_ProbNNmu>0.5",true,"test4.eps");
  myFitter1D.fit_normalization_Data("BDT>0&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","test5.eps");
  //myFitter1D.fit_Data("BDT>0&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>565 ","test6.eps");
  //myFitter1D.fit_resonant_Data("BDT>0.&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass<565&&D_DiMuon_Mass>525","test6.eps");

  */

  //myFitter1D.fit_Kpipipi_misID_fromHistogramm("BDT>0.5&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5",true,""); 
  //myFitter1D.fit_normalization_Data("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&h0_PIDK>0&&h1_PIDK>0&&mu1_MuonNShared==0&&mu0_MuonNShared==0","test3.eps");
  


  //myFitter1D.makeToyStudy("BDT>0.6&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass>525&&D_DiMuon_Mass<565","BDT>0.6&&mu1_ProbNNmu>0.5","test.eps",10,10,10,false);

  // myFitter1D.addNormalizationSWeights("mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","mu1_ProbNNmu>0.5");
  
  //std::cout << "misID "<< myFitter.getMisIDbkgExp("BDT>0&&mu0_ProbNNmu>0.3&&mu1_ProbNNmu>0.3&&nTracks%2==0&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875","test.eps"); 
  //EfficiencyCalculator myEfficiency("D2KKmumu");
  //std::cout<<"getMCRelativeSigToNormMisIDEfficiency "<< myEfficiency.getMCRelativeSigToNormMisIDEfficiency("BDT>0&&mu1_ProbNNmu>0.3","nTracks%2==0","D_DiMuon_Mass>525","D_DiMuon_Mass>675&&D_DiMuon_Mass<875") <<std::endl;
  //std::cout<<"getMisIDFractionQ2Range "<< myEfficiency.getMisIDFractionQ2Range("D_DiMuon_Mass>525") <<std::endl;
  
  //myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  //myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  //myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  //myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");

  //myFitter1D.makeToyStudy("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","BDT>0.5&&mu1_ProbNNmu>0.5","test.eps",50,50,5);
  //myFitter1D.fit_Data("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","test4.eps"); 
  //myFitter1D.getCombBkgFromDeltaM("BDT>0.5&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","");
  

  //StandardHypoTestInvDemo();  
  

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

  //////////////////////////////////////
  //                                  //
  //  Selection optimisation          //
  //  optimizes BDT and PID sim.      //
  //  defined in selectionOptimization//
  //                                  //
  //////////////////////////////////////

  // optimizeSelectionInBins("nTracksOdd","D_DiMuon_Mass<525");
  //optimizeSelectionInBins("nTracksOdd","D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
  //optimizeSelectionInBins("nTracksOdd","D_DiMuon_Mass>565");
  
  //optimizeSelectionInBins("nTracksEven","D_DiMuon_Mass<525");
  //optimizeSelectionInBins("nTracksEven","D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
  //optimizeSelectionInBins("nTracksEven","D_DiMuon_Mass>565");

  //optimizeSelectionInBins("nTracksEven","D_DiMuon_Mass>0");
  //optimizeSelectionInBins("nTracksOdd","D_DiMuon_Mass>0");

  
  //optimizeSelection2D("nTracksOdd");
  //optimizeSelection("nTracksEven"); 
  //optimizeSelection2D("nTracksEven"); 
  //draw_BDT_crosschecks();
  //check_peakingBackground("D2KKmumu",true); //last argument decides if PID cuts are aplied or not
  //check_peakingBackground("D2KKmumu",false);

  //studyKKMCEfficiency(0.6,0.5,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/EfficiencyStudies/offlineSelected.eps");
  //studyKKMCEfficiency(-10,-10,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/EfficiencyStudies/stripTrigSelected.eps");
  
  cout << "==============================================" << endl;
  cout << " Done " 
       << " \n Time since start " << (time(0) - startTime)/60.0
       << " min." << endl;
  cout << "==============================================" << endl;

   return 0;

}
