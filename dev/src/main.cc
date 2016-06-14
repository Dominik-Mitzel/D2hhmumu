#include <TChain.h>
#include <iostream>
#include <TString.h>
#include "LumiTupleReader.C"
#include "D2KKmumuReader.h"
#include "D2pipimumuReader.h"
#include "D2KpimumuReader.h"
#include "D2KKpipiReader.h"
#include "D2KpipipiReader.h"
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


  // D2pipimumuMC();      
  //D2pipimumuData();                                                                                                                
  //D2KpimumuMC();        
  //D2KpimumuData();
  //D2KKmumuMC();                                                                                                                     
  //D2KKmumuData();
  //D2KKpipiMC();                                                                                   
  //D2KpipipiMC();                         
  //D2KKpipiData();                                                                                                               
  //D2KpipipiData();                                                                                        
  //createGeneratorLevelMCTuple("D2KKmumu");
  //createGeneratorLevelMCTuple("D2Kpimumu");
  //createGeneratorLevelMCTuple("D2KKpipi");
  //createGeneratorLevelMCTuple("D2Kpipipi");
  
  EfficiencyCalculator myEfficiencies("D2KKmumu");
  myEfficiencies.getMCSignalMisIDEfficiency("BDT>0.5","nTracks%2!=0","D_DiMuon_Mass>565");
  myEfficiencies.getMCRelativeSigToNormMisIDEfficiency("BDT>0.5","nTracks%2!=0","D_DiMuon_Mass>565","D_DiMuon_Mass>675&&D_DiMuon_Mass<875");

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

  
  //D2hhmumuFitter_Applications myApplication("D2KKmumu","2012");
  //myApplication.runFull1DFits("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5 ","BDT>0.6&&mu0_ProbNNmu>0.5","D_DiMuon_Mass>0&&D_DiMuon_Mass<2000");
  //myApplication.runAllFull1DFits("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5 ","BDT>0.6&&mu0_ProbNNmu>0.5");
  //myApplication.performAllToyStudies();
  //myApplication.compare_misID_shapes();
  //myApplication.compare_misID_shapes_2D();
  //myApplication.saveAllModelConfig("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5 ","BDT>0.6&&mu0_ProbNNmu>0.5");
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

  
         
  //Fit part 1D
  /*
  D2hhmumuFitter1D myFitter;
  myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  myFitter.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  //myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root"); 
  myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpipipi_D2KKmumuBDT.root");
  myFitter.setPathToKKpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root");

  std::cout << "misID "<< myFitter.getMisIDbkgExp("BDT>0&&mu0_ProbNNmu>0.3&&mu1_ProbNNmu>0.3&&nTracks%2==0&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875","test.eps"); 
  EfficiencyCalculator myEfficiency("D2KKmumu");
  std::cout<<"getMCRelativeSigToNormMisIDEfficiency "<< myEfficiency.getMCRelativeSigToNormMisIDEfficiency("BDT>0&&mu1_ProbNNmu>0.3","nTracks%2==0","D_DiMuon_Mass>525","D_DiMuon_Mass>675&&D_DiMuon_Mass<875") <<std::endl;
  std::cout<<"getMisIDFractionQ2Range "<< myEfficiency.getMisIDFractionQ2Range("D_DiMuon_Mass>525") <<std::endl;
  */
  // myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  //myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  //myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  //myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");

  //myFitter1D.makeToyStudy("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","BDT>0.5&&mu1_ProbNNmu>0.5","test.eps",50,50,5);

  // myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  //myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  //myFitter1D.setPathToKpipipiHistoData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/zoz-5000.root");
  //myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  //myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");   


  //myFitter1D.fillModelConfig("BDT>0.7 && mu0_ProbNNmu>0.5 && mu1_ProbNNmu>0.5 ","Data_model1D.root");
  //myFitter1D.GausExpModel(100,100);
  //myFitter1D.ExtractLimit();
  
  //myFitter1D.getCombBkg("BDT>0.6&&mu0_ProbNNmu>0.4&&mu1_ProbNNmu>0.4","");

  //myFitter1D.fit_MC("BDT>0.6&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",true,"test1.eps");
  //myFitter1D.fit_Kpipipi_misID("BDT>0.6&&mu1_ProbNNmu>0.5",true,"test2.eps");
  //myFitter1D.fit_Kpipipi_misID_fromHistogramm("BDT>0.5&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5",true,""); 
  //myFitter1D.fit_normalization_Data("BDT>0.6&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","test3.eps");
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
  //D2KKmumuCrossapplication();


  //////////////////////////////////////
  //                                  //
  //  Selection optimisation          //
  //  optimizes BDT and PID sim.      //
  //  defined in selectionOptimization//
  //                                  //
  //////////////////////////////////////

  //optimizeSelectionInBins("nTracksOdd","D_DiMuon_Mass<525");
  //optimizeSelectionInBins("nTracksOdd","D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
  //optimizeSelectionInBins("nTracksOdd","D_DiMuon_Mass>565");

  //optimizeSelectionInBins("nTracksEven","D_DiMuon_Mass<525");
  //optimizeSelectionInBins("nTracksEven","D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
  //optimizeSelectionInBins("nTracksEven","D_DiMuon_Mass>565");
  
  //optimizeSelection2D("nTracksOdd");
  //optimizeSelection("nTracksEven"); 
  //optimizeSelection2D("nTracksEven"); 
  //draw_BDT_crosschecks();
  //check_peakingBackground("D2KKmumu",true); //last argument decides if PID cuts are aplied or not
  //check_peakingBackground("D2KKmumu",false);


  cout << "==============================================" << endl;
  cout << " Done " 
       << " \n Time since start " << (time(0) - startTime)/60.0
       << " min." << endl;
  cout << "==============================================" << endl;

   return 0;

}
