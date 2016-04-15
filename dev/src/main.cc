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
  //D2pipimumuData();                                                                                                                                       //D2KpimumuMC();                                                                                                                                          //D2KpimumuData();
  //D2KKmumuMC();                                                                                                                                           //D2KKmumuData();                                                                                                                                         //D2KKpipiMC();                                                                                                                                           //D2KpipipiMC();                                                                                                                                          //D2KKpipiData();                                                                                                                                       
  //D2KpipipiData();                                                                                                                                                          


  ////////////////////////////////////                                                                                                                      //                                //                                                                                                                      // data MC comparison with        //
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

  //Fit part
  D2hhmumuFitter myFitter;
  //myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root");
  //myFitter.fit_MC("",true);
  //myFitter.fit_PIDinverted_Data(true);
  //myFitter.getCombBkg("BDT>0.6 && mu0_ProbNNmu>0.4 && mu1_ProbNNmu>0.4","");
  //myFitter.fit_Kpipipi_misID("BDT>0.7 && mu1_ProbNNmu>0.5",true);
  //myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  myFitter.fit_normalization_Data("BDT>0.7 && mu0_ProbNNmu>0.5 && mu1_ProbNNmu>0.5 ");
  // myFitter.fit_Data("BDT>0.7 && mu0_ProbNNmu>0.5 && mu1_ProbNNmu>0.5 "); 
 
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

  // optimizeSelection();

  cout << "==============================================" << endl;
  cout << " Done " 
       << " \n Time since start " << (time(0) - startTime)/60.0
       << " min." << endl;
  cout << "==============================================" << endl;

   return 0;

}
