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



void D2KKmumuData(){


  TChain* Tree_D2KKmumu = new TChain("DstD2KKMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2KKmumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));
    Tree_D2KKmumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }

  bool isMC= false;
  D2KKmumuReader* KK_Reader = new D2KKmumuReader(Tree_D2KKmumu);
  KK_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_PreselectedSubsample1.root",40);
  //KK_Reader->fillHistograms("../rootFiles/test.root",isMC);                                                                          

}

void D2pipimumuData(){

  TChain* Tree_D2pipimumu = new TChain("DstD2PiPiMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2pipimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));
    Tree_D2pipimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }

  D2pipimumuReader* pipi_Reader = new D2pipimumuReader(Tree_D2pipimumu);
  pipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_PreselectedSubsample.root",10);
  //pipi_Reader->fillHistograms("../rootfiles/Data2012_Kinematical_Distrubutions_D2KK_noPreselection.root",false);                                                                           

}

void D2KpimumuData(){

  TChain* Tree_D2Kpimumu = new TChain("DstD2KPiMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2Kpimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));                                                                 
    Tree_D2Kpimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }

  D2KpimumuReader* Kpi_Reader = new D2KpimumuReader(Tree_D2Kpimumu);
  Kpi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_PreselectedSubsample.root",10);                                                                    
  //Kpi_Reader->createValidationSubsample("D2Kpimumu_ValidationSubsample.root");
}

void D2KKmumuMC(){


  TChain* Tree_MC_D2KKmumu = new TChain("MC12_DstD2KKMuMu/DecayTree");
  Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/highStat_magDw/MC12_DstD2KKmumu_magDw.root");                                                         
  Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/highStat_magUp/MC12_DstD2KKmumu_magUp.root");                                                                            

  Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/magDw/MC12_DstD2KKmumu_magDw.root");
  Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/magUp/MC12_DstD2KKmumu_magUp.root");  

  //Mai June TCK                                                                                                                                                                          
  //Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_MaiJune/magDw/MC12_DstD2KKmumu_MaiJune_magDw.root");
  //Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_MaiJune/magUp/MC12_DstD2KKmumu_MaiJune_magUp.root");

  D2KKmumuReader* KK_MC_Reader = new D2KKmumuReader(Tree_MC_D2KKmumu);
  KK_MC_Reader->InitMC();
  //KK_MC_Reader->fillHistograms("../rootfiles/MC2012_Kinematical_Distrubutions_D2KK.root",true);                                                                                          
  KK_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_highStat_MCtrainingSample.root");
  //KK_MC_Reader->studyTriggerEfficiency();                                                                                                                                                

}


void D2KKpipiMC(){


  TChain* Tree_MC_D2KKpipi = new TChain("MC12_DstD2KKpipi/DecayTree");

  for (int i=0; i<50; ++i) {
    Tree_MC_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKpipi_filtered/magDw/MC12_DstD2KKpipi_%i.root",i));
    Tree_MC_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKpipi_filtered/magUp/D2KKpipi_%i.root",i));
  }
 

  D2KKpipiReader* KK_MC_Reader = new D2KKpipiReader(Tree_MC_D2KKpipi);
  KK_MC_Reader->InitMC();
                                                                                                                                                                            
  KK_MC_Reader->addMisIdMasses("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_filteredt_MCSample1.root");
 

}


void D2KpipipiMC(){


  TChain* Tree_MC_D2Kpipipi = new TChain("MC12_DstD2Kpipipi/DecayTree");

  for (int i=0; i<50; ++i) {
    Tree_MC_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpipipi_filtered/magDw/MC12_DstD2Kpipipi_%i.root",i));
    Tree_MC_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpipipi_filtered/magUp/D2Kpipipi_%i.root",i));
  }


  D2KpipipiReader* Kpi_MC_Reader = new D2KpipipiReader(Tree_MC_D2Kpipipi);
  Kpi_MC_Reader->InitMC();
                                                                                                                                                                         
  Kpi_MC_Reader->addMisIdMasses("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_filteredt_MCSample1.root");


}



void D2KpimumuMC(){


  TChain* Tree_MC_D2Kpimumu = new TChain("MC12_DstD2KKpipi/DecayTree");
  Tree_MC_D2Kpimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu/magDw/MC12_DstD2Kpimumu_magDw.root");
  Tree_MC_D2Kpimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu/magUp/MC12_DstD2Kpimumu_magUp.root");

  D2KpimumuReader* Kpi_MC_Reader = new D2KpimumuReader(Tree_MC_D2Kpimumu);
  Kpi_MC_Reader->InitMC();
  //KK_MC_Reader->fillHistograms("../rootfiles/MC2012_Kinematical_Distrubutions_D2KK.root",true);                                                                                          
  Kpi_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root");
  //KK_MC_Reader->studyTriggerEfficiency();                                                                                                                                                

}

void D2pipimumuMC(){

  TChain* Tree_MC_D2pipimumu = new TChain("MC12_DstD2pipiMuMu/DecayTree");
  Tree_MC_D2pipimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu/magDw/MC12_DstD2pipimumu_magDw.root");
  Tree_MC_D2pipimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu/magUp/MC12_DstD2pipimumu_magUp.root");

  D2pipimumuReader* pipi_MC_Reader = new D2pipimumuReader(Tree_MC_D2pipimumu);
  pipi_MC_Reader->InitMC();
  //pipi_MC_Reader->fillHistograms("../rootfiles/MC2012_Kinematical_Distrubutions_D2KK.root",true);                                                                                        
  pipi_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_MCtrainingSample.root");
  //pipi_MC_Reader->studyTriggerEfficiency();                                                                                                                                              

}



int main() {



  ////////////////////////////////////
  //                                //
  //  create preselected samples    //
  //  split by run numbers          //
  //  create TMVA trainign samples  //
  //                                //
  ////////////////////////////////////



  //D2pipimumuMC();                                                                                                                                                                          
  //D2pipimumuData();                                                                                                                                                                         
  //D2KpimumuMC();                                                                                                                                                                           
  //D2KpimumuData();

  //D2KKmumuMC();                                                                                                                                                         
  //D2KKmumuData();                                                                                                                                                   
  //D2KKpipiMC();                                                                                                                                                        
  //D2KpipipiMC();                                                                                                                                                          


  ////////////////////////////////////                                                                                                                                                         //                                //                                                                                                                                                         // data MC comparison with        //
  //   sweights                     //                                                                                                                                
  //                                //
  ////////////////////////////////////

  //sweights part                     
  //calculateSweights("../Data_MC_Comparison/D2Kpimumu_BDT_selected.root");
  //data_MC_comparison();
  // checkLuminosity();      


  ////////////////////////////////////
  //                                //
  // Fitter                         //
  //                                //
  ////////////////////////////////////

  
  //Fit part
  D2hhmumuFitter myFitter;
  //myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root");
  myFitter.fit_MC(true);
  myFitter.fit_PIDinverted_Data(true);
  myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter.fit_Data("BDT>0.7 && mu0_ProbNNmu>0.5 && mu1_ProbNNmu>0.5"); 
 //myFitter.getCombBkg("BDT>0.6 && mu0_ProbNNmu>0.4 && mu1_ProbNNmu>0.4","");
 //myFitter.toyStudy();


  ////////////////////////////////////
  //                                //
  //  TMVA training and application //
  //                                //
  ////////////////////////////////////


  //TMVA
  //D2KKmumuCrosstraining();
  //D2KKmumuCrossapplication();


  ////////////////////////////////////
  //                                //
  //  Selection optimisation        //
  //                                //
  ////////////////////////////////////

  // optimizeSelection();

   return 0;

}
