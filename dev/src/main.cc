#include <TChain.h>
#include <iostream>
#include <TString.h>
#include "LumiTupleReader.C"
#include "D2KKmumuReader.h"
#include "D2pipimumuReader.h"
#include "D2KpimumuReader.h"
#include "D2hhmumuFitter.h"
#include "sWeights.h"

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
  KK_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_PreselectedSubsample.root",5);
  //KK_Reader->fillHistograms("../rootFiles/test.root",isMC);                                                                          

}

void D2pipimumuData(){

  TChain* Tree_D2pipimumu = new TChain("DstD2PiPiMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2pipimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));
    Tree_D2pipimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }

  D2pipimumuReader* pipi_Reader = new D2pipimumuReader(Tree_D2pipimumu);
  pipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_PreselectedSubsample.root",2);
  //pipi_Reader->fillHistograms("../rootfiles/Data2012_Kinematical_Distrubutions_D2KK_noPreselection.root",false);                                                                           

}

void D2KpimumuData(){

  TChain* Tree_D2Kpimumu = new TChain("DstD2KPiMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2Kpimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));                                                                 
    Tree_D2Kpimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }

  D2KpimumuReader* Kpi_Reader = new D2KpimumuReader(Tree_D2Kpimumu);
  Kpi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_PreselectedSubsample.root",2);                                                                    
  //Kpi_Reader->createValidationSubsample("D2Kpimumu_ValidationSubsample.root");
}

void D2KKmumuMC(){


  TChain* Tree_MC_D2KKmumu = new TChain("MC12_DstD2KKMuMu/DecayTree");
  Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/magDw/MC12_DstD2KKmumu_magDw.root");                                                                            
  Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/magUp/MC12_DstD2KKmumu_magUp.root");                                                                            

  //Mai June TCK                                                                                                                                                                          
  Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_MaiJune/magDw/MC12_DstD2KKmumu_MaiJune_magDw.root");
  Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_MaiJune/magUp/MC12_DstD2KKmumu_MaiJune_magUp.root");

  D2KKmumuReader* KK_MC_Reader = new D2KKmumuReader(Tree_MC_D2KKmumu);
  KK_MC_Reader->InitMC();
  //KK_MC_Reader->fillHistograms("../rootfiles/MC2012_Kinematical_Distrubutions_D2KK.root",true);                                                                                          
  KK_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_MaiJune_MCtrainingSample.root");
  //KK_MC_Reader->studyTriggerEfficiency();                                                                                                                                                

}

void D2KpimumuMC(){


  TChain* Tree_MC_D2Kpimumu = new TChain("MC12_DstD2KPiMuMu/DecayTree");
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

  //D2pipimumuMC();                                                                                                                                                                          
  //D2pipimumuData();                                                                                                                                                                         
  
  //D2KpimumuMC();                                                                                                                                                                           
  //D2KpimumuData();

  //D2KKmumuMC();                                                                                                                                                                            
  //D2KKmumuData();                                                                                                                                                                             
  //calculateSweights("../Data_MC_Comparison/D2Kpimumu_BDT_selected.root");
  //data_MC_comparison();
  // checkLuminosity();      
  //fit_MC();
  fit_Data();
  return 0;

}
