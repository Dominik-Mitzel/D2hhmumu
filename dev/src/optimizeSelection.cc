#include "optimizeSelection.h"
#include "D2hhmumuFitter.h"
#include "TEventList.h"
#include "TPaletteAxis.h"
#include "D2KKmumuReader.h"
#include "D2pipimumuReader.h"
#include "D2KpimumuReader.h"
#include "D2KKpipiReader.h"
#include "D2KpipipiReader.h"



void optimizeSelection() {

  TH2* h2= new TH2D("h2_cutEfficiency","h3_cutEfficiency",15,-0.5,1.0,10,0,1);
  TH2* h2_misIDBkg = new TH2D("h2_misIDBkg","h2_misIDBkg",15,-0.5,1,10,0,1);
  TH2* h2_combBkg = new TH2D("h2_combBkg","h2_combBkg",15,-0.5,1,10,0,1);
  TH2* h2_sigEff = new TH2D("h2_sigEff","h2_sigEff",15,-0.5,1,10,0,1);


  //double PIDcuts[20]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95};
  //double BDTcuts[20]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95};
  //double PIDcuts[14]={0,0.1,0.2,0.3,0.4,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9};
  //double BDTcuts[14]={0,0.1,0.2,0.3,0.4,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9};
  
  double PIDcuts[10]={0,0.1,0.2,0.3,0.4,0.5,0.60,0.7,0.8,0.9};                                                 
  //  double BDTcuts[10]={0,0.1,0.2,0.3,0.4,0.5,0.60,0.7,0.8,0.9};
  double BDTcuts[15]={-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  
  TString cut;
  TString namePlot;

  double misIDBkg;
  double combBkg;
  double misIDBKG_scaled;
  double EffRatio=1;
  double signalEff=1;

  D2hhmumuFitter myFitter;
  myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root");
  myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  myFitter.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter.setPathToInvData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/inverted_D2Kpimumu_BDT_selected.root");

  myFitter.fit_MC("",true); //fix MC signal shape
  myFitter.fit_PIDinverted_Data(true);
  double FOM;

  for(int i=0;i<15;++i){
    for(int j=0;j<10;++j){

      cut = TString::Format("BDT>%f&&mu0_ProbNNmu>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j],PIDcuts[j]);
      namePlot =  TString::Format("BDT%fPID%f",BDTcuts[i],PIDcuts[j]);

      //std::cout<<cut<<std::endl;
      misIDBkg = myFitter.getMisIDbkgExp(cut,namePlot);
      combBkg = myFitter.getCombBkg(cut,namePlot);
      signalEff = getMCSignalEfficiency(cut);  
      EffRatio=EffD2KKpipiToEffD2Kpipipi(cut);
      misIDBKG_scaled= misIDBkg/(33.14*2)*EffRatio;
      
      FOM=signalEff/( 5/2 + TMath::Sqrt(misIDBKG_scaled + combBkg) );

      h2->SetBinContent(i+1,j+1,FOM);
      h2_misIDBkg->SetBinContent(i+1,j+1,misIDBKG_scaled); 
      h2_combBkg ->SetBinContent(i+1,j+1,combBkg);
      h2_sigEff ->SetBinContent(i+1,j+1,signalEff);

      std::cout<<i<<" "<<j<<"  "<<FOM<<"  "<<signalEff << "  " <<misIDBKG_scaled <<"  " << combBkg <<std::endl;

    }
  }

 TCanvas* c = new TCanvas("canvas","canvas"); 
			  
 c->Divide(2,2);
			  
  Float_t newMargin1 = 0.13;
  Float_t newMargin2 = 0.15;

  c->SetGrid();
  c->SetTicks();
  
  gStyle->SetPalette( 1, 0 );

  gStyle->SetPaintTextFormat( "2g" );
  gStyle->SetOptStat(0);

  //h2->SetMarkerSize( 1.5 );
  h2->SetMarkerColor( 0 );
  //Float_t labelSize = 0.040;
  //h2->LabelsOption( "d" );
  c->cd(1);
  //h2->SetLabelOffset( 0.011 );// label offset on x axis                                                                                                                                 
  h2->Draw("colz"); // color pads                                                                                                                                               
  h2->Draw("textsame");  // add text      
  c->Update();
        
  c->cd(2);
  h2_misIDBkg->Draw("colz");
  h2_misIDBkg->Draw("textsame");
  c->cd(3);
  h2_combBkg->Draw("colz");
  h2_combBkg->Draw("textsame");
  c->cd(4);
  h2_sigEff->Draw("colz");
  h2_sigEff->Draw("textsame");
  c->Update();
  c->Print("../img/selection.eps");


}

double getMCSignalEfficiency(TString cut){

  TFile* file;
  TString pathToSignalMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root";
  file= new TFile(pathToSignalMC,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  double nTotal = tree->GetEntries();
  
  double nSelected = tree->GetEntries(cut);
  
  file->Close();
  
  return nSelected/nTotal;
}


double EffD2KKpipiToEffD2Kpipipi(TString cut){

  TFile* file1;
  TFile* file2;
  TString pathToKKpipiMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_filtered_D2KKmumuBDT.root";
  TString pathToKpipipiMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_filtered_D2KKmumuBDT.root";
 
  file1= new TFile(pathToKKpipiMC,"OPEN");
  file2= new TFile(pathToKpipipiMC,"OPEN");

  TTree* tree1 = (TTree*) file1->Get("BDT_Tree");
  TTree* tree2= (TTree*) file2->Get("BDT_Tree");

  double nSelected1 = tree1->GetEntries(cut);
  double nSelected2 = tree2->GetEntries(cut);


  file1->Close();  
  file2->Close();
  
  return nSelected1/nSelected2; //nKKpipi/nKpipipi

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
  Kpi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_PreselectedSubsample.root",10);                                  //Kpi_Reader->createValidationSubsample("D2Kpimumu_ValidationSubsample.root");
}

void D2KKmumuMC(){

  TChain* Tree_MC_D2KKmumu = new TChain("MC12_DstD2KKMuMu/DecayTree");
  Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/highStat_magDw/MC12_DstD2KKmumu_magDw.root");                                 
  Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/highStat_magUp/MC12_DstD2KKmumu_magUp.root");                                       Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/magDw/MC12_DstD2KKmumu_magDw.root");//"old one ""
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


void D2KpipipiData(){

  TChain* Tree_D2Kpipipi = new TChain("DstD2KPiPiPi/DecayTree");

  for (int i=0; i<450; ++i) {
    Tree_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magUp/2012Data_D2hhhh_%i.root",i));
    Tree_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magDw/2012Data_D2hhhh_%i.root",i));
  }


  D2KpipipiReader* Kpipipi_Reader = new D2KpipipiReader(Tree_D2Kpipipi);
  Kpipipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_PreselectedSubsample.root",10);
}


void D2KKpipiData(){

  TChain* Tree_D2KKpipi = new TChain("DstD2KKPiPi/DecayTree");

  for (int i=0; i<2000; ++i) {
    Tree_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/D2hhhhPIDline/magUp/2012Data_D2hhhh_%i.root",i));
    Tree_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/D2hhhhPIDline/magDw/2012Data_D2hhhh_%i.root",i));
  }

  D2KKpipiReader* KKpipi_Reader = new D2KKpipiReader(Tree_D2KKpipi);
  KKpipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_PreselectedSubsample.root",10);
}


void D2KKpipiMC(){


  TChain* Tree_MC_D2KKpipi = new TChain("MC12_DstD2KKpipi/DecayTree");

  for (int i=0; i<50; ++i) {
    Tree_MC_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKpipi_filtered/magDw/MC12_DstD2KKpipi_%i.root",i));
    Tree_MC_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKpipi_filtered/magUp/D2KKpipi_%i.root",i));
  }
 

  D2KKpipiReader* KK_MC_Reader = new D2KKpipiReader(Tree_MC_D2KKpipi);
  KK_MC_Reader->InitMC();
  KK_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_MCtrainingSample.root");
                                                                                                                                          
  //KK_MC_Reader->addMisIdMasses("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_filteredt_MCSample1.root");

}


void D2KpipipiMC(){


  TChain* Tree_MC_D2Kpipipi = new TChain("MC12_DstD2Kpipipi/DecayTree");

  for (int i=0; i<50; ++i) {
    Tree_MC_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpipipi_filtered/magDw/MC12_DstD2Kpipipi_%i.root",i));
    Tree_MC_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpipipi_filtered/magUp/D2Kpipipi_%i.root",i));
  }

  D2KpipipiReader* Kpi_MC_Reader = new D2KpipipiReader(Tree_MC_D2Kpipipi);
  Kpi_MC_Reader->InitMC();
  Kpi_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_MCtrainingSample.root");                            
}


