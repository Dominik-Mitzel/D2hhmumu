#include "optimizeSelection.h"
#include "D2hhmumuFitter.h"
#include "TEventList.h"
#include "TPaletteAxis.h"
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

  myFitter.fit_MC(true); //fix MC signal shape
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

