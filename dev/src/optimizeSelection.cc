#include "optimizeSelection.h"
#include "D2hhmumuFitter.h"
#include "TEventList.h"
#include "TPaletteAxis.h"
#include "D2KKmumuReader.h"
#include "D2pipimumuReader.h"
#include "D2KpimumuReader.h"
#include "D2KKpipiReader.h"
#include "D2KpipipiReader.h"
#include "D2pipipipiReader.h"
#include "TProfile.h"
#include "D2hhmumuFitter1D.h"
#include "TGenPhaseSpace.h"
#include "EfficiencyCalculator.h"
#include "EfficiencyStudies.h"
#include "TEfficiency.h"
#include "functions.C"
#include "D2hhmumuFitter1D.h"
#include "RooGenericPdf.h"

using namespace std;
using namespace RooFit ;

//scan PID BDT cut to find best working point by maximizing FOM
void optimizeSelection(TString cutNtracks="") { //do it for ntracks even and odd

  double factor = 2;
  TString nTrackSplitting="";
  if(cutNtracks=="nTracksEven") {nTrackSplitting="&&nTracks%2==0";factor=1;};
  if(cutNtracks=="nTracksOdd") {nTrackSplitting="&&nTracks%2!=0";factor=1;};

  std::cout<<"cutNtracks= " << cutNtracks << " nTrackSplitting =" << nTrackSplitting << " factor "<< factor<<std::endl;

  TFile *fout = new TFile("../rootFiles/PID_BDT_cut_scan_1D.root","recreate");
  TH2* h2= new TH2D("h2_cutEfficiency","h3_cutEfficiency",15,-0.5,1.0,10,0,1);
  TH2* h2_misIDBkg = new TH2D("h2_misIDBkg","h2_misIDBkg",15,-0.5,1,10,0,1);
  TH2* h2_combBkg = new TH2D("h2_combBkg","h2_combBkg",15,-0.5,1,10,0,1);
  TH2* h2_sigEff = new TH2D("h2_sigEff","h2_sigEff",15,-0.5,1,10,0,1);
  TH2* h2_relEff = new TH2D("h2_releff","h2_relEff",15,-0.5,1,10,0,1);


  //double PIDcuts[20]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95};
  //double BDTcuts[20]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95};
  //double PIDcuts[14]={0,0.1,0.2,0.3,0.4,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9};
  //double BDTcuts[14]={0,0.1,0.2,0.3,0.4,0.5,0.55,0.60,0.65,0.7,0.75,0.8,0.85,0.9};
  
  double PIDcuts[10]={0,0.1,0.2,0.3,0.4,0.5,0.60,0.7,0.8,0.9};                                                 
  //  double BDTcuts[10]={0,0.1,0.2,0.3,0.4,0.5,0.60,0.7,0.8,0.9};
  double BDTcuts[15]={-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  
  TString cut;
  TString cut_singleID;
  TString namePlot;

  double misIDBkg;
  double combBkg;
  double misIDBKG_scaled;
  double EffRatio=1;
  double signalEff=1;

  //set the files for the fitter. The fitter fixes the shapes for MC, misID bkg.
  D2hhmumuFitter1D myFitter;
  //myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root");
  //myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root");
  myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  myFitter.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  //try with MC /// myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root");   

  myFitter.fit_MC("deltaM>144.5 && deltaM<146.5"+nTrackSplitting,true,""); //fix MC signal shape
  double FOM;

  for(int i=0;i<15;++i){
     for(int j=0;j<10;++j){
  
      cut = TString::Format("BDT>%f&&mu0_ProbNNmu>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j],PIDcuts[j])+nTrackSplitting;
      namePlot =  TString::Format("BDT%fPID%f",BDTcuts[i],PIDcuts[j]);
      cut_singleID = TString::Format("BDT>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j]);


      myFitter.fit_Kpipipi_misID(cut_singleID,true,"");
      misIDBkg = myFitter.getMisIDbkgExp(cut,namePlot);
      combBkg = myFitter.getCombBkg(cut,namePlot);
      signalEff = getMCSignalEfficiency(cut)*2/factor;  //if splitting is done, factor = 1 and total efficiency has to be multiplied by 2
      EffRatio=EffD2KKpipiToEffD2Kpipipi(cut_singleID);
      misIDBKG_scaled= misIDBkg/(33.14)*EffRatio/factor*10/14; //factor 10/14 due to number of produced MC events 
      
      FOM=signalEff/( 5/2 + TMath::Sqrt(misIDBKG_scaled + combBkg) );

      h2->SetBinContent(i+1,j+1,FOM);
      h2_misIDBkg->SetBinContent(i+1,j+1,misIDBKG_scaled); 
      h2_combBkg ->SetBinContent(i+1,j+1,combBkg);
      h2_sigEff ->SetBinContent(i+1,j+1,signalEff);
      h2_relEff->SetBinContent(i+1,j+1,EffRatio);

      std::cout<<i<<" "<<j<<" FOM: "<<FOM<<"  signalEff:"<<signalEff << " msiID BKG: " <<misIDBKG_scaled <<" CombBkg: " << combBkg <<"  relEff  "<< EffRatio <<std::endl;

    }
  }

 TCanvas* c = new TCanvas("canvas","canvas"); 
	       	  
 c->Divide(2,2);
			  
  Float_t newMargin1 = 0.13;
  Float_t newMargin2 = 0.15;

  c->SetGrid();
  c->SetTicks();
  
  gStyle->SetPalette( 1, 0 );

  gStyle->SetPaintTextFormat( "3.3f" );
  gStyle->SetOptStat(0);

  h2->SetMarkerColor( 0 );
  c->cd(1);
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

  c->Write();
  h2->Write();
  h2_misIDBkg->Write();
  h2_combBkg ->Write();
  h2_sigEff ->Write();
  h2_relEff->Write();

  c->SaveAs("../img/selection1D.root");
  
  if(cutNtracks=="nTracksEven")c->Print("../img/selection1D_nTracksEven.eps");
  if(cutNtracks=="nTracksOdd") c->Print("../img/selection1D_nTracksOdd.eps");

  c->Print("../img/selection1D.eps");

  fout->Write();
 
}
void optimizeSelectionInBins_pipimumu(TString cutNtracks="",TString q2Range="") { //do it for ntracks even and odd

  double factor = 2;
  double combinatorialFactor =1/2; //because there are 2 possible ways to misidentify D2Kpipipi with respect to D2KKpipi
  TString cutQ2Range="";
  if(q2Range!="") cutQ2Range="&&"+q2Range; //just adds the && to the range
  TString nTrackSplitting="";
  TString splitting="";
  //translates NTracks even and odd to the corresponding mathematical expression
  if(cutNtracks=="nTracksEven") {nTrackSplitting="&&nTracks%2==0";splitting="nTracks%2==0";}; 
  if(cutNtracks=="nTracksOdd") {nTrackSplitting="&&nTracks%2!=0";splitting="nTracks%2!=0";};

  TFile *fout = new TFile("../rootFiles/pipimumu_PID_BDT_cut_scan_1D_"+q2Range+".root","recreate");
  TH2* h2= new TH2D("h2_cutEfficiency","h3_cutEfficiency",15,-0.5,1.0,10,0,1);
  TH2* h2_misIDBkg = new TH2D("h2_misIDBkg","h2_misIDBkg",15,-0.5,1,10,0,1);
  TH2* h2_combBkg = new TH2D("h2_combBkg","h2_combBkg",15,-0.5,1,10,0,1);
  TH2* h2_sigEff = new TH2D("h2_sigEff","h2_sigEff",15,-0.5,1,10,0,1);
  TH2* h2_relEff = new TH2D("h2_releff","h2_relEff",15,-0.5,1,10,0,1);

  double PIDcuts[10]={0,0.1,0.2,0.3,0.4,0.5,0.60,0.7,0.8,0.9};                                                 
  double BDTcuts[15]={-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  
  TString cut_doubleID;
  TString cutNormalizationChannel;
  TString cut_singleID;
  TString namePlot;

  double misIDBkg;
  double combBkg;
  double misIDBKG_scaled;
  double EffRatio=1;
  double signalEff=1;

  D2hhmumuFitter1D myFitter;
  myFitter.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2pipimumuBDT.root");
  myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_BDT.root");
  myFitter.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2pipimumuBDT.root");
  //myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  //here, the misID shape is taken from MC
  myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpipipi_D2pipimumuBDT.root");
  myFitter.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipipipi_PIDline_D2pipimumuBDT.root");   

  EfficiencyCalculator myEfficiency("D2pipimumu");

  myFitter.fit_normalization_MC("deltaM>144.5&&deltaM<146.5"+nTrackSplitting+cutQ2Range,true,""); //fix MC signal shape from fit
  double FOM;

  for(int i=0;i<15;++i){
   for(int j=0;j<10;++j){
  //for(int i=0;i<10;++i){
  //   for(int j=0;j<10;++j){
  
       ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
       /// perfom a scan in 2D space. Expected misID is calculated with respect to misID yield in norm channel. 
       /// Kpipipi misID is fitted in q2 range of normalzation channel and then extrapolated to corresponding bin in KKpipi 
       // nKKpipi = nKKpipi * rel.Efficiency * rel.BR * fractionOfq2Bin
       //////////////////////////////////////////////////////////////////////////////////

       
      cut_doubleID = TString::Format("BDT>%f&&mu0_ProbNNmu>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j],PIDcuts[j]);
      cut_singleID = TString::Format("BDT>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j]);
      cutNormalizationChannel = cut_doubleID+nTrackSplitting; //HERE: consider FULL dimuon spectrum, because BR of Kpipipi and KKpipi are for full known for full range
      std::cout<<" cut_doubleID " << cut_doubleID << " cut_singleID " << cut_singleID << "cutNormalizationChannel" << cutNormalizationChannel << "cutQ2Range " << cutQ2Range <<std::endl;
      
      namePlot ="newBDT_"+  TString::Format("BDT%fPID%f",BDTcuts[i],PIDcuts[j])+cutNtracks+q2Range;
      
      myFitter.fit_Kpipipi_misID(cut_singleID,true,""); //no trck splitting here, just to get shape
      misIDBkg = myFitter.getMisIDbkgExp(cutNormalizationChannel,namePlot) * factor/combinatorialFactor; //account for fact that only half of data is used, no PID for misID accounted for; 
      combBkg = myFitter.getCombBkg(cut_doubleID+nTrackSplitting+cutQ2Range,namePlot) * factor;
      misIDBKG_scaled = misIDBkg/24.52 * myEfficiency.getMCRelativeSigToNormMisIDEfficiency(TString::Format("BDT>%f",BDTcuts[i]),splitting,q2Range,"D_DiMuon_Mass>0"/*full range*/) * myEfficiency.getMisIDFractionQ2Range(q2Range);
      signalEff = myEfficiency.getMCSignalEfficiency(cut_doubleID,splitting,q2Range);
      
      FOM=signalEff/( 5/2 + TMath::Sqrt(misIDBKG_scaled + combBkg) );

      h2->SetBinContent(i+1,j+1,FOM);
      h2_misIDBkg->SetBinContent(i+1,j+1,misIDBKG_scaled); 
      h2_combBkg ->SetBinContent(i+1,j+1,combBkg);
      h2_sigEff ->SetBinContent(i+1,j+1,signalEff);
      h2_relEff->SetBinContent(i+1,j+1,EffRatio);

      std::cout<<i<<" "<<j<<" FOM: "<<FOM<<"  signalEff:"<<signalEff << " msiID BKG: " <<misIDBKG_scaled <<" CombBkg: " << combBkg  <<std::endl;

    }
  }

 TCanvas* c = new TCanvas("canvas","canvas"); 
	       	  
 c->Divide(2,2);
			  
  Float_t newMargin1 = 0.13;
  Float_t newMargin2 = 0.15;

  c->SetGrid();
  c->SetTicks();
  
  gStyle->SetPalette( 1, 0 );

  gStyle->SetPaintTextFormat( "3.3f" );
  gStyle->SetOptStat(0);

  h2->SetMarkerColor( 0 );
  c->cd(1);
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

  c->Write();
  h2->Write();
  h2_misIDBkg->Write();
  h2_combBkg ->Write();
  h2_sigEff ->Write();
  h2_relEff->Write();

  //c->SaveAs("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/selectionOptimisation/selection1D_"+q2Range+".root");
  //c->Print("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/selectionOptimisation/selection1D_"+q2Range+".eps");
  
  if(cutNtracks=="nTracksEven")c->Print("/work/mitzel/D2hhmumu/dev/D2pipimumu/img/selectionOptimisation/newBDT_selection1D_nTracksEven_"+q2Range+".eps");
  if(cutNtracks=="nTracksOdd") c->Print("/work/mitzel/D2hhmumu/dev/D2pipimumu/img/selectionOptimisation/newBDT_selection1D_nTracksOdd_"+q2Range+".eps");


  fout->Write();
 
}
void optimizeSelectionInBins(TString cutNtracks="",TString q2Range="") { //do it for ntracks even and odd

  double factor = 2;
  double combinatorialFactor =2; //because there are 2 possible ways to misidentify D2Kpipipi with respect to D2KKpipi
  TString cutQ2Range="";
  if(q2Range!="") cutQ2Range="&&"+q2Range; //just adds the && to the range
  TString nTrackSplitting="";
  TString splitting="";
  //translates NTracks even and odd to the corresponding mathematical expression
  if(cutNtracks=="nTracksEven") {nTrackSplitting="&&nTracks%2==0";splitting="nTracks%2==0";}; 
  if(cutNtracks=="nTracksOdd") {nTrackSplitting="&&nTracks%2!=0";splitting="nTracks%2!=0";};

  TFile *fout = new TFile("../rootFiles/PID_BDT_cut_scan_1D_"+q2Range+".root","recreate");
  TH2* h2= new TH2D("h2_cutEfficiency","h3_cutEfficiency",15,-0.5,1.0,10,0,1);
  TH2* h2_misIDBkg = new TH2D("h2_misIDBkg","h2_misIDBkg",15,-0.5,1,10,0,1);
  TH2* h2_combBkg = new TH2D("h2_combBkg","h2_combBkg",15,-0.5,1,10,0,1);
  TH2* h2_sigEff = new TH2D("h2_sigEff","h2_sigEff",15,-0.5,1,10,0,1);
  TH2* h2_relEff = new TH2D("h2_releff","h2_relEff",15,-0.5,1,10,0,1);

  double PIDcuts[10]={0,0.1,0.2,0.3,0.4,0.5,0.60,0.7,0.8,0.9};                                                 
  double BDTcuts[15]={-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  
  TString cut_doubleID;
  TString cutNormalizationChannel;
  TString cut_singleID;
  TString namePlot;

  double misIDBkg;
  double combBkg;
  double misIDBKG_scaled;
  double EffRatio=1;
  double signalEff=1;

  D2hhmumuFitter1D myFitter;
  myFitter.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  myFitter.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  //myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  //here, the misID shape is taken from MC
  myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpipipi_D2KKmumuBDT.root");
  myFitter.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root");   

  EfficiencyCalculator myEfficiency("D2KKmumu");

  myFitter.fit_normalization_MC("deltaM>144.5&&deltaM<146.5"+nTrackSplitting+cutQ2Range,true,""); //fix MC signal shape from fit
  double FOM;

  for(int i=0;i<15;++i){
   for(int j=0;j<10;++j){
  //for(int i=0;i<10;++i){
  //   for(int j=0;j<10;++j){
  
       ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
       /// perfom a scan in 2D space. Expected misID is calculated with respect to misID yield in norm channel. 
       /// Kpipipi misID is fitted in q2 range of normalzation channel and then extrapolated to corresponding bin in KKpipi 
       // nKKpipi = nKKpipi * rel.Efficiency * rel.BR * fractionOfq2Bin
       //////////////////////////////////////////////////////////////////////////////////

       
      cut_doubleID = TString::Format("BDT>%f&&mu0_ProbNNmu>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j],PIDcuts[j]);
      cut_singleID = TString::Format("BDT>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j]);
      cutNormalizationChannel = cut_doubleID+nTrackSplitting; //HERE: consider FULL dimuon spectrum, because BR of Kpipipi and KKpipi are for full known for full range
      std::cout<<" cut_doubleID " << cut_doubleID << " cut_singleID " << cut_singleID << "cutNormalizationChannel" << cutNormalizationChannel << "cutQ2Range " << cutQ2Range <<std::endl;
      
      namePlot ="newBDT_"+  TString::Format("BDT%fPID%f",BDTcuts[i],PIDcuts[j])+cutNtracks+q2Range;
      
      myFitter.fit_Kpipipi_misID(cut_singleID,true,""); //no trck splitting here, just to get shape
      misIDBkg = myFitter.getMisIDbkgExp(cutNormalizationChannel,namePlot) * factor/combinatorialFactor; //account for fact that only half of data is used, no PID for misID accounted for; 
      combBkg = myFitter.getCombBkg(cut_doubleID+nTrackSplitting+cutQ2Range,namePlot) * factor;
      misIDBKG_scaled = misIDBkg/33.14 * myEfficiency.getMCRelativeSigToNormMisIDEfficiency(TString::Format("BDT>%f",BDTcuts[i]),splitting,q2Range,"D_DiMuon_Mass>0"/*full range*/) * myEfficiency.getMisIDFractionQ2Range(q2Range);
      signalEff = myEfficiency.getMCSignalEfficiency(cut_doubleID,splitting,q2Range);
      
      FOM=signalEff/( 5/2 + TMath::Sqrt(misIDBKG_scaled + combBkg) );

      h2->SetBinContent(i+1,j+1,FOM);
      h2_misIDBkg->SetBinContent(i+1,j+1,misIDBKG_scaled); 
      h2_combBkg ->SetBinContent(i+1,j+1,combBkg);
      h2_sigEff ->SetBinContent(i+1,j+1,signalEff);
      h2_relEff->SetBinContent(i+1,j+1,EffRatio);

      std::cout<<i<<" "<<j<<" FOM: "<<FOM<<"  signalEff:"<<signalEff << " msiID BKG: " <<misIDBKG_scaled <<" CombBkg: " << combBkg  <<std::endl;

    }
  }

 TCanvas* c = new TCanvas("canvas","canvas"); 
	       	  
 c->Divide(2,2);
			  
  Float_t newMargin1 = 0.13;
  Float_t newMargin2 = 0.15;

  c->SetGrid();
  c->SetTicks();
  
  gStyle->SetPalette( 1, 0 );

  gStyle->SetPaintTextFormat( "3.3f" );
  gStyle->SetOptStat(0);

  h2->SetMarkerColor( 0 );
  c->cd(1);
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

  c->Write();
  h2->Write();
  h2_misIDBkg->Write();
  h2_combBkg ->Write();
  h2_sigEff ->Write();
  h2_relEff->Write();

  //c->SaveAs("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/selectionOptimisation/selection1D_"+q2Range+".root");
  //c->Print("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/selectionOptimisation/selection1D_"+q2Range+".eps");
  
  if(cutNtracks=="nTracksEven")c->Print("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/selectionOptimisation/newBDT_selection1D_nTracksEven_"+q2Range+".eps");
  if(cutNtracks=="nTracksOdd") c->Print("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/selectionOptimisation/newBDT_selection1D_nTracksOdd_"+q2Range+".eps");


  fout->Write();
 
}

void optimizeSelection2D(TString cutNtracks="") { //do it for ntracks even and odd

  double factor = 2;
  TString nTrackSplitting="";
  if(cutNtracks=="nTracksEven") {nTrackSplitting="&&nTracks%2==0";factor=1;};
  if(cutNtracks=="nTracksOdd") {nTrackSplitting="&&nTracks%2!=0";factor=1;};

  std::cout<<"cutNtracks= " << cutNtracks << " nTrackSplitting =" << nTrackSplitting << " factor "<< factor<<std::endl;

  TFile *fout = new TFile("../rootFiles/PID_BDT_cut_scan_2D.root","recreate");
  TH2* h2= new TH2D("h2_cutEfficiency","h3_cutEfficiency",15,-0.5,1.0,10,0,1);
  TH2* h2_misIDBkg = new TH2D("h2_misIDBkg","h2_misIDBkg",15,-0.5,1,10,0,1);
  TH2* h2_combBkg = new TH2D("h2_combBkg","h2_combBkg",15,-0.5,1,10,0,1);
  TH2* h2_sigEff = new TH2D("h2_sigEff","h2_sigEff",15,-0.5,1,10,0,1);
  TH2* h2_relEff = new TH2D("h2_releff","h2_relEff",15,-0.5,1,10,0,1);
  
  double PIDcuts[10]={0,0.1,0.2,0.3,0.4,0.5,0.60,0.7,0.8,0.9};                                                 
  double BDTcuts[15]={-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  
  TString cut;
  TString cut_singleID;
  TString namePlot;

  double misIDBkg;
  double combBkg;
  double misIDBKG_scaled;
  double EffRatio=1;
  double signalEff=1;

  D2hhmumuFitter myFitter; //2D fitter
  myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  myFitter.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  //try with MC /// myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples_backup/MC_D2Kpipipi_D2KKmumuBDT_temp.root");
  //myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");   
  
  std::cout<<"fixing signal MC shape for cut "<< cut << std::endl;
  myFitter.fit_MC("",true,""); //fix MC signal shape
  double FOM;

  for(int i=0;i<15;++i){
     for(int j=0;j<10;++j){
 
      cut = TString::Format("BDT>%f&&mu0_ProbNNmu>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j],PIDcuts[j])+nTrackSplitting;
      namePlot =  TString::Format("BDT%fPID%f",BDTcuts[i],PIDcuts[j]);
      cut_singleID = TString::Format("BDT>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j]);

      //      std::cout<<cut<<std::endl;
      std::cout<<"fixing misiID shape for cut "<< cut << std::endl;
      myFitter.fit_Kpipipi_misID(cut_singleID,true,"");
      std::cout<<"extracting misID bkg..." << std::endl;
      misIDBkg = myFitter.getMisIDbkgExp(cut,namePlot);
      std::cout<<"..and comb bkg... "<< std::endl;
      combBkg = myFitter.getCombBkg(cut,namePlot);
      std::cout<<"..computing signal efficiency... "<< std::endl;
      signalEff = getMCSignalEfficiency(cut)*2/factor;  //if splitting is done, factor = 1 and total efficiency has to be multiplied by 2
      std::cout<<"..and efficinecy ratio for hadronic modes... "<< std::endl;
      EffRatio=EffD2KKpipiToEffD2Kpipipi(cut_singleID);
      misIDBKG_scaled= misIDBkg/(33.14)*EffRatio/factor*10/14; //factor 10/14 due to number of produced MC events 
      
      FOM=signalEff/( 5/2 + TMath::Sqrt(misIDBKG_scaled + combBkg) );

      h2->SetBinContent(i+1,j+1,FOM);
      h2_misIDBkg->SetBinContent(i+1,j+1,misIDBKG_scaled); 
      h2_combBkg ->SetBinContent(i+1,j+1,combBkg);
      h2_sigEff ->SetBinContent(i+1,j+1,signalEff);
      h2_relEff->SetBinContent(i+1,j+1,EffRatio);

      std::cout<<i<<" "<<j<<" FOM: "<<FOM<<"  signalEff:"<<signalEff << " msiID BKG: " <<misIDBKG_scaled <<" CombBkg: " << combBkg <<"  relEff  "<< EffRatio <<std::endl;

    }
  }

 TCanvas* c = new TCanvas("canvas","canvas"); 
	       	  
 c->Divide(2,2);
			  
  Float_t newMargin1 = 0.13;
  Float_t newMargin2 = 0.15;

  c->SetGrid();
  c->SetTicks();
  
  gStyle->SetPalette( 1, 0 );

  gStyle->SetPaintTextFormat( "3.3f" );
  gStyle->SetOptStat(0);

  h2->SetMarkerColor( 0 );
  c->cd(1);
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

  c->Write();
  h2->Write();
  h2_misIDBkg->Write();
  h2_combBkg ->Write();
  h2_sigEff ->Write();
  h2_relEff->Write();

  c->SaveAs("../img/selection2D.root");
  
  if(cutNtracks=="nTracksEven")c->Print("../img/selection2D_nTracksEven.eps");
  if(cutNtracks=="nTracksOdd") c->Print("../img/selection2D_nTracksOdd.eps");

  c->Print("../img/selection2D.eps");

  fout->Write();

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
  TString pathToKKpipiMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKpipi_D2KKmumuBDT.root";
  TString pathToKpipipiMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpipipi_D2KKmumuBDT.root"; //!!!!!! set to temp file due to server error
 
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



//the following functions are used to merge tuples into one tuple, activate only relevant branches and add some useful variables. 

void D2KKmumuData(){


  TChain* Tree_D2KKmumu = new TChain("DstD2KKMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2KKmumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));
    Tree_D2KKmumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }
  
  D2KKmumuReader* KK_Reader = new D2KKmumuReader(Tree_D2KKmumu);
  KK_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_PreselectedSubsample.root",100);
}

void D2pipimumuData(){

  TChain* Tree_D2pipimumu = new TChain("DstD2PiPiMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2pipimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));
    Tree_D2pipimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }

  D2pipimumuReader* pipi_Reader = new D2pipimumuReader(Tree_D2pipimumu);
  pipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_PreselectedSubsample.root",100);

}

void D2KpimumuData(){

  TChain* Tree_D2Kpimumu = new TChain("DstD2KPiMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2Kpimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));                                                                 
    Tree_D2Kpimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }

  D2KpimumuReader* Kpi_Reader = new D2KpimumuReader(Tree_D2Kpimumu);
  Kpi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_PreselectedSubsample.root",0);
}

void D2KKmumuMC(){

  TChain* Tree_MC_D2KKmumu = new TChain("MC12_DstD2KKMuMu/DecayTree");

  
  for (int i=0; i<100; ++i) {
    Tree_MC_D2KKmumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_withGenLevel/magDw/MC12_DstD2KKmumu_%i.root",i));
    Tree_MC_D2KKmumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_withGenLevel/magUp/MC12_DstD2KKmumu_%i.root",i));
  }

  D2KKmumuReader* KK_MC_Reader = new D2KKmumuReader(Tree_MC_D2KKmumu);
  KK_MC_Reader->InitMC();
  KK_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_MCtrainingSample.root");

  createGeneratorLevelMCTuple("D2KKmumu");

}


void D2pipimumuMC(){

  TChain* Tree_MC_D2pipimumu = new TChain("MC12_DstD2PiPiMuMu/DecayTree");
  Tree_MC_D2pipimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu_withGenLevel/magDw/MC12_DstD2pipimumu_magDw.root");
  Tree_MC_D2pipimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu_withGenLevel/magUp/MC12_DstD2pipimumu_magUp.root");

  D2pipimumuReader* pipi_MC_Reader = new D2pipimumuReader(Tree_MC_D2pipimumu);
  pipi_MC_Reader->InitMC();
  pipi_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_MCtrainingSample.root");
  createGeneratorLevelMCTuple("D2pipimumu");
                                                                                                   
}

void D2KpimumuMC(){


  TChain* Tree_MC_D2Kpimumu = new TChain("MC12_DstD2KPiMuMu/DecayTree");
  //old MC
  //  Tree_MC_D2Kpimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu/magDw/MC12_DstD2Kpimumu_magDw.root");
  //Tree_MC_D2Kpimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu/magUp/MC12_DstD2Kpimumu_magUp.root");

  for (int i=0; i<100; ++i) {
    Tree_MC_D2Kpimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu_withGenLevel/magDw/MC12_DstD2Kpimumu_%i.root",i));
    Tree_MC_D2Kpimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu_withGenLevel/magUp/MC12_DstD2Kpimumu_%i.root",i));
  }

  //new hight Stat, no GenLevel info
  //Tree_MC_D2Kpimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu/highStat_magDw/MC12_DstD2Kpimumu_magDw.root");
  //Tree_MC_D2Kpimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu/highStat_magUp/MC12_DstD2Kpimumu_magUp.root");

  D2KpimumuReader* Kpi_MC_Reader = new D2KpimumuReader(Tree_MC_D2Kpimumu);
  Kpi_MC_Reader->InitMC();
  Kpi_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_MCtrainingSample.root");
  createGeneratorLevelMCTuple("D2Kpimumu");
}


void D2KpipipiData(){

  TChain* Tree_D2Kpipipi = new TChain("DstD2KPiPiPi/DecayTree");

  for (int i=0; i<1500; ++i) {
    Tree_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magUp/2012Data_D2hhhh_%i.root",i));
    Tree_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magDw/2012Data_D2hhhh_%i.root",i));
  }


  D2KpipipiReader* Kpipipi_Reader = new D2KpipipiReader(Tree_D2Kpipipi);
  Kpipipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_PreselectedSubsample.root",0);
}

void D2KpipipiRandomizedData(){

  TChain* Tree_D2Kpipipi = new TChain("DstD2KPiPiPi/DecayTree");

  for (int i=0; i<1500; ++i) {
    Tree_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magUp/2012Data_D2hhhh_%i.root",i));
    Tree_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magDw/2012Data_D2hhhh_%i.root",i));
  }

  //ADDED TRIGGER INFOS HERE 27.4.17
  D2KpipipiReader* Kpipipi_Reader = new D2KpipipiReader(Tree_D2Kpipipi);
  Kpipipi_Reader->createRandomizedSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_PreselectedRandomizedSubsample_withTrigger.root");
}


void D2KKpipiData(){

  TChain* Tree_D2KKpipi = new TChain("DstD2KKPiPi/DecayTree");

  for (int i=0; i<1500; ++i) {
    Tree_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magUp/2012Data_D2hhhh_%i.root",i));
    Tree_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magDw/2012Data_D2hhhh_%i.root",i));
  }

  D2KKpipiReader* KKpipi_Reader = new D2KKpipiReader(Tree_D2KKpipi);
  KKpipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_PreselectedSubsample.root",10);
}

void D2pipipipiData(){

  TChain* Tree_D2pipipipi = new TChain("DstD2PiPiPiPi/DecayTree");

  for (int i=0; i<1500; ++i) {
    Tree_D2pipipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magUp/2012Data_D2hhhh_%i.root",i));
    Tree_D2pipipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magDw/2012Data_D2hhhh_%i.root",i));
  }

  D2pipipipiReader* pipipipi_Reader = new D2pipipipiReader(Tree_D2pipipipi);
  pipipipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_PreselectedSubsample.root",10);
}

void D2pipipipiRandomizedData(){

  TChain* Tree_D2pipipipi = new TChain("DstD2PiPiPiPi/DecayTree");

  for (int i=0; i<1500; ++i) {
    Tree_D2pipipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magUp/2012Data_D2hhhh_%i.root",i));
    Tree_D2pipipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magDw/2012Data_D2hhhh_%i.root",i));
  }

  D2pipipipiReader* pipipipi_Reader = new D2pipipipiReader(Tree_D2pipipipi);
  pipipipi_Reader->createRandomizedSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_PreselectedRandomizedSubsample.root");

}


void D2pipipipiMC(){

  TChain* Tree_D2pipipipi = new TChain("MC12_DstD2pipipipi/DecayTree");

  Tree_D2pipipipi->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipipipi_withGenLevel/magUp/MC12_DstD2pipipipi_magUp.root");
  Tree_D2pipipipi->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipipipi_withGenLevel/magUp/MC12_DstD2pipipipi_magDw.root");

  D2pipipipiReader* pipipipi_Reader = new D2pipipipiReader(Tree_D2pipipipi);
  pipipipi_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_MCtrainingSample.root");

  createGeneratorLevelMCTuple("D2pipipipi");

}


void D2KKpipiMC(){


  TChain* Tree_MC_D2KKpipi = new TChain("MC12_DstD2KKpipi/DecayTree");

  for (int i=0; i<100; ++i) {
    Tree_MC_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKpipi_withGenLevel/magDw/MC12_DstD2KKpipi_%i.root",i));
    Tree_MC_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKpipi_withGenLevel/magUp/MC12_DstD2KKpipi_%i.root",i));
  }
 

  D2KKpipiReader* KK_MC_Reader = new D2KKpipiReader(Tree_MC_D2KKpipi);
  KK_MC_Reader->InitMC();
  KK_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_MCtrainingSample.root");
                                                                                                                                          
  //KK_MC_Reader->addMisIdMasses("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_filteredt_MCSample1.root");
  createGeneratorLevelMCTuple("D2KKpipi");
}


void D2KpipipiMC(){


  TChain* Tree_MC_D2Kpipipi = new TChain("MC12_DstD2Kpipipi/DecayTree");

  for (int i=0; i<100; ++i) {
    Tree_MC_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpipipi_withGenLevel/magDw/MC12_DstD2Kpipipi_%i.root",i));
    Tree_MC_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpipipi_withGenLevel/magUp/MC12_DstD2Kpipipi_%i.root",i));
  }

  D2KpipipiReader* Kpi_MC_Reader = new D2KpipipiReader(Tree_MC_D2Kpipipi);
  Kpi_MC_Reader->InitMC();
  Kpi_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_MCtrainingSample.root");                            
  
  //create also generator level tuple with same data 
  createGeneratorLevelMCTuple("D2Kpipipi");

}


//draw plots to check correlation of BDT with mass variables
void draw_BDT_crosschecks(){

  draw_BDT_crosschecks_forVariable("Dst_DTF_D0_M",1840,1950,0.7);
  draw_BDT_crosschecks_forVariable("deltaM",140,155,0.7);
  draw_BDT_crosschecks_forVariable("D_DiMuon_Mass",200,900,0.7);
  draw_BDT_crosschecks_forVariable("mHH",950,1450,0.7);

}

    


void draw_BDT_crosschecks_forVariable(TString variable, int xLow, int xHigh, double BDTcut){

  TH1::SetDefaultSumw2();

  
  TString outname="../img/BDTcrosschecks/"+variable+".root" ;
  TChain* Chain_MC_sig = new TChain("BDT_Tree");
  TChain* Chain_MC_norm = new TChain("BDT_Tree");
  TChain* Chain_data_sideband = new TChain("BDT_Tree");
  
  TString inputfile_MC_sig= "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root";
  TString inputfile_MC_Norm= "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root";
  TString inputfile_data_sideband = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT.root";

  Chain_MC_sig->AddFile(inputfile_MC_sig);
  Chain_MC_norm->AddFile(inputfile_MC_Norm);
  Chain_data_sideband->AddFile(inputfile_data_sideband);

  TFile *outfile = new TFile(outname,"recreate");

  TCanvas *c1_MC_sig = new TCanvas("c1_MC_sig","c1_MC_sig",200,10,700,500);
  TCanvas *c1_MC_norm = new TCanvas("c1_MC_norm","c1_MC_norm",200,10,700,500);
  TCanvas *c1_MC_sideband = new TCanvas("c1_data_sideband","c1_data_sidenand",200,10,700,500);
  TCanvas *c1_ratio = new TCanvas("eff ratio","effRatio",200,10,700,500);


  TH2* h2_MC_sig = new TH2D ("h2_MC_sig_"+variable,"BDT versus "+variable+" signal MC",100,xLow,xHigh,100,-1,1);
  TH2* h2_MC_norm = new TH2D ("h2_MC_norm_"+variable,"BDT versus "+variable+" snormalisation MC",100,xLow,xHigh,100,-1,1);
  TH2* h2_sideband = new TH2D ("h2_sideband_sig_"+variable,"BDT versus "+variable+" data sideband",100,xLow,xHigh,100,-1,1);

  TH1* h1_BDTeff_sig = new TH1D ("h1_MC_sig_"+variable,"BDT versus "+variable+" signal MC",20,xLow,xHigh);
  TH1* h1_BDTeff_norm = new TH1D ("h1_MC_norm_"+variable,"BDT versus "+variable+" normalisation MC",20,xLow,xHigh);
  TH1* h1_BDTeff_sideband = new TH1D ("h1_MC_sideband_"+variable,"BDT versus "+variable+" data sideband",20,xLow,xHigh);

  TH1* h1_BDTeff_sig_total = new TH1D ("h1_MC_sig_"+variable+"_total","BDT versus "+variable+" signal MC",20,xLow,xHigh);
  TH1* h1_BDTeff_norm_total = new TH1D ("h1_MC_norm_"+variable+"_total","BDT versus "+variable+" normalisation MC",20,xLow,xHigh);
  TH1* h1_BDTeff_sideband_total = new TH1D ("h1_MC_sideband_"+variable+"_total","BDT versus "+variable+" data sideband",20,xLow,xHigh);

  std::vector<TH2*> h2s;
  std::vector<TH1*> h1s;
  std::vector<TH1*> h1norms;
  std::vector<TCanvas*> canvases;
  std::vector<TChain*> chains;

  canvases.push_back(c1_MC_sig);
  canvases.push_back(c1_MC_norm);
  canvases.push_back(c1_MC_sideband);

  h2s.push_back(h2_MC_sig);
  h2s.push_back(h2_MC_norm);
  h2s.push_back(h2_sideband);

  h1s.push_back(h1_BDTeff_sig);
  h1s.push_back(h1_BDTeff_norm);
  h1s.push_back(h1_BDTeff_sideband);

  h1norms.push_back(h1_BDTeff_sig_total);
  h1norms.push_back(h1_BDTeff_norm_total);
  h1norms.push_back(h1_BDTeff_sideband_total);

  chains.push_back(Chain_MC_sig);
  chains.push_back(Chain_MC_norm);
  chains.push_back(Chain_data_sideband);
  
  double BDT, temp;

  for(int i=0;i<3;++i) {

      chains[i]->SetBranchAddress("BDT",&BDT);
      chains[i]->SetBranchAddress(variable,&temp);

      Int_t nentries = chains[i]->GetEntries();
      for (Int_t j=0;j<nentries;++j) { // loop over tree entries 
	chains[i]->GetEntry(j);
	h2s[i]->Fill(temp,BDT);
	h1norms[i]->Fill(temp);
	//	std::cout<<BDT<<std::endl;
	if(BDT>BDTcut) h1s[i]->Fill(temp); 
      }
      
      canvases[i]->Divide(2);
      canvases[i]->cd(1);
      h2s[i]->Draw();
      h2s[i]->ProfileX()->Draw("same");
      canvases[i]->cd(2);
      h1s[i]->Divide(h1norms[i]);
      h1s[i]->GetXaxis()->SetRangeUser(0,1);
      h1s[i]->Draw();	
      canvases[i]->Write();
    }   
    
    c1_ratio->cd();
    TH1D* effSig = (TH1D*)h1s[0]->Clone();
    TH1D* effNorm = (TH1D*)h1s[1]->Clone();
    effSig->Divide(effNorm);
    effSig->Draw();
    c1_ratio->Write();
    
    outfile->Write();
}

void check_peakingBackground(TString kind,bool PIDCut=false) {

  TChain* fChain = new TChain("BDT_Tree");

  if(kind!= "D2KKmumu") std::cout<<" Selection optimization: Error decay mode not implemeted"<<std::endl;
  if(kind== "D2KKmumu"){ fChain->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");}

  TH1* h_mKpipipi = new TH1D ("h_mKpipipi","h_mKpipipi",50,1650,1950);
  TH1* h_mpiKpipi = new TH1D ("h_mpiKpipi","h_mpiKpipi",50,1650,1950);
  TH2* h_mpiKpipi_mD = new TH2D ("h_mpiKpipi_vs_mD","h_mpiKpipi_vs_mD",20,1650,1950,20,1800,1950);

  TH1* h_mKpimumu = new TH1D ("h_mKpimumu","h_mKpimumu",50,1650,1950);
  TH1* h_mpiKmumu = new TH1D ("h_mpiKmumu","h_mpiKmumu",50,1650,1950);
  TH1* h_mKKmumu = new TH1D ("h_mKKmumu","h_mKKmumu",50,1750,2000);

  TH1* h_mKpipipi_asKKmumu = new TH1D ("h_mKpipipi_asKKmumu","h_mKpipipi_asKKmumu",50,1900,2350);
  TH1* h_mKpimumu_asKKmumu = new TH1D ("h_mKpimumu_asKKmumu","h_mKpiMmumu_asKKmumu",50,1900,2350);
  TH1* h_mKKmumu_asKpipipi = new TH1D ("h_mKKmumu_asKpipipi","h_mKKmumu_asKpipipi",50,1600,1850);
  TH1* h_mKpipipipi0_asKKmumu = new TH1D ("h_mKpipipipi0_asKKmumu","h_mKpipipipi0_asKKmumu",50,1200,2150);
 
 
  double h0_PX,h0_PY,h0_PZ;
  double h1_PX,h1_PY,h1_PZ;
  double mu0_PX,mu0_PY,mu0_PZ;
  double mu1_PX,mu1_PY,mu1_PZ;
  double deltaM, mD;
  
  double mu0_ProbNNmu,mu1_ProbNNmu,h1_ProbNNK,h0_ProbNNK;
  double BDT;

  fChain->SetBranchAddress("h1_PX",&h1_PX);
  fChain->SetBranchAddress("h1_PY",&h1_PY);
  fChain->SetBranchAddress("h1_PZ",&h1_PZ);
  fChain->SetBranchAddress("h0_PX",&h0_PX);
  fChain->SetBranchAddress("h0_PY",&h0_PY);
  fChain->SetBranchAddress("h0_PZ",&h0_PZ);
  fChain->SetBranchAddress("mu0_PX",&mu0_PX);
  fChain->SetBranchAddress("mu0_PY",&mu0_PY);
  fChain->SetBranchAddress("mu0_PZ",&mu0_PZ);
  fChain->SetBranchAddress("mu1_PX",&mu1_PX);
  fChain->SetBranchAddress("mu1_PY",&mu1_PY);
  fChain->SetBranchAddress("mu1_PZ",&mu1_PZ);
  fChain->SetBranchAddress("deltaM",&deltaM);
  fChain->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  fChain->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
  fChain->SetBranchAddress("BDT",&BDT);
  fChain->SetBranchAddress("Dst_DTF_D0_M",&mD);
  fChain->SetBranchAddress("h1_ProbNNk",&h1_ProbNNK);
  fChain->SetBranchAddress("h0_ProbNNk",&h0_ProbNNK);

  TLorentzVector h0_hypK,h0_hypPi;
  TLorentzVector h1_hypK,h1_hypPi;
  TLorentzVector mu0_hypK,mu0_hypPi,mu0_hypMu;
  TLorentzVector mu1_hypK,mu1_hypPi,mu1_hypMu;
  TString fileName = "/work/mitzel/D2hhmumu/dev/D2KKmumu/img/peakingBkgStudies/peakingBkg_studies_noPIDCut.eps";  

  for(int i=0; i<fChain->GetEntries();++i) {

    fChain->GetEntry(i);

    if(BDT<0.6) continue;
    if(PIDCut) {
      if(mu0_ProbNNmu<0.5) continue;
      if(mu1_ProbNNmu<0.5) continue;
      fileName = "/work/mitzel/D2hhmumu/dev/D2KKmumu/img/peakingBkgStudies/peakingBkg_studies.eps";  
      
    }
    //if(h1_ProbNNK<0.5 || h0_ProbNNK<0.5) continue;
    //if(deltaM<144.5 || deltaM >146.5) continue;
    if(deltaM>144. && deltaM <147) continue; //cut signal in deltaM
    //if(mD<1880&&mD>1840) continue;
    

    h0_hypK.SetXYZM(h0_PX,h0_PY,h0_PZ,Mass::K());
    h0_hypPi.SetXYZM(h0_PX,h0_PY,h0_PZ,Mass::Pi());
    h1_hypK.SetXYZM(h1_PX,h1_PY,h1_PZ,Mass::K());
    h1_hypPi.SetXYZM(h1_PX,h1_PY,h1_PZ,Mass::Pi());
    mu0_hypK.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::K());
    mu0_hypPi.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Pi());
    mu0_hypMu.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1_hypK.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::K());
    mu1_hypPi.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Pi());
    mu1_hypMu.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());

    h_mKpipipi->Fill((h0_hypK+h1_hypPi+ mu0_hypPi + mu1_hypPi) .M()); 
    h_mpiKpipi->Fill((h0_hypPi+h1_hypK+ mu0_hypPi + mu1_hypPi) .M());
    h_mKpimumu->Fill((h0_hypK+h1_hypPi+ mu0_hypMu + mu1_hypMu) .M());
    h_mpiKmumu->Fill((h0_hypPi+h1_hypK+ mu0_hypMu + mu1_hypMu) .M());
    h_mKKmumu->Fill((h0_hypK+h1_hypK+ mu0_hypMu + mu1_hypMu) .M());
    //if((h0_hypPi+h1_hypK+ mu0_hypPi + mu1_hypPi).M()<1800) { h_mKKmumu->Fill(mD); }

    h_mpiKpipi_mD->Fill((h0_hypPi+h1_hypK+ mu0_hypPi + mu1_hypPi).M(),mD);
   
  }

  TCanvas *canvas = new TCanvas("misID shapes","miIDshapes");
  canvas->Divide(3,3);
  canvas->cd(1);
  h_mKpipipi->Draw();
  canvas->cd(2);
  h_mpiKpipi->Draw();
  canvas->cd(3);
  h_mKpimumu->Draw();
  canvas->cd(4);
  h_mpiKmumu->Draw();
  canvas->cd(7);
  h_mpiKpipi_mD->Draw("colz");
  canvas->cd(9);
  h_mKKmumu->Draw();

  
  TLorentzVector pD(0.0, 0.0, 0.0, Mass::D0());
  Double_t masses1[4] = {Mass::K(),Mass::Pi(),Mass::Pi(),Mass::Pi()} ;
  Double_t masses2[4] = {Mass::K(),Mass::Pi(),Mass::Mu(),Mass::Mu()} ;
  Double_t masses3[4] = {Mass::K(),Mass::K(),Mass::Mu(),Mass::Mu()} ;
  Double_t masses4[5] = {Mass::K(),Mass::Pi(),Mass::Pi(),Mass::Pi(),Mass::Pi0()} ;
  

  TGenPhaseSpace event_Kpipipi;
  event_Kpipipi.SetDecay(pD, 4, masses1);
  TGenPhaseSpace event_Kpimumu;
  event_Kpimumu.SetDecay(pD, 4, masses2);
  TGenPhaseSpace event_KKmumu;
  event_KKmumu.SetDecay(pD, 4, masses3);
  TGenPhaseSpace event_Kpipipipi0;
  event_Kpipipipi0.SetDecay(pD, 5, masses4);

  for (Int_t n=0;n<100000;n++) {
    Double_t weight1 = event_Kpipipi.Generate();
    TLorentzVector *pK      = event_Kpipipi.GetDecay(0);
    TLorentzVector *pPi1    = event_Kpipipi.GetDecay(1);
    TLorentzVector *pPi2    = event_Kpipipi.GetDecay(2);
    TLorentzVector *pPi3    = event_Kpipipi.GetDecay(3);

    
    h0_hypK.SetXYZM(pK->Px(),pK->Py(),pK->Pz(),Mass::K());
    h1_hypK.SetXYZM(pPi1->Px(),pPi1->Py(),pPi1->Pz(),Mass::K());
    mu0_hypMu.SetXYZM(pPi2->Px(),pPi2->Py(),pPi2->Pz(),Mass::Mu());
    mu1_hypMu.SetXYZM(pPi3->Px(),pPi3->Py(),pPi3->Pz(),Mass::Mu());

    h_mKpipipi_asKKmumu ->Fill((h0_hypK+h1_hypK+mu0_hypMu+mu1_hypMu).M(),weight1);
    
    Double_t weight2 = event_Kpimumu.Generate();
    TLorentzVector *pK_2      = event_Kpimumu.GetDecay(0);
    TLorentzVector *pPi_2    = event_Kpimumu.GetDecay(1);
    TLorentzVector *pMu1_2    = event_Kpimumu.GetDecay(2);
    TLorentzVector *pMu2_2    = event_Kpimumu.GetDecay(3);

    h0_hypK.SetXYZM(pK_2->Px(),pK_2->Py(),pK_2->Pz(),Mass::K());
    h1_hypK.SetXYZM(pPi_2->Px(),pPi_2->Py(),pPi_2->Pz(),Mass::K());
    mu0_hypMu.SetXYZM(pMu1_2->Px(),pMu1_2->Py(),pMu1_2->Pz(),Mass::Mu());
    mu1_hypMu.SetXYZM(pMu2_2->Px(),pMu2_2->Py(),pMu2_2->Pz(),Mass::Mu());

    h_mKpimumu_asKKmumu ->Fill((h0_hypK+h1_hypK+mu0_hypMu+mu1_hypMu).M(),weight2);  

    Double_t weight3 = event_KKmumu.Generate();
    TLorentzVector *pK_3      = event_KKmumu.GetDecay(0);
    TLorentzVector *pK2_3      = event_KKmumu.GetDecay(1);
    TLorentzVector *pMu1_3    = event_KKmumu.GetDecay(2);
    TLorentzVector *pMu2_3    = event_KKmumu.GetDecay(3);

    h0_hypK.SetXYZM(pK_3->Px(),pK_3->Py(),pK_3->Pz(),Mass::K());
    h1_hypPi.SetXYZM(pK2_3->Px(),pK2_3->Py(),pK2_3->Pz(),Mass::Pi());
    mu0_hypPi.SetXYZM(pMu1_3->Px(),pMu1_3->Py(),pMu1_3->Pz(),Mass::Pi());
    mu1_hypPi.SetXYZM(pMu2_3->Px(),pMu2_3->Py(),pMu2_3->Pz(),Mass::Pi());

    h_mKKmumu_asKpipipi ->Fill((h0_hypK+h1_hypPi+mu0_hypPi+mu1_hypPi).M(),weight3);  

    Double_t weight4 = event_Kpipipipi0.Generate();
    TLorentzVector *pK_4      = event_Kpipipipi0.GetDecay(0);
    TLorentzVector *pPi_4     = event_Kpipipipi0.GetDecay(1);
    TLorentzVector *pPi1_4    = event_Kpipipipi0.GetDecay(2);
    TLorentzVector *pPi2_4    = event_Kpipipipi0.GetDecay(3);
    TLorentzVector *pPi0_4    = event_Kpipipipi0.GetDecay(4);

    h0_hypK.SetXYZM(pK_4->Px(),pK_4->Py(),pK_4->Pz(),Mass::K());
    h1_hypK.SetXYZM(pPi_4->Px(),pPi_4->Py(),pPi_4->Pz(),Mass::K());
    mu0_hypMu.SetXYZM(pPi1_4->Px(),pPi1_4->Py(),pPi1_4->Pz(),Mass::Mu());
    mu1_hypMu.SetXYZM(pPi2_4->Px(),pPi2_4->Py(),pPi2_4->Pz(),Mass::Mu());

    h_mKpipipipi0_asKKmumu ->Fill((h0_hypK+h1_hypPi+mu0_hypMu+mu1_hypMu).M(),weight4);  

    //TLorentzVector temp1; temp1.SetXYZM(pPi_4->Px(),pPi_4->Py(),pPi_4->Pz(),Mass::Pi());
    //TLorentzVector temp2; temp2.SetXYZM(pPi1_4->Px(),pPi1_4->Py(),pPi1_4->Pz(),Mass::Pi());
    //TLorentzVector temp3; temp3.SetXYZM(pPi2_4->Px(),pPi2_4->Py(),pPi2_4->Pz(),Mass::Pi());
    //TLorentzVector temp4; temp4.SetXYZM(pPi0_4->Px(),pPi0_4->Py(),pPi0_4->Pz(),Mass::Pi0());

    //std::cout<<"mass "<< (h0_hypK+temp1+temp2+temp3+temp4).M()  << std::endl;
        
  }

  canvas->cd(5);
  h_mKpipipi_asKKmumu->Draw();
  canvas->cd(6);
  h_mKpimumu_asKKmumu->Draw();
  canvas->cd(8);
  h_mKKmumu_asKpipipi->Draw();
  canvas->SaveAs(fileName);
  
  TCanvas* canvas2 = new TCanvas("b","b");
  h_mKpipipipi0_asKKmumu->Draw();
  canvas2->SaveAs("test.eps");

}

//make tuples with the generator level candidates. This is needed to get the right normlization on the efficiency calculation

void createGeneratorLevelMCTuple (TString kind) {

  TChain* fChain = new TChain("MCTruthTuple/MCTruthTuple");

  for(int i=0;i<50;++i) {
    fChain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/"+kind+"_withGenLevel/magUp/MC12_Dst"+kind+"_%i.root",i));
    fChain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/"+kind+"_withGenLevel/magDw/MC12_Dst"+kind+"_%i.root",i));
  }    

  double mu0_PX,mu0_PY,mu0_PZ;
  double mu1_PX,mu1_PY,mu1_PZ;

  if(kind=="D2KKmumu" || kind=="D2Kpimumu" || kind=="D2pipimumu"){
    fChain->SetBranchAddress("muminus_TRUEP_X",&mu1_PX);
    fChain->SetBranchAddress("muminus_TRUEP_Y",&mu1_PY);
    fChain->SetBranchAddress("muminus_TRUEP_Z",&mu1_PZ);
    fChain->SetBranchAddress("muplus_TRUEP_X",&mu0_PX);
    fChain->SetBranchAddress("muplus_TRUEP_Y",&mu0_PY);
    fChain->SetBranchAddress("muplus_TRUEP_Z",&mu0_PZ);
  }

  if(kind=="D2KKpipi"){
    fChain->SetBranchAddress("piminus_TRUEP_X",&mu1_PX);
    fChain->SetBranchAddress("piminus_TRUEP_Y",&mu1_PY);
    fChain->SetBranchAddress("piminus_TRUEP_Z",&mu1_PZ);
    fChain->SetBranchAddress("piplus_TRUEP_X",&mu0_PX);
    fChain->SetBranchAddress("piplus_TRUEP_Y",&mu0_PY);
    fChain->SetBranchAddress("piplus_TRUEP_Z",&mu0_PZ);
  }

  if(kind=="D2Kpipipi"){
    fChain->SetBranchAddress("piminus_TRUEP_X",&mu1_PX);
    fChain->SetBranchAddress("piminus_TRUEP_Y",&mu1_PY);
    fChain->SetBranchAddress("piminus_TRUEP_Z",&mu1_PZ);
    fChain->SetBranchAddress("piplus0_TRUEP_X",&mu0_PX);
    fChain->SetBranchAddress("piplus0_TRUEP_Y",&mu0_PY);
    fChain->SetBranchAddress("piplus0_TRUEP_Z",&mu0_PZ);
  }

  TLorentzVector mu0,mu1;  
  double D_DiMuon_Mass; 

  //location of tuples again set to HD location
  TString fName = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/"+kind+"_MCgeneratorLevelTuple.root";
  TFile *newfile = new TFile(fName,"recreate");
  TTree *newtree = new TTree("MCTruthTuple","MCTruthTuple");
  newtree->Branch("D_TRUE_DiMuon_Mass",&D_DiMuon_Mass);
  newtree->Branch("D_DiMuon_Mass",&D_DiMuon_Mass);

  for (Long64_t i=0;i<fChain->GetEntries(); i++) {
    fChain->GetEntry(i);
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass = (mu0+mu1).M();
    newtree->Fill();
  }
  newtree->AutoSave();
  delete newfile;
}

void studyKKMCEfficiency(double BDTCut, double PIDCut, TString fOut) {

  TH1* h_nSelKKmumu = new TH1D ("h_nSelKKmumu","h_nSelKKmumu",50,200,900);
  TH1* h_nGenKKmumu = new TH1D ("h_nGenKKmumu","h_nGenKKmumu",50,200,900);

  TH1* h_nSelKpimumu = new TH1D ("h_nSelKpimumu","h_nSelKpimumu",50,200,1200);
  TH1* h_nGenKpimumu = new TH1D ("h_nGenKpimumu","h_nGenKpimumu",50,200,1200);

  h_nSelKKmumu->Sumw2();
  h_nGenKKmumu->Sumw2();
  h_nSelKpimumu->Sumw2();
  h_nGenKpimumu->Sumw2();

  TChain* c_KKmumuSel = new TChain("BDT_Tree");
  TChain* c_KKmumuGen = new TChain("MCTruthTuple");
  TChain* c_KpimumuSel = new TChain("BDT_Tree");
  TChain* c_KpimumuGen = new TChain("MCTruthTuple");
  
  double D_DiMuon_Mass;
  double BDT;
  double mu0_ProbNNmu, mu1_ProbNNmu;
  double nSelNorm=0;
  double nTotNorm=0;  

  c_KKmumuSel->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root");
  c_KpimumuSel->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");

  c_KKmumuGen->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_MCgeneratorLevelTuple.root");
  c_KpimumuGen->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCgeneratorLevelTuple.root");


  c_KKmumuSel->SetBranchAddress("BDT",&BDT);
  c_KpimumuSel->SetBranchAddress("BDT",&BDT);
  c_KKmumuSel->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  c_KpimumuSel->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  c_KKmumuSel->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
  c_KpimumuSel->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
 
  c_KKmumuGen->SetBranchAddress("D_DiMuon_Mass",&D_DiMuon_Mass); //generator level DiMuon Mass is TRUE quantity
  c_KpimumuGen->SetBranchAddress("D_DiMuon_Mass",&D_DiMuon_Mass);
  c_KKmumuSel->SetBranchAddress("D_TRUE_DiMuon_Mass",&D_DiMuon_Mass); //also here true qunatities!
  c_KpimumuSel->SetBranchAddress("D_TRUE_DiMuon_Mass",&D_DiMuon_Mass);


  for(int i=0;i<c_KKmumuSel->GetEntries();++i) {
    c_KKmumuSel->GetEntry(i);
    if (BDT<BDTCut || mu0_ProbNNmu < PIDCut || mu1_ProbNNmu<PIDCut ) continue;
    h_nSelKKmumu->Fill(D_DiMuon_Mass);
  }
  for(int i=0;i<c_KKmumuGen->GetEntries();++i) {
    c_KKmumuGen->GetEntry(i);
    h_nGenKKmumu->Fill(D_DiMuon_Mass);
  }
  
  for(int i=0;i<c_KpimumuSel->GetEntries();++i) {
    c_KpimumuSel->GetEntry(i);
    if (BDT<BDTCut || mu0_ProbNNmu < PIDCut || mu1_ProbNNmu<PIDCut ) continue;
    h_nSelKpimumu->Fill(D_DiMuon_Mass);
    if(D_DiMuon_Mass>675&&D_DiMuon_Mass<875)nSelNorm+=1;
  }
  
  for(int i=0;i<c_KpimumuGen->GetEntries();++i) {
    c_KpimumuGen->GetEntry(i);
    h_nGenKpimumu->Fill(D_DiMuon_Mass);
    if(D_DiMuon_Mass>675&&D_DiMuon_Mass<875)nTotNorm+=1;
  }

  std::cout<<"nSelNorm/nTotNorm "<<nSelNorm/nTotNorm<<std::endl;

  TEfficiency* effKKmumu = new TEfficiency(*h_nSelKKmumu,*h_nGenKKmumu);  
  TEfficiency* effKpimumu= new TEfficiency(*h_nSelKpimumu,*h_nGenKpimumu);
  /*
  TH1* effKKmumu = (TH1*)h_nSelKKmumu->Clone();
  effKKmumu->Divide(h_nGenKKmumu);
  TH1* effKpimumu = (TH1*)h_nSelKpimumu->Clone();
  effKpimumu->Divide(h_nGenKpimumu);
  */

  gStyle->SetOptStat(0);

  TCanvas* a = new TCanvas("a","a");
  a->Divide(2,3);
  a->cd(1);
  h_nSelKKmumu->Draw();
  a->cd(3);
  h_nGenKKmumu->Draw();
  a->cd(2);
  h_nSelKpimumu->Draw();
  a->cd(4);
  h_nGenKpimumu->Draw();
  a->cd(5);
  effKKmumu->Draw();
  a->cd(6);
  effKpimumu->Draw();
  
  a->SaveAs(fOut);
}

void studypipiMCEfficiency(double BDTCut, double PIDCut, TString fOut) {

  TH1* h_nSelpipimumu = new TH1D ("h_nSelpipimumu","h_nSelpipimumu",50,200,1450);
  TH1* h_nGenpipimumu = new TH1D ("h_nGenpipimumu","h_nGenpipimumu",50,200,1450);

  TH1* h_nSelKpimumu = new TH1D ("h_nSelKpimumu","h_nSelKpimumu",50,200,1200);
  TH1* h_nGenKpimumu = new TH1D ("h_nGenKpimumu","h_nGenKpimumu",50,200,1200);

  h_nSelpipimumu->Sumw2();
  h_nGenpipimumu->Sumw2();
  h_nSelKpimumu->Sumw2();
  h_nGenKpimumu->Sumw2();

  TChain* c_pipimumuSel = new TChain("BDT_Tree");
  TChain* c_pipimumuGen = new TChain("MCTruthTuple");
  TChain* c_KpimumuSel = new TChain("BDT_Tree");
  TChain* c_KpimumuGen = new TChain("MCTruthTuple");
  
  double D_DiMuon_Mass;
  double BDT;
  double mu0_ProbNNmu, mu1_ProbNNmu;
  double nSelNorm=0;
  double nTotNorm=0;  

  c_pipimumuSel->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2pipimumu_BDT.root");
  c_KpimumuSel->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2pipimumuBDT.root");

  c_pipimumuGen->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_MCgeneratorLevelTuple.root");
  c_KpimumuGen->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCgeneratorLevelTuple.root");


  c_pipimumuSel->SetBranchAddress("BDT",&BDT);
  c_KpimumuSel->SetBranchAddress("BDT",&BDT);
  c_pipimumuSel->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  c_KpimumuSel->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  c_pipimumuSel->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
  c_KpimumuSel->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
 
  c_pipimumuGen->SetBranchAddress("D_DiMuon_Mass",&D_DiMuon_Mass); //generator level DiMuon Mass is TRUE quantity
  c_KpimumuGen->SetBranchAddress("D_DiMuon_Mass",&D_DiMuon_Mass);
  c_pipimumuSel->SetBranchAddress("D_TRUE_DiMuon_Mass",&D_DiMuon_Mass); //also here true qunatities!
  c_KpimumuSel->SetBranchAddress("D_TRUE_DiMuon_Mass",&D_DiMuon_Mass);


  for(int i=0;i<c_pipimumuSel->GetEntries();++i) {
    c_pipimumuSel->GetEntry(i);
    if (BDT<BDTCut || mu0_ProbNNmu < PIDCut || mu1_ProbNNmu<PIDCut ) continue;
    h_nSelpipimumu->Fill(D_DiMuon_Mass);
  }
  for(int i=0;i<c_pipimumuGen->GetEntries();++i) {
    c_pipimumuGen->GetEntry(i);
    h_nGenpipimumu->Fill(D_DiMuon_Mass);
  }
  
  for(int i=0;i<c_KpimumuSel->GetEntries();++i) {
    c_KpimumuSel->GetEntry(i);
    if (BDT<BDTCut || mu0_ProbNNmu < PIDCut || mu1_ProbNNmu<PIDCut ) continue;
    h_nSelKpimumu->Fill(D_DiMuon_Mass);
    if(D_DiMuon_Mass>675&&D_DiMuon_Mass<875)nSelNorm+=1;
  }
  
  for(int i=0;i<c_KpimumuGen->GetEntries();++i) {
    c_KpimumuGen->GetEntry(i);
    h_nGenKpimumu->Fill(D_DiMuon_Mass);
    if(D_DiMuon_Mass>675&&D_DiMuon_Mass<875)nTotNorm+=1;
  }

  std::cout<<"nSelNorm/nTotNorm "<<nSelNorm/nTotNorm<<std::endl;

  TEfficiency* effpipimumu = new TEfficiency(*h_nSelpipimumu,*h_nGenpipimumu);  
  TEfficiency* effKpimumu= new TEfficiency(*h_nSelKpimumu,*h_nGenKpimumu);
  /*
  TH1* effpipimumu = (TH1*)h_nSelpipimumu->Clone();
  effpipimumu->Divide(h_nGenpipimumu);
  TH1* effKpimumu = (TH1*)h_nSelKpimumu->Clone();
  effKpimumu->Divide(h_nGenKpimumu);
  */

  gStyle->SetOptStat(0);

  TCanvas* a = new TCanvas("a","a");
  a->Divide(2,3);
  a->cd(1);
  h_nSelpipimumu->Draw();
  a->cd(3);
  h_nGenpipimumu->Draw();
  a->cd(2);
  h_nSelKpimumu->Draw();
  a->cd(4);
  h_nGenKpimumu->Draw();
  a->cd(5);
  effpipimumu->Draw();
  a->cd(6);
  effKpimumu->Draw();
  
  a->SaveAs(fOut);
}





  int nPar=3;
  double mMin=0;
  double mMax=1600;

  Double_t myfunction(Double_t *x, Double_t *p)
  {
    Double_t f = p[0] * relativisticBreitWigner(x[0],p[2],p[1]);
    return f;
  }

//Breit-Wigner function
Double_t mybw(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  return par[0]*arg1*arg2/(arg3 + arg4);
}
    
void define_binning(TString kind) {


  RooRealVar m("m","m(#mu#mu)",0,1600);

  RooRealVar mRho("mRho","mRho",775.26);
  RooRealVar wRho("wRho","wRho",149.1);

  RooRealVar mOmega("mOmega","mOmega",782.65);
  RooRealVar wOmega("wOmega","wOmega",8.49);
 
  RooRealVar mPhi("mPhi","mPhi",1019.455);
  RooRealVar wPhi("wPhi","wPhi",4.26);

  RooRealVar mEta("mEta","mEta",547.862);
  RooRealVar wEta("wEta","wEta",0.00131);


  RooGenericPdf pdfRho("pdfRho","pdfRho"," 14.0/22.0 * @0*@0*@1*@1 /(  ((@2*@2) - (@1*@1))*((@2*@2) - (@1*@1))  +  @2*@2*@2*@2*((@0*@0)/(@1*@1))  )",RooArgSet(wRho,mRho,m));
  RooGenericPdf pdfOmega("pdfOmega","pdfOmega"," 14.0/22.0 * @0*@0*@1*@1 /(  ((@2*@2) - (@1*@1))*((@2*@2) - (@1*@1))  +  @2*@2*@2*@2*((@0*@0)/(@1*@1))  )",RooArgSet(wOmega,mOmega,m));
  RooGenericPdf pdfPhi("pdfPhi","pdfPhi"," 14.0/22.0 * @0*@0*@1*@1 /(  ((@2*@2) - (@1*@1))*((@2*@2) - (@1*@1))  +  @2*@2*@2*@2*((@0*@0)/(@1*@1))  )",RooArgSet(wPhi,mPhi,m));
  RooGenericPdf pdfEta("pdfEta","pdfEta"," 14.0/22.0 * @0*@0*@1*@1 /(  ((@2*@2) - (@1*@1))*((@2*@2) - (@1*@1))  +  @2*@2*@2*@2*((@0*@0)/(@1*@1))  )",RooArgSet(wEta,mEta,m));

  //RooAddPdf sum ("sum","sum",RooArgList(pdfRho,pdfOmega,pdfPhi));


  TCanvas a("a","a");
  gPad-> SetLogy() ;  

RooPlot* xframe = m.frame(Title("Interpreted expression pdf")) ;
  //sum.plotOn(xframe) ;
  

  pdfRho.plotOn(xframe) ;
  pdfOmega.plotOn(xframe) ;
  pdfPhi.plotOn(xframe) ;
  //pdfEta.plotOn(xframe) ;

  xframe->Draw();
  a.Print("test.eps");


}

void Draw_MC_L0Muon_Efficiencies(TString normSelection, TString variable, int nBins, int xLow, int xHigh, TString SelCut) {

  //TString outname="../"+variable+".root" ; //to de defined
  TChain* Chain_MC_sig = new TChain("BDT_Tree");
  //TChain* Chain_MC_norm = new TChain("BDT_Tree");
  
  TString inputfile_MC_sig= "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT_noCuts.root";
  //TString inputfile_MC_Norm= "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_noCuts_D2KKmumuBDT.root";

  double x;
  bool mu0_L0MuonDecision_TOS,mu1_L0MuonDecision_TOS;
  bool mu0_L0MuonDecision_TIS,mu1_L0MuonDecision_TIS;
  bool mu0_L0MuonDecision_Dec,mu1_L0MuonDecision_Dec;
  bool h0_L0HadronDecision_Dec,h0_L0HadronDecision_TIS,h0_L0HadronDecision_TOS;
  bool h1_L0HadronDecision_Dec,h1_L0HadronDecision_TIS,h1_L0HadronDecision_TOS;
  bool Dst_L0Global_TIS;
  bool normSelection_TOS,normSelection_TIS, normSelection_Dec;
  double mu0_PT, mu1_PT;

  Chain_MC_sig->SetBranchAddress(normSelection+"_TOS",&normSelection_TOS);
  Chain_MC_sig->SetBranchAddress(normSelection+"_TIS",&normSelection_TIS);
  Chain_MC_sig->SetBranchAddress(normSelection+"_Dec",&normSelection_Dec);
  Chain_MC_sig->SetBranchAddress("mu1_PT",&mu0_PT);
  Chain_MC_sig->SetBranchAddress("mu0_PT",&mu1_PT);
  
  Chain_MC_sig->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Chain_MC_sig->SetBranchAddress("mu0_L0MuonDecision_TIS",&mu0_L0MuonDecision_TIS);
  Chain_MC_sig->SetBranchAddress("mu0_L0MuonDecision_Dec",&mu0_L0MuonDecision_Dec);
  Chain_MC_sig->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Chain_MC_sig->SetBranchAddress("mu1_L0MuonDecision_TIS",&mu1_L0MuonDecision_TIS);
  Chain_MC_sig->SetBranchAddress("mu1_L0MuonDecision_Dec",&mu1_L0MuonDecision_Dec);
  Chain_MC_sig->SetBranchAddress(variable,&x);


  Chain_MC_sig->AddFile(inputfile_MC_sig);
  //Chain_MC_norm->AddFile(inputfile_MC_Norm);

  TTree* Chain_MC_sig_cut = Chain_MC_sig->CopyTree(SelCut+"&&mu0_MuonNShared==0&&mu1_MuonNShared==0");
  //TTree* Chain_MC_norm_cut= Chain_MC_norm->CopyTree(SelCut);

  //TFile *outfile = new TFile(outname,"recreate");
  TCanvas *c1_MC_sig = new TCanvas("c1_MC_sig","c1_MC_sig",200,10,700,500);
  TCanvas *c1_MC_norm = new TCanvas("c1_MC_norm","c1_MC_norm",200,10,700,500);

  TH1* h1_nMuonTOS = new TH1D("h1_nMuonTOS"+variable,"nTOS as function of "+variable,nBins,xLow,xHigh);
  TH1* h1_nMuonTIS = new TH1D("h1_nMuonTIS"+variable,"nTIS as function of "+variable,nBins,xLow,xHigh);
  TH1* h1_nMuonTISTOS = new TH1D("h1_nMuonTISTOS"+variable,"nTISTOS as function of "+variable,nBins,xLow,xHigh);
  TH1* h1_nMuonDec = new TH1D("h1_nMuonDec"+variable,"nDec as function of "+variable,nBins,xLow,xHigh);
  TH1* h1_nMuonSel = new TH1D("h1_nMuonSel"+variable,"nSel as function of "+variable,nBins,xLow,xHigh);
  
  
  for (int i=0; i<Chain_MC_sig_cut->GetEntries();++i){
    
    Chain_MC_sig_cut->GetEntry(i);

    bool L0TOS=mu1_L0MuonDecision_TOS||mu0_L0MuonDecision_TOS;
    bool L0Dec=mu1_L0MuonDecision_Dec||mu0_L0MuonDecision_Dec;
    bool L0TIS= normSelection_TIS;
    bool L0TISTOS= normSelection_TIS && L0TOS;
    
    std::cout<<i<<"  "<< x  <<"  "<<L0TOS<<"  "<<L0TIS<<"  "<<L0TISTOS<<std::endl;
    
    h1_nMuonSel->Fill(x);
    if(L0TOS) h1_nMuonTOS->Fill(x);
    if(L0TIS) h1_nMuonTIS->Fill(x);
    if(L0Dec) h1_nMuonDec->Fill(x);
    if(L0TISTOS) h1_nMuonTISTOS->Fill(x);
    
  }

  TEfficiency* Eff_TOS_TISTOS = new TEfficiency(*h1_nMuonTISTOS,*h1_nMuonTIS);
  TEfficiency* Eff_TOS_TRUE= new TEfficiency(*h1_nMuonTOS,*h1_nMuonSel);
  //TEfficiency* Eff_TIS_TISTOS= new TEfficiency(*h1_nMuonTISTOS,*h1_nMuonTOS);
  //TEfficiency* Eff_TIS_TRUE= new TEfficiency(*h1_nMuonTIS,*h1_nMuonSel);

  TCanvas *a = new TCanvas("a","a");
  a->Divide(2,3);
  a->cd(1);
  Eff_TOS_TISTOS->Draw();
  a->cd(2);
  Eff_TOS_TRUE->Draw();
  a->cd(3);
  h1_nMuonTISTOS->Draw();
  a->cd(4);
  h1_nMuonTOS->Draw();
  a->cd(5);
  h1_nMuonTIS->Draw();
  a->cd(6);
  h1_nMuonSel->Draw();
  a->Print("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/TISTOSEfficiency/"+normSelection+"_"+variable+"_"+SelCut+"_1.eps");
  
  TCanvas *b = new TCanvas("b","b");
  a->cd(7);
  Eff_TOS_TISTOS->SetLineColor(kRed);
  Eff_TOS_TISTOS->Draw();
  Eff_TOS_TRUE->Draw("SAME");
  b->Print("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/TISTOSEfficiency/"+normSelection+"_"+variable+"_"+SelCut+"_2.eps");

  delete Chain_MC_sig_cut;
  delete Chain_MC_sig;
  delete a;
  delete b;

}




void Draw_MC_L0Muon_Efficiencies_2D() {

  //TString outname="../"+variable+".root" ; //to de defined
  TChain* Chain_MC_sig = new TChain("BDT_Tree");
  //TChain* Chain_MC_norm = new TChain("BDT_Tree");
  
  TString inputfile_MC_sig= "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/noPreselectionCuts/MC_D2Kpimumu_D2KKmumuBDT_noCuts.root";
  //TString inputfile_MC_Norm= "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_noCuts_D2KKmumuBDT.root";

  double D_PT, D_PZ; 
  bool mu0_L0MuonDecision_TOS,mu1_L0MuonDecision_TOS;
  bool mu0_L0MuonDecision_TIS,mu1_L0MuonDecision_TIS;
  bool mu0_L0MuonDecision_Dec,mu1_L0MuonDecision_Dec;
  bool h0_L0HadronDecision_Dec,h0_L0HadronDecision_TIS,h0_L0HadronDecision_TOS;
  bool h1_L0HadronDecision_Dec,h1_L0HadronDecision_TIS,h1_L0HadronDecision_TOS;
  bool Dst_L0Global_TIS;
  double mu0_PT, mu1_PT;
  bool D_L0MuonDecision_TIS,D_L0MuonDecision_TOS,D_L0MuonDecision_Dec;
  bool D_L0HadronDecision_TIS,D_L0HadronDecision_TOS,D_L0HadronDecision_Dec;
  int nTracks;
  bool D_L0Global_TIS;
  double Dst_PT, Dst_PZ;

  Chain_MC_sig->SetBranchAddress("mu1_PT",&mu0_PT);
  Chain_MC_sig->SetBranchAddress("mu0_PT",&mu1_PT);
  Chain_MC_sig->SetBranchAddress("D_PT",&D_PT);
  Chain_MC_sig->SetBranchAddress("D_PZ",&D_PZ);
  Chain_MC_sig->SetBranchAddress("Dst_PT",&Dst_PT);
  Chain_MC_sig->SetBranchAddress("Dst_PZ",&Dst_PZ);
  Chain_MC_sig->SetBranchAddress("nTracks",&nTracks);
  
  Chain_MC_sig->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Chain_MC_sig->SetBranchAddress("mu0_L0MuonDecision_TIS",&mu0_L0MuonDecision_TIS);
  Chain_MC_sig->SetBranchAddress("mu0_L0MuonDecision_Dec",&mu0_L0MuonDecision_Dec);
  Chain_MC_sig->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Chain_MC_sig->SetBranchAddress("mu1_L0MuonDecision_TIS",&mu1_L0MuonDecision_TIS);

  Chain_MC_sig->SetBranchAddress("D_L0MuonDecision_TIS",&D_L0MuonDecision_TIS);
  Chain_MC_sig->SetBranchAddress("D_L0MuonDecision_TOS",&D_L0MuonDecision_TOS);
  Chain_MC_sig->SetBranchAddress("D_L0MuonDecision_Dec",&D_L0MuonDecision_Dec);
  Chain_MC_sig->SetBranchAddress("D_L0HadronDecision_TIS",&D_L0HadronDecision_TIS);
  Chain_MC_sig->SetBranchAddress("D_L0HadronDecision_TOS",&D_L0HadronDecision_TOS);
  Chain_MC_sig->SetBranchAddress("D_L0HadronDecision_Dec",&D_L0HadronDecision_Dec);

  Chain_MC_sig->SetBranchAddress("D_L0Global_TIS",&D_L0Global_TIS);

  Chain_MC_sig->AddFile(inputfile_MC_sig);
  //Chain_MC_norm->AddFile(inputfile_MC_Norm);

  TTree* Chain_MC_sig_cut = Chain_MC_sig->CopyTree("Dst_BKGCAT<11 && D_DiMuon_Mass>675&&D_DiMuon_Mass<875&&mu0_ProbNNghost<0.5&&mu1_ProbNNghost<0.5&&h0_ProbNNghost<0.5&&h1_ProbNNghost<0.5&&Slowpi_ProbNNghost<0.5 && h0_ProbNNk>0.2&&h1_ProbNNpi>0.2 && deltaM>144.5&&deltaM<146.5 && mu0_MuonNShared==0&&mu1_MuonNShared==0 &&BDT>0 && mu0_ProbNNmu>0.5 && mu1_ProbNNmu>0.5");

  TFile* fout = new TFile("/work/mitzel/D2hhmumu/dev/D2Kpimumu/img/TISTOSStudies/muonTIS_3bins_DstPT.root","RECREATE");
  //TFile* fout = new TFile("test.root","RECREATE");
  
  std::vector<TH1D*> histos_sel,histos_TISTOS, histos_TIS,histos_TOS;
  TH1D* temp1;
  TH1D* temp2;
  TH1D* temp3;
  TH1D* temp4;

  /////////////////////////////BINNING///////////////////
  ///////////////////////////////////////////////////////

  int nX=3;
  //Float_t xRanges[4] = {0,4000,6000,25000}; 
  //Float_t xRanges[3] = {0,7000,25000}; 
  //Float_t xRanges[3] = {2000,6000,80e3}; 
  Float_t xRanges[4] = {2e3,5e3,7e3,20e3}; // Dst_PT
  //Float_t xRanges[4] = {0,150,500,800};
   
  //std::vector<double> yRanges {0,50e3,100e3,150e3,2000e3};
  // std::vector<double> yRanges {0,50e3,90e3,250e3};
  //  std::vector<double> yRanges {0,80e3,250e3}; //D_PZ, relatively ok
  //   std::vector<double> yRanges {0,150,800};
   std::vector<double> yRanges {0,400e3};

  //int nX=2;
  //Float_t xRanges[3] = {0,8000,25000};                                                                                                                                                   
  //std::vector<double> yRanges {0,90e3,250e3};     

  
 // std::vector<int> yRanges {0,100,150,200,600}; nTracks
  
  //int nX= xRanges.size() ;
  int nY= yRanges.size() ;
  nY-=1;
  
  std::cout<<"nX= "<< nX << ",nY = "<<nY<<std::endl; 

  for(int i = 0;i<nY;++i) {
    temp1 = new TH1D(TString::Format("nSel_bin_%i",i),TString::Format("nSel_bin_%i",i),nX,xRanges);
    temp2 = new TH1D(TString::Format("nTISTOS_bin_%i",i),TString::Format("nTISTOS_bin_%i",i),nX,xRanges);
    temp3 = new TH1D(TString::Format("nTOS_bin_%i",i),TString::Format("nTOS_bin_%i",i),nX,xRanges);
    temp4 = new TH1D(TString::Format("nTIS_bin_%i",i),TString::Format("nTIS_bin_%i",i),nX,xRanges);
    histos_sel.push_back(temp1);
    histos_TISTOS.push_back(temp2);
    histos_TOS.push_back(temp3);
    histos_TIS.push_back(temp4);
   }  

  TCanvas *c1_MC_sig = new TCanvas("c1_MC_sig","c1_MC_sig",200,10,700,500);
  TCanvas *c1_MC_norm = new TCanvas("c1_MC_norm","c1_MC_norm",200,10,700,500);

  for(int i = 0;i<nY;++i) {

    std::cout<<"processing " <<i+1<< "/" <<nY<<std::endl;
    std::cout<<"in Range " <<yRanges[i] << " - " <<yRanges[i+1]<<std::endl;
    
    std::cout<<"entries " <<Chain_MC_sig_cut->GetEntries()<<std::endl;

    
    for (int j=0; j<Chain_MC_sig_cut->GetEntries();++j){
    
      Chain_MC_sig_cut->GetEntry(j);
      //if(D_PZ>yRanges[i] && D_PZ<yRanges[i+1]) {
      //if(nTracks>yRanges[i] && nTracks<yRanges[i+1]) {
      if(Dst_PZ>yRanges[i] && Dst_PZ<yRanges[i+1]) {

	bool L0TOS=D_L0MuonDecision_TOS;
	//bool L0Dec=D_L0MuonDecision_Dec;
	bool L0TIS= D_L0MuonDecision_TIS;
	//bool L0TIS= D_L0Global_TIS;
	//bool L0TIS= D_L0HadronDecision_TIS;
	
	bool L0TISTOS= L0TIS && L0TOS;
    
	histos_sel[i]->Fill(Dst_PT);
	
	//std::cout<<"j "<<j<<std::endl;
	//std::cout<<"D_PT "<<D_PT<<std::endl;
	//std::cout<<"D PZ "<<D_PZ<<std::endl;
	//std::cout<<"L0TOS "<<L0TOS<<std::endl;
	//std::cout<<"L0TIS "<<L0TIS<<std::endl;
	//std::cout<<"L0TISTOS "<<L0TISTOS<<std::endl;
	
	if(L0TOS) histos_TOS[i]->Fill(Dst_PT);
	if(L0TIS) histos_TIS[i]->Fill(Dst_PT);
	if(L0TISTOS) histos_TISTOS[i]->Fill(Dst_PT);
      }
      else continue; 
    }
  
  }
 
  std::vector<TH1D*> Eff_sel,Eff_TISTOS;  
  std::vector<double> errors,errors2;
  TH1D* Eff_sel_temp;
  TH1D* Eff_TISTOS_temp;
 
  for(int i = 0;i<nY;++i) {

      histos_TISTOS[i]->Write();
      histos_TIS[i]->Write();
      histos_TOS[i]->Write();
      histos_sel[i]->Write();
    
    //calculating the error
    for(int j = 1;j<nX+1;++j) {
      double k = histos_TOS[i]->GetBinContent(j);
      double n = histos_sel[i]->GetBinContent(j);
      double dE = 1/n * TMath::Sqrt(k *(1 - (k/n) ) );
      errors.push_back(dE);
      std::cout<<"true eff errors: " <<" j: "<<j<<" dE: "<<dE<<" n: "<<n<<" k: "<<k<< "k/n: "<<k/n<< std::endl;
    }
    //Dividing histogramms
    Eff_sel_temp = (TH1D*)histos_TOS[i]->Clone();
    Eff_sel_temp->Divide(histos_sel[i]);
    
    //setting the error
    for(int j = 1;j<nX+1;++j) {
      Eff_sel_temp->SetBinError(j, errors[j]);
    }
    Eff_sel.push_back(Eff_sel_temp);
    
    //calculating the errors for tistos
    for(int j = 1;j<nX+1;++j) {
      double k = histos_TISTOS[i]->GetBinContent(j);
      double n = histos_TIS[i]->GetBinContent(j);
      double dE = 1/n * TMath::Sqrt(k *(1 - (k/n) ) );
      errors2.push_back(dE);
      std::cout<<"TISTOS eff errors: " <<" j: "<<j<<" dE: "<<dE<<" n: "<<n<<" k: "<<k<< "k/n: "<<k/n<< std::endl;
    }
    
    Eff_TISTOS_temp = (TH1D*)histos_TISTOS[i]->Clone();
    Eff_TISTOS_temp->Divide(histos_TIS[i]);
    
    //setting the error
    for(int j = 1;j<nX+1;++j) {
      Eff_TISTOS_temp->SetBinError(j, errors2[j]);
    }
  
    Eff_TISTOS.push_back(Eff_TISTOS_temp);
    
  }



  ///////////////////////////////////////
  ///////////// DATA ////////////////////

  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/noPreselectionCuts/D2Kpimumu_D2KKmumuBDT_noCuts.root");
  myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpipipi_D2KKmumuBDT.root");

  TString preselection  = "D_DiMuon_Mass>675&&D_DiMuon_Mass<875&&mu0_ProbNNghost<0.5&&mu1_ProbNNghost<0.5&&h0_ProbNNghost<0.5&&h1_ProbNNghost<0.5&&Slowpi_ProbNNghost<0.5&&mu0_MuonNShared==0&&mu1_MuonNShared==0&&h0_ProbNNk>0.2&&h1_ProbNNpi>0.2";

  TString offlineSelection = "&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.";
  //TString L0TIS_data = "&&D_L0Global_TIS";
  //TString L0TISTOS_data = "&&D_L0Global_TIS&&D_L0MuonDecision_TOS";
  TString L0TIS_data = "&&D_L0MuonDecision_TIS";
  TString L0TISTOS_data = "&&D_L0MuonDecision_TIS&&D_L0MuonDecision_TOS";
  //TString L0TIS_data = "&&D_L0HadronDecision_TIS";
  //TString L0TISTOS_data = "&&D_L0HadronDecision_TIS&&D_L0MuonDecision_TOS";

  double nTIS_data=0;
  double nTISTOS_data=0;
  double nSel_data=0;

  TH1D* temp1_data;
  TH1D* temp2_data;
  TH1D* temp3_data;

  std::vector<TH1D*> histos_TISTOS_data, histos_TIS_data,histos_sel_data;

  for(int i = 0;i<nY;++i) {

    temp1_data = new TH1D(TString::Format("nTISTOS_data_bin_%i",i),TString::Format("nTISTOS_data_bin_%i",i),nX,xRanges);
    temp2_data = new TH1D(TString::Format("nTIS_data_bin_%i",i),TString::Format("nTIS_data_bin_%i",i),nX,xRanges);
    temp3_data = new TH1D(TString::Format("nSel_data_bin_%i",i),TString::Format("nSel_data_bin_%i",i),nX,xRanges);

    histos_TISTOS_data.push_back(temp1_data);
    histos_TIS_data.push_back(temp2_data);
    histos_sel_data.push_back(temp3_data);

  }

  myFitter1D.fit_normalization_MC(preselection+offlineSelection,true,"test1.eps");
  myFitter1D.fit_Kpipipi_misID(preselection+"&&mu1_ProbNNmu>0.&&BDT>0",true,"test2.eps");

  TString targetFolder = "/work/mitzel/D2hhmumu/dev/D2Kpimumu/img/TISTOSStudies/";

  for(int i = 0;i<nY;++i) {

    for(int j = 1;j<nX+1;++j){

      TString totalCut_TIS=preselection+offlineSelection+TString::Format("&&Dst_PT>%f&&Dst_PT<%f",xRanges[j-1],xRanges[j])+TString::Format("&&Dst_PZ>%f&&Dst_PZ<%f",yRanges[i],yRanges[i+1])+L0TIS_data;
      std::cout<<"CUT FOR TIS "<<i<<"  "<<j<<"  "<<totalCut_TIS<<std::endl;
      nTIS_data = myFitter1D.fit_normalization_Data(totalCut_TIS,targetFolder+TString::Format("data_TIS_Dst_PT>%f&&Dst_PT<%f",xRanges[j-1],xRanges[j])+TString::Format("&&Dst_PZ>%f&&Dst_PZ<%f.eps",yRanges[i],yRanges[i+1]));
      histos_TIS_data[i]->SetBinContent(j,nTIS_data);
      TString totalCut_TISTOS=preselection+offlineSelection+TString::Format("&&Dst_PT>%f&&Dst_PT<%f",xRanges[j-1],xRanges[j])+TString::Format("&&Dst_PZ>%f&&Dst_PZ<%f",yRanges[i],yRanges[i+1])+L0TISTOS_data;
      std::cout<<"CUT FOR TISTOS "<<i<<"  "<<j<<"  "<<totalCut_TISTOS<<std::endl;
      nTISTOS_data= myFitter1D.fit_normalization_Data(totalCut_TISTOS,targetFolder+TString::Format("data_TISTOS_Dst_PT>%f&&Dst_PT<%f",xRanges[j-1],xRanges[j])+TString::Format("&&Dst_PZ>%f&&Dst_PZ<%f.eps",yRanges[i],yRanges[i+1]));
      histos_TISTOS_data[i]->SetBinContent(j,nTISTOS_data);
      TString totalCut_sel=preselection+offlineSelection+TString::Format("&&Dst_PT>%f&&Dst_PT<%f",xRanges[j-1],xRanges[j])+TString::Format("&&Dst_PZ>%f&&Dst_PZ<%f",yRanges[i],yRanges[i+1]);
      std::cout<<"CUT FOR SELECTED "<<i<<"  "<<j<<"  "<<totalCut_TISTOS<<std::endl;
      nSel_data= myFitter1D.fit_normalization_Data(totalCut_sel,targetFolder+TString::Format("data_sel_Dst_PT>%f&&Dst_PT<%f",xRanges[j-1],xRanges[j])+TString::Format("&&Dst_PZ>%f&&Dst_PZ<%f.eps",yRanges[i],yRanges[i+1]));
      histos_sel_data[i]->SetBinContent(j,nSel_data);
    }
  }


  fout->cd();
  std::vector<TH1D*> Eff_TISTOS_data;
  std::vector<double> errors_data;
  TH1D* Eff_TISTOS_temp_data;

  for(int i = 0;i<nY;++i) {

      histos_TISTOS_data[i]->Write();
      histos_TIS_data[i]->Write();
      histos_sel_data[i]->Write();
    //calculating the errors for tistos                                                                                  
 for(int j = 1;j<nX+1;++j) {
      double k = histos_TISTOS_data[i]->GetBinContent(j);
      double n = histos_TIS_data[i]->GetBinContent(j);
      double dE = 1/n * TMath::Sqrt(k *(1 - (k/n) ) );
      errors_data.push_back(dE);
      std::cout<<"TISTOS eff errors: " <<" j: "<<j<<" dE: "<<dE<<" n: "<<n<<" k: "<<k<< "k/n: "<<k/n<< std::endl;
    }

    Eff_TISTOS_temp_data = (TH1D*)histos_TISTOS_data[i]->Clone();
    Eff_TISTOS_temp_data->Divide(histos_TIS_data[i]);

    //setting the error                                                                                                                                                                                                                                                                                                                                                             
    for(int j = 1;j<nX+1;++j) {
      Eff_TISTOS_temp_data->SetBinError(j, errors_data[j]);
    }

    Eff_TISTOS_data.push_back(Eff_TISTOS_temp_data);

  }


  //////////////////////////////////////////////////
  ////////// DRAWING ///////////////////////////////

  TCanvas *a = new TCanvas("a","a");
  a->Divide(3,3);
  for(int i = 0;i<nY;++i) {
   a->cd(i+1);
   Eff_TISTOS[i]->GetYaxis()->SetRangeUser(0,1);
   Eff_TISTOS[i]->Draw("E");
   Eff_sel[i]->SetLineColor(kRed);
   Eff_sel[i]->Draw("SAME E");
   Eff_TISTOS_data[i]->SetLineColor(kCyan);
   Eff_TISTOS_data[i]->Draw("SAME E");

   TLegend * leg = new TLegend(0.1,0.7,0.48,0.9);
   leg->AddEntry(Eff_TISTOS[i]);
   leg->AddEntry(Eff_sel[i]);
   leg->AddEntry(Eff_TISTOS_data[i]);
   leg->Draw();

   Eff_TISTOS[i]->Write();
   Eff_sel[i]->Write();
   Eff_TISTOS_data[i]->Write();
  }
  
  //a->Print("muonTIS_6_bins.eps");
  a->Write();


  TCanvas *b = new TCanvas("b","b");
  b->Divide(2,nY);
  for(int i = 1;i<nY+1;++i) {
    b->cd(2*i-1);
    //histos_TIS[i]->GetYaxis()->SetRangeUser(0,1);
    histos_TIS[i-1]->DrawNormalized("");
    histos_sel[i-1]->SetLineColor(kRed);
    histos_sel[i-1]->DrawNormalized("SAME E");
    histos_TISTOS[i-1]->SetLineColor(kCyan);
    histos_TISTOS[i-1]->DrawNormalized("SAME E");
    histos_TOS[i-1]->SetLineColor(kBlack);
    histos_TOS[i-1]->DrawNormalized("SAME E");

   TLegend * leg = new TLegend(0.1,0.7,0.48,0.9);
   leg->AddEntry(histos_TISTOS[i-1]);
   leg->AddEntry(histos_sel[i-1]);
   leg->AddEntry(histos_TIS[i-1]);
   leg->AddEntry(histos_TOS[i-1]);
   leg->Draw();

    b->cd(2*i);
      //histos_TIS_data[i]->SetLineColor(n);
    histos_TIS_data[i-1]->DrawNormalized("SAME E");
    histos_sel_data[i-1]->SetLineColor(kRed);
    histos_sel_data[i-1]->DrawNormalized("SAME E");
    histos_TISTOS_data[i-1]->SetLineColor(kBlack);
    histos_TISTOS_data[i-1]->DrawNormalized("SAME E");

    TLegend * leg2 = new TLegend(0.1,0.7,0.48,0.9);
    leg2->AddEntry(histos_TISTOS_data[i-1]);
    leg2->AddEntry(histos_sel_data[i-1]);
    leg2->AddEntry(histos_TIS_data[i-1]);
    leg2->Draw();

  }
  
  b->Write();

  double mean_TISTOS=0;
  double mean_true=0;
  double nSel=0;
  double nTOS=0;
  double nTISTOS=0;
  double nTIS=0;

  double mean_data=0;
  double mean_data_alternative=0;
  double nSel_data2=0;
  double nTIS_data2=0;
  double factor_data=1;
  double factor_data_alternative=1;

  double dY=1;
  double dX=1;
  double factor=1; 


  for(int i = 0;i<nY;++i) {
    //dY = yRanges[i+1] - yRanges[i];
    std::cout<<"dY "<<dY<<std::endl;

    for(int j = 1;j<nX+1;++j) {
      //dX = xRanges[j] - xRanges[j-1];
      factor = histos_sel[i]->GetBinContent(j); //spectrum of MC selected events
      std::cout<<"dX "<<dX<<std::endl;
      nSel+= histos_sel[i]->GetBinContent(j);
      nTOS+= histos_TOS[i]->GetBinContent(j);
      nTISTOS+=histos_TISTOS[i]->GetBinContent(j);
      nTIS+=histos_TIS[i]->GetBinContent(j);
      
      ///DATA
      factor_data=histos_sel_data[i]->GetBinContent(j);
      factor_data_alternative = histos_TIS_data[i]->GetBinContent(j); //take spectrum of TIS data sample
      mean_data+=Eff_TISTOS_data[i]->GetBinContent(j) * factor_data;
      mean_data_alternative+=Eff_TISTOS_data[i]->GetBinContent(j) * factor_data_alternative;
      nSel_data2+=histos_sel_data[i]->GetBinContent(j);
      nTIS_data2+=histos_TIS_data[i]->GetBinContent(j);

      mean_true+= Eff_sel[i]->GetBinContent(j) *factor_data; //this one is weighted by the PT spectrum from data 
      mean_TISTOS += Eff_TISTOS[i]->GetBinContent(j) * dY *dX* factor_data; //also this one 

      std::cout<<"mean TISTOS"<<mean_TISTOS<<"( "<<Eff_TISTOS[i]->GetBinContent(j) <<") * "<< factor << std::endl;
      std::cout<<"mean true"<<mean_true<<"( "<<Eff_sel[i]->GetBinContent(j) <<") * "<< factor << std::endl;
      std::cout<<"mean data"<<mean_data<<"( "<<Eff_TISTOS_data[i]->GetBinContent(j) <<") * "<< factor_data << std::endl;
      std::cout<<"mean data alternative"<<mean_data_alternative<<"( "<<Eff_TISTOS_data[i]->GetBinContent(j) <<") * "<< factor_data_alternative << std::endl;

    }  
  }
  //mean /= (xRanges[nX] - xRanges[0]);
  //mean /= (yRanges[nY] - yRanges[0]);
  mean_TISTOS/= nSel_data2;
  mean_true/=nSel_data2; // take MC eff weighted by data PT/PZ distrubution
  mean_data_alternative/=nTIS_data2;
  mean_data/=nSel_data2;
  

  std::cout << "mean TISTOS weighted" << mean_TISTOS<< "(true unweighted"<<nTOS/nSel<<", nTOS "<<nTOS <<",nSel "<<nSel   <<", nTIS "<<nTIS <<",nTISTOS "<<nTISTOS << ")"<< std::endl;
  std::cout << "difference mean TISTOS weighted and unweighted true eff " << mean_TISTOS - nTOS/nSel<< std::endl;
  std::cout << "true mean MC reweighted " <<  mean_true <<  std::endl;
  std::cout << "mean in data" << mean_data<<" and weighted to distibution in TIS "<<mean_data_alternative<< std::endl;
  std::cout << "diff weighted mean data and weighted MC true" << mean_data - mean_true <<" weighed MC true and MC weighted TISTOS "<<mean_true - mean_TISTOS<< std::endl;


  //  std::cout <<"range X "<< xRanges[nX] - xRanges[0] <<" range Y  "<< yRanges[nY] - yRanges[0] << std::endl;
  
fout->Write();  

}

void Draw_MC_L0Muon_Efficiencies_forVariable(TString normSelection,int nBins,TString SelCut) {
  //Draw_MC_L0Muon_Efficiencies(normSelection,"D_PT",nBins,1500,20000,SelCut);
  // Draw_MC_L0Muon_Efficiencies(normSelection,"Dst_PT",nBins,1500,20000,SelCut);
  //Draw_MC_L0Muon_Efficiencies(normSelection,"D_DiMuon_Mass",nBins,200,1200,SelCut);
  // Draw_MC_L0Muon_Efficiencies(normSelection,"mumu_P",nBins,0,200e3,SelCut);
  //Draw_MC_L0Muon_Efficiencies(normSelection,"mumu_PT",nBins,0,10000,SelCut);
  //Draw_MC_L0Muon_Efficiencies(normSelection,"max_muPT",nBins,0,6000,SelCut);
   Draw_MC_L0Muon_Efficiencies(normSelection,"mu0_PT",nBins,0,4000,SelCut);


   //   Draw_MC_L0Muon_Efficiencies("D_L0Global","nTracks",nBins,0,600,SelCut);
}


void Draw_MC_TriggerEfficiencies( TString variable, int nBins, int xLow, int xHigh, TString SelCut) {

  TH1::SetDefaultSumw2();

  //TString outname="../"+variable+".root" ; //to de defined
  TChain* Chain_MC_sig = new TChain("BDT_Tree");
  //TChain* Chain_MC_norm = new TChain("BDT_Tree");
  
  TString inputfile_MC_sig= "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT_noCuts.root";
  //  TString inputfile_MC_Norm= "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_noCuts_D2KKmumuBDT.root";

  double x;
  bool mu0_L0MuonDecision_TOS,mu1_L0MuonDecision_TOS;
  bool mu0_L0MuonDecision_TIS,mu1_L0MuonDecision_TIS;
  bool mu0_L0MuonDecision_Dec,mu1_L0MuonDecision_Dec;
  
  bool mu0_Hlt1TrackMuonDecision_TOS,mu1_Hlt1TrackMuonDecision_TOS,D_Hlt1TrackAllL0Decision_TOS;
  bool D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS,D_Hlt2CharmSemilepD02KKMuMuDecision_TOS;

  bool normSelection_TOS,normSelection_TIS, normSelection_Dec;

  
  Chain_MC_sig->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Chain_MC_sig->SetBranchAddress("mu0_L0MuonDecision_TIS",&mu0_L0MuonDecision_TIS);
  Chain_MC_sig->SetBranchAddress("mu0_L0MuonDecision_Dec",&mu0_L0MuonDecision_Dec);
  Chain_MC_sig->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Chain_MC_sig->SetBranchAddress("mu1_L0MuonDecision_TIS",&mu1_L0MuonDecision_TIS);
  Chain_MC_sig->SetBranchAddress("mu1_L0MuonDecision_Dec",&mu1_L0MuonDecision_Dec);
  Chain_MC_sig->SetBranchAddress(variable,&x);
  Chain_MC_sig->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  Chain_MC_sig->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  Chain_MC_sig->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);
  Chain_MC_sig->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS",&D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
  Chain_MC_sig->SetBranchAddress("D_Hlt2CharmSemilepD02KKMuMuDecision_TOS",&D_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
  


  Chain_MC_sig->AddFile(inputfile_MC_sig);
  //Chain_MC_norm->AddFile(inputfile_MC_Norm);

  TTree* Chain_MC_sig_cut = Chain_MC_sig->CopyTree(SelCut+"&&mu0_MuonNShared==0&&mu1_MuonNShared==0");
  //TTree* Chain_MC_norm_cut= Chain_MC_norm->CopyTree(SelCut);

  //TFile *outfile = new TFile(outname,"recreate");
  TCanvas *c1_MC_sig = new TCanvas("c1_MC_sig","c1_MC_sig",200,10,700,500);
  TCanvas *c1_MC_norm = new TCanvas("c1_MC_norm","c1_MC_norm",200,10,700,500);

  TH1* h1_L0TOS = new TH1D("h1_L0TOS"+variable,"nTOS L0 as function of "+variable,nBins,xLow,xHigh);
  TH1* h1_HLT1TOS = new TH1D("h1_HLT1TOS"+variable,"nTOS Hlt1 as function of "+variable,nBins,xLow,xHigh);
  TH1* h1_HLT2TOS = new TH1D("h1_HLT2TOS"+variable,"nTOS Hlt2 as function of "+variable,nBins,xLow,xHigh);
  TH1* h1_Sel = new TH1D("h1_Sel"+variable,"nSel as function of "+variable,nBins,xLow,xHigh);
  
  
  for (int i=0; i<Chain_MC_sig_cut->GetEntries();++i){
    
    Chain_MC_sig_cut->GetEntry(i);
    bool L0TOS=mu1_L0MuonDecision_TOS||mu0_L0MuonDecision_TOS;
    bool HLT1TOS=L0TOS && (mu0_Hlt1TrackMuonDecision_TOS || mu1_Hlt1TrackMuonDecision_TOS ||  D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
    bool HLT2TOS=L0TOS&& HLT1TOS && D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;

            
    h1_Sel->Fill(x);
    if(L0TOS) h1_L0TOS->Fill(x);
    if(HLT1TOS) h1_HLT1TOS->Fill(x);
    if(HLT2TOS) h1_HLT2TOS->Fill(x);
    
  }

  ///TEfficiency* Eff_L0TOS= new TEfficiency(*h1_L0TOS,*h1_Sel);
  //TEfficiency* Eff_HLT1TOS= new TEfficiency(*h1_HLT1TOS,*h1_L0TOS);
  //TEfficiency* Eff_HLT2TOS= new TEfficiency(*h1_HLT2TOS,*h1_HLT1TOS);

  TH1* Eff_L0TOS = (TH1D*)h1_L0TOS->Clone();
  Eff_L0TOS->Divide(h1_Sel);
  TH1* Eff_HLT1TOS = (TH1D*) h1_HLT1TOS->Clone();
  Eff_HLT1TOS->Divide(h1_L0TOS);
  TH1* Eff_HLT2TOS =  (TH1D*)h1_HLT2TOS->Clone();
  Eff_HLT2TOS->Divide(h1_HLT1TOS);

  TCanvas *a = new TCanvas("a","a");
  a->Divide(2,2);
  a->cd(1);
  Eff_L0TOS->Draw();
  a->cd(2);
  Eff_HLT1TOS->Draw();
  a->cd(3);
  Eff_HLT2TOS->Draw();
  a->cd(4);
  Eff_L0TOS->SetLineColor(kRed);
  Eff_HLT1TOS->SetLineColor(kBlue);
  Eff_HLT2TOS->SetLineColor(kCyan);
  Eff_L0TOS->GetYaxis()->SetRangeUser(0,1);
  Eff_L0TOS->Draw();
  Eff_HLT1TOS->Draw("SAME");
  Eff_HLT2TOS->Draw("SAME");
  a->Print("test2.eps");
  //a->Print("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/TISTOSEfficiency/"+normSelection+"_"+variable+"_"+SelCut+"_1.eps");
  

  delete Chain_MC_sig_cut;
  delete Chain_MC_sig;
  delete a;

}

void newCutOptimization(TString kind, TString polarity="", TString q2Range="") {

  
  TFile *fout = new TFile("/work/mitzel/D2hhmumu/dev/img/newSelectionOptimisation/newFeb17/"+kind+"_PID_BDT_cut_scan_1D_"+q2Range+polarity+".root","recreate");
  TString path="/work/mitzel/D2hhmumu/dev/img/newSelectionOptimisation/newFeb17/";

  /*
  TH2* h2_misIDBkg = new TH2D("h2_misIDBkg","h2_misIDBkg",15,-0.5,1,10,0,1);
  TH2* h2_misIDBkg_alt = new TH2D("h2_misIDBkg_alt","h2_misIDBkg_alt",15,-0.5,1,10,0,1);
  TH2* h2_combBkg = new TH2D("h2_combBkg","h2_combBkg",15,-0.5,1,10,0,1);
  TH2* h2_sigEff = new TH2D("h2_sigEff","h2_sigEff",15,-0.5,1,10,0,1);
  TH2* h2_relEff = new TH2D("h2_releff","h2_relEff",15,-0.5,1,10,0,1);
  TH2* h2= new TH2D("h2_cutEfficiency","h2_cutEfficiency",15,-0.5,1.0,10,0,1);
  */


  double PIDcuts[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.};
  //double BDTcuts[15]={-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  //double BDTcuts[21]={-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  double BDTcuts[11]={-1,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1};

  int nBinsX=10;
  int nBinsY=10;

  //int nBinsX = sizeof(BDTcuts)/sizeof(BDTcuts[0]);
  //int nBinsY = sizeof(PIDcuts)/sizeof(PIDcuts[0]); 
 
  TH2* h2_misIDBkg = new TH2D("h2_misIDBkg","h2_misIDBkg",nBinsX,BDTcuts[0],BDTcuts[nBinsX],nBinsY,PIDcuts[0],PIDcuts[nBinsY]);
  TH2* h2_misIDBkg_alt = new TH2D("h2_misIDBkg_alt","h2_misIDBkg_alt",nBinsX,BDTcuts[nBinsX],BDTcuts[nBinsX-1],nBinsY,PIDcuts[0],PIDcuts[nBinsY]);
  TH2* h2_combBkg = new TH2D("h2_combBkg","h2_combBkg",nBinsX,BDTcuts[0],BDTcuts[nBinsX],nBinsY,PIDcuts[0],PIDcuts[nBinsY]);
  TH2* h2_sigEff = new TH2D("h2_sigEff","h2_sigEff",nBinsX,BDTcuts[0],BDTcuts[nBinsX],nBinsY,PIDcuts[0],PIDcuts[nBinsY]);
  TH2* h2_relEff = new TH2D("h2_releff","h2_relEff",nBinsX,BDTcuts[0],BDTcuts[nBinsX],nBinsY,PIDcuts[0],PIDcuts[nBinsY]);
  TH2* h2= new TH2D("h2_cutEfficiency","h2_cutEfficiency",nBinsX,BDTcuts[0],BDTcuts[nBinsX],nBinsY,PIDcuts[0],PIDcuts[nBinsY]);

  double misIDBkg,misIDBKG_scaled;
  D2hhmumuFitter1D myFitter;

  myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_"+kind+"_BDT.root");
  myFitter.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+kind+"_BDT.root");
  myFitter.fit_MC("deltaM>144.5&&deltaM<146.5&&"+polarity+"&&"+q2Range,true,"");
  if(kind=="D2KKmumu") myFitter.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_"+kind+"BDT.root");
  else myFitter.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_"+kind+"BDT.root");
  

  TString cut_singleID;
  TString cut_doubleID;

  //take a single shappe for the misID
  cut_singleID="BDT>.4&&mu1_ProbNNmu>.3";
  myFitter.fit_HHpipi_misID(cut_singleID+"&&"+q2Range,true,"test.eps");

  //MC file to get efficiency
  TString fIn1="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_"+kind+"_BDT.root";

  //  TString fIn2=;

    for(int i=0;i<nBinsX;++i){ //BDT
      for(int j=0;j<nBinsY;++j){ //was 9 before //PID

      //for(int i=7;i<8;++i){ //BDT                                                                                                                         
      //for(int j=4;j<7;++j){ //was 9 before //PID  

      cut_doubleID = TString::Format("BDT>%f&&mu0_ProbNNmu>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j],PIDcuts[j]);
      cut_singleID = TString::Format("BDT>%f&&mu1_ProbNNmu>%f",BDTcuts[i],PIDcuts[j]);
      //cut_singleID = TString::Format("BDT>%f",BDTcuts[i]);

      cout<<i<<"  "<<j<<" "<<cut_doubleID<<"  "<<cut_singleID<<endl;
      
      std::pair<double,double> eff = getMCEfficiency(fIn1,"BDT_Tree",cut_doubleID,polarity+"&&"+q2Range);
      double signalEff = eff.first;
      double dSignalEff = eff.second;
      
      //std::pair<std::pair<double,double>,std::pair<double,double> >bkg = myFitter.getBkgFromBlindedFit(cut_doubleID+"&&"+polarity,q2Range,path+"Fits/"+kind+"/"+cut_doubleID+"_"+q2Range+"_"+polarity);

      std::pair<double,double> bkg = myFitter.getBkgFromBlindedFit(cut_doubleID+"&&"+polarity,q2Range,path+"Fits/"+kind+"/"+cut_doubleID+"_"+q2Range+"_"+polarity);

      //int roundedCombBkg = (int)bkg.second;
      //int roundedMisID = (int)bkg.first;
      
      //double combBkg = (double)roundedCombBkg;
      //double misIDBkg = (double)roundedMisID;
      
      /*
      double combBkg=bkg.second.first;
      double misIDBkg=bkg.first.first;
      
      double dcombBkg=bkg.second.second;
      double dmisIDBkg=bkg.first.second;
      
      double nBkg = combBkg+ misIDBkg;
      double dNBkg = TMath::Sqrt(dcombBkg*dcombBkg + dmisIDBkg*dmisIDBkg);
      */
      
      double nBkg= bkg.first;
      double dNBkg= bkg.second;

      //h2_misIDBkg->SetBinContent(i+1,j+1,misIDBkg); 
      //h2_combBkg->SetBinContent(i+1,j+1,combBkg);

      h2_sigEff->SetBinContent(i+1,j+1,signalEff);     

      double FOM=signalEff/( 2.5 + TMath::Sqrt( nBkg ) );
      double dFOM = TMath::Sqrt( TMath::Power(dSignalEff/( 2.5 + TMath::Sqrt( nBkg )),2) +  TMath::Power(signalEff*dNBkg/( (2.5 + TMath::Sqrt(nBkg))*(2.5 + TMath::Sqrt(nBkg))*2*TMath::Sqrt(nBkg)),2) ); 
      double dFOM2 = getDFOMwithToys(signalEff,dSignalEff,nBkg,dNBkg);
      
      double FOM2 = signalEff/( 2.5 + TMath::Sqrt( nBkg+dNBkg ) );
      double dFOM3 = TMath::Sqrt( TMath::Power(dSignalEff/( 2.5 + TMath::Sqrt( nBkg )),2) +  TMath::Power(FOM2-FOM,2) ); 
      
      std::cout<<"analytical solution "<<dFOM<<endl;
      std::cout<<"toy solution "<<dFOM2<<endl;
      std::cout<<"alternative solution "<<dFOM3<<endl;

      h2->SetBinContent(i+1,j+1,FOM);
      if(nBkg<1) 
	h2->SetBinError(i+1,j+1,dFOM3);
      else h2->SetBinError(i+1,j+1,dFOM);

      //cout<<"combBkg "<<combBkg<<" "<<misIDBkg<<" "<<FOM<<"  "<<signalEff<<endl;
      //cout<<"nBkg "<<nBkg<<"  "<<"dNBkg  "<<dNBkg<<"  "<<"dFOM  "<<dFOM<<std::endl;  
      
      }

   }

  
   fout->cd();

   TCanvas* c = new TCanvas("canvas","canvas"); 
 
  c->Divide(2,2);
			  
  Float_t newMargin1 = 0.13;
  Float_t newMargin2 = 0.15;


  c->SetGrid();
  c->SetTicks();
  
  gStyle->SetPalette( 1, 0 );

  gStyle->SetPaintTextFormat( "3.3f" );
  gStyle->SetOptStat(0);

  h2->SetMarkerColor( 0 );
  c->cd(1);
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

  c->Write();
  c->Print(path+kind+"_"+q2Range+"_"+polarity+".eps");
  
  h2_misIDBkg->Write();
  h2_misIDBkg_alt->Write();
  h2_combBkg->Write();
  h2_sigEff->Write();
  h2_relEff->Write();
  h2->Write();
  
  fout->Write();

}

double getDFOMwithToys(double signalEff, double dSignalEff, double nBkg,double dNBkg) {

  double FOM0= signalEff/( 2.5 + TMath::Sqrt( nBkg ) );

  double xlow = FOM0-FOM0;
  double xhigh = FOM0+FOM0;

  TH1* h1= new TH1D("temp","temp",100,xlow,xhigh);
  TRandom3 *generator = new TRandom3(0);
  double FOM, BkgVar, EffVar;
    
  int tries = 0;

  while(tries<4000) {
    

    EffVar=generator->Gaus(signalEff, dSignalEff);
    //BkgVar=TMath::Abs(generator->Gaus(nBkg,dNBkg));
    BkgVar=generator->Gaus(nBkg,dNBkg);
    //BkgVar=generator->PoissonD(dNBkg+0.51);

    FOM = EffVar/( 2.5 + TMath::Sqrt( BkgVar ) );
    //cout<<BkgVar<<"  "<<EffVar<<"  "<<TMath::Sqrt( BkgVar )<<"  "<< 2.5 <<"  "<<FOM<<endl;
    h1->Fill(FOM);
    
    tries++;
  }
  
  double FOM2 = EffVar/( 2.5 + TMath::Sqrt( nBkg+dNBkg ) );
  cout<<FOM0<<endl;
  cout<<FOM2<<endl;
  cout<<FOM0-FOM2<<endl;



  h1->Draw();
  return h1->GetRMS();

}
 
