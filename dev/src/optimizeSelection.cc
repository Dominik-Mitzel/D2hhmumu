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
#include "TEfficiency.h"

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
  
  bool isMC= false;
  D2KKmumuReader* KK_Reader = new D2KKmumuReader(Tree_D2KKmumu);
  KK_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_PreselectedSubsample.root",100);
  //KK_Reader->fillHistograms("../rootFiles/test.root",isMC);                                                                          

}

void D2pipimumuData(){

  TChain* Tree_D2pipimumu = new TChain("DstD2PiPiMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2pipimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));
    Tree_D2pipimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }

  D2pipimumuReader* pipi_Reader = new D2pipimumuReader(Tree_D2pipimumu);
  pipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_PreselectedSubsample.root",100);
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
  // Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/highStat_magDw/MC12_DstD2KKmumu_magDw.root");                                 
  //Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/highStat_magUp/MC12_DstD2KKmumu_magUp.root");  
  //Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/magDw/MC12_DstD2KKmumu_magDw.root");//"old one ""
  //Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu/magUp/MC12_DstD2KKmumu_magUp.root");  

  for (int i=0; i<100; ++i) {
    Tree_MC_D2KKmumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_withGenLevel/magDw/MC12_DstD2KKmumu_%i.root",i));
    Tree_MC_D2KKmumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_withGenLevel/magUp/MC12_DstD2KKmumu_%i.root",i));
  }

  //Mai June TCK                                                                                                                                                                
  //Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_MaiJune/magDw/MC12_DstD2KKmumu_MaiJune_magDw.root");
  //Tree_MC_D2KKmumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_MaiJune/magUp/MC12_DstD2KKmumu_MaiJune_magUp.root");

  D2KKmumuReader* KK_MC_Reader = new D2KKmumuReader(Tree_MC_D2KKmumu);
  KK_MC_Reader->InitMC();
  KK_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_highStat_MCtrainingSample.root");
  //KK_MC_Reader->studyTriggerEfficiency();                                                                                                            

  createGeneratorLevelMCTuple("D2KKmumu");

}


void D2pipimumuMC(){

  TChain* Tree_MC_D2pipimumu = new TChain("MC12_DstD2PiPiMuMu/DecayTree");
  Tree_MC_D2pipimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu_withGenLevel/magDw/MC12_DstD2pipimumu_magDw.root");
  Tree_MC_D2pipimumu->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu_withGenLevel/magUp/MC12_DstD2pipimumu_magUp.root");

  D2pipimumuReader* pipi_MC_Reader = new D2pipimumuReader(Tree_MC_D2pipimumu);
  pipi_MC_Reader->InitMC();
  //pipi_MC_Reader->fillHistograms("../rootfiles/MC2012_Kinematical_Distrubutions_D2KK.root",true);                                                                                        
  pipi_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_MCtrainingSample.root");
  //pipi_MC_Reader->studyTriggerEfficiency();                                           
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
  //KK_MC_Reader->fillHistograms("../rootfiles/MC2012_Kinematical_Distrubutions_D2KK.root",true);                                                   
  Kpi_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_highStat_MCtrainingSample.root");
  //KK_MC_Reader->studyTriggerEfficiency();                                                                                                                                       
  createGeneratorLevelMCTuple("D2Kpimumu");
}


void D2KpipipiData(){

  TChain* Tree_D2Kpipipi = new TChain("DstD2KPiPiPi/DecayTree");

  for (int i=0; i<1500; ++i) {
    Tree_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magUp/2012Data_D2hhhh_%i.root",i));
    Tree_D2Kpipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magDw/2012Data_D2hhhh_%i.root",i));
  }


  D2KpipipiReader* Kpipipi_Reader = new D2KpipipiReader(Tree_D2Kpipipi);
  Kpipipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_PreselectedSubsample.root",10);
}


void D2KKpipiData(){

  TChain* Tree_D2KKpipi = new TChain("DstD2KKPiPi/DecayTree");

  for (int i=0; i<1500; ++i) {
    Tree_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magUp/2012Data_D2hhhh_%i.root",i));
    Tree_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magDw/2012Data_D2hhhh_%i.root",i));
  }

  D2KKpipiReader* KKpipi_Reader = new D2KKpipiReader(Tree_D2KKpipi);
  KKpipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_PreselectedSubsample.root",10);
}

void D2pipipipiData(){

  TChain* Tree_D2pipipipi = new TChain("DstD2PiPiPiPi/DecayTree");

  for (int i=0; i<1500; ++i) {
    Tree_D2pipipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magUp/2012Data_D2hhhh_%i.root",i));
    Tree_D2pipipipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/newD2hhhhPIDLine/magDw/2012Data_D2hhhh_%i.root",i));
  }

  D2pipipipiReader* pipipipi_Reader = new D2pipipipiReader(Tree_D2pipipipi);
  pipipipi_Reader->createSubsample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipipipi_PIDline_PreselectedSubsample.root",10);
}


void D2KKpipiMC(){


  TChain* Tree_MC_D2KKpipi = new TChain("MC12_DstD2KKpipi/DecayTree");

  for (int i=0; i<100; ++i) {
    Tree_MC_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKpipi_withGenLevel/magDw/MC12_DstD2KKpipi_%i.root",i));
    Tree_MC_D2KKpipi->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/MC/D2KKpipi_withGenLevel/magUp/MC12_DstD2KKpipi_%i.root",i));
  }
 

  D2KKpipiReader* KK_MC_Reader = new D2KKpipiReader(Tree_MC_D2KKpipi);
  KK_MC_Reader->InitMC();
  KK_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_MCtrainingSample.root");
                                                                                                                                          
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
  Kpi_MC_Reader->createMCtrainingSample("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_MCtrainingSample.root");                            
  
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
  c_KKmumuSel->SetBranchAddress("D_TRUE_DiMuon_Mass",&D_DiMuon_Mass);
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

  TH1* effKKmumu = (TH1*)h_nSelKKmumu->Clone();
  effKKmumu->Divide(h_nGenKKmumu);
  TH1* effKpimumu = (TH1*)h_nSelKpimumu->Clone();
  effKpimumu->Divide(h_nGenKpimumu);
  
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

    
