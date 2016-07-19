#include "D2hhmumuFitter_Applications.h"
#include "D2hhmumuFitter1D.h"
#include "StandardHypoTestInvDemo.C"


D2hhmumuFitter_Applications::D2hhmumuFitter_Applications(TString m_kind,TString m_year){

  kind= m_kind;
  year= m_year;
  targetFolder="/work/mitzel/D2hhmumu/dev/D2KKmumu/";
  setDefaultPathToData(kind);  
  setQ2Ranges(kind);
  q2RangeNormalizationMode="D_DiMuon_Mass>675&&D_DiMuon_Mass<875";
 }

D2hhmumuFitter_Applications::~D2hhmumuFitter_Applications(){};

void D2hhmumuFitter_Applications::setDefaultPathToData(TString m_kind){

  if(m_kind!="D2KKmumu" && m_kind!="D2pipimumu") std::cout<<" D2hhmumuFitter_Applications: Error decay mode not implemeted"<<std::endl;
  if(m_kind== "D2KKmumu"){
    pathToSignalMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root";
    pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root";
    pathToSidebandData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT.root";
    pathToKpipipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root";//"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root";
    pathToKKpipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root";//"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root"; 
    pathToKpimumuMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root";
    pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root";
    pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root";
  }

  if(m_kind== "D2pipimumu"){
    pathToSignalMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2pipimumu_BDT.root";
    pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2pipimumuBDT.root";
    pathToSidebandData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2pipimumu_BDT.root";
    pathToKpipipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2pipimumuBDT.root";
    pathToKKpipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipipipi_PIDline_D2pipimumuBDT.root";
    pathToKpimumuMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2pipimumuBDT.root";
    pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_BDT.root";
    pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2pipimumuBDT.root";
  }


}


void D2hhmumuFitter_Applications::setQ2Ranges(TString m_kind){
 
  if(m_kind!= ("D2KKmumu" || "D2pipimumu")) std::cout<<" D2hhmumuFitter_Applications: Error decay mode not implemeted"<<std::endl;
  
  if(m_kind== "D2KKmumu"){
    q2Ranges.push_back("D_DiMuon_Mass<525");
    q2Ranges.push_back("D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
    q2Ranges.push_back("D_DiMuon_Mass>565");
  }   

  if(m_kind== "D2pipimumu"){
    q2Ranges.push_back("D_DiMuon_Mass<525");
    q2Ranges.push_back("D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
    q2Ranges.push_back("D_DiMuon_Mass>565&&D_DiMuon_Mass<950");
    q2Ranges.push_back("D_DiMuon_Mass>950&&D_DiMuon_Mass<1100");
    q2Ranges.push_back("D_DiMuon_Mass>1100");
  }   


  std::cout<<"q^2 intervals: "<<std::endl;
  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) { 
    std::cout<<*it<<" , ";
  }
  std::cout<<" "<<std::endl;
}

void D2hhmumuFitter_Applications::saveModelConfig(TString dataCut,TString misIDCut,TString q2Range){

  D2hhmumuFitter1D myFitter1D;                                                                                                                                                         
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root");

  TString fOut=targetFolder+"ModelConfigs/"+ dataCut + "_" + q2Range + ".root"; 
  //  TString totalDataCut=dataCut+"&&"+q2Range;
  //TString totalNormalizationCut=dataCut+"&&"+q2RangeNormalizationMode; //
  //TString totalMisIDCut=misIDCut+"&&"+q2Range;

  myFitter1D.fillModelConfig(dataCut,misIDCut,q2Range,fOut); 

}

void D2hhmumuFitter_Applications::saveAllModelConfig(TString dataCut,TString misIDCut){

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

    saveModelConfig(dataCut,misIDCut,*it);
  }
}
 
void D2hhmumuFitter_Applications::compare_1D_and_2D_fit(TString dataCut,TString normalizationCut,TString misIDCut){

  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.fit_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/comparison1D_2D_fit/MC"+dataCut+".eps");
  myFitter1D.fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/comparison1D_2D_fit/misID"+misIDCut+".eps");
  myFitter1D.fit_normalization_Data(normalizationCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/comparison1D_2D_fit/norm"+normalizationCut+".eps");

  D2hhmumuFitter myFitter; 
  myFitter.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter.fit_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/comparison1D_2D_fit/2D_MC"+dataCut+".eps");
  myFitter.fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/comparison1D_2D_fit/2D_misID"+misIDCut+".eps");
  myFitter.fit_normalization_Data(normalizationCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/comparison1D_2D_fit/norm"+normalizationCut+".eps");
}  

void D2hhmumuFitter_Applications::compare_misID_shapes(){

  TString PathToFolder = "/work/mitzel/D2hhmumu/dev/D2KKmumu/img/comparison_misID_shapes/";
  TString dataCut="BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5";
  TString pathToKpimumuMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root";
  TString pathToKpipipiMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpipipi_D2KKmumuBDT.root";

    
  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToSignalMC(pathToKpimumuMC);
  myFitter1D.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D.setPathToNormData(pathToNormData);
  //fix MC signal shape
  myFitter1D.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
  myFitter1D.ResolutionScale.setVal(1);myFitter1D.ResolutionScale.setConstant();
  myFitter1D.globalShift.setVal(1);myFitter1D.globalShift.setConstant();
  //single ID to pi-
  myFitter1D.fit_Kpipipi_misID("BDT>0.5&&mu1_ProbNNmu>0.5",true,PathToFolder+"Kipipi_singleMuonID_piminus.eps");
  myFitter1D.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_singleMuonID_piminus.eps");
 
  D2hhmumuFitter1D myFitter1D2;
  myFitter1D2.setPathToSignalMC(pathToKpimumuMC);
  myFitter1D2.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D2.setPathToNormData(pathToNormData);
  //fix MC signal shape
  myFitter1D2.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
  myFitter1D2.ResolutionScale.setVal(1);myFitter1D2.ResolutionScale.setConstant();
  myFitter1D2.globalShift.setVal(1);myFitter1D2.globalShift.setConstant();
  //lowPT pi ID  
  myFitter1D2.fit_Kpipipi_misID("BDT>0.5 && ( ( (mu1_PT < mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT > mu0_PT) && mu0_ProbNNmu>0.5))",true,PathToFolder+"Kipipi_lowPTMuonID.eps");
  myFitter1D2.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_lowPTMuonID.eps");

  D2hhmumuFitter1D myFitter1D3;
  myFitter1D3.setPathToSignalMC(pathToKpimumuMC);
  myFitter1D3.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D3.setPathToNormData(pathToNormData);
  //fix MC signal shape
  myFitter1D3.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
  myFitter1D3.ResolutionScale.setVal(1);myFitter1D3.ResolutionScale.setConstant();
  myFitter1D3.globalShift.setVal(1);myFitter1D3.globalShift.setConstant();
  //highPT pi ID  
  myFitter1D3.fit_Kpipipi_misID("BDT>0.5 && ( ( (mu1_PT > mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT < mu0_PT) && mu0_ProbNNmu>0.5))",true,PathToFolder+"Kipipi_highPTMuonID.eps");
  myFitter1D3.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_highPTMuonID.eps");
  
  D2hhmumuFitter1D myFitter1D4;
  myFitter1D4.setPathToSignalMC(pathToKpimumuMC);
  myFitter1D4.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D4.setPathToNormData(pathToNormData);
  myFitter1D4.setPathToKpipipiHistoData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/zoz-5000.root"); 
  //fix MC signal shape
  myFitter1D4.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
  myFitter1D4.ResolutionScale.setVal(1);myFitter1D4.ResolutionScale.setConstant();
  myFitter1D4.globalShift.setVal(1);myFitter1D4.globalShift.setConstant();
  //Benoit shape  
  myFitter1D4.fit_Kpipipi_misID_fromHistogramm("BDT>0.5&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5",true,PathToFolder+"Kipipi_BenoitShape.eps");
  myFitter1D4.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_BenoitShape.eps");

  D2hhmumuFitter1D myFitter1D5;
  myFitter1D5.setPathToSignalMC(pathToKpimumuMC);
  myFitter1D5.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D5.setPathToNormData(pathToNormData);
  //fix MC signal shape
  myFitter1D5.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
  myFitter1D5.ResolutionScale.setVal(1);myFitter1D5.ResolutionScale.setConstant();
  myFitter1D5.globalShift.setVal(1);myFitter1D5.globalShift.setConstant();
  //flat cobminatorial bkg  
  myFitter1D5.D0_M_chebyA.setVal(0);
  myFitter1D5.D0_M_chebyA.setConstant();
  myFitter1D5.fit_Kpipipi_misID("BDT>0.5 && mu1_ProbNNmu>0.5",true,PathToFolder+"Kipipi_flatCombBkg.eps");
  myFitter1D5.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_flatCombBkg.eps");
  
  D2hhmumuFitter1D myFitter1D7;
  myFitter1D7.setPathToSignalMC(pathToKpimumuMC);
  myFitter1D7.setPathToKpipipiData(pathToKpipipiMC);
  myFitter1D7.setPathToNormData(pathToNormData);
  //fix MC signal shape                                                                                                                                                                    
  myFitter1D7.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
  myFitter1D7.ResolutionScale.setVal(1);myFitter1D7.ResolutionScale.setConstant();
  myFitter1D7.globalShift.setVal(1);myFitter1D7.globalShift.setConstant();
  //Kpipipi MC                                                                                                                                                                           
  myFitter1D7.fit_Kpipipi_misID("BDT>0.5&&mu1_ProbNNmu>0.5",true,PathToFolder+"Kipipi_MC.eps");
  myFitter1D7.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_KpipipiMC.eps");
  
}  
void D2hhmumuFitter_Applications::compare_misID_shapes_2D(){

  TString PathToFolder = "/work/mitzel/D2hhmumu/dev/D2KKmumu/img/comparison_misID_shapes/2D/";
  TString dataCut="BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5";
  TString pathToKpimumuMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root";
  TString pathToKpipipiMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples_backup/MC_D2Kpipipi_D2KKmumuBDT_temp.root";

  /*  
  D2hhmumuFitter myFitter;
  myFitter.setPathToSignalMC(pathToKpimumuMC);
  myFitter.setPathToKpipipiData(pathToKpipipiData);
  myFitter.setPathToNormData(pathToNormData);
  //fix MC signal shape
  myFitter.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
   //single ID to pi-
  myFitter.fit_Kpipipi_misID("BDT>0.5&&mu1_ProbNNmu>0.5",true,PathToFolder+"Kipipi_singleMuonID_piminus.eps");
  myFitter.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_singleMuonID_piminus.eps");
  
  D2hhmumuFitter myFitter2;
  myFitter2.setPathToSignalMC(pathToKpimumuMC);
  myFitter2.setPathToKpipipiData(pathToKpipipiData);
  myFitter2.setPathToNormData(pathToNormData);
  //fix MC signal shape
  myFitter2.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
  //lowPT pi ID  
  myFitter2.fit_Kpipipi_misID("BDT>0.5 && ( ( (mu1_PT < mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT > mu0_PT) && mu0_ProbNNmu>0.5))",true,PathToFolder+"Kipipi_lowPTMuonID.eps");
  myFitter2.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_lowPTMuonID.eps");
  

  D2hhmumuFitter myFitter3;
  myFitter3.setPathToSignalMC(pathToKpimumuMC);
  myFitter3.setPathToKpipipiData(pathToKpipipiData);
  myFitter3.setPathToNormData(pathToNormData);
  //fix MC signal shape
  myFitter3.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
  //highPT pi ID  
  myFitter3.fit_Kpipipi_misID("BDT>0.5 && ( ( (mu1_PT > mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT < mu0_PT) && mu0_ProbNNmu>0.5))",true,PathToFolder+"Kipipi_highPTMuonID.eps");
  myFitter3.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_highPTMuonID.eps");

  
  D2hhmumuFitter myFitter5;
  myFitter5.setPathToSignalMC(pathToKpimumuMC);
  myFitter5.setPathToKpipipiData(pathToKpipipiData);
  myFitter5.setPathToNormData(pathToNormData);
  //fix MC signal shape
  myFitter5.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
  //flat cobminatorial bkg  
  myFitter5.D0_M_chebyA.setVal(0);
  myFitter5.D0_M_chebyA.setConstant();
  myFitter5.fit_Kpipipi_misID("BDT>0.5 && mu1_ProbNNmu>0.5",true,PathToFolder+"Kipipi_flatCombBkg.eps");
  myFitter5.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_flatCombBkg.eps");
  */

  D2hhmumuFitter myFitter7;
  myFitter7.setPathToSignalMC(pathToKpimumuMC);
  myFitter7.setPathToKpipipiData(pathToKpipipiMC);
  myFitter7.setPathToNormData(pathToNormData);
  //fix MC signal shape                                                                                                                                                                    
  myFitter7.fit_MC(dataCut,true,PathToFolder+"MCShape.eps");
  //Kpipipi MC                                                                                                                                                                           
  myFitter7.fit_Kpipipi_misID("BDT>0.5&&mu1_ProbNNmu>0.5",true,PathToFolder+"Kipipi_MC.eps");
  myFitter7.fit_normalization_Data("BDT>0.5&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5",PathToFolder+"norm_KpipipiMC.eps");
    
}  


void D2hhmumuFitter_Applications::ExtractExpectedLimit(TString fIn){

  // StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 0, 3, true, 10,1e-9,5e-8,20000);
   StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 2, 3, true, 20,1e-9,5e-8);                                                                    
}

void D2hhmumuFitter_Applications::ExtractAllExpectedLimit(TString dataCut){

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {
    TString fIn = targetFolder+"ModelConfigs/"+ dataCut + "_" + *it + ".root";
    ExtractExpectedLimit(fIn);
  }
}


void D2hhmumuFitter_Applications::runFull1DFits(TString dataCut,TString misIDCut,TString q2Cut){

  D2hhmumuFitter1D myFitter1D;

  if(kind=="D2KKmumu"){
    myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root");
    myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
    myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
    myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root");
    myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
    myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  }

  if(kind=="D2pipimumu"){
    myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2pipimumu_BDT.root");
    myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2pipimumuBDT.root");
    myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2pipimumuBDT.root");
    myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipipipi_PIDline_D2pipimumuBDT.root");
    myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2pipimumuBDT.root");
    myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_BDT.root");
  }

  std::cout<<dataCut+"&&"+q2Cut<<std::endl;
  myFitter1D.fit_MC(dataCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/MC"+dataCut+q2Cut+".eps");
  myFitter1D.fit_normalization_MC(dataCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/MC"+dataCut+q2Cut+".eps");
  
  myFitter1D.fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/misID_Kpi"+misIDCut+q2RangeNormalizationMode+".eps");
  myFitter1D.fit_HHpipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/misID_Kpi"+misIDCut+q2Cut+".eps");

  myFitter1D.fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/norm_noShapeConstraint_"+dataCut+q2RangeNormalizationMode+".eps");
  myFitter1D.fit_Data(dataCut+"&&"+q2Cut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/data_"+dataCut+q2Cut+".eps");

  //myFitter1D.fit_HHpipi_misID(misIDCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/Fits/misID_KK_"+misIDCut+q2Cut+".eps");
}

void D2hhmumuFitter_Applications::runFullResonant1DFits(TString dataCut,TString misIDCut,TString q2Cut){

  D2hhmumuFitter1D myFitter1D;

  if(kind=="D2KKmumu"){
    myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root");
    myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
    myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
    myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root");
    myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
    myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  }

  if(kind=="D2pipimumu"){
    myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2pipimumu_BDT.root");
    myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2pipimumuBDT.root");
    myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2pipimumuBDT.root");
    myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipipipi_PIDline_D2pipimumuBDT.root");
    myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2pipimumuBDT.root");
    myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_BDT.root");
  }

  std::cout<<dataCut+"&&"+q2Cut<<std::endl;
  myFitter1D.fit_MC(dataCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/MC_"+dataCut+q2Cut+".eps");
  myFitter1D.fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/MC_normMode_"+dataCut+".eps");
  std::cout<<"Monte Carlo fits done.."<<std::endl;
  
  myFitter1D.fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/misID_Kpi"+misIDCut+q2RangeNormalizationMode+".eps");
  myFitter1D.fit_HHpipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/misID_pipi"+misIDCut+q2Cut+".eps");
  std::cout<<"misID fits done.."<<std::endl;

  myFitter1D.fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/normMode_"+dataCut+q2RangeNormalizationMode+".eps");
  myFitter1D.fit_resonant_Data(dataCut+"&&"+q2Cut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/resonant_data_"+dataCut+q2Cut+".eps");
  std::cout<<"data fits done.."<<std::endl;


}

//perform fits in all q^2 bins
void D2hhmumuFitter_Applications::runAllFull1DFits(TString dataCut,TString misIDCut){

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {
    runFull1DFits(dataCut,misIDCut,*it);
  }

}

void D2hhmumuFitter_Applications::runAllResonantFull1DFits(TString dataCut,TString misIDCut){


    runFullResonant1DFits(dataCut,misIDCut,q2Ranges[1]);
    runFullResonant1DFits(dataCut,misIDCut,q2Ranges[2]);
    runFullResonant1DFits(dataCut,misIDCut,q2Ranges[3]);

}

void D2hhmumuFitter_Applications::studyResolutionScale(TString dataCut,TString misIDCut,TString q2Cut,TString channel){

  D2hhmumuFitter1D* myFitter1D = new D2hhmumuFitter1D();
  myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_"+channel+"_D2KKmumuBDT.root");
  myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/"+channel+"_PIDline_D2KKmumuBDT.root");
  if(channel=="D2Kpipipi") myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/"+channel+"_D2KKmumuBDT_even1.root");
  std::cout<<"study MC resolution "<<channel<<"  "<<dataCut+"&&"+q2Cut<<std::endl; 
  myFitter1D->fit_MC(dataCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyResolutionScale/MC"+channel+"_"+q2Cut+".eps");
  myFitter1D->fit_HHpipi(dataCut+"&&nTracks%6==0"+"&&"+q2Cut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyResolutionScale/Data"+channel+"_"+q2Cut+".eps");
  //delete myFitter1D;
}

//perform fits in all q^2 bins
void D2hhmumuFitter_Applications::runAllResolutionScaleStudies(TString dataCut,TString misIDCut){
  

  //  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {
  //  studyResolutionScale(dataCut,misIDCut,*it,"D2KKpipi");
  // }
  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {
    studyResolutionScale(dataCut,misIDCut,*it,"D2Kpipipi");
  }

}


void D2hhmumuFitter_Applications::performAllToyStudies(){

  TString dataCut="BDT>0&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5";
  TString misIDCut="BDT>0.95&&mu1_ProbNNmu>0.5";
  
  D2hhmumuFitter1D myFitter1D;
    
  myFitter1D.setPathToSignalMC(pathToKpimumuMC);
  myFitter1D.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D.setPathToHHpipiData(pathToKKpipiData);
  myFitter1D.setPathToSignalData(pathToSignalData);
  myFitter1D.setPathToNormData(pathToNormData);

  double nSig[4] = {2240,30,1,2}; //{norm,reso,eta,low}
  double nBkg[4] = {107,10,1,10};
  double nMisID[4] = {378,4,1,1};

  /*
  myFitter1D.makeToyStudy(dataCut,q2Ranges[1],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_EtaBin.eps",nSig[2],nBkg[2],nMisID[2],false,false);
  myFitter1D.makeToyStudy(dataCut,q2Ranges[0],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_lowMassBin.eps",nSig[3],nBkg[3],nMisID[3],false,false);  
  myFitter1D.makeToyStudy(dataCut,q2Ranges[2],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_resonantBin.eps",nSig[1],nBkg[1],nMisID[1],false,false);
  myFitter1D.makeToyStudy(dataCut,q2RangeNormalizationMode,misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_normalization.eps",nSig[0],nBkg[0],nMisID[0],false,false);//toys are thrown with the expected yields    
  */
  //shape of comb free, more unstable
  myFitter1D.makeToyStudy(dataCut,q2Ranges[1],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_EtaBin_combBkgShapeFree.eps",nSig[2],nBkg[2],nMisID[2],false,true);
  myFitter1D.makeToyStudy(dataCut,q2Ranges[0],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_lowMassBin_combBkgShapeFree.eps",nSig[3],nBkg[3],nMisID[3],false,true);
  myFitter1D.makeToyStudy(dataCut,q2Ranges[2],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_resonantBin_combBkgShapeFree.eps",nSig[1],nBkg[1],nMisID[1],false,true);
  myFitter1D.makeToyStudy(dataCut,q2RangeNormalizationMode,misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_normalization_combBkgShapeFree.eps",nSig[0],nBkg[0],nMisID[0],false,true);

  /* set misID bkg to 0 
  myFitter1D.makeToyStudy(dataCut,q2Ranges[1],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_EtaBin_misID20.eps",nSig[2],nBkg[2],nMisID[2],true,true);
  myFitter1D.makeToyStudy(dataCut,q2Ranges[0],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_lowMassBin_misID20.eps",nSig[3],nBkg[3],nMisID[3],true,true);
  myFitter1D.makeToyStudy(dataCut,q2Ranges[2],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_resonantBin_misID20.eps",nSig[1],nBkg[1],nMisID[1],true,true);
  myFitter1D.makeToyStudy(dataCut,q2RangeNormalizationMode,misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_normalization_misID20.eps",nSig[0],nBkg[0],nMisID[0],false,true);
  */
}

void D2hhmumuFitter_Applications::studyNormalizationFits(TString dataCut,TString misIDCut, bool doNSharedCut){

  TString path = "/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyNormalizationFits/"; //where to safe results
  TString cutNshared=""; 
  if(doNSharedCut) {
    cutNshared= "&&mu0_MuonNShared==0&&mu1_MuonNShared==0";
    path=path+"noSharedMuonHits_";
  }
  dataCut=dataCut+cutNshared;

  //fix shapes,scale factors
  
  D2hhmumuFitter1D myFitter1D;
    
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_D2KKmumuBDT_even1.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_D2KKmumuBDT_even1.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.fit_MC(dataCut,true,path+"MCfit_Kpimumu.eps");
  myFitter1D.fit_Kpipipi_misID(misIDCut,true,path+"misIDfit_Kpipipi.eps");
  myFitter1D.ResolutionScale.setVal(1.097);                                                                                                                                 
  myFitter1D.ResolutionScale.setConstant();
  myFitter1D.globalShift.setVal(1.00055);
  myFitter1D.globalShift.setConstant();
  myFitter1D.D0_M_chebyB.setVal(0);
  myFitter1D.D0_M_chebyB.setConstant();
  myFitter1D.fit_normalization_Data(dataCut,path+"fixedMisID_scaleFactorFromKpipipi.eps");

  //fix shapes,scale factors free
  
  D2hhmumuFitter1D myFitter1D_2;
    
  myFitter1D_2.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D_2.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_D2KKmumuBDT_even1.root");
  myFitter1D_2.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_D2KKmumuBDT_even1.root");
  myFitter1D_2.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D_2.fit_MC(dataCut+"&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875",true,"");
  myFitter1D_2.fit_Kpipipi_misID(misIDCut,true,"");
  myFitter1D_2.D0_M_chebyB.setVal(0);
  myFitter1D_2.D0_M_chebyB.setConstant();  
  myFitter1D_2.fit_normalization_Data(dataCut,path+"fixedMisID_scaleFactorFree.eps");

  //fix shapes,scale factors fix to results from D2Kpimumu fit 
   
  D2hhmumuFitter1D myFitter1D_3;
    
  myFitter1D_3.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D_3.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_D2KKmumuBDT_even1.root");
  myFitter1D_3.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_D2KKmumuBDT_even1.root");
  myFitter1D_3.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D_3.fit_MC(dataCut+"&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875",true,"");
  myFitter1D_3.fit_Kpipipi_misID(misIDCut,true,"");
  myFitter1D_3.ResolutionScale.setVal(1.095);                                                                                                                                 
  myFitter1D_3.ResolutionScale.setConstant();
  myFitter1D_3.globalShift.setVal(1.00037);
  myFitter1D_3.globalShift.setConstant();
  myFitter1D_3.D0_M_chebyB.setVal(0);
  myFitter1D_3.D0_M_chebyB.setConstant();  
  myFitter1D_3.fit_normalization_Data(dataCut,path+"fixedMisID_scaleFactorFromKpimumu.eps");

   //misID shape free,scale factors fix 
   
  D2hhmumuFitter1D myFitter1D_4;
    
  myFitter1D_4.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D_4.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_D2KKmumuBDT_even1.root");
  myFitter1D_4.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_D2KKmumuBDT_even1.root");
  myFitter1D_4.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D_4.fit_MC(dataCut+"&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875",true,"");
  myFitter1D_4.ResolutionScale.setVal(1.097);                                                                                                                                 
  myFitter1D_4.ResolutionScale.setConstant();
  myFitter1D_4.globalShift.setVal(1.00055);
  myFitter1D_4.globalShift.setConstant();
  myFitter1D_4.D0_M_chebyB.setVal(0);
  myFitter1D_4.D0_M_chebyB.setConstant();  
  myFitter1D_4.fit_normalization_Data(dataCut,path+"freeMisID_scaleFactorFromKpipipi.eps");

  //misID shape free,scale factors free

  D2hhmumuFitter1D myFitter1D_5;
    
  myFitter1D_5.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D_5.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_D2KKmumuBDT_even1.root");
  myFitter1D_5.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_D2KKmumuBDT_even1.root");
  myFitter1D_5.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D_5.fit_MC(dataCut+"&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875",true,"");
  myFitter1D_5.D0_M_chebyB.setVal(0);
  myFitter1D_5.D0_M_chebyB.setConstant();  
  myFitter1D_5.fit_normalization_Data(dataCut,path+"freeMisID_freeScaleFactor.eps");
  

}

