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

  if(m_kind!= "D2KKmumu") std::cout<<" D2hhmumuFitter_Applications: Error decay mode not implemeted"<<std::endl;
  if(m_kind== "D2KKmumu"){
    pathToSignalMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT.root";
    pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root";
    pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root";
    pathToSidebandData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT.root";
    pathToKpipipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root";
    pathToKKpipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root"; 
    pathToKpimumuMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root";
  }
}


void D2hhmumuFitter_Applications::setQ2Ranges(TString m_kind){
 
  if(m_kind!= "D2KKmumu") std::cout<<" D2hhmumuFitter_Applications: Error decay mode not implemeted"<<std::endl;
  
  if(m_kind== "D2KKmumu"){
    q2Ranges.push_back("D_DiMuon_Mass<525");
    q2Ranges.push_back("D_DiMuon_Mass>525&&D_DiMuon_Mass<565");
    q2Ranges.push_back("D_DiMuon_Mass>565");
  }   
  std::cout<<"q^2 intervals: "<<std::endl;
  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) { 
    std::cout<<*it<<" , ";
  }
  std::cout<<" "<<std::endl;
}

void D2hhmumuFitter_Applications::saveModelConfig(TString dataCut,TString normalizationCut,TString misIDCut,TString q2Range){

  D2hhmumuFitter1D myFitter1D;                                                                                                                                                         
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKKpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root");

  TString fOut=targetFolder+"ModelConfigs/"+ dataCut + "_" + q2Range + ".root"; 
  TString totalDataCut=dataCut+"&&"+q2Range;
  TString totalNormalizationCut=normalizationCut+"&&"+q2RangeNormalizationMode; //
  TString totalMisIDCut=misIDCut+"&&"+q2Range;

  myFitter1D.fillModelConfig(totalDataCut,totalNormalizationCut,totalMisIDCut,fOut); 

}

void D2hhmumuFitter_Applications::saveAllModelConfig(TString dataCut,TString nomalizationCut,TString misIDCut){

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

    saveModelConfig(dataCut,nomalizationCut,misIDCut,*it);
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

  /*  
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
  */
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


void D2hhmumuFitter_Applications::runFull1DFits(TString dataCut,TString normalizationCut,TString misIDCut,TString q2Cut){

  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter1D.setPathToKKpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_PIDline_D2KKmumuBDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");

  std::cout<<dataCut+"&&"+q2Cut<<std::endl;
  myFitter1D.fit_MC(dataCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/Fits/MC"+dataCut+q2Cut+".eps");
  myFitter1D.fit_Kpipipi_misID(misIDCut+"&&"+q2RangeNormalizationMode,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/Fits/misID_Kpi"+misIDCut+q2RangeNormalizationMode+".eps");
  myFitter1D.fit_normalization_Data(normalizationCut+"&&"+q2RangeNormalizationMode,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/Fits/norm"+normalizationCut+q2RangeNormalizationMode+".eps");

  myFitter1D.fit_HHpipi_misID(misIDCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/Fits/misID_KK_"+misIDCut+q2Cut+".eps");
  myFitter1D.fit_Data(dataCut+"&&"+q2Cut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/Fits/data_"+misIDCut+q2Cut+".eps");
}

//perform fits in all q^2 bins
void D2hhmumuFitter_Applications::runAllFull1DFits(TString dataCut,TString normalizationCut,TString misIDCut){

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {
    runFull1DFits(dataCut,normalizationCut,misIDCut,*it);
  }

}


void D2hhmumuFitter_Applications::performAllToyStudies(){

  TString dataCut="BDT>0.6&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5";
  TString normalizationCut="BDT>0.6&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5";
  TString misIDCut="BDT>0.6&&mu1_ProbNNmu>0.5";
  
  D2hhmumuFitter1D myFitter1D;
  
  
  myFitter1D.setPathToSignalMC(pathToKpimumuMC);
  myFitter1D.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D.setPathToKKpipiData(pathToKKpipiData);
  myFitter1D.setPathToSignalData(pathToSignalData);
  myFitter1D.setPathToNormData(pathToNormData);

  double nSig[4] = {2300,60,5,3}; //{norm,reso,eta,low}
  double nBkg[4] = {83,13,1,28};
  double nMisID[4] = {391,3,1,2};
  
  
  myFitter1D.makeToyStudy(dataCut+"&&"+q2RangeNormalizationMode,normalizationCut+"&&"+q2RangeNormalizationMode, //cut values are needed to fix the shapes for the toy studies
			  misIDCut+"&&"+q2RangeNormalizationMode,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/normalization.eps",nSig[0],nBkg[0],nMisID[0],false);//toys are thrown with the expected yields
    myFitter1D.makeToyStudy(dataCut+"&&"+q2Ranges[2],normalizationCut+"&&"+normalizationCut,
			    misIDCut+"&&"+q2Ranges[2],"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/resonantBinn.eps",nSig[1],nBkg[1],nMisID[1],false);
    myFitter1D.makeToyStudy(dataCut+"&&"+q2Ranges[1],normalizationCut+"&&"+normalizationCut, 
			    misIDCut+"&&"+q2Ranges[1],"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/etaBin.eps",nSig[2],nBkg[2],nMisID[2],false);
    myFitter1D.makeToyStudy(dataCut+"&&"+q2Ranges[0],normalizationCut+"&&"+normalizationCut, 
			    misIDCut+"&&"+q2Ranges[0],"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/lowmassBin.eps",nSig[3],nBkg[3],nMisID[3],false);

    //again with misID background set to 0, check bias on signal yield
    myFitter1D.makeToyStudy(dataCut+"&&"+q2RangeNormalizationMode,normalizationCut+"&&"+q2RangeNormalizationMode, //cut values are needed to fix the shapes for the toy studies
			  misIDCut+"&&"+q2RangeNormalizationMode,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/normalization_misID20.eps",nSig[0],nBkg[0],nMisID[0],true);//toys are thrown with the expected yields
    myFitter1D.makeToyStudy(dataCut+"&&"+q2Ranges[2],normalizationCut+"&&"+normalizationCut,
			    misIDCut+"&&"+q2Ranges[2],"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/resonantBinn_misID20.eps",nSig[1],nBkg[1],nMisID[1],true);
    myFitter1D.makeToyStudy(dataCut+"&&"+q2Ranges[1],normalizationCut+"&&"+normalizationCut, 
			    misIDCut+"&&"+q2Ranges[1],"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/etaBin_misID20.eps",nSig[2],nBkg[2],nMisID[2],true);
    myFitter1D.makeToyStudy(dataCut+"&&"+q2Ranges[0],normalizationCut+"&&"+normalizationCut, 
			    misIDCut+"&&"+q2Ranges[0],"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/lowmassBin_misID20.eps",nSig[3],nBkg[3],nMisID[3],true);

}

