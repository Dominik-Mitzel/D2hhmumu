#include "D2hhmumuFitter_Applications.h"
#include "D2hhmumuFitter1D.h"
#include "StandardHypoTestInvDemo.C"
#include "EfficiencyCalculator.h"
#include "dcastyle.h"


D2hhmumuFitter_Applications::D2hhmumuFitter_Applications(TString m_kind,TString m_year){

  dcastyle();
  kind= m_kind;
  year= m_year;
  if(m_kind== "D2KKmumu") 
    targetFolder="/work/mitzel/D2hhmumu/dev/D2KKmumu/";
  if(m_kind== "D2pipimumu")
    targetFolder="/work/mitzel/D2hhmumu/dev/D2pipimumu/";
  setDefaultPathToData(kind);  
  setQ2Ranges(kind);
  q2RangeNormalizationMode="D_DiMuon_Mass>675&&D_DiMuon_Mass<875";
 
}

D2hhmumuFitter_Applications::~D2hhmumuFitter_Applications(){};

void D2hhmumuFitter_Applications::setDefaultPathToData(TString m_kind){

  if(! (m_kind=="D2KKmumu" || m_kind=="D2pipimumu") ) std::cout<<" D2hhmumuFitter_Applications: Error decay mode not implemeted"<<std::endl;
  if(m_kind== "D2KKmumu"){
    pathToSignalMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT.root";
    pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root";
    pathToSidebandData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/sideband_D2KKmumu_BDT.root";
    pathToKpipipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_new.root";
    pathToKKpipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT_new.root";
    pathToKpimumuMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root";
    pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT.root";
    pathToNormMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root";
  }

  if(m_kind== "D2pipimumu"){
    pathToSignalMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root";
    pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root";
    pathToSidebandData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/sideband_D2pipimumu_BDT.root";
    pathToKpipipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_new.root";
    pathToKKpipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_new.root";
    pathToKpimumuMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root";
    pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root";
    pathToNormMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root";
  }


}


void D2hhmumuFitter_Applications::setQ2Ranges(TString m_kind){
 
  if(! (m_kind=="D2KKmumu" || m_kind=="D2pipimumu")) std::cout<<" D2hhmumuFitter_Applications: Error decay mode not implemeted"<<std::endl;
  
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

void D2hhmumuFitter_Applications::saveModelConfig(TString dataCut,TString misIDCut){

  D2hhmumuFitter1D* myFitter1D;

  for(int q2Bin=0; q2Bin<q2Ranges.size();++q2Bin) {

    myFitter1D = new D2hhmumuFitter1D();

    //omit the resonant bins
    if(q2Bin==2 || q2Bin==3) continue;
 
   
    TString fOut= targetFolder+"ModelConfigs/"+ kind + "_" + TString::Format("bin_%i",q2Bin)+".root"; 
    std::cout<<q2Bin<<"  "<<fOut<<std::endl;

    if(kind=="D2KKmumu"){
      myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT.root");
      myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root");
      myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_new.root");
      myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT_new.root");
      myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root");
      myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT.root");
    }
    
    if(kind=="D2pipimumu"){
      myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root");
      myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");
      myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_new.root");
      myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_new.root");
      myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root");
      myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root");
    }

    myFitter1D->fillModelConfig(kind,dataCut,misIDCut,q2Bin,fOut); 

  }

}

 
void D2hhmumuFitter_Applications::compare_1D_and_2D_fit(TString dataCut,TString normalizationCut,TString misIDCut){

  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2KKmumuBDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.fit_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/comparison1D_2D_fit/MC"+dataCut+".eps","","");
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

void D2hhmumuFitter_Applications::compare_misID_shapes(TString dataCut="BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5"){

  TString PathToFolder = "/work/mitzel/D2hhmumu/dev/"+kind+"/img/comparison_misID_shapes/freeCombinatorialBkgShape/";


  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToNormMC(pathToNormMC);
  myFitter1D.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D.setPathToNormData(pathToNormData);
  myFitter1D.setPathToSignalData(pathToSignalData);
  myFitter1D.setPathToSignalMC(pathToSignalMC);
  myFitter1D.setPathToHHpipiData(pathToKKpipiData);
  //fix MC signal shape
  myFitter1D.fit_normalization_MC(dataCut,true,PathToFolder+"MCShape.eps");
  //myFitter1D.ResolutionScale.setVal(1);myFitter1D.ResolutionScale.setConstant();
  //myFitter1D.globalShift.setVal(1);myFitter1D.globalShift.setConstant();
  //single ID to pi-
  myFitter1D.fit_HHpipi_misID("BDT>0.4&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"pipipipi_singleMuonID_piminus.eps");
  myFitter1D.fit_Kpipipi_misID("BDT>0.4&&mu1_ProbNNmu>0.5",true,PathToFolder+"Kpipipi_singleMuonID_piminus.eps");
  myFitter1D.fit_normalization_Data(dataCut,PathToFolder+"norm_singleMuonID_piminus.eps");
  myFitter1D.fit_MC(dataCut+"&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"MCSignalShape.eps");
  myFitter1D.fit_resonant_Data("D2pipimumu",dataCut,"D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"signal_singleMuonID_piminus.eps");

  D2hhmumuFitter1D myFitter1D_2;
  myFitter1D_2.setPathToNormMC(pathToNormMC);
  myFitter1D_2.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D_2.setPathToNormData(pathToNormData);
  myFitter1D_2.setPathToSignalData(pathToSignalData);
  myFitter1D_2.setPathToSignalMC(pathToSignalMC);
  myFitter1D_2.setPathToHHpipiData(pathToKKpipiData);
  //fix MC signal shape
  myFitter1D_2.fit_normalization_MC(dataCut,true,PathToFolder+"MCShape.eps");
  //myFitter1D.ResolutionScale.setVal(1);myFitter1D.ResolutionScale.setConstant();
  //myFitter1D.globalShift.setVal(1);myFitter1D.globalShift.setConstant();
  //single ID to pi+
  myFitter1D_2.fit_HHpipi_misID("BDT>0.4&&mu0_ProbNNmu>0.5&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"pipipipi_singleMuonID_piplus.eps");
  myFitter1D_2.fit_Kpipipi_misID("BDT>0.4&&mu0_ProbNNmu>0.5",true,PathToFolder+"Kpipipi_singleMuonID_piplus.eps");
  myFitter1D_2.fit_normalization_Data(dataCut,PathToFolder+"norm_singleMuonID_piplus.eps");
  myFitter1D_2.fit_MC(dataCut+"&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"MCSignalShape.eps");
  myFitter1D_2.fit_resonant_Data("D2pipimumu",dataCut,"D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"signal_singleMuonID_piplus.eps");

  D2hhmumuFitter1D myFitter1D_3;
  myFitter1D_3.setPathToNormMC(pathToNormMC);
  myFitter1D_3.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D_3.setPathToNormData(pathToNormData);
  myFitter1D_3.setPathToSignalData(pathToSignalData);
  myFitter1D_3.setPathToSignalMC(pathToSignalMC);
  myFitter1D_3.setPathToHHpipiData(pathToKKpipiData);
  //fix MC signal shape
  myFitter1D_3.fit_normalization_MC(dataCut,true,PathToFolder+"MCShape.eps");
  //myFitter1D.ResolutionScale.setVal(1);myFitter1D.ResolutionScale.setConstant();
  //myFitter1D.globalShift.setVal(1);myFitter1D.globalShift.setConstant();
  //low pT pi
  myFitter1D_3.fit_HHpipi_misID("BDT>0.4 && ( ( (mu1_PT < mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT > mu0_PT) && mu0_ProbNNmu>0.5))&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"pipipipi_lowPTMuonID.eps");
  myFitter1D_3.fit_Kpipipi_misID("BDT>0.4 && ( ( (mu1_PT < mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT > mu0_PT) && mu0_ProbNNmu>0.5))",true,PathToFolder+"Kpipipi_lowPTMuonID.eps");
  myFitter1D_3.fit_normalization_Data(dataCut,PathToFolder+"norm_lowPTMuonID.eps");
  myFitter1D_3.fit_MC(dataCut+"&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"MCSignalShape.eps");
  myFitter1D_3.fit_resonant_Data("D2pipimumu",dataCut,"D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"signal_lowPTMuonID.eps");


  D2hhmumuFitter1D myFitter1D_4;
  myFitter1D_4.setPathToNormMC(pathToNormMC);
  myFitter1D_4.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D_4.setPathToNormData(pathToNormData);
  myFitter1D_4.setPathToSignalData(pathToSignalData);
  myFitter1D_4.setPathToSignalMC(pathToSignalMC);
  myFitter1D_4.setPathToHHpipiData(pathToKKpipiData);
  //fix MC signal shape
  myFitter1D_4.fit_normalization_MC(dataCut,true,PathToFolder+"MCShape.eps");
  //myFitter1D.ResolutionScale.setVal(1);myFitter1D.ResolutionScale.setConstant();
  //myFitter1D.globalShift.setVal(1);myFitter1D.globalShift.setConstant();
  //low pT pi
  myFitter1D_4.fit_HHpipi_misID("BDT>0.4 && ( ( (mu1_PT > mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT < mu0_PT) && mu0_ProbNNmu>0.5))&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"pipipipi_highPTMuonID.eps"); 
  myFitter1D_4.fit_Kpipipi_misID("BDT>0.4 && ( ( (mu1_PT > mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT < mu0_PT) && mu0_ProbNNmu>0.5))",true,PathToFolder+"Kpipipi_lowPTMuonID.eps");
  myFitter1D_4.fit_normalization_Data(dataCut,PathToFolder+"norm_highPTMuonID.eps");
  myFitter1D_4.fit_MC(dataCut+"&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"MCSignalShape.eps");
  myFitter1D_4.fit_resonant_Data("D2pipimumu",dataCut,"D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"signal_highPTMuonID.eps");
  

  
  D2hhmumuFitter1D myFitter1D_5;
  myFitter1D_5.setPathToNormMC(pathToNormMC);
  myFitter1D_5.setPathToKpipipiData(pathToKpipipiData);
  myFitter1D_5.setPathToNormData(pathToNormData);
  myFitter1D_5.setPathToSignalData(pathToSignalData);
  myFitter1D_5.setPathToSignalMC(pathToSignalMC);
  myFitter1D_5.setPathToHHpipiData(pathToKKpipiData);
  //fix MC signal shape
  myFitter1D_5.fit_normalization_MC(dataCut,true,PathToFolder+"MCShape.eps");
  //myFitter1D_5.ResolutionScale.setVal(1);myFitter1D_5.ResolutionScale.setConstant();
  //myFitter1D_5.globalShift.setVal(1);myFitter1D_5.globalShift.setConstant();
  //single ID to pi-
  myFitter1D_5.fit_HHpipi_misID("BDT>0.4&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"pipipipi_noMounID.eps");
  myFitter1D_5.fit_Kpipipi_misID("BDT>0.4&&nTracks%15==0",true,PathToFolder+"Kpipipi_noMuonID.eps");
  myFitter1D_5.fit_normalization_Data(dataCut,PathToFolder+"norm_noMuonID.eps");
  myFitter1D_5.fit_MC(dataCut+"&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"MCSignalShape.eps");
  myFitter1D_5.fit_resonant_Data("D2pipimumu",dataCut,"D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"signal_noMuonID.eps");
  

  /*  
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
  */

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

void D2hhmumuFitter_Applications::ExtractExpectedLimit(){
  

  //StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 0, 3, true, 20,1e-9,5e-8,20000);
  //good range for KKmumu //StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 2, 3, true, 20,1e-9,7e-8);                                                           

  TString fIn;  

  for(int q2Bin=0; q2Bin<q2Ranges.size();++q2Bin) {

    fIn = targetFolder+"ModelConfigs/"+ kind + "_" + TString::Format("bin_%i",q2Bin)+".root";

    StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 2, 3, true, 20,1e-9,7e-8);                                                                   
  //StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 2, 3, true, 20,1e-8,5e-7);                 
  }

}



void D2hhmumuFitter_Applications::runFull1DFits(TString dataCut,TString misIDCut){

  D2hhmumuFitter1D* myFitter1D;
  int counter=0;

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

    myFitter1D= new D2hhmumuFitter1D();

    if(kind=="D2KKmumu"){
      myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT.root");
      myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root");
      myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_new.root");
      myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT_new.root");
      myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root");
      myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT.root");


      myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_KKmumu_bin_%i.eps",counter),"m(KK#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangesKK_low[counter],rangesKK_high[counter]));
      myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_norm_KKTrained.eps");
      std::cout<<"Monte Carlo fits done.."<<std::endl;
      myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpipipi_KKTrained.eps");
      myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_KKmumu_bin_%i.eps",counter),"m_{KK#mu#mu}(KK#pi#pi)",TString::Format("%.0fMeV<m(#pi#pi)<%.0fMeV",rangesKK_low[counter],rangesKK_high[counter])); 
      std::cout<<"misID fits done.."<<std::endl;
      myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_norm_KKTrained.eps");
      myFitter1D->fit_Data(kind,dataCut,+(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_blind_KKmumu_bin_%i.eps",counter),"m(KK#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangesKK_low[counter],rangesKK_high[counter]),true);
      std::cout<<"blind data fits done.."<<std::endl;

      ++counter;
    
    }

    if(kind=="D2pipimumu"){
      myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root");
      myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");
      myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_new.root");
      myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_new.root");
      myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root");
      myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root");

      myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_pipimumu_bin_%i.eps",counter),"m(#pi#pi#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangespipi_low[counter],rangespipi_high[counter]));
      myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_norm_pipiTrained.eps");
      std::cout<<"Monte Carlo fits done.."<<std::endl;
      myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpipipi_pipiTrained.eps");
      myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_pipimumu_bin_%i.eps",counter),"m_{#pi#pi#mu#mu}(#pi#pi#pi#pi)",TString::Format("%.0fMeV<m(#pi#pi)<%.0fMeV",rangespipi_low[counter],rangespipi_high[counter])); 
      std::cout<<"misID fits done.."<<std::endl;
      myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_norm_pipiTrained.eps");
      myFitter1D->fit_Data(kind,dataCut,+(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_blind_pipimumu_bin_%i.eps",counter),"m(#pi#pi#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangespipi_low[counter],rangespipi_high[counter]),true);
      std::cout<<"blind data fits done.."<<std::endl;
      
      ++counter;
    }
  }
}


void D2hhmumuFitter_Applications::runFullResonant1DFits(TString dataCut,TString misIDCut){


  D2hhmumuFitter1D* myFitter1D;
  int counter=0;

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

      myFitter1D= new D2hhmumuFitter1D();

      if(kind=="D2KKmumu" && (counter!=0 && counter!=1 && counter!=2 ) ){
	myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT.root");
	myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root");
	myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_new.root");
	myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT_new.root");
	myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root");
	myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT.root");
  
	myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_KKmumu_bin_%i.eps",counter),"m(KK#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangesKK_low[counter],rangesKK_high[counter]));
	myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_norm_KKTrained.eps");
	std::cout<<"Monte Carlo fits done.."<<std::endl;
	myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpipipi_KKTrained.eps");
	myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_KKmumu_bin_%i.eps",counter),"m_{KK#mu#mu}(KK#pi#pi)",TString::Format("%.0fMeV<m(#pi#pi)<%.0fMeV",rangesKK_low[counter],rangesKK_high[counter])); 
	std::cout<<"misID fits done.."<<std::endl;
	myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_norm_KKTrained.eps");
	myFitter1D->fit_resonant_Data(kind,dataCut,+(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_KKmumu_bin_%i.eps",counter),"m(KK#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangesKK_low[counter],rangesKK_high[counter]),true);
	std::cout<<"data fits done.."<<std::endl;

      }
 
      if(kind=="D2pipimumu" && (counter!=0 && counter!=1 && counter!=4)){
	myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root");
	myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");
	myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_new.root");
	myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_new.root");
	myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root");
	myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root");

	myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_pipimumu_bin_%i.eps",counter),"m(#pi#pi#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangespipi_low[counter],rangespipi_high[counter]));
	myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_norm_pipiTrained.eps");
	std::cout<<"Monte Carlo fits done.."<<std::endl;
	myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpipipi_pipiTrained.eps");
	myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_pipimumu_bin_%i.eps",counter),"m_{#pi#pi#mu#mu}(#pi#pi#pi#pi)",TString::Format("%.0fMeV<m(#pi#pi)<%.0fMeV",rangespipi_low[counter],rangespipi_high[counter])); 
	std::cout<<"misID fits done.."<<std::endl;
	myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_norm_pipiTrained.eps");
	myFitter1D->fit_resonant_Data(kind,dataCut,+(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_pipimumu_bin_%i.eps",counter),"m(#pi#pi#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangespipi_low[counter],rangespipi_high[counter]),true);
	std::cout<<" data fits done.."<<std::endl;
      
      }
      ++counter;
  }
  
  std::cout<<"data fits done.."<<std::endl;

}


void D2hhmumuFitter_Applications::constrainCombBkgShapes(TString dataCut,TString misIDCut){

  D2hhmumuFitter1D* myFitter1D;
  int counter=0;

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

    myFitter1D= new D2hhmumuFitter1D();

    if(kind=="D2KKmumu"){
      myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT.root");
      myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root");
      myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_new.root");
      myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT_new.root");
      myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root");
      myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT.root");


      myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/MC_KKmumu_bin_%i.eps",counter),"m(KK#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangesKK_low[counter],rangesKK_high[counter]));
      myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/MC_norm_KKTrained.eps");
      std::cout<<"Monte Carlo fits done.."<<std::endl;
      myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/misID_Kpipipi_KKTrained.eps");
      myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/misID_KKmumu_bin_%i.eps",counter),"m_{KK#mu#mu}(KK#pi#pi)",TString::Format("%.0fMeV<m(#pi#pi)<%.0fMeV",rangesKK_low[counter],rangesKK_high[counter])); 
      std::cout<<"misID fits done.."<<std::endl;
      myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/data_norm_KKTrained.eps");
      //inverderd BDT CUT!!
      myFitter1D->fit_Data(kind,"BDT<.4&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>.5",+(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/data_blind_invertedBDT_KKmumu_bin_%i.eps",counter),"m(KK#mu#mu)","",false);
      std::cout<<"blind data fits done.."<<std::endl;

      ++counter;
    
    }

    if(kind=="D2pipimumu"){
      myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root");
      myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");
      myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_new.root");
      myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_new.root");
      myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root");
      myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root");

      myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/MC_pipimumu_bin_%i.eps",counter),"m(#pi#pi#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangespipi_low[counter],rangespipi_high[counter]));
      myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/MC_norm_pipiTrained.eps");
      std::cout<<"Monte Carlo fits done.."<<std::endl;
      myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/misID_Kpipipi_pipiTrained.eps");
      myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/misID_pipimumu_bin_%i.eps",counter),"m_{#pi#pi#mu#mu}(#pi#pi#pi#pi)",TString::Format("%.0fMeV<m(#pi#pi)<%.0fMeV",rangespipi_low[counter],rangespipi_high[counter])); 
      std::cout<<"misID fits done.."<<std::endl;
      myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/data_norm_pipiTrained.eps");
      myFitter1D->fit_Data(kind,"BDT<.4&&BDT>-.5&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>.5",(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/constrainCombBkgShape/data_blind_invertedBDTpipimumu_bin_%i.eps",counter),"m(#pi#pi#mu#mu)","",false);
      std::cout<<"blind data fits done.."<<std::endl;
      
      ++counter;
    }
  }
}



void D2hhmumuFitter_Applications::plotRelativeEfficiencies(TString dataCut){

  const Int_t NBINS = q2Ranges.size();
  Double_t edges[q2Ranges.size()+1]; 

  if(kind=="D2KKmumu"){
    edges[0]=0;
    edges[1]=525;
    edges[2]=565;
    edges[3]=900;
  }
  if(kind=="D2pipimumu"){
    edges[0]=0;
    edges[1]=525;
    edges[2]=565;
    edges[3]=950;
    edges[4]=1100;
    edges[5]=1450;
  }
  
  TH1* h = new TH1D("relativeEff","relative Efficiency",NBINS,edges);  

  EfficiencyCalculator myEfficiency(kind);
  double relEff;
  int counter=1;

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {
    relEff = myEfficiency.getMCRelativeSigToNormEfficiency(dataCut,"",*it,"D_DiMuon_Mass>675&&D_DiMuon_Mass<875");
    h->SetBinContent(counter,relEff);
    counter+=1;
  }
  
  TCanvas a("a","a");
  a.cd();
  h->GetYaxis()->SetRangeUser(0,2);
  h->Draw();
  a.Print(targetFolder+"img/EfficiencyStudies/RelativeEfficiencyVariableBinning.eps");

}


void D2hhmumuFitter_Applications::performAllToyStudies(){
  /*

  TString dataCut="BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5";
  TString misIDCut="BDT>0.4&&mu0_ProbNNmu>0.5";
  
  D2hhmumuFitter1D myFitter1D;

  myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT.root");
  myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT.root");

  double nSig[4] = {2240,30,1,2}; //{norm,reso,eta,low}
  double nBkg[4] = {107,10,1,10};
  double nMisID[4] = {378,4,1,1};
  */
  /*
  myFitter1D.makeToyStudy(dataCut,q2Ranges[1],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_EtaBin.eps",nSig[2],nBkg[2],nMisID[2],false,false);
  myFitter1D.makeToyStudy(dataCut,q2Ranges[0],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_lowMassBin.eps",nSig[3],nBkg[3],nMisID[3],false,false);  
  myFitter1D.makeToyStudy(dataCut,q2Ranges[2],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_resonantBin.eps",nSig[1],nBkg[1],nMisID[1],false,false);
  myFitter1D.makeToyStudy(dataCut,q2RangeNormalizationMode,misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_normalization.eps",nSig[0],nBkg[0],nMisID[0],false,false);//toys are thrown with the expected yields    
  */
  
  //shape of comb free, more unstable
  /*
  myFitter1D.makeToyStudy(dataCut,q2Ranges[1],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_EtaBin_combBkgShapeFree.eps",nSig[2],nBkg[2],nMisID[2],false,true);
  myFitter1D.makeToyStudy(dataCut,q2Ranges[0],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_lowMassBin_combBkgShapeFree.eps",nSig[3],nBkg[3],nMisID[3],false,true);
  myFitter1D.makeToyStudy(dataCut,q2Ranges[2],misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_resonantBin_combBkgShapeFree.eps",nSig[1],nBkg[1],nMisID[1],false,true);
  myFitter1D.makeToyStudy(dataCut,q2RangeNormalizationMode,misIDCut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/fitterValidation/KKmumu_normalization_combBkgShapeFree.eps",nSig[0],nBkg[0],nMisID[0],false,true);
  */
  
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

void D2hhmumuFitter_Applications::drawMisIDShapes(TString cut){

  D2hhmumuFitter1D* myFitter1D;
  int counter=0;
  
  
  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

    myFitter1D = new D2hhmumuFitter1D();

    if(kind=="D2KKmumu"){
      myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT_new.root");
      myFitter1D->fit_HHpipi_misID(cut+"&&"+*it,true,TString::Format("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/misIDShapes/misID_KKmumu_bin_%i.eps",counter),"m_{KK#mu#mu}(KK#pi#pi)",TString::Format("%.0fMeV<m_{#mu#mu}(#pi#pi)<%.0fMeV",rangesKK_low[counter],rangesKK_high[counter]));
    }

    //if(counter==1)
    if(kind=="D2pipimumu"){
      myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_new.root");
      myFitter1D->fit_HHpipi_misID(cut+"&&"+*it,true,TString::Format("/work/mitzel/D2hhmumu/dev/D2pipimumu/img/misIDShapes/misID_pipimumu_bin_%i.eps",counter),"m_{#pi#pi#mu#mu}(#pi#pi#pi#pi)",TString::Format("%.0fMeV<m_{#mu#mu}(#pi#pi)<%.0fMeV",rangespipi_low[counter],rangespipi_high[counter]));
    }

    ++counter;			 
  }
  
  
  if(kind=="D2KKmumu"){
    myFitter1D = new D2hhmumuFitter1D();
    myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_new.root");
    myFitter1D->fit_HHpipi_misID(cut+"&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875",true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/misIDShapes/misID_Kpimumu_KKTrained.eps","m_{K#pi#mu#mu}(K#pi#pi#pi)");
  }
  
  if(kind=="D2pipimumu"){
    myFitter1D = new D2hhmumuFitter1D();
    myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_new.root");
    myFitter1D->fit_HHpipi_misID(cut+"&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875",true,"/work/mitzel/D2hhmumu/dev/D2pipimumu/img/misIDShapes/misID_Kpimumu_pipiTrained.eps","m_{K#pi#mu#mu}(K#pi#pi#pi)");
  }
  
} 
    


void D2hhmumuFitter_Applications::drawMCSignalShapes(TString cut){

  D2hhmumuFitter1D* myFitter1D;
  int counter=0;
  
  
  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

    myFitter1D = new D2hhmumuFitter1D();

    if(kind=="D2KKmumu"){
      myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT.root");
      myFitter1D->fit_MC(cut+"&&"+*it,true,TString::Format("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/signalShapes/MC_KKmumu_bin_%i.eps",counter),"m(KK#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangesKK_low[counter],rangesKK_high[counter]));
    }

    if(kind=="D2pipimumu"){
      myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root");
      myFitter1D->fit_MC(cut+"&&"+*it,true,TString::Format("/work/mitzel/D2hhmumu/dev/D2pipimumu/img/signalShapes/MC_pipimumu_bin_%i.eps",counter),"m(#pi#pi#mu#mu)",TString::Format("%.0fMeV<m(#mu#mu)<%.0fMeV",rangespipi_low[counter],rangespipi_high[counter]));
    }

    ++counter;			 
  }
  

  if(kind=="D2KKmumu"){
    myFitter1D = new D2hhmumuFitter1D();
    myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root");
    myFitter1D->fit_normalization_MC(cut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/signalShapes/MC_norm_KKTrained.eps");
  }

  if(kind=="D2pipimumu"){
    myFitter1D = new D2hhmumuFitter1D();
    myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");
    myFitter1D->fit_normalization_MC(cut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/signalShapes/MC_norm_pipiTrained.eps");
  }

} 
    
      
    
void D2hhmumuFitter_Applications::studyResolutionScale(TString cut){

  D2hhmumuFitter1D* myFitter1D;
  int counter=0;

  /*
  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

    myFitter1D = new D2hhmumuFitter1D();

    if(kind=="D2KKmumu"){
      myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKpipi_D2KKmumuBDT.root");
      myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT.root");
      myFitter1D->fit_MC(cut+"&&"+*it,true,TString::Format("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyResolutionScale/MC_D2KKpipi_bin_%i.eps",counter));
      myFitter1D->fit_HHpipi(cut+"&&"+*it,TString::Format("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyResolutionScale/Data_D2KKpipi_bin_%i.eps",counter));
    }

    if(kind=="D2pipimumu"){
      myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipipipi_D2pipimumuBDT.root");
      myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT.root");
      myFitter1D->fit_MC(cut+"&&"+*it,true,TString::Format("/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyResolutionScale/MC_D2pipipipi_bin_%i.eps",counter));
      myFitter1D->fit_HHpipi(cut+"&&"+*it,TString::Format("/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyResolutionScale/Data_D2pipipipi_bin_%i.eps",counter));
    }

    ++counter;
  }
  */

  myFitter1D = new D2hhmumuFitter1D();

  if(kind=="D2KKmumu"){
    myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2KKmumuBDT.root");
    myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT.root");
    myFitter1D->fit_MC(cut+"&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875",true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyResolutionScale/MC_D2Kpipipi.eps");
    myFitter1D->fit_HHpipi(cut+"&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyResolutionScale/Data_D2Kpipipi.eps");  
  }

  if(kind=="D2pipimumu"){
    myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2pipimumuBDT.root");
    myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT.root");
    myFitter1D->fit_MC(cut+"&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875",true,"/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyResolutionScale/MC_D2Kpipipi.eps");
    myFitter1D->fit_HHpipi(cut+"&&D_DiMuon_Mass>675&&D_DiMuon_Mass<875","/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyResolutionScale/Data_D2Kpipipi.eps");  
  }

  myFitter1D = new D2hhmumuFitter1D();

  if(kind=="D2pipimumu"){
    myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipipipi_D2pipimumuBDT.root");
    myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT.root");
    myFitter1D->fit_MC(cut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyResolutionScale/MC_D2pipipipi.eps");
    myFitter1D->fit_HHpipi(cut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyResolutionScale/Data_D2pipipipi.eps");
  }

  if(kind=="D2KKmumu"){
    myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKpipi_D2KKmumuBDT.root");
    myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT.root");
    myFitter1D->fit_MC(cut,true,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyResolutionScale/MC_D2KKpipi.eps");
    myFitter1D->fit_HHpipi(cut,"/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyResolutionScale/Data_D2KKpipi.eps");
  }


}



void D2hhmumuFitter_Applications::studyTriggerEfficiencyNormalizationMode(){

  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/noPreselectionCuts/D2Kpimumu_D2KKmumuBDT_noCuts.root");
  myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpipipi_D2KKmumuBDT.root");


  TString preselection  = "D_DiMuon_Mass>675&&D_DiMuon_Mass<875&&mu0_ProbNNghost<0.5&&mu1_ProbNNghost<0.5&&h0_ProbNNghost<0.5&&h1_ProbNNghost<0.5&&Slowpi_ProbNNghost<0.5&&mu0_MuonNShared==0&&mu1_MuonNShared==0&&h0_ProbNNk>0.2&&h1_ProbNNpi>0.2"; 
  TString offlineSelection = "&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.";

  // only very lose preselection
  
  myFitter1D.fit_normalization_MC(preselection,true,"test1.eps");
  myFitter1D.fit_Kpipipi_misID(preselection+"&&mu1_ProbNNmu>0.",true,"test2.eps");
  myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>0.&&mu0_ProbNNmu>0.","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection.eps");
/*
  //BDT and PID 
  myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel.eps");
  //BDT and PID && L0_Muon_TOS
  myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0MuonDecision_TOS","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0muon_TOS.eps");
  //BDT and PID && L0_global_TIS
  myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0Global_TIS","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0global_TIS.eps");
  //BDT and PID && L0_global_TIS && L0_Muon_TOS 
  myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0Global_TIS&&D_L0MuonDecision_TOS","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0global_TIS_and_L0Muon_TOS.eps");
  //BDT and PID && L0_Muon_TOS && HLT1                                                                                                                                                   
  myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0MuonDecision_TOS&&(mu0_Hlt1TrackMuonDecision_TOS==1 || mu1_Hlt1TrackMuonDecision_TOS==1 || D_Hlt1TrackAllL0Decision_TOS==1)","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0muonTOS_and_Hlt1_TOS.eps");
    //BDT and PID && L0_Muon_TOS && HLT1 && HLT2
 myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0MuonDecision_TOS&&(mu0_Hlt1TrackMuonDecision_TOS || mu1_Hlt1TrackMuonDecision_TOS || D_Hlt1TrackAllL0Decision_TOS)&&D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0muonTOS_and_Hlt1_TOS_and_Hlt2_TOS.eps");
  //BDT and PID && L0_Muon_TIS && HLT1 && HLT2 
 myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0Global_TIS&&(mu0_Hlt1TrackMuonDecision_TOS || mu1_Hlt1TrackMuonDecision_TOS || D_Hlt1TrackAllL0Decision_TOS)&&D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0globalTIS_and_Hlt1_TOS_and_Hlt2_TOS.eps");
    //BDT and PID && L0_Muon_TOS && L0_global_TIS && HLT1 && HLT2 
 myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0Global_TIS&&D_L0MuonDecision_TOS&&(mu0_Hlt1TrackMuonDecision_TOS || mu1_Hlt1TrackMuonDecision_TOS || D_Hlt1TrackAllL0Decision_TOS)&&D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0globalTISMuonTOS_and_Hlt1_TOS_and_Hlt2_TOS.eps");

 //BDT and PID && L0_Muon_TOS && L0_muon_TIS 
 myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0MuonDecision_TIS&&D_L0MuonDecision_TOS","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0MuonTISTOS.eps");
 //BDT and PID && L0_Muon_TOS && L0_muon_TIS 
 myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0MuonDecision_TIS","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0MuonTIS.eps");
*/

 //hadron TIS
 myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0HadronDecision_TIS&&D_L0MuonDecision_TOS","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0HadronTIS_L0MuonTOS.eps");
 //BDT and PID && L0_Muon_TOS && L0_muon_TIS 
 myFitter1D.fit_normalization_Data(preselection+"&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.&&D_L0HadronDecision_TIS","/work/mitzel/D2hhmumu/dev/D2KKmumu/img/studyTriggerEfficiencyNormMode/preselection_and_offlineSel_and_L0HadronTIS.eps");

  
}

void D2hhmumuFitter_Applications::studyTISEfficiencyNormalizationModeBinned(){

  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/noPreselectionCuts/D2Kpimumu_D2KKmumuBDT_noCuts.root");
  myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpipipi_D2KKmumuBDT.root");

  TString preselection  = "D_DiMuon_Mass>675&&D_DiMuon_Mass<875&&mu0_ProbNNghost<0.5&&mu1_ProbNNghost<0.5&&h0_ProbNNghost<0.5&&h1_ProbNNghost<0.5&&Slowpi_ProbNNghost<0.5&&mu0_MuonNShared==0&&mu1_MuonNShared==0&&h0_ProbNNk>0.2&&h1_ProbNNpi>0.2"; 

  TString offlineSelection = "&&mu1_ProbNNmu>.5&&mu0_ProbNNmu>0.5&&BDT>0.";
  TString L0TIS = "&&D_L0Global_TIS";
  TString L0TISTOS = "&&D_L0Global_TIS&&D_L0MuonDecision_TOS";
 
  double nTIS=0;
  double nTISTOS=0;
  
  ///// int nX=3;
  /////Float_t xRanges[4] = {0,5000,8000,25000};
  //std::vector<double> yRanges {0,50e3,100e3,150e3,2000e3};                                                                                                                         
  /////std::vector<double> yRanges {0,50e3,90e3,250e3};
 
  int nX=2;
  Float_t xRanges[3] = {0,8000,25000};
  //std::vector<double> yRanges {0,50e3,100e3,150e3,2000e3};                                                                                                                         
  std::vector<double> yRanges {0,90e3,250e3};
  

  int nY= yRanges.size() ;
  nY-=1;

  TH1D* temp1;
  TH1D* temp2;

  std::vector<TH1D*> histos_TISTOS, histos_TIS;

  for(int i = 0;i<nY;++i) {

    temp1 = new TH1D(TString::Format("nTISTOS_bin_%i",i),TString::Format("nTISTOS_bin_%i",i),nX,xRanges);
    temp2 = new TH1D(TString::Format("nTIS_bin_%i",i),TString::Format("nTIS_bin_%i",i),nX,xRanges);
    histos_TISTOS.push_back(temp1);
    histos_TIS.push_back(temp2);
  }

  myFitter1D.fit_normalization_MC(preselection+offlineSelection,true,"test1.eps");
  myFitter1D.fit_Kpipipi_misID(preselection+"&&mu1_ProbNNmu>0.&&BDT>0",true,"test2.eps");
    

  for(int i = 0;i<nY;++i) {
    
    for(int j = 1;j<nX+1;++j){

      TString totalCut_TIS=preselection+offlineSelection+TString::Format("&&D_PT>%f&&D_PT<%f",xRanges[j-1],xRanges[j])+TString::Format("&&D_PZ>%f&&D_PZ<%f",yRanges[i],yRanges[i+1])+L0TIS;
      std::cout<<"CUT FOR TIS "<<i<<"  "<<j<<"  "<<totalCut_TIS<<std::endl;
      nTIS = myFitter1D.fit_normalization_Data(totalCut_TIS,"test3.eps");
      histos_TIS[i]->SetBinContent(j,nTIS);
      TString totalCut_TISTOS=preselection+offlineSelection+TString::Format("&&D_PT>%f&&D_PT<%f",xRanges[j-1],xRanges[j])+TString::Format("&&D_PZ>%f&&D_PZ<%f",yRanges[i],yRanges[i+1])+L0TISTOS;
      std::cout<<"CUT FOR TISTOS "<<i<<"  "<<j<<"  "<<totalCut_TISTOS<<std::endl;
      nTISTOS= myFitter1D.fit_normalization_Data(totalCut_TISTOS,"test4.eps");
      histos_TISTOS[i]->SetBinContent(j,nTISTOS);;
      }
  }

    std::vector<TH1D*> Eff_TISTOS;
    std::vector<double> errors;
    TH1D* Eff_TISTOS_temp;

    for(int i = 0;i<nY;++i) {


      //calculating the errors for tistos                                                                                                                                                     
      for(int j = 1;j<nX+1;++j) {
	double k = histos_TISTOS[i]->GetBinContent(j);
	double n = histos_TIS[i]->GetBinContent(j);
	double dE = 1/n * TMath::Sqrt(k *(1 - (k/n) ) );
	errors.push_back(dE);
	std::cout<<"TISTOS eff errors: " <<" j: "<<j<<" dE: "<<dE<<" n: "<<n<<" k: "<<k<< "k/n: "<<k/n<< std::endl;
      }

      Eff_TISTOS_temp = (TH1D*)histos_TISTOS[i]->Clone();
      Eff_TISTOS_temp->Divide(histos_TIS[i]);

      //setting the error                                                                                                                                                                     
      for(int j = 1;j<nX+1;++j) {
	Eff_TISTOS_temp->SetBinError(j, errors[j]);
      }

      Eff_TISTOS.push_back(Eff_TISTOS_temp);

    }

    TCanvas *a = new TCanvas("a","a");
    a->Divide(3,3);
    for(int i = 0;i<nY;++i) {
      a->cd(i+1);
      Eff_TISTOS[i]->GetYaxis()->SetRangeUser(0,1);
      Eff_TISTOS[i]->Draw("E");
    }

    a->Print("test.eps");

}
