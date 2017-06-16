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
    pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT_noMultCand.root";
    pathToSidebandData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/sideband_D2KKmumu_BDT.root";
    pathToKpipipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_Randomized.root";
    pathToKKpipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT.root";
    pathToKpimumuMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root";
    pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT_noMultCand.root";
    pathToNormMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root";
  }

  if(m_kind== "D2pipimumu"){
    pathToSignalMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root";
    pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT_noMultCand.root";
    pathToSidebandData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/sideband_D2pipimumu_BDT.root";
    pathToKpipipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_Randomized.root";
    pathToKKpipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_Randomized.root";
    pathToKpimumuMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root";
    pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT_noMultCand.root";
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
    //if(q2Bin==2 || q2Bin==3) continue;
    
    TString fOut= targetFolder+"ModelConfigs/"+ kind + "_" + TString::Format("bin_%i",q2Bin)+".root"; 

    myFitter1D->setPathToSignalMC(pathToSignalMC);
    myFitter1D->setPathToNormMC(pathToNormMC);
    myFitter1D->setPathToKpipipiData(pathToKpipipiData);
    myFitter1D->setPathToHHpipiData(pathToKKpipiData);
    myFitter1D->setPathToNormData(pathToNormData);
    myFitter1D->setPathToSignalData(pathToSignalData);

    std::cout<<kind<<"  "<<dataCut<<"  "<<misIDCut<<"   "<<q2Bin<<"   "<<fOut<<std::endl;
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


void D2hhmumuFitter_Applications::draw_misID_shapes_singlePlot() {
 
  RooRealVar D0_M("Dst_DTF_D0_M", "m_{K#pi#mu#mu}(K#pi#pi#pi)", 1760., 1930.,"MeV/c^{2}");

  /*
  RooRealVar D0_M_xi_bkg_piplus("D0_M_xi_bkg_piplus","D0_M_xi_bkg_piplus",1835.72);
  RooRealVar  D0_M_lambda_bkg_piplus("D0_M_lambda_bkg_piplus","D0_M_lambda_bkg_piplus",28.30);
  RooRealVar  D0_M_gamma_bkg_piplus("D0_M_gamma_bkg_piplus","D0_M_gamma_bkg_piplus",-0.82);
  RooRealVar  D0_M_delta_bkg_piplus("D0_M_delta_bkg_piplus","D0_M_delta_bkg_piplus",1.01);

    RooRealVar D0_M_xi_bkg_piminus("D0_M_xi_bkg_piminus","D0_M_xi_bkg_piminus",1838.20);
  RooRealVar  D0_M_lambda_bkg_piminus("D0_M_lambda_bkg_piminus","D0_M_lambda_bkg_piminus",28.57);
  RooRealVar  D0_M_gamma_bkg_piminus("D0_M_gamma_bkg_piminus","D0_M_gamma_bkg_piminus",-0.63);
  RooRealVar  D0_M_delta_bkg_piminus("D0_M_delta_bkg_piminus","D0_M_delta_bkg_piminus",1.09);

  RooRealVar D0_M_xi_bkg_noPID("D0_M_xi_bkg_noPID","D0_M_xi_bkg_noPID",1842.01);
  RooRealVar  D0_M_lambda_bkg_noPID("D0_M_lambda_bkg_noPID","D0_M_lambda_bkg_noPID",14.48);
  RooRealVar  D0_M_gamma_bkg_noPID("D0_M_gamma_bkg_noPID","D0_M_gamma_bkg_noPID",-0.89);
  RooRealVar  D0_M_delta_bkg_noPID("D0_M_delta_bkg_noPID","D0_M_delta_bkg_noPID",0.69);

  RooRealVar D0_M_xi_bkg_lowPT("D0_M_xi_bkg_lowPT","D0_M_xi_bkg_lowPT",1836.44);
  RooRealVar  D0_M_lambda_bkg_lowPT("D0_M_lambda_bkg_lowPT","D0_M_lambda_bkg_lowPT",34.73);
  RooRealVar  D0_M_gamma_bkg_lowPT("D0_M_gamma_bkg_lowPT","D0_M_gamma_bkg_lowPT",-0.64);
  RooRealVar  D0_M_delta_bkg_lowPT("D0_M_delta_bkg_lowPT","D0_M_delta_bkg_lowPT",1.17);

  RooRealVar D0_M_xi_bkg_highPT("D0_M_xi_bkg_highPT","D0_M_xi_bkg_highPT",1837.19);
  RooRealVar  D0_M_lambda_bkg_highPT("D0_M_lambda_bkg_highPT","D0_M_lambda_bkg_highPT",25.85);
  RooRealVar  D0_M_gamma_bkg_highPT("D0_M_gamma_bkg_highPT","D0_M_gamma_bkg_highPT",-0.78);
  RooRealVar  D0_M_delta_bkg_highPT("D0_M_delta_bkg_highPT","D0_M_delta_bkg_highPT",0.99);
  */


  //non randomized shape
  RooRealVar D0_M_xi_bkg_piplus("D0_M_xi_bkg_piplus","D0_M_xi_bkg_piplus",1834.28);
  RooRealVar  D0_M_lambda_bkg_piplus("D0_M_lambda_bkg_piplus","D0_M_lambda_bkg_piplus",31.35);
  RooRealVar  D0_M_gamma_bkg_piplus("D0_M_gamma_bkg_piplus","D0_M_gamma_bkg_piplus",-0.85);
  RooRealVar  D0_M_delta_bkg_piplus("D0_M_delta_bkg_piplus","D0_M_delta_bkg_piplus",1.04);

  RooRealVar D0_M_xi_bkg_piminus("D0_M_xi_bkg_piminus","D0_M_xi_bkg_piminus",1839.20);
  RooRealVar  D0_M_lambda_bkg_piminus("D0_M_lambda_bkg_piminus","D0_M_lambda_bkg_piminus",28.10);
  RooRealVar  D0_M_gamma_bkg_piminus("D0_M_gamma_bkg_piminus","D0_M_gamma_bkg_piminus",-0.63);
  RooRealVar  D0_M_delta_bkg_piminus("D0_M_delta_bkg_piminus","D0_M_delta_bkg_piminus",1.12);

  RooRealVar D0_M_xi_bkg_noPID("D0_M_xi_bkg_noPID","D0_M_xi_bkg_noPID",1841.87);
  RooRealVar  D0_M_lambda_bkg_noPID("D0_M_lambda_bkg_noPID","D0_M_lambda_bkg_noPID",15.12);
  RooRealVar  D0_M_gamma_bkg_noPID("D0_M_gamma_bkg_noPID","D0_M_gamma_bkg_noPID",-0.88);
  RooRealVar  D0_M_delta_bkg_noPID("D0_M_delta_bkg_noPID","D0_M_delta_bkg_noPID",0.73);

  RooRealVar D0_M_xi_bkg_lowPT("D0_M_xi_bkg_lowPT","D0_M_xi_bkg_lowPT",1837.34);
  RooRealVar  D0_M_lambda_bkg_lowPT("D0_M_lambda_bkg_lowPT","D0_M_lambda_bkg_lowPT",33.88);
  RooRealVar  D0_M_gamma_bkg_lowPT("D0_M_gamma_bkg_lowPT","D0_M_gamma_bkg_lowPT",-0.67);
  RooRealVar  D0_M_delta_bkg_lowPT("D0_M_delta_bkg_lowPT","D0_M_delta_bkg_lowPT",1.17);

  RooRealVar D0_M_xi_bkg_highPT("D0_M_xi_bkg_highPT","D0_M_xi_bkg_highPT",1836.40);
  RooRealVar  D0_M_lambda_bkg_highPT("D0_M_lambda_bkg_highPT","D0_M_lambda_bkg_highPT",28.63);
  RooRealVar  D0_M_gamma_bkg_highPT("D0_M_gamma_bkg_highPT","D0_M_gamma_bkg_highPT",-0.79);
  RooRealVar  D0_M_delta_bkg_highPT("D0_M_delta_bkg_highPT","D0_M_delta_bkg_highPT",1.04);


  RooRealVar D0_M_xi_bkg_Hlt1TOS("D0_M_xi_bkg_Hlt1TOS","D0_M_xi_bkg_Hlt1TOS",1.8393e+03);
  RooRealVar  D0_M_lambda_bkg_Hlt1TOS("D0_M_lambda_bkg_Hlt1TOS","D0_M_lambda_bkg_Hlt1TOS",2.7318e+01);
  RooRealVar  D0_M_gamma_bkg_Hlt1TOS("D0_M_gamma_bkg_Hlt1TOS","D0_M_gamma_bkg_Hlt1TOS",-6.5558e-01);
  RooRealVar  D0_M_delta_bkg_Hlt1TOS("D0_M_delta_bkg_Hlt1TOS","D0_M_delta_bkg_Hlt1TOS",1.1003e+00);

  RooRealVar D0_M_xi_bkg_L0TOS("D0_M_xi_bkg_L0TOS","D0_M_xi_bkg_L0TOS",1.8405e+03);
  RooRealVar  D0_M_lambda_bkg_L0TOS("D0_M_lambda_bkg_L0TOS","D0_M_lambda_bkg_L0TOS",2.5788e+01);
  RooRealVar  D0_M_gamma_bkg_L0TOS("D0_M_gamma_bkg_L0TOS","D0_M_gamma_bkg_L0TOS",-6.0702e-01);
  RooRealVar  D0_M_delta_bkg_L0TOS("D0_M_delta_bkg_L0TOS","D0_M_delta_bkg_L0TOS",1.0987e+00);

  
  RooJohnsonSU Signal_piplus("Signal_piplus", "Signal D^{0} JSU",D0_M,D0_M_xi_bkg_piplus,D0_M_lambda_bkg_piplus,D0_M_gamma_bkg_piplus,D0_M_delta_bkg_piplus);
  RooJohnsonSU Signal_piminus("Signal_piminus", "Signal D^{0} JSU",D0_M,D0_M_xi_bkg_piminus,D0_M_lambda_bkg_piminus,D0_M_gamma_bkg_piminus,D0_M_delta_bkg_piminus);
  RooJohnsonSU Signal_noPID("Signal_noPID", "Signal D^{0} JSU",D0_M,D0_M_xi_bkg_noPID,D0_M_lambda_bkg_noPID,D0_M_gamma_bkg_noPID,D0_M_delta_bkg_noPID);
 
  RooJohnsonSU Signal_lowPT("Signal_lowPT", "Signal D^{0} JSU",D0_M,D0_M_xi_bkg_lowPT,D0_M_lambda_bkg_lowPT,D0_M_gamma_bkg_lowPT,D0_M_delta_bkg_lowPT);
  RooJohnsonSU Signal_highPT("Signal_highPT", "Signal D^{0} JSU",D0_M,D0_M_xi_bkg_highPT,D0_M_lambda_bkg_highPT,D0_M_gamma_bkg_highPT,D0_M_delta_bkg_highPT);

  RooJohnsonSU Signal_Hlt1TOS("Signal_Hlt1TOS", "Signal D^{0} JSU",D0_M,D0_M_xi_bkg_Hlt1TOS,D0_M_lambda_bkg_Hlt1TOS,D0_M_gamma_bkg_Hlt1TOS,D0_M_delta_bkg_Hlt1TOS);
  RooJohnsonSU Signal_L0TOS("Signal_L0TOS", "Signal D^{0} JSU",D0_M,D0_M_xi_bkg_L0TOS,D0_M_lambda_bkg_L0TOS,D0_M_gamma_bkg_L0TOS,D0_M_delta_bkg_L0TOS);



  TCanvas* c1= new TCanvas("");
  //c1->Divide(2,2);                                                                                                                                                              
  //CreateSubPad(c1,0.25);

  c1->cd(1);

  RooPlot* frame_m= D0_M.frame();                                                                                                                                                 

  frame_m->SetTitle("");                                                                                                                                                           Signal_piplus.plotOn(frame_m,LineColor(kRed),LineWidth(2));
  Signal_noPID.plotOn(frame_m,LineColor(kBlue),LineWidth(2));

  Signal_Hlt1TOS.plotOn(frame_m,LineColor(kBlack),LineWidth(2));
  Signal_L0TOS.plotOn(frame_m,LineColor(11),LineWidth(2));

  Signal_piminus.plotOn(frame_m,LineColor(kGreen),LineWidth(2),LineStyle(kDashed)); 

  Signal_lowPT.plotOn(frame_m,LineColor(kCyan),LineWidth(2));
  Signal_highPT.plotOn(frame_m,LineColor(kViolet),LineWidth(2));

  frame_m->GetYaxis()->SetTitle("PDF");

  frame_m->Draw();

  TLatex *leg2 = new TLatex();
  leg2->SetNDC();
  leg2->SetTextAlign(13);
  leg2->SetTextFont(63);
  leg2->SetTextSizePixels(18);
  leg2->DrawLatex(0.2,0.635+0.25,"#color[3]{PID on #pi^{-}(default)}");
  leg2->DrawLatex(0.2,0.595+0.25,"#color[2]{PID on #pi^{+}}");
  leg2->DrawLatex(0.2,0.555+0.25,"#color[4]{no PID}");
  leg2->DrawLatex(0.2,0.515+0.25,"#color[7]{PID high p_{T} #pi}");
  leg2->DrawLatex(0.2,0.475+0.25,"#color[6]{PID low p_{T} #pi}");
  leg2->DrawLatex(0.2,0.435+0.25,"#color[11]{PID #pi^{-} and L0 Trigger} ");
  leg2->DrawLatex(0.2,0.395+0.25,"#color[1]{PID #pi^{-} and Hlt1 Trigger}");

  c1->Print("DifferentMisIDShapes.eps");


}

void D2hhmumuFitter_Applications::compare_misID_shapes(TString dataCut="BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5"){

  TString PathToFolder = "/work/mitzel/D2hhmumu/dev/"+kind+"/img/comparison_misID_shapes/";

  //the non randomized shapes
  //TString PathToFolder = "/work/mitzel/D2hhmumu/dev/"+kind+"/img/comparison_misID_shapes/nonRandomizedSample/";
  //pathToKpipipiData="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT.root";
  //pathToKKpipiData="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT.root";

  TString nameExtension="";

  bool array[2] = {true,false};

  for(int i=0;i<2;++i){
    
    if(!array[i]) { 
      PathToFolder+="freeCombinatorialBkgShape/";
      nameExtension+="freeCombBkg";
    }

    
    D2hhmumuFitter1D myFitter1D("comparison_misID_shapes/singleMuonID_piminus"+nameExtension+".txt");
    //D2hhmumuFitter1D myFitter1D("comparison_misID_shapes/nonRandomizedSample/singleMuonID_piminus"+nameExtension+".txt");
    myFitter1D.setPathToNormMC(pathToNormMC);
    myFitter1D.setPathToKpipipiData(pathToKpipipiData);
    myFitter1D.setPathToNormData(pathToNormData);
    myFitter1D.setPathToSignalData(pathToSignalData);
    myFitter1D.setPathToSignalMC(pathToSignalMC);
    myFitter1D.setPathToHHpipiData(pathToKKpipiData);
    //fix MC signal shape
    myFitter1D.fit_MC(dataCut+"&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"MCSignalShape.eps");
    myFitter1D.fit_normalization_MC(dataCut,true,PathToFolder+"MCShape.eps");
    //single ID to pi-
    myFitter1D.fit_HHpipi_misID("BDT>0.4&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"pipipipi_singleMuonID_piminus.eps","","",array[i]);
    myFitter1D.fit_Kpipipi_misID("BDT>0.4&&mu1_ProbNNmu>0.5",true,PathToFolder+"Kpipipi_singleMuonID_piminus.eps",array[i]);
    myFitter1D.fit_normalization_Data(dataCut,PathToFolder+"norm_singleMuonID_piminus.eps");
    myFitter1D.fit_unblinded_Data("D2pipimumu",dataCut,"D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"signal_singleMuonID_piminus.eps");

    D2hhmumuFitter1D myFitter1D_2("comparison_misID_shapes/singleMuonID_piplus"+nameExtension+".txt");
    //D2hhmumuFitter1D myFitter1D_2("comparison_misID_shapes/nonRandomizedSample/singleMuonID_piplus"+nameExtension+".txt");
    myFitter1D_2.setPathToNormMC(pathToNormMC);
    myFitter1D_2.setPathToKpipipiData(pathToKpipipiData);
    myFitter1D_2.setPathToNormData(pathToNormData);
    myFitter1D_2.setPathToSignalData(pathToSignalData);
    myFitter1D_2.setPathToSignalMC(pathToSignalMC);
    myFitter1D_2.setPathToHHpipiData(pathToKKpipiData);
    //fix MC signal shape
    myFitter1D_2.fit_normalization_MC(dataCut,true,PathToFolder+"MCShape.eps");
    myFitter1D_2.fit_MC(dataCut+"&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"MCSignalShape.eps");
    //single ID to pi+
    myFitter1D_2.fit_HHpipi_misID("BDT>0.4&&mu0_ProbNNmu>0.5&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"pipipipi_singleMuonID_piplus.eps",array[i]);
    myFitter1D_2.fit_Kpipipi_misID("BDT>0.4&&mu0_ProbNNmu>0.5",true,PathToFolder+"Kpipipi_singleMuonID_piplus.eps",array[i]);
    myFitter1D_2.fit_normalization_Data(dataCut,PathToFolder+"norm_singleMuonID_piplus.eps");
    myFitter1D_2.fit_unblinded_Data("D2pipimumu",dataCut,"D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"signal_singleMuonID_piplus.eps");

    D2hhmumuFitter1D myFitter1D_3("comparison_misID_shapes/lowPTMuonID.txt");
    //D2hhmumuFitter1D myFitter1D_3("comparison_misID_shapes/nonRandomizedSample/lowPTMuonID.txt");
    myFitter1D_3.setPathToNormMC(pathToNormMC);
    myFitter1D_3.setPathToKpipipiData(pathToKpipipiData);
    myFitter1D_3.setPathToNormData(pathToNormData);
    myFitter1D_3.setPathToSignalData(pathToSignalData);
    myFitter1D_3.setPathToSignalMC(pathToSignalMC);
    myFitter1D_3.setPathToHHpipiData(pathToKKpipiData);
    //fix MC signal shape
    myFitter1D_3.fit_normalization_MC(dataCut,true,PathToFolder+"MCShape.eps");
    myFitter1D_3.fit_MC(dataCut+"&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"MCSignalShape.eps");
    //low pT pi
    myFitter1D_3.fit_HHpipi_misID("BDT>0.4 && ( ( (mu1_PT < mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT > mu0_PT) && mu0_ProbNNmu>0.5))&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"pipipipi_lowPTMuonID.eps",array[i]);
    myFitter1D_3.fit_Kpipipi_misID("BDT>0.4 && ( ( (mu1_PT < mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT > mu0_PT) && mu0_ProbNNmu>0.5))",true,PathToFolder+"Kpipipi_lowPTMuonID.eps",array[i]);
    myFitter1D_3.fit_normalization_Data(dataCut,PathToFolder+"norm_lowPTMuonID.eps");
    myFitter1D_3.fit_unblinded_Data("D2pipimumu",dataCut,"D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"signal_lowPTMuonID.eps");


    D2hhmumuFitter1D myFitter1D_4("comparison_misID_shapes/highPTMuonID"+nameExtension+".txt");
    //D2hhmumuFitter1D myFitter1D_4("comparison_misID_shapes/nonRandomizedSample/highPTMuonID"+nameExtension+".txt");
    myFitter1D_4.setPathToNormMC(pathToNormMC);
    myFitter1D_4.setPathToKpipipiData(pathToKpipipiData);
    myFitter1D_4.setPathToNormData(pathToNormData);
    myFitter1D_4.setPathToSignalData(pathToSignalData);
    myFitter1D_4.setPathToSignalMC(pathToSignalMC);
    myFitter1D_4.setPathToHHpipiData(pathToKKpipiData);
    //fix MC signal shape
    myFitter1D_4.fit_normalization_MC(dataCut,true,PathToFolder+"MCShape.eps");
    myFitter1D_4.fit_MC(dataCut+"&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"MCSignalShape.eps");
    //low pT pi
    myFitter1D_4.fit_HHpipi_misID("BDT>0.4 && ( ( (mu1_PT > mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT < mu0_PT) && mu0_ProbNNmu>0.5))&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"pipipipi_highPTMuonID.eps",array[i]); 
    myFitter1D_4.fit_Kpipipi_misID("BDT>0.4 && ( ( (mu1_PT > mu0_PT) && mu1_ProbNNmu>0.5) || ((mu1_PT < mu0_PT) && mu0_ProbNNmu>0.5))",true,PathToFolder+"Kpipipi_lowPTMuonID.eps",array[i]);
    myFitter1D_4.fit_normalization_Data(dataCut,PathToFolder+"norm_highPTMuonID.eps");
    myFitter1D_4.fit_unblinded_Data("D2pipimumu",dataCut,"D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"signal_highPTMuonID.eps");
  

    D2hhmumuFitter1D myFitter1D_5("comparison_misID_shapes/noMuonID"+nameExtension+".txt");
    //D2hhmumuFitter1D myFitter1D_5("comparison_misID_shapes/nonRandomizedSample/noMuonID"+nameExtension+".txt");
    myFitter1D_5.setPathToNormMC(pathToNormMC);
    myFitter1D_5.setPathToKpipipiData(pathToKpipipiData);
    myFitter1D_5.setPathToNormData(pathToNormData);
    myFitter1D_5.setPathToSignalData(pathToSignalData);
    myFitter1D_5.setPathToSignalMC(pathToSignalMC);
    myFitter1D_5.setPathToHHpipiData(pathToKKpipiData);
    //fix MC signal shape
    myFitter1D_5.fit_normalization_MC(dataCut,true,PathToFolder+"MCShape.eps");
    myFitter1D_5.fit_MC(dataCut+"&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"MCSignalShape.eps");
    //single ID to pi-
    myFitter1D_5.fit_HHpipi_misID("BDT>0.4&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"pipipipi_noMounID.eps",array[i]);
    myFitter1D_5.fit_Kpipipi_misID("BDT>0.4&&nTracks%15==0",true,PathToFolder+"Kpipipi_noMuonID.eps",array[i]);
    myFitter1D_5.fit_normalization_Data(dataCut,PathToFolder+"norm_noMuonID.eps");
    myFitter1D_5.fit_unblinded_Data("D2pipimumu",dataCut,"D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"signal_noMuonID.eps");
  }  


}  

void D2hhmumuFitter_Applications::ExtractExpectedLimit(){
  
  //StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 0, 3, true, 20,1e-9,5e-8,20000);
  //good range for KKmumu //StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 2, 3, true, 20,1e-9,7e-8);                                                          
 
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);                                                                                                                         
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1);                                                                                                                             
  ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-12); 

  TString fIn;  
  
  //StandardHypoTestInvDemo("/work/mitzel/D2hhmumu/dev/D2pipimumu/ModelConfigs/D2pipimumu_bin_1.root", "m_ws", "ModelConfig", "B_only_model", "data", 2, 3, true, 10000,1e-10,1e-7);
  //StandardHypoTestInvDemo("/work/mitzel/D2hhmumu/dev/D2pipimumu/ModelConfigs/D2pipimumu_bin_4.root", "m_ws", "ModelConfig", "B_only_model", "data", 2, 3, true, 10000,1e-10,1e-7);  
  //StandardHypoTestInvDemo("/work/mitzel/D2hhmumu/dev/D2pipimumu/ModelConfigs/D2pipimumu_bin_1.root", "m_ws", "ModelConfig", "B_only_model", "data", 0, 3, true, 40,1e-10,5e-8,500,false);
  //StandardHypoTestInvDemo("/work/mitzel/D2hhmumu/dev/D2pipimumu/ModelConfigs/D2pipimumu_bin_4.root", "m_ws", "ModelConfig", "B_only_model", "data", 0, 3, true, 40,1e-10,5e-8,500,false);
  //StandardHypoTestInvDemo("testConfig.root", "m_ws", "ModelConfig", "B_only_model", "data", 2, 3, true, 5000,1e-10,1e-7);
  
  for(int q2Bin=0; q2Bin<q2Ranges.size();++q2Bin) {

	fIn = targetFolder+"ModelConfigs/"+ kind + "_" + TString::Format("bin_%i",q2Bin)+".root";

	if(kind=="D2KKmumu" ){

	//StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 2, 3, true, 30000,1e-10,1e-7);                                                                   
	//StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 0, 3, true, 500,1e-10,8e-8,5000,false);                                                                
 	
	//eat the moment we agreed on taking the LEP test statistics
	  
	  if(q2Bin==0)StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 0, 0, true, 15,0.8e-8,7e-8,5000,false);                                            
	  if(q2Bin==1) continue;
	  if(q2Bin==2)StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 0, 0, true, 15,5e-8,0.17e-6,5000,false);                                                    
	}

      if(kind=="D2pipimumu" ){

	  if(q2Bin==0)StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 0, 0, true, 15,0.8e-8,1.5e-7,5000,false);                                            
	  if(q2Bin==1)StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 0, 0, true, 15,0.8e-8,3e-8,5000,false);
	  if(q2Bin==2) continue;
          if(q2Bin==3) continue;
	  if(q2Bin==4)StandardHypoTestInvDemo(fIn, "m_ws", "ModelConfig", "B_only_model", "data", 0, 0, true, 15,0.8e-8,4e-8,5000,false);                                                  	  
      }



  }


}


void D2hhmumuFitter_Applications::runFull1DFits(TString dataCut,TString misIDCut){

  D2hhmumuFitter1D* myFitter1D;
  int counter=0;

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

    myFitter1D= new D2hhmumuFitter1D(kind+"_"+*it+".txt");

    myFitter1D->setPathToSignalMC(pathToSignalMC);
    myFitter1D->setPathToNormMC(pathToNormMC);
    myFitter1D->setPathToKpipipiData(pathToKpipipiData);
    myFitter1D->setPathToHHpipiData(pathToKKpipiData);
    myFitter1D->setPathToNormData(pathToNormData);
    myFitter1D->setPathToSignalData(pathToSignalData);

    if(kind=="D2KKmumu"){


      myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_KKmumu_bin_%i.eps",counter),"m(K^{+}K^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangesKK_low[counter],rangesKK_high[counter]));
      myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_norm_KKTrained.eps");
      std::cout<<"Monte Carlo fits done.."<<std::endl;
      myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpipipi_KKTrained.eps");
      myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_KKmumu_bin_%i.eps",counter),"m_{K^{+}K^{-}#mu^{+}#mu^{-}}(KK#pi#pi)",TString::Format("%.0f<m(#pi#pi)<%.0fMeV/c^{2}",rangesKK_low[counter],rangesKK_high[counter])); 
      std::cout<<"misID fits done.."<<std::endl;
      myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_norm_KKTrained.eps");
      myFitter1D->fit_Data(kind,dataCut,+(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_blind_KKmumu_bin_%i.eps",counter),"m(K^{+}K^{-}#mu^{+}#mu^{-})",TString::Format("%.0fMeV<m(#mu^{+}#mu^{-})<%.0f/c^{2}",rangesKK_low[counter],rangesKK_high[counter]),true);
      std::cout<<"blind data fits done.."<<std::endl;

      ++counter;
    
    }

    if(kind=="D2pipimumu"){


      myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_pipimumu_bin_%i.eps",counter),"m(#pi^{+}#pi^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangespipi_low[counter],rangespipi_high[counter]));
      myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_norm_pipiTrained.eps");
      std::cout<<"Monte Carlo fits done.."<<std::endl;
      myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpipipi_pipiTrained.eps");
      myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_pipimumu_bin_%i.eps",counter),"m_{#pi^{+}#pi^{-}#mu^{+}#mu^{-}}(#pi#pi#pi#pi)",TString::Format("%.0f<m(#pi#pi)<%.0fMeV/c^{2}",rangespipi_low[counter],rangespipi_high[counter])); 
      std::cout<<"misID fits done.."<<std::endl;
      myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_norm_pipiTrained.eps");
      myFitter1D->fit_Data(kind,dataCut,+(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_blind_pipimumu_bin_%i.eps",counter),"m(#pi^{+}#pi^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangespipi_low[counter],rangespipi_high[counter]),true);
      std::cout<<"blind data fits done.."<<std::endl;
      
      ++counter;
    }
  }
}


void D2hhmumuFitter_Applications::runFullResonant1DFits(TString dataCut,TString misIDCut){


  D2hhmumuFitter1D* myFitter1D;
  int counter=0;

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

      myFitter1D= new D2hhmumuFitter1D(kind+"_unblinded_"+*it+".txt");

      myFitter1D->setPathToSignalMC(pathToSignalMC);
      myFitter1D->setPathToNormMC(pathToNormMC);
      myFitter1D->setPathToKpipipiData(pathToKpipipiData);
      myFitter1D->setPathToHHpipiData(pathToKKpipiData);
      myFitter1D->setPathToNormData(pathToNormData);
      myFitter1D->setPathToSignalData(pathToSignalData);


      //if(kind=="D2KKmumu" && (counter!=0 && counter!=1 && counter!=2 ) ){ //we're unblinded now!
      if(kind=="D2KKmumu" ){

	myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_KKmumu_bin_%i.eps",counter),"m(K^{+}K^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangesKK_low[counter],rangesKK_high[counter]));
	myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_norm_KKTrained.eps");
	std::cout<<"Monte Carlo fits done.."<<std::endl;
	myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpipipi_KKTrained.eps");
	myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_KKmumu_bin_%i.eps",counter),"m_{K^{+}K^{-}#mu^{+}#mu^{-}}(KK#pi#pi)",TString::Format("%.0f<m(#pi#pi)<%.0fMeV/c^{2}",rangesKK_low[counter],rangesKK_high[counter])); 
	std::cout<<"misID fits done.."<<std::endl;
	myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_norm_KKTrained.eps");
	myFitter1D->fit_unblinded_Data(kind,dataCut,+(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_KKmumu_bin_%i.eps",counter),"m(K^{+}K^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangesKK_low[counter],rangesKK_high[counter]),true);
	std::cout<<"data fits done.."<<std::endl;

      }
 
      //      if(kind=="D2pipimumu" && (counter!=0 && counter!=1 && counter!=4)){ // we are unblinded!
        if(kind=="D2pipimumu"){
	
	myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_pipimumu_bin_%i.eps",counter),"m(#pi^{+}#pi^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangespipi_low[counter],rangespipi_high[counter]));
	myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_norm_pipiTrained.eps");
	std::cout<<"Monte Carlo fits done.."<<std::endl;
	myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpipipi_pipiTrained.eps");
	myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_pipimumu_bin_%i.eps",counter),"m_{#pi^{+}#pi^{-}#mu^{+}#mu^{-}}(#pi#pi#pi#pi)",TString::Format("%.0f<m(#pi#pi)<%.0fMeV/c^{2}",rangespipi_low[counter],rangespipi_high[counter])); 
	std::cout<<"misID fits done.."<<std::endl;
	myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_norm_pipiTrained.eps");
	myFitter1D->fit_unblinded_Data(kind,dataCut,+(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_pipimumu_bin_%i.eps",counter),"m(#pi^{+}#pi^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangespipi_low[counter],rangespipi_high[counter]),true);
	std::cout<<" data fits done.."<<std::endl;
      
      }
      ++counter;
  }
  
  std::cout<<"data fits done.."<<std::endl;

}


void D2hhmumuFitter_Applications::addAllSignalSweights(TString dataCut,TString misIDCut){


  D2hhmumuFitter1D* myFitter1D;
  int counter=0;

  TString targetFolder = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/sWeights/";

  for (std::vector<TString>::iterator it = q2Ranges.begin() ; it != q2Ranges.end(); ++it) {

      myFitter1D= new D2hhmumuFitter1D(kind+"_sWeights_"+*it+".txt");

      myFitter1D->setPathToSignalMC(pathToSignalMC);
      myFitter1D->setPathToNormMC(pathToNormMC);
      myFitter1D->setPathToKpipipiData(pathToKpipipiData);
      myFitter1D->setPathToHHpipiData(pathToKKpipiData);
      myFitter1D->setPathToNormData(pathToNormData);
      myFitter1D->setPathToSignalData(pathToSignalData);


      //if(kind=="D2KKmumu" && (counter!=0 && counter!=1 && counter!=2 ) ){ //we're unblinded now!
      if(kind=="D2KKmumu" ){

	myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_KKmumu_bin_%i.eps",counter),"m(K^{+}K^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangesKK_low[counter],rangesKK_high[counter]));
	myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_norm_KKTrained.eps");
	std::cout<<"Monte Carlo fits done.."<<std::endl;
	myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpipipi_KKTrained.eps");
	myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_KKmumu_bin_%i.eps",counter),"m_{K^{+}K^{-}#mu^{+}#mu^{-}}(KK#pi#pi)",TString::Format("%.0f<m(#pi#pi)<%.0fMeV/c^{2}",rangesKK_low[counter],rangesKK_high[counter])); 
	std::cout<<"misID fits done.."<<std::endl;
	myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_norm_KKTrained.eps");
	myFitter1D->addSignalSWeights(kind,dataCut,(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/sWeightFit_KKmumu_bin_%i.eps",counter),"m(K^{+}K^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangesKK_low[counter],rangesKK_high[counter]),true,targetFolder+TString::Format("sWeights_KKmumu_bin_%i.root",counter));

	std::cout<<"data fits done.."<<std::endl;

      }
 
      //      if(kind=="D2pipimumu" && (counter!=0 && counter!=1 && counter!=4)){ // we are unblinded!
        if(kind=="D2pipimumu"){
	
	myFitter1D->fit_MC(dataCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_pipimumu_bin_%i.eps",counter),"m(#pi^{+}#pi^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangespipi_low[counter],rangespipi_high[counter]));
	myFitter1D->fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_norm_pipiTrained.eps");
	std::cout<<"Monte Carlo fits done.."<<std::endl;
	myFitter1D->fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpipipi_pipiTrained.eps");
	myFitter1D->fit_HHpipi_misID(misIDCut+"&&"+(*it),true,TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_pipimumu_bin_%i.eps",counter),"m_{#pi^{+}#pi^{-}#mu^{+}#mu^{-}}(#pi#pi#pi#pi)",TString::Format("%.0f<m(#pi#pi)<%.0fMeV/c^{2}",rangespipi_low[counter],rangespipi_high[counter])); 
	std::cout<<"misID fits done.."<<std::endl;
	myFitter1D->fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/data_norm_pipiTrained.eps");

	myFitter1D->addSignalSWeights(kind,dataCut,(*it),TString::Format("/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/sWeightFit_pipimumu_bin_%i.eps",counter),"m(#pi^{+}#pi^{-}#mu^{+}#mu^{-})",TString::Format("%.0f<m(#mu^{+}#mu^{-})<%.0fMeV/c^{2}",rangespipi_low[counter],rangespipi_high[counter]),true,targetFolder+TString::Format("sWeights_pipimumu_bin_%i.root",counter));

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

    myFitter1D->setPathToSignalMC(pathToSignalMC);
    myFitter1D->setPathToNormMC(pathToNormMC);
    myFitter1D->setPathToKpipipiData(pathToKpipipiData);
    myFitter1D->setPathToHHpipiData(pathToKKpipiData);
    myFitter1D->setPathToNormData(pathToNormData);
    myFitter1D->setPathToSignalData(pathToSignalData);
    

    if(kind=="D2KKmumu"){
    

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




void D2hhmumuFitter_Applications::performAllToyStudies(){

  //////////////////////////////////                                                                                                                                                          
  // TOY STUDY FOR FITTER VALIDATION                                                                                                                                                          
  /*
  D2hhmumuFitter1D* myFitter1D = new D2hhmumuFitter1D();
  myFitter1D->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");
  myFitter1D->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root");
  myFitter1D->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_Randomized.root");
  myFitter1D->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_Randomized.root");
  myFitter1D->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root");
  myFitter1D->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root");
  myFitter1D->makeToyStudy("D2pipimumu","BDT>.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass<525","BDT>0.4&&mu1_ProbNNmu>0.5","/work/mitzel/D2hhmumu/dev/D2pipimumu/fitValidation/1000Toys/D2pipimumu_FitterValidation_BRe_9_positiveBkg.root",0.5,24,20);
  myFitter1D->makeToyStudy("D2pipimumu","BDT>.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass<525","BDT>0.4&&mu1_ProbNNmu>0.5","/work/mitzel/D2hhmumu/dev/D2pipimumu/fitValidation/1000Toys/D2pipimumu_FitterValidation_BRe_8_positiveBkg.root",5,24,20);
  myFitter1D->makeToyStudy("D2pipimumu","BDT>.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass<525","BDT>0.4&&mu1_ProbNNmu>0.5","/work/mitzel/D2hhmumu/dev/D2pipimumu/fitValidation/1000Toys/D2pipimumu_FitterValidation_BRe_7_positiveBkg.root",50,24,20);
  myFitter1D->makeToyStudy("D2pipimumu","BDT>.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass<525","BDT>0.4&&mu1_ProbNNmu>0.5","/work/mitzel/D2hhmumu/dev/D2pipimumu/fitValidation/1000Toys/D2pipimumu_FitterValidation_BRe_6_positiveBkg.root",500,24,20);
  */

  D2hhmumuFitter1D* myFitter1D2 = new D2hhmumuFitter1D();
  myFitter1D2->setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D2->setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT.root");
  myFitter1D2->setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT.root");
  myFitter1D2->setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_Randomized.root");
  myFitter1D2->setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D2->setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT.root");
  
  myFitter1D2->makeToyStudy("D2KKmumu","BDT>.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass<525","BDT>0.4&&mu1_ProbNNmu>0.5","/work/mitzel/D2hhmumu/dev/D2KKmumu/fitValidation/1000Toys/D2KKmumu_FitterValidation_BRe_9_positiveBkg.root",0.5,7,1);                                                                                   
  myFitter1D2->makeToyStudy("D2KKmumu","BDT>.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass<525","BDT>0.4&&mu1_ProbNNmu>0.5","/work/mitzel/D2hhmumu/dev/D2KKmumu/fitValidation/1000Toys/D2KKmumu_FitterValidation_BRe_8_positiveBkg.root",5,7,1);
  myFitter1D2->makeToyStudy("D2KKmumu","BDT>.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass<525","BDT>0.4&&mu1_ProbNNmu>0.5","/work/mitzel/D2hhmumu/dev/D2KKmumu/fitValidation/1000Toys/D2KKmumu_FitterValidation_BRe_7_positiveBkg.root",50,7,1);
  myFitter1D2->makeToyStudy("D2KKmumu","BDT>.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass<525","BDT>0.4&&mu1_ProbNNmu>0.5","/work/mitzel/D2hhmumu/dev/D2KKmumu/fitValidation/1000Toys/D2KKmumu_FitterValidation_BRe_6_positiveBkg.root",500,7,1);
  myFitter1D2->makeToyStudy("D2KKmumu","BDT>.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","D_DiMuon_Mass<525","BDT>0.4&&mu1_ProbNNmu>0.5","/work/mitzel/D2hhmumu/dev/D2KKmumu/fitValidation/1000Toys/test_inv.root",0.5,8,1);
  
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

  TString PathToFolder = "/work/mitzel/D2hhmumu/dev/"+kind+"/img/studyResolutionScale/";

  D2hhmumuFitter1D* myFitter1D;
  int counter=0;

  myFitter1D = new D2hhmumuFitter1D("studyResolutionScale/results.txt"); 

  myFitter1D->setPathToSignalMC(pathToSignalMC);
  myFitter1D->setPathToNormMC(pathToNormMC);
  myFitter1D->setPathToKpipipiData(pathToKpipipiData);
  myFitter1D->setPathToHHpipiData(pathToKKpipiData);
  myFitter1D->setPathToNormData(pathToNormData);
  myFitter1D->setPathToSignalData(pathToSignalData);

  myFitter1D->fit_MC("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,"test.eps","m(KK#mu#mu)","950MeV<m(#mu#mu)<1100MeV");
  myFitter1D->fit_normalization_MC("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5",true,PathToFolder+"MC_Fit.eps");

  myFitter1D->fit_Kpipipi_misID("BDT>0.4&&mu1_ProbNNmu>0.5",true,PathToFolder+"Kpipipi_Fit.eps",true);
  myFitter1D->fit_HHpipi_misID("BDT>0.4&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,PathToFolder+"HHpipi_Fit.eps","m_{#pi#pi#mu#mu}(#pi#pi#pi#pi)","950MeV<m(#pi#pi)<1100MeV",true);

  myFitter1D->ResolutionScale.setRange(0.9,1.2);                                                                                                                              
  myFitter1D->globalShift.setRange(0.9,1.2);

  myFitter1D->fit_unblinded_Data("D2pipimumu","BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"freeFactor_Fit.eps","m(#pi#pi#mu#mu)","950MeV<m(#mu#mu)<1100MeV",false);
  myFitter1D->fit_normalization_Data("BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","test4.eps");
  myFitter1D->fit_unblinded_Data("D2pipimumu","BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",PathToFolder+"fixedFactor_Fit.eps","m(#pi#pi#mu#mu)","950MeV<m(#mu#mu)<1100MeV",false);


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

  */


}


void D2hhmumuFitter_Applications::studyCombBkgShape(){

  D2hhmumuFitter1D* myFitter1D;
  int counter=0;

  myFitter1D = new D2hhmumuFitter1D("studyCombBkgShape/results.txt"); 

  myFitter1D->setPathToSignalMC(pathToSignalMC);
  myFitter1D->setPathToNormMC(pathToNormMC);
  myFitter1D->setPathToKpipipiData(pathToKpipipiData);
  myFitter1D->setPathToHHpipiData(pathToKKpipiData);
  myFitter1D->setPathToNormData(pathToNormData);
  myFitter1D->setPathToSignalData(pathToSignalData);
 
  myFitter1D->fit_MC("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,"test.eps","m(KK#mu#mu)","950MeV<m(#mu#mu)<1100MeV");
  myFitter1D->fit_normalization_MC("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5",true,"test1.eps");

  myFitter1D->fit_Kpipipi_misID("BDT>0.4&&mu1_ProbNNmu>0.5",true,"test2.eps",true);
  myFitter1D->fit_HHpipi_misID("BDT>0.4&&mu1_ProbNNmu>0.5&&D_DiMuon_Mass>950&&D_DiMuon_Mass<1100",true,"test3.eps","m_{#pi#pi#mu#mu}(#pi#pi#pi#pi)","950MeV<m(#pi#pi)<1100MeV",true);
  myFitter1D->fit_normalization_Data("BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","test4.eps");
  //now all parameters fixed except for gamma

  myFitter1D->fit_invertedBDT_data("mu0_ProbNNmu>.5&&mu1_ProbNNmu>.5&&BDT<.4","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100","/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyCombBkgShape/v1.eps","","");
  myFitter1D->fit_invertedBDT_data("mu0_ProbNNmu>.5&&mu1_ProbNNmu>.5&&BDT<.4&&deltaM<165","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100","/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyCombBkgShape/v2.eps","","");
  myFitter1D->fit_invertedBDT_data("mu0_ProbNNmu>.5&&mu1_ProbNNmu>.5&&BDT<.4&&deltaM>165","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100","/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyCombBkgShape/v3.eps","","");

  myFitter1D->fit_invertedBDT_data("mu0_ProbNNmu>.5&&mu1_ProbNNmu>.5&&BDT<1","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100","/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyCombBkgShape/v4.eps","","");
  myFitter1D->fit_invertedBDT_data("mu0_ProbNNmu>.5&&mu1_ProbNNmu>.5&&BDT<1&&deltaM<165","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100","/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyCombBkgShape/v5.eps","","");
  myFitter1D->fit_invertedBDT_data("mu0_ProbNNmu>.5&&mu1_ProbNNmu>.5&&BDT<1&&deltaM>165","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100","/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyCombBkgShape/v6.eps","","");

  myFitter1D->fit_invertedBDT_data("mu0_ProbNNmu>.5&&mu1_ProbNNmu>.5&&BDT<.0","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100","","","");
  myFitter1D->fit_invertedBDT_data("mu0_ProbNNmu>.5&&mu1_ProbNNmu>.5&&BDT<.0&&deltaM<165","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100","/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyCombBkgShape/v7.eps","","");
  myFitter1D->fit_invertedBDT_data("mu0_ProbNNmu>.5&&mu1_ProbNNmu>.5&&BDT<.0&&deltaM>165","D_DiMuon_Mass>950&&D_DiMuon_Mass<1100","/work/mitzel/D2hhmumu/dev/D2pipimumu/img/studyCombBkgShape/v8.eps","","");
  
 
  //myFitter1D->fit_normalization_Data("BDT>0.4&&mu1_ProbNNmu>0.5&&mu0_ProbNNmu>0.5","test4.eps");
  //myFitter1D->fit_unblinded_Data("D2pipimumu","BDT>0.4&&mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5","D_DiMuon_Mass<950&&D_DiMuon_Mass>565","test5.eps","m(#pi#pi#mu#mu)","950MeV<m(#mu#mu)<1100MeV",true);

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
