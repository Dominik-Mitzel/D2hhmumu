#include "D2hhmumuFitter.h"

using namespace std;
using namespace RooFit ;

void fit_MC(){


  gStyle->SetOptStat(0);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelOffset(0.01,"X");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  ///Load file                                                                                                                                                                                
  TFile* file;
  TString pathToFile="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root"; //"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root";
  file= new TFile(pathToFile,"OPEN");
  TTree* tree = (TTree*) file->Get("DecayTree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);

  ///Fill all needed variables in RooDataSet                                                                                                                                                  
  ///---------------------------------------                                                                                                                                                  

  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1800., 1940.,"MeV");
  RooRealVar deltaM("deltaM","#delta m", 139.8,150,"MeV");
  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", tree, RooArgSet(D0_M,deltaM));
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));
  RooDataHist* data_binned = data_small->binnedClone();
 ///Define fit model                                                                                                                                                                         
  ///----------------------- delta M                                                                                                                                                                 
  RooRealVar deltaM_xi("deltaM_xi","deltaM_xi",145,144,146);
  RooRealVar deltaM_lambda("deltaM_lambda","deltaM_lambda",0.3,0.1,2);
  RooRealVar deltaM_gamma("deltaM_gamma","deltaM_gamma",0.,-2,2);
  RooRealVar deltaM_delta("deltaM_delta","deltaM_delta",0.3,0.,10);

  RooJohnsonSU deltaM_Johnson("deltaM_J","Johnson SU distribution",deltaM,deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta);

  // mD0

  RooRealVar D0_M_xi("D0_M_xi","D0_M_xi",1865,1860,1870);
  RooRealVar D0_M_lambda("D0_M_lambda","D0_M_lambda",7,0.1,20);
  RooRealVar D0_M_gamma("D0_M_gamma","D0_M_gamma",0.,-2,2);
  RooRealVar D0_M_delta("D0_M_delta","D0_M_delta",0.3,0.,10);

  RooJohnsonSU D0_M_Johnson("D0_M_J","Johnson SU distribution",D0_M,D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta);
 
  //2D fit function

 RooProdPdf signal_model("signal_model","signal_model",RooArgSet(D0_M_Johnson,deltaM_Johnson));

  RooFitResult *result;
  //result = Gauss2.fitTo(*data_binned,Save(kTRUE));                                                                                                                                          
  result = signal_model.fitTo(*data,Save(kTRUE),NumCPU(3));
  
  ///Plot                                                                                                                                                                                     
  ///----------                                                                                                                                                                               
  TCanvas* c1= new TCanvas("");
  c1->Divide(2);
  c1->cd(1);
  
  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));

  signal_model.plotOn(frame_m,LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  
  
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));

  signal_model.plotOn(frame_m2,LineColor(kRed),LineWidth(1));
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m2->Draw();

  c1->Draw();
  c1->Print("massFit_MC.eps");
}


void fit_Data(){

  gStyle->SetOptStat(0);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelOffset(0.01,"X");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  ///Load file                                                                                                                                                                                
  TFile* file;
  TString pathToFile="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_BDT_selected.root"; //"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root";
  file= new TFile(pathToFile,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);

  ///Fill all needed variables in RooDataSet                                                                                                                                                  
  ///---------------------------------------                                                                                                                                                  

  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1820., 1940.,"MeV");
  RooRealVar deltaM("deltaM","#delta m", 139.8,150,"MeV");
  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", tree, RooArgSet(D0_M,deltaM));
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));
  RooDataHist* data_binned = data_small->binnedClone();

  ///Define fit model                                                                                                                                                                         
  ///----------------------- delta M                                                                                                                                                                 
  RooRealVar deltaM_xi("deltaM_xi","deltaM_xi",145,144,146);
  RooRealVar deltaM_lambda("deltaM_lambda","deltaM_lambda",0.3,0.1,2);
  RooRealVar deltaM_gamma("deltaM_gamma","deltaM_gamma",0.,-2,2);
  RooRealVar deltaM_delta("deltaM_delta","deltaM_delta",0.3,0.,10);

  RooJohnsonSU deltaM_Johnson("deltaM_J","Johnson SU distribution",deltaM,deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta);

  /// mD0

  RooRealVar D0_M_xi("D0_M_xi","D0_M_xi",1865,1860,1870);
  RooRealVar D0_M_lambda("D0_M_lambda","D0_M_lambda",7,0.1,20);
  RooRealVar D0_M_gamma("D0_M_gamma","D0_M_gamma",0.,-2,2);
  RooRealVar D0_M_delta("D0_M_delta","D0_M_delta",0.3,0.,10);

  RooJohnsonSU D0_M_Johnson("D0_M_J","Johnson SU distribution",D0_M,D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta);


  //background models

  //purely combinatorial
  RooRealVar deltaM_threshold("deltaM_threshold","deltaM_threshold",139.57018);
  RooRealVar deltaM_alpha("deltaM_alpha","deltaM_alpha",2.1799e-02,0,10.);

  RooThreshold deltaM_combinatorial("deltaM_combinatorial","deltaM combinatorial",deltaM,deltaM_threshold,deltaM_alpha);

  //randompi
  RooRealVar deltaM_alpha_randompi("deltaM_alpha_randompi","deltaM_alpha random pi",2.1799e-02,0,10.);
  RooThreshold deltaM_randompi("deltaM_randompi","deltaM randompi",deltaM,deltaM_threshold,deltaM_alpha_randompi);

  RooJohnsonSU D0_M_Johnson_randompi("D0_M_J_randompi","Johnson SU distribution",D0_M,D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta);


  RooRealVar D0_M_chebyA("D0_M_chebyA","D0_M_chebyA",0,-1,1);
  RooRealVar D0_M_chebyB("D0_M_chebyB","D0_M_chebyB",0,-1,1);
  RooRealVar D0_M_chebyC("D0_M_chebyC","D0_M_chebyC",0,-1,1);

  RooChebychev D0_M_combinatorial("D0_M_combinatorial","D0_M_combinatorial",D0_M,RooArgList(D0_M_chebyA,D0_M_chebyB,D0_M_chebyC));
  RooRealVar f_D0_M_sig("f_sig","fraction of CB in D0 mass signal",0.1,0.,1.);
  RooRealVar f_D0_M_sig2("f_sig2","fraction of CB in D0 mass signal",0.1,0.,1.);


  //peaking , in mD0 jsut a Gaussian... to be improved!

  RooRealVar mean1("mu1", "mean1", 1840,1835.,1845.);                                                                                                           
  RooRealVar sigma1("sigma_{1}", "sigma1", 5.45,3.,25.);                                                                                                        
  RooGaussian Gauss1("Gauss1", "Gauss1", D0_M, mean1, sigma1);  
 
  RooJohnsonSU deltaM_Johnson_D2hhhh("deltaM_J_D2hhhh","Johnson SU distribution",deltaM,deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta);

  RooAddPdf total_D0_M("PDF_D0_M_sig","signal PDF in D0_M",RooArgSet(D0_M_combinatorial,D0_M_Johnson,Gauss1),RooArgSet(f_D0_M_sig,f_D0_M_sig2));

  //total pdf

  RooRealVar f_deltaM_sig("f_sig","fraction of CB in D0 mass signal",0.1,0.,1.);                                                                                
  RooAddPdf PDF_deltaM_sig("PDF_deltaM_sig","signal PDF in deltaM",RooArgSet(deltaM_combinatorial,deltaM_Johnson),f_deltaM_sig);   
  


  ///signal model                                                                                                                                                                         
  ///-------------------------                                                                                                                                                                
  RooProdPdf signal_model("signal_model","signal_model",RooArgSet(D0_M_Johnson,deltaM_Johnson));
  RooProdPdf combinatorial_model("combinatorial_model","model",RooArgSet(D0_M_combinatorial,deltaM_combinatorial));
  RooProdPdf randompi_model("randompi_model","model",RooArgSet(D0_M_Johnson_randompi,deltaM_randompi));
  RooProdPdf D2hhhh_model("D2hhhh_model","model",RooArgSet(Gauss1,deltaM_Johnson_D2hhhh));
  
  RooRealVar nSignal("nSignal","number of signal events",5000,0,20e3);
  RooExtendPdf eSignalPDF("eSignalPDF","extended signal PDF 2D",signal_model,nSignal);

  RooRealVar f1("f1","fraction of CB in D0 mass signal",0.1,0.,1.);
  RooRealVar f2("f2","fraction of CB in D0 mass signal",0.1,0.,1.);
  RooRealVar f3("f3","fraction of CB in D0 mass signal",0.1,0.,1.);

  RooRealVar nBkg("nBkg","number of bkg events",5000, 0,20e3);
  RooAddPdf bkg_model("bkg_model","total bkg pdf",RooArgSet(combinatorial_model,randompi_model,D2hhhh_model),RooArgSet(f1,f2));
  RooExtendPdf eBKGPDF("eBKGPDF","extended bkg PDF 2D",bkg_model,nBkg);

  //RooAddPdf total_bkg("total","total pdf",RooArgSet(combinatorial_model,randompi_model,D2hhhh_model),RooArgSet(f1,f2));
  ///RooAddPdf total("total","total pdf",RooArgSet(signal_model,combinatorial_model,randompi_model,D2hhhh_model),RooArgSet(f1,f2,f3));
  RooAddPdf total("total","total pdf",RooArgSet(eSignalPDF,eBKGPDF));


  ///Fit                                                                                                                                                                                     
  // RooFitResult *result;                                                                                                                                                                    
  //if(binned) result = pdf->fitTo(*data_binned,Save(kTRUE),Extended());                                                                                                                      
  //else result = pdf->fitTo(*data,Save(kTRUE),Extended(),NumCPU(3));                                                                                                                         

  RooFitResult *result;
  //result = Gauss2.fitTo(*data_binned,Save(kTRUE));                                                                                                                                          
  result = total.fitTo(*data,Save(kTRUE),NumCPU(3));
  //result = PDF_D0_M_sig.fitTo(*data_small,Save(kTRUE),NumCPU(3));                                                                                                          


  cout << "result is --------------- "<<endl;
  result->Print();

  ///Plot                                                                                                                                                                                     
  ///----------                                                                                                                                                                               
  TCanvas* c1= new TCanvas("");
  c1->Divide(2);
  c1->cd(1);
  
  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  total.plotOn(frame_m,LineColor(kRed),LineStyle(kDashed),LineWidth(1));
   total.plotOn(frame_m,Components(combinatorial_model),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
   total.plotOn(frame_m,Components(randompi_model),LineColor(kCyan),LineStyle(kDashed),LineWidth(1));
   total.plotOn(frame_m,Components(D2hhhh_model),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
   data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  
  
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  total.plotOn(frame_m2,LineColor(kRed),LineWidth(1));
   total.plotOn(frame_m2,Components(combinatorial_model),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
   total.plotOn(frame_m2,Components(randompi_model),LineColor(kCyan),LineStyle(kDashed),LineWidth(1));
   total.plotOn(frame_m2,Components(D2hhhh_model),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
   data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m2->Draw();
  
std::cout<<"nentries"<<tree->GetEntries()<<std::endl;	
  c1->Draw();
  c1->Print("massFit2.eps");
}
/*
function graveyard 

 

  RooRealVar mean1("mu1", "mean1", 145,144.,146.);
  RooRealVar sigma1("sigma_{1}", "sigma1", 0.45,0.1,2.);
  RooGaussian Gauss1("Gauss1", "Gauss1", deltaM, mean1, sigma1);
  RooRealVar mean2("mu2", "mean2", 145.,144.,147.);
  RooRealVar sigma2("sigma_{2}", "sigma2", 0.45,0.1,3.);
  RooGaussian Gauss2("Gauss2", "Gauss2", deltaM, mean1, sigma2);
  RooRealVar sigma3("sigma_{3}", "sigma3", 0.45,0.01,3.);
  RooGaussian Gauss3("Gauss3", "Gauss3", deltaM, mean1, sigma3);

  RooAddPdf PDF_deltaM_sig("PDF_deltaM_sig","signal PDF in deltaM",RooArgSet(Gauss1,Gauss2),RooArgSet(f_deltaM_sig));
  RooRealVar f_deltaM_sig("f_sig","fraction of CB in D0 mass sigal",0.8,0.,1.);
  RooRealVar f_deltaM_sig2("f_sig2","fraction of CB in D0 mass sigal",0.1,0.,1.);


    //better model                                                                                                                                                                              

  RooRealVar mean_mD0_sig("mean_CB_mD0_sig","mean of signal CB ",1866,1860,1880);
  RooRealVar sigma_mD0_sig("sigma_CB_mD0_sig","mean of signal CB ",7.7,5,10);
  RooRealVar alpha_mD0_sig1("alpha_mD0_sig1","alpha1 signal CB",2.0,1.,5.);
  RooRealVar alpha_mD0_sig2("alpha_mD0_sig2","alpha2 signal CB",1.5,1.,5.);
  RooRealVar n_mD0_sig1("n_mD0_sig1","n1 signal CB",1.8,0.,15.);
  RooRealVar n_mD0_sig2("n_mD0_sig2","n2 signal CB",20.,0.,30.);

  RooCBShape CB_mD0_sig1("CB_mD0_sig1","crystall ball 1 in mD0 signal PDF",D0_M,mean_mD0_sig,sigma_mD0_sig,alpha_mD0_sig1,n_mD0_sig1);
  RooCBShape CB_mD0_sig2("CB_mD0_sig2","crystall ball 2 in mD0 signal PDF",D0_M,mean_mD0_sig,sigma_mD0_sig,alpha_mD0_sig2,n_mD0_sig2);

  RooRealVar f_mD0_sig("f_sig","fraction of CB in D0 mass signal",0.8,0.,1.);
  RooAddPdf PDF_mD0_sig("PDF_mD0_sig","signal PDF in mD0",RooArgSet(CB_mD0_sig1,CB_mD0_sig2),f_mD0_sig);


  //delta M alternative with double crystal ball                                                                                                                                              

  RooRealVar mean_deltaM_sig("mean_CB_deltaM_sig","mean of signal CB ",144,145,146);
  RooRealVar sigma_deltaM_sig("sigma_CB_deltaM_sig","mean of signal CB ",0.1,0.1,3);
  RooRealVar alpha_deltaM_sig1("alpha_deltaM_sig1","alpha1 signal CB",2.0,0.,5.);
  RooRealVar alpha_deltaM_sig2("alpha_deltaM_sig2","alpha2 signal CB",1.5,0.,5.);
  RooRealVar n_deltaM_sig1("n_deltaM_sig1","n1 signal CB",0,0.,15.);
  RooRealVar n_deltaM_sig2("n_deltaM_sig2","n2 signal CB",0.,0.,30.);

  RooCBShape CB_deltaM_sig1("CB_deltaM_sig1","crystall ball 1 in deltaM signal PDF",deltaM,mean_deltaM_sig,sigma_deltaM_sig,alpha_deltaM_sig1,n_deltaM_sig1);
  RooCBShape CB_deltaM_sig2("CB_deltaM_sig2","crystall ball 2 in deltaM signal PDF",deltaM,mean_deltaM_sig,sigma_deltaM_sig,alpha_deltaM_sig2,n_deltaM_sig2);

  //  RooRealVar f_deltaM_sig("f_sig","fraction of CB in D0 mass signal",0.8,0.,1.);                                                                                                          
  //RooAddPdf PDF_deltaM_sig("PDF_deltaM_sig","signal PDF in deltaM",RooArgSet(CB_deltaM_sig1,CB_deltaM_sig2),f_deltaM_sig);                                                                 


  ///peaking bkg                                                                                                                                                                              
  RooRealVar mean_bkg("mu_{bkg}", "mean_bkg", 1840.,1820.,1860.);
  RooRealVar sigma_bkg("sigma_{bkg}", "sigma_bkg", 28.5, 5., 70.);
  RooGaussian bkg_Gauss("bkg_Gauss", "bkg_Gauss", D0_M, mean_bkg, sigma_bkg);

  ///Exponential                                                                                                                                                                              
  RooRealVar exp_par("lambda","exp",0.003535,-1.,1.);
  RooExponential bkg_exp("bkg_exp","exponential background",D0_M,exp_par);

  ///bkg pdf                                                                                                                                                                                  
  RooRealVar f_bkg("f_{bkg}", "background fraction", 0.11, 0., 1.);
  RooAddPdf bkg("bkg", "bkg", RooArgList(bkg_Gauss, bkg_exp), RooArgList(f_bkg));
  ///total pdf                                                                                                                                                                                
  ///----------------------                                                                                                                                                                   
  RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
  RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());

  RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(Gauss1,bkg), RooArgList(n_sig,n_bkg));

 
  */
