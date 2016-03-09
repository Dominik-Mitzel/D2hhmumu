#include "D2hhmumuFitter.h"
#include "RooAddPdf.h"
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
  RooRealVar deltaM_xi("deltaM_xi","deltaM_xi",1.45437e+02,144,146);
  RooRealVar deltaM_lambda("deltaM_lambda","deltaM_lambda",5.35039e-01,0.1,2);
  RooRealVar deltaM_gamma("deltaM_gamma","deltaM_gamma",1.00164e-01,-2,2);
  RooRealVar deltaM_delta("deltaM_delta","deltaM_delta",1.07829e+00,0.,10);

  RooRealVar D0_M_xi("D0_M_xi","D0_M_xi",1.86560e+03,1860,1870);
  RooRealVar D0_M_lambda("D0_M_lambda","D0_M_lambda",8.54839e+00,0.1,20);
  RooRealVar D0_M_gamma("D0_M_gamma","D0_M_gamma",-5.42758e-02,-2,2);
  RooRealVar D0_M_delta("D0_M_delta","D0_M_delta",6.09288e-01,0.,10);

  D2hhmumuModel myMCModel;
  myMCModel.Signal(D0_M,deltaM,
		   D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
		   deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta
		   );

  std::string components="Signal";
  myMCModel.Model(components);


  RooFitResult *result;
  //result = Gauss2.fitTo(*data_binned,Save(kTRUE));                                                                                                                                          
  result = myMCModel.GetWorkspace().pdf("D2hhmumuModel")->fitTo(*data,Save(kTRUE),NumCPU(3));
  
  ///Plot                                                                                                                                                                                     
  ///----------                                                                                                                                                                               
  TCanvas* c1= new TCanvas("");
  c1->Divide(2);
  c1->cd(1);
  
  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));

  myMCModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,LineColor(kRed),LineWidth(1));
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  
  
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));

  myMCModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,LineColor(kRed),LineWidth(1));
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


  RooRealVar deltaM_xi("deltaM_xi","deltaM_xi",1.45437e+02,144,146);
  RooRealVar deltaM_lambda("deltaM_lambda","deltaM_lambda",5.35039e-01,0.1,2);
  RooRealVar deltaM_gamma("deltaM_gamma","deltaM_gamma",1.00164e-01,-2,2);
  RooRealVar deltaM_delta("deltaM_delta","deltaM_delta",1.07829e+00,0.,10);

  RooRealVar D0_M_xi("D0_M_xi","D0_M_xi",1.86560e+03,1860,1870);
  RooRealVar D0_M_lambda("D0_M_lambda","D0_M_lambda",8.54839e+00,0.1,20);
  RooRealVar D0_M_gamma("D0_M_gamma","D0_M_gamma",-5.42758e-02,-2,2);
  RooRealVar D0_M_delta("D0_M_delta","D0_M_delta",6.09288e-01,0.,10);

  D2hhmumuModel myModel;
  RooAbsPdf *my_signal = myModel.Signal(D0_M,deltaM,
                   D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
                   deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta
                   );


  //background models

  //purely combinatorial
  RooRealVar deltaM_threshold("deltaM_threshold","deltaM_threshold",139.57018);
  RooRealVar deltaM_alpha("deltaM_alpha","deltaM_alpha",3.9761e+00,0,10.);

  RooRealVar D0_M_chebyA("D0_M_chebyA","D0_M_chebyA",-3.5906e-02,-1,1);
  RooRealVar D0_M_chebyB("D0_M_chebyB","D0_M_chebyB",-1.7004e-02,-1,1);
  RooRealVar D0_M_chebyC("D0_M_chebyC","D0_M_chebyC",-1.7882e-02,-1,1);


  myModel.CombinatoricBackground(D0_M,deltaM,
				 D0_M_chebyA,D0_M_chebyB,D0_M_chebyC,
				 deltaM_threshold,deltaM_alpha
				 );


  //randompi
  myModel.RandomPionBackground(D0_M,deltaM,
			       D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
			       deltaM_threshold,deltaM_alpha
                               );

  //peaking , in mD0 jsut a Gaussian... to be improved!

  RooRealVar mean1("mu1", "mean1", 1.8350e3,1835.,1845.);                                                                                             
  RooRealVar sigma1("sigma_{1}", "sigma1", 1.4187e+01,3.,25.);                                                                                      

  myModel.D2hhhhBackground(D0_M,deltaM,
			   mean1,sigma1,
			   deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta
			   );

  myModel.D2hhhhRandomPionBackground(D0_M,deltaM,
                           mean1,sigma1,
			   deltaM_threshold,deltaM_alpha
                           );

  ///Fit                                                                                                                                                                                     
  std::string components="Signal CombinatoricBkg RandomPionBkg D2hhhhBkg D2hhhhRandomPionBkg";
  RooAbsPdf* finalPDF = myModel.Model(components);

  RooFitResult *result;
  //result = myModel.GetWorkspace().pdf("D2hhmumuModel")->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));                                                           
  result = finalPDF->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));

  cout << "result is --------------- "<<endl;
  result->Print();

  ///Plot                                                                                                                                                                                     
  ///----------                                                                                                                                                                               
  TCanvas* c1= new TCanvas("");
  c1->Divide(2,2);
  c1->cd(1);
  
  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel.GetWorkspace().pdf("Signal")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel.GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel.GetWorkspace().pdf("RandomPionBkg")),LineColor(kCyan),LineStyle(kDashed),LineWidth(1));
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel.GetWorkspace().pdf("D2hhhhBkg")),LineColor(kYellow),LineStyle(kDashed),LineWidth(1));
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel.GetWorkspace().pdf("D2hhhhRandomPionBkg")),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1)); 
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel.GetWorkspace().pdf("Signal")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel.GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel.GetWorkspace().pdf("RandomPionBkg")),LineColor(kCyan),LineStyle(kDashed),LineWidth(1));
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel.GetWorkspace().pdf("D2hhhhBkg")),LineColor(kYellow),LineStyle(kDashed),LineWidth(1)); 
  finalPDF->plotOn(frame_m2,LineColor(kBlack),LineWidth(1)); 
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m2->Draw();
  
  RooHist* hpull = frame_m->pullHist() ;
  RooPlot* frame_m3 = D0_M.frame(Title("Pull Distribution")) ;
  frame_m3->addPlotable(hpull,"P") ;
  c1->cd(3);
  frame_m3->Draw();

  RooHist* hpull2 = frame_m2->pullHist() ;
  RooPlot* frame_m4 = deltaM.frame(Title("Pull Distribution")) ;
  frame_m4->addPlotable(hpull2,"P") ;
  c1->cd(4);
  frame_m4->Draw();

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
