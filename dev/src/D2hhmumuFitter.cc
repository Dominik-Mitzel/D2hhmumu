#include "D2hhmumuFitter.h"
#include "RooAddPdf.h"
using namespace std;
using namespace RooFit ;

D2hhmumuFitter::D2hhmumuFitter():

  D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1820., 1940.,"MeV"),
  deltaM("deltaM","#delta m", 139.8,150,"MeV"),
  //signal                                                                                                                                                                      
  deltaM_xi("deltaM_xi","deltaM_xi",1.45437e+02,144,146),
  deltaM_lambda("deltaM_lambda","deltaM_lambda",5.35039e-01,0.1,2),
  deltaM_gamma("deltaM_gamma","deltaM_gamma",1.00164e-01,-2,2),
  deltaM_delta("deltaM_delta","deltaM_delta",1.07829e+00,0.,10),
  D0_M_xi("D0_M_xi","D0_M_xi",1.86560e+03,1860,1870),
  D0_M_lambda("D0_M_lambda","D0_M_lambda",8.54839e+00,0.1,20),
  D0_M_gamma("D0_M_gamma","D0_M_gamma",-5.42758e-02,-2,2),
  D0_M_delta("D0_M_delta","D0_M_delta",6.09288e-01,0.,10),
  //purely combinatorial background                                                                                                                                             
  deltaM_threshold("deltaM_threshold","deltaM_threshold",139.57018),
  deltaM_alpha("deltaM_alpha","deltaM_alpha",3.9761e+00,0,10.),
  D0_M_chebyA("D0_M_chebyA","D0_M_chebyA",-3.5906e-02,-1,1),
  D0_M_chebyB("D0_M_chebyB","D0_M_chebyB",-1.7004e-02,-1,1),
  D0_M_chebyC("D0_M_chebyC","D0_M_chebyC",-1.7882e-02,-1,1),
  //peaking , in mD0 just a Gaussian... to be improved!                                                                                                                         
  mean1("mu1", "mean1", 1.8350e3,1835.,1845.),
  sigma1("sigma_{1}", "sigma1", 1.4187e+01,3.,25.)

{};

D2hhmumuFitter::~D2hhmumuFitter(){}


void::D2hhmumuFitter::setKpimumuStartParameters(){                                                                                                                             
                                                 
  //later like this...                                                                                                                                                         
  //deltaM_xi.SetValue()..
  
}


void D2hhmumuFitter::fit_MC(){

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
  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", tree, RooArgSet(D0_M,deltaM));
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));
  RooDataHist* data_binned = data_small->binnedClone();

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


void D2hhmumuFitter::fit_Data(){

  ///Load file                                                                                                                                                                                
  TFile* file;
  TString pathToFile="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_BDT_selected.root"; //"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root";
  file= new TFile(pathToFile,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);

  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", tree, RooArgSet(D0_M,deltaM));
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));
  RooDataHist* data_binned = data_small->binnedClone();

  D2hhmumuModel myModel;
  RooAbsPdf *my_signal = myModel.Signal(D0_M,deltaM,
                   D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
                   deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta
                   );


  //background models

  myModel.CombinatoricBackground(D0_M,deltaM,
				 D0_M_chebyA,D0_M_chebyB,D0_M_chebyC,
				 deltaM_threshold,deltaM_alpha
				 );


  //randompi
  myModel.RandomPionBackground(D0_M,deltaM,
			       D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
			       deltaM_threshold,deltaM_alpha
                               );

  //peaking , in mD0 just a Gaussian... to be improved!

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
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel.GetWorkspace().pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1));
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
  myModel.GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel.GetWorkspace().pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1)); 
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
  c1->Print("../img/massFit2.eps");
}


void D2hhmumuFitter::setStyle(){

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

}
