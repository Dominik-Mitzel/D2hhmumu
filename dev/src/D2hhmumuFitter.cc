#include "D2hhmumuFitter.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
using namespace std;
using namespace RooFit ;

//this class is is a fitter for Dst->D(hhmumu)pi decays.
//
//
//

D2hhmumuFitter::D2hhmumuFitter():

  D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1820., 1940.,"MeV"),
  //deltaM("deltaM","#delta m", 139.8,150,"MeV"),
  deltaM("deltaM","#delta m", 142,149,"MeV"),
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
  deltaM_xi_bkg("deltaM_xi_bkg","deltaM_xi_bkg",1.45437e+02,144,146),
  deltaM_lambda_bkg("deltaM_lambda_bkg","deltaM_lambda_bkg",5.35039e-01,0.1,2),
  deltaM_gamma_bkg("deltaM_gamma_bkg","deltaM_gamma_bkg",1.00164e-01,-2,2),
  deltaM_delta_bkg("deltaM_delta_bkg","deltaM_delta_bkg",1.07829e+00,0.,10),
  D0_M_xi_bkg("D0_M_xi_bkg","D0_M_xi_bkg",1.8312e+03,1.8212e+03,1.8412e+03),
  D0_M_lambda_bkg("D0_M_lambda_bkg","D0_M_lambda_bkg",4.0000e+01,1.0000e+01,10.0000e+01),
  D0_M_gamma_bkg("D0_M_gamma_bkg","D0_M_gamma_bkg",-5.8442e-01,-10.8442e-01,-5.8442e-01),
  D0_M_delta_bkg("D0_M_delta_bkg","D0_M_delta_bkg",7.5826e-01,1.5826e-01,17.5826e-01),

  mean1("mu1", "mean1", 1.8350e3,1835.,1845.),
  sigma1("sigma_{1}", "sigma1", 1.4187e+01,3.,25.)
{

  pathToSignalMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root";
  pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root";
  pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root";
  pathToInvData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/inverted_D2Kpimumu_BDT_selected.root";
  pathToSidebandData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT.root";


}

void D2hhmumuFitter::setPathToSignalMC(TString path){
  pathToSignalMC = path;
}
void D2hhmumuFitter::setPathToSignalData(TString path){
  pathToSignalData = path;
}
void D2hhmumuFitter::setPathToNormData(TString path){
  pathToNormData = path;
}
void D2hhmumuFitter::setPathToInvData(TString path){
  pathToInvData = path;
}
void D2hhmumuFitter::setPathToSidebandData(TString path){
  pathToSidebandData = path;
}


D2hhmumuFitter::~D2hhmumuFitter(){}


void D2hhmumuFitter::setKpimumuStartParameters(){                                                                                                                             
                                                 
  //later like this...                                                                                                                                                         
  //deltaM_xi.SetValue()..
  
}

void D2hhmumuFitter::setKKmumuStartParameters(){

  //later like this...                                                                                                                                                         
  //deltaM_xi.SetValue()..                                                                                                                                                      

}

void D2hhmumuFitter::setpipimumuStartParameters(){

  //later like this...                                                                                                                                                         
  //deltaM_xi.SetValue()..                                                                                                                                                      

}



void D2hhmumuFitter::fit_MC(bool fixShape = true){

  ///Load file                                                                                                                                                                                
  TFile* file;
  file= new TFile(pathToSignalMC,"OPEN");
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


  D2hhmumuModel* myModel = initializeModel();

  std::string components="Signal";
  RooAbsPdf* finalPDF = myModel->Model(components);


  RooFitResult *result;
  result = finalPDF->fitTo(*data,Save(kTRUE),NumCPU(3));
  
  ///Plot                                                                                                                                                                                     
  ///----------                                                                                                                                                                               
  TCanvas* c1= new TCanvas("");
  c1->Divide(2);
  c1->cd(1);
  
  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));

  finalPDF->plotOn(frame_m,LineColor(kRed),LineWidth(1));
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  
  
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));

  finalPDF->plotOn(frame_m2,LineColor(kRed),LineWidth(1));
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m2->Draw();

  c1->Draw();
  c1->Print("massFit_MC.eps");


  deltaM_xi.setVal(myModel->GetWorkspace().var("deltaM_xi")->getValV() );
  deltaM_lambda.setVal(myModel->GetWorkspace().var("deltaM_lambda")->getValV() );
  deltaM_gamma.setVal(myModel->GetWorkspace().var("deltaM_gamma")->getValV() );
  deltaM_delta.setVal(myModel->GetWorkspace().var("deltaM_delta")->getValV() );
  D0_M_xi.setVal(myModel->GetWorkspace().var("D0_M_xi")->getValV() );
  D0_M_lambda.setVal(myModel->GetWorkspace().var("D0_M_lambda")->getValV() );
  D0_M_gamma.setVal(myModel->GetWorkspace().var("D0_M_gamma")->getValV() );
  D0_M_delta.setVal(myModel->GetWorkspace().var("D0_M_delta")->getValV() );


  if(fixShape){

    deltaM_xi.setConstant();
    deltaM_lambda.setConstant();
    deltaM_gamma.setConstant();
    deltaM_delta.setConstant();
    D0_M_xi.setConstant();
    D0_M_lambda.setConstant();
    D0_M_gamma.setConstant();
    D0_M_delta.setConstant();
  }
  file->Close();

}


void D2hhmumuFitter::fit_Data(TString cut=""){

  ///Load file                                                                                                                                                                                
  TFile* file;
  file= new TFile(pathToSignalData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);  
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);

  TTree* cutTree = tree->CopyTree(cut);


  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M,deltaM));
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));
  RooDataHist* data_binned = data_small->binnedClone();

 
  ///Fit                                                                                                                                                                        
             
  D2hhmumuModel* myModel = initializeModel();
  std::string components="Signal CombinatoricBkg RandomPionBkg D2hhhhBkg";
  RooAbsPdf* finalPDF = myModel->Model(components);
  std::cout<<"ok.."<<std::endl;

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
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("Signal")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("RandomPionBkg")),LineColor(kCyan),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("D2hhhhRandomPionBkg")),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1)); 
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel->GetWorkspace().pdf("Signal")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel->GetWorkspace().pdf("RandomPionBkg")),LineColor(kCyan),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel->GetWorkspace().pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1)); 
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


  file->Close();


}


void D2hhmumuFitter::fit_PIDinverted_Data(bool fixShape){

  ///Load file                                                                                                                                                                 \
                                                                                                                                                                               
  TFile* file;
  file= new TFile(pathToInvData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);

  D0_M.setRange(1800,1940);

  D2hhmumuModel* myModel = initializeModel();

  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", tree, RooArgSet(D0_M,deltaM));
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));
  RooDataHist* data_binned = data_small->binnedClone();

  ///Fit                                                                                                          

  std::string components="D2hhhhBkg CombinatoricBkg";
  RooAbsPdf* finalPDF = myModel->Model(components);
  RooFitResult *result;
  result = finalPDF->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));

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
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("D2hhhhBkg")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1));
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();

  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel->GetWorkspace().pdf("D2hhhhBkg")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m2,LineColor(kBlack),LineWidth(1));
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m2->Draw();

  std::cout<<"nentries"<<tree->GetEntries()<<std::endl;
  c1->Draw();
  c1->Print("../img/massFit_invertedPID.eps");

  deltaM_xi_bkg.setVal(myModel->GetWorkspace().var("deltaM_xi_bkg")->getValV() );
  deltaM_lambda_bkg.setVal(myModel->GetWorkspace().var("deltaM_lambda_bkg")->getValV() );
  deltaM_gamma_bkg.setVal(myModel->GetWorkspace().var("deltaM_gamma_bkg")->getValV() );
  deltaM_delta_bkg.setVal(myModel->GetWorkspace().var("deltaM_delta_bkg")->getValV() );
  D0_M_xi_bkg.setVal(myModel->GetWorkspace().var("D0_M_xi_bkg")->getValV() );
  D0_M_lambda_bkg.setVal(myModel->GetWorkspace().var("D0_M_lambda_bkg")->getValV() );
  D0_M_gamma_bkg.setVal(myModel->GetWorkspace().var("D0_M_gamma_bkg")->getValV() );
  D0_M_delta_bkg.setVal(myModel->GetWorkspace().var("D0_M_delta_bkg")->getValV() );

  if(fixShape){

    deltaM_xi_bkg.setConstant();
    deltaM_lambda_bkg.setConstant();
    deltaM_gamma_bkg.setConstant();
    deltaM_delta_bkg.setConstant();
    //D0_M_xi_bkg.setConstant();                                                                                                                                                           
    //D0_M_lambda_bkg.setConstant();                                                                                                  
                                                    
    //D0_M_gamma_bkg.setConstant();                                                                                                                                                        
    //D0_M_delta_bkg.setConstant();                                                                                                                                                         
  }
  file->Close();
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

double D2hhmumuFitter::getMisIDbkgExp(TString cut, TString namePlot){


  TFile* file;
  file= new TFile(pathToNormData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  TTree* cutTree = tree->CopyTree(cut);

   
  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M,deltaM));
  cout <<cutTree->GetEntries() <<endl;

  ///Fit                                                                                                                                                                                   
  D2hhmumuModel* myModel = initializeModel();
  std::string components="Signal CombinatoricBkg D2hhhhBkg";
  RooAbsPdf* finalPDF = myModel->Model(components);

  RooFitResult *result;
  result = finalPDF->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));                                                                                                                   

  TCanvas* c1= new TCanvas("");
  c1->Divide(2);
  c1->cd(1);

  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("D2hhhhBkg")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1));
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel->GetWorkspace().pdf("D2hhhhBkg")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m2,LineColor(kBlack),LineWidth(1));
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m2->Draw();
  c1->Draw();
  c1->Print("../img/selectionOptimization/"+namePlot+"midID.eps");
  delete c1;

  double nhhhhBkg=myModel->GetWorkspace().var("nD2hhhhBkg")->getValV();

  file->Close();
  
  return nhhhhBkg;


}

double D2hhmumuFitter::getCombBkg(TString cut,TString namePlot){


  TFile* file;
  //file= new TFile(pathToSidebandData,"OPEN");
  file= new TFile(pathToSignalData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  TTree* cutTree = tree->CopyTree(cut);

  D0_M.setRange(1880,1940);
  deltaM.setRange(146.5,160);
  D0_M_chebyB.setVal(0);
  D0_M_chebyB.setConstant();

  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M,deltaM));
  cout <<cutTree->GetEntries() <<endl;

  ///Fit               

  D2hhmumuModel* myModel = initializeModel();
  //std::string components="Signal CombinatoricBkg  D2hhhhBkg ";
  
  std::string components="CombinatoricBkg ";                                                                                                         
  RooAbsPdf* finalPDF = myModel->Model(components);

  RooFitResult *result;
  result = finalPDF->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));                                                                                                                    
  TCanvas* c1= new TCanvas("");
  c1->Divide(2);
  c1->cd(1);

  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1));
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m2,LineColor(kBlack),LineWidth(1));
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m2->Draw();
  c1->Draw();
  c1->Print("../img/selectionOptimization/"+namePlot+"combBkg.eps");
  delete c1;


  //extrapolate background to signal region 
  D0_M.setRange(1800,1940);
  deltaM.setRange(143,160);
  RooRealVar dM_threshold("dM_threshold","deltaM_threshold",139.57018);
  RooRealVar dM_alpha("dM_alpha","deltaM_alpha",myModel->GetWorkspace().var("deltaM_alpha")->getValV());
  RooRealVar M_chebyA("M_chebyA","D0_M_chebyA",myModel->GetWorkspace().var("D0_M_chebyA")->getValV());
  RooRealVar M_chebyB("M_chebyB","D0_M_chebyB",myModel->GetWorkspace().var("D0_M_chebyB")->getValV());
  RooChebychev CombinatoricBkgM("CombinatoricBkgM", "Combinatoric Background (M)", D0_M, RooArgList(M_chebyA,M_chebyB));
  RooThreshold CombinatoricBkgDm("CombinatoricBkgDm", "Combinatoric Background (#Deltam)", deltaM, dM_threshold, dM_alpha);
  RooProdPdf CombinatoricBkg("CombinatoricBkg", "Combinatoric Background", RooArgList(CombinatoricBkgM, CombinatoricBkgDm));

  D0_M.setRange("signal",1840,1880);
  deltaM.setRange("signal",144,147);
  RooArgSet variables(D0_M,deltaM);
  RooAbsReal* fr_sig = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("signal"));  
  
  D0_M.setRange("sideband",1880,1940);
  deltaM.setRange("sideband",146.5,160);
  RooAbsReal* fr_sideband = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("sideband"));

  std::cout<<myModel->GetWorkspace().var("nCombinatoricBkg")->getValV()<<"  "<< fr_sig->getVal() << "  "<< fr_sideband->getVal()<<std::endl;                             
  std::cout<<myModel->GetWorkspace().var("nCombinatoricBkg")->getValV() *  fr_sig->getVal() /fr_sideband->getVal()<<std::endl;                                                                                     
  D0_M.setRange(1800,1940);//set back to nominal range
  deltaM.setRange(142,149);
  
  //double nComb=myModel->GetWorkspace().var("nCombinatoricBkg")->getValV();
  double nComb = myModel->GetWorkspace().var("nCombinatoricBkg")->getValV() *  fr_sig->getVal() /fr_sideband->getVal();

  file->Close();

  // return myModel->GetWorkspace().var("nCombinatoricBkg")->getValV()+myModel->GetWorkspace().var("nRandomPionBkg")->getValV();
  return nComb;
 

}



D2hhmumuModel* D2hhmumuFitter::initializeModel(){

  
  //funcition is called by the constructor and initializes a model defined in D2hhmumuModel.h
  // all the variables used in the fit are also initialized by the constructor and can be set indivudially
  //for all channels with the functions set<XX>mumuStartParameters(), XX=Kpi , KK, pipi

  std::cout<<"building a model for D2hhmumu dacays..."<<std::endl;


   D2hhmumuModel* myModel= new D2hhmumuModel();

   myModel->Signal(D0_M,deltaM,
                 D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
                 deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta
                 );


  //background models

  myModel->CombinatoricBackground(D0_M,deltaM,
				 D0_M_chebyA,D0_M_chebyB,
				 deltaM_threshold,deltaM_alpha
				 );


  //randompi
  myModel->RandomPionBackground(D0_M,deltaM,
			       D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
			       deltaM_threshold,deltaM_alpha
                               );

  //peaking , in mD0 just a Gaussian... to be improved!

  /*
  myModel->D2hhhhBackground(D0_M,deltaM,
			   mean1,sigma1,
			   deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta
			   );
  */

  myModel->D2hhhhBackground(D0_M,deltaM,
			    D0_M_xi_bkg,D0_M_lambda_bkg,D0_M_gamma_bkg,D0_M_delta_bkg,
			    deltaM_xi_bkg,deltaM_lambda_bkg,deltaM_gamma_bkg,deltaM_delta_bkg
			    );


  myModel->D2hhhhRandomPionBackground(D0_M,deltaM,
                           mean1,sigma1,
			   deltaM_threshold,deltaM_alpha
                           );



  std::cout<<"...done"<<std::endl;
  return myModel;

}

void D2hhmumuFitter::toyStudy() {

  D2hhmumuModel* myModel=initializeModel();
  std::string components="Signal CombinatoricBkg RandomPionBkg D2hhhhBkg D2hhhhRandomPionBkg";
  RooAbsPdf* finalPDF = myModel->Model(components);

  RooMCStudy* mcstudy = new RooMCStudy(*finalPDF,RooArgSet(D0_M,deltaM),Binned(kTRUE),Silence(),Extended(),
				       FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;


  mcstudy->generateAndFit(10);

  // E x p l o r e   r e s u l t s   o f   s t u d y 
  // ------------------------------------------------

  // Make plots of the distributions of mean, the error on mean and the pull of mean
  //RooPlot* frame1 = mcstudy->plotParam(mean,Bins(40)) ;
  //RooPlot* frame2 = mcstudy->plotError(mean,Bins(40)) ;
  //RooPlot* frame3 = mcstudy->plotPull(mean,Bins(40),FitGauss(kTRUE)) ;

  // Plot distribution of minimized likelihood
  RooPlot* frame4 = mcstudy->plotNLL(Bins(40)) ;

  // Make some histograms from the parameter dataset
  //TH1* hh_cor_a0_s1f = mcstudy->fitParDataSet().createHistogram("hh",a1,YVar(sig1frac)) ;
  //TH1* hh_cor_a0_a1  = mcstudy->fitParDataSet().createHistogram("hh",a0,YVar(a1)) ;

  // Access some of the saved fit results from individual toys
  TH2* corrHist000 = mcstudy->fitResult(0)->correlationHist("c000") ;
  TH2* corrHist127 = mcstudy->fitResult(3)->correlationHist("c127") ;
  TH2* corrHist953 = mcstudy->fitResult(6)->correlationHist("c953") ;

  // Draw all plots on a canvas
  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  TCanvas* c = new TCanvas("rf801_mcstudy","rf801_mcstudy",900,900) ;
  c->Divide(3,3) ;
  // c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
  //c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
  //c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;
  c->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;
  //c->cd(5) ; gPad->SetLeftMargin(0.15) ; hh_cor_a0_s1f->GetYaxis()->SetTitleOffset(1.4) ; hh_cor_a0_s1f->Draw("box") ;
  //c->cd(6) ; gPad->SetLeftMargin(0.15) ; hh_cor_a0_a1->GetYaxis()->SetTitleOffset(1.4) ; hh_cor_a0_a1->Draw("box") ;
  c->cd(7) ; gPad->SetLeftMargin(0.15) ; corrHist000->GetYaxis()->SetTitleOffset(1.4) ; corrHist000->Draw("colz") ;
  c->cd(8) ; gPad->SetLeftMargin(0.15) ; corrHist127->GetYaxis()->SetTitleOffset(1.4) ; corrHist127->Draw("colz") ;
  c->cd(9) ; gPad->SetLeftMargin(0.15) ; corrHist953->GetYaxis()->SetTitleOffset(1.4) ; corrHist953->Draw("colz") ; 

  c->Draw();
  c->Print("../img/toyStudy.eps");

}
  
