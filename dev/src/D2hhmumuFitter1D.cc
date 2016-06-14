#include "D2hhmumuFitter1D.h" 
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooStats/ModelConfig.h"
#include "RooUnblindPrecision.h"
#include "RooUnblindOffset.h"
#include "RooCategory.h"
#include "TBox.h"
#include "TLatex.h"
#include "RooRandom.h"
//#include "StandardHypoTestInvDemo.C"

using namespace std;
using namespace RooFit ;

//this class is is a fitter for Dst->D(hhmumu)pi decays.
//
//
//

D2hhmumuFitter1D::D2hhmumuFitter1D():


  EffRatio("EffRatio","Efficieny ratio",1,0,3),
  nNorm("nNorm","number events normalisation channel",2.7421e+03,100,10000),
  BFnorm("BFnorm","normalizaiton mode Branching fraction",1,0,1),
 
  //signal                                                                                                                                           
  D0_M_xi("D0_M_xi","D0_M_xi",1.86560e+03,1860,1870),
  D0_M_lambda("D0_M_lambda","D0_M_lambda",8.54839e+00,0.1,20),
  D0_M_gamma("D0_M_gamma","D0_M_gamma",-5.42758e-02,-2,40),
  D0_M_delta("D0_M_delta","D0_M_delta",6.09288e-01,0.,10),

  ResolutionScale("ResolutionScale","ResolutionScale",1,0.9,1.1),
  globalShift("globalShift","globalShift",1,0.9,1.1),

  D0_M_mean("D0_M_mean","D0_M_mean",1.86560e+03,1860,1870),
  D0_M_sigma("D0_M_sigma","D0_M_sigma",8.54839e+00,0.1,20),
  D0_M_alphaR("D0_M_alphaR","D0_M_alphaR",-1.,-5.,0.),
  D0_M_alphaL("D0_M_alphaL","D0_M_alphaL",-1.,-5.,-0.5),
  D0_M_nL("D0_M_nL","D0_M_nL",8e-01,-10.,20.),
  D0_M_nR("D0_M_nR","D0_M_nR",6.09288e-01,-10.,20.),

  D0_M_mean_bkg("D0_M_mean_bkg","D0_M_mean_bkg",1.84560e+03,1835,1865),
  D0_M_sigma_bkg("D0_M_sigma_bkg","D0_M_sigma_bkg",8.054839e+00,5.,15),
  D0_M_alphaR_bkg("D0_M_alphaR_bkg","D0_M_alphaR_bkg",-2.,-5.,0.),
  D0_M_alphaL_bkg("D0_M_alphaL_bkg","D0_M_alphaL_bkg",2,0,5),//,0.,3.),
  D0_M_nL_bkg("D0_M_nL_bkg","D0_M_nL_bkg",1,0.,40.),
  D0_M_nR_bkg("D0_M_nR_bkg","D0_M_nR_bkg",1,0.,10.), 

  //purely combinatorial background                                                                                                                 
  //D0_M_chebyA("D0_M_chebyA","D0_M_chebyA",-3.5906e-02,-1,1),
  D0_M_chebyA("D0_M_chebyA","D0_M_chebyA",0,-1,1),
  D0_M_chebyB("D0_M_chebyB","D0_M_chebyB",-1.7004e-02,-1,1),
  D0_M_chebyC("D0_M_chebyC","D0_M_chebyC",-1.7882e-02,-1,1),

  //peaking in mD0 and dM D2Kpipipi                                                                                                                       
  D0_M_xi_bkg_norm("D0_M_xi_bkg_norm","D0_M_xi_bkg_norm",1.8312e+03,1.8012e+03,1.912e+03),
  D0_M_lambda_bkg_norm("D0_M_lambda_bkg_norm","D0_M_lambda_bkg_norm",4.0000e+01,1.0000e+01,10.0000e+01),
  D0_M_gamma_bkg_norm("D0_M_gamma_bkg_norm","D0_M_gamma_bkg_norm",-5.8442e-01,-10.8442e-01,-1.8442e-01),
  D0_M_delta_bkg_norm("D0_M_delta_bkg_norm","D0_M_delta_bkg_norm",7.5826e-01,0,1.5826),

  //peaking in mD0 and dM signal shape                                                                                                                       
  D0_M_xi_bkg("D0_M_xi_bkg","D0_M_xi_bkg",1.8312e+03,1.8012e+03,1.912e+03),
  D0_M_lambda_bkg("D0_M_lambda_bkg","D0_M_lambda_bkg",4.0000e+01,1.0000e+01,10.0000e+01),
  D0_M_gamma_bkg("D0_M_gamma_bkg","D0_M_gamma_bkg",-5.8442e-01,-10.8442e-01,-1.8442e-01),
  D0_M_delta_bkg("D0_M_delta_bkg","D0_M_delta_bkg",7.5826e-01,0,1.5826),

  Nsig_blinding("Nsig_blinding","blinding offset for signal yield",0),
  BFsig_blinding("BFsig_blinding","blinding offset for BF measurement",0)
{

  //constructor

  pathToSignalMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root";
  //pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root",
  pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root";                                                      
  pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root";
  pathToInvData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/inverted_D2Kpimumu_BDT_selected.root";
  pathToSidebandData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT.root";
  pathToKpipipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_D2KKmumuBDT.root"; 
  pathToKKpipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKpipi_D2KKmumuBDT.root"; 
  pathToKpipipiHistoData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/zoz-5000.root";

  //blinding part 

  generator= new TRandom3(13031989);
  Nsig_offset = generator->Rndm()*1000;
  BFsig_offset =generator->Rndm()*1e-2; 

  Nsig_blinding.setVal(Nsig_offset);
  BFsig_blinding.setVal(BFsig_offset);

}

void D2hhmumuFitter1D::setPathToSignalMC(TString path){
  pathToSignalMC = path;
}
void D2hhmumuFitter1D::setPathToSignalData(TString path){
  pathToSignalData = path;
}
void D2hhmumuFitter1D::setPathToNormData(TString path){
  pathToNormData = path;
}
void D2hhmumuFitter1D::setPathToInvData(TString path){
  pathToInvData = path;
}
void D2hhmumuFitter1D::setPathToSidebandData(TString path){
  pathToSidebandData = path;
}
void D2hhmumuFitter1D::setPathToKpipipiHistoData(TString path){
  pathToKpipipiHistoData = path;
}
void D2hhmumuFitter1D::setPathToKpipipiData(TString path){
  pathToKpipipiData = path;
}
void D2hhmumuFitter1D::setPathToKKpipiData(TString path){
  pathToKKpipiData = path;
}
D2hhmumuFitter1D::~D2hhmumuFitter1D(){}


void D2hhmumuFitter1D::setKpimumuStartParameters(){                                                                                                                             
                                                 
  //later like this...                                                                                                                                                         
  //deltaM_xi.SetValue()..
  
}

void D2hhmumuFitter1D::setKKmumuStartParameters(){

  //later like this...                                                                                                                                                         
  //deltaM_xi.SetValue()..                                                                                                                                                      

}

void D2hhmumuFitter1D::setpipimumuStartParameters(){

  //later like this...                                                                                                                                                         
  //deltaM_xi.SetValue()..                                                                                                                                                      

}

void D2hhmumuFitter1D::fit_MC(TString cut="",bool fixShape = true,TString namePlot=""){

  //observables                                                                                                                                                                            
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1820., 1900.,"MeV");
  RooRealVar nSignal("nSignal","#signal events ",100000,0,100000,"MeV");

  ///create Model with desired components                                                                                                                                                  

  ResolutionScale.setRange(1,1);
  globalShift.setRange(1,1);

  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  initializeModel(myModel,D0_M); //initializes the model with parameters                                                                                                                    
  std::string components="Signal";
  myModel->Model(components); //acutally build the model with default coefficients                                                                                                          
  RooWorkspace m_ws = myModel->GetWorkspace(); // get workspace with  all PDFs and coefficients                                                                                            
  RooAbsPdf* finalPDF = m_ws.pdf("D2hhmumuModel");

  //get data to fit                                                                                                                                                                         
  TFile* file;
  file= new TFile(pathToSignalMC,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);
  tree->SetBranchStatus("nTracks",1);
  //apply cuts if needed                                                                                                                                                                    
  TTree* cutTree = tree->CopyTree(cut+"&& deltaM>144.5 && deltaM <146.5");
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M));
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));
  RooDataHist* data_binned = data_small->binnedClone();
  m_ws.import(*data);
  //do the fit                          
  RooFitResult *result;
  result = finalPDF->fitTo(*data,Save(kTRUE),NumCPU(3));
  m_ws.import(*result);

  ///Plot                                                                                                                                                                                   
  ///----------                                                                                                                               
                                                                                                                                                                                           
  TCanvas* c1= new TCanvas("");
  c1->Divide(2,2);
  c1->cd(1);

  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));

  finalPDF->plotOn(frame_m,LineColor(kRed),LineWidth(1));
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  //finalPDF->paramOn(frame_m);
  frame_m->Draw();

  c1->cd(2);
  RooPlot* frame_m2 = D0_M.frame(Title("fit parameters")) ;
  finalPDF->paramOn(frame_m2,Layout(.2,.80,.99));
  frame_m2->Draw();

  RooHist* hpull = frame_m->pullHist() ;
  RooPlot* frame_m3 = D0_M.frame(Title("Pull Distribution")) ;
  frame_m3->addPlotable(hpull,"P") ;
  c1->cd(3);
  frame_m3->Draw();

  c1->Draw();
  //c1->Print("../img/massFit1D_MC.eps");
  c1->Print("../D2KKmumu/img/MC_fit_"+cut+".eps");
  c1->Print(namePlot);

  //std::cout<<"nSig "<<(m_ws.var("nSignal")->getValV())<<std::endl;
  D0_M_xi.setVal(m_ws.var("D0_M_xi")->getValV() );
  D0_M_lambda.setVal(m_ws.var("D0_M_lambda")->getValV() );
  D0_M_gamma.setVal(m_ws.var("D0_M_gamma")->getValV() );
  D0_M_delta.setVal(m_ws.var("D0_M_delta")->getValV() );

  // write the workspace in the file                                                                                                                                                       
  //  TString fileName = "MC_model1D.root";
  // m_ws.writeToFile(fileName,true);
  m_ws.writeToFile("../D2KKmumu/img/MC_WS_"+cut+".root",true);
  //  cout << "model written to file " << fileName << endl;

  if(fixShape){
    D0_M_xi.setConstant();                                                                                                                                                                
    D0_M_lambda.setConstant();                                                                                                                                                            
    D0_M_gamma.setConstant();
    D0_M_delta.setConstant();
  }
  
  delete tree;
  delete cutTree;

  file->Close();

  delete file;
  
}


void D2hhmumuFitter1D::fit_Kpipipi_misID(TString cut="",bool fixShape=false,TString namePlot=""){


  RooRealVar D0_M("misID_mD_OS", "m(h h #mu #mu)", 1820., 1940.,"MeV");
                                       
  TFile* file;
  file= new TFile(pathToKpipipiData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);
  tree->SetBranchStatus("misID_mD_OS",1);
  tree->SetBranchStatus("mu0_PT",1);
  tree->SetBranchStatus("mu1_PT",1);

  //tree->SetBranchStatus("nTracks",1);

  D0_M_chebyB.setVal(0);
  D0_M_chebyB.setConstant();
  D0_M.setRange(1780,1900);


  TTree* cutTree = tree->CopyTree(cut+"&& deltaM > 144.5 && deltaM <146.5");
  RooArgList list =  RooArgList( D0_M );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M));
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));
  RooDataHist* data_binned = data_small->binnedClone();

  ///Fit                                                                                                          

  ///create Model with desired components                                                                                                                                    
  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  initializeNormalizationModel(myModel,D0_M);
  std::string components="D2hhhhBkg CombinatoricBkg";
  RooAbsPdf* finalPDF = myModel->Model(components);
  RooFitResult *result;
  result = finalPDF->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));

  cout << "result is --------------- "<<endl;
  result->Print();

  ///Plot                                                                                                                                                   ///----------                                                                                                                                                                
  TCanvas* c1= new TCanvas("");
  c1->Divide(2,2);
  c1->cd(1);
  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("D2hhhhBkg")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1));
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  //finalPDF->paramOn(frame_m);
  frame_m->Draw();
  
  c1->cd(2);
  RooPlot* frame_m2 = D0_M.frame(Title("fit parameters")) ;
  finalPDF->paramOn(frame_m2,Layout(.2,.80,.99));
  frame_m2->Draw();

  RooHist* hpull = frame_m->pullHist() ;
  RooPlot* frame_m3 = D0_M.frame(Title("Pull Distribution")) ;
  frame_m3->addPlotable(hpull,"P") ;
  c1->cd(3);
  frame_m3->Draw();

  std::cout<<"nentries"<<tree->GetEntries()<<std::endl;
  c1->Draw();
  c1->Print("../D2KKmumu/img/misID_fit_"+cut+".eps");
  c1->Print(namePlot);

   D0_M_xi_bkg_norm.setVal(myModel->GetWorkspace().var("D0_M_xi_bkg_norm")->getValV() );
   D0_M_lambda_bkg_norm.setVal(myModel->GetWorkspace().var("D0_M_lambda_bkg_norm")->getValV() );
   D0_M_gamma_bkg_norm.setVal(myModel->GetWorkspace().var("D0_M_gamma_bkg_norm")->getValV() );
   D0_M_delta_bkg_norm.setVal(myModel->GetWorkspace().var("D0_M_delta_bkg_norm")->getValV() );

  if(fixShape){

     D0_M_xi_bkg_norm.setConstant();                                                                                                                     
     D0_M_lambda_bkg_norm.setConstant();                                                                                                  
     D0_M_gamma_bkg_norm.setConstant();                                                                                                          
     D0_M_delta_bkg_norm.setConstant();                                                                                                                                                       
  
  }
  delete tree;
  delete cutTree;
  file->Close();
  delete file;


}
void D2hhmumuFitter1D::fit_HHpipi_misID(TString cut="",bool fixShape=false,TString namePlot=""){


  RooRealVar D0_M("misID_mD_OS", "m(h h #mu #mu)", 1820., 1940.,"MeV");
                                       
  TFile* file;
  file= new TFile(pathToKKpipiData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);
  tree->SetBranchStatus("misID_mD_OS",1);
  tree->SetBranchStatus("mu0_PT",1);
  tree->SetBranchStatus("mu1_PT",1);

  //tree->SetBranchStatus("nTracks",1);

  D0_M_chebyB.setVal(0);
  D0_M_chebyB.setConstant();
  D0_M.setRange(1780,1900);


  TTree* cutTree = tree->CopyTree(cut+"&& deltaM > 144.5 && deltaM <146.5");
  RooArgList list =  RooArgList( D0_M );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M));
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));

  ///Fit                                                                                                          

  ///create Model with desired components                                                                                                                                    
  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  initializeModel(myModel,D0_M);
  std::string components="D2hhhhBkg CombinatoricBkg";
  RooAbsPdf* finalPDF = myModel->Model(components);
  RooFitResult *result;
  result = finalPDF->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));

  cout << "result is --------------- "<<endl;
  result->Print();

  ///Plot                                                                                                                                                   ///----------                                                                                                                                                                
  TCanvas* c1= new TCanvas("");
  c1->Divide(2,2);
  c1->cd(1);
  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("D2hhhhBkg")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1));
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  //finalPDF->paramOn(frame_m);
  frame_m->Draw();
  
  c1->cd(2);
  RooPlot* frame_m2 = D0_M.frame(Title("fit parameters")) ;
  finalPDF->paramOn(frame_m2,Layout(.2,.80,.99));
  frame_m2->Draw();

  RooHist* hpull = frame_m->pullHist() ;
  RooPlot* frame_m3 = D0_M.frame(Title("Pull Distribution")) ;
  frame_m3->addPlotable(hpull,"P") ;
  c1->cd(3);
  frame_m3->Draw();

  std::cout<<"nentries"<<tree->GetEntries()<<std::endl;
  c1->Draw();
  c1->Print("../D2KKmumu/img/misID_fit_"+cut+".eps");
  c1->Print(namePlot);

   D0_M_xi_bkg.setVal(myModel->GetWorkspace().var("D0_M_xi_bkg")->getValV() );
   D0_M_lambda_bkg.setVal(myModel->GetWorkspace().var("D0_M_lambda_bkg")->getValV() );
   D0_M_gamma_bkg.setVal(myModel->GetWorkspace().var("D0_M_gamma_bkg")->getValV() );
   D0_M_delta_bkg.setVal(myModel->GetWorkspace().var("D0_M_delta_bkg")->getValV() );

  if(fixShape){

     D0_M_xi_bkg.setConstant();                                                                                                                     
     D0_M_lambda_bkg.setConstant();                                                                                                  
     D0_M_gamma_bkg.setConstant();                                                                                                          
     D0_M_delta_bkg.setConstant();                                                                                                                                                       
  
  }
  delete tree;
  delete cutTree;
  file->Close();
  delete file;


}
void D2hhmumuFitter1D::fit_Kpipipi_misID_fromHistogramm(TString cut="",bool fixShape=false,TString namePlot=""){


  RooRealVar D0_M("misID_mD_OS", "m(h h #mu #mu)", 1720., 1940.,"MeV");
                                       
  TFile* file;
  file= new TFile(pathToKpipipiHistoData,"OPEN");
  TH1* temp = (TH1D*)file->Get("float_model_aux__D_MM");
  TH1* toyData = new TH1D("toyData","toyData",500,1720,1940);
  
  for(int i =0;i<10000;++i) {
    toyData->Fill(temp->GetRandom());
  }

  RooDataHist histo("histo","histo",D0_M,Import(*toyData));

  D0_M_chebyB.setVal(0);
  D0_M_chebyB.setConstant();
  //D0_M.setRange(1780,1900);

  ///Fit                                                                                                          

  ///create Model with desired components                                                                                                                                    
  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  initializeNormalizationModel(myModel,D0_M);
  std::string components="D2hhhhBkg ";
  RooAbsPdf* finalPDF = myModel->Model(components);
  RooFitResult *result;
  result = finalPDF->fitTo(histo,Save(kTRUE),Extended(kTRUE),NumCPU(3));

  cout << "result is --------------- "<<endl;
  result->Print();

  ///Plot                                                                                                                                                   ///----------                                                                                                                                                                
  TCanvas* c1= new TCanvas("");
  c1->Divide(2,2);
  c1->cd(1);
  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  histo.plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("D2hhhhBkg")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1));
  histo.plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  //finalPDF->paramOn(frame_m);
  frame_m->Draw();
  
  c1->cd(2);
  RooPlot* frame_m2 = D0_M.frame(Title("fit parameters")) ;
  finalPDF->paramOn(frame_m2,Layout(.2,.80,.99));
  frame_m2->Draw();

  RooHist* hpull = frame_m->pullHist() ;
  RooPlot* frame_m3 = D0_M.frame(Title("Pull Distribution")) ;
  frame_m3->addPlotable(hpull,"P") ;
  c1->cd(3);
  frame_m3->Draw();

  c1->Draw();
  c1->Print("../D2KKmumu/img/template_misID_fit_"+cut+".eps");
  c1->Print(namePlot);

   D0_M_xi_bkg_norm.setVal(myModel->GetWorkspace().var("D0_M_xi_bkg_norm")->getValV() );
   D0_M_lambda_bkg_norm.setVal(myModel->GetWorkspace().var("D0_M_lambda_bkg_norm")->getValV() );
   D0_M_gamma_bkg_norm.setVal(myModel->GetWorkspace().var("D0_M_gamma_bkg_norm")->getValV() );
   D0_M_delta_bkg_norm.setVal(myModel->GetWorkspace().var("D0_M_delta_bkg_norm")->getValV() );

  if(fixShape){

     D0_M_xi_bkg_norm.setConstant();                                                                                                                     
     D0_M_lambda_bkg_norm.setConstant();                                                                                                  
     D0_M_gamma_bkg_norm.setConstant();                                                                                                          
     D0_M_delta_bkg_norm.setConstant();                                                                                                                                                         
  }
  file->Close();
  delete file;
}

double D2hhmumuFitter1D::fit_normalization_Data(TString cut="",TString namePlot=""){
 
  //observables                                                                                                                                                                             
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1770., 1940.,"MeV");
  TString dmRange = "&& deltaM>144.5 && deltaM<146.5";

  //allow for MC Data differences in resolution and global mass shift
  ResolutionScale.setRange(0.9,1.2); 
  globalShift.setRange(0.9,1.2);
  //D0_M_chebyA.setRange(0,2);  

  ///create Model with desired components                                                                                                                                                   
  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  RooWorkspace m_ws = initializeNormalizationModel(myModel,D0_M); //get a workspace with all PDFs       
  //set the yields and add the to ws
  RooRealVar nSignal("nSignal","number of signal events",2000,0,5000);                              
  RooRealVar nCombinatoricBkg("nCombinatoricBkg","number of combinatorial bkg events",100,0,10000);                              
  RooRealVar nD2hhhhBkg("nD2hhhhBkg","number of misidentified D2hhhh background",500,0,2000);                               

  m_ws.import(nSignal);
  m_ws.import(nCombinatoricBkg);
  m_ws.import(nD2hhhhBkg);

  //create the PDF
  m_ws.factory("SUM::tot(nSignal*Signal,nCombinatoricBkg*CombinatoricBkg,nD2hhhhBkg*D2hhhhBkg)");
  RooAbsPdf* finalPDF = m_ws.pdf("tot");

  //get data to fit
  TFile* file;
  file= new TFile(pathToNormData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);  
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);
  tree->SetBranchStatus("nTracks",1);

  //apply cuts if needed
  TTree* cutTree = tree->CopyTree(cut+dmRange);
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M));
  m_ws.import(*data);

  //do the fit 
  RooFitResult *result;
  result = finalPDF->fitTo(*data,Save(kTRUE),NumCPU(3));
  m_ws.import(*result);
  cout << "result is --------------- "<<endl;
  result->Print();

  // write the workspace in the file
  TString fileName = "normData_model.root";
  nNorm.setVal(m_ws.var("nSignal")->getVal());
  nNorm.setError(m_ws.var("nSignal")->getError());

  m_ws.writeToFile(fileName,true);
  cout << "model written to file " << fileName << endl;
  
  ///Plot                                                                                                                                                                 ///----------                                                                                                                                                                
  TCanvas* c1= new TCanvas("canvas");
  c1->Divide(2,2);
  c1->cd(1);

  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  m_ws.pdf("tot")->plotOn(frame_m,Components(*m_ws.pdf("Signal")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  m_ws.pdf("tot")->plotOn(frame_m,Components(*m_ws.pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  m_ws.pdf("tot")->plotOn(frame_m,Components(*m_ws.pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1)); 
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  //finalPDF->paramOn(frame_m);
  frame_m->Draw();

  c1->cd(2);
  RooPlot* frame_m2 = D0_M.frame(Title("fit parameters")) ;
  finalPDF->paramOn(frame_m2,Layout(.2,.80,.99));
  frame_m2->Draw();
 
  RooHist* hpull = frame_m->pullHist() ;
  RooPlot* frame_m3 = D0_M.frame(Title("Pull Distribution")) ;
  frame_m3->addPlotable(hpull,"P") ;
  c1->cd(3);
  frame_m3->Draw();
 
  std::cout<<"nentries"<<tree->GetEntries()<<std::endl;	
  c1->Draw();
  //c1->Print("../img/massFit1D_Normalization.eps");
  c1->Print("../D2KKmumu/img/normalization_fit_"+cut+".eps");
  c1->Print(namePlot);


  //fix scale factor for MC - data differences
  ResolutionScale.setVal(myModel->GetWorkspace().var("ResolutionScale")->getValV() );
  globalShift.setVal(myModel->GetWorkspace().var("globalShift")->getValV() );  
  ResolutionScale.setConstant();
  globalShift.setConstant();
  

  delete tree;
  delete cutTree;

  delete file;
 
  return m_ws.var("nSignal")->getVal();

}

void D2hhmumuFitter1D::fit_Data(TString cut="",TString namePlot=""){

  /*
 D0_M_delta_bkg.setVal(  8.6917e-01 );
 D0_M_gamma_bkg.setVal( -5.6508e-01 );
 D0_M_lambda_bkg.setVal( 1.7008e+01 );
 D0_M_xi_bkg.setVal(  1.8431e+03 );
 D0_M_delta_bkg.setConstant();
 D0_M_gamma_bkg.setConstant();
 D0_M_lambda_bkg.setConstant();
 D0_M_xi_bkg.setConstant();
 nNorm.setVal(2.7421e+03);
 nNorm.setConstant();
 D0_M_chebyB.setVal(0);
 D0_M_chebyB.setConstant();
  */

  //observables                                                                                                                                                                             
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1780., 1950.,"MeV");
  TString dmRange = "&& deltaM>144.5 && deltaM<146.5";

  //get data to fit                                                                                                                                                                          
  TFile* file;
  file= new TFile(pathToSignalData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);
  tree->SetBranchStatus("nTracks",1);
  //apply cuts if needed                                                                                                                                                                    
  TTree* cutTree = tree->CopyTree(cut+dmRange);
  RooArgList list =  RooArgList( D0_M );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M));

  ///////////////////////
  //                   //
  //   fit the yield   //
  //                   //
  ///////////////////////


  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  RooWorkspace m_ws = initializeModel(myModel,D0_M); //get a workspace with all PDFs                                                                                                        
  m_ws.SetName("m_ws");
  m_ws.import(*data);

  RooRealVar nSignal("nSignal","number of signal events",50,0,200);                                                                                       
  RooRealVar nCombinatoricBkg("nCombinatoricBkg","number of combinatorial bkg events",10,0,300);                                                   
  RooRealVar nD2hhhhBkg("nD2hhhhBkg","number of misidentified D2hhhh background",5,0,20);                                              
  RooRealVar nSignal_blind("nSignal_blind","number of signal events blind",nSignal.getVal()+Nsig_blinding.getVal(),nSignal.getMin()+Nsig_blinding.getVal(),nSignal.getMax()+Nsig_blinding.getVal());

  m_ws.import(nCombinatoricBkg);
  m_ws.import(nD2hhhhBkg);
  m_ws.import(nSignal_blind);                                                                                                                                   
  m_ws.import(Nsig_blinding);
  m_ws.factory("expr::yield('nSignal_blind-Nsig_blinding',nSignal_blind,Nsig_blinding)");                                                                                             
  m_ws.factory("SUM::tot(yield*Signal,nCombinatoricBkg*CombinatoricBkg,nD2hhhhBkg*D2hhhhBkg)");

  RooAbsPdf* finalPDF = m_ws.pdf("tot");                                                                                                                                                
  m_ws.import(*finalPDF);
 
  RooFitResult *result;                                                                                                                                                                  
  result = finalPDF->fitTo(*data,Save(kTRUE),NumCPU(3));                                                                                                                                 
  m_ws.import(*result);                                                                                                                                                                  
  cout << "result is --------------- "<<endl;                                                                                                                                            
  result->Print();                                                                                                                                                                       

  TCanvas* c1= new TCanvas("canvas");
  c1->Divide(2);
  c1->cd(1);

  D0_M.setRange("lowersideband",1780.,1830.);
  D0_M.setRange("uppersideband",1900.,1950.);
  double sidebands = data->sumEntries("1","lowersideband,uppersideband");
  std::cout<<"sidebands "<< sidebands << std::endl;

  c1->cd(1);
  RooPlot* frame_m= D0_M.frame();                                                                                                                                          
  frame_m->SetTitle("");                                                                                                                                                    
  TBox * blindingbox = new TBox();
  blindingbox->SetFillColor(0);
  data->plotOn(frame_m,Name("data"),CutRange("lowersideband,uppersideband"),MarkerSize(0.5),Binning(50));                
  finalPDF->plotOn(frame_m,Range("lowersideband,uppersideband"),LineColor(kRed),LineWidth(1),Normalization(sidebands,RooAbsReal::NumEvent));                                  
  frame_m->Draw();
  blindingbox->DrawBox(1830.,0.0,1900.,frame_m->GetMaximum()*0.97);

  c1->cd(2);                                                                                                                                    
  RooPlot* frame_m2 = D0_M.frame(Title("fit parameters")) ;
  finalPDF->paramOn(frame_m2,Layout(.2,.80,.99));
  frame_m2->Draw();

  c1->Draw();
  //c1->Print("../D2KKmumu/img/Data_fit_"+cut+".eps");
  c1->Print(namePlot);
  

  ///////////////////////////
  //                       //
  //   FIT BR              //
  //                       //
  ///////////////////////////

                                                                                                                                                                                        
  //redfintion of PDF in terms of BR                                                                                                                                      
  EffRatio.setVal(1); EffRatio.setConstant();                                                                                                              
  std::cout<<"nnorm "<< nNorm.getVal() <<std::endl;                                                                                                         
  nNorm.setVal(nNorm.getVal()*1);                                                                                                                                       
  nNorm.setConstant();                                                                                                         
  BFnorm.setVal(4.17e-6);BFnorm.setConstant();                                                                                                                                            
  RooRealVar BFsig("BFsig","signal Branching fraction",1e-7,1e-9,1e-6);                  
  RooRealVar BFblind("BFblind","BFblind",BFsig.getVal()+BFsig_blinding.getVal(),BFsig.getMin()+BFsig_blinding.getVal(),BFsig.getMax()+BFsig_blinding.getVal());
  
  std::cout<<"BFBLIND"<<BFblind.getVal()<<" "<< BFblind.getMin()  <<"  "<< BFblind.getMax()<<std::endl;

  m_ws.import(BFblind); // add nuissance parameters to the WS                                                                                                                      
  m_ws.import(BFsig);
  m_ws.import(BFsig_blinding);
  m_ws.import(EffRatio);                                                                                                                                                              
  m_ws.import(BFnorm);                                                                                                                                                                
  m_ws.import(nNorm);                                                                                                                                                                     
  
  m_ws.factory("expr::sig_yield('nNorm*((BFblind-BFsig_blinding)/BFnorm)*EffRatio',nNorm,BFblind,BFsig_blinding,BFnorm,EffRatio)");    
  m_ws.factory("SUM::tot_pdf(nCombinatoricBkg*CombinatoricBkg,sig_yield*Signal,nD2hhhhBkg*D2hhhhBkg)");                                               
        
  RooAbsPdf* finalPDF2 = m_ws.pdf("tot_pdf");                        
  RooFitResult *result2;                                                                                                                                                              
  result2 = finalPDF2->fitTo(*data,Save(kTRUE),NumCPU(3));                                                                                                
  result2->Print();                                                                                                                                                                          
}

void D2hhmumuFitter1D::fillModelConfig(TString dataCut,TString nomalizationCut,TString misIDCut, TString fileName){
  
  /*
  D0_M_delta_bkg.setVal(  8.6917e-01 );
  D0_M_gamma_bkg.setVal( -5.6508e-01 );
  D0_M_lambda_bkg.setVal( 1.7008e+01 );
  D0_M_xi_bkg.setVal(  1.8431e+03 );
  D0_M_delta_bkg.setConstant();
  D0_M_gamma_bkg.setConstant();
  D0_M_lambda_bkg.setConstant();
  D0_M_xi_bkg.setConstant();
  nNorm.setVal(2.7421e+03);
  nNorm.setConstant();
  D0_M_chebyB.setVal(0);
  D0_M_chebyB.setConstant();
  */

  //fix the mass shapes and get normalization 
  this->fit_MC(dataCut,true);
  this->fit_Kpipipi_misID(misIDCut,true);
  this->fit_normalization_Data(nomalizationCut);  //neeeded to get nNorm

  this->fit_HHpipi_misID(misIDCut,true); //fix D2HHmumu shape

  std::cout<< "nnorm " << nNorm.getVal()  <<std::endl;

  std::cout << "Filling workpace for cut: "<< dataCut << " to "<< fileName << std::endl;
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1780., 1950.,"MeV");
  TString dmRange = "&& deltaM>144.5 && deltaM<146.5";
 
  TFile* file;
  file= new TFile(pathToSignalData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);

  TTree* cutTree = tree->CopyTree(dataCut+dmRange);
  RooArgList list =  RooArgList( D0_M );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M));

  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  RooWorkspace m_ws = initializeModel(myModel,D0_M); //get a workspace with all PDFs                                                                                                                                                                                                                                                                             
  m_ws.SetName("m_ws");
  m_ws.import(*data);

  RooRealVar nCombinatoricBkg("nCombinatoricBkg","number of combinatorial bkg events",30,-300,300);
  RooRealVar nD2hhhhBkg("nD2hhhhBkg","number of misidentified D2hhhh background",2,-50,50);
  RooRealVar BFsig("BFsig","signal Branching fraction",1e-8,1e-9,1e-6);
  D0_M_chebyA.setVal(0.4);

  // add nuissance parameters to the WS                                                                                                                               
  
  EffRatio.setVal(1);
  BFnorm.setVal(4.17e-6);

  //nNorm.setConstant();
  BFnorm.setConstant();
  EffRatio.setConstant();

  m_ws.import(nCombinatoricBkg);
  m_ws.import(nD2hhhhBkg);
  m_ws.import(BFsig);
  m_ws.import(EffRatio);
  m_ws.import(BFnorm);
  m_ws.import(nNorm);
  
  m_ws.factory("expr::sig_yield('nNorm*BFsig/BFnorm*EffRatio',nNorm,BFsig,BFnorm,EffRatio)");
  m_ws.factory("SUM::tot_pdf(nCombinatoricBkg*CombinatoricBkg,sig_yield*Signal,nD2hhhhBkg*D2hhhhBkg)");
  m_ws.factory("Gaussian::constraintnNorm(nNorm0[1000,40000],nNorm,nNorm_err[1])");
  m_ws.factory("Gaussian::constraintEff(Eff0[0,1.5],EffRatio,EffRatio_err[1])");
  m_ws.factory("Gaussian::constraintBFnorm(BFnorm0[1e-6,1e-5],BFnorm,BFnorm_err[1])");
 
  m_ws.factory("PROD:model(tot_pdf,constraintnNorm,constraintEff,constraintBFnorm)"); 

  //set constraints to nuisance parameteres
  m_ws.var("nNorm0")->setVal(nNorm.getVal());
  m_ws.var("nNorm_err")->setVal(nNorm.getError());
  m_ws.var("nNorm0")->setConstant(true);
  std::cout<<"nNorm0 " << m_ws.var("nNorm0")->getVal() << "+-" << m_ws.var("nNorm_err")->getVal()  <<std::endl; 

  m_ws.var("Eff0")->setVal(EffRatio.getVal());
  m_ws.var("EffRatio_err")->setVal(0.05);
  m_ws.var("Eff0")->setConstant(true);

  m_ws.var("BFnorm0")->setVal(BFnorm.getVal());
  m_ws.var("BFnorm_err")->setVal(0.1*BFnorm.getVal());
  m_ws.var("BFnorm0")->setConstant(true);


  //RooAbsPdf* finalPDF2 = m_ws.pdf("model");
  //RooFitResult *result2;
  //finalPDF2->fitTo(*data,Save(kTRUE),NumCPU(3));
  //result2->Print();    

  m_ws.defineSet("obs","Dst_DTF_D0_M");                                                                                                                                          
  m_ws.defineSet("poi","BFsig");                                                                                                                                                 
  m_ws.defineSet("nuispar","nCombinatoricBkg,nNorm,EffRatio,BFnorm,nD2hhhhBkg,D0_M_chebyA");                                                                                          
  //m_ws.defineSet("nuispar","nCombinatoricBkg,nD2hhhhBkg,D0_M_chebyA");                                                                                          
  //m_ws.defineSet("gObs","nNorm0,eff0,BFnorm0");                                                                                                                                          
  //create the Model config which is passes to limit calculator                                                                                                                           
  m_ws.defineSet("gObs","nNorm0,Eff0,BFnorm0");

  RooStats::ModelConfig*  myModelConfig= new RooStats::ModelConfig("ModelConfig",&m_ws);                                                                                                     
  //myModelConfig->SetWorkspace(m_ws);                                                                                                                                            
  myModelConfig->SetPdf(*(m_ws.pdf("model")));  
  myModelConfig->SetObservables(*m_ws.set("obs"));                                                                                                                                
  myModelConfig->SetParametersOfInterest(*m_ws.set("poi"));                                                                                                                 
  myModelConfig->SetNuisanceParameters(*m_ws.set("nuispar"));                                                                                                    
  myModelConfig->SetSnapshot(*m_ws.set("poi"));                                                                                                                         
  myModelConfig->SetGlobalObservables(*m_ws.set("gObs"));  
 
  m_ws.import(*myModelConfig);             
  
  RooStats::ModelConfig* bModel = (RooStats::ModelConfig*) myModelConfig->Clone();                                                                            
  bModel->SetName("B_only_model");                                                                                                                               
  RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());                                                                                  
  if (!var) exit(1);                                                                                                                                                          
  double oldval = var->getVal();                                                                                                                                   
  var->setVal(0);                                                                                                                                                      
  bModel->SetSnapshot( RooArgSet(*var)  );                                                                                                                                
  var->setVal(oldval);                                                                                                                                                         
  m_ws.import(*bModel);                                                                                                                                                         

  // write the workspace in the file                                                                                                                                                      
  //TString fileName = "Data_model1D.root";                                                                                                                                                
  m_ws.writeToFile(fileName,true);                                                                                                                                             
  cout << "model written to file " << fileName << endl;                                                                                                                                  

}

void D2hhmumuFitter1D::makeToyStudy(TString dataCut,TString nomalizationCut,TString misIDCut,TString targetFile,double nSig_exp, double nCombBkg_exp, double nMisID_exp, bool fixMisID20=false){
   
  //fix the mass shapes and get normalization                                                                                                                                            
  std::cout<<"Beginning toy study: data cut "<< dataCut << " with expected yields: nsig = "<< nSig_exp <<  "ncombbkg = "<< nCombBkg_exp<<  "nmisID = "<< nMisID_exp <<std::endl;   

  TH1* pullSigYield = new TH1D("pullSigYield","pull of signal yield",50,-5,5); 
  TH1* pullMisIDYield = new TH1D("pullmsiIDYield","pull of misID yield",50,-5,5); 
  TH1* pullCombYield = new TH1D("pullCombYield","pull of comb bkg yield",50,-5,5); 
  TH1* SigYield = new TH1D("SigYield"," signal yield",50,nSig_exp - 5*TMath::Sqrt(nSig_exp),nSig_exp + 5*TMath::Sqrt(nSig_exp)); //+-5 sigma window
  TH1* MisIDYield = new TH1D("misIDYield"," misID yield",50,nMisID_exp - 5*TMath::Sqrt(nMisID_exp),nMisID_exp + 5*TMath::Sqrt(nMisID_exp)); 
  TH1* CombYield = new TH1D("CombYield"," comb bkg yield",50,nCombBkg_exp - 5*TMath::Sqrt(nCombBkg_exp),nCombBkg_exp + 5*TMath::Sqrt(nCombBkg_exp)); 
  TH1* errorSigYield = new TH1D("errorSigYield","error of signal yield",50,TMath::Sqrt(nSig_exp) - 5*TMath::Sqrt(TMath::Sqrt(nSig_exp)),TMath::Sqrt(nSig_exp) + 5*TMath::Sqrt(TMath::Sqrt(nSig_exp))); 
TH1* errorMisIDYield = new TH1D("errormsiIDYield","error of misID yield",50,TMath::Sqrt(nMisID_exp) - 5*TMath::Sqrt(TMath::Sqrt(nMisID_exp)),TMath::Sqrt(nMisID_exp) + 5*TMath::Sqrt(TMath::Sqrt(nMisID_exp)));
TH1* errorCombYield = new TH1D("errorCombYield","error of comb bkg yield",50,TMath::Sqrt(nCombBkg_exp) - 5*TMath::Sqrt(TMath::Sqrt(nCombBkg_exp)),TMath::Sqrt(nCombBkg_exp) + 5*TMath::Sqrt(TMath::Sqrt(nCombBkg_exp))); 
  
  this->fit_MC(dataCut,true);
  this->fit_Kpipipi_misID(misIDCut,true);
  this->fit_HHpipi_misID(misIDCut,true);
  this->fit_normalization_Data(nomalizationCut);  //all shapes fixed

  TRandom3 generatorSig(13031989);
  TRandom3 generatorComb(27051987);
  TRandom3 generatorMisID(27031959);

  // std::cout << "Doing toy study for cut: "<< dataCut << " to "<< fileName << std::endl;
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1780., 1950.,"MeV");

  RooRealVar nCombinatoricBkg("nCombinatoricBkg","number of combinatorial bkg events",nCombBkg_exp,nCombBkg_exp-500,nCombBkg_exp+500);
  RooRealVar nD2hhhhBkg("nD2hhhhBkg","number of misidentified D2hhhh background",nMisID_exp,nMisID_exp-500,nMisID_exp+500);
  RooRealVar nSigExp("nSigExp","number of expected signal events",nSig_exp,nSig_exp-1000,nSig_exp+1000);
  
  if(fixMisID20) {nD2hhhhBkg.setVal(0); nD2hhhhBkg.setConstant(); } //check bias if misID yield is set to zero

  int i=1;

  while(i<10001){//100000
  
    D2hhmumuModel1D* myModel= new D2hhmumuModel1D();

    RooWorkspace m_ws = initializeModel(myModel,D0_M); //get a workspace with all PDFs                                                                                                     
    m_ws.SetName("m_ws");

    m_ws.import(nCombinatoricBkg);
    m_ws.import(nD2hhhhBkg);
    m_ws.import(nSigExp);

    double nSigfluct = generatorSig.Poisson(nSig_exp);
    double nMisIDfluct = generatorMisID.Poisson(nMisID_exp) ;
    double nCombfluct = generatorComb.Poisson(nCombBkg_exp);
  
    std::cout<<"nSigFluc "<< nSigfluct<<::std::endl;

    // m_ws.factory("Gaussian:gaus1(Dst_DTF_D0_M,mu1[1860],sigma1[10])");
    //m_ws.factory("Gaussian:gaus2(Dst_DTF_D0_M,mu2[1840],sigma2[10])");
    //m_ws.factory("Gaussian:gaus3(Dst_DTF_D0_M,mu3[1880],sigma3[10])");
 

    // RooDataSet* ds_signal = m_ws.pdf("gaus1")->generate(D0_M,nSigfluct);
    //RooDataSet* ds_misID = m_ws.pdf("gaus2")->generate(D0_M,nMisIDfluct);
    //RooDataSet* ds_combinatorial = m_ws.pdf("gaus3")->generate(D0_M,nCombfluct);
  
    RooDataSet* ds_signal = m_ws.pdf("Signal")->generate(D0_M,nSigfluct);
    RooDataSet* ds_misID = m_ws.pdf("D2hhhhBkg")->generate(D0_M,nMisIDfluct);
    RooDataSet* ds_combinatorial = m_ws.pdf("CombinatoricBkg")->generate(D0_M,nCombfluct);

    ds_signal->append(*ds_misID);
    ds_signal->append(*ds_combinatorial);

    m_ws.factory("SUM::tot_pdf(nCombinatoricBkg*CombinatoricBkg,nSigExp*Signal,nD2hhhhBkg*D2hhhhBkg)");
    //RooAbsPdf* finalPDF = m_ws.pdf("tot_pdf"); 
  
    //m_ws.factory("SUM::tot_pdf(nSigExp*gaus1)");                                                                              

    RooFitResult *result = m_ws.pdf("tot_pdf")->fitTo(*ds_signal,Save(kTRUE),NumCPU(3));
  
    //pullSigYield->Fill( (m_ws.var("nSigExp")->getVal() - nSigfluct)/m_ws.var("nSigExp")->getError() );
   
    pullSigYield->Fill( (m_ws.var("nSigExp")->getVal() - nSigExp.getVal() )/m_ws.var("nSigExp")->getError() ); 
    pullMisIDYield->Fill( (m_ws.var("nD2hhhhBkg")->getVal() - nD2hhhhBkg.getVal() )/m_ws.var("nD2hhhhBkg")->getError() );
    pullCombYield->Fill( (m_ws.var("nCombinatoricBkg")->getVal() - nCombinatoricBkg.getVal() )/m_ws.var("nCombinatoricBkg")->getError() );

    SigYield->Fill( m_ws.var("nSigExp")->getVal() );
    MisIDYield->Fill( m_ws.var("nD2hhhhBkg")->getVal() );
    CombYield->Fill( m_ws.var("nCombinatoricBkg")->getVal());

    errorSigYield->Fill( m_ws.var("nSigExp")->getError() );
    errorMisIDYield->Fill( m_ws.var("nD2hhhhBkg")->getError() );
    errorCombYield->Fill( m_ws.var("nCombinatoricBkg")->getError() );


    //std::cout << "nSigExp " << m_ws.var("nSigExp")->getVal()  << "pull "<< (m_ws.var("nSigExp")->getVal() - nSigfluct)/m_ws.var("nSigExp")->getError() << "uncertainty "<< m_ws.var("nSigExp")->getError() <<  "generated "<< nSigfluct << std::endl;
   
    delete myModel;
    delete ds_signal;
    delete ds_misID;
    delete ds_combinatorial;
    //delete finalPDF;
    delete result;

    ++i;

  }  

  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  gStyle->SetOptFit(1) ;

  TCanvas* c = new TCanvas("rf801_mcstudy","rf801_mcstudy",900,900) ;
  c->Divide(3,3);
  c->cd(1);
  pullSigYield->Draw();
  pullSigYield->Fit("gaus");
  c->cd(2);
  SigYield->Draw();
  c->cd(3);
  errorSigYield->Draw();
  c->cd(4);
  pullMisIDYield->Draw();
  pullMisIDYield->Fit("gaus");
  c->cd(5);
  MisIDYield->Draw();
  c->cd(6);
  errorMisIDYield->Draw();
  c->cd(7);
  pullCombYield->Draw();
  pullCombYield->Fit("gaus");
  c->cd(8);
  CombYield->Draw();
  c->cd(9);
  errorCombYield->Draw();

  c->SaveAs(targetFile);

  std::cout<<"DONE"<<std::endl;
  
  /*
  RooMCStudy* myRooMCStudy = new RooMCStudy(*m_ws.pdf("tot_pdf"),D0_M,Binned(kFALSE),Extended(),FitOptions(Save(kFALSE),PrintEvalErrors(0)));
  myRooMCStudy->generateAndFit(5000);

  // E x p l o r e   r e s u l t s   o f   s t u d y 
  // ------------------------------------------------

  // Make plots of the distributions of mean, the error on mean and the pull of mean
  RooPlot* frame1 = myRooMCStudy->plotParam(nSigExp,Bins(40)) ;
  RooPlot* frame2 = myRooMCStudy->plotError(nSigExp,Bins(40)) ;
  RooPlot* frame3 = myRooMCStudy->plotPull(nSigExp,Bins(40),FitGauss(kTRUE)) ;

  RooPlot* frame4 = myRooMCStudy->plotParam(nCombinatoricBkg,Bins(40)) ;
  RooPlot* frame5 = myRooMCStudy->plotError(nCombinatoricBkg,Bins(40)) ;
  RooPlot* frame6 = myRooMCStudy->plotPull(nCombinatoricBkg,Bins(40),FitGauss(kTRUE)) ;

  RooPlot* frame7 = myRooMCStudy->plotParam(nD2hhhhBkg,Bins(40)) ;
  RooPlot* frame8 = myRooMCStudy->plotError(nD2hhhhBkg,Bins(40)) ;
  RooPlot* frame9 = myRooMCStudy->plotPull(nD2hhhhBkg,Bins(40),FitGauss(kTRUE)) ;


  // Draw all plots on a canvas
  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;
  TCanvas* c = new TCanvas("rf801_mcstudy","rf801_mcstudy",900,900) ;
  c->Divide(3,3) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;

  c->cd(4) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;
  c->cd(5) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame5->Draw() ;
  c->cd(6) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame6->Draw() ;

  c->cd(7) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame7->Draw() ;
  c->cd(8) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame8->Draw() ;
  c->cd(9) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame9->Draw() ;

  c->SaveAs("MCstudy.eps");

  //  cout << "model written to file " << fileName << endl;                                                                                                                                  */

}


void D2hhmumuFitter1D::setStyle(){

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

double D2hhmumuFitter1D::getMisIDbkgExp(TString cut, TString namePlot){
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1780., 1940.,"MeV");
  TString dmRange = "&& deltaM>144.5 && deltaM<146.5";

  TFile* file;
  file= new TFile(pathToNormData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  tree->SetBranchStatus("nTracks",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);
 
  TTree* cutTree = tree->CopyTree(cut+dmRange);

  RooArgList list =  RooArgList( D0_M );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M));
  cout <<cutTree->GetEntries() <<endl;
  ///Fit                                                                                                                                                                                   
  ///create Model with desired components                                                                                                                                                   
  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  initializeNormalizationModel(myModel,D0_M);
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
  c1->Draw();
  c1->Print("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/selectionOptimisation/fits/"+namePlot+"misID.eps");
  delete c1;

  double nhhhhBkg=myModel->GetWorkspace().var("nD2hhhhBkg")->getValV();

  delete tree;
  delete cutTree;

  file->Close();
  delete file;

  return nhhhhBkg;


}

double D2hhmumuFitter1D::getCombBkg(TString cut,TString namePlot){

  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1780., 1950.,"MeV");
  TString dmRange = "&& deltaM>144.5 && deltaM<146.5";

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
  tree->SetBranchStatus("nTracks",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);

  TTree* cutTree = tree->CopyTree(cut+dmRange);

  // D0_M.setRange(1880,1950);
  D0_M_chebyB.setVal(0);
  D0_M_chebyB.setConstant();
  //D0_M_chebyA.setRange(0,1.5);

  RooArgList list =  RooArgList( D0_M );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M));
  cout <<cutTree->GetEntries() <<endl;

  ///Fit                                                                                                                                                                                    
  ///create Model with desired components                                                                                                                                                   
  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  initializeModel(myModel,D0_M);
  std::string components="CombinatoricBkg";
  RooAbsPdf* finalPDF = myModel->Model(components);

  D0_M.setRange("R1",1880,1950);
  D0_M.setRange("R2",1780,1820);
  D0_M.setRange("signal",1820,1880);

  double sidebands = data->sumEntries("1","R1,R2");
  RooFitResult *result;
  result = finalPDF->fitTo(*data,Range("R1,R2"),Save(kTRUE),Extended(kTRUE),NumCPU(3));                                                                                            

  TCanvas* c1= new TCanvas("");
  c1->Divide(2);
  c1->cd(1);

  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),CutRange("R1,R2"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,Range("lowersideband,uppersideband"),LineColor(kRed),LineWidth(1),Normalization(sidebands,RooAbsReal::NumEvent));                                  
  data->plotOn(frame_m,Name("data"),CutRange("R1,R2"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  c1->Draw();
  c1->Print("/work/mitzel/D2hhmumu/dev/D2KKmumu/img/selectionOptimisation/fits/"+namePlot+"combBkg.eps");
  delete c1;


  //extrapolate background to signal region                                                                                                                                                 
  D0_M.setRange(1780,1950);
  RooRealVar M_chebyA("M_chebyA","D0_M_chebyA",myModel->GetWorkspace().var("D0_M_chebyA")->getValV());
  RooRealVar M_chebyB("M_chebyB","D0_M_chebyB",myModel->GetWorkspace().var("D0_M_chebyB")->getValV());
  RooChebychev CombinatoricBkg("CombinatoricBkg", "Combinatoric Background (M)", D0_M, RooArgList(M_chebyA,M_chebyB));


  RooArgSet variables(D0_M);
  RooAbsReal* fr_sig = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("signal"));
  RooAbsReal* fr_sideband1 = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("R1"));
  RooAbsReal* fr_sideband2 = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("R2"));
  RooRealVar ftot("ftot","ftot",(fr_sideband1->getVal()+fr_sideband2->getVal()) ); 

  //
  //D0_M.setRange("sideband",1880,1950);
  //RooAbsReal* fr_sideband = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("sideband"));
  //RooAbsReal* fr_sideband = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("R1","R2")); 

  std::cout<<myModel->GetWorkspace().var("nCombinatoricBkg")->getValV()<<"  "<< fr_sig->getVal() << "  "<< fr_sideband1->getVal()
	   << "  " << fr_sideband2->getVal() << "   "<< ftot.getVal() << std::endl;
  std::cout<<myModel->GetWorkspace().var("nCombinatoricBkg")->getValV() *  fr_sig->getVal() /ftot.getVal()<<std::endl;       

  // D0_M.setRange(1800,1950);//set back to nominal range                                                                                                                                      

  //double nComb=myModel->GetWorkspace().var("nCombinatoricBkg")->getValV();                                                                                                                
  double nComb = myModel->GetWorkspace().var("nCombinatoricBkg")->getValV() *  fr_sig->getVal() /ftot.getVal();
  
  delete tree;
  delete cutTree;

  file->Close();

  // return myModel->GetWorkspace().var("nCombinatoricBkg")->getValV()+myModel->GetWorkspace().var("nRandomPionBkg")->getValV();                                                           

  delete file;
  return nComb;

}

double D2hhmumuFitter1D::getCombBkgFromDeltaM(TString cut,TString namePlot){

  RooRealVar deltaM("deltaM", "#delta m", 142., 155.,"MeV");
  TString dmRange = "&& deltaM>144.5 && deltaM<146.5";

  D0_M_xi.setRange(144,146);
  D0_M_lambda.setRange(0.1,2);
  D0_M_gamma.setRange(-2,2);
  D0_M_delta.setRange(0.,10);
  D0_M_xi.setVal(145);
  D0_M_lambda.setVal(0.5);
  D0_M_gamma.setVal(-1);
  D0_M_delta.setVal(3);


  TFile* file;
  //file= new TFile(pathToSidebandData,"OPEN");                                                                                                                                             
  file= new TFile(pathToNormData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  tree->SetBranchStatus("nTracks",1);

  TTree* cutTree = tree->CopyTree(cut);

  // D0_M.setRange(1880,1950);
  //D0_M_chebyB.setVal(0);
  //D0_M_chebyB.setConstant();
  //D0_M_chebyA.setRange(0,1.5);

  RooArgList list =  RooArgList( deltaM );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(deltaM));
  cout <<cutTree->GetEntries() <<endl;

  ///Fit                                                                                                                                                                                    
  ///create Model with desired components                                                                                                                                                   
  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  initializeModel(myModel,deltaM);
  std::string components="Signal CombinatoricBkg";
  RooAbsPdf* finalPDF = myModel->Model(components);

  deltaM.setRange("R1",142,155);
  deltaM.setRange("signal",144.5,146.5);

  //double sidebands = data->sumEntries("1","R1,R2");
  RooFitResult *result;
  result = finalPDF->fitTo(*data,Range("R1"),Save(kTRUE),Extended(kTRUE),NumCPU(3));                                                                                            

  TCanvas* c1= new TCanvas("");
  c1->Divide(2);
  c1->cd(1);

  RooPlot* frame_m= deltaM.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,Range("R1"),LineColor(kRed),LineWidth(1));//,Normalization(sidebands,RooAbsReal::NumEvent));                                  
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  c1->Draw();
  //c1->Print("../img/selectionOptimization1D/"+namePlot+"combBkg.eps");
  c1->Print("test.eps");
  delete c1;

  //extrapolate background to signal region                                                                                                                                                 
  RooRealVar M_chebyA("M_chebyA","D0_M_chebyA",myModel->GetWorkspace().var("D0_M_chebyA")->getValV());
  RooRealVar M_chebyB("M_chebyB","D0_M_chebyB",myModel->GetWorkspace().var("D0_M_chebyB")->getValV());
  RooChebychev CombinatoricBkg("CombinatoricBkg", "Combinatoric Background (M)", deltaM, RooArgList(M_chebyA,M_chebyB));

  RooArgSet variables(deltaM);
  RooAbsReal* fr_sig = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("signal"));
  RooAbsReal* fr_sideband = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("R1")); 

  std::cout<<myModel->GetWorkspace().var("nCombinatoricBkg")->getValV()<<"  "<< fr_sig->getVal() << "  "<< fr_sideband->getVal() << std::endl;
  std::cout<<myModel->GetWorkspace().var("nCombinatoricBkg")->getValV() *  fr_sig->getVal() /fr_sideband->getVal()<<std::endl;       

  // D0_M.setRange(1800,1950);//set back to nominal range                                                                                                                                      

  //double nComb=myModel->GetWorkspace().var("nCombinatoricBkg")->getValV();                                                                                                                
  double nComb = myModel->GetWorkspace().var("nCombinatoricBkg")->getValV() *  fr_sig->getVal() /fr_sideband->getVal();
 
  delete tree;
  delete cutTree;
 
  file->Close();

  // return myModel->GetWorkspace().var("nCombinatoricBkg")->getValV()+myModel->GetWorkspace().var("nRandomPionBkg")->getValV();                                                             delete file;
  
  return nComb;

}



RooWorkspace D2hhmumuFitter1D::initializeModel(D2hhmumuModel1D* myModel, RooRealVar D0_M){
  
  //funcition is called by the constructor and initializes a model defined in D2hhmumuModel.h
  // all the variables used in the fit are also initialized by the constructor and can be set indivudially
  //for all channels with the functions set<XX>mumuStartParameters(), XX=Kpi , KK, pipi
  //the model can be built dirctly to fit, it also fills the workspace with all the PDFs.


  std::cout<<"building a model for D2hhmumu dacays..."<<std::endl;
  
  // myModel->Signal(D0_M,D0_M_mean,D0_M_sigma,D0_M_alphaR,D0_M_alphaL, D0_M_nR, D0_M_nL);
   

  myModel->Signal(D0_M,
		  D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
		 ResolutionScale,globalShift
		  );

 
  myModel->CombinatoricBackground(D0_M,
				 D0_M_chebyA,D0_M_chebyB
				 				 );



   myModel->D2hhhhBackground(D0_M,
  			    D0_M_xi_bkg,D0_M_lambda_bkg,D0_M_gamma_bkg,D0_M_delta_bkg
  			    			    );
  //myModel->D2hhhhBackground(D0_M,D0_M_mean_bkg,D0_M_sigma_bkg,D0_M_alphaR_bkg,D0_M_alphaL_bkg, D0_M_nR_bkg, D0_M_nL_bkg);
 
  //myModel->D2hhhhBackground(D0_M,D0_M_mean_bkg,D0_M_sigma_bkg,D0_M_alphaL_bkg, D0_M_nL_bkg);

 
   std::cout<<"return filled workspace with PFDs of all implemented components..done"<<std::endl;
  return myModel->GetWorkspace();
 
}
RooWorkspace D2hhmumuFitter1D::initializeNormalizationModel(D2hhmumuModel1D* myModel, RooRealVar D0_M){
  
  //funcition is called by the constructor and initializes a model defined in D2hhmumuModel.h
  // all the variables used in the fit are also initialized by the constructor and can be set indivudially
  //for all channels with the functions set<XX>mumuStartParameters(), XX=Kpi , KK, pipi
  //the model can be built dirctly to fit, it also fills the workspace with all the PDFs.


  std::cout<<"building a model for D2hhmumu dacays D2Kpimumu normalization channel..."<<std::endl;
  
  myModel->Signal(D0_M,
		  D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
		 ResolutionScale,globalShift
		  );

 
  myModel->CombinatoricBackground(D0_M,
				 D0_M_chebyA,D0_M_chebyB
				 				 );

   myModel->D2hhhhBackground(D0_M,
			     D0_M_xi_bkg_norm,D0_M_lambda_bkg_norm,D0_M_gamma_bkg_norm,D0_M_delta_bkg_norm
  			      ); 
  //myModel->D2hhhhBackground(D0_M,D0_M_mean_bkg,D0_M_sigma_bkg,D0_M_alphaR_bkg,D0_M_alphaL_bkg, D0_M_nR_bkg, D0_M_nL_bkg);
 
  //myModel->D2hhhhBackground(D0_M,D0_M_mean_bkg,D0_M_sigma_bkg,D0_M_alphaL_bkg, D0_M_nL_bkg);

 
   std::cout<<"return filled workspace with PFDs of all implemented components..done"<<std::endl;
  return myModel->GetWorkspace();
 
}


void D2hhmumuFitter1D::fit_PIDinverted_Data(bool fixShape,TString namePlot){

  ///Load file                                                                                                                                                                 
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1820., 1940.,"MeV");

  TFile* file;
  file= new TFile(pathToInvData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);

  D0_M.setRange(1800,1940);

  ///create Model with desired components                                                                                                                                    
  D2hhmumuModel1D* myModel= new D2hhmumuModel1D();
  initializeModel(myModel,D0_M);

  RooDataSet* data = new RooDataSet("data", "data", tree, RooArgSet(D0_M));
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


  std::cout<<"nentries"<<tree->GetEntries()<<std::endl;
  c1->Draw();
  c1->Print("../img/massFit_invertedPID.eps");

  D0_M_xi_bkg.setVal(myModel->GetWorkspace().var("D0_M_xi_bkg")->getValV() );
  D0_M_lambda_bkg.setVal(myModel->GetWorkspace().var("D0_M_lambda_bkg")->getValV() );
  D0_M_gamma_bkg.setVal(myModel->GetWorkspace().var("D0_M_gamma_bkg")->getValV() );
  D0_M_delta_bkg.setVal(myModel->GetWorkspace().var("D0_M_delta_bkg")->getValV() );

  if(fixShape){
    D0_M_xi_bkg.setConstant();                                                                                                                                                           
    D0_M_lambda_bkg.setConstant();                                                                                                  
                                                    
    D0_M_gamma_bkg.setConstant();                                                                                                                                                        
    D0_M_delta_bkg.setConstant();                                                                                                                                                         
  }
  file->Close();
}

void D2hhmumuFitter1D::GausExpModel(int nsig = 100,    // number of signal events                                                                                                           
				    int nbkg = 10000 )  // number of background events                                                                                                                         
{

  RooWorkspace w("w");
  w.factory("Exponential:bkg_pdf(x[1800,1900], a[-0.05,-1,0.])");
  //w.factory("Gaussian:sig_pdf(x, mass[2], sigma[0.3])");                                                                                                                                  

  RooRealVar sigma("sigma", "m(h h #mu #mu)",1865);
  RooRealVar mean("mean", "m(h h #mu #mu)",8);
 
  RooRealVar * x = w.var("x");  // the observable                                                                                                                                           
  RooJohnsonSU Signal("sig", "Signal D^{0} JSU", *x,D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta);                                                                                          
  //RooGaussian Signal ("sig", "Signal D^{0} JSU", *x,mean,sigma);

  std::cout<<"test1"<<std::endl;
  TRandom3* generator = new TRandom3(13031989);
  double offset = generator->Rndm()*1000;
  RooRealVar nSignal("nSignal","#signal events ",5,0,50,"MeV");
  RooRealVar blinding ("blinding","blinding",offset);
  blinding.setConstant();

  blinding.setVal(0);
  double nSiglow = 0;
  double nSighigh=50;

  std::cout<<"test2"<<std::endl;
  RooRealVar nSignal_blind("nSignal_blind","number of signal events blind",nSignal.getVal()+blinding.getVal(),nSiglow+blinding.getVal(),nSighigh+blinding.getVal());
  w.import(nSignal);
  w.import(blinding);
  w.import(nSignal_blind);
  w.import(Signal);
  std::cout<<"test3"<<std::endl;
  RooFormulaVar yield("yield","@0-@1",RooArgList(nSignal_blind,blinding));
  w.import(yield);
  RooUnblindOffset test("test","test","scheisse",100,nSignal) ;
  w.import(test);

  w.factory("expr::yield('nSignal_blind-blinding',nSignal_blind,blinding)");                                                                                                          
  w.factory("SUM:model(yield*sig, nbkg[0,100000]*bkg_pdf)");  // for extended model                                                                                                  
  //w.factory("SUM:model(test*sig, nbkg[0,10000]*bkg_pdf)");  // for extended model                                                                                                  

  // set the desired value of signal and background events                                                                                                                                  
  //w.var("nsig")->setVal(nsig);
    w.var("nbkg")->setVal(1000);
  
  RooAbsPdf * pdf = w.pdf("model");
  std::cout<<"test4"<<std::endl;

  // generate the data                                                                                                                                                                      

  // use fixed random numbers for reproducibility (use 0 for changing every time)                                                                                                           
  RooRandom::randomGenerator()->SetSeed(111);

  // fix number of bins to 50 to plot or to generate data (default is 100 bins)                                                                                                             
  x->setBins(50);

  RooDataSet * data = pdf->generate( *x);  // will generate accordint to total S+B events                                                                                                   
  //RooDataSet * data = pdf->generate( *x, AllBinned());  // will generate accordint to total S+B events                                                                                    
  data->SetName("data");
  w.import(*data);

  data->Print();

  RooPlot * plot = x->frame(Title("Gaussian Signal over Exponential Background"));
  data->plotOn(plot);
  plot->Draw();
  std::cout<<"test5"<<std::endl;

  //RooFitResult * r = pdf->fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
  //r->Print();

  std::cout<<"test6"<<std::endl;

  pdf->plotOn(plot);
  //draw the two separate pdf's                                                                                                                                                             
  //pdf->plotOn(plot, RooFit::Components("bkg_pdf"), RooFit::LineStyle(kDashed) );
  //pdf->plotOn(plot, RooFit::Components("sig"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed) );

  //pdf->paramOn(plot,Layout(0.5,0.9,0.85));
  TCanvas a("a","a");
  plot->Draw();
  a.SaveAs("test.eps");
 
  RooStats::ModelConfig mc("ModelConfig",&w);
  mc.SetPdf(*pdf);
  std::cout<<"test7"<<std::endl;

  mc.SetParametersOfInterest(*w.var("nSignal_blind"));
  mc.SetObservables(*w.var("x"));
  // define set of nuisance parameters                                                                                                                                                      
  std::cout<<"test8"<<std::endl;

  w.defineSet("nuisParams","a,nbkg,D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta");
   //w.defineSet("nuisParams","a,nbkg");
  mc.SetNuisanceParameters(*w.set("nuisParams"));

  // import model in the workspace                                                                                                                                                          
  w.import(mc);

  RooStats::ModelConfig*  sbModel = (RooStats::ModelConfig*) w.obj("ModelConfig");
  //RooRealVar* poi = (RooRealVar*) sbModel->GetParametersOfInterest()->first();
  RooRealVar* poi = dynamic_cast<RooRealVar*> (sbModel->GetParametersOfInterest()->first()); 
  // define the S+B snapshot (this is used for computing the expected significance)
  //poi->setVal(50);
  std::cout<<"test10"<<std::endl;
  sbModel->SetSnapshot(*poi);
  sbModel->SetName("modelconfig");
  // create B only model
  RooStats::ModelConfig * bModel = (RooStats::ModelConfig*) sbModel->Clone();
  bModel->SetName("B_only_model");      
  poi->setVal(0);
  bModel->SetSnapshot( *poi  );

  w.import(*bModel);
  w.import(*sbModel);


 // write the workspace in the file                                                                                                                                                        
  TString fileName = "GausExpModel.root";
  w.writeToFile(fileName,true);
  cout << "model written to file " << fileName << endl;
}




