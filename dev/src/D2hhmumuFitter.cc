#include "D2hhmumuFitter.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooStats/ModelConfig.h"
#include "RooUnblindPrecision.h"
#include "RooUnblindOffset.h"
#include "RooCategory.h"
#include "TBox.h"
#include "TLatex.h"

using namespace std;
using namespace RooFit ;

//this class is is a fitter for Dst->D(hhmumu)pi decays.
//
//
//

D2hhmumuFitter::D2hhmumuFitter():


  EffRatio("EffRatio","Efficieny ratio",1,0,3),
  nNorm("nNorm","number events normalisation channel",2500,100,10000),
  BFsig("BFsig","signal Branching fraction",1e-8,1e-9,1e-6),
  BFnorm("BFnorm","normalizaiton mode Branching fraction",1,0,1),
 
  ResolutionScale("ResolutionScale","ResolutionScale",1,0.9,1.1),
  globalShift("globalShift","globalShift",1,0.9,1.1),

  //signal                                                                                                                                           
  deltaM_xi("deltaM_xi","deltaM_xi",1.45437e+02,144,146),
  deltaM_lambda("deltaM_lambda","deltaM_lambda",5.35039e-01,0.1,2),
  deltaM_gamma("deltaM_gamma","deltaM_gamma",1.00164e-01,-2,2),
  deltaM_delta("deltaM_delta","deltaM_delta",1.07829e+00,0.,10),
  D0_M_xi("D0_M_xi","D0_M_xi",1.86560e+03,1860,1870),
  D0_M_lambda("D0_M_lambda","D0_M_lambda",8.54839e+00,0.1,20),
  D0_M_gamma("D0_M_gamma","D0_M_gamma",-5.42758e-02,-2,40),
  D0_M_delta("D0_M_delta","D0_M_delta",6.09288e-01,0.,10),

  //purely combinatorial background                                                                                                                 
  deltaM_threshold("deltaM_threshold","deltaM_threshold",139.57018),
  deltaM_alpha("deltaM_alpha","deltaM_alpha",3.9761e+00,0,10.),
  D0_M_chebyA("D0_M_chebyA","D0_M_chebyA",-3.5906e-02,-1,1),
  D0_M_chebyB("D0_M_chebyB","D0_M_chebyB",-1.7004e-02,-1,1),
  D0_M_chebyC("D0_M_chebyC","D0_M_chebyC",-1.7882e-02,-1,1),

  //peaking in mD0 and dM                                                                                                                        
  deltaM_xi_bkg("deltaM_xi_bkg","deltaM_xi_bkg",1.45437e+02,144,146),
  deltaM_lambda_bkg("deltaM_lambda_bkg","deltaM_lambda_bkg",5.35039e-01,0.1,2),
  deltaM_gamma_bkg("deltaM_gamma_bkg","deltaM_gamma_bkg",1.00164e-01,-2,2),
  deltaM_delta_bkg("deltaM_delta_bkg","deltaM_delta_bkg",1.07829e+00,0.,10),
  D0_M_xi_bkg("D0_M_xi_bkg","D0_M_xi_bkg",1.8312e+03,1.8012e+03,1.912e+03),
  D0_M_lambda_bkg("D0_M_lambda_bkg","D0_M_lambda_bkg",4.0000e+01,1.0000e+01,10.0000e+01),
  D0_M_gamma_bkg("D0_M_gamma_bkg","D0_M_gamma_bkg",-5.8442e-01,-10.8442e-01,-1.8442e-01),
  D0_M_delta_bkg("D0_M_delta_bkg","D0_M_delta_bkg",7.5826e-01,0,1.5826),

  mean1("mu1", "mean1", 1.8350e3,1835.,1845.),
  sigma1("sigma_{1}", "sigma1", 1.4187e+01,3.,25.)
{

  pathToSignalMC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root";
  //pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root",
  pathToSignalData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root";                                                      
  pathToNormData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root";
  pathToInvData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/inverted_D2Kpimumu_BDT_selected.root";
  pathToSidebandData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT.root";
  pathToKpipipiData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_D2KKmumuBDT.root"; 
  pathToKpipipiHistoData = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/zoz-5000.root";
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
void D2hhmumuFitter::setPathToKpipipiData(TString path){
  pathToKpipipiData = path;
}
void D2hhmumuFitter::setPathToKpipipiHistoData(TString path){
  pathToKpipipiHistoData = path;
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



void D2hhmumuFitter::fit_MC(TString cut="",bool fixShape = true, TString namePlot=""){

  
  //observables
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1820., 1940.,"MeV");
  RooRealVar deltaM("deltaM","#delta m", 143,148,"MeV");
  RooRealVar nSignal("nSignal","#signal events ",100000,0,100000,"MeV");

  ///create Model with desired components
  D2hhmumuModel* myModel= new D2hhmumuModel();
  initializeModel(myModel,D0_M,deltaM); //initializes the model with parameters 
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
  //apply cuts if needed
  TTree* cutTree = tree->CopyTree(cut);
  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M,deltaM));
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
  finalPDF->paramOn(frame_m);
  frame_m->Draw();  
  
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  //finalPDF->paramOn(frame_m2);
  finalPDF->plotOn(frame_m2,LineColor(kRed),LineWidth(1));
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

  c1->Draw();
  c1->Print("../D2KKmumu/img/2D_MC_fit_"+cut+".eps");
  c1->Print(namePlot);

  // c1->Print("../img/massFit_MC.eps");
  //c1->Print(fOut+".eps");                                                                                                                                                     

  std::cout<<"nsig"<<(m_ws.var("nSignal")->getValV())<<std::endl;
  deltaM_xi.setVal(m_ws.var("deltaM_xi")->getValV() );
  deltaM_lambda.setVal(m_ws.var("deltaM_lambda")->getValV() ); 
  deltaM_gamma.setVal(m_ws.var("deltaM_gamma")->getValV() );
  deltaM_delta.setVal(m_ws.var("deltaM_delta")->getValV() );
  D0_M_xi.setVal(m_ws.var("D0_M_xi")->getValV() );
  D0_M_lambda.setVal(m_ws.var("D0_M_lambda")->getValV() );
  D0_M_gamma.setVal(m_ws.var("D0_M_gamma")->getValV() );
  D0_M_delta.setVal(m_ws.var("D0_M_delta")->getValV() );

  // write the workspace in the file
  TString fileName = "MC_model.root";
  //m_ws.writeToFile(fOut+".root",true);
  cout << "model written to file " << fileName << endl;

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
  delete tree;
  delete cutTree;
  file->Close();
  delete file;
  }


void D2hhmumuFitter::fit_normalization_Data(TString cut="",TString namePlot=""){

  //observables
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1780., 1950.,"MeV");
  RooRealVar deltaM("deltaM","#delta m", 140,155,"MeV");
  
  ///create Model with desired components
  D2hhmumuModel* myModel= new D2hhmumuModel();
  RooWorkspace m_ws = initializeModel(myModel,D0_M,deltaM); //get a workspace with all PDFs 

  //set the yields and add the to ws
  RooRealVar nSignal("nSignal","number of signal events",2000,0,5000);                              
  RooRealVar nCombinatoricBkg("nCombinatoricBkg","number of combinatorial bkg events",100,0,10000);                              
  RooRealVar nRandomPionBkg("nRandomPionBkg","number of random pion associated to true DO",200,0,2000);                               
  RooRealVar nD2hhhhBkg("nD2hhhhBkg","number of misidentified D2hhhh background",500,0,2000);                               

  m_ws.import(nSignal);
  m_ws.import(nCombinatoricBkg);
  m_ws.import(nRandomPionBkg);
  m_ws.import(nD2hhhhBkg);

  //create the PDF
  m_ws.factory("SUM::tot(nSignal*Signal,nCombinatoricBkg*CombinatoricBkg,nRandomPionBkg*RandomPionBkg,nD2hhhhBkg*D2hhhhBkg)");
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
  tree->SetBranchStatus("mu0_ProbNNghost",1);
  tree->SetBranchStatus("mu1_ProbNNghost",1);
  tree->SetBranchStatus("h0_ProbNNghost",1);
  tree->SetBranchStatus("h1_ProbNNghost",1);
  tree->SetBranchStatus("h0_ProbNNK",1);
  tree->SetBranchStatus("h1_ProbNNpi",1);
  tree->SetBranchStatus("Slowpi_ProbNNghost",1);
  tree->SetBranchStatus("mu0_MuonNShared",1);
  tree->SetBranchStatus("mu1_MuonNShared",1);

  //apply cuts if needed
  TTree* cutTree = tree->CopyTree(cut);
  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M,deltaM));
  //RooDataSet* data_small = (RooDataSet*) data->reduce(Cut("!( (deltaM>144 && deltaM <147) || (Dst_DTF_D0_M>1840 && Dst_DTF_D0_M<1880) )") );
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
  m_ws.pdf("tot")->plotOn(frame_m,Components(*m_ws.pdf("RandomPionBkg")),LineColor(kViolet),LineStyle(kDashed),LineWidth(1));
  m_ws.pdf("tot")->plotOn(frame_m,Components(*m_ws.pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1)); 
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  finalPDF->paramOn(frame_m);
  frame_m->Draw();
 
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  m_ws.pdf("tot")->plotOn(frame_m2,Components(*m_ws.pdf("Signal")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  m_ws.pdf("tot")->plotOn(frame_m2,Components(*m_ws.pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  m_ws.pdf("tot")->plotOn(frame_m2,Components(*m_ws.pdf("RandomPionBkg")),LineColor(kViolet),LineStyle(kDashed),LineWidth(1));
  m_ws.pdf("tot")->plotOn(frame_m2,Components(*m_ws.pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1)); 
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
  c1->SaveAs("test.root");
  //c1->Print("../img/massFitNormalization.eps");
  c1->Print("../D2KKmumu/img/2D_normalization_fit_"+cut+".eps");
  c1->Print(namePlot);
  
  delete tree;
  delete cutTree;

  delete file;
  
}


void D2hhmumuFitter::fit_Data(TString cut="",TString namePlot=""){

  //observables
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1800., 1940.,"MeV");
  RooRealVar deltaM("deltaM","#delta m", 142,149,"MeV");
  D2hhmumuModel* myModel= new D2hhmumuModel();
  RooWorkspace m_ws = initializeModel(myModel,D0_M,deltaM); //get a workspace with all PDFs                                                               
  m_ws.SetName("m_ws");

  //set the yields and add the to ws                                                                                                                      
  // RooRealVar nSignal("nSignal","number of signal events",2000,0,5000);
  //RooRealVar nCombinatoricBkg("nCombinatoricBkg","number of combinatorial bkg events",100,0,2000);
  // RooRealVar nRandomPionBkg("nRandomPionBkg","number of random pion associated to true DO",200,0,2000);
  //RooRealVar nD2hhhhBkg("nD2hhhhBkg","number of misidentified D2hhhh background",500,0,2000);

  RooRealVar nSignal("nSignal","number of signal events",50,0,200);
  RooRealVar nCombinatoricBkg("nCombinatoricBkg","number of combinatorial bkg events",10,0,300);
  RooRealVar nRandomPionBkg("nRandomPionBkg","number of random pion associated to true DO",10,0,30);
  RooRealVar nD2hhhhBkg("nD2hhhhBkg","number of misidentified D2hhhh background",5,0,20);
  RooUnblindPrecision  nSignal_blind ("nSignal_blind","nSignal blind","scheisse",nSignal.getVal(),nSignal.getVal(),nSignal) ;                                       
  m_ws.import(nSignal);
  m_ws.import(nSignal_blind);
  m_ws.import(nCombinatoricBkg);
  m_ws.import(nRandomPionBkg);
  m_ws.import(nD2hhhhBkg);

  // m_ws.factory("blinding[Unblind=0, Blind=1]");
  //m_ws.cat("blinding")->setLabel("Unblind");

  m_ws.factory("SUM::tot(nSignal_blind*Signal,nCombinatoricBkg*CombinatoricBkg,nRandomPionBkg*RandomPionBkg,nD2hhhhBkg*D2hhhhBkg)");
  RooAbsPdf* finalPDF = m_ws.pdf("tot");
 
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
 
  //apply cuts if needed                                                                                                                                                            
  TTree* cutTree = tree->CopyTree(cut);
  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M,deltaM));
  RooDataSet* data_small = (RooDataSet*) data->reduce(Cut("! (deltaM>145 && deltaM <146) && !(Dst_DTF_D0_M>1850 && Dst_DTF_D0_M<1880) ") );                                      

  m_ws.import(*data);

  //do the fit 
  RooFitResult *result;
  result = finalPDF->fitTo(*data,Save(kTRUE),NumCPU(3));
  m_ws.import(*result);
              
  cout << "result is --------------- "<<endl;
  result->Print();


  //redfintion of PDF in terms of BR 
  BFsig.setVal(1e-7);
  EffRatio.setVal(1); EffRatio.setConstant();
  nNorm.setConstant();
  BFnorm.setVal(4.17e-6);BFnorm.setConstant();
  RooUnblindPrecision  BFblind ("BFblind","Signal BF blind","scheisse",BFsig.getVal(),BFsig.getVal(),BFsig) ;
  m_ws.import(BFblind); // add nuissance parameters to the WS
  m_ws.import(EffRatio);
  m_ws.import(BFnorm);
  m_ws.import(nNorm);

  m_ws.factory("expr::sig_yield('nNorm*(BFblind/BFnorm)*EffRatio',nNorm,BFblind,BFnorm,EffRatio)");  
  m_ws.factory("SUM::tot_pdf(nCombinatoricBkg*CombinatoricBkg,sig_yield*Signal,nRandomPionBkg*RandomPionBkg,nD2hhhhBkg*D2hhhhBkg)");

  /*
  //signal model, simple tests                                                                                                                                                      
  m_ws.defineSet("obs","Dst_DTF_D0_M,deltaM");                                                                                                                                      
  m_ws.defineSet("poi","nSignal");                                                                                                                                                  
  m_ws.defineSet("nuispar","nCombinatoricBkg,nNorm,EffRatio,BFnorm,nRandomPionBkg,nD2hhhhBkg");                                                                          
  m_ws.defineSet("nuispar","D0_M_chebyA,D0_M_xi,deltaM_alpha,deltaM_xi,nCombinatoricBkg");                                                                                
  
  RooStats::ModelConfig* mc = new RooStats::ModelConfig("modelconfig",&m_ws);
  mc->SetPdf(*m_ws.pdf("tot"));
  mc->SetParametersOfInterest(*m_ws.set("poi"));
  mc->SetObservables(*m_ws.set("obs"));
  mc->SetNuisanceParameters(*m_ws.set("nuispar"));
    RooRealVar* poi = dynamic_cast<RooRealVar*> (mc->GetParametersOfInterest()->first());
  // import model in the workspace                                                                                                                                                  
  m_ws.import(*mc);
  mc->SetSnapshot(*poi);

  // create B only model                                                                                                                                                            
  RooStats::ModelConfig * bModel = (RooStats::ModelConfig*) mc->Clone();
  bModel->SetName("B_only_model");
  poi->setVal(0);
  bModel->SetSnapshot( *poi  );
  m_ws.import(*bModel);
  

	      //write the workspace in the file                                                                                                                                                 
  TString fileName = "GausExpModel1.root";
  m_ws.writeToFile(fileName,true);
  cout << "model written to file " << fileName << endl;  */
    //add gaussian constraints to the nuissance parameters
  m_ws.factory("Gaussian::constraintEff(eff0[0,1.5],EffRatio,EffRatio_err[1])");
  m_ws.var("EffRatio")->setVal(EffRatio.getVal());
  m_ws.var("EffRatio_err")->setVal(EffRatio.getVal()*1e-9);
  m_ws.var("eff0")->setVal(EffRatio.getVal());
  m_ws.var("eff0")->setConstant(true);
  std::cout<<"test3"<<std::endl;


  m_ws.factory("Gaussian::constraintNormmodeBF(BFnorm0[74],BFnorm,BFnorm_err[1])");
  m_ws.var("BFnorm")->setVal(BFnorm.getVal());
  m_ws.var("BFnorm_err")->setVal(BFnorm.getVal()*1e-9);
  m_ws.var("BFnorm0")->setVal(BFnorm.getVal());
  m_ws.var("BFnorm0")->setConstant(true);

   m_ws.factory("Gaussian:constraintNormmodeYield(nNorm0[100],nNorm,nNorm_err[1])");
  m_ws.var("nNorm")->setVal(nNorm.getVal());
  m_ws.var("nNorm_err")->setVal(nNorm.getVal()*1e-9);
  m_ws.var("nNorm0")->setVal(nNorm.getVal());
  m_ws.var("nNorm0")->setConstant(true);
  
  //multiply PDF by gaussian constraints
  m_ws.factory("PROD:modelcf(tot_pdf,constraintNormmodeBF,constraintEff,constraintNormmodeYield)");

  RooAbsPdf* finalPDF2 = m_ws.pdf("tot_pdf");
  RooFitResult *result2;
  result2 = finalPDF2->fitTo(*data,Save(kTRUE),NumCPU(3));
  result2->Print();
  TCanvas* c1= new TCanvas("");
  c1->Divide(2,2);
  //c1->cd(1);
 
  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  //m_ws.pdf("tot_pdf")->plotOn(frame_m,Components(*m_ws.pdf("Signal")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  //m_ws.pdf("tot_pdf")->plotOn(frame_m,Components(*m_ws.pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  //m_ws.pdf("tot_pdf")->plotOn(frame_m,Components(*m_ws.pdf("RandomPionBkg")),LineColor(kCyan),LineStyle(kDashed),LineWidth(1));
  //m_ws.pdf("tot_pdf")->plotOn(frame_m,Components(*m_ws.pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1)); 
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  //frame_m->Draw();
  //TPad *newpad1=new TPad("newpad1","a transparent pad",0,0,1,1);
  //newpad1->SetFillStyle(4000);
  //newpad1->Draw();
  //newpad1->cd();
  //TBox *box1 = new TBox(0.33,0.11,0.66,0.88);
  //box1->Draw();  
  //box1->SetFillColor(0);
  //TLatex *t = new TLatex();
  //t->SetTextSize(0.12);
  //t->SetTextColor(kRed);
  //t->DrawLatex(0.35,0.5,"blinded ");
  //c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50));
  //  m_ws.pdf("tot_pdf")->plotOn(frame_m2,Components(*m_ws.pdf("Signal")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  //m_ws.pdf("tot_pdf")->plotOn(frame_m2,Components(*m_ws.pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  //m_ws.pdf("tot_pdf")->plotOn(frame_m2,Components(*m_ws.pdf("RandomPionBkg")),LineColor(kCyan),LineStyle(kDashed),LineWidth(1));
  //m_ws.pdf("tot_pdf")->plotOn(frame_m2,Components(*m_ws.pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1)); 
  finalPDF->plotOn(frame_m2,LineColor(kBlack),LineWidth(1)); 
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),Binning(50),DataError(RooAbsData::SumW2));
  //frame_m2->Draw();
  //TPad *newpad2=new TPad("newpad2","a transparent pad",0,0,1,1);
  //newpad2->SetFillStyle(4000);
  //newpad2->Draw();
  //newpad2->cd();
  //TBox *box2 = new TBox(0.33,0.11,0.66,0.88);
  //box2->Draw();
  //box2->SetFillColor(0);
  //t->DrawLatex(0.35,0.5,"blinded ");
  
  c1->cd(3);
  RooHist* hpull1 = frame_m->pullHist() ;
  RooPlot* frame_m3 = D0_M.frame(Title("Pull Distribution")) ;
  frame_m3->addPlotable(hpull1,"") ; 
  frame_m3->Draw();

  //RooPlot* frame_m3 = D0_M.frame(Title("Pull Distribution")) ;
  //frame_m3->addPlotable(hpull,"HIST") ;
  //c1->cd(3);
  //frame_m3->Draw();

  RooHist* hpull2 = frame_m2->pullHist() ;
  RooPlot* frame_m4 = deltaM.frame(Title("Pull Distribution")) ;
  frame_m4->addPlotable(hpull2,"H") ;
  c1->cd(4);
  frame_m4->Draw();

  c1->cd(1);                                                                                                                                                                      
  RooPlot* frame_m_blind= D0_M.frame();
  frame_m_blind->SetTitle("");
  data_small->plotOn(frame_m_blind,Name("data"),MarkerSize(0.5),Binning(50));
 //m_ws.pdf("tot_pdf")->plotOn(frame_m,Components(*m_ws.pdf("Signal")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));                                                            
 m_ws.pdf("tot_pdf")->plotOn(frame_m_blind,Components(*m_ws.pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));                                          
 m_ws.pdf("tot_pdf")->plotOn(frame_m_blind,Components(*m_ws.pdf("RandomPionBkg")),LineColor(kCyan),LineStyle(kDashed),LineWidth(1));                                              
 m_ws.pdf("tot_pdf")->plotOn(frame_m_blind,Components(*m_ws.pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1));                                        
 //finalPDF->plotOn(frame_m_blind,LineColor(kBlack),LineWidth(1));
 data_small->plotOn(frame_m_blind,Name("data"),MarkerSize(0.5),Binning(50)); 
 frame_m_blind->Draw();                                                                                                                                                              c1->cd(2);
 RooPlot* frame_dm_blind= deltaM.frame();
 frame_dm_blind->SetTitle("");
 data_small->plotOn(frame_dm_blind,Name("data"),MarkerSize(0.5),Binning(50));
 //m_ws.pdf("tot_pdf")->plotOn(frame_m,Components(*m_ws.pdf("Signal")),LineColor(kRed),LineStyle(kDashed),LineWidth(1));                                                            
 m_ws.pdf("tot_pdf")->plotOn(frame_dm_blind,Components(*m_ws.pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));                                          
 m_ws.pdf("tot_pdf")->plotOn(frame_dm_blind,Components(*m_ws.pdf("RandomPionBkg")),LineColor(kCyan),LineStyle(kDashed),LineWidth(1));                                             
 m_ws.pdf("tot_pdf")->plotOn(frame_dm_blind,Components(*m_ws.pdf("D2hhhhBkg")),LineColor(kGreen+3),LineStyle(kDashed),LineWidth(1));                                               
 //finalPDF->plotOn(frame_dm_blind,LineColor(kBlack),LineWidth(1));
 data_small->plotOn(frame_dm_blind,Name("data"),MarkerSize(0.5),Binning(50));
 frame_dm_blind->Draw();        

  std::cout<<"nentries"<<tree->GetEntries()<<std::endl;	
  c1->Draw();
  c1->Print("../img/massFitSignal.eps");
  c1->Print(namePlot);

  m_ws.defineSet("obs","Dst_DTF_D0_M,deltaM");
  m_ws.defineSet("poi","BFsig");
  m_ws.defineSet("nuispar","nCombinatoricBkg,nNorm,EffRatio,BFnorm,nRandomPionBkg,nD2hhhhBkg");
  m_ws.defineSet("gObs","nNorm0,eff0,BFnorm0");

  //create the Model config which is passes to limit calculator
  
  RooStats::ModelConfig*  myModelConfig= new RooStats::ModelConfig("modelconfig");
  myModelConfig->SetWorkspace(m_ws);  
  myModelConfig->SetPdf(*(m_ws.pdf("modelcf")));
  //myModelConfig->SetPdf(*(m_ws.pdf("tot_pdf")));
  myModelConfig->SetObservables(*m_ws.set("obs"));
  myModelConfig->SetParametersOfInterest(*m_ws.set("poi"));
  myModelConfig->SetSnapshot(*m_ws.set("poi"));
  myModelConfig->SetNuisanceParameters(*m_ws.set("nuispar"));
  myModelConfig->SetGlobalObservables(*m_ws.set("gObs"));

  m_ws.import(*myModelConfig);

  //background only model
  RooStats::ModelConfig* bModel = (RooStats::ModelConfig*) myModelConfig->Clone();
  bModel->SetName("modelconfig_withpoi0");
  RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
  if (!var) exit(1);
  double oldval = var->getVal();
  var->setVal(0);
  bModel->SetSnapshot( RooArgSet(*var)  );
  var->setVal(oldval);
  m_ws.import(*bModel);
  

  // write the workspace in the file
  TString fileName = "Data_model.root";
  m_ws.writeToFile(fileName,true);
  cout << "model written to file " << fileName << endl;
  
  
}

void D2hhmumuFitter::fit_Kpipipi_misID(TString cut="",bool fixShape=false,TString namePlot=""){


  RooRealVar D0_M("misID_mD_OS", "m(h h #mu #mu)", 1820., 1940.,"MeV");
  RooRealVar deltaM("deltaM","#delta m", 142,149,"MeV");
                                       
  TFile* file;
  file= new TFile(pathToKpipipiData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("BDT",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  tree->SetBranchStatus("misID_mD_OS",1);
  tree->SetBranchStatus("mu0_PT",1);
  tree->SetBranchStatus("mu1_PT",1);
  tree->SetBranchStatus("nTracks",1);

  D0_M_chebyB.setVal(0);
  D0_M_chebyB.setConstant();
  D0_M.setRange(1780,1920);
  
  TTree* cutTree = tree->CopyTree(cut);
  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M,deltaM));
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));
  RooDataHist* data_binned = data_small->binnedClone();

  ///Fit                                                                                                          

  ///create Model with desired components                                                                                                                                    
  D2hhmumuModel* myModel= new D2hhmumuModel();
  initializeModel(myModel,D0_M,deltaM);
  //  std::string components="D2hhhhBkg CombinatoricBkg";
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
  finalPDF->paramOn(frame_m);
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
  //c1->Print("../img/massFit_Kpipipi.eps");
  c1->Print("../D2KKmumu/img/2D_misID_fit_"+cut+".eps");
  c1->Print(namePlot);

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


  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1780., 1940.,"MeV");
  RooRealVar deltaM("deltaM","#delta m", 142,149,"MeV");
  RooRealVar nSignal("nSignal","#signal events ",100000,0,100000,"MeV");


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

  TTree* cutTree = tree->CopyTree(cut);
   
  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M,deltaM));
  cout <<cutTree->GetEntries() <<endl;

  D2hhmumuModel* myModel= new D2hhmumuModel();
  initializeModel(myModel,D0_M,deltaM); //initializes the model with parameters                                                                                                              
  std::string components="Signal CombinatoricBkg D2hhhhBkg";
  RooAbsPdf* finalPDF =myModel->Model(components); //acutally build the model with default coefficients                                                                                                           
  //RooWorkspace m_ws = myModel->GetWorkspace(); // get workspace with  all PDFs and coefficients                                                                                          
  //RooAbsPdf* finalPDF = m_ws.pdf("D2hhmumuModel");

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
  c1->Print("../img/selectionOptimization/"+namePlot+"misID.eps");
  delete c1;

  double nhhhhBkg=myModel->GetWorkspace().var("nD2hhhhBkg")->getValV();
  std::cout << "nHHHH bkg "<< nhhhhBkg <<std::endl;

  delete tree;
  delete cutTree;
  file->Close();
  delete file;
  
  return nhhhhBkg;
 
}

double D2hhmumuFitter::getCombBkg(TString cut,TString namePlot){

  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1820., 1940.,"MeV");
  RooRealVar deltaM("deltaM","#delta m", 142,149,"MeV");
  RooRealVar nSignal("nSignal","#signal events ",100000,0,100000,"MeV");

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

  TTree* cutTree = tree->CopyTree(cut);

  D0_M.setRange(1780,1950);
  deltaM.setRange(140,155);
  D0_M_chebyB.setVal(0);
  D0_M_chebyB.setConstant();

  RooArgList list =  RooArgList( D0_M,deltaM );
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D0_M,deltaM));
  cout <<cutTree->GetEntries() <<endl;

  ///Fit               

  D2hhmumuModel* myModel= new D2hhmumuModel();
  initializeModel(myModel,D0_M,deltaM); //initializes the model with parameters                                                                                                            
  std::string components="CombinatoricBkg";
  RooAbsPdf* finalPDF=myModel->Model(components); //acutally build the model with default coefficients                                                                                                           
  //RooWorkspace m_ws = myModel->GetWorkspace(); // get workspace with  all PDFs and coefficients                                                                                             
  //RooAbsPdf* finalPDF = m_ws.pdf("D2hhmumuModel");

  D0_M.setRange("signal",1820,1880);
  deltaM.setRange("signal",144.5,146.5);
  D0_M.setRange("R1",1880,1950);
  //deltaM.setRange("R1",144.5,155);
  D0_M.setRange("R2",1780,1820);
  //deltaM.setRange("R2",144.5,155);

  RooFitResult *result;
  result = finalPDF->fitTo(*data,Range("R1,R2"),Save(kTRUE),Extended(kTRUE),NumCPU(3));                                                                                                                   
  TCanvas* c1= new TCanvas("");
  c1->Divide(2);
  c1->cd(1);

  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");
  data->plotOn(frame_m,Name("data"),CutRange("R1,R2"),MarkerSize(0.5),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m,LineColor(kBlack),LineWidth(1),Range("R1,R2"));
  data->plotOn(frame_m,Name("data"),CutRange("R1,R2"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  c1->cd(2);
  RooPlot* frame_m2= deltaM.frame();
  frame_m2->SetTitle("");
  data->plotOn(frame_m2,Name("data"),MarkerSize(0.5),CutRange("R1,R2"),Binning(50));
  myModel->GetWorkspace().pdf("D2hhmumuModel")->plotOn(frame_m2,Components(*myModel->GetWorkspace().pdf("CombinatoricBkg")),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  finalPDF->plotOn(frame_m2,LineColor(kBlack),LineWidth(1),Range("R1,R2"));
  data->plotOn(frame_m2,Name("data"),CutRange("R1,R2"),MarkerSize(0.5),Binning(50));
  frame_m2->Draw();
  c1->Draw();
  c1->Print("../img/selectionOptimization/"+namePlot+"combBkg.eps");
  delete c1;

  //extrapolate background to signal region 
  //D0_M.setRange(1800,1950);
  //deltaM.setRange(143,160);
  
  RooRealVar dM_threshold("dM_threshold","deltaM_threshold",139.57018);
  RooRealVar dM_alpha("dM_alpha","deltaM_alpha",myModel->GetWorkspace().var("deltaM_alpha")->getValV());
  RooRealVar M_chebyA("M_chebyA","D0_M_chebyA",myModel->GetWorkspace().var("D0_M_chebyA")->getValV());
  RooRealVar M_chebyB("M_chebyB","D0_M_chebyB",myModel->GetWorkspace().var("D0_M_chebyB")->getValV());
  RooChebychev CombinatoricBkgM("CombinatoricBkgM", "Combinatoric Background (M)", D0_M, RooArgList(M_chebyA,M_chebyB));
  RooThreshold CombinatoricBkgDm("CombinatoricBkgDm", "Combinatoric Background (#Deltam)", deltaM, dM_threshold, dM_alpha);
  RooProdPdf CombinatoricBkg("CombinatoricBkg", "Combinatoric Background", RooArgList(CombinatoricBkgM, CombinatoricBkgDm));

  RooArgSet variables(D0_M,deltaM);
  RooAbsReal* fr_sig = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("signal"));  
  
  //D0_M.setRange("sideband",1880,1940);
  //deltaM.setRange("sideband",146.5,160);
  
  RooArgSet variables2(D0_M); //only look in mD sidebands, integrate over dM 
  RooAbsReal* fr_sideband1 = CombinatoricBkg.createIntegral(variables,NormSet(variables2),Range("R1"));
  RooAbsReal* fr_sideband2 = CombinatoricBkg.createIntegral(variables,NormSet(variables2),Range("R2"));
  RooRealVar ftot("ftot","ftot",(fr_sideband1->getVal()+fr_sideband2->getVal()) );

  ///RooAbsReal* fr_sideband = CombinatoricBkg.createIntegral(variables,NormSet(variables),Range("sideband"));

  ///std::cout<<myModel->GetWorkspace().var("nCombinatoricBkg")->getValV()<<"  "<< fr_sig->getVal() << "  "<< fr_sideband->getVal()<<std::endl;                             
  ///std::cout<<myModel->GetWorkspace().var("nCombinatoricBkg")->getValV() *  fr_sig->getVal() /fr_sideband->getVal()<<std::endl;                                                                 
  std::cout<<myModel->GetWorkspace().var("nCombinatoricBkg")->getValV()<<"  "<< fr_sig->getVal() << "  "<< fr_sideband1->getVal()
           << "  " << fr_sideband2->getVal() << "   "<< ftot.getVal() << std::endl;
  std::cout<<myModel->GetWorkspace().var("nCombinatoricBkg")->getValV() *  fr_sig->getVal() /ftot.getVal()<<std::endl;  

  ///D0_M.setRange(1800,1940);//set back to nominal range
  ///deltaM.setRange(142,149);
  
  //double nComb=myModel->GetWorkspace().var("nCombinatoricBkg")->getValV();
  ///double nComb = myModel->GetWorkspace().var("nCombinatoricBkg")->getValV() *  fr_sig->getVal() /fr_sideband->getVal();
  double nComb = myModel->GetWorkspace().var("nCombinatoricBkg")->getValV() *  fr_sig->getVal() /ftot.getVal();

  delete tree;
  delete cutTree;
  file->Close();
  delete file;

  // return myModel->GetWorkspace().var("nCombinatoricBkg")->getValV()+myModel->GetWorkspace().var("nRandomPionBkg")->getValV();
  return nComb;
 


}



RooWorkspace D2hhmumuFitter::initializeModel(D2hhmumuModel* myModel, RooRealVar D0_M,RooRealVar deltaM){
  
  //funcition is called by the constructor and initializes a model defined in D2hhmumuModel.h
  // all the variables used in the fit are also initialized by the constructor and can be set indivudially
  //for all channels with the functions set<XX>mumuStartParameters(), XX=Kpi , KK, pipi
  //the model can be built dirctly to fit, it also fills the workspace with all the PDFs.

  std::cout<<"building a model for D2hhmumu dacays..."<<std::endl;
   myModel->Signal(D0_M,deltaM,
                 D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
		   deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta,
		   ResolutionScale,globalShift
                 );

   /*
   myModel->Signal_forLimit(D0_M,deltaM,
		   EffRatio,nNorm,BFsig,BFnorm,
		   D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta,
		   deltaM_xi,deltaM_lambda,deltaM_gamma,deltaM_delta
		   );
   */
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



  myModel->D2hhhhBackground(D0_M,deltaM,
			    D0_M_xi_bkg,D0_M_lambda_bkg,D0_M_gamma_bkg,D0_M_delta_bkg,
			    deltaM_xi_bkg,deltaM_lambda_bkg,deltaM_gamma_bkg,deltaM_delta_bkg
			    );

 
   std::cout<<"return filled workspace with PFDs of all implemented components..done"<<std::endl;
  return myModel->GetWorkspace();
  ;

}


void D2hhmumuFitter::fit_PIDinverted_Data(bool fixShape,TString namePlot){

  ///Load file                                                                                                                                                                 
  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1820., 1940.,"MeV");
  RooRealVar deltaM("deltaM","#delta m", 142,149,"MeV");

  TFile* file;
  file= new TFile(pathToInvData,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("deltaM",1);

  D0_M.setRange(1800,1940);

  ///create Model with desired components                                                                                                                                    
  D2hhmumuModel* myModel= new D2hhmumuModel();
  initializeModel(myModel,D0_M,deltaM);

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
  c1->Print(namePlot);

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
    D0_M_xi_bkg.setConstant();                                                                                                                                                           
    D0_M_lambda_bkg.setConstant();                                                                                                  
                                                    
    D0_M_gamma_bkg.setConstant();                                                                                                                                                        
    D0_M_delta_bkg.setConstant();                                                                                                                                                         
  }
  file->Close();
}

