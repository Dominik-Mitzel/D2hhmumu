#include <iostream>
#include "sWeights.h"
#include <cmath>
#include <iostream>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include <sstream>
#include <RooDataSet.h>
#include "RooGaussModel.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddModel.h"
#include "RooPolynomial.h"
#include "RooTruthModel.h"
#include "RooFitResult.h"
#include "RooDecay.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooDstD0BG.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooCategory.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooHist.h"
#include "RooStats/SPlot.h"
#include "RooTreeDataStore.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooGlobalFunc.h"
#include "Tools.h"


using namespace std;
using namespace RooFit ;
using namespace RooStats;


void calculateSweights(TString pathToFile,TString cuts){

  bool binned=true;
  bool sWeight=true;

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
  file= new TFile(pathToFile,"OPEN");
  TTree* tree = (TTree*) file->Get("BDT_Tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);

  ///Fill all needed variables in RooDataSet                                                                                                                                            
  ///---------------------------------------                                                                                                                                            

  RooRealVar D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1770., 1940.,"MeV");
  RooArgList list =  RooArgList(D0_M);
  RooDataSet* data = new RooDataSet("data", "data", tree, list, "");
  RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(D0_M)));
  RooDataHist* data_binned = data_small->binnedClone();


  ///Define fit model                                                                                                                                                                   
  ///----------------                                                                                                                                                                   

  ///Signal model                                                                                                                                                                       
  ///-----------------------                                                                                                                                                            

  RooRealVar mean1("mu", "mean1", 1865,1863.,1868.);
  RooRealVar sigma1("sigma_{1}", "sigma1", 10.45,5.,25.);

  RooGaussian Gauss1("Gauss1", "Gauss1", D0_M, mean1, sigma1);


  ///Background model                                                                                                                                                                   
  ///-------------------------                                                                                                                                                          

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

  ///Fit                                                                                                                                                                                
  RooFitResult *result;
  if(binned) result = pdf->fitTo(*data_binned,Save(kTRUE),Extended());
  else result = pdf->fitTo(*data,Save(kTRUE),Extended(),NumCPU(3));
  cout << "result is --------------- "<<endl;
  result->Print();

  ///Plot                                                                                                                                                                               
  ///----------                                                                                                                                                                         
  TCanvas* c1= new TCanvas("");
  RooPlot* frame_m= D0_M.frame();
  frame_m->SetTitle("");


  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
  pdf->plotOn(frame_m,Components(Gauss1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
  pdf->plotOn(frame_m,Components(bkg),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
  pdf->paramOn(frame_m,Layout(0.6));
  data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
  frame_m->Draw();
  c1->Print("massFit.eps");

  //sweight part

  mean1.setConstant();
  sigma1.setConstant();
  mean_bkg.setConstant();
  sigma_bkg.setConstant();
  f_bkg.setConstant();
  exp_par.setConstant();

  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",*data,pdf,RooArgList(n_sig,n_bkg));
  gStyle->SetOptStat(0);

  ///Plot the sWeight distributions as a function of mass                                                                                                                       
  TCanvas* SwDs = new TCanvas("Bd sWeight","Bd sWeight distribution");
  TH2 * SwDsHist = (TH2*)data->createHistogram("Dst_DTF_D0_M,n_sig_sw");
  SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
  SwDsHist->SetTitle("");
  //SwDs->Write();                                                                                                                                                              
  SwDsHist->Draw();
  SwDs->Print("../Data_MC_Comparison/sWeight.eps");

  ///Create output file                                                                                                                                                         
  TFile* output = new TFile("../Data_MC_Comparison/sweighted_D2Kpimumu.root","RECREATE");
  TH1* histo = new TH1D("sWeights","sWeights",100,-1,2);
  tree->SetBranchStatus("*",1);
  TTree* new_tree = tree->CopyTree("");
  double w, Dst_DTF_D0_M;
  TBranch* Bra_sw = new_tree->Branch("n_sig_sw",&w);
 
  ///loop over events                                                                                                                                                         
  int numEvents = tree->GetEntries();
  std::cout<<"calculating sWeights for "<<numEvents <<"entries."<<std::endl;


  for(int i=0; i< numEvents; i++){
    if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
    tree->GetEntry(i);
    //if(Dst_D0_M<1820 || Dst_D0_M > 1940) continue;
    w=sData->GetSWeight(i,"n_sig_sw");
    histo->Fill(w);
    Bra_sw->Fill();
  }

  std::cout<<"calculates sWeights"<<std::endl;
  histo->Write();
  new_tree->Write();
  output->Close();

}



std::vector<TH1*> fillHistograms(TString f_input,TString f_output, bool isData=false) {

    TH1::SetDefaultSumw2();

    TChain* fChain;
    if(!isData) fChain = new TChain("DecayTree");
    if(isData) fChain = new TChain("BDT_Tree");

    fChain->AddFile(f_input);

    Int_t           Dst_BKGCAT;   //!                                                                                                                                                           
    Double_t        w;
    if(!isData) w=1;
    Double_t Dst_MAXDOCA,D_MAXDOCA, Dst_MINIP, D_FD_OWNPV, D_FDCHI2_OWNPV, D_DIRA_OWNPV;
    Double_t mu0_IPCHI2_OWNPV , mu1_IPCHI2_OWNPV , h0_IPCHI2_OWNPV , h1_IPCHI2_OWNPV , Slowpi_IPCHI2_OWNPV;
    Double_t Dst_IPCHI2_OWNPV,D_IPCHI2_OWNPV;
    Double_t D_MINIP, D_MINIPCHI2;
    Double_t Dst_PT, D_PT, mu1_PT, mu0_PT, h0_PT,h1_PT, Slowpi_PT;
    Double_t Dst_P, D_P, mu1_P, mu0_P, h0_P,h1_P;
    Double_t mu0_PX, mu0_PY, mu0_PZ, mu1_PX, mu1_PY, mu1_PZ,h0_PX, h0_PY, h0_PZ,h1_PX, h1_PY,h1_PZ;
    Double_t mu0_PIDmu, mu1_PIDmu;
    Double_t D_M, Dst_M, D_DiMuon_Mass;
    Double_t Slowpi_cosh,D_cosh, mu0_cosh, D_ENDVERTEX_CHI2,Slowpi_P;
    Double_t Slowpi_ProbNNghost, mu0_ProbNNghost, mu1_ProbNNghost,h0_ProbNNghost,h1_ProbNNghost;
    Double_t BDT;
    Double_t deltaM, Dst_DTF_D0_M;
    Double_t mu1_ProbNNmu, mu0_ProbNNmu;
    Double_t D_TAU;
    double D_Conemult,Dst_Conemult,D_Coneptasy,Dst_Coneptasy;
    Int_t nTracks_data,nPVs_data,nSPDHits;
    double Slowpi_TRACK_GhostProb,h0_TRACK_GhostProb,h1_TRACK_GhostProb,mu0_TRACK_GhostProb,mu1_TRACK_GhostProb;



    fChain->SetBranchAddress( "Dst_MAXDOCA", &Dst_MAXDOCA );
    fChain->SetBranchAddress( "D_MAXDOCA", &D_MAXDOCA );
    fChain->SetBranchAddress( "Dst_MINIP", & Dst_MINIP);
    fChain->SetBranchAddress( "D_FD_OWNPV", &D_FD_OWNPV );
    fChain->SetBranchAddress( "D_FDCHI2_OWNPV", &D_FDCHI2_OWNPV );
    fChain->SetBranchAddress( "D_DIRA_OWNPV", &D_DIRA_OWNPV );
    fChain->SetBranchAddress( "D_IPCHI2_OWNPV", &D_IPCHI2_OWNPV );
    fChain->SetBranchAddress( "Dst_IPCHI2_OWNPV", &Dst_IPCHI2_OWNPV );
    fChain->SetBranchAddress( "mu0_IPCHI2_OWNPV", &mu0_IPCHI2_OWNPV );
    fChain->SetBranchAddress( "mu1_IPCHI2_OWNPV", &mu1_IPCHI2_OWNPV );
    fChain->SetBranchAddress( "h0_IPCHI2_OWNPV", &h0_IPCHI2_OWNPV );
    fChain->SetBranchAddress( "h1_IPCHI2_OWNPV", &h1_IPCHI2_OWNPV );
    fChain->SetBranchAddress( "Slowpi_IPCHI2_OWNPV", &Slowpi_IPCHI2_OWNPV );
    fChain->SetBranchAddress( "D_MINIP", &D_MINIP );
    fChain->SetBranchAddress( "D_MINIPCHI2", &D_MINIPCHI2 );
    fChain->SetBranchAddress( "Dst_PT", &Dst_PT );
    fChain->SetBranchAddress( "D_PT", &D_PT );
    fChain->SetBranchAddress( "mu1_PT", &mu1_PT );
    fChain->SetBranchAddress( "mu0_PT", &mu0_PT );
    fChain->SetBranchAddress( "h0_PT", &h0_PT );
    fChain->SetBranchAddress( "h1_PT", &h1_PT );
    fChain->SetBranchAddress( "Slowpi_PT", &Slowpi_PT );

    fChain->SetBranchAddress( "Dst_P", &Dst_P );
    fChain->SetBranchAddress( "D_P", &D_P );
    fChain->SetBranchAddress( "mu1_P", &mu1_P );
    fChain->SetBranchAddress( "mu0_P", &mu0_P );
    fChain->SetBranchAddress( "h0_P", &h0_P );
    fChain->SetBranchAddress( "h1_P", &h1_P );
    fChain->SetBranchAddress( "Slowpi_P", &Slowpi_P );


    fChain->SetBranchAddress( "mu0_PX", &mu0_PX );
    fChain->SetBranchAddress( "mu0_PY", &mu0_PY );
    fChain->SetBranchAddress( "mu0_PZ", &mu0_PZ );
    fChain->SetBranchAddress( "mu1_PX", &mu1_PX );
    fChain->SetBranchAddress( "mu1_PY", &mu1_PY );
    fChain->SetBranchAddress( "mu1_PZ", &mu1_PZ );
    fChain->SetBranchAddress( "h0_PX", &h0_PX );
    fChain->SetBranchAddress( "h0_PY", &h0_PY );
    fChain->SetBranchAddress( "h0_PZ", &h0_PZ );
    fChain->SetBranchAddress( "mu0_PIDmu", &mu0_PIDmu );
    fChain->SetBranchAddress( "mu1_PIDmu", &mu1_PIDmu );
    fChain->SetBranchAddress( "mu0_ProbNNmu", &mu0_ProbNNmu );
    fChain->SetBranchAddress( "mu1_ProbNNmu", &mu1_ProbNNmu );
    fChain->SetBranchAddress( "D_M", &D_M );
    fChain->SetBranchAddress( "Dst_M", &Dst_M );
    fChain->SetBranchAddress( "D_DiMuon_Mass", &D_DiMuon_Mass );
    fChain->SetBranchAddress( "Slowpi_P", &Slowpi_P);
    fChain->SetBranchAddress( "Slowpi_cosh", &Slowpi_cosh);
    fChain->SetBranchAddress( "D_ENDVERTEX_CHI2", & D_ENDVERTEX_CHI2);
    fChain->SetBranchAddress( "mu0_ProbNNghost", & mu0_ProbNNghost);
    fChain->SetBranchAddress( "mu1_ProbNNghost", & mu1_ProbNNghost);
    fChain->SetBranchAddress( "h0_ProbNNghost", & h0_ProbNNghost);
    fChain->SetBranchAddress( "h1_ProbNNghost", & h1_ProbNNghost);
    fChain->SetBranchAddress( "Slowpi_ProbNNghost", & Slowpi_ProbNNghost);
    fChain->SetBranchAddress( "deltaM", & deltaM);
    fChain->SetBranchAddress( "Dst_DTF_D0_M", & Dst_DTF_D0_M);
    fChain->SetBranchAddress( "D_TAU",&D_TAU);
    fChain->SetBranchAddress( "D_cosh",&D_cosh);
    fChain->SetBranchAddress( "mu0_cosh",&mu0_cosh);

    fChain->SetBranchAddress("D_Conemult",&D_Conemult);
    fChain->SetBranchAddress("Dst_Conemult",&Dst_Conemult);
    fChain->SetBranchAddress("D_Coneptasy",&D_Coneptasy);
    fChain->SetBranchAddress("Dst_Coneptasy",&Dst_Coneptasy);
    fChain->SetBranchAddress("nSPDHits",&nSPDHits);
    fChain->SetBranchAddress("nTracks_data",&nTracks_data);
    fChain->SetBranchAddress("mu1_TRACK_GhostProb",&mu1_TRACK_GhostProb);
    fChain->SetBranchAddress("mu0_TRACK_GhostProb",&mu0_TRACK_GhostProb);
    fChain->SetBranchAddress("h0_TRACK_GhostProb",&h0_TRACK_GhostProb);
    fChain->SetBranchAddress("h1_TRACK_GhostProb",&h1_TRACK_GhostProb);
    fChain->SetBranchAddress("Slowpi_TRACK_GhostProb",&Slowpi_TRACK_GhostProb);
    fChain->SetBranchAddress( "D_IPCHI2_OWNPV", &D_IPCHI2_OWNPV );
    fChain->SetBranchAddress( "Dst_IPCHI2_OWNPV", &Dst_IPCHI2_OWNPV );



    if(isData)fChain->SetBranchAddress("n_sig_sw", &w);
    TFile* target = new TFile(f_output,"RECREATE");

    std::vector<TH1*> histograms;


    TH1* h_D_MAXDOCA =       new TH1D( "H_D_MAXDOCA","D max DOCA",100,0,0.3);
    TH1* h_D_DIRA_OWNPV =    new TH1D ("H_D_DIRA_OWNPV","D DIRA", 100,0.9998,1);
    TH1* h_D_ENDVERTEX_CHI2 =new TH1D ("H_D_ENDVERTEX_CHI2","D Decay Vertex CHI2",100,0,100);
    TH1* h_D_FDCHI2_OWNPV =  new TH1D( "h_D_FDCHI2_OWNPV","D FD CHI2",100,0,8000);

    TH1* h_Dst_IPCHI2_OWNPV = new TH1D("Dst_IPCHI2_OWNPV","Dst IPCHI2_OWNPV",100,-4,4);
    TH1* h_D_IPCHI2_OWNPV =   new TH1D("D_IPCHI2_OWNPV","D0 IPCHI2_OWNPV",100,-4,4);
    TH1* h_h0_IPCHI2_OWNPV =  new TH1D("K_IPCHI2_OWNPV","K IPCHI2_OWNPV",100,0,4);
    TH1* h_h1_IPCHI2_OWNPV =  new TH1D("pi_IPCHI2_OWNPV","pi IPCHI2_OWNP",100,0,4);
    TH1* h_mu0_IPCHI2_OWNPV = new TH1D("mu0_IPCHI2_OWNPV","Muon1 IPCHI2_OWNPV",100,0,4);
    TH1* h_mu1_IPCHI2_OWNPV = new TH1D("mu1_IPCHI2_OWNPV","Muon2 IPCHI2_OWNPV",100,0,4);
    TH1* h_pis_IPCHI2_OWNPV = new TH1D("pis_IPCHI2_OWNPV","slow pi IPCHI2_OWNPV",100,-4,4);

    TH1* h_Dst_pt =         new TH1D("Dst_PT","Dst pt",100,0,25000);
    TH1* h_D_pt =           new TH1D("D_PT","D0 pt",100,0,20000);
    TH1* h_h0_pt =          new TH1D("K_PT","K pt",100,0,10000);
    TH1* h_h1_pt =          new TH1D("pi_PT","pi pt",100,0,8000);
    TH1* h_mu0_pt =         new TH1D("mu0_PT","Muon1 pt",100,0,8000);
    TH1* h_mu1_pt =         new TH1D("mu1_PT","Muon2 pt",100,0,8000);
    TH1* h_pis_pt =         new TH1D("pis_PT","slow pion pt",100,0,2000);

    TH1* h_Dst_p =          new TH1D("Dst_P","Dst p",100,0,400e3);
    TH1* h_D_p =            new TH1D("D_P","D0 p",100,0,400e3);
    TH1* h_h0_p =           new TH1D("K_P","K p",100,0,150e3);
    TH1* h_h1_p =           new TH1D("pi_P","pi p",100,0,100e3);
    TH1* h_mu0_p =          new TH1D("mu0_P","Muon1 p",100,0,100e3);
    TH1* h_mu1_p =          new TH1D("mu1_P","Muon2 p",100,0,100e3);
    TH1* h_pis_p =          new TH1D("pis_P","slow pi p",100,0,30000);

    TH1* h_dimuonMass =     new TH1D("h_dimuonMass","dimuonMass",100,200,1200);
    TH1* h_D_ptasy =        new TH1D("h_D_ptasy","D pt asymmetry",100,-1,1);
    TH1* h_Dst_ptasy =      new TH1D("h_Dst_ptasy","Dst pt asymmetry",100,-1,1);

    TH1* h_Dst_M =          new TH1D("h_Dst_M","Dst mass",100,1900,2140);
    TH1* h_D_M =            new TH1D("h_D_M","D mass",100,1750,2000);
    TH1* h_dm =             new TH1D("h_dm","delta mass",100,140,155);

    TH1* h_mu0_ProbNN_mu =  new TH1D("mu0_ProbNN_mu","Muon0 ProbNN mu",100,0,1);
    TH1* h_mu1_ProbNN_mu =  new TH1D("mu1_ProbNN_mu","Muon1 ProbNN mu",100,0,1);

    TH1* h_D_cosh =             new TH1D("h_D_cosh","cosh D",100,-1,1);
    TH1* h_mu0_cosh =           new TH1D("h_mo0_cosh","cosh mu0",100,-1,1);
    TH1* h_Slowpi_cosh =        new TH1D("h_Slowpi_cosh","cosh slowpi",100,-1,1);


    /*
    TH1* h_Dst_OWNPV_X=     new TH1D("Dst_OWNPV_X","Dst_OWNPV X",100,0.45,0.8);
    TH1* h_Dst_OWNPV_Y=     new TH1D("h_Dst_OWNPV_Y","Dst_OWNPV Y",100,-.1,.3);
    TD2Kpimumu_BDT.rootH1* h_Dst_OWNPV_Z=     new TH1D("h_Dst_OWNPV_Z","Dst_OWNPV Z",100,-200,200);
    TH1* h_nVeloClusters=   new TH1D("nVeloClusters","nVelo clusters",100,0,6000);
    TH1* h_nVeloTracks=     new TH1D("nVeloTracks","nVelo Tracks",100,0,500);
    TH1* h_nITClusters=     new TH1D("h_nITClusters","nIT clusters",100,0,2500);
    TH1* h_nTTClusters=     new TH1D("h_nTTClusters","nTT clusters",100,0,2500);
    TH1* h_nOTClusters=     new TH1D("h_OTClusters","nOT clusters",100,0,13e3);
    TH1* h_mu1_TRACK_CHI2=  new TH1D("h_mu1_TRACK_CHI2","mu1 track CHI2",100,0,120);
    TH1* h_mu0_TRACK_CHI2=  new TH1D("h_mu0_TRACK_CHI2","m0 track CHI2",100,0,120);
    TH1* h_h1_TRACK_CHI2=   new TH1D("h_h1_TRACK_CHI2","p track CHI2",100,0,120);
    TH1* h_h0_TRACK_CHI2=   new TH1D("h_h0_TRACK_CHI2","K track CHI2",100,0,120);
    TH1* h_Slowpi_TRACK_CHI2= new TH1D("h_Slowpi_TRACK_CHI2","slow pi track CHI2",100,0,120);
    */
    TH1* h_mu1_TRACK_GhostProb=new TH1D("h_mu1_TRACK_GhostProb","mu1 Track ghostProb",100,0,0.25);
    TH1* h_mu0_TRACK_GhostProb=new TH1D("h_mu0_TRACK_GhostProb","mu0 Track ghsotProb",100,0,0.25);
    TH1* h_h1_TRACK_GhostProb=new TH1D("h_h1_TRACK_GhostProb","pi track ghostProb",100,0,0.25);
    TH1* h_h0_TRACK_GhostProb=new TH1D("h_h0_TRACK_GhostProb","K track ghostProb",100,0,0.25);
    TH1* h_Slowpi_TRACK_GhostProb=new TH1D("h_h0_TRACK_GhostProb","slow pion track ghostProb",100,0,0.25);


    //TH1* h_nPVs =           new TH1D("nPVs","nPVs",10,0,10);
    TH1* h_nSPDHits =       new TH1D("nSPDHits","nSPDHits",100,0,1000);
    TH1* h_nTracks =        new TH1D("nTracks","nTracks",100,0,700);


    TH1* h_Dst_cmult =      new TH1D("h_Dst_cmult","Dst cone multiplicity",100,0,120);
    TH1* h_D_cmult =        new TH1D("h_D_cmult","D cone multiplicity",100,0,120);

    /*
    TH1* h_Dst_eta =        new TH1D("Dst_eta","Dst eta",100,1.5,5.5);
    TH1* h_D_eta =          new TH1D("D_eta","D0 eta",100,1.5,5.5);
    TH1* h_h0_eta =         new TH1D("K_eta","K eta",100,1.5,5.5);
    TH1* h_h1_eta =         new TH1D("h1_eta","pi eta",100,1.5,5.5);
    TH1* h_mu0_eta =        new TH1D("mu0_eta","Muon1 eta",100,1.5,5.5);
    TH1* h_mu1_eta =        new TH1D("mu1_eta","Muon2 eta",100,1.5,5.5);
    TH1* h_pis_eta =        new TH1D("pis_eta","slow pi eta",100,1.5,5.5);

    TH1* h_h0_PID_K =       new TH1D("K_PID_K","Kaon PID K",100,-20,120);
    TH1* h_mu0_PID_mu =     new TH1D("mu0_PID_mu","Muon1 PID mu",100,-20,30);
    TH1* h_mu1_PID_mu =     new TH1D("mu1_PID_mu","Muon2 PID mu",100,-20,30);

    TH1* h_h0_ProbNN_K =    new TH1D("K_ProbNN_K","Kaon ProbNN K",100,0,1);
    TH1* h_h1_ProbNN_pi =   new TH1D("pi_ProbNN_pi","Pi ProbNN pi",100,0,1);
    TH1* h_slowpi_ProbNN_pi =new TH1D("slowpi_ProbNN_pi","Slow pi ProbNN pi",100,0,1);
    */

    TH1* h_h0_ProbNN_ghost =    new TH1D("K_ProbNN_ghost","Kaon ProbNN ghost",100,0,1);                                                                                                      
    TH1* h_h1_ProbNN_ghost =    new TH1D("pi_ProbNN_ghost","Pion ProbNN ghost",100,0,1);   
    TH1* h_mu0_ProbNN_ghost =    new TH1D("mu0_ProbNN_ghost","muon0 ProbNN ghost",100,0,1);
    TH1* h_mu1_ProbNN_ghost =    new TH1D("mu1_ProbNN_ghost","muon1 ProbNN ghost",100,0,1);
    TH1* h_Slowpi_ProbNN_ghost =    new TH1D("Slowpi_ProbNN_ghost","Pis ProbNN ghost",100,0,1);

    TLorentzVector pH0,pH1,pMu0, pMu1, pD, pDst,pPis;

    int nEntries=fChain->GetEntries();
    std::cout<<f_input<<" nEntries: "<<nEntries<<std::endl;
    for (int i =0 ;i<nEntries; ++i) {

      fChain->GetEntry(i);
      // if(!isData && Dst_BKGCAT > 20) continue;
      if(mu0_ProbNNmu<0.5 ||mu1_ProbNNmu<0.5|| Dst_DTF_D0_M <1820 || Dst_DTF_D0_M>1940) continue; //PID cut applied in data also to MC!

      /*
      pH1.SetXYZM(h1_PX,h1_PY,h1_PZ,Mass::Pi());
      pH0.SetXYZM(h0_PX,h0_PY,h0_PZ,Mass::K());
      pMu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
      pMu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());

      pD.SetXYZM(D_PX,D_PY,D_PZ,Mass::D0());
      pDst.SetXYZM(Dst_PX,Dst_PY,Dst_PZ,Mass::Ds());
      pPis.SetXYZM(Slowpi_PX,Slowpi_PY,Slowpi_PZ,Mass::Pi());
      */

      h_Dst_M->Fill(Dst_M,w)    ;
      h_D_M ->Fill(Dst_DTF_D0_M,w);
      h_dm ->Fill(deltaM,w);

      /*
      h_Dst_eta ->Fill(pDst.Eta(),w);
      h_D_eta ->Fill(pD.Eta(),w);
      h_h0_eta ->Fill(pH0.Eta(),w);
      h_h1_eta ->Fill(pH1.Eta(),w);
      h_mu0_eta ->Fill(pMu0.Eta(),w);
      h_mu1_eta ->Fill(pMu1.Eta(),w);
      h_pis_eta ->Fill(pPis.Eta(),w);
      */
      h_nSPDHits      ->Fill( nSPDHits ,w);
      h_nTracks       ->Fill( nTracks_data ,w);

      //if(isData) std::cout<<nPV<<"  "<<nTTClusters_data<<"  "<<nITClusters_data<<"  "<<std::endl;                                                                                 
      
      h_dimuonMass    ->Fill( D_DiMuon_Mass ,w);
      
      h_D_ptasy       ->Fill( D_Coneptasy ,w);
      h_Dst_ptasy     ->Fill( Dst_Coneptasy ,w);

      h_Dst_cmult     ->Fill( Dst_Conemult ,w);
      h_D_cmult       ->Fill( D_Conemult ,w);
      
      h_Dst_pt        ->Fill( Dst_PT ,w);
      h_D_pt         ->Fill( D_PT ,w);
      h_h0_pt         ->Fill( h0_PT ,w);
      h_h1_pt         ->Fill( h1_PT ,w);
      h_mu0_pt        ->Fill( mu0_PT ,w);
      h_mu1_pt        ->Fill( mu1_PT ,w);
      h_pis_pt        ->Fill( Slowpi_PT ,w);

      h_Dst_p        ->Fill( Dst_P ,w);
      h_D_p           ->Fill( D_P ,w);
      h_h0_p          ->Fill( h0_P ,w);
      h_h1_p          ->Fill( h1_P ,w);
      h_mu0_p         ->Fill( mu0_P ,w);
      h_mu1_p         ->Fill( mu1_P ,w);
      h_pis_p         ->Fill( Slowpi_P ,w);
      
      /*
      h_h0_PID_K      ->Fill( h0_PIDK ,w);
      //              h_h1_PID_pi     ->Fill( h1_PIDpi ,w);                                                                                                                                           
      h_mu0_PID_mu    ->Fill( mu0_PIDmu ,w);
      h_mu1_PID_mu    ->Fill( mu1_PIDmu ,w);

      h_h0_ProbNN_K  ->Fill(( h0_ProbNNk) ,w);
      h_h1_ProbNN_pi ->Fill(( h1_ProbNNpi) ,w);
      */
      h_mu0_ProbNN_mu->Fill((  mu0_ProbNNmu),w);
      h_mu1_ProbNN_mu ->Fill(( mu1_ProbNNmu ),w);
      //h_slowpi_ProbNN_pi->Fill(( Slowpi_ProbNNpi) ,w);

      h_Dst_IPCHI2_OWNPV  ->Fill( TMath::Log10(Dst_IPCHI2_OWNPV ),w);
      h_D_IPCHI2_OWNPV    ->Fill( TMath::Log10(D_IPCHI2_OWNPV) ,w);
      h_h0_IPCHI2_OWNPV   ->Fill( TMath::Log10(h0_IPCHI2_OWNPV) ,w);
      h_h1_IPCHI2_OWNPV   ->Fill( TMath::Log10(h1_IPCHI2_OWNPV) ,w);
      h_mu0_IPCHI2_OWNPV  ->Fill( TMath::Log10(mu0_IPCHI2_OWNPV) ,w);
      h_mu1_IPCHI2_OWNPV  ->Fill( TMath::Log10(mu1_IPCHI2_OWNPV) ,w);
      h_pis_IPCHI2_OWNPV  ->Fill( TMath::Log10(Slowpi_IPCHI2_OWNPV) ,w);
      

      h_h0_ProbNN_ghost ->Fill(h0_ProbNNghost,w);
      h_h1_ProbNN_ghost ->Fill(h1_ProbNNghost,w);
      h_mu0_ProbNN_ghost->Fill(mu0_ProbNNghost,w);
      h_mu1_ProbNN_ghost->Fill(mu1_ProbNNghost,w);
      h_Slowpi_ProbNN_ghost->Fill(Slowpi_ProbNNghost,w);

      h_D_cosh ->Fill(D_cosh,w);     
      h_mu0_cosh->Fill(mu0_cosh,w);   
      h_Slowpi_cosh ->Fill(Slowpi_cosh,w);


      /*
      h_Dst_OWNPV_X->Fill(Dst_OWNPV_X,w);
      h_Dst_OWNPV_Y->Fill(Dst_OWNPV_Y,w);
      h_Dst_OWNPV_Z->Fill(Dst_OWNPV_Z,w);
      h_nVeloClusters->Fill(nVeloClusters,w);

      if(!isData){
	h_nITClusters->Fill(nITClusters,w);
	h_nTTClusters->Fill(nTTClusters,w);
	h_nOTClusters->Fill(nOTClusters,w);
	h_nPVs          ->Fill( nPVs ,w);
	h_nVeloTracks->Fill(nVeloTracks,w);
      }
      else{
	h_nITClusters->Fill(nITClusters_data,w);
	h_nTTClusters->Fill(nTTClusters_data,w);
	h_nPVs          ->Fill( nPV ,w);
	h_nVeloTracks->Fill(nVeloTracks_data,w);
      }
      
      h_mu1_TRACK_CHI2->Fill(mu1_TRACK_CHI2,w);
      h_mu0_TRACK_CHI2->Fill(mu0_TRACK_CHI2,w);
      h_h1_TRACK_CHI2->Fill(h1_TRACK_CHI2,w);
      h_h0_TRACK_CHI2->Fill(h0_TRACK_CHI2,w);
      h_Slowpi_TRACK_CHI2->Fill(Slowpi_TRACK_CHI2,w);
      */
      h_mu1_TRACK_GhostProb->Fill((mu1_TRACK_GhostProb),w);
      h_mu0_TRACK_GhostProb->Fill((mu0_TRACK_GhostProb),w);
      h_h1_TRACK_GhostProb->Fill((h1_TRACK_GhostProb),w);
      h_h0_TRACK_GhostProb->Fill((h0_TRACK_GhostProb),w);
      h_Slowpi_TRACK_GhostProb->Fill((Slowpi_TRACK_GhostProb),w);
      

      //std::cout<<D_MAXDOCA<<"  "<<D_FDCHI2_OWNPV<<std::endl;                                                                                                                      

      h_D_MAXDOCA ->Fill(D_MAXDOCA,w);
      h_D_DIRA_OWNPV ->Fill(D_DIRA_OWNPV,w);
      h_D_ENDVERTEX_CHI2 ->Fill(D_ENDVERTEX_CHI2,w);
      h_D_FDCHI2_OWNPV ->Fill(D_FDCHI2_OWNPV,w);




    }



/*
    histograms.push_back(h_Dst_OWNPV_X);
    histograms.push_back(h_Dst_OWNPV_Y);
    histograms.push_back(h_Dst_OWNPV_Z);
    histograms.push_back(h_nVeloClusters);
    histograms.push_back(h_nVeloTracks);
    histograms.push_back(h_nITClusters);
    histograms.push_back(h_nTTClusters);
    histograms.push_back(h_nOTClusters);*/
//     histograms.push_back(h_mu1_TRACK_CHI2);
//     histograms.push_back(h_mu0_TRACK_CHI2);
//     histograms.push_back(h_h1_TRACK_CHI2);
//     histograms.push_back(h_h0_TRACK_CHI2);
//     histograms.push_back(h_Slowpi_TRACK_CHI2);
     histograms.push_back(h_mu1_TRACK_GhostProb);
     histograms.push_back(h_mu0_TRACK_GhostProb);
     histograms.push_back(h_h1_TRACK_GhostProb);
     histograms.push_back(h_h0_TRACK_GhostProb);
     histograms.push_back(h_Slowpi_TRACK_GhostProb);
//     histograms.push_back(h_nPVs  );
     histograms.push_back(h_nSPDHits  );
     histograms.push_back(h_nTracks  );
    histograms.push_back(h_dimuonMass  );
    histograms.push_back(h_D_ptasy  );
    histograms.push_back(h_Dst_ptasy  );
    histograms.push_back(h_Dst_cmult  );
    histograms.push_back(h_D_cmult  );
    histograms.push_back(h_Dst_M);
    histograms.push_back(h_D_M);
    histograms.push_back(h_dm );
    histograms.push_back(h_Dst_pt   );
    histograms.push_back(h_D_pt     );
    histograms.push_back(h_h0_pt    );
    histograms.push_back(h_h1_pt    );
    histograms.push_back(h_mu0_pt   );
    histograms.push_back(h_mu1_pt   );
    histograms.push_back(h_pis_pt   );
    histograms.push_back(h_Dst_p        );
    histograms.push_back(h_D_p    );
    histograms.push_back(h_h0_p    );
    histograms.push_back(h_h1_p     );
    histograms.push_back(h_mu0_p   );
    histograms.push_back(h_mu1_p    );
    histograms.push_back(h_pis_p    );
/*    histograms.push_back(h_h0_PID_K  );*/
//     histograms.push_back(h_mu0_PID_mu);
//     histograms.push_back(h_mu1_PID_mu  );
//     histograms.push_back(h_h0_ProbNN_K  );
//     histograms.push_back(h_h1_ProbNN_pi );
    histograms.push_back(h_mu0_ProbNN_mu);
    histograms.push_back(h_mu1_ProbNN_mu );
    //    histograms.push_back(h_slowpi_ProbNN_pi);
    histograms.push_back(h_Dst_IPCHI2_OWNPV  );
    histograms.push_back(h_D_IPCHI2_OWNPV    );
    histograms.push_back(h_h0_IPCHI2_OWNPV   );
    histograms.push_back(h_h1_IPCHI2_OWNPV   );
    histograms.push_back(h_mu0_IPCHI2_OWNPV  );
    histograms.push_back(h_mu1_IPCHI2_OWNPV  );
    histograms.push_back(h_pis_IPCHI2_OWNPV  );
//     histograms.push_back(h_Dst_eta );
//     histograms.push_back(h_D_eta);
//     histograms.push_back(h_h0_eta );
//     histograms.push_back(h_h1_eta );
//     histograms.push_back(h_mu0_eta );
//     histograms.push_back(h_mu1_eta );
//     histograms.push_back(h_pis_eta );

    histograms.push_back(h_D_MAXDOCA);
    histograms.push_back(h_D_DIRA_OWNPV);
    histograms.push_back(h_D_ENDVERTEX_CHI2);
    histograms.push_back(h_D_FDCHI2_OWNPV);

    histograms.push_back(h_h0_ProbNN_ghost);
    histograms.push_back(h_h1_ProbNN_ghost);
    histograms.push_back(h_mu0_ProbNN_ghost);
    histograms.push_back(h_mu1_ProbNN_ghost);
    histograms.push_back(h_Slowpi_ProbNN_ghost);

    histograms.push_back(h_D_cosh);
    histograms.push_back(h_mu0_cosh);
    histograms.push_back(h_Slowpi_cosh);


    return histograms;

  }

void CreateSubPad(TCanvas *canvas, Double_t vfrac= 0.25)
{
  canvas->SetCanvasSize(canvas->GetWw(),(1.+vfrac)*canvas->GetWh());

  Double_t xlow, ylow, xup, yup;
  canvas->GetPad(0)->GetPadPar(xlow,ylow,xup,yup);

  canvas->Divide(1,2);

  TVirtualPad *upPad = canvas->GetPad(1);
  upPad->SetPad(xlow,ylow+vfrac*(yup-ylow),xup,yup);

  TVirtualPad *dwPad = canvas->GetPad(2);
  dwPad->SetPad(xlow,ylow,xup,ylow+vfrac*(yup-ylow));

  canvas->Update();
  return;
}



void data_MC_comparison(){


  TString input_data = "../Data_MC_Comparison/sweighted_D2Kpimumu.root";
  TString output_data = "../Data_MC_Comparison/histos_data_D2Kpimumu.root";
  

  std::vector<TH1*> data_histos = fillHistograms(input_data,output_data,true);
  for (int i=0; i<data_histos.size();++i) {
   //InitHist(sim9_histos_DST[i]);                                                                                                                                                       
    data_histos[i]->SetLineColor(2);
  }

  TString input_MC = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_MCtrainingSample.root";
  TString output_MC = "../Data_MC_Comparison/histos_MC_D2Kpimumu.root";

  std::vector<TH1*> MC_histos = fillHistograms(input_MC,output_MC,false);
  for (int i=0; i<MC_histos.size();++i) {
    MC_histos[i]->SetLineColor(4);
  }


  //drawing                                                                                                                                                                                  
  for (int i=0; i<data_histos.size() ;++i) {
    const char* str=data_histos[i]->GetTitle();

    std::string s =str;
    TCanvas*a=new TCanvas(str,str);
     CreateSubPad(a);
    a->cd(1);
    data_histos[i]->GetYaxis()->SetTitle("normalized");
    data_histos[i]->GetXaxis()->SetTitle(data_histos[i]->GetTitle());

    data_histos[i]->DrawNormalized();
       //Data_histos_DST[i]->DrawNormalized("same");                                                                                                                                         

    MC_histos[i]->DrawNormalized("same");
    
    a->cd(2);
    
    TH1* histo3 = (TH1D*)data_histos[i]->Clone();
    histo3->Scale(1/histo3->Integral());
    MC_histos[i]->Scale(1/MC_histos[i]->Integral());

    histo3->Divide(MC_histos[i]);
    histo3->GetYaxis()->SetTitle("data/Sim");
    histo3->GetYaxis()->SetRangeUser(0.3,2.0);

    histo3->Draw();
    

    a->SaveAs(TString::Format("../Data_MC_Comparison/plots/new_%i.pdf",i));
  }



}

