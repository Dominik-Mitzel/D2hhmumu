#include "TMVA_applications.h"
#include "TMath.h"
#include <math.h>
using namespace TMVA;


void Classification_D2KKmumu(int part) {

  TMVA::Tools::Instance();

  TString signalTree_train;
  TString bkgTree_train;

  TString signalTree_test;
  TString bkgTree_test;

  TString targetFile;

  if(part==1) {

    signalTree_train="DecayTree_even";
    bkgTree_train="sideband/DecayTree_even";

    signalTree_test="DecayTree_odd";
    bkgTree_test="sideband/DecayTree_odd";

    targetFile = "training_D2KKmumu_evenTrained.root";
  }

  if(part==2) {
    signalTree_train="DecayTree_odd";
    bkgTree_train="sideband/DecayTree_odd";

    signalTree_test="DecayTree_even";
    bkgTree_test="sideband/DecayTree_even";

    targetFile = "training_D2KKmumu_oddTrained.root";
  }


  TFile* outputFile = TFile::Open(targetFile,"RECREATE" );

  TChain* signal_train = new TChain(signalTree_train);
  TChain* background_train= new TChain(bkgTree_train);

  TChain* signal_test = new TChain(signalTree_test);
  TChain* background_test= new TChain(bkgTree_test);

  signal_train->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_highStat_MCtrainingSample.root");
  background_train->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_PreselectedSubsample.root");

  signal_test->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_highStat_MCtrainingSample.root");
  background_test->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_PreselectedSubsample.root");

  // Factory                                                                                                                                                                                
  TString factoryOptions = "!V:!Silent:Color:Transformations=I;N:AnalysisType=Classification:DrawProgressBar";
  TMVA::Factory *factory = new TMVA::Factory( targetFile, outputFile , factoryOptions);

  //add variables used in training                                                                                                                                                          

  //pointing                                                                                                                                                                                
  factory->AddVariable("D_MAXDOCA",'F');////                                                                                                                                                
  factory->AddVariable("D_cosh",'F');////                                                                                                                                                   
  factory->AddVariable("mu0_cosh",'F');////                                                                                                                                                 

  factory->AddVariable("log_D_FDCHI2_OWNPV:=log10(D_FDCHI2_OWNPV)",'F');////                                                                                                                
  factory->AddVariable("log_D_DIRA_OWNPV:=log(D_DIRA_OWNPV)",'F');////                                                                                                                   
  //factory->AddVariable("D_DIRA_OWNPV",'F'a);////                                                                                                                    


  factory->AddVariable("D_ENDVERTEX_CHI2",'F');                                                                                                                                            \

  factory->AddVariable("log_Slowpi_IPCHI2_OWNPV:=log10(Slowpi_IPCHI2_OWNPV)",'F');///                                                                                                       

  //factory->AddVariable("D_MINIP",'F');////                                                                                                                                                
  factory->AddVariable("D_MINIPCHI2",'F');////                                                                                                                                              

  //kinematic                                                                                                                                                                               
  factory->AddVariable("Slowpi_P",'F');
  factory->AddVariable("min_mu_PT:=min(mu0_PT,mu1_PT)",'F');////                                                                                                                            
  factory->AddVariable("min_h_PT:=min(h0_PT,h1_PT)",'F');////                                                                                                                               

  factory->AddVariable("Dst_Coneptasy",'F');

  factory->AddVariable("Slowpi_PT",'F');                                                                                                                                                  
  factory->AddVariable("max_mu_PT:=max(mu0_PT,mu1_PT)",'F');////                                                                                                                          



  TCut cutsS=" Slowpi_ProbNNghost<0.5 && h0_ProbNNghost<0.5 && h1_ProbNNghost<0.5 && mu0_ProbNNghost <0.5 && mu1_ProbNNghost<0.5";//"D_DiMuon_Mass>525";                                    
  TCut cutsB="Slowpi_ProbNNghost<0.5 && h0_ProbNNghost<0.5 && h1_ProbNNghost<0.5 && mu0_ProbNNghost <0.5 && mu1_ProbNNghost<0.5";//"D_DiMuon_Mass>525";                                     


  factory->AddSignalTree(signal_train,1.0,TMVA::Types::kTraining);
  factory->AddBackgroundTree(background_train,1.0,TMVA::Types::kTraining);
  factory->AddSignalTree(signal_test,1.0,TMVA::Types::kTesting);
  factory->AddBackgroundTree(background_test,1.0,TMVA::Types::kTesting);

  factory->PrepareTrainingAndTestTree( cutsS, cutsB, "NormMode=None:!V" );



  //factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=250:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" \
		       );
  factory->BookMethod( TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );


  // Train and test MVAs                                                                                                                                                                    
  factory->TrainAllMethods();
  factory->TestAllMethods();

  // Evaluate performances                                                                                                                                                                  
  factory->EvaluateAllMethods();

  outputFile->Close();

  delete factory;

  // Launch the GUI for the root macros                                                                                                                                                     

//if (!gROOT->IsBatch()) TMVAGui( targetFile );                                                                                                                                          

}


void D2KKmumuCrosstraining() {

  Classification_D2KKmumu(1);
  Classification_D2KKmumu(2);


}


void Application_D2KKmumu(TString treeName, TString fileIn, TString fileOut, int part) {


  TString dir    = "weights/";
  TString prefix;
  TString weightfile;
  TString methodName="BDTG method";
  TString dataDir = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/";
  TString targetName =  dataDir + fileOut;

  //select even trained BDT for part 1 and odd for part 2                                                                                                                             
  if(part==1){
    prefix = "training_D2KKmumu_evenTrained.root";
    weightfile = dir + prefix + TString("_")  + TString("BDTG") + TString(".weights.xml");
  }

  if(part==2){
    prefix = "training_D2KKmumu_oddTrained.root";
    weightfile = dir + prefix + TString("_")  + TString("BDTG") + TString(".weights.xml");
  }


  // This loads the library                                                                                                                                                           
  TMVA::Tools::Instance();

  TChain* theTree = new TChain(treeName);
  theTree->AddFile(dataDir+fileIn);


  // --- Create the Reader object. Reader only takes floats!!!                        
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

  // Create a set of variables and declare them to the reader                                                                                                                         
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used                                                                                 

  Float_t Dst_MAXDOCA, D_MAXDOCA, Dst_MINIP, D_FD_OWNPV, D_FDCHI2_OWNPV, D_DIRA_OWNPV;
  Float_t mu0_IPCHI2_OWNPV,mu1_IPCHI2_OWNPV,h0_IPCHI2_OWNPV,h1_IPCHI2_OWNPV,Slowpi_IPCHI2_OWNPV;
  Float_t log_mu1_IPCHI2_OWNPV,log_mu0_IPCHI2_OWNPV,log_h1_IPCHI2_OWNPV,log_h0_IPCHI2_OWNPV,log_Slowpi_IPCHI2_OWNPV;
  Float_t log_D_FD_OWNPV, log_D_FDCHI2_OWNPV, log_D_DIRA_OWNPV;
  Float_t mu0_P, h0_P,h1_P,mu1_P;
  Float_t D_MINIP, D_MINIPCHI2;
  Float_t Dst_PT, D_PT, mu1_PT, mu0_PT, h0_PT,h1_PT, Slowpi_PT,min_mu_PT,min_h_PT;
  Float_t mu0_PX, mu0_PY, mu0_PZ, mu1_PX, mu1_PY, mu1_PZ,h0_PX, h0_PY, h0_PZ,h1_PX, h1_PY, h1_PZ, Slowpi_P;
  Float_t D_M, Dst_M, Slowpi_cosh,D_ENDVERTEX_CHI2,D_P;
  Float_t Dst_Coneptasy, log_D_TAU;
  Float_t log_Slowpi_ProbNNghost, log_mu0_ProbNNghost, log_mu1_ProbNNghost,log_h0_ProbNNghost,log_h1_ProbNNghost;
  Float_t D_cosh,mu0_cosh,max_mu_PT;

  reader->AddVariable("D_MAXDOCA",&D_MAXDOCA);
  reader->AddVariable("D_cosh",&D_cosh);
  reader->AddVariable("mu0_cosh",&mu0_cosh);
  reader->AddVariable("log_D_FDCHI2_OWNPV:=log10(D_FDCHI2_OWNPV)",&log_D_FDCHI2_OWNPV);
  reader->AddVariable("log_D_DIRA_OWNPV:=log(D_DIRA_OWNPV)",&log_D_DIRA_OWNPV);
  //reader->AddVariable("D_DIRA_OWNPV",&D_DIRA_OWNPV);                                                                                                            
  reader->AddVariable("D_ENDVERTEX_CHI2",&D_ENDVERTEX_CHI2);
  reader->AddVariable("log_Slowpi_IPCHI2_OWNPV:=log10(Slowpi_IPCHI2_OWNPV)",&log_Slowpi_IPCHI2_OWNPV);
  reader->AddVariable("D_MINIPCHI2",&D_MINIPCHI2);
  reader->AddVariable("Slowpi_P",&Slowpi_P);
  reader->AddVariable("min_mu_PT:=min(mu0_PT,mu1_PT)",&min_mu_PT);
  reader->AddVariable("min_h_PT:=min(h0_PT,h1_PT)",&min_h_PT);
  reader->AddVariable("Dst_Coneptasy",&Dst_Coneptasy);
  reader->AddVariable("Slowpi_PT",&Slowpi_PT);                                                                                                                                      
  reader->AddVariable("max_mu_PT:=max(mu0_PT,mu1_PT)",&max_mu_PT);                        

  reader->BookMVA( methodName, weightfile );

  // read variables of the data tree and create new tree with BDT variable                                                                                                            

  Double_t user_Dst_MAXDOCA,user_D_MAXDOCA, user_Dst_MINIP, user_D_FD_OWNPV, user_D_FDCHI2_OWNPV, user_D_DIRA_OWNPV;
  Double_t user_mu0_IPCHI2_OWNPV , user_mu1_IPCHI2_OWNPV , user_h0_IPCHI2_OWNPV , user_h1_IPCHI2_OWNPV , user_Slowpi_IPCHI2_OWNPV;
  Double_t user_D_MINIP, user_D_MINIPCHI2;
  Double_t user_Dst_PT, user_D_PT, user_mu1_PT, user_mu0_PT, user_h0_PT,user_h1_PT, user_Slowpi_PT;
  Double_t user_mu0_PX, user_mu0_PY, user_mu0_PZ, user_mu1_PX, user_mu1_PY, user_mu1_PZ,user_h0_PX, user_h0_PY, user_h0_PZ,user_h1_PX, user_h1_PY,user_h1_PZ;
  Double_t user_mu0_PIDmu, user_mu1_PIDmu;
  Double_t user_D_M, user_Dst_M, user_D_DiMuon_Mass;
  Double_t user_Slowpi_cosh,user_D_cosh, user_mu0_cosh, user_D_ENDVERTEX_CHI2,user_Slowpi_P;
  Double_t user_Dst_Coneptasy;
  Double_t user_Slowpi_ProbNNghost, user_mu0_ProbNNghost, user_mu1_ProbNNghost,user_h0_ProbNNghost,user_h1_ProbNNghost;
  Double_t user_BDT;
  Double_t user_deltaM, user_Dst_DTF_D0_M;
  Double_t user_mu1_ProbNNmu, user_mu0_ProbNNmu;
  Double_t user_D_TAU;
  Int_t user_eventNumber;

  theTree->SetBranchAddress( "Dst_MAXDOCA", &user_Dst_MAXDOCA );
  theTree->SetBranchAddress( "D_MAXDOCA", &user_D_MAXDOCA );
  theTree->SetBranchAddress( "Dst_MINIP", & user_Dst_MINIP);
  theTree->SetBranchAddress( "D_FD_OWNPV", &user_D_FD_OWNPV );
  theTree->SetBranchAddress( "D_FDCHI2_OWNPV", &user_D_FDCHI2_OWNPV );
  theTree->SetBranchAddress( "D_DIRA_OWNPV", &user_D_DIRA_OWNPV );
  theTree->SetBranchAddress( "mu0_IPCHI2_OWNPV", &user_mu0_IPCHI2_OWNPV );
  theTree->SetBranchAddress( "mu1_IPCHI2_OWNPV", &user_mu1_IPCHI2_OWNPV );
  theTree->SetBranchAddress( "h0_IPCHI2_OWNPV", &user_h0_IPCHI2_OWNPV );
  theTree->SetBranchAddress( "h1_IPCHI2_OWNPV", &user_h1_IPCHI2_OWNPV );
  theTree->SetBranchAddress( "Slowpi_IPCHI2_OWNPV", &user_Slowpi_IPCHI2_OWNPV );
  theTree->SetBranchAddress( "D_MINIP", &user_D_MINIP );
  theTree->SetBranchAddress( "D_MINIPCHI2", &user_D_MINIPCHI2 );
  theTree->SetBranchAddress( "Dst_PT", &user_Dst_PT );
  theTree->SetBranchAddress( "D_PT", &user_D_PT );
  theTree->SetBranchAddress( "mu1_PT", &user_mu1_PT );
  theTree->SetBranchAddress( "mu0_PT", &user_mu0_PT );
  theTree->SetBranchAddress( "h0_PT", &user_h0_PT );
  theTree->SetBranchAddress( "h1_PT", &user_h1_PT );
  theTree->SetBranchAddress( "Slowpi_PT", &user_Slowpi_PT );
  theTree->SetBranchAddress( "mu0_PX", &user_mu0_PX );
  theTree->SetBranchAddress( "mu0_PY", &user_mu0_PY );
  theTree->SetBranchAddress( "mu0_PZ", &user_mu0_PZ );
  theTree->SetBranchAddress( "mu1_PX", &user_mu1_PX );
  theTree->SetBranchAddress( "mu1_PY", &user_mu1_PY );
  theTree->SetBranchAddress( "mu1_PIDmu", &user_mu1_PIDmu );
  theTree->SetBranchAddress( "mu1_PZ", &user_mu1_PZ );
  theTree->SetBranchAddress( "h0_PX", &user_h0_PX );
  theTree->SetBranchAddress( "h0_PY", &user_h0_PY );
  theTree->SetBranchAddress( "h0_PZ", &user_h0_PZ );
  theTree->SetBranchAddress( "mu0_PIDmu", &user_mu0_PIDmu );
  theTree->SetBranchAddress( "mu0_ProbNNmu", &user_mu0_ProbNNmu );
  theTree->SetBranchAddress( "mu1_ProbNNmu", &user_mu1_ProbNNmu );
  theTree->SetBranchAddress( "D_M", &user_D_M );
  theTree->SetBranchAddress( "Dst_M", &user_Dst_M );
  theTree->SetBranchAddress( "D_DiMuon_Mass", &user_D_DiMuon_Mass );
  theTree->SetBranchAddress( "Slowpi_P", &user_Slowpi_P);
  theTree->SetBranchAddress( "Slowpi_cosh", &user_Slowpi_cosh);
  theTree->SetBranchAddress( "D_ENDVERTEX_CHI2", & user_D_ENDVERTEX_CHI2);
  theTree->SetBranchAddress( "Dst_Coneptasy", & user_Dst_Coneptasy);
  theTree->SetBranchAddress( "mu0_ProbNNghost", & user_mu0_ProbNNghost);
  theTree->SetBranchAddress( "mu1_ProbNNghost", & user_mu1_ProbNNghost);
  theTree->SetBranchAddress( "h0_ProbNNghost", & user_h0_ProbNNghost);
  theTree->SetBranchAddress( "h1_ProbNNghost", & user_h1_ProbNNghost);
  theTree->SetBranchAddress( "Slowpi_ProbNNghost", & user_Slowpi_ProbNNghost);
  theTree->SetBranchAddress( "deltaM", & user_deltaM);
  theTree->SetBranchAddress( "Dst_DTF_D0_M", & user_Dst_DTF_D0_M);
  theTree->SetBranchAddress( "D_TAU",&user_D_TAU);
  theTree->SetBranchAddress( "D_cosh",&user_D_cosh);
  theTree->SetBranchAddress( "mu0_cosh",&user_mu0_cosh);
  theTree->SetBranchAddress( "eventNumber",&user_eventNumber);

  TFile *target  = new TFile(targetName,"RECREATE" );
  TTree* Tree = new TTree("BDT_Tree","BDT_Tree");         // new tree with BDT variable and all relevant other variables                                                              

  Tree->Branch( "Dst_MAXDOCA", &user_Dst_MAXDOCA );
  Tree->Branch( "D_MAXDOCA", &user_D_MAXDOCA );
  Tree->Branch( "Dst_MINIP", &user_Dst_MINIP);
  Tree->Branch( "D_FD_OWNPV,", &user_D_FD_OWNPV );
  Tree->Branch( "D_FDCHI2_OWNPV,", &user_D_FDCHI2_OWNPV );
  Tree->Branch( "D_DIRA_OWNPV", &user_D_DIRA_OWNPV );
  Tree->Branch( "mu0_IPCHI2_OWNPV", &user_mu0_IPCHI2_OWNPV );
  Tree->Branch( "mu1_IPCHI2_OWNPV", &user_mu1_IPCHI2_OWNPV );
  Tree->Branch( "h0_IPCHI2_OWNPV", &user_h0_IPCHI2_OWNPV );
  Tree->Branch( "h1_IPCHI2_OWNPV", &user_h1_IPCHI2_OWNPV );
  Tree->Branch( "Slowpi_IPCHI2_OWNPV", &user_Slowpi_IPCHI2_OWNPV );
  Tree->Branch( "D_MINIP", &user_D_MINIP );
  Tree->Branch( "D_MINIPCHI2", &user_D_MINIPCHI2 );
  Tree->Branch( "Dst_PT", &user_Dst_PT );
  Tree->Branch( "D_PT", &user_D_PT );
  Tree->Branch( "mu1_PT", &user_mu1_PT );
  Tree->Branch( "mu0_PT", &user_mu0_PT );
  Tree->Branch( "h0_PT", &user_h0_PT );
  Tree->Branch( "h1_PT", &user_h1_PT );
  Tree->Branch( "Slowpi_PT", &user_Slowpi_PT );
  Tree->Branch( "Slowpi_P", &user_Slowpi_P );
  Tree->Branch( "mu0_PX", &user_mu0_PX );
  Tree->Branch( "mu0_PY", &user_mu0_PY );
  Tree->Branch( "mu0_PZ", &user_mu0_PZ );
  Tree->Branch( "mu1_PX", &user_mu1_PX );
  Tree->Branch( "mu1_PY", &user_mu1_PY );
  Tree->Branch( "m1_PZ", &user_mu1_PZ );
  Tree->Branch( "h0_PX", &user_h0_PX );
  Tree->Branch( "h0_PY", &user_h0_PY );
  Tree->Branch( "h0_PZ", &user_h0_PZ );
  Tree->Branch( "D_M", &user_D_M );
  Tree->Branch( "Dst_M", &user_Dst_M );
  Tree->Branch( "BDT", &user_BDT );
  Tree->Branch( "mu0_PIDmu", &user_mu0_PIDmu );
  Tree->Branch( "mu1_PIDmu", &user_mu1_PIDmu );
  Tree->Branch( "D_DiMuon_Mass", &user_D_DiMuon_Mass );
  Tree->Branch( "Dst_DTF_D0_M" , &user_Dst_DTF_D0_M);
  Tree->Branch( "deltaM" , &user_deltaM);
  Tree->Branch( "mu0_ProbNNmu", &user_mu0_ProbNNmu );
  Tree->Branch( "mu1_ProbNNmu", &user_mu1_ProbNNmu );
  Tree->Branch( "mu0_cosh", &user_mu0_cosh );
  Tree->Branch( "D_cosh", &user_D_cosh );
  Tree->Branch( "eventNumber", &user_eventNumber );
  Tree->Branch( "mu0_ProbNNghost", & user_mu0_ProbNNghost);
  Tree->Branch( "mu1_ProbNNghost", & user_mu1_ProbNNghost);
  Tree->Branch( "h0_ProbNNghost", & user_h0_ProbNNghost);
  Tree->Branch( "h1_ProbNNghost", & user_h1_ProbNNghost);
  Tree->Branch( "Slowpi_ProbNNghost", & user_Slowpi_ProbNNghost);

  std::cout << "number Events: " << theTree->GetEntries() << std::endl;

  for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

    if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
    //if (ievt>60000) break;                                                                                                                                                    
    theTree->GetEntry(ievt);
    if(user_Slowpi_ProbNNghost>0.5 || user_h1_ProbNNghost>0.5 ||
       user_h0_ProbNNghost>0.5 || user_mu1_ProbNNghost>0.5 ||  user_mu0_ProbNNghost > 0.5) continue;
    //              if(user_mu0_ProbNNmu<0.5 || user_mu1_ProbNNmu<0.5) continue;                                                                                                


    D_MAXDOCA= float(user_D_MAXDOCA);
    Dst_MINIP= float(user_Dst_MINIP);
    D_FDCHI2_OWNPV= float(user_D_FDCHI2_OWNPV);
    D_DIRA_OWNPV= float(user_D_DIRA_OWNPV);
    Slowpi_cosh = float(user_Slowpi_cosh);
    D_cosh = float(user_D_cosh);
    mu0_cosh = float(user_mu0_cosh);
    Dst_Coneptasy = float(user_Dst_Coneptasy);
    //transformations                                                                                                                                                           
    log_Slowpi_IPCHI2_OWNPV = float(TMath::Log10(user_Slowpi_IPCHI2_OWNPV)) ;
    log_D_FD_OWNPV=float( TMath::Log10(user_D_FD_OWNPV) );
    log_D_FDCHI2_OWNPV=float( TMath::Log10(user_D_FDCHI2_OWNPV));
    log_D_DIRA_OWNPV=float( TMath::Log10(user_D_DIRA_OWNPV));
    min_mu_PT = float(TMath::Min(user_mu1_PT,user_mu0_PT));
    min_h_PT = float(TMath::Min(user_h1_PT,user_h0_PT));
    max_mu_PT = float(TMath::Max(user_mu1_PT,user_mu0_PT));
    Slowpi_P = float(user_Slowpi_P);
    Slowpi_PT = float(user_Slowpi_PT);
    D_MINIP= float(user_D_MINIP);
    D_MINIPCHI2= float(user_D_MINIPCHI2);
    D_ENDVERTEX_CHI2 = float(user_D_ENDVERTEX_CHI2);

    user_BDT= reader->EvaluateMVA("BDTG method");
    //if(ievt>50000)break;                                                                                                                                                      

    Tree->Fill();
  }

  target->Write();

  delete reader;
  //delete theTree;                                                                                                                                                                   

  target->Close();
  //target->Delete();                                          

}



void D2KKmumuCrossapplication(){



  //data                                                                                 
  //Application takes name of tree in input file, name of the input file, the target name and part 1 or 2 for even or odd trained events, respectively. 

  /*
  Application_D2KKmumu("data/DecayTree_odd","D2KKmumu_PreselectedSubsample.root","D2KKmumu_BDT_odd.root",1);
  Application_D2KKmumu("data/DecayTree_even","D2KKmumu_PreselectedSubsample.root","D2KKmumu_BDT_even.root",2);

  TChain* myChain1= new TChain("BDT_Tree");
  myChain1->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT_odd.root");
  myChain1->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT_even.root");
  myChain1->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2KKmumu_BDT.root");


  //sideband                                                                                                                                                                          

  Application_D2KKmumu("sideband/DecayTree_odd","D2KKmumu_PreselectedSubsample.root","sideband_D2KKmumu_BDT_odd.root",1);
  Application_D2KKmumu("sideband/DecayTree_even","D2KKmumu_PreselectedSubsample.root","sideband_D2KKmumu_BDT_even.root",2);
  TChain* myChain2= new TChain("BDT_Tree");
  myChain2->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT_odd.root");
  myChain2->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT_even.root");
  myChain2->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT.root");
  */

  //MC                                                                                                                                                                                

  Application_D2KKmumu("DecayTree_odd","D2KKmumu_highStat_MCtrainingSample.root","MC_D2KKmumu_BDT_odd.root",1);
  Application_D2KKmumu("DecayTree_even","D2KKmumu_highStat_MCtrainingSample.root","MC_D2KKmumu_BDT_even.root",2);

  TChain* myChain3= new TChain("BDT_Tree");
  myChain3->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT_odd.root");
  myChain3->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2KKmumu_BDT_even.root");
  myChain3->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_2KKmumu_BDT.root");

  /*
  //normalization mode 
  Application_D2KKmumu("data/DecayTree_odd","D2Kpimumu_PreselectedSubsample.root","D2Kpimumu_D2KKmumuBDT_odd.root",1);
  Application_D2KKmumu("data/DecayTree_even","D2Kpimumu_PreselectedSubsample.root","D2Kpimumu_D2KKmumuBDT_even.root",2);

  TChain* myChain4= new TChain("BDT_Tree");
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT_odd.root");
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT_even.root");
  myChain4->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2KKmumuBDT.root");

  std::cout << "==> TMVAClassificationApplication is done!"  << std::endl;
  */
}


