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

  signal_train->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_MCtrainingSample.root");
  background_train->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_PreselectedSubsample.root");

  signal_test->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_MCtrainingSample.root");
  background_test->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_PreselectedSubsample.root");

  // Factory                                                                                                                                                                              
  TString factoryOptions = "!V:!Silent:Color:Transformations=I;N:AnalysisType=Classification:DrawProgressBar";
  TMVA::Factory *factory = new TMVA::Factory( targetFile, outputFile , factoryOptions);

  //add variables used in training                                                                                                                                                          

  //pointing                                                                                                                                                                                
  factory->AddVariable("D_MAXDOCA",'F');////                                                                                                                                       
  //factory->AddVariable("D_cosh",'F');////                                                                                                                                       
  //factory->AddVariable("mu0_cosh",'F');////                                                                                                                                       
  factory->AddVariable("log_D_FDCHI2_OWNPV:=log10(D_FDCHI2_OWNPV)",'F');////                                                                                                     
  factory->AddVariable("log_D_DIRA_OWNPV:=log(D_DIRA_OWNPV)",'F');////                                                                                                          
  //factory->AddVariable("D_DIRA_OWNPV",'F'a);////                                                                                                              
  factory->AddVariable("D_ENDVERTEX_CHI2",'F');                                                                                                                                   
  factory->AddVariable("log_Slowpi_IPCHI2_OWNPV:=log10(Slowpi_IPCHI2_OWNPV)",'F');///                                                                                              
  //factory->AddVariable("D_MINIP",'F');////                                                                                                                                       
  factory->AddVariable("D_MINIPCHI2",'F');////                                                                                                                                              

  //kinematic                                                                                                                                                                               
  factory->AddVariable("Slowpi_P",'F');
  //factory->AddVariable("min_mu_PT:=min(mu0_PT,mu1_PT)",'F');////remove newest version
  //factory->AddVariable("min_h_PT:=min(h0_PT,h1_PT)",'F');//// remove newst version
  factory->AddVariable("Dst_Coneptasy",'F');
  factory->AddVariable("Slowpi_PT",'F');                                                                                                                                           
  //factory->AddVariable("max_mu_PT:=max(mu0_PT,mu1_PT)",'F');////remove for newest version                                                                                                                           

  TCut cutsS="h0_ProbNNk>0.2 && h1_ProbNNk>0.2 && mu0_MuonNShared==0 && mu1_MuonNShared==0";
  TCut cutsB="h0_ProbNNk>0.2 && h1_ProbNNk>0.2 && mu0_MuonNShared==0 && mu1_MuonNShared==0";


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
void Classification_D2pipimumu(int part) {

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

    targetFile = "training_D2pipimumu_evenTrained.root";
  }

  if(part==2) {
    signalTree_train="DecayTree_odd";
    bkgTree_train="sideband/DecayTree_odd";

    signalTree_test="DecayTree_even";
    bkgTree_test="sideband/DecayTree_even";

    targetFile = "training_D2pipimumu_oddTrained.root";
  }


  TFile* outputFile = TFile::Open(targetFile,"RECREATE" );

  TChain* signal_train = new TChain(signalTree_train);
  TChain* background_train= new TChain(bkgTree_train);

  TChain* signal_test = new TChain(signalTree_test);
  TChain* background_test= new TChain(bkgTree_test);

  signal_train->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_MCtrainingSample.root");
  background_train->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_PreselectedSubsample.root");

  signal_test->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_MCtrainingSample.root");
  background_test->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_PreselectedSubsample.root");

  // Factory                                                                                                                                                                                
  TString factoryOptions = "!V:!Silent:Color:Transformations=I;N:AnalysisType=Classification:DrawProgressBar";
  TMVA::Factory *factory = new TMVA::Factory( targetFile, outputFile , factoryOptions);

  //add variables used in training                                                                                                                                                          

  //pointing                                                                                                                                                                                
  factory->AddVariable("D_MAXDOCA",'F');////                                                                                                                                       
  //factory->AddVariable("D_cosh",'F');////                                                                                                                                       
  //factory->AddVariable("mu0_cosh",'F');////                                                                                                                                       
  factory->AddVariable("log_D_FDCHI2_OWNPV:=log10(D_FDCHI2_OWNPV)",'F');////                                                                                                     
  factory->AddVariable("log_D_DIRA_OWNPV:=log(D_DIRA_OWNPV)",'F');////                                                                                                          
  //factory->AddVariable("D_DIRA_OWNPV",'F'a);////                                                                                                              
  factory->AddVariable("D_ENDVERTEX_CHI2",'F');                                                                                                                                   
  factory->AddVariable("log_Slowpi_IPCHI2_OWNPV:=log10(Slowpi_IPCHI2_OWNPV)",'F');///                                                                                              
  //factory->AddVariable("D_MINIP",'F');////                                                                                                                                       
  factory->AddVariable("D_MINIPCHI2",'F');////                                                                                                                                              

  //kinematic                                                                                                                                                                               
  factory->AddVariable("Slowpi_P",'F');
  //factory->AddVariable("min_mu_PT:=min(mu0_PT,mu1_PT)",'F');////remove newest version
  //factory->AddVariable("min_h_PT:=min(h0_PT,h1_PT)",'F');//// remove newst version
  factory->AddVariable("Dst_Coneptasy",'F');
  factory->AddVariable("Slowpi_PT",'F');                                                                                                                                           
  //factory->AddVariable("max_mu_PT:=max(mu0_PT,mu1_PT)",'F');////remove for newest version                                                                                                                           
  TCut cutsS="h0_ProbNNpi>0.2 && h1_ProbNNpi>0.2 && mu0_MuonNShared==0 && mu1_MuonNShared==0";
  TCut cutsB="h0_ProbNNpi>0.2 && h1_ProbNNpi>0.2 && mu0_MuonNShared==0 && mu1_MuonNShared==0";

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

void D2pipimumuCrosstraining() {

  Classification_D2pipimumu(1);
  Classification_D2pipimumu(2);


}


void Application_D2KKmumu(TString treeName, TString fileIn, TString fileOut, int part,bool isMC = false, bool isNormalizationMode=false, bool skipPID=false) {


  TString dir    = "weights/";
  TString prefix;
  TString weightfile;
  TString methodName="BDTG method";
  TString dataDir = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/";
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
  //reader->AddVariable("D_cosh",&D_cosh);
  //reader->AddVariable("mu0_cosh",&mu0_cosh);
  reader->AddVariable("log_D_FDCHI2_OWNPV:=log10(D_FDCHI2_OWNPV)",&log_D_FDCHI2_OWNPV);
  reader->AddVariable("log_D_DIRA_OWNPV:=log(D_DIRA_OWNPV)",&log_D_DIRA_OWNPV);
  //reader->AddVariable("D_DIRA_OWNPV",&D_DIRA_OWNPV);                                                                                                            
  reader->AddVariable("D_ENDVERTEX_CHI2",&D_ENDVERTEX_CHI2);
  reader->AddVariable("log_Slowpi_IPCHI2_OWNPV:=log10(Slowpi_IPCHI2_OWNPV)",&log_Slowpi_IPCHI2_OWNPV);
  reader->AddVariable("D_MINIPCHI2",&D_MINIPCHI2);
  reader->AddVariable("Slowpi_P",&Slowpi_P);
  //reader->AddVariable("min_mu_PT:=min(mu0_PT,mu1_PT)",&min_mu_PT);///
  //reader->AddVariable("min_h_PT:=min(h0_PT,h1_PT)",&min_h_PT);///
  reader->AddVariable("Dst_Coneptasy",&Dst_Coneptasy);
  reader->AddVariable("Slowpi_PT",&Slowpi_PT);                                                                                                                                      
  //reader->AddVariable("max_mu_PT:=max(mu0_PT,mu1_PT)",&max_mu_PT);///                        

  reader->BookMVA( methodName, weightfile );

  // read variables of the data tree and create new tree with BDT variable                                                                                                            

  Double_t user_Dst_MAXDOCA,user_D_MAXDOCA, user_Dst_MINIP, user_D_FD_OWNPV, user_D_FDCHI2_OWNPV, user_D_DIRA_OWNPV;
  Double_t user_mu0_IPCHI2_OWNPV , user_mu1_IPCHI2_OWNPV , user_h0_IPCHI2_OWNPV , user_h1_IPCHI2_OWNPV , user_Slowpi_IPCHI2_OWNPV;
  Double_t user_D_MINIP, user_D_MINIPCHI2;
  Double_t user_Dst_PT, user_D_PT, user_mu1_PT, user_mu0_PT, user_h0_PT,user_h1_PT, user_Slowpi_PT;
  Double_t user_D_P, user_Dst_P;
  Double_t user_D_PZ, user_Dst_PZ;
  Double_t user_mu0_PX, user_mu0_PY, user_mu0_PZ, user_mu1_PX, user_mu1_PY, user_mu1_PZ,user_h0_PX, user_h0_PY, user_h0_PZ,user_h1_PX, user_h1_PY,user_h1_PZ;
  Double_t user_mu0_P, user_mu1_P, user_h0_P, user_h1_P;
  Double_t user_mu0_PIDmu, user_mu1_PIDmu;
  Double_t user_D_M, user_Dst_M, user_D_DiMuon_Mass;
  Double_t user_Slowpi_cosh,user_D_cosh, user_mu0_cosh, user_D_ENDVERTEX_CHI2,user_Slowpi_P;
  Double_t user_Dst_Coneptasy;
  Double_t user_Slowpi_ProbNNghost, user_mu0_ProbNNghost, user_mu1_ProbNNghost,user_h0_ProbNNghost,user_h1_ProbNNghost;
  Double_t user_BDT;
  Double_t user_deltaM, user_Dst_DTF_D0_M;
  Double_t user_mu1_ProbNNmu, user_mu0_ProbNNmu;
  Double_t user_D_TAU;
  ULong64_t user_eventNumber;
  Int_t user_nTracks_data;
  Int_t user_nTracks_MC;
  Double_t misID_mD_OS=-1000;
  Double_t misID_dm_OS=-1000;
  Double_t mpiKpipi;
  Double_t h0_PIDK, h1_PIDK; 
  double min_muPT, max_muPT;
  double user_mu0_TRACK_Eta;
  double user_mu1_TRACK_Eta;
  double user_h0_TRACK_Eta;
  double user_h1_TRACK_Eta;

  double h0_TRACK_GhostProb;
  double h1_TRACK_GhostProb;
  double mu0_TRACK_GhostProb;
  double mu1_TRACK_GhostProb;
  double Slowpi_TRACK_GhostProb;

  double user_Slowpi_TRACK_Eta;  
  UInt_t          nCandidate;
  ULong64_t       totCandidates;
  UInt_t          runNumber;
  double h0_RichDLLp,h1_RichDLLp;
  double h1_ProbNNp,h0_ProbNNp;

  
  bool mu0_L0MuonDecision_TOS,mu1_L0MuonDecision_TOS,mu0_L0DiMuonDecision_TOS,mu1_L0DiMuonDecision_TOS,h1_L0MuonDecision_TOS,h1_L0DiMuonDecision_TOS,Dst_L0Global_TIS,D_L0Global_TIS,mu0_Hlt1TrackMuonDecision_TOS,mu1_Hlt1TrackMuonDecision_TOS,D_Hlt1TrackAllL0Decision_TOS,D_Hlt1DiMuonHighMassDecision_TOS,D_Hlt1DiMuonLowMassDecision_TOS,mu0_Hlt1SingleMuonNoIPDecision_TOS,mu1_Hlt1SingleMuonNoIPDecision_TOS,mu0_Hlt1SingleMuonHighPTDecision_TOS,mu1_Hlt1SingleMuonHighPTDecision_TOS,h1_0_Hlt1TrackMuonDecision_TOS,D_Hlt2CharmSemilepD02KKMuMuDecision_TOS,D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS,D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS,Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS,Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS,D_Hlt2DiMuonDetachedDecision_TOS,Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS,Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS,Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS,D_Hlt2CharmHadD02HHHH_K3piDecision_TOS,D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS,D_Hlt2CharmHadD02HHHH_4piDecision_TOS,h1_Hlt1TrackMuonDecision_TOS;

 bool mu0_L0MuonDecision_TIS,mu1_L0MuonDecision_TIS,mu0_L0DiMuonDecision_TIS,mu1_L0DiMuonDecision_TIS,h1_L0MuonDecision_TIS,h1_L0DiMuonDecision_TIS,mu0_Hlt1TrackMuonDecision_TIS,mu1_Hlt1TrackMuonDecision_TIS,D_Hlt1TrackAllL0Decision_TIS,D_Hlt1DiMuonHighMassDecision_TIS,D_Hlt1DiMuonLowMassDecision_TIS,mu0_Hlt1SingleMuonNoIPDecision_TIS,mu1_Hlt1SingleMuonNoIPDecision_TIS,mu0_Hlt1SingleMuonHighPTDecision_TIS,mu1_Hlt1SingleMuonHighPTDecision_TIS,h1_0_Hlt1TrackMuonDecision_TIS,D_Hlt2CharmSemilepD02KKMuMuDecision_TIS,D_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS,D_Hlt2CharmSemilepD02KPiMuMuDecision_TIS,Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TIS,Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS,D_Hlt2DiMuonDetachedDecision_TIS,Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TIS,Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS,Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS,D_Hlt2CharmHadD02HHHH_K3piDecision_TIS,D_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS,D_Hlt2CharmHadD02HHHH_4piDecision_TIS,h1_Hlt1TrackMuonDecision_TIS;

 bool mu0_L0MuonDecision_Dec,mu1_L0MuonDecision_Dec,mu0_L0DiMuonDecision_Dec,mu1_L0DiMuonDecision_Dec,h1_L0MuonDecision_Dec,h1_L0DiMuonDecision_Dec,mu0_Hlt1TrackMuonDecision_Dec,mu1_Hlt1TrackMuonDecision_Dec,D_Hlt1TrackAllL0Decision_Dec,D_Hlt1DiMuonHighMassDecision_Dec,D_Hlt1DiMuonLowMassDecision_Dec,mu0_Hlt1SingleMuonNoIPDecision_Dec,mu1_Hlt1SingleMuonNoIPDecision_Dec,mu0_Hlt1SingleMuonHighPTDecision_Dec,mu1_Hlt1SingleMuonHighPTDecision_Dec,h1_0_Hlt1TrackMuonDecision_Dec,D_Hlt2CharmSemilepD02KKMuMuDecision_Dec,D_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec,D_Hlt2CharmSemilepD02KPiMuMuDecision_Dec,Dst_Hlt2CharmHadD02HHXDst_hhXDecision_Dec,Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec,D_Hlt2DiMuonDetachedDecision_Dec,Dst_Hlt2CharmHadD02HHHHDst_4piDecision_Dec,Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec,Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec,D_Hlt2CharmHadD02HHHH_K3piDecision_Dec,D_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec,D_Hlt2CharmHadD02HHHH_4piDecision_Dec,h1_Hlt1TrackMuonDecision_Dec;

bool D_L0DiMuonDecision_TOS,D_L0MuonDecision_TOS,D_L0HadronDecision_TOS,D_L0PhotonDecision_TOS,D_L0ElectronDecision_TOS;
bool Dst_L0DiMuonDecision_TOS,Dst_L0MuonDecision_TOS,Dst_L0HadronDecision_TOS,Dst_L0PhotonDecision_TOS,Dst_L0ElectronDecision_TOS;

bool D_L0DiMuonDecision_TIS,D_L0MuonDecision_TIS,D_L0HadronDecision_TIS,D_L0PhotonDecision_TIS,D_L0ElectronDecision_TIS;
bool Dst_L0DiMuonDecision_TIS,Dst_L0MuonDecision_TIS,Dst_L0HadronDecision_TIS,Dst_L0PhotonDecision_TIS,Dst_L0ElectronDecision_TIS;

bool D_L0DiMuonDecision_Dec,D_L0MuonDecision_Dec,D_L0HadronDecision_Dec,D_L0PhotonDecision_Dec,D_L0ElectronDecision_Dec;
bool Dst_L0DiMuonDecision_Dec,Dst_L0MuonDecision_Dec,Dst_L0HadronDecision_Dec,Dst_L0PhotonDecision_Dec,Dst_L0ElectronDecision_Dec;


  //MC Variables                                                                                                                                                            
  Int_t           Dst_BKGCAT;   //!                                                                                                                                      
  Int_t           Dst_TRUEID;   //!                                                                                                                                      
  Int_t           Dst_MC_MOTHER_ID;
  Int_t           Dst_MC_MOTHER_KEY;
  Int_t           Dst_MC_GD_MOTHER_ID;
  Int_t           Dst_MC_GD_MOTHER_KEY;
  Int_t           Dst_MC_GD_GD_MOTHER_ID;
  Int_t           Dst_MC_GD_GD_MOTHER_KEY;
  Double_t        Dst_TRUEP_E;
  Double_t        Dst_TRUEP_X;
  Double_t        Dst_TRUEP_Y;
  Double_t        Dst_TRUEP_Z;
  Double_t        Dst_TRUEPT;
  Double_t        Dst_TRUEORIGINVERTEX_X;
  Double_t        Dst_TRUEORIGINVERTEX_Y;
  Double_t        Dst_TRUEORIGINVERTEX_Z;
  Double_t        Dst_TRUEENDVERTEX_X;
  Double_t        Dst_TRUEENDVERTEX_Y;
  Double_t        Dst_TRUEENDVERTEX_Z;
  Bool_t          Dst_TRUEISSTABLE;
  Double_t        Dst_TRUETAU;
  Int_t           D_BKGCAT;
  Int_t           D_TRUEID;
  Int_t           D_MC_MOTHER_ID;
  Int_t           D_MC_MOTHER_KEY;
  Double_t        D_TRUEP_E;
  Double_t        D_TRUEP_X;
  Double_t        D_TRUEP_Y;
  Double_t        D_TRUEP_Z;
  Double_t        D_TRUEPT;
  Double_t        D_TRUEORIGINVERTEX_X;
  Double_t        D_TRUEORIGINVERTEX_Y;
  Double_t        D_TRUEORIGINVERTEX_Z;
  Double_t        D_TRUEENDVERTEX_X;
  Double_t        D_TRUEENDVERTEX_Y;
  Double_t        D_TRUEENDVERTEX_Z;
  Bool_t          D_TRUEISSTABLE;
  Double_t        D_TRUETAU;
  Int_t           h0_TRUEID;
  Int_t           h0_MC_MOTHER_ID;
  Double_t        h0_TRUEP_E;
  Double_t        h0_TRUEP_X;
  Double_t        h0_TRUEP_Y;
  Double_t        h0_TRUEP_Z;
  Double_t        h0_TRUEPT;
  Double_t        h0_TRUEORIGINVERTEX_X;
  Double_t        h0_TRUEORIGINVERTEX_Y;
  Double_t        h0_TRUEORIGINVERTEX_Z;
  Bool_t          h0_TRUEISSTABLE;
  Double_t        h0_TRUETAU;
  Int_t           h1_TRUEID;
  Int_t           h1_MC_MOTHER_ID;
  Double_t        h1_TRUEP_E;
  Double_t        h1_TRUEP_X;
  Double_t        h1_TRUEP_Y;
  Double_t        h1_TRUEP_Z;
  Double_t        h1_TRUEPT;
  Double_t        h1_TRUEORIGINVERTEX_X;
  Double_t        h1_TRUEORIGINVERTEX_Y;
  Double_t        h1_TRUEORIGINVERTEX_Z;
  Bool_t          h1_TRUEISSTABLE;
  Double_t        h1_TRUETAU;
  Int_t           mu0_TRUEID;
  Int_t           mu0_MC_MOTHER_ID;
  Double_t        mu0_TRUEP_E;
  Double_t        mu0_TRUEP_X;
  Double_t        mu0_TRUEP_Y;
  Double_t        mu0_TRUEP_Z;
  Double_t        mu0_TRUEPT;
  Double_t        mu0_TRUEORIGINVERTEX_X;
  Double_t        mu0_TRUEORIGINVERTEX_Y;
  Double_t        mu0_TRUEORIGINVERTEX_Z;
  Double_t        mu0_TRUEENDVERTEX_X;
  Double_t        mu0_TRUEENDVERTEX_Y;
  Double_t        mu0_TRUEENDVERTEX_Z;
  Bool_t          mu0_TRUEISSTABLE;
  Double_t        mu0_TRUETAU;
  Int_t           mu1_TRUEID;
  Int_t           mu1_MC_MOTHER_ID;
  Int_t           mu1_MC_MOTHER_KEY;
  Double_t        mu1_TRUEP_E;
  Double_t        mu1_TRUEP_X;
  Double_t        mu1_TRUEP_Y;
  Double_t        mu1_TRUEP_Z;
  Double_t        mu1_TRUEPT;
  Double_t        mu1_TRUEORIGINVERTEX_X;
  Double_t        mu1_TRUEORIGINVERTEX_Y;
  Double_t        mu1_TRUEORIGINVERTEX_Z;
  Bool_t          mu1_TRUEISSTABLE;
  Double_t        mu1_TRUETAU;
  Int_t           Slowpi_TRUEID;
  Int_t           Slowpi_MC_MOTHER_ID;
  Double_t        Slowpi_TRUEP_E;
  Double_t        Slowpi_TRUEP_X;
  Double_t        Slowpi_TRUEP_Y;
  Double_t        Slowpi_TRUEP_Z;
  Double_t        Slowpi_TRUEPT;
  Double_t        Slowpi_TRUEORIGINVERTEX_X;
  Double_t        Slowpi_TRUEORIGINVERTEX_Y;
  Double_t        Slowpi_TRUEORIGINVERTEX_Z;
  Bool_t          Slowpi_TRUEISSTABLE;
  Double_t        Slowpi_TRUETAU;
  Double_t        user_mHH;
  Double_t user_h0_ProbNNk,user_h1_ProbNNk,user_h0_ProbNNpi,user_h1_ProbNNpi;
  Double_t D_TRUE_DiMuon_Mass;
  Int_t           mu1_MuonNShared;
  Int_t           mu0_MuonNShared;
  Bool_t          mu0_isMuon;
  Bool_t          mu1_isMuon;
  double mumu_P,mumu_PT;
  Short_t         Polarity;

  bool Dst_Hlt1Phys_TIS,D_Hlt1GlobalDecision_TIS,D_Hlt1L0AnyDecision_TIS;
  double mKPiMuMu,mKMuPiMu,mKMuMuPi;
  double m34,m24,m23;

  double DTF_h0_PX,DTF_h0_PY,DTF_h0_PZ;
  double DTF_h1_PX,DTF_h1_PY,DTF_h1_PZ;
  double DTF_mu0_PX,DTF_mu0_PY,DTF_mu0_PZ;
  double DTF_mu1_PX,DTF_mu1_PY,DTF_mu1_PZ;

  theTree->SetBranchAddress( "Dst_DTF_h0_PX", &DTF_h0_PX);
  theTree->SetBranchAddress( "Dst_DTF_h0_PY", &DTF_h0_PY );
  theTree->SetBranchAddress( "Dst_DTF_h0_PZ", &DTF_h0_PZ );
  theTree->SetBranchAddress( "Dst_DTF_h1_PX", &DTF_h1_PX);
  theTree->SetBranchAddress( "Dst_DTF_h1_PY", &DTF_h1_PY );
  theTree->SetBranchAddress( "Dst_DTF_h1_PZ", &DTF_h1_PZ );
  theTree->SetBranchAddress( "Dst_DTF_mu0_PX", &DTF_mu0_PX);
  theTree->SetBranchAddress( "Dst_DTF_mu0_PY", &DTF_mu0_PY );
  theTree->SetBranchAddress( "Dst_DTF_mu0_PZ", &DTF_mu0_PZ );
  theTree->SetBranchAddress( "Dst_DTF_mu1_PX", &DTF_mu1_PX);
  theTree->SetBranchAddress( "Dst_DTF_mu1_PY", &DTF_mu1_PY );
  theTree->SetBranchAddress( "Dst_DTF_mu1_PZ", &DTF_mu1_PZ );

  theTree->SetBranchAddress("D_Hlt1TrackAllL0Decision_TIS",&D_Hlt1TrackAllL0Decision_TIS);
  theTree->SetBranchAddress("D_Hlt1L0AnyDecision_TIS",&D_Hlt1L0AnyDecision_TIS);
  theTree->SetBranchAddress("D_Hlt1GlobalDecision_TIS",&D_Hlt1GlobalDecision_TIS);
  theTree->SetBranchAddress("Dst_Hlt1Phys_TIS",&Dst_Hlt1Phys_TIS);
  theTree->SetBranchAddress("h0_RichDLLp",&h0_RichDLLp);
  theTree->SetBranchAddress("h1_RichDLLp",&h1_RichDLLp);
  theTree->SetBranchAddress("h1_ProbNNp",&h1_ProbNNp);
  theTree->SetBranchAddress("h0_ProbNNp",&h0_ProbNNp);
  theTree->SetBranchAddress("h0_TRACK_GhostProb",&h0_TRACK_GhostProb);
  theTree->SetBranchAddress("h1_TRACK_GhostProb",&h1_TRACK_GhostProb);
  theTree->SetBranchAddress("mu0_TRACK_GhostProb",&mu0_TRACK_GhostProb);
  theTree->SetBranchAddress("mu1_TRACK_GhostProb",&mu1_TRACK_GhostProb);
  theTree->SetBranchAddress("Slowpi_TRACK_GhostProb",&Slowpi_TRACK_GhostProb);
  
  //trigger                                                                                                                                                                       
  theTree->SetBranchAddress("mu0_L0MuonDecision_Dec",&mu0_L0MuonDecision_Dec);
  theTree->SetBranchAddress("mu1_L0MuonDecision_Dec",&mu1_L0MuonDecision_Dec);
  theTree->SetBranchAddress("mu0_L0DiMuonDecision_Dec",&mu0_L0DiMuonDecision_Dec);
  theTree->SetBranchAddress("mu1_L0DiMuonDecision_Dec",&mu1_L0DiMuonDecision_Dec);
  theTree->SetBranchAddress("h1_L0MuonDecision_Dec",&h1_L0MuonDecision_Dec);
  theTree->SetBranchAddress("h1_L0DiMuonDecision_Dec",&h1_L0DiMuonDecision_Dec);
  theTree->SetBranchAddress("Dst_L0Global_TIS",&Dst_L0Global_TIS);
  theTree->SetBranchAddress("D_L0Global_TIS",&D_L0Global_TIS);

  theTree->SetBranchAddress("D_L0MuonDecision_Dec",&D_L0MuonDecision_Dec);
  theTree->SetBranchAddress("D_L0DiMuonDecision_Dec",&D_L0DiMuonDecision_Dec);
  theTree->SetBranchAddress("D_L0PhotonDecision_Dec",&D_L0PhotonDecision_Dec);
  theTree->SetBranchAddress("D_L0ElectronDecision_Dec",&D_L0ElectronDecision_Dec);
  theTree->SetBranchAddress("D_L0HadronDecision_Dec",&D_L0HadronDecision_Dec);
  theTree->SetBranchAddress("Dst_L0MuonDecision_Dec",&Dst_L0MuonDecision_Dec);
  theTree->SetBranchAddress("Dst_L0DiMuonDecision_Dec",&Dst_L0DiMuonDecision_Dec);
  theTree->SetBranchAddress("Dst_L0PhotonDecision_Dec",&Dst_L0PhotonDecision_Dec);
  theTree->SetBranchAddress("Dst_L0ElectronDecision_Dec",&Dst_L0ElectronDecision_Dec);
  theTree->SetBranchAddress("Dst_L0HadronDecision_Dec",&Dst_L0HadronDecision_Dec);
  theTree->SetBranchAddress("runNumber",&runNumber);

  theTree->SetBranchAddress("mu0_Hlt1TrackMuonDecision_Dec",&mu0_Hlt1TrackMuonDecision_Dec);
  theTree->SetBranchAddress("mu1_Hlt1TrackMuonDecision_Dec",&mu1_Hlt1TrackMuonDecision_Dec);
  theTree->SetBranchAddress("D_Hlt1TrackAllL0Decision_Dec",&D_Hlt1TrackAllL0Decision_Dec);
  theTree->SetBranchAddress("D_Hlt1DiMuonHighMassDecision_Dec",&D_Hlt1DiMuonHighMassDecision_Dec);
  theTree->SetBranchAddress("D_Hlt1DiMuonLowMassDecision_Dec",&D_Hlt1DiMuonLowMassDecision_Dec);
  theTree->SetBranchAddress("mu0_Hlt1SingleMuonNoIPDecision_Dec",&mu0_Hlt1SingleMuonNoIPDecision_Dec);
  theTree->SetBranchAddress("mu1_Hlt1SingleMuonNoIPDecision_Dec",&mu1_Hlt1SingleMuonNoIPDecision_Dec);
  theTree->SetBranchAddress("mu0_Hlt1SingleMuonHighPTDecision_Dec",&mu0_Hlt1SingleMuonHighPTDecision_Dec);
  theTree->SetBranchAddress("mu1_Hlt1SingleMuonHighPTDecision_Dec",&mu1_Hlt1SingleMuonHighPTDecision_Dec);
  theTree->SetBranchAddress("h1_0_Hlt1TrackMuonDecision_Dec",&h1_0_Hlt1TrackMuonDecision_Dec);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02KKMuMuDecision_Dec",&D_Hlt2CharmSemilepD02KKMuMuDecision_Dec);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec",&D_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_Dec",&D_Hlt2CharmSemilepD02KPiMuMuDecision_Dec);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_Dec",&Dst_Hlt2CharmHadD02HHXDst_hhXDecision_Dec);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec",&Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec);
  theTree->SetBranchAddress("D_Hlt2DiMuonDetachedDecision_Dec",&D_Hlt2DiMuonDetachedDecision_Dec);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_Dec",&Dst_Hlt2CharmHadD02HHHHDst_4piDecision_Dec);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec",&Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec",&Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_K3piDecision_Dec",&D_Hlt2CharmHadD02HHHH_K3piDecision_Dec);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec",&D_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_4piDecision_Dec",&D_Hlt2CharmHadD02HHHH_4piDecision_Dec);


  theTree->SetBranchAddress("mu0_L0MuonDecision_TIS",&mu0_L0MuonDecision_TIS);
  theTree->SetBranchAddress("mu1_L0MuonDecision_TIS",&mu1_L0MuonDecision_TIS);
  theTree->SetBranchAddress("mu0_L0DiMuonDecision_TIS",&mu0_L0DiMuonDecision_TIS);
  theTree->SetBranchAddress("mu1_L0DiMuonDecision_TIS",&mu1_L0DiMuonDecision_TIS);
  theTree->SetBranchAddress("h1_L0MuonDecision_TIS",&h1_L0MuonDecision_TIS);
  theTree->SetBranchAddress("h1_L0DiMuonDecision_TIS",&h1_L0DiMuonDecision_TIS);
  theTree->SetBranchAddress("Dst_L0Global_TIS",&Dst_L0Global_TIS);
  theTree->SetBranchAddress("D_L0Global_TIS",&D_L0Global_TIS);
  theTree->SetBranchAddress("D_L0MuonDecision_TIS",&D_L0MuonDecision_TIS);
  theTree->SetBranchAddress("D_L0DiMuonDecision_TIS",&D_L0DiMuonDecision_TIS);
  theTree->SetBranchAddress("D_L0PhotonDecision_TIS",&D_L0PhotonDecision_TIS);
  theTree->SetBranchAddress("D_L0ElectronDecision_TIS",&D_L0ElectronDecision_TIS);
  theTree->SetBranchAddress("D_L0HadronDecision_TIS",&D_L0HadronDecision_TIS);
  theTree->SetBranchAddress("Dst_L0MuonDecision_TIS",&Dst_L0MuonDecision_TIS);
  theTree->SetBranchAddress("Dst_L0DiMuonDecision_TIS",&Dst_L0DiMuonDecision_TIS);
  theTree->SetBranchAddress("Dst_L0PhotonDecision_TIS",&Dst_L0PhotonDecision_TIS);
  theTree->SetBranchAddress("Dst_L0ElectronDecision_TIS",&Dst_L0ElectronDecision_TIS);
  theTree->SetBranchAddress("Dst_L0HadronDecision_TIS",&Dst_L0HadronDecision_TIS);

  theTree->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TIS",&mu0_Hlt1TrackMuonDecision_TIS);
  theTree->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TIS",&mu1_Hlt1TrackMuonDecision_TIS);
  theTree->SetBranchAddress("D_Hlt1TrackAllL0Decision_TIS",&D_Hlt1TrackAllL0Decision_TIS);
  theTree->SetBranchAddress("D_Hlt1DiMuonHighMassDecision_TIS",&D_Hlt1DiMuonHighMassDecision_TIS);
  theTree->SetBranchAddress("D_Hlt1DiMuonLowMassDecision_TIS",&D_Hlt1DiMuonLowMassDecision_TIS);
  theTree->SetBranchAddress("mu0_Hlt1SingleMuonNoIPDecision_TIS",&mu0_Hlt1SingleMuonNoIPDecision_TIS);
  theTree->SetBranchAddress("mu1_Hlt1SingleMuonNoIPDecision_TIS",&mu1_Hlt1SingleMuonNoIPDecision_TIS);
  theTree->SetBranchAddress("mu0_Hlt1SingleMuonHighPTDecision_TIS",&mu0_Hlt1SingleMuonHighPTDecision_TIS);
  theTree->SetBranchAddress("mu1_Hlt1SingleMuonHighPTDecision_TIS",&mu1_Hlt1SingleMuonHighPTDecision_TIS);
  theTree->SetBranchAddress("h1_0_Hlt1TrackMuonDecision_TIS",&h1_0_Hlt1TrackMuonDecision_TIS);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02KKMuMuDecision_TIS",&D_Hlt2CharmSemilepD02KKMuMuDecision_TIS);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS",&D_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_TIS",&D_Hlt2CharmSemilepD02KPiMuMuDecision_TIS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TIS",&Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TIS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS",&Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS);
  theTree->SetBranchAddress("D_Hlt2DiMuonDetachedDecision_TIS",&D_Hlt2DiMuonDetachedDecision_TIS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TIS",&Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TIS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS",&Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS",&Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_K3piDecision_TIS",&D_Hlt2CharmHadD02HHHH_K3piDecision_TIS);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS",&D_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_4piDecision_TIS",&D_Hlt2CharmHadD02HHHH_4piDecision_TIS);

  theTree->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  theTree->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  theTree->SetBranchAddress("mu0_L0DiMuonDecision_TOS",&mu0_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("mu1_L0DiMuonDecision_TOS",&mu1_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("h1_L0MuonDecision_TOS",&h1_L0MuonDecision_TOS);
  theTree->SetBranchAddress("h1_L0DiMuonDecision_TOS",&h1_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("D_L0MuonDecision_TOS",&D_L0MuonDecision_TOS);
  theTree->SetBranchAddress("D_L0DiMuonDecision_TOS",&D_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("D_L0PhotonDecision_TOS",&D_L0PhotonDecision_TOS);
  theTree->SetBranchAddress("D_L0ElectronDecision_TOS",&D_L0ElectronDecision_TOS);
  theTree->SetBranchAddress("D_L0HadronDecision_TOS",&D_L0HadronDecision_TOS);
  theTree->SetBranchAddress("Dst_L0MuonDecision_TOS",&Dst_L0MuonDecision_TOS);
  theTree->SetBranchAddress("Dst_L0DiMuonDecision_TOS",&Dst_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("Dst_L0PhotonDecision_TOS",&Dst_L0PhotonDecision_TOS);
  theTree->SetBranchAddress("Dst_L0ElectronDecision_TOS",&Dst_L0ElectronDecision_TOS);
  theTree->SetBranchAddress("Dst_L0HadronDecision_TOS",&Dst_L0HadronDecision_TOS);

  theTree->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  theTree->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  theTree->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);
  theTree->SetBranchAddress("D_Hlt1DiMuonHighMassDecision_TOS",&D_Hlt1DiMuonHighMassDecision_TOS);
  theTree->SetBranchAddress("D_Hlt1DiMuonLowMassDecision_TOS",&D_Hlt1DiMuonLowMassDecision_TOS);
  theTree->SetBranchAddress("mu0_Hlt1SingleMuonNoIPDecision_TOS",&mu0_Hlt1SingleMuonNoIPDecision_TOS);
  theTree->SetBranchAddress("mu1_Hlt1SingleMuonNoIPDecision_TOS",&mu1_Hlt1SingleMuonNoIPDecision_TOS);
  theTree->SetBranchAddress("mu0_Hlt1SingleMuonHighPTDecision_TOS",&mu0_Hlt1SingleMuonHighPTDecision_TOS);
  theTree->SetBranchAddress("mu1_Hlt1SingleMuonHighPTDecision_TOS",&mu1_Hlt1SingleMuonHighPTDecision_TOS);
  theTree->SetBranchAddress("h1_0_Hlt1TrackMuonDecision_TOS",&h1_0_Hlt1TrackMuonDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02KKMuMuDecision_TOS",&D_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS",&D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS",&D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS",&Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS",&Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2DiMuonDetachedDecision_TOS",&D_Hlt2DiMuonDetachedDecision_TOS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_K3piDecision_TOS",&D_Hlt2CharmHadD02HHHH_K3piDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS",&D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_4piDecision_TOS",&D_Hlt2CharmHadD02HHHH_4piDecision_TOS);

  theTree->SetBranchAddress( "Polarity",&Polarity);
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
  theTree->SetBranchAddress( "Dst_P", &user_Dst_P );
  theTree->SetBranchAddress( "D_P", &user_D_P );
  theTree->SetBranchAddress( "Dst_PZ", &user_Dst_PZ );
  theTree->SetBranchAddress( "D_PZ", &user_D_PZ );
  theTree->SetBranchAddress( "mu1_PT", &user_mu1_PT );
  theTree->SetBranchAddress( "mu0_PT", &user_mu0_PT );
  theTree->SetBranchAddress( "h0_PT", &user_h0_PT );
  theTree->SetBranchAddress( "h1_PT", &user_h1_PT );
  theTree->SetBranchAddress( "Slowpi_PT", &user_Slowpi_PT );
  theTree->SetBranchAddress( "mu1_P", &user_mu1_P );
  theTree->SetBranchAddress( "mu0_P", &user_mu0_P );
  theTree->SetBranchAddress( "h0_P", &user_h0_P );
  theTree->SetBranchAddress( "h1_P", &user_h1_P );
  theTree->SetBranchAddress( "Slowpi_P", &user_Slowpi_P );
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
  theTree->SetBranchAddress( "h1_PX", &user_h1_PX );
  theTree->SetBranchAddress( "h1_PY", &user_h1_PY );
  theTree->SetBranchAddress( "h1_PZ", &user_h1_PZ );
  theTree->SetBranchAddress( "mu0_PIDmu", &user_mu0_PIDmu );
  theTree->SetBranchAddress( "mu0_ProbNNmu", &user_mu0_ProbNNmu );
  theTree->SetBranchAddress( "mu1_ProbNNmu", &user_mu1_ProbNNmu );
  theTree->SetBranchAddress( "h0_ProbNNk", &user_h0_ProbNNk );
  theTree->SetBranchAddress( "h1_ProbNNk", &user_h1_ProbNNk );
  theTree->SetBranchAddress( "h0_ProbNNpi", &user_h0_ProbNNpi );
  theTree->SetBranchAddress( "h1_ProbNNpi", &user_h1_ProbNNpi );
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
  theTree->SetBranchAddress( "misID_mD_OS",&misID_mD_OS);
  theTree->SetBranchAddress( "misID_dm_OS",&misID_dm_OS);
  theTree->SetBranchAddress("mHH", & user_mHH);
  theTree->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
  theTree->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);
  theTree->SetBranchAddress("mu0_isMuon",&mu0_isMuon);
  theTree->SetBranchAddress("h0_PIDK",&h0_PIDK);
  theTree->SetBranchAddress("h1_PIDK",&h1_PIDK);
  theTree->SetBranchAddress("nCandidate",&nCandidate);
  theTree->SetBranchAddress("totCandidates",&totCandidates);

  if(!isMC) theTree->SetBranchAddress("nTracks_data", &user_nTracks_data);
  //MC truth informations 

  if(isMC){

   theTree->SetBranchAddress("Dst_BKGCAT", &Dst_BKGCAT);
   theTree->SetBranchAddress("Dst_TRUEID", &Dst_TRUEID);
   theTree->SetBranchAddress("Dst_MC_MOTHER_ID", &Dst_MC_MOTHER_ID);
   theTree->SetBranchAddress("Dst_MC_MOTHER_KEY", &Dst_MC_MOTHER_KEY );
   theTree->SetBranchAddress("Dst_TRUEP_E", &Dst_TRUEP_E);
   theTree->SetBranchAddress("Dst_TRUEP_X", &Dst_TRUEP_X);
   theTree->SetBranchAddress("Dst_TRUEP_Y", &Dst_TRUEP_Y);
   theTree->SetBranchAddress("Dst_TRUEP_Z", &Dst_TRUEP_Z);
   theTree->SetBranchAddress("Dst_TRUEPT", &Dst_TRUEPT);
   theTree->SetBranchAddress("Dst_TRUEORIGINVERTEX_X", &Dst_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("Dst_TRUEORIGINVERTEX_Y", &Dst_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("Dst_TRUEORIGINVERTEX_Z", &Dst_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("Dst_TRUEENDVERTEX_X", &Dst_TRUEENDVERTEX_X);
   theTree->SetBranchAddress("Dst_TRUEENDVERTEX_Y", &Dst_TRUEENDVERTEX_Y);
   theTree->SetBranchAddress("Dst_TRUEENDVERTEX_Z", &Dst_TRUEENDVERTEX_Z);
   theTree->SetBranchAddress("Dst_TRUEISSTABLE", &Dst_TRUEISSTABLE);
   theTree->SetBranchAddress("Dst_TRUETAU", &Dst_TRUETAU);
   theTree->SetBranchAddress("D_BKGCAT", &D_BKGCAT);
   theTree->SetBranchAddress("D_TRUEID", &D_TRUEID);
   theTree->SetBranchAddress("D_MC_MOTHER_ID", &D_MC_MOTHER_ID);
   theTree->SetBranchAddress("D_TRUEP_E", &D_TRUEP_E);
   theTree->SetBranchAddress("D_TRUEP_X", &D_TRUEP_X);
   theTree->SetBranchAddress("D_TRUEP_Y", &D_TRUEP_Y);
   theTree->SetBranchAddress("D_TRUEP_Z", &D_TRUEP_Z);
   theTree->SetBranchAddress("D_TRUEORIGINVERTEX_X", &D_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("D_TRUEORIGINVERTEX_Y", &D_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("D_TRUEORIGINVERTEX_Z", &D_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("D_TRUEENDVERTEX_X", &D_TRUEENDVERTEX_X);
   theTree->SetBranchAddress("D_TRUEENDVERTEX_Y", &D_TRUEENDVERTEX_Y);
   theTree->SetBranchAddress("D_TRUEENDVERTEX_Z", &D_TRUEENDVERTEX_Z);
   theTree->SetBranchAddress("D_TRUEISSTABLE", &D_TRUEISSTABLE);
   theTree->SetBranchAddress("D_TRUETAU", &D_TRUETAU);
   theTree->SetBranchAddress("h0_TRUEID", &h0_TRUEID);
   theTree->SetBranchAddress("h0_TRUEP_E", &h0_TRUEP_E);
   theTree->SetBranchAddress("h0_TRUEP_X", &h0_TRUEP_X);
   theTree->SetBranchAddress("h0_TRUEP_Y", &h0_TRUEP_Y);
   theTree->SetBranchAddress("h0_TRUEP_Z", &h0_TRUEP_Z);
   theTree->SetBranchAddress("h0_TRUEPT", &h0_TRUEPT);
   theTree->SetBranchAddress("h0_TRUEORIGINVERTEX_X", &h0_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("h0_TRUEORIGINVERTEX_Y", &h0_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("h0_TRUEORIGINVERTEX_Z", &h0_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("h0_TRUEISSTABLE", &h0_TRUEISSTABLE);
   theTree->SetBranchAddress("h0_TRUETAU", &h0_TRUETAU);
   theTree->SetBranchAddress("h1_TRUEID", &h1_TRUEID);
   theTree->SetBranchAddress("h1_TRUEP_E", &h1_TRUEP_E);
   theTree->SetBranchAddress("h1_TRUEP_X", &h1_TRUEP_X);
   theTree->SetBranchAddress("h1_TRUEP_Y", &h1_TRUEP_Y);
   theTree->SetBranchAddress("h1_TRUEP_Z", &h1_TRUEP_Z);
   theTree->SetBranchAddress("h1_TRUEPT", &h1_TRUEPT);
   theTree->SetBranchAddress("h1_TRUEORIGINVERTEX_X", &h1_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("h1_TRUEORIGINVERTEX_Y", &h1_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("h1_TRUEORIGINVERTEX_Z", &h1_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("h1_TRUETAU", &h1_TRUETAU );
   theTree->SetBranchAddress("mu0_TRUEID", &mu0_TRUEID);
   theTree->SetBranchAddress("mu0_TRUEP_E", &mu0_TRUEP_E);
   theTree->SetBranchAddress("mu0_TRUEP_X", &mu0_TRUEP_X);
   theTree->SetBranchAddress("mu0_TRUEP_Y", &mu0_TRUEP_Y);
   theTree->SetBranchAddress("mu0_TRUEP_Z", &mu0_TRUEP_Z);
   theTree->SetBranchAddress("mu0_TRUEPT", &mu0_TRUEPT);
   theTree->SetBranchAddress("mu0_TRUEORIGINVERTEX_X", &mu0_TRUEORIGINVERTEX_X );
   theTree->SetBranchAddress("mu0_TRUEORIGINVERTEX_Y", &mu0_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("mu0_TRUEORIGINVERTEX_Z", &mu0_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("mu0_TRUETAU", &mu0_TRUETAU);
   theTree->SetBranchAddress("mu1_TRUEID", &mu1_TRUEID );
   theTree->SetBranchAddress("mu1_MC_MOTHER_ID", &mu1_MC_MOTHER_ID );
   theTree->SetBranchAddress("mu1_TRUEP_E", &mu1_TRUEP_E);
   theTree->SetBranchAddress("mu1_TRUEP_X", &mu1_TRUEP_X);
   theTree->SetBranchAddress("mu1_TRUEP_Y", &mu1_TRUEP_Y);
   theTree->SetBranchAddress("mu1_TRUEP_Z", &mu1_TRUEP_Z);
   theTree->SetBranchAddress("mu1_TRUEPT", &mu1_TRUEPT);
   theTree->SetBranchAddress("mu1_TRUEORIGINVERTEX_X", &mu1_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("mu1_TRUEORIGINVERTEX_Y", &mu1_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("mu1_TRUEORIGINVERTEX_Z", &mu1_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("mu1_TRUEISSTABLE", &mu1_TRUEISSTABLE);
   theTree->SetBranchAddress("mu1_TRUETAU", &mu1_TRUETAU);
   theTree->SetBranchAddress("Slowpi_TRUEID", &Slowpi_TRUEID);
   theTree->SetBranchAddress("Slowpi_MC_MOTHER_ID", &Slowpi_MC_MOTHER_ID);
   theTree->SetBranchAddress("Slowpi_TRUEP_E", &Slowpi_TRUEP_E);
   theTree->SetBranchAddress("Slowpi_TRUEP_X", &Slowpi_TRUEP_X);
   theTree->SetBranchAddress("Slowpi_TRUEP_Y", &Slowpi_TRUEP_Y);
   theTree->SetBranchAddress("Slowpi_TRUEP_Z", &Slowpi_TRUEP_Z);
   theTree->SetBranchAddress("Slowpi_TRUEPT", &Slowpi_TRUEPT);
   theTree->SetBranchAddress("Slowpi_TRUEORIGINVERTEX_X", &Slowpi_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("Slowpi_TRUEORIGINVERTEX_Y", &Slowpi_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("Slowpi_TRUEORIGINVERTEX_Z", &Slowpi_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("Slowpi_TRUEISSTABLE", &Slowpi_TRUEISSTABLE);
   theTree->SetBranchAddress("Slowpi_TRUETAU", &Slowpi_TRUETAU);  
   theTree->SetBranchAddress("nTracks", &user_nTracks_data);
   theTree->SetBranchAddress("mu0_TRACK_Eta", &user_mu0_TRACK_Eta);
   theTree->SetBranchAddress("mu1_TRACK_Eta", &user_mu1_TRACK_Eta);
   theTree->SetBranchAddress("h1_TRACK_Eta", &user_h1_TRACK_Eta);
   theTree->SetBranchAddress("h0_TRACK_Eta", &user_h0_TRACK_Eta);
   theTree->SetBranchAddress("Slowpi_TRACK_Eta", &user_Slowpi_TRACK_Eta);
   

  

  }


  TFile *target  = new TFile(targetName,"RECREATE" );
  TTree* Tree = new TTree("BDT_Tree","BDT_Tree");         // new tree with BDT variable and all relevant other variables                                                        
  Tree->Branch("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Tree->Branch("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Tree->Branch("mu0_L0DiMuonDecision_TOS",&mu0_L0DiMuonDecision_TOS);
  Tree->Branch("mu1_L0DiMuonDecision_TOS",&mu1_L0DiMuonDecision_TOS);
  Tree->Branch("h1_L0MuonDecision_TOS",&h1_L0MuonDecision_TOS);
  Tree->Branch("h1_L0DiMuonDecision_TOS",&h1_L0DiMuonDecision_TOS);
  Tree->Branch("Dst_L0Global_TIS",&Dst_L0Global_TIS);
  Tree->Branch("D_L0Global_TIS",&D_L0Global_TIS);
  Tree->Branch("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  Tree->Branch("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  Tree->Branch("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);
  Tree->Branch("D_Hlt1DiMuonHighMassDecision_TOS",&D_Hlt1DiMuonHighMassDecision_TOS);
  Tree->Branch("D_Hlt1DiMuonLowMassDecision_TOS",&D_Hlt1DiMuonLowMassDecision_TOS);
  Tree->Branch("mu0_Hlt1SingleMuonNoIPDecision_TOS",&mu0_Hlt1SingleMuonNoIPDecision_TOS);
  Tree->Branch("mu1_Hlt1SingleMuonNoIPDecision_TOS",&mu1_Hlt1SingleMuonNoIPDecision_TOS);
  Tree->Branch("mu0_Hlt1SingleMuonHighPTDecision_TOS",&mu0_Hlt1SingleMuonHighPTDecision_TOS);
  Tree->Branch("mu1_Hlt1SingleMuonHighPTDecision_TOS",&mu1_Hlt1SingleMuonHighPTDecision_TOS);
  Tree->Branch("h1_Hlt1TrackMuonDecision_TOS",&h1_Hlt1TrackMuonDecision_TOS);
  Tree->Branch("D_Hlt2CharmSemilepD02KKMuMuDecision_TOS",&D_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
  Tree->Branch("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS",&D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS);
  Tree->Branch("D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS",&D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS",&Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS",&Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS);
  Tree->Branch("D_Hlt2DiMuonDetachedDecision_TOS",&D_Hlt2DiMuonDetachedDecision_TOS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_K3piDecision_TOS",&D_Hlt2CharmHadD02HHHH_K3piDecision_TOS);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS",&D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_4piDecision_Dec",&D_Hlt2CharmHadD02HHHH_4piDecision_Dec);

  Tree->Branch("mu0_L0MuonDecision_Dec",&mu0_L0MuonDecision_Dec);
  Tree->Branch("mu1_L0MuonDecision_Dec",&mu1_L0MuonDecision_Dec);
  Tree->Branch("mu0_L0DiMuonDecision_Dec",&mu0_L0DiMuonDecision_Dec);
  Tree->Branch("mu1_L0DiMuonDecision_Dec",&mu1_L0DiMuonDecision_Dec);
  Tree->Branch("h1_L0MuonDecision_Dec",&h1_L0MuonDecision_Dec);
  Tree->Branch("h1_L0DiMuonDecision_Dec",&h1_L0DiMuonDecision_Dec);
  Tree->Branch("Dst_L0Global_TIS",&Dst_L0Global_TIS);
  Tree->Branch("D_L0Global_TIS",&D_L0Global_TIS);
  Tree->Branch("mu0_Hlt1TrackMuonDecision_Dec",&mu0_Hlt1TrackMuonDecision_Dec);
  Tree->Branch("mu1_Hlt1TrackMuonDecision_Dec",&mu1_Hlt1TrackMuonDecision_Dec);
  Tree->Branch("D_Hlt1TrackAllL0Decision_Dec",&D_Hlt1TrackAllL0Decision_Dec);
  Tree->Branch("D_Hlt1DiMuonHighMassDecision_Dec",&D_Hlt1DiMuonHighMassDecision_Dec);
  Tree->Branch("D_Hlt1DiMuonLowMassDecision_Dec",&D_Hlt1DiMuonLowMassDecision_Dec);
  Tree->Branch("mu0_Hlt1SingleMuonNoIPDecision_Dec",&mu0_Hlt1SingleMuonNoIPDecision_Dec);
  Tree->Branch("mu1_Hlt1SingleMuonNoIPDecision_Dec",&mu1_Hlt1SingleMuonNoIPDecision_Dec);
  Tree->Branch("mu0_Hlt1SingleMuonHighPTDecision_Dec",&mu0_Hlt1SingleMuonHighPTDecision_Dec);
  Tree->Branch("mu1_Hlt1SingleMuonHighPTDecision_Dec",&mu1_Hlt1SingleMuonHighPTDecision_Dec);
  Tree->Branch("h1_Hlt1TrackMuonDecision_Dec",&h1_Hlt1TrackMuonDecision_Dec);
  Tree->Branch("D_Hlt2CharmSemilepD02KKMuMuDecision_Dec",&D_Hlt2CharmSemilepD02KKMuMuDecision_Dec);
  Tree->Branch("D_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec",&D_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec);
  Tree->Branch("D_Hlt2CharmSemilepD02KPiMuMuDecision_Dec",&D_Hlt2CharmSemilepD02KPiMuMuDecision_Dec);
  Tree->Branch("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_Dec",&Dst_Hlt2CharmHadD02HHXDst_hhXDecision_Dec);
  Tree->Branch("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec",&Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec);
  Tree->Branch("D_Hlt2DiMuonDetachedDecision_Dec",&D_Hlt2DiMuonDetachedDecision_Dec);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_Dec",&Dst_Hlt2CharmHadD02HHHHDst_4piDecision_Dec);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec",&Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec",&Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_K3piDecision_Dec",&D_Hlt2CharmHadD02HHHH_K3piDecision_Dec);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec",&D_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_4piDecision_Dec",&D_Hlt2CharmHadD02HHHH_4piDecision_Dec);

  Tree->Branch("mu0_L0MuonDecision_TIS",&mu0_L0MuonDecision_TIS);
  Tree->Branch("mu1_L0MuonDecision_TIS",&mu1_L0MuonDecision_TIS);
  Tree->Branch("mu0_L0DiMuonDecision_TIS",&mu0_L0DiMuonDecision_TIS);
  Tree->Branch("mu1_L0DiMuonDecision_TIS",&mu1_L0DiMuonDecision_TIS);
  Tree->Branch("h1_L0MuonDecision_TIS",&h1_L0MuonDecision_TIS);
  Tree->Branch("h1_L0DiMuonDecision_TIS",&h1_L0DiMuonDecision_TIS);
  Tree->Branch("Dst_L0Global_TIS",&Dst_L0Global_TIS);
  Tree->Branch("D_L0Global_TIS",&D_L0Global_TIS);
  Tree->Branch("mu0_Hlt1TrackMuonDecision_TIS",&mu0_Hlt1TrackMuonDecision_TIS);
  Tree->Branch("mu1_Hlt1TrackMuonDecision_TIS",&mu1_Hlt1TrackMuonDecision_TIS);
  Tree->Branch("D_Hlt1TrackAllL0Decision_TIS",&D_Hlt1TrackAllL0Decision_TIS);
  Tree->Branch("D_Hlt1DiMuonHighMassDecision_TIS",&D_Hlt1DiMuonHighMassDecision_TIS);
  Tree->Branch("D_Hlt1DiMuonLowMassDecision_TIS",&D_Hlt1DiMuonLowMassDecision_TIS);
  Tree->Branch("mu0_Hlt1SingleMuonNoIPDecision_TIS",&mu0_Hlt1SingleMuonNoIPDecision_TIS);
  Tree->Branch("mu1_Hlt1SingleMuonNoIPDecision_TIS",&mu1_Hlt1SingleMuonNoIPDecision_TIS);
  Tree->Branch("mu0_Hlt1SingleMuonHighPTDecision_TIS",&mu0_Hlt1SingleMuonHighPTDecision_TIS);
  Tree->Branch("mu1_Hlt1SingleMuonHighPTDecision_TIS",&mu1_Hlt1SingleMuonHighPTDecision_TIS);
  Tree->Branch("h1_Hlt1TrackMuonDecision_TIS",&h1_Hlt1TrackMuonDecision_TIS);
  Tree->Branch("D_Hlt2CharmSemilepD02KKMuMuDecision_TIS",&D_Hlt2CharmSemilepD02KKMuMuDecision_TIS);
  Tree->Branch("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS",&D_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS);
  Tree->Branch("D_Hlt2CharmSemilepD02KPiMuMuDecision_TIS",&D_Hlt2CharmSemilepD02KPiMuMuDecision_TIS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TIS",&Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TIS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS",&Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS);
  Tree->Branch("D_Hlt2DiMuonDetachedDecision_TIS",&D_Hlt2DiMuonDetachedDecision_TIS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TIS",&Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TIS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS",&Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS",&Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_K3piDecision_TIS",&D_Hlt2CharmHadD02HHHH_K3piDecision_TIS);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS",&D_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_4piDecision_TIS",&D_Hlt2CharmHadD02HHHH_4piDecision_TIS);

  Tree->Branch("D_L0MuonDecision_TIS",&D_L0MuonDecision_TIS);
  Tree->Branch("D_L0DiMuonDecision_TIS",&D_L0DiMuonDecision_TIS);
  Tree->Branch("D_L0PhotonDecision_TIS",&D_L0PhotonDecision_TIS);
  Tree->Branch("D_L0ElectronDecision_TIS",&D_L0ElectronDecision_TIS);
  Tree->Branch("D_L0HadronDecision_TIS",&D_L0HadronDecision_TIS);
  Tree->Branch("Dst_L0MuonDecision_TIS",&Dst_L0MuonDecision_TIS);
  Tree->Branch("Dst_L0DiMuonDecision_TIS",&Dst_L0DiMuonDecision_TIS);
  Tree->Branch("Dst_L0PhotonDecision_TIS",&Dst_L0PhotonDecision_TIS);
  Tree->Branch("Dst_L0ElectronDecision_TIS",&Dst_L0ElectronDecision_TIS);
  Tree->Branch("Dst_L0HadronDecision_TIS",&Dst_L0HadronDecision_TIS);

  Tree->Branch("D_L0MuonDecision_Dec",&D_L0MuonDecision_Dec);
  Tree->Branch("D_L0DiMuonDecision_Dec",&D_L0DiMuonDecision_Dec);
  Tree->Branch("D_L0PhotonDecision_Dec",&D_L0PhotonDecision_Dec);
  Tree->Branch("D_L0ElectronDecision_Dec",&D_L0ElectronDecision_Dec);
  Tree->Branch("D_L0HadronDecision_Dec",&D_L0HadronDecision_Dec);
  Tree->Branch("Dst_L0MuonDecision_Dec",&Dst_L0MuonDecision_Dec);
  Tree->Branch("Dst_L0DiMuonDecision_Dec",&Dst_L0DiMuonDecision_Dec);
  Tree->Branch("Dst_L0PhotonDecision_Dec",&Dst_L0PhotonDecision_Dec);
  Tree->Branch("Dst_L0ElectronDecision_Dec",&Dst_L0ElectronDecision_Dec);
  Tree->Branch("Dst_L0HadronDecision_Dec",&Dst_L0HadronDecision_Dec);

  Tree->Branch("D_L0MuonDecision_TOS",&D_L0MuonDecision_TOS);
  Tree->Branch("D_L0DiMuonDecision_TOS",&D_L0DiMuonDecision_TOS);
  Tree->Branch("D_L0PhotonDecision_TOS",&D_L0PhotonDecision_TOS);
  Tree->Branch("D_L0ElectronDecision_TOS",&D_L0ElectronDecision_TOS);
  Tree->Branch("D_L0HadronDecision_TOS",&D_L0HadronDecision_TOS);
  Tree->Branch("Dst_L0MuonDecision_TOS",&Dst_L0MuonDecision_TOS);
  Tree->Branch("Dst_L0DiMuonDecision_TOS",&Dst_L0DiMuonDecision_TOS);
  Tree->Branch("Dst_L0PhotonDecision_TOS",&Dst_L0PhotonDecision_TOS);
  Tree->Branch("Dst_L0ElectronDecision_TOS",&Dst_L0ElectronDecision_TOS);
  Tree->Branch("Dst_L0HadronDecision_TOS",&Dst_L0HadronDecision_TOS);
  Tree->Branch( "runNumber",&runNumber);

  Tree->Branch( "Polarity",&Polarity);
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
  Tree->Branch( "Dst_P", &user_Dst_P );
  Tree->Branch( "D_P", &user_D_P );
  Tree->Branch( "Dst_PZ", &user_Dst_PZ );
  Tree->Branch( "D_PZ", &user_D_PZ );
  Tree->Branch( "mu1_PT", &user_mu1_PT );
  Tree->Branch( "mu0_PT", &user_mu0_PT );
  Tree->Branch( "h0_PT", &user_h0_PT );
  Tree->Branch( "h1_PT", &user_h1_PT );
  Tree->Branch( "mu1_P", &user_mu1_P );
  Tree->Branch( "mu0_P", &user_mu0_P );
  Tree->Branch( "h0_P", &user_h0_P );
  Tree->Branch( "h1_P", &user_h1_P );
  Tree->Branch( "Slowpi_P", &user_Slowpi_P );
  Tree->Branch( "mu0_PX", &user_mu0_PX );
  Tree->Branch( "mu0_PY", &user_mu0_PY );
  Tree->Branch( "mu0_PZ", &user_mu0_PZ );
  Tree->Branch( "mu1_PX", &user_mu1_PX );
  Tree->Branch( "mu1_PY", &user_mu1_PY );
  Tree->Branch( "mu1_PZ", &user_mu1_PZ );
  Tree->Branch( "h0_PX", &user_h0_PX );
  Tree->Branch( "h0_PY", &user_h0_PY );
  Tree->Branch( "h0_PZ", &user_h0_PZ );
  Tree->Branch( "h1_PX", &user_h1_PX );
  Tree->Branch( "h1_PY", &user_h1_PY );
  Tree->Branch( "h1_PZ", &user_h1_PZ );
  Tree->Branch( "min_muPT", &min_muPT );
  Tree->Branch( "max_muPT", &max_muPT );
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
  Tree->Branch( "misID_dm_OS", &misID_dm_OS );
  Tree->Branch( "misID_mD_OS", &misID_mD_OS );
  Tree->Branch( "mHH", & user_mHH);
  Tree->Branch( "Dst_Coneptasy", & user_Dst_Coneptasy);
  Tree->Branch( "h0_ProbNNk", &user_h0_ProbNNk );
  Tree->Branch( "h1_ProbNNk", &user_h1_ProbNNk );
  Tree->Branch( "h0_ProbNNpi", &user_h0_ProbNNpi );
  Tree->Branch( "h1_ProbNNpi", &user_h1_ProbNNpi );
  Tree->Branch( "nTracks", &user_nTracks_data );
  Tree->Branch("mu1_MuonNShared",&mu1_MuonNShared);
  Tree->Branch("mu0_MuonNShared",&mu0_MuonNShared);
  Tree->Branch("mu0_isMuon",&mu0_isMuon);
  Tree->Branch("mu1_isMuon",&mu1_isMuon);
  Tree->Branch( "D_ENDVERTEX_CHI2",&user_D_ENDVERTEX_CHI2);
  Tree->Branch( "D_FDCHI2_OWNPV",&user_D_FDCHI2_OWNPV );
  Tree->Branch("h0_PIDK", &h0_PIDK);  
  Tree->Branch("h1_PIDK", &h1_PIDK);  
  Tree->Branch("mpiKpipi", &mpiKpipi);  
  Tree->Branch("mumu_P", &mumu_P);
  Tree->Branch("mumu_PT", &mumu_PT);
  Tree->Branch("nCandidate",&nCandidate);
  Tree->Branch("totCandidates",&totCandidates);

  Tree->Branch("D_Hlt1TrackAllL0Decision_TIS",&D_Hlt1TrackAllL0Decision_TIS);
  Tree->Branch("D_Hlt1L0AnyDecision_TIS,",&D_Hlt1L0AnyDecision_TIS);
  Tree->Branch("D_Hlt1GlobalDecision_TIS",&D_Hlt1GlobalDecision_TIS);
  Tree->Branch("Dst_Hlt1Phys_TIS",&Dst_Hlt1Phys_TIS);
  Tree->Branch("h0_RichDLLp",&h0_RichDLLp);
  Tree->Branch("h1_RichDLLp",&h1_RichDLLp);
  Tree->Branch("h1_ProbNNp",&h1_ProbNNp);
  Tree->Branch("h0_ProbNNp",&h0_ProbNNp);
  Tree->Branch("h0_TRACK_GhostProb",&h0_TRACK_GhostProb);
  Tree->Branch("h1_TRACK_GhostProb",&h1_TRACK_GhostProb);
  Tree->Branch("mu0_TRACK_GhostProb",&mu0_TRACK_GhostProb);
  Tree->Branch("mu1_TRACK_GhostProb",&mu1_TRACK_GhostProb);
  Tree->Branch("Slowpi_TRACK_GhostProb",&Slowpi_TRACK_GhostProb);
  Tree->Branch("mKPiMuMu",&mKPiMuMu);
  Tree->Branch("mKMuPiMu",&mKMuPiMu);
  Tree->Branch("mKMuMuPi",&mKMuMuPi);
  
  Tree->Branch("m34",&m34);
  Tree->Branch("m24",&m24);
  Tree->Branch("m23",&m23);


  if(isMC) {

   Tree->Branch("Dst_BKGCAT", &Dst_BKGCAT);
   Tree->Branch("Dst_TRUEID", &Dst_TRUEID);
   Tree->Branch("Dst_MC_MOTHER_ID", &Dst_MC_MOTHER_ID);
   Tree->Branch("Dst_MC_MOTHER_KEY", &Dst_MC_MOTHER_KEY );
   Tree->Branch("Dst_TRUEP_E", &Dst_TRUEP_E);
   Tree->Branch("Dst_TRUEP_X", &Dst_TRUEP_X);
   Tree->Branch("Dst_TRUEP_Y", &Dst_TRUEP_Y);
   Tree->Branch("Dst_TRUEP_Z", &Dst_TRUEP_Z);
   Tree->Branch("Dst_TRUEPT", &Dst_TRUEPT);
   Tree->Branch("Dst_TRUEORIGINVERTEX_X", &Dst_TRUEORIGINVERTEX_X);
   Tree->Branch("Dst_TRUEORIGINVERTEX_Y", &Dst_TRUEORIGINVERTEX_Y);
   Tree->Branch("Dst_TRUEORIGINVERTEX_Z", &Dst_TRUEORIGINVERTEX_Z);
   Tree->Branch("Dst_TRUEENDVERTEX_X", &Dst_TRUEENDVERTEX_X);
   Tree->Branch("Dst_TRUEENDVERTEX_Y", &Dst_TRUEENDVERTEX_Y);
   Tree->Branch("Dst_TRUEENDVERTEX_Z", &Dst_TRUEENDVERTEX_Z);
   Tree->Branch("Dst_TRUEISSTABLE", &Dst_TRUEISSTABLE);
   Tree->Branch("Dst_TRUETAU", &Dst_TRUETAU);
   Tree->Branch("D_BKGCAT", &D_BKGCAT);
   Tree->Branch("D_TRUEID", &D_TRUEID);
   Tree->Branch("D_MC_MOTHER_ID", &D_MC_MOTHER_ID);
   Tree->Branch("D_TRUEP_E", &D_TRUEP_E);
   Tree->Branch("D_TRUEP_X", &D_TRUEP_X);
   Tree->Branch("D_TRUEP_Y", &D_TRUEP_Y);
   Tree->Branch("D_TRUEP_Z", &D_TRUEP_Z);
   Tree->Branch("D_TRUE_DiMuon_Mass", &D_TRUE_DiMuon_Mass);
   Tree->Branch("D_TRUEORIGINVERTEX_X", &D_TRUEORIGINVERTEX_X);
   Tree->Branch("D_TRUEORIGINVERTEX_Y", &D_TRUEORIGINVERTEX_Y);
   Tree->Branch("D_TRUEORIGINVERTEX_Z", &D_TRUEORIGINVERTEX_Z);
   Tree->Branch("D_TRUEENDVERTEX_X", &D_TRUEENDVERTEX_X);
   Tree->Branch("D_TRUEENDVERTEX_Y", &D_TRUEENDVERTEX_Y);
   Tree->Branch("D_TRUEENDVERTEX_Z", &D_TRUEENDVERTEX_Z);
   Tree->Branch("D_TRUEISSTABLE", &D_TRUEISSTABLE);
   Tree->Branch("D_TRUETAU", &D_TRUETAU);
   Tree->Branch("h0_TRUEID", &h0_TRUEID);
   Tree->Branch("h0_TRUEP_E", &h0_TRUEP_E);
   Tree->Branch("h0_TRUEP_X", &h0_TRUEP_X);
   Tree->Branch("h0_TRUEP_Y", &h0_TRUEP_Y);
   Tree->Branch("h0_TRUEP_Z", &h0_TRUEP_Z);
   Tree->Branch("h0_TRUEPT", &h0_TRUEPT);
   Tree->Branch("h0_TRUEORIGINVERTEX_X", &h0_TRUEORIGINVERTEX_X);
   Tree->Branch("h0_TRUEORIGINVERTEX_Y", &h0_TRUEORIGINVERTEX_Y);
   Tree->Branch("h0_TRUEORIGINVERTEX_Z", &h0_TRUEORIGINVERTEX_Z);
   Tree->Branch("h0_TRUEISSTABLE", &h0_TRUEISSTABLE);
   Tree->Branch("h0_TRUETAU", &h0_TRUETAU);
   Tree->Branch("h1_TRUEID", &h1_TRUEID);
   Tree->Branch("h1_TRUEP_E", &h1_TRUEP_E);
   Tree->Branch("h1_TRUEP_X", &h1_TRUEP_X);
   Tree->Branch("h1_TRUEP_Y", &h1_TRUEP_Y);
   Tree->Branch("h1_TRUEP_Z", &h1_TRUEP_Z);
   Tree->Branch("h1_TRUEPT", &h1_TRUEPT);
   Tree->Branch("h1_TRUEORIGINVERTEX_X", &h1_TRUEORIGINVERTEX_X);
   Tree->Branch("h1_TRUEORIGINVERTEX_Y", &h1_TRUEORIGINVERTEX_Y);
   Tree->Branch("h1_TRUEORIGINVERTEX_Z", &h1_TRUEORIGINVERTEX_Z);
   Tree->Branch("h1_TRUETAU", &h1_TRUETAU);
   Tree->Branch("mu0_TRUEID", &mu0_TRUEID);
   Tree->Branch("mu0_TRUEP_E", &mu0_TRUEP_E);
   Tree->Branch("mu0_TRUEP_X", &mu0_TRUEP_X);
   Tree->Branch("mu0_TRUEP_Y", &mu0_TRUEP_Y);
   Tree->Branch("mu0_TRUEP_Z", &mu0_TRUEP_Z);
   Tree->Branch("mu0_TRUEPT", &mu0_TRUEPT);
   Tree->Branch("mu0_TRUEORIGINVERTEX_X", &mu0_TRUEORIGINVERTEX_X);
   Tree->Branch("mu0_TRUEORIGINVERTEX_Y", &mu0_TRUEORIGINVERTEX_Y);
   Tree->Branch("mu0_TRUEORIGINVERTEX_Z", &mu0_TRUEORIGINVERTEX_Z);
   Tree->Branch("mu0_TRUETAU", &mu0_TRUETAU);
   Tree->Branch("mu1_TRUEID", &mu1_TRUEID );
   Tree->Branch("mu1_MC_MOTHER_ID", &mu1_MC_MOTHER_ID );
   Tree->Branch("mu1_TRUEP_E", &mu1_TRUEP_E);
   Tree->Branch("mu1_TRUEP_X", &mu1_TRUEP_X);
   Tree->Branch("mu1_TRUEP_Y", &mu1_TRUEP_Y);
   Tree->Branch("mu1_TRUEP_Z", &mu1_TRUEP_Z);
   Tree->Branch("mu1_TRUEPT", &mu1_TRUEPT);
   Tree->Branch("mu1_TRUEORIGINVERTEX_X", &mu1_TRUEORIGINVERTEX_X);
   Tree->Branch("mu1_TRUEORIGINVERTEX_Y", &mu1_TRUEORIGINVERTEX_Y);
   Tree->Branch("mu1_TRUEORIGINVERTEX_Z", &mu1_TRUEORIGINVERTEX_Z);
   Tree->Branch("mu1_TRUEISSTABLE", &mu1_TRUEISSTABLE);
   Tree->Branch("mu1_TRUETAU", &mu1_TRUETAU);
   Tree->Branch("Slowpi_TRUEID", &Slowpi_TRUEID);
   Tree->Branch("Slowpi_MC_MOTHER_ID", &Slowpi_MC_MOTHER_ID);
   Tree->Branch("Slowpi_TRUEP_E", &Slowpi_TRUEP_E);
   Tree->Branch("Slowpi_TRUEP_X", &Slowpi_TRUEP_X);
   Tree->Branch("Slowpi_TRUEP_Y", &Slowpi_TRUEP_Y);
   Tree->Branch("Slowpi_TRUEP_Z", &Slowpi_TRUEP_Z);
   Tree->Branch("Slowpi_TRUEPT", &Slowpi_TRUEPT);
   Tree->Branch("Slowpi_TRUEORIGINVERTEX_X", &Slowpi_TRUEORIGINVERTEX_X);
   Tree->Branch("Slowpi_TRUEORIGINVERTEX_Y", &Slowpi_TRUEORIGINVERTEX_Y);
   Tree->Branch("Slowpi_TRUEORIGINVERTEX_Z", &Slowpi_TRUEORIGINVERTEX_Z);
   Tree->Branch("Slowpi_TRUEISSTABLE", &Slowpi_TRUEISSTABLE);
   Tree->Branch("Slowpi_TRUETAU", &Slowpi_TRUETAU);  
   Tree->Branch("mu0_TRACK_Eta", &user_mu0_TRACK_Eta);
   Tree->Branch("mu1_TRACK_Eta", &user_mu1_TRACK_Eta);
   Tree->Branch("h1_TRACK_Eta", &user_h1_TRACK_Eta);
   Tree->Branch("h0_TRACK_Eta", &user_h0_TRACK_Eta);
   Tree->Branch("Slowpi_TRACK_Eta", &user_Slowpi_TRACK_Eta);

  }

  TLorentzVector h0_hypK,h0_hypPi,h0_hypMu;
  TLorentzVector h1_hypK,h1_hypPi,h1_hypMu;
  TLorentzVector mu0_hypK,mu0_hypPi,mu0_hypMu;
  TLorentzVector mu1_hypK,mu1_hypPi,mu1_hypMu;


  std::cout << "number Events: " << theTree->GetEntries() << std::endl;

  for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

    if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
    //if (ievt>60000) break;                                                                                                                                                    
    theTree->GetEntry(ievt);
    ////ALL CUTS REMOVED RIGHT NOW!!!////////////
    
    //if(user_Slowpi_ProbNNghost>0.5 || user_h1_ProbNNghost>0.5 ||
    // user_h0_ProbNNghost>0.5 || user_mu1_ProbNNghost>0.5 ||  user_mu0_ProbNNghost > 0.5) continue;
    
    //cut on nMuonShared
    if( mu0_MuonNShared!=0 || mu1_MuonNShared!=0 ) continue;
        
    if(skipPID)  goto skipped; //..again sorry for the goto... but is simplest solution
      
    //// HADRON PID, different for signal and normalization modes /////                                                                                                                     
    if(isNormalizationMode && ( user_h0_ProbNNk<0.2 || user_h1_ProbNNpi<0.2) ) continue;
    if(!isNormalizationMode && ( user_h0_ProbNNk<0.2 || user_h1_ProbNNk<0.2)) continue;

  skipped:

    mu0_hypMu.SetXYZM(mu0_TRUEP_X,mu0_TRUEP_Y,mu0_TRUEP_Z,Mass::Mu());
    mu1_hypMu.SetXYZM(mu1_TRUEP_X,mu1_TRUEP_Y,mu1_TRUEP_Z,Mass::Mu());
    D_TRUE_DiMuon_Mass=(mu0_hypMu+mu1_hypMu).M();

    mumu_P=(mu0_hypMu+mu1_hypMu).P();
    mumu_PT=(mu0_hypMu+mu1_hypMu).Pt();

    h1_hypK.SetXYZM(DTF_h1_PX,DTF_h1_PY,DTF_h1_PZ,Mass::K());
    h1_hypPi.SetXYZM(DTF_h1_PX,DTF_h1_PY,DTF_h1_PZ,Mass::Pi());
    h1_hypMu.SetXYZM(DTF_h1_PX,DTF_h1_PY,DTF_h1_PZ,Mass::Mu());
    h0_hypPi.SetXYZM(DTF_h0_PX,DTF_h0_PY,DTF_h0_PZ,Mass::Pi());
    h0_hypK.SetXYZM(DTF_h0_PX,DTF_h0_PY,DTF_h0_PZ,Mass::K());
    h0_hypMu.SetXYZM(DTF_h0_PX,DTF_h0_PY,DTF_h0_PZ,Mass::Mu());
    mu0_hypPi.SetXYZM(DTF_mu0_PX,DTF_mu0_PY,DTF_mu0_PZ,Mass::Pi());
    mu1_hypPi.SetXYZM(DTF_mu1_PX,DTF_mu1_PY,DTF_mu1_PZ,Mass::Pi());

    mpiKpipi=(h0_hypPi+h1_hypK+ mu0_hypPi + mu1_hypPi).M();

    //recompute Di_DiuonMass (needed especially for misID)                                                                                                                         
    mu0_hypMu.SetXYZM(DTF_mu0_PX,DTF_mu0_PY,DTF_mu0_PZ,Mass::Mu());
    mu1_hypMu.SetXYZM(DTF_mu1_PX,DTF_mu1_PY,DTF_mu1_PZ,Mass::Mu());
    user_D_DiMuon_Mass = (mu0_hypMu+mu1_hypMu).M();

    mKPiMuMu = ( h0_hypK + h1_hypPi + mu0_hypMu + mu1_hypMu).M();
    mKMuPiMu = ( h0_hypK + h1_hypMu + mu0_hypPi + mu1_hypMu).M();
    mKMuMuPi = ( h0_hypK + h1_hypMu + mu0_hypMu + mu1_hypPi).M();
    
    m34 = (mu0_hypMu + mu1_hypMu).M();
    m24 = (h1_hypMu + mu1_hypMu).M();
    m23 = (h1_hypMu + mu0_hypMu).M();

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
    log_D_DIRA_OWNPV=float( TMath::Log(user_D_DIRA_OWNPV));
    min_mu_PT = float(TMath::Min(user_mu1_PT,user_mu0_PT));
    min_h_PT = float(TMath::Min(user_h1_PT,user_h0_PT));
    max_mu_PT = float(TMath::Max(user_mu1_PT,user_mu0_PT));
    Slowpi_P = float(user_Slowpi_P);
    Slowpi_PT = float(user_Slowpi_PT);
    D_MINIP= float(user_D_MINIP);
    D_MINIPCHI2= float(user_D_MINIPCHI2);
    D_ENDVERTEX_CHI2 = float(user_D_ENDVERTEX_CHI2);
    min_muPT = TMath::Min(user_mu1_PT,user_mu0_PT);
    max_muPT = TMath::Max(user_mu1_PT,user_mu0_PT);


    user_BDT= reader->EvaluateMVA("BDTG method");
    //if(BDT<0) continue;
    //if(ievt>50000)break;                                                                                                                                                      

    Tree->Fill();
  }

  target->Write();

  delete reader;
  //delete theTree;                                                                                                                                                                   

  target->Close();
  //target->Delete();                                          

}
void Application_D2pipimumu(TString treeName, TString fileIn, TString fileOut, int part,bool isMC = false, bool isNormalizationMode=false, bool skipPID=false) {


  TString dir    = "weights/";
  TString prefix;
  TString weightfile;
  TString methodName="BDTG method";
  TString dataDir = "/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/";
  TString targetName =  dataDir + fileOut;

  //select even trained BDT for part 1 and odd for part 2                                                                                                                             
  if(part==1){
    prefix = "training_D2pipimumu_evenTrained.root";
    weightfile = dir + prefix + TString("_")  + TString("BDTG") + TString(".weights.xml");
  }

  if(part==2){
    prefix = "training_D2pipimumu_oddTrained.root";
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
  //reader->AddVariable("D_cosh",&D_cosh);
  //reader->AddVariable("mu0_cosh",&mu0_cosh);
  reader->AddVariable("log_D_FDCHI2_OWNPV:=log10(D_FDCHI2_OWNPV)",&log_D_FDCHI2_OWNPV);
  reader->AddVariable("log_D_DIRA_OWNPV:=log(D_DIRA_OWNPV)",&log_D_DIRA_OWNPV);
  //reader->AddVariable("D_DIRA_OWNPV",&D_DIRA_OWNPV);                                                                                                            
  reader->AddVariable("D_ENDVERTEX_CHI2",&D_ENDVERTEX_CHI2);
  reader->AddVariable("log_Slowpi_IPCHI2_OWNPV:=log10(Slowpi_IPCHI2_OWNPV)",&log_Slowpi_IPCHI2_OWNPV);
  reader->AddVariable("D_MINIPCHI2",&D_MINIPCHI2);
  reader->AddVariable("Slowpi_P",&Slowpi_P);
  //reader->AddVariable("min_mu_PT:=min(mu0_PT,mu1_PT)",&min_mu_PT);///
  //reader->AddVariable("min_h_PT:=min(h0_PT,h1_PT)",&min_h_PT);///
  reader->AddVariable("Dst_Coneptasy",&Dst_Coneptasy);
  reader->AddVariable("Slowpi_PT",&Slowpi_PT);                                                                                                                                      
  //reader->AddVariable("max_mu_PT:=max(mu0_PT,mu1_PT)",&max_mu_PT);///                        

  reader->BookMVA( methodName, weightfile );

  // read variables of the data tree and create new tree with BDT variable                                                                                                            

  Double_t user_Dst_MAXDOCA,user_D_MAXDOCA, user_Dst_MINIP, user_D_FD_OWNPV, user_D_FDCHI2_OWNPV, user_D_DIRA_OWNPV;
  Double_t user_mu0_IPCHI2_OWNPV , user_mu1_IPCHI2_OWNPV , user_h0_IPCHI2_OWNPV , user_h1_IPCHI2_OWNPV , user_Slowpi_IPCHI2_OWNPV;
  Double_t user_D_MINIP, user_D_MINIPCHI2;
  Double_t user_Dst_PT, user_D_PT, user_mu1_PT, user_mu0_PT, user_h0_PT,user_h1_PT, user_Slowpi_PT;
  Double_t user_mu0_PX, user_mu0_PY, user_mu0_PZ, user_mu1_PX, user_mu1_PY, user_mu1_PZ,user_h0_PX, user_h0_PY, user_h0_PZ,user_h1_PX, user_h1_PY,user_h1_PZ;
  Double_t user_mu0_P, user_mu1_P, user_h0_P,user_h1_P;
  Double_t user_mu0_PIDmu, user_mu1_PIDmu;
  Double_t user_D_M, user_Dst_M, user_D_DiMuon_Mass;
  Double_t user_Slowpi_cosh,user_D_cosh, user_mu0_cosh, user_D_ENDVERTEX_CHI2,user_Slowpi_P;
  Double_t user_Dst_Coneptasy;
  Double_t user_Slowpi_ProbNNghost, user_mu0_ProbNNghost, user_mu1_ProbNNghost,user_h0_ProbNNghost,user_h1_ProbNNghost;
  Double_t user_BDT;
  Double_t user_deltaM, user_Dst_DTF_D0_M;
  Double_t user_mu1_ProbNNmu, user_mu0_ProbNNmu;
  Double_t user_D_TAU;
  ULong64_t user_eventNumber;
  Int_t user_nTracks_data;
  Int_t user_nTracks_MC;
  Double_t misID_mD_OS=-1000;
  Double_t misID_dm_OS=-1000;
  Double_t h0_PIDK, h1_PIDK; 
  double user_mu0_TRACK_Eta;
  double user_mu1_TRACK_Eta;
  double user_h0_TRACK_Eta;
  double user_h1_TRACK_Eta;
  double user_Slowpi_TRACK_Eta;  
  UInt_t          nCandidate;
  ULong64_t       totCandidates;
  UInt_t          runNumber;
  double h0_RichDLLp,h1_RichDLLp;
  double h1_ProbNNp,h0_ProbNNp;

  double h0_TRACK_GhostProb;
  double h1_TRACK_GhostProb;
  double mu0_TRACK_GhostProb;
  double mu1_TRACK_GhostProb;
  double Slowpi_TRACK_GhostProb;
  double mKPiMuMu,mKMuPiMu,mKMuMuPi;
  double mPiPiMuMu,mMuMuPiPi,mPiMuPiMu,mPiMuMuPi,mMuPiPiMu,mMuPiMuPi,mPiPiPiPiPi;
  double m34,m24,m23,m12,m13,m14;

  
  bool D_Hlt1TrackAllL0Decision_TIS, D_Hlt1L0AnyDecision_TIS, D_Hlt1GlobalDecision_TIS,Dst_Hlt1Phys_TIS;

  bool mu0_L0MuonDecision_TOS,mu1_L0MuonDecision_TOS,mu0_L0DiMuonDecision_TOS,mu1_L0DiMuonDecision_TOS,h1_L0MuonDecision_TOS,h1_L0DiMuonDecision_TOS,Dst_L0Global_TIS,D_L0Global_TIS,mu0_Hlt1TrackMuonDecision_TOS,mu1_Hlt1TrackMuonDecision_TOS,D_Hlt1TrackAllL0Decision_TOS,D_Hlt1DiMuonHighMassDecision_TOS,D_Hlt1DiMuonLowMassDecision_TOS,mu0_Hlt1SingleMuonNoIPDecision_TOS,mu1_Hlt1SingleMuonNoIPDecision_TOS,mu0_Hlt1SingleMuonHighPTDecision_TOS,mu1_Hlt1SingleMuonHighPTDecision_TOS,h1_0_Hlt1TrackMuonDecision_TOS,D_Hlt2CharmSemilepD02KKMuMuDecision_TOS,D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS,D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS,Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS,Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS,D_Hlt2DiMuonDetachedDecision_TOS,Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS,Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS,Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS,D_Hlt2CharmHadD02HHHH_K3piDecision_TOS,D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS,D_Hlt2CharmHadD02HHHH_4piDecision_TOS,h1_Hlt1TrackMuonDecision_TOS;

  bool D_L0MuonDecision_TIS,D_L0DiMuonDecision_TIS,D_L0PhotonDecision_TIS,D_L0ElectronDecision_TIS,D_L0HadronDecision_TOS,D_L0MuonDecision_TOS,D_L0DiMuonDecision_TOS,D_L0PhotonDecision_TOS,D_L0ElectronDecision_TOS,mu0_L0HadronDecision_TIS,mu0_L0MuonDecision_TIS,mu0_L0DiMuonDecision_TIS,mu0_L0PhotonDecision_TIS,mu0_L0ElectronDecision_TIS,mu0_L0HadronDecision_TOS,mu0_L0PhotonDecision_TOS,mu0_L0ElectronDecision_TOS,mu1_L0HadronDecision_TIS,mu1_L0MuonDecision_TIS,mu1_L0DiMuonDecision_TIS,mu1_L0PhotonDecision_TIS,mu1_L0ElectronDecision_TIS,mu1_L0HadronDecision_TOS,mu1_L0PhotonDecision_TOS,mu1_L0ElectronDecision_TOS,D_L0HadronDecision_TIS;

  //MC Variables                                                                                                                                                            
  Int_t           Dst_BKGCAT;   //!                                                                                                                                      
  Int_t           Dst_TRUEID;   //!                                                                                                                                      
  Int_t           Dst_MC_MOTHER_ID;
  Int_t           Dst_MC_MOTHER_KEY;
  Int_t           Dst_MC_GD_MOTHER_ID;
  Int_t           Dst_MC_GD_MOTHER_KEY;
  Int_t           Dst_MC_GD_GD_MOTHER_ID;
  Int_t           Dst_MC_GD_GD_MOTHER_KEY;
  Double_t        Dst_TRUEP_E;
  Double_t        Dst_TRUEP_X;
  Double_t        Dst_TRUEP_Y;
  Double_t        Dst_TRUEP_Z;
  Double_t        Dst_TRUEPT;
  Double_t        Dst_TRUEORIGINVERTEX_X;
  Double_t        Dst_TRUEORIGINVERTEX_Y;
  Double_t        Dst_TRUEORIGINVERTEX_Z;
  Double_t        Dst_TRUEENDVERTEX_X;
  Double_t        Dst_TRUEENDVERTEX_Y;
  Double_t        Dst_TRUEENDVERTEX_Z;
  Bool_t          Dst_TRUEISSTABLE;
  Double_t        Dst_TRUETAU;
  Int_t           D_BKGCAT;
  Int_t           D_TRUEID;
  Int_t           D_MC_MOTHER_ID;
  Int_t           D_MC_MOTHER_KEY;
  Double_t        D_TRUEP_E;
  Double_t        D_TRUEP_X;
  Double_t        D_TRUEP_Y;
  Double_t        D_TRUEP_Z;
  Double_t        D_TRUEPT;
  Double_t        D_TRUEORIGINVERTEX_X;
  Double_t        D_TRUEORIGINVERTEX_Y;
  Double_t        D_TRUEORIGINVERTEX_Z;
  Double_t        D_TRUEENDVERTEX_X;
  Double_t        D_TRUEENDVERTEX_Y;
  Double_t        D_TRUEENDVERTEX_Z;
  Bool_t          D_TRUEISSTABLE;
  Double_t        D_TRUETAU;
  Int_t           h0_TRUEID;
  Int_t           h0_MC_MOTHER_ID;
  Double_t        h0_TRUEP_E;
  Double_t        h0_TRUEP_X;
  Double_t        h0_TRUEP_Y;
  Double_t        h0_TRUEP_Z;
  Double_t        h0_TRUEPT;
  Double_t        h0_TRUEORIGINVERTEX_X;
  Double_t        h0_TRUEORIGINVERTEX_Y;
  Double_t        h0_TRUEORIGINVERTEX_Z;
  Bool_t          h0_TRUEISSTABLE;
  Double_t        h0_TRUETAU;
  Int_t           h1_TRUEID;
  Int_t           h1_MC_MOTHER_ID;
  Double_t        h1_TRUEP_E;
  Double_t        h1_TRUEP_X;
  Double_t        h1_TRUEP_Y;
  Double_t        h1_TRUEP_Z;
  Double_t        h1_TRUEPT;
  Double_t        h1_TRUEORIGINVERTEX_X;
  Double_t        h1_TRUEORIGINVERTEX_Y;
  Double_t        h1_TRUEORIGINVERTEX_Z;
  Bool_t          h1_TRUEISSTABLE;
  Double_t        h1_TRUETAU;
  Int_t           mu0_TRUEID;
  Int_t           mu0_MC_MOTHER_ID;
  Double_t        mu0_TRUEP_E;
  Double_t        mu0_TRUEP_X;
  Double_t        mu0_TRUEP_Y;
  Double_t        mu0_TRUEP_Z;
  Double_t        mu0_TRUEPT;
  Double_t        mu0_TRUEORIGINVERTEX_X;
  Double_t        mu0_TRUEORIGINVERTEX_Y;
  Double_t        mu0_TRUEORIGINVERTEX_Z;
  Double_t        mu0_TRUEENDVERTEX_X;
  Double_t        mu0_TRUEENDVERTEX_Y;
  Double_t        mu0_TRUEENDVERTEX_Z;
  Bool_t          mu0_TRUEISSTABLE;
  Double_t        mu0_TRUETAU;
  Int_t           mu1_TRUEID;
  Int_t           mu1_MC_MOTHER_ID;
  Int_t           mu1_MC_MOTHER_KEY;
  Double_t        mu1_TRUEP_E;
  Double_t        mu1_TRUEP_X;
  Double_t        mu1_TRUEP_Y;
  Double_t        mu1_TRUEP_Z;
  Double_t        mu1_TRUEPT;
  Double_t        mu1_TRUEORIGINVERTEX_X;
  Double_t        mu1_TRUEORIGINVERTEX_Y;
  Double_t        mu1_TRUEORIGINVERTEX_Z;
  Bool_t          mu1_TRUEISSTABLE;
  Double_t        mu1_TRUETAU;
  Int_t           Slowpi_TRUEID;
  Int_t           Slowpi_MC_MOTHER_ID;
  Double_t        Slowpi_TRUEP_E;
  Double_t        Slowpi_TRUEP_X;
  Double_t        Slowpi_TRUEP_Y;
  Double_t        Slowpi_TRUEP_Z;
  Double_t        Slowpi_TRUEPT;
  Double_t        Slowpi_TRUEORIGINVERTEX_X;
  Double_t        Slowpi_TRUEORIGINVERTEX_Y;
  Double_t        Slowpi_TRUEORIGINVERTEX_Z;
  Bool_t          Slowpi_TRUEISSTABLE;
  Double_t        Slowpi_TRUETAU;
  Double_t        user_mHH;
  Double_t user_h0_ProbNNk,user_h1_ProbNNk,user_h0_ProbNNpi,user_h1_ProbNNpi,user_Slowpi_ProbNNpi;
  Double_t D_TRUE_DiMuon_Mass;
  Int_t           mu1_MuonNShared;
  Int_t           mu0_MuonNShared;
  Bool_t          mu0_isMuon;
  Bool_t          mu1_isMuon;
  Short_t         Polarity;
  //double          h0_ProbNNp;
  //double          h1_ProbNNp;
  double DTF_h0_PX,DTF_h0_PY,DTF_h0_PZ;
  double DTF_h1_PX,DTF_h1_PY,DTF_h1_PZ;
  double DTF_mu0_PX,DTF_mu0_PY,DTF_mu0_PZ;
  double DTF_mu1_PX,DTF_mu1_PY,DTF_mu1_PZ;
  double user_Slowpi_PX,user_Slowpi_PY,user_Slowpi_PZ;
  double DTF_Pis_PX,DTF_Pis_PY,DTF_Pis_PZ;

  theTree->SetBranchAddress( "Dst_DTF_h0_PX", &DTF_h0_PX);
  theTree->SetBranchAddress( "Dst_DTF_h0_PY", &DTF_h0_PY );
  theTree->SetBranchAddress( "Dst_DTF_h0_PZ", &DTF_h0_PZ );
  theTree->SetBranchAddress( "Dst_DTF_h1_PX", &DTF_h1_PX);
  theTree->SetBranchAddress( "Dst_DTF_h1_PY", &DTF_h1_PY );
  theTree->SetBranchAddress( "Dst_DTF_h1_PZ", &DTF_h1_PZ );
  theTree->SetBranchAddress( "Dst_DTF_mu0_PX", &DTF_mu0_PX);
  theTree->SetBranchAddress( "Dst_DTF_mu0_PY", &DTF_mu0_PY );
  theTree->SetBranchAddress( "Dst_DTF_mu0_PZ", &DTF_mu0_PZ );
  theTree->SetBranchAddress( "Dst_DTF_mu1_PX", &DTF_mu1_PX);
  theTree->SetBranchAddress( "Dst_DTF_mu1_PY", &DTF_mu1_PY );
  theTree->SetBranchAddress( "Dst_DTF_mu1_PZ", &DTF_mu1_PZ );
  theTree->SetBranchAddress( "Dst_DTF_Pis_PX", &DTF_Pis_PX);
  theTree->SetBranchAddress( "Dst_DTF_Pis_PY", &DTF_Pis_PY );
  theTree->SetBranchAddress( "Dst_DTF_Pis_PZ", &DTF_Pis_PZ );

  theTree->SetBranchAddress("D_Hlt1TrackAllL0Decision_TIS",&D_Hlt1TrackAllL0Decision_TIS);
  theTree->SetBranchAddress("D_Hlt1L0AnyDecision_TIS,",&D_Hlt1L0AnyDecision_TIS);
  theTree->SetBranchAddress("D_Hlt1GlobalDecision_TIS",&D_Hlt1GlobalDecision_TIS);
  theTree->SetBranchAddress("Dst_Hlt1Phys_TIS",&Dst_Hlt1Phys_TIS);
  theTree->SetBranchAddress("h0_RichDLLp",&h0_RichDLLp);
  theTree->SetBranchAddress("h1_RichDLLp",&h1_RichDLLp);
  theTree->SetBranchAddress("h1_ProbNNp",&h1_ProbNNp);
  theTree->SetBranchAddress("h0_ProbNNp",&h0_ProbNNp);
  theTree->SetBranchAddress("h0_TRACK_GhostProb",&h0_TRACK_GhostProb);
  theTree->SetBranchAddress("h1_TRACK_GhostProb",&h1_TRACK_GhostProb);
  theTree->SetBranchAddress("mu0_TRACK_GhostProb",&mu0_TRACK_GhostProb);
  theTree->SetBranchAddress("mu1_TRACK_GhostProb",&mu1_TRACK_GhostProb);
  theTree->SetBranchAddress("Slowpi_TRACK_GhostProb",&Slowpi_TRACK_GhostProb);


  //trigger                                                                                                                                                                       

  theTree->SetBranchAddress("Dst_L0Global_TIS",&Dst_L0Global_TIS);

  theTree->SetBranchAddress("D_L0HadronDecision_TIS",&D_L0HadronDecision_TIS);
  theTree->SetBranchAddress("D_L0MuonDecision_TIS",&D_L0MuonDecision_TIS);
  theTree->SetBranchAddress("D_L0DiMuonDecision_TIS",&D_L0DiMuonDecision_TIS);
  theTree->SetBranchAddress("D_L0PhotonDecision_TIS",&D_L0PhotonDecision_TIS);
  theTree->SetBranchAddress("D_L0ElectronDecision_TIS",&D_L0ElectronDecision_TIS);
  theTree->SetBranchAddress("D_L0HadronDecision_TOS",&D_L0HadronDecision_TOS);
  theTree->SetBranchAddress("D_L0MuonDecision_TOS",&D_L0MuonDecision_TOS);
  theTree->SetBranchAddress("D_L0DiMuonDecision_TOS",&D_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("D_L0PhotonDecision_TOS",&D_L0PhotonDecision_TOS);
  theTree->SetBranchAddress("D_L0ElectronDecision_TOS",&D_L0ElectronDecision_TOS);
 
  theTree->SetBranchAddress("mu0_L0HadronDecision_TIS",&mu0_L0HadronDecision_TIS);
  theTree->SetBranchAddress("mu0_L0MuonDecision_TIS",&mu0_L0MuonDecision_TIS);
  theTree->SetBranchAddress("mu0_L0DiMuonDecision_TIS",&mu0_L0DiMuonDecision_TIS);
  theTree->SetBranchAddress("mu0_L0PhotonDecision_TIS",&mu0_L0PhotonDecision_TIS);
  theTree->SetBranchAddress("mu0_L0ElectronDecision_TIS",&mu0_L0ElectronDecision_TIS);
  theTree->SetBranchAddress("mu0_L0HadronDecision_TOS",&mu0_L0HadronDecision_TOS);
  theTree->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  theTree->SetBranchAddress("mu0_L0DiMuonDecision_TOS",&mu0_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("mu0_L0PhotonDecision_TOS",&mu0_L0PhotonDecision_TOS);
  theTree->SetBranchAddress("mu0_L0ElectronDecision_TOS",&mu0_L0ElectronDecision_TOS);

  theTree->SetBranchAddress("mu1_L0HadronDecision_TIS",&mu1_L0HadronDecision_TIS);
  theTree->SetBranchAddress("mu1_L0MuonDecision_TIS",&mu1_L0MuonDecision_TIS);
  theTree->SetBranchAddress("mu1_L0DiMuonDecision_TIS",&mu1_L0DiMuonDecision_TIS);
  theTree->SetBranchAddress("mu1_L0PhotonDecision_TIS",&mu1_L0PhotonDecision_TIS);
  theTree->SetBranchAddress("mu1_L0ElectronDecision_TIS",&mu1_L0ElectronDecision_TIS);
  theTree->SetBranchAddress("mu1_L0HadronDecision_TOS",&mu1_L0HadronDecision_TOS);
  theTree->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  theTree->SetBranchAddress("mu1_L0DiMuonDecision_TOS",&mu1_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("mu1_L0PhotonDecision_TOS",&mu1_L0PhotonDecision_TOS);
  theTree->SetBranchAddress("mu1_L0ElectronDecision_TOS",&mu1_L0ElectronDecision_TOS);

  theTree->SetBranchAddress("h0_L0HadronDecision_TIS",&mu0_L0HadronDecision_TIS);
  theTree->SetBranchAddress("h0_L0MuonDecision_TIS",&mu0_L0MuonDecision_TIS);
  theTree->SetBranchAddress("h0_L0DiMuonDecision_TIS",&mu0_L0DiMuonDecision_TIS);
  theTree->SetBranchAddress("h0_L0PhotonDecision_TIS",&mu0_L0PhotonDecision_TIS);
  theTree->SetBranchAddress("h0_L0ElectronDecision_TIS",&mu0_L0ElectronDecision_TIS);
  theTree->SetBranchAddress("h0_L0HadronDecision_TOS",&mu0_L0HadronDecision_TOS);
  theTree->SetBranchAddress("h0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  theTree->SetBranchAddress("h0_L0DiMuonDecision_TOS",&mu0_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("h0_L0PhotonDecision_TOS",&mu0_L0PhotonDecision_TOS);
  theTree->SetBranchAddress("h0_L0ElectronDecision_TOS",&mu0_L0ElectronDecision_TOS);

  theTree->SetBranchAddress("h1_L0HadronDecision_TIS",&mu1_L0HadronDecision_TIS);
  theTree->SetBranchAddress("h1_L0MuonDecision_TIS",&mu1_L0MuonDecision_TIS);
  theTree->SetBranchAddress("h1_L0DiMuonDecision_TIS",&mu1_L0DiMuonDecision_TIS);
  theTree->SetBranchAddress("h1_L0PhotonDecision_TIS",&mu1_L0PhotonDecision_TIS);
  theTree->SetBranchAddress("h1_L0ElectronDecision_TIS",&mu1_L0ElectronDecision_TIS);
  theTree->SetBranchAddress("h1_L0HadronDecision_TOS",&mu1_L0HadronDecision_TOS);
  theTree->SetBranchAddress("h1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  theTree->SetBranchAddress("h1_L0DiMuonDecision_TOS",&mu1_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("h1_L0PhotonDecision_TOS",&mu1_L0PhotonDecision_TOS);
  theTree->SetBranchAddress("h1_L0ElectronDecision_TOS",&mu1_L0ElectronDecision_TOS);


  theTree->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  theTree->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  theTree->SetBranchAddress("mu0_L0DiMuonDecision_TOS",&mu0_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("mu1_L0DiMuonDecision_TOS",&mu1_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("h1_L0MuonDecision_TOS",&h1_L0MuonDecision_TOS);
  theTree->SetBranchAddress("h1_L0DiMuonDecision_TOS",&h1_L0DiMuonDecision_TOS);
  theTree->SetBranchAddress("Dst_L0Global_TIS",&Dst_L0Global_TIS);
  
  theTree->SetBranchAddress("D_L0Global_TIS",&D_L0Global_TIS);
  theTree->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  theTree->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  theTree->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);
  theTree->SetBranchAddress("D_Hlt1DiMuonHighMassDecision_TOS",&D_Hlt1DiMuonHighMassDecision_TOS);
  theTree->SetBranchAddress("D_Hlt1DiMuonLowMassDecision_TOS",&D_Hlt1DiMuonLowMassDecision_TOS);
  theTree->SetBranchAddress("mu0_Hlt1SingleMuonNoIPDecision_TOS",&mu0_Hlt1SingleMuonNoIPDecision_TOS);
  theTree->SetBranchAddress("mu1_Hlt1SingleMuonNoIPDecision_TOS",&mu1_Hlt1SingleMuonNoIPDecision_TOS);
  theTree->SetBranchAddress("mu0_Hlt1SingleMuonHighPTDecision_TOS",&mu0_Hlt1SingleMuonHighPTDecision_TOS);
  theTree->SetBranchAddress("mu1_Hlt1SingleMuonHighPTDecision_TOS",&mu1_Hlt1SingleMuonHighPTDecision_TOS);
  theTree->SetBranchAddress("h1_0_Hlt1TrackMuonDecision_TOS",&h1_0_Hlt1TrackMuonDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02KKMuMuDecision_TOS",&D_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS",&D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS",&D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS",&Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS",&Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2DiMuonDetachedDecision_TOS",&D_Hlt2DiMuonDetachedDecision_TOS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS);
  theTree->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_K3piDecision_TOS",&D_Hlt2CharmHadD02HHHH_K3piDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS",&D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS);
  theTree->SetBranchAddress("D_Hlt2CharmHadD02HHHH_4piDecision_TOS",&D_Hlt2CharmHadD02HHHH_4piDecision_TOS);
  theTree->SetBranchAddress( "runNumber",&runNumber);

  theTree->SetBranchAddress( "Polarity",&Polarity);
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
  theTree->SetBranchAddress( "mu1_P", &user_mu1_P );
  theTree->SetBranchAddress( "mu0_P", &user_mu0_P );
  theTree->SetBranchAddress( "h0_P", &user_h0_P );
  theTree->SetBranchAddress( "h1_P", &user_h1_P );
  theTree->SetBranchAddress( "Slowpi_PT", &user_Slowpi_PT );
  theTree->SetBranchAddress( "mu0_PX", &user_mu0_PX );
  theTree->SetBranchAddress( "mu0_PY", &user_mu0_PY );
  theTree->SetBranchAddress( "mu0_PZ", &user_mu0_PZ );
  theTree->SetBranchAddress( "Slowpi_PX", &user_Slowpi_PX );
  theTree->SetBranchAddress( "Slowpi_PY", &user_Slowpi_PY );
  theTree->SetBranchAddress( "Slowpi_PZ", &user_Slowpi_PZ );
  theTree->SetBranchAddress( "mu1_PX", &user_mu1_PX );
  theTree->SetBranchAddress( "mu1_PY", &user_mu1_PY );
  theTree->SetBranchAddress( "mu1_PIDmu", &user_mu1_PIDmu );
  theTree->SetBranchAddress( "mu1_PZ", &user_mu1_PZ );
  theTree->SetBranchAddress( "h0_PX", &user_h0_PX );
  theTree->SetBranchAddress( "h0_PY", &user_h0_PY );
  theTree->SetBranchAddress( "h0_PZ", &user_h0_PZ );
  theTree->SetBranchAddress( "h1_PX", &user_h1_PX );
  theTree->SetBranchAddress( "h1_PY", &user_h1_PY );
  theTree->SetBranchAddress( "h1_PZ", &user_h1_PZ );
  theTree->SetBranchAddress( "mu0_PIDmu", &user_mu0_PIDmu );
  theTree->SetBranchAddress( "mu0_ProbNNmu", &user_mu0_ProbNNmu );
  theTree->SetBranchAddress( "mu1_ProbNNmu", &user_mu1_ProbNNmu );
  theTree->SetBranchAddress( "h0_ProbNNk", &user_h0_ProbNNk );
  theTree->SetBranchAddress( "h1_ProbNNk", &user_h1_ProbNNk );
  theTree->SetBranchAddress( "h0_ProbNNpi", &user_h0_ProbNNpi );
  theTree->SetBranchAddress( "h1_ProbNNpi", &user_h1_ProbNNpi );
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
  theTree->SetBranchAddress( "Slowpi_ProbNNpi", & user_Slowpi_ProbNNpi);
  theTree->SetBranchAddress( "deltaM", & user_deltaM);
  theTree->SetBranchAddress( "Dst_DTF_D0_M", & user_Dst_DTF_D0_M);
  theTree->SetBranchAddress( "D_TAU",&user_D_TAU);
  theTree->SetBranchAddress( "D_cosh",&user_D_cosh);
  theTree->SetBranchAddress( "mu0_cosh",&user_mu0_cosh);
  theTree->SetBranchAddress( "eventNumber",&user_eventNumber);
  theTree->SetBranchAddress( "misID_mD_OS",&misID_mD_OS);
  theTree->SetBranchAddress( "misID_dm_OS",&misID_dm_OS);
  theTree->SetBranchAddress("mHH", & user_mHH);
  theTree->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
  theTree->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);
  theTree->SetBranchAddress("mu0_isMuon",&mu0_isMuon);
  theTree->SetBranchAddress("h0_PIDK",&h0_PIDK);
  theTree->SetBranchAddress("h1_PIDK",&h1_PIDK);
  theTree->SetBranchAddress("nCandidate",&nCandidate);
  theTree->SetBranchAddress("totCandidates",&totCandidates);

  theTree->SetBranchAddress("mu0_TRACK_Eta", &user_mu0_TRACK_Eta);
  theTree->SetBranchAddress("mu1_TRACK_Eta", &user_mu1_TRACK_Eta);
  theTree->SetBranchAddress("h1_TRACK_Eta", &user_h1_TRACK_Eta);
  theTree->SetBranchAddress("h0_TRACK_Eta", &user_h0_TRACK_Eta);
  theTree->SetBranchAddress("Slowpi_TRACK_Eta", &user_Slowpi_TRACK_Eta);

  if(!isMC) theTree->SetBranchAddress("nTracks_data", &user_nTracks_data);
  //MC truth informations 

  if(isMC){

   theTree->SetBranchAddress("Dst_BKGCAT", &Dst_BKGCAT);
   theTree->SetBranchAddress("Dst_TRUEID", &Dst_TRUEID);
   theTree->SetBranchAddress("Dst_MC_MOTHER_ID", &Dst_MC_MOTHER_ID);
   theTree->SetBranchAddress("Dst_MC_MOTHER_KEY", &Dst_MC_MOTHER_KEY );
   theTree->SetBranchAddress("Dst_TRUEP_E", &Dst_TRUEP_E);
   theTree->SetBranchAddress("Dst_TRUEP_X", &Dst_TRUEP_X);
   theTree->SetBranchAddress("Dst_TRUEP_Y", &Dst_TRUEP_Y);
   theTree->SetBranchAddress("Dst_TRUEP_Z", &Dst_TRUEP_Z);
   theTree->SetBranchAddress("Dst_TRUEPT", &Dst_TRUEPT);
   theTree->SetBranchAddress("Dst_TRUEORIGINVERTEX_X", &Dst_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("Dst_TRUEORIGINVERTEX_Y", &Dst_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("Dst_TRUEORIGINVERTEX_Z", &Dst_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("Dst_TRUEENDVERTEX_X", &Dst_TRUEENDVERTEX_X);
   theTree->SetBranchAddress("Dst_TRUEENDVERTEX_Y", &Dst_TRUEENDVERTEX_Y);
   theTree->SetBranchAddress("Dst_TRUEENDVERTEX_Z", &Dst_TRUEENDVERTEX_Z);
   theTree->SetBranchAddress("Dst_TRUEISSTABLE", &Dst_TRUEISSTABLE);
   theTree->SetBranchAddress("Dst_TRUETAU", &Dst_TRUETAU);
   theTree->SetBranchAddress("D_BKGCAT", &D_BKGCAT);
   theTree->SetBranchAddress("D_TRUEID", &D_TRUEID);
   theTree->SetBranchAddress("D_MC_MOTHER_ID", &D_MC_MOTHER_ID);
   theTree->SetBranchAddress("D_TRUEP_E", &D_TRUEP_E);
   theTree->SetBranchAddress("D_TRUEP_X", &D_TRUEP_X);
   theTree->SetBranchAddress("D_TRUEP_Y", &D_TRUEP_Y);
   theTree->SetBranchAddress("D_TRUEP_Z", &D_TRUEP_Z);
   theTree->SetBranchAddress("D_TRUEORIGINVERTEX_X", &D_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("D_TRUEORIGINVERTEX_Y", &D_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("D_TRUEORIGINVERTEX_Z", &D_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("D_TRUEENDVERTEX_X", &D_TRUEENDVERTEX_X);
   theTree->SetBranchAddress("D_TRUEENDVERTEX_Y", &D_TRUEENDVERTEX_Y);
   theTree->SetBranchAddress("D_TRUEENDVERTEX_Z", &D_TRUEENDVERTEX_Z);
   theTree->SetBranchAddress("D_TRUEISSTABLE", &D_TRUEISSTABLE);
   theTree->SetBranchAddress("D_TRUETAU", &D_TRUETAU);
   theTree->SetBranchAddress("h0_TRUEID", &h0_TRUEID);
   theTree->SetBranchAddress("h0_TRUEP_E", &h0_TRUEP_E);
   theTree->SetBranchAddress("h0_TRUEP_X", &h0_TRUEP_X);
   theTree->SetBranchAddress("h0_TRUEP_Y", &h0_TRUEP_Y);
   theTree->SetBranchAddress("h0_TRUEP_Z", &h0_TRUEP_Z);
   theTree->SetBranchAddress("h0_TRUEPT", &h0_TRUEPT);
   theTree->SetBranchAddress("h0_TRUEORIGINVERTEX_X", &h0_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("h0_TRUEORIGINVERTEX_Y", &h0_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("h0_TRUEORIGINVERTEX_Z", &h0_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("h0_TRUEISSTABLE", &h0_TRUEISSTABLE);
   theTree->SetBranchAddress("h0_TRUETAU", &h0_TRUETAU);
   theTree->SetBranchAddress("h1_TRUEID", &h1_TRUEID);
   theTree->SetBranchAddress("h1_TRUEP_E", &h1_TRUEP_E);
   theTree->SetBranchAddress("h1_TRUEP_X", &h1_TRUEP_X);
   theTree->SetBranchAddress("h1_TRUEP_Y", &h1_TRUEP_Y);
   theTree->SetBranchAddress("h1_TRUEP_Z", &h1_TRUEP_Z);
   theTree->SetBranchAddress("h1_TRUEPT", &h1_TRUEPT);
   theTree->SetBranchAddress("h1_TRUEORIGINVERTEX_X", &h1_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("h1_TRUEORIGINVERTEX_Y", &h1_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("h1_TRUEORIGINVERTEX_Z", &h1_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("h1_TRUETAU", &h1_TRUETAU );
   theTree->SetBranchAddress("mu0_TRUEID", &mu0_TRUEID);
   theTree->SetBranchAddress("mu0_TRUEP_E", &mu0_TRUEP_E);
   theTree->SetBranchAddress("mu0_TRUEP_X", &mu0_TRUEP_X);
   theTree->SetBranchAddress("mu0_TRUEP_Y", &mu0_TRUEP_Y);
   theTree->SetBranchAddress("mu0_TRUEP_Z", &mu0_TRUEP_Z);
   theTree->SetBranchAddress("mu0_TRUEPT", &mu0_TRUEPT);
   theTree->SetBranchAddress("mu0_TRUEORIGINVERTEX_X", &mu0_TRUEORIGINVERTEX_X );
   theTree->SetBranchAddress("mu0_TRUEORIGINVERTEX_Y", &mu0_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("mu0_TRUEORIGINVERTEX_Z", &mu0_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("mu0_TRUETAU", &mu0_TRUETAU);
   theTree->SetBranchAddress("mu1_TRUEID", &mu1_TRUEID );
   theTree->SetBranchAddress("mu1_MC_MOTHER_ID", &mu1_MC_MOTHER_ID );
   theTree->SetBranchAddress("mu1_TRUEP_E", &mu1_TRUEP_E);
   theTree->SetBranchAddress("mu1_TRUEP_X", &mu1_TRUEP_X);
   theTree->SetBranchAddress("mu1_TRUEP_Y", &mu1_TRUEP_Y);
   theTree->SetBranchAddress("mu1_TRUEP_Z", &mu1_TRUEP_Z);
   theTree->SetBranchAddress("mu1_TRUEPT", &mu1_TRUEPT);
   theTree->SetBranchAddress("mu1_TRUEORIGINVERTEX_X", &mu1_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("mu1_TRUEORIGINVERTEX_Y", &mu1_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("mu1_TRUEORIGINVERTEX_Z", &mu1_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("mu1_TRUEISSTABLE", &mu1_TRUEISSTABLE);
   theTree->SetBranchAddress("mu1_TRUETAU", &mu1_TRUETAU);
   theTree->SetBranchAddress("Slowpi_TRUEID", &Slowpi_TRUEID);
   theTree->SetBranchAddress("Slowpi_MC_MOTHER_ID", &Slowpi_MC_MOTHER_ID);
   theTree->SetBranchAddress("Slowpi_TRUEP_E", &Slowpi_TRUEP_E);
   theTree->SetBranchAddress("Slowpi_TRUEP_X", &Slowpi_TRUEP_X);
   theTree->SetBranchAddress("Slowpi_TRUEP_Y", &Slowpi_TRUEP_Y);
   theTree->SetBranchAddress("Slowpi_TRUEP_Z", &Slowpi_TRUEP_Z);
   theTree->SetBranchAddress("Slowpi_TRUEPT", &Slowpi_TRUEPT);
   theTree->SetBranchAddress("Slowpi_TRUEORIGINVERTEX_X", &Slowpi_TRUEORIGINVERTEX_X);
   theTree->SetBranchAddress("Slowpi_TRUEORIGINVERTEX_Y", &Slowpi_TRUEORIGINVERTEX_Y);
   theTree->SetBranchAddress("Slowpi_TRUEORIGINVERTEX_Z", &Slowpi_TRUEORIGINVERTEX_Z);
   theTree->SetBranchAddress("Slowpi_TRUEISSTABLE", &Slowpi_TRUEISSTABLE);
   theTree->SetBranchAddress("Slowpi_TRUETAU", &Slowpi_TRUETAU);  
   theTree->SetBranchAddress("nTracks", &user_nTracks_data);
  }


  TFile *target  = new TFile(targetName,"RECREATE" );
  TTree* Tree = new TTree("BDT_Tree","BDT_Tree");         // new tree with BDT variable and all relevant other variables                                                        

  Tree->Branch("D_L0HadronDecision_TIS",&D_L0HadronDecision_TIS);
  Tree->Branch("D_L0MuonDecision_TIS",&D_L0MuonDecision_TIS);
  Tree->Branch("D_L0DiMuonDecision_TIS",&D_L0DiMuonDecision_TIS);
  Tree->Branch("D_L0PhotonDecision_TIS",&D_L0PhotonDecision_TIS);
  Tree->Branch("D_L0ElectronDecision_TIS",&D_L0ElectronDecision_TIS);
  Tree->Branch("D_L0HadronDecision_TOS",&D_L0HadronDecision_TOS);
  Tree->Branch("D_L0MuonDecision_TOS",&D_L0MuonDecision_TOS);
  Tree->Branch("D_L0DiMuonDecision_TOS",&D_L0DiMuonDecision_TOS);
  Tree->Branch("D_L0PhotonDecision_TOS",&D_L0PhotonDecision_TOS);
  Tree->Branch("D_L0ElectronDecision_TOS",&D_L0ElectronDecision_TOS);
 
  Tree->Branch("mu0_L0HadronDecision_TIS",&mu0_L0HadronDecision_TIS);
  Tree->Branch("mu0_L0MuonDecision_TIS",&mu0_L0MuonDecision_TIS);
  Tree->Branch("mu0_L0DiMuonDecision_TIS",&mu0_L0DiMuonDecision_TIS);
  Tree->Branch("mu0_L0PhotonDecision_TIS",&mu0_L0PhotonDecision_TIS);
  Tree->Branch("mu0_L0ElectronDecision_TIS",&mu0_L0ElectronDecision_TIS);
  Tree->Branch("mu0_L0HadronDecision_TOS",&mu0_L0HadronDecision_TOS);
  Tree->Branch("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Tree->Branch("mu0_L0DiMuonDecision_TOS",&mu0_L0DiMuonDecision_TOS);
  Tree->Branch("mu0_L0PhotonDecision_TOS",&mu0_L0PhotonDecision_TOS);
  Tree->Branch("mu0_L0ElectronDecision_TOS",&mu0_L0ElectronDecision_TOS);

  Tree->Branch("mu1_L0HadronDecision_TIS",&mu1_L0HadronDecision_TIS);
  Tree->Branch("mu1_L0MuonDecision_TIS",&mu1_L0MuonDecision_TIS);
  Tree->Branch("mu1_L0DiMuonDecision_TIS",&mu1_L0DiMuonDecision_TIS);
  Tree->Branch("mu1_L0PhotonDecision_TIS",&mu1_L0PhotonDecision_TIS);
  Tree->Branch("mu1_L0ElectronDecision_TIS",&mu1_L0ElectronDecision_TIS);
  Tree->Branch("mu1_L0HadronDecision_TOS",&mu1_L0HadronDecision_TOS);
  Tree->Branch("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Tree->Branch("mu1_L0DiMuonDecision_TOS",&mu1_L0DiMuonDecision_TOS);
  Tree->Branch("mu1_L0PhotonDecision_TOS",&mu1_L0PhotonDecision_TOS);
  Tree->Branch("mu1_L0ElectronDecision_TOS",&mu1_L0ElectronDecision_TOS);

  Tree->Branch("h0_L0HadronDecision_TIS",&mu0_L0HadronDecision_TIS);
  Tree->Branch("h0_L0MuonDecision_TIS",&mu0_L0MuonDecision_TIS);
  Tree->Branch("h0_L0DiMuonDecision_TIS",&mu0_L0DiMuonDecision_TIS);
  Tree->Branch("h0_L0PhotonDecision_TIS",&mu0_L0PhotonDecision_TIS);
  Tree->Branch("h0_L0ElectronDecision_TIS",&mu0_L0ElectronDecision_TIS);
  Tree->Branch("h0_L0HadronDecision_TOS",&mu0_L0HadronDecision_TOS);
  Tree->Branch("h0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Tree->Branch("h0_L0DiMuonDecision_TOS",&mu0_L0DiMuonDecision_TOS);
  Tree->Branch("h0_L0PhotonDecision_TOS",&mu0_L0PhotonDecision_TOS);
  Tree->Branch("h0_L0ElectronDecision_TOS",&mu0_L0ElectronDecision_TOS);

  Tree->Branch("h1_L0HadronDecision_TIS",&mu1_L0HadronDecision_TIS);
  Tree->Branch("h1_L0MuonDecision_TIS",&mu1_L0MuonDecision_TIS);
  Tree->Branch("h1_L0DiMuonDecision_TIS",&mu1_L0DiMuonDecision_TIS);
  Tree->Branch("h1_L0PhotonDecision_TIS",&mu1_L0PhotonDecision_TIS);
  Tree->Branch("h1_L0ElectronDecision_TIS",&mu1_L0ElectronDecision_TIS);
  Tree->Branch("h1_L0HadronDecision_TOS",&mu1_L0HadronDecision_TOS);
  Tree->Branch("h1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Tree->Branch("h1_L0DiMuonDecision_TOS",&mu1_L0DiMuonDecision_TOS);
  Tree->Branch("h1_L0PhotonDecision_TOS",&mu1_L0PhotonDecision_TOS);
  Tree->Branch("h1_L0ElectronDecision_TOS",&mu1_L0ElectronDecision_TOS);


  Tree->Branch("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Tree->Branch("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Tree->Branch("mu0_L0DiMuonDecision_TOS",&mu0_L0DiMuonDecision_TOS);
  Tree->Branch("mu1_L0DiMuonDecision_TOS",&mu1_L0DiMuonDecision_TOS);
  Tree->Branch("h1_L0MuonDecision_TOS",&h1_L0MuonDecision_TOS);
  Tree->Branch("h1_L0DiMuonDecision_TOS",&h1_L0DiMuonDecision_TOS);
  Tree->Branch("Dst_L0Global_TIS",&Dst_L0Global_TIS);
  Tree->Branch("D_L0Global_TIS",&D_L0Global_TIS);
  Tree->Branch("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  Tree->Branch("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  Tree->Branch("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);
  Tree->Branch("D_Hlt1DiMuonHighMassDecision_TOS",&D_Hlt1DiMuonHighMassDecision_TOS);
  Tree->Branch("D_Hlt1DiMuonLowMassDecision_TOS",&D_Hlt1DiMuonLowMassDecision_TOS);
  Tree->Branch("mu0_Hlt1SingleMuonNoIPDecision_TOS",&mu0_Hlt1SingleMuonNoIPDecision_TOS);
  Tree->Branch("mu1_Hlt1SingleMuonNoIPDecision_TOS",&mu1_Hlt1SingleMuonNoIPDecision_TOS);
  Tree->Branch("mu0_Hlt1SingleMuonHighPTDecision_TOS",&mu0_Hlt1SingleMuonHighPTDecision_TOS);
  Tree->Branch("mu1_Hlt1SingleMuonHighPTDecision_TOS",&mu1_Hlt1SingleMuonHighPTDecision_TOS);
  Tree->Branch("h1_Hlt1TrackMuonDecision_TOS",&h1_Hlt1TrackMuonDecision_TOS);
  Tree->Branch("D_Hlt2CharmSemilepD02KKMuMuDecision_TOS",&D_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
  Tree->Branch("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS",&D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS);
  Tree->Branch("D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS",&D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS",&Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS",&Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS);
  Tree->Branch("D_Hlt2DiMuonDetachedDecision_TOS",&D_Hlt2DiMuonDetachedDecision_TOS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS);
  Tree->Branch("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS",&Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_K3piDecision_TOS",&D_Hlt2CharmHadD02HHHH_K3piDecision_TOS);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS",&D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS);
  Tree->Branch("D_Hlt2CharmHadD02HHHH_4piDecision_TOS",&D_Hlt2CharmHadD02HHHH_4piDecision_TOS);

  Tree->Branch( "runNumber",&runNumber);
  Tree->Branch( "Polarity",&Polarity);
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
  Tree->Branch( "mu1_P", &user_mu1_P );
  Tree->Branch( "mu0_P", &user_mu0_P );
  Tree->Branch( "h0_P", &user_h0_P );
  Tree->Branch( "h1_P", &user_h1_P );
  Tree->Branch( "Slowpi_PT", &user_Slowpi_PT );
  Tree->Branch( "Slowpi_P", &user_Slowpi_P );
  Tree->Branch( "mu0_PX", &user_mu0_PX );
  Tree->Branch( "mu0_PY", &user_mu0_PY );
  Tree->Branch( "mu0_PZ", &user_mu0_PZ );
  Tree->Branch( "mu1_PX", &user_mu1_PX );
  Tree->Branch( "mu1_PY", &user_mu1_PY );
  Tree->Branch( "mu1_PZ", &user_mu1_PZ );
  Tree->Branch( "Slowpi_PX", &user_Slowpi_PX );
  Tree->Branch( "Slowpi_PY", &user_Slowpi_PY );
  Tree->Branch( "Slowpi_PZ", &user_Slowpi_PZ );
  Tree->Branch( "h0_PX", &user_h0_PX );
  Tree->Branch( "h0_PY", &user_h0_PY );
  Tree->Branch( "h0_PZ", &user_h0_PZ );
  Tree->Branch( "h1_PX", &user_h1_PX );
  Tree->Branch( "h1_PY", &user_h1_PY );
  Tree->Branch( "h1_PZ", &user_h1_PZ );
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
  Tree->Branch( "misID_dm_OS", &misID_dm_OS );
  Tree->Branch( "misID_mD_OS", &misID_mD_OS );
  Tree->Branch( "mHH", & user_mHH);
  Tree->Branch( "Dst_Coneptasy", & user_Dst_Coneptasy);
  Tree->Branch( "h0_ProbNNk", &user_h0_ProbNNk );
  Tree->Branch( "h1_ProbNNk", &user_h1_ProbNNk );
  Tree->Branch( "h0_ProbNNpi", &user_h0_ProbNNpi );
  Tree->Branch( "h1_ProbNNpi", &user_h1_ProbNNpi );
  Tree->Branch( "Slowpi_ProbNNpi", &user_Slowpi_ProbNNpi );
  Tree->Branch( "nTracks", &user_nTracks_data );
  Tree->Branch("mu1_MuonNShared",&mu1_MuonNShared);
  Tree->Branch("mu0_MuonNShared",&mu0_MuonNShared);
  Tree->Branch("mu0_isMuon",&mu0_isMuon);
  Tree->Branch("mu1_isMuon",&mu1_isMuon);
  Tree->Branch( "D_ENDVERTEX_CHI2",&user_D_ENDVERTEX_CHI2);
  Tree->Branch( "D_FDCHI2_OWNPV",&user_D_FDCHI2_OWNPV );
  Tree->Branch("h0_PIDK", &h0_PIDK);  
  Tree->Branch("h1_PIDK", &h1_PIDK);  
  Tree->Branch("nCandidate",&nCandidate);
  Tree->Branch("totCandidates",&totCandidates);

  Tree->Branch("D_Hlt1TrackAllL0Decision_TIS",&D_Hlt1TrackAllL0Decision_TIS);
  Tree->Branch("D_Hlt1L0AnyDecision_TIS",&D_Hlt1L0AnyDecision_TIS);
  Tree->Branch("D_Hlt1GlobalDecision_TIS",&D_Hlt1GlobalDecision_TIS);
  Tree->Branch("Dst_Hlt1Phys_TIS",&Dst_Hlt1Phys_TIS);
  Tree->Branch("h0_RichDLLp",&h0_RichDLLp);
  Tree->Branch("h1_RichDLLp",&h1_RichDLLp);
  Tree->Branch("h1_ProbNNp",&h1_ProbNNp);
  Tree->Branch("h0_ProbNNp",&h0_ProbNNp);
  Tree->Branch("h0_TRACK_GhostProb",&h0_TRACK_GhostProb);
  Tree->Branch("h1_TRACK_GhostProb",&h1_TRACK_GhostProb);
  Tree->Branch("mu0_TRACK_GhostProb",&mu0_TRACK_GhostProb);
  Tree->Branch("mu1_TRACK_GhostProb",&mu1_TRACK_GhostProb);
  Tree->Branch("Slowpi_TRACK_GhostProb",&Slowpi_TRACK_GhostProb);

  Tree->Branch("mKPiMuMu",&mKPiMuMu);
  Tree->Branch("mKMuPiMu",&mKMuPiMu);
  Tree->Branch("mKMuMuPi",&mKMuMuPi);

  Tree->Branch("mPiPiMuMu",&mPiPiMuMu);
  Tree->Branch("mPiMuPiMu",&mPiMuPiMu);
  Tree->Branch("mMuMuPiPi",&mMuMuPiPi);
  Tree->Branch("mMuPiMuPi",&mMuPiMuPi);
  Tree->Branch("mMuPiPiMu",&mMuPiPiMu);
  Tree->Branch("mPiMuMuPi",&mPiMuMuPi);
  Tree->Branch("mPiPiPiPiPi",&mPiPiPiPiPi);

  Tree->Branch("m34",&m34);
  Tree->Branch("m24",&m24);
  Tree->Branch("m23",&m23);
  Tree->Branch("m12",&m12);
  Tree->Branch("m13",&m13);
  Tree->Branch("m14",&m14);

  if(isMC) {

   Tree->Branch("Dst_BKGCAT", &Dst_BKGCAT);
   Tree->Branch("Dst_TRUEID", &Dst_TRUEID);
   Tree->Branch("Dst_MC_MOTHER_ID", &Dst_MC_MOTHER_ID);
   Tree->Branch("Dst_MC_MOTHER_KEY", &Dst_MC_MOTHER_KEY );
   Tree->Branch("Dst_TRUEP_E", &Dst_TRUEP_E);
   Tree->Branch("Dst_TRUEP_X", &Dst_TRUEP_X);
   Tree->Branch("Dst_TRUEP_Y", &Dst_TRUEP_Y);
   Tree->Branch("Dst_TRUEP_Z", &Dst_TRUEP_Z);
   Tree->Branch("Dst_TRUEPT", &Dst_TRUEPT);
   Tree->Branch("Dst_TRUEORIGINVERTEX_X", &Dst_TRUEORIGINVERTEX_X);
   Tree->Branch("Dst_TRUEORIGINVERTEX_Y", &Dst_TRUEORIGINVERTEX_Y);
   Tree->Branch("Dst_TRUEORIGINVERTEX_Z", &Dst_TRUEORIGINVERTEX_Z);
   Tree->Branch("Dst_TRUEENDVERTEX_X", &Dst_TRUEENDVERTEX_X);
   Tree->Branch("Dst_TRUEENDVERTEX_Y", &Dst_TRUEENDVERTEX_Y);
   Tree->Branch("Dst_TRUEENDVERTEX_Z", &Dst_TRUEENDVERTEX_Z);
   Tree->Branch("Dst_TRUEISSTABLE", &Dst_TRUEISSTABLE);
   Tree->Branch("Dst_TRUETAU", &Dst_TRUETAU);
   Tree->Branch("D_BKGCAT", &D_BKGCAT);
   Tree->Branch("D_TRUEID", &D_TRUEID);
   Tree->Branch("D_MC_MOTHER_ID", &D_MC_MOTHER_ID);
   Tree->Branch("D_TRUEP_E", &D_TRUEP_E);
   Tree->Branch("D_TRUEP_X", &D_TRUEP_X);
   Tree->Branch("D_TRUEP_Y", &D_TRUEP_Y);
   Tree->Branch("D_TRUEP_Z", &D_TRUEP_Z);
   Tree->Branch("D_TRUE_DiMuon_Mass", &D_TRUE_DiMuon_Mass);
   Tree->Branch("D_TRUEORIGINVERTEX_X", &D_TRUEORIGINVERTEX_X);
   Tree->Branch("D_TRUEORIGINVERTEX_Y", &D_TRUEORIGINVERTEX_Y);
   Tree->Branch("D_TRUEORIGINVERTEX_Z", &D_TRUEORIGINVERTEX_Z);
   Tree->Branch("D_TRUEENDVERTEX_X", &D_TRUEENDVERTEX_X);
   Tree->Branch("D_TRUEENDVERTEX_Y", &D_TRUEENDVERTEX_Y);
   Tree->Branch("D_TRUEENDVERTEX_Z", &D_TRUEENDVERTEX_Z);
   Tree->Branch("D_TRUEISSTABLE", &D_TRUEISSTABLE);
   Tree->Branch("D_TRUETAU", &D_TRUETAU);
   Tree->Branch("h0_TRUEID", &h0_TRUEID);
   Tree->Branch("h0_TRUEP_E", &h0_TRUEP_E);
   Tree->Branch("h0_TRUEP_X", &h0_TRUEP_X);
   Tree->Branch("h0_TRUEP_Y", &h0_TRUEP_Y);
   Tree->Branch("h0_TRUEP_Z", &h0_TRUEP_Z);
   Tree->Branch("h0_TRUEPT", &h0_TRUEPT);
   Tree->Branch("h0_TRUEORIGINVERTEX_X", &h0_TRUEORIGINVERTEX_X);
   Tree->Branch("h0_TRUEORIGINVERTEX_Y", &h0_TRUEORIGINVERTEX_Y);
   Tree->Branch("h0_TRUEORIGINVERTEX_Z", &h0_TRUEORIGINVERTEX_Z);
   Tree->Branch("h0_TRUEISSTABLE", &h0_TRUEISSTABLE);
   Tree->Branch("h0_TRUETAU", &h0_TRUETAU);
   Tree->Branch("h1_TRUEID", &h1_TRUEID);
   Tree->Branch("h1_TRUEP_E", &h1_TRUEP_E);
   Tree->Branch("h1_TRUEP_X", &h1_TRUEP_X);
   Tree->Branch("h1_TRUEP_Y", &h1_TRUEP_Y);
   Tree->Branch("h1_TRUEP_Z", &h1_TRUEP_Z);
   Tree->Branch("h1_TRUEPT", &h1_TRUEPT);
   Tree->Branch("h1_TRUEORIGINVERTEX_X", &h1_TRUEORIGINVERTEX_X);
   Tree->Branch("h1_TRUEORIGINVERTEX_Y", &h1_TRUEORIGINVERTEX_Y);
   Tree->Branch("h1_TRUEORIGINVERTEX_Z", &h1_TRUEORIGINVERTEX_Z);
   Tree->Branch("h1_TRUETAU", &h1_TRUETAU);
   Tree->Branch("mu0_TRUEID", &mu0_TRUEID);
   Tree->Branch("mu0_TRUEP_E", &mu0_TRUEP_E);
   Tree->Branch("mu0_TRUEP_X", &mu0_TRUEP_X);
   Tree->Branch("mu0_TRUEP_Y", &mu0_TRUEP_Y);
   Tree->Branch("mu0_TRUEP_Z", &mu0_TRUEP_Z);
   Tree->Branch("mu0_TRUEPT", &mu0_TRUEPT);
   Tree->Branch("mu0_TRUEORIGINVERTEX_X", &mu0_TRUEORIGINVERTEX_X);
   Tree->Branch("mu0_TRUEORIGINVERTEX_Y", &mu0_TRUEORIGINVERTEX_Y);
   Tree->Branch("mu0_TRUEORIGINVERTEX_Z", &mu0_TRUEORIGINVERTEX_Z);
   Tree->Branch("mu0_TRUETAU", &mu0_TRUETAU);
   Tree->Branch("mu1_TRUEID", &mu1_TRUEID );
   Tree->Branch("mu1_MC_MOTHER_ID", &mu1_MC_MOTHER_ID );
   Tree->Branch("mu1_TRUEP_E", &mu1_TRUEP_E);
   Tree->Branch("mu1_TRUEP_X", &mu1_TRUEP_X);
   Tree->Branch("mu1_TRUEP_Y", &mu1_TRUEP_Y);
   Tree->Branch("mu1_TRUEP_Z", &mu1_TRUEP_Z);
   Tree->Branch("mu1_TRUEPT", &mu1_TRUEPT);
   Tree->Branch("mu1_TRUEORIGINVERTEX_X", &mu1_TRUEORIGINVERTEX_X);
   Tree->Branch("mu1_TRUEORIGINVERTEX_Y", &mu1_TRUEORIGINVERTEX_Y);
   Tree->Branch("mu1_TRUEORIGINVERTEX_Z", &mu1_TRUEORIGINVERTEX_Z);
   Tree->Branch("mu1_TRUEISSTABLE", &mu1_TRUEISSTABLE);
   Tree->Branch("mu1_TRUETAU", &mu1_TRUETAU);
   Tree->Branch("Slowpi_TRUEID", &Slowpi_TRUEID);
   Tree->Branch("Slowpi_MC_MOTHER_ID", &Slowpi_MC_MOTHER_ID);
   Tree->Branch("Slowpi_TRUEP_E", &Slowpi_TRUEP_E);
   Tree->Branch("Slowpi_TRUEP_X", &Slowpi_TRUEP_X);
   Tree->Branch("Slowpi_TRUEP_Y", &Slowpi_TRUEP_Y);
   Tree->Branch("Slowpi_TRUEP_Z", &Slowpi_TRUEP_Z);
   Tree->Branch("Slowpi_TRUEPT", &Slowpi_TRUEPT);
   Tree->Branch("Slowpi_TRUEORIGINVERTEX_X", &Slowpi_TRUEORIGINVERTEX_X);
   Tree->Branch("Slowpi_TRUEORIGINVERTEX_Y", &Slowpi_TRUEORIGINVERTEX_Y);
   Tree->Branch("Slowpi_TRUEORIGINVERTEX_Z", &Slowpi_TRUEORIGINVERTEX_Z);
   Tree->Branch("Slowpi_TRUEISSTABLE", &Slowpi_TRUEISSTABLE);
   Tree->Branch("Slowpi_TRUETAU", &Slowpi_TRUETAU);  

   Tree->Branch("mu0_TRACK_Eta", &user_mu0_TRACK_Eta);
   Tree->Branch("mu1_TRACK_Eta", &user_mu1_TRACK_Eta);
   Tree->Branch("h1_TRACK_Eta", &user_h1_TRACK_Eta);
   Tree->Branch("h0_TRACK_Eta", &user_h0_TRACK_Eta);
   Tree->Branch("Slowpi_TRACK_Eta", &user_Slowpi_TRACK_Eta);

  }

  TLorentzVector mu0_hypMu;
  TLorentzVector mu1_hypMu;
  TLorentzVector h0_hypK,h0_hypPi,h0_hypMu;
  TLorentzVector h1_hypK,h1_hypPi,h1_hypMu;
  TLorentzVector mu0_hypK,mu0_hypPi;
  TLorentzVector mu1_hypK,mu1_hypPi;
  TLorentzVector Slowpi_hypPi;


  std::cout << "number Events: " << theTree->GetEntries() << std::endl;

  for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) { 

    if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
    //if (ievt>60000) break;                                                                                                                                                    
    theTree->GetEntry(ievt);
    ///////////SOME CUTS REMOVED RIGHT NOW!!!!!!!////////////////
    // IF NEEDED CUTS ARE ALREADY APLIED IN CREATESUBSAMPLE
    ////////////////////////////////////////////////////////////

    //if(user_Slowpi_ProbNNghost>0.5 || user_h1_ProbNNghost>0.5 ||
       //   user_h0_ProbNNghost>0.5 || user_mu1_ProbNNghost>0.5 ||  user_mu0_ProbNNghost > 0.5) continue;
    
    //cut on nMuonShared, included here as muon n shared cut is not included for hadron 4 body decays                                                                                        
    
    if( mu0_MuonNShared!=0 || mu1_MuonNShared!=0 ) continue;
    
    if(skipPID) goto skipped; //..sorry for the goto, but skip hadron PID if needed

    //// HADRON PID, different for signal and normalization modes /////
    if(isNormalizationMode && ( user_h0_ProbNNk<0.2 || user_h1_ProbNNpi<0.2) ) continue;
    if(!isNormalizationMode && ( user_h0_ProbNNpi<0.2 || user_h1_ProbNNpi<0.2)) continue;
    
  skipped:
    
    //compute true Di_DiuonMass
    mu0_hypMu.SetXYZM(mu0_TRUEP_X,mu0_TRUEP_Y,mu0_TRUEP_Z,Mass::Mu());
    mu1_hypMu.SetXYZM(mu1_TRUEP_X,mu1_TRUEP_Y,mu1_TRUEP_Z,Mass::Mu());
    D_TRUE_DiMuon_Mass=(mu0_hypMu+mu1_hypMu).M();
    
    //recompute Di_DiuonMass (needed especially for misID)
    mu0_hypMu.SetXYZM(DTF_mu0_PX,DTF_mu0_PY,DTF_mu0_PZ,Mass::Mu());
    mu1_hypMu.SetXYZM(DTF_mu1_PX,DTF_mu1_PY,DTF_mu1_PZ,Mass::Mu());
    user_D_DiMuon_Mass = (mu0_hypMu+mu1_hypMu).M();

    h1_hypK.SetXYZM(DTF_h1_PX,DTF_h1_PY,DTF_h1_PZ,Mass::K());
    h1_hypPi.SetXYZM(DTF_h1_PX,DTF_h1_PY,DTF_h1_PZ,Mass::Pi());
    h1_hypMu.SetXYZM(DTF_h1_PX,DTF_h1_PY,DTF_h1_PZ,Mass::Mu());
    h0_hypMu.SetXYZM(DTF_h0_PX,DTF_h0_PY,DTF_h0_PZ,Mass::Mu());
    h0_hypPi.SetXYZM(DTF_h0_PX,DTF_h0_PY,DTF_h0_PZ,Mass::Pi());
    h0_hypK.SetXYZM(DTF_h0_PX,DTF_h0_PY,DTF_h0_PZ,Mass::K());
    mu0_hypPi.SetXYZM(DTF_mu0_PX,DTF_mu0_PY,DTF_mu0_PZ,Mass::Pi());
    mu1_hypPi.SetXYZM(DTF_mu1_PX,DTF_mu1_PY,DTF_mu1_PZ,Mass::Pi());
    Slowpi_hypPi.SetXYZM(DTF_Pis_PX,DTF_Pis_PY,DTF_Pis_PZ,Mass::Pi());

    mKPiMuMu = ( h0_hypK + h1_hypPi + mu0_hypMu + mu1_hypMu).M();
    mKMuPiMu = ( h0_hypK + h1_hypMu + mu0_hypPi + mu1_hypMu).M();
    mKMuMuPi = ( h0_hypK + h1_hypMu + mu0_hypMu + mu1_hypPi).M();

    mPiPiMuMu=( h0_hypPi + h1_hypPi + mu0_hypMu + mu1_hypMu).M();
    mMuMuPiPi=( h0_hypMu + h1_hypMu + mu0_hypPi + mu1_hypPi).M();
    mMuPiMuPi=( h0_hypMu + h1_hypPi + mu0_hypMu + mu1_hypPi).M();
    mMuPiPiMu=( h0_hypMu + h1_hypPi + mu0_hypPi + mu1_hypMu).M();
    mPiMuMuPi=( h0_hypPi + h1_hypMu + mu0_hypMu + mu1_hypPi).M();
    mPiMuPiMu=( h0_hypPi + h1_hypMu + mu0_hypPi + mu1_hypMu).M();
    mPiPiPiPiPi=( h0_hypPi + h1_hypPi + mu0_hypPi + mu1_hypPi + Slowpi_hypPi).M();

    m34 = (mu0_hypMu + mu1_hypMu).M();
    m24 = (h1_hypMu + mu1_hypMu).M();
    m23 = (h1_hypMu + mu0_hypMu).M();
    m12 = (h0_hypMu + h1_hypMu).M();
    m13 = (h0_hypMu + mu0_hypMu).M();
    m14 = (h0_hypMu + mu1_hypMu).M();

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
    log_D_DIRA_OWNPV=float( TMath::Log(user_D_DIRA_OWNPV));
    min_mu_PT = float(TMath::Min(user_mu1_PT,user_mu0_PT));
    min_h_PT = float(TMath::Min(user_h1_PT,user_h0_PT));
    max_mu_PT = float(TMath::Max(user_mu1_PT,user_mu0_PT));
    Slowpi_P = float(user_Slowpi_P);
    Slowpi_PT = float(user_Slowpi_PT);
    D_MINIP= float(user_D_MINIP);
    D_MINIPCHI2= float(user_D_MINIPCHI2);
    D_ENDVERTEX_CHI2 = float(user_D_ENDVERTEX_CHI2);

    user_BDT= reader->EvaluateMVA("BDTG method");
    //if(BDT<0) continue;
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



  //signal data                                                                                 
  //Application takes name of tree in input file, name of the input file, the target name and part 1 or 2 for even or odd trained events, respectively. 
  

  /////////////////////////////                                                                                                                            
  //     Signal  Mode        //                                                                                                                        

  //////////////////////////// 
  
  /*
  
  Application_D2KKmumu("data/DecayTree_odd","D2KKmumu_PreselectedSubsample.root","D2KKmumu_BDT_odd.root",1);
  Application_D2KKmumu("data/DecayTree_even","D2KKmumu_PreselectedSubsample.root","D2KKmumu_BDT_even.root",2);

  TChain* myChain1= new TChain("BDT_Tree");
  myChain1->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT_odd.root");
  myChain1->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT_even.root");
  myChain1->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKmumu_BDT.root");
  
    
  //signal sideband                                                                                                                                                                          
  Application_D2KKmumu("sideband/DecayTree_odd","D2KKmumu_PreselectedSubsample.root","sideband_D2KKmumu_BDT_odd.root",1);
  Application_D2KKmumu("sideband/DecayTree_even","D2KKmumu_PreselectedSubsample.root","sideband_D2KKmumu_BDT_even.root",2);
  TChain* myChain2= new TChain("BDT_Tree");
  myChain2->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT_odd.root");
  myChain2->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT_even.root");
  myChain2->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2KKmumu_BDT.root");
  
  
  //signal MC                                                                                                                                                                                
  
  Application_D2KKmumu("DecayTree_odd","D2KKmumu_MCtrainingSample.root","MC_D2KKmumu_BDT_odd.root",1,true);
  Application_D2KKmumu("DecayTree_even","D2KKmumu_MCtrainingSample.root","MC_D2KKmumu_BDT_even.root",2,true);

  TChain* myChain3= new TChain("BDT_Tree");
  myChain3->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT_odd.root");
  myChain3->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT_even.root");
  myChain3->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKmumu_BDT.root");
 
  */ 
    
  
  /////////////////////////////

  // Normalization Mode      //
  
  ////////////////////////////
  
  /*
  //normalization mode data
  Application_D2KKmumu("data/DecayTree_odd","D2Kpimumu_PreselectedSubsample.root","D2Kpimumu_D2KKmumuBDT_odd.root",1,false,true);
  Application_D2KKmumu("data/DecayTree_even","D2Kpimumu_PreselectedSubsample.root","D2Kpimumu_D2KKmumuBDT_even.root",2,false,true);

  TChain* myChain4= new TChain("BDT_Tree");
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT_odd.root");
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT_even.root");
  myChain4->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2KKmumuBDT.root");
  

  */
    
  /*
  //normalization mode sideband                                                                                                                                                                           
  Application_D2KKmumu("sideband/DecayTree_odd","D2Kpimumu_PreselectedSubsample.root","sideband_D2Kpimumu_D2KKmumuBDT_odd.root",1);
  Application_D2KKmumu("sideband/DecayTree_even","D2Kpimumu_PreselectedSubsample.root","sideband_D2Kpimumu_D2KKmumuBDT_even.root",2);
  TChain* myChain7= new TChain("BDT_Tree");
  myChain7->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2Kpimumu_D2KKmumuBDT_odd.root");
  myChain7->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2Kpimumu_D2KKmumuBDT_even.root");
  myChain7->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2Kpimumu_D2KKmumuBDT.root");

  
  
  //nomrmalization MC                                                                                                                                          
                                                                                                                                                          
  Application_D2KKmumu("DecayTree_odd","D2Kpimumu_MCtrainingSample.root","MC_D2Kpimumu_D2KKmumuBDT_odd.root",1,true,true);
  Application_D2KKmumu("DecayTree_even","D2Kpimumu_MCtrainingSample.root","MC_D2Kpimumu_D2KKmumuBDT_even.root",2,true,true);

  TChain* myChain10= new TChain("BDT_Tree");
  myChain10->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT_odd.root");
  myChain10->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT_even.root");
  myChain10->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root");
  */ 


  
  /////////////////////////////                                                                                                                            

  // misID D2hhhh  Mode      //                                                                                                                           

  //////////////////////////// 
  
    
  //mis ID bkg Data                                                                                                                                                            
  Application_D2KKmumu("data/DecayTree_odd","D2Kpipipi_PreselectedRandomizedSubsample.root","D2Kpipipi_D2KKmumuBDT_odd.root",1,false,true,false); 
  Application_D2KKmumu("data/DecayTree_even","D2Kpipipi_PreselectedRandomizedSubsample.root","D2Kpipipi_D2KKmumuBDT_even.root",2,false,true,false);                                 
  TChain* myChain8= new TChain("BDT_Tree");                                                                                                              
  myChain8->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_odd.root");        
  myChain8->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_even.root");                                                    
  myChain8->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2KKmumuBDT_new.root");    
  
  
  //Application_D2KKmumu("data/DecayTree_odd","D2KKpipi_PreselectedSubsample.root","D2KKpipi_D2KKmumuBDT_odd.root",1,false,false,false);
  //Application_D2KKmumu("data/DecayTree_even","D2KKpipi_PreselectedSubsample.root","D2KKpipi_D2KKmumuBDT_even.root",2,false,false,false);
  //TChain* myChain11= new TChain("BDT_Tree");
  //myChain11->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT_odd.root");
  //myChain11->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT_even.root");
  //myChain11->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2KKpipi_D2KKmumuBDT_new.root");
    
  
  /*
  //mis ID bkg MC 

  Application_D2KKmumu("DecayTree_odd","D2KKpipi_MCtrainingSample.root","MC_D2KKpipi_D2KKmumuBDT_odd.root",1,true,false,false);
  Application_D2KKmumu("DecayTree_even","D2KKpipi_MCtrainingSample.root","MC_D2KKpipi_D2KKmumuBDT_even.root",2,true,false,false);
  TChain* myChain5= new TChain("BDT_Tree");
  myChain5->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKpipi_D2KKmumuBDT_odd.root");
  myChain5->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKpipi_D2KKmumuBDT_even.root");
  myChain5->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2KKpipi_D2KKmumuBDT.root");

  TChain* myChain6= new TChain("BDT_Tree");  
  Application_D2KKmumu("DecayTree_odd","D2Kpipipi_MCtrainingSample.root","MC_D2Kpipipi_D2KKmumuBDT_odd.root",1,true,true,false); 
  Application_D2KKmumu("DecayTree_even","D2Kpipipi_MCtrainingSample.root","MC_D2Kpipipi_D2KKmumuBDT_even.root",2,true,true,false);  
  myChain6->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2KKmumuBDT_odd.root");
  myChain6->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2KKmumuBDT_even.root");
  myChain6->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2KKmumuBDT.root");  
  */
 
  std::cout << "==> TMVAClassificationApplication is done!"  << std::endl;
  
  }


void D2pipimumuCrossapplication(){



  //signal data                                                                                 
  //Application takes name of tree in input file, name of the input file, the target name and part 1 or 2 for even or odd trained events, respectively. 
  

  /////////////////////////////                                                                                                                            
  //     Signal  Mode        //                                                                                                                        

  //////////////////////////// 
  
  /*  
    
  Application_D2pipimumu("data/DecayTree_odd","D2pipimumu_PreselectedSubsample.root","D2pipimumu_BDT_odd.root",1);
  Application_D2pipimumu("data/DecayTree_even","D2pipimumu_PreselectedSubsample.root","D2pipimumu_BDT_even.root",2);

  TChain* myChain1= new TChain("BDT_Tree");
  myChain1->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT_odd.root");
  myChain1->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT_even.root");
  myChain1->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root");
  
  */
  /*
  ///NO CUTS                                                                                                                                                                                  
  Application_D2pipimumu("data/DecayTree_odd","D2pipimumu_PreselectedSubsample_noCuts.root","D2pipimumu_D2pipimumuBDT_noCuts_odd.root",1,false,true);
  Application_D2pipimumu("data/DecayTree_even","D2pipimumu_PreselectedSubsample_noCuts.root","D2pipimumu_D2pipimumuBDT_noCuts_even.root",2,false,true);

  TChain* myChain5= new TChain("BDT_Tree");
  myChain5->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_D2pipimumuBDT_noCuts_odd.root");
  myChain5->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_D2pipimumuBDT_noCuts_even.root");
  myChain5->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipimumu_D2pipimumuBDT_noCuts.root");
  */


  /*
  //signal sideband                                                                                                                                                                          
  Application_D2pipimumu("sideband/DecayTree_odd","D2pipimumu_PreselectedSubsample.root","sideband_D2pipimumu_BDT_odd.root",1);
  Application_D2pipimumu("sideband/DecayTree_even","D2pipimumu_PreselectedSubsample.root","sideband_D2pipimumu_BDT_even.root",2);
  TChain* myChain2= new TChain("BDT_Tree");
  myChain2->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2pipimumu_BDT_odd.root");
  myChain2->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2pipimumu_BDT_even.root");
  myChain2->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2pipimumu_BDT.root");
  
  */
  //signal MC                                                                                                                                                                         
  /*
  Application_D2pipimumu("DecayTree_odd","D2pipimumu_MCtrainingSample.root","MC_D2pipimumu_BDT_odd.root",1,true);
  Application_D2pipimumu("DecayTree_even","D2pipimumu_MCtrainingSample.root","MC_D2pipimumu_BDT_even.root",2,true);

  TChain* myChain3= new TChain("BDT_Tree");
  myChain3->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT_odd.root");
  myChain3->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT_even.root");
  myChain3->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root");
  */

  /////////////////////////////
  // Normalization Mode      //
  ////////////////////////////
  
  /*
  //normalization mode data
  
  Application_D2pipimumu("data/DecayTree_odd","D2Kpimumu_PreselectedSubsample.root","D2Kpimumu_D2pipimumuBDT_odd.root",1,false,true);
  Application_D2pipimumu("data/DecayTree_even","D2Kpimumu_PreselectedSubsample.root","D2Kpimumu_D2pipimumuBDT_even.root",2,false,true);

  TChain* myChain4= new TChain("BDT_Tree");
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT_odd.root");
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT_even.root");
  myChain4->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root");
  */
  
  
  /*
  ///NO CUTS                                                                                                                                                                                  
  Application_D2pipimumu("data/DecayTree_odd","D2Kpimumu_PreselectedSubsample_noCuts.root","D2Kpimumu_D2pipimumuBDT_noCuts_odd.root",1,false,true);
  Application_D2pipimumu("data/DecayTree_even","D2Kpimumu_PreselectedSubsample_noCuts.root","D2Kpimumu_D2pipimumuBDT_noCuts_even.root",2,false,true);

  TChain* myChain4= new TChain("BDT_Tree");
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2pipimumuBDT_noCuts_odd.root");
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2pipimumuBDT_noCuts_even.root");
  myChain4->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpimumu_D2pipimumuBDT_noCuts.root");
  */
  /* 
  //normalization mode sideband                                                                                                                                                                           
  Application_D2KKmumu("sideband/DecayTree_odd","D2Kpimumu_PreselectedSubsample.root","sideband_D2Kpimumu_D2KKmumuBDT_odd.root",1);
  Application_D2KKmumu("sideband/DecayTree_even","D2Kpimumu_PreselectedSubsample.root","sideband_D2Kpimumu_D2KKmumuBDT_even.root",2);
  TChain* myChain7= new TChain("BDT_Tree");
  myChain7->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2Kpimumu_D2KKmumuBDT_odd.root");
  myChain7->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2Kpimumu_D2KKmumuBDT_even.root");
  myChain7->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/sideband_D2Kpimumu_D2KKmumuBDT.root");
    

  //nomrmalization MC                                                                                                                                          
                                                                                                                                                          
  Application_D2pipimumu("DecayTree_odd","D2Kpimumu_MCtrainingSample.root","MC_D2Kpimumu_D2pipimumuBDT_odd.root",1,true,true);
  Application_D2pipimumu("DecayTree_even","D2Kpimumu_MCtrainingSample.root","MC_D2Kpimumu_D2pipimumuBDT_even.root",2,true,true);

  TChain* myChain10= new TChain("BDT_Tree");
  myChain10->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT_odd.root");
  myChain10->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT_even.root");
  myChain10->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");
  
  */

  /////////////////////////////                                                                                                                            

  // misID D2hhhh  Mode      //                                                                                                                           

  //////////////////////////// 

  
  Application_D2pipimumu("data/DecayTree_odd","D2pipipipi_PreselectedRandomizedSubsample.root","D2pipipipi_D2pipimumuBDT_odd.root",1,false,false,false);
  Application_D2pipimumu("data/DecayTree_even","D2pipipipi_PreselectedRandomizedSubsample.root","D2pipipipi_D2pipimumuBDT_even.root",2,false,false,false);
  TChain* myChain11= new TChain("BDT_Tree");
  myChain11->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_odd.root");
  myChain11->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_even.root");
  myChain11->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT_new.root");

   
  //mis ID bkg Data                                                                                                                                                            
  Application_D2pipimumu("data/DecayTree_odd","D2Kpipipi_PreselectedRandomizedSubsample.root","D2Kpipipi_D2pipimumuBDT_odd.root",1,false,true,false); 
  Application_D2pipimumu("data/DecayTree_even","D2Kpipipi_PreselectedRandomizedSubsample.root","D2Kpipipi_D2pipimumuBDT_even.root",2,false,true,false);   
  TChain* myChain8= new TChain("BDT_Tree");                                                                                                              
  myChain8->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_odd.root");        
  myChain8->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_even.root");                                                    
  myChain8->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpipipi_D2pipimumuBDT_new.root");     
  
  
  /*    
  //mis ID bkg MC 

  Application_D2pipimumu("DecayTree_odd","D2pipipipi_MCtrainingSample.root","MC_D2pipipipi_D2pipimumuBDT_odd.root",1,true,false,false);
  Application_D2pipimumu("DecayTree_even","D2pipipipi_MCtrainingSample.root","MC_D2pipipipi_D2pipimumuBDT_even.root",2,true,false,false);
  TChain* myChain5= new TChain("BDT_Tree");
  myChain5->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipipipi_D2pipimumuBDT_odd.root");
  myChain5->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipipipi_D2pipimumuBDT_even.root");
  myChain5->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipipipi_D2pipimumuBDT.root");

  
  
  Application_D2pipimumu("DecayTree_odd","D2Kpipipi_MCtrainingSample.root","MC_D2Kpipipi_D2pipimumuBDT_odd.root",1,true,true,false); 
  Application_D2pipimumu("DecayTree_even","D2Kpipipi_MCtrainingSample.root","MC_D2Kpipipi_D2pipimumuBDT_even.root",2,true,true,false);  
  TChain* myChain6= new TChain("BDT_Tree");                                                                                                              
  myChain6->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2pipimumuBDT_odd.root");
  myChain6->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2pipimumuBDT_even.root");
  myChain6->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2pipimumuBDT.root");  
  
  */

  std::cout << "==> TMVAClassificationApplication is done!"  << std::endl;
  
  }



void CrossapplicationForEfficiencyStudies(){



  /////////////////////////////                                                                                                                            
  //     Signal  Mode        //                                                                                                                        

  //////////////////////////// 
  
  std::vector<double> rangesKpi_low = {675};
  std::vector<double> rangesKpi_high = {875};
  std::vector<double> rangespipi_low = {200,525,565,950,1100};
  std::vector<double> rangespipi_high = {525,565,950,1100,1600};
  std::vector<double> rangesKK_low = {200,525,565};
  std::vector<double> rangesKK_high = {525,565,900};
 
  //std::vector<double> rangesKpi_low = {675};
  //std::vector<double> rangesKpi_high = {875};
  //std::vector<double> rangespipi_low = {0};
  //std::vector<double> rangespipi_high = {1500};
  //std::vector<double> rangesKK_low = {0,525};                                                                                                                                        
  // std::vector<double> rangesKK_high = {525,565};//!!//950!!                                                                                                                         
  //std::vector<double> rangesKK_low = {0};
  //std::vector<double> rangesKK_high = {950};//!!//950!!

 

  //D2Kpimumu
  
  for(int i=0; i<rangesKpi_low.size();++i){
    //magUp
    Application_D2pipimumu("DecayTree_odd",TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_odd_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),1,true,true,true);
    Application_D2pipimumu("DecayTree_even",TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_event_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),2,true,true,true);

  TChain* myChain1= new TChain("BDT_Tree");
  myChain1->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_odd_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain1->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_event_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain1->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));

    //magDw
  Application_D2pipimumu("DecayTree_odd",TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_odd_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),1,true,true,true);
  Application_D2pipimumu("DecayTree_even",TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_event_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),2,true,true,true);

  TChain* myChain2= new TChain("BDT_Tree");
  myChain2->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_odd_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain2->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_event_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain2->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
 
  //D2KKmumu trained

  //magUp
  Application_D2KKmumu("DecayTree_odd",TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_odd_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),1,true,true,true);
  Application_D2KKmumu("DecayTree_even",TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_event_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),2,true,true,true);

  TChain* myChain3= new TChain("BDT_Tree");
  myChain3->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_odd_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain3->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_event_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain3->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));

    //magDw
   Application_D2KKmumu("DecayTree_odd",TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_odd_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),1,true,true,true);
   Application_D2KKmumu("DecayTree_even",TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_event_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),2,true,true,true);

  TChain* myChain4= new TChain("BDT_Tree");
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_odd_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_event_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain4->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));

 }
  
  
  //D2pipimumu

  for(int i=0; i<rangespipi_low.size();++i){
    //magUp
    Application_D2pipimumu("DecayTree_odd",TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_magUp.root",rangespipi_low[i],rangespipi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_odd_magUp.root",rangespipi_low[i],rangespipi_high[i]),1,true,true,true);
    Application_D2pipimumu("DecayTree_even",TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_magUp.root",rangespipi_low[i],rangespipi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_event_magUp.root",rangespipi_low[i],rangespipi_high[i]),2,true,true,true);

  TChain* myChain5= new TChain("BDT_Tree");
  myChain5->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_odd_magUp.root",rangespipi_low[i],rangespipi_high[i]));
  myChain5->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_event_magUp.root",rangespipi_low[i],rangespipi_high[i]));
  myChain5->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]));

    //magDw
  Application_D2pipimumu("DecayTree_odd",TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_magDw.root",rangespipi_low[i],rangespipi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_odd_magDw.root",rangespipi_low[i],rangespipi_high[i]),1,true,true,true);
  Application_D2pipimumu("DecayTree_even",TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_magDw.root",rangespipi_low[i],rangespipi_high[i]),TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_event_magDw.root",rangespipi_low[i],rangespipi_high[i]),2,true,true,true);

  TChain* myChain6= new TChain("BDT_Tree");
  myChain6->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_odd_magDw.root",rangespipi_low[i],rangespipi_high[i]));
  myChain6->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_event_magDw.root",rangespipi_low[i],rangespipi_high[i]));
  myChain6->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]));

 }
 

  //KKmumu
  for(int i=0; i<rangesKK_low.size();++i){
    //magUp
     Application_D2KKmumu("DecayTree_odd",TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_magUp.root",rangesKK_low[i],rangesKK_high[i]),TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_odd_magUp.root",rangesKK_low[i],rangesKK_high[i]),1,true,true,true);
    Application_D2KKmumu("DecayTree_even",TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_magUp.root",rangesKK_low[i],rangesKK_high[i]),TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_event_magUp.root",rangesKK_low[i],rangesKK_high[i]),2,true,true,true);

  TChain* myChain7= new TChain("BDT_Tree");
  myChain7->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_odd_magUp.root",rangesKK_low[i],rangesKK_high[i]));
  myChain7->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_event_magUp.root",rangesKK_low[i],rangesKK_high[i]));
  myChain7->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]));

    //magDw
   Application_D2KKmumu("DecayTree_odd",TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_magDw.root",rangesKK_low[i],rangesKK_high[i]),TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_odd_magDw.root",rangesKK_low[i],rangesKK_high[i]),1,true,true,true);
   Application_D2KKmumu("DecayTree_even",TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_magDw.root",rangesKK_low[i],rangesKK_high[i]),TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_event_magDw.root",rangesKK_low[i],rangesKK_high[i]),2,true,true,true);

  TChain* myChain8= new TChain("BDT_Tree");
  myChain8->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_odd_magDw.root",rangesKK_low[i],rangesKK_high[i]));
  myChain8->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_event_magDw.root",rangesKK_low[i],rangesKK_high[i]));
  myChain8->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]));

 }
  
  std::cout << "==> TMVAClassificationApplication is done!"  << std::endl;
  
}



void CrossapplicationForEfficiencyStudiesNoTruthmatching(){


  /////////////////////////////                                                                                                                            
  //     Signal  Mode        //                                                                                                                        

  //////////////////////////// 
  
  std::vector<double> rangesKpi_low = {200};
  std::vector<double> rangesKpi_high = {1600};
  std::vector<double> rangespipi_low = {200};
  std::vector<double> rangespipi_high = {1600};
  std::vector<double> rangesKK_low = {200};
  std::vector<double> rangesKK_high = {1600};
 
 

  //D2Kpimumu
  
  for(int i=0; i<rangesKpi_low.size();++i){
    //magUp
    Application_D2pipimumu("DecayTree_odd",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_odd_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),1,true,true,true);
    Application_D2pipimumu("DecayTree_even",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_event_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),2,true,true,true);

  TChain* myChain1= new TChain("BDT_Tree");
  myChain1->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_odd_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain1->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_event_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain1->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));

  

    //magDw
  Application_D2pipimumu("DecayTree_odd",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_odd_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),1,true,true,true);
  Application_D2pipimumu("DecayTree_even",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_event_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),2,true,true,true);

  TChain* myChain2= new TChain("BDT_Tree");
  myChain2->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_odd_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain2->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_event_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain2->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
 
  //D2KKmumu trained

  //magUp
  Application_D2KKmumu("DecayTree_odd",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_odd_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),1,true,true,true);
  Application_D2KKmumu("DecayTree_even",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_event_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),2,true,true,true);

  TChain* myChain3= new TChain("BDT_Tree");
  myChain3->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_odd_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain3->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_event_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain3->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));

    //magDw
   Application_D2KKmumu("DecayTree_odd",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_odd_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),1,true,true,true);
   Application_D2KKmumu("DecayTree_even",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_event_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),2,true,true,true);

  TChain* myChain4= new TChain("BDT_Tree");
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_odd_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain4->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_event_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
  myChain4->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));

 }
  
  
  //D2pipimumu

  for(int i=0; i<rangespipi_low.size();++i){
    //magUp
    Application_D2pipimumu("DecayTree_odd",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_magUp.root",rangespipi_low[i],rangespipi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_odd_magUp.root",rangespipi_low[i],rangespipi_high[i]),1,true,true,true);
    Application_D2pipimumu("DecayTree_even",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_magUp.root",rangespipi_low[i],rangespipi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_event_magUp.root",rangespipi_low[i],rangespipi_high[i]),2,true,true,true);

  TChain* myChain5= new TChain("BDT_Tree");
  myChain5->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_odd_magUp.root",rangespipi_low[i],rangespipi_high[i]));
  myChain5->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_event_magUp.root",rangespipi_low[i],rangespipi_high[i]));
  myChain5->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]));

    //magDw
  Application_D2pipimumu("DecayTree_odd",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_magDw.root",rangespipi_low[i],rangespipi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_odd_magDw.root",rangespipi_low[i],rangespipi_high[i]),1,true,true,true);
  Application_D2pipimumu("DecayTree_even",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_magDw.root",rangespipi_low[i],rangespipi_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_event_magDw.root",rangespipi_low[i],rangespipi_high[i]),2,true,true,true);

  TChain* myChain6= new TChain("BDT_Tree");
  myChain6->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_odd_magDw.root",rangespipi_low[i],rangespipi_high[i]));
  myChain6->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_event_magDw.root",rangespipi_low[i],rangespipi_high[i]));
  myChain6->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]));

 }
 

  //KKmumu
  for(int i=0; i<rangesKK_low.size();++i){
    //magUp
     Application_D2KKmumu("DecayTree_odd",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_magUp.root",rangesKK_low[i],rangesKK_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_odd_magUp.root",rangesKK_low[i],rangesKK_high[i]),1,true,true,true);
    Application_D2KKmumu("DecayTree_even",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_magUp.root",rangesKK_low[i],rangesKK_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_event_magUp.root",rangesKK_low[i],rangesKK_high[i]),2,true,true,true);

  TChain* myChain7= new TChain("BDT_Tree");
  myChain7->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_odd_magUp.root",rangesKK_low[i],rangesKK_high[i]));
  myChain7->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_event_magUp.root",rangesKK_low[i],rangesKK_high[i]));
  myChain7->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]));

    //magDw
   Application_D2KKmumu("DecayTree_odd",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_magDw.root",rangesKK_low[i],rangesKK_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_odd_magDw.root",rangesKK_low[i],rangesKK_high[i]),1,true,true,true);
   Application_D2KKmumu("DecayTree_even",TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_magDw.root",rangesKK_low[i],rangesKK_high[i]),TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_event_magDw.root",rangesKK_low[i],rangesKK_high[i]),2,true,true,true);

  TChain* myChain8= new TChain("BDT_Tree");
  myChain8->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_odd_magDw.root",rangesKK_low[i],rangesKK_high[i]));
  myChain8->Add("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_event_magDw.root",rangesKK_low[i],rangesKK_high[i]));
  myChain8->Merge("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/"+TString::Format("MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]));

  

 }
  
 std::cout << "==> TMVAClassificationApplication is done!"  << std::endl;
  
}


