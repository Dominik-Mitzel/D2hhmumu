#include "D2KpipipiReader.h"
#include <iostream>
#include "TGenPhaseSpace.h"

double D2KpipipiReader::get_mMuMu_noDCF()
{
  initializeMomenta();
  double mMuMu_noDCF = (pDTFMu1+pDTFMu0+pDTFH1+pDTFH0).M();
  
  return mMuMu_noDCF;
}


double D2KpipipiReader::get_mMuMu_doubleDCF()
{
  initializeMomenta();

  TVector3 b0=(pDTFPi0).BoostVector();
  TVector3 b1=(pDTFPi0).BoostVector();

  TLorentzVector pcm0(pDTFPi0);
  TLorentzVector pcm1(pDTFPi1);

  pcm0.Boost(-b0);
  pcm1.Boost(-b1);

  Double_t masses[2] = {Mass::Mu(), 0. };
  TGenPhaseSpace event1;
  event1.SetDecay(pcm0, 2, masses);
  event1.Generate();
  TLorentzVector* pPiDecayed0  = event1.GetDecay(0);
  TGenPhaseSpace event2;
  event2.SetDecay(pcm1, 2, masses);
  event2.Generate();
  TLorentzVector *pPiDecayed1  = event2.GetDecay(0);

  pPiDecayed0->Boost(b0);
  pPiDecayed1->Boost(b1);

  double mMuMu_doubleDCF = ( *pPiDecayed0+ *pPiDecayed1 +pDTFH1+pDTFH0).M();

  return mMuMu_doubleDCF;
}


double D2KpipipiReader::get_mMuMu_DCF_lowP()
{
  initializeMomenta();

  Double_t masses[2] = {Mass::Mu(), 0. };
  TGenPhaseSpace event;
  
  if(pDTFPi0.P() < pDTFPi1.P() ) {
    event.SetDecay(pDTFPi0, 2, masses);
  }
  else {
    event.SetDecay(pDTFPi1, 2, masses);
  }

  event.Generate();
  TLorentzVector* pPiDecayed  = event.GetDecay(0);
  double mMuMu_DCF_lowP;

  if(pDTFPi0.P() < pDTFPi1.P() ) {
    mMuMu_DCF_lowP=(*pPiDecayed+pDTFMu1+pDTFH1+pDTFH0).M();
  }
  else {
    mMuMu_DCF_lowP=(*pPiDecayed+pDTFMu0+pDTFH1+pDTFH0).M();
  }  

  return mMuMu_DCF_lowP;
}

double D2KpipipiReader::get_mMuMu_DCF_highP()
{
  initializeMomenta();

  Double_t masses[2] = {Mass::Mu(), 0. };
  TGenPhaseSpace event;

  if(pDTFPi0.P() > pDTFPi1.P() ) {
    event.SetDecay(pDTFPi0, 2, masses);
  }
  else {
    event.SetDecay(pDTFPi1, 2, masses);
  }

  event.Generate();
  TLorentzVector *pPiDecayed  = event.GetDecay(0);
  double mMuMu_DCF_highP;

  if(pDTFPi0.P() > pDTFPi1.P() ) {
    mMuMu_DCF_highP=(*pPiDecayed+pDTFMu1+pDTFH0+pDTFH1).M();
  }
  else {
    mMuMu_DCF_highP=(*pPiDecayed+pDTFMu0+pDTFH1+pDTFH0).M();
  }

  return mMuMu_DCF_highP;
}



    

D2KpipipiReader::D2KpipipiReader(TTree *tree) 
{

  Init(tree);
}


D2KpipipiReader::~D2KpipipiReader()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}


void D2KpipipiReader::initializeMomenta(){                                                                                                                                                  


  pH1.SetXYZM(h1_PX,h1_PY,h1_PZ,Mass::Pi());
  pH0.SetXYZM(h0_PX,h0_PY,h0_PZ,Mass::K());
  pMu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
  pMu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());

  pPi0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Pi());
  pPi1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Pi());

  pDTFH1.SetXYZM(Dst_DTF_h1_PX,Dst_DTF_h1_PY,Dst_DTF_h1_PZ,Mass::Pi());
  pDTFH0.SetXYZM(Dst_DTF_h0_PX,Dst_DTF_h0_PY,Dst_DTF_h0_PZ,Mass::K());
  pDTFMu1.SetXYZM(Dst_DTF_mu1_PX,Dst_DTF_mu1_PY,Dst_DTF_mu1_PZ,Mass::Mu());
  pDTFMu0.SetXYZM(Dst_DTF_mu0_PX,Dst_DTF_mu0_PY,Dst_DTF_mu0_PZ,Mass::Mu());

  pDTFPi1.SetXYZM(Dst_DTF_mu1_PX,Dst_DTF_mu1_PY,Dst_DTF_mu1_PZ,Mass::Pi());
  pDTFPi0.SetXYZM(Dst_DTF_mu0_PX,Dst_DTF_mu0_PY,Dst_DTF_mu0_PZ,Mass::Pi());

  pD.SetXYZM(D_PX,D_PY,D_PZ,Mass::D0());
  pDst.SetXYZM(Dst_PX,Dst_PY,Dst_PZ,Mass::Ds());
  pPis.SetXYZM(Slowpi_PX,Slowpi_PY,Slowpi_PZ,Mass::Pi());

  pDTFDst.SetXYZM(Dst_DTF_Dstarplus_PX,Dst_DTF_Dstarplus_PY,Dst_DTF_Dstarplus_PZ,Mass::Ds());
  pDTFD.SetXYZM(Dst_DTF_D0_PX,Dst_DTF_D0_PY,Dst_DTF_D0_PZ,Mass::D0());
  pDTFPis.SetXYZM(Dst_DTF_Pis_PX,Dst_DTF_Pis_PY,Dst_DTF_Pis_PZ,Mass::Pi());



}


bool D2KpipipiReader::isHlt2Selected(){

  // if(D_Hlt2CharmSemilepD02KKMuMuDecision_TOS==1) return true;
  //else return false;
  // at the moment, only HLT2 filtered MC of the D2hhhh modes is available, therefore do not cut on HLT2 
  return true;
}



void D2KpipipiReader::addMisIdMasses(TString name) {

  InitMC();
  //activateRelevantBranches();                                                                                                                    

  Long64_t nentries = fChain->GetEntries();
  std::cout<<"Found tree with "<<nentries <<" entries..."<<std::endl;

  fChain->GetEntry(0);
  double mMuMu_doubleDCF, mMuMu_noDCF, mMuMu_DCF_highP, mMuMu_DCF_lowP;
  double Slowpi_cosh,mu0_cosh,D_cosh, deltaM;
  double D_Conemult,Dst_Conemult,D_Coneptasy,Dst_Coneptasy;
  double mHH;

  //Create a new file + a clone of old tree in new file                                                                                            
  TFile *newfile = new TFile(name,"recreate");
  TTree *newtree = fChain->CloneTree(0);
  
  newtree->Branch("mMuMu_noDCF", &mMuMu_noDCF );
  newtree->Branch("mMuMu_doubleDCF", &mMuMu_doubleDCF );
  newtree->Branch("mMuMu_DCF_highP", &mMuMu_DCF_highP );
  newtree->Branch("mMuMu_DCF_lowP", &mMuMu_DCF_lowP );
  newtree->Branch("Slowpi_cosh",&Slowpi_cosh);
  newtree->Branch("D_cosh", & D_cosh);
  newtree->Branch("mu0_cosh",&mu0_cosh);
  newtree->Branch("deltaM",&deltaM);
  newtree->Branch("D_Conemult",&D_Conemult);
  newtree->Branch("Dst_Conemult",&Dst_Conemult);
  newtree->Branch("D_Coneptasy",&D_Coneptasy);
  newtree->Branch("Dst_Coneptasy",&Dst_Coneptasy);
  newtree->Branch("mHH", & mHH);
  

  //  for (Long64_t i=0;i<nentries; i++) {
  for (Long64_t i=0;i<200000; i++) {

    fChain->GetEntry(i);
    //aply trigger selection criteria and MC truth matching                                                                                       
    if(!MCTruthmatched()) continue;
    if(!isL0Selected() || !isHlt1Selected() || !isHlt2Selected() )continue;

    initializeMomenta();
    Slowpi_cosh=slowpi_helicityAngle();
    D_cosh = D0_helicityAngle();
    mu0_cosh=muon_helicityAngle();
    deltaM = Dst_DTF_Dstarplus_M - Dst_DTF_D0_M;
    D_Conemult = D_cmult_1_50;
    Dst_Conemult = Dst_cmult_1_50;
    D_Coneptasy = D_ptasy_1_50;
    Dst_Coneptasy = Dst_ptasy_1_50;
    mHH = (pH1+pH0).M();
    mMuMu_noDCF=get_mMuMu_noDCF();
    mMuMu_doubleDCF=get_mMuMu_doubleDCF();
    mMuMu_DCF_highP=get_mMuMu_DCF_highP();
    mMuMu_DCF_lowP=get_mMuMu_DCF_lowP();
    newtree->Fill();
  }

  std::cout<<"Created samples "<< name <<" with "<<newtree->GetEntries()<<" entries."<<std::endl;
  //newtree->Print();                                                                                                                             
  newtree->AutoSave();

  delete newfile;

}

