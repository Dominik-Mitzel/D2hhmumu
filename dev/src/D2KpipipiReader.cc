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
  pDTFMu2.SetXYZM(Dst_DTF_h1_PX,Dst_DTF_h1_PY,Dst_DTF_h1_PZ,Mass::Mu());

  pDTFPi1.SetXYZM(Dst_DTF_mu1_PX,Dst_DTF_mu1_PY,Dst_DTF_mu1_PZ,Mass::Pi());
  pDTFPi0.SetXYZM(Dst_DTF_mu0_PX,Dst_DTF_mu0_PY,Dst_DTF_mu0_PZ,Mass::Pi());

  pD.SetXYZM(D_PX,D_PY,D_PZ,Mass::D0());
  pDst.SetXYZM(Dst_PX,Dst_PY,Dst_PZ,Mass::Ds());
  pPis.SetXYZM(Slowpi_PX,Slowpi_PY,Slowpi_PZ,Mass::Pi());

  pDTFDst.SetXYZM(Dst_DTF_Dstarplus_PX,Dst_DTF_Dstarplus_PY,Dst_DTF_Dstarplus_PZ,Mass::Ds());
  pDTFD.SetXYZM(Dst_DTF_D0_PX,Dst_DTF_D0_PY,Dst_DTF_D0_PZ,Mass::D0());
  pDTFPis.SetXYZM(Dst_DTF_Pis_PX,Dst_DTF_Pis_PY,Dst_DTF_Pis_PZ,Mass::Pi());



}

void D2KpipipiReader::createSubsample(TString name, double percentage) {//specify percentage of full data sample written in sideband subsample                                     

  Long64_t nentries = fChain->GetEntries();
  std::cout<<"This programm creates ("<<percentage<<"%) subsample in mass sideband and applies trigger selection to the rest"<<std::endl;
  std::cout<<"Found tree with "<<nentries <<" entries....start random sampling ("<<percentage<<"%)..."<<std::endl;

  activateRelevantBranches();
  fChain->GetEntry(0);

  TRandom3* generator = new TRandom3(runNumber); //Seed is RunNumber of first event                                                                                                        

  //Create a new file + a clone of old tree in new file                                                                                                                                   
 
  TFile *newfile = new TFile(name,"recreate");
  TDirectory *f_sideband = newfile->mkdir("sideband");
  TTree *newtree_data_odd = fChain->CloneTree(0);
  TTree *newtree_data_even = fChain->CloneTree(0);
  newtree_data_odd->SetName("DecayTree_odd");
  newtree_data_even->SetName("DecayTree_even");

  TDirectory *f_data = newfile->mkdir("data");
  TTree *newtree_bkg_even = fChain->CloneTree(0);
  TTree *newtree_bkg_odd = fChain->CloneTree(0);
  newtree_bkg_odd->SetName("DecayTree_odd");
  newtree_bkg_even->SetName("DecayTree_even");
  double Slowpi_cosh,mu0_cosh, D_cosh, deltaM,nVeloTracks;
  double D_Conemult,Dst_Conemult,D_Coneptasy,Dst_Coneptasy;
  Int_t nTracks_data,nPVs_data,nSPDHits;
  double mHH;
  bool isSideband;
  
  double mKpiOS,mKpiSS,mpipiOS,mpipiSS,misID_mD_OS,misID_mD_SS,misID_dm_SS,misID_dm_OS; 
  double pKpiOS, pKpiSS, ppipiSS, ppipiOS;

  std::vector<TTree*> trees;
  trees.push_back(newtree_data_odd);
  trees.push_back(newtree_data_even);
  trees.push_back(newtree_bkg_odd);
  trees.push_back(newtree_bkg_even);

  for (std::vector<TTree*>::iterator it = trees.begin() ; it != trees.end(); ++it) {

    (*it)->Branch("isBkgSideband",&isSideband);
    (*it)->Branch("Slowpi_cosh",&Slowpi_cosh);
    (*it)->Branch("D_cosh", & D_cosh);
    (*it)->Branch("mu0_cosh",&mu0_cosh);
    (*it)->Branch("deltaM",&deltaM);
    (*it)->Branch("D_Conemult",&D_Conemult);
    (*it)->Branch("Dst_Conemult",&Dst_Conemult);
    (*it)->Branch("D_Coneptasy",&D_Coneptasy);
    (*it)->Branch("Dst_Coneptasy",&Dst_Coneptasy);
    (*it)->Branch("nSPDHits",&nSPDHits);
    (*it)->Branch("nVeloTracks",&nVeloTracks);
    (*it)->Branch("nTracks_data",&nTracks_data);
    (*it)->Branch("mKpiOS", & mKpiOS);
    (*it)->Branch("mKpiSS", & mKpiSS);
    (*it)->Branch("mpipiOS", & mpipiOS);
    (*it)->Branch("mpipiSS", & mpipiSS);
    (*it)->Branch("misID_mD_OS", & misID_mD_OS);
    (*it)->Branch("misID_mD_SS", &misID_mD_SS);   
    (*it)->Branch("misID_dm_SS", & misID_dm_SS);
    (*it)->Branch("misID_dm_OS", & misID_dm_OS);
    (*it)->Branch("pKpiOS", &pKpiOS );
    (*it)->Branch("pKpiSS", & pKpiSS);
    (*it)->Branch("ppipiOS", & ppipiOS);
    (*it)->Branch("ppipiSS", & ppipiSS);


  }
  for (Long64_t i=0;i<nentries; i++) {

    if(i%10000 == 0) std::cout<<i<<" events processed...." << std::endl;

    fChain->GetEntry(i);
    //select only events that pass the trigger seelction criteria                                                                                                                         

    if( Dst_DTF_Dstarplus_M - Dst_DTF_D0_M < 140 ||  Dst_DTF_Dstarplus_M - Dst_DTF_D0_M > 154 ) continue;
    if( Dst_DTF_D0_M < 1760 ||  Dst_DTF_D0_M > 1980 ) continue;
    if(!passGhostProbCut(0.5)) continue;


    //if(!isL0Selected() || !isHlt1Selected() || !isHlt2Selected()  )continue;
    initializeMomenta();

    //additional variables                                                                                                                                                                  
    Slowpi_cosh=slowpi_helicityAngle();
    D_cosh = D0_helicityAngle();
    mu0_cosh=muon_helicityAngle();
    deltaM = Dst_DTF_Dstarplus_M - Dst_DTF_D0_M;
    D_Conemult = Dst_CONEMULT_D;
    Dst_Conemult = Dst_CONEMULT_Dstar;
    D_Coneptasy = Dst_CONEPTASYM_D;
    Dst_Coneptasy = Dst_CONEPTASYM_Dstar;
    nSPDHits = (int)nSpdDigits;

    nVeloTracks=nVELO;
    nTracks_data=(int)nTracks;
    nPVs_data=(int)nPVs;
    mHH = (pH1+pH0).M();

    mKpiOS = (pDTFH0+pDTFH1).M();
    mKpiSS = (pDTFH0+pDTFPi1).M();
    mpipiOS =(pDTFPi0+pDTFPi1).M();
    mpipiSS =(pDTFH1+pDTFPi0).M();

    misID_mD_OS=( pDTFH1+pDTFH0+pDTFMu0+pDTFMu1).M();
    misID_mD_SS=( pDTFPi1 + pDTFH0+pDTFMu0 + pDTFMu2).M();
    misID_dm_SS=(pDTFPi1 + pDTFH0+pDTFMu0 + pDTFMu2 + pDTFPis).M()- (pDTFPi1 + pDTFH0+pDTFMu0 + pDTFMu2).M();
    misID_dm_OS=( pDTFH1+pDTFH0+pDTFMu0+pDTFMu1+pDTFPis).M() -  ( pDTFH1+pDTFH0+pDTFMu0+pDTFMu1).M();

    pKpiOS=(pDTFH0+pDTFH1).P();
    pKpiSS=(pDTFH0+pDTFPi1).P();
    ppipiOS=(pDTFPi0+pDTFPi1).P();
    ppipiSS=(pDTFH1+pDTFPi0).P();

    if(eventNumber%2!=0){ //odd events                                                                                                                                                      
      trees[0]->Fill(); //odd data                                                                                                                                                          
      if(generator->Rndm()<percentage/100 && isBkgSideband()) trees[2]->Fill();
      }
     
    if(eventNumber%2==0){ //even events                                                                                                                                                    
      trees[1]->Fill(); //even data                                                                                                                                                        
      if(generator->Rndm()<percentage/100 && isBkgSideband()) trees[3]->Fill();
      }
      
  }

   
  f_data->cd();
  trees[0]->Write();
  trees[1]->Write();
  
  f_sideband->cd();
  trees[2]->Write();
  trees[3]->Write();

  delete newfile;

}

void D2KpipipiReader::createMCtrainingSample(TString name) {

  InitMC();
  //activateRelevantBranches();                                                                                                                                                               

  Long64_t nentries = fChain->GetEntries();
  std::cout<<"Found tree with "<<nentries <<" entries..."<<std::endl;

  fChain->GetEntry(0);

  //Create a new file + a clone of old tree in new file                                                                                                                                       
  TFile *newfile = new TFile(name,"recreate");
  TTree *newtree_even = fChain->CloneTree(0);
  TTree *newtree_odd = fChain->CloneTree(0);
  newtree_odd->SetName("DecayTree_odd");
  newtree_even->SetName("DecayTree_even ");


  double Slowpi_cosh,mu0_cosh,D_cosh, deltaM;
  double D_Conemult,Dst_Conemult,D_Coneptasy,Dst_Coneptasy;
  double mHH;
  double mKpiOS,mKpiSS,mpipiOS,mpipiSS,misID_mD_OS,misID_mD_SS,misID_dm_SS,misID_dm_OS;
  double pKpiOS, pKpiSS, ppipiSS, ppipiOS;


  newtree_even->Branch("Slowpi_cosh",&Slowpi_cosh);
  newtree_even->Branch("D_cosh", & D_cosh);
  newtree_even->Branch("mu0_cosh",&mu0_cosh);
  newtree_even->Branch("deltaM",&deltaM);
  newtree_even->Branch("D_Conemult",&D_Conemult);
  newtree_even->Branch("Dst_Conemult",&Dst_Conemult);
  newtree_even->Branch("D_Coneptasy",&D_Coneptasy);
  newtree_even->Branch("Dst_Coneptasy",&Dst_Coneptasy);
  newtree_even->Branch("mHH", & mHH);
  newtree_even->Branch("mKpiOS", & mKpiOS);
  newtree_even->Branch("mKpiSS", & mKpiSS);
  newtree_even->Branch("mpipiOS", & mpipiOS);
  newtree_even->Branch("mpipiSS", & mpipiSS);
  newtree_even->Branch("misID_mD_OS", & misID_mD_OS);
  newtree_even->Branch("misID_mD_SS", &misID_mD_SS);
  newtree_even->Branch("misID_dm_SS", & misID_dm_SS);
  newtree_even->Branch("misID_dm_OS", & misID_dm_OS);
  newtree_even->Branch("pKpiOS", &pKpiOS );
  newtree_even->Branch("pKpiSS", & pKpiSS);
  newtree_even->Branch("ppipiOS", & ppipiOS);
  newtree_even->Branch("ppipiSS", & ppipiSS);


  newtree_odd->Branch("Slowpi_cosh",&Slowpi_cosh);
  newtree_odd->Branch("D_cosh", & D_cosh);
  newtree_odd->Branch("mu0_cosh",&mu0_cosh);
  newtree_odd->Branch("deltaM",&deltaM);
  newtree_odd->Branch("D_Conemult",&D_Conemult);
  newtree_odd->Branch("Dst_Conemult",&Dst_Conemult);
  newtree_odd->Branch("D_Coneptasy",&D_Coneptasy);
  newtree_odd->Branch("Dst_Coneptasy",&Dst_Coneptasy);
  newtree_odd->Branch("mHH", & mHH);
  newtree_odd->Branch("mKpiOS", & mKpiOS);
  newtree_odd->Branch("mKpiSS", & mKpiSS);
  newtree_odd->Branch("mpipiOS", & mpipiOS);
  newtree_odd->Branch("mpipiSS", & mpipiSS);
  newtree_odd->Branch("misID_mD_OS", & misID_mD_OS);
  newtree_odd->Branch("misID_mD_SS", &misID_mD_SS);
  newtree_odd->Branch("misID_dm_SS", & misID_dm_SS);
  newtree_odd->Branch("misID_dm_OS", & misID_dm_OS);
  newtree_odd->Branch("pKpiOS", &pKpiOS );
  newtree_odd->Branch("pKpiSS", & pKpiSS);
  newtree_odd->Branch("ppipiOS", & ppipiOS);
  newtree_odd->Branch("ppipiSS", & ppipiSS);


  //  for (Long64_t i=0;i<nentries; i++) {
  for (Long64_t i=0;i<nentries; i++) {
    fChain->GetEntry(i);
    //aply trigger selection criteria and MC truth matching                                                                                                                                  
    if(!MCTruthmatched()) continue;
    if(!passGhostProbCut(0.5)) continue;

    //    if(!isL0Selected() || !isHlt1Selected() || !isHlt2Selected() )continue;

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
    mKpiOS = (pDTFH0+pDTFH1).M();
    mKpiSS = (pDTFH0+pDTFPi1).M();
    mpipiOS =(pDTFPi0+pDTFPi1).M();
    mpipiSS =(pDTFH1+pDTFPi0).M();
    misID_mD_OS=( pDTFH1+pDTFH0+pDTFMu0+pDTFMu1).M();
    misID_mD_SS=( pDTFPi1 + pDTFH0+pDTFMu0 + pDTFMu2).M();
    misID_dm_SS=(pDTFPi1 + pDTFH0+pDTFMu0 + pDTFMu2 + pDTFPis).M()- (pDTFPi1 + pDTFH0+pDTFMu0 + pDTFMu2).M();
    misID_dm_OS=( pDTFH1+pDTFH0+pDTFMu0+pDTFMu1+pDTFPis).M() -  ( pDTFH1+pDTFH0+pDTFMu0+pDTFMu1).M();
    pKpiOS=(pDTFH0+pDTFH1).P();
    pKpiSS=(pDTFH0+pDTFPi1).P();
    ppipiOS=(pDTFPi0+pDTFPi1).P();
    ppipiSS=(pDTFH1+pDTFPi0).P();


    if(eventNumber%2==0) newtree_even->Fill();
    else newtree_odd->Fill();
  }

  std::cout<<"Created MC taining samples "<< name <<" with "<<newtree_even->GetEntries()<<" entries (even sample) and "<< newtree_odd->GetEntries() <<"(odd.)"<<std::endl;

  //newtree->Print();                                                                                                                                                                        
  newtree_even->AutoSave();
  newtree_odd->AutoSave();

  delete newfile;

}

bool D2KpipipiReader::isHlt2Selected(){

  // if(D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS==1) return true;
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
    //if(!isL0Selected() || !isHlt1Selected() || !isHlt2Selected() )continue;

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

