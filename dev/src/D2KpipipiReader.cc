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
  double h0_PIDK,h1_PIDK;
  Int_t mu0_MuonNShared,mu1_MuonNShared;
  
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
    (*it)->Branch("mu0_MuonNShared", &mu0_MuonNShared );
    (*it)->Branch("mu1_MuonNShared", &mu1_MuonNShared );
    (*it)->Branch("h0_PIDK", &h0_PIDK );
    (*it)->Branch("h1_PIDK", &h1_PIDK );


  }
  for (Long64_t i=0;i<nentries; i++) {

    if(i%10000 == 0) std::cout<<i<<" events processed...." << std::endl;

    fChain->GetEntry(i);
    //select only events that pass the trigger seelction criteria                                                                                                                         

    //if( Dst_DTF_Dstarplus_M - Dst_DTF_D0_M < 140 ||  Dst_DTF_Dstarplus_M - Dst_DTF_D0_M > 154 ) continue; //tighten here small window for 1D fit                 
    if( Dst_DTF_Dstarplus_M - Dst_DTF_D0_M < 144.5 ||  Dst_DTF_Dstarplus_M - Dst_DTF_D0_M > 146.5 ) continue; //tighten here small window for 1D fit
    if( Dst_DTF_D0_M < 1700 ||  Dst_DTF_D0_M > 2000 ) continue;
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


void D2KpipipiReader::createRandomizedSubsample(TString name) {//randomize pion pairs                                                                                                    

  Long64_t nentries = fChain->GetEntries();

  activateRelevantBranches();
  fChain->GetEntry(0);

  //Set trigger variables status to 0 bacuase they are not swapped                                                                                                                        
  fChain->SetBranchStatus("mu0_L0HadronDecision_TIS",0);
  fChain->SetBranchStatus("mu0_L0MuonDecision_TIS",0);
  fChain->SetBranchStatus("mu0_L0DiMuonDecision_TIS",0);
  fChain->SetBranchStatus("mu0_L0PhotonDecision_TIS",0);
  fChain->SetBranchStatus("mu0_L0ElectronDecision_TIS",0);

  fChain->SetBranchStatus("mu0_L0HadronDecision_TOS",0);
  fChain->SetBranchStatus("mu0_L0MuonDecision_TOS",0);
  fChain->SetBranchStatus("mu0_L0DiMuonDecision_TOS",0);
  fChain->SetBranchStatus("mu0_L0PhotonDecision_TOS",0);
  fChain->SetBranchStatus("mu0_L0ElectronDecision_TOS",0);

  fChain->SetBranchStatus("mu1_L0HadronDecision_TIS",0);
  fChain->SetBranchStatus("mu1_L0MuonDecision_TIS",0);
  fChain->SetBranchStatus("mu1_L0DiMuonDecision_TIS",0);
  fChain->SetBranchStatus("mu1_L0PhotonDecision_TIS",0);
  fChain->SetBranchStatus("mu1_L0ElectronDecision_TIS",0);

  fChain->SetBranchStatus("mu1_L0HadronDecision_TOS",0);
  fChain->SetBranchStatus("mu1_L0MuonDecision_TOS",0);
  fChain->SetBranchStatus("mu1_L0DiMuonDecision_TOS",0);
  fChain->SetBranchStatus("mu1_L0PhotonDecision_TOS",0);
  fChain->SetBranchStatus("mu1_L0ElectronDecision_TOS",0);

  fChain->SetBranchStatus("h1_L0HadronDecision_TIS",0);
  fChain->SetBranchStatus("h1_L0MuonDecision_TIS",0);
  fChain->SetBranchStatus("h1_L0DiMuonDecision_TIS",0);
  fChain->SetBranchStatus("h1_L0PhotonDecision_TIS",0);
  fChain->SetBranchStatus("h1_L0ElectronDecision_TIS",0);

  fChain->SetBranchStatus("h0_L0HadronDecision_TOS",0);
  fChain->SetBranchStatus("h0_L0MuonDecision_TOS",0);
  fChain->SetBranchStatus("h0_L0DiMuonDecision_TOS",0);
  fChain->SetBranchStatus("h0_L0PhotonDecision_TOS",0);
  fChain->SetBranchStatus("h0_L0ElectronDecision_TOS",0);

  fChain->SetBranchStatus("h0_L0HadronDecision_TIS",0);
  fChain->SetBranchStatus("h0_L0MuonDecision_TIS",0);
  fChain->SetBranchStatus("h0_L0DiMuonDecision_TIS",0);
  fChain->SetBranchStatus("h0_L0PhotonDecision_TIS",0);
  fChain->SetBranchStatus("h0_L0ElectronDecision_TIS",0);


  fChain->SetBranchStatus("h0_L0HadronDecision_TOS",0);
  fChain->SetBranchStatus("h0_L0MuonDecision_TOS",0);
  fChain->SetBranchStatus("h0_L0DiMuonDecision_TOS",0);
  fChain->SetBranchStatus("h0_L0PhotonDecision_TOS",0);
  fChain->SetBranchStatus("h0_L0ElectronDecision_TOS",0);

  fChain->SetBranchStatus("mu0_Hlt1TrackMuonDecision_Dec",0);
  fChain->SetBranchStatus("mu1_Hlt1TrackMuonDecision_Dec",0);

  fChain->SetBranchStatus("mu0_L0MuonDecision_TIS",0);
  fChain->SetBranchStatus("mu1_L0MuonDecision_TIS",0);
  fChain->SetBranchStatus("mu0_L0DiMuonDecision_TIS",0);
  fChain->SetBranchStatus("mu1_L0DiMuonDecision_TIS",0);
  fChain->SetBranchStatus("h1_L0MuonDecision_TIS",0);
  fChain->SetBranchStatus("h1_L0DiMuonDecision_TIS",0);

  fChain->SetBranchStatus("mu0_L0MuonDecision_TOS",0);
  fChain->SetBranchStatus("mu1_L0MuonDecision_TOS",0);
  fChain->SetBranchStatus("mu0_L0DiMuonDecision_TOS",0);
  fChain->SetBranchStatus("mu1_L0DiMuonDecision_TOS",0);
  fChain->SetBranchStatus("h1_L0MuonDecision_TOS",0);
  fChain->SetBranchStatus("h1_L0DiMuonDecision_TOS",0);

  TRandom3* generator = new TRandom3(runNumber); //Seed is RunNumber of first event                                                                                                       

  //Create a new file + a clone of old the in new file                                                                                                                                    

  TFile *newfile = new TFile(name,"recreate");
  TTree *newtree_data_odd = new TTree("DecayTree_odd","DecayTree_odd");
  TTree *newtree_data_even =  new TTree("DecayTree_even","DecayTree_even");

  TDirectory *f_data = newfile->mkdir("data");

  double Slowpi_cosh,mu0_cosh, D_cosh, deltaM,nVeloTracks;
  double D_Conemult,Dst_Conemult,D_Coneptasy,Dst_Coneptasy;
  Int_t nTracks_data,nPVs_data,nSPDHits;
  double mHH;
  bool isSideband;
  double mKpiOS,mKpiSS,mpipiOS,mpipiSS,misID_mD_OS,misID_mD_SS,misID_dm_SS,misID_dm_OS;
  double pKpiOS, pKpiSS, ppipiSS, ppipiOS;
  double h0_PIDK,h1_PIDK;

  std::vector<TTree*> trees;
  trees.push_back(newtree_data_odd);
  trees.push_back(newtree_data_even);
  
  double mup_PT;
  double mup_PX;
  double mup_PY;
  double mup_PZ;
  double mup_P;
  double mup_PIDe;
  double mup_PIDmu;
  double mup_PIDK;
  double mup_PIDp;
  double mup_ProbNNe;
  double mup_ProbNNk;
  double mup_ProbNNp;
  double mup_ProbNNpi;
  double mup_ProbNNmu;
  double mup_ProbNNghost;
  double mup_ID;
  double mup_isMuon;
  double mup_isMuonLoose;
  double mup_IP_ORIVX;
  double mup_IPCHI2_ORIVX;
  double mup_IP_OWNPV;
  double mup_IPCHI2_OWNPV;
  double mup_MINIP;
  double mup_MINIPCHI2;
  double mup_NShared;
  double mup_MuonNShared;
  double mup_TRACK_GhostProb;
  double mup_TRACK_CHI2NDOF;

  double hp_PT;
  double hp_PX;
  double hp_PY;
  double hp_PZ;
  double hp_P;
  double hp_PIDe;
  double hp_PIDmu;
  double hp_PIDK;
  double hp_PIDp;
  double hp_ProbNNe;
  double hp_ProbNNk;
  double hp_ProbNNp;
  double hp_ProbNNpi;
  double hp_ProbNNmu;
  double hp_ProbNNghost;
  double hp_ID;
  double hp_isMuon;
  double hp_isMuonLoose;
  double hp_IP_ORIVX;
  double hp_IPCHI2_OWNPV;
  double hp_MINIP;
  double hp_MINIPCHI2;
  double hp_NShared;
  double hp_MuonNShared;
  double hp_TRACK_GhostProb;
  double hp_TRACK_CHI2NDOF;
  double hp_IPCHI2_ORIVX;
  double hp_IP_OWNPV;

  double DTF_hp_PX;
  double DTF_hp_PY;
  double DTF_hp_PZ;
  double DTF_mup_PX;
  double DTF_mup_PY;
  double DTF_mup_PZ;

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
    (*it)->Branch("mpipiOS", & mpipiOS);
    (*it)->Branch("misID_mD_OS", & misID_mD_OS);
    (*it)->Branch("misID_dm_OS", & misID_dm_OS);
    (*it)->Branch("pKpiOS", &pKpiOS );
    (*it)->Branch("ppipiOS", & ppipiOS);
    (*it)->Branch("D_MAXDOCA",&D_MAXDOCA);
    (*it)->Branch("D_FDCHI2_OWNPV",&D_FDCHI2_OWNPV);
    (*it)->Branch("D_DIRA_OWNPV",&D_DIRA_OWNPV);
    (*it)->Branch("D_ENDVERTEX_CHI2",&D_ENDVERTEX_CHI2);
    (*it)->Branch("Slowpi_IPCHI2_OWNPV",&Slowpi_IPCHI2_OWNPV);
    (*it)->Branch("D_MINIPCHI2",&D_MINIPCHI2);
    (*it)->Branch("Slowpi_P",&Slowpi_P);
    (*it)->Branch("Slowpi_PT",&Slowpi_PT);
    (*it)->Branch("Dst_DTF_D0_M",&Dst_DTF_D0_M);
    (*it)->Branch("Dst_DTF_Dstarplus_M",&Dst_DTF_Dstarplus_M);

    (*it)->Branch("mu0_PT",&mup_PT);
    (*it)->Branch("mu0_PX",&mup_PX);
    (*it)->Branch("mu0_PY",&mup_PY);
    (*it)->Branch("mu0_PZ",&mup_PZ);
    (*it)->Branch("mu0_P",&mup_P);
    (*it)->Branch("mu0_PIDe",&mup_PIDe);
    (*it)->Branch("mu0_PIDmu",&mup_PIDmu);
    (*it)->Branch("mu0_PIDK",&mup_PIDK);
    (*it)->Branch("mu0_PIDp",&mup_PIDp);
    (*it)->Branch("mu0_ProbNNe",&mup_ProbNNe);
    (*it)->Branch("mu0_ProbNNk",&mup_ProbNNk);
    (*it)->Branch("mu0_ProbNNp",&mup_ProbNNp);
    (*it)->Branch("mu0_ProbNNpi",&mup_ProbNNpi);
    (*it)->Branch("mu0_ProbNNmu",&mup_ProbNNmu);
    (*it)->Branch("mu0_ProbNNghost",&mup_ProbNNghost);
    (*it)->Branch("mu0_ID",&mup_ID);
    (*it)->Branch("mu0_isMuon",&mup_isMuon);
    (*it)->Branch("mu0_isMuonLoose",&mup_isMuonLoose);
    (*it)->Branch("mu0_IP_ORIVX",&mup_IP_ORIVX);
    (*it)->Branch("mu0_IPCHI2_ORIVX",&mup_IPCHI2_ORIVX);
    (*it)->Branch("mu0_IP_OWNPV",&mup_IP_OWNPV);
    (*it)->Branch("mu0_IPCHI2_OWNPV",&mup_IPCHI2_OWNPV);
    (*it)->Branch("mu0_MINIP",&mup_MINIP);
    (*it)->Branch("mu0_MINIPCHI2",&mup_MINIPCHI2);
    (*it)->Branch("mu0_NShared",&mup_NShared);
    (*it)->Branch("mu0_MuonNShared",&mup_MuonNShared);
    (*it)->Branch("mu0_TRACK_GhostProb",&mup_TRACK_GhostProb);
    (*it)->Branch("mu0_TRACK_CHI2NDOF",&mup_TRACK_CHI2NDOF);

    (*it)->Branch("mu1_PT",&mu1_PT);
    (*it)->Branch("mu1_PX",&mu1_PX);
    (*it)->Branch("mu1_PY",&mu1_PY);
    (*it)->Branch("mu1_PZ",&mu1_PZ);
    (*it)->Branch("mu1_P",&mu1_P);
    (*it)->Branch("mu1_PIDe",&mu1_PIDe);
    (*it)->Branch("mu1_PIDmu",&mu1_PIDmu);
    (*it)->Branch("mu1_PIDK",&mu1_PIDK);
    (*it)->Branch("mu1_PIDp",&mu1_PIDp);
    (*it)->Branch("mu1_ProbNNe",&mu1_ProbNNe);
    (*it)->Branch("mu1_ProbNNk",&mu1_ProbNNk);
    (*it)->Branch("mu1_ProbNNp",&mu1_ProbNNp);
    (*it)->Branch("mu1_ProbNNpi",&mu1_ProbNNpi);
    (*it)->Branch("mu1_ProbNNmu",&mu1_ProbNNmu);
    (*it)->Branch("mu1_ProbNNghost",&mu1_ProbNNghost);
    (*it)->Branch("mu1_ID",&mu1_ID);
    (*it)->Branch("mu1_isMuon",&mu1_isMuon);
    (*it)->Branch("mu1_isMuonLoose",&mu1_isMuonLoose);
    (*it)->Branch("mu1_IP_ORIVX",&mu1_IP_ORIVX);
    (*it)->Branch("mu1_IPCHI2_ORIVX",&mu1_IPCHI2_ORIVX);
    (*it)->Branch("mu1_IP_OWNPV",&mu1_IP_OWNPV);
    (*it)->Branch("mu1_IPCHI2_OWNPV",&mu1_IPCHI2_OWNPV);
    (*it)->Branch("mu1_MINIP",&mu1_MINIP);
    (*it)->Branch("mu1_MINIPCHI2",&mu1_MINIPCHI2);
    (*it)->Branch("mu1_NShared",&mu1_NShared);
    (*it)->Branch("mu1_MuonNShared",&mu1_MuonNShared);
    (*it)->Branch("mu1_TRACK_GhostProb",&mu1_TRACK_GhostProb);
    (*it)->Branch("mu1_TRACK_CHI2NDOF",&mu1_TRACK_CHI2NDOF);

    (*it)->Branch("h0_PT",&h0_PT);
    (*it)->Branch("h0_PX",&h0_PX);
    (*it)->Branch("h0_PY",&h0_PY);
    (*it)->Branch("h0_PZ",&h0_PZ);
    (*it)->Branch("h0_P",&h0_P);
    (*it)->Branch("h0_PIDe",&h0_PIDe);
    (*it)->Branch("h0_PIDmu",&h0_PIDmu);
    (*it)->Branch("h0_PIDK",&h0_PIDK);
    (*it)->Branch("h0_PIDp",&h0_PIDp);
    (*it)->Branch("h0_ProbNNe",&h0_ProbNNe);
    (*it)->Branch("h0_ProbNNk",&h0_ProbNNk);
    (*it)->Branch("h0_ProbNNp",&h0_ProbNNp);
    (*it)->Branch("h0_ProbNNe",&h0_ProbNNe);
    (*it)->Branch("h0_ProbNNk",&h0_ProbNNk);
    (*it)->Branch("h0_ProbNNp",&h0_ProbNNp);
    (*it)->Branch("h0_ProbNNpi",&h0_ProbNNpi);
    (*it)->Branch("h0_ProbNNmu",&h0_ProbNNmu);
    (*it)->Branch("h0_ProbNNghost",&h0_ProbNNghost);
    (*it)->Branch("h0_ID",&h0_ID);
    (*it)->Branch("h0_isMuon",&h0_isMuon);
    (*it)->Branch("h0_isMuonLoose",&h0_isMuonLoose);
    (*it)->Branch("h0_IP_ORIVX",&h0_IP_ORIVX);
    (*it)->Branch("h0_IPCHI2_ORIVX",&h0_IPCHI2_ORIVX);
    (*it)->Branch("h0_IP_OWNPV",&h0_IP_OWNPV);
    (*it)->Branch("h0_IPCHI2_OWNPV",&h0_IPCHI2_OWNPV);
    (*it)->Branch("h0_MINIP",&h0_MINIP);
    (*it)->Branch("h0_MINIPCHI2",&h0_MINIPCHI2);
    (*it)->Branch("h0_NShared",&h0_NShared);
    (*it)->Branch("h0_MuonNShared",&h0_MuonNShared);
    (*it)->Branch("h0_TRACK_GhostProb",&h0_TRACK_GhostProb);
    (*it)->Branch("h0_TRACK_CHI2NDOF",&h0_TRACK_CHI2NDOF);

    (*it)->Branch("h1_PT",&hp_PT);
    (*it)->Branch("h1_PX",&hp_PX);
    (*it)->Branch("h1_PY",&hp_PY);
    (*it)->Branch("h1_PZ",&hp_PZ);
    (*it)->Branch("h1_P",&hp_P);
    (*it)->Branch("h1_PIDe",&hp_PIDe);
    (*it)->Branch("h1_PIDmu",&hp_PIDmu);
    (*it)->Branch("h1_PIDK",&hp_PIDK);
    (*it)->Branch("h1_PIDp",&hp_PIDp);
    (*it)->Branch("h1_ProbNNe",&hp_ProbNNe);
    (*it)->Branch("h1_ProbNNk",&hp_ProbNNk);
    (*it)->Branch("h1_ProbNNp",&hp_ProbNNp);
    (*it)->Branch("h1_ProbNNpi",&hp_ProbNNpi);
    (*it)->Branch("h1_ProbNNmu",&hp_ProbNNmu);
    (*it)->Branch("h1_ProbNNghost",&hp_ProbNNghost);
    (*it)->Branch("h1_ID",&hp_ID);
    (*it)->Branch("h1_isMuon",&hp_isMuon);
    (*it)->Branch("h1_isMuonLoose",&hp_isMuonLoose);
    (*it)->Branch("h1_IP_ORIVX",&hp_IP_ORIVX);
    (*it)->Branch("h1_IPCHI2_ORIVX",&hp_IPCHI2_ORIVX);
    (*it)->Branch("h1_IP_OWNPV",&hp_IP_OWNPV);
    (*it)->Branch("h1_IPCHI2_OWNPV",&hp_IPCHI2_OWNPV);
    (*it)->Branch("h1_MINIP",&hp_MINIP);
    (*it)->Branch("h1_MINIPCHI2",&hp_MINIPCHI2);
    (*it)->Branch("h1_NShared",&hp_NShared);
    (*it)->Branch("h1_MuonNShared",&hp_MuonNShared);
    (*it)->Branch("h1_TRACK_GhostProb",&hp_TRACK_GhostProb);
    (*it)->Branch("h1_TRACK_CHI2NDOF",&hp_TRACK_CHI2NDOF);

    (*it)->Branch("Dst_DTF_h1_PX",&DTF_hp_PX);
    (*it)->Branch("Dst_DTF_h1_PY",&DTF_hp_PY);
    (*it)->Branch("Dst_DTF_h1_PZ",&DTF_hp_PZ);
    (*it)->Branch("Dst_DTF_h0_PX",&Dst_DTF_h0_PX);
    (*it)->Branch("Dst_DTF_h0_PY",&Dst_DTF_h0_PY);
    (*it)->Branch("Dst_DTF_h0_PZ",&Dst_DTF_h0_PZ);
    (*it)->Branch("Dst_DTF_mu0_PX",&DTF_mup_PX);
    (*it)->Branch("Dst_DTF_mu0_PY",&DTF_mup_PY);
    (*it)->Branch("Dst_DTF_mu0_PZ",&DTF_mup_PZ);
    (*it)->Branch("Dst_DTF_mu1_PX",&Dst_DTF_mu1_PX);
    (*it)->Branch("Dst_DTF_mu1_PY",&Dst_DTF_mu1_PY);
    (*it)->Branch("Dst_DTF_mu1_PZ",&Dst_DTF_mu1_PZ);

    (*it)->Branch("D_DiMuon_Mass",&D_DiMuon_Mass);

  }


    for (Long64_t i=0;i<nentries; i++) {

      if(i%10000 == 0) std::cout<<i<<" events processed...." << std::endl;

      fChain->GetEntry(i);
      //select only events that pass the trigger seelction criteria                                                                                                                         
      //if( Dst_DTF_Dstarplus_M - Dst_DTF_D0_M < 140 ||  Dst_DTF_Dstarplus_M - Dst_DTF_D0_M > 154 ) continue; //tighten here small window for 1D fit                                        
      if( Dst_DTF_Dstarplus_M - Dst_DTF_D0_M < 144.5 ||  Dst_DTF_Dstarplus_M - Dst_DTF_D0_M > 146.5 ) continue; //tighten here small window for 1D fit                                      
      if( Dst_DTF_D0_M < 1700 ||  Dst_DTF_D0_M > 2000 ) continue;
      if(!passGhostProbCut(0.5)) continue;

      double randomNr1 = generator->Rndm();
      //assign randomly one of the positively charged pions as muon and the other one as pion                                                                                               
      if(randomNr1>0.5) {

	//std::cout<<"mup = mu0"<<std::endl;                                                                                                                                                
	mup_PT = mu0_PT;
	mup_PX= mu0_PX;
	mup_PY= mu0_PY;
	mup_PZ=mu0_PZ;
	mup_P=mu0_P;
	mup_PIDe=mu0_PIDe;
	mup_PIDmu=mu0_PIDmu;
	mup_PIDK=mu0_PIDK;
	mup_PIDp=mu0_PIDp;
	mup_ProbNNe=mu0_ProbNNe;
	mup_ProbNNk=mu0_ProbNNk;
	mup_ProbNNp=mu0_ProbNNp;
	mup_ProbNNpi=mu0_ProbNNpi;
	mup_ProbNNmu=mu0_ProbNNmu;
	mup_ProbNNghost=mu0_ProbNNghost;
	mup_ID=mu0_ID;
	mup_isMuon= mu0_isMuon;
	mup_isMuonLoose=mu0_isMuonLoose;
	mup_IP_ORIVX=mu0_IP_ORIVX;
	mup_IPCHI2_ORIVX=mu0_IPCHI2_ORIVX;
	mup_IP_OWNPV=mu0_IP_OWNPV;
	mup_IPCHI2_OWNPV=mu0_IPCHI2_OWNPV;
	mup_MINIP=mu0_MINIP;
	mup_MINIPCHI2=mu0_MINIPCHI2;
	mup_NShared=mu0_NShared;
	mup_MuonNShared=mu0_MuonNShared;
	mup_TRACK_GhostProb= mu0_TRACK_GhostProb;
	mup_TRACK_CHI2NDOF=mu0_TRACK_CHI2NDOF;

	hp_PT = h1_PT;
	hp_PX= h1_PX;
	hp_PY= h1_PY;
	hp_PZ=h1_PZ;
	hp_P=h1_P;
	hp_PIDe=h1_PIDe;
	hp_PIDmu=h1_PIDmu;
	hp_PIDK=h1_PIDK;
	hp_PIDp=h1_PIDp;
	hp_ProbNNe=h1_ProbNNe;
	hp_ProbNNk=h1_ProbNNk;
	hp_ProbNNp=h1_ProbNNp;
	hp_ProbNNpi=h1_ProbNNpi;
	hp_ProbNNmu=h1_ProbNNmu;
	hp_ProbNNghost=h1_ProbNNghost;
	hp_ID=h1_ID;
	hp_isMuon= h1_isMuon;
	hp_isMuonLoose=h1_isMuonLoose;
	hp_IP_ORIVX=h1_IP_ORIVX;
	hp_IPCHI2_ORIVX=h1_IPCHI2_ORIVX;
	hp_IP_OWNPV=h1_IP_OWNPV;
	hp_IPCHI2_OWNPV=h1_IPCHI2_OWNPV;
	hp_MINIP=h1_MINIP;
	hp_MINIPCHI2=h1_MINIPCHI2;
	hp_NShared=h1_NShared;
	hp_MuonNShared=h1_MuonNShared;
	hp_TRACK_GhostProb= h1_TRACK_GhostProb;
	hp_TRACK_CHI2NDOF=h1_TRACK_CHI2NDOF;

	DTF_hp_PX = Dst_DTF_h1_PX;
	DTF_hp_PY = Dst_DTF_h1_PY;
	DTF_hp_PZ=Dst_DTF_h1_PZ;
	DTF_mup_PX =Dst_DTF_mu0_PX;
	DTF_mup_PY=Dst_DTF_mu0_PY;
	DTF_mup_PZ=Dst_DTF_mu0_PZ;

      }

      else {

	mup_PT = h1_PT;
	mup_PX= h1_PX;
	mup_PY= h1_PY;
	mup_PZ=h1_PZ;
	mup_P=h1_P;
	mup_PIDe=h1_PIDe;
	mup_PIDmu=h1_PIDmu;
	mup_PIDK=h1_PIDK;
	mup_PIDp=h1_PIDp;
	mup_ProbNNe=h1_ProbNNe;
	mup_ProbNNk=h1_ProbNNk;
	mup_ProbNNp=h1_ProbNNp;
	mup_ProbNNpi=h1_ProbNNpi;
	mup_ProbNNmu=h1_ProbNNmu;
	mup_ProbNNghost=h1_ProbNNghost;
	mup_ID=h1_ID;
	mup_isMuon= h1_isMuon;
	mup_isMuonLoose=h1_isMuonLoose;
	mup_IP_ORIVX=h1_IP_ORIVX;
	mup_IPCHI2_ORIVX=h1_IPCHI2_ORIVX;
	mup_IP_OWNPV=h1_IP_OWNPV;
	mup_IPCHI2_OWNPV=h1_IPCHI2_OWNPV;
	mup_MINIP=h1_MINIP;
	mup_MINIPCHI2=h1_MINIPCHI2;
	mup_NShared=h1_NShared;
	mup_MuonNShared=h1_MuonNShared;
	mup_TRACK_GhostProb= h1_TRACK_GhostProb;
	mup_TRACK_CHI2NDOF=h1_TRACK_CHI2NDOF;

	hp_PT = mu0_PT;
	hp_PX= mu0_PX;
	hp_PY= mu0_PY;
	hp_PZ=mu0_PZ;
	hp_P=mu0_P;
	hp_PIDe=mu0_PIDe;
	hp_PIDmu=mu0_PIDmu;
	hp_PIDK=mu0_PIDK;
	hp_PIDp=mu0_PIDp;
	hp_ProbNNe=mu0_ProbNNe;
	hp_ProbNNk=mu0_ProbNNk;
	hp_ProbNNp=mu0_ProbNNp;
	hp_ProbNNpi=mu0_ProbNNpi;
	hp_ProbNNmu=mu0_ProbNNmu;
	hp_ProbNNghost=mu0_ProbNNghost;
	hp_ID=mu0_ID;
	hp_isMuon= mu0_isMuon;
	hp_isMuonLoose=mu0_isMuonLoose;
	hp_IP_ORIVX=mu0_IP_ORIVX;
	hp_IPCHI2_ORIVX=mu0_IPCHI2_ORIVX;
	hp_IP_OWNPV=mu0_IP_OWNPV;
	hp_IPCHI2_OWNPV=mu0_IPCHI2_OWNPV;
	hp_MINIP=mu0_MINIP;
	hp_MINIPCHI2=mu0_MINIPCHI2;
	hp_NShared=mu0_NShared;
	hp_MuonNShared=mu0_MuonNShared;
	hp_TRACK_GhostProb= mu0_TRACK_GhostProb;
	hp_TRACK_CHI2NDOF=mu0_TRACK_CHI2NDOF;

	DTF_hp_PX = Dst_DTF_mu0_PX;
	DTF_hp_PY = Dst_DTF_mu0_PY;
	DTF_hp_PZ=Dst_DTF_mu0_PZ;
	DTF_mup_PX =Dst_DTF_h1_PX;
	DTF_mup_PY=Dst_DTF_h1_PY;
	DTF_mup_PZ=Dst_DTF_h1_PZ;

      }

      //std::cout<<mup_PT<<" "<<mu0_PT<<"  "<<h0_PT<<std::endl;                                                                                                                             

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

      pDTFH1.SetXYZM(DTF_hp_PX,DTF_hp_PY,DTF_hp_PZ,Mass::Pi());
      pDTFH0.SetXYZM(Dst_DTF_h0_PX,Dst_DTF_h0_PY,Dst_DTF_h0_PZ,Mass::K());
      pDTFMu1.SetXYZM(Dst_DTF_mu1_PX,Dst_DTF_mu1_PY,Dst_DTF_mu1_PZ,Mass::Mu());
      pDTFMu0.SetXYZM(DTF_mup_PX,DTF_mup_PY,DTF_mup_PZ,Mass::Mu());
      pDTFPis.SetXYZM(Dst_DTF_Pis_PX,Dst_DTF_Pis_PY,Dst_DTF_Pis_PZ,Mass::Pi());

      pDTFPi1.SetXYZM(Dst_DTF_mu1_PX,Dst_DTF_mu1_PY,Dst_DTF_mu1_PZ,Mass::Pi());
      pDTFPi0.SetXYZM(DTF_mup_PX,DTF_mup_PY,DTF_mup_PZ,Mass::Pi());

      D_DiMuon_Mass = (pDTFMu1+pDTFMu0).M();
      mHH = (pH1+pH0).M();

      mpipiOS =(pDTFPi0+pDTFPi1).M();

      misID_mD_OS=( pDTFH1+pDTFH0+pDTFMu0+pDTFMu1).M();
      misID_dm_OS=( pDTFH1+pDTFH0+pDTFMu0+pDTFMu1+pDTFPis).M() -  ( pDTFH1+pDTFH0+pDTFMu0+pDTFMu1).M();

      pKpiOS=(pDTFH0+pDTFH1).P();
      ppipiOS=(pDTFPi0+pDTFPi1).P();

      if(eventNumber%2!=0){ //odd events                                                                                                                                                    
	trees[0]->Fill(); //odd data                                                                                                                                                        
      }
      if(eventNumber%2==0){ //even events                                                                                                                                                   
	trees[1]->Fill(); //even data                                                                                                                                                       
      }
    }

    f_data->cd();
    trees[0]->Write();
    trees[1]->Write();

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
  double h0_PIDK,h1_PIDK;
  Int_t mu0_MuonNShared,mu1_MuonNShared;


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
  newtree_even->Branch("mu0_MuonNShared", &mu0_MuonNShared );
  newtree_even->Branch("mu1_MuonNShared", &mu1_MuonNShared );
  newtree_even->Branch("h0_PIDK", &h0_PIDK );
  newtree_even->Branch("h1_PIDK", &h1_PIDK );


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
  newtree_odd->Branch("mu0_MuonNShared", &mu0_MuonNShared );
  newtree_odd->Branch("mu1_MuonNShared", &mu1_MuonNShared );
  newtree_odd->Branch("h0_PIDK", &h0_PIDK );
  newtree_odd->Branch("h1_PIDK", &h1_PIDK );


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

