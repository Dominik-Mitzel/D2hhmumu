#include "D2hhmumuReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include "Tools.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>

bool D2hhmumuReader::isBkgSideband(){

  if (Dst_DTF_D0_M < 1840 )return false; //cut these events also in event loop pr preselection... 

  //  if(  !((D_M > 1840 && D_M<1880) && (DTFdm() > 144.5 && DTFdm() < 146.5)) ) return true;
  if(  !((Dst_DTF_D0_M > 1840 && Dst_DTF_D0_M<1880) && (DTFdm() > 144.5 && DTFdm() < 146.5)) ) return true;                                                                                                   
  else return false; 	 
}

double D2hhmumuReader::slowpi_helicityAngle() {

  initializeMomenta();
  
  TVector3 b=(pDTFD+pDTFPis).BoostVector();
  TLorentzVector pcm(pDTFPis);
  pcm.Boost(-b);
  double Pis_cosh=pcm.Vect().Dot(b)/(pcm.P()*b.Mag());
  return Pis_cosh;
}  


double D2hhmumuReader::muon_helicityAngle() {

  initializeMomenta();
  
  TVector3 b=(pDTFD).BoostVector();
  TLorentzVector pcm(pDTFMu0);
  pcm.Boost(-b);
  double mu0_cosh=pcm.Vect().Dot(b)/(pcm.P()*b.Mag());
  return mu0_cosh;
}   

double D2hhmumuReader::D0_helicityAngle() {

  initializeMomenta();
  
  TVector3 b=(pDTFDst).BoostVector();
  TLorentzVector pcm(pDTFD);
  pcm.Boost(-b);
  double D_cosh=pcm.Vect().Dot(b)/(pcm.P()*b.Mag());
  return D_cosh;
}   
 



bool D2hhmumuReader::isInMassRange(){

if(Dst_DTF_D0_M < 1800 && Dst_DTF_D0_M > 1930 ) return false;
return true;
}

bool D2hhmumuReader::isL0Selected(){

  if(mu0_L0MuonDecision_TOS==1 || mu1_L0MuonDecision_TOS ==1 || mu0_L0DiMuonDecision_TOS==1 || mu1_L0DiMuonDecision_TOS ==1 || Dst_L0Global_TIS==1) return true;
  else return false; 	 
}

bool D2hhmumuReader::isHlt1Selected(){

  if(mu0_Hlt1TrackMuonDecision_TOS==1 || mu1_Hlt1TrackMuonDecision_TOS==1 || D_Hlt1TrackAllL0Decision_TOS==1) return true;
  else return false; 	 
}


bool D2hhmumuReader::MCTruthmatched(){

  if(Dst_BKGCAT > 10) return false;

  return true;
}

double D2hhmumuReader::DTFdm(){
  
  initializeMomenta();
  return  (pDTFH1+pDTFH0+pDTFMu0+pDTFMu1+pDTFPis).M()-(pDTFH1+pDTFH0+pDTFMu0+pDTFMu1).M();
}

double D2hhmumuReader::dm(){
  
  initializeMomenta();
  return (pH1+pH0+pMu0+pMu1+pPis).M()-(pH1+pH0+pMu0+pMu1).M();
}

void D2hhmumuReader::fillHistograms(TString fname, bool isMC=false) {

  TFile *newfile = new TFile(fname,"recreate");

  h1_DTFdm = new TH1D("h1_DTFdm","delta mass (DTF)",100,140,150);
  //TH1* h1_DTFdmtest = new TH1D("h1_DTFdmtest","delta mass (DTF)",100,140,150);
  h1_dm = new TH1D("h1_dm","delta mass ",100,140,150);
  h1_DTFmMuMu = new TH1D("h1_DTFMuMu","Dimuon mass (DTF)",100,100,900);
  h1_mMuMu = new TH1D("h1_mMuMu","Dimuon mass ",100,100,1200);
  h1_DTFmD = new TH1D("h1_DTFmD","reconstructed D mass (DTF)",50,1820,1900); 
  h1_mD = new TH1D("h1_mD","reconstucted mass",50,1820,1900);
  h2_mD_deltaM = new TH2D("h1_md_deltaM","deltaM vs. D mass",30,1840,1900,30,140,150);

     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(isMC && !MCTruthmatched()) continue;
      if(mu1_ProbNNmu<0.5 || mu0_ProbNNmu<0.5) continue;
      initializeMomenta(); 
      activateRelevantBranches();
      //if(!isL0Selected())continue;    
      //if(!isHlt1Selected())continue;
      //if(!isHlt2Selected())continue;       
      if(DTFdm()<144.5 || DTFdm() >146) continue;
 
   
      h1_DTFdm->Fill(DTFdm());
      h1_dm->Fill(dm());
      h1_DTFmMuMu->Fill((pDTFMu0+pDTFMu1).M());
      h1_mMuMu->Fill((pMu0+pMu1).M());	
      h1_DTFmD->Fill(Dst_DTF_D0_M );
      h1_mD->Fill((pH1+pH0+pMu0+pMu1).M());
      h2_mD_deltaM->Fill((pH1+pH0+pMu0+pMu1).M(),DTFdm());
      //h1_DTFdmtest->Fill(Dst_DTF_Dstarplus_M-Dst_DTF_D0_M);
 }
   
   h1_DTFdm->Write();
   h1_dm->Write();
   h1_DTFmMuMu->Write();
   h1_mMuMu->Write();
   h1_DTFmD->Write();
   h1_mD->Write();
   h2_mD_deltaM->Write();
   //   h1_DTFdmtest->Write();

   delete newfile;

}



 
void D2hhmumuReader::createMCtrainingSample(TString name) {

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

  newtree_even->Branch("Slowpi_cosh",&Slowpi_cosh);
  newtree_even->Branch("D_cosh", & D_cosh);
  newtree_even->Branch("mu0_cosh",&mu0_cosh);
  newtree_even->Branch("deltaM",&deltaM);
  newtree_even->Branch("D_Conemult",&D_Conemult);
  newtree_even->Branch("Dst_Conemult",&Dst_Conemult);
  newtree_even->Branch("D_Coneptasy",&D_Coneptasy);
  newtree_even->Branch("Dst_Coneptasy",&Dst_Coneptasy);
  newtree_even->Branch("mHH", & mHH);

  newtree_odd->Branch("Slowpi_cosh",&Slowpi_cosh);
  newtree_odd->Branch("D_cosh", & D_cosh);
  newtree_odd->Branch("mu0_cosh",&mu0_cosh);
  newtree_odd->Branch("deltaM",&deltaM);
  newtree_odd->Branch("D_Conemult",&D_Conemult);
  newtree_odd->Branch("Dst_Conemult",&Dst_Conemult);
  newtree_odd->Branch("D_Coneptasy",&D_Coneptasy);
  newtree_odd->Branch("Dst_Coneptasy",&Dst_Coneptasy);
  newtree_odd->Branch("mHH", & mHH);


  for (Long64_t i=0;i<nentries; i++) {
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

      if(eventNumber%2==0) newtree_even->Fill(); 
      else newtree_odd->Fill();
   }
   
  std::cout<<"Created MC taining samples "<< name <<" with "<<newtree_even->GetEntries()<<" entries (even sample) and "<< newtree_odd->GetEntries() <<"(odd.)"<<std::endl;

   //newtree->Print();
   newtree_even->AutoSave();
   newtree_odd->AutoSave();
 
   delete newfile;

} 

void D2hhmumuReader::createSubsample(TString name, double percentage) {//specify percentage of full data sample written in sideband subsample

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
     (*it)->Branch("mHH", & mHH);

   }

   for (Long64_t i=0;i<nentries; i++) {
     
     if(i%10000 == 0) std::cout<<i<<" events processed...." << std::endl; 
     
     fChain->GetEntry(i);
      //select only events that pass the trigger seelction criteria 
      if(!isL0Selected() || !isHlt1Selected() || !isHlt2Selected()  )continue;   
      
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
        
      if(eventNumber%2!=0){ //odd events

	trees[0]->Fill(); //odd data
	if(generator->Rndm()<percentage/100 && isBkgSideband()) trees[2]->Fill();
      }
	
      if(eventNumber%2==0){ //even events

	trees[1]->Fill(); //even data
	if(generator->Rndm()<percentage/100 && isBkgSideband()) trees[3]->Fill();
      }
	  
   }
   
   //   std::cout<<"Created random subsample "<< name <<" with "<<newtree->GetEntries()<<" entries "<<std::endl;

 
  //   newtree->AutoSave();
   f_data->cd();
   trees[0]->Write();
   trees[1]->Write();

   //newtree->Print();
   f_sideband->cd();
   trees[2]->Write();
   trees[3]->Write();
  
   delete newfile;

}



void D2hhmumuReader::createValidationSubsample(TString name) {//specify percentage of full data sample written in subsample                                                                   

   Long64_t nentries = fChain->GetEntries();

   activateRelevantBranches();
   fChain->GetEntry(0);

   //Create a new file + a clone of old tree in new file                                                                                                                                      
   TFile *newfile = new TFile(name,"recreate");
   TDirectory *f_data = newfile->mkdir("data");
   TTree *newtree2 = fChain->CloneTree(0);

   double Slowpi_cosh, deltaM,nVeloTracks;
   double D_Conemult,Dst_Conemult,D_Coneptasy,Dst_Coneptasy;
   Int_t nTracks_data,nPVs_data,nSPDHits;

   newtree2->Branch("Slowpi_cosh",&Slowpi_cosh);
   newtree2->Branch("deltaM",&deltaM);
   newtree2->Branch("D_Conemult",&D_Conemult);
   newtree2->Branch("Dst_Conemult",&Dst_Conemult);
   newtree2->Branch("D_Coneptasy",&D_Coneptasy);
   newtree2->Branch("Dst_Coneptasy",&Dst_Coneptasy);
   newtree2->Branch("nSPDHits",&nSPDHits);
   newtree2->Branch("nVeloTracks",&nVeloTracks);
   newtree2->Branch("nTracks_data",&nTracks_data);

   for (Long64_t i=0;i<nentries; i++) {

     if(i%10000 == 0) std::cout<<i<<" events processed...." << std::endl;
     fChain->GetEntry(i);

      //additional variables                                                                                                                                                                  
      Slowpi_cosh=slowpi_helicityAngle();
      deltaM = Dst_DTF_Dstarplus_M - Dst_DTF_D0_M;
      D_Conemult = Dst_CONEMULT_D;
      Dst_Conemult = Dst_CONEMULT_Dstar;
      D_Coneptasy = Dst_CONEPTASYM_D;
      Dst_Coneptasy = Dst_CONEPTASYM_Dstar;
      nSPDHits = (int)nSpdDigits;
      nVeloTracks=nVELO;
      nTracks_data=(int)nTracks;
      nPVs_data=(int)nPVs;

      if( !(Dst_DTF_D0_M>1850 && Dst_DTF_D0_M<1880 && mu0_ProbNNmu>0.5 && mu1_ProbNNmu>0.5 && deltaM >144.5 && deltaM < 146.5 ) ) continue;
      newtree2->Fill();
   }

   std::cout<<"Created subsample "<< name <<" with "<<newtree2->GetEntries()<<" entries "<<std::endl;

   //newtree->Print();                                                                                                                                                                        
   f_data->cd();
   //newtree2->AutoSave();                                                                                                                                                                    
   newtree2->Write();

   delete newfile;

}



void D2hhmumuReader::studyTriggerEfficiency(){

  InitMC(); //only meant to be used for for MC 
 
  std::cout<<"Start evaluating MC trigger efficiencies..."<<  std::cout<<"Done:Histograms are saved in ../rootfiles/triggerEfficienciesMC.root "<<std::endl;

  TFile *newfile = new TFile("../rootFiles/triggerEfficienciesMC.root","recreate");

  TH1D* h1_L0Efficiency= new TH1D("h1_L0Efficiency","Efficiency of L0 selection",4,-.5,3.5);
  TH1D* h1_Hlt1Efficiency= new TH1D("h1_Hlt1Efficiency","Efficiency of Hlt1 selection",7,-.5,6.5);
  TH1D* h1_Hlt2Efficiency= new TH1D("h1_Hlt2Efficiency","Efficiency of Hlt2 selection",4,-.5,3.5);
  TH1D* h1_triggerEfficiency= new TH1D("h1_triggerEfficiency","Efficiency of trigger selection",3,-.5,2.5);

  h1_L0Efficiency->GetXaxis()->SetBinLabel(1,"L0MuonDecision_TOS==1");
  h1_L0Efficiency->GetXaxis()->SetBinLabel(2,"L0DiMuonDecision_TOS==1");
  h1_L0Efficiency->GetXaxis()->SetBinLabel(3,"Ds_L0Global_TIS==1");
  h1_L0Efficiency->GetXaxis()->SetBinLabel(4,"baselineL0");

  h1_Hlt1Efficiency->GetXaxis()->SetBinLabel(1,"Hlt1TrackMuonDecision_TOS==1");
  h1_Hlt1Efficiency->GetXaxis()->SetBinLabel(2,"D_Hlt1TrackAllL0Decision_TOS==1");
  h1_Hlt1Efficiency->GetXaxis()->SetBinLabel(3,"D_Hlt1DiMuonHighMassDecision_TOS ==1");
  h1_Hlt1Efficiency->GetXaxis()->SetBinLabel(4,"D_Hlt1DiMuonLowMassDecision_TOS==1");
  h1_Hlt1Efficiency->GetXaxis()->SetBinLabel(5,"mu0_Hlt1SingleMuonNoIPDecision_TOS==1");
  h1_Hlt1Efficiency->GetXaxis()->SetBinLabel(6,"mu0_Hlt1SingleMuonHighPTDecision_TOS==1");
  h1_Hlt1Efficiency->GetXaxis()->SetBinLabel(7,"baselineHLT1");

  h1_Hlt2Efficiency->GetXaxis()->SetBinLabel(1,"D_Hlt2CharmSemilepD02KKMuMuDecision_TOS==1");
  h1_Hlt2Efficiency->GetXaxis()->SetBinLabel(2,"Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS==1");
  h1_Hlt2Efficiency->GetXaxis()->SetBinLabel(3,"Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS==1");
  h1_Hlt2Efficiency->GetXaxis()->SetBinLabel(4,"logical or");
  
  h1_triggerEfficiency->GetXaxis()->SetBinLabel(1,"L0");
  h1_triggerEfficiency->GetXaxis()->SetBinLabel(2,"L0 && Hlt1");
  h1_triggerEfficiency->GetXaxis()->SetBinLabel(3,"L0 && Hlt1 && Hlt2");


  if (fChain == 0) return;

  
  Long64_t nentries = fChain->GetEntriesFast();
  int nTruthmatched = 0;
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(!MCTruthmatched()) continue;
      nTruthmatched = nTruthmatched +1;

      //Level0 
      if(mu0_L0MuonDecision_TOS==1 || mu1_L0MuonDecision_TOS ==1) h1_L0Efficiency->Fill(0);
      if(mu0_L0DiMuonDecision_TOS==1 || mu1_L0DiMuonDecision_TOS ==1) h1_L0Efficiency->Fill(1);
      if(Dst_L0Global_TIS==1) h1_L0Efficiency->Fill(2);
      if(mu0_L0MuonDecision_TOS==1 || mu1_L0MuonDecision_TOS ==1 || Dst_L0Global_TIS==1 || mu0_L0DiMuonDecision_TOS==1 || mu1_L0DiMuonDecision_TOS ==1 )  h1_L0Efficiency->Fill(3);

      //HLT1  
      if(mu0_Hlt1TrackMuonDecision_TOS==1 || mu1_Hlt1TrackMuonDecision_TOS==1) h1_Hlt1Efficiency->Fill(0);
      if(D_Hlt1TrackAllL0Decision_TOS==1)  h1_Hlt1Efficiency->Fill(1);
      if(D_Hlt1DiMuonHighMassDecision_TOS == 1 )  h1_Hlt1Efficiency->Fill(2);
      if(D_Hlt1DiMuonLowMassDecision_TOS == 1 )  h1_Hlt1Efficiency->Fill(3);
      if( mu0_Hlt1SingleMuonNoIPDecision_TOS == 1 || mu1_Hlt1SingleMuonNoIPDecision_TOS == 1)  h1_Hlt1Efficiency->Fill(4);
      if( mu0_Hlt1SingleMuonHighPTDecision_TOS == 1 || mu1_Hlt1SingleMuonHighPTDecision_TOS == 1 )  h1_Hlt1Efficiency->Fill(5);
      if(mu0_Hlt1TrackMuonDecision_TOS==1 || mu1_Hlt1TrackMuonDecision_TOS==1 || D_Hlt1TrackAllL0Decision_TOS==1 ) h1_Hlt1Efficiency->Fill(6);

      //HLT2
      if(D_Hlt2CharmSemilepD02KKMuMuDecision_TOS==1) h1_Hlt2Efficiency->Fill(0);
      if(Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS==1 /*D_Hlt2DiMuonDetachedDecision_TOS ==1*/) h1_Hlt2Efficiency->Fill(1);
      if(Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS==1) h1_Hlt2Efficiency->Fill(2);
      if(Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS==1 ||Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS==1||D_Hlt2CharmSemilepD02KKMuMuDecision_TOS==1  ) h1_Hlt2Efficiency->Fill(3);


      //total (Standard selection)
      if(isL0Selected()) h1_triggerEfficiency->Fill(0);
      if(isL0Selected() && isHlt1Selected() ) h1_triggerEfficiency->Fill(1);
      if(isL0Selected() && isHlt1Selected() && (isHlt2Selected() || Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS==1 || Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS==1)) h1_triggerEfficiency->Fill(2);

  }

  for(int i =1; i < h1_L0Efficiency->GetNbinsX()+1 ;++i) h1_L0Efficiency->SetBinContent(i,(h1_L0Efficiency->GetBinContent(i)/nTruthmatched));
  for(int i =1; i < h1_Hlt1Efficiency->GetNbinsX()+1 ;++i) h1_Hlt1Efficiency->SetBinContent(i,(h1_Hlt1Efficiency->GetBinContent(i)/nTruthmatched));
  for(int i =1; i < h1_Hlt2Efficiency->GetNbinsX()+1 ;++i) h1_Hlt2Efficiency->SetBinContent(i,(h1_Hlt2Efficiency->GetBinContent(i)/nTruthmatched));
  for(int i =1; i < h1_triggerEfficiency->GetNbinsX()+1 ;++i) h1_triggerEfficiency->SetBinContent(i,(h1_triggerEfficiency->GetBinContent(i)/nTruthmatched));

  TCanvas* canv_trigger = new TCanvas("trigger","trigger");
  canv_trigger->Divide(2,2);
  canv_trigger->cd(1);
  h1_L0Efficiency->Draw();
  canv_trigger->cd(2);
  h1_Hlt2Efficiency->Draw();
  canv_trigger->cd(3);
  h1_Hlt1Efficiency->Draw();
  canv_trigger->cd(4);
  h1_triggerEfficiency->Draw();
  canv_trigger->SaveAs("MCtriggerEfficiency.pdf");

  h1_L0Efficiency->Write();
  h1_Hlt2Efficiency->Write();
  h1_Hlt1Efficiency->Write();
  h1_triggerEfficiency->Write();

  std::cout<<"...total number of truthmatched events :"<<nTruthmatched <<std::endl;
  std::cout<<"After baseline trigger selection:"<<nTruthmatched*h1_triggerEfficiency->GetBinContent(3) <<std::endl;  
  std::cout<<"Efficicncy: "<<h1_triggerEfficiency->GetBinContent(3) <<std::endl;  
  std::cout<<"Done" <<std::endl;
  
  delete newfile;

}  



D2hhmumuReader::D2hhmumuReader(TTree *tree) : fChain(0)  
{
   if (!tree) return;
   Init(tree);
}

D2hhmumuReader::~D2hhmumuReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t D2hhmumuReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t D2hhmumuReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void D2hhmumuReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Dst_MINIP", &Dst_MINIP, &b_Dst_MINIP);
   fChain->SetBranchAddress("Dst_MINIPCHI2", &Dst_MINIPCHI2, &b_Dst_MINIPCHI2);
   fChain->SetBranchAddress("Dst_MINIPNEXTBEST", &Dst_MINIPNEXTBEST, &b_Dst_MINIPNEXTBEST);
   fChain->SetBranchAddress("Dst_MINIPCHI2NEXTBEST", &Dst_MINIPCHI2NEXTBEST, &b_Dst_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("Dst_ENDVERTEX_X", &Dst_ENDVERTEX_X, &b_Dst_ENDVERTEX_X);
   fChain->SetBranchAddress("Dst_ENDVERTEX_Y", &Dst_ENDVERTEX_Y, &b_Dst_ENDVERTEX_Y);
   fChain->SetBranchAddress("Dst_ENDVERTEX_Z", &Dst_ENDVERTEX_Z, &b_Dst_ENDVERTEX_Z);
   fChain->SetBranchAddress("Dst_ENDVERTEX_XERR", &Dst_ENDVERTEX_XERR, &b_Dst_ENDVERTEX_XERR);
   fChain->SetBranchAddress("Dst_ENDVERTEX_YERR", &Dst_ENDVERTEX_YERR, &b_Dst_ENDVERTEX_YERR);
   fChain->SetBranchAddress("Dst_ENDVERTEX_ZERR", &Dst_ENDVERTEX_ZERR, &b_Dst_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("Dst_ENDVERTEX_CHI2", &Dst_ENDVERTEX_CHI2, &b_Dst_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("Dst_ENDVERTEX_NDOF", &Dst_ENDVERTEX_NDOF, &b_Dst_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("Dst_ENDVERTEX_COV_", Dst_ENDVERTEX_COV_, &b_Dst_ENDVERTEX_COV_);
   fChain->SetBranchAddress("Dst_OWNPV_X", &Dst_OWNPV_X, &b_Dst_OWNPV_X);
   fChain->SetBranchAddress("Dst_OWNPV_Y", &Dst_OWNPV_Y, &b_Dst_OWNPV_Y);
   fChain->SetBranchAddress("Dst_OWNPV_Z", &Dst_OWNPV_Z, &b_Dst_OWNPV_Z);
   fChain->SetBranchAddress("Dst_OWNPV_XERR", &Dst_OWNPV_XERR, &b_Dst_OWNPV_XERR);
   fChain->SetBranchAddress("Dst_OWNPV_YERR", &Dst_OWNPV_YERR, &b_Dst_OWNPV_YERR);
   fChain->SetBranchAddress("Dst_OWNPV_ZERR", &Dst_OWNPV_ZERR, &b_Dst_OWNPV_ZERR);
   fChain->SetBranchAddress("Dst_OWNPV_CHI2", &Dst_OWNPV_CHI2, &b_Dst_OWNPV_CHI2);
   fChain->SetBranchAddress("Dst_OWNPV_NDOF", &Dst_OWNPV_NDOF, &b_Dst_OWNPV_NDOF);
   fChain->SetBranchAddress("Dst_OWNPV_COV_", Dst_OWNPV_COV_, &b_Dst_OWNPV_COV_);
   fChain->SetBranchAddress("Dst_IP_OWNPV", &Dst_IP_OWNPV, &b_Dst_IP_OWNPV);
   fChain->SetBranchAddress("Dst_IPCHI2_OWNPV", &Dst_IPCHI2_OWNPV, &b_Dst_IPCHI2_OWNPV);
   fChain->SetBranchAddress("Dst_FD_OWNPV", &Dst_FD_OWNPV, &b_Dst_FD_OWNPV);
   fChain->SetBranchAddress("Dst_FDCHI2_OWNPV", &Dst_FDCHI2_OWNPV, &b_Dst_FDCHI2_OWNPV);
   fChain->SetBranchAddress("Dst_DIRA_OWNPV", &Dst_DIRA_OWNPV, &b_Dst_DIRA_OWNPV);
   fChain->SetBranchAddress("Dst_TOPPV_X", &Dst_TOPPV_X, &b_Dst_TOPPV_X);
   fChain->SetBranchAddress("Dst_TOPPV_Y", &Dst_TOPPV_Y, &b_Dst_TOPPV_Y);
   fChain->SetBranchAddress("Dst_TOPPV_Z", &Dst_TOPPV_Z, &b_Dst_TOPPV_Z);
   fChain->SetBranchAddress("Dst_TOPPV_XERR", &Dst_TOPPV_XERR, &b_Dst_TOPPV_XERR);
   fChain->SetBranchAddress("Dst_TOPPV_YERR", &Dst_TOPPV_YERR, &b_Dst_TOPPV_YERR);
   fChain->SetBranchAddress("Dst_TOPPV_ZERR", &Dst_TOPPV_ZERR, &b_Dst_TOPPV_ZERR);
   fChain->SetBranchAddress("Dst_TOPPV_CHI2", &Dst_TOPPV_CHI2, &b_Dst_TOPPV_CHI2);
   fChain->SetBranchAddress("Dst_TOPPV_NDOF", &Dst_TOPPV_NDOF, &b_Dst_TOPPV_NDOF);
   fChain->SetBranchAddress("Dst_TOPPV_COV_", Dst_TOPPV_COV_, &b_Dst_TOPPV_COV_);
   fChain->SetBranchAddress("Dst_IP_TOPPV", &Dst_IP_TOPPV, &b_Dst_IP_TOPPV);
   fChain->SetBranchAddress("Dst_IPCHI2_TOPPV", &Dst_IPCHI2_TOPPV, &b_Dst_IPCHI2_TOPPV);
   fChain->SetBranchAddress("Dst_FD_TOPPV", &Dst_FD_TOPPV, &b_Dst_FD_TOPPV);
   fChain->SetBranchAddress("Dst_FDCHI2_TOPPV", &Dst_FDCHI2_TOPPV, &b_Dst_FDCHI2_TOPPV);
   fChain->SetBranchAddress("Dst_DIRA_TOPPV", &Dst_DIRA_TOPPV, &b_Dst_DIRA_TOPPV);
   fChain->SetBranchAddress("Dst_P", &Dst_P, &b_Dst_P);
   fChain->SetBranchAddress("Dst_PT", &Dst_PT, &b_Dst_PT);
   fChain->SetBranchAddress("Dst_PE", &Dst_PE, &b_Dst_PE);
   fChain->SetBranchAddress("Dst_PX", &Dst_PX, &b_Dst_PX);
   fChain->SetBranchAddress("Dst_PY", &Dst_PY, &b_Dst_PY);
   fChain->SetBranchAddress("Dst_PZ", &Dst_PZ, &b_Dst_PZ);
   fChain->SetBranchAddress("Dst_MM", &Dst_MM, &b_Dst_MM);
   fChain->SetBranchAddress("Dst_MMERR", &Dst_MMERR, &b_Dst_MMERR);
   fChain->SetBranchAddress("Dst_M", &Dst_M, &b_Dst_M);
   fChain->SetBranchAddress("Dst_ID", &Dst_ID, &b_Dst_ID);
   fChain->SetBranchAddress("Dst_TAU", &Dst_TAU, &b_Dst_TAU);
   fChain->SetBranchAddress("Dst_TAUERR", &Dst_TAUERR, &b_Dst_TAUERR);
   fChain->SetBranchAddress("Dst_TAUCHI2", &Dst_TAUCHI2, &b_Dst_TAUCHI2);
   fChain->SetBranchAddress("Dst_L0Global_Dec", &Dst_L0Global_Dec, &b_Dst_L0Global_Dec);
   fChain->SetBranchAddress("Dst_L0Global_TIS", &Dst_L0Global_TIS, &b_Dst_L0Global_TIS);
   fChain->SetBranchAddress("Dst_L0Global_TOS", &Dst_L0Global_TOS, &b_Dst_L0Global_TOS);
   fChain->SetBranchAddress("Dst_Hlt1Global_Dec", &Dst_Hlt1Global_Dec, &b_Dst_Hlt1Global_Dec);
   fChain->SetBranchAddress("Dst_Hlt1Global_TIS", &Dst_Hlt1Global_TIS, &b_Dst_Hlt1Global_TIS);
   fChain->SetBranchAddress("Dst_Hlt1Global_TOS", &Dst_Hlt1Global_TOS, &b_Dst_Hlt1Global_TOS);
   fChain->SetBranchAddress("Dst_Hlt1Phys_Dec", &Dst_Hlt1Phys_Dec, &b_Dst_Hlt1Phys_Dec);
   fChain->SetBranchAddress("Dst_Hlt1Phys_TIS", &Dst_Hlt1Phys_TIS, &b_Dst_Hlt1Phys_TIS);
   fChain->SetBranchAddress("Dst_Hlt1Phys_TOS", &Dst_Hlt1Phys_TOS, &b_Dst_Hlt1Phys_TOS);
   fChain->SetBranchAddress("Dst_Hlt2Global_Dec", &Dst_Hlt2Global_Dec, &b_Dst_Hlt2Global_Dec);
   fChain->SetBranchAddress("Dst_Hlt2Global_TIS", &Dst_Hlt2Global_TIS, &b_Dst_Hlt2Global_TIS);
   fChain->SetBranchAddress("Dst_Hlt2Global_TOS", &Dst_Hlt2Global_TOS, &b_Dst_Hlt2Global_TOS);
   fChain->SetBranchAddress("Dst_Hlt2Phys_Dec", &Dst_Hlt2Phys_Dec, &b_Dst_Hlt2Phys_Dec);
   fChain->SetBranchAddress("Dst_Hlt2Phys_TIS", &Dst_Hlt2Phys_TIS, &b_Dst_Hlt2Phys_TIS);
   fChain->SetBranchAddress("Dst_Hlt2Phys_TOS", &Dst_Hlt2Phys_TOS, &b_Dst_Hlt2Phys_TOS);
   fChain->SetBranchAddress("Dst_DTF_D0_E", &Dst_DTF_D0_E, &b_Dst_DTF_D0_E);
   fChain->SetBranchAddress("Dst_DTF_D0_M", &Dst_DTF_D0_M, &b_Dst_DTF_D0_M);
   fChain->SetBranchAddress("Dst_DTF_D0_P", &Dst_DTF_D0_P, &b_Dst_DTF_D0_P);
   fChain->SetBranchAddress("Dst_DTF_D0_PT", &Dst_DTF_D0_PT, &b_Dst_DTF_D0_PT);
   fChain->SetBranchAddress("Dst_DTF_D0_PX", &Dst_DTF_D0_PX, &b_Dst_DTF_D0_PX);
   fChain->SetBranchAddress("Dst_DTF_D0_PY", &Dst_DTF_D0_PY, &b_Dst_DTF_D0_PY);
   fChain->SetBranchAddress("Dst_DTF_D0_PZ", &Dst_DTF_D0_PZ, &b_Dst_DTF_D0_PZ);
   fChain->SetBranchAddress("Dst_DTF_Dstarplus_E", &Dst_DTF_Dstarplus_E, &b_Dst_DTF_Dstarplus_E);
   fChain->SetBranchAddress("Dst_DTF_Dstarplus_M", &Dst_DTF_Dstarplus_M, &b_Dst_DTF_Dstarplus_M);
   fChain->SetBranchAddress("Dst_DTF_Dstarplus_P", &Dst_DTF_Dstarplus_P, &b_Dst_DTF_Dstarplus_P);
   fChain->SetBranchAddress("Dst_DTF_Dstarplus_PT", &Dst_DTF_Dstarplus_PT, &b_Dst_DTF_Dstarplus_PT);
   fChain->SetBranchAddress("Dst_DTF_Dstarplus_PX", &Dst_DTF_Dstarplus_PX, &b_Dst_DTF_Dstarplus_PX);
   fChain->SetBranchAddress("Dst_DTF_Dstarplus_PY", &Dst_DTF_Dstarplus_PY, &b_Dst_DTF_Dstarplus_PY);
   fChain->SetBranchAddress("Dst_DTF_Dstarplus_PZ", &Dst_DTF_Dstarplus_PZ, &b_Dst_DTF_Dstarplus_PZ);
   fChain->SetBranchAddress("Dst_DTF_NDOF", &Dst_DTF_NDOF, &b_Dst_DTF_NDOF);
   fChain->SetBranchAddress("Dst_DTF_Pis_BPVIPCHI2", &Dst_DTF_Pis_BPVIPCHI2, &b_Dst_DTF_Pis_BPVIPCHI2);
   fChain->SetBranchAddress("Dst_DTF_Pis_E", &Dst_DTF_Pis_E, &b_Dst_DTF_Pis_E);
   fChain->SetBranchAddress("Dst_DTF_Pis_M", &Dst_DTF_Pis_M, &b_Dst_DTF_Pis_M);
   fChain->SetBranchAddress("Dst_DTF_Pis_P", &Dst_DTF_Pis_P, &b_Dst_DTF_Pis_P);
   fChain->SetBranchAddress("Dst_DTF_Pis_PT", &Dst_DTF_Pis_PT, &b_Dst_DTF_Pis_PT);
   fChain->SetBranchAddress("Dst_DTF_Pis_PX", &Dst_DTF_Pis_PX, &b_Dst_DTF_Pis_PX);
   fChain->SetBranchAddress("Dst_DTF_Pis_PY", &Dst_DTF_Pis_PY, &b_Dst_DTF_Pis_PY);
   fChain->SetBranchAddress("Dst_DTF_Pis_PZ", &Dst_DTF_Pis_PZ, &b_Dst_DTF_Pis_PZ);
   fChain->SetBranchAddress("Dst_DTF_h0_E", &Dst_DTF_h0_E, &b_Dst_DTF_h0_E);
   fChain->SetBranchAddress("Dst_DTF_h0_P", &Dst_DTF_h0_P, &b_Dst_DTF_h0_P);
   fChain->SetBranchAddress("Dst_DTF_h0_PT", &Dst_DTF_h0_PT, &b_Dst_DTF_h0_PT);
   fChain->SetBranchAddress("Dst_DTF_h0_PX", &Dst_DTF_h0_PX, &b_Dst_DTF_h0_PX);
   fChain->SetBranchAddress("Dst_DTF_h0_PY", &Dst_DTF_h0_PY, &b_Dst_DTF_h0_PY);
   fChain->SetBranchAddress("Dst_DTF_h0_PZ", &Dst_DTF_h0_PZ, &b_Dst_DTF_h0_PZ);
   fChain->SetBranchAddress("Dst_DTF_h1_E", &Dst_DTF_h1_E, &b_Dst_DTF_h1_E);
   fChain->SetBranchAddress("Dst_DTF_h1_P", &Dst_DTF_h1_P, &b_Dst_DTF_h1_P);
   fChain->SetBranchAddress("Dst_DTF_h1_PT", &Dst_DTF_h1_PT, &b_Dst_DTF_h1_PT);
   fChain->SetBranchAddress("Dst_DTF_h1_PX", &Dst_DTF_h1_PX, &b_Dst_DTF_h1_PX);
   fChain->SetBranchAddress("Dst_DTF_h1_PY", &Dst_DTF_h1_PY, &b_Dst_DTF_h1_PY);
   fChain->SetBranchAddress("Dst_DTF_h1_PZ", &Dst_DTF_h1_PZ, &b_Dst_DTF_h1_PZ);
   fChain->SetBranchAddress("Dst_DTF_mu0_E", &Dst_DTF_mu0_E, &b_Dst_DTF_mu0_E);
   fChain->SetBranchAddress("Dst_DTF_mu0_P", &Dst_DTF_mu0_P, &b_Dst_DTF_mu0_P);
   fChain->SetBranchAddress("Dst_DTF_mu0_PT", &Dst_DTF_mu0_PT, &b_Dst_DTF_mu0_PT);
   fChain->SetBranchAddress("Dst_DTF_mu0_PX", &Dst_DTF_mu0_PX, &b_Dst_DTF_mu0_PX);
   fChain->SetBranchAddress("Dst_DTF_mu0_PY", &Dst_DTF_mu0_PY, &b_Dst_DTF_mu0_PY);
   fChain->SetBranchAddress("Dst_DTF_mu0_PZ", &Dst_DTF_mu0_PZ, &b_Dst_DTF_mu0_PZ);
   fChain->SetBranchAddress("Dst_DTF_mu1_E", &Dst_DTF_mu1_E, &b_Dst_DTF_mu1_E);
   fChain->SetBranchAddress("Dst_DTF_mu1_P", &Dst_DTF_mu1_P, &b_Dst_DTF_mu1_P);
   fChain->SetBranchAddress("Dst_DTF_mu1_PT", &Dst_DTF_mu1_PT, &b_Dst_DTF_mu1_PT);
   fChain->SetBranchAddress("Dst_DTF_mu1_PX", &Dst_DTF_mu1_PX, &b_Dst_DTF_mu1_PX);
   fChain->SetBranchAddress("Dst_DTF_mu1_PY", &Dst_DTF_mu1_PY, &b_Dst_DTF_mu1_PY);
   fChain->SetBranchAddress("Dst_DTF_mu1_PZ", &Dst_DTF_mu1_PZ, &b_Dst_DTF_mu1_PZ);
   fChain->SetBranchAddress("Dst_CONEANGLE_D", &Dst_CONEANGLE_D, &b_Dst_CONEANGLE_D);
   fChain->SetBranchAddress("Dst_CONEANGLE_Dstar", &Dst_CONEANGLE_Dstar, &b_Dst_CONEANGLE_Dstar);
   fChain->SetBranchAddress("Dst_CONEMULT_D", &Dst_CONEMULT_D, &b_Dst_CONEMULT_D);
   fChain->SetBranchAddress("Dst_CONEMULT_Dstar", &Dst_CONEMULT_Dstar, &b_Dst_CONEMULT_Dstar);
   fChain->SetBranchAddress("Dst_CONEPTASYM_D", &Dst_CONEPTASYM_D, &b_Dst_CONEPTASYM_D);
   fChain->SetBranchAddress("Dst_CONEPTASYM_Dstar", &Dst_CONEPTASYM_Dstar, &b_Dst_CONEPTASYM_Dstar);
   fChain->SetBranchAddress("Dst_MAXDOCA", &Dst_MAXDOCA, &b_Dst_MAXDOCA);
   fChain->SetBranchAddress("Dst_L0HadronDecision_Dec", &Dst_L0HadronDecision_Dec, &b_Dst_L0HadronDecision_Dec);
   fChain->SetBranchAddress("Dst_L0HadronDecision_TIS", &Dst_L0HadronDecision_TIS, &b_Dst_L0HadronDecision_TIS);
   fChain->SetBranchAddress("Dst_L0HadronDecision_TOS", &Dst_L0HadronDecision_TOS, &b_Dst_L0HadronDecision_TOS);
   fChain->SetBranchAddress("Dst_L0MuonDecision_Dec", &Dst_L0MuonDecision_Dec, &b_Dst_L0MuonDecision_Dec);
   fChain->SetBranchAddress("Dst_L0MuonDecision_TIS", &Dst_L0MuonDecision_TIS, &b_Dst_L0MuonDecision_TIS);
   fChain->SetBranchAddress("Dst_L0MuonDecision_TOS", &Dst_L0MuonDecision_TOS, &b_Dst_L0MuonDecision_TOS);
   fChain->SetBranchAddress("Dst_L0DiMuonDecision_Dec", &Dst_L0DiMuonDecision_Dec, &b_Dst_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("Dst_L0DiMuonDecision_TIS", &Dst_L0DiMuonDecision_TIS, &b_Dst_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("Dst_L0DiMuonDecision_TOS", &Dst_L0DiMuonDecision_TOS, &b_Dst_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("Dst_L0ElectronDecision_Dec", &Dst_L0ElectronDecision_Dec, &b_Dst_L0ElectronDecision_Dec);
   fChain->SetBranchAddress("Dst_L0ElectronDecision_TIS", &Dst_L0ElectronDecision_TIS, &b_Dst_L0ElectronDecision_TIS);
   fChain->SetBranchAddress("Dst_L0ElectronDecision_TOS", &Dst_L0ElectronDecision_TOS, &b_Dst_L0ElectronDecision_TOS);
   fChain->SetBranchAddress("Dst_L0PhotonDecision_Dec", &Dst_L0PhotonDecision_Dec, &b_Dst_L0PhotonDecision_Dec);
   fChain->SetBranchAddress("Dst_L0PhotonDecision_TIS", &Dst_L0PhotonDecision_TIS, &b_Dst_L0PhotonDecision_TIS);
   fChain->SetBranchAddress("Dst_L0PhotonDecision_TOS", &Dst_L0PhotonDecision_TOS, &b_Dst_L0PhotonDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt1DiMuonHighMassDecision_Dec", &Dst_Hlt1DiMuonHighMassDecision_Dec, &b_Dst_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt1DiMuonHighMassDecision_TIS", &Dst_Hlt1DiMuonHighMassDecision_TIS, &b_Dst_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt1DiMuonHighMassDecision_TOS", &Dst_Hlt1DiMuonHighMassDecision_TOS, &b_Dst_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt1DiMuonLowMassDecision_Dec", &Dst_Hlt1DiMuonLowMassDecision_Dec, &b_Dst_Hlt1DiMuonLowMassDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt1DiMuonLowMassDecision_TIS", &Dst_Hlt1DiMuonLowMassDecision_TIS, &b_Dst_Hlt1DiMuonLowMassDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt1DiMuonLowMassDecision_TOS", &Dst_Hlt1DiMuonLowMassDecision_TOS, &b_Dst_Hlt1DiMuonLowMassDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt1SingleMuonNoIPDecision_Dec", &Dst_Hlt1SingleMuonNoIPDecision_Dec, &b_Dst_Hlt1SingleMuonNoIPDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt1SingleMuonNoIPDecision_TIS", &Dst_Hlt1SingleMuonNoIPDecision_TIS, &b_Dst_Hlt1SingleMuonNoIPDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt1SingleMuonNoIPDecision_TOS", &Dst_Hlt1SingleMuonNoIPDecision_TOS, &b_Dst_Hlt1SingleMuonNoIPDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt1SingleMuonHighPTDecision_Dec", &Dst_Hlt1SingleMuonHighPTDecision_Dec, &b_Dst_Hlt1SingleMuonHighPTDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt1SingleMuonHighPTDecision_TIS", &Dst_Hlt1SingleMuonHighPTDecision_TIS, &b_Dst_Hlt1SingleMuonHighPTDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt1SingleMuonHighPTDecision_TOS", &Dst_Hlt1SingleMuonHighPTDecision_TOS, &b_Dst_Hlt1SingleMuonHighPTDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt1TrackAllL0Decision_Dec", &Dst_Hlt1TrackAllL0Decision_Dec, &b_Dst_Hlt1TrackAllL0Decision_Dec);
   fChain->SetBranchAddress("Dst_Hlt1TrackAllL0Decision_TIS", &Dst_Hlt1TrackAllL0Decision_TIS, &b_Dst_Hlt1TrackAllL0Decision_TIS);
   fChain->SetBranchAddress("Dst_Hlt1TrackAllL0Decision_TOS", &Dst_Hlt1TrackAllL0Decision_TOS, &b_Dst_Hlt1TrackAllL0Decision_TOS);
   fChain->SetBranchAddress("Dst_Hlt1TrackMuonDecision_Dec", &Dst_Hlt1TrackMuonDecision_Dec, &b_Dst_Hlt1TrackMuonDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt1TrackMuonDecision_TIS", &Dst_Hlt1TrackMuonDecision_TIS, &b_Dst_Hlt1TrackMuonDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt1TrackMuonDecision_TOS", &Dst_Hlt1TrackMuonDecision_TOS, &b_Dst_Hlt1TrackMuonDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt1TrackPhotonDecision_Dec", &Dst_Hlt1TrackPhotonDecision_Dec, &b_Dst_Hlt1TrackPhotonDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt1TrackPhotonDecision_TIS", &Dst_Hlt1TrackPhotonDecision_TIS, &b_Dst_Hlt1TrackPhotonDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt1TrackPhotonDecision_TOS", &Dst_Hlt1TrackPhotonDecision_TOS, &b_Dst_Hlt1TrackPhotonDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt1L0AnyDecision_Dec", &Dst_Hlt1L0AnyDecision_Dec, &b_Dst_Hlt1L0AnyDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt1L0AnyDecision_TIS", &Dst_Hlt1L0AnyDecision_TIS, &b_Dst_Hlt1L0AnyDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt1L0AnyDecision_TOS", &Dst_Hlt1L0AnyDecision_TOS, &b_Dst_Hlt1L0AnyDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt1GlobalDecision_Dec", &Dst_Hlt1GlobalDecision_Dec, &b_Dst_Hlt1GlobalDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt1GlobalDecision_TIS", &Dst_Hlt1GlobalDecision_TIS, &b_Dst_Hlt1GlobalDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt1GlobalDecision_TOS", &Dst_Hlt1GlobalDecision_TOS, &b_Dst_Hlt1GlobalDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2SingleMuonDecision_Dec", &Dst_Hlt2SingleMuonDecision_Dec, &b_Dst_Hlt2SingleMuonDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2SingleMuonDecision_TIS", &Dst_Hlt2SingleMuonDecision_TIS, &b_Dst_Hlt2SingleMuonDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2SingleMuonDecision_TOS", &Dst_Hlt2SingleMuonDecision_TOS, &b_Dst_Hlt2SingleMuonDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2DiMuonDetachedDecision_Dec", &Dst_Hlt2DiMuonDetachedDecision_Dec, &b_Dst_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2DiMuonDetachedDecision_TIS", &Dst_Hlt2DiMuonDetachedDecision_TIS, &b_Dst_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2DiMuonDetachedDecision_TOS", &Dst_Hlt2DiMuonDetachedDecision_TOS, &b_Dst_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD2HMuMuDecision_Dec", &Dst_Hlt2CharmSemilepD2HMuMuDecision_Dec, &b_Dst_Hlt2CharmSemilepD2HMuMuDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD2HMuMuDecision_TIS", &Dst_Hlt2CharmSemilepD2HMuMuDecision_TIS, &b_Dst_Hlt2CharmSemilepD2HMuMuDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD2HMuMuDecision_TOS", &Dst_Hlt2CharmSemilepD2HMuMuDecision_TOS, &b_Dst_Hlt2CharmSemilepD2HMuMuDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec", &Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec, &b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS", &Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS, &b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS", &Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS, &b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec", &Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec, &b_Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS", &Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS, &b_Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS", &Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS, &b_Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec", &Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec, &b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS", &Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS, &b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS", &Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS, &b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec", &Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec, &b_Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS", &Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS, &b_Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS", &Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS, &b_Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec", &Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec, &b_Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS", &Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS, &b_Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS", &Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS, &b_Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD02KKMuMuDecision_Dec", &Dst_Hlt2CharmSemilepD02KKMuMuDecision_Dec, &b_Dst_Hlt2CharmSemilepD02KKMuMuDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD02KKMuMuDecision_TIS", &Dst_Hlt2CharmSemilepD02KKMuMuDecision_TIS, &b_Dst_Hlt2CharmSemilepD02KKMuMuDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD02KKMuMuDecision_TOS", &Dst_Hlt2CharmSemilepD02KKMuMuDecision_TOS, &b_Dst_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD02KPiMuMuDecision_Dec", &Dst_Hlt2CharmSemilepD02KPiMuMuDecision_Dec, &b_Dst_Hlt2CharmSemilepD02KPiMuMuDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD02KPiMuMuDecision_TIS", &Dst_Hlt2CharmSemilepD02KPiMuMuDecision_TIS, &b_Dst_Hlt2CharmSemilepD02KPiMuMuDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmSemilepD02KPiMuMuDecision_TOS", &Dst_Hlt2CharmSemilepD02KPiMuMuDecision_TOS, &b_Dst_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_Dec", &Dst_Hlt2CharmHadD02HHHHDst_4piDecision_Dec, &b_Dst_Hlt2CharmHadD02HHHHDst_4piDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TIS", &Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TIS, &b_Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS", &Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS, &b_Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec", &Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec, &b_Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS", &Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS, &b_Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS", &Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS, &b_Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec", &Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec, &b_Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS", &Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS, &b_Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS", &Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS, &b_Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_Dec", &Dst_Hlt2CharmHadD02HHXDst_hhXDecision_Dec, &b_Dst_Hlt2CharmHadD02HHXDst_hhXDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TIS", &Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TIS, &b_Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS", &Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS, &b_Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec", &Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec, &b_Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS", &Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS, &b_Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS", &Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS, &b_Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec", &Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec, &b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS", &Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS, &b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS", &Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS, &b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec", &Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec, &b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS", &Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS, &b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS", &Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS, &b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec", &Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec, &b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS", &Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS, &b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS", &Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS, &b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec", &Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec, &b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS", &Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS, &b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS", &Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS, &b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_K3piDecision_Dec", &Dst_Hlt2CharmHadD02HHHH_K3piDecision_Dec, &b_Dst_Hlt2CharmHadD02HHHH_K3piDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_K3piDecision_TIS", &Dst_Hlt2CharmHadD02HHHH_K3piDecision_TIS, &b_Dst_Hlt2CharmHadD02HHHH_K3piDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_K3piDecision_TOS", &Dst_Hlt2CharmHadD02HHHH_K3piDecision_TOS, &b_Dst_Hlt2CharmHadD02HHHH_K3piDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec", &Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec, &b_Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS", &Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS, &b_Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS", &Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS, &b_Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec", &Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec, &b_Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS", &Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS, &b_Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS", &Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS, &b_Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec", &Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec, &b_Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS", &Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS, &b_Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS", &Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS, &b_Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_4piDecision_Dec", &Dst_Hlt2CharmHadD02HHHH_4piDecision_Dec, &b_Dst_Hlt2CharmHadD02HHHH_4piDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_4piDecision_TIS", &Dst_Hlt2CharmHadD02HHHH_4piDecision_TIS, &b_Dst_Hlt2CharmHadD02HHHH_4piDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_4piDecision_TOS", &Dst_Hlt2CharmHadD02HHHH_4piDecision_TOS, &b_Dst_Hlt2CharmHadD02HHHH_4piDecision_TOS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec", &Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec, &b_Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS", &Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS, &b_Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS);
   fChain->SetBranchAddress("Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS", &Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS, &b_Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS);
   fChain->SetBranchAddress("Dst_cpx_1.00", &Dst_cpx_1_00, &b_Dst_cpx_1_00);
   fChain->SetBranchAddress("Dst_cpy_1.00", &Dst_cpy_1_00, &b_Dst_cpy_1_00);
   fChain->SetBranchAddress("Dst_cpz_1.00", &Dst_cpz_1_00, &b_Dst_cpz_1_00);
   fChain->SetBranchAddress("Dst_cpt_1.00", &Dst_cpt_1_00, &b_Dst_cpt_1_00);
   fChain->SetBranchAddress("Dst_cp_1.00", &Dst_cp_1_00, &b_Dst_cp_1_00);
   fChain->SetBranchAddress("Dst_cmult_1.00", &Dst_cmult_1_00, &b_Dst_cmult_1_00);
   fChain->SetBranchAddress("Dst_deltaEta_1.00", &Dst_deltaEta_1_00, &b_Dst_deltaEta_1_00);
   fChain->SetBranchAddress("Dst_deltaPhi_1.00", &Dst_deltaPhi_1_00, &b_Dst_deltaPhi_1_00);
   fChain->SetBranchAddress("Dst_pxasy_1.00", &Dst_pxasy_1_00, &b_Dst_pxasy_1_00);
   fChain->SetBranchAddress("Dst_pyasy_1.00", &Dst_pyasy_1_00, &b_Dst_pyasy_1_00);
   fChain->SetBranchAddress("Dst_pzasy_1.00", &Dst_pzasy_1_00, &b_Dst_pzasy_1_00);
   fChain->SetBranchAddress("Dst_pasy_1.00", &Dst_pasy_1_00, &b_Dst_pasy_1_00);
   fChain->SetBranchAddress("Dst_ptasy_1.00", &Dst_ptasy_1_00, &b_Dst_ptasy_1_00);
   fChain->SetBranchAddress("Dst_cpx_1.10", &Dst_cpx_1_10, &b_Dst_cpx_1_10);
   fChain->SetBranchAddress("Dst_cpy_1.10", &Dst_cpy_1_10, &b_Dst_cpy_1_10);
   fChain->SetBranchAddress("Dst_cpz_1.10", &Dst_cpz_1_10, &b_Dst_cpz_1_10);
   fChain->SetBranchAddress("Dst_cpt_1.10", &Dst_cpt_1_10, &b_Dst_cpt_1_10);
   fChain->SetBranchAddress("Dst_cp_1.10", &Dst_cp_1_10, &b_Dst_cp_1_10);
   fChain->SetBranchAddress("Dst_cmult_1.10", &Dst_cmult_1_10, &b_Dst_cmult_1_10);
   fChain->SetBranchAddress("Dst_deltaEta_1.10", &Dst_deltaEta_1_10, &b_Dst_deltaEta_1_10);
   fChain->SetBranchAddress("Dst_deltaPhi_1.10", &Dst_deltaPhi_1_10, &b_Dst_deltaPhi_1_10);
   fChain->SetBranchAddress("Dst_pxasy_1.10", &Dst_pxasy_1_10, &b_Dst_pxasy_1_10);
   fChain->SetBranchAddress("Dst_pyasy_1.10", &Dst_pyasy_1_10, &b_Dst_pyasy_1_10);
   fChain->SetBranchAddress("Dst_pzasy_1.10", &Dst_pzasy_1_10, &b_Dst_pzasy_1_10);
   fChain->SetBranchAddress("Dst_pasy_1.10", &Dst_pasy_1_10, &b_Dst_pasy_1_10);
   fChain->SetBranchAddress("Dst_ptasy_1.10", &Dst_ptasy_1_10, &b_Dst_ptasy_1_10);
   fChain->SetBranchAddress("Dst_cpx_1.20", &Dst_cpx_1_20, &b_Dst_cpx_1_20);
   fChain->SetBranchAddress("Dst_cpy_1.20", &Dst_cpy_1_20, &b_Dst_cpy_1_20);
   fChain->SetBranchAddress("Dst_cpz_1.20", &Dst_cpz_1_20, &b_Dst_cpz_1_20);
   fChain->SetBranchAddress("Dst_cpt_1.20", &Dst_cpt_1_20, &b_Dst_cpt_1_20);
   fChain->SetBranchAddress("Dst_cp_1.20", &Dst_cp_1_20, &b_Dst_cp_1_20);
   fChain->SetBranchAddress("Dst_cmult_1.20", &Dst_cmult_1_20, &b_Dst_cmult_1_20);
   fChain->SetBranchAddress("Dst_deltaEta_1.20", &Dst_deltaEta_1_20, &b_Dst_deltaEta_1_20);
   fChain->SetBranchAddress("Dst_deltaPhi_1.20", &Dst_deltaPhi_1_20, &b_Dst_deltaPhi_1_20);
   fChain->SetBranchAddress("Dst_pxasy_1.20", &Dst_pxasy_1_20, &b_Dst_pxasy_1_20);
   fChain->SetBranchAddress("Dst_pyasy_1.20", &Dst_pyasy_1_20, &b_Dst_pyasy_1_20);
   fChain->SetBranchAddress("Dst_pzasy_1.20", &Dst_pzasy_1_20, &b_Dst_pzasy_1_20);
   fChain->SetBranchAddress("Dst_pasy_1.20", &Dst_pasy_1_20, &b_Dst_pasy_1_20);
   fChain->SetBranchAddress("Dst_ptasy_1.20", &Dst_ptasy_1_20, &b_Dst_ptasy_1_20);
   fChain->SetBranchAddress("Dst_cpx_1.30", &Dst_cpx_1_30, &b_Dst_cpx_1_30);
   fChain->SetBranchAddress("Dst_cpy_1.30", &Dst_cpy_1_30, &b_Dst_cpy_1_30);
   fChain->SetBranchAddress("Dst_cpz_1.30", &Dst_cpz_1_30, &b_Dst_cpz_1_30);
   fChain->SetBranchAddress("Dst_cpt_1.30", &Dst_cpt_1_30, &b_Dst_cpt_1_30);
   fChain->SetBranchAddress("Dst_cp_1.30", &Dst_cp_1_30, &b_Dst_cp_1_30);
   fChain->SetBranchAddress("Dst_cmult_1.30", &Dst_cmult_1_30, &b_Dst_cmult_1_30);
   fChain->SetBranchAddress("Dst_deltaEta_1.30", &Dst_deltaEta_1_30, &b_Dst_deltaEta_1_30);
   fChain->SetBranchAddress("Dst_deltaPhi_1.30", &Dst_deltaPhi_1_30, &b_Dst_deltaPhi_1_30);
   fChain->SetBranchAddress("Dst_pxasy_1.30", &Dst_pxasy_1_30, &b_Dst_pxasy_1_30);
   fChain->SetBranchAddress("Dst_pyasy_1.30", &Dst_pyasy_1_30, &b_Dst_pyasy_1_30);
   fChain->SetBranchAddress("Dst_pzasy_1.30", &Dst_pzasy_1_30, &b_Dst_pzasy_1_30);
   fChain->SetBranchAddress("Dst_pasy_1.30", &Dst_pasy_1_30, &b_Dst_pasy_1_30);
   fChain->SetBranchAddress("Dst_ptasy_1.30", &Dst_ptasy_1_30, &b_Dst_ptasy_1_30);
   fChain->SetBranchAddress("Dst_cpx_1.40", &Dst_cpx_1_40, &b_Dst_cpx_1_40);
   fChain->SetBranchAddress("Dst_cpy_1.40", &Dst_cpy_1_40, &b_Dst_cpy_1_40);
   fChain->SetBranchAddress("Dst_cpz_1.40", &Dst_cpz_1_40, &b_Dst_cpz_1_40);
   fChain->SetBranchAddress("Dst_cpt_1.40", &Dst_cpt_1_40, &b_Dst_cpt_1_40);
   fChain->SetBranchAddress("Dst_cp_1.40", &Dst_cp_1_40, &b_Dst_cp_1_40);
   fChain->SetBranchAddress("Dst_cmult_1.40", &Dst_cmult_1_40, &b_Dst_cmult_1_40);
   fChain->SetBranchAddress("Dst_deltaEta_1.40", &Dst_deltaEta_1_40, &b_Dst_deltaEta_1_40);
   fChain->SetBranchAddress("Dst_deltaPhi_1.40", &Dst_deltaPhi_1_40, &b_Dst_deltaPhi_1_40);
   fChain->SetBranchAddress("Dst_pxasy_1.40", &Dst_pxasy_1_40, &b_Dst_pxasy_1_40);
   fChain->SetBranchAddress("Dst_pyasy_1.40", &Dst_pyasy_1_40, &b_Dst_pyasy_1_40);
   fChain->SetBranchAddress("Dst_pzasy_1.40", &Dst_pzasy_1_40, &b_Dst_pzasy_1_40);
   fChain->SetBranchAddress("Dst_pasy_1.40", &Dst_pasy_1_40, &b_Dst_pasy_1_40);
   fChain->SetBranchAddress("Dst_ptasy_1.40", &Dst_ptasy_1_40, &b_Dst_ptasy_1_40);
   fChain->SetBranchAddress("Dst_cpx_1.50", &Dst_cpx_1_50, &b_Dst_cpx_1_50);
   fChain->SetBranchAddress("Dst_cpy_1.50", &Dst_cpy_1_50, &b_Dst_cpy_1_50);
   fChain->SetBranchAddress("Dst_cpz_1.50", &Dst_cpz_1_50, &b_Dst_cpz_1_50);
   fChain->SetBranchAddress("Dst_cpt_1.50", &Dst_cpt_1_50, &b_Dst_cpt_1_50);
   fChain->SetBranchAddress("Dst_cp_1.50", &Dst_cp_1_50, &b_Dst_cp_1_50);
   fChain->SetBranchAddress("Dst_cmult_1.50", &Dst_cmult_1_50, &b_Dst_cmult_1_50);
   fChain->SetBranchAddress("Dst_deltaEta_1.50", &Dst_deltaEta_1_50, &b_Dst_deltaEta_1_50);
   fChain->SetBranchAddress("Dst_deltaPhi_1.50", &Dst_deltaPhi_1_50, &b_Dst_deltaPhi_1_50);
   fChain->SetBranchAddress("Dst_pxasy_1.50", &Dst_pxasy_1_50, &b_Dst_pxasy_1_50);
   fChain->SetBranchAddress("Dst_pyasy_1.50", &Dst_pyasy_1_50, &b_Dst_pyasy_1_50);
   fChain->SetBranchAddress("Dst_pzasy_1.50", &Dst_pzasy_1_50, &b_Dst_pzasy_1_50);
   fChain->SetBranchAddress("Dst_pasy_1.50", &Dst_pasy_1_50, &b_Dst_pasy_1_50);
   fChain->SetBranchAddress("Dst_ptasy_1.50", &Dst_ptasy_1_50, &b_Dst_ptasy_1_50);
   fChain->SetBranchAddress("Dst_cpx_1.60", &Dst_cpx_1_60, &b_Dst_cpx_1_60);
   fChain->SetBranchAddress("Dst_cpy_1.60", &Dst_cpy_1_60, &b_Dst_cpy_1_60);
   fChain->SetBranchAddress("Dst_cpz_1.60", &Dst_cpz_1_60, &b_Dst_cpz_1_60);
   fChain->SetBranchAddress("Dst_cpt_1.60", &Dst_cpt_1_60, &b_Dst_cpt_1_60);
   fChain->SetBranchAddress("Dst_cp_1.60", &Dst_cp_1_60, &b_Dst_cp_1_60);
   fChain->SetBranchAddress("Dst_cmult_1.60", &Dst_cmult_1_60, &b_Dst_cmult_1_60);
   fChain->SetBranchAddress("Dst_deltaEta_1.60", &Dst_deltaEta_1_60, &b_Dst_deltaEta_1_60);
   fChain->SetBranchAddress("Dst_deltaPhi_1.60", &Dst_deltaPhi_1_60, &b_Dst_deltaPhi_1_60);
   fChain->SetBranchAddress("Dst_pxasy_1.60", &Dst_pxasy_1_60, &b_Dst_pxasy_1_60);
   fChain->SetBranchAddress("Dst_pyasy_1.60", &Dst_pyasy_1_60, &b_Dst_pyasy_1_60);
   fChain->SetBranchAddress("Dst_pzasy_1.60", &Dst_pzasy_1_60, &b_Dst_pzasy_1_60);
   fChain->SetBranchAddress("Dst_pasy_1.60", &Dst_pasy_1_60, &b_Dst_pasy_1_60);
   fChain->SetBranchAddress("Dst_ptasy_1.60", &Dst_ptasy_1_60, &b_Dst_ptasy_1_60);
   fChain->SetBranchAddress("Dst_cpx_1.70", &Dst_cpx_1_70, &b_Dst_cpx_1_70);
   fChain->SetBranchAddress("Dst_cpy_1.70", &Dst_cpy_1_70, &b_Dst_cpy_1_70);
   fChain->SetBranchAddress("Dst_cpz_1.70", &Dst_cpz_1_70, &b_Dst_cpz_1_70);
   fChain->SetBranchAddress("Dst_cpt_1.70", &Dst_cpt_1_70, &b_Dst_cpt_1_70);
   fChain->SetBranchAddress("Dst_cp_1.70", &Dst_cp_1_70, &b_Dst_cp_1_70);
   fChain->SetBranchAddress("Dst_cmult_1.70", &Dst_cmult_1_70, &b_Dst_cmult_1_70);
   fChain->SetBranchAddress("Dst_deltaEta_1.70", &Dst_deltaEta_1_70, &b_Dst_deltaEta_1_70);
   fChain->SetBranchAddress("Dst_deltaPhi_1.70", &Dst_deltaPhi_1_70, &b_Dst_deltaPhi_1_70);
   fChain->SetBranchAddress("Dst_pxasy_1.70", &Dst_pxasy_1_70, &b_Dst_pxasy_1_70);
   fChain->SetBranchAddress("Dst_pyasy_1.70", &Dst_pyasy_1_70, &b_Dst_pyasy_1_70);
   fChain->SetBranchAddress("Dst_pzasy_1.70", &Dst_pzasy_1_70, &b_Dst_pzasy_1_70);
   fChain->SetBranchAddress("Dst_pasy_1.70", &Dst_pasy_1_70, &b_Dst_pasy_1_70);
   fChain->SetBranchAddress("Dst_ptasy_1.70", &Dst_ptasy_1_70, &b_Dst_ptasy_1_70);
   fChain->SetBranchAddress("D_MINIP", &D_MINIP, &b_D_MINIP);
   fChain->SetBranchAddress("D_MINIPCHI2", &D_MINIPCHI2, &b_D_MINIPCHI2);
   fChain->SetBranchAddress("D_MINIPNEXTBEST", &D_MINIPNEXTBEST, &b_D_MINIPNEXTBEST);
   fChain->SetBranchAddress("D_MINIPCHI2NEXTBEST", &D_MINIPCHI2NEXTBEST, &b_D_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("D_ENDVERTEX_X", &D_ENDVERTEX_X, &b_D_ENDVERTEX_X);
   fChain->SetBranchAddress("D_ENDVERTEX_Y", &D_ENDVERTEX_Y, &b_D_ENDVERTEX_Y);
   fChain->SetBranchAddress("D_ENDVERTEX_Z", &D_ENDVERTEX_Z, &b_D_ENDVERTEX_Z);
   fChain->SetBranchAddress("D_ENDVERTEX_XERR", &D_ENDVERTEX_XERR, &b_D_ENDVERTEX_XERR);
   fChain->SetBranchAddress("D_ENDVERTEX_YERR", &D_ENDVERTEX_YERR, &b_D_ENDVERTEX_YERR);
   fChain->SetBranchAddress("D_ENDVERTEX_ZERR", &D_ENDVERTEX_ZERR, &b_D_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("D_ENDVERTEX_CHI2", &D_ENDVERTEX_CHI2, &b_D_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("D_ENDVERTEX_NDOF", &D_ENDVERTEX_NDOF, &b_D_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("D_ENDVERTEX_COV_", D_ENDVERTEX_COV_, &b_D_ENDVERTEX_COV_);
   fChain->SetBranchAddress("D_OWNPV_X", &D_OWNPV_X, &b_D_OWNPV_X);
   fChain->SetBranchAddress("D_OWNPV_Y", &D_OWNPV_Y, &b_D_OWNPV_Y);
   fChain->SetBranchAddress("D_OWNPV_Z", &D_OWNPV_Z, &b_D_OWNPV_Z);
   fChain->SetBranchAddress("D_OWNPV_XERR", &D_OWNPV_XERR, &b_D_OWNPV_XERR);
   fChain->SetBranchAddress("D_OWNPV_YERR", &D_OWNPV_YERR, &b_D_OWNPV_YERR);
   fChain->SetBranchAddress("D_OWNPV_ZERR", &D_OWNPV_ZERR, &b_D_OWNPV_ZERR);
   fChain->SetBranchAddress("D_OWNPV_CHI2", &D_OWNPV_CHI2, &b_D_OWNPV_CHI2);
   fChain->SetBranchAddress("D_OWNPV_NDOF", &D_OWNPV_NDOF, &b_D_OWNPV_NDOF);
   fChain->SetBranchAddress("D_OWNPV_COV_", D_OWNPV_COV_, &b_D_OWNPV_COV_);
   fChain->SetBranchAddress("D_IP_OWNPV", &D_IP_OWNPV, &b_D_IP_OWNPV);
   fChain->SetBranchAddress("D_IPCHI2_OWNPV", &D_IPCHI2_OWNPV, &b_D_IPCHI2_OWNPV);
   fChain->SetBranchAddress("D_FD_OWNPV", &D_FD_OWNPV, &b_D_FD_OWNPV);
   fChain->SetBranchAddress("D_FDCHI2_OWNPV", &D_FDCHI2_OWNPV, &b_D_FDCHI2_OWNPV);
   fChain->SetBranchAddress("D_DIRA_OWNPV", &D_DIRA_OWNPV, &b_D_DIRA_OWNPV);
   fChain->SetBranchAddress("D_TOPPV_X", &D_TOPPV_X, &b_D_TOPPV_X);
   fChain->SetBranchAddress("D_TOPPV_Y", &D_TOPPV_Y, &b_D_TOPPV_Y);
   fChain->SetBranchAddress("D_TOPPV_Z", &D_TOPPV_Z, &b_D_TOPPV_Z);
   fChain->SetBranchAddress("D_TOPPV_XERR", &D_TOPPV_XERR, &b_D_TOPPV_XERR);
   fChain->SetBranchAddress("D_TOPPV_YERR", &D_TOPPV_YERR, &b_D_TOPPV_YERR);
   fChain->SetBranchAddress("D_TOPPV_ZERR", &D_TOPPV_ZERR, &b_D_TOPPV_ZERR);
   fChain->SetBranchAddress("D_TOPPV_CHI2", &D_TOPPV_CHI2, &b_D_TOPPV_CHI2);
   fChain->SetBranchAddress("D_TOPPV_NDOF", &D_TOPPV_NDOF, &b_D_TOPPV_NDOF);
   fChain->SetBranchAddress("D_TOPPV_COV_", D_TOPPV_COV_, &b_D_TOPPV_COV_);
   fChain->SetBranchAddress("D_IP_TOPPV", &D_IP_TOPPV, &b_D_IP_TOPPV);
   fChain->SetBranchAddress("D_IPCHI2_TOPPV", &D_IPCHI2_TOPPV, &b_D_IPCHI2_TOPPV);
   fChain->SetBranchAddress("D_FD_TOPPV", &D_FD_TOPPV, &b_D_FD_TOPPV);
   fChain->SetBranchAddress("D_FDCHI2_TOPPV", &D_FDCHI2_TOPPV, &b_D_FDCHI2_TOPPV);
   fChain->SetBranchAddress("D_DIRA_TOPPV", &D_DIRA_TOPPV, &b_D_DIRA_TOPPV);
   fChain->SetBranchAddress("D_ORIVX_X", &D_ORIVX_X, &b_D_ORIVX_X);
   fChain->SetBranchAddress("D_ORIVX_Y", &D_ORIVX_Y, &b_D_ORIVX_Y);
   fChain->SetBranchAddress("D_ORIVX_Z", &D_ORIVX_Z, &b_D_ORIVX_Z);
   fChain->SetBranchAddress("D_ORIVX_XERR", &D_ORIVX_XERR, &b_D_ORIVX_XERR);
   fChain->SetBranchAddress("D_ORIVX_YERR", &D_ORIVX_YERR, &b_D_ORIVX_YERR);
   fChain->SetBranchAddress("D_ORIVX_ZERR", &D_ORIVX_ZERR, &b_D_ORIVX_ZERR);
   fChain->SetBranchAddress("D_ORIVX_CHI2", &D_ORIVX_CHI2, &b_D_ORIVX_CHI2);
   fChain->SetBranchAddress("D_ORIVX_NDOF", &D_ORIVX_NDOF, &b_D_ORIVX_NDOF);
   fChain->SetBranchAddress("D_ORIVX_COV_", D_ORIVX_COV_, &b_D_ORIVX_COV_);
   fChain->SetBranchAddress("D_IP_ORIVX", &D_IP_ORIVX, &b_D_IP_ORIVX);
   fChain->SetBranchAddress("D_IPCHI2_ORIVX", &D_IPCHI2_ORIVX, &b_D_IPCHI2_ORIVX);
   fChain->SetBranchAddress("D_FD_ORIVX", &D_FD_ORIVX, &b_D_FD_ORIVX);
   fChain->SetBranchAddress("D_FDCHI2_ORIVX", &D_FDCHI2_ORIVX, &b_D_FDCHI2_ORIVX);
   fChain->SetBranchAddress("D_DIRA_ORIVX", &D_DIRA_ORIVX, &b_D_DIRA_ORIVX);
   fChain->SetBranchAddress("D_P", &D_P, &b_D_P);
   fChain->SetBranchAddress("D_PT", &D_PT, &b_D_PT);
   fChain->SetBranchAddress("D_PE", &D_PE, &b_D_PE);
   fChain->SetBranchAddress("D_PX", &D_PX, &b_D_PX);
   fChain->SetBranchAddress("D_PY", &D_PY, &b_D_PY);
   fChain->SetBranchAddress("D_PZ", &D_PZ, &b_D_PZ);
   fChain->SetBranchAddress("D_MM", &D_MM, &b_D_MM);
   fChain->SetBranchAddress("D_MMERR", &D_MMERR, &b_D_MMERR);
   fChain->SetBranchAddress("D_M", &D_M, &b_D_M);
   fChain->SetBranchAddress("D_ID", &D_ID, &b_D_ID);
   fChain->SetBranchAddress("D_TAU", &D_TAU, &b_D_TAU);
   fChain->SetBranchAddress("D_TAUERR", &D_TAUERR, &b_D_TAUERR);
   fChain->SetBranchAddress("D_TAUCHI2", &D_TAUCHI2, &b_D_TAUCHI2);
   fChain->SetBranchAddress("D_L0Global_Dec", &D_L0Global_Dec, &b_D_L0Global_Dec);
   fChain->SetBranchAddress("D_L0Global_TIS", &D_L0Global_TIS, &b_D_L0Global_TIS);
   fChain->SetBranchAddress("D_L0Global_TOS", &D_L0Global_TOS, &b_D_L0Global_TOS);
   fChain->SetBranchAddress("D_Hlt1Global_Dec", &D_Hlt1Global_Dec, &b_D_Hlt1Global_Dec);
   fChain->SetBranchAddress("D_Hlt1Global_TIS", &D_Hlt1Global_TIS, &b_D_Hlt1Global_TIS);
   fChain->SetBranchAddress("D_Hlt1Global_TOS", &D_Hlt1Global_TOS, &b_D_Hlt1Global_TOS);
   fChain->SetBranchAddress("D_Hlt1Phys_Dec", &D_Hlt1Phys_Dec, &b_D_Hlt1Phys_Dec);
   fChain->SetBranchAddress("D_Hlt1Phys_TIS", &D_Hlt1Phys_TIS, &b_D_Hlt1Phys_TIS);
   fChain->SetBranchAddress("D_Hlt1Phys_TOS", &D_Hlt1Phys_TOS, &b_D_Hlt1Phys_TOS);
   fChain->SetBranchAddress("D_Hlt2Global_Dec", &D_Hlt2Global_Dec, &b_D_Hlt2Global_Dec);
   fChain->SetBranchAddress("D_Hlt2Global_TIS", &D_Hlt2Global_TIS, &b_D_Hlt2Global_TIS);
   fChain->SetBranchAddress("D_Hlt2Global_TOS", &D_Hlt2Global_TOS, &b_D_Hlt2Global_TOS);
   fChain->SetBranchAddress("D_Hlt2Phys_Dec", &D_Hlt2Phys_Dec, &b_D_Hlt2Phys_Dec);
   fChain->SetBranchAddress("D_Hlt2Phys_TIS", &D_Hlt2Phys_TIS, &b_D_Hlt2Phys_TIS);
   fChain->SetBranchAddress("D_Hlt2Phys_TOS", &D_Hlt2Phys_TOS, &b_D_Hlt2Phys_TOS);
   fChain->SetBranchAddress("D_L0HadronDecision_Dec", &D_L0HadronDecision_Dec, &b_D_L0HadronDecision_Dec);
   fChain->SetBranchAddress("D_L0HadronDecision_TIS", &D_L0HadronDecision_TIS, &b_D_L0HadronDecision_TIS);
   fChain->SetBranchAddress("D_L0HadronDecision_TOS", &D_L0HadronDecision_TOS, &b_D_L0HadronDecision_TOS);
   fChain->SetBranchAddress("D_L0MuonDecision_Dec", &D_L0MuonDecision_Dec, &b_D_L0MuonDecision_Dec);
   fChain->SetBranchAddress("D_L0MuonDecision_TIS", &D_L0MuonDecision_TIS, &b_D_L0MuonDecision_TIS);
   fChain->SetBranchAddress("D_L0MuonDecision_TOS", &D_L0MuonDecision_TOS, &b_D_L0MuonDecision_TOS);
   fChain->SetBranchAddress("D_L0DiMuonDecision_Dec", &D_L0DiMuonDecision_Dec, &b_D_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("D_L0DiMuonDecision_TIS", &D_L0DiMuonDecision_TIS, &b_D_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("D_L0DiMuonDecision_TOS", &D_L0DiMuonDecision_TOS, &b_D_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("D_L0ElectronDecision_Dec", &D_L0ElectronDecision_Dec, &b_D_L0ElectronDecision_Dec);
   fChain->SetBranchAddress("D_L0ElectronDecision_TIS", &D_L0ElectronDecision_TIS, &b_D_L0ElectronDecision_TIS);
   fChain->SetBranchAddress("D_L0ElectronDecision_TOS", &D_L0ElectronDecision_TOS, &b_D_L0ElectronDecision_TOS);
   fChain->SetBranchAddress("D_L0PhotonDecision_Dec", &D_L0PhotonDecision_Dec, &b_D_L0PhotonDecision_Dec);
   fChain->SetBranchAddress("D_L0PhotonDecision_TIS", &D_L0PhotonDecision_TIS, &b_D_L0PhotonDecision_TIS);
   fChain->SetBranchAddress("D_L0PhotonDecision_TOS", &D_L0PhotonDecision_TOS, &b_D_L0PhotonDecision_TOS);
   fChain->SetBranchAddress("D_Hlt1DiMuonHighMassDecision_Dec", &D_Hlt1DiMuonHighMassDecision_Dec, &b_D_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("D_Hlt1DiMuonHighMassDecision_TIS", &D_Hlt1DiMuonHighMassDecision_TIS, &b_D_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("D_Hlt1DiMuonHighMassDecision_TOS", &D_Hlt1DiMuonHighMassDecision_TOS, &b_D_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("D_Hlt1DiMuonLowMassDecision_Dec", &D_Hlt1DiMuonLowMassDecision_Dec, &b_D_Hlt1DiMuonLowMassDecision_Dec);
   fChain->SetBranchAddress("D_Hlt1DiMuonLowMassDecision_TIS", &D_Hlt1DiMuonLowMassDecision_TIS, &b_D_Hlt1DiMuonLowMassDecision_TIS);
   fChain->SetBranchAddress("D_Hlt1DiMuonLowMassDecision_TOS", &D_Hlt1DiMuonLowMassDecision_TOS, &b_D_Hlt1DiMuonLowMassDecision_TOS);
   fChain->SetBranchAddress("D_Hlt1SingleMuonNoIPDecision_Dec", &D_Hlt1SingleMuonNoIPDecision_Dec, &b_D_Hlt1SingleMuonNoIPDecision_Dec);
   fChain->SetBranchAddress("D_Hlt1SingleMuonNoIPDecision_TIS", &D_Hlt1SingleMuonNoIPDecision_TIS, &b_D_Hlt1SingleMuonNoIPDecision_TIS);
   fChain->SetBranchAddress("D_Hlt1SingleMuonNoIPDecision_TOS", &D_Hlt1SingleMuonNoIPDecision_TOS, &b_D_Hlt1SingleMuonNoIPDecision_TOS);
   fChain->SetBranchAddress("D_Hlt1SingleMuonHighPTDecision_Dec", &D_Hlt1SingleMuonHighPTDecision_Dec, &b_D_Hlt1SingleMuonHighPTDecision_Dec);
   fChain->SetBranchAddress("D_Hlt1SingleMuonHighPTDecision_TIS", &D_Hlt1SingleMuonHighPTDecision_TIS, &b_D_Hlt1SingleMuonHighPTDecision_TIS);
   fChain->SetBranchAddress("D_Hlt1SingleMuonHighPTDecision_TOS", &D_Hlt1SingleMuonHighPTDecision_TOS, &b_D_Hlt1SingleMuonHighPTDecision_TOS);
   fChain->SetBranchAddress("D_Hlt1TrackAllL0Decision_Dec", &D_Hlt1TrackAllL0Decision_Dec, &b_D_Hlt1TrackAllL0Decision_Dec);
   fChain->SetBranchAddress("D_Hlt1TrackAllL0Decision_TIS", &D_Hlt1TrackAllL0Decision_TIS, &b_D_Hlt1TrackAllL0Decision_TIS);
   fChain->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS", &D_Hlt1TrackAllL0Decision_TOS, &b_D_Hlt1TrackAllL0Decision_TOS);
   fChain->SetBranchAddress("D_Hlt1TrackMuonDecision_Dec", &D_Hlt1TrackMuonDecision_Dec, &b_D_Hlt1TrackMuonDecision_Dec);
   fChain->SetBranchAddress("D_Hlt1TrackMuonDecision_TIS", &D_Hlt1TrackMuonDecision_TIS, &b_D_Hlt1TrackMuonDecision_TIS);
   fChain->SetBranchAddress("D_Hlt1TrackMuonDecision_TOS", &D_Hlt1TrackMuonDecision_TOS, &b_D_Hlt1TrackMuonDecision_TOS);
   fChain->SetBranchAddress("D_Hlt1TrackPhotonDecision_Dec", &D_Hlt1TrackPhotonDecision_Dec, &b_D_Hlt1TrackPhotonDecision_Dec);
   fChain->SetBranchAddress("D_Hlt1TrackPhotonDecision_TIS", &D_Hlt1TrackPhotonDecision_TIS, &b_D_Hlt1TrackPhotonDecision_TIS);
   fChain->SetBranchAddress("D_Hlt1TrackPhotonDecision_TOS", &D_Hlt1TrackPhotonDecision_TOS, &b_D_Hlt1TrackPhotonDecision_TOS);
   fChain->SetBranchAddress("D_Hlt1L0AnyDecision_Dec", &D_Hlt1L0AnyDecision_Dec, &b_D_Hlt1L0AnyDecision_Dec);
   fChain->SetBranchAddress("D_Hlt1L0AnyDecision_TIS", &D_Hlt1L0AnyDecision_TIS, &b_D_Hlt1L0AnyDecision_TIS);
   fChain->SetBranchAddress("D_Hlt1L0AnyDecision_TOS", &D_Hlt1L0AnyDecision_TOS, &b_D_Hlt1L0AnyDecision_TOS);
   fChain->SetBranchAddress("D_Hlt1GlobalDecision_Dec", &D_Hlt1GlobalDecision_Dec, &b_D_Hlt1GlobalDecision_Dec);
   fChain->SetBranchAddress("D_Hlt1GlobalDecision_TIS", &D_Hlt1GlobalDecision_TIS, &b_D_Hlt1GlobalDecision_TIS);
   fChain->SetBranchAddress("D_Hlt1GlobalDecision_TOS", &D_Hlt1GlobalDecision_TOS, &b_D_Hlt1GlobalDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2SingleMuonDecision_Dec", &D_Hlt2SingleMuonDecision_Dec, &b_D_Hlt2SingleMuonDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2SingleMuonDecision_TIS", &D_Hlt2SingleMuonDecision_TIS, &b_D_Hlt2SingleMuonDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2SingleMuonDecision_TOS", &D_Hlt2SingleMuonDecision_TOS, &b_D_Hlt2SingleMuonDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2DiMuonDetachedDecision_Dec", &D_Hlt2DiMuonDetachedDecision_Dec, &b_D_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2DiMuonDetachedDecision_TIS", &D_Hlt2DiMuonDetachedDecision_TIS, &b_D_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2DiMuonDetachedDecision_TOS", &D_Hlt2DiMuonDetachedDecision_TOS, &b_D_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD2HMuMuDecision_Dec", &D_Hlt2CharmSemilepD2HMuMuDecision_Dec, &b_D_Hlt2CharmSemilepD2HMuMuDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD2HMuMuDecision_TIS", &D_Hlt2CharmSemilepD2HMuMuDecision_TIS, &b_D_Hlt2CharmSemilepD2HMuMuDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD2HMuMuDecision_TOS", &D_Hlt2CharmSemilepD2HMuMuDecision_TOS, &b_D_Hlt2CharmSemilepD2HMuMuDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec", &D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec, &b_D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS", &D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS, &b_D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS", &D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS, &b_D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec", &D_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec, &b_D_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS", &D_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS, &b_D_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS", &D_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS, &b_D_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec", &D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec, &b_D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS", &D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS, &b_D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS", &D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS, &b_D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec", &D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec, &b_D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS", &D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS, &b_D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS", &D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS, &b_D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec", &D_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec, &b_D_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS", &D_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS, &b_D_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS", &D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS, &b_D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD02KKMuMuDecision_Dec", &D_Hlt2CharmSemilepD02KKMuMuDecision_Dec, &b_D_Hlt2CharmSemilepD02KKMuMuDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD02KKMuMuDecision_TIS", &D_Hlt2CharmSemilepD02KKMuMuDecision_TIS, &b_D_Hlt2CharmSemilepD02KKMuMuDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD02KKMuMuDecision_TOS", &D_Hlt2CharmSemilepD02KKMuMuDecision_TOS, &b_D_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_Dec", &D_Hlt2CharmSemilepD02KPiMuMuDecision_Dec, &b_D_Hlt2CharmSemilepD02KPiMuMuDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_TIS", &D_Hlt2CharmSemilepD02KPiMuMuDecision_TIS, &b_D_Hlt2CharmSemilepD02KPiMuMuDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS", &D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS, &b_D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHHDst_4piDecision_Dec", &D_Hlt2CharmHadD02HHHHDst_4piDecision_Dec, &b_D_Hlt2CharmHadD02HHHHDst_4piDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHHDst_4piDecision_TIS", &D_Hlt2CharmHadD02HHHHDst_4piDecision_TIS, &b_D_Hlt2CharmHadD02HHHHDst_4piDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHHDst_4piDecision_TOS", &D_Hlt2CharmHadD02HHHHDst_4piDecision_TOS, &b_D_Hlt2CharmHadD02HHHHDst_4piDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec", &D_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec, &b_D_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS", &D_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS, &b_D_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS", &D_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS, &b_D_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec", &D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec, &b_D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS", &D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS, &b_D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS", &D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS, &b_D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_hhXDecision_Dec", &D_Hlt2CharmHadD02HHXDst_hhXDecision_Dec, &b_D_Hlt2CharmHadD02HHXDst_hhXDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_hhXDecision_TIS", &D_Hlt2CharmHadD02HHXDst_hhXDecision_TIS, &b_D_Hlt2CharmHadD02HHXDst_hhXDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_hhXDecision_TOS", &D_Hlt2CharmHadD02HHXDst_hhXDecision_TOS, &b_D_Hlt2CharmHadD02HHXDst_hhXDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec", &D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec, &b_D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS", &D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS, &b_D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS", &D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS, &b_D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec", &D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec, &b_D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS", &D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS, &b_D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS", &D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS, &b_D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec", &D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec, &b_D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS", &D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS, &b_D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS", &D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS, &b_D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec", &D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec, &b_D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS", &D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS, &b_D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS", &D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS, &b_D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec", &D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec, &b_D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS", &D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS, &b_D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS", &D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS, &b_D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_K3piDecision_Dec", &D_Hlt2CharmHadD02HHHH_K3piDecision_Dec, &b_D_Hlt2CharmHadD02HHHH_K3piDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_K3piDecision_TIS", &D_Hlt2CharmHadD02HHHH_K3piDecision_TIS, &b_D_Hlt2CharmHadD02HHHH_K3piDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_K3piDecision_TOS", &D_Hlt2CharmHadD02HHHH_K3piDecision_TOS, &b_D_Hlt2CharmHadD02HHHH_K3piDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec", &D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec, &b_D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS", &D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS, &b_D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS", &D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS, &b_D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec", &D_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec, &b_D_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS", &D_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS, &b_D_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS", &D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS, &b_D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec", &D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec, &b_D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS", &D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS, &b_D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS", &D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS, &b_D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_4piDecision_Dec", &D_Hlt2CharmHadD02HHHH_4piDecision_Dec, &b_D_Hlt2CharmHadD02HHHH_4piDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_4piDecision_TIS", &D_Hlt2CharmHadD02HHHH_4piDecision_TIS, &b_D_Hlt2CharmHadD02HHHH_4piDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_4piDecision_TOS", &D_Hlt2CharmHadD02HHHH_4piDecision_TOS, &b_D_Hlt2CharmHadD02HHHH_4piDecision_TOS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec", &D_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec, &b_D_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS", &D_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS, &b_D_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS);
   fChain->SetBranchAddress("D_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS", &D_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS, &b_D_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS);
   fChain->SetBranchAddress("D_cpx_1.00", &D_cpx_1_00, &b_D_cpx_1_00);
   fChain->SetBranchAddress("D_cpy_1.00", &D_cpy_1_00, &b_D_cpy_1_00);
   fChain->SetBranchAddress("D_cpz_1.00", &D_cpz_1_00, &b_D_cpz_1_00);
   fChain->SetBranchAddress("D_cpt_1.00", &D_cpt_1_00, &b_D_cpt_1_00);
   fChain->SetBranchAddress("D_cp_1.00", &D_cp_1_00, &b_D_cp_1_00);
   fChain->SetBranchAddress("D_cmult_1.00", &D_cmult_1_00, &b_D_cmult_1_00);
   fChain->SetBranchAddress("D_deltaEta_1.00", &D_deltaEta_1_00, &b_D_deltaEta_1_00);
   fChain->SetBranchAddress("D_deltaPhi_1.00", &D_deltaPhi_1_00, &b_D_deltaPhi_1_00);
   fChain->SetBranchAddress("D_pxasy_1.00", &D_pxasy_1_00, &b_D_pxasy_1_00);
   fChain->SetBranchAddress("D_pyasy_1.00", &D_pyasy_1_00, &b_D_pyasy_1_00);
   fChain->SetBranchAddress("D_pzasy_1.00", &D_pzasy_1_00, &b_D_pzasy_1_00);
   fChain->SetBranchAddress("D_pasy_1.00", &D_pasy_1_00, &b_D_pasy_1_00);
   fChain->SetBranchAddress("D_ptasy_1.00", &D_ptasy_1_00, &b_D_ptasy_1_00);
   fChain->SetBranchAddress("D_cpx_1.10", &D_cpx_1_10, &b_D_cpx_1_10);
   fChain->SetBranchAddress("D_cpy_1.10", &D_cpy_1_10, &b_D_cpy_1_10);
   fChain->SetBranchAddress("D_cpz_1.10", &D_cpz_1_10, &b_D_cpz_1_10);
   fChain->SetBranchAddress("D_cpt_1.10", &D_cpt_1_10, &b_D_cpt_1_10);
   fChain->SetBranchAddress("D_cp_1.10", &D_cp_1_10, &b_D_cp_1_10);
   fChain->SetBranchAddress("D_cmult_1.10", &D_cmult_1_10, &b_D_cmult_1_10);
   fChain->SetBranchAddress("D_deltaEta_1.10", &D_deltaEta_1_10, &b_D_deltaEta_1_10);
   fChain->SetBranchAddress("D_deltaPhi_1.10", &D_deltaPhi_1_10, &b_D_deltaPhi_1_10);
   fChain->SetBranchAddress("D_pxasy_1.10", &D_pxasy_1_10, &b_D_pxasy_1_10);
   fChain->SetBranchAddress("D_pyasy_1.10", &D_pyasy_1_10, &b_D_pyasy_1_10);
   fChain->SetBranchAddress("D_pzasy_1.10", &D_pzasy_1_10, &b_D_pzasy_1_10);
   fChain->SetBranchAddress("D_pasy_1.10", &D_pasy_1_10, &b_D_pasy_1_10);
   fChain->SetBranchAddress("D_ptasy_1.10", &D_ptasy_1_10, &b_D_ptasy_1_10);
   fChain->SetBranchAddress("D_cpx_1.20", &D_cpx_1_20, &b_D_cpx_1_20);
   fChain->SetBranchAddress("D_cpy_1.20", &D_cpy_1_20, &b_D_cpy_1_20);
   fChain->SetBranchAddress("D_cpz_1.20", &D_cpz_1_20, &b_D_cpz_1_20);
   fChain->SetBranchAddress("D_cpt_1.20", &D_cpt_1_20, &b_D_cpt_1_20);
   fChain->SetBranchAddress("D_cp_1.20", &D_cp_1_20, &b_D_cp_1_20);
   fChain->SetBranchAddress("D_cmult_1.20", &D_cmult_1_20, &b_D_cmult_1_20);
   fChain->SetBranchAddress("D_deltaEta_1.20", &D_deltaEta_1_20, &b_D_deltaEta_1_20);
   fChain->SetBranchAddress("D_deltaPhi_1.20", &D_deltaPhi_1_20, &b_D_deltaPhi_1_20);
   fChain->SetBranchAddress("D_pxasy_1.20", &D_pxasy_1_20, &b_D_pxasy_1_20);
   fChain->SetBranchAddress("D_pyasy_1.20", &D_pyasy_1_20, &b_D_pyasy_1_20);
   fChain->SetBranchAddress("D_pzasy_1.20", &D_pzasy_1_20, &b_D_pzasy_1_20);
   fChain->SetBranchAddress("D_pasy_1.20", &D_pasy_1_20, &b_D_pasy_1_20);
   fChain->SetBranchAddress("D_ptasy_1.20", &D_ptasy_1_20, &b_D_ptasy_1_20);
   fChain->SetBranchAddress("D_cpx_1.30", &D_cpx_1_30, &b_D_cpx_1_30);
   fChain->SetBranchAddress("D_cpy_1.30", &D_cpy_1_30, &b_D_cpy_1_30);
   fChain->SetBranchAddress("D_cpz_1.30", &D_cpz_1_30, &b_D_cpz_1_30);
   fChain->SetBranchAddress("D_cpt_1.30", &D_cpt_1_30, &b_D_cpt_1_30);
   fChain->SetBranchAddress("D_cp_1.30", &D_cp_1_30, &b_D_cp_1_30);
   fChain->SetBranchAddress("D_cmult_1.30", &D_cmult_1_30, &b_D_cmult_1_30);
   fChain->SetBranchAddress("D_deltaEta_1.30", &D_deltaEta_1_30, &b_D_deltaEta_1_30);
   fChain->SetBranchAddress("D_deltaPhi_1.30", &D_deltaPhi_1_30, &b_D_deltaPhi_1_30);
   fChain->SetBranchAddress("D_pxasy_1.30", &D_pxasy_1_30, &b_D_pxasy_1_30);
   fChain->SetBranchAddress("D_pyasy_1.30", &D_pyasy_1_30, &b_D_pyasy_1_30);
   fChain->SetBranchAddress("D_pzasy_1.30", &D_pzasy_1_30, &b_D_pzasy_1_30);
   fChain->SetBranchAddress("D_pasy_1.30", &D_pasy_1_30, &b_D_pasy_1_30);
   fChain->SetBranchAddress("D_ptasy_1.30", &D_ptasy_1_30, &b_D_ptasy_1_30);
   fChain->SetBranchAddress("D_cpx_1.40", &D_cpx_1_40, &b_D_cpx_1_40);
   fChain->SetBranchAddress("D_cpy_1.40", &D_cpy_1_40, &b_D_cpy_1_40);
   fChain->SetBranchAddress("D_cpz_1.40", &D_cpz_1_40, &b_D_cpz_1_40);
   fChain->SetBranchAddress("D_cpt_1.40", &D_cpt_1_40, &b_D_cpt_1_40);
   fChain->SetBranchAddress("D_cp_1.40", &D_cp_1_40, &b_D_cp_1_40);
   fChain->SetBranchAddress("D_cmult_1.40", &D_cmult_1_40, &b_D_cmult_1_40);
   fChain->SetBranchAddress("D_deltaEta_1.40", &D_deltaEta_1_40, &b_D_deltaEta_1_40);
   fChain->SetBranchAddress("D_deltaPhi_1.40", &D_deltaPhi_1_40, &b_D_deltaPhi_1_40);
   fChain->SetBranchAddress("D_pxasy_1.40", &D_pxasy_1_40, &b_D_pxasy_1_40);
   fChain->SetBranchAddress("D_pyasy_1.40", &D_pyasy_1_40, &b_D_pyasy_1_40);
   fChain->SetBranchAddress("D_pzasy_1.40", &D_pzasy_1_40, &b_D_pzasy_1_40);
   fChain->SetBranchAddress("D_pasy_1.40", &D_pasy_1_40, &b_D_pasy_1_40);
   fChain->SetBranchAddress("D_ptasy_1.40", &D_ptasy_1_40, &b_D_ptasy_1_40);
   fChain->SetBranchAddress("D_cpx_1.50", &D_cpx_1_50, &b_D_cpx_1_50);
   fChain->SetBranchAddress("D_cpy_1.50", &D_cpy_1_50, &b_D_cpy_1_50);
   fChain->SetBranchAddress("D_cpz_1.50", &D_cpz_1_50, &b_D_cpz_1_50);
   fChain->SetBranchAddress("D_cpt_1.50", &D_cpt_1_50, &b_D_cpt_1_50);
   fChain->SetBranchAddress("D_cp_1.50", &D_cp_1_50, &b_D_cp_1_50);
   fChain->SetBranchAddress("D_cmult_1.50", &D_cmult_1_50, &b_D_cmult_1_50);
   fChain->SetBranchAddress("D_deltaEta_1.50", &D_deltaEta_1_50, &b_D_deltaEta_1_50);
   fChain->SetBranchAddress("D_deltaPhi_1.50", &D_deltaPhi_1_50, &b_D_deltaPhi_1_50);
   fChain->SetBranchAddress("D_pxasy_1.50", &D_pxasy_1_50, &b_D_pxasy_1_50);
   fChain->SetBranchAddress("D_pyasy_1.50", &D_pyasy_1_50, &b_D_pyasy_1_50);
   fChain->SetBranchAddress("D_pzasy_1.50", &D_pzasy_1_50, &b_D_pzasy_1_50);
   fChain->SetBranchAddress("D_pasy_1.50", &D_pasy_1_50, &b_D_pasy_1_50);
   fChain->SetBranchAddress("D_ptasy_1.50", &D_ptasy_1_50, &b_D_ptasy_1_50);
   fChain->SetBranchAddress("D_cpx_1.60", &D_cpx_1_60, &b_D_cpx_1_60);
   fChain->SetBranchAddress("D_cpy_1.60", &D_cpy_1_60, &b_D_cpy_1_60);
   fChain->SetBranchAddress("D_cpz_1.60", &D_cpz_1_60, &b_D_cpz_1_60);
   fChain->SetBranchAddress("D_cpt_1.60", &D_cpt_1_60, &b_D_cpt_1_60);
   fChain->SetBranchAddress("D_cp_1.60", &D_cp_1_60, &b_D_cp_1_60);
   fChain->SetBranchAddress("D_cmult_1.60", &D_cmult_1_60, &b_D_cmult_1_60);
   fChain->SetBranchAddress("D_deltaEta_1.60", &D_deltaEta_1_60, &b_D_deltaEta_1_60);
   fChain->SetBranchAddress("D_deltaPhi_1.60", &D_deltaPhi_1_60, &b_D_deltaPhi_1_60);
   fChain->SetBranchAddress("D_pxasy_1.60", &D_pxasy_1_60, &b_D_pxasy_1_60);
   fChain->SetBranchAddress("D_pyasy_1.60", &D_pyasy_1_60, &b_D_pyasy_1_60);
   fChain->SetBranchAddress("D_pzasy_1.60", &D_pzasy_1_60, &b_D_pzasy_1_60);
   fChain->SetBranchAddress("D_pasy_1.60", &D_pasy_1_60, &b_D_pasy_1_60);
   fChain->SetBranchAddress("D_ptasy_1.60", &D_ptasy_1_60, &b_D_ptasy_1_60);
   fChain->SetBranchAddress("D_cpx_1.70", &D_cpx_1_70, &b_D_cpx_1_70);
   fChain->SetBranchAddress("D_cpy_1.70", &D_cpy_1_70, &b_D_cpy_1_70);
   fChain->SetBranchAddress("D_cpz_1.70", &D_cpz_1_70, &b_D_cpz_1_70);
   fChain->SetBranchAddress("D_cpt_1.70", &D_cpt_1_70, &b_D_cpt_1_70);
   fChain->SetBranchAddress("D_cp_1.70", &D_cp_1_70, &b_D_cp_1_70);
   fChain->SetBranchAddress("D_cmult_1.70", &D_cmult_1_70, &b_D_cmult_1_70);
   fChain->SetBranchAddress("D_deltaEta_1.70", &D_deltaEta_1_70, &b_D_deltaEta_1_70);
   fChain->SetBranchAddress("D_deltaPhi_1.70", &D_deltaPhi_1_70, &b_D_deltaPhi_1_70);
   fChain->SetBranchAddress("D_pxasy_1.70", &D_pxasy_1_70, &b_D_pxasy_1_70);
   fChain->SetBranchAddress("D_pyasy_1.70", &D_pyasy_1_70, &b_D_pyasy_1_70);
   fChain->SetBranchAddress("D_pzasy_1.70", &D_pzasy_1_70, &b_D_pzasy_1_70);
   fChain->SetBranchAddress("D_pasy_1.70", &D_pasy_1_70, &b_D_pasy_1_70);
   fChain->SetBranchAddress("D_ptasy_1.70", &D_ptasy_1_70, &b_D_ptasy_1_70);
   fChain->SetBranchAddress("h0_MINIP", &h0_MINIP, &b_h0_MINIP);
   fChain->SetBranchAddress("h0_MINIPCHI2", &h0_MINIPCHI2, &b_h0_MINIPCHI2);
   fChain->SetBranchAddress("h0_MINIPNEXTBEST", &h0_MINIPNEXTBEST, &b_h0_MINIPNEXTBEST);
   fChain->SetBranchAddress("h0_MINIPCHI2NEXTBEST", &h0_MINIPCHI2NEXTBEST, &b_h0_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("h0_OWNPV_X", &h0_OWNPV_X, &b_h0_OWNPV_X);
   fChain->SetBranchAddress("h0_OWNPV_Y", &h0_OWNPV_Y, &b_h0_OWNPV_Y);
   fChain->SetBranchAddress("h0_OWNPV_Z", &h0_OWNPV_Z, &b_h0_OWNPV_Z);
   fChain->SetBranchAddress("h0_OWNPV_XERR", &h0_OWNPV_XERR, &b_h0_OWNPV_XERR);
   fChain->SetBranchAddress("h0_OWNPV_YERR", &h0_OWNPV_YERR, &b_h0_OWNPV_YERR);
   fChain->SetBranchAddress("h0_OWNPV_ZERR", &h0_OWNPV_ZERR, &b_h0_OWNPV_ZERR);
   fChain->SetBranchAddress("h0_OWNPV_CHI2", &h0_OWNPV_CHI2, &b_h0_OWNPV_CHI2);
   fChain->SetBranchAddress("h0_OWNPV_NDOF", &h0_OWNPV_NDOF, &b_h0_OWNPV_NDOF);
   fChain->SetBranchAddress("h0_OWNPV_COV_", h0_OWNPV_COV_, &b_h0_OWNPV_COV_);
   fChain->SetBranchAddress("h0_IP_OWNPV", &h0_IP_OWNPV, &b_h0_IP_OWNPV);
   fChain->SetBranchAddress("h0_IPCHI2_OWNPV", &h0_IPCHI2_OWNPV, &b_h0_IPCHI2_OWNPV);
   fChain->SetBranchAddress("h0_TOPPV_X", &h0_TOPPV_X, &b_h0_TOPPV_X);
   fChain->SetBranchAddress("h0_TOPPV_Y", &h0_TOPPV_Y, &b_h0_TOPPV_Y);
   fChain->SetBranchAddress("h0_TOPPV_Z", &h0_TOPPV_Z, &b_h0_TOPPV_Z);
   fChain->SetBranchAddress("h0_TOPPV_XERR", &h0_TOPPV_XERR, &b_h0_TOPPV_XERR);
   fChain->SetBranchAddress("h0_TOPPV_YERR", &h0_TOPPV_YERR, &b_h0_TOPPV_YERR);
   fChain->SetBranchAddress("h0_TOPPV_ZERR", &h0_TOPPV_ZERR, &b_h0_TOPPV_ZERR);
   fChain->SetBranchAddress("h0_TOPPV_CHI2", &h0_TOPPV_CHI2, &b_h0_TOPPV_CHI2);
   fChain->SetBranchAddress("h0_TOPPV_NDOF", &h0_TOPPV_NDOF, &b_h0_TOPPV_NDOF);
   fChain->SetBranchAddress("h0_TOPPV_COV_", h0_TOPPV_COV_, &b_h0_TOPPV_COV_);
   fChain->SetBranchAddress("h0_IP_TOPPV", &h0_IP_TOPPV, &b_h0_IP_TOPPV);
   fChain->SetBranchAddress("h0_IPCHI2_TOPPV", &h0_IPCHI2_TOPPV, &b_h0_IPCHI2_TOPPV);
   fChain->SetBranchAddress("h0_ORIVX_X", &h0_ORIVX_X, &b_h0_ORIVX_X);
   fChain->SetBranchAddress("h0_ORIVX_Y", &h0_ORIVX_Y, &b_h0_ORIVX_Y);
   fChain->SetBranchAddress("h0_ORIVX_Z", &h0_ORIVX_Z, &b_h0_ORIVX_Z);
   fChain->SetBranchAddress("h0_ORIVX_XERR", &h0_ORIVX_XERR, &b_h0_ORIVX_XERR);
   fChain->SetBranchAddress("h0_ORIVX_YERR", &h0_ORIVX_YERR, &b_h0_ORIVX_YERR);
   fChain->SetBranchAddress("h0_ORIVX_ZERR", &h0_ORIVX_ZERR, &b_h0_ORIVX_ZERR);
   fChain->SetBranchAddress("h0_ORIVX_CHI2", &h0_ORIVX_CHI2, &b_h0_ORIVX_CHI2);
   fChain->SetBranchAddress("h0_ORIVX_NDOF", &h0_ORIVX_NDOF, &b_h0_ORIVX_NDOF);
   fChain->SetBranchAddress("h0_ORIVX_COV_", h0_ORIVX_COV_, &b_h0_ORIVX_COV_);
   fChain->SetBranchAddress("h0_IP_ORIVX", &h0_IP_ORIVX, &b_h0_IP_ORIVX);
   fChain->SetBranchAddress("h0_IPCHI2_ORIVX", &h0_IPCHI2_ORIVX, &b_h0_IPCHI2_ORIVX);
   fChain->SetBranchAddress("h0_P", &h0_P, &b_h0_P);
   fChain->SetBranchAddress("h0_PT", &h0_PT, &b_h0_PT);
   fChain->SetBranchAddress("h0_PE", &h0_PE, &b_h0_PE);
   fChain->SetBranchAddress("h0_PX", &h0_PX, &b_h0_PX);
   fChain->SetBranchAddress("h0_PY", &h0_PY, &b_h0_PY);
   fChain->SetBranchAddress("h0_PZ", &h0_PZ, &b_h0_PZ);
   fChain->SetBranchAddress("h0_M", &h0_M, &b_h0_M);
   fChain->SetBranchAddress("h0_L0Calo_HCAL_realET", &h0_L0Calo_HCAL_realET, &b_h0_L0Calo_HCAL_realET);
   fChain->SetBranchAddress("h0_L0Calo_HCAL_xProjection", &h0_L0Calo_HCAL_xProjection, &b_h0_L0Calo_HCAL_xProjection);
   fChain->SetBranchAddress("h0_L0Calo_HCAL_yProjection", &h0_L0Calo_HCAL_yProjection, &b_h0_L0Calo_HCAL_yProjection);
   fChain->SetBranchAddress("h0_L0Calo_HCAL_region", &h0_L0Calo_HCAL_region, &b_h0_L0Calo_HCAL_region);
   fChain->SetBranchAddress("h0_L0Calo_HCAL_TriggerET", &h0_L0Calo_HCAL_TriggerET, &b_h0_L0Calo_HCAL_TriggerET);
   fChain->SetBranchAddress("h0_L0Calo_HCAL_TriggerHCALET", &h0_L0Calo_HCAL_TriggerHCALET, &b_h0_L0Calo_HCAL_TriggerHCALET);
   fChain->SetBranchAddress("h0_L0Calo_HCAL_xTrigger", &h0_L0Calo_HCAL_xTrigger, &b_h0_L0Calo_HCAL_xTrigger);
   fChain->SetBranchAddress("h0_L0Calo_HCAL_yTrigger", &h0_L0Calo_HCAL_yTrigger, &b_h0_L0Calo_HCAL_yTrigger);
   fChain->SetBranchAddress("h0_ID", &h0_ID, &b_h0_ID);
   fChain->SetBranchAddress("h0_CombDLLMu", &h0_CombDLLMu, &b_h0_CombDLLMu);
   fChain->SetBranchAddress("h0_ProbNNmu", &h0_ProbNNmu, &b_h0_ProbNNmu);
   fChain->SetBranchAddress("h0_ProbNNghost", &h0_ProbNNghost, &b_h0_ProbNNghost);
   fChain->SetBranchAddress("h0_InMuonAcc", &h0_InMuonAcc, &b_h0_InMuonAcc);
   fChain->SetBranchAddress("h0_MuonDist2", &h0_MuonDist2, &b_h0_MuonDist2);
   fChain->SetBranchAddress("h0_regionInM2", &h0_regionInM2, &b_h0_regionInM2);
   fChain->SetBranchAddress("h0_hasMuon", &h0_hasMuon, &b_h0_hasMuon);
   fChain->SetBranchAddress("h0_isMuon", &h0_isMuon, &b_h0_isMuon);
   fChain->SetBranchAddress("h0_isMuonLoose", &h0_isMuonLoose, &b_h0_isMuonLoose);
   fChain->SetBranchAddress("h0_NShared", &h0_NShared, &b_h0_NShared);
   fChain->SetBranchAddress("h0_MuonLLmu", &h0_MuonLLmu, &b_h0_MuonLLmu);
   fChain->SetBranchAddress("h0_MuonLLbg", &h0_MuonLLbg, &b_h0_MuonLLbg);
   fChain->SetBranchAddress("h0_isMuonFromProto", &h0_isMuonFromProto, &b_h0_isMuonFromProto);
   fChain->SetBranchAddress("h0_PIDe", &h0_PIDe, &b_h0_PIDe);
   fChain->SetBranchAddress("h0_PIDmu", &h0_PIDmu, &b_h0_PIDmu);
   fChain->SetBranchAddress("h0_PIDK", &h0_PIDK, &b_h0_PIDK);
   fChain->SetBranchAddress("h0_PIDp", &h0_PIDp, &b_h0_PIDp);
   fChain->SetBranchAddress("h0_ProbNNe", &h0_ProbNNe, &b_h0_ProbNNe);
   fChain->SetBranchAddress("h0_ProbNNk", &h0_ProbNNk, &b_h0_ProbNNk);
   fChain->SetBranchAddress("h0_ProbNNp", &h0_ProbNNp, &b_h0_ProbNNp);
   fChain->SetBranchAddress("h0_ProbNNpi", &h0_ProbNNpi, &b_h0_ProbNNpi);
   fChain->SetBranchAddress("h0_hasRich", &h0_hasRich, &b_h0_hasRich);
   fChain->SetBranchAddress("h0_hasCalo", &h0_hasCalo, &b_h0_hasCalo);
   fChain->SetBranchAddress("h0_UsedRichAerogel", &h0_UsedRichAerogel, &b_h0_UsedRichAerogel);
   fChain->SetBranchAddress("h0_UsedRich1Gas", &h0_UsedRich1Gas, &b_h0_UsedRich1Gas);
   fChain->SetBranchAddress("h0_UsedRich2Gas", &h0_UsedRich2Gas, &b_h0_UsedRich2Gas);
   fChain->SetBranchAddress("h0_RichAboveElThres", &h0_RichAboveElThres, &b_h0_RichAboveElThres);
   fChain->SetBranchAddress("h0_RichAboveMuThres", &h0_RichAboveMuThres, &b_h0_RichAboveMuThres);
   fChain->SetBranchAddress("h0_RichAbovePiThres", &h0_RichAbovePiThres, &b_h0_RichAbovePiThres);
   fChain->SetBranchAddress("h0_RichAboveKaThres", &h0_RichAboveKaThres, &b_h0_RichAboveKaThres);
   fChain->SetBranchAddress("h0_RichAbovePrThres", &h0_RichAbovePrThres, &b_h0_RichAbovePrThres);
   fChain->SetBranchAddress("h0_RichDLLe", &h0_RichDLLe, &b_h0_RichDLLe);
   fChain->SetBranchAddress("h0_RichDLLmu", &h0_RichDLLmu, &b_h0_RichDLLmu);
   fChain->SetBranchAddress("h0_RichDLLpi", &h0_RichDLLpi, &b_h0_RichDLLpi);
   fChain->SetBranchAddress("h0_RichDLLk", &h0_RichDLLk, &b_h0_RichDLLk);
   fChain->SetBranchAddress("h0_RichDLLp", &h0_RichDLLp, &b_h0_RichDLLp);
   fChain->SetBranchAddress("h0_RichDLLbt", &h0_RichDLLbt, &b_h0_RichDLLbt);
   fChain->SetBranchAddress("h0_InAccMuon", &h0_InAccMuon, &b_h0_InAccMuon);
   fChain->SetBranchAddress("h0_MuonMuLL", &h0_MuonMuLL, &b_h0_MuonMuLL);
   fChain->SetBranchAddress("h0_MuonBkgLL", &h0_MuonBkgLL, &b_h0_MuonBkgLL);
   fChain->SetBranchAddress("h0_MuonNShared", &h0_MuonNShared, &b_h0_MuonNShared);
   fChain->SetBranchAddress("h0_InAccEcal", &h0_InAccEcal, &b_h0_InAccEcal);
   fChain->SetBranchAddress("h0_CaloEcalE", &h0_CaloEcalE, &b_h0_CaloEcalE);
   fChain->SetBranchAddress("h0_EcalPIDe", &h0_EcalPIDe, &b_h0_EcalPIDe);
   fChain->SetBranchAddress("h0_EcalPIDmu", &h0_EcalPIDmu, &b_h0_EcalPIDmu);
   fChain->SetBranchAddress("h0_InAccHcal", &h0_InAccHcal, &b_h0_InAccHcal);
   fChain->SetBranchAddress("h0_CaloHcalE", &h0_CaloHcalE, &b_h0_CaloHcalE);
   fChain->SetBranchAddress("h0_HcalPIDe", &h0_HcalPIDe, &b_h0_HcalPIDe);
   fChain->SetBranchAddress("h0_HcalPIDmu", &h0_HcalPIDmu, &b_h0_HcalPIDmu);
   fChain->SetBranchAddress("h0_InAccPrs", &h0_InAccPrs, &b_h0_InAccPrs);
   fChain->SetBranchAddress("h0_PrsPIDe", &h0_PrsPIDe, &b_h0_PrsPIDe);
   fChain->SetBranchAddress("h0_CaloPrsE", &h0_CaloPrsE, &b_h0_CaloPrsE);
   fChain->SetBranchAddress("h0_InAccSpd", &h0_InAccSpd, &b_h0_InAccSpd);
   fChain->SetBranchAddress("h0_CaloSpdE", &h0_CaloSpdE, &b_h0_CaloSpdE);
   fChain->SetBranchAddress("h0_InAccBrem", &h0_InAccBrem, &b_h0_InAccBrem);
   fChain->SetBranchAddress("h0_BremPIDe", &h0_BremPIDe, &b_h0_BremPIDe);
   fChain->SetBranchAddress("h0_VeloCharge", &h0_VeloCharge, &b_h0_VeloCharge);
   fChain->SetBranchAddress("h0_RICHDLLe", &h0_RICHDLLe, &b_h0_RICHDLLe);
   fChain->SetBranchAddress("h0_RICHDLLmu", &h0_RICHDLLmu, &b_h0_RICHDLLmu);
   fChain->SetBranchAddress("h0_RICHDLLpi", &h0_RICHDLLpi, &b_h0_RICHDLLpi);
   fChain->SetBranchAddress("h0_RICHDLLK", &h0_RICHDLLK, &b_h0_RICHDLLK);
   fChain->SetBranchAddress("h0_RICHDLLp", &h0_RICHDLLp, &b_h0_RICHDLLp);
   fChain->SetBranchAddress("h0_RICHDLLbt", &h0_RICHDLLbt, &b_h0_RICHDLLbt);
   fChain->SetBranchAddress("h0_RICHBestID", &h0_RICHBestID, &b_h0_RICHBestID);
   fChain->SetBranchAddress("h0_RICHThreshold", &h0_RICHThreshold, &b_h0_RICHThreshold);
   fChain->SetBranchAddress("h0_RICHThresholdEl", &h0_RICHThresholdEl, &b_h0_RICHThresholdEl);
   fChain->SetBranchAddress("h0_RICHThresholdMu", &h0_RICHThresholdMu, &b_h0_RICHThresholdMu);
   fChain->SetBranchAddress("h0_RICHThresholdPi", &h0_RICHThresholdPi, &b_h0_RICHThresholdPi);
   fChain->SetBranchAddress("h0_RICHThresholdKa", &h0_RICHThresholdKa, &b_h0_RICHThresholdKa);
   fChain->SetBranchAddress("h0_RICHThresholdPr", &h0_RICHThresholdPr, &b_h0_RICHThresholdPr);
   fChain->SetBranchAddress("h0_RICHAerogelUsed", &h0_RICHAerogelUsed, &b_h0_RICHAerogelUsed);
   fChain->SetBranchAddress("h0_RICH1GasUsed", &h0_RICH1GasUsed, &b_h0_RICH1GasUsed);
   fChain->SetBranchAddress("h0_RICH2GasUsed", &h0_RICH2GasUsed, &b_h0_RICH2GasUsed);
   fChain->SetBranchAddress("h0_TRACK_Eta", &h0_TRACK_Eta, &b_h0_TRACK_Eta);
   fChain->SetBranchAddress("h0_TRACK_Phi", &h0_TRACK_Phi, &b_h0_TRACK_Phi);
   fChain->SetBranchAddress("h0_Aerogel_X", &h0_Aerogel_X, &b_h0_Aerogel_X);
   fChain->SetBranchAddress("h0_Aerogel_Y", &h0_Aerogel_Y, &b_h0_Aerogel_Y);
   fChain->SetBranchAddress("h0_Aerogel_Z", &h0_Aerogel_Z, &b_h0_Aerogel_Z);
   fChain->SetBranchAddress("h0_Aerogel_Rho", &h0_Aerogel_Rho, &b_h0_Aerogel_Rho);
   fChain->SetBranchAddress("h0_Aerogel_Phi", &h0_Aerogel_Phi, &b_h0_Aerogel_Phi);
   fChain->SetBranchAddress("h0_Rich1Gas_X", &h0_Rich1Gas_X, &b_h0_Rich1Gas_X);
   fChain->SetBranchAddress("h0_Rich1Gas_Y", &h0_Rich1Gas_Y, &b_h0_Rich1Gas_Y);
   fChain->SetBranchAddress("h0_Rich1Gas_Z", &h0_Rich1Gas_Z, &b_h0_Rich1Gas_Z);
   fChain->SetBranchAddress("h0_Rich1Gas_Rho", &h0_Rich1Gas_Rho, &b_h0_Rich1Gas_Rho);
   fChain->SetBranchAddress("h0_Rich1Gas_Phi", &h0_Rich1Gas_Phi, &b_h0_Rich1Gas_Phi);
   fChain->SetBranchAddress("h0_Rich2Gas_X", &h0_Rich2Gas_X, &b_h0_Rich2Gas_X);
   fChain->SetBranchAddress("h0_Rich2Gas_Y", &h0_Rich2Gas_Y, &b_h0_Rich2Gas_Y);
   fChain->SetBranchAddress("h0_Rich2Gas_Z", &h0_Rich2Gas_Z, &b_h0_Rich2Gas_Z);
   fChain->SetBranchAddress("h0_Rich2Gas_Rho", &h0_Rich2Gas_Rho, &b_h0_Rich2Gas_Rho);
   fChain->SetBranchAddress("h0_Rich2Gas_Phi", &h0_Rich2Gas_Phi, &b_h0_Rich2Gas_Phi);
   fChain->SetBranchAddress("h0_L0Global_Dec", &h0_L0Global_Dec, &b_h0_L0Global_Dec);
   fChain->SetBranchAddress("h0_L0Global_TIS", &h0_L0Global_TIS, &b_h0_L0Global_TIS);
   fChain->SetBranchAddress("h0_L0Global_TOS", &h0_L0Global_TOS, &b_h0_L0Global_TOS);
   fChain->SetBranchAddress("h0_Hlt1Global_Dec", &h0_Hlt1Global_Dec, &b_h0_Hlt1Global_Dec);
   fChain->SetBranchAddress("h0_Hlt1Global_TIS", &h0_Hlt1Global_TIS, &b_h0_Hlt1Global_TIS);
   fChain->SetBranchAddress("h0_Hlt1Global_TOS", &h0_Hlt1Global_TOS, &b_h0_Hlt1Global_TOS);
   fChain->SetBranchAddress("h0_Hlt1Phys_Dec", &h0_Hlt1Phys_Dec, &b_h0_Hlt1Phys_Dec);
   fChain->SetBranchAddress("h0_Hlt1Phys_TIS", &h0_Hlt1Phys_TIS, &b_h0_Hlt1Phys_TIS);
   fChain->SetBranchAddress("h0_Hlt1Phys_TOS", &h0_Hlt1Phys_TOS, &b_h0_Hlt1Phys_TOS);
   fChain->SetBranchAddress("h0_Hlt2Global_Dec", &h0_Hlt2Global_Dec, &b_h0_Hlt2Global_Dec);
   fChain->SetBranchAddress("h0_Hlt2Global_TIS", &h0_Hlt2Global_TIS, &b_h0_Hlt2Global_TIS);
   fChain->SetBranchAddress("h0_Hlt2Global_TOS", &h0_Hlt2Global_TOS, &b_h0_Hlt2Global_TOS);
   fChain->SetBranchAddress("h0_Hlt2Phys_Dec", &h0_Hlt2Phys_Dec, &b_h0_Hlt2Phys_Dec);
   fChain->SetBranchAddress("h0_Hlt2Phys_TIS", &h0_Hlt2Phys_TIS, &b_h0_Hlt2Phys_TIS);
   fChain->SetBranchAddress("h0_Hlt2Phys_TOS", &h0_Hlt2Phys_TOS, &b_h0_Hlt2Phys_TOS);
   fChain->SetBranchAddress("h0_TRACK_Type", &h0_TRACK_Type, &b_h0_TRACK_Type);
   fChain->SetBranchAddress("h0_TRACK_Key", &h0_TRACK_Key, &b_h0_TRACK_Key);
   fChain->SetBranchAddress("h0_TRACK_CHI2", &h0_TRACK_CHI2, &b_h0_TRACK_CHI2);
   fChain->SetBranchAddress("h0_TRACK_NDOF", &h0_TRACK_NDOF, &b_h0_TRACK_NDOF);
   fChain->SetBranchAddress("h0_TRACK_CHI2NDOF", &h0_TRACK_CHI2NDOF, &b_h0_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("h0_TRACK_PCHI2", &h0_TRACK_PCHI2, &b_h0_TRACK_PCHI2);
   fChain->SetBranchAddress("h0_TRACK_VeloCHI2NDOF", &h0_TRACK_VeloCHI2NDOF, &b_h0_TRACK_VeloCHI2NDOF);
   fChain->SetBranchAddress("h0_TRACK_TCHI2NDOF", &h0_TRACK_TCHI2NDOF, &b_h0_TRACK_TCHI2NDOF);
   fChain->SetBranchAddress("h0_VELO_UTID", &h0_VELO_UTID, &b_h0_VELO_UTID);
   fChain->SetBranchAddress("h0_TRACK_FirstMeasurementX", &h0_TRACK_FirstMeasurementX, &b_h0_TRACK_FirstMeasurementX);
   fChain->SetBranchAddress("h0_TRACK_FirstMeasurementY", &h0_TRACK_FirstMeasurementY, &b_h0_TRACK_FirstMeasurementY);
   fChain->SetBranchAddress("h0_TRACK_FirstMeasurementZ", &h0_TRACK_FirstMeasurementZ, &b_h0_TRACK_FirstMeasurementZ);
   fChain->SetBranchAddress("h0_TRACK_MatchCHI2", &h0_TRACK_MatchCHI2, &b_h0_TRACK_MatchCHI2);
   fChain->SetBranchAddress("h0_TRACK_GhostProb", &h0_TRACK_GhostProb, &b_h0_TRACK_GhostProb);
   fChain->SetBranchAddress("h0_TRACK_CloneDist", &h0_TRACK_CloneDist, &b_h0_TRACK_CloneDist);
   fChain->SetBranchAddress("h0_TRACK_Likelihood", &h0_TRACK_Likelihood, &b_h0_TRACK_Likelihood);
   fChain->SetBranchAddress("h0_L0HadronDecision_Dec", &h0_L0HadronDecision_Dec, &b_h0_L0HadronDecision_Dec);
   fChain->SetBranchAddress("h0_L0HadronDecision_TIS", &h0_L0HadronDecision_TIS, &b_h0_L0HadronDecision_TIS);
   fChain->SetBranchAddress("h0_L0HadronDecision_TOS", &h0_L0HadronDecision_TOS, &b_h0_L0HadronDecision_TOS);
   fChain->SetBranchAddress("h0_L0MuonDecision_Dec", &h0_L0MuonDecision_Dec, &b_h0_L0MuonDecision_Dec);
   fChain->SetBranchAddress("h0_L0MuonDecision_TIS", &h0_L0MuonDecision_TIS, &b_h0_L0MuonDecision_TIS);
   fChain->SetBranchAddress("h0_L0MuonDecision_TOS", &h0_L0MuonDecision_TOS, &b_h0_L0MuonDecision_TOS);
   fChain->SetBranchAddress("h0_L0DiMuonDecision_Dec", &h0_L0DiMuonDecision_Dec, &b_h0_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("h0_L0DiMuonDecision_TIS", &h0_L0DiMuonDecision_TIS, &b_h0_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("h0_L0DiMuonDecision_TOS", &h0_L0DiMuonDecision_TOS, &b_h0_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("h0_L0ElectronDecision_Dec", &h0_L0ElectronDecision_Dec, &b_h0_L0ElectronDecision_Dec);
   fChain->SetBranchAddress("h0_L0ElectronDecision_TIS", &h0_L0ElectronDecision_TIS, &b_h0_L0ElectronDecision_TIS);
   fChain->SetBranchAddress("h0_L0ElectronDecision_TOS", &h0_L0ElectronDecision_TOS, &b_h0_L0ElectronDecision_TOS);
   fChain->SetBranchAddress("h0_L0PhotonDecision_Dec", &h0_L0PhotonDecision_Dec, &b_h0_L0PhotonDecision_Dec);
   fChain->SetBranchAddress("h0_L0PhotonDecision_TIS", &h0_L0PhotonDecision_TIS, &b_h0_L0PhotonDecision_TIS);
   fChain->SetBranchAddress("h0_L0PhotonDecision_TOS", &h0_L0PhotonDecision_TOS, &b_h0_L0PhotonDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt1DiMuonHighMassDecision_Dec", &h0_Hlt1DiMuonHighMassDecision_Dec, &b_h0_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt1DiMuonHighMassDecision_TIS", &h0_Hlt1DiMuonHighMassDecision_TIS, &b_h0_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt1DiMuonHighMassDecision_TOS", &h0_Hlt1DiMuonHighMassDecision_TOS, &b_h0_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt1DiMuonLowMassDecision_Dec", &h0_Hlt1DiMuonLowMassDecision_Dec, &b_h0_Hlt1DiMuonLowMassDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt1DiMuonLowMassDecision_TIS", &h0_Hlt1DiMuonLowMassDecision_TIS, &b_h0_Hlt1DiMuonLowMassDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt1DiMuonLowMassDecision_TOS", &h0_Hlt1DiMuonLowMassDecision_TOS, &b_h0_Hlt1DiMuonLowMassDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt1SingleMuonNoIPDecision_Dec", &h0_Hlt1SingleMuonNoIPDecision_Dec, &b_h0_Hlt1SingleMuonNoIPDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt1SingleMuonNoIPDecision_TIS", &h0_Hlt1SingleMuonNoIPDecision_TIS, &b_h0_Hlt1SingleMuonNoIPDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt1SingleMuonNoIPDecision_TOS", &h0_Hlt1SingleMuonNoIPDecision_TOS, &b_h0_Hlt1SingleMuonNoIPDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt1SingleMuonHighPTDecision_Dec", &h0_Hlt1SingleMuonHighPTDecision_Dec, &b_h0_Hlt1SingleMuonHighPTDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt1SingleMuonHighPTDecision_TIS", &h0_Hlt1SingleMuonHighPTDecision_TIS, &b_h0_Hlt1SingleMuonHighPTDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt1SingleMuonHighPTDecision_TOS", &h0_Hlt1SingleMuonHighPTDecision_TOS, &b_h0_Hlt1SingleMuonHighPTDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt1TrackAllL0Decision_Dec", &h0_Hlt1TrackAllL0Decision_Dec, &b_h0_Hlt1TrackAllL0Decision_Dec);
   fChain->SetBranchAddress("h0_Hlt1TrackAllL0Decision_TIS", &h0_Hlt1TrackAllL0Decision_TIS, &b_h0_Hlt1TrackAllL0Decision_TIS);
   fChain->SetBranchAddress("h0_Hlt1TrackAllL0Decision_TOS", &h0_Hlt1TrackAllL0Decision_TOS, &b_h0_Hlt1TrackAllL0Decision_TOS);
   fChain->SetBranchAddress("h0_Hlt1TrackMuonDecision_Dec", &h0_Hlt1TrackMuonDecision_Dec, &b_h0_Hlt1TrackMuonDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt1TrackMuonDecision_TIS", &h0_Hlt1TrackMuonDecision_TIS, &b_h0_Hlt1TrackMuonDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt1TrackMuonDecision_TOS", &h0_Hlt1TrackMuonDecision_TOS, &b_h0_Hlt1TrackMuonDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt1TrackPhotonDecision_Dec", &h0_Hlt1TrackPhotonDecision_Dec, &b_h0_Hlt1TrackPhotonDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt1TrackPhotonDecision_TIS", &h0_Hlt1TrackPhotonDecision_TIS, &b_h0_Hlt1TrackPhotonDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt1TrackPhotonDecision_TOS", &h0_Hlt1TrackPhotonDecision_TOS, &b_h0_Hlt1TrackPhotonDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt1L0AnyDecision_Dec", &h0_Hlt1L0AnyDecision_Dec, &b_h0_Hlt1L0AnyDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt1L0AnyDecision_TIS", &h0_Hlt1L0AnyDecision_TIS, &b_h0_Hlt1L0AnyDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt1L0AnyDecision_TOS", &h0_Hlt1L0AnyDecision_TOS, &b_h0_Hlt1L0AnyDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt1GlobalDecision_Dec", &h0_Hlt1GlobalDecision_Dec, &b_h0_Hlt1GlobalDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt1GlobalDecision_TIS", &h0_Hlt1GlobalDecision_TIS, &b_h0_Hlt1GlobalDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt1GlobalDecision_TOS", &h0_Hlt1GlobalDecision_TOS, &b_h0_Hlt1GlobalDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2SingleMuonDecision_Dec", &h0_Hlt2SingleMuonDecision_Dec, &b_h0_Hlt2SingleMuonDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2SingleMuonDecision_TIS", &h0_Hlt2SingleMuonDecision_TIS, &b_h0_Hlt2SingleMuonDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2SingleMuonDecision_TOS", &h0_Hlt2SingleMuonDecision_TOS, &b_h0_Hlt2SingleMuonDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2DiMuonDetachedDecision_Dec", &h0_Hlt2DiMuonDetachedDecision_Dec, &b_h0_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2DiMuonDetachedDecision_TIS", &h0_Hlt2DiMuonDetachedDecision_TIS, &b_h0_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2DiMuonDetachedDecision_TOS", &h0_Hlt2DiMuonDetachedDecision_TOS, &b_h0_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD2HMuMuDecision_Dec", &h0_Hlt2CharmSemilepD2HMuMuDecision_Dec, &b_h0_Hlt2CharmSemilepD2HMuMuDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD2HMuMuDecision_TIS", &h0_Hlt2CharmSemilepD2HMuMuDecision_TIS, &b_h0_Hlt2CharmSemilepD2HMuMuDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD2HMuMuDecision_TOS", &h0_Hlt2CharmSemilepD2HMuMuDecision_TOS, &b_h0_Hlt2CharmSemilepD2HMuMuDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec", &h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec, &b_h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS", &h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS, &b_h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS", &h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS, &b_h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec", &h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec, &b_h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS", &h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS, &b_h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS", &h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS, &b_h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec", &h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec, &b_h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS", &h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS, &b_h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS", &h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS, &b_h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec", &h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec, &b_h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS", &h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS, &b_h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS", &h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS, &b_h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec", &h0_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec, &b_h0_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS", &h0_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS, &b_h0_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS", &h0_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS, &b_h0_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD02KKMuMuDecision_Dec", &h0_Hlt2CharmSemilepD02KKMuMuDecision_Dec, &b_h0_Hlt2CharmSemilepD02KKMuMuDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD02KKMuMuDecision_TIS", &h0_Hlt2CharmSemilepD02KKMuMuDecision_TIS, &b_h0_Hlt2CharmSemilepD02KKMuMuDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD02KKMuMuDecision_TOS", &h0_Hlt2CharmSemilepD02KKMuMuDecision_TOS, &b_h0_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD02KPiMuMuDecision_Dec", &h0_Hlt2CharmSemilepD02KPiMuMuDecision_Dec, &b_h0_Hlt2CharmSemilepD02KPiMuMuDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD02KPiMuMuDecision_TIS", &h0_Hlt2CharmSemilepD02KPiMuMuDecision_TIS, &b_h0_Hlt2CharmSemilepD02KPiMuMuDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmSemilepD02KPiMuMuDecision_TOS", &h0_Hlt2CharmSemilepD02KPiMuMuDecision_TOS, &b_h0_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHHDst_4piDecision_Dec", &h0_Hlt2CharmHadD02HHHHDst_4piDecision_Dec, &b_h0_Hlt2CharmHadD02HHHHDst_4piDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHHDst_4piDecision_TIS", &h0_Hlt2CharmHadD02HHHHDst_4piDecision_TIS, &b_h0_Hlt2CharmHadD02HHHHDst_4piDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHHDst_4piDecision_TOS", &h0_Hlt2CharmHadD02HHHHDst_4piDecision_TOS, &b_h0_Hlt2CharmHadD02HHHHDst_4piDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec", &h0_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec, &b_h0_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS", &h0_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS, &b_h0_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS", &h0_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS, &b_h0_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec", &h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec, &b_h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS", &h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS, &b_h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS", &h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS, &b_h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_hhXDecision_Dec", &h0_Hlt2CharmHadD02HHXDst_hhXDecision_Dec, &b_h0_Hlt2CharmHadD02HHXDst_hhXDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_hhXDecision_TIS", &h0_Hlt2CharmHadD02HHXDst_hhXDecision_TIS, &b_h0_Hlt2CharmHadD02HHXDst_hhXDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_hhXDecision_TOS", &h0_Hlt2CharmHadD02HHXDst_hhXDecision_TOS, &b_h0_Hlt2CharmHadD02HHXDst_hhXDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec", &h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec, &b_h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS", &h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS, &b_h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS", &h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS, &b_h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec", &h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec, &b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS", &h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS, &b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS", &h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS, &b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec", &h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec, &b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS", &h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS, &b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS", &h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS, &b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec", &h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec, &b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS", &h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS, &b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS", &h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS, &b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec", &h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec, &b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS", &h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS, &b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS", &h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS, &b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_K3piDecision_Dec", &h0_Hlt2CharmHadD02HHHH_K3piDecision_Dec, &b_h0_Hlt2CharmHadD02HHHH_K3piDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_K3piDecision_TIS", &h0_Hlt2CharmHadD02HHHH_K3piDecision_TIS, &b_h0_Hlt2CharmHadD02HHHH_K3piDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_K3piDecision_TOS", &h0_Hlt2CharmHadD02HHHH_K3piDecision_TOS, &b_h0_Hlt2CharmHadD02HHHH_K3piDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec", &h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec, &b_h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS", &h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS, &b_h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS", &h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS, &b_h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec", &h0_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec, &b_h0_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS", &h0_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS, &b_h0_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS", &h0_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS, &b_h0_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec", &h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec, &b_h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS", &h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS, &b_h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS", &h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS, &b_h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_4piDecision_Dec", &h0_Hlt2CharmHadD02HHHH_4piDecision_Dec, &b_h0_Hlt2CharmHadD02HHHH_4piDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_4piDecision_TIS", &h0_Hlt2CharmHadD02HHHH_4piDecision_TIS, &b_h0_Hlt2CharmHadD02HHHH_4piDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_4piDecision_TOS", &h0_Hlt2CharmHadD02HHHH_4piDecision_TOS, &b_h0_Hlt2CharmHadD02HHHH_4piDecision_TOS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec", &h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec, &b_h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS", &h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS, &b_h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS);
   fChain->SetBranchAddress("h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS", &h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS, &b_h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS);
   fChain->SetBranchAddress("h1_MINIP", &h1_MINIP, &b_h1_MINIP);
   fChain->SetBranchAddress("h1_MINIPCHI2", &h1_MINIPCHI2, &b_h1_MINIPCHI2);
   fChain->SetBranchAddress("h1_MINIPNEXTBEST", &h1_MINIPNEXTBEST, &b_h1_MINIPNEXTBEST);
   fChain->SetBranchAddress("h1_MINIPCHI2NEXTBEST", &h1_MINIPCHI2NEXTBEST, &b_h1_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("h1_OWNPV_X", &h1_OWNPV_X, &b_h1_OWNPV_X);
   fChain->SetBranchAddress("h1_OWNPV_Y", &h1_OWNPV_Y, &b_h1_OWNPV_Y);
   fChain->SetBranchAddress("h1_OWNPV_Z", &h1_OWNPV_Z, &b_h1_OWNPV_Z);
   fChain->SetBranchAddress("h1_OWNPV_XERR", &h1_OWNPV_XERR, &b_h1_OWNPV_XERR);
   fChain->SetBranchAddress("h1_OWNPV_YERR", &h1_OWNPV_YERR, &b_h1_OWNPV_YERR);
   fChain->SetBranchAddress("h1_OWNPV_ZERR", &h1_OWNPV_ZERR, &b_h1_OWNPV_ZERR);
   fChain->SetBranchAddress("h1_OWNPV_CHI2", &h1_OWNPV_CHI2, &b_h1_OWNPV_CHI2);
   fChain->SetBranchAddress("h1_OWNPV_NDOF", &h1_OWNPV_NDOF, &b_h1_OWNPV_NDOF);
   fChain->SetBranchAddress("h1_OWNPV_COV_", h1_OWNPV_COV_, &b_h1_OWNPV_COV_);
   fChain->SetBranchAddress("h1_IP_OWNPV", &h1_IP_OWNPV, &b_h1_IP_OWNPV);
   fChain->SetBranchAddress("h1_IPCHI2_OWNPV", &h1_IPCHI2_OWNPV, &b_h1_IPCHI2_OWNPV);
   fChain->SetBranchAddress("h1_TOPPV_X", &h1_TOPPV_X, &b_h1_TOPPV_X);
   fChain->SetBranchAddress("h1_TOPPV_Y", &h1_TOPPV_Y, &b_h1_TOPPV_Y);
   fChain->SetBranchAddress("h1_TOPPV_Z", &h1_TOPPV_Z, &b_h1_TOPPV_Z);
   fChain->SetBranchAddress("h1_TOPPV_XERR", &h1_TOPPV_XERR, &b_h1_TOPPV_XERR);
   fChain->SetBranchAddress("h1_TOPPV_YERR", &h1_TOPPV_YERR, &b_h1_TOPPV_YERR);
   fChain->SetBranchAddress("h1_TOPPV_ZERR", &h1_TOPPV_ZERR, &b_h1_TOPPV_ZERR);
   fChain->SetBranchAddress("h1_TOPPV_CHI2", &h1_TOPPV_CHI2, &b_h1_TOPPV_CHI2);
   fChain->SetBranchAddress("h1_TOPPV_NDOF", &h1_TOPPV_NDOF, &b_h1_TOPPV_NDOF);
   fChain->SetBranchAddress("h1_TOPPV_COV_", h1_TOPPV_COV_, &b_h1_TOPPV_COV_);
   fChain->SetBranchAddress("h1_IP_TOPPV", &h1_IP_TOPPV, &b_h1_IP_TOPPV);
   fChain->SetBranchAddress("h1_IPCHI2_TOPPV", &h1_IPCHI2_TOPPV, &b_h1_IPCHI2_TOPPV);
   fChain->SetBranchAddress("h1_ORIVX_X", &h1_ORIVX_X, &b_h1_ORIVX_X);
   fChain->SetBranchAddress("h1_ORIVX_Y", &h1_ORIVX_Y, &b_h1_ORIVX_Y);
   fChain->SetBranchAddress("h1_ORIVX_Z", &h1_ORIVX_Z, &b_h1_ORIVX_Z);
   fChain->SetBranchAddress("h1_ORIVX_XERR", &h1_ORIVX_XERR, &b_h1_ORIVX_XERR);
   fChain->SetBranchAddress("h1_ORIVX_YERR", &h1_ORIVX_YERR, &b_h1_ORIVX_YERR);
   fChain->SetBranchAddress("h1_ORIVX_ZERR", &h1_ORIVX_ZERR, &b_h1_ORIVX_ZERR);
   fChain->SetBranchAddress("h1_ORIVX_CHI2", &h1_ORIVX_CHI2, &b_h1_ORIVX_CHI2);
   fChain->SetBranchAddress("h1_ORIVX_NDOF", &h1_ORIVX_NDOF, &b_h1_ORIVX_NDOF);
   fChain->SetBranchAddress("h1_ORIVX_COV_", h1_ORIVX_COV_, &b_h1_ORIVX_COV_);
   fChain->SetBranchAddress("h1_IP_ORIVX", &h1_IP_ORIVX, &b_h1_IP_ORIVX);
   fChain->SetBranchAddress("h1_IPCHI2_ORIVX", &h1_IPCHI2_ORIVX, &b_h1_IPCHI2_ORIVX);
   fChain->SetBranchAddress("h1_P", &h1_P, &b_h1_P);
   fChain->SetBranchAddress("h1_PT", &h1_PT, &b_h1_PT);
   fChain->SetBranchAddress("h1_PE", &h1_PE, &b_h1_PE);
   fChain->SetBranchAddress("h1_PX", &h1_PX, &b_h1_PX);
   fChain->SetBranchAddress("h1_PY", &h1_PY, &b_h1_PY);
   fChain->SetBranchAddress("h1_PZ", &h1_PZ, &b_h1_PZ);
   fChain->SetBranchAddress("h1_M", &h1_M, &b_h1_M);
   fChain->SetBranchAddress("h1_L0Calo_HCAL_realET", &h1_L0Calo_HCAL_realET, &b_h1_L0Calo_HCAL_realET);
   fChain->SetBranchAddress("h1_L0Calo_HCAL_xProjection", &h1_L0Calo_HCAL_xProjection, &b_h1_L0Calo_HCAL_xProjection);
   fChain->SetBranchAddress("h1_L0Calo_HCAL_yProjection", &h1_L0Calo_HCAL_yProjection, &b_h1_L0Calo_HCAL_yProjection);
   fChain->SetBranchAddress("h1_L0Calo_HCAL_region", &h1_L0Calo_HCAL_region, &b_h1_L0Calo_HCAL_region);
   fChain->SetBranchAddress("h1_L0Calo_HCAL_TriggerET", &h1_L0Calo_HCAL_TriggerET, &b_h1_L0Calo_HCAL_TriggerET);
   fChain->SetBranchAddress("h1_L0Calo_HCAL_TriggerHCALET", &h1_L0Calo_HCAL_TriggerHCALET, &b_h1_L0Calo_HCAL_TriggerHCALET);
   fChain->SetBranchAddress("h1_L0Calo_HCAL_xTrigger", &h1_L0Calo_HCAL_xTrigger, &b_h1_L0Calo_HCAL_xTrigger);
   fChain->SetBranchAddress("h1_L0Calo_HCAL_yTrigger", &h1_L0Calo_HCAL_yTrigger, &b_h1_L0Calo_HCAL_yTrigger);
   fChain->SetBranchAddress("h1_ID", &h1_ID, &b_h1_ID);
   fChain->SetBranchAddress("h1_CombDLLMu", &h1_CombDLLMu, &b_h1_CombDLLMu);
   fChain->SetBranchAddress("h1_ProbNNmu", &h1_ProbNNmu, &b_h1_ProbNNmu);
   fChain->SetBranchAddress("h1_ProbNNghost", &h1_ProbNNghost, &b_h1_ProbNNghost);
   fChain->SetBranchAddress("h1_InMuonAcc", &h1_InMuonAcc, &b_h1_InMuonAcc);
   fChain->SetBranchAddress("h1_MuonDist2", &h1_MuonDist2, &b_h1_MuonDist2);
   fChain->SetBranchAddress("h1_regionInM2", &h1_regionInM2, &b_h1_regionInM2);
   fChain->SetBranchAddress("h1_hasMuon", &h1_hasMuon, &b_h1_hasMuon);
   fChain->SetBranchAddress("h1_isMuon", &h1_isMuon, &b_h1_isMuon);
   fChain->SetBranchAddress("h1_isMuonLoose", &h1_isMuonLoose, &b_h1_isMuonLoose);
   fChain->SetBranchAddress("h1_NShared", &h1_NShared, &b_h1_NShared);
   fChain->SetBranchAddress("h1_MuonLLmu", &h1_MuonLLmu, &b_h1_MuonLLmu);
   fChain->SetBranchAddress("h1_MuonLLbg", &h1_MuonLLbg, &b_h1_MuonLLbg);
   fChain->SetBranchAddress("h1_isMuonFromProto", &h1_isMuonFromProto, &b_h1_isMuonFromProto);
   fChain->SetBranchAddress("h1_PIDe", &h1_PIDe, &b_h1_PIDe);
   fChain->SetBranchAddress("h1_PIDmu", &h1_PIDmu, &b_h1_PIDmu);
   fChain->SetBranchAddress("h1_PIDK", &h1_PIDK, &b_h1_PIDK);
   fChain->SetBranchAddress("h1_PIDp", &h1_PIDp, &b_h1_PIDp);
   fChain->SetBranchAddress("h1_ProbNNe", &h1_ProbNNe, &b_h1_ProbNNe);
   fChain->SetBranchAddress("h1_ProbNNk", &h1_ProbNNk, &b_h1_ProbNNk);
   fChain->SetBranchAddress("h1_ProbNNp", &h1_ProbNNp, &b_h1_ProbNNp);
   fChain->SetBranchAddress("h1_ProbNNpi", &h1_ProbNNpi, &b_h1_ProbNNpi);
   fChain->SetBranchAddress("h1_hasRich", &h1_hasRich, &b_h1_hasRich);
   fChain->SetBranchAddress("h1_hasCalo", &h1_hasCalo, &b_h1_hasCalo);
   fChain->SetBranchAddress("h1_UsedRichAerogel", &h1_UsedRichAerogel, &b_h1_UsedRichAerogel);
   fChain->SetBranchAddress("h1_UsedRich1Gas", &h1_UsedRich1Gas, &b_h1_UsedRich1Gas);
   fChain->SetBranchAddress("h1_UsedRich2Gas", &h1_UsedRich2Gas, &b_h1_UsedRich2Gas);
   fChain->SetBranchAddress("h1_RichAboveElThres", &h1_RichAboveElThres, &b_h1_RichAboveElThres);
   fChain->SetBranchAddress("h1_RichAboveMuThres", &h1_RichAboveMuThres, &b_h1_RichAboveMuThres);
   fChain->SetBranchAddress("h1_RichAbovePiThres", &h1_RichAbovePiThres, &b_h1_RichAbovePiThres);
   fChain->SetBranchAddress("h1_RichAboveKaThres", &h1_RichAboveKaThres, &b_h1_RichAboveKaThres);
   fChain->SetBranchAddress("h1_RichAbovePrThres", &h1_RichAbovePrThres, &b_h1_RichAbovePrThres);
   fChain->SetBranchAddress("h1_RichDLLe", &h1_RichDLLe, &b_h1_RichDLLe);
   fChain->SetBranchAddress("h1_RichDLLmu", &h1_RichDLLmu, &b_h1_RichDLLmu);
   fChain->SetBranchAddress("h1_RichDLLpi", &h1_RichDLLpi, &b_h1_RichDLLpi);
   fChain->SetBranchAddress("h1_RichDLLk", &h1_RichDLLk, &b_h1_RichDLLk);
   fChain->SetBranchAddress("h1_RichDLLp", &h1_RichDLLp, &b_h1_RichDLLp);
   fChain->SetBranchAddress("h1_RichDLLbt", &h1_RichDLLbt, &b_h1_RichDLLbt);
   fChain->SetBranchAddress("h1_InAccMuon", &h1_InAccMuon, &b_h1_InAccMuon);
   fChain->SetBranchAddress("h1_MuonMuLL", &h1_MuonMuLL, &b_h1_MuonMuLL);
   fChain->SetBranchAddress("h1_MuonBkgLL", &h1_MuonBkgLL, &b_h1_MuonBkgLL);
   fChain->SetBranchAddress("h1_MuonNShared", &h1_MuonNShared, &b_h1_MuonNShared);
   fChain->SetBranchAddress("h1_InAccEcal", &h1_InAccEcal, &b_h1_InAccEcal);
   fChain->SetBranchAddress("h1_CaloEcalE", &h1_CaloEcalE, &b_h1_CaloEcalE);
   fChain->SetBranchAddress("h1_EcalPIDe", &h1_EcalPIDe, &b_h1_EcalPIDe);
   fChain->SetBranchAddress("h1_EcalPIDmu", &h1_EcalPIDmu, &b_h1_EcalPIDmu);
   fChain->SetBranchAddress("h1_InAccHcal", &h1_InAccHcal, &b_h1_InAccHcal);
   fChain->SetBranchAddress("h1_CaloHcalE", &h1_CaloHcalE, &b_h1_CaloHcalE);
   fChain->SetBranchAddress("h1_HcalPIDe", &h1_HcalPIDe, &b_h1_HcalPIDe);
   fChain->SetBranchAddress("h1_HcalPIDmu", &h1_HcalPIDmu, &b_h1_HcalPIDmu);
   fChain->SetBranchAddress("h1_InAccPrs", &h1_InAccPrs, &b_h1_InAccPrs);
   fChain->SetBranchAddress("h1_PrsPIDe", &h1_PrsPIDe, &b_h1_PrsPIDe);
   fChain->SetBranchAddress("h1_CaloPrsE", &h1_CaloPrsE, &b_h1_CaloPrsE);
   fChain->SetBranchAddress("h1_InAccSpd", &h1_InAccSpd, &b_h1_InAccSpd);
   fChain->SetBranchAddress("h1_CaloSpdE", &h1_CaloSpdE, &b_h1_CaloSpdE);
   fChain->SetBranchAddress("h1_InAccBrem", &h1_InAccBrem, &b_h1_InAccBrem);
   fChain->SetBranchAddress("h1_BremPIDe", &h1_BremPIDe, &b_h1_BremPIDe);
   fChain->SetBranchAddress("h1_VeloCharge", &h1_VeloCharge, &b_h1_VeloCharge);
   fChain->SetBranchAddress("h1_RICHDLLe", &h1_RICHDLLe, &b_h1_RICHDLLe);
   fChain->SetBranchAddress("h1_RICHDLLmu", &h1_RICHDLLmu, &b_h1_RICHDLLmu);
   fChain->SetBranchAddress("h1_RICHDLLpi", &h1_RICHDLLpi, &b_h1_RICHDLLpi);
   fChain->SetBranchAddress("h1_RICHDLLK", &h1_RICHDLLK, &b_h1_RICHDLLK);
   fChain->SetBranchAddress("h1_RICHDLLp", &h1_RICHDLLp, &b_h1_RICHDLLp);
   fChain->SetBranchAddress("h1_RICHDLLbt", &h1_RICHDLLbt, &b_h1_RICHDLLbt);
   fChain->SetBranchAddress("h1_RICHBestID", &h1_RICHBestID, &b_h1_RICHBestID);
   fChain->SetBranchAddress("h1_RICHThreshold", &h1_RICHThreshold, &b_h1_RICHThreshold);
   fChain->SetBranchAddress("h1_RICHThresholdEl", &h1_RICHThresholdEl, &b_h1_RICHThresholdEl);
   fChain->SetBranchAddress("h1_RICHThresholdMu", &h1_RICHThresholdMu, &b_h1_RICHThresholdMu);
   fChain->SetBranchAddress("h1_RICHThresholdPi", &h1_RICHThresholdPi, &b_h1_RICHThresholdPi);
   fChain->SetBranchAddress("h1_RICHThresholdKa", &h1_RICHThresholdKa, &b_h1_RICHThresholdKa);
   fChain->SetBranchAddress("h1_RICHThresholdPr", &h1_RICHThresholdPr, &b_h1_RICHThresholdPr);
   fChain->SetBranchAddress("h1_RICHAerogelUsed", &h1_RICHAerogelUsed, &b_h1_RICHAerogelUsed);
   fChain->SetBranchAddress("h1_RICH1GasUsed", &h1_RICH1GasUsed, &b_h1_RICH1GasUsed);
   fChain->SetBranchAddress("h1_RICH2GasUsed", &h1_RICH2GasUsed, &b_h1_RICH2GasUsed);
   fChain->SetBranchAddress("h1_TRACK_Eta", &h1_TRACK_Eta, &b_h1_TRACK_Eta);
   fChain->SetBranchAddress("h1_TRACK_Phi", &h1_TRACK_Phi, &b_h1_TRACK_Phi);
   fChain->SetBranchAddress("h1_Aerogel_X", &h1_Aerogel_X, &b_h1_Aerogel_X);
   fChain->SetBranchAddress("h1_Aerogel_Y", &h1_Aerogel_Y, &b_h1_Aerogel_Y);
   fChain->SetBranchAddress("h1_Aerogel_Z", &h1_Aerogel_Z, &b_h1_Aerogel_Z);
   fChain->SetBranchAddress("h1_Aerogel_Rho", &h1_Aerogel_Rho, &b_h1_Aerogel_Rho);
   fChain->SetBranchAddress("h1_Aerogel_Phi", &h1_Aerogel_Phi, &b_h1_Aerogel_Phi);
   fChain->SetBranchAddress("h1_Rich1Gas_X", &h1_Rich1Gas_X, &b_h1_Rich1Gas_X);
   fChain->SetBranchAddress("h1_Rich1Gas_Y", &h1_Rich1Gas_Y, &b_h1_Rich1Gas_Y);
   fChain->SetBranchAddress("h1_Rich1Gas_Z", &h1_Rich1Gas_Z, &b_h1_Rich1Gas_Z);
   fChain->SetBranchAddress("h1_Rich1Gas_Rho", &h1_Rich1Gas_Rho, &b_h1_Rich1Gas_Rho);
   fChain->SetBranchAddress("h1_Rich1Gas_Phi", &h1_Rich1Gas_Phi, &b_h1_Rich1Gas_Phi);
   fChain->SetBranchAddress("h1_Rich2Gas_X", &h1_Rich2Gas_X, &b_h1_Rich2Gas_X);
   fChain->SetBranchAddress("h1_Rich2Gas_Y", &h1_Rich2Gas_Y, &b_h1_Rich2Gas_Y);
   fChain->SetBranchAddress("h1_Rich2Gas_Z", &h1_Rich2Gas_Z, &b_h1_Rich2Gas_Z);
   fChain->SetBranchAddress("h1_Rich2Gas_Rho", &h1_Rich2Gas_Rho, &b_h1_Rich2Gas_Rho);
   fChain->SetBranchAddress("h1_Rich2Gas_Phi", &h1_Rich2Gas_Phi, &b_h1_Rich2Gas_Phi);
   fChain->SetBranchAddress("h1_L0Global_Dec", &h1_L0Global_Dec, &b_h1_L0Global_Dec);
   fChain->SetBranchAddress("h1_L0Global_TIS", &h1_L0Global_TIS, &b_h1_L0Global_TIS);
   fChain->SetBranchAddress("h1_L0Global_TOS", &h1_L0Global_TOS, &b_h1_L0Global_TOS);
   fChain->SetBranchAddress("h1_Hlt1Global_Dec", &h1_Hlt1Global_Dec, &b_h1_Hlt1Global_Dec);
   fChain->SetBranchAddress("h1_Hlt1Global_TIS", &h1_Hlt1Global_TIS, &b_h1_Hlt1Global_TIS);
   fChain->SetBranchAddress("h1_Hlt1Global_TOS", &h1_Hlt1Global_TOS, &b_h1_Hlt1Global_TOS);
   fChain->SetBranchAddress("h1_Hlt1Phys_Dec", &h1_Hlt1Phys_Dec, &b_h1_Hlt1Phys_Dec);
   fChain->SetBranchAddress("h1_Hlt1Phys_TIS", &h1_Hlt1Phys_TIS, &b_h1_Hlt1Phys_TIS);
   fChain->SetBranchAddress("h1_Hlt1Phys_TOS", &h1_Hlt1Phys_TOS, &b_h1_Hlt1Phys_TOS);
   fChain->SetBranchAddress("h1_Hlt2Global_Dec", &h1_Hlt2Global_Dec, &b_h1_Hlt2Global_Dec);
   fChain->SetBranchAddress("h1_Hlt2Global_TIS", &h1_Hlt2Global_TIS, &b_h1_Hlt2Global_TIS);
   fChain->SetBranchAddress("h1_Hlt2Global_TOS", &h1_Hlt2Global_TOS, &b_h1_Hlt2Global_TOS);
   fChain->SetBranchAddress("h1_Hlt2Phys_Dec", &h1_Hlt2Phys_Dec, &b_h1_Hlt2Phys_Dec);
   fChain->SetBranchAddress("h1_Hlt2Phys_TIS", &h1_Hlt2Phys_TIS, &b_h1_Hlt2Phys_TIS);
   fChain->SetBranchAddress("h1_Hlt2Phys_TOS", &h1_Hlt2Phys_TOS, &b_h1_Hlt2Phys_TOS);
   fChain->SetBranchAddress("h1_TRACK_Type", &h1_TRACK_Type, &b_h1_TRACK_Type);
   fChain->SetBranchAddress("h1_TRACK_Key", &h1_TRACK_Key, &b_h1_TRACK_Key);
   fChain->SetBranchAddress("h1_TRACK_CHI2", &h1_TRACK_CHI2, &b_h1_TRACK_CHI2);
   fChain->SetBranchAddress("h1_TRACK_NDOF", &h1_TRACK_NDOF, &b_h1_TRACK_NDOF);
   fChain->SetBranchAddress("h1_TRACK_CHI2NDOF", &h1_TRACK_CHI2NDOF, &b_h1_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("h1_TRACK_PCHI2", &h1_TRACK_PCHI2, &b_h1_TRACK_PCHI2);
   fChain->SetBranchAddress("h1_TRACK_VeloCHI2NDOF", &h1_TRACK_VeloCHI2NDOF, &b_h1_TRACK_VeloCHI2NDOF);
   fChain->SetBranchAddress("h1_TRACK_TCHI2NDOF", &h1_TRACK_TCHI2NDOF, &b_h1_TRACK_TCHI2NDOF);
   fChain->SetBranchAddress("h1_VELO_UTID", &h1_VELO_UTID, &b_h1_VELO_UTID);
   fChain->SetBranchAddress("h1_TRACK_FirstMeasurementX", &h1_TRACK_FirstMeasurementX, &b_h1_TRACK_FirstMeasurementX);
   fChain->SetBranchAddress("h1_TRACK_FirstMeasurementY", &h1_TRACK_FirstMeasurementY, &b_h1_TRACK_FirstMeasurementY);
   fChain->SetBranchAddress("h1_TRACK_FirstMeasurementZ", &h1_TRACK_FirstMeasurementZ, &b_h1_TRACK_FirstMeasurementZ);
   fChain->SetBranchAddress("h1_TRACK_MatchCHI2", &h1_TRACK_MatchCHI2, &b_h1_TRACK_MatchCHI2);
   fChain->SetBranchAddress("h1_TRACK_GhostProb", &h1_TRACK_GhostProb, &b_h1_TRACK_GhostProb);
   fChain->SetBranchAddress("h1_TRACK_CloneDist", &h1_TRACK_CloneDist, &b_h1_TRACK_CloneDist);
   fChain->SetBranchAddress("h1_TRACK_Likelihood", &h1_TRACK_Likelihood, &b_h1_TRACK_Likelihood);
   fChain->SetBranchAddress("h1_L0HadronDecision_Dec", &h1_L0HadronDecision_Dec, &b_h1_L0HadronDecision_Dec);
   fChain->SetBranchAddress("h1_L0HadronDecision_TIS", &h1_L0HadronDecision_TIS, &b_h1_L0HadronDecision_TIS);
   fChain->SetBranchAddress("h1_L0HadronDecision_TOS", &h1_L0HadronDecision_TOS, &b_h1_L0HadronDecision_TOS);
   fChain->SetBranchAddress("h1_L0MuonDecision_Dec", &h1_L0MuonDecision_Dec, &b_h1_L0MuonDecision_Dec);
   fChain->SetBranchAddress("h1_L0MuonDecision_TIS", &h1_L0MuonDecision_TIS, &b_h1_L0MuonDecision_TIS);
   fChain->SetBranchAddress("h1_L0MuonDecision_TOS", &h1_L0MuonDecision_TOS, &b_h1_L0MuonDecision_TOS);
   fChain->SetBranchAddress("h1_L0DiMuonDecision_Dec", &h1_L0DiMuonDecision_Dec, &b_h1_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("h1_L0DiMuonDecision_TIS", &h1_L0DiMuonDecision_TIS, &b_h1_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("h1_L0DiMuonDecision_TOS", &h1_L0DiMuonDecision_TOS, &b_h1_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("h1_L0ElectronDecision_Dec", &h1_L0ElectronDecision_Dec, &b_h1_L0ElectronDecision_Dec);
   fChain->SetBranchAddress("h1_L0ElectronDecision_TIS", &h1_L0ElectronDecision_TIS, &b_h1_L0ElectronDecision_TIS);
   fChain->SetBranchAddress("h1_L0ElectronDecision_TOS", &h1_L0ElectronDecision_TOS, &b_h1_L0ElectronDecision_TOS);
   fChain->SetBranchAddress("h1_L0PhotonDecision_Dec", &h1_L0PhotonDecision_Dec, &b_h1_L0PhotonDecision_Dec);
   fChain->SetBranchAddress("h1_L0PhotonDecision_TIS", &h1_L0PhotonDecision_TIS, &b_h1_L0PhotonDecision_TIS);
   fChain->SetBranchAddress("h1_L0PhotonDecision_TOS", &h1_L0PhotonDecision_TOS, &b_h1_L0PhotonDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt1DiMuonHighMassDecision_Dec", &h1_Hlt1DiMuonHighMassDecision_Dec, &b_h1_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt1DiMuonHighMassDecision_TIS", &h1_Hlt1DiMuonHighMassDecision_TIS, &b_h1_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt1DiMuonHighMassDecision_TOS", &h1_Hlt1DiMuonHighMassDecision_TOS, &b_h1_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt1DiMuonLowMassDecision_Dec", &h1_Hlt1DiMuonLowMassDecision_Dec, &b_h1_Hlt1DiMuonLowMassDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt1DiMuonLowMassDecision_TIS", &h1_Hlt1DiMuonLowMassDecision_TIS, &b_h1_Hlt1DiMuonLowMassDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt1DiMuonLowMassDecision_TOS", &h1_Hlt1DiMuonLowMassDecision_TOS, &b_h1_Hlt1DiMuonLowMassDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt1SingleMuonNoIPDecision_Dec", &h1_Hlt1SingleMuonNoIPDecision_Dec, &b_h1_Hlt1SingleMuonNoIPDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt1SingleMuonNoIPDecision_TIS", &h1_Hlt1SingleMuonNoIPDecision_TIS, &b_h1_Hlt1SingleMuonNoIPDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt1SingleMuonNoIPDecision_TOS", &h1_Hlt1SingleMuonNoIPDecision_TOS, &b_h1_Hlt1SingleMuonNoIPDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt1SingleMuonHighPTDecision_Dec", &h1_Hlt1SingleMuonHighPTDecision_Dec, &b_h1_Hlt1SingleMuonHighPTDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt1SingleMuonHighPTDecision_TIS", &h1_Hlt1SingleMuonHighPTDecision_TIS, &b_h1_Hlt1SingleMuonHighPTDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt1SingleMuonHighPTDecision_TOS", &h1_Hlt1SingleMuonHighPTDecision_TOS, &b_h1_Hlt1SingleMuonHighPTDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt1TrackAllL0Decision_Dec", &h1_Hlt1TrackAllL0Decision_Dec, &b_h1_Hlt1TrackAllL0Decision_Dec);
   fChain->SetBranchAddress("h1_Hlt1TrackAllL0Decision_TIS", &h1_Hlt1TrackAllL0Decision_TIS, &b_h1_Hlt1TrackAllL0Decision_TIS);
   fChain->SetBranchAddress("h1_Hlt1TrackAllL0Decision_TOS", &h1_Hlt1TrackAllL0Decision_TOS, &b_h1_Hlt1TrackAllL0Decision_TOS);
   fChain->SetBranchAddress("h1_Hlt1TrackMuonDecision_Dec", &h1_Hlt1TrackMuonDecision_Dec, &b_h1_Hlt1TrackMuonDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt1TrackMuonDecision_TIS", &h1_Hlt1TrackMuonDecision_TIS, &b_h1_Hlt1TrackMuonDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt1TrackMuonDecision_TOS", &h1_Hlt1TrackMuonDecision_TOS, &b_h1_Hlt1TrackMuonDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt1TrackPhotonDecision_Dec", &h1_Hlt1TrackPhotonDecision_Dec, &b_h1_Hlt1TrackPhotonDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt1TrackPhotonDecision_TIS", &h1_Hlt1TrackPhotonDecision_TIS, &b_h1_Hlt1TrackPhotonDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt1TrackPhotonDecision_TOS", &h1_Hlt1TrackPhotonDecision_TOS, &b_h1_Hlt1TrackPhotonDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt1L0AnyDecision_Dec", &h1_Hlt1L0AnyDecision_Dec, &b_h1_Hlt1L0AnyDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt1L0AnyDecision_TIS", &h1_Hlt1L0AnyDecision_TIS, &b_h1_Hlt1L0AnyDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt1L0AnyDecision_TOS", &h1_Hlt1L0AnyDecision_TOS, &b_h1_Hlt1L0AnyDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt1GlobalDecision_Dec", &h1_Hlt1GlobalDecision_Dec, &b_h1_Hlt1GlobalDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt1GlobalDecision_TIS", &h1_Hlt1GlobalDecision_TIS, &b_h1_Hlt1GlobalDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt1GlobalDecision_TOS", &h1_Hlt1GlobalDecision_TOS, &b_h1_Hlt1GlobalDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2SingleMuonDecision_Dec", &h1_Hlt2SingleMuonDecision_Dec, &b_h1_Hlt2SingleMuonDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2SingleMuonDecision_TIS", &h1_Hlt2SingleMuonDecision_TIS, &b_h1_Hlt2SingleMuonDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2SingleMuonDecision_TOS", &h1_Hlt2SingleMuonDecision_TOS, &b_h1_Hlt2SingleMuonDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2DiMuonDetachedDecision_Dec", &h1_Hlt2DiMuonDetachedDecision_Dec, &b_h1_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2DiMuonDetachedDecision_TIS", &h1_Hlt2DiMuonDetachedDecision_TIS, &b_h1_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2DiMuonDetachedDecision_TOS", &h1_Hlt2DiMuonDetachedDecision_TOS, &b_h1_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD2HMuMuDecision_Dec", &h1_Hlt2CharmSemilepD2HMuMuDecision_Dec, &b_h1_Hlt2CharmSemilepD2HMuMuDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD2HMuMuDecision_TIS", &h1_Hlt2CharmSemilepD2HMuMuDecision_TIS, &b_h1_Hlt2CharmSemilepD2HMuMuDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD2HMuMuDecision_TOS", &h1_Hlt2CharmSemilepD2HMuMuDecision_TOS, &b_h1_Hlt2CharmSemilepD2HMuMuDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec", &h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec, &b_h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS", &h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS, &b_h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS", &h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS, &b_h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec", &h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec, &b_h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS", &h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS, &b_h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS", &h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS, &b_h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec", &h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec, &b_h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS", &h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS, &b_h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS", &h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS, &b_h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec", &h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec, &b_h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS", &h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS, &b_h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS", &h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS, &b_h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec", &h1_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec, &b_h1_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS", &h1_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS, &b_h1_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS", &h1_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS, &b_h1_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD02KKMuMuDecision_Dec", &h1_Hlt2CharmSemilepD02KKMuMuDecision_Dec, &b_h1_Hlt2CharmSemilepD02KKMuMuDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD02KKMuMuDecision_TIS", &h1_Hlt2CharmSemilepD02KKMuMuDecision_TIS, &b_h1_Hlt2CharmSemilepD02KKMuMuDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD02KKMuMuDecision_TOS", &h1_Hlt2CharmSemilepD02KKMuMuDecision_TOS, &b_h1_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD02KPiMuMuDecision_Dec", &h1_Hlt2CharmSemilepD02KPiMuMuDecision_Dec, &b_h1_Hlt2CharmSemilepD02KPiMuMuDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD02KPiMuMuDecision_TIS", &h1_Hlt2CharmSemilepD02KPiMuMuDecision_TIS, &b_h1_Hlt2CharmSemilepD02KPiMuMuDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmSemilepD02KPiMuMuDecision_TOS", &h1_Hlt2CharmSemilepD02KPiMuMuDecision_TOS, &b_h1_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHHDst_4piDecision_Dec", &h1_Hlt2CharmHadD02HHHHDst_4piDecision_Dec, &b_h1_Hlt2CharmHadD02HHHHDst_4piDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHHDst_4piDecision_TIS", &h1_Hlt2CharmHadD02HHHHDst_4piDecision_TIS, &b_h1_Hlt2CharmHadD02HHHHDst_4piDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHHDst_4piDecision_TOS", &h1_Hlt2CharmHadD02HHHHDst_4piDecision_TOS, &b_h1_Hlt2CharmHadD02HHHHDst_4piDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec", &h1_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec, &b_h1_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS", &h1_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS, &b_h1_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS", &h1_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS, &b_h1_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec", &h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec, &b_h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS", &h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS, &b_h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS", &h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS, &b_h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_hhXDecision_Dec", &h1_Hlt2CharmHadD02HHXDst_hhXDecision_Dec, &b_h1_Hlt2CharmHadD02HHXDst_hhXDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_hhXDecision_TIS", &h1_Hlt2CharmHadD02HHXDst_hhXDecision_TIS, &b_h1_Hlt2CharmHadD02HHXDst_hhXDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_hhXDecision_TOS", &h1_Hlt2CharmHadD02HHXDst_hhXDecision_TOS, &b_h1_Hlt2CharmHadD02HHXDst_hhXDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec", &h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec, &b_h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS", &h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS, &b_h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS", &h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS, &b_h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec", &h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec, &b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS", &h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS, &b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS", &h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS, &b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec", &h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec, &b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS", &h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS, &b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS", &h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS, &b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec", &h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec, &b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS", &h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS, &b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS", &h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS, &b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec", &h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec, &b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS", &h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS, &b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS", &h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS, &b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_K3piDecision_Dec", &h1_Hlt2CharmHadD02HHHH_K3piDecision_Dec, &b_h1_Hlt2CharmHadD02HHHH_K3piDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_K3piDecision_TIS", &h1_Hlt2CharmHadD02HHHH_K3piDecision_TIS, &b_h1_Hlt2CharmHadD02HHHH_K3piDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_K3piDecision_TOS", &h1_Hlt2CharmHadD02HHHH_K3piDecision_TOS, &b_h1_Hlt2CharmHadD02HHHH_K3piDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec", &h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec, &b_h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS", &h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS, &b_h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS", &h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS, &b_h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec", &h1_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec, &b_h1_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS", &h1_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS, &b_h1_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS", &h1_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS, &b_h1_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec", &h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec, &b_h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS", &h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS, &b_h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS", &h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS, &b_h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_4piDecision_Dec", &h1_Hlt2CharmHadD02HHHH_4piDecision_Dec, &b_h1_Hlt2CharmHadD02HHHH_4piDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_4piDecision_TIS", &h1_Hlt2CharmHadD02HHHH_4piDecision_TIS, &b_h1_Hlt2CharmHadD02HHHH_4piDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_4piDecision_TOS", &h1_Hlt2CharmHadD02HHHH_4piDecision_TOS, &b_h1_Hlt2CharmHadD02HHHH_4piDecision_TOS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec", &h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec, &b_h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS", &h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS, &b_h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS);
   fChain->SetBranchAddress("h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS", &h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS, &b_h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS);
   fChain->SetBranchAddress("mu0_MINIP", &mu0_MINIP, &b_mu0_MINIP);
   fChain->SetBranchAddress("mu0_MINIPCHI2", &mu0_MINIPCHI2, &b_mu0_MINIPCHI2);
   fChain->SetBranchAddress("mu0_MINIPNEXTBEST", &mu0_MINIPNEXTBEST, &b_mu0_MINIPNEXTBEST);
   fChain->SetBranchAddress("mu0_MINIPCHI2NEXTBEST", &mu0_MINIPCHI2NEXTBEST, &b_mu0_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("mu0_OWNPV_X", &mu0_OWNPV_X, &b_mu0_OWNPV_X);
   fChain->SetBranchAddress("mu0_OWNPV_Y", &mu0_OWNPV_Y, &b_mu0_OWNPV_Y);
   fChain->SetBranchAddress("mu0_OWNPV_Z", &mu0_OWNPV_Z, &b_mu0_OWNPV_Z);
   fChain->SetBranchAddress("mu0_OWNPV_XERR", &mu0_OWNPV_XERR, &b_mu0_OWNPV_XERR);
   fChain->SetBranchAddress("mu0_OWNPV_YERR", &mu0_OWNPV_YERR, &b_mu0_OWNPV_YERR);
   fChain->SetBranchAddress("mu0_OWNPV_ZERR", &mu0_OWNPV_ZERR, &b_mu0_OWNPV_ZERR);
   fChain->SetBranchAddress("mu0_OWNPV_CHI2", &mu0_OWNPV_CHI2, &b_mu0_OWNPV_CHI2);
   fChain->SetBranchAddress("mu0_OWNPV_NDOF", &mu0_OWNPV_NDOF, &b_mu0_OWNPV_NDOF);
   fChain->SetBranchAddress("mu0_OWNPV_COV_", mu0_OWNPV_COV_, &b_mu0_OWNPV_COV_);
   fChain->SetBranchAddress("mu0_IP_OWNPV", &mu0_IP_OWNPV, &b_mu0_IP_OWNPV);
   fChain->SetBranchAddress("mu0_IPCHI2_OWNPV", &mu0_IPCHI2_OWNPV, &b_mu0_IPCHI2_OWNPV);
   fChain->SetBranchAddress("mu0_TOPPV_X", &mu0_TOPPV_X, &b_mu0_TOPPV_X);
   fChain->SetBranchAddress("mu0_TOPPV_Y", &mu0_TOPPV_Y, &b_mu0_TOPPV_Y);
   fChain->SetBranchAddress("mu0_TOPPV_Z", &mu0_TOPPV_Z, &b_mu0_TOPPV_Z);
   fChain->SetBranchAddress("mu0_TOPPV_XERR", &mu0_TOPPV_XERR, &b_mu0_TOPPV_XERR);
   fChain->SetBranchAddress("mu0_TOPPV_YERR", &mu0_TOPPV_YERR, &b_mu0_TOPPV_YERR);
   fChain->SetBranchAddress("mu0_TOPPV_ZERR", &mu0_TOPPV_ZERR, &b_mu0_TOPPV_ZERR);
   fChain->SetBranchAddress("mu0_TOPPV_CHI2", &mu0_TOPPV_CHI2, &b_mu0_TOPPV_CHI2);
   fChain->SetBranchAddress("mu0_TOPPV_NDOF", &mu0_TOPPV_NDOF, &b_mu0_TOPPV_NDOF);
   fChain->SetBranchAddress("mu0_TOPPV_COV_", mu0_TOPPV_COV_, &b_mu0_TOPPV_COV_);
   fChain->SetBranchAddress("mu0_IP_TOPPV", &mu0_IP_TOPPV, &b_mu0_IP_TOPPV);
   fChain->SetBranchAddress("mu0_IPCHI2_TOPPV", &mu0_IPCHI2_TOPPV, &b_mu0_IPCHI2_TOPPV);
   fChain->SetBranchAddress("mu0_ORIVX_X", &mu0_ORIVX_X, &b_mu0_ORIVX_X);
   fChain->SetBranchAddress("mu0_ORIVX_Y", &mu0_ORIVX_Y, &b_mu0_ORIVX_Y);
   fChain->SetBranchAddress("mu0_ORIVX_Z", &mu0_ORIVX_Z, &b_mu0_ORIVX_Z);
   fChain->SetBranchAddress("mu0_ORIVX_XERR", &mu0_ORIVX_XERR, &b_mu0_ORIVX_XERR);
   fChain->SetBranchAddress("mu0_ORIVX_YERR", &mu0_ORIVX_YERR, &b_mu0_ORIVX_YERR);
   fChain->SetBranchAddress("mu0_ORIVX_ZERR", &mu0_ORIVX_ZERR, &b_mu0_ORIVX_ZERR);
   fChain->SetBranchAddress("mu0_ORIVX_CHI2", &mu0_ORIVX_CHI2, &b_mu0_ORIVX_CHI2);
   fChain->SetBranchAddress("mu0_ORIVX_NDOF", &mu0_ORIVX_NDOF, &b_mu0_ORIVX_NDOF);
   fChain->SetBranchAddress("mu0_ORIVX_COV_", mu0_ORIVX_COV_, &b_mu0_ORIVX_COV_);
   fChain->SetBranchAddress("mu0_IP_ORIVX", &mu0_IP_ORIVX, &b_mu0_IP_ORIVX);
   fChain->SetBranchAddress("mu0_IPCHI2_ORIVX", &mu0_IPCHI2_ORIVX, &b_mu0_IPCHI2_ORIVX);
   fChain->SetBranchAddress("mu0_P", &mu0_P, &b_mu0_P);
   fChain->SetBranchAddress("mu0_PT", &mu0_PT, &b_mu0_PT);
   fChain->SetBranchAddress("mu0_PE", &mu0_PE, &b_mu0_PE);
   fChain->SetBranchAddress("mu0_PX", &mu0_PX, &b_mu0_PX);
   fChain->SetBranchAddress("mu0_PY", &mu0_PY, &b_mu0_PY);
   fChain->SetBranchAddress("mu0_PZ", &mu0_PZ, &b_mu0_PZ);
   fChain->SetBranchAddress("mu0_M", &mu0_M, &b_mu0_M);
   fChain->SetBranchAddress("mu0_L0Calo_HCAL_realET", &mu0_L0Calo_HCAL_realET, &b_mu0_L0Calo_HCAL_realET);
   fChain->SetBranchAddress("mu0_L0Calo_HCAL_xProjection", &mu0_L0Calo_HCAL_xProjection, &b_mu0_L0Calo_HCAL_xProjection);
   fChain->SetBranchAddress("mu0_L0Calo_HCAL_yProjection", &mu0_L0Calo_HCAL_yProjection, &b_mu0_L0Calo_HCAL_yProjection);
   fChain->SetBranchAddress("mu0_L0Calo_HCAL_region", &mu0_L0Calo_HCAL_region, &b_mu0_L0Calo_HCAL_region);
   fChain->SetBranchAddress("mu0_L0Calo_HCAL_TriggerET", &mu0_L0Calo_HCAL_TriggerET, &b_mu0_L0Calo_HCAL_TriggerET);
   fChain->SetBranchAddress("mu0_L0Calo_HCAL_TriggerHCALET", &mu0_L0Calo_HCAL_TriggerHCALET, &b_mu0_L0Calo_HCAL_TriggerHCALET);
   fChain->SetBranchAddress("mu0_L0Calo_HCAL_xTrigger", &mu0_L0Calo_HCAL_xTrigger, &b_mu0_L0Calo_HCAL_xTrigger);
   fChain->SetBranchAddress("mu0_L0Calo_HCAL_yTrigger", &mu0_L0Calo_HCAL_yTrigger, &b_mu0_L0Calo_HCAL_yTrigger);
   fChain->SetBranchAddress("mu0_ID", &mu0_ID, &b_mu0_ID);
   fChain->SetBranchAddress("mu0_CombDLLMu", &mu0_CombDLLMu, &b_mu0_CombDLLMu);
   fChain->SetBranchAddress("mu0_ProbNNmu", &mu0_ProbNNmu, &b_mu0_ProbNNmu);
   fChain->SetBranchAddress("mu0_ProbNNghost", &mu0_ProbNNghost, &b_mu0_ProbNNghost);
   fChain->SetBranchAddress("mu0_InMuonAcc", &mu0_InMuonAcc, &b_mu0_InMuonAcc);
   fChain->SetBranchAddress("mu0_MuonDist2", &mu0_MuonDist2, &b_mu0_MuonDist2);
   fChain->SetBranchAddress("mu0_regionInM2", &mu0_regionInM2, &b_mu0_regionInM2);
   fChain->SetBranchAddress("mu0_hasMuon", &mu0_hasMuon, &b_mu0_hasMuon);
   fChain->SetBranchAddress("mu0_isMuon", &mu0_isMuon, &b_mu0_isMuon);
   fChain->SetBranchAddress("mu0_isMuonLoose", &mu0_isMuonLoose, &b_mu0_isMuonLoose);
   fChain->SetBranchAddress("mu0_NShared", &mu0_NShared, &b_mu0_NShared);
   fChain->SetBranchAddress("mu0_MuonLLmu", &mu0_MuonLLmu, &b_mu0_MuonLLmu);
   fChain->SetBranchAddress("mu0_MuonLLbg", &mu0_MuonLLbg, &b_mu0_MuonLLbg);
   fChain->SetBranchAddress("mu0_isMuonFromProto", &mu0_isMuonFromProto, &b_mu0_isMuonFromProto);
   fChain->SetBranchAddress("mu0_PIDe", &mu0_PIDe, &b_mu0_PIDe);
   fChain->SetBranchAddress("mu0_PIDmu", &mu0_PIDmu, &b_mu0_PIDmu);
   fChain->SetBranchAddress("mu0_PIDK", &mu0_PIDK, &b_mu0_PIDK);
   fChain->SetBranchAddress("mu0_PIDp", &mu0_PIDp, &b_mu0_PIDp);
   fChain->SetBranchAddress("mu0_ProbNNe", &mu0_ProbNNe, &b_mu0_ProbNNe);
   fChain->SetBranchAddress("mu0_ProbNNk", &mu0_ProbNNk, &b_mu0_ProbNNk);
   fChain->SetBranchAddress("mu0_ProbNNp", &mu0_ProbNNp, &b_mu0_ProbNNp);
   fChain->SetBranchAddress("mu0_ProbNNpi", &mu0_ProbNNpi, &b_mu0_ProbNNpi);
   fChain->SetBranchAddress("mu0_hasRich", &mu0_hasRich, &b_mu0_hasRich);
   fChain->SetBranchAddress("mu0_hasCalo", &mu0_hasCalo, &b_mu0_hasCalo);
   fChain->SetBranchAddress("mu0_UsedRichAerogel", &mu0_UsedRichAerogel, &b_mu0_UsedRichAerogel);
   fChain->SetBranchAddress("mu0_UsedRich1Gas", &mu0_UsedRich1Gas, &b_mu0_UsedRich1Gas);
   fChain->SetBranchAddress("mu0_UsedRich2Gas", &mu0_UsedRich2Gas, &b_mu0_UsedRich2Gas);
   fChain->SetBranchAddress("mu0_RichAboveElThres", &mu0_RichAboveElThres, &b_mu0_RichAboveElThres);
   fChain->SetBranchAddress("mu0_RichAboveMuThres", &mu0_RichAboveMuThres, &b_mu0_RichAboveMuThres);
   fChain->SetBranchAddress("mu0_RichAbovePiThres", &mu0_RichAbovePiThres, &b_mu0_RichAbovePiThres);
   fChain->SetBranchAddress("mu0_RichAboveKaThres", &mu0_RichAboveKaThres, &b_mu0_RichAboveKaThres);
   fChain->SetBranchAddress("mu0_RichAbovePrThres", &mu0_RichAbovePrThres, &b_mu0_RichAbovePrThres);
   fChain->SetBranchAddress("mu0_RichDLLe", &mu0_RichDLLe, &b_mu0_RichDLLe);
   fChain->SetBranchAddress("mu0_RichDLLmu", &mu0_RichDLLmu, &b_mu0_RichDLLmu);
   fChain->SetBranchAddress("mu0_RichDLLpi", &mu0_RichDLLpi, &b_mu0_RichDLLpi);
   fChain->SetBranchAddress("mu0_RichDLLk", &mu0_RichDLLk, &b_mu0_RichDLLk);
   fChain->SetBranchAddress("mu0_RichDLLp", &mu0_RichDLLp, &b_mu0_RichDLLp);
   fChain->SetBranchAddress("mu0_RichDLLbt", &mu0_RichDLLbt, &b_mu0_RichDLLbt);
   fChain->SetBranchAddress("mu0_InAccMuon", &mu0_InAccMuon, &b_mu0_InAccMuon);
   fChain->SetBranchAddress("mu0_MuonMuLL", &mu0_MuonMuLL, &b_mu0_MuonMuLL);
   fChain->SetBranchAddress("mu0_MuonBkgLL", &mu0_MuonBkgLL, &b_mu0_MuonBkgLL);
   fChain->SetBranchAddress("mu0_MuonNShared", &mu0_MuonNShared, &b_mu0_MuonNShared);
   fChain->SetBranchAddress("mu0_InAccEcal", &mu0_InAccEcal, &b_mu0_InAccEcal);
   fChain->SetBranchAddress("mu0_CaloEcalE", &mu0_CaloEcalE, &b_mu0_CaloEcalE);
   fChain->SetBranchAddress("mu0_EcalPIDe", &mu0_EcalPIDe, &b_mu0_EcalPIDe);
   fChain->SetBranchAddress("mu0_EcalPIDmu", &mu0_EcalPIDmu, &b_mu0_EcalPIDmu);
   fChain->SetBranchAddress("mu0_InAccHcal", &mu0_InAccHcal, &b_mu0_InAccHcal);
   fChain->SetBranchAddress("mu0_CaloHcalE", &mu0_CaloHcalE, &b_mu0_CaloHcalE);
   fChain->SetBranchAddress("mu0_HcalPIDe", &mu0_HcalPIDe, &b_mu0_HcalPIDe);
   fChain->SetBranchAddress("mu0_HcalPIDmu", &mu0_HcalPIDmu, &b_mu0_HcalPIDmu);
   fChain->SetBranchAddress("mu0_InAccPrs", &mu0_InAccPrs, &b_mu0_InAccPrs);
   fChain->SetBranchAddress("mu0_PrsPIDe", &mu0_PrsPIDe, &b_mu0_PrsPIDe);
   fChain->SetBranchAddress("mu0_CaloPrsE", &mu0_CaloPrsE, &b_mu0_CaloPrsE);
   fChain->SetBranchAddress("mu0_InAccSpd", &mu0_InAccSpd, &b_mu0_InAccSpd);
   fChain->SetBranchAddress("mu0_CaloSpdE", &mu0_CaloSpdE, &b_mu0_CaloSpdE);
   fChain->SetBranchAddress("mu0_InAccBrem", &mu0_InAccBrem, &b_mu0_InAccBrem);
   fChain->SetBranchAddress("mu0_BremPIDe", &mu0_BremPIDe, &b_mu0_BremPIDe);
   fChain->SetBranchAddress("mu0_VeloCharge", &mu0_VeloCharge, &b_mu0_VeloCharge);
   fChain->SetBranchAddress("mu0_RICHDLLe", &mu0_RICHDLLe, &b_mu0_RICHDLLe);
   fChain->SetBranchAddress("mu0_RICHDLLmu", &mu0_RICHDLLmu, &b_mu0_RICHDLLmu);
   fChain->SetBranchAddress("mu0_RICHDLLpi", &mu0_RICHDLLpi, &b_mu0_RICHDLLpi);
   fChain->SetBranchAddress("mu0_RICHDLLK", &mu0_RICHDLLK, &b_mu0_RICHDLLK);
   fChain->SetBranchAddress("mu0_RICHDLLp", &mu0_RICHDLLp, &b_mu0_RICHDLLp);
   fChain->SetBranchAddress("mu0_RICHDLLbt", &mu0_RICHDLLbt, &b_mu0_RICHDLLbt);
   fChain->SetBranchAddress("mu0_RICHBestID", &mu0_RICHBestID, &b_mu0_RICHBestID);
   fChain->SetBranchAddress("mu0_RICHThreshold", &mu0_RICHThreshold, &b_mu0_RICHThreshold);
   fChain->SetBranchAddress("mu0_RICHThresholdEl", &mu0_RICHThresholdEl, &b_mu0_RICHThresholdEl);
   fChain->SetBranchAddress("mu0_RICHThresholdMu", &mu0_RICHThresholdMu, &b_mu0_RICHThresholdMu);
   fChain->SetBranchAddress("mu0_RICHThresholdPi", &mu0_RICHThresholdPi, &b_mu0_RICHThresholdPi);
   fChain->SetBranchAddress("mu0_RICHThresholdKa", &mu0_RICHThresholdKa, &b_mu0_RICHThresholdKa);
   fChain->SetBranchAddress("mu0_RICHThresholdPr", &mu0_RICHThresholdPr, &b_mu0_RICHThresholdPr);
   fChain->SetBranchAddress("mu0_RICHAerogelUsed", &mu0_RICHAerogelUsed, &b_mu0_RICHAerogelUsed);
   fChain->SetBranchAddress("mu0_RICH1GasUsed", &mu0_RICH1GasUsed, &b_mu0_RICH1GasUsed);
   fChain->SetBranchAddress("mu0_RICH2GasUsed", &mu0_RICH2GasUsed, &b_mu0_RICH2GasUsed);
   fChain->SetBranchAddress("mu0_TRACK_Eta", &mu0_TRACK_Eta, &b_mu0_TRACK_Eta);
   fChain->SetBranchAddress("mu0_TRACK_Phi", &mu0_TRACK_Phi, &b_mu0_TRACK_Phi);
   fChain->SetBranchAddress("mu0_Aerogel_X", &mu0_Aerogel_X, &b_mu0_Aerogel_X);
   fChain->SetBranchAddress("mu0_Aerogel_Y", &mu0_Aerogel_Y, &b_mu0_Aerogel_Y);
   fChain->SetBranchAddress("mu0_Aerogel_Z", &mu0_Aerogel_Z, &b_mu0_Aerogel_Z);
   fChain->SetBranchAddress("mu0_Aerogel_Rho", &mu0_Aerogel_Rho, &b_mu0_Aerogel_Rho);
   fChain->SetBranchAddress("mu0_Aerogel_Phi", &mu0_Aerogel_Phi, &b_mu0_Aerogel_Phi);
   fChain->SetBranchAddress("mu0_Rich1Gas_X", &mu0_Rich1Gas_X, &b_mu0_Rich1Gas_X);
   fChain->SetBranchAddress("mu0_Rich1Gas_Y", &mu0_Rich1Gas_Y, &b_mu0_Rich1Gas_Y);
   fChain->SetBranchAddress("mu0_Rich1Gas_Z", &mu0_Rich1Gas_Z, &b_mu0_Rich1Gas_Z);
   fChain->SetBranchAddress("mu0_Rich1Gas_Rho", &mu0_Rich1Gas_Rho, &b_mu0_Rich1Gas_Rho);
   fChain->SetBranchAddress("mu0_Rich1Gas_Phi", &mu0_Rich1Gas_Phi, &b_mu0_Rich1Gas_Phi);
   fChain->SetBranchAddress("mu0_Rich2Gas_X", &mu0_Rich2Gas_X, &b_mu0_Rich2Gas_X);
   fChain->SetBranchAddress("mu0_Rich2Gas_Y", &mu0_Rich2Gas_Y, &b_mu0_Rich2Gas_Y);
   fChain->SetBranchAddress("mu0_Rich2Gas_Z", &mu0_Rich2Gas_Z, &b_mu0_Rich2Gas_Z);
   fChain->SetBranchAddress("mu0_Rich2Gas_Rho", &mu0_Rich2Gas_Rho, &b_mu0_Rich2Gas_Rho);
   fChain->SetBranchAddress("mu0_Rich2Gas_Phi", &mu0_Rich2Gas_Phi, &b_mu0_Rich2Gas_Phi);
   fChain->SetBranchAddress("mu0_L0Global_Dec", &mu0_L0Global_Dec, &b_mu0_L0Global_Dec);
   fChain->SetBranchAddress("mu0_L0Global_TIS", &mu0_L0Global_TIS, &b_mu0_L0Global_TIS);
   fChain->SetBranchAddress("mu0_L0Global_TOS", &mu0_L0Global_TOS, &b_mu0_L0Global_TOS);
   fChain->SetBranchAddress("mu0_Hlt1Global_Dec", &mu0_Hlt1Global_Dec, &b_mu0_Hlt1Global_Dec);
   fChain->SetBranchAddress("mu0_Hlt1Global_TIS", &mu0_Hlt1Global_TIS, &b_mu0_Hlt1Global_TIS);
   fChain->SetBranchAddress("mu0_Hlt1Global_TOS", &mu0_Hlt1Global_TOS, &b_mu0_Hlt1Global_TOS);
   fChain->SetBranchAddress("mu0_Hlt1Phys_Dec", &mu0_Hlt1Phys_Dec, &b_mu0_Hlt1Phys_Dec);
   fChain->SetBranchAddress("mu0_Hlt1Phys_TIS", &mu0_Hlt1Phys_TIS, &b_mu0_Hlt1Phys_TIS);
   fChain->SetBranchAddress("mu0_Hlt1Phys_TOS", &mu0_Hlt1Phys_TOS, &b_mu0_Hlt1Phys_TOS);
   fChain->SetBranchAddress("mu0_Hlt2Global_Dec", &mu0_Hlt2Global_Dec, &b_mu0_Hlt2Global_Dec);
   fChain->SetBranchAddress("mu0_Hlt2Global_TIS", &mu0_Hlt2Global_TIS, &b_mu0_Hlt2Global_TIS);
   fChain->SetBranchAddress("mu0_Hlt2Global_TOS", &mu0_Hlt2Global_TOS, &b_mu0_Hlt2Global_TOS);
   fChain->SetBranchAddress("mu0_Hlt2Phys_Dec", &mu0_Hlt2Phys_Dec, &b_mu0_Hlt2Phys_Dec);
   fChain->SetBranchAddress("mu0_Hlt2Phys_TIS", &mu0_Hlt2Phys_TIS, &b_mu0_Hlt2Phys_TIS);
   fChain->SetBranchAddress("mu0_Hlt2Phys_TOS", &mu0_Hlt2Phys_TOS, &b_mu0_Hlt2Phys_TOS);
   fChain->SetBranchAddress("mu0_TRACK_Type", &mu0_TRACK_Type, &b_mu0_TRACK_Type);
   fChain->SetBranchAddress("mu0_TRACK_Key", &mu0_TRACK_Key, &b_mu0_TRACK_Key);
   fChain->SetBranchAddress("mu0_TRACK_CHI2", &mu0_TRACK_CHI2, &b_mu0_TRACK_CHI2);
   fChain->SetBranchAddress("mu0_TRACK_NDOF", &mu0_TRACK_NDOF, &b_mu0_TRACK_NDOF);
   fChain->SetBranchAddress("mu0_TRACK_CHI2NDOF", &mu0_TRACK_CHI2NDOF, &b_mu0_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("mu0_TRACK_PCHI2", &mu0_TRACK_PCHI2, &b_mu0_TRACK_PCHI2);
   fChain->SetBranchAddress("mu0_TRACK_VeloCHI2NDOF", &mu0_TRACK_VeloCHI2NDOF, &b_mu0_TRACK_VeloCHI2NDOF);
   fChain->SetBranchAddress("mu0_TRACK_TCHI2NDOF", &mu0_TRACK_TCHI2NDOF, &b_mu0_TRACK_TCHI2NDOF);
   fChain->SetBranchAddress("mu0_VELO_UTID", &mu0_VELO_UTID, &b_mu0_VELO_UTID);
   fChain->SetBranchAddress("mu0_TRACK_FirstMeasurementX", &mu0_TRACK_FirstMeasurementX, &b_mu0_TRACK_FirstMeasurementX);
   fChain->SetBranchAddress("mu0_TRACK_FirstMeasurementY", &mu0_TRACK_FirstMeasurementY, &b_mu0_TRACK_FirstMeasurementY);
   fChain->SetBranchAddress("mu0_TRACK_FirstMeasurementZ", &mu0_TRACK_FirstMeasurementZ, &b_mu0_TRACK_FirstMeasurementZ);
   fChain->SetBranchAddress("mu0_TRACK_MatchCHI2", &mu0_TRACK_MatchCHI2, &b_mu0_TRACK_MatchCHI2);
   fChain->SetBranchAddress("mu0_TRACK_GhostProb", &mu0_TRACK_GhostProb, &b_mu0_TRACK_GhostProb);
   fChain->SetBranchAddress("mu0_TRACK_CloneDist", &mu0_TRACK_CloneDist, &b_mu0_TRACK_CloneDist);
   fChain->SetBranchAddress("mu0_TRACK_Likelihood", &mu0_TRACK_Likelihood, &b_mu0_TRACK_Likelihood);
   fChain->SetBranchAddress("mu0_L0HadronDecision_Dec", &mu0_L0HadronDecision_Dec, &b_mu0_L0HadronDecision_Dec);
   fChain->SetBranchAddress("mu0_L0HadronDecision_TIS", &mu0_L0HadronDecision_TIS, &b_mu0_L0HadronDecision_TIS);
   fChain->SetBranchAddress("mu0_L0HadronDecision_TOS", &mu0_L0HadronDecision_TOS, &b_mu0_L0HadronDecision_TOS);
   fChain->SetBranchAddress("mu0_L0MuonDecision_Dec", &mu0_L0MuonDecision_Dec, &b_mu0_L0MuonDecision_Dec);
   fChain->SetBranchAddress("mu0_L0MuonDecision_TIS", &mu0_L0MuonDecision_TIS, &b_mu0_L0MuonDecision_TIS);
   fChain->SetBranchAddress("mu0_L0MuonDecision_TOS", &mu0_L0MuonDecision_TOS, &b_mu0_L0MuonDecision_TOS);
   fChain->SetBranchAddress("mu0_L0DiMuonDecision_Dec", &mu0_L0DiMuonDecision_Dec, &b_mu0_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("mu0_L0DiMuonDecision_TIS", &mu0_L0DiMuonDecision_TIS, &b_mu0_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("mu0_L0DiMuonDecision_TOS", &mu0_L0DiMuonDecision_TOS, &b_mu0_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("mu0_L0ElectronDecision_Dec", &mu0_L0ElectronDecision_Dec, &b_mu0_L0ElectronDecision_Dec);
   fChain->SetBranchAddress("mu0_L0ElectronDecision_TIS", &mu0_L0ElectronDecision_TIS, &b_mu0_L0ElectronDecision_TIS);
   fChain->SetBranchAddress("mu0_L0ElectronDecision_TOS", &mu0_L0ElectronDecision_TOS, &b_mu0_L0ElectronDecision_TOS);
   fChain->SetBranchAddress("mu0_L0PhotonDecision_Dec", &mu0_L0PhotonDecision_Dec, &b_mu0_L0PhotonDecision_Dec);
   fChain->SetBranchAddress("mu0_L0PhotonDecision_TIS", &mu0_L0PhotonDecision_TIS, &b_mu0_L0PhotonDecision_TIS);
   fChain->SetBranchAddress("mu0_L0PhotonDecision_TOS", &mu0_L0PhotonDecision_TOS, &b_mu0_L0PhotonDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt1DiMuonHighMassDecision_Dec", &mu0_Hlt1DiMuonHighMassDecision_Dec, &b_mu0_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt1DiMuonHighMassDecision_TIS", &mu0_Hlt1DiMuonHighMassDecision_TIS, &b_mu0_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt1DiMuonHighMassDecision_TOS", &mu0_Hlt1DiMuonHighMassDecision_TOS, &b_mu0_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt1DiMuonLowMassDecision_Dec", &mu0_Hlt1DiMuonLowMassDecision_Dec, &b_mu0_Hlt1DiMuonLowMassDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt1DiMuonLowMassDecision_TIS", &mu0_Hlt1DiMuonLowMassDecision_TIS, &b_mu0_Hlt1DiMuonLowMassDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt1DiMuonLowMassDecision_TOS", &mu0_Hlt1DiMuonLowMassDecision_TOS, &b_mu0_Hlt1DiMuonLowMassDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt1SingleMuonNoIPDecision_Dec", &mu0_Hlt1SingleMuonNoIPDecision_Dec, &b_mu0_Hlt1SingleMuonNoIPDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt1SingleMuonNoIPDecision_TIS", &mu0_Hlt1SingleMuonNoIPDecision_TIS, &b_mu0_Hlt1SingleMuonNoIPDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt1SingleMuonNoIPDecision_TOS", &mu0_Hlt1SingleMuonNoIPDecision_TOS, &b_mu0_Hlt1SingleMuonNoIPDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt1SingleMuonHighPTDecision_Dec", &mu0_Hlt1SingleMuonHighPTDecision_Dec, &b_mu0_Hlt1SingleMuonHighPTDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt1SingleMuonHighPTDecision_TIS", &mu0_Hlt1SingleMuonHighPTDecision_TIS, &b_mu0_Hlt1SingleMuonHighPTDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt1SingleMuonHighPTDecision_TOS", &mu0_Hlt1SingleMuonHighPTDecision_TOS, &b_mu0_Hlt1SingleMuonHighPTDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt1TrackAllL0Decision_Dec", &mu0_Hlt1TrackAllL0Decision_Dec, &b_mu0_Hlt1TrackAllL0Decision_Dec);
   fChain->SetBranchAddress("mu0_Hlt1TrackAllL0Decision_TIS", &mu0_Hlt1TrackAllL0Decision_TIS, &b_mu0_Hlt1TrackAllL0Decision_TIS);
   fChain->SetBranchAddress("mu0_Hlt1TrackAllL0Decision_TOS", &mu0_Hlt1TrackAllL0Decision_TOS, &b_mu0_Hlt1TrackAllL0Decision_TOS);
   fChain->SetBranchAddress("mu0_Hlt1TrackMuonDecision_Dec", &mu0_Hlt1TrackMuonDecision_Dec, &b_mu0_Hlt1TrackMuonDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TIS", &mu0_Hlt1TrackMuonDecision_TIS, &b_mu0_Hlt1TrackMuonDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS", &mu0_Hlt1TrackMuonDecision_TOS, &b_mu0_Hlt1TrackMuonDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt1TrackPhotonDecision_Dec", &mu0_Hlt1TrackPhotonDecision_Dec, &b_mu0_Hlt1TrackPhotonDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt1TrackPhotonDecision_TIS", &mu0_Hlt1TrackPhotonDecision_TIS, &b_mu0_Hlt1TrackPhotonDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt1TrackPhotonDecision_TOS", &mu0_Hlt1TrackPhotonDecision_TOS, &b_mu0_Hlt1TrackPhotonDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt1L0AnyDecision_Dec", &mu0_Hlt1L0AnyDecision_Dec, &b_mu0_Hlt1L0AnyDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt1L0AnyDecision_TIS", &mu0_Hlt1L0AnyDecision_TIS, &b_mu0_Hlt1L0AnyDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt1L0AnyDecision_TOS", &mu0_Hlt1L0AnyDecision_TOS, &b_mu0_Hlt1L0AnyDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt1GlobalDecision_Dec", &mu0_Hlt1GlobalDecision_Dec, &b_mu0_Hlt1GlobalDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt1GlobalDecision_TIS", &mu0_Hlt1GlobalDecision_TIS, &b_mu0_Hlt1GlobalDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt1GlobalDecision_TOS", &mu0_Hlt1GlobalDecision_TOS, &b_mu0_Hlt1GlobalDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2SingleMuonDecision_Dec", &mu0_Hlt2SingleMuonDecision_Dec, &b_mu0_Hlt2SingleMuonDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2SingleMuonDecision_TIS", &mu0_Hlt2SingleMuonDecision_TIS, &b_mu0_Hlt2SingleMuonDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2SingleMuonDecision_TOS", &mu0_Hlt2SingleMuonDecision_TOS, &b_mu0_Hlt2SingleMuonDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2DiMuonDetachedDecision_Dec", &mu0_Hlt2DiMuonDetachedDecision_Dec, &b_mu0_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2DiMuonDetachedDecision_TIS", &mu0_Hlt2DiMuonDetachedDecision_TIS, &b_mu0_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2DiMuonDetachedDecision_TOS", &mu0_Hlt2DiMuonDetachedDecision_TOS, &b_mu0_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD2HMuMuDecision_Dec", &mu0_Hlt2CharmSemilepD2HMuMuDecision_Dec, &b_mu0_Hlt2CharmSemilepD2HMuMuDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD2HMuMuDecision_TIS", &mu0_Hlt2CharmSemilepD2HMuMuDecision_TIS, &b_mu0_Hlt2CharmSemilepD2HMuMuDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD2HMuMuDecision_TOS", &mu0_Hlt2CharmSemilepD2HMuMuDecision_TOS, &b_mu0_Hlt2CharmSemilepD2HMuMuDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec", &mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec, &b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS", &mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS, &b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS", &mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS, &b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec", &mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec, &b_mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS", &mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS, &b_mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS", &mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS, &b_mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec", &mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec, &b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS", &mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS, &b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS", &mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS, &b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec", &mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec, &b_mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS", &mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS, &b_mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS", &mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS, &b_mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec", &mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec, &b_mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS", &mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS, &b_mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS", &mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS, &b_mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD02KKMuMuDecision_Dec", &mu0_Hlt2CharmSemilepD02KKMuMuDecision_Dec, &b_mu0_Hlt2CharmSemilepD02KKMuMuDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD02KKMuMuDecision_TIS", &mu0_Hlt2CharmSemilepD02KKMuMuDecision_TIS, &b_mu0_Hlt2CharmSemilepD02KKMuMuDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD02KKMuMuDecision_TOS", &mu0_Hlt2CharmSemilepD02KKMuMuDecision_TOS, &b_mu0_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD02KPiMuMuDecision_Dec", &mu0_Hlt2CharmSemilepD02KPiMuMuDecision_Dec, &b_mu0_Hlt2CharmSemilepD02KPiMuMuDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD02KPiMuMuDecision_TIS", &mu0_Hlt2CharmSemilepD02KPiMuMuDecision_TIS, &b_mu0_Hlt2CharmSemilepD02KPiMuMuDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmSemilepD02KPiMuMuDecision_TOS", &mu0_Hlt2CharmSemilepD02KPiMuMuDecision_TOS, &b_mu0_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHHDst_4piDecision_Dec", &mu0_Hlt2CharmHadD02HHHHDst_4piDecision_Dec, &b_mu0_Hlt2CharmHadD02HHHHDst_4piDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHHDst_4piDecision_TIS", &mu0_Hlt2CharmHadD02HHHHDst_4piDecision_TIS, &b_mu0_Hlt2CharmHadD02HHHHDst_4piDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHHDst_4piDecision_TOS", &mu0_Hlt2CharmHadD02HHHHDst_4piDecision_TOS, &b_mu0_Hlt2CharmHadD02HHHHDst_4piDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec", &mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec, &b_mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS", &mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS, &b_mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS", &mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS, &b_mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec", &mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec, &b_mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS", &mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS, &b_mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS", &mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS, &b_mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_hhXDecision_Dec", &mu0_Hlt2CharmHadD02HHXDst_hhXDecision_Dec, &b_mu0_Hlt2CharmHadD02HHXDst_hhXDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_hhXDecision_TIS", &mu0_Hlt2CharmHadD02HHXDst_hhXDecision_TIS, &b_mu0_Hlt2CharmHadD02HHXDst_hhXDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_hhXDecision_TOS", &mu0_Hlt2CharmHadD02HHXDst_hhXDecision_TOS, &b_mu0_Hlt2CharmHadD02HHXDst_hhXDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec", &mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec, &b_mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS", &mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS, &b_mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS", &mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS, &b_mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec", &mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec, &b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS", &mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS, &b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS", &mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS, &b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec", &mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec, &b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS", &mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS, &b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS", &mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS, &b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec", &mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec, &b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS", &mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS, &b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS", &mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS, &b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec", &mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec, &b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS", &mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS, &b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS", &mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS, &b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_K3piDecision_Dec", &mu0_Hlt2CharmHadD02HHHH_K3piDecision_Dec, &b_mu0_Hlt2CharmHadD02HHHH_K3piDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_K3piDecision_TIS", &mu0_Hlt2CharmHadD02HHHH_K3piDecision_TIS, &b_mu0_Hlt2CharmHadD02HHHH_K3piDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_K3piDecision_TOS", &mu0_Hlt2CharmHadD02HHHH_K3piDecision_TOS, &b_mu0_Hlt2CharmHadD02HHHH_K3piDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec", &mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec, &b_mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS", &mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS, &b_mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS", &mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS, &b_mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec", &mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec, &b_mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS", &mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS, &b_mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS", &mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS, &b_mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec", &mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec, &b_mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS", &mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS, &b_mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS", &mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS, &b_mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_4piDecision_Dec", &mu0_Hlt2CharmHadD02HHHH_4piDecision_Dec, &b_mu0_Hlt2CharmHadD02HHHH_4piDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_4piDecision_TIS", &mu0_Hlt2CharmHadD02HHHH_4piDecision_TIS, &b_mu0_Hlt2CharmHadD02HHHH_4piDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_4piDecision_TOS", &mu0_Hlt2CharmHadD02HHHH_4piDecision_TOS, &b_mu0_Hlt2CharmHadD02HHHH_4piDecision_TOS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec", &mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec, &b_mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS", &mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS, &b_mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS);
   fChain->SetBranchAddress("mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS", &mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS, &b_mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS);
   fChain->SetBranchAddress("mu1_MINIP", &mu1_MINIP, &b_mu1_MINIP);
   fChain->SetBranchAddress("mu1_MINIPCHI2", &mu1_MINIPCHI2, &b_mu1_MINIPCHI2);
   fChain->SetBranchAddress("mu1_MINIPNEXTBEST", &mu1_MINIPNEXTBEST, &b_mu1_MINIPNEXTBEST);
   fChain->SetBranchAddress("mu1_MINIPCHI2NEXTBEST", &mu1_MINIPCHI2NEXTBEST, &b_mu1_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("mu1_OWNPV_X", &mu1_OWNPV_X, &b_mu1_OWNPV_X);
   fChain->SetBranchAddress("mu1_OWNPV_Y", &mu1_OWNPV_Y, &b_mu1_OWNPV_Y);
   fChain->SetBranchAddress("mu1_OWNPV_Z", &mu1_OWNPV_Z, &b_mu1_OWNPV_Z);
   fChain->SetBranchAddress("mu1_OWNPV_XERR", &mu1_OWNPV_XERR, &b_mu1_OWNPV_XERR);
   fChain->SetBranchAddress("mu1_OWNPV_YERR", &mu1_OWNPV_YERR, &b_mu1_OWNPV_YERR);
   fChain->SetBranchAddress("mu1_OWNPV_ZERR", &mu1_OWNPV_ZERR, &b_mu1_OWNPV_ZERR);
   fChain->SetBranchAddress("mu1_OWNPV_CHI2", &mu1_OWNPV_CHI2, &b_mu1_OWNPV_CHI2);
   fChain->SetBranchAddress("mu1_OWNPV_NDOF", &mu1_OWNPV_NDOF, &b_mu1_OWNPV_NDOF);
   fChain->SetBranchAddress("mu1_OWNPV_COV_", mu1_OWNPV_COV_, &b_mu1_OWNPV_COV_);
   fChain->SetBranchAddress("mu1_IP_OWNPV", &mu1_IP_OWNPV, &b_mu1_IP_OWNPV);
   fChain->SetBranchAddress("mu1_IPCHI2_OWNPV", &mu1_IPCHI2_OWNPV, &b_mu1_IPCHI2_OWNPV);
   fChain->SetBranchAddress("mu1_TOPPV_X", &mu1_TOPPV_X, &b_mu1_TOPPV_X);
   fChain->SetBranchAddress("mu1_TOPPV_Y", &mu1_TOPPV_Y, &b_mu1_TOPPV_Y);
   fChain->SetBranchAddress("mu1_TOPPV_Z", &mu1_TOPPV_Z, &b_mu1_TOPPV_Z);
   fChain->SetBranchAddress("mu1_TOPPV_XERR", &mu1_TOPPV_XERR, &b_mu1_TOPPV_XERR);
   fChain->SetBranchAddress("mu1_TOPPV_YERR", &mu1_TOPPV_YERR, &b_mu1_TOPPV_YERR);
   fChain->SetBranchAddress("mu1_TOPPV_ZERR", &mu1_TOPPV_ZERR, &b_mu1_TOPPV_ZERR);
   fChain->SetBranchAddress("mu1_TOPPV_CHI2", &mu1_TOPPV_CHI2, &b_mu1_TOPPV_CHI2);
   fChain->SetBranchAddress("mu1_TOPPV_NDOF", &mu1_TOPPV_NDOF, &b_mu1_TOPPV_NDOF);
   fChain->SetBranchAddress("mu1_TOPPV_COV_", mu1_TOPPV_COV_, &b_mu1_TOPPV_COV_);
   fChain->SetBranchAddress("mu1_IP_TOPPV", &mu1_IP_TOPPV, &b_mu1_IP_TOPPV);
   fChain->SetBranchAddress("mu1_IPCHI2_TOPPV", &mu1_IPCHI2_TOPPV, &b_mu1_IPCHI2_TOPPV);
   fChain->SetBranchAddress("mu1_ORIVX_X", &mu1_ORIVX_X, &b_mu1_ORIVX_X);
   fChain->SetBranchAddress("mu1_ORIVX_Y", &mu1_ORIVX_Y, &b_mu1_ORIVX_Y);
   fChain->SetBranchAddress("mu1_ORIVX_Z", &mu1_ORIVX_Z, &b_mu1_ORIVX_Z);
   fChain->SetBranchAddress("mu1_ORIVX_XERR", &mu1_ORIVX_XERR, &b_mu1_ORIVX_XERR);
   fChain->SetBranchAddress("mu1_ORIVX_YERR", &mu1_ORIVX_YERR, &b_mu1_ORIVX_YERR);
   fChain->SetBranchAddress("mu1_ORIVX_ZERR", &mu1_ORIVX_ZERR, &b_mu1_ORIVX_ZERR);
   fChain->SetBranchAddress("mu1_ORIVX_CHI2", &mu1_ORIVX_CHI2, &b_mu1_ORIVX_CHI2);
   fChain->SetBranchAddress("mu1_ORIVX_NDOF", &mu1_ORIVX_NDOF, &b_mu1_ORIVX_NDOF);
   fChain->SetBranchAddress("mu1_ORIVX_COV_", mu1_ORIVX_COV_, &b_mu1_ORIVX_COV_);
   fChain->SetBranchAddress("mu1_IP_ORIVX", &mu1_IP_ORIVX, &b_mu1_IP_ORIVX);
   fChain->SetBranchAddress("mu1_IPCHI2_ORIVX", &mu1_IPCHI2_ORIVX, &b_mu1_IPCHI2_ORIVX);
   fChain->SetBranchAddress("mu1_P", &mu1_P, &b_mu1_P);
   fChain->SetBranchAddress("mu1_PT", &mu1_PT, &b_mu1_PT);
   fChain->SetBranchAddress("mu1_PE", &mu1_PE, &b_mu1_PE);
   fChain->SetBranchAddress("mu1_PX", &mu1_PX, &b_mu1_PX);
   fChain->SetBranchAddress("mu1_PY", &mu1_PY, &b_mu1_PY);
   fChain->SetBranchAddress("mu1_PZ", &mu1_PZ, &b_mu1_PZ);
   fChain->SetBranchAddress("mu1_M", &mu1_M, &b_mu1_M);
   fChain->SetBranchAddress("mu1_L0Calo_HCAL_realET", &mu1_L0Calo_HCAL_realET, &b_mu1_L0Calo_HCAL_realET);
   fChain->SetBranchAddress("mu1_L0Calo_HCAL_xProjection", &mu1_L0Calo_HCAL_xProjection, &b_mu1_L0Calo_HCAL_xProjection);
   fChain->SetBranchAddress("mu1_L0Calo_HCAL_yProjection", &mu1_L0Calo_HCAL_yProjection, &b_mu1_L0Calo_HCAL_yProjection);
   fChain->SetBranchAddress("mu1_L0Calo_HCAL_region", &mu1_L0Calo_HCAL_region, &b_mu1_L0Calo_HCAL_region);
   fChain->SetBranchAddress("mu1_L0Calo_HCAL_TriggerET", &mu1_L0Calo_HCAL_TriggerET, &b_mu1_L0Calo_HCAL_TriggerET);
   fChain->SetBranchAddress("mu1_L0Calo_HCAL_TriggerHCALET", &mu1_L0Calo_HCAL_TriggerHCALET, &b_mu1_L0Calo_HCAL_TriggerHCALET);
   fChain->SetBranchAddress("mu1_L0Calo_HCAL_xTrigger", &mu1_L0Calo_HCAL_xTrigger, &b_mu1_L0Calo_HCAL_xTrigger);
   fChain->SetBranchAddress("mu1_L0Calo_HCAL_yTrigger", &mu1_L0Calo_HCAL_yTrigger, &b_mu1_L0Calo_HCAL_yTrigger);
   fChain->SetBranchAddress("mu1_ID", &mu1_ID, &b_mu1_ID);
   fChain->SetBranchAddress("mu1_CombDLLMu", &mu1_CombDLLMu, &b_mu1_CombDLLMu);
   fChain->SetBranchAddress("mu1_ProbNNmu", &mu1_ProbNNmu, &b_mu1_ProbNNmu);
   fChain->SetBranchAddress("mu1_ProbNNghost", &mu1_ProbNNghost, &b_mu1_ProbNNghost);
   fChain->SetBranchAddress("mu1_InMuonAcc", &mu1_InMuonAcc, &b_mu1_InMuonAcc);
   fChain->SetBranchAddress("mu1_MuonDist2", &mu1_MuonDist2, &b_mu1_MuonDist2);
   fChain->SetBranchAddress("mu1_regionInM2", &mu1_regionInM2, &b_mu1_regionInM2);
   fChain->SetBranchAddress("mu1_hasMuon", &mu1_hasMuon, &b_mu1_hasMuon);
   fChain->SetBranchAddress("mu1_isMuon", &mu1_isMuon, &b_mu1_isMuon);
   fChain->SetBranchAddress("mu1_isMuonLoose", &mu1_isMuonLoose, &b_mu1_isMuonLoose);
   fChain->SetBranchAddress("mu1_NShared", &mu1_NShared, &b_mu1_NShared);
   fChain->SetBranchAddress("mu1_MuonLLmu", &mu1_MuonLLmu, &b_mu1_MuonLLmu);
   fChain->SetBranchAddress("mu1_MuonLLbg", &mu1_MuonLLbg, &b_mu1_MuonLLbg);
   fChain->SetBranchAddress("mu1_isMuonFromProto", &mu1_isMuonFromProto, &b_mu1_isMuonFromProto);
   fChain->SetBranchAddress("mu1_PIDe", &mu1_PIDe, &b_mu1_PIDe);
   fChain->SetBranchAddress("mu1_PIDmu", &mu1_PIDmu, &b_mu1_PIDmu);
   fChain->SetBranchAddress("mu1_PIDK", &mu1_PIDK, &b_mu1_PIDK);
   fChain->SetBranchAddress("mu1_PIDp", &mu1_PIDp, &b_mu1_PIDp);
   fChain->SetBranchAddress("mu1_ProbNNe", &mu1_ProbNNe, &b_mu1_ProbNNe);
   fChain->SetBranchAddress("mu1_ProbNNk", &mu1_ProbNNk, &b_mu1_ProbNNk);
   fChain->SetBranchAddress("mu1_ProbNNp", &mu1_ProbNNp, &b_mu1_ProbNNp);
   fChain->SetBranchAddress("mu1_ProbNNpi", &mu1_ProbNNpi, &b_mu1_ProbNNpi);
   fChain->SetBranchAddress("mu1_hasRich", &mu1_hasRich, &b_mu1_hasRich);
   fChain->SetBranchAddress("mu1_hasCalo", &mu1_hasCalo, &b_mu1_hasCalo);
   fChain->SetBranchAddress("mu1_UsedRichAerogel", &mu1_UsedRichAerogel, &b_mu1_UsedRichAerogel);
   fChain->SetBranchAddress("mu1_UsedRich1Gas", &mu1_UsedRich1Gas, &b_mu1_UsedRich1Gas);
   fChain->SetBranchAddress("mu1_UsedRich2Gas", &mu1_UsedRich2Gas, &b_mu1_UsedRich2Gas);
   fChain->SetBranchAddress("mu1_RichAboveElThres", &mu1_RichAboveElThres, &b_mu1_RichAboveElThres);
   fChain->SetBranchAddress("mu1_RichAboveMuThres", &mu1_RichAboveMuThres, &b_mu1_RichAboveMuThres);
   fChain->SetBranchAddress("mu1_RichAbovePiThres", &mu1_RichAbovePiThres, &b_mu1_RichAbovePiThres);
   fChain->SetBranchAddress("mu1_RichAboveKaThres", &mu1_RichAboveKaThres, &b_mu1_RichAboveKaThres);
   fChain->SetBranchAddress("mu1_RichAbovePrThres", &mu1_RichAbovePrThres, &b_mu1_RichAbovePrThres);
   fChain->SetBranchAddress("mu1_RichDLLe", &mu1_RichDLLe, &b_mu1_RichDLLe);
   fChain->SetBranchAddress("mu1_RichDLLmu", &mu1_RichDLLmu, &b_mu1_RichDLLmu);
   fChain->SetBranchAddress("mu1_RichDLLpi", &mu1_RichDLLpi, &b_mu1_RichDLLpi);
   fChain->SetBranchAddress("mu1_RichDLLk", &mu1_RichDLLk, &b_mu1_RichDLLk);
   fChain->SetBranchAddress("mu1_RichDLLp", &mu1_RichDLLp, &b_mu1_RichDLLp);
   fChain->SetBranchAddress("mu1_RichDLLbt", &mu1_RichDLLbt, &b_mu1_RichDLLbt);
   fChain->SetBranchAddress("mu1_InAccMuon", &mu1_InAccMuon, &b_mu1_InAccMuon);
   fChain->SetBranchAddress("mu1_MuonMuLL", &mu1_MuonMuLL, &b_mu1_MuonMuLL);
   fChain->SetBranchAddress("mu1_MuonBkgLL", &mu1_MuonBkgLL, &b_mu1_MuonBkgLL);
   fChain->SetBranchAddress("mu1_MuonNShared", &mu1_MuonNShared, &b_mu1_MuonNShared);
   fChain->SetBranchAddress("mu1_InAccEcal", &mu1_InAccEcal, &b_mu1_InAccEcal);
   fChain->SetBranchAddress("mu1_CaloEcalE", &mu1_CaloEcalE, &b_mu1_CaloEcalE);
   fChain->SetBranchAddress("mu1_EcalPIDe", &mu1_EcalPIDe, &b_mu1_EcalPIDe);
   fChain->SetBranchAddress("mu1_EcalPIDmu", &mu1_EcalPIDmu, &b_mu1_EcalPIDmu);
   fChain->SetBranchAddress("mu1_InAccHcal", &mu1_InAccHcal, &b_mu1_InAccHcal);
   fChain->SetBranchAddress("mu1_CaloHcalE", &mu1_CaloHcalE, &b_mu1_CaloHcalE);
   fChain->SetBranchAddress("mu1_HcalPIDe", &mu1_HcalPIDe, &b_mu1_HcalPIDe);
   fChain->SetBranchAddress("mu1_HcalPIDmu", &mu1_HcalPIDmu, &b_mu1_HcalPIDmu);
   fChain->SetBranchAddress("mu1_InAccPrs", &mu1_InAccPrs, &b_mu1_InAccPrs);
   fChain->SetBranchAddress("mu1_PrsPIDe", &mu1_PrsPIDe, &b_mu1_PrsPIDe);
   fChain->SetBranchAddress("mu1_CaloPrsE", &mu1_CaloPrsE, &b_mu1_CaloPrsE);
   fChain->SetBranchAddress("mu1_InAccSpd", &mu1_InAccSpd, &b_mu1_InAccSpd);
   fChain->SetBranchAddress("mu1_CaloSpdE", &mu1_CaloSpdE, &b_mu1_CaloSpdE);
   fChain->SetBranchAddress("mu1_InAccBrem", &mu1_InAccBrem, &b_mu1_InAccBrem);
   fChain->SetBranchAddress("mu1_BremPIDe", &mu1_BremPIDe, &b_mu1_BremPIDe);
   fChain->SetBranchAddress("mu1_VeloCharge", &mu1_VeloCharge, &b_mu1_VeloCharge);
   fChain->SetBranchAddress("mu1_RICHDLLe", &mu1_RICHDLLe, &b_mu1_RICHDLLe);
   fChain->SetBranchAddress("mu1_RICHDLLmu", &mu1_RICHDLLmu, &b_mu1_RICHDLLmu);
   fChain->SetBranchAddress("mu1_RICHDLLpi", &mu1_RICHDLLpi, &b_mu1_RICHDLLpi);
   fChain->SetBranchAddress("mu1_RICHDLLK", &mu1_RICHDLLK, &b_mu1_RICHDLLK);
   fChain->SetBranchAddress("mu1_RICHDLLp", &mu1_RICHDLLp, &b_mu1_RICHDLLp);
   fChain->SetBranchAddress("mu1_RICHDLLbt", &mu1_RICHDLLbt, &b_mu1_RICHDLLbt);
   fChain->SetBranchAddress("mu1_RICHBestID", &mu1_RICHBestID, &b_mu1_RICHBestID);
   fChain->SetBranchAddress("mu1_RICHThreshold", &mu1_RICHThreshold, &b_mu1_RICHThreshold);
   fChain->SetBranchAddress("mu1_RICHThresholdEl", &mu1_RICHThresholdEl, &b_mu1_RICHThresholdEl);
   fChain->SetBranchAddress("mu1_RICHThresholdMu", &mu1_RICHThresholdMu, &b_mu1_RICHThresholdMu);
   fChain->SetBranchAddress("mu1_RICHThresholdPi", &mu1_RICHThresholdPi, &b_mu1_RICHThresholdPi);
   fChain->SetBranchAddress("mu1_RICHThresholdKa", &mu1_RICHThresholdKa, &b_mu1_RICHThresholdKa);
   fChain->SetBranchAddress("mu1_RICHThresholdPr", &mu1_RICHThresholdPr, &b_mu1_RICHThresholdPr);
   fChain->SetBranchAddress("mu1_RICHAerogelUsed", &mu1_RICHAerogelUsed, &b_mu1_RICHAerogelUsed);
   fChain->SetBranchAddress("mu1_RICH1GasUsed", &mu1_RICH1GasUsed, &b_mu1_RICH1GasUsed);
   fChain->SetBranchAddress("mu1_RICH2GasUsed", &mu1_RICH2GasUsed, &b_mu1_RICH2GasUsed);
   fChain->SetBranchAddress("mu1_TRACK_Eta", &mu1_TRACK_Eta, &b_mu1_TRACK_Eta);
   fChain->SetBranchAddress("mu1_TRACK_Phi", &mu1_TRACK_Phi, &b_mu1_TRACK_Phi);
   fChain->SetBranchAddress("mu1_Aerogel_X", &mu1_Aerogel_X, &b_mu1_Aerogel_X);
   fChain->SetBranchAddress("mu1_Aerogel_Y", &mu1_Aerogel_Y, &b_mu1_Aerogel_Y);
   fChain->SetBranchAddress("mu1_Aerogel_Z", &mu1_Aerogel_Z, &b_mu1_Aerogel_Z);
   fChain->SetBranchAddress("mu1_Aerogel_Rho", &mu1_Aerogel_Rho, &b_mu1_Aerogel_Rho);
   fChain->SetBranchAddress("mu1_Aerogel_Phi", &mu1_Aerogel_Phi, &b_mu1_Aerogel_Phi);
   fChain->SetBranchAddress("mu1_Rich1Gas_X", &mu1_Rich1Gas_X, &b_mu1_Rich1Gas_X);
   fChain->SetBranchAddress("mu1_Rich1Gas_Y", &mu1_Rich1Gas_Y, &b_mu1_Rich1Gas_Y);
   fChain->SetBranchAddress("mu1_Rich1Gas_Z", &mu1_Rich1Gas_Z, &b_mu1_Rich1Gas_Z);
   fChain->SetBranchAddress("mu1_Rich1Gas_Rho", &mu1_Rich1Gas_Rho, &b_mu1_Rich1Gas_Rho);
   fChain->SetBranchAddress("mu1_Rich1Gas_Phi", &mu1_Rich1Gas_Phi, &b_mu1_Rich1Gas_Phi);
   fChain->SetBranchAddress("mu1_Rich2Gas_X", &mu1_Rich2Gas_X, &b_mu1_Rich2Gas_X);
   fChain->SetBranchAddress("mu1_Rich2Gas_Y", &mu1_Rich2Gas_Y, &b_mu1_Rich2Gas_Y);
   fChain->SetBranchAddress("mu1_Rich2Gas_Z", &mu1_Rich2Gas_Z, &b_mu1_Rich2Gas_Z);
   fChain->SetBranchAddress("mu1_Rich2Gas_Rho", &mu1_Rich2Gas_Rho, &b_mu1_Rich2Gas_Rho);
   fChain->SetBranchAddress("mu1_Rich2Gas_Phi", &mu1_Rich2Gas_Phi, &b_mu1_Rich2Gas_Phi);
   fChain->SetBranchAddress("mu1_L0Global_Dec", &mu1_L0Global_Dec, &b_mu1_L0Global_Dec);
   fChain->SetBranchAddress("mu1_L0Global_TIS", &mu1_L0Global_TIS, &b_mu1_L0Global_TIS);
   fChain->SetBranchAddress("mu1_L0Global_TOS", &mu1_L0Global_TOS, &b_mu1_L0Global_TOS);
   fChain->SetBranchAddress("mu1_Hlt1Global_Dec", &mu1_Hlt1Global_Dec, &b_mu1_Hlt1Global_Dec);
   fChain->SetBranchAddress("mu1_Hlt1Global_TIS", &mu1_Hlt1Global_TIS, &b_mu1_Hlt1Global_TIS);
   fChain->SetBranchAddress("mu1_Hlt1Global_TOS", &mu1_Hlt1Global_TOS, &b_mu1_Hlt1Global_TOS);
   fChain->SetBranchAddress("mu1_Hlt1Phys_Dec", &mu1_Hlt1Phys_Dec, &b_mu1_Hlt1Phys_Dec);
   fChain->SetBranchAddress("mu1_Hlt1Phys_TIS", &mu1_Hlt1Phys_TIS, &b_mu1_Hlt1Phys_TIS);
   fChain->SetBranchAddress("mu1_Hlt1Phys_TOS", &mu1_Hlt1Phys_TOS, &b_mu1_Hlt1Phys_TOS);
   fChain->SetBranchAddress("mu1_Hlt2Global_Dec", &mu1_Hlt2Global_Dec, &b_mu1_Hlt2Global_Dec);
   fChain->SetBranchAddress("mu1_Hlt2Global_TIS", &mu1_Hlt2Global_TIS, &b_mu1_Hlt2Global_TIS);
   fChain->SetBranchAddress("mu1_Hlt2Global_TOS", &mu1_Hlt2Global_TOS, &b_mu1_Hlt2Global_TOS);
   fChain->SetBranchAddress("mu1_Hlt2Phys_Dec", &mu1_Hlt2Phys_Dec, &b_mu1_Hlt2Phys_Dec);
   fChain->SetBranchAddress("mu1_Hlt2Phys_TIS", &mu1_Hlt2Phys_TIS, &b_mu1_Hlt2Phys_TIS);
   fChain->SetBranchAddress("mu1_Hlt2Phys_TOS", &mu1_Hlt2Phys_TOS, &b_mu1_Hlt2Phys_TOS);
   fChain->SetBranchAddress("mu1_TRACK_Type", &mu1_TRACK_Type, &b_mu1_TRACK_Type);
   fChain->SetBranchAddress("mu1_TRACK_Key", &mu1_TRACK_Key, &b_mu1_TRACK_Key);
   fChain->SetBranchAddress("mu1_TRACK_CHI2", &mu1_TRACK_CHI2, &b_mu1_TRACK_CHI2);
   fChain->SetBranchAddress("mu1_TRACK_NDOF", &mu1_TRACK_NDOF, &b_mu1_TRACK_NDOF);
   fChain->SetBranchAddress("mu1_TRACK_CHI2NDOF", &mu1_TRACK_CHI2NDOF, &b_mu1_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("mu1_TRACK_PCHI2", &mu1_TRACK_PCHI2, &b_mu1_TRACK_PCHI2);
   fChain->SetBranchAddress("mu1_TRACK_VeloCHI2NDOF", &mu1_TRACK_VeloCHI2NDOF, &b_mu1_TRACK_VeloCHI2NDOF);
   fChain->SetBranchAddress("mu1_TRACK_TCHI2NDOF", &mu1_TRACK_TCHI2NDOF, &b_mu1_TRACK_TCHI2NDOF);
   fChain->SetBranchAddress("mu1_VELO_UTID", &mu1_VELO_UTID, &b_mu1_VELO_UTID);
   fChain->SetBranchAddress("mu1_TRACK_FirstMeasurementX", &mu1_TRACK_FirstMeasurementX, &b_mu1_TRACK_FirstMeasurementX);
   fChain->SetBranchAddress("mu1_TRACK_FirstMeasurementY", &mu1_TRACK_FirstMeasurementY, &b_mu1_TRACK_FirstMeasurementY);
   fChain->SetBranchAddress("mu1_TRACK_FirstMeasurementZ", &mu1_TRACK_FirstMeasurementZ, &b_mu1_TRACK_FirstMeasurementZ);
   fChain->SetBranchAddress("mu1_TRACK_MatchCHI2", &mu1_TRACK_MatchCHI2, &b_mu1_TRACK_MatchCHI2);
   fChain->SetBranchAddress("mu1_TRACK_GhostProb", &mu1_TRACK_GhostProb, &b_mu1_TRACK_GhostProb);
   fChain->SetBranchAddress("mu1_TRACK_CloneDist", &mu1_TRACK_CloneDist, &b_mu1_TRACK_CloneDist);
   fChain->SetBranchAddress("mu1_TRACK_Likelihood", &mu1_TRACK_Likelihood, &b_mu1_TRACK_Likelihood);
   fChain->SetBranchAddress("mu1_L0HadronDecision_Dec", &mu1_L0HadronDecision_Dec, &b_mu1_L0HadronDecision_Dec);
   fChain->SetBranchAddress("mu1_L0HadronDecision_TIS", &mu1_L0HadronDecision_TIS, &b_mu1_L0HadronDecision_TIS);
   fChain->SetBranchAddress("mu1_L0HadronDecision_TOS", &mu1_L0HadronDecision_TOS, &b_mu1_L0HadronDecision_TOS);
   fChain->SetBranchAddress("mu1_L0MuonDecision_Dec", &mu1_L0MuonDecision_Dec, &b_mu1_L0MuonDecision_Dec);
   fChain->SetBranchAddress("mu1_L0MuonDecision_TIS", &mu1_L0MuonDecision_TIS, &b_mu1_L0MuonDecision_TIS);
   fChain->SetBranchAddress("mu1_L0MuonDecision_TOS", &mu1_L0MuonDecision_TOS, &b_mu1_L0MuonDecision_TOS);
   fChain->SetBranchAddress("mu1_L0DiMuonDecision_Dec", &mu1_L0DiMuonDecision_Dec, &b_mu1_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("mu1_L0DiMuonDecision_TIS", &mu1_L0DiMuonDecision_TIS, &b_mu1_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("mu1_L0DiMuonDecision_TOS", &mu1_L0DiMuonDecision_TOS, &b_mu1_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("mu1_L0ElectronDecision_Dec", &mu1_L0ElectronDecision_Dec, &b_mu1_L0ElectronDecision_Dec);
   fChain->SetBranchAddress("mu1_L0ElectronDecision_TIS", &mu1_L0ElectronDecision_TIS, &b_mu1_L0ElectronDecision_TIS);
   fChain->SetBranchAddress("mu1_L0ElectronDecision_TOS", &mu1_L0ElectronDecision_TOS, &b_mu1_L0ElectronDecision_TOS);
   fChain->SetBranchAddress("mu1_L0PhotonDecision_Dec", &mu1_L0PhotonDecision_Dec, &b_mu1_L0PhotonDecision_Dec);
   fChain->SetBranchAddress("mu1_L0PhotonDecision_TIS", &mu1_L0PhotonDecision_TIS, &b_mu1_L0PhotonDecision_TIS);
   fChain->SetBranchAddress("mu1_L0PhotonDecision_TOS", &mu1_L0PhotonDecision_TOS, &b_mu1_L0PhotonDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt1DiMuonHighMassDecision_Dec", &mu1_Hlt1DiMuonHighMassDecision_Dec, &b_mu1_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt1DiMuonHighMassDecision_TIS", &mu1_Hlt1DiMuonHighMassDecision_TIS, &b_mu1_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt1DiMuonHighMassDecision_TOS", &mu1_Hlt1DiMuonHighMassDecision_TOS, &b_mu1_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt1DiMuonLowMassDecision_Dec", &mu1_Hlt1DiMuonLowMassDecision_Dec, &b_mu1_Hlt1DiMuonLowMassDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt1DiMuonLowMassDecision_TIS", &mu1_Hlt1DiMuonLowMassDecision_TIS, &b_mu1_Hlt1DiMuonLowMassDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt1DiMuonLowMassDecision_TOS", &mu1_Hlt1DiMuonLowMassDecision_TOS, &b_mu1_Hlt1DiMuonLowMassDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt1SingleMuonNoIPDecision_Dec", &mu1_Hlt1SingleMuonNoIPDecision_Dec, &b_mu1_Hlt1SingleMuonNoIPDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt1SingleMuonNoIPDecision_TIS", &mu1_Hlt1SingleMuonNoIPDecision_TIS, &b_mu1_Hlt1SingleMuonNoIPDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt1SingleMuonNoIPDecision_TOS", &mu1_Hlt1SingleMuonNoIPDecision_TOS, &b_mu1_Hlt1SingleMuonNoIPDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt1SingleMuonHighPTDecision_Dec", &mu1_Hlt1SingleMuonHighPTDecision_Dec, &b_mu1_Hlt1SingleMuonHighPTDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt1SingleMuonHighPTDecision_TIS", &mu1_Hlt1SingleMuonHighPTDecision_TIS, &b_mu1_Hlt1SingleMuonHighPTDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt1SingleMuonHighPTDecision_TOS", &mu1_Hlt1SingleMuonHighPTDecision_TOS, &b_mu1_Hlt1SingleMuonHighPTDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt1TrackAllL0Decision_Dec", &mu1_Hlt1TrackAllL0Decision_Dec, &b_mu1_Hlt1TrackAllL0Decision_Dec);
   fChain->SetBranchAddress("mu1_Hlt1TrackAllL0Decision_TIS", &mu1_Hlt1TrackAllL0Decision_TIS, &b_mu1_Hlt1TrackAllL0Decision_TIS);
   fChain->SetBranchAddress("mu1_Hlt1TrackAllL0Decision_TOS", &mu1_Hlt1TrackAllL0Decision_TOS, &b_mu1_Hlt1TrackAllL0Decision_TOS);
   fChain->SetBranchAddress("mu1_Hlt1TrackMuonDecision_Dec", &mu1_Hlt1TrackMuonDecision_Dec, &b_mu1_Hlt1TrackMuonDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TIS", &mu1_Hlt1TrackMuonDecision_TIS, &b_mu1_Hlt1TrackMuonDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS", &mu1_Hlt1TrackMuonDecision_TOS, &b_mu1_Hlt1TrackMuonDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt1TrackPhotonDecision_Dec", &mu1_Hlt1TrackPhotonDecision_Dec, &b_mu1_Hlt1TrackPhotonDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt1TrackPhotonDecision_TIS", &mu1_Hlt1TrackPhotonDecision_TIS, &b_mu1_Hlt1TrackPhotonDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt1TrackPhotonDecision_TOS", &mu1_Hlt1TrackPhotonDecision_TOS, &b_mu1_Hlt1TrackPhotonDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt1L0AnyDecision_Dec", &mu1_Hlt1L0AnyDecision_Dec, &b_mu1_Hlt1L0AnyDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt1L0AnyDecision_TIS", &mu1_Hlt1L0AnyDecision_TIS, &b_mu1_Hlt1L0AnyDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt1L0AnyDecision_TOS", &mu1_Hlt1L0AnyDecision_TOS, &b_mu1_Hlt1L0AnyDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt1GlobalDecision_Dec", &mu1_Hlt1GlobalDecision_Dec, &b_mu1_Hlt1GlobalDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt1GlobalDecision_TIS", &mu1_Hlt1GlobalDecision_TIS, &b_mu1_Hlt1GlobalDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt1GlobalDecision_TOS", &mu1_Hlt1GlobalDecision_TOS, &b_mu1_Hlt1GlobalDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2SingleMuonDecision_Dec", &mu1_Hlt2SingleMuonDecision_Dec, &b_mu1_Hlt2SingleMuonDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2SingleMuonDecision_TIS", &mu1_Hlt2SingleMuonDecision_TIS, &b_mu1_Hlt2SingleMuonDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2SingleMuonDecision_TOS", &mu1_Hlt2SingleMuonDecision_TOS, &b_mu1_Hlt2SingleMuonDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2DiMuonDetachedDecision_Dec", &mu1_Hlt2DiMuonDetachedDecision_Dec, &b_mu1_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2DiMuonDetachedDecision_TIS", &mu1_Hlt2DiMuonDetachedDecision_TIS, &b_mu1_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2DiMuonDetachedDecision_TOS", &mu1_Hlt2DiMuonDetachedDecision_TOS, &b_mu1_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD2HMuMuDecision_Dec", &mu1_Hlt2CharmSemilepD2HMuMuDecision_Dec, &b_mu1_Hlt2CharmSemilepD2HMuMuDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD2HMuMuDecision_TIS", &mu1_Hlt2CharmSemilepD2HMuMuDecision_TIS, &b_mu1_Hlt2CharmSemilepD2HMuMuDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD2HMuMuDecision_TOS", &mu1_Hlt2CharmSemilepD2HMuMuDecision_TOS, &b_mu1_Hlt2CharmSemilepD2HMuMuDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec", &mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec, &b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS", &mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS, &b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS", &mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS, &b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec", &mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec, &b_mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS", &mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS, &b_mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS", &mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS, &b_mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec", &mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec, &b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS", &mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS, &b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS", &mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS, &b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec", &mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec, &b_mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS", &mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS, &b_mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS", &mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS, &b_mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec", &mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec, &b_mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS", &mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS, &b_mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS", &mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS, &b_mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD02KKMuMuDecision_Dec", &mu1_Hlt2CharmSemilepD02KKMuMuDecision_Dec, &b_mu1_Hlt2CharmSemilepD02KKMuMuDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD02KKMuMuDecision_TIS", &mu1_Hlt2CharmSemilepD02KKMuMuDecision_TIS, &b_mu1_Hlt2CharmSemilepD02KKMuMuDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD02KKMuMuDecision_TOS", &mu1_Hlt2CharmSemilepD02KKMuMuDecision_TOS, &b_mu1_Hlt2CharmSemilepD02KKMuMuDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD02KPiMuMuDecision_Dec", &mu1_Hlt2CharmSemilepD02KPiMuMuDecision_Dec, &b_mu1_Hlt2CharmSemilepD02KPiMuMuDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD02KPiMuMuDecision_TIS", &mu1_Hlt2CharmSemilepD02KPiMuMuDecision_TIS, &b_mu1_Hlt2CharmSemilepD02KPiMuMuDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmSemilepD02KPiMuMuDecision_TOS", &mu1_Hlt2CharmSemilepD02KPiMuMuDecision_TOS, &b_mu1_Hlt2CharmSemilepD02KPiMuMuDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHHDst_4piDecision_Dec", &mu1_Hlt2CharmHadD02HHHHDst_4piDecision_Dec, &b_mu1_Hlt2CharmHadD02HHHHDst_4piDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHHDst_4piDecision_TIS", &mu1_Hlt2CharmHadD02HHHHDst_4piDecision_TIS, &b_mu1_Hlt2CharmHadD02HHHHDst_4piDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHHDst_4piDecision_TOS", &mu1_Hlt2CharmHadD02HHHHDst_4piDecision_TOS, &b_mu1_Hlt2CharmHadD02HHHHDst_4piDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec", &mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec, &b_mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS", &mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS, &b_mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS", &mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS, &b_mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec", &mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec, &b_mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS", &mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS, &b_mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS", &mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS, &b_mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_hhXDecision_Dec", &mu1_Hlt2CharmHadD02HHXDst_hhXDecision_Dec, &b_mu1_Hlt2CharmHadD02HHXDst_hhXDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_hhXDecision_TIS", &mu1_Hlt2CharmHadD02HHXDst_hhXDecision_TIS, &b_mu1_Hlt2CharmHadD02HHXDst_hhXDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_hhXDecision_TOS", &mu1_Hlt2CharmHadD02HHXDst_hhXDecision_TOS, &b_mu1_Hlt2CharmHadD02HHXDst_hhXDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec", &mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec, &b_mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS", &mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS, &b_mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS", &mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS, &b_mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec", &mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec, &b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS", &mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS, &b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS", &mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS, &b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec", &mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec, &b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS", &mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS, &b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS", &mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS, &b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec", &mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec, &b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS", &mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS, &b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS", &mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS, &b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec", &mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec, &b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS", &mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS, &b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS", &mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS, &b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_K3piDecision_Dec", &mu1_Hlt2CharmHadD02HHHH_K3piDecision_Dec, &b_mu1_Hlt2CharmHadD02HHHH_K3piDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_K3piDecision_TIS", &mu1_Hlt2CharmHadD02HHHH_K3piDecision_TIS, &b_mu1_Hlt2CharmHadD02HHHH_K3piDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_K3piDecision_TOS", &mu1_Hlt2CharmHadD02HHHH_K3piDecision_TOS, &b_mu1_Hlt2CharmHadD02HHHH_K3piDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec", &mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec, &b_mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS", &mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS, &b_mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS", &mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS, &b_mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec", &mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec, &b_mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS", &mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS, &b_mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS", &mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS, &b_mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec", &mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec, &b_mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS", &mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS, &b_mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS", &mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS, &b_mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_4piDecision_Dec", &mu1_Hlt2CharmHadD02HHHH_4piDecision_Dec, &b_mu1_Hlt2CharmHadD02HHHH_4piDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_4piDecision_TIS", &mu1_Hlt2CharmHadD02HHHH_4piDecision_TIS, &b_mu1_Hlt2CharmHadD02HHHH_4piDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_4piDecision_TOS", &mu1_Hlt2CharmHadD02HHHH_4piDecision_TOS, &b_mu1_Hlt2CharmHadD02HHHH_4piDecision_TOS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec", &mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec, &b_mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS", &mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS, &b_mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS);
   fChain->SetBranchAddress("mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS", &mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS, &b_mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS);
   fChain->SetBranchAddress("Slowpi_MINIP", &Slowpi_MINIP, &b_Slowpi_MINIP);
   fChain->SetBranchAddress("Slowpi_MINIPCHI2", &Slowpi_MINIPCHI2, &b_Slowpi_MINIPCHI2);
   fChain->SetBranchAddress("Slowpi_MINIPNEXTBEST", &Slowpi_MINIPNEXTBEST, &b_Slowpi_MINIPNEXTBEST);
   fChain->SetBranchAddress("Slowpi_MINIPCHI2NEXTBEST", &Slowpi_MINIPCHI2NEXTBEST, &b_Slowpi_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("Slowpi_OWNPV_X", &Slowpi_OWNPV_X, &b_Slowpi_OWNPV_X);
   fChain->SetBranchAddress("Slowpi_OWNPV_Y", &Slowpi_OWNPV_Y, &b_Slowpi_OWNPV_Y);
   fChain->SetBranchAddress("Slowpi_OWNPV_Z", &Slowpi_OWNPV_Z, &b_Slowpi_OWNPV_Z);
   fChain->SetBranchAddress("Slowpi_OWNPV_XERR", &Slowpi_OWNPV_XERR, &b_Slowpi_OWNPV_XERR);
   fChain->SetBranchAddress("Slowpi_OWNPV_YERR", &Slowpi_OWNPV_YERR, &b_Slowpi_OWNPV_YERR);
   fChain->SetBranchAddress("Slowpi_OWNPV_ZERR", &Slowpi_OWNPV_ZERR, &b_Slowpi_OWNPV_ZERR);
   fChain->SetBranchAddress("Slowpi_OWNPV_CHI2", &Slowpi_OWNPV_CHI2, &b_Slowpi_OWNPV_CHI2);
   fChain->SetBranchAddress("Slowpi_OWNPV_NDOF", &Slowpi_OWNPV_NDOF, &b_Slowpi_OWNPV_NDOF);
   fChain->SetBranchAddress("Slowpi_OWNPV_COV_", Slowpi_OWNPV_COV_, &b_Slowpi_OWNPV_COV_);
   fChain->SetBranchAddress("Slowpi_IP_OWNPV", &Slowpi_IP_OWNPV, &b_Slowpi_IP_OWNPV);
   fChain->SetBranchAddress("Slowpi_IPCHI2_OWNPV", &Slowpi_IPCHI2_OWNPV, &b_Slowpi_IPCHI2_OWNPV);
   fChain->SetBranchAddress("Slowpi_TOPPV_X", &Slowpi_TOPPV_X, &b_Slowpi_TOPPV_X);
   fChain->SetBranchAddress("Slowpi_TOPPV_Y", &Slowpi_TOPPV_Y, &b_Slowpi_TOPPV_Y);
   fChain->SetBranchAddress("Slowpi_TOPPV_Z", &Slowpi_TOPPV_Z, &b_Slowpi_TOPPV_Z);
   fChain->SetBranchAddress("Slowpi_TOPPV_XERR", &Slowpi_TOPPV_XERR, &b_Slowpi_TOPPV_XERR);
   fChain->SetBranchAddress("Slowpi_TOPPV_YERR", &Slowpi_TOPPV_YERR, &b_Slowpi_TOPPV_YERR);
   fChain->SetBranchAddress("Slowpi_TOPPV_ZERR", &Slowpi_TOPPV_ZERR, &b_Slowpi_TOPPV_ZERR);
   fChain->SetBranchAddress("Slowpi_TOPPV_CHI2", &Slowpi_TOPPV_CHI2, &b_Slowpi_TOPPV_CHI2);
   fChain->SetBranchAddress("Slowpi_TOPPV_NDOF", &Slowpi_TOPPV_NDOF, &b_Slowpi_TOPPV_NDOF);
   fChain->SetBranchAddress("Slowpi_TOPPV_COV_", Slowpi_TOPPV_COV_, &b_Slowpi_TOPPV_COV_);
   fChain->SetBranchAddress("Slowpi_IP_TOPPV", &Slowpi_IP_TOPPV, &b_Slowpi_IP_TOPPV);
   fChain->SetBranchAddress("Slowpi_IPCHI2_TOPPV", &Slowpi_IPCHI2_TOPPV, &b_Slowpi_IPCHI2_TOPPV);
   fChain->SetBranchAddress("Slowpi_ORIVX_X", &Slowpi_ORIVX_X, &b_Slowpi_ORIVX_X);
   fChain->SetBranchAddress("Slowpi_ORIVX_Y", &Slowpi_ORIVX_Y, &b_Slowpi_ORIVX_Y);
   fChain->SetBranchAddress("Slowpi_ORIVX_Z", &Slowpi_ORIVX_Z, &b_Slowpi_ORIVX_Z);
   fChain->SetBranchAddress("Slowpi_ORIVX_XERR", &Slowpi_ORIVX_XERR, &b_Slowpi_ORIVX_XERR);
   fChain->SetBranchAddress("Slowpi_ORIVX_YERR", &Slowpi_ORIVX_YERR, &b_Slowpi_ORIVX_YERR);
   fChain->SetBranchAddress("Slowpi_ORIVX_ZERR", &Slowpi_ORIVX_ZERR, &b_Slowpi_ORIVX_ZERR);
   fChain->SetBranchAddress("Slowpi_ORIVX_CHI2", &Slowpi_ORIVX_CHI2, &b_Slowpi_ORIVX_CHI2);
   fChain->SetBranchAddress("Slowpi_ORIVX_NDOF", &Slowpi_ORIVX_NDOF, &b_Slowpi_ORIVX_NDOF);
   fChain->SetBranchAddress("Slowpi_ORIVX_COV_", Slowpi_ORIVX_COV_, &b_Slowpi_ORIVX_COV_);
   fChain->SetBranchAddress("Slowpi_IP_ORIVX", &Slowpi_IP_ORIVX, &b_Slowpi_IP_ORIVX);
   fChain->SetBranchAddress("Slowpi_IPCHI2_ORIVX", &Slowpi_IPCHI2_ORIVX, &b_Slowpi_IPCHI2_ORIVX);
   fChain->SetBranchAddress("Slowpi_P", &Slowpi_P, &b_Slowpi_P);
   fChain->SetBranchAddress("Slowpi_PT", &Slowpi_PT, &b_Slowpi_PT);
   fChain->SetBranchAddress("Slowpi_PE", &Slowpi_PE, &b_Slowpi_PE);
   fChain->SetBranchAddress("Slowpi_PX", &Slowpi_PX, &b_Slowpi_PX);
   fChain->SetBranchAddress("Slowpi_PY", &Slowpi_PY, &b_Slowpi_PY);
   fChain->SetBranchAddress("Slowpi_PZ", &Slowpi_PZ, &b_Slowpi_PZ);
   fChain->SetBranchAddress("Slowpi_M", &Slowpi_M, &b_Slowpi_M);
   fChain->SetBranchAddress("Slowpi_L0Calo_HCAL_realET", &Slowpi_L0Calo_HCAL_realET, &b_Slowpi_L0Calo_HCAL_realET);
   fChain->SetBranchAddress("Slowpi_L0Calo_HCAL_xProjection", &Slowpi_L0Calo_HCAL_xProjection, &b_Slowpi_L0Calo_HCAL_xProjection);
   fChain->SetBranchAddress("Slowpi_L0Calo_HCAL_yProjection", &Slowpi_L0Calo_HCAL_yProjection, &b_Slowpi_L0Calo_HCAL_yProjection);
   fChain->SetBranchAddress("Slowpi_L0Calo_HCAL_region", &Slowpi_L0Calo_HCAL_region, &b_Slowpi_L0Calo_HCAL_region);
   fChain->SetBranchAddress("Slowpi_L0Calo_HCAL_TriggerET", &Slowpi_L0Calo_HCAL_TriggerET, &b_Slowpi_L0Calo_HCAL_TriggerET);
   fChain->SetBranchAddress("Slowpi_L0Calo_HCAL_TriggerHCALET", &Slowpi_L0Calo_HCAL_TriggerHCALET, &b_Slowpi_L0Calo_HCAL_TriggerHCALET);
   fChain->SetBranchAddress("Slowpi_L0Calo_HCAL_xTrigger", &Slowpi_L0Calo_HCAL_xTrigger, &b_Slowpi_L0Calo_HCAL_xTrigger);
   fChain->SetBranchAddress("Slowpi_L0Calo_HCAL_yTrigger", &Slowpi_L0Calo_HCAL_yTrigger, &b_Slowpi_L0Calo_HCAL_yTrigger);
   fChain->SetBranchAddress("Slowpi_ID", &Slowpi_ID, &b_Slowpi_ID);
   fChain->SetBranchAddress("Slowpi_CombDLLMu", &Slowpi_CombDLLMu, &b_Slowpi_CombDLLMu);
   fChain->SetBranchAddress("Slowpi_ProbNNmu", &Slowpi_ProbNNmu, &b_Slowpi_ProbNNmu);
   fChain->SetBranchAddress("Slowpi_ProbNNghost", &Slowpi_ProbNNghost, &b_Slowpi_ProbNNghost);
   fChain->SetBranchAddress("Slowpi_InMuonAcc", &Slowpi_InMuonAcc, &b_Slowpi_InMuonAcc);
   fChain->SetBranchAddress("Slowpi_MuonDist2", &Slowpi_MuonDist2, &b_Slowpi_MuonDist2);
   fChain->SetBranchAddress("Slowpi_regionInM2", &Slowpi_regionInM2, &b_Slowpi_regionInM2);
   fChain->SetBranchAddress("Slowpi_hasMuon", &Slowpi_hasMuon, &b_Slowpi_hasMuon);
   fChain->SetBranchAddress("Slowpi_isMuon", &Slowpi_isMuon, &b_Slowpi_isMuon);
   fChain->SetBranchAddress("Slowpi_isMuonLoose", &Slowpi_isMuonLoose, &b_Slowpi_isMuonLoose);
   fChain->SetBranchAddress("Slowpi_NShared", &Slowpi_NShared, &b_Slowpi_NShared);
   fChain->SetBranchAddress("Slowpi_MuonLLmu", &Slowpi_MuonLLmu, &b_Slowpi_MuonLLmu);
   fChain->SetBranchAddress("Slowpi_MuonLLbg", &Slowpi_MuonLLbg, &b_Slowpi_MuonLLbg);
   fChain->SetBranchAddress("Slowpi_isMuonFromProto", &Slowpi_isMuonFromProto, &b_Slowpi_isMuonFromProto);
   fChain->SetBranchAddress("Slowpi_PIDe", &Slowpi_PIDe, &b_Slowpi_PIDe);
   fChain->SetBranchAddress("Slowpi_PIDmu", &Slowpi_PIDmu, &b_Slowpi_PIDmu);
   fChain->SetBranchAddress("Slowpi_PIDK", &Slowpi_PIDK, &b_Slowpi_PIDK);
   fChain->SetBranchAddress("Slowpi_PIDp", &Slowpi_PIDp, &b_Slowpi_PIDp);
   fChain->SetBranchAddress("Slowpi_ProbNNe", &Slowpi_ProbNNe, &b_Slowpi_ProbNNe);
   fChain->SetBranchAddress("Slowpi_ProbNNk", &Slowpi_ProbNNk, &b_Slowpi_ProbNNk);
   fChain->SetBranchAddress("Slowpi_ProbNNp", &Slowpi_ProbNNp, &b_Slowpi_ProbNNp);
   fChain->SetBranchAddress("Slowpi_ProbNNpi", &Slowpi_ProbNNpi, &b_Slowpi_ProbNNpi);
   fChain->SetBranchAddress("Slowpi_hasRich", &Slowpi_hasRich, &b_Slowpi_hasRich);
   fChain->SetBranchAddress("Slowpi_hasCalo", &Slowpi_hasCalo, &b_Slowpi_hasCalo);
   fChain->SetBranchAddress("Slowpi_UsedRichAerogel", &Slowpi_UsedRichAerogel, &b_Slowpi_UsedRichAerogel);
   fChain->SetBranchAddress("Slowpi_UsedRich1Gas", &Slowpi_UsedRich1Gas, &b_Slowpi_UsedRich1Gas);
   fChain->SetBranchAddress("Slowpi_UsedRich2Gas", &Slowpi_UsedRich2Gas, &b_Slowpi_UsedRich2Gas);
   fChain->SetBranchAddress("Slowpi_RichAboveElThres", &Slowpi_RichAboveElThres, &b_Slowpi_RichAboveElThres);
   fChain->SetBranchAddress("Slowpi_RichAboveMuThres", &Slowpi_RichAboveMuThres, &b_Slowpi_RichAboveMuThres);
   fChain->SetBranchAddress("Slowpi_RichAbovePiThres", &Slowpi_RichAbovePiThres, &b_Slowpi_RichAbovePiThres);
   fChain->SetBranchAddress("Slowpi_RichAboveKaThres", &Slowpi_RichAboveKaThres, &b_Slowpi_RichAboveKaThres);
   fChain->SetBranchAddress("Slowpi_RichAbovePrThres", &Slowpi_RichAbovePrThres, &b_Slowpi_RichAbovePrThres);
   fChain->SetBranchAddress("Slowpi_RichDLLe", &Slowpi_RichDLLe, &b_Slowpi_RichDLLe);
   fChain->SetBranchAddress("Slowpi_RichDLLmu", &Slowpi_RichDLLmu, &b_Slowpi_RichDLLmu);
   fChain->SetBranchAddress("Slowpi_RichDLLpi", &Slowpi_RichDLLpi, &b_Slowpi_RichDLLpi);
   fChain->SetBranchAddress("Slowpi_RichDLLk", &Slowpi_RichDLLk, &b_Slowpi_RichDLLk);
   fChain->SetBranchAddress("Slowpi_RichDLLp", &Slowpi_RichDLLp, &b_Slowpi_RichDLLp);
   fChain->SetBranchAddress("Slowpi_RichDLLbt", &Slowpi_RichDLLbt, &b_Slowpi_RichDLLbt);
   fChain->SetBranchAddress("Slowpi_InAccMuon", &Slowpi_InAccMuon, &b_Slowpi_InAccMuon);
   fChain->SetBranchAddress("Slowpi_MuonMuLL", &Slowpi_MuonMuLL, &b_Slowpi_MuonMuLL);
   fChain->SetBranchAddress("Slowpi_MuonBkgLL", &Slowpi_MuonBkgLL, &b_Slowpi_MuonBkgLL);
   fChain->SetBranchAddress("Slowpi_MuonNShared", &Slowpi_MuonNShared, &b_Slowpi_MuonNShared);
   fChain->SetBranchAddress("Slowpi_InAccEcal", &Slowpi_InAccEcal, &b_Slowpi_InAccEcal);
   fChain->SetBranchAddress("Slowpi_CaloEcalE", &Slowpi_CaloEcalE, &b_Slowpi_CaloEcalE);
   fChain->SetBranchAddress("Slowpi_EcalPIDe", &Slowpi_EcalPIDe, &b_Slowpi_EcalPIDe);
   fChain->SetBranchAddress("Slowpi_EcalPIDmu", &Slowpi_EcalPIDmu, &b_Slowpi_EcalPIDmu);
   fChain->SetBranchAddress("Slowpi_InAccHcal", &Slowpi_InAccHcal, &b_Slowpi_InAccHcal);
   fChain->SetBranchAddress("Slowpi_CaloHcalE", &Slowpi_CaloHcalE, &b_Slowpi_CaloHcalE);
   fChain->SetBranchAddress("Slowpi_HcalPIDe", &Slowpi_HcalPIDe, &b_Slowpi_HcalPIDe);
   fChain->SetBranchAddress("Slowpi_HcalPIDmu", &Slowpi_HcalPIDmu, &b_Slowpi_HcalPIDmu);
   fChain->SetBranchAddress("Slowpi_InAccPrs", &Slowpi_InAccPrs, &b_Slowpi_InAccPrs);
   fChain->SetBranchAddress("Slowpi_PrsPIDe", &Slowpi_PrsPIDe, &b_Slowpi_PrsPIDe);
   fChain->SetBranchAddress("Slowpi_CaloPrsE", &Slowpi_CaloPrsE, &b_Slowpi_CaloPrsE);
   fChain->SetBranchAddress("Slowpi_InAccSpd", &Slowpi_InAccSpd, &b_Slowpi_InAccSpd);
   fChain->SetBranchAddress("Slowpi_CaloSpdE", &Slowpi_CaloSpdE, &b_Slowpi_CaloSpdE);
   fChain->SetBranchAddress("Slowpi_InAccBrem", &Slowpi_InAccBrem, &b_Slowpi_InAccBrem);
   fChain->SetBranchAddress("Slowpi_BremPIDe", &Slowpi_BremPIDe, &b_Slowpi_BremPIDe);
   fChain->SetBranchAddress("Slowpi_VeloCharge", &Slowpi_VeloCharge, &b_Slowpi_VeloCharge);
   fChain->SetBranchAddress("Slowpi_RICHDLLe", &Slowpi_RICHDLLe, &b_Slowpi_RICHDLLe);
   fChain->SetBranchAddress("Slowpi_RICHDLLmu", &Slowpi_RICHDLLmu, &b_Slowpi_RICHDLLmu);
   fChain->SetBranchAddress("Slowpi_RICHDLLpi", &Slowpi_RICHDLLpi, &b_Slowpi_RICHDLLpi);
   fChain->SetBranchAddress("Slowpi_RICHDLLK", &Slowpi_RICHDLLK, &b_Slowpi_RICHDLLK);
   fChain->SetBranchAddress("Slowpi_RICHDLLp", &Slowpi_RICHDLLp, &b_Slowpi_RICHDLLp);
   fChain->SetBranchAddress("Slowpi_RICHDLLbt", &Slowpi_RICHDLLbt, &b_Slowpi_RICHDLLbt);
   fChain->SetBranchAddress("Slowpi_RICHBestID", &Slowpi_RICHBestID, &b_Slowpi_RICHBestID);
   fChain->SetBranchAddress("Slowpi_RICHThreshold", &Slowpi_RICHThreshold, &b_Slowpi_RICHThreshold);
   fChain->SetBranchAddress("Slowpi_RICHThresholdEl", &Slowpi_RICHThresholdEl, &b_Slowpi_RICHThresholdEl);
   fChain->SetBranchAddress("Slowpi_RICHThresholdMu", &Slowpi_RICHThresholdMu, &b_Slowpi_RICHThresholdMu);
   fChain->SetBranchAddress("Slowpi_RICHThresholdPi", &Slowpi_RICHThresholdPi, &b_Slowpi_RICHThresholdPi);
   fChain->SetBranchAddress("Slowpi_RICHThresholdKa", &Slowpi_RICHThresholdKa, &b_Slowpi_RICHThresholdKa);
   fChain->SetBranchAddress("Slowpi_RICHThresholdPr", &Slowpi_RICHThresholdPr, &b_Slowpi_RICHThresholdPr);
   fChain->SetBranchAddress("Slowpi_RICHAerogelUsed", &Slowpi_RICHAerogelUsed, &b_Slowpi_RICHAerogelUsed);
   fChain->SetBranchAddress("Slowpi_RICH1GasUsed", &Slowpi_RICH1GasUsed, &b_Slowpi_RICH1GasUsed);
   fChain->SetBranchAddress("Slowpi_RICH2GasUsed", &Slowpi_RICH2GasUsed, &b_Slowpi_RICH2GasUsed);
   fChain->SetBranchAddress("Slowpi_TRACK_Eta", &Slowpi_TRACK_Eta, &b_Slowpi_TRACK_Eta);
   fChain->SetBranchAddress("Slowpi_TRACK_Phi", &Slowpi_TRACK_Phi, &b_Slowpi_TRACK_Phi);
   fChain->SetBranchAddress("Slowpi_Aerogel_X", &Slowpi_Aerogel_X, &b_Slowpi_Aerogel_X);
   fChain->SetBranchAddress("Slowpi_Aerogel_Y", &Slowpi_Aerogel_Y, &b_Slowpi_Aerogel_Y);
   fChain->SetBranchAddress("Slowpi_Aerogel_Z", &Slowpi_Aerogel_Z, &b_Slowpi_Aerogel_Z);
   fChain->SetBranchAddress("Slowpi_Aerogel_Rho", &Slowpi_Aerogel_Rho, &b_Slowpi_Aerogel_Rho);
   fChain->SetBranchAddress("Slowpi_Aerogel_Phi", &Slowpi_Aerogel_Phi, &b_Slowpi_Aerogel_Phi);
   fChain->SetBranchAddress("Slowpi_Rich1Gas_X", &Slowpi_Rich1Gas_X, &b_Slowpi_Rich1Gas_X);
   fChain->SetBranchAddress("Slowpi_Rich1Gas_Y", &Slowpi_Rich1Gas_Y, &b_Slowpi_Rich1Gas_Y);
   fChain->SetBranchAddress("Slowpi_Rich1Gas_Z", &Slowpi_Rich1Gas_Z, &b_Slowpi_Rich1Gas_Z);
   fChain->SetBranchAddress("Slowpi_Rich1Gas_Rho", &Slowpi_Rich1Gas_Rho, &b_Slowpi_Rich1Gas_Rho);
   fChain->SetBranchAddress("Slowpi_Rich1Gas_Phi", &Slowpi_Rich1Gas_Phi, &b_Slowpi_Rich1Gas_Phi);
   fChain->SetBranchAddress("Slowpi_Rich2Gas_X", &Slowpi_Rich2Gas_X, &b_Slowpi_Rich2Gas_X);
   fChain->SetBranchAddress("Slowpi_Rich2Gas_Y", &Slowpi_Rich2Gas_Y, &b_Slowpi_Rich2Gas_Y);
   fChain->SetBranchAddress("Slowpi_Rich2Gas_Z", &Slowpi_Rich2Gas_Z, &b_Slowpi_Rich2Gas_Z);
   fChain->SetBranchAddress("Slowpi_Rich2Gas_Rho", &Slowpi_Rich2Gas_Rho, &b_Slowpi_Rich2Gas_Rho);
   fChain->SetBranchAddress("Slowpi_Rich2Gas_Phi", &Slowpi_Rich2Gas_Phi, &b_Slowpi_Rich2Gas_Phi);
   fChain->SetBranchAddress("Slowpi_L0Global_Dec", &Slowpi_L0Global_Dec, &b_Slowpi_L0Global_Dec);
   fChain->SetBranchAddress("Slowpi_L0Global_TIS", &Slowpi_L0Global_TIS, &b_Slowpi_L0Global_TIS);
   fChain->SetBranchAddress("Slowpi_L0Global_TOS", &Slowpi_L0Global_TOS, &b_Slowpi_L0Global_TOS);
   fChain->SetBranchAddress("Slowpi_Hlt1Global_Dec", &Slowpi_Hlt1Global_Dec, &b_Slowpi_Hlt1Global_Dec);
   fChain->SetBranchAddress("Slowpi_Hlt1Global_TIS", &Slowpi_Hlt1Global_TIS, &b_Slowpi_Hlt1Global_TIS);
   fChain->SetBranchAddress("Slowpi_Hlt1Global_TOS", &Slowpi_Hlt1Global_TOS, &b_Slowpi_Hlt1Global_TOS);
   fChain->SetBranchAddress("Slowpi_Hlt1Phys_Dec", &Slowpi_Hlt1Phys_Dec, &b_Slowpi_Hlt1Phys_Dec);
   fChain->SetBranchAddress("Slowpi_Hlt1Phys_TIS", &Slowpi_Hlt1Phys_TIS, &b_Slowpi_Hlt1Phys_TIS);
   fChain->SetBranchAddress("Slowpi_Hlt1Phys_TOS", &Slowpi_Hlt1Phys_TOS, &b_Slowpi_Hlt1Phys_TOS);
   fChain->SetBranchAddress("Slowpi_Hlt2Global_Dec", &Slowpi_Hlt2Global_Dec, &b_Slowpi_Hlt2Global_Dec);
   fChain->SetBranchAddress("Slowpi_Hlt2Global_TIS", &Slowpi_Hlt2Global_TIS, &b_Slowpi_Hlt2Global_TIS);
   fChain->SetBranchAddress("Slowpi_Hlt2Global_TOS", &Slowpi_Hlt2Global_TOS, &b_Slowpi_Hlt2Global_TOS);
   fChain->SetBranchAddress("Slowpi_Hlt2Phys_Dec", &Slowpi_Hlt2Phys_Dec, &b_Slowpi_Hlt2Phys_Dec);
   fChain->SetBranchAddress("Slowpi_Hlt2Phys_TIS", &Slowpi_Hlt2Phys_TIS, &b_Slowpi_Hlt2Phys_TIS);
   fChain->SetBranchAddress("Slowpi_Hlt2Phys_TOS", &Slowpi_Hlt2Phys_TOS, &b_Slowpi_Hlt2Phys_TOS);
   fChain->SetBranchAddress("Slowpi_TRACK_Type", &Slowpi_TRACK_Type, &b_Slowpi_TRACK_Type);
   fChain->SetBranchAddress("Slowpi_TRACK_Key", &Slowpi_TRACK_Key, &b_Slowpi_TRACK_Key);
   fChain->SetBranchAddress("Slowpi_TRACK_CHI2", &Slowpi_TRACK_CHI2, &b_Slowpi_TRACK_CHI2);
   fChain->SetBranchAddress("Slowpi_TRACK_NDOF", &Slowpi_TRACK_NDOF, &b_Slowpi_TRACK_NDOF);
   fChain->SetBranchAddress("Slowpi_TRACK_CHI2NDOF", &Slowpi_TRACK_CHI2NDOF, &b_Slowpi_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("Slowpi_TRACK_PCHI2", &Slowpi_TRACK_PCHI2, &b_Slowpi_TRACK_PCHI2);
   fChain->SetBranchAddress("Slowpi_TRACK_VeloCHI2NDOF", &Slowpi_TRACK_VeloCHI2NDOF, &b_Slowpi_TRACK_VeloCHI2NDOF);
   fChain->SetBranchAddress("Slowpi_TRACK_TCHI2NDOF", &Slowpi_TRACK_TCHI2NDOF, &b_Slowpi_TRACK_TCHI2NDOF);
   fChain->SetBranchAddress("Slowpi_VELO_UTID", &Slowpi_VELO_UTID, &b_Slowpi_VELO_UTID);
   fChain->SetBranchAddress("Slowpi_TRACK_FirstMeasurementX", &Slowpi_TRACK_FirstMeasurementX, &b_Slowpi_TRACK_FirstMeasurementX);
   fChain->SetBranchAddress("Slowpi_TRACK_FirstMeasurementY", &Slowpi_TRACK_FirstMeasurementY, &b_Slowpi_TRACK_FirstMeasurementY);
   fChain->SetBranchAddress("Slowpi_TRACK_FirstMeasurementZ", &Slowpi_TRACK_FirstMeasurementZ, &b_Slowpi_TRACK_FirstMeasurementZ);
   fChain->SetBranchAddress("Slowpi_TRACK_MatchCHI2", &Slowpi_TRACK_MatchCHI2, &b_Slowpi_TRACK_MatchCHI2);
   fChain->SetBranchAddress("Slowpi_TRACK_GhostProb", &Slowpi_TRACK_GhostProb, &b_Slowpi_TRACK_GhostProb);
   fChain->SetBranchAddress("Slowpi_TRACK_CloneDist", &Slowpi_TRACK_CloneDist, &b_Slowpi_TRACK_CloneDist);
   fChain->SetBranchAddress("Slowpi_TRACK_Likelihood", &Slowpi_TRACK_Likelihood, &b_Slowpi_TRACK_Likelihood);
   fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
   fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
   fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
   fChain->SetBranchAddress("nBackward", &nBackward, &b_nBackward);
   fChain->SetBranchAddress("nDownstream", &nDownstream, &b_nDownstream);
   fChain->SetBranchAddress("nITClusters", &nITClusters, &b_nITClusters);
   fChain->SetBranchAddress("nLong", &nLong, &b_nLong);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("nPVs", &nPVs, &b_nPVs);
   fChain->SetBranchAddress("nSpdDigits", &nSpdDigits, &b_nSpdDigits);
   fChain->SetBranchAddress("nTTClusters", &nTTClusters, &b_nTTClusters);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("nUpstream", &nUpstream, &b_nUpstream);
   fChain->SetBranchAddress("nVELO", &nVELO, &b_nVELO);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
   fChain->SetBranchAddress("BCType", &BCType, &b_BCType);
   fChain->SetBranchAddress("OdinTCK", &OdinTCK, &b_OdinTCK);
   fChain->SetBranchAddress("L0DUTCK", &L0DUTCK, &b_L0DUTCK);
   fChain->SetBranchAddress("HLT1TCK", &HLT1TCK, &b_HLT1TCK);
   fChain->SetBranchAddress("HLT2TCK", &HLT2TCK, &b_HLT2TCK);
   fChain->SetBranchAddress("GpsTime", &GpsTime, &b_GpsTime);
   fChain->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVX", PVX, &b_PVX);
   fChain->SetBranchAddress("PVY", PVY, &b_PVY);
   fChain->SetBranchAddress("PVZ", PVZ, &b_PVZ);
   fChain->SetBranchAddress("PVXERR", PVXERR, &b_PVXERR);
   fChain->SetBranchAddress("PVYERR", PVYERR, &b_PVYERR);
   fChain->SetBranchAddress("PVZERR", PVZERR, &b_PVZERR);
   fChain->SetBranchAddress("PVCHI2", PVCHI2, &b_PVCHI2);
   fChain->SetBranchAddress("PVNDOF", PVNDOF, &b_PVNDOF);
   fChain->SetBranchAddress("PVNTRACKS", PVNTRACKS, &b_PVNTRACKS);
   Notify();
}

void D2hhmumuReader::activateRelevantBranches()
{
 
  fChain->SetBranchStatus("*",0);
  
 //kinematic

  //DTF

  fChain->SetBranchStatus("Dst_DTF_Pis_PX",1);
  fChain->SetBranchStatus("Dst_DTF_Pis_PY",1);
  fChain->SetBranchStatus("Dst_DTF_Pis_PZ",1);
  fChain->SetBranchStatus("Dst_DTF_h1_PX",1);
  fChain->SetBranchStatus("Dst_DTF_h1_PY",1);
  fChain->SetBranchStatus("Dst_DTF_h1_PZ",1);
  fChain->SetBranchStatus("Dst_DTF_h0_PX",1);
  fChain->SetBranchStatus("Dst_DTF_h0_PY",1);
  fChain->SetBranchStatus("Dst_DTF_h0_PZ",1);
  fChain->SetBranchStatus("Dst_DTF_mu0_PX",1); 
  fChain->SetBranchStatus("Dst_DTF_mu0_PY",1);
  fChain->SetBranchStatus("Dst_DTF_mu0_PZ",1);
  fChain->SetBranchStatus("Dst_DTF_mu1_PX",1); 
  fChain->SetBranchStatus("Dst_DTF_mu1_PY",1);
  fChain->SetBranchStatus("Dst_DTF_mu1_PZ",1);
  fChain->SetBranchStatus("Dst_DTF_mu1_P",1);
  fChain->SetBranchStatus("Dst_DTF_mu1_PT",1);
  fChain->SetBranchStatus("Dst_DTF_mu0_P",1);
  fChain->SetBranchStatus("Dst_DTF_mu0_PT",1);
  fChain->SetBranchStatus("Dst_DTF_h1_P",1);
  fChain->SetBranchStatus("Dst_DTF_h1_PT",1);
  fChain->SetBranchStatus("Dst_DTF_h0_P",1);
  fChain->SetBranchStatus("Dst_DTF_h0_PT",1);
  fChain->SetBranchStatus("Dst_DTF_D0_PX",1);
  fChain->SetBranchStatus("Dst_DTF_D0_PY",1);
  fChain->SetBranchStatus("Dst_DTF_D0_PZ",1);
  fChain->SetBranchStatus("Dst_DTF_D0_P",1);
  fChain->SetBranchStatus("Dst_DTF_D0_PT",1);
  fChain->SetBranchStatus("Dst_DTF_D0_M",1);
  fChain->SetBranchStatus("Dst_DTF_D0_BPVIPCHI2",1);
  fChain->SetBranchStatus("Dst_DTF_CHI2",1);
  fChain->SetBranchStatus("Dst_DTF_NDOF",1);
  fChain->SetBranchStatus("Dst_DTF_Dstarplus_M",1);
  fChain->SetBranchStatus("Dst_DTF_Dstarplus_P",1);
  fChain->SetBranchStatus("Dst_DTF_Dstarplus_PT",1);
  fChain->SetBranchStatus("Dst_DTF_Dstarplus_PX",1);
  fChain->SetBranchStatus("Dst_DTF_Dstarplus_PY",1);
  fChain->SetBranchStatus("Dst_DTF_Dstarplus_PZ",1);

  fChain->SetBranchStatus("Dst_CONEANGLE_D",1);
  fChain->SetBranchStatus("Dst_CONEANGLE_Dstar",1);
  fChain->SetBranchStatus("Dst_CONEMULT_D",1);
  fChain->SetBranchStatus("Dst_CONEMULT_Dstar",1);
  fChain->SetBranchStatus("Dst_CONEPTASYM_D",1);
  fChain->SetBranchStatus("Dst_CONEPTASYM_Dstar",1);

  fChain->SetBranchStatus("Dst_PT",1);
  fChain->SetBranchStatus("Dst_PX",1);
  fChain->SetBranchStatus("Dst_PY",1);
  fChain->SetBranchStatus("Dst_PZ",1);
  fChain->SetBranchStatus("Dst_P",1);
  fChain->SetBranchStatus("Dst_M",1);

  fChain->SetBranchStatus("D_PT",1);
  fChain->SetBranchStatus("D_PX",1);
  fChain->SetBranchStatus("D_PY",1);
  fChain->SetBranchStatus("D_PZ",1);
  fChain->SetBranchStatus("D_P",1);
  fChain->SetBranchStatus("D_DiMuon_Mass",1);
  fChain->SetBranchStatus("D_M",1);

  fChain->SetBranchStatus("D_TAU",1);
  fChain->SetBranchStatus("D_TAUCHI2",1);

  fChain->SetBranchStatus("mu0_PT",1);
  fChain->SetBranchStatus("mu0_PX",1);
  fChain->SetBranchStatus("mu0_PY",1);
  fChain->SetBranchStatus("mu0_PZ",1);
  fChain->SetBranchStatus("mu0_P",1);

  fChain->SetBranchStatus("mu1_PT",1);
  fChain->SetBranchStatus("mu1_PX",1);
  fChain->SetBranchStatus("mu1_PY",1);
  fChain->SetBranchStatus("mu1_PZ",1);
  fChain->SetBranchStatus("mu1_P",1);

  fChain->SetBranchStatus("h0_PT",1);
  fChain->SetBranchStatus("h0_PX",1);
  fChain->SetBranchStatus("h0_PY",1);
  fChain->SetBranchStatus("h0_PZ",1);
  fChain->SetBranchStatus("h0_P",1);

  fChain->SetBranchStatus("h1_PT",1);
  fChain->SetBranchStatus("h1_PX",1);
  fChain->SetBranchStatus("h1_PY",1);
  fChain->SetBranchStatus("h1_PZ",1);
  fChain->SetBranchStatus("h1_P",1);

  fChain->SetBranchStatus("Slowpi_PT",1);
  fChain->SetBranchStatus("Slowpi_PX",1);
  fChain->SetBranchStatus("Slowpi_PY",1);
  fChain->SetBranchStatus("Slowpi_PZ",1);
  fChain->SetBranchStatus("Slowpi_P",1);

  //PID

  fChain->SetBranchStatus("Slowpi_PIDe",1);
  fChain->SetBranchStatus("Slowpi_PIDmu",1);
  fChain->SetBranchStatus("Slowpi_PIDK",1);
  fChain->SetBranchStatus("Slowpi_PIDp",1);
  fChain->SetBranchStatus("Slowpi_ProbNNe",1);
  fChain->SetBranchStatus("Slowpi_ProbNNk",1);
  fChain->SetBranchStatus("Slowpi_ProbNNp",1);
  fChain->SetBranchStatus("Slowpi_ProbNNpi",1);
  fChain->SetBranchStatus("Slowpi_ProbNNghost",1);

  fChain->SetBranchStatus("mu0_PIDe",1);
  fChain->SetBranchStatus("mu0_PIDmu",1);
  fChain->SetBranchStatus("mu0_PIDK",1);
  fChain->SetBranchStatus("mu0_PIDp",1);
  fChain->SetBranchStatus("mu0_ProbNNe",1);
  fChain->SetBranchStatus("mu0_ProbNNk",1);
  fChain->SetBranchStatus("mu0_ProbNNp",1);
  fChain->SetBranchStatus("mu0_ProbNNpi",1);
  fChain->SetBranchStatus("mu0_ProbNNmu",1);
  fChain->SetBranchStatus("mu0_ProbNNghost",1);

  fChain->SetBranchStatus("mu1_PIDe",1);
  fChain->SetBranchStatus("mu1_PIDmu",1);
  fChain->SetBranchStatus("mu1_PIDK",1);
  fChain->SetBranchStatus("mu1_PIDp",1);
  fChain->SetBranchStatus("mu1_ProbNNe",1);
  fChain->SetBranchStatus("mu1_ProbNNk",1);
  fChain->SetBranchStatus("mu1_ProbNNp",1);
  fChain->SetBranchStatus("mu1_ProbNNpi",1);
  fChain->SetBranchStatus("mu1_ProbNNmu",1);
  fChain->SetBranchStatus("mu1_ProbNNghost",1);
   
  fChain->SetBranchStatus("h1_PIDe",1);
  fChain->SetBranchStatus("h1_PIDmu",1);
  fChain->SetBranchStatus("h1_PIDK",1);
  fChain->SetBranchStatus("h1_PIDp",1);
  fChain->SetBranchStatus("h1_ProbNNe",1);
  fChain->SetBranchStatus("h1_ProbNNk",1);
  fChain->SetBranchStatus("h1_ProbNNp",1);
  fChain->SetBranchStatus("h1_ProbNNpi",1);
  fChain->SetBranchStatus("h1_ProbNNmu",1);
  fChain->SetBranchStatus("h1_ProbNNghost",1);

  fChain->SetBranchStatus("h0_PIDe",1);
  fChain->SetBranchStatus("h0_PIDmu",1);
  fChain->SetBranchStatus("h0_PIDK",1);
  fChain->SetBranchStatus("h0_PIDp",1);
  fChain->SetBranchStatus("h0_ProbNNe",1);
  fChain->SetBranchStatus("h0_ProbNNk",1);
  fChain->SetBranchStatus("h0_ProbNNp",1);
  fChain->SetBranchStatus("h0_ProbNNpi",1);
  fChain->SetBranchStatus("h0_ProbNNmu",1);
  fChain->SetBranchStatus("h0_ProbNNghost",1);


  //trigger 

  fChain->SetBranchStatus("mu0_L0MuonDecision_TOS",1);
  fChain->SetBranchStatus("mu1_L0MuonDecision_TOS",1);
  fChain->SetBranchStatus("mu0_L0DiMuonDecision_TOS",1);
  fChain->SetBranchStatus("mu1_L0DiMuonDecision_TOS",1);
  fChain->SetBranchStatus("Dst_L0Global_TIS",1);
  fChain->SetBranchStatus("D_L0Global_TIS",1);

  fChain->SetBranchStatus("mu0_Hlt1TrackMuonDecision_TOS",1);
  fChain->SetBranchStatus("mu1_Hlt1TrackMuonDecision_TOS",1);
  fChain->SetBranchStatus("D_Hlt1TrackAllL0Decision_TOS",1);
  fChain->SetBranchStatus("D_Hlt1DiMuonHighMassDecision_TOS",1);
  fChain->SetBranchStatus("D_Hlt1DiMuonLowMassDecision_TOS",1);
  fChain->SetBranchStatus("mu0_Hlt1SingleMuonNoIPDecision_TOS",1);
  fChain->SetBranchStatus("mu1_Hlt1SingleMuonNoIPDecision_TOS",1);
  fChain->SetBranchStatus("mu0_Hlt1SingleMuonHighPTDecision_TOS",1);
  fChain->SetBranchStatus("mu1_Hlt1SingleMuonHighPTDecision_TOS",1);
  fChain->SetBranchStatus("D_Hlt2CharmSemilepD02KKMuMuDecision_TOS",1);
  fChain->SetBranchStatus("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS",1);
  fChain->SetBranchStatus("D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS",1);
  fChain->SetBranchStatus("Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS",1);
  fChain->SetBranchStatus("Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS",1);
  fChain->SetBranchStatus("D_Hlt2DiMuonDetachedDecision_TOS",1);

  fChain->SetBranchStatus("mu1_TRACK_GhostProb",1);
  fChain->SetBranchStatus("mu0_TRACK_GhostProb",1);
  fChain->SetBranchStatus("h0_TRACK_GhostProb",1);
  fChain->SetBranchStatus("h1_TRACK_GhostProb",1);
  fChain->SetBranchStatus("Slowpi_TRACK_GhostProb",1);


  //general 

  fChain->SetBranchStatus("nPVs",1);
  fChain->SetBranchStatus("nPV",1);
  fChain->SetBranchStatus("PVNTRACKS",1);
  fChain->SetBranchStatus("PVCHI2",1);
  fChain->SetBranchStatus("PVNDOF",1);
  fChain->SetBranchStatus("nVELO",1);
  fChain->SetBranchStatus("nSpdDigits",1);
  fChain->SetBranchStatus("nTracks",1);
  fChain->SetBranchStatus("nTTClusters",1);
  fChain->SetBranchStatus("nITClusters",1);
  //fChain->SetBranchStatus("nOTClusters",1); //not in trees!

  fChain->SetBranchStatus("Dst_OWNPV_X",1);
  fChain->SetBranchStatus("Dst_OWNPV_Y",1);
  fChain->SetBranchStatus("Dst_OWNPV_Z",1);


  //IP, pointing

  fChain->SetBranchStatus("Dst_ENDVERTEX_CHI2",1);
  fChain->SetBranchStatus("D_ENDVERTEX_CHI2",1);

  fChain->SetBranchStatus("Dst_MINIP",1);
  fChain->SetBranchStatus("Dst_MINIPCHI2",1);
  
  fChain->SetBranchStatus("D_DIRA_OWNPV",1);
  fChain->SetBranchStatus("D_DIRA_ORIVX",1);
  fChain->SetBranchStatus("Dst_DIRA_OWNPV",1);
  
  fChain->SetBranchStatus("D_MAXDOCA",1);
  fChain->SetBranchStatus("Dst_MAXDOCA",1);

  fChain->SetBranchStatus("D_IP_ORIVX",1);
  fChain->SetBranchStatus("D_IPCHI2_ORIVX",1);
  fChain->SetBranchStatus("D_IP_OWNPV",1);
  fChain->SetBranchStatus("D_IPCHI2_OWNPV",1);
  fChain->SetBranchStatus("D_MINIP",1); 
  fChain->SetBranchStatus("D_MINIPCHI2",1);
  
  fChain->SetBranchStatus("D_FD_ORIVX",1); 
  fChain->SetBranchStatus("D_FDCHI2_ORIVX",1);
  fChain->SetBranchStatus("D_FD_OWNPV",1); 
  fChain->SetBranchStatus("D_FDCHI2_OWNPV",1);
  
  fChain->SetBranchStatus("Dst_FD_OWNPV",1); 
  fChain->SetBranchStatus("Dst_FDCHI2_OWNPV",1);

  fChain->SetBranchStatus("Dst_IP_OWNPV",1);
  fChain->SetBranchStatus("Dst_IPCHI2_OWNPV",1);
  fChain->SetBranchStatus("Dst_MINIP",1); 
  fChain->SetBranchStatus("Dst_MINIPCHI2",1);
  
  fChain->SetBranchStatus("Slowpi_isMuon",1);
  fChain->SetBranchStatus("Slowpi_isMuonLoose",1);
  fChain->SetBranchStatus("Slowpi_IP_ORIVX",1);
  fChain->SetBranchStatus("Slowpi_IPCHI2_ORIVX",1);
  fChain->SetBranchStatus("Slowpi_IP_OWNPV",1);
  fChain->SetBranchStatus("Slowpi_IPCHI2_OWNPV",1);
  fChain->SetBranchStatus("Slowpi_MINIP",1); 
  fChain->SetBranchStatus("Slowpi_MINIPCHI2",1);
  fChain->SetBranchStatus("Slowpi_NShared",1);
 
  fChain->SetBranchStatus("mu0_isMuon",1);
  fChain->SetBranchStatus("mu0_isMuonLoose",1);
  fChain->SetBranchStatus("mu0_IP_ORIVX",1);
  fChain->SetBranchStatus("mu0_IPCHI2_ORIVX",1);
  fChain->SetBranchStatus("mu0_IP_OWNPV",1);
  fChain->SetBranchStatus("mu0_IPCHI2_OWNPV",1);
  fChain->SetBranchStatus("mu0_MINIP",1); 
  fChain->SetBranchStatus("mu0_MINIPCHI2",1);
  fChain->SetBranchStatus("mu0_NShared",1);
  fChain->SetBranchStatus("mu0_MuonNShared",1);

  fChain->SetBranchStatus("mu1_isMuon",1);
  fChain->SetBranchStatus("mu1_isMuonLoose",1);
  fChain->SetBranchStatus("mu1_IP_ORIVX",1);
  fChain->SetBranchStatus("mu1_IPCHI2_ORIVX",1);
  fChain->SetBranchStatus("mu1_IP_OWNPV",1);
  fChain->SetBranchStatus("mu1_IPCHI2_OWNPV",1);
  fChain->SetBranchStatus("mu1_MINIP",1); 
  fChain->SetBranchStatus("mu1_MINIPCHI2",1);
  fChain->SetBranchStatus("mu1_NShared",1);
  fChain->SetBranchStatus("mu1_MuonNShared",1);

  fChain->SetBranchStatus("h1_isMuon",1);
  fChain->SetBranchStatus("h1_isMuonLoose",1);
  fChain->SetBranchStatus("h1_IP_ORIVX",1);
  fChain->SetBranchStatus("h1_IPCHI2_ORIVX",1);
  fChain->SetBranchStatus("h1_IP_OWNPV",1);
  fChain->SetBranchStatus("h1_IPCHI2_OWNPV",1);
  fChain->SetBranchStatus("h1_MINIP",1); 
  fChain->SetBranchStatus("h1_MINIPCHI2",1);
  fChain->SetBranchStatus("h1_NShared",1);
 
  fChain->SetBranchStatus("h0_isMuon",1);
  fChain->SetBranchStatus("h0_isMuonLoose",1);
  fChain->SetBranchStatus("h0_IP_ORIVX",1);
  fChain->SetBranchStatus("h0_IPCHI2_ORIVX",1);
  fChain->SetBranchStatus("h0_IP_OWNPV",1);
  fChain->SetBranchStatus("h0_IPCHI2_OWNPV",1);
  fChain->SetBranchStatus("h0_MINIP",1); 
  fChain->SetBranchStatus("h0_MINIPCHI2",1);
  fChain->SetBranchStatus("h0_NShared",1);

  fChain->SetBranchStatus("Dst_CONEANGLE_D",1);
  fChain->SetBranchStatus("Dst_CONEANGLE_Dstar",1);
  fChain->SetBranchStatus("Dst_CONEMULT_D",1);
  fChain->SetBranchStatus("Dst_CONEMULT_Dstar",1);
  fChain->SetBranchStatus("Dst_CONEPTASYM_D",1);
  fChain->SetBranchStatus("Dst_CONEPTASYM_Dstar",1);

  fChain->SetBranchStatus("mu1_TRACK_CHI2NDOF",1);
  fChain->SetBranchStatus("mu0_TRACK_CHI2NDOF",1);
  fChain->SetBranchStatus("h1_TRACK_CHI2NDOF",1);
  fChain->SetBranchStatus("h0_TRACK_CHI2NDOF",1);
  fChain->SetBranchStatus("Slowpi_TRACK_CHI2NDOF",1);
  fChain->SetBranchStatus("eventNumber",1);

  Notify();
  

}

void D2hhmumuReader::InitMC()
{
 
   fChain->SetBranchAddress("Dst_BKGCAT", &Dst_BKGCAT, &b_Dst_BKGCAT);
   fChain->SetBranchAddress("Dst_TRUEID", &Dst_TRUEID, &b_Dst_TRUEID);
   fChain->SetBranchAddress("Dst_MC_MOTHER_ID", &Dst_MC_MOTHER_ID, &b_Dst_MC_MOTHER_ID);
   fChain->SetBranchAddress("Dst_MC_MOTHER_KEY", &Dst_MC_MOTHER_KEY, &b_Dst_MC_MOTHER_KEY);
   fChain->SetBranchAddress("Dst_MC_GD_MOTHER_ID", &Dst_MC_GD_MOTHER_ID, &b_Dst_MC_GD_MOTHER_ID);
   fChain->SetBranchAddress("Dst_MC_GD_MOTHER_KEY", &Dst_MC_GD_MOTHER_KEY, &b_Dst_MC_GD_MOTHER_KEY);
   fChain->SetBranchAddress("Dst_MC_GD_GD_MOTHER_ID", &Dst_MC_GD_GD_MOTHER_ID, &b_Dst_MC_GD_GD_MOTHER_ID);
   fChain->SetBranchAddress("Dst_MC_GD_GD_MOTHER_KEY", &Dst_MC_GD_GD_MOTHER_KEY, &b_Dst_MC_GD_GD_MOTHER_KEY);
   fChain->SetBranchAddress("Dst_TRUEP_E", &Dst_TRUEP_E, &b_Dst_TRUEP_E);
   fChain->SetBranchAddress("Dst_TRUEP_X", &Dst_TRUEP_X, &b_Dst_TRUEP_X);
   fChain->SetBranchAddress("Dst_TRUEP_Y", &Dst_TRUEP_Y, &b_Dst_TRUEP_Y);
   fChain->SetBranchAddress("Dst_TRUEP_Z", &Dst_TRUEP_Z, &b_Dst_TRUEP_Z);
   fChain->SetBranchAddress("Dst_TRUEPT", &Dst_TRUEPT, &b_Dst_TRUEPT);
   fChain->SetBranchAddress("Dst_TRUEORIGINVERTEX_X", &Dst_TRUEORIGINVERTEX_X, &b_Dst_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("Dst_TRUEORIGINVERTEX_Y", &Dst_TRUEORIGINVERTEX_Y, &b_Dst_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("Dst_TRUEORIGINVERTEX_Z", &Dst_TRUEORIGINVERTEX_Z, &b_Dst_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("Dst_TRUEENDVERTEX_X", &Dst_TRUEENDVERTEX_X, &b_Dst_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("Dst_TRUEENDVERTEX_Y", &Dst_TRUEENDVERTEX_Y, &b_Dst_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("Dst_TRUEENDVERTEX_Z", &Dst_TRUEENDVERTEX_Z, &b_Dst_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("Dst_TRUEISSTABLE", &Dst_TRUEISSTABLE, &b_Dst_TRUEISSTABLE);
   fChain->SetBranchAddress("Dst_TRUETAU", &Dst_TRUETAU, &b_Dst_TRUETAU);
   fChain->SetBranchAddress("D_BKGCAT", &D_BKGCAT, &b_D_BKGCAT);
   fChain->SetBranchAddress("D_TRUEID", &D_TRUEID, &b_D_TRUEID);
   fChain->SetBranchAddress("D_MC_MOTHER_ID", &D_MC_MOTHER_ID, &b_D_MC_MOTHER_ID);
   fChain->SetBranchAddress("D_MC_MOTHER_KEY", &D_MC_MOTHER_KEY, &b_D_MC_MOTHER_KEY);
   fChain->SetBranchAddress("D_MC_GD_MOTHER_ID", &D_MC_GD_MOTHER_ID, &b_D_MC_GD_MOTHER_ID);
   fChain->SetBranchAddress("D_MC_GD_MOTHER_KEY", &D_MC_GD_MOTHER_KEY, &b_D_MC_GD_MOTHER_KEY);
   fChain->SetBranchAddress("D_MC_GD_GD_MOTHER_ID", &D_MC_GD_GD_MOTHER_ID, &b_D_MC_GD_GD_MOTHER_ID);
   fChain->SetBranchAddress("D_MC_GD_GD_MOTHER_KEY", &D_MC_GD_GD_MOTHER_KEY, &b_D_MC_GD_GD_MOTHER_KEY);
   fChain->SetBranchAddress("D_TRUEP_E", &D_TRUEP_E, &b_D_TRUEP_E);
   fChain->SetBranchAddress("D_TRUEP_X", &D_TRUEP_X, &b_D_TRUEP_X);
   fChain->SetBranchAddress("D_TRUEP_Y", &D_TRUEP_Y, &b_D_TRUEP_Y);
   fChain->SetBranchAddress("D_TRUEP_Z", &D_TRUEP_Z, &b_D_TRUEP_Z);
   fChain->SetBranchAddress("D_TRUEPT", &D_TRUEPT, &b_D_TRUEPT);
   fChain->SetBranchAddress("D_TRUEORIGINVERTEX_X", &D_TRUEORIGINVERTEX_X, &b_D_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("D_TRUEORIGINVERTEX_Y", &D_TRUEORIGINVERTEX_Y, &b_D_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("D_TRUEORIGINVERTEX_Z", &D_TRUEORIGINVERTEX_Z, &b_D_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("D_TRUEENDVERTEX_X", &D_TRUEENDVERTEX_X, &b_D_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("D_TRUEENDVERTEX_Y", &D_TRUEENDVERTEX_Y, &b_D_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("D_TRUEENDVERTEX_Z", &D_TRUEENDVERTEX_Z, &b_D_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("D_TRUEISSTABLE", &D_TRUEISSTABLE, &b_D_TRUEISSTABLE);
   fChain->SetBranchAddress("D_TRUETAU", &D_TRUETAU, &b_D_TRUETAU);
   fChain->SetBranchAddress("h0_TRUEID", &h0_TRUEID, &b_h0_TRUEID);
   fChain->SetBranchAddress("h0_MC_MOTHER_ID", &h0_MC_MOTHER_ID, &b_h0_MC_MOTHER_ID);
   fChain->SetBranchAddress("h0_MC_MOTHER_KEY", &h0_MC_MOTHER_KEY, &b_h0_MC_MOTHER_KEY);
   fChain->SetBranchAddress("h0_MC_GD_MOTHER_ID", &h0_MC_GD_MOTHER_ID, &b_h0_MC_GD_MOTHER_ID);
   fChain->SetBranchAddress("h0_MC_GD_MOTHER_KEY", &h0_MC_GD_MOTHER_KEY, &b_h0_MC_GD_MOTHER_KEY);
   fChain->SetBranchAddress("h0_MC_GD_GD_MOTHER_ID", &h0_MC_GD_GD_MOTHER_ID, &b_h0_MC_GD_GD_MOTHER_ID);
   fChain->SetBranchAddress("h0_MC_GD_GD_MOTHER_KEY", &h0_MC_GD_GD_MOTHER_KEY, &b_h0_MC_GD_GD_MOTHER_KEY);
   fChain->SetBranchAddress("h0_TRUEP_E", &h0_TRUEP_E, &b_h0_TRUEP_E);
   fChain->SetBranchAddress("h0_TRUEP_X", &h0_TRUEP_X, &b_h0_TRUEP_X);
   fChain->SetBranchAddress("h0_TRUEP_Y", &h0_TRUEP_Y, &b_h0_TRUEP_Y);
   fChain->SetBranchAddress("h0_TRUEP_Z", &h0_TRUEP_Z, &b_h0_TRUEP_Z);
   fChain->SetBranchAddress("h0_TRUEPT", &h0_TRUEPT, &b_h0_TRUEPT);
   fChain->SetBranchAddress("h0_TRUEORIGINVERTEX_X", &h0_TRUEORIGINVERTEX_X, &b_h0_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("h0_TRUEORIGINVERTEX_Y", &h0_TRUEORIGINVERTEX_Y, &b_h0_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("h0_TRUEORIGINVERTEX_Z", &h0_TRUEORIGINVERTEX_Z, &b_h0_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("h0_TRUEENDVERTEX_X", &h0_TRUEENDVERTEX_X, &b_h0_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("h0_TRUEENDVERTEX_Y", &h0_TRUEENDVERTEX_Y, &b_h0_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("h0_TRUEENDVERTEX_Z", &h0_TRUEENDVERTEX_Z, &b_h0_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("h0_TRUEISSTABLE", &h0_TRUEISSTABLE, &b_h0_TRUEISSTABLE);
   fChain->SetBranchAddress("h0_TRUETAU", &h0_TRUETAU, &b_h0_TRUETAU);
   fChain->SetBranchAddress("h1_TRUEID", &h1_TRUEID, &b_h1_TRUEID);
   fChain->SetBranchAddress("h1_MC_MOTHER_ID", &h1_MC_MOTHER_ID, &b_h1_MC_MOTHER_ID);
   fChain->SetBranchAddress("h1_MC_MOTHER_KEY", &h1_MC_MOTHER_KEY, &b_h1_MC_MOTHER_KEY);
   fChain->SetBranchAddress("h1_MC_GD_MOTHER_ID", &h1_MC_GD_MOTHER_ID, &b_h1_MC_GD_MOTHER_ID);
   fChain->SetBranchAddress("h1_MC_GD_MOTHER_KEY", &h1_MC_GD_MOTHER_KEY, &b_h1_MC_GD_MOTHER_KEY);
   fChain->SetBranchAddress("h1_MC_GD_GD_MOTHER_ID", &h1_MC_GD_GD_MOTHER_ID, &b_h1_MC_GD_GD_MOTHER_ID);
   fChain->SetBranchAddress("h1_MC_GD_GD_MOTHER_KEY", &h1_MC_GD_GD_MOTHER_KEY, &b_h1_MC_GD_GD_MOTHER_KEY);
   fChain->SetBranchAddress("h1_TRUEP_E", &h1_TRUEP_E, &b_h1_TRUEP_E);
   fChain->SetBranchAddress("h1_TRUEP_X", &h1_TRUEP_X, &b_h1_TRUEP_X);
   fChain->SetBranchAddress("h1_TRUEP_Y", &h1_TRUEP_Y, &b_h1_TRUEP_Y);
   fChain->SetBranchAddress("h1_TRUEP_Z", &h1_TRUEP_Z, &b_h1_TRUEP_Z);
   fChain->SetBranchAddress("h1_TRUEPT", &h1_TRUEPT, &b_h1_TRUEPT);
   fChain->SetBranchAddress("h1_TRUEORIGINVERTEX_X", &h1_TRUEORIGINVERTEX_X, &b_h1_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("h1_TRUEORIGINVERTEX_Y", &h1_TRUEORIGINVERTEX_Y, &b_h1_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("h1_TRUEORIGINVERTEX_Z", &h1_TRUEORIGINVERTEX_Z, &b_h1_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("h1_TRUEENDVERTEX_X", &h1_TRUEENDVERTEX_X, &b_h1_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("h1_TRUEENDVERTEX_Y", &h1_TRUEENDVERTEX_Y, &b_h1_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("h1_TRUEENDVERTEX_Z", &h1_TRUEENDVERTEX_Z, &b_h1_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("h1_TRUEISSTABLE", &h1_TRUEISSTABLE, &b_h1_TRUEISSTABLE);
   fChain->SetBranchAddress("h1_TRUETAU", &h1_TRUETAU, &b_h1_TRUETAU);
   fChain->SetBranchAddress("mu0_TRUEID", &mu0_TRUEID, &b_mu0_TRUEID);
   fChain->SetBranchAddress("mu0_MC_MOTHER_ID", &mu0_MC_MOTHER_ID, &b_mu0_MC_MOTHER_ID);
   fChain->SetBranchAddress("mu0_MC_MOTHER_KEY", &mu0_MC_MOTHER_KEY, &b_mu0_MC_MOTHER_KEY);
   fChain->SetBranchAddress("mu0_MC_GD_MOTHER_ID", &mu0_MC_GD_MOTHER_ID, &b_mu0_MC_GD_MOTHER_ID);
   fChain->SetBranchAddress("mu0_MC_GD_MOTHER_KEY", &mu0_MC_GD_MOTHER_KEY, &b_mu0_MC_GD_MOTHER_KEY);
   fChain->SetBranchAddress("mu0_MC_GD_GD_MOTHER_ID", &mu0_MC_GD_GD_MOTHER_ID, &b_mu0_MC_GD_GD_MOTHER_ID);
   fChain->SetBranchAddress("mu0_MC_GD_GD_MOTHER_KEY", &mu0_MC_GD_GD_MOTHER_KEY, &b_mu0_MC_GD_GD_MOTHER_KEY);
   fChain->SetBranchAddress("mu0_TRUEP_E", &mu0_TRUEP_E, &b_mu0_TRUEP_E);
   fChain->SetBranchAddress("mu0_TRUEP_X", &mu0_TRUEP_X, &b_mu0_TRUEP_X);
   fChain->SetBranchAddress("mu0_TRUEP_Y", &mu0_TRUEP_Y, &b_mu0_TRUEP_Y);
   fChain->SetBranchAddress("mu0_TRUEP_Z", &mu0_TRUEP_Z, &b_mu0_TRUEP_Z);
   fChain->SetBranchAddress("mu0_TRUEPT", &mu0_TRUEPT, &b_mu0_TRUEPT);
   fChain->SetBranchAddress("mu0_TRUEORIGINVERTEX_X", &mu0_TRUEORIGINVERTEX_X, &b_mu0_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("mu0_TRUEORIGINVERTEX_Y", &mu0_TRUEORIGINVERTEX_Y, &b_mu0_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("mu0_TRUEORIGINVERTEX_Z", &mu0_TRUEORIGINVERTEX_Z, &b_mu0_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("mu0_TRUEENDVERTEX_X", &mu0_TRUEENDVERTEX_X, &b_mu0_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("mu0_TRUEENDVERTEX_Y", &mu0_TRUEENDVERTEX_Y, &b_mu0_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("mu0_TRUEENDVERTEX_Z", &mu0_TRUEENDVERTEX_Z, &b_mu0_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("mu0_TRUEISSTABLE", &mu0_TRUEISSTABLE, &b_mu0_TRUEISSTABLE);
   fChain->SetBranchAddress("mu0_TRUETAU", &mu0_TRUETAU, &b_mu0_TRUETAU);
   fChain->SetBranchAddress("mu1_TRUEID", &mu1_TRUEID, &b_mu1_TRUEID);
   fChain->SetBranchAddress("mu1_MC_MOTHER_ID", &mu1_MC_MOTHER_ID, &b_mu1_MC_MOTHER_ID);
   fChain->SetBranchAddress("mu1_MC_MOTHER_KEY", &mu1_MC_MOTHER_KEY, &b_mu1_MC_MOTHER_KEY);
   fChain->SetBranchAddress("mu1_MC_GD_MOTHER_ID", &mu1_MC_GD_MOTHER_ID, &b_mu1_MC_GD_MOTHER_ID);
   fChain->SetBranchAddress("mu1_MC_GD_MOTHER_KEY", &mu1_MC_GD_MOTHER_KEY, &b_mu1_MC_GD_MOTHER_KEY);
   fChain->SetBranchAddress("mu1_MC_GD_GD_MOTHER_ID", &mu1_MC_GD_GD_MOTHER_ID, &b_mu1_MC_GD_GD_MOTHER_ID);
   fChain->SetBranchAddress("mu1_MC_GD_GD_MOTHER_KEY", &mu1_MC_GD_GD_MOTHER_KEY, &b_mu1_MC_GD_GD_MOTHER_KEY);
   fChain->SetBranchAddress("mu1_TRUEP_E", &mu1_TRUEP_E, &b_mu1_TRUEP_E);
   fChain->SetBranchAddress("mu1_TRUEP_X", &mu1_TRUEP_X, &b_mu1_TRUEP_X);
   fChain->SetBranchAddress("mu1_TRUEP_Y", &mu1_TRUEP_Y, &b_mu1_TRUEP_Y);
   fChain->SetBranchAddress("mu1_TRUEP_Z", &mu1_TRUEP_Z, &b_mu1_TRUEP_Z);
   fChain->SetBranchAddress("mu1_TRUEPT", &mu1_TRUEPT, &b_mu1_TRUEPT);
   fChain->SetBranchAddress("mu1_TRUEORIGINVERTEX_X", &mu1_TRUEORIGINVERTEX_X, &b_mu1_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("mu1_TRUEORIGINVERTEX_Y", &mu1_TRUEORIGINVERTEX_Y, &b_mu1_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("mu1_TRUEORIGINVERTEX_Z", &mu1_TRUEORIGINVERTEX_Z, &b_mu1_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("mu1_TRUEENDVERTEX_X", &mu1_TRUEENDVERTEX_X, &b_mu1_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("mu1_TRUEENDVERTEX_Y", &mu1_TRUEENDVERTEX_Y, &b_mu1_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("mu1_TRUEENDVERTEX_Z", &mu1_TRUEENDVERTEX_Z, &b_mu1_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("mu1_TRUEISSTABLE", &mu1_TRUEISSTABLE, &b_mu1_TRUEISSTABLE);
   fChain->SetBranchAddress("mu1_TRUETAU", &mu1_TRUETAU, &b_mu1_TRUETAU);
   fChain->SetBranchAddress("Slowpi_TRUEID", &Slowpi_TRUEID, &b_Slowpi_TRUEID);
   fChain->SetBranchAddress("Slowpi_MC_MOTHER_ID", &Slowpi_MC_MOTHER_ID, &b_Slowpi_MC_MOTHER_ID);
   fChain->SetBranchAddress("Slowpi_MC_MOTHER_KEY", &Slowpi_MC_MOTHER_KEY, &b_Slowpi_MC_MOTHER_KEY);
   fChain->SetBranchAddress("Slowpi_MC_GD_MOTHER_ID", &Slowpi_MC_GD_MOTHER_ID, &b_Slowpi_MC_GD_MOTHER_ID);
   fChain->SetBranchAddress("Slowpi_MC_GD_MOTHER_KEY", &Slowpi_MC_GD_MOTHER_KEY, &b_Slowpi_MC_GD_MOTHER_KEY);
   fChain->SetBranchAddress("Slowpi_MC_GD_GD_MOTHER_ID", &Slowpi_MC_GD_GD_MOTHER_ID, &b_Slowpi_MC_GD_GD_MOTHER_ID);
   fChain->SetBranchAddress("Slowpi_MC_GD_GD_MOTHER_KEY", &Slowpi_MC_GD_GD_MOTHER_KEY, &b_Slowpi_MC_GD_GD_MOTHER_KEY);
   fChain->SetBranchAddress("Slowpi_TRUEP_E", &Slowpi_TRUEP_E, &b_Slowpi_TRUEP_E);
   fChain->SetBranchAddress("Slowpi_TRUEP_X", &Slowpi_TRUEP_X, &b_Slowpi_TRUEP_X);
   fChain->SetBranchAddress("Slowpi_TRUEP_Y", &Slowpi_TRUEP_Y, &b_Slowpi_TRUEP_Y);
   fChain->SetBranchAddress("Slowpi_TRUEP_Z", &Slowpi_TRUEP_Z, &b_Slowpi_TRUEP_Z);
   fChain->SetBranchAddress("Slowpi_TRUEPT", &Slowpi_TRUEPT, &b_Slowpi_TRUEPT);
   fChain->SetBranchAddress("Slowpi_TRUEORIGINVERTEX_X", &Slowpi_TRUEORIGINVERTEX_X, &b_Slowpi_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("Slowpi_TRUEORIGINVERTEX_Y", &Slowpi_TRUEORIGINVERTEX_Y, &b_Slowpi_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("Slowpi_TRUEORIGINVERTEX_Z", &Slowpi_TRUEORIGINVERTEX_Z, &b_Slowpi_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("Slowpi_TRUEENDVERTEX_X", &Slowpi_TRUEENDVERTEX_X, &b_Slowpi_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("Slowpi_TRUEENDVERTEX_Y", &Slowpi_TRUEENDVERTEX_Y, &b_Slowpi_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("Slowpi_TRUEENDVERTEX_Z", &Slowpi_TRUEENDVERTEX_Z, &b_Slowpi_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("Slowpi_TRUEISSTABLE", &Slowpi_TRUEISSTABLE, &b_Slowpi_TRUEISSTABLE);
   fChain->SetBranchAddress("Slowpi_TRUETAU", &Slowpi_TRUETAU, &b_Slowpi_TRUETAU);
}


Bool_t D2hhmumuReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void D2hhmumuReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t D2hhmumuReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//#endif // #ifdef D2hhmumuReader_cxx


