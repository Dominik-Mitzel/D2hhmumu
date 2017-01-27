#include "optimizeSelection.h"
#include "D2hhmumuFitter.h"
#include "TEventList.h"
#include "TPaletteAxis.h"
#include "D2KKmumuReader.h"
#include "D2pipimumuReader.h"
#include "D2KpimumuReader.h"
#include "D2KKpipiReader.h"
#include "D2KpipipiReader.h"
#include "D2pipipipiReader.h"
#include "TProfile.h"
#include "D2hhmumuFitter1D.h"
#include "TGenPhaseSpace.h"
#include "EfficiencyCalculator.h"
#include "TEfficiency.h"
#include "D2hhmumuFitter1D.h"
#include "RooGenericPdf.h"
#include "EfficiencyStudies.h"
#include <utility> 
#include <fstream>
#include "TList.h"
#include "TMVA_applications.h"
#include "dcastyle.h"

using namespace std;
using namespace RooFit ;

//set the binning scheme global


  std::vector<double> rangesKpi_low = {675};
  std::vector<double> rangesKpi_high = {875};
  std::vector<double> rangespipi_low = {200,525,565,950,1100};
  std::vector<double> rangespipi_high = {525,565,950,1100,1600};
  std::vector<double> rangesKK_low = {200,525,565};
  std::vector<double> rangesKK_high = {525,565,900};//!!//THIS IS THE NEW BINNING



  double binsKpi[2]={675,875};
  double binsKK[4]={200,525,565,900};
  double binspipi[6]={200,525,565,950,1100,1600};

/*

std::vector<double> rangesKpi_low = {675};
std::vector<double> rangesKpi_high = {875};
std::vector<double> rangespipi_low = {0,525,565,950,1100};
std::vector<double> rangespipi_high = {525,565,950,1100,1500};
recoandstd::vector<double> rangesKK_low = {0,525,565};
std::vector<double> rangesKK_high = {525,565,950};
 


double binsKpi[2]={675,875};
double binsKK[4]={0,525,565,950};
double binspipi[6]={0,525,565,950,1100,1500};

*/

/*
  std::vector<double> rangesKpi_low = {675};
  std::vector<double> rangesKpi_high = {875};
  std::vector<double> rangespipi_low = {0};
  std::vector<double> rangespipi_high = {1500};
//std::vector<double> rangesKK_low = {0,525};
// std::vector<double> rangesKK_high = {525,565};//!!//950!!
  std::vector<double> rangesKK_low = {0};
  std::vector<double> rangesKK_high = {950};//!!//950!!

  double binsKpi[2]={675,875};
  double binsKK[2]={0,950};
  double binspipi[2]={0,1500};
*/

///std::vector<double> rangesPt_low = {800,1200,1600,2000,4000};
///std::vector<double> rangesPt_high = {1200,1600,2000,4000,8000};
///double binsPt[6]={800,1200,1600,2000,4000,8000};
/*///!!!fineBinning!!///*/std::vector<double> rangesPt_low = {0,800,1100,1400,1700,2200,3000};
/*///!!!fineBinning!!///*/std::vector<double> rangesPt_high = {800,1100,1400,1700,2200,3000,8000};
/*///!!!fineBinning!!///*/double binsPt[]={0,800,1100,1400,1700,2200,3000,8000};


/*///!!!fineBinning!!///*/  std::vector<double> rangesPt_low_fine = {0,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500,4000,5000};
/*///!!!fineBinning!!///*/  std::vector<double> rangesPt_high_fine = {200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500,4000,5000,8000};
/*///!!!fineBinning!!///*/double binsPt_fine[]={0,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3500,4000,5000,8000};


//std::vector<double> rangesPt_low = {0,800,1400,2200,3000};
//std::vector<double> rangesPt_high = {800,1400,2200,3000,8000};
//double binsPt[]={0,800,1400,2200,3000,8000};

std::vector<double> rangesPtBothMuons_low = {0,800};
std::vector<double> rangesPtBothMuons_high = {800,8000};
double binsPtBothMuons[]={0,800,8000};

//full range including the 0 
//std::vector<double> full_rangesPt_low = {0,800,1100,1400,1700,2200,3000};
//std::vector<double> full_rangesPt_high = {800,1100,1400,1700,2200,3000,8000};
//double full_binsPt[]={0,800,1100,1400,1700,2200,3000,8000};

void performFullPIDEfficiencyStudy(){

/////////////////////////////////
/*
  Here, the PID efficiency is determined using data driven methods.
  The PID calib muon reference sample does not cover the phase space
  of the low pt momentum of signal and reference decay (pt cut >800MeV)

  PREPARATIONS
  
  check if the kinmetics of the PID calib reference sample and the signal samples match
*/
//  plotPIDCalibMuonSampleAndSignalMC();
//  plotPIDCalibKaonSampleAndSignalMC();
//  plotPIDCalibPionSampleAndSignalMC();
 /*
  split the MC tuples in the dimuon mass bins with:
*/
  //createMCTuplesForEffStudies(); 
  //CrossapplicationForEfficiencyStudies(); //apply BDT to the splitted tuples 
/*
  then, split the tuples (already BDT variable added to them) by the pt of the muon:
*/
  //splitMCtuplesInPt(true); //bins in pt ot mu0,argument refers to cut un muonNshared
  //here, the dataset is split accoring to rangesPTBothMuons, which is basically just cutting the low pt regime
  //splitMCtuplesPtCutOnSingleMuon(true,0);
  //splitMCtuplesPtCutOnSingleMuon(true,1);
  //splitMCtuplesInPtOfBothMuons(true); //switched of now, has to be done once

/*
  evaluate the muon PID eff in bins of muon pt for each polarity. Beforehand, the PID calibs scripts located in ~/cmtuser/../PIDCALIB/../D2hhmumu/ have to be executed.
  ///////////////////////////////////////////////////////
  HOWTO EXECUTE THE PIDCALIB SCRIPTS
  
  Setup environment: LbLogin -c x86_64-slc6-gcc48-opt
                     SetupUrania v5r0
		    
  needed scripts: in /home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu
                  create the tables for muon and hadron PID cut:
		  
		  make the tables for a certain PID cut:

		  ./createTablesMuonPID.sh
		  ./createTablesHadronPID.sh
	  
		  perform the scripts using the tables and signal MC
		  ./performDoubeMuonPID.py  ( for the momentum range 800-8000MeV )
		  ./performMuonSingePID.py  ( for the momentum range 800-8000MeV )
		  ./performMuonPIDInPTBins.py (in the defined pt bins for single muon, defined in performMuonPIDInPTBins.py. The binning
		                               scheme must match the one defined here in EfficiencyStudy.cpp)

		  for hadrons:
		  ./performHadronPID.py
		  
		  the PID calib results are saved in folders in /home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu and will be read out by the
		  following programms

  //////////////////////////////////////////////////
		  
  Note: in evaluatePIDCalibMuonEfficiencyForPTBins, also the extrapolated low pt efficiency is computed by a linear fit! It is done in bin of dimuon mass bins and only for muon index 0
*/

  evaluatePIDCalibMuonEfficiencyForPTBins("magDw");
  evaluatePIDCalibMuonEfficiencyForPTBins("magUp");
 /*
   Evaluate the single and double muon efficiencies for both polarities and muon indices in pt range > 800MeV
   PID calib has to be performed induvidually on magUp and Dw
 */
  

  evaluatePIDCalibSingleMuonEfficiency("magUp",0,"default"); //default is the standard PID calib binning scheme. 
  evaluatePIDCalibSingleMuonEfficiency("magDw",0,"default");
  evaluatePIDCalibSingleMuonEfficiency("magUp",1,"default");
  evaluatePIDCalibSingleMuonEfficiency("magDw",1,"default");
  
  evaluatePIDCalibDoubleMuonEfficiency("magDw");
  evaluatePIDCalibDoubleMuonEfficiency("magUp");

 // also get the MC single and double muon PID efficiencies. No splitting by magUp an Dw done
 
  evaluateMCSingleMuonEfficiency(0);
  evaluateMCSingleMuonEfficiency(1);
  evaluateMCDoubleMuonEfficiency();

 /*   
   as PID calib has to be done separately by polarity, the average is calculated  
  
 */

 //drawBothPolaritiesPIDCalibMuonEfficiencyForPTBins(); //write a NEW one that onl does the things we need right now!

 //computeGlobalLowPtEfficiency(); //compute the global low pt eff, averaged over polarity, channel, m(mumu)

 //average the PID calib results over polarity
 drawBothPolaritiesPIDCalibDoubleMuonEfficiency();

 drawBothPolaritiesPIDCalibSingleMuonEfficiency(0,"default");//again mouon index 0 and 1
 drawBothPolaritiesPIDCalibSingleMuonEfficiency(1,"default");
 
 //takes PID calib and MC PID efficiencies and corrects the low pt regime, therfor the mu_Pt spectra are needed
 //from MC
 
 //drawPtSpectra(0); //pt spectra for mu0 and mu1
 //drawPtSpectra(1);
 drawLowPtCorrectedEfficiency("default");
  
 //checkMuonPIDFactorization();
  
 /*                                                                                                                                                        
   for Hadron efficiency, no splitting by pt is needed, otherwise similar strategy. Phase space fully coverd by PID calib reference sample
   -evaluate PID calib efficiency for both polarities
   -average
   -evaluate MC eff
 */
 
  evaluatePIDCalibHadronEfficiency("magUp","default");
  evaluatePIDCalibHadronEfficiency("magDw","default");
  drawBothPolaritiesPIDCalibHadronEfficiency("default");
 
 //the total PID efficiency is the product of Hadron and muon PID
 
  drawTotalPIDEfficiency("default",0.015,0.01);

}

void evaluatePIDEfficiencySystematicUncertainty(){

  //evaluate the uncertainty due to the binning scheme in PID calib

  //muon PID
  
  evaluatePIDCalibSingleMuonEfficiency("magUp",0,"alternative");
  evaluatePIDCalibSingleMuonEfficiency("magDw",0,"alternative");
  evaluatePIDCalibSingleMuonEfficiency("magUp",1,"alternative");
  evaluatePIDCalibSingleMuonEfficiency("magDw",1,"alternative");

  evaluatePIDCalibSingleMuonEfficiency("magUp",0,"alternative2");
  evaluatePIDCalibSingleMuonEfficiency("magDw",0,"alternative2");
  evaluatePIDCalibSingleMuonEfficiency("magUp",1,"alternative2");
  evaluatePIDCalibSingleMuonEfficiency("magDw",1,"alternative2");

  drawBothPolaritiesPIDCalibSingleMuonEfficiency(0,"alternative");//again mouon index 0 and 1
  drawBothPolaritiesPIDCalibSingleMuonEfficiency(1,"alternative");

  drawBothPolaritiesPIDCalibSingleMuonEfficiency(0,"alternative2");//again mouon index 0 and 1
  drawBothPolaritiesPIDCalibSingleMuonEfficiency(1,"alternative2");

  drawLowPtCorrectedEfficiency("alternative");
  drawLowPtCorrectedEfficiency("alternative2");
    
  //Hadron PID

  evaluatePIDCalibHadronEfficiency("magUp","alternative");
  evaluatePIDCalibHadronEfficiency("magDw","alternative");

  evaluatePIDCalibHadronEfficiency("magUp","alternative2");
  evaluatePIDCalibHadronEfficiency("magDw","alternative2");

  drawBothPolaritiesPIDCalibHadronEfficiency("alternative");
  drawBothPolaritiesPIDCalibHadronEfficiency("alternative2");

  drawTotalPIDEfficiency("alternative",0.015,0.01);
  drawTotalPIDEfficiency("alternative2",0.015,0.01);
  
  compareHadonAndMuonEfficienciesWithDifferentBinning();

}

void compareHadonAndMuonEfficienciesWithDifferentBinning(){

  TFile* fDefaultMuon = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalMuonPIDEfficiency_default.root");
  TFile* fAlternativeMuon = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalMuonPIDEfficiency_alternative.root");
  TFile* fAlternative2Muon = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalMuonPIDEfficiency_alternative2.root");

  TFile* fDefaultHadron = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalHadronPIDEfficiency_default.root","OPEN");
  TFile* fAlternativeHadron = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalHadronPIDEfficiency_alternative.root","OPEN");
  TFile* fAlternative2Hadron = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalHadronPIDEfficiency_alternative2.root","OPEN");

  TH1D * relativeMuonPID_KK_default = (TH1D*)fDefaultMuon->Get("relEff_KKmumu");
  TH1D * relativeMuonPID_KK_alternative = (TH1D*)fAlternativeMuon->Get("relEff_KKmumu");
  TH1D * relativeMuonPID_KK_alternative2 = (TH1D*)fAlternative2Muon->Get("relEff_KKmumu");

  TH1D * relativeMuonPID_pipi_default = (TH1D*)fDefaultMuon->Get("relEff_pipimumu");
  TH1D * relativeMuonPID_pipi_alternative = (TH1D*)fAlternativeMuon->Get("relEff_pipimumu");
  TH1D * relativeMuonPID_pipi_alternative2 = (TH1D*)fAlternative2Muon->Get("relEff_pipimumu");

  TH1D * relativeHadronPID_KK_default = (TH1D*)fDefaultHadron->Get("relEff_KKmumu_average");
  TH1D * relativeHadronPID_KK_alternative = (TH1D*)fAlternativeHadron->Get("relEff_KKmumu_average");
  TH1D * relativeHadronPID_KK_alternative2 = (TH1D*)fAlternative2Hadron->Get("relEff_KKmumu_average");

  TH1D * relativeHadronPID_pipi_default = (TH1D*)fDefaultHadron->Get("relEff_pipimumu_average");
  TH1D * relativeHadronPID_pipi_alternative = (TH1D*)fAlternativeHadron->Get("relEff_pipimumu_average");
  TH1D * relativeHadronPID_pipi_alternative2 = (TH1D*)fAlternative2Hadron->Get("relEff_pipimumu_average");

  dcastyle();

    
  TCanvas* a = new TCanvas("a","a");
  a->Divide(1,2);
  a->cd(1);
  relativeMuonPID_KK_default->GetYaxis()->SetRangeUser(.9,1.2);
  relativeMuonPID_KK_default->Draw();
  relativeMuonPID_KK_alternative->SetLineColor(2);
  relativeMuonPID_KK_alternative->Draw("SAME");
  relativeMuonPID_KK_alternative2->SetLineColor(3);
  relativeMuonPID_KK_alternative2->Draw("SAME");
  cout<<"test1"<<endl;

  a->cd(2);
  relativeMuonPID_pipi_default->GetYaxis()->SetRangeUser(.9,1.2);
  relativeMuonPID_pipi_default->Draw();
  relativeMuonPID_pipi_alternative->SetLineColor(2);
  relativeMuonPID_pipi_alternative->Draw("SAME");
  relativeMuonPID_pipi_alternative2->SetLineColor(3);
  relativeMuonPID_pipi_alternative2->Draw("SAME");
  cout<<"test1"<<endl;
  a->Print("test1.eps");
    
  TCanvas*  b= new TCanvas("b","b");
  b->Divide(1,2);
  b->cd(1);
  relativeHadronPID_KK_default->GetYaxis()->SetRangeUser(.9,1.2);
  relativeHadronPID_KK_default->Draw();
  relativeHadronPID_KK_alternative->SetLineColor(2);
  relativeHadronPID_KK_alternative->Draw("SAME");
  relativeHadronPID_KK_alternative2->SetLineColor(3);
  relativeHadronPID_KK_alternative2->Draw("SAME");
  cout<<"test1"<<endl;

  b->cd(2);
  relativeHadronPID_pipi_default->GetYaxis()->SetRangeUser(.9,1.2);
  relativeHadronPID_pipi_default->Draw();
  relativeHadronPID_pipi_alternative->SetLineColor(2);
  relativeHadronPID_pipi_alternative->Draw("SAME");
  relativeHadronPID_pipi_alternative2->SetLineColor(3);
  relativeHadronPID_pipi_alternative2->Draw("SAME");
  cout<<"test1"<<endl;
  
  b->Print("test.eps");
    

}


class Cand {
public :
  ULong64_t eventNumber;
  UInt_t    runNumber;
  Double_t  mass;
  UInt_t    nCandidate;
  int       Dst_BKGCAT;
  Double_t  id() { return ((double) runNumber)/((double) eventNumber); }
};

class TrueCand {
public :
  ULong64_t eventNumber;
  UInt_t    runNumber;
  double mu0_PX, mu0_PY, mu0_PZ;
  double mu1_PX, mu1_PY, mu1_PZ;
  double h0_PX, h0_PY, h0_PZ;
  double h1_PX, h1_PY, h1_PZ;
  double Slowpi_PX,Slowpi_PY,Slowpi_PZ;
  Double_t  id() { return ((double) runNumber)/((double) eventNumber); }
};

 
//multiple candidate stuff

 
void writeMultipleCanidatesToFile(bool withGhosts){
  
  //write all the multiple canidated to a file with totCandidates>1 (after dm preselction and nShared cut for muons)
  

  TString nameKpi;
  TString nameKpi_D2KKmumuBDT;
  TString nameKK;
  TString namepipi;

  if(withGhosts) {
    nameKpi="candidateLists/candidatesKpi_2012_filtered_D2pipimumuBDT_withGhosts.txt";
    nameKpi_D2KKmumuBDT="candidateLists/candidatesKpi_2012_filtered_D2KKmumuBDT_withGhosts.txt";
    nameKK="candidateLists/candidatesKK_2012_filtered_withGhosts.txt";
    namepipi="candidateLists/candidatespipi_2012_filtered_withGhosts.txt";
  }

  else {
    nameKpi="candidateLists/candidatesKpi_2012_filtered_D2pipimumuBDT_noGhosts.txt";
    nameKpi_D2KKmumuBDT="candidateLists/candidatesKpi_2012_filtered_D2KKmumuBDT_noGhosts.txt";
    nameKK="candidateLists/candidatesKK_2012_filtered_noGhosts.txt";
    namepipi="candidateLists/candidatespipi_2012_filtered_noGhosts.txt";
  }

  
  std::ofstream fileKpi(nameKpi);
  std::ofstream fileKpi_D2KKmumuBDT(nameKpi_D2KKmumuBDT);
  std::ofstream fileKK(nameKK);
  std::ofstream filepipi(namepipi);

  double cutBDT=0.4;
  double cutPID=0.5; 

  //TChain* Chain_KK = new TChain("DecayTree");
  //TChain* Chain_Kpi=  new TChain("DecayTree");
  //TChain* Chain_pipi=  new TChain("DecayTree");
  TChain* Chain_KK = new TChain("BDT_Tree");
  TChain* Chain_Kpi=  new TChain("BDT_Tree");
  TChain* Chain_Kpi_D2KKmumuBDT=  new TChain("BDT_Tree");
  TChain* Chain_pipi=  new TChain("BDT_Tree");
  
  // TChain* Chain_pipi=  new TChain("MC12_DstD2pipiMuMu/DecayTree");

  //Chain_KK->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/MC12_DstD2KKmumu.root");
  //Chain_Kpi->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/MC12_DstD2Kpimumu.root");
  //Chain_pipi->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MC12_DstD2pipimumu.root");
  
  Chain_KK->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_200.0_1600.0_D2KKmumuBDT.root");
  Chain_Kpi->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2pipimumuBDT.root");
  Chain_pipi->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_200.0_1600.0_D2pipimumuBDT.root"); 
  Chain_Kpi_D2KKmumuBDT->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2KKmumuBDT.root");

  double m_D;
  ULong64_t eventNumber;
  UInt_t runNumber;
  UInt_t nCandidate;
  int Dst_BKGCAT;
  ULong64_t totCandidates;

  double mu0_PX, mu0_PY, mu0_PZ;
  double mu1_PX, mu1_PY, mu1_PZ;
  double h0_PX, h0_PY, h0_PZ;
  double h1_PX, h1_PY, h1_PZ;
  double Slowpi_PX,Slowpi_PY,Slowpi_PZ;
  double deltaM; 
  int mu1_MuonNShared,mu0_MuonNShared;

  double BDT, mu0_ProbNNmu, mu1_ProbNNmu;
  bool mu1_L0MuonDecision_TOS,mu0_L0MuonDecision_TOS;
  bool mu1_Hlt1TrackMuonDecision_TOS, mu0_Hlt1TrackMuonDecision_TOS, D_Hlt1TrackAllL0Decision_TOS;
  bool Hlt2_TOS;

  double h0_ProbNNghost, h1_ProbNNghost,mu0_ProbNNghost, mu1_ProbNNghost, Slowpi_ProbNNghost;
  double h0_ProbNNk,h1_ProbNNk,h0_ProbNNpi,h1_ProbNNpi;

  Chain_KK->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);
  Chain_KK->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
  Chain_KK->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
  Chain_KK->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
  Chain_KK->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
  Chain_KK->SetBranchAddress("h0_ProbNNk",&h0_ProbNNk);
  Chain_KK->SetBranchAddress("h1_ProbNNk",&h1_ProbNNk);

  Chain_KK->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS); 
  Chain_KK->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS); 
  Chain_KK->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  Chain_KK->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  Chain_KK->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);

  Chain_KK->SetBranchAddress("D_Hlt2CharmSemilepD02KKMuMuDecision_TOS",&Hlt2_TOS);

  Chain_KK->SetBranchAddress("totCandidates",&totCandidates);
  Chain_KK->SetBranchAddress("nCandidate",&nCandidate);
  Chain_KK->SetBranchAddress("eventNumber",&eventNumber);
  Chain_KK->SetBranchAddress("runNumber",&runNumber);
  Chain_KK->SetBranchAddress("Dst_DTF_D0_M",&m_D);
  Chain_KK->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  Chain_KK->SetBranchAddress("deltaM",&deltaM);

  Chain_KK->SetBranchAddress("BDT",&BDT);
  Chain_KK->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  Chain_KK->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);

  Chain_KK->SetBranchAddress("mu0_PX",&mu0_PX);
  Chain_KK->SetBranchAddress("mu0_PY",&mu0_PY);
  Chain_KK->SetBranchAddress("mu0_PZ",&mu0_PZ);
  Chain_KK->SetBranchAddress("mu1_PX",&mu1_PX);
  Chain_KK->SetBranchAddress("mu1_PY",&mu1_PY);
  Chain_KK->SetBranchAddress("mu1_PZ",&mu1_PZ);

  Chain_KK->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
  Chain_KK->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

  Chain_KK->SetBranchAddress("h0_PX",&h0_PX);
  Chain_KK->SetBranchAddress("h0_PY",&h0_PY);
  Chain_KK->SetBranchAddress("h0_PZ",&h0_PZ);
  Chain_KK->SetBranchAddress("h1_PX",&h1_PX);
  Chain_KK->SetBranchAddress("h1_PY",&h1_PY);
  Chain_KK->SetBranchAddress("h1_PZ",&h1_PZ);

  Chain_Kpi->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);
  Chain_Kpi->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
  Chain_Kpi->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
  Chain_Kpi->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
  Chain_Kpi->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
  Chain_Kpi->SetBranchAddress("h0_ProbNNk",&h0_ProbNNk);
  Chain_Kpi->SetBranchAddress("h1_ProbNNpi",&h1_ProbNNpi);

  Chain_Kpi->SetBranchAddress("nCandidate",&nCandidate);
  Chain_Kpi->SetBranchAddress("eventNumber",&eventNumber);
  Chain_Kpi->SetBranchAddress("runNumber",&runNumber);
  Chain_Kpi->SetBranchAddress("Dst_DTF_D0_M",&m_D);
  Chain_Kpi->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  Chain_Kpi->SetBranchAddress("totCandidates",&totCandidates); 
  Chain_Kpi->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS",&Hlt2_TOS);

  Chain_Kpi->SetBranchAddress("deltaM",&deltaM);

  Chain_Kpi->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Chain_Kpi->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Chain_Kpi->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  Chain_Kpi->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  Chain_Kpi->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);

  Chain_Kpi->SetBranchAddress("BDT",&BDT);
  Chain_Kpi->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  Chain_Kpi->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);

  Chain_Kpi->SetBranchAddress("mu0_PX",&mu0_PX);
  Chain_Kpi->SetBranchAddress("mu0_PY",&mu0_PY);
  Chain_Kpi->SetBranchAddress("mu0_PZ",&mu0_PZ);
  Chain_Kpi->SetBranchAddress("mu1_PX",&mu1_PX);
  Chain_Kpi->SetBranchAddress("mu1_PY",&mu1_PY);
  Chain_Kpi->SetBranchAddress("mu1_PZ",&mu1_PZ);

  Chain_Kpi->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
  Chain_Kpi->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

  Chain_Kpi->SetBranchAddress("h0_PX",&h0_PX);
  Chain_Kpi->SetBranchAddress("h0_PY",&h0_PY);
  Chain_Kpi->SetBranchAddress("h0_PZ",&h0_PZ);
  Chain_Kpi->SetBranchAddress("h1_PX",&h1_PX);
  Chain_Kpi->SetBranchAddress("h1_PY",&h1_PY);
  Chain_Kpi->SetBranchAddress("h1_PZ",&h1_PZ);

  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("h0_ProbNNk",&h0_ProbNNk);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("h1_ProbNNpi",&h1_ProbNNpi);

  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("nCandidate",&nCandidate);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("eventNumber",&eventNumber);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("runNumber",&runNumber);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("Dst_DTF_D0_M",&m_D);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("totCandidates",&totCandidates); 
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS",&Hlt2_TOS);

  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("deltaM",&deltaM);

  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);

  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("BDT",&BDT);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);

  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu0_PX",&mu0_PX);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu0_PY",&mu0_PY);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu0_PZ",&mu0_PZ);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu1_PX",&mu1_PX);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu1_PY",&mu1_PY);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu1_PZ",&mu1_PZ);

  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("h0_PX",&h0_PX);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("h0_PY",&h0_PY);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("h0_PZ",&h0_PZ);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("h1_PX",&h1_PX);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("h1_PY",&h1_PY);
  Chain_Kpi_D2KKmumuBDT->SetBranchAddress("h1_PZ",&h1_PZ);


  Chain_pipi->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);
  Chain_pipi->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
  Chain_pipi->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
  Chain_pipi->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
  Chain_pipi->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
  Chain_pipi->SetBranchAddress("h0_ProbNNpi",&h0_ProbNNpi);
  Chain_pipi->SetBranchAddress("h1_ProbNNpi",&h1_ProbNNpi);

  Chain_pipi->SetBranchAddress("nCandidate",&nCandidate);
  Chain_pipi->SetBranchAddress("eventNumber",&eventNumber);
  Chain_pipi->SetBranchAddress("runNumber",&runNumber);
  Chain_pipi->SetBranchAddress("Dst_DTF_D0_M",&m_D);
  Chain_pipi->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  Chain_pipi->SetBranchAddress("totCandidates",&totCandidates);
  Chain_pipi->SetBranchAddress("deltaM",&deltaM);

  Chain_pipi->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Chain_pipi->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Chain_pipi->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  Chain_pipi->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  Chain_pipi->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);

  Chain_pipi->SetBranchAddress("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS",&Hlt2_TOS);

  Chain_pipi->SetBranchAddress("BDT",&BDT);
  Chain_pipi->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  Chain_pipi->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);

  Chain_pipi->SetBranchAddress("mu0_PX",&mu0_PX);
  Chain_pipi->SetBranchAddress("mu0_PY",&mu0_PY);
  Chain_pipi->SetBranchAddress("mu0_PZ",&mu0_PZ);
  Chain_pipi->SetBranchAddress("mu1_PX",&mu1_PX);
  Chain_pipi->SetBranchAddress("mu1_PY",&mu1_PY);
  Chain_pipi->SetBranchAddress("mu1_PZ",&mu1_PZ);

  Chain_pipi->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
  Chain_pipi->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

  Chain_pipi->SetBranchAddress("h0_PX",&h0_PX);
  Chain_pipi->SetBranchAddress("h0_PY",&h0_PY);
  Chain_pipi->SetBranchAddress("h0_PZ",&h0_PZ);
  Chain_pipi->SetBranchAddress("h1_PX",&h1_PX);
  Chain_pipi->SetBranchAddress("h1_PY",&h1_PY);
  Chain_pipi->SetBranchAddress("h1_PZ",&h1_PZ);

  //double Dst_DTF_D0_M,Dst_DTF_Dstarplus_M;
  //Chain_pipi->SetBranchAddress("Dst_DTF_D0_M",&Dst_DTF_D0_M);
  //Chain_pipi->SetBranchAddress("Dst_DTF_Dstarplus_M",&Dst_DTF_Dstarplus_M);

  
  //Loop                                                                                                                                                                              
  Long64_t nentriesKK = Chain_KK->GetEntries();
  Long64_t nentriesKpi = Chain_Kpi->GetEntries();
  Long64_t nentriesKpi_D2KKmumuBDT = Chain_Kpi_D2KKmumuBDT->GetEntries();
  Long64_t nentriespipi = Chain_pipi->GetEntries();

  int counter1=0;
  int counter2=0;

  for (Long64_t jentry=0; jentry<nentriesKK;jentry++) {

      Chain_KK->GetEntry(jentry);
      
      if(deltaM>146.5 || deltaM < 144.5) continue;
      if(mu0_MuonNShared!=0 || mu1_MuonNShared!=0) continue;
      
      if(withGhosts && !(Dst_BKGCAT<11||Dst_BKGCAT==60))continue;
      //if(withGhosts && !(Dst_BKGCAT<100))continue;
      if(!withGhosts && !(Dst_BKGCAT<11))continue;        
  
      if(BDT<cutBDT) continue;
      if(mu0_ProbNNmu<cutPID||mu1_ProbNNmu<cutPID) continue;
      if(!(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)) continue;
      if(!(mu0_Hlt1TrackMuonDecision_TOS || mu1_Hlt1TrackMuonDecision_TOS|| D_Hlt1TrackAllL0Decision_TOS))continue;
      if(!Hlt2_TOS) continue;
      if(mu0_ProbNNghost>0.5 || mu1_ProbNNghost>0.5 || h0_ProbNNghost>0.5 || h1_ProbNNghost>0.5 || Slowpi_ProbNNghost>0.5) continue;
      if(h0_ProbNNk<0.2 || h1_ProbNNk<0.2 ) continue;
      

      if(totCandidates==1) {counter1++; continue; }

      counter2++;
      
      fileKK << eventNumber << "\t" << runNumber << "\t" << m_D << "\t" << nCandidate << "\t" << Dst_BKGCAT << std::endl;
      /*
      fileKK << eventNumber << "\t" << runNumber << "\t" << m_D  << "\t" << nCandidate << "\t" << Dst_BKGCAT 
	     << "\t" << mu0_PX << "\t" << mu0_PY
	     << "\t" <<mu0_PZ<< "\t" << mu1_PX<< "\t" << mu1_PY<< "\t" << mu1_PZ
      	     << "\t" << h0_PX << "\t" << h0_PY
	     << "\t" << h0_PZ<< "\t" << h1_PX<< "\t" << h1_PY<< "\t" << h1_PZ<< "\t" <<  Slowpi_PX<< "\t" << Slowpi_PY<< "\t" << Slowpi_PZ<< std::endl;
      */
           }
  std::cout<<"KK: single candidate events "<<counter1<<" candidates in multiple cand. "<<counter2<<" sum "<< counter1+counter2<<std::endl; 

  counter1=0;
  counter2=0;

  for (Long64_t jentry=0; jentry<nentriesKpi;jentry++) {

      Chain_Kpi->GetEntry(jentry);
      
      
      if(deltaM>146.5 || deltaM < 144.5) continue;
      if(mu0_MuonNShared!=0 || mu1_MuonNShared!=0) continue;
      
      if(withGhosts && !(Dst_BKGCAT<11||Dst_BKGCAT==60))continue;
      //if(withGhosts && !(Dst_BKGCAT<100))continue;
      if(!withGhosts && !(Dst_BKGCAT<11))continue;        

      if(BDT<cutBDT) continue;
      if(mu0_ProbNNmu<cutPID||mu1_ProbNNmu<cutPID) continue;
      if(!(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)) continue;
      if(!(mu0_Hlt1TrackMuonDecision_TOS || mu1_Hlt1TrackMuonDecision_TOS|| D_Hlt1TrackAllL0Decision_TOS))continue;
      if(!Hlt2_TOS) continue;
      if(mu0_ProbNNghost>0.5 || mu1_ProbNNghost>0.5 || h0_ProbNNghost>0.5 || h1_ProbNNghost>0.5 || Slowpi_ProbNNghost>0.5) continue;
      if(h0_ProbNNk<0.2 || h1_ProbNNpi<0.2 ) continue;
      
      if(totCandidates==1) {counter1++; continue; }

      counter2++;

      fileKpi << eventNumber << "\t" << runNumber << "\t" << m_D << "\t" << nCandidate << "\t" << Dst_BKGCAT  <<std::endl;

      /*
      fileKpi << eventNumber << "\t" << runNumber << "\t" << m_D  << "\t" << nCandidate << "\t" << Dst_BKGCAT
	     << "\t" << mu0_PX << "\t" << mu0_PY
	     << "\t" <<mu0_PZ<< "\t" << mu1_PX<< "\t" << mu1_PY<< "\t" << mu1_PZ
      	     << "\t" << h0_PX << "\t" << h0_PY
           << "\t" << h0_PZ<< "\t" << h1_PX<< "\t" << h1_PY<< "\t" << h1_PZ<< "\t" <<  Slowpi_PX<< "\t" << Slowpi_PY<< "\t" << Slowpi_PZ<< std::endl;
      */

    }

  std::cout<<"Kpi: single candidate events "<<counter1<<" candidates in multiple cand. "<<counter2<<" sum "<< counter1+counter2<<std::endl;

  counter1=0;
  counter2=0;


  for (Long64_t jentry=0; jentry<nentriesKpi_D2KKmumuBDT;jentry++) {

      Chain_Kpi_D2KKmumuBDT->GetEntry(jentry);
      
      
      if(deltaM>146.5 || deltaM < 144.5) continue;
      if(mu0_MuonNShared!=0 || mu1_MuonNShared!=0) continue;
  
      if(withGhosts && !(Dst_BKGCAT<11||Dst_BKGCAT==60))continue;
      //if(withGhosts && !(Dst_BKGCAT<100))continue;
      if(!withGhosts && !(Dst_BKGCAT<11))continue;        

      if(BDT<cutBDT) continue;
      if(mu0_ProbNNmu<cutPID||mu1_ProbNNmu<cutPID) continue;
      if(!(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)) continue;
      if(!(mu0_Hlt1TrackMuonDecision_TOS || mu1_Hlt1TrackMuonDecision_TOS|| D_Hlt1TrackAllL0Decision_TOS))continue;
      if(!Hlt2_TOS) continue;
      if(mu0_ProbNNghost>0.5 || mu1_ProbNNghost>0.5 || h0_ProbNNghost>0.5 || h1_ProbNNghost>0.5 || Slowpi_ProbNNghost>0.5) continue;
      if(h0_ProbNNk<0.2 || h1_ProbNNpi<0.2 ) continue;
      
      if(totCandidates==1) {counter1++; continue; }

      counter2++;

      fileKpi_D2KKmumuBDT << eventNumber << "\t" << runNumber << "\t" << m_D << "\t" << nCandidate << "\t" << Dst_BKGCAT  <<std::endl;

      /*
      fileKpi << eventNumber << "\t" << runNumber << "\t" << m_D  << "\t" << nCandidate << "\t" << Dst_BKGCAT
	     << "\t" << mu0_PX << "\t" << mu0_PY
	     << "\t" <<mu0_PZ<< "\t" << mu1_PX<< "\t" << mu1_PY<< "\t" << mu1_PZ
      	     << "\t" << h0_PX << "\t" << h0_PY
           << "\t" << h0_PZ<< "\t" << h1_PX<< "\t" << h1_PY<< "\t" << h1_PZ<< "\t" <<  Slowpi_PX<< "\t" << Slowpi_PY<< "\t" << Slowpi_PZ<< std::endl;
      */

    }

  std::cout<<"Kpi: single candidate events "<<counter1<<" candidates in multiple cand. "<<counter2<<" sum "<< counter1+counter2<<std::endl;



  counter1=0;
  counter2=0;

  for (Long64_t jentry=0; jentry<nentriespipi;jentry++) {

      Chain_pipi->GetEntry(jentry);
  
      if(deltaM>146.5 || deltaM < 144.5) continue;
      //if(Dst_DTF_Dstarplus_M-Dst_DTF_D0_M < 144.5) continue;
      //if(Dst_DTF_Dstarplus_M-Dst_DTF_D0_M > 146.5) continue;

      if(mu0_MuonNShared!=0 || mu1_MuonNShared!=0) continue;
      
      if(withGhosts && !(Dst_BKGCAT<11||Dst_BKGCAT==60))continue;
      //if(withGhosts && !(Dst_BKGCAT<100))continue;
      if(!withGhosts && !(Dst_BKGCAT<11))continue;        

      if(!(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)) continue;
      if(!(mu0_Hlt1TrackMuonDecision_TOS || mu1_Hlt1TrackMuonDecision_TOS|| D_Hlt1TrackAllL0Decision_TOS))continue;
      if(BDT<cutBDT) continue;
      if(mu0_ProbNNmu<cutPID||mu1_ProbNNmu<cutPID) continue;
      if(!Hlt2_TOS) continue;
      if(mu0_ProbNNghost>0.5 || mu1_ProbNNghost>0.5 || h0_ProbNNghost>0.5 || h1_ProbNNghost>0.5 || Slowpi_ProbNNghost>0.5) continue;
      if(h0_ProbNNpi<0.2 || h1_ProbNNk<0.2 ) continue;
      

      if(totCandidates==1) {counter1++; continue; }

      counter2++;

      filepipi << eventNumber << "\t" << runNumber << "\t" << m_D << "\t" << nCandidate << "\t" << Dst_BKGCAT << std::endl;

      /*
      filepipi << eventNumber << "\t" << runNumber << "\t" << m_D  << "\t" << nCandidate << "\t" << Dst_BKGCAT
	     << "\t" << mu0_PX << "\t" << mu0_PY
	     << "\t" <<mu0_PZ<< "\t" << mu1_PX<< "\t" << mu1_PY<< "\t" << mu1_PZ
      	     << "\t" << h0_PX << "\t" << h0_PY
           << "\t" << h0_PZ<< "\t" << h1_PX<< "\t" << h1_PY<< "\t" << h1_PZ<< "\t" <<  Slowpi_PX<< "\t" << Slowpi_PY<< "\t" << Slowpi_PZ<< std::endl;
      */
  }
  std::cout<<"pi pi: single candidate events "<<counter1<<" candidates in multiple cand. "<<counter2<<" sum "<< counter1+counter2<<std::endl;

    fileKK.close();
    fileKpi.close();
    filepipi.close();

  }
 


bool MatchToMultipleCandidate(vector<Cand> chosen_cand_list, ULong64_t eventNumber, UInt_t runNumber,UInt_t nCandidate){

  for (vector<Cand>::iterator it = chosen_cand_list.begin(); it < chosen_cand_list.end(); ++it) {
    if (eventNumber != (*it).eventNumber) continue;
    if (runNumber != (*it).runNumber) continue;
    if (nCandidate != (*it).nCandidate) continue;
    chosen_cand_list.erase(it);
    return true;
  }
  return false;
}



bool matchToTrueCandidate(vector<TrueCand> trueCand_list, TrueCand recoCand){

    //set the matching criteria 
    double reso_muPX = 0.15;
    double reso_muPY = 0.15;
    double reso_muPZ = 0.05;

    double reso_pisPX = 0.2;
    double reso_pisPY = 0.2;
    double reso_pisPZ = 0.05;

    double reso_hPX = 0.2;
    double reso_hPY = 0.2;
    double reso_hPZ = 0.05;

    int counter_mu0;
    int counter_mu1;
    int counter_h1;
    int counter_h0;
    int counter_slowpi;
    

    for (vector<TrueCand>::iterator it = trueCand_list.begin(); it < trueCand_list.end(); ++it) {

      counter_mu0=0;
      counter_mu1=0;
      counter_h1=0;
      counter_h0=0;
      counter_slowpi=0;

      if (recoCand.eventNumber != (*it).eventNumber) continue;
      if (recoCand.runNumber != (*it).runNumber) continue;
      
      /*
      cout<<recoCand.eventNumber<<"  "<<recoCand.runNumber<<endl;

      cout<<"mu0_PX "<<(recoCand.mu0_PX-(*it).mu0_PX)/(*it).mu0_PX<<"  "<<
	"mu0_PY "<<(recoCand.mu0_PY-(*it).mu0_PY)/(*it).mu0_PY<<"  "<<
	"mu0_PZ "<<(recoCand.mu0_PZ-(*it).mu0_PZ)/(*it).mu0_PZ<<"  "<<       
	"mu1_PX "<<(recoCand.mu1_PX-(*it).mu1_PX)/(*it).mu1_PX<<"  "<<
	"mu1_PY "<<(recoCand.mu1_PY-(*it).mu1_PY)/(*it).mu1_PY<<"  "<<
	"mu1_PZ "<<(recoCand.mu1_PZ-(*it).mu1_PZ)/(*it).mu1_PZ<<"  "<<endl;

      cout<<"h0_PX "<<(recoCand.h0_PX-(*it).h0_PX)/(*it).h0_PX<<"  "<<
	"h0_PY "<<(recoCand.h0_PY-(*it).h0_PY)/(*it).h0_PY<<"  "<<
	"h0_PZ "<<(recoCand.h0_PZ-(*it).h0_PZ)/(*it).h0_PZ<<"  "<<       
	"h1_PX "<<(recoCand.h1_PX-(*it).h1_PX)/(*it).h1_PX<<"  "<<
	"h1_PY "<<(recoCand.h1_PY-(*it).h1_PY)/(*it).h1_PY<<"  "<<
	"h1_PZ "<<(recoCand.h1_PZ-(*it).h1_PZ)/(*it).h1_PZ<<"  "<<endl;

      cout<<"Slowpi_PX "<<(recoCand.Slowpi_PX-(*it).Slowpi_PX)/(*it).Slowpi_PX<<"  "<<
	"Slowpi_PY "<<(recoCand.Slowpi_PY-(*it).Slowpi_PY)/(*it).Slowpi_PY<<"  "<<
	"Slowpi_PZ "<<(recoCand.Slowpi_PZ-(*it).Slowpi_PZ)/(*it).Slowpi_PZ<<"  "<<       
	endl;
      */

      if ( TMath::Abs( (recoCand.mu0_PX-(*it).mu0_PX)/(*it).mu0_PX )<reso_muPX) counter_mu0+=1;
      if ( TMath::Abs( (recoCand.mu0_PY-(*it).mu0_PY)/(*it).mu0_PY )<reso_muPY) counter_mu0+=1;
      if ( TMath::Abs( (recoCand.mu0_PZ-(*it).mu0_PZ)/(*it).mu0_PZ )<reso_muPZ) counter_mu0+=1;

      if ( TMath::Abs( (recoCand.mu1_PX-(*it).mu1_PX)/(*it).mu1_PX )<reso_muPX) counter_mu1+=1;
      if ( TMath::Abs( (recoCand.mu1_PY-(*it).mu1_PY)/(*it).mu1_PY )<reso_muPY) counter_mu1+=1;
      if ( TMath::Abs( (recoCand.mu1_PZ-(*it).mu1_PZ)/(*it).mu1_PZ )<reso_muPZ) counter_mu1+=1;

      if ( TMath::Abs( (recoCand.h0_PX-(*it).h0_PX)/(*it).h0_PX )<reso_hPX) counter_h0+=1;
      if ( TMath::Abs( (recoCand.h0_PY-(*it).h0_PY)/(*it).h0_PY )<reso_hPY) counter_h0+=1;
      if ( TMath::Abs( (recoCand.h0_PZ-(*it).h0_PZ)/(*it).h0_PZ )<reso_hPZ) counter_h0+=1;

      if ( TMath::Abs( (recoCand.h1_PX-(*it).h1_PX)/(*it).h1_PX )<reso_hPX) counter_h1+=1;
      if ( TMath::Abs( (recoCand.h1_PY-(*it).h1_PY)/(*it).h1_PY )<reso_hPY) counter_h1+=1;
      if ( TMath::Abs( (recoCand.h1_PZ-(*it).h1_PZ)/(*it).h1_PZ )<reso_hPZ) counter_h1+=1;

      if ( TMath::Abs( (recoCand.Slowpi_PX-(*it).Slowpi_PX)/(*it).Slowpi_PX )<reso_pisPX) counter_slowpi+=1;
      if ( TMath::Abs( (recoCand.Slowpi_PY-(*it).Slowpi_PY)/(*it).Slowpi_PY )<reso_pisPY) counter_slowpi+=1;
      if ( TMath::Abs( (recoCand.Slowpi_PZ-(*it).Slowpi_PZ)/(*it).Slowpi_PZ )<reso_pisPZ) counter_slowpi+=1;

      //cout<<recoCand.eventNumber<<"  "<<counter_mu0<<"  "<<counter_mu1<<"  "<<counter_h0<<"  "<<counter_h1<<"  "<<counter_slowpi<<endl;
      //2 momentum components have to match
      if(counter_mu0<2) continue;
      if(counter_mu1<2) continue;
      if(counter_h0<2) continue;
      if(counter_h1<2) continue;
      if(counter_slowpi<2) continue;
      //cout<<"match"<<endl;
      /*
      if ( TMath::Abs( (recoCand.mu0_PX-(*it).mu0_PX)/(*it).mu0_PX )>reso_muPX) continue;
      if ( TMath::Abs( (recoCand.mu0_PY-(*it).mu0_PY)/(*it).mu0_PY )>reso_muPY) continue;
      if ( TMath::Abs( (recoCand.mu0_PZ-(*it).mu0_PZ)/(*it).mu0_PZ )>reso_muPZ) continue;

      if ( TMath::Abs( (recoCand.mu1_PX-(*it).mu1_PX)/(*it).mu1_PX )>reso_muPX) continue;
      if ( TMath::Abs( (recoCand.mu1_PY-(*it).mu1_PY)/(*it).mu1_PY )>reso_muPY) continue;
      if ( TMath::Abs( (recoCand.mu1_PZ-(*it).mu1_PZ)/(*it).mu1_PZ )>reso_muPZ) continue;

      if ( TMath::Abs( (recoCand.h0_PX-(*it).h0_PX)/(*it).h0_PX )>reso_hPX) continue;
      if ( TMath::Abs( (recoCand.h0_PY-(*it).h0_PY)/(*it).h0_PY )>reso_hPY) continue;
      if ( TMath::Abs( (recoCand.h0_PZ-(*it).h0_PZ)/(*it).h0_PZ )>reso_hPZ) continue;

      if ( TMath::Abs( (recoCand.h1_PX-(*it).h1_PX)/(*it).h1_PX )>reso_hPX) continue;
      if ( TMath::Abs( (recoCand.h1_PY-(*it).h1_PY)/(*it).h1_PY )>reso_hPY) continue;
      if ( TMath::Abs( (recoCand.h1_PZ-(*it).h1_PZ)/(*it).h1_PZ )>reso_hPZ) continue;

      if ( TMath::Abs( (recoCand.Slowpi_PX-(*it).Slowpi_PX)/(*it).Slowpi_PX )>reso_pisPX) continue;
      if ( TMath::Abs( (recoCand.Slowpi_PY-(*it).Slowpi_PY)/(*it).Slowpi_PY )>reso_pisPY) continue;
      if ( TMath::Abs( (recoCand.Slowpi_PZ-(*it).Slowpi_PZ)/(*it).Slowpi_PZ )>reso_pisPZ) continue;
      */
      trueCand_list.erase(it);

      return true;
    }
    return false;
  }
			  
			  

 
void writeAllTrueCanidatesToFile(){
  
  //write all the multiple canidated to a file with totCandidates>1 (after dm preselction and nShared cut for muons)

  std::ofstream fileKpi("AllTrue_candidatesKpi_2012.txt");
  std::ofstream fileKK("AllTruel_candidatesKK_2012.txt");
  std::ofstream filepipi("AllTrue_candidatespipi_2012.txt");


  TChain* Chain_KK = new TChain("MCTruthTuple");
  TChain* Chain_Kpi=  new TChain("MCTruthTuple");
  TChain* Chain_pipi=  new TChain("MCTruthTuple");

  Chain_KK->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/MCTruthTuple_DstD2KKmumu.root");
  Chain_Kpi->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/MCTruthTuple_DstD2Kpimumu.root");
  Chain_pipi->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MCTruthTuple_DstD2pipimumu.root");
  
  double m_D;
  ULong64_t eventNumber;
  UInt_t runNumber;
  UInt_t nCandidate;
  int Dst_BKGCAT;
  ULong64_t totCandidates;

  double mu0_TRUEP_X, mu0_TRUEP_Y, mu0_TRUEP_Z;
  double mu1_TRUEP_X, mu1_TRUEP_Y, mu1_TRUEP_Z;
  double h0_TRUEP_X, h0_TRUEP_Y, h0_TRUEP_Z;
  double h1_TRUEP_X, h1_TRUEP_Y, h1_TRUEP_Z;
  double Slowpi_TRUEP_X,Slowpi_TRUEP_Y,Slowpi_TRUEP_Z;
  double deltaM; 

  Chain_KK->SetBranchAddress("eventNumber",&eventNumber);
  Chain_KK->SetBranchAddress("runNumber",&runNumber);

  Chain_KK->SetBranchAddress("muplus_TRUEP_X",&mu0_TRUEP_X);
  Chain_KK->SetBranchAddress("muplus_TRUEP_Y",&mu0_TRUEP_Y);
  Chain_KK->SetBranchAddress("muplus_TRUEP_Z",&mu0_TRUEP_Z);
  Chain_KK->SetBranchAddress("muminus_TRUEP_X",&mu1_TRUEP_X);
  Chain_KK->SetBranchAddress("muminus_TRUEP_Y",&mu1_TRUEP_Y);
  Chain_KK->SetBranchAddress("muminus_TRUEP_Z",&mu1_TRUEP_Z);

  Chain_KK->SetBranchAddress("Kminus_TRUEP_X",&h0_TRUEP_X);
  Chain_KK->SetBranchAddress("Kminus_TRUEP_Y",&h0_TRUEP_Y);
  Chain_KK->SetBranchAddress("Kminus_TRUEP_Z",&h0_TRUEP_Z);
  Chain_KK->SetBranchAddress("Kplus_TRUEP_X",&h1_TRUEP_X);
  Chain_KK->SetBranchAddress("Kplus_TRUEP_Y",&h1_TRUEP_Y);
  Chain_KK->SetBranchAddress("Kplus_TRUEP_Z",&h1_TRUEP_Z);

  Chain_KK->SetBranchAddress("piplus_TRUEP_X",&Slowpi_TRUEP_X);
  Chain_KK->SetBranchAddress("piplus_TRUEP_Y",&Slowpi_TRUEP_Y);
  Chain_KK->SetBranchAddress("piplus_TRUEP_Z",&Slowpi_TRUEP_Z);

  Chain_Kpi->SetBranchAddress("eventNumber",&eventNumber);
  Chain_Kpi->SetBranchAddress("runNumber",&runNumber);

  Chain_Kpi->SetBranchAddress("muplus_TRUEP_X",&mu0_TRUEP_X);
  Chain_Kpi->SetBranchAddress("muplus_TRUEP_Y",&mu0_TRUEP_Y);
  Chain_Kpi->SetBranchAddress("muplus_TRUEP_Z",&mu0_TRUEP_Z);
  Chain_Kpi->SetBranchAddress("muminus_TRUEP_X",&mu1_TRUEP_X);
  Chain_Kpi->SetBranchAddress("muminus_TRUEP_Y",&mu1_TRUEP_Y);
  Chain_Kpi->SetBranchAddress("muminus_TRUEP_Z",&mu1_TRUEP_Z);

  Chain_Kpi->SetBranchAddress("Kminus_TRUEP_X",&h0_TRUEP_X);
  Chain_Kpi->SetBranchAddress("Kminus_TRUEP_Y",&h0_TRUEP_Y);
  Chain_Kpi->SetBranchAddress("Kminus_TRUEP_Z",&h0_TRUEP_Z);
  Chain_Kpi->SetBranchAddress("piplus_TRUEP_X",&h1_TRUEP_X);
  Chain_Kpi->SetBranchAddress("piplus_TRUEP_Y",&h1_TRUEP_Y);
  Chain_Kpi->SetBranchAddress("piplus_TRUEP_Z",&h1_TRUEP_Z);

  Chain_Kpi->SetBranchAddress("piplus0_TRUEP_X",&Slowpi_TRUEP_X);
  Chain_Kpi->SetBranchAddress("piplus0_TRUEP_Y",&Slowpi_TRUEP_Y);
  Chain_Kpi->SetBranchAddress("piplus0_TRUEP_Z",&Slowpi_TRUEP_Z);

  Chain_pipi->SetBranchAddress("eventNumber",&eventNumber);
  Chain_pipi->SetBranchAddress("runNumber",&runNumber);

  Chain_pipi->SetBranchAddress("muplus_TRUEP_X",&mu0_TRUEP_X);
  Chain_pipi->SetBranchAddress("muplus_TRUEP_Y",&mu0_TRUEP_Y);
  Chain_pipi->SetBranchAddress("muplus_TRUEP_Z",&mu0_TRUEP_Z);
  Chain_pipi->SetBranchAddress("muminus_TRUEP_X",&mu1_TRUEP_X);
  Chain_pipi->SetBranchAddress("muminus_TRUEP_Y",&mu1_TRUEP_Y);
  Chain_pipi->SetBranchAddress("muminus_TRUEP_Z",&mu1_TRUEP_Z);

  Chain_pipi->SetBranchAddress("piminus_TRUEP_X",&h0_TRUEP_X);
  Chain_pipi->SetBranchAddress("piminus_TRUEP_Y",&h0_TRUEP_Y);
  Chain_pipi->SetBranchAddress("piminus_TRUEP_Z",&h0_TRUEP_Z);
  Chain_pipi->SetBranchAddress("piplus_TRUEP_X",&h1_TRUEP_X);
  Chain_pipi->SetBranchAddress("piplus_TRUEP_Y",&h1_TRUEP_Y);
  Chain_pipi->SetBranchAddress("piplus_TRUEP_Z",&h1_TRUEP_Z);
 
  Chain_pipi->SetBranchAddress("piplus0_TRUEP_X",&Slowpi_TRUEP_X);
  Chain_pipi->SetBranchAddress("piplus0_TRUEP_Y",&Slowpi_TRUEP_Y);
  Chain_pipi->SetBranchAddress("piplus0_TRUEP_Z",&Slowpi_TRUEP_Z);

  
  
  //Loop                                                                                                                                                                              
  Long64_t nentriesKK = Chain_KK->GetEntries();
  Long64_t nentriesKpi = Chain_Kpi->GetEntries();
  Long64_t nentriespipi = Chain_pipi->GetEntries();

  for (Long64_t jentry=0; jentry<nentriesKK;jentry++) {

      Chain_KK->GetEntry(jentry);
        
      fileKK << eventNumber << "\t" << runNumber 
	     << "\t" << mu0_TRUEP_X << "\t" << mu0_TRUEP_Y <<"\t" <<mu0_TRUEP_Z
	     << "\t" << mu1_TRUEP_X<< "\t" << mu1_TRUEP_Y<< "\t" << mu1_TRUEP_Z
      	     << "\t" << h0_TRUEP_X << "\t" << h0_TRUEP_Y << "\t" << h0_TRUEP_Z
	     << "\t" << h1_TRUEP_X<< "\t" << h1_TRUEP_Y<< "\t" << h1_TRUEP_Z
	     << "\t" <<  Slowpi_TRUEP_X<< "\t" << Slowpi_TRUEP_Y<< "\t" << Slowpi_TRUEP_Z<< std::endl;
  }

  for (Long64_t jentry=0; jentry<nentriesKpi;jentry++) {

      Chain_Kpi->GetEntry(jentry);

      fileKpi << eventNumber << "\t" << runNumber 
	     << "\t" << mu0_TRUEP_X << "\t" << mu0_TRUEP_Y <<"\t" <<mu0_TRUEP_Z
	     << "\t" << mu1_TRUEP_X<< "\t" << mu1_TRUEP_Y<< "\t" << mu1_TRUEP_Z
      	     << "\t" << h0_TRUEP_X << "\t" << h0_TRUEP_Y << "\t" << h0_TRUEP_Z
	     << "\t" << h1_TRUEP_X<< "\t" << h1_TRUEP_Y<< "\t" << h1_TRUEP_Z
	     << "\t" <<  Slowpi_TRUEP_X<< "\t" << Slowpi_TRUEP_Y<< "\t" << Slowpi_TRUEP_Z<< std::endl;

    }


  for (Long64_t jentry=0; jentry<nentriespipi;jentry++) {

      Chain_pipi->GetEntry(jentry);

      filepipi << eventNumber << "\t" << runNumber 
	     << "\t" << mu0_TRUEP_X << "\t" << mu0_TRUEP_Y <<"\t" <<mu0_TRUEP_Z
	     << "\t" << mu1_TRUEP_X<< "\t" << mu1_TRUEP_Y<< "\t" << mu1_TRUEP_Z
      	     << "\t" << h0_TRUEP_X << "\t" << h0_TRUEP_Y << "\t" << h0_TRUEP_Z
	     << "\t" << h1_TRUEP_X<< "\t" << h1_TRUEP_Y<< "\t" << h1_TRUEP_Z
	     << "\t" <<  Slowpi_TRUEP_X<< "\t" << Slowpi_TRUEP_Y<< "\t" << Slowpi_TRUEP_Z<< std::endl;

   }

    fileKK.close();
    fileKpi.close();
    filepipi.close();

  }
 



void chose_multiple_events(
			   const char *input="/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MC12_DstD2pipimumu.root",
			   TString channel="D2pipimumu",
			   bool withGhosts=false
			   ){
  //take a MC file and identify the multiple candidates writte to candidates**_2012.txt and chose the one with lowest Dst_BKGCAT. chosen candidated are written to chosen_candidates**_2012.txt

  double cutBDT=0.4;
  double cutPID=0.5;

  vector<Cand> events;
  vector<Cand> chosen_events;

  Cand tmp;
  TString pathToFiles="candidateLists/";

  TString filename = pathToFiles;
  
  if(channel=="D2pipimumu") filename+="candidatespipi_2012_filtered";
  if(channel=="D2Kpimumu_D2KKmumuBDT") filename+="candidatesKpi_2012_filtered_D2KKmumuBDT";
  if(channel=="D2Kpimumu_D2pipimumuBDT") filename+="candidatesKpi_2012_filtered_D2pipimumuBDT";
  if(channel=="D2KKmumu") filename+="candidatesKK_2012_filtered";
  
  if(withGhosts) filename+="_withGhosts.txt";
  else filename+="_noGhosts.txt";

  TString outputfilename;
  if(channel=="D2pipimumu")outputfilename = pathToFiles+"chosen_candidatespipi_2012_filtered";
  if(channel=="D2Kpimumu_D2KKmumuBDT")outputfilename = pathToFiles+"chosen_candidatesKpi_2012_filtered_D2KKmumuBDT";
  if(channel=="D2Kpimumu_D2pipimumuBDT")outputfilename = pathToFiles+"chosen_candidatesKpi_2012_filtered_D2pipimumuBDT";
  if(channel=="D2KKmumu") outputfilename= pathToFiles+"chosen_candidatesKK_2012_filtered";

  if(withGhosts) outputfilename+="_withGhosts.txt";
  else outputfilename+="_noGhosts.txt";

  /// read files with candidates that had originally more than one candidate in the event                                                                                                                                                                                      /// and create a vector with all the matches candidates from the file                                                                                                                                                                                                                                                                                                                                                                                     
  cout << "Readingfile "<<filename<<" ... to " <<outputfilename<< endl;
  ifstream rs(filename);
  if (!rs) {
    cout << "Unable to open " << filename << endl;
    return;
  }


  while (rs >> tmp.eventNumber >> tmp.runNumber  >>  tmp.mass >> tmp.nCandidate >> tmp.Dst_BKGCAT ) {
    events.push_back(tmp);
  }
  rs.close();
  cout << "Done: " << events.size() << " events found." << endl;



  ///___________________________________________________________________________________                                                                                                                                                                                                                                                                                                                                                                                                           
  TString inputfile_data;
  TChain* chain = new TChain("BDT_Tree");
  //TChain* chain = new TChain("MC12_DstD2pipiMuMu/DecayTree");                                                                                                                                      

  chain->AddFile(input);

  //____________________________________________________________________________________                                                   

  double mD;
  ULong64_t loop_eventNumber;
  UInt_t loop_runNumber;
  UInt_t nCandidate;
  int Dst_BKGCAT;
  ULong64_t totCandidates;
 
  Cand loopCand;
  vector<Cand> chosenCandidates;
  vector<Cand> multCandidates;

  double h0_ProbNNghost, h1_ProbNNghost,mu0_ProbNNghost, mu1_ProbNNghost, Slowpi_ProbNNghost;
  double h0_ProbNN,h1_ProbNN;
  double BDT, mu0_ProbNNmu, mu1_ProbNNmu;
  bool mu1_L0MuonDecision_TOS,mu0_L0MuonDecision_TOS;
  bool mu1_Hlt1TrackMuonDecision_TOS, mu0_Hlt1TrackMuonDecision_TOS, D_Hlt1TrackAllL0Decision_TOS;
  bool Hlt2_TOS;
  int mu1_MuonNShared,mu0_MuonNShared;
  double deltaM;

  chain->SetBranchAddress("nCandidate",&nCandidate);
  chain->SetBranchAddress("eventNumber",&loop_eventNumber);
  chain->SetBranchAddress("runNumber",&loop_runNumber);
  chain->SetBranchAddress("Dst_DTF_D0_M",&mD);
  chain->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);

  chain->SetBranchAddress("totCandidates",&totCandidates);
  chain->SetBranchAddress("deltaM",&deltaM);
  chain->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);
  chain->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
  chain->SetBranchAddress("BDT",&BDT);
  chain->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  chain->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
  chain->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  chain->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  chain->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  chain->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  chain->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);
  chain->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
  chain->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
  chain->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);
  chain->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
  chain->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);

  if(channel=="D2KKmumu") {
    chain->SetBranchAddress("h0_ProbNNk",&h0_ProbNN);
    chain->SetBranchAddress("h1_ProbNNk",&h1_ProbNN);
    chain->SetBranchAddress("D_Hlt2CharmSemilepD02KKMuMuDecision_TOS",&Hlt2_TOS);
  }

  if(channel=="D2Kpimumu") {
    chain->SetBranchAddress("h0_ProbNNk",&h0_ProbNN);
    chain->SetBranchAddress("h1_ProbNNpi",&h1_ProbNN);
    chain->SetBranchAddress("D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS",&Hlt2_TOS);
  }


  if(channel=="D2pipimumu") {
    chain->SetBranchAddress("h0_ProbNNpi",&h0_ProbNN);
    chain->SetBranchAddress("h1_ProbNNpi",&h1_ProbNN);
    chain->SetBranchAddress("D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS",&Hlt2_TOS);
  }

  //Loop                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
  Long64_t nentries = chain->GetEntries();
  
  int counter1=0;int counter2=0;int counter3=0;


  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    chain->GetEntry(jentry);

    
    if(deltaM>146.5 || deltaM < 144.5) continue;
    if(mu0_MuonNShared!=0 || mu1_MuonNShared!=0) continue;
    if(BDT<cutBDT) continue;
    if(mu0_ProbNNmu<cutPID||mu1_ProbNNmu<cutPID) continue;
    if(!(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)) continue;
    if(!(mu0_Hlt1TrackMuonDecision_TOS || mu1_Hlt1TrackMuonDecision_TOS|| D_Hlt1TrackAllL0Decision_TOS))continue;
    if(!Hlt2_TOS) continue;
    if(mu0_ProbNNghost>0.5 || mu1_ProbNNghost>0.5 || h0_ProbNNghost>0.5 || h1_ProbNNghost>0.5 || Slowpi_ProbNNghost>0.5) continue;
    if(h0_ProbNN<0.2 || h1_ProbNN<0.2 ) continue;

    if(withGhosts && !(Dst_BKGCAT<11||Dst_BKGCAT==60)) continue;
    //if(withGhosts && !(Dst_BKGCAT<100)) continue;
    if(!withGhosts && !(Dst_BKGCAT<11)) continue;
    

    if(totCandidates==1) continue;

    /// create temporary candidate to campare with candidates from file                                                                                                                                                                                                                                                                                                                                  
    loopCand.mass=mD;
    loopCand.eventNumber=loop_eventNumber;
    loopCand.runNumber=loop_runNumber;
    loopCand.nCandidate=nCandidate;
    loopCand.Dst_BKGCAT=Dst_BKGCAT;

   ///every candidate is put into a vector which will be added by the candidates of the same event in the following                                                                                                                                                                                                                                                                                         
    multCandidates.push_back(loopCand);

    for (int i=1;i<events.size();++i) { ///compare to every candidate from file                                                                                                                                                                                                 
      if(events[i].eventNumber==0) continue; ///checks if candidate has already been matched (to avoid double counting)                 

      if(loop_eventNumber==events[i].eventNumber && loop_runNumber==events[i].runNumber && nCandidate!=events[i].nCandidate ) {          ///find partner, but not the cand itself! (nCand)
      //if( (TMath::Abs((double)loop_eventNumber-(double)events[i].eventNumber)<1) && (TMath::Abs((double)loop_runNumber-(double)events[i].runNumber)<1) && (TMath::Abs((double)nCandidate-(double)events[i].nCandidate)>1) ) {  
	//cout<<jentry<<"  "<<i<<"  "<<nCandidate<<"  "<<events[i].nCandidate <<endl;               
	multCandidates.push_back(events[i]); /// if partner is found put into vector                                                    
	events[i].eventNumber=0; ///candidate has already found its partners, so set event number to 0                               
      }
      if(loop_eventNumber==events[i].eventNumber && loop_runNumber==events[i].runNumber && nCandidate==events[i].nCandidate ) {  ///find itself and "delete" from list by setting eventnumber to 0     
	events[i].eventNumber=0;
      }

    }
    //cout<<"event "<<jentry<<endl;
    int nCand=multCandidates.size();
    //cout<<"nCand "<<nCand<<endl;
    counter1+=nCand;
    if(nCand>1) { ///multiple candidate found if in the vector multCand has more than one entry, chose one randomly                                                   
      //cout<<"in Loop"<<endl;
      counter2+=1;
      TRandom3 generator(loop_eventNumber); ///seed is event number                                                                                             
      double randomNr=generator.Rndm();
      int index=int(nCand*randomNr);
      //cout<<"chosen index "<<index<<" out of "<< nCand<<endl;                                                                                                               
      int littleCounter=0;
      for(int k=0;k<multCandidates.size();++k){
	if(k!=index) {chosenCandidates.push_back(multCandidates[index]); ///chosen candidates contains the candidates to be removed from the events in the loop   
	  //std::cout<<"reject "<<k<<std::endl;
	  counter3+=1;
	  littleCounter+=1;}
      }
      //cout<<"removed "<<littleCounter<<"  "<<littleCounter+1-nCand<<endl; 
    }

    multCandidates.clear(); ///clear vector for next event   


    //std::vector<int> categories;
    //double nCand=multCandidates.size();
    //cout<<"candidates "<<multCandidates.size()<<endl;
    //if(nCand>1) { ///multiple candidate found if in the vector multCand has more than one entry, chose one with lowest BDKCAT  
      //cout<<"new event with mult cand: "<<nCand<<endl;
      //for(vector<Cand>::iterator it = multCandidates.begin();it!=multCandidates.end();++it){
      //categories.push_back((*it).Dst_BKGCAT);
	//cout<<(*it).Dst_BKGCAT<<endl;
      //}
      //int randomIndex = rand() % multCandidates.size();
      ///chosen candidates contains the candidates to be removed from the events in the loop                              
      //int index=distance(categories.begin(),min_element(categories.begin(),categories.end()));
      //chosenCandidates.push_back(multCandidates[index]); 
      //chosenCandidates.push_back(multCandidates[randomIndex]);
      //if(nCand>1) cout<<"random index: "<<randomIndex<<endl;
       //}
    //multCandidates.clear(); ///clear vector for next event                                                                                
  }                                                                                                                                     
  ///write the chosen candidates to an output file                                                                              

  cout << "Writing output file..." << endl;
  ofstream fout(outputfilename);
  for (vector<Cand>::iterator it=chosenCandidates.begin(); it<chosenCandidates.end(); ++it) {
    fout << (*it).eventNumber << "\t" << (*it).runNumber << "\t" << (*it).nCandidate << "\t" << (*it).mass<< "\t" << (*it).Dst_BKGCAT<< endl;
  }
  fout.close();

  cout << "Done. Chosen candidates: " << chosenCandidates.size()<<endl;
  //cout << counter1 <<"  "<<counter2<<"  "<<counter3<< endl;
}


void addMultipleCandidateBranch(TChain *chain, TString fOut,const char *f_chosenCand){
   
  //take the flaggeg MC tuples and add a flag for selected multiple candidates

  vector<Cand> matched_list;
  Cand event;
  std::ifstream file(f_chosenCand);

  if ( !file ) {
    std::cout << "Unable to open " << "matched_cands_list_2012.txt" << std::endl;
  } else
    std::cout << "Reading matched candidates from " << f_chosenCand  << "..." << std::endl;


  while ( (file >> event.eventNumber >> event.runNumber >>  event.nCandidate >> event.mass >> event.Dst_BKGCAT )) {
    matched_list.push_back(event);
  }
  file.close();
  std::cout << "Done." << std::endl;

  double mD;
  ULong64_t eventNumber;
  UInt_t runNumber;
  UInt_t nCandidate;
  int Dst_BKGCAT;
  ULong64_t totCandidates;
  
  chain->SetBranchAddress("nCandidate",&nCandidate);
  chain->SetBranchAddress("eventNumber",&eventNumber);
  chain->SetBranchAddress("runNumber",&runNumber);
  chain->SetBranchAddress("Dst_DTF_D0_M",&mD);
  chain->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain->SetBranchAddress("totCandidates",&totCandidates);
    
   TFile* targetFile = new TFile(fOut,"RECREATE");

  TTree* new_tree = chain->CopyTree("");
  bool isRejectedMultipleCandidate;
  TBranch* Bra = new_tree->Branch("isRejectedMultipleCandidate",&isRejectedMultipleCandidate);
 
  int numEvents = chain ->GetEntries();
  
  for(int i=0; i< numEvents; i++){
    
    chain->GetEntry(i);
    if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;

    if( MatchToMultipleCandidate(matched_list, eventNumber, runNumber, nCandidate) ) isRejectedMultipleCandidate=true;  // if is found on list, reject!
    else isRejectedMultipleCandidate = false;
    Bra->Fill();
  } 
  
  new_tree->Write();
  targetFile->Write();
  targetFile->Close();
  //delete new_tree;
  delete targetFile;

}


void addCustomTruthMatchingBranch(TChain *chain, TString fOut,const char *f_chosenCand){
   
  //take the flaggeg MC tuples and add a flag for selected multiple candidates

  vector<TrueCand> matched_list;
  TrueCand event;
  std::ifstream file(f_chosenCand);

  if ( !file ) {
    std::cout << "Unable to open " << "matched_cands_list_2012.txt" << std::endl;
  } else
    std::cout << "Reading matched candidates from " << f_chosenCand  << "..." << std::endl;


  while ( (file >> event.eventNumber >> event.runNumber 
	   >>event.mu0_PX >> event.mu0_PY  >> event.mu0_PZ 
	   >>event.mu1_PX >> event.mu1_PY  >> event.mu1_PZ 
	   >>event.h0_PX >> event.h0_PY  >> event.h0_PZ 
	   >>event.h1_PX >> event.h1_PY  >> event.h1_PZ 
	   >>event.Slowpi_PX >> event.Slowpi_PY  >> event.Slowpi_PZ )){
    matched_list.push_back(event);
  }

  file.close();
  std::cout << "Done." << std::endl;

  double mD;
  ULong64_t eventNumber;
  UInt_t runNumber;
  UInt_t nCandidate;
  int Dst_BKGCAT;
  ULong64_t totCandidates;
  double mu0_PX, mu0_PY, mu0_PZ;
  double mu1_PX, mu1_PY, mu1_PZ;
  double h0_PX, h0_PY, h0_PZ;
  double h1_PX, h1_PY, h1_PZ;
  double Slowpi_PX,Slowpi_PY,Slowpi_PZ;

  TrueCand recoCand; //sound abit weird. Is a candidate with the reconstructed momenta of the candidate
  
  chain->SetBranchAddress("eventNumber",&eventNumber);
  chain->SetBranchAddress("runNumber",&runNumber);

  chain->SetBranchAddress("mu0_PX",&mu0_PX);
  chain->SetBranchAddress("mu0_PY",&mu0_PY);
  chain->SetBranchAddress("mu0_PZ",&mu0_PZ);

  chain->SetBranchAddress("mu1_PX",&mu1_PX);
  chain->SetBranchAddress("mu1_PY",&mu1_PY);
  chain->SetBranchAddress("mu1_PZ",&mu1_PZ);

  chain->SetBranchAddress("h0_PX",&h0_PX);
  chain->SetBranchAddress("h0_PY",&h0_PY);
  chain->SetBranchAddress("h0_PZ",&h0_PZ);

  chain->SetBranchAddress("h1_PX",&h1_PX);
  chain->SetBranchAddress("h1_PY",&h1_PY);
  chain->SetBranchAddress("h1_PZ",&h1_PZ);

  chain->SetBranchAddress("Slowpi_PX",&Slowpi_PX);
  chain->SetBranchAddress("Slowpi_PY",&Slowpi_PY);
  chain->SetBranchAddress("Slowpi_PZ",&Slowpi_PZ);
    
  TFile* targetFile = new TFile(fOut,"RECREATE");

  TTree* new_tree = chain->CopyTree("");
  bool isMatchedCandidate;
  TBranch* Bra = new_tree->Branch("isMatchedCandidate",&isMatchedCandidate);
  int index;
  TBranch* Bra_inex = new_tree->Branch("index",&index);
 
  int numEvents = chain ->GetEntries();
  
  for(int i=0; i< numEvents; i++){
  //for(int i=0; i< 10000; i++){
    
    chain->GetEntry(i);
    if (0ul == (i % 20ul)) cout << "Read event " << i << "/" << numEvents << endl;
    
    //set reco candidate;
    recoCand.eventNumber = eventNumber;
    recoCand.runNumber = runNumber;

    recoCand.mu0_PX = mu0_PX;
    recoCand.mu0_PY = mu0_PY;
    recoCand.mu0_PZ = mu0_PZ;

    recoCand.mu1_PX = mu1_PX;
    recoCand.mu1_PY = mu1_PY;
    recoCand.mu1_PZ = mu1_PZ;

    recoCand.h0_PX = h0_PX;
    recoCand.h0_PY = h0_PY;
    recoCand.h0_PZ = h0_PZ;

    recoCand.h1_PX = h1_PX;
    recoCand.h1_PY = h1_PY;
    recoCand.h1_PZ = h1_PZ;

    recoCand.Slowpi_PX = Slowpi_PX;
    recoCand.Slowpi_PY = Slowpi_PY;
    recoCand.Slowpi_PZ = Slowpi_PZ;

    isMatchedCandidate = matchToTrueCandidate(matched_list, recoCand);
    index=i;
    Bra->Fill();
    Bra_inex->Fill();
    //    recoCand.Clear();
  } 
  
  new_tree->Write();
  targetFile->Write();
  targetFile->Close();
  //delete new_tree;
  delete targetFile;

}



 
void createMCTuplesForEffStudies() {
 

  TChain* Tree_D2Kpimumu_magUp = new TChain("MC12_DstD2KPiMuMu/DecayTree");
  Tree_D2Kpimumu_magUp->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu_withGenLevel/magUp/MC12_DstD2Kpimumu_magUp.root");
  TChain* Tree_D2Kpimumu_magDw = new TChain("MC12_DstD2KPiMuMu/DecayTree");
  Tree_D2Kpimumu_magDw->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu_withGenLevel/magDw/MC12_DstD2Kpimumu_magDw.root");

  TChain* Tree_D2KKmumu_magUp = new TChain("MC12_DstD2KKMuMu/DecayTree");
  Tree_D2KKmumu_magUp->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_withGenLevel/magUp/MC12_DstD2KKmumu_magUp.root");
  TChain* Tree_D2KKmumu_magDw = new TChain("MC12_DstD2KKMuMu/DecayTree");
  Tree_D2KKmumu_magDw->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_withGenLevel/magDw/MC12_DstD2KKmumu_magDw.root");
 
  TChain* Tree_D2pipimumu_magUp = new TChain("MC12_DstD2PiPiMuMu/DecayTree");
  Tree_D2pipimumu_magUp->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu_withGenLevel/magUp/MC12_DstD2pipimumu_magUp.root");
  TChain* Tree_D2pipimumu_magDw = new TChain("MC12_DstD2PiPiMuMu/DecayTree");
  Tree_D2pipimumu_magDw->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu_withGenLevel/magDw/MC12_DstD2pipimumu_magDw.root");
   

  D2KpimumuReader* Kpi_Reader_magUp = new D2KpimumuReader(Tree_D2Kpimumu_magUp);
  D2KpimumuReader* Kpi_Reader_magDw = new D2KpimumuReader(Tree_D2Kpimumu_magDw);
  D2KKmumuReader* KK_Reader_magUp = new D2KKmumuReader(Tree_D2KKmumu_magUp);
  D2KKmumuReader* KK_Reader_magDw = new D2KKmumuReader(Tree_D2KKmumu_magDw);
  D2pipimumuReader* pipi_Reader_magUp = new D2pipimumuReader(Tree_D2pipimumu_magUp);
  D2pipimumuReader* pipi_Reader_magDw = new D2pipimumuReader(Tree_D2pipimumu_magDw);
 
  
  for(int i=0; i<rangesKpi_low.size();++i){
    Kpi_Reader_magUp->createMCEfficiencyStudySample(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),rangesKpi_low[i],rangesKpi_high[i]);
    Kpi_Reader_magDw->createMCEfficiencyStudySample(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),rangesKpi_low[i],rangesKpi_high[i]);
  }

  
  
  for(int i=0; i<rangesKK_low.size();++i){
    KK_Reader_magUp->createMCEfficiencyStudySample(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_magUp.root",rangesKK_low[i],rangesKK_high[i]),rangesKK_low[i],rangesKK_high[i]);
    KK_Reader_magDw->createMCEfficiencyStudySample(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_magDw.root",rangesKK_low[i],rangesKK_high[i]),rangesKK_low[i],rangesKK_high[i]);
  }
  
  for(int i=0; i<rangespipi_low.size();++i){
  
   
  pipi_Reader_magUp->createMCEfficiencyStudySample(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_magUp.root",rangespipi_low[i],rangespipi_high[i]),rangespipi_low[i],rangespipi_high[i]);
    pipi_Reader_magDw->createMCEfficiencyStudySample(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_magDw.root",rangespipi_low[i],rangespipi_high[i]),rangespipi_low[i],rangespipi_high[i]);
  }
 

}

 


void createMCTuplesForEffStudiesNoTruthmatching() {
 

  TChain* Tree_D2Kpimumu_magUp = new TChain("MC12_DstD2KPiMuMu/DecayTree");
  Tree_D2Kpimumu_magUp->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu_withGenLevel/magUp/MC12_DstD2Kpimumu_magUp.root");
  TChain* Tree_D2Kpimumu_magDw = new TChain("MC12_DstD2KPiMuMu/DecayTree");
  Tree_D2Kpimumu_magDw->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu_withGenLevel/magDw/MC12_DstD2Kpimumu_magDw.root");

  TChain* Tree_D2KKmumu_magUp = new TChain("MC12_DstD2KKMuMu/DecayTree");
  Tree_D2KKmumu_magUp->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_withGenLevel/magUp/MC12_DstD2KKmumu_magUp.root");
  TChain* Tree_D2KKmumu_magDw = new TChain("MC12_DstD2KKMuMu/DecayTree");
  Tree_D2KKmumu_magDw->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_withGenLevel/magDw/MC12_DstD2KKmumu_magDw.root");
 
  TChain* Tree_D2pipimumu_magUp = new TChain("MC12_DstD2PiPiMuMu/DecayTree");
  Tree_D2pipimumu_magUp->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu_withGenLevel/magUp/MC12_DstD2pipimumu_magUp.root");
  TChain* Tree_D2pipimumu_magDw = new TChain("MC12_DstD2PiPiMuMu/DecayTree");
  Tree_D2pipimumu_magDw->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu_withGenLevel/magDw/MC12_DstD2pipimumu_magDw.root");
   

  D2KpimumuReader* Kpi_Reader_magUp = new D2KpimumuReader(Tree_D2Kpimumu_magUp);
  D2KpimumuReader* Kpi_Reader_magDw = new D2KpimumuReader(Tree_D2Kpimumu_magDw);
  D2KKmumuReader* KK_Reader_magUp = new D2KKmumuReader(Tree_D2KKmumu_magUp);
  D2KKmumuReader* KK_Reader_magDw = new D2KKmumuReader(Tree_D2KKmumu_magDw);
  D2pipimumuReader* pipi_Reader_magUp = new D2pipimumuReader(Tree_D2pipimumu_magUp);
  D2pipimumuReader* pipi_Reader_magDw = new D2pipimumuReader(Tree_D2pipimumu_magDw);
  
  std::cout<<TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magUp.root",200.,1600.)<<std::endl;

  Kpi_Reader_magUp->createMCEfficiencyStudySampleNoTruthmatching(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magUp.root",200.,1600.),200.,1600.);
    Kpi_Reader_magDw->createMCEfficiencyStudySampleNoTruthmatching(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_%.1f_%.1f_magDw.root",200.,1600.),200.,1600.);

    

    KK_Reader_magUp->createMCEfficiencyStudySampleNoTruthmatching(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_magUp.root",200.,1600.),200.,1600.);
    KK_Reader_magDw->createMCEfficiencyStudySampleNoTruthmatching(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_%.1f_%.1f_magDw.root",200.,1600.),200.,1600.);
  
  pipi_Reader_magUp->createMCEfficiencyStudySampleNoTruthmatching(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_magUp.root",200.,1600.),200.,1600.);
    pipi_Reader_magDw->createMCEfficiencyStudySampleNoTruthmatching(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_%.1f_%.1f_magDw.root",200.,1600.),200.,1600.);



}




void addDiMuonMassBranch(TChain *chain, TString fOut, bool isReco){
   
  cout<<"add true Dimuon mass Branch for "<<fOut<<endl;
  double mu0_PX, mu0_PY, mu0_PZ;
  double mu1_PX, mu1_PY, mu1_PZ;
  double D_TRUE_DiMuon_Mass;
  double Dst_DTF_D0_M,Dst_DTF_Dstarplus_M;
  double deltaM;

  if(!isReco) {    
    chain->SetBranchAddress("muminus_TRUEP_X",&mu0_PX);
    chain->SetBranchAddress("muminus_TRUEP_Y",&mu0_PY);
    chain->SetBranchAddress("muminus_TRUEP_Z",&mu0_PZ);
  
    chain->SetBranchAddress("muplus_TRUEP_X",&mu1_PX);
    chain->SetBranchAddress("muplus_TRUEP_Y",&mu1_PY);
    chain->SetBranchAddress("muplus_TRUEP_Z",&mu1_PZ);
  }
  else{
    chain->SetBranchAddress("mu0_TRUEP_X",&mu0_PX);
    chain->SetBranchAddress("mu0_TRUEP_Y",&mu0_PY);
    chain->SetBranchAddress("mu0_TRUEP_Z",&mu0_PZ);

    chain->SetBranchAddress("mu1_TRUEP_X",&mu1_PX);
    chain->SetBranchAddress("mu1_TRUEP_Y",&mu1_PY);
    chain->SetBranchAddress("mu1_TRUEP_Z",&mu1_PZ);

    chain->SetBranchAddress("Dst_DTF_D0_M",&Dst_DTF_D0_M);
    chain->SetBranchAddress("Dst_DTF_Dstarplus_M",&Dst_DTF_Dstarplus_M);
    

  }
 
  TFile* targetFile = new TFile(fOut,"RECREATE");

  TTree* new_tree = chain->CopyTree("");
  TBranch* Bra = new_tree->Branch("D_TRUE_DiMuon_Mass",&D_TRUE_DiMuon_Mass);
  TBranch* Bra_dm;
  if(isReco) Bra_dm=new_tree->Branch("deltaM",&deltaM);
  TLorentzVector mu0, mu1;

  int numEvents = chain ->GetEntries();
  
  for(int i=0; i< numEvents; i++){
    
    chain->GetEntry(i);
    if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());    
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_TRUE_DiMuon_Mass=(mu0+mu1).M();
    
    std::cout<<D_TRUE_DiMuon_Mass<<"  "<<mu0_PX<<"  "<<mu0_PY<<"  "<<mu0_PZ<<"  "<<mu1_PX<<"  "<<mu1_PY<<"  "<<mu1_PZ<<
    
    Bra->Fill();
    if(isReco) {deltaM = Dst_DTF_Dstarplus_M-Dst_DTF_D0_M;Bra_dm->Fill();}
  } 
  
  new_tree->Write();
  targetFile->Write();
  targetFile->Close();
  //delete new_tree;
  delete targetFile;

}





void createTuplesForRecoEfficiency(){


  //also merge up and down!

  TChain *chain_KK_truth = new TChain("MCTruthTuple/MCTruthTuple");
  TChain *chain_Kpi_truth = new TChain("MCTruthTuple/MCTruthTuple");
  TChain *chain_pipi_truth = new TChain("MCTruthTuple/MCTruthTuple");
  
  TChain *chain_KK_reco = new TChain("MC12_DstD2KKMuMu/DecayTree");
  TChain *chain_Kpi_reco = new TChain("MC12_DstD2KPiMuMu/DecayTree");
  TChain *chain_pipi_reco = new TChain("MC12_DstD2pipiMuMu/DecayTree");

  //chain_pipi_truth->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/magUp/MC12_DstD2pipimumu_magUp.root");
  chain_pipi_reco->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/magUp/MC12_DstD2pipimumu_magUp.root");
  //chain_pipi_truth->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/magDw/MC12_DstD2pipimumu_magDw.root");
  chain_pipi_reco->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/magDw/MC12_DstD2pipimumu_magDw.root");

  //addDiMuonMassBranch(chain_pipi_truth,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MCTruthTuple_DstD2pipimumu.root",false);
  addDiMuonMassBranch(chain_pipi_reco,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MC12_DstD2pipimumu.root",true);
  
/*
  chain_Kpi_truth->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/magUp/MC12_DstD2Kpimumu_magUp.root");
  chain_Kpi_reco->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/magUp/MC12_DstD2Kpimumu_magUp.root");
  chain_Kpi_truth->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/magDw/MC12_DstD2Kpimumu_magDw.root");
  chain_Kpi_reco->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/magDw/MC12_DstD2Kpimumu_magDw.root");

  addDiMuonMassBranch(chain_Kpi_truth,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/MCTruthTuple_DstD2Kpimumu.root",false);
  addDiMuonMassBranch(chain_Kpi_reco,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/MC12_DstD2Kpimumu.root",true);

  chain_KK_truth->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/magUp/MC12_DstD2KKmumu_magUp.root");
  chain_KK_reco->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/magUp/MC12_DstD2KKmumu_magUp.root");
  chain_KK_truth->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/magDw/MC12_DstD2KKmumu_magDw.root");
  chain_KK_reco->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/magDw/MC12_DstD2KKmumu_magDw.root");

  addDiMuonMassBranch(chain_KK_truth,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/MCTruthTuple_DstD2KKmumu.root",false);
  addDiMuonMassBranch(chain_KK_reco,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/MC12_DstD2KKmumu.root",true);
*/ 

}

void createTuplesWithSelectedMultipleCand(){
  
  TChain *chain_KK_reco = new TChain("BDT_Tree");
  TChain *chain_Kpi_reco = new TChain("BDT_Tree");
  TChain *chain_Kpi_reco_D2KKmumuBDT = new TChain("BDT_Tree");
  TChain *chain_pipi_reco = new TChain("BDT_Tree");

  chain_pipi_reco->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_200.0_1600.0_D2pipimumuBDT.root");
  chain_KK_reco->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_200.0_1600.0_D2KKmumuBDT.root");
  chain_Kpi_reco->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2pipimumuBDT.root");
  chain_Kpi_reco_D2KKmumuBDT->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2KKmumuBDT.root");

  //addMultipleCandidateBranch(chain_KK_reco,"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_200.0_1600.0_D2KKmumuBDT_multipleCand_withGhosts.root","candidateLists/chosen_candidatesKK_2012_filtered_withGhosts.txt");
  //addMultipleCandidateBranch(chain_pipi_reco,"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_200.0_1600.0_D2pipimumuBDT_multipleCand_withGhosts.root","candidateLists/chosen_candidatespipi_2012_filtered_withGhosts.txt");
  //addMultipleCandidateBranch(chain_Kpi_reco,"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2pipimumuBDT_multipleCand_withGhosts.root","candidateLists/chosen_candidatesKpi_2012_filtered_D2pipimumuBDT_withGhosts.txt");
  //addMultipleCandidateBranch(chain_Kpi_reco_D2KKmumuBDT,"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2KKmumuBDT_multipleCand_withGhosts.root","candidateLists/chosen_candidatesKpi_2012_filtered_D2KKmumuBDT_withGhosts.txt");

  
  addMultipleCandidateBranch(chain_KK_reco,"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_200.0_1600.0_D2KKmumuBDT_multipleCand_noGhosts.root","candidateLists/chosen_candidatesKK_2012_filtered_noGhosts.txt");
  addMultipleCandidateBranch(chain_pipi_reco,"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_200.0_1600.0_D2pipimumuBDT_multipleCand_noGhosts.root","candidateLists/chosen_candidatespipi_2012_filtered_noGhosts.txt");
  addMultipleCandidateBranch(chain_Kpi_reco,"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2pipimumuBDT_multipleCand_noGhosts.root","candidateLists/chosen_candidatesKpi_2012_filtered_D2pipimumuBDT_noGhosts.txt");
  addMultipleCandidateBranch(chain_Kpi_reco_D2KKmumuBDT,"/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2KKmumuBDT_multipleCand_noGhosts.root","candidateLists/chosen_candidatesKpi_2012_filtered_D2KKmumuBDT_noGhosts.txt");
 
  
  //addMultipleCandidateBranch(chain_pipi_reco,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MC12_DstD2pipimumu_matched_andMultCand.root","chosen_candidatespipi_2012.txt");
  //addMultipleCandidateBranch(chain_Kpi_reco,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/MC12_DstD2Kpimumu_matched_andMultCand.root","chosen_candidatesKpi_2012.txt");

  // addCustomTruthMatchingBranch(chain_KK_reco,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/MC12_DstD2KKmumu_matched.root","AllTruel_candidatesKK_2012.txt");
  //addCustomTruthMatchingBranch(chain_Kpi_reco,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/MC12_DstD2Kpimumu_matched.root","AllTrue_candidatesKpi_2012.txt");
  //addCustomTruthMatchingBranch(chain_pipi_reco,"/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MC12_DstD2pipimumu_matched.root","AllTrue_candidatespipi_2012.txt");
  
}


std::pair<double,double> getMCEfficiency(TString fIn,TString nameTree, TString cut_sel,TString cut_norm ) {
  
  TFile* file= new TFile(fIn,"OPEN");
  TTree* tuple = (TTree*) file->Get(nameTree);
  
  double nNorm = (double)tuple->GetEntries(cut_norm);
  double nSel = (double)tuple->GetEntries(cut_sel+"&&"+cut_norm);

  double eff = (double)nSel/nNorm;
  double err = (double) 1/nNorm* TMath::Sqrt(nSel* (1- (nSel/nNorm) ) );

  tuple->Delete();
  file->Close();
  file->Delete();

  std::cout<<fIn<<" selection cut "<<cut_sel<<" norm cut "<<cut_norm<<" eff "<<(double)nSel/nNorm<<"  "<<nSel <<"  "<<nNorm<<std::endl;
  return std::make_pair(eff , err);

}

std::pair<double,double> getMCEfficiencyFrom2Files(TString fIn1,TString fIn2,TString nameTree, TString cut_sel,TString cut_norm ) {
  
  //  TFile* file1= new TFile(fIn1,"OPEN");
  //TTree* tuple1 = (TTree*) file1->Get(nameTree);

  //TFile* file2= new TFile(fIn2,"OPEN");
  //TTree* tuple2 = (TTree*) file2->Get(nameTree);
  
  TChain* myChain = new TChain(nameTree);
  myChain->AddFile(fIn1);
  myChain->AddFile(fIn2);

  //double nNorm = (double)tuple1->GetEntries(cut_norm) + (double)tuple2->GetEntries(cut_norm);
  //double nSel = (double)tuple1->GetEntries(cut_sel+"&&"+cut_norm) + (double)tuple2->GetEntries(cut_sel+"&&"+cut_norm);
  double nNorm = (double)myChain->GetEntries(cut_norm) ;
  double nSel = (double)myChain->GetEntries(cut_sel+"&&"+cut_norm);

  double eff = (double)nSel/nNorm;
  double err = (double) 1/nNorm* TMath::Sqrt(nSel* (1- (nSel/nNorm) ) );
  std::cout<<fIn1<<" selection cut "<<cut_sel<<" norm cut "<<cut_norm<<" eff "<<(double)nSel/nNorm<<"  "<<nSel <<"  "<<nNorm<<std::endl;

  //tuple1->Delete();
  //file1->Close();
  //file1->Delete();
  
  //tuple2->Delete();
  //file2->Close();
  //file2->Delete();
  myChain->Delete();

  return std::make_pair(eff , err);

}




std::pair<double,double> evaluatePIDCalibEfficiency(TString fIn) {

  TFile* file= new TFile(fIn,"OPEN");
  TTree* tuple = (TTree*) file->Get("CalibTool_PIDCalibTree");
  float weight=0;
  float err=0;
  double totalEff=0;
  double totalErr=0;

  double nEntries = (double)tuple->GetEntries();
  
  tuple->SetBranchAddress("Event_PIDCalibEff",&weight);
  tuple->SetBranchAddress("Event_PIDCalibErr",&err);

  for(int i=0;i<nEntries;++i) {
    tuple->GetEntry(i);
    totalEff+=(double)weight;
    totalErr+=(double)TMath::Power(err,2);;
    }

  totalEff=totalEff/nEntries;
  totalErr=TMath::Sqrt(totalErr)/nEntries;

  //std::cout<<"PID Calib Efficiency= "<<totalEff<<" for file "<<fIn<<std::endl; 
  tuple->Delete();
  file->Close();
  file->Delete();
  return std::make_pair(totalEff,totalErr);
}


void computeGlobalLowPtEfficiency(){

  dcastyle();

  TFile *fUp = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_magUp_ptBins.root","READ");
  TFile *fDw = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_magDw_ptBins.root","READ");

  TH1D * PIDCalib_extrapolatedEff_KKmumu_up =(TH1D*)fUp->Get("PIDCalib_extrapolatedEff_KKmumu");
  TH1D * PIDCalib_extrapolatedEff_pipimumu_up =(TH1D*)fUp->Get("PIDCalib_extrapolatedEff_pipimumu");
  TH1D * PIDCalib_extrapolatedEff_Kpimumu_up=(TH1D*)fUp->Get("PIDCalib_extrapolatedEff_Kpimumu");

  TH1D * PIDCalib_extrapolatedEff_KKmumu_dw =(TH1D*)fDw->Get("PIDCalib_extrapolatedEff_KKmumu");
  TH1D * PIDCalib_extrapolatedEff_pipimumu_dw =(TH1D*)fDw->Get("PIDCalib_extrapolatedEff_pipimumu");
  TH1D * PIDCalib_extrapolatedEff_Kpimumu_dw=(TH1D*)fDw->Get("PIDCalib_extrapolatedEff_Kpimumu");

  TFile *fout = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/globalLowPtEfficiency.root","RECREATE");
  TH1D * PIDCalib_extrapolatedEff_average = new TH1D("PIDCalib_extrapolatedEff_average","PIDCalib_extrapolatedEff_average",4,0,4);


  double average=0;
  double dAverage=0;
  double averageExtrapolEff;
  double dAverageExtrapolEff;

  for(int i=0; i<rangesKK_low.size();++i){

    //average extrapolated efficiency over m(mumu) and polarity                                                                                             
    averageExtrapolEff+=PIDCalib_extrapolatedEff_KKmumu_up->GetBinContent(i+1);
    dAverageExtrapolEff+=PIDCalib_extrapolatedEff_KKmumu_up->GetBinError(i+1);
    averageExtrapolEff+=PIDCalib_extrapolatedEff_KKmumu_dw->GetBinContent(i+1);
    dAverageExtrapolEff+=PIDCalib_extrapolatedEff_KKmumu_dw->GetBinError(i+1);
  }

  averageExtrapolEff/=double(2*rangesKK_low.size());
  dAverageExtrapolEff/=double(2*rangesKK_low.size());
  PIDCalib_extrapolatedEff_average->SetBinContent(1,averageExtrapolEff);
  PIDCalib_extrapolatedEff_average->SetBinError(1,dAverageExtrapolEff);
  averageExtrapolEff=0;
  dAverageExtrapolEff=0;

  //also average magUp and Dw for Kpimumu                                                                                                                                                     
  averageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_up->GetBinContent(1);
  dAverageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_up->GetBinError(1);
  averageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_dw->GetBinContent(1);
  dAverageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_dw->GetBinError(1);
  averageExtrapolEff/=2.;
  dAverageExtrapolEff/=2.;
  PIDCalib_extrapolatedEff_average->SetBinContent(3,averageExtrapolEff);
  PIDCalib_extrapolatedEff_average->SetBinError(3,dAverageExtrapolEff);
  averageExtrapolEff=0;
  dAverageExtrapolEff=0;


  for(int i=0; i<rangespipi_low.size();++i){

    //average extrapolated efficiency over m(mumu)                                                                                                                                            
    averageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_up->GetBinContent(i+1);
    dAverageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_up->GetBinError(i+1);
    averageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_dw->GetBinContent(i+1);
    dAverageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_dw->GetBinError(i+1);

  }

  averageExtrapolEff/=double(2*rangespipi_low.size());
  dAverageExtrapolEff/=double(2*rangespipi_low.size());
  PIDCalib_extrapolatedEff_average->SetBinContent(2,averageExtrapolEff);
  PIDCalib_extrapolatedEff_average->SetBinError(2,dAverageExtrapolEff);
  averageExtrapolEff=0;
  dAverageExtrapolEff=0;


  //average over all channels extrapolated low pt eff                                                                                                                                         
  averageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinContent(1);
  averageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinContent(2);
  averageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinContent(3);
  dAverageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinError(1);
  dAverageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinError(2);
  dAverageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinError(3);
  averageExtrapolEff/=3.;
  dAverageExtrapolEff/=3.;
  PIDCalib_extrapolatedEff_average->SetBinContent(4,averageExtrapolEff);
  PIDCalib_extrapolatedEff_average->SetBinError(4,dAverageExtrapolEff);


  TCanvas *f = new TCanvas("f","f");
  f->Divide(3,2);
  f->cd(1);
  PIDCalib_extrapolatedEff_KKmumu_up->Draw();
  PIDCalib_extrapolatedEff_KKmumu_dw->SetLineColor(kRed);
  PIDCalib_extrapolatedEff_KKmumu_dw->Draw("SAME");

  f->cd(2);
  PIDCalib_extrapolatedEff_pipimumu_up->Draw();
  PIDCalib_extrapolatedEff_pipimumu_dw->SetLineColor(kRed);
  PIDCalib_extrapolatedEff_pipimumu_dw->Draw("SAME");

  f->cd(3);
  PIDCalib_extrapolatedEff_Kpimumu_up->Draw();
  PIDCalib_extrapolatedEff_Kpimumu_dw->SetLineColor(kRed);
  PIDCalib_extrapolatedEff_Kpimumu_dw->Draw("SAME");

  f->cd(4);
  PIDCalib_extrapolatedEff_average->GetXaxis()->SetBinLabel(1,"D->KK#mu#mu");
  PIDCalib_extrapolatedEff_average->GetXaxis()->SetBinLabel(2,"D->#pi#pi#mu#mu");
  PIDCalib_extrapolatedEff_average->GetXaxis()->SetBinLabel(3,"D->K#pi#mu#mu");
  PIDCalib_extrapolatedEff_average->GetXaxis()->SetBinLabel(4,"average");
  PIDCalib_extrapolatedEff_average->Draw();
  PIDCalib_extrapolatedEff_average->Write();

  f->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/globalLowPtEfficiency.eps");
  fout->Write();
  fout->Close();


}

void drawBothPolaritiesPIDCalibMuonEfficiencyForPTBins(){

  //here, also a lot of things that are not needed (any more) are done. 

  dcastyle();

  TFile *fUp = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_magUp_ptBins.root","READ");
  TFile *fDw = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_magDw_ptBins.root","READ");

  TFile *fMC = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon0_PID_MC.root","READ");

  TF1 * f1= new TF1("f1","[0]",8000,80000);
  f1->SetLineColor(kBlue);
  TF1 * f2= new TF1("f2","[0]",8000,80000);
  f2->SetLineColor(kRed);
  
  std::vector<TH1D*> MC_PIDCalib_EffRatio_pipimumu_magUp;
  std::vector<TH1D*> MC_PIDCalib_EffRatio_KKmumu_magUp;
  std::vector<TH1D*> MC_PIDCalib_EffRatio_pipimumu_magDw;
  std::vector<TH1D*> MC_PIDCalib_EffRatio_KKmumu_magDw;

  std::vector<TH1D*> PIDCalib_Eff_pipimumu_magUp;
  std::vector<TH1D*> PIDCalib_Eff_KKmumu_magUp;
  std::vector<TH1D*> PIDCalib_Eff_pipimumu_magDw;
  std::vector<TH1D*> PIDCalib_Eff_KKmumu_magDw;

  TH1D * PIDCalib_extrapolatedEff_KKmumu_up =(TH1D*)fUp->Get("PIDCalib_extrapolatedEff_KKmumu");
  TH1D * PIDCalib_extrapolatedEff_pipimumu_up =(TH1D*)fUp->Get("PIDCalib_extrapolatedEff_pipimumu");
  TH1D * PIDCalib_extrapolatedEff_Kpimumu_up=(TH1D*)fUp->Get("PIDCalib_extrapolatedEff_Kpimumu");
 
  TH1D * PIDCalib_extrapolatedEff_KKmumu_dw =(TH1D*)fDw->Get("PIDCalib_extrapolatedEff_KKmumu");
  TH1D * PIDCalib_extrapolatedEff_pipimumu_dw =(TH1D*)fDw->Get("PIDCalib_extrapolatedEff_pipimumu");
  TH1D * PIDCalib_extrapolatedEff_Kpimumu_dw=(TH1D*)fDw->Get("PIDCalib_extrapolatedEff_Kpimumu");

  //also draw MC as comparison
  TH1* relMCEff_KKmumu = (TH1D*)fMC->Get("relMCEff_KKmumu_mu0");
  TH1* relMCEff_pipimumu =(TH1D*)fMC->Get("relMCEff_pipimumu_mu0");

  TFile *fout = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/muonPID_correctionFactors.root","RECREATE");

  TH1D * PIDCalib_extrapolatedEff_average = new TH1D("PIDCalib_extrapolatedEff_average","PIDCalib_extrapolatedEff_average",4,0,4);

 
  //average over pt and polarity
  TH1D * PIDCalib_average_KKmumu = new TH1D("PIDCalib_average_KKmumu","PIDCalib Muon Eff averaged",sizeof(binsKK)/sizeof(double)-1,binsKK); 
  TH1D * PIDCalib_average_pipimumu = new TH1D("PIDCalib_average_pipimumu","PIDCalib Muon Eff averaged",sizeof(binspipi)/sizeof(double)-1,binspipi); 

  TH1D * correctionFactors_KKmumu = new TH1D("correctionFactors","correctionFactors_KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK); 
  TH1D * correctionFactors_pipimumu = new TH1D("correctionFactors","correctionFactors_pipimumu",sizeof(binspipi)/sizeof(double)-1,binspipi); 
  TH1D * correctionFactors_KKmumu_magUp = new TH1D("correctionFactors_KKmumu_magUp","correctionFactors_KKmumu_magUp",sizeof(binsKK)/sizeof(double)-1,binsKK); 
  TH1D * correctionFactors_pipimumu_magUp = new TH1D("correctionFactors_pipimumu_magUp","correctionFactors_pipimumu_magUp",sizeof(binspipi)/sizeof(double)-1,binspipi); 
  TH1D * correctionFactors_KKmumu_magDw = new TH1D("correctionFactors_KKmumu_magDw","correctionFactors_KKmumu_magDw",sizeof(binsKK)/sizeof(double)-1,binsKK); 
  TH1D * correctionFactors_pipimumu_magDw = new TH1D("correctionFactors_pipimumu_magDw","correctionFactors_pipimumu_magDw",sizeof(binspipi)/sizeof(double)-1,binspipi); 

  TH1D* temp;

  for(int i=0; i<rangesKK_low.size();++i){
    temp = (TH1D*)fUp->Get(TString::Format("MC_data_PIDCalib_relEff_KKmumu_%i",i));
    MC_PIDCalib_EffRatio_KKmumu_magUp.push_back(temp);
    temp = (TH1D*)fDw->Get(TString::Format("MC_data_PIDCalib_relEff_KKmumu_%i",i));
    MC_PIDCalib_EffRatio_KKmumu_magDw.push_back(temp);
    //directly the PID calib eff ratios
    temp = (TH1D*)fUp->Get(TString::Format("relEff_KKmumu_%i",i));
    PIDCalib_Eff_KKmumu_magUp.push_back(temp);
    temp = (TH1D*)fDw->Get(TString::Format("relEff_KKmumu_%i",i));
    PIDCalib_Eff_KKmumu_magDw.push_back(temp);
  }

  for(int i=0; i<rangespipi_low.size();++i){
    temp = (TH1D*)fUp->Get(TString::Format("MC_data_PIDCalib_relEff_pipimumu_%i",i));
    MC_PIDCalib_EffRatio_pipimumu_magUp.push_back(temp);
    temp = (TH1D*)fDw->Get(TString::Format("MC_data_PIDCalib_relEff_pipimumu_%i",i));
    MC_PIDCalib_EffRatio_pipimumu_magDw.push_back(temp);
    //directly the PID calib eff ratios
    temp = (TH1D*)fUp->Get(TString::Format("relEff_pipimumu_%i",i));
    PIDCalib_Eff_pipimumu_magUp.push_back(temp);
    temp = (TH1D*)fDw->Get(TString::Format("relEff_pipimumu_%i",i));
    PIDCalib_Eff_pipimumu_magDw.push_back(temp);
  }


  TCanvas *a = new TCanvas("a","a");

  double average=0;
  double dAverage=0;
  double averageExtrapolEff;
  double dAverageExtrapolEff;
  

  a->Divide(2,2);

  for(int i=0; i<rangesKK_low.size();++i){
    a->cd(i+1);
    MC_PIDCalib_EffRatio_KKmumu_magUp[i]->Draw();
    MC_PIDCalib_EffRatio_KKmumu_magUp[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    MC_PIDCalib_EffRatio_KKmumu_magUp[i]->GetXaxis()->SetTitle("P_{t}(#mu)");
    MC_PIDCalib_EffRatio_KKmumu_magUp[i]->GetYaxis()->SetTitle("#epsilon_{R}^{PIDCalib}/#epsilon_{R}^{MC}");
    MC_PIDCalib_EffRatio_KKmumu_magUp[i]->Fit("f1","N");
    average+= f1->GetParameter(0);
    dAverage=+ f1->GetParError(0)*f1->GetParError(0)/2;

    correctionFactors_KKmumu_magUp->SetBinContent(i+1,f1->GetParameter(0));
    correctionFactors_KKmumu_magUp->SetBinError(i+1,f1->GetParError(0));

    MC_PIDCalib_EffRatio_KKmumu_magDw[i]->SetLineColor(kRed);
    MC_PIDCalib_EffRatio_KKmumu_magDw[i]->Draw("SAME");
    MC_PIDCalib_EffRatio_KKmumu_magDw[i]->Fit("f2","N");
    average+= f2->GetParameter(0);
    dAverage=+ f2->GetParError(0)*f2->GetParError(0)/2;
    average/=2; 
    dAverage = TMath::Sqrt(dAverage);

    correctionFactors_KKmumu_magDw->SetBinContent(i+1,f2->GetParameter(0));
    correctionFactors_KKmumu_magDw->SetBinError(i+1,f2->GetParError(0));

    correctionFactors_KKmumu->SetBinContent(i+1,average);
    correctionFactors_KKmumu->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;

    //average extrapolated efficiency over m(mumu) and polarity
    averageExtrapolEff+=PIDCalib_extrapolatedEff_KKmumu_up->GetBinContent(i+1);
    dAverageExtrapolEff+=PIDCalib_extrapolatedEff_KKmumu_up->GetBinError(i+1);
    averageExtrapolEff+=PIDCalib_extrapolatedEff_KKmumu_dw->GetBinContent(i+1);
    dAverageExtrapolEff+=PIDCalib_extrapolatedEff_KKmumu_dw->GetBinError(i+1);
  }
  averageExtrapolEff/=double(2*rangesKK_low.size());
  dAverageExtrapolEff/=double(2*rangesKK_low.size());
  PIDCalib_extrapolatedEff_average->SetBinContent(1,averageExtrapolEff);
  PIDCalib_extrapolatedEff_average->SetBinError(1,dAverageExtrapolEff);
  averageExtrapolEff=0;
  dAverageExtrapolEff=0;
  
  //also average magUp and Dw for Kpimumu
  averageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_up->GetBinContent(1);
  dAverageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_up->GetBinError(1);
  averageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_dw->GetBinContent(1);
  dAverageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_dw->GetBinError(1);
  averageExtrapolEff/=2.;
  dAverageExtrapolEff/=2.;
  PIDCalib_extrapolatedEff_average->SetBinContent(3,averageExtrapolEff);
  PIDCalib_extrapolatedEff_average->SetBinError(3,dAverageExtrapolEff);
  averageExtrapolEff=0;
  dAverageExtrapolEff=0;

  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/MC_PIDCalib_EffRatio_KKmumu.eps");
  a->Write();

  TCanvas *b = new TCanvas("b","b");

  b->Divide(3,2);
  for(int i=0; i<rangespipi_low.size();++i){
    b->cd(i+1);
    MC_PIDCalib_EffRatio_pipimumu_magUp[i]->Draw();
    MC_PIDCalib_EffRatio_pipimumu_magUp[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    MC_PIDCalib_EffRatio_pipimumu_magUp[i]->GetXaxis()->SetTitle("P_{t}(#mu)");
    MC_PIDCalib_EffRatio_pipimumu_magUp[i]->GetYaxis()->SetTitle("#epsilon_{R}^{PIDCalib}/#epsilon_{R}^{MC}");
    MC_PIDCalib_EffRatio_pipimumu_magUp[i]->Fit("f1","N");
    average+= f1->GetParameter(0);
    dAverage=+ f1->GetParError(0)*f1->GetParError(0)/2;

    correctionFactors_pipimumu_magUp->SetBinContent(i+1,f1->GetParameter(0));
    correctionFactors_pipimumu_magUp->SetBinError(i+1,f1->GetParError(0));

    MC_PIDCalib_EffRatio_pipimumu_magDw[i]->SetLineColor(kRed);
    MC_PIDCalib_EffRatio_pipimumu_magDw[i]->Draw("SAME");
    MC_PIDCalib_EffRatio_pipimumu_magDw[i]->Fit("f2","N");
    average+= f2->GetParameter(0);
    dAverage=+ f2->GetParError(0)*f2->GetParError(0)/2;
    average/=2; 
    dAverage = TMath::Sqrt(dAverage);

    correctionFactors_pipimumu_magDw->SetBinContent(i+1,f2->GetParameter(0));
    correctionFactors_pipimumu_magDw->SetBinError(i+1,f2->GetParError(0));

    correctionFactors_pipimumu->SetBinContent(i+1,average);
    correctionFactors_pipimumu->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;

    //average extrapolated efficiency over m(mumu)
    averageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_up->GetBinContent(i+1);
    dAverageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_up->GetBinError(i+1);
    averageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_dw->GetBinContent(i+1);
    dAverageExtrapolEff+=PIDCalib_extrapolatedEff_pipimumu_dw->GetBinError(i+1);
 
  }

  averageExtrapolEff/=double(2*rangespipi_low.size());
  dAverageExtrapolEff/=double(2*rangespipi_low.size());
  PIDCalib_extrapolatedEff_average->SetBinContent(2,averageExtrapolEff);
  PIDCalib_extrapolatedEff_average->SetBinError(2,dAverageExtrapolEff);
  averageExtrapolEff=0;
  dAverageExtrapolEff=0;
 
  b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/MC_PIDCalib_EffRatio_pipimumu.eps");
  b->Write();

  TCanvas *c = new TCanvas("c","c");
  c->Divide(2,2);
  c->cd(1);
  correctionFactors_KKmumu->GetYaxis()->SetRangeUser(0.95,1.05);
  correctionFactors_KKmumu->GetXaxis()->SetTitle("m(#mu #mu)[MeV]");
  correctionFactors_KKmumu->GetYaxis()->SetTitle("#epsilon_{R}^{PIDCalib}/#epsilon_{R}^{MC}");
  correctionFactors_KKmumu->Draw();
  correctionFactors_KKmumu->Write();
  c->cd(2);
  correctionFactors_pipimumu->GetYaxis()->SetRangeUser(0.95,1.05);
  correctionFactors_pipimumu->GetXaxis()->SetTitle("m(#mu #mu)[MeV]");
  correctionFactors_pipimumu->GetYaxis()->SetTitle("#epsilon_{R}^{PIDCalib}/#epsilon_{R}^{MC}");
  correctionFactors_pipimumu->Draw();
  correctionFactors_pipimumu->Write();

  c->cd(3);
  correctionFactors_KKmumu_magUp->GetYaxis()->SetRangeUser(0.95,1.05);
  correctionFactors_KKmumu_magUp->GetXaxis()->SetTitle("m(#mu #mu)[MeV]");
  correctionFactors_KKmumu_magUp->GetYaxis()->SetTitle("#epsilon_{R}^{PIDCalib}/#epsilon_{R}^{MC}");
  correctionFactors_KKmumu_magUp->Draw();
  correctionFactors_KKmumu_magUp->Write();
  correctionFactors_KKmumu_magDw->SetLineColor(kRed);
  correctionFactors_KKmumu_magDw->Draw("SAME");
  correctionFactors_KKmumu_magDw->Write();
 
  c->cd(4);
  correctionFactors_pipimumu_magUp->GetYaxis()->SetRangeUser(0.95,1.05);
  correctionFactors_pipimumu_magUp->GetXaxis()->SetTitle("m(#mu #mu)[MeV]");
  correctionFactors_pipimumu_magUp->GetYaxis()->SetTitle("#epsilon_{R}^{PIDCalib}/#epsilon_{R}^{MC}");
  correctionFactors_pipimumu_magUp->Draw();
  correctionFactors_pipimumu_magUp->Write();
  correctionFactors_pipimumu_magDw->SetLineColor(kRed);
  correctionFactors_pipimumu_magDw->Draw("SAME");
  correctionFactors_pipimumu_magDw->Write();
  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/MC_PIDCalib_correctionFactors.eps");
  c->Write();
   
  TCanvas *d = new TCanvas("d","d");
  d->Divide(2,2);

  for(int i=0; i<rangesKK_low.size();++i){
    d->cd(i+1);
    PIDCalib_Eff_KKmumu_magUp[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    PIDCalib_Eff_KKmumu_magUp[i]->GetXaxis()->SetTitle("m(#mu#mu)");
    PIDCalib_Eff_KKmumu_magUp[i]->GetYaxis()->SetTitle("#epsilon_{R}^{PIDCalib}");
    PIDCalib_Eff_KKmumu_magUp[i]->Draw();
    PIDCalib_Eff_KKmumu_magUp[i]->Fit("f1","N");
    average+= f1->GetParameter(0);
    dAverage=+ f1->GetParError(0)*f1->GetParError(0)/2;

    PIDCalib_Eff_KKmumu_magDw[i]->SetLineColor(kRed);
    PIDCalib_Eff_KKmumu_magDw[i]->Draw("SAME");
    PIDCalib_Eff_KKmumu_magDw[i]->Fit("f2","N");
    average+= f2->GetParameter(0);
    dAverage=+ f2->GetParError(0)*f2->GetParError(0)/2;
    average/=2; 
    dAverage = TMath::Sqrt(dAverage);

    PIDCalib_average_KKmumu->SetBinContent(i+1,average);
    PIDCalib_average_KKmumu->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
  } 
  d->cd(4);
  PIDCalib_average_KKmumu->GetYaxis()->SetRangeUser(0.97,1.03);
  relMCEff_KKmumu->SetLineColor(kCyan);
  PIDCalib_average_KKmumu->Draw();
  relMCEff_KKmumu->Draw("SAME");
  PIDCalib_average_KKmumu->Write();
 
  d->Print("test.eps");

  TCanvas *e = new TCanvas("e","e");
  e->Divide(3,2);

  for(int i=0; i<rangespipi_low.size();++i){
    e->cd(i+1);
    PIDCalib_Eff_pipimumu_magUp[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    PIDCalib_Eff_pipimumu_magUp[i]->GetXaxis()->SetTitle("p_{T}(#mu)");
    PIDCalib_Eff_pipimumu_magUp[i]->GetYaxis()->SetTitle("#epsilon_{R}^{PIDCalib}");
    PIDCalib_Eff_pipimumu_magUp[i]->Draw();
    PIDCalib_Eff_pipimumu_magUp[i]->Fit("f1","N");
    average+= f1->GetParameter(0);
    dAverage=+ f1->GetParError(0)*f1->GetParError(0)/2;

    PIDCalib_Eff_pipimumu_magDw[i]->SetLineColor(kRed);
    PIDCalib_Eff_pipimumu_magDw[i]->Draw("SAME");
    PIDCalib_Eff_pipimumu_magDw[i]->Fit("f2","N");
    average+= f2->GetParameter(0);
    dAverage=+ f2->GetParError(0)*f2->GetParError(0)/2;
    average/=2; 
    dAverage = TMath::Sqrt(dAverage);

    PIDCalib_average_pipimumu->SetBinContent(i+1,average);
    PIDCalib_average_pipimumu->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
  } 
  e->cd(6);
  PIDCalib_average_pipimumu->GetYaxis()->SetRangeUser(0.97,1.03);
  PIDCalib_average_pipimumu->Draw();
  relMCEff_pipimumu->SetLineColor(kCyan);
  relMCEff_pipimumu->Draw("SAME");
  PIDCalib_average_pipimumu->Write();
  e->Print("test2.eps");


  //average over all channels extrapolated low pt eff
  averageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinContent(1);
  averageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinContent(2);
  averageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinContent(3);
  dAverageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinError(1);
  dAverageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinError(2);
  dAverageExtrapolEff+=PIDCalib_extrapolatedEff_average->GetBinError(3);
  averageExtrapolEff/=3.;
  dAverageExtrapolEff/=3.;
  PIDCalib_extrapolatedEff_average->SetBinContent(4,averageExtrapolEff);
  PIDCalib_extrapolatedEff_average->SetBinError(4,dAverageExtrapolEff);

  TCanvas *f = new TCanvas("f","f");
  f->Divide(3,2);
  f->cd(1);
  PIDCalib_extrapolatedEff_KKmumu_up->Draw();
  PIDCalib_extrapolatedEff_KKmumu_dw->SetLineColor(kRed);
  PIDCalib_extrapolatedEff_KKmumu_dw->Draw("SAME");

  f->cd(2);
  PIDCalib_extrapolatedEff_pipimumu_up->Draw();
  PIDCalib_extrapolatedEff_pipimumu_dw->SetLineColor(kRed);
  PIDCalib_extrapolatedEff_pipimumu_dw->Draw("SAME");

  f->cd(3);
  PIDCalib_extrapolatedEff_Kpimumu_up->Draw();
  PIDCalib_extrapolatedEff_Kpimumu_dw->SetLineColor(kRed);
  PIDCalib_extrapolatedEff_Kpimumu_dw->Draw("SAME");
  
  f->cd(4);
  PIDCalib_extrapolatedEff_average->Draw();
  PIDCalib_extrapolatedEff_average->Write();

  f->Print("test3.eps");
  fout->Write();
  fout->Close();
 

 
}

double computeUncertaintyOnR(double a, double b, double c, double d, double e, double f, double g, double h, double i, double dB, double dC, double dE, double dG,double dI){

  double dRdC = -((1-d)*(a*b+(1-a)*c)*(c*(1-f)+f*g))/(TMath::Power((c*(1-d)+d*e),2)*(c*(1-h)+h*i))-((1-h)*(a*b+(1-a)*c)*(c*(1-f)+f*g))/((c*(1-d)+d*e)*TMath::Power((c*(1-h)+h*i),2))+((1-f)*(a*b+(1-a)*c))/((c*(1-d)+d*e)*(c*(1-h)+h*i))+((1-a)*(c*(1-f)+f*g))/((c*(1-d)+d*e)*(c*(1-h)+h*i));  
  double dRdB = (a*(c*(-f)+c+f*g))/((c*(d-1)-d*e)*(c*(h-1)-h*i));
  double dRdE = -(d*(a*b+(1-a)*c)*(c*(1-f)+f*g))/(TMath::Power((c*(1-d)+d*e),2)*(c*(1-h)+h*i));
  double dRdG = (f*(a*(b-c)+c))/((c*(d-1)-d*e)*(c*(h-1)-h*i));
  double dRdI = -(h*(a*b+(1-a)*c)*(c*(1-f)+f*g))/((c*(1-d)+d*e)*TMath::Power((c*(1-h)+h*i),2));
  

  double dR=TMath::Sqrt( TMath::Power((dRdC*dC),2) + TMath::Power((dRdB*dB),2) + TMath::Power((dRdE*dE),2)  + TMath::Power((dRdG*dG),2) + TMath::Power((dRdI*dI),2) );

  return dR;
}

void drawLowPtCorrectedEfficiency(TString binningScheme){

  TFile *fDataSinglePID_mu0 = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_bothPolarities_mu0_"+binningScheme+".root","READ");
  TFile *fDataSinglePID_mu1 = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_bothPolarities_mu1_"+binningScheme+".root","READ");
  TFile *fDataLowPt = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/globalLowPtEfficiency.root","READ");

  dcastyle();
  
  TH1* dataSingleEff_mu0_KKmumu = (TH1D*)fDataSinglePID_mu0->Get("AbsEff_KKmumu_average");
  TH1* dataSingleEff_mu0_pipimumu =(TH1D*)fDataSinglePID_mu0->Get("AbsEff_pipimumu_average");
  TH1* dataSingleEff_mu1_KKmumu = (TH1D*)fDataSinglePID_mu1->Get("AbsEff_KKmumu_average");
  TH1* dataSingleEff_mu1_pipimumu =(TH1D*)fDataSinglePID_mu1->Get("AbsEff_pipimumu_average");
  TH1* dataSingleEff_mu0_Kpimumu =(TH1D*)fDataSinglePID_mu0->Get("AbsEff_Kpimumu_average");
  TH1* dataSingleEff_mu1_Kpimumu =(TH1D*)fDataSinglePID_mu1->Get("AbsEff_Kpimumu_average");
  
  TH1*  extrapolatedLowPTEff = (TH1D*)fDataLowPt->Get("PIDCalib_extrapolatedEff_average");


  //these ones are corrected in low pt region!
  TH1* relEff_pipimumu = new TH1D("relEff_pipimumu","relEff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relEff_KKmumu = new TH1D("relEff_KKmumu","relEff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1* absEff_pipimumu = new TH1D("absEff_pipimumu","absEff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* absEff_KKmumu = new TH1D("absEff_KKmumu","absEff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* absEff_Kpimumu = new TH1D("absEff_Kpimumu","absEff_Kpimumu in bins of dimuon mass",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1* absSinglePIDEff_pipimumu = new TH1D("absSinglePIDEff_pipimumu","absSinglePIDEff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* absSinglePIDEff_KKmumu = new TH1D("absSinglePIDEff_KKmumu","absSinglePIDEff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* absSinglePIDEff_Kpimumu = new TH1D("absSinglePIDEff_Kpimumu","absSinglePIDEff_Kpimumu in bins of dimuon mass",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
 

  TFile *fSinglePID_mu0 = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon0_PID_MC.root","READ");
  TFile *fSinglePID_mu1 = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon1_PID_MC.root","READ");
  /*
  TH1* dataSingleEff_mu0_KKmumu = (TH1D*)fSinglePID_mu0->Get("AbsMCEff_KKmumu_lowPtCut_mu0");
  TH1* dataSingleEff_mu0_pipimumu =(TH1D*)fSinglePID_mu0->Get("AbsMCEff_pipimumu_lowPtCut_mu0");
  TH1* dataSingleEff_mu1_KKmumu = (TH1D*)fSinglePID_mu1->Get("AbsMCEff_KKmumu_lowPtCut_mu1");
  TH1* dataSingleEff_mu1_pipimumu =(TH1D*)fSinglePID_mu1->Get("AbsMCEff_pipimumu_lowPtCut_mu1");
  TH1* dataSingleEff_mu0_Kpimumu =(TH1D*)fSinglePID_mu0->Get("AbsMCEff_Kpimumu_lowPtCut_mu0");
  TH1* dataSingleEff_mu1_Kpimumu =(TH1D*)fSinglePID_mu1->Get("AbsMCEff_Kpimumu_lowPtCut_mu1");
  */
  TH1* MCSingleEff_mu0_KKmumu_ptCut = (TH1D*)fSinglePID_mu0->Get("AbsMCEff_KKmumu_highPtCut_mu0");
  TH1* MCSingleEff_mu0_pipimumu_ptCut =(TH1D*)fSinglePID_mu0->Get("AbsMCEff_pipimumu_highPtCut_mu0");
  TH1* MCSingleEff_mu1_KKmumu_ptCut = (TH1D*)fSinglePID_mu1->Get("AbsMCEff_KKmumu_highPtCut_mu1");
  TH1* MCSingleEff_mu1_pipimumu_ptCut =(TH1D*)fSinglePID_mu1->Get("AbsMCEff_pipimumu_highPtCut_mu1");

  TH1* MCSingleEff_mu0_Kpimumu_ptCut =(TH1D*)fSinglePID_mu0->Get("AbsMCEff_Kpimumu_highPtCut_mu0");
  TH1* MCSingleEff_mu1_Kpimumu_ptCut =(TH1D*)fSinglePID_mu1->Get("AbsMCEff_Kpimumu_highPtCut_mu1");

  TFile *fPtFractions_mu0 = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PtSpectrum_Muon0_MC.root","READ");
  TFile *fPtFractions_mu1 = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PtSpectrum_Muon1_MC.root","READ");
  
  std::vector<TH1D*> ptSpectrum_KKmumu_2Bins_mu0  ;
  std::vector<TH1D*> ptSpectrum_Kpimumu_2Bins_mu0;
  std::vector<TH1D*> ptSpectrum_pipimumu_2Bins_mu0;

  std::vector<TH1D*> ptSpectrum_KKmumu_2Bins_mu1  ;
  std::vector<TH1D*> ptSpectrum_Kpimumu_2Bins_mu1;
  std::vector<TH1D*> ptSpectrum_pipimumu_2Bins_mu1;
  TH1D* temp;


  for(int i=0; i<rangesKpi_low.size();++i){
    temp=(TH1D*)fPtFractions_mu0->Get(TString::Format("ptSpectrum_2Bins_Kpimumu_q2bin_%i_mu0",i));
    ptSpectrum_Kpimumu_2Bins_mu0.push_back(temp);
    temp=(TH1D*)fPtFractions_mu1->Get(TString::Format("ptSpectrum_2Bins_Kpimumu_q2bin_%i_mu1",i));
    ptSpectrum_Kpimumu_2Bins_mu1.push_back(temp);
  }

  for(int i=0; i<rangesKK_low.size();++i){ 
    temp=(TH1D*)fPtFractions_mu0->Get(TString::Format("ptSpectrum_2Bins_KKmumu_q2bin_%i_mu0",i));
    ptSpectrum_KKmumu_2Bins_mu0.push_back(temp);
    temp=(TH1D*)fPtFractions_mu1->Get(TString::Format("ptSpectrum_2Bins_KKmumu_q2bin_%i_mu1",i));
    ptSpectrum_KKmumu_2Bins_mu1.push_back(temp);
  }

  for(int i=0; i<rangespipi_low.size();++i){
    temp=(TH1D*)fPtFractions_mu0->Get(TString::Format("ptSpectrum_2Bins_pipimumu_q2bin_%i_mu0",i));
    ptSpectrum_pipimumu_2Bins_mu0.push_back(temp);
    temp=(TH1D*)fPtFractions_mu1->Get(TString::Format("ptSpectrum_2Bins_pipimumu_q2bin_%i_mu1",i));
    ptSpectrum_pipimumu_2Bins_mu1.push_back(temp);
  }
  std::cout<<"test1"<<std::endl;


  double Rmu1_tot,R_mu1_data_highPt,R_mu1_MC_lowPt;
  double Rmu0_tot,R_mu0_data_highPt,R_mu0_MC_lowPt;
  double fracPt_mu1, fracPt_mu0;
  double R_tot;
  double fsig_mu0,fsig_mu1;
  double fnorm_mu0,fnorm_mu1;
  double Effpipi_cor_mu0,EffKpi_cor_mu0;
  double Effpipi_cor_mu1,EffKpi_cor_mu1;
  double EffKK_cor_mu0;
  double EffKK_cor_mu1;
  
  for(int i=0; i<rangespipi_low.size();++i){
    

    fsig_mu0 = ptSpectrum_pipimumu_2Bins_mu0[i]->Integral(0,1)/ptSpectrum_pipimumu_2Bins_mu0[i]->Integral();
    fnorm_mu0 = ptSpectrum_Kpimumu_2Bins_mu0[0]->Integral(0,1)/ptSpectrum_Kpimumu_2Bins_mu0[0]->Integral();

    
    Effpipi_cor_mu0 = dataSingleEff_mu0_pipimumu->GetBinContent(i+1)*(1-fsig_mu0)+extrapolatedLowPTEff->GetBinContent(4)*fsig_mu0;
    EffKpi_cor_mu0 = dataSingleEff_mu0_Kpimumu->GetBinContent(1)*(1-fnorm_mu0)+extrapolatedLowPTEff->GetBinContent(4)*fnorm_mu0;
    
    Rmu0_tot=Effpipi_cor_mu0/EffKpi_cor_mu0;

    fsig_mu1 = ptSpectrum_pipimumu_2Bins_mu1[i]->Integral(0,1)/ptSpectrum_pipimumu_2Bins_mu1[i]->Integral();
    fnorm_mu1 = ptSpectrum_Kpimumu_2Bins_mu1[0]->Integral(0,1)/ptSpectrum_Kpimumu_2Bins_mu1[0]->Integral();

    Effpipi_cor_mu1 = dataSingleEff_mu1_pipimumu->GetBinContent(i+1)*(1-fsig_mu1)+extrapolatedLowPTEff->GetBinContent(4)*fsig_mu1;
    EffKpi_cor_mu1 = dataSingleEff_mu1_Kpimumu->GetBinContent(1)*(1-fnorm_mu1)+extrapolatedLowPTEff->GetBinContent(4)*fnorm_mu1;
     
    Rmu1_tot=Effpipi_cor_mu1/EffKpi_cor_mu1;
    R_tot=Rmu1_tot*Rmu0_tot;
    relEff_pipimumu->SetBinContent(i+1,R_tot);

    double EffData_mu0= dataSingleEff_mu0_pipimumu->GetBinContent(i+1);
    double EffNorm_mu0=dataSingleEff_mu0_Kpimumu->GetBinContent(1);
    double EffData_mu1= dataSingleEff_mu1_pipimumu->GetBinContent(i+1);
    double EffNorm_mu1=dataSingleEff_mu1_Kpimumu->GetBinContent(1);
    double EffExtr = extrapolatedLowPTEff->GetBinContent(4);

    double dEffData_mu0= dataSingleEff_mu0_pipimumu->GetBinError(i+1);
    double dEffNorm_mu0=dataSingleEff_mu0_Kpimumu->GetBinError(1);
    double dEffData_mu1= dataSingleEff_mu1_pipimumu->GetBinError(i+1);
    double dEffNorm_mu1=dataSingleEff_mu1_Kpimumu->GetBinError(1);
    double dEffExtr = extrapolatedLowPTEff->GetBinError(4);

    double dR = computeUncertaintyOnR(fsig_mu0,EffData_mu0,EffExtr,fnorm_mu0,EffNorm_mu0,fsig_mu1,EffData_mu1,fnorm_mu1,EffNorm_mu1,dEffData_mu0,dEffExtr,dEffNorm_mu0,dEffData_mu1,dEffNorm_mu1);
    relEff_pipimumu->SetBinError(i+1,dR);

    absEff_pipimumu->SetBinContent(i+1,Effpipi_cor_mu0*Effpipi_cor_mu1);
    absEff_Kpimumu->SetBinContent(1,EffKpi_cor_mu1);
    cout<<i<<" "<<Effpipi_cor_mu0<<"  "<<fsig_mu0<<"  "<< dataSingleEff_mu0_pipimumu->GetBinContent(i+1)<<"  "<<extrapolatedLowPTEff->GetBinContent(4)<<endl;

    absSinglePIDEff_pipimumu->SetBinContent(i+1,Effpipi_cor_mu0);
    absSinglePIDEff_Kpimumu->SetBinContent(1,EffKpi_cor_mu0);
 
    std::cout<<"dR "<<dR<<std::endl;

      }  

  
  for(int i=0; i<rangesKK_low.size();++i){
    

    fsig_mu0 = ptSpectrum_KKmumu_2Bins_mu0[i]->Integral(0,1)/ptSpectrum_KKmumu_2Bins_mu0[i]->Integral();
    fnorm_mu0 = ptSpectrum_Kpimumu_2Bins_mu0[0]->Integral(0,1)/ptSpectrum_Kpimumu_2Bins_mu0[0]->Integral();

    EffKK_cor_mu0 = dataSingleEff_mu0_KKmumu->GetBinContent(i+1)*(1-fsig_mu0)+extrapolatedLowPTEff->GetBinContent(4)*fsig_mu0;
    EffKpi_cor_mu0 = dataSingleEff_mu0_Kpimumu->GetBinContent(1)*(1-fnorm_mu0)+extrapolatedLowPTEff->GetBinContent(4)*fnorm_mu0;
    
    Rmu0_tot=EffKK_cor_mu0/EffKpi_cor_mu0;

    fsig_mu1 = ptSpectrum_KKmumu_2Bins_mu1[i]->Integral(0,1)/ptSpectrum_KKmumu_2Bins_mu1[i]->Integral();
    fnorm_mu1 = ptSpectrum_Kpimumu_2Bins_mu1[0]->Integral(0,1)/ptSpectrum_Kpimumu_2Bins_mu1[0]->Integral();

    EffKK_cor_mu1 = dataSingleEff_mu1_KKmumu->GetBinContent(i+1)*(1-fsig_mu1)+extrapolatedLowPTEff->GetBinContent(4)*fsig_mu1;
    EffKpi_cor_mu1 = dataSingleEff_mu1_Kpimumu->GetBinContent(1)*(1-fnorm_mu1)+extrapolatedLowPTEff->GetBinContent(4)*fnorm_mu1;
     
    Rmu1_tot=EffKK_cor_mu1/EffKpi_cor_mu1;
    R_tot=Rmu1_tot*Rmu0_tot;
    relEff_KKmumu->SetBinContent(i+1,R_tot);

    double EffData_mu0= dataSingleEff_mu0_pipimumu->GetBinContent(i+1);
    double EffNorm_mu0=dataSingleEff_mu0_Kpimumu->GetBinContent(1);
    double EffData_mu1= dataSingleEff_mu1_pipimumu->GetBinContent(i+1);
    double EffNorm_mu1=dataSingleEff_mu1_Kpimumu->GetBinContent(1);
    double EffExtr = extrapolatedLowPTEff->GetBinContent(4);

    double dEffData_mu0= dataSingleEff_mu0_pipimumu->GetBinError(i+1);
    double dEffNorm_mu0=dataSingleEff_mu0_Kpimumu->GetBinError(1);
    double dEffData_mu1= dataSingleEff_mu1_pipimumu->GetBinError(i+1);
    double dEffNorm_mu1=dataSingleEff_mu1_Kpimumu->GetBinError(1);
    double dEffExtr = extrapolatedLowPTEff->GetBinError(4);

    double dR = computeUncertaintyOnR(fsig_mu0,EffData_mu0,EffExtr,fnorm_mu0,EffNorm_mu0,fsig_mu1,EffData_mu1,fnorm_mu1,EffNorm_mu1,dEffData_mu0,dEffExtr,dEffNorm_mu0,dEffData_mu1,dEffNorm_mu1);
    relEff_KKmumu->SetBinError(i+1,dR);
    
    absEff_KKmumu->SetBinContent(i+1,EffKK_cor_mu0*EffKK_cor_mu1);
    absSinglePIDEff_KKmumu->SetBinContent(i+1,EffKK_cor_mu0);

  }  

  TFile* fout = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalMuonPIDEfficiency_"+binningScheme+".root","RECREATE");
  TCanvas* a = new TCanvas("a","a");
  a->Divide(1,2);
  a->cd(1);
  relEff_pipimumu->GetYaxis()->SetRangeUser(0.9,1.2);
  relEff_pipimumu->Draw();
  relEff_pipimumu->Write();
  a->cd(2);
  relEff_KKmumu->GetYaxis()->SetRangeUser(0.9,1.2);
  relEff_KKmumu->Draw();
  relEff_KKmumu->Write();
  if(binningScheme=="default")
    a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalMuonEfficiency.eps");
  
  absEff_KKmumu->Write();
  absEff_Kpimumu->Write();
  absEff_pipimumu->Write();
  absSinglePIDEff_pipimumu->Write();
  absSinglePIDEff_KKmumu->Write();
  absSinglePIDEff_Kpimumu->Write();


  fout->Write();
  fout->Close();
} 



void checkMuonPIDFactorization(){

  //check if Eff(m1&&mu2) = Eff(mu1)* Eff(mu2) in MC and compare to PID calib result

  dcastyle();

  TFile *fSinglePID_mu0 = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon0_PID_MC.root","READ");
  TFile *fSinglePID_mu1 = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon1_PID_MC.root","READ");
  TFile *fDoublePID = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC.root","READ");

  TFile *fDataSinglePID_mu0 = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_bothPolarities_mu0.root","READ");
  TFile *fDataSinglePID_mu1 = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_bothPolarities_mu1.root","READ");
  TFile *fDataDoublePID = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_and_PIDCalib_magDw.root","READ");

  TFile *flowPtCoroected = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalMuonPIDEfficiency.root","READ");
  TH1* h_corrected_pipi=(TH1D*)flowPtCoroected->Get("relEff_pipimumu");
  TH1* h_corrected_KK=(TH1D*)flowPtCoroected->Get("relEff_KKmumu");

  //TFile *fDataDoublePID = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/Muon_PID_MC_and_PIDCalib_bothPolarities.root","READ");
  
  /*
  TH1* MCSingleEff_mu0_KKmumu = (TH1D*)fSinglePID_mu0->Get("relMCEff_KKmumu_mu0");
  TH1* MCSingleEff_mu0_pipimumu =(TH1D*)fSinglePID_mu0->Get("relMCEff_pipimumu_mu0");
  TH1* MCSingleEff_mu1_KKmumu = (TH1D*)fSinglePID_mu1->Get("relMCEff_KKmumu_mu1");
  TH1* MCSingleEff_mu1_pipimumu =(TH1D*)fSinglePID_mu1->Get("relMCEff_pipimumu_mu1");

  TH1* MCDoubleEff_KKmumu = (TH1D*)fDoublePID->Get("relMCEff_KKmumu");
  TH1* MCDoubleEff_pipimumu =(TH1D*)fDoublePID->Get("relMCEff_pipimumu");
  */

  TH1* MCSingleEff_mu0_KKmumu_ptCut = (TH1D*)fSinglePID_mu0->Get("relMCEff_KKmumu_highPtCut_mu0");
  TH1* MCSingleEff_mu0_pipimumu_ptCut =(TH1D*)fSinglePID_mu0->Get("relMCEff_pipimumu_highPtCut_mu0");
  TH1* MCSingleEff_mu1_KKmumu_ptCut = (TH1D*)fSinglePID_mu1->Get("relMCEff_KKmumu_highPtCut_mu1");
  TH1* MCSingleEff_mu1_pipimumu_ptCut =(TH1D*)fSinglePID_mu1->Get("relMCEff_pipimumu_highPtCut_mu1");

  TH1* MCSingleEff_mu0_KKmumu = (TH1D*)fSinglePID_mu0->Get("relMCEff_KKmumu_mu0");
  TH1* MCSingleEff_mu0_pipimumu =(TH1D*)fSinglePID_mu0->Get("relMCEff_pipimumu_mu0");
  TH1* MCSingleEff_mu1_KKmumu = (TH1D*)fSinglePID_mu1->Get("relMCEff_KKmumu_mu1");
  TH1* MCSingleEff_mu1_pipimumu =(TH1D*)fSinglePID_mu1->Get("relMCEff_pipimumu_mu1");


  TH1* MCDoubleEff_KKmumu_ptCut = (TH1D*)fDoublePID->Get("relMCEff_KKmumu_highPtCut");
  TH1* MCDoubleEff_pipimumu_ptCut =(TH1D*)fDoublePID->Get("relMCEff_pipimumu_highPtCut");


  TH1* MCDoubleEff_KKmumu = (TH1D*)fDoublePID->Get("relMCEff_KKmumu");
  TH1* MCDoubleEff_pipimumu =(TH1D*)fDoublePID->Get("relMCEff_pipimumu");

 
  TH1* dataSingleEff_mu0_KKmumu = (TH1D*)fDataSinglePID_mu0->Get("relEff_KKmumu_average");
  TH1* dataSingleEff_mu0_pipimumu =(TH1D*)fDataSinglePID_mu0->Get("relEff_pipimumu_average");

  TH1* dataSingleEff_mu1_KKmumu = (TH1D*)fDataSinglePID_mu1->Get("relEff_KKmumu_average");
  TH1* dataSingleEff_mu1_pipimumu =(TH1D*)fDataSinglePID_mu1->Get("relEff_pipimumu_average");

 
  //TH1* dataDoubleEff_KKmumu = (TH1D*)fDataDoublePID->Get("relEff_KKmumu_average");
  //TH1* dataDoubleEff_pipimumu =(TH1D*)fDataDoublePID->Get("relEff_pipimumu_average");
  TH1* dataDoubleEff_KKmumu = (TH1D*)fDataDoublePID->Get("relEff_KKmumu");
  TH1* dataDoubleEff_pipimumu =(TH1D*)fDataDoublePID->Get("relEff_pipimumu");
 
  //TH1D * dataDoubleEff_KKmumu = new TH1D("dataDoubleEff_KKmumu","dataDoubleEff_KKmumu (Eff1^2)",sizeof(binsKK)/sizeof(double)-1,binsKK); 
  //TH1D * dataDoubleEff_KKmumu = new TH1D("dataDoubleEff_KKmumu","dataDoubleEff_KKmumu (Eff1^2)",sizeof(binsKK)/sizeof(double)-1,binsKK); 

  TH1D * MC_product_Eff_KKmumu = new TH1D("MC_product_Eff_KKmumu","MC_product_Eff_KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK); 
  TH1D * MC_product_Eff_pipimumu = new TH1D("MC_product_Eff_pipimumu","MC_product_Eff_pipimumu",sizeof(binspipi)/sizeof(double)-1,binspipi); 

  TH1D * MC_product_Eff_KKmumu_ptCut = new TH1D("MC_product_Eff_KKmumu_ptCut","MC_product_Eff_KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK); 
  TH1D * MC_product_Eff_pipimumu_ptCut = new TH1D("MC_product_Eff_pipimumu_ptCut","MC_product_Eff_pipimumu",sizeof(binspipi)/sizeof(double)-1,binspipi); 

  TH1D * data_product_Eff_KKmumu = new TH1D("data_product_Eff_KKmumu","data_product_Eff_KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK); 
  TH1D * data_product_Eff_pipimumu = new TH1D("data_product_Eff_pipimumu","data_product_Eff_pipimumu",sizeof(binspipi)/sizeof(double)-1,binspipi); 
    
  
  TFile *fout = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/checkMuonPIDFactorization.root","RECREATE");


  double prodEff;
  double dProdEff;
  

  for(int i=0; i<rangesKK_low.size();++i){

    //MC
    prodEff=MCSingleEff_mu0_KKmumu->GetBinContent(i+1) * MCSingleEff_mu1_KKmumu->GetBinContent(i+1);
    MC_product_Eff_KKmumu->SetBinContent(i+1,prodEff);
    dProdEff=TMath::Sqrt(TMath::Power(MCSingleEff_mu0_KKmumu->GetBinContent(i+1)*MCSingleEff_mu1_KKmumu->GetBinError(i+1),2)+TMath::Power(MCSingleEff_mu1_KKmumu->GetBinContent(i+1)*MCSingleEff_mu0_KKmumu->GetBinError(i+1),2));
    MC_product_Eff_KKmumu->SetBinError(i+1,dProdEff);
    dProdEff=0;
    prodEff=0;

    //data
    prodEff=dataSingleEff_mu0_KKmumu->GetBinContent(i+1)*dataSingleEff_mu1_KKmumu->GetBinContent(i+1);
    data_product_Eff_KKmumu->SetBinContent(i+1,prodEff);
    dProdEff=TMath::Sqrt(TMath::Power(dataSingleEff_mu0_KKmumu->GetBinContent(i+1)*dataSingleEff_mu1_KKmumu->GetBinError(i+1),2)+TMath::Power(dataSingleEff_mu1_KKmumu->GetBinContent(i+1)*dataSingleEff_mu0_KKmumu->GetBinError(i+1),2));
    data_product_Eff_KKmumu->SetBinError(i+1,dProdEff);
    dProdEff=0;
    prodEff=0;

  }  
  TCanvas *a = new TCanvas("a","a"); 
  MC_product_Eff_KKmumu->GetYaxis()->SetRangeUser(0.9,1.);
  MC_product_Eff_KKmumu->Draw();
  MCDoubleEff_KKmumu->SetLineColor(kBlue);
  MCDoubleEff_KKmumu->Draw("SAME");
  dataDoubleEff_KKmumu->SetLineColor(kOrange);
  //dataDoubleEff_KKmumu->Draw("SAME");
  data_product_Eff_KKmumu->SetLineColor(kViolet);
  data_product_Eff_KKmumu->Draw("SAME");
  h_corrected_KK->SetLineColor(kRed);
  h_corrected_KK->Draw("SAME");

  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/checkMuonPIDFactorization1.eps");

  for(int i=0; i<rangespipi_low.size();++i){
    prodEff=MCSingleEff_mu0_pipimumu->GetBinContent(i+1) * MCSingleEff_mu1_pipimumu->GetBinContent(i+1);
    MC_product_Eff_pipimumu->SetBinContent(i+1,prodEff);
    dProdEff=TMath::Sqrt(TMath::Power(MCSingleEff_mu0_pipimumu->GetBinContent(i+1)*MCSingleEff_mu1_pipimumu->GetBinError(i+1),2)+TMath::Power(MCSingleEff_mu1_pipimumu->GetBinContent(i+1)*MCSingleEff_mu0_pipimumu->GetBinError(i+1),2));
    MC_product_Eff_pipimumu->SetBinError(i+1,dProdEff);
    dProdEff=0;
    prodEff=0;
    //data
    prodEff=dataSingleEff_mu0_pipimumu->GetBinContent(i+1)*dataSingleEff_mu1_pipimumu->GetBinContent(i+1);
    //data_product_Eff_KKmumu->SetBinContent(i+1,prodEff);?????
    data_product_Eff_pipimumu->SetBinContent(i+1,prodEff);
    dProdEff=TMath::Sqrt(TMath::Power(dataSingleEff_mu0_pipimumu->GetBinContent(i+1)*dataSingleEff_mu1_pipimumu->GetBinError(i+1),2)+TMath::Power(dataSingleEff_mu1_pipimumu->GetBinContent(i+1)*dataSingleEff_mu0_pipimumu->GetBinError(i+1),2));
    data_product_Eff_pipimumu->SetBinError(i+1,dProdEff);
    dProdEff=0;
    prodEff=0;
    //MC pt Cut
    prodEff=MCSingleEff_mu0_pipimumu_ptCut->GetBinContent(i+1) * MCSingleEff_mu1_pipimumu_ptCut->GetBinContent(i+1);
    MC_product_Eff_pipimumu_ptCut->SetBinContent(i+1,prodEff);
    dProdEff=TMath::Sqrt(TMath::Power(MCSingleEff_mu0_pipimumu_ptCut->GetBinContent(i+1)*MCSingleEff_mu1_pipimumu_ptCut->GetBinError(i+1),2)+TMath::Power(MCSingleEff_mu1_pipimumu_ptCut->GetBinContent(i+1)*MCSingleEff_mu0_pipimumu_ptCut->GetBinError(i+1),2));
    MC_product_Eff_pipimumu_ptCut->SetBinError(i+1,dProdEff);
    dProdEff=0;
    prodEff=0;

  }  
  TCanvas *b = new TCanvas("b","b");
  MC_product_Eff_pipimumu->GetYaxis()->SetRangeUser(0.9,1.2);
  MC_product_Eff_pipimumu->Draw();
  MCDoubleEff_pipimumu->SetLineColor(kBlue);
  MCDoubleEff_pipimumu->Draw("SAME");
  //dataDoubleEff_pipimumu->SetLineColor(kOrange);
  //dataDoubleEff_pipimumu->Draw("SAME");
  data_product_Eff_pipimumu->SetLineColor(kViolet);
  data_product_Eff_pipimumu->Draw("SAME");
  //MC_product_Eff_pipimumu_ptCut->SetLineColor(kGreen);
  //MC_product_Eff_pipimumu_ptCut->Draw("SAME");
  //MCDoubleEff_pipimumu_ptCut->SetLineColor(kGreen+2);;
  //MCDoubleEff_pipimumu_ptCut->Draw("SAME");
  h_corrected_pipi->SetLineColor(kRed);
  h_corrected_pipi->Draw("SAME");
  
  b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/checkMuonPIDFactorization2.eps");
}


void drawBothPolaritiesPIDCalibHadronEfficiency(TString binningScheme="default"){

  dcastyle();

  TFile *fUp = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Hadron_PID_MC_and_PIDCalib_"+binningScheme+"_magUp.root","READ");
  TFile *fDw = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Hadron_PID_MC_and_PIDCalib_"+binningScheme+"_magDw.root","READ");

  TH1* relEff_KKmumu_up =(TH1D*)fUp->Get("relEff_KKmumu"); 
  TH1* relMCEff_KKmumu_up= (TH1D*)fUp->Get("relMCEff_KKmumu");
  TH1* MC_PIDCalib_EffRatio_KKmumu_up=(TH1D*)fUp->Get("MC_PIDCalib_EffRatio_KKmumu");
  TH1* relEff_pipimumu_up=(TH1D*)fUp->Get("relEff_pipimumu");
  TH1* relMCEff_pipimumu_up=(TH1D*)fUp->Get("relMCEff_pipimumu");
  TH1* MC_PIDCalib_EffRatio_pipimumu_up=(TH1D*)fUp->Get("MC_PIDCalib_EffRatio_pipimumu");

  TH1* relEff_KKmumu_dw =(TH1D*)fDw->Get("relEff_KKmumu"); 
  TH1* relMCEff_KKmumu_dw= (TH1D*)fDw->Get("relMCEff_KKmumu");
  TH1* MC_PIDCalib_EffRatio_KKmumu_dw=(TH1D*)fDw->Get("MC_PIDCalib_EffRatio_KKmumu");
  TH1* relEff_pipimumu_dw=(TH1D*)fDw->Get("relEff_pipimumu");
  TH1* relMCEff_pipimumu_dw=(TH1D*)fDw->Get("relMCEff_pipimumu");
  TH1* MC_PIDCalib_EffRatio_pipimumu_dw=(TH1D*)fDw->Get("MC_PIDCalib_EffRatio_pipimumu");
  
  TH1* relEff_KKmumu_average = new TH1D("relEff_KKmumu_average","relEff_KKmumu in bins of dimuon mass averaged both polarities",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relEff_pipimumu_average = new TH1D("relEff_pipimumu_average","relEff_pipimumu in bins of dimuon mass averaged both polarities",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1* relMCEff_KKmumu_average = new TH1D("relMCEff_KKmumu_average","relEff_KKmumu in bins of dimuon mass averaged both polarities",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_pipimumu_average = new TH1D("relMCEff_pipimumu_average","relEff_pipimumu in bins of dimuon mass averaged both polarities",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TFile *fout = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalHadronPIDEfficiency_"+binningScheme+".root","RECREATE");

  TCanvas *a = new TCanvas("a","a");

  double average=0;
  double dAverage=0;

  a->Divide(2,2);
  a->cd(1);
  relEff_KKmumu_up->Draw();
  relEff_KKmumu_dw->SetLineColor(kRed);
  relEff_KKmumu_dw->Draw("SAME");
  a->cd(2);
  relMCEff_KKmumu_up->GetYaxis()->SetRangeUser(0.9,1);
  relMCEff_KKmumu_up->Draw();
  relMCEff_KKmumu_dw->SetLineColor(kRed);
  relMCEff_KKmumu_dw->Draw("SAME");
  a->cd(3);
  relEff_pipimumu_up->Draw();
  relEff_pipimumu_dw->SetLineColor(kRed);
  relEff_pipimumu_dw->Draw("SAME");
  a->cd(4);
  relMCEff_pipimumu_up->GetYaxis()->SetRangeUser(1.02,1.14);
  relMCEff_pipimumu_up->Draw();
  relMCEff_pipimumu_dw->SetLineColor(kRed);
  relMCEff_pipimumu_dw->Draw("SAME");
  if(binningScheme=="default")
    a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Hadron_PID_MC_and_PIDCalib_bothPolarities1.eps");

  TCanvas *b = new TCanvas("b","b");
  b->Divide(1,2);
  b->cd(1);
  MC_PIDCalib_EffRatio_KKmumu_up->Draw();
  MC_PIDCalib_EffRatio_KKmumu_dw->SetLineColor(kRed);
  MC_PIDCalib_EffRatio_KKmumu_dw->Draw("SAME");
  b->cd(2);
  MC_PIDCalib_EffRatio_pipimumu_up->Draw();
  MC_PIDCalib_EffRatio_pipimumu_dw->SetLineColor(kRed);
  MC_PIDCalib_EffRatio_pipimumu_dw->Draw("SAME");
  if(binningScheme=="default")
    b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Hadron_PID_MC_and_PIDCalib_bothPolarities2.eps");

  //average over polarities
  for(int i=0; i<rangesKK_low.size();++i){
    average+=relEff_KKmumu_up->GetBinContent(i+1);
    dAverage=+relEff_KKmumu_up->GetBinError(i+1)*relEff_KKmumu_up->GetBinError(i+1)/2;
    average+=relEff_KKmumu_dw->GetBinContent(i+1);
    dAverage=+relEff_KKmumu_dw->GetBinError(i+1)*relEff_KKmumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                               
    dAverage = TMath::Sqrt(dAverage);
    relEff_KKmumu_average->SetBinContent(i+1,average);
    relEff_KKmumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
    //MC
    average+=relMCEff_KKmumu_up->GetBinContent(i+1);
    dAverage=+relMCEff_KKmumu_up->GetBinError(i+1)*relMCEff_KKmumu_up->GetBinError(i+1)/2;
    average+=relMCEff_KKmumu_dw->GetBinContent(i+1);
    dAverage=+relMCEff_KKmumu_dw->GetBinError(i+1)*relMCEff_KKmumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                
    dAverage = TMath::Sqrt(dAverage);
    relMCEff_KKmumu_average->SetBinContent(i+1,average);
    relMCEff_KKmumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;

  }    

  for(int i=0; i<rangespipi_low.size();++i){
    average+=relEff_pipimumu_up->GetBinContent(i+1);
    dAverage=+relEff_pipimumu_up->GetBinError(i+1)*relEff_pipimumu_up->GetBinError(i+1)/2;
    average+=relEff_pipimumu_dw->GetBinContent(i+1);
    dAverage=+relEff_pipimumu_dw->GetBinError(i+1)*relEff_pipimumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                               
    dAverage = TMath::Sqrt(dAverage);
    relEff_pipimumu_average->SetBinContent(i+1,average);
    relEff_pipimumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
    //MC
    average+=relMCEff_pipimumu_up->GetBinContent(i+1);
    dAverage=+relMCEff_pipimumu_up->GetBinError(i+1)*relMCEff_pipimumu_up->GetBinError(i+1)/2;
    average+=relMCEff_pipimumu_dw->GetBinContent(i+1);
    dAverage=+relMCEff_pipimumu_dw->GetBinError(i+1)*relMCEff_pipimumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                               
    dAverage = TMath::Sqrt(dAverage);
    relMCEff_pipimumu_average->SetBinContent(i+1,average);
    relMCEff_pipimumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;

  }    

  TCanvas *c = new TCanvas("c","c");
  c->Divide(1,2);
  c->cd(1);
  relEff_KKmumu_average->GetYaxis()->SetRangeUser(0.9,1.);
  relEff_KKmumu_average->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_KKmumu_average->Draw(); 
  relEff_KKmumu_average->Write(); 
  relMCEff_KKmumu_average->SetLineColor(kCyan); 
  relMCEff_KKmumu_average->Draw("SAME"); 
  relMCEff_KKmumu_average->Write(); 
  c->cd(2);
  relEff_pipimumu_average->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_pipimumu_average->GetYaxis()->SetRangeUser(1.0,1.1);
  relEff_pipimumu_average->Draw();
  relEff_pipimumu_average->Write();
  relMCEff_pipimumu_average->SetLineColor(kCyan); 
  relMCEff_pipimumu_average->Draw("SAME"); 
  relMCEff_pipimumu_average->Write(); 
  if(binningScheme=="default")
    c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalHadronEfficiency.eps");

  fout->Write();
  fout->Close();
}

void drawTotalPIDEfficiency(TString binningScheme, double dEffMu_sys, double dEffHad_sys) {

  
  TFile *fMuon = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalMuonPIDEfficiency_"+binningScheme+".root","READ");
  TFile *fHadron = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalHadronPIDEfficiency_"+binningScheme+".root","READ");
  TFile *fMCMuon = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC.root","READ");
		  
  TH1* muonPID_pipi=(TH1D*)fMuon->Get("relEff_pipimumu");
  TH1* muonPID_KK=(TH1D*)fMuon->Get("relEff_KKmumu");

  TH1* hadronPID_pipi=(TH1D*)fHadron->Get("relEff_pipimumu_average");
  TH1* hadronPID_KK=(TH1D*)fHadron->Get("relEff_KKmumu_average");

  TH1* MChadronPID_pipi=(TH1D*)fHadron->Get("relMCEff_pipimumu_average");
  TH1* MChadronPID_KK=(TH1D*)fHadron->Get("relMCEff_KKmumu_average");
  
  TH1* MCmuonPID_pipi=(TH1D*)fMCMuon->Get("relMCEff_pipimumu");
  TH1* MCmuonPID_KK=(TH1D*)fMCMuon->Get("relMCEff_KKmumu");

  TFile* fOut = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/totalPIDEfficiency_"+binningScheme+".root","RECREATE"); 

  TH1* totalRelPIDEff_KKmumu = new TH1D("totalRelPIDEff_KKmumu","rel PID Eff KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* totalRelPIDEff_pipimumu = new TH1D("totalRelPIDEff_pipimumu","rel PID Eff pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
    
  TH1* totalRelMCPIDEff_KKmumu = new TH1D("totalRelMCPIDEff_KKmumu","rel PID Eff KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* totalRelMCPIDEff_pipimumu = new TH1D("totalRelMCPIDEff_pipimumu","rel PID Eff pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);

  for(int i=0; i<rangesKK_low.size();++i){
    
    double PIDrelEff,dPIDrelEff;
    PIDrelEff=muonPID_KK->GetBinContent(i+1) * hadronPID_KK->GetBinContent(i+1);
    cout<<muonPID_KK->GetBinContent(i+1)<<"  "<<hadronPID_KK->GetBinContent(i+1)<<"  "<<PIDrelEff<<endl;
 
    dPIDrelEff=PIDrelEff*TMath::Sqrt( TMath::Power(muonPID_KK->GetBinError(i+1)/muonPID_KK->GetBinContent(i+1),2)+  TMath::Power(dEffMu_sys/muonPID_KK->GetBinContent(i+1),2) +
     				      TMath::Power(hadronPID_KK->GetBinError(i+1)/hadronPID_KK->GetBinContent(i+1),2)+  TMath::Power(dEffHad_sys/hadronPID_KK->GetBinContent(i+1),2) );

    //dPIDrelEff=TMath::Sqrt(TMath::Power(muonPID_KK->GetBinError(i+1),2) * TMath::Power(hadronPID_KK->GetBinContent(i+1),2) );
    totalRelPIDEff_KKmumu->SetBinContent(i+1,PIDrelEff);
    totalRelPIDEff_KKmumu->SetBinError(i+1,dPIDrelEff);
    //MC
    PIDrelEff=MCmuonPID_KK->GetBinContent(i+1) * MChadronPID_KK->GetBinContent(i+1);
    dPIDrelEff=TMath::Sqrt(TMath::Power(MCmuonPID_KK->GetBinError(i+1),2) * TMath::Power(MChadronPID_KK->GetBinContent(i+1),2) );
    totalRelMCPIDEff_KKmumu->SetBinContent(i+1,PIDrelEff);
    totalRelMCPIDEff_KKmumu->SetBinError(i+1,dPIDrelEff);

  }

  for(int i=0; i<rangespipi_low.size();++i){
    
    double PIDrelEff,dPIDrelEff;
    PIDrelEff=muonPID_pipi->GetBinContent(i+1) * hadronPID_pipi->GetBinContent(i+1);

    dPIDrelEff=PIDrelEff*TMath::Sqrt( TMath::Power(muonPID_pipi->GetBinError(i+1)/muonPID_pipi->GetBinContent(i+1),2)+  TMath::Power(dEffMu_sys/muonPID_pipi->GetBinContent(i+1),2) +
     				      TMath::Power(hadronPID_pipi->GetBinError(i+1)/hadronPID_pipi->GetBinContent(i+1),2)+  TMath::Power(dEffHad_sys/hadronPID_pipi->GetBinContent(i+1),2) );

    //std::cout<<TMath::Power(muonPID_pipi->GetBinError(i+1)/muonPID_pipi->GetBinContent(i+1),2)<<"  "<<TMath::Power(dEffMu_sys/muonPID_pipi->GetBinContent(i+1),2)<<"  "<<TMath::Power(hadronPID_pipi->GetBinError(i+1)/hadronPID_pipi->GetBinContent(i+1),2)<<"  "<<TMath::Power(dEffHad_sys/hadronPID_pipi->GetBinContent(i+1),2)<<endl;
    //dPIDrelEff=TMath::Sqrt(TMath::Power(muonPID_pipi->GetBinError(i+1),2) * TMath::Power(hadronPID_pipi->GetBinContent(i+1),2) );
    totalRelPIDEff_pipimumu->SetBinContent(i+1,PIDrelEff);
    totalRelPIDEff_pipimumu->SetBinError(i+1,dPIDrelEff);
    //MC
    PIDrelEff=MCmuonPID_pipi->GetBinContent(i+1) * MChadronPID_pipi->GetBinContent(i+1);
    dPIDrelEff=TMath::Sqrt(TMath::Power(MCmuonPID_pipi->GetBinError(i+1),2) * TMath::Power(MChadronPID_pipi->GetBinContent(i+1),2) );
    totalRelMCPIDEff_pipimumu->SetBinContent(i+1,PIDrelEff);
    totalRelMCPIDEff_pipimumu->SetBinError(i+1,dPIDrelEff);
  }
  
  dcastyle();

  TCanvas *c = new TCanvas("c","c");
  c->Divide(1,2);
  c->cd(1);
  totalRelPIDEff_KKmumu->GetYaxis()->SetRangeUser(0.8,1.0);
  totalRelPIDEff_KKmumu->GetYaxis()->SetTitle("#epsilon(D->KK#mu#mu)/#epsilon(D->K#pi#mu#mu)");
  totalRelPIDEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  totalRelPIDEff_KKmumu->Draw();
  totalRelMCPIDEff_KKmumu->SetLineColor(kCyan);
  totalRelMCPIDEff_KKmumu->Draw("SAME");

  c->cd(2);
  totalRelPIDEff_pipimumu->GetYaxis()->SetRangeUser(.95,1.3);
  totalRelPIDEff_pipimumu->GetYaxis()->SetTitle("#epsilon(D->#pi#pi#mu#mu)/#epsilon(D->K#pi#mu#mu)");
  totalRelPIDEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  totalRelPIDEff_pipimumu->Draw();
  totalRelMCPIDEff_pipimumu->SetLineColor(kCyan);
  totalRelMCPIDEff_pipimumu->Draw("SAME");

  totalRelPIDEff_KKmumu->Write();
  totalRelPIDEff_pipimumu->Write();
  
  totalRelMCPIDEff_KKmumu->Write();
  totalRelMCPIDEff_pipimumu->Write();

  c->Write();
  if(binningScheme=="default")
    c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/totalPIDefficiency.eps");

  fOut->Write();
  fOut->Close();

  //    c->cd(3);
  //TH1F* ratio = DrawRatio(totalRelMCPIDEff_pipimumu,totalRelPIDEff_pipimumu);
  //ratio->Draw();

}




void drawBothPolaritiesPIDCalibDoubleMuonEfficiency(){

  dcastyle();

  TFile *fUp = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_and_PIDCalib_magUp.root","READ");
  TFile *fDw = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_and_PIDCalib_magDw.root","READ");

  TH1* relEff_KKmumu_up =(TH1D*)fUp->Get("relEff_KKmumu"); 
  TH1* relMCEff_KKmumu_up= (TH1D*)fUp->Get("relMCEff_KKmumu");
  TH1* MC_PIDCalib_EffRatio_KKmumu_up=(TH1D*)fUp->Get("MC_PIDCalib_EffRatio_KKmumu");
  TH1* relEff_pipimumu_up=(TH1D*)fUp->Get("relEff_pipimumu");
  TH1* relMCEff_pipimumu_up=(TH1D*)fUp->Get("relMCEff_pipimumu");
  TH1* MC_PIDCalib_EffRatio_pipimumu_up=(TH1D*)fUp->Get("MC_PIDCalib_EffRatio_pipimumu");

  TH1* relEff_KKmumu_dw =(TH1D*)fDw->Get("relEff_KKmumu"); 
  TH1* relMCEff_KKmumu_dw= (TH1D*)fDw->Get("relMCEff_KKmumu");
  TH1* MC_PIDCalib_EffRatio_KKmumu_dw=(TH1D*)fDw->Get("MC_PIDCalib_EffRatio_KKmumu");
  TH1* relEff_pipimumu_dw=(TH1D*)fDw->Get("relEff_pipimumu");
  TH1* relMCEff_pipimumu_dw=(TH1D*)fDw->Get("relMCEff_pipimumu");
  TH1* MC_PIDCalib_EffRatio_pipimumu_dw=(TH1D*)fDw->Get("MC_PIDCalib_EffRatio_pipimumu");
  
  TH1* relEff_KKmumu_average = new TH1D("relEff_KKmumu_average","relEff_KKmumu in bins of dimuon mass averaged both polarities",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relEff_pipimumu_average = new TH1D("relEff_pipimumu_average","relEff_pipimumu in bins of dimuon mass averaged both polarities",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1* relMCEff_KKmumu_average = new TH1D("relMCEff_KKmumu_average","relEff_KKmumu in bins of dimuon mass averaged both polarities",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_pipimumu_average = new TH1D("relMCEff_pipimumu_average","relEff_pipimumu in bins of dimuon mass averaged both polarities",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TFile *fout = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies//PIDEfficiencydoubleMuon_PID_MC_and_PIDCalib_bothPolarities.root","RECREATE");

  TCanvas *a = new TCanvas("a","a");

  double average=0;
  double dAverage=0;

  a->Divide(2,2);
  a->cd(1);
  relEff_KKmumu_up->Draw();
  relEff_KKmumu_dw->SetLineColor(kRed);
  relEff_KKmumu_dw->Draw("SAME");
  a->cd(2);
  relMCEff_KKmumu_up->GetYaxis()->SetRangeUser(0.9,1);
  relMCEff_KKmumu_up->Draw();
  relMCEff_KKmumu_dw->SetLineColor(kRed);
  relMCEff_KKmumu_dw->Draw("SAME");
  a->cd(3);
  relEff_pipimumu_up->Draw();
  relEff_pipimumu_dw->SetLineColor(kRed);
  relEff_pipimumu_dw->Draw("SAME");
  a->cd(4);
  relMCEff_pipimumu_up->GetYaxis()->SetRangeUser(1.02,1.14);
  relMCEff_pipimumu_up->Draw();
  relMCEff_pipimumu_dw->SetLineColor(kRed);
  relMCEff_pipimumu_dw->Draw("SAME");
  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/douleMuon_PID_MC_and_PIDCalib_bothPolarities1.eps");

  TCanvas *b = new TCanvas("b","b");
  b->Divide(1,2);
  b->cd(1);
  MC_PIDCalib_EffRatio_KKmumu_up->Draw();
  MC_PIDCalib_EffRatio_KKmumu_dw->SetLineColor(kRed);
  MC_PIDCalib_EffRatio_KKmumu_dw->Draw("SAME");
  b->cd(2);
  MC_PIDCalib_EffRatio_pipimumu_up->Draw();
  MC_PIDCalib_EffRatio_pipimumu_dw->SetLineColor(kRed);
  MC_PIDCalib_EffRatio_pipimumu_dw->Draw("SAME");
  b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_and_PIDCalib_bothPolarities2.eps");

  //average over polarities
  for(int i=0; i<rangesKK_low.size();++i){
    average+=relEff_KKmumu_up->GetBinContent(i+1);
    dAverage=+relEff_KKmumu_up->GetBinError(i+1)*relEff_KKmumu_up->GetBinError(i+1)/2;
    average+=relEff_KKmumu_dw->GetBinContent(i+1);
    dAverage=+relEff_KKmumu_dw->GetBinError(i+1)*relEff_KKmumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                dAverage = TMath::Sqrt(dAverage);
    relEff_KKmumu_average->SetBinContent(i+1,average);
    relEff_KKmumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
    //MC
    average+=relMCEff_KKmumu_up->GetBinContent(i+1);
    dAverage=+relMCEff_KKmumu_up->GetBinError(i+1)*relMCEff_KKmumu_up->GetBinError(i+1)/2;
    average+=relMCEff_KKmumu_dw->GetBinContent(i+1);
    dAverage=+relMCEff_KKmumu_dw->GetBinError(i+1)*relMCEff_KKmumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                dAverage = TMath::Sqrt(dAverage);
    relMCEff_KKmumu_average->SetBinContent(i+1,average);
    relMCEff_KKmumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;

  }    

  for(int i=0; i<rangespipi_low.size();++i){
    average+=relEff_pipimumu_up->GetBinContent(i+1);
    dAverage=+relEff_pipimumu_up->GetBinError(i+1)*relEff_pipimumu_up->GetBinError(i+1)/2;
    average+=relEff_pipimumu_dw->GetBinContent(i+1);
    dAverage=+relEff_pipimumu_dw->GetBinError(i+1)*relEff_pipimumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                dAverage = TMath::Sqrt(dAverage);
    relEff_pipimumu_average->SetBinContent(i+1,average);
    relEff_pipimumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
    //MC
    average+=relMCEff_pipimumu_up->GetBinContent(i+1);
    dAverage=+relMCEff_pipimumu_up->GetBinError(i+1)*relMCEff_pipimumu_up->GetBinError(i+1)/2;
    average+=relMCEff_pipimumu_dw->GetBinContent(i+1);
    dAverage=+relMCEff_pipimumu_dw->GetBinError(i+1)*relMCEff_pipimumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                dAverage = TMath::Sqrt(dAverage);
    relMCEff_pipimumu_average->SetBinContent(i+1,average);
    relMCEff_pipimumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;

  }    

  TCanvas *c = new TCanvas("c","c");
  c->Divide(1,2);
  c->cd(1);
  relEff_KKmumu_average->GetYaxis()->SetRangeUser(0.8,1.);
  relEff_KKmumu_average->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_KKmumu_average->Draw(); 
  relEff_KKmumu_average->Write(); 
  relMCEff_KKmumu_average->SetLineColor(kCyan); 
  relMCEff_KKmumu_average->Draw("SAME");
  relMCEff_KKmumu_average->Write();
  c->cd(2);
  relEff_pipimumu_average->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_pipimumu_average->GetYaxis()->SetRangeUser(.9,1.2);
  relEff_pipimumu_average->Draw();
  relEff_pipimumu_average->Write();
  relMCEff_pipimumu_average->SetLineColor(kCyan); 
  relMCEff_pipimumu_average->Draw("SAME"); 
  relMCEff_pipimumu_average->Write();
  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_and_PIDCalib_bothPolarities3.eps");
  fout->Write();
  fout->Close();
}

void drawBothPolaritiesPIDCalibSingleMuonEfficiency(int muonIndex, TString binningScheme="default"){

  dcastyle();

  TString fileUp;
  TString fileDw;

  TFile *fUp = new TFile(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_mu%i_",muonIndex)+binningScheme+"_magUp.root","READ");
  TFile *fDw = new TFile(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_mu%i_",muonIndex)+binningScheme+"_magDw.root","READ");

  TH1* relEff_KKmumu_up =(TH1D*)fUp->Get("relEff_KKmumu"); 
  TH1* relMCEff_KKmumu_up= (TH1D*)fUp->Get("relMCEff_KKmumu");
  TH1* MC_PIDCalib_EffRatio_KKmumu_up=(TH1D*)fUp->Get("MC_PIDCalib_EffRatio_KKmumu");
  TH1* relEff_pipimumu_up=(TH1D*)fUp->Get("relEff_pipimumu");
  TH1* relMCEff_pipimumu_up=(TH1D*)fUp->Get("relMCEff_pipimumu");
  TH1* MC_PIDCalib_EffRatio_pipimumu_up=(TH1D*)fUp->Get("MC_PIDCalib_EffRatio_pipimumu");

  TH1* relEff_KKmumu_dw =(TH1D*)fDw->Get("relEff_KKmumu"); 
  TH1* relMCEff_KKmumu_dw= (TH1D*)fDw->Get("relMCEff_KKmumu");
  TH1* MC_PIDCalib_EffRatio_KKmumu_dw=(TH1D*)fDw->Get("MC_PIDCalib_EffRatio_KKmumu");
  TH1* relEff_pipimumu_dw=(TH1D*)fDw->Get("relEff_pipimumu");
  TH1* relMCEff_pipimumu_dw=(TH1D*)fDw->Get("relMCEff_pipimumu");
  TH1* MC_PIDCalib_EffRatio_pipimumu_dw=(TH1D*)fDw->Get("MC_PIDCalib_EffRatio_pipimumu");
  
  TH1* relEff_KKmumu_average = new TH1D("relEff_KKmumu_average","relEff_KKmumu in bins of dimuon mass averaged both polarities",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relEff_pipimumu_average = new TH1D("relEff_pipimumu_average","relEff_pipimumu in bins of dimuon mass averaged both polarities",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1* relMCEff_KKmumu_average = new TH1D("relMCEff_KKmumu_average","relEff_KKmumu in bins of dimuon mass averaged both polarities",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_pipimumu_average = new TH1D("relMCEff_pipimumu_average","relEff_pipimumu in bins of dimuon mass averaged both polarities",sizeof(binspipi)/sizeof(double)-1,binspipi);

  //absolute efficiencies

  TH1* AbsEff_KKmumu_up =(TH1D*)fUp->Get("AbsEff_KKmumu"); 
  TH1* AbsMCEff_KKmumu_up= (TH1D*)fUp->Get("AbsMCEff_KKmumu");
  TH1* AbsEff_pipimumu_up=(TH1D*)fUp->Get("AbsEff_pipimumu");
  TH1* AbsMCEff_pipimumu_up=(TH1D*)fUp->Get("AbsMCEff_pipimumu");
  TH1* AbsEff_Kpimumu_up=(TH1D*)fUp->Get("AbsEff_Kpimumu");
  TH1* AbsMCEff_Kpimumu_up=(TH1D*)fUp->Get("AbsMCEff_Kpimumu");

  TH1* AbsEff_KKmumu_dw =(TH1D*)fDw->Get("AbsEff_KKmumu"); 
  TH1* AbsMCEff_KKmumu_dw= (TH1D*)fDw->Get("AbsMCEff_KKmumu");
  TH1* AbsEff_pipimumu_dw=(TH1D*)fDw->Get("AbsEff_pipimumu");
  TH1* AbsMCEff_pipimumu_dw=(TH1D*)fDw->Get("AbsMCEff_pipimumu");
  TH1* AbsEff_Kpimumu_dw=(TH1D*)fDw->Get("AbsEff_Kpimumu");
  TH1* AbsMCEff_Kpimumu_dw=(TH1D*)fDw->Get("AbsMCEff_Kpimumu");
  
  TH1* AbsEff_KKmumu_average = new TH1D("AbsEff_KKmumu_average","AbsEff_KKmumu in bins of dimuon mass averaged both polarities",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsEff_pipimumu_average = new TH1D("AbsEff_pipimumu_average","AbsEff_pipimumu in bins of dimuon mass averaged both polarities",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* AbsEff_Kpimumu_average = new TH1D("AbsEff_Kpimumu_average","AbsEff_Kpimumu in bins of dimuon mass averaged both polarities",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1* AbsMCEff_KKmumu_average = new TH1D("AbsMCEff_KKmumu_average","AbsEff_KKmumu in bins of dimuon mass averaged both polarities",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsMCEff_pipimumu_average = new TH1D("AbsMCEff_pipimumu_average","AbsEff_pipimumu in bins of dimuon mass averaged both polarities",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* AbsMCEff_Kpimumu_average = new TH1D("AbsMCEff_Kpimumu_average","AbsEff_Kpimumu in bins of dimuon mass averaged both polarities",sizeof(binsKpi)/sizeof(double)-1,binsKpi);


  TFile *fout = new TFile(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_bothPolarities_mu%i_",muonIndex)+binningScheme+".root","RECREATE");

  TCanvas *a = new TCanvas("a","a");

  double average=0;
  double dAverage=0;

  a->Divide(2,2);
  a->cd(1);
  relEff_KKmumu_up->Draw();
  relEff_KKmumu_dw->SetLineColor(kRed);
  relEff_KKmumu_dw->Draw("SAME");
  a->cd(2);
  relMCEff_KKmumu_up->GetYaxis()->SetRangeUser(0.9,1);
  relMCEff_KKmumu_up->Draw();
  relMCEff_KKmumu_dw->SetLineColor(kRed);
  relMCEff_KKmumu_dw->Draw("SAME");
  a->cd(3);
  relEff_pipimumu_up->Draw();
  relEff_pipimumu_dw->SetLineColor(kRed);
  relEff_pipimumu_dw->Draw("SAME");
  a->cd(4);
  relMCEff_pipimumu_up->GetYaxis()->SetRangeUser(1.02,1.14);
  relMCEff_pipimumu_up->Draw();
  relMCEff_pipimumu_dw->SetLineColor(kRed);
  relMCEff_pipimumu_dw->Draw("SAME");
  a->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon_PID_MC_and_PIDCalib_mu%i_bothPolarities1.eps",muonIndex));

  TCanvas *b = new TCanvas("b","b");
  b->Divide(1,2);
  b->cd(1);
  MC_PIDCalib_EffRatio_KKmumu_up->Draw();
  MC_PIDCalib_EffRatio_KKmumu_dw->SetLineColor(kRed);
  MC_PIDCalib_EffRatio_KKmumu_dw->Draw("SAME");
  b->cd(2);
  MC_PIDCalib_EffRatio_pipimumu_up->Draw();
  MC_PIDCalib_EffRatio_pipimumu_dw->SetLineColor(kRed);
  MC_PIDCalib_EffRatio_pipimumu_dw->Draw("SAME");
  if(binningScheme=="default")
    b->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon_PID_MC_and_PIDCalib_mu%i_bothPolarities2.eps",muonIndex));

  //average over polarities
  for(int i=0; i<rangesKK_low.size();++i){
    average+=relEff_KKmumu_up->GetBinContent(i+1);
    dAverage=+relEff_KKmumu_up->GetBinError(i+1)*relEff_KKmumu_up->GetBinError(i+1)/2;
    average+=relEff_KKmumu_dw->GetBinContent(i+1);
    dAverage=+relEff_KKmumu_dw->GetBinError(i+1)*relEff_KKmumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                dAverage = TMath::Sqrt(dAverage);
    relEff_KKmumu_average->SetBinContent(i+1,average);
    relEff_KKmumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
    //MC
    average+=relMCEff_KKmumu_up->GetBinContent(i+1);
    dAverage=+relMCEff_KKmumu_up->GetBinError(i+1)*relMCEff_KKmumu_up->GetBinError(i+1)/2;
    average+=relMCEff_KKmumu_dw->GetBinContent(i+1);
    dAverage=+relMCEff_KKmumu_dw->GetBinError(i+1)*relMCEff_KKmumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                dAverage = TMath::Sqrt(dAverage);
    relMCEff_KKmumu_average->SetBinContent(i+1,average);
    relMCEff_KKmumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;

    average+=AbsEff_KKmumu_up->GetBinContent(i+1);
    dAverage=+AbsEff_KKmumu_up->GetBinError(i+1)*AbsEff_KKmumu_up->GetBinError(i+1)/2;
    average+=AbsEff_KKmumu_dw->GetBinContent(i+1);
    dAverage=+AbsEff_KKmumu_dw->GetBinError(i+1)*AbsEff_KKmumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                dAverage = TMath::Sqrt(dAverage);
    AbsEff_KKmumu_average->SetBinContent(i+1,average);
    AbsEff_KKmumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
    //MC
    average+=AbsMCEff_KKmumu_up->GetBinContent(i+1);
    dAverage=+AbsMCEff_KKmumu_up->GetBinError(i+1)*AbsMCEff_KKmumu_up->GetBinError(i+1)/2;
    average+=AbsMCEff_KKmumu_dw->GetBinContent(i+1);
    dAverage=+AbsMCEff_KKmumu_dw->GetBinError(i+1)*AbsMCEff_KKmumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                dAverage = TMath::Sqrt(dAverage);
    AbsMCEff_KKmumu_average->SetBinContent(i+1,average);
    AbsMCEff_KKmumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
    
  }    

  for(int i=0; i<rangespipi_low.size();++i){
    average+=relEff_pipimumu_up->GetBinContent(i+1);
    dAverage=+relEff_pipimumu_up->GetBinError(i+1)*relEff_pipimumu_up->GetBinError(i+1)/2;
    average+=relEff_pipimumu_dw->GetBinContent(i+1);
    dAverage=+relEff_pipimumu_dw->GetBinError(i+1)*relEff_pipimumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                
    dAverage = TMath::Sqrt(dAverage);
    relEff_pipimumu_average->SetBinContent(i+1,average);
    relEff_pipimumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
    //MC
    average+=relMCEff_pipimumu_up->GetBinContent(i+1);
    dAverage=+relMCEff_pipimumu_up->GetBinError(i+1)*relMCEff_pipimumu_up->GetBinError(i+1)/2;
    average+=relMCEff_pipimumu_dw->GetBinContent(i+1);
    dAverage=+relMCEff_pipimumu_dw->GetBinError(i+1)*relMCEff_pipimumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                
    dAverage = TMath::Sqrt(dAverage);
    relMCEff_pipimumu_average->SetBinContent(i+1,average);
    relMCEff_pipimumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;

    average+=AbsEff_pipimumu_up->GetBinContent(i+1);
    dAverage=+AbsEff_pipimumu_up->GetBinError(i+1)*AbsEff_pipimumu_up->GetBinError(i+1)/2;
    average+=AbsEff_pipimumu_dw->GetBinContent(i+1);
    dAverage=+AbsEff_pipimumu_dw->GetBinError(i+1)*AbsEff_pipimumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                
    dAverage = TMath::Sqrt(dAverage);
    AbsEff_pipimumu_average->SetBinContent(i+1,average);
    AbsEff_pipimumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;
    //MC
    average+=AbsMCEff_pipimumu_up->GetBinContent(i+1);
    dAverage=+AbsMCEff_pipimumu_up->GetBinError(i+1)*AbsMCEff_pipimumu_up->GetBinError(i+1)/2;
    average+=AbsMCEff_pipimumu_dw->GetBinContent(i+1);
    dAverage=+AbsMCEff_pipimumu_dw->GetBinError(i+1)*AbsMCEff_pipimumu_dw->GetBinError(i+1)/2;
    average/=2;                                                                                                                                                                
    dAverage = TMath::Sqrt(dAverage);
    AbsMCEff_pipimumu_average->SetBinContent(i+1,average);
    AbsMCEff_pipimumu_average->SetBinError(i+1,dAverage);
    average=0;
    dAverage=0;

  }    

  //average Kpimumu polarities
  average+=AbsEff_Kpimumu_up->GetBinContent(1);
  dAverage=+AbsEff_Kpimumu_up->GetBinError(1)*AbsEff_Kpimumu_up->GetBinError(1)/2;
  average+=AbsEff_Kpimumu_dw->GetBinContent(1);
  dAverage=+AbsEff_Kpimumu_dw->GetBinError(1)*AbsEff_Kpimumu_dw->GetBinError(1)/2;
  average/=2;                                                                                                                                                                
  dAverage = TMath::Sqrt(dAverage);
  AbsEff_Kpimumu_average->SetBinContent(1,average);
  AbsEff_Kpimumu_average->SetBinError(1,dAverage);
  AbsEff_Kpimumu_average->Write();
  average=0;
  dAverage=0;
  //MC
  average+=AbsMCEff_Kpimumu_up->GetBinContent(1);
  dAverage=+AbsMCEff_Kpimumu_up->GetBinError(1)*AbsMCEff_Kpimumu_up->GetBinError(1)/2;
  average+=AbsMCEff_Kpimumu_dw->GetBinContent(1);
  dAverage=+AbsMCEff_Kpimumu_dw->GetBinError(1)*AbsMCEff_Kpimumu_dw->GetBinError(1)/2;
  average/=2;                                                                                                                                                                
  dAverage = TMath::Sqrt(dAverage);
  AbsMCEff_Kpimumu_average->SetBinContent(1,average);
  AbsMCEff_Kpimumu_average->SetBinError(1,dAverage);
  average=0;
  dAverage=0;
  AbsMCEff_Kpimumu_average->Write();

  
  TCanvas *c = new TCanvas("c","c");
  c->Divide(1,2);
  c->cd(1);
  relEff_KKmumu_average->GetYaxis()->SetRangeUser(0.9,1.);
  relEff_KKmumu_average->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_KKmumu_average->Draw(); 
  relEff_KKmumu_average->Write(); 
  relMCEff_KKmumu_average->SetLineColor(kCyan); 
  relMCEff_KKmumu_average->Draw("SAME");
  relMCEff_KKmumu_average->Write();
  c->cd(2);
  relEff_pipimumu_average->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_pipimumu_average->GetYaxis()->SetRangeUser(1.0,1.1);
  relEff_pipimumu_average->Draw();
  relEff_pipimumu_average->Write();
  relMCEff_pipimumu_average->SetLineColor(kCyan); 
  relMCEff_pipimumu_average->Draw("SAME"); 
  relMCEff_pipimumu_average->Write();
  if(binningScheme=="default")
    c->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon_PID_MC_and_PIDCalib_mu%i_bothPolarities3.eps",muonIndex));


  TCanvas *d = new TCanvas("d","d");
  d->Divide(1,2);
  d->cd(1);
  //AbsEff_KKmumu_average->GetYaxis()->SetRangeUser(0.9,1.);
  AbsEff_KKmumu_average->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsEff_KKmumu_average->Draw(); 
  AbsEff_KKmumu_average->Write(); 
  AbsMCEff_KKmumu_average->SetLineColor(kCyan); 
  AbsMCEff_KKmumu_average->Draw("SAME");
  AbsMCEff_KKmumu_average->Write();
  d->cd(2);
  AbsEff_pipimumu_average->GetXaxis()->SetTitle("m(#mu#mu)");
  //AbsEff_pipimumu_average->GetYaxis()->SetRangeUser(1.0,1.1);
  AbsEff_pipimumu_average->Draw();
  AbsEff_pipimumu_average->Write();
  AbsMCEff_pipimumu_average->SetLineColor(kCyan); 
  AbsMCEff_pipimumu_average->Draw("SAME"); 
  AbsMCEff_pipimumu_average->Write();
  if(binningScheme=="default")
    d->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon_PID_MC_and_PIDCalib_bothPolarities4.eps");

  TCanvas *e = new TCanvas("e","e");
  e->Divide(2,2);
  e->cd(1);
  AbsEff_KKmumu_up->Draw();
  AbsEff_KKmumu_dw->SetLineColor(kRed);
  AbsEff_KKmumu_dw->Draw("SAME");
  e->cd(2);
  //  AbsMCEff_KKmumu_up->GetYaxis()->SetRangeUser(0.9,1);
  AbsMCEff_KKmumu_up->Draw();
  AbsMCEff_KKmumu_dw->SetLineColor(kRed);
  AbsMCEff_KKmumu_dw->Draw("SAME");
  e->cd(3);
  AbsEff_pipimumu_up->Draw();
  AbsEff_pipimumu_dw->SetLineColor(kRed);
  AbsEff_pipimumu_dw->Draw("SAME");
  e->cd(4);
  //AbsMCEff_pipimumu_up->GetYaxis()->SetRangeUser(1.02,1.14);
  AbsMCEff_pipimumu_up->Draw();
  AbsMCEff_pipimumu_dw->SetLineColor(kRed);
  AbsMCEff_pipimumu_dw->Draw("SAME");
  if(binningScheme=="default")
    e->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon_PID_MC_and_PIDCalib_mu%i_bothPolarities4.eps",muonIndex));
  
  fout->Write();
  fout->Close();
}
  
 
void  evaluatePIDCalibMuonEfficiencyForPTBins(TString polarity) {

  dcastyle();
  
  TH1* AbsEff_Kpimumu = new TH1D("AbsEff_Kpimumu","AbsEff_Kpimumu in bins of Pt",sizeof(binsPt)/sizeof(double)-1,binsPt);
  TH1* AbsMCEff_Kpimumu = new TH1D("AbsMCEff_Kpimumu","Abs MC Eff_Kpimumu in bins of Pt",sizeof(binsPt)/sizeof(double)-1,binsPt);


  std::vector<TH1D*> AbsEff_KKmumu;
  std::vector<TH1D*> relEff_KKmumu;
  std::vector<TH1D*> AbsMCEff_KKmumu;
  std::vector<TH1D*> relMCEff_KKmumu;
  TH1D* AbsEff_KKmumu_temp;
  TH1D* relEff_KKmumu_temp;

  std::vector<TH1D*> AbsEff_pipimumu;
  std::vector<TH1D*> relEff_pipimumu;
  std::vector<TH1D*> AbsMCEff_pipimumu;
  std::vector<TH1D*> relMCEff_pipimumu;

  std::vector<TH1D*> MC_PIDCalib_EffRatio_pipimumu;
  std::vector<TH1D*> MC_PIDCalib_EffRatio_KKmumu;

  std::vector<TH1D*> MC_PIDCalib_AbsEffRatio_pipimumu;
  std::vector<TH1D*> MC_PIDCalib_AbsEffRatio_KKmumu;
  TH1D* MC_PIDCalib_AbsEffRatio_Kpimumu =  new TH1D(TString::Format("MC_data_PIDCalib_relAbsEff_Kpimumu_%i",0),TString::Format("ratio MC and PIDcalib eff ratios pipimumu m(mumu) bin %i",0),sizeof(binsPt)/sizeof(double)-1,binsPt);

  TF1 * f1= new TF1("f1","[0]+x*[1]",0,3000);
  f1->SetLineColor(kRed);
  TF1 * f2= new TF1("f2","[0]+x*[1]",0,3000);
  f2->SetLineColor(kRed);
  TF1 * f3= new TF1("f3","[0]+x*[1]",0,3000);
  f3->SetLineColor(kRed);


  TH1D * PIDCalib_extrapolatedEff_KKmumu = new TH1D("PIDCalib_extrapolatedEff_KKmumu","PIDCalib Eff extrapolated to low pt KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * PIDCalib_extrapolatedEff_pipimumu = new TH1D("PIDCalib_extrapolatedEff_pipimumu","PIDCalib Eff extrapolated to low pt pipimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * PIDCalib_extrapolatedEff_Kpimumu = new TH1D("PIDCalib_extrapolatedEff_Kpimumu","PIDCalib Eff extrapolated to low pt Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D* AbsEff_pipimumu_temp;
  TH1D* relEff_pipimumu_temp;
 
  double tempEff, tempEff_norm,tempEff_err, tempEff_norm_err;
  double tempMCEff, tempMCEff_norm,tempMCEff_norm_err,tempMCEff_err;
  double relMCEff, relMCEff_err,relEff_err;
  double ratioEffRatios,ratioEffRatios_err;

  TString MC_file;
  TString pathToPIDCalibFile="/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/resultsPtBins/"+polarity+"/";

  // std::vector<TH1D*> MC_PIDCalib_AbsEffRatio_pipimumu;
  //std::vector<TH1D*> MC_PIDCalib_AbsEffRatio_KKmumu;
  //std::vector<TH1D*> MC_PIDCalib_AbsEffRatio_Kpimumu;

  //perform the study in each bin of dimuon mass indivudially
  for(int i=0; i<rangesKK_low.size();++i){
    AbsEff_KKmumu_temp = new TH1D(TString::Format("AbsEff_KKmumu_%i",i),TString::Format("AbsEff_KKmumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    AbsEff_KKmumu.push_back(AbsEff_KKmumu_temp);
    relEff_KKmumu_temp = new TH1D(TString::Format("relEff_KKmumu_%i",i),TString::Format("relEff_KKmumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    relEff_KKmumu.push_back(relEff_KKmumu_temp);
    AbsEff_KKmumu_temp = new TH1D(TString::Format("AbsMCEff_KKmumu_%i",i),TString::Format("AbsMCEff_KKmumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    AbsMCEff_KKmumu.push_back(AbsEff_KKmumu_temp);
    relEff_KKmumu_temp = new TH1D(TString::Format("relMCEff_KKmumu_%i",i),TString::Format("relMCEff_KKmumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    relMCEff_KKmumu.push_back(relEff_KKmumu_temp);
    //histo for ratio of eff ratios data/MC
    relEff_KKmumu_temp = new TH1D(TString::Format("MC_data_PIDCalib_relEff_KKmumu_%i",i),TString::Format("ratio MC and PIDcalib eff ratios KKmumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    MC_PIDCalib_EffRatio_KKmumu.push_back(relEff_KKmumu_temp);
    relEff_KKmumu_temp = new TH1D(TString::Format("MC_data_PIDCalib_AbsEff_KKmumu_%i",i),TString::Format("ratio MC and PIDcalib eff ratios KKmumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    MC_PIDCalib_AbsEffRatio_KKmumu.push_back(relEff_KKmumu_temp);
    
  }    

  for(int i=0; i<rangespipi_low.size();++i){
    AbsEff_pipimumu_temp = new TH1D(TString::Format("AbsEff_pipimumu_%i",i),TString::Format("AbsEff_pipimumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    AbsEff_pipimumu.push_back(AbsEff_pipimumu_temp);
    relEff_pipimumu_temp = new TH1D(TString::Format("relEff_pipimumu_%i",i),TString::Format("relEff_pipimumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    relEff_pipimumu.push_back(relEff_pipimumu_temp);
    AbsEff_pipimumu_temp = new TH1D(TString::Format("AbsMCEff_pipimumu_%i",i),TString::Format("AbsMCEff_pipimumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    AbsMCEff_pipimumu.push_back(AbsEff_pipimumu_temp);
    relEff_pipimumu_temp = new TH1D(TString::Format("relMCEff_pipimumu_%i",i),TString::Format("relMCEff_pipimumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    relMCEff_pipimumu.push_back(relEff_pipimumu_temp);
    //histo for ratio of eff ratios data/MC
    relEff_pipimumu_temp = new TH1D(TString::Format("MC_data_PIDCalib_relEff_pipimumu_%i",i),TString::Format("ratio MC and PIDcalib eff ratios pipimumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    MC_PIDCalib_EffRatio_pipimumu.push_back(relEff_pipimumu_temp);
    relEff_pipimumu_temp = new TH1D(TString::Format("MC_data_PIDCalib_AbsEff_pipimumu_%i",i),TString::Format("ratio MC and PIDcalib eff ratios pipimumu m(mumu) bin %i",i),sizeof(binsPt)/sizeof(double)-1,binsPt);
    MC_PIDCalib_AbsEffRatio_pipimumu.push_back(relEff_pipimumu_temp);
  }    



  TString fIn ;

  /////////// LOOP OVER PT BINS//////////////////////

  for(int j=0; j<rangesPt_low.size();++j){                   

    // std::cout<<"j "<<j<<std::endl;

    ////////////////////////
    // normalization channel 
    ////////////////////////

    if(j!=0) {
      fIn = pathToPIDCalibFile+TString::Format("D2Kpimumu_%.1f_%.1f_%.1f_%.1f.root",rangesKpi_low[0],rangesKpi_high[0],rangesPt_low[j],rangesPt_high[j]);  
      tempEff_norm = evaluatePIDCalibEfficiency(fIn).first;
      tempEff_norm_err = evaluatePIDCalibEfficiency(fIn).second;
      AbsEff_Kpimumu->SetBinContent(j+1,tempEff_norm);
      AbsEff_Kpimumu->SetBinError(j+1,tempEff_norm_err);
    } 
    else std::cout<<"first bin not possible for data"<<std::endl;

    //MC
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBins/MC_D2Kpimumu_%.1f_%.1f_%.1f_%.1f_D2KKmumuBDT_",rangesKpi_low[0],rangesKpi_high[0],rangesPt_low[j],rangesPt_high[j])+polarity+".root"; //take MC file of ONE pt bin
    tempMCEff_norm = getMCEfficiency(MC_file,"BDT_Tree","mu0_ProbNNmu>0.5","D_DiMuon_Mass>0").first;
    tempMCEff_norm_err = getMCEfficiency(MC_file,"BDT_Tree","mu0_ProbNNmu>0.5","D_DiMuon_Mass>0").second;
    AbsMCEff_Kpimumu->SetBinContent(j+1,tempMCEff_norm);
    AbsMCEff_Kpimumu->SetBinError(j+1,tempMCEff_norm_err);
    //ratio PID calib - MC abolute eff
    if(j!=0){
      ratioEffRatios=AbsMCEff_Kpimumu->GetBinContent(j+1)/AbsEff_Kpimumu->GetBinContent(j+1);
      MC_PIDCalib_AbsEffRatio_Kpimumu->SetBinContent(j+1,ratioEffRatios);
      std::cout<<"ratioEffRatios "<<ratioEffRatios<<std::endl;
    }
    //std::cout<<"temp eff norm "<<"  "<<j<<"  "<<tempEff_norm<<"   "<<tempMCEff_norm<<"  "<<AbsEff_Kpimumu->GetBinContent(j+1)<<std::endl;
 
    ////////////////////////
    // loop over signal channels
    ////////////////////////

    //KKmumu mode
  for(int i=0; i<rangesKK_low.size();++i){
    if(j!=0) {  //avoid fitst bin pt<800 for data as no calibration data is available
      fIn = pathToPIDCalibFile+TString::Format("D2KKmumu_%.1f_%.1f_%.1f_%.1f.root",rangesKK_low[i],rangesKK_high[i],rangesPt_low[j],rangesPt_high[j]);
      tempEff=evaluatePIDCalibEfficiency(fIn).first;
      tempEff_err=evaluatePIDCalibEfficiency(fIn).second;
      AbsEff_KKmumu[i]->SetBinContent(j+1,tempEff);
      AbsEff_KKmumu[i]->SetBinError(j+1,tempEff_err);
      //std::cout<<"relative eff "<<j<<"  "<<i<<"  "<<tempEff<<"  "<<AbsEff_Kpimumu->GetBinContent(j+1)<<"  "<<tempEff/AbsEff_Kpimumu->GetBinContent(j+1)<<endl;
      relEff_KKmumu[i]->SetBinContent(j+1,tempEff/AbsEff_Kpimumu->GetBinContent(j+1));
      relEff_err=TMath::Sqrt(TMath::Power(tempEff_err/AbsEff_Kpimumu->GetBinContent(j+1),2) + TMath::Power(tempEff*AbsEff_Kpimumu->GetBinError(j+1)/(AbsEff_Kpimumu->GetBinContent(j+1)*AbsEff_Kpimumu->GetBinContent(j+1)),2));
      relEff_KKmumu[i]->SetBinError(j+1,relEff_err);
    }
     else std::cout<<"first bin not possible for data"<<std::endl;
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBins/MC_D2KKmumu_%.1f_%.1f_%.1f_%.1f_D2KKmumuBDT_",rangesKK_low[i],rangesKK_high[i],rangesPt_low[j],rangesPt_high[j])+polarity+".root";
    tempMCEff = getMCEfficiency(MC_file,"BDT_Tree","mu0_ProbNNmu>0.5","D_DiMuon_Mass>0").first;
    tempMCEff_err = getMCEfficiency(MC_file,"BDT_Tree","mu0_ProbNNmu>0.5","D_DiMuon_Mass>0").second;
    AbsMCEff_KKmumu[i]->SetBinContent(j+1,tempMCEff);
    AbsMCEff_KKmumu[i]->SetBinError(j+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(j+1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(j+1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(j+1)/(AbsMCEff_Kpimumu->GetBinContent(j+1)*AbsMCEff_Kpimumu->GetBinContent(j+1)),2) );
    relMCEff_KKmumu[i]->SetBinContent(j+1,relMCEff);
    relMCEff_KKmumu[i]->SetBinError(j+1,relMCEff_err);
    //divide MC and PID calib Efficiency Ratios!
    if(j!=0) {
      ratioEffRatios=relMCEff_KKmumu[i]->GetBinContent(j+1)/relEff_KKmumu[i]->GetBinContent(j+1);
      ratioEffRatios_err=TMath::Sqrt(TMath::Power(relMCEff_KKmumu[i]->GetBinError(j+1)/relEff_KKmumu[i]->GetBinContent(j+1),2) + TMath::Power(relMCEff_KKmumu[i]->GetBinContent(j+1)*relEff_KKmumu[i]->GetBinError(j+1)/(relEff_KKmumu[i]->GetBinContent(j+1)*relEff_KKmumu[i]->GetBinContent(j+1)),2));
      MC_PIDCalib_EffRatio_KKmumu[i]->SetBinContent(j+1,ratioEffRatios);
      MC_PIDCalib_EffRatio_KKmumu[i]->SetBinError(j+1,ratioEffRatios_err);
      //also for absolute eff!
      ratioEffRatios=AbsMCEff_KKmumu[i]->GetBinContent(j+1)/AbsEff_KKmumu[i]->GetBinContent(j+1);
      ratioEffRatios_err=TMath::Sqrt(TMath::Power(AbsMCEff_KKmumu[i]->GetBinError(j+1)/AbsEff_KKmumu[i]->GetBinContent(j+1),2) + TMath::Power(AbsMCEff_KKmumu[i]->GetBinContent(j+1)*AbsEff_KKmumu[i]->GetBinError(j+1)/(AbsEff_KKmumu[i]->GetBinContent(j+1)*AbsEff_KKmumu[i]->GetBinContent(j+1)),2));
      MC_PIDCalib_AbsEffRatio_KKmumu[i]->SetBinContent(j+1,ratioEffRatios);
      MC_PIDCalib_AbsEffRatio_KKmumu[i]->SetBinError(j+1,relMCEff_err);      
    }
   else std::cout<<"first bin not possible for data"<<std::endl;
    
  }

   //pipimumu mode
  for(int i=0; i<rangespipi_low.size();++i){
    if(j!=0) {  //avoid fitst bin pt<800 for data as no calibration data is available
      fIn = pathToPIDCalibFile+TString::Format("D2pipimumu_%.1f_%.1f_%.1f_%.1f.root",rangespipi_low[i],rangespipi_high[i],rangesPt_low[j],rangesPt_high[j]);
      tempEff=evaluatePIDCalibEfficiency(fIn).first;
      tempEff_err=evaluatePIDCalibEfficiency(fIn).second;
      AbsEff_pipimumu[i]->SetBinContent(j+1,tempEff);
      AbsEff_pipimumu[i]->SetBinError(j+1,tempEff_err);
      //std::cout<<"relative eff "<<j<<"  "<<i<<"  "<<tempEff<<"  "<<AbsEff_Kpimumu->GetBinContent(j+1)<<"  "<<tempEff/AbsEff_Kpimumu->GetBinContent(j+1)<<endl;
      relEff_pipimumu[i]->SetBinContent(j+1,tempEff/AbsEff_Kpimumu->GetBinContent(j+1));
      relEff_err=TMath::Sqrt(TMath::Power(tempEff_err/AbsEff_Kpimumu->GetBinContent(j+1),2) + TMath::Power(tempEff*AbsEff_Kpimumu->GetBinError(j+1)/(AbsEff_Kpimumu->GetBinContent(j+1)*AbsEff_Kpimumu->GetBinContent(j+1)),2));
      relEff_pipimumu[i]->SetBinError(j+1,relEff_err);
    }
    else std::cout<<"first bin not possible for data"<<std::endl;
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBins/MC_D2pipimumu_%.1f_%.1f_%.1f_%.1f_D2pipimumuBDT_",rangespipi_low[i],rangespipi_high[i],rangesPt_low[j],rangesPt_high[j])+polarity+".root";
    tempMCEff = getMCEfficiency(MC_file,"BDT_Tree","mu0_ProbNNmu>0.5","D_DiMuon_Mass>0").first;
    tempMCEff_err = getMCEfficiency(MC_file,"BDT_Tree","mu0_ProbNNmu>0.5","D_DiMuon_Mass>0").second;
    AbsMCEff_pipimumu[i]->SetBinContent(j+1,tempMCEff);
    AbsMCEff_pipimumu[i]->SetBinError(j+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(j+1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(j+1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(j+1)/(AbsMCEff_Kpimumu->GetBinContent(j+1)*AbsMCEff_Kpimumu->GetBinContent(j+1)),2) );
    relMCEff_pipimumu[i]->SetBinContent(j+1,relMCEff);
    relMCEff_pipimumu[i]->SetBinError(j+1,relMCEff_err);
    //divide MC and PID calib Efficiency Ratios!
    if(j!=0) {
      ratioEffRatios=relMCEff_pipimumu[i]->GetBinContent(j+1)/relEff_pipimumu[i]->GetBinContent(j+1);
      ratioEffRatios_err=TMath::Sqrt(TMath::Power(relMCEff_pipimumu[i]->GetBinError(j+1)/relEff_pipimumu[i]->GetBinContent(j+1),2) + TMath::Power(relMCEff_pipimumu[i]->GetBinContent(j+1)*relEff_pipimumu[i]->GetBinError(j+1)/(relEff_pipimumu[i]->GetBinContent(j+1)*relEff_pipimumu[i]->GetBinContent(j+1)),2));
      MC_PIDCalib_EffRatio_pipimumu[i]->SetBinContent(j+1,ratioEffRatios);
      MC_PIDCalib_EffRatio_pipimumu[i]->SetBinError(j+1,ratioEffRatios_err);
      //also for absolute eff!
ratioEffRatios=AbsMCEff_pipimumu[i]->GetBinContent(j+1)/AbsEff_pipimumu[i]->GetBinContent(j+1);
      ratioEffRatios_err=TMath::Sqrt(TMath::Power(AbsMCEff_pipimumu[i]->GetBinError(j+1)/AbsEff_pipimumu[i]->GetBinContent(j+1),2) + TMath::Power(AbsMCEff_pipimumu[i]->GetBinContent(j+1)*AbsEff_pipimumu[i]->GetBinError(j+1)/(AbsEff_pipimumu[i]->GetBinContent(j+1)*AbsEff_pipimumu[i]->GetBinContent(j+1)),2));
      MC_PIDCalib_AbsEffRatio_pipimumu[i]->SetBinContent(j+1,ratioEffRatios);
      MC_PIDCalib_AbsEffRatio_pipimumu[i]->SetBinError(j+1,relMCEff_err);      

   }
    else std::cout<<"first bin not possible for data"<<std::endl;
  }



 }


  //FIT and also extrapolate the Efficiency to low ptbin which is not covered by PID calib
 
  TFile *fout = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_"+polarity+"_ptBins.root","RECREATE");
  

  TCanvas *a = new TCanvas("a","a"); 
   
  a->Divide(2,2);
  for(int i=0; i<rangesKK_low.size();++i){
    a->cd(i+1);
  a->SetGridx();
  a->SetGridy();
    
  AbsEff_KKmumu[i]->GetYaxis()->SetRangeUser(0.6,1.1);
  AbsEff_KKmumu[i]->GetXaxis()->SetTitle("P_{t}(#mu)");
  AbsEff_KKmumu[i]->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}^{PIDCalib}");
  AbsEff_KKmumu[i]->Draw();
  AbsEff_KKmumu[i]->Write();
  AbsEff_KKmumu[i]->Fit("f1","R");
  //std::cout<<"extrapolated "<<f1->GetParameter(0)+f1->GetParameter(1)*400<<std::endl;
  PIDCalib_extrapolatedEff_KKmumu->SetBinContent(i+1,f1->GetParameter(0)+f1->GetParameter(1)*400);
  PIDCalib_extrapolatedEff_KKmumu->SetBinError(i+1,f1->GetParameter(1)*400);
  AbsMCEff_KKmumu[i]->SetLineColor(kCyan);
  AbsMCEff_KKmumu[i]->GetXaxis()->SetTitle("P_{t}(#mu)[MeV]");
  AbsMCEff_KKmumu[i]->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}^{MC}");
  AbsMCEff_KKmumu[i]->Draw("SAME");
  AbsMCEff_KKmumu[i]->Write();
  }

  a->cd(4);
  AbsEff_Kpimumu->GetYaxis()->SetRangeUser(0.6,1.1);
  AbsEff_Kpimumu->GetXaxis()->SetTitle("P_{t}(#mu)");
  AbsEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}^{PIDCalib}");
  AbsEff_Kpimumu->Draw();
  AbsEff_Kpimumu->Fit("f2","R");
  PIDCalib_extrapolatedEff_Kpimumu->SetBinContent(1,f2->GetParameter(0)+f2->GetParameter(1)*400);
  PIDCalib_extrapolatedEff_Kpimumu->SetBinError(1,f2->GetParameter(1)*400);
  AbsMCEff_Kpimumu->SetLineColor(kCyan);
  AbsEff_Kpimumu->GetXaxis()->SetTitle("P_{t}(#mu)[MeV]");
  AbsEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}^{MC}");
  AbsMCEff_Kpimumu->Draw("SAME");
  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PID_MC_and_PIDCalib_"+polarity+"_ptBins_1.eps");

  TCanvas *b = new TCanvas("b","b"); 
  b->Divide(2,2);
  for(int i=0; i<rangesKK_low.size();++i){
    b->cd(i+1);
    relEff_KKmumu[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    relEff_KKmumu[i]->GetXaxis()->SetTitle("P_{t}(#mu)[MeV]");
    relEff_KKmumu[i]->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}^{PIDCalib}/#epsilon_{D->K#pi#mu#mu}^{PIDCalib}");
    relEff_KKmumu[i]->Draw();
    relEff_KKmumu[i]->Write();
    relMCEff_KKmumu[i]->SetLineColor(kCyan);
    relMCEff_KKmumu[i]->GetXaxis()->SetTitle("P_{t}(#mu)[MeV]");
    relMCEff_KKmumu[i]->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}^{MC}/#epsilon_{D->K#pi#mu#mu}^{MC}");
    relMCEff_KKmumu[i]->Draw("SAME");
    relMCEff_KKmumu[i]->Write();
  }
  b->cd(4);
  //  relEff_KKmumu[0]->Draw();
  //relEff_KKmumu[1]->SetLineColor(kRed);
  //relEff_KKmumu[1]->Draw("SAME");
  //relEff_KKmumu[2]->SetLineColor(kBlack);
  //relEff_KKmumu[2]->Draw("SAME");
  
  b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PID_MC_and_PIDCalib_"+polarity+"_ptBins_2.eps");

  TCanvas *c = new TCanvas("c","c"); 
  c->Divide(3,2);
  for(int i=0; i<rangespipi_low.size();++i){
    c->cd(i+1);
    AbsEff_pipimumu[i]->GetYaxis()->SetRangeUser(0.5,1.1);
    AbsEff_pipimumu[i]->GetXaxis()->SetTitle("P_{t}(#mu)");
    AbsEff_pipimumu[i]->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}^{PIDCalib}");
    AbsEff_pipimumu[i]->Draw();
    AbsEff_pipimumu[i]->Write();
    AbsEff_pipimumu[i]->Fit("f3","R");
    PIDCalib_extrapolatedEff_pipimumu->SetBinContent(i+1,f3->GetParameter(0)+f3->GetParameter(1)*400);
    PIDCalib_extrapolatedEff_pipimumu->SetBinError(i+1,f3->GetParameter(1)*400);
    AbsMCEff_pipimumu[i]->SetLineColor(kCyan);
    AbsMCEff_pipimumu[i]->GetXaxis()->SetTitle("P_{t}(#mu)[MeV]");
    AbsMCEff_pipimumu[i]->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}^{MC}");
    AbsMCEff_pipimumu[i]->Draw("SAME");
    AbsMCEff_pipimumu[i]->Write();
  }
  c->cd(6);
  AbsEff_Kpimumu->GetYaxis()->SetRangeUser(0.5,1.1);
  //AbsEff_Kpimumu->Draw();
  AbsEff_Kpimumu->Draw("");
  AbsEff_Kpimumu->Write();
  AbsMCEff_Kpimumu->SetLineColor(kCyan);
  AbsMCEff_Kpimumu->Draw("SAME");
  AbsMCEff_Kpimumu->Write();

  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_"+polarity+"_ptBins_3.eps");

  TCanvas *d = new TCanvas("d","d"); 
  d->Divide(3,2);
  for(int i=0; i<rangespipi_low.size();++i){
    d->cd(i+1);
    relEff_pipimumu[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    relEff_pipimumu[i]->GetXaxis()->SetTitle("P_{t}(#mu)");
    relEff_pipimumu[i]->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}^{PIDCalib}/#epsilon_{D->K#pi#mu#mu}^{PIDCalib}");
    relEff_pipimumu[i]->Draw();
    relEff_pipimumu[i]->Write();
    relMCEff_pipimumu[i]->SetLineColor(kCyan);
    relEff_pipimumu[i]->GetXaxis()->SetTitle("P_{t}(#mu)[MeV]");
    relEff_pipimumu[i]->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}^{MC}/#epsilon_{D->K#pi#mu#mu}^{MC}");
    relMCEff_pipimumu[i]->Draw("SAME");
    relMCEff_pipimumu[i]->Write();
  }
  /*
  d->cd(6);
  relEff_pipimumu[0]->Draw();
  relEff_pipimumu[1]->SetLineColor(kRed);
  relEff_pipimumu[1]->Draw("SAME");
  relEff_pipimumu[2]->SetLineColor(kBlack);
  relEff_pipimumu[2]->Draw("SAME");
  relEff_pipimumu[3]->SetLineColor(kViolet);
  relEff_pipimumu[3]->Draw("SAME");
  relEff_pipimumu[4]->SetLineColor(kOrange);
  relEff_pipimumu[4]->Draw("SAME");
  */
  d->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_"+polarity+"_ptBins_4.eps");
  
  a->Write();
  b->Write();
  c->Write();
  d->Write();

 TCanvas *e = new TCanvas("e","e"); 

 e->Divide(2,2);
  for(int i=0; i<rangesKK_low.size();++i){
    e->cd(i+1); 
    AbsEff_KKmumu[i]->Divide(AbsMCEff_KKmumu[i]);
    AbsEff_KKmumu[i]->Draw();
  }

  e->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_."+polarity+"_ptBins_5.eps");
  e->Write();

 TCanvas *f = new TCanvas("f","f"); 

  f->Divide(2,2);
  for(int i=0; i<rangesKK_low.size();++i){
    f->cd(i+1); 
    MC_PIDCalib_EffRatio_KKmumu[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    MC_PIDCalib_EffRatio_KKmumu[i]->Draw();  
    MC_PIDCalib_EffRatio_KKmumu[i]->Write();
}

  f->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_."+polarity+"_ptBins_6.eps");
  f->Write();


 TCanvas *g = new TCanvas("g","g"); 

  g->Divide(3,2);
  for(int i=0; i<rangespipi_low.size();++i){
    g->cd(i+1); 
    MC_PIDCalib_EffRatio_pipimumu[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    MC_PIDCalib_EffRatio_pipimumu[i]->Draw();
    MC_PIDCalib_EffRatio_pipimumu[i]->Write();
  }

  g->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_."+polarity+"_ptBins_7.eps");
  g->Write();


  //    MC_PIDCalib_EffRatio_KKmumu_magDw[i]->Fit("f2","N");
  // average+= f2->GetParameter(0);
  //  dAverage=+ f2->GetParError(0)*f2->GetParError(0)/2;

 TCanvas *h = new TCanvas("h","h"); 

 h->Divide(2,2);
  h->cd(1);
  MC_PIDCalib_AbsEffRatio_Kpimumu->GetYaxis()->SetRangeUser(0.9,1.1);
  MC_PIDCalib_AbsEffRatio_Kpimumu->Draw();
  h->cd(2);
  MC_PIDCalib_AbsEffRatio_KKmumu[0]->GetYaxis()->SetRangeUser(0.9,1.1);
  MC_PIDCalib_AbsEffRatio_KKmumu[0]->Draw();
  h->cd(3);
  MC_PIDCalib_AbsEffRatio_pipimumu[0]->GetYaxis()->SetRangeUser(0.9,1.1);
  MC_PIDCalib_AbsEffRatio_pipimumu[0]->Draw();
  h->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_."+polarity+"_ptBins_8.eps");
  h->Write();


 TCanvas *i = new TCanvas("i","i"); 

  i->Divide(2,2);
  i->cd(1);
  PIDCalib_extrapolatedEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}^{extrapolalated}(p_{T}<800MeV)");
  PIDCalib_extrapolatedEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  PIDCalib_extrapolatedEff_KKmumu->Draw();
  i->cd(2);
  PIDCalib_extrapolatedEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}^{extrapolalated}(p_{T}<800MeV)");
  PIDCalib_extrapolatedEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  PIDCalib_extrapolatedEff_pipimumu->Draw();
  i->cd(3);
  PIDCalib_extrapolatedEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}^{extrapolalated}(p_{T}<800MeV)");
  PIDCalib_extrapolatedEff_Kpimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  PIDCalib_extrapolatedEff_Kpimumu->Draw();

  i->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PID_MC_and_PIDCalib_extrapolatedEff_"+polarity+"_ptBins.eps");

  PIDCalib_extrapolatedEff_KKmumu->Write();
  PIDCalib_extrapolatedEff_Kpimumu->Write();
  PIDCalib_extrapolatedEff_pipimumu->Write();
}


void  evaluatePIDCalibHadronEfficiency(TString polarity, TString binningScheme) {


  dcastyle();

  TString pathToPIDCalibFile;

  pathToPIDCalibFile = "/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/resultsHadronPID/"+polarity+"/";
  if(binningScheme=="alternative") pathToPIDCalibFile = "/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/resultsHadronPID/alternativeBinning/"+polarity+"/";
  if(binningScheme=="alternative2") pathToPIDCalibFile = "/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/resultsHadronPID/alternativeBinning2/"+polarity+"/";


  TFile *fout = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Hadron_PID_MC_and_PIDCalib_"+binningScheme+"_"+polarity+".root","RECREATE");
  
  TH1* AbsEff_Kpimumu = new TH1D("AbsEff_Kpimumu","AbsEff_Kpimumu in bins of dimuon mass",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1* AbsMCEff_Kpimumu = new TH1D("AbsMCEff_Kpimumu","Abs MC Eff_Kpimumu in bins of dimuon mass",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1* AbsEff_KKmumu = new TH1D("AbsEff_KKmumu","AbsEff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsMCEff_KKmumu = new TH1D("AbsMCEff_KKmumu","Abs MC Eff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relEff_KKmumu = new TH1D("relEff_KKmumu","relEff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_KKmumu = new TH1D("relMCEff_KKmumu","rel MC Eff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* MC_PIDCalib_EffRatio_KKmumu= new TH1D("MC_PIDCalib_EffRatio_KKmumu","data and MC ratio of relative signal/nomralization EffRatio",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1* AbsEff_pipimumu = new TH1D("AbsEff_pipimumu","AbsEff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* AbsMCEff_pipimumu = new TH1D("AbsMCEff_pipimumu","Abs MC Eff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relEff_pipimumu = new TH1D("relEff_pipimumu","relEff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relMCEff_pipimumu = new TH1D("relMCEff_pipimumu","rel MC Eff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* MC_PIDCalib_EffRatio_pipimumu= new TH1D("MC_PIDCalib_EffRatio_pipimumu","data and MC ratio of relative signal/nomralization EffRatio",sizeof(binspipi)/sizeof(double)-1,binspipi);

 
  double tempEff, tempEff_norm,tempEff_err, tempEff_norm_err;
  double tempMCEff, tempMCEff_norm,tempMCEff_norm_err,tempMCEff_err;
  double relMCEff, relMCEff_err,relEff_err;
  double ratioEffRatios,ratioEffRatios_err;

  TString MC_file;

  TString hadronPIDcut_Kpi="h0_ProbNNk>0.2&&h0_ProbNNghost<0.5&&h1_ProbNNpi>0.2&&h1_ProbNNghost<0.5&&Slowpi_ProbNNghost<0.5";
  TString hadronPIDcut_KK="h0_ProbNNk>0.2&&h0_ProbNNghost<0.5&&h1_ProbNNk>0.2&&h1_ProbNNghost<0.5&&Slowpi_ProbNNghost<0.5";
  TString hadronPIDcut_pipi="h0_ProbNNpi>0.2&&h0_ProbNNghost<0.5&&h1_ProbNNpi>0.2&&h1_ProbNNghost<0.5&&Slowpi_ProbNNghost<0.5";
  
  
  //norm mode
  for(int i=0; i<rangesKpi_low.size();++i){
    TString fIn = pathToPIDCalibFile+TString::Format("D2Kpimumu_%.1f_%.1f.root",rangesKpi_low[0],rangesKpi_high[0]);  
    tempEff_norm = evaluatePIDCalibEfficiency(fIn).first;
    tempEff_norm_err = evaluatePIDCalibEfficiency(fIn).second;
    AbsEff_Kpimumu->SetBinContent(i+1,tempEff_norm);
    AbsEff_Kpimumu->SetBinError(i+1,tempEff_norm_err);
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_",rangesKpi_low[0],rangesKpi_high[0])+polarity+".root";
    tempMCEff_norm = getMCEfficiency(MC_file,"BDT_Tree",hadronPIDcut_Kpi,"D_DiMuon_Mass>0").first;
    tempMCEff_norm_err = getMCEfficiency(MC_file,"BDT_Tree",hadronPIDcut_Kpi,"D_DiMuon_Mass>0").second;
    AbsMCEff_Kpimumu->SetBinContent(i+1,tempMCEff_norm);
    AbsMCEff_Kpimumu->SetBinError(i+1,tempMCEff_norm_err);
  }

  //KKmumu mode
  for(int i=0; i<rangesKK_low.size();++i){
    TString fIn = pathToPIDCalibFile+TString::Format("D2KKmumu_%.1f_%.1f.root",rangesKK_low[i],rangesKK_high[i]);  
    tempEff = evaluatePIDCalibEfficiency(fIn).first;
    tempEff_err = evaluatePIDCalibEfficiency(fIn).second;
    AbsEff_KKmumu->SetBinContent(i+1,tempEff);
    AbsEff_KKmumu->SetBinError(i+1,tempEff_err);
    relEff_KKmumu->SetBinContent(i+1,tempEff/AbsEff_Kpimumu->GetBinContent(1));
    std::cout<<rangesKK_low[i]<<"-"<<rangesKK_high[i]<<"  "<<tempEff<<"  "<<AbsEff_Kpimumu->GetBinContent(1)<<std::endl;
    relEff_err=TMath::Sqrt(TMath::Power(tempEff_err/AbsEff_Kpimumu->GetBinContent(1),2) + TMath::Power(tempEff*AbsEff_Kpimumu->GetBinError(1)/(AbsEff_Kpimumu->GetBinContent(1)*AbsEff_Kpimumu->GetBinContent(1)),2));
    relEff_KKmumu->SetBinError(i+1,relEff_err);
    //MC 
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_",rangesKK_low[i],rangesKK_high[i])+polarity+".root";
    tempMCEff = getMCEfficiency(MC_file,"BDT_Tree",hadronPIDcut_KK,"D_DiMuon_Mass>0").first;
    tempMCEff_err = getMCEfficiency(MC_file,"BDT_Tree",hadronPIDcut_KK,"D_DiMuon_Mass>0").second;
    AbsMCEff_KKmumu->SetBinContent(i+1,tempMCEff);
    AbsMCEff_KKmumu->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(1)/(AbsMCEff_Kpimumu->GetBinContent(1)*AbsMCEff_Kpimumu->GetBinContent(1)),2) );
    relMCEff_KKmumu->SetBinContent(i+1,relMCEff);
    relMCEff_KKmumu->SetBinError(i+1,relMCEff_err);
    //relative MC and PID calib
    ratioEffRatios=relMCEff_KKmumu->GetBinContent(i+1)/relEff_KKmumu->GetBinContent(i+1);
    ratioEffRatios_err=TMath::Sqrt(TMath::Power(relMCEff_KKmumu->GetBinError(i+1)/relEff_KKmumu->GetBinContent(i+1),2) + TMath::Power(relMCEff_KKmumu->GetBinContent(i+1)*relEff_KKmumu->GetBinError(i+1)/(relEff_KKmumu->GetBinContent(i+1)*relEff_KKmumu->GetBinContent(i+1)),2));
    MC_PIDCalib_EffRatio_KKmumu->SetBinContent(i+1,ratioEffRatios);
    MC_PIDCalib_EffRatio_KKmumu->SetBinError(i+1,ratioEffRatios_err);
  }

  //pipimumu mode
  for(int i=0; i<rangespipi_low.size();++i){
    TString fIn = pathToPIDCalibFile+TString::Format("D2pipimumu_%.1f_%.1f.root",rangespipi_low[i],rangespipi_high[i]);  
    tempEff = evaluatePIDCalibEfficiency(fIn).first;
    tempEff_err = evaluatePIDCalibEfficiency(fIn).second;
    AbsEff_pipimumu->SetBinContent(i+1,tempEff);
    AbsEff_pipimumu->SetBinError(i+1,tempEff_err);
    relEff_pipimumu->SetBinContent(i+1,tempEff/AbsEff_Kpimumu->GetBinContent(1));
    relEff_err=TMath::Sqrt(TMath::Power(tempEff_err/AbsEff_Kpimumu->GetBinContent(1),2) + TMath::Power(tempEff*AbsEff_Kpimumu->GetBinError(1)/(AbsEff_Kpimumu->GetBinContent(1)*AbsEff_Kpimumu->GetBinContent(1)),2));
    relEff_pipimumu->SetBinError(i+1,relEff_err);
    //MC 
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_",rangespipi_low[i],rangespipi_high[i])+polarity+".root";
    tempMCEff = getMCEfficiency(MC_file,"BDT_Tree",hadronPIDcut_pipi,"D_DiMuon_Mass>0").first;
    tempMCEff_err = getMCEfficiency(MC_file,"BDT_Tree",hadronPIDcut_pipi,"D_DiMuon_Mass>0").second;
    AbsMCEff_pipimumu->SetBinContent(i+1,tempMCEff);
    AbsMCEff_pipimumu->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(1)/(AbsMCEff_Kpimumu->GetBinContent(1)*AbsMCEff_Kpimumu->GetBinContent(1)),2) );
    relMCEff_pipimumu->SetBinContent(i+1,relMCEff);
    relMCEff_pipimumu->SetBinError(i+1,relMCEff_err);
    //relative MC and PID calib
    ratioEffRatios=relMCEff_pipimumu->GetBinContent(i+1)/relEff_pipimumu->GetBinContent(i+1);
    ratioEffRatios_err=TMath::Sqrt(TMath::Power(relMCEff_pipimumu->GetBinError(i+1)/relEff_pipimumu->GetBinContent(i+1),2) + TMath::Power(relMCEff_pipimumu->GetBinContent(i+1)*relEff_pipimumu->GetBinError(i+1)/(relEff_pipimumu->GetBinContent(i+1)*relEff_pipimumu->GetBinContent(i+1)),2));
    MC_PIDCalib_EffRatio_pipimumu->SetBinContent(i+1,ratioEffRatios);
    MC_PIDCalib_EffRatio_pipimumu->SetBinError(i+1,ratioEffRatios_err);
  }

 
  fout->cd();

  TCanvas *a = new TCanvas("a","a"); 
  a->Divide(1,2);
  a->cd(1);
  AbsEff_KKmumu->GetYaxis()->SetRangeUser(0.8,.9);
  AbsEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}");
  AbsEff_KKmumu->Draw();
  AbsEff_KKmumu->Write();
  AbsMCEff_KKmumu->SetLineColor(kCyan);
  AbsMCEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsMCEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}");
  AbsMCEff_KKmumu->Draw("SAME");
  AbsMCEff_KKmumu->Write();
  
  a->cd(2);
  relEff_KKmumu->GetYaxis()->SetRangeUser(0.9,1.0);
  relEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relEff_KKmumu->Draw();
  relEff_KKmumu->Write();
  relMCEff_KKmumu->SetLineColor(kCyan);
  relMCEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  relMCEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relMCEff_KKmumu->Draw("SAME");
  relMCEff_KKmumu->Write();
  if(binningScheme=="default")
    a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Hadron_PID_MC_and_PIDCalib_"+polarity+"_1.eps");

  
  TCanvas *b = new TCanvas("b","b"); 
 
  AbsEff_Kpimumu->GetYaxis()->SetRangeUser(0.85,.95);
  AbsEff_Kpimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}");
  AbsEff_Kpimumu->Draw();
  AbsEff_Kpimumu->Write();
  AbsMCEff_Kpimumu->SetLineColor(kCyan);
  AbsEff_Kpimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}");
  AbsMCEff_Kpimumu->Draw("SAME"); 
  AbsMCEff_Kpimumu->Write();
  if(binningScheme=="default")
    b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Hadron_PID_MC_and_PIDCalib_"+polarity+"_2.eps");


  TCanvas *c = new TCanvas("c","c"); 
  c->Divide(1,2);
  c->cd(1);
  AbsEff_pipimumu->GetYaxis()->SetRangeUser(0.9,1.);
  AbsEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}");
  AbsEff_pipimumu->Draw();
  AbsEff_pipimumu->Write();
  AbsMCEff_pipimumu->SetLineColor(kCyan);
  AbsMCEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsMCEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}");
  AbsMCEff_pipimumu->Draw("SAME");
  AbsMCEff_pipimumu->Write();
  
  c->cd(2);  
  relEff_pipimumu->GetYaxis()->SetRangeUser(1.02,1.14);
  relEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relEff_pipimumu->Draw();
  relEff_pipimumu->Write();
  relMCEff_pipimumu->SetLineColor(kCyan);
  relEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relMCEff_pipimumu->Draw("SAME");
  relMCEff_pipimumu->Write();
  c->Write();
  if(binningScheme=="default")
    c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Hadron_PID_MC_and_PIDCalib_"+polarity+"_3.eps");
 
  
  TCanvas *d = new TCanvas("d","d"); 
  d->Divide(1,2);
  
  d->cd(1); 
  MC_PIDCalib_EffRatio_KKmumu->GetYaxis()->SetRangeUser(0.95,1.04);
  MC_PIDCalib_EffRatio_KKmumu->Draw();  
  MC_PIDCalib_EffRatio_KKmumu->Write();

  d->cd(2); 
  MC_PIDCalib_EffRatio_pipimumu->GetYaxis()->SetRangeUser(0.95,1.04);
  MC_PIDCalib_EffRatio_pipimumu->Draw();
  MC_PIDCalib_EffRatio_pipimumu->Write();

  if(binningScheme=="default")
    d->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Hadron_PID_MC_and_PIDCalib_"+polarity+"_4.eps");
  d->Write();
  
  a->Write();
  b->Write();
  c->Write();

  fout->Write();
  fout->Close();

}

void  evaluatePIDCalibDoubleMuonEfficiency(TString polarity) {
  //both muons, integtated Pt

  dcastyle();

  TFile *fout = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_and_PIDCalib_"+polarity+".root","RECREATE");
  
  TH1* AbsEff_Kpimumu = new TH1D("AbsEff_Kpimumu","AbsEff_Kpimumu in bins of dimuon mass",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1* AbsMCEff_Kpimumu = new TH1D("AbsMCEff_Kpimumu","Abs MC Eff_Kpimumu in bins of dimuon mass",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1* AbsEff_KKmumu = new TH1D("AbsEff_KKmumu","AbsEff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsMCEff_KKmumu = new TH1D("AbsMCEff_KKmumu","Abs MC Eff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relEff_KKmumu = new TH1D("relEff_KKmumu","relEff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_KKmumu = new TH1D("relMCEff_KKmumu","rel MC Eff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* MC_PIDCalib_EffRatio_KKmumu= new TH1D("MC_PIDCalib_EffRatio_KKmumu","data and MC ratio of relative signal/nomralization EffRatio",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1* AbsEff_pipimumu = new TH1D("AbsEff_pipimumu","AbsEff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* AbsMCEff_pipimumu = new TH1D("AbsMCEff_pipimumu","Abs MC Eff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relEff_pipimumu = new TH1D("relEff_pipimumu","relEff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relMCEff_pipimumu = new TH1D("relMCEff_pipimumu","rel MC Eff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* MC_PIDCalib_EffRatio_pipimumu= new TH1D("MC_PIDCalib_EffRatio_pipimumu","data and MC ratio of relative signal/nomralization EffRatio",sizeof(binspipi)/sizeof(double)-1,binspipi);

 
  double tempEff, tempEff_norm,tempEff_err, tempEff_norm_err;
  double tempMCEff, tempMCEff_norm,tempMCEff_norm_err,tempMCEff_err;
  double relMCEff, relMCEff_err,relEff_err;
  double ratioEffRatios,ratioEffRatios_err;

  TString MC_file;
  TString pathToPIDCalibFile="/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/resultsDoubleMuonPID/"+polarity+"/";

  TString muonPIDCut="mu0_ProbNNmu>0.5&&mu0_ProbNNghost<0.5&&mu1_ProbNNmu>0.5&&mu1_ProbNNghost<0.5" ;
  
  
  //norm mode
  for(int i=0; i<rangesKpi_low.size();++i){
    TString fIn = pathToPIDCalibFile+TString::Format("D2Kpimumu_%.1f_%.1f_800.0_8000.0.root",rangesKpi_low[0],rangesKpi_high[0]);  
    tempEff_norm = evaluatePIDCalibEfficiency(fIn).first;
    tempEff_norm_err = evaluatePIDCalibEfficiency(fIn).second;
    AbsEff_Kpimumu->SetBinContent(i+1,tempEff_norm);
    AbsEff_Kpimumu->SetBinError(i+1,tempEff_norm_err);
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_",rangesKpi_low[0],rangesKpi_high[0])+polarity+".root";
    tempMCEff_norm = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"D_DiMuon_Mass>0").first;
    tempMCEff_norm_err = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"D_DiMuon_Mass>0").second;
    AbsMCEff_Kpimumu->SetBinContent(i+1,tempMCEff_norm);
    AbsMCEff_Kpimumu->SetBinError(i+1,tempMCEff_norm_err);
  }

  //KKmumu mode
  for(int i=0; i<rangesKK_low.size();++i){
    TString fIn = pathToPIDCalibFile+TString::Format("D2KKmumu_%.1f_%.1f_800.0_8000.0.root",rangesKK_low[i],rangesKK_high[i]);  
    tempEff = evaluatePIDCalibEfficiency(fIn).first;
    tempEff_err = evaluatePIDCalibEfficiency(fIn).second;
    AbsEff_KKmumu->SetBinContent(i+1,tempEff);
    AbsEff_KKmumu->SetBinError(i+1,tempEff_err);
    relEff_KKmumu->SetBinContent(i+1,tempEff/AbsEff_Kpimumu->GetBinContent(1));
    std::cout<<rangesKK_low[i]<<"-"<<rangesKK_high[i]<<"  "<<tempEff<<"  "<<AbsEff_Kpimumu->GetBinContent(1)<<std::endl;
    relEff_err=TMath::Sqrt(TMath::Power(tempEff_err/AbsEff_Kpimumu->GetBinContent(1),2) + TMath::Power(tempEff*AbsEff_Kpimumu->GetBinError(1)/(AbsEff_Kpimumu->GetBinContent(1)*AbsEff_Kpimumu->GetBinContent(1)),2));
    relEff_KKmumu->SetBinError(i+1,relEff_err);
    //MC 
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_",rangesKK_low[i],rangesKK_high[i])+polarity+".root";
    tempMCEff = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"D_DiMuon_Mass>0").first;
    tempMCEff_err = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"D_DiMuon_Mass>0").second;
    AbsMCEff_KKmumu->SetBinContent(i+1,tempMCEff);
    AbsMCEff_KKmumu->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(1)/(AbsMCEff_Kpimumu->GetBinContent(1)*AbsMCEff_Kpimumu->GetBinContent(1)),2) );
    relMCEff_KKmumu->SetBinContent(i+1,relMCEff);
    relMCEff_KKmumu->SetBinError(i+1,relMCEff_err);
    //relative MC and PID calib
    ratioEffRatios=relMCEff_KKmumu->GetBinContent(i+1)/relEff_KKmumu->GetBinContent(i+1);
    ratioEffRatios_err=TMath::Sqrt(TMath::Power(relMCEff_KKmumu->GetBinError(i+1)/relEff_KKmumu->GetBinContent(i+1),2) + TMath::Power(relMCEff_KKmumu->GetBinContent(i+1)*relEff_KKmumu->GetBinError(i+1)/(relEff_KKmumu->GetBinContent(i+1)*relEff_KKmumu->GetBinContent(i+1)),2));
    MC_PIDCalib_EffRatio_KKmumu->SetBinContent(i+1,ratioEffRatios);
    MC_PIDCalib_EffRatio_KKmumu->SetBinError(i+1,ratioEffRatios_err);
  }

  //pipimumu mode
  for(int i=0; i<rangespipi_low.size();++i){
    TString fIn = pathToPIDCalibFile+TString::Format("D2pipimumu_%.1f_%.1f_800.0_8000.0.root",rangespipi_low[i],rangespipi_high[i]);  
    tempEff = evaluatePIDCalibEfficiency(fIn).first;
    tempEff_err = evaluatePIDCalibEfficiency(fIn).second;
    AbsEff_pipimumu->SetBinContent(i+1,tempEff);
    AbsEff_pipimumu->SetBinError(i+1,tempEff_err);
    relEff_pipimumu->SetBinContent(i+1,tempEff/AbsEff_Kpimumu->GetBinContent(1));
    relEff_err=TMath::Sqrt(TMath::Power(tempEff_err/AbsEff_Kpimumu->GetBinContent(1),2) + TMath::Power(tempEff*AbsEff_Kpimumu->GetBinError(1)/(AbsEff_Kpimumu->GetBinContent(1)*AbsEff_Kpimumu->GetBinContent(1)),2));
    relEff_pipimumu->SetBinError(i+1,relEff_err);
    //MC 
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_",rangespipi_low[i],rangespipi_high[i])+polarity+".root";
    tempMCEff = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"D_DiMuon_Mass>0").first;
    tempMCEff_err = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"D_DiMuon_Mass>0").second;
    AbsMCEff_pipimumu->SetBinContent(i+1,tempMCEff);
    AbsMCEff_pipimumu->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(1)/(AbsMCEff_Kpimumu->GetBinContent(1)*AbsMCEff_Kpimumu->GetBinContent(1)),2) );
    relMCEff_pipimumu->SetBinContent(i+1,relMCEff);
    relMCEff_pipimumu->SetBinError(i+1,relMCEff_err);
    //relative MC and PID calib
    ratioEffRatios=relMCEff_pipimumu->GetBinContent(i+1)/relEff_pipimumu->GetBinContent(i+1);
    ratioEffRatios_err=TMath::Sqrt(TMath::Power(relMCEff_pipimumu->GetBinError(i+1)/relEff_pipimumu->GetBinContent(i+1),2) + TMath::Power(relMCEff_pipimumu->GetBinContent(i+1)*relEff_pipimumu->GetBinError(i+1)/(relEff_pipimumu->GetBinContent(i+1)*relEff_pipimumu->GetBinContent(i+1)),2));
    MC_PIDCalib_EffRatio_pipimumu->SetBinContent(i+1,ratioEffRatios);
    MC_PIDCalib_EffRatio_pipimumu->SetBinError(i+1,ratioEffRatios_err);
  }


 
  fout->cd();

  TCanvas *a = new TCanvas("a","a"); 
  a->Divide(1,2);
  a->cd(1);
  AbsEff_KKmumu->GetYaxis()->SetRangeUser(0.2,.8);
  AbsEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}");
  AbsEff_KKmumu->Draw();
  AbsEff_KKmumu->Write();
  AbsMCEff_KKmumu->SetLineColor(kCyan);
  AbsMCEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsMCEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}");
  AbsMCEff_KKmumu->Draw("SAME");
  AbsMCEff_KKmumu->Write();
  
  a->cd(2);
  relEff_KKmumu->GetYaxis()->SetRangeUser(0.85,1.1);
  relEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relEff_KKmumu->Draw();
  relEff_KKmumu->Write();
  relMCEff_KKmumu->SetLineColor(kCyan);
  relMCEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  relMCEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relMCEff_KKmumu->Draw("SAME");
  relMCEff_KKmumu->Write();
  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_and_PIDCalib_."+polarity+"_1.eps");

  
  TCanvas *b = new TCanvas("b","b"); 
 
  AbsEff_Kpimumu->GetYaxis()->SetRangeUser(0.2,.9);
  AbsEff_Kpimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}");
  AbsEff_Kpimumu->Draw();
  AbsEff_Kpimumu->Write();
  AbsMCEff_Kpimumu->SetLineColor(kCyan);
  AbsEff_Kpimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}");
  AbsMCEff_Kpimumu->Draw("SAME"); 
  AbsMCEff_Kpimumu->Write();

  b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_and_PIDCalib_."+polarity+"_2.eps");


  TCanvas *c = new TCanvas("c","c"); 
  c->Divide(1,2);
  c->cd(1);
  AbsEff_pipimumu->GetYaxis()->SetRangeUser(0.2,.8);
  AbsEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}");
  AbsEff_pipimumu->Draw();
  AbsEff_pipimumu->Write();
  AbsMCEff_pipimumu->SetLineColor(kCyan);
  AbsMCEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsMCEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}");
  AbsMCEff_pipimumu->Draw("SAME");
  AbsMCEff_pipimumu->Write();
  
  c->cd(2);  
  relEff_pipimumu->GetYaxis()->SetRangeUser(.85,1.2);
  relEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relEff_pipimumu->Draw();
  relEff_pipimumu->Write();
  relMCEff_pipimumu->SetLineColor(kCyan);
  relEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  relEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relMCEff_pipimumu->Draw("SAME");
  relMCEff_pipimumu->Write();
  c->Write();
  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_and_PIDCalib_."+polarity+"_3.eps");
 
  
  TCanvas *d = new TCanvas("d","d"); 
  d->Divide(1,2);
  
  d->cd(1); 
  MC_PIDCalib_EffRatio_KKmumu->GetYaxis()->SetRangeUser(0.9,1.2);
  MC_PIDCalib_EffRatio_KKmumu->Draw();  
  MC_PIDCalib_EffRatio_KKmumu->Write();

  d->cd(2); 
  MC_PIDCalib_EffRatio_pipimumu->GetYaxis()->SetRangeUser(0.9,1.2);
  MC_PIDCalib_EffRatio_pipimumu->Draw();
  MC_PIDCalib_EffRatio_pipimumu->Write();
 
  d->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_and_PIDCalib_."+polarity+"_4.eps");
  d->Write();
  
  a->Write();
  b->Write();
  c->Write();

  fout->Close();
}

void  evaluatePIDCalibSingleMuonEfficiency(TString polarity,int muonIndex, TString binningScheme="default") {
  //both muons, integtated Pt
  //BUT:take same pt range also for MC efficienciencies(800-8000)


  TString pathToPIDCalibFile;
  
  pathToPIDCalibFile = "/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/resultsSingleMuon/"+polarity+"/";
  if(binningScheme=="alternative") pathToPIDCalibFile = "/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/resultsSingleMuon/alternativeBinning/"+polarity+"/";
  if(binningScheme=="alternative2") pathToPIDCalibFile = "/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/resultsSingleMuon/alternativeBinning2/"+polarity+"/";


  dcastyle();

  TFile *fout = new TFile(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_mu%i_",muonIndex)+binningScheme+"_"+polarity+".root","RECREATE");
  
  TH1* AbsEff_Kpimumu = new TH1D("AbsEff_Kpimumu","AbsEff_Kpimumu in bins of dimuon mass",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1* AbsMCEff_Kpimumu = new TH1D("AbsMCEff_Kpimumu","Abs MC Eff_Kpimumu in bins of dimuon mass",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1* AbsEff_KKmumu = new TH1D("AbsEff_KKmumu","AbsEff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsMCEff_KKmumu = new TH1D("AbsMCEff_KKmumu","Abs MC Eff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relEff_KKmumu = new TH1D("relEff_KKmumu","relEff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_KKmumu = new TH1D("relMCEff_KKmumu","rel MC Eff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* MC_PIDCalib_EffRatio_KKmumu= new TH1D("MC_PIDCalib_EffRatio_KKmumu","data and MC ratio of relative signal/nomralization EffRatio",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1* AbsEff_pipimumu = new TH1D("AbsEff_pipimumu","AbsEff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* AbsMCEff_pipimumu = new TH1D("AbsMCEff_pipimumu","Abs MC Eff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relEff_pipimumu = new TH1D("relEff_pipimumu","relEff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relMCEff_pipimumu = new TH1D("relMCEff_pipimumu","rel MC Eff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* MC_PIDCalib_EffRatio_pipimumu= new TH1D("MC_PIDCalib_EffRatio_pipimumu","data and MC ratio of relative signal/nomralization EffRatio",sizeof(binspipi)/sizeof(double)-1,binspipi);

 
  double tempEff, tempEff_norm,tempEff_err, tempEff_norm_err;
  double tempMCEff, tempMCEff_norm,tempMCEff_norm_err,tempMCEff_err;
  double relMCEff, relMCEff_err,relEff_err;
  double ratioEffRatios,ratioEffRatios_err;

  TString MC_file;

  TString muonPIDCut=TString::Format("mu%i_ProbNNmu>0.5&&mu%i_ProbNNghost<0.5",muonIndex,muonIndex);
  TString cutPT= TString::Format("&&mu%i_PT>800&&mu%i_PT<8000",muonIndex,muonIndex);

  //norm mode
  for(int i=0; i<rangesKpi_low.size();++i){
    TString fIn = pathToPIDCalibFile+TString::Format("D2Kpimumu_%.1f_%.1f_800.0_8000.0_mu%i.root",rangesKpi_low[0],rangesKpi_high[0],muonIndex);  
    tempEff_norm = evaluatePIDCalibEfficiency(fIn).first;
    tempEff_norm_err = evaluatePIDCalibEfficiency(fIn).second;
    AbsEff_Kpimumu->SetBinContent(i+1,tempEff_norm);
    AbsEff_Kpimumu->SetBinError(i+1,tempEff_norm_err);
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_",rangesKpi_low[0],rangesKpi_high[0])+polarity+".root";
    tempMCEff_norm = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"mu0_MuonNShared==0&&mu1_MuonNShared==0"+cutPT).first;
    tempMCEff_norm_err = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"mu0_MuonNShared==0&&mu1_MuonNShared==0"+cutPT).second;
    AbsMCEff_Kpimumu->SetBinContent(i+1,tempMCEff_norm);
    AbsMCEff_Kpimumu->SetBinError(i+1,tempMCEff_norm_err);
  }

  //KKmumu mode
  for(int i=0; i<rangesKK_low.size();++i){
    TString fIn = pathToPIDCalibFile+TString::Format("D2KKmumu_%.1f_%.1f_800.0_8000.0_mu%i.root",rangesKK_low[i],rangesKK_high[i],muonIndex);  
    tempEff = evaluatePIDCalibEfficiency(fIn).first;
    tempEff_err = evaluatePIDCalibEfficiency(fIn).second;
    AbsEff_KKmumu->SetBinContent(i+1,tempEff);
    AbsEff_KKmumu->SetBinError(i+1,tempEff_err);
    relEff_KKmumu->SetBinContent(i+1,tempEff/AbsEff_Kpimumu->GetBinContent(1));
    std::cout<<rangesKK_low[i]<<"-"<<rangesKK_high[i]<<"  "<<tempEff<<"  "<<AbsEff_Kpimumu->GetBinContent(1)<<std::endl;
    relEff_err=TMath::Sqrt(TMath::Power(tempEff_err/AbsEff_Kpimumu->GetBinContent(1),2) + TMath::Power(tempEff*AbsEff_Kpimumu->GetBinError(1)/(AbsEff_Kpimumu->GetBinContent(1)*AbsEff_Kpimumu->GetBinContent(1)),2));
    relEff_KKmumu->SetBinError(i+1,relEff_err);
    //MC 
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_",rangesKK_low[i],rangesKK_high[i])+polarity+".root";
    tempMCEff = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"mu0_MuonNShared==0&&mu1_MuonNShared==0"+cutPT).first;
    tempMCEff_err = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"mu0_MuonNShared==0&&mu1_MuonNShared==0"+cutPT).second;
    AbsMCEff_KKmumu->SetBinContent(i+1,tempMCEff);
    AbsMCEff_KKmumu->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(1)/(AbsMCEff_Kpimumu->GetBinContent(1)*AbsMCEff_Kpimumu->GetBinContent(1)),2) );
    relMCEff_KKmumu->SetBinContent(i+1,relMCEff);
    relMCEff_KKmumu->SetBinError(i+1,relMCEff_err);
    //relative MC and PID calib
    ratioEffRatios=relMCEff_KKmumu->GetBinContent(i+1)/relEff_KKmumu->GetBinContent(i+1);
    ratioEffRatios_err=TMath::Sqrt(TMath::Power(relMCEff_KKmumu->GetBinError(i+1)/relEff_KKmumu->GetBinContent(i+1),2) + TMath::Power(relMCEff_KKmumu->GetBinContent(i+1)*relEff_KKmumu->GetBinError(i+1)/(relEff_KKmumu->GetBinContent(i+1)*relEff_KKmumu->GetBinContent(i+1)),2));
    MC_PIDCalib_EffRatio_KKmumu->SetBinContent(i+1,ratioEffRatios);
    MC_PIDCalib_EffRatio_KKmumu->SetBinError(i+1,ratioEffRatios_err);
  }

  //pipimumu mode
  for(int i=0; i<rangespipi_low.size();++i){
    TString fIn = pathToPIDCalibFile+TString::Format("D2pipimumu_%.1f_%.1f_800.0_8000.0_mu%i.root",rangespipi_low[i],rangespipi_high[i],muonIndex);  
    tempEff = evaluatePIDCalibEfficiency(fIn).first;
    tempEff_err = evaluatePIDCalibEfficiency(fIn).second;
    AbsEff_pipimumu->SetBinContent(i+1,tempEff);
    AbsEff_pipimumu->SetBinError(i+1,tempEff_err);
    relEff_pipimumu->SetBinContent(i+1,tempEff/AbsEff_Kpimumu->GetBinContent(1));
    relEff_err=TMath::Sqrt(TMath::Power(tempEff_err/AbsEff_Kpimumu->GetBinContent(1),2) + TMath::Power(tempEff*AbsEff_Kpimumu->GetBinError(1)/(AbsEff_Kpimumu->GetBinContent(1)*AbsEff_Kpimumu->GetBinContent(1)),2));
    relEff_pipimumu->SetBinError(i+1,relEff_err);
    //MC 
    MC_file = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_",rangespipi_low[i],rangespipi_high[i])+polarity+".root";
    tempMCEff = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"mu0_MuonNShared==0&&mu1_MuonNShared==0"+cutPT).first;
    tempMCEff_err = getMCEfficiency(MC_file,"BDT_Tree",muonPIDCut,"mu0_MuonNShared==0&&mu1_MuonNShared==0"+cutPT).second;
    AbsMCEff_pipimumu->SetBinContent(i+1,tempMCEff);
    AbsMCEff_pipimumu->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(1)/(AbsMCEff_Kpimumu->GetBinContent(1)*AbsMCEff_Kpimumu->GetBinContent(1)),2) );
    relMCEff_pipimumu->SetBinContent(i+1,relMCEff);
    relMCEff_pipimumu->SetBinError(i+1,relMCEff_err);
    //relative MC and PID calib
    ratioEffRatios=relMCEff_pipimumu->GetBinContent(i+1)/relEff_pipimumu->GetBinContent(i+1);
    ratioEffRatios_err=TMath::Sqrt(TMath::Power(relMCEff_pipimumu->GetBinError(i+1)/relEff_pipimumu->GetBinContent(i+1),2) + TMath::Power(relMCEff_pipimumu->GetBinContent(i+1)*relEff_pipimumu->GetBinError(i+1)/(relEff_pipimumu->GetBinContent(i+1)*relEff_pipimumu->GetBinContent(i+1)),2));
    MC_PIDCalib_EffRatio_pipimumu->SetBinContent(i+1,ratioEffRatios);
    MC_PIDCalib_EffRatio_pipimumu->SetBinError(i+1,ratioEffRatios_err);
  }


 
  fout->cd();

  TCanvas *a = new TCanvas("a","a"); 
  a->Divide(1,2);
  a->cd(1);
  AbsEff_KKmumu->GetYaxis()->SetRangeUser(0.2,.8);
  AbsEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  AbsEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}");
  AbsEff_KKmumu->Draw();
  AbsEff_KKmumu->Write();
  AbsMCEff_KKmumu->SetLineColor(kCyan);
  AbsMCEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsMCEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}");
  AbsMCEff_KKmumu->Draw("SAME");
  AbsMCEff_KKmumu->Write();
  
  a->cd(2);
  relEff_KKmumu->GetYaxis()->SetRangeUser(0.85,1.1);
  relEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  relEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relEff_KKmumu->Draw();
  relEff_KKmumu->Write();
  relMCEff_KKmumu->SetLineColor(kCyan);
  relMCEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  relMCEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relMCEff_KKmumu->Draw("SAME");
  relMCEff_KKmumu->Write();
  if(binningScheme=="default")
    a->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_mu%i_",muonIndex)+polarity+"_1.eps");

  
  TCanvas *b = new TCanvas("b","b"); 
 
  AbsEff_Kpimumu->GetYaxis()->SetRangeUser(0.2,.9);
  AbsEff_Kpimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  AbsEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}");
  AbsEff_Kpimumu->Draw();
  AbsEff_Kpimumu->Write();
  AbsMCEff_Kpimumu->SetLineColor(kCyan);
  AbsEff_Kpimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  AbsEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}");
  AbsMCEff_Kpimumu->Draw("SAME"); 
  AbsMCEff_Kpimumu->Write();
  if(binningScheme=="default")
    b->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_mu%i_",muonIndex)+polarity+"_2.eps");



  TCanvas *c = new TCanvas("c","c"); 
  c->Divide(1,2);
  c->cd(1);
  AbsEff_pipimumu->GetYaxis()->SetRangeUser(0.2,.9);
  AbsEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  AbsEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}");
  AbsEff_pipimumu->Draw();
  AbsEff_pipimumu->Write();
  AbsMCEff_pipimumu->SetLineColor(kCyan);
  AbsMCEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  AbsMCEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}");
  AbsMCEff_pipimumu->Draw("SAME");
  AbsMCEff_pipimumu->Write();
  
  c->cd(2);  
  relEff_pipimumu->GetYaxis()->SetRangeUser(.85,1.2);
  relEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  relEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relEff_pipimumu->Draw();
  relEff_pipimumu->Write();
  relMCEff_pipimumu->SetLineColor(kCyan);
  relEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  relEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relMCEff_pipimumu->Draw("SAME");
  relMCEff_pipimumu->Write();
  c->Write();
  if(binningScheme=="default")
    c->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_mu%i_",muonIndex)+polarity+"_3.eps");
  
  TCanvas *d = new TCanvas("d","d"); 
  d->Divide(1,2);
  
  d->cd(1); 
  MC_PIDCalib_EffRatio_KKmumu->GetYaxis()->SetRangeUser(0.9,1.2);
  MC_PIDCalib_EffRatio_KKmumu->Draw();  
  MC_PIDCalib_EffRatio_KKmumu->Write();

  d->cd(2); 
  MC_PIDCalib_EffRatio_pipimumu->GetYaxis()->SetRangeUser(0.9,1.2);
  MC_PIDCalib_EffRatio_pipimumu->Draw();
  MC_PIDCalib_EffRatio_pipimumu->Write();
 
  d->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/SingleMuon_PID_MC_and_PIDCalib_mu%i_",muonIndex)+polarity+"_4.eps");

  d->Write();
  
  a->Write();
  b->Write();
  c->Write();

  fout->Close();
}


void  evaluateMCSingleMuonEfficiency(int muonIndex) {


  dcastyle();

  TFile *fout = new TFile(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon%i_PID_MC.root",muonIndex),"RECREATE");
  
  TH1* AbsMCEff_Kpimumu = new TH1D(TString::Format("AbsMCEff_Kpimumu_mu%i",muonIndex),"Abs MC Eff_Kpimumu in bins of dimuon mass",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1* AbsMCEff_KKmumu = new TH1D(TString::Format("AbsMCEff_KKmumu_mu%i",muonIndex),"Abs MC Eff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_KKmumu = new TH1D(TString::Format("relMCEff_KKmumu_mu%i",muonIndex),"rel MC Eff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsMCEff_pipimumu = new TH1D(TString::Format("AbsMCEff_pipimumu_mu%i",muonIndex),"Abs MC Eff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relMCEff_pipimumu = new TH1D(TString::Format("relMCEff_pipimumu_mu%i",muonIndex),"rel MC Eff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1* AbsMCEff_Kpimumu_highPtCut = new TH1D(TString::Format("AbsMCEff_Kpimumu_highPtCut_mu%i",muonIndex),"Abs MC Eff_Kpimumu in bins of dimuon mass pt<800",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1* AbsMCEff_KKmumu_highPtCut = new TH1D(TString::Format("AbsMCEff_KKmumu_highPtCut_mu%i",muonIndex),"Abs MC Eff_KKmumu in bins of dimuon mass pt<800",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_KKmumu_highPtCut = new TH1D(TString::Format("relMCEff_KKmumu_highPtCut_mu%i",muonIndex),"rel MC Eff_KKmumu in bins of dimuon mass pt<800",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsMCEff_pipimumu_highPtCut = new TH1D(TString::Format("AbsMCEff_pipimumu_highPtCut_mu%i",muonIndex),"Abs MC Eff_pipimumu in bins of dimuon mass pt<800",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relMCEff_pipimumu_highPtCut = new TH1D(TString::Format("relMCEff_pipimumu_highPtCut_mu%i",muonIndex),"rel MC Eff_pipimumu in bins of dimuon mass pt<800",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1* AbsMCEff_Kpimumu_lowPtCut = new TH1D(TString::Format("AbsMCEff_Kpimumu_lowPtCut_mu%i",muonIndex),"Abs MC Eff_Kpimumu in bins of dimuon mass pt<800",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1* AbsMCEff_KKmumu_lowPtCut = new TH1D(TString::Format("AbsMCEff_KKmumu_lowPtCut_mu%i",muonIndex),"Abs MC Eff_KKmumu in bins of dimuon mass pt<800",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_KKmumu_lowPtCut = new TH1D(TString::Format("relMCEff_KKmumu_lowPtCut_mu%i",muonIndex),"rel MC Eff_KKmumu in bins of dimuon mass pt<800",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsMCEff_pipimumu_lowPtCut = new TH1D(TString::Format("AbsMCEff_pipimumu_lowPtCut_mu%i",muonIndex),"Abs MC Eff_pipimumu in bins of dimuon mass pt<800",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relMCEff_pipimumu_lowPtCut = new TH1D(TString::Format("relMCEff_pipimumu_lowPtCut_mu%i",muonIndex),"rel MC Eff_pipimumu in bins of dimuon mass pt<800",sizeof(binspipi)/sizeof(double)-1,binspipi);

  double tempEff, tempEff_norm,tempEff_err, tempEff_norm_err;
  double tempMCEff, tempMCEff_norm,tempMCEff_norm_err,tempMCEff_err;
  double relMCEff, relMCEff_err,relEff_err;
  double ratioEffRatios,ratioEffRatios_err;

  TString MC_file_up,MC_file_dw; 

  TString muonPIDcut=TString::Format("mu%i_ProbNNmu>0.5&&mu%i_ProbNNghost<0.5",muonIndex,muonIndex);
  //also look at efficiency for events with pt higher and lower the 800 MeV cut
  TString ptCut= TString::Format("mu%i_PT<800&&mu0_MuonNShared==0&&mu1_MuonNShared==0",muonIndex);
  TString ptCut_low= TString::Format("mu%i_PT>800&&mu0_MuonNShared==0&&mu1_MuonNShared==0",muonIndex);
  
  ///////////////////////////////////////
  // 
  //  Fill the MC efficiencies 
  //
  ///////////////////////////////////////

  //norm mode
  for(int i=0; i<rangesKpi_low.size();++i){
    MC_file_up = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_",rangesKpi_low[0],rangesKpi_high[0])+"magUp.root";
    MC_file_dw = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_",rangesKpi_low[0],rangesKpi_high[0])+"magDw.root";
    //std::cout<<"noCut!"<<std::endl;
    tempMCEff_norm = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").first;
    tempMCEff_norm_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").second;
    AbsMCEff_Kpimumu->SetBinContent(i+1,tempMCEff_norm);
    AbsMCEff_Kpimumu->SetBinError(i+1,tempMCEff_norm_err);
    //pt<800
    //std::cout<<"first Cut!"<<std::endl;
    tempMCEff_norm = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).first;
    tempMCEff_norm_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).second;
    AbsMCEff_Kpimumu_highPtCut->SetBinContent(i+1,tempMCEff_norm);
    AbsMCEff_Kpimumu_highPtCut->SetBinError(i+1,tempMCEff_norm_err);
    ///>>800
    //std::cout<<"second Cut!"<<std::endl;
    tempMCEff_norm = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).first;
    tempMCEff_norm_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).second;
    AbsMCEff_Kpimumu_lowPtCut->SetBinContent(i+1,tempMCEff_norm);
    AbsMCEff_Kpimumu_lowPtCut->SetBinError(i+1,tempMCEff_norm_err);
    
  }

  //KKmumu mode
  for(int i=0; i<rangesKK_low.size();++i){
    //MC 
    MC_file_up = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_",rangesKK_low[i],rangesKK_high[i])+"magUp.root";
    MC_file_dw = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_",rangesKK_low[i],rangesKK_high[i])+"magDw.root";
    //std::cout<<"noCut!"<<std::endl;
    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").second;
    AbsMCEff_KKmumu->SetBinContent(i+1,tempMCEff);
    AbsMCEff_KKmumu->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(1)/(AbsMCEff_Kpimumu->GetBinContent(1)*AbsMCEff_Kpimumu->GetBinContent(1)),2) );
    relMCEff_KKmumu->SetBinContent(i+1,relMCEff);
    relMCEff_KKmumu->SetBinError(i+1,relMCEff_err);
    //std::cout<<"nocut "<<tempMCEff<<"/"<<AbsMCEff_Kpimumu->GetBinContent(1)<<"="<<tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1)<<" "<<relMCEff_KKmumu->GetBinContent(i+1)<<std::endl;
    //pt<800
    //std::cout<<"fitst Cut!"<<std::endl;
    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).second;
    AbsMCEff_KKmumu_highPtCut->SetBinContent(i+1,tempMCEff);
    AbsMCEff_KKmumu_highPtCut->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu_highPtCut->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu_highPtCut->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu_highPtCut->GetBinError(1)/(AbsMCEff_Kpimumu_highPtCut->GetBinContent(1)*AbsMCEff_Kpimumu_highPtCut->GetBinContent(1)),2) );
    relMCEff_KKmumu_highPtCut->SetBinContent(i+1,relMCEff);
    relMCEff_KKmumu_highPtCut->SetBinError(i+1,relMCEff_err);
    //std::cout<<"pt>800: " <<tempMCEff<<"/"<<AbsMCEff_Kpimumu_highPtCut->GetBinContent(1)<<"="<<tempMCEff/AbsMCEff_Kpimumu_highPtCut->GetBinContent(1)<<" "<<relMCEff_KKmumu_highPtCut->GetBinContent(i+1)<<std::endl;
    //std::cout<<"second Cut!"<<std::endl;
    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).second;
    AbsMCEff_KKmumu_lowPtCut->SetBinContent(i+1,tempMCEff);
    AbsMCEff_KKmumu_lowPtCut->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu_lowPtCut->GetBinError(1)/(AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1)*AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1)),2) );
    relMCEff_KKmumu_lowPtCut->SetBinContent(i+1,relMCEff);
    relMCEff_KKmumu_lowPtCut->SetBinError(i+1,relMCEff_err);
    //std::cout<<"pt<800: " <<tempMCEff<<"/"<<AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1)<<"="<<tempMCEff/AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1)<<" "<<relMCEff_KKmumu_lowPtCut->GetBinContent(i+1)<<std::endl;
  
 }

  //pipimumu mode
  for(int i=0; i<rangespipi_low.size();++i){
    //MC 
    std::cout<<"no Cut!"<<std::endl;
    MC_file_up = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_",rangespipi_low[i],rangespipi_high[i])+"magUp.root";
    MC_file_dw = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_",rangespipi_low[i],rangespipi_high[i])+"magDw.root";
    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").second;
    //tempMCEff = getMCEfficiency(MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").first;
    //tempMCEff_err = getMCEfficiency(MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").second;
    //std::cout<<" Q2 bin" <<i<<" TEMP EFF: "<<tempMCEff<<std::endl;
    AbsMCEff_pipimumu->SetBinContent(i+1,tempMCEff);
    AbsMCEff_pipimumu->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(1)/(AbsMCEff_Kpimumu->GetBinContent(1)*AbsMCEff_Kpimumu->GetBinContent(1)),2) );
    relMCEff_pipimumu->SetBinContent(i+1,relMCEff);
    relMCEff_pipimumu->SetBinError(i+1,relMCEff_err);
    //pt<800
    //std::cout<<"fist Cut!"<<std::endl;
    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).second;
    AbsMCEff_pipimumu_highPtCut->SetBinContent(i+1,tempMCEff);
    AbsMCEff_pipimumu_highPtCut->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu_highPtCut->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu_highPtCut->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu_highPtCut->GetBinError(1)/(AbsMCEff_Kpimumu_highPtCut->GetBinContent(1)*AbsMCEff_Kpimumu_highPtCut->GetBinContent(1)),2) );
    relMCEff_pipimumu_highPtCut->SetBinContent(i+1,relMCEff);
    relMCEff_pipimumu_highPtCut->SetBinError(i+1,relMCEff_err);
    //std::cout<<"second Cut!"<<std::endl;
    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).second;
    AbsMCEff_pipimumu_lowPtCut->SetBinContent(i+1,tempMCEff);
    AbsMCEff_pipimumu_lowPtCut->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu_lowPtCut->GetBinError(1)/(AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1)*AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1)),2) );
    relMCEff_pipimumu_lowPtCut->SetBinContent(i+1,relMCEff);
    relMCEff_pipimumu_lowPtCut->SetBinError(i+1,relMCEff_err);

  }


  fout->cd();
  TCanvas *a = new TCanvas("a","a"); 
  a->Divide(1,2);
  a->cd(1);
  AbsMCEff_KKmumu->SetLineColor(kCyan);
  AbsMCEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsMCEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}");
  AbsMCEff_KKmumu->GetYaxis()->SetRangeUser(0.5,1.5);
  AbsMCEff_KKmumu->Draw();
  AbsMCEff_KKmumu_highPtCut->SetLineColor(kViolet);
  AbsMCEff_KKmumu_highPtCut->Draw("SAME");  
  AbsMCEff_KKmumu_lowPtCut->SetLineColor(kBlue);
  AbsMCEff_KKmumu_lowPtCut->Draw("SAME");  
  AbsMCEff_KKmumu->Write();
  
  a->cd(2);
  relMCEff_KKmumu->SetLineColor(kCyan);
  relMCEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)");
  relMCEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relMCEff_KKmumu->GetYaxis()->SetRangeUser(0.95,1.05);
  relMCEff_KKmumu->Draw();
  relMCEff_KKmumu_highPtCut->SetLineColor(kViolet);
  relMCEff_KKmumu_highPtCut->Draw("SAME");  
  relMCEff_KKmumu_lowPtCut->SetLineColor(kBlue);
  relMCEff_KKmumu_lowPtCut->Draw("SAME");  
  relMCEff_KKmumu->Write();
  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon_PID_MC_1.eps");

  
  TCanvas *b = new TCanvas("b","b"); 
  AbsMCEff_Kpimumu->SetLineColor(kCyan);
  AbsMCEff_Kpimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsMCEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}");
  AbsMCEff_Kpimumu->GetYaxis()->SetRangeUser(0.5,1.5);
  AbsMCEff_Kpimumu->Draw(); 
  AbsMCEff_Kpimumu_highPtCut->SetLineColor(kViolet);
  AbsMCEff_Kpimumu_highPtCut->Draw("SAME");  
  AbsMCEff_Kpimumu_lowPtCut->SetLineColor(kBlue);
  AbsMCEff_Kpimumu_lowPtCut->Draw("SAME");  
  AbsMCEff_Kpimumu->Write();

  b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon_PID_MC_2.eps");


  TCanvas *c = new TCanvas("c","c"); 
  c->Divide(1,2);
  c->cd(1);
  AbsMCEff_pipimumu->SetLineColor(kCyan);
  AbsMCEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsMCEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}");
  AbsMCEff_pipimumu->GetYaxis()->SetRangeUser(0.5,1.5);
  AbsMCEff_pipimumu->Draw();
  AbsMCEff_pipimumu_highPtCut->SetLineColor(kViolet);
  AbsMCEff_pipimumu_highPtCut->Draw("SAME");  
  AbsMCEff_pipimumu_lowPtCut->SetLineColor(kBlue);
  AbsMCEff_pipimumu_lowPtCut->Draw("SAME");  
  AbsMCEff_pipimumu->Write();
  
  c->cd(2);  
  relMCEff_pipimumu->SetLineColor(kCyan);
  relMCEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  relMCEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relMCEff_pipimumu->GetYaxis()->SetRangeUser(0.95,1.1);
  relMCEff_pipimumu->Draw();
  relMCEff_pipimumu_highPtCut->SetLineColor(kViolet);
  relMCEff_pipimumu_highPtCut->Draw("SAME");  
  relMCEff_pipimumu_lowPtCut->SetLineColor(kBlue);
  relMCEff_pipimumu_lowPtCut->Draw("SAME");  
  relMCEff_pipimumu->Write();
  //c->Write();
  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/Muon_PID_MC_3.eps");
 
  fout->Write();
  //a->Write();
  //b->Write();
  //c->Write();
  

}
void evaluateMCDoubleMuonEfficiency() {

  //evaluates the double muon PID efficiency in q2 bins, magUp and Down merged.

  dcastyle();

  TFile *fout = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC.root","RECREATE");
  
  TH1* AbsMCEff_Kpimumu = new TH1D("AbsMCEff_Kpimumu","Abs MC Eff_Kpimumu in bins of dimuon mass",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1* AbsMCEff_KKmumu = new TH1D("AbsMCEff_KKmumu","Abs MC Eff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_KKmumu = new TH1D("relMCEff_KKmumu","rel MC Eff_KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsMCEff_pipimumu = new TH1D("AbsMCEff_pipimumu","Abs MC Eff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relMCEff_pipimumu = new TH1D("relMCEff_pipimumu","rel MC Eff_pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1* AbsMCEff_Kpimumu_highPtCut = new TH1D("AbsMCEff_Kpimumu_highPtCut","Abs MC Eff_Kpimumu in bins of dimuon mass pt<800",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1* AbsMCEff_KKmumu_highPtCut = new TH1D("AbsMCEff_KKmumu_highPtCut","Abs MC Eff_KKmumu in bins of dimuon mass pt<800",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_KKmumu_highPtCut = new TH1D("relMCEff_KKmumu_highPtCut","rel MC Eff_KKmumu in bins of dimuon mass pt<800",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsMCEff_pipimumu_highPtCut = new TH1D("AbsMCEff_pipimumu_highPtCut","Abs MC Eff_pipimumu in bins of dimuon mass pt<800",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relMCEff_pipimumu_highPtCut = new TH1D("relMCEff_pipimumu_highPtCut","rel MC Eff_pipimumu in bins of dimuon mass pt<800",sizeof(binspipi)/sizeof(double)-1,binspipi);


  TH1* AbsMCEff_Kpimumu_lowPtCut = new TH1D("AbsMCEff_Kpimumu_lowPtCut","Abs MC Eff_Kpimumu in bins of dimuon mass pt<800",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1* AbsMCEff_KKmumu_lowPtCut = new TH1D("AbsMCEff_KKmumu_lowPtCut","Abs MC Eff_KKmumu in bins of dimuon mass pt<800",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* relMCEff_KKmumu_lowPtCut = new TH1D("relMCEff_KKmumu_lowPtCut","rel MC Eff_KKmumu in bins of dimuon mass pt<800",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* AbsMCEff_pipimumu_lowPtCut = new TH1D("AbsMCEff_pipimumu_lowPtCut","Abs MC Eff_pipimumu in bins of dimuon mass pt<800",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1* relMCEff_pipimumu_lowPtCut = new TH1D("relMCEff_pipimumu_lowPtCut","rel MC Eff_pipimumu in bins of dimuon mass pt<800",sizeof(binspipi)/sizeof(double)-1,binspipi);

  double tempEff, tempEff_norm,tempEff_err, tempEff_norm_err;
  double tempMCEff, tempMCEff_norm,tempMCEff_norm_err,tempMCEff_err;
  double relMCEff, relMCEff_err,relEff_err;
  double ratioEffRatios,ratioEffRatios_err;

  TString MC_file_up,MC_file_dw;

  TString muonPIDcut="mu0_ProbNNmu>0.5&&mu0_ProbNNghost<0.5&&mu1_ProbNNmu>0.5&&mu1_ProbNNghost<0.5";
  TString ptCut = "mu0_PT<800&&mu1_PT<800";
  TString ptCut_low = "mu0_PT>800&&mu1_PT>800";

  //norm mode
  for(int i=0; i<rangesKpi_low.size();++i){
    MC_file_up = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_",rangesKpi_low[0],rangesKpi_high[0])+"magUp.root";
    MC_file_dw = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_",rangesKpi_low[0],rangesKpi_high[0])+"magDw.root";
    tempMCEff_norm = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").first;
    tempMCEff_norm_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").second;
    AbsMCEff_Kpimumu->SetBinContent(i+1,tempMCEff_norm);
    AbsMCEff_Kpimumu->SetBinError(i+1,tempMCEff_norm_err);
    //pt<800
    tempMCEff_norm = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).first;
    tempMCEff_norm_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).second;
    AbsMCEff_Kpimumu_highPtCut->SetBinContent(i+1,tempMCEff_norm);
    AbsMCEff_Kpimumu_highPtCut->SetBinError(i+1,tempMCEff_norm_err);

    tempMCEff_norm = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).first;
    tempMCEff_norm_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).second;
    AbsMCEff_Kpimumu_lowPtCut->SetBinContent(i+1,tempMCEff_norm);
    AbsMCEff_Kpimumu_lowPtCut->SetBinError(i+1,tempMCEff_norm_err);

  }

  //KKmumu mode
  for(int i=0; i<rangesKK_low.size();++i){
    //MC 
    MC_file_up = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_",rangesKK_low[i],rangesKK_high[i])+"magUp.root";
    MC_file_dw = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_",rangesKK_low[i],rangesKK_high[i])+"magDw.root";
    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").second;
    AbsMCEff_KKmumu->SetBinContent(i+1,tempMCEff);
    AbsMCEff_KKmumu->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(1)/(AbsMCEff_Kpimumu->GetBinContent(1)*AbsMCEff_Kpimumu->GetBinContent(1)),2) );
    relMCEff_KKmumu->SetBinContent(i+1,relMCEff);
    relMCEff_KKmumu->SetBinError(i+1,relMCEff_err);
    //pt<800
    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).second;
    AbsMCEff_KKmumu_highPtCut->SetBinContent(i+1,tempMCEff);
    AbsMCEff_KKmumu_highPtCut->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu_highPtCut->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu_highPtCut->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu_highPtCut->GetBinError(1)/(AbsMCEff_Kpimumu_highPtCut->GetBinContent(1)*AbsMCEff_Kpimumu_highPtCut->GetBinContent(1)),2) );
    relMCEff_KKmumu_highPtCut->SetBinContent(i+1,relMCEff);
    relMCEff_KKmumu_highPtCut->SetBinError(i+1,relMCEff_err);

    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).second;
    AbsMCEff_KKmumu_lowPtCut->SetBinContent(i+1,tempMCEff);
    AbsMCEff_KKmumu_lowPtCut->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu_lowPtCut->GetBinError(1)/(AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1)*AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1)),2) );
    relMCEff_KKmumu_lowPtCut->SetBinContent(i+1,relMCEff);
    relMCEff_KKmumu_lowPtCut->SetBinError(i+1,relMCEff_err);
  }

  //pipimumu mode
  for(int i=0; i<rangespipi_low.size();++i){
    //MC 
    MC_file_up = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_",rangespipi_low[i],rangespipi_high[i])+"magUp.root";
    MC_file_dw = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_",rangespipi_low[i],rangespipi_high[i])+"magDw.root";
    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,"mu0_MuonNShared==0&&mu1_MuonNShared==0").second;
    AbsMCEff_pipimumu->SetBinContent(i+1,tempMCEff);
    AbsMCEff_pipimumu->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu->GetBinError(1)/(AbsMCEff_Kpimumu->GetBinContent(1)*AbsMCEff_Kpimumu->GetBinContent(1)),2) );
    relMCEff_pipimumu->SetBinContent(i+1,relMCEff);
    relMCEff_pipimumu->SetBinError(i+1,relMCEff_err);
    //pt<800
    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut).second;
    AbsMCEff_pipimumu_highPtCut->SetBinContent(i+1,tempMCEff);
    AbsMCEff_pipimumu_highPtCut->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu_highPtCut->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu_highPtCut->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu_highPtCut->GetBinError(1)/(AbsMCEff_Kpimumu_highPtCut->GetBinContent(1)*AbsMCEff_Kpimumu_highPtCut->GetBinContent(1)),2) );
    relMCEff_pipimumu_highPtCut->SetBinContent(i+1,relMCEff);
    relMCEff_pipimumu_highPtCut->SetBinError(i+1,relMCEff_err);

    tempMCEff = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).first;
    tempMCEff_err = getMCEfficiencyFrom2Files(MC_file_up,MC_file_dw,"BDT_Tree",muonPIDcut,ptCut_low).second;
    AbsMCEff_pipimumu_lowPtCut->SetBinContent(i+1,tempMCEff);
    AbsMCEff_pipimumu_lowPtCut->SetBinError(i+1,tempMCEff_err);
    relMCEff=tempMCEff/AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1);
    relMCEff_err=TMath::Sqrt( TMath::Power(tempMCEff_err/AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1),2)+TMath::Power(tempMCEff*AbsMCEff_Kpimumu_lowPtCut->GetBinError(1)/(AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1)*AbsMCEff_Kpimumu_lowPtCut->GetBinContent(1)),2) );
    relMCEff_pipimumu_lowPtCut->SetBinContent(i+1,relMCEff);
    relMCEff_pipimumu_lowPtCut->SetBinError(i+1,relMCEff_err);

  }


  fout->cd();
  TCanvas *a = new TCanvas("a","a"); 
  a->Divide(1,2);
  a->cd(1);
  AbsMCEff_KKmumu->SetLineColor(kCyan);
  AbsMCEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  AbsMCEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}");
  AbsMCEff_KKmumu->GetYaxis()->SetRangeUser(0.4,0.8);
  AbsMCEff_KKmumu->Draw();
  AbsMCEff_KKmumu_highPtCut->SetLineColor(kViolet);
  AbsMCEff_KKmumu_highPtCut->Draw("SAME");  
  AbsMCEff_KKmumu_lowPtCut->SetLineColor(kBlue);
  AbsMCEff_KKmumu_lowPtCut->Draw("SAME");  
  AbsMCEff_KKmumu->Write();
  
  a->cd(2);
  relMCEff_KKmumu->SetLineColor(kCyan);
  relMCEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  relMCEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{D->KK#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relMCEff_KKmumu->GetYaxis()->SetRangeUser(0.9,1.2);
  relMCEff_KKmumu->Draw();
  relMCEff_KKmumu_highPtCut->SetLineColor(kViolet);
  relMCEff_KKmumu_highPtCut->Draw("SAME");  
  relMCEff_KKmumu_lowPtCut->SetLineColor(kBlue);
  relMCEff_KKmumu_lowPtCut->Draw("SAME");  
  relMCEff_KKmumu->Write();
  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_1.eps");

  
  TCanvas *b = new TCanvas("b","b"); 
  AbsMCEff_Kpimumu->SetLineColor(kCyan);
  AbsMCEff_Kpimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  AbsMCEff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{D->K#pi#mu#mu}");
  AbsMCEff_Kpimumu->GetYaxis()->SetRangeUser(0.3,0.9);
  AbsMCEff_Kpimumu->Draw(); 
  AbsMCEff_Kpimumu_highPtCut->SetLineColor(kViolet);
  AbsMCEff_Kpimumu_highPtCut->Draw("SAME");  
  AbsMCEff_Kpimumu_lowPtCut->SetLineColor(kBlue);
  AbsMCEff_Kpimumu_lowPtCut->Draw("SAME");  
  AbsMCEff_Kpimumu->Write();

  b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_2.eps");


  TCanvas *c = new TCanvas("c","c"); 
  c->Divide(1,2);
  c->cd(1);
  AbsMCEff_pipimumu->SetLineColor(kCyan);
  AbsMCEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)");
  AbsMCEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}");
  AbsMCEff_pipimumu->GetYaxis()->SetRangeUser(0.3,0.9);
  AbsMCEff_pipimumu->Draw();
  AbsMCEff_pipimumu_highPtCut->SetLineColor(kViolet);
  AbsMCEff_pipimumu_highPtCut->Draw("SAME");  
  AbsMCEff_pipimumu_lowPtCut->SetLineColor(kBlue);
  AbsMCEff_pipimumu_lowPtCut->Draw("SAME");  
  AbsMCEff_pipimumu->Write();
  
  c->cd(2);  
  relMCEff_pipimumu->SetLineColor(kCyan);
  relMCEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  relMCEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{D->#pi#pi#mu#mu}/#epsilon_{D->K#pi#mu#mu}");
  relMCEff_pipimumu->GetYaxis()->SetRangeUser(0.9,1.2);
  relMCEff_pipimumu->Draw();
  relMCEff_pipimumu_highPtCut->SetLineColor(kViolet);
  relMCEff_pipimumu_highPtCut->Draw("SAME");  
  relMCEff_pipimumu_lowPtCut->SetLineColor(kBlue);
  relMCEff_pipimumu_lowPtCut->Draw("SAME");  
  relMCEff_pipimumu->Write();
  //c->Write();
  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/doubleMuon_PID_MC_3.eps");
 
  fout->Write();
  //a->Write();
  //b->Write();
  //c->Write();
  
  fout->Close();
}


void copyMCTuplesWithCut(TString fileIn,TString target, TString nameTree,TString cut) {
  
  //takes an input tree and saves a cipz with cut string applied
  TFile* file= new TFile(fileIn,"OPEN");
  TTree* tuple = (TTree*) file->Get(nameTree);
  TFile* output = new TFile(target,"RECREATE");
  TTree* sel_tuple = tuple->CopyTree(cut);
  sel_tuple->Write();
  output->Close();
  file->Close();


}

void splitMCtuplesInPt(bool cutNshared){
//split in bins of pt 
//cuts on pt of single muon m0
  TString Nshared;
  if(cutNshared)  Nshared = "&&mu1_MuonNShared==0&&mu0_MuonNShared==0";
  else Nshared = "&&mu1_MuonNShared<100"; //cut effectively does nothing at all

 //Kpimumu KKmumu trained
for(int i=0; i<rangesKpi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPt_low.size();++j){                                                                                                                                       
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBins/MC_D2Kpimumu_%.1f_%.1f_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i],rangesPt_low[j],rangesPt_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f",rangesPt_low[j],rangesPt_high[j])+Nshared;
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

//KKmumumu  
for(int i=0; i<rangesKK_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPt_low.size();++j){                                                                                                                                        
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBins/MC_D2KKmumu_%.1f_%.1f_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i],rangesPt_low[j],rangesPt_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f",rangesPt_low[j],rangesPt_high[j])+Nshared;
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

//pipimumumu  
for(int i=0; i<rangespipi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPt_low.size();++j){                                                                                                                                        
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBins/MC_D2pipimumu_%.1f_%.1f_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i],rangesPt_low[j],rangesPt_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f",rangesPt_low[j],rangesPt_high[j]);
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }


for(int i=0; i<rangesKpi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPt_low.size();++j){                                                                                                                                       
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBins/MC_D2Kpimumu_%.1f_%.1f_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i],rangesPt_low[j],rangesPt_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f",rangesPt_low[j],rangesPt_high[j])+Nshared;
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

//KKmumumu  
for(int i=0; i<rangesKK_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPt_low.size();++j){                                                                                                                                        
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBins/MC_D2KKmumu_%.1f_%.1f_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i],rangesPt_low[j],rangesPt_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f",rangesPt_low[j],rangesPt_high[j])+Nshared;
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

//pipimumumu  
for(int i=0; i<rangespipi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPt_low.size();++j){                                                                                                                                        
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBins/MC_D2pipimumu_%.1f_%.1f_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i],rangesPt_low[j],rangesPt_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f",rangesPt_low[j],rangesPt_high[j])+Nshared;
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

}

void splitMCtuplesInPtOfBothMuons(bool cutNshared){
//split in bins of pt 

  TString Nshared;
  if(cutNshared)  Nshared = "&&mu1_MuonNShared==0&&mu0_MuonNShared==0";
  else Nshared = "&&mu1_MuonNShared<100"; //cut effectively does nothing at all

 //Kpimumu KKmumu trained
for(int i=0; i<rangesKpi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                                       
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBinsBothMuons/MC_D2Kpimumu_%.1f_%.1f_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f&&mu1_PT>%f&&mu1_PT<%f",rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j])+Nshared;
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

//KKmumumu  
for(int i=0; i<rangesKK_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                                        
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBinsBothMuons/MC_D2KKmumu_%.1f_%.1f_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f&&mu1_PT>%f&&mu1_PT<%f",rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j])+Nshared;
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

//pipimumumu  
for(int i=0; i<rangespipi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                                        
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBinsBothMuons/MC_D2pipimumu_%.1f_%.1f_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f&&mu1_PT>%f&&mu1_PT<%f",rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j]);
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }


for(int i=0; i<rangesKpi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                                       
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBinsBothMuons/MC_D2Kpimumu_%.1f_%.1f_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f&&mu1_PT>%f&&mu1_PT<%f",rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j])+Nshared;
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

//KKmumumu  
for(int i=0; i<rangesKK_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                                        
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBinsBothMuons/MC_D2KKmumu_%.1f_%.1f_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f&&mu1_PT>%f&&mu1_PT<%f",rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j])+Nshared;
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

//pipimumumu  
for(int i=0; i<rangespipi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                                        
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/PTBinsBothMuons/MC_D2pipimumu_%.1f_%.1f_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j]);
    TString cut =  TString::Format("mu0_PT>%f&&mu0_PT<%f&&mu1_PT>%f&&mu1_PT<%f",rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j])+Nshared;
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",cut);
    std::cout<<cut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

}

void splitMCtuplesPtCutOnSingleMuon(bool cutNshared,int index){
//appy cut on single muon  pt>800 

  TString Nshared;
  if(cutNshared)  Nshared = "&&mu1_MuonNShared==0&&mu0_MuonNShared==0";
  else Nshared = "&&mu1_MuonNShared<100"; //cut effectively does nothing at all
  TString ptCut;
  
 //Kpimumu KKmumu trained
for(int i=0; i<rangesKpi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                           ptCut=TString::Format("mu%i_PT>%f&&mu%i_PT<%f",index,rangesPtBothMuons_low[j],index,rangesPtBothMuons_high[j])+Nshared;            
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/BinsSingleMuon/MC_D2Kpimumu_%.1f_%.1f_%.1f_%.1f_mu%i_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],index);
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",ptCut);
    std::cout<<ptCut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }
for(int i=0; i<rangesKK_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                           ptCut=TString::Format("mu%i_PT>%f&&mu%i_PT<%f",index,rangesPtBothMuons_low[j],index,rangesPtBothMuons_high[j])+Nshared;            
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/BinsSingleMuon/MC_D2KKmumu_%.1f_%.1f_%.1f_%.1f_mu%i_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],index);
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",ptCut);
    std::cout<<ptCut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }
for(int i=0; i<rangespipi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                           ptCut=TString::Format("mu%i_PT>%f&&mu%i_PT<%f",index,rangesPtBothMuons_low[j],index,rangesPtBothMuons_high[j])+Nshared;            
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/BinsSingleMuon/MC_D2pipimumu_%.1f_%.1f_%.1f_%.1f_mu%i_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],index);
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",ptCut);
    std::cout<<ptCut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }
//MagDw
 //Kpimumu KKmumu trained
for(int i=0; i<rangesKpi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                           ptCut=TString::Format("mu%i_PT>%f&&mu%i_PT<%f",index,rangesPtBothMuons_low[j],index,rangesPtBothMuons_high[j])+Nshared;            
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/BinsSingleMuon/MC_D2Kpimumu_%.1f_%.1f_%.1f_%.1f_mu%i_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],index);
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",ptCut);
    std::cout<<ptCut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }
for(int i=0; i<rangesKK_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                           ptCut=TString::Format("mu%i_PT>%f&&mu%i_PT<%f",index,rangesPtBothMuons_low[j],index,rangesPtBothMuons_high[j])+Nshared;            
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/BinsSingleMuon/MC_D2KKmumu_%.1f_%.1f_%.1f_%.1f_mu%i_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],index);
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",ptCut);
    std::cout<<ptCut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }
for(int i=0; i<rangespipi_low.size();++i){                                                                                                                                       
  for(int j=0; j<rangesPtBothMuons_low.size();++j){                                                                                                                           ptCut=TString::Format("mu%i_PT>%f&&mu%i_PT<%f",index,rangesPtBothMuons_low[j],index,rangesPtBothMuons_high[j])+Nshared;            
    TString fIn = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]);  
    TString fOut = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/BinsSingleMuon/MC_D2pipimumu_%.1f_%.1f_%.1f_%.1f_mu%i_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i],rangesPtBothMuons_low[j],rangesPtBothMuons_high[j],index);
    copyMCTuplesWithCut(fIn,fOut,"BDT_Tree",ptCut);
    std::cout<<ptCut<<"on file "<<fIn<<" target is "<<fOut<<std::endl;
  }  
 }

}


void MC_PID_efficiency(double ghostProbCut,double hadronPID,double muonPID){

  
   TString cut_hadronPID_Kpimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_pipimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNpi>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_KKmumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNk>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 

  TString cut_MuonPID=TString::Format("mu0_ProbNNghost<%.1f&&mu1_ProbNNghost<%.1f&&mu0_ProbNNmu>%.1f&&mu1_ProbNNmu>%.1f",ghostProbCut,ghostProbCut,muonPID,muonPID); 
  
  TString cut_PID_Kpimumu = cut_hadronPID_Kpimumu + "&&" + cut_MuonPID;
  TString cut_PID_pipimumu = cut_hadronPID_pipimumu + "&&" + cut_MuonPID;
  TString cut_PID_KKmumu = cut_hadronPID_KKmumu + "&&" + cut_MuonPID;

  TFile* fOut= new TFile("PIDEfficiencies.root","RECREATE"); 

  TH1D* effHadronPIDKpimumu_magUp;
  TH1D* effHadronPIDKpimumu_magDw;
  TH1D* effMuonPIDKpimumu_magUp;
  TH1D* effMuonPIDKpimumu_magDw;
  TH1D* effPIDKpimumu_magUp;
  TH1D* effPIDKpimumu_magDw;

  
  effHadronPIDKpimumu_magUp= new TH1D("effHadronPIDKpimumu_magUp","efficiency hadron MC PID D->Kpimumu magUp",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effHadronPIDKpimumu_magDw= new TH1D("effHadronPIDKpimumu_magDw","efficiency hadron MC PID D->Kpimumu magDw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effMuonPIDKpimumu_magUp= new TH1D("effMuonPIDKpimumu_magUp","efficiency Muon MC PID D->Kpimumu magUp",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effMuonPIDKpimumu_magDw= new TH1D("effMuonPIDKpimumu_magDw","efficiency Muon MC PID D->Kpimumu magDw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effPIDKpimumu_magUp= new TH1D("effPIDKpimumu_magUp","efficiency MC PID D->Kpimumu magUp",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effPIDKpimumu_magDw= new TH1D("effPIDKpimumu_magDw","efficiency MC PID D->Kpimumu magDw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D* effHadronPIDpipimumu_magUp;
  TH1D* effHadronPIDpipimumu_magDw;
  TH1D* effMuonPIDpipimumu_magUp;
  TH1D* effMuonPIDpipimumu_magDw;
  TH1D* effPIDpipimumu_magUp;
  TH1D* effPIDpipimumu_magDw;

  effHadronPIDpipimumu_magUp= new TH1D("effHadronPIDpipimumu_magUp","efficiency hadron MC PID D->pipimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effHadronPIDpipimumu_magDw= new TH1D("effHadronPIDpipimumu_magDw","efficiency hadron MC PID D->pipimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effMuonPIDpipimumu_magUp= new TH1D("effMuonPIDpipimumu_magUp","efficiency Muon MC PID D->pipimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effMuonPIDpipimumu_magDw= new TH1D("effMuonPIDpipimumu_magDw","efficiency Muon MC PID D->pipimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effPIDpipimumu_magUp= new TH1D("effPIDpipimumu_magUp","efficiency MC PID D->pipimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effPIDpipimumu_magDw= new TH1D("effPIDpipimumu_magDw","efficiency MC PID D->pipimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1D* effHadronPIDKKmumu_magUp;
  TH1D* effHadronPIDKKmumu_magDw;
  TH1D* effMuonPIDKKmumu_magUp;
  TH1D* effMuonPIDKKmumu_magDw;
  TH1D* effPIDKKmumu_magUp;
  TH1D* effPIDKKmumu_magDw;

  effHadronPIDKKmumu_magUp= new TH1D("effHadronPIDKKmumu_magUp","efficiency hadron MC PID D->KKmumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  effHadronPIDKKmumu_magDw= new TH1D("effHadronPIDKKmumu_magDw","efficiency hadron MC PID D->KKmumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);
  effMuonPIDKKmumu_magUp= new TH1D("effMuonPIDKKmumu_magUp","efficiency Muon MC PID D->KKmumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  effMuonPIDKKmumu_magDw= new TH1D("effMuonPIDKKmumu_magDw","efficiency Muon MC PID D->KKmumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);
  effPIDKKmumu_magUp= new TH1D("effPIDKKmumu_magUp","efficiency MC PID D->KKmumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  effPIDKKmumu_magDw= new TH1D("effPIDKKmumu_magDw","efficiency MC PID D->KKmumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);
  
  TH1D* relEffHadronPIDKKmumu_magUp;
  TH1D* relEffHadronPIDKKmumu_magDw;
  TH1D* relEffMuonPIDKKmumu_magUp;
  TH1D* relEffMuonPIDKKmumu_magDw;
  TH1D* relEffPIDKKmumu_magUp;
  TH1D* relEffPIDKKmumu_magDw;

  relEffHadronPIDKKmumu_magUp= new TH1D("relEffHadronPIDKKmumu_magUp","relative efficiency hadron MC PID D->KKmumu/D->Kpimumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffHadronPIDKKmumu_magDw= new TH1D("relEffHadronPIDKKmumu_magDw","efficiency hadron MC PID D->KKmumu/D->Kpimumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffMuonPIDKKmumu_magUp= new TH1D("relEffMuonPIDKKmumu_magUp","efficiency Muon MC PID D->KKmumu/D->Kpimumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffMuonPIDKKmumu_magDw= new TH1D("relEffMuonPIDKKmumu_magDw","efficiency Muon MC PID D->KKmumu/D->Kpimumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffPIDKKmumu_magUp= new TH1D("relEffPIDKKmumu_magUp","efficiency MC PID D->KKmumu/D->Kpimumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffPIDKKmumu_magDw= new TH1D("relEffPIDKKmumu_magDw","efficiency MC PID D->KKmumu/D->Kpimumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1D* relEffHadronPIDpipimumu_magUp;
  TH1D* relEffHadronPIDpipimumu_magDw;
  TH1D* relEffMuonPIDpipimumu_magUp;
  TH1D* relEffMuonPIDpipimumu_magDw;
  TH1D* relEffPIDpipimumu_magUp;
  TH1D* relEffPIDpipimumu_magDw;

  relEffHadronPIDpipimumu_magUp= new TH1D("relEffHadronPIDpipimumu_magUp","relative efficiency hadron MC PID D->pipimumu/D->Kpimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffHadronPIDpipimumu_magDw= new TH1D("relEffHadronPIDpipimumu_magDw","efficiency hadron MC PID D->pipimumu/D->Kpimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffMuonPIDpipimumu_magUp= new TH1D("relEffMuonPIDpipimumu_magUp","efficiency Muon MC PID D->pipimumu/D->Kpimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffMuonPIDpipimumu_magDw= new TH1D("relEffMuonPIDpipimumu_magDw","efficiency Muon MC PID D->pipimumu/D->Kpimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffPIDpipimumu_magUp= new TH1D("relEffPIDpipimumu_magUp","efficiency MC PID D->pipimumu/D->Kpimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffPIDpipimumu_magDw= new TH1D("relEffPIDpipimumu_magDw","efficiency MC PID D->pipimumu/D->Kpimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  

  for(int i=0; i<rangesKpi_low.size();++i){
    
    TChain* tempChain_magUp= new TChain("BDT_Tree");
    tempChain_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
    TChain* tempChain_magDw= new TChain("BDT_Tree");
    tempChain_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));
    //hadron PID 
    //mag up
    double nNorm = tempChain_magUp->GetEntries();
    double  nSel = tempChain_magUp->GetEntries(cut_hadronPID_Kpimumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effHadronPIDKpimumu_magUp->SetBinContent(i+1,tempEff);
    effHadronPIDKpimumu_magUp->SetBinError(i+1,dTempEff);

    //mag dw
    nNorm = tempChain_magDw->GetEntries();
    nSel = tempChain_magDw->GetEntries(cut_hadronPID_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effHadronPIDKpimumu_magDw->SetBinContent(i+1,tempEff);
    effHadronPIDKpimumu_magDw->SetBinError(i+1,dTempEff);
    
    
    //muon PID
    nNorm = tempChain_magUp->GetEntries();
    nSel = tempChain_magUp->GetEntries(cut_MuonPID);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effMuonPIDKpimumu_magUp->SetBinContent(i+1,tempEff);
    effMuonPIDKpimumu_magUp->SetBinError(i+1,dTempEff);
    //mag dw
    nNorm = tempChain_magDw->GetEntries();
    nSel = tempChain_magDw->GetEntries(cut_MuonPID);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effMuonPIDKpimumu_magDw->SetBinContent(i+1,tempEff);
    effMuonPIDKpimumu_magDw->SetBinError(i+1,dTempEff);
     
    //total PID
    nNorm = tempChain_magUp->GetEntries();
    nSel = tempChain_magUp->GetEntries(cut_PID_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effPIDKpimumu_magUp->SetBinContent(i+1,tempEff);
    effPIDKpimumu_magUp->SetBinError(i+1,dTempEff);
    //mag dw
    nNorm = tempChain_magDw->GetEntries();
    nSel = tempChain_magDw->GetEntries(cut_PID_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effPIDKpimumu_magDw->SetBinContent(i+1,tempEff);
    effPIDKpimumu_magDw->SetBinError(i+1,dTempEff);
     
  }    

  for(int i=0; i<rangesKK_low.size();++i){
    
     TChain* tempChain_magUp= new TChain("BDT_Tree");
    tempChain_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]));
    TChain* tempChain_magDw= new TChain("BDT_Tree");
    tempChain_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]));
    //hadron PID 
    //mag up
    double nNorm = tempChain_magUp->GetEntries();
    double nSel = tempChain_magUp->GetEntries(cut_hadronPID_KKmumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effHadronPIDKKmumu_magUp->SetBinContent(i+1,tempEff);
    effHadronPIDKKmumu_magUp->SetBinError(i+1,dTempEff);
    double relEff = tempEff/effHadronPIDKKmumu_magUp->GetBinContent(1);
    double dRelEff = TMath::Sqrt( (dTempEff/effHadronPIDKKmumu_magUp->GetBinContent(1))*(dTempEff/effHadronPIDKKmumu_magUp->GetBinContent(1)) 
				  + (tempEff*effHadronPIDKKmumu_magUp->GetBinError(1)/(effHadronPIDKKmumu_magUp->GetBinContent(1)*effHadronPIDKKmumu_magUp->GetBinContent(1)) * 
				     (tempEff*effHadronPIDKKmumu_magUp->GetBinError(1)/(effHadronPIDKKmumu_magUp->GetBinContent(1)*effHadronPIDKKmumu_magUp->GetBinContent(1))
				      ) ));
    relEffHadronPIDKKmumu_magUp->SetBinContent(i+1,relEff);				     
    relEffHadronPIDKKmumu_magUp->SetBinError(i+1,dRelEff);			     
    
    //mag dw
    nNorm = tempChain_magDw->GetEntries();
    nSel = tempChain_magDw->GetEntries(cut_hadronPID_KKmumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effHadronPIDKKmumu_magDw->SetBinContent(i+1,tempEff);
    effHadronPIDKKmumu_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effHadronPIDKKmumu_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effHadronPIDKKmumu_magDw->GetBinContent(1))*(dTempEff/effHadronPIDKKmumu_magDw->GetBinContent(1)) 
				  + (tempEff*effHadronPIDKKmumu_magDw->GetBinError(1)/(effHadronPIDKKmumu_magDw->GetBinContent(1)*effHadronPIDKKmumu_magDw->GetBinContent(1)) * 
				     (tempEff*effHadronPIDKKmumu_magDw->GetBinError(1)/(effHadronPIDKKmumu_magDw->GetBinContent(1)*effHadronPIDKKmumu_magDw->GetBinContent(1))
				      ) ));
    relEffHadronPIDKKmumu_magDw->SetBinContent(i+1,relEff);				     
    relEffHadronPIDKKmumu_magDw->SetBinError(i+1,dRelEff);			     

    //muon PID
    nNorm = tempChain_magUp->GetEntries();
    nSel = tempChain_magUp->GetEntries(cut_MuonPID);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effMuonPIDKKmumu_magUp->SetBinContent(i+1,tempEff);
    effMuonPIDKKmumu_magUp->SetBinError(i+1,dTempEff);
    relEff = tempEff/effMuonPIDKKmumu_magUp->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effMuonPIDKKmumu_magUp->GetBinContent(1))*(dTempEff/effMuonPIDKKmumu_magUp->GetBinContent(1)) 
				  + (tempEff*effMuonPIDKKmumu_magUp->GetBinError(1)/(effMuonPIDKKmumu_magUp->GetBinContent(1)*effMuonPIDKKmumu_magUp->GetBinContent(1)) * 
				     (tempEff*effMuonPIDKKmumu_magUp->GetBinError(1)/(effMuonPIDKKmumu_magUp->GetBinContent(1)*effMuonPIDKKmumu_magUp->GetBinContent(1))
				      )));
    relEffMuonPIDKKmumu_magUp->SetBinContent(i+1,relEff);				     
    relEffMuonPIDKKmumu_magUp->SetBinError(i+1,dRelEff);			     
    //mag dw
    nNorm = tempChain_magDw->GetEntries();
    nSel = tempChain_magDw->GetEntries(cut_MuonPID);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effMuonPIDKKmumu_magDw->SetBinContent(i+1,tempEff);
    effMuonPIDKKmumu_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effMuonPIDKKmumu_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effMuonPIDKKmumu_magDw->GetBinContent(1))*(dTempEff/effMuonPIDKKmumu_magDw->GetBinContent(1)) 
				  + (tempEff*effMuonPIDKKmumu_magDw->GetBinError(1)/(effMuonPIDKKmumu_magDw->GetBinContent(1)*effMuonPIDKKmumu_magDw->GetBinContent(1)) * 
				     (tempEff*effMuonPIDKKmumu_magDw->GetBinError(1)/(effMuonPIDKKmumu_magDw->GetBinContent(1)*effMuonPIDKKmumu_magDw->GetBinContent(1))
				      )));
    relEffMuonPIDKKmumu_magDw->SetBinContent(i+1,relEff);				     
    relEffMuonPIDKKmumu_magDw->SetBinError(i+1,dRelEff);			     
     
    //total PID
    nNorm = tempChain_magUp->GetEntries();
    nSel = tempChain_magUp->GetEntries(cut_PID_KKmumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effPIDKKmumu_magUp->SetBinContent(i+1,tempEff);
    effPIDKKmumu_magUp->SetBinError(i+1,dTempEff);
    relEff = tempEff/effPIDKKmumu_magUp->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effPIDKKmumu_magUp->GetBinContent(1))*(dTempEff/effPIDKKmumu_magUp->GetBinContent(1)) 
				  + (tempEff*effPIDKKmumu_magUp->GetBinError(1)/(effPIDKKmumu_magUp->GetBinContent(1)*effPIDKKmumu_magUp->GetBinContent(1)) * 
				     (tempEff*effPIDKKmumu_magUp->GetBinError(1)/(effPIDKKmumu_magUp->GetBinContent(1)*effPIDKKmumu_magUp->GetBinContent(1))
				      )));
    relEffPIDKKmumu_magUp->SetBinContent(i+1,relEff);				     
    relEffPIDKKmumu_magUp->SetBinError(i+1,dRelEff);			     
    //mag dw
    nNorm = tempChain_magDw->GetEntries();
    nSel = tempChain_magDw->GetEntries(cut_PID_KKmumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effPIDKKmumu_magDw->SetBinContent(i+1,tempEff);
    effPIDKKmumu_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effPIDKKmumu_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effPIDKKmumu_magDw->GetBinContent(1))*(dTempEff/effPIDKKmumu_magDw->GetBinContent(1)) 
				  + (tempEff*effPIDKKmumu_magDw->GetBinError(1)/(effPIDKKmumu_magDw->GetBinContent(1)*effPIDKKmumu_magDw->GetBinContent(1)) * 
				     (tempEff*effPIDKKmumu_magDw->GetBinError(1)/(effPIDKKmumu_magDw->GetBinContent(1)*effPIDKKmumu_magDw->GetBinContent(1))
				      )));
    relEffPIDKKmumu_magDw->SetBinContent(i+1,relEff);				     
    relEffPIDKKmumu_magDw->SetBinError(i+1,dRelEff);			     
      

  }    

  for(int i=0; i<rangespipi_low.size();++i){
    
     TChain* tempChain_magUp= new TChain("BDT_Tree");
    tempChain_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]));
    TChain* tempChain_magDw= new TChain("BDT_Tree");
    tempChain_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]));
    //hadron PID 
    //mag up
    double nNorm = tempChain_magUp->GetEntries();
    double nSel = tempChain_magUp->GetEntries(cut_hadronPID_pipimumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effHadronPIDpipimumu_magUp->SetBinContent(i+1,tempEff);
    effHadronPIDpipimumu_magUp->SetBinError(i+1,dTempEff);
    double relEff = tempEff/effHadronPIDpipimumu_magUp->GetBinContent(1);
    double dRelEff = TMath::Sqrt( (dTempEff/effHadronPIDpipimumu_magUp->GetBinContent(1))*(dTempEff/effHadronPIDpipimumu_magUp->GetBinContent(1)) 
				  + (tempEff*effHadronPIDpipimumu_magUp->GetBinError(1)/(effHadronPIDpipimumu_magUp->GetBinContent(1)*effHadronPIDpipimumu_magUp->GetBinContent(1)) * 
				     (tempEff*effHadronPIDpipimumu_magUp->GetBinError(1)/(effHadronPIDpipimumu_magUp->GetBinContent(1)*effHadronPIDpipimumu_magUp->GetBinContent(1))
				      )));
    
    relEffHadronPIDpipimumu_magUp->SetBinContent(i+1,relEff);				     
    relEffHadronPIDpipimumu_magUp->SetBinError(i+1,dRelEff);			     

    std::cout<<nNorm <<"  "<<nSel <<"  "<<tempEff<<"  "<<dTempEff<<"  "<<std::endl; 
    std::cout<<cut_hadronPID_Kpimumu <<std::endl; 

    //mag dw
    nNorm = tempChain_magDw->GetEntries();
    nSel = tempChain_magDw->GetEntries(cut_hadronPID_pipimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effHadronPIDpipimumu_magDw->SetBinContent(i+1,tempEff);
    effHadronPIDpipimumu_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effHadronPIDpipimumu_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effHadronPIDpipimumu_magDw->GetBinContent(1))*(dTempEff/effHadronPIDpipimumu_magDw->GetBinContent(1)) 
				  + (tempEff*effHadronPIDpipimumu_magDw->GetBinError(1)/(effHadronPIDpipimumu_magDw->GetBinContent(1)*effHadronPIDpipimumu_magDw->GetBinContent(1)) * 
				     (tempEff*effHadronPIDpipimumu_magDw->GetBinError(1)/(effHadronPIDpipimumu_magDw->GetBinContent(1)*effHadronPIDpipimumu_magDw->GetBinContent(1))
				      )));
    relEffHadronPIDpipimumu_magDw->SetBinContent(i+1,relEff);				     
    relEffHadronPIDpipimumu_magDw->SetBinError(i+1,dRelEff);			     

    //muon PID
    nNorm = tempChain_magUp->GetEntries();
    nSel = tempChain_magUp->GetEntries(cut_MuonPID);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effMuonPIDpipimumu_magUp->SetBinContent(i+1,tempEff);
    effMuonPIDpipimumu_magUp->SetBinError(i+1,dTempEff);
    relEff = tempEff/effMuonPIDpipimumu_magUp->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effMuonPIDpipimumu_magUp->GetBinContent(1))*(dTempEff/effMuonPIDpipimumu_magUp->GetBinContent(1)) 
				  + (tempEff*effMuonPIDpipimumu_magUp->GetBinError(1)/(effMuonPIDpipimumu_magUp->GetBinContent(1)*effMuonPIDpipimumu_magUp->GetBinContent(1)) * 
				     (tempEff*effMuonPIDpipimumu_magUp->GetBinError(1)/(effMuonPIDpipimumu_magUp->GetBinContent(1)*effMuonPIDpipimumu_magUp->GetBinContent(1))
				      )));
    relEffMuonPIDpipimumu_magUp->SetBinContent(i+1,relEff);				     
    relEffMuonPIDpipimumu_magUp->SetBinError(i+1,dRelEff);			     
    //mag dw
    nNorm = tempChain_magDw->GetEntries();
    nSel = tempChain_magDw->GetEntries(cut_MuonPID);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effMuonPIDpipimumu_magDw->SetBinContent(i+1,tempEff);
    effMuonPIDpipimumu_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effMuonPIDpipimumu_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effMuonPIDpipimumu_magDw->GetBinContent(1))*(dTempEff/effMuonPIDpipimumu_magDw->GetBinContent(1)) 
				  + (tempEff*effMuonPIDpipimumu_magDw->GetBinError(1)/(effMuonPIDpipimumu_magDw->GetBinContent(1)*effMuonPIDpipimumu_magDw->GetBinContent(1)) * 
				     (tempEff*effMuonPIDpipimumu_magDw->GetBinError(1)/(effMuonPIDpipimumu_magDw->GetBinContent(1)*effMuonPIDpipimumu_magDw->GetBinContent(1))
				      )));
    relEffMuonPIDpipimumu_magDw->SetBinContent(i+1,relEff);				     
    relEffMuonPIDpipimumu_magDw->SetBinError(i+1,dRelEff);			     
     
    //total PID
    nNorm = tempChain_magUp->GetEntries();
    nSel = tempChain_magUp->GetEntries(cut_PID_pipimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effPIDpipimumu_magUp->SetBinContent(i+1,tempEff);
    effPIDpipimumu_magUp->SetBinError(i+1,dTempEff);
    relEff = tempEff/effPIDpipimumu_magUp->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effPIDpipimumu_magUp->GetBinContent(1))*(dTempEff/effPIDpipimumu_magUp->GetBinContent(1)) 
				  + (tempEff*effPIDpipimumu_magUp->GetBinError(1)/(effPIDpipimumu_magUp->GetBinContent(1)*effPIDpipimumu_magUp->GetBinContent(1)) * 
				     (tempEff*effPIDpipimumu_magUp->GetBinError(1)/(effPIDpipimumu_magUp->GetBinContent(1)*effPIDpipimumu_magUp->GetBinContent(1))
				      )));
    relEffPIDpipimumu_magUp->SetBinContent(i+1,relEff);				     
    relEffPIDpipimumu_magUp->SetBinError(i+1,dRelEff);			     
    //mag dw
    nNorm = tempChain_magDw->GetEntries();
    nSel = tempChain_magDw->GetEntries(cut_PID_pipimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effPIDpipimumu_magDw->SetBinContent(i+1,tempEff);
    effPIDpipimumu_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effPIDpipimumu_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effPIDpipimumu_magDw->GetBinContent(1))*(dTempEff/effPIDpipimumu_magDw->GetBinContent(1)) 
				  + (tempEff*effPIDpipimumu_magDw->GetBinError(1)/(effPIDpipimumu_magDw->GetBinContent(1)*effPIDpipimumu_magDw->GetBinContent(1)) * 
				     (tempEff*effPIDpipimumu_magDw->GetBinError(1)/(effPIDpipimumu_magDw->GetBinContent(1)*effPIDpipimumu_magDw->GetBinContent(1)) 
				      )));
    relEffPIDpipimumu_magDw->SetBinContent(i+1,relEff);				     
    relEffPIDpipimumu_magDw->SetBinError(i+1,dRelEff);			     
				     

  }    

  
  //Drawing part 

  TCanvas *a  = new TCanvas("pipimumu","Eff pipimumu");
  a->Divide(2,2);
  a->cd(1);
  effHadronPIDpipimumu_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effHadronPIDpipimumu_magUp->Draw();
  effHadronPIDpipimumu_magDw->SetLineColor(kRed);
  effHadronPIDpipimumu_magDw->Draw("SAME");

  a->cd(2);
  effMuonPIDpipimumu_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effMuonPIDpipimumu_magUp->Draw();
  effMuonPIDpipimumu_magDw->SetLineColor(kRed);
  effMuonPIDpipimumu_magDw->Draw("SAME");

  a->cd(3);
  effPIDpipimumu_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effPIDpipimumu_magUp->Draw();
  effPIDpipimumu_magDw->SetLineColor(kRed);
  effPIDpipimumu_magDw->Draw("SAME");

  a->cd(4);
  relEffHadronPIDpipimumu_magUp->GetYaxis()->SetRangeUser(0.5,1.3);
  relEffHadronPIDpipimumu_magUp->Draw();
  relEffMuonPIDpipimumu_magUp->SetLineColor(kRed);
  relEffMuonPIDpipimumu_magUp->Draw("SAME");
  relEffPIDpipimumu_magUp->SetLineColor(kViolet);
  relEffPIDpipimumu_magUp->Draw("SAME");

  a->Write();


  TCanvas *b  = new TCanvas("KKmumu","Eff KKmumu");
  b->Divide(2,2);
  b->cd(1);
  effHadronPIDKKmumu_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effHadronPIDKKmumu_magUp->Draw();
  effHadronPIDKKmumu_magDw->SetLineColor(kRed);
  effHadronPIDKKmumu_magDw->Draw("SAME");

  b->cd(2);
  effMuonPIDKKmumu_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effMuonPIDKKmumu_magUp->Draw();
  effMuonPIDKKmumu_magDw->SetLineColor(kRed);
  effMuonPIDKKmumu_magDw->Draw("SAME");

  b->cd(3);
  effPIDKKmumu_magUp->GetYaxis()->SetRangeUser(0.3,1);
  effPIDKKmumu_magUp->Draw();
  effPIDKKmumu_magDw->SetLineColor(kRed);
  effPIDKKmumu_magDw->Draw("SAME");

  b->cd(4);
  relEffHadronPIDKKmumu_magUp->GetYaxis()->SetRangeUser(0.5,1.3);
  relEffHadronPIDKKmumu_magUp->Draw();
  relEffMuonPIDKKmumu_magUp->SetLineColor(kRed);
  relEffMuonPIDKKmumu_magUp->Draw("SAME");
  relEffPIDKKmumu_magUp->SetLineColor(kViolet);
  relEffPIDKKmumu_magUp->Draw("SAME");

  b->Write();
  
  TCanvas *c  = new TCanvas("Kpimumu","Eff Kpimumu");
  c->Divide(2,2);
  c->cd(1);
  effHadronPIDKpimumu_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effHadronPIDKpimumu_magUp->Draw();
  effHadronPIDKpimumu_magDw->SetLineColor(kRed);
  effHadronPIDKpimumu_magDw->Draw("SAME");

  c->cd(2);
  effMuonPIDKpimumu_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effMuonPIDKpimumu_magUp->Draw();
  effMuonPIDKpimumu_magDw->SetLineColor(kRed);
  effMuonPIDKpimumu_magDw->Draw("SAME");

  c->cd(3);
  effPIDKpimumu_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effPIDKpimumu_magUp->Draw();
  effPIDKpimumu_magDw->SetLineColor(kRed);
  effPIDKpimumu_magDw->Draw("SAME");

  a->Print("test1.eps");
  c->Write();
  
}
			   
void MC_BDT_efficiency(double ghostProbCut,double hadronPID,double muonPID,double BDT,bool applyPID, bool applyTrigger){

  dcastyle;
  
  TString cutNshared="mu0_MuonNShared==0&&mu1_MuonNShared==0";
  TString cut_hadronPID_Kpimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_pipimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNpi>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_KKmumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNk>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_SlowpiPID = TString::Format("Slowpi_ProbNNghost<%.1f",ghostProbCut); 
  
  TString cutTrigger_KKmumu = "(mu0_L0MuonDecision_TOS==1 || mu1_L0MuonDecision_TOS ==1) && (mu0_Hlt1TrackMuonDecision_TOS==1 || mu1_Hlt1TrackMuonDecision_TOS==1 || D_Hlt1TrackAllL0Decision_TOS==1)&&D_Hlt2CharmSemilepD02KKMuMuDecision_TOS==1";
  TString cutTrigger_Kpimumu = "(mu0_L0MuonDecision_TOS==1 || mu1_L0MuonDecision_TOS ==1) && (mu0_Hlt1TrackMuonDecision_TOS==1 || mu1_Hlt1TrackMuonDecision_TOS==1 || D_Hlt1TrackAllL0Decision_TOS==1)&&D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS==1";
  TString cutTrigger_pipimumu = "(mu0_L0MuonDecision_TOS==1 || mu1_L0MuonDecision_TOS ==1) && (mu0_Hlt1TrackMuonDecision_TOS==1 || mu1_Hlt1TrackMuonDecision_TOS==1 || D_Hlt1TrackAllL0Decision_TOS==1)&&D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS==1";


  TString cut_MuonPID=TString::Format("mu0_ProbNNghost<%.1f&&mu1_ProbNNghost<%.1f&&mu0_ProbNNmu>%.1f&&mu1_ProbNNmu>%.1f",ghostProbCut,ghostProbCut,muonPID,muonPID); 
 
  TString cut_BDT=TString::Format("BDT>%.1f",BDT);

  TString cut_PID_Kpimumu; 
  TString cut_PID_pipimumu;
  TString cut_PID_KKmumu ;

  TString cut_Kpimumu ;
  TString cut_pipimumu;
  TString cut_KKmumu ; 

  TString nameTarget ;
  nameTarget="/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/";
 
  //cut *PID* will the one be aplied to the normalization yield, cut_*mumu is the actual selection cut

  //if(applyTrigger && !applyPID) cout<<"trigger but no PID cut is not implemented right now. Either PID, PID and Trigger or neither of them."<<endl;

  
  if(applyPID && applyTrigger) {
    
    cut_PID_Kpimumu = cut_hadronPID_Kpimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID + "&&" + cutNshared  + "&&" + cutTrigger_Kpimumu;
    cut_PID_pipimumu = cut_hadronPID_pipimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID+ "&&" + cutNshared + "&&" + cutTrigger_pipimumu;
    cut_PID_KKmumu = cut_hadronPID_KKmumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID+ "&&" + cutNshared + "&&" + cutTrigger_KKmumu;
    
    nameTarget+= "MCBDTEfficiencies_afterPID_andTrigger.root";

  }

  if(applyPID && !applyTrigger) {
    
    cut_PID_Kpimumu = cut_hadronPID_Kpimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID + "&&" + cutNshared ; 
    cut_PID_pipimumu = cut_hadronPID_pipimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID+ "&&" + cutNshared; 
    cut_PID_KKmumu = cut_hadronPID_KKmumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID+ "&&" + cutNshared ;
        
    nameTarget+= "MCBDTEfficiencies_afterPID.root";
  }

    if(!applyPID && !applyTrigger) {

    cut_PID_Kpimumu = cutNshared ; 
    cut_PID_pipimumu = cutNshared; 
    cut_PID_KKmumu =  cutNshared ;
        
    nameTarget+= "MCBDTEfficiencies.root";
  }

  //build the selection string
  cut_Kpimumu = cut_PID_Kpimumu + "&&" + cut_BDT;
  cut_pipimumu =cut_PID_pipimumu + "&&" + cut_BDT;
  cut_KKmumu = cut_PID_KKmumu + "&&" + cut_BDT;
   

  TFile* fOut= new TFile(nameTarget,"RECREATE"); 

  TH1D* effBDTKpimumu_KKmumuTrained_magUp;
  TH1D* effBDTKpimumu_KKmumuTrained_magDw;
  TH1D* effBDTKpimumu_KKmumuTrained;

  effBDTKpimumu_KKmumuTrained_magUp= new TH1D("effBDTKpimumu_KKmumuTrained_magUp","efficiency MC BDT D->Kpimumu magUp",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effBDTKpimumu_KKmumuTrained_magDw= new TH1D("effBDTKpimumu_KKmumuTrained_magDw","efficiency MC BDT D->Kpimumu magDw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effBDTKpimumu_KKmumuTrained= new TH1D("effBDTKpimumu_KKmumuTrained","efficiency MC BDT D->Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D* effBDTKpimumu_pipimumuTrained_magUp;
  TH1D* effBDTKpimumu_pipimumuTrained_magDw; 
  TH1D* effBDTKpimumu_pipimumuTrained; 
  TH1D* effBDTKpimumu_pipimumuTrained_data;

  effBDTKpimumu_pipimumuTrained_magUp= new TH1D("effBDTKpimumu_pipimumuTrained_magUp","efficiency MC BDT D->Kpimumu magUp",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effBDTKpimumu_pipimumuTrained_magDw= new TH1D("effBDTKpimumu_pipimumuTrained_magDw","efficiency MC BDT D->Kpimumu magDw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effBDTKpimumu_pipimumuTrained= new TH1D("effBDTKpimumu_pipimumuTrained","efficiency MC BDT D->Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effBDTKpimumu_pipimumuTrained_data= new TH1D("effBDTKpimumu_pipimumuTrained_data","efficiency MC BDT D->Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
 
  TH1D* effBDTKKmumu_KKmumuTrained_magUp;
  TH1D* effBDTKKmumu_KKmumuTrained_magDw;
  TH1D* effBDTKKmumu_KKmumuTrained;

  effBDTKKmumu_KKmumuTrained_magUp= new TH1D("effBDTKKmumu_KKmumuTrained_magUp","efficiency MC BDT D->KKmumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  effBDTKKmumu_KKmumuTrained_magDw= new TH1D("effBDTKKmumu_KKmumuTrained_magDw","efficiency MC BDT D->KKmumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);
  effBDTKKmumu_KKmumuTrained= new TH1D("effBDTKKmumu_KKmumuTrained","efficiency MC BDT D->KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1D* effBDTpipimumu_pipimumuTrained_magUp;
  TH1D* effBDTpipimumu_pipimumuTrained_magDw;
  TH1D* effBDTpipimumu_pipimumuTrained;
  TH1D* effBDTpipimumu_pipimumuTrained_data;

  effBDTpipimumu_pipimumuTrained_magUp= new TH1D("effBDTpipimumu_pipimumuTrained_magUp","efficiency MC BDT D->pipimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effBDTpipimumu_pipimumuTrained_magDw= new TH1D("effBDTKpipimumu_pipimumuTrained_magDw","efficiency MC BDT D->pipimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effBDTpipimumu_pipimumuTrained= new TH1D("effBDTKpipimumu_pipimumuTrained","efficiency MC BDT D->pipimumu ",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effBDTpipimumu_pipimumuTrained_data= new TH1D("effBDTpipipimumu_pipimumuTrained_data","efficiency MC BDT D->pipimumu ",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1D* relEffBDTKKmumu_magUp;
  TH1D* relEffBDTKKmumu_magDw;
  TH1D* relEffBDTKKmumu;
  TH1D* relEffBDTpipimumu_magUp;
  TH1D* relEffBDTpipimumu_magDw;
  TH1D* relEffBDTpipimumu;
  TH1D* relEffBDTpipimumu_data;
  

  relEffBDTKKmumu_magUp= new TH1D("relEffBDTDKKmumu_magUp","relative efficiency MC BDT D->KKmumu/D->Kpimumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffBDTKKmumu_magDw= new TH1D("relEffBDTDKKmumu_magDw","relative efficiency MC BDT D->KKmumu/D->Kpimumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffBDTKKmumu= new TH1D("relEffBDTDKKmumu","relative efficiency MC BDT D->KKmumu/D->Kpimumu ",sizeof(binsKK)/sizeof(double)-1,binsKK);

  relEffBDTpipimumu_magUp= new TH1D("relEffBDTDpipimumu_magUp","relative efficiency MC BDT D->pipimumu/D->Kpimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi); 
  relEffBDTpipimumu_magDw= new TH1D("relEffBDTDpipimumu_magDw","relative efficiency MC BDT D->pipimumu/D->Kpimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffBDTpipimumu= new TH1D("relEffBDTDpipimumu","relative efficiency MC BDT D->pipimumu/D->Kpimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffBDTpipimumu_data= new TH1D("relEffBDTDpipimumu_data","relative efficiency MC BDT D->pipimumu/D->Kpimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);
 

  for(int i=0; i<rangesKpi_low.size();++i){
    
    TChain* tempChain_KKmumuTrained_magUp= new TChain("BDT_Tree");
    tempChain_KKmumuTrained_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
    TChain* tempChain_KKmumuTrained_magDw= new TChain("BDT_Tree");
    tempChain_KKmumuTrained_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));

    TChain* tempChain_pipimumuTrained_magUp= new TChain("BDT_Tree");
    tempChain_pipimumuTrained_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
    TChain* tempChain_pipimumuTrained_magDw= new TChain("BDT_Tree");
    tempChain_pipimumuTrained_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));


    //mag up KKmumu trained
    double nNorm = (double)tempChain_KKmumuTrained_magUp->GetEntries(cut_PID_Kpimumu);
    double  nSel =(double)tempChain_KKmumuTrained_magUp->GetEntries(cut_Kpimumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKpimumu_KKmumuTrained_magUp->SetBinContent(i+1,tempEff);
    effBDTKpimumu_KKmumuTrained_magUp->SetBinError(i+1,dTempEff);

    //mag Dw KKmumu trained
    nNorm = (double)tempChain_KKmumuTrained_magDw->GetEntries(cut_PID_Kpimumu);
    nSel = (double)tempChain_KKmumuTrained_magDw->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKpimumu_KKmumuTrained_magDw->SetBinContent(i+1,tempEff);
    effBDTKpimumu_KKmumuTrained_magDw->SetBinError(i+1,dTempEff);

    //both polarities
    //mag up KKmumu trained                                                                                                                                                                  
    //mag Dw KKmumu trained
    nNorm = (double)tempChain_KKmumuTrained_magDw->GetEntries(cut_PID_Kpimumu) + tempChain_KKmumuTrained_magUp->GetEntries(cut_PID_Kpimumu);
    nSel = (double)tempChain_KKmumuTrained_magDw->GetEntries(cut_Kpimumu) + tempChain_KKmumuTrained_magUp->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKpimumu_KKmumuTrained->SetBinContent(i+1,tempEff);
    effBDTKpimumu_KKmumuTrained->SetBinError(i+1,dTempEff);

    //mag up pipimumu trained
    nNorm = (double)tempChain_pipimumuTrained_magUp->GetEntries(cut_PID_Kpimumu);
    nSel = (double)tempChain_pipimumuTrained_magUp->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKpimumu_pipimumuTrained_magUp->SetBinContent(i+1,tempEff);
    effBDTKpimumu_pipimumuTrained_magUp->SetBinError(i+1,dTempEff);

    //mag Dw pipimumu trained
    nNorm = (double)tempChain_pipimumuTrained_magDw->GetEntries(cut_PID_Kpimumu);
    nSel = (double)tempChain_pipimumuTrained_magDw->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKpimumu_pipimumuTrained_magDw->SetBinContent(i+1,tempEff);
    effBDTKpimumu_pipimumuTrained_magDw->SetBinError(i+1,dTempEff);


    //both polarities
    //mag up KKmumu trained                                                                                                                                                                  
    //mag Dw KKmumu trained
    nNorm = (double)tempChain_pipimumuTrained_magDw->GetEntries(cut_PID_Kpimumu) + (double)tempChain_pipimumuTrained_magUp->GetEntries(cut_PID_Kpimumu);
    nSel = (double)tempChain_pipimumuTrained_magDw->GetEntries(cut_Kpimumu) + (double)tempChain_pipimumuTrained_magUp->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKpimumu_pipimumuTrained->SetBinContent(i+1,tempEff);
    effBDTKpimumu_pipimumuTrained->SetBinError(i+1,dTempEff);

     
  }    

  for(int i=0; i<rangesKK_low.size();++i){
    
     TChain* tempChain_magUp= new TChain("BDT_Tree");
    tempChain_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]));
    TChain* tempChain_magDw= new TChain("BDT_Tree");
    tempChain_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]));
    //hadron PID 
    //mag up
    double nNorm = (double)tempChain_magUp->GetEntries(cut_PID_KKmumu);
    double nSel = (double)tempChain_magUp->GetEntries(cut_KKmumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    //cout<<"bin "<<i<<" "<<nNorm<<"  "<<nSel<<"  "<<tempEff<<"+-"<<dTempEff<<cut_PID_KKmumu<<" "<<cut_KKmumu<< endl;
    effBDTKKmumu_KKmumuTrained_magUp->SetBinContent(i+1,tempEff);
    effBDTKKmumu_KKmumuTrained_magUp->SetBinError(i+1,dTempEff);
    double relEff = tempEff/effBDTKpimumu_KKmumuTrained_magUp->GetBinContent(1);
    double dRelEff = TMath::Sqrt( (dTempEff/effBDTKKmumu_KKmumuTrained_magUp->GetBinContent(1))*(dTempEff/effBDTKKmumu_KKmumuTrained_magUp->GetBinContent(1)) 
				  + (tempEff*effBDTKKmumu_KKmumuTrained_magUp->GetBinError(1)/(effBDTKKmumu_KKmumuTrained_magUp->GetBinContent(1)*effBDTKKmumu_KKmumuTrained_magUp->GetBinContent(1)) * 
				     (tempEff*effBDTKKmumu_KKmumuTrained_magUp->GetBinError(1)/(effBDTKKmumu_KKmumuTrained_magUp->GetBinContent(1)*effBDTKKmumu_KKmumuTrained_magUp->GetBinContent(1))
				      ) ));
    relEffBDTKKmumu_magUp->SetBinContent(i+1,relEff);				     
    relEffBDTKKmumu_magUp->SetBinError(i+1,dRelEff);			     
    
    //mag Dw
    nNorm =(double) tempChain_magDw->GetEntries(cut_PID_KKmumu);
    nSel =(double) tempChain_magDw->GetEntries(cut_KKmumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKKmumu_KKmumuTrained_magDw->SetBinContent(i+1,tempEff);
    effBDTKKmumu_KKmumuTrained_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effBDTKpimumu_KKmumuTrained_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effBDTKKmumu_KKmumuTrained_magDw->GetBinContent(1))*(dTempEff/effBDTKKmumu_KKmumuTrained_magDw->GetBinContent(1)) 
				  + (tempEff*effBDTKKmumu_KKmumuTrained_magDw->GetBinError(1)/(effBDTKKmumu_KKmumuTrained_magDw->GetBinContent(1)*effBDTKKmumu_KKmumuTrained_magDw->GetBinContent(1)) * (tempEff*effBDTKKmumu_KKmumuTrained_magDw->GetBinError(1)/(effBDTKKmumu_KKmumuTrained_magDw->GetBinContent(1)*effBDTKKmumu_KKmumuTrained_magDw->GetBinContent(1))
				      ) ));
    relEffBDTKKmumu_magDw->SetBinContent(i+1,relEff);				     
    relEffBDTKKmumu_magDw->SetBinError(i+1,dRelEff);			     
      
    //both polarities
    nNorm = (double)tempChain_magDw->GetEntries(cut_PID_KKmumu)+(double)tempChain_magUp->GetEntries(cut_PID_KKmumu);
    nSel = (double)tempChain_magDw->GetEntries(cut_KKmumu)+(double)tempChain_magUp->GetEntries(cut_KKmumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKKmumu_KKmumuTrained->SetBinContent(i+1,tempEff);
    effBDTKKmumu_KKmumuTrained->SetBinError(i+1,dTempEff);
    relEff = tempEff/effBDTKpimumu_KKmumuTrained->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effBDTKKmumu_KKmumuTrained->GetBinContent(1))*(dTempEff/effBDTKKmumu_KKmumuTrained->GetBinContent(1)) 
				  + (tempEff*effBDTKKmumu_KKmumuTrained->GetBinError(1)/(effBDTKKmumu_KKmumuTrained->GetBinContent(1)*effBDTKKmumu_KKmumuTrained->GetBinContent(1)) * (tempEff*effBDTKKmumu_KKmumuTrained->GetBinError(1)/(effBDTKKmumu_KKmumuTrained->GetBinContent(1)*effBDTKKmumu_KKmumuTrained->GetBinContent(1))
				      ) ));
    relEffBDTKKmumu->SetBinContent(i+1,relEff);				     
    relEffBDTKKmumu->SetBinError(i+1,dRelEff);			     
    
   }    

  for(int i=0; i<rangespipi_low.size();++i){
    
     TChain* tempChain_magUp= new TChain("BDT_Tree");
    tempChain_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]));
    TChain* tempChain_magDw= new TChain("BDT_Tree");
    tempChain_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]));
    //hadron PID 
    //mag up
    double nNorm = (double)tempChain_magUp->GetEntries(cut_PID_pipimumu);
    double nSel = (double)tempChain_magUp->GetEntries(cut_pipimumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTpipimumu_pipimumuTrained_magUp->SetBinContent(i+1,tempEff);
    effBDTpipimumu_pipimumuTrained_magUp->SetBinError(i+1,dTempEff);
    double relEff = tempEff/effBDTKpimumu_pipimumuTrained_magUp->GetBinContent(1);
    double dRelEff = TMath::Sqrt( (dTempEff/effBDTpipimumu_pipimumuTrained_magUp->GetBinContent(1))*(dTempEff/effBDTpipimumu_pipimumuTrained_magUp->GetBinContent(1)) 
				  + (tempEff*effBDTpipimumu_pipimumuTrained_magUp->GetBinError(1)/(effBDTpipimumu_pipimumuTrained_magUp->GetBinContent(1)*effBDTpipimumu_pipimumuTrained_magUp->GetBinContent(1)) * 
				     (tempEff*effBDTpipimumu_pipimumuTrained_magUp->GetBinError(1)/(effBDTpipimumu_pipimumuTrained_magUp->GetBinContent(1)*effBDTpipimumu_pipimumuTrained_magUp->GetBinContent(1))
				      ) ));
    relEffBDTpipimumu_magUp->SetBinContent(i+1,relEff);				     
    relEffBDTpipimumu_magUp->SetBinError(i+1,dRelEff);			     
    
    //mag Dw
    nNorm = (double)tempChain_magDw->GetEntries(cut_PID_pipimumu);
    nSel = (double)tempChain_magDw->GetEntries(cut_pipimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTpipimumu_pipimumuTrained_magDw->SetBinContent(i+1,tempEff);
    effBDTpipimumu_pipimumuTrained_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effBDTKpimumu_pipimumuTrained_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effBDTpipimumu_pipimumuTrained_magDw->GetBinContent(1))*(dTempEff/effBDTpipimumu_pipimumuTrained_magDw->GetBinContent(1)) 
				  + (tempEff*effBDTpipimumu_pipimumuTrained_magDw->GetBinError(1)/(effBDTpipimumu_pipimumuTrained_magDw->GetBinContent(1)*effBDTpipimumu_pipimumuTrained_magDw->GetBinContent(1)) * 
				     (tempEff*effBDTpipimumu_pipimumuTrained_magDw->GetBinError(1)/(effBDTpipimumu_pipimumuTrained_magDw->GetBinContent(1)*effBDTpipimumu_pipimumuTrained_magDw->GetBinContent(1))
				      ) ));
    relEffBDTpipimumu_magDw->SetBinContent(i+1,relEff);				     
    relEffBDTpipimumu_magDw->SetBinError(i+1,dRelEff);			     
      
   //both polarities
    nNorm = (double)tempChain_magDw->GetEntries(cut_PID_pipimumu)+(double)tempChain_magUp->GetEntries(cut_PID_pipimumu);
    nSel = (double)tempChain_magDw->GetEntries(cut_pipimumu)+(double)tempChain_magUp->GetEntries(cut_pipimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTpipimumu_pipimumuTrained->SetBinContent(i+1,tempEff);
    effBDTpipimumu_pipimumuTrained->SetBinError(i+1,dTempEff);
    relEff = tempEff/effBDTKpimumu_pipimumuTrained->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effBDTpipimumu_pipimumuTrained->GetBinContent(1))*(dTempEff/effBDTpipimumu_pipimumuTrained->GetBinContent(1)) 
				  + (tempEff*effBDTpipimumu_pipimumuTrained->GetBinError(1)/(effBDTpipimumu_pipimumuTrained->GetBinContent(1)*effBDTpipimumu_pipimumuTrained->GetBinContent(1)) * (tempEff*effBDTpipimumu_pipimumuTrained->GetBinError(1)/(effBDTpipimumu_pipimumuTrained->GetBinContent(1)*effBDTpipimumu_pipimumuTrained->GetBinContent(1))
				      ) ));
    relEffBDTpipimumu->SetBinContent(i+1,relEff);				     
    relEffBDTpipimumu->SetBinError(i+1,dRelEff);			     
 
  }    
  

  ////////////////////////////////////////////////
  ////////// data part

  //ghost prob and nshared cut already aplied to the *noCuts* files
  TString kind ="D2pipimumu";
  TString dataCut = "mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&BDT>0";
  TString dataCutNoBDT="mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5";
  TString misIDCut = "mu0_ProbNNmu>0.5&&BDT>0";
  TString q2Cut= "D_DiMuon_Mass>950&&D_DiMuon_Mass<1100";
  TString q2RangeNormalizationMode="D_DiMuon_Mass>675&&D_DiMuon_Mass<875";

  cut_MuonPID=TString::Format("mu0_ProbNNghost<%.1f&&mu1_ProbNNghost<%.1f&&mu0_ProbNNmu>%.1f&&mu1_ProbNNmu>%.1f",ghostProbCut,ghostProbCut,muonPID,muonPID);

  cut_PID_Kpimumu = cut_hadronPID_Kpimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID + "&&" + cutNshared;
  cut_PID_pipimumu = cut_hadronPID_pipimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID+ "&&" + cutNshared;
  cut_PID_KKmumu = cut_hadronPID_KKmumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID+ "&&" + cutNshared;

  cut_Kpimumu = cut_PID_Kpimumu + "&&" + cut_BDT;
  cut_pipimumu =cut_PID_pipimumu + "&&" + cut_BDT;
   
  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2pipimumuBDT.root");
  myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2pipimumuBDT.root");
//  myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2Kpimumu_D2pipimumuBDT.root");
  myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipimumu_BDT.root");


  //didnt create extra files for misID and MC without preselection
  myFitter1D.fit_MC(dataCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_"+dataCut+q2Cut+".eps");
  myFitter1D.fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_normMode_"+dataCut+".eps");
  std::cout<<"Monte Carlo fits done.."<<std::endl;

  myFitter1D.fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/misID_Kpi"+misIDCut+q2RangeNormalizationMode+".eps");
  myFitter1D.fit_HHpipi_misID(misIDCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/misID_pipi"+misIDCut+q2Cut+".eps"); ////q2 cut missing? check!        
  std::cout<<"misID fits done.."<<std::endl;
  /*
  double nSel_norm = myFitter1D.fit_normalization_Data(dataCut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/normMode_"+dataCut+q2RangeNormalizationMode+".eps");
  double nSel_pipi=myFitter1D.fit_resonant_Data(dataCut+"&&"+q2Cut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/resonant_data_"+dataCut+q2Cut+".eps");

  double nTot_norm = myFitter1D.fit_normalization_Data(dataCutNoBDT,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/normMode_"+dataCutNoBDT+q2RangeNormalizationMode+".eps");
  double nTot_pipi=myFitter1D.fit_resonant_Data(dataCutNoBDT+"&&"+q2Cut,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/resonant_data_"+dataCutNoBDT+q2Cut+".eps");
  */

  
  double nSel_norm = myFitter1D.fit_normalization_Data(cut_Kpimumu,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff_Kpimumu_fit1.eps");
  double nSel_pipi=myFitter1D.fit_resonant_Data("D2pipimumu",cut_pipimumu,q2Cut,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff_pipimumu_fit1.eps");

  double nTot_norm = myFitter1D.fit_normalization_Data(cut_PID_Kpimumu,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff_Kpimumu_fit2.eps");
  double nTot_pipi=myFitter1D.fit_resonant_Data("D2pipimumu",cut_PID_pipimumu,q2Cut,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff_pipimumu_fit2.eps");

  double Eff_norm = nSel_norm/nTot_norm;
  double Eff_pipi = nSel_pipi/nTot_pipi;

  double dEff_norm = 1/nTot_norm * TMath::Sqrt(nSel_norm*(1-(nSel_norm/nTot_norm) ) );
  double dEff_pipi = 1/nTot_pipi * TMath::Sqrt(nSel_pipi*(1-(nSel_pipi/nTot_pipi) ) );

  effBDTpipimumu_pipimumuTrained_data->SetBinContent(4,Eff_pipi);
  effBDTKpimumu_pipimumuTrained_data->SetBinContent(1,Eff_norm);
  effBDTpipimumu_pipimumuTrained_data->SetBinError(4,dEff_pipi);
  effBDTKpimumu_pipimumuTrained_data->SetBinError(1,dEff_norm);

  double relEff = Eff_pipi/Eff_norm;
  double dRelEff = TMath::Sqrt( (dEff_pipi/Eff_norm)*(dEff_pipi/Eff_norm) + (dEff_norm*Eff_pipi/(Eff_norm*Eff_norm))*(dEff_norm*Eff_pipi/(Eff_norm*Eff_norm)) );
  relEffBDTpipimumu_data->SetBinContent(4,relEff);
  relEffBDTpipimumu_data->SetBinError(4,dRelEff);

  std::cout<<"nSel_norm "<< nSel_norm <<std::endl;
  std::cout<<"nSel_pipi "<< nSel_pipi <<std::endl;
  std::cout<<"nTot_norm "<<  nTot_norm<<std::endl;
  std::cout<<"nTot_pipi "<<  nTot_pipi<<std::endl;
  std::cout<<"Eff_norm  "<<  Eff_norm<<"+-"<< dEff_norm <<std::endl;
  std::cout<<"Eff_pipi  "<<  Eff_pipi<<"+-"<<  dEff_pipi<<std::endl;
  std::cout<<"relEff  "<<  relEff<<"+-"<< dRelEff <<std::endl;
  std::cout<<"difference in ratio: "<<relEff-relEffBDTpipimumu->GetBinContent(4);


  std::cout<<"data fits done.."<<std::endl;

  

  ////////////////////
  //Draing part 

  dcastyle();

  TCanvas *a  = new TCanvas("pipimumu","Eff pipimumu");
  a->Divide(2,2);
  a->cd(1);
  effBDTpipimumu_pipimumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effBDTpipimumu_pipimumuTrained_magUp->Draw();
  effBDTpipimumu_pipimumuTrained_magDw->SetLineColor(kRed);
  effBDTpipimumu_pipimumuTrained_magDw->Draw("SAME");
  effBDTpipimumu_pipimumuTrained->SetLineColor(kBlack);
  effBDTpipimumu_pipimumuTrained->Draw("SAME");
  effBDTpipimumu_pipimumuTrained->Write();
  
  a->cd(2);
  relEffBDTpipimumu_magUp->GetYaxis()->SetRangeUser(0.5,1.3);
  relEffBDTpipimumu_magUp->Draw();
  relEffBDTpipimumu_magDw->SetLineColor(kRed);
  relEffBDTpipimumu_magDw->Draw("SAME");
  relEffBDTpipimumu->SetLineColor(kBlack);
  relEffBDTpipimumu->Draw("SAME");
  relEffBDTpipimumu->Write();

  a->cd(3);
  effBDTKKmumu_KKmumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effBDTKKmumu_KKmumuTrained_magUp->Draw();
  effBDTKKmumu_KKmumuTrained_magDw->SetLineColor(kRed);
  effBDTKKmumu_KKmumuTrained_magDw->Draw("SAME");
  effBDTKKmumu_KKmumuTrained->SetLineColor(kBlack);
  effBDTKKmumu_KKmumuTrained->Draw("SAME");
  effBDTKKmumu_KKmumuTrained->Write();

  a->cd(4);
  relEffBDTKKmumu_magUp->GetYaxis()->SetRangeUser(0.5,1.3);
  relEffBDTKKmumu_magUp->Draw();
  relEffBDTKKmumu_magDw->SetLineColor(kRed);
  relEffBDTKKmumu_magDw->Draw("SAME");
  relEffBDTKKmumu->SetLineColor(kBlack);
  relEffBDTKKmumu->Draw("SAME");
  relEffBDTKKmumu->Write();
  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff1.eps");
 

 a->Write();
 

  
  TCanvas *c  = new TCanvas("Kpimumu","Eff Kpimumu");

  c->Divide(2);
  c->cd(1);
  effBDTKpimumu_KKmumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effBDTKpimumu_KKmumuTrained_magUp->Draw();
  effBDTKpimumu_KKmumuTrained_magDw->SetLineColor(kRed);
  effBDTKpimumu_KKmumuTrained_magDw->Draw("SAME");
  effBDTKpimumu_KKmumuTrained->SetLineColor(kBlack);
  effBDTKpimumu_KKmumuTrained->Draw("SAME");
  c->cd(2);
  effBDTKpimumu_pipimumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effBDTKpimumu_pipimumuTrained_magUp->Draw();
  effBDTKpimumu_pipimumuTrained_magDw->SetLineColor(kRed);
  effBDTKpimumu_pipimumuTrained_magDw->Draw("SAME");
  effBDTKpimumu_pipimumuTrained->SetLineColor(kBlack);
  effBDTKpimumu_pipimumuTrained->Draw("SAME");


  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff2.eps");
  c->Write();

  TCanvas *d  = new TCanvas("Kpimumu","Eff Kpimumu");

  d->Divide(1,3);
  d->cd(1);
  effBDTKpimumu_pipimumuTrained_data->GetYaxis()->SetRangeUser(0.5,1);
  effBDTKpimumu_pipimumuTrained_data->Draw();
  effBDTKpimumu_pipimumuTrained->Draw("SAME");
  d->cd(2);
  effBDTpipimumu_pipimumuTrained_data->GetYaxis()->SetRangeUser(0.5,1);
  effBDTpipimumu_pipimumuTrained_data->Draw();
  effBDTpipimumu_pipimumuTrained->Draw("SAME");
  d->cd(3);
  relEffBDTpipimumu_data->GetYaxis()->SetRangeUser(0.9,1.1);
  relEffBDTpipimumu_data->Draw();
  relEffBDTpipimumu->Draw("SAME");
  d->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff3.eps");
  

  TCanvas *e  = new TCanvas("e","e");
  e->Divide(1,2);
  e->cd(1);
  relEffBDTKKmumu->GetYaxis()->SetRangeUser(0.9,1.1);
  relEffBDTKKmumu->Draw("");

  e->cd(2);
  relEffBDTpipimumu_data->GetYaxis()->SetRangeUser(0.9,1.1);
  relEffBDTpipimumu_data->Draw();
  relEffBDTpipimumu->Draw("SAME");
  e->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff4.eps");


  
  fOut->Write();
  fOut->Close();
  

}
			   

void MC_Combined_Hlt_BDT_efficiency(double ghostProbCut,double hadronPID,double muonPID,double BDT,bool applyPID,bool withGhosts){

  dcastyle;
  
  TString cutNshared="mu0_MuonNShared==0&&mu1_MuonNShared==0&&deltaM<146.6&&deltaM>144.5";
  TString cut_hadronPID_Kpimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_pipimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNpi>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_KKmumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNk>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_SlowpiPID = TString::Format("Slowpi_ProbNNghost<%.1f",ghostProbCut); 
  
  TString cutTriggerL0 = "(mu0_L0MuonDecision_TOS==1 || mu1_L0MuonDecision_TOS ==1)";

  TString cutTrigger_KKmumu = "(mu0_Hlt1TrackMuonDecision_TOS==1 || mu1_Hlt1TrackMuonDecision_TOS==1 || D_Hlt1TrackAllL0Decision_TOS==1)&&D_Hlt2CharmSemilepD02KKMuMuDecision_TOS==1";
  TString cutTrigger_Kpimumu = "(mu0_Hlt1TrackMuonDecision_TOS==1 || mu1_Hlt1TrackMuonDecision_TOS==1 || D_Hlt1TrackAllL0Decision_TOS==1)&&D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS==1";
  TString cutTrigger_pipimumu = "(mu0_Hlt1TrackMuonDecision_TOS==1 || mu1_Hlt1TrackMuonDecision_TOS==1 || D_Hlt1TrackAllL0Decision_TOS==1)&&D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS==1";


  TString cut_MuonPID=TString::Format("mu0_ProbNNghost<%.1f&&mu1_ProbNNghost<%.1f&&mu0_ProbNNmu>%.1f&&mu1_ProbNNmu>%.1f",ghostProbCut,ghostProbCut,muonPID,muonPID); 
 
  TString cut_BDT=TString::Format("BDT>%.1f",BDT);
 
  TString cut_norm_Kpimumu; 
  TString cut_norm_pipimumu;
  TString cut_norm_KKmumu ;
 
  TString cut_norm_Kpimumu_givenHlt; 
  TString cut_norm_pipimumu_givenHlt;
  TString cut_norm_KKmumu_givenHlt ;

  TString nameTarget ;
  nameTarget="/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/";
 
  //cut *PID* will the one be aplied to the normalization yield, cut_*mumu is the actual selection cut

  //if(applyTrigger && !applyPID) cout<<"trigger but no PID cut is not implemented right now. Either PID, PID and Trigger or neither of them."<<endl;

  
  TChain* tempChain_Kpi_KKmumuTrained = new TChain("BDT_Tree");
  TChain* tempChain_Kpi_pipimumuTrained = new TChain("BDT_Tree");
  TChain* tempChain_pipi_pipimumuTrained = new TChain("BDT_Tree");
  TChain* tempChain_KK_KKmumuTrained = new TChain("BDT_Tree");




  if(!applyPID) {
    //Form the string which are the "normalisation" cut later
    cut_norm_Kpimumu = cut_hadronPID_Kpimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID + "&&" + cutNshared +"&&"+ cutTriggerL0 ; 
    cut_norm_pipimumu = cut_hadronPID_pipimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID+ "&&" + cutNshared +"&&"+ cutTriggerL0; 
    cut_norm_KKmumu = cut_hadronPID_KKmumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID+ "&&" + cutNshared +"&&"+ cutTriggerL0;
        
    nameTarget+= "MCBDT_And_Hlt_Efficiencies_afterPID";
  }
 
  else {
    
    cut_norm_Kpimumu = cutNshared +"&&"+ cutTriggerL0; 
    cut_norm_pipimumu = cutNshared+"&&"+ cutTriggerL0; 
    cut_norm_KKmumu =  cutNshared +"&&"+ cutTriggerL0;
        
    nameTarget+= "MCBDT_and_Hlt_Efficiencies";
  }

  if(withGhosts) {
    cut_norm_Kpimumu+="&&(Dst_BKGCAT<11||Dst_BKGCAT==60)"; 
    cut_norm_pipimumu+="&&(Dst_BKGCAT<11||Dst_BKGCAT==60)";
    cut_norm_KKmumu+="&&(Dst_BKGCAT<11||Dst_BKGCAT==60)";
  
    tempChain_Kpi_KKmumuTrained->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2KKmumuBDT_multipleCand_withGhosts.root");
    tempChain_Kpi_pipimumuTrained->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2pipimumuBDT_multipleCand_withGhosts.root");//!!!!
    tempChain_pipi_pipimumuTrained->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_200.0_1600.0_D2pipimumuBDT_multipleCand_withGhosts.root");
    tempChain_KK_KKmumuTrained->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_200.0_1600.0_D2KKmumuBDT_multipleCand_withGhosts.root");
    nameTarget+= "withGhosts.root";

}
  else  {
    cut_norm_Kpimumu+="&&(Dst_BKGCAT<11)"; 
    cut_norm_pipimumu+="&&(Dst_BKGCAT<11)";
    cut_norm_KKmumu+="&&(Dst_BKGCAT<11)";
  
    tempChain_Kpi_KKmumuTrained->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2KKmumuBDT_multipleCand_noGhosts.root");
    tempChain_Kpi_pipimumuTrained->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2Kpimumu_200.0_1600.0_D2pipimumuBDT_multipleCand_noGhosts.root");//!!!!
    tempChain_pipi_pipimumuTrained->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2pipimumu_200.0_1600.0_D2pipimumuBDT_multipleCand_noGhosts.root");
    tempChain_KK_KKmumuTrained->AddFile("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudiesNoTruthmatching/MC_D2KKmumu_200.0_1600.0_D2KKmumuBDT_multipleCand_noGhosts.root");
    nameTarget+= "noGhosts.root";

  }

  //for BDT get conditional Eff after trigger Hlt cut
  cut_norm_Kpimumu_givenHlt=cut_norm_Kpimumu+"&&"+cutTrigger_Kpimumu;
  cut_norm_KKmumu_givenHlt= cut_norm_KKmumu+"&&"+cutTrigger_KKmumu;
  cut_norm_pipimumu_givenHlt=cut_norm_pipimumu+"&&"+cutTrigger_pipimumu;
  
  TString cut_Kpimumu_Hlt ;
  TString cut_pipimumu_Hlt;
  TString cut_KKmumu_Hlt ; 

  TString cut_Kpimumu_BDT ;
  TString cut_pipimumu_BDT;
  TString cut_KKmumu_BDT ; 

  //build the selection string for Hlt Eff
  cut_Kpimumu_Hlt=cut_norm_Kpimumu + "&&" + cutTrigger_Kpimumu;
  cut_KKmumu_Hlt=cut_norm_KKmumu + "&&"+ cutTrigger_KKmumu;
  cut_pipimumu_Hlt=cut_norm_pipimumu + "&&"+ cutTrigger_pipimumu;

  //build the selection string for BDT Eff
  cut_Kpimumu_BDT=cut_norm_Kpimumu_givenHlt + "&&" + cut_BDT  +"&&!isRejectedMultipleCandidate";
  cut_KKmumu_BDT=cut_norm_KKmumu_givenHlt + "&&" + cut_BDT    +"&&!isRejectedMultipleCandidate";
  cut_pipimumu_BDT=cut_norm_pipimumu_givenHlt + "&&" + cut_BDT+"&&!isRejectedMultipleCandidate";


  TFile* fOut= new TFile(nameTarget,"RECREATE"); 

  TH1D* effBDTKpimumu_KKmumuTrained;
  effBDTKpimumu_KKmumuTrained= new TH1D("effBDTKpimumu_KKmumuTrained","efficiency MC BDT D->Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D* effBDTKpimumu_pipimumuTrained; 
  effBDTKpimumu_pipimumuTrained= new TH1D("effBDTKpimumu_pipimumuTrained","efficiency MC BDT D->Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D* effBDTKKmumu_KKmumuTrained;
  effBDTKKmumu_KKmumuTrained= new TH1D("effBDTKKmumu_KKmumuTrained","efficiency MC BDT D->KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D* effBDTpipimumu_pipimumuTrained;
  effBDTpipimumu_pipimumuTrained= new TH1D("effBDTKpipimumu_pipimumuTrained","efficiency MC BDT D->pipimumu ",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1D* relEffBDTKKmumu;
  TH1D* relEffBDTpipimumu;
 
  relEffBDTKKmumu= new TH1D("relEffBDTDKKmumu","relative efficiency MC BDT D->KKmumu/D->Kpimumu ",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffBDTpipimumu= new TH1D("relEffBDTDpipimumu","relative efficiency MC BDT D->pipimumu/D->Kpimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);


  TH1D* effHltKpimumu_KKmumuTrained;
  effHltKpimumu_KKmumuTrained= new TH1D("effHltKpimumu_KKmumuTrained","efficiency MC Hlt D->Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D* effHltKpimumu_pipimumuTrained; 
  effHltKpimumu_pipimumuTrained= new TH1D("effHltKpimumu_pipimumuTrained","efficiency MC Hlt D->Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D* effHltKKmumu_KKmumuTrained;
  effHltKKmumu_KKmumuTrained= new TH1D("effHltKKmumu_KKmumuTrained","efficiency MC Hlt D->KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D* effHltpipimumu_pipimumuTrained;
  effHltpipimumu_pipimumuTrained= new TH1D("effHltKpipimumu_pipimumuTrained","efficiency MC Hlt D->pipimumu ",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1D* relEffHltKKmumu;
  TH1D* relEffHltpipimumu;
 
  relEffHltKKmumu= new TH1D("relEffHltDKKmumu","relative efficiency MC Hlt D->KKmumu/D->Kpimumu ",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffHltpipimumu= new TH1D("relEffHltDpipimumu","relative efficiency MC Hlt D->pipimumu/D->Kpimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1D* effCombinedHltBDTKpimumu_KKmumuTrained;
  effCombinedHltBDTKpimumu_KKmumuTrained= new TH1D("effCombinedHltBDTKpimumu_KKmumuTrained","efficiency MC CombinedHltBDT D->Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D* effCombinedHltBDTKpimumu_pipimumuTrained; 
  effCombinedHltBDTKpimumu_pipimumuTrained= new TH1D("effCombinedHltBDTKpimumu_pipimumuTrained","efficiency MC CombinedHltBDT D->Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D* effCombinedHltBDTKKmumu_KKmumuTrained;
  effCombinedHltBDTKKmumu_KKmumuTrained= new TH1D("effCombinedHltBDTKKmumu_KKmumuTrained","efficiency MC CombinedHltBDT D->KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D* effCombinedHltBDTpipimumu_pipimumuTrained;
  effCombinedHltBDTpipimumu_pipimumuTrained= new TH1D("effCombinedHltBDTKpipimumu_pipimumuTrained","efficiency MC CombinedHltBDT D->pipimumu ",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1D* relEffCombinedHltBDTKKmumu;
  TH1D* relEffCombinedHltBDTpipimumu;
 
  relEffCombinedHltBDTKKmumu= new TH1D("relEffCombinedHltBDTDKKmumu","relative efficiency MC CombinedHltBDT D->KKmumu/D->Kpimumu ",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffCombinedHltBDTpipimumu= new TH1D("relEffCombinedHltBDTDpipimumu","relative efficiency MC CombinedHltBDT D->pipimumu/D->Kpimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);



  
  for(int i=0; i<rangesKpi_low.size();++i){
   
    TString DimuonMassRange=TString::Format("&&D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKpi_low[i],rangesKpi_high[i]);
    TString temp_cut_Kpimumu_Hlt = cut_Kpimumu_Hlt+DimuonMassRange;
    TString temp_cut_Kpimumu_BDT =cut_Kpimumu_BDT+ DimuonMassRange;
    TString temp_cut_norm_Kpimumu_givenHlt=cut_norm_Kpimumu_givenHlt+ DimuonMassRange;
    TString temp_cut_norm_Kpimumu=cut_norm_Kpimumu+DimuonMassRange;


    //*********************Hlt*******************************************************

    //KKmumu trained
    double nNorm = (double)tempChain_Kpi_KKmumuTrained->GetEntries(cut_norm_Kpimumu);
    double  nSel =(double)tempChain_Kpi_KKmumuTrained->GetEntries(cut_Kpimumu_Hlt);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effHltKpimumu_KKmumuTrained->SetBinContent(i+1,tempEff);
    effHltKpimumu_KKmumuTrained->SetBinError(i+1,dTempEff);

    
    //pipimumu trained
    nNorm = (double)tempChain_Kpi_pipimumuTrained->GetEntries(cut_norm_Kpimumu);
    nSel = (double)tempChain_Kpi_pipimumuTrained->GetEntries(cut_Kpimumu_Hlt);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effHltKpimumu_pipimumuTrained->SetBinContent(i+1,tempEff);
    effHltKpimumu_pipimumuTrained->SetBinError(i+1,dTempEff);

    //********************* BDT *******************************************************

   
    //KKmumu trained
    nNorm = (double)tempChain_Kpi_KKmumuTrained->GetEntries(cut_norm_Kpimumu_givenHlt);
    nSel =(double)tempChain_Kpi_KKmumuTrained->GetEntries(cut_Kpimumu_BDT);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKpimumu_KKmumuTrained->SetBinContent(i+1,tempEff);
    effBDTKpimumu_KKmumuTrained->SetBinError(i+1,dTempEff);

  
    //pipimumu trained
    nNorm = (double)tempChain_Kpi_pipimumuTrained->GetEntries(cut_norm_Kpimumu_givenHlt);
    nSel = (double)tempChain_Kpi_pipimumuTrained->GetEntries(cut_Kpimumu_BDT);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKpimumu_pipimumuTrained->SetBinContent(i+1,tempEff);
    effBDTKpimumu_pipimumuTrained->SetBinError(i+1,dTempEff);

    //************************ combined ******************************************//

    
    //KKmumu trained
    nNorm = (double)tempChain_Kpi_KKmumuTrained->GetEntries(cut_norm_Kpimumu);
    nSel =(double)tempChain_Kpi_KKmumuTrained->GetEntries(cut_Kpimumu_BDT);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effCombinedHltBDTKpimumu_KKmumuTrained->SetBinContent(i+1,tempEff);
    effCombinedHltBDTKpimumu_KKmumuTrained->SetBinError(i+1,dTempEff);

    //pipimumu trained
    nNorm = (double)tempChain_Kpi_pipimumuTrained->GetEntries(cut_norm_Kpimumu);
    nSel = (double)tempChain_Kpi_pipimumuTrained->GetEntries(cut_Kpimumu_BDT);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effCombinedHltBDTKpimumu_pipimumuTrained->SetBinContent(i+1,tempEff);
    effCombinedHltBDTKpimumu_pipimumuTrained->SetBinError(i+1,dTempEff);
     
  }    


  
  for(int i=0; i<rangesKK_low.size();++i){
   
    TString DimuonMassRange=TString::Format("&&D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKK_low[i],rangesKK_high[i]);
    TString temp_cut_KKmumu_Hlt=cut_KKmumu_Hlt+DimuonMassRange;
    TString temp_cut_KKmumu_BDT=cut_KKmumu_BDT+DimuonMassRange;
    TString temp_cut_norm_KKmumu_givenHlt=cut_norm_KKmumu_givenHlt+DimuonMassRange;
    TString temp_cut_norm_KKmumu=cut_norm_KKmumu+DimuonMassRange;


    //*********************Hlt*******************************************************

    double nNorm = (double)tempChain_KK_KKmumuTrained->GetEntries(temp_cut_norm_KKmumu);
    double  nSel =(double)tempChain_KK_KKmumuTrained->GetEntries(temp_cut_KKmumu_Hlt);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effHltKKmumu_KKmumuTrained->SetBinContent(i+1,tempEff);
    effHltKKmumu_KKmumuTrained->SetBinError(i+1,dTempEff);

    std::cout<< temp_cut_norm_KKmumu  <<"  "<<temp_cut_KKmumu_Hlt <<"  "<<nNorm<<"  "<<nSel<<std::endl;


    double relEff=tempEff/effHltKpimumu_KKmumuTrained->GetBinContent(1);
    double dRelEff = TMath::Sqrt( TMath::Power(dTempEff/effHltKpimumu_KKmumuTrained->GetBinContent(1),2) + 
				  TMath::Power(effHltKpimumu_KKmumuTrained->GetBinError(1)*tempEff/(effHltKpimumu_KKmumuTrained->GetBinContent(1)*effHltKpimumu_KKmumuTrained->GetBinContent(1)) ,2) );
			  
    relEffHltKKmumu->SetBinContent(i+1,relEff);
    relEffHltKKmumu->SetBinError(i+1,dRelEff);

    //********************* BDT *******************************************************

   
    nNorm = (double)tempChain_KK_KKmumuTrained->GetEntries(temp_cut_norm_KKmumu_givenHlt);
    nSel =(double)tempChain_KK_KKmumuTrained->GetEntries(temp_cut_KKmumu_BDT);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTKKmumu_KKmumuTrained->SetBinContent(i+1,tempEff);
    effBDTKKmumu_KKmumuTrained->SetBinError(i+1,dTempEff);

    relEff=tempEff/effBDTKpimumu_KKmumuTrained->GetBinContent(1);
    dRelEff = TMath::Sqrt( TMath::Power(dTempEff/effBDTKpimumu_KKmumuTrained->GetBinContent(1),2) + 
				  TMath::Power(effBDTKpimumu_KKmumuTrained->GetBinError(1)*tempEff/(effBDTKpimumu_KKmumuTrained->GetBinContent(1)*effBDTKpimumu_KKmumuTrained->GetBinContent(1)) ,2) );

    relEffBDTKKmumu->SetBinContent(i+1,relEff);
    relEffBDTKKmumu->SetBinError(i+1,dRelEff);


    //************************ combined ******************************************//

    
    nNorm = (double)tempChain_KK_KKmumuTrained->GetEntries(temp_cut_norm_KKmumu);
    nSel =(double)tempChain_KK_KKmumuTrained->GetEntries(temp_cut_KKmumu_BDT);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effCombinedHltBDTKKmumu_KKmumuTrained->SetBinContent(i+1,tempEff);
    effCombinedHltBDTKKmumu_KKmumuTrained->SetBinError(i+1,dTempEff);
   
    relEff=tempEff/effCombinedHltBDTKpimumu_KKmumuTrained->GetBinContent(1);
    dRelEff = TMath::Sqrt( TMath::Power(dTempEff/effCombinedHltBDTKpimumu_KKmumuTrained->GetBinContent(1),2) + 
				  TMath::Power(effCombinedHltBDTKpimumu_KKmumuTrained->GetBinError(1)*tempEff/(effCombinedHltBDTKpimumu_KKmumuTrained->GetBinContent(1)*effCombinedHltBDTKpimumu_KKmumuTrained->GetBinContent(1)) ,2) );

    relEffCombinedHltBDTKKmumu->SetBinContent(i+1,relEff);
    relEffCombinedHltBDTKKmumu->SetBinError(i+1,dRelEff);
    
  
  }    



  for(int i=0; i<rangespipi_low.size();++i){

    TString DimuonMassRange=TString::Format("&&D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangespipi_low[i],rangespipi_high[i]);
    TString temp_cut_pipimumu_Hlt=cut_pipimumu_Hlt+DimuonMassRange;
    TString temp_cut_pipimumu_BDT=cut_pipimumu_BDT+DimuonMassRange;
    TString temp_cut_norm_pipimumu_givenHlt=cut_norm_pipimumu_givenHlt+DimuonMassRange;
    TString temp_cut_norm_pipimumu=cut_norm_pipimumu+DimuonMassRange;

    //*********************Hlt*******************************************************

    double nNorm = (double)tempChain_pipi_pipimumuTrained->GetEntries(temp_cut_norm_pipimumu);
    double  nSel =(double)tempChain_pipi_pipimumuTrained->GetEntries(temp_cut_pipimumu_Hlt);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effHltpipimumu_pipimumuTrained->SetBinContent(i+1,tempEff);
    effHltpipimumu_pipimumuTrained->SetBinError(i+1,dTempEff);

    double relEff=tempEff/effHltKpimumu_pipimumuTrained->GetBinContent(1);
    double dRelEff = TMath::Sqrt( TMath::Power(dTempEff/effHltKpimumu_pipimumuTrained->GetBinContent(1),2) + 
				  TMath::Power(effHltKpimumu_pipimumuTrained->GetBinError(1)*tempEff/(effHltKpimumu_pipimumuTrained->GetBinContent(1)*effHltKpimumu_pipimumuTrained->GetBinContent(1)) ,2) );
			  
    relEffHltpipimumu->SetBinContent(i+1,relEff);
    relEffHltpipimumu->SetBinError(i+1,dRelEff);

    //********************* BDT *******************************************************

   
    nNorm = (double)tempChain_pipi_pipimumuTrained->GetEntries(temp_cut_norm_pipimumu_givenHlt);
    nSel =(double)tempChain_pipi_pipimumuTrained->GetEntries(temp_cut_pipimumu_BDT);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effBDTpipimumu_pipimumuTrained->SetBinContent(i+1,tempEff);
    effBDTpipimumu_pipimumuTrained->SetBinError(i+1,dTempEff);

    relEff=tempEff/effBDTKpimumu_pipimumuTrained->GetBinContent(1);
    dRelEff = TMath::Sqrt( TMath::Power(dTempEff/effBDTKpimumu_pipimumuTrained->GetBinContent(1),2) + 
				  TMath::Power(effBDTKpimumu_pipimumuTrained->GetBinError(1)*tempEff/(effBDTKpimumu_pipimumuTrained->GetBinContent(1)*effBDTKpimumu_pipimumuTrained->GetBinContent(1)) ,2) );

    relEffBDTpipimumu->SetBinContent(i+1,relEff);
    relEffBDTpipimumu->SetBinError(i+1,dRelEff);


    //************************ combined ******************************************//

    
    nNorm = (double)tempChain_pipi_pipimumuTrained->GetEntries(temp_cut_norm_pipimumu);
    nSel =(double)tempChain_pipi_pipimumuTrained->GetEntries(temp_cut_pipimumu_BDT);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effCombinedHltBDTpipimumu_pipimumuTrained->SetBinContent(i+1,tempEff);
    effCombinedHltBDTpipimumu_pipimumuTrained->SetBinError(i+1,dTempEff);
    
    relEff=tempEff/effCombinedHltBDTKpimumu_pipimumuTrained->GetBinContent(1);
    dRelEff = TMath::Sqrt( TMath::Power(dTempEff/effCombinedHltBDTKpimumu_pipimumuTrained->GetBinContent(1),2) + 
				  TMath::Power(effCombinedHltBDTKpimumu_pipimumuTrained->GetBinError(1)*tempEff/(effCombinedHltBDTKpimumu_pipimumuTrained->GetBinContent(1)*effCombinedHltBDTKpimumu_pipimumuTrained->GetBinContent(1)) ,2) );

    relEffCombinedHltBDTpipimumu->SetBinContent(i+1,relEff);
    relEffCombinedHltBDTpipimumu->SetBinError(i+1,dRelEff);
    

 
  }    
  


  ////////////////////
  //Draing part 

  dcastyle();

  TCanvas *a  = new TCanvas("pipimumu","Eff pipimumu");
  a->Divide(2,2);
  a->cd(1);
  effBDTpipimumu_pipimumuTrained->SetLineColor(kBlack);
  effBDTpipimumu_pipimumuTrained->Draw("SAME");
  effBDTpipimumu_pipimumuTrained->Write();
  
  a->cd(2);
  relEffBDTpipimumu->SetLineColor(kBlack);
  relEffBDTpipimumu->Draw("SAME");
  relEffBDTpipimumu->Write();

  a->cd(3);
  effBDTKKmumu_KKmumuTrained->SetLineColor(kBlack);
  effBDTKKmumu_KKmumuTrained->Draw("SAME");
  effBDTKKmumu_KKmumuTrained->Write();

  a->cd(4);
  relEffBDTKKmumu->SetLineColor(kBlack);
  relEffBDTKKmumu->Draw("SAME");
  relEffBDTKKmumu->Write();
  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff1.eps");
 

 a->Write();
 

  
  TCanvas *c  = new TCanvas("Kpimumu","Eff Kpimumu");

  c->Divide(2);
  c->cd(1);
  effBDTKpimumu_KKmumuTrained->SetLineColor(kBlack);
  effBDTKpimumu_KKmumuTrained->Draw("SAME");
  c->cd(2);
  effBDTKpimumu_pipimumuTrained->SetLineColor(kBlack);
  effBDTKpimumu_pipimumuTrained->Draw("SAME");


  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff2.eps");
  c->Write();

  TCanvas *d  = new TCanvas("Kpimumu","Eff Kpimumu");

  d->Divide(1,3);
  d->cd(1);
  effBDTKpimumu_pipimumuTrained->Draw("SAME");
  d->cd(2);
  effBDTpipimumu_pipimumuTrained->Draw("SAME");
  d->cd(3);
  relEffBDTpipimumu->Draw("SAME");
  d->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff3.eps");
  

  TCanvas *e  = new TCanvas("e","e");
  e->Divide(1,2);
  e->cd(1);
  relEffBDTKKmumu->GetYaxis()->SetRangeUser(0.9,1.1);
  relEffBDTKKmumu->Draw("");

  e->cd(2);
  relEffBDTpipimumu->Draw("SAME");
  e->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MC_BDT_Eff4.eps");
  
  fOut->Write();
  fOut->Close();
  

}
			   
void MC_L0Trigger_efficiency(double ghostProbCut,double hadronPID,double muonPID,double BDT,bool applyPID, bool applyBDT, TString cut_Trigger_TIS){


  TString cutNshared="mu0_MuonNShared==0&&mu1_MuonNShared==0&&deltaM>144.5&&deltaM<146.5&&(mHH>850&&mHH<950)";

  TString cut_hadronPID_Kpimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_pipimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNpi>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_KKmumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNk>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_SlowpiPID = TString::Format("Slowpi_ProbNNghost<%.1f",ghostProbCut); 

  TString cut_MuonPID=TString::Format("mu0_ProbNNghost<%.1f&&mu1_ProbNNghost<%.1f&&mu0_ProbNNmu>%.1f&&mu1_ProbNNmu>%.1f",ghostProbCut,ghostProbCut,muonPID,muonPID); 
  
  TString cut_BDT=TString::Format("BDT>%.1f",BDT);

  //TString cut_Trigger="(D_L0MuonDecision_TOS)";
  TString cut_Trigger="(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)";
  //TString cut_Trigger="mu0_L0MuonDecision_TOS";
  //TString cut_Trigger="( ( (mu0_P>mu1_P)&&mu0_L0MuonDecision_TOS) || ( (mu1_P>mu0_P)&&mu1_L0MuonDecision_TOS) )";
  
  //  TString cut_Trigger_TIS="D_L0MuonDecision_TIS";
  
  

  TString cut_PID_Kpimumu; 
  TString cut_PID_pipimumu;
  TString cut_PID_KKmumu ;

  TString cut_Kpimumu ;
  TString cut_pipimumu;
  TString cut_KKmumu ; 

  TString cut_norm_Kpimumu ;
  TString cut_norm_pipimumu;
  TString cut_norm_KKmumu ; 
  
  TString nameTarget = "MCL0Trigger_Efficiencies";
 
    cut_PID_Kpimumu = cut_hadronPID_Kpimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID + "&&" + cutNshared;
    cut_PID_pipimumu = cut_hadronPID_pipimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID + "&&" + cutNshared;
    cut_PID_KKmumu = cut_hadronPID_KKmumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID + "&&" + cutNshared;

  if(applyPID) {
 
    cut_Kpimumu=cut_PID_Kpimumu;
    cut_pipimumu=cut_PID_pipimumu;  
    cut_KKmumu=cut_PID_KKmumu;
    nameTarget = nameTarget+"_afterPID";
  }
  if(applyBDT) {

    if(applyPID){
      cut_Kpimumu = cut_Kpimumu+"&&"+cut_BDT;
      cut_pipimumu = cut_pipimumu+"&&"+cut_BDT;
      cut_KKmumu =   cut_pipimumu+"&&"+cut_BDT;
    }
    else{
    cut_Kpimumu = cut_BDT;
    cut_pipimumu =cut_BDT;
    cut_KKmumu =  cut_BDT;;
    }
    nameTarget = nameTarget+"_afterBDT";
  }


  //again norm is the 'loser cut' used for the yield used in normalisation
  cut_norm_Kpimumu = cut_Kpimumu;
  cut_norm_pipimumu =cut_pipimumu;
  cut_norm_KKmumu = cut_KKmumu;

  cut_Kpimumu = cut_Kpimumu + "&&" + cut_Trigger;
  cut_pipimumu =cut_pipimumu + "&&" + cut_Trigger;
  cut_KKmumu = cut_KKmumu + "&&" + cut_Trigger;

  if(!applyBDT && !applyPID) {
    
    cut_Kpimumu= cutNshared;
    cut_pipimumu= cutNshared;
    cut_KKmumu= cutNshared;
      
    cut_norm_Kpimumu = cut_Kpimumu;
    cut_norm_pipimumu =cut_pipimumu;
    cut_norm_KKmumu = cut_KKmumu;
    cut_Kpimumu = cut_Kpimumu + "&&" + cut_Trigger;
    cut_pipimumu =cut_pipimumu + "&&" + cut_Trigger;
    cut_KKmumu = cut_KKmumu + "&&" + cut_Trigger;

  }

  nameTarget+=".root";
  
  cout<<"Cuts:S "<<endl;
  cout<<cut_Kpimumu+"&&"+cut_Trigger_TIS<<endl;
  cout<<cut_norm_Kpimumu+"&&"+cut_Trigger_TIS<<endl;


  TFile* fOut= new TFile(nameTarget,"RECREATE"); 

  TH1D* effTriggerKpimumu_KKmumuTrained_magUp;
  TH1D* effTriggerKpimumu_KKmumuTrained_magDw;

  effTriggerKpimumu_KKmumuTrained_magUp= new TH1D("effTriggerKpimumu_KKmumuTrained_magUp","efficiency MC Trigger D->Kpimumu magUp",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effTriggerKpimumu_KKmumuTrained_magDw= new TH1D("effTriggerKpimumu_KKmumuTrained_magDw","efficiency MC Trigger D->Kpimumu magDw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D* effTriggerKpimumu_pipimumuTrained_magUp;
  TH1D* effTriggerKpimumu_pipimumuTrained_TISTOS_magUp;
  TH1D* effTriggerKpimumu_pipimumuTrained_magDw;
  TH1D* effTriggerKpimumu_pipimumuTrained_data;

  effTriggerKpimumu_pipimumuTrained_TISTOS_magUp= new TH1D("effTriggerKpimumu_pipimumuTrained_TISTOS_magUp","efficiency MC Trigger D->Kpimumu magUp",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effTriggerKpimumu_pipimumuTrained_magUp= new TH1D("effTriggerKpimumu_pipimumuTrained_magUp","efficiency MC Trigger D->Kpimumu magUp",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effTriggerKpimumu_pipimumuTrained_magDw= new TH1D("effTriggerKpimumu_pipimumuTrained_magDw","efficiency MC Trigger D->Kpimumu magDw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effTriggerKpimumu_pipimumuTrained_data= new TH1D("effTriggerKpimumu_pipimumuTrained_data","efficiency MC Trigger D->Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D* effTriggerKKmumu_KKmumuTrained_magUp;
  TH1D* effTriggerKKmumu_KKmumuTrained_magDw;

  effTriggerKKmumu_KKmumuTrained_magUp= new TH1D("effTriggerKKmumu_KKmumuTrained_magUp","efficiency MC Trigger D->KKmumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  effTriggerKKmumu_KKmumuTrained_magDw= new TH1D("effTriggerKKmumu_KKmumuTrained_magDw","efficiency MC Trigger D->KKmumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1D* effTriggerpipimumu_pipimumuTrained_magUp;
  TH1D* effTriggerpipimumu_pipimumuTrained_TISTOS_magUp;
  TH1D* effTriggerpipimumu_pipimumuTrained_magDw;
  TH1D* effTriggerpipimumu_pipimumuTrained_data;

  effTriggerpipimumu_pipimumuTrained_magUp= new TH1D("effTriggerpipimumu_pipimumuTrained_magUp","efficiency MC Trigger D->pipimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effTriggerpipimumu_pipimumuTrained_TISTOS_magUp= new TH1D("effTriggerpipimumu_pipimumuTrained_TISTOS_magUp","efficiency MC Trigger D->pipimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effTriggerpipimumu_pipimumuTrained_magDw= new TH1D("effTriggerKpipimumu_pipimumuTrained_magDw","efficiency MC Trigger D->pipimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effTriggerpipimumu_pipimumuTrained_data= new TH1D("effTriggerpipipimumu_pipimumuTrained_data","efficiency MC Trigger D->pipimumu ",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1D* relEffTriggerKKmumu_magUp;
  TH1D* relEffTriggerKKmumu_magDw;
  TH1D* relEffTriggerpipimumu_magUp;
  TH1D* relEffTriggerpipimumu_magDw;
  TH1D* relEffTriggerpipimumu_data;
  TH1D* relEffTriggerpipimumu_TISTOS_magUp;

  relEffTriggerpipimumu_data= new TH1D("relEffTriggerDpipimumu_data","relative efficiency MC Trigger D->pipimumu/D->Kpimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffTriggerKKmumu_magUp= new TH1D("relEffTriggerDKKmumu_magUp","relative efficiency MC Trigger D->KKmumu/D->Kpimumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffTriggerKKmumu_magDw= new TH1D("relEffTriggerDKKmumu_magDw","relative efficiency MC Trigger D->KKmumu/D->Kpimumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffTriggerpipimumu_magUp= new TH1D("relEffTriggerDpipimumu_magUp","relative efficiency MC Trigger D->pipimumu/D->Kpimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffTriggerpipimumu_TISTOS_magUp= new TH1D("relEffTriggerpipimumu_TISTOS_magUp","relative efficiency MC Trigger D->pipimumu/D->Kpimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffTriggerpipimumu_magDw= new TH1D("relEffTriggerDpipimumu_magDw","relative efficiency MC Trigger D->pipimumu/D->Kpimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  
  ofstream myfile;
  TString textFile = "../img/EfficiencyStudies/triggerTests/"+cut_Trigger_TIS;
  if(applyBDT) textFile +="_afterBDT";
  if(applyPID) textFile +="_afterPID";

  myfile.open(textFile+".txt");

  for(int i=0; i<rangesKpi_low.size();++i){
    
    TChain* tempChain_KKmumuTrained_magUp= new TChain("BDT_Tree");
    tempChain_KKmumuTrained_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
    TChain* tempChain_KKmumuTrained_magDw= new TChain("BDT_Tree");
    tempChain_KKmumuTrained_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));

    TChain* tempChain_pipimumuTrained_magUp= new TChain("BDT_Tree");
    tempChain_pipimumuTrained_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
    TChain* tempChain_pipimumuTrained_magDw= new TChain("BDT_Tree");
    tempChain_pipimumuTrained_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));


    //mag up KKmumu trained
    double nNorm = (double)tempChain_KKmumuTrained_magUp->GetEntries(cut_norm_Kpimumu);
    double  nSel = (double)tempChain_KKmumuTrained_magUp->GetEntries(cut_Kpimumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTriggerKpimumu_KKmumuTrained_magUp->SetBinContent(i+1,tempEff);
    effTriggerKpimumu_KKmumuTrained_magUp->SetBinError(i+1,dTempEff);

    //mag Dw KKmumu trained
    nNorm = tempChain_KKmumuTrained_magDw->GetEntries(cut_norm_Kpimumu);
    nSel = tempChain_KKmumuTrained_magDw->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTriggerKpimumu_KKmumuTrained_magDw->SetBinContent(i+1,tempEff);
    effTriggerKpimumu_KKmumuTrained_magDw->SetBinError(i+1,dTempEff);

    //mag up KKmumu trained
    nNorm = tempChain_pipimumuTrained_magUp->GetEntries(cut_norm_Kpimumu);
    nSel = tempChain_pipimumuTrained_magUp->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTriggerKpimumu_pipimumuTrained_magUp->SetBinContent(i+1,tempEff);
    effTriggerKpimumu_pipimumuTrained_magUp->SetBinError(i+1,dTempEff);

    //mag dw pipimumu trained
    nNorm = tempChain_pipimumuTrained_magDw->GetEntries(cut_norm_Kpimumu);
    nSel = tempChain_pipimumuTrained_magDw->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTriggerKpimumu_pipimumuTrained_magDw->SetBinContent(i+1,tempEff);
    effTriggerKpimumu_pipimumuTrained_magDw->SetBinError(i+1,dTempEff);

    //mag up
    nNorm = tempChain_pipimumuTrained_magUp->GetEntries(cut_norm_Kpimumu);
    nSel = tempChain_pipimumuTrained_magUp->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    myfile<<"MC Kpimumu (up): nNorm= "<<nNorm<<" nSel= "<< nSel <<" Eff "<<tempEff <<"+-"<<dTempEff<<std::endl; 
    effTriggerKpimumu_pipimumuTrained_magUp->SetBinContent(i+1,tempEff);
    effTriggerKpimumu_pipimumuTrained_magUp->SetBinError(i+1,dTempEff);
    
    //TISTOS magUp
    nNorm = tempChain_pipimumuTrained_magUp->GetEntries(cut_norm_Kpimumu+"&&"+cut_Trigger_TIS);
    nSel = tempChain_pipimumuTrained_magUp->GetEntries(cut_Kpimumu+"&&"+cut_Trigger_TIS);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    myfile<<"MC Kpimumu TISOS: nNorm= "<<nNorm<<" nSel= "<< nSel <<" Eff "<<tempEff <<"+-"<<dTempEff<<std::endl; 
    effTriggerKpimumu_pipimumuTrained_TISTOS_magUp->SetBinContent(i+1,tempEff);
    effTriggerKpimumu_pipimumuTrained_TISTOS_magUp->SetBinError(i+1,dTempEff);

  }    

  for(int i=0; i<rangesKK_low.size();++i){
    
    TChain* tempChain_magUp= new TChain("BDT_Tree");
    tempChain_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]));
    TChain* tempChain_magDw= new TChain("BDT_Tree");
    tempChain_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]));

    //mag up
    double nNorm = (double)tempChain_magUp->GetEntries(cut_norm_KKmumu);
    double nSel = (double)tempChain_magUp->GetEntries(cut_KKmumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTriggerKKmumu_KKmumuTrained_magUp->SetBinContent(i+1,tempEff);
    effTriggerKKmumu_KKmumuTrained_magUp->SetBinError(i+1,dTempEff);
    double relEff = tempEff/effTriggerKpimumu_KKmumuTrained_magUp->GetBinContent(1);
    double dRelEff = TMath::Sqrt( (dTempEff/effTriggerKKmumu_KKmumuTrained_magUp->GetBinContent(1))*(dTempEff/effTriggerKKmumu_KKmumuTrained_magUp->GetBinContent(1)) 
				  + (tempEff*effTriggerKKmumu_KKmumuTrained_magUp->GetBinError(1)/(effTriggerKKmumu_KKmumuTrained_magUp->GetBinContent(1)*effTriggerKKmumu_KKmumuTrained_magUp->GetBinContent(1)) * 
				     (tempEff*effTriggerKKmumu_KKmumuTrained_magUp->GetBinError(1)/(effTriggerKKmumu_KKmumuTrained_magUp->GetBinContent(1)*effTriggerKKmumu_KKmumuTrained_magUp->GetBinContent(1))
				      ) ));
    relEffTriggerKKmumu_magUp->SetBinContent(i+1,relEff);				     
    relEffTriggerKKmumu_magUp->SetBinError(i+1,dRelEff);			     
    
    //mag Dw
    nNorm = tempChain_magDw->GetEntries(cut_norm_KKmumu);
    nSel = tempChain_magDw->GetEntries(cut_KKmumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTriggerKKmumu_KKmumuTrained_magDw->SetBinContent(i+1,tempEff);
    effTriggerKKmumu_KKmumuTrained_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effTriggerKpimumu_KKmumuTrained_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effTriggerKKmumu_KKmumuTrained_magDw->GetBinContent(1))*(dTempEff/effTriggerKKmumu_KKmumuTrained_magDw->GetBinContent(1)) 
				  + (tempEff*effTriggerKKmumu_KKmumuTrained_magDw->GetBinError(1)/(effTriggerKKmumu_KKmumuTrained_magDw->GetBinContent(1)*effTriggerKKmumu_KKmumuTrained_magDw->GetBinContent(1)) * 
				     (tempEff*effTriggerKKmumu_KKmumuTrained_magDw->GetBinError(1)/(effTriggerKKmumu_KKmumuTrained_magDw->GetBinContent(1)*effTriggerKKmumu_KKmumuTrained_magDw->GetBinContent(1))
				      ) ));
    relEffTriggerKKmumu_magDw->SetBinContent(i+1,relEff);				     
    relEffTriggerKKmumu_magDw->SetBinError(i+1,dRelEff);			     
      

  }    

  for(int i=0; i<rangespipi_low.size();++i){
    
    TChain* tempChain_magUp= new TChain("BDT_Tree");
    tempChain_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]));
    TChain* tempChain_magDw= new TChain("BDT_Tree");
    tempChain_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]));

    //mag up
    double nNorm = (double)tempChain_magUp->GetEntries(cut_norm_pipimumu);
    double nSel = (double)tempChain_magUp->GetEntries(cut_pipimumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    myfile<<"MC pipimumu (up): nNorm= "<<nNorm<<" nSel= "<< nSel <<" Eff "<<tempEff <<"+-"<<dTempEff<<std::endl; 
    effTriggerpipimumu_pipimumuTrained_magUp->SetBinContent(i+1,tempEff);
    effTriggerpipimumu_pipimumuTrained_magUp->SetBinError(i+1,dTempEff);
    double relEff = tempEff/effTriggerKpimumu_pipimumuTrained_magUp->GetBinContent(1);
    double dRelEff = TMath::Sqrt( (dTempEff/effTriggerpipimumu_pipimumuTrained_magUp->GetBinContent(1))*(dTempEff/effTriggerpipimumu_pipimumuTrained_magUp->GetBinContent(1)) 
				  + (tempEff*effTriggerpipimumu_pipimumuTrained_magUp->GetBinError(1)/(effTriggerpipimumu_pipimumuTrained_magUp->GetBinContent(1)*effTriggerpipimumu_pipimumuTrained_magUp->GetBinContent(1)) * 
				     (tempEff*effTriggerpipimumu_pipimumuTrained_magUp->GetBinError(1)/(effTriggerpipimumu_pipimumuTrained_magUp->GetBinContent(1)*effTriggerpipimumu_pipimumuTrained_magUp->GetBinContent(1))
				      ) ));
    relEffTriggerpipimumu_magUp->SetBinContent(i+1,relEff);				     
    relEffTriggerpipimumu_magUp->SetBinError(i+1,dRelEff);			     

    //TISTOS
    nNorm = (double)tempChain_magUp->GetEntries(cut_norm_pipimumu+"&&"+cut_Trigger_TIS);
    nSel = (double)tempChain_magUp->GetEntries(cut_pipimumu+"&&"+cut_Trigger_TIS);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->SetBinContent(i+1,tempEff);
    effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->SetBinError(i+1,dTempEff);
    myfile<<"MC pipimumu (TISTOS): nNorm= "<<nNorm<<" nSel= "<< nSel <<" Eff "<<tempEff <<"+-"<<dTempEff<<std::endl; 
    relEff = tempEff/effTriggerKpimumu_pipimumuTrained_TISTOS_magUp->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->GetBinContent(1))*(dTempEff/effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->GetBinContent(1)) 
				  + (tempEff*effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->GetBinError(1)/(effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->GetBinContent(1)*effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->GetBinContent(1)) * 
				     (tempEff*effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->GetBinError(1)/(effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->GetBinContent(1)*effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->GetBinContent(1))
				      ) ));
    relEffTriggerpipimumu_TISTOS_magUp->SetBinContent(i+1,relEff);				     
    relEffTriggerpipimumu_TISTOS_magUp->SetBinError(i+1,dRelEff);			     
    

    //mag Dw
    nNorm = (double)tempChain_magDw->GetEntries(cut_norm_pipimumu);
    nSel = (double)tempChain_magDw->GetEntries(cut_pipimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTriggerpipimumu_pipimumuTrained_magDw->SetBinContent(i+1,tempEff);
    effTriggerpipimumu_pipimumuTrained_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effTriggerKpimumu_pipimumuTrained_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effTriggerpipimumu_pipimumuTrained_magDw->GetBinContent(1))*(dTempEff/effTriggerpipimumu_pipimumuTrained_magDw->GetBinContent(1)) 
				  + (tempEff*effTriggerpipimumu_pipimumuTrained_magDw->GetBinError(1)/(effTriggerpipimumu_pipimumuTrained_magDw->GetBinContent(1)*effTriggerpipimumu_pipimumuTrained_magDw->GetBinContent(1)) * 
				     (tempEff*effTriggerpipimumu_pipimumuTrained_magDw->GetBinError(1)/(effTriggerpipimumu_pipimumuTrained_magDw->GetBinContent(1)*effTriggerpipimumu_pipimumuTrained_magDw->GetBinContent(1))
				      ) ));
    relEffTriggerpipimumu_magDw->SetBinContent(i+1,relEff);				     
    relEffTriggerpipimumu_magDw->SetBinError(i+1,dRelEff);			     
      

  }    


  ////////////////////////////////////////////////
  ////////// data part

  //ghost prob and nshared cut already aplied to the *noCuts* files
  TString kind = "D2pipimumu";
  TString dataCut = "mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&BDT>0";
  TString dataCutNoBDT="mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5";
  TString misIDCut = "mu0_ProbNNmu>0.5&&BDT>0";
  TString q2Cut= "D_DiMuon_Mass>565&&D_DiMuon_Mass<950";
  //TString q2Cut= "D_DiMuon_Mass>950&&D_DiMuon_Mass<1100";

  TString q2RangeNormalizationMode="D_DiMuon_Mass>675&&D_DiMuon_Mass<875";

  
  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MC_D2Kpimumu_D2KKmumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2Kpipipi_PIDline_D2pipimumuBDT.root");
  myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/D2pipipipi_PIDline_D2pipimumuBDT.root");
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/noPreselectionCuts/D2Kpimumu_D2pipimumuBDT_noCuts.root");
  myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/noPreselectionCuts/D2pipimumu_D2pipimumuBDT_noCuts.root");

  
  //didnt create extra files for misID and MC without preselection
  myFitter1D.fit_MC(dataCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_"+dataCut+q2Cut+".eps");
  myFitter1D.fit_normalization_MC(dataCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/MC_normMode_"+dataCut+".eps");
  std::cout<<"Monte Carlo fits done.."<<std::endl;

  myFitter1D.fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_Kpi"+misIDCut+q2RangeNormalizationMode+".eps");
  myFitter1D.fit_HHpipi_misID(misIDCut+"&&"+q2Cut,true,"/work/mitzel/D2hhmumu/dev/"+kind+"/img/Fits/full1DFits/misID_pipi"+misIDCut+q2Cut+".eps"); ////q2 cut missing? check!        
  std::cout<<"misID fits done.."<<std::endl;


  //use TIS TOS method
  double nSel_norm = myFitter1D.fit_normalization_Data(cut_Kpimumu+"&&"+q2RangeNormalizationMode+"&&"+cut_Trigger_TIS,textFile+"1.eps"); //TIS&&TOS
  double nSel_pipi=myFitter1D.fit_resonant_Data("D2pipimumu",cut_pipimumu,q2Cut+"&&"+cut_Trigger_TIS,textFile+"2.eps");//TIS&&TOS

  double nTot_norm = myFitter1D.fit_normalization_Data(cut_norm_Kpimumu+"&&"+q2RangeNormalizationMode+"&&"+cut_Trigger_TIS,textFile+"3.eps");  //TIS
  double nTot_pipi=myFitter1D.fit_resonant_Data("D2pipimumu",cut_norm_pipimumu,q2Cut+"&&"+cut_Trigger_TIS,textFile+"4.eps"); //TIS
  
  double Eff_norm = nSel_norm/nTot_norm;
  double Eff_pipi = nSel_pipi/nTot_pipi;

  double dEff_norm = 1/nTot_norm * TMath::Sqrt(nSel_norm*(1-(nSel_norm/nTot_norm) ) );
  double dEff_pipi = 1/nTot_pipi * TMath::Sqrt(nSel_pipi*(1-(nSel_pipi/nTot_pipi) ) );

  effTriggerpipimumu_pipimumuTrained_data->SetBinContent(4,Eff_pipi);///
  effTriggerKpimumu_pipimumuTrained_data->SetBinContent(1,Eff_norm);///
  effTriggerpipimumu_pipimumuTrained_data->SetBinError(4,dEff_pipi);///
  effTriggerKpimumu_pipimumuTrained_data->SetBinError(1,dEff_norm);///

  double relEff = Eff_pipi/Eff_norm;
  double dRelEff = TMath::Sqrt( (dEff_pipi/Eff_norm)*(dEff_pipi/Eff_norm) + (dEff_norm*Eff_pipi/(Eff_norm*Eff_norm))*(dEff_norm*Eff_pipi/(Eff_norm*Eff_norm)) );
  relEffTriggerpipimumu_data->SetBinContent(4,relEff);////
  relEffTriggerpipimumu_data->SetBinError(4,dRelEff);////


  myfile<<"Cuts:S "<<endl;
  myfile<<cut_Kpimumu+"&&"+cut_Trigger_TIS<<endl;
  myfile<<cut_norm_Kpimumu+"&&"+cut_Trigger_TIS<<endl;

  myfile<<"nSel_norm "<< nSel_norm <<endl;
  myfile<<"nSel_pipi "<< nSel_pipi <<endl;
  myfile<<"nTot_norm "<<  nTot_norm<<endl;
  myfile<<"nTot_pipi "<<  nTot_pipi<<endl;
  myfile<<"data Eff_norm  "<<  Eff_norm<<"+-"<< dEff_norm <<endl;
  myfile<<"data Eff_pipi  "<<  Eff_pipi<<"+-"<<  dEff_pipi<<endl;
  myfile<<"data relEff  "<<  relEff<<"+-"<< dRelEff <<endl;
  myfile<<"MC (norm): "<< effTriggerKpimumu_pipimumuTrained_magUp->GetBinContent(1) <<endl;
  myfile<<"and Sig " << effTriggerpipimumu_pipimumuTrained_magUp->GetBinContent(4) << endl;///
  myfile<<"MC ratio: " <<relEffTriggerpipimumu_magUp->GetBinContent(4)<<endl; ///
  myfile<<"difference in ratio: "<<relEff-relEffTriggerpipimumu_magUp->GetBinContent(4)<<endl;///
  myfile<<"TISTOS MC (norm): "<< effTriggerKpimumu_pipimumuTrained_TISTOS_magUp->GetBinContent(1) <<endl;
  myfile<<"TISTOS MC (sig) " << effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->GetBinContent(4) <<endl;///
  myfile<<"TISTOS MC ratio:" <<relEffTriggerpipimumu_TISTOS_magUp->GetBinContent(4)<<endl;///
  
  myfile.close();

  std::cout<<"data fits done.."<<std::endl;




  //Draing part 

  TCanvas *a  = new TCanvas("pipimumu","Eff pipimumu");
  a->Divide(3,2);
  a->cd(1);
  effTriggerpipimumu_pipimumuTrained_magUp->GetYaxis()->SetRangeUser(0.3,1);
  effTriggerpipimumu_pipimumuTrained_magUp->Draw();
  //effTriggerpipimumu_pipimumuTrained_magDw->SetLineColor(kRed);
  //effTriggerpipimumu_pipimumuTrained_magDw->Draw("SAME");
  effTriggerpipimumu_pipimumuTrained_data->SetLineColor(kOrange);
  effTriggerpipimumu_pipimumuTrained_data->Draw("SAME");

  effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->SetLineColor(kViolet);
  effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->Draw("SAME");

  a->cd(2);
  relEffTriggerpipimumu_magUp->GetYaxis()->SetRangeUser(0.8,1.5);
  relEffTriggerpipimumu_magUp->Draw();
  //relEffTriggerpipimumu_magDw->SetLineColor(kRed);
  //relEffTriggerpipimumu_magDw->Draw("SAME");
  relEffTriggerpipimumu_data->SetLineColor(kOrange);
  relEffTriggerpipimumu_data->Draw("SAME");
  relEffTriggerpipimumu_TISTOS_magUp->SetLineColor(kViolet);
  relEffTriggerpipimumu_TISTOS_magUp->Draw("SAME");
  
  a->cd(3);
  effTriggerKKmumu_KKmumuTrained_magUp->GetYaxis()->SetRangeUser(0.3,1);
  effTriggerKKmumu_KKmumuTrained_magUp->Draw();
  effTriggerKKmumu_KKmumuTrained_magDw->SetLineColor(kRed);
  effTriggerKKmumu_KKmumuTrained_magDw->Draw("SAME");
  a->cd(4);
  relEffTriggerKKmumu_magUp->GetYaxis()->SetRangeUser(0.5,1.5);
  relEffTriggerKKmumu_magUp->Draw();
  relEffTriggerKKmumu_magDw->SetLineColor(kRed);
  relEffTriggerKKmumu_magDw->Draw("SAME");
  a->cd(5);
  effTriggerKpimumu_KKmumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effTriggerKpimumu_KKmumuTrained_magUp->Draw();
  effTriggerKpimumu_KKmumuTrained_magDw->SetLineColor(kRed);
  effTriggerKpimumu_KKmumuTrained_magDw->Draw("SAME");
  a->cd(6);
  effTriggerKpimumu_pipimumuTrained_magUp->GetYaxis()->SetRangeUser(0.3,1.0);
  effTriggerKpimumu_pipimumuTrained_magUp->Draw();
  //effTriggerKpimumu_pipimumuTrained_magDw->SetLineColor(kRed);
  //effTriggerKpimumu_pipimumuTrained_magDw->Draw("SAME");
  
  effTriggerKpimumu_pipimumuTrained_data->SetLineColor(kOrange);
  effTriggerKpimumu_pipimumuTrained_data->Draw("SAME");

  effTriggerKpimumu_pipimumuTrained_TISTOS_magUp->SetLineColor(kViolet);
  effTriggerKpimumu_pipimumuTrained_TISTOS_magUp->Draw("SAME");

  a->Print(textFile+"5.eps");

  a->Write();
 

  
  TCanvas *c  = new TCanvas("Kpimumu","Eff Kpimumu");

  c->Divide(2);
  c->cd(1);
  effTriggerKpimumu_KKmumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effTriggerKpimumu_KKmumuTrained_magUp->Draw();
  effTriggerKpimumu_KKmumuTrained_magDw->SetLineColor(kRed);
  effTriggerKpimumu_KKmumuTrained_magDw->Draw("SAME");
  c->cd(2);
  effTriggerKpimumu_pipimumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effTriggerKpimumu_pipimumuTrained_magUp->Draw();
  effTriggerKpimumu_pipimumuTrained_magDw->SetLineColor(kRed);
  effTriggerKpimumu_pipimumuTrained_magDw->Draw("SAME");
  effTriggerKpimumu_pipimumuTrained_data->Draw("SAME");

  c->Print(textFile+"6.eps");

  c->Write();





}
			   
			   

			 
void MC_total_efficiency(double ghostProbCut,double hadronPID,double muonPID,double BDT){

  
  TString cut_hadronPID_Kpimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_pipimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNpi>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_KKmumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNk>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 

  TString cut_MuonPID=TString::Format("mu0_ProbNNghost<%.1f&&mu1_ProbNNghost<%.1f&&mu0_ProbNNmu>%.1f&&mu1_ProbNNmu>%.1f",ghostProbCut,ghostProbCut,muonPID,muonPID); 
  
  TString cut_BDT=TString::Format("BDT>%.1f",BDT);
  //TString cut_Trigger="D_L0MuonDecision_TOS";
  TString cut_Trigger="(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)";
  
  TString cut_PID_Kpimumu; 
  TString cut_PID_pipimumu;
  TString cut_PID_KKmumu ;

  TString cut_Kpimumu ;
  TString cut_pipimumu;
  TString cut_KKmumu ; 

  TString cut_norm_Kpimumu ;
  TString cut_norm_pipimumu;
  TString cut_norm_KKmumu ; 


  TString nameTarget = "MCL0Total_Efficiencies";
 
  cut_PID_Kpimumu = cut_hadronPID_Kpimumu + "&&" + cut_MuonPID;
  cut_PID_pipimumu = cut_hadronPID_pipimumu + "&&" + cut_MuonPID;
  cut_PID_KKmumu = cut_hadronPID_KKmumu + "&&" + cut_MuonPID;

  cut_Kpimumu=cut_PID_Kpimumu;
  cut_pipimumu=cut_PID_pipimumu;  
  cut_KKmumu=cut_PID_KKmumu; 

  nameTarget = nameTarget+"_afterPID";
  cut_Kpimumu = cut_Kpimumu+"&&"+cut_BDT;
  cut_pipimumu = cut_pipimumu+"&&"+cut_BDT;
  cut_KKmumu =   cut_KKmumu+"&&"+cut_BDT;
    
    
  nameTarget = nameTarget+"_afterBDT";
 
  cut_norm_Kpimumu = "";
  cut_norm_pipimumu ="";
  cut_norm_KKmumu = "";

  cut_Kpimumu = cut_Kpimumu + "&&" + cut_Trigger;
  cut_pipimumu =cut_pipimumu + "&&" + cut_Trigger;
  cut_KKmumu = cut_KKmumu + "&&" + cut_Trigger;

  nameTarget+=".root";
  

  TFile* fOut= new TFile(nameTarget,"RECREATE"); 

  TH1D* effTotalKpimumu_KKmumuTrained_magUp;
  TH1D* effTotalKpimumu_KKmumuTrained_magDw;

  effTotalKpimumu_KKmumuTrained_magUp= new TH1D("effTotalKpimumu_KKmumuTrained_magUp","efficiency MC Total D->Kpimumu magUp",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effTotalKpimumu_KKmumuTrained_magDw= new TH1D("effTotalKpimumu_KKmumuTrained_magDw","efficiency MC Total D->Kpimumu magDw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D* effTotalKpimumu_pipimumuTrained_magUp;
  TH1D* effTotalKpimumu_pipimumuTrained_magDw;

  effTotalKpimumu_pipimumuTrained_magUp= new TH1D("effTotalKpimumu_pipimumuTrained_magUp","efficiency MC Total D->Kpimumu magUp",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  effTotalKpimumu_pipimumuTrained_magDw= new TH1D("effTotalKpimumu_pipimumuTrained_magDw","efficiency MC Total D->Kpimumu magDw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D* effTotalKKmumu_KKmumuTrained_magUp;
  TH1D* effTotalKKmumu_KKmumuTrained_magDw;

  effTotalKKmumu_KKmumuTrained_magUp= new TH1D("effTotalKKmumu_KKmumuTrained_magUp","efficiency MC Total D->KKmumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  effTotalKKmumu_KKmumuTrained_magDw= new TH1D("effTotalKKmumu_KKmumuTrained_magDw","efficiency MC Total D->KKmumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1D* effTotalpipimumu_pipimumuTrained_magUp;
  TH1D* effTotalpipimumu_pipimumuTrained_magDw;

  effTotalpipimumu_pipimumuTrained_magUp= new TH1D("effTotalpipimumu_pipimumuTrained_magUp","efficiency MC Total D->pipimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  effTotalpipimumu_pipimumuTrained_magDw= new TH1D("effTotalKpipimumu_pipimumuTrained_magDw","efficiency MC Total D->pipimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1D* relEffTotalKKmumu_magUp;
  TH1D* relEffTotalKKmumu_magDw;
  TH1D* relEffTotalpipimumu_magUp;
  TH1D* relEffTotalpipimumu_magDw;

  relEffTotalKKmumu_magUp= new TH1D("relEffTotalDKKmumu_magUp","relative efficiency MC Total D->KKmumu/D->Kpimumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffTotalKKmumu_magDw= new TH1D("relEffTotalDKKmumu_magDw","relative efficiency MC Total D->KKmumu/D->Kpimumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);

  relEffTotalpipimumu_magUp= new TH1D("relEffTotalDpipimumu_magUp","relative efficiency MC Total D->pipimumu/D->Kpimumu magUp",sizeof(binspipi)/sizeof(double)-1,binspipi);
  relEffTotalpipimumu_magDw= new TH1D("relEffTotalDpipimumu_magDw","relative efficiency MC Total D->pipimumu/D->Kpimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  

  for(int i=0; i<rangesKpi_low.size();++i){
    
    TChain* tempChain_KKmumuTrained_magUp= new TChain("BDT_Tree");
    tempChain_KKmumuTrained_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
    TChain* tempChain_KKmumuTrained_magDw= new TChain("BDT_Tree");
    tempChain_KKmumuTrained_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));

    TChain* tempChain_pipimumuTrained_magUp= new TChain("BDT_Tree");
    tempChain_pipimumuTrained_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
    TChain* tempChain_pipimumuTrained_magDw= new TChain("BDT_Tree");
    tempChain_pipimumuTrained_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));


    //mag up KKmumu trained
    double nNorm = tempChain_KKmumuTrained_magUp->GetEntries(cut_norm_Kpimumu);
    double  nSel = tempChain_KKmumuTrained_magUp->GetEntries(cut_Kpimumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTotalKpimumu_KKmumuTrained_magUp->SetBinContent(i+1,tempEff);
    effTotalKpimumu_KKmumuTrained_magUp->SetBinError(i+1,dTempEff);

    //mag Dw KKmumu trained
    nNorm = tempChain_KKmumuTrained_magDw->GetEntries(cut_norm_Kpimumu);
    nSel = tempChain_KKmumuTrained_magDw->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTotalKpimumu_KKmumuTrained_magDw->SetBinContent(i+1,tempEff);
    effTotalKpimumu_KKmumuTrained_magDw->SetBinError(i+1,dTempEff);

    //mag up KKmumu trained
    nNorm = tempChain_pipimumuTrained_magUp->GetEntries(cut_norm_Kpimumu);
    nSel = tempChain_pipimumuTrained_magUp->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTotalKpimumu_pipimumuTrained_magUp->SetBinContent(i+1,tempEff);
    effTotalKpimumu_pipimumuTrained_magUp->SetBinError(i+1,dTempEff);

    //mag Dw pipimumu trained
    nNorm = tempChain_pipimumuTrained_magDw->GetEntries(cut_norm_Kpimumu);
    nSel = tempChain_pipimumuTrained_magDw->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTotalKpimumu_pipimumuTrained_magDw->SetBinContent(i+1,tempEff);
    effTotalKpimumu_pipimumuTrained_magDw->SetBinError(i+1,dTempEff);

    //mag dw
    nNorm = tempChain_pipimumuTrained_magUp->GetEntries(cut_norm_Kpimumu);
    nSel = tempChain_pipimumuTrained_magUp->GetEntries(cut_Kpimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTotalKpimumu_pipimumuTrained_magUp->SetBinContent(i+1,tempEff);
    effTotalKpimumu_pipimumuTrained_magUp->SetBinError(i+1,dTempEff);
    
     
  }    

  for(int i=0; i<rangesKK_low.size();++i){
    
    TChain* tempChain_magUp= new TChain("BDT_Tree");
    tempChain_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]));
    TChain* tempChain_magDw= new TChain("BDT_Tree");
    tempChain_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]));

    //mag up
    double nNorm = tempChain_magUp->GetEntries(cut_norm_KKmumu);
    double nSel = tempChain_magUp->GetEntries(cut_KKmumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTotalKKmumu_KKmumuTrained_magUp->SetBinContent(i+1,tempEff);
    effTotalKKmumu_KKmumuTrained_magUp->SetBinError(i+1,dTempEff);
    double relEff = tempEff/effTotalKpimumu_KKmumuTrained_magUp->GetBinContent(1);
    double dRelEff = TMath::Sqrt( (dTempEff/effTotalKKmumu_KKmumuTrained_magUp->GetBinContent(1))*(dTempEff/effTotalKKmumu_KKmumuTrained_magUp->GetBinContent(1)) 
				  + (tempEff*effTotalKKmumu_KKmumuTrained_magUp->GetBinError(1)/(effTotalKKmumu_KKmumuTrained_magUp->GetBinContent(1)*effTotalKKmumu_KKmumuTrained_magUp->GetBinContent(1)) * 
				     (tempEff*effTotalKKmumu_KKmumuTrained_magUp->GetBinError(1)/(effTotalKKmumu_KKmumuTrained_magUp->GetBinContent(1)*effTotalKKmumu_KKmumuTrained_magUp->GetBinContent(1))
				      ) ));
    relEffTotalKKmumu_magUp->SetBinContent(i+1,relEff);				     
    relEffTotalKKmumu_magUp->SetBinError(i+1,dRelEff);			     
    
    //mag Dw
    nNorm = tempChain_magDw->GetEntries(cut_norm_KKmumu);
    nSel = tempChain_magDw->GetEntries(cut_KKmumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTotalKKmumu_KKmumuTrained_magDw->SetBinContent(i+1,tempEff);
    effTotalKKmumu_KKmumuTrained_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effTotalKpimumu_KKmumuTrained_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effTotalKKmumu_KKmumuTrained_magDw->GetBinContent(1))*(dTempEff/effTotalKKmumu_KKmumuTrained_magDw->GetBinContent(1)) 
				  + (tempEff*effTotalKKmumu_KKmumuTrained_magDw->GetBinError(1)/(effTotalKKmumu_KKmumuTrained_magDw->GetBinContent(1)*effTotalKKmumu_KKmumuTrained_magDw->GetBinContent(1)) * 
				     (tempEff*effTotalKKmumu_KKmumuTrained_magDw->GetBinError(1)/(effTotalKKmumu_KKmumuTrained_magDw->GetBinContent(1)*effTotalKKmumu_KKmumuTrained_magDw->GetBinContent(1))
				      ) ));
    relEffTotalKKmumu_magDw->SetBinContent(i+1,relEff);				     
    relEffTotalKKmumu_magDw->SetBinError(i+1,dRelEff);			     
      

  }    

  for(int i=0; i<rangespipi_low.size();++i){
    
    TChain* tempChain_magUp= new TChain("BDT_Tree");
    tempChain_magUp ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]));
    TChain* tempChain_magDw= new TChain("BDT_Tree");
    tempChain_magDw ->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]));

    //mag up
    double nNorm = tempChain_magUp->GetEntries(cut_norm_pipimumu);
    double nSel = tempChain_magUp->GetEntries(cut_pipimumu);
    double tempEff = nSel/nNorm;
    double dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTotalpipimumu_pipimumuTrained_magUp->SetBinContent(i+1,tempEff);
    effTotalpipimumu_pipimumuTrained_magUp->SetBinError(i+1,dTempEff);
    double relEff = tempEff/effTotalKpimumu_pipimumuTrained_magUp->GetBinContent(1);
    double dRelEff = TMath::Sqrt( (dTempEff/effTotalpipimumu_pipimumuTrained_magUp->GetBinContent(1))*(dTempEff/effTotalpipimumu_pipimumuTrained_magUp->GetBinContent(1)) 
				  + (tempEff*effTotalpipimumu_pipimumuTrained_magUp->GetBinError(1)/(effTotalpipimumu_pipimumuTrained_magUp->GetBinContent(1)*effTotalpipimumu_pipimumuTrained_magUp->GetBinContent(1)) * 
				     (tempEff*effTotalpipimumu_pipimumuTrained_magUp->GetBinError(1)/(effTotalpipimumu_pipimumuTrained_magUp->GetBinContent(1)*effTotalpipimumu_pipimumuTrained_magUp->GetBinContent(1))
				      ) ));
    relEffTotalpipimumu_magUp->SetBinContent(i+1,relEff);				     
    relEffTotalpipimumu_magUp->SetBinError(i+1,dRelEff);			     
    
    //mag Dw
    nNorm = tempChain_magDw->GetEntries(cut_norm_pipimumu);
    nSel = tempChain_magDw->GetEntries(cut_pipimumu);
    tempEff = nSel/nNorm;
    dTempEff = 1/nNorm * TMath::Sqrt(nSel*(1-(nSel/nNorm) ) );
    effTotalpipimumu_pipimumuTrained_magDw->SetBinContent(i+1,tempEff);
    effTotalpipimumu_pipimumuTrained_magDw->SetBinError(i+1,dTempEff);
    relEff = tempEff/effTotalKpimumu_pipimumuTrained_magDw->GetBinContent(1);
    dRelEff = TMath::Sqrt( (dTempEff/effTotalpipimumu_pipimumuTrained_magDw->GetBinContent(1))*(dTempEff/effTotalpipimumu_pipimumuTrained_magDw->GetBinContent(1)) 
				  + (tempEff*effTotalpipimumu_pipimumuTrained_magDw->GetBinError(1)/(effTotalpipimumu_pipimumuTrained_magDw->GetBinContent(1)*effTotalpipimumu_pipimumuTrained_magDw->GetBinContent(1)) * 
				     (tempEff*effTotalpipimumu_pipimumuTrained_magDw->GetBinError(1)/(effTotalpipimumu_pipimumuTrained_magDw->GetBinContent(1)*effTotalpipimumu_pipimumuTrained_magDw->GetBinContent(1))
				      ) ));
    relEffTotalpipimumu_magDw->SetBinContent(i+1,relEff);				     
    relEffTotalpipimumu_magDw->SetBinError(i+1,dRelEff);			     
      

  }    


  //Draing part 

  TCanvas *a  = new TCanvas("pipimumu","Eff pipimumu");
  a->Divide(2,2);
  a->cd(1);
  effTotalpipimumu_pipimumuTrained_magUp->GetYaxis()->SetRangeUser(0.1,1);
  effTotalpipimumu_pipimumuTrained_magUp->Draw();
  effTotalpipimumu_pipimumuTrained_magDw->SetLineColor(kRed);
  effTotalpipimumu_pipimumuTrained_magDw->Draw("SAME");
  a->cd(2);
  relEffTotalpipimumu_magUp->GetYaxis()->SetRangeUser(0.7,1.8);
  relEffTotalpipimumu_magUp->Draw();
  relEffTotalpipimumu_magDw->SetLineColor(kRed);
  relEffTotalpipimumu_magDw->Draw("SAME");
    a->cd(3);
  effTotalKKmumu_KKmumuTrained_magUp->GetYaxis()->SetRangeUser(0.1,.3);
  effTotalKKmumu_KKmumuTrained_magUp->Draw();
  effTotalKKmumu_KKmumuTrained_magDw->SetLineColor(kRed);
  effTotalKKmumu_KKmumuTrained_magDw->Draw("SAME");
  a->cd(4);
  relEffTotalKKmumu_magUp->GetYaxis()->SetRangeUser(0.5,1.3);
  relEffTotalKKmumu_magUp->Draw();
  relEffTotalKKmumu_magDw->SetLineColor(kRed);
  relEffTotalKKmumu_magDw->Draw("SAME");

 
  a->Write();
  a->Print("test4.eps");
 

  
  TCanvas *c  = new TCanvas("Kpimumu","Eff Kpimumu");

  c->Divide(2);
  c->cd(1);
  effTotalKpimumu_KKmumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effTotalKpimumu_KKmumuTrained_magUp->Draw();
  effTotalKpimumu_KKmumuTrained_magDw->SetLineColor(kRed);
  effTotalKpimumu_KKmumuTrained_magDw->Draw("SAME");
  c->cd(2);
  effTotalKpimumu_pipimumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effTotalKpimumu_pipimumuTrained_magUp->Draw();
  effTotalKpimumu_pipimumuTrained_magDw->SetLineColor(kRed);
  effTotalKpimumu_pipimumuTrained_magDw->Draw("SAME");


  c->Write();

}
			   
void plotPIDCalibMuonSampleAndSignalMC() {
  
  TChain* chain_PIDCalib = new TChain("JpsiCalib");
  chain_PIDCalib->AddFile("/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/tree.root");
  
  TChain* chain_MC = new TChain("MC12_DstD2KKMuMu/DecayTree");
  chain_MC->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_withGenLevel/magUp/MC12_DstD2KKmumu_magUp.root");
  
  TFile* fOut = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/MC_PIDCalib_P_ETA_MuonSample.root","RECREATE");

  double nsig_sw;
  int Dst_BKGCAT;
  double p,pt;
  double eta; 
  bool isMuon;
  bool inMuonAcc;
  bool D_L0MuonDecision;
  chain_PIDCalib->SetBranchAddress("nsig_sw",&nsig_sw);
  chain_PIDCalib->SetBranchAddress("Mu_P",&p);
  chain_PIDCalib->SetBranchAddress("Mu_PT",&pt);
  chain_PIDCalib->SetBranchAddress("Mu_Eta",&eta);
  //chain_PIDCalib->SetBranchAddress("Mu_IsMuon",&isMuon);
  //chain_PIDCalib->SetBranchAddress("Mu_InMuonAcc",&inMuonAcc);

  chain_MC->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_MC->SetBranchAddress("mu0_P",&p);
  chain_MC->SetBranchAddress("mu0_PT",&pt);
  chain_MC->SetBranchAddress("mu0_TRACK_Eta",&eta);
  chain_MC->SetBranchAddress("D_L0MuonDecision_TOS",&D_L0MuonDecision);


  TH2* P_ETA_PIDCalib = new TH2D("P_ETA_PIDCalib","P_ETA_PIDCalib",30,0,100000,30,1.5,5);
  TH2* P_ETA_MC = new TH2D("P_ETA_MC","P_ETA_MC",30,0,100000,30,1.5,5);

  TH2* P_PT_PIDCalib = new TH2D("P_PT_PIDCalib","P_PT_PIDCalib",30,0,100000,30,0,8000);
  TH2* P_PT_MC = new TH2D("P_PT_MC","PT_ETA_MC",30,0,100000,30,0,8000);

  std::vector<TH2D*> P_ETA_MC_PtBins;
  std::vector<TH2D*> P_ETA_PIDCalib_PtBins;

  TH2D* temp;

  for(int j=0; j<rangesPt_low.size();++j){                                                                      
    temp = new TH2D(TString::Format("MC_p_eta_%.1f_1f.root",rangesPt_low[j],rangesPt_high[j]),TString::Format("p_eta_%.1f_%.1f.root",rangesPt_low[j],rangesPt_high[j]),30,0,100000,30,1.5,5);  
    P_ETA_MC_PtBins.push_back(temp);
    temp = new TH2D(TString::Format("PIDCalib_p_eta_%.1f_1f.root",rangesPt_low[j],rangesPt_high[j]),TString::Format("p_eta_%.1f_%1f.root",rangesPt_low[j],rangesPt_high[j]),30,0,100000,30,1.5,5);  
    P_ETA_PIDCalib_PtBins.push_back(temp);
  }

 
  for (int i=0;i<chain_PIDCalib->GetEntries();++i) {

    chain_PIDCalib->GetEntry(i);
    //if(!isMuon) continue;
    //if(!inMuonAcc) continue;

    P_ETA_PIDCalib->Fill(p,eta,nsig_sw);
    P_PT_PIDCalib->Fill(p,pt,nsig_sw);

    for(int j=0; j<rangesPt_low.size();++j){                                                                      
      if(pt>rangesPt_low[j] && pt<rangesPt_high[j]) P_ETA_PIDCalib_PtBins[j]->Fill(p,eta,nsig_sw);
    }

  }

    // std::cout<<p<<"  "<<eta<<"  "<<std::endl;
  
  std::cout<<"part 1 done"<<std::endl;
  for (int i=0;i<chain_MC->GetEntries();++i) {

    chain_MC->GetEntry(i);    
    //if(i>40000) break;
    //if(!D_L0MuonDecision) continue;
    if(Dst_BKGCAT<11) { 
     P_ETA_MC->Fill(p,eta);
     P_PT_MC->Fill(p,pt);

     for(int j=0; j<rangesPt_low.size();++j){                                                                      
       if(pt>rangesPt_low[j] && pt<rangesPt_high[j]) P_ETA_MC_PtBins[j]->Fill(p,eta);
     }


    }  
     //    std::cout<<i<<"  "<<p<<"  "<<eta<<"  "<<std::endl;
}

  dcastyle();

  TCanvas *b  = new TCanvas("b","b");
  //b->Divide(3,2);
  b->cd(1);

  int colors[]={1,2,800,8,6,kCyan};

  for(int j=0; j<rangesPt_low.size();++j){                                                                      
    P_ETA_MC_PtBins[j]->SetMarkerColor(colors[j]);
    P_ETA_MC_PtBins[j]->SetMarkerSize(2);
    P_ETA_MC_PtBins[j]->SetLineWidth(2);
    P_ETA_MC_PtBins[j]->Draw("same"); 
    std::cout<<j<<" "<<rangesPt_low[j]<<" -  "<<rangesPt_high[j] <<"  "<<colors[j]<<"  "<<P_ETA_MC_PtBins[j]->GetEntries() <<std::endl;
   }
  P_ETA_PIDCalib->SetLineWidth(2);
  P_ETA_PIDCalib->SetLineColor(kBlue);
  P_ETA_PIDCalib->Draw("box same"); 
  
  b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/MC_PIDCalib_P_ETA_PTBINS_MuonSample.eps");

  TCanvas *c  = new TCanvas("c","c");
  c->Divide(3,2);
  c->cd(1);

  //P_ETA_PIDCalib->Draw("colz");
  //c->cd(2);
  //P_ETA_MC->Draw("");
  //c->cd(3);
  P_ETA_PIDCalib->SetLineWidth(2);
  P_ETA_PIDCalib->SetLineColor(kBlue);
  P_ETA_MC->Draw("");
  P_ETA_PIDCalib->Draw("box same");
  c->cd(2);
  P_ETA_MC->ProjectionX()->DrawNormalized("");
  P_ETA_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  c->cd(3);
  P_ETA_MC->ProjectionY()->DrawNormalized("");
  P_ETA_PIDCalib->ProjectionY()->DrawNormalized("SAME");

  c->cd(4);
  P_PT_PIDCalib->SetLineWidth(2);
  P_PT_PIDCalib->SetLineColor(kBlue);
  P_PT_MC->Draw("");
  P_PT_PIDCalib->Draw("box same");
  c->cd(5);
  P_PT_MC->ProjectionX()->DrawNormalized();
  P_PT_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  c->cd(6);
  P_PT_MC->ProjectionY()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionY()->DrawNormalized("SAME");
  
  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/MC_PIDCalib_P_ETA_MuonSample.eps");
  fOut->Write();
}


			 


void plotPIDCalibKaonSampleAndSignalMC() {
  
  TChain* chain_PIDCalib = new TChain("RSDStCalib");
  chain_PIDCalib->AddFile("/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/PIDCalibTree_Kaon.root");
  
  TChain* chain_MC = new TChain("MC12_DstD2KKMuMu/DecayTree");
  chain_MC->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2KKmumu_withGenLevel/magUp/MC12_DstD2KKmumu_magUp.root");
  TChain* chain_MC_norm = new TChain("MC12_DstD2KPiMuMu/DecayTree");
  chain_MC_norm->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu_withGenLevel/magUp/MC12_DstD2Kpimumu_magUp.root");

  TFile* fOut = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/MC_PIDCalib_P_ETA_KaonSample.root","RECREATE");
  

  double nsig_sw;
  int Dst_BKGCAT;
  double p,pt;
  double eta; 
  bool isMuon;
  bool inMuonAcc;
  bool D_L0MuonDecision;
  chain_PIDCalib->SetBranchAddress("nsig_sw",&nsig_sw);
  chain_PIDCalib->SetBranchAddress("K_P",&p);
  chain_PIDCalib->SetBranchAddress("K_PT",&pt);
  chain_PIDCalib->SetBranchAddress("K_Eta",&eta);
  //chain_PIDCalib->SetBranchAddress("K_IsKon",&isKon);
  //chain_PIDCalib->SetBranchAddress("Mu_InMuonAcc",&inMuonAcc);

  chain_MC->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_MC->SetBranchAddress("h0_P",&p);
  chain_MC->SetBranchAddress("h0_PT",&pt);
  chain_MC->SetBranchAddress("h0_TRACK_Eta",&eta);
  chain_MC->SetBranchAddress("D_L0MuonDecision_TOS",&D_L0MuonDecision);

  chain_MC_norm->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_MC_norm->SetBranchAddress("h0_P",&p);
  chain_MC_norm->SetBranchAddress("h0_PT",&pt);
  chain_MC_norm->SetBranchAddress("h0_TRACK_Eta",&eta);
  chain_MC_norm->SetBranchAddress("D_L0MuonDecision_TOS",&D_L0MuonDecision);


  TH2* P_ETA_PIDCalib = new TH2D("P_ETA_PIDCalib","P_ETA_PIDCalib",30,0,100000,30,1.5,5);
  TH2* P_ETA_MC = new TH2D("P_ETA_MC","P_ETA_MC",30,0,100000,30,1.5,5);
  TH2* P_ETA_MC_norm = new TH2D("P_ETA_MC_norm","P_ETA_MC_norm",30,0,100000,30,1.5,5);

  TH2* P_PT_PIDCalib = new TH2D("P_PT_PIDCalib","P_PT_PIDCalib",30,0,100000,30,0,8000);
  TH2* P_PT_MC = new TH2D("P_PT_MC","PT_ETA_MC",30,0,100000,30,0,8000);
  TH2* P_PT_MC_norm = new TH2D("P_PT_MC_norm","PT_ETA_MC_norm",30,0,100000,30,0,8000);

 
  for (int i=0;i<chain_PIDCalib->GetEntries();++i) {

    chain_PIDCalib->GetEntry(i);
    //if(!isMuon) continue;
    //if(!inMuonAcc) continue;

    P_ETA_PIDCalib->Fill(p,eta,nsig_sw);
    P_PT_PIDCalib->Fill(p,pt,nsig_sw);

  }

    // std::cout<<p<<"  "<<eta<<"  "<<std::endl;
  
  std::cout<<"part 1 done"<<std::endl;
  for (int i=0;i<chain_MC->GetEntries();++i) {

    chain_MC->GetEntry(i);    
    //if(i>40000) break;
    //if(!D_L0MuonDecision) continue;
    if(Dst_BKGCAT<11) { 
     P_ETA_MC->Fill(p,eta);
     P_PT_MC->Fill(p,pt);
    }  
  }
  for (int i=0;i<chain_MC_norm->GetEntries();++i) {

    chain_MC_norm->GetEntry(i);    
    //if(i>40000) break;
    //if(!D_L0MuonDecision) continue;
    if(Dst_BKGCAT<11) { 
     P_ETA_MC_norm->Fill(p,eta);
     P_PT_MC_norm->Fill(p,pt);
    }  
  }

  dcastyle();

  TCanvas *c  = new TCanvas("c","c");
  c->Divide(3,2);
  c->cd(1);

  P_ETA_PIDCalib->SetLineWidth(2);
  P_ETA_PIDCalib->SetLineColor(kBlue);
  P_ETA_MC->Draw("");
  P_ETA_PIDCalib->Draw("box same");

  c->cd(2);
  P_ETA_MC->ProjectionX()->DrawNormalized("");
  P_ETA_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  c->cd(3);
  P_ETA_MC->ProjectionY()->DrawNormalized("");
  P_ETA_PIDCalib->ProjectionY()->DrawNormalized("SAME");

  c->cd(4);
  P_ETA_PIDCalib->SetLineWidth(2);
  P_ETA_PIDCalib->SetLineColor(kBlue);
  P_ETA_MC_norm->Draw("");
  P_ETA_PIDCalib->Draw("box same");
  c->cd(5);
  P_ETA_MC_norm->ProjectionX()->DrawNormalized("");
  P_ETA_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  c->cd(6);
  P_ETA_MC_norm->ProjectionY()->DrawNormalized("");
  P_ETA_PIDCalib->ProjectionY()->DrawNormalized("SAME");
  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/MC_PIDCalib_P_ETA_KaonSample.eps");

  TCanvas *d  = new TCanvas("d","d");
  d->Divide(3,2);
  d->cd(1);

  P_PT_PIDCalib->SetLineWidth(2);
  P_PT_PIDCalib->SetLineColor(kBlue);
  P_PT_MC->Draw("");
  P_PT_PIDCalib->Draw("box same");

  d->cd(2);
  P_PT_MC->ProjectionX()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  c->cd(3);
  P_PT_MC->ProjectionY()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionY()->DrawNormalized("SAME");

  d->cd(4);
  P_PT_PIDCalib->SetLineWidth(2);
  P_PT_PIDCalib->SetLineColor(kBlue);
  P_PT_MC_norm->Draw("");
  P_PT_PIDCalib->Draw("box same");
  d->cd(5);
  P_PT_MC_norm->ProjectionX()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  d->cd(6);
  P_PT_MC_norm->ProjectionY()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionY()->DrawNormalized("SAME");


  
  d->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/norm_MC_PIDCalib_P_PT_KaonSample.eps");
  fOut->Write();
  //d->Print("test2.eps");

}

void plotPIDCalibPionSampleAndSignalMC() {
  
  TChain* chain_PIDCalib = new TChain("RSDStCalib");
  chain_PIDCalib->AddFile("/home/he/mitzel/cmtuser/Urania_v5r0/PIDCalib/PIDPerfScripts/scripts/D2hhmumu/PIDCalibTree_Pion.root");
  
  TChain* chain_MC = new TChain("MC12_DstD2PiPiMuMu/DecayTree");
  chain_MC->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2pipimumu_withGenLevel/magUp/MC12_DstD2pipimumu_magUp.root");
  TChain* chain_MC_norm = new TChain("MC12_DstD2KPiMuMu/DecayTree");
  chain_MC_norm->AddFile("/auto/data/mitzel/D2hhmumu/new/MC/D2Kpimumu_withGenLevel/magUp/MC12_DstD2Kpimumu_magUp.root");
  
  TFile* fOut = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/MC_PIDCalib_P_ETA_PionSample.root","RECREATE");

  double nsig_sw;
  int Dst_BKGCAT;
  double p,pt,slow_pt, slow_p, slow_eta;
  double eta; 
  bool isMuon;
  bool inMuonAcc;
  bool D_L0MuonDecision;
  chain_PIDCalib->SetBranchAddress("nsig_sw",&nsig_sw);
  chain_PIDCalib->SetBranchAddress("Pi_P",&p);
  chain_PIDCalib->SetBranchAddress("Pi_PT",&pt);
  chain_PIDCalib->SetBranchAddress("Pi_Eta",&eta);
  //chain_PIDCalib->SetBranchAddress("K_IsKon",&isKon);
  //chain_PIDCalib->SetBranchAddress("Mu_InMuonAcc",&inMuonAcc);

  chain_MC->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_MC->SetBranchAddress("h0_P",&p);
  chain_MC->SetBranchAddress("h0_PT",&pt);
  chain_MC->SetBranchAddress("h0_TRACK_Eta",&eta);
  chain_MC->SetBranchAddress("D_L0MuonDecision_TOS",&D_L0MuonDecision);

  chain_MC_norm->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_MC_norm->SetBranchAddress("h1_P",&p);
  chain_MC_norm->SetBranchAddress("h1_PT",&pt);
  chain_MC_norm->SetBranchAddress("h1_TRACK_Eta",&eta);
  chain_MC_norm->SetBranchAddress("D_L0MuonDecision_TOS",&D_L0MuonDecision);

  chain_MC_norm->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_MC_norm->SetBranchAddress("Slowpi_P",&slow_p);
  chain_MC_norm->SetBranchAddress("Slowpi_PT",&slow_pt);
  chain_MC_norm->SetBranchAddress("Slowpi_TRACK_Eta",&slow_eta);
  chain_MC_norm->SetBranchAddress("D_L0MuonDecision_TOS",&D_L0MuonDecision);

  TH2* P_ETA_PIDCalib = new TH2D("P_ETA_PIDCalib","P_ETA_PIDCalib",30,0,100000,30,1.5,5);
  TH2* P_ETA_MC = new TH2D("P_ETA_MC","P_ETA_MC",30,0,100000,30,1.5,5);
  TH2* P_ETA_MC_norm = new TH2D("P_ETA_MC_norm","P_ETA_MC_norm",30,0,100000,30,1.5,5);
  TH2* P_ETA_MC_slowpi = new TH2D("P_ETA_MC_slowpi","P_ETA_MC_slowpi",30,0,100000,30,1.5,5);

  TH2* P_PT_PIDCalib = new TH2D("P_PT_PIDCalib","P_PT_PIDCalib",30,0,100000,30,0,8000);
  TH2* P_PT_MC = new TH2D("P_PT_MC","PT_ETA_MC",30,0,100000,30,0,8000);
  TH2* P_PT_MC_norm = new TH2D("P_PT_MC_norm","PT_ETA_MC_norm",30,0,100000,30,0,8000);
  TH2* P_PT_MC_slowpi = new TH2D("P_PT_MC_slowpi","PT_ETA_MC_slowpi",30,0,100000,30,0,8000);

 
  for (int i=0;i<chain_PIDCalib->GetEntries();++i) {

    chain_PIDCalib->GetEntry(i);
    //if(!isMuon) continue;
    //if(!inMuonAcc) continue;

    P_ETA_PIDCalib->Fill(p,eta,nsig_sw);
    P_PT_PIDCalib->Fill(p,pt,nsig_sw);

  }

    // std::cout<<p<<"  "<<eta<<"  "<<std::endl;
  
  std::cout<<"part 1 done"<<std::endl;
  for (int i=0;i<chain_MC->GetEntries();++i) {

    chain_MC->GetEntry(i);    
    //if(i>40000) break;
    //if(!D_L0MuonDecision) continue;
    if(Dst_BKGCAT<11) { 
     P_ETA_MC->Fill(p,eta);
     P_PT_MC->Fill(p,pt);
    }  
  }
  for (int i=0;i<chain_MC_norm->GetEntries();++i) {

    chain_MC_norm->GetEntry(i);    
    //if(i>40000) break;
    //if(!D_L0MuonDecision) continue;
    if(Dst_BKGCAT<11) { 
     P_ETA_MC_norm->Fill(p,eta);
     P_ETA_MC_slowpi->Fill(slow_p,slow_eta);
     P_PT_MC_norm->Fill(p,pt);
     P_PT_MC_slowpi->Fill(slow_p,slow_pt);
    }  
  }

  dcastyle();

  TCanvas *c  = new TCanvas("c","c");
  c->Divide(3,2);
  c->cd(1);

  P_ETA_PIDCalib->SetLineWidth(2);
  P_ETA_PIDCalib->SetLineColor(kBlue);
  P_ETA_MC->Draw("");
  P_ETA_PIDCalib->Draw("box same");

  c->cd(2);
  P_ETA_MC->ProjectionX()->DrawNormalized("");
  P_ETA_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  c->cd(3);
  P_ETA_MC->ProjectionY()->DrawNormalized("");
  P_ETA_PIDCalib->ProjectionY()->DrawNormalized("SAME");

  c->cd(4);
  P_ETA_PIDCalib->SetLineWidth(2);
  P_ETA_PIDCalib->SetLineColor(kBlue);
  P_ETA_MC_norm->Draw("");
  P_ETA_PIDCalib->Draw("box same");
  c->cd(5);
  P_ETA_MC_norm->ProjectionX()->DrawNormalized("");
  P_ETA_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  c->cd(6);
  P_ETA_MC_norm->ProjectionY()->DrawNormalized("");
  P_ETA_PIDCalib->ProjectionY()->DrawNormalized("SAME");
  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/MC_PIDCalib_P_ETA_PionSample.eps");

  TCanvas *d  = new TCanvas("d","d");
  d->Divide(3,2);
  d->cd(1);

  P_PT_PIDCalib->SetLineWidth(2);
  P_PT_PIDCalib->SetLineColor(kBlue);
  P_PT_MC->Draw("");
  P_PT_PIDCalib->Draw("box same");

  d->cd(2);
  P_PT_MC->ProjectionX()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  c->cd(3);
  P_PT_MC->ProjectionY()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionY()->DrawNormalized("SAME");

  d->cd(4);
  P_PT_PIDCalib->SetLineWidth(2);
  P_PT_PIDCalib->SetLineColor(kBlue);
  P_PT_MC_norm->Draw("");
  P_PT_PIDCalib->Draw("box same");
  d->cd(5);
  P_PT_MC_norm->ProjectionX()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  d->cd(6);
  P_PT_MC_norm->ProjectionY()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionY()->DrawNormalized("SAME");
  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/norm_MC_PIDCalib_P_ETA_PionSample.eps");
  

  TCanvas *e  = new TCanvas("e","e");
  e->Divide(3,2);
  e->cd(1);

  P_PT_PIDCalib->SetLineWidth(2);
  P_PT_PIDCalib->SetLineColor(kBlue);
  P_PT_MC->Draw("");
  P_PT_PIDCalib->Draw("box same");

  e->cd(2);
  P_PT_MC->ProjectionX()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  e->cd(3);
  P_PT_MC->ProjectionY()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionY()->DrawNormalized("SAME");

  e->cd(4);
  P_PT_PIDCalib->SetLineWidth(2);
  P_PT_PIDCalib->SetLineColor(kBlue);
  P_PT_MC_slowpi->Draw("");
  P_PT_PIDCalib->Draw("box same");
  e->cd(5);
  P_PT_MC_slowpi->ProjectionX()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionX()->DrawNormalized("SAME");
  e->cd(6);
  P_PT_MC_slowpi->ProjectionY()->DrawNormalized("");
  P_PT_PIDCalib->ProjectionY()->DrawNormalized("SAME");

  e->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/slowpi_MC_PIDCalib_P_PT_KaonSample.eps");
  fOut->Write();
  
}

 
TH1D*  fillPtSpectrum(int muonIndex,TString fIn,TString fIn2,TString Tree,TString nameHisto) {
  
  double pt;
  
  TChain * chain = new TChain(Tree);
  chain->AddFile(fIn);
  chain->AddFile(fIn2);
  chain->SetBranchAddress(TString::Format("mu%i_PT",muonIndex),&pt);
  
  TH1D* h_pt  = new TH1D(nameHisto,"pt",50,0,8000);
  for(int i=0;i<chain->GetEntries();++i){
    chain->GetEntry(i);
    h_pt->Fill(pt);
    //  std::cout<<i<<"  "<<pt<<std::endl;
  }
  return h_pt;
}

TH1D*  fillPtSpectrum2Bins(int muonIndex,TString fIn,TString fIn2,TString Tree,TString nameHisto) {
  
  double pt;
  
  TChain * chain = new TChain(Tree);
  chain->AddFile(fIn);
  chain->AddFile(fIn2);
  chain->SetBranchAddress(TString::Format("mu%i_PT",muonIndex),&pt);

  double bins[]={0,800,8000};
  TH1D* h_pt  = new TH1D(nameHisto,"pt",2,bins);
  for(int i=0;i<chain->GetEntries();++i){
    chain->GetEntry(i);
    h_pt->Fill(pt);
  }
  return h_pt;
}

			 
void  drawPtSpectra(int muonIndex) {

  dcastyle();
  TFile *fout = new TFile(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/PtSpectrum_Muon%i_MC.root",muonIndex),"RECREATE");

  std::vector<TH1D*> ptSpectrum_KKmumu;
  std::vector<TH1D*> ptSpectrum_Kpimumu;
  std::vector<TH1D*> ptSpectrum_pipimumu;

  std::vector<TH1D*> ptSpectrum_KKmumu_2Bins;
  std::vector<TH1D*> ptSpectrum_Kpimumu_2Bins;
  std::vector<TH1D*> ptSpectrum_pipimumu_2Bins;

  TString MC_file_dw;
  TString MC_file_up;
  
  //norm mode
  for(int i=0; i<rangesKpi_low.size();++i){
    MC_file_up = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_",rangesKpi_low[0],rangesKpi_high[0])+"magUp.root";
    MC_file_dw = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_",rangesKpi_low[0],rangesKpi_high[0])+"magDw.root";
    TH1D* temp = fillPtSpectrum(muonIndex,MC_file_dw,MC_file_up,"BDT_Tree",TString::Format("ptSpectrum_Kpimumu_q2bin_%i_mu%i",i,muonIndex));
    ptSpectrum_Kpimumu.push_back(temp);
    ptSpectrum_Kpimumu_2Bins.push_back(fillPtSpectrum2Bins(muonIndex,MC_file_dw,MC_file_up,"BDT_Tree",TString::Format("ptSpectrum_2Bins_Kpimumu_q2bin_%i_mu%i",i,muonIndex))); 
  }

  //KKmumu mode
  for(int i=0; i<rangesKK_low.size();++i){
    //MC 
    MC_file_up = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_",rangesKK_low[i],rangesKK_high[i])+"magUp.root";
    MC_file_dw = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_",rangesKK_low[i],rangesKK_high[i])+"magDw.root";
    ptSpectrum_KKmumu.push_back(fillPtSpectrum(0,MC_file_dw,MC_file_up,"BDT_Tree",TString::Format("ptSpectrum_KKmumu_q2bin_%i_mu%i",i,muonIndex)));
    ptSpectrum_KKmumu_2Bins.push_back(fillPtSpectrum2Bins(0,MC_file_dw,MC_file_up,"BDT_Tree",TString::Format("ptSpectrum_2Bins_KKmumu_q2bin_%i_mu%i",i,muonIndex))); 
  }

  //pipimumu mode
  for(int i=0; i<rangespipi_low.size();++i){
    //MC 
    MC_file_up = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_",rangespipi_low[i],rangespipi_high[i])+"magUp.root";
    MC_file_dw = TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_",rangespipi_low[i],rangespipi_high[i])+"magDw.root";
    ptSpectrum_pipimumu.push_back(fillPtSpectrum(muonIndex,MC_file_dw,MC_file_up,"BDT_Tree",TString::Format("ptSpectrum_pipimumu_q2bin_%i_mu%i",i,muonIndex)));
    ptSpectrum_pipimumu_2Bins.push_back(fillPtSpectrum2Bins(muonIndex,MC_file_dw,MC_file_up,"BDT_Tree",TString::Format("ptSpectrum_2Bins_pipimumu_q2bin_%i_mu%i",i,muonIndex))); 
   }
  TCanvas *d = new TCanvas("d","d");
  d->Divide(2);

    
  for(int i=0; i<rangesKK_low.size();++i){                                                                                                                                     
    d->cd(1);                                                                                                                                                              
    ptSpectrum_KKmumu[i]->SetLineColor(i+1);
    ptSpectrum_KKmumu[i]->DrawNormalized("SAME");
    ptSpectrum_KKmumu[i]->Write();
    d->cd(2);                                                                                                                                                              
    ptSpectrum_KKmumu_2Bins[i]->SetLineColor(i+1);
    ptSpectrum_KKmumu_2Bins[i]->DrawNormalized("SAME");
    ptSpectrum_KKmumu_2Bins[i]->Write();
  }
  
  TCanvas *e = new TCanvas("e","e");
  e->Divide(2);
  for(int i=0; i<rangesKpi_low.size();++i){                                                                                                                                     
    e->cd(1);                                                                                                                                                            
    ptSpectrum_Kpimumu[i]->SetLineColor(i+1);
    ptSpectrum_Kpimumu[i]->DrawNormalized("SAME");
    e->cd(2);                                                                                                                                                              
    ptSpectrum_KKmumu_2Bins[i]->SetLineColor(i+1);
    ptSpectrum_KKmumu_2Bins[i]->DrawNormalized("SAME");
  }
  TCanvas *f = new TCanvas("f","f");
  f->Divide(2);

  for(int i=0; i<rangespipi_low.size();++i){                                                                                                                                     
    f->cd(1);                                                                                                                                                             
    ptSpectrum_pipimumu[i]->SetLineColor(i+1);
    ptSpectrum_pipimumu[i]->DrawNormalized("SAME");
    ptSpectrum_pipimumu[i]->Write();
    f->cd(2);                                                                                                                                                              
    ptSpectrum_pipimumu_2Bins[i]->SetLineColor(i+1);
    ptSpectrum_pipimumu_2Bins[i]->DrawNormalized("SAME");
    ptSpectrum_pipimumu_2Bins[i]->Write();
  }
  ptSpectrum_Kpimumu_2Bins[0]->Write();
  
  fout->Close();
  d->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/ptSectrum_mu%i_KKmumu.eps",muonIndex));
  e->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/ptSectrum_mu%i_Kpimumu.eps",muonIndex));
  f->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/ptSectrum_mu%i_pipimumu.eps",muonIndex));


}


void drawGenLevelCutEfficiencies(){

 
  TChain *chain_KK_beforeCut_up = new TChain("MCTree_beforecut");
  TChain *chain_KK_afterCut_up = new TChain("MCTree_aftercut");
  TChain *chain_KK_beforeCut_dw = new TChain("MCTree_beforecut");
  TChain *chain_KK_afterCut_dw = new TChain("MCTree_aftercut");

  chain_KK_beforeCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/KKmumu/genLevelTuple_27215002_up.root");
  chain_KK_afterCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/KKmumu/genLevelTuple_27215002_up.root");
  chain_KK_beforeCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/KKmumu/genLevelTuple_27215002_dw.root");
  chain_KK_afterCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/KKmumu/genLevelTuple_27215002_dw.root");

  chain_KK_beforeCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/KKmumu/genLevelTuple_27215002_up_2.root");
  chain_KK_afterCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/KKmumu/genLevelTuple_27215002_up_2.root");
  chain_KK_beforeCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/KKmumu/genLevelTuple_27215002_dw_2.root");
  chain_KK_afterCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/KKmumu/genLevelTuple_27215002_dw_2.root");

  TChain *chain_pipi_beforeCut_up = new TChain("MCTree_beforecut");
  TChain *chain_pipi_afterCut_up = new TChain("MCTree_aftercut");
  TChain *chain_pipi_beforeCut_dw = new TChain("MCTree_beforecut");
  TChain *chain_pipi_afterCut_dw = new TChain("MCTree_aftercut");

  chain_pipi_beforeCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/pipimumu/genLevelTuple_27215003_up.root");
  chain_pipi_afterCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/pipimumu/genLevelTuple_27215003_up.root");
  chain_pipi_beforeCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/pipimumu/genLevelTuple_27215003_dw.root");
  chain_pipi_afterCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/pipimumu/genLevelTuple_27215003_dw.root");

  chain_pipi_beforeCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/pipimumu/genLevelTuple_27215003_up_2.root");
  chain_pipi_afterCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/pipimumu/genLevelTuple_27215003_up_2.root");
  chain_pipi_beforeCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/pipimumu/genLevelTuple_27215003_dw_2.root");
  chain_pipi_afterCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/pipimumu/genLevelTuple_27215003_dw_2.root");

  TChain *chain_Kpi_beforeCut_up = new TChain("MCTree_beforecut");
  TChain *chain_Kpi_afterCut_up = new TChain("MCTree_aftercut");
  TChain *chain_Kpi_beforeCut_dw = new TChain("MCTree_beforecut");
  TChain *chain_Kpi_afterCut_dw = new TChain("MCTree_aftercut");

  chain_Kpi_beforeCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/Kpimumu/genLevelTuple_27215000_up.root");
  chain_Kpi_afterCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/Kpimumu/genLevelTuple_27215000_up.root");
  chain_Kpi_beforeCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/Kpimumu/genLevelTuple_27215000_dw.root");
  chain_Kpi_afterCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/genLevelMC/Kpimumu/genLevelTuple_27215000_dw.root");

  dcastyle();

  TFile *fOut= new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/GenLevelCutEfficiency/generatorLevelCutEfficiency.root","RECREATE");
  
  TH1D * genLevelEff_KK_up = new TH1D("genLevelEff_KK_up","genLevelEff_KK_up",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * genLevelEff_KK_dw = new TH1D("genLevelEff_KK_dw","genLevelEff_KK_up",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * genLevelEff_KK = new TH1D("genLevelEff_KK","genLevelEff_KK",sizeof(binsKK)/sizeof(double)-1,binsKK);
  
  TH1D * genLevelEff_pipi_up = new TH1D("genLevelEff_pipi_up","genLevelEff_pipi_up",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * genLevelEff_pipi_dw = new TH1D("genLevelEff_pipi_dw","genLevelEff_pipi_up",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * genLevelEff_pipi = new TH1D("genLevelEff_pipi","genLevelEff_pipi",sizeof(binspipi)/sizeof(double)-1,binspipi);
  
  TH1D * genLevelEff_Kpi_up = new TH1D("genLevelEff_Kpi_up","genLevelEff_Kpi_up",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D * genLevelEff_Kpi_dw = new TH1D("genLevelEff_Kpi_dw","genLevelEff_Kpi_up",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D * genLevelEff_Kpi = new TH1D("genLevelEff_Kpi","genLevelEff_Kpi",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D * genLevelRelativeEff_KK = new TH1D("genLevelRelativeEff_KK","genLevelRelativeEff_KK",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * genLevelRelativeEff_pipi = new TH1D("genLevelRelativeEff_pipi","genLevelRelativeEff_pipi",sizeof(binspipi)/sizeof(double)-1,binspipi);


  //D2KKmumu

  //magUp

  double D_DiMuon_Mass;
  double Eff, dEff;

  double nTot_up, nSel_up; 
  double nTot_dw, nSel_dw; 
  double nTot, nSel; 

  for (int j=0;j<sizeof(binsKpi)/sizeof(double)-1;++j) {
    
    nTot_up = (double)chain_Kpi_beforeCut_up->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKpi_low[j],rangesKpi_high[j] ) );
    nSel_up = (double)chain_Kpi_afterCut_up->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKpi_low[j],rangesKpi_high[j] ) );
    Eff=nSel_up/nTot_up;
    dEff=(1/nTot_up) * TMath::Sqrt(nSel_up * (1 - (nSel_up/nTot_up) ) );
    genLevelEff_Kpi_up->SetBinContent(j+1,Eff);
    genLevelEff_Kpi_up->SetBinError(j+1,dEff);

    nTot_dw = (double)chain_Kpi_beforeCut_dw->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKpi_low[j],rangesKpi_high[j] ) );
    nSel_dw = (double)chain_Kpi_afterCut_dw->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKpi_low[j],rangesKpi_high[j] ) );
    Eff=nSel_dw/nTot_dw;
    dEff=(1/nTot_dw) * TMath::Sqrt(nSel_dw * (1 - (nSel_dw/nTot_dw) ) );
    genLevelEff_Kpi_dw->SetBinContent(j+1,Eff);
    genLevelEff_Kpi_dw->SetBinError(j+1,dEff);

    nTot = nTot_up + nTot_dw;
    nSel = nSel_up + nSel_dw;
    Eff=nSel/nTot;
    dEff=(1/nTot) * TMath::Sqrt(nSel * (1 - (nSel/nTot) ) );
    genLevelEff_Kpi->SetBinContent(j+1,Eff);
    genLevelEff_Kpi->SetBinError(j+1,dEff);
  }

  for (int j=0;j<sizeof(binsKK)/sizeof(double)-1;++j) {
    
    nTot_up = (double)chain_KK_beforeCut_up->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKK_low[j],rangesKK_high[j] ) );
    nSel_up = (double)chain_KK_afterCut_up->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKK_low[j],rangesKK_high[j] ) );
    Eff=nSel_up/nTot_up;
    dEff=(1/nTot_up) * TMath::Sqrt(nSel_up * (1 - (nSel_up/nTot_up) ) );
    genLevelEff_KK_up->SetBinContent(j+1,Eff);
    genLevelEff_KK_up->SetBinError(j+1,dEff);

    nTot_dw = (double)chain_KK_beforeCut_dw->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKK_low[j],rangesKK_high[j] ) );
    nSel_dw = (double)chain_KK_afterCut_dw->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKK_low[j],rangesKK_high[j] ) );
    Eff=nSel_dw/nTot_dw;
    dEff=(1/nTot_dw) * TMath::Sqrt(nSel_dw * (1 - (nSel_dw/nTot_dw) ) );
    genLevelEff_KK_dw->SetBinContent(j+1,Eff);
    genLevelEff_KK_dw->SetBinError(j+1,dEff);

    nTot = nTot_up + nTot_dw;
    nSel = nSel_up + nSel_dw;
    Eff=nSel/nTot;
    dEff=(1/nTot) * TMath::Sqrt(nSel * (1 - (nSel/nTot) ) );
    genLevelEff_KK->SetBinContent(j+1,Eff);
    genLevelEff_KK->SetBinError(j+1,dEff);

    genLevelRelativeEff_KK->SetBinContent(j+1,Eff/genLevelEff_Kpi->GetBinContent(1));
    genLevelRelativeEff_KK->SetBinError(j+1,TMath::Sqrt( TMath::Power( (dEff/genLevelEff_Kpi->GetBinContent(1)),2) + 
							 TMath::Power( (Eff*genLevelEff_Kpi->GetBinError(1))/(genLevelEff_Kpi->GetBinContent(1)*genLevelEff_Kpi->GetBinContent(1)),2) ) );
    double dRelEff=TMath::Sqrt( TMath::Power( (dEff/genLevelEff_Kpi->GetBinContent(1)),2) +TMath::Power( (Eff*genLevelEff_Kpi->GetBinError(1))/(genLevelEff_Kpi->GetBinContent(1)*genLevelEff_Kpi->GetBinContent(1)),2) ) ;

  std::cout<<"Eff "<<Eff<<"+-"<<dEff<<" EffNorm= " << genLevelEff_Kpi->GetBinContent(1)<<"+-"<<genLevelEff_Kpi->GetBinError(1)<< std::endl;
std::cout<<"rel Eff "<<Eff/genLevelEff_Kpi->GetBinContent(1)<<"+-"<< dRelEff  << std::endl;
  }
  


  for (int j=0;j<sizeof(binspipi)/sizeof(double)-1;++j) {
    
    nTot_up = (double)chain_pipi_beforeCut_up->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangespipi_low[j],rangespipi_high[j] ) );
    nSel_up = (double)chain_pipi_afterCut_up->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangespipi_low[j],rangespipi_high[j] ) );
    Eff=nSel_up/nTot_up;
    dEff=(1/nTot_up) * TMath::Sqrt(nSel_up * (1 - (nSel_up/nTot_up) ) );
    genLevelEff_pipi_up->SetBinContent(j+1,Eff);
    genLevelEff_pipi_up->SetBinError(j+1,dEff);

    nTot_dw = (double)chain_pipi_beforeCut_dw->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangespipi_low[j],rangespipi_high[j] ) );
    nSel_dw = (double)chain_pipi_afterCut_dw->GetEntries(TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangespipi_low[j],rangespipi_high[j] ) );
    Eff=nSel_dw/nTot_dw;
    dEff=(1/nTot_dw) * TMath::Sqrt(nSel_dw * (1 - (nSel_dw/nTot_dw) ) );
    genLevelEff_pipi_dw->SetBinContent(j+1,Eff);
    genLevelEff_pipi_dw->SetBinError(j+1,dEff);

    nTot = nTot_up + nTot_dw;
    nSel = nSel_up + nSel_dw;
    Eff=nSel/nTot;
    dEff=(1/nTot) * TMath::Sqrt(nSel * (1 - (nSel/nTot) ) );
    genLevelEff_pipi->SetBinContent(j+1,Eff);
    genLevelEff_pipi->SetBinError(j+1,dEff);

    genLevelRelativeEff_pipi->SetBinContent(j+1,Eff/genLevelEff_Kpi->GetBinContent(1));
    genLevelRelativeEff_pipi->SetBinError(j+1,TMath::Sqrt( TMath::Power( (dEff/genLevelEff_Kpi->GetBinContent(1)),2) + 
							   TMath::Power( (Eff*genLevelEff_Kpi->GetBinError(1))/(genLevelEff_Kpi->GetBinContent(1)*genLevelEff_Kpi->GetBinContent(1)),2) ) );

    double dRelEff=TMath::Sqrt( TMath::Power( (dEff/genLevelEff_Kpi->GetBinContent(1)),2) +TMath::Power( (Eff*genLevelEff_Kpi->GetBinError(1))/(genLevelEff_Kpi->GetBinContent(1)*genLevelEff_Kpi->GetBinContent(1)),2) ) ;

  std::cout<<"Eff "<<Eff<<"+-"<<dEff<<" EffNorm= " << genLevelEff_Kpi->GetBinContent(1)<<"+-"<<genLevelEff_Kpi->GetBinError(1)<< std::endl;
  std::cout<<"rel Eff "<<Eff/genLevelEff_Kpi->GetBinContent(1)<<"+-"<< dRelEff  << std::endl;
    

  }
  

 
  TCanvas* a = new TCanvas("a","a");
  a->Divide(3,1);
  a->cd(1);
  
  genLevelEff_KK_up->GetYaxis()->SetRangeUser(0.,0.3);
  genLevelEff_KK_up->GetYaxis()->SetTitle("#epsilon_{genLevel}(D->KK#mu#mu)");
  genLevelEff_KK_up->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  genLevelEff_KK_up->Draw();
  genLevelEff_KK_dw->SetLineColor(kRed);
  genLevelEff_KK_dw->Draw("SAME");
  a->cd(2);
  genLevelEff_Kpi_up->GetYaxis()->SetRangeUser(0.,0.3);
  genLevelEff_Kpi_up->GetYaxis()->SetTitle("#epsilon_{genLevel}(D->K#pi#mu#mu)");
  genLevelEff_Kpi_up->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  genLevelEff_Kpi_up->Draw();
  genLevelEff_Kpi_dw->SetLineColor(kRed);
  genLevelEff_Kpi_dw->Draw("SAME");
  a->cd(3);
  genLevelEff_pipi_up->GetYaxis()->SetRangeUser(0.,0.3);
  genLevelEff_pipi_up->GetYaxis()->SetTitle("#epsilon_{genLevel}(D->#pi#pi#mu#mu)");
  genLevelEff_pipi_up->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  genLevelEff_pipi_up->Draw();
  genLevelEff_pipi_dw->SetLineColor(kRed);
  genLevelEff_pipi_dw->Draw("SAME");
  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/GenLevelCutEfficiency/generatorLevelCutEfficiency1.eps");

  TCanvas* b = new TCanvas("b","b");
  b->Divide(2,1);
  b->cd(1);
  genLevelRelativeEff_KK->GetYaxis()->SetRangeUser(0.7,1.2);
  genLevelRelativeEff_KK->GetYaxis()->SetTitle("#epsilon_{genLevel}(D->KK#mu#mu)/#epsilon_{genLevel}(D->K#pi#mu#mu)");
  genLevelRelativeEff_KK->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  genLevelRelativeEff_KK->Draw();
  b->cd(2);
  genLevelRelativeEff_pipi->GetYaxis()->SetRangeUser(0.7,1.2);
  genLevelRelativeEff_pipi->GetYaxis()->SetTitle("#epsilon_{genLevel}(D->#pi#pi#mu#mu)/#epsilon_{genLevel}(D->K#pi#mu#mu)");
  genLevelRelativeEff_pipi->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  genLevelRelativeEff_pipi->Draw();
  b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/GenLevelCutEfficiency/generatorLevelCutEfficiency2.eps");

  
  genLevelEff_pipi->Write();
  genLevelEff_KK->Write();
  genLevelEff_Kpi->Write();
  genLevelRelativeEff_KK->Write();
  genLevelRelativeEff_pipi->Write();

  TCanvas* niceC1 = new TCanvas("niceC1");
  genLevelEff_KK_up->Draw();
  genLevelEff_KK_dw->SetLineColor(kRed);
  genLevelEff_KK_dw->Draw("SAME");
  niceC1->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/GenLevelCutEfficiency/genLevelEffKK.eps");

  TCanvas* niceC2 = new TCanvas("niceC2");
  genLevelEff_Kpi_up->Draw();
  genLevelEff_Kpi_dw->SetLineColor(kRed);
  genLevelEff_Kpi_dw->Draw("SAME");
  niceC2->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/GenLevelCutEfficiency/genLevelEffKpi.eps");

  TCanvas* niceC3 = new TCanvas("niceC3");
  genLevelEff_pipi_up->Draw();
  genLevelEff_pipi_dw->SetLineColor(kRed);
  genLevelEff_pipi_dw->Draw("SAME");
  niceC3->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/GenLevelCutEfficiency/genLevelEffpipi.eps");

  TCanvas* niceC4 = new TCanvas("niceC4");
  genLevelRelativeEff_KK->Draw();
  niceC4->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/GenLevelCutEfficiency/genLevelRelativeEffKK.eps");

  TCanvas* niceC5 = new TCanvas("niceC5");
  genLevelRelativeEff_pipi->Draw();
  genLevelRelativeEff_pipi->GetYaxis()->SetTitleOffset(1.4);
  niceC5->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/GenLevelCutEfficiency/genLevelRelativeEffpipi.eps");

  fOut->Write();
  fOut->Close();

}

std::pair<double,double> fitMC(TString file,TString treeName, TString cut, TString namePlot){


  TFile* fileIn;
  fileIn= new TFile(file,"OPEN");
  TTree* tree = (TTree*) fileIn->Get(treeName);

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);
  tree->SetBranchStatus("D_TRUE_DiMuon_Mass",1);
  tree->SetBranchStatus("mu0_MuonNShared",1);
  tree->SetBranchStatus("mu1_MuonNShared",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("Dst_BKGCAT",1);
  tree->SetBranchStatus("totCandidates",1);
  tree->SetBranchStatus("isSelectedMultipleCandidate",1);
  tree->SetBranchStatus("isMatchedCandidate",1);


  RooRealVar Dst_DTF_D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1780., 1940.,"MeV");


  TFile* fOut = new TFile("temp.root","RECREATE");
  TTree* cutTree = tree->CopyTree(cut);
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(Dst_DTF_D0_M));
  TCanvas* c1= new TCanvas("test");


  RooRealVar D0_M_xi("D0_M_xi_","D0_M_xi",1.86560e+03,1860,1870);
  RooRealVar D0_M_lambda("D0_M_lambda_","D0_M_lambda",8.54839e+00,0.1,20);
  RooRealVar D0_M_gamma("D0_M_gamma_","D0_M_gamma",-5.42758e-02,-2,40);
  RooRealVar D0_M_delta("D0_M_delta_","D0_M_delta",6.09288e-01,0.,10);
  RooRealVar D0_M_chebyA("D0_M_chebyA_","D0_M_chebyA",0,-8,8);
  RooRealVar D0_M_chebyB("D0_M_chebyB_","D0_M_chebyB",-1.7004e-02,-1,1);
  RooRealVar D0_M_chebyC("D0_M_chebyC_","D0_M_chebyC",-1.7882e-02,-1,1);
  RooRealVar lambda("lambda_","lambda",0,-1,1);

  RooJohnsonSU Signal("Signal", "Signal D^{0} JSU", Dst_DTF_D0_M, D0_M_xi,D0_M_lambda,D0_M_gamma,D0_M_delta);
  // RooChebychev CombinatoricBkg("CombinatoricBkg", "Combinatoric Background (M)", Dst_DTF_D0_M, RooArgList(D0_M_chebyA));
  RooExponential CombinatoricBkg("CombinatoricBkg", "Combinatoric Background (M)", Dst_DTF_D0_M,lambda);

  RooWorkspace* w = new RooWorkspace(namePlot,kTRUE) ;
  w->import(Dst_DTF_D0_M);
  w->import(Signal);
  w->import(CombinatoricBkg);
 
  w->factory("SUM::sum(nSig[1000,0,10000]*Signal,nbkg[500,0,20000]*CombinatoricBkg)") ;

  // --- Perform extended ML fit of composite PDF to data ---                                                                                                                                 
  RooFitResult *result = w->pdf("sum")->fitTo(*data) ;
  //result->Print();                                                                                                                                                                          

  // --- Plot toy data and composite PDF overlaid ---                                                                                                                                         
  RooPlot* mesframe = w->var("Dst_DTF_D0_M")->frame() ;
  data->plotOn(mesframe) ;
  w->pdf("sum")->plotOn(mesframe) ;
  w->pdf("sum")->plotOn(mesframe,Components(*w->pdf("Signal")),LineStyle(kDashed),LineColor(kRed)) ;
  w->pdf("sum")->plotOn(mesframe,Components(*w->pdf("CombinatoricBkg")),LineStyle(kDashed),LineColor(kBlack)) ;
  w->pdf("sum")->paramOn(mesframe,Parameters(RooArgSet(*w->var("nSig"),*w->var("nbkg"))),Layout(.6,.93,.90));


  mesframe->Draw()  ;
  // c1->Print("../img/EfficiencyStudies/TriggerEffTagAndProbe/calibrationFits/"+namePlot+".eps");
  TString path = "../img/EfficiencyStudies/strippingEfficiency/withFit/";
  c1->Print(path+namePlot+".eps");

  fOut->Close();
  std::pair<double,double> signalYield =std::make_pair(w->var("nSig")->getVal(),w->var("nSig")->getError());

  delete fOut;
  delete result;
  delete tree;
  delete c1;
  w->Delete();

  return signalYield;

} 


std::pair<double,double> fitMC_alternativeModel(TString file,TString treeName, TString cut, TString namePlot){


  TFile* fileIn;
  fileIn= new TFile(file,"OPEN");
  TTree* tree = (TTree*) fileIn->Get(treeName);

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Dst_DTF_D0_M",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);
  tree->SetBranchStatus("D_TRUE_DiMuon_Mass",1);
  tree->SetBranchStatus("mu0_MuonNShared",1);
  tree->SetBranchStatus("mu1_MuonNShared",1);
  tree->SetBranchStatus("deltaM",1);
  tree->SetBranchStatus("Dst_BKGCAT",1);
  tree->SetBranchStatus("totCandidates",1);
  tree->SetBranchStatus("isSelectedMultipleCandidate",1);
  tree->SetBranchStatus("isMatchedCandidate",1);


  RooRealVar Dst_DTF_D0_M("Dst_DTF_D0_M", "m(h h #mu #mu)", 1780., 1940.,"MeV");


  TFile* fOut = new TFile("temp.root","RECREATE");
  TTree* cutTree = tree->CopyTree(cut);
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(Dst_DTF_D0_M));
  TCanvas* c1= new TCanvas("test");

  RooRealVar lambda("lambda_","lambda",0,-1,1);


  RooWorkspace* w = new RooWorkspace(namePlot,kTRUE) ;
  w->import(Dst_DTF_D0_M);
  w->factory("Gaussian::gauss1(Dst_DTF_D0_M,mean[1870,1760,1980],width[15,7,30])");
  w->factory("Gaussian::gauss2(Dst_DTF_D0_M,mean,width2[3,2,50])");
  w->factory("SUM::SigPDF(f[0.5,0,1]*gauss1,gauss2)") ;

  w->factory("Exponential::expBkg(Dst_DTF_D0_M,p[0,0,0])");
  w->factory("SUM::sum(nSig[10000,0,100000]*SigPDF,nbkg[5000,0,200000]*expBkg)") ;


  // --- Perform extended ML fit of composite PDF to data ---                                                                                                                                 
  RooFitResult *result = w->pdf("sum")->fitTo(*data) ;
  //result->Print();                                                                                                                                                                          

  // --- Plot toy data and composite PDF overlaid ---                                                                                                                                         
  RooPlot* mesframe = w->var("Dst_DTF_D0_M")->frame() ;
  data->plotOn(mesframe) ;
  w->pdf("sum")->plotOn(mesframe) ;
  w->pdf("sum")->plotOn(mesframe,Components(*w->pdf("SigPDF")),LineStyle(kDashed),LineColor(kGreen+3)) ;
  w->pdf("sum")->plotOn(mesframe,Components(*w->pdf("expBkg")),LineStyle(kDashed),LineColor(kViolet)) ;
  //RooArgSet set (nSig,nbkg);
  w->pdf("sum")->paramOn(mesframe,Parameters(RooArgSet(*w->var("nSig"),*w->var("nbkg"))),Layout(.6,.93,.90));

  mesframe->Draw()  ;
  // c1->Print("../img/EfficiencyStudies/TriggerEffTagAndProbe/calibrationFits/"+namePlot+".eps");
  TString path = "../img/EfficiencyStudies/strippingEfficiency/withFit/";
  c1->Print(path+namePlot+".eps");

  fOut->Close();
  std::pair<double,double> signalYield =std::make_pair(w->var("nSig")->getVal(),w->var("nSig")->getError());

  delete fOut;
  delete result;
  delete tree;
  delete c1;
  w->Delete();

  return signalYield;

} 



std::pair<double,double> getMCYield(TString fIn,TString nameTree, TString cut_sel) {

  TFile* file= new TFile(fIn,"OPEN");
  TTree* tuple = (TTree*) file->Get(nameTree);

  double nSel = (double)tuple->GetEntries(cut_sel);

  tuple->Delete();
  file->Close();
  file->Delete();

  return std::make_pair(nSel , TMath::Sqrt(nSel));

}


void drawRecoAndStrippingEfficiencyWithFit(bool useFit=true,TString additionalCut=""){

  //if useFit is set to 0, also efficiency with just counting is possible where a BKGCAT cut can be specified with additionalCut
  dcastyle();

  TString cutSel,cutNorm;
  std::pair<double,double> nNorm;
  std::pair<double,double> nSel;
  double Eff,dEff;
  TString data_file;
  TString qualityCut = "mu0_MuonNShared==0&&mu1_MuonNShared==0&&deltaM>144.5&&deltaM<146.5";
  if(additionalCut!="") qualityCut+="&&"+additionalCut;

  TString nameTarget="../img/EfficiencyStudies/strippingEfficiency/RecoAndStrippingEfficiency";
  if(useFit) nameTarget+="_fitted";
  if(additionalCut!="") nameTarget+="_"+additionalCut;


  TFile* fOut= new TFile(nameTarget+".root","RECREATE");

  TH1D * strippingEff_KK = new TH1D("strippingEff_KK","strippingEff_KK",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * strippingEff_pipi = new TH1D("strippingEff_pipi","strippingEff_pipi",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * strippingEff_Kpi = new TH1D("strippingEff_Kpi","strippingEff_Kpi",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D * strippingRelativeEff_KK = new TH1D("strippingRelativeEff_KK","strippingRelativeEff_KK",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * strippingRelativeEff_pipi = new TH1D("strippingRelativeEff_pipi","strippingRelativeEff_pipi",sizeof(binspipi)/sizeof(double)-1,binspipi);


  for(int j=0; j<rangesKpi_low.size();++j){

    data_file ="/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/MC12_DstD2Kpimumu_matched_andMultCand.root";
    // cutSel= qualityCut+"&&"+TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangesKpi_low[j],rangesKpi_high[j]);
    cutSel= qualityCut+"&&"+TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKpi_low[j],rangesKpi_high[j]);

    nNorm=getMCYield("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/MCTruthTuple_DstD2Kpimumu.root","MCTruthTuple",TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangesKpi_low[j],rangesKpi_high[j]));

    if(useFit) nSel=fitMC(data_file,"DecayTree",cutSel,"fitMC_Kpimumu_"+TString::Format("D_DiMuon_Mass_%.0f_D_DiMuon_Mass_%.0f",rangesKpi_low[j],rangesKpi_high[j]));
    else nSel=getMCYield(data_file,"DecayTree",cutSel);

    Eff=nSel.first/nNorm.first;
    dEff=1/nNorm.first* TMath::Sqrt(nSel.first* (1- (nSel.first/nNorm.first) ) );                                                                                                           

    //dEff= TMath::Sqrt((nSel.second/nNorm.first)*(nSel.second/nNorm.first) + (nSel.first*nNorm.second/(nNorm.first*nNorm.first))*(nSel.first*nNorm.second/(nNorm.first*nNorm.first)) );
    strippingEff_Kpi->SetBinContent(j+1,Eff);
    strippingEff_Kpi->SetBinError(j+1,dEff);
            
  }

  for(int j=0; j<rangesKK_low.size();++j){

    data_file ="/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/MC12_DstD2KKmumu_matched_andMultCand.root";
    //cutSel= qualityCut+"&&"+TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangesKK_low[j],rangesKK_high[j]);
    cutSel= qualityCut+"&&"+TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKK_low[j],rangesKK_high[j]);

    nNorm=getMCYield("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/MCTruthTuple_DstD2KKmumu.root","MCTruthTuple",TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangesKK_low[j],rangesKK_high[j]));
    if(useFit)nSel=fitMC(data_file,"DecayTree",cutSel,"fitMC_KKmumu_"+TString::Format("D_DiMuon_Mass_%.0f_D_DiMuon_Mass_%.0f",rangesKK_low[j],rangesKK_high[j]));
    else nSel=getMCYield(data_file,"DecayTree",cutSel);

    Eff=nSel.first/nNorm.first;
    dEff=1/nNorm.first* TMath::Sqrt(nSel.first* (1- (nSel.first/nNorm.first) ) );                                                                                                           
    cout<<"result "<<j<<"  "<<nNorm.first<<"  "<<nSel.first<<"  "<<Eff<<"  "<<dEff<<endl;

    //dEff= TMath::Sqrt((nSel.second/nNorm.first)*(nSel.second/nNorm.first) + (nSel.first*nNorm.second/(nNorm.first*nNorm.first))*(nSel.first*nNorm.second/(nNorm.first*nNorm.first)) );
    strippingEff_KK->SetBinContent(j+1,Eff);
    strippingEff_KK->SetBinError(j+1,dEff);

    strippingRelativeEff_KK->SetBinContent(j+1,Eff/strippingEff_Kpi->GetBinContent(1));
    strippingRelativeEff_KK->SetBinError(j+1,TMath::Sqrt( TMath::Power( (dEff/strippingEff_Kpi->GetBinContent(1)),2) +
							    TMath::Power( (Eff*strippingEff_Kpi->GetBinError(1))/(strippingEff_Kpi->GetBinContent(1)*strippingEff_Kpi->GetBinContent(1)),2) ));  
  
            
  }

  for(int j=0; j<rangespipi_low.size();++j){

    data_file ="/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MC12_DstD2pipimumu_matched_andMultCand.root";
    //cutSel= qualityCut+"&&"+TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangespipi_low[j],rangespipi_high[j]);
    cutSel= qualityCut+"&&"+TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangespipi_low[j],rangespipi_high[j]);
    
    nNorm=getMCYield("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MCTruthTuple_DstD2pipimumu.root","MCTruthTuple",TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangespipi_low[j],rangespipi_high[j]));
    if(useFit) nSel=fitMC(data_file,"DecayTree",cutSel,"fitMC_pipimumu_"+TString::Format("D_DiMuon_Mass_%.0f_D_DiMuon_Mass_%.0f",rangespipi_low[j],rangespipi_high[j]));
    else nSel=getMCYield(data_file,"DecayTree",cutSel);

    Eff=nSel.first/nNorm.first;
    dEff=1/nNorm.first* TMath::Sqrt(nSel.first* (1- (nSel.first/nNorm.first) ) );                                                                                                           

    //dEff= TMath::Sqrt((nSel.second/nNorm.first)*(nSel.second/nNorm.first) + (nSel.first*nNorm.second/(nNorm.first*nNorm.first))*(nSel.first*nNorm.second/(nNorm.first*nNorm.first)) );
    strippingEff_pipi->SetBinContent(j+1,Eff);
    strippingEff_pipi->SetBinError(j+1,dEff);
  
    strippingRelativeEff_pipi->SetBinContent(j+1,Eff/strippingEff_Kpi->GetBinContent(1));
    strippingRelativeEff_pipi->SetBinError(j+1,TMath::Sqrt( TMath::Power( (dEff/strippingEff_Kpi->GetBinContent(1)),2) +
							    TMath::Power( (Eff*strippingEff_Kpi->GetBinError(1))/(strippingEff_Kpi->GetBinContent(1)*strippingEff_Kpi->GetBinContent(1)),2) ));          
  }
  fOut->cd();
 
  strippingEff_Kpi->Write();
  strippingEff_KK->Write();
  strippingEff_pipi->Write();
  strippingRelativeEff_pipi->Write();
  strippingRelativeEff_pipi->Write();
  fOut->Write();


  }

void drawRecoAndStrippingEfficiencyWithFit_alternativeFitModel(bool useFit=true,TString additionalCut=""){ //for systematic studies

  //if useFit is set to 0, also efficiency with just counting is possible where a BKGCAT cut can be specified with additionalCut
  dcastyle();


  TString cutSel,cutNorm;
  std::pair<double,double> nNorm;
  std::pair<double,double> nSel;
  double Eff,dEff;
  TString data_file;
  TString qualityCut = "mu0_MuonNShared==0&&mu1_MuonNShared==0&&deltaM>144.5&&deltaM<146.5";
  if(additionalCut!="") qualityCut+="&&"+additionalCut;

  TString nameTarget="../img/EfficiencyStudies/strippingEfficiency/RecoAndStrippingEfficiency";
  if(useFit) nameTarget+="_fitted";
  if(additionalCut!="") nameTarget+="_"+additionalCut;


  TFile* fOut= new TFile(nameTarget+"_alternativeFitModel.root","RECREATE");

  TH1D * strippingEff_KK = new TH1D("strippingEff_KK","strippingEff_KK",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * strippingEff_pipi = new TH1D("strippingEff_pipi","strippingEff_pipi",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * strippingEff_Kpi = new TH1D("strippingEff_Kpi","strippingEff_Kpi",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D * strippingRelativeEff_KK = new TH1D("strippingRelativeEff_KK","strippingRelativeEff_KK",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * strippingRelativeEff_pipi = new TH1D("strippingRelativeEff_pipi","strippingRelativeEff_pipi",sizeof(binspipi)/sizeof(double)-1,binspipi);


  for(int j=0; j<rangesKpi_low.size();++j){

    data_file ="/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/MC12_DstD2Kpimumu_matched_andMultCand.root";
    // cutSel= qualityCut+"&&"+TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangesKpi_low[j],rangesKpi_high[j]);
    cutSel= qualityCut+"&&"+TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKpi_low[j],rangesKpi_high[j]);

    nNorm=getMCYield("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/MCTruthTuple_DstD2Kpimumu.root","MCTruthTuple",TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangesKpi_low[j],rangesKpi_high[j]));

    if(useFit) nSel=fitMC_alternativeModel(data_file,"DecayTree",cutSel,"fitMC_alternativeModel_Kpimumu_"+TString::Format("D_DiMuon_Mass_%.0f_D_DiMuon_Mass_%.0f",rangesKpi_low[j],rangesKpi_high[j]));
    else nSel=getMCYield(data_file,"DecayTree",cutSel);

    Eff=nSel.first/nNorm.first;
    dEff=1/nNorm.first* TMath::Sqrt(nSel.first* (1- (nSel.first/nNorm.first) ) );                                                                                                           

    //dEff= TMath::Sqrt((nSel.second/nNorm.first)*(nSel.second/nNorm.first) + (nSel.first*nNorm.second/(nNorm.first*nNorm.first))*(nSel.first*nNorm.second/(nNorm.first*nNorm.first)) );
    strippingEff_Kpi->SetBinContent(j+1,Eff);
    strippingEff_Kpi->SetBinError(j+1,dEff);
            
  }

  for(int j=0; j<rangesKK_low.size();++j){

    data_file ="/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/MC12_DstD2KKmumu_matched_andMultCand.root";
    //cutSel= qualityCut+"&&"+TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangesKK_low[j],rangesKK_high[j]);
    cutSel= qualityCut+"&&"+TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangesKK_low[j],rangesKK_high[j]);

    nNorm=getMCYield("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/MCTruthTuple_DstD2KKmumu.root","MCTruthTuple",TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangesKK_low[j],rangesKK_high[j]));
    if(useFit)nSel=fitMC_alternativeModel(data_file,"DecayTree",cutSel,"fitMC_alternativeModel_KKmumu_"+TString::Format("D_DiMuon_Mass_%.0f_D_DiMuon_Mass_%.0f",rangesKK_low[j],rangesKK_high[j]));
    else nSel=getMCYield(data_file,"DecayTree",cutSel);

    Eff=nSel.first/nNorm.first;
    dEff=1/nNorm.first* TMath::Sqrt(nSel.first* (1- (nSel.first/nNorm.first) ) );                                                                                                           
    cout<<"result "<<j<<"  "<<nNorm.first<<"  "<<nSel.first<<"  "<<Eff<<"  "<<dEff<<endl;

    //dEff= TMath::Sqrt((nSel.second/nNorm.first)*(nSel.second/nNorm.first) + (nSel.first*nNorm.second/(nNorm.first*nNorm.first))*(nSel.first*nNorm.second/(nNorm.first*nNorm.first)) );
    strippingEff_KK->SetBinContent(j+1,Eff);
    strippingEff_KK->SetBinError(j+1,dEff);

    strippingRelativeEff_KK->SetBinContent(j+1,Eff/strippingEff_Kpi->GetBinContent(1));
    strippingRelativeEff_KK->SetBinError(j+1,TMath::Sqrt( TMath::Power( (dEff/strippingEff_Kpi->GetBinContent(1)),2) +
							    TMath::Power( (Eff*strippingEff_Kpi->GetBinError(1))/(strippingEff_Kpi->GetBinContent(1)*strippingEff_Kpi->GetBinContent(1)),2) ));  
  
            
  }

  for(int j=0; j<rangespipi_low.size();++j){

    data_file ="/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MC12_DstD2pipimumu_matched_andMultCand.root";
    //cutSel= qualityCut+"&&"+TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangespipi_low[j],rangespipi_high[j]);
    cutSel= qualityCut+"&&"+TString::Format("D_DiMuon_Mass>%f&&D_DiMuon_Mass<%f",rangespipi_low[j],rangespipi_high[j]);
    
    nNorm=getMCYield("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/MCTruthTuple_DstD2pipimumu.root","MCTruthTuple",TString::Format("D_TRUE_DiMuon_Mass>%f&&D_TRUE_DiMuon_Mass<%f",rangespipi_low[j],rangespipi_high[j]));
    if(useFit) nSel=fitMC_alternativeModel(data_file,"DecayTree",cutSel,"fitMC_alternativeModel_pipimumu_"+TString::Format("D_DiMuon_Mass_%.0f_D_DiMuon_Mass_%.0f",rangespipi_low[j],rangespipi_high[j]));
    else nSel=getMCYield(data_file,"DecayTree",cutSel);

    Eff=nSel.first/nNorm.first;
    dEff=1/nNorm.first* TMath::Sqrt(nSel.first* (1- (nSel.first/nNorm.first) ) );                                                                                                           

    //dEff= TMath::Sqrt((nSel.second/nNorm.first)*(nSel.second/nNorm.first) + (nSel.first*nNorm.second/(nNorm.first*nNorm.first))*(nSel.first*nNorm.second/(nNorm.first*nNorm.first)) );
    strippingEff_pipi->SetBinContent(j+1,Eff);
    strippingEff_pipi->SetBinError(j+1,dEff);
  
    strippingRelativeEff_pipi->SetBinContent(j+1,Eff/strippingEff_Kpi->GetBinContent(1));
    strippingRelativeEff_pipi->SetBinError(j+1,TMath::Sqrt( TMath::Power( (dEff/strippingEff_Kpi->GetBinContent(1)),2) +
							    TMath::Power( (Eff*strippingEff_Kpi->GetBinError(1))/(strippingEff_Kpi->GetBinContent(1)*strippingEff_Kpi->GetBinContent(1)),2) ));          
  }
  fOut->cd();
 
  strippingEff_Kpi->Write();
  strippingEff_KK->Write();
  strippingEff_pipi->Write();
  strippingRelativeEff_pipi->Write();
  strippingRelativeEff_pipi->Write();
  fOut->Write();


  }




void compareDifferentRecoAndStrippingEfficiencies(){

  TString pathToFiles = "/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/strippingEfficiency/";

  TFile* fitted = new TFile(pathToFiles+"RecoAndStrippingEfficiency_fitted.root","OPEN");
  TFile* fitted_andSelectedMultCand = new TFile(pathToFiles+"RecoAndStrippingEfficiency_fitted_isSelectedMultipleCandidate.root","OPEN");
  TFile* bkgcat = new TFile(pathToFiles+"RecoAndStrippingEfficiency_Dst_BKGCAT<11.root","OPEN");
  TFile* bkgcat_andSelectedMultCand = new TFile(pathToFiles+"RecoAndStrippingEfficiency_Dst_BKGCAT<11&&isSelectedMultipleCandidate.root","OPEN");
  TFile* fitted_andMatched = new TFile(pathToFiles+"RecoAndStrippingEfficiency_fitted_isMatchedCandidate.root","OPEN");
  TFile* fitted_matchedAndSelectedMultCand = new TFile(pathToFiles+"RecoAndStrippingEfficiency_fitted_isMatchedCandidate&&isSelectedMultipleCandidate.root","OPEN");
  TFile* fitted_withGhosts = new TFile(pathToFiles+"RecoAndStrippingEfficiency_fitted_(Dst_BKGCAT<11||Dst_BKGCAT==60).root","OPEN");
  //TFile* fitted_withGhosts = new TFile(pathToFiles+"RecoAndStrippingEfficiency_fitted_withGhosts.root","OPEN");


  TH1D * fitted_KK_withGhosts= (TH1D*)fitted_withGhosts->Get("strippingRelativeEff_KK");
  TH1D * fitted_pipi_withGhosts= (TH1D*)fitted_withGhosts->Get("strippingRelativeEff_pipi");

  TH1D * fitted_KK= (TH1D*)fitted->Get("strippingRelativeEff_KK");
  TH1D * fitted_pipi= (TH1D*)fitted->Get("strippingRelativeEff_pipi");

  TH1D * fitted_andSelectedMultCand_KK= (TH1D*)fitted_andSelectedMultCand->Get("strippingRelativeEff_KK");
  TH1D * fitted_andSelectedMultCand_pipi= (TH1D*)fitted_andSelectedMultCand->Get("strippingRelativeEff_pipi");

  TH1D * bkgcat_KK= (TH1D*)bkgcat->Get("strippingRelativeEff_KK");
  TH1D * bkgcat_pipi= (TH1D*)bkgcat->Get("strippingRelativeEff_pipi");

  TH1D * bkgcat_andSelectedMultCand_KK= (TH1D*)bkgcat_andSelectedMultCand->Get("strippingRelativeEff_KK");
  TH1D * bkgcat_andSelectedMultCand_pipi= (TH1D*)bkgcat_andSelectedMultCand->Get("strippingRelativeEff_pipi");

  TH1D * fitted_andMatched_KK= (TH1D*)fitted_andMatched->Get("strippingRelativeEff_KK");
  TH1D * fitted_andMatched_pipi= (TH1D*)fitted_andMatched->Get("strippingRelativeEff_pipi");

  TH1D * fitted_matchedAndSelectedMultCand_KK= (TH1D*)fitted_matchedAndSelectedMultCand->Get("strippingRelativeEff_KK");
  TH1D * fitted_matchedAndSelectedMultCand_pipi= (TH1D*)fitted_matchedAndSelectedMultCand->Get("strippingRelativeEff_pipi");

  dcastyle();

  TCanvas* a = new TCanvas("a","a");
  a->Divide(1,2);
  a->cd(1);
  fitted_KK->GetYaxis()->SetRangeUser(0.5,0.8);
  fitted_KK->Draw();
  bkgcat_KK->SetLineColor(kRed);
  bkgcat_KK->Draw("SAME");
  fitted_KK_withGhosts->SetLineColor(kBlack);
  fitted_KK_withGhosts->Draw("SAME");

  a->cd(2);
  fitted_pipi->GetYaxis()->SetRangeUser(0.6,1.3);
  fitted_pipi->Draw();
  bkgcat_pipi->SetLineColor(kRed);
  bkgcat_pipi->Draw("SAME");
  fitted_pipi_withGhosts->SetLineColor(kBlack);
  fitted_pipi_withGhosts->Draw("SAME");

  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/strippingEfficiency/differentMethods1.eps");

  TCanvas* b = new TCanvas("b","b");
  b->Divide(1,2);
  b->cd(1);

  fitted_andSelectedMultCand_KK->GetYaxis()->SetRangeUser(0.5,0.8);
  fitted_andSelectedMultCand_KK->Draw();
  bkgcat_KK->SetLineColor(kRed);
  bkgcat_KK->Draw("SAME");
  bkgcat_andSelectedMultCand_KK->SetLineColor(kGreen);
  bkgcat_andSelectedMultCand_KK->Draw("SAME");
  fitted_andMatched_KK->SetLineColor(kOrange);
  fitted_andMatched_KK->Draw("SAME");
  fitted_matchedAndSelectedMultCand_KK->SetLineColor(kViolet);
  fitted_matchedAndSelectedMultCand_KK->Draw("SAME");
  fitted_KK_withGhosts->Draw("SAME");

  b->cd(2);
  fitted_andSelectedMultCand_pipi->GetYaxis()->SetRangeUser(0.7,1.25);
  fitted_andSelectedMultCand_pipi->Draw();
  bkgcat_pipi->SetLineColor(kRed);
  bkgcat_pipi->Draw("SAME");
  bkgcat_andSelectedMultCand_pipi->SetLineColor(kGreen);
  bkgcat_andSelectedMultCand_pipi->Draw("SAME");
  fitted_andMatched_pipi->SetLineColor(kOrange);
  fitted_andMatched_pipi->Draw("SAME");
  fitted_matchedAndSelectedMultCand_pipi->SetLineColor(kViolet);
  fitted_matchedAndSelectedMultCand_pipi->Draw("SAME");
  fitted_pipi_withGhosts->Draw("SAME");
  b->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/strippingEfficiency/differentMethods2.eps");

  TCanvas* c = new TCanvas("c","c");
  c->Divide(1,2);
  c->cd(1);

  bkgcat_KK->GetYaxis()->SetRangeUser(0.5,0.8);
  bkgcat_KK->SetLineColor(kRed);
  bkgcat_KK->Draw("SAME");
  fitted_andMatched_KK->SetLineColor(kOrange);
  fitted_andMatched_KK->Draw("SAME");
  fitted_KK_withGhosts->Draw("SAME");

  c->cd(2);
  bkgcat_pipi->GetYaxis()->SetRangeUser(0.7,1.25);
  bkgcat_pipi->SetLineColor(kRed);
  bkgcat_pipi->Draw("SAME");
  fitted_andMatched_pipi->SetLineColor(kOrange);
  fitted_andMatched_pipi->Draw("SAME");
  fitted_pipi_withGhosts->Draw("SAME");
  c->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/strippingEfficiency/differentMethods2.eps");

}




void drawRecoAndStrippingEfficiency(int cat){

  
  TChain *chain_KK_beforeCut_up = new TChain("MCTruthTuple/MCTruthTuple");
  TChain *chain_KK_afterCut_up = new TChain("MC12_DstD2KKMuMu/DecayTree");
  TChain *chain_KK_beforeCut_dw = new TChain("MCTruthTuple/MCTruthTuple");
  TChain *chain_KK_afterCut_dw = new TChain("MC12_DstD2KKMuMu/DecayTree");
  
  TChain *chain_pipi_beforeCut_up = new TChain("MCTruthTuple/MCTruthTuple");
  TChain *chain_pipi_afterCut_up = new TChain("MC12_DstD2pipiMuMu/DecayTree");
  TChain *chain_pipi_beforeCut_dw = new TChain("MCTruthTuple/MCTruthTuple");
  TChain *chain_pipi_afterCut_dw = new TChain("MC12_DstD2pipiMuMu/DecayTree");
  
  TChain *chain_Kpi_beforeCut_up = new TChain("MCTruthTuple/MCTruthTuple");
  TChain *chain_Kpi_afterCut_up = new TChain("MC12_DstD2KPiMuMu/DecayTree");
  TChain *chain_Kpi_beforeCut_dw = new TChain("MCTruthTuple/MCTruthTuple");
  TChain *chain_Kpi_afterCut_dw = new TChain("MC12_DstD2KPiMuMu/DecayTree");

 
  chain_Kpi_beforeCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/magUp/MC12_DstD2Kpimumu_magUp.root");
  chain_Kpi_afterCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/magUp/MC12_DstD2Kpimumu_magUp.root");
  chain_Kpi_beforeCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/magDw/MC12_DstD2Kpimumu_magDw.root");
  chain_Kpi_afterCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2Kpimumu/magDw/MC12_DstD2Kpimumu_magDw.root");

  chain_KK_beforeCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/magUp/MC12_DstD2KKmumu_magUp.root");
  chain_KK_afterCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/magUp/MC12_DstD2KKmumu_magUp.root");
  chain_KK_beforeCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/magDw/MC12_DstD2KKmumu_magDw.root");
  chain_KK_afterCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2KKmumu/magDw/MC12_DstD2KKmumu_magDw.root");

  
  chain_pipi_beforeCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/magUp/MC12_DstD2pipimumu_magUp.root");
  chain_pipi_afterCut_up->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/magUp/MC12_DstD2pipimumu_magUp.root");
  chain_pipi_beforeCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/magDw/MC12_DstD2pipimumu_magDw.root");
  chain_pipi_afterCut_dw->AddFile("/auto/data/mitzel/D2hhmumu/new/flaggedMC/D2pipimumu/magDw/MC12_DstD2pipimumu_magDw.root");
  

  double mu0_PX, mu0_PY, mu0_PZ;
  double mu1_PX, mu1_PY, mu1_PZ;
  int Dst_BKGCAT;
  
  chain_pipi_beforeCut_up->SetBranchAddress("muminus_TRUEP_X",&mu0_PX);  
  chain_pipi_beforeCut_up->SetBranchAddress("muminus_TRUEP_Y",&mu0_PY);  
  chain_pipi_beforeCut_up->SetBranchAddress("muminus_TRUEP_Z",&mu0_PZ);  
  chain_pipi_beforeCut_up->SetBranchAddress("muplus_TRUEP_X",&mu1_PX);  
  chain_pipi_beforeCut_up->SetBranchAddress("muplus_TRUEP_Y",&mu1_PY);  
  chain_pipi_beforeCut_up->SetBranchAddress("muplus_TRUEP_Z",&mu1_PZ);

  chain_pipi_beforeCut_dw->SetBranchAddress("muminus_TRUEP_X",&mu0_PX);  
  chain_pipi_beforeCut_dw->SetBranchAddress("muminus_TRUEP_Y",&mu0_PY);  
  chain_pipi_beforeCut_dw->SetBranchAddress("muminus_TRUEP_Z",&mu0_PZ);  
  chain_pipi_beforeCut_dw->SetBranchAddress("muplus_TRUEP_X",&mu1_PX);  
  chain_pipi_beforeCut_dw->SetBranchAddress("muplus_TRUEP_Y",&mu1_PY);  
  chain_pipi_beforeCut_dw->SetBranchAddress("muplus_TRUEP_Z",&mu1_PZ);
  
  chain_KK_beforeCut_up->SetBranchAddress("muminus_TRUEP_X",&mu0_PX);  
  chain_KK_beforeCut_up->SetBranchAddress("muminus_TRUEP_Y",&mu0_PY);  
  chain_KK_beforeCut_up->SetBranchAddress("muminus_TRUEP_Z",&mu0_PZ);  
  chain_KK_beforeCut_up->SetBranchAddress("muplus_TRUEP_X",&mu1_PX);  
  chain_KK_beforeCut_up->SetBranchAddress("muplus_TRUEP_Y",&mu1_PY);  
  chain_KK_beforeCut_up->SetBranchAddress("muplus_TRUEP_Z",&mu1_PZ);

  chain_KK_beforeCut_dw->SetBranchAddress("muminus_TRUEP_X",&mu0_PX);  
  chain_KK_beforeCut_dw->SetBranchAddress("muminus_TRUEP_Y",&mu0_PY);  
  chain_KK_beforeCut_dw->SetBranchAddress("muminus_TRUEP_Z",&mu0_PZ);  
  chain_KK_beforeCut_dw->SetBranchAddress("muplus_TRUEP_X",&mu1_PX);  
  chain_KK_beforeCut_dw->SetBranchAddress("muplus_TRUEP_Y",&mu1_PY);  
  chain_KK_beforeCut_dw->SetBranchAddress("muplus_TRUEP_Z",&mu1_PZ);

  chain_Kpi_beforeCut_up->SetBranchAddress("muminus_TRUEP_X",&mu0_PX);  
  chain_Kpi_beforeCut_up->SetBranchAddress("muminus_TRUEP_Y",&mu0_PY);  
  chain_Kpi_beforeCut_up->SetBranchAddress("muminus_TRUEP_Z",&mu0_PZ);  
  chain_Kpi_beforeCut_up->SetBranchAddress("muplus_TRUEP_X",&mu1_PX);  
  chain_Kpi_beforeCut_up->SetBranchAddress("muplus_TRUEP_Y",&mu1_PY);  
  chain_Kpi_beforeCut_up->SetBranchAddress("muplus_TRUEP_Z",&mu1_PZ);

  chain_Kpi_beforeCut_dw->SetBranchAddress("muminus_TRUEP_X",&mu0_PX);  
  chain_Kpi_beforeCut_dw->SetBranchAddress("muminus_TRUEP_Y",&mu0_PY);  
  chain_Kpi_beforeCut_dw->SetBranchAddress("muminus_TRUEP_Z",&mu0_PZ);  
  chain_Kpi_beforeCut_dw->SetBranchAddress("muplus_TRUEP_X",&mu1_PX);  
  chain_Kpi_beforeCut_dw->SetBranchAddress("muplus_TRUEP_Y",&mu1_PY);  
  chain_Kpi_beforeCut_dw->SetBranchAddress("muplus_TRUEP_Z",&mu1_PZ);
  
  chain_pipi_afterCut_up->SetBranchAddress("mu0_TRUEP_X",&mu0_PX);  
  chain_pipi_afterCut_up->SetBranchAddress("mu0_TRUEP_Y",&mu0_PY);  
  chain_pipi_afterCut_up->SetBranchAddress("mu0_TRUEP_Z",&mu0_PZ);  
  chain_pipi_afterCut_up->SetBranchAddress("mu1_TRUEP_X",&mu1_PX);  
  chain_pipi_afterCut_up->SetBranchAddress("mu1_TRUEP_Y",&mu1_PY);  
  chain_pipi_afterCut_up->SetBranchAddress("mu1_TRUEP_Z",&mu1_PZ);

  chain_pipi_afterCut_dw->SetBranchAddress("mu0_TRUEP_X",&mu0_PX);  
  chain_pipi_afterCut_dw->SetBranchAddress("mu0_TRUEP_Y",&mu0_PY);  
  chain_pipi_afterCut_dw->SetBranchAddress("mu0_TRUEP_Z",&mu0_PZ);  
  chain_pipi_afterCut_dw->SetBranchAddress("mu1_TRUEP_X",&mu1_PX);  
  chain_pipi_afterCut_dw->SetBranchAddress("mu1_TRUEP_Y",&mu1_PY);  
  chain_pipi_afterCut_dw->SetBranchAddress("mu1_TRUEP_Z",&mu1_PZ);
  
  chain_KK_afterCut_up->SetBranchAddress("mu0_TRUEP_X",&mu0_PX);  
  chain_KK_afterCut_up->SetBranchAddress("mu0_TRUEP_Y",&mu0_PY);  
  chain_KK_afterCut_up->SetBranchAddress("mu0_TRUEP_Z",&mu0_PZ);  
  chain_KK_afterCut_up->SetBranchAddress("mu1_TRUEP_X",&mu1_PX);  
  chain_KK_afterCut_up->SetBranchAddress("mu1_TRUEP_Y",&mu1_PY);  
  chain_KK_afterCut_up->SetBranchAddress("mu1_TRUEP_Z",&mu1_PZ);

  chain_KK_afterCut_dw->SetBranchAddress("mu0_TRUEP_X",&mu0_PX);  
  chain_KK_afterCut_dw->SetBranchAddress("mu0_TRUEP_Y",&mu0_PY);  
  chain_KK_afterCut_dw->SetBranchAddress("mu0_TRUEP_Z",&mu0_PZ);  
  chain_KK_afterCut_dw->SetBranchAddress("mu1_TRUEP_X",&mu1_PX);  
  chain_KK_afterCut_dw->SetBranchAddress("mu1_TRUEP_Y",&mu1_PY);  
  chain_KK_afterCut_dw->SetBranchAddress("mu1_TRUEP_Z",&mu1_PZ);

  chain_Kpi_afterCut_up->SetBranchAddress("mu0_TRUEP_X",&mu0_PX);  
  chain_Kpi_afterCut_up->SetBranchAddress("mu0_TRUEP_Y",&mu0_PY);  
  chain_Kpi_afterCut_up->SetBranchAddress("mu0_TRUEP_Z",&mu0_PZ);  
  chain_Kpi_afterCut_up->SetBranchAddress("mu1_TRUEP_X",&mu1_PX);  
  chain_Kpi_afterCut_up->SetBranchAddress("mu1_TRUEP_Y",&mu1_PY);  
  chain_Kpi_afterCut_up->SetBranchAddress("mu1_TRUEP_Z",&mu1_PZ);

  chain_Kpi_afterCut_dw->SetBranchAddress("mu0_TRUEP_X",&mu0_PX);  
  chain_Kpi_afterCut_dw->SetBranchAddress("mu0_TRUEP_Y",&mu0_PY);  
  chain_Kpi_afterCut_dw->SetBranchAddress("mu0_TRUEP_Z",&mu0_PZ);  
  chain_Kpi_afterCut_dw->SetBranchAddress("mu1_TRUEP_X",&mu1_PX);  
  chain_Kpi_afterCut_dw->SetBranchAddress("mu1_TRUEP_Y",&mu1_PY);  
  chain_Kpi_afterCut_dw->SetBranchAddress("mu1_TRUEP_Z",&mu1_PZ);

  chain_Kpi_afterCut_dw->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_Kpi_afterCut_up->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_KK_afterCut_dw->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_KK_afterCut_up->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_pipi_afterCut_dw->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);
  chain_pipi_afterCut_up->SetBranchAddress("Dst_BKGCAT",&Dst_BKGCAT);



  TFile* fOut= new TFile(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/strippingEfficiency_BKDCAT%i.root",cat),"RECREATE");

  TH1D * strippingEff_KK_up = new TH1D("strippingEff_KK_up","strippingEff_KK_up",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * strippingEff_KK_dw = new TH1D("strippingEff_KK_dw","strippingEff_KK_up",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * strippingEff_KK = new TH1D("strippingEff_KK","strippingEff_KK",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * nSel_KK_up = new TH1D("nSel_KK_up","nSel_KK_up",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * nTot_KK_up = new TH1D("nTot_KK_up","nTot_KK_up",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * nSel_KK_dw = new TH1D("nSel_KK_dw","nSel_KK_dw",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * nTot_KK_dw = new TH1D("nTot_KK_dw","nTot_KK_dw",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1D * strippingEff_pipi_up = new TH1D("strippingEff_pipi_up","strippingEff_pipi_up",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * strippingEff_pipi_dw = new TH1D("strippingEff_pipi_dw","strippingEff_pipi_up",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * strippingEff_pipi = new TH1D("strippingEff_pipi","strippingEff_pipi",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * nSel_pipi_up = new TH1D("nSel_pipi_up","nSel_pipi_up",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * nTot_pipi_up = new TH1D("nTot_pipi_up","nTot_pipi_up",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * nSel_pipi_dw = new TH1D("nSel_pipi_dw","nSel_pipi_dw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * nTot_pipi_dw = new TH1D("nTot_pipi_dw","nTot_pipi_dw",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1D * strippingEff_Kpi_up = new TH1D("strippingEff_Kpi_up","strippingEff_Kpi_up",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D * strippingEff_Kpi_dw = new TH1D("strippingEff_Kpi_dw","strippingEff_Kpi_up",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D * strippingEff_Kpi = new TH1D("strippingEff_Kpi","strippingEff_Kpi",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D * nSel_Kpi_up = new TH1D("nSel_Kpi_up","nSel_Kpi_up",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D * nTot_Kpi_up = new TH1D("nTot_Kpi_up","nTot_Kpi_up",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D * nSel_Kpi_dw = new TH1D("nSel_Kpi_dw","nSel_Kpi_dw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  TH1D * nTot_Kpi_dw = new TH1D("nTot_Kpi_dw","nTot_Kpi_dw",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D * strippingRelativeEff_KK = new TH1D("strippingRelativeEff_KK","strippingRelativeEff_KK",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * strippingRelativeEff_pipi = new TH1D("strippingRelativeEff_pipi","strippingRelativeEff_pipi",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TLorentzVector mu0, mu1;
  double D_DiMuon_Mass;

   std::cout<<"Filling histograms for pipimumu.."<<std::endl;

   // int cat=60;

  for(int i=0;i<chain_pipi_beforeCut_up->GetEntries();++i) {
    chain_pipi_beforeCut_up->GetEntry(i);
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nTot_pipi_up->Fill(D_DiMuon_Mass);
  }

  for(int i=0;i<chain_pipi_afterCut_up->GetEntries();++i) {
    chain_pipi_afterCut_up->GetEntry(i);
    if(Dst_BKGCAT>cat) continue;
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nSel_pipi_up->Fill(D_DiMuon_Mass);
  }

  for(int i=0;i<chain_pipi_beforeCut_dw->GetEntries();++i) {
    chain_pipi_beforeCut_dw->GetEntry(i);
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nTot_pipi_dw->Fill(D_DiMuon_Mass);
  }

  for(int i=0;i<chain_pipi_afterCut_dw->GetEntries();++i) {
    chain_pipi_afterCut_dw->GetEntry(i);
    if(Dst_BKGCAT>cat) continue;
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nSel_pipi_dw->Fill(D_DiMuon_Mass);
  }

  

   std::cout<<"Filling histograms for KKmumu.."<<std::endl;
 
  for(int i=0;i<chain_KK_beforeCut_up->GetEntries();++i) {
    chain_KK_beforeCut_up->GetEntry(i);
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nTot_KK_up->Fill(D_DiMuon_Mass);
  }

  for(int i=0;i<chain_KK_afterCut_up->GetEntries();++i) {
    chain_KK_afterCut_up->GetEntry(i);
    if(Dst_BKGCAT>cat) continue;
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nSel_KK_up->Fill(D_DiMuon_Mass);
  }

  for(int i=0;i<chain_KK_beforeCut_dw->GetEntries();++i) {
    chain_KK_beforeCut_dw->GetEntry(i);
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nTot_KK_dw->Fill(D_DiMuon_Mass);
  }

  for(int i=0;i<chain_KK_afterCut_dw->GetEntries();++i) {
    chain_KK_afterCut_dw->GetEntry(i);
    if(Dst_BKGCAT>cat) continue;
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nSel_KK_dw->Fill(D_DiMuon_Mass);
  }
  
   std::cout<<"Filling histograms for Kpimumu.."<<std::endl;

   std::cout<<"chain_Kpi_beforeCut_up->GetEntries() "<<chain_Kpi_beforeCut_up->GetEntries()<<std::endl;
   for(int i=0;i<chain_Kpi_beforeCut_up->GetEntries();++i) {
    chain_Kpi_beforeCut_up->GetEntry(i);
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nTot_Kpi_up->Fill(D_DiMuon_Mass);
  }
  std::cout<<"chain_Kpi_afterCut_up->GetEntries() "<<chain_Kpi_afterCut_up->GetEntries()<<std::endl;
  for(int i=0;i<chain_Kpi_afterCut_up->GetEntries();++i) {
    chain_Kpi_afterCut_up->GetEntry(i);
     if(Dst_BKGCAT>cat) continue;
     mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nSel_Kpi_up->Fill(D_DiMuon_Mass);
  }
  std::cout<<"chain_Kpi_beforeCut_dw->GetEntries() "<<chain_Kpi_beforeCut_dw->GetEntries()<<std::endl;
  for(int i=0;i<chain_Kpi_beforeCut_dw->GetEntries();++i) {
    chain_Kpi_beforeCut_dw->GetEntry(i);
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nTot_Kpi_dw->Fill(D_DiMuon_Mass);
  }
  std::cout<<"chain_Kpi_afterCut_dw->GetEntries() "<<chain_Kpi_afterCut_dw->GetEntries()<<std::endl;
  for(int i=0;i<chain_Kpi_afterCut_dw->GetEntries();++i) {
    chain_Kpi_afterCut_dw->GetEntry(i);
    if(Dst_BKGCAT>cat) continue;
    mu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
    mu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());
    D_DiMuon_Mass=(mu0+mu1).M();
    nSel_Kpi_dw->Fill(D_DiMuon_Mass);
  }

  double nTot_up,nTot_dw, nTot;
  double nSel_up,nSel_dw, nSel;
  double Eff, dEff;

  for (int j=0;j<sizeof(binsKpi)/sizeof(double)-1;++j) {

    std::cout<<"Computing efficiencies histograms for Kpimumu"<<std::endl;
    nTot_up = (double)nTot_Kpi_up->GetBinContent(j+1);
    nSel_up = (double)nSel_Kpi_up->GetBinContent(j+1);
    Eff=nSel_up/nTot_up;
    dEff=(1/nTot_up) * TMath::Sqrt(nSel_up * (1 - (nSel_up/nTot_up) ) );
    strippingEff_Kpi_up->SetBinContent(j+1,Eff);
    strippingEff_Kpi_up->SetBinError(j+1,dEff);

    nTot_dw = (double)nTot_Kpi_dw->GetBinContent(j+1);
    nSel_dw = (double)nSel_Kpi_dw->GetBinContent(j+1);
    Eff=nSel_dw/nTot_dw;
    dEff=(1/nTot_dw) * TMath::Sqrt(nSel_dw * (1 - (nSel_dw/nTot_dw) ) );
    strippingEff_Kpi_dw->SetBinContent(j+1,Eff);
    strippingEff_Kpi_dw->SetBinError(j+1,dEff);

    nTot = nTot_up + nTot_dw;
    nSel = nSel_up + nSel_dw;
    Eff=nSel/nTot;
    dEff=(1/nTot) * TMath::Sqrt(nSel * (1 - (nSel/nTot) ) );
    strippingEff_Kpi->SetBinContent(j+1,Eff);
    strippingEff_Kpi->SetBinError(j+1,dEff);
    std::cout<<"Kpi eff "<<j<<"  "<<Eff<<"+-"<<dEff<<std::endl; 
    std::cout<<"nTot "<<nTot<<" nSel "<<nSel<<" up "<<nTot_up<<","<<nSel_up<<" dw "<<nTot_dw<<","<<nSel_dw<<std::endl;  

 
  }
 

  for (int j=0;j<sizeof(binsKK)/sizeof(double)-1;++j) {

    std::cout<<"Computing efficiencies histograms for KKmumu"<<std::endl;

    nTot_up = (double)nTot_KK_up->GetBinContent(j+1);
    nSel_up = (double)nSel_KK_up->GetBinContent(j+1);
    Eff=nSel_up/nTot_up;
    dEff=(1/nTot_up) * TMath::Sqrt(nSel_up * (1 - (nSel_up/nTot_up) ) );
    strippingEff_KK_up->SetBinContent(j+1,Eff);
    strippingEff_KK_up->SetBinError(j+1,dEff);

    nTot_dw = (double)nTot_KK_dw->GetBinContent(j+1);
    nSel_dw = (double)nSel_KK_dw->GetBinContent(j+1);
    Eff=nSel_dw/nTot_dw;
    dEff=(1/nTot_dw) * TMath::Sqrt(nSel_dw * (1 - (nSel_dw/nTot_dw) ) );
    strippingEff_KK_dw->SetBinContent(j+1,Eff);
    strippingEff_KK_dw->SetBinError(j+1,dEff);


    nTot = nTot_up + nTot_dw;
    nSel = nSel_up + nSel_dw;
    Eff=nSel/nTot;
    dEff=(1/nTot) * TMath::Sqrt(nSel * (1 - (nSel/nTot) ) );
    strippingEff_KK->SetBinContent(j+1,Eff);
    strippingEff_KK->SetBinError(j+1,dEff);

    std::cout<<"nTot "<<nTot<<" nSel "<<nSel<<" up "<<nTot_up<<","<<nSel_up<<" dw "<<nTot_dw<<","<<nSel_dw<<std::endl;  
    std::cout<<"KK eff "<<j<<"  "<<Eff<<"+-"<<dEff<<std::endl; 

    strippingRelativeEff_KK->SetBinContent(j+1,Eff/strippingEff_Kpi->GetBinContent(1));
    strippingRelativeEff_KK->SetBinError(j+1,TMath::Sqrt( TMath::Power( (dEff/strippingEff_Kpi->GetBinContent(1)),2) + 
							  TMath::Power( (Eff*strippingEff_Kpi->GetBinError(1))/(strippingEff_Kpi->GetBinContent(1)*strippingEff_Kpi->GetBinContent(1)),2) ));

  }

  
  for (int j=0;j<sizeof(binspipi)/sizeof(double)-1;++j) {

    std::cout<<"Computing efficiencies histograms for pipimumu"<<std::endl;

    nTot_up = (double)nTot_pipi_up->GetBinContent(j+1);
    nSel_up = (double)nSel_pipi_up->GetBinContent(j+1);
    Eff=nSel_up/nTot_up;
    dEff=(1/nTot_up) * TMath::Sqrt(nSel_up * (1 - (nSel_up/nTot_up) ) );
    strippingEff_pipi_up->SetBinContent(j+1,Eff);
    strippingEff_pipi_up->SetBinError(j+1,dEff);

    nTot_dw = (double)nTot_pipi_dw->GetBinContent(j+1);
    nSel_dw = (double)nSel_pipi_dw->GetBinContent(j+1);
    Eff=nSel_dw/nTot_dw;
    dEff=(1/nTot_dw) * TMath::Sqrt(nSel_dw * (1 - (nSel_dw/nTot_dw) ) );
    strippingEff_pipi_dw->SetBinContent(j+1,Eff);
    strippingEff_pipi_dw->SetBinError(j+1,dEff);

    nTot = nTot_up + nTot_dw;
    nSel = nSel_up + nSel_dw;
    Eff=nSel/nTot;
    dEff=(1/nTot) * TMath::Sqrt(nSel * (1 - (nSel/nTot) ) );
    strippingEff_pipi->SetBinContent(j+1,Eff);
    strippingEff_pipi->SetBinError(j+1,dEff);

    strippingRelativeEff_pipi->SetBinContent(j+1,Eff/strippingEff_Kpi->GetBinContent(1));
    strippingRelativeEff_pipi->SetBinError(j+1,TMath::Sqrt( TMath::Power( (dEff/strippingEff_Kpi->GetBinContent(1)),2) + 
							  TMath::Power( (Eff*strippingEff_Kpi->GetBinError(1))/(strippingEff_Kpi->GetBinContent(1)*strippingEff_Kpi->GetBinContent(1)),2) ));

  }
  nTot_pipi_up->Write();  
  nTot_pipi_dw->Write();  
  nSel_pipi_up->Write();  
  nSel_pipi_dw->Write();  
  nTot_KK_up->Write();  
  nTot_KK_dw->Write();  
  nSel_KK_up->Write();  
  nSel_KK_dw->Write();  
  nTot_Kpi_up->Write();  
  nTot_Kpi_dw->Write();  
  nSel_Kpi_up->Write();  
  nSel_Kpi_dw->Write();  
 
  strippingEff_pipi_dw->Write();
  strippingEff_pipi_up->Write();
  strippingEff_KK_dw->Write();
  strippingEff_KK_up->Write();
  strippingEff_Kpi_dw->Write();
  strippingEff_Kpi_up->Write();

  strippingRelativeEff_pipi->Write();
  strippingRelativeEff_KK->Write();

  TCanvas* a = new TCanvas("a","a");
  a->Divide(1,3);
  a->cd(1);
  
  //strippingEff_KK_up->GetYaxis()->SetRangeUser(0.,0.3);
  strippingEff_KK_up->Draw();
  strippingEff_KK_dw->SetLineColor(kRed);
  strippingEff_KK_dw->Draw("SAME");
  a->cd(2);
  //strippingEff_Kpi_up->GetYaxis()->SetRangeUser(0.,0.3);
  strippingEff_Kpi_up->Draw();
  strippingEff_Kpi_dw->SetLineColor(kRed);
  strippingEff_Kpi_dw->Draw("SAME");
  a->cd(3);
  //strippingEff_pipi_up->GetYaxis()->SetRangeUser(0.,0.3);565.0_950
  strippingEff_pipi_up->Draw();
  strippingEff_pipi_dw->SetLineColor(kRed);
  strippingEff_pipi_dw->Draw("SAME");
  a->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/strippingEfficiency_BKDCAT%i_1.eps",cat));

  TCanvas* b = new TCanvas("b","b");
  b->Divide(1,2);
  b->cd(1);
  //strippingRelativeEff_KK->GetYaxis()->SetRangeUser(0.3,1.4);
  strippingRelativeEff_KK->Draw();
  b->cd(2);
  //strippingRelativeEff_pipi->GetYaxis()->SetRangeUser(0.3,1.4);
  strippingRelativeEff_pipi->Draw();
  b->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/strippingEfficiency_BKDCAT%i_2.eps",cat));
  fOut->Write();

}


std::pair<double,double> fitDsPhiPi(TString file,TString treeName, TString cut, TString namePlot){


  std::cout<<"Calibrating for "<<cut<<std::endl;
  dcastyle();

  TFile* fileIn;
  fileIn= new TFile(file,"OPEN");
  TTree* tree = (TTree*) fileIn->Get(treeName);
  
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("D_M",1);
  tree->SetBranchStatus("mu0_ProbNNmu",1);
  tree->SetBranchStatus("mu1_ProbNNmu",1);
  tree->SetBranchStatus("D_DiMuon_Mass",1);
  tree->SetBranchStatus("nTracks",1);
  tree->SetBranchStatus("mu0_MuonNShared",1);
  tree->SetBranchStatus("mu1_MuonNShared",1);
  tree->SetBranchStatus("mu0_L0MuonDecision_TOS",1);
  tree->SetBranchStatus("mu1_L0MuonDecision_TOS",1);

  tree->SetBranchStatus("mu0_PT",1);
  tree->SetBranchStatus("mu1_PT",1);

  tree->SetBranchStatus("mu0_ProbNNghost",1);
  tree->SetBranchStatus("mu1_ProbNNghost",1);
  tree->SetBranchStatus("h0_ProbNNghost",1);
  tree->SetBranchStatus("Slowpi_ProbNNghost",1);
  tree->SetBranchStatus("h0_ProbNNk",1);
 
 
  RooRealVar D_M("D_M", "m(#phi #pi)", 1920., 2050.,"MeV");


  TFile* fOut = new TFile("temp.root","RECREATE");
  TTree* cutTree = tree->CopyTree(cut);
  RooDataSet* data = new RooDataSet("data", "data", cutTree, RooArgSet(D_M));
  TCanvas* c1= new TCanvas("test");
  CreateSubPad(c1,0.25);

  RooWorkspace* w = new RooWorkspace(namePlot,kTRUE) ;
  w->import(D_M);
  w->factory("Gaussian::gauss1(D_M,mean[1970,1860,1980],width[15,7,20])");
  w->factory("Gaussian::gauss2(D_M,mean,width2[3,2,20])");
  w->factory("SUM::SigPDF(f[0.5,0,1]*gauss1,gauss2)") ;

  // w.factory("Gaussian::gauss2(D_M,mean2[1865,1860,1890],width2[10,5,15])");
  w->factory("Exponential::expBkg(D_M,p[0,-2,2])");
  w->factory("SUM::sum(nSig[10000,0,100000]*SigPDF,nbkg[5000,0,200000]*expBkg)") ;


  // --- Perform extended ML fit of composite PDF to data ---
  RooFitResult *result = w->pdf("sum")->fitTo(*data) ;
  //result->Print();

  c1->cd(1);
  // --- Plot toy data and composite PDF overlaid ---
  RooPlot* mesframe = w->var("D_M")->frame() ;
  data->plotOn(mesframe) ;
  w->pdf("sum")->plotOn(mesframe,Components(*w->pdf("SigPDF")),LineStyle(kDashed),LineColor(kRed)) ;
  w->pdf("sum")->plotOn(mesframe) ;
  
  mesframe->Draw()  ;

  c1->cd(2);
  RooHist* hpull = mesframe->pullHist() ;
  RooPlot* mesframe2 = w->var("D_M")->frame() ;
  mesframe2->addPlotable(hpull,"P") ;
  mesframe2->Draw();

  c1->Print("../img/EfficiencyStudies/TriggerEffTagAndProbe/calibrationFits/"+namePlot+".eps");

  fOut->Close();
  std::pair<double,double> signalYield =std::make_pair(w->var("nSig")->getVal(),w->var("nSig")->getError());

  delete fOut;
  delete result;
  delete tree;
  delete c1;
  w->Delete();
    
  return signalYield;

}


void calibrateTagAndProbeL0EfficiencyWithDsPhiPi(){

  TFile* fOut= new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/calibrationFile.root","RECREATE");
  TH1* AbsMCEff_Kpimumu = new TH1D("AbsMCEff_Kpimumu","Abs MC Eff_Kpimumu in bins of Pt",sizeof(binsPt)/sizeof(double)-1,binsPt);

  TString cutSel,cutNorm;
  std::pair<double,double> nNorm;
  std::pair<double,double> nSel;
  double Eff,dEff;
  TString data_file;
  TString qualityCut = "mu0_ProbNNmu>0.3&&mu0_L0MuonDecision_TOS&&mu0_MuonNShared==0&&mu1_MuonNShared==0";

  for(int j=0; j<rangesPt_low.size();++j){

    data_file ="/auto/data/mitzel/Ds2phipi/Ds2phipi.root";
    cutNorm= qualityCut+"&&"+TString::Format("mu1_PT>%f&&mu1_PT<%f",rangesPt_low[j],rangesPt_high[j]);
    cutSel = cutNorm+"&&"+"mu1_L0MuonDecision_TOS";  

    nNorm=fitDsPhiPi(data_file,"DecayTree",cutNorm,"normFit_"+TString::Format("mu1_PT_%.0f_mu1_PT_%.0f",rangesPt_low[j],rangesPt_high[j]));
    nSel=fitDsPhiPi(data_file,"DecayTree",cutSel,"selFit_"+TString::Format("mu1_PT_%.0f_mu1_PT_%.0f",rangesPt_low[j],rangesPt_high[j]));
       
    Eff=nSel.first/nNorm.first;
    dEff=1/nNorm.first* TMath::Sqrt(nSel.first* (1- (nSel.first/nNorm.first) ) );

    //dEff= TMath::Sqrt((nSel.second/nNorm.first)*(nSel.second/nNorm.first) + (nSel.first*nNorm.second/(nNorm.first*nNorm.first))*(nSel.first*nNorm.second/(nNorm.first*nNorm.first)) );

    AbsMCEff_Kpimumu->SetBinContent(j+1,Eff);
    AbsMCEff_Kpimumu->SetBinError(j+1,dEff);
  }
  fOut->cd();
  AbsMCEff_Kpimumu->Write();
  fOut->Write();

}


void calibrateTagAndProbeL0EfficiencyWithMC(){

  TFile* fOut= new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/calibrationMCFile.root","RECREATE");
  TH1* AbsMCEff_Kpimumu = new TH1D("AbsMCEff_Kpimumu","Abs MC Eff_Kpimumu in bins of Pt",sizeof(binsPt)/sizeof(double)-1,binsPt);

  TString MC_file, MC_file2;
  TString cut,cutNorm;
  std::pair<double,double> tempPair;

  for(int j=0; j<rangesPt_low.size();++j){

    MC_file ="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_675.0_875.0_D2pipimumuBDT_magUp.root";
    MC_file2 ="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/MCEfficiencyStudies/MC_D2Kpimumu_675.0_875.0_D2pipimumuBDT_magDw.root";
    //MC_file ="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_565.0_900.0_D2KKmumuBDT_magDw.root";
    //MC_file2 ="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_565.0_900.0_D2KKmumuBDT_magUp.root";
    //MC_file ="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_950.0_1100.0_D2pipimumuBDT_magUp.root";
    //MC_file2 ="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_950.0_1100.0_D2pipimumuBDT_magDw.root";
    cutNorm= TString::Format("mu0_L0MuonDecision_TOS&&mu1_PT>%f&&mu1_PT<%f",rangesPt_low[j],rangesPt_high[j]);
 
    cut = "mu1_L0MuonDecision_TOS";  
    tempPair = getMCEfficiencyFrom2Files(MC_file,MC_file2,"BDT_Tree",cut,cutNorm);
    AbsMCEff_Kpimumu->SetBinContent(j+1,tempPair.first);
    AbsMCEff_Kpimumu->SetBinError(j+1,tempPair.second);
  }
  fOut->cd();
  AbsMCEff_Kpimumu->Write();
  fOut->Write();

}

					    

void mergeDs2PhiPiTuples(){


  TChain* chain_up = new TChain("DstD2PhiPi/DecayTree");
  for(int i=1; i<400;++i) {
    chain_up->AddFile(TString::Format("/auto/data/mitzel/Ds2phipi/magUp/Ds2phipi_%i.root",i));
  }  
  TFile* fOut_up= new TFile("/auto/data/mitzel/Ds2phipi/magUp/Ds2phipi_magUp.root","RECREATE");
  TTree* cutTree_up=(TTree*)chain_up->CopyTree("abs(D_DiMuon_Mass-1019.5)<20 && D_M>(1968.5-60) && D_M<(1968.5+80)");
  cutTree_up->Write();  
  fOut_up->Write();
  fOut_up->Close();

  
  TChain* chain_dw = new TChain("DstD2PhiPi/DecayTree");
  for(int i=1; i<400;++i) {
    chain_dw->AddFile(TString::Format("/auto/data/mitzel/Ds2phipi/magDw/Ds2phipi_%i.root",i));
  }  
  TFile* fOut_dw= new TFile("/auto/data/mitzel/Ds2phipi/magDw/Ds2phipi_magDw.root","RECREATE");
  TTree* cutTree_dw=(TTree*)chain_dw->CopyTree("D_DiMuon_Mass>999.5&&D_DiMuon_Mass<1049.5");
  cutTree_dw->Write();  
  fOut_dw->Write();
  fOut_dw->Close();

  TFile* fOut_tot = new TFile("/auto/data/mitzel/Ds2phipi/Ds2phipi.root","RECREATE");
  TList* list = new TList();
  list->Add(cutTree_up);
  list->Add(cutTree_dw);
  
  TTree* mergedTree = TTree::MergeTrees(list);
  mergedTree->Write();
  fOut_tot->Close();
  


}



void applyTagAndProbeL0Efficiency(bool useData=true){
  
  //the calibration file 
  TString calibrationFile;
  if(useData) calibrationFile="/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/calibrationFile.root";
  else calibrationFile="/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/calibrationMCFile.root";

  TFile* fIn= new TFile(calibrationFile,"OPEN");
  TH1D* AbsMCEff_Kpimumu=(TH1D*)fIn->Get("AbsMCEff_Kpimumu");

  TString target;
  if(useData) target = "/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/L0TriggerEfficiencies";
  else target = "/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/L0TriggerEfficienciesMCClosureTest2";

  dcastyle();

  TFile* fOut = new TFile(target+".root","RECREATE");
  TH1D * L0Trigger_Eff_KKmumu = new TH1D("L0Trigger_Eff_KKmumu","L0Trigger_Eff_KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * L0Trigger_Eff_pipimumu = new TH1D("L0Trigger_Eff_pipimumu","L0Trigger_Eff_pipimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * L0Trigger_Eff_Kpimumu = new TH1D("L0Trigger_Eff_Kpimumu","L0Trigger_Eff_Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D * MC_L0Trigger_Eff_KKmumu = new TH1D("MC_L0Trigger_Eff_KKmumu","MC_L0Trigger_Eff_KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * MC_L0Trigger_Eff_pipimumu = new TH1D("MC_L0Trigger_Eff_pipimumu","MC_L0Trigger_Eff_pipimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * MC_L0Trigger_Eff_Kpimumu = new TH1D("MC_L0Trigger_Eff_Kpimumu","MC_L0Trigger_Eff_Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D * L0Trigger_RelativeEff_KKmumu = new TH1D("L0Trigger_RelativeEff_KKmumu","L0Trigger_RelativeEff_KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * L0Trigger_RelativeEff_pipimumu = new TH1D("L0Trigger_RelativeEff_pipimumu","L0Trigger_RelativeEff_pipimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);
 
  TH1D * MC_L0Trigger_RelativeEff_KKmumu = new TH1D("MC_L0Trigger_RelativeEff_KKmumu","MC_L0Trigger_RelativeEff_KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * MC_L0Trigger_RelativeEff_pipimumu = new TH1D("MC_L0Trigger_RelativeEff_pipimumu","MC_L0Trigger_RelativeEff_pipimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);
  


  double mu1_PT,mu0_PT,totEff;
  bool mu1_L0MuonDecision_TOS,mu0_L0MuonDecision_TOS;
  double BDT, mu0_ProbNNmu, mu1_ProbNNmu;
  int mu1_MuonNShared, mu0_MuonNShared;
  double Slowpi_ProbNNghost, h0_ProbNNghost, h1_ProbNNghost, mu0_ProbNNghost, mu1_ProbNNghost;
  double h0_ProbNNX, h1_ProbNNX;

  //KPimumu

  for(int i=0; i<rangesKpi_low.size();++i){
    
    TChain* chain;
    chain = new TChain("BDT_Tree");
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));

    chain->SetBranchAddress("mu1_PT",&mu1_PT);
    chain->SetBranchAddress("mu0_PT",&mu0_PT);
    chain->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
    chain->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
    chain->SetBranchAddress("BDT",&BDT);
    chain->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
    chain->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
    chain->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
    chain->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

    chain->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
    chain->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
    chain->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
    chain->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
    chain->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);

    chain->SetBranchAddress("h0_ProbNNk",&h0_ProbNNX);
    chain->SetBranchAddress("h1_ProbNNpi",&h1_ProbNNX);


    double tempEff1, tempEff0, tempEffProduct;
    double Eff_i,dEff_i,dTotEff;
    double dTempEff1, dTempEff0;
    int norm=0;

    for(int j=0; j<chain->GetEntries();++j){

      chain->GetEntry(j);
      //get trigger efficiency after PID
      //ALSO ADDHDRON PID HERE 
      
      if(mu0_ProbNNmu<0.5) continue;
      if(mu1_ProbNNmu<0.5) continue;
      //if(BDT<0) continue;
      if(mu1_MuonNShared>0)continue;
      if(mu0_MuonNShared>0)continue;
      
      //hadron PID
      if(mu0_ProbNNghost>0.5) continue;
      if(mu1_ProbNNghost>0.5) continue;
      if(h0_ProbNNghost>0.5) continue;
      if(h1_ProbNNghost>0.5) continue;
      if(Slowpi_ProbNNghost>0.5) continue;

      if(h0_ProbNNX<0.2) continue;
      if(h1_ProbNNX<0.2) continue;
      

      norm+=1;
      tempEff1=AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu1_PT));
      dTempEff1=AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu1_PT));
      tempEff0=AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu0_PT));
      dTempEff0=AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu0_PT));
      
      Eff_i=tempEff1+tempEff0-tempEff1*tempEff0;
      dEff_i = TMath::Sqrt( (1-tempEff1)*dTempEff0*(1-tempEff1)*dTempEff0 + (1-tempEff0)*dTempEff1*(1-tempEff0)*dTempEff1 );
      totEff+=Eff_i;
      dTotEff+=dEff_i*dEff_i;
      //cout<<j<<"  "<<dTempEff0<<"  "<<dTempEff1<<"  "<<dEff_i<<"  "<< dEff_i*dEff_i<<endl;
    }
    
    std::cout<<i<<"  "<<totEff<<"  "<< (double)norm <<" "<<dTotEff<<std::endl;
    totEff/=(double)norm;
    dTotEff=TMath::Sqrt(dTotEff)/(double)norm;
    std::cout<<i<<"  "<<totEff<<" +- "<<dTotEff<<std::endl;

    L0Trigger_Eff_Kpimumu->SetBinContent(i+1,totEff);
    L0Trigger_Eff_Kpimumu->SetBinError(i+1,dTotEff);

    TString cutHadronPID = "&&mu0_ProbNNghost<0.5&&mu1_ProbNNghost<0.5&&h0_ProbNNghost<0.5&&h1_ProbNNghost<0.5&&h0_ProbNNk>0.2&&h1_ProbNNpi>0.2&&Slowpi_ProbNNghost<0.5";
    std::pair<double,double> trueEff = getMCEfficiencyFrom2Files(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]),TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]),"BDT_Tree","(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)","mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&mu1_MuonNShared==0&&mu0_MuonNShared==0"+cutHadronPID);

    MC_L0Trigger_Eff_Kpimumu->SetBinContent(i+1,trueEff.first);
    MC_L0Trigger_Eff_Kpimumu->SetBinError(i+1,trueEff.second);

    chain->Clear();
  }

 
  for(int i=0; i<rangespipi_low.size();++i){

    TChain* chain;
    chain = new TChain("BDT_Tree");
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]));
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]));

    chain->SetBranchAddress("mu1_PT",&mu1_PT);
    chain->SetBranchAddress("mu0_PT",&mu0_PT);
    chain->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
    chain->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
    chain->SetBranchAddress("BDT",&BDT);
    chain->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
    chain->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
    chain->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
    chain->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

    chain->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
    chain->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
    chain->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
    chain->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
    chain->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);

    chain->SetBranchAddress("h0_ProbNNpi",&h0_ProbNNX);
    chain->SetBranchAddress("h1_ProbNNpi",&h1_ProbNNX);

    double tempEff1, tempEff0, tempEffProduct;
    double Eff_i,dEff_i,dTotEff;
    double dTempEff1, dTempEff0;
    int norm=0;

    for(int j=0; j<chain->GetEntries();++j){

      chain->GetEntry(j);
      //get trigger efficiency after PID
      //ALSO ADDHDRON PID HERE
      if(mu0_ProbNNmu<0.5) continue;
      if(mu1_ProbNNmu<0.5) continue;
      //if(BDT<0) continue;
      if(mu1_MuonNShared>0)continue;
      if(mu0_MuonNShared>0)continue;

      //hadron PID
      if(mu0_ProbNNghost>0.5) continue;
      if(mu1_ProbNNghost>0.5) continue;
      if(h0_ProbNNghost>0.5) continue;
      if(h1_ProbNNghost>0.5) continue;
      if(Slowpi_ProbNNghost>0.5) continue;

      if(h0_ProbNNX<0.2) continue;
      if(h1_ProbNNX<0.2) continue;
      
      norm+=1;
      tempEff1=AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu1_PT));
      dTempEff1=AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu1_PT));
      tempEff0=AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu0_PT));
      dTempEff0=AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu0_PT));
      
      Eff_i=tempEff1+tempEff0-tempEff1*tempEff0;
      dEff_i = TMath::Sqrt( (1-tempEff1)*dTempEff0*(1-tempEff1)*dTempEff0 + (1-tempEff0)*dTempEff1*(1-tempEff0)*dTempEff1 );
      totEff+=Eff_i;
      dTotEff+=dEff_i*dEff_i;

    }
    totEff/=(double)norm;
    dTotEff=TMath::Sqrt(dTotEff)/(double)norm;
    std::cout<<i<<"  "<<totEff<<" +- "<<dTotEff<<std::endl;

    L0Trigger_Eff_pipimumu->SetBinContent(i+1,totEff);
    L0Trigger_Eff_pipimumu->SetBinError(i+1,dTotEff);

    TString cutHadronPID = "&&mu0_ProbNNghost<0.5&&mu1_ProbNNghost<0.5&&h0_ProbNNghost<0.5&&h1_ProbNNghost<0.5&&h0_ProbNNpi>0.2&&h1_ProbNNpi>0.2&&Slowpi_ProbNNghost<0.5";
    std::pair<double,double> trueEff = getMCEfficiencyFrom2Files(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]),TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]),"BDT_Tree","(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)","mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&mu1_MuonNShared==0&&mu0_MuonNShared==0"+cutHadronPID);

    MC_L0Trigger_Eff_pipimumu->SetBinContent(i+1,trueEff.first);
    MC_L0Trigger_Eff_pipimumu->SetBinError(i+1,trueEff.second);
 
    //ratio

    MC_L0Trigger_RelativeEff_pipimumu->SetBinContent(i+1,trueEff.first/MC_L0Trigger_Eff_Kpimumu->GetBinContent(1));
    L0Trigger_RelativeEff_pipimumu->SetBinContent(i+1,totEff/L0Trigger_Eff_Kpimumu->GetBinContent(1));

    //error on ratio, assume uncorrelated numerator and denominator
    double Esig = totEff;
    double Enorm = L0Trigger_Eff_Kpimumu->GetBinContent(1);
    double dEsig= dTotEff;
    double dEnorm= L0Trigger_Eff_Kpimumu->GetBinError(1);
    
    double dR = TMath::Sqrt(TMath::Power( (dEsig/Enorm),2) + TMath::Power( (dEnorm*Esig/(Enorm*Enorm) ) ,2));
    L0Trigger_RelativeEff_pipimumu->SetBinError(i+1,dR);

    double MCEsig = trueEff.first;
    double MCEnorm = MC_L0Trigger_Eff_Kpimumu->GetBinContent(1);
    double MCdEsig= trueEff.second;
    double MCdEnorm= MC_L0Trigger_Eff_Kpimumu->GetBinError(1);

    double MCdR = TMath::Sqrt(TMath::Power( (MCdEsig/MCEnorm),2) + TMath::Power( (MCdEnorm*MCEsig/(MCEnorm*MCEnorm) ) ,2));
    MC_L0Trigger_RelativeEff_pipimumu->SetBinError(i+1,MCdR);


    chain->Clear();
  
    
  }

  for(int i=0; i<rangesKK_low.size();++i){
    
    TChain* chain;
    chain = new TChain("BDT_Tree");
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]));
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]));

    chain->SetBranchAddress("mu1_PT",&mu1_PT);
    chain->SetBranchAddress("mu0_PT",&mu0_PT);
    chain->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
    chain->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
    chain->SetBranchAddress("BDT",&BDT);
    chain->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
    chain->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
    chain->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
    chain->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

    chain->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
    chain->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
    chain->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
    chain->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
    chain->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);

    chain->SetBranchAddress("h0_ProbNNk",&h0_ProbNNX);
    chain->SetBranchAddress("h1_ProbNNk",&h1_ProbNNX);

    double tempEff1, tempEff0, tempEffProduct;
    double Eff_i,dEff_i,dTotEff;
    double dTempEff1, dTempEff0;
    int norm=0;

    for(int j=0; j<chain->GetEntries();++j){

      chain->GetEntry(j);
      //get trigger efficiency after PID
      //ALSO ADDHDRON PID HERE 
      if(mu0_ProbNNmu<0.5) continue;
      if(mu1_ProbNNmu<0.5) continue;
      //if(BDT<0) continue;
      if(mu1_MuonNShared>0)continue;
      if(mu0_MuonNShared>0)continue;

      //hadron PID
      if(mu0_ProbNNghost>0.5) continue;
      if(mu1_ProbNNghost>0.5) continue;
      if(h0_ProbNNghost>0.5) continue;
      if(h1_ProbNNghost>0.5) continue;
      if(Slowpi_ProbNNghost>0.5) continue;

      if(h0_ProbNNX<0.2) continue;
      if(h1_ProbNNX<0.2) continue;

      norm+=1;
      tempEff1=AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu1_PT));
      dTempEff1=AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu1_PT));
      tempEff0=AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu0_PT));
      dTempEff0=AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu0_PT));
      
      Eff_i=tempEff1+tempEff0-tempEff1*tempEff0;
      dEff_i = TMath::Sqrt( (1-tempEff1)*dTempEff0*(1-tempEff1)*dTempEff0 + (1-tempEff0)*dTempEff1*(1-tempEff0)*dTempEff1 );
      totEff+=Eff_i;
      dTotEff+=dEff_i*dEff_i;    

    }
    std::cout<<i<<"  "<<totEff<<" +- "<<dTotEff<<"  "<< norm <<std::endl;

    totEff/=(double)norm;
    dTotEff=TMath::Sqrt(dTotEff)/(double)norm;

    L0Trigger_Eff_KKmumu->SetBinContent(i+1,totEff);
    L0Trigger_Eff_KKmumu->SetBinError(i+1,dTotEff);

    std::cout<<i<<"  "<<totEff<<" +- "<<dTotEff<<"  "<< norm <<std::endl;

    TString cutHadronPID = "&&mu0_ProbNNghost<0.5&&mu1_ProbNNghost<0.5&&h0_ProbNNghost<0.5&&h1_ProbNNghost<0.5&&h0_ProbNNk>0.2&&h1_ProbNNk>0.2&&Slowpi_ProbNNghost<0.5";
 
    std::pair<double,double> trueEff = getMCEfficiencyFrom2Files(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]),TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]),"BDT_Tree","(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)","mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&mu1_MuonNShared==0&&mu0_MuonNShared==0"+cutHadronPID);

    MC_L0Trigger_Eff_KKmumu->SetBinContent(i+1,trueEff.first);
    MC_L0Trigger_Eff_KKmumu->SetBinError(i+1,trueEff.second);
 
    //ratio
    MC_L0Trigger_RelativeEff_KKmumu->SetBinContent(i+1,trueEff.first/MC_L0Trigger_Eff_Kpimumu->GetBinContent(1));
    L0Trigger_RelativeEff_KKmumu->SetBinContent(i+1,totEff/L0Trigger_Eff_Kpimumu->GetBinContent(1));

    //error on ratio, assume uncorrelated numerator and denominator
    double Esig = totEff;
    double Enorm = L0Trigger_Eff_Kpimumu->GetBinContent(1);
    double dEsig= dTotEff;
    double dEnorm= L0Trigger_Eff_Kpimumu->GetBinError(1);
    
    double dR = TMath::Sqrt(TMath::Power( (dEsig/Enorm),2) + TMath::Power( (dEnorm*Esig/(Enorm*Enorm) ) ,2));
    cout<<"KK "<<i<<" "<<Esig<<"  "<<Enorm<<"  "<<Esig/Enorm<<"  "<< dR<<endl;
    L0Trigger_RelativeEff_KKmumu->SetBinError(i+1,dR);

    double MCEsig = trueEff.first;
    double MCEnorm = MC_L0Trigger_Eff_Kpimumu->GetBinContent(1);
    double MCdEsig= trueEff.second;
    double MCdEnorm= MC_L0Trigger_Eff_Kpimumu->GetBinError(1);

    double MCdR = TMath::Sqrt(TMath::Power( (MCdEsig/MCEnorm),2) + TMath::Power( (MCdEnorm*MCEsig/(MCEnorm*MCEnorm) ) ,2));
    MC_L0Trigger_RelativeEff_KKmumu->SetBinError(i+1,MCdR);
    //cout<<"MC "<<i<<" "<<MCEsig<<"  "<<MCEnorm<<"  "<<MCEsig/MCEnorm<<" "<<MCdR<<endl;
    

    chain->Clear();
  }




  MC_L0Trigger_Eff_KKmumu->Write();
  MC_L0Trigger_Eff_pipimumu->Write();
  L0Trigger_Eff_KKmumu->Write();
  L0Trigger_Eff_pipimumu->Write();
  L0Trigger_RelativeEff_KKmumu->Write();
  L0Trigger_RelativeEff_pipimumu->Write();
  
  TCanvas* a = new TCanvas("a","a");
  a->Divide(4,3);
  a->cd(1);
  MC_L0Trigger_Eff_KKmumu->SetLineColor(kCyan);
  MC_L0Trigger_Eff_KKmumu->GetYaxis()->SetRangeUser(0,1);
  MC_L0Trigger_Eff_KKmumu->GetYaxis()->SetTitle("#epsilon_{L0}(D->KK#mu#mu)");
  MC_L0Trigger_Eff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  MC_L0Trigger_Eff_KKmumu->Draw();
  L0Trigger_Eff_KKmumu->Draw("SAME");
  a->cd(2);
  TH1* ratio1 = (TH1D*)MC_L0Trigger_Eff_KKmumu->Clone("DiffAbsEffKK");
  ratio1->Add(L0Trigger_Eff_KKmumu,-1); 
  ratio1->Draw();
  a->cd(3);
  MC_L0Trigger_RelativeEff_KKmumu->SetLineColor(kCyan);
  MC_L0Trigger_RelativeEff_KKmumu->GetYaxis()->SetTitle("#epsilon_{L0}(D->KK#mu#mu)/#epsilon_{L0}(D->KK#mu#mu)");
  MC_L0Trigger_RelativeEff_KKmumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  MC_L0Trigger_RelativeEff_KKmumu->GetYaxis()->SetRangeUser(0.7,.9);
  MC_L0Trigger_RelativeEff_KKmumu->Draw();
  L0Trigger_RelativeEff_KKmumu->Draw("SAME");
 
  a->cd(4);
  TH1* clone1 = (TH1D*)L0Trigger_RelativeEff_KKmumu->Clone("DiffRelativeEffKK");
  clone1->Add(MC_L0Trigger_RelativeEff_KKmumu,-1);
  clone1->Draw();
 
  a->cd(5);
  MC_L0Trigger_Eff_pipimumu->SetLineColor(kCyan);
  MC_L0Trigger_Eff_pipimumu->GetYaxis()->SetRangeUser(0,1);
  MC_L0Trigger_Eff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  MC_L0Trigger_Eff_pipimumu->GetYaxis()->SetTitle("#epsilon_{L0}(D->#pi#pi#mu#mu)");
  MC_L0Trigger_Eff_pipimumu->Draw();
  L0Trigger_Eff_pipimumu->Draw("SAME");
  a->cd(6);
  TH1* ratio2 = (TH1D*)MC_L0Trigger_Eff_pipimumu->Clone("DiffAbsEffpipi");
  ratio2->Add(L0Trigger_Eff_pipimumu,-1); 
  ratio2->Draw();

  a->cd(7);
  MC_L0Trigger_RelativeEff_pipimumu->SetLineColor(kCyan);
  MC_L0Trigger_RelativeEff_pipimumu->GetYaxis()->SetTitle("#epsilon_{L0}(D->#pi#pi#mu#mu)/#epsilon_{L0}(D->KK#mu#mu)");
  MC_L0Trigger_RelativeEff_pipimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  MC_L0Trigger_RelativeEff_pipimumu->GetYaxis()->SetRangeUser(0.8,1.4);
  MC_L0Trigger_RelativeEff_pipimumu->Draw();
  L0Trigger_RelativeEff_pipimumu->Draw("SAME");

  a->cd(8);
  TH1* clone2 = (TH1D*)L0Trigger_RelativeEff_pipimumu->Clone("DiffRelativeEffpipi");
  clone2->Add(MC_L0Trigger_RelativeEff_pipimumu,-1);
  clone2->Draw();

  a->cd(9);
  MC_L0Trigger_Eff_Kpimumu->SetLineColor(kCyan);
  MC_L0Trigger_Eff_Kpimumu->GetYaxis()->SetTitle("#epsilon_{L0}(D->K#pi#mu#mu)");
  MC_L0Trigger_Eff_Kpimumu->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  MC_L0Trigger_Eff_Kpimumu->GetYaxis()->SetRangeUser(0.,1.);
  MC_L0Trigger_Eff_Kpimumu->Draw();
  L0Trigger_Eff_Kpimumu->Draw("SAME");
  a->cd(10);
  TH1* ratio3 = (TH1D*)MC_L0Trigger_Eff_Kpimumu->Clone("DiffAbsEffKpi");
  ratio3->Add(L0Trigger_Eff_Kpimumu,-1); 
  ratio3->Draw();

  //a->Print(TString::Format("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/.eps",cat));
  a->Print(target+".eps");

  fOut->Write();

}


//finsih at some point...
/*
void applyTagAndProbeL0EfficiencyToyStudy(bool useData=true, int nToys=100){
  
  //the calibration file 
  TString calibrationFile;
  if(useData) calibrationFile="/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/calibrationFile.root";
  else calibrationFile="/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/calibrationMCFile.root";

  TFile* fIn= new TFile(calibrationFile,"OPEN");
  TH1D* AbsMCEff_Kpimumu=(TH1D*)fIn->Get("AbsMCEff_Kpimumu");

  TString target;
  if(useData) target = "/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/L0TriggerEfficienciesToyStudy";
  else target = "/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/L0TriggerEfficienciesMCClosureTestToyStudy";

  dcastyle();

  TFile* fOut = new TFile(target+".root","RECREATE");
  TH1D * L0Trigger_Eff_KKmumu = new TH1D("L0Trigger_Eff_KKmumu","L0Trigger_Eff_KKmumu",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1D * L0Trigger_Eff_pipimumu = new TH1D("L0Trigger_Eff_pipimumu","L0Trigger_Eff_pipimumu",sizeof(binspipi)/sizeof(double)-1,binspipi);
  TH1D * L0Trigger_Eff_Kpimumu = new TH1D("L0Trigger_Eff_Kpimumu","L0Trigger_Eff_Kpimumu",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  
  std::vector<std::pair<double,double>> Kpi_pt;
  std::vector<std::pair<double,double>> pipi_pt[5];

  std::vector<TH1*>  reso_pipi;
  TH1D* temp;
  for(int i=0; i<rangespipi_low.size();++i){
    temp = new TH1D(TString::Format("reso_pipi_bin_%i",i),TString::Format("reso_pipi_bin_%i",i),100,0.7,1.4);
    reso_pipi.push_back(temp);
  }

  double mu1_PT,mu0_PT,totEff;
  bool mu1_L0MuonDecision_TOS,mu0_L0MuonDecision_TOS;
  double BDT, mu0_ProbNNmu, mu1_ProbNNmu;
  int mu1_MuonNShared, mu0_MuonNShared;
  double Slowpi_ProbNNghost, h0_ProbNNghost, h1_ProbNNghost, mu0_ProbNNghost, mu1_ProbNNghost;
  double h0_ProbNNX, h1_ProbNNX;

  //KPimumu
  TRandom3* generator = new TRandom3(0);

  for(int i=0; i<rangesKpi_low.size();++i){
    
    TChain* chain;
    chain = new TChain("BDT_Tree");
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]));
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangesKpi_low[i],rangesKpi_high[i]));

    chain->SetBranchAddress("mu1_PT",&mu1_PT);
    chain->SetBranchAddress("mu0_PT",&mu0_PT);
    chain->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
    chain->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
    chain->SetBranchAddress("BDT",&BDT);
    chain->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
    chain->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
    chain->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
    chain->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

    chain->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
    chain->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
    chain->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
    chain->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
    chain->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);

    chain->SetBranchAddress("h0_ProbNNk",&h0_ProbNNX);
    chain->SetBranchAddress("h1_ProbNNpi",&h1_ProbNNX);


    double tempEff1, tempEff0, tempEffProduct;
    double Eff_i,dEff_i,dTotEff;
    double dTempEff1, dTempEff0;
    int norm=0;

    for(int j=0; j<chain->GetEntries();++j){

      chain->GetEntry(j);
      //get trigger efficiency after PID
      //ALSO ADDHDRON PID HERE 
      
      if(mu0_ProbNNmu<0.5) continue;
      if(mu1_ProbNNmu<0.5) continue;
      //if(BDT<0) continue;
      if(mu1_MuonNShared>0)continue;
      if(mu0_MuonNShared>0)continue;
      
      //hadron PID
      if(mu0_ProbNNghost>0.5) continue;
      if(mu1_ProbNNghost>0.5) continue;
      if(h0_ProbNNghost>0.5) continue;
      if(h1_ProbNNghost>0.5) continue;
      if(Slowpi_ProbNNghost>0.5) continue;

      if(h0_ProbNNX<0.2) continue;
      if(h1_ProbNNX<0.2) continue;
      
      norm+=1;
      Kpi_pt.push_back(std::make_pair(mu0_PT,mu1_PT));
              
      //cout<<j<<"  "<<dTempEff0<<"  "<<dTempEff1<<"  "<<dEff_i<<"  "<< dEff_i*dEff_i<<endl;
    }
    
    chain->Clear();
  }
    
 
  for(int i=0; i<rangespipi_low.size();++i){

    TChain* chain;

    std::vector<std::pair<double,double>> temp;

    chain = new TChain("BDT_Tree");
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]));
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magDw.root",rangespipi_low[i],rangespipi_high[i]));

    chain->SetBranchAddress("mu1_PT",&mu1_PT);
    chain->SetBranchAddress("mu0_PT",&mu0_PT);
    chain->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
    chain->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
    chain->SetBranchAddress("BDT",&BDT);
    chain->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
    chain->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
    chain->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
    chain->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

    chain->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
    chain->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
    chain->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
    chain->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
    chain->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);

    chain->SetBranchAddress("h0_ProbNNpi",&h0_ProbNNX);
    chain->SetBranchAddress("h1_ProbNNpi",&h1_ProbNNX);

    double tempEff1, tempEff0, tempEffProduct;
    double Eff_i,dEff_i,dTotEff;
    double dTempEff1, dTempEff0;
    int norm=0;

    for(int j=0; j<chain->GetEntries();++j){

      chain->GetEntry(j);
      //get trigger efficiency after PID
      //ALSO ADDHDRON PID HERE
      if(mu0_ProbNNmu<0.5) continue;
      if(mu1_ProbNNmu<0.5) continue;
      //if(BDT<0) continue;
      if(mu1_MuonNShared>0)continue;
      if(mu0_MuonNShared>0)continue;

      //hadron PID
      if(mu0_ProbNNghost>0.5) continue;
      if(mu1_ProbNNghost>0.5) continue;
      if(h0_ProbNNghost>0.5) continue;
      if(h1_ProbNNghost>0.5) continue;
      if(Slowpi_ProbNNghost>0.5) continue;

      if(h0_ProbNNX<0.2) continue;
      if(h1_ProbNNX<0.2) continue;
      
      norm+=1;

      temp.push_back(std::make_pair(mu0_PT,mu1_PT));
		  
    }
    pipi_pt[i]=temp;
    chain->Clear();
  
  } //end loop over bins

  for(int k=0;k<nToys;++k){

    for(int l=0;l<)
    tempEff1=generator->Gaus(AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu1_PT)),AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu1_PT)));
    tempEff0=generator->Gaus(AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu0_PT)),AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu0_PT)))
      
  /*
  
    reso_pipi[i]->Fill(totEff/L0Trigger_Eff_Kpimumu->GetBinContent(1));



      tempEff1=generator->Gaus(AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu1_PT)),AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu1_PT)));
      tempEff0=generator->Gaus(AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu0_PT)),AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu0_PT)));
            
      Eff_i=tempEff1+tempEff0-tempEff1*tempEff0;
      totEff+=Eff_i;


      tempEff1=generator->Gaus(AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu1_PT)),AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu1_PT)));
      tempEff0=generator->Gaus(AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu0_PT)),AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu0_PT)));
      Eff_i=tempEff1+tempEff0-tempEff1*tempEff0;
      totEff+=Eff_i;
	
      totEff/=(double)norm;
    
  *////////

  /*
  for(int i=0; i<rangesKK_low.size();++i){
    
    TChain* chain;
    chain = new TChain("BDT_Tree");
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]));
    chain->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]));

    chain->SetBranchAddress("mu1_PT",&mu1_PT);
    chain->SetBranchAddress("mu0_PT",&mu0_PT);
    chain->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
    chain->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
    chain->SetBranchAddress("BDT",&BDT);
    chain->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
    chain->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
    chain->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
    chain->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

    chain->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
    chain->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
    chain->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
    chain->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
    chain->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);

    chain->SetBranchAddress("h0_ProbNNk",&h0_ProbNNX);
    chain->SetBranchAddress("h1_ProbNNk",&h1_ProbNNX);

    double tempEff1, tempEff0, tempEffProduct;
    double Eff_i,dEff_i,dTotEff;
    double dTempEff1, dTempEff0;
    int norm=0;

    for(int j=0; j<chain->GetEntries();++j){

      chain->GetEntry(j);
      //get trigger efficiency after PID
      //ALSO ADDHDRON PID HERE 
      if(mu0_ProbNNmu<0.5) continue;
      if(mu1_ProbNNmu<0.5) continue;
      //if(BDT<0) continue;
      if(mu1_MuonNShared>0)continue;
      if(mu0_MuonNShared>0)continue;

      //hadron PID
      if(mu0_ProbNNghost>0.5) continue;
      if(mu1_ProbNNghost>0.5) continue;
      if(h0_ProbNNghost>0.5) continue;
      if(h1_ProbNNghost>0.5) continue;
      if(Slowpi_ProbNNghost>0.5) continue;

      if(h0_ProbNNX<0.2) continue;
      if(h1_ProbNNX<0.2) continue;

      norm+=1;
      tempEff1=AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu1_PT));
      dTempEff1=AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu1_PT));
      tempEff0=AbsMCEff_Kpimumu->GetBinContent(AbsMCEff_Kpimumu->FindBin(mu0_PT));
      dTempEff0=AbsMCEff_Kpimumu->GetBinError(AbsMCEff_Kpimumu->FindBin(mu0_PT));
      
      Eff_i=tempEff1+tempEff0-tempEff1*tempEff0;
      dEff_i = TMath::Sqrt( (1-tempEff1)*dTempEff0*(1-tempEff1)*dTempEff0 + (1-tempEff0)*dTempEff1*(1-tempEff0)*dTempEff1 );
      totEff+=Eff_i;
      dTotEff+=dEff_i*dEff_i;    

    }
    std::cout<<i<<"  "<<totEff<<" +- "<<dTotEff<<"  "<< norm <<std::endl;

    totEff/=(double)norm;
    dTotEff=TMath::Sqrt(dTotEff)/(double)norm;

    L0Trigger_Eff_KKmumu->SetBinContent(i+1,totEff);
    L0Trigger_Eff_KKmumu->SetBinError(i+1,dTotEff);

    std::cout<<i<<"  "<<totEff<<" +- "<<dTotEff<<"  "<< norm <<std::endl;

    TString cutHadronPID = "&&mu0_ProbNNghost<0.5&&mu1_ProbNNghost<0.5&&h0_ProbNNghost<0.5&&h1_ProbNNghost<0.5&&h0_ProbNNk>0.2&&h1_ProbNNk>0.2&&Slowpi_ProbNNghost<0.5";
 
    std::pair<double,double> trueEff = getMCEfficiencyFrom2Files(TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]),TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magDw.root",rangesKK_low[i],rangesKK_high[i]),"BDT_Tree","(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)","mu0_ProbNNmu>0.5&&mu1_ProbNNmu>0.5&&mu1_MuonNShared==0&&mu0_MuonNShared==0"+cutHadronPID);

    MC_L0Trigger_Eff_KKmumu->SetBinContent(i+1,trueEff.first);
    MC_L0Trigger_Eff_KKmumu->SetBinError(i+1,trueEff.second);
 
    //ratio
    MC_L0Trigger_RelativeEff_KKmumu->SetBinContent(i+1,trueEff.first/MC_L0Trigger_Eff_Kpimumu->GetBinContent(1));
    L0Trigger_RelativeEff_KKmumu->SetBinContent(i+1,totEff/L0Trigger_Eff_Kpimumu->GetBinContent(1));

    //error on ratio, assume uncorrelated numerator and denominator
    double Esig = totEff;
    double Enorm = L0Trigger_Eff_Kpimumu->GetBinContent(1);
    double dEsig= dTotEff;
    double dEnorm= L0Trigger_Eff_Kpimumu->GetBinError(1);
    
    double dR = TMath::Sqrt(TMath::Power( (dEsig/Enorm),2) + TMath::Power( (dEnorm*Esig/(Enorm*Enorm) ) ,2));
    cout<<"KK "<<i<<" "<<Esig<<"  "<<Enorm<<"  "<<Esig/Enorm<<"  "<< dR<<endl;
    L0Trigger_RelativeEff_KKmumu->SetBinError(i+1,dR);

    double MCEsig = trueEff.first;
    double MCEnorm = MC_L0Trigger_Eff_Kpimumu->GetBinContent(1);
    double MCdEsig= trueEff.second;
    double MCdEnorm= MC_L0Trigger_Eff_Kpimumu->GetBinError(1);

    double MCdR = TMath::Sqrt(TMath::Power( (MCdEsig/MCEnorm),2) + TMath::Power( (MCdEnorm*MCEsig/(MCEnorm*MCEnorm) ) ,2));
    MC_L0Trigger_RelativeEff_KKmumu->SetBinError(i+1,MCdR);
    //cout<<"MC "<<i<<" "<<MCEsig<<"  "<<MCEnorm<<"  "<<MCEsig/MCEnorm<<" "<<MCdR<<endl;
    

    chain->Clear();
    }

      for(int i=0; i<rangespipi_low.size();++i){
    reso_pipi[i]->Write();
  }

  fOut->Write();

}

*/

void createDataTuplesForEffStudies(){
  
  TChain* Tree_D2Kpimumu = new TChain("DstD2KPiMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2Kpimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));
    Tree_D2Kpimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }

  D2KpimumuReader* Kpi_Reader = new D2KpimumuReader(Tree_D2Kpimumu);
  Kpi_Reader->createSubsampleNoTriggerCuts("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/noTriggerCuts/D2Kpimumu_noTriggerCuts.root");


  TChain* Tree_D2pipimumu = new TChain("DstD2PiPiMuMu/DecayTree");

  for (int i=0; i<1600; ++i) {
    Tree_D2pipimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magUp/2012Data_D2hhmumu_%i.root",i));
    Tree_D2pipimumu->AddFile(TString::Format("/auto/data/mitzel/D2hhmumu/new/2012Data/magDw/2012Data_D2hhmumu_%i.root",i));
  }

  D2pipimumuReader* pipi_Reader = new D2pipimumuReader(Tree_D2Kpimumu);
  pipi_Reader->createSubsampleNoTriggerCuts("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/noTriggerCuts/D2pipimumu_noTriggerCuts.root");

}

void MC_HltTrigger_efficiency(double ghostProbCut,double hadronPID,double muonPID,double BDT,bool applyPID, bool applyBDT, TString triggerLevel){

  

  TString cutHlt1 = "(mu0_Hlt1TrackMuonDecision_TOS || mu1_Hlt1TrackMuonDecision_TOS || D_Hlt1TrackAllL0Decision_TOS )";
  //TString cutHlt1 = "D_Hlt1TrackAllL0Decision_TOS";

  TString cutHlt2_Kpi =  "D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS";
  TString cutHlt2_KK =  "D_Hlt2CharmSemilepD02KKMuMuDecision_TOS";
  TString cutHlt2_pipi =  "D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS";

  TString selectionCut_Trigger_Kpi, selectionCut_Trigger_KK,selectionCut_Trigger_pipi;
  TString normalizationCut_Kpi,normalizationCut_KK,normalizationCut_pipi;

  TString nameTarget =  "/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/HltMCEfficiencies/MCHltTrigger_Efficiencies";

  normalizationCut_Kpi="(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)";
  normalizationCut_KK="(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)";
  normalizationCut_pipi="(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)";

  
  if(triggerLevel=="Hlt1") {
    selectionCut_Trigger_Kpi = cutHlt1 ;
    selectionCut_Trigger_KK = cutHlt1 ;
    selectionCut_Trigger_pipi = cutHlt1 ;
    nameTarget+="_Hlt1" ;
 }
  if(triggerLevel=="Hlt2") {
    selectionCut_Trigger_Kpi = cutHlt2_Kpi;
    selectionCut_Trigger_KK = cutHlt2_KK;
    selectionCut_Trigger_pipi = cutHlt2_pipi;
    normalizationCut_Kpi+="&&"+cutHlt1;
    normalizationCut_KK+="&&"+cutHlt1;
    normalizationCut_pipi+="&&"+cutHlt1;
    nameTarget+="_Hlt2";

  }
  
  if(triggerLevel=="FullHlt") { 
    selectionCut_Trigger_Kpi = cutHlt1+"&&"+cutHlt2_Kpi;
    selectionCut_Trigger_KK = cutHlt1+"&&"+cutHlt2_KK;
    selectionCut_Trigger_pipi = cutHlt1+"&&"+cutHlt2_pipi;
    nameTarget+="_Hlt1_and_Hlt2";
  } 
  
  //form PID strings
  TString cutNshared="mu0_MuonNShared==0&&mu1_MuonNShared==0&&deltaM>144.5&&deltaM<146.5"; //preselction cuts
  TString cut_hadronPID_Kpimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_pipimumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNpi>%.1f&&h1_ProbNNpi>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_hadronPID_KKmumu=TString::Format("h0_ProbNNghost<%.1f&&h1_ProbNNghost<%.1f&&h0_ProbNNk>%.1f&&h1_ProbNNk>%.1f",ghostProbCut,ghostProbCut,hadronPID,hadronPID); 
  TString cut_SlowpiPID = TString::Format("Slowpi_ProbNNghost<%.1f",ghostProbCut); 
  TString cut_MuonPID=TString::Format("mu0_ProbNNghost<%.1f&&mu1_ProbNNghost<%.1f&&mu0_ProbNNmu>%.1f&&mu1_ProbNNmu>%.1f",ghostProbCut,ghostProbCut,muonPID,muonPID); 
  TString cut_BDT=TString::Format("BDT>%.1f",BDT);

  TString cut_PID_Kpimumu = cut_hadronPID_Kpimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID ;
  TString cut_PID_pipimumu = cut_hadronPID_pipimumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID;
  TString cut_PID_KKmumu = cut_hadronPID_KKmumu + "&&" + cut_MuonPID + "&&" + cut_SlowpiPID;

  //append preselection 
  normalizationCut_Kpi+="&&"+cutNshared;
  normalizationCut_KK+="&&"+cutNshared;
  normalizationCut_pipi+="&&"+cutNshared;

  if(applyPID) {
    normalizationCut_Kpi+="&&"+cut_PID_Kpimumu;
    normalizationCut_pipi+="&&"+cut_PID_pipimumu;
    normalizationCut_KK+="&&"+cut_PID_KKmumu;
  }  
  if(applyBDT) {
    normalizationCut_Kpi+="&&"+cut_BDT;
    normalizationCut_KK+="&&"+cut_BDT;
    normalizationCut_pipi+="&&"+cut_BDT;
  }

 
   dcastyle();


  TFile* fOut= new TFile(nameTarget+".root","RECREATE"); 

  //histograms 

  TH1D* effTriggerKpimumu_KKmumuTrained;
  effTriggerKpimumu_KKmumuTrained= new TH1D("effTriggerKpimumu_KKmumuTrained","efficiency MC Trigger D->Kpimumu ",sizeof(binsKpi)/sizeof(double)-1,binsKpi);
  
  TH1D* effTriggerKpimumu_pipimumuTrained;
  effTriggerKpimumu_pipimumuTrained= new TH1D("effTriggerKpimumu_pipimumuTrained","efficiency MC Trigger D->Kpimumu ",sizeof(binsKpi)/sizeof(double)-1,binsKpi);

  TH1D* effTriggerKKmumu_KKmumuTrained;
  effTriggerKKmumu_KKmumuTrained= new TH1D("effTriggerKKmumu_KKmumuTrained","efficiency MC Trigger D->KKmumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1D* effTriggerpipimumu_pipimumuTrained;
  effTriggerpipimumu_pipimumuTrained= new TH1D("effTriggerKpipimumu_pipimumuTrained","efficiency MC Trigger D->pipimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1D* relEffTriggerKKmumu_KKmumuTrained;
  relEffTriggerKKmumu_KKmumuTrained= new TH1D("relEffTriggerKKmumu_KKmumuTrained","relEfficiency MC Trigger D->KKmumu magUp",sizeof(binsKK)/sizeof(double)-1,binsKK);

  TH1D* relEffTriggerpipimumu_pipimumuTrained;
  relEffTriggerpipimumu_pipimumuTrained= new TH1D("relEffTriggerpipimumu_pipimumuTrained","relEfficiency MC Trigger D->pipimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);


  /*
  CHECK
    TH1D* relEffTriggerKKmumu;
  TH1D* relEffTriggerpipimumu;

  relEffTriggerKKmumu= new TH1D("relEffTriggerDKKmumu","relative efficiency MC Trigger D->KKmumu/D->Kpimumu magDw",sizeof(binsKK)/sizeof(double)-1,binsKK);
  relEffTriggerpipimumu= new TH1D("relEffTriggerDpipimumu","relative efficiency MC Trigger D->pipimumu/D->Kpimumu magDw",sizeof(binspipi)/sizeof(double)-1,binspipi);
  */

  cout<<"CUTS Kpi: "<<selectionCut_Trigger_Kpi<<"  "<<normalizationCut_Kpi<<endl;
  cout<<"CUTS pipi: "<<selectionCut_Trigger_pipi<<"  "<<normalizationCut_pipi<<endl;

  for(int i=0; i<rangesKpi_low.size();++i){
    
    //mag up KKmumu trained
    TString f1=TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]);
    TString f2=TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]);
    std::pair<double,double> tempEff = getMCEfficiencyFrom2Files(f1,f2,"BDT_Tree",selectionCut_Trigger_Kpi,normalizationCut_Kpi);
    
    effTriggerKpimumu_KKmumuTrained->SetBinContent(i+1,tempEff.first);
    effTriggerKpimumu_KKmumuTrained->SetBinError(i+1,tempEff.second);
 
    //mag up pipimumu trained
    f1=TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]);
    f2=TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2Kpimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangesKpi_low[i],rangesKpi_high[i]);
    tempEff = getMCEfficiencyFrom2Files(f1,f2,"BDT_Tree",selectionCut_Trigger_Kpi,normalizationCut_Kpi);

    effTriggerKpimumu_pipimumuTrained->SetBinContent(i+1,tempEff.first);
    effTriggerKpimumu_pipimumuTrained->SetBinError(i+1,tempEff.second);

  }    

  for(int i=0; i<rangesKK_low.size();++i){
    
    TString f1=TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]);
    TString f2=TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2KKmumu_%.1f_%.1f_D2KKmumuBDT_magUp.root",rangesKK_low[i],rangesKK_high[i]);
    std::pair<double,double> tempEff = getMCEfficiencyFrom2Files(f1,f2,"BDT_Tree",selectionCut_Trigger_KK,normalizationCut_KK);
    
    effTriggerKKmumu_KKmumuTrained->SetBinContent(i+1,tempEff.first);
    effTriggerKKmumu_KKmumuTrained->SetBinError(i+1,tempEff.second);
    
    double relEff = tempEff.first/effTriggerKpimumu_KKmumuTrained->GetBinContent(1);
    double dRelEff = TMath::Sqrt( TMath::Power( tempEff.second/effTriggerKpimumu_KKmumuTrained->GetBinContent(1),2) + TMath::Power( (tempEff.first*effTriggerKpimumu_KKmumuTrained->GetBinError(1))/(effTriggerKpimumu_KKmumuTrained->GetBinContent(1)*effTriggerKpimumu_KKmumuTrained->GetBinContent(1)) ,2) );
    relEffTriggerKKmumu_KKmumuTrained->SetBinContent(i+1,relEff);
    relEffTriggerKKmumu_KKmumuTrained->SetBinError(i+1,dRelEff);

  }    

  for(int i=0; i<rangespipi_low.size();++i){

    TString f1=TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]);
    TString f2=TString::Format("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MCEfficiencyStudies/MC_D2pipimumu_%.1f_%.1f_D2pipimumuBDT_magUp.root",rangespipi_low[i],rangespipi_high[i]);
    std::pair<double,double> tempEff = getMCEfficiencyFrom2Files(f1,f2,"BDT_Tree",selectionCut_Trigger_pipi,normalizationCut_pipi);
    
    effTriggerpipimumu_pipimumuTrained->SetBinContent(i+1,tempEff.first);
    effTriggerpipimumu_pipimumuTrained->SetBinError(i+1,tempEff.second);

    double relEff = tempEff.first/effTriggerKpimumu_pipimumuTrained->GetBinContent(1);
    double dRelEff = TMath::Sqrt( TMath::Power( tempEff.second/effTriggerKpimumu_pipimumuTrained->GetBinContent(1),2) + TMath::Power( (tempEff.first*effTriggerKpimumu_pipimumuTrained->GetBinError(1))/(effTriggerKpimumu_pipimumuTrained->GetBinContent(1)*effTriggerKpimumu_pipimumuTrained->GetBinContent(1)) ,2) );
    relEffTriggerpipimumu_pipimumuTrained->SetBinContent(i+1,relEff);
    relEffTriggerpipimumu_pipimumuTrained->SetBinError(i+1,dRelEff);

  }    
  
  /*
  ////////////////////////////////////////////////                                              
  // data part
  // 
  //
  
  
  TString q2Cut= "D_DiMuon_Mass>950&&D_DiMuon_Mass<1100";
  TString q2RangeNormalizationMode="D_DiMuon_Mass>675&&D_DiMuon_Mass<875";

  TString dataCut_Kpi=  cutNshared+ "&&" + cut_PID_Kpimumu + "&&(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)&&deltaM>144.5&&deltaM<146.5" ;
  TString dataCut_pipi=  cutNshared+ "&&" + q2Cut + "&&" + cut_PID_pipimumu + "&&(mu0_L0MuonDecision_TOS||mu1_L0MuonDecision_TOS)&&deltaM>144.5&&deltaM<146.5";
  TString misIDCut = "mu0_ProbNNmu>0.5&&mu0_ProbNNghost<0.5";

  //configure the fitter
  D2hhmumuFitter1D myFitter1D;
  myFitter1D.setPathToSignalMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2pipimumu_BDT.root");
  myFitter1D.setPathToNormMC("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpimumu_D2pipimumuBDT.root");
  myFitter1D.setPathToKpipipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2pipimumuBDT.root");
  myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/MC_D2Kpipipi_D2pipimumuBDT.root");
  //  myFitter1D.setPathToHHpipiData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/finalPreselected/D2pipipipi_D2pipimumuBDT.root");                                              
  myFitter1D.setPathToNormData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/noPreselectionCuts/D2Kpimumu_D2pipimumuBDT_noCuts.root");
  myFitter1D.setPathToSignalData("/auto/data/mitzel/D2hhmumu/new/preselectedSamples/noPreselectionCuts/D2pipimumu_D2pipimumuBDT_noCuts.root");
  
  //didnt create extra files for misID and MC without preselection                                                                                             
  myFitter1D.fit_MC(dataCut_pipi,true,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/HLTEfficiency/MC_"+dataCut_pipi+q2Cut+".eps");
  myFitter1D.fit_normalization_MC(dataCut_Kpi,true,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/HLTEfficiency/MC_normMode_"+dataCut_Kpi+".eps");
  std::cout<<"Monte Carlo fits done.."<<std::endl;

  myFitter1D.fit_Kpipipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/HLTEfficiency/misID_Kpi"+misIDCut+q2RangeNormalizationMode+".eps");
  myFitter1D.fit_HHpipi_misID(misIDCut,true,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/HLTEfficiency/misID_pipi"+misIDCut+q2Cut+".eps"); ////q2 cut missing? check!   
  std::cout<<"misID fits done.."<<std::endl;

  double nSel_norm = myFitter1D.fit_normalization_Data(dataCut_Kpi+"&&"+selectionCut_Trigger_Kpi,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/HLTEfficiency/HLT_Eff_Kpimumu_fit1.eps");
  double nSel_pipi=myFitter1D.fit_resonant_Data(dataCut_pipi+"&&"+selectionCut_Trigger_pipi,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/HLTEfficiency/HLT_Eff_pipimumu_fit1.eps");

  double nTot_norm = myFitter1D.fit_normalization_Data(dataCut_Kpi,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/HLTEfficiency/HLT_Eff_Kpimumu_fit2.eps");
  double nTot_pipi=myFitter1D.fit_resonant_Data(dataCut_pipi,"/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/HLTEfficiency/HLT_Eff_pipimumu_fit2.eps");

  double Eff_norm = nSel_norm/nTot_norm;
  double Eff_pipi = nSel_pipi/nTot_pipi;

  double dEff_norm = 1/nTot_norm * TMath::Sqrt(nSel_norm*(1-(nSel_norm/nTot_norm) ) );
  double dEff_pipi = 1/nTot_pipi * TMath::Sqrt(nSel_pipi*(1-(nSel_pipi/nTot_pipi) ) );

  std::cout<<"Eff_norm  "<<  Eff_norm<<"+-"<< dEff_norm <<std::endl;
  std::cout<<"Eff_pipi  "<<  Eff_pipi<<"+-"<<  dEff_pipi<<std::endl;
  std::cout<<"ratio  "<<  Eff_pipi/Eff_norm <<std::endl;

  std::cout<<"Eff_norm MC "<< effTriggerKpimumu_pipimumuTrained->GetBinContent(1)   <<std::endl;
  std::cout<<"Eff_pipi MC "<< effTriggerpipimumu_pipimumuTrained->GetBinContent(4)   <<std::endl;
  std::cout<<"ratio  "<<   effTriggerpipimumu_pipimumuTrained->GetBinContent(4)/effTriggerKpimumu_pipimumuTrained->GetBinContent(1) <<std::endl;

  */


  //Drawing part 

  TCanvas *a  = new TCanvas("pipimumu","Eff pipimumu");
  a->Divide(3,2);
  a->cd(1);
  effTriggerKpimumu_KKmumuTrained->GetYaxis()->SetTitle("#epsilon_{HLT}(D->K#pi#mu#mu)");
  effTriggerKpimumu_KKmumuTrained->GetYaxis()->SetRangeUser(0.5,0.8);
  effTriggerKpimumu_KKmumuTrained->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  effTriggerKpimumu_KKmumuTrained->Draw();
  a->cd(2);
  effTriggerKpimumu_pipimumuTrained->GetYaxis()->SetTitle("#epsilon_{HLT}(D->K#pi#mu#mu)");
  effTriggerKpimumu_pipimumuTrained->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  effTriggerKpimumu_pipimumuTrained->GetYaxis()->SetRangeUser(0.5,0.8);
  effTriggerKpimumu_pipimumuTrained->Draw();
  a->cd(3);
  effTriggerKKmumu_KKmumuTrained->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  effTriggerKKmumu_KKmumuTrained->GetYaxis()->SetTitle("#epsilon_{HLT}(D->KK#mu#mu)");
  //effTriggerKKmumu_KKmumuTrained->GetYaxis()->SetRangeUser(0.4,0.7);
  effTriggerKKmumu_KKmumuTrained->Draw();
  a->cd(4);
  effTriggerpipimumu_pipimumuTrained->Draw();
  effTriggerpipimumu_pipimumuTrained->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  effTriggerpipimumu_pipimumuTrained->GetYaxis()->SetTitle("#epsilon_{HLT}(D->#pi#pi#mu#mu)");
  effTriggerpipimumu_pipimumuTrained->GetYaxis()->SetRangeUser(0.4,0.7);
  a->cd(5);
  relEffTriggerKKmumu_KKmumuTrained->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  relEffTriggerKKmumu_KKmumuTrained->GetYaxis()->SetTitle("#epsilon_{HLT}(D->KK#mu#mu)/eEpsilon_{HLT}(D->K#pi#mu#mu)");
  relEffTriggerKKmumu_KKmumuTrained->GetYaxis()->SetRangeUser(0.8,1.1);
  relEffTriggerKKmumu_KKmumuTrained->Draw();
  a->cd(6);
  relEffTriggerpipimumu_pipimumuTrained->GetXaxis()->SetTitle("m(#mu#mu)[MeV]");
  relEffTriggerpipimumu_pipimumuTrained->GetYaxis()->SetTitle("#epsilon_{HLT}(D->#pi#pi#mu#mu)/#epsilon_{HLT}(D->K#pi#mu#mu)");
  relEffTriggerpipimumu_pipimumuTrained->GetYaxis()->SetRangeUser(0.8,1.1);
  relEffTriggerpipimumu_pipimumuTrained->Draw();
  a->Print(nameTarget+".eps");

  TCanvas* b = new TCanvas("b","b");
  b->Divide(1,2);
  b->cd(1);
  relEffTriggerpipimumu_pipimumuTrained->Draw();
  b->cd(2);
  relEffTriggerKKmumu_KKmumuTrained->Draw();
  b->Print(nameTarget+"_2.eps");


  effTriggerKpimumu_KKmumuTrained->Write();
  effTriggerKpimumu_pipimumuTrained->Write();
  effTriggerKKmumu_KKmumuTrained->Write();
  effTriggerpipimumu_pipimumuTrained->Write();
  
  relEffTriggerpipimumu_pipimumuTrained->Write();
  relEffTriggerKKmumu_KKmumuTrained->Write();

  fOut->Write();

  /*
  effTriggerpipimumu_pipimumuTrained_magUp->GetYaxis()->SetRangeUser(0.3,1);
  effTriggerpipimumu_pipimumuTrained_magUp->Draw();
  //effTriggerpipimumu_pipimumuTrained_magDw->SetLineColor(kRed);
  //effTriggerpipimumu_pipimumuTrained_magDw->Draw("SAME");
  effTriggerpipimumu_pipimumuTrained_data->SetLineColor(kOrange);
  effTriggerpipimumu_pipimumuTrained_data->Draw("SAME");

  effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->SetLineColor(kViolet);
  effTriggerpipimumu_pipimumuTrained_TISTOS_magUp->Draw("SAME");

  a->cd(2);
  relEffTriggerpipimumu_magUp->GetYaxis()->SetRangeUser(0.8,1.5);
  relEffTriggerpipimumu_magUp->Draw();
  //relEffTriggerpipimumu_magDw->SetLineColor(kRed);
  //relEffTriggerpipimumu_magDw->Draw("SAME");
  relEffTriggerpipimumu_data->SetLineColor(kOrange);
  relEffTriggerpipimumu_data->Draw("SAME");
  relEffTriggerpipimumu_TISTOS_magUp->SetLineColor(kViolet);
  relEffTriggerpipimumu_TISTOS_magUp->Draw("SAME");
  
  a->cd(3);
  effTriggerKKmumu_KKmumuTrained_magUp->GetYaxis()->SetRangeUser(0.3,1);
  effTriggerKKmumu_KKmumuTrained_magUp->Draw();
  effTriggerKKmumu_KKmumuTrained_magDw->SetLineColor(kRed);
  effTriggerKKmumu_KKmumuTrained_magDw->Draw("SAME");
  a->cd(4);
  relEffTriggerKKmumu_magUp->GetYaxis()->SetRangeUser(0.5,1.5);
  relEffTriggerKKmumu_magUp->Draw();
  relEffTriggerKKmumu_magDw->SetLineColor(kRed);
  relEffTriggerKKmumu_magDw->Draw("SAME");
  a->cd(5);
  effTriggerKpimumu_KKmumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effTriggerKpimumu_KKmumuTrained_magUp->Draw();
  effTriggerKpimumu_KKmumuTrained_magDw->SetLineColor(kRed);
  effTriggerKpimumu_KKmumuTrained_magDw->Draw("SAME");
  a->cd(6);
  effTriggerKpimumu_pipimumuTrained_magUp->GetYaxis()->SetRangeUser(0.3,1.0);
  effTriggerKpimumu_pipimumuTrained_magUp->Draw();
  //effTriggerKpimumu_pipimumuTrained_magDw->SetLineColor(kRed);
  //effTriggerKpimumu_pipimumuTrained_magDw->Draw("SAME");
  
  effTriggerKpimumu_pipimumuTrained_data->SetLineColor(kOrange);
  effTriggerKpimumu_pipimumuTrained_data->Draw("SAME");

  effTriggerKpimumu_pipimumuTrained_TISTOS_magUp->SetLineColor(kViolet);
  effTriggerKpimumu_pipimumuTrained_TISTOS_magUp->Draw("SAME");

  a->Print(textFile+"5.eps");

  a->Write();
 

  
  TCanvas *c  = new TCanvas("Kpimumu","Eff Kpimumu");

  c->Divide(2);
  c->cd(1);
  effTriggerKpimumu_KKmumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effTriggerKpimumu_KKmumuTrained_magUp->Draw();
  effTriggerKpimumu_KKmumuTrained_magDw->SetLineColor(kRed);
  effTriggerKpimumu_KKmumuTrained_magDw->Draw("SAME");
  c->cd(2);
  effTriggerKpimumu_pipimumuTrained_magUp->GetYaxis()->SetRangeUser(0.5,1);
  effTriggerKpimumu_pipimumuTrained_magUp->Draw();
  effTriggerKpimumu_pipimumuTrained_magDw->SetLineColor(kRed);
  effTriggerKpimumu_pipimumuTrained_magDw->Draw("SAME");
  effTriggerKpimumu_pipimumuTrained_data->Draw("SAME");

  c->Print(textFile+"6.eps");

  c->Write();
  */
}

void drawTotalEfficiency(){

  TFile* fGenLevel = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/GenLevelCutEfficiency/generatorLevelCutEfficiency.root","OPEN");
  TFile* fStripping= new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/strippingEfficiency/RecoAndStrippingEfficiency_fitted_(Dst_BKGCAT<11||Dst_BKGCAT==60).root","OPEN");
  //TFile* fPID = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/totalPIDEfficiency_default.root","OPEN");
  TFile *fMuonPID = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalMuonPIDEfficiency_default.root","READ");
  TFile *fHadronPID = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/PIDEfficiency/finalHadronPIDEfficiency_default.root","READ");
  TFile* fTrigger = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/TriggerEffTagAndProbe/L0TriggerEfficiencies.root","OPEN");
  //TFile* fHLT = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/HltMCEfficiencies/MCHltTrigger_Efficiencies_Hlt1_and_Hlt2.root","OPEN");
  //TFile* fBDT = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MCBDTEfficiencies_afterPID_andTrigger.root","OPEN");
  TFile* fHLT = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MCBDT_and_Hlt_EfficienciesnoGhosts.root","OPEN");
  TFile* fBDT = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/BDTEfficiency/MCBDT_and_Hlt_EfficienciesnoGhosts.root","OPEN");


  double sysL0 = 0.013;
  double sysBDT = 0.014;
  double sysMonPID=0.006 ;
  double sysHadPID =0.005 ;
  double sysReco = 0.009;
    
  TString bins[5]={"<525","525-565","565-950","950-1100",">1100"};
  

  dcastyle();
  
  //TH1* totalRelPIDEff_KKmumu= (TH1D*)fPID->Get("totalRelPIDEff_KKmumu");
  //TH1* totalRelPIDEff_pipimumu= (TH1D*)fPID->Get("totalRelPIDEff_pipimumu");

  TH1* totalRelMuonPIDEff_KKmumu= (TH1D*)fMuonPID->Get("relEff_KKmumu");
  TH1* totalRelMuonPIDEff_pipimumu= (TH1D*)fMuonPID->Get("relEff_pipimumu");

  TH1* totalRelHadronPIDEff_KKmumu= (TH1D*)fHadronPID->Get("relEff_KKmumu_average");
  TH1* totalRelHadronPIDEff_pipimumu= (TH1D*)fHadronPID->Get("relEff_pipimumu_average");

  TH1* genLevelRelativeEff_KK= (TH1D*)fGenLevel->Get("genLevelRelativeEff_KK");
  TH1* genLevelRelativeEff_pipi= (TH1D*)fGenLevel->Get("genLevelRelativeEff_pipi");

  TH1* strippingRelativeEff_KK= (TH1D*)fStripping->Get("strippingRelativeEff_KK");
  TH1* strippingRelativeEff_pipi= (TH1D*)fStripping->Get("strippingRelativeEff_pipi");

  TH1* L0Trigger_RelativeEff_KKmumu = (TH1D*)fTrigger->Get("L0Trigger_RelativeEff_KKmumu");
  TH1* L0Trigger_RelativeEff_pipimumu = (TH1D*)fTrigger->Get("L0Trigger_RelativeEff_pipimumu");
  
  TH1* relEffBDTKKmumu =(TH1D*)fBDT->Get("relEffBDTDKKmumu");
  TH1* relEffBDTpipimumu =(TH1D*)fBDT->Get("relEffBDTDpipimumu");
 
  TH1* relEffTriggerKKmumu=(TH1D*)fHLT->Get("relEffHltDKKmumu");
  TH1* relEffTriggerpipimumu=(TH1D*)fHLT->Get("relEffHltDpipimumu");
  
  TH1* relEffHltAndBDTKKmumu=(TH1D*)fHLT->Get("relEffCombinedHltBDTDKKmumu");
  TH1* relEffHltAndBDTpipimumu=(TH1D*)fHLT->Get("relEffCombinedHltBDTDpipimumu");
 
  TH1* totalRelEff_KKmumu = new TH1D("totalRelEff_KKmumu","rel PID Eff KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* totalRelEff_pipimumu = new TH1D("totalRelEff_pipimumu","rel PID Eff pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);

  TH1* totalRelEff_KKmumu_onlyStatError = new TH1D("totalRelEff_KKmumu_onlyStatError","rel PID Eff KKmumu in bins of dimuon mass",sizeof(binsKK)/sizeof(double)-1,binsKK);
  TH1* totalRelEff_pipimumu_onlyStatError = new TH1D("totalRelEff_pipimumu_onlyStatError","rel PID Eff pipimumu in bins of dimuon mass",sizeof(binspipi)/sizeof(double)-1,binspipi);


  TCanvas* cPID = new TCanvas("cMuonPID","cMuonID");
  cPID->Divide(1,2);
  cPID->cd(1);
  totalRelMuonPIDEff_KKmumu->Draw();
  cPID->cd(2);
  totalRelMuonPIDEff_pipimumu->Draw();
  cPID->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/totalEfficiency/MuonPID.eps");

  TCanvas* cHPID = new TCanvas("cHadronPID","cHadronPID");
  cHPID->Divide(1,2);
  cHPID->cd(1);
  totalRelHadronPIDEff_KKmumu->Draw();
  cHPID->cd(2);
  totalRelHadronPIDEff_pipimumu->Draw();
  cHPID->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/totalEfficiency/HadronPID.eps");

  TCanvas* cgenLevel = new TCanvas("cgenLevel","cgenLevel");
  cgenLevel->Divide(1,2);
  cgenLevel->cd(1);
  genLevelRelativeEff_KK->Draw();
  cgenLevel->cd(2);
  genLevelRelativeEff_pipi->Draw();
  cgenLevel->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/totalEfficiency/genLevel.eps");

  TCanvas* cStrip = new TCanvas("cStrip","cStrip");
  cStrip->Divide(1,2);
  cStrip->cd(1);
  strippingRelativeEff_KK->Draw();
  cStrip->cd(2);
  strippingRelativeEff_pipi->Draw();
  cStrip->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/totalEfficiency/Strip.eps");

  TCanvas* cL0Trigger = new TCanvas("cL0Trigger","cL0Trigger");
  cL0Trigger->Divide(1,2);
  cL0Trigger->cd(1);
  L0Trigger_RelativeEff_KKmumu->Draw();
  cL0Trigger->cd(2);
  L0Trigger_RelativeEff_pipimumu->Draw();
  cL0Trigger->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/totalEfficiency/L0Trigger.eps");

  TCanvas* cHlt = new TCanvas("cHlt","cHlt");
  cHlt->Divide(1,2);
  cHlt->cd(1);
  relEffTriggerKKmumu->Draw();
  cHlt->cd(2);
  relEffTriggerpipimumu->Draw();
  cHlt->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/totalEfficiency/Hlt.eps");

  TCanvas* cBDT = new TCanvas("cBDT","cBDT");
  cBDT->Divide(1,2);
  cBDT->cd(1);
  relEffBDTKKmumu->Draw();
  cBDT->cd(2);
  relEffBDTpipimumu->Draw();
  cBDT->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/totalEfficiency/BDT.eps");


  ofstream myfile;
  TString textFile = "../img/EfficiencyStudies/totalEfficiency/summary";
  myfile.open(textFile+".txt");


  myfile<<"D->KKmumu"<<endl;

  for(int i=0; i<rangesKK_low.size();++i){

    double totalEff=0;
    double dTotalEff=0;
    double dTotalEff_sys=0;

    totalEff = totalRelMuonPIDEff_KKmumu->GetBinContent(i+1);
    totalEff *= totalRelHadronPIDEff_KKmumu->GetBinContent(i+1);
    totalEff *= genLevelRelativeEff_KK->GetBinContent(i+1);
    totalEff *= strippingRelativeEff_KK->GetBinContent(i+1);
    totalEff *= L0Trigger_RelativeEff_KKmumu->GetBinContent(i+1);
    //totalEff *= relEffBDTKKmumu->GetBinContent(i+1);
    //totalEff *= relEffTriggerKKmumu->GetBinContent(i+1);
    totalEff *= relEffHltAndBDTKKmumu->GetBinContent(i+1);
    myfile<<bins[i]<<std::setprecision(3)<<" &  "<<totalEff;

    dTotalEff = TMath::Power(totalRelMuonPIDEff_KKmumu->GetBinError(i+1)/totalRelMuonPIDEff_KKmumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(totalRelHadronPIDEff_KKmumu->GetBinError(i+1)/totalRelHadronPIDEff_KKmumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(genLevelRelativeEff_KK->GetBinError(i+1)/genLevelRelativeEff_KK->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(strippingRelativeEff_KK->GetBinError(i+1)/strippingRelativeEff_KK->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(L0Trigger_RelativeEff_KKmumu->GetBinError(i+1)/L0Trigger_RelativeEff_KKmumu->GetBinContent(i+1),2);
    //dTotalEff += TMath::Power(relEffBDTKKmumu->GetBinError(i+1)/relEffBDTKKmumu->GetBinContent(i+1),2);
    //dTotalEff += TMath::Power(relEffTriggerKKmumu->GetBinError(i+1)/relEffTriggerKKmumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(relEffHltAndBDTKKmumu->GetBinError(i+1)/relEffHltAndBDTKKmumu->GetBinContent(i+1),2);
    
    myfile<<" & "<<TMath::Sqrt(dTotalEff)*100;
 
    totalRelEff_KKmumu_onlyStatError->SetBinContent(i+1,totalEff);
    totalRelEff_KKmumu_onlyStatError->SetBinError(i+1,TMath::Sqrt(dTotalEff)*totalEff);

    //systematics
    dTotalEff += TMath::Power(sysL0/L0Trigger_RelativeEff_KKmumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(sysMonPID/totalRelMuonPIDEff_KKmumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(sysHadPID/totalRelHadronPIDEff_KKmumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(sysBDT/relEffBDTKKmumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(sysReco/strippingRelativeEff_KK->GetBinContent(i+1),2);

    dTotalEff_sys += TMath::Power(sysL0/L0Trigger_RelativeEff_KKmumu->GetBinContent(i+1),2);
    dTotalEff_sys += TMath::Power(sysMonPID/totalRelMuonPIDEff_KKmumu->GetBinContent(i+1),2);
    dTotalEff_sys += TMath::Power(sysHadPID/totalRelHadronPIDEff_KKmumu->GetBinContent(i+1),2);
    dTotalEff_sys += TMath::Power(sysBDT/relEffBDTKKmumu->GetBinContent(i+1),2);
    dTotalEff_sys += TMath::Power(sysReco/strippingRelativeEff_KK->GetBinContent(i+1),2);
    
    myfile<<" & "<<TMath::Sqrt(dTotalEff_sys)*100;
    myfile<<" & "<<TMath::Sqrt(dTotalEff)*100<<endl;
  
    dTotalEff = TMath::Sqrt(dTotalEff)*totalEff;

    totalRelEff_KKmumu->SetBinContent(i+1,totalEff); 
    totalRelEff_KKmumu->SetBinError(i+1,dTotalEff);
 
    myfile<<"bin "<<i<<" muon PID "<<totalRelMuonPIDEff_KKmumu->GetBinContent(i+1)<<"+-"<<totalRelMuonPIDEff_KKmumu->GetBinError(i+1)
	  <<" dR/R stat ="<<totalRelMuonPIDEff_KKmumu->GetBinError(i+1)/totalRelMuonPIDEff_KKmumu->GetBinContent(i+1)
          <<" dR/R sys ="<<sysMonPID/totalRelMuonPIDEff_KKmumu->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" hadron PID "<<totalRelHadronPIDEff_KKmumu->GetBinContent(i+1)<<"+-"<<totalRelHadronPIDEff_KKmumu->GetBinError(i+1)
	  <<" dR/R stat ="<<totalRelHadronPIDEff_KKmumu->GetBinError(i+1)/totalRelHadronPIDEff_KKmumu->GetBinContent(i+1)
	  <<" dR/R sys ="<<sysHadPID/totalRelHadronPIDEff_KKmumu->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" genlevel "<<genLevelRelativeEff_KK->GetBinContent(i+1)<<"+-"<<genLevelRelativeEff_KK->GetBinError(i+1)
	  <<" dR/R ="<<genLevelRelativeEff_KK->GetBinError(i+1)/genLevelRelativeEff_KK->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" strip "<<strippingRelativeEff_KK->GetBinContent(i+1)<<"+-"<<strippingRelativeEff_KK->GetBinError(i+1)
	  <<" dR/R ="<<strippingRelativeEff_KK->GetBinError(i+1)/strippingRelativeEff_KK->GetBinContent(i+1)
	  <<" dR/R sys ="<<sysReco/strippingRelativeEff_KK->GetBinContent(i+1)<<endl;    
    myfile<<"bin "<<i<<" L0 "<<L0Trigger_RelativeEff_KKmumu->GetBinContent(i+1)<<"+-"<<L0Trigger_RelativeEff_KKmumu->GetBinError(i+1)
	  <<" dR/R stat="<<L0Trigger_RelativeEff_KKmumu->GetBinError(i+1)/L0Trigger_RelativeEff_KKmumu->GetBinContent(i+1)
	  <<" dR/R sys ="<<sysL0/L0Trigger_RelativeEff_KKmumu->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" BDT "<<relEffBDTKKmumu->GetBinContent(i+1)<<"+-"<<relEffBDTKKmumu->GetBinError(i+1)
	  <<" dR/R stat="<<relEffBDTKKmumu->GetBinError(i+1)/relEffBDTKKmumu->GetBinContent(i+1)
	  <<" dR/R sys="<<sysBDT/relEffBDTKKmumu->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" HLT "<<relEffTriggerKKmumu->GetBinContent(i+1)<<"+-"<<relEffTriggerKKmumu->GetBinError(i+1)
	  <<" dR/R ="<<relEffTriggerKKmumu->GetBinError(i+1)/relEffTriggerKKmumu->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" total "<<totalRelEff_KKmumu->GetBinContent(i+1)<<"+-"<<totalRelEff_KKmumu->GetBinError(i+1)
	  <<" dR/R ="<<totalRelEff_KKmumu->GetBinError(i+1)/totalRelEff_KKmumu->GetBinContent(i+1)<<std::endl;
      //<<" dR/R add. sys = "<<0.022/totalRelEff_KKmumu->GetBinContent(i+1)<<endl;

   
  }



  myfile<<"D->pipimumu"<<endl;

  for(int i=0; i<rangespipi_low.size();++i){

    double totalEff=0;
    double dTotalEff=0;
    double dTotalEff_sys=0;

    totalEff = totalRelMuonPIDEff_pipimumu->GetBinContent(i+1);
    totalEff *= totalRelHadronPIDEff_pipimumu->GetBinContent(i+1);
    totalEff *= genLevelRelativeEff_pipi->GetBinContent(i+1);
    totalEff *= strippingRelativeEff_pipi->GetBinContent(i+1);
    totalEff *= L0Trigger_RelativeEff_pipimumu->GetBinContent(i+1);
    //totalEff *= relEffBDTpipimumu->GetBinContent(i+1);
    //totalEff *= relEffTriggerpipimumu->GetBinContent(i+1);
    totalEff *= relEffHltAndBDTpipimumu->GetBinContent(i+1);
    myfile<<bins[i]<<std::setprecision(3)<<" &  "<<totalEff;

    dTotalEff = TMath::Power(totalRelMuonPIDEff_pipimumu->GetBinError(i+1)/totalRelMuonPIDEff_pipimumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(totalRelHadronPIDEff_pipimumu->GetBinError(i+1)/totalRelHadronPIDEff_pipimumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(genLevelRelativeEff_pipi->GetBinError(i+1)/genLevelRelativeEff_pipi->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(strippingRelativeEff_pipi->GetBinError(i+1)/strippingRelativeEff_pipi->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(L0Trigger_RelativeEff_pipimumu->GetBinError(i+1)/L0Trigger_RelativeEff_pipimumu->GetBinContent(i+1),2);
    //dTotalEff += TMath::Power(relEffBDTpipimumu->GetBinError(i+1)/relEffBDTpipimumu->GetBinContent(i+1),2);
    //dTotalEff += TMath::Power(relEffTriggerpipimumu->GetBinError(i+1)/relEffTriggerpipimumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(relEffHltAndBDTpipimumu->GetBinError(i+1)/relEffHltAndBDTpipimumu->GetBinContent(i+1),2);
    
    myfile<<" & "<<TMath::Sqrt(dTotalEff)*100;
 
    totalRelEff_pipimumu_onlyStatError->SetBinContent(i+1,totalEff);
    totalRelEff_pipimumu_onlyStatError->SetBinError(i+1,TMath::Sqrt(dTotalEff)*totalEff);

    //systematics
    dTotalEff += TMath::Power(sysL0/L0Trigger_RelativeEff_pipimumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(sysMonPID/totalRelMuonPIDEff_pipimumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(sysHadPID/totalRelHadronPIDEff_pipimumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(sysBDT/relEffBDTpipimumu->GetBinContent(i+1),2);
    dTotalEff += TMath::Power(sysReco/strippingRelativeEff_pipi->GetBinContent(i+1),2);

    dTotalEff_sys += TMath::Power(sysL0/L0Trigger_RelativeEff_pipimumu->GetBinContent(i+1),2);
    dTotalEff_sys += TMath::Power(sysMonPID/totalRelMuonPIDEff_pipimumu->GetBinContent(i+1),2);
    dTotalEff_sys += TMath::Power(sysHadPID/totalRelHadronPIDEff_pipimumu->GetBinContent(i+1),2);
    dTotalEff_sys += TMath::Power(sysBDT/relEffBDTpipimumu->GetBinContent(i+1),2);
    dTotalEff_sys += TMath::Power(sysReco/strippingRelativeEff_pipi->GetBinContent(i+1),2);

    myfile<<" & "<<TMath::Sqrt(dTotalEff_sys)*100;
    myfile<<" & "<<TMath::Sqrt(dTotalEff)*100<<endl;

    dTotalEff = TMath::Sqrt(dTotalEff)*totalEff;

    totalRelEff_pipimumu->SetBinContent(i+1,totalEff); 
    totalRelEff_pipimumu->SetBinError(i+1,dTotalEff);
    
    myfile<<"bin "<<i<<" muon PID "<<totalRelMuonPIDEff_pipimumu->GetBinContent(i+1)<<"+-"<<totalRelMuonPIDEff_pipimumu->GetBinError(i+1)
	  <<" dR/R stat ="<<totalRelMuonPIDEff_pipimumu->GetBinError(i+1)/totalRelMuonPIDEff_pipimumu->GetBinContent(i+1)
          <<" dR/R sys ="<<sysMonPID/totalRelMuonPIDEff_pipimumu->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" hadron PID "<<totalRelHadronPIDEff_pipimumu->GetBinContent(i+1)<<"+-"<<totalRelHadronPIDEff_pipimumu->GetBinError(i+1)
	  <<" dR/R stat ="<<totalRelHadronPIDEff_pipimumu->GetBinError(i+1)/totalRelHadronPIDEff_pipimumu->GetBinContent(i+1)
	  <<" dR/R sys ="<<sysHadPID/totalRelHadronPIDEff_pipimumu->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" genlevel "<<genLevelRelativeEff_pipi->GetBinContent(i+1)<<"+-"<<genLevelRelativeEff_pipi->GetBinError(i+1)
	  <<" dR/R ="<<genLevelRelativeEff_pipi->GetBinError(i+1)/genLevelRelativeEff_pipi->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" strip "<<strippingRelativeEff_pipi->GetBinContent(i+1)<<"+-"<<strippingRelativeEff_pipi->GetBinError(i+1)
	  <<" dR/R ="<<strippingRelativeEff_pipi->GetBinError(i+1)/strippingRelativeEff_pipi->GetBinContent(i+1)
	  <<" dR/R sys ="<<sysReco/strippingRelativeEff_pipi->GetBinContent(i+1)<<endl;    
    myfile<<"bin "<<i<<" L0 "<<L0Trigger_RelativeEff_pipimumu->GetBinContent(i+1)<<"+-"<<L0Trigger_RelativeEff_pipimumu->GetBinError(i+1)
	  <<" dR/R stat="<<L0Trigger_RelativeEff_pipimumu->GetBinError(i+1)/L0Trigger_RelativeEff_pipimumu->GetBinContent(i+1)
	  <<" dR/R sys ="<<sysL0/L0Trigger_RelativeEff_pipimumu->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" BDT "<<relEffBDTpipimumu->GetBinContent(i+1)<<"+-"<<relEffBDTpipimumu->GetBinError(i+1)
	  <<" dR/R stat="<<relEffBDTpipimumu->GetBinError(i+1)/relEffBDTpipimumu->GetBinContent(i+1)
	  <<" dR/R sys="<<sysBDT/relEffBDTpipimumu->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" HLT "<<relEffTriggerpipimumu->GetBinContent(i+1)<<"+-"<<relEffTriggerpipimumu->GetBinError(i+1)
	  <<" dR/R ="<<relEffTriggerpipimumu->GetBinError(i+1)/relEffTriggerpipimumu->GetBinContent(i+1)<<endl;
    myfile<<"bin "<<i<<" total "<<totalRelEff_pipimumu->GetBinContent(i+1)<<"+-"<<totalRelEff_pipimumu->GetBinError(i+1)
	  <<" dR/R ="<<totalRelEff_pipimumu->GetBinError(i+1)/totalRelEff_pipimumu->GetBinContent(i+1)<<std::endl;
      //<<" dR/R add. sys = "<<0.022/totalRelEff_pipimumu->GetBinContent(i+1)<<endl;

   
  }


  myfile.close();
  
  TCanvas* a = new TCanvas("a","a");
  a->Divide(1,2);
  a->cd(1);
  totalRelEff_KKmumu->GetYaxis()->SetRangeUser(0.4,0.6);
  totalRelEff_KKmumu->Draw();

  a->cd(2);
  totalRelEff_pipimumu->GetYaxis()->SetRangeUser(0.6,1.7);
  totalRelEff_pipimumu->Draw();
  a->Print("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/totalEfficiency/total.eps");

  TFile* fOut_final  = new TFile("/work/mitzel/D2hhmumu/dev/img/EfficiencyStudies/totalEfficiency/totalRelativeEfficiency.root","RECREATE");
  fOut_final->cd();
  totalRelEff_KKmumu->Write();
  totalRelEff_pipimumu->Write();

  totalRelEff_pipimumu_onlyStatError->Write();
  totalRelEff_KKmumu_onlyStatError->Write();

  fOut_final->Write();
  fOut_final->Close();
 
  

}
