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
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

using namespace std;
using namespace RooFit ;
using namespace RooStats;

int preselect(bool MC=false , bool preselection =true) {
        
    ///Load file	
    TChain* tree=new TChain("Bu2kpipimumuTuple/DecayTree");
    tree->Add("/auto/data/dargent/Bu2psiKpipi/data/11D/*.root");
    tree->Add("/auto/data/dargent/Bu2psiKpipi/data/11U/*.root");
    tree->Add("/auto/data/dargent/Bu2psiKpipi/data/12D/*.root");
    tree->Add("/auto/data/dargent/Bu2psiKpipi/data/12U/*.root");

    int N = tree->GetEntries();
    cout << "Old file contains " << N << " events" <<  endl;
    
    //Disable all branches
    tree->SetBranchStatus("*",0);
        
    //activate needed branches
    tree->SetBranchStatus("Bplus_L0MuonDecision_*",1) ;
    tree->SetBranchStatus("Bplus_L0DiMuonDecision_*",1) ;
    
	tree->SetBranchStatus("Bplus_Hlt1TrackAllL0Decision_*",1) ;
	tree->SetBranchStatus("Bplus_Hlt1TrackMuonDecision_*",1) ;
    tree->SetBranchStatus("Bplus_Hlt1DiMuonHighMassDecision_*",1) ;

	tree->SetBranchStatus("Bplus_Hlt2Topo2BodyBBDTDecision_*",1) ;
	tree->SetBranchStatus("Bplus_Hlt2Topo3BodyBBDTDecision_*",1) ;
	tree->SetBranchStatus("Bplus_Hlt2Topo4BodyBBDTDecision_*",1) ;
	tree->SetBranchStatus("Bplus_Hlt2TopoMu2BodyBBDTDecision_*",1) ;
	tree->SetBranchStatus("Bplus_Hlt2TopoMu3BodyBBDTDecision_*",1) ;
	tree->SetBranchStatus("Bplus_Hlt2TopoMu4BodyBBDTDecision_*",1) ;
	tree->SetBranchStatus("Bplus_Hlt2SingleMuonDecision_*",1) ;
	tree->SetBranchStatus("Bplus_Hlt2DiMuonDetached*",1) ;
    tree->SetBranchStatus("Bplus_Hlt2DiMuon*Psi2S*",1) ;
    tree->SetBranchStatus("Bplus_Hlt2DiMuon*JPsi*",1) ;

    tree->SetBranchStatus("*PID*",1) ;
    tree->SetBranchStatus("*ProbNN*",1) ;
    //tree->SetBranchStatus("*GhostProb",1) ;
    tree->SetBranchStatus("mu*_isMuon",1) ;

	tree->SetBranchStatus("nCandidate",1) ;
	tree->SetBranchStatus("nTracks",1) ;
	tree->SetBranchStatus("nPV",1) ;
	tree->SetBranchStatus("eventNumber",1) ;
	tree->SetBranchStatus("runNumber",1) ;
	tree->SetBranchStatus("EventInSequence",1) ;
	tree->SetBranchStatus("totCandidates",1) ;
	tree->SetBranchStatus("Polarity",1) ;
	
    tree->SetBranchStatus("*M",1) ;
	tree->SetBranchStatus("*MM",1) ;
    tree->SetBranchStatus( "*P", 1 );
	tree->SetBranchStatus( "*PX", 1 );
	tree->SetBranchStatus( "*PY", 1);
	tree->SetBranchStatus( "*PZ", 1);
	tree->SetBranchStatus( "*PE", 1);
	tree->SetBranchStatus( "*PT", 1 );
	tree->SetBranchStatus( "*TAU", 1 );
    tree->SetBranchStatus( "*ETA", 1 );
    tree->SetBranchStatus( "*FD*", 1 );
	tree->SetBranchStatus( "*IP*", 1 );

	tree->SetBranchStatus( "*IPCHI2_OWNPV", 1 );
	tree->SetBranchStatus( "*FDCHI2_OWNPV",1 );
	tree->SetBranchStatus( "*DIRA_OWNPV",1);
	tree->SetBranchStatus( "*ENDVERTEX_CHI2",1 );
    tree->SetBranchStatus( "*ENDVERTEX_NDOF",1 );
	tree->SetBranchStatus("*TRACK_CHI2NDOF",1) ;
    
    tree->SetBranchStatus( "*chi2",1 );
    tree->SetBranchStatus( "*nDOF",1 );
    tree->SetBranchStatus( "*status",1 );
    tree->SetBranchStatus( "*ctau",1 );

	if(!MC)tree->SetBranchStatus("*TRUE*",0) ;
    else tree->SetBranchStatus("*TRUE*",1) ;
    
	///Create output file
	TFile* output = new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected.root","RECREATE");
		
	///cuts
	string cuts;
	string mass_cuts = "Bplus_MM>5200 && Bplus_MM<5600 && abs(J_psi_1S_MM - 3686.109) < 60 ";
        //string kin_cut = "Bplus_TAU>0.00035" ;
        string kin_cut = "Bplus_PT>2000 && Bplus_TAU>0.0" ;

    string trigger_cuts = "Bplus_L0MuonDecision_TOS && (Bplus_Hlt1TrackAllL0Decision_TOS || Bplus_Hlt1TrackMuonDecision_TOS) && (Bplus_Hlt2Topo2BodyBBDTDecision_TOS || Bplus_Hlt2Topo3BodyBBDTDecision_TOS || Bplus_Hlt2Topo4BodyBBDTDecision_TOS || Bplus_Hlt2TopoMu2BodyBBDTDecision_TOS || Bplus_Hlt2TopoMu3BodyBBDTDecision_TOS || Bplus_Hlt2TopoMu4BodyBBDTDecision_TOS || Bplus_Hlt2SingleMuonDecision_TOS || Bplus_Hlt2DiMuonDetachedDecision_TOS)";
    
    string Ztrigger_cuts = "(Bplus_L0MuonDecision_TOS || Bplus_L0DiMuonDecision_TOS) && (Bplus_Hlt1TrackAllL0Decision_TOS || Bplus_Hlt1TrackMuonDecision_TOS || Bplus_Hlt1DiMuonHighMassDecision_TOS) && (Bplus_Hlt2Topo2BodyBBDTDecision_TOS || Bplus_Hlt2Topo3BodyBBDTDecision_TOS || Bplus_Hlt2Topo4BodyBBDTDecision_TOS ||Bplus_Hlt2DiMuonJPsiDecision_TOS || Bplus_Hlt2DiMuonJPsiHighPTDecision_TOS || Bplus_Hlt2DiMuonDetachedHeavyDecision_TOS || Bplus_Hlt2DiMuonPsi2SDecision_TOS || Bplus_Hlt2DiMuonPsi2SHighPTDecision_TOS || Bplus_Hlt2DiMuonDetachedDecision_TOS)";
    
    string Newtrigger_cuts = "(Bplus_L0MuonDecision_TOS || Bplus_L0DiMuonDecision_TOS) && (Bplus_Hlt1TrackAllL0Decision_TOS || Bplus_Hlt1TrackMuonDecision_TOS || Bplus_Hlt1DiMuonHighMassDecision_TOS) && (Bplus_Hlt2Topo2BodyBBDTDecision_TOS || Bplus_Hlt2Topo3BodyBBDTDecision_TOS || Bplus_Hlt2Topo4BodyBBDTDecision_TOS || Bplus_Hlt2TopoMu2BodyBBDTDecision_TOS || Bplus_Hlt2TopoMu3BodyBBDTDecision_TOS || Bplus_Hlt2TopoMu4BodyBBDTDecision_TOS || Bplus_Hlt2SingleMuonDecision_TOS || Bplus_Hlt2DiMuonDetachedDecision_TOS || Bplus_Hlt2DiMuonJPsiDecision_TOS || Bplus_Hlt2DiMuonJPsiHighPTDecision_TOS || Bplus_Hlt2DiMuonDetachedHeavyDecision_TOS || Bplus_Hlt2DiMuonPsi2SDecision_TOS || Bplus_Hlt2DiMuonPsi2SHighPTDecision_TOS || Bplus_Hlt2DiMuonDetachedDecision_TOS)";
    
    string quality_cuts = "(Bplus_ENDVERTEX_CHI2)/(Bplus_ENDVERTEX_NDOF) < 5 && piminus_TRACK_CHI2NDOF < 5 && piplus_TRACK_CHI2NDOF <5 && muminus_TRACK_CHI2NDOF < 5 && muplus_TRACK_CHI2NDOF < 5 && Kplus_TRACK_CHI2NDOF < 5";// && (Bplus_DTF_chi2 / Bplus_DTF_nDOF) < 10  && Bplus_DTF_status  == 0" ;
    string B0_cut = "abs(sqrt((J_psi_1S_PE + Kplus_PE+ piminus_PE)*(J_psi_1S_PE + Kplus_PE+ piminus_PE)-((J_psi_1S_PX + Kplus_PX+ piminus_PX)*(J_psi_1S_PX + Kplus_PX+ piminus_PX) + (J_psi_1S_PY + Kplus_PY+ piminus_PY)*(J_psi_1S_PY + Kplus_PY+ piminus_PY) + (J_psi_1S_PZ + Kplus_PZ + piminus_PZ)*(J_psi_1S_PZ + Kplus_PZ + piminus_PZ)))-5279.55)>60.";
    string pid_cuts = "Kplus_PIDK > 3.5 && piplus_PIDK < 14.5 && piminus_PIDK < 14.5 && (Kplus_PIDK - piplus_PIDK) > 10 && muplus_PIDmu > 0 &&  muminus_PIDmu > 0 && muplus_isMuon  && muminus_isMuon " ;
    
    string momenta_cutsZ = "Bplus_PT > 2000 && Kplus_PT > 260 && piplus_PT > 260 && piminus_PT > 260 && J_psi_1S_PT > 2000 && muplus_PT > 1000 && muminus_PT > 1000 " ;

    string momenta_cuts = "Bplus_PT > 2000 && Kplus_PT > 200 && piplus_PT > 200 && piminus_PT > 200" ;
    
    string MCTRUE_cut = "abs(Bplus_TRUEID)==521 && abs(J_psi_1S_TRUEID)==443 && abs(Kplus_TRUEID)==321 && abs(piplus_TRUEID)==211 && abs(piminus_TRUEID)==211" ;
    
    
	cuts.append(mass_cuts);
    if (preselection) {
        cuts.append("&&");
        cuts.append(kin_cut);
        cuts.append("&&");
        cuts.append(Newtrigger_cuts);
        cuts.append("&&");
        cuts.append(quality_cuts);
        cuts.append("&&");
        cuts.append(B0_cut);
        cuts.append("&&");
        cuts.append(pid_cuts);
        //cuts.append("&&");
        //cuts.append(momenta_cuts);
    }
    
    if (MC) {
        cuts.append("&&");
        cuts.append(MCTRUE_cut);
    }
    
    TTree* tree_sel = tree->CopyTree(cuts.c_str());    
    cout << "New file contains " << tree_sel->GetEntries() << " events" <<  endl;
    if(MC)cout << "Signal Eff = " << (double) ((double) tree_sel->GetEntries() / (double) tree->GetEntries()) << endl;
    
    output->WriteTObject(tree_sel);
	output->Close();
    delete output;
     
}    

void chooseBestPV(string fit="psiDTF", string input="/auto/data/dargent/Bu2psiKpipi/data/data_preselected.root"){
 
    double unit=1000000.;

    ///Read file
    TChain* tree=new TChain("DecayTree");
    tree->Add(input.c_str());
    
    ///momenta
    int nPV;
    float chi2[100], ndof[100];

    float B_M[100] ,K1_M[100] , B_TAU[100];
    float K_PX[100], K_PY[100], K_PZ[100], K_PE[100];
    float pip_PX[100], pip_PY[100], pip_PZ[100], pip_PE[100];
    float pim_PX[100], pim_PY[100], pim_PZ[100], pim_PE[100];
    float mup_PX[100], mup_PY[100], mup_PZ[100], mup_PE[100];
    float mum_PX[100], mum_PY[100], mum_PZ[100], mum_PE[100];

    tree->SetBranchAddress(("Bplus_"+fit+"_nPV").c_str(),&nPV);
    tree->SetBranchAddress(("Bplus_"+fit+"_chi2").c_str(),&chi2);
    tree->SetBranchAddress(("Bplus_"+fit+"_nDOF").c_str(),&ndof);

    tree->SetBranchAddress(("Bplus_"+fit+"_M").c_str(),&B_M);
    tree->SetBranchAddress(("Bplus_"+fit+"_ctau").c_str(),&B_TAU);

    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_Kplus_PX").c_str(),&K_PX);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_Kplus_PY").c_str(),&K_PY);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_Kplus_PZ").c_str(),&K_PZ); 
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_Kplus_PE").c_str(),&K_PE); 

    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_PX").c_str(),&pip_PX);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_PY").c_str(),&pip_PY);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_PZ").c_str(),&pip_PZ); 
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_PE").c_str(),&pip_PE); 

    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_0_PX").c_str(),&pim_PX);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_0_PY").c_str(),&pim_PY);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_0_PZ").c_str(),&pim_PZ); 
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_0_PE").c_str(),&pim_PE); 

    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_0_PX").c_str(),&mup_PX);
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_0_PY").c_str(),&mup_PY);
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_0_PZ").c_str(),&mup_PZ); 
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_0_PE").c_str(),&mup_PE); 

    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_PX").c_str(),&mum_PX);
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_PY").c_str(),&mum_PY);
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_PZ").c_str(),&mum_PZ); 
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_PE").c_str(),&mum_PE); 

    ///create new tree
    TFile* f = new TFile(("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_bestPV_"+fit+".root").c_str(),"RECREATE");
    TTree* new_tree = tree->CloneTree();//CopyTree();    
    double DTF_Bplus_M, DTF_psi_M, DTF_chi2, DTF_Bplus_TAU; 
    double DTF_Bplus_PX, DTF_Bplus_PY, DTF_Bplus_PZ, DTF_Bplus_PE;      
    double DTF_Kplus_PX, DTF_Kplus_PY, DTF_Kplus_PZ, DTF_Kplus_PE;      
    double DTF_piplus_PX, DTF_piplus_PY, DTF_piplus_PZ, DTF_piplus_PE;      
    double DTF_piminus_PX, DTF_piminus_PY, DTF_piminus_PZ, DTF_piminus_PE;      
    double DTF_muplus_PX, DTF_muplus_PY, DTF_muplus_PZ, DTF_muplus_PE;      
    double DTF_muminus_PX, DTF_muminus_PY, DTF_muminus_PZ, DTF_muminus_PE; 
    double DTF_Jpsi_1S_PX, DTF_Jpsi_1S_PY, DTF_Jpsi_1S_PZ, DTF_Jpsi_1S_PE;      

    TBranch* Bra_DTF_Bplus_M = new_tree->Branch((fit+"_Bplus_M").c_str(), &DTF_Bplus_M, (fit+"_Bplus_M/D").c_str());
    TBranch* Bra_DTF_Bplus_TAU = new_tree->Branch((fit+"_TAU").c_str(), &DTF_Bplus_TAU, (fit+"_TAU/D").c_str());
    TBranch* Bra_DTF_psi_M = new_tree->Branch((fit+"_psi_M").c_str(), &DTF_psi_M,  (fit+"_psi_M/D").c_str());
    TBranch* Bra_DTF_chi2 = new_tree->Branch((fit+"_chi2").c_str(), &DTF_chi2,  (fit+"_chi2/D").c_str());

    TBranch* Bra_DTF_Bplus_PX = new_tree->Branch((fit+"_Bplus_PX").c_str(), &DTF_Bplus_PX, (fit+"_Bplus_PX/D").c_str());
    TBranch* Bra_DTF_Bplus_PY = new_tree->Branch((fit+"_Bplus_PY").c_str(), &DTF_Bplus_PY, (fit+"_Bplus_PY/D").c_str());
    TBranch* Bra_DTF_Bplus_PZ = new_tree->Branch((fit+"_Bplus_PZ").c_str(), &DTF_Bplus_PZ, (fit+"_Bplus_PZ/D").c_str());
    TBranch* Bra_DTF_Bplus_PE = new_tree->Branch((fit+"_Bplus_PE").c_str(), &DTF_Bplus_PE, (fit+"_Bplus_PE/D").c_str());

    TBranch* Bra_DTF_Kplus_PX = new_tree->Branch((fit+"_Kplus_PX").c_str(), &DTF_Kplus_PX, (fit+"_Kplus_PX/D").c_str());
    TBranch* Bra_DTF_Kplus_PY = new_tree->Branch((fit+"_Kplus_PY").c_str(), &DTF_Kplus_PY, (fit+"_Kplus_PY/D").c_str());
    TBranch* Bra_DTF_Kplus_PZ = new_tree->Branch((fit+"_Kplus_PZ").c_str(), &DTF_Kplus_PZ, (fit+"_Kplus_PZ/D").c_str());
    TBranch* Bra_DTF_Kplus_PE = new_tree->Branch((fit+"_Kplus_PE").c_str(), &DTF_Kplus_PE, (fit+"_Kplus_PE/D").c_str());

    TBranch* Bra_DTF_piplus_PX = new_tree->Branch((fit+"_piplus_PX").c_str(), &DTF_piplus_PX, (fit+"_piplus_PX/D").c_str());
    TBranch* Bra_DTF_piplus_PY = new_tree->Branch((fit+"_piplus_PY").c_str(), &DTF_piplus_PY, (fit+"_piplus_PY/D").c_str());
    TBranch* Bra_DTF_piplus_PZ = new_tree->Branch((fit+"_piplus_PZ").c_str(), &DTF_piplus_PZ, (fit+"_piplus_PZ/D").c_str());
    TBranch* Bra_DTF_piplus_PE = new_tree->Branch((fit+"_piplus_PE").c_str(), &DTF_piplus_PE, (fit+"_piplus_PE/D").c_str());

    TBranch* Bra_DTF_piminus_PX = new_tree->Branch((fit+"_piminus_PX").c_str(), &DTF_piminus_PX, (fit+"_piminus_PX/D").c_str());
    TBranch* Bra_DTF_piminus_PY = new_tree->Branch((fit+"_piminus_PY").c_str(), &DTF_piminus_PY, (fit+"_piminus_PY/D").c_str());
    TBranch* Bra_DTF_piminus_PZ = new_tree->Branch((fit+"_piminus_PZ").c_str(), &DTF_piminus_PZ, (fit+"_piminus_PZ/D").c_str());
    TBranch* Bra_DTF_piminus_PE = new_tree->Branch((fit+"_piminus_PE").c_str(), &DTF_piminus_PE, (fit+"_piminus_PE/D").c_str());

    TBranch* Bra_DTF_muplus_PX = new_tree->Branch((fit+"_muplus_PX").c_str(), &DTF_muplus_PX, (fit+"_muplus_PX/D").c_str());
    TBranch* Bra_DTF_muplus_PY = new_tree->Branch((fit+"_muplus_PY").c_str(), &DTF_muplus_PY, (fit+"_muplus_PY/D").c_str());
    TBranch* Bra_DTF_muplus_PZ = new_tree->Branch((fit+"_muplus_PZ").c_str(), &DTF_muplus_PZ, (fit+"_muplus_PZ/D").c_str());
    TBranch* Bra_DTF_muplus_PE = new_tree->Branch((fit+"_muplus_PE").c_str(), &DTF_muplus_PE, (fit+"_muplus_PE/D").c_str());

    TBranch* Bra_DTF_muminus_PX = new_tree->Branch((fit+"_muminus_PX").c_str(), &DTF_muminus_PX, (fit+"_muminus_PX/D").c_str());
    TBranch* Bra_DTF_muminus_PY = new_tree->Branch((fit+"_muminus_PY").c_str(), &DTF_muminus_PY, (fit+"_muminus_PY/D").c_str());
    TBranch* Bra_DTF_muminus_PZ = new_tree->Branch((fit+"_muminus_PZ").c_str(), &DTF_muminus_PZ, (fit+"_muminus_PZ/D").c_str());
    TBranch* Bra_DTF_muminus_PE = new_tree->Branch((fit+"_muminus_PE").c_str(), &DTF_muminus_PE, (fit+"_muminus_PE/D").c_str());

    TBranch* Bra_DTF_Jpsi_1S_PX = new_tree->Branch((fit+"_Jpsi_1S_PX").c_str(), &DTF_Jpsi_1S_PX, (fit+"_Jpsi_1S_PX/D").c_str());
    TBranch* Bra_DTF_Jpsi_1S_PY = new_tree->Branch((fit+"_Jpsi_1S_PY").c_str(), &DTF_Jpsi_1S_PY, (fit+"_Jpsi_1S_PY/D").c_str());
    TBranch* Bra_DTF_Jpsi_1S_PZ = new_tree->Branch((fit+"_Jpsi_1S_PZ").c_str(), &DTF_Jpsi_1S_PZ, (fit+"_Jpsi_1S_PZ/D").c_str());
    TBranch* Bra_DTF_Jpsi_1S_PE = new_tree->Branch((fit+"_Jpsi_1S_PE").c_str(), &DTF_Jpsi_1S_PE, (fit+"_Jpsi_1S_PE/D").c_str());

   ///loop over events
    int numEvents = tree->GetEntries();
    
    for(int i=0; i< numEvents; i++)
    {	
	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
	tree->GetEntry(i);

	///find best PV
	int best=0;
	double tmp=chi2[0]/ndof[0];
	if (0ul == (i % 10000ul))cout << "nPV= " << nPV << endl;
	if (0ul == (i % 10000ul))cout << "first= " << tmp << endl;

	for(int j=1; j<nPV ;j++){
                if (0ul == (i % 10000ul))cout << j << "_chi2= " << chi2[j]/ndof[j] << endl;
		if(chi2[j]/ndof[j]< tmp){
			 tmp = chi2[j]/ndof[j];
			 best=j;
		}
	}

	if (0ul == (i % 10000ul)){
		cout << "best= " << tmp << endl;
		cout << endl ;
		cout << endl ;
	}

        TLorentzVector K_P(K_PX[best],K_PY[best],K_PZ[best],K_PE[best]);
        TLorentzVector pip_P(pip_PX[best],pip_PY[best],pip_PZ[best],pip_PE[best]);
        TLorentzVector pim_P(pim_PX[best],pim_PY[best],pim_PZ[best],pim_PE[best]);
        TLorentzVector mup_P(mup_PX[best],mup_PY[best],mup_PZ[best],mup_PE[best]);
        TLorentzVector mum_P(mum_PX[best],mum_PY[best],mum_PZ[best],mum_PE[best]);

	TLorentzVector B_P = K_P + pip_P + pim_P + mup_P + mum_P;
	TLorentzVector psi_P = mup_P + mum_P;
	

	///Fill branches
	DTF_Bplus_M = B_M[best];
	if(DTF_Bplus_M<5200. || DTF_Bplus_M>5600.)continue;
	DTF_Bplus_TAU = B_TAU[best];
	DTF_chi2 = tmp;
	DTF_psi_M = psi_P.M();

	DTF_Bplus_PX = B_P.X();
	DTF_Bplus_PY = B_P.Y();
	DTF_Bplus_PZ = B_P.Z();
	DTF_Bplus_PE = B_P.E();

	DTF_Jpsi_1S_PX = psi_P.X();
	DTF_Jpsi_1S_PY = psi_P.Y();
	DTF_Jpsi_1S_PZ = psi_P.Z();
	DTF_Jpsi_1S_PE = psi_P.E();

	DTF_Kplus_PX = K_PX[best];
	DTF_Kplus_PY = K_PY[best];
	DTF_Kplus_PZ = K_PZ[best];
	DTF_Kplus_PE = K_PE[best];

	DTF_piplus_PX = pip_PX[best];
	DTF_piplus_PY = pip_PY[best];
	DTF_piplus_PZ = pip_PZ[best];
	DTF_piplus_PE = pip_PE[best];

	DTF_piminus_PX = pim_PX[best];
	DTF_piminus_PY = pim_PY[best];
	DTF_piminus_PZ = pim_PZ[best];
	DTF_piminus_PE = pim_PE[best];

	DTF_muplus_PX = mup_PX[best];
	DTF_muplus_PY = mup_PY[best];
	DTF_muplus_PZ = mup_PZ[best];
	DTF_muplus_PE = mup_PE[best];

	DTF_muminus_PX = mum_PX[best];
	DTF_muminus_PY = mum_PY[best];
	DTF_muminus_PZ = mum_PZ[best];
	DTF_muminus_PE = mum_PE[best];

        Bra_DTF_Bplus_M->Fill();
        Bra_DTF_Bplus_TAU->Fill();
	Bra_DTF_psi_M->Fill();
        Bra_DTF_chi2->Fill();

        Bra_DTF_Bplus_PX->Fill();
        Bra_DTF_Bplus_PY->Fill();
        Bra_DTF_Bplus_PZ->Fill();
        Bra_DTF_Bplus_PE->Fill();

        Bra_DTF_Jpsi_1S_PX->Fill();
        Bra_DTF_Jpsi_1S_PY->Fill();
        Bra_DTF_Jpsi_1S_PZ->Fill();
        Bra_DTF_Jpsi_1S_PE->Fill();

        Bra_DTF_Kplus_PX->Fill();
        Bra_DTF_Kplus_PY->Fill();
        Bra_DTF_Kplus_PZ->Fill();
        Bra_DTF_Kplus_PE->Fill();

        Bra_DTF_piplus_PX->Fill();
        Bra_DTF_piplus_PY->Fill();
        Bra_DTF_piplus_PZ->Fill();
        Bra_DTF_piplus_PE->Fill();

        Bra_DTF_piminus_PX->Fill();
        Bra_DTF_piminus_PY->Fill();
        Bra_DTF_piminus_PZ->Fill();
        Bra_DTF_piminus_PE->Fill();

        Bra_DTF_muplus_PX->Fill();
        Bra_DTF_muplus_PY->Fill();
        Bra_DTF_muplus_PZ->Fill();
        Bra_DTF_muplus_PE->Fill();

        Bra_DTF_muminus_PX->Fill();
        Bra_DTF_muminus_PY->Fill();
        Bra_DTF_muminus_PZ->Fill();
        Bra_DTF_muminus_PE->Fill();

     }

    new_tree->Write();
    f->Close();
}


void fitPreselected(){

	bool binned=false;
	bool sWeight=true;

   	 gStyle->SetOptStat(0);
    	//gStyle->SetTitleXOffset(1.1);
    	//gStyle->SetTitleYOffset(1.3);
    	gStyle->SetTitleXSize(0.05);
    	gStyle->SetTitleYSize(0.05);
    	gStyle->SetTitleFont(42,"X");
    	gStyle->SetTitleFont(42,"Y");
    	//gStyle->SetLabelSize(0.033,"X");
    	//gStyle->SetLabelSize(0.033,"Y");
    	gStyle->SetLabelFont(42,"X");
    	gStyle->SetLabelFont(42,"Y");
   	gStyle->SetLabelOffset(0.01,"X");
    	gStyle->SetPadTickX(1);
    	gStyle->SetPadTickY(1);
    
    	TH1::SetDefaultSumw2();
    	TH2::SetDefaultSumw2();
	///Load file
	TFile* file;
	file= new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_bestPV_psiDTF.root");	
	//file= new TFile("data_test.root");
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bplus_MM",1);
	tree->SetBranchStatus("Bplus_TAU",1);
	tree->SetBranchStatus("psiFit_Bplus_M",1);
   	tree->SetBranchStatus("*PT",1);

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///B+
        RooRealVar Bplus_MM("Bplus_MM", "m(#psi K #pi #pi)", 5200., 5600.,"MeV");
	RooRealVar Bplus_PT("Bplus_PT","Bplus_PT",0.);
	RooRealVar Bplus_TAU("Bplus_TAU","Bplus_TAU",0.);
  	RooRealVar J_psi_1S_PT("J_psi_1S_PT", "J_psi_1S_PT", 0.);
  	RooRealVar muminus_PT("muminus_PT", "muminus_PT", 0.);
	RooRealVar muplus_PT("muplus_PT", "muplus_PT", 0.);
  	RooRealVar piplus_PT("piplus_PT","piplus_PT",0.);
	RooRealVar piminus_PT("piminus_PT","piminus_PT",0.);
  	RooRealVar Kplus_PT("Kplus_PT", "Kplus_PT", 0.);

	RooArgList list =  RooArgList(Bplus_MM,Bplus_PT,Bplus_TAU,J_psi_1S_PT,Kplus_PT,piplus_PT,piminus_PT,muplus_PT,muminus_PT);
	RooDataSet* data = new RooDataSet("data", "data", tree, list, "");
	 //"Bplus_PT>2000 && Kplus_PT>200 && piplus_PT>200 && piminus_PT>200 && Bplus_TAU >0.00035");
	
	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(Bplus_MM)));
	RooDataHist* data_binned = data_small->binnedClone();

	///Define fit model
	///----------------

	///Signal model
	///-----------------------

	RooRealVar mean1("mu", "mean1", 5281.71,5150.,5350.); 
	RooRealVar sigma1("sigma_{1}", "sigma1", 15.45,10.,20.);	
	RooRealVar sigma2("sigma_{2}", "sigma2", 29.,20.,40.);
	RooRealVar scale("s", "scale", 1.924,1.,3.);
	RooFormulaVar sigma3("sigma3", "s*sigma_{1}", RooArgSet(scale,sigma1));

	//cout<<sigma3.getVal()<<endl;

	///Gaussian
	RooGaussian Gauss1("Gauss1", "Gauss1", Bplus_MM, mean1, sigma1);
	RooGaussian Gauss2("Gauss2", "Gauss2", Bplus_MM, mean1, sigma2);
	//RooGaussian Gauss3("Gauss3", "Gauss3", DeltaM, mean2, sigma3);
	//RooRealVar f_sig("f_{sig}", "signal fraction", 0.56, 0., 1.);

	/// Crystal Ball
	RooRealVar alpha("alpha","alpha",1.41,0.5,10.);
	RooRealVar n("n","n",34.4);
	RooCBShape crystal1("crystal1","CB PDF1",Bplus_MM,mean1,sigma1, alpha, n);
        RooCBShape crystal2("crystal2","CB PDF2",Bplus_MM,mean1,sigma3, alpha, n);	

	///signal pdf
	RooRealVar f_sig("f_{sig}", "signal fraction", 0.812,0.,0.95);
	//RooAddPdf signal("signal", "signal", RooArgList(Gauss1,Gauss2),RooArgList(f_sig));
	RooAddPdf signal("signal", "signal", RooArgList(crystal1,crystal2),RooArgList(f_sig));

	///Background model
	///-------------------------

	///Exponential
	RooRealVar exp_par("lambda","exp",-0.003535,-1.,0.);	
	RooExponential bkg_exp("bkg_exp","exponential background",Bplus_MM,exp_par);
	
	///Polynomial
	RooRealVar c_0("c_{0}","c_0",0.);
	//c_0.setConstant(0.);
	RooRealVar c_1("c_{1}","c_1",-.356,-1.,0.);
	RooRealVar c_2("c_{2}","c_2",-0.01,-1.,1.);

	RooChebychev bkg_Chebychev("bkg_Chebychev","bkg_Chebychev", Bplus_MM, RooArgList(c_1));
	RooPolynomial bkg_Poly("bkg_Poly","bkg_Poly",Bplus_MM,RooArgList(c_0));

	///Gaussian
	RooRealVar mean_bkg("mu_{bkg}", "mean_bkg", 5063.,5000.,5100.); 
	RooRealVar sigma_bkg("sigma_{bkg}", "sigma_bkg", 48.5, 40., 80.);
	RooRealVar sigma2_bkg("sigma_{2}", "sigma2", 19.2, 0.1, 100.5);
	RooGaussian bkg_Gauss("bkg_Gauss", "bkg_Gauss", Bplus_MM, mean_bkg, sigma_bkg);
	
	///bkg pdf
	RooRealVar f_bkg("f_{bkg}", "background fraction", 0.11, 0., 1.);
	RooRealVar f2_bkg("f2_{bkg}", "background fraction2", 0.11, 0., 1.);

	//RooAddPdf bkg("bkg", "bkg", RooArgList(bkg_exp, bkg_Chebychev), RooArgList(f_bkg));
	RooAddPdf bkg("bkg", "bkg", RooArgList(bkg_Gauss, bkg_Chebychev), RooArgList(f_bkg));

	///total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(signal, bkg_Chebychev), RooArgList(n_sig,n_bkg));

	///Fit
	RooFitResult *result;
	if(binned) result = pdf->fitTo(*data_binned,Save(kTRUE),Extended());
	else result = pdf->fitTo(*data,Save(kTRUE),Extended(),NumCPU(3));
	cout << "result is --------------- "<<endl;
	result->Print(); 

	cout << "sigma2 =" << sigma1.getVal()*scale.getVal() << " pm " << sigma1.getError()*scale.getVal() << endl;

	///calculate # (signal)background events in signal region
	cout << endl;
	cout << endl;

	Bplus_MM.setRange("signal",mean1.getVal()-60.,mean1.getVal()+60.);
	
	RooAbsReal* S_fr= signal.createIntegral(Bplus_MM,NormSet(Bplus_MM),Range("signal"));
	Double_t S = S_fr->getVal()*n_sig.getVal();
	cout<<"S= " << S << endl;
	RooAbsReal* B_fr= bkg.createIntegral(Bplus_MM,NormSet(Bplus_MM),Range("signal"));
	Double_t B = B_fr->getVal()*n_bkg.getVal();
	cout<<"B= " << B << endl;
	
	cout<< "S/sqrt(S+B)= " << S/sqrt(S+B) << endl;
   	cout<<"S/B= " << S/B<< endl;

	cout << endl;
	cout << endl;

	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bplus_MM.frame();
	frame_m->SetTitle("");

	data->plotOn(frame_m,Name("data"),MarkerSize(0.1),Binning(100),DrawOption("e"),DataError(RooAbsData::SumW2));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack));
	pdf->plotOn(frame_m,Components(signal),LineColor(kBlue),LineStyle(kDashed));
	pdf->plotOn(frame_m,Components(bkg_Chebychev),LineColor(kRed),LineStyle(kDashed));
	pdf->paramOn(frame_m,Layout(0.6));
	frame_m->Draw();
	c1->Print("Preselection/BmassFit.eps");

	RooPlot* frame_m2= Bplus_MM.frame();
	frame_m2->SetTitle("");
	//data->plotOn(frame_m2,MarkerSize(0.2),Binning(100),DrawOption("e"),DataError(RooAbsData::SumW2));
	data->plotOn(frame_m2,MarkerSize(0.5),Binning(50));
	pdf->plotOn(frame_m2,LineColor(kBlack),LineWidth(2));
	pdf->plotOn(frame_m2,Components(signal),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m2,Components(bkg_Chebychev),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	data->plotOn(frame_m2,MarkerSize(0.5),Binning(50));
	frame_m2->Draw();
	c1->Print("Preselection/BmassFit2.eps");

	RooPlot* frame_m3= Bplus_MM.frame("");
	data->plotOn(frame_m3);
	pdf->plotOn(frame_m3);
	gPad->SetLogy();
	frame_m3->Draw();
	c1->Print("Preselection/BmassFit_log.eps");


	gPad->SetLogy(0);
	///fit results
	double chi2 = frame_m->chiSquare("pdf","data",9);
	double covmatr = result->covQual();
	double edm = result->edm();
	cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl;

  	// Construct a histogram with the pulls of the data w.r.t the curve
	RooHist* hresid = frame_m->residHist("data","pdf") ;
  	RooHist* hpull = frame_m->pullHist("data","pdf") ;

	// Create a new frame to draw the residual distribution and add the distribution to the frame
	RooPlot* frame2 = Bplus_MM.frame(Title("Residual Distribution")) ;
	frame2->SetTitle("");
	frame2->addPlotable(hresid,"P") ;
	frame2->Draw();
	c1->Print("Preselection/residual.eps");

	// Create a new frame to draw the pull distribution and add the distribution to the frame
	RooPlot* frame3 = Bplus_MM.frame(Title("Pull Distribution")) ;
	frame3->SetTitle("");	
	frame3->SetLabelFont(62,"Y");
	frame3->addPlotable(hpull,"P") ;
	frame3->Draw();
	c1->Print("Preselection/pull.eps");

	if(sWeight){
		mean1.setConstant();
		sigma1.setConstant();
		sigma2.setConstant();
		f_sig.setConstant();
		f_bkg.setConstant();
		exp_par.setConstant();
		c_0.setConstant();
		alpha.setConstant();
		n.setConstant();
		//sigma3.setConstant();
		c_1.setConstant();
		scale.setConstant();
	
		SPlot* sData = new SPlot("sData","An SPlot",*data,pdf,RooArgList(n_sig,n_bkg)); 
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs = new TCanvas("Bd sWeight","Bd sWeight distribution");
		TH2 * SwDsHist = (TH2*)data->createHistogram("Bplus_MM,n_sig_sw");
		SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist->SetTitle("");
		//SwDs->Write();
		SwDsHist->Draw();
		SwDs->Print("Preselection/Bd_sWeight.eps");

    		///Create output file
   		 TFile* output = new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_sweight.root","RECREATE");

		 tree->SetBranchStatus("*",1);
   		 TTree* new_tree = tree->CloneTree();//CopyTree();    
    		 double w;
    		 TBranch* Bra_sw = new_tree->Branch("n_sig_sw", &w, "n_sig_sw/D");

  		  ///loop over events
    		  int numEvents = tree->GetEntries();
    
    		  for(int i=0; i< numEvents; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
			tree->GetEntry(i);
			w=sData->GetSWeight(i,"n_sig_sw");
			Bra_sw->Fill();
  		  }
   		 new_tree->Write();
   		 output->Close();
	}
}



void fitPreselected_psiConstrained(){

	bool binned=false;
	bool sWeight=true;

   	 gStyle->SetOptStat(0);
    	//gStyle->SetTitleXOffset(1.1);
    	//gStyle->SetTitleYOffset(1.3);
    	gStyle->SetTitleXSize(0.05);
    	gStyle->SetTitleYSize(0.05);
    	gStyle->SetTitleFont(42,"X");
    	gStyle->SetTitleFont(42,"Y");
    	//gStyle->SetLabelSize(0.033,"X");
    	//gStyle->SetLabelSize(0.033,"Y");
    	gStyle->SetLabelFont(42,"X");
    	gStyle->SetLabelFont(42,"Y");
   	gStyle->SetLabelOffset(0.01,"X");
    	gStyle->SetPadTickX(1);
    	gStyle->SetPadTickY(1);
    
    	TH1::SetDefaultSumw2();
    	TH2::SetDefaultSumw2();
	///Load file
	TFile* file;
	file= new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_bestPV_psiDTF.root");	
	//file= new TFile("data_test.root");
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bplus_MM",1);
	tree->SetBranchStatus("Bplus_TAU",1);
	tree->SetBranchStatus("psiFit_Bplus_M",1);
   	tree->SetBranchStatus("*PT",1);

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///B+
        RooRealVar Bplus_MM("psiFit_Bplus_M", "m(#psi K #pi #pi)", 5200., 5600.,"MeV");
	RooRealVar Bplus_PT("Bplus_PT","Bplus_PT",0.);
	RooRealVar Bplus_TAU("Bplus_TAU","Bplus_TAU",0.);
  	RooRealVar J_psi_1S_PT("J_psi_1S_PT", "J_psi_1S_PT", 0.);
  	RooRealVar muminus_PT("muminus_PT", "muminus_PT", 0.);
	RooRealVar muplus_PT("muplus_PT", "muplus_PT", 0.);
  	RooRealVar piplus_PT("piplus_PT","piplus_PT",0.);
	RooRealVar piminus_PT("piminus_PT","piminus_PT",0.);
  	RooRealVar Kplus_PT("Kplus_PT", "Kplus_PT", 0.);

	RooArgList list =  RooArgList(Bplus_MM,Bplus_PT,Bplus_TAU,J_psi_1S_PT,Kplus_PT,piplus_PT,piminus_PT,muplus_PT,muminus_PT);
	RooDataSet* data = new RooDataSet("data", "data", tree, list, "");
	 //"Bplus_PT>2000 && Kplus_PT>200 && piplus_PT>200 && piminus_PT>200 && Bplus_TAU >0.00035");
	
	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(Bplus_MM)));
	RooDataHist* data_binned = data_small->binnedClone();

	///Define fit model
	///----------------

	///Signal model
	///-----------------------

	RooRealVar mean1("mu", "mean1", 5281.71,5150.,5350.); 
	RooRealVar sigma1("sigma_{1}", "sigma1", 15.45,0.,50.);	
	RooRealVar sigma2("sigma_{2}", "sigma2", 29.,20.,40.);
	RooRealVar scale("s", "scale", 1.924,1.,3.);
	RooFormulaVar sigma3("sigma3", "s*sigma_{1}", RooArgSet(scale,sigma1));

	//cout<<sigma3.getVal()<<endl;

	///Gaussian
	RooGaussian Gauss1("Gauss1", "Gauss1", Bplus_MM, mean1, sigma1);
	RooGaussian Gauss2("Gauss2", "Gauss2", Bplus_MM, mean1, sigma2);
	//RooGaussian Gauss3("Gauss3", "Gauss3", DeltaM, mean2, sigma3);
	//RooRealVar f_sig("f_{sig}", "signal fraction", 0.56, 0., 1.);

	/// Crystal Ball
	RooRealVar alpha("alpha","alpha",1.41,0.5,10.);
	RooRealVar n("n","n",34.4,0.,100.);
	RooCBShape crystal1("crystal1","CB PDF1",Bplus_MM,mean1,sigma1, alpha, n);
        RooCBShape crystal2("crystal2","CB PDF2",Bplus_MM,mean1,sigma3, alpha, n);	

	///signal pdf
	RooRealVar f_sig("f_{sig}", "signal fraction", 0.812,0.,0.95);
	//RooAddPdf signal("signal", "signal", RooArgList(Gauss1,Gauss2),RooArgList(f_sig));
	RooAddPdf signal("signal", "signal", RooArgList(crystal1,crystal2),RooArgList(f_sig));

	///Background model
	///-------------------------

	///Exponential
	RooRealVar exp_par("lambda","exp",-0.003535,-1.,0.);	
	RooExponential bkg_exp("bkg_exp","exponential background",Bplus_MM,exp_par);
	
	///Polynomial
	RooRealVar c_0("c_{0}","c_0",0.);
	//c_0.setConstant(0.);
	RooRealVar c_1("c_{1}","c_1",-.356,-1.,1.);
	RooRealVar c_2("c_{2}","c_2",-0.01,-1.,1.);

	RooChebychev bkg_Chebychev("bkg_Chebychev","bkg_Chebychev", Bplus_MM, RooArgList(c_1));
	RooPolynomial bkg_Poly("bkg_Poly","bkg_Poly",Bplus_MM,RooArgList(c_0));

	///Gaussian
	RooRealVar mean_bkg("mu_{bkg}", "mean_bkg", 5063.,5000.,5100.); 
	RooRealVar sigma_bkg("sigma_{bkg}", "sigma_bkg", 48.5, 40., 80.);
	RooRealVar sigma2_bkg("sigma_{2}", "sigma2", 19.2, 0.1, 100.5);
	RooGaussian bkg_Gauss("bkg_Gauss", "bkg_Gauss", Bplus_MM, mean_bkg, sigma_bkg);
	
	///bkg pdf
	RooRealVar f_bkg("f_{bkg}", "background fraction", 0.11, 0., 1.);
	RooRealVar f2_bkg("f2_{bkg}", "background fraction2", 0.11, 0., 1.);

	//RooAddPdf bkg("bkg", "bkg", RooArgList(bkg_exp, bkg_Chebychev), RooArgList(f_bkg));
	RooAddPdf bkg("bkg", "bkg", RooArgList(bkg_Gauss, bkg_Chebychev), RooArgList(f_bkg));

	///total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(Gauss1, bkg_Chebychev), RooArgList(n_sig,n_bkg));

	///Fit
	RooFitResult *result;
	if(binned) result = pdf->fitTo(*data_binned,Save(kTRUE),Extended());
	else result = pdf->fitTo(*data,Save(kTRUE),Extended(),NumCPU(3));
	cout << "result is --------------- "<<endl;
	result->Print(); 

	cout << "sigma2 =" << sigma1.getVal()*scale.getVal() << " pm " << sigma1.getError()*scale.getVal() << endl;

	///calculate # (signal)background events in signal region
	cout << endl;
	cout << endl;

	Bplus_MM.setRange("signal",mean1.getVal()-60.,mean1.getVal()+60.);
	
	RooAbsReal* S_fr= signal.createIntegral(Bplus_MM,NormSet(Bplus_MM),Range("signal"));
	Double_t S = S_fr->getVal()*n_sig.getVal();
	cout<<"S= " << S << endl;
	RooAbsReal* B_fr= bkg.createIntegral(Bplus_MM,NormSet(Bplus_MM),Range("signal"));
	Double_t B = B_fr->getVal()*n_bkg.getVal();
	cout<<"B= " << B << endl;
	
	cout<< "S/sqrt(S+B)= " << S/sqrt(S+B) << endl;
   	cout<<"S/B= " << S/B<< endl;

	cout << endl;
	cout << endl;

	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bplus_MM.frame();
	frame_m->SetTitle("");

	data->plotOn(frame_m,Name("data"),MarkerSize(0.1),Binning(100),DrawOption("e"),DataError(RooAbsData::SumW2));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack));
	pdf->plotOn(frame_m,Components(signal),LineColor(kBlue),LineStyle(kDashed));
	pdf->plotOn(frame_m,Components(bkg_Chebychev),LineColor(kRed),LineStyle(kDashed));
	pdf->paramOn(frame_m,Layout(0.6));
	frame_m->Draw();
	c1->Print("Preselection2/BmassFit.eps");

	RooPlot* frame_m2= Bplus_MM.frame();
	frame_m2->SetTitle("");
	//data->plotOn(frame_m2,MarkerSize(0.2),Binning(100),DrawOption("e"),DataError(RooAbsData::SumW2));
	data->plotOn(frame_m2,MarkerSize(0.5),Binning(50));
	pdf->plotOn(frame_m2,LineColor(kBlack),LineWidth(2));
	pdf->plotOn(frame_m2,Components(signal),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m2,Components(bkg_Chebychev),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	data->plotOn(frame_m2,MarkerSize(0.5),Binning(50));
	frame_m2->Draw();
	c1->Print("Preselection2/BmassFit2.eps");

	RooPlot* frame_m3= Bplus_MM.frame("");
	data->plotOn(frame_m3);
	pdf->plotOn(frame_m3);
	gPad->SetLogy();
	frame_m3->Draw();
	c1->Print("Preselection2/BmassFit_log.eps");


	gPad->SetLogy(0);
	///fit results
	double chi2 = frame_m->chiSquare("pdf","data",9);
	double covmatr = result->covQual();
	double edm = result->edm();
	cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl;

  	// Construct a histogram with the pulls of the data w.r.t the curve
	RooHist* hresid = frame_m->residHist("data","pdf") ;
  	RooHist* hpull = frame_m->pullHist("data","pdf") ;

	// Create a new frame to draw the residual distribution and add the distribution to the frame
	RooPlot* frame2 = Bplus_MM.frame(Title("Residual Distribution")) ;
	frame2->SetTitle("");
	frame2->addPlotable(hresid,"P") ;
	frame2->Draw();
	c1->Print("Preselection2/residual.eps");

	// Create a new frame to draw the pull distribution and add the distribution to the frame
	RooPlot* frame3 = Bplus_MM.frame(Title("Pull Distribution")) ;
	frame3->SetTitle("");	
	frame3->SetLabelFont(62,"Y");
	frame3->addPlotable(hpull,"P") ;
	frame3->Draw();
	c1->Print("Preselection2/pull.eps");

	if(sWeight){
		mean1.setConstant();
		sigma1.setConstant();
		sigma2.setConstant();
		f_sig.setConstant();
		f_bkg.setConstant();
		exp_par.setConstant();
		c_2.setConstant();
		alpha.setConstant();
		n.setConstant();
		//sigma3.setConstant();
		c_1.setConstant();
		scale.setConstant();
	
		SPlot* sData = new SPlot("sData","An SPlot",*data,pdf,RooArgList(n_sig,n_bkg)); 
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs = new TCanvas("Bd sWeight","Bd sWeight distribution");
		TH2 * SwDsHist = (TH2*)data->createHistogram("psiFit_Bplus_M,n_sig_sw");
		SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist->SetTitle("");
		//SwDs->Write();
		SwDsHist->Draw();
		SwDs->Print("Preselection2/Bd_sWeight.eps");

    		///Create output file
   		 TFile* output = new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_psiConstrained_sweight.root","RECREATE");

		 tree->SetBranchStatus("*",1);
   		 TTree* new_tree = tree->CloneTree();//CopyTree();    
    		 double w;
    		 TBranch* Bra_sw = new_tree->Branch("n_sig_sw", &w, "n_sig_sw/D");

  		  ///loop over events
    		  int numEvents = tree->GetEntries();
    
    		  for(int i=0; i< numEvents; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
			tree->GetEntry(i);
			w=sData->GetSWeight(i,"n_sig_sw");
			Bra_sw->Fill();
  		  }
   		 new_tree->Write();
   		 output->Close();
	}
}

void addVarsForBDT(){
	
    ///Load file
    TFile* file;
    file= new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_sweight.root");	
    TTree* tree = (TTree*) file->Get("DecayTree");	

    ///Reconstructed momenta
    double K_rec[4]; 
    double pip_rec[4]; 
    double pim_rec[4]; 
    double psi_rec[4];

    tree->SetBranchAddress("Kplus_PX",&K_rec[0]);
    tree->SetBranchAddress("Kplus_PY",&K_rec[1]);
    tree->SetBranchAddress("Kplus_PZ",&K_rec[2]); 
    tree->SetBranchAddress("Kplus_PE",&K_rec[3]); 

    tree->SetBranchAddress("piplus_PX",&pip_rec[0]);
    tree->SetBranchAddress("piplus_PY",&pip_rec[1]);
    tree->SetBranchAddress("piplus_PZ",&pip_rec[2]); 
    tree->SetBranchAddress("piplus_PE",&pip_rec[3]); 

    tree->SetBranchAddress("piminus_PX",&pim_rec[0]);
    tree->SetBranchAddress("piminus_PY",&pim_rec[1]);
    tree->SetBranchAddress("piminus_PZ",&pim_rec[2]); 
    tree->SetBranchAddress("piminus_PE",&pim_rec[3]); 
    
    tree->SetBranchAddress("J_psi_1S_PX",&psi_rec[0]);
    tree->SetBranchAddress("J_psi_1S_PY",&psi_rec[1]);
    tree->SetBranchAddress("J_psi_1S_PZ",&psi_rec[2]); 
    tree->SetBranchAddress("J_psi_1S_PE",&psi_rec[3]);

    ///Refitted momenta
    double K_dtf[4]; 
    double pip_dtf[4]; 
    double pim_dtf[4]; 
    double psi_dtf[4]; 

    tree->SetBranchAddress("psiDTF_Kplus_PX",&K_dtf[0]);
    tree->SetBranchAddress("psiDTF_Kplus_PY",&K_dtf[1]);
    tree->SetBranchAddress("psiDTF_Kplus_PZ",&K_dtf[2]); 
    tree->SetBranchAddress("psiDTF_Kplus_PE",&K_dtf[3]); 

    tree->SetBranchAddress("psiDTF_piplus_PX",&pip_dtf[0]);
    tree->SetBranchAddress("psiDTF_piplus_PY",&pip_dtf[1]);
    tree->SetBranchAddress("psiDTF_piplus_PZ",&pip_dtf[2]); 
    tree->SetBranchAddress("psiDTF_piplus_PE",&pip_dtf[3]); 

    tree->SetBranchAddress("psiDTF_piminus_PX",&pim_dtf[0]);
    tree->SetBranchAddress("psiDTF_piminus_PY",&pim_dtf[1]);
    tree->SetBranchAddress("psiDTF_piminus_PZ",&pim_dtf[2]); 
    tree->SetBranchAddress("psiDTF_piminus_PE",&pim_dtf[3]); 

    tree->SetBranchAddress("psiDTF_Jpsi_1S_PX",&psi_dtf[0]);
    tree->SetBranchAddress("psiDTF_Jpsi_1S_PY",&psi_dtf[1]);
    tree->SetBranchAddress("psiDTF_Jpsi_1S_PZ",&psi_dtf[2]); 
    tree->SetBranchAddress("psiDTF_Jpsi_1S_PE",&psi_dtf[3]); 
   		
    ///Create output file
    TFile* output = new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_sweight_bdt.root","RECREATE");

    TTree* new_tree = tree->CloneTree();//CopyTree();    
    int sample;
    double angK,angPip, angPim;
    double mKpi, mpipi, mKpipi, mJpsipi, mJpsipipi;
    double DTF_mKpi, DTF_mpipi, DTF_mKpipi, DTF_mJpsipi, DTF_mJpsipipi;

    TBranch* Bra_s = new_tree->Branch("sample", &sample, "sample/I");
    TBranch* Bra_angK = new_tree->Branch("angK", &angK, "angK/D");
    TBranch* Bra_angPip = new_tree->Branch("angPip", &angPip, "angPip/D");
    TBranch* Bra_angPim = new_tree->Branch("angPim", &angPim, "angPim/D");
    TBranch* Bra_mKpi = new_tree->Branch("mKpi", &mKpi, "mKpi/D");
    TBranch* Bra_mpipi = new_tree->Branch("mpipi", &mpipi, "mpipi/D");
    TBranch* Bra_mKpipi = new_tree->Branch("mKpipi", &mKpipi, "mKpipi/D");
    TBranch* Bra_mJpsipi = new_tree->Branch("mJpsipi", &mJpsipi, "mJpsipi/D");
    TBranch* Bra_mJpsipipi = new_tree->Branch("mJpsipipi", &mJpsipipi, "mJpsipipi/D");
    TBranch* Bra_DTF_mKpi = new_tree->Branch("DTF_mKpi", &DTF_mKpi, "DTF_mKpi/D");
    TBranch* Bra_DTF_mpipi = new_tree->Branch("DTF_mpipi", &DTF_mpipi, "DTF_mpipi/D");
    TBranch* Bra_DTF_mKpipi = new_tree->Branch("DTF_mKpipi", &DTF_mKpipi, "DTF_mKpipi/D");
    TBranch* Bra_DTF_mJpsipi = new_tree->Branch("DTF_mJpsipi", &DTF_mJpsipi, "DTF_mJpsipi/D");
    TBranch* Bra_DTF_mJpsipipi = new_tree->Branch("DTF_mJpsipipi", &DTF_mJpsipipi, "DTF_mJpsipipi/D");

    TRandom3 r;
    ///loop over events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
			tree->GetEntry(i);
			if(r.Rndm()<0.5)sample=1;
			else sample=2;
			Bra_s->Fill();
			
			TVector3 psi(psi_rec[0],psi_rec[1],0.);
			TVector3 K(K_rec[0],K_rec[1],0.);
			TVector3 pip(pip_rec[0],pip_rec[1],0.);
			TVector3 pim(pim_rec[0],pim_rec[1],0.);
			angK= psi.Angle(K);
			angPip= psi.Angle(pip);
			angPim= psi.Angle(pim);
			Bra_angK->Fill();
			Bra_angPip->Fill();
			Bra_angPim->Fill();

    		        TLorentzVector K_recP(K_rec[0],K_rec[1],K_rec[2],K_rec[3]);
        		TLorentzVector pip_recP(pip_rec[0],pip_rec[1],pip_rec[2],pip_rec[3]);
			TLorentzVector pim_recP(pim_rec[0],pim_rec[1],pim_rec[2],pim_rec[3]);
        		TLorentzVector psi_recP(psi_rec[0],psi_rec[1],psi_rec[2],psi_rec[3]);

			mKpi= (K_recP+pim_recP).M2();
			mpipi= (pip_recP+pim_recP).M2();
			mKpipi= (K_recP+pim_recP+pip_recP).M2();
			mJpsipi= (pip_recP+psi_recP).M2();
			mJpsipipi= (psi_recP+pim_recP+pip_recP).M2();

			Bra_mKpi->Fill();
			Bra_mpipi->Fill();
			Bra_mKpipi->Fill();
			Bra_mJpsipi->Fill();
			Bra_mJpsipipi->Fill();

        		TLorentzVector K_dtfP(K_dtf[0],K_dtf[1],K_dtf[2],K_dtf[3]);
       			TLorentzVector pip_dtfP(pip_dtf[0],pip_dtf[1],pip_dtf[2],pip_dtf[3]);
			TLorentzVector pim_dtfP(pim_dtf[0],pim_dtf[1],pim_dtf[2],pim_dtf[3]);
        		TLorentzVector psi_dtfP(psi_dtf[0],psi_dtf[1],psi_dtf[2],psi_dtf[3]);

			DTF_mKpi= (K_dtfP+pim_dtfP).M2();
			DTF_mpipi= (pip_dtfP+pim_dtfP).M2();
			DTF_mKpipi= (K_dtfP+pim_dtfP+pip_dtfP).M2();
			DTF_mJpsipi= (pip_dtfP+psi_dtfP).M2();
			DTF_mJpsipipi= (psi_dtfP+pim_dtfP+pip_dtfP).M2();

			Bra_DTF_mKpi->Fill();
			Bra_DTF_mpipi->Fill();
			Bra_DTF_mKpipi->Fill();
			Bra_DTF_mJpsipi->Fill();
			Bra_DTF_mJpsipipi->Fill();
			
			
     }
     new_tree->Write();
     output->Close();	
}


void applyBDTcut(){
   ///Load file
   TFile* file= new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_bdt.root");
   TTree* tree = (TTree*) file->Get("DecayTree");	
   tree->SetBranchStatus("*",0);
   tree->SetBranchStatus("Bplus_MM",1);
   tree->SetBranchStatus("*PX",1);
   tree->SetBranchStatus("*PY",1);
   tree->SetBranchStatus("*PZ",1);
   tree->SetBranchStatus("*PE",1);
   tree->SetBranchStatus("m*pi",1);
   tree->SetBranchStatus("DTF_m*pi",1);
   tree->SetBranchStatus("BDT*",1);
   tree->SetBranchStatus("sample",1);
   tree->SetBranchStatus("nCandidate",1);
   tree->SetBranchStatus("eventNumber",1);
   tree->SetBranchStatus("runNumber",1);

   TFile* output=new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_final.root","RECREATE");
   //TTree* new_tree = tree->CopyTree("(sample==1 && BDT_response2>-0.0786)||(sample==2 && BDT_response1>-0.0888)");
   TTree* new_tree = tree->CopyTree("BDT_response2>0.");
   new_tree->Write();
   
   file->Close();
   output->Close();	
}

void chooseMultCand(){
   ///Load file
   TFile* file= new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_final_multCand_PIDK.root");
   TTree* tree = (TTree*) file->Get("DecayTree");	
   ///Create output
   TFile* output=new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_final_chosenCand_PIDK_new.root","RECREATE");
   TTree* new_tree = tree->CopyTree("isSelectedMultipleCand==1 && Bplus_MM<5500");
   new_tree->Write();
   file->Close();
   output->Close();	
}

void fitBDT(){

	bool binned=false;
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
	file= new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_final.root");	
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bplus_MM",1);
	tree->SetBranchStatus("isSelectedMultipleCand",1);

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///B+
        RooRealVar Bplus_MM("Bplus_MM", "m(J/#psi K #pi #pi)", 5200., 5600.,"MeV");
	RooArgList list =  RooArgList(Bplus_MM);
	RooDataSet* data = new RooDataSet("data", "data", tree, list, "");
	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(Bplus_MM)));
	RooDataHist* data_binned = data_small->binnedClone();

	///Define fit model
	///----------------

	///Signal model
	///-----------------------

	RooRealVar mean1("mu", "mean1", 5281.71,5150.,5350.); 
	RooRealVar sigma1("sigma_{1}", "sigma1", 15.45,10.,50.);	
	RooRealVar sigma2("sigma_{2}", "sigma2", 29.,20.,40.);
	RooRealVar scale("s", "scale", 1.924,0.,10.);
	RooFormulaVar sigma3("sigma3", "s*sigma_{1}", RooArgSet(scale,sigma1));

	//cout<<sigma3.getVal()<<endl;

	///Gaussian
	RooGaussian Gauss1("Gauss1", "Gauss1", Bplus_MM, mean1, sigma1);
	RooGaussian Gauss2("Gauss2", "Gauss2", Bplus_MM, mean1, sigma2);
	//RooGaussian Gauss3("Gauss3", "Gauss3", DeltaM, mean2, sigma3);
	//RooRealVar f_sig("f_{sig}", "signal fraction", 0.56, 0., 1.);

	/// Crystal Ball
	RooRealVar alpha("alpha","alpha",1.41,0.5,10.);
	RooRealVar n("n","n",4.4,0.,100.);
	RooCBShape crystal1("crystal1","CB PDF1",Bplus_MM,mean1,sigma1, alpha, n);
        RooCBShape crystal2("crystal2","CB PDF2",Bplus_MM,mean1,sigma3, alpha, n);	

	///signal pdf
	RooRealVar f_sig("f_{sig}", "signal fraction", 0.812,0.,0.9);
	//RooAddPdf signal("signal", "signal", RooArgList(Gauss1,Gauss2),RooArgList(f_sig));
	RooAddPdf signal("signal", "signal", RooArgList(crystal1,crystal2),RooArgList(f_sig));

	///Background model
	///-------------------------

	///Exponential
	RooRealVar exp_par("lambda","exp",-0.003535,-1.,0.);	
	RooExponential bkg_exp("bkg_exp","exponential background",Bplus_MM,exp_par);
	
	///Polynomial
	RooRealVar c_0("c_{0}","c_0",0.);
	//c_0.setConstant(0.);
	RooRealVar c_1("c_{1}","c_1",-.356,-1.,0.);
	RooRealVar c_2("c_{2}","c_2",-0.01,-1.,1.);

	RooChebychev bkg_Chebychev("bkg_Chebychev","bkg_Chebychev", Bplus_MM, RooArgList(c_1));
	RooPolynomial bkg_Poly("bkg_Poly","bkg_Poly",Bplus_MM,RooArgList(c_0));

	///Gaussian
	RooRealVar mean_bkg("mu_{bkg}", "mean_bkg", 5063.,5000.,5100.); 
	RooRealVar sigma_bkg("sigma_{bkg}", "sigma_bkg", 48.5, 40., 80.);
	RooRealVar sigma2_bkg("sigma_{2}", "sigma2", 19.2, 0.1, 100.5);
	RooGaussian bkg_Gauss("bkg_Gauss", "bkg_Gauss", Bplus_MM, mean_bkg, sigma_bkg);
	
	///bkg pdf
	RooRealVar f_bkg("f_{bkg}", "background fraction", 0.11, 0., 1.);
	RooRealVar f2_bkg("f2_{bkg}", "background fraction2", 0.11, 0., 1.);

	//RooAddPdf bkg("bkg", "bkg", RooArgList(bkg_exp, bkg_Chebychev), RooArgList(f_bkg));
	RooAddPdf bkg("bkg", "bkg", RooArgList(bkg_Gauss, bkg_Chebychev), RooArgList(f_bkg));

	///total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(signal, bkg_Chebychev), RooArgList(n_sig,n_bkg));

	///Fit
	RooFitResult *result;
	if(binned) result = pdf->fitTo(*data_binned,Save(kTRUE),Extended());
	else result = pdf->fitTo(*data,Save(kTRUE),Extended(),NumCPU(3));
	cout << "result is --------------- "<<endl;
	result->Print(); 

	cout << "sigma2 =" << sigma1.getVal()*scale.getVal() << " pm " << sigma1.getError()*scale.getVal() << endl;

	///calculate # (signal)background events in signal region
	cout << endl;
	cout << endl;

	Bplus_MM.setRange("signal",mean1.getVal()-60.,mean1.getVal()+60.);
	
	RooAbsReal* S_fr= signal.createIntegral(Bplus_MM,NormSet(Bplus_MM),Range("signal"));
	Double_t S = S_fr->getVal()*n_sig.getVal();
	cout<<"S= " << S << endl;
	RooAbsReal* B_fr= bkg.createIntegral(Bplus_MM,NormSet(Bplus_MM),Range("signal"));
	Double_t B = B_fr->getVal()*n_bkg.getVal();
	cout<<"B= " << B << endl;
	
	cout<< "S/sqrt(S+B)= " << S/sqrt(S+B) << endl;
   	cout<<"S/B= " << S/B<< endl;

	cout << endl;
	cout << endl;

	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bplus_MM.frame();
	frame_m->SetTitle("");


	data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
	pdf->plotOn(frame_m,Components(signal),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	//pdf->plotOn(frame_m,Components(crystal1),LineColor(kBlue),LineStyle(kDashed));
	//pdf->plotOn(frame_m,Components(crystal2),LineColor(kBlue),LineStyle(kDashed));
	pdf->plotOn(frame_m,Components(bkg_Chebychev),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	pdf->paramOn(frame_m,Layout(0.6));
	data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
	frame_m->Draw();
	c1->Print("Final/BmassFit.eps");

	RooPlot* frame_m2= Bplus_MM.frame();
	frame_m2->SetTitle("");
	data->plotOn(frame_m2,MarkerSize(0.5),Binning(50));
	pdf->plotOn(frame_m2,LineColor(kBlack),LineWidth(2));
	pdf->plotOn(frame_m2,Components(signal),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m2,Components(bkg_Chebychev),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	data->plotOn(frame_m2,MarkerSize(0.5),Binning(50));
	frame_m2->Draw();
	c1->Print("Final/BmassFit2.eps");

	RooPlot* frame_m3= Bplus_MM.frame("");
	data->plotOn(frame_m3);
	pdf->plotOn(frame_m3);
	gPad->SetLogy();
	frame_m3->Draw();
	c1->Print("Final/BmassFit_log.eps");


	gPad->SetLogy(0);
	///fit results
	double chi2 = frame_m->chiSquare("pdf","data",9);
	double covmatr = result->covQual();
	double edm = result->edm();
	cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl;

  	// Construct a histogram with the pulls of the data w.r.t the curve
	RooHist* hresid = frame_m->residHist("data","pdf") ;
  	RooHist* hpull = frame_m->pullHist("data","pdf") ;

	// Create a new frame to draw the residual distribution and add the distribution to the frame
	RooPlot* frame2 = Bplus_MM.frame(Title("Residual Distribution")) ;
	frame2->SetTitle("");
	frame2->addPlotable(hresid,"P") ;
	frame2->Draw();
	c1->Print("Final/residual.eps");

	// Create a new frame to draw the pull distribution and add the distribution to the frame
	RooPlot* frame3 = Bplus_MM.frame(Title("Pull Distribution")) ;
	frame3->SetTitle("");	
	frame3->SetLabelFont(62,"Y");
	frame3->addPlotable(hpull,"P") ;
	frame3->Draw();
	c1->Print("Final/pull.eps");

	if(sWeight){
		mean1.setConstant();
		sigma1.setConstant();
		sigma2.setConstant();
		f_sig.setConstant();
		f_bkg.setConstant();
		exp_par.setConstant();
		c_0.setConstant();
		alpha.setConstant();
		n.setConstant();
		//sigma3.setConstant();
		c_1.setConstant();
		scale.setConstant();
	
		SPlot* sData = new SPlot("sData","An SPlot",*data,pdf,RooArgList(n_sig,n_bkg)); 
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs = new TCanvas("Bd sWeight","Bd sWeight distribution");
		TH2 * SwDsHist = (TH2*)data->createHistogram("Bplus_MM,n_sig_sw");
		SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist->SetTitle("");
		//SwDs->Write();
		SwDsHist->Draw();
		SwDs->Print("Final/Bd_sWeight.eps");

    		///Create output file
   		 TFile* output = new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_final_sweight.root","RECREATE");
		 tree->SetBranchStatus("*",1);
   		 TTree* new_tree = tree->CopyTree("");    
    		 double w;
    		 TBranch* Bra_sw = new_tree->Branch("n_sig_sw", &w, "n_sig_sw/D");

  		  ///loop over events
    		  int numEvents = tree->GetEntries();
    
    		  for(int i=0; i< numEvents; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
			tree->GetEntry(i);
			w=sData->GetSWeight(i,"n_sig_sw");
			Bra_sw->Fill();
  		  }
   		 new_tree->Write();
   		 output->Close();
	}
}

void makePlots(){
	gStyle->SetOptStat(0);
    	gStyle->SetTitleXOffset(0.9);
    	//gStyle->SetTitleYOffset(1.3);
    	gStyle->SetTitleXSize(0.05);
    	gStyle->SetTitleYSize(0.05);
    	gStyle->SetTitleFont(42,"X");
    	gStyle->SetTitleFont(42,"Y");
    	//gStyle->SetLabelSize(0.033,"X");
    	//gStyle->SetLabelSize(0.033,"Y");
    	gStyle->SetLabelFont(42,"X");
    	gStyle->SetLabelFont(42,"Y");
   	gStyle->SetLabelOffset(0.01,"X");
    	gStyle->SetPadTickX(1);
    	gStyle->SetPadTickY(1);
    
    	TH1::SetDefaultSumw2();
    	TH2::SetDefaultSumw2();
    
    ///Load file
    TFile* file;
    file= new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_sweight_bdt.root");	
    TTree* tree = (TTree*) file->Get("DecayTree");	

    double DTF_mJpsipipi, mJpsipipi;
    tree->SetBranchAddress("DTF_mJpsipipi",&DTF_mJpsipipi);
    tree->SetBranchAddress("mJpsipipi",&mJpsipipi);
    TH1D* hJpsipipi = new TH1D("hJpsipipi","",100,12,23);
    TH1D* hJpsipipi_DTF = new TH1D("hJpsipipi_DTF",";m^{2}(J/#psi #pi^{+} #pi^{-}) [GeV^{2}];Yield [norm.]",100,12,23);

    double DTF_mKpipi, mKpipi, DTF_mKpi, mKpi, mPsi;
    double mB, sw;
    tree->SetBranchAddress("DTF_mKpipi",&DTF_mKpipi);
    tree->SetBranchAddress("mKpipi",&mKpipi);
    tree->SetBranchAddress("DTF_mKpi",&DTF_mKpi);
    tree->SetBranchAddress("mKpi",&mKpi);
    tree->SetBranchAddress("Bplus_MM",&mB);
    tree->SetBranchAddress("n_sig_sw",&sw);
    tree->SetBranchAddress("J_psi_1S_MM",&mPsi);

    TH1D* hKpi = new TH1D("hKpi","",100,0.,4.5);
    TH1D* hKpi_DTF = new TH1D("hKpi_DTF",";m^{2}(K #pi) [GeV^{2}]; Yield [norm.]",100,0.,4.5);
    
    TH1D* hKpipi = new TH1D("hKpipi","",100,0.5,6);
    TH1D* hKpipi_DTF = new TH1D("hKpipi_DTF",";m^{2}(K #pi #pi) [GeV^{2}]; Yield [norm.]",100,0.5,6);
    TH1D* hKpipi_s = new TH1D("hKpipi_s",";m^{2}(K^{+} #pi^{+} #pi^{-}) [GeV^{2}]; Candidates ",100,0.5,6);
    TH1D* hKpipi_DTF_s = new TH1D("hKpipi_DTF_s",";m^{2}(K^{+} #pi^{+} #pi^{-}) [GeV^{2}]; Candidates",100,0.5,6);
    TH1D* hKpipi_b = new TH1D("hKpipi_b",";m^{2}(K #pi #pi) [GeV^{2}]; Yield [norm.]",100,0.5,6);
    TH1D* hKpipi_DTF_b = new TH1D("hKpipi_DTF_b",";m^{2}(K #pi #pi) [GeV^{2}]; Yield [norm.]",100,0.5,6);
    
    TH2D* mB_mKpipi = new TH2D("mB_mKpipi",";m^{2}(K^{+} #pi^{+} #pi^{-}) [GeV^{2}];m(B^{+}) [MeV]",50,0,6.5,50,5200,5600);
    TH1D* hpsi = new TH1D("mPsi","",100,3640,3740);

    ///loop over events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++){	
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
		tree->GetEntry(i);
		hJpsipipi->Fill(mJpsipipi/1000000,sw);
		hJpsipipi_DTF->Fill(DTF_mJpsipipi/1000000,sw);
		hKpipi->Fill(mKpipi/1000000,sw);
		hKpipi_DTF->Fill(DTF_mKpipi/1000000,sw);
		hKpi->Fill(mKpi/1000000,sw);
		hKpi_DTF->Fill(DTF_mKpi/1000000,sw);
		mB_mKpipi->Fill(mKpipi/1000000,mB);
		hpsi->Fill(mPsi,sw);
		if(abs(mB-5283.3)<60.){
			hKpipi_s->Fill(mKpipi/1000000);
			hKpipi_DTF_s->Fill(DTF_mKpipi/1000000);
		}
		if(mB>5400.){
			hKpipi_b->Fill(mKpipi/1000000);
			hKpipi_DTF_b->Fill(DTF_mKpipi/1000000);
		}
    }

    TCanvas* c= new TCanvas();

    hpsi->Draw("e");
    c->Print("plots/mpsi.eps");

    hKpi_DTF->SetLineColor(kRed);
    hKpi_DTF->Draw("");   
    hKpi->Draw("esame");
    c->Print("plots/mKpi.eps");

    hKpipi_DTF->SetLineColor(kRed);
    hKpipi_DTF->Draw("");   
    hKpipi->Draw("esame");
    c->Print("plots/mKpipi.eps");

    hKpipi_b->SetLineColor(kRed);
    hKpipi_s->Draw("");   
    hKpipi_b->DrawNormalized("esame",51495);
    c->Print("plots/mKpipi_sb.eps");

    hKpipi_DTF_b->SetLineColor(kRed);
    hKpipi_DTF_s->Draw("");   
    hKpipi_DTF_b->DrawNormalized("esame",51495);
    c->Print("plots/mKpipi_DTF_sb.eps");

    mB_mKpipi->Draw();
    c->Print("plots/mB_mKpipi.eps");
    
    hJpsipipi_DTF->SetLineColor(kBlue);
    hJpsipipi_DTF->DrawNormalized("hist",1);   
    hJpsipipi->SetLineColor(kRed);
    hJpsipipi->DrawNormalized("histsame",1);
    c->Print("plots/mJpsipipi.eps");
    gPad->SetLogy(1);
    c->Print("plots/mJpsipipi_log.eps");
}


 
int main(){
    time_t startTime = time(0);

    //preselect();
    //chooseBestPV("BFit","/auto/data/dargent/Bu2psiKpipi/data/data_preselected.root");
    //chooseBestPV("psiFit","/auto/data/dargent/Bu2psiKpipi/data/data_preselected_bestPV_BFit.root");
    //chooseBestPV("psiDTF","/auto/data/dargent/Bu2psiKpipi/data/data_preselected_bestPV_psiFit.root");
    //fitPreselected();  
    //fitPreselected_psiConstrained();      
    //addVarsForBDT();
    //applyBDTcut();
    //chooseMultCand();
    fitBDT();
    //makePlots();

    cout << "==============================================" << endl;
    cout << " Done " 
    << " \n Time since start " << (time(0) - startTime)/60.0
    << " min." << endl;
    cout << "==============================================" << endl;

return 0;
}

