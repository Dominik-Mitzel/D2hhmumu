#ifndef D2hhmumuReader_h
#define D2hhmumuReader_h
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>

// Header file for the classes stored in the TTree if any.

class D2hhmumuReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain


   D2hhmumuReader(TTree *tree=0);
   virtual ~D2hhmumuReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     InitMC();
   virtual void     activateRelevantBranches();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     createSubsample(TString name, double percentage); 
   virtual void     createMCtrainingSample(TString name);
   virtual void     createValidationSubsample(TString name); //used for sim09 validation purposes
   virtual void     fillHistograms(TString fname, bool isMC);
   virtual void     studyTriggerEfficiency();
   virtual double   slowpi_helicityAngle();
   virtual double   muon_helicityAngle();    
   virtual double   D0_helicityAngle();    


   double DTFdm();
   double dm();
   bool passGhostProbCut(double cut);
   bool isBkgSideband();
   bool isInMassRange();
   bool MCTruthmatched();
  
   virtual bool isL0Selected();
   virtual bool isHlt1Selected();
   //these functions have to be defined for each decay channel seperatley, here purely virtual defined
   virtual bool isHlt2Selected()=0;
   virtual void initializeMomenta()=0;
   
   TLorentzVector pDst;
   TLorentzVector pD  ;
   TLorentzVector pPis;
   TLorentzVector pMu0;
   TLorentzVector pMu1;
   TLorentzVector pH0;
   TLorentzVector pH1;

   TLorentzVector pDTFDst;
   TLorentzVector pDTFD ;
   TLorentzVector pDTFPis ;
   TLorentzVector pDTFMu0;
   TLorentzVector pDTFMu1;
   TLorentzVector pDTFH0;
   TLorentzVector pDTFH1;


   TH1* h1_DTFdm; 
   TH1* h1_dm; 
   TH1* h1_DTFmMuMu;
   TH1* h1_mMuMu;
   TH1* h1_DTFmD;
   TH1* h1_mD;
   TH2* h2_mD_deltaM;


// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxDst_ENDVERTEX_COV = 1;
   const Int_t kMaxDst_OWNPV_COV = 1;
   const Int_t kMaxDst_TOPPV_COV = 1;
   const Int_t kMaxD_ENDVERTEX_COV = 1;
   const Int_t kMaxD_OWNPV_COV = 1;
   const Int_t kMaxD_TOPPV_COV = 1;
   const Int_t kMaxD_ORIVX_COV = 1;
   const Int_t kMaxh0_OWNPV_COV = 1;
   const Int_t kMaxh0_TOPPV_COV = 1;
   const Int_t kMaxh0_ORIVX_COV = 1;
   const Int_t kMaxh1_OWNPV_COV = 1;
   const Int_t kMaxh1_TOPPV_COV = 1;
   const Int_t kMaxh1_ORIVX_COV = 1;
   const Int_t kMaxmu0_OWNPV_COV = 1;
   const Int_t kMaxmu0_TOPPV_COV = 1;
   const Int_t kMaxmu0_ORIVX_COV = 1;
   const Int_t kMaxmu1_OWNPV_COV = 1;
   const Int_t kMaxmu1_TOPPV_COV = 1;
   const Int_t kMaxmu1_ORIVX_COV = 1;
   const Int_t kMaxSlowpi_OWNPV_COV = 1;
   const Int_t kMaxSlowpi_TOPPV_COV = 1;
   const Int_t kMaxSlowpi_ORIVX_COV = 1;

   // Declaration of leaf types
   Double_t        Dst_MINIP;
   Double_t        Dst_MINIPCHI2;
   Double_t        Dst_MINIPNEXTBEST;
   Double_t        Dst_MINIPCHI2NEXTBEST;
   Double_t        Dst_ENDVERTEX_X;
   Double_t        Dst_ENDVERTEX_Y;
   Double_t        Dst_ENDVERTEX_Z;
   Double_t        Dst_ENDVERTEX_XERR;
   Double_t        Dst_ENDVERTEX_YERR;
   Double_t        Dst_ENDVERTEX_ZERR;
   Double_t        Dst_ENDVERTEX_CHI2;
   Int_t           Dst_ENDVERTEX_NDOF;
   Float_t         Dst_ENDVERTEX_COV_[3][3];
   Double_t        Dst_OWNPV_X;
   Double_t        Dst_OWNPV_Y;
   Double_t        Dst_OWNPV_Z;
   Double_t        Dst_OWNPV_XERR;
   Double_t        Dst_OWNPV_YERR;
   Double_t        Dst_OWNPV_ZERR;
   Double_t        Dst_OWNPV_CHI2;
   Int_t           Dst_OWNPV_NDOF;
   Float_t         Dst_OWNPV_COV_[3][3];
   Double_t        Dst_IP_OWNPV;
   Double_t        Dst_IPCHI2_OWNPV;
   Double_t        Dst_FD_OWNPV;
   Double_t        Dst_FDCHI2_OWNPV;
   Double_t        Dst_DIRA_OWNPV;
   Double_t        Dst_TOPPV_X;
   Double_t        Dst_TOPPV_Y;
   Double_t        Dst_TOPPV_Z;
   Double_t        Dst_TOPPV_XERR;
   Double_t        Dst_TOPPV_YERR;
   Double_t        Dst_TOPPV_ZERR;
   Double_t        Dst_TOPPV_CHI2;
   Int_t           Dst_TOPPV_NDOF;
   Float_t         Dst_TOPPV_COV_[3][3];
   Double_t        Dst_IP_TOPPV;
   Double_t        Dst_IPCHI2_TOPPV;
   Double_t        Dst_FD_TOPPV;
   Double_t        Dst_FDCHI2_TOPPV;
   Double_t        Dst_DIRA_TOPPV;
   Double_t        Dst_P;
   Double_t        Dst_PT;
   Double_t        Dst_PE;
   Double_t        Dst_PX;
   Double_t        Dst_PY;
   Double_t        Dst_PZ;
   Double_t        Dst_MM;
   Double_t        Dst_MMERR;
   Double_t        Dst_M;
   Int_t           Dst_ID;
   Double_t        Dst_TAU;
   Double_t        Dst_TAUERR;
   Double_t        Dst_TAUCHI2;
   Bool_t          Dst_L0Global_Dec;
   Bool_t          Dst_L0Global_TIS;
   Bool_t          Dst_L0Global_TOS;
   Bool_t          Dst_Hlt1Global_Dec;
   Bool_t          Dst_Hlt1Global_TIS;
   Bool_t          Dst_Hlt1Global_TOS;
   Bool_t          Dst_Hlt1Phys_Dec;
   Bool_t          Dst_Hlt1Phys_TIS;
   Bool_t          Dst_Hlt1Phys_TOS;
   Bool_t          Dst_Hlt2Global_Dec;
   Bool_t          Dst_Hlt2Global_TIS;
   Bool_t          Dst_Hlt2Global_TOS;
   Bool_t          Dst_Hlt2Phys_Dec;
   Bool_t          Dst_Hlt2Phys_TIS;
   Bool_t          Dst_Hlt2Phys_TOS;
   Double_t        Dst_DTF_CHI2;
   Double_t        Dst_DTF_D0_BPVIPCHI2;
   Double_t        Dst_DTF_D0_E;
   Double_t        Dst_DTF_D0_M;
   Double_t        Dst_DTF_D0_P;
   Double_t        Dst_DTF_D0_PT;
   Double_t        Dst_DTF_D0_PX;
   Double_t        Dst_DTF_D0_PY;
   Double_t        Dst_DTF_D0_PZ;
   Double_t        Dst_DTF_Dstarplus_E;
   Double_t        Dst_DTF_Dstarplus_M;
   Double_t        Dst_DTF_Dstarplus_P;
   Double_t        Dst_DTF_Dstarplus_PT;
   Double_t        Dst_DTF_Dstarplus_PX;
   Double_t        Dst_DTF_Dstarplus_PY;
   Double_t        Dst_DTF_Dstarplus_PZ;
   Double_t        Dst_DTF_NDOF;
   Double_t        Dst_DTF_Pis_BPVIPCHI2;
   Double_t        Dst_DTF_Pis_E;
   Double_t        Dst_DTF_Pis_M;
   Double_t        Dst_DTF_Pis_P;
   Double_t        Dst_DTF_Pis_PT;
   Double_t        Dst_DTF_Pis_PX;
   Double_t        Dst_DTF_Pis_PY;
   Double_t        Dst_DTF_Pis_PZ;
   Double_t        Dst_DTF_h0_E;
   Double_t        Dst_DTF_h0_P;
   Double_t        Dst_DTF_h0_PT;
   Double_t        Dst_DTF_h0_PX;
   Double_t        Dst_DTF_h0_PY;
   Double_t        Dst_DTF_h0_PZ;
   Double_t        Dst_DTF_h1_E;
   Double_t        Dst_DTF_h1_P;
   Double_t        Dst_DTF_h1_PT;
   Double_t        Dst_DTF_h1_PX;
   Double_t        Dst_DTF_h1_PY;
   Double_t        Dst_DTF_h1_PZ;
   Double_t        Dst_DTF_mu0_E;
   Double_t        Dst_DTF_mu0_P;
   Double_t        Dst_DTF_mu0_PT;
   Double_t        Dst_DTF_mu0_PX;
   Double_t        Dst_DTF_mu0_PY;
   Double_t        Dst_DTF_mu0_PZ;
   Double_t        Dst_DTF_mu1_E;
   Double_t        Dst_DTF_mu1_P;
   Double_t        Dst_DTF_mu1_PT;
   Double_t        Dst_DTF_mu1_PX;
   Double_t        Dst_DTF_mu1_PY;
   Double_t        Dst_DTF_mu1_PZ;
   Double_t        Dst_CONEANGLE_D;
   Double_t        Dst_CONEANGLE_Dstar;
   Double_t        Dst_CONEMULT_D;
   Double_t        Dst_CONEMULT_Dstar;
   Double_t        Dst_CONEPTASYM_D;
   Double_t        Dst_CONEPTASYM_Dstar;
   Double_t        Dst_MAXDOCA;
   Bool_t          Dst_L0HadronDecision_Dec;
   Bool_t          Dst_L0HadronDecision_TIS;
   Bool_t          Dst_L0HadronDecision_TOS;
   Bool_t          Dst_L0MuonDecision_Dec;
   Bool_t          Dst_L0MuonDecision_TIS;
   Bool_t          Dst_L0MuonDecision_TOS;
   Bool_t          Dst_L0DiMuonDecision_Dec;
   Bool_t          Dst_L0DiMuonDecision_TIS;
   Bool_t          Dst_L0DiMuonDecision_TOS;
   Bool_t          Dst_L0ElectronDecision_Dec;
   Bool_t          Dst_L0ElectronDecision_TIS;
   Bool_t          Dst_L0ElectronDecision_TOS;
   Bool_t          Dst_L0PhotonDecision_Dec;
   Bool_t          Dst_L0PhotonDecision_TIS;
   Bool_t          Dst_L0PhotonDecision_TOS;
   Bool_t          Dst_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          Dst_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          Dst_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          Dst_Hlt1DiMuonLowMassDecision_Dec;
   Bool_t          Dst_Hlt1DiMuonLowMassDecision_TIS;
   Bool_t          Dst_Hlt1DiMuonLowMassDecision_TOS;
   Bool_t          Dst_Hlt1SingleMuonNoIPDecision_Dec;
   Bool_t          Dst_Hlt1SingleMuonNoIPDecision_TIS;
   Bool_t          Dst_Hlt1SingleMuonNoIPDecision_TOS;
   Bool_t          Dst_Hlt1SingleMuonHighPTDecision_Dec;
   Bool_t          Dst_Hlt1SingleMuonHighPTDecision_TIS;
   Bool_t          Dst_Hlt1SingleMuonHighPTDecision_TOS;
   Bool_t          Dst_Hlt1TrackAllL0Decision_Dec;
   Bool_t          Dst_Hlt1TrackAllL0Decision_TIS;
   Bool_t          Dst_Hlt1TrackAllL0Decision_TOS;
   Bool_t          Dst_Hlt1TrackMuonDecision_Dec;
   Bool_t          Dst_Hlt1TrackMuonDecision_TIS;
   Bool_t          Dst_Hlt1TrackMuonDecision_TOS;
   Bool_t          Dst_Hlt1TrackPhotonDecision_Dec;
   Bool_t          Dst_Hlt1TrackPhotonDecision_TIS;
   Bool_t          Dst_Hlt1TrackPhotonDecision_TOS;
   Bool_t          Dst_Hlt1L0AnyDecision_Dec;
   Bool_t          Dst_Hlt1L0AnyDecision_TIS;
   Bool_t          Dst_Hlt1L0AnyDecision_TOS;
   Bool_t          Dst_Hlt1GlobalDecision_Dec;
   Bool_t          Dst_Hlt1GlobalDecision_TIS;
   Bool_t          Dst_Hlt1GlobalDecision_TOS;
   Bool_t          Dst_Hlt2SingleMuonDecision_Dec;
   Bool_t          Dst_Hlt2SingleMuonDecision_TIS;
   Bool_t          Dst_Hlt2SingleMuonDecision_TOS;
   Bool_t          Dst_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          Dst_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          Dst_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          Dst_Hlt2CharmSemilepD2HMuMuDecision_Dec;
   Bool_t          Dst_Hlt2CharmSemilepD2HMuMuDecision_TIS;
   Bool_t          Dst_Hlt2CharmSemilepD2HMuMuDecision_TOS;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;
   Bool_t          Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;
   Bool_t          Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;
   Bool_t          Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;
   Bool_t          Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;
   Bool_t          Dst_Hlt2CharmSemilepD02KKMuMuDecision_Dec;
   Bool_t          Dst_Hlt2CharmSemilepD02KKMuMuDecision_TIS;
   Bool_t          Dst_Hlt2CharmSemilepD02KKMuMuDecision_TOS;
   Bool_t          Dst_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;
   Bool_t          Dst_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;
   Bool_t          Dst_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_K3piDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_K3piDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_K3piDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_4piDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_4piDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_4piDecision_TOS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;
   Bool_t          Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;
   Double_t        Dst_cpx_1_00;
   Double_t        Dst_cpy_1_00;
   Double_t        Dst_cpz_1_00;
   Double_t        Dst_cpt_1_00;
   Double_t        Dst_cp_1_00;
   Int_t           Dst_cmult_1_00;
   Double_t        Dst_deltaEta_1_00;
   Double_t        Dst_deltaPhi_1_00;
   Double_t        Dst_pxasy_1_00;
   Double_t        Dst_pyasy_1_00;
   Double_t        Dst_pzasy_1_00;
   Double_t        Dst_pasy_1_00;
   Double_t        Dst_ptasy_1_00;
   Double_t        Dst_cpx_1_10;
   Double_t        Dst_cpy_1_10;
   Double_t        Dst_cpz_1_10;
   Double_t        Dst_cpt_1_10;
   Double_t        Dst_cp_1_10;
   Int_t           Dst_cmult_1_10;
   Double_t        Dst_deltaEta_1_10;
   Double_t        Dst_deltaPhi_1_10;
   Double_t        Dst_pxasy_1_10;
   Double_t        Dst_pyasy_1_10;
   Double_t        Dst_pzasy_1_10;
   Double_t        Dst_pasy_1_10;
   Double_t        Dst_ptasy_1_10;
   Double_t        Dst_cpx_1_20;
   Double_t        Dst_cpy_1_20;
   Double_t        Dst_cpz_1_20;
   Double_t        Dst_cpt_1_20;
   Double_t        Dst_cp_1_20;
   Int_t           Dst_cmult_1_20;
   Double_t        Dst_deltaEta_1_20;
   Double_t        Dst_deltaPhi_1_20;
   Double_t        Dst_pxasy_1_20;
   Double_t        Dst_pyasy_1_20;
   Double_t        Dst_pzasy_1_20;
   Double_t        Dst_pasy_1_20;
   Double_t        Dst_ptasy_1_20;
   Double_t        Dst_cpx_1_30;
   Double_t        Dst_cpy_1_30;
   Double_t        Dst_cpz_1_30;
   Double_t        Dst_cpt_1_30;
   Double_t        Dst_cp_1_30;
   Int_t           Dst_cmult_1_30;
   Double_t        Dst_deltaEta_1_30;
   Double_t        Dst_deltaPhi_1_30;
   Double_t        Dst_pxasy_1_30;
   Double_t        Dst_pyasy_1_30;
   Double_t        Dst_pzasy_1_30;
   Double_t        Dst_pasy_1_30;
   Double_t        Dst_ptasy_1_30;
   Double_t        Dst_cpx_1_40;
   Double_t        Dst_cpy_1_40;
   Double_t        Dst_cpz_1_40;
   Double_t        Dst_cpt_1_40;
   Double_t        Dst_cp_1_40;
   Int_t           Dst_cmult_1_40;
   Double_t        Dst_deltaEta_1_40;
   Double_t        Dst_deltaPhi_1_40;
   Double_t        Dst_pxasy_1_40;
   Double_t        Dst_pyasy_1_40;
   Double_t        Dst_pzasy_1_40;
   Double_t        Dst_pasy_1_40;
   Double_t        Dst_ptasy_1_40;
   Double_t        Dst_cpx_1_50;
   Double_t        Dst_cpy_1_50;
   Double_t        Dst_cpz_1_50;
   Double_t        Dst_cpt_1_50;
   Double_t        Dst_cp_1_50;
   Int_t           Dst_cmult_1_50;
   Double_t        Dst_deltaEta_1_50;
   Double_t        Dst_deltaPhi_1_50;
   Double_t        Dst_pxasy_1_50;
   Double_t        Dst_pyasy_1_50;
   Double_t        Dst_pzasy_1_50;
   Double_t        Dst_pasy_1_50;
   Double_t        Dst_ptasy_1_50;
   Double_t        Dst_cpx_1_60;
   Double_t        Dst_cpy_1_60;
   Double_t        Dst_cpz_1_60;
   Double_t        Dst_cpt_1_60;
   Double_t        Dst_cp_1_60;
   Int_t           Dst_cmult_1_60;
   Double_t        Dst_deltaEta_1_60;
   Double_t        Dst_deltaPhi_1_60;
   Double_t        Dst_pxasy_1_60;
   Double_t        Dst_pyasy_1_60;
   Double_t        Dst_pzasy_1_60;
   Double_t        Dst_pasy_1_60;
   Double_t        Dst_ptasy_1_60;
   Double_t        Dst_cpx_1_70;
   Double_t        Dst_cpy_1_70;
   Double_t        Dst_cpz_1_70;
   Double_t        Dst_cpt_1_70;
   Double_t        Dst_cp_1_70;
   Int_t           Dst_cmult_1_70;
   Double_t        Dst_deltaEta_1_70;
   Double_t        Dst_deltaPhi_1_70;
   Double_t        Dst_pxasy_1_70;
   Double_t        Dst_pyasy_1_70;
   Double_t        Dst_pzasy_1_70;
   Double_t        Dst_pasy_1_70;
   Double_t        Dst_ptasy_1_70;
   Double_t        D_MINIP;
   Double_t        D_MINIPCHI2;
   Double_t        D_MINIPNEXTBEST;
   Double_t        D_MINIPCHI2NEXTBEST;
   Double_t        D_ENDVERTEX_X;
   Double_t        D_ENDVERTEX_Y;
   Double_t        D_ENDVERTEX_Z;
   Double_t        D_ENDVERTEX_XERR;
   Double_t        D_ENDVERTEX_YERR;
   Double_t        D_ENDVERTEX_ZERR;
   Double_t        D_ENDVERTEX_CHI2;
   Int_t           D_ENDVERTEX_NDOF;
   Float_t         D_ENDVERTEX_COV_[3][3];
   Double_t        D_OWNPV_X;
   Double_t        D_OWNPV_Y;
   Double_t        D_OWNPV_Z;
   Double_t        D_OWNPV_XERR;
   Double_t        D_OWNPV_YERR;
   Double_t        D_OWNPV_ZERR;
   Double_t        D_OWNPV_CHI2;
   Int_t           D_OWNPV_NDOF;
   Float_t         D_OWNPV_COV_[3][3];
   Double_t        D_IP_OWNPV;
   Double_t        D_IPCHI2_OWNPV;
   Double_t        D_FD_OWNPV;
   Double_t        D_FDCHI2_OWNPV;
   Double_t        D_DIRA_OWNPV;
   Double_t        D_TOPPV_X;
   Double_t        D_TOPPV_Y;
   Double_t        D_TOPPV_Z;
   Double_t        D_TOPPV_XERR;
   Double_t        D_TOPPV_YERR;
   Double_t        D_TOPPV_ZERR;
   Double_t        D_TOPPV_CHI2;
   Int_t           D_TOPPV_NDOF;
   Float_t         D_TOPPV_COV_[3][3];
   Double_t        D_IP_TOPPV;
   Double_t        D_IPCHI2_TOPPV;
   Double_t        D_FD_TOPPV;
   Double_t        D_FDCHI2_TOPPV;
   Double_t        D_DIRA_TOPPV;
   Double_t        D_ORIVX_X;
   Double_t        D_ORIVX_Y;
   Double_t        D_ORIVX_Z;
   Double_t        D_ORIVX_XERR;
   Double_t        D_ORIVX_YERR;
   Double_t        D_ORIVX_ZERR;
   Double_t        D_ORIVX_CHI2;
   Int_t           D_ORIVX_NDOF;
   Float_t         D_ORIVX_COV_[3][3];
   Double_t        D_IP_ORIVX;
   Double_t        D_IPCHI2_ORIVX;
   Double_t        D_FD_ORIVX;
   Double_t        D_FDCHI2_ORIVX;
   Double_t        D_DIRA_ORIVX;
   Double_t        D_P;
   Double_t        D_PT;
   Double_t        D_PE;
   Double_t        D_PX;
   Double_t        D_PY;
   Double_t        D_PZ;
   Double_t        D_MM;
   Double_t        D_MMERR;
   Double_t        D_M;
   Int_t           D_ID;
   Double_t        D_TAU;
   Double_t        D_TAUERR;
   Double_t        D_TAUCHI2;
   Bool_t          D_L0Global_Dec;
   Bool_t          D_L0Global_TIS;
   Bool_t          D_L0Global_TOS;
   Bool_t          D_Hlt1Global_Dec;
   Bool_t          D_Hlt1Global_TIS;
   Bool_t          D_Hlt1Global_TOS;
   Bool_t          D_Hlt1Phys_Dec;
   Bool_t          D_Hlt1Phys_TIS;
   Bool_t          D_Hlt1Phys_TOS;
   Bool_t          D_Hlt2Global_Dec;
   Bool_t          D_Hlt2Global_TIS;
   Bool_t          D_Hlt2Global_TOS;
   Bool_t          D_Hlt2Phys_Dec;
   Bool_t          D_Hlt2Phys_TIS;
   Bool_t          D_Hlt2Phys_TOS;
   Double_t        D_DTF_CHI2;
   Double_t        D_DTF_NDOF;
   Double_t        D_D_DTF_M;
   Double_t        D_D_DTF_M_False;
   Double_t        D_DiMuon_Mass;
   Double_t        D_LV01;
   Double_t        D_LV02;
   Double_t        D_LV03;
   Double_t        D_LV04;
   Double_t        D_MAXDOCA;
   Bool_t          D_L0HadronDecision_Dec;
   Bool_t          D_L0HadronDecision_TIS;
   Bool_t          D_L0HadronDecision_TOS;
   Bool_t          D_L0MuonDecision_Dec;
   Bool_t          D_L0MuonDecision_TIS;
   Bool_t          D_L0MuonDecision_TOS;
   Bool_t          D_L0DiMuonDecision_Dec;
   Bool_t          D_L0DiMuonDecision_TIS;
   Bool_t          D_L0DiMuonDecision_TOS;
   Bool_t          D_L0ElectronDecision_Dec;
   Bool_t          D_L0ElectronDecision_TIS;
   Bool_t          D_L0ElectronDecision_TOS;
   Bool_t          D_L0PhotonDecision_Dec;
   Bool_t          D_L0PhotonDecision_TIS;
   Bool_t          D_L0PhotonDecision_TOS;
   Bool_t          D_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          D_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          D_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          D_Hlt1DiMuonLowMassDecision_Dec;
   Bool_t          D_Hlt1DiMuonLowMassDecision_TIS;
   Bool_t          D_Hlt1DiMuonLowMassDecision_TOS;
   Bool_t          D_Hlt1SingleMuonNoIPDecision_Dec;
   Bool_t          D_Hlt1SingleMuonNoIPDecision_TIS;
   Bool_t          D_Hlt1SingleMuonNoIPDecision_TOS;
   Bool_t          D_Hlt1SingleMuonHighPTDecision_Dec;
   Bool_t          D_Hlt1SingleMuonHighPTDecision_TIS;
   Bool_t          D_Hlt1SingleMuonHighPTDecision_TOS;
   Bool_t          D_Hlt1TrackAllL0Decision_Dec;
   Bool_t          D_Hlt1TrackAllL0Decision_TIS;
   Bool_t          D_Hlt1TrackAllL0Decision_TOS;
   Bool_t          D_Hlt1TrackMuonDecision_Dec;
   Bool_t          D_Hlt1TrackMuonDecision_TIS;
   Bool_t          D_Hlt1TrackMuonDecision_TOS;
   Bool_t          D_Hlt1TrackPhotonDecision_Dec;
   Bool_t          D_Hlt1TrackPhotonDecision_TIS;
   Bool_t          D_Hlt1TrackPhotonDecision_TOS;
   Bool_t          D_Hlt1L0AnyDecision_Dec;
   Bool_t          D_Hlt1L0AnyDecision_TIS;
   Bool_t          D_Hlt1L0AnyDecision_TOS;
   Bool_t          D_Hlt1GlobalDecision_Dec;
   Bool_t          D_Hlt1GlobalDecision_TIS;
   Bool_t          D_Hlt1GlobalDecision_TOS;
   Bool_t          D_Hlt2SingleMuonDecision_Dec;
   Bool_t          D_Hlt2SingleMuonDecision_TIS;
   Bool_t          D_Hlt2SingleMuonDecision_TOS;
   Bool_t          D_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          D_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          D_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          D_Hlt2CharmSemilepD2HMuMuDecision_Dec;
   Bool_t          D_Hlt2CharmSemilepD2HMuMuDecision_TIS;
   Bool_t          D_Hlt2CharmSemilepD2HMuMuDecision_TOS;
   Bool_t          D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;
   Bool_t          D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;
   Bool_t          D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;
   Bool_t          D_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;
   Bool_t          D_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;
   Bool_t          D_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;
   Bool_t          D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;
   Bool_t          D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;
   Bool_t          D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;
   Bool_t          D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;
   Bool_t          D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;
   Bool_t          D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;
   Bool_t          D_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;
   Bool_t          D_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;
   Bool_t          D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;
   Bool_t          D_Hlt2CharmSemilepD02KKMuMuDecision_Dec;
   Bool_t          D_Hlt2CharmSemilepD02KKMuMuDecision_TIS;
   Bool_t          D_Hlt2CharmSemilepD02KKMuMuDecision_TOS;
   Bool_t          D_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;
   Bool_t          D_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;
   Bool_t          D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHHH_K3piDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHHH_K3piDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHHH_K3piDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHHH_4piDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHHH_4piDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHHH_4piDecision_TOS;
   Bool_t          D_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;
   Bool_t          D_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;
   Bool_t          D_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;
   Double_t        D_cpx_1_00;
   Double_t        D_cpy_1_00;
   Double_t        D_cpz_1_00;
   Double_t        D_cpt_1_00;
   Double_t        D_cp_1_00;
   Int_t           D_cmult_1_00;
   Double_t        D_deltaEta_1_00;
   Double_t        D_deltaPhi_1_00;
   Double_t        D_pxasy_1_00;
   Double_t        D_pyasy_1_00;
   Double_t        D_pzasy_1_00;
   Double_t        D_pasy_1_00;
   Double_t        D_ptasy_1_00;
   Double_t        D_cpx_1_10;
   Double_t        D_cpy_1_10;
   Double_t        D_cpz_1_10;
   Double_t        D_cpt_1_10;
   Double_t        D_cp_1_10;
   Int_t           D_cmult_1_10;
   Double_t        D_deltaEta_1_10;
   Double_t        D_deltaPhi_1_10;
   Double_t        D_pxasy_1_10;
   Double_t        D_pyasy_1_10;
   Double_t        D_pzasy_1_10;
   Double_t        D_pasy_1_10;
   Double_t        D_ptasy_1_10;
   Double_t        D_cpx_1_20;
   Double_t        D_cpy_1_20;
   Double_t        D_cpz_1_20;
   Double_t        D_cpt_1_20;
   Double_t        D_cp_1_20;
   Int_t           D_cmult_1_20;
   Double_t        D_deltaEta_1_20;
   Double_t        D_deltaPhi_1_20;
   Double_t        D_pxasy_1_20;
   Double_t        D_pyasy_1_20;
   Double_t        D_pzasy_1_20;
   Double_t        D_pasy_1_20;
   Double_t        D_ptasy_1_20;
   Double_t        D_cpx_1_30;
   Double_t        D_cpy_1_30;
   Double_t        D_cpz_1_30;
   Double_t        D_cpt_1_30;
   Double_t        D_cp_1_30;
   Int_t           D_cmult_1_30;
   Double_t        D_deltaEta_1_30;
   Double_t        D_deltaPhi_1_30;
   Double_t        D_pxasy_1_30;
   Double_t        D_pyasy_1_30;
   Double_t        D_pzasy_1_30;
   Double_t        D_pasy_1_30;
   Double_t        D_ptasy_1_30;
   Double_t        D_cpx_1_40;
   Double_t        D_cpy_1_40;
   Double_t        D_cpz_1_40;
   Double_t        D_cpt_1_40;
   Double_t        D_cp_1_40;
   Int_t           D_cmult_1_40;
   Double_t        D_deltaEta_1_40;
   Double_t        D_deltaPhi_1_40;
   Double_t        D_pxasy_1_40;
   Double_t        D_pyasy_1_40;
   Double_t        D_pzasy_1_40;
   Double_t        D_pasy_1_40;
   Double_t        D_ptasy_1_40;
   Double_t        D_cpx_1_50;
   Double_t        D_cpy_1_50;
   Double_t        D_cpz_1_50;
   Double_t        D_cpt_1_50;
   Double_t        D_cp_1_50;
   Int_t           D_cmult_1_50;
   Double_t        D_deltaEta_1_50;
   Double_t        D_deltaPhi_1_50;
   Double_t        D_pxasy_1_50;
   Double_t        D_pyasy_1_50;
   Double_t        D_pzasy_1_50;
   Double_t        D_pasy_1_50;
   Double_t        D_ptasy_1_50;
   Double_t        D_cpx_1_60;
   Double_t        D_cpy_1_60;
   Double_t        D_cpz_1_60;
   Double_t        D_cpt_1_60;
   Double_t        D_cp_1_60;
   Int_t           D_cmult_1_60;
   Double_t        D_deltaEta_1_60;
   Double_t        D_deltaPhi_1_60;
   Double_t        D_pxasy_1_60;
   Double_t        D_pyasy_1_60;
   Double_t        D_pzasy_1_60;
   Double_t        D_pasy_1_60;
   Double_t        D_ptasy_1_60;
   Double_t        D_cpx_1_70;
   Double_t        D_cpy_1_70;
   Double_t        D_cpz_1_70;
   Double_t        D_cpt_1_70;
   Double_t        D_cp_1_70;
   Int_t           D_cmult_1_70;
   Double_t        D_deltaEta_1_70;
   Double_t        D_deltaPhi_1_70;
   Double_t        D_pxasy_1_70;
   Double_t        D_pyasy_1_70;
   Double_t        D_pzasy_1_70;
   Double_t        D_pasy_1_70;
   Double_t        D_ptasy_1_70;
   Double_t        h0_MINIP;
   Double_t        h0_MINIPCHI2;
   Double_t        h0_MINIPNEXTBEST;
   Double_t        h0_MINIPCHI2NEXTBEST;
   Double_t        h0_OWNPV_X;
   Double_t        h0_OWNPV_Y;
   Double_t        h0_OWNPV_Z;
   Double_t        h0_OWNPV_XERR;
   Double_t        h0_OWNPV_YERR;
   Double_t        h0_OWNPV_ZERR;
   Double_t        h0_OWNPV_CHI2;
   Int_t           h0_OWNPV_NDOF;
   Float_t         h0_OWNPV_COV_[3][3];
   Double_t        h0_IP_OWNPV;
   Double_t        h0_IPCHI2_OWNPV;
   Double_t        h0_TOPPV_X;
   Double_t        h0_TOPPV_Y;
   Double_t        h0_TOPPV_Z;
   Double_t        h0_TOPPV_XERR;
   Double_t        h0_TOPPV_YERR;
   Double_t        h0_TOPPV_ZERR;
   Double_t        h0_TOPPV_CHI2;
   Int_t           h0_TOPPV_NDOF;
   Float_t         h0_TOPPV_COV_[3][3];
   Double_t        h0_IP_TOPPV;
   Double_t        h0_IPCHI2_TOPPV;
   Double_t        h0_ORIVX_X;
   Double_t        h0_ORIVX_Y;
   Double_t        h0_ORIVX_Z;
   Double_t        h0_ORIVX_XERR;
   Double_t        h0_ORIVX_YERR;
   Double_t        h0_ORIVX_ZERR;
   Double_t        h0_ORIVX_CHI2;
   Int_t           h0_ORIVX_NDOF;
   Float_t         h0_ORIVX_COV_[3][3];
   Double_t        h0_IP_ORIVX;
   Double_t        h0_IPCHI2_ORIVX;
   Double_t        h0_P;
   Double_t        h0_PT;
   Double_t        h0_PE;
   Double_t        h0_PX;
   Double_t        h0_PY;
   Double_t        h0_PZ;
   Double_t        h0_M;
   Double_t        h0_L0Calo_HCAL_realET;
   Double_t        h0_L0Calo_HCAL_xProjection;
   Double_t        h0_L0Calo_HCAL_yProjection;
   Int_t           h0_L0Calo_HCAL_region;
   Double_t        h0_L0Calo_HCAL_TriggerET;
   Double_t        h0_L0Calo_HCAL_TriggerHCALET;
   Double_t        h0_L0Calo_HCAL_xTrigger;
   Double_t        h0_L0Calo_HCAL_yTrigger;
   Int_t           h0_ID;
   Double_t        h0_CombDLLMu;
   Double_t        h0_ProbNNmu;
   Double_t        h0_ProbNNghost;
   Double_t        h0_InMuonAcc;
   Double_t        h0_MuonDist2;
   Int_t           h0_regionInM2;
   Bool_t          h0_hasMuon;
   Bool_t          h0_isMuon;
   Bool_t          h0_isMuonLoose;
   Int_t           h0_NShared;
   Double_t        h0_MuonLLmu;
   Double_t        h0_MuonLLbg;
   Int_t           h0_isMuonFromProto;
   Double_t        h0_PIDe;
   Double_t        h0_PIDmu;
   Double_t        h0_PIDK;
   Double_t        h0_PIDp;
   Double_t        h0_ProbNNe;
   Double_t        h0_ProbNNk;
   Double_t        h0_ProbNNp;
   Double_t        h0_ProbNNpi;
   Bool_t          h0_hasRich;
   Bool_t          h0_hasCalo;
   Bool_t          h0_UsedRichAerogel;
   Bool_t          h0_UsedRich1Gas;
   Bool_t          h0_UsedRich2Gas;
   Bool_t          h0_RichAboveElThres;
   Bool_t          h0_RichAboveMuThres;
   Bool_t          h0_RichAbovePiThres;
   Bool_t          h0_RichAboveKaThres;
   Bool_t          h0_RichAbovePrThres;
   Double_t        h0_RichDLLe;
   Double_t        h0_RichDLLmu;
   Double_t        h0_RichDLLpi;
   Double_t        h0_RichDLLk;
   Double_t        h0_RichDLLp;
   Double_t        h0_RichDLLbt;
   Bool_t          h0_InAccMuon;
   Double_t        h0_MuonMuLL;
   Double_t        h0_MuonBkgLL;
   Int_t           h0_MuonNShared;
   Bool_t          h0_InAccEcal;
   Double_t        h0_CaloEcalE;
   Double_t        h0_EcalPIDe;
   Double_t        h0_EcalPIDmu;
   Bool_t          h0_InAccHcal;
   Double_t        h0_CaloHcalE;
   Double_t        h0_HcalPIDe;
   Double_t        h0_HcalPIDmu;
   Bool_t          h0_InAccPrs;
   Double_t        h0_PrsPIDe;
   Double_t        h0_CaloPrsE;
   Bool_t          h0_InAccSpd;
   Double_t        h0_CaloSpdE;
   Bool_t          h0_InAccBrem;
   Double_t        h0_BremPIDe;
   Double_t        h0_VeloCharge;
   Double_t        h0_RICHDLLe;
   Double_t        h0_RICHDLLmu;
   Double_t        h0_RICHDLLpi;
   Double_t        h0_RICHDLLK;
   Double_t        h0_RICHDLLp;
   Double_t        h0_RICHDLLbt;
   Int_t           h0_RICHBestID;
   Int_t           h0_RICHThreshold;
   Int_t           h0_RICHThresholdEl;
   Int_t           h0_RICHThresholdMu;
   Int_t           h0_RICHThresholdPi;
   Int_t           h0_RICHThresholdKa;
   Int_t           h0_RICHThresholdPr;
   Int_t           h0_RICHAerogelUsed;
   Int_t           h0_RICH1GasUsed;
   Int_t           h0_RICH2GasUsed;
   Double_t        h0_TRACK_Eta;
   Double_t        h0_TRACK_Phi;
   Double_t        h0_Aerogel_X;
   Double_t        h0_Aerogel_Y;
   Double_t        h0_Aerogel_Z;
   Double_t        h0_Aerogel_Rho;
   Double_t        h0_Aerogel_Phi;
   Double_t        h0_Rich1Gas_X;
   Double_t        h0_Rich1Gas_Y;
   Double_t        h0_Rich1Gas_Z;
   Double_t        h0_Rich1Gas_Rho;
   Double_t        h0_Rich1Gas_Phi;
   Double_t        h0_Rich2Gas_X;
   Double_t        h0_Rich2Gas_Y;
   Double_t        h0_Rich2Gas_Z;
   Double_t        h0_Rich2Gas_Rho;
   Double_t        h0_Rich2Gas_Phi;
   Bool_t          h0_L0Global_Dec;
   Bool_t          h0_L0Global_TIS;
   Bool_t          h0_L0Global_TOS;
   Bool_t          h0_Hlt1Global_Dec;
   Bool_t          h0_Hlt1Global_TIS;
   Bool_t          h0_Hlt1Global_TOS;
   Bool_t          h0_Hlt1Phys_Dec;
   Bool_t          h0_Hlt1Phys_TIS;
   Bool_t          h0_Hlt1Phys_TOS;
   Bool_t          h0_Hlt2Global_Dec;
   Bool_t          h0_Hlt2Global_TIS;
   Bool_t          h0_Hlt2Global_TOS;
   Bool_t          h0_Hlt2Phys_Dec;
   Bool_t          h0_Hlt2Phys_TIS;
   Bool_t          h0_Hlt2Phys_TOS;
   Int_t           h0_TRACK_Type;
   Int_t           h0_TRACK_Key;
   Double_t        h0_TRACK_CHI2;
   Int_t           h0_TRACK_NDOF;
   Double_t        h0_TRACK_CHI2NDOF;
   Double_t        h0_TRACK_PCHI2;
   Double_t        h0_TRACK_VeloCHI2NDOF;
   Double_t        h0_TRACK_TCHI2NDOF;
   Double_t        h0_VELO_UTID;
   Double_t        h0_TRACK_FirstMeasurementX;
   Double_t        h0_TRACK_FirstMeasurementY;
   Double_t        h0_TRACK_FirstMeasurementZ;
   Double_t        h0_TRACK_MatchCHI2;
   Double_t        h0_TRACK_GhostProb;
   Double_t        h0_TRACK_CloneDist;
   Double_t        h0_TRACK_Likelihood;
   Bool_t          h0_L0HadronDecision_Dec;
   Bool_t          h0_L0HadronDecision_TIS;
   Bool_t          h0_L0HadronDecision_TOS;
   Bool_t          h0_L0MuonDecision_Dec;
   Bool_t          h0_L0MuonDecision_TIS;
   Bool_t          h0_L0MuonDecision_TOS;
   Bool_t          h0_L0DiMuonDecision_Dec;
   Bool_t          h0_L0DiMuonDecision_TIS;
   Bool_t          h0_L0DiMuonDecision_TOS;
   Bool_t          h0_L0ElectronDecision_Dec;
   Bool_t          h0_L0ElectronDecision_TIS;
   Bool_t          h0_L0ElectronDecision_TOS;
   Bool_t          h0_L0PhotonDecision_Dec;
   Bool_t          h0_L0PhotonDecision_TIS;
   Bool_t          h0_L0PhotonDecision_TOS;
   Bool_t          h0_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          h0_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          h0_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          h0_Hlt1DiMuonLowMassDecision_Dec;
   Bool_t          h0_Hlt1DiMuonLowMassDecision_TIS;
   Bool_t          h0_Hlt1DiMuonLowMassDecision_TOS;
   Bool_t          h0_Hlt1SingleMuonNoIPDecision_Dec;
   Bool_t          h0_Hlt1SingleMuonNoIPDecision_TIS;
   Bool_t          h0_Hlt1SingleMuonNoIPDecision_TOS;
   Bool_t          h0_Hlt1SingleMuonHighPTDecision_Dec;
   Bool_t          h0_Hlt1SingleMuonHighPTDecision_TIS;
   Bool_t          h0_Hlt1SingleMuonHighPTDecision_TOS;
   Bool_t          h0_Hlt1TrackAllL0Decision_Dec;
   Bool_t          h0_Hlt1TrackAllL0Decision_TIS;
   Bool_t          h0_Hlt1TrackAllL0Decision_TOS;
   Bool_t          h0_Hlt1TrackMuonDecision_Dec;
   Bool_t          h0_Hlt1TrackMuonDecision_TIS;
   Bool_t          h0_Hlt1TrackMuonDecision_TOS;
   Bool_t          h0_Hlt1TrackPhotonDecision_Dec;
   Bool_t          h0_Hlt1TrackPhotonDecision_TIS;
   Bool_t          h0_Hlt1TrackPhotonDecision_TOS;
   Bool_t          h0_Hlt1L0AnyDecision_Dec;
   Bool_t          h0_Hlt1L0AnyDecision_TIS;
   Bool_t          h0_Hlt1L0AnyDecision_TOS;
   Bool_t          h0_Hlt1GlobalDecision_Dec;
   Bool_t          h0_Hlt1GlobalDecision_TIS;
   Bool_t          h0_Hlt1GlobalDecision_TOS;
   Bool_t          h0_Hlt2SingleMuonDecision_Dec;
   Bool_t          h0_Hlt2SingleMuonDecision_TIS;
   Bool_t          h0_Hlt2SingleMuonDecision_TOS;
   Bool_t          h0_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          h0_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          h0_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          h0_Hlt2CharmSemilepD2HMuMuDecision_Dec;
   Bool_t          h0_Hlt2CharmSemilepD2HMuMuDecision_TIS;
   Bool_t          h0_Hlt2CharmSemilepD2HMuMuDecision_TOS;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;
   Bool_t          h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;
   Bool_t          h0_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;
   Bool_t          h0_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;
   Bool_t          h0_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;
   Bool_t          h0_Hlt2CharmSemilepD02KKMuMuDecision_Dec;
   Bool_t          h0_Hlt2CharmSemilepD02KKMuMuDecision_TIS;
   Bool_t          h0_Hlt2CharmSemilepD02KKMuMuDecision_TOS;
   Bool_t          h0_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;
   Bool_t          h0_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;
   Bool_t          h0_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_K3piDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHHH_K3piDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_K3piDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_4piDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHHH_4piDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_4piDecision_TOS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;
   Bool_t          h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;
   Bool_t          h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;
   Double_t        h1_MINIP;
   Double_t        h1_MINIPCHI2;
   Double_t        h1_MINIPNEXTBEST;
   Double_t        h1_MINIPCHI2NEXTBEST;
   Double_t        h1_OWNPV_X;
   Double_t        h1_OWNPV_Y;
   Double_t        h1_OWNPV_Z;
   Double_t        h1_OWNPV_XERR;
   Double_t        h1_OWNPV_YERR;
   Double_t        h1_OWNPV_ZERR;
   Double_t        h1_OWNPV_CHI2;
   Int_t           h1_OWNPV_NDOF;
   Float_t         h1_OWNPV_COV_[3][3];
   Double_t        h1_IP_OWNPV;
   Double_t        h1_IPCHI2_OWNPV;
   Double_t        h1_TOPPV_X;
   Double_t        h1_TOPPV_Y;
   Double_t        h1_TOPPV_Z;
   Double_t        h1_TOPPV_XERR;
   Double_t        h1_TOPPV_YERR;
   Double_t        h1_TOPPV_ZERR;
   Double_t        h1_TOPPV_CHI2;
   Int_t           h1_TOPPV_NDOF;
   Float_t         h1_TOPPV_COV_[3][3];
   Double_t        h1_IP_TOPPV;
   Double_t        h1_IPCHI2_TOPPV;
   Double_t        h1_ORIVX_X;
   Double_t        h1_ORIVX_Y;
   Double_t        h1_ORIVX_Z;
   Double_t        h1_ORIVX_XERR;
   Double_t        h1_ORIVX_YERR;
   Double_t        h1_ORIVX_ZERR;
   Double_t        h1_ORIVX_CHI2;
   Int_t           h1_ORIVX_NDOF;
   Float_t         h1_ORIVX_COV_[3][3];
   Double_t        h1_IP_ORIVX;
   Double_t        h1_IPCHI2_ORIVX;
   Double_t        h1_P;
   Double_t        h1_PT;
   Double_t        h1_PE;
   Double_t        h1_PX;
   Double_t        h1_PY;
   Double_t        h1_PZ;
   Double_t        h1_M;
   Double_t        h1_L0Calo_HCAL_realET;
   Double_t        h1_L0Calo_HCAL_xProjection;
   Double_t        h1_L0Calo_HCAL_yProjection;
   Int_t           h1_L0Calo_HCAL_region;
   Double_t        h1_L0Calo_HCAL_TriggerET;
   Double_t        h1_L0Calo_HCAL_TriggerHCALET;
   Double_t        h1_L0Calo_HCAL_xTrigger;
   Double_t        h1_L0Calo_HCAL_yTrigger;
   Int_t           h1_ID;
   Double_t        h1_CombDLLMu;
   Double_t        h1_ProbNNmu;
   Double_t        h1_ProbNNghost;
   Double_t        h1_InMuonAcc;
   Double_t        h1_MuonDist2;
   Int_t           h1_regionInM2;
   Bool_t          h1_hasMuon;
   Bool_t          h1_isMuon;
   Bool_t          h1_isMuonLoose;
   Int_t           h1_NShared;
   Double_t        h1_MuonLLmu;
   Double_t        h1_MuonLLbg;
   Int_t           h1_isMuonFromProto;
   Double_t        h1_PIDe;
   Double_t        h1_PIDmu;
   Double_t        h1_PIDK;
   Double_t        h1_PIDp;
   Double_t        h1_ProbNNe;
   Double_t        h1_ProbNNk;
   Double_t        h1_ProbNNp;
   Double_t        h1_ProbNNpi;
   Bool_t          h1_hasRich;
   Bool_t          h1_hasCalo;
   Bool_t          h1_UsedRichAerogel;
   Bool_t          h1_UsedRich1Gas;
   Bool_t          h1_UsedRich2Gas;
   Bool_t          h1_RichAboveElThres;
   Bool_t          h1_RichAboveMuThres;
   Bool_t          h1_RichAbovePiThres;
   Bool_t          h1_RichAboveKaThres;
   Bool_t          h1_RichAbovePrThres;
   Double_t        h1_RichDLLe;
   Double_t        h1_RichDLLmu;
   Double_t        h1_RichDLLpi;
   Double_t        h1_RichDLLk;
   Double_t        h1_RichDLLp;
   Double_t        h1_RichDLLbt;
   Bool_t          h1_InAccMuon;
   Double_t        h1_MuonMuLL;
   Double_t        h1_MuonBkgLL;
   Int_t           h1_MuonNShared;
   Bool_t          h1_InAccEcal;
   Double_t        h1_CaloEcalE;
   Double_t        h1_EcalPIDe;
   Double_t        h1_EcalPIDmu;
   Bool_t          h1_InAccHcal;
   Double_t        h1_CaloHcalE;
   Double_t        h1_HcalPIDe;
   Double_t        h1_HcalPIDmu;
   Bool_t          h1_InAccPrs;
   Double_t        h1_PrsPIDe;
   Double_t        h1_CaloPrsE;
   Bool_t          h1_InAccSpd;
   Double_t        h1_CaloSpdE;
   Bool_t          h1_InAccBrem;
   Double_t        h1_BremPIDe;
   Double_t        h1_VeloCharge;
   Double_t        h1_RICHDLLe;
   Double_t        h1_RICHDLLmu;
   Double_t        h1_RICHDLLpi;
   Double_t        h1_RICHDLLK;
   Double_t        h1_RICHDLLp;
   Double_t        h1_RICHDLLbt;
   Int_t           h1_RICHBestID;
   Int_t           h1_RICHThreshold;
   Int_t           h1_RICHThresholdEl;
   Int_t           h1_RICHThresholdMu;
   Int_t           h1_RICHThresholdPi;
   Int_t           h1_RICHThresholdKa;
   Int_t           h1_RICHThresholdPr;
   Int_t           h1_RICHAerogelUsed;
   Int_t           h1_RICH1GasUsed;
   Int_t           h1_RICH2GasUsed;
   Double_t        h1_TRACK_Eta;
   Double_t        h1_TRACK_Phi;
   Double_t        h1_Aerogel_X;
   Double_t        h1_Aerogel_Y;
   Double_t        h1_Aerogel_Z;
   Double_t        h1_Aerogel_Rho;
   Double_t        h1_Aerogel_Phi;
   Double_t        h1_Rich1Gas_X;
   Double_t        h1_Rich1Gas_Y;
   Double_t        h1_Rich1Gas_Z;
   Double_t        h1_Rich1Gas_Rho;
   Double_t        h1_Rich1Gas_Phi;
   Double_t        h1_Rich2Gas_X;
   Double_t        h1_Rich2Gas_Y;
   Double_t        h1_Rich2Gas_Z;
   Double_t        h1_Rich2Gas_Rho;
   Double_t        h1_Rich2Gas_Phi;
   Bool_t          h1_L0Global_Dec;
   Bool_t          h1_L0Global_TIS;
   Bool_t          h1_L0Global_TOS;
   Bool_t          h1_Hlt1Global_Dec;
   Bool_t          h1_Hlt1Global_TIS;
   Bool_t          h1_Hlt1Global_TOS;
   Bool_t          h1_Hlt1Phys_Dec;
   Bool_t          h1_Hlt1Phys_TIS;
   Bool_t          h1_Hlt1Phys_TOS;
   Bool_t          h1_Hlt2Global_Dec;
   Bool_t          h1_Hlt2Global_TIS;
   Bool_t          h1_Hlt2Global_TOS;
   Bool_t          h1_Hlt2Phys_Dec;
   Bool_t          h1_Hlt2Phys_TIS;
   Bool_t          h1_Hlt2Phys_TOS;
   Int_t           h1_TRACK_Type;
   Int_t           h1_TRACK_Key;
   Double_t        h1_TRACK_CHI2;
   Int_t           h1_TRACK_NDOF;
   Double_t        h1_TRACK_CHI2NDOF;
   Double_t        h1_TRACK_PCHI2;
   Double_t        h1_TRACK_VeloCHI2NDOF;
   Double_t        h1_TRACK_TCHI2NDOF;
   Double_t        h1_VELO_UTID;
   Double_t        h1_TRACK_FirstMeasurementX;
   Double_t        h1_TRACK_FirstMeasurementY;
   Double_t        h1_TRACK_FirstMeasurementZ;
   Double_t        h1_TRACK_MatchCHI2;
   Double_t        h1_TRACK_GhostProb;
   Double_t        h1_TRACK_CloneDist;
   Double_t        h1_TRACK_Likelihood;
   Bool_t          h1_L0HadronDecision_Dec;
   Bool_t          h1_L0HadronDecision_TIS;
   Bool_t          h1_L0HadronDecision_TOS;
   Bool_t          h1_L0MuonDecision_Dec;
   Bool_t          h1_L0MuonDecision_TIS;
   Bool_t          h1_L0MuonDecision_TOS;
   Bool_t          h1_L0DiMuonDecision_Dec;
   Bool_t          h1_L0DiMuonDecision_TIS;
   Bool_t          h1_L0DiMuonDecision_TOS;
   Bool_t          h1_L0ElectronDecision_Dec;
   Bool_t          h1_L0ElectronDecision_TIS;
   Bool_t          h1_L0ElectronDecision_TOS;
   Bool_t          h1_L0PhotonDecision_Dec;
   Bool_t          h1_L0PhotonDecision_TIS;
   Bool_t          h1_L0PhotonDecision_TOS;
   Bool_t          h1_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          h1_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          h1_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          h1_Hlt1DiMuonLowMassDecision_Dec;
   Bool_t          h1_Hlt1DiMuonLowMassDecision_TIS;
   Bool_t          h1_Hlt1DiMuonLowMassDecision_TOS;
   Bool_t          h1_Hlt1SingleMuonNoIPDecision_Dec;
   Bool_t          h1_Hlt1SingleMuonNoIPDecision_TIS;
   Bool_t          h1_Hlt1SingleMuonNoIPDecision_TOS;
   Bool_t          h1_Hlt1SingleMuonHighPTDecision_Dec;
   Bool_t          h1_Hlt1SingleMuonHighPTDecision_TIS;
   Bool_t          h1_Hlt1SingleMuonHighPTDecision_TOS;
   Bool_t          h1_Hlt1TrackAllL0Decision_Dec;
   Bool_t          h1_Hlt1TrackAllL0Decision_TIS;
   Bool_t          h1_Hlt1TrackAllL0Decision_TOS;
   Bool_t          h1_Hlt1TrackMuonDecision_Dec;
   Bool_t          h1_Hlt1TrackMuonDecision_TIS;
   Bool_t          h1_Hlt1TrackMuonDecision_TOS;
   Bool_t          h1_Hlt1TrackPhotonDecision_Dec;
   Bool_t          h1_Hlt1TrackPhotonDecision_TIS;
   Bool_t          h1_Hlt1TrackPhotonDecision_TOS;
   Bool_t          h1_Hlt1L0AnyDecision_Dec;
   Bool_t          h1_Hlt1L0AnyDecision_TIS;
   Bool_t          h1_Hlt1L0AnyDecision_TOS;
   Bool_t          h1_Hlt1GlobalDecision_Dec;
   Bool_t          h1_Hlt1GlobalDecision_TIS;
   Bool_t          h1_Hlt1GlobalDecision_TOS;
   Bool_t          h1_Hlt2SingleMuonDecision_Dec;
   Bool_t          h1_Hlt2SingleMuonDecision_TIS;
   Bool_t          h1_Hlt2SingleMuonDecision_TOS;
   Bool_t          h1_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          h1_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          h1_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          h1_Hlt2CharmSemilepD2HMuMuDecision_Dec;
   Bool_t          h1_Hlt2CharmSemilepD2HMuMuDecision_TIS;
   Bool_t          h1_Hlt2CharmSemilepD2HMuMuDecision_TOS;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;
   Bool_t          h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;
   Bool_t          h1_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;
   Bool_t          h1_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;
   Bool_t          h1_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;
   Bool_t          h1_Hlt2CharmSemilepD02KKMuMuDecision_Dec;
   Bool_t          h1_Hlt2CharmSemilepD02KKMuMuDecision_TIS;
   Bool_t          h1_Hlt2CharmSemilepD02KKMuMuDecision_TOS;
   Bool_t          h1_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;
   Bool_t          h1_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;
   Bool_t          h1_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_K3piDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHHH_K3piDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_K3piDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_4piDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHHH_4piDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_4piDecision_TOS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;
   Bool_t          h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;
   Bool_t          h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;
   Double_t        mu0_MINIP;
   Double_t        mu0_MINIPCHI2;
   Double_t        mu0_MINIPNEXTBEST;
   Double_t        mu0_MINIPCHI2NEXTBEST;
   Double_t        mu0_OWNPV_X;
   Double_t        mu0_OWNPV_Y;
   Double_t        mu0_OWNPV_Z;
   Double_t        mu0_OWNPV_XERR;
   Double_t        mu0_OWNPV_YERR;
   Double_t        mu0_OWNPV_ZERR;
   Double_t        mu0_OWNPV_CHI2;
   Int_t           mu0_OWNPV_NDOF;
   Float_t         mu0_OWNPV_COV_[3][3];
   Double_t        mu0_IP_OWNPV;
   Double_t        mu0_IPCHI2_OWNPV;
   Double_t        mu0_TOPPV_X;
   Double_t        mu0_TOPPV_Y;
   Double_t        mu0_TOPPV_Z;
   Double_t        mu0_TOPPV_XERR;
   Double_t        mu0_TOPPV_YERR;
   Double_t        mu0_TOPPV_ZERR;
   Double_t        mu0_TOPPV_CHI2;
   Int_t           mu0_TOPPV_NDOF;
   Float_t         mu0_TOPPV_COV_[3][3];
   Double_t        mu0_IP_TOPPV;
   Double_t        mu0_IPCHI2_TOPPV;
   Double_t        mu0_ORIVX_X;
   Double_t        mu0_ORIVX_Y;
   Double_t        mu0_ORIVX_Z;
   Double_t        mu0_ORIVX_XERR;
   Double_t        mu0_ORIVX_YERR;
   Double_t        mu0_ORIVX_ZERR;
   Double_t        mu0_ORIVX_CHI2;
   Int_t           mu0_ORIVX_NDOF;
   Float_t         mu0_ORIVX_COV_[3][3];
   Double_t        mu0_IP_ORIVX;
   Double_t        mu0_IPCHI2_ORIVX;
   Double_t        mu0_P;
   Double_t        mu0_PT;
   Double_t        mu0_PE;
   Double_t        mu0_PX;
   Double_t        mu0_PY;
   Double_t        mu0_PZ;
   Double_t        mu0_M;
   Double_t        mu0_L0Calo_HCAL_realET;
   Double_t        mu0_L0Calo_HCAL_xProjection;
   Double_t        mu0_L0Calo_HCAL_yProjection;
   Int_t           mu0_L0Calo_HCAL_region;
   Double_t        mu0_L0Calo_HCAL_TriggerET;
   Double_t        mu0_L0Calo_HCAL_TriggerHCALET;
   Double_t        mu0_L0Calo_HCAL_xTrigger;
   Double_t        mu0_L0Calo_HCAL_yTrigger;
   Int_t           mu0_ID;
   Double_t        mu0_CombDLLMu;
   Double_t        mu0_ProbNNmu;
   Double_t        mu0_ProbNNghost;
   Double_t        mu0_InMuonAcc;
   Double_t        mu0_MuonDist2;
   Int_t           mu0_regionInM2;
   Bool_t          mu0_hasMuon;
   Bool_t          mu0_isMuon;
   Bool_t          mu0_isMuonLoose;
   Int_t           mu0_NShared;
   Double_t        mu0_MuonLLmu;
   Double_t        mu0_MuonLLbg;
   Int_t           mu0_isMuonFromProto;
   Double_t        mu0_PIDe;
   Double_t        mu0_PIDmu;
   Double_t        mu0_PIDK;
   Double_t        mu0_PIDp;
   Double_t        mu0_ProbNNe;
   Double_t        mu0_ProbNNk;
   Double_t        mu0_ProbNNp;
   Double_t        mu0_ProbNNpi;
   Bool_t          mu0_hasRich;
   Bool_t          mu0_hasCalo;
   Bool_t          mu0_UsedRichAerogel;
   Bool_t          mu0_UsedRich1Gas;
   Bool_t          mu0_UsedRich2Gas;
   Bool_t          mu0_RichAboveElThres;
   Bool_t          mu0_RichAboveMuThres;
   Bool_t          mu0_RichAbovePiThres;
   Bool_t          mu0_RichAboveKaThres;
   Bool_t          mu0_RichAbovePrThres;
   Double_t        mu0_RichDLLe;
   Double_t        mu0_RichDLLmu;
   Double_t        mu0_RichDLLpi;
   Double_t        mu0_RichDLLk;
   Double_t        mu0_RichDLLp;
   Double_t        mu0_RichDLLbt;
   Bool_t          mu0_InAccMuon;
   Double_t        mu0_MuonMuLL;
   Double_t        mu0_MuonBkgLL;
   Int_t           mu0_MuonNShared;
   Bool_t          mu0_InAccEcal;
   Double_t        mu0_CaloEcalE;
   Double_t        mu0_EcalPIDe;
   Double_t        mu0_EcalPIDmu;
   Bool_t          mu0_InAccHcal;
   Double_t        mu0_CaloHcalE;
   Double_t        mu0_HcalPIDe;
   Double_t        mu0_HcalPIDmu;
   Bool_t          mu0_InAccPrs;
   Double_t        mu0_PrsPIDe;
   Double_t        mu0_CaloPrsE;
   Bool_t          mu0_InAccSpd;
   Double_t        mu0_CaloSpdE;
   Bool_t          mu0_InAccBrem;
   Double_t        mu0_BremPIDe;
   Double_t        mu0_VeloCharge;
   Double_t        mu0_RICHDLLe;
   Double_t        mu0_RICHDLLmu;
   Double_t        mu0_RICHDLLpi;
   Double_t        mu0_RICHDLLK;
   Double_t        mu0_RICHDLLp;
   Double_t        mu0_RICHDLLbt;
   Int_t           mu0_RICHBestID;
   Int_t           mu0_RICHThreshold;
   Int_t           mu0_RICHThresholdEl;
   Int_t           mu0_RICHThresholdMu;
   Int_t           mu0_RICHThresholdPi;
   Int_t           mu0_RICHThresholdKa;
   Int_t           mu0_RICHThresholdPr;
   Int_t           mu0_RICHAerogelUsed;
   Int_t           mu0_RICH1GasUsed;
   Int_t           mu0_RICH2GasUsed;
   Double_t        mu0_TRACK_Eta;
   Double_t        mu0_TRACK_Phi;
   Double_t        mu0_Aerogel_X;
   Double_t        mu0_Aerogel_Y;
   Double_t        mu0_Aerogel_Z;
   Double_t        mu0_Aerogel_Rho;
   Double_t        mu0_Aerogel_Phi;
   Double_t        mu0_Rich1Gas_X;
   Double_t        mu0_Rich1Gas_Y;
   Double_t        mu0_Rich1Gas_Z;
   Double_t        mu0_Rich1Gas_Rho;
   Double_t        mu0_Rich1Gas_Phi;
   Double_t        mu0_Rich2Gas_X;
   Double_t        mu0_Rich2Gas_Y;
   Double_t        mu0_Rich2Gas_Z;
   Double_t        mu0_Rich2Gas_Rho;
   Double_t        mu0_Rich2Gas_Phi;
   Bool_t          mu0_L0Global_Dec;
   Bool_t          mu0_L0Global_TIS;
   Bool_t          mu0_L0Global_TOS;
   Bool_t          mu0_Hlt1Global_Dec;
   Bool_t          mu0_Hlt1Global_TIS;
   Bool_t          mu0_Hlt1Global_TOS;
   Bool_t          mu0_Hlt1Phys_Dec;
   Bool_t          mu0_Hlt1Phys_TIS;
   Bool_t          mu0_Hlt1Phys_TOS;
   Bool_t          mu0_Hlt2Global_Dec;
   Bool_t          mu0_Hlt2Global_TIS;
   Bool_t          mu0_Hlt2Global_TOS;
   Bool_t          mu0_Hlt2Phys_Dec;
   Bool_t          mu0_Hlt2Phys_TIS;
   Bool_t          mu0_Hlt2Phys_TOS;
   Int_t           mu0_TRACK_Type;
   Int_t           mu0_TRACK_Key;
   Double_t        mu0_TRACK_CHI2;
   Int_t           mu0_TRACK_NDOF;
   Double_t        mu0_TRACK_CHI2NDOF;
   Double_t        mu0_TRACK_PCHI2;
   Double_t        mu0_TRACK_VeloCHI2NDOF;
   Double_t        mu0_TRACK_TCHI2NDOF;
   Double_t        mu0_VELO_UTID;
   Double_t        mu0_TRACK_FirstMeasurementX;
   Double_t        mu0_TRACK_FirstMeasurementY;
   Double_t        mu0_TRACK_FirstMeasurementZ;
   Double_t        mu0_TRACK_MatchCHI2;
   Double_t        mu0_TRACK_GhostProb;
   Double_t        mu0_TRACK_CloneDist;
   Double_t        mu0_TRACK_Likelihood;
   Bool_t          mu0_L0HadronDecision_Dec;
   Bool_t          mu0_L0HadronDecision_TIS;
   Bool_t          mu0_L0HadronDecision_TOS;
   Bool_t          mu0_L0MuonDecision_Dec;
   Bool_t          mu0_L0MuonDecision_TIS;
   Bool_t          mu0_L0MuonDecision_TOS;
   Bool_t          mu0_L0DiMuonDecision_Dec;
   Bool_t          mu0_L0DiMuonDecision_TIS;
   Bool_t          mu0_L0DiMuonDecision_TOS;
   Bool_t          mu0_L0ElectronDecision_Dec;
   Bool_t          mu0_L0ElectronDecision_TIS;
   Bool_t          mu0_L0ElectronDecision_TOS;
   Bool_t          mu0_L0PhotonDecision_Dec;
   Bool_t          mu0_L0PhotonDecision_TIS;
   Bool_t          mu0_L0PhotonDecision_TOS;
   Bool_t          mu0_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          mu0_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          mu0_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          mu0_Hlt1DiMuonLowMassDecision_Dec;
   Bool_t          mu0_Hlt1DiMuonLowMassDecision_TIS;
   Bool_t          mu0_Hlt1DiMuonLowMassDecision_TOS;
   Bool_t          mu0_Hlt1SingleMuonNoIPDecision_Dec;
   Bool_t          mu0_Hlt1SingleMuonNoIPDecision_TIS;
   Bool_t          mu0_Hlt1SingleMuonNoIPDecision_TOS;
   Bool_t          mu0_Hlt1SingleMuonHighPTDecision_Dec;
   Bool_t          mu0_Hlt1SingleMuonHighPTDecision_TIS;
   Bool_t          mu0_Hlt1SingleMuonHighPTDecision_TOS;
   Bool_t          mu0_Hlt1TrackAllL0Decision_Dec;
   Bool_t          mu0_Hlt1TrackAllL0Decision_TIS;
   Bool_t          mu0_Hlt1TrackAllL0Decision_TOS;
   Bool_t          mu0_Hlt1TrackMuonDecision_Dec;
   Bool_t          mu0_Hlt1TrackMuonDecision_TIS;
   Bool_t          mu0_Hlt1TrackMuonDecision_TOS;
   Bool_t          mu0_Hlt1TrackPhotonDecision_Dec;
   Bool_t          mu0_Hlt1TrackPhotonDecision_TIS;
   Bool_t          mu0_Hlt1TrackPhotonDecision_TOS;
   Bool_t          mu0_Hlt1L0AnyDecision_Dec;
   Bool_t          mu0_Hlt1L0AnyDecision_TIS;
   Bool_t          mu0_Hlt1L0AnyDecision_TOS;
   Bool_t          mu0_Hlt1GlobalDecision_Dec;
   Bool_t          mu0_Hlt1GlobalDecision_TIS;
   Bool_t          mu0_Hlt1GlobalDecision_TOS;
   Bool_t          mu0_Hlt2SingleMuonDecision_Dec;
   Bool_t          mu0_Hlt2SingleMuonDecision_TIS;
   Bool_t          mu0_Hlt2SingleMuonDecision_TOS;
   Bool_t          mu0_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          mu0_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          mu0_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          mu0_Hlt2CharmSemilepD2HMuMuDecision_Dec;
   Bool_t          mu0_Hlt2CharmSemilepD2HMuMuDecision_TIS;
   Bool_t          mu0_Hlt2CharmSemilepD2HMuMuDecision_TOS;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;
   Bool_t          mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;
   Bool_t          mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;
   Bool_t          mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;
   Bool_t          mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;
   Bool_t          mu0_Hlt2CharmSemilepD02KKMuMuDecision_Dec;
   Bool_t          mu0_Hlt2CharmSemilepD02KKMuMuDecision_TIS;
   Bool_t          mu0_Hlt2CharmSemilepD02KKMuMuDecision_TOS;
   Bool_t          mu0_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;
   Bool_t          mu0_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;
   Bool_t          mu0_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_K3piDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_K3piDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_K3piDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_4piDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_4piDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_4piDecision_TOS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;
   Bool_t          mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;
   Double_t        mu1_MINIP;
   Double_t        mu1_MINIPCHI2;
   Double_t        mu1_MINIPNEXTBEST;
   Double_t        mu1_MINIPCHI2NEXTBEST;
   Double_t        mu1_OWNPV_X;
   Double_t        mu1_OWNPV_Y;
   Double_t        mu1_OWNPV_Z;
   Double_t        mu1_OWNPV_XERR;
   Double_t        mu1_OWNPV_YERR;
   Double_t        mu1_OWNPV_ZERR;
   Double_t        mu1_OWNPV_CHI2;
   Int_t           mu1_OWNPV_NDOF;
   Float_t         mu1_OWNPV_COV_[3][3];
   Double_t        mu1_IP_OWNPV;
   Double_t        mu1_IPCHI2_OWNPV;
   Double_t        mu1_TOPPV_X;
   Double_t        mu1_TOPPV_Y;
   Double_t        mu1_TOPPV_Z;
   Double_t        mu1_TOPPV_XERR;
   Double_t        mu1_TOPPV_YERR;
   Double_t        mu1_TOPPV_ZERR;
   Double_t        mu1_TOPPV_CHI2;
   Int_t           mu1_TOPPV_NDOF;
   Float_t         mu1_TOPPV_COV_[3][3];
   Double_t        mu1_IP_TOPPV;
   Double_t        mu1_IPCHI2_TOPPV;
   Double_t        mu1_ORIVX_X;
   Double_t        mu1_ORIVX_Y;
   Double_t        mu1_ORIVX_Z;
   Double_t        mu1_ORIVX_XERR;
   Double_t        mu1_ORIVX_YERR;
   Double_t        mu1_ORIVX_ZERR;
   Double_t        mu1_ORIVX_CHI2;
   Int_t           mu1_ORIVX_NDOF;
   Float_t         mu1_ORIVX_COV_[3][3];
   Double_t        mu1_IP_ORIVX;
   Double_t        mu1_IPCHI2_ORIVX;
   Double_t        mu1_P;
   Double_t        mu1_PT;
   Double_t        mu1_PE;
   Double_t        mu1_PX;
   Double_t        mu1_PY;
   Double_t        mu1_PZ;
   Double_t        mu1_M;
   Double_t        mu1_L0Calo_HCAL_realET;
   Double_t        mu1_L0Calo_HCAL_xProjection;
   Double_t        mu1_L0Calo_HCAL_yProjection;
   Int_t           mu1_L0Calo_HCAL_region;
   Double_t        mu1_L0Calo_HCAL_TriggerET;
   Double_t        mu1_L0Calo_HCAL_TriggerHCALET;
   Double_t        mu1_L0Calo_HCAL_xTrigger;
   Double_t        mu1_L0Calo_HCAL_yTrigger;
   Int_t           mu1_ID;
   Double_t        mu1_CombDLLMu;
   Double_t        mu1_ProbNNmu;
   Double_t        mu1_ProbNNghost;
   Double_t        mu1_InMuonAcc;
   Double_t        mu1_MuonDist2;
   Int_t           mu1_regionInM2;
   Bool_t          mu1_hasMuon;
   Bool_t          mu1_isMuon;
   Bool_t          mu1_isMuonLoose;
   Int_t           mu1_NShared;
   Double_t        mu1_MuonLLmu;
   Double_t        mu1_MuonLLbg;
   Int_t           mu1_isMuonFromProto;
   Double_t        mu1_PIDe;
   Double_t        mu1_PIDmu;
   Double_t        mu1_PIDK;
   Double_t        mu1_PIDp;
   Double_t        mu1_ProbNNe;
   Double_t        mu1_ProbNNk;
   Double_t        mu1_ProbNNp;
   Double_t        mu1_ProbNNpi;
   Bool_t          mu1_hasRich;
   Bool_t          mu1_hasCalo;
   Bool_t          mu1_UsedRichAerogel;
   Bool_t          mu1_UsedRich1Gas;
   Bool_t          mu1_UsedRich2Gas;
   Bool_t          mu1_RichAboveElThres;
   Bool_t          mu1_RichAboveMuThres;
   Bool_t          mu1_RichAbovePiThres;
   Bool_t          mu1_RichAboveKaThres;
   Bool_t          mu1_RichAbovePrThres;
   Double_t        mu1_RichDLLe;
   Double_t        mu1_RichDLLmu;
   Double_t        mu1_RichDLLpi;
   Double_t        mu1_RichDLLk;
   Double_t        mu1_RichDLLp;
   Double_t        mu1_RichDLLbt;
   Bool_t          mu1_InAccMuon;
   Double_t        mu1_MuonMuLL;
   Double_t        mu1_MuonBkgLL;
   Int_t           mu1_MuonNShared;
   Bool_t          mu1_InAccEcal;
   Double_t        mu1_CaloEcalE;
   Double_t        mu1_EcalPIDe;
   Double_t        mu1_EcalPIDmu;
   Bool_t          mu1_InAccHcal;
   Double_t        mu1_CaloHcalE;
   Double_t        mu1_HcalPIDe;
   Double_t        mu1_HcalPIDmu;
   Bool_t          mu1_InAccPrs;
   Double_t        mu1_PrsPIDe;
   Double_t        mu1_CaloPrsE;
   Bool_t          mu1_InAccSpd;
   Double_t        mu1_CaloSpdE;
   Bool_t          mu1_InAccBrem;
   Double_t        mu1_BremPIDe;
   Double_t        mu1_VeloCharge;
   Double_t        mu1_RICHDLLe;
   Double_t        mu1_RICHDLLmu;
   Double_t        mu1_RICHDLLpi;
   Double_t        mu1_RICHDLLK;
   Double_t        mu1_RICHDLLp;
   Double_t        mu1_RICHDLLbt;
   Int_t           mu1_RICHBestID;
   Int_t           mu1_RICHThreshold;
   Int_t           mu1_RICHThresholdEl;
   Int_t           mu1_RICHThresholdMu;
   Int_t           mu1_RICHThresholdPi;
   Int_t           mu1_RICHThresholdKa;
   Int_t           mu1_RICHThresholdPr;
   Int_t           mu1_RICHAerogelUsed;
   Int_t           mu1_RICH1GasUsed;
   Int_t           mu1_RICH2GasUsed;
   Double_t        mu1_TRACK_Eta;
   Double_t        mu1_TRACK_Phi;
   Double_t        mu1_Aerogel_X;
   Double_t        mu1_Aerogel_Y;
   Double_t        mu1_Aerogel_Z;
   Double_t        mu1_Aerogel_Rho;
   Double_t        mu1_Aerogel_Phi;
   Double_t        mu1_Rich1Gas_X;
   Double_t        mu1_Rich1Gas_Y;
   Double_t        mu1_Rich1Gas_Z;
   Double_t        mu1_Rich1Gas_Rho;
   Double_t        mu1_Rich1Gas_Phi;
   Double_t        mu1_Rich2Gas_X;
   Double_t        mu1_Rich2Gas_Y;
   Double_t        mu1_Rich2Gas_Z;
   Double_t        mu1_Rich2Gas_Rho;
   Double_t        mu1_Rich2Gas_Phi;
   Bool_t          mu1_L0Global_Dec;
   Bool_t          mu1_L0Global_TIS;
   Bool_t          mu1_L0Global_TOS;
   Bool_t          mu1_Hlt1Global_Dec;
   Bool_t          mu1_Hlt1Global_TIS;
   Bool_t          mu1_Hlt1Global_TOS;
   Bool_t          mu1_Hlt1Phys_Dec;
   Bool_t          mu1_Hlt1Phys_TIS;
   Bool_t          mu1_Hlt1Phys_TOS;
   Bool_t          mu1_Hlt2Global_Dec;
   Bool_t          mu1_Hlt2Global_TIS;
   Bool_t          mu1_Hlt2Global_TOS;
   Bool_t          mu1_Hlt2Phys_Dec;
   Bool_t          mu1_Hlt2Phys_TIS;
   Bool_t          mu1_Hlt2Phys_TOS;
   Int_t           mu1_TRACK_Type;
   Int_t           mu1_TRACK_Key;
   Double_t        mu1_TRACK_CHI2;
   Int_t           mu1_TRACK_NDOF;
   Double_t        mu1_TRACK_CHI2NDOF;
   Double_t        mu1_TRACK_PCHI2;
   Double_t        mu1_TRACK_VeloCHI2NDOF;
   Double_t        mu1_TRACK_TCHI2NDOF;
   Double_t        mu1_VELO_UTID;
   Double_t        mu1_TRACK_FirstMeasurementX;
   Double_t        mu1_TRACK_FirstMeasurementY;
   Double_t        mu1_TRACK_FirstMeasurementZ;
   Double_t        mu1_TRACK_MatchCHI2;
   Double_t        mu1_TRACK_GhostProb;
   Double_t        mu1_TRACK_CloneDist;
   Double_t        mu1_TRACK_Likelihood;
   Bool_t          mu1_L0HadronDecision_Dec;
   Bool_t          mu1_L0HadronDecision_TIS;
   Bool_t          mu1_L0HadronDecision_TOS;
   Bool_t          mu1_L0MuonDecision_Dec;
   Bool_t          mu1_L0MuonDecision_TIS;
   Bool_t          mu1_L0MuonDecision_TOS;
   Bool_t          mu1_L0DiMuonDecision_Dec;
   Bool_t          mu1_L0DiMuonDecision_TIS;
   Bool_t          mu1_L0DiMuonDecision_TOS;
   Bool_t          mu1_L0ElectronDecision_Dec;
   Bool_t          mu1_L0ElectronDecision_TIS;
   Bool_t          mu1_L0ElectronDecision_TOS;
   Bool_t          mu1_L0PhotonDecision_Dec;
   Bool_t          mu1_L0PhotonDecision_TIS;
   Bool_t          mu1_L0PhotonDecision_TOS;
   Bool_t          mu1_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          mu1_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          mu1_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          mu1_Hlt1DiMuonLowMassDecision_Dec;
   Bool_t          mu1_Hlt1DiMuonLowMassDecision_TIS;
   Bool_t          mu1_Hlt1DiMuonLowMassDecision_TOS;
   Bool_t          mu1_Hlt1SingleMuonNoIPDecision_Dec;
   Bool_t          mu1_Hlt1SingleMuonNoIPDecision_TIS;
   Bool_t          mu1_Hlt1SingleMuonNoIPDecision_TOS;
   Bool_t          mu1_Hlt1SingleMuonHighPTDecision_Dec;
   Bool_t          mu1_Hlt1SingleMuonHighPTDecision_TIS;
   Bool_t          mu1_Hlt1SingleMuonHighPTDecision_TOS;
   Bool_t          mu1_Hlt1TrackAllL0Decision_Dec;
   Bool_t          mu1_Hlt1TrackAllL0Decision_TIS;
   Bool_t          mu1_Hlt1TrackAllL0Decision_TOS;
   Bool_t          mu1_Hlt1TrackMuonDecision_Dec;
   Bool_t          mu1_Hlt1TrackMuonDecision_TIS;
   Bool_t          mu1_Hlt1TrackMuonDecision_TOS;
   Bool_t          mu1_Hlt1TrackPhotonDecision_Dec;
   Bool_t          mu1_Hlt1TrackPhotonDecision_TIS;
   Bool_t          mu1_Hlt1TrackPhotonDecision_TOS;
   Bool_t          mu1_Hlt1L0AnyDecision_Dec;
   Bool_t          mu1_Hlt1L0AnyDecision_TIS;
   Bool_t          mu1_Hlt1L0AnyDecision_TOS;
   Bool_t          mu1_Hlt1GlobalDecision_Dec;
   Bool_t          mu1_Hlt1GlobalDecision_TIS;
   Bool_t          mu1_Hlt1GlobalDecision_TOS;
   Bool_t          mu1_Hlt2SingleMuonDecision_Dec;
   Bool_t          mu1_Hlt2SingleMuonDecision_TIS;
   Bool_t          mu1_Hlt2SingleMuonDecision_TOS;
   Bool_t          mu1_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          mu1_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          mu1_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          mu1_Hlt2CharmSemilepD2HMuMuDecision_Dec;
   Bool_t          mu1_Hlt2CharmSemilepD2HMuMuDecision_TIS;
   Bool_t          mu1_Hlt2CharmSemilepD2HMuMuDecision_TOS;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;
   Bool_t          mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;
   Bool_t          mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;
   Bool_t          mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;
   Bool_t          mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;
   Bool_t          mu1_Hlt2CharmSemilepD02KKMuMuDecision_Dec;
   Bool_t          mu1_Hlt2CharmSemilepD02KKMuMuDecision_TIS;
   Bool_t          mu1_Hlt2CharmSemilepD02KKMuMuDecision_TOS;
   Bool_t          mu1_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;
   Bool_t          mu1_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;
   Bool_t          mu1_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_K3piDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_K3piDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_K3piDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_4piDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_4piDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_4piDecision_TOS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;
   Bool_t          mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;
   Double_t        Slowpi_MINIP;
   Double_t        Slowpi_MINIPCHI2;
   Double_t        Slowpi_MINIPNEXTBEST;
   Double_t        Slowpi_MINIPCHI2NEXTBEST;
   Double_t        Slowpi_OWNPV_X;
   Double_t        Slowpi_OWNPV_Y;
   Double_t        Slowpi_OWNPV_Z;
   Double_t        Slowpi_OWNPV_XERR;
   Double_t        Slowpi_OWNPV_YERR;
   Double_t        Slowpi_OWNPV_ZERR;
   Double_t        Slowpi_OWNPV_CHI2;
   Int_t           Slowpi_OWNPV_NDOF;
   Float_t         Slowpi_OWNPV_COV_[3][3];
   Double_t        Slowpi_IP_OWNPV;
   Double_t        Slowpi_IPCHI2_OWNPV;
   Double_t        Slowpi_TOPPV_X;
   Double_t        Slowpi_TOPPV_Y;
   Double_t        Slowpi_TOPPV_Z;
   Double_t        Slowpi_TOPPV_XERR;
   Double_t        Slowpi_TOPPV_YERR;
   Double_t        Slowpi_TOPPV_ZERR;
   Double_t        Slowpi_TOPPV_CHI2;
   Int_t           Slowpi_TOPPV_NDOF;
   Float_t         Slowpi_TOPPV_COV_[3][3];
   Double_t        Slowpi_IP_TOPPV;
   Double_t        Slowpi_IPCHI2_TOPPV;
   Double_t        Slowpi_ORIVX_X;
   Double_t        Slowpi_ORIVX_Y;
   Double_t        Slowpi_ORIVX_Z;
   Double_t        Slowpi_ORIVX_XERR;
   Double_t        Slowpi_ORIVX_YERR;
   Double_t        Slowpi_ORIVX_ZERR;
   Double_t        Slowpi_ORIVX_CHI2;
   Int_t           Slowpi_ORIVX_NDOF;
   Float_t         Slowpi_ORIVX_COV_[3][3];
   Double_t        Slowpi_IP_ORIVX;
   Double_t        Slowpi_IPCHI2_ORIVX;
   Double_t        Slowpi_P;
   Double_t        Slowpi_PT;
   Double_t        Slowpi_PE;
   Double_t        Slowpi_PX;
   Double_t        Slowpi_PY;
   Double_t        Slowpi_PZ;
   Double_t        Slowpi_M;
   Double_t        Slowpi_L0Calo_HCAL_realET;
   Double_t        Slowpi_L0Calo_HCAL_xProjection;
   Double_t        Slowpi_L0Calo_HCAL_yProjection;
   Int_t           Slowpi_L0Calo_HCAL_region;
   Double_t        Slowpi_L0Calo_HCAL_TriggerET;
   Double_t        Slowpi_L0Calo_HCAL_TriggerHCALET;
   Double_t        Slowpi_L0Calo_HCAL_xTrigger;
   Double_t        Slowpi_L0Calo_HCAL_yTrigger;
   Int_t           Slowpi_ID;
   Double_t        Slowpi_CombDLLMu;
   Double_t        Slowpi_ProbNNmu;
   Double_t        Slowpi_ProbNNghost;
   Double_t        Slowpi_InMuonAcc;
   Double_t        Slowpi_MuonDist2;
   Int_t           Slowpi_regionInM2;
   Bool_t          Slowpi_hasMuon;
   Bool_t          Slowpi_isMuon;
   Bool_t          Slowpi_isMuonLoose;
   Int_t           Slowpi_NShared;
   Double_t        Slowpi_MuonLLmu;
   Double_t        Slowpi_MuonLLbg;
   Int_t           Slowpi_isMuonFromProto;
   Double_t        Slowpi_PIDe;
   Double_t        Slowpi_PIDmu;
   Double_t        Slowpi_PIDK;
   Double_t        Slowpi_PIDp;
   Double_t        Slowpi_ProbNNe;
   Double_t        Slowpi_ProbNNk;
   Double_t        Slowpi_ProbNNp;
   Double_t        Slowpi_ProbNNpi;
   Bool_t          Slowpi_hasRich;
   Bool_t          Slowpi_hasCalo;
   Bool_t          Slowpi_UsedRichAerogel;
   Bool_t          Slowpi_UsedRich1Gas;
   Bool_t          Slowpi_UsedRich2Gas;
   Bool_t          Slowpi_RichAboveElThres;
   Bool_t          Slowpi_RichAboveMuThres;
   Bool_t          Slowpi_RichAbovePiThres;
   Bool_t          Slowpi_RichAboveKaThres;
   Bool_t          Slowpi_RichAbovePrThres;
   Double_t        Slowpi_RichDLLe;
   Double_t        Slowpi_RichDLLmu;
   Double_t        Slowpi_RichDLLpi;
   Double_t        Slowpi_RichDLLk;
   Double_t        Slowpi_RichDLLp;
   Double_t        Slowpi_RichDLLbt;
   Bool_t          Slowpi_InAccMuon;
   Double_t        Slowpi_MuonMuLL;
   Double_t        Slowpi_MuonBkgLL;
   Int_t           Slowpi_MuonNShared;
   Bool_t          Slowpi_InAccEcal;
   Double_t        Slowpi_CaloEcalE;
   Double_t        Slowpi_EcalPIDe;
   Double_t        Slowpi_EcalPIDmu;
   Bool_t          Slowpi_InAccHcal;
   Double_t        Slowpi_CaloHcalE;
   Double_t        Slowpi_HcalPIDe;
   Double_t        Slowpi_HcalPIDmu;
   Bool_t          Slowpi_InAccPrs;
   Double_t        Slowpi_PrsPIDe;
   Double_t        Slowpi_CaloPrsE;
   Bool_t          Slowpi_InAccSpd;
   Double_t        Slowpi_CaloSpdE;
   Bool_t          Slowpi_InAccBrem;
   Double_t        Slowpi_BremPIDe;
   Double_t        Slowpi_VeloCharge;
   Double_t        Slowpi_RICHDLLe;
   Double_t        Slowpi_RICHDLLmu;
   Double_t        Slowpi_RICHDLLpi;
   Double_t        Slowpi_RICHDLLK;
   Double_t        Slowpi_RICHDLLp;
   Double_t        Slowpi_RICHDLLbt;
   Int_t           Slowpi_RICHBestID;
   Int_t           Slowpi_RICHThreshold;
   Int_t           Slowpi_RICHThresholdEl;
   Int_t           Slowpi_RICHThresholdMu;
   Int_t           Slowpi_RICHThresholdPi;
   Int_t           Slowpi_RICHThresholdKa;
   Int_t           Slowpi_RICHThresholdPr;
   Int_t           Slowpi_RICHAerogelUsed;
   Int_t           Slowpi_RICH1GasUsed;
   Int_t           Slowpi_RICH2GasUsed;
   Double_t        Slowpi_TRACK_Eta;
   Double_t        Slowpi_TRACK_Phi;
   Double_t        Slowpi_Aerogel_X;
   Double_t        Slowpi_Aerogel_Y;
   Double_t        Slowpi_Aerogel_Z;
   Double_t        Slowpi_Aerogel_Rho;
   Double_t        Slowpi_Aerogel_Phi;
   Double_t        Slowpi_Rich1Gas_X;
   Double_t        Slowpi_Rich1Gas_Y;
   Double_t        Slowpi_Rich1Gas_Z;
   Double_t        Slowpi_Rich1Gas_Rho;
   Double_t        Slowpi_Rich1Gas_Phi;
   Double_t        Slowpi_Rich2Gas_X;
   Double_t        Slowpi_Rich2Gas_Y;
   Double_t        Slowpi_Rich2Gas_Z;
   Double_t        Slowpi_Rich2Gas_Rho;
   Double_t        Slowpi_Rich2Gas_Phi;
   Bool_t          Slowpi_L0Global_Dec;
   Bool_t          Slowpi_L0Global_TIS;
   Bool_t          Slowpi_L0Global_TOS;
   Bool_t          Slowpi_Hlt1Global_Dec;
   Bool_t          Slowpi_Hlt1Global_TIS;
   Bool_t          Slowpi_Hlt1Global_TOS;
   Bool_t          Slowpi_Hlt1Phys_Dec;
   Bool_t          Slowpi_Hlt1Phys_TIS;
   Bool_t          Slowpi_Hlt1Phys_TOS;
   Bool_t          Slowpi_Hlt2Global_Dec;
   Bool_t          Slowpi_Hlt2Global_TIS;
   Bool_t          Slowpi_Hlt2Global_TOS;
   Bool_t          Slowpi_Hlt2Phys_Dec;
   Bool_t          Slowpi_Hlt2Phys_TIS;
   Bool_t          Slowpi_Hlt2Phys_TOS;
   Int_t           Slowpi_TRACK_Type;
   Int_t           Slowpi_TRACK_Key;
   Double_t        Slowpi_TRACK_CHI2;
   Int_t           Slowpi_TRACK_NDOF;
   Double_t        Slowpi_TRACK_CHI2NDOF;
   Double_t        Slowpi_TRACK_PCHI2;
   Double_t        Slowpi_TRACK_VeloCHI2NDOF;
   Double_t        Slowpi_TRACK_TCHI2NDOF;
   Double_t        Slowpi_VELO_UTID;
   Double_t        Slowpi_TRACK_FirstMeasurementX;
   Double_t        Slowpi_TRACK_FirstMeasurementY;
   Double_t        Slowpi_TRACK_FirstMeasurementZ;
   Double_t        Slowpi_TRACK_MatchCHI2;
   Double_t        Slowpi_TRACK_GhostProb;
   Double_t        Slowpi_TRACK_CloneDist;
   Double_t        Slowpi_TRACK_Likelihood;
   UInt_t          nCandidate;
   ULong64_t       totCandidates;
   ULong64_t       EventInSequence;
   Double_t        nBackward;
   Double_t        nDownstream;
   Double_t        nITClusters;
   Double_t        nLong;
   Double_t        nMuon;
   Double_t        nPVs;
   Double_t        nSpdDigits;
   Double_t        nTTClusters;
   Double_t        nTracks;
   Double_t        nUpstream;
   Double_t        nVELO;
   UInt_t          runNumber;
   ULong64_t       eventNumber;
   UInt_t          BCID;
   Int_t           BCType;
   UInt_t          OdinTCK;
   UInt_t          L0DUTCK;
   UInt_t          HLT1TCK;
   UInt_t          HLT2TCK;
   ULong64_t       GpsTime;
   Short_t         Polarity;
   Int_t           nPV;
   Float_t         PVX[100];   //[nPV]
   Float_t         PVY[100];   //[nPV]
   Float_t         PVZ[100];   //[nPV]
   Float_t         PVXERR[100];   //[nPV]
   Float_t         PVYERR[100];   //[nPV]
   Float_t         PVZERR[100];   //[nPV]
   Float_t         PVCHI2[100];   //[nPV]
   Float_t         PVNDOF[100];   //[nPV]
   Float_t         PVNTRACKS[100];   //[nPV]

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
   Int_t           D_MC_GD_MOTHER_ID;
   Int_t           D_MC_GD_MOTHER_KEY;
   Int_t           D_MC_GD_GD_MOTHER_ID;
   Int_t           D_MC_GD_GD_MOTHER_KEY;
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
   Int_t           h0_MC_MOTHER_KEY;
   Int_t           h0_MC_GD_MOTHER_ID;
   Int_t           h0_MC_GD_MOTHER_KEY;
   Int_t           h0_MC_GD_GD_MOTHER_ID;
   Int_t           h0_MC_GD_GD_MOTHER_KEY;
   Double_t        h0_TRUEP_E;
   Double_t        h0_TRUEP_X;
   Double_t        h0_TRUEP_Y;
   Double_t        h0_TRUEP_Z;
   Double_t        h0_TRUEPT;
   Double_t        h0_TRUEORIGINVERTEX_X;
   Double_t        h0_TRUEORIGINVERTEX_Y;
   Double_t        h0_TRUEORIGINVERTEX_Z;
   Double_t        h0_TRUEENDVERTEX_X;
   Double_t        h0_TRUEENDVERTEX_Y;
   Double_t        h0_TRUEENDVERTEX_Z;
   Bool_t          h0_TRUEISSTABLE;
   Double_t        h0_TRUETAU;
   Int_t           h1_TRUEID;
   Int_t           h1_MC_MOTHER_ID;
   Int_t           h1_MC_MOTHER_KEY;
   Int_t           h1_MC_GD_MOTHER_ID;
   Int_t           h1_MC_GD_MOTHER_KEY;
   Int_t           h1_MC_GD_GD_MOTHER_ID;
   Int_t           h1_MC_GD_GD_MOTHER_KEY;
   Double_t        h1_TRUEP_E;
   Double_t        h1_TRUEP_X;
   Double_t        h1_TRUEP_Y;
   Double_t        h1_TRUEP_Z;
   Double_t        h1_TRUEPT;
   Double_t        h1_TRUEORIGINVERTEX_X;
   Double_t        h1_TRUEORIGINVERTEX_Y;
   Double_t        h1_TRUEORIGINVERTEX_Z;
   Double_t        h1_TRUEENDVERTEX_X;
   Double_t        h1_TRUEENDVERTEX_Y;
   Double_t        h1_TRUEENDVERTEX_Z;
   Bool_t          h1_TRUEISSTABLE;
   Double_t        h1_TRUETAU;
   Int_t           mu0_TRUEID;
   Int_t           mu0_MC_MOTHER_ID;
   Int_t           mu0_MC_MOTHER_KEY;
   Int_t           mu0_MC_GD_MOTHER_ID;
   Int_t           mu0_MC_GD_MOTHER_KEY;
   Int_t           mu0_MC_GD_GD_MOTHER_ID;
   Int_t           mu0_MC_GD_GD_MOTHER_KEY;
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
   Int_t           mu1_MC_GD_MOTHER_ID;
   Int_t           mu1_MC_GD_MOTHER_KEY;
   Int_t           mu1_MC_GD_GD_MOTHER_ID;
   Int_t           mu1_MC_GD_GD_MOTHER_KEY;
   Double_t        mu1_TRUEP_E;
   Double_t        mu1_TRUEP_X;
   Double_t        mu1_TRUEP_Y;
   Double_t        mu1_TRUEP_Z;
   Double_t        mu1_TRUEPT;
   Double_t        mu1_TRUEORIGINVERTEX_X;
   Double_t        mu1_TRUEORIGINVERTEX_Y;
   Double_t        mu1_TRUEORIGINVERTEX_Z;
   Double_t        mu1_TRUEENDVERTEX_X;
   Double_t        mu1_TRUEENDVERTEX_Y;
   Double_t        mu1_TRUEENDVERTEX_Z;
   Bool_t          mu1_TRUEISSTABLE;
   Double_t        mu1_TRUETAU;
   Int_t           Slowpi_TRUEID;
   Int_t           Slowpi_MC_MOTHER_ID;
   Int_t           Slowpi_MC_MOTHER_KEY;
   Int_t           Slowpi_MC_GD_MOTHER_ID;
   Int_t           Slowpi_MC_GD_MOTHER_KEY;
   Int_t           Slowpi_MC_GD_GD_MOTHER_ID;
   Int_t           Slowpi_MC_GD_GD_MOTHER_KEY;
   Double_t        Slowpi_TRUEP_E;
   Double_t        Slowpi_TRUEP_X;
   Double_t        Slowpi_TRUEP_Y;
   Double_t        Slowpi_TRUEP_Z;
   Double_t        Slowpi_TRUEPT;
   Double_t        Slowpi_TRUEORIGINVERTEX_X;
   Double_t        Slowpi_TRUEORIGINVERTEX_Y;
   Double_t        Slowpi_TRUEORIGINVERTEX_Z;
   Double_t        Slowpi_TRUEENDVERTEX_X;
   Double_t        Slowpi_TRUEENDVERTEX_Y;
   Double_t        Slowpi_TRUEENDVERTEX_Z;
   Bool_t          Slowpi_TRUEISSTABLE;
   Double_t        Slowpi_TRUETAU;
  



   // List of branches
   TBranch        *b_Dst_MINIP;   //!
   TBranch        *b_Dst_MINIPCHI2;   //!
   TBranch        *b_Dst_MINIPNEXTBEST;   //!
   TBranch        *b_Dst_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_Dst_ENDVERTEX_X;   //!
   TBranch        *b_Dst_ENDVERTEX_Y;   //!
   TBranch        *b_Dst_ENDVERTEX_Z;   //!
   TBranch        *b_Dst_ENDVERTEX_XERR;   //!
   TBranch        *b_Dst_ENDVERTEX_YERR;   //!
   TBranch        *b_Dst_ENDVERTEX_ZERR;   //!
   TBranch        *b_Dst_ENDVERTEX_CHI2;   //!
   TBranch        *b_Dst_ENDVERTEX_NDOF;   //!
   TBranch        *b_Dst_ENDVERTEX_COV_;   //!
   TBranch        *b_Dst_OWNPV_X;   //!
   TBranch        *b_Dst_OWNPV_Y;   //!
   TBranch        *b_Dst_OWNPV_Z;   //!
   TBranch        *b_Dst_OWNPV_XERR;   //!
   TBranch        *b_Dst_OWNPV_YERR;   //!
   TBranch        *b_Dst_OWNPV_ZERR;   //!
   TBranch        *b_Dst_OWNPV_CHI2;   //!
   TBranch        *b_Dst_OWNPV_NDOF;   //!
   TBranch        *b_Dst_OWNPV_COV_;   //!
   TBranch        *b_Dst_IP_OWNPV;   //!
   TBranch        *b_Dst_IPCHI2_OWNPV;   //!
   TBranch        *b_Dst_FD_OWNPV;   //!
   TBranch        *b_Dst_FDCHI2_OWNPV;   //!
   TBranch        *b_Dst_DIRA_OWNPV;   //!
   TBranch        *b_Dst_TOPPV_X;   //!
   TBranch        *b_Dst_TOPPV_Y;   //!
   TBranch        *b_Dst_TOPPV_Z;   //!
   TBranch        *b_Dst_TOPPV_XERR;   //!
   TBranch        *b_Dst_TOPPV_YERR;   //!
   TBranch        *b_Dst_TOPPV_ZERR;   //!
   TBranch        *b_Dst_TOPPV_CHI2;   //!
   TBranch        *b_Dst_TOPPV_NDOF;   //!
   TBranch        *b_Dst_TOPPV_COV_;   //!
   TBranch        *b_Dst_IP_TOPPV;   //!
   TBranch        *b_Dst_IPCHI2_TOPPV;   //!
   TBranch        *b_Dst_FD_TOPPV;   //!
   TBranch        *b_Dst_FDCHI2_TOPPV;   //!
   TBranch        *b_Dst_DIRA_TOPPV;   //!
   TBranch        *b_Dst_P;   //!
   TBranch        *b_Dst_PT;   //!
   TBranch        *b_Dst_PE;   //!
   TBranch        *b_Dst_PX;   //!
   TBranch        *b_Dst_PY;   //!
   TBranch        *b_Dst_PZ;   //!
   TBranch        *b_Dst_MM;   //!
   TBranch        *b_Dst_MMERR;   //!
   TBranch        *b_Dst_M;   //!
   TBranch        *b_Dst_ID;   //!
   TBranch        *b_Dst_TAU;   //!
   TBranch        *b_Dst_TAUERR;   //!
   TBranch        *b_Dst_TAUCHI2;   //!
   TBranch        *b_Dst_L0Global_Dec;   //!
   TBranch        *b_Dst_L0Global_TIS;   //!
   TBranch        *b_Dst_L0Global_TOS;   //!
   TBranch        *b_Dst_Hlt1Global_Dec;   //!
   TBranch        *b_Dst_Hlt1Global_TIS;   //!
   TBranch        *b_Dst_Hlt1Global_TOS;   //!
   TBranch        *b_Dst_Hlt1Phys_Dec;   //!
   TBranch        *b_Dst_Hlt1Phys_TIS;   //!
   TBranch        *b_Dst_Hlt1Phys_TOS;   //!
   TBranch        *b_Dst_Hlt2Global_Dec;   //!
   TBranch        *b_Dst_Hlt2Global_TIS;   //!
   TBranch        *b_Dst_Hlt2Global_TOS;   //!
   TBranch        *b_Dst_Hlt2Phys_Dec;   //!
   TBranch        *b_Dst_Hlt2Phys_TIS;   //!
   TBranch        *b_Dst_Hlt2Phys_TOS;   //!
   TBranch        *b_Dst_DTF_CHI2;   //!
   TBranch        *b_Dst_DTF_D0_BPVIPCHI2;   //!
   TBranch        *b_Dst_DTF_D0_E;   //!
   TBranch        *b_Dst_DTF_D0_M;   //!
   TBranch        *b_Dst_DTF_D0_P;   //!
   TBranch        *b_Dst_DTF_D0_PT;   //!
   TBranch        *b_Dst_DTF_D0_PX;   //!
   TBranch        *b_Dst_DTF_D0_PY;   //!
   TBranch        *b_Dst_DTF_D0_PZ;   //!
   TBranch        *b_Dst_DTF_Dstarplus_E;   //!
   TBranch        *b_Dst_DTF_Dstarplus_M;   //!
   TBranch        *b_Dst_DTF_Dstarplus_P;   //!
   TBranch        *b_Dst_DTF_Dstarplus_PT;   //!
   TBranch        *b_Dst_DTF_Dstarplus_PX;   //!
   TBranch        *b_Dst_DTF_Dstarplus_PY;   //!
   TBranch        *b_Dst_DTF_Dstarplus_PZ;   //!
   TBranch        *b_Dst_DTF_NDOF;   //!
   TBranch        *b_Dst_DTF_Pis_BPVIPCHI2;   //!
   TBranch        *b_Dst_DTF_Pis_E;   //!
   TBranch        *b_Dst_DTF_Pis_M;   //!
   TBranch        *b_Dst_DTF_Pis_P;   //!
   TBranch        *b_Dst_DTF_Pis_PT;   //!
   TBranch        *b_Dst_DTF_Pis_PX;   //!
   TBranch        *b_Dst_DTF_Pis_PY;   //!
   TBranch        *b_Dst_DTF_Pis_PZ;   //!
   TBranch        *b_Dst_DTF_h0_E;   //!
   TBranch        *b_Dst_DTF_h0_P;   //!
   TBranch        *b_Dst_DTF_h0_PT;   //!
   TBranch        *b_Dst_DTF_h0_PX;   //!
   TBranch        *b_Dst_DTF_h0_PY;   //!
   TBranch        *b_Dst_DTF_h0_PZ;   //!
   TBranch        *b_Dst_DTF_h1_E;   //!
   TBranch        *b_Dst_DTF_h1_P;   //!
   TBranch        *b_Dst_DTF_h1_PT;   //!
   TBranch        *b_Dst_DTF_h1_PX;   //!
   TBranch        *b_Dst_DTF_h1_PY;   //!
   TBranch        *b_Dst_DTF_h1_PZ;   //!
   TBranch        *b_Dst_DTF_mu0_E;   //!
   TBranch        *b_Dst_DTF_mu0_P;   //!
   TBranch        *b_Dst_DTF_mu0_PT;   //!
   TBranch        *b_Dst_DTF_mu0_PX;   //!
   TBranch        *b_Dst_DTF_mu0_PY;   //!
   TBranch        *b_Dst_DTF_mu0_PZ;   //!
   TBranch        *b_Dst_DTF_mu1_E;   //!
   TBranch        *b_Dst_DTF_mu1_P;   //!
   TBranch        *b_Dst_DTF_mu1_PT;   //!
   TBranch        *b_Dst_DTF_mu1_PX;   //!
   TBranch        *b_Dst_DTF_mu1_PY;   //!
   TBranch        *b_Dst_DTF_mu1_PZ;   //!
   TBranch        *b_Dst_MAXDOCA;   //!
   TBranch        *b_Dst_CONEANGLE_D;   //!
   TBranch        *b_Dst_CONEANGLE_Dstar;   //!
   TBranch        *b_Dst_CONEMULT_D;   //!
   TBranch        *b_Dst_CONEMULT_Dstar;   //!
   TBranch        *b_Dst_CONEPTASYM_D;   //!
   TBranch        *b_Dst_CONEPTASYM_Dstar; 
   TBranch        *b_Dst_L0HadronDecision_Dec;   //!
   TBranch        *b_Dst_L0HadronDecision_TIS;   //!
   TBranch        *b_Dst_L0HadronDecision_TOS;   //!
   TBranch        *b_Dst_L0MuonDecision_Dec;   //!
   TBranch        *b_Dst_L0MuonDecision_TIS;   //!
   TBranch        *b_Dst_L0MuonDecision_TOS;   //!
   TBranch        *b_Dst_L0DiMuonDecision_Dec;   //!
   TBranch        *b_Dst_L0DiMuonDecision_TIS;   //!
   TBranch        *b_Dst_L0DiMuonDecision_TOS;   //!
   TBranch        *b_Dst_L0ElectronDecision_Dec;   //!
   TBranch        *b_Dst_L0ElectronDecision_TIS;   //!
   TBranch        *b_Dst_L0ElectronDecision_TOS;   //!
   TBranch        *b_Dst_L0PhotonDecision_Dec;   //!
   TBranch        *b_Dst_L0PhotonDecision_TIS;   //!
   TBranch        *b_Dst_L0PhotonDecision_TOS;   //!
   TBranch        *b_Dst_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_Dst_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_Dst_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_Dst_Hlt1DiMuonLowMassDecision_Dec;   //!
   TBranch        *b_Dst_Hlt1DiMuonLowMassDecision_TIS;   //!
   TBranch        *b_Dst_Hlt1DiMuonLowMassDecision_TOS;   //!
   TBranch        *b_Dst_Hlt1SingleMuonNoIPDecision_Dec;   //!
   TBranch        *b_Dst_Hlt1SingleMuonNoIPDecision_TIS;   //!
   TBranch        *b_Dst_Hlt1SingleMuonNoIPDecision_TOS;   //!
   TBranch        *b_Dst_Hlt1SingleMuonHighPTDecision_Dec;   //!
   TBranch        *b_Dst_Hlt1SingleMuonHighPTDecision_TIS;   //!
   TBranch        *b_Dst_Hlt1SingleMuonHighPTDecision_TOS;   //!
   TBranch        *b_Dst_Hlt1TrackAllL0Decision_Dec;   //!
   TBranch        *b_Dst_Hlt1TrackAllL0Decision_TIS;   //!
   TBranch        *b_Dst_Hlt1TrackAllL0Decision_TOS;   //!
   TBranch        *b_Dst_Hlt1TrackMuonDecision_Dec;   //!
   TBranch        *b_Dst_Hlt1TrackMuonDecision_TIS;   //!
   TBranch        *b_Dst_Hlt1TrackMuonDecision_TOS;   //!
   TBranch        *b_Dst_Hlt1TrackPhotonDecision_Dec;   //!
   TBranch        *b_Dst_Hlt1TrackPhotonDecision_TIS;   //!
   TBranch        *b_Dst_Hlt1TrackPhotonDecision_TOS;   //!
   TBranch        *b_Dst_Hlt1L0AnyDecision_Dec;   //!
   TBranch        *b_Dst_Hlt1L0AnyDecision_TIS;   //!
   TBranch        *b_Dst_Hlt1L0AnyDecision_TOS;   //!
   TBranch        *b_Dst_Hlt1GlobalDecision_Dec;   //!
   TBranch        *b_Dst_Hlt1GlobalDecision_TIS;   //!
   TBranch        *b_Dst_Hlt1GlobalDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2SingleMuonDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2SingleMuonDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2SingleMuonDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD2HMuMuDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD2HMuMuDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD2HMuMuDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD02KKMuMuDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD02KKMuMuDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD02KKMuMuDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_K3piDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_K3piDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_K3piDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_4piDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_4piDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_4piDecision_TOS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;   //!
   TBranch        *b_Dst_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;   //!
   TBranch        *b_Dst_cpx_1_00;   //!
   TBranch        *b_Dst_cpy_1_00;   //!
   TBranch        *b_Dst_cpz_1_00;   //!
   TBranch        *b_Dst_cpt_1_00;   //!
   TBranch        *b_Dst_cp_1_00;   //!
   TBranch        *b_Dst_cmult_1_00;   //!
   TBranch        *b_Dst_deltaEta_1_00;   //!
   TBranch        *b_Dst_deltaPhi_1_00;   //!
   TBranch        *b_Dst_pxasy_1_00;   //!
   TBranch        *b_Dst_pyasy_1_00;   //!
   TBranch        *b_Dst_pzasy_1_00;   //!
   TBranch        *b_Dst_pasy_1_00;   //!
   TBranch        *b_Dst_ptasy_1_00;   //!
   TBranch        *b_Dst_cpx_1_10;   //!
   TBranch        *b_Dst_cpy_1_10;   //!
   TBranch        *b_Dst_cpz_1_10;   //!
   TBranch        *b_Dst_cpt_1_10;   //!
   TBranch        *b_Dst_cp_1_10;   //!
   TBranch        *b_Dst_cmult_1_10;   //!
   TBranch        *b_Dst_deltaEta_1_10;   //!
   TBranch        *b_Dst_deltaPhi_1_10;   //!
   TBranch        *b_Dst_pxasy_1_10;   //!
   TBranch        *b_Dst_pyasy_1_10;   //!
   TBranch        *b_Dst_pzasy_1_10;   //!
   TBranch        *b_Dst_pasy_1_10;   //!
   TBranch        *b_Dst_ptasy_1_10;   //!
   TBranch        *b_Dst_cpx_1_20;   //!
   TBranch        *b_Dst_cpy_1_20;   //!
   TBranch        *b_Dst_cpz_1_20;   //!
   TBranch        *b_Dst_cpt_1_20;   //!
   TBranch        *b_Dst_cp_1_20;   //!
   TBranch        *b_Dst_cmult_1_20;   //!
   TBranch        *b_Dst_deltaEta_1_20;   //!
   TBranch        *b_Dst_deltaPhi_1_20;   //!
   TBranch        *b_Dst_pxasy_1_20;   //!
   TBranch        *b_Dst_pyasy_1_20;   //!
   TBranch        *b_Dst_pzasy_1_20;   //!
   TBranch        *b_Dst_pasy_1_20;   //!
   TBranch        *b_Dst_ptasy_1_20;   //!
   TBranch        *b_Dst_cpx_1_30;   //!
   TBranch        *b_Dst_cpy_1_30;   //!
   TBranch        *b_Dst_cpz_1_30;   //!
   TBranch        *b_Dst_cpt_1_30;   //!
   TBranch        *b_Dst_cp_1_30;   //!
   TBranch        *b_Dst_cmult_1_30;   //!
   TBranch        *b_Dst_deltaEta_1_30;   //!
   TBranch        *b_Dst_deltaPhi_1_30;   //!
   TBranch        *b_Dst_pxasy_1_30;   //!
   TBranch        *b_Dst_pyasy_1_30;   //!
   TBranch        *b_Dst_pzasy_1_30;   //!
   TBranch        *b_Dst_pasy_1_30;   //!
   TBranch        *b_Dst_ptasy_1_30;   //!
   TBranch        *b_Dst_cpx_1_40;   //!
   TBranch        *b_Dst_cpy_1_40;   //!
   TBranch        *b_Dst_cpz_1_40;   //!
   TBranch        *b_Dst_cpt_1_40;   //!
   TBranch        *b_Dst_cp_1_40;   //!
   TBranch        *b_Dst_cmult_1_40;   //!
   TBranch        *b_Dst_deltaEta_1_40;   //!
   TBranch        *b_Dst_deltaPhi_1_40;   //!
   TBranch        *b_Dst_pxasy_1_40;   //!
   TBranch        *b_Dst_pyasy_1_40;   //!
   TBranch        *b_Dst_pzasy_1_40;   //!
   TBranch        *b_Dst_pasy_1_40;   //!
   TBranch        *b_Dst_ptasy_1_40;   //!
   TBranch        *b_Dst_cpx_1_50;   //!
   TBranch        *b_Dst_cpy_1_50;   //!
   TBranch        *b_Dst_cpz_1_50;   //!
   TBranch        *b_Dst_cpt_1_50;   //!
   TBranch        *b_Dst_cp_1_50;   //!
   TBranch        *b_Dst_cmult_1_50;   //!
   TBranch        *b_Dst_deltaEta_1_50;   //!
   TBranch        *b_Dst_deltaPhi_1_50;   //!
   TBranch        *b_Dst_pxasy_1_50;   //!
   TBranch        *b_Dst_pyasy_1_50;   //!
   TBranch        *b_Dst_pzasy_1_50;   //!
   TBranch        *b_Dst_pasy_1_50;   //!
   TBranch        *b_Dst_ptasy_1_50;   //!
   TBranch        *b_Dst_cpx_1_60;   //!
   TBranch        *b_Dst_cpy_1_60;   //!
   TBranch        *b_Dst_cpz_1_60;   //!
   TBranch        *b_Dst_cpt_1_60;   //!
   TBranch        *b_Dst_cp_1_60;   //!
   TBranch        *b_Dst_cmult_1_60;   //!
   TBranch        *b_Dst_deltaEta_1_60;   //!
   TBranch        *b_Dst_deltaPhi_1_60;   //!
   TBranch        *b_Dst_pxasy_1_60;   //!
   TBranch        *b_Dst_pyasy_1_60;   //!
   TBranch        *b_Dst_pzasy_1_60;   //!
   TBranch        *b_Dst_pasy_1_60;   //!
   TBranch        *b_Dst_ptasy_1_60;   //!
   TBranch        *b_Dst_cpx_1_70;   //!
   TBranch        *b_Dst_cpy_1_70;   //!
   TBranch        *b_Dst_cpz_1_70;   //!
   TBranch        *b_Dst_cpt_1_70;   //!
   TBranch        *b_Dst_cp_1_70;   //!
   TBranch        *b_Dst_cmult_1_70;   //!
   TBranch        *b_Dst_deltaEta_1_70;   //!
   TBranch        *b_Dst_deltaPhi_1_70;   //!
   TBranch        *b_Dst_pxasy_1_70;   //!
   TBranch        *b_Dst_pyasy_1_70;   //!
   TBranch        *b_Dst_pzasy_1_70;   //!
   TBranch        *b_Dst_pasy_1_70;   //!
   TBranch        *b_Dst_ptasy_1_70;   //!
   TBranch        *b_D_MINIP;   //!
   TBranch        *b_D_MINIPCHI2;   //!
   TBranch        *b_D_MINIPNEXTBEST;   //!
   TBranch        *b_D_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_D_ENDVERTEX_X;   //!
   TBranch        *b_D_ENDVERTEX_Y;   //!
   TBranch        *b_D_ENDVERTEX_Z;   //!
   TBranch        *b_D_ENDVERTEX_XERR;   //!
   TBranch        *b_D_ENDVERTEX_YERR;   //!
   TBranch        *b_D_ENDVERTEX_ZERR;   //!
   TBranch        *b_D_ENDVERTEX_CHI2;   //!
   TBranch        *b_D_ENDVERTEX_NDOF;   //!
   TBranch        *b_D_ENDVERTEX_COV_;   //!
   TBranch        *b_D_OWNPV_X;   //!
   TBranch        *b_D_OWNPV_Y;   //!
   TBranch        *b_D_OWNPV_Z;   //!
   TBranch        *b_D_OWNPV_XERR;   //!
   TBranch        *b_D_OWNPV_YERR;   //!
   TBranch        *b_D_OWNPV_ZERR;   //!
   TBranch        *b_D_OWNPV_CHI2;   //!
   TBranch        *b_D_OWNPV_NDOF;   //!
   TBranch        *b_D_OWNPV_COV_;   //!
   TBranch        *b_D_IP_OWNPV;   //!
   TBranch        *b_D_IPCHI2_OWNPV;   //!
   TBranch        *b_D_FD_OWNPV;   //!
   TBranch        *b_D_FDCHI2_OWNPV;   //!
   TBranch        *b_D_DIRA_OWNPV;   //!
   TBranch        *b_D_TOPPV_X;   //!
   TBranch        *b_D_TOPPV_Y;   //!
   TBranch        *b_D_TOPPV_Z;   //!
   TBranch        *b_D_TOPPV_XERR;   //!
   TBranch        *b_D_TOPPV_YERR;   //!
   TBranch        *b_D_TOPPV_ZERR;   //!
   TBranch        *b_D_TOPPV_CHI2;   //!
   TBranch        *b_D_TOPPV_NDOF;   //!
   TBranch        *b_D_TOPPV_COV_;   //!
   TBranch        *b_D_IP_TOPPV;   //!
   TBranch        *b_D_IPCHI2_TOPPV;   //!
   TBranch        *b_D_FD_TOPPV;   //!
   TBranch        *b_D_FDCHI2_TOPPV;   //!
   TBranch        *b_D_DIRA_TOPPV;   //!
   TBranch        *b_D_ORIVX_X;   //!
   TBranch        *b_D_ORIVX_Y;   //!
   TBranch        *b_D_ORIVX_Z;   //!
   TBranch        *b_D_ORIVX_XERR;   //!
   TBranch        *b_D_ORIVX_YERR;   //!
   TBranch        *b_D_ORIVX_ZERR;   //!
   TBranch        *b_D_ORIVX_CHI2;   //!
   TBranch        *b_D_ORIVX_NDOF;   //!
   TBranch        *b_D_ORIVX_COV_;   //!
   TBranch        *b_D_IP_ORIVX;   //!
   TBranch        *b_D_IPCHI2_ORIVX;   //!
   TBranch        *b_D_FD_ORIVX;   //!
   TBranch        *b_D_FDCHI2_ORIVX;   //!
   TBranch        *b_D_DIRA_ORIVX;   //!
   TBranch        *b_D_P;   //!
   TBranch        *b_D_PT;   //!
   TBranch        *b_D_PE;   //!
   TBranch        *b_D_PX;   //!
   TBranch        *b_D_PY;   //!
   TBranch        *b_D_PZ;   //!
   TBranch        *b_D_MM;   //!
   TBranch        *b_D_MMERR;   //!
   TBranch        *b_D_M;   //!
   TBranch        *b_D_ID;   //!
   TBranch        *b_D_TAU;   //!
   TBranch        *b_D_TAUERR;   //!
   TBranch        *b_D_TAUCHI2;   //!
   TBranch        *b_D_L0Global_Dec;   //!
   TBranch        *b_D_L0Global_TIS;   //!
   TBranch        *b_D_L0Global_TOS;   //!
   TBranch        *b_D_Hlt1Global_Dec;   //!
   TBranch        *b_D_Hlt1Global_TIS;   //!
   TBranch        *b_D_Hlt1Global_TOS;   //!
   TBranch        *b_D_Hlt1Phys_Dec;   //!
   TBranch        *b_D_Hlt1Phys_TIS;   //!
   TBranch        *b_D_Hlt1Phys_TOS;   //!
   TBranch        *b_D_Hlt2Global_Dec;   //!
   TBranch        *b_D_Hlt2Global_TIS;   //!
   TBranch        *b_D_Hlt2Global_TOS;   //!
   TBranch        *b_D_Hlt2Phys_Dec;   //!
   TBranch        *b_D_Hlt2Phys_TIS;   //!
   TBranch        *b_D_Hlt2Phys_TOS;   //!
   TBranch        *b_D_DTF_CHI2;   //!
   TBranch        *b_D_DTF_NDOF;   //!
   TBranch        *b_D_D_DTF_M;   //!
   TBranch        *b_D_D_DTF_M_False;   //!
   TBranch        *b_D_DiMuon_Mass;   //!
   TBranch        *b_D_LV01;   //!
   TBranch        *b_D_LV02;   //!
   TBranch        *b_D_LV03;   //!
   TBranch        *b_D_LV04;   //!
   TBranch        *b_D_MAXDOCA;   //!
   TBranch        *b_D_h0_DTF_P;   //!
   TBranch        *b_D_h0_DTF_PX;   //!
   TBranch        *b_D_h0_DTF_PY;   //!
   TBranch        *b_D_h0_DTF_PZ;   //!
   TBranch        *b_D_h1_DTF_P;   //!
   TBranch        *b_D_h1_DTF_PX;   //!
   TBranch        *b_D_h1_DTF_PY;   //!
   TBranch        *b_D_h1_DTF_PZ;   //!
   TBranch        *b_D_mu0_DTF_P;   //!
   TBranch        *b_D_mu0_DTF_PX;   //!
   TBranch        *b_D_mu0_DTF_PY;   //!
   TBranch        *b_D_mu0_DTF_PZ;   //!
   TBranch        *b_D_mu1_DTF_P;   //!
   TBranch        *b_D_mu1_DTF_PX;   //!
   TBranch        *b_D_mu1_DTF_PY;   //!
   TBranch        *b_D_mu1_DTF_PZ;   //!
   TBranch        *b_D_L0HadronDecision_Dec;   //!
   TBranch        *b_D_L0HadronDecision_TIS;   //!
   TBranch        *b_D_L0HadronDecision_TOS;   //!
   TBranch        *b_D_L0MuonDecision_Dec;   //!
   TBranch        *b_D_L0MuonDecision_TIS;   //!
   TBranch        *b_D_L0MuonDecision_TOS;   //!
   TBranch        *b_D_L0DiMuonDecision_Dec;   //!
   TBranch        *b_D_L0DiMuonDecision_TIS;   //!
   TBranch        *b_D_L0DiMuonDecision_TOS;   //!
   TBranch        *b_D_L0ElectronDecision_Dec;   //!
   TBranch        *b_D_L0ElectronDecision_TIS;   //!
   TBranch        *b_D_L0ElectronDecision_TOS;   //!
   TBranch        *b_D_L0PhotonDecision_Dec;   //!
   TBranch        *b_D_L0PhotonDecision_TIS;   //!
   TBranch        *b_D_L0PhotonDecision_TOS;   //!
   TBranch        *b_D_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_D_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_D_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_D_Hlt1DiMuonLowMassDecision_Dec;   //!
   TBranch        *b_D_Hlt1DiMuonLowMassDecision_TIS;   //!
   TBranch        *b_D_Hlt1DiMuonLowMassDecision_TOS;   //!
   TBranch        *b_D_Hlt1SingleMuonNoIPDecision_Dec;   //!
   TBranch        *b_D_Hlt1SingleMuonNoIPDecision_TIS;   //!
   TBranch        *b_D_Hlt1SingleMuonNoIPDecision_TOS;   //!
   TBranch        *b_D_Hlt1SingleMuonHighPTDecision_Dec;   //!
   TBranch        *b_D_Hlt1SingleMuonHighPTDecision_TIS;   //!
   TBranch        *b_D_Hlt1SingleMuonHighPTDecision_TOS;   //!
   TBranch        *b_D_Hlt1TrackAllL0Decision_Dec;   //!
   TBranch        *b_D_Hlt1TrackAllL0Decision_TIS;   //!
   TBranch        *b_D_Hlt1TrackAllL0Decision_TOS;   //!
   TBranch        *b_D_Hlt1TrackMuonDecision_Dec;   //!
   TBranch        *b_D_Hlt1TrackMuonDecision_TIS;   //!
   TBranch        *b_D_Hlt1TrackMuonDecision_TOS;   //!
   TBranch        *b_D_Hlt1TrackPhotonDecision_Dec;   //!
   TBranch        *b_D_Hlt1TrackPhotonDecision_TIS;   //!
   TBranch        *b_D_Hlt1TrackPhotonDecision_TOS;   //!
   TBranch        *b_D_Hlt1L0AnyDecision_Dec;   //!
   TBranch        *b_D_Hlt1L0AnyDecision_TIS;   //!
   TBranch        *b_D_Hlt1L0AnyDecision_TOS;   //!
   TBranch        *b_D_Hlt1GlobalDecision_Dec;   //!
   TBranch        *b_D_Hlt1GlobalDecision_TIS;   //!
   TBranch        *b_D_Hlt1GlobalDecision_TOS;   //!
   TBranch        *b_D_Hlt2SingleMuonDecision_Dec;   //!
   TBranch        *b_D_Hlt2SingleMuonDecision_TIS;   //!
   TBranch        *b_D_Hlt2SingleMuonDecision_TOS;   //!
   TBranch        *b_D_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_D_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_D_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmSemilepD2HMuMuDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmSemilepD2HMuMuDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmSemilepD2HMuMuDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmSemilepD02KKMuMuDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmSemilepD02KKMuMuDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmSemilepD02KKMuMuDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_K3piDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_K3piDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_K3piDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_4piDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_4piDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_4piDecision_TOS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;   //!
   TBranch        *b_D_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;   //!
   TBranch        *b_D_cpx_1_00;   //!
   TBranch        *b_D_cpy_1_00;   //!
   TBranch        *b_D_cpz_1_00;   //!
   TBranch        *b_D_cpt_1_00;   //!
   TBranch        *b_D_cp_1_00;   //!
   TBranch        *b_D_cmult_1_00;   //!
   TBranch        *b_D_deltaEta_1_00;   //!
   TBranch        *b_D_deltaPhi_1_00;   //!
   TBranch        *b_D_pxasy_1_00;   //!
   TBranch        *b_D_pyasy_1_00;   //!
   TBranch        *b_D_pzasy_1_00;   //!
   TBranch        *b_D_pasy_1_00;   //!
   TBranch        *b_D_ptasy_1_00;   //!
   TBranch        *b_D_cpx_1_10;   //!
   TBranch        *b_D_cpy_1_10;   //!
   TBranch        *b_D_cpz_1_10;   //!
   TBranch        *b_D_cpt_1_10;   //!
   TBranch        *b_D_cp_1_10;   //!
   TBranch        *b_D_cmult_1_10;   //!
   TBranch        *b_D_deltaEta_1_10;   //!
   TBranch        *b_D_deltaPhi_1_10;   //!
   TBranch        *b_D_pxasy_1_10;   //!
   TBranch        *b_D_pyasy_1_10;   //!
   TBranch        *b_D_pzasy_1_10;   //!
   TBranch        *b_D_pasy_1_10;   //!
   TBranch        *b_D_ptasy_1_10;   //!
   TBranch        *b_D_cpx_1_20;   //!
   TBranch        *b_D_cpy_1_20;   //!
   TBranch        *b_D_cpz_1_20;   //!
   TBranch        *b_D_cpt_1_20;   //!
   TBranch        *b_D_cp_1_20;   //!
   TBranch        *b_D_cmult_1_20;   //!
   TBranch        *b_D_deltaEta_1_20;   //!
   TBranch        *b_D_deltaPhi_1_20;   //!
   TBranch        *b_D_pxasy_1_20;   //!
   TBranch        *b_D_pyasy_1_20;   //!
   TBranch        *b_D_pzasy_1_20;   //!
   TBranch        *b_D_pasy_1_20;   //!
   TBranch        *b_D_ptasy_1_20;   //!
   TBranch        *b_D_cpx_1_30;   //!
   TBranch        *b_D_cpy_1_30;   //!
   TBranch        *b_D_cpz_1_30;   //!
   TBranch        *b_D_cpt_1_30;   //!
   TBranch        *b_D_cp_1_30;   //!
   TBranch        *b_D_cmult_1_30;   //!
   TBranch        *b_D_deltaEta_1_30;   //!
   TBranch        *b_D_deltaPhi_1_30;   //!
   TBranch        *b_D_pxasy_1_30;   //!
   TBranch        *b_D_pyasy_1_30;   //!
   TBranch        *b_D_pzasy_1_30;   //!
   TBranch        *b_D_pasy_1_30;   //!
   TBranch        *b_D_ptasy_1_30;   //!
   TBranch        *b_D_cpx_1_40;   //!
   TBranch        *b_D_cpy_1_40;   //!
   TBranch        *b_D_cpz_1_40;   //!
   TBranch        *b_D_cpt_1_40;   //!
   TBranch        *b_D_cp_1_40;   //!
   TBranch        *b_D_cmult_1_40;   //!
   TBranch        *b_D_deltaEta_1_40;   //!
   TBranch        *b_D_deltaPhi_1_40;   //!
   TBranch        *b_D_pxasy_1_40;   //!
   TBranch        *b_D_pyasy_1_40;   //!
   TBranch        *b_D_pzasy_1_40;   //!
   TBranch        *b_D_pasy_1_40;   //!
   TBranch        *b_D_ptasy_1_40;   //!
   TBranch        *b_D_cpx_1_50;   //!
   TBranch        *b_D_cpy_1_50;   //!
   TBranch        *b_D_cpz_1_50;   //!
   TBranch        *b_D_cpt_1_50;   //!
   TBranch        *b_D_cp_1_50;   //!
   TBranch        *b_D_cmult_1_50;   //!
   TBranch        *b_D_deltaEta_1_50;   //!
   TBranch        *b_D_deltaPhi_1_50;   //!
   TBranch        *b_D_pxasy_1_50;   //!
   TBranch        *b_D_pyasy_1_50;   //!
   TBranch        *b_D_pzasy_1_50;   //!
   TBranch        *b_D_pasy_1_50;   //!
   TBranch        *b_D_ptasy_1_50;   //!
   TBranch        *b_D_cpx_1_60;   //!
   TBranch        *b_D_cpy_1_60;   //!
   TBranch        *b_D_cpz_1_60;   //!
   TBranch        *b_D_cpt_1_60;   //!
   TBranch        *b_D_cp_1_60;   //!
   TBranch        *b_D_cmult_1_60;   //!
   TBranch        *b_D_deltaEta_1_60;   //!
   TBranch        *b_D_deltaPhi_1_60;   //!
   TBranch        *b_D_pxasy_1_60;   //!
   TBranch        *b_D_pyasy_1_60;   //!
   TBranch        *b_D_pzasy_1_60;   //!
   TBranch        *b_D_pasy_1_60;   //!
   TBranch        *b_D_ptasy_1_60;   //!
   TBranch        *b_D_cpx_1_70;   //!
   TBranch        *b_D_cpy_1_70;   //!
   TBranch        *b_D_cpz_1_70;   //!
   TBranch        *b_D_cpt_1_70;   //!
   TBranch        *b_D_cp_1_70;   //!
   TBranch        *b_D_cmult_1_70;   //!
   TBranch        *b_D_deltaEta_1_70;   //!
   TBranch        *b_D_deltaPhi_1_70;   //!
   TBranch        *b_D_pxasy_1_70;   //!
   TBranch        *b_D_pyasy_1_70;   //!
   TBranch        *b_D_pzasy_1_70;   //!
   TBranch        *b_D_pasy_1_70;   //!
   TBranch        *b_D_ptasy_1_70;   //!
   TBranch        *b_h0_MINIP;   //!
   TBranch        *b_h0_MINIPCHI2;   //!
   TBranch        *b_h0_MINIPNEXTBEST;   //!
   TBranch        *b_h0_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_h0_OWNPV_X;   //!
   TBranch        *b_h0_OWNPV_Y;   //!
   TBranch        *b_h0_OWNPV_Z;   //!
   TBranch        *b_h0_OWNPV_XERR;   //!
   TBranch        *b_h0_OWNPV_YERR;   //!
   TBranch        *b_h0_OWNPV_ZERR;   //!
   TBranch        *b_h0_OWNPV_CHI2;   //!
   TBranch        *b_h0_OWNPV_NDOF;   //!
   TBranch        *b_h0_OWNPV_COV_;   //!
   TBranch        *b_h0_IP_OWNPV;   //!
   TBranch        *b_h0_IPCHI2_OWNPV;   //!
   TBranch        *b_h0_TOPPV_X;   //!
   TBranch        *b_h0_TOPPV_Y;   //!
   TBranch        *b_h0_TOPPV_Z;   //!
   TBranch        *b_h0_TOPPV_XERR;   //!
   TBranch        *b_h0_TOPPV_YERR;   //!
   TBranch        *b_h0_TOPPV_ZERR;   //!
   TBranch        *b_h0_TOPPV_CHI2;   //!
   TBranch        *b_h0_TOPPV_NDOF;   //!
   TBranch        *b_h0_TOPPV_COV_;   //!
   TBranch        *b_h0_IP_TOPPV;   //!
   TBranch        *b_h0_IPCHI2_TOPPV;   //!
   TBranch        *b_h0_ORIVX_X;   //!
   TBranch        *b_h0_ORIVX_Y;   //!
   TBranch        *b_h0_ORIVX_Z;   //!
   TBranch        *b_h0_ORIVX_XERR;   //!
   TBranch        *b_h0_ORIVX_YERR;   //!
   TBranch        *b_h0_ORIVX_ZERR;   //!
   TBranch        *b_h0_ORIVX_CHI2;   //!
   TBranch        *b_h0_ORIVX_NDOF;   //!
   TBranch        *b_h0_ORIVX_COV_;   //!
   TBranch        *b_h0_IP_ORIVX;   //!
   TBranch        *b_h0_IPCHI2_ORIVX;   //!
   TBranch        *b_h0_P;   //!
   TBranch        *b_h0_PT;   //!
   TBranch        *b_h0_PE;   //!
   TBranch        *b_h0_PX;   //!
   TBranch        *b_h0_PY;   //!
   TBranch        *b_h0_PZ;   //!
   TBranch        *b_h0_M;   //!
   TBranch        *b_h0_L0Calo_HCAL_realET;   //!
   TBranch        *b_h0_L0Calo_HCAL_xProjection;   //!
   TBranch        *b_h0_L0Calo_HCAL_yProjection;   //!
   TBranch        *b_h0_L0Calo_HCAL_region;   //!
   TBranch        *b_h0_L0Calo_HCAL_TriggerET;   //!
   TBranch        *b_h0_L0Calo_HCAL_TriggerHCALET;   //!
   TBranch        *b_h0_L0Calo_HCAL_xTrigger;   //!
   TBranch        *b_h0_L0Calo_HCAL_yTrigger;   //!
   TBranch        *b_h0_ID;   //!
   TBranch        *b_h0_CombDLLMu;   //!
   TBranch        *b_h0_ProbNNmu;   //!
   TBranch        *b_h0_ProbNNghost;   //!
   TBranch        *b_h0_InMuonAcc;   //!
   TBranch        *b_h0_MuonDist2;   //!
   TBranch        *b_h0_regionInM2;   //!
   TBranch        *b_h0_hasMuon;   //!
   TBranch        *b_h0_isMuon;   //!
   TBranch        *b_h0_isMuonLoose;   //!
   TBranch        *b_h0_NShared;   //!
   TBranch        *b_h0_MuonLLmu;   //!
   TBranch        *b_h0_MuonLLbg;   //!
   TBranch        *b_h0_isMuonFromProto;   //!
   TBranch        *b_h0_PIDe;   //!
   TBranch        *b_h0_PIDmu;   //!
   TBranch        *b_h0_PIDK;   //!
   TBranch        *b_h0_PIDp;   //!
   TBranch        *b_h0_ProbNNe;   //!
   TBranch        *b_h0_ProbNNk;   //!
   TBranch        *b_h0_ProbNNp;   //!
   TBranch        *b_h0_ProbNNpi;   //!
   TBranch        *b_h0_hasRich;   //!
   TBranch        *b_h0_hasCalo;   //!
   TBranch        *b_h0_UsedRichAerogel;   //!
   TBranch        *b_h0_UsedRich1Gas;   //!
   TBranch        *b_h0_UsedRich2Gas;   //!
   TBranch        *b_h0_RichAboveElThres;   //!
   TBranch        *b_h0_RichAboveMuThres;   //!
   TBranch        *b_h0_RichAbovePiThres;   //!
   TBranch        *b_h0_RichAboveKaThres;   //!
   TBranch        *b_h0_RichAbovePrThres;   //!
   TBranch        *b_h0_RichDLLe;   //!
   TBranch        *b_h0_RichDLLmu;   //!
   TBranch        *b_h0_RichDLLpi;   //!
   TBranch        *b_h0_RichDLLk;   //!
   TBranch        *b_h0_RichDLLp;   //!
   TBranch        *b_h0_RichDLLbt;   //!
   TBranch        *b_h0_InAccMuon;   //!
   TBranch        *b_h0_MuonMuLL;   //!
   TBranch        *b_h0_MuonBkgLL;   //!
   TBranch        *b_h0_MuonNShared;   //!
   TBranch        *b_h0_InAccEcal;   //!
   TBranch        *b_h0_CaloEcalE;   //!
   TBranch        *b_h0_EcalPIDe;   //!
   TBranch        *b_h0_EcalPIDmu;   //!
   TBranch        *b_h0_InAccHcal;   //!
   TBranch        *b_h0_CaloHcalE;   //!
   TBranch        *b_h0_HcalPIDe;   //!
   TBranch        *b_h0_HcalPIDmu;   //!
   TBranch        *b_h0_InAccPrs;   //!
   TBranch        *b_h0_PrsPIDe;   //!
   TBranch        *b_h0_CaloPrsE;   //!
   TBranch        *b_h0_InAccSpd;   //!
   TBranch        *b_h0_CaloSpdE;   //!
   TBranch        *b_h0_InAccBrem;   //!
   TBranch        *b_h0_BremPIDe;   //!
   TBranch        *b_h0_VeloCharge;   //!
   TBranch        *b_h0_RICHDLLe;   //!
   TBranch        *b_h0_RICHDLLmu;   //!
   TBranch        *b_h0_RICHDLLpi;   //!
   TBranch        *b_h0_RICHDLLK;   //!
   TBranch        *b_h0_RICHDLLp;   //!
   TBranch        *b_h0_RICHDLLbt;   //!
   TBranch        *b_h0_RICHBestID;   //!
   TBranch        *b_h0_RICHThreshold;   //!
   TBranch        *b_h0_RICHThresholdEl;   //!
   TBranch        *b_h0_RICHThresholdMu;   //!
   TBranch        *b_h0_RICHThresholdPi;   //!
   TBranch        *b_h0_RICHThresholdKa;   //!
   TBranch        *b_h0_RICHThresholdPr;   //!
   TBranch        *b_h0_RICHAerogelUsed;   //!
   TBranch        *b_h0_RICH1GasUsed;   //!
   TBranch        *b_h0_RICH2GasUsed;   //!
   TBranch        *b_h0_TRACK_Eta;   //!
   TBranch        *b_h0_TRACK_Phi;   //!
   TBranch        *b_h0_Aerogel_X;   //!
   TBranch        *b_h0_Aerogel_Y;   //!
   TBranch        *b_h0_Aerogel_Z;   //!
   TBranch        *b_h0_Aerogel_Rho;   //!
   TBranch        *b_h0_Aerogel_Phi;   //!
   TBranch        *b_h0_Rich1Gas_X;   //!
   TBranch        *b_h0_Rich1Gas_Y;   //!
   TBranch        *b_h0_Rich1Gas_Z;   //!
   TBranch        *b_h0_Rich1Gas_Rho;   //!
   TBranch        *b_h0_Rich1Gas_Phi;   //!
   TBranch        *b_h0_Rich2Gas_X;   //!
   TBranch        *b_h0_Rich2Gas_Y;   //!
   TBranch        *b_h0_Rich2Gas_Z;   //!
   TBranch        *b_h0_Rich2Gas_Rho;   //!
   TBranch        *b_h0_Rich2Gas_Phi;   //!
   TBranch        *b_h0_L0Global_Dec;   //!
   TBranch        *b_h0_L0Global_TIS;   //!
   TBranch        *b_h0_L0Global_TOS;   //!
   TBranch        *b_h0_Hlt1Global_Dec;   //!
   TBranch        *b_h0_Hlt1Global_TIS;   //!
   TBranch        *b_h0_Hlt1Global_TOS;   //!
   TBranch        *b_h0_Hlt1Phys_Dec;   //!
   TBranch        *b_h0_Hlt1Phys_TIS;   //!
   TBranch        *b_h0_Hlt1Phys_TOS;   //!
   TBranch        *b_h0_Hlt2Global_Dec;   //!
   TBranch        *b_h0_Hlt2Global_TIS;   //!
   TBranch        *b_h0_Hlt2Global_TOS;   //!
   TBranch        *b_h0_Hlt2Phys_Dec;   //!
   TBranch        *b_h0_Hlt2Phys_TIS;   //!
   TBranch        *b_h0_Hlt2Phys_TOS;   //!
   TBranch        *b_h0_TRACK_Type;   //!
   TBranch        *b_h0_TRACK_Key;   //!
   TBranch        *b_h0_TRACK_CHI2;   //!
   TBranch        *b_h0_TRACK_NDOF;   //!
   TBranch        *b_h0_TRACK_CHI2NDOF;   //!
   TBranch        *b_h0_TRACK_PCHI2;   //!
   TBranch        *b_h0_TRACK_VeloCHI2NDOF;   //!
   TBranch        *b_h0_TRACK_TCHI2NDOF;   //!
   TBranch        *b_h0_VELO_UTID;   //!
   TBranch        *b_h0_TRACK_FirstMeasurementX;   //!
   TBranch        *b_h0_TRACK_FirstMeasurementY;   //!
   TBranch        *b_h0_TRACK_FirstMeasurementZ;   //!
   TBranch        *b_h0_TRACK_MatchCHI2;   //!
   TBranch        *b_h0_TRACK_GhostProb;   //!
   TBranch        *b_h0_TRACK_CloneDist;   //!
   TBranch        *b_h0_TRACK_Likelihood;   //!
   TBranch        *b_h0_L0HadronDecision_Dec;   //!
   TBranch        *b_h0_L0HadronDecision_TIS;   //!
   TBranch        *b_h0_L0HadronDecision_TOS;   //!
   TBranch        *b_h0_L0MuonDecision_Dec;   //!
   TBranch        *b_h0_L0MuonDecision_TIS;   //!
   TBranch        *b_h0_L0MuonDecision_TOS;   //!
   TBranch        *b_h0_L0DiMuonDecision_Dec;   //!
   TBranch        *b_h0_L0DiMuonDecision_TIS;   //!
   TBranch        *b_h0_L0DiMuonDecision_TOS;   //!
   TBranch        *b_h0_L0ElectronDecision_Dec;   //!
   TBranch        *b_h0_L0ElectronDecision_TIS;   //!
   TBranch        *b_h0_L0ElectronDecision_TOS;   //!
   TBranch        *b_h0_L0PhotonDecision_Dec;   //!
   TBranch        *b_h0_L0PhotonDecision_TIS;   //!
   TBranch        *b_h0_L0PhotonDecision_TOS;   //!
   TBranch        *b_h0_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_h0_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_h0_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_h0_Hlt1DiMuonLowMassDecision_Dec;   //!
   TBranch        *b_h0_Hlt1DiMuonLowMassDecision_TIS;   //!
   TBranch        *b_h0_Hlt1DiMuonLowMassDecision_TOS;   //!
   TBranch        *b_h0_Hlt1SingleMuonNoIPDecision_Dec;   //!
   TBranch        *b_h0_Hlt1SingleMuonNoIPDecision_TIS;   //!
   TBranch        *b_h0_Hlt1SingleMuonNoIPDecision_TOS;   //!
   TBranch        *b_h0_Hlt1SingleMuonHighPTDecision_Dec;   //!
   TBranch        *b_h0_Hlt1SingleMuonHighPTDecision_TIS;   //!
   TBranch        *b_h0_Hlt1SingleMuonHighPTDecision_TOS;   //!
   TBranch        *b_h0_Hlt1TrackAllL0Decision_Dec;   //!
   TBranch        *b_h0_Hlt1TrackAllL0Decision_TIS;   //!
   TBranch        *b_h0_Hlt1TrackAllL0Decision_TOS;   //!
   TBranch        *b_h0_Hlt1TrackMuonDecision_Dec;   //!
   TBranch        *b_h0_Hlt1TrackMuonDecision_TIS;   //!
   TBranch        *b_h0_Hlt1TrackMuonDecision_TOS;   //!
   TBranch        *b_h0_Hlt1TrackPhotonDecision_Dec;   //!
   TBranch        *b_h0_Hlt1TrackPhotonDecision_TIS;   //!
   TBranch        *b_h0_Hlt1TrackPhotonDecision_TOS;   //!
   TBranch        *b_h0_Hlt1L0AnyDecision_Dec;   //!
   TBranch        *b_h0_Hlt1L0AnyDecision_TIS;   //!
   TBranch        *b_h0_Hlt1L0AnyDecision_TOS;   //!
   TBranch        *b_h0_Hlt1GlobalDecision_Dec;   //!
   TBranch        *b_h0_Hlt1GlobalDecision_TIS;   //!
   TBranch        *b_h0_Hlt1GlobalDecision_TOS;   //!
   TBranch        *b_h0_Hlt2SingleMuonDecision_Dec;   //!
   TBranch        *b_h0_Hlt2SingleMuonDecision_TIS;   //!
   TBranch        *b_h0_Hlt2SingleMuonDecision_TOS;   //!
   TBranch        *b_h0_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_h0_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_h0_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD2HMuMuDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD2HMuMuDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD2HMuMuDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD02KKMuMuDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD02KKMuMuDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD02KKMuMuDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_K3piDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_K3piDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_K3piDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_4piDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_4piDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_4piDecision_TOS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;   //!
   TBranch        *b_h0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;   //!
   TBranch        *b_h1_MINIP;   //!
   TBranch        *b_h1_MINIPCHI2;   //!
   TBranch        *b_h1_MINIPNEXTBEST;   //!
   TBranch        *b_h1_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_h1_OWNPV_X;   //!
   TBranch        *b_h1_OWNPV_Y;   //!
   TBranch        *b_h1_OWNPV_Z;   //!
   TBranch        *b_h1_OWNPV_XERR;   //!
   TBranch        *b_h1_OWNPV_YERR;   //!
   TBranch        *b_h1_OWNPV_ZERR;   //!
   TBranch        *b_h1_OWNPV_CHI2;   //!
   TBranch        *b_h1_OWNPV_NDOF;   //!
   TBranch        *b_h1_OWNPV_COV_;   //!
   TBranch        *b_h1_IP_OWNPV;   //!
   TBranch        *b_h1_IPCHI2_OWNPV;   //!
   TBranch        *b_h1_TOPPV_X;   //!
   TBranch        *b_h1_TOPPV_Y;   //!
   TBranch        *b_h1_TOPPV_Z;   //!
   TBranch        *b_h1_TOPPV_XERR;   //!
   TBranch        *b_h1_TOPPV_YERR;   //!
   TBranch        *b_h1_TOPPV_ZERR;   //!
   TBranch        *b_h1_TOPPV_CHI2;   //!
   TBranch        *b_h1_TOPPV_NDOF;   //!
   TBranch        *b_h1_TOPPV_COV_;   //!
   TBranch        *b_h1_IP_TOPPV;   //!
   TBranch        *b_h1_IPCHI2_TOPPV;   //!
   TBranch        *b_h1_ORIVX_X;   //!
   TBranch        *b_h1_ORIVX_Y;   //!
   TBranch        *b_h1_ORIVX_Z;   //!
   TBranch        *b_h1_ORIVX_XERR;   //!
   TBranch        *b_h1_ORIVX_YERR;   //!
   TBranch        *b_h1_ORIVX_ZERR;   //!
   TBranch        *b_h1_ORIVX_CHI2;   //!
   TBranch        *b_h1_ORIVX_NDOF;   //!
   TBranch        *b_h1_ORIVX_COV_;   //!
   TBranch        *b_h1_IP_ORIVX;   //!
   TBranch        *b_h1_IPCHI2_ORIVX;   //!
   TBranch        *b_h1_P;   //!
   TBranch        *b_h1_PT;   //!
   TBranch        *b_h1_PE;   //!
   TBranch        *b_h1_PX;   //!
   TBranch        *b_h1_PY;   //!
   TBranch        *b_h1_PZ;   //!
   TBranch        *b_h1_M;   //!
   TBranch        *b_h1_L0Calo_HCAL_realET;   //!
   TBranch        *b_h1_L0Calo_HCAL_xProjection;   //!
   TBranch        *b_h1_L0Calo_HCAL_yProjection;   //!
   TBranch        *b_h1_L0Calo_HCAL_region;   //!
   TBranch        *b_h1_L0Calo_HCAL_TriggerET;   //!
   TBranch        *b_h1_L0Calo_HCAL_TriggerHCALET;   //!
   TBranch        *b_h1_L0Calo_HCAL_xTrigger;   //!
   TBranch        *b_h1_L0Calo_HCAL_yTrigger;   //!
   TBranch        *b_h1_ID;   //!
   TBranch        *b_h1_CombDLLMu;   //!
   TBranch        *b_h1_ProbNNmu;   //!
   TBranch        *b_h1_ProbNNghost;   //!
   TBranch        *b_h1_InMuonAcc;   //!
   TBranch        *b_h1_MuonDist2;   //!
   TBranch        *b_h1_regionInM2;   //!
   TBranch        *b_h1_hasMuon;   //!
   TBranch        *b_h1_isMuon;   //!
   TBranch        *b_h1_isMuonLoose;   //!
   TBranch        *b_h1_NShared;   //!
   TBranch        *b_h1_MuonLLmu;   //!
   TBranch        *b_h1_MuonLLbg;   //!
   TBranch        *b_h1_isMuonFromProto;   //!
   TBranch        *b_h1_PIDe;   //!
   TBranch        *b_h1_PIDmu;   //!
   TBranch        *b_h1_PIDK;   //!
   TBranch        *b_h1_PIDp;   //!
   TBranch        *b_h1_ProbNNe;   //!
   TBranch        *b_h1_ProbNNk;   //!
   TBranch        *b_h1_ProbNNp;   //!
   TBranch        *b_h1_ProbNNpi;   //!
   TBranch        *b_h1_hasRich;   //!
   TBranch        *b_h1_hasCalo;   //!
   TBranch        *b_h1_UsedRichAerogel;   //!
   TBranch        *b_h1_UsedRich1Gas;   //!
   TBranch        *b_h1_UsedRich2Gas;   //!
   TBranch        *b_h1_RichAboveElThres;   //!
   TBranch        *b_h1_RichAboveMuThres;   //!
   TBranch        *b_h1_RichAbovePiThres;   //!
   TBranch        *b_h1_RichAboveKaThres;   //!
   TBranch        *b_h1_RichAbovePrThres;   //!
   TBranch        *b_h1_RichDLLe;   //!
   TBranch        *b_h1_RichDLLmu;   //!
   TBranch        *b_h1_RichDLLpi;   //!
   TBranch        *b_h1_RichDLLk;   //!
   TBranch        *b_h1_RichDLLp;   //!
   TBranch        *b_h1_RichDLLbt;   //!
   TBranch        *b_h1_InAccMuon;   //!
   TBranch        *b_h1_MuonMuLL;   //!
   TBranch        *b_h1_MuonBkgLL;   //!
   TBranch        *b_h1_MuonNShared;   //!
   TBranch        *b_h1_InAccEcal;   //!
   TBranch        *b_h1_CaloEcalE;   //!
   TBranch        *b_h1_EcalPIDe;   //!
   TBranch        *b_h1_EcalPIDmu;   //!
   TBranch        *b_h1_InAccHcal;   //!
   TBranch        *b_h1_CaloHcalE;   //!
   TBranch        *b_h1_HcalPIDe;   //!
   TBranch        *b_h1_HcalPIDmu;   //!
   TBranch        *b_h1_InAccPrs;   //!
   TBranch        *b_h1_PrsPIDe;   //!
   TBranch        *b_h1_CaloPrsE;   //!
   TBranch        *b_h1_InAccSpd;   //!
   TBranch        *b_h1_CaloSpdE;   //!
   TBranch        *b_h1_InAccBrem;   //!
   TBranch        *b_h1_BremPIDe;   //!
   TBranch        *b_h1_VeloCharge;   //!
   TBranch        *b_h1_RICHDLLe;   //!
   TBranch        *b_h1_RICHDLLmu;   //!
   TBranch        *b_h1_RICHDLLpi;   //!
   TBranch        *b_h1_RICHDLLK;   //!
   TBranch        *b_h1_RICHDLLp;   //!
   TBranch        *b_h1_RICHDLLbt;   //!
   TBranch        *b_h1_RICHBestID;   //!
   TBranch        *b_h1_RICHThreshold;   //!
   TBranch        *b_h1_RICHThresholdEl;   //!
   TBranch        *b_h1_RICHThresholdMu;   //!
   TBranch        *b_h1_RICHThresholdPi;   //!
   TBranch        *b_h1_RICHThresholdKa;   //!
   TBranch        *b_h1_RICHThresholdPr;   //!
   TBranch        *b_h1_RICHAerogelUsed;   //!
   TBranch        *b_h1_RICH1GasUsed;   //!
   TBranch        *b_h1_RICH2GasUsed;   //!
   TBranch        *b_h1_TRACK_Eta;   //!
   TBranch        *b_h1_TRACK_Phi;   //!
   TBranch        *b_h1_Aerogel_X;   //!
   TBranch        *b_h1_Aerogel_Y;   //!
   TBranch        *b_h1_Aerogel_Z;   //!
   TBranch        *b_h1_Aerogel_Rho;   //!
   TBranch        *b_h1_Aerogel_Phi;   //!
   TBranch        *b_h1_Rich1Gas_X;   //!
   TBranch        *b_h1_Rich1Gas_Y;   //!
   TBranch        *b_h1_Rich1Gas_Z;   //!
   TBranch        *b_h1_Rich1Gas_Rho;   //!
   TBranch        *b_h1_Rich1Gas_Phi;   //!
   TBranch        *b_h1_Rich2Gas_X;   //!
   TBranch        *b_h1_Rich2Gas_Y;   //!
   TBranch        *b_h1_Rich2Gas_Z;   //!
   TBranch        *b_h1_Rich2Gas_Rho;   //!
   TBranch        *b_h1_Rich2Gas_Phi;   //!
   TBranch        *b_h1_L0Global_Dec;   //!
   TBranch        *b_h1_L0Global_TIS;   //!
   TBranch        *b_h1_L0Global_TOS;   //!
   TBranch        *b_h1_Hlt1Global_Dec;   //!
   TBranch        *b_h1_Hlt1Global_TIS;   //!
   TBranch        *b_h1_Hlt1Global_TOS;   //!
   TBranch        *b_h1_Hlt1Phys_Dec;   //!
   TBranch        *b_h1_Hlt1Phys_TIS;   //!
   TBranch        *b_h1_Hlt1Phys_TOS;   //!
   TBranch        *b_h1_Hlt2Global_Dec;   //!
   TBranch        *b_h1_Hlt2Global_TIS;   //!
   TBranch        *b_h1_Hlt2Global_TOS;   //!
   TBranch        *b_h1_Hlt2Phys_Dec;   //!
   TBranch        *b_h1_Hlt2Phys_TIS;   //!
   TBranch        *b_h1_Hlt2Phys_TOS;   //!
   TBranch        *b_h1_TRACK_Type;   //!
   TBranch        *b_h1_TRACK_Key;   //!
   TBranch        *b_h1_TRACK_CHI2;   //!
   TBranch        *b_h1_TRACK_NDOF;   //!
   TBranch        *b_h1_TRACK_CHI2NDOF;   //!
   TBranch        *b_h1_TRACK_PCHI2;   //!
   TBranch        *b_h1_TRACK_VeloCHI2NDOF;   //!
   TBranch        *b_h1_TRACK_TCHI2NDOF;   //!
   TBranch        *b_h1_VELO_UTID;   //!
   TBranch        *b_h1_TRACK_FirstMeasurementX;   //!
   TBranch        *b_h1_TRACK_FirstMeasurementY;   //!
   TBranch        *b_h1_TRACK_FirstMeasurementZ;   //!
   TBranch        *b_h1_TRACK_MatchCHI2;   //!
   TBranch        *b_h1_TRACK_GhostProb;   //!
   TBranch        *b_h1_TRACK_CloneDist;   //!
   TBranch        *b_h1_TRACK_Likelihood;   //!
   TBranch        *b_h1_L0HadronDecision_Dec;   //!
   TBranch        *b_h1_L0HadronDecision_TIS;   //!
   TBranch        *b_h1_L0HadronDecision_TOS;   //!
   TBranch        *b_h1_L0MuonDecision_Dec;   //!
   TBranch        *b_h1_L0MuonDecision_TIS;   //!
   TBranch        *b_h1_L0MuonDecision_TOS;   //!
   TBranch        *b_h1_L0DiMuonDecision_Dec;   //!
   TBranch        *b_h1_L0DiMuonDecision_TIS;   //!
   TBranch        *b_h1_L0DiMuonDecision_TOS;   //!
   TBranch        *b_h1_L0ElectronDecision_Dec;   //!
   TBranch        *b_h1_L0ElectronDecision_TIS;   //!
   TBranch        *b_h1_L0ElectronDecision_TOS;   //!
   TBranch        *b_h1_L0PhotonDecision_Dec;   //!
   TBranch        *b_h1_L0PhotonDecision_TIS;   //!
   TBranch        *b_h1_L0PhotonDecision_TOS;   //!
   TBranch        *b_h1_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_h1_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_h1_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_h1_Hlt1DiMuonLowMassDecision_Dec;   //!
   TBranch        *b_h1_Hlt1DiMuonLowMassDecision_TIS;   //!
   TBranch        *b_h1_Hlt1DiMuonLowMassDecision_TOS;   //!
   TBranch        *b_h1_Hlt1SingleMuonNoIPDecision_Dec;   //!
   TBranch        *b_h1_Hlt1SingleMuonNoIPDecision_TIS;   //!
   TBranch        *b_h1_Hlt1SingleMuonNoIPDecision_TOS;   //!
   TBranch        *b_h1_Hlt1SingleMuonHighPTDecision_Dec;   //!
   TBranch        *b_h1_Hlt1SingleMuonHighPTDecision_TIS;   //!
   TBranch        *b_h1_Hlt1SingleMuonHighPTDecision_TOS;   //!
   TBranch        *b_h1_Hlt1TrackAllL0Decision_Dec;   //!
   TBranch        *b_h1_Hlt1TrackAllL0Decision_TIS;   //!
   TBranch        *b_h1_Hlt1TrackAllL0Decision_TOS;   //!
   TBranch        *b_h1_Hlt1TrackMuonDecision_Dec;   //!
   TBranch        *b_h1_Hlt1TrackMuonDecision_TIS;   //!
   TBranch        *b_h1_Hlt1TrackMuonDecision_TOS;   //!
   TBranch        *b_h1_Hlt1TrackPhotonDecision_Dec;   //!
   TBranch        *b_h1_Hlt1TrackPhotonDecision_TIS;   //!
   TBranch        *b_h1_Hlt1TrackPhotonDecision_TOS;   //!
   TBranch        *b_h1_Hlt1L0AnyDecision_Dec;   //!
   TBranch        *b_h1_Hlt1L0AnyDecision_TIS;   //!
   TBranch        *b_h1_Hlt1L0AnyDecision_TOS;   //!
   TBranch        *b_h1_Hlt1GlobalDecision_Dec;   //!
   TBranch        *b_h1_Hlt1GlobalDecision_TIS;   //!
   TBranch        *b_h1_Hlt1GlobalDecision_TOS;   //!
   TBranch        *b_h1_Hlt2SingleMuonDecision_Dec;   //!
   TBranch        *b_h1_Hlt2SingleMuonDecision_TIS;   //!
   TBranch        *b_h1_Hlt2SingleMuonDecision_TOS;   //!
   TBranch        *b_h1_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_h1_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_h1_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD2HMuMuDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD2HMuMuDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD2HMuMuDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD02KKMuMuDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD02KKMuMuDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD02KKMuMuDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_K3piDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_K3piDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_K3piDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_4piDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_4piDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_4piDecision_TOS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;   //!
   TBranch        *b_h1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;   //!
   TBranch        *b_mu0_MINIP;   //!
   TBranch        *b_mu0_MINIPCHI2;   //!
   TBranch        *b_mu0_MINIPNEXTBEST;   //!
   TBranch        *b_mu0_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_mu0_OWNPV_X;   //!
   TBranch        *b_mu0_OWNPV_Y;   //!
   TBranch        *b_mu0_OWNPV_Z;   //!
   TBranch        *b_mu0_OWNPV_XERR;   //!
   TBranch        *b_mu0_OWNPV_YERR;   //!
   TBranch        *b_mu0_OWNPV_ZERR;   //!
   TBranch        *b_mu0_OWNPV_CHI2;   //!
   TBranch        *b_mu0_OWNPV_NDOF;   //!
   TBranch        *b_mu0_OWNPV_COV_;   //!
   TBranch        *b_mu0_IP_OWNPV;   //!
   TBranch        *b_mu0_IPCHI2_OWNPV;   //!
   TBranch        *b_mu0_TOPPV_X;   //!
   TBranch        *b_mu0_TOPPV_Y;   //!
   TBranch        *b_mu0_TOPPV_Z;   //!
   TBranch        *b_mu0_TOPPV_XERR;   //!
   TBranch        *b_mu0_TOPPV_YERR;   //!
   TBranch        *b_mu0_TOPPV_ZERR;   //!
   TBranch        *b_mu0_TOPPV_CHI2;   //!
   TBranch        *b_mu0_TOPPV_NDOF;   //!
   TBranch        *b_mu0_TOPPV_COV_;   //!
   TBranch        *b_mu0_IP_TOPPV;   //!
   TBranch        *b_mu0_IPCHI2_TOPPV;   //!
   TBranch        *b_mu0_ORIVX_X;   //!
   TBranch        *b_mu0_ORIVX_Y;   //!
   TBranch        *b_mu0_ORIVX_Z;   //!
   TBranch        *b_mu0_ORIVX_XERR;   //!
   TBranch        *b_mu0_ORIVX_YERR;   //!
   TBranch        *b_mu0_ORIVX_ZERR;   //!
   TBranch        *b_mu0_ORIVX_CHI2;   //!
   TBranch        *b_mu0_ORIVX_NDOF;   //!
   TBranch        *b_mu0_ORIVX_COV_;   //!
   TBranch        *b_mu0_IP_ORIVX;   //!
   TBranch        *b_mu0_IPCHI2_ORIVX;   //!
   TBranch        *b_mu0_P;   //!
   TBranch        *b_mu0_PT;   //!
   TBranch        *b_mu0_PE;   //!
   TBranch        *b_mu0_PX;   //!
   TBranch        *b_mu0_PY;   //!
   TBranch        *b_mu0_PZ;   //!
   TBranch        *b_mu0_M;   //!
   TBranch        *b_mu0_L0Calo_HCAL_realET;   //!
   TBranch        *b_mu0_L0Calo_HCAL_xProjection;   //!
   TBranch        *b_mu0_L0Calo_HCAL_yProjection;   //!
   TBranch        *b_mu0_L0Calo_HCAL_region;   //!
   TBranch        *b_mu0_L0Calo_HCAL_TriggerET;   //!
   TBranch        *b_mu0_L0Calo_HCAL_TriggerHCALET;   //!
   TBranch        *b_mu0_L0Calo_HCAL_xTrigger;   //!
   TBranch        *b_mu0_L0Calo_HCAL_yTrigger;   //!
   TBranch        *b_mu0_ID;   //!
   TBranch        *b_mu0_CombDLLMu;   //!
   TBranch        *b_mu0_ProbNNmu;   //!
   TBranch        *b_mu0_ProbNNghost;   //!
   TBranch        *b_mu0_InMuonAcc;   //!
   TBranch        *b_mu0_MuonDist2;   //!
   TBranch        *b_mu0_regionInM2;   //!
   TBranch        *b_mu0_hasMuon;   //!
   TBranch        *b_mu0_isMuon;   //!
   TBranch        *b_mu0_isMuonLoose;   //!
   TBranch        *b_mu0_NShared;   //!
   TBranch        *b_mu0_MuonLLmu;   //!
   TBranch        *b_mu0_MuonLLbg;   //!
   TBranch        *b_mu0_isMuonFromProto;   //!
   TBranch        *b_mu0_PIDe;   //!
   TBranch        *b_mu0_PIDmu;   //!
   TBranch        *b_mu0_PIDK;   //!
   TBranch        *b_mu0_PIDp;   //!
   TBranch        *b_mu0_ProbNNe;   //!
   TBranch        *b_mu0_ProbNNk;   //!
   TBranch        *b_mu0_ProbNNp;   //!
   TBranch        *b_mu0_ProbNNpi;   //!
   TBranch        *b_mu0_hasRich;   //!
   TBranch        *b_mu0_hasCalo;   //!
   TBranch        *b_mu0_UsedRichAerogel;   //!
   TBranch        *b_mu0_UsedRich1Gas;   //!
   TBranch        *b_mu0_UsedRich2Gas;   //!
   TBranch        *b_mu0_RichAboveElThres;   //!
   TBranch        *b_mu0_RichAboveMuThres;   //!
   TBranch        *b_mu0_RichAbovePiThres;   //!
   TBranch        *b_mu0_RichAboveKaThres;   //!
   TBranch        *b_mu0_RichAbovePrThres;   //!
   TBranch        *b_mu0_RichDLLe;   //!
   TBranch        *b_mu0_RichDLLmu;   //!
   TBranch        *b_mu0_RichDLLpi;   //!
   TBranch        *b_mu0_RichDLLk;   //!
   TBranch        *b_mu0_RichDLLp;   //!
   TBranch        *b_mu0_RichDLLbt;   //!
   TBranch        *b_mu0_InAccMuon;   //!
   TBranch        *b_mu0_MuonMuLL;   //!
   TBranch        *b_mu0_MuonBkgLL;   //!
   TBranch        *b_mu0_MuonNShared;   //!
   TBranch        *b_mu0_InAccEcal;   //!
   TBranch        *b_mu0_CaloEcalE;   //!
   TBranch        *b_mu0_EcalPIDe;   //!
   TBranch        *b_mu0_EcalPIDmu;   //!
   TBranch        *b_mu0_InAccHcal;   //!
   TBranch        *b_mu0_CaloHcalE;   //!
   TBranch        *b_mu0_HcalPIDe;   //!
   TBranch        *b_mu0_HcalPIDmu;   //!
   TBranch        *b_mu0_InAccPrs;   //!
   TBranch        *b_mu0_PrsPIDe;   //!
   TBranch        *b_mu0_CaloPrsE;   //!
   TBranch        *b_mu0_InAccSpd;   //!
   TBranch        *b_mu0_CaloSpdE;   //!
   TBranch        *b_mu0_InAccBrem;   //!
   TBranch        *b_mu0_BremPIDe;   //!
   TBranch        *b_mu0_VeloCharge;   //!
   TBranch        *b_mu0_RICHDLLe;   //!
   TBranch        *b_mu0_RICHDLLmu;   //!
   TBranch        *b_mu0_RICHDLLpi;   //!
   TBranch        *b_mu0_RICHDLLK;   //!
   TBranch        *b_mu0_RICHDLLp;   //!
   TBranch        *b_mu0_RICHDLLbt;   //!
   TBranch        *b_mu0_RICHBestID;   //!
   TBranch        *b_mu0_RICHThreshold;   //!
   TBranch        *b_mu0_RICHThresholdEl;   //!
   TBranch        *b_mu0_RICHThresholdMu;   //!
   TBranch        *b_mu0_RICHThresholdPi;   //!
   TBranch        *b_mu0_RICHThresholdKa;   //!
   TBranch        *b_mu0_RICHThresholdPr;   //!
   TBranch        *b_mu0_RICHAerogelUsed;   //!
   TBranch        *b_mu0_RICH1GasUsed;   //!
   TBranch        *b_mu0_RICH2GasUsed;   //!
   TBranch        *b_mu0_TRACK_Eta;   //!
   TBranch        *b_mu0_TRACK_Phi;   //!
   TBranch        *b_mu0_Aerogel_X;   //!
   TBranch        *b_mu0_Aerogel_Y;   //!
   TBranch        *b_mu0_Aerogel_Z;   //!
   TBranch        *b_mu0_Aerogel_Rho;   //!
   TBranch        *b_mu0_Aerogel_Phi;   //!
   TBranch        *b_mu0_Rich1Gas_X;   //!
   TBranch        *b_mu0_Rich1Gas_Y;   //!
   TBranch        *b_mu0_Rich1Gas_Z;   //!
   TBranch        *b_mu0_Rich1Gas_Rho;   //!
   TBranch        *b_mu0_Rich1Gas_Phi;   //!
   TBranch        *b_mu0_Rich2Gas_X;   //!
   TBranch        *b_mu0_Rich2Gas_Y;   //!
   TBranch        *b_mu0_Rich2Gas_Z;   //!
   TBranch        *b_mu0_Rich2Gas_Rho;   //!
   TBranch        *b_mu0_Rich2Gas_Phi;   //!
   TBranch        *b_mu0_L0Global_Dec;   //!
   TBranch        *b_mu0_L0Global_TIS;   //!
   TBranch        *b_mu0_L0Global_TOS;   //!
   TBranch        *b_mu0_Hlt1Global_Dec;   //!
   TBranch        *b_mu0_Hlt1Global_TIS;   //!
   TBranch        *b_mu0_Hlt1Global_TOS;   //!
   TBranch        *b_mu0_Hlt1Phys_Dec;   //!
   TBranch        *b_mu0_Hlt1Phys_TIS;   //!
   TBranch        *b_mu0_Hlt1Phys_TOS;   //!
   TBranch        *b_mu0_Hlt2Global_Dec;   //!
   TBranch        *b_mu0_Hlt2Global_TIS;   //!
   TBranch        *b_mu0_Hlt2Global_TOS;   //!
   TBranch        *b_mu0_Hlt2Phys_Dec;   //!
   TBranch        *b_mu0_Hlt2Phys_TIS;   //!
   TBranch        *b_mu0_Hlt2Phys_TOS;   //!
   TBranch        *b_mu0_TRACK_Type;   //!
   TBranch        *b_mu0_TRACK_Key;   //!
   TBranch        *b_mu0_TRACK_CHI2;   //!
   TBranch        *b_mu0_TRACK_NDOF;   //!
   TBranch        *b_mu0_TRACK_CHI2NDOF;   //!
   TBranch        *b_mu0_TRACK_PCHI2;   //!
   TBranch        *b_mu0_TRACK_VeloCHI2NDOF;   //!
   TBranch        *b_mu0_TRACK_TCHI2NDOF;   //!
   TBranch        *b_mu0_VELO_UTID;   //!
   TBranch        *b_mu0_TRACK_FirstMeasurementX;   //!
   TBranch        *b_mu0_TRACK_FirstMeasurementY;   //!
   TBranch        *b_mu0_TRACK_FirstMeasurementZ;   //!
   TBranch        *b_mu0_TRACK_MatchCHI2;   //!
   TBranch        *b_mu0_TRACK_GhostProb;   //!
   TBranch        *b_mu0_TRACK_CloneDist;   //!
   TBranch        *b_mu0_TRACK_Likelihood;   //!
   TBranch        *b_mu0_L0HadronDecision_Dec;   //!
   TBranch        *b_mu0_L0HadronDecision_TIS;   //!
   TBranch        *b_mu0_L0HadronDecision_TOS;   //!
   TBranch        *b_mu0_L0MuonDecision_Dec;   //!
   TBranch        *b_mu0_L0MuonDecision_TIS;   //!
   TBranch        *b_mu0_L0MuonDecision_TOS;   //!
   TBranch        *b_mu0_L0DiMuonDecision_Dec;   //!
   TBranch        *b_mu0_L0DiMuonDecision_TIS;   //!
   TBranch        *b_mu0_L0DiMuonDecision_TOS;   //!
   TBranch        *b_mu0_L0ElectronDecision_Dec;   //!
   TBranch        *b_mu0_L0ElectronDecision_TIS;   //!
   TBranch        *b_mu0_L0ElectronDecision_TOS;   //!
   TBranch        *b_mu0_L0PhotonDecision_Dec;   //!
   TBranch        *b_mu0_L0PhotonDecision_TIS;   //!
   TBranch        *b_mu0_L0PhotonDecision_TOS;   //!
   TBranch        *b_mu0_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_mu0_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_mu0_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_mu0_Hlt1DiMuonLowMassDecision_Dec;   //!
   TBranch        *b_mu0_Hlt1DiMuonLowMassDecision_TIS;   //!
   TBranch        *b_mu0_Hlt1DiMuonLowMassDecision_TOS;   //!
   TBranch        *b_mu0_Hlt1SingleMuonNoIPDecision_Dec;   //!
   TBranch        *b_mu0_Hlt1SingleMuonNoIPDecision_TIS;   //!
   TBranch        *b_mu0_Hlt1SingleMuonNoIPDecision_TOS;   //!
   TBranch        *b_mu0_Hlt1SingleMuonHighPTDecision_Dec;   //!
   TBranch        *b_mu0_Hlt1SingleMuonHighPTDecision_TIS;   //!
   TBranch        *b_mu0_Hlt1SingleMuonHighPTDecision_TOS;   //!
   TBranch        *b_mu0_Hlt1TrackAllL0Decision_Dec;   //!
   TBranch        *b_mu0_Hlt1TrackAllL0Decision_TIS;   //!
   TBranch        *b_mu0_Hlt1TrackAllL0Decision_TOS;   //!
   TBranch        *b_mu0_Hlt1TrackMuonDecision_Dec;   //!
   TBranch        *b_mu0_Hlt1TrackMuonDecision_TIS;   //!
   TBranch        *b_mu0_Hlt1TrackMuonDecision_TOS;   //!
   TBranch        *b_mu0_Hlt1TrackPhotonDecision_Dec;   //!
   TBranch        *b_mu0_Hlt1TrackPhotonDecision_TIS;   //!
   TBranch        *b_mu0_Hlt1TrackPhotonDecision_TOS;   //!
   TBranch        *b_mu0_Hlt1L0AnyDecision_Dec;   //!
   TBranch        *b_mu0_Hlt1L0AnyDecision_TIS;   //!
   TBranch        *b_mu0_Hlt1L0AnyDecision_TOS;   //!
   TBranch        *b_mu0_Hlt1GlobalDecision_Dec;   //!
   TBranch        *b_mu0_Hlt1GlobalDecision_TIS;   //!
   TBranch        *b_mu0_Hlt1GlobalDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2SingleMuonDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2SingleMuonDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2SingleMuonDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD2HMuMuDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD2HMuMuDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD2HMuMuDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD02KKMuMuDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD02KKMuMuDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD02KKMuMuDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_K3piDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_K3piDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_K3piDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_4piDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_4piDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_4piDecision_TOS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;   //!
   TBranch        *b_mu0_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;   //!
   TBranch        *b_mu1_MINIP;   //!
   TBranch        *b_mu1_MINIPCHI2;   //!
   TBranch        *b_mu1_MINIPNEXTBEST;   //!
   TBranch        *b_mu1_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_mu1_OWNPV_X;   //!
   TBranch        *b_mu1_OWNPV_Y;   //!
   TBranch        *b_mu1_OWNPV_Z;   //!
   TBranch        *b_mu1_OWNPV_XERR;   //!
   TBranch        *b_mu1_OWNPV_YERR;   //!
   TBranch        *b_mu1_OWNPV_ZERR;   //!
   TBranch        *b_mu1_OWNPV_CHI2;   //!
   TBranch        *b_mu1_OWNPV_NDOF;   //!
   TBranch        *b_mu1_OWNPV_COV_;   //!
   TBranch        *b_mu1_IP_OWNPV;   //!
   TBranch        *b_mu1_IPCHI2_OWNPV;   //!
   TBranch        *b_mu1_TOPPV_X;   //!
   TBranch        *b_mu1_TOPPV_Y;   //!
   TBranch        *b_mu1_TOPPV_Z;   //!
   TBranch        *b_mu1_TOPPV_XERR;   //!
   TBranch        *b_mu1_TOPPV_YERR;   //!
   TBranch        *b_mu1_TOPPV_ZERR;   //!
   TBranch        *b_mu1_TOPPV_CHI2;   //!
   TBranch        *b_mu1_TOPPV_NDOF;   //!
   TBranch        *b_mu1_TOPPV_COV_;   //!
   TBranch        *b_mu1_IP_TOPPV;   //!
   TBranch        *b_mu1_IPCHI2_TOPPV;   //!
   TBranch        *b_mu1_ORIVX_X;   //!
   TBranch        *b_mu1_ORIVX_Y;   //!
   TBranch        *b_mu1_ORIVX_Z;   //!
   TBranch        *b_mu1_ORIVX_XERR;   //!
   TBranch        *b_mu1_ORIVX_YERR;   //!
   TBranch        *b_mu1_ORIVX_ZERR;   //!
   TBranch        *b_mu1_ORIVX_CHI2;   //!
   TBranch        *b_mu1_ORIVX_NDOF;   //!
   TBranch        *b_mu1_ORIVX_COV_;   //!
   TBranch        *b_mu1_IP_ORIVX;   //!
   TBranch        *b_mu1_IPCHI2_ORIVX;   //!
   TBranch        *b_mu1_P;   //!
   TBranch        *b_mu1_PT;   //!
   TBranch        *b_mu1_PE;   //!
   TBranch        *b_mu1_PX;   //!
   TBranch        *b_mu1_PY;   //!
   TBranch        *b_mu1_PZ;   //!
   TBranch        *b_mu1_M;   //!
   TBranch        *b_mu1_L0Calo_HCAL_realET;   //!
   TBranch        *b_mu1_L0Calo_HCAL_xProjection;   //!
   TBranch        *b_mu1_L0Calo_HCAL_yProjection;   //!
   TBranch        *b_mu1_L0Calo_HCAL_region;   //!
   TBranch        *b_mu1_L0Calo_HCAL_TriggerET;   //!
   TBranch        *b_mu1_L0Calo_HCAL_TriggerHCALET;   //!
   TBranch        *b_mu1_L0Calo_HCAL_xTrigger;   //!
   TBranch        *b_mu1_L0Calo_HCAL_yTrigger;   //!
   TBranch        *b_mu1_ID;   //!
   TBranch        *b_mu1_CombDLLMu;   //!
   TBranch        *b_mu1_ProbNNmu;   //!
   TBranch        *b_mu1_ProbNNghost;   //!
   TBranch        *b_mu1_InMuonAcc;   //!
   TBranch        *b_mu1_MuonDist2;   //!
   TBranch        *b_mu1_regionInM2;   //!
   TBranch        *b_mu1_hasMuon;   //!
   TBranch        *b_mu1_isMuon;   //!
   TBranch        *b_mu1_isMuonLoose;   //!
   TBranch        *b_mu1_NShared;   //!
   TBranch        *b_mu1_MuonLLmu;   //!
   TBranch        *b_mu1_MuonLLbg;   //!
   TBranch        *b_mu1_isMuonFromProto;   //!
   TBranch        *b_mu1_PIDe;   //!
   TBranch        *b_mu1_PIDmu;   //!
   TBranch        *b_mu1_PIDK;   //!
   TBranch        *b_mu1_PIDp;   //!
   TBranch        *b_mu1_ProbNNe;   //!
   TBranch        *b_mu1_ProbNNk;   //!
   TBranch        *b_mu1_ProbNNp;   //!
   TBranch        *b_mu1_ProbNNpi;   //!
   TBranch        *b_mu1_hasRich;   //!
   TBranch        *b_mu1_hasCalo;   //!
   TBranch        *b_mu1_UsedRichAerogel;   //!
   TBranch        *b_mu1_UsedRich1Gas;   //!
   TBranch        *b_mu1_UsedRich2Gas;   //!
   TBranch        *b_mu1_RichAboveElThres;   //!
   TBranch        *b_mu1_RichAboveMuThres;   //!
   TBranch        *b_mu1_RichAbovePiThres;   //!
   TBranch        *b_mu1_RichAboveKaThres;   //!
   TBranch        *b_mu1_RichAbovePrThres;   //!
   TBranch        *b_mu1_RichDLLe;   //!
   TBranch        *b_mu1_RichDLLmu;   //!
   TBranch        *b_mu1_RichDLLpi;   //!
   TBranch        *b_mu1_RichDLLk;   //!
   TBranch        *b_mu1_RichDLLp;   //!
   TBranch        *b_mu1_RichDLLbt;   //!
   TBranch        *b_mu1_InAccMuon;   //!
   TBranch        *b_mu1_MuonMuLL;   //!
   TBranch        *b_mu1_MuonBkgLL;   //!
   TBranch        *b_mu1_MuonNShared;   //!
   TBranch        *b_mu1_InAccEcal;   //!
   TBranch        *b_mu1_CaloEcalE;   //!
   TBranch        *b_mu1_EcalPIDe;   //!
   TBranch        *b_mu1_EcalPIDmu;   //!
   TBranch        *b_mu1_InAccHcal;   //!
   TBranch        *b_mu1_CaloHcalE;   //!
   TBranch        *b_mu1_HcalPIDe;   //!
   TBranch        *b_mu1_HcalPIDmu;   //!
   TBranch        *b_mu1_InAccPrs;   //!
   TBranch        *b_mu1_PrsPIDe;   //!
   TBranch        *b_mu1_CaloPrsE;   //!
   TBranch        *b_mu1_InAccSpd;   //!
   TBranch        *b_mu1_CaloSpdE;   //!
   TBranch        *b_mu1_InAccBrem;   //!
   TBranch        *b_mu1_BremPIDe;   //!
   TBranch        *b_mu1_VeloCharge;   //!
   TBranch        *b_mu1_RICHDLLe;   //!
   TBranch        *b_mu1_RICHDLLmu;   //!
   TBranch        *b_mu1_RICHDLLpi;   //!
   TBranch        *b_mu1_RICHDLLK;   //!
   TBranch        *b_mu1_RICHDLLp;   //!
   TBranch        *b_mu1_RICHDLLbt;   //!
   TBranch        *b_mu1_RICHBestID;   //!
   TBranch        *b_mu1_RICHThreshold;   //!
   TBranch        *b_mu1_RICHThresholdEl;   //!
   TBranch        *b_mu1_RICHThresholdMu;   //!
   TBranch        *b_mu1_RICHThresholdPi;   //!
   TBranch        *b_mu1_RICHThresholdKa;   //!
   TBranch        *b_mu1_RICHThresholdPr;   //!
   TBranch        *b_mu1_RICHAerogelUsed;   //!
   TBranch        *b_mu1_RICH1GasUsed;   //!
   TBranch        *b_mu1_RICH2GasUsed;   //!
   TBranch        *b_mu1_TRACK_Eta;   //!
   TBranch        *b_mu1_TRACK_Phi;   //!
   TBranch        *b_mu1_Aerogel_X;   //!
   TBranch        *b_mu1_Aerogel_Y;   //!
   TBranch        *b_mu1_Aerogel_Z;   //!
   TBranch        *b_mu1_Aerogel_Rho;   //!
   TBranch        *b_mu1_Aerogel_Phi;   //!
   TBranch        *b_mu1_Rich1Gas_X;   //!
   TBranch        *b_mu1_Rich1Gas_Y;   //!
   TBranch        *b_mu1_Rich1Gas_Z;   //!
   TBranch        *b_mu1_Rich1Gas_Rho;   //!
   TBranch        *b_mu1_Rich1Gas_Phi;   //!
   TBranch        *b_mu1_Rich2Gas_X;   //!
   TBranch        *b_mu1_Rich2Gas_Y;   //!
   TBranch        *b_mu1_Rich2Gas_Z;   //!
   TBranch        *b_mu1_Rich2Gas_Rho;   //!
   TBranch        *b_mu1_Rich2Gas_Phi;   //!
   TBranch        *b_mu1_L0Global_Dec;   //!
   TBranch        *b_mu1_L0Global_TIS;   //!
   TBranch        *b_mu1_L0Global_TOS;   //!
   TBranch        *b_mu1_Hlt1Global_Dec;   //!
   TBranch        *b_mu1_Hlt1Global_TIS;   //!
   TBranch        *b_mu1_Hlt1Global_TOS;   //!
   TBranch        *b_mu1_Hlt1Phys_Dec;   //!
   TBranch        *b_mu1_Hlt1Phys_TIS;   //!
   TBranch        *b_mu1_Hlt1Phys_TOS;   //!
   TBranch        *b_mu1_Hlt2Global_Dec;   //!
   TBranch        *b_mu1_Hlt2Global_TIS;   //!
   TBranch        *b_mu1_Hlt2Global_TOS;   //!
   TBranch        *b_mu1_Hlt2Phys_Dec;   //!
   TBranch        *b_mu1_Hlt2Phys_TIS;   //!
   TBranch        *b_mu1_Hlt2Phys_TOS;   //!
   TBranch        *b_mu1_TRACK_Type;   //!
   TBranch        *b_mu1_TRACK_Key;   //!
   TBranch        *b_mu1_TRACK_CHI2;   //!
   TBranch        *b_mu1_TRACK_NDOF;   //!
   TBranch        *b_mu1_TRACK_CHI2NDOF;   //!
   TBranch        *b_mu1_TRACK_PCHI2;   //!
   TBranch        *b_mu1_TRACK_VeloCHI2NDOF;   //!
   TBranch        *b_mu1_TRACK_TCHI2NDOF;   //!
   TBranch        *b_mu1_VELO_UTID;   //!
   TBranch        *b_mu1_TRACK_FirstMeasurementX;   //!
   TBranch        *b_mu1_TRACK_FirstMeasurementY;   //!
   TBranch        *b_mu1_TRACK_FirstMeasurementZ;   //!
   TBranch        *b_mu1_TRACK_MatchCHI2;   //!
   TBranch        *b_mu1_TRACK_GhostProb;   //!
   TBranch        *b_mu1_TRACK_CloneDist;   //!
   TBranch        *b_mu1_TRACK_Likelihood;   //!
   TBranch        *b_mu1_L0HadronDecision_Dec;   //!
   TBranch        *b_mu1_L0HadronDecision_TIS;   //!
   TBranch        *b_mu1_L0HadronDecision_TOS;   //!
   TBranch        *b_mu1_L0MuonDecision_Dec;   //!
   TBranch        *b_mu1_L0MuonDecision_TIS;   //!
   TBranch        *b_mu1_L0MuonDecision_TOS;   //!
   TBranch        *b_mu1_L0DiMuonDecision_Dec;   //!
   TBranch        *b_mu1_L0DiMuonDecision_TIS;   //!
   TBranch        *b_mu1_L0DiMuonDecision_TOS;   //!
   TBranch        *b_mu1_L0ElectronDecision_Dec;   //!
   TBranch        *b_mu1_L0ElectronDecision_TIS;   //!
   TBranch        *b_mu1_L0ElectronDecision_TOS;   //!
   TBranch        *b_mu1_L0PhotonDecision_Dec;   //!
   TBranch        *b_mu1_L0PhotonDecision_TIS;   //!
   TBranch        *b_mu1_L0PhotonDecision_TOS;   //!
   TBranch        *b_mu1_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_mu1_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_mu1_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_mu1_Hlt1DiMuonLowMassDecision_Dec;   //!
   TBranch        *b_mu1_Hlt1DiMuonLowMassDecision_TIS;   //!
   TBranch        *b_mu1_Hlt1DiMuonLowMassDecision_TOS;   //!
   TBranch        *b_mu1_Hlt1SingleMuonNoIPDecision_Dec;   //!
   TBranch        *b_mu1_Hlt1SingleMuonNoIPDecision_TIS;   //!
   TBranch        *b_mu1_Hlt1SingleMuonNoIPDecision_TOS;   //!
   TBranch        *b_mu1_Hlt1SingleMuonHighPTDecision_Dec;   //!
   TBranch        *b_mu1_Hlt1SingleMuonHighPTDecision_TIS;   //!
   TBranch        *b_mu1_Hlt1SingleMuonHighPTDecision_TOS;   //!
   TBranch        *b_mu1_Hlt1TrackAllL0Decision_Dec;   //!
   TBranch        *b_mu1_Hlt1TrackAllL0Decision_TIS;   //!
   TBranch        *b_mu1_Hlt1TrackAllL0Decision_TOS;   //!
   TBranch        *b_mu1_Hlt1TrackMuonDecision_Dec;   //!
   TBranch        *b_mu1_Hlt1TrackMuonDecision_TIS;   //!
   TBranch        *b_mu1_Hlt1TrackMuonDecision_TOS;   //!
   TBranch        *b_mu1_Hlt1TrackPhotonDecision_Dec;   //!
   TBranch        *b_mu1_Hlt1TrackPhotonDecision_TIS;   //!
   TBranch        *b_mu1_Hlt1TrackPhotonDecision_TOS;   //!
   TBranch        *b_mu1_Hlt1L0AnyDecision_Dec;   //!
   TBranch        *b_mu1_Hlt1L0AnyDecision_TIS;   //!
   TBranch        *b_mu1_Hlt1L0AnyDecision_TOS;   //!
   TBranch        *b_mu1_Hlt1GlobalDecision_Dec;   //!
   TBranch        *b_mu1_Hlt1GlobalDecision_TIS;   //!
   TBranch        *b_mu1_Hlt1GlobalDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2SingleMuonDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2SingleMuonDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2SingleMuonDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD2HMuMuDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD2HMuMuDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD2HMuMuDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2KMuMuDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2PiMuMuSSDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilep3bodyD2KMuMuSSDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD02PiPiMuMuDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD02KKMuMuDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD02KKMuMuDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD02KKMuMuDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD02KPiMuMuDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD02KPiMuMuDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmSemilepD02KPiMuMuDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHHDst_4piDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHHDst_4piDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHHDst_4piDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHHDst_K3piDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHHDst_KKpipiDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_hhXDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_hhXDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_hhXDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_hhXWideMassDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_BaryonhhXWideMassDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHXDst_LeptonhhXWideMassDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_K3piDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_K3piDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_K3piDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_K3piWideMassDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_KKpipiDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_KKpipiWideMassDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_4piDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_4piDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_4piDecision_TOS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_Dec;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TIS;   //!
   TBranch        *b_mu1_Hlt2CharmHadD02HHHH_4piWideMassDecision_TOS;   //!
   TBranch        *b_Slowpi_MINIP;   //!
   TBranch        *b_Slowpi_MINIPCHI2;   //!
   TBranch        *b_Slowpi_MINIPNEXTBEST;   //!
   TBranch        *b_Slowpi_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_Slowpi_OWNPV_X;   //!
   TBranch        *b_Slowpi_OWNPV_Y;   //!
   TBranch        *b_Slowpi_OWNPV_Z;   //!
   TBranch        *b_Slowpi_OWNPV_XERR;   //!
   TBranch        *b_Slowpi_OWNPV_YERR;   //!
   TBranch        *b_Slowpi_OWNPV_ZERR;   //!
   TBranch        *b_Slowpi_OWNPV_CHI2;   //!
   TBranch        *b_Slowpi_OWNPV_NDOF;   //!
   TBranch        *b_Slowpi_OWNPV_COV_;   //!
   TBranch        *b_Slowpi_IP_OWNPV;   //!
   TBranch        *b_Slowpi_IPCHI2_OWNPV;   //!
   TBranch        *b_Slowpi_TOPPV_X;   //!
   TBranch        *b_Slowpi_TOPPV_Y;   //!
   TBranch        *b_Slowpi_TOPPV_Z;   //!
   TBranch        *b_Slowpi_TOPPV_XERR;   //!
   TBranch        *b_Slowpi_TOPPV_YERR;   //!
   TBranch        *b_Slowpi_TOPPV_ZERR;   //!
   TBranch        *b_Slowpi_TOPPV_CHI2;   //!
   TBranch        *b_Slowpi_TOPPV_NDOF;   //!
   TBranch        *b_Slowpi_TOPPV_COV_;   //!
   TBranch        *b_Slowpi_IP_TOPPV;   //!
   TBranch        *b_Slowpi_IPCHI2_TOPPV;   //!
   TBranch        *b_Slowpi_ORIVX_X;   //!
   TBranch        *b_Slowpi_ORIVX_Y;   //!
   TBranch        *b_Slowpi_ORIVX_Z;   //!
   TBranch        *b_Slowpi_ORIVX_XERR;   //!
   TBranch        *b_Slowpi_ORIVX_YERR;   //!
   TBranch        *b_Slowpi_ORIVX_ZERR;   //!
   TBranch        *b_Slowpi_ORIVX_CHI2;   //!
   TBranch        *b_Slowpi_ORIVX_NDOF;   //!
   TBranch        *b_Slowpi_ORIVX_COV_;   //!
   TBranch        *b_Slowpi_IP_ORIVX;   //!
   TBranch        *b_Slowpi_IPCHI2_ORIVX;   //!
   TBranch        *b_Slowpi_P;   //!
   TBranch        *b_Slowpi_PT;   //!
   TBranch        *b_Slowpi_PE;   //!
   TBranch        *b_Slowpi_PX;   //!
   TBranch        *b_Slowpi_PY;   //!
   TBranch        *b_Slowpi_PZ;   //!
   TBranch        *b_Slowpi_M;   //!
   TBranch        *b_Slowpi_L0Calo_HCAL_realET;   //!
   TBranch        *b_Slowpi_L0Calo_HCAL_xProjection;   //!
   TBranch        *b_Slowpi_L0Calo_HCAL_yProjection;   //!
   TBranch        *b_Slowpi_L0Calo_HCAL_region;   //!
   TBranch        *b_Slowpi_L0Calo_HCAL_TriggerET;   //!
   TBranch        *b_Slowpi_L0Calo_HCAL_TriggerHCALET;   //!
   TBranch        *b_Slowpi_L0Calo_HCAL_xTrigger;   //!
   TBranch        *b_Slowpi_L0Calo_HCAL_yTrigger;   //!
   TBranch        *b_Slowpi_ID;   //!
   TBranch        *b_Slowpi_CombDLLMu;   //!
   TBranch        *b_Slowpi_ProbNNmu;   //!
   TBranch        *b_Slowpi_ProbNNghost;   //!
   TBranch        *b_Slowpi_InMuonAcc;   //!
   TBranch        *b_Slowpi_MuonDist2;   //!
   TBranch        *b_Slowpi_regionInM2;   //!
   TBranch        *b_Slowpi_hasMuon;   //!
   TBranch        *b_Slowpi_isMuon;   //!
   TBranch        *b_Slowpi_isMuonLoose;   //!
   TBranch        *b_Slowpi_NShared;   //!
   TBranch        *b_Slowpi_MuonLLmu;   //!
   TBranch        *b_Slowpi_MuonLLbg;   //!
   TBranch        *b_Slowpi_isMuonFromProto;   //!
   TBranch        *b_Slowpi_PIDe;   //!
   TBranch        *b_Slowpi_PIDmu;   //!
   TBranch        *b_Slowpi_PIDK;   //!
   TBranch        *b_Slowpi_PIDp;   //!
   TBranch        *b_Slowpi_ProbNNe;   //!
   TBranch        *b_Slowpi_ProbNNk;   //!
   TBranch        *b_Slowpi_ProbNNp;   //!
   TBranch        *b_Slowpi_ProbNNpi;   //!
   TBranch        *b_Slowpi_hasRich;   //!
   TBranch        *b_Slowpi_hasCalo;   //!
   TBranch        *b_Slowpi_UsedRichAerogel;   //!
   TBranch        *b_Slowpi_UsedRich1Gas;   //!
   TBranch        *b_Slowpi_UsedRich2Gas;   //!
   TBranch        *b_Slowpi_RichAboveElThres;   //!
   TBranch        *b_Slowpi_RichAboveMuThres;   //!
   TBranch        *b_Slowpi_RichAbovePiThres;   //!
   TBranch        *b_Slowpi_RichAboveKaThres;   //!
   TBranch        *b_Slowpi_RichAbovePrThres;   //!
   TBranch        *b_Slowpi_RichDLLe;   //!
   TBranch        *b_Slowpi_RichDLLmu;   //!
   TBranch        *b_Slowpi_RichDLLpi;   //!
   TBranch        *b_Slowpi_RichDLLk;   //!
   TBranch        *b_Slowpi_RichDLLp;   //!
   TBranch        *b_Slowpi_RichDLLbt;   //!
   TBranch        *b_Slowpi_InAccMuon;   //!
   TBranch        *b_Slowpi_MuonMuLL;   //!
   TBranch        *b_Slowpi_MuonBkgLL;   //!
   TBranch        *b_Slowpi_MuonNShared;   //!
   TBranch        *b_Slowpi_InAccEcal;   //!
   TBranch        *b_Slowpi_CaloEcalE;   //!
   TBranch        *b_Slowpi_EcalPIDe;   //!
   TBranch        *b_Slowpi_EcalPIDmu;   //!
   TBranch        *b_Slowpi_InAccHcal;   //!
   TBranch        *b_Slowpi_CaloHcalE;   //!
   TBranch        *b_Slowpi_HcalPIDe;   //!
   TBranch        *b_Slowpi_HcalPIDmu;   //!
   TBranch        *b_Slowpi_InAccPrs;   //!
   TBranch        *b_Slowpi_PrsPIDe;   //!
   TBranch        *b_Slowpi_CaloPrsE;   //!
   TBranch        *b_Slowpi_InAccSpd;   //!
   TBranch        *b_Slowpi_CaloSpdE;   //!
   TBranch        *b_Slowpi_InAccBrem;   //!
   TBranch        *b_Slowpi_BremPIDe;   //!
   TBranch        *b_Slowpi_VeloCharge;   //!
   TBranch        *b_Slowpi_RICHDLLe;   //!
   TBranch        *b_Slowpi_RICHDLLmu;   //!
   TBranch        *b_Slowpi_RICHDLLpi;   //!
   TBranch        *b_Slowpi_RICHDLLK;   //!
   TBranch        *b_Slowpi_RICHDLLp;   //!
   TBranch        *b_Slowpi_RICHDLLbt;   //!
   TBranch        *b_Slowpi_RICHBestID;   //!
   TBranch        *b_Slowpi_RICHThreshold;   //!
   TBranch        *b_Slowpi_RICHThresholdEl;   //!
   TBranch        *b_Slowpi_RICHThresholdMu;   //!
   TBranch        *b_Slowpi_RICHThresholdPi;   //!
   TBranch        *b_Slowpi_RICHThresholdKa;   //!
   TBranch        *b_Slowpi_RICHThresholdPr;   //!
   TBranch        *b_Slowpi_RICHAerogelUsed;   //!
   TBranch        *b_Slowpi_RICH1GasUsed;   //!
   TBranch        *b_Slowpi_RICH2GasUsed;   //!
   TBranch        *b_Slowpi_TRACK_Eta;   //!
   TBranch        *b_Slowpi_TRACK_Phi;   //!
   TBranch        *b_Slowpi_Aerogel_X;   //!
   TBranch        *b_Slowpi_Aerogel_Y;   //!
   TBranch        *b_Slowpi_Aerogel_Z;   //!
   TBranch        *b_Slowpi_Aerogel_Rho;   //!
   TBranch        *b_Slowpi_Aerogel_Phi;   //!
   TBranch        *b_Slowpi_Rich1Gas_X;   //!
   TBranch        *b_Slowpi_Rich1Gas_Y;   //!
   TBranch        *b_Slowpi_Rich1Gas_Z;   //!
   TBranch        *b_Slowpi_Rich1Gas_Rho;   //!
   TBranch        *b_Slowpi_Rich1Gas_Phi;   //!
   TBranch        *b_Slowpi_Rich2Gas_X;   //!
   TBranch        *b_Slowpi_Rich2Gas_Y;   //!
   TBranch        *b_Slowpi_Rich2Gas_Z;   //!
   TBranch        *b_Slowpi_Rich2Gas_Rho;   //!
   TBranch        *b_Slowpi_Rich2Gas_Phi;   //!
   TBranch        *b_Slowpi_L0Global_Dec;   //!
   TBranch        *b_Slowpi_L0Global_TIS;   //!
   TBranch        *b_Slowpi_L0Global_TOS;   //!
   TBranch        *b_Slowpi_Hlt1Global_Dec;   //!
   TBranch        *b_Slowpi_Hlt1Global_TIS;   //!
   TBranch        *b_Slowpi_Hlt1Global_TOS;   //!
   TBranch        *b_Slowpi_Hlt1Phys_Dec;   //!
   TBranch        *b_Slowpi_Hlt1Phys_TIS;   //!
   TBranch        *b_Slowpi_Hlt1Phys_TOS;   //!
   TBranch        *b_Slowpi_Hlt2Global_Dec;   //!
   TBranch        *b_Slowpi_Hlt2Global_TIS;   //!
   TBranch        *b_Slowpi_Hlt2Global_TOS;   //!
   TBranch        *b_Slowpi_Hlt2Phys_Dec;   //!
   TBranch        *b_Slowpi_Hlt2Phys_TIS;   //!
   TBranch        *b_Slowpi_Hlt2Phys_TOS;   //!
   TBranch        *b_Slowpi_TRACK_Type;   //!
   TBranch        *b_Slowpi_TRACK_Key;   //!
   TBranch        *b_Slowpi_TRACK_CHI2;   //!
   TBranch        *b_Slowpi_TRACK_NDOF;   //!
   TBranch        *b_Slowpi_TRACK_CHI2NDOF;   //!
   TBranch        *b_Slowpi_TRACK_PCHI2;   //!
   TBranch        *b_Slowpi_TRACK_VeloCHI2NDOF;   //!
   TBranch        *b_Slowpi_TRACK_TCHI2NDOF;   //!
   TBranch        *b_Slowpi_VELO_UTID;   //!
   TBranch        *b_Slowpi_TRACK_FirstMeasurementX;   //!
   TBranch        *b_Slowpi_TRACK_FirstMeasurementY;   //!
   TBranch        *b_Slowpi_TRACK_FirstMeasurementZ;   //!
   TBranch        *b_Slowpi_TRACK_MatchCHI2;   //!
   TBranch        *b_Slowpi_TRACK_GhostProb;   //!
   TBranch        *b_Slowpi_TRACK_CloneDist;   //!
   TBranch        *b_Slowpi_TRACK_Likelihood;   //!
   TBranch        *b_nCandidate;   //!
   TBranch        *b_totCandidates;   //!
   TBranch        *b_EventInSequence;   //!
   TBranch        *b_nBackward;   //!
   TBranch        *b_nDownstream;   //!
   TBranch        *b_nITClusters;   //!
   TBranch        *b_nLong;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_nPVs;   //!
   TBranch        *b_nSpdDigits;   //!
   TBranch        *b_nTTClusters;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_nUpstream;   //!
   TBranch        *b_nVELO;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_BCID;   //!
   TBranch        *b_BCType;   //!
   TBranch        *b_OdinTCK;   //!
   TBranch        *b_L0DUTCK;   //!
   TBranch        *b_HLT1TCK;   //!
   TBranch        *b_HLT2TCK;   //!
   TBranch        *b_GpsTime;   //!
   TBranch        *b_Polarity;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVX;   //!
   TBranch        *b_PVY;   //!
   TBranch        *b_PVZ;   //!
   TBranch        *b_PVXERR;   //!
   TBranch        *b_PVYERR;   //!
   TBranch        *b_PVZERR;   //!
   TBranch        *b_PVCHI2;   //!
   TBranch        *b_PVNDOF;   //!
   TBranch        *b_PVNTRACKS;   //!
  
   //MC Branches
   TBranch        *b_Dst_BKGCAT;   //!
   TBranch        *b_Dst_TRUEID;   //!
   TBranch        *b_Dst_MC_MOTHER_ID;   //!
   TBranch        *b_Dst_MC_MOTHER_KEY;   //!
   TBranch        *b_Dst_MC_GD_MOTHER_ID;   //!
   TBranch        *b_Dst_MC_GD_MOTHER_KEY;   //!
   TBranch        *b_Dst_MC_GD_GD_MOTHER_ID;   //!
   TBranch        *b_Dst_MC_GD_GD_MOTHER_KEY;   //!
   TBranch        *b_Dst_TRUEP_E;   //!
   TBranch        *b_Dst_TRUEP_X;   //!
   TBranch        *b_Dst_TRUEP_Y;   //!
   TBranch        *b_Dst_TRUEP_Z;   //!
   TBranch        *b_Dst_TRUEPT;   //!
   TBranch        *b_Dst_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_Dst_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_Dst_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_Dst_TRUEENDVERTEX_X;   //!
   TBranch        *b_Dst_TRUEENDVERTEX_Y;   //!
   TBranch        *b_Dst_TRUEENDVERTEX_Z;   //!
   TBranch        *b_Dst_TRUEISSTABLE;   //!
   TBranch        *b_Dst_TRUETAU;   //!
   TBranch        *b_D_BKGCAT;   //!
   TBranch        *b_D_TRUEID;   //!
   TBranch        *b_D_MC_MOTHER_ID;   //!
   TBranch        *b_D_MC_MOTHER_KEY;   //!
   TBranch        *b_D_MC_GD_MOTHER_ID;   //!
   TBranch        *b_D_MC_GD_MOTHER_KEY;   //!
   TBranch        *b_D_MC_GD_GD_MOTHER_ID;   //!
   TBranch        *b_D_MC_GD_GD_MOTHER_KEY;   //!
   TBranch        *b_D_TRUEP_E;   //!
   TBranch        *b_D_TRUEP_X;   //!
   TBranch        *b_D_TRUEP_Y;   //!
   TBranch        *b_D_TRUEP_Z;   //!
   TBranch        *b_D_TRUEPT;   //!
   TBranch        *b_D_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_D_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_D_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_D_TRUEENDVERTEX_X;   //!
   TBranch        *b_D_TRUEENDVERTEX_Y;   //!
   TBranch        *b_D_TRUEENDVERTEX_Z;   //!
   TBranch        *b_D_TRUEISSTABLE;   //!
   TBranch        *b_D_TRUETAU;   //!
   TBranch        *b_h0_TRUEID;   //!
   TBranch        *b_h0_MC_MOTHER_ID;   //!
   TBranch        *b_h0_MC_MOTHER_KEY;   //!
   TBranch        *b_h0_MC_GD_MOTHER_ID;   //!
   TBranch        *b_h0_MC_GD_MOTHER_KEY;   //!
   TBranch        *b_h0_MC_GD_GD_MOTHER_ID;   //!
   TBranch        *b_h0_MC_GD_GD_MOTHER_KEY;   //!
   TBranch        *b_h0_TRUEP_E;   //!
   TBranch        *b_h0_TRUEP_X;   //!
   TBranch        *b_h0_TRUEP_Y;   //!
   TBranch        *b_h0_TRUEP_Z;   //!
   TBranch        *b_h0_TRUEPT;   //!
   TBranch        *b_h0_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_h0_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_h0_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_h0_TRUEENDVERTEX_X;   //!
   TBranch        *b_h0_TRUEENDVERTEX_Y;   //!
   TBranch        *b_h0_TRUEENDVERTEX_Z;   //!
   TBranch        *b_h0_TRUEISSTABLE;   //!
   TBranch        *b_h0_TRUETAU;   //!
   TBranch        *b_h1_TRUEID;   //!
   TBranch        *b_h1_MC_MOTHER_ID;   //!
   TBranch        *b_h1_MC_MOTHER_KEY;   //!
   TBranch        *b_h1_MC_GD_MOTHER_ID;   //!
   TBranch        *b_h1_MC_GD_MOTHER_KEY;   //!
   TBranch        *b_h1_MC_GD_GD_MOTHER_ID;   //!
   TBranch        *b_h1_MC_GD_GD_MOTHER_KEY;   //!
   TBranch        *b_h1_TRUEP_E;   //!
   TBranch        *b_h1_TRUEP_X;   //!
   TBranch        *b_h1_TRUEP_Y;   //!
   TBranch        *b_h1_TRUEP_Z;   //!
   TBranch        *b_h1_TRUEPT;   //!
   TBranch        *b_h1_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_h1_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_h1_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_h1_TRUEENDVERTEX_X;   //!
   TBranch        *b_h1_TRUEENDVERTEX_Y;   //!
   TBranch        *b_h1_TRUEENDVERTEX_Z;   //!
   TBranch        *b_h1_TRUEISSTABLE;   //!
   TBranch        *b_h1_TRUETAU;   //!
   TBranch        *b_mu0_TRUEID;   //!
   TBranch        *b_mu0_MC_MOTHER_ID;   //!
   TBranch        *b_mu0_MC_MOTHER_KEY;   //!
   TBranch        *b_mu0_MC_GD_MOTHER_ID;   //!
   TBranch        *b_mu0_MC_GD_MOTHER_KEY;   //!
   TBranch        *b_mu0_MC_GD_GD_MOTHER_ID;   //!
   TBranch        *b_mu0_MC_GD_GD_MOTHER_KEY;   //!
   TBranch        *b_mu0_TRUEP_E;   //!
   TBranch        *b_mu0_TRUEP_X;   //!
   TBranch        *b_mu0_TRUEP_Y;   //!
   TBranch        *b_mu0_TRUEP_Z;   //!
   TBranch        *b_mu0_TRUEPT;   //!
   TBranch        *b_mu0_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_mu0_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_mu0_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_mu0_TRUEENDVERTEX_X;   //!
   TBranch        *b_mu0_TRUEENDVERTEX_Y;   //!
   TBranch        *b_mu0_TRUEENDVERTEX_Z;   //!
   TBranch        *b_mu0_TRUEISSTABLE;   //!
   TBranch        *b_mu0_TRUETAU;   //!
   TBranch        *b_mu1_TRUEID;   //!
   TBranch        *b_mu1_MC_MOTHER_ID;   //!
   TBranch        *b_mu1_MC_MOTHER_KEY;   //!
   TBranch        *b_mu1_MC_GD_MOTHER_ID;   //!
   TBranch        *b_mu1_MC_GD_MOTHER_KEY;   //!
   TBranch        *b_mu1_MC_GD_GD_MOTHER_ID;   //!
   TBranch        *b_mu1_MC_GD_GD_MOTHER_KEY;   //!
   TBranch        *b_mu1_TRUEP_E;   //!
   TBranch        *b_mu1_TRUEP_X;   //!
   TBranch        *b_mu1_TRUEP_Y;   //!
   TBranch        *b_mu1_TRUEP_Z;   //!
   TBranch        *b_mu1_TRUEPT;   //!
   TBranch        *b_mu1_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_mu1_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_mu1_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_mu1_TRUEENDVERTEX_X;   //!
   TBranch        *b_mu1_TRUEENDVERTEX_Y;   //!
   TBranch        *b_mu1_TRUEENDVERTEX_Z;   //!
   TBranch        *b_mu1_TRUEISSTABLE;   //!
   TBranch        *b_mu1_TRUETAU;   //!
   TBranch        *b_Slowpi_TRUEID;   //!
   TBranch        *b_Slowpi_MC_MOTHER_ID;   //!
   TBranch        *b_Slowpi_MC_MOTHER_KEY;   //!
   TBranch        *b_Slowpi_MC_GD_MOTHER_ID;   //!
   TBranch        *b_Slowpi_MC_GD_MOTHER_KEY;   //!
   TBranch        *b_Slowpi_MC_GD_GD_MOTHER_ID;   //!
   TBranch        *b_Slowpi_MC_GD_GD_MOTHER_KEY;   //!
   TBranch        *b_Slowpi_TRUEP_E;   //!
   TBranch        *b_Slowpi_TRUEP_X;   //!
   TBranch        *b_Slowpi_TRUEP_Y;   //!
   TBranch        *b_Slowpi_TRUEP_Z;   //!
   TBranch        *b_Slowpi_TRUEPT;   //!
   TBranch        *b_Slowpi_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_Slowpi_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_Slowpi_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_Slowpi_TRUEENDVERTEX_X;   //!
   TBranch        *b_Slowpi_TRUEENDVERTEX_Y;   //!
   TBranch        *b_Slowpi_TRUEENDVERTEX_Z;   //!
   TBranch        *b_Slowpi_TRUEISSTABLE;   //!
   TBranch        *b_Slowpi_TRUETAU;
   

};

#endif


