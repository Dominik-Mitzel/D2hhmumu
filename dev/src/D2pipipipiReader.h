#ifndef D2pipipipiReader_h
#define D2pipipipiReader_h

#include "D2hhmumuReader.h"
#include "Tools.h"
#include "TRandom3.h"


class D2pipipipiReader : public D2hhmumuReader {

 public :

  D2pipipipiReader(TTree *tree);
  ~D2pipipipiReader();
  void initializeMomenta();
  bool isHlt2Selected();

  //pion 4 momenta under correct mass hypothesis 
  TLorentzVector pDTFPi0;
  TLorentzVector pDTFPi1;
  TLorentzVector pPi0;
  TLorentzVector pPi1;

  void   createSubsample(TString name, double percentage);
  void   createMCtrainingSample(TString name);
  //get masses m(mumu) for decay in flight
  void   addMisIdMasses(TString name);
  void   set_mMuMu_misID();
  void   createRandomizedSubsample(TString name);
  
  double mMuMu_noDCF;
  double mMuMu_doubleDCF;
  double mMuMu_DCF0;
  double mMuMu_DCF1;

  double decayedPion0_PT;
  double decayedPion1_PT;

};

#endif

