#ifndef D2KpipipiReader_h
#define D2KpipipiReader_h

#include "D2hhmumuReader.h"
#include "Tools.h"
#include "TRandom3.h"

class D2KpipipiReader : public D2hhmumuReader {

 public :

  D2KpipipiReader(TTree *tree);
  ~D2KpipipiReader();
  void initializeMomenta();
  bool isHlt2Selected();

  //pion 4 momenta under correct mass hypothesis 
  TLorentzVector pDTFPi0;
  TLorentzVector pDTFPi1;
  TLorentzVector pPi0;
  TLorentzVector pPi1;
  TLorentzVector pDTFMu2;

  //get masses m(mumu) for decay in flight
  void   addMisIdMasses(TString name);
  void   createSubsample(TString name, double percentage);
  void   createMCtrainingSample(TString name);
  double get_mMuMu_noDCF();
  double get_mMuMu_doubleDCF();
  double get_mMuMu_DCF_lowP();
  double get_mMuMu_DCF_highP();
  void   createRandomizedSubsample(TString name);
};

#endif

