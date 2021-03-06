//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 28 11:16:52 2015 by ROOT version 6.02/08
// from TTree LumiTuple/LumiTuple
// found on file: /auto/data/mitzel/D2hhmumu/2012/D2KKmumu/magUp/2012Data_D2KKmumu_13.root
//////////////////////////////////////////////////////////

#ifndef LumiTupleReader_h
#define LumiTupleReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class LumiTupleReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        IntegratedLuminosity;
   Double_t        IntegratedLuminosityErr;

   // List of branches
   TBranch        *b_IntegratedLuminosity;   //!
   TBranch        *b_IntegratedLuminosityErr;   //!

   LumiTupleReader(TTree *tree=0);
   virtual ~LumiTupleReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual  double   getTotalLumi();
};

#endif

#ifdef LumiTupleReader_cxx
LumiTupleReader::LumiTupleReader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/auto/data/mitzel/D2hhmumu/2012/D2KKmumu/magUp/2012Data_D2KKmumu_13.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/auto/data/mitzel/D2hhmumu/2012/D2KKmumu/magUp/2012Data_D2KKmumu_13.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/auto/data/mitzel/D2hhmumu/2012/D2KKmumu/magUp/2012Data_D2KKmumu_13.root:/GetIntegratedLuminosity");
      dir->GetObject("LumiTuple",tree);

   }
   Init(tree);
}

LumiTupleReader::~LumiTupleReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t LumiTupleReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t LumiTupleReader::LoadTree(Long64_t entry)
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

void LumiTupleReader::Init(TTree *tree)
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

   fChain->SetBranchAddress("IntegratedLuminosity", &IntegratedLuminosity, &b_IntegratedLuminosity);
   fChain->SetBranchAddress("IntegratedLuminosityErr", &IntegratedLuminosityErr, &b_IntegratedLuminosityErr);
   Notify();
}

Bool_t LumiTupleReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LumiTupleReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t LumiTupleReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef LumiTupleReader_cxx
