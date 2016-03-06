#include <TH1.h>
#include <TCanvas.h>
#include <TString.h>

void calculateSweights(TString pathToFile);
void data_MC_comparison();
void CreateSubPad        (TCanvas *c, Double_t vfrac );
std::vector<TH1*> fillHistograms(TString f_input,TString f_output, bool isData);
