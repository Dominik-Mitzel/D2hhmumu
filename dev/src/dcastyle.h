#ifndef DCASTYLE_H
#define DCASTYLE_H
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH2.h>
#include <TGraph.h>
#include <TF2.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TGaxis.h>
#include <TMath.h>


// Style
void dcastyle();
void TestStyle(const char *tmpstyle = "dcastyle");
void InitGraph(TGraph *g, EColor color = kBlack);
void InitHist (TH1 *histo, EColor color = kBlack);
void InitHist (TH2 *histo, EColor color = kBlack);

// Plotting
TLatex     *PlotTextLabel       (const char *text, Double_t x = 0.93, Double_t y = 0.88, TString align="right", bool setNDC=true);
void        CreateSubPad        (TCanvas *c, Double_t vfrac = 0.25);
void        DrawWithColorPalette(TH2* h, TPad *pad = 0);
void        DrawWithSecondScale (TH1* h, TPad *pad = 0);
TH1F       *DrawResiduals       (TH1 *h,  TF1 *func);
TH1F       *DrawResiduals       (TH1 *h1, TH1 *h2);
TH1F       *DrawRatio           (TH1 *h1, TH1 *h2);
TPaveStats *PlotStats           (TH1 *h, TString option = "");
TPaveStats *PlotFitResults      (TF1 *f, TString option = "", TPaveStats *p = 0);

// Histos
Bool_t HistComp(TH1 *h1, TH1 *h2);
#endif // DCA_STYLE_H


