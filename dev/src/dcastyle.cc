#include "dcastyle.h"

// Use times new roman, precision 2                                                                                                           

/* main function that create the style. Running the macro with the optional argument value set to 1 will produce a test plot */
void dcastyle()
{
	TStyle *dcastyle = new TStyle("dcastyle","The real perfect style");

	Int_t lhcbFont     = 132;
	Double_t lhcbWidth = 2.00;
	Double_t lhcbTSize = 0.055;  

	// canvas
	dcastyle->SetCanvasColor(0);
	dcastyle->SetCanvasBorderSize(10);
	dcastyle->SetCanvasBorderMode(0);
	dcastyle->SetCanvasDefW(600);
	dcastyle->SetCanvasDefH(600);

    // paper
    dcastyle->SetPaperSize(20,26);
    
	// pads
	dcastyle->SetPadColor(0);
	dcastyle->SetPadBorderSize(10);
	dcastyle->SetPadBorderMode(0);
	dcastyle->SetPadLeftMargin(.18);
	dcastyle->SetPadRightMargin(.05);
	dcastyle->SetPadBottomMargin(.16);
	dcastyle->SetPadTopMargin(.07);
	dcastyle->SetPadGridX(0);
	dcastyle->SetPadGridY(0);
	dcastyle->SetPadTickX(1);
	dcastyle->SetPadTickY(1);

	// frame
	dcastyle->SetFrameBorderMode(0);
	dcastyle->SetFrameBorderSize(10);
	dcastyle->SetFrameFillStyle(0);
	dcastyle->SetFrameFillColor(0);
	dcastyle->SetFrameLineColor(0);
	dcastyle->SetFrameLineStyle(0);
	dcastyle->SetFrameLineWidth(0);

	// histogram
	dcastyle->SetHistFillColor(0);
	dcastyle->SetHistFillStyle(1001);// solid
	dcastyle->SetHistLineColor(1);
	dcastyle->SetHistLineStyle(0);
	dcastyle->SetHistLineWidth(lhcbWidth);
	dcastyle->SetHistMinimumZero();

    // look of the statistics box:
	dcastyle->SetOptStat(0);
	dcastyle->SetOptFit(0);
	dcastyle->SetStatColor(0);
    dcastyle->SetStatBorderSize(0);
    dcastyle->SetStatFont(lhcbFont);
    dcastyle->SetStatFontSize(0.05);
    dcastyle->SetStatX(0.9);
    dcastyle->SetStatY(0.9);
    dcastyle->SetStatW(0.25);
    dcastyle->SetStatH(0.15);
	
	// graph
	dcastyle->SetEndErrorSize(0);
	dcastyle->SetErrorX(0.5);

	// marker
	dcastyle->SetMarkerStyle(20);
	dcastyle->SetMarkerColor(kBlack);
	dcastyle->SetMarkerSize(0.5);

	// title 
	dcastyle->SetOptTitle(0);//1
	dcastyle->SetTitleFillColor(0);
	dcastyle->SetTitleBorderSize(0);
	dcastyle->SetTitleStyle(0);
    dcastyle->SetTitleX(0.0);
    dcastyle->SetTitleY(1.0);
    dcastyle->SetTitleW(1.0);
    dcastyle->SetTitleH(0.05);
    dcastyle->SetTitleFont(lhcbFont);

	// axes
	dcastyle->SetNdivisions(505,"X");
	dcastyle->SetNdivisions(510,"Y");
	
	dcastyle->SetTitleSize(1.1*lhcbTSize,"X");
    dcastyle->SetTitleFont(lhcbFont,"X");
	dcastyle->SetTitleOffset(.95,"X");
	dcastyle->SetLabelOffset(0.01,"X");
	dcastyle->SetLabelSize(lhcbTSize,"X");
	dcastyle->SetLabelFont(lhcbFont,"X");

	dcastyle->SetTitleSize(1.1*lhcbTSize,"Y");
	dcastyle->SetTitleFont(lhcbFont,"Y");
	dcastyle->SetTitleOffset(1.4,"Y");
	dcastyle->SetLabelOffset(0.01,"Y");
	dcastyle->SetLabelSize(lhcbTSize,"Y");
	dcastyle->SetLabelFont(lhcbFont,"Y");

	dcastyle->SetStripDecimals(kTRUE);

	dcastyle->SetTitleSize(1.1*lhcbTSize,"Z");
	dcastyle->SetTitleFont(lhcbFont,"Z");
	dcastyle->SetTitleOffset(1.800,"Z");
	dcastyle->SetLabelOffset(0.008,"Z");
	dcastyle->SetLabelSize(lhcbTSize,"Z");
	dcastyle->SetLabelFont(lhcbFont,"Z");
    
    dcastyle->SetTextSize(lhcbTSize);
	dcastyle->SetTextFont(lhcbFont);

	// function
	dcastyle->SetFuncColor(kBlue);
	dcastyle->SetFuncStyle(0);
	dcastyle->SetFuncWidth(lhcbWidth);

	// legend
    dcastyle->SetLegendBorderSize(0);
    dcastyle->SetLegendFillColor(0);
    dcastyle->SetLegendFont(42);

    // use medium bold lines
    dcastyle->SetLineWidth(lhcbWidth);
    dcastyle->SetGridWidth(lhcbWidth);
    dcastyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
    
    // If you want the usual gradient palette (blue -> red)
	dcastyle->SetPalette(1);
	dcastyle->SetNumberContours(20);
    // If you want colors that correspond to gray scale in black and white:
    //int colors[8] = {0,5,7,3,6,2,4,1};
    //dcastyle->SetPalette(8,colors);
    
	// set dcastyle as current style
	gROOT->SetStyle("dcastyle");
    gROOT->ForceStyle();

    //if(test) TestStyle();
}

/* function that tests the dcastyle */
void TestStyle(const char *tmpstyle)
{
	TString oldstyle = gStyle->GetName();
	gROOT->SetStyle(tmpstyle);

	TCanvas *c1 = new TCanvas("c1","c1");
    c1->cd();
	TH1F *h1 = new TH1F("h1","",100,-6.,6.);
	TH1F *h2 = new TH1F("h2","",100,-6.,6.);
    h1->SetXTitle("Variable");
    h1->SetYTitle("Events per bin");
	InitHist(h1,kBlack);
	InitHist(h2,kRed);

	TF1 *func = new TF1("func","gaus",-6.,6.);
    func->SetParameters(1.,1.,1.);
    h1->FillRandom("gaus",10000);
    h1->Draw("e");
    h2->FillRandom("func",5000);
    h2->Draw("esame");
	h1->Fit(func,"QR","same");

	TLatex *label = PlotTextLabel(Form("#chi^{2}/ndf = %.0f/%d",func->GetChisquare(),func->GetNDF()),0.92,0.88,"right");
    label->SetTextColor(kBlue);
    
    TCanvas *c2 = new TCanvas("c2","c2");
    c2->cd();
    TH2F *h3 = new TH2F("h3","",100,0.,5.,100,0.,5.);
    h3->SetXTitle("Variable 1");
    h3->SetYTitle("Variable 2");
	InitHist(h3);
    
    TF2 *func2 = new TF2("func2","abs(sin(x)*sin(y)/(x*y))",0.,5.,0.,5.);
    func2->SetLineColor(kBlack);
    for (int i=0; i<100000; i++) {
        double x,y;
        func2->GetRandom2(x,y);
        h3->Fill(x,y);
    }
    h3->Draw("col");
    func2->Draw("same");
    
	gROOT->SetStyle(oldstyle);

	return;
}

/* functions that initialize an histogram (or a graph) with the dcastyle */
void InitGraph(TGraph *histo, EColor color)
{
	histo->UseCurrentStyle();
	histo->SetLineColor(color);
	histo->SetMarkerColor(color);
}

void InitHist(TH1 *histo, EColor color)
{
	histo->UseCurrentStyle();
	histo->SetLineColor(color);
	histo->SetMarkerColor(color);
}

void InitHist(TH2 *histo, EColor color)
{
	histo->UseCurrentStyle();
	histo->SetLineColor(color);
	histo->SetMarkerStyle(1);
	histo->SetMarkerColor(color);
}

/* function that plots a text in (x,y) position */
TLatex *PlotTextLabel(const char *text, Double_t x, Double_t y, TString align, bool setNDC)
{
	Int_t lhcbFont     = 132;
	Double_t lhcbWidth = 2.00;
	Double_t lhcbTSize = 0.055;  

	gPad->cd();

	TLatex *pre_text = new TLatex(x,y,text);
	if (setNDC) pre_text->SetNDC();
    pre_text->SetTextFont(lhcbFont);
	pre_text->SetTextSize(lhcbTSize);
	if(align.Contains("right"))
		pre_text->SetTextAlign(33);
	else if(align.Contains("left"))
		pre_text->SetTextAlign(13);
	else
		pre_text->SetTextAlign(23);
	pre_text->Draw();

	return pre_text;
}

/* function to draw an histogram with visible color palette at the right of pad. */
void DrawWithColorPalette(TH2 *h, TPad *pad)
{
	if(!pad) pad = (TPad*) gPad;
	if(!pad) {
		new TCanvas();
		pad = (TPad*) gPad;
	}
	pad->cd();
	h->Draw("col");
	
	pad->SetRightMargin(.18);
	pad->Update();
	
	h->Draw("colz");
	return;
}

/* function to superimpose on a pad an histogram with different scale on the y-axis */
void DrawWithSecondScale(TH1* h, TPad *pad)
{
	if(!pad) pad = (TPad*) gPad;
	pad->cd();
	
	pad->SetRightMargin(pad->GetLeftMargin());
	pad->Update();
	
	Double_t rightmax = 1.10*h->GetMaximum();
	Double_t scale = (gPad->GetUymax()-gPad->GetUymin())/rightmax;
    
	h->Scale(scale);
	h->Draw("same");
    
	TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,505,"+L");
	axis->SetVertical();
	
	Color_t color = h->GetMarkerColor();
	
	axis->SetLineColor (color);
	axis->SetLabelColor(color);
	axis->SetTitleColor(color);
    
	TAxis *yaxis = h->GetYaxis();
	axis->SetTitle      (yaxis->GetTitle());
	axis->SetNdivisions (yaxis->GetNdivisions());
	axis->SetTitleSize  (yaxis->GetTitleSize());
	axis->SetTitleOffset(yaxis->GetTitleOffset());
	axis->SetTitleFont  (yaxis->GetTitleFont());
	axis->SetLabelOffset(yaxis->GetLabelOffset());
	axis->SetLabelSize  (yaxis->GetLabelSize());
	axis->SetLabelFont  (yaxis->GetLabelFont());
    
	axis->Draw();
    
	h->Draw("same");
    
	return;
}

/* function that creates a sub pad at the bottom of the canvas */
/*
void CreateSubPad(TCanvas *canvas, Double_t vfrac)
{
	canvas->SetCanvasSize(canvas->GetWw(),(1.+vfrac)*canvas->GetWh());
    
	Double_t xlow, ylow, xup, yup;
	canvas->GetPad(0)->GetPadPar(xlow,ylow,xup,yup);
    
	canvas->Divide(1,2);
	
	TVirtualPad *upPad = canvas->GetPad(1);
	upPad->SetPad(xlow,ylow+vfrac*(yup-ylow),xup,yup);
	
	TVirtualPad *dwPad = canvas->GetPad(2);
	dwPad->SetPad(xlow,ylow,xup,ylow+vfrac*(yup-ylow));
	
	canvas->Update();
	return;
}
*/
/* function that calculates the normalized residuals between an histograms and a function,
 and draws them in the current pad */
TH1F *DrawResiduals(TH1 *histo, TF1 *func)
{
	TAxis *xaxis = histo->GetXaxis();
	Int_t nbins  = histo->GetNbinsX();
	TH1F *histo_res = new TH1F(Form("res_%s",histo->GetName()),"",nbins,xaxis->GetXmin(),xaxis->GetXmax());
    
	histo_res->SetFillColor(kBlue);
	histo_res->GetYaxis()->SetTitle("#Delta/#sigma");
	histo_res->GetYaxis()->SetTitleSize(.15);
	histo_res->GetYaxis()->SetTitleOffset(.63);
	histo_res->GetYaxis()->SetLabelSize(.15);
	histo_res->GetXaxis()->SetLabelSize(.0);
    
	for(Int_t i=1; i<nbins; i++){
		Double_t delta = histo->GetBinContent(i)-func->Eval(histo->GetBinCenter(i));
		Double_t sigma = histo->GetBinError(i);
		if( sigma!=0 && histo->GetBinCenter(i)<func->GetXmax() && histo->GetBinCenter(i)>func->GetXmin() )
			histo_res->SetBinContent(i,delta/sigma);
	}
	histo_res->Draw();
	return histo_res;
}

/* function that calculates the normalized residuals between two histograms,
 and draws them in the current pad */
TH1F *DrawResiduals(TH1 *h1, TH1 *h2)
{
	if(!HistComp(h1,h2))
		return NULL;
	Int_t nbins1    = h1->GetNbinsX();
	TAxis *xaxis    = h1->GetXaxis();
	TH1F *histo_res = new TH1F(Form("res_%s_%s",h1->GetName(),h2->GetName()),"",nbins1,xaxis->GetXmin(),xaxis->GetXmax());
    
	histo_res->SetFillColor(kBlue);
	histo_res->GetYaxis()->SetTitle("#Delta/#sigma");
	histo_res->GetYaxis()->SetTitleSize(.15);
	histo_res->GetYaxis()->SetTitleOffset(.63);
	histo_res->GetYaxis()->SetLabelSize(.15);
	histo_res->GetXaxis()->SetLabelSize(.0);
    
	for(Int_t i=1; i<nbins1; i++){
		Double_t delta  = h1->GetBinContent(i)-h2->GetBinContent(i);
		Double_t sigma1 = h1->GetBinError(i);
		Double_t sigma2 = h2->GetBinError(i);
		Double_t sigma  = TMath::Sqrt(sigma1*sigma1+sigma2*sigma2);
		if( sigma!=0 )
			histo_res->SetBinContent(i,delta/sigma);
	}
	histo_res->Draw();
	return histo_res;
}

/* function that calculates the bin per bin ratio between two histograms,
 and draws it in the current pad */
TH1F *DrawRatio(TH1 *h1, TH1 *h2)
{
	if(!HistComp(h1,h2))
		return NULL;
	Int_t nbins1    = h1->GetNbinsX();
	TAxis *xaxis    = h1->GetXaxis();
	TH1F *histo_rat = new TH1F(Form("rat_%s_%s",h1->GetName(),h2->GetName()),"",nbins1,xaxis->GetXmin(),xaxis->GetXmax());
    
	histo_rat->GetYaxis()->SetTitle("Ratio");
	histo_rat->GetYaxis()->SetTitleSize(.15);
	histo_rat->GetYaxis()->SetTitleOffset(.63);
	histo_rat->GetYaxis()->SetLabelSize(.15);
	histo_rat->GetXaxis()->SetLabelSize(.0);
    
	for(Int_t i=1; i<=nbins1; i++){
		Double_t c1 = h1->GetBinContent(i);
		Double_t c2 = h2->GetBinContent(i);
		if(c1==0. || c2==0.) continue;
		Double_t ratio  = c1/c2;
		Double_t sigma1 = h1->GetBinError(i);
		Double_t sigma2 = h2->GetBinError(i);
		Double_t sigma  = ratio*TMath::Sqrt(sigma1*sigma1/(c1*c1)+sigma2*sigma2/(c2*c2));
		histo_rat->SetBinContent(i,ratio);
		histo_rat->SetBinError(i,sigma);
	}
	histo_rat->Draw();
	return histo_rat;
}

/* function that plot histogram statistical informations to the current pad */
TPaveStats *PlotStats(TH1 *h, TString option)
{
	TPaveStats *ptstats = new TPaveStats(0.659396,0.6573427,0.942953,0.8898601,"brNDC");
	ptstats->SetName("stats");
	ptstats->SetBorderSize(1);
	ptstats->SetFillColor(0);
	ptstats->SetTextAlign(12);
	TText *text = ptstats->AddText(Form("Entries = %g",h->GetEntries()));
	if(option.Contains("I"))
		text = ptstats->AddText(Form("Integral = %g",h->Integral()));
    
	Int_t N = h->GetDimension();
	if (N==1) {
		text = ptstats->AddText(Form("Mean  = %g",h->GetMean()));
		text = ptstats->AddText(Form("RMS   = %g",h->GetRMS()));
		if(option.Contains("S"))
			text = ptstats->AddText(Form("Skewness = %g",h->GetSkewness()));
		if(option.Contains("K"))
            text = ptstats->AddText(Form("Kurtosis = %g",h->GetKurtosis()));
		if(option.Contains("U")) {
			Int_t under = h->GetBinContent(0);
			text = ptstats->AddText(Form("Underflow = %d",under));
		}
		if(option.Contains("O")) {
			Int_t over = h->GetBinContent( h->GetNbinsX()+1 );
			text = ptstats->AddText(Form("Overflow  = %d",over));
		}
	} else {
        text = ptstats->AddText(Form("Mean x = %g",h->GetMean(1)));
        text = ptstats->AddText(Form("Mean y = %g",h->GetMean(2)));
        text = ptstats->AddText(Form("RMS x  = %g",h->GetRMS(1)));
        text = ptstats->AddText(Form("RMS y = %g",h->GetRMS(2)));
		if(option.Contains("S")) {
            text = ptstats->AddText(Form("Skewness x = %g",h->GetSkewness(1)));
            text = ptstats->AddText(Form("Skewness y = %g",h->GetSkewness(2)));
		}
		if(option.Contains("K")) {
            text = ptstats->AddText(Form("Kurtosis x = %g",h->GetKurtosis(1)));
            text = ptstats->AddText(Form("Kurtosis y = %g",h->GetKurtosis(2)));
		}
  		if(option.Contains("U") || option.Contains("O")) {
  			Int_t nbinsX = h->GetNbinsX();
  			Int_t nbinsY = h->GetNbinsY();
  			Double_t bl = h->GetBinContent(0,0);
  			Double_t bc = 0; for(int i=1; i<=nbinsX; i++) bc += h->GetBinContent(i,0);
  			Double_t br = h->GetBinContent(nbinsX+1,0);
  			Double_t cl = 0; for(int i=1; i<=nbinsY; i++) cl += h->GetBinContent(0,i);
			Double_t cc = h->GetEffectiveEntries();
  			Double_t cr = 0; for(int i=1; i<=nbinsX; i++) cr += h->GetBinContent(nbinsX+1,i);
  			Double_t tl = h->GetBinContent(0,nbinsY+1);
  			Double_t tc = 0; for(int i=1; i<=nbinsX; i++) tc += h->GetBinContent(i,nbinsY+1);
  			Double_t tr = h->GetBinContent(nbinsX+1,nbinsY+1);
			text = ptstats->AddText(Form("%g| %g| %g",tl,tc,tr));
			text = ptstats->AddText(Form("%g| %g| %g",cl,cc,cr));
			text = ptstats->AddText(Form("%g| %g| %g",bl,bc,br));
		}
	}
	ptstats->Draw();
	
	return ptstats;
}

/* function that plot fit results to the current pad */
TPaveStats *PlotFitResults(TF1 *f, TString option, TPaveStats *ptstats)
{
	if (!ptstats) {
		ptstats = new TPaveStats(0.659396,0.6573427,0.942953,0.8898601,"brNDC");
		ptstats->SetName("fitres");
		ptstats->SetBorderSize(0);
		ptstats->SetFillColor(0);
		ptstats->SetFillStyle(0);
		ptstats->SetTextAlign(12);
	}
	
	double chi2 = f->GetChisquare();
	TText *text = (chi2<2.) ? ptstats->AddText(Form("#chi^{2}/ndf = %.2f/%d",chi2,f->GetNDF()))
    : ptstats->AddText(Form("#chi^{2}/ndf = %.0f/%d",chi2,f->GetNDF()));
    
	if(option.Contains("P"))
		text = ptstats->AddText(Form("Prob = %g",f->GetProb()));
    
	Int_t npar = f->GetNpar();
	for (Int_t i=0; i<npar; ++i) {
		if (option.Contains("%")) {
			const char *format = option.Remove(0,option.First("%"));
			const char *line = Form("%s  = %s #pm %s","%s",format,format);
			text = ptstats->AddText(Form(line,f->GetParName(i),f->GetParameter(i),f->GetParError(i)));
		} else
			text = ptstats->AddText(Form("%s  = %g #pm %g",f->GetParName(i),f->GetParameter(i),f->GetParError(i)));
	}
	ptstats->Draw();
	
	return ptstats;
}

/* function that checks if two histograms are compatible */
void PrintError(const char *macro, const char *message)
{
	std::cout << "Error: " << macro << " - " << message << std::endl;
	return;
}

Bool_t HistComp(TH1 *h1, TH1 *h2)
{
    if (!h1 || !h2) {
        PrintError("HistComp","Nonexistent histogram!");
        return kFALSE;
    }
    
	Int_t dimension = h1->GetDimension();
	if( dimension != h2->GetDimension() ) {
		PrintError("HistComp","Histograms have different dimensions!");
		return kFALSE;
	}
    
	if ( h1->GetNbinsX() != h2->GetNbinsX() ) {
		PrintError("HistComp","Histograms have different number of bins!");
		return kFALSE;
	}
    
	if( (h1->GetXaxis()->GetXmin() != h2->GetXaxis()->GetXmin()) ||
       (h1->GetXaxis()->GetXmax() != h2->GetXaxis()->GetXmax()) ){
		PrintError("HistComp","Histograms have different limits!");
        return kFALSE;
    }
    
	if ( dimension > 1 ) {
		if ( h1->GetNbinsY() != h2->GetNbinsY() ) {
			PrintError("HistComp","Histograms have different number of bins!");
			return kFALSE;
		}
        
		if( (h1->GetYaxis()->GetXmin() != h2->GetYaxis()->GetXmin()) ||
           (h1->GetYaxis()->GetXmax() != h2->GetYaxis()->GetXmax()) ){
			PrintError("HistComp","Histograms have different limits!");
            return kFALSE;
        }
    }
    
	if ( dimension > 2 ) {
		if ( h1->GetNbinsZ() != h2->GetNbinsZ() ) {
			PrintError("HistComp","Histograms have different number of bins!");
			return kFALSE;
		}
        
		if( (h1->GetZaxis()->GetXmin() != h2->GetZaxis()->GetXmin()) ||
           (h1->GetZaxis()->GetXmax() != h2->GetZaxis()->GetXmax()) ){
			PrintError("HistComp","Histograms have different limits!");
            return kFALSE;
        }
    }
    
    return kTRUE;
}

