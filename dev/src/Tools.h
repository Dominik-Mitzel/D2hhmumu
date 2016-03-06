#ifndef TOOLS_H
#define TOOLS_H

// -----------------------------------------
// Particle masses [MeV]
// -----------------------------------------
namespace Mass
{
	inline Double_t E  () { return 0.000510998910e3; }
	inline Double_t Mu () { return 0.105658367e3; }
	inline Double_t Tau() { return 1.77684e3; }
	inline Double_t K  () { return 0.493677e3; }
	inline Double_t Pi () { return 0.13957018e3; }
	inline Double_t Pi0() { return 0.1349766e3; }
	inline Double_t D0 () { return 1.86484e3; }
	inline Double_t DS () { return 2.01027e3; }
	inline Double_t D  () { return 1.86962e3; }
	inline Double_t Ds () { return 1.96849e3; }
	inline Double_t P  () { return 0.938272013e3; }
	inline Double_t Phi() { return 1.019455e3; }
	inline Double_t B  () { return 5.27917e3; }
	inline Double_t B0 () { return 5.2795e3; }
	inline Double_t Bs () { return 5.3663e3; }
	inline Double_t K0 () { return 0.497614e3; }
	inline Double_t Lm () { return 1.115683e3; }
}

// -----------------------------------------
// Particle MC Codes
// -----------------------------------------
namespace PdgCode
{
	inline Int_t Gamma(){ return 22; }
	inline Int_t E   () { return 11; }
	inline Int_t NuE () { return 12; }
	inline Int_t Mu  () { return 13; }
	inline Int_t NuMu() { return 14; }
	inline Int_t Tau()  { return 15; }
	inline Int_t NuTau(){ return 16; }
	inline Int_t K   () { return 321; }
	inline Int_t Pi  () { return 211; }
	inline Int_t Pi0 () { return 111; }
	inline Int_t D0  () { return 421; }
	inline Int_t DS  () { return 413; }
	inline Int_t D   () { return 411; }
	inline Int_t Ds  () { return 431; }
	inline Int_t P   () { return 2212; }
	inline Int_t Phi () { return 333; }
	inline Int_t B   () { return 521; }
	inline Int_t B0  () { return 511; }
	inline Int_t Bs  () { return 531; }
	inline Int_t K0  () { return 311; }
	inline Int_t K0s () { return 310; }
	inline Int_t K0l () { return 130; }
	inline Int_t Lm  () { return 3122; }
}

#endif // TOOLS_H

