#include "D2hhmumuFitter.h"
#include "D2KKmumuReader.h"
#include "D2pipimumuReader.h"
#include "D2KpimumuReader.h"
#include "D2KKpipiReader.h"
#include "D2KpipipiReader.h"



void optimizeSelection();
double getMCSignalEfficiency(TString cut);
double EffD2KKpipiToEffD2Kpipipi(TString cut);

void D2pipimumuMC();
void D2pipimumuData();
void D2KpimumuMC();
void D2KpimumuData();
void D2KKmumuMC(); 
void D2KKmumuData(); 
void D2KKpipiMC();
void D2KpipipiMC();
void D2KKpipiData();
void D2KpipipiData();
