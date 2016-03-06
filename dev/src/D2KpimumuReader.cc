#include "D2KpimumuReader.h"

D2KpimumuReader::D2KpimumuReader(TTree *tree) 
{

  Init(tree);
}

D2KpimumuReader::~D2KpimumuReader()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}



void D2KpimumuReader::initializeMomenta(){                                                                                                                                                     
 
  pH1.SetXYZM(h1_PX,h1_PY,h1_PZ,Mass::K());
  pH0.SetXYZM(h0_PX,h0_PY,h0_PZ,Mass::Pi());
  pMu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
  pMu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());

  pD.SetXYZM(D_PX,D_PY,D_PZ,Mass::D0());
  pDst.SetXYZM(Dst_PX,Dst_PY,Dst_PZ,Mass::Ds());
  pPis.SetXYZM(Slowpi_PX,Slowpi_PY,Slowpi_PZ,Mass::Pi());

  pDTFDst.SetXYZM(Dst_DTF_Dstarplus_PX,Dst_DTF_Dstarplus_PY,Dst_DTF_Dstarplus_PZ,Mass::Ds());
  pDTFD.SetXYZM(Dst_DTF_D0_PX,Dst_DTF_D0_PY,Dst_DTF_D0_PZ,Mass::D0());
  pDTFPis.SetXYZM(Dst_DTF_Pis_PX,Dst_DTF_Pis_PY,Dst_DTF_Pis_PZ,Mass::Pi());

  pDTFH1.SetXYZM(Dst_DTF_h1_PX,Dst_DTF_h1_PY,Dst_DTF_h1_PZ,Mass::K());
  pDTFH0.SetXYZM(Dst_DTF_h0_PX,Dst_DTF_h0_PY,Dst_DTF_h0_PZ,Mass::Pi());
  pDTFMu1.SetXYZM(Dst_DTF_mu1_PX,Dst_DTF_mu1_PY,Dst_DTF_mu1_PZ,Mass::Mu());
  pDTFMu0.SetXYZM(Dst_DTF_mu0_PX,Dst_DTF_mu0_PY,Dst_DTF_mu0_PZ,Mass::Mu());

}


bool D2KpimumuReader::isHlt2Selected(){
  
  if(D_Hlt2CharmSemilepD02KPiMuMuDecision_TOS==1) return true;
  else return false;
}

