#include "D2KKmumuReader.h"
#include <iostream>

D2KKmumuReader::D2KKmumuReader(TTree *tree) 
{

  Init(tree);
}


D2KKmumuReader::~D2KKmumuReader()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}


void D2KKmumuReader::initializeMomenta(){                                                                                                                                                  


  pH1.SetXYZM(h1_PX,h1_PY,h1_PZ,Mass::K());
  pH0.SetXYZM(h0_PX,h0_PY,h0_PZ,Mass::K());
  pMu0.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Mu());
  pMu1.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Mu());

  pD.SetXYZM(D_PX,D_PY,D_PZ,Mass::D0());
  pDst.SetXYZM(Dst_PX,Dst_PY,Dst_PZ,Mass::Ds());
  pPis.SetXYZM(Slowpi_PX,Slowpi_PY,Slowpi_PZ,Mass::Pi());

  pDTFDst.SetXYZM(Dst_DTF_Dstarplus_PX,Dst_DTF_Dstarplus_PY,Dst_DTF_Dstarplus_PZ,Mass::Ds());
  pDTFD.SetXYZM(Dst_DTF_D0_PX,Dst_DTF_D0_PY,Dst_DTF_D0_PZ,Mass::D0());
  pDTFPis.SetXYZM(Dst_DTF_Pis_PX,Dst_DTF_Pis_PY,Dst_DTF_Pis_PZ,Mass::Pi());

  pDTFH1.SetXYZM(Dst_DTF_h1_PX,Dst_DTF_h1_PY,Dst_DTF_h1_PZ,Mass::K());
  pDTFH0.SetXYZM(Dst_DTF_h0_PX,Dst_DTF_h0_PY,Dst_DTF_h0_PZ,Mass::K());
  pDTFMu1.SetXYZM(Dst_DTF_mu1_PX,Dst_DTF_mu1_PY,Dst_DTF_mu1_PZ,Mass::Mu());
  pDTFMu0.SetXYZM(Dst_DTF_mu0_PX,Dst_DTF_mu0_PY,Dst_DTF_mu0_PZ,Mass::Mu());

  pMu0_hypPi.SetXYZM(mu0_PX,mu0_PY,mu0_PZ,Mass::Pi());
  pMu1_hypPi.SetXYZM(mu1_PX,mu1_PY,mu1_PZ,Mass::Pi());
  pH0_hypMu.SetXYZM(h0_PX,h0_PY,h0_PZ,Mass::Mu());
  pH1_hypMu.SetXYZM(h1_PX,h1_PY,h1_PZ,Mass::Mu());
  pH0_hypKaon.SetXYZM(h0_PX,h0_PY,h0_PZ,Mass::K());
  pH1_hypKaon.SetXYZM(h1_PX,h1_PY,h1_PZ,Mass::K());
  pH0_hypProton.SetXYZM(h0_PX,h0_PY,h0_PZ,Mass::P());
  pH1_hypProton.SetXYZM(h1_PX,h1_PY,h1_PZ,Mass::P());
  pH0_hypPi.SetXYZM(h0_PX,h0_PY,h0_PZ,Mass::Pi());
  pH1_hypPi.SetXYZM(h1_PX,h1_PY,h1_PZ,Mass::Pi());


}


bool D2KKmumuReader::isHlt2Selected(){

  if(D_Hlt2CharmSemilepD02KKMuMuDecision_TOS==1) return true;
  else return false;
}
