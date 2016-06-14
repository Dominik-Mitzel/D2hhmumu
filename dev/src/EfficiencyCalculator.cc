#include "EfficiencyCalculator.h"

EfficiencyCalculator::EfficiencyCalculator(TString kind)
{

  tree_recoSignal=new TChain("BDT_Tree");
  tree_recoNorm =new TChain("BDT_Tree");
  tree_recoSigMisID=new TChain("BDT_Tree");
  tree_recoNormMisID=new TChain("BDT_Tree");

  tree_genSignal=new TChain("MCTruthTuple");
  tree_genNorm=new TChain("MCTruthTuple");
  tree_genSigMisID=new TChain("MCTruthTuple");
  tree_genNormMisID=new TChain("MCTruthTuple");
 
  m_kind=kind;
  //path to files is hardcoded here to location in HD
  pathToFiles="/auto/data/mitzel/D2hhmumu/new/preselectedSamples/";

  if(m_kind!="D2KKmumu" && m_kind!="D2pipimumu") {
    std::cout<<"ERROR:EfficienyCalculator : Channel specified is not implemented"<<std::endl; 
  }

  //add all the files(signal, normalization and misID) to the corresponding trees
  if(m_kind=="D2KKmumu") {
    tree_recoSignal->AddFile(pathToFiles+"MC_"+m_kind+"_BDT.root");
    tree_recoNorm->AddFile(pathToFiles+"MC_D2Kpimumu_"+m_kind+"BDT.root");
    tree_recoSigMisID->AddFile(pathToFiles+"MC_D2KKpipi_"+m_kind+"BDT.root");
    tree_recoNormMisID->AddFile(pathToFiles+"MC_D2Kpipipi_"+m_kind+"BDT.root");

    tree_genSignal->AddFile(pathToFiles+m_kind+"_MCgeneratorLevelTuple.root");
    tree_genNorm->AddFile(pathToFiles+"D2Kpimumu_MCgeneratorLevelTuple.root");
    tree_genSigMisID->AddFile(pathToFiles+"D2KKpipi_MCgeneratorLevelTuple.root");
    tree_genNormMisID->AddFile(pathToFiles+"D2Kpipipi_MCgeneratorLevelTuple.root");
  }
   
  if(m_kind=="D2pipimumu") {
    //PATHS HAVE TO BE SET PROPERLY LATER
    tree_recoSignal->AddFile(pathToFiles+"MC_"+m_kind+"_BDT.root");
    tree_recoNorm->AddFile(pathToFiles+"MC_D2Kpimumu_"+m_kind+"BDT.root");
    tree_recoSigMisID->AddFile(pathToFiles+"MC_D2KKpipi_"+m_kind+"BDT.root");
    tree_recoNormMisID->AddFile(pathToFiles+"MC_D2Kpipipi_"+m_kind+"BDT.root");

    tree_genSignal->AddFile(pathToFiles+m_kind+"_MCgeneratorLevelTuple.root");
    tree_genNorm->AddFile(pathToFiles+"D2Kpimumu_MCgeneratorLevelTuple.root");
    tree_genSigMisID->AddFile(pathToFiles+"D2KKpipi_MCgeneratorLevelTuple.root");
    tree_genNormMisID->AddFile(pathToFiles+"D2Kpipipi_MCgeneratorLevelTuple.root");
  }
   
}

double EfficiencyCalculator::getMCSignalEfficiency(TString selectionCut, TString splitting, TString q2Range){

  //factor is set to 2 of the data set is splitted to account for the fact, that only half of the data is used
  double factor = 1; 
  TString totalCut = selectionCut;
  TString normCut = q2Range;
  
  if(splitting!="") {
    totalCut+="&&"+splitting;
    //normCut+="&&"+splitting; //doesnt work right now as nTracks is not added to tuples
    factor =2;
  }
  if(q2Range!="") totalCut+="&&"+q2Range;
  
  double nSel = tree_recoSignal->GetEntries(totalCut);
  double nTot = tree_genSignal->GetEntries(normCut)/factor;
  std::cout<<"MCSignalEfficiency: "<<nSel/nTot<<std::endl;
  
  return nSel/nTot;

}

double EfficiencyCalculator::getMCSignalMisIDEfficiency(TString selectionCut, TString splitting, TString q2Range){
  
  double factor = 1;
  TString totalCut = selectionCut;
  TString normCut = q2Range;

  if(splitting!="") {
    totalCut+="&&"+splitting;
    //normCut+="&&"+splitting;
    factor =2;
  }
  if(q2Range!="") totalCut+="&&"+q2Range;
  
  double nSel = tree_recoSigMisID->GetEntries(totalCut);
  double nTot = tree_genSigMisID->GetEntries(normCut)/factor;
  std::cout<<"MCSignalMisIDEfficiency:"<<nSel/nTot<< "(nSel/nTot) " << nSel << " / " << nTot <<std::endl;
  
  return nSel/nTot;

}

double EfficiencyCalculator::getMCNormalizationEfficiency(TString selectionCut, TString splitting, TString q2Range){

  double factor = 1;
  if (q2Range =="") q2Range = "D_DiMuon_Mass>675&&D_DiMuon_Mass<875"; //fix if not otherwise stated
  TString totalCut = selectionCut+"&&"+q2Range;
  TString normCut = q2Range;
  
  if(splitting!="") {
    totalCut+="&&"+splitting;
    //normCut+="&&"+splitting;
    factor =2;
  }
  
  double nSel = tree_recoNorm->GetEntries(totalCut);
  double nTot = tree_genNorm->GetEntries(normCut)/factor;
  std::cout<<"MCNormalizationEfficiency: "<<nSel/nTot<< "(nSel/nTot) " << nSel << " / " << nTot << std::endl;
  
  return nSel/nTot;

}
double EfficiencyCalculator::getMCNormalizationMisIDEfficiency(TString selectionCut, TString splitting, TString q2Range){

  double factor = 1;
  if (q2Range =="") q2Range = "D_DiMuon_Mass>675&&D_DiMuon_Mass<875"; //fix if not otherwise stated
  TString totalCut = selectionCut+"&&"+q2Range;
  TString normCut= q2Range;
  
  if(splitting!="") {
    totalCut+="&&"+splitting;
    //normCut+="&&"+splitting;
    factor =2;
  }
  
  double nSel = tree_recoNormMisID->GetEntries(totalCut);
  double nTot = tree_genNormMisID->GetEntries(normCut)/factor;
  std::cout<<"MCNormalizationMisIDEfficiency: "<<nSel/nTot<< "(nSel/nTot) " << nSel << " / " << nTot << std::endl;
  
  return nSel/nTot;

}

double EfficiencyCalculator::getMCRelativeSigToNormEfficiency(TString selectionCut, TString splitting, TString q2RangeSig, TString q2RangeNorm){

  double relEff = getMCSignalEfficiency(selectionCut,splitting,q2RangeSig)/getMCNormalizationEfficiency(selectionCut,splitting,q2RangeNorm);
  std::cout<<"MCRelativeSigToNormEfficiency "<<relEff<<std::endl;
  return relEff;
}

double EfficiencyCalculator::getMCRelativeSigToNormMisIDEfficiency(TString selectionCut, TString splitting, TString q2RangeSig, TString q2RangeNorm){

  double relEff = getMCSignalMisIDEfficiency(selectionCut,splitting,q2RangeSig)/getMCNormalizationMisIDEfficiency(selectionCut,splitting,q2RangeNorm);
  std::cout<<"MCRelativeSigToNormMisIDEfficiency "<<relEff<<std::endl;
  return relEff;
}

double EfficiencyCalculator::getMisIDFractionQ2Range(TString q2Range){

  double nSel = tree_genSigMisID->GetEntries(q2Range);
  double nNorm = tree_genSigMisID->GetEntries();
  
  return nSel/nNorm;
}

EfficiencyCalculator::~EfficiencyCalculator(){};
