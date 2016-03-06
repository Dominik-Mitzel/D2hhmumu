#ifndef D2KpimumuReader_h
#define D2KpimumuReader_h

#include "D2hhmumuReader.h"
#include "Tools.h"
//#include "D2hhmumuReader.C"


class D2KpimumuReader : public D2hhmumuReader {

 public :

  D2KpimumuReader(TTree *tree);
  ~D2KpimumuReader();
  void initializeMomenta();
  bool isHlt2Selected();
};

#endif


