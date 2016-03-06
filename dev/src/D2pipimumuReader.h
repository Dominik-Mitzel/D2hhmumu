#ifndef D2pipimumuReader_h
#define D2pipimumuReader_h

#include "D2hhmumuReader.h"
#include "Tools.h"


class D2pipimumuReader : public D2hhmumuReader {

 public :

  D2pipimumuReader(TTree *tree);
  ~D2pipimumuReader();
  void initializeMomenta();
  bool isHlt2Selected();
};

#endif

