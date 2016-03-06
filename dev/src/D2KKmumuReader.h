#ifndef D2KKmumuReader_h
#define D2KKmumuReader_h

#include "D2hhmumuReader.h"
#include "Tools.h"


class D2KKmumuReader : public D2hhmumuReader {

 public :

  D2KKmumuReader(TTree *tree);
  ~D2KKmumuReader();
  void initializeMomenta();
  bool isHlt2Selected();
};

#endif

