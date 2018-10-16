#ifndef histio_c
#define histio_c

#include "TClass.h"
#include "TList.h"
#include "TFile.h"
#include "TIterator.h"
#include "TRegexp.h"
#include "TDirectory.h"
#include "TObject.h"
#include "TSystem.h"
#include "TKey.h"
#include "TH1.h"

void saveHist(const char* filename, const char* pat, bool delete_hists = false ) ;

void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;

#endif

