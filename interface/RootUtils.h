#ifndef _rootutils_h_
#define _rootutils_h_

#include "TH1.h"

namespace utils
{
  namespace root
  {
    //
    void fixExtremities(TH1* h,bool addOverflow, bool addUnderflow)
    {
      if(h==0) return;
      
      if(addUnderflow){
	double fbin  = h->GetBinContent(0) + h->GetBinContent(1);
	double fbine = sqrt(h->GetBinError(0)*h->GetBinError(0)
			    + h->GetBinError(1)*h->GetBinError(1));
	h->SetBinContent(1,fbin);
	h->SetBinError(1,fbine);
	h->SetBinContent(0,0);
	h->SetBinError(0,0);
      }
      
      if(addOverflow){  
	int nbins = h->GetNbinsX();
	double fbin  = h->GetBinContent(nbins) + h->GetBinContent(nbins+1);
	double fbine = sqrt(h->GetBinError(nbins)*h->GetBinError(nbins) 
			    + h->GetBinError(nbins+1)*h->GetBinError(nbins+1));
	h->SetBinContent(nbins,fbin);
	h->SetBinError(nbins,fbine);
	h->SetBinContent(nbins+1,0);
	h->SetBinError(nbins+1,0);
      }
    }
  }
}

#endif
