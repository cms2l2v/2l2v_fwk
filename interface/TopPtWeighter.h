#ifndef _toptweighter_hh_
#define _toptweighter_hh_

#include "TH1.h"
#include "TGraph.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>
#include <vector>

class TopPtWeighter
{
 public:

  /**
     @short CTOR
  */
  TopPtWeighter(TString procTag, TString outDir, TString shapesDir, TTree *evTree) :
    wgtGr_(0), wgtGrUp_(0), wgtGrDown_(0), curWgt_(1.0), curWgtUp_(1.0), curWgtDown_(1.0)
    { 
      //first look for the weights in the outDir
      TString wgtsFileUrl = outDir + "/" + procTag + "_toppt.root";
      gSystem->ExpandPathName( wgtsFileUrl );
      TFile *wgtsFile = TFile::Open(wgtsFileUrl);
      if(wgtsFile && !wgtsFile->IsZombie())
	{
	  std::cout << "[TopPtWeighter] will use weights previously stored @ " << wgtsFileUrl << std::endl;
	  wgtGr_    =(TGraph *)wgtsFile->Get("topptwgt");
	  wgtGrUp_  =(TGraph *)wgtsFile->Get("topptwgt_up");
	  wgtGrDown_=(TGraph *)wgtsFile->Get("topptwgt_down");
	  wgtsFile->Close();
	}
      if( wgtGr_!=0 ) return;
      
      
      //no weights file was found compute again...
      TObjArray *tkns=procTag.Tokenize("_");
      float mtop(172.5);
      for(int itkn=0; itkn<tkns->GetEntriesFast(); itkn++)
	{
	  TString mtopstr(tkns->At(itkn)->GetName());
	  mtopstr.ReplaceAll("v",".");
	  float imtop=mtopstr.Atof();
	  if(imtop<100) continue;
	  mtop=imtop;
	}
      std::cout << "[TopPtWeighter] will compute the top pT weights for m=" << mtop << " and store them in a local file" << std::endl;
      
      TString topPtShapeUrl(shapesDir); topPtShapeUrl += "/toppt_approxnnlo_8TeV.root";
      gSystem->ExpandPathName(topPtShapeUrl);
      TFile *topPtF=TFile::Open(topPtShapeUrl);
      char buf[100];
      sprintf(buf,"m%d",int(mtop*10));
      TH1 *targetH=(TH1 *)topPtF->Get(buf);
      if(targetH==0){
	std::cout << "[TopPtWeighter] error: could not find target shape " << buf << std::endl; 
	return;
      }
      targetH->SetDirectory(0);
      topPtF->Close();
      
      //init the reference histograms
      TH1 *genH=(TH1 *)targetH->Clone("gen");
      genH->SetDirectory(0);
      genH->Reset("ICE");
      
      //determine the weights that normalize the shape by running over the tree
      if(evTree==0){
	std::cout << "[TopPtWeighter] error: need an events tree to compute generated shape " << std::endl;
	return;
      }
      const Int_t totalEntries=evTree->GetEntriesFast();      
      Int_t mcn, mc_id[1000];
      Float_t  mc_px[1000],mc_py[1000];
      evTree->SetBranchAddress("mcn",   &mcn);
      evTree->SetBranchAddress("mc_id", mc_id);
      evTree->SetBranchAddress("mc_px", mc_px);
      evTree->SetBranchAddress("mc_py", mc_py);
      for(int inum=0; inum<totalEntries; ++inum)
	{
	  evTree->GetEntry(inum);
	  for(Int_t igen=0; igen<mcn; ++igen){
	    if(abs(mc_id[igen])!=6) continue;
	    float pt=sqrt(pow(mc_px[igen],2)+pow(mc_py[igen],2));
	    genH->Fill(pt);
	  }
	}
      
      //scale to the same number of events
      targetH->Scale(genH->Integral()/targetH->Integral());
      //std::cout << "[TopPtWeighter] average pT[GeV] in target is " << targetH->GetMean() << " generator is: " << genH->GetMean() << std::endl;

      //compute weights
      TH1 *ratio=(TH1 *) targetH->Clone("ratio"); ratio->Divide(genH); wgtGr_=new TGraph(ratio); delete ratio;
      wgtGrUp_ = new TGraph; 
      wgtGrDown_ = new TGraph;
      for(Int_t ip=0; ip<wgtGr_->GetN(); ip++)
	{
	  Double_t x,y;
	  wgtGr_->GetPoint(ip,x,y);
	  wgtGrUp_->SetPoint(ip,x,y);
	  wgtGrDown_->SetPoint(ip,x,1.0);
	}
      
      //determine the final normalization weights
      float totalExpected(0), wgtSum(0), wgtSumUp(0), wgtSumDown(0);
      for(int inum=0; inum<totalEntries; ++inum)
	{
	  evTree->GetEntry(inum);
	  
	  float pttop(-1), ptantitop(-1);
	  for(Int_t igen=0; igen<mcn; ++igen){
	    if(abs(mc_id[igen])!=6) continue;
	    float pt=sqrt(pow(mc_px[igen],2)+pow(mc_py[igen],2));
	    if(mc_id[igen]==6) pttop=pt;
	    else               ptantitop=pt;
	    if(pttop>0 && ptantitop>0) break;
	  }
	  if(pttop<0 || ptantitop<0) continue;
	  computeWeight(pttop,ptantitop);
	  totalExpected += 1.0;
	  wgtSum        += curWgt_;
	  wgtSumUp      += curWgtUp_;
	  wgtSumDown    += curWgtDown_;
	}
      //std::cout << "[TopPtWeighter] normalization computed for " << totalExpected << " events" << std::endl;
      
      //finalize the normalization
      if(totalExpected!=0){
	for(Int_t ip=0; ip<wgtGr_->GetN(); ip++)
	  {
	    Double_t x,y;
	    wgtGr_->GetPoint(ip,x,y);
	    wgtGr_->SetPoint(ip,x,y*totalExpected/wgtSum);
	    
	    wgtGrUp_->GetPoint(ip,x,y);
	    wgtGrUp_->SetPoint(ip,x,y*totalExpected/wgtSumUp);
	    
	    wgtGrDown_->GetPoint(ip,x,y);
	    wgtGrDown_->SetPoint(ip,x,y*totalExpected/wgtSumDown);
	  }
      }
      else{
	std::cout << "[TopPtWeighter] error: no ttbar events have been found in the events tree" << std::endl;
      }
      
      wgtsFile = TFile::Open(wgtsFileUrl,"RECREATE");
      wgtGr_->Write("topptwgt");
      wgtGrUp_->Write("topptwgt_up");
      wgtGrDown_->Write("topptwgt_down");
      wgtsFile->Close();
      //std::cout << "[TopPtWeighter] normalized weights saved to " << wgtsFileUrl << std::endl;
    }
    
    // FWLite
    TopPtWeighter(TString procTag, TString outDir, TString shapesDir, fwlite::ChainEvent& ev) :
      wgtGr_(0), wgtGrUp_(0), wgtGrDown_(0), curWgt_(1.0), curWgtUp_(1.0), curWgtDown_(1.0)
      { 
	//first look for the weights in the outDir
	TString wgtsFileUrl = outDir + "/" + procTag + "_toppt.root";
	gSystem->ExpandPathName( wgtsFileUrl );
	TFile *wgtsFile = TFile::Open(wgtsFileUrl);
	if(wgtsFile && !wgtsFile->IsZombie())
	  {
	    std::cout << "[TopPtWeighter] will use weights previously stored @ " << wgtsFileUrl << std::endl;
	    wgtGr_    =(TGraph *)wgtsFile->Get("topptwgt");
	    wgtGrUp_  =(TGraph *)wgtsFile->Get("topptwgt_up");
	    wgtGrDown_=(TGraph *)wgtsFile->Get("topptwgt_down");
	    wgtsFile->Close();
	  }
	if( wgtGr_!=0 ) return;
	
	
	//no weights file was found compute again...
	TObjArray *tkns=procTag.Tokenize("_");
	float mtop(172.5);
	for(int itkn=0; itkn<tkns->GetEntriesFast(); itkn++)
	  {
	    TString mtopstr(tkns->At(itkn)->GetName());
	    mtopstr.ReplaceAll("v",".");
	    float imtop=mtopstr.Atof();
	    if(imtop<100) continue;
	    mtop=imtop;
	  }
	std::cout << "[TopPtWeighter] will compute the top pT weights for m=" << mtop << " and store them in a local file" << std::endl;
	
	TString topPtShapeUrl(shapesDir); topPtShapeUrl += "/toppt_approxnnlo_8TeV.root";
	gSystem->ExpandPathName(topPtShapeUrl);
	TFile *topPtF=TFile::Open(topPtShapeUrl);
	char buf[100];
	sprintf(buf,"m%d",int(mtop*10));
	TH1 *targetH=(TH1 *)topPtF->Get(buf);
	if(targetH==0){
	  std::cout << "[TopPtWeighter] error: could not find target shape " << buf << std::endl; 
	  return;
	}
	targetH->SetDirectory(0);
	topPtF->Close();
	
	//init the reference histograms
	TH1 *genH=(TH1 *)targetH->Clone("gen");
	genH->SetDirectory(0);
	genH->Reset("ICE");
	
	//determine the weights that normalize the shape by running over the tree
	if(ev.size()==0){
	  std::cout << "[TopPtWeighter] error: need an events tree to compute generated shape " << std::endl;
	  return;
	}
	const Int_t totalEntries=ev.size();      
	Int_t mcn, mc_id[1000];
	Float_t  mc_px[1000],mc_py[1000];
	//evTree->SetBranchAddress("mcn",   &mcn);
	//evTree->SetBranchAddress("mc_id", mc_id);
	//evTree->SetBranchAddress("mc_px", mc_px);
	//evTree->SetBranchAddress("mc_py", mc_py);
	for(int inum=0; inum<totalEntries; ++inum)
	  {
	    ev.to(inum);
	    fwlite::Handle< llvvGenParticleCollection > genPartCollHandle;
	    genPartCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
	    if(!genPartCollHandle.isValid()){printf("llvvGenParticleCollection Object NotFound\n");continue;}
	    llvvGenParticleCollection gen = *genPartCollHandle;
	    for(size_t igen=0; igen<gen.size(); ++igen){
	      if(abs(gen[igen].id)!=6) continue;
	      float pt=sqrt(pow(gen[igen].px(),2)+pow(gen[igen].py(),2));
	      genH->Fill(pt);
	    }
	  }
	
	//scale to the same number of events
	targetH->Scale(genH->Integral()/targetH->Integral());
	//std::cout << "[TopPtWeighter] average pT[GeV] in target is " << targetH->GetMean() << " generator is: " << genH->GetMean() << std::endl;
	
	//compute weights
	TH1 *ratio=(TH1 *) targetH->Clone("ratio"); ratio->Divide(genH); wgtGr_=new TGraph(ratio); delete ratio;
	wgtGrUp_ = new TGraph; 
	wgtGrDown_ = new TGraph;
	for(Int_t ip=0; ip<wgtGr_->GetN(); ip++)
	  {
	    Double_t x,y;
	    wgtGr_->GetPoint(ip,x,y);
	    wgtGrUp_->SetPoint(ip,x,y);
	    wgtGrDown_->SetPoint(ip,x,1.0);
	  }
	
	//determine the final normalization weights
	float totalExpected(0), wgtSum(0), wgtSumUp(0), wgtSumDown(0);
	for(int inum=0; inum<totalEntries; ++inum)
	  {
	    ev.to(inum);
	    
	    float pttop(-1), ptantitop(-1);
	    fwlite::Handle< llvvGenParticleCollection > genPartCollHandle;
	    genPartCollHandle.getByLabel(ev, "llvvObjectProducersUsed");
	    if(!genPartCollHandle.isValid()){printf("llvvGenParticleCollection Object NotFound\n");continue;}
	    llvvGenParticleCollection gen = *genPartCollHandle;
	    for(size_t igen=0; igen<gen.size(); ++igen){
	      if(abs(gen[igen].id)!=6) continue;
	      float pt=sqrt(pow(gen[igen].px(),2)+pow(gen[igen].py(),2));
	      if(gen[igen].id==6) pttop=pt;
	      else               ptantitop=pt;
	      if(pttop>0 && ptantitop>0) break;
	    }
	    if(pttop<0 || ptantitop<0) continue;
	    computeWeight(pttop,ptantitop);
	    totalExpected += 1.0;
	    wgtSum        += curWgt_;
	    wgtSumUp      += curWgtUp_;
	    wgtSumDown    += curWgtDown_;
	  }
	//std::cout << "[TopPtWeighter] normalization computed for " << totalExpected << " events" << std::endl;
	
	//finalize the normalization
	if(totalExpected!=0){
	  for(Int_t ip=0; ip<wgtGr_->GetN(); ip++)
	    {
	      Double_t x,y;
	      wgtGr_->GetPoint(ip,x,y);
	      wgtGr_->SetPoint(ip,x,y*totalExpected/wgtSum);
	      
	      wgtGrUp_->GetPoint(ip,x,y);
	      wgtGrUp_->SetPoint(ip,x,y*totalExpected/wgtSumUp);
	      
	      wgtGrDown_->GetPoint(ip,x,y);
	      wgtGrDown_->SetPoint(ip,x,y*totalExpected/wgtSumDown);
	    }
	}
	else{
	  std::cout << "[TopPtWeighter] error: no ttbar events have been found in the events tree" << std::endl;
	}
	
	wgtsFile = TFile::Open(wgtsFileUrl,"RECREATE");
	wgtGr_->Write("topptwgt");
	wgtGrUp_->Write("topptwgt_up");
	wgtGrDown_->Write("topptwgt_down");
	wgtsFile->Close();
	//std::cout << "[TopPtWeighter] normalized weights saved to " << wgtsFileUrl << std::endl;
      }
      
      
    /**
       @short computes the ttbar pT weights for this event
    */
   inline void computeWeight(float pttop,float ptantitop)
     {
	if(wgtGr_==0) return;
	
	//saturate at 400 GeV (Kidonakis curves go to 500 but it's prone to stat variations also in MC)
	if(pttop>400)     pttop=400;
	if(ptantitop>400) ptantitop=400;

	float wgt(wgtGr_->Eval(pttop) * wgtGr_->Eval(ptantitop));  if(wgt<0) wgt=0;
	curWgt_     = sqrt( wgt );

	wgt=wgtGrUp_->Eval(pttop) * wgtGrUp_->Eval(ptantitop);     if(wgt<0) wgt=0;
	curWgtUp_   = wgt;

	wgt=wgtGrDown_->Eval(pttop) * wgtGrDown_->Eval(ptantitop); if(wgt<0) wgt=0;
	curWgtDown_ = wgt;
      }
    
    /**
       @short get the last events computed
     */
    inline void getEventWeight(float &nom, float &up, float &down)
    {
      nom=curWgt_;
      up=curWgtUp_;
      down=curWgtDown_;
    }
			 
    
    /**
       @short DTOR
     */
    ~TopPtWeighter()
      {
      }
  
 private:

  TGraph *wgtGr_, *wgtGrUp_, *wgtGrDown_;
  float curWgt_, curWgtUp_, curWgtDown_;
};

#endif
