#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "TH1F.h"
#include "TSystem.h"

namespace utils
{
  namespace cmssw
  {
    //
    FactorizedJetCorrector* getJetCorrector(TString baseDir, bool isMC)
    {
      gSystem->ExpandPathName(baseDir);
      TString pf(isMC ? "MC" : "DATA");
      
      //order matters: L1 -> L2 -> L3 (-> Residuals)
      std::vector<std::string> jetCorFiles;
      std::cout << baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt" << std::endl;
      jetCorFiles.push_back((baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt").Data());
      jetCorFiles.push_back((baseDir+"/"+pf+"_L2Relative_AK4PFchs.txt").Data());
      jetCorFiles.push_back((baseDir+"/"+pf+"_L3Absolute_AK4PFchs.txt").Data());
      if(!isMC) jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
     
      //init the parameters for correction
      std::vector<JetCorrectorParameters> corSteps;
      for(size_t i=0; i<jetCorFiles.size(); i++) corSteps.push_back(JetCorrectorParameters(jetCorFiles[i]));
      
      //return the corrector
      return new FactorizedJetCorrector(corSteps);
    }
    
     std::vector<double> smearJER(double pt, double eta, double genPt){
         std::vector<double> toReturn(3,pt);
         if(genPt<=0) return toReturn;
         
         // FIXME: These are the 8 TeV values.
         //
         eta=fabs(eta);
         double ptSF(1.0), ptSF_err(0.06);
         if(eta<0.8)                  { ptSF=1.061; ptSF_err=sqrt(pow(0.012,2)+pow(0.023,2)); }
         else if(eta>=0.8 && eta<1.3) { ptSF=1.088; ptSF_err=sqrt(pow(0.012,2)+pow(0.029,2)); }
         else if(eta>=1.3 && eta<1.9) { ptSF=1.106; ptSF_err=sqrt(pow(0.017,2)+pow(0.030,2)); }
         else if(eta>=1.9 && eta<2.5) { ptSF=1.126; ptSF_err=sqrt(pow(0.035,2)+pow(0.094,2)); }
         else if(eta>=2.5 && eta<3.0) { ptSF=1.343; ptSF_err=sqrt(pow(0.127,2)+pow(0.123,2)); }
         else if(eta>=3.0 && eta<3.2) { ptSF=1.303; ptSF_err=sqrt(pow(0.127,2)+pow(1.303,2)); }
         else if(eta>=3.2 && eta<5.0) { ptSF=1.320; ptSF_err=sqrt(pow(0.127,2)+pow(1.320,2)); }
         
         toReturn[0]=TMath::Max(0.,((genPt+ptSF*(pt-genPt)))/pt);
         toReturn[1]=TMath::Max(0.,((genPt+(ptSF+ptSF_err)*(pt-genPt)))/pt);
         toReturn[2]=TMath::Max(0.,((genPt+(ptSF-ptSF_err)*(pt-genPt)))/pt);
         return toReturn;
     }
     
     //
     std::vector<float> smearJES(double pt, double eta, JetCorrectionUncertainty *jecUnc){
         jecUnc->setJetEta(eta);
         jecUnc->setJetPt(pt);
         double relShift=fabs(jecUnc->getUncertainty(true));
         std::vector<float> toRet;
         toRet.push_back((1.0+relShift));
         toRet.push_back((1.0-relShift));
         return toRet;
     }
     
     void updateJEC(pat::JetCollection& jets, FactorizedJetCorrector *jesCor, JetCorrectionUncertainty *totalJESUnc, float rho, int nvtx,bool isMC){
         for(size_t ijet=0; ijet<jets.size(); ijet++){
             pat::Jet& jet = jets[ijet];
             
             //correct JES
             LorentzVector rawJet = jet.correctedP4("Uncorrected");

             //double toRawSF=jet.correctedJet("Uncorrected").pt()/jet.pt();
             //LorentzVector rawJet(jet*toRawSF);
             jesCor->setJetEta(rawJet.eta());
             jesCor->setJetPt(rawJet.pt());
             jesCor->setJetA(jet.jetArea());
             jesCor->setRho(rho);
             jesCor->setNPV(nvtx);
             jet.setP4(rawJet*jesCor->getCorrection());

             //smear JER
             //double newJERSF(1.0);
             if(isMC){
                 const reco::GenJet* genJet=jet.genJet();
                 if(genJet){
                   double genjetpt( genJet ? genJet->pt(): 0.);                    
                    std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
                    jet.setP4(jet.p4()*smearJER[0]);
                
                    //printf("jet pt=%f gen pt = %f smearing %f %f %f\n", jet.pt(), genjetpt, smearJER[0], smearJER[1], smearJER[2]);
                    // //set the JER up/down alternatives
                    jet.addUserFloat("jerup", smearJER[1]);  //kept for backward compatibility
                    jet.addUserFloat("jerdown", smearJER[2] ); //kept for backward compatibility
                    jet.addUserFloat("_res_jup", smearJER[1]);
                    jet.addUserFloat("_res_jdown", smearJER[2] );
                 }else{
                    jet.addUserFloat("jerup", 1.0); //kept for backward compatibility
                    jet.addUserFloat("jerdown", 1.0);  //kept for backward compatibility
                    jet.addUserFloat("_res_jup", 1.0);
                    jet.addUserFloat("_res_jdown", 1.0 );
                 }
             }

             if(isMC){
                ////set the JES up/down pT alternatives
                std::vector<float> ptUnc=utils::cmssw::smearJES(jet.pt(),jet.eta(), totalJESUnc);
                jet.addUserFloat("jesup",    ptUnc[0] );  //kept for backward compatibility
                jet.addUserFloat("jesdown",  ptUnc[1] );  //kept for backward compatibility
                jet.addUserFloat("_scale_jup",    ptUnc[0] );
                jet.addUserFloat("_scale_jdown",  ptUnc[1] );
             }
             
             // FIXME: this is not to be re-set. Check that this is a desired non-feature.
             // i.e. check that the uncorrectedJet remains the same even when the corrected momentum is changed by this routine.
             //to get the raw jet again
             //jets[ijet].setVal("torawsf",1./(newJECSF*newJERSF));
         }
     }
     
//    //
//    std::vector<LorentzVector> getMETvariations(LorentzVector &rawMETP4, pat::JetCollection &jets, std::vector<patUtils::GenericLepton> &leptons,bool isMC)
//    {
//      std::vector<LorentzVector> newMetsP4(9,rawMETP4);
//      if(!isMC) return newMetsP4;
//      
//      LorentzVector nullP4(0,0,0,0);
//      
//      //recompute the clustered and unclustered fluxes with energy variations
//      for(size_t ivar=1; ivar<=8; ivar++)
//        {
//          
//          //leptonic flux
//          LorentzVector leptonFlux(nullP4), lepDiff(nullP4);
//          for(size_t ilep=0; ilep<leptons.size(); ilep++) {
//            double varSign( (ivar==LESUP ? 1.0 : (ivar==LESDOWN ? -1.0 : 0.0) ) );
//            int id( abs(leptons[ilep].get("id")) );
//            double sf(1.0);
//            if(id==13) sf=(1.0+varSign*0.01);
//            if(id==11) {
//              if(fabs(leptons[ilep].eta())<1.442) sf=(1.0+varSign*0.02);
//              else                                sf=(1.0-varSign*0.05);
//            }
//            leptonFlux += leptons[ilep];
//            lepDiff += (sf-1)*leptons[ilep];
//          }
//      
//          //clustered flux
//          LorentzVector jetDiff(nullP4), clusteredFlux(nullP4);
//          for(size_t ijet=0; ijet<jets.size(); ijet++)
//            {
//              if(jets[ijet].pt()==0) continue;
//              double jetsf(1.0);
//              if(ivar==JERUP)   jetsf=jets[ijet].getVal("jerup")/jets[ijet].pt();
//              if(ivar==JERDOWN) jetsf=jets[ijet].getVal("jerdown")/jets[ijet].pt();
//              if(ivar==JESUP)   jetsf=jets[ijet].getVal("jesup")/jets[ijet].pt();
//              if(ivar==JESDOWN) jetsf=jets[ijet].getVal("jesdown")/jets[ijet].pt();
//              LorentzVector newJet( jets[ijet] ); newJet *= jetsf;
//              jetDiff       += (newJet-jets[ijet]);
//              clusteredFlux += jets[ijet];
//            }
//          LorentzVector iMet=rawMETP4-jetDiff-lepDiff;
//
//          //unclustered flux
//          if(ivar==UMETUP || ivar==UMETDOWN)
//            {
//              LorentzVector unclusteredFlux=-(iMet+clusteredFlux+leptonFlux);
//              unclusteredFlux *= (ivar==UMETUP ? 1.1 : 0.9); 
//              iMet = -clusteredFlux -leptonFlux - unclusteredFlux;
//            }
//      
//          //save new met
//          newMetsP4[ivar]=iMet;
//        }
//  
//      //all done here
//      return newMetsP4;
//    }
//    
    //
    const reco::Candidate *getGeneratorFinalStateFor(const reco::Candidate *p, bool isSherpa)
    {
      if(p==0) return 0;
      
      const reco::Candidate *prevState=p;
      do{	
	const reco::Candidate *nextState=0;
	int nDaughters = prevState->numberOfDaughters();
	for(int iDaughter=0; iDaughter<nDaughters; iDaughter++)
	  {
	    const reco::Candidate *dau = prevState->daughter(iDaughter);
	    if(dau==0) continue;
	    if(dau->pdgId()!= p->pdgId()) continue;
	    nextState=dau;	   
	    break;
	  }
	if(nextState==0) break;
	if(nextState==prevState) break;
	prevState=nextState;
      }while(1);
      return prevState;
    }

    //
    bool isBhadron(int pdgId)
    {
      int absid=abs(pdgId);
      return ( (absid>= 5122 && absid<= 5554) ||    //baryons
	       (absid>=20513 && absid<=20543) ||    //mesons
	       (absid>=10511 && absid<=10543) || 
	       (absid>=  511 && absid<=  545) ||
	       (absid>=  551 && absid<=  557) ||    //bbar mesons
	       (absid>=10551 && absid<=10557) ||
	       (absid>=100551 && absid<=100557) ||
	       (absid>=110551 && absid<=110557) ||
	       (absid>=200551 && absid<=200557) ||
	       (absid>=210551 && absid<=210557) );
    }

    //cf. https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    std::vector<float> smearJER(float pt, float eta, float genPt)
    {
      std::vector<float> toReturn(3,pt);
      if(genPt<=0) return toReturn;
      
      //
      eta=fabs(eta);
      double ptSF(1.0), ptSF_err(0.06);
      if(eta<0.5)                  { ptSF=1.052; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2)); }
      else if(eta>=0.5 && eta<1.1) { ptSF=1.057; ptSF_err=sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2)); }
      else if(eta>=1.1 && eta<1.7) { ptSF=1.096; ptSF_err=sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2)); }
      else if(eta>=1.7 && eta<2.3) { ptSF=1.134; ptSF_err=sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2)); }
      else if(eta>=2.3 && eta<5.0) { ptSF=1.288; ptSF_err=sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2)); }
      
      toReturn[0]=TMath::Max(0.,(genPt+ptSF*(pt-genPt)));
      toReturn[1]=TMath::Max(0.,(genPt+(ptSF+ptSF_err)*(pt-genPt)));
      toReturn[2]=TMath::Max(0.,(genPt+(ptSF-ptSF_err)*(pt-genPt)));
      return toReturn;
    }

    //
    std::vector<float> smearJES(float pt, float eta, JetCorrectionUncertainty *jecUnc)
    {
      jecUnc->setJetEta(eta);
      jecUnc->setJetPt(pt);
      float relShift=fabs(jecUnc->getUncertainty(true));
      std::vector<float> toRet;
      toRet.push_back((1.0+relShift)*pt);
      toRet.push_back((1.0-relShift)*pt);
      return toRet;
    }


    //
    PuShifter_t getPUshifters(std::vector< float > &Lumi_distr, float puUnc)
    {
      Int_t NBins = Lumi_distr.size();
      TH1F *pu=new TH1F("putmp","",NBins,-0.5,float(NBins)-0.5);
      TH1F *puup=(TH1F *)pu->Clone("puuptmp");
      TH1F *pudown=(TH1F *)pu->Clone("pudowntmp");
      for(size_t i=0; i<Lumi_distr.size(); i++)  pu->SetBinContent(i+1,Lumi_distr[i]);
      
      for(int ibin=1; ibin<=pu->GetXaxis()->GetNbins(); ibin++)
	{
	  Double_t xval=pu->GetBinCenter(ibin);
	  TGraph *gr = new TGraph;
	  for(int ishift=-3; ishift<3; ishift++)
	    {
	      if(ibin+ishift<0) continue;
	      if(ibin+ishift>pu->GetXaxis()->GetNbins()) continue;
	      
	      gr->SetPoint(gr->GetN(),xval+ishift,pu->GetBinContent(ibin+ishift));
	    }
	  if(gr->GetN()>1)
	    {
	      Double_t newval(gr->Eval(xval*(1+puUnc)));
	      pudown->SetBinContent(ibin,newval>0?newval:0.0);
	      newval=gr->Eval(xval*(1-puUnc));
	      puup->SetBinContent(ibin,newval>0?newval:0.0);
	    }
	  delete gr;
	}
      puup->Scale(pu->Integral()/puup->Integral());
      pudown->Scale(pu->Integral()/pudown->Integral());
      std::cout << "getPUshifts will shift average PU by " << puup->GetMean()-pu->GetMean() << " / " << pudown->GetMean()-pu->GetMean() << std::endl; 
      
      puup->Divide(pu);    TGraph *puupWgt = new TGraph(puup);
      pudown->Divide(pu);  TGraph *pudownWgt = new TGraph(pudown);
      delete puup;
      delete pudown;  
      delete pu;
      
      PuShifter_t res(2);
      res[PUDOWN] = pudownWgt;
      res[PUUP]   = puupWgt;
      return res;
    }
    

    //
    Float_t getEffectiveArea(int id,float eta,int cone,TString isoSum)
    {
      Float_t Aeff(1.0);
      if(abs(id)==11){ // electron 
	// PHYS14  https://indico.cern.ch/event/367861/contribution/2/material/slides/0.pdf 
	if(fabs(eta)<0.8)                         Aeff=(cone==3? 0.1013 : 0.180);
	else if(fabs(eta)>0.8 && fabs(eta)<1.3)   Aeff=(cone==3? 0.0988 : 0.200);
	else if(fabs(eta)>1.3 && fabs(eta)<2.0)   Aeff=(cone==3? 0.0572 : 0.150);
	else if(fabs(eta)>2.0 && fabs(eta)<2.2)   Aeff=(cone==3? 0.0842 : 0.190);
	else if(fabs(eta)>2.2 && fabs(eta)<2.5)   Aeff=(cone==3? 0.1530 : 0.210);
      }

      else if(abs(id)==13){ // muon 
	// PHYS14  https://indico.cern.ch/event/367861/contribution/2/material/slides/0.pdf 
	if(fabs(eta)<0.8)                         Aeff=(cone==3? 0.0913 : 0.180);
	else if(fabs(eta)>0.8 && fabs(eta)<1.3)   Aeff=(cone==3? 0.0765 : 0.200);
	else if(fabs(eta)>1.3 && fabs(eta)<2.0)   Aeff=(cone==3? 0.0546 : 0.150);
	else if(fabs(eta)>2.0 && fabs(eta)<2.2)   Aeff=(cone==3? 0.0728 : 0.190);
	else if(fabs(eta)>2.2 && fabs(eta)<2.5)   Aeff=(cone==3? 0.1177 : 0.210);
      }
      
      else if(abs(id)==22){ // photon 
	//https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Recipe_for_regular_users_for_74X
	// Effective areas for the PHYS14, conditions: PU20 bx25 
	if(isoSum=="chIso"){
	  if(fabs(eta)<1.0)                         Aeff=0.0234;
          else if(fabs(eta)>1.0 && fabs(eta)<1.479) Aeff=0.0189;
          else if(fabs(eta)>1.479 && fabs(eta)<2.0) Aeff=0.0171;
          else if(fabs(eta)>2.0 && fabs(eta)<2.2)   Aeff=0.0129;
          else if(fabs(eta)>2.2 && fabs(eta)<2.3)   Aeff=0.0110;
          else if(fabs(eta)>2.3 && fabs(eta)<2.4)   Aeff=0.0074;
          else                                      Aeff=0.0035;
	}
	if(isoSum=="nhIso"){
	  if(fabs(eta)<1.0)                         Aeff=0.0053;
          else if(fabs(eta)>1.0 && fabs(eta)<1.479) Aeff=0.0103;
          else if(fabs(eta)>1.479 && fabs(eta)<2.0) Aeff=0.0057;
          else if(fabs(eta)>2.0 && fabs(eta)<2.2)   Aeff=0.0070;
          else if(fabs(eta)>2.2 && fabs(eta)<2.3)   Aeff=0.0152;
          else if(fabs(eta)>2.3 && fabs(eta)<2.4)   Aeff=0.0232;
          else                                      Aeff=0.1709;
	}
	if(isoSum=="gIso"){
	  if(fabs(eta)<1.0)                         Aeff=0.0780;
          else if(fabs(eta)>1.0 && fabs(eta)<1.479) Aeff=0.0629;
          else if(fabs(eta)>1.479 && fabs(eta)<2.0) Aeff=0.0264;
          else if(fabs(eta)>2.0 && fabs(eta)<2.2)   Aeff=0.0462;
          else if(fabs(eta)>2.2 && fabs(eta)<2.3)   Aeff=0.0740;
          else if(fabs(eta)>2.3 && fabs(eta)<2.4)   Aeff=0.0924;
          else                                      Aeff=0.1484;
	}
      }
      return Aeff;
    }


//   double relIso(llvvLepton lep, double rho){
//      if(abs(lep.id)==11){
//          return (TMath::Max(lep.nhIso03+lep.gIso03-rho*utils::cmssw::getEffectiveArea(11,lep.electronInfoRef->sceta),double(0.))+lep.chIso03)/lep.pt();
//      }else if(abs(lep.id)==13){
//          return (TMath::Max(lep.nhIso04+lep.gIso04-0.5*lep.puchIso04,double(0.))+lep.chIso04)/lep.pt();
//      }else{
//          return -1;
//      }
//   }
    
    
    void getSingleMuTrigEff(const double& pt, const double& abseta, double& muontriggerefficiency){
      // Muon trigger/ID/Iso scale factors for efficiency are taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs                                                                                                                                           
      if(abseta>=0. && abseta <0.9){ // ABCD
	if(pt>=140. /*&& pt<500.*/) muontriggerefficiency=0.98041749810533507;
	if(pt>=25.  && pt<30.)      muontriggerefficiency=0.98372524384334614;
	if(pt>=30.  && pt<35.)      muontriggerefficiency=0.98406344315477012;
	if(pt>=35.  && pt<40.)      muontriggerefficiency=0.98391658181685537;
	if(pt>=40.  && pt<50.)      muontriggerefficiency=0.98345252700570363;
	if(pt>=50.  && pt<60.)      muontriggerefficiency=0.98429177039157478;
	if(pt>=60.  && pt<90.)      muontriggerefficiency=0.98467201842489449;
	if(pt>=90.  && pt<140.)     muontriggerefficiency=0.98091711658069591;
      }
      if(abseta>=0.9 && abseta <1.2){ // ABCD
	if(pt>=140. /*&& pt<500.*/) muontriggerefficiency=0.97127896196175556;
	if(pt>=25.  && pt<30.)      muontriggerefficiency=0.96838127559931908;
	if(pt>=30.  && pt<35.)      muontriggerefficiency=0.96538054889610103;
	if(pt>=35.  && pt<40.)      muontriggerefficiency=0.96696514151670487;
	if(pt>=40.  && pt<50.)      muontriggerefficiency=0.96667958160832501;
	if(pt>=50.  && pt<60.)      muontriggerefficiency=0.96273957552501865;
	if(pt>=60.  && pt<90.)      muontriggerefficiency=0.95952416834753307;
	if(pt>=90.  && pt<140.)     muontriggerefficiency=0.96444182461126438;
      }
      if(abseta>=1.2 && abseta <2.1){ // ABCD
	if(pt>=140. /*&& pt<500.*/) muontriggerefficiency=0.99416866829048334;
	if(pt>=25.  && pt<30.)      muontriggerefficiency=1.0051991254438037;
	if(pt>=30.  && pt<35.)      muontriggerefficiency=1.0013781590159485;
	if(pt>=35.  && pt<40.)      muontriggerefficiency=0.99616640424792002;
	if(pt>=40.  && pt<50.)      muontriggerefficiency=0.99425410141043047;
	if(pt>=50.  && pt<60.)      muontriggerefficiency=0.99054467301217797;
	if(pt>=60.  && pt<90.)      muontriggerefficiency=0.98829374192885855;
	if(pt>=90.  && pt<140.)     muontriggerefficiency=0.98187598993908232;
      }
    }
  
  }


  //
  std::string toLatexRounded(double value, double error, double systError,bool doPowers)
  {
    using namespace std;

    bool ValueWasNull = false;

    if(value==0.0 && error==0.0)return string("");
    if(value==0.0){value=error; ValueWasNull=true;}
    
    if(!doPowers){
      char tmpchar[255];
      if(systError<0)
	sprintf(tmpchar,"$%.0f\\pm%.0f$",value,error);
      else
	sprintf(tmpchar,"$%.0f\\pm%.0f\\pm%.0f$",value,error,systError);
      return string(tmpchar);
    }
    
    double power = floor(log10(value));
    if(power<=-3)     {power=power+3;}
    else if(power>=2) {power=power-2;}
    else              {power=0;}
    
    value = value / pow(10,power);
    error = error / pow(10,power);
    if(systError>=0)systError = systError / pow(10,power);
    int ValueFloating;
    if(systError<0){
      ValueFloating = 1 + std::max(-1*log10(error),0.0);
    }else{
      ValueFloating = 1 + std::max(-1*log10(systError), std::max(-1*log10(error),0.0));
    }
    int ErrorFloating = ValueFloating;
    
 
    if(ValueWasNull){value=0.0;}

    char tmpchar[255];
    if(value<=1E-4){
        //sum in quadrature errors
        double erroSum = 0;
        if(error>0){erroSum+=error*error;}
        if(systError>0){erroSum+=systError*systError;}
        sprintf(tmpchar,"$<%.*f$",ErrorFloating,sqrt(erroSum));
    }else{
       if(power!=0){
         if(systError<0){
           sprintf(tmpchar,"$(%.*f\\pm%.*f)\\times 10^{%g}$",ValueFloating,value,ErrorFloating,error,power);
         }else{
           sprintf(tmpchar,"$(%.*f\\pm%.*f\\pm%.*f)\\times 10^{%g}$",ValueFloating,value,ErrorFloating,error,ErrorFloating,systError,power);
         }
         
       }else{
         if(systError<0){
           sprintf(tmpchar,"$%.*f\\pm%.*f$",ValueFloating,value,ErrorFloating,error);
         }else{
           sprintf(tmpchar,"$%.*f\\pm%.*f\\pm%.*f$",ValueFloating,value,ErrorFloating,error,ErrorFloating,systError);
         }
       }
    }
    return string(tmpchar);
  }

  //
  void TLatexToTex(TString &expr)
  {
    expr = "$"+expr;
    expr += "$";
    expr.ReplaceAll("mu","\\mu"); 
    expr.ReplaceAll("_"," "); 
    expr.ReplaceAll("#","\\");
  }







	// loop on all the lumi blocks for an EDM file in order to count the number of events that are in a sample
	// this is useful to determine how to normalize the events (compute weight)
	unsigned long getMergeableCounterValue(const std::vector<std::string>& urls, std::string counter)
	{
	   unsigned long Total = 0;
	   for(unsigned int f=0;f<urls.size();f++){
	      TFile *file = TFile::Open(urls[f].c_str());      
	      fwlite::LuminosityBlock ls( file );
	      for(ls.toBegin(); !ls.atEnd(); ++ls){
		 fwlite::Handle<edm::MergeableCounter> nEventsTotalCounter;
		 nEventsTotalCounter.getByLabel(ls,counter.c_str());
		 if(!nEventsTotalCounter.isValid()){printf("Invalid nEventsTotalCounterH\n");continue;}
		 Total+= nEventsTotalCounter->value;
	      }
	   }
	   return Total;
	}


  double getTotalNumberOfEvents(std::vector<std::string>& urls, bool fast, bool weightSum)
  {
    double toReturn = 0;
    for(unsigned int f=0;f<urls.size();f++){
       TFile* file = TFile::Open(urls[f].c_str() );
       fwlite::Event ev(file);
       if(!fast){
          for(ev.toBegin(); !ev.atEnd(); ++ev){
             fwlite::Handle< GenEventInfoProduct > genEventInfoHandle;
             genEventInfoHandle.getByLabel(ev, "generator");
             if(!genEventInfoHandle.isValid()){fast=true; break;} //if this object is missing, it's likely missing for the entire sample, move to the fast method
             if(weightSum){toReturn+=genEventInfoHandle->weight();
             }else{                 
               if(genEventInfoHandle->weight()<0){toReturn--;}else{toReturn++;}
             }
          }
       }

       if(fast){
          toReturn += ev.size();          
       }
       delete file;
     }
     return toReturn;
  }


  void getMCPileupDistributionFromMiniAOD(std::vector<std::string>& urls, unsigned int Npu, std::vector<float>& mcpileup)
  {
    mcpileup.clear();
    mcpileup.resize(Npu);
    for(unsigned int f=0;f<urls.size();f++){
       TFile* file = TFile::Open(urls[f].c_str() );
       fwlite::Event ev(file);
       for(ev.toBegin(); !ev.atEnd(); ++ev){
          fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
          puInfoH.getByLabel(ev, "slimmedAddPileupInfo");
          if(!puInfoH.isValid()){printf("collection PileupSummaryInfos with name addPileupInfo does not exist\n"); exit(0);}
          unsigned int ngenITpu = 0;
          for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
             if(it->getBunchCrossing()==0)      { ngenITpu += it->getTrueNumInteractions(); }
          }
          if(ngenITpu>=Npu){printf("ngenITpu is larger than vector size... vector is being resized, but you should check that all is ok!"); mcpileup.resize(ngenITpu+1);}
          mcpileup[ngenITpu]++;
       }
       delete file;
     }
  }


  double getMCPileupDistributionAndTotalEventFromMiniAOD(std::vector<std::string>& urls, unsigned int Npu, std::vector<float>& mcpileup)
  {
    double toReturn=0;
    mcpileup.clear();
    mcpileup.resize(Npu);
    for(unsigned int f=0;f<urls.size();f++){
       TFile* file = TFile::Open(urls[f].c_str() );
       fwlite::Event ev(file);
       for(ev.toBegin(); !ev.atEnd(); ++ev){
          fwlite::Handle< GenEventInfoProduct > genEventInfoHandle;
          genEventInfoHandle.getByLabel(ev, "generator");
          if(!genEventInfoHandle.isValid()){printf("collection generator is not found\n");} //if this object is missing, it's likely missing for the entire sample, move to the fast method
          toReturn+=genEventInfoHandle->weight();


          fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
          puInfoH.getByLabel(ev, "slimmedAddPileupInfo");
          if(!puInfoH.isValid()){printf("collection PileupSummaryInfos with name addPileupInfo does not exist\n"); exit(0);}
          unsigned int ngenITpu = 0;
          for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
             if(it->getBunchCrossing()==0)      { ngenITpu += it->getTrueNumInteractions(); }
          }
          if(ngenITpu>=Npu){printf("ngenITpu is larger than vector size... vector is being resized, but you should check that all is ok!"); mcpileup.resize(ngenITpu+1);}
          mcpileup[ngenITpu]++;
       }
       delete file;
     }
    return toReturn;
  }


  bool isGoodVertex(reco::Vertex& vtx)
  {

    if(vtx.chi2()==0 && vtx.ndof()==0) return false; // Corresponds to the AOD method vtx->isFake()  

    if(vtx.ndof() < 4)            return false;
    if(abs(vtx.z())>24.)          return false;
    if(vtx.position().Rho() >2.0) return false;
    // else
    return true;
  }
  


  void getPileupNormalization(std::vector<float>& mcpileup, double* PUNorm, edm::LumiReWeighting* LumiWeights, utils::cmssw::PuShifter_t PuShifters){
    PUNorm[0]=0; PUNorm[1]=0; PUNorm[2]=0;
    double NEvents=0;
    for(unsigned int i=0;i<mcpileup.size();i++){
      NEvents+=mcpileup[i];
      double puWeight = LumiWeights->weight((int)i);
      PUNorm[0]+=mcpileup[i]*puWeight;
      PUNorm[1]+=mcpileup[i]*puWeight*PuShifters[utils::cmssw::PUDOWN]->Eval(i);
      PUNorm[2]+=mcpileup[i]*puWeight*PuShifters[utils::cmssw::PUUP  ]->Eval(i);
    }
    PUNorm[0]/=NEvents;
    PUNorm[1]/=NEvents;
    PUNorm[2]/=NEvents;
  }



  bool passTriggerPatternsAndGetName(edm::TriggerResultsByName& tr, std::string& pathName, std::string pattern){
     if(edm::is_glob(pattern)){
        std::vector< std::vector<std::string>::const_iterator > matches = edm::regexMatch(tr.triggerNames(), pattern);
        for(size_t t=0;t<matches.size();t++){
           if(tr.accept( matches[t]->c_str() ) ){pathName = *matches[t]; return true;}
        }
     }else{
        if(tr.accept( pattern.c_str() ) ) { pathName = pattern; return true;}
     }
     return false;
  }


  bool passTriggerPatterns(edm::TriggerResultsByName& tr, std::string pattern){
     if(edm::is_glob(pattern)){
        std::vector< std::vector<std::string>::const_iterator > matches = edm::regexMatch(tr.triggerNames(), pattern);
        for(size_t t=0;t<matches.size();t++){
           if(tr.accept( matches[t]->c_str() ) )return true;
        }
     }else{
        if(tr.accept( pattern.c_str() ) ) return true;
     }
     return false;
  }

  bool passTriggerPatterns(edm::TriggerResultsByName& tr, std::string pattern1, std::string pattern2, std::string pattern3, std::string pattern4){
     if(pattern1!="" && passTriggerPatterns(tr, pattern1))return true;
     if(pattern2!="" && passTriggerPatterns(tr, pattern2))return true;
     if(pattern3!="" && passTriggerPatterns(tr, pattern3))return true;
     if(pattern4!="" && passTriggerPatterns(tr, pattern4))return true;
     return false;
  }

  bool passTriggerPatterns(edm::TriggerResultsByName& tr, std::vector<std::string>& patterns){
     for(size_t p=0;p<patterns.size();p++){
        if(passTriggerPatterns(tr, patterns[p]))return true;
     }
     return false;
  }
  
  
   void getHiggsLineshapeFromMiniAOD(std::vector<std::string>& urls, TH1D* hGen){
      if(!hGen)return;
      hGen->Reset();
      for(unsigned int f=0;f<urls.size();f++){
       TFile* file = TFile::Open(urls[f].c_str() );
       fwlite::Event ev(file);
       for(ev.toBegin(); !ev.atEnd(); ++ev){
          reco::GenParticleCollection gen;
          fwlite::Handle< reco::GenParticleCollection > genHandle;
          genHandle.getByLabel(ev, "prunedGenParticles");
          if(genHandle.isValid()){ gen = *genHandle;}

          LorentzVector higgs(0,0,0,0);
	  for(size_t igen=0; igen<gen.size(); igen++){
	     if(!gen[igen].isHardProcess()) continue;
	     if(abs(gen[igen].pdgId())>=11 && abs(gen[igen].pdgId())<=16){ higgs += gen[igen].p4(); }
	  }         
          hGen->Fill(higgs.mass());
       }
       delete file;
     }
  }
 
}
