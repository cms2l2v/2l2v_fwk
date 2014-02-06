#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "TGraphErrors.h"
#include "HiggsAnalysis/CombinedLimit/interface/th1fmorph.h"

namespace higgs{

  namespace utils{

    //
    EventCategory::EventCategory(int mode):mode_(mode) { }

    EventCategory::~EventCategory() { }
    
  
    //
    TString EventCategory::GetCategory(data::PhysicsObjectCollection_t &jets, LorentzVector &boson)
    {
      //jet multiplicity
      int NJets(0);
      for(size_t ijet=0; ijet<jets.size(); ijet++){
	if(jets[ijet].pt()<=30)continue;
	NJets++;
      }
      
      //VBF tag
      bool isVBF(false);
      if(NJets>=2){
	LorentzVector VBFSyst = jets[0] + jets[1];
	double j1eta=jets[0].eta() ;
	double j2eta=jets[1].eta();
	double dEta = fabs(j1eta-j2eta);
	
	int NCentralJet(0), NCentralBoson(0);
	double MaxEta, MinEta;
	if(j1eta<j2eta) { MinEta=j1eta; MaxEta=j2eta;}
	else            { MinEta=j2eta; MaxEta=j1eta;}
	for(size_t ijet=2; ijet<jets.size(); ijet++){
	  float jpt=jets[ijet].pt();
	  float jeta=jets[ijet].eta();
	  if(jpt<30)continue; 
	  if(jeta>MinEta && jeta<MaxEta) NCentralJet++;  
	}
	
	if(boson.eta()>MinEta && boson.eta()<MaxEta) NCentralBoson=1;
	isVBF=( (dEta>4.0) && (VBFSyst.M()>500) && (NCentralJet==0) && (NCentralBoson==1) );
      }
      
      //build classification
      TString cat("");
      switch(mode_)
	{
	case EXCLUSIVEVBF:
	  {
	    cat= isVBF ? "vbf":"novbf";
	    break;
	  }
	case EXCLUSIVE3JETS:
	  {
	    if(NJets==0)      cat="eq0jets";
	    else if(NJets==1) cat="eq1jets";
	    else              cat="geq2jets";
	    break;
	  }
	case EXCLUSIVE3JETSVBF:
	  {
	    if(isVBF)         cat="vbf";
	    else if(NJets==0) cat="eq0jets";
	    else if(NJets==1) cat="eq1jets";
	    else              cat="geq2jets";
	    break;
	  }
	case EXCLUSIVE2JETS:
	  {
	    if(NJets==0)      cat="eq0jets";
	    else              cat="geq1jets";
	    break;
	  }
	case  EXCLUSIVE2JETSVBF:
	  {
	    if(isVBF)         cat="vbf";
	    else if(NJets==0) cat="eq0jets";
	    else              cat="geq1jets";
	    break;
	  }
	default:
	  break;
	}
	
      return cat;
    }
    
    //
    double transverseMass(LorentzVector &visible, LorentzVector &invisible, bool assumeSameMass){
      if(assumeSameMass){
	LorentzVector sum=visible+invisible;
	double tMass = TMath::Power(TMath::Sqrt(TMath::Power(visible.pt(),2)+pow(visible.mass(),2))+TMath::Sqrt(TMath::Power(invisible.pt(),2)+pow(visible.mass(),2)),2);
	tMass-=TMath::Power(sum.pt(),2);
	return TMath::Sqrt(tMass);
      }else{
	double dphi=fabs(deltaPhi(invisible.phi(),visible.phi()));
	return TMath::Sqrt(2*invisible.pt()*visible.pt()*(1-TMath::Cos(dphi)));
      }
      return -1;
    }
    
    //    
    double weightToH125Interference(double mass,double width,TFile *intFile, TString var)
    {
      if(width==0 || intFile==0) return 1.0;
      TString name("weights_ceq"); name+=Int_t(width); 
      if(var!="") name += "_"+var;
      TGraphErrors *gr=(TGraphErrors *)intFile->Get(name);
      if(gr==0) return 1.0;
      return gr->Eval(mass);
    }

    //    
    TGraph* weightNarrowResonnance(std::string SampleName, double m_gen, double mass, double Cprime, double BRnew, TGraph* hLineShapeNominal, TF1 *decayProbPdf, TFile *nrLineShapesFile,TString pf){
      if((Cprime<0 || BRnew<0) || (Cprime==0 && BRnew==0)){
         TGraph* g = new TGraph(2);
         g->SetPoint(0,    0, 1.0);
         g->SetPoint(1, 9999, 1.0);
         return g;
      }


      //We only have lineshape for Cprime=X and BRnew=0, so we need to find the lineshape equivalent to (Cprime, BRnew) to (Csecond, 0)
      //that is easy because the signal width for (Cprime, Brnew) is SM width * cprimeÂ / (1-BRnew), so the equivalent with can be taken from the pair (Cprime/sqrt(1-BRnew), 0) in order to get the same width
      //some care is needed for the cross-section because it does not scale the same way, but this does not matter since we have the xsection normalized to SM afterward
      //do not allow Csecond to be larger than 1 though
      double Csecond = std::min(1.0, Cprime / sqrt(1-BRnew));
      double BRnew2  = 0;

      if(BRnew!=0)printf("BRnew is different than 0, so apply the change of variable (cprime=%f, Brnew=%f) --> (csecond=%f, 0)\n", Cprime, BRnew, Csecond);

      //1st check if is in the file
      if(nrLineShapesFile)
	{
          //very likely Csecond is not a multiple of 0.1 anymore,
          //we need to morph between to 0.1 multiple
          //first check if Csecond is multiple of 0.1
          if( int(Csecond*1000)%100 ==0 ){
 	     char nrShapeBuf[100];
  	     sprintf(nrShapeBuf,"NR_%04d_%02d_%02d_weights%s",int(m_gen),int(Csecond*10),int(BRnew2*10),pf.Data());
	     TGraph *nrGr=(TGraph *)nrLineShapesFile->Get(nrShapeBuf);               
	     if(nrGr){
               for(int i=0;i<nrGr->GetN();i++){double x, y; nrGr->GetPoint(i, x, y);  if(y<0 || y>1000){y=0;} nrGr->SetPoint(i, x, y);}

               return nrGr;
//	       float weight=nrGr->Eval(mass);
//	       if(weight<0) weight=0;
//	       return weight;
	       //float targetNorm=nrGr->Integral();
	       //std::cout << nrShapeBuf<< " " << targetNorm << std::endl;
	       //if(targetNorm==0) targetNorm=1.0; 
	       //float targetProb=nrGr->Eval(mass);
	       //if(targetProb<0) targetProb=0;
	       //float nominalProb=hLineShapeNominal->Eval(mass);
	       //if(nominalProb==0) return 0;
	       //else return targetProb/(targetNorm*nominalProb);
  	     }
          }else{               //Csecond is NOT a multiple of 0.1
             printf("Csecond %f is NOT a multiple of 0.1\n", Csecond);

             //identify neighboring values that are multiple of 0.1
             double CsecondL = (int(Csecond*1000)/100)/10.0;
             double CsecondR = CsecondL+0.1;
             printf("morph the lineshape between %f and %f\n", CsecondL, CsecondR);

             char nrShapeBufL[100];   sprintf(nrShapeBufL,"NR_%04d_%02d_%02d_weights%s",int(m_gen),int(CsecondL*10),int(BRnew2*10),pf.Data());
             TGraph *nrGrL=(TGraph *)nrLineShapesFile->Get(nrShapeBufL);
             char nrShapeBufR[100];   sprintf(nrShapeBufR,"NR_%04d_%02d_%02d_weights%s",int(m_gen),int(CsecondR*10),int(BRnew2*10),pf.Data());
             TGraph *nrGrR=(TGraph *)nrLineShapesFile->Get(nrShapeBufR);
             if(nrGrL && nrGrR){
                for(int i=0;i<nrGrL->GetN();i++){double x, y; nrGrL->GetPoint(i, x, y);  if(y<0 || y>1000){y=0;}  nrGrL->SetPoint(i, x, y);}
                for(int i=0;i<nrGrR->GetN();i++){double x, y; nrGrR->GetPoint(i, x, y);  if(y<0 || y>1000){y=0;}  nrGrR->SetPoint(i, x, y);}

                TH1F* hL = new TH1F("hL", "hL", 1000, 0, std::max(nrGrL->GetX()[nrGrL->GetN()-1], nrGrR->GetX()[nrGrR->GetN()-1]) );
                TH1F* hR = new TH1F("hR", "hR", 1000, 0, std::max(nrGrL->GetX()[nrGrL->GetN()-1], nrGrR->GetX()[nrGrR->GetN()-1]) );
                for(int i=0;i<hL->GetXaxis()->GetNbins();i++){
                   if(nrGrL->Eval(hL->GetXaxis()->GetBinCenter(i))>1000 || nrGrR->Eval(hR->GetXaxis()->GetBinCenter(i))>1000)printf("check AB %f --> %f - %f\n", hL->GetXaxis()->GetBinCenter(i), nrGrL->Eval(hL->GetXaxis()->GetBinCenter(i)), nrGrR->Eval(hR->GetXaxis()->GetBinCenter(i)));
                   float valL = std::max(0.0, nrGrL->Eval(hL->GetXaxis()->GetBinCenter(i)));  if(valL>1000)valL=0;
                   float valR = std::max(0.0, nrGrR->Eval(hR->GetXaxis()->GetBinCenter(i)));  if(valR>1000)valR=0;                   
                   hL->SetBinContent(i, valL);
                   hR->SetBinContent(i, valR);
                }
                TH1F* hC = th1fmorph("hC","hC", hL, hR, CsecondL, CsecondR, Csecond, 1.0, 0);

                TGraph* nrGrC = new TGraph(hC->GetXaxis()->GetNbins());
                for(int i=0;i<hC->GetXaxis()->GetNbins();i++){
                   if(hC->GetBinContent(i)>1000) printf("check C %f --> %f\n", hC->GetXaxis()->GetBinCenter(i), hC->GetBinContent(i));
                   nrGrC->SetPoint(i, hC->GetXaxis()->GetBinCenter(i), hC->GetBinContent(i));
                }

                delete hL;
                delete hR;
                delete hC;

                return nrGrC;
             }
          }
      }
      printf("weight for narrow resonnance not in file\n");
      exit(0);

/*
      //if not found than use relativistic breit-wigner
      double decay_width = -1;
      if(m_gen == 130){    decay_width =   0.00487; 
      }else if(m_gen == 140){    decay_width =   0.00812; 
      }else if(m_gen == 150){    decay_width =   0.01730; 
      }else if(m_gen == 160){    decay_width =   0.08290; 
      }else if(m_gen == 170){    decay_width =   0.38000; 
      }else if(m_gen == 180){    decay_width =   0.63000; 
      }else if(m_gen == 190){    decay_width =   1.04000; 
      }else if(m_gen == 200){    decay_width =   1.43000; 
      }else if(m_gen == 250){    decay_width =   4.04000; 
      }else if(m_gen == 300){    decay_width =   8.43000; 
      }else if(m_gen == 350){    decay_width =  15.20000; 
      }else if(m_gen == 400){    decay_width =  29.20000; 
      }else if(m_gen == 450){    decay_width =  46.95000; 
      }else if(m_gen == 500){    decay_width =  68.00000;
      }else if(m_gen == 550){    decay_width =  93.15000;
      }else if(m_gen == 600){    decay_width = 123.00000;
      }else if(m_gen == 650){    decay_width = 158.00000;
      }else if(m_gen == 700){    decay_width = 199.00000;
      }else if(m_gen == 750){    decay_width = 247.00000;
      }else if(m_gen == 800){    decay_width = 304.00000;
      }else if(m_gen == 850){    decay_width = 371.00000;
      }else if(m_gen == 900){    decay_width = 449.00000;
      }else if(m_gen == 950){    decay_width = 540.00000;
      }else if(m_gen ==1000){    decay_width = 647.00000;
      }
	 
      double OverallXSectionScaleFactor = 1.0;// *hLineShapeNominal->Integral();
      decay_width = decay_width * pow(Csecond,2) / (1-BRnew);      
      //      OverallXSectionScaleFactor *= pow(Csecond,2) * (1-BRnew);
    
      //the CPS shape
      double weight_CPS = hLineShapeNominal->Eval(mass);
    
      //standard BW
      double toReturn(-1);
      if(decayProbPdf==0)
	{
	  double weight_BW = TMath::BreitWigner(mass, m_gen, decay_width);
	  toReturn = OverallXSectionScaleFactor * (weight_BW / weight_CPS);   
	}
      else
	{
	  //other prob function where [0] is the mass of resonance and [1] is the new width
	  // e.g. relativistic BW   
	  //TF1 *f=new TF1("relbw","(2*sqrt(2)*[0]*[1]*sqrt(pow([0],2)*(pow([0],2)+pow([1],2)))/(TMath::Pi()*sqrt(pow([0],2)+sqrt(pow([0],2)*(pow([0],2)+pow([1],2))))))/(pow(pow(x,2)-pow([0],2),2)+pow([0]*[1],2))",0,2000);
	  decayProbPdf->SetParameter(0,m_gen);
	  decayProbPdf->SetParameter(1,decay_width);
	  double weight_BW = decayProbPdf->Eval( mass );
	  toReturn = OverallXSectionScaleFactor * ( weight_BW / weight_CPS);
	}
    
      //   printf("%E / %E = %E\n",weight_BW, hLineShapeNominal->Eval(mass), (weight_BW / hLineShapeNominal->Eval(mass)) );
      if(toReturn<0){toReturn=1.0;}
      if(weight_CPS<=0)toReturn = 0.0;
      return toReturn;
*/
    }
    
  }
}
