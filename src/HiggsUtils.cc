#include "UserCode/llvv_fwk/interface/HiggsUtils.h"

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
    double weightNarrowResonnance(std::string SampleName, double m_gen, double mass, double Cprime, double BRnew, TGraph* hLineShapeNominal, TF1 *decayProbPdf){
      if((Cprime<0 || BRnew<0) || (Cprime==0 && BRnew==0)) return 1.0;
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
	 
      double OverallXSectionScaleFactor = 1.0;//*hLineShapeNominal->Integral();
      decay_width = decay_width * pow(Cprime,2) / (1-BRnew);      
      //      OverallXSectionScaleFactor *= pow(Cprime,2) * (1-BRnew);
    
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
    }
    
  }
}
