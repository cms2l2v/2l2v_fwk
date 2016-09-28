#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "TGraphErrors.h"
#include "UserCode/llvv_fwk/interface/th1fmorph.h"

namespace higgs{

  namespace utils{

    //
    EventCategory::EventCategory(int mode):mode_(mode) { }

    EventCategory::~EventCategory() { }
    
    TString EventCategory::GetCategory(pat::JetCollection &jets, LorentzVector &boson)
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
	LorentzVector VBFSyst = jets[0].p4() + jets[1].p4();
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
    double transverseMass(const LorentzVector &visible, const LorentzVector &invisible, bool assumeSameMass){
      if(assumeSameMass){
	LorentzVector sum=visible+invisible;
	double tMass = pow(sqrt(pow(visible.pt(),2)+pow(visible.mass(),2))+sqrt(pow(invisible.pt(),2)+pow(91.188,2)),2);
	tMass-=pow(sum.pt(),2);
	return sqrt(tMass);
      }else{
	double dphi=fabs(deltaPhi(invisible.phi(),visible.phi()));
	return sqrt(2*invisible.pt()*visible.pt()*(1-cos(dphi)));
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


    TGraph* getWeightGraphFromShapes(TH1D* newLineShape, TH1D* originalLineShape, double mH){
      double RMS = originalLineShape->GetRMS();
      int RebinFactor = 1; while((originalLineShape->GetBinWidth(1)*RebinFactor) < RMS/ 8.0){RebinFactor*=2;}

      TH1D* weightsHDenominator = (TH1D*)originalLineShape->Clone("weightsHDenominator");    weightsHDenominator->Rebin(RebinFactor);
//    TH1D* weightsH   = mH>=400?(TH1D*)hSI_nnlo  ->Clone("weightsH")  : (TH1D*)hS_nnlo->Clone("weightsH");   weightsH  ->Rebin(RebinFactor);  weightsH->  Divide(weightsHDenominator);
      TH1D* weightsH   = (TH1D*)newLineShape->Clone("weightsH");   weightsH  ->Rebin(RebinFactor);  weightsH->  Divide(weightsHDenominator);
      for(int x=1;x<=weightsH->GetNbinsX();x++){
         weightsH  ->SetBinContent(x, std::max(0.0,weightsH  ->GetBinContent(x))  );   weightsH  ->SetBinError(x, 0.0  ); //truncate to positive values
         if(fabs(weightsH->GetXaxis()->GetBinCenter(x) - mH)>3*RMS){ weightsH  ->SetBinContent(x, 0.0  );   weightsH  ->SetBinError(x, 0.0  );} //truncate to 3sigma
      }
      TGraph* toReturn=new TGraph(weightsH);
      delete weightsH; 
      delete weightsHDenominator;
      return toReturn;     
   }

    //    
    TH1D* getHistoFromNRfile(std::string histoName, double mass, double Cprime, double BRnew, TFile *nrLineShapesFile){
      //We only have lineshape for Cprime=X and BRnew=0, so we need to find the lineshape equivalent to (Cprime, BRnew) to (Csecond, 0)
      //that is easy because the signal width for (Cprime, Brnew) is SM width * cprimeÂ / (1-BRnew), so the equivalent with can be taken from the pair (Cprime/sqrt(1-BRnew), 0) in order to get the same width
      //some care is needed for the cross-section because it does not scale the same way, but this does not matter since we have the xsection normalized to SM afterward
      //do not allow Csecond to be larger than 1 though

      char nrShapeBuf[100];


      double Csecond = Cprime;
      double BRnew2 = BRnew;
      if(BRnew!=0){
         Csecond = BRnew<=0?Cprime:std::min(1.0, Cprime / sqrt(1-BRnew));
         BRnew2  = 0;        
         printf("BRnew is different than 0, so apply the change of variable (cprime=%f, Brnew=%f) --> (csecond=%f, 0)\n", Cprime, BRnew, Csecond);
      }

      //1st check if is in the file
      if(!nrLineShapesFile){
         printf("LineShapeFile not ok for narrow resonnance\n");  fflush(stdout);
         return NULL;
      }

      //very likely Csecond is not a multiple of 0.1 anymore,
      //we need to morph between to 0.1 multiple
      //first check if Csecond is multiple of 0.1  (in steps of 0.01)
      if( (int(Csecond*10000)/10)%100==0 ){  
  	     sprintf(nrShapeBuf,"%d_%3.1f_%3.1f/%s",int(mass),Csecond,BRnew2, histoName.c_str());
             return (TH1D*)nrLineShapesFile->Get(nrShapeBuf);

     //Csecond is NOT a multiple of 0.1 --> need to morph
      }else{ 
         printf("Csecond %f is NOT a multiple of 0.1 --> %i --> %i\n", Csecond, int(Csecond*1000), int(Csecond*1000)%100);

         //identify neighboring values that are multiple of 0.1
         double CsecondL = (int(Csecond*1000)/100)/10.0;
         double CsecondR = CsecondL+0.1;
         printf("morph the lineshape between %f and %f\n", CsecondL, CsecondR);
         TH1D* HL=getHistoFromNRfile(histoName, mass, CsecondL, BRnew2, nrLineShapesFile);
         TH1D* HR=getHistoFromNRfile(histoName, mass, CsecondR, BRnew2, nrLineShapesFile);
         if(!HL || !HR){printf("Left and Right NR shapes can not be found in file\n"); return NULL;}          
         return th1fmorph("hC","hC", HL, HR, CsecondL, CsecondR, Csecond, 1.0, 0);
      }

      return NULL;
  }


    //    
    TGraph* weightNarrowResonnance(bool isVBF, double mass, double Cprime, double BRnew, TFile *nrLineShapesFile, double& Norm, TString pf){
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
      double Csecond = BRnew<=0?Cprime:std::min(1.0, Cprime / sqrt(1-BRnew));
      double BRnew2  = 0;

      if(BRnew!=0)printf("BRnew is different than 0, so apply the change of variable (cprime=%f, Brnew=%f) --> (csecond=%f, 0)\n", Cprime, BRnew, Csecond);

      //1st check if is in the file
      if(!nrLineShapesFile){
         printf("LineShapeFile not ok for narrow resonnance\n");  fflush(stdout);
         return NULL;
      }

      TH1D* original  = getHistoFromNRfile("mH_S_NR"                  , mass, 1.0, 0.0, nrLineShapesFile);  // SM signal only
//      TH1D* lineshape = getHistoFromNRfile(TString("mH_SI_NR_nnlo")+pf, mass, Csecond, BRnew2, nrLineShapesFile);  // (signal + interference) * NNLOKFactors   
      TH1D* lineshape = NULL;
      if(isVBF){
         //Signal+Interference[h2-h1, h2-Bckg, h1-Bckg] (LO)
         lineshape = getHistoFromNRfile((TString("mH_SI_NR")+pf).Data(), mass, Csecond, BRnew2, nrLineShapesFile);              
      } else {
	 //Signal+Interference[h2-h1, h2-Bckg, h1-Bckg] (LO*KFactor NNLO)
         lineshape = getHistoFromNRfile((TString("mH_SI_NR_nnlo")+pf).Data(), mass, Csecond, BRnew2, nrLineShapesFile);     
      }
      if(!original || !lineshape)return NULL;

      Norm = lineshape->Integral()/original->Integral();

      TGraph* nrGr = getWeightGraphFromShapes(lineshape, original, mass);
      for(int i=0;i<nrGr->GetN();i++){double x, y; nrGr->GetPoint(i, x, y);  if(y<0 || y>1000){y=0;} nrGr->SetPoint(i, x, y);} //make sure weights are not crazy

      //add 20% uncertainty on VBF lineshape
      if(isVBF){
         if(pf.Contains("up"  )){        for(int i=0;i<nrGr->GetN();i++){double x, y; nrGr->GetPoint(i, x, y);  nrGr->SetPoint(i, x, y*1.2);}      }
         if(pf.Contains("down")){        for(int i=0;i<nrGr->GetN();i++){double x, y; nrGr->GetPoint(i, x, y);  nrGr->SetPoint(i, x, y*0.8);}      }
      }


      return nrGr;
   }


    TGraph* weightGGZZContinuum(TFile *nrLineShapesFile, double& Norm, TString pf){
      //1st check if is in the file
      if(!nrLineShapesFile){
         printf("LineShapeFile not ok for ggZZ Continuum\n");  fflush(stdout);
         return NULL;
      }

      char nrShapeBuf[100];
      sprintf(nrShapeBuf,"%d_%3.1f_%3.1f/%s%s",int(1000) ,1.0, 0.0, "kFactors", pf.Data());
      TGraph* nrGr = (TGraph*)nrLineShapesFile->Get(nrShapeBuf);
      for(int i=0;i<nrGr->GetN();i++){double x, y; nrGr->GetPoint(i, x, y);  if(y<0 || y>1000){y=0;} nrGr->SetPoint(i, x, y);} //make sure weights are not crazy

      Norm = 1.0;

      return nrGr;
   }


  }
}
