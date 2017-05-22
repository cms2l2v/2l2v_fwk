#include "UserCode/llvv_fwk/interface/HiggsUtils.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "UserCode/llvv_fwk/interface/th1fmorph.h"

#include <sstream>

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

      Cprime = 1.0;
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

   TGraph* Get_NNLO_kFactors(){

     //double sF=1;
     TString nnlosf_FileUrl(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/weights/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
     gSystem->ExpandPathName(nnlosf_FileUrl);
     TFile *nnlosf_File = TFile::Open(nnlosf_FileUrl);
     if(!nnlosf_File){
         printf("nnlo k-Factors for Signal is are not found\n");  fflush(stdout);
     }
     TGraph *sFGr = (TGraph*) nnlosf_File->Get("kfactor_Nominal");
     //sF = sFGr->Eval(mass);
     return sFGr;
     nnlosf_File->Close();

   }

   TGraph* Get_CPS_weights( double nominal_mass){

	std::vector<double> Mass;
	std::vector<double> Wgt;
	std::ostringstream strs;
	strs << nominal_mass;
	std::string mass_str = strs.str();
	string File="CPS_weights_M"+mass_str+"GeV.txt";
	TString cps_FileUrl(string(std::getenv("CMSSW_BASE"))+"/src/UserCode/llvv_fwk/data/CPS_weights/"+File);
	gSystem->ExpandPathName(cps_FileUrl);

     	fstream file( cps_FileUrl.Data() );
        string line;
        if ( file.is_open() ){
                 int itrLine=0;
                 while ( !file.eof() ){
			itrLine++;
                        getline(file,line);
                        stringstream file_line(line,ios_base::in); 
			double mass, weight;
			file_line >> mass >> weight;
			Mass.push_back(mass); Wgt.push_back(weight);	
		 }
	}

	double MassArr[Mass.size()]; double WgtArr[Mass.size()];

	for(unsigned int k=0; k<Mass.size(); k++){ MassArr[k]=Mass[k]; WgtArr[k]=Wgt[k];}

     	TGraph *Cps_Gr = new TGraph(Mass.size(),MassArr,WgtArr);
   
	return Cps_Gr;
   }

   float ComputeInterfWeight( Mela& mela, bool isVBF, TString MelaMode, double heavyMass, double heavyWidth, SimpleParticleCollection_t& daughters, SimpleParticleCollection_t& associated, SimpleParticleCollection_t& mothers){

	float Interf_weight=0;

	if(isVBF){
		mela.setInputEvent(&daughters, &associated, &mothers, true);
                if( MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("Continuum") ){

			float Bckg_wg=0; float Sigh2_wg=0; float All_wg=0;

			//Bckg Only
                	mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::JJVBF_S);
               		mela.computeProdDecP( Bckg_wg, false);

			//Heavy Boson Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);	
			mela.computeProdDecP( Sigh2_wg, false);	


			//Bckg+Heavy Boson
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);	
                        mela.computeProdDecP( All_wg, false); 

			Interf_weight = All_wg - Sigh2_wg - Bckg_wg;

		} else if( MelaMode.Contains("Interf") && MelaMode.Contains("h1") && MelaMode.Contains("Continuum") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float All_wg=0;

			//Bckg Only
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::JJVBF_S);
                        mela.computeProdDecP( Bckg_wg, false); 

                        //Light Higgs Only 
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( Sigh1_wg, false); 

                        //Bckg+Light Higgs 
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S); 
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( All_wg, false); 
                        
                        Interf_weight = All_wg - Sigh1_wg - Bckg_wg;
                
		} else if( MelaMode.Contains("Interf") && MelaMode.Contains("Full") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float Sigh2_wg=0; float All_wg=0;

			//Bckg Only
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::JJVBF_S);
                        mela.computeProdDecP( Bckg_wg, false); 

			//Light Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( Sigh1_wg, false); 

			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S); 
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        mela.computeProdDecP( Sigh2_wg, false); 

			//Bckg+Heavy Higgs+Light Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S); 
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);
			mela.computeProdDecP( All_wg, false);

			Interf_weight = All_wg - Sigh1_wg - Sigh2_wg - Bckg_wg;

                } else if( MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("h1") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float Sigh2_wg=0; float All_wg=0;
                        float All_h2Bckg=0; float All_h1Bckg=0; float h2Bckg=0; float h1Bckg=0;

			//Bckg Only
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::JJVBF_S);
                        mela.computeProdDecP( Bckg_wg, false);

			//Light Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S); 
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( Sigh1_wg, false);

			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        mela.computeProdDecP( Sigh2_wg, false);

			//Bckg+Light Higgs+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);	
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1); 
                        mela.computeProdDecP( All_wg, false);

			//Bckg+Light Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( All_h1Bckg, false);

			//Bckg+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        mela.computeProdDecP( All_h2Bckg, false);

                        h2Bckg = All_h2Bckg - Sigh2_wg - Bckg_wg;
                        h1Bckg = All_h1Bckg - Sigh1_wg - Bckg_wg;
                        Interf_weight = All_wg - Sigh1_wg - Sigh2_wg - Bckg_wg - h2Bckg - h1Bckg;

                } else if(MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("h1") && MelaMode.Contains("Continuum")){

			float All_wg=0; float All_h1Bckg=0; float Sigh2_wg=0;	
			//Bckg+Light Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeProdDecP( All_h1Bckg, false);

			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);       
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        mela.computeProdDecP( Sigh2_wg, false);

			//Light Higgs+Heavy Higgs+Continuum
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S); 
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);
                        mela.computeProdDecP( All_wg, false);

			Interf_weight = All_wg - All_h1Bckg - Sigh2_wg;
		}

	} else {
                mela.setInputEvent(&daughters, 0, &mothers, true);
                if( MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("Continuum") ){

                        float Bckg_wg=0; float Sigh2_wg=0; float All_wg=0;
                       
			//Bckg Only 
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
                        mela.computeP( Bckg_wg, false);
                       
			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0); 
			mela.computeP( Sigh2_wg, false);
			  
			//Bckg+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0); 
                        mela.computeP( All_wg, false); 

                        Interf_weight = All_wg - Sigh2_wg - Bckg_wg;

                } else if( MelaMode.Contains("Interf") && MelaMode.Contains("h1") && MelaMode.Contains("Continuum") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float All_wg=0;
                       
			//Bckg Only 
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
                        mela.computeP( Bckg_wg, false);
                       
			//Light Higgs Only 
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0); 
                        mela.computeP( Sigh1_wg, false);
                       
			//Bckg+Light Higgs Only 
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0); 
                        mela.computeP( All_wg, false); 
                        
                        Interf_weight = All_wg - Sigh1_wg - Bckg_wg;

                } else if( MelaMode.Contains("Interf") && MelaMode.Contains("Full") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float Sigh2_wg=0; float All_wg=0;

			//Bckg Only
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
                        mela.computeP( Bckg_wg, false);

			//Light Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0); 
                        mela.computeP( Sigh1_wg, false);

			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0); 
                        mela.computeP( Sigh2_wg, false);

			//Bckg+Light Higgs+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);  
                        mela.computeP( All_wg, false);

                        Interf_weight = All_wg - Sigh1_wg - Sigh2_wg - Bckg_wg;

                } else if( MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("h1") ){

                        float Bckg_wg=0; float Sigh1_wg=0; float Sigh2_wg=0; float All_wg=0; 
			float All_h2Bckg=0; float All_h1Bckg=0; float h2Bckg=0; float h1Bckg=0;

			//Bckg Only
                        mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
                        mela.computeP( Bckg_wg, false);

			//Light Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeP( Sigh1_wg, false); 

			//Heavy Higgs Only
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0); 
                        mela.computeP( Sigh2_wg, false); 

			//Bckg+Light Higgs+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);           
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);
                        mela.computeP( All_wg, false);

			//Bckg+Light Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeP( All_h1Bckg, false);

			//Bckg+Heavy Higgs
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        mela.computeP( All_h2Bckg, false);

			h2Bckg = All_h2Bckg - Sigh2_wg - Bckg_wg;
			h1Bckg = All_h1Bckg - Sigh1_wg - Bckg_wg;
                        Interf_weight = All_wg - Sigh1_wg - Sigh2_wg - Bckg_wg - h2Bckg - h1Bckg;

		} else if(MelaMode.Contains("Interf") && MelaMode.Contains("h2") && MelaMode.Contains("h1") && MelaMode.Contains("Continuum")){

                        float All_wg=0; float All_h1Bckg=0; float Sigh2_wg=0;

			//Light higgs+Continuum
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG); 
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.computeP( All_h1Bckg, false);

			//Heavy Higgs
                        mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                        mela.computeP( Sigh2_wg, false);

			//Light Higgs+Heavy Higgs+Continuum
                        mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1); 
                        mela.computeP( All_wg, false);

                        Interf_weight = All_wg - All_h1Bckg - Sigh2_wg;

		} 

	}

	return Interf_weight;

   }

   double weightNarrowResonnance_MELA( Mela& mela, bool isVBF, TString MelaMode, double CP, double heavyMass, fwlite::Event& eV){

        //Mela mela( 13, heavyMass, TVar::DEBUG);

        fwlite::Handle< LHEEventProduct > lheEv;
        lheEv.getByLabel(eV, "externalLHEProducer");

        //Weight to reweight the MELA shape to the real cross-section
        double continuum_weight=1.;
        continuum_weight = weightContinuum_MELA(isVBF,CP,heavyMass);
	
	//Fill a Map with Mass and Width SM Like
	double heavyWidth=0; float weightSM=1; float weightMELA=1; float finalweight=1; float cpsweight=1;
	float  propFixedW=1; float propCPSW=1; 

	std::map< double, double>  SM_Info;
	SM_Info[200]=1.43;  SM_Info[300]=8.43;  SM_Info[400]=29.3; 
	SM_Info[500]=68;    SM_Info[600]=123;   SM_Info[700]=199;
	SM_Info[800]=304;   SM_Info[900]=499;   SM_Info[1000]=647;
	SM_Info[1500]=1500; SM_Info[2000]=2000; SM_Info[2500]=2500;
	SM_Info[3000]=3000;

	heavyWidth=SM_Info[heavyMass];

        SimpleParticleCollection_t daughters, mothers, associated; // associated;
        TLorentzVector Higgs;
	std::vector< TLorentzVector> Partons, AssPartons, Lep;
	TLorentzVector AssGluon;

	bool isVBF_gqInitial(false);

	//Loop on particles and fill SimpleParticleCollection_t 
        for(int k=0; k<lheEv->hepeup().NUP; k++){

	    //if( isVBF && k==5 ) continue;
	    int PdgId=0.; int Status=0.;
	    PdgId=lheEv->hepeup().IDUP.at(k); 
	    Status=lheEv->hepeup().ISTUP.at(k);
	    double Px=lheEv->hepeup().PUP.at(k)[0]; double Py=lheEv->hepeup().PUP.at(k)[1]; 
	    double Pz=lheEv->hepeup().PUP.at(k)[2]; double  E=lheEv->hepeup().PUP.at(k)[3]; 
	    TLorentzVector check( Px, Py, Pz, E); 
            if( (abs(PdgId)<7.0 || PdgId==21.0) && Status<0.0 ){	
                TLorentzVector partons( Px, Py, Pz, E);
		if (abs(PdgId)<7.0 && isVBF) mothers.push_back( SimpleParticle_t( PdgId, partons)); //Filling Infos
                else mothers.push_back(SimpleParticle_t(0, partons)); //Else fill gluons as 0 (unknown parton) in case the initial state is qg in ggF, or qg or gg in VBF
		if( isVBF && PdgId==21.0 && Status<0.0) isVBF_gqInitial = true;
	    } else if ( (abs(PdgId)<7.0 || PdgId==21.0) && Status>0.0){
                TLorentzVector extra_partons( Px, Py, Pz, E);
		if (abs(PdgId)<7.0 && isVBF) AssPartons.push_back( extra_partons );
		else if( PdgId==21.0 && isVBF ) AssGluon.SetPxPyPzE( Px, Py, Pz, E);
                //if (abs(PdgId)<7.0 && isVBF) associated.push_back( SimpleParticle_t( PdgId, extra_partons));
                //else if(abs(PdgId)==21.0 && isVBF) associated.push_back(SimpleParticle_t(0, extra_partons)); 
	    } else if ( abs(PdgId)==11.0 || abs(PdgId)==12.0 || abs(PdgId)==13.0 || abs(PdgId)==14.0 || abs(PdgId)==15.0 || abs(PdgId)==16.0 ){
		TLorentzVector lepP( Px, Py, Pz, E);
		daughters.push_back( SimpleParticle_t( PdgId, lepP)); //Filling Infos
	    } else if (  abs(PdgId)==25.0 ){
                Higgs.SetPxPyPzE( Px, Py, Pz, E);
	    }

	}

	unsigned int position = 10000;
	if(isVBF){
		float mindR = 100;
		for(unsigned int s=0; s<AssPartons.size(); s++){
			float dR = AssPartons[s].DeltaR(AssGluon);
			if( dR < mindR){
				mindR = dR;
				position = s;
			}	
		}
	}

	if(isVBF){
		for(unsigned int ns=0; ns<AssPartons.size(); ns++){
			if(position==ns){
				TLorentzVector NewAss =  AssPartons[ns]+AssGluon;
				associated.push_back( SimpleParticle_t( 0, NewAss));
			}else{
				associated.push_back( SimpleParticle_t( 0, AssPartons[ns]));
			}
		}	
	}
	
	std::sort( associated.begin(), associated.end(), utils::sort_CandidatesByPt_V2);

        mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ); //Mela Candidate mode initialized

	if(isVBF){ 
		mela.setInputEvent(&daughters, &associated, &mothers, true);
		mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
        	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
		mela.computeProdDecP( weightSM, false);
		//if(weightSM==0.)TUtil::PrintCandidateSummary(mela.getCurrentCandidate());	
	}else{ 
		mela.setInputEvent(&daughters, 0, &mothers, true);
		mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
        	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
                mela.computeP( weightSM, false);
	} 

	//CPS pole scheme reweight
	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
	mela.getXPropagator(TVar::FixedWidth, propFixedW);
	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0);
	mela.getXPropagator(TVar::CPS, propCPSW);
	cpsweight= propFixedW/propCPSW;
	
        //BSM reweighiting
        mela.resetInputEvent();

	heavyWidth=CP;
	//heavyWidth=heavyWidth*CP*CP;

	if( !MelaMode.Contains("Interf") ){
        	if(isVBF){
                	mela.setInputEvent(&daughters, &associated, &mothers, true);
                	if(MelaMode.Contains("Continuum")){
				mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::JJVBF_S); 
			} else if(MelaMode.Contains("Bckg")){
                                mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S);  
                                mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);	
			} else if(MelaMode.Contains("Sigh2")){	
                        	mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S); 
        			mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0); 
			} else if(MelaMode.Contains("Sigh1")){
                        	mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF_S);
                        	mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);		
			} else if(MelaMode.Contains("All")){
                        	mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::JJVBF_S); 
                        	mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1);	
			}
       			mela.computeProdDecP( weightMELA, false);
	 	}else{
                	mela.setInputEvent(&daughters, 0, &mothers, true);	
			if(MelaMode.Contains("Continuum")){
                		mela.setProcess( TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
			} else if(MelaMode.Contains("Bckg")){
                                mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                                mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);	
			} else if(MelaMode.Contains("Sigh2")){ 
                        	mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 0); 
			} else if(MelaMode.Contains("Sigh1")){ 
                        	mela.setProcess( TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
                        	mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
			} else if(MelaMode.Contains("All")){	
                        	mela.setProcess( TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
                        	mela.setMelaHiggsMassWidth( 125, 4.07e-3, 0);
                        	mela.setMelaHiggsMassWidth( heavyMass, heavyWidth, 1); 
                	}
			mela.computeP( weightMELA, false);
        	}
		
	} else if( MelaMode.Contains("Interf") ){
		weightMELA = ComputeInterfWeight( mela, isVBF, MelaMode, heavyMass, heavyWidth, daughters, associated, mothers);
	}
 
	mela.resetInputEvent();

	if(weightSM==0){ finalweight=0; }
	else if(isVBF_gqInitial){ finalweight=0; }
        else{ finalweight=(weightMELA/weightSM)*cpsweight*continuum_weight;}

	/*if(isnan(finalweight)){
		printf(" \n");
		printf("Particle Size: %5i WeightMELA: %20.18f WeightSM: %20.18f Continuum: %20.18f \n", lheEv->hepeup().NUP, weightMELA, weightSM, continuum_weight);
        	for(int k=0; k<lheEv->hepeup().NUP; k++){

            		//if( isVBF && k==5 ) continue; 
            		int PdgId=0.; int Status=0.;
            		PdgId=lheEv->hepeup().IDUP.at(k);
            		Status=lheEv->hepeup().ISTUP.at(k);
            		double Px=lheEv->hepeup().PUP.at(k)[0]; double Py=lheEv->hepeup().PUP.at(k)[1];
            		double Pz=lheEv->hepeup().PUP.at(k)[2]; double  E=lheEv->hepeup().PUP.at(k)[3];
            		TLorentzVector check( Px, Py, Pz, E);
            		printf("Particle: %4i Mass: %10.5f Status: %4i Px: %10.5f Py: %10.5f Pz: %10.5f E: %10.5f Eta: %6.3f \n", PdgId, check.M(), Status, check.Px(), check.Py(), check.Pz(), check.E(), check.Eta());
	  	}
	}*/

	
        return finalweight;

    }

 double weightContinuum_MELA( bool isVBF, double CP, double heavyMass){

	double continuumWeight=1;
	std::map< double, std::map<double,double> > cpWeight_MapVBF;
        std::map< double, std::map<double,double> > cpWeight_MapggH;

	//GGH Continuum weights, Wider Resonance
        cpWeight_MapggH[100.0][200] = 0.0087132425794885122; 
        cpWeight_MapggH[100.0][300] = 0.0373777202252024363;
        cpWeight_MapggH[100.0][400] = 0.0359384479550159119;
        cpWeight_MapggH[100.0][500] = 0.0155390105379753559;
        cpWeight_MapggH[100.0][600] = 0.0057398053414657946;
        cpWeight_MapggH[100.0][700] = 0.0028350558580945489;
        cpWeight_MapggH[100.0][800] = 0.0013594444780210208;
        cpWeight_MapggH[100.0][900] = 0.0009418577837818124;
        cpWeight_MapggH[100.0][1000] = 0.0005459666328769485;
        cpWeight_MapggH[100.0][1500] = 0.0000742221719463059;
        cpWeight_MapggH[100.0][2000] = 0.0000273049684146098;
        cpWeight_MapggH[100.0][2500] = 0.0000110560511182888;
        cpWeight_MapggH[100.0][3000] = 0.0000051082914848699;

        cpWeight_MapggH[10.0][200] = 0.0087132425794885122;
        cpWeight_MapggH[10.0][300] = 0.0373777202252024363;
        cpWeight_MapggH[10.0][400] = 0.0359384479550159119;
        cpWeight_MapggH[10.0][500] = 0.0155390105379753559; 
        cpWeight_MapggH[10.0][600] = 0.0057398053213052885;
        cpWeight_MapggH[10.0][700] = 0.0028350558580945489;
        cpWeight_MapggH[10.0][800] = 0.0013594444780210208;
        cpWeight_MapggH[10.0][900] = 0.0009418577837818124;
        cpWeight_MapggH[10.0][1000] = 0.0005459666328769485;
        cpWeight_MapggH[10.0][1500] = 0.0000742221719463059;
        cpWeight_MapggH[10.0][2000] = 0.0000273049684146098;
        cpWeight_MapggH[10.0][2500] = 0.0000110560511182888;
        cpWeight_MapggH[10.0][3000] = 0.0000051082914848699;

        cpWeight_MapggH[5.0][200] = 0.0087132425794885122;
        cpWeight_MapggH[5.0][300] = 0.0373777202252024363;
        cpWeight_MapggH[5.0][400] = 0.0359384479550159119;
        cpWeight_MapggH[5.0][500] = 0.0155390105379753559;
        cpWeight_MapggH[5.0][600] = 0.0057398053213052885;
        cpWeight_MapggH[5.0][700] = 0.0028350558580945489;
        cpWeight_MapggH[5.0][800] = 0.0013594444780210208;
        cpWeight_MapggH[5.0][900] = 0.0009418577837818124;
        cpWeight_MapggH[5.0][1000] = 0.0005459666328769485;
        cpWeight_MapggH[5.0][1500] = 0.0000742221719463059;
        cpWeight_MapggH[5.0][2000] = 0.0000273049684146098;
        cpWeight_MapggH[5.0][2500] = 0.0000110560511182888;
        cpWeight_MapggH[5.0][3000] = 0.0000051082914848699;

	//Filling ContinuumWeights VBF
        cpWeight_MapVBF[100.0][200] = 0.0000988541685841062; 
        cpWeight_MapVBF[100.0][300] = 0.0156497585728107742;
        cpWeight_MapVBF[100.0][400] = 0.0065885676955840574;
        cpWeight_MapVBF[100.0][500] = 0.0029627096904723373;
        cpWeight_MapVBF[100.0][600] = 0.0028489759402046014;
        cpWeight_MapVBF[100.0][700] = 0.0017261912468847672;
        cpWeight_MapVBF[100.0][800] = 0.0011086243396082442;
        cpWeight_MapVBF[100.0][900] = 0.0007438900405574646;
        cpWeight_MapVBF[100.0][1000] = 0.0007343592746045676;
        cpWeight_MapVBF[100.0][1500] = 0.0002057277715874602;
        cpWeight_MapVBF[100.0][2000] = 0.0001122712310492367;
        cpWeight_MapVBF[100.0][2500] = 0.0000927441982705444;
        cpWeight_MapVBF[100.0][3000] = 0.0000593464387931730;

        cpWeight_MapVBF[10.0][200] = 0.0000988541685841062;
        cpWeight_MapVBF[10.0][300] = 0.0156497585728107742;
        cpWeight_MapVBF[10.0][400] = 0.0065885676955840574;
        cpWeight_MapVBF[10.0][500] = 0.0029627096904723373;
        cpWeight_MapVBF[10.0][600] = 0.0028489759402046014;
        cpWeight_MapVBF[10.0][700] = 0.0017261912468847672;
        cpWeight_MapVBF[10.0][800] = 0.0011086243396082442;
        cpWeight_MapVBF[10.0][900] = 0.0007438900405574646;
        cpWeight_MapVBF[10.0][1000] = 0.0007343592746045676;
        cpWeight_MapVBF[10.0][1500] = 0.0002057277715874602;
        cpWeight_MapVBF[10.0][2000] = 0.0001122712310492367;
        cpWeight_MapVBF[10.0][2500] = 0.0000927441982705444;
        cpWeight_MapVBF[10.0][3000] = 0.0000593464387931730;

        cpWeight_MapVBF[5.0][200] = 0.0000988541685841062;
        cpWeight_MapVBF[5.0][300] = 0.0156497585728107742;
        cpWeight_MapVBF[5.0][400] = 0.0065885676955840574;
        cpWeight_MapVBF[5.0][500] = 0.0029627096904723373;
        cpWeight_MapVBF[5.0][600] = 0.0028489759402046014;
        cpWeight_MapVBF[5.0][700] = 0.0017261912468847672;
        cpWeight_MapVBF[5.0][800] = 0.0011086243396082442;
        cpWeight_MapVBF[5.0][900] = 0.0007438900405574646;
        cpWeight_MapVBF[5.0][1000] = 0.0007343592746045676;
        cpWeight_MapVBF[5.0][1500] = 0.0002057277715874602;
        cpWeight_MapVBF[5.0][2000] = 0.0001122712310492367;
        cpWeight_MapVBF[5.0][2500] = 0.0000927441982705444;
        cpWeight_MapVBF[5.0][3000] = 0.0000593464387931730;

	if(isVBF){ continuumWeight = cpWeight_MapVBF[CP][heavyMass]; } else { continuumWeight = cpWeight_MapggH[CP][heavyMass]; }
	return continuumWeight;

 }

  }
}
