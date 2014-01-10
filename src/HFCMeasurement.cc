#include "UserCode/llvv_fwk/interface/HFCMeasurement.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/ProfileInspector.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/RooStatsUtils.h"

#include "RooNumIntConfig.h"
#include "RooNLLVar.h"
#include "RooConstVar.h"
#include "RooExtendPdf.h"

#include "TGraphErrors.h"
#include "TSystem.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;


HFCMeasurement::HFCMeasurement(int fitType,TString fitConfig, TString wpConfig, int maxJets)
{ 

  fitType_=fitType;
  switch(fitType_)
   {
   case FIT_EB:              fitTypeTitle_="#varepsilon_{b}";                          fitTypeName_="effb";                        break;
   case FIT_R_AND_EB:        fitTypeTitle_="R;#varepsilon_{b}";                        fitTypeName_="rvseffb";                     break;
   case FIT_R_AND_MU:        fitTypeTitle_="R;#mu=#sigma/#sigma_{th}";                 fitTypeName_="rvsmu";                       break;
   case FIT_EB_AND_MU:        fitTypeTitle_="#varepsilon_{b};#mu=#sigma/#sigma_{th}";  fitTypeName_="ebvsmu";                      break;
   case FIT_VTB:             fitTypeTitle_="|V_{tb}|";                                 fitTypeName_="vtb";                         break;
   case FIT_GAMMAT:          fitTypeTitle_="#Gamma_{t} [GeV]";                         fitTypeName_="gammat";                      break;
   default:                  fitTypeTitle_="R=B(t#rightarrow Wb)/B(t#rightarrow Wq)";  fitTypeName_="r";           fitType_=FIT_R; break;
    }

  maxJets_=maxJets;

  ws_ = new RooWorkspace("w");
  mc_ = 0;
  data_=0;

  //
  //parameters of interest
  //
  RooArgSet poi;
  char expBuf[500];
  float cenR(1.0), minR(0.0), maxR(2.0);
  if(fitType_==FIT_R_AND_EB)                            { minR=0.95;  maxR=1.05;            }
  if(fitType_==FIT_EB || fitType_==FIT_EB_AND_MU)       { cenR=1.0;   minR=1.0;   maxR=1.0; }
  if(fitType_==FIT_VTB)                                 { minR=0.0;   maxR=1.08;            }
  if(fitType_==FIT_GAMMAT)                              { cenR=1.34;  minR=0.0;   maxR=5.0; }
  sprintf(expBuf,"rb[%f,%f,%f]",cenR,minR,maxR);
  ws_->factory(expBuf);
  if(minR!=maxR) poi.add( *ws_->var("rb") );
  if(fitType_==FIT_VTB)  
    {
      ws_->factory("FormulaVar::r('@0*@0',{rb})");
    }
  else if(fitType_==FIT_GAMMAT)
    {
      instantiateSingleTopContribution(ws_);
    }
  else
    {
      ws_->factory("FormulaVar::r('@0',{rb})");
    }

  float minSFb(1.0),maxSFb(1.0);
  if(fitType_==FIT_R_AND_EB)                       { minSFb=0.9;  maxSFb=1.1; }
  if(fitType_== FIT_EB || fitType_==FIT_EB_AND_MU) { minSFb=0.01; maxSFb=1.4; }
  sprintf(expBuf,"sfeb[1.0,%f,%f]",minSFb,maxSFb);
  ws_->factory(expBuf);
  if(minSFb!=maxSFb) poi.add( *ws_->var("sfeb") );

  float minMu(1.0),maxMu(1.0);
  if(fitType_==FIT_R_AND_MU || fitType_==FIT_EB_AND_MU) { minMu=0; maxMu=2.0; }
  sprintf(expBuf,"mu[1.0,%f,%f]",minMu,maxMu);
  ws_->factory(expBuf);
  if(minMu!=maxMu) poi.add( *ws_->var("mu") );
  
  ws_->defineSet("poi",poi);  
  parseFitConfig(fitConfig);
  parseFitConfig(wpConfig);
  initHFCModel();
}


//
void HFCMeasurement::parseFitConfig(TString url)
{
  //check if path is ok
  gSystem->ExpandPathName(url);
  if(gSystem->AccessPathName(url)) return;

  char expBuf[500];  //free buffer	  
  RooArgSet nuis,constr,globalObs;
  
  //check if Vtd/Vts information is already available
  if(ws_->var("vtdvts")==0)
    {
      RooGaussian *nuisConstr=(RooGaussian *)ws_->factory("Gaussian::vtdvts_constr(vtdvts0[0,-5,5],vtdvts[0,-5,5],1.0)");
      constr.add( *nuisConstr );
      
      RooRealVar *globalObsVar=ws_->var("vtdvts0");
      globalObsVar->setConstant(true);
      globalObsVar->SetTitle("|V_{td}|/|V_{ts}|(0)");
      globalObs.add( *globalObsVar );
      
      RooRealVar *nuisVar = ws_->var("vtdvts");
      nuisVar->SetTitle("|V_{td}|/|V_{ts}|");
      nuis.add( *nuisVar );

      ws_->factory("FormulaVar::vtdvts2('pow(0.211*(1+0.0061*@0),2)',{vtdvts})");
    }
  
  //check if tW production correction hasn't been instantiated and if R is already available
  if(ws_->function("tWprod")==0 && ws_->function("r")!=0)
    {
      ws_->factory("FormulaVar::tWprod('(1.0-max(min(@0,1.0),0.0))/(1+@1)*(@1*13.37+2.23)+max(min(@0,1.0),0.0)',{r,vtdvts2})");
    }

  //read parameters from file
  JSONWrapper::Object jsonF(url.Data(), true);

  std::vector<JSONWrapper::Object> params=jsonF.daughters();
  for(size_t iparam=0; iparam<params.size(); iparam++)
    {
      //special case for the tagger name
      if(jsonF.key[iparam]=="tagger")
	{
	  wp_=jsonF["tagger"].toString();
	  continue;
	}
      else if(jsonF.key[iparam]=="sample")
	{
	  sampleType_=jsonF["sample"].toString();
	  continue;
	}

      std::vector<JSONWrapper::Object> &dau = params[iparam].daughters();
      for(size_t icat=0; icat<dau.size(); icat++)
	{
	  JSONWrapper::Object &descript=dau[icat];
	  string param=jsonF.key[iparam];
	  if(params[iparam].key[icat]!="")
	    {
	      param+="_"+params[iparam].key[icat];
	      if(params[iparam].key[icat].find("_")==string::npos) sampleCats_.insert(params[iparam].key[icat]);
	    }
	  std::vector<std::string> uncs = descript.key;

	  //for each parameter instantiate a formula, per category of the type: x\prod_i(1+theta_i)
	  //where theta_i are the nuisances which will be varied with a given prior
	  RooRealVar *cenVal=ws_->var((param+"_cen").c_str());
	  if(cenVal==0)
	    {
	      string formula("@0");
	      RooArgList varsInFormula;

	      //central value
	      float val=descript["val"].toDouble();
	      sprintf(expBuf,"%s_cen[%f,%f,%f]",param.c_str(),val,val,val);
	      cenVal=(RooRealVar *)ws_->factory(expBuf);
	      cenVal->setConstant(true);
	      varsInFormula.add(*cenVal);

	      int varCntr(0);
	
	      //tW correction
	      if(jsonF.key[iparam]=="kst")
		{
		  cout << "Attaching tW correction factor" << endl;
		  varCntr++; 
		  formula += "*@1";
		  varsInFormula.add( *(ws_->function("tWprod")) );
		}

	      //signal strength correction
	      if(jsonF.key[iparam]=="ftt")
		{
		  cout << "Attaching signal strength correction factor" << endl;
		  varCntr++; 
		  formula += "*@1";
		  varsInFormula.add( *(ws_->var("mu")) );
		}
	      
	      //modifiers
	      for(size_t iunc=0; iunc<uncs.size(); iunc++)
		{
		  if(uncs[iunc]=="val") continue;
		  if(uncs[iunc]=="epsb" && (fitType_==FIT_EB || fitType_==FIT_EB_AND_MU)) continue;

		  //add new constraint if required
		  RooRealVar *nuisVar=ws_->var(uncs[iunc].c_str());
		  RooRealVar *modVar=ws_->var((uncs[iunc]+"_sigma").c_str());

		  if(nuisVar==0)
		    {
		      sprintf(expBuf,"Gaussian::%s_constr(%s0[0,-5,5],%s[0,-5,5],1.0)",uncs[iunc].c_str(), uncs[iunc].c_str(), uncs[iunc].c_str());
		      //is this a bug?
		      //sprintf(expBuf,"Gaussian::%s_constr(%s0[0,-5,5],%s[0,-5,5].0,1.0)",uncs[iunc].c_str(), uncs[iunc].c_str(), uncs[iunc].c_str());
		      RooGaussian *nuisConstr=(RooGaussian *)ws_->factory(expBuf);
		      constr.add( *nuisConstr );
		      
		      //this is a global observable will be used for toys
		      //cf. https://indico.desy.de/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=6083
		      //cf. https://twiki.cern.ch/twiki/bin/view/RooStats/RooStatsTutorialsAugust2012#Create_Poisson_Counting_model
		      RooRealVar *globalObsVar = ws_->var( (uncs[iunc]+"0").c_str() );         
		      globalObsVar->setConstant(true);
		      globalObsVar->SetTitle( (getNuisanceTitle(uncs[iunc])+"(0)").c_str() );
		      globalObs.add( *globalObsVar );

		      nuisVar = ws_->var( uncs[iunc].c_str() );
		      nuisVar->SetTitle( getNuisanceTitle(uncs[iunc]).c_str() ); //this is just for fancy labelling
		      nuis.add( *nuisVar );
		    }
		    
		  if(modVar==0)
		    {
		      float sigmaVal=descript[uncs[iunc].c_str()].toDouble();
		      sprintf(expBuf,"%s_sigma[%f,%f,%f]",uncs[iunc].c_str(),sigmaVal,sigmaVal,sigmaVal);
		      modVar=(RooRealVar *)ws_->factory(expBuf);
		      modVar->setConstant(true);
		    }

		  //add to formula
		  sprintf(expBuf,"*max((1+@%d*@%d),0.)",varCntr+1,varCntr+2);
		  formula+=expBuf;
		  varCntr+=2;
		  varsInFormula.add( *modVar );
		  varsInFormula.add( *nuisVar );
		}

	      //final formula
	      RooFormulaVar formulaVar(param.c_str(),formula.c_str(),varsInFormula);
	      ws_->import(formulaVar);
	    }

	  else
	    {
	      //if parameter has already been parsed, update the values only
	      cenVal->setVal(descript["val"].toDouble());
	      for(size_t iunc=0; iunc<uncs.size(); iunc++)
		{
		  if(uncs[iunc]=="val") continue;
		  
		  //add new constraint if required
		  RooRealVar *nuisVar=ws_->var((uncs[iunc]+"_sigma").c_str());
		  if(nuisVar==0) cout << "[Warning] failed to update uncertainty for: " << uncs[iunc] << " ...neglecting this source" << endl; 
		  else    nuisVar->setVal(descript[ uncs[iunc].c_str() ].toDouble());
		}
	    }
	}
    }

  if(ws_->set("nuisances")!=0) nuis.add( *(ws_->set("nuisances")) );
  ws_->defineSet("nuisances",nuis);

  if(ws_->set("constr")!=0) constr.add( *(ws_->set("constr")) );
  ws_->defineSet("constr",constr);
  
  if(ws_->set("globalobservables")!=0) globalObs.add( *ws_->set("globalobservables"));
  ws_->defineSet("globalobservables",globalObs);
}

//
void HFCMeasurement::initHFCModel()
{
  if(ws_==0) return;
  if(ws_->var("bmultobs")!=0) return;

  char expBuf[200];

  //observables
  ws_->factory("bmultobs[0.,99999999999.]");
  sprintf(expBuf,"bmult[0,%f]",float(maxJets_+1));
  ws_->factory(expBuf);
  ws_->var("bmult")->setBins(maxJets_+1);
  RooCategory sample("sample","sample");

  RooFormulaVar *r = (RooFormulaVar *)ws_->function("r");
  RooRealVar *sfeb  = ws_->var("sfeb");


  //
  // model
  //
  //acceptance corrections
  RooRealVar *acc1  = (RooRealVar *) ws_->factory("a1[1.0]");
  RooRealVar *acc05 = (RooRealVar *) ws_->factory("a05[1.0]");
  RooRealVar *acc0  = (RooRealVar *) ws_->factory("a0[1.0]");

  std::map<string,string> basePdfsPerCat; //auxiliary
  for(std::set<string>::iterator cIt=sampleCats_.begin(); cIt!=sampleCats_.end(); cIt++)
    {
      string tag=*cIt;
      int njets(2);
      if(tag.find("3jets")!=string::npos) njets=3;
      if(tag.find("4jets")!=string::npos) njets=4;

      //sample composition base parameters
      RooAbsReal *fcor = ws_->function(("fcor_"+tag).c_str());
      RooAbsReal *ftt  = ws_->function(("ftt_"+tag).c_str());
      RooAbsReal *kst  = ws_->function(("kst_"+tag).c_str());
      if(fcor==0 || ftt==0 || kst==0)
	{
	  cout << "[Warning] unable to instantiate model for " << tag << " check the configuration file for sample composition parameters" << endl;
	  continue;
	}

      //probability to select/reconstruct N(t->Wq)
      sprintf(expBuf,"FormulaVar::alpha_%s('@1*(2.+@2)>0 ? min((2*%d*@0)/(@1*(2.+@2)),1.0) : 1.0',{%s,%s,%s})",tag.c_str(),njets,fcor->GetName(),ftt->GetName(),kst->GetName());
      RooAbsReal *alpha = (RooAbsReal *) ws_->factory(expBuf);

      //sample composition in terms of N(t->Wq)
      sprintf(expBuf,"FormulaVar::alpha2_%s('pow(@0,2)*@1',{%s,%s})",tag.c_str(),alpha->GetName(),ftt->GetName());
      RooAbsReal *alpha2=(RooAbsReal *)ws_->factory(expBuf);
      sprintf(expBuf,"FormulaVar::alpha1_%s('(2*@0*(1-@0)+@0*@2)*@1',{%s,%s,%s})",tag.c_str(),alpha->GetName(),ftt->GetName(),kst->GetName());
      RooAbsReal *alpha1=(RooAbsReal *)ws_->factory(expBuf);
      sprintf(expBuf,"FormulaVar::alpha0_%s('1-@0-@1',{%s,%s})",tag.c_str(),alpha2->GetName(),alpha1->GetName());
      RooAbsReal *alpha0=(RooAbsReal *)ws_->factory(expBuf);


      //tag efficiencies for each sub-category in terms of number of jets from t->Wq reconstructed and selected
      std::vector<RooAbsReal *> eb,eq,eqstar;
      for(int ia=0; ia<=2; ia++)
	{
	  char pfbuf[10];
	  sprintf(pfbuf,"_%d",ia);

	  //b-tagging base parameters
	  RooAbsReal *abseb     = ws_->function(("abseb_"+tag+pfbuf).c_str());
	  RooAbsReal *abseq     = ws_->function(("abseq_"+tag+pfbuf).c_str());
	  if(abseb!=0 && abseq!=0)
	    {
	      //modify locally the b-tag rate with a correction factors (for SF_b fitting)
	      sprintf(expBuf,"FormulaVar::eb_%s_%d('max(min(@0*@1,1.0),0.0)',{%s,%s})",tag.c_str(),ia,abseb->GetName(),sfeb->GetName());
	      eb.push_back( (RooAbsReal *)ws_->factory(expBuf) );
	      eq.push_back( abseq );
	    }
	  else
	    {
	      if(ia!=0) cout << "[Warning] could not find eb or eq for " << tag << " with " << ia << " jets from top decay selected" << endl;
	      eb.push_back(0);
	      eq.push_back(0);
	    }

	  //extra jet contamination (non-matched jets)
	  RooAbsReal *abseqstar = ws_->function(("abseqstar_"+tag+pfbuf).c_str());
	  if(abseqstar!=0)
	    {
	      eqstar.push_back(abseqstar);
	    }
	  else
	    {
	      if(ia!=2 && njets!=2) cout << "[Warning] could not find eqstar for " << tag << " with " << ia << " jets from top decay selected" << endl;
	      eqstar.push_back(0);
	    }
	}


      //KERNEL PDFs
      //new version
     
      //alpha_0         @0 
      RooArgList a0args(*(eqstar[0]));
      RooFormulaVar pdf0a2j0t(("pdf0a2j0t_"+tag).c_str(), "pow(1-@0,2)", a0args);
      RooFormulaVar pdf0a2j1t(("pdf0a2j1t_"+tag).c_str(), "2*(1-@0)*@0", a0args);	
      RooFormulaVar pdf0a2j2t(("pdf0a2j2t_"+tag).c_str(), "pow(@0,2)",   a0args);

      //alpha_1:        @0 @1       @2       @3           @4    @5     @6             
      RooArgList a1args(*r,*(eb[1]),*(eq[1]),*(eqstar[1]),*acc1,*acc05,*acc0);
      RooFormulaVar pdf1a2j0t(("pdf1a2j0t_"+tag).c_str(), "pow(@0,2)*(1-@1)*(1-@3)*@4         + @0*(1-@0)*((1-@1)+(1-@2))*(1-@3)*@5        + pow(1-@0,2)*(1-@2)*(1-@3)*@6",         a1args); 
      RooFormulaVar pdf1a2j1t(("pdf1a2j1t_"+tag).c_str(), "pow(@0,2)*(@1*(1-@3)+(1-@1)*@3)*@4 + @0*(1-@0)*((@1+@2)*(1-@3)+(2-@1-@2)*@3)*@5 + pow(1-@0,2)*(@2*(1-@3)+(1-@2)*@3)*@6", a1args);
      RooFormulaVar pdf1a2j2t(("pdf1a2j2t_"+tag).c_str(), "pow(@0,2)*@1*@3*@4                 + @0*(1-@0)*(@1+@2)*@3*@5                    + pow(1-@0,2)*@2*@3*@6",                 a1args);

      //alpha_2 :       @0 @1       @2       @3    @4     @5
      RooArgList a2args(*r,*(eb[2]),*(eq[2]),*acc1,*acc05,*acc0);
      RooFormulaVar pdf2a2j0t(("pdf2a2j0t_"+tag).c_str(), "pow(@0*(1-@1),2)*@3      + 2*@0*(1-@0)*(1-@1)*(1-@2)*@4         + pow((1-@0)*(1-@2),2)*@5",      a2args);
      RooFormulaVar pdf2a2j1t(("pdf2a2j1t_"+tag).c_str(), "2*pow(@0,2)*(1-@1)*@1*@3 + 2*@0*(1-@0)*((1-@1)*@2+@1*(1-@2))*@4 + 2*pow((1-@0),2)*(1-@2)*@2*@5", a2args);
      RooFormulaVar pdf2a2j2t(("pdf2a2j2t_"+tag).c_str(), "pow(@0*@1,2)*@3          + 2*@0*(1-@0)*@1*@2*@4                 + pow((1-@0)*@2,2)*@5",          a2args);

      //define the specialized PDFs per jet multiplicity
      if(njets==2)
	{
	  //add according to alpha coefficients
	  RooGenericPdf pdf2j0t(("pdf2j0t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a2j0t,pdf1a2j0t,pdf0a2j0t,*alpha2,*alpha1,*alpha0));
	  RooGenericPdf pdf2j1t(("pdf2j1t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a2j1t,pdf1a2j1t,pdf0a2j1t,*alpha2,*alpha1,*alpha0));
	  RooGenericPdf pdf2j2t(("pdf2j2t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a2j2t,pdf1a2j2t,pdf0a2j2t,*alpha2,*alpha1,*alpha0));

	  //import to workspace
	  ws_->import( pdf2j0t );
	  ws_->import( pdf2j1t );
	  ws_->import( pdf2j2t );
	}
      else
	{
	  //extra terms for "non-matched" tags 
	  RooFormulaVar pdf0a0t(("pdf0a0t_"+tag).c_str(), "(1-@0)", RooArgList(*(eqstar[0])) );
	  RooFormulaVar pdf0a1t(("pdf0a1t_"+tag).c_str(), "@0",     RooArgList(*(eqstar[0])) );
	  RooFormulaVar pdf1a0t(("pdf1a0t_"+tag).c_str(), "(1-@0)", RooArgList(*(eqstar[1])) );
	  RooFormulaVar pdf1a1t(("pdf1a1t_"+tag).c_str(), "@0",     RooArgList(*(eqstar[1])) );
	  RooFormulaVar pdf2a0t(("pdf2a0t_"+tag).c_str(), "(1-@0)", RooArgList(*(eqstar[2])) );
	  RooFormulaVar pdf2a1t(("pdf2a1t_"+tag).c_str(), "@0",     RooArgList(*(eqstar[2])) );

	  if(njets==3)
	    {
	      RooFormulaVar pdf0a3j0t(("pdf0a3j0t_"+tag).c_str(),"@0*@1",       RooArgList(pdf0a2j0t,pdf0a0t));
	      RooFormulaVar pdf0a3j1t(("pdf0a3j1t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf0a2j0t,pdf0a1t,pdf0a2j1t,pdf0a0t));
	      RooFormulaVar pdf0a3j2t(("pdf0a3j2t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf0a2j1t,pdf0a1t,pdf0a2j2t,pdf0a0t));
	      RooFormulaVar pdf0a3j3t(("pdf0a3j3t_"+tag).c_str(),"@0*@1",       RooArgList(pdf0a2j2t,pdf0a1t));

	      RooFormulaVar pdf1a3j0t(("pdf1a3j0t_"+tag).c_str(),"@0*@1",       RooArgList(pdf1a2j0t,pdf1a0t));
	      RooFormulaVar pdf1a3j1t(("pdf1a3j1t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf1a2j0t,pdf1a1t,pdf1a2j1t,pdf1a0t));
	      RooFormulaVar pdf1a3j2t(("pdf1a3j2t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf1a2j1t,pdf1a1t,pdf1a2j2t,pdf1a0t));
	      RooFormulaVar pdf1a3j3t(("pdf1a3j3t_"+tag).c_str(),"@0*@1",       RooArgList(pdf1a2j2t,pdf1a1t));

	      RooFormulaVar pdf2a3j0t(("pdf2a3j0t_"+tag).c_str(),"@0*@1",       RooArgList(pdf2a2j0t,pdf2a0t));
	      RooFormulaVar pdf2a3j1t(("pdf2a3j1t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf2a2j0t,pdf2a1t,pdf2a2j1t,pdf2a0t));
	      RooFormulaVar pdf2a3j2t(("pdf2a3j2t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf2a2j1t,pdf2a1t,pdf2a2j2t,pdf2a0t));
	      RooFormulaVar pdf2a3j3t(("pdf2a3j3t_"+tag).c_str(),"@0*@1",       RooArgList(pdf2a2j2t,pdf2a1t));

	      //add according to alpha coefficients
	      RooGenericPdf pdf3j0t(("pdf3j0t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a3j0t,pdf1a3j0t,pdf0a3j0t,*alpha2,*alpha1,*alpha0));
	      RooGenericPdf pdf3j1t(("pdf3j1t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a3j1t,pdf1a3j1t,pdf0a3j1t,*alpha2,*alpha1,*alpha0));
	      RooGenericPdf pdf3j2t(("pdf3j2t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a3j2t,pdf1a3j2t,pdf0a3j2t,*alpha2,*alpha1,*alpha0));
	      RooGenericPdf pdf3j3t(("pdf3j3t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a3j3t,pdf1a3j3t,pdf0a3j3t,*alpha2,*alpha1,*alpha0));

	      //import to workspace
	      ws_->import( pdf3j0t );
	      ws_->import( pdf3j1t,RecycleConflictNodes() );
	      ws_->import( pdf3j2t,RecycleConflictNodes() );
	      ws_->import( pdf3j3t,RecycleConflictNodes() );
	    }
	  if(njets==4)
	    {
	      RooFormulaVar pdf0a4j0t(("pdf0a4j0t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf0a2j0t,pdf0a0t));
	      RooFormulaVar pdf0a4j1t(("pdf0a4j1t_"+tag).c_str(),"2*@0*@2*@3+@1*@3*@3",            RooArgList(pdf0a2j0t,pdf0a2j1t,pdf0a1t,pdf0a0t));
	      RooFormulaVar pdf0a4j2t(("pdf0a4j2t_"+tag).c_str(),"@0*@3*@3+2*@1*@3*@4+@2*@4*@4",   RooArgList(pdf0a2j0t,pdf0a2j1t,pdf0a2j2t,pdf0a1t,pdf0a0t));
	      RooFormulaVar pdf0a4j3t(("pdf0a4j3t_"+tag).c_str(),"@0*@2*@2+2*@1*@2*@3",            RooArgList(pdf0a2j1t,pdf0a2j2t,pdf0a1t,pdf0a0t));
	      RooFormulaVar pdf0a4j4t(("pdf0a4j4t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf0a2j2t,pdf0a1t));

	      RooFormulaVar pdf1a4j0t(("pdf1a4j0t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf1a2j0t,pdf1a0t));
	      RooFormulaVar pdf1a4j1t(("pdf1a4j1t_"+tag).c_str(),"2*@0*@2*@3+@1*@3*@3",            RooArgList(pdf1a2j0t,pdf1a2j1t,pdf1a1t,pdf1a0t));
	      RooFormulaVar pdf1a4j2t(("pdf1a4j2t_"+tag).c_str(),"@0*@3*@3+2*@1*@3*@4+@2*@4*@4",   RooArgList(pdf1a2j0t,pdf1a2j1t,pdf1a2j2t,pdf1a1t,pdf1a0t));
	      RooFormulaVar pdf1a4j3t(("pdf1a4j3t_"+tag).c_str(),"@0*@2*@2+2*@1*@2*@3",            RooArgList(pdf1a2j1t,pdf1a2j2t,pdf1a1t,pdf1a0t));
	      RooFormulaVar pdf1a4j4t(("pdf1a4j4t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf1a2j2t,pdf1a1t));

	      RooFormulaVar pdf2a4j0t(("pdf2a4j0t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf2a2j0t,pdf2a0t));
	      RooFormulaVar pdf2a4j1t(("pdf2a4j1t_"+tag).c_str(),"2*@0*@2*@3+@1*@3*@3",            RooArgList(pdf2a2j0t,pdf2a2j1t,pdf2a1t,pdf2a0t));
	      RooFormulaVar pdf2a4j2t(("pdf2a4j2t_"+tag).c_str(),"@0*@3*@3+2*@1*@3*@4+@2*@4*@4",   RooArgList(pdf2a2j0t,pdf2a2j1t,pdf2a2j2t,pdf2a1t,pdf2a0t));
	      RooFormulaVar pdf2a4j3t(("pdf2a4j3t_"+tag).c_str(),"@0*@2*@2+2*@1*@2*@3",            RooArgList(pdf2a2j1t,pdf2a2j2t,pdf2a1t,pdf2a0t));
	      RooFormulaVar pdf2a4j4t(("pdf2a4j4t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf2a2j2t,pdf2a1t));

	      //add according to alpha coefficients
	      RooGenericPdf pdf4j0t(("pdf4j0t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a4j0t,pdf1a4j0t,pdf0a4j0t,*alpha2,*alpha1,*alpha0));
	      RooGenericPdf pdf4j1t(("pdf4j1t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a4j1t,pdf1a4j1t,pdf0a4j1t,*alpha2,*alpha1,*alpha0));
	      RooGenericPdf pdf4j2t(("pdf4j2t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a4j2t,pdf1a4j2t,pdf0a4j2t,*alpha2,*alpha1,*alpha0));
	      RooGenericPdf pdf4j3t(("pdf4j3t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a4j3t,pdf1a4j3t,pdf0a4j3t,*alpha2,*alpha1,*alpha0));
	      RooGenericPdf pdf4j4t(("pdf4j4t_"+tag).c_str(),"@0*@3+@1*@4+@2*@5",RooArgList(pdf2a4j4t,pdf1a4j4t,pdf0a4j4t,*alpha2,*alpha1,*alpha0));

	      //import to workspace
	      ws_->import( pdf4j0t );
	      ws_->import( pdf4j1t, RecycleConflictNodes() );
	      ws_->import( pdf4j2t, RecycleConflictNodes() );
	      ws_->import( pdf4j3t, RecycleConflictNodes() );
	      ws_->import( pdf4j4t, RecycleConflictNodes() );
	    }
	}

      


      /* pre-approval version

      //alpha_0         @0 
      RooArgList a0args(*(eqstar[0]));
      RooGenericPdf pdf0a2j0t(("pdf0a2j0t_"+tag).c_str(), "pow(1-@0,2)", a0args);
      RooGenericPdf pdf0a2j1t(("pdf0a2j1t_"+tag).c_str(), "2*(1-@0)*@0", a0args);	
      RooGenericPdf pdf0a2j2t(("pdf0a2j2t_"+tag).c_str(), "pow(@0,2)",   a0args);

      //alpha_1:        @0 @1       @2       @3           @4    @5     @6             
      RooArgList a1args(*r,*(eb[1]),*(eq[1]),*(eqstar[1]),*acc1,*acc05,*acc0);
      RooGenericPdf pdf1a2j0t(("pdf1a2j0t_"+tag).c_str(), "pow(@0,2)*(1-@1)*(1-@3)*@4         + @0*(1-@0)*((1-@1)+(1-@2))*(1-@3)*@5        + pow(1-@0,2)*(1-@2)*(1-@3)*@6",         a1args); 
      RooGenericPdf pdf1a2j1t(("pdf1a2j1t_"+tag).c_str(), "pow(@0,2)*(@1*(1-@3)+(1-@1)*@3)*@4 + @0*(1-@0)*((@1+@2)*(1-@3)+(2-@1-@2)*@3)*@5 + pow(1-@0,2)*(@2*(1-@3)+(1-@2)*@3)*@6", a1args);
      RooGenericPdf pdf1a2j2t(("pdf1a2j2t_"+tag).c_str(), "pow(@0,2)*@1*@3*@4                 + @0*(1-@0)*(@1+@2)*@3*@5                    + pow(1-@0,2)*@2*@3*@6",                 a1args);

      //alpha_2 :       @0 @1       @2       @3    @4     @5
      RooArgList a2args(*r,*(eb[2]),*(eq[2]),*acc1,*acc05,*acc0);
      RooGenericPdf pdf2a2j0t(("pdf2a2j0t_"+tag).c_str(), "pow(@0*(1-@1),2)*@3      + 2*@0*(1-@0)*(1-@1)*(1-@2)*@4         + pow((1-@0)*(1-@2),2)*@5",      a2args);
      RooGenericPdf pdf2a2j1t(("pdf2a2j1t_"+tag).c_str(), "2*pow(@0,2)*(1-@1)*@1*@3 + 2*@0*(1-@0)*((1-@1)*@2+@1*(1-@2))*@4 + 2*pow((1-@0),2)*(1-@2)*@2*@5", a2args);
      RooGenericPdf pdf2a2j2t(("pdf2a2j2t_"+tag).c_str(), "pow(@0*@1,2)*@3          + 2*@0*(1-@0)*@1*@2*@4                 + pow((1-@0)*@2,2)*@5",          a2args);

      //define the specialized PDFs per jet multiplicity
      RooArgList acoeffs(*alpha2,*alpha1,*alpha0);
      if(njets==2)
	{
	  //add according to alpha coefficients
	  RooAddPdf pdf2j0t(("pdf2j0t_"+tag).c_str(),"",RooArgList(pdf2a2j0t,pdf1a2j0t,pdf0a2j0t),acoeffs);
	  RooAddPdf pdf2j1t(("pdf2j1t_"+tag).c_str(),"",RooArgList(pdf2a2j1t,pdf1a2j1t,pdf0a2j1t),acoeffs);
	  RooAddPdf pdf2j2t(("pdf2j2t_"+tag).c_str(),"",RooArgList(pdf2a2j2t,pdf1a2j2t,pdf0a2j2t),acoeffs);

	  //import to workspace
	  ws_->import( pdf2j0t );
	  ws_->import( pdf2j1t );
	  ws_->import( pdf2j2t );
	}
      else
	{
	  //extra terms for "non-matched" tags 
	  RooGenericPdf pdf0a0t(("pdf0a0t_"+tag).c_str(), "(1-@0)", RooArgList(*(eqstar[0])) );
	  RooGenericPdf pdf0a1t(("pdf0a1t_"+tag).c_str(), "@0",     RooArgList(*(eqstar[0])) );
	  RooGenericPdf pdf1a0t(("pdf1a0t_"+tag).c_str(), "(1-@0)", RooArgList(*(eqstar[1])) );
	  RooGenericPdf pdf1a1t(("pdf1a1t_"+tag).c_str(), "@0",     RooArgList(*(eqstar[1])) );
	  RooGenericPdf pdf2a0t(("pdf2a0t_"+tag).c_str(), "(1-@0)", RooArgList(*(eqstar[2])) );
	  RooGenericPdf pdf2a1t(("pdf2a1t_"+tag).c_str(), "@0",     RooArgList(*(eqstar[2])) );

	  if(njets==3)
	    {
	      RooGenericPdf pdf0a3j0t(("pdf0a3j0t_"+tag).c_str(),"@0*@1",       RooArgList(pdf0a2j0t,pdf0a0t));
	      RooGenericPdf pdf0a3j1t(("pdf0a3j1t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf0a2j0t,pdf0a1t,pdf0a2j1t,pdf0a0t));
	      RooGenericPdf pdf0a3j2t(("pdf0a3j2t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf0a2j1t,pdf0a1t,pdf0a2j2t,pdf0a0t));
	      RooGenericPdf pdf0a3j3t(("pdf0a3j3t_"+tag).c_str(),"@0*@1",       RooArgList(pdf0a2j2t,pdf0a1t));

	      RooGenericPdf pdf1a3j0t(("pdf1a3j0t_"+tag).c_str(),"@0*@1",       RooArgList(pdf1a2j0t,pdf1a0t));
	      RooGenericPdf pdf1a3j1t(("pdf1a3j1t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf1a2j0t,pdf1a1t,pdf1a2j1t,pdf1a0t));
	      RooGenericPdf pdf1a3j2t(("pdf1a3j2t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf1a2j1t,pdf1a1t,pdf1a2j2t,pdf1a0t));
	      RooGenericPdf pdf1a3j3t(("pdf1a3j3t_"+tag).c_str(),"@0*@1",       RooArgList(pdf1a2j2t,pdf1a1t));

	      RooGenericPdf pdf2a3j0t(("pdf2a3j0t_"+tag).c_str(),"@0*@1",       RooArgList(pdf2a2j0t,pdf2a0t));
	      RooGenericPdf pdf2a3j1t(("pdf2a3j1t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf2a2j0t,pdf2a1t,pdf2a2j1t,pdf2a0t));
	      RooGenericPdf pdf2a3j2t(("pdf2a3j2t_"+tag).c_str(),"@0*@1+@2*@3", RooArgList(pdf2a2j1t,pdf2a1t,pdf2a2j2t,pdf2a0t));
	      RooGenericPdf pdf2a3j3t(("pdf2a3j3t_"+tag).c_str(),"@0*@1",       RooArgList(pdf2a2j2t,pdf2a1t));

	      //add according to alpha coefficients
	      RooAddPdf pdf3j0t(("pdf3j0t_"+tag).c_str(),"",RooArgList(pdf2a3j0t,pdf1a3j0t,pdf0a3j0t),acoeffs);
	      RooAddPdf pdf3j1t(("pdf3j1t_"+tag).c_str(),"",RooArgList(pdf2a3j1t,pdf1a3j1t,pdf0a3j1t),acoeffs);
	      RooAddPdf pdf3j2t(("pdf3j2t_"+tag).c_str(),"",RooArgList(pdf2a3j2t,pdf1a3j2t,pdf0a3j2t),acoeffs);
	      RooAddPdf pdf3j3t(("pdf3j3t_"+tag).c_str(),"",RooArgList(pdf2a3j3t,pdf1a3j3t,pdf0a3j3t),acoeffs);

	      //import to workspace
	      ws_->import( pdf3j0t );
	      ws_->import( pdf3j1t,RecycleConflictNodes() );
	      ws_->import( pdf3j2t,RecycleConflictNodes() );
	      ws_->import( pdf3j3t,RecycleConflictNodes() );
	    }
	  if(njets==4)
	    {
	      RooGenericPdf pdf0a4j0t(("pdf0a4j0t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf0a2j0t,pdf0a0t));
	      RooGenericPdf pdf0a4j1t(("pdf0a4j1t_"+tag).c_str(),"2*@0*@2*@3+@1*@3*@3",            RooArgList(pdf0a2j0t,pdf0a2j1t,pdf0a1t,pdf0a0t));
	      RooGenericPdf pdf0a4j2t(("pdf0a4j2t_"+tag).c_str(),"@0*@3*@3+2*@1*@3*@4+@2*@4*@4",   RooArgList(pdf0a2j0t,pdf0a2j1t,pdf0a2j2t,pdf0a1t,pdf0a0t));
	      RooGenericPdf pdf0a4j3t(("pdf0a4j3t_"+tag).c_str(),"@0*@2*@2+2*@1*@2*@3",            RooArgList(pdf0a2j1t,pdf0a2j2t,pdf0a1t,pdf0a0t));
	      RooGenericPdf pdf0a4j4t(("pdf0a4j4t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf0a2j2t,pdf0a1t));

	      RooGenericPdf pdf1a4j0t(("pdf1a4j0t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf1a2j0t,pdf1a0t));
	      RooGenericPdf pdf1a4j1t(("pdf1a4j1t_"+tag).c_str(),"2*@0*@2*@3+@1*@3*@3",            RooArgList(pdf1a2j0t,pdf1a2j1t,pdf1a1t,pdf1a0t));
	      RooGenericPdf pdf1a4j2t(("pdf1a4j2t_"+tag).c_str(),"@0*@3*@3+2*@1*@3*@4+@2*@4*@4",   RooArgList(pdf1a2j0t,pdf1a2j1t,pdf1a2j2t,pdf1a1t,pdf1a0t));
	      RooGenericPdf pdf1a4j3t(("pdf1a4j3t_"+tag).c_str(),"@0*@2*@2+2*@1*@2*@3",            RooArgList(pdf1a2j1t,pdf1a2j2t,pdf1a1t,pdf1a0t));
	      RooGenericPdf pdf1a4j4t(("pdf1a4j4t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf1a2j2t,pdf1a1t));

	      RooGenericPdf pdf2a4j0t(("pdf2a4j0t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf2a2j0t,pdf2a0t));
	      RooGenericPdf pdf2a4j1t(("pdf2a4j1t_"+tag).c_str(),"2*@0*@2*@3+@1*@3*@3",            RooArgList(pdf2a2j0t,pdf2a2j1t,pdf2a1t,pdf2a0t));
	      RooGenericPdf pdf2a4j2t(("pdf2a4j2t_"+tag).c_str(),"@0*@3*@3+2*@1*@3*@4+@2*@4*@4",   RooArgList(pdf2a2j0t,pdf2a2j1t,pdf2a2j2t,pdf2a1t,pdf2a0t));
	      RooGenericPdf pdf2a4j3t(("pdf2a4j3t_"+tag).c_str(),"@0*@2*@2+2*@1*@2*@3",            RooArgList(pdf2a2j1t,pdf2a2j2t,pdf2a1t,pdf2a0t));
	      RooGenericPdf pdf2a4j4t(("pdf2a4j4t_"+tag).c_str(),"@0*@1*@1",                       RooArgList(pdf2a2j2t,pdf2a1t));

	      //add according to alpha coefficients
	      RooAddPdf pdf4j0t(("pdf4j0t_"+tag).c_str(),"",RooArgList(pdf2a4j0t,pdf1a4j0t,pdf0a4j0t),acoeffs);
	      RooAddPdf pdf4j1t(("pdf4j1t_"+tag).c_str(),"",RooArgList(pdf2a4j1t,pdf1a4j1t,pdf0a4j1t),acoeffs);
	      RooAddPdf pdf4j2t(("pdf4j2t_"+tag).c_str(),"",RooArgList(pdf2a4j2t,pdf1a4j2t,pdf0a4j2t),acoeffs);
	      RooAddPdf pdf4j3t(("pdf4j3t_"+tag).c_str(),"",RooArgList(pdf2a4j3t,pdf1a4j3t,pdf0a4j3t),acoeffs);
	      RooAddPdf pdf4j4t(("pdf4j4t_"+tag).c_str(),"",RooArgList(pdf2a4j4t,pdf1a4j4t,pdf0a4j4t),acoeffs);

	      //import to workspace
	      ws_->import( pdf4j0t );
	      ws_->import( pdf4j1t, RecycleConflictNodes() );
	      ws_->import( pdf4j2t, RecycleConflictNodes() );
	      ws_->import( pdf4j3t, RecycleConflictNodes() );
	      ws_->import( pdf4j4t, RecycleConflictNodes() );
	    }
	}

      */

      //each n-jets n-b-tags is a category
      for(int itag=0; itag<=njets; itag++)
      	{
	  char catname[100];
	  sprintf(catname,"n%dbtags_%s",itag,tag.c_str());
	  sample.defineType(catname);

	  char pdfname[100];	  
	  sprintf(pdfname,"pdf%dj%dt_%s",njets,itag,tag.c_str());
	  basePdfsPerCat[catname]=pdfname;
	}
    }

  //base model is a simultaneous pdf for each category
  RooSimultaneous *basemodel = new RooSimultaneous("basemodel","basemodel",sample);
  for(std::map<string, string>::iterator it=basePdfsPerCat.begin(); it!=basePdfsPerCat.end(); it++)
    {
      RooAbsPdf *pdf=ws_->pdf(it->second.c_str());
      if(pdf==0)
	{
	  cout << "[Warning] failed to link " << it->second << " to category " << it->first << endl;
	  continue;
	}
      basemodel->addPdf( *pdf, it->first.c_str() );
    }

  //add the event yields so it can be extended
  RooRealVar *yields = (RooRealVar *)ws_->factory("yields[1,0,999999999.]");

  //add the constraints if any
  RooAbsPdf *model=0;
  if(ws_->set("constr")!=0)
    {
      RooExtendPdf *extmodel = new RooExtendPdf("extmodel","extmodel",*basemodel,*yields);
      RooArgSet modelFactors(*extmodel); 
      modelFactors.add( *(ws_->set("constr")) );
      model = (RooAbsPdf *)(new RooProdPdf("model","model",modelFactors));
    }
  else
    {
      model = (RooAbsPdf *)(new RooExtendPdf("model","model",*basemodel,*yields));
    }

  //import the model
  ws_->import(*model);

  //finalize workspace
  ws_->defineSet("observables","bmult,bmultobs,sample");
  RooArgSet baseGlobalObs( *(ws_->var("bmult")), *(ws_->var("bmultobs")) ); 
  baseGlobalObs.add( *(ws_->set("globalobservables") ) );
  ws_->defineSet("globalobservables", baseGlobalObs );
  RooArgSet *nullParams = (RooArgSet *)ws_->allVars().snapshot();
  ws_->saveSnapshot("default",*nullParams,kTRUE);

  printConfiguration(cout);

  //instantiate a model configurator for RooStats
  mc_ = new ModelConfig("mc",ws_);
  mc_->SetPdf(*model);
  mc_->SetParametersOfInterest(*(ws_->set("poi")));
  mc_->SetObservables(*(ws_->set("observables")));
  mc_->SetGlobalObservables(*(ws_->set("globalobservables")));
  mc_->SetNuisanceParameters(*(ws_->set("nuisances")));
}


//
void HFCMeasurement::resetModelValues()
{
  if(ws_==0) return;
  ws_->loadSnapshot("default");
}

//
void HFCMeasurement::fitHFCfrom(TH1 *btagObs, bool debug)
{
  if( !createDataset(btagObs) ) return;;
  if(debug) { 
    cout << "Saving workspace to HeavyFlavorWorkspace.root" << endl; 
    //     TIterator *nuis_params_itr = mc_->GetNuisanceParameters()->createIterator();
    //     TObject *nuis_params_obj;
    //     while((nuis_params_obj=nuis_params_itr->Next())){
    //       RooRealVar *nuiVar=(RooRealVar *)nuis_params_obj;
    //       nuiVar->setRange(-1,1); //limit to +/-1-sigma for FC
    //     }
    saveWorkspace("HeavyFlavorWorkspace.root");
  }
  fit(debug);
}

//
HFCMeasurement::FitResult_t HFCMeasurement::plrFit(TH1F *h)
{
  createDataset(h);
  return plrFit(data_,mc_,false);
}




//
TH1F *HFCMeasurement::generateBtagObs()
{
  TH1F *btagObs = new TH1F("btagsextended","btagsextended",3*3*5,0.,3*3*5);
  btagObs->SetDirectory(0);

  for(int ich=0; ich<=2; ich++)
    {
      TString ch("ee");
      if(ich==1) ch="mumu";
      if(ich==2) ch="emu";
      
      for(int njets=2; njets<=maxJets_; njets++)
	{
	  int ibin(0);
	  if(njets==3) ibin +=5;
	  if(njets==4) ibin +=10;
	  if(ich==1)   ibin +=15;
	  if(ich==2)   ibin += 2*15;
	  TString tag(ch); tag+=njets; tag+="jets";

	  //get the number of observed events after projecting the data
	  //TString cut("sample==sample::n1btags_"+tag+" || sample==sample::n2btags_"+tag);
	  TString cut("sample==sample::n0btags_"+tag+" || sample==sample::n1btags_"+tag+" || sample==sample::n2btags_"+tag);
	  if(njets>2) cut+=" || sample==sample::n3btags_"+tag;
	  if(njets>3) cut+=" || sample==sample::n4btags_"+tag;
	  
	  RooDataSet *data = (RooDataSet *) data_->reduce(cut);
	  Double_t cts=data->sumEntries();
	  
	  //project the model with these counts
	  TH1F *icatH=getProjectedModel(tag,cts,"tmp","tmp");
	  for(int jbin=1; jbin<=icatH->GetXaxis()->GetNbins(); jbin++)
	    {
	      int idx=ibin+jbin;
	      btagObs->SetBinContent(idx,icatH->GetBinContent(jbin));
	      btagObs->SetBinError(idx,sqrt(icatH->GetBinContent(jbin)));
	    }
	  delete icatH;
	}
    }

  //all done
  return btagObs;
}


//
bool HFCMeasurement::createDataset(TH1 *btagObs)
{
  if(btagObs==0) return false;

  //re-create a RooDataSet with categories
  RooRealVar *bmult    = ws_->var("bmult");
  RooRealVar *bmultObs = ws_->var("bmultobs");
  RooCategory *sample  = ws_->cat("sample");
  if(data_==0) data_ = new RooDataSet("data","data",RooArgSet( *bmult, *bmultObs, *sample), RooFit::WeightVar(*bmultObs) );
  else         data_->reset();
  for(int ibin=1; ibin<btagObs->GetXaxis()->GetNbins(); ibin++)
    {
      float counts(btagObs->GetBinContent(ibin));
      if(counts==0) continue;

      //build the category name
      int chCode((ibin-1)/(3*5));
      int njets(((ibin-1)/5)%3+2);
      int ntags((ibin-1)%5);
      TString ch("ee");
      if(chCode==1) ch="mumu";
      if(chCode==2) ch="emu";
      char catName[200];
      sprintf(catName,"n%dbtags_%s%djets",ntags,ch.Data(),njets);

      //add entry
      bmult->setVal(ntags);
      bmultObs->setVal(counts);  
      sample->setLabel(catName);
      data_->add( RooArgSet(*bmult,*bmultObs,*sample), counts );
    }


  return true;
}

//
bool HFCMeasurement::fit(bool debug)
{
  bool fitResult(false);
  
  //check inputs
  RooAbsPdf *pdf   = ws_->pdf("model");
  if(data_==0 || mc_==0 || pdf==0) return fitResult;

  //reset results
  curResults_.clear();

  //exclusive categories
  std::map<string,RooDataSet *> dataSlice;
  std::map<string, TH1F *> preFitModel, postExclusiveFitModel, postInclusiveFitModel;
  std::map<string, std::vector<TH1F *> > preFitModelFunc;
  std::map<string, TString> combChannelCuts;
  std::map<int, TString> combMultCuts; 
  TString allCuts("");
  for(std::set<string>::iterator cIt = sampleCats_.begin(); cIt != sampleCats_.end(); cIt++)
    {
      TString tag=*cIt;
      int njets(2);
      if(tag.Contains("3")) njets=3;
      if(tag.Contains("4")) njets=4;
      
      //      if(tag.Contains("ee") && njets==4) continue;

      ////project data
      //TString cut("sample==sample::n1btags_"+tag+" || sample==sample::n2btags_"+tag);
      TString cut("sample==sample::n0btags_"+tag+" || sample==sample::n1btags_"+tag+" || sample==sample::n2btags_"+tag);
      if(njets>2) cut+=" || sample==sample::n3btags_"+tag;
      if(njets>3) cut+=" || sample==sample::n4btags_"+tag;
      RooDataSet *data = (RooDataSet *) data_->reduce(cut);
      dataSlice[tag.Data()]=data;
      Double_t cts=data->sumEntries();

      //add to the combined channel cuts
      string key("ee");
      if(tag.Contains("mumu")) key="mumu";
      if(tag.Contains("emu")) key="emu";
      if(combChannelCuts.find(key)==combChannelCuts.end()) combChannelCuts[key]=cut;
      else                                                 combChannelCuts[key]=combChannelCuts[key] + " || " + cut; 
      if(allCuts=="") allCuts=cut;
      else            allCuts=allCuts + " || " + cut;

      //add to the combined multiplicity cuts
      if(combMultCuts.find(njets)==combMultCuts.end()) combMultCuts[njets]=cut;
      else                                             combMultCuts[njets]=combMultCuts[njets]+ " || " + cut; 

      //reset to default values
      resetModelValues();

      //save the pre-fit snapshot
      preFitModel[tag.Data()]=getProjectedModel(tag,cts,tag+"_prefitprojpdf","Pre-fit");
      preFitModelFunc[tag.Data()]=getProbabilityModelFunction(tag,tag+"_prefitprojpdffunc");

      //profile likelihood method 
      curResults_[tag.Data()]=plrFit(data,mc_,debug);
      
      //save the pre-fit snapshot
      postExclusiveFitModel[tag.Data()]=getProjectedModel(tag,cts,tag+"_postexcfitprojpdf","Post-fit (exc.)");
    }

  //exclusive channel results
  for(std::map<std::string,TString>::iterator cIt=combChannelCuts.begin(); cIt!=combChannelCuts.end(); cIt++)
    {
      resetModelValues();
      RooDataSet *data = (RooDataSet *) data_->reduce(cIt->second);
      dataSlice[cIt->first]   = data;
      curResults_[cIt->first] = plrFit(data,mc_,debug);
    } 

  //exclusive multiplicity results
  for(std::map<int,TString>::iterator cIt=combMultCuts.begin(); cIt!=combMultCuts.end(); cIt++)
    {
      char tag[200]; sprintf(tag,"=%d jets",cIt->first);
      resetModelValues();
      RooDataSet *data = (RooDataSet *) data_->reduce(cIt->second);
      dataSlice[tag]   = data;
      curResults_[tag] = plrFit(data,mc_,debug);
    } 


  //the inclusive result <- this is the final result to be used
  resetModelValues();
  //curResults_["inclusive"]=plrFit(data_,mc_,debug);
  RooDataSet *data = (RooDataSet *) data_->reduce(allCuts);
  curResults_["inclusive"]=plrFit(data,mc_,debug,debug);
  for(std::set<string>::iterator cIt = sampleCats_.begin(); cIt != sampleCats_.end(); cIt++)
    {
      string tag=*cIt;
      postInclusiveFitModel[tag]=getProjectedModel(tag,postExclusiveFitModel[tag]->Integral(),tag+"_postincfitprojpdf","Post-fit (inc.)");
    }
  
  //show results of the fit if required
  if(!debug) return fitResult;

  RooRealVar *bmult    = ws_->var("bmult");
  TCanvas *c = new TCanvas("c","c",1200,1200);
  c->Divide(maxJets_-1,3);
  TCanvas *modelc = new TCanvas("modelc","modelc",1200,1200);
  modelc->Divide(maxJets_-1,3);
  TCanvas *pfnc = new TCanvas("pfnc","pfnc",1200,1200);
  pfnc->Divide(maxJets_-1,3);
  int icat=0;
  TH1F *dataExtended=new TH1F("dataextended",";b-tag multiplicity;Events", 3*3*5, 0.,3*3*5); dataExtended->SetDirectory(0);
  TH1F *modelExtended=(TH1F *)dataExtended->Clone("modelextended");                          modelExtended->SetDirectory(0);        modelExtended->SetLineColor(30);
  TH1F *postFitModelExtended=(TH1F *) dataExtended->Clone("pfmodelextended");                postFitModelExtended->SetDirectory(0); postFitModelExtended->SetLineColor(kBlue);
  for(std::set<string>::iterator cIt = sampleCats_.begin(); cIt != sampleCats_.end(); cIt++, icat++)
    {
      string tag=*cIt;

      //select only fully exclusive categories
      if(tag.find("jets")==string::npos) { icat--; continue; }

      c->cd(icat+1);
      RooPlot *frame = bmult->frame();
      RooDataSet *data = dataSlice[tag];
      data->plotOn(frame,Name("data"));
      frame->Draw();
      preFitModel[tag]->Draw("histsame");           preFitModel[tag]->SetLineStyle(7); preFitModel[tag]->SetLineColor(30);
      postExclusiveFitModel[tag]->Draw("histsame"); postExclusiveFitModel[tag]->SetLineColor(30);
      postInclusiveFitModel[tag]->Draw("histsame");

      for(int ibin=1; ibin<=modelExtended->GetXaxis()->GetNbins(); ibin++)
	{
	  Int_t ncounts(0);
	  if(ibin<=data->numEntries())
	    {
	      data->get(ibin-1);
	      ncounts = data->weight();
	    }
	  dataExtended->SetBinContent(ibin+icat*5,ncounts);
	  dataExtended->SetBinError(ibin+icat*5,sqrt(ncounts));
	  modelExtended->SetBinContent(ibin+icat*5,preFitModel[tag]->GetBinContent(ibin));
	  modelExtended->SetBinError(ibin+icat*5,preFitModel[tag]->GetBinError(ibin));
	  postFitModelExtended->SetBinContent(ibin+icat*5,postInclusiveFitModel[tag]->GetBinContent(ibin));
	  postFitModelExtended->SetBinError(ibin+icat*5,postInclusiveFitModel[tag]->GetBinError(ibin));
	}

      frame->GetXaxis()->SetNdivisions(maxJets_+1);
      frame->GetXaxis()->SetTitle("b-tag multiplicity");
      frame->GetYaxis()->SetTitle("Events");
      for(int ibin=1; ibin<=frame->GetXaxis()->GetNbins(); ibin++)
	{
	  TString label("="); label+=(ibin-1);
	  frame->GetXaxis()->SetBinLabel(ibin,label);
	}
      frame->GetYaxis()->SetRangeUser(0.01,data->sumEntries()*0.55);
      frame->GetYaxis()->SetTitleOffset(1.4);
	  
      //category
      TPaveText *pt = new TPaveText(0.6,0.85,0.85,0.94,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      TString caption=tag.c_str();
      caption=caption.ReplaceAll("mu","#mu");
      caption=caption.ReplaceAll("2",",=2 ");
      caption=caption.ReplaceAll("3",",=3 ");
      caption=caption.ReplaceAll("4",",=4 ");
      pt->AddText(caption);
      pt->Draw();

      //overall header
      if(icat==0)
	{
	  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(12);
	  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
	  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
	  pt->Draw();

	  //build the legend
	  TLegend *leg=new TLegend(0.2,0.7,0.5,0.94);
     	  leg->SetBorderSize(0);
	  leg->SetFillStyle(0);
	  leg->SetTextFont(42);
	  leg->AddEntry("data","data","p");
	  leg->AddEntry(preFitModel[tag],preFitModel[tag]->GetTitle(),"l");
	  leg->AddEntry(postExclusiveFitModel[tag],postExclusiveFitModel[tag]->GetTitle(),"l");
	  leg->AddEntry(postInclusiveFitModel[tag],postInclusiveFitModel[tag]->GetTitle(),"l");
	  leg->Draw();
	}

      //model
      modelc->cd(icat+1);
      TLegend *leg=0;
      if(icat==maxJets_-2)
	{
	  leg=new TLegend(0.2,0.7,0.5,0.94);
     	  leg->SetBorderSize(0);
	  leg->SetFillStyle(0);
	  leg->SetTextFont(42);
	}
      THStack *stack=new THStack(("modelstack_"+tag).c_str(),("modelstack_"+tag).c_str());
      for(size_t ih=0; ih<preFitModelFunc[tag].size(); ih++)
	{
	  TH1 *preH= preFitModelFunc[tag][ih];
	  preH->SetFillStyle(1001);
	  preH->SetFillColor(preH->GetLineColor());
	  stack->Add(preH,"HIST");
	  if(leg) leg->AddEntry(preH,preH->GetTitle(),"f");
	}
      stack->Draw("");
      stack->GetXaxis()->SetTitle( fitTypeTitle_ );
      stack->GetYaxis()->SetTitle("Probability");
      stack->GetYaxis()->SetTitleOffset(1.2);

      if(leg) leg->Draw();
      
      pt = new TPaveText(0.6,0.85,0.85,0.94,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->AddText(caption);
      pt->Draw();

      //overall header
      if(icat==0)
	{
	  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(12);
	  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
	  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
	  pt->Draw();
	}
      
      //post-fit nuisance parameters
      pfnc->cd(icat+1)->SetBottomMargin(0.3);
      if(curResults_[tag].postFitNuisGr)
	{
	  curResults_[tag].postFitNuisGr->Draw("b");
	  curResults_[tag].postFitNuisGr->GetXaxis()->CenterTitle();
	  curResults_[tag].postFitNuisGr->GetYaxis()->SetRangeUser(-2.5,2.5);
	  curResults_[tag].postFitNuisGr->LabelsDeflate("X");
	  curResults_[tag].postFitNuisGr->LabelsOption("v");
	}

      pt = new TPaveText(0.6,0.85,0.85,0.94,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->AddText(caption);
      pt->Draw();
      
      if(icat==0)
	{
	  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(12);
	  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
	  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
	  pt->Draw();
	}
    }

  //extended model canvas
  TCanvas *extc = new TCanvas("extc","extc",970,600);
  extc->cd();
  {
    TLegend *leg=new TLegend(0.6,0.96,0.95,1.0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetNColumns(3);
    modelExtended->Draw("hist");             leg->AddEntry(modelExtended,"Pre-fit model","f");
    postFitModelExtended->Draw("histsame");  leg->AddEntry(postFitModelExtended,"Post-fit model","f");
    dataExtended->Draw("e1same");            leg->AddEntry(dataExtended,"data","f");
    Float_t maxY(max(modelExtended->GetMaximum(),dataExtended->GetMaximum())*1.1);
    modelExtended->GetYaxis()->SetRangeUser(1e-1,maxY);
    leg->Draw();
    TPaveText *pt = new TPaveText(0.12,0.96,0.6,1.0,"brNDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
    else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
    pt->Draw();
    extc->Modified();
    extc->Update();
    
    //auxiliary lines and labels
    TLine *l=new TLine;
    l->SetLineStyle(7);
    l->DrawLine(15,0,15,maxY);
    l->DrawLine(30,0,30,maxY);
    for(int i=0; i<9; i++)
      {
	TString label("="); label += (i%3+2); label+=" jets";
	TPaveText *labelt=new TPaveText(2+5*i,maxY*0.95,3+5*i,maxY,"");
	labelt->SetBorderSize(0);
	labelt->SetFillStyle(0);
	labelt->AddText(label);
	labelt->SetTextFont(52);
	labelt->SetTextSize(0.03);
	labelt->Draw();
      }
    for(int i=0; i<3; i++)
      {
	TString label("ee events");
	if(i==1) label="e#mu events";
	else if(i==2) label="#mu#mu events";
	TPaveText *labelt=new TPaveText(6+15*i,maxY*0.9,7+15*i,maxY*0.95,"");
	labelt->SetBorderSize(0);
	labelt->SetFillStyle(0);
	labelt->SetTextFont(42);
	labelt->SetTextAlign(12);
	labelt->SetTextSize(0.03);
	labelt->AddText(label);
	labelt->Draw();
      }
  }

  saveGraphicalResult(extc,"DataSlicesExtended");
  saveGraphicalResult(c,"DataSlicesFit");
  saveGraphicalResult(modelc,"DataSlicesModel");
  saveGraphicalResult(pfnc,"ExclusivePostFitNuisances");

  bool is2Dfit( mc_->GetParametersOfInterest()->getSize()>1 );
  TString plrXtitle( fitTypeTitle_ );
  TString plrYtitle( "-log #lambda" );
  if(is2Dfit){
    TObjArray *tkns=fitTypeTitle_.Tokenize(";");
    plrXtitle=tkns->At(0)->GetName();
    plrYtitle=tkns->At(1)->GetName();
  }
  TCanvas *excllc = new TCanvas("excllc","excllc",1200,400);
  excllc->Divide(3,1);
  TCanvas *excpfnc = new TCanvas("excpfnc","excpfnc",800,1200);
  excpfnc->Divide(1,3);
  for(int ich=0; ich<3; ich++)
    {
      string ch("ee"),chTitle("ee"); 
      if(ich==1) { ch="mumu"; chTitle="#mu#mu"; }
      if(ich==2) { ch="emu";  chTitle="e#mu";   }

      //exclusive and combined PLR plots
      excllc->cd(ich+1);

      TLegend *leg=0;
      if(ich==0) 
	{
	  leg=new TLegend(0.2,0.7,0.5,0.94);
     	  leg->SetBorderSize(0);
	  leg->SetFillStyle(0);
	  leg->SetTextFont(42);
	  leg->SetNColumns(2);
	}
      
      if(curResults_[ch].plrGr)
	{
	  curResults_[ch].plrGr->Draw("al");
	  curResults_[ch].plrGr->SetTitle("combined");
	  curResults_[ch].plrGr->SetName((ch+"_plr").c_str());
	  if(!is2Dfit) curResults_[ch].plrGr->GetYaxis()->SetRangeUser(0,10);
	  if(leg) leg->AddEntry(curResults_[ch].plrGr,"combined","f");
	  curResults_[ch].plrGr->GetXaxis()->SetTitle(plrXtitle);
	  curResults_[ch].plrGr->GetYaxis()->SetTitle(plrYtitle);
	}

      int igr(0);
      int colors[]={2,8,9};
      for(std::map<std::string,FitResult_t>::iterator rIt=curResults_.begin(); rIt!=curResults_.end(); rIt++)
	{
	  if(rIt->first.find(ch)==string::npos || rIt->first==ch) continue;

	  TGraph *iplr=rIt->second.plrGr;
	  if(iplr==0) continue;
	  TString subCatTitle("");
	  if(rIt->first.find("2")!=string::npos) subCatTitle += "=2 jets";
	  if(rIt->first.find("3")!=string::npos) subCatTitle += "=3 jets";
	  if(rIt->first.find("4")!=string::npos) subCatTitle += "=4 jets";
	  iplr->SetName((rIt->first+"_plr").c_str());
	  iplr->SetLineStyle(7);
	  iplr->SetLineColor(colors[igr]);
	  iplr->Draw("l");
	  if(leg) leg->AddEntry(iplr,subCatTitle,"f");
	  igr++;
	}

      if(leg)leg->Draw();

      TPaveText *pt = new TPaveText(0.6,0.85,0.85,0.94,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->AddText(chTitle.c_str());
      pt->Draw();

      if(ich==0)
	{
	  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(12);
	  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
	  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
	  pt->Draw();
	}


      //post-fit nuisance parameters
      excpfnc->cd(ich+1)->SetBottomMargin(0.3);
      if(curResults_[ch].postFitNuisGr)
	{
	  curResults_[ch].postFitNuisGr->Draw("b");
	  curResults_[ch].postFitNuisGr->GetXaxis()->CenterTitle();
	  curResults_[ch].postFitNuisGr->GetYaxis()->SetRangeUser(-2.5,2.5);
	  curResults_[ch].postFitNuisGr->LabelsDeflate("X");
	  curResults_[ch].postFitNuisGr->LabelsOption("v");
	}

      pt = new TPaveText(0.6,0.85,0.85,0.94,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->AddText(chTitle.c_str());
      pt->Draw();

      if(ich==0)
	{
	  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(0);
	  pt->SetFillStyle(0);
	  pt->SetTextAlign(12);
	  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
	  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
	  pt->Draw();
	}

    }
  saveGraphicalResult(excllc,"ExclusiveFitPLR");
  saveGraphicalResult(excpfnc,"PostFitNuisances");

  //inclusive likelihood plot
  //the result is printed out
  float poiFit    = curResults_["inclusive"].poiFit;
  float poiFitElo = poiFit-curResults_["inclusive"].poiFitLoLim;
  float poiFitEHi = curResults_["inclusive"].poiFitUpLim-poiFit;

  TCanvas *incllc = new TCanvas("incllc","incllc",800,800);
  incllc->SetTopMargin(0.15);
  TLegend *leg=new TLegend(0.2,0.85,0.6,0.95);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42); 
  leg->SetNColumns(2);
  leg->SetTextSize(0.05);
  if(curResults_["inclusive"].plrGr)
    {
      curResults_["inclusive"].plrGr->Draw("al");
      curResults_["inclusive"].plrGr->SetTitle("combined");
      curResults_["inclusive"].plrGr->GetXaxis()->SetTitle(plrXtitle);
      curResults_["inclusive"].plrGr->GetYaxis()->SetTitle(plrYtitle);
      if(!is2Dfit) curResults_["inclusive"].plrGr->GetYaxis()->SetRangeUser(0,10);
      leg->AddEntry(curResults_["inclusive"].plrGr,"combined","f");
    }

  int colors[]={2,8,9};
  for(int ich=0; ich<3; ich++)
    {
      string ch("ee"),chTitle("ee"); 
      if(ich==1) { ch="mumu"; chTitle="#mu#mu"; }
      if(ich==2) { ch="emu";  chTitle="e#mu";   }
      if(curResults_[ch].plrGr==0) continue; 

      TGraph *plrGr=(TGraph *)curResults_[ch].plrGr->Clone((ch+"_combplr").c_str());
      plrGr->SetLineStyle(7);
      plrGr->SetLineColor(colors[ich]);
      plrGr->Draw("l");
      leg->AddEntry(plrGr,chTitle.c_str(),"f");
    }
  leg->Draw();
  
  TPaveText *pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.05);
  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
  pt->Draw();

  //draw the data and the model sum in a sub-pad
  TPad *npad = new TPad("llpad","ll", 0.4, 0.64, 0.8, 0.94);
  npad->Draw();
  npad->cd();
  RooPlot *frame = bmult->frame();
  data_->plotOn(frame,DrawOption("pz"));
  frame->Draw();
  frame->GetXaxis()->SetNdivisions(maxJets_+1);
  frame->GetXaxis()->SetTitle("b-tag multiplicity");
  frame->GetYaxis()->SetTitle("Events");
  for(int ibin=1; ibin<=frame->GetXaxis()->GetNbins(); ibin++)
    {
      TString label("="); label+=(ibin-1);
      frame->GetXaxis()->SetBinLabel(ibin,label);
    }
  frame->GetYaxis()->SetRangeUser(0.01,data_->sumEntries()*0.75);
  frame->GetYaxis()->SetTitleOffset(1.4);
  
  TH1 *modelTotalH=0;
  for(std::map<string, TH1F *>::iterator hIt=postInclusiveFitModel.begin(); hIt!=postInclusiveFitModel.end(); hIt++)
    {
      if(modelTotalH==0) { modelTotalH = (TH1 *)hIt->second->Clone("totalmodel");  modelTotalH->SetDirectory(0); }
      else               { modelTotalH->Add(hIt->second); }
    }
  modelTotalH->Draw("histsame");

  pt = new TPaveText(0.1,0.8,0.9,0.92,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  char buf[200];
  sprintf(buf,"%3.3f < %s < %3.3f",poiFit-poiFitElo,fitTypeTitle_.Data(),poiFit+poiFitEHi);
  pt->AddText(buf);
  pt->Draw("same");
  saveGraphicalResult(incllc,"InclusiveFitPLR");
	
  //post-fit nuisance parameters
  TCanvas *incpfnc = new TCanvas("incpfnc","incpfnc",800,400);
  incpfnc->SetBottomMargin(0.3);
  if(curResults_["inclusive"].postFitNuisGr)
    {
      curResults_["inclusive"].postFitNuisGr->Draw("b");
      curResults_["inclusive"].postFitNuisGr->GetXaxis()->CenterTitle();
      curResults_["inclusive"].postFitNuisGr->GetYaxis()->SetRangeUser(-2.5,2.5);
      curResults_["inclusive"].postFitNuisGr->LabelsDeflate("X");
      curResults_["inclusive"].postFitNuisGr->LabelsOption("v");
    }

  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
  pt->Draw();
  saveGraphicalResult(incpfnc,"InclusivePostFitNuisances");

  //the inclusive model
  TCanvas *incmodelc = new TCanvas("incmodelc","incmodelc",600,600);
  leg=new TLegend(0.2,0.2,0.5,0.4);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  std::vector<TH1F *> totalModel;
  float totalModelNorm(0);
  for(std::map<string,std::vector<TH1F *> >::iterator cIt=preFitModelFunc.begin(); cIt!=preFitModelFunc.end(); cIt++)
    {
      if(cIt->first.find("jets")==string::npos) continue;
      for(size_t ih=0; ih<cIt->second.size(); ih++)
	{
	  if(ih==0) totalModelNorm+=1.0;
	  if(totalModel.size()<=ih)
	    {
	      totalModel.push_back( (TH1F *)cIt->second[ih]->Clone( TString("inc")+cIt->second[ih]->GetName() ) );
	      totalModel[ih]->SetDirectory(0);
	    }
	  else
	    totalModel[ih]->Add(cIt->second[ih] );
	}
    }
  THStack *stack=new THStack("modelstack","modelstack");
  for(size_t ih=0; ih<totalModel.size(); ih++)
    {
      TH1 *preH=totalModel[ih];
      preH->Scale(1./totalModelNorm);
      stack->Add(preH,"HIST");
      leg->AddEntry(preH,preH->GetTitle(),"f");
    }
  stack->Draw("");
  stack->GetXaxis()->SetTitle( fitTypeTitle_ );
  stack->GetYaxis()->SetTitle("Probability");
  stack->GetYaxis()->SetTitleOffset(1.2);
  leg->Draw();
  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
  pt->Draw();
  saveGraphicalResult(incmodelc,"InclusiveModel");

  //finaly a LEP style plot with all the different measurements and the overall combination
  TCanvas *sc=new TCanvas("sc","sc",600,600);
  TGraphAsymmErrors *exclusiveFitsGr=new TGraphAsymmErrors;     exclusiveFitsGr->SetFillStyle(0);      exclusiveFitsGr->SetLineColor(kRed);   exclusiveFitsGr->SetMarkerStyle(1);      exclusiveFitsGr->SetMarkerColor(0);
  TGraphAsymmErrors *exclusiveFitsStatGr=new TGraphAsymmErrors; exclusiveFitsStatGr->SetFillStyle(0);  exclusiveFitsStatGr->SetLineColor(1);  exclusiveFitsStatGr->SetMarkerStyle(20); exclusiveFitsStatGr->SetMarkerColor(1);
  std::vector<TPaveText *>  exclusiveFitNames;
  for(std::map<std::string,FitResult_t>::iterator rIt=curResults_.begin(); rIt!=curResults_.end(); rIt++)
    {
      if(rIt->first=="inclusive") continue;

      float val=rIt->second.poiFit;
      float eLo=val-rIt->second.poiFitLoLim;
      float eHi=rIt->second.poiFitUpLim-val;
      float eStatLo=rIt->second.poiStatErr; //val-rIt->second.poiFitStatLoLim;
      float eStatHi=rIt->second.poiStatErr; //rIt->second.poiFitStatUpLim-val;

      Int_t np=exclusiveFitsGr->GetN();
      exclusiveFitsGr->SetPoint(np,val,np*2);
      exclusiveFitsGr->SetPointError(np,eLo,eHi,0.1,0.1);
      
      exclusiveFitsStatGr->SetPoint(np,val,np*2);
      exclusiveFitsStatGr->SetPointError(np,eStatLo,eStatHi,0,0);

      TPaveText *pt = new TPaveText(1.08,np*2+0.2,1.1,np*2+0.6,"br");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->SetTextSize(0.03);
      //pt->SetTextAlign(42);
      string caption(""), captionSuf("");
      if(rIt->first.find("ee") != string::npos)   caption="ee";
      if(rIt->first.find("mumu") != string::npos) caption="#mu#mu";
      if(rIt->first.find("emu") != string::npos)  caption="e#mu";
      if(rIt->first.find("2")!=string::npos)      captionSuf+="=2 jets";
      else if(rIt->first.find("3")!=string::npos) captionSuf+="=3 jets";
      else if(rIt->first.find("4")!=string::npos) captionSuf+="=4 jets";
      else                                        pt->SetTextFont(62); 
      if(caption=="")         caption=captionSuf;
      else if(captionSuf!="") caption = " "+captionSuf; 
      pt->AddText(caption.c_str());
      exclusiveFitNames.push_back(pt);
    }

  exclusiveFitsGr->Draw("ae2p");
  exclusiveFitsGr->GetXaxis()->SetTitle( fitTypeTitle_ );
  exclusiveFitsGr->GetYaxis()->SetNdivisions(0);
  exclusiveFitsStatGr->Draw("p");

  TGraph *inclusiveFitGr =new TGraph; inclusiveFitGr->SetLineColor(kGray);
  inclusiveFitGr->SetPoint(0,poiFit-poiFitElo,exclusiveFitsGr->GetYaxis()->GetXmin());
  inclusiveFitGr->SetPoint(1,poiFit+poiFitEHi,exclusiveFitsGr->GetYaxis()->GetXmin());
  inclusiveFitGr->SetPoint(2,poiFit+poiFitEHi,exclusiveFitsGr->GetYaxis()->GetXmax());
  inclusiveFitGr->SetPoint(3,poiFit-poiFitElo,exclusiveFitsGr->GetYaxis()->GetXmax());
  inclusiveFitGr->SetPoint(4,poiFit-poiFitElo,exclusiveFitsGr->GetYaxis()->GetXmin());
  inclusiveFitGr->Draw("l");

  for(size_t ntxt=0; ntxt<exclusiveFitNames.size(); ntxt++) exclusiveFitNames[ntxt]->Draw("same");

  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
  pt->Draw();

  TLine *auxLine=new TLine(0.92,21.25,1.1,21.25);
  auxLine->SetLineStyle(7);
  auxLine->DrawLine(0.92,13.25,1.1,13.25);
  auxLine->DrawLine(0.92,5.25,1.1,5.25);

  pt=new TPaveText(1.05,-2.9,1.07,0.59,"br");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextColor(9);
  pt->SetTextAngle(90);
  pt->SetTextFont(0.035);
  pt->AddText("inclusive");
  pt->Draw();

  saveGraphicalResult(sc,"FitResultsSummary");

  return fitResult;
}

void HFCMeasurement::saveGraphicalResult(TCanvas *c,string name)
{
  if(c==0) return;
  c->SaveAs((name+"_"+sampleType_+".png").c_str());
  c->SaveAs((name+"_"+sampleType_+".pdf").c_str());
  c->SaveAs((name+"_"+sampleType_+".C").c_str());
}


//
TH1F *HFCMeasurement::getProjectedModel(TString tag,Double_t cts,TString name, TString title)
{
  int njets(2);
  if(tag.Contains("3")) njets=3;
  if(tag.Contains("4")) njets=4;
  
  TH1F *pdfGr = new TH1F(name,title,maxJets_+1,0,maxJets_+1);
  pdfGr->SetLineColor(kBlue);
  pdfGr->SetLineWidth(2);
  pdfGr->SetDirectory(0);
  for(int itag=0; itag<=njets; itag++)
    {
      char catName[200];
      sprintf(catName,"n%dbtags_%s",itag,tag.Data());	      
      RooAbsPdf *pdf = ((RooSimultaneous *) ws_->pdf("basemodel"))->getPdf(catName);
      if(pdf==0) pdfGr->Fill(itag,0);
      else       pdfGr->Fill(itag, pdf->getVal()*cts);
    }
  
  return pdfGr;
}

//
std::vector<TH1F *> HFCMeasurement::getProbabilityModelFunction(TString tag,TString baseName)
{
  std::vector<TH1F *> res;

  int njets(2);
  if(tag.Contains("3")) njets=3;
  if(tag.Contains("4")) njets=4;
  
  RooRealVar* firstPOI = (RooRealVar*) mc_->GetParametersOfInterest()->first();
  int colors[]={kGray,kGreen,kBlue,kOrange,kRed+2,kGreen+3};
  for(int itag=0; itag<=njets; itag++)
    {
      char catName[200];
      sprintf(catName,"n%dbtags_%s",itag,tag.Data());	      
      RooAbsPdf *pdf = ((RooSimultaneous *) ws_->pdf("basemodel"))->getPdf(catName);

      TString suf("tags"); suf+=itag;
      TString title("=");  title+=itag; title += " b-tags";
      TH1F *pdfGr = new TH1F(baseName+suf,title,200,0,1.2);
      pdfGr->SetLineColor(colors[itag]);
      pdfGr->SetLineWidth(2);
      pdfGr->SetDirectory(0);
      if(firstPOI)
	{
	  for(int xbin=1; xbin<=pdfGr->GetXaxis()->GetNbins(); xbin++)
	    {
	      float x=pdfGr->GetBinCenter(xbin);
	      firstPOI->setVal(x);
	      if(pdf==0) pdfGr->Fill(x,0);
	      else       pdfGr->Fill(x, pdf->getVal());
	    }
	}
      res.push_back(pdfGr);
    }

  //all done here
  return res;
}


//
HFCMeasurement::FitResult_t HFCMeasurement::plrFit(RooDataSet *data, ModelConfig *mc,bool debug,bool systBreakup)
{
  FitResult_t res;
  res.status=false;

  //check inputs
  if(data==0||mc==0) return res;

  //for debugging purposes
  //freezeNuisances(mc,CORRELATEDNUISANCES,true);
  //freezeNuisances(mc,UNCORRELATEDNUISANCES,true);

  RooRealVar* firstPOI = (RooRealVar*) mc_->GetParametersOfInterest()->first();
  
  //the following is based on http://root.cern.ch/root/html/tutorials/roostats/StandardProfileLikelihoodDemo.C.html
  ProfileLikelihoodCalculator pl(*data,*mc_);
  pl.SetConfidenceLevel(0.68); 	
  LikelihoodInterval* interval = pl.GetInterval();
  if(interval==0) { res.status=false; return res; }
  
  res.status      = true;
  res.poiFit      = firstPOI->getVal(); 
  res.poiFitLoLim = interval->LowerLimit(*firstPOI); 
  res.poiFitUpLim = interval->UpperLimit(*firstPOI);
  res.poiErr      = 0.5*(res.poiFitUpLim-res.poiFitLoLim);


  //get post fit nuisance pulls
  TIterator *nuis_params_itr = mc->GetNuisanceParameters()->createIterator();
  TObject *nuis_params_obj;
  while((nuis_params_obj=nuis_params_itr->Next())){
    RooRealVar *nuiVar=(RooRealVar *)nuis_params_obj;
    if(nuiVar==0) continue;
    if(fabs(nuiVar->getVal())<0.01 || isnan(float(nuiVar->getVal())) ) continue;
    res.postFitNuis[nuiVar->GetTitle()]=nuiVar->getVal();
  }
  
  //save post fit pulls as plots
  res.postFitNuisGr=new TH1F(TString("pfnpull_")+data->GetName(),";;Pull (n x #sigma);",res.postFitNuis.size(),0,res.postFitNuis.size());
  res.postFitNuisGr->SetDirectory(0);
  res.postFitNuisGr->SetFillColor(38);
  res.postFitNuisGr->SetStats(0);
  int ibin(1);
  for(std::map<std::string,Double_t>::iterator nIt=res.postFitNuis.begin(); nIt!=res.postFitNuis.end(); nIt++,ibin++)
    {
      res.postFitNuisGr->GetXaxis()->SetBinLabel(ibin,nIt->first.c_str());
      res.postFitNuisGr->SetBinContent(ibin,nIt->second);
    }	

  if(sampleType_!="mc"){
    //statistical only uncertainty
    freezeNuisances(mc,ALLNUISANCES,true);
    ProfileLikelihoodCalculator statpl(*data,*mc);
    statpl.SetConfidenceLevel(0.68); 	
    statpl.GetInterval();
    res.poiStatErr = firstPOI->getError();
    res.uncBreakup["stat"]=res.poiStatErr;
  }

  //loop over the nuisances to estimate breakup (only for data)
  if(systBreakup && sampleType_!="mc")
    {
      cout << "Breaking up systematic uncertainties" << endl;
      nuis_params_itr = mc->GetNuisanceParameters()->createIterator();
      while((nuis_params_obj=nuis_params_itr->Next())){
	RooRealVar *nVar=(RooRealVar *)nuis_params_obj;
	float postfitVal=nVar->getVal();
	nVar->setVal(0.0);
	nVar->setConstant(kFALSE);
	ProfileLikelihoodCalculator ipl(*data,*mc);
	ipl.SetConfidenceLevel(0.68);
	ipl.GetInterval();
	float iunc=sqrt(pow(firstPOI->getError(),2)-pow(res.poiStatErr,2));
	res.uncBreakup[ nVar->GetName() ] = iunc;
	nVar->setVal(postfitVal);
	nVar->setConstant(kTRUE);
      }
    }

 
  //get the likelihood graph if required 
  res.plrGr=0; 
  freezeNuisances(mc,ALLNUISANCES,false);
  if(debug)
    {     
      res.plrGr=new TGraph;
      res.plrGr->SetName("plr");

      //this is the most stupid thing ... to get the likelihood curve
      //from the list of primitives in the canvas, convert the RooCurve to a TGraph
      //even if the class is derived from it (otherwise it crashes)
      TCanvas *c= new TCanvas("tmpc","tmpc",600,600);
      if(mc_->GetParametersOfInterest()->getSize()>1)
	{
	  cout << "More than on POI found: generating a 90% CL contour" << endl;
	  ProfileLikelihoodCalculator plr95(*data,*mc);
	  plr95.SetConfidenceLevel(0.95); 	
	  LikelihoodIntervalPlot plot( plr95.GetInterval());
	  plot.SetRange(0.5,2.0,0.5,2.0);
	  plot.SetNPoints(200);
	  plot.Draw(""); 
	}
      else
	{
	  LikelihoodIntervalPlot plot(interval);
	  //       float rmin=max(firstPOI->getVal()-10*firstPOI->getError(),0.);
	  //       float rmax=min(firstPOI->getVal()+10*firstPOI->getError(),2.0);
	  //       plot.SetRange(rmin,rmax);
	  plot.SetRange(0.9,1.1);
	  if(fitType_==FIT_GAMMAT) plot.SetRange(0.1,3.0);
	  plot.SetNPoints(100);
	  plot.Draw(""); 
	}

      TIter nextpobj(c->GetListOfPrimitives());
      TObject *pobj;
      while ((pobj = nextpobj()))
	{
	  if(pobj==0) break;
	  TString pobjName(pobj->GetName());
	  if(pobjName.BeginsWith("nll")){
	    RooCurve *nllCurve=(RooCurve *)pobj;
	    for(int ipt=0; ipt<nllCurve->GetN(); ipt++)
	      {
		Double_t ix,iy;
		nllCurve->GetPoint(ipt,ix,iy);
		if(fabs(iy)>10 || iy<0) continue;
		res.plrGr->SetPoint(res.plrGr->GetN(),ix,iy);
	      }		
	  }
	  else if(pobjName.BeginsWith("Graph_of_Likelihood")){
	    res.plrGr=(TGraph *) pobj;
	  }
	}

      //final format
      res.plrGr->SetFillStyle(0);
      res.plrGr->SetFillColor(0);
      res.plrGr->SetMarkerStyle(1);
      res.plrGr->SetMarkerColor(kBlue);
      res.plrGr->SetLineWidth(2);
      res.plrGr->SetLineColor(kBlue);
      delete c;
    }
  
  return res;
}

//
void HFCMeasurement::freezeNuisances(ModelConfig *mc,int mode,bool setConstant)
{
  
  TIterator *nuis_params_itr = mc->GetNuisanceParameters()->createIterator();
  TObject *nuis_params_obj;
  while((nuis_params_obj=nuis_params_itr->Next())){
    RooRealVar *nuiVar=(RooRealVar *)nuis_params_obj;
    TString parName(nuiVar->GetName());
    if(mode==CORRELATEDNUISANCES)    if(parName.Contains("_stat") && !parName.Contains("kst_")) continue;
    if(mode==UNCORRELATEDNUISANCES)  if(!parName.Contains("_stat") || parName.Contains("kst_")) continue;
    nuiVar->setConstant(setConstant);
  }
}
   


//
std::string HFCMeasurement::getNuisanceTitle(std::string &name)
{
  if(name=="dy")           return "DY";
  if(name=="gsplit")       return "Extra HF";
  if(name=="br")           return "BR(W#rightarrow l#nu)";
  if(name=="pu")           return "PU";
  if(name=="jer")          return "JER";
  if(name=="jes")          return "JES";
  if(name=="pdf"   )       return "PDF";
  if(name=="ue")           return "UE";
  if(name=="cr")           return "CR";
  if(name=="meps")         return "ME-PS";
  if(name=="q2")           return "Q^{2}";
  if(name=="signal")       return "Signal";
  if(name=="lumi")         return "Lumi";
  if(name=="fakes")        return "Fakes";
  if(name=="seleff")       return "Sel.eff.";
  if(name=="mcstat")       return "MC stats";
  if(name=="kst_stat")     return "single t frac";
  if(name=="pttop")        return "p_{T}(top)";
  if(name=="epsb")         return "#varepsilon_{b}";
  if(name=="epsq")         return "#varepsilon_{q}";
  if(name=="epsqstar")     return "#varepsilon_{q*}";

  string compTitle(name);
  if(name.find("tt_")!=string::npos)   compTitle="#sigma_{t#bar{t}}^{stat}("; 
  if(name.find("fcor_")!=string::npos) compTitle="f_{correct}^{stat}("; 
  if(name.find("ee")!=string::npos)    compTitle+="ee"; 
  if(name.find("emu")!=string::npos)   compTitle+="e#mu"; 
  if(name.find("mumu")!=string::npos)  compTitle+="#mu#mu"; 
  if(name.find("2")!=string::npos)     compTitle+=",=2j)"; 
  if(name.find("3")!=string::npos)     compTitle+=",=3j)"; 
  if(name.find("4")!=string::npos)     compTitle+=",=4j)";
  return compTitle;
}


//
void HFCMeasurement::printConfiguration(std::ostream &os)
{
  //debug
  cout << "******************** HFC model summary  *************************" << endl;
  ws_->Print("v");
  cout << "Categories defined: ";
  for(std::set<std::string>::iterator it = sampleCats_.begin(); it!=sampleCats_.end(); it++) cout << *it << " "; 
  cout << endl;
  cout << "Tagger: " << wp_ << endl;
  cout << "Sample: " << sampleType_ << endl;
  cout << "*****************************************************************" << endl;
}

//
void HFCMeasurement::instantiateSingleTopContribution(RooWorkspace *wspace)
{
  
  //
  // THEORY TOP QUARK WIDTH
  // is a simple constant
  //

  TString Gamma_NLO("(GF*pow(mt,3)/(8*TMath::Pi()*sqrt(2)))*pow(1-pow(mW/mt,2),2)*(1+2*pow(mW/mt,2))*(1-(2*as)/(3*TMath::Pi())*(2.*pow(TMath::Pi(),2)/3.-5./2.))");
  Gamma_NLO.ReplaceAll("mt","x");
  Gamma_NLO.ReplaceAll("GF","[0]");
  Gamma_NLO.ReplaceAll("mW","[1]");
  Gamma_NLO.ReplaceAll("as","[2]");
  TF1 *f_Gamma_NLO=new TF1("Gamma_NLO",Gamma_NLO,166.5,178.5);
  f_Gamma_NLO->SetParameter(0,1.16637E-5);
  f_Gamma_NLO->SetParameter(1,80.399);
  f_Gamma_NLO->SetParameter(2,0.118);
  
  char pdf_gt_th[100];
  sprintf(pdf_gt_th,"Gamma_b[%f]",f_Gamma_NLO->Eval(172.5));
  wspace->factory(pdf_gt_th);
  cout << pdf_gt_th << endl;

  //
  // EXPERIMENTAL CROSS SECTION
  // JHEP 1212 (2012) 035
  //\sigma(t-ch) = 67.2 +/- 6.1 pb = 67.2 +/- 3.7 (stat.) +/- 3.0 (syst.) +/- 3.5 (theor.) +/- 1.5 (lum.) pb
  parseFitConfig("$CMSSW_BASE/src/UserCode/llvv_fwk/test/top/tch_2012_data_cfg.json");

  //
  // THEORY CROSS SECTION
  // parameterized as gaussian(PDF) x Rect(xsec | mu_r,mu_f)  
  //

  //Phys. Rev. D 83 (2011) 091503(R).
  //\sigma(t-ch) = 64.6 +/- 1.3 (scale) +1.4 -0.7 (PDF) pb
  float sigma_st   = 64.6;
  float sigma_st_h = sigma_st+1.3;
  float sigma_st_l = sigma_st-1.3;
  float delta_pdf  = 1.4;
  char pdf_st_th[1000]; 
  sprintf(pdf_st_th,"st_th[%f,%f,%f]",sigma_st,0.,4*sigma_st);
  wspace->factory(pdf_st_th);
  sprintf(pdf_st_th,"GenericPdf::pdf_st_th('(1/(2*(%f-%f)))*(TMath::Erf((%f-@0)/(sqrt(2)*%f))-TMath::Erf((%f-@0)/(sqrt(2)*%f)))',{st_th})",
 	  sigma_st_h,sigma_st_l,sigma_st_h,delta_pdf,sigma_st_l,delta_pdf);
  wspace->factory(pdf_st_th);

  RooArgSet constr;
  constr.add( *(ws_->pdf("pdf_st_th")) );
  if(ws_->set("constr")!=0) constr.add( *(ws_->set("constr")) );
  ws_->defineSet("constr",constr);

  ws_->factory("FormulaVar::r('(@0*@1)/(@2*@3)',{Gamma_b,st_exp,st_th,rb})");

}
