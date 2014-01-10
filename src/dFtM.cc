#include "UserCode/llvv_fwk/interface/dFtM.h"

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
#include "RooStats/ProofConfig.h"

#include "RooNumIntConfig.h"
#include "RooNLLVar.h"
#include "RooConstVar.h"
#include "RooExtendPdf.h"

#include "TGraphErrors.h"
#include "TSystem.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;


//
dFtM::dFtM(int fitType, TString flavCfg, TString effCfg, TString btagsUrl) : 
  fitType_(fitType),
  ws_(new RooWorkspace("w")),
  mc_(0),
  data_(0)
{ 
  RooArgSet poi;
  ws_->factory("musignal[1.0]");
  //ws_->factory("musignal[1.0,0.5,2]");
  //poi.add( *ws_->var("musignal") );

  //parse configuration
  parseFitConfig(flavCfg);
  parseFitConfig(effCfg);

  //read b-tag histograms for expectations
  readoutExpectationsFrom(btagsUrl);

  //define the parameters of interest according to the type of the fit and the number of categories
  switch(fitType_)
   {
   case FIT_SFb:         fitTypeTitle_="SF_{b}";        fitTypeName_="sfb";      break;
   case FIT_SFq:         fitTypeTitle_="SF_{q}";        fitTypeName_="sfq";      break;
   case FIT_SFb_AND_SFq: fitTypeTitle_="SF_{b};SF_{q}"; fitTypeName_="sfbvssfq"; break;
   }
  float min_SFl(1.0), max_SFl(1.0),min_SFb(1.0),max_SFb(1.0),min_SFc(1.0), max_SFc(1.0);;
  if(fitType_==FIT_SFq || fitType_==FIT_SFb_AND_SFq) { min_SFl=0.5; max_SFl=2.0; }
  //if(fitType_==FIT_SFc || fitType_==FIT_SFc_AND_SFq) { min_SFc=0.5; max_SFc=2.0; }
  if(fitType_==FIT_SFb || fitType_==FIT_SFb_AND_SFq) { min_SFb=0.5; max_SFb=2.0; }
  std::set<Int_t> jetKin;
  for(std::map<std::string, std::pair<Int_t,Int_t> >::iterator cIt=catToJetKinematics_.begin();
      cIt!=catToJetKinematics_.end();
      cIt++)
    {
      jetKin.insert(cIt->second.first);
      jetKin.insert(cIt->second.second);
    }
  char expBuf[500];
  for(std::set<Int_t>::iterator jkIt=jetKin.begin(); jkIt!=jetKin.end(); jkIt++)
    {
      sprintf(expBuf,"SFbk%d[1.0,%f,%f]",*jkIt,min_SFb,max_SFb);
      RooAbsArg *var=ws_->factory(expBuf);    
      if(min_SFb!=max_SFb)  poi.add(*var);
      //     else                  var->setConstant(true);
      
      sprintf(expBuf,"SFck%d[1.0,%f,%f]",*jkIt,min_SFc,max_SFc);
      var=ws_->factory(expBuf);    
      //if(min_SFc!=max_SFc)  poi.add(*var);
      //var->setConstant(true);

      sprintf(expBuf,"SFlk%d[1.0,%f,%f]",*jkIt,min_SFl,max_SFl);
      ws_->factory(expBuf);    
      var=ws_->factory(expBuf);    
      if(min_SFl!=max_SFl)  poi.add(*var);
      //else                  var->setConstant(true);
    }
  ws_->defineSet("poi",poi);  

  //initiate the fit model 
  initModel();

  //define the dataset
  readoutDataFrom(btagsUrl);
}

//
void dFtM::parseFitConfig(TString url)
{
  //check if path is ok
  gSystem->ExpandPathName(url);
  if(gSystem->AccessPathName(url)) return;

  char expBuf[500]; // free buffer	  
  RooArgSet nuis,constr,globalObs;

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
      else if(jsonF.key[iparam]=="catmap")
	{
	  std::vector<std::string> catKeys = params[iparam].key;
	  for(size_t icat=0; icat<catKeys.size(); icat++)
	    {
	      std::string val=params[iparam][ catKeys[icat] ].toString();
	      Int_t j1(-1),j2(-1);
	      sscanf(val.c_str(),"(%d,%d)",&j1,&j2);
	      if(j1<0 || j2<0) continue;
	      catToJetKinematics_[ catKeys[icat] ] = std::pair<int,int>(j1,j2);
	    }
	  continue;
	}
      
      std::vector<JSONWrapper::Object> &dau = params[iparam].daughters();
      for(size_t icat=0; icat<dau.size(); icat++)
	{
	  JSONWrapper::Object &descript=dau[icat];
	  string param=jsonF.key[iparam];
	  //bool skipModifiers(param=="absepsb" || param=="absepsq" || param=="absepsbc");
	  bool skipModifiers(true);
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

	      //modifiers
	      int varCntr(0);    
	      for(size_t iunc=0; iunc<uncs.size(); iunc++)
		{
		  if(uncs[iunc]=="val") continue;

		  //add new constraint if required
		  RooRealVar *nuisVar=ws_->var(uncs[iunc].c_str());
		  RooRealVar *modVar=ws_->var((uncs[iunc]+"_sigma").c_str());

		  if(nuisVar==0)
		    {
		      sprintf(expBuf,"Gaussian::%s_constr(0.0,%s[0,-5,5],1.0)",uncs[iunc].c_str(),uncs[iunc].c_str());
		      RooGaussian *nuisConstr=(RooGaussian *)ws_->factory(expBuf);
		      constr.add( *nuisConstr );
		      
		      //this is a global observable will be used for toys
		      //cf. https://indico.desy.de/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=6083
		      //cf. https://twiki.cern.ch/twiki/bin/view/RooStats/RooStatsTutorialsAugust2012#Create_Poisson_Counting_model
		      //RooRealVar *globalObsVar = ws_->var( (uncs[iunc]+"0").c_str() );         
		      //globalObsVar->setConstant(true);
		      //globalObs.add( *globalObsVar );

		      nuisVar = ws_->var( uncs[iunc].c_str() );
		      if(skipModifiers) nuisVar->setConstant(true);
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
void dFtM::readoutDataFrom(TString url)
{
  TFile *inF=TFile::Open(url);
  if(inF==0 || inF->IsZombie()){
    cout << "[dFtM::readoutDataFrom][Warning] unable to open file @ " << url << endl;
    return;
  }

  //instantiate the dataset                                                                                                                                                                      
  RooRealVar *bmult    = ws_->var("bmult");
  RooRealVar *bmultObs = ws_->var("bmultobs");
  RooCategory *sample  = ws_->cat("sample");
  data_ = new RooDataSet("data","data",RooArgSet( *bmult, *bmultObs, *sample), RooFit::WeightVar(*bmultObs) );
  
  for(std::map<std::string, std::pair<Int_t,Int_t> >::iterator it =catToJetKinematics_.begin();
      it!=catToJetKinematics_.end();
      it++)
    {
      TString hname(it->first.c_str());  hname.ReplaceAll("k","btvkin"); hname+= wp_;
      TH1F *data   = (TH1F *)inF->Get(hname);
      if(data==0) continue;

      for(int xbin=1; xbin<=data->GetXaxis()->GetNbins(); xbin++)
	{
	  TString ch("ee");
	  if(xbin>3) ch="mumu";
	  if(xbin>6) ch="emu";

	  Int_t nbtags((xbin-1)%3);

	  TString tag(ch); tag+=it->first.c_str(); tag+="_"; tag+=nbtags; tag+="t";

	  //add entry for data
	  Float_t counts(data->GetBinContent(xbin));
	  bmult->setVal(nbtags);
	  bmultObs->setVal(counts);
	  sample->setLabel("n_"+tag);

	  data_->add( RooArgSet(*bmult,*bmultObs,*sample), counts );
	}
    }

  inF->Close();
}


//
void dFtM::readoutExpectationsFrom(TString url)
{
  TFile *inF=TFile::Open(url);
  if(inF==0 || inF->IsZombie()){
    cout << "[dFtM::readoutDataFrom][Warning] unable to open file @ " << url << endl;
    return;
  }

  //read histograms from file
  char expBuf[500];
  RooArgSet nuis,constr,globalObs;
  for(std::map<std::string, std::pair<Int_t,Int_t> >::iterator it =catToJetKinematics_.begin();
      it!=catToJetKinematics_.end();
      it++)
    {
      TString hname(it->first.c_str());  hname.ReplaceAll("k","btvkin"); hname+= wp_;
      TH1F *ttbar  = (TH1F *)inF->Get("ttbarmc_"+hname);
      TH1F *dy     = (TH1F *)inF->Get("dydatadriven_"+hname);
      TH1F *others = (TH1F *)inF->Get("othersmc_"+hname);
      if(ttbar==0 || dy==0 || others==0) continue;
      
      TH1F *totalBckg=(TH1F *) dy->Clone("totalbckg_"+hname);
      totalBckg->Add(others);

      for(int xbin=1; xbin<=ttbar->GetXaxis()->GetNbins(); xbin+=3)
	{
	  TString ch("ee");
	  if(xbin>3) ch="mumu";
	  if(xbin>6) ch="emu";

	  TString tag(ch); tag+=it->first.c_str();
	  
	  //compute expected signal fraction and uncertainty
	  Int_t firstBin((xbin/3)*3+1), lastBin((xbin/3)*3+3); 
	  Double_t ttbarExp,ttbarExpUnc;
	  ttbarExp=ttbar->IntegralAndError(firstBin,lastBin,ttbarExpUnc);
	  Double_t othersExp,othersExpUnc;
	  othersExp=totalBckg->IntegralAndError(firstBin,lastBin,othersExpUnc);
	  Double_t totalExp=ttbarExp+othersExp;
	
	  //build signal fraction expression
	  TString varName("fsignal_"+tag);
	  sprintf(expBuf,"%s_0[%f]",varName.Data(),ttbarExp/totalExp);
	  RooRealVar *nomFracVar=(RooRealVar *)ws_->factory(expBuf);
	  sprintf(expBuf,"%s_sigma[%f]",varName.Data(),ttbarExpUnc/totalExp);
	  RooRealVar *nomFracUnc=(RooRealVar *)ws_->factory(expBuf);
	  nomFracUnc->setConstant(true);
	  sprintf(expBuf,"%s_others_sigma[%f]",varName.Data(),othersExpUnc/totalExp);
	  RooRealVar *nomFracBckgUnc=(RooRealVar *)ws_->factory(expBuf);
	  nomFracBckgUnc->setConstant(true);
	  sprintf(expBuf,"Gaussian::%s_bckg_constr(0.0,nuis_%s_bckg[0,-5,5],1.0)",varName.Data(),varName.Data());
	  RooGaussian *nomFracBckgConstr=(RooGaussian *)ws_->factory(expBuf);
	  constr.add( *nomFracBckgConstr );
	  //RooRealVar *nomFracBckgGlobalObsVar = ws_->var( "nuis_"+varName+"_bckg_0");
	  //nomFracBckgGlobalObsVar->setConstant(true);
	  //globalObs.add( *nomFracBckgGlobalObsVar );
	  RooRealVar *nomFracBckgNuis = ws_->var( "nuis_"+varName +"_bckg" );
	  nuis.add( *nomFracBckgNuis );
	  //	  RooFormulaVar sigFracFormula(varName,
	  //				       "@0*@1*max((1+0*@2*@3)*(1-@4*@5),0.)",
	  //				       RooArgList(*(ws_->var("musignal")),*nomFracVar,
	  //					  *nomFracUnc,*(ws_->var("mcstat")),
	  //					  *nomFracBckgUnc,*nomFracBckgNuis));	
	  //RooFormulaVar sigFracFormula(varName,
	  //			       "@0*@1*max(1-@2*@3,0.)",
	  //			       RooArgList(*(ws_->var("musignal")),*nomFracVar,
	  //					  *nomFracBckgUnc,*nomFracBckgNuis));	
	  RooFormulaVar sigFracFormula(varName,
				       "@0*@1",
				       RooArgList(*(ws_->var("musignal")),*nomFracVar));
	  ws_->import(sigFracFormula);

	  //build background pdf (at this point fixed may consider some variation in the future)
	  varName="bckgmodel_"+tag;
	  Float_t f0= othersExp>0 ? totalBckg->GetBinContent(firstBin)/othersExp :  1e-3;
	  sprintf(expBuf,"%s_0t[%f]",varName.Data(),f0);
	  ws_->factory(expBuf);
	  Float_t f1= othersExp>0 ? totalBckg->GetBinContent(firstBin+1)/othersExp :  1e-3;
	  sprintf(expBuf,"%s_1t[%f]",varName.Data(),f1);
	  ws_->factory(expBuf);
	  Float_t f2= othersExp>0 ? totalBckg->GetBinContent(firstBin+2)/othersExp :  1e-3;
	  sprintf(expBuf,"%s_2t[%f]",varName.Data(),f2);
	  ws_->factory(expBuf);
	}
    }
  

  inF->Close();

  //add to the global sets
  if(ws_->set("nuisances")!=0) nuis.add( *(ws_->set("nuisances")) );
  ws_->defineSet("nuisances",nuis);
  
  if(ws_->set("constr")!=0) constr.add( *(ws_->set("constr")) );
  ws_->defineSet("constr",constr);
  
  if(ws_->set("globalobservables")!=0) globalObs.add( *ws_->set("globalobservables"));
  ws_->defineSet("globalobservables",globalObs);
}


//
void dFtM::initModel()
{
  if(ws_==0) return;
  if(ws_->var("bmultobs")!=0) return;
  
  char expBuf[200];

  //observables
  ws_->factory("bmultobs[0.,99999999999.]");
  sprintf(expBuf,"bmult[0,3]");
  ws_->factory(expBuf);
  ws_->var("bmult")->setBins(3);

  RooCategory sample("sample","sample");
  std::map<string,string> pdfsPerCat; 
  for(std::set<std::string>::iterator cIt = sampleCats_.begin(); cIt!=sampleCats_.end(); cIt++)
    {
      std::string tag=*cIt;
      std::string kinTag = tag.substr ( tag.find("k") );       
      std::pair<Int_t,Int_t> jetKinematics=catToJetKinematics_[kinTag];

      //flavour content variables
      std::map<TString,RooAbsReal *> F;
      F["bb"]=ws_->function(("fbb_"+tag).c_str());
      F["bc"]=ws_->function(("fbc_"+tag).c_str());
      F["bl"]=ws_->function(("fbl_"+tag).c_str());
      F["cb"]=ws_->function(("fcb_"+tag).c_str());
      F["cc"]=ws_->function(("fcc_"+tag).c_str());
      F["cl"]=ws_->function(("fcl_"+tag).c_str());
      F["lb"]=ws_->function(("flb_"+tag).c_str());
      F["lc"]=ws_->function(("flc_"+tag).c_str());
      F["ll"]=ws_->function(("fll_"+tag).c_str());

      //b-tagging variables
      std::vector<RooAbsReal *> SFb,Eb, SFc, Ec, SFl, El;
      for(size_t ijet=0; ijet<2; ijet++)
	{
	  TString ijetKin("k"); 
	  if(ijet==0) ijetKin += jetKinematics.first;
	  if(ijet==1) ijetKin += jetKinematics.second;
	  SFb.push_back( ws_->var("SFb"+ijetKin) );
	  SFc.push_back( ws_->var("SFc"+ijetKin) );
	  SFl.push_back( ws_->var("SFl"+ijetKin) );
	  ijetKin="k_"; 
	  if(ijet==0) ijetKin += jetKinematics.first;
	  if(ijet==1) ijetKin += jetKinematics.second;
	  Eb.push_back( ws_->function("absepsb_"+ijetKin) );
	  Ec.push_back( ws_->function("absepsc_"+ijetKin) );
	  El.push_back( ws_->function("absepsl_"+ijetKin) );
	}

      //the probability function
      //2 b-tags is writen below (1 and 0 just follow the complementary probabilty i.e. SF_i e_i -> (1- SF_i e_i) 
      //"f_bb*SF_b[0]*e_b[0]*SF_b[1]*e_b[1]*+f_bc*SF_b[0]*e_b[0]*SF_c[1]*e_c[1]+f_bl*SF_b[0]*e_b[0]*SF_l[1]*e_l[1]+"
      //"f_cb*SF_c[0]*e_c[0]*SF_b[1]*e_b[1]*+f_cc*SF_c[0]*e_c[0]*SF_c[1]*e_c[1]+f_cl*SF_c[0]*e_c[0]*SF_l[1]*e_l[1]+"
      //"f_lb*SF_l[0]*e_l[0]*SF_l[1]*e_b[1]*+f_lc*SF_l[0]*e_l[0]*SF_c[1]*e_c[1]+f_ll*SF_l[0]*e_l[0]*SF_l[1]*e_l[1]"
      RooArgList sigpdfArgs;
      /*0*/  sigpdfArgs.add(*F["bb"]); sigpdfArgs.add(*(SFb[0])); sigpdfArgs.add(*(Eb[0]));  sigpdfArgs.add(*(SFb[1])); sigpdfArgs.add(*(Eb[1]));
      /*5*/  sigpdfArgs.add(*F["bc"]);                                                       sigpdfArgs.add(*(SFc[1])); sigpdfArgs.add(*(Ec[1])); 
      /*8*/  sigpdfArgs.add(*F["bl"]);                                                       sigpdfArgs.add(*(SFl[1])); sigpdfArgs.add(*(El[1]));
      /*11*/ sigpdfArgs.add(*F["cb"]); sigpdfArgs.add(*(SFc[0])); sigpdfArgs.add(*(Ec[0]));
      /*14*/ sigpdfArgs.add(*F["cc"]); 
      /*15*/ sigpdfArgs.add(*F["cl"]);
      /*16*/ sigpdfArgs.add(*F["lb"]); sigpdfArgs.add(*(SFl[0])); sigpdfArgs.add(*(El[0]));
      /*19*/ sigpdfArgs.add(*F["lc"]); 
      /*20*/ sigpdfArgs.add(*F["ll"]); 

      //2 b-tags
      RooFormulaVar sig2pdf(("sigmodel_"+tag+"_2t").c_str(),
 			    "(@0 *@3*@4+@5 *@6*@7+@8 *@9*@10)*@1*@2+"
			    "(@11*@3*@4+@14*@6*@7+@15*@9*@10)*@12*@13+"
			    "(@16*@3*@4+@19*@6*@7+@20*@9*@10)*@17*@18",
			    sigpdfArgs);
      RooGenericPdf btag2pdf(("model_"+tag+"_2t").c_str(),
      			     "@0*@1+(1-@2)*@3",
      			     RooArgList(*(ws_->function(("fsignal_"+tag).c_str())),sig2pdf,*(ws_->var(("fsignal_"+tag+"_0").c_str())),*(ws_->var(("bckgmodel_"+tag+"_2t").c_str()))));
      sample.defineType(("n_"+tag+"_2t").c_str());
      pdfsPerCat["n_"+tag+"_2t"]=btag2pdf.GetName();
      ws_->import( btag2pdf,RecycleConflictNodes());

      //1 b-tag
      RooFormulaVar sig1pdf(("sigmodel_"+tag+"_1t").c_str(),
 			    "(@0 *@3*@4+@5 *@6*@7+@8 *@9*@10)*(1-@1*@2)+"
 			    "(@0 *(1-@3*@4)+@5 *(1-@6*@7)+@8 *(1-@9*@10))*@1*@2+"
			    "(@11*@3*@4+@14*@6*@7+@15*@9*@10)*(1-@12*@13)+"
			    "(@11*(1-@3*@4)+@14*(1-@6*@7)+@15*(1-@9*@10))*@12*@13+"
			    "(@16*@3*@4+@19*@6*@7+@20*@9*@10)*(1-@17*@18)+"
			    "(@16*(1-@3*@4)+@19*(1-@6*@7)+@20*(1-@9*@10))*@17*@18",
			    sigpdfArgs);
      RooGenericPdf btag1pdf(("model_"+tag+"_1t").c_str(),
			     "@0*@1+(1-@2)*@3",
      			     RooArgList(*(ws_->function(("fsignal_"+tag).c_str())),sig1pdf,*(ws_->var(("fsignal_"+tag+"_0").c_str())),*(ws_->var(("bckgmodel_"+tag+"_1t").c_str()))));
      sample.defineType(("n_"+tag+"_1t").c_str());
      pdfsPerCat["n_"+tag+"_1t"]=btag1pdf.GetName();
      ws_->import( btag1pdf,RecycleConflictNodes() ); 

      //0 b-tags
      RooFormulaVar sig0pdf(("sigmodel_"+tag+"_0t").c_str(),
 			    "(@0 *(1-@3*@4)+@5 *(1-@6*@7)+@8 *(1-@9*@10))*(1-@1*@2)+"
			    "(@11*(1-@3*@4)+@14*(1-@6*@7)+@15*(1-@9*@10))*(1-@12*@13)+"
			    "(@16*(1-@3*@4)+@19*(1-@6*@7)+@20*(1-@9*@10))*(1-@17*@18)",
			    sigpdfArgs);
      RooGenericPdf btag0pdf(("model_"+tag+"_0t").c_str(),
			     "@0*@1+(1-@2)*@3",
      			     RooArgList(*(ws_->function(("fsignal_"+tag).c_str())),sig0pdf,*(ws_->var(("fsignal_"+tag+"_0").c_str())),*(ws_->var(("bckgmodel_"+tag+"_0t").c_str()))));
      sample.defineType(("n_"+tag+"_0t").c_str());
      pdfsPerCat["n_"+tag+"_0t"]=btag0pdf.GetName();
      ws_->import( btag0pdf,RecycleConflictNodes());
    }

  //base model is a simultaneous pdf for each category
  RooSimultaneous *dFtMmodel = new RooSimultaneous("dFtMbase","dFtM",sample);
  for(std::map<string, string>::iterator it=pdfsPerCat.begin(); it!=pdfsPerCat.end(); it++)
    {
      RooAbsPdf *pdf=ws_->pdf(it->second.c_str());
      if(pdf==0)
 	{
 	  cout << "[Warning] failed to link " << it->second << " to category " << it->first << endl;
 	  continue;
 	}
      dFtMmodel->addPdf( *pdf, it->first.c_str() );
    }
  
  //add the event yields so it can be extended                              
  RooRealVar *yields = (RooRealVar *)ws_->factory("yields[1,0,999999999.]");

  //add the constraints if any
  RooAbsPdf *model=(RooAbsPdf *)dFtMmodel;
  if(ws_->set("constr")!=0)
    {
      RooExtendPdf *extmodel = new RooExtendPdf("extmodel","extmodel",*dFtMmodel,*yields);
      RooArgSet modelFactors(*extmodel);
      modelFactors.add( *(ws_->set("constr")) );
      model = (RooAbsPdf *)(new RooProdPdf("dFtM","dFtM",modelFactors));
    }
  else
    {
      model = (RooAbsPdf *)(new RooExtendPdf("model","model",*dFtMmodel,*yields));
    }

  //import the model
  ws_->import(*model);

  //finalize workspace
  ws_->defineSet("observables","bmult,bmultobs,sample");
  RooArgSet baseGlobalObs( *(ws_->var("bmult")), *(ws_->var("bmultobs")) ); 
  baseGlobalObs.add( *(ws_->set("globalobservables") ) );
  ws_->defineSet("globalobservables", baseGlobalObs );
  RooArgSet *nullParams = (RooArgSet *)ws_->allVars().snapshot();
  ws_->saveSnapshot("prefit",*nullParams,kTRUE);
  
  //instantiate a model configurator for RooStats (not necessarily needed...) 
  mc_ = new ModelConfig("mc",ws_);
  mc_->SetPdf(*model);
  mc_->SetParametersOfInterest(*(ws_->set("poi")));
  mc_->SetObservables(*(ws_->set("observables")));
  mc_->SetGlobalObservables(*(ws_->set("globalobservables")));
  mc_->SetNuisanceParameters(*(ws_->set("nuisances")));
}

//
void dFtM::reset()
{
  if(ws_==0) return;
  ws_->loadSnapshot("prefit");
}

//
std::map<TString,TH1F *> dFtM::getExtendedBtagMultiplicityHistograms()
{
  std::map<TString, TH1F *> btagsPerCh;
  size_t nCategsPerChannel(sampleCats_.size()/3);

  //RooRealVar *yields=ws_->var("yields");

  for(size_t ich=0; ich<3; ich++)
    {
      TString ch("ee");
      if(ich==1) ch="mumu";
      if(ich==2) ch="emu";
      TString title(ch); title.ReplaceAll("mu","#mu");

      //create the histograms
      TH1F *btagObs = new TH1F(ch+"btags",title+";b-tag multiplicity;Events",3*nCategsPerChannel,0,3*nCategsPerChannel);
      btagObs->SetDirectory(0);
      btagObs->SetMarkerStyle(20);
      btagObs->SetMarkerColor(1);
      btagObs->SetLineColor(1);
      btagObs->SetLineWidth(2);
      btagsPerCh[ch]=btagObs;
      TH1F *btagExp = (TH1F *)btagObs->Clone(ch+"btagsexp");
      btagExp->SetDirectory(0);
      btagExp->SetMarkerStyle(1);
      btagExp->SetMarkerColor(kBlue);
      btagExp->SetLineColor(kBlue);
      btagExp->SetLineWidth(2);
      btagsPerCh[ch+"_exp"]=btagExp;
      TH1F *btagBckgExp = (TH1F *)btagObs->Clone(ch+"btagsbckgexp");
      btagBckgExp->SetDirectory(0);
      btagBckgExp->SetMarkerStyle(1);
      btagBckgExp->SetMarkerColor(kCyan-5);
      btagBckgExp->SetLineColor(kCyan-5);
      btagBckgExp->SetLineWidth(2);
      btagBckgExp->SetFillStyle(1001);
      btagBckgExp->SetFillColor(kCyan-5);
      btagsPerCh[ch+"_bckgexp"]=btagBckgExp;

      //fill the histograms
      for(size_t icat=1; icat<=nCategsPerChannel; icat++)
	{
	  TString cat(ch+"k"); cat += icat;
 	  
	  for(Int_t ibtags=0; ibtags<=2; ibtags++)
	    {
	      TString binLabel(""); binLabel += ibtags; binLabel+="t";

	      //data
	      TString cut("sample==sample::n_"+cat+"_"+binLabel);
	      RooDataSet *data = (RooDataSet *) data_->reduce(cut);
	      Double_t cts=data->sumEntries();
	      Int_t xbin=(icat-1)*3+ibtags+1;
	      btagObs->SetBinContent(xbin,cts);
	      btagObs->SetBinError(xbin,sqrt(cts));
	      btagObs->GetXaxis()->SetBinLabel(xbin,binLabel);
	      
	      //expectation
	      RooAbsPdf *pdf = (RooAbsPdf *) ws_->pdf("model_"+cat+"_"+binLabel);
	      
	      cut=("sample==sample::n_"+cat+"_0t || sample==sample::n_"+cat+"_1t || sample==sample::n_"+cat+"_2t" );
              cts=((RooDataSet *) data_->reduce(cut))->sumEntries();

	      Double_t expcts=pdf->getVal()*cts;
	      btagExp->SetBinContent(xbin,expcts);
	      btagExp->SetBinError(xbin,0);
	      btagExp->GetXaxis()->SetBinLabel(xbin,binLabel);

	      //background only
	      RooAbsReal *bckgExp=ws_->var("bckgmodel_"+cat+"_"+binLabel);
	      RooAbsReal *fsignal=ws_->function("fsignal_"+cat);
	      expcts=bckgExp->getVal()*cts*(1-fsignal->getVal());
	     
	      btagBckgExp->SetBinContent(xbin,expcts);
              btagBckgExp->SetBinError(xbin,0);
              btagBckgExp->GetXaxis()->SetBinLabel(xbin,binLabel);
	    }
	}
    }

  return btagsPerCh;
}



dFtM::FitResult_t dFtM::fit()
{
  dFtM::FitResult_t ret;
  
  //check inputs
  RooAbsPdf *pdf   =   ws_->pdf("dFtM");
  if(data_==0 || mc_==0 || pdf==0) return ret;

  // TString cut("sample==sample::n_eek1_0t || sample==sample::n_eek1_1t || sample==sample::n_eek1_2t"
  //	      "|| sample==sample::n_mumuk1_0t || sample==sample::n_mumuk1_1t || sample==sample::n_mumuk1_2t"
  //.	      "|| sample==sample::n_emuk1_0t || sample==sample::n_emuk1_1t || sample==sample::n_emuk1_2t");
  //(RooDataSet *) data_->reduce(cut);
  RooDataSet *data = data_;

  //  ProofConfig pc(*ws_, 4, "workers=4", kFALSE);

  ProfileLikelihoodCalculator pl(*data,*mc_); 
  pl.SetConfidenceLevel(0.68);
  LikelihoodInterval* interval = pl.GetInterval();  
  
  RooArgSet *postFitParams = (RooArgSet *)ws_->allVars().snapshot();
  ws_->saveSnapshot("postfit",*postFitParams,kTRUE);


  return ret;
}
//   //reset results
//   curResults_.clear();

//   //exclusive categories
//   std::map<string,RooDataSet *> dataSlice;
//   std::map<string, TH1F *> preFitModel, postExclusiveFitModel, postInclusiveFitModel;
//   std::map<string, std::vector<TH1F *> > preFitModelFunc;
//   std::map<string, TString> combChannelCuts;
//   std::map<int, TString> combMultCuts; 
//   TString allCuts("");
//   for(std::set<string>::iterator cIt = sampleCats_.begin(); cIt != sampleCats_.end(); cIt++)
//     {
//       TString tag=*cIt;
//       int njets(2);
//       if(tag.Contains("3")) njets=3;
//       if(tag.Contains("4")) njets=4;
      
//       //      if(tag.Contains("ee") && njets==4) continue;

//       ////project data
//       //TString cut("sample==sample::n1btags_"+tag+" || sample==sample::n2btags_"+tag);
//       TString cut("sample==sample::n0btags_"+tag+" || sample==sample::n1btags_"+tag+" || sample==sample::n2btags_"+tag);
//       if(njets>2) cut+=" || sample==sample::n3btags_"+tag;
//       if(njets>3) cut+=" || sample==sample::n4btags_"+tag;
//       RooDataSet *data = (RooDataSet *) data_->reduce(cut);
//       dataSlice[tag.Data()]=data;
//       Double_t cts=data->sumEntries();

//       //add to the combined channel cuts
//       string key("ee");
//       if(tag.Contains("mumu")) key="mumu";
//       if(tag.Contains("emu")) key="emu";
//       if(combChannelCuts.find(key)==combChannelCuts.end()) combChannelCuts[key]=cut;
//       else                                                 combChannelCuts[key]=combChannelCuts[key] + " || " + cut; 
//       if(allCuts=="") allCuts=cut;
//       else            allCuts=allCuts + " || " + cut;

//       //add to the combined multiplicity cuts
//       if(combMultCuts.find(njets)==combMultCuts.end()) combMultCuts[njets]=cut;
//       else                                             combMultCuts[njets]=combMultCuts[njets]+ " || " + cut; 

//       //reset to default values
//       resetModelValues();

//       //save the pre-fit snapshot
//       preFitModel[tag.Data()]=getProjectedModel(tag,cts,tag+"_prefitprojpdf","Pre-fit");
//       preFitModelFunc[tag.Data()]=getProbabilityModelFunction(tag,tag+"_prefitprojpdffunc");

//       //profile likelihood method 
//       curResults_[tag.Data()]=plrFit(data,mc_,debug);
      
//       //save the pre-fit snapshot
//       postExclusiveFitModel[tag.Data()]=getProjectedModel(tag,cts,tag+"_postexcfitprojpdf","Post-fit (exc.)");
//     }

//   //exclusive channel results
//   for(std::map<std::string,TString>::iterator cIt=combChannelCuts.begin(); cIt!=combChannelCuts.end(); cIt++)
//     {
//       resetModelValues();
//       RooDataSet *data = (RooDataSet *) data_->reduce(cIt->second);
//       dataSlice[cIt->first]   = data;
//       curResults_[cIt->first] = plrFit(data,mc_,debug);
//     } 

//   //exclusive multiplicity results
//   for(std::map<int,TString>::iterator cIt=combMultCuts.begin(); cIt!=combMultCuts.end(); cIt++)
//     {
//       char tag[200]; sprintf(tag,"=%d jets",cIt->first);
//       resetModelValues();
//       RooDataSet *data = (RooDataSet *) data_->reduce(cIt->second);
//       dataSlice[tag]   = data;
//       curResults_[tag] = plrFit(data,mc_,debug);
//     } 


//   //the inclusive result <- this is the final result to be used
//   resetModelValues();
//   //curResults_["inclusive"]=plrFit(data_,mc_,debug);
//   RooDataSet *data = (RooDataSet *) data_->reduce(allCuts);
//   curResults_["inclusive"]=plrFit(data,mc_,debug,debug);
//   for(std::set<string>::iterator cIt = sampleCats_.begin(); cIt != sampleCats_.end(); cIt++)
//     {
//       string tag=*cIt;
//       postInclusiveFitModel[tag]=getProjectedModel(tag,postExclusiveFitModel[tag]->Integral(),tag+"_postincfitprojpdf","Post-fit (inc.)");
//     }
  
//   //show results of the fit if required
//   if(!debug) return fitResult;

//   RooRealVar *bmult    = ws_->var("bmult");
//   TCanvas *c = new TCanvas("c","c",1200,1200);
//   c->Divide(maxJets_-1,3);
//   TCanvas *modelc = new TCanvas("modelc","modelc",1200,1200);
//   modelc->Divide(maxJets_-1,3);
//   TCanvas *pfnc = new TCanvas("pfnc","pfnc",1200,1200);
//   pfnc->Divide(maxJets_-1,3);
//   int icat=0;
//   TH1F *dataExtended=new TH1F("dataextended",";b-tag multiplicity;Events", 3*3*5, 0.,3*3*5); dataExtended->SetDirectory(0);
//   TH1F *modelExtended=(TH1F *)dataExtended->Clone("modelextended");                          modelExtended->SetDirectory(0);        modelExtended->SetLineColor(30);
//   TH1F *postFitModelExtended=(TH1F *) dataExtended->Clone("pfmodelextended");                postFitModelExtended->SetDirectory(0); postFitModelExtended->SetLineColor(kBlue);
//   for(std::set<string>::iterator cIt = sampleCats_.begin(); cIt != sampleCats_.end(); cIt++, icat++)
//     {
//       string tag=*cIt;

//       //select only fully exclusive categories
//       if(tag.find("jets")==string::npos) { icat--; continue; }

//       c->cd(icat+1);
//       RooPlot *frame = bmult->frame();
//       RooDataSet *data = dataSlice[tag];
//       data->plotOn(frame,Name("data"));
//       frame->Draw();
//       preFitModel[tag]->Draw("histsame");           preFitModel[tag]->SetLineStyle(7); preFitModel[tag]->SetLineColor(30);
//       postExclusiveFitModel[tag]->Draw("histsame"); postExclusiveFitModel[tag]->SetLineColor(30);
//       postInclusiveFitModel[tag]->Draw("histsame");

//       for(int ibin=1; ibin<=modelExtended->GetXaxis()->GetNbins(); ibin++)
// 	{
// 	  Int_t ncounts(0);
// 	  if(ibin<=data->numEntries())
// 	    {
// 	      data->get(ibin-1);
// 	      ncounts = data->weight();
// 	    }
// 	  dataExtended->SetBinContent(ibin+icat*5,ncounts);
// 	  dataExtended->SetBinError(ibin+icat*5,sqrt(ncounts));
// 	  modelExtended->SetBinContent(ibin+icat*5,preFitModel[tag]->GetBinContent(ibin));
// 	  modelExtended->SetBinError(ibin+icat*5,preFitModel[tag]->GetBinError(ibin));
// 	  postFitModelExtended->SetBinContent(ibin+icat*5,postInclusiveFitModel[tag]->GetBinContent(ibin));
// 	  postFitModelExtended->SetBinError(ibin+icat*5,postInclusiveFitModel[tag]->GetBinError(ibin));
// 	}

//       frame->GetXaxis()->SetNdivisions(maxJets_+1);
//       frame->GetXaxis()->SetTitle("b-tag multiplicity");
//       frame->GetYaxis()->SetTitle("Events");
//       for(int ibin=1; ibin<=frame->GetXaxis()->GetNbins(); ibin++)
// 	{
// 	  TString label("="); label+=(ibin-1);
// 	  frame->GetXaxis()->SetBinLabel(ibin,label);
// 	}
//       frame->GetYaxis()->SetRangeUser(0.01,data->sumEntries()*0.55);
//       frame->GetYaxis()->SetTitleOffset(1.4);
	  
//       //category
//       TPaveText *pt = new TPaveText(0.6,0.85,0.85,0.94,"brNDC");
//       pt->SetBorderSize(0);
//       pt->SetFillColor(0);
//       pt->SetFillStyle(0);
//       pt->SetTextFont(42);
//       TString caption=tag.c_str();
//       caption=caption.ReplaceAll("mu","#mu");
//       caption=caption.ReplaceAll("2",",=2 ");
//       caption=caption.ReplaceAll("3",",=3 ");
//       caption=caption.ReplaceAll("4",",=4 ");
//       pt->AddText(caption);
//       pt->Draw();

//       //overall header
//       if(icat==0)
// 	{
// 	  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
// 	  pt->SetBorderSize(0);
// 	  pt->SetFillColor(0);
// 	  pt->SetFillStyle(0);
// 	  pt->SetTextAlign(12);
// 	  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
// 	  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
// 	  pt->Draw();

// 	  //build the legend
// 	  TLegend *leg=new TLegend(0.2,0.7,0.5,0.94);
//      	  leg->SetBorderSize(0);
// 	  leg->SetFillStyle(0);
// 	  leg->SetTextFont(42);
// 	  leg->AddEntry("data","data","p");
// 	  leg->AddEntry(preFitModel[tag],preFitModel[tag]->GetTitle(),"l");
// 	  leg->AddEntry(postExclusiveFitModel[tag],postExclusiveFitModel[tag]->GetTitle(),"l");
// 	  leg->AddEntry(postInclusiveFitModel[tag],postInclusiveFitModel[tag]->GetTitle(),"l");
// 	  leg->Draw();
// 	}

//       //model
//       modelc->cd(icat+1);
//       TLegend *leg=0;
//       if(icat==maxJets_-2)
// 	{
// 	  leg=new TLegend(0.2,0.7,0.5,0.94);
//      	  leg->SetBorderSize(0);
// 	  leg->SetFillStyle(0);
// 	  leg->SetTextFont(42);
// 	}
//       THStack *stack=new THStack(("modelstack_"+tag).c_str(),("modelstack_"+tag).c_str());
//       for(size_t ih=0; ih<preFitModelFunc[tag].size(); ih++)
// 	{
// 	  TH1 *preH= preFitModelFunc[tag][ih];
// 	  preH->SetFillStyle(1001);
// 	  preH->SetFillColor(preH->GetLineColor());
// 	  stack->Add(preH,"HIST");
// 	  if(leg) leg->AddEntry(preH,preH->GetTitle(),"f");
// 	}
//       stack->Draw("");
//       stack->GetXaxis()->SetTitle( fitTypeTitle_ );
//       stack->GetYaxis()->SetTitle("Probability");
//       stack->GetYaxis()->SetTitleOffset(1.2);

//       if(leg) leg->Draw();
      
//       pt = new TPaveText(0.6,0.85,0.85,0.94,"brNDC");
//       pt->SetBorderSize(0);
//       pt->SetFillColor(0);
//       pt->SetFillStyle(0);
//       pt->SetTextFont(42);
//       pt->AddText(caption);
//       pt->Draw();

//       //overall header
//       if(icat==0)
// 	{
// 	  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
// 	  pt->SetBorderSize(0);
// 	  pt->SetFillColor(0);
// 	  pt->SetFillStyle(0);
// 	  pt->SetTextAlign(12);
// 	  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
// 	  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
// 	  pt->Draw();
// 	}
      
//       //post-fit nuisance parameters
//       pfnc->cd(icat+1)->SetBottomMargin(0.3);
//       if(curResults_[tag].postFitNuisGr)
// 	{
// 	  curResults_[tag].postFitNuisGr->Draw("b");
// 	  curResults_[tag].postFitNuisGr->GetXaxis()->CenterTitle();
// 	  curResults_[tag].postFitNuisGr->GetYaxis()->SetRangeUser(-2.5,2.5);
// 	  curResults_[tag].postFitNuisGr->LabelsDeflate("X");
// 	  curResults_[tag].postFitNuisGr->LabelsOption("v");
// 	}

//       pt = new TPaveText(0.6,0.85,0.85,0.94,"brNDC");
//       pt->SetBorderSize(0);
//       pt->SetFillColor(0);
//       pt->SetFillStyle(0);
//       pt->SetTextFont(42);
//       pt->AddText(caption);
//       pt->Draw();
      
//       if(icat==0)
// 	{
// 	  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
// 	  pt->SetBorderSize(0);
// 	  pt->SetFillColor(0);
// 	  pt->SetFillStyle(0);
// 	  pt->SetTextAlign(12);
// 	  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
// 	  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
// 	  pt->Draw();
// 	}
//     }

//   //extended model canvas
//   TCanvas *extc = new TCanvas("extc","extc",970,600);
//   extc->cd();
//   {
//     TLegend *leg=new TLegend(0.6,0.96,0.95,1.0);
//     leg->SetBorderSize(0);
//     leg->SetFillStyle(0);
//     leg->SetTextFont(42);
//     leg->SetNColumns(3);
//     modelExtended->Draw("hist");             leg->AddEntry(modelExtended,"Pre-fit model","f");
//     postFitModelExtended->Draw("histsame");  leg->AddEntry(postFitModelExtended,"Post-fit model","f");
//     dataExtended->Draw("e1same");            leg->AddEntry(dataExtended,"data","f");
//     Float_t maxY(max(modelExtended->GetMaximum(),dataExtended->GetMaximum())*1.1);
//     modelExtended->GetYaxis()->SetRangeUser(1e-1,maxY);
//     leg->Draw();
//     TPaveText *pt = new TPaveText(0.12,0.96,0.6,1.0,"brNDC");
//     pt->SetBorderSize(0);
//     pt->SetFillColor(0);
//     pt->SetFillStyle(0);
//     pt->SetTextAlign(12);
//     if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
//     else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
//     pt->Draw();
//     extc->Modified();
//     extc->Update();
    
//     //auxiliary lines and labels
//     TLine *l=new TLine;
//     l->SetLineStyle(7);
//     l->DrawLine(15,0,15,maxY);
//     l->DrawLine(30,0,30,maxY);
//     for(int i=0; i<9; i++)
//       {
// 	TString label("="); label += (i%3+2); label+=" jets";
// 	TPaveText *labelt=new TPaveText(2+5*i,maxY*0.95,3+5*i,maxY,"");
// 	labelt->SetBorderSize(0);
// 	labelt->SetFillStyle(0);
// 	labelt->AddText(label);
// 	labelt->SetTextFont(52);
// 	labelt->SetTextSize(0.03);
// 	labelt->Draw();
//       }
//     for(int i=0; i<3; i++)
//       {
// 	TString label("ee events");
// 	if(i==1) label="e#mu events";
// 	else if(i==2) label="#mu#mu events";
// 	TPaveText *labelt=new TPaveText(6+15*i,maxY*0.9,7+15*i,maxY*0.95,"");
// 	labelt->SetBorderSize(0);
// 	labelt->SetFillStyle(0);
// 	labelt->SetTextFont(42);
// 	labelt->SetTextAlign(12);
// 	labelt->SetTextSize(0.03);
// 	labelt->AddText(label);
// 	labelt->Draw();
//       }
//   }

//   saveGraphicalResult(extc,"DataSlicesExtended");
//   saveGraphicalResult(c,"DataSlicesFit");
//   saveGraphicalResult(modelc,"DataSlicesModel");
//   saveGraphicalResult(pfnc,"ExclusivePostFitNuisances");

//   bool is2Dfit( mc_->GetParametersOfInterest()->getSize()>1 );
//   TString plrXtitle( fitTypeTitle_ );
//   TString plrYtitle( "-log #lambda" );
//   if(is2Dfit){
//     TObjArray *tkns=fitTypeTitle_.Tokenize(";");
//     plrXtitle=tkns->At(0)->GetName();
//     plrYtitle=tkns->At(1)->GetName();
//   }
//   TCanvas *excllc = new TCanvas("excllc","excllc",1200,400);
//   excllc->Divide(3,1);
//   TCanvas *excpfnc = new TCanvas("excpfnc","excpfnc",800,1200);
//   excpfnc->Divide(1,3);
//   for(int ich=0; ich<3; ich++)
//     {
//       string ch("ee"),chTitle("ee"); 
//       if(ich==1) { ch="mumu"; chTitle="#mu#mu"; }
//       if(ich==2) { ch="emu";  chTitle="e#mu";   }

//       //exclusive and combined PLR plots
//       excllc->cd(ich+1);

//       TLegend *leg=0;
//       if(ich==0) 
// 	{
// 	  leg=new TLegend(0.2,0.7,0.5,0.94);
//      	  leg->SetBorderSize(0);
// 	  leg->SetFillStyle(0);
// 	  leg->SetTextFont(42);
// 	  leg->SetNColumns(2);
// 	}
      
//       if(curResults_[ch].plrGr)
// 	{
// 	  curResults_[ch].plrGr->Draw("al");
// 	  curResults_[ch].plrGr->SetTitle("combined");
// 	  curResults_[ch].plrGr->SetName((ch+"_plr").c_str());
// 	  if(!is2Dfit) curResults_[ch].plrGr->GetYaxis()->SetRangeUser(0,10);
// 	  if(leg) leg->AddEntry(curResults_[ch].plrGr,"combined","f");
// 	  curResults_[ch].plrGr->GetXaxis()->SetTitle(plrXtitle);
// 	  curResults_[ch].plrGr->GetYaxis()->SetTitle(plrYtitle);
// 	}

//       int igr(0);
//       int colors[]={2,8,9};
//       for(std::map<std::string,FitResult_t>::iterator rIt=curResults_.begin(); rIt!=curResults_.end(); rIt++)
// 	{
// 	  if(rIt->first.find(ch)==string::npos || rIt->first==ch) continue;

// 	  TGraph *iplr=rIt->second.plrGr;
// 	  if(iplr==0) continue;
// 	  TString subCatTitle("");
// 	  if(rIt->first.find("2")!=string::npos) subCatTitle += "=2 jets";
// 	  if(rIt->first.find("3")!=string::npos) subCatTitle += "=3 jets";
// 	  if(rIt->first.find("4")!=string::npos) subCatTitle += "=4 jets";
// 	  iplr->SetName((rIt->first+"_plr").c_str());
// 	  iplr->SetLineStyle(7);
// 	  iplr->SetLineColor(colors[igr]);
// 	  iplr->Draw("l");
// 	  if(leg) leg->AddEntry(iplr,subCatTitle,"f");
// 	  igr++;
// 	}

//       if(leg)leg->Draw();

//       TPaveText *pt = new TPaveText(0.6,0.85,0.85,0.94,"brNDC");
//       pt->SetBorderSize(0);
//       pt->SetFillColor(0);
//       pt->SetFillStyle(0);
//       pt->SetTextFont(42);
//       pt->AddText(chTitle.c_str());
//       pt->Draw();

//       if(ich==0)
// 	{
// 	  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
// 	  pt->SetBorderSize(0);
// 	  pt->SetFillColor(0);
// 	  pt->SetFillStyle(0);
// 	  pt->SetTextAlign(12);
// 	  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
// 	  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
// 	  pt->Draw();
// 	}


//       //post-fit nuisance parameters
//       excpfnc->cd(ich+1)->SetBottomMargin(0.3);
//       if(curResults_[ch].postFitNuisGr)
// 	{
// 	  curResults_[ch].postFitNuisGr->Draw("b");
// 	  curResults_[ch].postFitNuisGr->GetXaxis()->CenterTitle();
// 	  curResults_[ch].postFitNuisGr->GetYaxis()->SetRangeUser(-2.5,2.5);
// 	  curResults_[ch].postFitNuisGr->LabelsDeflate("X");
// 	  curResults_[ch].postFitNuisGr->LabelsOption("v");
// 	}

//       pt = new TPaveText(0.6,0.85,0.85,0.94,"brNDC");
//       pt->SetBorderSize(0);
//       pt->SetFillColor(0);
//       pt->SetFillStyle(0);
//       pt->SetTextFont(42);
//       pt->AddText(chTitle.c_str());
//       pt->Draw();

//       if(ich==0)
// 	{
// 	  pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
// 	  pt->SetBorderSize(0);
// 	  pt->SetFillColor(0);
// 	  pt->SetFillStyle(0);
// 	  pt->SetTextAlign(12);
// 	  if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
// 	  else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
// 	  pt->Draw();
// 	}

//     }
//   saveGraphicalResult(excllc,"ExclusiveFitPLR");
//   saveGraphicalResult(excpfnc,"PostFitNuisances");

//   //inclusive likelihood plot
//   //the result is printed out
//   float poiFit    = curResults_["inclusive"].poiFit;
//   float poiFitElo = poiFit-curResults_["inclusive"].poiFitLoLim;
//   float poiFitEHi = curResults_["inclusive"].poiFitUpLim-poiFit;

//   TCanvas *incllc = new TCanvas("incllc","incllc",800,800);
//   incllc->SetTopMargin(0.15);
//   TLegend *leg=new TLegend(0.2,0.85,0.6,0.95);
//   leg->SetBorderSize(0);
//   leg->SetFillStyle(0);
//   leg->SetTextFont(42); 
//   leg->SetNColumns(2);
//   leg->SetTextSize(0.05);
//   if(curResults_["inclusive"].plrGr)
//     {
//       curResults_["inclusive"].plrGr->Draw("al");
//       curResults_["inclusive"].plrGr->SetTitle("combined");
//       curResults_["inclusive"].plrGr->GetXaxis()->SetTitle(plrXtitle);
//       curResults_["inclusive"].plrGr->GetYaxis()->SetTitle(plrYtitle);
//       if(!is2Dfit) curResults_["inclusive"].plrGr->GetYaxis()->SetRangeUser(0,10);
//       leg->AddEntry(curResults_["inclusive"].plrGr,"combined","f");
//     }

//   int colors[]={2,8,9};
//   for(int ich=0; ich<3; ich++)
//     {
//       string ch("ee"),chTitle("ee"); 
//       if(ich==1) { ch="mumu"; chTitle="#mu#mu"; }
//       if(ich==2) { ch="emu";  chTitle="e#mu";   }
//       if(curResults_[ch].plrGr==0) continue; 

//       TGraph *plrGr=(TGraph *)curResults_[ch].plrGr->Clone((ch+"_combplr").c_str());
//       plrGr->SetLineStyle(7);
//       plrGr->SetLineColor(colors[ich]);
//       plrGr->Draw("l");
//       leg->AddEntry(plrGr,chTitle.c_str(),"f");
//     }
//   leg->Draw();
  
//   TPaveText *pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
//   pt->SetBorderSize(0);
//   pt->SetFillColor(0);
//   pt->SetFillStyle(0);
//   pt->SetTextAlign(12);
//   pt->SetTextSize(0.05);
//   if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
//   else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
//   pt->Draw();

//   //draw the data and the model sum in a sub-pad
//   TPad *npad = new TPad("llpad","ll", 0.4, 0.64, 0.8, 0.94);
//   npad->Draw();
//   npad->cd();
//   RooPlot *frame = bmult->frame();
//   data_->plotOn(frame,DrawOption("pz"));
//   frame->Draw();
//   frame->GetXaxis()->SetNdivisions(maxJets_+1);
//   frame->GetXaxis()->SetTitle("b-tag multiplicity");
//   frame->GetYaxis()->SetTitle("Events");
//   for(int ibin=1; ibin<=frame->GetXaxis()->GetNbins(); ibin++)
//     {
//       TString label("="); label+=(ibin-1);
//       frame->GetXaxis()->SetBinLabel(ibin,label);
//     }
//   frame->GetYaxis()->SetRangeUser(0.01,data_->sumEntries()*0.75);
//   frame->GetYaxis()->SetTitleOffset(1.4);
  
//   TH1 *modelTotalH=0;
//   for(std::map<string, TH1F *>::iterator hIt=postInclusiveFitModel.begin(); hIt!=postInclusiveFitModel.end(); hIt++)
//     {
//       if(modelTotalH==0) { modelTotalH = (TH1 *)hIt->second->Clone("totalmodel");  modelTotalH->SetDirectory(0); }
//       else               { modelTotalH->Add(hIt->second); }
//     }
//   modelTotalH->Draw("histsame");

//   pt = new TPaveText(0.1,0.8,0.9,0.92,"brNDC");
//   pt->SetBorderSize(0);
//   pt->SetFillColor(0);
//   pt->SetFillStyle(0);
//   pt->SetTextFont(42);
//   char buf[200];
//   sprintf(buf,"%3.3f < %s < %3.3f",poiFit-poiFitElo,fitTypeTitle_.Data(),poiFit+poiFitEHi);
//   pt->AddText(buf);
//   pt->Draw("same");
//   saveGraphicalResult(incllc,"InclusiveFitPLR");
	
//   //post-fit nuisance parameters
//   TCanvas *incpfnc = new TCanvas("incpfnc","incpfnc",800,400);
//   incpfnc->SetBottomMargin(0.3);
//   if(curResults_["inclusive"].postFitNuisGr)
//     {
//       curResults_["inclusive"].postFitNuisGr->Draw("b");
//       curResults_["inclusive"].postFitNuisGr->GetXaxis()->CenterTitle();
//       curResults_["inclusive"].postFitNuisGr->GetYaxis()->SetRangeUser(-2.5,2.5);
//       curResults_["inclusive"].postFitNuisGr->LabelsDeflate("X");
//       curResults_["inclusive"].postFitNuisGr->LabelsOption("v");
//     }

//   pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
//   pt->SetBorderSize(0);
//   pt->SetFillColor(0);
//   pt->SetFillStyle(0);
//   pt->SetTextAlign(12);
//   if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
//   else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
//   pt->Draw();
//   saveGraphicalResult(incpfnc,"InclusivePostFitNuisances");

//   //the inclusive model
//   TCanvas *incmodelc = new TCanvas("incmodelc","incmodelc",600,600);
//   leg=new TLegend(0.2,0.2,0.5,0.4);
//   leg->SetBorderSize(0);
//   leg->SetFillStyle(0);
//   leg->SetTextFont(42);
//   std::vector<TH1F *> totalModel;
//   float totalModelNorm(0);
//   for(std::map<string,std::vector<TH1F *> >::iterator cIt=preFitModelFunc.begin(); cIt!=preFitModelFunc.end(); cIt++)
//     {
//       if(cIt->first.find("jets")==string::npos) continue;
//       for(size_t ih=0; ih<cIt->second.size(); ih++)
// 	{
// 	  if(ih==0) totalModelNorm+=1.0;
// 	  if(totalModel.size()<=ih)
// 	    {
// 	      totalModel.push_back( (TH1F *)cIt->second[ih]->Clone( TString("inc")+cIt->second[ih]->GetName() ) );
// 	      totalModel[ih]->SetDirectory(0);
// 	    }
// 	  else
// 	    totalModel[ih]->Add(cIt->second[ih] );
// 	}
//     }
//   THStack *stack=new THStack("modelstack","modelstack");
//   for(size_t ih=0; ih<totalModel.size(); ih++)
//     {
//       TH1 *preH=totalModel[ih];
//       preH->Scale(1./totalModelNorm);
//       stack->Add(preH,"HIST");
//       leg->AddEntry(preH,preH->GetTitle(),"f");
//     }
//   stack->Draw("");
//   stack->GetXaxis()->SetTitle( fitTypeTitle_ );
//   stack->GetYaxis()->SetTitle("Probability");
//   stack->GetYaxis()->SetTitleOffset(1.2);
//   leg->Draw();
//   pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
//   pt->SetBorderSize(0);
//   pt->SetFillColor(0);
//   pt->SetFillStyle(0);
//   pt->SetTextAlign(12);
//   if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
//   else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
//   pt->Draw();
//   saveGraphicalResult(incmodelc,"InclusiveModel");

//   //finaly a LEP style plot with all the different measurements and the overall combination
//   TCanvas *sc=new TCanvas("sc","sc",600,600);
//   TGraphAsymmErrors *exclusiveFitsGr=new TGraphAsymmErrors;     exclusiveFitsGr->SetFillStyle(0);      exclusiveFitsGr->SetLineColor(kRed);   exclusiveFitsGr->SetMarkerStyle(1);      exclusiveFitsGr->SetMarkerColor(0);
//   TGraphAsymmErrors *exclusiveFitsStatGr=new TGraphAsymmErrors; exclusiveFitsStatGr->SetFillStyle(0);  exclusiveFitsStatGr->SetLineColor(1);  exclusiveFitsStatGr->SetMarkerStyle(20); exclusiveFitsStatGr->SetMarkerColor(1);
//   std::vector<TPaveText *>  exclusiveFitNames;
//   for(std::map<std::string,FitResult_t>::iterator rIt=curResults_.begin(); rIt!=curResults_.end(); rIt++)
//     {
//       if(rIt->first=="inclusive") continue;

//       float val=rIt->second.poiFit;
//       float eLo=val-rIt->second.poiFitLoLim;
//       float eHi=rIt->second.poiFitUpLim-val;
//       float eStatLo=rIt->second.poiStatErr; //val-rIt->second.poiFitStatLoLim;
//       float eStatHi=rIt->second.poiStatErr; //rIt->second.poiFitStatUpLim-val;

//       Int_t np=exclusiveFitsGr->GetN();
//       exclusiveFitsGr->SetPoint(np,val,np*2);
//       exclusiveFitsGr->SetPointError(np,eLo,eHi,0.1,0.1);
      
//       exclusiveFitsStatGr->SetPoint(np,val,np*2);
//       exclusiveFitsStatGr->SetPointError(np,eStatLo,eStatHi,0,0);

//       TPaveText *pt = new TPaveText(1.08,np*2+0.2,1.1,np*2+0.6,"br");
//       pt->SetBorderSize(0);
//       pt->SetFillColor(0);
//       pt->SetFillStyle(0);
//       pt->SetTextFont(42);
//       pt->SetTextSize(0.03);
//       //pt->SetTextAlign(42);
//       string caption(""), captionSuf("");
//       if(rIt->first.find("ee") != string::npos)   caption="ee";
//       if(rIt->first.find("mumu") != string::npos) caption="#mu#mu";
//       if(rIt->first.find("emu") != string::npos)  caption="e#mu";
//       if(rIt->first.find("2")!=string::npos)      captionSuf+="=2 jets";
//       else if(rIt->first.find("3")!=string::npos) captionSuf+="=3 jets";
//       else if(rIt->first.find("4")!=string::npos) captionSuf+="=4 jets";
//       else                                        pt->SetTextFont(62); 
//       if(caption=="")         caption=captionSuf;
//       else if(captionSuf!="") caption = " "+captionSuf; 
//       pt->AddText(caption.c_str());
//       exclusiveFitNames.push_back(pt);
//     }

//   exclusiveFitsGr->Draw("ae2p");
//   exclusiveFitsGr->GetXaxis()->SetTitle( fitTypeTitle_ );
//   exclusiveFitsGr->GetYaxis()->SetNdivisions(0);
//   exclusiveFitsStatGr->Draw("p");

//   TGraph *inclusiveFitGr =new TGraph; inclusiveFitGr->SetLineColor(kGray);
//   inclusiveFitGr->SetPoint(0,poiFit-poiFitElo,exclusiveFitsGr->GetYaxis()->GetXmin());
//   inclusiveFitGr->SetPoint(1,poiFit+poiFitEHi,exclusiveFitsGr->GetYaxis()->GetXmin());
//   inclusiveFitGr->SetPoint(2,poiFit+poiFitEHi,exclusiveFitsGr->GetYaxis()->GetXmax());
//   inclusiveFitGr->SetPoint(3,poiFit-poiFitElo,exclusiveFitsGr->GetYaxis()->GetXmax());
//   inclusiveFitGr->SetPoint(4,poiFit-poiFitElo,exclusiveFitsGr->GetYaxis()->GetXmin());
//   inclusiveFitGr->Draw("l");

//   for(size_t ntxt=0; ntxt<exclusiveFitNames.size(); ntxt++) exclusiveFitNames[ntxt]->Draw("same");

//   pt = new TPaveText(0.12,0.96,0.9,1.0,"brNDC");
//   pt->SetBorderSize(0);
//   pt->SetFillColor(0);
//   pt->SetFillStyle(0);
//   pt->SetTextAlign(12);
//   if(sampleType_!="data")  pt->AddText("CMS simulation, #sqrt{s}=8 TeV");
//   else                     pt->AddText("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
//   pt->Draw();

//   TLine *auxLine=new TLine(0.92,21.25,1.1,21.25);
//   auxLine->SetLineStyle(7);
//   auxLine->DrawLine(0.92,13.25,1.1,13.25);
//   auxLine->DrawLine(0.92,5.25,1.1,5.25);

//   pt=new TPaveText(1.05,-2.9,1.07,0.59,"br");
//   pt->SetBorderSize(0);
//   pt->SetFillColor(0);
//   pt->SetFillStyle(0);
//   pt->SetTextColor(9);
//   pt->SetTextAngle(90);
//   pt->SetTextFont(0.035);
//   pt->AddText("inclusive");
//   pt->Draw();

//   saveGraphicalResult(sc,"FitResultsSummary");

//   return fitResult;
// }

// void dFtM::saveGraphicalResult(TCanvas *c,string name)
// {
//   if(c==0) return;
//   c->SaveAs((name+"_"+sampleType_+".png").c_str());
//   c->SaveAs((name+"_"+sampleType_+".pdf").c_str());
//   c->SaveAs((name+"_"+sampleType_+".C").c_str());
// }


// //
// TH1F *dFtM::getProjectedModel(TString tag,Double_t cts,TString name, TString title)
// {
//   int njets(2);
//   if(tag.Contains("3")) njets=3;
//   if(tag.Contains("4")) njets=4;
  
//   TH1F *pdfGr = new TH1F(name,title,maxJets_+1,0,maxJets_+1);
//   pdfGr->SetLineColor(kBlue);
//   pdfGr->SetLineWidth(2);
//   pdfGr->SetDirectory(0);
//   for(int itag=0; itag<=njets; itag++)
//     {
//       char catName[200];
//       sprintf(catName,"n%dbtags_%s",itag,tag.Data());	      
//       RooAbsPdf *pdf = ((RooSimultaneous *) ws_->pdf("basemodel"))->getPdf(catName);
//       if(pdf==0) pdfGr->Fill(itag,0);
//       else       pdfGr->Fill(itag, pdf->getVal()*cts);
//     }
  
//   return pdfGr;
// }

// //
// std::vector<TH1F *> dFtM::getProbabilityModelFunction(TString tag,TString baseName)
// {
//   std::vector<TH1F *> res;

//   int njets(2);
//   if(tag.Contains("3")) njets=3;
//   if(tag.Contains("4")) njets=4;
  
//   RooRealVar* firstPOI = (RooRealVar*) mc_->GetParametersOfInterest()->first();
//   int colors[]={kGray,kGreen,kBlue,kOrange,kRed+2,kGreen+3};
//   for(int itag=0; itag<=njets; itag++)
//     {
//       char catName[200];
//       sprintf(catName,"n%dbtags_%s",itag,tag.Data());	      
//       RooAbsPdf *pdf = ((RooSimultaneous *) ws_->pdf("basemodel"))->getPdf(catName);

//       TString suf("tags"); suf+=itag;
//       TString title("=");  title+=itag; title += " b-tags";
//       TH1F *pdfGr = new TH1F(baseName+suf,title,200,0,1.2);
//       pdfGr->SetLineColor(colors[itag]);
//       pdfGr->SetLineWidth(2);
//       pdfGr->SetDirectory(0);
//       if(firstPOI)
// 	{
// 	  for(int xbin=1; xbin<=pdfGr->GetXaxis()->GetNbins(); xbin++)
// 	    {
// 	      float x=pdfGr->GetBinCenter(xbin);
// 	      firstPOI->setVal(x);
// 	      if(pdf==0) pdfGr->Fill(x,0);
// 	      else       pdfGr->Fill(x, pdf->getVal());
// 	    }
// 	}
//       res.push_back(pdfGr);
//     }

//   //all done here
//   return res;
// }


// //
// dFtM::FitResult_t dFtM::plrFit(RooDataSet *data, ModelConfig *mc,bool debug,bool systBreakup)
// {
//   FitResult_t res;
//   res.status=false;

//   //check inputs
//   if(data==0||mc==0) return res;

//   //for debugging purposes
//   //freezeNuisances(mc,CORRELATEDNUISANCES,true);
//   //freezeNuisances(mc,UNCORRELATEDNUISANCES,true);

//   RooRealVar* firstPOI = (RooRealVar*) mc_->GetParametersOfInterest()->first();
  
//   //the following is based on http://root.cern.ch/root/html/tutorials/roostats/StandardProfileLikelihoodDemo.C.html
//   ProfileLikelihoodCalculator pl(*data,*mc_);
//   pl.SetConfidenceLevel(0.68); 	
//   LikelihoodInterval* interval = pl.GetInterval();
//   if(interval==0) { res.status=false; return res; }
  
//   res.status      = true;
//   res.poiFit      = firstPOI->getVal(); 
//   res.poiFitLoLim = interval->LowerLimit(*firstPOI); 
//   res.poiFitUpLim = interval->UpperLimit(*firstPOI);
//   res.poiErr      = 0.5*(res.poiFitUpLim-res.poiFitLoLim);


//   //get post fit nuisance pulls
//   TIterator *nuis_params_itr = mc->GetNuisanceParameters()->createIterator();
//   TObject *nuis_params_obj;
//   while((nuis_params_obj=nuis_params_itr->Next())){
//     RooRealVar *nuiVar=(RooRealVar *)nuis_params_obj;
//     if(nuiVar==0) continue;
//     if(fabs(nuiVar->getVal())<0.01 || isnan(float(nuiVar->getVal())) ) continue;
//     res.postFitNuis[nuiVar->GetTitle()]=nuiVar->getVal();
//   }
  
//   //save post fit pulls as plots
//   res.postFitNuisGr=new TH1F(TString("pfnpull_")+data->GetName(),";;Pull (n x #sigma);",res.postFitNuis.size(),0,res.postFitNuis.size());
//   res.postFitNuisGr->SetDirectory(0);
//   res.postFitNuisGr->SetFillColor(38);
//   res.postFitNuisGr->SetStats(0);
//   int ibin(1);
//   for(std::map<std::string,Double_t>::iterator nIt=res.postFitNuis.begin(); nIt!=res.postFitNuis.end(); nIt++,ibin++)
//     {
//       res.postFitNuisGr->GetXaxis()->SetBinLabel(ibin,nIt->first.c_str());
//       res.postFitNuisGr->SetBinContent(ibin,nIt->second);
//     }	

//   if(sampleType_!="mc"){
//     //statistical only uncertainty
//     freezeNuisances(mc,ALLNUISANCES,true);
//     ProfileLikelihoodCalculator statpl(*data,*mc);
//     statpl.SetConfidenceLevel(0.68); 	
//     statpl.GetInterval();
//     res.poiStatErr = firstPOI->getError();
//     res.uncBreakup["stat"]=res.poiStatErr;
//   }

//   //loop over the nuisances to estimate breakup (only for data)
//   if(systBreakup && sampleType_!="mc")
//     {
//       cout << "Breaking up systematic uncertainties" << endl;
//       nuis_params_itr = mc->GetNuisanceParameters()->createIterator();
//       while((nuis_params_obj=nuis_params_itr->Next())){
// 	RooRealVar *nVar=(RooRealVar *)nuis_params_obj;
// 	float postfitVal=nVar->getVal();
// 	nVar->setVal(0.0);
// 	nVar->setConstant(kFALSE);
// 	ProfileLikelihoodCalculator ipl(*data,*mc);
// 	ipl.SetConfidenceLevel(0.68);
// 	ipl.GetInterval();
// 	float iunc=sqrt(pow(firstPOI->getError(),2)-pow(res.poiStatErr,2));
// 	res.uncBreakup[ nVar->GetName() ] = iunc;
// 	nVar->setVal(postfitVal);
// 	nVar->setConstant(kTRUE);
//       }
//     }

 
//   //get the likelihood graph if required 
//   res.plrGr=0; 
//   freezeNuisances(mc,ALLNUISANCES,false);
//   if(debug)
//     {     
//       res.plrGr=new TGraph;
//       res.plrGr->SetName("plr");

//       //this is the most stupid thing ... to get the likelihood curve
//       //from the list of primitives in the canvas, convert the RooCurve to a TGraph
//       //even if the class is derived from it (otherwise it crashes)
//       TCanvas *c= new TCanvas("tmpc","tmpc",600,600);
//       if(mc_->GetParametersOfInterest()->getSize()>1)
// 	{
// 	  cout << "More than on POI found: generating a 90% CL contour" << endl;
// 	  ProfileLikelihoodCalculator plr95(*data,*mc);
// 	  plr95.SetConfidenceLevel(0.95); 	
// 	  LikelihoodIntervalPlot plot( plr95.GetInterval());
// 	  plot.SetRange(0.5,2.0,0.5,2.0);
// 	  plot.SetNPoints(200);
// 	  plot.Draw(""); 
// 	}
//       else
// 	{
// 	  LikelihoodIntervalPlot plot(interval);
// 	  //       float rmin=max(firstPOI->getVal()-10*firstPOI->getError(),0.);
// 	  //       float rmax=min(firstPOI->getVal()+10*firstPOI->getError(),2.0);
// 	  //       plot.SetRange(rmin,rmax);
// 	  plot.SetRange(0.9,1.1);
// 	  if(fitType_==FIT_GAMMAT) plot.SetRange(0.1,3.0);
// 	  plot.SetNPoints(100);
// 	  plot.Draw(""); 
// 	}

//       TIter nextpobj(c->GetListOfPrimitives());
//       TObject *pobj;
//       while ((pobj = nextpobj()))
// 	{
// 	  if(pobj==0) break;
// 	  TString pobjName(pobj->GetName());
// 	  if(pobjName.BeginsWith("nll")){
// 	    RooCurve *nllCurve=(RooCurve *)pobj;
// 	    for(int ipt=0; ipt<nllCurve->GetN(); ipt++)
// 	      {
// 		Double_t ix,iy;
// 		nllCurve->GetPoint(ipt,ix,iy);
// 		if(fabs(iy)>10 || iy<0) continue;
// 		res.plrGr->SetPoint(res.plrGr->GetN(),ix,iy);
// 	      }		
// 	  }
// 	  else if(pobjName.BeginsWith("Graph_of_Likelihood")){
// 	    res.plrGr=(TGraph *) pobj;
// 	  }
// 	}

//       //final format
//       res.plrGr->SetFillStyle(0);
//       res.plrGr->SetFillColor(0);
//       res.plrGr->SetMarkerStyle(1);
//       res.plrGr->SetMarkerColor(kBlue);
//       res.plrGr->SetLineWidth(2);
//       res.plrGr->SetLineColor(kBlue);
//       delete c;
//     }
  
//   return res;
// }

// //
// void dFtM::freezeNuisances(ModelConfig *mc,int mode,bool setConstant)
// {
  
//   TIterator *nuis_params_itr = mc->GetNuisanceParameters()->createIterator();
//   TObject *nuis_params_obj;
//   while((nuis_params_obj=nuis_params_itr->Next())){
//     RooRealVar *nuiVar=(RooRealVar *)nuis_params_obj;
//     TString parName(nuiVar->GetName());
//     if(mode==CORRELATEDNUISANCES)    if(parName.Contains("_stat") && !parName.Contains("kst_")) continue;
//     if(mode==UNCORRELATEDNUISANCES)  if(!parName.Contains("_stat") || parName.Contains("kst_")) continue;
//     nuiVar->setConstant(setConstant);
//   }
// }
       

//
void dFtM::printConfiguration(std::ostream &os)
{
  //debug
  cout << "******************** HFC model summary  *************************" << endl;
  ws_->Print("v");
  cout << sampleCats_.size() << " categories defined: ";
  for(std::set<std::string>::iterator it = sampleCats_.begin(); it!=sampleCats_.end(); it++) cout << *it << " "; 
  cout << endl;
  cout << "Tagger: " << wp_ << endl;
  cout << "*****************************************************************" << endl;
}

