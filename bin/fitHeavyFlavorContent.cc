#include "TList.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSystem.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLegend.h"

#include "UserCode/llvv_fwk/interface/tdrstyle.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/HFCMeasurement.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <vector>

using namespace std;
using namespace RooFit;

typedef std::vector<TString> Proc_t;
typedef std::vector<Float_t> ProcYields_t;
typedef std::vector<TH1D *> ProcHistos_t;


stringstream report;
 
TH1 *getObservedBtagFrom(TString url, TString wp, TString sampleType, TH1F* mcSignalTemplate=0);
TH1F *generateBtagsFor(float R, TString url);
void reweightPredictionToExpectation(TH1F *targetBtag, TH1F *modelBtag);

//check me
void runCalibration(TString url, int fitType, int nuisanceType, int runMode, TString btagWP, 
		    std::map<TString,Double_t> &fitPars, int maxPE, Double_t dataLumi, TString syst,  bool freezeResults);


//
void reweightPredictionToExpectation(TH1 *targetBtag, TH1 *modelBtag)
{
  if(targetBtag==0 || modelBtag==0) return;
  for(Int_t ibin=1; ibin<=targetBtag->GetXaxis()->GetNbins()/5; ibin++)
    {
      Float_t targetTotal(0),modelTotal(0);
      for(Int_t jbin=1; jbin<=5; jbin++){
	Int_t idx=jbin+(ibin-1)*5;
	targetTotal += targetBtag->GetBinContent(idx);
	modelTotal  += modelBtag->GetBinContent(idx);
      }

      //if not counts in model, assume the target
      if(modelTotal==0) {
	for(Int_t jbin=1; jbin<=5; jbin++){
	  Int_t idx=jbin+(ibin-1)*5;
	  modelBtag->SetBinContent(idx,targetBtag->GetBinContent(idx));
	  modelBtag->SetBinError(idx,targetBtag->GetBinError(idx));
	}
      }
      else{
	Float_t sf=targetTotal/modelTotal;
	cout << ibin << " " << sf << endl;
	for(Int_t jbin=1; jbin<=5; jbin++){
	  Int_t idx=jbin+(ibin-1)*5;
	  Float_t cts=modelBtag->GetBinContent(idx);
	  Float_t err=modelBtag->GetBinError(idx);
	  modelBtag->SetBinContent(idx,cts*sf);
	  modelBtag->SetBinError(idx,err*sf);
	}
      }
    }
}


//
void printHelp()
{
  printf("--help       --> print this\n");
  printf("--par        --> fitter configuration file\n");
  printf("--btag       --> tag efficiencies configuration file\n");
  printf("--in         --> input file with b-tag multiplicity distribution\n");
  printf("--syst       --> extra input file with signal b-tag multiplicity distributions for different R scenarios\n");
  printf("--rtogen     --> value of R to generate\n");
  printf("--fit        --> fit type code                         (default: R=0,eb=1,ebvsR=2,vtb=3)\n");
  printf("--ws         --> input file with the workspace\n");
  printf("--npe        --> number of pseudo-experiments to throw (default:1)\n");
  printf("--seed       --> seed for randon number generator (default:0 = TUUID based)\n");
  printf("--study      --> study to perform: lin, syst (default:lin)\n");
  printf("--nosampling --> don't re-sample the shapes for MC studies (use full stats)\n");
  printf("--fast       --> perform a fast preview of the fit results\n");
  printf("command line examples: fitHeavyFlavorContent --in plotter.root --par hfcParams_cfg.json --btag wpParams_cfg.json\n");
  printf("                       fitHeavyFlavorContent --ws WS.root --npe 100\n");
}

//
std::map<TString,Double_t> parseParametersFrom(TString parfileURL)
{
  std::map<TString,Double_t> fitPars;
  ifstream in;
  in.open(parfileURL);
  TString line;
  while (1) {
    in >> line;
    if (!in.good()) break;
    TObjArray *tokens = line.Tokenize(":");
    if(tokens->GetEntriesFast()<2) continue;
    TString key = ((TObjString *)tokens->At(0))->GetString();
    TString val = ((TObjString *)tokens->At(1))->GetString();
    fitPars[key]=val.Atof();
  }
  in.close();
  return fitPars;
}

//
TH1F *generateBtagsFor(float R, TString url)
{

  TFile *inF=TFile::Open(url);
  TH1F *r1=(TH1F*)inF->Get("t#bar{t}systr1/csvLbtagsextended")->Clone("r1");       r1->SetDirectory(0);   r1->Scale(R*R/r1->Integral());
  TH1F *r05a=(TH1F*)inF->Get("t#bar{t}systr05a/csvLbtagsextended")->Clone("r05a"); r05a->SetDirectory(0); r05a->Scale(R*(1-R)/r05a->Integral());
  TH1F *r05b=(TH1F*)inF->Get("t#bar{t}systr05b/csvLbtagsextended")->Clone("r05b"); r05b->SetDirectory(0); r05b->Scale(R*(1-R)/r05b->Integral());
  TH1F *r05=(TH1F *)r05a->Clone("r05"); r05->Add(r05b); r05->SetDirectory(0);
  TH1F *r0=(TH1F*)inF->Get("t#bar{t}systr0/csvLbtagsextended")->Clone("r0");       r0->SetDirectory(0);   r0->Scale((1-R)*(1-R)/r0->Integral());
  inF->Close();

  r1->SetFillStyle(1001); r1->SetFillColor(614);
  r05->SetFillStyle(1001); r05->SetFillColor(824);
  r0->SetFillStyle(1001); r0->SetFillColor(592);


  THStack *total=new THStack;
  total->Add(r0,"hist");
  total->Add(r05,"hist");
  total->Add(r1,"hist");

  TCanvas *c=new TCanvas("c","c",600,600);
  total->Draw();
  total->GetXaxis()->SetTitle(r1->GetXaxis()->GetTitle());
  total->GetYaxis()->SetTitle("Events (a.u.)");
  total->GetYaxis()->SetRangeUser(0,0.2);


  for(int i=0; i<3; i++)
    {
      TString ch("ee");
      if(i==1) ch="#mu#mu";
      if(i==2) ch="e#mu";
      
      TPaveText *pt=new TPaveText(3*5*i+7,0.12,3*5*i+12,0.15);
      pt->SetFillStyle(0);
      pt->SetTextAlign(12);
      pt->SetBorderSize(0);
      pt->SetTextFont(52);
      pt->AddText(ch);
      pt->SetTextSize(0.025);
      pt->Draw();
      
      TLine *l=new TLine(3*5*(1+i),0,3*5*(1+i),0.15);
      l->SetLineStyle(7);
      l->Draw();
    }

  TPaveText *pt=new TPaveText(0.1,0.95,0.6,0.99,"brNDC");
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetBorderSize(0);
  pt->SetTextSize(0.03);
  char buf[50];
  //sprintf(buf,"|V_{tb}|=%3.2f, R=%3.2f",sqrt(R),R);
  sprintf(buf,"R=%3.2f",R);
  pt->AddText("CMS simulation, #sqrt{s}=8 TeV, " + TString(buf));
  pt->Draw();

  TLegend *leg=new TLegend(0.65,0.95,0.96,0.99);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(r0,"R=0","f");
  leg->AddEntry(r05,"R=0.5","f");
  leg->AddEntry(r1,"R=1.0","f");
  leg->SetNColumns(3);
  leg->Draw();

  TString name("csvLbtagsextended_r");
  name += (int)(R*100);

  c->Modified();
  c->Update();
  c->SaveAs(name+".pdf");



  TH1F *h=(TH1F *)total->GetStack()->At( total->GetStack()->GetEntriesFast()-1 )->Clone(name);
  h->SetDirectory(0);
  return h;
}


//
int main(int argc, char* argv[])
{
  // load framework libraries                                                                                          
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  TString url(""),systUrl(""),wsurl("");
  TString fitParFile(""),btagParFile("");
  Int_t fitType(HFCMeasurement::FIT_R);
  TString syst("");
  int maxPE(1);
  int rgSeed(0);
  float rToGen(-1.0);
  TString study("lin");
  bool fastScan(false);
  bool noSampling(false);
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos)              { printHelp();	  return 0; }
      if(arg.find("--in")!=string::npos && i+1<argc)    { url=argv[i+1];         gSystem->ExpandPathName(url);         i++;  printf("in      = %s\n", url.Data());              }
      if(arg.find("--syst")!=string::npos && i+1<argc)  { systUrl=argv[i+1];     gSystem->ExpandPathName(url);         i++;  printf("signal templates in = %s\n", systUrl.Data());              }
      if(arg.find("--rtogen")!=string::npos && i+1<argc){ sscanf(argv[i+1],"%f",&rToGen);                              i++;  printf("R scenario =%f\n", rToGen);                }
      if(arg.find("--ws")!=string::npos && i+1<argc)    { wsurl=argv[i+1];       gSystem->ExpandPathName(wsurl);       i++;  printf("ws      = %s\n", wsurl.Data());            }
      if(arg.find("--par")!=string::npos && i+1<argc)   { fitParFile=argv[i+1];  gSystem->ExpandPathName(fitParFile);  i++;  printf("parFile = %s\n", fitParFile.Data());       }
      if(arg.find("--btag")!=string::npos && i+1<argc)  { btagParFile=argv[i+1]; gSystem->ExpandPathName(btagParFile); i++;  printf("btagFile  = %s\n", btagParFile.Data());    }
      if(arg.find("--npe")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&maxPE);                               i++;  printf("N_{PE}  = %d\n", maxPE);                   }
      if(arg.find("--seed")!=string::npos && i+1<argc)  { sscanf(argv[i+1],"%d",&rgSeed);                              i++;  printf("seed  = %d\n", rgSeed);                   }
      if(arg.find("--fit")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&fitType);                             i++;  printf("fitType  = %d\n", fitType);                }
      if(arg.find("--study")!=string::npos && i+1<argc) { study=argv[i+1];                                             i++;  printf("study    = %s\n",study.Data());            }
      if(arg.find("--fast")!=string::npos)              { fastScan=true;                                                     printf("Fast scan will be run (no plots/tables)\n"); }
      if(arg.find("--nosampling")!=string::npos)        { noSampling=true;                                                   printf("Won't sample the expected b-tag multiplicity\n"); }
    } 
  if(url=="" && wsurl=="")      { printHelp(); return 0; }
  if(url!="" && fitParFile=="") { printHelp(); return 0; }
 
  setTDRStyle();
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);

  //keep RooFit quiet                                                                                                                                                          
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().getStream(0).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);
  

  //
  // STANDARD FIT OR MC CLOSURE
  //
  TH1F *mcSignalTemplate=0;
  if(systUrl!="" && rToGen>=0 && rToGen<=1) mcSignalTemplate=generateBtagsFor(rToGen,systUrl);
  if(rToGen==-1) rToGen=1.0;

  if(wsurl=="")
    {
      //run pseudo-experiments
      HFCMeasurement *fitter=new HFCMeasurement(fitType, fitParFile, btagParFile);
      TString sampleType = fitter->getSampleType(); 
      TString wp         = fitter->getWP();
      
      //get observed/expected b-tag multiplicity from plotter
      TH1 *btagObs=getObservedBtagFrom(url,wp,sampleType,mcSignalTemplate);
      if(btagObs==0) {
	cout << "Unable to find " << sampleType << " b-tag distribution for " << wp << " in " << url << endl;
	return 0;
      }
      
      //generate a random seed: 0 ensures TUUID independet seeds, but may be tricky when running in parallel (better to use job number)
      gRandom->SetSeed(rgSeed);
      
      if(sampleType!="data")
	{
	  //fitter.fitHFCfrom(btagObs,true);
	  
	  std::map<std::string, TH1F *> closureTests;
	  std::map<std::string, TH1F *> closureTestsUnc;
	  TH1 *btagObsSample=0;
	  for(int ipe=1; ipe<=maxPE; ipe++)
	    {
	      
	      //instantiate a new fitter
	      //fitter =new HFCMeasurement(HFCMeasurement::FIT_R, fitParFile,btagParFile);
	      fitter =new HFCMeasurement(fitType, fitParFile,btagParFile);
	      if(btagObsSample==0){
		btagObsSample=(TH1 *)btagObs->Clone("btagobssample");
		btagObsSample->SetDirectory(0);
	      }
	      
	      //sample distribution
	      if(!noSampling){
		btagObsSample->Reset("ICE");
		btagObsSample->FillRandom(btagObs,btagObs->Integral());
	      }
	      
	      //fit it
	      fitter->fitHFCfrom(btagObsSample,!fastScan);
	      
	      //fill closure test histograms
	      std::map<std::string,HFCMeasurement::FitResult_t> &res=fitter->getResults();
	      for(std::map<std::string,HFCMeasurement::FitResult_t>::iterator rit = res.begin(); rit!= res.end(); rit++)
		{
		  //if(rit->first.find("jets")==string::npos) continue;
		  if(rit->second.status==false) continue;
		  TH1F *biasH=0,*uncH=0; 
		  if(closureTests.find(rit->first)==closureTests.end())
		    {
		      biasH = new TH1F(("bias_"+rit->first).c_str(),";bias;Pseudo-experiments",100,-0.0505,0.0495);               
		      biasH->SetDirectory(0);
		      closureTests[rit->first]=biasH;

		      uncH = new TH1F(("unc_"+rit->first).c_str(),";Uncertainty;Pseudo-experiments",100,0,0.06);               
		      uncH->SetDirectory(0);
		      closureTestsUnc[rit->first]=uncH;
		    }
		  else
		    {
		      biasH=closureTests[rit->first];
		      uncH=closureTestsUnc[rit->first];
		    }
		  if(biasH) biasH->Fill(rit->second.poiFit-rToGen);
		  if(uncH)  uncH->Fill(rit->second.poiErr);
		}
	      
	      delete fitter;
	    }

	  //save the histograms to a file
	  TFile *fOut=TFile::Open("HFCClosureTest.root","RECREATE");
	  for(std::map<std::string, TH1F *>::iterator rit=closureTests.begin(); rit!=closureTests.end(); rit++)
	    {
	      rit->second->Write();
	      closureTestsUnc[rit->first]->Write();
	    }
	  fOut->Close();
	} 
      else
	{
	  fitter->fitHFCfrom(btagObs,!fastScan);
	  
	  std::map<std::string,HFCMeasurement::FitResult_t> &res=fitter->getResults();
	  
	  //display results it a tex table
	  string sample[]={"inclusive","ee","mumu","emu"};
	  string sampleTitle[]={"inclusive","$ee$","$\\mu\\mu$","$e\\mu$"};
	  
	  //table header
	  cout << endl
	       << "\\begin{table}" << endl
	       << "\\begin{center}" << endl
	       << "\\caption{}" << endl
	       << "\\label{tab:hfcsummarytable}" << endl
	       << "\\begin{tabular}{lcccc}\\hline\\hline" << endl
	       << "Sample" << flush;
	  for(int isample=0; isample<4; isample++) cout << " & " << sampleTitle[isample];
	  cout << "\\\\\\hline" << endl;
	  
	  //result
	  cout << "$\\mathcal{" << fitter->title() << "}$" << flush; 
	  for(int isample=0; isample<4; isample++) 
	    {
	      float val=res[ sample[isample] ].poiFit;
	      float errHi=res[ sample[isample] ].poiFitUpLim-val;
	      float errLo=val-res[ sample[isample] ].poiFitLoLim;
	      cout << " & " << val << "$^{+" << errHi << "}_{-" << errLo << "}$";
	      //cout << " & " << utils::toLatexRounded(val,0.5*(errHi+errLo));
	    }
	  cout << "\\\\\\hline" << endl;
	  
	  //syst uncertainties
	  std::map<std::string,Double_t> &nuiList=res[sample[0]].uncBreakup;
	  for(std::map<std::string,Double_t>::iterator it=nuiList.begin(); it!=nuiList.end(); it++)
	    {
	      if(it->first=="stat")      cout << "Stat. unc." << flush;
	      else if(it->first=="syst") cout << "Syst. unc." << flush;
	      else                       cout << "~~~~" << it->first << flush; 
	      for(int isample=0; isample<4; isample++) 
		{
		  if(res[sample[isample]].uncBreakup.find(it->first)==res[sample[isample]].uncBreakup.end())
		    cout << " & ";
		  else
		    cout << " & " << res[ sample[isample] ].uncBreakup[it->first];
		}
	      
	      if(it->first=="stat" || it->first=="syst")  cout << "\\\\\\hline" << endl;
	      else                                        cout << "\\\\" << endl;
	    }
	  cout << "\\hline\\end{tabular}" << endl
	       << "\\end{center}" << endl
	       << "\\end{table}" << endl << endl;
	}
    }
  //
  // UNCERTAINTY BREAKUP
  //
  else if (study=="lin")
    {
      //linearity check
      TProfile *linFit=new TProfile("linfit",";Generated R;<Fitter R>;",100,0,1.1,0,1.1);
      Int_t np=10;
      for(int ip=0; ip<=np; ip++)
	{
	  float rgen=ip*1./np;

	  TH1F *templateH=0, *ensembleH=0;
	  for(int ipe=0; ipe<=maxPE; ipe++)
	    {
	      TFile *inF = TFile::Open(wsurl);
	      RooWorkspace *ws=(RooWorkspace *) inF->Get("w");
	      
	      HFCMeasurement *fitter=new HFCMeasurement(ws,4,HFCMeasurement::FIT_R);
	      RooRealVar* firstPOI = (RooRealVar*) ws->set("poi")->first();
	      ws->loadSnapshot("default");
	      firstPOI->setVal(rgen);

	      if(templateH==0) { 
		templateH=fitter->generateBtagObs();
		ensembleH=(TH1F *)templateH->Clone("ensemble");
		ensembleH->SetDirectory(0);
		delete fitter;
		continue;
	      }
	      
	      ensembleH->Reset("ICE");
	      ensembleH->FillRandom(templateH,templateH->Integral());
	      HFCMeasurement::FitResult_t res=fitter->plrFit(ensembleH);
	      linFit->Fill(rgen,res.poiFit,1);
	      
	      inF->Close();
	    }
	  
	  delete templateH;
	  delete ensembleH;
	}

      TCanvas *c=new TCanvas("c","c",600,600);
      linFit->SetMarkerStyle(20);
      linFit->Draw("e1");
      TLine *lin=new TLine(0,0,1,1);
      lin->SetLineColor(kGray);
      lin->Draw();
      TPaveText *pt = new TPaveText(0.1,0.96,0.9,1.0,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetFillStyle(0);
      pt->AddText("CMS simulation");
      pt->Draw();
      c->SaveAs("LinearityTest.png");
      c->SaveAs("LinearityTest.pdf");
      c->SaveAs("LinearityTest.C");
    }
  else
    {

      std::map<string,HFCMeasurement::FitResult_t> res;

      //get workspace and model from file
      TFile *inF = TFile::Open(wsurl);
      RooWorkspace *ws     = (RooWorkspace *) inF->Get("w");
      RooStats::ModelConfig *mc      = (RooStats::ModelConfig *) ws->obj("mc");
      RooDataSet *data     = (RooDataSet *) ws->data("data");

      //central fit
      HFCMeasurement *fitter=new HFCMeasurement(ws,4,HFCMeasurement::FIT_R);
      res["central"] = fitter->plrFit(data,mc,true);
      RooArgSet *nullParams = (RooArgSet *)ws->allVars().snapshot();
      ws->saveSnapshot("bestfit",*nullParams,kTRUE);

      //stat only
      //ws->loadSnapshot("bestfit");
      ws->loadSnapshot("default");
      TIterator *nuisItr_k = ws->set("nuisances")->createIterator();
      RooRealVar *nuis_k=0;
      while((nuis_k=(RooRealVar *)nuisItr_k->Next())) nuis_k->setConstant(kTRUE);
      res["stat"] = fitter->plrFit(data,mc,false);      
      
      //individual syst fit
      nuisItr_k = ws->set("nuisances")->createIterator();
      nuis_k=0;
      while((nuis_k=(RooRealVar *)nuisItr_k->Next())){
	
	//reset the nuisances
	//	ws->loadSnapshot("bestfit");
	ws->loadSnapshot("default");

	//fix all except the one being tested
	TIterator *nuisItr_j = ws->set("nuisances")->createIterator();
	RooRealVar *nuis_j=0;
	while((nuis_j=(RooRealVar *)nuisItr_j->Next())){
	  if(nuis_j==0) continue;
	  if(nuis_j->GetName()!=nuis_k->GetName()) nuis_j->setConstant(kTRUE);
	  else                                     nuis_j->setVal(0);
	}

	res[nuis_k->GetTitle()] = fitter->plrFit(data,mc,false);      
      }
      
      inF->Close();
      
      //printout table
      float sumDiff(0),sumUnc(0);
      for(std::map<string,HFCMeasurement::FitResult_t>::iterator resIt=res.begin();  resIt!=res.end(); resIt++)
	{
	  if(resIt->first=="central") continue;
	  cout << resIt->first << " " << resIt->second.poiFit-res["central"].poiFit << " " << resIt->second.poiFitLoLim << " - " << resIt->second.poiFitUpLim << endl; 
	  if(resIt->first=="stat") continue;
	  sumDiff += pow(resIt->second.poiFit-res["central"].poiFit,2);
	  sumUnc += pow(0.5*(resIt->second.poiFitLoLim-resIt->second.poiFitUpLim),2);
	}
      cout << "--------------------------" << endl;
      cout << "total syst" << sqrt(sumDiff) << " " << sqrt(sumUnc) << endl;
      cout << "--------------------------" << endl;
      cout << "central" << res["central"].poiFit << " +" << res["central"].poiFitUpLim-res["central"].poiFit << " -" << res["central"].poiFit-res["central"].poiFitLoLim << endl;  
    }
}

//
TH1 *getObservedBtagFrom(TString url, TString wp, TString sampleType, TH1F* mcSignalTemplate)
{
  
  TH1 *btagObs=0;
  TFile *fIn = TFile::Open(url);
  if(fIn==0) return btagObs;

  TList *dirs=fIn->GetListOfKeys();
  for(int iproc=0; iproc<dirs->GetEntries(); iproc++)
    {
      TString iDir(dirs->At(iproc)->GetName());
      if(sampleType!="data" && iDir=="data") continue;
      if(sampleType=="data" && iDir!="data") continue;
      TH1 *h    = (TH1 *)fIn->Get(iDir+"/"+wp+"btagsextended");
      if(h==0) continue;
      cout << "Accepting " << iDir << " for " << wp << " and sample type is " << sampleType << endl;
      if(sampleType!="data" && iDir=="t#bar{t}" && mcSignalTemplate!=0){
	cout << "Replacing nominal sample with artificial signal template. Re-weighting yields per categories." << endl;
	//mcSignalTemplate->Scale(h->Integral()/mcSignalTemplate->Integral());
	reweightPredictionToExpectation(h,mcSignalTemplate);
	h=mcSignalTemplate;
      }
      if(btagObs==0) { btagObs = (TH1 *)h->Clone("btagobs"); btagObs->SetDirectory(0); }
      else           { btagObs->Add(h); }
    }
  fIn->Close();
  
  return btagObs;
}  
