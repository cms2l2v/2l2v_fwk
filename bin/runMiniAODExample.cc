#include <iostream>
#include <boost/shared_ptr.hpp>

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h" //for svfit

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/llvvObjects.h"
#include "UserCode/llvv_fwk/interface/TMVAUtils.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/PDFInfo.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"
#include "UserCode/llvv_fwk/interface/GammaWeightsHandler.h"
//#include "UserCode/llvv_fwk/interface/i2HDMUtils.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"
#include <Math/VectorUtil.h>

using namespace std;

int main(int argc, char* argv[])
{
	//##############################################
	//########    GLOBAL INITIALIZATION     ########
	//##############################################

	// check arguments
	if(argc<2){ std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl; exit(0); }

	// load framework libraries
	gSystem->Load( "libFWCoreFWLite" );
	AutoLibraryLoader::enable();

	// configure the process
	const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

	bool isMC = runProcess.getParameter<bool>("isMC");
	double xsec = runProcess.getParameter<double>("xsec");
	int mctruthmode=runProcess.getParameter<int>("mctruthmode");

	std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
	TString baseDir    = runProcess.getParameter<std::string>("dirName");
	TString url = TString(argv[1]);
	TString outFileUrl(gSystem->BaseName(url));
	outFileUrl.ReplaceAll(".py","");
	if(mctruthmode!=0) { outFileUrl += "_filt"; outFileUrl += mctruthmode; }
	TString outdir=runProcess.getParameter<std::string>("outdir");
	TString outUrl( outdir );
	gSystem->Exec("mkdir -p " + outUrl);
//	bool filterOnlyEE(false), filterOnlyMUMU(false);
//	if(!isMC)
//	{
//		if(url.Contains("DoubleEle")) filterOnlyEE=true;
//		if(url.Contains("DoubleMu"))  filterOnlyMUMU=true;
//	}
//	bool isSingleMuPD(!isMC && url.Contains("SingleMu"));  
//	bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
//	bool isSignal(isMC && (url.Contains("VBFNLO") || url.Contains("lljj")) );

	TString outTxtUrl= outUrl + "/" + outFileUrl + ".txt";
	FILE* outTxtFile = NULL;
	if(!isMC)outTxtFile = fopen(outTxtUrl.Data(), "w");
	printf("TextFile URL = %s\n",outTxtUrl.Data());

	//##############################################
	//########    INITIATING HISTOGRAMS     ########
	//##############################################
	SmartSelectionMonitor mon;

	//pileup control
	mon.addHistogram( new TH1F( "nvtx",";Vertices;Events",50,-0.5,49.5) ); 
	mon.addHistogram( new TH1F( "nvtxraw",";Vertices;Events",50,-0.5,49.5) ); 

	//##############################################
	//######## GET READY FOR THE EVENT LOOP ########
	//##############################################

	fwlite::ChainEvent ev(urls);
	const unsigned int totalEntries= ev.size();

	//MC normalization (to 1/pb)
	double xsecWeight = xsec/totalEntries;
	if(!isMC) xsecWeight=1.0;

	//pileup weighting: based on vtx for now...
	edm::LumiReWeighting* LumiWeights = NULL;
	utils::cmssw::PuShifter_t PuShifters;
	double PUNorm[] = {1,1,1};
	if(isMC){
		std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
		std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
		std::vector<float> mcPileupDistribution;
		utils::getMCPileupDistributionFromMiniAOD(ev,dataPileupDistribution.size(), mcPileupDistribution);
		while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
		while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
                printf("MCPU:");for(unsigned int i=0;i<  mcPileupDistribution.size();i++){printf("%6.2E ",  mcPileupDistribution[i]);}printf("\n");
                printf("MCPU:");for(unsigned int i=0;i<dataPileupDistribution.size();i++){printf("%6.2E ",dataPileupDistribution[i]);}printf("\n");
       	        gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
		LumiWeights = new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
		PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);
		utils::getPileupNormalization(mcPileupDistribution, PUNorm, LumiWeights, PuShifters);
	}

	//##############################################
	//########           EVENT LOOP         ########
	//##############################################
	//loop on all the events
	printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
	printf("Scanning the ntuple :");
	unsigned int step(totalEntries/50);
	for( unsigned int iev=0; iev<totalEntries; iev++ ){
		if(iev%step==0){printf(".");fflush(stdout);}

		//##############################################   EVENT LOOP STARTS   ##############################################
		ev.to(iev); //load the event content from the EDM file

                reco::VertexCollection vtx;
		fwlite::Handle< reco::VertexCollection > vtxHandle;
		vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
		if(vtxHandle.isValid()){ vtx = *vtxHandle;}
                int nvtx = vtx.size();

                fwlite::Handle< std::vector<PileupSummaryInfo> > puInfoH;
                puInfoH.getByLabel(ev, "addPileupInfo");
                if(!puInfoH.isValid()){printf("collection PileupSummaryInfos with name addPileupInfo does not exist\n"); exit(0);}
                int ngenITpu = 0;
                for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++){
                   if(it->getBunchCrossing()==0)      { ngenITpu += it->getPU_NumInteractions(); }
                }

		//pileup weight
		float weight = 1.0;
		//double TotalWeight_plus = 1.0;
		//double TotalWeight_minus = 1.0;
		float puWeight(1.0);

		if(isMC){          
                        //std::cout << ngenITpu << std::endl;     
			puWeight          = LumiWeights->weight(ngenITpu) * PUNorm[0];
			weight            = xsecWeight*puWeight;
			//TotalWeight_plus  = PuShifters[utils::cmssw::PUUP  ]->Eval(genEv.ngenITpu) * (PUNorm[2]/PUNorm[0]);
			//TotalWeight_minus = PuShifters[utils::cmssw::PUDOWN]->Eval(genEv.ngenITpu) * (PUNorm[1]/PUNorm[0]);
		}


                std::vector<TString> chTags;
                chTags.push_back("all");

   	        mon.fillHisto("nvtx"        ,   chTags, nvtx,      weight);
		mon.fillHisto("nvtxraw"     ,   chTags, nvtx,      weight/puWeight);

	}//end of the event loop  
	printf("\n"); 

	//##############################################
	//########     SAVING HISTO TO FILE     ########
	//##############################################
	//save control plots to file
	outUrl += "/";
	outUrl += outFileUrl + ".root";
	printf("Results save in %s\n", outUrl.Data());

	//save all to the file
	TFile *ofile=TFile::Open(outUrl, "recreate");
	mon.Write();
	ofile->Close();

	if(outTxtFile)fclose(outTxtFile);

}  




//TEST


