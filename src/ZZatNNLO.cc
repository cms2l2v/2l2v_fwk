#include "UserCode/llvv_fwk/interface/ZZatNNLO.h"

namespace ZZatNNLO
{
	//Read correction table
	std::vector<std::vector<float>> readFile_and_loadTable(TString dtag){
		ifstream myReadFile;
		std::vector<float> Table_line;
		std::vector<std::vector<float>> Table;
		TString name;
		TString cmssw_path;
		cmssw_path = getenv("CMSSW_BASE");
		TString path = cmssw_path+"/src/UserCode/llvv_fwk/src/";
		
		if(dtag.Contains("ZZ")) name = path+"Corrections/ZZ_NNLOQCD.dat";
		myReadFile.open(name);
		if(!myReadFile.is_open()) cout<<"WARNING: "+name+" NOT FOUND"<<endl;
		int Start=0;
		while (!myReadFile.eof()){
			Start++;
			string output;
			myReadFile >> output;
			if(Start%2!=0) Table_line.push_back(atof(output.c_str()));
			if(Start%2==0){
				Table_line.push_back(atof(output.c_str()));
				Table.push_back(Table_line);
				Table_line.clear();
			}
		}
		myReadFile.close();
		return Table;
	}


	//Find the right correction in the file	
	double findCorrection(const std::vector<std::vector<float>> & Table, double dphi){
		
		double kFactor = 1.;
		if(dphi >= 3.1416) return 1.;
		for( int i = 0 ; i <= 29; i++){

			if(dphi == i*0.1){
				return Table[i][1];
			}
			else if(dphi < i*0.1){
				return Table[i-1][1];
			}
		}
		if(dphi >= 2.9) kFactor = Table[29][1];
		
		return kFactor ;
	}


	//The main function, will return the kfactor
	double getNNLOCorrections(TString dtag, const reco::GenParticleCollection & genParticles, const std::vector<std::vector<float>> & Table){
		double kFactor = 1.;

		reco::GenParticleCollection genLeptons;
		reco::GenParticleCollection genNeutrinos;

		for (unsigned int i =0; i < genParticles.size(); i++){
			reco::GenParticle genParticle = genParticles[i];

			if(genParticle.status()==1 && (fabs(genParticle.pdgId())==11 || fabs(genParticle.pdgId())==13 || fabs(genParticle.pdgId())==15)) genLeptons.push_back(genParticle); //status 1 : final state
			if(genParticle.status()==1 && (fabs(genParticle.pdgId())==12 || fabs(genParticle.pdgId())==14 || fabs(genParticle.pdgId())==16)) genNeutrinos.push_back(genParticle);
		}

		std::sort(genLeptons.begin(), genLeptons.end(), utils::sort_CandidatesByPt);
		std::sort(genNeutrinos.begin(), genNeutrinos.end(), utils::sort_CandidatesByPt);

		if( genLeptons.size() < 2 || genNeutrinos.size() < 2) return 1; //no corrections can be applied if we don't find our two Z's

		LorentzVector Z1 = genLeptons[0].p4() + genLeptons[1].p4(); //First Z : charged leptons
		LorentzVector Z2 = genNeutrinos[0].p4() + genNeutrinos[1].p4(); //Second Z : neutrinos

    double dphi=fabs(deltaPhi(Z1.phi(),Z2.phi()));


		kFactor = findCorrection( Table, dphi ); //Extract the corrections for the values of s and t computed

		return kFactor;
	}

}
