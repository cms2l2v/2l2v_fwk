#include "UserCode/llvv_fwk/interface/EwkCorrections.h"

namespace EwkCorrections
{
	//Read correction table
	std::vector<std::vector<float>> readFile_and_loadEwkTable(TString url){
		ifstream myReadFile;
		std::vector<float> Table_line;
		std::vector<std::vector<float>> Table_EWK;
		TString name;
		TString cmssw_path;
		cmssw_path = getenv("CMSSW_BASE");
		TString path = cmssw_path+"/src/UserCode/llvv_fwk/src/";
		
		if(url.Contains("ZZ")) name = path+"EwkCorrections/ZZ_EwkCorrections.dat";
		myReadFile.open(name);
		if(!myReadFile.is_open()) cout<<"WARNING: "+name+" NOT FOUND"<<endl;
		int Start=0;
		while (!myReadFile.eof()){
			Start++;
			string output;
			myReadFile >> output;
			if(Start%5!=0) Table_line.push_back(atof(output.c_str()));
			if(Start%5==0){
				Table_line.push_back(atof(output.c_str()));
				Table_EWK.push_back(Table_line);
				Table_line.clear();
			}
		}
		myReadFile.close();
		return Table_EWK;
	}


	//Find the right correction in the file	
	std::vector<float> findCorrection(const std::vector<std::vector<float>> & Table_EWK, float sqrt_s_hat, float t_hat){
		//find the range of sqrt s hat (each 200 lines it changes)
		unsigned int j = 0;
		float best = 0.8E+04; //highest value of sqrt s hat in the table
		if( sqrt_s_hat > best) j = 39800; //in the very rare case where we have bigger s than our table (table is for 8TeV and we run at 13TeV)
		else{
			for(unsigned int i = 0 ; i < 40000 ; i = i+200){
				if(fabs(sqrt_s_hat - Table_EWK[i][0]) < best){
					best = fabs(sqrt_s_hat - Table_EWK[i][0]);
					j = i;
				}
				else break ;
			}
		}
		best = Table_EWK[j+199][1];
		if(t_hat > best) j = j+199; //in the very rare case where we have bigger t than our table
		else{
			best = 0.1E+09;
			for(unsigned int k = j ; k < j + 200 ; k++){
				if(fabs(t_hat - Table_EWK[k][1]) < best){
					best = fabs(t_hat - Table_EWK[k][1]);
					j = k;
				}
				else break ;
			}
		}
		std::vector<float> EWK_w2_vec;
		EWK_w2_vec.push_back(Table_EWK[j][2]); //ewk corrections for quark u/c
		EWK_w2_vec.push_back(Table_EWK[j][3]); //ewk corrections for quark d/s
		EWK_w2_vec.push_back(Table_EWK[j][4]); //ewk corrections for quark b
		return EWK_w2_vec ;
	}


	//The main function, will return the kfactor
	double getEwkCorrections(TString url, reco::GenParticleCollection genParticles, std::vector<std::vector<float>> Table, GenEventInfoProduct eventInfo, SmartSelectionMonitor& mon){
		double kFactor = 1.;

		reco::GenParticleCollection genIncomingQuarks;
		reco::GenParticleCollection genIncomingGluons;
		reco::GenParticleCollection genLeptons;
		reco::GenParticleCollection genNeutrinos;

		for (unsigned int i =0; i < genParticles.size(); i++){
			reco::GenParticle genParticle = genParticles[i];

			if(fabs(genParticle.pdgId()) >= 1 && fabs(genParticle.pdgId()) <= 5 && genParticle.status() == 21) genIncomingQuarks.push_back(genParticle); //status 21 : incoming particles of hardest subprocess
			if(fabs(genParticle.pdgId()) == 21 && genParticle.status() == 21) genIncomingGluons.push_back(genParticle);
			if(genParticle.status()==1 && (fabs(genParticle.pdgId())==11 || fabs(genParticle.pdgId())==13 || fabs(genParticle.pdgId())==15)) genLeptons.push_back(genParticle); //status 1 : final state
			if(genParticle.status()==1 && (fabs(genParticle.pdgId())==12 || fabs(genParticle.pdgId())==14 || fabs(genParticle.pdgId())==16)) genNeutrinos.push_back(genParticle);
		}

		std::sort(genIncomingQuarks.begin(), genIncomingQuarks.end(), utils::sort_CandidatesByPt);
		std::sort(genIncomingGluons.begin(), genIncomingGluons.end(), utils::sort_CandidatesByPt);
		std::sort(genLeptons.begin(), genLeptons.end(), utils::sort_CandidatesByPt);
		std::sort(genNeutrinos.begin(), genNeutrinos.end(), utils::sort_CandidatesByPt);

		if(!eventInfo.pdf()) return 1; //no corrections can be applied because we need x1 and x2 
		if( genLeptons.size() < 2 || genNeutrinos.size() < 2) return 1; //no corrections can be applied if we don't find our two Z's
		double x1 = eventInfo.pdf()->x.first; 
		double x2 = eventInfo.pdf()->x.second; 
		int id1=0;
		int id2=0;

		LorentzVector Z1 = genLeptons[0].p4() + genLeptons[1].p4(); //First Z : charged leptons
		LorentzVector Z2 = genNeutrinos[0].p4() + genNeutrinos[1].p4(); //Second Z : neutrinos
		LorentzVector ZZ = Z1+Z2;
		TLorentzVector ZZ_t(ZZ.X(),ZZ.Y(),ZZ.Z(),ZZ.T()); //Need TLorentzVectors for several methods (boosts)
		TLorentzVector Z1_t(Z1.X(),Z1.Y(),Z1.Z(),Z1.T());
		TLorentzVector Z2_t(Z2.X(),Z2.Y(),Z2.Z(),Z2.T());

		double s_hat = pow(ZZ.M(),2); // s_hat = center-of-mass energy of 2 Z system

		//Boost quarks and Z1
		TLorentzVector Z1_b = Z1_t;
		TLorentzVector p1_b, p2_b;
		double energy = 6500. ; //13 TeV in total
		p1_b.SetXYZT(0.,0.,x1*energy,x1*energy); //x1 = fraction of momentum taken by the particle initiating the hard process
		p2_b.SetXYZT(0.,0.,-x2*energy,x2*energy);
		Z1_b.Boost( -ZZ_t.BoostVector()); //Inverse Lorentz transformation, to get to the center-of-mass frame
		p1_b.Boost( -ZZ_t.BoostVector());
		p2_b.Boost( -ZZ_t.BoostVector());

		//Unitary vectors
		TLorentzVector Z1_b_u = Z1_b*(1/Z1_b.P()); //Normalized to 1
		TLorentzVector p1_b_u = p1_b*(1/p1_b.P());
		TLorentzVector p2_b_u = p2_b*(1/p2_b.P());

		//Effective beam axis
		TLorentzVector diff_p = p1_b_u - p2_b_u;
		TLorentzVector eff_beam_axis = diff_p*(1./diff_p.P());
		double cos_theta = eff_beam_axis.X()*Z1_b_u.X() + eff_beam_axis.Y()*Z1_b_u.Y() + eff_beam_axis.Z()*Z1_b_u.Z();

		double m_z = 91.1876; //Z bosons assumed to be on-shell
		double t_hat = m_z*m_z - 0.5*s_hat + cos_theta * sqrt( 0.25*s_hat*s_hat - m_z*m_z*s_hat );

		int quark_type = 0; //Flavour of incident quark
		if(genIncomingQuarks.size() > 0) quark_type = fabs(genIncomingQuarks[0].pdgId()); //Works unless if gg->ZZ process : it shouldn't be the case as we're using POWHEG

		std::vector<float> Correction_vec = findCorrection( Table, sqrt(s_hat), t_hat ); //Extract the corrections for the values of s and t computed

		if(quark_type==1) kFactor = 1. + Correction_vec[1]; //d
		if(quark_type==2) kFactor = 1. + Correction_vec[0]; //u
		if(quark_type==3) kFactor = 1. + Correction_vec[1]; //s as d
		if(quark_type==4) kFactor = 1. + Correction_vec[0]; //c as u
		if(quark_type==5) kFactor = 1. + Correction_vec[2]; //b

		if(sqrt(s_hat)< 2*m_z) kFactor = 1.; //Off-shell cases, not corrected to avoid non-defined values for t.

  		//Fill control histograms
		if(sqrt(s_hat) > 2*m_z){
			mon.fillHisto("Nevent_vs_ZpT", "ll_LO", Z1.Pt(), 1.);
			mon.fillHisto("Nevent_vs_ZpT", "ll_NLO", Z1.Pt(), kFactor);
			mon.fillHisto("Nevent_vs_ZpT", "vv_LO", Z2.Pt(), 1.);
			mon.fillHisto("Nevent_vs_ZpT", "vv_NLO", Z2.Pt(), kFactor);
		}
		mon.fillHisto("Nevent_vs_Mzz", "LO", ZZ.mass(), 1.);
		mon.fillHisto("Nevent_vs_Mzz", "NLO", ZZ.mass(), kFactor);

		if(genIncomingQuarks.size() == 2 && genIncomingGluons.size() ==0){
			id2 = quark_type;
			id1 = id2;
		}
		else if( genIncomingQuarks.size() == 1 && genIncomingGluons.size() ==1){
			id1 = quark_type;
			id2 = 21;
		}

		if(id1 == 2 || id1 == 4){ //u or c quark
  		mon.fillHisto("s_vs_t", "uc", t_hat, sqrt(s_hat));
  		mon.fillHisto("k_vs_s", "uc", sqrt(s_hat), kFactor);
  		mon.fillHisto("k_vs_t", "uc", t_hat, kFactor);
		}
		else if(id1 ==1 || id1 == 3){ //d or s quark
 		 	mon.fillHisto("s_vs_t", "ds", t_hat, sqrt(s_hat));
  		mon.fillHisto("k_vs_s", "ds", sqrt(s_hat), kFactor);
  		mon.fillHisto("k_vs_t", "ds", t_hat, kFactor);
		}
		else if( id1 == 5){ //b quark
		 	mon.fillHisto("s_vs_t", "b", t_hat, sqrt(s_hat));
  		mon.fillHisto("k_vs_s", "b", sqrt(s_hat), kFactor);
  		mon.fillHisto("k_vs_t", "b", t_hat, kFactor);
		}

		//qq, uu, cc, dd, ss, bb, qg, ug, cg, dg, sg, bg, other
		if( id1 >0 && id1 <9 && id1 == id2 ){ //qq case
			mon.fillHisto( "count_quarks_type", "study", 0, 1.); //qq fill
			switch ( id1 ){
				case 1:	//dd
					mon.fillHisto( "count_quarks_type", "study", 3, 1.);
					break;
				case 2:	//uu
					mon.fillHisto( "count_quarks_type", "study", 1, 1.);
					break;
				case 3:	//ss
					mon.fillHisto( "count_quarks_type", "study", 4, 1.);
					break;
				case 4:	//cc
					mon.fillHisto( "count_quarks_type", "study", 2, 1.);
					break;
				case 5:	//bb
					mon.fillHisto( "count_quarks_type", "study", 5, 1.);
					break;
				default:
					mon.fillHisto( "count_quarks_type", "study", 12, 1.);
					break;
			}
		}
		else if ( (id1 == 21 && id2 >0 && id2 < 9) || (id2 == 21 && id1 >0 && id1 < 9) ){ //qg case
			mon.fillHisto( "count_quarks_type", "study", 6, 1.); //qg fill
			switch ( min(id1, id2) ){
				case 1:	//dg
					mon.fillHisto( "count_quarks_type", "study", 9, 1.);
					break;
				case 2:	//ug
					mon.fillHisto( "count_quarks_type", "study", 7, 1.);
					break;
				case 3:	//sg
					mon.fillHisto( "count_quarks_type", "study", 10, 1.);
					break;
				case 4:	//cg
					mon.fillHisto( "count_quarks_type", "study", 8, 1.);
					break;
				case 5:	//bg
					mon.fillHisto( "count_quarks_type", "study", 11, 1.);
					break;
				default:
					mon.fillHisto( "count_quarks_type", "study", 12, 1.);
					break;
			}

		}
		else{
			mon.fillHisto( "count_quarks_type", "study", 12, 1.);
		}

		return kFactor;
	}

}
