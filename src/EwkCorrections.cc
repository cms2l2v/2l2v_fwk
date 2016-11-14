#include "UserCode/llvv_fwk/interface/EwkCorrections.h"

namespace EwkCorrections
{
	//Read correction table
	std::vector<std::vector<float>> readFile_and_loadEwkTable(TString dtag){
		ifstream myReadFile;
		std::vector<float> Table_line;
		std::vector<std::vector<float>> Table_EWK;
		TString name;
		TString cmssw_path;
		cmssw_path = getenv("CMSSW_BASE");
		TString path = cmssw_path+"/src/UserCode/llvv_fwk/src/";
		
		if(dtag.Contains("ZZ")) name = path+"Corrections/ZZ_EwkCorrections.dat";
		if(dtag.Contains("WZ")) name = path+"Corrections/WZ_EwkCorrections.dat";
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
	double getEwkCorrections(TString dtag, const reco::GenParticleCollection & genParticles, const std::vector<std::vector<float>> & Table, const GenEventInfoProduct & eventInfo, double & ewkCorrections_error){
		double kFactor = 1.;

		bool isWZ(dtag.Contains("WZ"));
		bool isZZ(dtag.Contains("ZZ"));

		reco::GenParticleCollection genIncomingQuarks;
		reco::GenParticleCollection genIncomingGluons;
		reco::GenParticleCollection genWBosons;
		reco::GenParticleCollection genZBosons;
		reco::GenParticleCollection genLeptonsFromZ;
		reco::GenParticleCollection genLeptonsFromW;
		reco::GenParticleCollection genNeutrinosFromZ;
		reco::GenParticleCollection genNeutrinosFromW;

		for (unsigned int i =0; i < genParticles.size(); i++){
			reco::GenParticle genParticle = genParticles[i];

			if(fabs(genParticle.pdgId()) >= 1 && fabs(genParticle.pdgId()) <= 5 && genParticle.status() == 21) genIncomingQuarks.push_back(genParticle); //status 21 : incoming particles of hardest subprocess
			if(fabs(genParticle.pdgId()) == 21 && genParticle.status() == 21) genIncomingGluons.push_back(genParticle);
			if(fabs(genParticle.pdgId()) == 24 && genParticle.status() == 22) genWBosons.push_back(genParticle);
			if(fabs(genParticle.pdgId()) == 23 && genParticle.status() == 22) genZBosons.push_back(genParticle);
			if(fabs(genParticle.pdgId())==11 || fabs(genParticle.pdgId())==13 || fabs(genParticle.pdgId())==15){
				if(fabs(genParticle.mother()->pdgId())==23 && genParticle.mother()->status()==62) genLeptonsFromZ.push_back(genParticle); //lepton originating directly from Z boson. Status 62 seems to correspond to status just before the decay of the boson.
				if(fabs(genParticle.mother()->pdgId())==24 && genParticle.mother()->status()==62) genLeptonsFromW.push_back(genParticle); //lepton originating directly from W boson
				}
			if(fabs(genParticle.pdgId())==12 || fabs(genParticle.pdgId())==14 || fabs(genParticle.pdgId())==16){
				if(fabs(genParticle.mother()->pdgId())==23 && genParticle.mother()->status()==62) genNeutrinosFromZ.push_back(genParticle); //neutrino originating directly from Z boson
				if(fabs(genParticle.mother()->pdgId())==24 && genParticle.mother()->status()==62) genNeutrinosFromW.push_back(genParticle); //neutrino originating directly from W boson
				}
		}

		std::sort(genIncomingQuarks.begin(), genIncomingQuarks.end(), utils::sort_CandidatesByPt);
		std::sort(genIncomingGluons.begin(), genIncomingGluons.end(), utils::sort_CandidatesByPt);
		std::sort(genWBosons.begin(), genWBosons.end(), utils::sort_CandidatesByPt);
		std::sort(genZBosons.begin(), genZBosons.end(), utils::sort_CandidatesByPt);
		std::sort(genLeptonsFromW.begin(), genLeptonsFromW.end(), utils::sort_CandidatesByPt);
		std::sort(genNeutrinosFromW.begin(), genNeutrinosFromW.end(), utils::sort_CandidatesByPt);
		std::sort(genLeptonsFromZ.begin(), genLeptonsFromZ.end(), utils::sort_CandidatesByPt);
		std::sort(genNeutrinosFromZ.begin(), genNeutrinosFromZ.end(), utils::sort_CandidatesByPt);

		if(!eventInfo.pdf()) return 1; //no corrections can be applied because we need x1 and x2 
		if( genZBosons.size() < 2 && isZZ) return 1; //no corrections can be applied if we don't find our two Z's
		if( (genZBosons.size() < 1 || genWBosons.size() < 1) && isWZ) return 1; //no corrections can be applied if we don't find at least one W and one Z
		double x1 = eventInfo.pdf()->x.first; 
		double x2 = eventInfo.pdf()->x.second; 

		LorentzVector V1;
		LorentzVector V2;

		if (isZZ){
			V1 = genLeptonsFromZ[0].p4() + genLeptonsFromZ[1].p4(); //First Z : charged leptons
			V2 = genNeutrinosFromZ[0].p4() + genNeutrinosFromZ[1].p4(); //Second Z : neutrinos
		}
		if (isWZ){
			V1 = genLeptonsFromZ[0].p4() + genLeptonsFromZ[1].p4(); //Z decaying in 2 leptons
			V2 = genLeptonsFromW[0].p4() + genNeutrinosFromW[0].p4(); //W decaying in 1 lepton and 1 neutrino
		}

		LorentzVector VV = V1+V2;
		TLorentzVector VV_t(VV.X(),VV.Y(),VV.Z(),VV.T()); //Need TLorentzVectors for several methods (boosts)
		TLorentzVector V1_t(V1.X(),V1.Y(),V1.Z(),V1.T());
		TLorentzVector V2_t(V2.X(),V2.Y(),V2.Z(),V2.T());

		double s_hat = pow(VV.M(),2); // s_hat = center-of-mass energy of 2 vector boson system

		//Boost quarks and V1
		TLorentzVector V1_b = V1_t;
		TLorentzVector p1_b, p2_b;
		double energy = 6500. ; //13 TeV in total
		p1_b.SetXYZT(0.,0.,x1*energy,x1*energy); //x1 = fraction of momentum taken by the particle initiating the hard process
		p2_b.SetXYZT(0.,0.,-x2*energy,x2*energy);
		V1_b.Boost( -VV_t.BoostVector()); //Inverse Lorentz transformation, to get to the center-of-mass frame
		p1_b.Boost( -VV_t.BoostVector());
		p2_b.Boost( -VV_t.BoostVector());

		//Unitary vectors
		TLorentzVector V1_b_u = V1_b*(1/V1_b.P()); //Normalized to 1
		TLorentzVector p1_b_u = p1_b*(1/p1_b.P());
		TLorentzVector p2_b_u = p2_b*(1/p2_b.P());

		//Effective beam axis
		TLorentzVector diff_p = p1_b_u - p2_b_u;
		TLorentzVector eff_beam_axis = diff_p*(1./diff_p.P());
		double cos_theta = eff_beam_axis.X()*V1_b_u.X() + eff_beam_axis.Y()*V1_b_u.Y() + eff_beam_axis.Z()*V1_b_u.Z();

		double m_z = 91.1876; //Z bosons assumed to be on-shell
		double m_w = 80.385;
		double t_hat = 0.;

		if(isZZ) t_hat = m_z*m_z - 0.5*s_hat + cos_theta * sqrt( 0.25*s_hat*s_hat - m_z*m_z*s_hat );
		if(isWZ){
			double b = 1./2./sqrt(s_hat) * sqrt(pow(s_hat-m_z*m_z-m_w*m_w,2) - 4*m_w*m_w*m_z*m_z);
			double a = sqrt(b*b + m_z*m_z);
			t_hat = m_z*m_z - sqrt(s_hat) * (a - b * cos_theta); //awful calculation, needed to put ourselves to the center-of-mass frame with the 2 particles having a different mass !
			}

		int quark_type = 0; //Flavour of incident quark
		if(genIncomingQuarks.size() > 0) quark_type = fabs(genIncomingQuarks[0].pdgId()); //Works unless if gg->ZZ process : it shouldn't be the case as we're using POWHEG

		std::vector<float> Correction_vec = findCorrection( Table, sqrt(s_hat), t_hat ); //Extract the corrections for the values of s and t computed

		if(quark_type==1) kFactor = 1. + Correction_vec[1]; //d
		if(quark_type==2) kFactor = 1. + Correction_vec[0]; //u
		if(quark_type==3) kFactor = 1. + Correction_vec[1]; //s as d
		if(quark_type==4) kFactor = 1. + Correction_vec[0]; //c as u
		if(quark_type==5) kFactor = 1. + Correction_vec[2]; //b  //Notice that the quark types are irrelevant for the case of WZ (same numbers in the last 3 columns).

		if(sqrt(s_hat)< 2*m_z && isZZ) kFactor = 1.; //Off-shell cases, not corrected to avoid non-defined values for t.
		if(sqrt(s_hat)< m_z + m_w && isWZ) kFactor = 1.;

		//Computing the associated error:
		//Warning, several methods could be used.
		//In Run 1, CMS used (kFactor-1)*(kFactor_QCD -1) for all rho
		//And ATLAS used : 0 for rho < 0.3 and 1 for rho >0.3
		//
		//Here is an implementation that is using a mix of the two. It may change in the future (but the change won't be critical)
		double kFactor_QCD = 1.;
		if(isZZ) kFactor_QCD = 15.99/9.89; //From arXiv1105.0020
		if(isWZ && genWBosons[0].pdgId() > 0) kFactor_QCD = 28.55/15.51; //for W+Z
		if(isWZ && genWBosons[0].pdgId() < 0) kFactor_QCD = 18.19/9.53; //for W-Z

		
		//Definition of rho
		double rho = 0.;
		if(isZZ){
			double rho = (genLeptonsFromZ[0].p4() + genLeptonsFromZ[1].p4() + genNeutrinosFromZ[0].p4() + genNeutrinosFromZ[1].p4()).pt();
			rho = rho/(genLeptonsFromZ[0].pt() + genLeptonsFromZ[1].pt() + genNeutrinosFromZ[0].pt() + genNeutrinosFromZ[1].pt());
		}
		if(isWZ){
			double rho = (genLeptonsFromZ[0].p4() + genLeptonsFromZ[1].p4() + genLeptonsFromW[0].p4() + genNeutrinosFromW[0].p4()).pt();
			rho = rho/(genLeptonsFromZ[0].pt() + genLeptonsFromZ[1].pt() + genLeptonsFromW[0].pt() + genNeutrinosFromW[0].pt());
		}

		if(rho<0.3) ewkCorrections_error = fabs((kFactor-1)*(kFactor_QCD -1));
		else ewkCorrections_error = fabs(1-kFactor);

		//At this point, we have the relative error on the delta_ewk ( = k_ewk -1 )
		//Let's - instead - return the absolute error on k: we do delta_ewk* the_relative_errir_on_it. This gives absolute error on delta, and so on k
		ewkCorrections_error = fabs(ewkCorrections_error*kFactor);

		//For WZ, contribution from gamma-induced processes
		double gamma_induced_uncertainty = 0.;
		if(isWZ){
			//We multiply the kFactor from virtual processes (computed above) with the one from gamma-induced processes (based on a fit form a separate study)
			if(genWBosons[0].pdgId() > 0) kFactor = kFactor*(1 + 0.00559445 - 5.17082e-6 * sqrt(s_hat) + 3.63331e-8 * s_hat); //W+Z
			if(genWBosons[0].pdgId() < 0) kFactor = kFactor*(1 + 0.00174737 + 1.70668e-5 * sqrt(s_hat) + 2.26398e-8 * s_hat); //W-Z

			//Uncertainty due to WZ gamma-induced contribution is set to 0 with a very good approximation (less than 0.1%) when using the LUXqed photon PDF (instead of NNPDF23)
			if(genWBosons[0].pdgId() > 0) gamma_induced_uncertainty = 0.; //0.00286804 - 8.4624e-6 * sqrt(s_hat) + 3.90611e-8 * s_hat; //W+Z
			if(genWBosons[0].pdgId() < 0) gamma_induced_uncertainty = 0.; //0.00417376 - 1.51319e-5 * sqrt(s_hat) + 5.68576e-8 * s_hat; //W-Z
			ewkCorrections_error = sqrt(pow(ewkCorrections_error,2) + pow(gamma_induced_uncertainty,2));
		}
		
		return kFactor;
	}

}
