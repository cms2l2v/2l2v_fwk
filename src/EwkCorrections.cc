#include "UserCode/llvv_fwk/interface/EwkCorrections.h"
#include <stdlib.h>

namespace EwkCorrections
{
	//Let's first define some functions

	//Read correction table
	std::vector<std::vector<float>> readFile_and_loadEwkTable(TString url){
		ifstream myReadFile;
		std::vector<float> Table_line;
		std::vector<std::vector<float>> Table_EWK;
		TString name;
		TString cmssw_path;
		cmssw_path = getenv("CMSSW_BASE");
		TString path = cmssw_path+"/src/UserCode/llvv_fwk/src/";
		if(url.Contains("ZZ")){
			name = path+"EwkCorrections/ZZ_EwkCorrections.dat";
		}
		//else if (url.Contains("WZ")){ //Ready for WZ but don't have the table yet
		//	name = path+"EwkCorrections/WZ_EwkCorrections.dat";
		//}
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
		int j = 0;
		float best = 0.8E+04; //highest value of sqrt s hat in the table
		if( sqrt_s_hat > best) j = 39800; //in the very rare case where we have bigger s than our table (table is for 8TeV and we run at 13TeV)
		else{
			for(int i = 0 ; i < 40000 ; i = i+200){
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
			for(int k = j ; k < j + 200 ; k++){
				if(fabs(t_hat - Table_EWK[k][1]) < fabs(best)){
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
		cout<< "getting ewk corrections..." << endl;
		double kFactor = 1.;

		//create leptons and photons and neutrinos
		reco::GenParticleCollection genIncomingQuarks;
		reco::GenParticleCollection genIncomingGluons;
		reco::GenParticleCollection genLeptons;
		//reco::GenParticleCollection genPhotons;
		reco::GenParticleCollection genNeutrinos;
		reco::GenParticleCollection genIncomingProtons;

		for (unsigned int i =0; i < genParticles.size(); i++){
			reco::GenParticle genParticle = genParticles[i];

			if(genParticle.status() == 2) continue; //drop the shower part of the history

			if(genParticle.pdgId() == 2212 && genParticle.status() == 4) genIncomingProtons.push_back(genParticle); //4 : incoming beam particle
			if(fabs(genParticle.pdgId()) >= 1 && fabs(genParticle.pdgId()) <= 5 && genParticle.status() == 21) genIncomingQuarks.push_back(genParticle); //21 : incoming particles of hardest subprocess
			if(fabs(genParticle.pdgId()) == 21 && genParticle.status() == 21) genIncomingGluons.push_back(genParticle);


			//leptons
			if(genParticle.status()==1 && (fabs(genParticle.pdgId())==11 || fabs(genParticle.pdgId())==13 || fabs(genParticle.pdgId())==15)){
				genLeptons.push_back(genParticle);
			}
			//photons
			//else if(genParticle.status()==1 && genParticle.pdgId()==22){
			//  genPhotons.push_back(genParticle);
			//}
			//neutrinos
			else if(genParticle.status()==1 && (fabs(genParticle.pdgId())==12 || fabs(genParticle.pdgId())==14 || fabs(genParticle.pdgId())==16)){
				genNeutrinos.push_back(genParticle);
			}


		}

		std::sort(genIncomingQuarks.begin(), genIncomingQuarks.end(), utils::sort_CandidatesByPt);
		std::sort(genIncomingGluons.begin(), genIncomingGluons.end(), utils::sort_CandidatesByPt);
		std::sort(genLeptons.begin(), genLeptons.end(), utils::sort_CandidatesByPt);
		//std::sort(genPhotons.begin(), genPhotons.end(), utils::sort_CandidatesByPt);
		std::sort(genNeutrinos.begin(), genNeutrinos.end(), utils::sort_CandidatesByPt);
		std::sort(genIncomingProtons.begin(), genIncomingProtons.end(), utils::sort_CandidatesByPt);

		//This addition is from GenEventInfo. It contains a lot of info on PDFs, x1, x2, id1, id2... Maybe we could (should?) use this instead of our previous way to look at things? At least we should ask for a correspondance (same id...)

		if(!eventInfo.pdf()) return 1; //no corrections can be applied because we need x1 and x2 
		cout << "Q scale : " << eventInfo.pdf()->scalePDF << endl; 
		double x1 = eventInfo.pdf()->x.first; 
		double x2 = eventInfo.pdf()->x.second; 
		int id1 = fabs(eventInfo.pdf()->id.first); 
		int id2 = fabs(eventInfo.pdf()->id.second); 





		//cout << "test table : en (6,2) on a " << Table[6][2] << endl;
		//Here : generated particles are put in vectors

		//All what follows is valable only for ZZ->2l2nu. Another code is needed for WZ.

		//We don't compute rho for the moment (necessary after ?)

		LorentzVector Z1 = genLeptons[0].p4() + genLeptons[1].p4(); //First Z : charged leptons
		LorentzVector Z2 = genNeutrinos[0].p4() + genNeutrinos[1].p4(); //Second Z : neutrinos
		LorentzVector ZZ = Z1+Z2;
		//Mainpulation to have TLorentzVectors instead of LorentzVectors
		TLorentzVector ZZ_t;
		TLorentzVector Z1_t;
		TLorentzVector Z2_t;
		ZZ_t.SetXYZT(ZZ.X(),ZZ.Y(),ZZ.Z(),ZZ.T());
		Z1_t.SetXYZT(Z1.X(),Z1.Y(),Z1.Z(),Z1.T());
		Z2_t.SetXYZT(Z2.X(),Z2.Y(),Z2.Z(),Z2.T());

		double s_hat = pow(ZZ.M(),2); // s_hat = center-of-mass energy of 2 Z system

		//In order to determine t_hat, we have to compute theta : angle between considered Z boson and the direction of the incident quarks in the rest frame of the 2 Z.
		//At the end, we'll have cos(theta) = (p_q1,b - p_q2,b)/(|p_q1,b - p_q2,b|) . p_Z1,b  , where p_qi,b is the *unitary* momentum *vector* of the i quark after the Lorentz boost.
		//(cf Nicolas's Master Thesis, p.15)
		//Let's first boost the quarks 1 and 2, as well as the Z1.
		TLorentzVector Z1_b = Z1_t;
		TLorentzVector p1_b, p2_b;
		double energy = 6500. ; //13 TeV in total
		p1_b.SetXYZT(0.,0.,x1*energy,x1*energy); //x1 = fraction of momentum taken by the particle initiating the hard process
		p2_b.SetXYZT(0.,0.,-x2*energy,x2*energy); //goes the other way !
		Z1_b.Boost( -ZZ_t.BoostVector()); //Inverse Lorentz transformation, to get to the center-of-mass frame
		p1_b.Boost( -ZZ_t.BoostVector());
		p2_b.Boost( -ZZ_t.BoostVector());
		//Unitary vectors
		TLorentzVector Z1_b_u = Z1_b*(1/Z1_b.P()); // Still LorentzVectors, but we look only at spatial components. The vector is normalized now.
		TLorentzVector p1_b_u = p1_b*(1/p1_b.P());
		TLorentzVector p2_b_u = p2_b*(1/p2_b.P());
		//Effective beam axis
		TLorentzVector diff_p = p1_b_u - p2_b_u;
		TLorentzVector eff_beam_axis = diff_p*(1./diff_p.P());
		//And finally the expression for cos(theta). Does anyone know a better way of doing things ?
		double cos_theta = eff_beam_axis.X()*Z1_b_u.X() + eff_beam_axis.Y()*Z1_b_u.Y() + eff_beam_axis.Z()*Z1_b_u.Z();
		//Now, compute t_hat according to the formula above.
		//double t_hat = pow(Z1.M(),2) - 0.5*s_hat + cos_theta * sqrt(0.25*s_hat*s_hat - pow(Z1.M(),2)*s_hat); //old version ; actually we assume the Z to be on-shell.
		double m_z = 91.1876;
		double t_hat = m_z*m_z - 0.5*s_hat + cos_theta * sqrt( 0.25*s_hat*s_hat - m_z*m_z*s_hat );

		//Find the quark type
		int quark_type = 0;
		if(genIncomingQuarks.size() > 0) quark_type = fabs(genIncomingQuarks[0].pdgId()); //works unless if gg->ZZ process : it shouldn't be the case as we're using POWHEG

		//Extract the corrections for the values of s and t just computed
		std::vector<float> Correction_vec = findCorrection( Table, sqrt(s_hat), t_hat );

		//Final correction factor
		if(quark_type==1) kFactor = 1. + Correction_vec[1]; //d
		if(quark_type==2) kFactor = 1. + Correction_vec[0]; //u
		if(quark_type==3) kFactor = 1. + Correction_vec[1]; //s as d
		if(quark_type==4) kFactor = 1. + Correction_vec[0]; //c as u
		if(quark_type==5) kFactor = 1. + Correction_vec[2]; //b
		//else, it stays at 1.

		cout<< "kFactor found !" << endl;

		cout << "x1 = 	" << x1 << endl;
		cout << "x2 = 	" << x2 << endl;
		cout << "m_Z = 	" << Z1.M() <<  endl;
		cout << "cos_theta = 	" << cos_theta << endl;
		cout << "den = 	" << diff_p.P() << endl;
		cout << "nb of quarks :		" << genIncomingQuarks.size() << endl;
		cout << "quark_type :	 " << quark_type << endl;
		cout << "id1 :	 " << id1 << endl;
		cout << "id2 :	 " << id2 << endl;
		cout << "sqrt(s_hat) =	 " << sqrt(s_hat) << endl;
		cout << "t_hat =	 " << t_hat << endl;
		cout << "kFactor =	 " << kFactor << endl;
		cout << endl;

  	//Fill the histograms of control
		std::vector <TString > quarks_type = {"uORc", "dORs", "b"};
  	mon.fillHisto("s_vs_t", quarks_type, t_hat, sqrt(s_hat));
  	mon.fillHisto("k_vs_s", quarks_type, sqrt(s_hat), kFactor);
  	mon.fillHisto("k_vs_t", quarks_type, t_hat, kFactor);

		mon.fillHisto("Nevent_vs_ZpT", "ll_LO", Z1.Pt(), 1.);
		mon.fillHisto("Nevent_vs_ZpT", "ll_NLO", Z1.Pt(), kFactor);
		mon.fillHisto("Nevent_vs_ZpT", "vv_LO", Z2.Pt(), 1.);
		mon.fillHisto("Nevent_vs_ZpT", "vv_NLO", Z2.Pt(), kFactor);

		mon.fillHisto("Nevent_vs_Mzz", "LO", ZZ.mass(), 1.);
		mon.fillHisto("Nevent_vs_Mzz", "NLO", ZZ.mass(), kFactor);

		//qq, uu, cc, dd, ss, bb, qg, ug, cg, dg, sg, bg, other
		if( id1 >0 && id1 <9 && id1 == id2 ){ //qq case (same flavor!)
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

	//
	//
	//  	 			if(nombre_leptons < 4) continue;
	//  	 			float gen_mz = 91.1876;
	//  	 			int gen_flag = 0; // 1 : 4 leptons de même saveur
	//
	//		 	 //Calcul de rho. Cela n'engage à rien de le calculer ; simplement, par la suite, on verra si on l'utilise ou pas...
	//  	 	 float num = sqrt(gen_ZZ.Px()*gen_ZZ.Px()+gen_ZZ.Py()*gen_ZZ.Py());
	//  	 	 float den = sqrt(gen_l1.Px()*gen_l1.Px()+gen_l1.Py()*gen_l1.Py()) + sqrt(gen_l2.Px()*gen_l2.Px()+gen_l2.Py()*gen_l2.Py()) + sqrt(gen_l3.Px()*gen_l3.Px()+gen_l3.Py()*gen_l3.Py()) + sqrt(gen_l4.Px()*gen_l4.Px()+gen_l4.Py()*gen_l4.Py());
	//  	 	 float gen_rho = num/den;
	//  	 	 //Si on a choisi de couper sur rho, cette coupure s'applique ici ! : (QCD veto)
	//  	 	 //En fait, il ne faut pas appliquer cette correction au niveau généré mais reconstruit. Voir plus loin dans la suite de l'analyse : cette coupure devient une coupure de sélection comme les autres.
	//  	 	 //if(include_rho_corrections && gen_rho > 0.3) continue;
	//  	 	 rho_genere->Fill(gen_rho);
	//
	//  	 	 //Here we extract the correction
	//  	 	 //gen_ZZ
	//  	 	 //Using the table
	//  	 	 //p1, p2, p3 and p4
	//
	//  	 	 float M_12 = 91.1876, M_22 = 91.1876;
	//  	 	 float Energy = 4000;
	//  	 	 TLorentzVector p1; p1.SetPxPyPzE(0.,0.,Energy * x1,Energy * x1);
	//  	 	 TLorentzVector p2; p2.SetPxPyPzE(0.,0.,-Energy * x2,Energy * x2);
	//  	 	 if( fabs(x1)==fabs(x2) )  cout<<"WARNING: x1=x2 "<<x1<<" "<<x2<<" with "<<id1<<" "<<id2<<endl;
	//  	 	 //S-HAT
	//  	 	 float s_hat = pow(gen_ZZ.E(),2)-pow(gen_ZZ.Px(),2)-pow(gen_ZZ.Py(),2)-pow(gen_ZZ.Pz(),2); // ScalarProd(gen_ZZ) (p1+p2)^2 = (p3+p4)^2 ~ +2*p1*p2 
	//  	 	 //T_HAT
	//  	 	 float t_hat2 = /*TypeFirst == 1 ?*/ /* p1.M()*p1.M() + M_12*M_12 - 2*( p1.E()*gen_Z1.E() - p1.Px()*gen_Z1.Px() - p1.Py()*gen_Z1.Py() - p1.Pz()*gen_Z1.Pz() ) :*/
	//     	 p1.M()*p1.M() + M_22*M_22 - 2*( p1.E()*gen_Z2.E() - p1.Px()*gen_Z2.Px() - p1.Py()*gen_Z2.Py() - p1.Pz()*gen_Z2.Pz() ); //T_HAT LO
	//  	 	 float la1 = sqrt( pow(s_hat,2) );            //la = sqrt( pow(a,2)+pow(b,2)+pow(c,2)-2*(a*b+a*c+b*c) ); 
	//  	 	 float la2 = sqrt( pow(s_hat,2) + pow(M_12,4) + pow(M_22,4) - 2*(s_hat*M_12*M_12 + s_hat*M_22*M_22 + M_12*M_12*M_22*M_22) );
	//  	 	 //  Boost: boost ext. momenta in CM frame of gen_Z1,gen_Z2
	//  	 	 TLorentzVector gen_Z1_b = gen_Z1, gen_Z2_b = gen_Z2;
	//  	 	 TLorentzVector p1_b = p1, p2_b = p2;
	//  	 	 gen_Z1_b.Boost( -gen_ZZ.BoostVector() );
	//  	 	 gen_Z2_b.Boost( -gen_ZZ.BoostVector() );
	//  	 	 p1_b.Boost( -gen_ZZ.BoostVector() );
	//  	 	 p2_b.Boost( -gen_ZZ.BoostVector() );
	//  	 	 //  Uni-vector
	//  	 	 TLorentzVector ee1 = p1_b*(1./sqrt( pow(p1_b.X(),2)+pow(p1_b.Y(),2)+pow(p1_b.Z(),2) ));
	//  	 	 TLorentzVector ee2 = p2_b*(1./sqrt( pow(p2_b.X(),2)+pow(p2_b.Y(),2)+pow(p2_b.Z(),2) ));
	//  	 	 TLorentzVector z1  = gen_Z1_b*(1./sqrt( pow(gen_Z1_b.X(),2)+pow(gen_Z1_b.Y(),2)+pow(gen_Z1_b.Z(),2) ));
	//  	 	 TLorentzVector z2  = gen_Z2_b*(1./sqrt( pow(gen_Z2_b.X(),2)+pow(gen_Z2_b.Y(),2)+pow(gen_Z2_b.Z(),2) ));
	//  	 	 //  "effective" beam axis
	//  	 	 float abse = sqrt( pow(ee1.X()-ee2.X(),2) + pow(ee1.Y()-ee2.Y(),2) + pow(ee1.Z()-ee2.Z(),2) );
	//  	 	 TLorentzVector ee = (ee1-ee2) * (1. / abse);
	//  	 	 //  "effective" scattering angle
	//  	 	 float costh = ee.X()*z1.X()+ee.Y()*z1.Y()+ee.Z()*z1.Z();
	//  	 	 //  final T_HAT
	//  	 	 float t_hat= M_12*M_12 - (1./2.)*(s_hat+M_12*M_12-M_22*M_22) + (1/(2*s_hat))*la1*la2*costh; //Mz-1/2*s+1/(2*s)*la(s,0,0)*la(s,mZ,mZ)*costh  or: (p1-p3)^2 = (p2-p4)^2 ~ -2*p1*p3
	//  	 	 //Quark Type
	//  	 	 int quarkType = -1.;
	//  	 	 if(      fabs(id1)==2 && fabs(id2)==2 ) quarkType=0; //delta_uub
	//  	 	 else if( fabs(id1)==1 && fabs(id2)==1 ) quarkType=1; //delta_ddb
	//  	 	 else if( fabs(id1)==5 && fabs(id2)==5 ) quarkType=2; //delta_bbb
	//  	 	 else if( fabs(id1)==4 && fabs(id2)==4 ) quarkType=0; // cc as delta_buu
	//  	 	 else if( fabs(id1)==3 && fabs(id2)==3 ) quarkType=1; // ss as delta_bdd
	//
	//  	 	 //Pour traiter les gluons !
	//  	 	 //Cas à 2 gluons (normalement, il n'y en a pas... mais sait-on jamais ?)
	//  	 	 if(id1 == 0 && id2 == 0){
	//     	 int nombredeZ = 0;
	//     	 for(int ii = 0; ii < mcn; ii++){
	//     	 if(mc_id[ii] == 23) nombredeZ ++;
	//     	 if(nombredeZ == 2 && mc_id[ii] != 21){
	//     	 if(fabs(mc_id[ii]) == 1) quarkType = 1;
	//     	 if(fabs(mc_id[ii]) == 2) quarkType = 0;
	//     	 if(fabs(mc_id[ii]) == 3) quarkType = 1;
	//     	 if(fabs(mc_id[ii]) == 4) quarkType = 0;
	//     	 if(fabs(mc_id[ii]) == 5) quarkType = 2;
	//     	 }
	//     	 if(quarkType != -1) break;
	//     	 }
	//  	 	 }
	//  	 	 //Cas à 1 gluon : dans ce cas, il suffit de regarder le type de l'autre quark.
	//  	 	 if(id1 == 0 && id2 != 0){
	//     	 if(fabs(id2) == 1) quarkType = 1;
	//     	 if(fabs(id2) == 2) quarkType = 0;
	//     	 if(fabs(id2) == 3) quarkType = 1;
	//     	 if(fabs(id2) == 4) quarkType = 0;
	//     	 if(fabs(id2) == 5) quarkType = 2;
	//  	 	 }
	//  	 	 if(id2 == 0 && id1 != 0){
	//     	 if(fabs(id1) == 1) quarkType = 1;
	//     	 if(fabs(id1) == 2) quarkType = 0;
	//     	 if(fabs(id1) == 3) quarkType = 1;
	//     	 if(fabs(id1) == 4) quarkType = 0;
	//     	 if(fabs(id1) == 5) quarkType = 2;
	//  	 	 }
	//
	//  	 	 //Extraction des corrections
	//  	 	 float sqrt_s_hat = sqrt(s_hat);
	//  	 	 std::vector<float> EWK_w2_vec = findCorrection( Table_EWK, sqrt_s_hat, t_hat );
	//  	 	 float EWK_w2 = 1. + EWK_w2_vec[quarkType];
	//  	 	 kFactor = EWK_w2;
	//		 	 }
	//		 	 if(!simulation_powheg) kFactor = 1.;
	//		 	 //Fin de l'extraction du facteur de correction. Maintenant, on peut passer tranquillement au reste de l'analyse sans être davantage perturbé par ce facteur.
	//
	//

}
