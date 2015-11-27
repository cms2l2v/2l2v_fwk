#include "UserCode/llvv_fwk/interface/EwkCorrections.h"

namespace EwkCorrections
{


	//Dans ma fonction main faudra ajouter une recherche d'url dans le code principal pour tomber sur les bonnes corrections en fonction des simus
	//La structure sera la suivante :
	//Dans le code central :
	// - il y aura un bool isMC avant la fonction
	// - elle prendra le nom du fichier sur lequel on tourne en argument
	// - elle prendra la collection de genParticules
	// - elle retournera le k-factor
	//Dans le src :
	// 1) Initialisation - Regarder le nom et choisir sur quelles corrections tourner
	// 2)




	//Let's first define some functions

	//Read correction table
	std::vector<std::vector<float>> readFile_and_loadEwkTable(TString url){
		ifstream myReadFile;
		std::vector<float> Table_line;
  	std::vector<std::vector<float>> Table_EWK;
		TString name;
		TString path = "$CMSSW_BASE/src/UserCode/llvv_fwk/src/";
		if(url.Contains("ZZ")){
		 	name = path+"EwkCorrections/ZZ_EwkCorrections.dat";
    }
	  else if (url.Contains("WZ")){ //Ready for WZ but don't have the table yet
	  	name = path+"EwkCorrections/WZ_EwkCorrections.dat";
	  }
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
 		cout << "File read" << endl;
 		cout << "test table : en (6,2) on a " << Table_EWK[6][2] << endl;
 		return Table_EWK;
	}






	//Find the right correction in the file	

	std::vector<float> findCorrection(const std::vector<std::vector<float>> & Table_EWK, float sqrt_s_hat, float t_hat){
		//find the range of sqrt s har (each 200 lines it changes)
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
		if( t_hat > best) j = j+199; //in the very rare case where we have bigger t than our table
		else{
			for(int k = j ; k < j + 200 ; k++){
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
	double getEwkCorrections(TString url, reco::GenParticleCollection genParticles){

		bool debug = true;
		double kFactor = 1.;

		if(debug) cout<< "Url name : "<< url<<endl;
    if(url.Contains("ZZ") || url.Contains("WZ")){

		//We can only apply corrections to specific diagrams in HZZ, so here we have to ask for some specific criteria:
		// - isMC-->put in the central code
		// - have at least 2 incident quarks and 2 Z
		// - rho < 0.3


		//I need to :
		// - be able to find information about the incident quark (status flag I guess?)
		// - find information on the 4 outgoing leptons
		//
		//



    //create leptons and photons and neutrinos

		reco::GenParticleCollection genIncomingQuarks;
		reco::GenParticleCollection genLeptons;
		//reco::GenParticleCollection genPhotons;
    reco::GenParticleCollection genNeutrinos;
    reco::GenParticleCollection genIncomingProtons;

    for (unsigned int i =0; i < genParticles.size(); i++){
      reco::GenParticle genParticle = genParticles[i];

      if(genParticle.status() == 2) continue; //drop the shower part of the history

			if(genParticle.status() == 3){ //his includes the two incoming colliding particles and partons produced in hard interaction
			 	if(genParticle.pdgId() == 2212) genIncomingProtons.push_back(genParticle);
			 	if(genParticle.pdgId() >= 1 && genParticle.pdgId() <= 5) genIncomingQuarks.push_back(genParticle);
			}

      //if(genParticle.pt() < 5) continue; //ask a minimal pt to remove garbage
      //let's be even more strict. I look at jet of 20GeV, so let's ask for no particles with a pT below 15
      //if(genParticle.pt() <15) continue;

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
    std::sort(genLeptons.begin(), genLeptons.end(), utils::sort_CandidatesByPt);
    //std::sort(genPhotons.begin(), genPhotons.end(), utils::sort_CandidatesByPt);
    std::sort(genNeutrinos.begin(), genNeutrinos.end(), utils::sort_CandidatesByPt);
    std::sort(genIncomingProtons.begin(), genIncomingProtons.end(), utils::sort_CandidatesByPt);

		if(debug){
			cout<<"# of incoming quarks : "<< genIncomingQuarks.size() <<endl;
			cout<<"# of incoming protons : "<< genIncomingProtons.size() <<endl;
			cout<<"# of leptons : "<< genLeptons.size() <<endl;
			cout<<"# of neutrinos : "<< genNeutrinos.size() <<endl;
		}


		//Read the electroweak corrections file and store it in the table
		std::vector<std::vector<float>> Table_EWK;
		Table_EWK = readFile_and_loadEwkTable(url);

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
