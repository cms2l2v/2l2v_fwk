#ifndef PhotonEfficiencySF_h
#define PhotonEfficiencySF_h


#include "UserCode/llvv_fwk/interface/PatUtils.h"
// cf.
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2 
//
class PhotonEfficiencySF
{
 public:
  //
  PhotonEfficiencySF() { }

  //
  ~PhotonEfficiencySF() {}

  //
  std::pair<float,float> getPhotonEfficiency(float pt, float eta, std::string wp, int cutVersion){
      float Abseta=fabs(eta);
      
      std::pair<float,float> eff(1.0,0.04);


        switch(cutVersion){
        case patUtils::CutVersion::Spring15Cut25ns :
          if(wp=="loose"){
              }else if(wp=="medium"){
              }else if(wp=="tight"){
              }
	  break;

        case patUtils::CutVersion::ICHEP16Cut :
          if(wp=="loose"){

        	 if(eta >= -2.5 && eta < -1.566){
        		if( pt < 30){ eff.first=0.910053; eff.second=0.0237691;
        		} else if( pt < 40){ eff.first=0.933894; eff.second=0.0131092;
        		} else if( pt < 50){ eff.first=0.941379; eff.second=0.00703792;
        		} else { eff.first=0.949772; eff.second=0.00992851;
        		}
        	}else if(eta >= -1.566 && eta < -1.444){
        		if( pt < 30){ eff.first=1.00741; eff.second=0.0659558;
        		} else if( pt < 40){ eff.first=0.969697; eff.second=0.0202936;
        		} else if( pt < 50){ eff.first=0.977247; eff.second=0.014855;
        		} else { eff.first=0.978035; eff.second=0.0278135;
        		}
        	}else if(eta >= -1.444 && eta < -1){
        		if( pt < 30){ eff.first=0.95267; eff.second=0.037875;
        		} else if( pt < 40){ eff.first=0.966819; eff.second=0.015697;
        		} else if( pt < 50){ eff.first=0.967991; eff.second=0.00459183;
        		} else { eff.first=0.972315; eff.second=0.00680787;
        		}
        	}else if(eta >= -1 && eta < 0){
        		if( pt < 30){ eff.first=0.936774; eff.second=0.0660392;
        		} else if( pt < 40){ eff.first=0.950762; eff.second=0.0187609;
        		} else if( pt < 50){ eff.first=0.957159; eff.second=0.00415637;
        		} else { eff.first=0.959551; eff.second=0.00810563;
        		}
        	}else if(eta >= 0 && eta < 1){
        		if( pt < 30){ eff.first=0.960154; eff.second=0.0660559;
        		} else if( pt < 40){ eff.first=0.962175; eff.second=0.0187609;
        		} else if( pt < 50){ eff.first=0.968218; eff.second=0.00415637;
        		} else { eff.first=0.968362; eff.second=0.00810563;
        		}
        	}else if(eta >= 1 && eta < 1.444){
        		if( pt < 30){ eff.first=0.969512; eff.second=0.037875;
        		} else if( pt < 40){ eff.first=0.965636; eff.second=0.015697;
        		} else if( pt < 50){ eff.first=0.972315; eff.second=0.00459183;
        		} else { eff.first=0.974501; eff.second=0.00680787;
        		}
        	}else if(eta >= 1.444 && eta < 1.566){
        		if( pt < 30){ eff.first=1.0177; eff.second=0.0659556;
        		} else if( pt < 40){ eff.first=0.979294; eff.second=0.0202936;
        		} else if( pt < 50){ eff.first=0.971591; eff.second=0.0148551;
        		} else { eff.first=0.95397; eff.second=0.0278135;
        		}
        	}else if(eta < 2.5){
        		if( pt < 30){ eff.first=0.915254; eff.second=0.0237685;
        		} else if( pt < 40){ eff.first=0.951249; eff.second=0.0131092;
        		} else if( pt < 50){ eff.first=0.960046; eff.second=0.00703792;
        		} else { eff.first=0.96376; eff.second=0.00992851;
        		}
        	}


          } else if(wp=="medium"){

        	 if(eta >= -2.5 && eta < -1.566){
        		if( pt < 30){ eff.first=0.862992; eff.second=0.0225504;
        		} else if( pt < 40){ eff.first=0.894958; eff.second=0.0129084;
        		} else if( pt < 50){ eff.first=0.902503; eff.second=0.0102019;
        		} else { eff.first=0.914286; eff.second=0.0139271;
        		}
        	}else if(eta >= -1.566 && eta < -1.444){
        		if( pt < 30){ eff.first=0.96846; eff.second=0.0814805;
        		} else if( pt < 40){ eff.first=0.951754; eff.second=0.0185324;
        		} else if( pt < 50){ eff.first=0.947861; eff.second=0.00903424;
        		} else { eff.first=0.937063; eff.second=0.0334198;
        		}
        	}else if(eta >= -1.444 && eta < -1){
        		if( pt < 30){ eff.first=0.926829; eff.second=0.0316868;
        		} else if( pt < 40){ eff.first=0.941567; eff.second=0.0143972;
        		} else if( pt < 50){ eff.first=0.940657; eff.second=0.0078592;
        		} else { eff.first=0.943325; eff.second=0.00791051;
        		}
        	}else if(eta >= -1 && eta < 0){
        		if( pt < 30){ eff.first=0.924499; eff.second=0.0522058;
        		} else if( pt < 40){ eff.first=0.924966; eff.second=0.0167213;
        		} else if( pt < 50){ eff.first=0.931701; eff.second=0.00803215;
        		} else { eff.first=0.933673; eff.second=0.00862579;
        		}
        	}else if(eta >= 0 && eta < 1){
        		if( pt < 30){ eff.first=0.941628; eff.second=0.0522058;
        		} else if( pt < 40){ eff.first=0.936726; eff.second=0.0167213;
        		} else if( pt < 50){ eff.first=0.942783; eff.second=0.00803215;
        		} else { eff.first=0.939744; eff.second=0.00862579;
        		}
        	}else if(eta >= 1 && eta < 1.444){
        		if( pt < 30){ eff.first=0.947522; eff.second=0.0316869;
        		} else if( pt < 40){ eff.first=0.938503; eff.second=0.0143972;
        		} else if( pt < 50){ eff.first=0.946633; eff.second=0.0078592;
        		} else { eff.first=0.942138; eff.second=0.00791051;
        		}
        	}else if(eta >= 1.444 && eta < 1.566){
        		if( pt < 30){ eff.first=0.985047; eff.second=0.0814803;
        		} else if( pt < 40){ eff.first=0.955752; eff.second=0.0185322;
        		} else if( pt < 50){ eff.first=0.952894; eff.second=0.0090344;
        		} else { eff.first=0.942611; eff.second=0.0334198;
        		}
        	}else if(eta < 2.5){
        		if( pt < 30){ eff.first=0.876161; eff.second=0.0225498;
        		} else if( pt < 40){ eff.first=0.916094; eff.second=0.0129084;
        		} else if( pt < 50){ eff.first=0.923277; eff.second=0.0102019;
        		} else { eff.first=0.934615; eff.second=0.0139271;
        		}
        	}


          } else if(wp=="tight"){


         	 if(eta >= -2.5 && eta < -1.566){
         		if( pt < 30){ eff.first=0.851852; eff.second=0.0286009;
         		} else if( pt < 40){ eff.first=0.880719; eff.second=0.0149598;
         		} else if( pt < 50){ eff.first=0.890076; eff.second=0.0123463;
         		} else { eff.first=0.902549; eff.second=0.016774;
         		}
         	}else if(eta >= -1.566 && eta < -1.444){
         		if( pt < 30){ eff.first=0.94577; eff.second=0.0837116;
         		} else if( pt < 40){ eff.first=0.945392; eff.second=0.0246238;
         		} else if( pt < 50){ eff.first=0.932412; eff.second=0.0125848;
         		} else { eff.first=0.917208; eff.second=0.0386943;
         		}
         	}else if(eta >= -1.444 && eta < -1){
         		if( pt < 30){ eff.first=0.897521; eff.second=0.0313796;
         		} else if( pt < 40){ eff.first=0.923781; eff.second=0.0163149;
         		} else if( pt < 50){ eff.first=0.919886; eff.second=0.0110466;
         		} else { eff.first=0.926031; eff.second=0.00913601;
         		}
         	}else if(eta >= -1 && eta < 0){
         		if( pt < 30){ eff.first=0.918033; eff.second=0.0556331;
         		} else if( pt < 40){ eff.first=0.905363; eff.second=0.0154303;
         		} else if( pt < 50){ eff.first=0.911635; eff.second=0.0103968;
         		} else { eff.first=0.914493; eff.second=0.0101519;
         		}
         	}else if(eta >= 0 && eta < 1){
         		if( pt < 30){ eff.first=0.92922; eff.second=0.0556331;
         		} else if( pt < 40){ eff.first=0.920128; eff.second=0.0154303;
         		} else if( pt < 50){ eff.first=0.928465; eff.second=0.0103968;
         		} else { eff.first=0.926686; eff.second=0.0101519;
         		}
         	}else if(eta >= 1 && eta < 1.444){
         		if( pt < 30){ eff.first=0.931973; eff.second=0.0313797;
         		} else if( pt < 40){ eff.first=0.921538; eff.second=0.0163149;
         		} else if( pt < 50){ eff.first=0.929191; eff.second=0.0110467;
         		} else { eff.first=0.926031; eff.second=0.00913601;
         		}
         	}else if(eta >= 1.444 && eta < 1.566){
         		if( pt < 30){ eff.first=0.975556; eff.second=0.0837115;
         		} else if( pt < 40){ eff.first=0.946918; eff.second=0.0246234;
         		} else if( pt < 50){ eff.first=0.944186; eff.second=0.0125839;
         		} else { eff.first=0.941374; eff.second=0.0386943;
         		}
         	}else if(eta < 2.5){
         		if( pt < 30){ eff.first=0.867754; eff.second=0.0286008;
         		} else if( pt < 40){ eff.first=0.900958; eff.second=0.0149598;
         		} else if( pt < 50){ eff.first=0.90991; eff.second=0.0123463;
         		} else { eff.first=0.920354; eff.second=0.016774;
         		}
         	}
         

	  }
        break;

	case patUtils::CutVersion::Moriond17Cut :
          if(wp=="loose"){

        	 if(eta >= -2.5 && eta < -1.566){
        		if( pt < 30){ eff.first=0.910053; eff.second=0.0237691;
        		} else if( pt < 40){ eff.first=0.933894; eff.second=0.0131092;
        		} else if( pt < 50){ eff.first=0.941379; eff.second=0.00703792;
        		} else { eff.first=0.949772; eff.second=0.00992851;
        		}
        	}else if(eta >= -1.566 && eta < -1.444){
        		if( pt < 30){ eff.first=1.00741; eff.second=0.0659558;
        		} else if( pt < 40){ eff.first=0.969697; eff.second=0.0202936;
        		} else if( pt < 50){ eff.first=0.977247; eff.second=0.014855;
        		} else { eff.first=0.978035; eff.second=0.0278135;
        		}
        	}else if(eta >= -1.444 && eta < -1){
        		if( pt < 30){ eff.first=0.95267; eff.second=0.037875;
        		} else if( pt < 40){ eff.first=0.966819; eff.second=0.015697;
        		} else if( pt < 50){ eff.first=0.967991; eff.second=0.00459183;
        		} else { eff.first=0.972315; eff.second=0.00680787;
        		}
        	}else if(eta >= -1 && eta < 0){
        		if( pt < 30){ eff.first=0.936774; eff.second=0.0660392;
        		} else if( pt < 40){ eff.first=0.950762; eff.second=0.0187609;
        		} else if( pt < 50){ eff.first=0.957159; eff.second=0.00415637;
        		} else { eff.first=0.959551; eff.second=0.00810563;
        		}
        	}else if(eta >= 0 && eta < 1){
        		if( pt < 30){ eff.first=0.960154; eff.second=0.0660559;
        		} else if( pt < 40){ eff.first=0.962175; eff.second=0.0187609;
        		} else if( pt < 50){ eff.first=0.968218; eff.second=0.00415637;
        		} else { eff.first=0.968362; eff.second=0.00810563;
        		}
        	}else if(eta >= 1 && eta < 1.444){
        		if( pt < 30){ eff.first=0.969512; eff.second=0.037875;
        		} else if( pt < 40){ eff.first=0.965636; eff.second=0.015697;
        		} else if( pt < 50){ eff.first=0.972315; eff.second=0.00459183;
        		} else { eff.first=0.974501; eff.second=0.00680787;
        		}
        	}else if(eta >= 1.444 && eta < 1.566){
        		if( pt < 30){ eff.first=1.0177; eff.second=0.0659556;
        		} else if( pt < 40){ eff.first=0.979294; eff.second=0.0202936;
        		} else if( pt < 50){ eff.first=0.971591; eff.second=0.0148551;
        		} else { eff.first=0.95397; eff.second=0.0278135;
        		}
        	}else if(eta < 2.5){
        		if( pt < 30){ eff.first=0.915254; eff.second=0.0237685;
        		} else if( pt < 40){ eff.first=0.951249; eff.second=0.0131092;
        		} else if( pt < 50){ eff.first=0.960046; eff.second=0.00703792;
        		} else { eff.first=0.96376; eff.second=0.00992851;
        		}
        	}


          } else if(wp=="medium"){

        	 if(eta >= -2.5 && eta < -1.566){
        		if( pt < 30){ eff.first=0.862992; eff.second=0.0225504;
        		} else if( pt < 40){ eff.first=0.894958; eff.second=0.0129084;
        		} else if( pt < 50){ eff.first=0.902503; eff.second=0.0102019;
        		} else { eff.first=0.914286; eff.second=0.0139271;
        		}
        	}else if(eta >= -1.566 && eta < -1.444){
        		if( pt < 30){ eff.first=0.96846; eff.second=0.0814805;
        		} else if( pt < 40){ eff.first=0.951754; eff.second=0.0185324;
        		} else if( pt < 50){ eff.first=0.947861; eff.second=0.00903424;
        		} else { eff.first=0.937063; eff.second=0.0334198;
        		}
        	}else if(eta >= -1.444 && eta < -1){
        		if( pt < 30){ eff.first=0.926829; eff.second=0.0316868;
        		} else if( pt < 40){ eff.first=0.941567; eff.second=0.0143972;
        		} else if( pt < 50){ eff.first=0.940657; eff.second=0.0078592;
        		} else { eff.first=0.943325; eff.second=0.00791051;
        		}
        	}else if(eta >= -1 && eta < 0){
        		if( pt < 30){ eff.first=0.924499; eff.second=0.0522058;
        		} else if( pt < 40){ eff.first=0.924966; eff.second=0.0167213;
        		} else if( pt < 50){ eff.first=0.931701; eff.second=0.00803215;
        		} else { eff.first=0.933673; eff.second=0.00862579;
        		}
        	}else if(eta >= 0 && eta < 1){
        		if( pt < 30){ eff.first=0.941628; eff.second=0.0522058;
        		} else if( pt < 40){ eff.first=0.936726; eff.second=0.0167213;
        		} else if( pt < 50){ eff.first=0.942783; eff.second=0.00803215;
        		} else { eff.first=0.939744; eff.second=0.00862579;
        		}
        	}else if(eta >= 1 && eta < 1.444){
        		if( pt < 30){ eff.first=0.947522; eff.second=0.0316869;
        		} else if( pt < 40){ eff.first=0.938503; eff.second=0.0143972;
        		} else if( pt < 50){ eff.first=0.946633; eff.second=0.0078592;
        		} else { eff.first=0.942138; eff.second=0.00791051;
        		}
        	}else if(eta >= 1.444 && eta < 1.566){
        		if( pt < 30){ eff.first=0.985047; eff.second=0.0814803;
        		} else if( pt < 40){ eff.first=0.955752; eff.second=0.0185322;
        		} else if( pt < 50){ eff.first=0.952894; eff.second=0.0090344;
        		} else { eff.first=0.942611; eff.second=0.0334198;
        		}
        	}else if(eta < 2.5){
        		if( pt < 30){ eff.first=0.876161; eff.second=0.0225498;
        		} else if( pt < 40){ eff.first=0.916094; eff.second=0.0129084;
        		} else if( pt < 50){ eff.first=0.923277; eff.second=0.0102019;
        		} else { eff.first=0.934615; eff.second=0.0139271;
        		}
        	}


          } else if(wp=="tight"){


         	 if(eta >= -2.5 && eta < -2.0){
         		if( pt < 35){ eff.first=0.903; eff.second=0.010;
         		} else if( pt < 50){ eff.first=0.916; eff.second=0.019;
         		} else if( pt < 90){ eff.first=0.921; eff.second=0.029;
         		} else if( pt < 150){ eff.first=0.940; eff.second=0.023;
         		} else { eff.first=0.915; eff.second=0.136;
         		}
         	}else if(eta >= -2.0 && eta < -1.566){
         		if( pt < 35){ eff.first=0.903; eff.second=0.010;
         		} else if( pt < 50){ eff.first=0.916; eff.second=0.019;
         		} else if( pt < 90){ eff.first=0.921; eff.second=0.029;
         		} else if( pt < 150){ eff.first=0.940; eff.second=0.023;
         		} else { eff.first=0.915; eff.second=0.136;
			}	
         	}else if(eta >= -1.566 && eta < -1.444){
         		if( pt < 35){ eff.first=1.0; eff.second=1.414;
         		} else if( pt < 50){ eff.first=1.0; eff.second=1.414;
         		} else if( pt < 90){ eff.first=1.0; eff.second=1.414;
         		} else if( pt < 150){ eff.first=1.0; eff.second=1.414;
         		} else { eff.first=1.0; eff.second=1.414;
         		}
         	}else if(eta >= -1.444 && eta < -0.8){
         		if( pt < 35){ eff.first=0.980; eff.second=0.013;
         		} else if( pt < 50){ eff.first=0.977; eff.second=0.011;
         		} else if( pt < 90){ eff.first=0.968; eff.second=0.026;
         		} else if( pt < 150){ eff.first=1.011; eff.second=0.016;
         		} else { eff.first=1.012; eff.second=0.025;
         		}
         	}else if(eta >= -0.8 && eta < 0.0){
         		if( pt < 35){ eff.first=0.969; eff.second=0.013;
         		} else if( pt < 50){ eff.first=0.963; eff.second=0.011;
         		} else if( pt < 90){ eff.first=0.953; eff.second=0.026;
         		} else if( pt < 150){ eff.first=0.999; eff.second=0.011;
         		} else { eff.first=0.989; eff.second=0.021;
         		}
         	}else if(eta >= 0.0 && eta < 0.8){
         		if( pt < 35){ eff.first=0.986; eff.second=0.028;
         		} else if( pt < 50){ eff.first=0.981; eff.second=0.011;
         		} else if( pt < 90){ eff.first=0.971; eff.second=0.026;
         		} else if( pt < 150){ eff.first=1.007; eff.second=0.011;
         		} else { eff.first=1.014; eff.second=0.021;
         		}
         	}else if(eta >= 0.8 && eta < 1.444){
         		if( pt < 35){ eff.first=0.990; eff.second=0.013;
         		} else if( pt < 50){ eff.first=0.980; eff.second=0.011;
         		} else if( pt < 90){ eff.first=0.971; eff.second=0.026;
         		} else if( pt < 150){ eff.first=1.026; eff.second=0.016;
         		} else { eff.first=0.977; eff.second=0.025;
         		}
         	}else if(eta >= 1.444 && eta < 1.566){
         		if( pt < 35){ eff.first=1.0; eff.second=1.414;
         		} else if( pt < 50){ eff.first=1.0; eff.second=1.414;
         		} else if( pt < 90){ eff.first=1.0; eff.second=1.414;
         		} else if( pt < 150){ eff.first=1.0; eff.second=1.414;
         		} else { eff.first=1.0; eff.second=1.414;
         		}
         	}else if(eta >= 1.566 && eta < 2.0){
         		if( pt < 35){ eff.first=0.937; eff.second=0.005;
         		} else if( pt < 50){ eff.first=0.957; eff.second=0.009;
         		} else if( pt < 90){ eff.first=0.948; eff.second=0.033;
         		} else if( pt < 150){ eff.first=0.983; eff.second=0.018;
         		} else { eff.first=1.006; eff.second=0.054;
         		}
         	}else if(eta < 2.5){
         		if( pt < 35){ eff.first=0.904; eff.second=0.010;
         		} else if( pt < 50){ eff.first=0.927; eff.second=0.019;
         		} else if( pt < 90){ eff.first=0.932; eff.second=0.029;
         		} else if( pt < 150){ eff.first=0.983; eff.second=0.023;
         		} else { eff.first=0.976; eff.second=0.136;
         		}
         	}
         

	  }
          break;
        }
      return eff;
    }

    
private:
    
};


#endif
