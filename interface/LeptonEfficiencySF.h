#ifndef LeptonEfficiencySF_h
#define LeptonEfficiencySF_h


#include "UserCode/llvv_fwk/interface/PatUtils.h"
// cf.
// https://twiki.cern.ch/twiki/bin/view/Main/EGammaScaleFactors2012#2012_8_TeV_data_53X
//
class LeptonEfficiencySF
{
 public:
  //
  LeptonEfficiencySF() { }

  //
  ~LeptonEfficiencySF() {}

  //
  std::pair<float,float> getLeptonEfficiency(float pt, float eta, int id, std::string wp, int cutVersion){
      float Abseta=fabs(eta);
      id=abs(id);
      
      std::pair<float,float> eff(1.0,0.04);


      switch(id){
      case 11:

        switch(cutVersion){
        case patUtils::CutVersion::Spring15Cut25ns :
          if(wp=="loose"){
                if( eta>= -2.5 && eta<-2.0){
                    if(pt<30)      { eff.first=0.8775;  eff.second=0.0067; 
                    }else if(pt<40) { eff.first=1.0271;  eff.second=0.0045; 
                    }else if(pt<50) { eff.first=1.0186;  eff.second=0.0033; 
                    }else           { eff.first=1.0114;  eff.second=0.0073;
                    }                 
                }else if(eta>=-2.0 && eta<-1.566){
                    if(pt<30)      { eff.first=0.9905;  eff.second=0.0067; 
                    }else if(pt<40) { eff.first=1.0024;  eff.second=0.0043; 
                    }else if(pt<50) { eff.first=1.0034;  eff.second=0.0032; 
                    }else           { eff.first=1.0044;  eff.second=0.0074;
                    }                  
                }else if(eta>=-1.566 && eta<-1.444){
                    if(pt<30)      { eff.first=1.0083;  eff.second=0.0185; 
                    }else if(pt<40) { eff.first=1.0133;  eff.second=0.0133; 
                    }else if(pt<50) { eff.first=0.9938;  eff.second=0.0079; 
                    }else           { eff.first=0.9805;  eff.second=0.0129;
                    }                  
                }else if(eta>=-1.444 && eta<-0.8){
                    if(pt<30)      { eff.first=1.0199;  eff.second=0.0072; 
                    }else if(pt<40) { eff.first=0.9919;  eff.second=0.0026; 
                    }else if(pt<50) { eff.first=0.9881;  eff.second=0.0015; 
                    }else           { eff.first=0.9798;  eff.second=0.0030;
                    }                  
                }else if(eta>-0.8 && eta<0){
                    if(pt<30)      { eff.first=0.9871;  eff.second=0.0069; 
                    }else if(pt<40) { eff.first=0.9752;  eff.second=0.0016; 
                    }else if(pt<50) { eff.first=0.9775;  eff.second=0.0015; 
                    }else           { eff.first=0.9778;  eff.second=0.0030;
                    }                  
                }else if(eta>0 && eta<0.8){
                    if(pt<30)      { eff.first=0.7954;  eff.second=0.0043; 
                    }else if(pt<40) { eff.first=0.9774;  eff.second=0.0016; 
                    }else if(pt<50) { eff.first=0.9786;  eff.second=0.0015; 
                    }else           { eff.first=0.9789;  eff.second=0.0030;
                    }                  
                }else if(eta>0.8 && eta<1.444){
                    if(pt<30)      { eff.first=1.0146;  eff.second=0.0072; 
                    }else if(pt<40) { eff.first=0.9988;  eff.second=0.0016; 
                    }else if(pt<50) { eff.first=0.9881;  eff.second=0.0015; 
                    }else           { eff.first=0.9861;  eff.second=0.0030;
                   }                  
                }else if(eta>1.444 && eta<1.556){
                    if(pt<30)      { eff.first=1.0393;  eff.second=0.0544; 
                    }else if(pt<40) { eff.first=0.9985;  eff.second=0.0133; 
                    }else if(pt<50) { eff.first=0.9739;  eff.second=0.0079; 
                    }else           { eff.first=0.9854;  eff.second=0.0169;
                    }                  
                  }else if(eta>1.556 && eta<2.0){
                    if(pt<30)      { eff.first=1.0233;  eff.second=0.0130; 
                    }else if(pt<40) { eff.first=0.9917;  eff.second=0.0043; 
                    }else if(pt<50) { eff.first=1.0034;  eff.second=0.0032; 
                    }else           { eff.first=1.0077;  eff.second=0.0063;
                    } 
                }else{
                    if(pt<30)      { eff.first=1.0306;  eff.second=0.0082; 
                    }else if(pt<40) { eff.first=1.0246;  eff.second=0.0045; 
                    }else if(pt<50) { eff.first=1.0140;  eff.second=0.0042; 
                    }else           { eff.first=1.0149;  eff.second=0.0083;
                    }
                }
              }else if(wp=="medium"){
                if(eta>=-2.5 && eta<-2.0){
                    if(pt<30)      { eff.first=1.0536;  eff.second=0.0099; 
                    }else if(pt<40) { eff.first=1.0123;  eff.second=0.0050; 
                    }else if(pt<50) { eff.first=1.0087;  eff.second=0.0045; 
                    }else           { eff.first=1.0060;  eff.second=0.0077;
                    }                  
                  }else if(eta>=-2.0 && eta<-1.566){
                    if(pt<30)      { eff.first=0.9695;  eff.second=0.0117; 
                    }else if(pt<40) { eff.first=0.9713;  eff.second=0.0047; 
                    }else if(pt<50) { eff.first=0.9833;  eff.second=0.0033; 
                    }else           { eff.first=0.9851;  eff.second=0.0057;
                    }                  
                }else if(eta>=-1.566 && eta<-1.444){
                    if(pt<30)      { eff.first=1.0404;  eff.second=0.0231; 
                    }else if(pt<40) { eff.first=1.0248;  eff.second=0.0120; 
                    }else if(pt<50) { eff.first=0.9902;  eff.second=0.0090; 
                    }else           { eff.first=0.9489;  eff.second=0.0143;
                    }                  
                }else if(eta>=-1.444 && eta<-0.8){
                    if(pt<30)      { eff.first=1.0594;  eff.second=0.0087; 
                    }else if(pt<40) { eff.first=0.9909;  eff.second=0.0029; 
                    }else if(pt<50) { eff.first=0.9741;  eff.second=0.0016; 
                    }else           { eff.first=0.9581;  eff.second=0.0040;
                    }                  
                }else if(eta>=-0.8 && eta<0){
                    if(pt<30)      { eff.first=1.0226;  eff.second=0.0081; 
                    }else if(pt<40) { eff.first=0.9786;  eff.second=0.0028; 
                    }else if(pt<50) { eff.first=0.9698;  eff.second=0.0016; 
                    }else           { eff.first=0.9675;  eff.second=0.0031;
                    }                  
                }else if(eta>=0 && eta<0.8){
                    if(pt<30)      { eff.first=1.0180;  eff.second=0.0067; 
                    }else if(pt<40) { eff.first=0.9862;  eff.second=0.0028; 
                    }else if(pt<50) { eff.first=0.9744;  eff.second=0.0016; 
                    }else           { eff.first=0.9708;  eff.second=0.0031;
                    }                  
                }else if(eta>=0.8 && eta<1.444){
                    if(pt<30)      { eff.first=1.0774;  eff.second=0.0276; 
                    }else if(pt<40) { eff.first=1.0065;  eff.second=0.0029; 
                    }else if(pt<50) { eff.first=0.9764;  eff.second=0.0017; 
                    }else           { eff.first=0.9648;  eff.second=0.0041;
                    }                  
                }else if(eta>=1.444 && eta<1.566){
                    if(pt<30)      { eff.first=1.1880;  eff.second=0.0264; 
                    }else if(pt<40) { eff.first=0.9910;  eff.second=0.0154; 
                    }else if(pt<50) { eff.first=0.9562;  eff.second=0.0400; 
                    }else           { eff.first=0.9744;  eff.second=0.0180;
                    }                  
                  }else if(eta>=1.566 && eta<2.0){
                    if(pt<30)      { eff.first=0.9804;  eff.second=0.0073; 
                    }else if(pt<40) { eff.first=0.9619;  eff.second=0.0047; 
                    }else if(pt<50) { eff.first=0.9785;  eff.second=0.0033; 
                    }else           { eff.first=0.9943;  eff.second=0.0057;
                    }                  
                }else {
                    if(pt<30)      { eff.first=1.0265;  eff.second=0.0097; 
                    }else if(pt<40) { eff.first=1.0096;  eff.second=0.0050; 
                    }else if(pt<50) { eff.first=1.0013;  eff.second=0.0045; 
                    }else           { eff.first=1.0048;  eff.second=0.0087;
                    }
                }
              }else if(wp=="tight"){
                if(eta>=-2.5 && eta<-2.0){
                    if(pt<30)      { eff.first=1.0586;  eff.second=0.0114; 
                    }else if(pt<40) { eff.first=1.0246;  eff.second=0.0060; 
                    }else if(pt<50) { eff.first=1.0101;  eff.second=0.0052; 
                    }else           { eff.first=0.9987;  eff.second=0.0097;
                    }                  
                }else if(eta>=-2.0 && eta<-1.566){
                    if(pt<30)      { eff.first=0.9838;  eff.second=0.0128; 
                    }else if(pt<40) { eff.first=0.9686;  eff.second=0.0056; 
                    }else if(pt<50) { eff.first=0.9819;  eff.second=0.0050; 
                    }else           { eff.first=0.9806;  eff.second=0.0075;
                    }                  
                  }else if(eta>=-1.566 && eta<-1.444){
                    if(pt<30)      { eff.first=1.1862;  eff.second=0.0498; 
                    }else if(pt<40) { eff.first=1.0152;  eff.second=0.0146; 
                    }else if(pt<50) { eff.first=0.9934;  eff.second=0.0156; 
                    }else           { eff.first=0.9449;  eff.second=0.0176;
                    }                  
                }else if(eta>=-1.444 && eta<-0.8){
                    if(pt<30)      { eff.first=1.0502;  eff.second=0.0109; 
                    }else if(pt<40) { eff.first=0.9892;  eff.second=0.0034; 
                    }else if(pt<50) { eff.first=0.9704;  eff.second=0.0030; 
                    }else           { eff.first=0.9562;  eff.second=0.0044;
                    }                  
                }else if(eta>=-0.8 && eta<0){
                    if(pt<30)       { eff.first=1.0188;  eff.second=0.0078; 
                    }else if(pt<40) { eff.first=0.9644;  eff.second=0.0033; 
                    }else if(pt<50) { eff.first=0.9565;  eff.second=0.0018; 
                    }else           { eff.first=0.9544;  eff.second=0.0034;
                    }                  
                }else if(eta>=0 && eta<0.8){
                    if(pt<30)       { eff.first=0.9944;  eff.second=0.0059; 
                    }else if(pt<40) { eff.first=0.9674;  eff.second=0.0021; 
                    }else if(pt<50) { eff.first=0.9617;  eff.second=0.0018; 
                    }else           { eff.first=0.9617;  eff.second=0.0034;
                    }                  
                  }else if(eta>=0.8 && eta<1.444){
                    if(pt<30)       { eff.first=1.0567;  eff.second=0.0148; 
                    }else if(pt<40) { eff.first=0.9969;  eff.second=0.0035; 
                    }else if(pt<50) { eff.first=0.9730;  eff.second=0.0030; 
                    }else           { eff.first=0.9561;  eff.second=0.0045;
                    }                  
                }else if(eta>=1.444 && eta<1.556){
                    if(pt<30)       { eff.first=1.0601;  eff.second=0.0439; 
                    }else if(pt<40) { eff.first=0.9868;  eff.second=0.1167; 
                    }else if(pt<50) { eff.first=0.9568;  eff.second=0.0110; 
                    }else           { eff.first=0.9493;  eff.second=0.0164;
                    }                  
                }else if(eta>=1.556 && eta<2.0){
                    if(pt<30)       { eff.first=1.0225;  eff.second=0.0130; 
                    }else if(pt<40) { eff.first=0.9670;  eff.second=0.0056; 
                    }else if(pt<50) { eff.first=0.9847;  eff.second=0.0050; 
                    }else           { eff.first=0.9922;  eff.second=0.0075; 
                    }                  
                }else{
                    if(pt<30)       { eff.first=1.0560;  eff.second=0.0113; 
                    }else if(pt<40) { eff.first=1.0229;  eff.second=0.0459; 
                    }else if(pt<50) { eff.first=1.0043;  eff.second=0.0052; 
                    }else           { eff.first=1.0013;  eff.second=0.0097; 
                    }
                  }
              }
	  break;

        case patUtils::CutVersion::ICHEP16Cut :
          if(wp=="loose"){
              if( eta >= -2.5 && eta < -2.0){
              		if( pt < 20.0){ eff.first=0.917582; eff.second=0.019551;
              	} else if( pt < 30.0){ eff.first=0.950617; eff.second=0.014836;
              	} else if( pt < 40.0){ eff.first=0.991474; eff.second=0.007577;
              	} else if( pt < 50.0){ eff.first=1.001157; eff.second=0.005660;
              	} else { eff.first=1.004561; eff.second=0.007909;
              	}
              }else if( eta >= -2.0 && eta < -1.6){
              		if( pt < 20.0){ eff.first=0.817717; eff.second=0.040757;
              	} else if( pt < 30.0){ eff.first=0.905138; eff.second=0.020285;
              	} else if( pt < 40.0){ eff.first=0.962485; eff.second=0.009990;
              	} else if( pt < 50.0){ eff.first=0.983259; eff.second=0.003530;
              	} else { eff.first=0.990099; eff.second=0.004044;
              	}
              }else if( eta >= -1.6 && eta < -1.4){
              		if( pt < 20.0){ eff.first=0.898678; eff.second=0.096118;
              	} else if( pt < 30.0){ eff.first=0.943609; eff.second=0.076008;
              	} else if( pt < 40.0){ eff.first=0.972740; eff.second=0.021814;
              	} else if( pt < 50.0){ eff.first=0.973171; eff.second=0.006222;
              	} else { eff.first=0.980746; eff.second=0.008877;
              	}
              }else if( eta >= -1.4 && eta < -0.8){
              		if( pt < 20.0){ eff.first=0.965135; eff.second=0.019147;
              	} else if( pt < 30.0){ eff.first=0.957031; eff.second=0.022926;
              	} else if( pt < 40.0){ eff.first=0.976932; eff.second=0.008641;
              	} else if( pt < 50.0){ eff.first=0.979303; eff.second=0.002928;
              	} else { eff.first=0.981740; eff.second=0.002213;
              	}
              }else if( eta >= -0.8 && eta < 0.0){
              		if( pt < 20.0){ eff.first=0.961353; eff.second=0.028838;
              	} else if( pt < 30.0){ eff.first=0.956975; eff.second=0.023614;
              	} else if( pt < 40.0){ eff.first=0.966705; eff.second=0.005886;
              	} else if( pt < 50.0){ eff.first=0.971678; eff.second=0.002726;
              	} else { eff.first=0.973205; eff.second=0.002008;
              	}
              }else if( eta >= 0.0 && eta < 0.8){
              		if( pt < 20.0){ eff.first=0.919431; eff.second=0.028838;
              	} else if( pt < 30.0){ eff.first=0.971466; eff.second=0.023614;
              	} else if( pt < 40.0){ eff.first=0.977012; eff.second=0.005886;
              	} else if( pt < 50.0){ eff.first=0.981441; eff.second=0.002726;
              	} else { eff.first=0.984946; eff.second=0.002008;
              	}
              }else if( eta >= 0.8 && eta < 1.4){
              		if( pt < 20.0){ eff.first=0.957075; eff.second=0.019147;
              	} else if( pt < 30.0){ eff.first=0.960938; eff.second=0.022926;
              	} else if( pt < 40.0){ eff.first=0.971231; eff.second=0.008641;
              	} else if( pt < 50.0){ eff.first=0.977199; eff.second=0.002928;
              	} else { eff.first=0.982833; eff.second=0.002213;
              	}
              }else if( eta >= 1.4 && eta < 1.6){
              		if( pt < 20.0){ eff.first=0.978520; eff.second=0.096117;
              	} else if( pt < 30.0){ eff.first=0.954802; eff.second=0.076008;
              	} else if( pt < 40.0){ eff.first=0.964029; eff.second=0.021814;
              	} else if( pt < 50.0){ eff.first=0.964591; eff.second=0.006222;
              	} else { eff.first=0.976077; eff.second=0.008877;
              	}
              }else if( eta >= 1.6 && eta < 2.0){
              		if( pt < 20.0){ eff.first=0.864111; eff.second=0.040757;
              	} else if( pt < 30.0){ eff.first=0.915436; eff.second=0.020285;
              	} else if( pt < 40.0){ eff.first=0.968310; eff.second=0.009990;
              	} else if( pt < 50.0){ eff.first=0.985523; eff.second=0.003530;
              	} else { eff.first=1.002203; eff.second=0.004044;
              	}
              }else {
              		if( pt < 20.0){ eff.first=0.896959; eff.second=0.019551;
              	} else if( pt < 30.0){ eff.first=0.945551; eff.second=0.014836;
              	} else if( pt < 40.0){ eff.first=0.976134; eff.second=0.007577;
              	} else if( pt < 50.0){ eff.first=0.990868; eff.second=0.005660;
              	} else { eff.first=0.996606; eff.second=0.007909;
              	}
              }

          } else if(wp=="medium"){
              if( eta >= -2.5 && eta < -2.0){
              		if( pt < 20.0){ eff.first=0.858209; eff.second=0.025860;
              	} else if( pt < 30.0){ eff.first=0.898527; eff.second=0.016977;
              	} else if( pt < 40.0){ eff.first=0.951219; eff.second=0.008878;
              	} else if( pt < 50.0){ eff.first=0.971429; eff.second=0.008115;
              	} else { eff.first=0.977273; eff.second=0.011493;
              	}
              }else if( eta >= -2.0 && eta < -1.6){
              		if( pt < 20.0){ eff.first=0.748848; eff.second=0.052968;
              	} else if( pt < 30.0){ eff.first=0.854003; eff.second=0.020792;
              	} else if( pt < 40.0){ eff.first=0.923077; eff.second=0.010527;
              	} else if( pt < 50.0){ eff.first=0.953405; eff.second=0.005876;
              	} else { eff.first=0.966705; eff.second=0.005664;
              	}
              }else if( eta >= -1.6 && eta < -1.4){
              		if( pt < 20.0){ eff.first=0.879630; eff.second=0.093792;
              	} else if( pt < 30.0){ eff.first=0.952267; eff.second=0.061282;
              	} else if( pt < 40.0){ eff.first=0.957045; eff.second=0.019637;
              	} else if( pt < 50.0){ eff.first=0.951253; eff.second=0.007580;
              	} else { eff.first=0.967347; eff.second=0.011131;
              	}
              }else if( eta >= -1.4 && eta < -0.8){
              		if( pt < 20.0){ eff.first=0.947581; eff.second=0.023546;
              	} else if( pt < 30.0){ eff.first=0.944870; eff.second=0.020984;
              	} else if( pt < 40.0){ eff.first=0.963636; eff.second=0.007988;
              	} else if( pt < 50.0){ eff.first=0.965144; eff.second=0.005550;
              	} else { eff.first=0.963869; eff.second=0.005485;
              	}
              }else if( eta >= -0.8 && eta < 0.0){
              		if( pt < 20.0){ eff.first=0.919386; eff.second=0.030269;
              	} else if( pt < 30.0){ eff.first=0.934718; eff.second=0.020054;
              	} else if( pt < 40.0){ eff.first=0.946701; eff.second=0.006414;
              	} else if( pt < 50.0){ eff.first=0.951537; eff.second=0.004149;
              	} else { eff.first=0.952928; eff.second=0.005486;
              	}
              }else if( eta >= 0.0 && eta < 0.8){
              		if( pt < 20.0){ eff.first=0.902622; eff.second=0.030269;
              	} else if( pt < 30.0){ eff.first=0.952733; eff.second=0.020054;
              	} else if( pt < 40.0){ eff.first=0.966879; eff.second=0.006413;
              	} else if( pt < 50.0){ eff.first=0.970273; eff.second=0.004149;
              	} else { eff.first=0.972350; eff.second=0.005486;
              	}
              }else if( eta >= 0.8 && eta < 1.4){
              		if( pt < 20.0){ eff.first=0.932000; eff.second=0.023546;
              	} else if( pt < 30.0){ eff.first=0.933131; eff.second=0.020984;
              	} else if( pt < 40.0){ eff.first=0.952196; eff.second=0.007988;
              	} else if( pt < 50.0){ eff.first=0.957092; eff.second=0.005550;
              	} else { eff.first=0.964037; eff.second=0.005485;
              	}
              }else if( eta >= 1.4 && eta < 1.6){
              		if( pt < 20.0){ eff.first=1.010601; eff.second=0.093791;
              	} else if( pt < 30.0){ eff.first=0.920673; eff.second=0.061282;
              	} else if( pt < 40.0){ eff.first=0.946274; eff.second=0.019637;
              	} else if( pt < 50.0){ eff.first=0.934266; eff.second=0.007580;
              	} else { eff.first=0.939678; eff.second=0.011132;
              	}
              }else if( eta >= 1.6 && eta < 2.0){
              		if( pt < 20.0){ eff.first=0.814554; eff.second=0.052968;
              	} else if( pt < 30.0){ eff.first=0.878205; eff.second=0.020792;
              	} else if( pt < 40.0){ eff.first=0.937255; eff.second=0.010527;
              	} else if( pt < 50.0){ eff.first=0.961905; eff.second=0.005876;
              	} else { eff.first=0.986175; eff.second=0.005664;
              	}
              }else {
              		if( pt < 20.0){ eff.first=0.858427; eff.second=0.025861;
              	} else if( pt < 30.0){ eff.first=0.911076; eff.second=0.016977;
              	} else if( pt < 40.0){ eff.first=0.946124; eff.second=0.008878;
              	} else if( pt < 50.0){ eff.first=0.968331; eff.second=0.008115;
              	} else { eff.first=0.978698; eff.second=0.011493;
              	}
              }

          } else if(wp=="tight"){
              if( eta >= -2.5 && eta < -2.0){
              		if( pt < 20.0){ eff.first=0.813559; eff.second=0.045692;
              	} else if( pt < 30.0){ eff.first=0.872690; eff.second=0.020570;
              	} else if( pt < 40.0){ eff.first=0.930195; eff.second=0.014515;
              	} else if( pt < 50.0){ eff.first=0.942611; eff.second=0.009042;
              	} else { eff.first=0.954116; eff.second=0.013293;
              	}
              }else if( eta >= -2.0 && eta < -1.6){
              		if( pt < 20.0){ eff.first=0.773700; eff.second=0.049056;
              	} else if( pt < 30.0){ eff.first=0.837255; eff.second=0.023665;
              	} else if( pt < 40.0){ eff.first=0.901563; eff.second=0.014323;
              	} else if( pt < 50.0){ eff.first=0.933426; eff.second=0.005379;
              	} else { eff.first=0.946753; eff.second=0.010596;
              	}
              }else if( eta >= -1.6 && eta < -1.4){
              		if( pt < 20.0){ eff.first=0.892430; eff.second=0.090915;
              	} else if( pt < 30.0){ eff.first=0.937685; eff.second=0.062067;
              	} else if( pt < 40.0){ eff.first=0.946392; eff.second=0.021237;
              	} else if( pt < 50.0){ eff.first=0.945073; eff.second=0.009198;
              	} else { eff.first=0.956790; eff.second=0.013983;
              	}
              }else if( eta >= -1.4 && eta < -0.8){
              		if( pt < 20.0){ eff.first=0.940860; eff.second=0.027953;
              	} else if( pt < 30.0){ eff.first=0.934363; eff.second=0.024353;
              	} else if( pt < 40.0){ eff.first=0.954829; eff.second=0.010925;
              	} else if( pt < 50.0){ eff.first=0.953103; eff.second=0.009542;
              	} else { eff.first=0.953428; eff.second=0.007131;
              	}
              }else if( eta >= -0.8 && eta < 0.0){
              		if( pt < 20.0){ eff.first=0.901515; eff.second=0.026808;
              	} else if( pt < 30.0){ eff.first=0.925373; eff.second=0.027809;
              	} else if( pt < 40.0){ eff.first=0.932127; eff.second=0.006671;
              	} else if( pt < 50.0){ eff.first=0.932615; eff.second=0.007123;
              	} else { eff.first=0.939009; eff.second=0.011901;
              	}
              }else if( eta >= 0.0 && eta < 0.8){
              		if( pt < 20.0){ eff.first=0.916256; eff.second=0.027175;
              	} else if( pt < 30.0){ eff.first=0.951941; eff.second=0.027809;
              	} else if( pt < 40.0){ eff.first=0.948563; eff.second=0.006671;
              	} else if( pt < 50.0){ eff.first=0.953867; eff.second=0.007123;
              	} else { eff.first=0.954140; eff.second=0.011901;
              	}
              }else if( eta >= 0.8 && eta < 1.4){
              		if( pt < 20.0){ eff.first=0.933333; eff.second=0.027953;
              	} else if( pt < 30.0){ eff.first=0.934741; eff.second=0.024253;
              	} else if( pt < 40.0){ eff.first=0.939628; eff.second=0.010925;
              	} else if( pt < 50.0){ eff.first=0.945280; eff.second=0.009542;
              	} else { eff.first=0.953608; eff.second=0.007131;
              	}
              }else if( eta >= 1.4 && eta < 1.6){
              		if( pt < 20.0){ eff.first=1.014563; eff.second=0.089153;
              	} else if( pt < 30.0){ eff.first=0.909366; eff.second=0.062481;
              	} else if( pt < 40.0){ eff.first=0.945833; eff.second=0.021237;
              	} else if( pt < 50.0){ eff.first=0.927066; eff.second=0.009366;
              	} else { eff.first=0.936267; eff.second=0.013983;
              	}
              }else if( eta >= 1.6 && eta < 2.0){
              		if( pt < 20.0){ eff.first=0.796923; eff.second=0.049056;
              	} else if( pt < 30.0){ eff.first=0.852071; eff.second=0.023558;
              	} else if( pt < 40.0){ eff.first=0.914861; eff.second=0.014254;
              	} else if( pt < 50.0){ eff.first=0.938441; eff.second=0.005379;
              	} else { eff.first=0.967742; eff.second=0.010596;
              	}
              }else {
              		if( pt < 20.0){ eff.first=0.829341; eff.second=0.045692;
              	} else if( pt < 30.0){ eff.first=0.880077; eff.second=0.020558;
              	} else if( pt < 40.0){ eff.first=0.916279; eff.second=0.014515;
              	} else if( pt < 50.0){ eff.first=0.942737; eff.second=0.008896;
              	} else { eff.first=0.949735; eff.second=0.013293;
              	}
              }

	  }
        }
        break;
      case 13:
        // taken from https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2#Results_for_2015_data
        switch(cutVersion){
        case patUtils::CutVersion::Spring15Cut25ns :

	  if(wp=="loose"){
                if( Abseta < 0.9 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){ eff.first=0.9855313897132874; eff.second=0.0025857207589580927;
                   }else if( pt < 30.0 ){ eff.first=0.9932104945182800; eff.second=0.0011055782087966714;
                   }else if( pt < 40.0 ){ eff.first=0.9978406429290771; eff.second=0.0003008564730347503;
                   }else if( pt < 50.0 ){ eff.first=0.9986774325370789; eff.second=0.0020699669694531743;
                   }else if( pt < 60.0 ){ eff.first=0.9962952136993408; eff.second=0.0005875918970945443;
                   }else                { eff.first=0.9991089105606079; eff.second=0.0013970545733599567;
                   }
                }else if( Abseta < 1.2 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){ eff.first=0.9967131614685059; eff.second=0.0038540071044748155;
                   }else if( pt < 30.0 ){ eff.first=0.9930336475372314; eff.second=0.002030384666413669;
                   }else if( pt < 40.0 ){ eff.first=0.9970052242279053; eff.second=0.0005914373468172853;
                   }else if( pt < 50.0 ){ eff.first=0.9979928135871887; eff.second=0.00033041698173930095;
                   }else if( pt < 60.0 ){ eff.first=0.9988929033279419; eff.second=0.0009884578713132953;
                   }else                { eff.first=0.9999018907546997; eff.second=0.0025239172940766615;
                   }
                }else if( Abseta < 2.1 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){ eff.first=1.0008281469345093; eff.second=0.00041628854762998813;
                   }else if( pt < 30.0 ){ eff.first=0.9959472417831421; eff.second=0.0011703408377733398;
                   }else if( pt < 40.0 ){ eff.first=0.9986226558685303; eff.second=0.0003983177208424227;
                   }else if( pt < 50.0 ){ eff.first=0.9991121292114258; eff.second=0.00022909925366406707;
                   }else if( pt < 60.0 ){ eff.first=0.9961471557617188; eff.second=0.004028990188838358;
                   }else                { eff.first=1.0007489919662476; eff.second=0.0016929141245890773;
                   }              
                }else if( Abseta < 2.4 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){ eff.first=1.000418782234192; eff.second=0.004228367068434378;
                   }else if( pt < 30.0 ){ eff.first=0.9931564927101135; eff.second=0.002443964514214304;
                   }else if( pt < 40.0 ){ eff.first=0.9958829879760742; eff.second=0.0010036017286035094;
                   }else if( pt < 50.0 ){ eff.first=0.9963051080703735; eff.second=0.0008026731752558733;
                   }else if( pt < 60.0 ){ eff.first=0.9898096323013306; eff.second=0.002625542653614785;
                   }else                { eff.first=1.0001412630081177; eff.second=0.010127272207232552;
                   }             
               }

          }else if(wp=="medium"){
                if( Abseta < 0.9 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9794142842292786; eff.second=0.002741508660046282;
                   }else if( pt < 30.0 ){eff.first=0.9848145842552185; eff.second=0.001290945415659455;
                   }else if( pt < 40.0 ){eff.first=0.9895508289337158; eff.second=0.0004331525546993884;
                   }else if( pt < 50.0 ){eff.first=0.9913147687911987; eff.second=0.00031971235747509097;
                   }else if( pt < 60.0 ){eff.first=0.9875446557998657; eff.second=0.0008303723453658396;
                   }else                {eff.first=0.9913212656974792; eff.second=0.0016589422525152388;
                   }
                }else if( Abseta < 1.2 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9879737496376038; eff.second=0.004043517161003101;
                   }else if( pt < 30.0 ){eff.first=0.9857050776481628; eff.second=0.0022352052935438993;
                   }else if( pt < 40.0 ){eff.first=0.9917845129966736; eff.second=0.0007618769816042719;
                   }else if( pt < 50.0 ){eff.first=0.9923620223999023; eff.second=0.0005198001949915494;
                   }else if( pt < 60.0 ){eff.first=0.9917807579040527; eff.second=0.0013459997816068128;
                   }else                {eff.first=0.9922217726707458; eff.second=0.0028402481587021917;
                   }
                }else if( Abseta < 2.1 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9957935214042664; eff.second=0.0021833565538039194;
                   }else if( pt < 30.0 ){eff.first=0.9913045763969421; eff.second=0.0012557686620512445;
                   }else if( pt < 40.0 ){eff.first=0.9932653903961182; eff.second=0.00047076468154621745;
                   }else if( pt < 50.0 ){eff.first=0.9944870471954346; eff.second=0.0003082350400714759;
                   }else if( pt < 60.0 ){eff.first=0.9912225008010864; eff.second=0.0009381904855814878;
                   }else                {eff.first=0.9955874085426331; eff.second=0.002450189400646501;
                   }
                }else if( Abseta < 2.4 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9767289161682129; eff.second=0.004725463342112523;
                   }else if( pt < 30.0 ){eff.first=0.9677559733390808; eff.second=0.003017063942248565;
                   }else if( pt < 40.0 ){eff.first=0.9679496884346008; eff.second=0.001433179514184055;
                   }else if( pt < 50.0 ){eff.first=0.9651092290878296; eff.second=0.0013090371257890464;
                   }else if( pt < 60.0 ){eff.first=0.9578791856765747; eff.second=0.003330288308070361;
                   }else                {eff.first=0.9747697710990906; eff.second=0.01065126483901303;
                   }
               }

           }else if(wp=="tight"){
                if( Abseta < 0.9 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9752116203308105; eff.second=0.0030660638813280626;
                   }else if( pt < 30.0 ){eff.first=0.9848297238349915; eff.second=0.0016307213764927449;
                   }else if( pt < 40.0 ){eff.first=0.9861794114112854; eff.second=0.0006187187412138267;
                   }else if( pt < 50.0 ){eff.first=0.987443208694458; eff.second=0.000494159746725046;
                   }else if( pt < 60.0 ){eff.first=0.9834294319152832; eff.second=0.0011818999573518245;
                   }else                {eff.first=0.9863178730010986; eff.second=0.002073330940717176;
                   }
                }else if( Abseta < 1.2 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;                      
                   }else if( pt < 25.0 ){eff.first=0.9738101959228516; eff.second=0.004502934246978295;
                   }else if( pt < 30.0 ){eff.first=0.978645384311676; eff.second=0.0027064755458685794;
                   }else if( pt < 40.0 ){eff.first=0.9798933267593384; eff.second=0.001057081371390319;
                   }else if( pt < 50.0 ){eff.first=0.980233907699585; eff.second=0.000819615406448897;
                   }else if( pt < 60.0 ){eff.first=0.9773300886154175; eff.second=0.001955436343316424;
                   }else                {eff.first=0.9795225858688354; eff.second=0.0035622593553725837;
                   }
                }else if( Abseta < 2.1 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9983288645744324; eff.second=0.002331323348626783;
                   }else if( pt < 30.0 ){eff.first=0.9905462265014648; eff.second=0.001402578599690647;
                   }else if( pt < 40.0 ){eff.first=0.9923668503761292; eff.second=0.0005653311393042486;
                   }else if( pt < 50.0 ){eff.first=0.9927627444267273; eff.second=0.0004155573807947332;
                   }else if( pt < 60.0 ){eff.first=0.9886322021484375; eff.second=0.0011254961157344963;
                   }else                {eff.first=0.9950451850891113; eff.second=0.002673833447209764;
                   }
                }else if( Abseta < 2.4 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9877836108207703; eff.second=0.004915740433340289;
                   }else if( pt < 30.0 ){eff.first=0.9802553653717041; eff.second=0.003173276637083633;
                   }else if( pt < 40.0 ){eff.first=0.9785045385360718; eff.second=0.0015542030446523895;
                   }else if( pt < 50.0 ){eff.first=0.9778544902801514; eff.second=0.001456799997296391;
                   }else if( pt < 60.0 ){eff.first=0.9654409885406494; eff.second=0.003709169009223743;
                   }else if( pt < 120.0 ){eff.first=0.9689615368843079; eff.second=0.011084748199568817;
                   }
               }

           }else if(wp=="tightiso"){ //with respect to tight Id
                //coming from Hugues Brun

                if( Abseta < 0.9 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;                            
                   }else if( pt < 25.0 ){ eff.first=1.0043761730194092; eff.second=0.003959090391076143;
                   }else if( pt < 30.0 ){ eff.first=0.9995378255844116; eff.second=0.0022512071035640673;
                   }else if( pt < 40.0 ){ eff.first=1.000901222229004; eff.second=0.0007979481788689052;
                   }else if( pt < 50.0 ){ eff.first=0.9986253976821899; eff.second=0.0004518361024064332;
                   }else if( pt < 60.0 ){ eff.first=1.0002487897872925; eff.second=0.000772847340102783;
                   }else                { eff.first=0.9986850023269653; eff.second=0.0008907575174433545;
                   }
                }else if( Abseta < 1.2 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;                            
                   }else if( pt < 25.0 ){ eff.first=1.004357933998108; eff.second=0.006125539530136138;
                   }else if( pt < 30.0 ){ eff.first=1.002331256866455; eff.second=0.004003683572512011;
                   }else if( pt < 40.0 ){ eff.first=1.004658579826355; eff.second=0.0014502638048416372;
                   }else if( pt < 50.0 ){ eff.first=1.0013608932495117; eff.second=0.0004888604573095644;
                   }else if( pt < 60.0 ){ eff.first=0.9986217021942139; eff.second=0.0012396364566794034;
                   }else                { eff.first=1.0054655075073242; eff.second=0.001589130019220112;
                   }
                }else if( Abseta < 2.1 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;                            
                   }else if( pt < 25.0 ){ eff.first=0.9970762133598328; eff.second=0.003109125287470401;
                   }else if( pt < 30.0 ){ eff.first=1.0006532669067383; eff.second=0.002067755362435184;
                   }else if( pt < 40.0 ){ eff.first=1.0023553371429443; eff.second=0.0008445520691793605;
                   }else if( pt < 50.0 ){ eff.first=0.999933660030365; eff.second=0.0004309914887707696;
                   }else if( pt < 60.0 ){ eff.first=1.0002963542938232; eff.second=0.0007614160360063238;
                   }else                { eff.first=1.0004935264587402; eff.second=0.0009382223143922724;
                   }

                }else if( Abseta < 2.4 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;                            
                   }else if( pt < 25.0 ){ eff.first=0.9957730770111084; eff.second=0.006137193387970902;
                   }else if( pt < 30.0 ){ eff.first=0.9939026832580566; eff.second=0.004261971076013437;
                   }else if( pt < 40.0 ){ eff.first=0.997478187084198; eff.second=0.001781225374381486;
                   }else if( pt < 50.0 ){ eff.first=1.002805233001709; eff.second=0.001100242856214239;
                   }else if( pt < 60.0 ){ eff.first=1.0043764114379883; eff.second=0.001806526581100641;
                   }else                { eff.first=1.0010104179382324; eff.second=0.0022795762936220253;
                   }
               }
           }
	   break;

         case patUtils::CutVersion::ICHEP16Cut :

	  if(wp=="loose"){
              if( Abseta < 0.9){
              		if( pt < 25.0){ eff.first=0.977091; eff.second=0.002205;
              	} else if( pt < 30.0){ eff.first=0.991379; eff.second=0.000923;
              	} else if( pt < 40.0){ eff.first=0.996729; eff.second=0.000253;
              	} else if( pt < 50.0){ eff.first=0.996736; eff.second=0.000194;
              	} else if( pt < 60.0){ eff.first=0.992798; eff.second=0.000631;
              	} else if( pt < 100.0){ eff.first=0.991625; eff.second=0.001656;
              	} else { eff.first=1.000993; eff.second=0.014993;
              	}
              }else if( Abseta < 1.2){
              		if( pt < 25.0){ eff.first=0.995499; eff.second=0.003316;
              	} else if( pt < 30.0){ eff.first=0.992138; eff.second=0.001565;
              	} else if( pt < 40.0){ eff.first=0.996141; eff.second=0.000482;
              	} else if( pt < 50.0){ eff.first=0.995362; eff.second=0.000296;
              	} else if( pt < 60.0){ eff.first=0.996281; eff.second=0.001891;
              	} else if( pt < 100.0){ eff.first=0.991941; eff.second=0.002883;
              	} else { eff.first=1.000000; eff.second=0.006808;
              	}
              }else if( Abseta < 2.1){
              		if( pt < 25.0){ eff.first=1.000469; eff.second=0.001280;
              	} else if( pt < 30.0){ eff.first=0.999010; eff.second=0.000927;
              	} else if( pt < 40.0){ eff.first=0.999120; eff.second=0.000329;
              	} else if( pt < 50.0){ eff.first=0.998493; eff.second=0.000211;
              	} else if( pt < 60.0){ eff.first=0.997555; eff.second=0.000946;
              	} else if( pt < 100.0){ eff.first=0.994120; eff.second=0.002451;
              	} else { eff.first=1.003575; eff.second=0.007436;
              	}
              }else {
              		if( pt < 25.0){ eff.first=1.002567; eff.second=0.000836;
              	} else if( pt < 30.0){ eff.first=0.997522; eff.second=0.001657;
              	} else if( pt < 40.0){ eff.first=0.995137; eff.second=0.000677;
              	} else if( pt < 50.0){ eff.first=0.994269; eff.second=0.000181;
              	} else if( pt < 60.0){ eff.first=0.998696; eff.second=0.002594;
              	} else if( pt < 100.0){ eff.first=0.982880; eff.second=0.006669;
              	} else { eff.first=0.952418; eff.second=0.068688;
              	}
              }
          }else if(wp=="medium"){

              if( Abseta < 0.9){
              		if( pt < 25.0){ eff.first=0.971891; eff.second=0.002295;
              	} else if( pt < 30.0){ eff.first=0.985390; eff.second=0.001010;
              	} else if( pt < 40.0){ eff.first=0.990564; eff.second=0.000320;
              	} else if( pt < 50.0){ eff.first=0.990956; eff.second=0.000139;
              	} else if( pt < 60.0){ eff.first=0.987879; eff.second=0.000662;
              	} else if( pt < 100.0){ eff.first=0.986944; eff.second=0.001750;
              	} else { eff.first=0.996012; eff.second=0.015247;
              	}
              }else if( Abseta < 1.2){
              		if( pt < 25.0){ eff.first=0.988667; eff.second=0.003403;
              	} else if( pt < 30.0){ eff.first=0.986373; eff.second=0.001691;
              	} else if( pt < 40.0){ eff.first=0.991130; eff.second=0.000577;
              	} else if( pt < 50.0){ eff.first=0.990543; eff.second=0.000395;
              	} else if( pt < 60.0){ eff.first=0.990838; eff.second=0.001251;
              	} else if( pt < 100.0){ eff.first=0.987564; eff.second=0.003048;
              	} else { eff.first=1.000000; eff.second=0.009091;
              	}
              }else if( Abseta < 2.1){
              		if( pt < 25.0){ eff.first=0.995891; eff.second=0.001844;
              	} else if( pt < 30.0){ eff.first=0.994749; eff.second=0.000973;
              	} else if( pt < 40.0){ eff.first=0.994626; eff.second=0.000366;
              	} else if( pt < 50.0){ eff.first=0.994182; eff.second=0.000256;
              	} else if( pt < 60.0){ eff.first=0.994234; eff.second=0.001026;
              	} else if( pt < 100.0){ eff.first=0.990053; eff.second=0.002490;
              	} else { eff.first=1.004801; eff.second=0.007700;
              	}
              }else {
              		if( pt < 25.0){ eff.first=0.980939; eff.second=0.003398;
              	} else if( pt < 30.0){ eff.first=0.972391; eff.second=0.002003;
              	} else if( pt < 40.0){ eff.first=0.971297; eff.second=0.000949;
              	} else if( pt < 50.0){ eff.first=0.968591; eff.second=0.000430;
              	} else if( pt < 60.0){ eff.first=0.975689; eff.second=0.004733;
              	} else if( pt < 100.0){ eff.first=0.960116; eff.second=0.007030;
              	} else { eff.first=0.921766; eff.second=0.078150;
              	}
              }

          }else if(wp=="tight"){
              if( Abseta < 0.9){
              		if( pt < 25.0){ eff.first=0.958376; eff.second=0.002567;
              	} else if( pt < 30.0){ eff.first=0.970763; eff.second=0.001265;
              	} else if( pt < 40.0){ eff.first=0.975881; eff.second=0.000481;
              	} else if( pt < 50.0){ eff.first=0.976182; eff.second=0.000387;
              	} else if( pt < 60.0){ eff.first=0.972266; eff.second=0.001000;
              	} else if( pt < 100.0){ eff.first=0.973523; eff.second=0.002016;
              	} else { eff.first=0.987640; eff.second=0.016082;
              	}
              }else if( Abseta < 1.2){
              		if( pt < 25.0){ eff.first=0.969183; eff.second=0.003758;
              	} else if( pt < 30.0){ eff.first=0.966866; eff.second=0.002066;
              	} else if( pt < 40.0){ eff.first=0.971573; eff.second=0.000830;
              	} else if( pt < 50.0){ eff.first=0.970849; eff.second=0.000627;
              	} else if( pt < 60.0){ eff.first=0.971497; eff.second=0.001676;
              	} else if( pt < 100.0){ eff.first=0.968216; eff.second=0.003472;
              	} else { eff.first=1.027140; eff.second=0.022479;
              	}
              }else if( Abseta < 2.1){
              		if( pt < 25.0){ eff.first=0.991418; eff.second=0.001973;
              	} else if( pt < 30.0){ eff.first=0.989645; eff.second=0.001097;
              	} else if( pt < 40.0){ eff.first=0.990198; eff.second=0.000450;
              	} else if( pt < 50.0){ eff.first=0.990848; eff.second=0.000200;
              	} else if( pt < 60.0){ eff.first=0.991589; eff.second=0.001150;
              	} else if( pt < 100.0){ eff.first=0.989786; eff.second=0.002626;
              	} else { eff.first=1.013486; eff.second=0.009446;
              	}
              }else {
              		if( pt < 25.0){ eff.first=0.985294; eff.second=0.003684;
              	} else if( pt < 30.0){ eff.first=0.975853; eff.second=0.002301;
              	} else if( pt < 40.0){ eff.first=0.973918; eff.second=0.001138;
              	} else if( pt < 50.0){ eff.first=0.969894; eff.second=0.001089;
              	} else if( pt < 60.0){ eff.first=0.981064; eff.second=0.003547;
              	} else if( pt < 100.0){ eff.first=0.975218; eff.second=0.007635;
              	} else { eff.first=0.918080; eff.second=0.079635;
              	}
              }

          }else if(wp=="tightiso"){ //with respect to tight Id
              if( Abseta < 0.9){
              		if( pt < 25.0){ eff.first=0.982323; eff.second=0.002986;
              	} else if( pt < 30.0){ eff.first=0.988605; eff.second=0.001663;
              	} else if( pt < 40.0){ eff.first=0.992043; eff.second=0.000581;
              	} else if( pt < 50.0){ eff.first=0.994218; eff.second=0.000328;
              	} else if( pt < 60.0){ eff.first=0.996457; eff.second=0.000556;
              	} else if( pt < 100.0){ eff.first=0.999023; eff.second=0.000704;
              	} else { eff.first=1.000072; eff.second=0.002222;
              	}
              }else if( Abseta < 1.2){
              		if( pt < 25.0){ eff.first=0.986009; eff.second=0.004569;
              	} else if( pt < 30.0){ eff.first=0.994709; eff.second=0.002892;
              	} else if( pt < 40.0){ eff.first=0.998090; eff.second=0.001061;
              	} else if( pt < 50.0){ eff.first=0.997873; eff.second=0.000448;
              	} else if( pt < 60.0){ eff.first=0.999352; eff.second=0.000927;
              	} else if( pt < 100.0){ eff.first=0.999509; eff.second=0.001104;
              	} else { eff.first=0.995709; eff.second=0.004246;
              	}
              }else if( Abseta < 2.1){
              		if( pt < 25.0){ eff.first=0.986859; eff.second=0.002204;
              	} else if( pt < 30.0){ eff.first=0.993867; eff.second=0.001401;
              	} else if( pt < 40.0){ eff.first=0.997524; eff.second=0.000582;
              	} else if( pt < 50.0){ eff.first=0.997950; eff.second=0.000247;
              	} else if( pt < 60.0){ eff.first=0.999048; eff.second=0.000521;
              	} else if( pt < 100.0){ eff.first=0.998433; eff.second=0.000670;
              	} else { eff.first=1.005706; eff.second=0.002814;
              	}
              }else {
              		if( pt < 25.0){ eff.first=0.983397; eff.second=0.003817;
              	} else if( pt < 30.0){ eff.first=0.987529; eff.second=0.002376;
              	} else if( pt < 40.0){ eff.first=0.993539; eff.second=0.001018;
              	} else if( pt < 50.0){ eff.first=0.996392; eff.second=0.000394;
              	} else if( pt < 60.0){ eff.first=0.997946; eff.second=0.001083;
              	} else if( pt < 100.0){ eff.first=1.000734; eff.second=0.001536;
              	} else { eff.first=0.998580; eff.second=0.005873;
              	}
              }

	  }
	  break;
	}
	break;
      }
      return eff;
    }


  std::pair<float,float> getTriggerEfficiencySF(float pt1, float eta1, float pt2, float eta2, int id, bool is2016){
     double etaCen = std::min(fabs(eta1), fabs(eta2));
     double etaFwd = std::max(fabs(eta1), fabs(eta2));
     id   = abs(id);
     std::pair<float,float> eff(1.0,0.00);

     if(id==121){

	if(is2016){
	if(etaCen<0.8){
                 if(etaFwd<0.8){ eff.first=0.989; eff.second=0.004;
           }else if(etaFwd<1.444){ eff.first=0.989; eff.second=0.005;
           }else if(etaFwd<1.566){ eff.first=0.981; eff.second=0.021;
           }else if(etaFwd<2.0){ eff.first=0.982; eff.second=0.010;
           }else if(etaFwd<2.5){ eff.first=0.959; eff.second=0.016;
           }
	}else if(etaCen<1.444){
                 if(etaFwd<0.8){ eff.first=0.989; eff.second=0.005;
           }else if(etaFwd<1.444){ eff.first=0.990; eff.second=0.007;
           }else if(etaFwd<1.566){ eff.first=0.985; eff.second=0.022;
           }else if(etaFwd<2.0){ eff.first=0.982; eff.second=0.009;
           }else if(etaFwd<2.5){ eff.first=0.964; eff.second=0.010;
           }
	}else if(etaCen<1.566){	//This block belongs to transition region 
                 if(etaFwd<0.8){ eff.first=0.981; eff.second=0.021;
           }else if(etaFwd<1.444){ eff.first=0.983; eff.second=0.023;
           }else if(etaFwd<1.566){ eff.first=0.986; eff.second=0.077;
           }else if(etaFwd<2.0){ eff.first=0.983; eff.second=0.028;
           }else if(etaFwd<2.5){ eff.first=0.959; eff.second=0.030;
           }
	}else if(etaCen<2.0){
                 if(etaFwd<0.8){ eff.first=0.982; eff.second=0.010;
           }else if(etaFwd<1.444){ eff.first=0.984; eff.second=0.009;
           }else if(etaFwd<1.566){ eff.first=0.982; eff.second=0.029;
           }else if(etaFwd<2.0){ eff.first=0.977; eff.second=0.011;
           }else if(etaFwd<2.5){ eff.first=0.965; eff.second=0.011;
           }
	}else if(etaCen<2.5){
                 if(etaFwd<0.8){ eff.first=0.959; eff.second=0.016;
           }else if(etaFwd<1.444){ eff.first=0.966; eff.second=0.010;
           }else if(etaFwd<1.566){ eff.first=0.962; eff.second=0.030;
           }else if(etaFwd<2.0){ eff.first=0.965; eff.second=0.011;
           }else if(etaFwd<2.5){ eff.first=0.943; eff.second=0.011;
           }
 	} 

        eff.second = sqrt(eff.second* eff.second + 0.02*0.02); //add 2% syst uncertainty        
	}
	
	else{
	if(etaCen<0.8){
                 if(etaFwd<0.8){ eff.first=0.998; eff.second=0.025;
           }else if(etaFwd<1.444){ eff.first=1.004; eff.second=0.025;
           }else if(etaFwd<1.566){ eff.first=0.996; eff.second=0.071;
           }else if(etaFwd<2.0){ eff.first=1.006; eff.second=0.044;
           }else if(etaFwd<2.5){ eff.first=1.004; eff.second=0.051;
           }
	}else if(etaCen<1.444){
                 if(etaFwd<0.8){ eff.first=0.994; eff.second=0.019;
           }else if(etaFwd<1.444){ eff.first=1.000; eff.second=0.030;
           }else if(etaFwd<1.566){ eff.first=0.987; eff.second=0.091;
           }else if(etaFwd<2.0){ eff.first=1.006; eff.second=0.035;
           }else if(etaFwd<2.5){ eff.first=1.004; eff.second=0.072;
           }
	}else if(etaCen<1.566){	//This block belongs to transition region 
                 if(etaFwd<0.8){ eff.first=1.003; eff.second=0.070;
           }else if(etaFwd<1.444){ eff.first=1.000; eff.second=0.049;
           }else if(etaFwd<1.566){ eff.first=0.973; eff.second=0.111;
           }else if(etaFwd<2.0){ eff.first=1.012; eff.second=0.119;
           }else if(etaFwd<2.5){ eff.first=1.000; eff.second=0.062;
           }
	}else if(etaCen<2.0){
                 if(etaFwd<0.8){ eff.first=1.000; eff.second=0.024;
           }else if(etaFwd<1.444){ eff.first=1.026; eff.second=0.050;
           }else if(etaFwd<1.566){ eff.first=0.995; eff.second=0.077;
           }else if(etaFwd<2.0){ eff.first=1.018; eff.second=0.030;
           }else if(etaFwd<2.5){ eff.first=0.998; eff.second=0.031;
           }
	}else if(etaCen<2.5){
                 if(etaFwd<0.8){ eff.first=0.991; eff.second=0.046;
           }else if(etaFwd<1.444){ eff.first=1.009; eff.second=0.058;
           }else if(etaFwd<1.566){ eff.first=0.994; eff.second=0.057;
           }else if(etaFwd<2.0){ eff.first=1.000; eff.second=0.032;
           }else if(etaFwd<2.5){ eff.first=1.000; eff.second=0.034;
           }
 	}
	}
        
     }else if(id==169){
        if(etaCen<0.9){
                 if(etaFwd<0.9){ eff.first=0.9490; eff.second=0.001;
           }else if(etaFwd<1.2){ eff.first=0.9628; eff.second=0.005;
           }else if(etaFwd<2.1){ eff.first=0.9477; eff.second=0.005;
           }else if(etaFwd<2.4){ eff.first=0.9420; eff.second=0.010;
           }
        }else if(etaCen<1.2){
                 if(etaFwd<1.2){ eff.first=0.9576; eff.second=0.005;
           }else if(etaFwd<2.1){ eff.first=0.9493; eff.second=0.010;
           }else if(etaFwd<2.4){ eff.first=0.9142; eff.second=0.015;          
           }
        }else if(etaCen<2.1){
                 if(etaFwd<2.1){ eff.first=0.9216; eff.second=0.010;
           }else if(etaFwd<2.4){ eff.first=0.9470; eff.second=0.015;
           }
        }else if(etaCen<2.4){
           if(etaFwd<2.4){ eff.first=0.9252; eff.second=0.025;
           }
        }
        eff.second = sqrt(eff.second* eff.second + 0.04*0.04); //add 4% syst uncertainty (conservative)
     }
   return eff;
 }

 private:

};


#endif
