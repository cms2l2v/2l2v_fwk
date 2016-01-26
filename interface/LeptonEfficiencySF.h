#ifndef LeptonEfficiencySF_h
#define LeptonEfficiencySF_h

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
  std::pair<float,float> getLeptonEfficiency(float pt, float eta, int id, std::string wp)
    {
      eta=fabs(eta);
      id=abs(id);
      
      std::pair<float,float> eff(1.0,0.04);
      switch(id){
      case 11:
	{
	  if(wp=="loose")
	    {
		if( eta>= -2.5 && eta<-2.0)
		  {
		    if(pt<30)      { eff.first=0.8775;  eff.second=0.0067; }
		    else if(pt<40) { eff.first=1.0271;  eff.second=0.0045; }
		    else if(pt<50) { eff.first=1.0186;  eff.second=0.0033; }
		    else           { eff.first=1.0114;  eff.second=0.0073; }

		  }
		else if(eta>=-2.0 && eta<-1.566)
		  {
		    if(pt<30)      { eff.first=0.9905;  eff.second=0.0067; }
		    else if(pt<40) { eff.first=1.0024;  eff.second=0.0043; }
		    else if(pt<50) { eff.first=1.0034;  eff.second=0.0032; }
		    else           { eff.first=1.0044;  eff.second=0.0074; }
		  }
		else if(eta>=-1.566 && eta<-1.444)
		  {
		    if(pt<30)      { eff.first=1.0083;  eff.second=0.0185; }
		    else if(pt<40) { eff.first=1.0133;  eff.second=0.0133; }
		    else if(pt<50) { eff.first=0.9938;  eff.second=0.0079; }
		    else           { eff.first=0.9805;  eff.second=0.0129; }
		  }
		else if(eta>=-1.444 && eta<-0.8)
		  {
		    if(pt<30)      { eff.first=1.0199;  eff.second=0.0072; }
		    else if(pt<40) { eff.first=0.9919;  eff.second=0.0026; }
		    else if(pt<50) { eff.first=0.9881;  eff.second=0.0015; }
		    else           { eff.first=0.9798;  eff.second=0.0030; }
		  }
		else if(eta>-0.8 && eta<0)
		  {
		    if(pt<30)      { eff.first=0.9871;  eff.second=0.0069; }
		    else if(pt<40) { eff.first=0.9752;  eff.second=0.0016; }
		    else if(pt<50) { eff.first=0.9775;  eff.second=0.0015; }
		    else           { eff.first=0.9778;  eff.second=0.0030; }
		  }
                else if(eta>0 && eta<0.8)
                  {
                    if(pt<30)      { eff.first=0.7954;  eff.second=0.0043; }
                    else if(pt<40) { eff.first=0.9774;  eff.second=0.0016; }
                    else if(pt<50) { eff.first=0.9786;  eff.second=0.0015; }
                    else           { eff.first=0.9789;  eff.second=0.0030; }
                  }
                else if(eta>0.8 && eta<1.444)
                  {
                    if(pt<30)      { eff.first=1.0146;  eff.second=0.0072; }
                    else if(pt<40) { eff.first=0.9988;  eff.second=0.0016; }
                    else if(pt<50) { eff.first=0.9881;  eff.second=0.0015; }
                    else           { eff.first=0.9861;  eff.second=0.0030; }
                  }
                else if(eta>1.444 && eta<1.556)
                  {
                    if(pt<30)      { eff.first=1.0393;  eff.second=0.0544; }
                    else if(pt<40) { eff.first=0.9985;  eff.second=0.0133; }
                    else if(pt<50) { eff.first=0.9739;  eff.second=0.0079; }
                    else           { eff.first=0.9854;  eff.second=0.0169; }
                  }
                  else if(eta>1.556 && eta<2.0)
                  {
                    if(pt<30)      { eff.first=1.0233;  eff.second=0.0130; }
                    else if(pt<40) { eff.first=0.9917;  eff.second=0.0043; }
                    else if(pt<50) { eff.first=1.0034;  eff.second=0.0032; }
                    else           { eff.first=1.0077;  eff.second=0.0063; }
                  }
                  
                else
                  {
                    if(pt<30)      { eff.first=1.0306;  eff.second=0.0082; }
                    else if(pt<40) { eff.first=1.0246;  eff.second=0.0045; }
                    else if(pt<50) { eff.first=1.0140;  eff.second=0.0042; }
                    else           { eff.first=1.0149;  eff.second=0.0083; }
                  }

	      }
	  if(wp=="medium")
	      {
		if(eta>=-2.5 && eta<-2.0)
		  {
		    if(pt<30)      { eff.first=1.0536;  eff.second=0.0099; }
		    else if(pt<40) { eff.first=1.0123;  eff.second=0.0050; }
		    else if(pt<50) { eff.first=1.0087;  eff.second=0.0045; }
		    else           { eff.first=1.0060;  eff.second=0.0077; }

		  }
		  else if(eta>=-2.0 && eta<-1.566)
		  {
		    if(pt<30)      { eff.first=0.9695;  eff.second=0.0117; }
		    else if(pt<40) { eff.first=0.9713;  eff.second=0.0047; }
		    else if(pt<50) { eff.first=0.9833;  eff.second=0.0033; }
		    else           { eff.first=0.9851;  eff.second=0.0057; }
		  }
		else if(eta>=-1.566 && eta<-1.444)
		  {
		    if(pt<30)      { eff.first=1.0404;  eff.second=0.0231; }
		    else if(pt<40) { eff.first=1.0248;  eff.second=0.0120; }
		    else if(pt<50) { eff.first=0.9902;  eff.second=0.0090; }
		    else           { eff.first=0.9489;  eff.second=0.0143; }
		  }
		else if(eta>=-1.444 && eta<-0.8)
		  {
		    if(pt<30)      { eff.first=1.0594;  eff.second=0.0087; }
		    else if(pt<40) { eff.first=0.9909;  eff.second=0.0029; }
		    else if(pt<50) { eff.first=0.9741;  eff.second=0.0016; }
		    else           { eff.first=0.9581;  eff.second=0.0040; }
		  }
		else if(eta>=-0.8 && eta<0)
		  {
		    if(pt<30)      { eff.first=1.0226;  eff.second=0.0081; }
		    else if(pt<40) { eff.first=0.9786;  eff.second=0.0028; }
		    else if(pt<50) { eff.first=0.9698;  eff.second=0.0016; }
		    else           { eff.first=0.9675;  eff.second=0.0031; }
		  }
		else if(eta>=0 && eta<0.8)
		  {
		    if(pt<30)      { eff.first=1.0180;  eff.second=0.0067; }
		    else if(pt<40) { eff.first=0.9862;  eff.second=0.0028; }
		    else if(pt<50) { eff.first=0.9744;  eff.second=0.0016; }
		    else           { eff.first=0.9708;  eff.second=0.0031; }
		  }
                else if(eta>=0.8 && eta<1.444)
                  {
                    if(pt<30)      { eff.first=1.0774;  eff.second=0.0276; }
                    else if(pt<40) { eff.first=1.0065;  eff.second=0.0029; }
                    else if(pt<50) { eff.first=0.9764;  eff.second=0.0017; }
                    else           { eff.first=0.9648;  eff.second=0.0041; }
                  }
                else if(eta>=1.444 && eta<1.566)
                  {
                    if(pt<30)      { eff.first=1.1880;  eff.second=0.0264; }
                    else if(pt<40) { eff.first=0.9910;  eff.second=0.0154; }
                    else if(pt<50) { eff.first=0.9562;  eff.second=0.0400; }
                    else           { eff.first=0.9744;  eff.second=0.0180; }
                  }
                  else if(eta>=1.566 && eta<2.0)
                  {
                    if(pt<30)      { eff.first=0.9804;  eff.second=0.0073; }
                    else if(pt<40) { eff.first=0.9619;  eff.second=0.0047; }
                    else if(pt<50) { eff.first=0.9785;  eff.second=0.0033; }
                    else           { eff.first=0.9943;  eff.second=0.0057; }
                  }
                else 
                  {
                    if(pt<30)      { eff.first=1.0265;  eff.second=0.0097; }
                    else if(pt<40) { eff.first=1.0096;  eff.second=0.0050; }
                    else if(pt<50) { eff.first=1.0013;  eff.second=0.0045; }
                    else           { eff.first=1.0048;  eff.second=0.0087; }
                  }

	      }
	  if(wp=="tight")
	      {
		if(eta>=-2.5 && eta<-2.0)
		  {
		    if(pt<30)      { eff.first=1.0586;  eff.second=0.0114; }
		    else if(pt<40) { eff.first=1.0246;  eff.second=0.0060; }
		    else if(pt<50) { eff.first=1.0101;  eff.second=0.0052; }
		    else           { eff.first=0.9987;  eff.second=0.0097; }

		  }
		else if(eta>=-2.0 && eta<-1.566)
		  {
		    if(pt<30)      { eff.first=0.9838;  eff.second=0.0128; }
		    else if(pt<40) { eff.first=0.9686;  eff.second=0.0056; }
		    else if(pt<50) { eff.first=0.9819;  eff.second=0.0050; }
		    else           { eff.first=0.9806;  eff.second=0.0075; }
		  }
		  else if(eta>=-1.566 && eta<-1.444)
		  {
		    if(pt<30)      { eff.first=1.1862;  eff.second=0.0498; }
		    else if(pt<40) { eff.first=1.0152;  eff.second=0.0146; }
		    else if(pt<50) { eff.first=0.9934;  eff.second=0.0156; }
		    else           { eff.first=0.9449;  eff.second=0.0176; }
		  }
		else if(eta>=-1.444 && eta<-0.8)
		  {
		    if(pt<30)      { eff.first=1.0502;  eff.second=0.0109; }
		    else if(pt<40) { eff.first=0.9892;  eff.second=0.0034; }
		    else if(pt<50) { eff.first=0.9704;  eff.second=0.0030; }
		    else           { eff.first=0.9562;  eff.second=0.0044; }
		  }
		else if(eta>=-0.8 && eta<0)
		  {
		    if(pt<30)      { eff.first=1.0188;  eff.second=0.0078; }
		    else if(pt<40) { eff.first=0.9644;  eff.second=0.0033; }
		    else if(pt<50) { eff.first=0.9565;  eff.second=0.0018; }
		    else           { eff.first=0.9544;  eff.second=0.0034; }
		  }
		else if(eta>=0 && eta<0.8)
		  {
		    if(pt<30)      { eff.first=0.9944;  eff.second=0.0059; }
		    else if(pt<40) { eff.first=0.9674;  eff.second=0.0021; }
		    else if(pt<50) { eff.first=0.9617;  eff.second=0.0018; }
		    else           { eff.first=0.9617;  eff.second=0.0034; }
		  }
		  else if(eta>=0.8 && eta<1.444)
                  {
                    if(pt<30)      { eff.first=1.0567;  eff.second=0.0148; }
                    else if(pt<40) { eff.first=0.9969;  eff.second=0.0035; }
                    else if(pt<50) { eff.first=0.9730;  eff.second=0.0030; }
                    else           { eff.first=0.9561;  eff.second=0.0045; }
                  }
                else if(eta>=1.444 && eta<1.556)
                  {
                    if(pt<30)      { eff.first=1.0601;  eff.second=0.0439; }
                    else if(pt<40) { eff.first=0.9868;  eff.second=0.1167; }
                    else if(pt<50) { eff.first=0.9568;  eff.second=0.0110; }
                    else           { eff.first=0.9493;  eff.second=0.0164; }
                  }
                else if(eta>=1.556 && eta<2.0)
                  {
                    if(pt<30)      { eff.first=1.0225;  eff.second=0.0130; }
                    else if(pt<40) { eff.first=0.9670;  eff.second=0.0056; }
                    else if(pt<50) { eff.first=0.9847;  eff.second=0.0050; }
                    else           { eff.first=0.9922;  eff.second=0.0075; }
                  }
                else
                  {
                    if(pt<30)      { eff.first=1.0560;  eff.second=0.0113; }
                    else if(pt<40) { eff.first=1.0229;  eff.second=0.0459; }
                    else if(pt<50) { eff.first=1.0043;  eff.second=0.0052; }
                    else           { eff.first=1.0013;  eff.second=0.0097; }
                  }
	      }
	}
	break;
      case 13:
	{
	  if(wp=="loose")
	    {
	      if(eta>=-2.4 && eta<-2.1){
		if(pt<20) { eff.first=1.; eff.second=1.; }
		else if(pt<25) { eff.first=1.07601; eff.second=0.16711; }
		else if(pt<30) { eff.first=1.09843; eff.second=0.12987; }
		else if(pt<35) { eff.first=1.02712; eff.second=0.10229; }
		else if(pt<40) { eff.first=0.99682; eff.second=0.09248; }
		else if(pt<50) { eff.first=1.01843; eff.second=0.06904; }
		else if(pt<60) { eff.first=1.01043; eff.second=0.12138; }
		else if(pt<90) { eff.first=1.00829; eff.second=0.19043; }
		else if(pt<140) { eff.first=1.0119; eff.second=0.44558; }
		else { eff.first=0.919903; eff.second=0.48963; }
	      }
	      else if(eta>=-2.1 && eta<-1.2){
		if(pt<20) { eff.first=1.; eff.second=1.; }
		else if(pt<25) { eff.first=1.04955; eff.second=0.08135; }
		else if(pt<30) { eff.first=1.02081; eff.second=0.06493; }
		else if(pt<35) { eff.first=1.01885; eff.second=0.05321; }
		else if(pt<40) { eff.first=1.01087; eff.second=0.04669; }
		else if(pt<50) { eff.first=1.01436; eff.second=0.03044; }
		else if(pt<60) { eff.first=1.00420; eff.second=0.06407; }
		else if(pt<90) { eff.first=1.00640; eff.second=0.09006; }
		else if(pt<140) { eff.first=0.96182; eff.second=0.19338; }
		else { eff.first=0.91990; eff.second=0.48963; }
	      }
	      else if(eta>=-1.2 && eta<-0.9){
		if(pt<20) { eff.first=1.; eff.second=1.; }
		else if(pt<25) { eff.first=0.91527; eff.second=0.12052; }
		else if(pt<30) { eff.first=0.99605; eff.second=0.11257; }
		else if(pt<35) { eff.first=1.02798; eff.second=0.08473; }
		else if(pt<40) { eff.first=1.00918; eff.second=0.06603; }
		else if(pt<50) { eff.first=1.00147; eff.second=0.04819; }
		else if(pt<60) { eff.first=1.00580; eff.second=0.09549; }
		else if(pt<90) { eff.first=0.99196; eff.second=0.14038; }
		else if(pt<140) { eff.first=0.97827; eff.second=0.30764; }
		else { eff.first=1.0444; eff.second=1.49927; }
	      }
	      else if(eta>=-0.9 && eta<0){
                if(pt<20) { eff.first=1.; eff.second=1.; }
                else if(pt<25) { eff.first=1.00446; eff.second=0.07225; }
                else if(pt<30) { eff.first=1.05353; eff.second=0.05499; }
                else if(pt<35) { eff.first=1.02960; eff.second=0.04229; }
                else if(pt<40) { eff.first=1.00097; eff.second=0.03549; }
                else if(pt<50) { eff.first=1.01108; eff.second=0.02664; }
                else if(pt<60) { eff.first=1.01003; eff.second=0.05316; }
                else if(pt<90) { eff.first=0.99317; eff.second=0.07422; }
                else if(pt<140) { eff.first=0.99001; eff.second=0.18719; }
                else { eff.first=1.03346; eff.second=0.43489; }
	      }
              else if(eta>=0 && eta<0.9){
                if(pt<20) { eff.first=1.; eff.second=1.; }
                else if(pt<25) { eff.first=1.02303; eff.second=0.07427; }
                else if(pt<30) { eff.first=0.99755; eff.second=0.05357; }
                else if(pt<35) { eff.first=1.01741; eff.second=0.04169; }
                else if(pt<40) { eff.first=1.02791; eff.second=0.03655; }
                else if(pt<50) { eff.first=1.01122; eff.second=0.02686; }
                else if(pt<60) { eff.first=1.00585; eff.second=0.05625; }
                else if(pt<90) { eff.first=0.99676; eff.second=0.07499; }
                else if(pt<140) { eff.first=0.99643; eff.second=0.17439; }
                else { eff.first=1.03247; eff.second=0.38604; }
              }
              else if(eta>0.9 && eta<1.2){
                if(pt<20) { eff.first=1.; eff.second=1.; }
                else if(pt<25) { eff.first=1.01873; eff.second=0.13685; }
                else if(pt<30) { eff.first=0.98156; eff.second=0.10515; }
                else if(pt<35) { eff.first=1.03577; eff.second=0.08535; }
                else if(pt<40) { eff.first=1.01255; eff.second=0.06999; }
                else if(pt<50) { eff.first=1.01294; eff.second=0.04773; }
                else if(pt<60) { eff.first=0.99653; eff.second=0.09673; }
                else if(pt<90) { eff.first=1.00852; eff.second=0.140352; }
                else if(pt<140) { eff.first=0.93710; eff.second=0.38227; }
                else { eff.first=1.01719; eff.second=0.66627; }
              }
              else if(eta>=1.2 && eta<2.1){
                if(pt<20) { eff.first=1.; eff.second=1.; }
                else if(pt<25) { eff.first=1.05358; eff.second=0.08379; }
                else if(pt<30) { eff.first=1.02293; eff.second=0.06561; }
                else if(pt<35) { eff.first=1.03277; eff.second=0.05539; }
                else if(pt<40) { eff.first=1.00034; eff.second=0.04529; }
                else if(pt<50) { eff.first=1.01050; eff.second=0.03045; }
                else if(pt<60) { eff.first=1.00007; eff.second=0.06254; }
                else if(pt<90) { eff.first=1.00413; eff.second=0.09162; }
                else if(pt<140) { eff.first=1.01646; eff.second=0.22425; }
                else { eff.first=0.92201; eff.second=0.43435; }
              }
              else{
                if(pt<20) { eff.first=1.; eff.second=1.; }
                else if(pt<25) { eff.first=1.06812; eff.second=0.15351; }
                else if(pt<30) { eff.first=1.03892; eff.second=0.12380; }
                else if(pt<35) { eff.first=1.07406; eff.second=0.11594; }
                else if(pt<40) { eff.first=1.00818; eff.second=0.09379; }
                else if(pt<50) { eff.first=1.01656; eff.second=0.06991; }
                else if(pt<60) { eff.first=0.98877; eff.second=0.14096; }
                else if(pt<90) { eff.first=0.98193; eff.second=0.18826; }
                else if(pt<140) { eff.first=1.02903; eff.second=0.47887; }
                else { eff.first=1.06485; eff.second=1.55669; }
              }
	    }

          if(wp=="medium")
           {
            if(eta>=-2.4 && eta<-2.1){
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.08065; eff.second=0.16783; }
              else if(pt<30) { eff.first=1.10371; eff.second=0.13051; }
              else if(pt<35) { eff.first=1.02229; eff.second=0.10207; }
              else if(pt<40) { eff.first=0.99576; eff.second=0.09261; }
              else if(pt<50) { eff.first=1.02214; eff.second=0.06937; }
              else if(pt<60) { eff.first=1.01638; eff.second=0.12209; }
              else if(pt<90) { eff.first=0.97848; eff.second=0.18639; }
              else if(pt<140) { eff.first=1.01412; eff.second=0.44659; }
              else { eff.first=0.92259; eff.second=0.49109; }
            }
            else if(eta>=-2.1 && eta<-1.2){
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.05110; eff.second=0.08148; }
              else if(pt<30) { eff.first=1.02055; eff.second=0.06495; }
              else if(pt<35) { eff.first=1.01467; eff.second=0.05306; }
              else if(pt<40) { eff.first=1.00707; eff.second=0.04658; }
              else if(pt<50) { eff.first=1.01100; eff.second=0.03037; }
              else if(pt<60) { eff.first=1.00401; eff.second=0.06409; }
              else if(pt<90) { eff.first=1.00853; eff.second=0.09025; }
              else if(pt<140) { eff.first=0.94422; eff.second=0.19076; }
              else { eff.first=0.92259; eff.second=0.49109; }
            }
            else if(eta>=-1.2 && eta<-0.9){
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=0.91035; eff.second=0.12021; }
              else if(pt<30) { eff.first=0.99329; eff.second=0.11246; }
              else if(pt<35) { eff.first=1.02539; eff.second=0.08468; }
              else if(pt<40) { eff.first=1.00856; eff.second=0.06606; }
              else if(pt<50) { eff.first=0.99907; eff.second=0.04813; }
              else if(pt<60) { eff.first=1.00121; eff.second=0.09528; }
              else if(pt<90) { eff.first=0.984833; eff.second=0.13971; }
              else if(pt<140) { eff.first=0.984494; eff.second=0.30964; }
              else { eff.first=1.04289; eff.second=1.49702; }
            }
            else if(eta>=-0.9 && eta<0.){
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.00739; eff.second=0.07252; }
              else if(pt<30) { eff.first=1.05238; eff.second=0.05504; }
              else if(pt<35) { eff.first=1.02493; eff.second=0.04220; }
              else if(pt<40) { eff.first=0.99887; eff.second=0.03549; }
              else if(pt<50) { eff.first=1.01009; eff.second=0.02666; }
              else if(pt<60) { eff.first=1.01087; eff.second=0.05328; }
              else if(pt<90) { eff.first=0.99322; eff.second=0.07433; }
              else if(pt<140) { eff.first=0.99655; eff.second=0.18844; }
              else { eff.first=1.04556; eff.second=0.44009; }
            }
            else if(eta>=0. && eta<0.9){
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.02619; eff.second=0.07457; }
              else if(pt<30) { eff.first=0.99752; eff.second=0.05366; }
              else if(pt<35) { eff.first=1.01785; eff.second=0.04177; }
              else if(pt<40) { eff.first=1.03051; eff.second=0.03667; }
              else if(pt<50) { eff.first=1.00945; eff.second=0.02686; }
              else if(pt<60) { eff.first=1.00479; eff.second=0.05629; }
              else if(pt<90) { eff.first=0.99433; eff.second=0.07497; }
              else if(pt<140) { eff.first=0.97057; eff.second=0.17111; }
              else { eff.first=1.04192; eff.second=0.38963; }
            }
            else if(eta>=0.9 && eta<1.2){
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.02499; eff.second=0.13769; }
              else if(pt<30) { eff.first=0.98629; eff.second=0.10566; }
              else if(pt<35) { eff.first=1.03725; eff.second=0.08555; }
              else if(pt<40) { eff.first=1.01218; eff.second=0.07006; }
              else if(pt<50) { eff.first=1.00944; eff.second=0.04766; }
              else if(pt<60) { eff.first=0.98667; eff.second=0.09611; }
              else if(pt<90) { eff.first=1.01357; eff.second=0.14106; }
              else if(pt<140) { eff.first=0.94121; eff.second=0.38397; }
              else { eff.first=1.04371; eff.second=0.68413; }
            }
            else if(eta>=1.2 && eta<2.1){
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.05232; eff.second=0.08376; }
              else if(pt<30) { eff.first=1.02059; eff.second=0.06554; }
              else if(pt<35) { eff.first=1.03308; eff.second=0.05543; }
              else if(pt<40) { eff.first=0.99802; eff.second=0.04524; }
              else if(pt<50) { eff.first=1.00905; eff.second=0.03043; }
              else if(pt<60) { eff.first=0.99832; eff.second=0.06249; }
              else if(pt<90) { eff.first=1.00147; eff.second=0.09146; }
              else if(pt<140) { eff.first=1.01792; eff.second=0.22458; }
              else { eff.first=0.92642; eff.second=0.43647; }
            }
            else{
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.05854; eff.second=0.15258; }
              else if(pt<30) { eff.first=1.03646; eff.second=0.12376; }
              else if(pt<35) { eff.first=1.08026; eff.second=0.11661; }
              else if(pt<40) { eff.first=1.01182; eff.second=0.09426; }
              else if(pt<50) { eff.first=1.02050; eff.second=0.07022; }
              else if(pt<60) { eff.first=0.98515; eff.second=0.14080; }
              else if(pt<90) { eff.first=0.98766; eff.second=0.18937; }
              else if(pt<140) { eff.first=1.03328; eff.second=0.4809; }
              else { eff.first=1.12581; eff.second=1.64865; }
            }
          }
	  if(wp=="tight")
           {
	    if(eta>=-2.4 && eta<-2.1){
	      if(pt<20) { eff.first=1.; eff.second=1.; }
	      else if(pt<25) { eff.first=1.06252; eff.second=0.16618; }
	      else if(pt<30) { eff.first=1.07752; eff.second=0.12862; }
	      else if(pt<35) { eff.first=1.01529; eff.second=0.10192; }
	      else if(pt<40) { eff.first=0.97559; eff.second=0.09148; }
	      else if(pt<50) { eff.first=1.00140; eff.second=0.06849; }
	      else if(pt<60) { eff.first=1.03034; eff.second=0.12403; }
	      else if(pt<90) { eff.first=1.01559; eff.second=0.19266; }
	      else if(pt<140) { eff.first=0.97144; eff.second=0.43821; }
	      else { eff.first=0.93187; eff.second=0.49604; }
	    }
	    else if(eta>=-2.1 && eta<-1.2){
	      if(pt<20) { eff.first=1.; eff.second=1.; }
	      else if(pt<25) { eff.first=1.03835; eff.second=0.08082; }
	      else if(pt<30) { eff.first=1.01140; eff.second=0.06458; }
	      else if(pt<35) { eff.first=1.01547; eff.second=0.05316; }
	      else if(pt<40) { eff.first=1.00008; eff.second=0.04639; }
	      else if(pt<50) { eff.first=1.00883; eff.second=0.03036; }
	      else if(pt<60) { eff.first=0.99369; eff.second=0.06372; }
	      else if(pt<90) { eff.first=1.00013; eff.second=0.08985; }
	      else if(pt<140) { eff.first=0.96954; eff.second=0.19494; }
	      else { eff.first=0.93187; eff.second=0.49605; }
	    }
	    else if(eta>=-1.2 && eta<-0.9){
	      if(pt<20) { eff.first=1.; eff.second=1.; }
	      else if(pt<25) { eff.first=0.91631; eff.second=0.12100; }
	      else if(pt<30) { eff.first=0.98860; eff.second=0.11235; }
	      else if(pt<35) { eff.first=1.02415; eff.second=0.08481; }
	      else if(pt<40) { eff.first=1.00118; eff.second=0.06588; }
	      else if(pt<50) { eff.first=0.99158; eff.second=0.04798; }
	      else if(pt<60) { eff.first=0.99014; eff.second=0.09476; }
	      else if(pt<90) { eff.first=0.99166; eff.second=0.14106; }
	      else if(pt<140) { eff.first=0.99929; eff.second=0.31436; }
	      else { eff.first=1.07497; eff.second=0.45277; }
	    }
            else if(eta>=-0.9 && eta<0){
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.00811; eff.second=0.07281; }
              else if(pt<30) { eff.first=1.04239; eff.second=0.05484; }
              else if(pt<35) { eff.first=1.03106; eff.second=0.04267; }
              else if(pt<40) { eff.first=1.00223; eff.second=0.03574; }
              else if(pt<50) { eff.first=1.02094; eff.second=0.02699; }
              else if(pt<60) { eff.first=1.01323; eff.second=0.05359; }
              else if(pt<90) { eff.first=0.99029; eff.second=0.07443; }
              else if(pt<140) { eff.first=0.99089; eff.second=0.18818; }
              else { eff.first=1.07497; eff.second=0.45277; }
            }
            else if(eta>=0 && eta<0.9){
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.03552; eff.second=0.07537; }
              else if(pt<30) { eff.first=0.99649; eff.second=0.05383; }
              else if(pt<35) { eff.first=1.01973; eff.second=0.04201; }
              else if(pt<40) { eff.first=1.03240; eff.second=0.03687; }
              else if(pt<50) { eff.first=1.01339; eff.second=0.02707; }
              else if(pt<60) { eff.first=1.01172; eff.second=0.05682; }
              else if(pt<90) { eff.first=0.99082; eff.second=0.07502; }
              else if(pt<140) { eff.first=0.97279; eff.second=0.17217; }
              else { eff.first=1.06743; eff.second=0.39929; }
            }
            else if(eta>=0.9 && eta<1.2){
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.02001; eff.second=0.13739; }
              else if(pt<30) { eff.first=0.96877; eff.second=0.10449; }
              else if(pt<35) { eff.first=1.02938; eff.second=0.08529; }
              else if(pt<40) { eff.first=0.98874; eff.second=0.06905; }
              else if(pt<50) { eff.first=1.00059; eff.second=0.04748; }
              else if(pt<60) { eff.first=0.95079; eff.second=0.09375; }
              else if(pt<90) { eff.first=1.02529; eff.second=0.14271; }
              else if(pt<140) { eff.first=0.86955; eff.second=0.36273; }
              else { eff.first=1.03529; eff.second=0.6782; }
            }
            else if(eta>=1.2 && eta<2.1){
              if(pt<20) { eff.first=1.03871; eff.second=0.08304; }
              else if(pt<25) { eff.first=1.02465; eff.second=0.06579; }
              else if(pt<30) { eff.first=1.03115; eff.second=0.05541; }
              else if(pt<35) { eff.first=0.99583; eff.second=0.04522; }
              else if(pt<40) { eff.first=1.00707; eff.second=0.03043; }
              else if(pt<50) { eff.first=0.98830; eff.second=0.06214; }
              else if(pt<60) { eff.first=0.99667; eff.second=0.09130; }
              else if(pt<90) { eff.first=0.93281; eff.second=0.21080; }
              else if(pt<140) { eff.first=0.83726; eff.second=0.40669; }
              else { eff.first=1.05929; eff.second=0.15316; }
            }
            else {
              if(pt<20) { eff.first=1.; eff.second=1.; }
              else if(pt<25) { eff.first=1.05929; eff.second=0.15316; }
              else if(pt<30) { eff.first=1.00336; eff.second=0.12103; }
              else if(pt<35) { eff.first=1.08700; eff.second=0.11752; }
              else if(pt<40) { eff.first=0.98945; eff.second=0.09291; }
              else if(pt<50) { eff.first=1.00016; eff.second=0.06941; }
              else if(pt<60) { eff.first=0.97174; eff.second=0.13999; }
              else if(pt<90) { eff.first=0.84562; eff.second=0.16955; }
              else if(pt<140) { eff.first=0.74416; eff.second=0.37917; }
              else { eff.first=1.09096; eff.second=1.5959; }
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
