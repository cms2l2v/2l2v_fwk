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
		    if(pt<30)      { eff.first=0.94346;  eff.second=0.08928; }
		    else if(pt<40) { eff.first=0.96368;  eff.second=0.0635; }
		    else if(pt<50) { eff.first=0.99027;  eff.second=0.06856; }
		    else           { eff.first=0.92276;  eff.second=0.11803; }

		  }
		else if(eta>=-2.0 && eta<-1.4)
		  {
		    if(pt<30)      { eff.first=0.95742;  eff.second=0.07349; }
		    else if(pt<40) { eff.first=0.98631;  eff.second=0.05149; }
		    else if(pt<50) { eff.first=0.96946;  eff.second=0.05609; }
		    else           { eff.first=0.96212;  eff.second=0.09426; }
		  }
		else if(eta>=-1.4 && eta<-0.8)
		  {
		    if(pt<30)      { eff.first=0.97916;  eff.second=0.06576; }
		    else if(pt<40) { eff.first=0.98246;  eff.second=0.04341; }
		    else if(pt<50) { eff.first=0.98483;  eff.second=0.04064; }
		    else           { eff.first=0.98467;  eff.second=0.07062; }
		  }
		else if(eta>=-0.8 && eta<0)
		  {
		    if(pt<30)      { eff.first=0.98226;  eff.second=0.05211; }
		    else if(pt<40) { eff.first=0.98214;  eff.second=0.0338; }
		    else if(pt<50) { eff.first=0.98527;  eff.second=0.03314; }
		    else           { eff.first=0.98424;  eff.second=0.05669; }
		  }
		else if(eta>0 && eta<0.8)
		  {
		    if(pt<30)      { eff.first=0.97089;  eff.second=0.05216; }
		    else if(pt<40) { eff.first=0.99258;  eff.second=0.03381; }
		    else if(pt<50) { eff.first=0.97086;  eff.second=0.03177; }
		    else           { eff.first=0.97901;  eff.second=0.05379; }
		  }
                else if(eta>0.8 && eta<1.4)
                  {
                    if(pt<30)      { eff.first=0.95247;  eff.second=0.06332; }
                    else if(pt<40) { eff.first=0.97602;  eff.second=0.04223; }
                    else if(pt<50) { eff.first=0.99021;  eff.second=0.04053; }
                    else           { eff.first=0.97189;  eff.second=0.06631; }
                  }
                else if(eta>1.4 && eta<2.0)
                  {
                    if(pt<30)      { eff.first=0.93671;  eff.second=0.07281; }
                    else if(pt<40) { eff.first=0.98022;  eff.second=0.05307; }
                    else if(pt<50) { eff.first=1.00077;  eff.second=0.05669; }
                    else           { eff.first=1.00785;  eff.second=0.09184; }
                  }
                else
                  {
                    if(pt<30)      { eff.first=0.98114;  eff.second=0.07752; }
                    else if(pt<40) { eff.first=0.99129;  eff.second=0.05948; }
                    else if(pt<50) { eff.first=1.00664;  eff.second=0.05948; }
                    else           { eff.first=1.00287;  eff.second=0.11191; }
                  }

	      }
	  if(wp=="medium")
	      {
		if(eta>=-2.5 && eta<-2.0)
		  {
		    if(pt<30)      { eff.first=0.84581;  eff.second=0.08471; }
		    else if(pt<40) { eff.first=0.89007;  eff.second=0.06096; }
		    else if(pt<50) { eff.first=0.92463;  eff.second=0.06599; }
		    else           { eff.first=0.88587;  eff.second=0.11588; }

		  }
		else if(eta>=-2.0 && eta<-1.4)
		  {
		    if(pt<30)      { eff.first=0.89439;  eff.second=0.07229; }
		    else if(pt<40) { eff.first=0.96171;  eff.second=0.05152; }
		    else if(pt<50) { eff.first=0.94212;  eff.second=0.05558; }
		    else           { eff.first=0.91704;  eff.second=0.09181; }
		  }
		else if(eta>=-1.4 && eta<-0.8)
		  {
		    if(pt<30)      { eff.first=0.96571;  eff.second=0.06729; }
		    else if(pt<40) { eff.first=0.94892;  eff.second=0.04324; }
		    else if(pt<50) { eff.first=0.95404;  eff.second=0.04033; }
		    else           { eff.first=0.94187;  eff.second=0.06918; }
		  }
		else if(eta>=-0.8 && eta<0)
		  {
		    if(pt<30)      { eff.first=0.99296;  eff.second=0.05415; }
		    else if(pt<40) { eff.first=0.97392;  eff.second=0.03435; }
		    else if(pt<50) { eff.first=0.97599;  eff.second=0.03342; }
		    else           { eff.first=0.98989;  eff.second=0.05762; }
		  }
		else if(eta>=0 && eta<0.8)
		  {
		    if(pt<30)      { eff.first=0.97001;  eff.second=0.05372; }
		    else if(pt<40) { eff.first=0.98839;  eff.second=0.03443; }
		    else if(pt<50) { eff.first=0.95987;  eff.second=0.03199; }
		    else           { eff.first=0.95920;  eff.second=0.05368; }
		  }
                else if(eta>=0.8 && eta<1.4)
                  {
                    if(pt<30)      { eff.first=0.93096;  eff.second=0.06442; }
                    else if(pt<40) { eff.first=0.95574;  eff.second=0.04254; }
                    else if(pt<50) { eff.first=0.96301;  eff.second=0.04038; }
                    else           { eff.first=0.96160;  eff.second=0.06669; }
                  }
                else if(eta>=1.4 && eta<2.0)
                  {
                    if(pt<30)      { eff.first=0.93805;  eff.second=0.07532; }
                    else if(pt<40) { eff.first=0.96064;  eff.second=0.05330; }
                    else if(pt<50) { eff.first=0.98193;  eff.second=0.05658; }
                    else           { eff.first=0.98035;  eff.second=0.09074; }
                  }
                else 
                  {
                    if(pt<30)      { eff.first=0.95528;  eff.second=0.07815; }
                    else if(pt<40) { eff.first=0.96873;  eff.second=0.05945; }
                    else if(pt<50) { eff.first=0.97990;  eff.second=0.06492; }
                    else           { eff.first=0.97498;  eff.second=0.11061; }
                  }

	      }
	  if(wp=="tight")
	      {
		if(eta>=-2.5 && eta<-2.0)
		  {
		    if(pt<30)      { eff.first=0.83898;  eff.second=0.08719; }
		    else if(pt<40) { eff.first=0.84466;  eff.second=0.05998; }
		    else if(pt<50) { eff.first=0.89877;  eff.second=0.06556; }
		    else           { eff.first=0.83868;  eff.second=0.11265; }

		  }
		else if(eta>=-2.0 && eta<-1.4)
		  {
		    if(pt<30)      { eff.first=0.83047;  eff.second=0.07218; }
		    else if(pt<40) { eff.first=0.96514;  eff.second=0.05336; }
		    else if(pt<50) { eff.first=0.89141;  eff.second=0.05452; }
		    else           { eff.first=0.87330;  eff.second=0.08977; }
		  }
		else if(eta>=-1.4 && eta<-0.8)
		  {
		    if(pt<30)      { eff.first=0.93803;  eff.second=0.06887; }
		    else if(pt<40) { eff.first=0.91022;  eff.second=0.04327; }
		    else if(pt<50) { eff.first=0.92954;  eff.second=0.04048; }
		    else           { eff.first=0.90445;  eff.second=0.06822; }
		  }
		else if(eta>=-0.8 && eta<0)
		  {
		    if(pt<30)      { eff.first=0.96447;  eff.second=0.05418; }
		    else if(pt<40) { eff.first=0.95349;  eff.second=0.03441; }
		    else if(pt<50) { eff.first=0.96792;  eff.second=0.03369; }
		    else           { eff.first=0.97049;  eff.second=0.05739; }
		  }
		else if(eta>=0 && eta<0.8)
		  {
		    if(pt<30)      { eff.first=0.97272;  eff.second=0.05497; }
		    else if(pt<40) { eff.first=0.96903;  eff.second=0.03456; }
		    else if(pt<50) { eff.first=0.95245;  eff.second=0.03226; }
		    else           { eff.first=0.93636;  eff.second=0.05335; }
		  }
                else if(eta>=0.8 && eta<1.4)
                  {
                    if(pt<30)      { eff.first=0.93304;  eff.second=0.06752; }
                    else if(pt<40) { eff.first=0.94285;  eff.second=0.04336; }
                    else if(pt<50) { eff.first=0.92849;  eff.second=0.04019; }
                    else           { eff.first=0.89527;  eff.second=0.06418; }
                  }
                else if(eta>=1.4 && eta<2.0)
                  {
                    if(pt<30)      { eff.first=0.92503;  eff.second=0.07850; }
                    else if(pt<40) { eff.first=0.91792;  eff.second=0.05323; }
                    else if(pt<50) { eff.first=0.93845;  eff.second=0.05589; }
                    else           { eff.first=0.96444;  eff.second=0.91054; }
                  }
                else
                  {
                    if(pt<30)      { eff.first=0.93980;  eff.second=0.07994; }
                    else if(pt<40) { eff.first=0.98190;  eff.second=0.06128; }
                    else if(pt<50) { eff.first=0.95518;  eff.second=0.06467; }
                    else           { eff.first=0.96954;  eff.second=0.11139; }
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
