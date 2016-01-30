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
  std::pair<float,float> getLeptonEfficiency(float pt, float eta, int id, std::string wp){
      eta=fabs(eta);
      id=abs(id);
      
      std::pair<float,float> eff(1.0,0.04);
      switch(id){
      case 11:
	{
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
        }
        break;
      case 13:
        // taken from https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2#Results_for_2015_data
	{
	  if(wp=="loose"){
                if( eta < 0.9 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){ eff.first=0.9855313897132874; eff.second=0.0025857207589580927;
                   }else if( pt < 30.0 ){ eff.first=0.9932104945182800; eff.second=0.0011055782087966714;
                   }else if( pt < 40.0 ){ eff.first=0.9978406429290771; eff.second=0.0003008564730347503;
                   }else if( pt < 50.0 ){ eff.first=0.9986774325370789; eff.second=0.0020699669694531743;
                   }else if( pt < 60.0 ){ eff.first=0.9962952136993408; eff.second=0.0005875918970945443;
                   }else                { eff.first=0.9991089105606079; eff.second=0.0013970545733599567;
                   }
                }else if( eta < 1.2 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){ eff.first=0.9967131614685059; eff.second=0.0038540071044748155;
                   }else if( pt < 30.0 ){ eff.first=0.9930336475372314; eff.second=0.002030384666413669;
                   }else if( pt < 40.0 ){ eff.first=0.9970052242279053; eff.second=0.0005914373468172853;
                   }else if( pt < 50.0 ){ eff.first=0.9979928135871887; eff.second=0.00033041698173930095;
                   }else if( pt < 60.0 ){ eff.first=0.9988929033279419; eff.second=0.0009884578713132953;
                   }else                { eff.first=0.9999018907546997; eff.second=0.0025239172940766615;
                   }
                }else if( eta < 2.1 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){ eff.first=1.0008281469345093; eff.second=0.00041628854762998813;
                   }else if( pt < 30.0 ){ eff.first=0.9959472417831421; eff.second=0.0011703408377733398;
                   }else if( pt < 40.0 ){ eff.first=0.9986226558685303; eff.second=0.0003983177208424227;
                   }else if( pt < 50.0 ){ eff.first=0.9991121292114258; eff.second=0.00022909925366406707;
                   }else if( pt < 60.0 ){ eff.first=0.9961471557617188; eff.second=0.004028990188838358;
                   }else                { eff.first=1.0007489919662476; eff.second=0.0016929141245890773;
                   }              
                }else if( eta < 2.4 ){
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
                if( eta < 0.9 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9794142842292786; eff.second=0.002741508660046282;
                   }else if( pt < 30.0 ){eff.first=0.9848145842552185; eff.second=0.001290945415659455;
                   }else if( pt < 40.0 ){eff.first=0.9895508289337158; eff.second=0.0004331525546993884;
                   }else if( pt < 50.0 ){eff.first=0.9913147687911987; eff.second=0.00031971235747509097;
                   }else if( pt < 60.0 ){eff.first=0.9875446557998657; eff.second=0.0008303723453658396;
                   }else                {eff.first=0.9913212656974792; eff.second=0.0016589422525152388;
                   }
                }else if( eta < 1.2 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9879737496376038; eff.second=0.004043517161003101;
                   }else if( pt < 30.0 ){eff.first=0.9857050776481628; eff.second=0.0022352052935438993;
                   }else if( pt < 40.0 ){eff.first=0.9917845129966736; eff.second=0.0007618769816042719;
                   }else if( pt < 50.0 ){eff.first=0.9923620223999023; eff.second=0.0005198001949915494;
                   }else if( pt < 60.0 ){eff.first=0.9917807579040527; eff.second=0.0013459997816068128;
                   }else                {eff.first=0.9922217726707458; eff.second=0.0028402481587021917;
                   }
                }else if( eta < 2.1 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9957935214042664; eff.second=0.0021833565538039194;
                   }else if( pt < 30.0 ){eff.first=0.9913045763969421; eff.second=0.0012557686620512445;
                   }else if( pt < 40.0 ){eff.first=0.9932653903961182; eff.second=0.00047076468154621745;
                   }else if( pt < 50.0 ){eff.first=0.9944870471954346; eff.second=0.0003082350400714759;
                   }else if( pt < 60.0 ){eff.first=0.9912225008010864; eff.second=0.0009381904855814878;
                   }else                {eff.first=0.9955874085426331; eff.second=0.002450189400646501;
                   }
                }else if( eta < 2.4 ){
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
                if( eta < 0.9 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9752116203308105; eff.second=0.0030660638813280626;
                   }else if( pt < 30.0 ){eff.first=0.9848297238349915; eff.second=0.0016307213764927449;
                   }else if( pt < 40.0 ){eff.first=0.9861794114112854; eff.second=0.0006187187412138267;
                   }else if( pt < 50.0 ){eff.first=0.987443208694458; eff.second=0.000494159746725046;
                   }else if( pt < 60.0 ){eff.first=0.9834294319152832; eff.second=0.0011818999573518245;
                   }else                {eff.first=0.9863178730010986; eff.second=0.002073330940717176;
                   }
                }else if( eta < 1.2 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;                      
                   }else if( pt < 25.0 ){eff.first=0.9738101959228516; eff.second=0.004502934246978295;
                   }else if( pt < 30.0 ){eff.first=0.978645384311676; eff.second=0.0027064755458685794;
                   }else if( pt < 40.0 ){eff.first=0.9798933267593384; eff.second=0.001057081371390319;
                   }else if( pt < 50.0 ){eff.first=0.980233907699585; eff.second=0.000819615406448897;
                   }else if( pt < 60.0 ){eff.first=0.9773300886154175; eff.second=0.001955436343316424;
                   }else                {eff.first=0.9795225858688354; eff.second=0.0035622593553725837;
                   }
                }else if( eta < 2.1 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;
                   }else if( pt < 25.0 ){eff.first=0.9983288645744324; eff.second=0.002331323348626783;
                   }else if( pt < 30.0 ){eff.first=0.9905462265014648; eff.second=0.001402578599690647;
                   }else if( pt < 40.0 ){eff.first=0.9923668503761292; eff.second=0.0005653311393042486;
                   }else if( pt < 50.0 ){eff.first=0.9927627444267273; eff.second=0.0004155573807947332;
                   }else if( pt < 60.0 ){eff.first=0.9886322021484375; eff.second=0.0011254961157344963;
                   }else                {eff.first=0.9950451850891113; eff.second=0.002673833447209764;
                   }
                }else if( eta < 2.4 ){
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

                if( eta < 0.9 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;                            
                   }else if( pt < 25.0 ){ eff.first=1.0043761730194092; eff.second=0.003959090391076143;
                   }else if( pt < 30.0 ){ eff.first=0.9995378255844116; eff.second=0.0022512071035640673;
                   }else if( pt < 40.0 ){ eff.first=1.000901222229004; eff.second=0.0007979481788689052;
                   }else if( pt < 50.0 ){ eff.first=0.9986253976821899; eff.second=0.0004518361024064332;
                   }else if( pt < 60.0 ){ eff.first=1.0002487897872925; eff.second=0.000772847340102783;
                   }else                { eff.first=0.9986850023269653; eff.second=0.0008907575174433545;
                   }
                }else if( eta < 1.2 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;                            
                   }else if( pt < 25.0 ){ eff.first=1.004357933998108; eff.second=0.006125539530136138;
                   }else if( pt < 30.0 ){ eff.first=1.002331256866455; eff.second=0.004003683572512011;
                   }else if( pt < 40.0 ){ eff.first=1.004658579826355; eff.second=0.0014502638048416372;
                   }else if( pt < 50.0 ){ eff.first=1.0013608932495117; eff.second=0.0004888604573095644;
                   }else if( pt < 60.0 ){ eff.first=0.9986217021942139; eff.second=0.0012396364566794034;
                   }else                { eff.first=1.0054655075073242; eff.second=0.001589130019220112;
                   }
                }else if( eta < 2.1 ){
                         if( pt < 20.0 ){ eff.first=1.0000000000000000; eff.second=0.0000000000000000000;                            
                   }else if( pt < 25.0 ){ eff.first=0.9970762133598328; eff.second=0.003109125287470401;
                   }else if( pt < 30.0 ){ eff.first=1.0006532669067383; eff.second=0.002067755362435184;
                   }else if( pt < 40.0 ){ eff.first=1.0023553371429443; eff.second=0.0008445520691793605;
                   }else if( pt < 50.0 ){ eff.first=0.999933660030365; eff.second=0.0004309914887707696;
                   }else if( pt < 60.0 ){ eff.first=1.0002963542938232; eff.second=0.0007614160360063238;
                   }else                { eff.first=1.0004935264587402; eff.second=0.0009382223143922724;
                   }

                }else if( eta < 2.4 ){
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
 

	}
	break;
      }
      return eff;
    }


  std::pair<float,float> getTriggerEfficiencySF(float pt1, float eta1, float pt2, float eta2, int id){
     double etaCen = std::min(fabs(eta1), fabs(eta2));
     double etaFwd = std::max(fabs(eta1), fabs(eta2));
     id   = abs(id);
     std::pair<float,float> eff(1.0,0.00);

     if(id==121){
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
