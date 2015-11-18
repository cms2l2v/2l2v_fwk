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

		if( eta >= -2.5 && eta < -2 ){
		       if( pt < 30 ){ eff.first=0.997125; eff.second=0.0233165; }  
		       else if( pt < 40 ){ eff.first=1.01092; eff.second=0.0168794; }  
		       else if( pt < 50 ){ eff.first=1.02161; eff.second=0.0163831; }  
		       else { eff.first=1.01311; eff.second=0.0297784; }  
		} 
		else if( eta >= -2 && eta < -1.4 ){
		       if( pt < 30 ){ eff.first=0.976505; eff.second=0.0212189; }  
		       else if( pt < 40 ){ eff.first=0.992568; eff.second=0.0144407; }  
		       else if( pt < 50 ){ eff.first=1.00484; eff.second=0.0131707; }  
		       else { eff.first=0.996669; eff.second=0.0236046; }  
		} 
		else if( eta >= -1.4 && eta < -0.8 ){
		       if( pt < 30 ){ eff.first=0.988313; eff.second=0.0178967; }  
		       else if( pt < 40 ){ eff.first=0.990623; eff.second=0.0109575; }  
		       else if( pt < 50 ){ eff.first=0.990852; eff.second=0.0100034; }  
		       else { eff.first=0.990867; eff.second=0.0180543; }  
		} 
		else if( eta >= -0.8 && eta < 0 ){
		       if( pt < 30 ){ eff.first=0.996828; eff.second=0.0138775; }  
		       else if( pt < 40 ){ eff.first=0.993903; eff.second=0.00860885; }  
		       else if( pt < 50 ){ eff.first=0.996103; eff.second=0.00827937; }  
		       else { eff.first=0.996272; eff.second=0.0149864; }  
		} 
		else if( eta >= 0 && eta < 0.8 ){
		       if( pt < 30 ){ eff.first=0.995547; eff.second=0.0137682; }  
		       else if( pt < 40 ){ eff.first=0.994457; eff.second=0.00863161; }  
		       else if( pt < 50 ){ eff.first=0.995634; eff.second=0.00817332; }  
		       else { eff.first=0.994405; eff.second=0.0148898; }  
		} 
		else if( eta >= 0.8 && eta < 1.4 ){
		       if( pt < 30 ){ eff.first=0.989891; eff.second=0.0179614; }  
		       else if( pt < 40 ){ eff.first=0.993066; eff.second=0.0109637; }  
		       else if( pt < 50 ){ eff.first=0.992857; eff.second=0.0100039; }  
		       else { eff.first=0.985771; eff.second=0.0178857; }  
		} 
		else if( eta >= 1.4 && eta < 2 ){
		       if( pt < 30 ){ eff.first=0.97154; eff.second=0.0212848; }  
		       else if( pt < 40 ){ eff.first=0.982701; eff.second=0.0143374; }  
		       else if( pt < 50 ){ eff.first=1.00196; eff.second=0.0132433; }  
		       else { eff.first=1.00692; eff.second=0.023856; }  
		} 
		else {
		       if( pt < 30 ){ eff.first=0.989345; eff.second=0.0224178; }  
		       else if( pt < 40 ){ eff.first=1.00027; eff.second=0.0165292; }  
		       else if( pt < 50 ){ eff.first=1.01094; eff.second=0.0161848; }  
		       else { eff.first=1.00576; eff.second=0.0298075; }  
		} 

	    }
	  if(wp=="medium")
	    {

		if( eta >= -2.5 && eta < -2 ){
		       if( pt < 30 ){ eff.first=0.986325; eff.second=0.0235428; }  
		       else if( pt < 40 ){ eff.first=1.00215; eff.second=0.0169098; }  
		       else if( pt < 50 ){ eff.first=1.01334; eff.second=0.0163589; }  
		       else { eff.first=1.0077; eff.second=0.0297907; }  
		} 
		else if( eta >= -2 && eta < -1.4 ){
		       if( pt < 30 ){ eff.first=0.93943; eff.second=0.0211193; }  
		       else if( pt < 40 ){ eff.first=0.96489; eff.second=0.0143102; }  
		       else if( pt < 50 ){ eff.first=0.986757; eff.second=0.0130904; }  
		       else { eff.first=0.98219; eff.second=0.0234694; }  
		} 
		else if( eta >= -1.4 && eta < -0.8 ){
		       if( pt < 30 ){ eff.first=0.970436; eff.second=0.0181068; }  
		       else if( pt < 40 ){ eff.first=0.977071; eff.second=0.0110629; }  
		       else if( pt < 50 ){ eff.first=0.973989; eff.second=0.0100267; }  
		       else { eff.first=0.983606; eff.second=0.0182678; }  
		} 
		else if( eta >= -0.8 && eta < 0 ){
		       if( pt < 30 ){ eff.first=0.983966; eff.second=0.0140767; }  
		       else if( pt < 40 ){ eff.first=0.984833; eff.second=0.00871137; }  
		       else if( pt < 50 ){ eff.first=0.986977; eff.second=0.00835061; }  
		       else { eff.first=0.98556; eff.second=0.0150574; }  
		} 
		else if( eta >= 0 && eta < 0.8 ){
		       if( pt < 30 ){ eff.first=0.988178; eff.second=0.0140277; }  
		       else if( pt < 40 ){ eff.first=0.988077; eff.second=0.00875592; }  
		       else if( pt < 50 ){ eff.first=0.987754; eff.second=0.00825274; }  
		       else { eff.first=0.986672; eff.second=0.014995; }  
		} 
		else if( eta >= 0.8 && eta < 1.4 ){
		       if( pt < 30 ){ eff.first=0.980004; eff.second=0.018298; }  
		       else if( pt < 40 ){ eff.first=0.979701; eff.second=0.0110658; }  
		       else if( pt < 50 ){ eff.first=0.979428; eff.second=0.0100688; }  
		       else { eff.first=0.963679; eff.second=0.0177785; }  
		} 
		else if( eta >= 1.4 && eta < 2 ){
		       if( pt < 30 ){ eff.first=0.93626; eff.second=0.0212019; }  
		       else if( pt < 40 ){ eff.first=0.955496; eff.second=0.0142203; }  
		       else if( pt < 50 ){ eff.first=0.985745; eff.second=0.0131967; }  
		       else { eff.first=0.998763; eff.second=0.0238728; }  
		} 
		else {
		       if( pt < 30 ){ eff.first=0.973818; eff.second=0.0224985; }  
		       else if( pt < 40 ){ eff.first=0.991322; eff.second=0.0165683; }  
		       else if( pt < 50 ){ eff.first=1.00246; eff.second=0.0161585; }  
		       else { eff.first=0.9997; eff.second=0.0297967; }  
		} 

	    }
	  if(wp=="tight")
	    {

		if( eta >= -2.5 && eta < -2 ){
		       if( pt < 30 ){ eff.first=0.978767; eff.second=0.0248665; }  
		       else if( pt < 40 ){ eff.first=0.97874; eff.second=0.0172434; }  
		       else if( pt < 50 ){ eff.first=0.991385; eff.second=0.0166254; }  
		       else { eff.first=0.979362; eff.second=0.0298182; }  
		} 
		else if( eta >= -2 && eta < -1.4 ){
		       if( pt < 30 ){ eff.first=0.901559; eff.second=0.0213775; }  
		       else if( pt < 40 ){ eff.first=0.935026; eff.second=0.0144915; }  
		       else if( pt < 50 ){ eff.first=0.963243; eff.second=0.0132384; }  
		       else { eff.first=0.958387; eff.second=0.0235538; }  
		} 
		else if( eta >= -1.4 && eta < -0.8 ){
		       if( pt < 30 ){ eff.first=0.970829; eff.second=0.0184402; }  
		       else if( pt < 40 ){ eff.first=0.975505; eff.second=0.0111909; }  
		       else if( pt < 50 ){ eff.first=0.974722; eff.second=0.0101348; }  
		       else { eff.first=0.983324; eff.second=0.01841; }  
		} 
		else if( eta >= -0.8 && eta < 0 ){
		       if( pt < 30 ){ eff.first=0.978812; eff.second=0.0142061; }  
		       else if( pt < 40 ){ eff.first=0.982944; eff.second=0.00878671; }  
		       else if( pt < 50 ){ eff.first=0.985259; eff.second=0.00840972; }  
		       else { eff.first=0.980211; eff.second=0.0150757; }  
		} 
		else if( eta >= 0 && eta < 0.8 ){
		       if( pt < 30 ){ eff.first=0.989706; eff.second=0.0142538; }  
		       else if( pt < 40 ){ eff.first=0.985887; eff.second=0.00882929; }  
		       else if( pt < 50 ){ eff.first=0.986271; eff.second=0.00831021; }  
		       else { eff.first=0.987067; eff.second=0.0151122; }  
		} 
		else if( eta >= 0.8 && eta < 1.4 ){
		       if( pt < 30 ){ eff.first=0.976055; eff.second=0.0185183; }  
		       else if( pt < 40 ){ eff.first=0.976944; eff.second=0.0111706; }  
		       else if( pt < 50 ){ eff.first=0.978238; eff.second=0.0101587; }  
		       else { eff.first=0.959685; eff.second=0.0178474; }  
		} 
		else if( eta >= 1.4 && eta < 2 ){
		       if( pt < 30 ){ eff.first=0.918795; eff.second=0.0219825; }  
		       else if( pt < 40 ){ eff.first=0.93171; eff.second=0.0144731; }  
		       else if( pt < 50 ){ eff.first=0.961052; eff.second=0.0133106; }  
		       else { eff.first=0.978459; eff.second=0.0241004; }  
		} 
		else {
		       if( pt < 30 ){ eff.first=0.966533; eff.second=0.0236321; }  
		       else if( pt < 40 ){ eff.first=0.980386; eff.second=0.0171141; }  
		       else if( pt < 50 ){ eff.first=0.984021; eff.second=0.0164513; }  
		       else { eff.first=0.983836; eff.second=0.0302232; }  
		} 

	    }
	}
	break;
      case 13:
	{
	  if(wp=="loose")
	    {

		if( eta >= -2.4 && eta < -2.1 ){
		       if( pt < 25 ){ eff.first=0.965616; eff.second=0.0345805; }  
		       else if( pt < 30 ){ eff.first=1.0071; eff.second=0.0292463; }  
		       else if( pt < 35 ){ eff.first=0.997079; eff.second=0.0244509; }  
		       else if( pt < 40 ){ eff.first=1.00335; eff.second=0.0213907; }  
		       else if( pt < 50 ){ eff.first=0.995354; eff.second=0.0153824; }  
		       else if( pt < 60 ){ eff.first=1.00054; eff.second=0.0316722; }  
		       else if( pt < 90 ){ eff.first=1.00162; eff.second=0.0512166; }  
		       else if( pt < 140 ){ eff.first=0.996104; eff.second=0.147148; }  
		       else { eff.first=1.12137; eff.second=0.500173; }  
		} 
		else if( eta >= -2.1 && eta < -1.2 ){
		       if( pt < 25 ){ eff.first=0.998579; eff.second=0.0191078; }  
		       else if( pt < 30 ){ eff.first=1.00219; eff.second=0.015411; }  
		       else if( pt < 35 ){ eff.first=0.99643; eff.second=0.0128988; }  
		       else if( pt < 40 ){ eff.first=0.998464; eff.second=0.0106316; }  
		       else if( pt < 50 ){ eff.first=0.998774; eff.second=0.00710623; }  
		       else if( pt < 60 ){ eff.first=0.997989; eff.second=0.0149045; }  
		       else if( pt < 90 ){ eff.first=0.992338; eff.second=0.0232768; }  
		       else if( pt < 140 ){ eff.first=1.00031; eff.second=0.0614911; }  
		       else { eff.first=0.996679; eff.second=0.148159; }  
		} 
		else if( eta >= -1.2 && eta < -0.9 ){
		       if( pt < 25 ){ eff.first=1.0049; eff.second=0.032519; }  
		       else if( pt < 30 ){ eff.first=1.00113; eff.second=0.0251281; }  
		       else if( pt < 35 ){ eff.first=0.989505; eff.second=0.0197657; }  
		       else if( pt < 40 ){ eff.first=1.00032; eff.second=0.0161035; }  
		       else if( pt < 50 ){ eff.first=1.00258; eff.second=0.0111377; }  
		       else if( pt < 60 ){ eff.first=1.00004; eff.second=0.0235085; }  
		       else if( pt < 90 ){ eff.first=0.995198; eff.second=0.036159; }  
		       else if( pt < 140 ){ eff.first=0.990185; eff.second=0.0940945; }  
		       else { eff.first=0.973794; eff.second=0.266703; }  
		} 
		else if( eta >= -0.9 && eta < 0 ){
		       if( pt < 25 ){ eff.first=0.990469; eff.second=0.0176405; }  
		       else if( pt < 30 ){ eff.first=0.998884; eff.second=0.0127258; }  
		       else if( pt < 35 ){ eff.first=0.99787; eff.second=0.010034; }  
		       else if( pt < 40 ){ eff.first=0.997114; eff.second=0.00850312; }  
		       else if( pt < 50 ){ eff.first=0.99701; eff.second=0.00621451; }
			else if( pt < 60 ){ eff.first=0.993928; eff.second=0.013169; }  
		       else if( pt < 90 ){ eff.first=0.994234; eff.second=0.0201856; }  
		       else if( pt < 140 ){ eff.first=0.985573; eff.second=0.052201; }  
		       else { eff.first=0.983287; eff.second=0.120356; }  
		} 
		else if( eta >= 0 && eta < 0.9 ){
		       if( pt < 25 ){ eff.first=1.00449; eff.second=0.0179545; }  
		       else if( pt < 30 ){ eff.first=0.997415; eff.second=0.0128337; }  
		       else if( pt < 35 ){ eff.first=0.999627; eff.second=0.0100688; }  
		       else if( pt < 40 ){ eff.first=0.997356; eff.second=0.00849888; }  
		       else if( pt < 50 ){ eff.first=0.999919; eff.second=0.00623494; }  
		       else if( pt < 60 ){ eff.first=0.998088; eff.second=0.0133187; }  
		       else if( pt < 90 ){ eff.first=0.997032; eff.second=0.0200465; }  
		       else if( pt < 140 ){ eff.first=0.997212; eff.second=0.0515662; }  
		       else { eff.first=0.993146; eff.second=0.107855; }  
		} 
		else if( eta >= 0.9 && eta < 1.2 ){
		       if( pt < 25 ){ eff.first=0.988306; eff.second=0.0315506; }  
		       else if( pt < 30 ){ eff.first=1.01033; eff.second=0.0254511; }  
		       else if( pt < 35 ){ eff.first=1.00265; eff.second=0.0202679; }  
		       else if( pt < 40 ){ eff.first=1.00365; eff.second=0.0160612; }  
		       else if( pt < 50 ){ eff.first=0.995975; eff.second=0.011078; }  
		       else if( pt < 60 ){ eff.first=1.00131; eff.second=0.0237556; }  
		       else if( pt < 90 ){ eff.first=0.985025; eff.second=0.0353872; }  
		       else if( pt < 140 ){ eff.first=0.981008; eff.second=0.0948365; }  
		       else { eff.first=0.981689; eff.second=0.207412; }  
		} 
		else if( eta >= 1.2 && eta < 2.1 ){
		       if( pt < 25 ){ eff.first=0.991561; eff.second=0.0187536; }  
		       else if( pt < 30 ){ eff.first=0.995135; eff.second=0.0153864; }  
		       else if( pt < 35 ){ eff.first=0.996955; eff.second=0.0128287; }  
		       else if( pt < 40 ){ eff.first=0.999252; eff.second=0.0106582; }  
		       else if( pt < 50 ){ eff.first=0.998156; eff.second=0.00706966; }  
		       else if( pt < 60 ){ eff.first=0.998532; eff.second=0.0149931; }  
		       else if( pt < 90 ){ eff.first=0.999152; eff.second=0.0237321; }  
		       else if( pt < 140 ){ eff.first=0.992795; eff.second=0.059789; }  
		       else { eff.first=0.993138; eff.second=0.145493; }  
		} 
		else {
		       if( pt < 25 ){ eff.first=0.972257; eff.second=0.0352844; }  
		       else if( pt < 30 ){ eff.first=0.993037; eff.second=0.0291596; }  
		       else if( pt < 35 ){ eff.first=1.00289; eff.second=0.0245242; }  
		       else if( pt < 40 ){ eff.first=0.997007; eff.second=0.0212478; }  
		       else if( pt < 50 ){ eff.first=0.998648; eff.second=0.0154273; }  
		       else if( pt < 60 ){ eff.first=0.99261; eff.second=0.0320107; }  
		       else if( pt < 90 ){ eff.first=0.993148; eff.second=0.0532669; }  
		       else if( pt < 140 ){ eff.first=1.00947; eff.second=0.143595; }  
		       else { eff.first=1.00433; eff.second=0.48198; }  
		} 

	    }
          if(wp=="medium")
           {

		if( eta >= -2.4 && eta < -2.1 ){
		       if( pt < 25 ){ eff.first=1.02571; eff.second=0.0378829; }  
		       else if( pt < 30 ){ eff.first=1.03564; eff.second=0.0308776; }  
		       else if( pt < 35 ){ eff.first=1.03845; eff.second=0.0261473; }  
		       else if( pt < 40 ){ eff.first=1.05754; eff.second=0.0231476; }  
		       else if( pt < 50 ){ eff.first=1.04513; eff.second=0.0165309; }  
		       else if( pt < 60 ){ eff.first=1.05099; eff.second=0.0340864; }  
		       else if( pt < 90 ){ eff.first=1.04946; eff.second=0.0551666; }  
		       else if( pt < 140 ){ eff.first=1.06673; eff.second=0.162059; }  
		       else { eff.first=1.15678; eff.second=0.540309; }  
		} 
		else if( eta >= -2.1 && eta < -1.2 ){
		       if( pt < 25 ){ eff.first=0.996322; eff.second=0.0193362; }  
		       else if( pt < 30 ){ eff.first=1.00114; eff.second=0.0155933; }  
		       else if( pt < 35 ){ eff.first=0.990151; eff.second=0.0129817; }  
		       else if( pt < 40 ){ eff.first=0.994275; eff.second=0.0107205; }  
		       else if( pt < 50 ){ eff.first=0.992173; eff.second=0.00714342; }  
		       else if( pt < 60 ){ eff.first=0.991019; eff.second=0.0149738; }  
		       else if( pt < 90 ){ eff.first=0.983683; eff.second=0.0233376; }  
		       else if( pt < 140 ){ eff.first=1.00059; eff.second=0.0620693; }  
		       else { eff.first=0.967963; eff.second=0.145003; }  
		} 
		else if( eta >= -1.2 && eta < -0.9 ){
		       if( pt < 25 ){ eff.first=1.00588; eff.second=0.0337015; }  
		       else if( pt < 30 ){ eff.first=0.975118; eff.second=0.025322; }  
		       else if( pt < 35 ){ eff.first=0.982748; eff.second=0.0203172; }  
		       else if( pt < 40 ){ eff.first=0.997692; eff.second=0.0165881; }  
		       else if( pt < 50 ){ eff.first=0.995606; eff.second=0.0114128; }  
		       else if( pt < 60 ){ eff.first=0.996507; eff.second=0.0241764; }  
		       else if( pt < 90 ){ eff.first=0.998382; eff.second=0.0373693; }  
		       else if( pt < 140 ){ eff.first=0.959569; eff.second=0.0936826; }  
		       else { eff.first=0.990039; eff.second=0.282097; }  
		} 
		else if( eta >= -0.9 && eta < 0 ){
		       if( pt < 25 ){ eff.first=0.991846; eff.second=0.018262; }  
		       else if( pt < 30 ){ eff.first=1.003; eff.second=0.0131838; }  
		       else if( pt < 35 ){ eff.first=1.00718; eff.second=0.0104463; }  
		       else if( pt < 40 ){ eff.first=1.00049; eff.second=0.0087867; }  
		       else if( pt < 50 ){ eff.first=1.00169; eff.second=0.00642563; }  
		       else if( pt < 60 ){ eff.first=0.990352; eff.second=0.0134724; }  
		       else if( pt < 90 ){ eff.first=0.993943; eff.second=0.0207645; }  
		       else if( pt < 140 ){ eff.first=0.994379; eff.second=0.0541797; }  
		       else { eff.first=1.01263; eff.second=0.130008; }  
		} 
		else if( eta >= 0 && eta < 0.9 ){
		       if( pt < 25 ){ eff.first=1.01484; eff.second=0.0187417; }  
		       else if( pt < 30 ){ eff.first=1.00679; eff.second=0.0133779; }  
		       else if( pt < 35 ){ eff.first=1.01778; eff.second=0.010568; }  
		       else if( pt < 40 ){ eff.first=1.01085; eff.second=0.00886548; }  
		       else if( pt < 50 ){ eff.first=1.01476; eff.second=0.0065113; }  
		       else if( pt < 60 ){ eff.first=1.00806; eff.second=0.0138288; }  
		       else if( pt < 90 ){ eff.first=1.00594; eff.second=0.0207906; }  
		       else if( pt < 140 ){ eff.first=1.00605; eff.second=0.0534435; }  
		       else { eff.first=0.950105; eff.second=0.105648; }  
		} 
		else if( eta >= 0.9 && eta < 1.2 ){
		       if( pt < 25 ){ eff.first=0.995476; eff.second=0.0329939; }  
		       else if( pt < 30 ){ eff.first=1.02444; eff.second=0.0267259; }  
		       else if( pt < 35 ){ eff.first=1.01852; eff.second=0.0213166; }  
		       else if( pt < 40 ){ eff.first=1.01207; eff.second=0.01672; }  
		       else if( pt < 50 ){ eff.first=1.00434; eff.second=0.0115271; }  
		       else if( pt < 60 ){ eff.first=1.0085; eff.second=0.0246806; }  
		       else if( pt < 90 ){ eff.first=0.993844; eff.second=0.036828; }  
		       else if( pt < 140 ){ eff.first=1.00346; eff.second=0.100301; }  
		       else {eff.first=0.989995; eff.second=0.216885; }  
		} 
		else if( eta >= 1.2 && eta < 2.1 ){
		       if( pt < 25 ){ eff.first=0.994996; eff.second=0.0190843; }  
		       else if( pt < 30 ){ eff.first=0.993814; eff.second=0.0155744; }  
		       else if( pt < 35 ){ eff.first=1.0003; eff.second=0.0130556; }  
		       else if( pt < 40 ){ eff.first=0.998499; eff.second=0.0107898; }  
		       else if( pt < 50 ){ eff.first=0.999478; eff.second=0.00716629; }  
		       else if( pt < 60 ){ eff.first=0.997264; eff.second=0.015153; }  
		       else if( pt < 90 ){ eff.first=1.00114; eff.second=0.0240556; }  
		       else if( pt < 140 ){ eff.first=0.986345; eff.second=0.060025; }  
		       else { eff.first=0.999144; eff.second=0.14824; }  
		} 
		else {
		       if( pt < 25 ){ eff.first=0.989057; eff.second=0.0371012; }  
		       else if( pt < 30 ){ eff.first=1.0124; eff.second=0.0307447; }  
		       else if( pt < 35 ){ eff.first=1.02534; eff.second=0.0258689; }  
		       else if( pt < 40 ){ eff.first=1.01988; eff.second=0.0223336; }  
		       else if( pt < 50 ){ eff.first=1.02514; eff.second=0.0162777; }  
		       else if( pt < 60 ){ eff.first=1.01158; eff.second=0.0335998; }  
		       else if( pt < 90 ){ eff.first=1.00015; eff.second=0.0551329; }  
		       else if( pt < 140 ){ eff.first=0.997006; eff.second=0.145836; }  
		       else { eff.first=1.09872; eff.second=0.554848; }  
		} 

           }
	  if(wp=="tight")
           {

		if( eta >= -2.4 && eta < -2.1 ){
		       if( pt < 25 ){ eff.first=0.957794; eff.second=0.0344595; }  
		       else if( pt < 30 ){ eff.first=1.00071; eff.second=0.0292252; }  
		       else if( pt < 35 ){ eff.first=0.986722; eff.second=0.024328; }  
		       else if( pt < 40 ){ eff.first=1.001; eff.second=0.021475; }  
		       else if( pt < 50 ){ eff.first=0.986578; eff.second=0.0153291; }  
		       else if( pt < 60 ){ eff.first=0.990299; eff.second=0.0315689; }  
		       else if( pt < 90 ){ eff.first=0.991747; eff.second=0.0510851; }  
		       else if( pt < 140 ){ eff.first=1.00194; eff.second=0.149376; }  
		       else { eff.first=1.12206; eff.second=0.508404; }  
		} 
		else if( eta >= -2.1 && eta < -1.2 ){
		       if( pt < 25 ){ eff.first=0.994641; eff.second=0.0190657; }  
		       else if( pt < 30 ){ eff.first=0.998472; eff.second=0.0153811; }  
		       else if( pt < 35 ){ eff.first=0.992344; eff.second=0.0128737; }  
		       else if( pt < 40 ){ eff.first=0.994021; eff.second=0.0106073; }  
		       else if( pt < 50 ){ eff.first=0.994611; eff.second=0.00709287; }  
		       else if( pt < 60 ){ eff.first=0.994233; eff.second=0.0148889; }  
		       else if( pt < 90 ){ eff.first=0.985544; eff.second=0.0231665; }  
		       else if( pt < 140 ){ eff.first=0.98566; eff.second=0.0607459; }  
		       else { eff.first=0.991678; eff.second=0.14838; }  
		} 
		else if( eta >= -1.2 && eta < -0.9 ){
		       if( pt < 25 ){ eff.first=0.996489; eff.second=0.0323875; }  
		       else if( pt < 30 ){ eff.first=0.992843; eff.second=0.0250381; }  
		       else if( pt < 35 ){ eff.first=0.983373; eff.second=0.0197413; }  
		       else if( pt < 40 ){ eff.first=0.994228; eff.second=0.0160745; }  
		       else if( pt < 50 ){ eff.first=0.996489; eff.second=0.0111176; }  
		       else if( pt < 60 ){ eff.first=0.990453; eff.second=0.0233958; }  
		       else if( pt < 90 ){ eff.first=0.982464; eff.second=0.0358538; }  
		       else if( pt < 140 ){ eff.first=0.985887; eff.second=0.0943882; }  
		       else { eff.first=0.946561; eff.second=0.262006; }  
		} 
		else if( eta >= -0.9 && eta < 0 ){
		       if( pt < 25 ){ eff.first=0.990533; eff.second=0.0177305; }  
		       else if( pt < 30 ){ eff.first=1.00012; eff.second=0.0128033; }  
		       else if( pt < 35 ){ eff.first=1.00091; eff.second=0.010119; }  
		       else if( pt < 40 ){ eff.first=1.00001; eff.second=0.00857493; }  
		       else if( pt < 50 ){ eff.first=1.0019; eff.second=0.00628073; }  
		       else if( pt < 60 ){ eff.first=0.994783; eff.second=0.0132564; }  
		       else if( pt < 90 ){ eff.first=0.989928; eff.second=0.0202109; }  
		       else if( pt < 140 ){ eff.first=0.972618; eff.second=0.0516637; }  
		       else { eff.first=0.959046; eff.second=0.117673; }  
		} 
		else if( eta >= 0 && eta < 0.9 ){
		       if( pt < 25 ){ eff.first=1.00463; eff.second=0.018064; }  
		       else if( pt < 30 ){ eff.first=0.997174; eff.second=0.0128978; }  
		       else if( pt < 35 ){ eff.first=0.998746; eff.second=0.0101131; }  
		       else if( pt < 40 ){ eff.first=0.996968; eff.second=0.00854337; }  
		       else if( pt < 50 ){ eff.first=1.00225; eff.second=0.00628758; }  
		       else if( pt < 60 ){ eff.first=0.995751; eff.second=0.0133663; }  
		       else if( pt < 90 ){ eff.first=0.992303; eff.second=0.0200658; }  
		       else if( pt < 140 ){ eff.first=0.98722; eff.second=0.0512357; }  
		       else { eff.first=0.989818; eff.second=0.108148; }  
		} 
		else if( eta >= 0.9 && eta < 1.2 ){
		       if( pt < 25 ){ eff.first=0.973148; eff.second=0.0312361; }  
		       else if( pt < 30 ){ eff.first=0.992351; eff.second=0.0251352; }  
		       else if( pt < 35 ){ eff.first=0.986834; eff.second=0.0200666; }  
		       else if( pt < 40 ){ eff.first=0.986431; eff.second=0.0158823; }  
		       else if( pt < 50 ){ eff.first=0.979025; eff.second=0.0109596; }  
		       else if( pt < 60 ){ eff.first=0.982174; eff.second=0.0234495; }  
		       else if( pt < 90 ){ eff.first=0.969787; eff.second=0.0350701; }  
		       else if( pt < 140 ){ eff.first=0.970168; eff.second=0.0945557; }  
		       else { eff.first=0.926858; eff.second=0.19651; }  
		} 
		else if( eta >= 1.2 && eta < 2.1 ){
		       if( pt < 25 ){ eff.first=0.987597; eff.second=0.0187088; }  
		       else if( pt < 30 ){ eff.first=0.992433; eff.second=0.0153711; }  
		       else if( pt < 35 ){ eff.first=0.995579; eff.second=0.0128378; }  
		       else if( pt < 40 ){ eff.first=0.996672; eff.second=0.0106535; }  
		       else if( pt < 50 ){ eff.first=0.995054; eff.second=0.00706222; }  
		       else if( pt < 60 ){ eff.first=0.996634; eff.second=0.0150105; }  
		       else if( pt < 90 ){ eff.first=0.996727; eff.second=0.0237315; }  
		       else if( pt < 140 ){ eff.first=0.994457; eff.second=0.060213; }  
		       else { eff.first=1.0061; eff.second=0.148525; }  
		} 
		else {
		       if( pt < 25 ){ eff.first=0.961479; eff.second=0.0350677; }  
		       else if( pt < 30 ){ eff.first=0.988415; eff.second=0.0291594; }  
		       else if( pt < 35 ){ eff.first=0.997842; eff.second=0.0245271; }  
		       else if( pt < 40 ){ eff.first=0.991395; eff.second=0.0212372; }  
		       else if( pt < 50 ){ eff.first=0.993476; eff.second=0.0154276; }  
		       else if( pt < 60 ){ eff.first=0.988009; eff.second=0.0320619; }  
		       else if( pt < 90 ){ eff.first=0.988104; eff.second=0.053417; }  
		       else if( pt < 140 ){ eff.first=1.01729; eff.second=0.146525; }  
		       else { eff.first=0.932599; eff.second=0.451738; }  
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
