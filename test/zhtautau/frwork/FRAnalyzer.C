#include "FRAnalyzer.h"

FRAnalyzer::FRAnalyzer(){
}
FRAnalyzer::~FRAnalyzer(){
}

FRAnalyzer::FRAnalyzer(string File1,
		string File2,
		int Index,
		int Ntoy,
		bool Data,
		string Dir){

	m_File1     = File1;
	m_File2     = File2;
	m_Index     = Index;
	m_Ntoy      = Ntoy;
	m_Data      = Data;
	m_Dir       = Dir;

	TH1::SetDefaultSumw2(kTRUE);
	TH2::SetDefaultSumw2(kTRUE);
	gStyle->SetOptStat(111);
	gStyle->SetOptFit(1011);
	CreateDir();
}

void FRAnalyzer::CreateDir(){
	struct stat st;
	if(stat(m_Dir.c_str(),&st) == 0){
		std::cout<<"Out Directory "<< m_Dir <<" already present!"<<std::endl;
	}else{
		std::cout<<"Creating Directory "<< m_Dir <<" ... "<<std::endl;
		int outD = system( ("mkdir "+m_Dir).c_str() );
		if(outD!=0)
			std::cout<<"Directory "<< m_Dir <<" could not be created!"<< std::endl;
	}
}

void FRAnalyzer::ReadFiles(){
	file1 = new TFile( (m_File1).c_str() );
	file2 = new TFile( (m_File2).c_str() );
}

void FRAnalyzer::Initialize(){
	debug=false;
	sample = "";
	finalStates.clear();
	channels.clear();
	shapeName.clear();
	tauIso.clear();
	muoIso.clear();
	eleIso.clear();
	hlabel.clear();

	if(m_Data){ 
	  sample = "data";
	}else{
	  sample = "Z#rightarrow ll";
	}
	
	cidx<<m_Index;
	cidx>>Idx;
	hlabel.push_back("");
	hlabel.push_back("MC_");
	cutIndexLabel.push_back("L_{T}>0GeV");
	cutIndexLabel.push_back("L_{T}>20GeV");
	cutIndexLabel.push_back("L_{T}>40GeV");
	cutIndexLabel.push_back("L_{T}>60GeV");
	cutIndexLabel.push_back("L_{T}>80GeV");
	cutIndexLabel.push_back("L_{T}>100GeV");
	cutIndexLabel.push_back("L_{T}>120GeV");
	cutIndexLabel.push_back("L_{T}>140GeV");
	cutIndexLabel.push_back("L_{T}>160GeV");
	cutIndexLabel.push_back("L_{T}>180GeV");
	cutIndexLabel.push_back("L_{T}>200GeV");
	shapeName.push_back("_Hsvfit_shapes");
	shapeName.push_back("_Asvfit_shapes");
	finalStates.push_back("ee_OSelmu");
	finalStates.push_back("ee_OSelha");
	finalStates.push_back("ee_OSmuha");
	finalStates.push_back("ee_OShaha");
	finalStates.push_back("mumu_OSelmu");
	finalStates.push_back("mumu_OSelha");
	finalStates.push_back("mumu_OSmuha");
	finalStates.push_back("mumu_OShaha");
	finalStates.push_back("ee_SSelmu");
	finalStates.push_back("ee_SSelha");
	finalStates.push_back("ee_SSmuha");
	finalStates.push_back("ee_SShaha");
	finalStates.push_back("mumu_SSelmu");
	finalStates.push_back("mumu_SSelha");
	finalStates.push_back("mumu_SSmuha");
	finalStates.push_back("mumu_SShaha");
	forName.push_back("ele01J");
	forName.push_back("ele02J");
	forName.push_back("ele03J");
	forName.push_back("muo01J");
	forName.push_back("muo02J");
	forName.push_back("muo03J");
	forName.push_back("tauLJ");
	forName.push_back("tauMJ");
	forName.push_back("ele01J_B");
	forName.push_back("ele02J_B");
	forName.push_back("ele03J_B");
	forName.push_back("muo01J_B");
	forName.push_back("muo02J_B");
	forName.push_back("muo03J_B");
	forName.push_back("tauLJ_B");
	forName.push_back("tauMJ_B");
	forName.push_back("ele01J_E");
	forName.push_back("ele02J_E");
	forName.push_back("ele03J_E");
	forName.push_back("muo01J_E");
	forName.push_back("muo02J_E");
	forName.push_back("muo03J_E");
	forName.push_back("tauLJ_E");
	forName.push_back("tauMJ_E");
	channels.push_back("eee#mu");
	channels.push_back("eee#tau_{h}");
	channels.push_back("ee#mu#tau_{h}");
	channels.push_back("ee#tau_{h}#tau_{h}");
	channels.push_back("#mu#mue#mu");
	channels.push_back("#mu#mue#tau_{h}");
	channels.push_back("#mu#mu#mu#tau_{h}");
	channels.push_back("#mu#mu#tau_{h}#tau_{h}");
	dir_data         = "data";               
	dir_wz           = "WZ";               
	dir_ww           = "WW";               
	dir_vg           = "V#gamma";          
	dir_ttZ          = "ttZ";       
	dir_top          = "tT,tTV,t,T";       
	dir_qcd          = "QCD";
	dir_qcdmu        = "QCD Mu";
	dir_wjets        = "W+jets";           
	dir_dy           = "Z#rightarrow ll";  
	dir_zz2l2tau     = "ZZ#rightarrow ll#tau#tau";
	dir_zz4l         = "ZZ#rightarrow 4l";
	dir_zz2l2nu      = "ZZ#rightarrow ll#nu#nu";
	dir_wwz          = "WWZ";
	dir_wzz          = "WZZ";
	dir_zzz          = "ZZZ";
	dir_www          = "WWW";
	dir_zh           = "Zh(125GeV)#rightarrow ll#tau#tau";

	double bins[]={5.0,32.4,70.2,110.2,189.,306.,550.8,1785.}; //OLD BINNING
	//double bins[]={5.0,32.4,70.2,110.2,170.,189.,306.,550.8,1785.}; //OLD BINNING
        nbins = 7;

	hBins = new TH2F( "bins2","",nbins,bins,nbins,bins);
	for(int i=1;i<=hBins->GetXaxis()->GetNbins();i++){
		for(int j=1;j<=hBins->GetXaxis()->GetNbins();j++){
			hBins->SetBinContent(i,j, hBins->GetBin(i,j) );
		}
	}

	h2ref      = (TH2F*)file1->Get("data/ee_OSelmu_Hsvfit_shapes");
	nbx        = h2ref->GetXaxis()->GetNbins();
	maxx       = h2ref->GetXaxis()->GetXmax();
	minx       = h2ref->GetXaxis()->GetXmin();

	DefineLikelihoodFunction();
        fitFun   = new TF1("fitFun","[0]*ROOT::Math::lognormal_pdf(x,[1],[2])");
}

//|******************|
//| Maths            |
//|******************|

TH1F* FRAnalyzer::Subtract(TH1F* h1, TH1F* h2){
	TH1F* h = (TH1F*)h1->Clone();
	h->Add(h2,-1);
	return h;
}
//|******************|
//|******************|

TPaveText* FRAnalyzer::GetPaveTextCMS(){

	TPaveText *T = new TPaveText(0.1,0.9,0.29,0.96,"brNDC");
	T->SetBorderSize(0);
	T->SetFillColor(0);
	T->SetFillStyle(0);
	T->SetLineColor(0);
	T->SetTextFont(62);
	T->SetTextAlign(12);
	T->SetTextSize(0.06);
	T->SetBorderSize(0);
	T->AddText("CMS");
	return T;
}

TPaveText* FRAnalyzer::GetPaveTextSim(bool simulation){
	TPaveText *T = new TPaveText(0.2332215,0.8916084,0.4228188,0.951049,"brNDC");
	T->SetBorderSize(0);
	T->SetFillColor(0);
	T->SetFillStyle(0);
	T->SetLineColor(0);
	T->SetTextAlign(12);
	T->SetTextFont(52);
	T->SetTextSize(0.04);
	if(simulation) T->AddText("Simulation");
	else T->AddText("Preliminary");
	return T;
}

TPaveText* FRAnalyzer::GetPaveTextLumi(){
	TPaveText *T = new TPaveText(0.6392617,0.8968531,0.8389262,0.958042,"brNDC");
	T->SetBorderSize(0);
	T->SetFillColor(0);
	T->SetFillStyle(0);
	T->SetLineColor(0);
	T->SetTextAlign(12);
	T->SetTextSize(0.04);
	T->SetTextFont(42);
	T->AddText("19.8 fb^{-1} (8 TeV)");
	return T;
}

TLatex* FRAnalyzer::ClosureTestLatexText(double x, double y, string channel){
	string name;
	if(channel.compare("em")==0) name = "eee#mu+#mu#mue#mu";
	if(channel.compare("et")==0) name = "eee#tau_{h}+#mu#mue#tau_{h}";
	if(channel.compare("mt")==0) name = "ee#mu#tau_{h}+#mu#mu#mu#tau_{h}";
	if(channel.compare("tt")==0) name = "ee#tau_{h}#tau_{h}+#mu#mu#tau_{h}#tau_{h}";
	TLatex *t = new TLatex(x,y,name.c_str());        
	t->SetNDC();
	t->SetTextColor(kBlue);
	t->SetTextFont(62);
	t->SetTextAlign(12);
	t->SetTextSize(0.045);
	return t;
}

vector<TLatex*> FRAnalyzer::FRlatexText(int i, double x, double y){
	string name;
	string cut;
	vector<TLatex*> vt;
	if(i==0 || i==1 || i==2)    name = "Electron Fake Rate";
	if(i==3 || i==4 || i==5)    name = "Muon Fake Rate";
	if(i==6 || i==7)            name = "Tau Fake Rate";
	if(i==8 || i==9 || i==10)   name = "Electron Prompt Ratio";
	if(i==11 || i==12 || i==13) name = "Muon Prompt Ratio";
	if(i==14 || i==15)          name = "Tau Prompt Ratio";
	if(i==0 || i==3)            cut  = "I_{rel}<0.1 (Barrel+Endcap)";
	if(i==1 || i==4)            cut  = "I_{rel}<0.2 (Barrel+Endcap)";
	if(i==2 || i==5)            cut  = "I_{rel}<0.3 (Barrel+Endcap)";
	if(i==6)                    cut  = "Loose-#Delta#beta 3hits (Barrel+Endcap)";
	if(i==7)                    cut  = "Medium-#Delta#beta 3hits (Barrel+Endcap)";
	if(i==8 || i==11)           cut  = "I_{rel}<0.1 (Barrel)";
	if(i==9 || i==12)           cut  = "I_{rel}<0.2 (Barrel)";
	if(i==10 || i==13)          cut  = "I_{rel}<0.3 (Barrel)";
	if(i==14)                   cut  = "Loose-#Delta#beta 3hits (Barrel)";
	if(i==15)                   cut  = "Medium-#Delta#beta 3hits (Barrel)";
	if(i==16 || i==19)          cut  = "I_{rel}<0.1 (Endcap)";
	if(i==17 || i==20)          cut  = "I_{rel}<0.2 (Endcap)";
	if(i==18 || i==21)          cut  = "I_{rel}<0.3 (Endcap)";
	if(i==22)                   cut  = "Loose-#Delta#beta 3hits (Endcap)";
	if(i==23)                   cut  = "Medium-#Delta#beta 3hits (Endcap)";
	TLatex *t1 = new TLatex(x,y,name.c_str());        
	TLatex *t2 = new TLatex(x,y+0.05,cut.c_str());        
	t1->SetNDC(); t2->SetNDC();
	t1->SetTextColor(kBlue); t2->SetTextColor(kBlue);
	t1->SetTextFont(62); t2->SetTextFont(42);
	t1->SetTextAlign(12); t2->SetTextAlign(12);
	t1->SetTextSize(0.045); t2->SetTextSize(0.045);
	vt.push_back(t1);
	vt.push_back(t2);
	return vt;
}


//|******************|
vector<TH1F*> FRAnalyzer::GetFRHistograms(bool numerator, string dir){
	TH1F* h1, *h2, *h3, *h4;
	TH1F* h5, *h6, *h7, *h8;
	TH1F* h9, *h10, *h11, *h12;
	TH1F* h13, *h14, *h15, *h16;
	if(numerator){
		h1 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleNumerator01Jet").c_str()) ;
		h2 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleNumerator02Jet").c_str()) ;
		h3 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleNumerator03Jet").c_str()) ;
		h4 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuNumerator01Jet").c_str()) ;
		h5 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuNumerator02Jet").c_str()) ;
		h6 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuNumerator03Jet").c_str()) ;
		h7 = (TH1F*)file1->Get((dir+"/FR_closestJetToTauNumeratorLJet").c_str()) ;
		h8 = (TH1F*)file1->Get((dir+"/FR_closestJetToTauNumeratorMJet").c_str()) ;

		h9 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleNumerator01Lep").c_str()) ;
		h10 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleNumerator02Lep").c_str()) ;
		h11 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleNumerator03Lep").c_str()) ;
		h12 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuNumerator01Lep").c_str()) ;
		h13 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuNumerator02Lep").c_str()) ;
		h14 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuNumerator03Lep").c_str()) ;
		h15 = (TH1F*)file1->Get((dir+"/FR_closestJetToTauNumeratorLLep").c_str()) ;
		h16 = (TH1F*)file1->Get((dir+"/FR_closestJetToTauNumeratorMLep").c_str()) ;
	}
	else{

		h1 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleDenominatorJet").c_str()) ;
		h2 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleDenominatorJet").c_str()) ;
		h3 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleDenominatorJet").c_str()) ;
		h4 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuDenominatorJet").c_str()) ;
		h5 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuDenominatorJet").c_str()) ;
		h6 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuDenominatorJet").c_str()) ;
		h7 = (TH1F*)file1->Get((dir+"/FR_closestJetToTauDenominatorJet").c_str()) ;
		h8 = (TH1F*)file1->Get((dir+"/FR_closestJetToTauDenominatorJet").c_str()) ;

		h9 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleDenominatorLep").c_str()) ;
		h10 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleDenominatorLep").c_str()) ;
		h11 = (TH1F*)file1->Get((dir+"/FR_closestJetToEleDenominatorLep").c_str()) ;
		h12 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuDenominatorLep").c_str()) ;
		h13 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuDenominatorLep").c_str()) ;
		h14 = (TH1F*)file1->Get((dir+"/FR_closestJetToMuDenominatorLep").c_str()) ;
		h15 = (TH1F*)file1->Get((dir+"/FR_closestJetToTauDenominatorLep").c_str()) ;
		h16 = (TH1F*)file1->Get((dir+"/FR_closestJetToTauDenominatorLep").c_str()) ;
	}

	vector<TH1F*> vec;
	vec.push_back(h1); vec.push_back(h2); vec.push_back(h3); vec.push_back(h4);
	vec.push_back(h5); vec.push_back(h6); vec.push_back(h7); vec.push_back(h8);
	vec.push_back(h9); vec.push_back(h10); vec.push_back(h11); vec.push_back(h12);
	vec.push_back(h13); vec.push_back(h14); vec.push_back(h15); vec.push_back(h16);

	return vec;
}
		
TF1* FRAnalyzer::GetLikelihoodFunction(int count){
  
  if(count<=1) return f01;         // = new TF1("f1","TMath::Poisson([0],x)",0,10);
  else if(count<=5)    return f02; // = new TF1("f1","TMath::Poisson([0],x)",0,20);
  else if(count<=10)   return f03; // = new TF1("f1","TMath::Poisson([0],x)",0,30);
  else if(count<=30){
    f04->SetRange(0,count+30);
    return f04; // = new TF1("f1","TMath::Poisson([0],x)",0,count+30);
  }
  else if(count<=50){
    f05->SetRange(0,count+40);
    return f05; // = new TF1("f1","TMath::Poisson([0],x)",0,count+40);
  }
  else if(count<=70){
    f06->SetRange(count-30,count+50);
    return f06; // = new TF1("f1","TMath::Poisson([0],x)",count-30,count+50);
  }
  else if(count<=100){
    f07->SetRange(count-40,count+60);
    return f07; // = new TF1("f1","TMath::Poisson([0],x)",count-40,count+60);
  }
  else if(count<=200){
    f08->SetRange(count-50,count+70);
    return f08; // = new TF1("f1","TMath::Poisson([0],x)",count-50,count+70);
  }
  else if(count<=600){
    f09->SetRange(count-100,count+100);
    return f09; // = new TF1("f1","TMath::Poisson([0],x)",count-100,count+100);
  }
  else if(count<=2000){
    f10->SetRange(count-200,count+200);
    return f10; // = new TF1("f1","TMath::Poisson([0],x)",count-200,count+200);
  }
  else{
    f11->SetRange(count-300,count+300);
    return f11; // = new TF1("f1","TMath::Poisson([0],x)",count-300,count+300);
  }
  //f1->SetParameter(0,count);
  //return f1;
}

TF1* FRAnalyzer::GetFitFunction(){
	return fitFun;
}

void FRAnalyzer::DefineLikelihoodFunction(){

  f01 = new TF1("f1","TMath::Poisson([0],x)",0,10);
  f02 = new TF1("f1","TMath::Poisson([0],x)",0,20);
  f03 = new TF1("f1","TMath::Poisson([0],x)",0,30);
  f04 = new TF1("f1","TMath::Poisson([0],x)"); //,0,count+30);
  f05 = new TF1("f1","TMath::Poisson([0],x)"); //,0,count+40);
  f06 = new TF1("f1","TMath::Poisson([0],x)"); //,count-30,count+50);
  f07 = new TF1("f1","TMath::Poisson([0],x)"); //,count-40,count+60);
  f08 = new TF1("f1","TMath::Poisson([0],x)"); //,count-50,count+70);
  f09 = new TF1("f1","TMath::Poisson([0],x)"); //,count-100,count+100);
  f10 = new TF1("f1","TMath::Poisson([0],x)"); //,count-200,count+200);
  f11 = new TF1("f1","TMath::Poisson([0],x)"); //,count-300,count+300);

}

vector<TH1F*> FRAnalyzer::GetRatiosData(){

	vector<TH1F*> num = GetFRHistograms(true,sample);
	vector<TH1F*> den = GetFRHistograms(false,sample);
	vector<TH1F*> ratios;
	for(int i=0;i<(int)num.size();i++){
		TH1F *tmpNum = (TH1F*)num.at(i)->Clone(); 
		TH1F *tmpDen = (TH1F*)den.at(i)->Clone();
		tmpNum->Divide(tmpDen); 
		ratios.push_back(tmpNum);
	}
	return ratios;
}

vector<TH1F*> FRAnalyzer::GetRatios(bool rebin){

	vector<TH1F*> num = GetFRHistograms(true,sample);
	vector<TH1F*> den = GetFRHistograms(false,sample);
	vector<TH1F*> ratios;

	if(m_Data){ 
		vector<TH1F*> num2 = GetFRHistograms(true,dir_wz);
		vector<TH1F*> den2 = GetFRHistograms(false,dir_wz);
		vector<TH1F*> num3 = GetFRHistograms(true,dir_zz4l);
		vector<TH1F*> den3 = GetFRHistograms(false,dir_zz4l);
		vector<TH1F*> num4 = GetFRHistograms(true,dir_zz2l2tau);
		vector<TH1F*> den4 = GetFRHistograms(false,dir_zz2l2tau);

		for(int i=0;i<(int)num.size();i++){
			TH1F* tmpNum2 = Subtract(num.at(i),num2.at(i));
			TH1F* tmpNum3 = Subtract(tmpNum2,num3.at(i));
			TH1F* tmpNum4 = Subtract(tmpNum3,num4.at(i));
			TH1F* tmpDen2 = Subtract(den.at(i),den2.at(i));
			TH1F* tmpDen3 = Subtract(tmpDen2,den3.at(i));
			TH1F* tmpDen4 = Subtract(tmpDen3,den4.at(i));
			if(rebin) tmpNum4->Rebin(num.at(i)->GetXaxis()->GetNbins());
			if(rebin) tmpDen4->Rebin(den.at(i)->GetXaxis()->GetNbins());
			tmpNum4->Divide(tmpDen4);
			ratios.push_back(tmpNum4);
		}
	}else{
		vector<TH1F*> num2  = GetFRHistograms(true,dir_wz);      vector<TH1F*> den2  = GetFRHistograms(false,dir_wz);
		vector<TH1F*> num3  = GetFRHistograms(true,dir_ww);      vector<TH1F*> den3  = GetFRHistograms(false,dir_ww);
		vector<TH1F*> num4  = GetFRHistograms(true,dir_vg);      vector<TH1F*> den4  = GetFRHistograms(false,dir_vg);
		vector<TH1F*> num5  = GetFRHistograms(true,dir_top);     vector<TH1F*> den5  = GetFRHistograms(false,dir_top);
		vector<TH1F*> num6  = GetFRHistograms(true,dir_qcd);     vector<TH1F*> den6  = GetFRHistograms(false,dir_qcd);
		vector<TH1F*> num7  = GetFRHistograms(true,dir_qcdmu);   vector<TH1F*> den7  = GetFRHistograms(false,dir_qcdmu);
		vector<TH1F*> num8  = GetFRHistograms(true,dir_wjets);   vector<TH1F*> den8  = GetFRHistograms(false,dir_wjets);
		vector<TH1F*> num9  = GetFRHistograms(true,dir_zz2l2nu); vector<TH1F*> den9  = GetFRHistograms(false,dir_zz2l2nu);
		vector<TH1F*> num10 = GetFRHistograms(true,dir_www);     vector<TH1F*> den10 = GetFRHistograms(false,dir_www);
		//vector<TH1F*> num10 = GetFRHistograms(true,dir_vvv);     vector<TH1F*> den10 = GetFRHistograms(false,dir_vvv);

		for(int i=0;i<(int)num.size();i++){
			TH1F *tmpNum = (TH1F*)num.at(i)->Clone(); 
			TH1F *tmpDen = (TH1F*)den.at(i)->Clone(); 
			tmpNum->Add(num2.at(i)); 
			tmpNum->Add(num3.at(i)); tmpNum->Add(num4.at(i)); tmpNum->Add(num5.at(i));
			tmpNum->Add(num6.at(i)); tmpNum->Add(num7.at(i)); tmpNum->Add(num8.at(i)); tmpNum->Add(num9.at(i));
			tmpNum->Add(num10.at(i));
			tmpDen->Add(den2.at(i)); 
			tmpDen->Add(den3.at(i)); tmpDen->Add(den4.at(i)); tmpDen->Add(den5.at(i));
			tmpDen->Add(den6.at(i)); tmpDen->Add(den7.at(i)); tmpDen->Add(den8.at(i)); tmpDen->Add(den9.at(i));
			tmpDen->Add(den10.at(i));
			if(rebin) tmpNum->Rebin(num.at(i)->GetXaxis()->GetNbins());
			if(rebin) tmpDen->Rebin(den.at(i)->GetXaxis()->GetNbins());
			tmpNum->Divide(tmpDen);
			ratios.push_back(tmpNum);
		}
	}

	return ratios;
}

vector<TH2F*> FRAnalyzer::BuildSummaryHistograms(string s){

	TH2F * h1 = new TH2F(("h"+s+"_ee_OSelmu"+Idx).c_str(),"s",8,0,8,100,0,30);
	TH2F * h2 = new TH2F(("h"+s+"_ee_OSelha"+Idx).c_str(),"s",8,0,8,100,0,30);
	TH2F * h3 = new TH2F(("h"+s+"_ee_OSmuha"+Idx).c_str(),"s",8,0,8,100,0,30);
	TH2F * h4 = new TH2F(("h"+s+"_ee_OShaha"+Idx).c_str(),"s",8,0,8,100,0,30);
	TH2F * h5 = new TH2F(("h"+s+"_mumu_OSelmu"+Idx).c_str(),"s",8,0,8,100,0,30);
	TH2F * h6 = new TH2F(("h"+s+"_mumu_OSelha"+Idx).c_str(),"s",8,0,8,100,0,30);
	TH2F * h7 = new TH2F(("h"+s+"_mumu_OSmuha"+Idx).c_str(),"s",8,0,8,100,0,30);
	TH2F * h8 = new TH2F(("h"+s+"_mumu_OShaha"+Idx).c_str(),"s",8,0,8,100,0,30);

	vector<TH2F*> h2D;
	h2D.push_back(h1); h2D.push_back(h2); h2D.push_back(h3); h2D.push_back(h4);   
	h2D.push_back(h5); h2D.push_back(h6); h2D.push_back(h7); h2D.push_back(h8);   
	
	return h2D;
}

vector<TH2F*> FRAnalyzer::BuildShapeHistograms(string s, string cat){
	string name="";
	double bins[]={5.0,32.4,70.2,110.2,189.,306.,550.8,1785.}; //OLD BINNING

	if(s.compare("H")==0) name = shapeName[0];
	if(s.compare("A")==0) name = shapeName[1];
	h_eeemOS = new TH2F(("ee_OSelmu"+name+cat).c_str(),TString("ee_OSelmu")+name,nbx,minx,maxx,nbins,bins);
	h_eeetOS = new TH2F(("ee_OSelha"+name+cat).c_str(),TString("ee_OSelha")+name,nbx,minx,maxx,nbins,bins);
	h_eemtOS = new TH2F(("ee_OSmuha"+name+cat).c_str(),TString("ee_OSmuha")+name,nbx,minx,maxx,nbins,bins);
	h_eettOS = new TH2F(("ee_OShaha"+name+cat).c_str(),TString("ee_OShaha")+name,nbx,minx,maxx,nbins,bins);
	h_mmemOS = new TH2F(("mumu_OSelmu"+name+cat).c_str(),TString("mumu_OSelmu")+name,nbx,minx,maxx,nbins,bins);
	h_mmetOS = new TH2F(("mumu_OSelha"+name+cat).c_str(),TString("mumu_OSelha")+name,nbx,minx,maxx,nbins,bins);
	h_mmmtOS = new TH2F(("mumu_OSmuha"+name+cat).c_str(),TString("mumu_OSmuha")+name,nbx,minx,maxx,nbins,bins);
	h_mmttOS = new TH2F(("mumu_OShaha"+name+cat).c_str(),TString("mumu_OShaha")+name,nbx,minx,maxx,nbins,bins);
	h_eeemSS = new TH2F(("ee_SSelmu"+name+cat).c_str(),TString("ee_SSelmu")+name,nbx,minx,maxx,nbins,bins);
	h_eeetSS = new TH2F(("ee_SSelha"+name+cat).c_str(),TString("ee_SSelha")+name,nbx,minx,maxx,nbins,bins);
	h_eemtSS = new TH2F(("ee_SSmuha"+name+cat).c_str(),TString("ee_SSmuha")+name,nbx,minx,maxx,nbins,bins);
	h_eettSS = new TH2F(("ee_SShaha"+name+cat).c_str(),TString("ee_SShaha")+name,nbx,minx,maxx,nbins,bins);
	h_mmemSS = new TH2F(("mumu_SSelmu"+name+cat).c_str(),TString("mumu_SSelmu")+name,nbx,minx,maxx,nbins,bins);
	h_mmetSS = new TH2F(("mumu_SSelha"+name+cat).c_str(),TString("mumu_SSelha")+name,nbx,minx,maxx,nbins,bins);
	h_mmmtSS = new TH2F(("mumu_SSmuha"+name+cat).c_str(),TString("mumu_SSmuha")+name,nbx,minx,maxx,nbins,bins);
	h_mmttSS = new TH2F(("mumu_SShaha"+name+cat).c_str(),TString("mumu_SShaha")+name,nbx,minx,maxx,nbins,bins);

	vector<TH2F*> h2D;
	h2D.push_back(h_eeemOS); h2D.push_back(h_eeetOS); h2D.push_back(h_eemtOS); h2D.push_back(h_eettOS);
	h2D.push_back(h_mmemOS); h2D.push_back(h_mmetOS); h2D.push_back(h_mmmtOS); h2D.push_back(h_mmttOS);
	h2D.push_back(h_eeemSS); h2D.push_back(h_eeetSS); h2D.push_back(h_eemtSS); h2D.push_back(h_eettSS);
	h2D.push_back(h_mmemSS); h2D.push_back(h_mmetSS); h2D.push_back(h_mmmtSS); h2D.push_back(h_mmttSS);

	return h2D;
}

void FRAnalyzer::PrintFinalNumbers(int i,int j, vector<TH2F*> shape,vector<TH2F*> shape1, vector<TH2F*> shape2){
	TH1D *tmp  = shape[i]->ProjectionY(TString::Format("p_%s", shape[i]->GetName()),j+1,j+1);
	TH1D *tmp1 = shape1[i]->ProjectionY(TString::Format("p_%s", shape[i]->GetName()),j+1,j+1);
	TH1D *tmp2 = shape2[i]->ProjectionY(TString::Format("p_%s", shape[i]->GetName()),j+1,j+1);
	double integral = tmp->GetBinContent(0); double integral1 = tmp1->GetBinContent(0); double integral2 = tmp2->GetBinContent(0);
	double error = tmp->GetBinError(0); double error1 = tmp1->GetBinError(0); double error2 = tmp2->GetBinError(0);
	if(i==0) cout << "$eee\\mu$ OS & " 		      << setprecision(3) <<integral<<"$\\pm$"<<error<<"&"<<integral1<<"$\\pm$"<<error1<<"&"<<integral2<<"$\\pm$"<<error2<<"\\\\\\hline" <<endl;
	if(i==1) cout << "$eee\\tau_{h}$ OS & " 	      << setprecision(3) <<integral<<"$\\pm$"<<error<<"&"<<integral1<<"$\\pm$"<<error1<<"&"<<integral2<<"$\\pm$"<<error2<<"\\\\\\hline" <<endl;
	if(i==2) cout << "$ee\\mu\\tau_{h}$ OS & " 	      << setprecision(3) <<integral<<"$\\pm$"<<error<<"&"<<integral1<<"$\\pm$"<<error1<<"&"<<integral2<<"$\\pm$"<<error2<<"\\\\\\hline" <<endl;
	if(i==3) cout << "$ee\\tau_{h}\\tau_{h}$ OS & "       << setprecision(3) <<integral<<"$\\pm$"<<error<<"&"<<integral1<<"$\\pm$"<<error1<<"&"<<integral2<<"$\\pm$"<<error2<<"\\\\\\hline" <<endl;
	if(i==4) cout << "$\\mu\\mu e\\mu$ OS &" 	      << setprecision(3) <<integral<<"$\\pm$"<<error<<"&"<<integral1<<"$\\pm$"<<error1<<"&"<<integral2<<"$\\pm$"<<error2<<"\\\\\\hline" <<endl;
	if(i==5) cout << "$\\mu\\mu e\\tau_{h}$ OS & " 	      << setprecision(3) <<integral<<"$\\pm$"<<error<<"&"<<integral1<<"$\\pm$"<<error1<<"&"<<integral2<<"$\\pm$"<<error2<<"\\\\\\hline" <<endl;
	if(i==6) cout << "$\\mu\\mu\\mu\\tau_{h}$ OS & "      << setprecision(3) <<integral<<"$\\pm$"<<error<<"&"<<integral1<<"$\\pm$"<<error1<<"&"<<integral2<<"$\\pm$"<<error2<<"\\\\\\hline" <<endl;
	if(i==7) cout << "$\\mu\\mu\\tau_{h}\\tau_{h}$ OS & " << setprecision(3) <<integral<<"$\\pm$"<<error<<"&"<<integral1<<"$\\pm$"<<error1<<"&"<<integral2<<"$\\pm$"<<error2<<"\\\\\\hline" <<endl;

} 

void FRAnalyzer::Printing(int i,int j, double integral01, double integral10, double integral00, double uCR01, double uCR10, double uCR00){

	double finalError = TMath::Sqrt(uCR01*uCR01+uCR10*uCR10+uCR00*uCR00);
	if((integral01+integral10)>=integral00) 
		if(i==0){ cout << "$eee\\mu$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral01+integral10-integral00<<"$\\pm$"<< finalError <<"\\\\" << endl;}

	if((integral01+integral10)<integral00) 
		if(i==0){ cout << "$eee\\mu$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral00<<"$\\pm$" << uCR00 <<"\\\\" << endl;}

	if((integral01+integral10)>=integral00) 
		if(i==1){ cout << "$eee\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral01+integral10-integral00<<"$\\pm$"<<finalError <<"\\\\" << endl;}

	if((integral01+integral10)<integral00) 
		if(i==1){ cout << "$eee\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral00<<"$\\pm$" << uCR00 <<"\\\\" << endl;}

	if((integral01+integral10)>=integral00) 
		if(i==2){ cout << "$ee\\mu\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral01+integral10-integral00<<"$\\pm$"<<finalError <<"\\\\" << endl;}

	if((integral01+integral10)<integral00) 
		if(i==2){ cout << "$ee\\mu\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral00<<"$\\pm$" << uCR00 <<"\\\\" << endl;}

	if((integral01+integral10)>=integral00) 
		if(i==3){ cout << "$ee\\tau_{h}\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral01+integral10-integral00<<"$\\pm$"<<finalError <<"\\\\" << endl;}

	if((integral01+integral10)<integral00) 
		if(i==3){ cout << "$ee\\tau_{h}\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral00<<"$\\pm$" << uCR00 <<"\\\\" << endl;}

	if((integral01+integral10)>=integral00) 
		if(i==4){ cout << "$\\mu\\mu e\\mu$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral01+integral10-integral00<<"$\\pm$"<<finalError <<"\\\\" << endl;}

	if((integral01+integral10)<integral00) 
		if(i==4){ cout << "$\\mu\\mu e\\mu$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral00<<"$\\pm$" << uCR00 <<"\\\\" << endl;}

	if((integral01+integral10)>=integral00) 
		if(i==5){ cout << "$\\mu\\mu e\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral01+integral10-integral00<<"$\\pm$"<<finalError <<"\\\\" << endl;}

	if((integral01+integral10)<integral00) 
		if(i==5){ cout << "$\\mu\\mu e\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral00<<"$\\pm$" << uCR00 <<"\\\\" << endl;}

	if((integral01+integral10)>=integral00) 
		if(i==6){ cout << "$\\mu\\mu\\mu\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral01+integral10-integral00<<"$\\pm$"<<finalError <<"\\\\" << endl;}

	if((integral01+integral10)<integral00) 
		if(i==6){ cout << "$\\mu\\mu\\mu\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral00<<"$\\pm$" << uCR00 <<"\\\\" << endl;}

	if((integral01+integral10)>=integral00) 
		if(i==7){ cout << "$\\mu\\mu\\tau_{h}\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral01+integral10-integral00<<"$\\pm$"<<finalError <<"\\\\" << endl;}

	if((integral01+integral10)<integral00) 
		if(i==7){ cout << "$\\mu\\mu\\tau_{h}\\tau_{h}$ OS & "  << integral01<<"$\\pm$" << uCR01<<" & " <<integral10<<"$\\pm$" << uCR10<<" & "  << integral00<<"$\\pm$" << uCR00 <<" & "  << integral00<<"$\\pm$" << uCR00 << endl;}

}

vector<TH2F*> FRAnalyzer::GetStraightMC(){
	vector<string> dirs;
	vector<TH2F*> shapes  = BuildShapeHistograms("A","");
	dirs.push_back(dir_ww); 
	dirs.push_back(dir_vg); 
	dirs.push_back(dir_top); 
	dirs.push_back(dir_qcd); 
	dirs.push_back(dir_qcdmu); 
	dirs.push_back(dir_wjets); 
	dirs.push_back(dir_wz); 
	dirs.push_back(dir_zz2l2nu);
	dirs.push_back(dir_www); 
	dirs.push_back(dir_dy);

	for(int i=0;i<(int)dirs.size();i++){

		TDirectoryFile *h = (TDirectoryFile*)file1->GetDirectory(dirs.at(i).c_str());
		TKey *key1 = h->GetKey((finalStates.at(0)+shapeName.at(1)).c_str());
		TKey *key2 = h->GetKey((finalStates.at(1)+shapeName.at(1)).c_str());
		TKey *key3 = h->GetKey((finalStates.at(2)+shapeName.at(1)).c_str());
		TKey *key4 = h->GetKey((finalStates.at(3)+shapeName.at(1)).c_str());
		TKey *key5 = h->GetKey((finalStates.at(4)+shapeName.at(1)).c_str());
		TKey *key6 = h->GetKey((finalStates.at(5)+shapeName.at(1)).c_str());
		TKey *key7 = h->GetKey((finalStates.at(6)+shapeName.at(1)).c_str());
		TKey *key8 = h->GetKey((finalStates.at(7)+shapeName.at(1)).c_str());
		TKey *key9 = h->GetKey((finalStates.at(8)+shapeName.at(1)).c_str());
		TKey *key10 = h->GetKey((finalStates.at(9)+shapeName.at(1)).c_str());
		TKey *key11 = h->GetKey((finalStates.at(10)+shapeName.at(1)).c_str());
		TKey *key12 = h->GetKey((finalStates.at(11)+shapeName.at(1)).c_str());
		TKey *key13 = h->GetKey((finalStates.at(12)+shapeName.at(1)).c_str());
		TKey *key14 = h->GetKey((finalStates.at(13)+shapeName.at(1)).c_str());
		TKey *key15 = h->GetKey((finalStates.at(14)+shapeName.at(1)).c_str());
		TKey *key16 = h->GetKey((finalStates.at(15)+shapeName.at(1)).c_str());


		if(key1!=0) shapes[0]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(0)+shapeName.at(1)).c_str()));
		if(key2!=0) shapes[1]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(1)+shapeName.at(1)).c_str()));
		if(key3!=0) shapes[2]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(2)+shapeName.at(1)).c_str()));
		if(key4!=0) shapes[3]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(3)+shapeName.at(1)).c_str()));
		if(key5!=0) shapes[4]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(4)+shapeName.at(1)).c_str()));
		if(key6!=0) shapes[5]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(5)+shapeName.at(1)).c_str()));
		if(key7!=0) shapes[6]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(6)+shapeName.at(1)).c_str()));
		if(key8!=0) shapes[7]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(7)+shapeName.at(1)).c_str()));
		if(key9!=0) shapes[8]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(8)+shapeName.at(1)).c_str()));
		if(key10!=0) shapes[9]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(9)+shapeName.at(1)).c_str()));
		if(key11!=0) shapes[10]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(10)+shapeName.at(1)).c_str()));
		if(key12!=0) shapes[11]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(11)+shapeName.at(1)).c_str()));
		if(key13!=0) shapes[12]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(12)+shapeName.at(1)).c_str()));
		if(key14!=0) shapes[13]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(13)+shapeName.at(1)).c_str()));
		if(key15!=0) shapes[14]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(14)+shapeName.at(1)).c_str()));
		if(key16!=0) shapes[15]->Add((TH2F*)file1->Get((dirs.at(i)+"/"+finalStates.at(15)+shapeName.at(1)).c_str()));
	}
	return shapes;
}

pair<double,double> FRAnalyzer::LogNormalFit(TH1 *h1){

	double xmin = h1->GetXaxis()->GetXmin();
	double xmax = h1->GetXaxis()->GetXmax();
	tmpFitFun = NULL;
	tmpFitFun = GetFitFunction();
	//TF1 * f1 = new TF1("f1","[0]*ROOT::Math::lognormal_pdf(x,[1],[2])",xmin,xmax);
	//setting initial parameters, extracting those from the histogram 
	tmpFitFun->SetRange(xmin,xmax);
	double prob[] = {0.5}; 
	double q[1]; 
	h1->GetQuantiles(1,q,prob);
        double median = q[0];
        double par0 = h1->Integral();
        double par1 = log(median);
        double par2 = log(1+(h1->GetRMS()/median));
		tmpFitFun->SetParameter(0,par0); 
		tmpFitFun->SetParameter(1,par1); 
		tmpFitFun->SetParameter(2,par2); 
		h1->Fit(tmpFitFun,"R");
        pair<double,double> estimate;
        estimate=make_pair( exp(tmpFitFun->GetParameter(1)), (exp(tmpFitFun->GetParameter(1)))*(exp(tmpFitFun->GetParameter(2)) - 1) ); 
        return estimate;
}

//perform the optimal binning
double FRAnalyzer::AutomatizeBinning(double content, double error, bool max){
	double x = (max==true)?(content+10*error):(content-10*error);
	return x;
}

TFile* FRAnalyzer::GetToyExperiments(vector<TH2F*>sCR01,vector<TH2F*>sCR10,vector<TH2F*>sCR00,vector<TH2F*>hCR01,vector<TH2F*>hCR10,vector<TH2F*>hCR00,
		vector<TH2F*>sDataFit,vector<TH2F*>sDataHis,vector<TH2F*>sMC,bool useBayesian){

	TFile *fileToys = NULL; 
	vector<TH2F*> hSB_OS = BuildSummaryHistograms("");
	
	for(int ch=0;ch<(int)finalStates.size()/2;ch++){
		vector<double> countOSCC, shapeOSCC, errorOSCC;
		vector<double> countSSCC, shapeSSCC, errorSSCC;
		double MCmeanOS, MCerrOS, MCmeanSS, MCerrSS;
		double we = (GetFRWeights())[0]; double wm = (GetFRWeights())[1]; double wt = (GetFRWeights())[2];
		countOSCC.push_back(hCR10[ch]->GetBinContent(m_Index+1,0)); 
		countOSCC.push_back(hCR01[ch]->GetBinContent(m_Index+1,0)); 
		countOSCC.push_back(hCR00[ch]->GetBinContent(m_Index+1,0));
		shapeOSCC.push_back(sCR10[ch]->GetBinContent(m_Index+1,0)); 
		shapeOSCC.push_back(sCR01[ch]->GetBinContent(m_Index+1,0)); 
		shapeOSCC.push_back(sCR00[ch]->GetBinContent(m_Index+1,0));
		errorOSCC.push_back(sCR10[ch]->GetBinError(m_Index+1,0)); 
		errorOSCC.push_back(sCR01[ch]->GetBinError(m_Index+1,0)); 
		errorOSCC.push_back(sCR00[ch]->GetBinError(m_Index+1,0));
		countSSCC.push_back(hCR10[ch]->GetBinContent(m_Index+1,9));
		countSSCC.push_back(hCR01[ch]->GetBinContent(m_Index+1,9));
		countSSCC.push_back(hCR00[ch]->GetBinContent(m_Index+1,9));
		shapeSSCC.push_back(sCR10[ch]->GetBinContent(m_Index+1,9));
		shapeSSCC.push_back(sCR01[ch]->GetBinContent(m_Index+1,9));
		shapeSSCC.push_back(sCR00[ch]->GetBinContent(m_Index+1,9));
		errorSSCC.push_back(sCR10[ch]->GetBinError(m_Index+1,9)); 
		errorSSCC.push_back(sCR01[ch]->GetBinError(m_Index+1,9)); 
		errorSSCC.push_back(sCR00[ch]->GetBinError(m_Index+1,9));
		MCmeanOS = sMC[ch]->GetBinContent(m_Index+1,0); MCerrOS = sMC[ch]->GetBinError(m_Index+1,0);
		MCmeanSS = sMC[ch]->GetBinContent(m_Index+1,9); MCerrSS = sMC[ch]->GetBinError(m_Index+1,9);


                cout << "INPUTS " << finalStates[ch] << " CR10: " << countOSCC[0] << " CR01:" << countOSCC[1] << " CR00:" << countOSCC[2] << " MC: " << MCmeanOS << endl; 
		if((ch==0||ch==4)) cout << "EM WEIGHTS CR10: " << wm << " CR01 " << we << " CR00 " << we*wm << endl;
		if((ch==1||ch==5)) cout << "ET WEIGHTS CR10: " << wt << " CR01 " << we << " CR00 " << we*wt << endl;
		if((ch==2||ch==6)) cout << "MT WEIGHTS CR10: " << wt << " CR01 " << wm << " CR00 " << wm*wt << endl;
	        if((ch==3||ch==7)) cout << "TT WEIGHTS CR10: " << wt << " CR01 " << wt << " CR00 " << wt*wt << endl;

                //double cSSCC = shapeSSCC[0]+shapeSSCC[1]-shapeSSCC[2]-MCmeanSS;
                //double eOSCC = sqrt(pow(errorOSCC[0],2)+pow(errorOSCC[1],2)+pow(errorOSCC[2],2)+pow(MCerrOS,2));
                //double eSSCC = sqrt(pow(errorSSCC[0],2)+pow(errorSSCC[1],2)+pow(errorSSCC[2],2)+pow(MCerrSS,2));
		//double xmaxOS = AutomatizeBinning(cOSCC,eOSCC,true);
		//double xmaxSS = AutomatizeBinning(cSSCC,eSSCC,true);

        	TH1F * hCC_OS  = new TH1F(("hOS"+finalStates[ch]+Idx).c_str(),"hfinal",100,0,30); 
		TH1F * hCC_SS  = new TH1F(("hSS"+finalStates[ch]+Idx).c_str(),"hfinal",100,0,30);

		for(int toy=0;toy<m_Ntoy;toy++){
			vector<double> muCC; vector<double> weightsCC; //vector<int> NCC;                 
			double mcxOSCC = (double)gRandom->Gaus(MCmeanOS,MCerrOS);
			double mcxSSCC = (double)gRandom->Gaus(MCmeanSS,MCerrSS);
			for(int j=0;j<(int)countOSCC.size();j++){
				f1 = NULL;
                                f1 = GetLikelihoodFunction(countOSCC.at(j));
				f1->SetParameter(0,countOSCC.at(j));
				double mux;
                                if(countOSCC.at(j)>0){(useBayesian)?(mux=f1->GetRandom()):(mux=gRandom->Poisson(countOSCC.at(j)));}
                                else if(countOSCC.at(j)==0) mux=f1->GetRandom();
				muCC.push_back(mux);  
				if(countOSCC.at(j)==0){
					if((ch==0||ch==4) && j==0) weightsCC.push_back(wm); 
					if((ch==0||ch==4) && j==1) weightsCC.push_back(we); 
					if((ch==0||ch==4) && j==2) weightsCC.push_back(we*wm); 
					if((ch==1||ch==5) && j==0) weightsCC.push_back(wt); 
					if((ch==1||ch==5) && j==1) weightsCC.push_back(we); 
					if((ch==1||ch==5) && j==2) weightsCC.push_back(we*wt); 
					if((ch==2||ch==6) && j==0) weightsCC.push_back(wt); 
					if((ch==2||ch==6) && j==1) weightsCC.push_back(wm); 
					if((ch==2||ch==6) && j==2) weightsCC.push_back(wm*wt); 
					if((ch==3||ch==7) && j<2)  weightsCC.push_back(wt); 
					if((ch==3||ch==7) && j==2) weightsCC.push_back(wt*wt); 
				}
				else{ 
					//muCC.push_back((double)gRandom->Poisson(countOSCC[j])); 
					weightsCC.push_back(shapeOSCC[j]/countOSCC[j]); 
				} 
			}
			
			for(int j=0;j<(int)countSSCC.size();j++){
				f1 = NULL;
				//double mux = hPDF0->GetRandom();
                                f1 = GetLikelihoodFunction(countSSCC.at(j));
				f1->SetParameter(0,countSSCC.at(j));
				double mux;
                                if(countSSCC.at(j)>0){(useBayesian)?(mux=f1->GetRandom()):(mux=gRandom->Poisson(countSSCC.at(j)));}
                                else if(countSSCC.at(j)==0) mux=f1->GetRandom();
				muCC.push_back(mux); //NCC.push_back(gRandom->Poisson(mux)); 
				if(countSSCC.at(j)==0){ 
					if((ch==0||ch==4) && j==0) weightsCC.push_back(wm); 
					if((ch==0||ch==4) && j==1) weightsCC.push_back(we); 
					if((ch==0||ch==4) && j==2) weightsCC.push_back(we*wm); 
					if((ch==1||ch==5) && j==0) weightsCC.push_back(wt); 
					if((ch==1||ch==5) && j==1) weightsCC.push_back(we); 
					if((ch==1||ch==5) && j==2) weightsCC.push_back(we*wt); 
					if((ch==2||ch==6) && j==0) weightsCC.push_back(wt); 
					if((ch==2||ch==6) && j==1) weightsCC.push_back(wm); 
					if((ch==2||ch==6) && j==2) weightsCC.push_back(wm*wt); 
					if((ch==3||ch==7) && j<2)  weightsCC.push_back(wt); 
					if((ch==3||ch==7) && j==2) weightsCC.push_back(wt*wt); 
				}
				else{ 
					//muCC.push_back((double)gRandom->Poisson(countSSCC[j])); 
					weightsCC.push_back(shapeSSCC[j]/countSSCC[j]);
				}
			}
			double cc_OS = muCC[0]*weightsCC[0]+muCC[1]*weightsCC[1]-muCC[2]*weightsCC[2]-mcxOSCC;
			double cc_SS = muCC[3]*weightsCC[3]+muCC[4]*weightsCC[4]-muCC[5]*weightsCC[5]-mcxSSCC;
			//cout << " cut & count - FINAL STATE " << finalStates[ch] << endl; 
			//cout << NCC[0] <<" * "<< weightsCC[0] << " + " << NCC[1] << " * " << weightsCC[1] << " - " << NCC[2] << " * " << weightsCC[2] << " - " << mcxOSCC << " = " << cc_NOS << endl;
			if(cc_OS<0) cc_OS=0;
			if(cc_SS<0) cc_SS=0;
			hCC_OS->Fill(cc_OS); 
			hCC_SS->Fill(cc_SS); 

			//cout << "WEIGHTTTTTTTTTTTTTTT " << weight << endl;

			for(int bin=0;bin<8;bin++){
				stringstream cbin; string sbin;
				cbin<<bin; cbin>>sbin;
				vector<double> muSB; vector<double> NSB; vector<double>weightsSB;                  
				vector<double> countOSSB, shapeOSSB, errorOSSB;
				vector<double> countSSSB, shapeSSSB, errorSSSB;
				countOSSB.push_back(hCR10[ch]->GetBinContent(m_Index+1,bin+1)); 
				countOSSB.push_back(hCR01[ch]->GetBinContent(m_Index+1,bin+1)); 
				countOSSB.push_back(hCR00[ch]->GetBinContent(m_Index+1,bin+1));
				shapeOSSB.push_back(sCR10[ch]->GetBinContent(m_Index+1,bin+1)); 
				shapeOSSB.push_back(sCR01[ch]->GetBinContent(m_Index+1,bin+1)); 
				shapeOSSB.push_back(sCR00[ch]->GetBinContent(m_Index+1,bin+1));
				countSSSB.push_back(hCR10[ch+8]->GetBinContent(m_Index+1,bin+1)); 
				countSSSB.push_back(hCR01[ch+8]->GetBinContent(m_Index+1,bin+1)); 
				countSSSB.push_back(hCR00[ch+8]->GetBinContent(m_Index+1,bin+1));
				shapeSSSB.push_back(sCR10[ch+8]->GetBinContent(m_Index+1,bin+1)); 
				shapeSSSB.push_back(sCR01[ch+8]->GetBinContent(m_Index+1,bin+1)); 
				shapeSSSB.push_back(sCR00[ch+8]->GetBinContent(m_Index+1,bin+1));
				double mcxOSSB = (double)gRandom->Gaus(sMC[ch]->GetBinContent(m_Index+1,bin+1));
				double mcxSSSB = (double)gRandom->Gaus(sMC[ch+8]->GetBinContent(m_Index+1,bin+1));
				//cout << " shape based - bin " << bin << " OS - FINAL STATE " << finalStates[ch] << countOSSB[0]<<"/"<<countOSSB[1]<<"/"<<countOSSB[2]<<endl; 
				//cout << " shape based - bin " << bin << " SS - FINAL STATE " << finalStates[ch] << countSSSB[0]<<"/"<<countSSSB[1]<<"/"<<countSSSB[2]<<endl; 
				for(int j=0;j<(int)countOSSB.size();j++){
					//double mux = hPDF0->GetRandom();
					f1 = NULL;
					f1 = GetLikelihoodFunction(countOSSB.at(j));
					f1->SetParameter(0,countOSSB.at(j));
					double mux;
					if(countOSSB.at(j)>0){(useBayesian)?(mux=f1->GetRandom()):(mux=gRandom->Poisson(countOSSB.at(j)));}
					else if(countOSSB.at(j)==0) mux=f1->GetRandom();
					muSB.push_back(mux); 
					if(countOSSB[j]==0){ 
						if((ch==0||ch==4) && j==0) weightsSB.push_back(wm); 
						if((ch==0||ch==4) && j==1) weightsSB.push_back(we); 
						if((ch==0||ch==4) && j==2) weightsSB.push_back(we*wm); 
						if((ch==1||ch==5) && j==0) weightsSB.push_back(wt); 
						if((ch==1||ch==5) && j==1) weightsSB.push_back(we); 
						if((ch==1||ch==5) && j==2) weightsSB.push_back(we*wt); 
						if((ch==2||ch==6) && j==0) weightsSB.push_back(wt); 
						if((ch==2||ch==6) && j==1) weightsSB.push_back(wm); 
						if((ch==2||ch==6) && j==2) weightsSB.push_back(wm*wt); 
						if((ch==3||ch==7) && j<2)  weightsSB.push_back(wt); 
						if((ch==3||ch==7) && j==2) weightsSB.push_back(wt*wt); 
					}
					else{ 
						weightsSB.push_back(shapeOSSB[j]/countOSSB[j]);
					}
				}
				for(int j=0;j<(int)countSSSB.size();j++){
					f1 = NULL;
					f1 = GetLikelihoodFunction(countSSSB.at(j));
					f1->SetParameter(0,countSSSB.at(j));
					double mux;
					if(countSSSB.at(j)>0){(useBayesian)?(mux=f1->GetRandom()):(mux=gRandom->Poisson(countSSSB.at(j)));}
					else if(countSSSB.at(j)==0) mux=f1->GetRandom();
					muSB.push_back(mux); 
					if(countSSSB[j]==0){ 
						if((ch==0||ch==4) && j==0) weightsSB.push_back(wm); 
						if((ch==0||ch==4) && j==1) weightsSB.push_back(we); 
						if((ch==0||ch==4) && j==2) weightsSB.push_back(we*wm); 
						if((ch==1||ch==5) && j==0) weightsSB.push_back(wt); 
						if((ch==1||ch==5) && j==1) weightsSB.push_back(we); 
						if((ch==1||ch==5) && j==2) weightsSB.push_back(we*wt); 
						if((ch==2||ch==6) && j==0) weightsSB.push_back(wt); 
						if((ch==2||ch==6) && j==1) weightsSB.push_back(wm); 
						if((ch==2||ch==6) && j==2) weightsSB.push_back(wm*wt); 
						if((ch==3||ch==7) && j<2)  weightsSB.push_back(wt); 
						if((ch==3||ch==7) && j==2) weightsSB.push_back(wt*wt); 
					}
					else{ 
						weightsSB.push_back(shapeSSSB[j]/countSSSB[j]);
					}
				}
				double sbTot = (cc_OS/(cc_OS+cc_SS)) * 
					( (muSB[0]*weightsSB[0])+
					 (muSB[3]*weightsSB[3])+
					 (muSB[1]*weightsSB[1])+
					 (muSB[4]*weightsSB[4])-
					 (muSB[2]*weightsSB[2])-
					 (muSB[5]*weightsSB[5])-
					 mcxOSSB-mcxSSSB );

				if(sbTot<0) sbTot=0;
				hSB_OS[ch]->Fill(bin,sbTot);
			}//close the loop on the other bins
		}//close toy loop

		pair<double,double> ccOS = LogNormalFit(hCC_OS); 
	        pair<double,double> ccSS = LogNormalFit(hCC_SS);
                vector<pair<double,double> > sbOS;

                string name="";
                if(useBayesian) name += "_bayesian";
                	
		fileToys = new TFile((m_Dir+"/toysLogNormal"+Idx+"_"+finalStates[ch]+name+".root").c_str(), "RECREATE");

		for(int bin=0;bin<8;bin++){
			TH1D *tmp  = (TH1D*)hSB_OS[ch]->ProjectionY(TString::Format("p_%s",hSB_OS[ch]->GetName())+Idx+finalStates[ch],bin+1,bin+1);
			sbOS.push_back(LogNormalFit(tmp));
			stringstream cbin; string Bin;
			cbin<<bin;
			cbin>>Bin;
			tmp->Write(("f"+Idx+"_LogNormalMu_SB_"+finalStates[ch]+"_OS_bin"+Bin).c_str());
                }
		hCC_OS->Write(("f"+Idx+"_LogNormalMu_CC_"+finalStates[ch]+"_OS").c_str());
		hCC_SS->Write(("f"+Idx+"_LogNormalMu_CC_"+finalStates[ch]+"_SS").c_str());

		sDataHis[ch]->SetBinContent(m_Index+1,0,hCC_OS->GetMean()); sDataHis[ch]->SetBinError(m_Index+1,0,hCC_OS->GetRMS());
		sDataFit[ch]->SetBinContent(m_Index+1,0,ccOS.first); sDataFit[ch]->SetBinError(m_Index+1,0,ccOS.second);
		sDataHis[ch]->SetBinContent(m_Index+1,9,hCC_SS->GetMean()); sDataHis[ch]->SetBinError(m_Index+1,9,hCC_SS->GetRMS());
		sDataFit[ch]->SetBinContent(m_Index+1,9,ccSS.first); sDataFit[ch]->SetBinError(m_Index+1,9,ccSS.second);

		double SF = hCC_OS->GetMean() / (hCC_OS->GetMean()+hCC_SS->GetMean());

		for(int bin=0;bin<8;bin++){
			TH1D *tmp  = (TH1D*)hSB_OS[ch]->ProjectionY(TString::Format("p_%s",hSB_OS[ch]->GetName())+Idx+finalStates[ch],bin+1,bin+1);
			sDataHis[ch]->SetBinContent(m_Index+1,bin+1,tmp->GetMean()); sDataHis[ch]->SetBinError(m_Index+1,bin+1,tmp->GetRMS());
			sDataFit[ch]->SetBinContent(m_Index+1,bin+1,(sbOS.at(bin)).first); sDataFit[ch]->SetBinError(m_Index+1,bin+1,(sbOS.at(bin)).second);
		}
	}//close the loop on the channels

        return fileToys;
	//hFinal_N->Scale(1/hFinal_N->Integral());         hFinal_mu->Scale(1/hFinal_mu->Integral());
	//hFinal_N_dyn->Scale(1/hFinal_N_dyn->Integral()); hFinal_mu_dyn->Scale(1/hFinal_mu_dyn->Integral());
	//TFile *fileToys  = new TFile((m_Dir+"/toys"+Idx+"_"+finalStates[ch]+".root").c_str(), "RECREATE");
	//hFinal_N->Write(("f"+Idx+"_N_"+finalStates[ch]).c_str());
	//hFinal_N_dyn->Write(("f"+Idx+"_Nd_"+finalStates[ch]).c_str());
	//hFinal_mu->Write(("f"+Idx+"_mu_"+finalStates[ch]).c_str());
	//hFinal_mu_dyn->Write(("f"+Idx+"_mud_"+finalStates[ch]).c_str());
}


void FRAnalyzer::Estimate(string dir,vector<TH1F*>ratios, 
		double elIso, double muIso, double tauIso, double sumPt, int cut_index,
		vector<TH2F*> h2DAc_CR01, vector<TH2F*> h2DAc_CR10, vector<TH2F*> h2DAc_CR00, vector<TH2F*> h2DAc_CR11,
		vector<TH2F*> h2DA_CR01,  vector<TH2F*> h2DA_CR10,  vector<TH2F*> h2DA_CR00,  vector<TH2F*> h2DA_CR11){

	TTree* tree                           = (TTree*)file2->Get((dir+"/CandTree").c_str());
	TTree* treeF                          = (TTree*)file2->Get((dir+"/CandTree_PWeight").c_str());

	UInt_t          treeEventId             = 0;   tree->SetBranchAddress("eventId", &treeEventId );
	UInt_t          treeLumiId              = 0;   tree->SetBranchAddress("lumiId", &treeLumiId );
	int             treeHiggsId             = 0;   tree->SetBranchAddress("higgsId", &treeHiggsId );
	int             treeZId                 = 0;   tree->SetBranchAddress("zId"    , &treeZId     );
	int             treeLeg1Id              = 0;   tree->SetBranchAddress("leg1Id" , &treeLeg1Id  );
	float           treeLeg1Iso             = 0;   tree->SetBranchAddress("leg1Iso" , &treeLeg1Iso);
	float           treeLeg1ClosestJetPt    = 0;   tree->SetBranchAddress("leg1ClosestJetPt" , &treeLeg1ClosestJetPt  );
	float           treeLeg1ClosestJetEta   = 0;   tree->SetBranchAddress("leg1ClosestJetEta" , &treeLeg1ClosestJetEta  );
	float           treeLeg1Pt              = 0;   tree->SetBranchAddress("leg1Pt" , &treeLeg1Pt  );
	float           treeLeg1Eta             = 0;   tree->SetBranchAddress("leg1Eta", &treeLeg1Eta );
	int             treeLeg1LepIDloose      = 0;   tree->SetBranchAddress("leg1LepIDloose" , &treeLeg1LepIDloose  );
	int             treeLeg1LepIDtight      = 0;   tree->SetBranchAddress("leg1LepIDtight" , &treeLeg1LepIDtight  );
	float           treeLeg1DeltaR          = 0;   tree->SetBranchAddress("leg1DeltaR" , &treeLeg1DeltaR  );
	int             treeLeg2Id              = 0;   tree->SetBranchAddress("leg2Id" , &treeLeg2Id  );
	float           treeLeg2Iso             = 0;   tree->SetBranchAddress("leg2Iso" , &treeLeg2Iso  );
	float           treeLeg2ClosestJetPt    = 0;   tree->SetBranchAddress("leg2ClosestJetPt" , &treeLeg2ClosestJetPt  );
	float           treeLeg2ClosestJetEta   = 0;   tree->SetBranchAddress("leg2ClosestJetEta" , &treeLeg2ClosestJetEta  );
	float           treeLeg2Pt              = 0;   tree->SetBranchAddress("leg2Pt" , &treeLeg2Pt  );
	float           treeLeg2Eta             = 0;   tree->SetBranchAddress("leg2Eta", &treeLeg2Eta );
	int             treeLeg2LepIDloose      = 0;   tree->SetBranchAddress("leg2LepIDloose" , &treeLeg2LepIDloose  );
	int             treeLeg2LepIDtight      = 0;   tree->SetBranchAddress("leg2LepIDtight" , &treeLeg2LepIDtight  );
	float           treeLeg2DeltaR          = 0;   tree->SetBranchAddress("leg2DeltaR" , &treeLeg2DeltaR  );
	UInt_t          treeTau1Loose3          = 0;   tree->SetBranchAddress("tau1Loose3", &treeTau1Loose3);
	UInt_t          treeTau1Medium3         = 0;   tree->SetBranchAddress("tau1Medium3", &treeTau1Medium3);
	UInt_t          treeTau1LooseMVA        = 0;   tree->SetBranchAddress("tau1LooseMVA", &treeTau1LooseMVA);
	UInt_t          treeTau1MediumMVA       = 0;   tree->SetBranchAddress("tau1MediumMVA", &treeTau1MediumMVA);
	UInt_t          treeTau2Loose3          = 0;   tree->SetBranchAddress("tau2Loose3", &treeTau2Loose3);
	UInt_t          treeTau2Medium3         = 0;   tree->SetBranchAddress("tau2Medium3", &treeTau2Medium3);
	UInt_t          treeTau2LooseMVA        = 0;   tree->SetBranchAddress("tau2LooseMVA", &treeTau2LooseMVA);
	UInt_t          treeTau2MediumMVA       = 0;   tree->SetBranchAddress("tau2MediumMVA", &treeTau2MediumMVA);
	float           treeSumPt               = 0;   tree->SetBranchAddress("sumPt", &treeSumPt );
	float           treeWeight              = 0;   tree->SetBranchAddress("weight", &treeWeight );
	float           treeLepEff1             = 0;   tree->SetBranchAddress("lepEff1", &treeLepEff1 );
	float           treeLepEff2             = 0;   tree->SetBranchAddress("lepEff2", &treeLepEff2 );
	float           treeHvisMass            = 0;   tree->SetBranchAddress("HvisMass", &treeHvisMass );
	float           treeAvisMass            = 0;   tree->SetBranchAddress("AvisMass", &treeAvisMass );
	float           treeHsvfMass            = 0;   tree->SetBranchAddress("HsvfMass", &treeHsvfMass );
	float           treeAsvfMass            = 0;   tree->SetBranchAddress("AsvfMass", &treeAsvfMass );

	float           treeFWeight             = 0;   treeF->SetBranchAddress("plotterWeight", &treeFWeight );

	tree->AddFriend(treeF);
	int  sign=+1;
	TH1F *eleFR=NULL,  *muoFR=NULL, *tauFR=NULL; 
	if(elIso==0.1)  eleFR=(TH1F*)ratios[0]->Clone();
	if(elIso==0.2)  eleFR=(TH1F*)ratios[1]->Clone(); 
	if(elIso==0.3)  eleFR=(TH1F*)ratios[2]->Clone();
	if(muIso==0.1)  muoFR=(TH1F*)ratios[3]->Clone();
	if(muIso==0.2)  muoFR=(TH1F*)ratios[4]->Clone();
	if(muIso==0.3)  muoFR=(TH1F*)ratios[5]->Clone();
	if(tauIso==0.0) tauFR=(TH1F*)ratios[6]->Clone();
	if(tauIso==1.0) tauFR=(TH1F*)ratios[7]->Clone();

	for(Int_t ev=0;ev<tree->GetEntries();ev++){
		tree->GetEntry(ev);
		double weight;
		double lepEff1=treeLepEff1;
		double lepEff2=treeLepEff2;
		if(dir.find("data") != string::npos) weight=1;
		else weight=treeWeight*treeFWeight;
		double f1jOS=0; double f1jSS=0;
		double f2jOS=0; double f2jSS=0;

		if(treeSumPt<sumPt) continue;
		if(debug) cout << "***** LT cut is --> " << sumPt << " - the event has LT > --> " << treeSumPt << endl; 
		if(debug) cout << "***** The event " << treeEventId << " has passed the LT cut " << endl; 
		bool TAU_LT=false;
		bool TAU1_TT=false;
		bool TAU2_TT=false;

		if(tauIso==0.0){
			if((abs(treeHiggsId)==165 || abs(treeHiggsId)==195) && (abs(treeLeg1Id)<15) && treeTau2Loose3==1) TAU_LT=true;
			if((abs(treeHiggsId)==165 || abs(treeHiggsId)==195) && (abs(treeLeg2Id)<15) && treeTau1Loose3==1) TAU_LT=true;
			if(abs(treeHiggsId)==225){if(treeLeg1Pt>treeLeg2Pt && treeTau1Loose3==1) TAU1_TT=true; }
			if(abs(treeHiggsId)==225){if(treeLeg1Pt<treeLeg2Pt && treeTau2Loose3==1) TAU1_TT=true; }
			if(abs(treeHiggsId)==225){if(treeLeg1Pt>treeLeg2Pt && treeTau2Loose3==1) TAU2_TT=true; }
			if(abs(treeHiggsId)==225){if(treeLeg1Pt<treeLeg2Pt && treeTau1Loose3==1) TAU2_TT=true; }
		}
		if(tauIso==1.0){
			if((abs(treeHiggsId)==165 || abs(treeHiggsId)==195) && (abs(treeLeg1Id)<15) && treeTau2Medium3==1) TAU_LT=true;
			if((abs(treeHiggsId)==165 || abs(treeHiggsId)==195) && (abs(treeLeg2Id)<15) && treeTau1Medium3==1) TAU_LT=true;
			if(abs(treeHiggsId)==225){if(treeLeg1Pt>treeLeg2Pt && treeTau1Medium3==1) TAU1_TT=true; }
			if(abs(treeHiggsId)==225){if(treeLeg1Pt<treeLeg2Pt && treeTau2Medium3==1) TAU1_TT=true; }
			if(abs(treeHiggsId)==225){if(treeLeg1Pt>treeLeg2Pt && treeTau2Medium3==1) TAU2_TT=true; }
			if(abs(treeHiggsId)==225){if(treeLeg1Pt<treeLeg2Pt && treeTau1Medium3==1) TAU2_TT=true; }
		}
		if(debug) cout << "Higgs ID: " << treeHiggsId << endl; 
		if(abs(treeHiggsId)==143){
			if(debug) cout << "Leg1ID: " << treeLeg1Id << " Leg1Iso: " << treeLeg1Iso << " Leg1ID: " << treeLeg1LepIDloose <<
				"- Leg2Id: " << treeLeg2Id << " Leg2Iso: " << treeLeg2Iso << " Leg2ID " << treeLeg2LepIDloose << endl;
			if(abs(treeZId)==121){
				if(treeHiggsId == -143 && abs(treeLeg1Id)==11 && (treeLeg1Iso<=elIso && treeLeg1LepIDloose==1) && (treeLeg2Iso <= muIso && treeLeg2LepIDloose==1) ){
					weight*=lepEff1;
					weight*=lepEff2;
					if(debug) cout << "--> EEEM channel, CR11 control region, OS" << endl;
					h2DAc_CR11[0]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[0]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 143 && abs(treeLeg1Id)==11 && (treeLeg1Iso<=elIso && treeLeg1LepIDloose==1) && (treeLeg2Iso <= muIso && treeLeg2LepIDloose==1) ){
					weight*=lepEff1;
					weight*=lepEff2;
					if(debug) cout << "--> EEEM channel, CR11 control region, SS" << endl;
					h2DAc_CR11[8]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[8]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeZId)==169){
				if(treeHiggsId == -143 && abs(treeLeg1Id)==11 && (treeLeg1Iso<=elIso && treeLeg1LepIDloose==1) && (treeLeg2Iso <= muIso && treeLeg2LepIDloose==1) ){
					weight*=lepEff1;
					weight*=lepEff2;
					if(debug) cout << "--> MMEM channel, CR11 control region, OS" << endl;
					h2DAc_CR11[4]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[4]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 143 && abs(treeLeg1Id)==11 && (treeLeg1Iso<=elIso && treeLeg1LepIDloose==1) && (treeLeg2Iso <= muIso && treeLeg2LepIDloose==1) ){
					weight*=lepEff1;
					weight*=lepEff2;
					if(debug) cout << "--> MMEM channel, CR11 control region, SS" << endl;
					h2DAc_CR11[12]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[12]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeLeg1Id)==11 && (treeLeg1Iso > elIso || treeLeg1LepIDloose==0) && (treeLeg2Iso <= muIso && treeLeg2LepIDloose==1) ){
				if(treeHiggsId == -143){
					f1jOS=eleFR->GetBinContent(eleFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR01 control region, OS" << endl;
						h2DAc_CR01[0]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[0]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMEM channel, CR01 control region, OS" << endl;
						h2DAc_CR01[4]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[4]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
				}
				if(treeHiggsId == 143){
					f1jSS=eleFR->GetBinContent(eleFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR01 control region, SS" << endl;
						h2DAc_CR01[8]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[8]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMEM channel, CR01 control region, SS" << endl;
						h2DAc_CR01[12]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[12]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
				}
			}
			if(abs(treeLeg1Id)==11 && (treeLeg1Iso <= elIso && treeLeg1LepIDloose==1) && (treeLeg2Iso > muIso || treeLeg2LepIDloose==0)){
				if(treeHiggsId == -143){
					f2jOS=muoFR->GetBinContent(muoFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR10 control region, OS" << endl;
						h2DAc_CR10[0]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[0]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
					if(abs(treeZId)==169){ 
						if(debug) cout << "--> MMEM channel, CR10 control region, OS" << endl;
						h2DAc_CR10[4]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[4]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
				}
				if(treeHiggsId == 143){
					f2jSS=muoFR->GetBinContent(muoFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR10 control region, SS" << endl;
						h2DAc_CR10[8]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[8]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMEM channel, CR10 control region, SS" << endl;
						h2DAc_CR10[12]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[12]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
				}
			}
			if(abs(treeLeg1Id)==11 && (treeLeg1Iso > elIso || treeLeg1LepIDloose==0) && (treeLeg2Iso > muIso || treeLeg2LepIDloose==0)){
				if(treeHiggsId == -143){
					f1jOS=eleFR->GetBinContent(eleFR->FindBin(treeLeg1ClosestJetPt));
					f2jOS=muoFR->GetBinContent(muoFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){ 
						if(debug) cout << "--> EEEM channel, CR00 control region, OS" << endl;
						h2DAc_CR00[0]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[0]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMEM channel, CR00 control region, OS" << endl;
						h2DAc_CR00[4]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[4]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
				}
				if(treeHiggsId == 143){
					f1jSS=eleFR->GetBinContent(eleFR->FindBin(treeLeg1ClosestJetPt));
					f2jSS=muoFR->GetBinContent(muoFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR00 control region, SS" << endl;
						h2DAc_CR00[8]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[8]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMEM channel, CR00 control region, SS" << endl;
						h2DAc_CR00[12]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[12]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
				}
			}
			if(abs(treeZId)==121){
				if(treeHiggsId == -143 && abs(treeLeg2Id)==11 && (treeLeg2Iso<=elIso && treeLeg2LepIDloose==1) && (treeLeg1Iso <= muIso && treeLeg1LepIDloose==1) ){
					weight*=lepEff1;
					weight*=lepEff2;
					if(debug) cout << "--> EEEM channel, CR11 control region, OS" << endl;
					h2DAc_CR11[0]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[0]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 143 && abs(treeLeg2Id)==11 && (treeLeg2Iso<=elIso && treeLeg2LepIDloose==1) && (treeLeg1Iso <= muIso && treeLeg1LepIDloose==1) ){
					weight*=lepEff1;
					weight*=lepEff2;
					if(debug) cout << "--> EEEM channel, CR11 control region, SS" << endl;
					h2DAc_CR11[8]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[8]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeZId)==169){
				if(treeHiggsId == -143 && abs(treeLeg2Id)==11 && (treeLeg2Iso<=elIso && treeLeg2LepIDloose==1) && (treeLeg1Iso <= muIso && treeLeg1LepIDloose==1) ){
					weight*=lepEff1;
					weight*=lepEff2;
					if(debug) cout << "--> MMEM channel, CR11 control region, OS" << endl;
					h2DAc_CR11[4]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[4]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 143 && abs(treeLeg2Id)==11 && (treeLeg2Iso<=elIso && treeLeg2LepIDloose==1) && (treeLeg1Iso <= muIso && treeLeg1LepIDloose==1) ){
					weight*=lepEff1;
					weight*=lepEff2;
					if(debug) cout << "--> MMEM channel, CR11 control region, SS" << endl;
					h2DAc_CR11[12]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[12]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeLeg2Id)==11 && (treeLeg2Iso > elIso || treeLeg2LepIDloose==0) && (treeLeg1Iso <= muIso && treeLeg1LepIDloose==1)){
				if(treeHiggsId == -143){
					f1jOS=eleFR->GetBinContent(eleFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR01 control region, OS" << endl;
						h2DAc_CR01[0]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[0]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMEM channel, CR01 control region, OS" << endl;
						h2DAc_CR01[4]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[4]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
				}
				if(treeHiggsId ==143){
					f1jSS=eleFR->GetBinContent(eleFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR01 control region, SS" << endl;
						h2DAc_CR01[8]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[8]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMEM channel, CR01 control region, SS" << endl;
						h2DAc_CR01[12]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[12]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
				}
			}
			if(abs(treeLeg2Id)==11 && (treeLeg2Iso <= elIso && treeLeg2LepIDloose==1) && (treeLeg1Iso > muIso || treeLeg1LepIDloose==0)){
				if(treeHiggsId == -143){
					f2jOS=muoFR->GetBinContent(muoFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR10 control region, OS" << endl;
						h2DAc_CR10[0]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[0]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMEM channel, CR10 control region, OS" << endl;
						h2DAc_CR10[4]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[4]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
				}
				if(treeHiggsId ==143){
					f2jSS=muoFR->GetBinContent(muoFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR10 control region, SS" << endl;
						h2DAc_CR10[8]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[8]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMEM channel, CR10 control region, SS" << endl;
						h2DAc_CR10[12]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[12]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
				}
			}
			if(abs(treeLeg2Id)==11 && (treeLeg2Iso > elIso || treeLeg2LepIDloose==0) && (treeLeg1Iso > muIso || treeLeg1LepIDloose==0)){
				if(treeHiggsId == -143){
					f1jOS=eleFR->GetBinContent(eleFR->FindBin(treeLeg2ClosestJetPt));
					f2jOS=muoFR->GetBinContent(muoFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR00 control region, OS" << endl;
						h2DAc_CR00[0]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[0]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
					if(abs(treeZId)==169){ 
						if(debug) cout << "--> MMEM channel, CR00 control region, OS" << endl;
						h2DAc_CR00[4]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[4]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
				}
				if(treeHiggsId ==143){
					f1jSS=eleFR->GetBinContent(eleFR->FindBin(treeLeg2ClosestJetPt));
					f2jSS=muoFR->GetBinContent(muoFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEEM channel, CR00 control region, SS" << endl;
						h2DAc_CR00[8]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[8]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMEM channel, CR00 control region, SS" << endl;
						h2DAc_CR00[12]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[12]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
				}
			}
		} //******** END OF THE EMU CHANNEL
		if(abs(treeHiggsId)==165){
			if(debug){
				if(abs(treeLeg1Id)==11) cout << "Leg1ID: " << treeLeg1Id << " Leg1Iso: " << treeLeg1Iso << " Leg1ID: " << treeLeg1LepIDtight << "- TAU_LT: " << TAU_LT << endl;
				if(abs(treeLeg2Id)==11) cout << "Leg2ID: " << treeLeg2Id << " Leg2Iso: " << treeLeg2Iso << " Leg2ID: " << treeLeg2LepIDtight << "- TAU_LT: " << TAU_LT << endl;
			}
			if(abs(treeZId)==121){
				if(treeHiggsId == -165 && abs(treeLeg1Id)==11 && (treeLeg1Iso<=elIso && treeLeg1LepIDtight==1) && TAU_LT ){
					weight*=lepEff1;
					h2DAc_CR11[1]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[1]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 165 && abs(treeLeg1Id)==11 && (treeLeg1Iso<=elIso && treeLeg1LepIDtight==1) && TAU_LT){
					weight*=lepEff1;
					h2DAc_CR11[9]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[9]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeZId)==169){
				if(treeHiggsId == -165 && abs(treeLeg1Id)==11 && (treeLeg1Iso<=elIso && treeLeg1LepIDtight==1) && TAU_LT ){
					weight*=lepEff1;
					h2DAc_CR11[5]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[5]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 165 && abs(treeLeg1Id)==11 && (treeLeg1Iso<=elIso && treeLeg1LepIDtight==1) && TAU_LT){
					weight*=lepEff1;
					h2DAc_CR11[13]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[13]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeLeg1Id)==11 && (treeLeg1Iso > elIso || treeLeg1LepIDtight==0) && TAU_LT){
				if(treeHiggsId == -165){
					f1jOS=eleFR->GetBinContent(eleFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR01 control region, OS" << endl;
						h2DAc_CR01[1]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[1]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR01 control region, OS" << endl;
						h2DAc_CR01[5]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[5]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
				}
				if(treeHiggsId ==165){
					f1jSS=eleFR->GetBinContent(eleFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR01 control region, SS" << endl;
						h2DAc_CR01[9]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[9]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR01 control region, SS" << endl;
						h2DAc_CR01[13]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[13]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
				}
			}
			if(abs(treeLeg1Id)==11 && treeLeg1Iso <= elIso && treeLeg1LepIDtight==1 && !TAU_LT){
				if(treeHiggsId == -165){
					f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					//cout << "PT " << treeLeg2ClosestJetPt << " value FR (tau): " << f2jOS << endl;
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR10 control region, OS" << endl;
						h2DAc_CR10[1]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[1]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR10 control region, OS" << endl;
						h2DAc_CR10[5]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[5]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
				}
				if(treeHiggsId ==165){
					f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR10 control region, SS" << endl;
						h2DAc_CR10[9]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[9]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
					if(abs(treeZId)==169){ 
						if(debug) cout << "--> MMET channel, CR10 control region, SS" << endl;
						h2DAc_CR10[13]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[13]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
				}
			}
			if(abs(treeLeg1Id)==11 && (treeLeg1Iso > elIso || treeLeg1LepIDtight==0) && !TAU_LT){
				if(treeHiggsId == -165){
					f1jOS=eleFR->GetBinContent(eleFR->FindBin(treeLeg1ClosestJetPt));
					f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR00 control region, OS" << endl;
						h2DAc_CR00[1]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[1]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR00 control region, OS" << endl;
						h2DAc_CR00[5]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[5]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
				}
				if(treeHiggsId ==165){
					f1jSS=eleFR->GetBinContent(eleFR->FindBin(treeLeg1ClosestJetPt));
					f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR00 control region, SS" << endl;
						h2DAc_CR00[9]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[9]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR00 control region, SS" << endl;
						h2DAc_CR00[13]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[13]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
				}
			}
			if(abs(treeZId)==121){
				if(treeHiggsId == -165 && abs(treeLeg2Id)==11 && (treeLeg2Iso<=elIso && treeLeg2LepIDtight==1) && TAU_LT ){
					weight*=lepEff2;
					h2DAc_CR11[1]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[1]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 165 && abs(treeLeg2Id)==11 && (treeLeg2Iso<=elIso && treeLeg2LepIDtight==1) && TAU_LT){
					weight*=lepEff2;
					h2DAc_CR11[9]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[9]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeZId)==169){
				if(treeHiggsId == -165 && abs(treeLeg2Id)==11 && (treeLeg2Iso<=elIso && treeLeg2LepIDtight==1) && TAU_LT ){
					weight*=lepEff2;
					h2DAc_CR11[5]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[5]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 165 && abs(treeLeg2Id)==11 && (treeLeg2Iso<=elIso && treeLeg2LepIDtight==1) && TAU_LT){
					weight*=lepEff2;
					h2DAc_CR11[13]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[13]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeLeg2Id)==11 && (treeLeg2Iso > elIso || treeLeg2LepIDtight==0) && TAU_LT){
				if(treeHiggsId == -165){
					f1jOS=eleFR->GetBinContent(eleFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR01 control region, OS" << endl;
						h2DAc_CR01[1]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[1]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR01 control region, OS" << endl;
						h2DAc_CR01[5]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[5]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
				}
				if(treeHiggsId ==165){
					f1jSS=eleFR->GetBinContent(eleFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR01 control region, SS" << endl;
						h2DAc_CR01[9]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[9]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR01 control region, SS" << endl;
						h2DAc_CR01[13]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[13]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
				}
			}
			if(abs(treeLeg2Id)==11 && treeLeg2Iso <= elIso && treeLeg2LepIDtight==1 && !TAU_LT){
				if(treeHiggsId == -165){
					f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR10 control region, OS" << endl;
						h2DAc_CR10[1]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[1]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR10 control region, OS" << endl;
						h2DAc_CR10[5]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[5]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
				}
				if(treeHiggsId ==165){
					f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR10 control region, SS" << endl;
						h2DAc_CR10[9]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[9]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR10 control region, SS" << endl;
						h2DAc_CR10[13]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[13]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
				}
			}
			if(abs(treeLeg2Id)==11 && (treeLeg2Iso > elIso || treeLeg2LepIDtight==0) && !TAU_LT){
				if(treeHiggsId == -165){
					f1jOS=eleFR->GetBinContent(eleFR->FindBin(treeLeg2ClosestJetPt));
					f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR00 control region, OS" << endl;
						h2DAc_CR00[1]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[1]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR00 control region, OS" << endl;
						h2DAc_CR00[5]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[5]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
				}
				if(treeHiggsId ==165){
					f1jSS=eleFR->GetBinContent(eleFR->FindBin(treeLeg2ClosestJetPt));
					f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEET channel, CR00 control region, SS" << endl;
						h2DAc_CR00[9]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[9]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMET channel, CR00 control region, SS" << endl;
						h2DAc_CR00[13]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[13]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
				}
			}
		} //******** END OF THE ETAU CHANNEL
		if(abs(treeHiggsId)==195){
			if(debug){
				if(abs(treeLeg1Id)==13) cout << "Leg1ID: " << treeLeg1Id << " Leg1Iso: " << treeLeg1Iso << " Leg1ID: " << treeLeg1LepIDloose << "- TAU_LT: " << TAU_LT << endl;
				if(abs(treeLeg2Id)==13) cout << "Leg2ID: " << treeLeg2Id << " Leg2Iso: " << treeLeg2Iso << " Leg2ID: " << treeLeg2LepIDloose << "- TAU_LT: " << TAU_LT << endl;
			}
			if(abs(treeZId)==121){
				if(treeHiggsId == -195 && abs(treeLeg1Id)==13 && (treeLeg1Iso<=muIso && treeLeg1LepIDloose==1) && TAU_LT ){
					weight*=lepEff1;
					h2DAc_CR11[2]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[2]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 195 && abs(treeLeg1Id)==13 && (treeLeg1Iso<=muIso && treeLeg1LepIDloose==1) && TAU_LT){
					weight*=lepEff1;
					h2DAc_CR11[10]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[10]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeZId)==169){
				if(treeHiggsId == -195 && abs(treeLeg1Id)==13 && (treeLeg1Iso<=muIso && treeLeg1LepIDloose==1) && TAU_LT ){
					weight*=lepEff1;
					h2DAc_CR11[6]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[6]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 195 && abs(treeLeg1Id)==13 && (treeLeg1Iso<=muIso && treeLeg1LepIDloose==1) && TAU_LT){
					weight*=lepEff1;
					h2DAc_CR11[14]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[14]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeLeg1Id)==13 && (treeLeg1Iso > muIso || treeLeg1LepIDloose==0) && TAU_LT){
				if(treeHiggsId == -195){
					f1jOS=muoFR->GetBinContent(muoFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR01 control region, OS" << endl;
						h2DAc_CR01[2]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[2]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMMT channel, CR01 control region, OS" << endl;
						h2DAc_CR01[6]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[6]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
				}
				if(treeHiggsId ==195){
					f1jSS=muoFR->GetBinContent(muoFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR01 control region, SS" << endl;
						h2DAc_CR01[10]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[10]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMMT channel, CR01 control region, SS" << endl;
						h2DAc_CR01[14]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[14]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
				}
			}
			if(abs(treeLeg1Id)==13 && (treeLeg1Iso <= muIso && treeLeg1LepIDloose==1) && !TAU_LT){
				if(treeHiggsId == -195){
					f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR10 control region, OS" << endl;
						h2DAc_CR10[2]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[2]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
					if(abs(treeZId)==169){ 
						if(debug) cout << "--> MMMT channel, CR10 control region, OS" << endl;
						h2DAc_CR10[6]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[6]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
				}
				if(treeHiggsId ==195){
					f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR10 control region, SS" << endl;
						h2DAc_CR10[10]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[10]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMMT channel, CR10 control region, SS" << endl;
						h2DAc_CR10[14]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[14]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
				}
			}
			if(abs(treeLeg1Id)==13 && (treeLeg1Iso > muIso || treeLeg1LepIDloose==0) && !TAU_LT){
				if(treeHiggsId == -195){
					f1jOS=muoFR->GetBinContent(muoFR->FindBin(treeLeg1ClosestJetPt));
					f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR00 control region, OS" << endl;
						h2DAc_CR00[2]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[2]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMMT channel, CR00 control region, OS" << endl;
						h2DAc_CR00[6]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[6]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
				}
				if(treeHiggsId ==195){
					f1jSS=muoFR->GetBinContent(muoFR->FindBin(treeLeg1ClosestJetPt));
					f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR00 control region, SS" << endl;
						h2DAc_CR00[10]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[10]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMMT channel, CR00 control region, SS" << endl;
						h2DAc_CR00[14]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[14]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
				}
			}
			if(abs(treeZId)==121){
				if(treeHiggsId == -195 && abs(treeLeg2Id)==13 && (treeLeg2Iso<=muIso && treeLeg2LepIDloose==1) && TAU_LT ){
					weight*=lepEff2;
					h2DAc_CR11[2]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[2]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 195 && abs(treeLeg2Id)==13 && (treeLeg2Iso<=muIso && treeLeg2LepIDloose==1) && TAU_LT){
					weight*=lepEff2;
					h2DAc_CR11[10]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[10]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeZId)==169){
				if(treeHiggsId == -195 && abs(treeLeg2Id)==13 && (treeLeg2Iso<=muIso && treeLeg2LepIDloose==1) && TAU_LT ){
					weight*=lepEff2;
					h2DAc_CR11[6]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[6]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 195 && abs(treeLeg2Id)==13 && (treeLeg2Iso<=muIso && treeLeg2LepIDloose==1) && TAU_LT){
					weight*=lepEff2;
					h2DAc_CR11[14]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[14]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeLeg2Id)==13 && (treeLeg2Iso > muIso || treeLeg2LepIDloose==0) && TAU_LT){
				if(treeHiggsId == -195){
					f1jOS=muoFR->GetBinContent(muoFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR01 control region, OS" << endl;
						h2DAc_CR01[2]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[2]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMMT channel, CR01 control region, OS" << endl;
						h2DAc_CR01[6]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[6]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
				}
				if(treeHiggsId ==195){
					f1jSS=muoFR->GetBinContent(muoFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR01 control region, SS" << endl;
						h2DAc_CR01[10]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[10]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
					if(abs(treeZId)==169){ 
						if(debug) cout << "--> MMMT channel, CR01 control region, SS" << endl;
						h2DAc_CR01[14]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[14]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
				}
			}
			if(abs(treeLeg2Id)==13 && (treeLeg2Iso <= muIso && treeLeg2LepIDloose==1) && !TAU_LT){
				if(treeHiggsId == -195){
					f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR10 control region, OS" << endl;
						h2DAc_CR10[2]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[2]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMMT channel, CR10 control region, OS" << endl;
						h2DAc_CR10[6]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[6]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
				}
				if(treeHiggsId ==195){
					f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR10 control region, SS" << endl;
						h2DAc_CR10[10]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[10]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
					if(abs(treeZId)==169){ 
						if(debug) cout << "--> MMMT channel, CR10 control region, SS" << endl;
						h2DAc_CR10[14]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[14]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
				}
			}
			if(abs(treeLeg2Id)==13 && (treeLeg2Iso > muIso || treeLeg2LepIDloose==0) && !TAU_LT){
				if(treeHiggsId == -195){
					f1jOS=muoFR->GetBinContent(muoFR->FindBin(treeLeg2ClosestJetPt));
					f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR00 control region, OS" << endl;
						h2DAc_CR00[2]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[2]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMMT channel, CR00 control region, OS" << endl;
						h2DAc_CR00[6]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[6]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
				}
				if(treeHiggsId ==195){
					f1jSS=muoFR->GetBinContent(muoFR->FindBin(treeLeg2ClosestJetPt));
					f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EEMT channel, CR00 control region, SS" << endl;
						h2DAc_CR00[10]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[10]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMMT channel, CR00 control region, SS" << endl;
						h2DAc_CR00[14]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[14]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
				}
			}
		}//******** END OF THE MUTAU CHANNEL
		if(abs(treeHiggsId)==225){
			if(debug) cout << "TAU1_TT: " << TAU1_TT  << " - TAU2_TT: " << TAU2_TT << endl;
			if(abs(treeZId)==121){
				if(treeHiggsId == -225 && TAU1_TT && TAU2_TT){
					if(debug) cout << "--> EETT channel, CR11 control region, OS" << endl;
					h2DAc_CR11[3]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[3]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 225 && TAU1_TT && TAU2_TT){
					if(debug) cout << "--> EETT channel, CR11 control region, SS" << endl;
					h2DAc_CR11[11]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[11]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(abs(treeZId)==169){
				if(treeHiggsId == -225 && TAU1_TT && TAU2_TT){
					if(debug) cout << "--> MMTT channel, CR11 control region, OS" << endl;
					h2DAc_CR11[7]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[7]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
				if(treeHiggsId == 225 && TAU1_TT && TAU2_TT){
					if(debug) cout << "--> MMTT channel, CR11 control region, SS" << endl;
					h2DAc_CR11[15]->Fill(cut_index, treeAsvfMass);
					h2DA_CR11[15]->Fill(cut_index, treeAsvfMass, weight*sign);
				}
			}
			if(TAU1_TT && !TAU2_TT){
				if(treeHiggsId == -225){
					if(treeLeg1Pt>treeLeg2Pt) f1jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(treeLeg2Pt>treeLeg1Pt) f1jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EETT channel, CR01 control region, OS" << endl;
						h2DAc_CR01[3]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[3]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMTT channel, CR01 control region, OS" << endl;
						h2DAc_CR01[7]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[7]->Fill(cut_index, treeAsvfMass, (f1jOS*weight*sign)/(1-f1jOS));
					}
				}
				if(treeHiggsId ==225){
					if(treeLeg1Pt>treeLeg2Pt) f1jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(treeLeg2Pt>treeLeg1Pt) f1jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EETT channel, CR01 control region, SS" << endl;
						h2DAc_CR01[11]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[11]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMTT channel, CR01 control region, SS" << endl;
						h2DAc_CR01[15]->Fill(cut_index, treeAsvfMass);
						h2DA_CR01[15]->Fill(cut_index, treeAsvfMass, (f1jSS*weight*sign)/(1-f1jSS));
					}
				}
			}
			else if(!TAU1_TT && TAU2_TT){
				if(treeHiggsId == -225){
					if(treeLeg1Pt>treeLeg2Pt) f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(treeLeg2Pt>treeLeg1Pt) f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EETT channel, CR10 control region, OS" << endl;
						h2DAc_CR10[3]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[3]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMTT channel, CR10 control region, OS" << endl;
						h2DAc_CR10[7]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[7]->Fill(cut_index, treeAsvfMass, (f2jOS*weight*sign)/(1-f2jOS));
					}
				}
				if(treeHiggsId == 225){
					if(treeLeg1Pt>treeLeg2Pt) f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(treeLeg2Pt>treeLeg1Pt) f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EETT channel, CR10 control region, SS" << endl;
						h2DAc_CR10[11]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[11]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
					if(abs(treeZId)==169){ 
						if(debug) cout << "--> MMTT channel, CR10 control region, SS" << endl;
						h2DAc_CR10[15]->Fill(cut_index, treeAsvfMass);
						h2DA_CR10[15]->Fill(cut_index, treeAsvfMass, (f2jSS*weight*sign)/(1-f2jSS));
					}
				}
			}
			else if(!TAU1_TT && !TAU2_TT){
				if(treeHiggsId == -225){
					f1jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					f2jOS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EETT channel, CR00 control region, OS" << endl;
						h2DAc_CR00[3]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[3]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
					if(abs(treeZId)==169){ 
						if(debug) cout << "--> MMTT channel, CR00 control region, OS" << endl;
						h2DAc_CR00[7]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[7]->Fill(cut_index, treeAsvfMass, (f1jOS*f2jOS*weight*sign)/((1-f1jOS)*(1-f2jOS)));
					}
				}
				if(treeHiggsId == 225){
					f1jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg2ClosestJetPt));
					f2jSS=tauFR->GetBinContent(tauFR->FindBin(treeLeg1ClosestJetPt));
					if(abs(treeZId)==121){
						if(debug) cout << "--> EETT channel, CR00 control region, SS" << endl;
						h2DAc_CR00[11]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[11]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
					if(abs(treeZId)==169){
						if(debug) cout << "--> MMTT channel, CR00 control region, SS" << endl;
						h2DAc_CR00[15]->Fill(cut_index, treeAsvfMass);
						h2DA_CR00[15]->Fill(cut_index, treeAsvfMass, (f1jSS*f2jSS*weight*sign)/((1-f1jSS)*(1-f2jSS)));
					}
				}
			}
		}//******** END OF THE TAUTAU CHANNEL
	}
}

vector<string> FRAnalyzer::GetStringVector(string s){
	vector<string> svec;
	svec.push_back(("ee_OSelmu"+s).c_str());
	svec.push_back(("ee_OSelha"+s).c_str());
	svec.push_back(("ee_OSmuha"+s).c_str());
	svec.push_back(("ee_OShaha"+s).c_str());
	svec.push_back(("mumu_OSelmu"+s).c_str());
	svec.push_back(("mumu_OSelha"+s).c_str());
	svec.push_back(("mumu_OSmuha"+s).c_str());
	svec.push_back(("mumu_OShaha"+s).c_str());
	svec.push_back(("ee_SSelmu"+s).c_str());
	svec.push_back(("ee_SSelha"+s).c_str());
	svec.push_back(("ee_SSmuha"+s).c_str());
	svec.push_back(("ee_SShaha"+s).c_str());
	svec.push_back(("mumu_SSelmu"+s).c_str());
	svec.push_back(("mumu_SSelha"+s).c_str());
	svec.push_back(("mumu_SSmuha"+s).c_str());
	svec.push_back(("mumu_SShaha"+s).c_str());

	return svec;
}


vector<TFile*> FRAnalyzer::Overall(){

	vector<TH1F*> ratios = GetRatios(false);
	vector<TFile*> finalFile;

	vector<TH2F*> h2DAc_CR01  = BuildShapeHistograms("H","_CR01"); vector<TH2F*> h2DAc_CR01_bkg  = BuildShapeHistograms("H","_MC_CR01"); 
	vector<TH2F*> h2DAc_CR10  = BuildShapeHistograms("H","_CR10"); vector<TH2F*> h2DAc_CR10_bkg  = BuildShapeHistograms("H","_MC_CR10"); 
	vector<TH2F*> h2DAc_CR00  = BuildShapeHistograms("H","_CR00"); vector<TH2F*> h2DAc_CR00_bkg  = BuildShapeHistograms("H","_MC_CR00"); 
	vector<TH2F*> h2DAc_CR11  = BuildShapeHistograms("H","_CR11"); vector<TH2F*> h2DAc_CR11_bkg  = BuildShapeHistograms("H","_MC_CR11"); 

	vector<TH2F*> h2DA_CR01  = BuildShapeHistograms("A","_CR01"); vector<TH2F*> h2DA_CR01_bkg  = BuildShapeHistograms("A","_MC_CR01"); 
	vector<TH2F*> h2DA_CR10  = BuildShapeHistograms("A","_CR10"); vector<TH2F*> h2DA_CR10_bkg  = BuildShapeHistograms("A","_MC_CR10");
	vector<TH2F*> h2DA_CR00  = BuildShapeHistograms("A","_CR00"); vector<TH2F*> h2DA_CR00_bkg  = BuildShapeHistograms("A","_MC_CR00");
	vector<TH2F*> h2DA_CR11  = BuildShapeHistograms("A","_CR11"); vector<TH2F*> h2DA_CR11_bkg  = BuildShapeHistograms("A","_MC_CR11");

	tauIso.push_back(0.0);
	eleIso.push_back(0.3);
	muoIso.push_back(0.3);

	int cut_index=0;

	if(m_Data){
		for(double elIso=0;elIso<eleIso.size();elIso++){
			for(double muIso=0;muIso<muoIso.size();muIso++){
				for(double taIso=0.;taIso<tauIso.size();taIso++){
					for(double sumPt=0;sumPt<=200;sumPt+=20){
						Estimate(sample,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01,h2DAc_CR10,h2DAc_CR00,h2DAc_CR11,
								h2DA_CR01,h2DA_CR10,h2DA_CR00,h2DA_CR11); //main
						Estimate(dir_zz2l2tau,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01_bkg,h2DAc_CR10_bkg,h2DAc_CR00_bkg,h2DAc_CR11_bkg,
								h2DA_CR01_bkg,h2DA_CR10_bkg,h2DA_CR00_bkg,h2DA_CR11_bkg); //MC
						Estimate(dir_zz4l,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01_bkg,h2DAc_CR10_bkg,h2DAc_CR00_bkg,h2DAc_CR11_bkg,
								h2DA_CR01_bkg,h2DA_CR10_bkg,h2DA_CR00_bkg,h2DA_CR11_bkg); //MC
						Estimate(dir_zh,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01_bkg,h2DAc_CR10_bkg,h2DAc_CR00_bkg,h2DAc_CR11_bkg,
								h2DA_CR01_bkg,h2DA_CR10_bkg,h2DA_CR00_bkg,h2DA_CR11_bkg); //MC
						Estimate(dir_wwz,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01_bkg,h2DAc_CR10_bkg,h2DAc_CR00_bkg,h2DAc_CR11_bkg,
								h2DA_CR01_bkg,h2DA_CR10_bkg,h2DA_CR00_bkg,h2DA_CR11_bkg); //MC
						Estimate(dir_wzz,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01_bkg,h2DAc_CR10_bkg,h2DAc_CR00_bkg,h2DAc_CR11_bkg,
								h2DA_CR01_bkg,h2DA_CR10_bkg,h2DA_CR00_bkg,h2DA_CR11_bkg); //MC
						Estimate(dir_ttZ,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01_bkg,h2DAc_CR10_bkg,h2DAc_CR00_bkg,h2DAc_CR11_bkg,
								h2DA_CR01_bkg,h2DA_CR10_bkg,h2DA_CR00_bkg,h2DA_CR11_bkg); //MC
						Estimate(dir_zzz,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01_bkg,h2DAc_CR10_bkg,h2DAc_CR00_bkg,h2DAc_CR11_bkg,
								h2DA_CR01_bkg,h2DA_CR10_bkg,h2DA_CR00_bkg,h2DA_CR11_bkg); //MC
						cut_index++;
					}
				}
			}
		}

		TFile *fileMain  = new TFile((m_Dir+"/frMain.root").c_str(), "RECREATE");
		for(uint i=0;i<h2DAc_CR00.size();i++){
			h2DAc_CR01[i]->Write((finalStates.at(i)+"_counts_CR01").c_str());
			h2DAc_CR10[i]->Write((finalStates.at(i)+"_counts_CR10").c_str());
			h2DAc_CR00[i]->Write((finalStates.at(i)+"_counts_CR00").c_str());
			h2DAc_CR11[i]->Write((finalStates.at(i)+"_counts_CR11").c_str());
			h2DA_CR01[i]->Write((finalStates.at(i)+"_Asvfit_shapes_CR01").c_str());
			h2DA_CR10[i]->Write((finalStates.at(i)+"_Asvfit_shapes_CR10").c_str());
			h2DA_CR00[i]->Write((finalStates.at(i)+"_Asvfit_shapes_CR00").c_str());
			h2DA_CR11[i]->Write((finalStates.at(i)+"_Asvfit_shapes_CR11").c_str());
		}
		TFile *fileMC  = new TFile((m_Dir+"/frMC.root").c_str(), "RECREATE");
		for(uint i=0;i<h2DAc_CR00_bkg.size();i++){
			h2DAc_CR01_bkg[i]->Write((finalStates.at(i)+"_counts_MC_CR01").c_str());
			h2DAc_CR10_bkg[i]->Write((finalStates.at(i)+"_counts_MC_CR10").c_str());
			h2DAc_CR00_bkg[i]->Write((finalStates.at(i)+"_counts_MC_CR00").c_str());
			h2DAc_CR11_bkg[i]->Write((finalStates.at(i)+"_counts_MC_CR11").c_str());
			h2DA_CR01_bkg[i]->Write((finalStates.at(i)+"_Asvfit_shapes_MC_CR01").c_str());
			h2DA_CR10_bkg[i]->Write((finalStates.at(i)+"_Asvfit_shapes_MC_CR10").c_str());
			h2DA_CR00_bkg[i]->Write((finalStates.at(i)+"_Asvfit_shapes_MC_CR00").c_str());
			h2DA_CR11_bkg[i]->Write((finalStates.at(i)+"_Asvfit_shapes_MC_CR11").c_str());
		}
		finalFile.push_back(fileMain);
		finalFile.push_back(fileMC);
	}else{
		for(double elIso=0;elIso<eleIso.size();elIso++){
			for(double muIso=0;muIso<muoIso.size();muIso++){
				for(double taIso=0.;taIso<tauIso.size();taIso++){
					for(double sumPt=0;sumPt<=200;sumPt+=20){
						Estimate(sample,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01,h2DAc_CR10,h2DAc_CR00,h2DAc_CR11,
								h2DA_CR01,h2DA_CR10,h2DA_CR00,h2DA_CR11); //main
						Estimate(dir_wz,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01,h2DAc_CR10,h2DAc_CR00,h2DAc_CR11,
								h2DA_CR01,h2DA_CR10,h2DA_CR00,h2DA_CR11); //main
						Estimate(dir_ww,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01,h2DAc_CR10,h2DAc_CR00,h2DAc_CR11,
								h2DA_CR01,h2DA_CR10,h2DA_CR00,h2DA_CR11); //main
						Estimate(dir_vg,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01,h2DAc_CR10,h2DAc_CR00,h2DAc_CR11,
								h2DA_CR01,h2DA_CR10,h2DA_CR00,h2DA_CR11); //main
						Estimate(dir_top,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01,h2DAc_CR10,h2DAc_CR00,h2DAc_CR11,
								h2DA_CR01,h2DA_CR10,h2DA_CR00,h2DA_CR11); //main
						Estimate(dir_qcd,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01,h2DAc_CR10,h2DAc_CR00,h2DAc_CR11,
								h2DA_CR01,h2DA_CR10,h2DA_CR00,h2DA_CR11); //main
						Estimate(dir_qcdmu,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01,h2DAc_CR10,h2DAc_CR00,h2DAc_CR11,
								h2DA_CR01,h2DA_CR10,h2DA_CR00,h2DA_CR11); //main
						Estimate(dir_wjets,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01,h2DAc_CR10,h2DAc_CR00,h2DAc_CR11,
								h2DA_CR01,h2DA_CR10,h2DA_CR00,h2DA_CR11); //main
						Estimate(dir_www,ratios,eleIso[elIso],muoIso[muIso],tauIso[taIso],sumPt,cut_index,
								h2DAc_CR01,h2DAc_CR10,h2DAc_CR00,h2DAc_CR11,
								h2DA_CR01,h2DA_CR10,h2DA_CR00,h2DA_CR11); //main
						cut_index++;
					}
				}
			}
		}

		TFile *fileMain  = new TFile((m_Dir+"/frMain.root").c_str(), "RECREATE");
		for(uint i=0;i<h2DAc_CR00.size();i++){
			h2DAc_CR01[i]->Write((finalStates.at(i)+"_counts_CR01").c_str());
			h2DAc_CR10[i]->Write((finalStates.at(i)+"_counts_CR10").c_str());
			h2DAc_CR00[i]->Write((finalStates.at(i)+"_counts_CR00").c_str());
			h2DAc_CR11[i]->Write((finalStates.at(i)+"_counts_CR11").c_str());
			h2DA_CR01[i]->Write((finalStates.at(i)+"_Asvfit_shapes_CR01").c_str());
			h2DA_CR10[i]->Write((finalStates.at(i)+"_Asvfit_shapes_CR10").c_str());
			h2DA_CR00[i]->Write((finalStates.at(i)+"_Asvfit_shapes_CR00").c_str());
			h2DA_CR11[i]->Write((finalStates.at(i)+"_Asvfit_shapes_CR11").c_str());
		}
		TFile *fileMC  = new TFile((m_Dir+"/frMC.root").c_str(), "RECREATE");
		for(uint i=0;i<h2DAc_CR00_bkg.size();i++){
			h2DAc_CR01_bkg[i]->Write((finalStates.at(i)+"_counts_MC_CR01").c_str());
			h2DAc_CR10_bkg[i]->Write((finalStates.at(i)+"_counts_MC_CR10").c_str());
			h2DAc_CR00_bkg[i]->Write((finalStates.at(i)+"_counts_MC_CR00").c_str());
			h2DAc_CR11_bkg[i]->Write((finalStates.at(i)+"_counts_MC_CR11").c_str());
			h2DA_CR01_bkg[i]->Write((finalStates.at(i)+"_Asvfit_shapes_MC_CR01").c_str());
			h2DA_CR10_bkg[i]->Write((finalStates.at(i)+"_Asvfit_shapes_MC_CR10").c_str());
			h2DA_CR00_bkg[i]->Write((finalStates.at(i)+"_Asvfit_shapes_MC_CR00").c_str());
			h2DA_CR11_bkg[i]->Write((finalStates.at(i)+"_Asvfit_shapes_MC_CR11").c_str());
		}
		finalFile.push_back(fileMain);
		finalFile.push_back(fileMC);
	}
	return finalFile;
}

void FRAnalyzer::SetExtremes(vector<TH2F*> hCR){
	if(debug) cout << "SETTING underflow  and overflow bins..." << endl;
	for(int i=0;i<(int)(finalStates.size()/2);i++){
		double int_OS(0.0),int_SS(0.0);
		double err_OS(0.0),err_SS(0.0);
		for(int k=0;k<(int)(hCR[0]->GetYaxis()->GetNbins());k++){
			int_OS += (hCR[i]->GetBinContent(m_Index+1,k+1));
			err_OS += TMath::Power(hCR[i]->GetBinError(m_Index+1,k+1),2);
			int_SS += (hCR[i+7]->GetBinContent(m_Index+1,k+1));
			//int_SS += (hCR[i+8]->GetBinContent(m_Index+1,k+1));
			err_SS += TMath::Power(hCR[i+7]->GetBinError(m_Index+1,k+1),2);
			//err_SS += TMath::Power(hCR[i+8]->GetBinError(m_Index+1,k+1),2);
		}
		//set in un the underflow the int +- its error for OS
		if(debug) cout << "Idx "<<m_Index<<" FS "<<i<< " INTEGRALS/ERR OS "  << int_OS <<"/"<< TMath::Sqrt(err_OS)<< endl;
		if(debug) cout << "Idx "<<m_Index<<" FS "<<i<< " INTEGRALS/ERR SS "  << int_SS <<"/" <<TMath::Sqrt(err_SS)<< endl;
		hCR[i]->SetBinContent(m_Index+1,0,int_OS); hCR[i]->SetBinError(m_Index+1,0,TMath::Sqrt(err_OS));
		hCR[i]->SetBinContent(m_Index+1,8,int_SS); hCR[i]->SetBinError(m_Index+1,8,TMath::Sqrt(err_SS));
		//hCR[i]->SetBinContent(m_Index+1,9,int_SS); hCR[i]->SetBinError(m_Index+1,9,TMath::Sqrt(err_SS));
		//set in un the underflow 0 +- 0 for SS
		hCR[i+(finalStates.size()/2)]->SetBinContent(m_Index+1,0,0); hCR[i+(finalStates.size()/2)]->SetBinError(m_Index+1,0,0);
		hCR[i+(finalStates.size()/2)]->SetBinContent(m_Index+1,8,0); hCR[i+(finalStates.size()/2)]->SetBinError(m_Index+1,8,0);
		//hCR[i+(finalStates.size()/2)]->SetBinContent(m_Index+1,9,0); hCR[i+(finalStates.size()/2)]->SetBinError(m_Index+1,9,0);
		if(debug) cout << "Idx "<<m_Index<<" FS "<<i<< " uf. OS "  << hCR[i]->GetBinContent(m_Index+1,0) << "/pm"<< hCR[i]->GetBinError(m_Index+1,0) << endl;
		if(debug) cout << "Idx "<<m_Index<<" FS "<<i<< " uf. SS "  << hCR[i+7]->GetBinContent(m_Index+1,0) <<"/pm"<<hCR[i+7]->GetBinError(m_Index+1,0) << endl;
		//if(debug) cout << "Idx "<<m_Index<<" FS "<<i<< " uf. SS "  << hCR[i+8]->GetBinContent(m_Index+1,0) <<"/pm"<<hCR[i+8]->GetBinError(m_Index+1,0) << endl;
		if(debug) cout << "Idx "<<m_Index<<" FS "<<i<< " of. OS "  << hCR[i]->GetBinContent(m_Index+1,8) << "/pm" <<hCR[i]->GetBinError(m_Index+1,8) << endl;
		//if(debug) cout << "Idx "<<m_Index<<" FS "<<i<< " of. OS "  << hCR[i]->GetBinContent(m_Index+1,9) << "/pm" <<hCR[i]->GetBinError(m_Index+1,9) << endl;
		if(debug) cout << "Idx "<<m_Index<<" FS "<<i<< " of. SS "  << hCR[i+7]->GetBinContent(m_Index+1,7) <<"/pm" <<hCR[i+7]->GetBinError(m_Index+1,7) <<endl;
		//if(debug) cout << "Idx "<<m_Index<<" FS "<<i<< " of. SS "  << hCR[i+8]->GetBinContent(m_Index+1,9) <<"/pm" <<hCR[i+8]->GetBinError(m_Index+1,9) <<endl;
	}
}

vector<TH2F*> FRAnalyzer::SumUpOSandSS(vector<TH2F*> hCR){
	vector<TH2F*> htCR;
	for(int i=0;i<(int)(finalStates.size()/2);i++){
		TH2F * tmp = (TH2F*)hCR[i]->Clone();
		tmp->Add(hCR[i+(finalStates.size()/2)]);
		htCR.push_back(tmp);	
	}
	if(debug){
		cout << "New vector size "<< htCR.size() <<endl;
		for(int j=0;j<nbx;j++){
			for(int i=0;i<(int)(finalStates.size()/2);i++){
				for(int k=0;k<=(int)(hCR[0]->GetYaxis()->GetNbins());k++){
					cout << "Idx "<<j<<" FS "<<i<<" bin " <<k<< " OS/SS/NEW:" <<
						hCR[i]->GetBinContent(j+1,k)<<"pm"<< hCR[i]->GetBinError(j+1,k) <<"/" <<
						hCR[i+7]->GetBinContent(j+1,k)<<"pm" << hCR[i+7]->GetBinError(j+1,k) <<"/" <<
						//hCR[i+8]->GetBinContent(j+1,k)<<"pm" << hCR[i+8]->GetBinError(j+1,k) <<"/" <<
						htCR[i]->GetBinContent(j+1,k)<<"pm"<<htCR[i]->GetBinError(j+1,k)<<endl;
				}
			}
		}
	}
	return htCR;
}

vector<double> FRAnalyzer::GetFRWeights(){
	vector<TH1F*> ratios = GetRatios(true);
	TH1F *eleFR=NULL,  *muoFR=NULL, *tauFR=NULL; 
	eleFR=(TH1F*)ratios[2]->Clone();
	muoFR=(TH1F*)ratios[5]->Clone();
	tauFR=(TH1F*)ratios[6]->Clone();
	double we = eleFR->GetBinContent(1)/(1-eleFR->GetBinContent(1));
	double wm = muoFR->GetBinContent(1)/(1-muoFR->GetBinContent(1));
	double wt = tauFR->GetBinContent(1)/(1-tauFR->GetBinContent(1));
	vector<double> weights;
	weights.push_back(we); weights.push_back(wm); weights.push_back(wt);
	return weights;
}

void FRAnalyzer::SetErrorForEmptyBins(vector<TH2F*> hCR01, vector<TH2F*> hCR10, vector<TH2F*> hCR00){
	double we = (GetFRWeights())[0];
	double wm = (GetFRWeights())[1];
	double wt = (GetFRWeights())[2];
	cout << "WE " << we << " WM " << wm << " WT " << wt << endl; 
	for(int j=0;j<nbx;j++){
		for(int i=0;i<(int)(finalStates.size()/2);i++){
			for(int k=0;k<(int)(hCR01[0]->GetYaxis()->GetNbins())+2;k++){
				if(i==0||i==4){
					if(hCR01[i]->GetBinContent(j+1,k)==0) hCR01[i]->SetBinError(j+1,k,1.8*we); 
					if(hCR10[i]->GetBinContent(j+1,k)==0) hCR10[i]->SetBinError(j+1,k,1.8*wm);
					if(hCR00[i]->GetBinContent(j+1,k)==0) hCR00[i]->SetBinError(j+1,k,1.8*we*wm);
				}
				if(i==1||i==5){
					if(hCR01[i]->GetBinContent(j+1,k)==0) hCR01[i]->SetBinError(j+1,k,1.8*we); 
					if(hCR10[i]->GetBinContent(j+1,k)==0) hCR10[i]->SetBinError(j+1,k,1.8*wt);
					if(hCR00[i]->GetBinContent(j+1,k)==0) hCR00[i]->SetBinError(j+1,k,1.8*we*wt);
				}
				if(i==2||i==6){
					if(hCR01[i]->GetBinContent(j+1,k)==0) hCR01[i]->SetBinError(j+1,k,1.8*wm); 
					if(hCR10[i]->GetBinContent(j+1,k)==0) hCR10[i]->SetBinError(j+1,k,1.8*wt);
					if(hCR00[i]->GetBinContent(j+1,k)==0) hCR00[i]->SetBinError(j+1,k,1.8*wm*wt);
				}
				if(i==3||i==7){
					if(hCR01[i]->GetBinContent(j+1,k)==0) hCR01[i]->SetBinError(j+1,k,1.8*wt); 
					if(hCR10[i]->GetBinContent(j+1,k)==0) hCR10[i]->SetBinError(j+1,k,1.8*wt);
					if(hCR00[i]->GetBinContent(j+1,k)==0) hCR00[i]->SetBinError(j+1,k,1.8*wt*wt);
				}
				if(debug) cout << "AFTER Idx "<<j<<" FS "<<i<<" bin " <<k<< " 01/10/00:" <<
					hCR01[i]->GetBinError(j+1,k)<<"/"<<hCR10[i]->GetBinError(j+1,k)<<"/"<<hCR00[i]->GetBinError(j+1,k)<<endl;
			}
		}
	}
}
vector<TH2F*> FRAnalyzer::GetShapes(vector<TH2F*> hCR01, vector<TH2F*> hCR10, vector<TH2F*> hCR00){
	vector<TH2F*> hSh;
	for(int i=0;i<(int)(finalStates.size());i++){
		TH2F * tmp = (TH2F*)hCR01[i]->Clone();
		tmp->Add(hCR10[i]);
		tmp->Add(hCR00[i],-1);
		hSh.push_back(tmp);
	}
	if(debug){
		for(int j=0;j<nbx;j++){
			for(int i=0;i<(int)(finalStates.size()/2);i++){
				for(int k=0;k<(int)(hCR01[0]->GetYaxis()->GetNbins())+2;k++){
					cout << "AFTER Idx "<<j<<" FS "<<i<<" bin " <<k<< " 01+10-00:" <<
						hSh[i]->GetBinContent(j+1,k)<<" pm "<<hSh[i]->GetBinError(j+1,k)<<endl;
				}
			}
		}}
	return hSh;
}

vector<TH2F*> FRAnalyzer::SubtractShapes(vector<TH2F*> hData, vector<TH2F*> hMC){
	vector<TH2F*> hSh;
	for(int i=0;i<(int)(finalStates.size()/2);i++){
		TH2F * tmp = (TH2F*)hData[i]->Clone();
		tmp->Add(hMC[i],-1);
		hSh.push_back(tmp);
	}
	if(debug){
		for(int j=0;j<nbx;j++){
			for(int i=0;i<(int)(finalStates.size()/2);i++){
				for(int k=0;k<(int)(hData[0]->GetYaxis()->GetNbins())+2;k++){
					cout << "AFTER Idx "<<j<<" FS "<<i<<" bin " <<k<< " 01+10-00:" <<
						hSh[i]->GetBinContent(j+1,k)<<" pm "<<hSh[i]->GetBinError(j+1,k)<<endl;
				}
			}
		}}
	return hSh;
}

void FRAnalyzer::RescaleShape( vector<TH2F*> shape ){	
	string shapeName = shape[0]->GetName();
	for(int i=0;i<(int)(finalStates.size()/2);i++){
		double intOS = shape[i]->GetBinContent(m_Index+1,0);
		double intSS = shape[i]->GetBinContent(m_Index+1,9);
		double errOS = shape[i]->GetBinError(m_Index+1,0);
		double SF = intOS/(intOS+intSS);
		if(debug) cout << "intOS/intSS/SF: " << intOS<<"/"<<intSS<<"/"<<SF<<endl;
		for(int k=0;k<(shape[0]->GetYaxis()->GetNbins());k++){
			shape[i]->SetBinContent(m_Index+1,k+1,SF*shape[i]->GetBinContent(m_Index+1,k+1));
			shape[i]->SetBinError(m_Index+1,k+1,SF*shape[i]->GetBinError(m_Index+1,k+1));
		}
		shape[i]->SetBinContent(m_Index+1,0,intOS);
		shape[i]->SetBinError(m_Index+1,0,errOS);
		shape[i]->SetBinContent(m_Index+1,9,0);
		shape[i]->SetBinError(m_Index+1,9,0);
	}
}
//WriteFinalNumbers(shA_toyN,shA_toyNFC,shA_toyMu,shA_toyMuFC);
void FRAnalyzer::WriteFinalNumbers(vector<TH2F*> h1, vector<TH2F*> h2, vector<TH2F*> h3, vector<TH2F*> h4, vector<TH2F*> h5, vector<TH2F*> h6){
	ofstream file;
	file.open((m_Dir+"/FinalNumbers.txt").c_str());
	//for(int j=0;j<nbx;j++){
	for(int i=0;i<(int)(finalStates.size()/2);i++){
		for(int k=0;k<(int)(h1[0]->GetYaxis()->GetNbins())+1;k++){
			file << finalStates.at(i) << "," << cutIndexLabel.at(m_Index) << " bin " << k <<
				//At_CR01,At_CR10,At_CR00
				" //toy N: " << h1[i]->GetBinContent(m_Index+1,k) << " pm "<< h1[i]->GetBinError(m_Index+1,k) << 
				" //toy N FC: " << h2[i]->GetBinContent(m_Index+1,k) << " pm "<< h2[i]->GetBinError(m_Index+1,k) << 
				//" //toy N tr. Neg.: " << h3[i]->GetBinContent(m_Index+1,k) << " pm "<< h3[i]->GetBinError(m_Index+1,k) << 
				" //toy Mu: " << h4[i]->GetBinContent(m_Index+1,k) << " pm "<< h4[i]->GetBinError(m_Index+1,k) << 
				" //toy Mu FC: " << h5[i]->GetBinContent(m_Index+1,k) << " pm "<< h5[i]->GetBinError(m_Index+1,k) << "\n";
			//" //toy Mu tr. Neg.: " << h6[i]->GetBinContent(j+1,k) << " pm "<< h6[i]->GetBinError(j+1,k) << "\n";
		}
	}
	//}

	file.close();
}

void FRAnalyzer::GetFCIntervals(vector<TH2F*> shape){
	for(int i=0;i<(int)shape.size();i++){
		if(debug) cout << " SHAPE N : " << i+1 << endl; 
		if(debug) cout << " N bins : " << shape[0]->GetYaxis()->GetNbins() << endl;
		for(int k=0;k<(shape[0]->GetYaxis()->GetNbins())+1;k++){
			if(shape[i]->GetBinContent(m_Index+1,k)>=0) continue;
			double content = shape[i]->GetBinContent(m_Index+1,k);
			double error   = shape[i]->GetBinError(m_Index+1,k);
			shape[i]->SetBinContent(m_Index+1,k,0);
			if((TMath::Abs((content/error))-0.0)<0.05) {
				shape[i]->SetBinError(m_Index+1,k,error);
				if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.9,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
				if((TMath::Abs((content/error))-0.0)>0.05&&(TMath::Abs((content/error))-0.1)<0.05) {
					shape[i]->SetBinError(m_Index+1,k,0.9*error);
					if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.9,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
					if((TMath::Abs((content/error))-0.1)>0.05&&(TMath::Abs((content/error))-0.2)<0.05) {
						shape[i]->SetBinError(m_Index+1,k,0.81*error);
						if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.81,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
						if((TMath::Abs((content/error))-0.2)>0.05&&(TMath::Abs((content/error))-0.3)<0.05) {
							shape[i]->SetBinError(m_Index+1,k,0.72*error);
							if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.72,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
							if((TMath::Abs((content/error))-0.3)>0.05&&(TMath::Abs((content/error))-0.4)<0.05) {
								shape[i]->SetBinError(m_Index+1,k,0.64*error);
								if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.64,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
								if((TMath::Abs((content/error))-0.4)>0.05&&(TMath::Abs((content/error))-0.5)<0.05) {
									shape[i]->SetBinError(m_Index+1,k,0.56*error);
									if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.56,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
									if((TMath::Abs((content/error))-0.5)>0.05&&(TMath::Abs((content/error))-0.6)<0.05) {
										shape[i]->SetBinError(m_Index+1,k,0.49*error);
										if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.49,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
										if((TMath::Abs((content/error))-0.6)>0.05&&(TMath::Abs((content/error))-0.7)<0.05) {
											shape[i]->SetBinError(m_Index+1,k,0.43*error);
											if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.43,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
											if((TMath::Abs((content/error))-0.7)>0.05&&(TMath::Abs((content/error))-0.8)<0.05) {
												shape[i]->SetBinError(m_Index+1,k,0.37*error);
												if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.37,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
												if((TMath::Abs((content/error))-0.8)>0.05&&(TMath::Abs((content/error))-0.9)<0.05) {
													shape[i]->SetBinError(m_Index+1,k,0.32*error);
													if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.32,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
													if((TMath::Abs((content/error))-0.9)>0.05&&(TMath::Abs((content/error))-1.0)<0.05) {
														shape[i]->SetBinError(m_Index+1,k,0.27*error);
														if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.27,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
														if((TMath::Abs((content/error))-1.0)>0.05&&(TMath::Abs((content/error))-1.1)<0.05) {
															shape[i]->SetBinError(m_Index+1,k,0.23*error);
															if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.23,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
															if((TMath::Abs((content/error))-1.1)>0.05&&(TMath::Abs((content/error))-1.2)<0.05) {
																shape[i]->SetBinError(m_Index+1,k,0.20*error);
																if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.20,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																if((TMath::Abs((content/error))-1.2)>0.05&&(TMath::Abs((content/error))-1.3)<0.05) {
																	shape[i]->SetBinError(m_Index+1,k,0.17*error);
																	if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.17,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																	if((TMath::Abs((content/error))-1.3)>0.05&&(TMath::Abs((content/error))-1.4)<0.05) {
																		shape[i]->SetBinError(m_Index+1,k,0.15*error);
																		if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.15,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																		if((TMath::Abs((content/error))-1.4)>0.05&&(TMath::Abs((content/error))-1.5)<0.05) {
																			shape[i]->SetBinError(m_Index+1,k,0.13*error);
																			if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.13,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																			if((TMath::Abs((content/error))-1.5)>0.05&&(TMath::Abs((content/error))-1.6)<0.05) {
																				shape[i]->SetBinError(m_Index+1,k,0.11*error);
																				if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.11,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																				if((TMath::Abs((content/error))-1.6)>0.05&&(TMath::Abs((content/error))-1.7)<0.05) {
																					shape[i]->SetBinError(m_Index+1,k,0.10*error);
																					if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.10,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																					if((TMath::Abs((content/error))-1.7)>0.05&&(TMath::Abs((content/error))-1.8)<0.05) {
																						shape[i]->SetBinError(m_Index+1,k,0.09*error);
																						if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.09,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																						if((TMath::Abs((content/error))-1.8)>0.05&&(TMath::Abs((content/error))-1.9)<0.05) {
																							shape[i]->SetBinError(m_Index+1,k,0.08*error);
																							if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.08,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																							if((TMath::Abs((content/error))-1.9)>0.05&&(TMath::Abs((content/error))-2.0)<0.05) {
																								shape[i]->SetBinError(m_Index+1,k,0.07*error);
																								if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.07,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																								if((TMath::Abs((content/error))-2.0)>0.05&&(TMath::Abs((content/error))-2.1)<0.05) {
																									shape[i]->SetBinError(m_Index+1,k,0.06*error);
																									if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.06,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																									if((TMath::Abs((content/error))-2.1)>0.05&&(TMath::Abs((content/error))-2.2)<0.05) {
																										shape[i]->SetBinError(m_Index+1,k,0.06*error);
																										if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.06,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																										if((TMath::Abs((content/error))-2.2)>0.05&&(TMath::Abs((content/error))-2.3)<0.05) {
																											shape[i]->SetBinError(m_Index+1,k,0.05*error);
																											if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.05,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																											if((TMath::Abs((content/error))-2.3)>0.05&&(TMath::Abs((content/error))-2.4)<0.05) {
																												shape[i]->SetBinError(m_Index+1,k,0.05*error);
																												if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.05,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																												if((TMath::Abs((content/error))-2.4)>0.05&&(TMath::Abs((content/error))-2.5)<0.05) {
																													shape[i]->SetBinError(m_Index+1,k,0.05*error);
																													if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.05,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																													if((TMath::Abs((content/error))-2.5)>0.05&&(TMath::Abs((content/error))-2.6)<0.05) {
																														shape[i]->SetBinError(m_Index+1,k,0.05*error);
																														if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.05,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
																														if((TMath::Abs((content/error))-2.6)>0.05) {
																															shape[i]->SetBinError(m_Index+1,k,0.04*error);
																															if(debug) cout << "m_Index: "<<m_Index<<",k: "<<k<<",content: "<<content<<",error: "<<error<<",correction:0.04,error corrected: "<<shape[i]->GetBinError(m_Index+1,k)<<endl;}
		}
	}
}


vector<TFile*> FRAnalyzer::FinalEstimate(){
	vector<TFile*> files = Overall();

	vector<TH2F*> H_CR01, H_CR01_MC;
	vector<TH2F*> H_CR10, H_CR10_MC;
	vector<TH2F*> H_CR00, H_CR00_MC;
	vector<TH2F*> H_CR11, H_CR11_MC;
	vector<TH2F*> A_CR01, A_CR01_MC;
	vector<TH2F*> A_CR10, A_CR10_MC;
	vector<TH2F*> A_CR00, A_CR00_MC;
	vector<TH2F*> A_CR11, A_CR11_MC;

	for(int i=0;i<(int)files.size();i++){
		for(int j=0;j<(int)finalStates.size();j++){
			TH2F* hcr01      = (TH2F*)files.at(i)->Get((finalStates.at(j)+"_counts_"+hlabel.at(i)+"CR01").c_str());
			TH2F* hcr10      = (TH2F*)files.at(i)->Get((finalStates.at(j)+"_counts_"+hlabel.at(i)+"CR10").c_str());
			TH2F* hcr00      = (TH2F*)files.at(i)->Get((finalStates.at(j)+"_counts_"+hlabel.at(i)+"CR00").c_str());
			TH2F* hcr11      = (TH2F*)files.at(i)->Get((finalStates.at(j)+"_counts_"+hlabel.at(i)+"CR11").c_str());
			TH2F* acr01      = (TH2F*)files.at(i)->Get((finalStates.at(j)+"_Asvfit_shapes_"+hlabel.at(i)+"CR01").c_str());
			TH2F* acr10      = (TH2F*)files.at(i)->Get((finalStates.at(j)+"_Asvfit_shapes_"+hlabel.at(i)+"CR10").c_str());
			TH2F* acr00      = (TH2F*)files.at(i)->Get((finalStates.at(j)+"_Asvfit_shapes_"+hlabel.at(i)+"CR00").c_str());
			TH2F* acr11      = (TH2F*)files.at(i)->Get((finalStates.at(j)+"_Asvfit_shapes_"+hlabel.at(i)+"CR11").c_str());


			if(i==0){
				H_CR01.push_back(hcr01); A_CR01.push_back(acr01);
				H_CR10.push_back(hcr10); A_CR10.push_back(acr10);
				H_CR00.push_back(hcr00); A_CR00.push_back(acr00);
				H_CR11.push_back(hcr11); A_CR11.push_back(acr11);
			}if(i==1){
				H_CR01_MC.push_back(hcr01); A_CR01_MC.push_back(acr01);
				H_CR10_MC.push_back(hcr10); A_CR10_MC.push_back(acr10); 
				H_CR00_MC.push_back(hcr00); A_CR00_MC.push_back(acr00); 
				H_CR11_MC.push_back(hcr11); A_CR11_MC.push_back(acr11); 
			}
		}

	}

	SetExtremes(H_CR01);    SetExtremes(H_CR10);    SetExtremes(H_CR00);     
	SetExtremes(A_CR01);    SetExtremes(A_CR10);    SetExtremes(A_CR00);     
	SetExtremes(H_CR01_MC); SetExtremes(H_CR10_MC); SetExtremes(H_CR00_MC);     
	SetExtremes(A_CR01_MC); SetExtremes(A_CR10_MC); SetExtremes(A_CR00_MC);     

	//#####################################################################################################################
	//Get prepared for toys ###############################################################################################
	vector<TH2F*> shBayFit  = BuildShapeHistograms("A","BayFit");
	vector<TH2F*> shBayHis  = BuildShapeHistograms("A","BayHis");
	vector<TH2F*> shFit     = BuildShapeHistograms("A","Fit");
	vector<TH2F*> shHis     = BuildShapeHistograms("A","His");
	vector<TH2F*> shA_toyMC = GetShapes(A_CR01_MC, A_CR10_MC, A_CR00_MC);   

	TFile *LNDist_Bay = GetToyExperiments(A_CR01,A_CR10,A_CR00,H_CR01,H_CR10,H_CR00,shBayFit,shBayHis,shA_toyMC,true);
	TFile *LNDist     = GetToyExperiments(A_CR01,A_CR10,A_CR00,H_CR01,H_CR10,H_CR00,shFit,shHis,shA_toyMC,false);

	TFile *shapeToyBayFit  = new TFile((m_Dir+"/shapeAtoy_DD_"+Idx+"_Bayesian_Fit.root").c_str(), "RECREATE");
	for(uint i=0;i<shBayFit.size()/2;i++){
		shBayFit[i]->Write((finalStates.at(i)+"_Asvfit_shapes").c_str());
	}
	TFile *shapeToyBayHis = new TFile((m_Dir+"/shapeAtoy_DD_"+Idx+"_Bayesian_His.root").c_str(), "RECREATE");
	for(uint i=0;i<shBayHis.size()/2;i++){
		shBayHis[i]->Write((finalStates.at(i)+"_Asvfit_shapes").c_str());
	}
	TFile *shapeToyFit  = new TFile((m_Dir+"/shapeAtoy_DD_"+Idx+"_Fit.root").c_str(), "RECREATE");
	for(uint i=0;i<shFit.size()/2;i++){
		shFit[i]->Write((finalStates.at(i)+"_Asvfit_shapes").c_str());
	}
	TFile *shapeToyHis = new TFile((m_Dir+"/shapeAtoy_DD_"+Idx+"_His.root").c_str(), "RECREATE");
	for(uint i=0;i<shHis.size()/2;i++){
		shHis[i]->Write((finalStates.at(i)+"_Asvfit_shapes").c_str());
	}

	vector<TFile*> file;
	//file.push_back(fileA); 
	file.push_back(LNDist_Bay); 
	file.push_back(LNDist); 
	file.push_back(shapeToyBayFit); 
	file.push_back(shapeToyBayHis); 
	file.push_back(shapeToyFit); 
	file.push_back(shapeToyHis); 
	return file;

}

pair<double,double> FRAnalyzer::GetClosureTestInfo(vector<TH2F*> shapeDD, vector<TH2F*> shapeMC, int ch1, int ch2){
	pair<double,double> info;

	TH1D *tmpDDee=NULL;
	TH1D *tmpDDmm=NULL;
	TH1D *tmpMCee=NULL;
	TH1D *tmpMCmm=NULL;
	TH1D *tmpDD=NULL;
	TH1D *tmpMC=NULL;

	tmpDDee   = shapeDD[ch1]->ProjectionY(TString::Format("p1ee_%s", shapeDD[ch1]->GetName()),m_Index+1,m_Index+1);
	tmpMCee   = shapeMC[ch1]->ProjectionY(TString::Format("p2ee_%s", shapeMC[ch1]->GetName()),m_Index+1,m_Index+1);
	tmpDDmm   = shapeDD[ch2]->ProjectionY(TString::Format("p1mm_%s", shapeDD[ch2]->GetName()),m_Index+1,m_Index+1);
	tmpMCmm   = shapeMC[ch2]->ProjectionY(TString::Format("p2mm_%s", shapeMC[ch2]->GetName()),m_Index+1,m_Index+1);
	tmpDD     = (TH1D*)tmpDDee->Clone(); tmpDD->Add(tmpDDmm);
	tmpMC     = (TH1D*)tmpMCee->Clone(); tmpMC->Add(tmpMCmm);

	//C&C result for DD is saved in the underflow bin with the full correct treatment of the errors
	double integralDD = tmpDD->GetBinContent(0);
	double tmpUncDD   = tmpDD->GetBinError(0);
	double integralMC = tmpMC->GetBinContent(0);
	double tmpUncMC   = tmpMC->GetBinError(0);
	double diff=TMath::Abs(integralDD-integralMC)/integralDD;
	double err=TMath::Sqrt( tmpUncDD * TMath::Power(integralMC/(integralDD*integralDD),2) + tmpUncMC * TMath::Power(1/integralDD,2)) ;
	info.first  = diff;
	info.second = err;

	return info;
}

TGraphErrors* FRAnalyzer::GetClosureTestGraph(pair<double,double> info){

        double x[1]; double ex[1];
        double y[1]; double ey[1];
	x[0]  = m_Index+1;
	ex[0] = 0;
	y[0]  = info.first;
	ey[0] = info.second;

        TGraphErrors* graph = new TGraphErrors(1, x, y, ex, ey);
        return graph;
}

TFile* FRAnalyzer::ClosureTest(){
	vector<TFile*> files=FinalEstimate();
	vector<TH2F*> shapeDD;

	for(int i=0;i<(int)finalStates.size()/2;i++){
		TH2F* dd  = (TH2F*)files[2]->Get((finalStates.at(i)+"_Asvfit_shapes").c_str());
		shapeDD.push_back(dd);
	}
 
	vector<TH2F*> shapeMC  = GetStraightMC();
	SetExtremes(shapeMC);
	vector<TH2F*> shapeMC2 = SumUpOSandSS(shapeMC);
	RescaleShape(shapeMC2);

        cout << "MC SHAPE SIZE " << shapeMC2.size() << endl;
	pair<double,double> emInfo   = GetClosureTestInfo(shapeDD,shapeMC2,0,4);
	pair<double,double> etInfo   = GetClosureTestInfo(shapeDD,shapeMC2,1,5);
	pair<double,double> mtInfo   = GetClosureTestInfo(shapeDD,shapeMC2,2,6);
	pair<double,double> ttInfo   = GetClosureTestInfo(shapeDD,shapeMC2,3,7);

        TGraphErrors* emGraph   = GetClosureTestGraph(emInfo);
        TGraphErrors* etGraph   = GetClosureTestGraph(etInfo);
        TGraphErrors* mtGraph   = GetClosureTestGraph(mtInfo);
        TGraphErrors* ttGraph   = GetClosureTestGraph(ttInfo);

	TFile* fileClosure = new TFile((m_Dir+"/ClosureTest"+Idx+".root").c_str(), "RECREATE");
        emGraph->Write("emGraph");
        etGraph->Write("etGraph");
        mtGraph->Write("mtGraph");
        ttGraph->Write("ttGraph");
        return fileClosure;
}




void FRAnalyzer::RetrieveResults(){
	vector<TFile*> file=FinalEstimate();
}
void FRAnalyzer::RetrieveCTResults(){
        TFile* fileC = ClosureTest();
}
