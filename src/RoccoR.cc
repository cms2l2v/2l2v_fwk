#include "TRandom3.h"
#include "TMath.h"
#include "UserCode/llvv_fwk/interface/RoccoR.h"
//#include "__utils__/RoccoR.h"


int RocRes::getBin(double x, const int NN, const double *b){
    for(int i=0; i<NN; ++i) if(x<b[i+1]) return i;
    return NN-1;
}

RocRes::RocRes():NETA(1),NTRK(1),NMIN(1){
    for(int H=0; H<NMAXETA; ++H){
	BETA[H]=0;
	kDat[H]=1.0;
	kRes[H]=1.0; //this is important
	for(int F=0; F<NMAXTRK+1; ++F){
	    ntrk[H][F]=0;
	    dtrk[H][F]=0;
	}
	for(int F=0; F<NMAXTRK; ++F){
	    width[H][F]=1;
	    alpha[H][F]=10;
	    power[H][F]=10;
	    cb[H][F].init(0.0, width[H][F], alpha[H][F], power[H][F]);
	}
    }
    BETA[NMAXETA]=0;
}

int RocRes::getEtaBin(double feta){
    return getBin(feta,NETA,BETA);
}

int RocRes::getNBinDT(double v, int H){
    return getBin(v,NTRK,dtrk[H]);
}

int RocRes::getNBinMC(double v, int H){
    return getBin(v,NTRK,ntrk[H]);
}

void RocRes::dumpParams(){
    cout << NMIN << endl;
    cout << NTRK << endl;
    cout << NETA << endl;
    for(int H=0; H<NETA+1; ++H) cout << BETA[H] << " ";
    cout << endl;
    for(int H=0; H<NETA; ++H){
	for(int F=0; F<NTRK; ++F){
	    cout << Form("%8.4f %8.4f %8.4f | ", width[H][F], alpha[H][F], power[H][F]);
	}
	cout << endl;
    }
    for(int H=0; H<NETA; ++H){
	for(int F=0; F<NTRK+1; ++F){
	    cout << Form("%8.4f %8.4f| ", ntrk[H][F], dtrk[H][F]);
	}
	cout << endl;
    }
    for(int H=0; H<NETA; ++H){
	for(int F=0; F<NTRK; ++F){
	    cb[H][F].init(0.0, width[H][F], alpha[H][F], power[H][F]);
	    cout << Form("%8.4f %8.4f %8.4f | ", rmsA[H][F], rmsB[H][F], rmsC[H][F]);
	}
	cout << endl;
    }
}


	
void RocRes::init(string filename){
  std::ifstream in(filename);
  // char buffer[256];
  char tag[4];
  int type, sys, mem, isdt, var, bin;	
  std::string s;
  while(std::getline(in, s)){
    stringstream ss(s); 
    if(s.substr(0,4)=="RMIN")       ss >> tag >> NMIN;
    else if(s.substr(0,4)=="RTRK")  ss >> tag >> NTRK;
    else if(s.substr(0,4)=="RETA")  {
      ss >> tag >> NETA;
      for(int i=0; i< NETA+1; ++i) ss >> BETA[i];
    }
    else if(s.substr(0,1)=="R")  {
      ss >> tag >> type >> sys >> mem >> isdt >> var >> bin; 
      if(var==0) for(int i=0; i<NTRK; ++i) ss >> rmsA[bin][i];  
      if(var==1) for(int i=0; i<NTRK; ++i) ss >> rmsB[bin][i];  
      if(var==2) for(int i=0; i<NTRK; ++i) {
	  ss >> rmsC[bin][i];  
	  rmsC[bin][i]/=100;
	}
      if(var==3) for(int i=0; i<NTRK; ++i) ss >> width[bin][i];  
      if(var==4) for(int i=0; i<NTRK; ++i) ss >> alpha[bin][i];  
      if(var==5) for(int i=0; i<NTRK; ++i) ss >> power[bin][i];  
    }
    else if(s.substr(0,1)=="T")  {
      ss >> tag >> type >> sys >> mem >> isdt >> var >> bin; 
      if(isdt==0) for(int i=0; i<NTRK+1; ++i) ss >> ntrk[bin][i];  
      if(isdt==1) for(int i=0; i<NTRK+1; ++i) ss >> dtrk[bin][i];  
    }
    else if(s.substr(0,1)=="F")  {
      ss >> tag >> type >> sys >> mem >> isdt >> var >> bin; 
      if(var==0){
	if(isdt==0) for(int i=0; i<NETA; ++i) ss >> kRes[i];  
	if(isdt==1) for(int i=0; i<NETA; ++i) ss >> kDat[i];
      }
    }
  }

  for(int H=0; H<NETA; ++H){
    for(int F=0; F<NTRK; ++F){
      cb[H][F].init(0.0, width[H][F], alpha[H][F], power[H][F]);
    }
  }
  in.close();
}

double RocRes::Sigma(double pt, int H, int F){
    double dpt=pt-45;
    return rmsA[H][F] + rmsB[H][F]*dpt + rmsC[H][F]*dpt*dpt;
}

double RocRes::kSpreadDet(double gpt, double rpt, double eta, int nlayers, double w){
    int     H = getBin(fabs(eta), NETA, BETA);
    int     F = nlayers-NMIN;
    double  v = ntrk[H][F]+(ntrk[H][F+1]-ntrk[H][F])*w;
    int     D = getBin(v, NTRK, dtrk[H]);
    double  kold = gpt / rpt;
    double  u = cb[H][F].cdf( (kold-1.0)/kRes[H]/Sigma(gpt,H,F) ); 
    double  knew = 1.0 + kDat[H]*Sigma(gpt,H,D)*cb[H][D].invcdf(u);
    if(knew<0) return 1.0;
    return kold/knew;
}

double RocRes::kSmearDet(double pt, double eta, TYPE type, double v, double u){
    int H = getBin(fabs(eta), NETA, BETA);
    const double (&trk) [NMAXTRK+1] = type==Data ? dtrk[H] : ntrk[H];
    int F = getBin(v, NTRK, trk);
    double K = type==Data ? kDat[H] : kRes[H]; 
    double x = K*Sigma(pt, H, F)*cb[H][F].invcdf(u);
    return 1.0/(1.0+x);
}

double RocRes::kExtraDet(double pt, double eta, int nlayers, double u, double w){
    int H = getBin(fabs(eta), NETA, BETA);
    int F = nlayers-NMIN;
    double  v = ntrk[H][F]+(ntrk[H][F+1]-ntrk[H][F])*w;
    int     D = getBin(v, NTRK, dtrk[H]);
    double RD = kDat[H]*Sigma(pt, H, D);
    double RM = kRes[H]*Sigma(pt, H, F);
    double x = RD>RM ? sqrt(RD*RD-RM*RM)*cb[H][F].invcdf(u) : 0;
    if(x<=-1) return 1.0;
    return 1.0/(1.0 + x); 
}

void RocRes::fillFitData(int &H, int &F, int &D, double &xmc, double &xdt, double &Rmc, double &Rdt, double pt, double eta){
    H = getBin(fabs(eta), NETA, BETA);
    double u=random.Rndm();
    F = getBin(u, NTRK, ntrk[H]); 
    D = getBin(u, NTRK, dtrk[H]);
    Rmc = Sigma(pt, H, F); 
    Rdt = Sigma(pt, H, D); 
    double v=random.Rndm();
    xmc = Rmc*cb[H][F].invcdf(v); 
    xdt = Rdt*cb[H][D].invcdf(v); 
    return;
}











const double RocOne::MPHI=-TMath::Pi();

int RocOne::getBin(double x, const int NN, const double *b){
    for(int i=0; i<NN; ++i) if(x<b[i+1]) return i;
    return NN-1;
}

int RocOne::getBin(double x, const int nmax, const double xmin, const double dx){
    int ibin=(x-xmin)/dx;
    if(ibin<0) return 0; 
    if(ibin>=nmax) return nmax-1;
    return ibin;
}

RocOne::RocOne():NETA(1),NPHI(1){
    DPHI=2*TMath::Pi()/NPHI;
    for(int H=0; H<NMAXETA; ++H){
	BETA[H]=0;
	D[MC][H]=1.0;
	D[DT][H]=1.0;
	for(int F=0; F<NMAXPHI; ++F){
	    for(int T=0; T<2; ++T){
		M[T][H][F]=1;
		A[T][H][F]=0;
	    }
	}
    }
    BETA[NMAXETA]=0;
}

RocOne::RocOne(std::string filename, int iTYPE, int iSYS, int iMEM){
    init(filename, iTYPE, iSYS, iMEM);
}


bool RocOne::checkSYS(int iSYS, int iMEM, int kSYS, int kMEM){
    if(iSYS==0 && iMEM==0)	      return true;
    if(iSYS==kSYS && iMEM==kMEM)      return true;
    return false;
}

bool RocOne::checkTIGHT(int iTYPE, int iSYS, int iMEM, int kTYPE, int kSYS, int kMEM){
    if(iTYPE!=kTYPE) return false;
    if(iSYS!=kSYS)   return false;
    if(iMEM!=kMEM)   return false;
    return true;
}

void RocOne::init(std::string filename, int iTYPE, int iSYS, int iMEM){

    RR.init(filename);

    std::ifstream in(filename);
    char tag[4];
    int type, sys, mem, isdt, var, bin;	

    std::string s;
    while(std::getline(in, s)){
	stringstream ss(s); 
	if(s.substr(0,4)=="CPHI")       {
	    ss >> tag >> NPHI;
	    DPHI=2*TMath::Pi()/NPHI;
	}
	else if(s.substr(0,4)=="CETA")  {
	    ss >> tag >> NETA;
	    for(int i=0; i< NETA+1; ++i) ss >> BETA[i];
	}
	else if(s.substr(0,1)=="C")  {
	    ss >> tag >> type >> sys >> mem >> isdt >> var >> bin; 
	    if(!checkTIGHT(type,sys,mem,iTYPE,iSYS,iMEM)) continue;
	    if(var==0) for(int i=0; i<NPHI; ++i) { ss >> M[isdt][bin][i]; M[isdt][bin][i]/=100; M[isdt][bin][i]+=1.0;} 
	    if(var==1) for(int i=0; i<NPHI; ++i) { ss >> A[isdt][bin][i]; A[isdt][bin][i]/=100;} 

	}
	else if(s.substr(0,1)=="F")  {
	    ss >> tag >> type >> sys >> mem >> isdt >> var >> bin; 
	    if(var==1){
		for(int i=0; i<NETA; ++i) { 
		    ss >> D[isdt][i];  
		    D[isdt][i]/=10000;
		    D[isdt][i]+=1.0;
		}
	    }
	}
    }
    in.close();
}

double RocOne::kScaleDT(int Q, double pt, double eta, double phi){
    int H=getBin(eta, NETA, BETA);
    int F=getBin(phi, NPHI, MPHI, DPHI);
    double m=M[DT][H][F];
    double a=A[DT][H][F];
    double d=D[DT][H];
    double k=d/(m+Q*a*pt);
    return k;
}


double RocOne::kScaleMC(int Q, double pt, double eta, double phi, double kSMR){
    int H=getBin(eta, NETA, BETA);
    int F=getBin(phi, NPHI, MPHI, DPHI);
    double m=M[MC][H][F];
    double a=A[MC][H][F];
    double d=D[MC][H];
    double k=d/(m+Q*a*pt);
    return k*kSMR;
}

double RocOne::kScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers, double u, double w){
    double k=kScaleMC(Q, pt, eta, phi);
    return k*RR.kExtraDet(k*pt, eta, nlayers, u, w);
}


double RocOne::kScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers, double w){
    double k=kScaleMC(Q, pt, eta, phi);
    return k*RR.kSpreadDet(gpt, k*pt, eta, nlayers, w);
}


double RocOne::kGenSmearDet(double pt, double eta, double v, double u){
    return RR.kSmearDet(pt, eta, RocRes::Data, v, u);
}








const int RoccoR::Nmem[RoccoR::Nset]={1, 100, 1, 41, 21};

RoccoR::RoccoR(){
    for(int s=0; s<Nset; ++s){
	for(int m=0; m<Nmem[s]; ++m){
	    RC[s][m]=0;
	}
    }
}


RoccoR::RoccoR(std::string dirname){
    for(int s=0; s<Nset; ++s){
	for(int m=0; m<Nmem[s]; ++m){
	    std::string inputfile=Form("%s/rc_%d_%d.txt", dirname.c_str(), s, m);
	    RC[s][m]=new RocOne(inputfile, 0, s, m);
	}
    }
}

RoccoR::~RoccoR(){
    for(int s=0; s<Nset; ++s){
	for(int m=0; m<Nmem[s]; ++m){
	    if(RC[s][m]) delete RC[s][m];
	}
    }
}


	
double RoccoR::kScaleDT(int Q, double pt, double eta, double phi, TYPE T, int m){
    return RC[T][m]->kScaleDT(Q, pt, eta, phi);
}


double RoccoR::kScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers){
    return RC[Default][0]->kScaleAndSmearMCDet(Q, pt, eta, phi, nlayers, random.Rndm(), random.Rndm());

}
	
double RoccoR::kScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers){
    return RC[Default][0]->kScaleFromGenMCDet(Q, pt, eta, phi, gpt, nlayers, random.Rndm());
}

double RoccoR::kScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers, double u, double w, TYPE T, int m){
    return RC[T][m]->kScaleAndSmearMCDet(Q, pt, eta, phi, nlayers, u, w);
}

double RoccoR::kScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers, double w, TYPE T, int m){
    return RC[T][m]->kScaleFromGenMCDet(Q, pt, eta, phi, gpt, nlayers, w);
}

std::vector<std::vector<double> > RoccoR::kkScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers){

    double u=random.Rndm();
    double w=random.Rndm();

    std::vector<std::vector<double> > result;
    for(int s=0; s<Nset; ++s){
	std::vector<double> d;
	for(int m=0; m<Nmem[s]; ++m){
	    double k=kScaleAndSmearMCDet(Q, pt, eta, phi, nlayers, u,  w, TYPE(s), m);
	    d.push_back(k);
	}
	result.push_back(d);
    }
    return result;
}

std::vector<std::vector<double> > RoccoR::kkScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers){
    // double u=random.Rndm();
    double w=random.Rndm();

    std::vector<std::vector<double> > result;
    for(int s=0; s<Nset; ++s){
	std::vector<double> d;
	for(int m=0; m<Nmem[s]; ++m){
	    double k=kScaleFromGenMCDet(Q, pt, eta, phi, gpt, nlayers, w, TYPE(s),  m);
	    d.push_back(k);
	}
	result.push_back(d);
    }
    return result;
}






