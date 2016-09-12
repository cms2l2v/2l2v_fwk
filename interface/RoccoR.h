#ifndef ElectroWeakAnalysis_RoccoR
#define ElectroWeakAnalysis_RoccoR
#include "TRandom3.h"
#include "TMath.h"

#include "iostream"
#include "fstream"
#include "sstream"

using namespace std;
using std::cout;
using std::endl;
using std::vector;

struct CrystalBall_2016{
    const double pi    = TMath::Pi();
    const double SPiO2 = sqrt(TMath::Pi()/2.0);
    const double S2    = sqrt(2.0);


    double m;
    double s;
    double a;
    double n;

    double B;
    double C;
    double D;
    double N;

    double NA;
    double Ns;
    double NC;
    double F;
    double G;
    double k;

    double cdfMa;
    double cdfPa;

    CrystalBall_2016(){
	init(0, 1, 10, 10);
    }
    CrystalBall_2016(double m_, double s_, double a_, double n_){
	init(m_, s_, a_, n_);
    }

    void init(double m_, double s_, double a_, double n_){
	m=m_;
	s=s_;
	a=a_;
	n=n_;

	double fa   = fabs(a);
	double expa = exp(-fa*fa/2);
	double A    = pow(n/fa, n)*expa;
	double C1   = n/fa/(n-1)*expa; 
	double D1   = 2*SPiO2*erf(fa/S2);

	B  = n/fa-fa;
	C  = (D1+2*C1)/C1;   
	D  = (D1+2*C1)/2;   

	N  = 1.0/s/(D1+2*C1); 
	k  = 1.0/(n-1);  

	NA = N*A;       
	Ns = N*s;       
	NC = Ns*C1;     
	F  = 1-fa*fa/n; 
	G  = s*n/fa;    

	cdfMa=cdf(m-a*s);
	cdfPa=cdf(m+a*s);
    }

    double pdf(double x){ 
	double d=(x-m)/s;
	if(d<-a) return NA*pow(B-d, -n);
	if(d> a) return NA*pow(B+d, -n);
	return N*exp(-d*d/2);
    }

    double cdf(double x){
	double d = (x-m)/s;
	if(d<-a) return NC / pow(F-s*d/G, n-1);
	if(d> a) return NC * (C - pow(F+s*d/G, 1-n) );
	return Ns*(D-SPiO2*erf(-d/S2));
    }

    double invcdf(double u){
	if(u<cdfMa) return m + G*(F - pow(NC/u,    k) );
	if(u>cdfPa) return m - G*(F - pow(C-u/NC, -k) );
	return m - S2*s*TMath::ErfInverse((D - u/Ns ) / SPiO2);
    }
};

class RocRes{
    private:
	static const int NMAXETA=12;
	static const int NMAXTRK=12;

	int NETA;
	int NTRK;
	int NMIN;

	double BETA[NMAXETA+1];
	double ntrk[NMAXETA][NMAXTRK+1];
	double dtrk[NMAXETA][NMAXTRK+1];

	double width[NMAXETA][NMAXTRK];
	double alpha[NMAXETA][NMAXTRK];
	double power[NMAXETA][NMAXTRK];

	double rmsA[NMAXETA][NMAXTRK];
	double rmsB[NMAXETA][NMAXTRK];
	double rmsC[NMAXETA][NMAXTRK];

	double kDat[NMAXETA];
	double kRes[NMAXETA];

	TRandom3 random; //only used while deriving corrections now
	
	int getBin(double x, const int NN, const double *b);


    public:
	enum TYPE {MC, Data, Extra};

	CrystalBall_2016  cb[NMAXETA][NMAXTRK]; 

	RocRes();
	int getEtaBin(double feta);
	int getNBinDT(double v, int H);
	int getNBinMC(double v, int H);
	void dumpParams();
	void init(string filename);

	~RocRes(){}

	double Sigma(double pt, int H, int F);
	double kSpreadDet(double gpt, double rpt, double eta, int nlayers, double w);
	double kSmearDet(double pt, double eta, TYPE type, double v, double u);
	double kExtraDet(double pt, double eta, int nlayers, double u, double w);
	void fillFitData(int &H, int &F, int &D, double &xmc, double &xdt, double &Rmc, double &Rdt, double pt, double eta);
};


class RocOne{
    private:
	enum TYPE{MC, DT};
	static const int NMAXETA=24;
	static const int NMAXPHI=16;
	static const double MPHI;

	int NETA;
	int NPHI;

	double BETA[NMAXETA+1];
	double DPHI;

	double M[2][NMAXETA][NMAXPHI];
	double A[2][NMAXETA][NMAXPHI];
	double D[2][NMAXETA];

	RocRes RR;

	int getBin(double x, const int NN, const double *b);
	int getBin(double x, const int nmax, const double xmin, const double dx);

    public:
	RocOne();
	~RocOne(){}
	RocOne(string filename, int iTYPE=0, int iSYS=0, int iMEM=0);
	bool checkSYS(int iSYS, int iMEM, int kSYS=0, int kMEM=0);
	bool checkTIGHT(int iTYPE, int iSYS, int iMEM, int kTYPE=0, int kSYS=0, int kMEM=0);
	void init(string filename, int iTYPE=0, int iSYS=0, int iMEM=0);

	double kScaleDT(int Q, double pt, double eta, double phi);
	double kScaleMC(int Q, double pt, double eta, double phi, double kSMR=1);
	double kScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers, double u, double w);
	double kScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers, double w);
	double kGenSmearDet(double pt, double eta, double v, double u);
};


class RoccoR{
    public:
	static const int Nset=5;
	enum TYPE{Default, Stat, ModPt, CorDM, FitDM};
	static const int Nmem[Nset];

	RoccoR(); //temporary, will change 
	RoccoR(std::string dirname); //temporary, will change 
	~RoccoR();
	
	double kScaleDT(int Q, double pt, double eta, double phi, TYPE T, int m);

	double kScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers);	      //only for default, non-reproducible
	double kScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers); //only for default, non-reproducible

	double kScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers, double u, double w, TYPE T, int m);  
	double kScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers, double w, TYPE T, int m); 

	std::vector<std::vector<double> > kkScaleAndSmearMCDet(int Q, double pt, double eta, double phi, int nlayers);
	std::vector<std::vector<double> > kkScaleFromGenMCDet(int Q, double pt, double eta, double phi, double gpt, int nlayers);

    private:
	TRandom3 random;
	RocOne *RC[Nset][100];

};


#endif

