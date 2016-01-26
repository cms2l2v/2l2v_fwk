#ifndef ElectroWeakAnalysis_MuResolution
#define ElectroWeakAnalysis_MuResolution

#include "TRandom3.h"
#include "TMath.h"

struct CrystalBall{ 
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
  
  CrystalBall():m(0),s(1),a(10),n(10){}
    CrystalBall(double m_, double s_, double a_, double n_){
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
    
    NA   = N*A;       
    Ns   = N*s;       
    NC   = Ns*C1;     
    F    = 1-fa*fa/n; 
    G    = s*n/fa;    
    
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


class muresolution {
  
 public:
  
  static const int NETA=12;
  static const int NTRK=12;
  static const int NMIN= 6; //minimum number of layers;
  
  enum TYPE {MC, Data, Extra};
  
  muresolution(){
    for(int H=0; H<NETA; ++H){
      for(int F=0; F<NTRK; ++F){
	cb[H][F].init(0.0, width[H][F], alpha[H][F], power[H][F]);
      }
    }
  }
  
  ~muresolution(){}
  
  double Sigma(double pt, int H, int F){
    double dpt=pt-45;
    return rmsA[H][F] + rmsB[H][F]*dpt + rmsC[H][F]*dpt*dpt;
  }
  
  double kSmear(double pt, double eta, int nlayers=0, TYPE type=Extra){
    int H = getBin(fabs(eta), NETA, BETA);
    int F = nlayers-NMIN;
    const double (&trk) [NTRK+1] = type==Data ? dtrk[H] : mtrk[H];
    if (F<0 || F>=NTRK) F = getBin(random.Rndm(), NTRK, trk);
    
    double x=-2;
    while(1+x<=0){
      double u=random.Rndm();
      if(type==Extra){ 
	double  v = random.Uniform(mtrk[H][F], mtrk[H][F+1]);
	int     D = getBin(v, NTRK, dtrk[H]);
	double RD = kDat[H]*Sigma(pt, H, D);
	double RM = kRes[H]*Sigma(pt, H, F);
	if(RD>RM) x = sqrt(RD*RD-RM*RM)*cb[H][F].invcdf(u);
	else      x = 0;
      }
      else if(type==Data) x = kDat[H]*Sigma(pt, H, F)*cb[H][F].invcdf(u); 
      else		  x = kRes[H]*Sigma(pt, H, F)*cb[H][F].invcdf(u);
    }
    
    return 1.0/(1.0 + x);
  }
  
  double kSpread(double gpt, double rpt, double eta, int nlayers){
    int     H = getBin(fabs(eta), NETA, BETA);
    int     F = nlayers-NMIN;
    double  v = random.Uniform(mtrk[H][F], mtrk[H][F+1]);
    int     D = getBin(v, NTRK, dtrk[H]);
    double  u = random.Rndm();
    double xd = kDat[H] * Sigma(gpt, H, D) * cb[H][D].invcdf(u);
    double xm = kRes[H] * Sigma(gpt, H, F) * cb[H][F].invcdf(u);
    double kold = gpt / rpt;
    double knew = 1.0 + (kold-1.0)*xd/xm; 
    if(knew<0) return kSmear(rpt, eta, nlayers, Extra);
    return kold/knew;
  }
  
  //static const double BETA[NETA+1];
  const double BETA[NETA+1] = {0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.85, 2.10, 2.25, 2.40};

  static const double mtrk[NETA][NTRK+1];
  static const double dtrk[NETA][NTRK+1];
  
  static const double width[NETA][NTRK];
  static const double alpha[NETA][NTRK];
  static const double power[NETA][NTRK];
  
  static const double rmsA[NETA][NTRK];
  static const double rmsB[NETA][NTRK];
  static const double rmsC[NETA][NTRK];
  
  //static const double kDat[NETA];
  //static const double kRes[NETA];
  const double kDat[NETA]={0.957838,0.998524,0.991463,0.949789,0.995753,0.9626,0.91684,0.891457,0.961592,1.03829*0.948057,0.969482*0.97541,1.02167*0.952802};
  const double kRes[NETA]={1.01101,0.928442,0.887687,0.938112,0.86722,0.86134,0.893858,0.89179,0.923499,0.878776*1.03576,0.85004*1.02033,0.928083*0.975536};
  
  CrystalBall  cb[NETA][NTRK]; //base
  
  TRandom3 random;
  
  int getBin(double x, const int NN, const double *b){
    for(int i=0; i<NN; ++i) if(x<b[i+1]) return i;
    return NN-1;
  }
  
};

#endif
