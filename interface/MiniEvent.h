#ifndef _minievent_h_
#define _minievent_h_

#include "TTree.h"

struct MiniEvent_t
{
  Bool_t isData;
  Int_t run,event,lumi;

  Double_t weight; // It includes pileup and NNLO weight, when applicable.

  Int_t ttbar_nw, ttbar_allmepartons, ttbar_matchmepartons,ttbar_genId;
  Double_t ttbar_w[500];
  Bool_t isFiducial;

  Int_t nvtx,pu,putrue;
  Double_t rho;
  
  Bool_t muTrigger,elTrigger;
  
  Int_t nl;
  Int_t isPromptFinalState[10], isDirectPromptTauDecayProductFinalState[10], l_tightId[10], l_tightIso[10];
  Int_t l_id[10],l_charge[10];
  Double_t l_pt[10],l_eta[10],l_phi[10], l_mass[10];

  Int_t nt;
  Int_t t_charge[10];
  Double_t t_pt[10], t_eta[10], t_phi[10], t_mass[10], t_leadChHadP[10], t_leadChHadPt[10], t_energy[10], t_et[10];

  Int_t nj,ngenj;
  Double_t j_pt[20],j_eta[20],j_phi[20],j_mass[20],j_area[20];
  Double_t genj_pt[20],genj_eta[20],genj_phi[20],genj_mass[20];
  Double_t j_csv[20];
  Int_t j_flav[20],j_pid[20],j_hadflav[20], j_isbtag[20];

  Double_t met_pt,met_phi;

  Int_t me_id,me_np,me_pid[25];
  Double_t me_px[25],me_py[25],me_pz[25],me_mass[25];
};

void createMiniEventTree(TTree *t,MiniEvent_t &ev);
void attachToMiniEventTree(TTree *t, MiniEvent_t &ev);

#endif
