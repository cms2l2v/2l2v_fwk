#include "UserCode/llvv_fwk/interface/MiniEvent.h"

//
void createMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  t->Branch("isData",     &ev.isData,     "isData/O");

  t->Branch("ttbar_nw",        &ev.ttbar_nw,        "ttbar_nw/I");
  t->Branch("ttbar_allmepartons",        &ev.ttbar_allmepartons,        "ttbar_allmepartons/I");
  t->Branch("ttbar_matchmepartons",        &ev.ttbar_matchmepartons,        "ttbar_matchmepartons/I");
  t->Branch("ttbar_w",        ev.ttbar_w,        "ttbar_w[ttbar_nw]/D");
  t->Branch("ttbar_genId",    &ev.ttbar_genId,    "ttbar_genId/I");

  t->Branch("run",       &ev.run,       "run/I");
  t->Branch("event",     &ev.event,     "event/I");
  t->Branch("lumi",      &ev.lumi,      "lumi/I");

  t->Branch("weight", &ev.weight, "weight/D");

  t->Branch("isFiducial",       &ev.isFiducial,       "isFiducial/O");
  t->Branch("muTrigger",        &ev.muTrigger,        "muTrigger/O");
  t->Branch("elTrigger",        &ev.elTrigger,        "elTrigger/O");

  t->Branch("nvtx",      &ev.nvtx,      "nvtx/I");
  t->Branch("pu",      &ev.pu,      "pu/I");
  t->Branch("putrue",      &ev.putrue,      "putrue/I");

  t->Branch("nl",        &ev.nl,        "nl/I");
  t->Branch("isPromptFinalState",        ev.isPromptFinalState,        "isPromptFinalState[nl]/I");
  t->Branch("isDirectPromptTauDecayProductFinalState",        ev.isDirectPromptTauDecayProductFinalState,        "isDirectPromptTauDecayProductFinalState[nl]/I");
  t->Branch("l_id",      ev.l_id,      "l_id[nl]/I");
  t->Branch("l_tightId", ev.l_tightId, "l_tightId[nl]/I");
  t->Branch("l_tightIso",ev.l_tightIso,"l_tightIso[nl]/I");
  t->Branch("l_charge",  ev.l_charge,  "l_charge[nl]/I");
  t->Branch("l_pt",      ev.l_pt,      "l_pt[nl]/D");
  t->Branch("l_eta",     ev.l_eta,     "l_eta[nl]/D");
  t->Branch("l_phi",     ev.l_phi,     "l_phi[nl]/D");
  t->Branch("l_mass",    ev.l_mass,    "l_mass[nl]/D");

  t->Branch("nt"           ,  &ev.nt          , "nt/I");
  t->Branch("t_charge"     ,  ev.t_charge     , "t_charge[nt]/I");
  t->Branch("t_pt"         ,  ev.t_pt         , "t_pt[nt]/D");
  t->Branch("t_eta"        ,  ev.t_eta        , "t_eta[nt]/D");
  t->Branch("t_phi"        ,  ev.t_phi        , "t_phi[nt]/D");
  t->Branch("t_mass"       ,  ev.t_mass       , "t_mass[nt]/D");
  t->Branch("t_leadChHadP" ,  ev.t_leadChHadP , "t_leadChHadP[nt]/D");
  t->Branch("t_leadChHadPt",  ev.t_leadChHadPt, "t_leadChHadPt[nt]/D");
  t->Branch("t_energy"     ,  ev.t_energy     , "t_energy[nt]/D");
  t->Branch("t_et"         ,  ev.t_et         , "t_et[nt]/D");

  t->Branch("t_dmfinding"       ,  ev.t_dmfinding        , "t_dmfinding/D");        
  t->Branch("t_isodb3hits"      ,  ev.t_isodb3hits       , "t_isodb3hits/D");       
  t->Branch("t_againstMu3Loose" ,  ev.t_againstMu3Loose  , "t_againstMu3Loose/D");  
  t->Branch("t_againstMu3Tight" ,  ev.t_againstMu3Tight  , "t_againstMu3Tight/D");  
  t->Branch("t_againstEl5VLoose",  ev.t_againstEl5VLoose , "t_againstEl5VLoose/D");       
  t->Branch("t_againstEl5Loose" ,  ev.t_againstEl5Loose  , "t_againstEl5Loose/D");       
  t->Branch("t_againstEl5Medium",  ev.t_againstEl5Medium , "t_againstEl5Medium/D");        

  t->Branch("me_id",        &ev.me_id,        "me_id/I");
  t->Branch("me_np",        &ev.me_np,        "me_np/I");
  t->Branch("me_pid",       ev.me_pid,        "me_pid/I");
  t->Branch("me_px",       ev.me_px,        "me_px/D");
  t->Branch("me_py",       ev.me_py,        "me_py/D");
  t->Branch("me_pz",       ev.me_pz,        "me_pz/D");
  t->Branch("me_mass",       ev.me_mass,        "me_mass/D");

  t->Branch("nj",        &ev.nj,        "nj/I");
  t->Branch("ngenj",        &ev.ngenj,        "ngenj/I");
  t->Branch("j_area",       ev.j_area,      "j_area[nj]/D");
  t->Branch("j_pt",       ev.j_pt,      "j_pt[nj]/D");
  t->Branch("j_eta",      ev.j_eta,     "j_eta[nj]/D");
  t->Branch("j_phi",      ev.j_phi,     "j_phi[nj]/D");
  t->Branch("j_mass",      ev.j_mass,     "j_mass[nj]/D");
  t->Branch("genj_pt",       ev.genj_pt,      "genj_pt[nj]/D");
  t->Branch("genj_eta",      ev.genj_eta,     "genj_eta[nj]/D");
  t->Branch("genj_phi",      ev.genj_phi,     "genj_phi[nj]/D");
  t->Branch("genj_mass",       ev.genj_mass,      "genj_mass[nj]/D");
  t->Branch("j_csv",      ev.j_csv,     "j_csv[nj]/D");
  t->Branch("j_isbtag",   ev.j_isbtag,  "j_isbtag[nj]/I");
  t->Branch("j_flav",     ev.j_flav,    "j_flav[nj]/I");
  t->Branch("j_hadflav",     ev.j_hadflav,    "j_hadflav[nj]/I");
  t->Branch("j_pid",      ev.j_pid,     "j_pid[nj]/I");

  t->Branch("met_pt",    &ev.met_pt,    "met_pt/D");
  t->Branch("met_phi",   &ev.met_phi,   "met_phi/D");
}

//
void attachToMiniEventTree(TTree *t,MiniEvent_t &ev)
{
  t->SetBranchAddress("isData",     &ev.isData);

  t->SetBranchAddress("ttbar_nw",        &ev.ttbar_nw);
  t->SetBranchAddress("ttbar_allmepartons",        &ev.ttbar_allmepartons);
  t->SetBranchAddress("ttbar_matchmepartons",        &ev.ttbar_matchmepartons);
  t->SetBranchAddress("ttbar_w",        ev.ttbar_w);
  t->SetBranchAddress("ttbar_genId",    &ev.ttbar_genId);

  t->SetBranchAddress("me_id",        &ev.me_id);
  t->SetBranchAddress("me_np",        &ev.me_np);
  t->SetBranchAddress("me_pid",       ev.me_pid);
  t->SetBranchAddress("me_px",       ev.me_px);
  t->SetBranchAddress("me_py",       ev.me_py);
  t->SetBranchAddress("me_pz",       ev.me_pz);
  t->SetBranchAddress("me_mass",       ev.me_mass);

  t->SetBranchAddress("run",       &ev.run);
  t->SetBranchAddress("event",     &ev.event);
  t->SetBranchAddress("lumi",      &ev.lumi);

  t->SetBranchAddress("weight",    &ev.weight);

  t->SetBranchAddress("isFiducial",       &ev.isFiducial);
  t->SetBranchAddress("muTrigger",        &ev.muTrigger);
  t->SetBranchAddress("elTrigger",        &ev.elTrigger);

  t->SetBranchAddress("nvtx",      &ev.nvtx);
  t->SetBranchAddress("pu",      &ev.pu);
  t->SetBranchAddress("putrue",      &ev.putrue);

  t->SetBranchAddress("nl", &ev.nl);
  t->SetBranchAddress("isPromptFinalState",        ev.isPromptFinalState);
  t->SetBranchAddress("isDirectPromptTauDecayProductFinalState",        ev.isDirectPromptTauDecayProductFinalState);
  t->SetBranchAddress("l_id",      ev.l_id);
  t->SetBranchAddress("l_tightId", ev.l_tightId);
  t->SetBranchAddress("l_tightIso",ev.l_tightIso);
  t->SetBranchAddress("l_charge",  ev.l_charge);
  t->SetBranchAddress("l_pt",      ev.l_pt);
  t->SetBranchAddress("l_eta",     ev.l_eta);
  t->SetBranchAddress("l_phi",     ev.l_phi);
  t->SetBranchAddress("l_mass",    ev.l_mass);

  t->SetBranchAddress("nt"                ,  &ev.nt          );
  t->SetBranchAddress("t_charge"          ,  ev.t_charge     );
  t->SetBranchAddress("t_pt"              ,  ev.t_pt         );
  t->SetBranchAddress("t_eta"             ,  ev.t_eta        );
  t->SetBranchAddress("t_phi"             ,  ev.t_phi        );
  t->SetBranchAddress("t_mass"            ,  ev.t_mass       );
  t->SetBranchAddress("t_leadChHadP"      ,  ev.t_leadChHadP );
  t->SetBranchAddress("t_leadChHadPt"     ,  ev.t_leadChHadPt);
  t->SetBranchAddress("t_energy"          ,  ev.t_energy     );
  t->SetBranchAddress("t_et"              ,  ev.t_et         );
  t->SetBranchAddress("t_dmfinding"       ,  ev.t_dmfinding       );
  t->SetBranchAddress("t_isodb3hits"      ,  ev.t_isodb3hits      );
  t->SetBranchAddress("t_againstMu3Loose" ,  ev.t_againstMu3Loose );
  t->SetBranchAddress("t_againstMu3Tight" ,  ev.t_againstMu3Tight );
  t->SetBranchAddress("t_againstEl5VLoose",  ev.t_againstEl5VLoose);
  t->SetBranchAddress("t_againstEl5Loose" ,  ev.t_againstEl5Loose );
  t->SetBranchAddress("t_againstEl5Medium",  ev.t_againstEl5Medium);

  t->SetBranchAddress("ngenj",        &ev.ngenj);
  t->SetBranchAddress("nj",        &ev.nj);
  t->SetBranchAddress("j_area",       ev.j_area);
  t->SetBranchAddress("j_pt",       ev.j_pt);
  t->SetBranchAddress("j_eta",      ev.j_eta);
  t->SetBranchAddress("j_phi",      ev.j_phi);
  t->SetBranchAddress("j_mass",      ev.j_mass);
  t->SetBranchAddress("genj_pt",       ev.genj_pt);
  t->SetBranchAddress("genj_eta",      ev.genj_eta);
  t->SetBranchAddress("genj_phi",      ev.genj_phi);
  t->SetBranchAddress("genj_mass",       ev.genj_mass);
  t->SetBranchAddress("j_csv",      ev.j_csv);
  t->SetBranchAddress("j_isbtag",   ev.j_isbtag);
  t->SetBranchAddress("j_flav",     ev.j_flav);
  t->SetBranchAddress("j_hadflav",     ev.j_hadflav);
  t->SetBranchAddress("j_pid",      ev.j_pid);

  t->SetBranchAddress("met_pt",    &ev.met_pt);
  t->SetBranchAddress("met_phi",   &ev.met_phi);
}

