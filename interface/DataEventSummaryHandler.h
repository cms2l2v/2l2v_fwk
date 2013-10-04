#ifndef dataeventsummaryhandler_h
#define dataeventsummaryhandler_h

#include "UserCode/llvv_fwk/interface/DataEventSummary.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"

#include "TTree.h"

namespace data
{
  class PhysicsObject_t : public LorentzVector
  {
  public:

    PhysicsObject_t()
      {
	SetPxPyPzE(0,0,0,0);
      }

    PhysicsObject_t(Float_t px, Float_t py, Float_t pz, Float_t en)
      {
	SetPxPyPzE(px,py,pz,en);
      } 

    PhysicsObject_t(const PhysicsObject_t & orig)
      {
	SetPxPyPzE(orig.px(),orig.py(),orig.pz(),orig.energy());
	info=orig.info;
	vals=orig.vals;
	flags=orig.flags;
	objs=orig.objs;
      }
    
    inline void setVal(TString key, Double_t val)                  { vals[key]=val; }
    inline const Double_t &getVal(TString key)                     { return vals.find(key)->second; }
    inline void setFlag(TString key, Bool_t val)                   { flags[key]=val; }
    inline const Bool_t &getFlag(TString key)                      { return flags.find(key)->second; }
    inline void set(TString key, Int_t val)                        { info[key]=val;}
    inline const Int_t &get(TString key)                           { return info.find(key)->second; } 
    inline void setObject(TString key, const PhysicsObject_t &obj) { objs.insert( std::pair<TString, PhysicsObject_t>(key,obj) ); }
    inline const PhysicsObject_t &getObject(TString key)           { return objs.find(key)->second; }

    static bool sortByPt(const data::PhysicsObject_t &a, const data::PhysicsObject_t &b)  { return a.pt()>b.pt(); }
    static bool sortByCSV(const data::PhysicsObject_t &a, const data::PhysicsObject_t &b) { return a.vals.find("csv")->second>b.vals.find("csv")->second; }

    std::map<TString,Int_t > info;
    std::map<TString,Double_t> vals;
    std::map<TString,Bool_t> flags;
    std::map<TString, PhysicsObject_t> objs;
  };

  typedef std::vector<PhysicsObject_t> PhysicsObjectCollection_t;

}


class DataEventSummaryHandler{

 public:

  //c/dtor
  DataEventSummaryHandler()  { };
  ~DataEventSummaryHandler() { };
  
  //initialize
  bool init(TTree *t,bool needsToRecreate=true);
  bool attach(TTree *t,bool readPFbranch=true);

  //r/w mode
  void fill() { if(t_) t_->Fill(); }
  int getEntries() { return (t_ ? t_->GetEntriesFast() : 0); }
  void getEntry(int ientry) { evSummary_.reset(); if(t_) t_->GetEntry(ientry);  }

  //getters
  TTree *getTree() { return t_; }
  DataEventSummary &getEvent() { return evSummary_; }
  
  enum PhysicsCode {LEPTONS, PHOTONS, JETS, MET, PFCANDIDATES, GENPARTICLES};
  data::PhysicsObjectCollection_t getPhysicsObject(int code);

  //the tree and the data holders
  TTree *t_;
  DataEventSummary evSummary_;
};


namespace utils
{
  namespace cmssw
  {
    //based on http://arxiv.org/pdf/1001.5027v3.pdf
    enum JetPullTypes {toJ1_ALL,toJ1_CH,toJ1_NEUT,toJ2_ALL,toJ2_CH,toJ2_NEUT};
    std::vector<float> pullDeltaTheta(data::PhysicsObject_t &jet1, data::PhysicsObject_t &jet2, data::PhysicsObjectCollection_t &pf);
    void setJetPulls(data::PhysicsObject_t &jet, data::PhysicsObjectCollection_t &pf);
    void setJetDirections(data::PhysicsObject_t &jet, data::PhysicsObjectCollection_t &pf);

    //set new jet energy corrections
    void updateJEC(data::PhysicsObjectCollection_t &jets, FactorizedJetCorrector *jesCor, JetCorrectionUncertainty *totalJESUnc, float rho, int nvtx,bool isMC);

    //apply MET variations
    enum METvariations { NOMINAL, JERUP, JERDOWN, JESUP, JESDOWN, UMETUP, UMETDOWN, LESUP, LESDOWN };
    std::vector<LorentzVector> getMETvariations(data::PhysicsObject_t &rawMETP4, data::PhysicsObjectCollection_t &jets, data::PhysicsObjectCollection_t &leptons, bool isMC);
  }
}

#endif
