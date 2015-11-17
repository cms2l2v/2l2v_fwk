#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "UserCode/llvv_fwk/interface/MiniEvent.h"
#include "UserCode/llvv_fwk/interface/ReadTree.hh"
#include "UserCode/llvv_fwk/interface/BtagUncertaintyComputer.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"

using namespace std;

Int_t getBtagCatBinToFill(Int_t nBtags,Int_t nJets)
{
  Int_t nJetsBin(nJets>4 ? 4 : nJets);

  Int_t binToFill(nBtags>=2?2:nBtags);
  binToFill+=3*(nJetsBin-1);
  return binToFill;
}

Double_t computeMT(TLorentzVector& a, TLorentzVector& b)
{
  return TMath::Sqrt(2*a.Pt()*b.Pt()*(1-TMath::Cos(a.DeltaPhi(b))));
}

void ReadTree(TString filename,
	      TString outname,
	      Int_t channelSelection, 
	      Int_t chargeSelection, 
	      TH1D *normH, 
	      Bool_t isTTbar,
	      Bool_t debug,
	      FlavourSplitting flavourSplitting,
	      GenWeightMode genWgtMode,
	      TGraph* puWgtGr, TGraph* puUpWgtGr, TGraph* puDownWgtGr)
{

  TString systs[]={"nom",
		   "puUp","puDown",
		   "muEffUp","muEffDown",
		   "eEffUp","eEffDown",
		   "qcdScaleDown","qcdScaleUp",
		   "umetDown", "umetUp",
		   "jesDown","jesUp",
		   "jerDown","jerUp",
		   "beffDown","beffUp",
		   "mistagDown","mistagUp"};



  //book histograms
  std::map<TString, TH1*> allPlots;

  std::map<Int_t,Double_t> lumiMap=lumiPerRun();
  for(Int_t ij=1; ij<=4; ij++)
    {
      for(Int_t itag=-1; itag<=2; itag++)
	{
	  if(itag>ij) continue;
	  TString tag(itag<0 ? Form("%dj",ij) : Form("%dj%dt",ij,itag));
	  allPlots["ratevsrun_"+tag] = new TH1D("ratevsrun_"+tag,";Run number; Events/pb",lumiMap.size(),0,lumiMap.size());
	  Int_t runCtr(0);
	  for(std::map<Int_t,Double_t>::iterator it=lumiMap.begin(); it!=lumiMap.end(); it++,runCtr++)
	    allPlots["ratevsrun_"+tag]->GetXaxis()->SetBinLabel(runCtr+1,Form("%d",it->first));
	  allPlots["lpt_"+tag]  = new TH1D("lpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["lsip3d_"+tag]  = new TH1D("lsip3d_"+tag,";3d impact parameter significance;Events" ,40,0.,20.);
	  allPlots["lchiso_"+tag]  = new TH1D("lchiso_"+tag,";Charged hadron isolation [GeV];Events" ,25,0.,50.);
	  allPlots["lchreliso_"+tag]  = new TH1D("lchreliso_"+tag,";Charged hadron relative isolation;Events" ,25,0.,0.2);
	  allPlots["leta_"+tag] = new TH1D("leta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["jpt_"+tag]  = new TH1D("jpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
	  allPlots["jeta_"+tag] = new TH1D("jeta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
	  allPlots["ht_"+tag]   = new TH1D("ht_"+tag,";H_{T} [GeV];Events",40,0,800);
	  allPlots["csv_"+tag]  = new TH1D("csv_"+tag,";CSV discriminator;Events",100,0,1.0);
	  allPlots["nvtx_"+tag] = new TH1D("nvtx_"+tag,";Vertex multiplicity;Events" ,50,0.,50.);
	  allPlots["met_"+tag]  = new TH1D("metpt_"+tag,";Missing transverse energy [GeV];Events" ,20,0.,300.);
	  allPlots["metphi_"+tag] = new TH1D("metphi_" + tag,";MET #phi [rad];Events" ,50,-3.2,3.2);
	  for(size_t isyst=0; isyst<sizeof(systs)/sizeof(TString); isyst++)
	    {
	      allPlots["mt_"+systs[isyst]+"_"+tag]     = new TH1D("mt_"+systs[isyst]+"_"+tag,";Transverse Mass [GeV];Events" ,20,0.,200.);
	      allPlots["minmlb_"+systs[isyst]+"_"+tag] = new TH1D("minmlb_"+systs[isyst]+"_"+tag,";min Mass(lepton,b) [GeV];Events" ,25,0.,250.);
	    }
	}
    }

  //category counting
  for(size_t i=0; i<sizeof(systs)/sizeof(TString); i++)
    {
      allPlots["njetsnbtags_"+systs[i]] = new TH1D("njetsnbtags_"+systs[i],";Category;Events" ,12,0.,12.);
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(1, "1j,=0b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(2, "1j,=1b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(3, "");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(4, "2j,=0b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(5, "2j,=1b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(6, "2j,#geq2b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(7, "3j,=0b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(8, "3j,=1b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(9, "3j,#geq2b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(10,"4j,=0b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(11,"4j,=1b");
      allPlots["njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(12,"4j,#geq2b");


      allPlots["eventflowold_"+systs[i]] = new TH1D("eventflowold_"+systs[i],";;Events", 5, 0., 5.); 
      allPlots["eventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(1, "1 #it{l}, #geq 2 jets");
      allPlots["eventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(2, "E_{T}^{miss}");
      allPlots["eventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(3, "#geq 1 b-tag");
      allPlots["eventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(4, "1 #tau");
      allPlots["eventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(5, "op. sign");
      
      allPlots["eventflownocat_"+systs[i]] = new TH1D("eventflownocat_"+systs[i], ";;Events", 4, 0., 4.);
      allPlots["eventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(1, "1 #it{l}, #geq 2 jets");
      allPlots["eventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(2, "E_{T}^{miss}");
      allPlots["eventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(3, "1 #tau");
      allPlots["eventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(4, "op. sign");

      allPlots["eventflowcat_"+systs[i]] = new TH1D("eventflowcat_"+systs[i], ";;Events", 7, 0., 7.);
      allPlots["eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(1, "1 #it{l}, #geq 2 jets");
      allPlots["eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(2, "E_{T}^{miss}");
      allPlots["eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(3, "1 #tau");
      allPlots["eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(4, "op. sign");
      allPlots["eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(5, "= 0 b-tag");
      allPlots["eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(6, "= 1 b-tag");
      allPlots["eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(7, "#geq 2 b-tag");



      allPlots["raweventflowold_"+systs[i]] = new TH1D("raweventflowold_"+systs[i],";;Events", 7, 0., 7.); 
      allPlots["raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(1, "Init");
      allPlots["raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(2, "1 iso lepton");
      allPlots["raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(3, "#geq 2 jets");
      allPlots["raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(4, "E_{T}^{miss}");
      allPlots["raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(5, "#geq 1 b-tag");
      allPlots["raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(6, "1 #tau");
      allPlots["raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(7, "op. sign");
      
      allPlots["raweventflownocat_"+systs[i]] = new TH1D("raweventflownocat_"+systs[i], ";;Events", 6, 0., 6.);
      allPlots["raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(1, "Init");
      allPlots["raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(2, "1 iso lepton");
      allPlots["raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(3, "#geq 2 jets");
      allPlots["raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(4, "E_{T}^{miss}");
      allPlots["raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(5, "1 #tau");
      allPlots["raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(6, "op. sign");

      allPlots["raweventflowcat_"+systs[i]] = new TH1D("raweventflowcat_"+systs[i], ";;Events", 9, 0., 9.);
      allPlots["raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(1, "Init");
      allPlots["raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(2, "1 iso lepton");
      allPlots["raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(3, "#geq 2 jets");
      allPlots["raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(4, "E_{T}^{miss}");
      allPlots["raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(5, "1 #tau");
      allPlots["raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(6, "op. sign");
      allPlots["raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(7, "= 0 b-tag");
      allPlots["raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(8, "= 1 b-tag");
      allPlots["raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(9, "#geq 2 b-tag");


      
      allPlots["btagcat_"+systs[i]] = new TH1D("btagcat_"+systs[i], ";;Events", 3, 0., 3.);
      allPlots["btagcat_"+systs[i]]->GetXaxis()->SetBinLabel(1, "= 0 b-tags");
      allPlots["btagcat_"+systs[i]]->GetXaxis()->SetBinLabel(2, "= 1 b-tags");
      allPlots["btagcat_"+systs[i]]->GetXaxis()->SetBinLabel(3, "#geq 2 b-tags");

    }

  for (auto& it : allPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  //jet uncertainty parameterization
  TString jecUncUrl("${CMSSW_BASE}/src/UserCode/llvv_fwk/data/jec/25ns/DATA_Uncertainty_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  JetCorrectionUncertainty* jecUnc = new JetCorrectionUncertainty(jecUncUrl.Data());
  
  // setup b-tag calibration readers
  TString btagUncUrl("${CMSSW_BASE}/src/UserCode/llvv_fwk/data/weights/CSVv2.csv");
  gSystem->ExpandPathName(btagUncUrl);
  BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
  BTagCalibrationReader btagSFbReader    (&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central" );
  BTagCalibrationReader btagSFbupReader  (&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "up"      );
  BTagCalibrationReader btagSFbdownReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "down"    ); 
  BTagCalibrationReader btagSFlReader    (&btvcalib, BTagEntry::OP_MEDIUM, "comb"  , "central" );
  BTagCalibrationReader btagSFlupReader  (&btvcalib, BTagEntry::OP_MEDIUM, "comb"  , "up"      );
  BTagCalibrationReader btagSFldownReader(&btvcalib, BTagEntry::OP_MEDIUM, "comb"  , "down"    ); 

  TString btagEffExpUrl("${CMSSW_BASE}/src/UserCode/llvv_fwk/data/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  TFile* beffIn=TFile::Open(btagEffExpUrl);
  TGraphAsymmErrors* expEff_b   =(TGraphAsymmErrors*) beffIn->Get("b"   );
  TGraphAsymmErrors* expEff_c   =(TGraphAsymmErrors*) beffIn->Get("c"   );
  TGraphAsymmErrors* expEff_udsg=(TGraphAsymmErrors*) beffIn->Get("udsg");
  BTagSFUtil myBTagSFUtil;

  //read tree from file
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TTree *t = (TTree*)f->Get("minievents");
  attachToMiniEventTree(t,ev);

  cout << "loop over events" << endl;
  //loop over events
  Int_t nentries(t->GetEntriesFast());
  cout << "...producing " << outname << " from " << nentries << " events" << endl;
  for (Int_t i=0;i<nentries;i++)
    {
      t->GetEntry(i);
      printf ("\r [%3.0f/100] done",100.*(double)(i)/(double)(nentries));

      
      size_t
        nMainLeptons(0),
        nVetoLeptons(0);

      size_t theLep(0);
      double theLepPt(-1.);
      // Select leptons

      std::vector<TLorentzVector> selLeptons; selLeptons.clear();
      
      if(debug) cout << "Event: " << i << endl;
      for (int ilep=0; ilep<ev.nl;++ilep)
	{
          
          bool
            passKin(true),     passId (true),    passIso(true),
            passVetoKin(true), passVetoId(true), passVetoIso(true);

          // MuScle fit corrections already applied. Consider the idea of applying them only now

	  // Check kinematics
	  TLorentzVector ilp4;
	  ilp4.SetPtEtaPhiM(ev.l_pt[ilep],ev.l_eta[ilep],ev.l_phi[ilep],ev.l_mass[ilep]);
          
          //unsigned int lid(fabs(ev.l_id[ilep]));
          double leta(fabs(ilp4.Eta()));

          if(ilp4.Pt() < 30.){ passKin=false; passVetoKin=false;}
          if(leta>2.4)       { passKin=false; passVetoKin=false;}
         
          if(ilp4.Pt() < 20.) passVetoKin=false;
          
          if(debug) cout << "Lepton " << ilep << " has id " << ev.l_tightId[ilep] << ", iso " << ev.l_tightIso[ilep] << endl;

          if(ev.l_tightId[ilep] == 0)  passId=false;
          if(ev.l_tightIso[ilep]== 0) passIso=false;
          // Veto ID and Iso are already applied

          if(passKin && passId && passIso){        
            nMainLeptons++;
            if(ilp4.Pt()>theLepPt) // Redundant: they should be already sorted, so it is OK to take the first one passing the requirements.
              {
                theLep = ilep;
                theLepPt = ilp4.Pt();
              }
            selLeptons.push_back(ilp4);
          }
          else if(passVetoKin && passVetoId && passVetoIso)
            nVetoLeptons++; // "else if" excludes main leptons from veto leptons
	  
        } // End loop on leptons
      
      // For the moment I do not need to sort the collection

      // For the moment, optimized for xsec. One tight lepton, no additional leptons
      if(debug) cout << "Event has " << nMainLeptons << " leptons, from size " << selLeptons.size() << ", and " << nVetoLeptons << " veto leptons." << endl;


      TLorentzVector lp4;
      lp4.SetPtEtaPhiM(ev.l_pt[theLep], ev.l_eta[theLep], ev.l_phi[theLep], ev.l_mass[theLep]);

      // Select according to the lepton id (charge is not needed right now)
      //      if(channelSelection!=0)
      //	{
      //	  if(abs(ev.l_id[theLep])!=abs(channelSelection)) continue;
      //	}
      
      if(debug) cout << "Event has lepton id " << ev.l_id[theLep] << " and flags: muTrigger " << ev.muTrigger << " and elTrigger " << ev.elTrigger << endl;
      //apply trigger requirement
      if(ev.muTrigger && abs(ev.l_id[theLep]) != 13) continue;
      if(ev.elTrigger && abs(ev.l_id[theLep]) != 11) continue;

      if(debug) cout << "\t\t Event passed trigger" << endl;
      
      // FILL HISTFLOW 1

      // Select taus
      size_t theTau(0);
      double theTauPt(-1.);
      std::vector<TLorentzVector> selTaus; selTaus.clear();
      for (int itau=0; itau<ev.nt;++itau)
	{
          // Kin (pre-applied, just a cross check)
          if(ev.t_pt[itau] < 20 || fabs(ev.t_eta[itau]) >2.3) continue;
          
          // Discriminators
          if(ev.t_isodb3hits[itau]      < 0.5) continue;
          if(ev.t_againstMu3Tight[itau] < 0.5) continue; // Or "Loose"
          if(ev.t_againstEl5Medium[itau]< 0.5) continue; // Or "VLoose" or "Loose"

          TLorentzVector itp4;
	  itp4.SetPtEtaPhiM(ev.t_pt[itau],ev.t_eta[itau],ev.t_phi[itau],ev.t_mass[itau]);
          // Lepton cleaning
          for(size_t ilep=0; ilep<selLeptons.size(); ++ilep)
            if(itp4.DeltaR(selLeptons[ilep])<0.4) continue;
          if(itp4.Pt()>theTauPt) // Redundant: they should be already sorted, so it is OK to take the first one passing the requirements.
            {
              theTau = itau;
              theTauPt = itp4.Pt();
            }
          
          selTaus.push_back(itp4);
        }      
      
      
      // Select jets
      uint32_t nJets(0),  nJetsJESLo(0),   nJetsJESHi(0),   nJetsJERLo(0),     nJetsJERHi(0);
      TLorentzVector jetSum(0,0,0,0), jetSumJESup(jetSum), jetSumJESdown(jetSum), jetSumJERup(jetSum), jetSumJERdown(jetSum);
      Double_t htsum(0);
      Int_t nudsgJets(0),ncJets(0), nbJets(0);
      uint32_t nBtags(0), nBtagsBeffLo(0), nBtagsBeffHi(0), nBtagsMistagLo(0), nBtagsMistagHi(0);
      Double_t minMlb(99999.);
      Double_t minMlbJESup(minMlb),  minMlbJESdown(minMlb);
      Double_t minMlbJERup(minMlb),  minMlbJERdown(minMlb);
      Double_t minMlbBeffHi(minMlb), minMlbBeffLo(minMlb);
      Double_t minMlbLeffHi(minMlb), minMlbLeffLo(minMlb);
      std::vector<int> selJetsIdx;           
      for (int ijet=0; ijet<ev.nj;ijet++)
	{
	  //check kinematics
	  TLorentzVector jp4;
	  jp4.SetPtEtaPhiM(ev.j_pt[ijet],ev.j_eta[ijet],ev.j_phi[ijet],ev.j_mass[ijet]);
          // Lepton cleaning
	  if(jp4.DeltaR(lp4)<0.4) continue;
          // Tau cleaning
          bool skipJet(false);
          for(size_t itau=0; itau<selTaus.size(); ++itau)
            if(jp4.DeltaR(selTaus[itau])<0.4){
              skipJet=true;
              break;
            }
          if(skipJet)
            continue;

	  double csv = ev.j_csv[ijet];	  

	  if(fabs(jp4.Eta()) > 2.4) continue;

	  //jet energy scale variations
	  jecUnc->setJetEta(fabs(jp4.Eta()));
	  jecUnc->setJetPt(jp4.Pt());
	  double unc = jecUnc->getUncertainty(true);    
	  
	  //jet energy resolution
	  std::vector<double> jerScale(3,1.0);
	  double genJet_pt(ev.genj_pt[ijet]);
	  if(!ev.isData && genJet_pt>0) jerScale=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt);
	  
	  //readout the b-tagging scale factors for this jet
	  bool isBTagged(csv>0.890),isBTaggedUp(isBTagged),isBTaggedDown(isBTagged);
	  if(!ev.isData)
	    {
              // No random update for now
              bool tempIsBTagged(isBTagged);
              
	      double jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
	      double expEff(1.0), jetBtagSF(1.0), jetBtagSFUp(1.0), jetBtagSFDown(1.0);
	      if(abs(ev.j_hadflav[ijet])==4) 
		{ 
		  expEff=expEff_c->Eval(jptForBtag); 
		  jetBtagSF=btagSFbReader.eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSFUp=btagSFbupReader.eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		  jetBtagSFDown=btagSFbdownReader.eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
		}
	      else if(abs(ev.j_hadflav[ijet])==5) 
		{ 
		  expEff=expEff_b->Eval(jptForBtag); 
		  jetBtagSF=btagSFbReader.eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSFUp=btagSFbupReader.eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		  jetBtagSFDown=btagSFbdownReader.eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
		}
	      else
		{
		  expEff=expEff_udsg->Eval(jptForBtag);
                  jetBtagSF=btagSFlReader.eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
                  jetBtagSFUp=btagSFlupReader.eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
                  jetBtagSFDown=btagSFldownReader.eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
		}
	      
              // No random update, for now (must find updated stuff)
	      myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
	      myBTagSFUtil.modifyBTagsWithSF(isBTaggedUp,  jetBtagSFUp,   expEff);
	      myBTagSFUtil.modifyBTagsWithSF(isBTaggedDown,jetBtagSFDown, expEff);
              isBTagged=tempIsBTagged;
	    }
	  
	  //select the jet
	  if(jp4.Pt()>30)          
	    { 
	      selJetsIdx.push_back(ijet);
	      nJets++;      
	      jetSum += jp4;       
	      htsum += jp4.Pt();
	      if(abs(ev.j_hadflav[ijet])==4)      ncJets++;
	      else if(abs(ev.j_hadflav[ijet])==5) nbJets++;
	      else nudsgJets++;
	      nBtags += isBTagged;
	      if(ev.isData)
		{
		  nBtagsBeffLo += isBTagged;
		  nBtagsBeffHi += isBTagged;
		  nBtagsMistagLo += isBTagged;
		  nBtagsMistagHi += isBTagged;
		}
	      else
		{
		  if(abs(ev.j_hadflav[ijet])==4|| abs(ev.j_hadflav[ijet])==5)
		    {
		      nBtagsBeffLo += isBTaggedDown;
		      nBtagsBeffHi += isBTaggedUp;
		      nBtagsMistagLo += isBTagged;
		      nBtagsMistagHi += isBTagged;
		    }
		  else
		    {
		      nBtagsBeffLo += isBTagged;
		      nBtagsBeffHi += isBTagged;
		      nBtagsMistagLo += isBTaggedDown;
		      nBtagsMistagHi += isBTaggedUp;
		    }
		}

	      if(isBTagged) minMlb = TMath::Min(minMlb,(Double_t)(jp4+lp4).M()); 
	      if(abs(ev.j_hadflav[ijet])==4 || abs(ev.j_hadflav[ijet])==5)
		{
		  if(isBTaggedDown) minMlbBeffLo=TMath::Min(minMlbBeffLo,(Double_t)(jp4+lp4).M());
		  if(isBTaggedUp)   minMlbBeffHi=TMath::Min(minMlbBeffHi,(Double_t)(jp4+lp4).M());
		} 
	      else
		{
		  if(isBTaggedDown) minMlbLeffLo=TMath::Min(minMlbLeffLo,(Double_t)(jp4+lp4).M());
		  if(isBTaggedUp)   minMlbLeffHi=TMath::Min(minMlbLeffHi,(Double_t)(jp4+lp4).M());
		}
	    }
	  if((jp4.Pt())*(1+unc)>30) 
	    { 
	      nJetsJESHi++;
	      jetSumJESup += (1+unc)*jp4;  
	      if(isBTagged) minMlbJESup   = TMath::Min(minMlbJESup,(Double_t)((1+unc)*jp4+lp4).M());
	    }
	  if((jp4.Pt())*(1-unc)>30) 
	    { 
	      nJetsJESLo++; 
	      jetSumJESdown += (1-unc)*jp4;
	      if(isBTagged) minMlbJESdown = TMath::Min(minMlbJESdown,(Double_t)((1-unc)*jp4+lp4).M());
	    } 
	  if(jerScale[1]*jp4.Pt()>30) 
	    { 
	      nJetsJERLo++; 
	      jetSumJERdown += jerScale[1]*jp4;
	      if(isBTagged)  minMlbJERdown = TMath::Min(minMlbJERdown,(Double_t)( jerScale[1]*jp4+lp4).M()); 
	    }
	  if(jerScale[2]*jp4.Pt()>30) 
	    { 
	      nJetsJERHi++; 
	      jetSumJERup += jerScale[2]*jp4;
	      if(isBTagged)  minMlbJERup   = TMath::Min(minMlbJERup,(Double_t)( jerScale[2]*jp4+lp4).M());
	    }
	}
      

      if(debug) cout << "\t\t event passed jets loop " << endl; 
    

      //varied MET
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);
      TLorentzVector metJESup( met+(jetSum-jetSumJESup)),metJESdown(met+(jetSum-jetSumJESdown));
      TLorentzVector metJERup( met+(jetSum-jetSumJERup)),metJERdown(met+(jetSum-jetSumJERdown));
      TLorentzVector metUMetdown(0.9*met-0.1*(jetSum+lp4)),metUMetup(1.1*met+0.1*(jetSum+lp4));
      
      //varied MT
      double mt( computeMT(lp4,met) );
      double mtJESup( computeMT(lp4,metJESup) ), mtJESdown( computeMT(lp4,metJESdown) );
      double mtJERup( computeMT(lp4,metJESup) ), mtJERdown( computeMT(lp4,metJERdown) );
      double mtUMetup( computeMT(lp4,metUMetup) ), mtUMetdown( computeMT(lp4,metUMetdown) );
      
      //check if flavour splitting was required
      if(flavourSplitting!=FlavourSplitting::NOFLAVOURSPLITTING)
	{
	  if(flavourSplitting==FlavourSplitting::BSPLITTING)         { if(nbJets==0)    continue; }
	  else if(flavourSplitting==FlavourSplitting::CSPLITTING)    { if(ncJets==0 || nbJets!=0)    continue; }
	  else if(flavourSplitting==FlavourSplitting::UDSGSPLITTING) { if(nudsgJets==0 || ncJets!=0 || nbJets!=0) continue; }
	}

      if(debug) cout << "\t\t event passed MET and MT variations " << endl;

      //generator level weights to apply
      std::vector<double> lepSF=getLeptonSelectionScaleFactor(ev.l_id[theLep],ev.l_pt[theLep], ev.l_eta[theLep],ev.isData);
      std::vector<double> puWeight(3,1.0);
      if(debug) cout << "\t\t lepton scale factors acquired " << endl;
      if(!ev.isData && puWgtGr)
	{
	  puWeight[0]=puWgtGr->Eval(ev.putrue);  puWeight[1]=puUpWgtGr->Eval(ev.putrue); puWeight[2]=puDownWgtGr->Eval(ev.putrue);
	}
      if(debug) cout << "\t\t pileup acquired" << endl;

      double norm=normH ? normH->GetBinContent(1) : 1.0;
      double wgt          (norm*lepSF[0]                                *puWeight[0]);
      if(debug) cout << "\t\t wgt = norm*lepSF*puweight[0] , i.e. " << wgt << "=" << norm << "*" << lepSF[0] << "*" <<puWeight[0] << endl;
      double wgtPuUp      (norm*lepSF[0]                                *puWeight[1]);
      double wgtPuDown    (norm*lepSF[0]                                *puWeight[2]);
      double wgtMuEffUp   (norm*(abs(ev.l_id[theLep])==13 ? lepSF[1] : lepSF[0])*puWeight[0]);
      double wgtMuEffDown (norm*(abs(ev.l_id[theLep])==13 ? lepSF[2] : lepSF[0])*puWeight[0]);
      double wgtElEffUp   (norm*(abs(ev.l_id[theLep])==11 ? lepSF[1] : lepSF[0])*puWeight[0]);
      double wgtElEffDown (norm*(abs(ev.l_id[theLep])==11 ? lepSF[2] : lepSF[0])*puWeight[0]);
      double wgtQCDScaleLo(wgt),wgtQCDScaleHi(wgt);
      if(debug) cout << "\t\t only ttbar weight is missing " << endl;
      if(genWgtMode!=NOGENWGT && !ev.isData) 
	{
	  wgt           *= ev.ttbar_w[0];
	  wgtPuUp       *= ev.ttbar_w[0];
	  wgtPuDown     *= ev.ttbar_w[0];
	  wgtMuEffUp    *= ev.ttbar_w[0];
	  wgtMuEffDown  *= ev.ttbar_w[0];
	  wgtElEffUp    *= ev.ttbar_w[0];
	  wgtElEffDown  *= ev.ttbar_w[0];
	  wgtQCDScaleLo *= ev.ttbar_w[0];
	  wgtQCDScaleHi *= ev.ttbar_w[0];
	  if(debug) cout << "\t\t weight has been multiplied by ev.ttbar_w[0]=" << ev.ttbar_w[0] << ", yielding wgt=" << wgt << endl;
	}
      if(debug) cout << "\t\t acquiring normalization from normH histogram, which has pointer: " << normH << endl;
      if(isTTbar)
	{
	  wgtQCDScaleLo   = wgt*(normH->GetBinContent(10)/norm)*(ev.ttbar_w[9]/ev.ttbar_w[0]);
	  wgtQCDScaleHi   = wgt*(normH->GetBinContent(6)/norm)*(ev.ttbar_w[5]/ev.ttbar_w[0]);	 
	}

      if(debug) cout << "\t\t event passed genevent stuff" << endl;


      allPlots["eventflowold_nom"]   ->Fill(0., wgt);
      allPlots["eventflownocat_nom"] ->Fill(0., wgt);
      allPlots["eventflowcat_nom"]   ->Fill(0., wgt);
      
      if(nMainLeptons!=1 || nVetoLeptons>0) continue;
      


      //nominal selection
      if(nJets>=1)
	{
	  int nJetsCat=TMath::Min((int)nJets,(int)4);
	  int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  
	  allPlots["njetsnbtags_nom"]->Fill(binToFill,wgt);
	  allPlots["njetsnbtags_qcdScaleDown"]->Fill(binToFill,wgtQCDScaleLo);
	  allPlots["njetsnbtags_qcdScaleUp"]->Fill(binToFill,wgtQCDScaleHi);
	  allPlots["njetsnbtags_puUp"]->Fill(binToFill,wgtPuUp);
	  allPlots["njetsnbtags_puDown"]->Fill(binToFill,wgtPuDown);
	  allPlots["njetsnbtags_muEffUp"]->Fill(binToFill,wgtMuEffUp);
	  allPlots["njetsnbtags_muEffDown"]->Fill(binToFill,wgtMuEffDown);
	  allPlots["njetsnbtags_umetUp"]->Fill(binToFill,wgt);
	  allPlots["njetsnbtags_umetDown"]->Fill(binToFill,wgt);
	  
	  std::vector<TString> catsToFill(2,Form("%dj",nJetsCat));
	  catsToFill[1]+= Form("%dt",nBtagsCat);
	  
	  std::map<Int_t,Double_t>::iterator rIt=lumiMap.find(ev.run);
	  for(size_t icat=0; icat<2; icat++)
	    {
	      TString tag=catsToFill[icat];
	      if(rIt!=lumiMap.end()) {
		Int_t runCtr=std::distance(lumiMap.begin(),rIt);
		allPlots["ratevsrun_"+tag]->Fill(runCtr,1.e+6/rIt->second);
	      }
	      allPlots["lpt_"+tag]->Fill(ev.l_pt[theLep],wgt);
	      allPlots["leta_"+tag]->Fill(ev.l_eta[theLep],wgt);
	      allPlots["jpt_"+tag]->Fill(ev.j_pt[ selJetsIdx[0] ],wgt);
	      allPlots["jeta_"+tag]->Fill(fabs(ev.j_eta[ selJetsIdx[0] ]),wgt);
	      allPlots["csv_"+tag]->Fill(ev.j_csv[ selJetsIdx[0] ],wgt);
	      allPlots["ht_"+tag]->Fill(htsum,wgt);
	      allPlots["nvtx_"+tag]->Fill(ev.nvtx,wgt);
	      allPlots["met_"+tag]->Fill(ev.met_pt,wgt);
	      allPlots["metphi_"+tag]->Fill(ev.met_phi,wgt);
	      allPlots["mt_nom_"+tag]->Fill(mt,wgt);
	      allPlots["mt_qcdScaleDown_"+tag]->Fill(mt,wgtQCDScaleLo);
	      allPlots["mt_qcdScaleUp_"+tag]->Fill(mt,wgtQCDScaleHi);
	      allPlots["mt_puUp_"+tag]->Fill(mt,wgtPuUp);
	      allPlots["mt_puDown_"+tag]->Fill(mt,wgtPuDown);
	      allPlots["mt_muEffUp_"+tag]->Fill(mt,wgtMuEffUp);
	      allPlots["mt_muEffDown_"+tag]->Fill(mt,wgtMuEffDown);
	      allPlots["mt_eEffUp_"+tag]->Fill(mt,wgtElEffUp);
	      allPlots["mt_eEffDown_"+tag]->Fill(mt,wgtElEffDown);	      
	      allPlots["mt_umetUp_"+tag]->Fill(mtUMetup,wgtElEffUp);
	      allPlots["mt_umetDown_"+tag]->Fill(mtUMetdown,wgtElEffDown);	      

	      if(nBtagsCat>0)
		{
		  allPlots["minmlb_nom_"+tag]->Fill(minMlb,wgt);
		  allPlots["minmlb_qcdScaleDown_"+tag]->Fill(minMlb,wgtQCDScaleLo);
		  allPlots["minmlb_qcdScaleUp_"+tag]->Fill(minMlb,wgtQCDScaleHi);
		  allPlots["minmlb_puUp_"+tag]->Fill(minMlb,wgtPuUp);
		  allPlots["minmlb_puDown_"+tag]->Fill(minMlb,wgtPuDown);
		  allPlots["minmlb_muEffUp_"+tag]->Fill(minMlb,wgtMuEffUp);
		  allPlots["minmlb_muEffDown_"+tag]->Fill(minMlb,wgtMuEffDown);
		  allPlots["minmlb_eEffUp_"+tag]->Fill(minMlb,wgtElEffUp);
		  allPlots["minmlb_eEffDown_"+tag]->Fill(minMlb,wgtElEffDown);
		  allPlots["minmlb_umetUp_"+tag]->Fill(minMlb,wgt);
		  allPlots["minmlb_umetDown_"+tag]->Fill(minMlb,wgt);
		}	  
	    }
	}

      if(nJetsJESHi>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJESHi,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_jesUp"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots["mt_jesUp_"+tag]->Fill(mtJESup,wgt);
	  allPlots["minmlb_jesUp_"+tag]->Fill(minMlbJESup,wgt);
	}
      if(nJetsJESLo>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJESLo,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_jesDown"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots["mt_jesDown_"+tag]->Fill(mtJESdown,wgt);
	  allPlots["minmlb_jesDown_"+tag]->Fill(minMlbJESdown,wgt);
	}
      if(nJetsJERHi>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJERHi,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_jerUp"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots["mt_jerUp_"+tag]->Fill(mtJERup,wgt);
	  allPlots["minmlb_jerUp_"+tag]->Fill(minMlbJERup,wgt);
	}
      if(nJetsJERLo>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJERLo,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_jerDown"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots["mt_jerDown_"+tag]->Fill(mtJERdown,wgt);
	  allPlots["minmlb_jerDown_"+tag]->Fill(minMlbJERdown,wgt);
	}
      if(nJets>=1)
	{
	  int nJetsCat=TMath::Min((int)nJets,(int)4);
          
	  int nBtagsCat=TMath::Min((int)nBtagsBeffLo,(int)2);
          int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots["njetsnbtags_beffDown"]->Fill(binToFill,wgt); 
	  
	  nBtagsCat=TMath::Min((int)nBtagsBeffHi,(int)2);
	  binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots["njetsnbtags_beffUp"]->Fill(binToFill,wgt); 

	  nBtagsCat=TMath::Min((int)nBtagsMistagLo,(int)2);
	  binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots["njetsnbtags_mistagDown"]->Fill(binToFill,wgt); 

	  nBtagsCat=TMath::Min((int)nBtagsMistagHi,(int)2);
	  binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots["njetsnbtags_mistagUp"]->Fill(binToFill,wgt); 
	}

      if(debug) cout << "\t\t Now nominal selection starts" << endl;
      // Nominal selection

      // One lepton requirement already applied

      bool
        pass2Jets   (nJets>1),
        passMet     (met.Pt()>40.),
        pass1btag   (nBtags>0),
        pass1tau    (selTaus.size()>0),
        passOS      (selTaus.size()>0 ? (ev.t_charge[theTau] * ev.l_charge[theLep] < 0): false  ),
        passCat0btag(nBtags==0),
        passCat1btag(nBtags==1),
        passCat2btag(nBtags>1);

      if(debug) cout << "\t\t Event categorization: pass2Jets " << pass2Jets << ", passMet " << passMet << ", pass1btag " << pass1btag << ", pass1tau " << pass1tau << ", passOS " << passOS << ", passCat0btag " << passCat0btag << ", passCat1btag " << passCat1btag << ", passCat2btag " << passCat2btag << endl;

      if(debug) cout << "\t\t Event weight: " << wgt << endl;

      
                                                                    allPlots["raweventflowold_nom"]->Fill(1., wgt);
      if(pass2Jets)                                               { allPlots["raweventflowold_nom"]->Fill(2., wgt); allPlots["eventflowold_nom"]->Fill(0., wgt); }
      if(pass2Jets && passMet)                                    { allPlots["raweventflowold_nom"]->Fill(3., wgt); allPlots["eventflowold_nom"]->Fill(1., wgt); }
      if(pass2Jets && passMet && pass1btag)                       { allPlots["raweventflowold_nom"]->Fill(4., wgt); allPlots["eventflowold_nom"]->Fill(2., wgt); }
      if(pass2Jets && passMet && pass1btag && pass1tau)           { allPlots["raweventflowold_nom"]->Fill(5., wgt); allPlots["eventflowold_nom"]->Fill(3., wgt); }
      if(pass2Jets && passMet && pass1btag && pass1tau && passOS) { allPlots["raweventflowold_nom"]->Fill(6., wgt); allPlots["eventflowold_nom"]->Fill(4., wgt); }
      
                                                       allPlots["raweventflownocat_nom"]->Fill(1., wgt);
      if(pass2Jets)                                  { allPlots["raweventflownocat_nom"]->Fill(2., wgt); allPlots["eventflownocat_nom"]->Fill(0., wgt); }
      if(pass2Jets && passMet)                       { allPlots["raweventflownocat_nom"]->Fill(3., wgt); allPlots["eventflownocat_nom"]->Fill(1., wgt); }
      if(pass2Jets && passMet && pass1tau)           { allPlots["raweventflownocat_nom"]->Fill(4., wgt); allPlots["eventflownocat_nom"]->Fill(2., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS) { allPlots["raweventflownocat_nom"]->Fill(5., wgt); allPlots["eventflownocat_nom"]->Fill(3., wgt); }

                                                                       allPlots["raweventflowcat_nom"]->Fill(1., wgt);
      if(pass2Jets)                                                  { allPlots["raweventflowcat_nom"]->Fill(2., wgt); allPlots["eventflowcat_nom"]->Fill(0., wgt); }
      if(pass2Jets && passMet)                                       { allPlots["raweventflowcat_nom"]->Fill(3., wgt); allPlots["eventflowcat_nom"]->Fill(1., wgt); }
      if(pass2Jets && passMet && pass1tau)                           { allPlots["raweventflowcat_nom"]->Fill(4., wgt); allPlots["eventflowcat_nom"]->Fill(2., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS)                 { allPlots["raweventflowcat_nom"]->Fill(5., wgt); allPlots["eventflowcat_nom"]->Fill(3., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && passCat0btag) { allPlots["raweventflowcat_nom"]->Fill(6., wgt); allPlots["eventflowcat_nom"]->Fill(4., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && passCat1btag) { allPlots["raweventflowcat_nom"]->Fill(7., wgt); allPlots["eventflowcat_nom"]->Fill(5., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && passCat2btag) { allPlots["raweventflowcat_nom"]->Fill(8., wgt); allPlots["eventflowcat_nom"]->Fill(6., wgt); }
     
      if(pass2Jets && passMet && pass1tau && passOS && passCat0btag) allPlots["btagcat_nom"]->Fill(0., wgt);
      if(pass2Jets && passMet && pass1tau && passOS && passCat1btag) allPlots["btagcat_nom"]->Fill(1., wgt);
      if(pass2Jets && passMet && pass1tau && passOS && passCat2btag) allPlots["btagcat_nom"]->Fill(2., wgt);

      if(debug) cout << "Event analyzed fully " << endl;

    }
  
  //close input file
  f->Close();

  //save histos to file  
  TString selPrefix("");  
  if(flavourSplitting!=NOFLAVOURSPLITTING) selPrefix=Form("%d_",flavourSplitting);
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+selPrefix+baseName,"RECREATE");
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

//
std::map<Int_t,Double_t> lumiPerRun()
{
  std::map<Int_t,Double_t> toReturn;
  toReturn[ 256630 ]=  948417.609 ;
  toReturn[ 256673 ]=   5534.408  ;
  toReturn[ 256674 ]=  92567.485  ;
  toReturn[ 256675 ]= 7099162.193 ;
  toReturn[ 256676 ]= 9172872.886 ;
  toReturn[ 256677 ]= 15581756.928;
  toReturn[ 256801 ]= 8830347.675 ;
  toReturn[ 256842 ]=  16510.582  ;
  toReturn[ 256843 ]= 37131085.338;
  toReturn[ 256866 ]=  58250.406  ;
  toReturn[ 256867 ]= 4546508.033 ;
  toReturn[ 256868 ]= 22542014.201;
  toReturn[ 256869 ]= 1539580.832 ;
  toReturn[ 256926 ]= 1499855.808 ;
  toReturn[ 256941 ]= 8741029.381 ;
  toReturn[ 257461 ]= 3057928.782 ;
  toReturn[ 257531 ]= 8418598.194 ;
  toReturn[ 257599 ]= 4940876.751 ;
  toReturn[ 257613 ]= 75519819.209;
  toReturn[ 257614 ]=  850778.922 ;
  toReturn[ 257645 ]= 62388624.946;
  toReturn[ 257682 ]= 13053256.987;
  toReturn[ 257722 ]=  719350.314 ;
  toReturn[ 257723 ]= 5941442.106 ;
  toReturn[ 257735 ]=  521278.124 ;
  toReturn[ 257751 ]= 27029514.967;
  toReturn[ 257804 ]=  210956.374 ;
  toReturn[ 257805 ]= 17038078.687;
  toReturn[ 257816 ]= 24328019.178;
  toReturn[ 257819 ]= 15147148.510;
  toReturn[ 257968 ]= 16769109.914;
  toReturn[ 257969 ]= 39179793.996;
  toReturn[ 258129 ]= 5813530.480 ;
  toReturn[ 258136 ]= 3617731.160 ;
  toReturn[ 258157 ]= 3866329.715 ;
  toReturn[ 258158 ]=105571609.093;
  toReturn[ 258159 ]= 25531210.007;
  toReturn[ 258177 ]=101938042.657;
  toReturn[ 258211 ]= 6371404.543 ;
  toReturn[ 258213 ]= 11238447.671;
  toReturn[ 258214 ]= 15172855.551;
  toReturn[ 258215 ]=  411364.505 ;
  toReturn[ 258287 ]= 12966493.641;
  toReturn[ 258403 ]= 15612888.766;
  toReturn[ 258425 ]= 10170405.992;
  toReturn[ 258426 ]=  751812.067 ;
  toReturn[ 258427 ]= 7901302.746 ;
  toReturn[ 258428 ]= 11578208.046;
  toReturn[ 258432 ]=  279944.749 ;
  toReturn[ 258434 ]= 30536738.787;
  toReturn[ 258440 ]= 44624917.855;
  toReturn[ 258444 ]= 2079367.112 ;
  toReturn[ 258445 ]= 16403375.902;
  toReturn[ 258446 ]= 7474593.995 ;
  toReturn[ 258448 ]= 35705866.510;
  toReturn[ 258655 ]=  383658.834 ;
  toReturn[ 258656 ]= 25581933.798;
  toReturn[ 258694 ]= 15219679.888;
  toReturn[ 258702 ]= 29093248.654;
  toReturn[ 258703 ]= 31138680.065;
  toReturn[ 258705 ]= 7604760.367 ;
  toReturn[ 258706 ]= 52122692.407;
  toReturn[ 258712 ]= 34495799.123;
  toReturn[ 258713 ]= 10164347.291;
  toReturn[ 258714 ]= 4168945.356 ;
  toReturn[ 258741 ]= 4446752.908 ;
  toReturn[ 258742 ]= 59299810.293;
  toReturn[ 258745 ]= 20966777.757;
  toReturn[ 258749 ]= 44752319.500;
  toReturn[ 258750 ]= 14330984.460;
  return toReturn;
};

//
std::vector<double> getLeptonSelectionScaleFactor(int id, double pt, double eta, bool isData)
{
  std::vector<double> lepSelSF(3,1.0);
  if(isData) return lepSelSF;

  std::pair<double,double>res(1.0,0.0);

  //electrons
  if(abs(id)==11)
    {
      if (fabs(eta)<0.8)
	{
	  if (pt<30)      { res.first=0.927; res.second=0.073; }
	  else if (pt<40) { res.first=0.975; res.second=0.018; }
	  else if (pt<50) { res.first=0.962; res.second=0.036; }
	  else            { res.first=0.955; res.second=0.022; }
	}
      else if (fabs(eta)<1.5)
	{
	  if (pt<30)      { res.first=0.891; res.second=0.074; }
	  else if (pt<40) { res.first=0.965; res.second=0.020; }
	  else if (pt<50) { res.first=0.968; res.second=0.018; }
	  else            { res.first=0.955; res.second=0.018; }
	}
      else
	{
	  if (pt<30)      { res.first=0.956; res.second=0.059; }
	  else if (pt<40) { res.first=0.995; res.second=0.018; }
	  else if (pt<50) { res.first=0.993; res.second=0.019; }
	  else            { res.first=0.985; res.second=0.023; }
	}
    }

  //muons
  if (abs(id)==13)
    {
      if (fabs(eta)<0.9)
	{
	  if (pt<30)      { res.first=1.003; res.second=0.019; }
	  else if (pt<40) { res.first=1.014; res.second=0.015; }
	  else if (pt<50) { res.first=1.001; res.second=0.014; }
	  else            { res.first=0.983; res.second=0.014; }
	}
      else if(fabs(eta)<1.2)
	{
	  if (pt<30)      { res.first=0.993; res.second=0.019; }
	  else if (pt<40) { res.first=0.994; res.second=0.015; }
	  else if (pt<50) { res.first=0.980; res.second=0.014; }
	  else            { res.first=0.987; res.second=0.015; }
	}
      else
	{
	  if (pt<30)      { res.first=1.023; res.second=0.028; }
	  else if (pt<40) { res.first=0.994; res.second=0.014; }
	  else if (pt<50) { res.first=0.996; res.second=0.014; }
	  else            { res.first=0.979; res.second=0.014; }
	}
    }

  lepSelSF[0]=res.first;
  lepSelSF[1]=res.first+res.second;
  lepSelSF[2]=res.first-res.second;
  return lepSelSF;
}

//Sources
//  Assuming nominal JER but uncertainties from Run I
//  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
std::vector<double> getJetResolutionScales(double pt, double eta, double genjpt)
{
  std::vector<double> res(3,1.0);

  double ptSF(1.0), ptSF_err(0.06);
  if(TMath::Abs(eta)<0.5) 
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2));
    }
  else if(TMath::Abs(eta)<1.1)
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2));
    }
  else if(TMath::Abs(eta)<1.7)
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2));
    }
  else if(TMath::Abs(eta)<2.3)
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    }
  else
    {
      ptSF=1.0;
      ptSF_err = TMath::Sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2));
    }

  res[0] = TMath::Max((Double_t)0.,(Double_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;
  res[1] = TMath::Max((Double_t)0.,(Double_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;
  res[2] = TMath::Max((Double_t)0.,(Double_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;
  
  return res;
}
