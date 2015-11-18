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


//void FillCutPlots(std::map<TString, TH1*> allPlots, TString& chan, TString cut, TString syst, MiniEvent_t& ev, size_t& theLep, std::vector<int> selJetsIdx, Double_t& wgt)
//{
//
//}


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

  std::vector<TString> chanList; chanList.clear();
  chanList.push_back("mu_");
  chanList.push_back("el_");

  std::map<Int_t,Double_t> lumiMap=lumiPerRun();

  for(std::vector<TString>::iterator ichan=chanList.begin(); ichan!=chanList.end(); ++ichan) // Loop on channels
    {
      TString chan(*ichan);
      for(Int_t ij=1; ij<=4; ij++)
        {
          for(Int_t itag=-1; itag<=2; itag++)
            {
              if(itag>ij) continue;
              TString tag(itag<0 ? Form("%dj",ij) : Form("%dj%dt",ij,itag));
              allPlots[chan+"ratevsrun_"+tag] = new TH1D(chan+"ratevsrun_"+tag,";Run number; Events/pb",lumiMap.size(),0,lumiMap.size());
              Int_t runCtr(0);
              for(std::map<Int_t,Double_t>::iterator it=lumiMap.begin(); it!=lumiMap.end(); it++,runCtr++)
                allPlots[chan+"ratevsrun_"+tag]->GetXaxis()->SetBinLabel(runCtr+1,Form("%d",it->first));
              allPlots[chan+"lpt_"+tag]  = new TH1D(chan+"lpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
              allPlots[chan+"lsip3d_"+tag]  = new TH1D(chan+"lsip3d_"+tag,";3d impact parameter significance;Events" ,40,0.,20.);
              allPlots[chan+"lchiso_"+tag]  = new TH1D(chan+"lchiso_"+tag,";Charged hadron isolation [GeV];Events" ,25,0.,50.);
              allPlots[chan+"lchreliso_"+tag]  = new TH1D(chan+"lchreliso_"+tag,";Charged hadron relative isolation;Events" ,25,0.,0.2);
              allPlots[chan+"leta_"+tag] = new TH1D(chan+"leta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
              allPlots[chan+"jpt_"+tag]  = new TH1D(chan+"jpt_"+tag,";Transverse momentum [GeV];Events" ,20,0.,300.);
              allPlots[chan+"jeta_"+tag] = new TH1D(chan+"jeta_"+tag,";Pseudo-rapidity;Events" ,12,0.,3.);
              allPlots[chan+"ht_"+tag]   = new TH1D(chan+"ht_"+tag,";H_{T} [GeV];Events",40,0,800);
              allPlots[chan+"csv_"+tag]  = new TH1D(chan+"csv_"+tag,";CSV discriminator;Events",100,0,1.0);
              allPlots[chan+"nvtx_"+tag] = new TH1D(chan+"nvtx_"+tag,";Vertex multiplicity;Events" ,50,0.,50.);
              allPlots[chan+"met_"+tag]  = new TH1D(chan+"metpt_"+tag,";Missing transverse energy [GeV];Events" ,20,0.,300.);
              allPlots[chan+"metphi_"+tag] = new TH1D(chan+"metphi_" + tag,";MET #phi [rad];Events" ,50,-3.2,3.2);
              for(size_t isyst=0; isyst<sizeof(systs)/sizeof(TString); isyst++)
                {
                  allPlots[chan+"mt_"+systs[isyst]+"_"+tag]     = new TH1D(chan+"mt_"+systs[isyst]+"_"+tag,";Transverse Mass [GeV];Events" ,20,0.,200.);
                  allPlots[chan+"minmlb_"+systs[isyst]+"_"+tag] = new TH1D(chan+"minmlb_"+systs[isyst]+"_"+tag,";min Mass(lepton,b) [GeV];Events" ,25,0.,250.);
                }
            }
        }
      
      //category counting
      for(size_t i=0; i<sizeof(systs)/sizeof(TString); i++)
        {
          allPlots[chan+"njetsnbtags_"+systs[i]] = new TH1D(chan+"njetsnbtags_"+systs[i],";Category;Events" ,12,0.,12.);
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(1, "1j,=0b");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(2, "1j,=1b");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(3, "");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(4, "2j,=0b");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(5, "2j,=1b");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(6, "2j,#geq2b");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(7, "3j,=0b");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(8, "3j,=1b");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(9, "3j,#geq2b");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(10,"4j,=0b");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(11,"4j,=1b");
          allPlots[chan+"njetsnbtags_"+systs[i]]->GetXaxis()->SetBinLabel(12,"4j,#geq2b");
          
          
          allPlots[chan+"eventflowold_"+systs[i]] = new TH1D(chan+"eventflowold_"+systs[i],";;Events", 5, 0., 5.); 
          allPlots[chan+"eventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(1, "1 #it{l}, #geq 2 jets");
          allPlots[chan+"eventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(2, "E_{T}^{miss}");
          allPlots[chan+"eventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(3, "1 #tau");
          allPlots[chan+"eventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(4, "op. sign");
          allPlots[chan+"eventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(5, "#geq 1 b-tag");
          
          allPlots[chan+"eventflownocat_"+systs[i]] = new TH1D(chan+"eventflownocat_"+systs[i], ";;Events", 4, 0., 4.);
          allPlots[chan+"eventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(1, "1 #it{l}, #geq 2 jets");
          allPlots[chan+"eventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(2, "E_{T}^{miss}");
          allPlots[chan+"eventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(3, "1 #tau");
          allPlots[chan+"eventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(4, "op. sign");
          
          allPlots[chan+"eventflowcat_"+systs[i]] = new TH1D(chan+"eventflowcat_"+systs[i], ";;Events", 7, 0., 7.);
          allPlots[chan+"eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(1, "1 #it{l}, #geq 2 jets");
          allPlots[chan+"eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(2, "E_{T}^{miss}");
          allPlots[chan+"eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(3, "1 #tau");
          allPlots[chan+"eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(4, "op. sign");
          allPlots[chan+"eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(5, "= 0 b-tag");
          allPlots[chan+"eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(6, "= 1 b-tag");
          allPlots[chan+"eventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(7, "#geq 2 b-tag");
          
          
          
          allPlots[chan+"raweventflowold_"+systs[i]] = new TH1D(chan+"raweventflowold_"+systs[i],";;Events", 7, 0., 7.); 
          allPlots[chan+"raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(1, "Init");
          allPlots[chan+"raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(2, "1 iso lepton");
          allPlots[chan+"raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(3, "#geq 2 jets");
          allPlots[chan+"raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(4, "E_{T}^{miss}");
          allPlots[chan+"raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(5, "1 #tau");
          allPlots[chan+"raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(6, "op. sign");
          allPlots[chan+"raweventflowold_"+systs[i]]->GetXaxis()->SetBinLabel(7, "#geq 1 b-tag");
          
          allPlots[chan+"raweventflownocat_"+systs[i]] = new TH1D(chan+"raweventflownocat_"+systs[i], ";;Events", 6, 0., 6.);
          allPlots[chan+"raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(1, "Init");
          allPlots[chan+"raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(2, "1 iso lepton");
          allPlots[chan+"raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(3, "#geq 2 jets");
          allPlots[chan+"raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(4, "E_{T}^{miss}");
          allPlots[chan+"raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(5, "1 #tau");
          allPlots[chan+"raweventflownocat_"+systs[i]]->GetXaxis()->SetBinLabel(6, "op. sign");
          
          allPlots[chan+"raweventflowcat_"+systs[i]] = new TH1D(chan+"raweventflowcat_"+systs[i], ";;Events", 9, 0., 9.);
          allPlots[chan+"raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(1, "Init");
          allPlots[chan+"raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(2, "1 iso lepton");
          allPlots[chan+"raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(3, "#geq 2 jets");
          allPlots[chan+"raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(4, "E_{T}^{miss}");
          allPlots[chan+"raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(5, "1 #tau");
          allPlots[chan+"raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(6, "op. sign");
          allPlots[chan+"raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(7, "= 0 b-tag");
          allPlots[chan+"raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(8, "= 1 b-tag");
          allPlots[chan+"raweventflowcat_"+systs[i]]->GetXaxis()->SetBinLabel(9, "#geq 2 b-tag");
          
          
          
          allPlots[chan+"btagcat_"+systs[i]] = new TH1D(chan+"btagcat_"+systs[i], ";;Events", 3, 0., 3.);
          allPlots[chan+"btagcat_"+systs[i]]->GetXaxis()->SetBinLabel(1, "= 0 b-tags");
          allPlots[chan+"btagcat_"+systs[i]]->GetXaxis()->SetBinLabel(2, "= 1 b-tags");
          allPlots[chan+"btagcat_"+systs[i]]->GetXaxis()->SetBinLabel(3, "#geq 2 b-tags");

          TString cuts[] = {
            "1l2j_",
            "met_",
            "1tau_",
            "os_",
            "1bi_",
            "0b_",
            "1b_",
            "2b_"
          };
          
          for(size_t icut=0; icut<sizeof(cuts)/sizeof(TString); ++icut)
            {
              TString cut(cuts[icut]);
              if(systs[i] != "nom") continue;
              allPlots[chan+cut+"lpt_"+systs[i]]       = new TH1D(chan+cut+"lpt_"+systs[i],";Transverse momentum [GeV];Events" ,50,0.,500.);
              allPlots[chan+cut+"leta_"+systs[i]]      = new TH1D(chan+cut+"leta_"+systs[i],";Pseudo-rapidity;Events" ,26,-2.6,2.6);
              allPlots[chan+cut+"jpt_"+systs[i]]       = new TH1D(chan+cut+"jpt_"+systs[i],";Transverse momentum [GeV];Events" ,50,0.,500.);
              allPlots[chan+cut+"jeta_"+systs[i]]      = new TH1D(chan+cut+"jeta_"+systs[i],";Pseudo-rapidity;Events" ,26,-2.6,2.6);
              allPlots[chan+cut+"tpt_"+systs[i]]       = new TH1D(chan+cut+"tpt_"+systs[i],";Transverse momentum [GeV];Events" ,30,0.,300.);
              allPlots[chan+cut+"teta_"+systs[i]]      = new TH1D(chan+cut+"teta_"+systs[i],";Pseudo-rapidity;Events" ,26,-2.6,2.6);
              //allPlots[chan+cut+"ht_"+systs[i]]        = new TH1D(chan+cut+"ht_"+systs[i],";H_{T} [GeV];Events",40,0,800);
              allPlots[chan+cut+"csv_"+systs[i]]       = new TH1D(chan+cut+"csv_"+systs[i],";CSV discriminator;Events",100,0,1.0);
              allPlots[chan+cut+"nvtx_"+systs[i]]      = new TH1D(chan+cut+"nvtx_"+systs[i],";Vertex multiplicity;Events" ,100,0.,100.);
              allPlots[chan+cut+"nvtxraw_"+systs[i]]   = new TH1D(chan+cut+"nvtxraw_"+systs[i],";Vertex multiplicity;Events" ,100,0.,100.);
              allPlots[chan+cut+"met_"+systs[i]]       = new TH1D(chan+cut+"metpt_"+systs[i],";Missing transverse energy [GeV];Events" ,50,0.,1000.);
              allPlots[chan+cut+"metphi_"+systs[i]]    = new TH1D(chan+cut+"metphi_" + systs[i],";MET #phi [rad];Events" ,64,-3.2,3.2);
              allPlots[chan+cut+"nbjets_"+systs[i]]    = new TH1D(chan+cut+"nbjets_" + systs[i],";b-tags multiplicity;Events",5, 0.,5.);
              
              

              allPlots[chan+cut+"nbjets_"+systs[i]]->GetXaxis()->SetBinLabel(1, "=0 b-tag");
              allPlots[chan+cut+"nbjets_"+systs[i]]->GetXaxis()->SetBinLabel(2, "=1 b-tag");
              allPlots[chan+cut+"nbjets_"+systs[i]]->GetXaxis()->SetBinLabel(3, "=2 b-tag");
              allPlots[chan+cut+"nbjets_"+systs[i]]->GetXaxis()->SetBinLabel(4, "=3 b-tag");
              allPlots[chan+cut+"nbjets_"+systs[i]]->GetXaxis()->SetBinLabel(5, "#geq4 b-tag");

            }

          
        }
      
    } // End loop on channels

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

      TString channel("");
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

      if(abs(ev.l_id[theLep]) == 13) channel="mu_";
      if(abs(ev.l_id[theLep]) == 11) channel="el_";
      if(channel!="mu_" && channel!="el_") cout << "CHANNEL IS " << channel << endl; 

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


      allPlots[channel+"raweventflowold_nom"]   ->Fill(0., wgt);
      allPlots[channel+"raweventflownocat_nom"] ->Fill(0., wgt);
      allPlots[channel+"raweventflowcat_nom"]   ->Fill(0., wgt);
      
      if(nMainLeptons!=1 || nVetoLeptons>0) continue;
      


      //nominal selection
      if(nJets>=1)
	{
	  int nJetsCat=TMath::Min((int)nJets,(int)4);
	  int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  
	  allPlots[channel+"njetsnbtags_nom"]->Fill(binToFill,wgt);
	  allPlots[channel+"njetsnbtags_qcdScaleDown"]->Fill(binToFill,wgtQCDScaleLo);
	  allPlots[channel+"njetsnbtags_qcdScaleUp"]->Fill(binToFill,wgtQCDScaleHi);
	  allPlots[channel+"njetsnbtags_puUp"]->Fill(binToFill,wgtPuUp);
	  allPlots[channel+"njetsnbtags_puDown"]->Fill(binToFill,wgtPuDown);
	  allPlots[channel+"njetsnbtags_muEffUp"]->Fill(binToFill,wgtMuEffUp);
	  allPlots[channel+"njetsnbtags_muEffDown"]->Fill(binToFill,wgtMuEffDown);
	  allPlots[channel+"njetsnbtags_umetUp"]->Fill(binToFill,wgt);
	  allPlots[channel+"njetsnbtags_umetDown"]->Fill(binToFill,wgt);
	  
	  std::vector<TString> catsToFill(2,Form("%dj",nJetsCat));
	  catsToFill[1]+= Form("%dt",nBtagsCat);
	  
	  std::map<Int_t,Double_t>::iterator rIt=lumiMap.find(ev.run);
	  for(size_t icat=0; icat<2; icat++)
	    {
	      TString tag=catsToFill[icat];
	      if(rIt!=lumiMap.end()) {
		Int_t runCtr=std::distance(lumiMap.begin(),rIt);
		allPlots[channel+"ratevsrun_"+tag]->Fill(runCtr,1.e+6/rIt->second);
	      }
	      allPlots[channel+"lpt_"+tag]->Fill(ev.l_pt[theLep],wgt);
	      allPlots[channel+"leta_"+tag]->Fill(ev.l_eta[theLep],wgt);
	      allPlots[channel+"jpt_"+tag]->Fill(ev.j_pt[ selJetsIdx[0] ],wgt);
	      allPlots[channel+"jeta_"+tag]->Fill(fabs(ev.j_eta[ selJetsIdx[0] ]),wgt);
	      allPlots[channel+"csv_"+tag]->Fill(ev.j_csv[ selJetsIdx[0] ],wgt);
	      allPlots[channel+"ht_"+tag]->Fill(htsum,wgt);
	      allPlots[channel+"nvtx_"+tag]->Fill(ev.nvtx,wgt);
	      allPlots[channel+"met_"+tag]->Fill(ev.met_pt,wgt);
	      allPlots[channel+"metphi_"+tag]->Fill(ev.met_phi,wgt);
	      allPlots[channel+"mt_nom_"+tag]->Fill(mt,wgt);
	      allPlots[channel+"mt_qcdScaleDown_"+tag]->Fill(mt,wgtQCDScaleLo);
	      allPlots[channel+"mt_qcdScaleUp_"+tag]->Fill(mt,wgtQCDScaleHi);
	      allPlots[channel+"mt_puUp_"+tag]->Fill(mt,wgtPuUp);
	      allPlots[channel+"mt_puDown_"+tag]->Fill(mt,wgtPuDown);
	      allPlots[channel+"mt_muEffUp_"+tag]->Fill(mt,wgtMuEffUp);
	      allPlots[channel+"mt_muEffDown_"+tag]->Fill(mt,wgtMuEffDown);
	      allPlots[channel+"mt_eEffUp_"+tag]->Fill(mt,wgtElEffUp);
	      allPlots[channel+"mt_eEffDown_"+tag]->Fill(mt,wgtElEffDown);	      
	      allPlots[channel+"mt_umetUp_"+tag]->Fill(mtUMetup,wgtElEffUp);
	      allPlots[channel+"mt_umetDown_"+tag]->Fill(mtUMetdown,wgtElEffDown);	      

	      if(nBtagsCat>0)
		{
		  allPlots[channel+"minmlb_nom_"+tag]->Fill(minMlb,wgt);
		  allPlots[channel+"minmlb_qcdScaleDown_"+tag]->Fill(minMlb,wgtQCDScaleLo);
		  allPlots[channel+"minmlb_qcdScaleUp_"+tag]->Fill(minMlb,wgtQCDScaleHi);
		  allPlots[channel+"minmlb_puUp_"+tag]->Fill(minMlb,wgtPuUp);
		  allPlots[channel+"minmlb_puDown_"+tag]->Fill(minMlb,wgtPuDown);
		  allPlots[channel+"minmlb_muEffUp_"+tag]->Fill(minMlb,wgtMuEffUp);
		  allPlots[channel+"minmlb_muEffDown_"+tag]->Fill(minMlb,wgtMuEffDown);
		  allPlots[channel+"minmlb_eEffUp_"+tag]->Fill(minMlb,wgtElEffUp);
		  allPlots[channel+"minmlb_eEffDown_"+tag]->Fill(minMlb,wgtElEffDown);
		  allPlots[channel+"minmlb_umetUp_"+tag]->Fill(minMlb,wgt);
		  allPlots[channel+"minmlb_umetDown_"+tag]->Fill(minMlb,wgt);
		}	  
	    }
	}

      if(nJetsJESHi>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJESHi,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots[channel+"njetsnbtags_jesUp"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots[channel+"mt_jesUp_"+tag]->Fill(mtJESup,wgt);
	  allPlots[channel+"minmlb_jesUp_"+tag]->Fill(minMlbJESup,wgt);
	}
      if(nJetsJESLo>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJESLo,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots[channel+"njetsnbtags_jesDown"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots[channel+"mt_jesDown_"+tag]->Fill(mtJESdown,wgt);
	  allPlots[channel+"minmlb_jesDown_"+tag]->Fill(minMlbJESdown,wgt);
	}
      if(nJetsJERHi>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJERHi,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots[channel+"njetsnbtags_jerUp"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots[channel+"mt_jerUp_"+tag]->Fill(mtJERup,wgt);
	  allPlots[channel+"minmlb_jerUp_"+tag]->Fill(minMlbJERup,wgt);
	}
      if(nJetsJERLo>=1)
	{
	  int nJetsCat=TMath::Min((int)nJetsJERLo,(int)4);
          int nBtagsCat=TMath::Min((int)nBtags,(int)2);
	  if(nBtagsCat>nJetsCat) nBtagsCat=nJetsCat;
	  int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots[channel+"njetsnbtags_jerDown"]->Fill(binToFill,wgt);
	  TString tag(Form("%dj%dt",nJetsCat,nBtagsCat));
	  allPlots[channel+"mt_jerDown_"+tag]->Fill(mtJERdown,wgt);
	  allPlots[channel+"minmlb_jerDown_"+tag]->Fill(minMlbJERdown,wgt);
	}
      if(nJets>=1)
	{
	  int nJetsCat=TMath::Min((int)nJets,(int)4);
          
	  int nBtagsCat=TMath::Min((int)nBtagsBeffLo,(int)2);
          int binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);
	  allPlots[channel+"njetsnbtags_beffDown"]->Fill(binToFill,wgt); 
	  
	  nBtagsCat=TMath::Min((int)nBtagsBeffHi,(int)2);
	  binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots[channel+"njetsnbtags_beffUp"]->Fill(binToFill,wgt); 

	  nBtagsCat=TMath::Min((int)nBtagsMistagLo,(int)2);
	  binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots[channel+"njetsnbtags_mistagDown"]->Fill(binToFill,wgt); 

	  nBtagsCat=TMath::Min((int)nBtagsMistagHi,(int)2);
	  binToFill=getBtagCatBinToFill(nBtagsCat,nJetsCat);          
	  allPlots[channel+"njetsnbtags_mistagUp"]->Fill(binToFill,wgt); 
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

                                                                    allPlots[channel+"raweventflowold_nom"]->Fill(1., wgt); 
      if(pass2Jets)                                               { allPlots[channel+"raweventflowold_nom"]->Fill(2., wgt); allPlots[channel+"eventflowold_nom"]->Fill(0., wgt); }
      if(pass2Jets && passMet)                                    { allPlots[channel+"raweventflowold_nom"]->Fill(3., wgt); allPlots[channel+"eventflowold_nom"]->Fill(1., wgt); }
      if(pass2Jets && passMet && pass1tau)                        { allPlots[channel+"raweventflowold_nom"]->Fill(4., wgt); allPlots[channel+"eventflowold_nom"]->Fill(2., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS)              { allPlots[channel+"raweventflowold_nom"]->Fill(5., wgt); allPlots[channel+"eventflowold_nom"]->Fill(3., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && pass1btag) { allPlots[channel+"raweventflowold_nom"]->Fill(6., wgt); allPlots[channel+"eventflowold_nom"]->Fill(4., wgt); }

                                                       allPlots[channel+"raweventflownocat_nom"]->Fill(1., wgt);
      if(pass2Jets)                                  { allPlots[channel+"raweventflownocat_nom"]->Fill(2., wgt); allPlots[channel+"eventflownocat_nom"]->Fill(0., wgt); }
      if(pass2Jets && passMet)                       { allPlots[channel+"raweventflownocat_nom"]->Fill(3., wgt); allPlots[channel+"eventflownocat_nom"]->Fill(1., wgt); }
      if(pass2Jets && passMet && pass1tau)           { allPlots[channel+"raweventflownocat_nom"]->Fill(4., wgt); allPlots[channel+"eventflownocat_nom"]->Fill(2., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS) { allPlots[channel+"raweventflownocat_nom"]->Fill(5., wgt); allPlots[channel+"eventflownocat_nom"]->Fill(3., wgt); }

                                                                       allPlots[channel+"raweventflowcat_nom"]->Fill(1., wgt);
      if(pass2Jets)                                                  { allPlots[channel+"raweventflowcat_nom"]->Fill(2., wgt); allPlots[channel+"eventflowcat_nom"]->Fill(0., wgt); }
      if(pass2Jets && passMet)                                       { allPlots[channel+"raweventflowcat_nom"]->Fill(3., wgt); allPlots[channel+"eventflowcat_nom"]->Fill(1., wgt); }
      if(pass2Jets && passMet && pass1tau)                           { allPlots[channel+"raweventflowcat_nom"]->Fill(4., wgt); allPlots[channel+"eventflowcat_nom"]->Fill(2., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS)                 { allPlots[channel+"raweventflowcat_nom"]->Fill(5., wgt); allPlots[channel+"eventflowcat_nom"]->Fill(3., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && passCat0btag) { allPlots[channel+"raweventflowcat_nom"]->Fill(6., wgt); allPlots[channel+"eventflowcat_nom"]->Fill(4., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && passCat1btag) { allPlots[channel+"raweventflowcat_nom"]->Fill(7., wgt); allPlots[channel+"eventflowcat_nom"]->Fill(5., wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && passCat2btag) { allPlots[channel+"raweventflowcat_nom"]->Fill(8., wgt); allPlots[channel+"eventflowcat_nom"]->Fill(6., wgt); }


      std::vector<TString> cuts; cuts.clear();

      if(pass2Jets)                                                  { cuts.push_back("1l2j_");} //FillCutPlots(allPlots, channel, cuts[0], "nom", ev, theLep, selJetsIdx, wgt); }
      if(pass2Jets && passMet)                                       { cuts.push_back("met_" );} //FillCutPlots(allPlots, channel, cuts[1], "nom", ev, theLep, selJetsIdx, wgt); }
      if(pass2Jets && passMet && pass1tau)                           { cuts.push_back("1tau_");} //FillCutPlots(allPlots, channel, cuts[2], "nom", ev, theLep, selJetsIdx, wgt); }
      if(pass2Jets && passMet && pass1tau && passOS)                 { cuts.push_back("os_"  );} //FillCutPlots(allPlots, channel, cuts[3], "nom", ev, theLep, selJetsIdx, wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && pass1btag)    { cuts.push_back("1bi_" );} //FillCutPlots(allPlots, channel, cuts[4], "nom", ev, theLep, selJetsIdx, wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && passCat0btag) { cuts.push_back("0b_"  );} //FillCutPlots(allPlots, channel, cuts[5], "nom", ev, theLep, selJetsIdx, wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && passCat1btag) { cuts.push_back("1b_"  );} //FillCutPlots(allPlots, channel, cuts[6], "nom", ev, theLep, selJetsIdx, wgt); }
      if(pass2Jets && passMet && pass1tau && passOS && passCat2btag) { cuts.push_back("2b_"  );}  //FillCutPlots(allPlots, channel, cuts[7], "nom", ev, theLep, selJetsIdx, wgt); }

      for(std::vector<TString>::iterator icut=cuts.begin(); icut!=cuts.end(); ++icut)
        {
          TString cut(*icut);
          TString syst("nom");
          allPlots[channel+cut+"lpt_"+syst]    ->Fill(ev.l_pt[theLep]                 ,wgt);
          allPlots[channel+cut+"leta_"+syst]   ->Fill(ev.l_eta[theLep]                ,wgt);
          allPlots[channel+cut+"jpt_"+syst]    ->Fill(ev.j_pt[ selJetsIdx[0] ]        ,wgt);
          allPlots[channel+cut+"jeta_"+syst]   ->Fill(fabs(ev.j_eta[ selJetsIdx[0] ]) ,wgt);
          allPlots[channel+cut+"tpt_"+syst]    ->Fill(ev.t_pt[theTau]                 ,wgt);
          allPlots[channel+cut+"teta_"+syst]   ->Fill(ev.t_eta[theTau]                ,wgt);
          allPlots[channel+cut+"csv_"+syst]    ->Fill(ev.j_csv[ selJetsIdx[0] ]       ,wgt);
          //allPlots[channel+cut+"ht_"+syst]     ->Fill(htsum                           ,wgt);
          allPlots[channel+cut+"nvtx_"+syst]   ->Fill(ev.nvtx                         ,wgt);
          allPlots[channel+cut+"nvtxraw_"+syst]->Fill(ev.nvtx                         ,wgt/puWeight[0]);
          allPlots[channel+cut+"met_"+syst]    ->Fill(ev.met_pt                       ,wgt);
          allPlots[channel+cut+"metphi_"+syst] ->Fill(ev.met_phi                      ,wgt);
          allPlots[channel+cut+"nbjets_"+syst] ->Fill(nBtags<4?nBtags:4               ,wgt);

        }
     
      if(pass2Jets && passMet && pass1tau && passOS && passCat0btag) allPlots[channel+"btagcat_nom"]->Fill(0., wgt);
      if(pass2Jets && passMet && pass1tau && passOS && passCat1btag) allPlots[channel+"btagcat_nom"]->Fill(1., wgt);
      if(pass2Jets && passMet && pass1tau && passOS && passCat2btag) allPlots[channel+"btagcat_nom"]->Fill(2., wgt);






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
