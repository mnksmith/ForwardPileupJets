/**************************************************************************
 **
 **   File:         Analysis_ForwardPileupJets.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      M. Swiatlowski
 **
 **************************************************************************/

#define Analysis_ForwardPileupJets_cxx

#include "Analysis_ForwardPileupJets.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "TKey.h"
#include "TObjString.h"
#include "TObjArray.h"


///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_ForwardPileupJets::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_ForwardPileupJets: DEBUG In WorkerBegin()" << endl;

  Analysis_JetMET_Base::WorkerBegin();
  fDoLeptonSelection = false;
  ChainCfg()->Get("doLeptonSelection",  fDoLeptonSelection);


  if (Debug()) cout << "Analysis_ForwardPileupJets: DEBUG Finish WorkerBegin()" << endl;  
} 

///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_ForwardPileupJets::ProcessEvent()
{

  if (Debug()) cout << "Analysis_ForwardPileupJets: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();

  // Event Selection: don't fill histograms if event selection is false
  //                  Set bool EventSelection for use in other (chained) analyses
  if( EventSelection() )       {Set("EventSelection", true );} 
  else                         {Set("EventSelection", false); return true;}
  
  // Jet Collections: remove overlap between muons and jets 
  //                  this makes a vector of jets with name jetsAntiKt4LCTopoGood
  MakeJetCollections("AntiKt4LCTopo");

  if(Debug()) cout << "Making jets and jet Plots!" << endl;

  static const MomKey JetKey4LC("AntiKt4LCTopoGood");
  static const MomKey JetKey4Tru("AntiKt4Truth");
  static const MomKey JetKey4IT("InTimeAntiKt4Truth");
  //Show();
  //cout << "NPV = " << Int("NPV") << endl;
  AnaKey JetKeyA(JetKey4LC.Data());
  Fill(JetKeyA+"_NPV", Int("NPV"), Weight(), 100, 0., 100.);
  Fill(JetKeyA+"_mu", Float("averageIntPerXing"), Weight(), 100, 0., 100.); 
  MakeJetPlots(JetKey4LC,JetKey4Tru,JetKey4IT);
  
  return true;
}

void Analysis_ForwardPileupJets::MakeJetPlots(const MomKey JetKey1, const MomKey JetKey2, const MomKey JetKey3){
  AnaKey JetKeyS(JetKey1.Data());
  double delR = 0.;
  double delRit = 0.;

  AddVec("jetsHardScatter");
  AddVec("jetsPileup");

  for(int iJet = 0; iJet < jets(JetKey1); iJet++){   
    int inDelRcut = 1; // for Pileup: stays 1 if there are no truth jets with pt > 4GeV within delR > 0.6
    jet(iJet, JetKey1).Set("isHardScatter", 0);
    jet(iJet, JetKey1).Set("isPileup", 0);
    if(isMC()){
      jet(iJet, JetKey1).Set("isQCDPileup1", 0);
      jet(iJet, JetKey1).Set("isQCDPileup2", 0);
      jet(iJet, JetKey1).Set("isQCDPileup3", 0);
      jet(iJet, JetKey1).Set("isQCDPileup4", 0);
    }
    Fill(JetKeyS+"_allj_pt", jet(iJet, JetKey1).p.Perp(), Weight(), 100, 0., 300.);
    Fill(JetKeyS+"_allj_m",  jet(iJet, JetKey1).p.M()   , Weight(), 100, 0., 300.);
    if(isMC()){
      for(int jJet = 0; jJet < jets(JetKey2); jJet++){
        delR = jet(iJet, JetKey1).p.DeltaR(jet(jJet, JetKey2).p);
        Fill("DeltaR_LCTopo_Truth", delR, Weight(), 100, 0., 10.);
    
        if(jet(jJet, JetKey2).p.Perp() > 10. && delR < 0.3){
          jet(iJet, JetKey1).Set("isHardScatter", 1);
          Add("jetsHardScatter", &(jet(iJet, JetKey1)));
        }
        if(jet(jJet, JetKey2).p.Perp() > 4. && delR < 0.6){
          inDelRcut = 0;
        }
      }
       

      if(inDelRcut == 1){
        jet(iJet, JetKey1).Set("isPileup", 1);
        Add("jetsPileup", &(jet(iJet, JetKey1)));
      
        for(int kJet=0; kJet < jets(JetKey3); kJet++){
          delRit = jet(iJet, JetKey1).p.DeltaR(jet(kJet, JetKey3).p);
          if(jet(kJet, JetKey3).p.Perp() > 15. && delRit < 0.4){          //criteria for a jet to be QCD pileup
            jet(iJet, JetKey1).Set("isQCDPileup1", 1);
          }
          if(jet(kJet, JetKey3).p.Perp()/jet(iJet, JetKey1).p.Perp() > 0.6 && delRit < 0.4){          //criteria for a jet to be QCD pileup
            jet(iJet, JetKey1).Set("isQCDPileup2", 1);
          }
          if(jet(kJet, JetKey3).p.Perp()/jet(iJet, JetKey1).p.Perp() > 0.7 && delRit < 0.4){          //criteria for a jet to be QCD pileup
            jet(iJet, JetKey1).Set("isQCDPileup3", 1);
          }
          if(jet(kJet, JetKey3).p.Perp()/jet(iJet, JetKey1).p.Perp() > 0.8 && delRit < 0.4){          //criteria for a jet to be QCD pileup
            jet(iJet, JetKey1).Set("isQCDPileup4", 1);
          }
        }
      }
    } 

    double sumNum = 0;
    double sumDenom = 0;
    double delRClust = 0;
    double sumPt1 = 0;
    double sumPt2 = 0;
    double sumPt3 = 0;
    double sumPt4 = 0;
    for(int iC = 0; iC < jet(iJet, JetKey1).Objs("constituents"); ++iC){
      Particle *constituent = (Particle*) jet(iJet, JetKey1).Obj("constituents", iC);
      delRClust = jet(iJet, JetKey1).p.DeltaR(constituent->p);
      sumNum = sumNum + pow(delRClust*constituent->p.Pt(), 2);
      sumDenom = sumDenom + pow(constituent->p.Pt(), 2);
    
      if(delRClust < 0.1){
        sumPt1 = sumPt1 + constituent->p.Pt();
      }
      if(delRClust > 0.1 && delRClust < 0.2){
        sumPt2 = sumPt2 + constituent->p.Pt();
      }
      if(delRClust > 0.2 && delRClust < 0.3){
        sumPt3 = sumPt3 + constituent->p.Pt();
      }
      if(delRClust > 0.3 && delRClust < 0.4){
        sumPt4 = sumPt4 + constituent->p.Pt();
      }
    }

    if(jet(iJet, JetKey1).Objs("constituents")>0){
      Fill(JetKeyS+"_allj_width", jet(iJet, JetKey1).Float("WIDTH"), Weight(), 100, 0., 1.);
      
      Fill(JetKeyS+"_allj_delRsqr", sumNum/sumDenom, Weight(), 100, 0., 0.5);
      jet(iJet, JetKey1).Set("delRsqr", sumNum/sumDenom);
      
      Fill(JetKeyS+"_allj_0_delR_1", sumPt1/jet(iJet, JetKey1).p.Perp(), Weight(), 100, 0., 1.);
      Fill(JetKeyS+"_allj_1_delR_2", sumPt2/jet(iJet, JetKey1).p.Perp(), Weight(), 100, 0., 1.);
      Fill(JetKeyS+"_allj_2_delR_3", sumPt3/jet(iJet, JetKey1).p.Perp(), Weight(), 100, 0., 1.);
      Fill(JetKeyS+"_allj_3_delR_4", sumPt4/jet(iJet, JetKey1).p.Perp(), Weight(), 100, 0., 1.);
      jet(iJet, JetKey1).Set("delR_01", sumPt1/jet(iJet, JetKey1).p.Perp());
      jet(iJet, JetKey1).Set("delR_12", sumPt2/jet(iJet, JetKey1).p.Perp());
      jet(iJet, JetKey1).Set("delR_23", sumPt3/jet(iJet, JetKey1).p.Perp());
      jet(iJet, JetKey1).Set("delR_34", sumPt4/jet(iJet, JetKey1).p.Perp());

    }
    if(jet(iJet, JetKey1).Objs("constituents")<=0){
      jet(iJet, JetKey1).Set("delRsqr", -1.);
      jet(iJet, JetKey1).Set("delR_01", -1.);
      jet(iJet, JetKey1).Set("delR_12", -1.);
      jet(iJet, JetKey1).Set("delR_23", -1.);
      jet(iJet, JetKey1).Set("delR_34", -1.);
    }

  }
    
}

// --------------------------------------------------
// Make jet collections with overlap removal
// ----------------------------------------------
void Analysis_ForwardPileupJets::MakeJetCollections(const MomKey JetKey){
  if(Debug()) cout <<"Analysis_ForwardPileupJets::MakeJetCollections Start" << endl;
    AddVec("jets"+JetKey+"Good");

    for(int iJet = 0; iJet < jets(JetKey); iJet++){
        Particle  *myjet        = &(jet(iJet, JetKey));
        
        bool overlap(false);
        if(Exists("muonsgood")){
            for(int iMu=0; iMu < muons("good"); iMu++){ // Muon - Jet Overlap ---------
                float dR = muon(iMu,"good").p.DeltaR(myjet->p);
                if(dR<0.4) overlap=true;
            }
        }
        if(overlap) continue;                          // ---------------------------

        // make collections
        Add("jets"+JetKey+"Good", myjet); 
    }
  if(Debug()) cout <<"Analysis_ForwardPileupJets::MakeJetCollections End" << endl;
}
///========================================
/// Lepton Selection 
///========================================

bool Analysis_ForwardPileupJets::EventSelection(){
  if(Debug()) cout <<"Analysis_ForwardPileupJets::EventSelection Begin" << endl;

  if (fDoLeptonSelection){
      if(Debug()) cout << "Analysis_ForwardPileupJets: doLeptonSelection"     << endl;
      if( ! (Objs("muonsgood")==2))                                                  {return false; } 
      if( ! Bool("TopSelection_GRL"))                                                {return false; }
      if( ! Bool("TopSelection_LarError"))                                           {return false; }
      if( ! Bool("TopSelection_NoBadLooseJets"))                                     {return false; }
      if( ! (muon(0, "good").Float("charge")*muon(1, "good").Float("charge") == -1)) {return false; }
      if( ! Bool("1muontrigger"))                                                    {return false;}

      // reconstruct Z
      AddVec("recosZCandidate");
      Particle* Z = new Particle();
      Z->p = muon(0, "good").p + muon(1, "good").p;
      Add("recosZCandidate", Z);
      if( ! ((reco(0, "ZCandidate").p.M() > 71) && 
             (reco(0, "ZCandidate").p.M()  < 111)))           return false;  
  }

  if(Debug()) cout <<"Analysis_ForwardPileupJets::EventSelection End" << endl;
  return true;
}

///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_ForwardPileupJets::WorkerTerminate()
{

  // Nothing

}

///==========================================
/// MakeJetMassCut: return true if lead jet of type specified is outside of mass window
///==========================================
bool Analysis_ForwardPileupJets::MakeJetMassCut(MomKey JetKey){
  if(jets(JetKey) < 1)
    return false;
  float m = jet(0,JetKey).p.M();
  if(m < 50. || m > 110.)
    return false; 
  return true;
}

///==========================================
/// MakeJetPtCut: return true if lead jet of type specified is outside of pT window
///==========================================
bool Analysis_ForwardPileupJets::MakeJetPtCut(MomKey JetKey){
  if(jets(JetKey) < 1)
    return false;
  float pt = jet(0,JetKey).p.Perp();
  if(pt <= 200. || pt > 350.)
    return false; 
  return true;
}
