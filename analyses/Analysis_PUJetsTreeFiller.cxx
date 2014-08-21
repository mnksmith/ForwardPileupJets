/**************************************************************************
 **
 **   File:         Analysis_PUJetsTreeFiller.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      P. Nef
 **
 **************************************************************************/

#define Analysis_PUJetsTreeFiller_cxx

#include "Analysis_PUJetsTreeFiller.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>      // std::setprecision
#include <cstdlib>
#include <sstream>
#include "TKey.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"


///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_PUJetsTreeFiller::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_PUJetsTreeFiller: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();


  // trees -------------------
  fEventTree = new TTree("EventTree", "Tree with event-by-event variables");
  AddBranches(fEventTree);
} 


///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_PUJetsTreeFiller::ProcessEvent()
{

  if (Debug()) cout << "Analysis_PUJetsTreeFiller: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		            << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();
  
  // Event Selection goes here... ---------------------------
  // only fill tree if event selection is good. this is set in Analysis_ForwardPileupJets.cxx
  if(Bool("EventSelection")==false) return true;  
  

  // Fill Tree-----------------------------------------------
  FillTree("AntiKt4LCTopoGood");


  // end----------------------------------------------------------------------
  if (Debug()) cout << "Analysis_PUJetsTreeFiller: DEBUG End ProcessEvent(): RunNumber = " << RunNumber() 
		            << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;
  return true;
}



///=========================================
/// WorkerTerminate: clean up
///=========================================
void Analysis_PUJetsTreeFiller::WorkerTerminate()
{
    fEventTree->Write();


  // Nothing more

}

///===========================================
/// Fill Tree
///========================================
void Analysis_PUJetsTreeFiller::FillTree(const MomKey JetKey){
  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::FillTree Start" << endl;

  for(int iJ=0; iJ<jets(JetKey); ++iJ){
    Particle *myjet = &(jet(iJ, JetKey)); 
    if(myjet->p.Pt()>20){
      ResetBranches(fEventTree);
      fTNJets++;
      fTJetIndex = iJ;
      FillEventVars(fEventTree, JetKey, myjet);
      fEventTree->Fill();
    }
  
  }  
  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::FillTree End" << endl;
  return;
}

///=============================================
/// Add Branches To Tree
///=============================================
void Analysis_PUJetsTreeFiller::AddBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::AddBranches Start" << endl;
  
    // Event Info
    tree->Branch("EventNumber",               &fTEventNumber,            "EventNumber/I");
    tree->Branch("RunNumber",                 &fTRunNumber,              "RunNumber/I");
    tree->Branch("Weight" ,                   &fTWeight,                 "Weight/F");
    tree->Branch("Mu" ,                       &fTMu,                     "Mu/F");
    tree->Branch("NPVtruth" ,                 &fTNPVtruth,               "NPVtruth/I");
    tree->Branch("NPV" ,                      &fTNPV,                    "NPV/I");
  
    // Jet vars --------------------------------------------------------------------------------------
    tree->Branch("JetIndex",                  &fTJetIndex,               "JetIndex/I");
    tree->Branch("NJets",                     &fTNJets,                  "NJets/I");
    tree->Branch("Jpt",                       &fTJPt,                    "Jpt/F");
    tree->Branch("Jeta",                      &fTJEta,                   "Jeta/F");
    tree->Branch("Jphi",                      &fTJPhi,                   "Jphi/F");
    tree->Branch("Jm",                        &fTJM,                     "Jm/F");
    tree->Branch("Jwidth",                    &fTJWidth,                 "Jwidth/F");
    tree->Branch("delRsqr",                   &fTdelRsqr,                "delRsqr/F");
    tree->Branch("delR_01",                   &fTdelR_01,                "delR_01/F");
    tree->Branch("delR_12",                   &fTdelR_12,                "delR_12/F");
    tree->Branch("delR_23",                   &fTdelR_23,                "delR_23/F");
    tree->Branch("delR_34",                   &fTdelR_34,                "delR_34/F");
    tree->Branch("isHardScatter",             &fTisHardScatter,          "isHardScatter/I");
    tree->Branch("isPileup",                  &fTisPileup,               "isPileup/I");
    tree->Branch("isQCDPileup1",              &fTisQCDPileup1,           "isQCDPileup1/I");
    tree->Branch("isQCDPileup2",              &fTisQCDPileup2,           "isQCDPileup2/I");
    tree->Branch("isQCDPileup3",              &fTisQCDPileup3,           "isQCDPileup3/I");
    tree->Branch("isQCDPileup4",              &fTisQCDPileup4,           "isQCDPileup4/I");
    tree->Branch("Timing",                    &fTTiming,                 "Timing/F");
    tree->Branch("NumTowers",                 &fTNumTowers,              "NumTowers/I");
    tree->Branch("NClus",                     &fTNClus,                  "NClus/I");
    tree->Branch("ClusPt",                    &fTClusPt,                 "ClusPt[NClus]/F");
    tree->Branch("ClusEta",                   &fTClusEta,                "ClusEta[NClus]/F");
    tree->Branch("ClusPhi",                   &fTClusPhi,                "ClusPhi[NClus]/F");
    tree->Branch("centlam",                   &fTcentlam,                "centlam[NClus]/F");
    tree->Branch("Edens",                     &fTEdens,                  "Edens[NClus]/F");
    tree->Branch("cellmaxfrac",               &fTcellmaxfrac,            "cellmaxfrac[NClus]/F");
    tree->Branch("long",                      &fTlong,                   "long[NClus]/F");
    tree->Branch("lat",                       &fTlat,                    "lat[NClus]/F");
    tree->Branch("secondLam",                 &fTsecondLam,              "secondLam[NClus]/F");
    tree->Branch("secondR",                   &fTsecondR,                "secondR[NClus]/F");
    tree->Branch("iso",                       &fTiso,                    "iso[NClus]/F");
    tree->Branch("sig",                       &fTsig,                    "sig[NClus]/F");
    tree->Branch("Epos",                      &fTEpos,                   "Epos[NClus]/F");
    tree->Branch("delTheta",                  &fTdelTheta,               "delTheta[NClus]/F");
    tree->Branch("delPhi",                    &fTdelPhi,                 "delPhi[NClus]/F");
    tree->Branch("centmag",                   &fTcentmag,                "centmag[NClus]/F");

  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::AddBranches End" << endl;
    return;
}

void Analysis_PUJetsTreeFiller::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::ResetBranches Start" << endl;
  
    // Event Info
    fTEventNumber           = -999;
    fTRunNumber             = -999;
    fTWeight                = -999;
    fTNPVtruth              = -999;
    fTNPV                   = -999;
    fTMu                    = -999;
    
    //Jet Info
    fTJetIndex           = -999;
    fTJPt              = -999.99;
    fTJEta             = -999.99;
    fTJPhi             = -999.99;
    fTJM               = -999.99;
    fTJWidth           = -999.99;
    fTdelRsqr          = -999.99;
    fTdelR_01          = -999.99;
    fTdelR_12          = -999.99;
    fTdelR_23          = -999.99;
    fTdelR_34          = -999.99;
    fTisHardScatter    = -999;
    fTisPileup         = -999;
    fTisQCDPileup1     = -999;
    fTisQCDPileup2     = -999;
    fTisQCDPileup3     = -999;
    fTisQCDPileup4     = -999;
    fTNumTowers        = -999;
    fTTiming          = -999.99;
    fTNClus            = 0;
    for(int iC=0; iC<MaxNCluster; ++iC){
	fTClusPt[iC]     = -999.99;
        fTClusEta[iC]    = -999.99;
        fTClusPhi[iC]    = -999.99;
        fTcentlam[iC]    = -999.99;
        fTEdens[iC]      = -999.99;
        fTcellmaxfrac[iC]= -999.99;
	fTlong[iC]       = -999.99;
 	fTlat[iC]        = -999.99;
        fTsecondLam[iC]  = -999.99;
        fTsecondR[iC]    = -999.99;
	fTiso[iC]        = -999.99;
 	fTsig[iC]        = -999.99;
        fTEpos[iC]       = -999.99;
 	fTdelTheta[iC]   = -999.99;
        fTdelPhi[iC]     = -999.99;
        fTcentmag[iC]    = -999.99;
    }

  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::ResetBranches End" << endl;
    return;
}


void Analysis_PUJetsTreeFiller::FillEventVars(TTree *tree, const MomKey JetKey, Particle *myjet){
  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::FillEventVars Begin" << endl;
    
    // Event Info -----------------------------------
    fTEventNumber                 = Int("EventNumber");
    fTRunNumber                   = Int("RunNumber");
    fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
    fTNPV                         = Exists("NPV")?      Int("NPV"):-1;
    fTWeight                      = DefaultWeight();
    fTMu                          = Float("averageIntPerXing");                  

    // Jet Info ----------------------
    fTJPt         = myjet->p.Pt();
    fTJEta        = myjet->p.Eta();
    fTJPhi        = myjet->p.Phi();
    fTJM          = myjet->p.M();
    fTJWidth      = myjet->Float("WIDTH");
    fTdelRsqr     = myjet->Float("delRsqr");
    fTdelR_01     = myjet->Float("delR_01");
    fTdelR_12     = myjet->Float("delR_12");
    fTdelR_23     = myjet->Float("delR_23");
    fTdelR_34     = myjet->Float("delR_34");
    fTisHardScatter = myjet->Int("isHardScatter");
    fTisPileup      = myjet->Int("isPileup");
    if(isMC()){
      fTisQCDPileup1   = myjet->Int("isQCDPileup1");
      fTisQCDPileup2   = myjet->Int("isQCDPileup2");
      fTisQCDPileup3   = myjet->Int("isQCDPileup3");
      fTisQCDPileup4   = myjet->Int("isQCDPileup4");
    }
    fTNumTowers     = myjet->Float("NumTowers");   
    fTTiming     = myjet->Float("Timing");

    for(int iC=0; iC<myjet->Objs("constituents"); ++iC){
	Particle *con = (Particle*) myjet->Obj("constituents",iC);
	if(fTNClus == MaxNCluster) continue;
	fTClusPt[iC] = con->p.Pt();
        fTClusEta[iC] = con->p.Eta();
        fTClusPhi[iC] = con->p.Phi();
        if(fTClusPt[iC]>10){
          fTcentlam[iC] = con->Float("centerlambda");
          fTEdens[iC] = con->Float("firstEdens");
	  fTcellmaxfrac[iC] = con->Float("cellmaxfrac");
  	  fTlong[iC] = con->Float("longitudinal");
          fTlat[iC] = con->Float("lateral");
 	  fTsecondLam[iC] = con->Float("secondlambda");
	  fTsecondR[iC] = con->Float("secondR");
	  //fTiso[iC] = con->Float("isolation");
	  //fTsig[iC] = con->Float("significance");
	  //fTEpos[iC] = con->Float("eng_pos");
	  fTdelTheta[iC] = con->Float("deltaTheta");
	  fTdelPhi[iC] = con->Float("deltaPhi");
	  fTcentmag[iC] = con->Float("centermag");
	}
        fTNClus++;
    }
   

  if(Debug()) cout <<"Analysis_PUJetsTreeFiller::FillEventVars End" << endl;
  return;
}

