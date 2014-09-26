#include "../scripts/runLocal.C"
#include "../scripts/runProof.C"
#include "../scripts/helperFunc.C"
#include "../scripts/helperJetMETCommon.C"
#include "../scripts/loadLibraries.C"
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

void runForwardPU_Zmumu(TString mode       = "local",         // local, lite, or cluster
TString identifier = "ForwardPileupJets",                      // tag 
//TString dataset = "Zmumu_PowhegPythia8_MC12_COMMON.jetmet2012",
TString dataset = "Zmumu_PowhegPythia8_MC12_COMMON_short.jetmet2012",

//TString dataset = "MuonsSMWZPeriodAB_short.jetmet2012",
//TString dataset = "MuonsSMWZPeriodAB_short2.jetmet2012",
//TString dataset = "MuonsSMWZPeriodAB.jetmet2012",
//TString dataset = "MuonsSMWZPeriodCD.jetmet2012",
//TString dataset = "MuonsSMWZPeriodE1.jetmet2012",
//TString dataset = "MuonsSMWZPeriodE2.jetmet2012",
//TString dataset = "MuonsSMWZPeriodE3.jetmet2012",
//TString dataset = "MuonsSMWZPeriodE4.jetmet2012",
//TString dataset = "MuonsSMWZPeriodE5.jetmet2012",
//TString dataset = "MuonsSMWZPeriodE6.jetmet2012",
//TString dataset = "MuonsSMWZPeriodG1.jetmet2012",
//TString dataset = "MuonsSMWZPeriodG2.jetmet2012",
//TString dataset = "MuonsSMWZPeriodG3.jetmet2012",
//TString dataset = "MuonsSMWZPeriodH.jetmet2012",
//TString dataset = "MuonsSMWZPeriodI.jetmet2012",
//TString dataset = "MuonsSMWZPeriodJ.jetmet2012",
//TString dataset = "MuonsSMWZPeriodL.jetmet2012",

TString username   = "mnks",                               // username (e.g. swiatlow, fizisist)
bool mcweights     = true,                                 // use mc weights?
bool debug         = false,                                // debug mode
Long64_t nentries  = 10                              // nevents
    ) 
{ 
    
    TString date(currentDateTime().c_str());
    identifier = date+"_"+identifier;

    ///----------------------------------------------------------------
    /// Load libraries , set the config file, treenam, and cluster info
    ///----------------------------------------------------------------

    cout << "trying to load libraries" << endl; 
    loadLibraries();

    cout << " Libraries loaded " << endl;

    // SetConfig 
    TString configfile("../config/pileupstudies.config");

    
    // Best to leave alone  
    TString pathLite("");
    TString pathCluster("root://atlprf01.slac.stanford.edu:2094//atlas/output/");
    pathCluster.Append(username);
    pathCluster.Append("/");

    // Determine eventbuilder from dataset name
    TString eventbuilder(dataset);
    eventbuilder.Remove(0,eventbuilder.Last('.')+1); 
    TString treename("physics");
   
 
    ///----------------------------------------------------------------
    /// Filename paths, URLs for PROOF running
    ///----------------------------------------------------------------
    bool runCluster(false);
    TString url(mode);
    TString path("");
    if(mode.CompareTo("lite")==0) {
        url = "lite://";
        path = pathLite;
    }
    else if(mode.CompareTo("cluster")==0) {
        url = TString(username+"@atlprf01.slac.stanford.edu");
        path = pathCluster;
        runCluster = true;
    }
    
    // Make an options file, edit as needed
    TFile* options = new TFile("options.root","RECREATE");

    ///----------------------------------------------------------------
    /// Overall Configuration
    ///----------------------------------------------------------------
    bool doBasic          = true;
    bool doJetCalibrations= true;
    bool doJetTriggers    = false;
    bool doPRW            = true;
    bool doParentChild    = false;
    bool doTrack          = true;
    bool doVertex         = true;
    bool doLCCluster      = true; 
    bool doEMCluster      = true;
    bool doTruth          = true;
    bool doEMJets         = false;
    bool doLCJets         = true;
    bool doJet4           = true;
    bool doJet6           = false; 
    bool doVectorJVFLinks = true;
    bool doTruthJets      = true;
    bool doOOTtruthJet4   = false;
    bool doTruthLinks     = false;
    bool doPhotons        = false;
    bool doElectrons      = false;
    bool doConstitLinks   = true;
    bool doInTimeTruthJet4= true;
    bool doMuons          = true;
    bool doMuonTriggers   = true;
    bool doMuonTriggersMatch = false;
    bool doElectrons      = true;
    bool doMETRefFinal    = true;
    bool doMuSmear        = false;
    bool doLeptonSelection= true;
    
    bool doCOMMON         = false;
    if(dataset.Contains("COMMON")){
      doCOMMON = true;
    }

    bool doSMWZ           = false;
    if(dataset.Contains("SMWZ")){
      doSMWZ = true;
    }

    int  counterMax       = -1;
    float LUMI            =1;
    TString prwTypes      = "PeriodAB_lumi"; 



    ///----------------------------------------------------------------
    /// Nominal Configuration
    /// set up the actual analyses
    ///----------------------------------------------------------------
    
    Config* TopSelection = new Config("TopSelection",configfile);
    TopSelection->Set("ANALYSIS","TopCommonSelection");
    TopSelection->Set("DEBUG",debug);

    Config* ForwardPileupJets = new Config("ForwardPileupJets",configfile);
    ForwardPileupJets->Set("ANALYSIS","ForwardPileupJets");
    ForwardPileupJets->Set("DEBUG",debug);

    Config* PUJetsTreeFiller = new Config("PUJetsTreeFiller",configfile);
    PUJetsTreeFiller->Set("ANALYSIS","PUJetsTreeFiller");
    PUJetsTreeFiller->Set("DEBUG",debug);
    
    Config* chain = new Config("chain",configfile);
    chain->AddVec("ANALYSIS");
    chain->Add("ANALYSIS",TopSelection);
    chain->Add("ANALYSIS",ForwardPileupJets);
    chain->Add("ANALYSIS",PUJetsTreeFiller);


    // set up configurations, this overwrites configs from configfile
    chain->Set("DOJETCALIBRATIONS",doJetCalibrations);
    chain->Set("DOJETTRIGGERS"   , doJetTriggers);
    chain->Set("DOJET4"          , doJet4          );
    chain->Set("DOJET6"          , doJet6          );
    chain->Set("DOVECTORJVFLINKS", doVectorJVFLinks);
    chain->Set("DOLCJETS"        , doLCJets        );
    chain->Set("DOEMJETS"        , doEMJets        );
    chain->Set("COUNTERMAX"      , counterMax      );
    chain->Set("DEBUG"           , debug           );
    chain->Set("MCWEIGHTS"       , mcweights       );
    chain->Set("PILE"            , doPRW           );
    chain->Set("DOBASIC"         , doBasic         );
    chain->Set("DOTRUTHLINKS"    , doTruthLinks    );
    chain->Set("DOCONSTITLINKS"  , doConstitLinks  );
    chain->Set("DOPARENTCHILD"   , doParentChild   );
    chain->Set("DOTRACK"         , doTrack         );
    chain->Set("DOLCCLUSTER"     , doLCCluster     );
    chain->Set("DOEMCLUSTER"     , doEMCluster     );
    chain->Set("DOTRUTH"         , doTruth         );
    chain->Set("DOTRUTHJETS"     , doTruthJets     );
    chain->Set("DOOOTTRUTHJET4"  , doOOTtruthJet4  );
    chain->Set("DOINTIMETRUTHJET4",doInTimeTruthJet4);
    chain->Set("DOVTX"           , doVertex        );
    chain->Set("DOPHOTON"        , doPhotons       );
    chain->Set("DOELECTRONS"     , doElectrons     );
    chain->Set("DOMUONS"         , doMuons         );
    chain->Set("DOMUONTRIGGERS"  , doMuonTriggers  );
    chain->Set("DOMUONTRIGGERSMATCH", doMuonTriggersMatch );
    chain->Set("DOMUSMEAR",        doMuSmear);
    chain->Set("DOCOMMON"        , doCOMMON        );
    chain->Set("JETTYPES"        , ""              );
    chain->Set("BTAGS"           , ""              ); 
    chain->Set("PRWTYPES"        , prwTypes        );
    chain->Set("LUMI",             LUMI);
    chain->Set("runCluster",       runCluster);
    chain->Set("DOMETREFFINAL"   , doMETRefFinal   );
    chain->Set("doLeptonSelection" , doLeptonSelection);
    chain->Set("DOSMWZ",           doSMWZ          );  

    chain->Write();


    
    WriteGRLObject("data12_8TeV.periodAB_HEAD_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml");
    WriteZllPRWO(options, "PeriodAB_lumi");
    WriteMuonUtilities(options);
    WriteBTagCalibObject(options,"MV1","0_7892");
    WriteJetCalibrationObjects(options);

    ///----------------------------------------------------------------
    /// ProofAna global Config object
    ///----------------------------------------------------------------
    Config* confProofAna = new Config("ProofAna");
    
    confProofAna->Set("DEBUG"          , false        );  // "false", 0, "0" etc. also works
    confProofAna->Set("SAVETIMERS"     , false        );  // ProofAna timer histos in output file   
    confProofAna->Set("IDENTIFIER"     , identifier   );
    confProofAna->Set("DATASET"        , dataset      );
    confProofAna->Set("OUTPUTPATH"     , path         );
    confProofAna->Set("EVENTBUILDER"   , eventbuilder );
    confProofAna->Set("MERGE",true);     // enable dataset mode
   
    cout << "set eventbuilder to " << eventbuilder << endl;
 
    ///----------------------------------------------------------------
    /// Read information used in MC weighting, multi-dataset jobs
    ///----------------------------------------------------------------
    ReadDatasetInfo(dataset  , confProofAna  ); 
    confProofAna->Write();
    options->Close();
    delete options;
 
 
    cout << "All setup, ready to go " << endl; 
    int runNevents=1; 
 

    // Decide to run local or on the cluster
    if(mode.CompareTo("local")==0) runLocal(dataset,treename,nentries);
    else{
        runProof(url,dataset,-1,treename);
    }
    gSystem->Unlink("options.root");

}


const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y%m%d.%H.%M", &tstruct);
    return buf;
}
