#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/baselineDef.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "TStopwatch.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "CS.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"

using namespace std;

// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 5 arguments "<<"SubsampleName"<<" Input Template" <<" MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Closure TTbarInc 1000 0 1" << std::endl;
      return -1;
    }
  const char *subsamplename = argv[1];
  const char *Maxevent = argv[2];
  const  char *Stratfile = argv[3];
  const  char *Filerun = argv[4];
  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);
  bool isData = false;
  TChain *fChain = 0;
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(subsamplename, startfile);
  const string condorSpec = argc==6 ? argv[5]: "";

  if(!FillChain(fChain, subsamplename, condorSpec, startfile, filerun))
       {
	 std::cerr << "Cannot get the tree " << std::endl;
       }
  const int maxevent = std::atoi(Maxevent);  
  
  TString sampleString(subsamplename);
  if(sampleString.Contains("Data")){Lumiscale = 1.0; isData = true;}
  if(!isData)
    {
      std::cout<<"!..................!"<<std::endl;
      std::cerr<<"Please run ./MC for this sample"<<std::endl;
      std::cout<<"!..................!"<<std::endl;
    } 
 //Searchbin                                                                                                            
  SearchBins SB("SB_69_2016");
  //Use BaselineVessel class for baseline variables and selections
  std::string spec = "Data";
  //  std::string filterevent = "SingleMuon_csc2015.txt";
  ExpBaselineVessel = new BaselineVessel(spec);
  AnaFunctions::prepareForNtupleReader();
  AnaFunctions::prepareTopTagger();
  NTupleReader *tr =0;
  
  tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames_DataOnly);
  tr->registerFunction(&passBaselineFuncExp);
// Add cleanJets function
  stopFunctions::cjh.setMuonIso("mini");
  stopFunctions::cjh.setElecIso("mini");
  stopFunctions::cjh.setRemove(false);
  tr->registerFunction(&stopFunctions::cleanJets);
 
  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<Lumiscale<<std::endl; 
  cout<<"maxevent: "<<maxevent<<endl;
  // Loop over the events (tree entries)
  while(tr->getNextEvent()){
    if(maxevent>=0 && tr->getEvtNum() > maxevent ) break;
    // Add print out of the progress of looping
    if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;

    const vector<TLorentzVector> &muonsLVec = tr->getVec<TLorentzVector>("muonsLVec");
    const vector<double> &muonsRelIso = tr->getVec<double>("muonsRelIso");
    const vector<double> &muonsMiniIso = tr->getVec<double>("muonsMiniIso");
    const vector<double> &muonsMtw = tr->getVec<double>("muonsMtw");    
    const vector<TLorentzVector> &jetsLVec = tr->getVec<TLorentzVector>("jetsLVec");
    const vector<double> &recoJetsBtag_0 = tr->getVec<double>("recoJetsBtag_0");
    const vector<int> &looseisoTrksMatchedJetIdx = tr->getVec<int>("looseisoTrksMatchedJetIdx");
    const vector<TLorentzVector> &loose_isoTrksLVec = tr->getVec<TLorentzVector>("loose_isoTrksLVec");
    const vector<double> &loose_isoTrks_iso = tr->getVec<double>("loose_isoTrks_iso");
    const vector<double> &loose_isoTrks_mtw = tr->getVec<double>("loose_isoTrks_mtw");
    const vector<int> &loose_isoTrks_pdgId = tr->getVec<int>("loose_isoTrks_pdgId");
    const std::vector<double> muonspfActivity = tr->getVec<double>("muonspfActivity");
    const std::vector<int> & muonsFlagIDVec = tr->getVec<int>("muonsFlagMedium");
    const std::vector<double>& recoJetschargedHadronEnergyFraction = tr->getVec<double>("recoJetschargedHadronEnergyFraction");
    const std::vector<double>& recoJetschargedEmEnergyFraction = tr->getVec<double>("recoJetschargedEmEnergyFraction");
    const std::vector<std::string> & TriggerNames = tr->getVec<std::string>("TriggerNames");
    const std::vector<int> & PassTrigger = tr->getVec<int>("PassTrigger");

    double met=tr->getVar<double>("met");
    double metphi=tr->getVar<double>("metphi");

    TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);
    int run = tr->getVar<int>("run");
    int lumi = tr->getVar<int>("lumi");
    int event = tr->getVar<int>("event");
    int nElectrons = tr->getVar<int>("nElectrons_CUT"+spec);
    int nMuons = tr->getVar<int>("nMuons_CUT"+spec);
    bool passNoiseEventFilter = tr->getVar<bool>("passNoiseEventFilter"+spec);
    const double EvtWt = tr->getVar<double>("evtWeight");
    //change event weight for MC sample
    EventWeight = EvtWt;
    Lumiscale = Lumiscale * EventWeight;
    bool passBaseline_tru = tr->getVar<bool>("passBaseline"+spec);
    bool passMuonVeto_tru = tr->getVar<bool>("passMuonVeto"+spec);
    bool passEleVeto_tru = tr->getVar<bool>("passEleVeto"+spec);
    bool passIsoTrkVeto_tru = tr->getVar<bool>("passIsoTrkVeto"+spec);
    bool passLeptVeto_tru = tr->getVar<bool>("passLeptVeto"+spec);
    bool passnJets_tru = tr->getVar<bool>("passnJets"+spec);
    bool passdPhis_tru = tr->getVar<bool>("passdPhis"+spec);
    bool passMET_tru =  tr->getVar<bool>("passMET"+spec);
    bool passBJets_tru = tr->getVar<bool>("passBJets"+spec);
    bool passTagger_tru = tr->getVar<bool>("passTagger"+spec);
    bool passHT_tru =  tr->getVar<bool>("passHT"+spec);
    bool passMT2_tru = tr->getVar<bool>("passMT2" + spec);
    const int nJets_tru = tr->getVar<int>("cntNJetsPt30Eta24"+spec);
    const int nbJets_tru = tr->getVar<int>("cntCSVS"+spec);
    const int nTops_tru = tr->getVar<int>("nTopCandSortedCnt"+spec);
    const double MT2_tru = tr->getVar<double>("best_had_brJet_MT2"+spec);
    const vector<double> dPhiVec_tru = tr->getVec<double>("dPhiVec"+spec);
    const double ht_tru = tr->getVar<double>("HT"+spec);
    
      //mht calculation
      TLorentzVector Mht_truLVec;	      
	for(unsigned int ij=0; ij<jetsLVec.size(); ij++){
	  if( !AnaFunctions::jetPassCuts(jetsLVec[ij], AnaConsts::pt30Arr) ) continue;
	  Mht_truLVec -= jetsLVec[ij];
	}
	const double Mht_tru = Mht_truLVec.Pt();

      bool passBaselineFull = passMuonVeto_tru && passEleVeto_tru && passIsoTrkVeto_tru && passnJets_tru && passdPhis_tru && passMET_tru && passBJets_tru && passTagger_tru && passHT_tru && passMT2_tru;
      bool passTrigger = true; 
      if(isData){
	bool foundTrigger = false;
	for(unsigned it=0; it<TriggerNames.size(); it++){
	  if( sampleString.Contains("SingleMuon") ){
	    if( TriggerNames[it].find("HLT_PFHT300_PFMET100_v") != string::npos){
	      if( PassTrigger[it] ) foundTrigger = true;
	    }
	  }
	}
	if( !foundTrigger ) passTrigger = false;
      }
      
      // Dist.
      if(passBaselineFull && passNoiseEventFilter && passTrigger){
	int jSR = SB.find_Binning_Index(nbJets_tru, nTops_tru, MT2_tru, met);
	if( jSR!= -1 ) {
	  myBaseHistgram.hYields->Fill(jSR, Lumiscale);
	}

	FillDouble(myBaseHistgram.hMET, met, Lumiscale);
	FillDouble(myBaseHistgram.hMT2, MT2_tru, Lumiscale);
	FillInt(myBaseHistgram.hNbJets, nbJets_tru, Lumiscale);
	FillInt(myBaseHistgram.hNTops, nTops_tru, Lumiscale);	

	FillInt(myBaseHistgram.hNJets, nJets_tru, Lumiscale);
	FillDouble(myBaseHistgram.hHT, ht_tru, Lumiscale);
	FillDouble(myBaseHistgram.hdPhi0, dPhiVec_tru[0], Lumiscale);
	FillDouble(myBaseHistgram.hdPhi1, dPhiVec_tru[1], Lumiscale);
	FillDouble(myBaseHistgram.hdPhi2, dPhiVec_tru[2], Lumiscale);
      }

  }	//event loop
  // --- Save the Histograms to File -----------------------------------
  (myBaseHistgram.oFile)->Write();
  
  return 0;
  
}


