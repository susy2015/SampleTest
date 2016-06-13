#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/samples.h"

static const int nSB = 69; //We use nSB serach bins depending on Nbjet, Ntop, met and MT2 value.

using namespace std;

static BaselineVessel *ExpBaselineVessel;
void passBaselineFuncExp(NTupleReader& tr)
{
  (*ExpBaselineVessel)(tr);
}
double Lumiscale = 1.0;
double EventWeight = 1.0;
class BaseHistgram
{
 public:
  void BookHistgram(const char *, const int&);
  TFile *oFile;
  TH1D *hMET;
  TH1D *hNbJets;
  TH1D *hNTops;
  TH1D *hMT2;
  
  TH1D *hNJets;
  TH1D *hHT;
  TH1D *hYields;
  TH1D *hdPhi0;
  TH1D *hdPhi1;
  TH1D *hdPhi2;

  const TString title = "LL MC";

};

void BaseHistgram::BookHistgram(const char *outFileName, const int& filerun)
{
  TString filename(outFileName);
  TString index(std::to_string(filerun));
  filename+= "_LLMC"+index+".root";
  oFile = new TFile(filename, "recreate");
 
  hMET = new TH1D("hMET",title+";met [GeV];Events",24,200.,800.);
  hMET->Sumw2();
  hNbJets = new TH1D("hNbJets",title+";N_{bjets};Events",4, 1, 5);
  hNbJets->Sumw2();
  hNTops = new TH1D("hNTops",title+";N_{tops};Events",4, 1, 5);
  hNTops->Sumw2();
  hMT2 = new TH1D("hMT2",title+";M_{T2}[GeV];Events",12,200,500);
  hMT2->Sumw2();
  hYields = new TH1D("hYields", title+";search bin;Events",nSB,0,nSB);
  hYields->Sumw2();
 
  hNJets = new TH1D("hNJets",title+";N_{jets};Events",6 ,4,10);
  hNJets->Sumw2();
  hHT = new TH1D("hHT",title+";H_{T} [GeV];Events",20,500.,1000.);
  hHT->Sumw2();  
  hdPhi0 = new TH1D("hdPhi0", title+";dPhi0;Events", 16, 0, 3.2);
  hdPhi0->Sumw2();
  hdPhi1 = new TH1D("hdPhi1", title+";dPhi1;Events", 16, 0, 3.2);
  hdPhi1->Sumw2();
  hdPhi2 = new TH1D("hdPhi2", title+";dPhi2;Events", 16, 0, 3.2);
  hdPhi2->Sumw2();
  
}


bool FillChain(TChain* &chain, const char *subsample, const string condorSpec, const int& startfile, const int& filerun){
  
  AnaSamples::SampleSet        allSamples = condorSpec.empty()? AnaSamples::SampleSet():AnaSamples::SampleSet(condorSpec);
  AnaSamples::SampleCollection allCollections(allSamples);
  bool find = false;  
  TString subsamplename(subsample);
  
  chain = new TChain(allSamples[subsample].treePath.c_str());
  if(allSamples[subsample] != allSamples.null())
    {
      allSamples[subsample].addFilesToChain(chain, startfile, filerun);
      find = true;
      Lumiscale = allSamples[subsample].getWeight();
    }
    return find;
}

void FillDouble(TH1 *hist, const double &a, const double &w){
  int nbin = hist->GetNbinsX();
  double low = hist->GetBinLowEdge(nbin);
  double high = hist->GetBinLowEdge(nbin + 1);
  double copy = a;
  if(copy >= high) copy = low;
  hist->Fill(copy, w);
}
void FillInt(TH1 *hist, const int &a, const double &w){
  int nbin = hist->GetNbinsX();
  int low = (int)hist->GetBinLowEdge(nbin);
  int high = (int)hist->GetBinLowEdge(nbin + 1);
  int copy = a;
  if(copy >= high) copy = low;
  hist->Fill(copy, w);
}
void Fill2D(TH2 *hist, const int &a, const double &b, const double &w){
  int nbinx = hist->GetNbinsX();
  int nbiny = hist->GetNbinsY();
  int lowx = hist->GetXaxis()->GetBinLowEdge(nbinx);
  int highx = hist->GetXaxis()->GetBinLowEdge(nbinx + 1);
  double lowy = hist->GetYaxis()->GetBinLowEdge(nbiny);
  double highy = hist->GetYaxis()->GetBinLowEdge(nbiny + 1);
  int copyx = a;
  if(copyx >= highx) copyx = lowx;
  double copyy = b;
  if(copyy >= highy) copyy = lowy;
  hist->Fill(copyx, copyy, w);
}

