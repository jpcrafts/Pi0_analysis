// ClusT_Comparison.cxx
//
// This program reads a slimmed ROOT file containing NPS replay data (with branch "NPS.cal.newClusT")
// and generates timing histograms:
//   - hRawClusterTime: raw cluster times (from branch "NPS.cal.clusT")
//   - hNewClusterTime: new energy-weighted cluster times (from branch "NPS.cal.newClusT")
//   - hDiff: the difference (new - raw)
// It applies basic HMS cuts on the following branches:
//   T.hms.hEDTM_tdcTimeRaw, H.gtr.dp, H.cal.etotnorm, H.cer.npeSum, H.gtr.th
// and then saves a canvas as a PNG file in a subdirectory called "Plots" with a filename
// that includes the run number (e.g., "4196_ClusT_comparison.png").
// The 4th pad displays a summary of the HMS cut information.
//
// Usage (from the command line):
//   ./ClusT_Comparison <input_root_file>
//
// Example:
//   ./ClusT_Comparison /my/output/dir/4196_newClusT.root

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TString.h"
#include "TPad.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

// HMS cut constants.
const double edtmtdcCut = 0.1;
const double hdeltaLow   = -8.5;
const double hdeltaHigh  = 8.5;
const double hcaltotCut  = 0.6;
const double hcernpeCut  = 1.0;
const double gtrthCut    = 0.09;

int main(int argc, char* argv[]){
  if(argc < 2){
    cout << "Usage: " << argv[0] << " <input_root_file>" << endl;
    return 1;
  }
  
  TString inputFile = argv[1];
  TFile* fin = TFile::Open(inputFile, "READ");
  if(!fin || fin->IsZombie()){
    cout << "Error opening file " << inputFile << endl;
    return 1;
  }
  
  // Get the slimmed tree.
  TTree* tree = (TTree*) fin->Get("T");
  if(!tree){
    cout << "Error: TTree 'T' not found in file " << inputFile << endl;
    fin->Close();
    return 1;
  }
  
  // Define histograms.
  TH1F* hRawClusterTime = new TH1F("hRawClusterTime", "Raw Cluster Times;Time (ns);Counts", 200, 100, 200);
  TH1F* hNewClusterTime = new TH1F("hNewClusterTime", "New Cluster Times;Time (ns);Counts", 200, 100, 200);
  TH1F* hDiff = new TH1F("hDiff", "New - Raw Cluster Time;Time Difference (ns);Counts", 100, -5, 5);
  
  // Set branch addresses.
  // Assume maximum clusters per event is 100.
  const int maxClusters = 100;
  Double_t clusT_raw[maxClusters];
  tree->SetBranchAddress("NPS.cal.clusT", clusT_raw);
  
  // The new cluster times are stored as a vector<double>
  std::vector<double>* newClusT = nullptr;
  tree->SetBranchAddress("NPS.cal.newClusT", &newClusT);
  
  // Also get number of clusters from branch "NPS.cal.nclust"
  Double_t nclust_d;
  tree->SetBranchAddress("NPS.cal.nclust", &nclust_d);
  
  // Set branch addresses for HMS cuts.
  Double_t edtmtdc;
  Double_t hdelta;
  Double_t hcaltot;
  Double_t hcernpe;
  Double_t gtrth;
  tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
  tree->SetBranchAddress("H.gtr.dp", &hdelta);
  tree->SetBranchAddress("H.cal.etotnorm", &hcaltot);
  tree->SetBranchAddress("H.cer.npeSum", &hcernpe);
  tree->SetBranchAddress("H.gtr.th", &gtrth);
  
  Long64_t nEntries = tree->GetEntries();
  cout << "Total entries in tree: " << nEntries << endl;
  
  // Loop over events.
  for(Long64_t ievt = 0; ievt < nEntries; ievt++){
    tree->GetEntry(ievt);
    
    // Apply HMS cuts.
    if(edtmtdc >= edtmtdcCut || hdelta <= hdeltaLow || hdelta >= hdeltaHigh ||
       hcaltot <= hcaltotCut || hcernpe <= hcernpeCut)
      continue;
    if(fabs(gtrth) > gtrthCut)
      continue;
    
    int nclust = static_cast<int>(nclust_d);
    if(nclust > maxClusters) nclust = maxClusters; // safety check
    
    for (int i = 0; i < nclust; i++){
      hRawClusterTime->Fill(clusT_raw[i]);
      if(newClusT && i < (int)newClusT->size()){
        double newT = (*newClusT)[i];
        hNewClusterTime->Fill(newT);
        hDiff->Fill(newT - clusT_raw[i]);
      }
    }
    
    if ((ievt+1) % 10000 == 0)
      cout << "Processed " << (ievt+1) << " events." << endl;
  }
  
  // Create a canvas and draw histograms.
  TCanvas* c = new TCanvas("c", "Timing Histograms", 1200, 800);
  c->Divide(2,2);
  
  c->cd(1);
  hRawClusterTime->Draw();
  
  c->cd(2);
  hNewClusterTime->Draw();
  
  c->cd(3);
  hDiff->Draw();
  
  // Pad 4: Display HMS cut information.
  c->cd(4);
  TPaveText *ptCuts = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
  ptCuts->AddText("HMS Cuts:");
  ptCuts->AddText(Form("T.hms.hEDTM_tdcTimeRaw < %.2f", edtmtdcCut));
  ptCuts->AddText(Form("H.gtr.dp between %.1f and %.1f", hdeltaLow, hdeltaHigh));
  ptCuts->AddText(Form("H.cal.etotnorm > %.2f", hcaltotCut));
  ptCuts->AddText(Form("H.cer.npeSum > %.1f", hcernpeCut));
  ptCuts->AddText(Form("H.gtr.th within ±%.2f", gtrthCut));
  ptCuts->SetFillColor(0);
  ptCuts->SetBorderSize(0);
  ptCuts->Draw();
  
  c->Update();
  
  // Create output directory "Plots" if it doesn't exist.
  system("mkdir -p Plots");
  
  // Extract run number from input file name.
  TString runStr = inputFile;
  Ssiz_t pos = runStr.Last('/');
  if(pos != kNPOS)
    runStr.Remove(0, pos+1);
  runStr.ReplaceAll("_newClusT.root", "");
  
  // Construct output PNG file name.
  TString outPNG = Form("Plots/%s_ClusT_comparison.png", runStr.Data());
  c->SaveAs(outPNG.Data());
  
  cout << "PNG file saved as: " << outPNG.Data() << endl;
  
  fin->Close();
  delete fin;
  
  return 0;
}
