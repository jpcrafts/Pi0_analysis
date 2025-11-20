// ClusT_Comparison.cxx
//
// This program reads a slimmed ROOT file containing NPS replay data with the following branches:
//   - NPS.cal.clusT      : Raw cluster times
//   - NPS.cal.newEWClusT : Energy-weighted new cluster times
//   - NPS.cal.newClusT   : Simple average new cluster times
//   - NPS.cal.nclust     : Number of clusters
//   - HMS-level branches: 
//         T.hms.hEDTM_tdcTimeRaw, H.gtr.dp, H.cal.etotnorm, 
//         H.cer.npeSum, H.gtr.th, H.gtr.ph
//
// It applies the basic HMS cuts and then fills histograms for:
//   1. Raw cluster times
//   2. Energy-weighted new cluster times
//   3. Simple average new cluster times
//   4. Difference: (Energy-weighted new - raw)
//   5. Difference: (Simple average new - raw)
//   6. A text pad displaying the HMS cut criteria
//
// The canvas is divided into 6 pads (2 rows x 3 columns) and saved as a PNG in the "Plots" directory.
//
// Usage (from the command line):
//   ./ClusT_Comparison <input_root_file>
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
const double gtrphCut    = 0.09;

void compareClusT(const char* inputFileName, const char* outputPNG) {
  // Open the input file.
  TFile *fin = TFile::Open(inputFileName, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "Error: Cannot open file " << inputFileName << endl;
    return;
  }
  
  // Get the TTree named "T".
  TTree *tree = (TTree*) fin->Get("T");
  if (!tree) {
    cout << "Error: TTree 'T' not found in file " << inputFileName << endl;
    fin->Close();
    return;
  }
  
  // Set branch addresses.
  // Raw cluster times.
  const int maxClusters = 100;
  Double_t clusT_raw[maxClusters];
  tree->SetBranchAddress("NPS.cal.clusT", clusT_raw);
  
  // Number of clusters.
  Double_t nclust_d;
  tree->SetBranchAddress("NPS.cal.nclust", &nclust_d);
  
  // HMS-level branches.
  Double_t edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph;
  tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
  tree->SetBranchAddress("H.gtr.dp", &hdelta);
  tree->SetBranchAddress("H.cal.etotnorm", &hcaltot);
  tree->SetBranchAddress("H.cer.npeSum", &hcernpe);
  tree->SetBranchAddress("H.gtr.th", &gtrth);
  tree->SetBranchAddress("H.gtr.ph", &gtrph);
  
  // New cluster time branches.
  std::vector<double>* vecNewEWClusT = nullptr;
  tree->SetBranchAddress("NPS.cal.newEWClusT", &vecNewEWClusT);
  
  std::vector<double>* vecNewClusT = nullptr;
  tree->SetBranchAddress("NPS.cal.newClusT", &vecNewClusT);
  
  // Create histograms.
  TH1F* hRawClusterTime = new TH1F("hRawClusterTime", "Raw Cluster Times;Time (ns);Counts", 200, 100, 200);
  TH1F* hNewEWClusT     = new TH1F("hNewEWClusT", "Energy-Weighted New Cluster Times;Time (ns);Counts", 200, 100, 200);
  TH1F* hNewClusT       = new TH1F("hNewClusT", "Simple Average New Cluster Times;Time (ns);Counts", 200, 100, 200);
  TH1F* hDiffEW         = new TH1F("hDiffEW", "Difference (EW - Raw);#Delta Time (ns);Counts", 100, -5, 5);
  TH1F* hDiffSimple     = new TH1F("hDiffSimple", "Difference (Simple - Raw);#Delta Time (ns);Counts", 100, -5, 5);
  
  Long64_t nEntries = tree->GetEntries();
  cout << "Total entries in tree: " << nEntries << endl;
  
  // Loop over events.
  for (Long64_t ievt = 0; ievt < nEntries; ievt++){
    tree->GetEntry(ievt);
    
    // Apply HMS cuts.
    if (edtmtdc >= edtmtdcCut || hdelta <= hdeltaLow || hdelta >= hdeltaHigh ||
        hcaltot <= hcaltotCut || hcernpe <= hcernpeCut)
      continue;
    if (fabs(gtrth) > gtrthCut || fabs(gtrph) > gtrphCut)
      continue;
    
    int nclust = static_cast<int>(nclust_d);
    if(nclust > maxClusters) nclust = maxClusters;
    
    for (int i = 0; i < nclust; i++){
      // Fill raw cluster times histogram.
      hRawClusterTime->Fill(clusT_raw[i]);
      
      // Fill energy-weighted new cluster times and difference.
      if(vecNewEWClusT && i < (int)vecNewEWClusT->size()){
         double newEW = (*vecNewEWClusT)[i];
         hNewEWClusT->Fill(newEW);
         hDiffEW->Fill(newEW - clusT_raw[i]);
      }
      
      // Fill simple average new cluster times and difference.
      if(vecNewClusT && i < (int)vecNewClusT->size()){
         double newSimple = (*vecNewClusT)[i];
         hNewClusT->Fill(newSimple);
         hDiffSimple->Fill(newSimple - clusT_raw[i]);
      }
    }
    
    if ((ievt+1) % 10000 == 0)
      cout << "Processed " << (ievt+1) << " events." << endl;
  }
  
  // Create a canvas with 6 pads (2 rows x 3 columns).
  TCanvas* c = new TCanvas("c", "Timing Histograms", 1200, 800);
  c->Divide(3,2);
  
  // Pad 1: Raw cluster times.
  c->cd(1);
  hRawClusterTime->Draw();
  
  // Pad 2: Energy-weighted new cluster times.
  c->cd(2);
  hNewEWClusT->Draw();
  
  // Pad 3: Simple average new cluster times.
  c->cd(3);
  hNewClusT->Draw();
  
  // Pad 4: Difference plot for energy-weighted times.
  c->cd(4);
  hDiffEW->Draw();
  
  // Pad 5: Difference plot for simple average times.
  c->cd(5);
  hDiffSimple->Draw();
  
  // Pad 6: Display HMS cut information.
  c->cd(6);
  TPaveText *ptCuts = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
  ptCuts->AddText("HMS Cuts:");
  ptCuts->AddText(Form("T.hms.hEDTM_tdcTimeRaw < %.2f", edtmtdcCut));
  ptCuts->AddText(Form("H.gtr.dp between %.1f and %.1f", hdeltaLow, hdeltaHigh));
  ptCuts->AddText(Form("H.cal.etotnorm > %.2f", hcaltotCut));
  ptCuts->AddText(Form("H.cer.npeSum > %.1f", hcernpeCut));
  ptCuts->AddText(Form("H.gtr.th within ±%.2f", gtrthCut));
  ptCuts->AddText(Form("H.gtr.ph within ±%.2f", gtrphCut));
  ptCuts->SetFillColor(0);
  ptCuts->SetBorderSize(0);
  ptCuts->Draw();
  
  c->Update();
  
  // Create output directory "Plots" if it doesn't exist.
  system("mkdir -p Plots");
  
  // Extract run number from inputFileName.
  TString runStr(inputFileName);
  Ssiz_t pos = runStr.Last('/');
  if(pos != kNPOS)
    runStr.Remove(0, pos+1);
  runStr.ReplaceAll("_newClusT.root", "");
  
  TString outPNG = Form("Plots/%s_ClusT_comparison.png", runStr.Data());
  c->SaveAs(outPNG.Data());
  
  cout << "PNG file saved as: " << outPNG.Data() << endl;
  
  fin->Close();
  delete fin;
}

#ifndef __CINT__
int main(int argc, char* argv[]){
  if(argc < 2){
    cout << "Usage: " << argv[0] << " <input_root_file>" << endl;
    return 1;
  }
  compareClusT(argv[1], "output.png");
  return 0;
}
#endif
