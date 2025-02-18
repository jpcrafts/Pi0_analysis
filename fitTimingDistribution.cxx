// fitTimingDistribution.cxx
//
// This macro reads a slimmed ROOT file containing NPS replay data,
// applies basic HMS cuts, fills a histogram of raw cluster times ("NPS.cal.clusT"),
// and fits it with a function that consists of a constant background plus a 
// cosine modulation (to model periodic accidental background with ~4 ns periodicity)
// plus a large Gaussian representing the main signal.
// 
// The fit function is defined as:
//   f(x) = [0] + [1]*cos(2*TMath::Pi()*(x-[2])/4) + [3]*exp(-0.5*((x-[4])/[5])^2)
// where:
//   [0] = constant background (~18000),
//   [1] = cosine modulation amplitude (~5000),
//   [2] = phase offset (initially 0),
//   [3] = signal amplitude (~160000),
//   [4] = signal mean (~150 ns),
//   [5] = signal sigma (~2 ns).
//
// HMS cuts applied (from HMS-level branches):
//   T.hms.hEDTM_tdcTimeRaw < 0.1,
//   H.gtr.dp between -8.5 and 8.5,
//   H.cal.etotnorm > 0.6,
//   H.cer.npeSum > 1.0,
//   |H.gtr.th| < 0.09
//
// The resulting canvas is saved as a PNG in the "Plots" directory.
// Usage:
//    ./fitTimingDistribution <input_root_file> <output_png>
// Example:
//    ./fitTimingDistribution /my/output/dir/4196_newClusT.root 4196_timingFit.png

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include <iostream>
#include <cstdlib>

using namespace std;

void fitTimingDistribution(const char* inputFileName, const char* outputPNG) {
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
  
  // Set branch addresses for cluster timing and HMS cuts.
  // Cluster-level.
  const int maxClusters = 100;
  Double_t clusT_raw[maxClusters];
  tree->SetBranchAddress("NPS.cal.clusT", clusT_raw);
  // Number of clusters.
  Double_t nclust_d;
  tree->SetBranchAddress("NPS.cal.nclust", &nclust_d);
  
  // HMS-level branches.
  Double_t edtmtdc, hdelta, hcaltot, hcernpe, gtrth;
  tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
  tree->SetBranchAddress("H.gtr.dp", &hdelta);
  tree->SetBranchAddress("H.cal.etotnorm", &hcaltot);
  tree->SetBranchAddress("H.cer.npeSum", &hcernpe);
  tree->SetBranchAddress("H.gtr.th", &gtrth);
  
  // Create a histogram for the timing distribution.
  TH1F *hTiming = new TH1F("hTiming", "Cluster Timing Distribution;Cluster Time (ns);Counts", 200, 100, 200);
  
  Long64_t nEntries = tree->GetEntries();
  cout << "Total entries in tree: " << nEntries << endl;
  
  // Loop over events and fill histogram only for events passing HMS cuts.
  for (Long64_t ievt = 0; ievt < nEntries; ievt++){
    tree->GetEntry(ievt);
    
    // Apply HMS cuts.
    if (edtmtdc >= 0.1 || hdelta <= -8.5 || hdelta >= 8.5 ||
        hcaltot <= 0.6 || hcernpe <= 1.0)
      continue;
    if (fabs(gtrth) > 0.09)
      continue;
    
    int nclust = static_cast<int>(nclust_d);
    if (nclust > maxClusters) nclust = maxClusters;
    for (int i = 0; i < nclust; i++){
      hTiming->Fill(clusT_raw[i]);
    }
    
    if ((ievt+1) % 10000 == 0)
      cout << "Processed " << (ievt+1) << " events." << endl;
  }
  
  // Define the combined fit function.
  // f(x) = [0] + [1]*cos(2*pi*(x-[2])/4) + [3]*exp(-0.5*((x-[4])/[5])^2)
  TF1 *fFit = new TF1("fFit",
    "[0] + [1]*cos(2*TMath::Pi()*(x-[2])/2) + [3]*exp(-0.5*((x-[4])/[5])^2)",
    110, 190);
  
  // Set initial parameter guesses.
  // [0]: constant background ~18000
  // [1]: cosine modulation amplitude ~5000
  // [2]: cosine phase offset (initially 0)
  // [3]: signal amplitude ~160000
  // [4]: signal mean ~150 ns
  // [5]: signal sigma ~2 ns
  fFit->SetParameters(18000, 5000, 0, 160000, 150, 2);
  
  // Optionally, set parameter limits.
  fFit->SetParLimits(0, 0, 1e6);   // Constant background between 0 and 1e6.
  fFit->SetParLimits(2, -10, 10);   // Phase offset within a reasonable range.
  fFit->SetParLimits(3, 0, 1e6);    // Signal amplitude between 0 and 1e6.
  fFit->SetParLimits(5, 0.1, 10);    // Signal sigma between 0.1 and 10 ns.
  
  // Fit the histogram with the combined function.
  hTiming->Fit(fFit, "R");
  
  // Create a canvas and draw histogram and fit.
  TCanvas *c = new TCanvas("c", "Timing Fit", 800, 600);
  hTiming->Draw();
  fFit->Draw("same");
  
  // Create output directory "Plots" if it doesn't exist.
  gSystem->Exec("mkdir -p Plots");
  
  // Optionally, extract run number from input file name.
  TString runStr(inputFileName);
  Ssiz_t pos = runStr.Last('/');
  if (pos != kNPOS)
    runStr.Remove(0, pos+1);
  runStr.ReplaceAll("_newClusT.root", "");
  
  TString outPNG = Form("Plots/%s_timingFit.png", runStr.Data());
  c->SaveAs(outPNG.Data());
  
  cout << "Saved PNG file: " << outPNG.Data() << endl;
  
  fin->Close();
  delete fin;
}

#ifndef __CINT__
int main(int argc, char* argv[]){
  if(argc < 3){
    cout << "Usage: " << argv[0] << " <input_root_file> <output_png>" << endl;
    return 1;
  }
  fitTimingDistribution(argv[1], argv[2]);
  return 0;
}
#endif
