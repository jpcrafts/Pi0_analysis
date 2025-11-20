// fitTimingDistribution.cxx
//
// This macro reads a slimmed ROOT file containing NPS replay data,
// applies basic HMS cuts, fills a histogram of energy‐weighted cluster times
// (from branch "NPS.cal.newEWClusT"), and fits the background in a sideband region
// (110–140 ns) with a function defined as:
//
//   B(x) = [0] + sum_{n=-N}^{N} [1]*exp(-0.5*((x - ([2] + n*[4]))/[3])^2)
// where:
//   [0] = constant background (C),
//   [1] = amplitude of each Gaussian (A),
//   [2] = center of the first peak (x0),
//   [3] = width (sigma),
//   [4] = period (T, ~2 ns),
//   [5] = number of peaks on each side (N)
//
// HMS cuts applied (from HMS-level branches):
//   T.hms.hEDTM_tdcTimeRaw < 0.1,
//   H.gtr.dp between -8.5 and 8.5,
//   H.cal.etotnorm > 0.6,
//   H.cer.npeSum > 1.0,
//   |H.gtr.th| < 0.09,
//   |H.gtr.ph| < 0.09
//
// The histogram is rebinned with 400 bins (from 100 to 200 ns) and its display is zoomed
// to 105–145 ns. The canvas is set to 1600×800. Two TPaveText boxes are added on the canvas:
//   - One (bottom right) shows the number of surviving entries (only entries that pass all cuts
//     and lie in the range 105–145 ns).
//   - One (bottom left) shows the HMS cut criteria.
// The output PNG is saved to "/Plots/bkgnd_plots/".
//
// Usage:
//    ./fitTimingDistribution <input_root_file> [output_png]
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
#include <vector>

using namespace std;

// HMS cut constants.
const double edtmtdcCut = 0.1;
const double hdeltaLow = -8.5;
const double hdeltaHigh = 8.5;
const double hcaltotCut = 0.6;
const double hcernpeCut = 1.0;
const double gtrthCut = 0.09;
const double gtrphCut = 0.09;

// Define the new background function.
// par[0] = C (constant background)
// par[1] = A (amplitude of each Gaussian)
// par[2] = x0 (center of the first peak)
// par[3] = sigma (width)
// par[4] = T (period, ~2 ns)
// par[5] = N (number of peaks on each side)
Double_t bgFunc(Double_t *x, Double_t *par)
{
  double C = par[0];
  double A = par[1];
  double x0 = par[2];
  double sigma = par[3];
  double T = par[4];
  int N = int(par[5]);

  double sum = C;
  for (int n = 0; n <= N; n++)
  {
    sum += A * exp(-0.5 * pow((x[0] - (x0 + n * T)) / sigma, 2));
  }
  return sum;
}

void fitTimingDistribution(const char *inputFileName, const char *outputPNG)
{
  // Open the input file.
  TFile *fin = TFile::Open(inputFileName, "READ");
  if (!fin || fin->IsZombie())
  {
    cout << "Error: Cannot open file " << inputFileName << endl;
    return;
  }

  // Get the TTree named "T".
  TTree *tree = (TTree *)fin->Get("T");
  if (!tree)
  {
    cout << "Error: TTree 'T' not found in file " << inputFileName << endl;
    fin->Close();
    return;
  }

  // Set branch addresses for energy-weighted cluster times and number of clusters.
  std::vector<double> *vecNewEWClusT = nullptr;
  tree->SetBranchAddress("NPS.cal.newEWClusT", &vecNewEWClusT);

  Double_t nclust_d;
  tree->SetBranchAddress("NPS.cal.nclust", &nclust_d);

  // Set branch addresses for HMS-level branches.
  Double_t edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph;
  tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
  tree->SetBranchAddress("H.gtr.dp", &hdelta);
  tree->SetBranchAddress("H.cal.etotnorm", &hcaltot);
  tree->SetBranchAddress("H.cer.npeSum", &hcernpe);
  tree->SetBranchAddress("H.gtr.th", &gtrth);
  tree->SetBranchAddress("H.gtr.ph", &gtrph);

  // Create a histogram for the energy-weighted timing distribution with 400 bins.
  TH1F *hTiming = new TH1F("hTiming", "Energy-Weighted Cluster Timing Distribution;Cluster Time (ns);Counts", 650, 100, 200);

  Long64_t nEntries = tree->GetEntries();
  cout << "Total entries in tree: " << nEntries << endl;

  int nSurviving = 0;
  // Loop over events.
  for (Long64_t ievt = 0; ievt < nEntries; ievt++)
  {
    tree->GetEntry(ievt);

    // Apply HMS cuts.
    if (edtmtdc >= edtmtdcCut || hdelta <= hdeltaLow || hdelta >= hdeltaHigh ||
        hcaltot <= hcaltotCut || hcernpe <= hcernpeCut)
      continue;
    if (fabs(gtrth) > gtrthCut || fabs(gtrph) > gtrphCut)
      continue;

    int nclust = static_cast<int>(nclust_d);
    for (int i = 0; i < nclust; i++)
    {
      if (vecNewEWClusT && i < (int)vecNewEWClusT->size())
      {
        double clusTime = (*vecNewEWClusT)[i];
        // Only count and fill if in the range 105-145 ns.
        if (clusTime >= 105 && clusTime <= 145)
        {
          hTiming->Fill(clusTime);
          nSurviving++;
        }
      }
    }

    if ((ievt + 1) % 10000 == 0)
      cout << "Processed " << (ievt + 1) << " events." << endl;
  }

  // Define the background fit function using the sum-of-Gaussians model.
  // We'll fit the sideband region from 110 to 140 ns.
  TF1 *fBg = new TF1("fBg", bgFunc, 110, 140, 6);
  // Initial parameters: C = 4000, A = 10000, x0 = 112 ns, sigma = 0.5 ns, T = 2.0 ns, N = 3.
  fBg->SetParameters(3500, 7000, 113.5, 3, 2.0, 0);
  fBg->SetParLimits(0, 3000, 3750);   // C: background level between 0 and 1e6
  //fBg->SetParLimits(1, 6000, 8000);   // A: amplitude between 0 and 1e6 (adjust as needed)
  //fBg->SetParLimits(2, 112, 114);     // x0: first peak must be in [110, 140] ns
  //fBg->SetParLimits(3, 0.1, 10);      // sigma: peak width between 0.1 and 10 ns
  //fBg->SetParLimits(4, 1.5, 2.5);     // T: period between 1.5 and 2.5 ns (if you expect ~2 ns)
  //fBg->SetParLimits(5, 10, 15);       // N: number of peaks from 1 to 10 (integer, but stored as double)

  // Fit the histogram in the sideband region.
  hTiming->Fit(fBg, "R");

  // Create a canvas that is 1600x800.
  TCanvas *c = new TCanvas("c", "Timing Fit", 1600, 800);

  // Set the x-axis range to 105-145 ns.
  hTiming->GetXaxis()->SetRangeUser(105, 145);
  // Draw the histogram and the background fit.
  hTiming->Draw();
  fBg->Draw("same");

  // Add a TPaveText for the surviving entries (bottom right).
  TPaveText *ptEntries = new TPaveText(0.65, 0.02, 0.95, 0.10, "NDC");
  ptEntries->AddText(Form("Entries after cuts: %d", nSurviving));
  ptEntries->SetFillColor(0);
  ptEntries->SetBorderSize(0);
  // ptEntries->Draw("same");

  // Add a TPaveText for HMS cut info (bottom left).
  // Add a TPaveText for HMS cut info (bottom left) with the previous coordinates.
  TPaveText *ptCuts = new TPaveText(0.05, 0.005, 0.95, 0.06, "NDC");
  ptCuts->AddText(Form("HMS Cuts: T.hms.hEDTM_tdcTimeRaw < %.2f, H.gtr.dp: [%.1f, %.1f], H.cal.etotnorm > %.2f, H.cer.npeSum > %.1f, |H.gtr.th| < %.2f, |H.gtr.ph| < %.2f",
                       edtmtdcCut, hdeltaLow, hdeltaHigh, hcaltotCut, hcernpeCut, gtrthCut, gtrphCut));
  ptCuts->SetFillColor(0);
  ptCuts->SetBorderSize(0);
  ptCuts->SetTextSize(0.03); // Adjust font size as needed.
  ptCuts->Draw("same");

  c->Update();

  // Create output directory "/Plots/bkgnd_plots" if needed.
  system("mkdir -p Plots/bkgnd_plots");

  // Generate output PNG filename from the input file name.
  TString runStr(inputFileName);
  Ssiz_t pos = runStr.Last('/');
  if (pos != kNPOS)
    runStr.Remove(0, pos + 1);
  runStr.ReplaceAll(".root", "");

  TString outPNG = Form("Plots/bkgnd_plots/%s_timingFit.png", runStr.Data());
  c->SaveAs(outPNG.Data());

  cout << "Saved PNG file: " << outPNG.Data() << endl;

  fin->Close();
  delete fin;
}

#ifndef __CINT__
int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    cout << "Usage: " << argv[0] << " <input_root_file> [output_png]" << endl;
    return 1;
  }
  const char *inputFileName = argv[1];
  TString outPNG;
  if (argc >= 3)
  {
    outPNG = argv[2];
  }
  else
  {
    TString runStr(inputFileName);
    Ssiz_t pos = runStr.Last('/');
    if (pos != kNPOS)
      runStr.Remove(0, pos + 1);
    runStr.ReplaceAll(".root", "");
    outPNG = Form("Plots/bkgnd_plots/%s_timingFit.png", runStr.Data());
  }
  fitTimingDistribution(inputFileName, outPNG.Data());
  return 0;
}
#endif
