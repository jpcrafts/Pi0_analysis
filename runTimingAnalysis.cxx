// File: runTimingAnalysis.cpp
// Description:
//   1) Reads a ROOT file from a specified directory,
//   2) Applies HMS cuts and fills two 1D timing histograms using TTreeReader,
//      and also fills two 2D maps using the detector geometry (36 rows x 30 columns):
//         - A Block Contribution Map (TH2F) and an Alternative Block Deviation Map (TH2F),
//           (these 2D maps are disabled in this version).
//   3) Additionally, the script fits:
//         - The Raw Cluster Time histogram (background subtracted) with a pure Gaussian,
//         - The Good Cluster Block Time histogram with a Gaussian plus constant background,
//         - The In-Time ADC-TDC Diff Time histogram (hIntime) with a Gaussian plus constant background,
//         - The In-Time Cluster Time histogram (hIntimeClusters) now filled using block-level ADC-TDC diff times,
//           but only from the two clusters whose clusT values are closest to 155.5 ns.
//   4) A combined canvas is drawn for four 1D histograms arranged in a 2×2 grid.
//   5) For the "intime" selection, an event must pass all HMS cuts and have at least two clusters
//      with clusE > 0.6 and clusT in [150,160] ns such that at least one pair is within ±3 ns.
//      For such events, block ADC-TDC diff times fill hIntime,
//      and hIntimeClusters is filled only for blocks associated with the two clusters closest to 155.5 ns.
//   6) All histograms are defined with a range from 100 to 200 ns.
//   
// Compile with:
//   g++ -O2 -o runTimingAnalysis runTimingAnalysis.cpp `root-config --cflags --libs`
// Run with:
//   ./runTimingAnalysis <nrun>
//
// Note: Requires ROOT 6 or later.

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TPaveText.h" // Used only for pads 1 and 2
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TColor.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <limits>
#include <algorithm> // For std::sort and fabs

using namespace std;

static const int ROWS = 36;
static const int COLS = 30;
static const int N_BLOCKS = ROWS * COLS;

int main(int argc, char *argv[])
{
   if (argc < 2)
   {
      cout << "Usage: " << argv[0] << " <nrun>" << endl;
      return 1;
   }

   // Enable full fit stats.
   gStyle->SetOptFit(1111);

   int nrun = atoi(argv[1]);
   TString inFileName = Form("/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%d_1_-1.root", nrun);
   TFile *f = TFile::Open(inFileName, "READ");
   if (!f || f->IsZombie())
   {
      cout << "Error opening file " << inFileName << endl;
      return 1;
   }

   // Get the TTree and disable unused branches.
   TTree *tree = (TTree *)f->Get("T");
   if (!tree)
   {
      cout << "Error: TTree 'T' not found in file " << inFileName << endl;
      return 1;
   }
   tree->SetBranchStatus("*", 0);
   tree->SetBranchStatus("NPS.cal.nclust", 1);
   tree->SetBranchStatus("NPS.cal.clusT", 1);
   tree->SetBranchStatus("NPS.cal.clusE", 1);
   tree->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
   tree->SetBranchStatus("H.gtr.dp", 1);
   tree->SetBranchStatus("H.cal.etotnorm", 1);
   tree->SetBranchStatus("H.cer.npeSum", 1);
   tree->SetBranchStatus("NPS.cal.fly.goodAdcTdcDiffTime", 1);
   tree->SetBranchStatus("NPS.cal.fly.e", 1);
   tree->SetBranchStatus("NPS.cal.fly.block_clusterID", 1);

   // Create a TTreeReader.
   TTreeReader reader(tree);
   TTreeReaderValue<Double_t> nclust(reader, "NPS.cal.nclust");
   TTreeReaderArray<Double_t> clusT(reader, "NPS.cal.clusT");
   TTreeReaderArray<Double_t> clusE(reader, "NPS.cal.clusE");
   TTreeReaderValue<Double_t> edtmtdc(reader, "T.hms.hEDTM_tdcTimeRaw");
   TTreeReaderValue<Double_t> hdelta(reader, "H.gtr.dp");
   TTreeReaderValue<Double_t> hcaltot(reader, "H.cal.etotnorm");
   TTreeReaderValue<Double_t> hcernpe(reader, "H.cer.npeSum");
   TTreeReaderArray<Double_t> block_t(reader, "NPS.cal.fly.goodAdcTdcDiffTime");
   TTreeReaderArray<Double_t> block_e(reader, "NPS.cal.fly.e");
   TTreeReaderArray<Double_t> block_clusterID(reader, "NPS.cal.fly.block_clusterID");

   // Create 1D histograms with range 100 to 200 ns.
   TH1F *hClusterTime = new TH1F("hClusterTime", "Raw Cluster Time Histogram (HMS Cuts);Cluster Time (ns);Counts", 200, 100, 200);
   TH1F *hGoodBlockTime = new TH1F("hGoodBlockTime", "Good Cluster Block Time Histogram;Block Time (ns);Counts", 200, 100, 200);
   TH1F *hIntime = new TH1F("hIntime", "In-Time ADC-TDC Diff Time;ADC-TDC Diff Time (ns);Counts", 200, 100, 200);
   hIntime->SetLineColor(kBlue); // Dark blue.
   // hIntimeClusters will now be filled using block data,
   // but only from blocks associated with the two clusters whose clusT are closest to 155.5 ns.
   TH1F *hIntimeClusters = new TH1F("hIntimeClusters", "In-Time Cluster Time Histogram (selected clusters);ADC-TDC Diff Time (ns);Counts", 200, 100, 200);

   // (2D maps and hAvgDiff are disabled as before.)

   // Loop over events.
   int processedEvents = 0;
   while (reader.Next())
   {
      // Apply HMS cuts.
      if (*edtmtdc >= 0.1 || *hdelta <= -12 || *hdelta >= 12 || *hcaltot <= 0.6 || *hcernpe <= 1.0)
         continue;

      processedEvents++;
      int numClus = static_cast<int>(*nclust);
      vector<bool> goodCluster(numClus, false);

      // Fill the raw cluster time histogram and mark "good" clusters.
      for (int i = 0; i < numClus; i++)
      {
         double tVal = clusT[i];
         if (tVal > 100 && tVal < 200)
            hClusterTime->Fill(tVal);
         // Mark a cluster as good if it passes the energy cut and a time window (130-170 ns).
         if (clusE[i] > 0.6 && tVal > 130 && tVal < 170)
            goodCluster[i] = true;
      }

      // For "intime" selection: require clusters with clusE > 0.6 and clusT in [150,160] ns.
      vector<double> narrowTimes;
      narrowTimes.reserve(numClus);  // Preallocate for efficiency.
      for (int i = 0; i < numClus; i++)
      {
         if (clusE[i] > 0.6 && clusT[i] >= 150 && clusT[i] <= 160)
            narrowTimes.push_back(clusT[i]);
      }

      // Use a sorted approach to check if any two cluster times are within ±3 ns.
      bool intimeEvent = false;
      if (narrowTimes.size() >= 2)
      {
         sort(narrowTimes.begin(), narrowTimes.end());
         for (size_t i = 0; i < narrowTimes.size() - 1; i++)
         {
            if ((narrowTimes[i+1] - narrowTimes[i]) <= 3)
            {
               intimeEvent = true;
               break;
            }
         }
      }

      // Fill the 1D Good Cluster Block Time histogram from block data.
      int nBlocks_evt = block_t.GetSize();
      for (int ib = 0; ib < nBlocks_evt; ib++)
      {
         int cid = static_cast<int>(block_clusterID[ib]);
         if (cid >= 0 && cid < numClus && goodCluster[cid])
         {
            if (block_e[ib] > 0)
            {
               double bTime = block_t[ib];
               if (bTime > 100 && bTime < 200)
                  hGoodBlockTime->Fill(bTime);
            }
         }
      }

      // For the fourth histogram, select the two clusters (from those with clusE > 0.6 and clusT in [150,160])
      // that are closest to 155.5 ns.
      int bestIndex = -1;
      int secondBestIndex = -1;
      double bestDiff = std::numeric_limits<double>::max();
      double secondBestDiff = std::numeric_limits<double>::max();
      for (int i = 0; i < numClus; i++)
      {
         if (clusE[i] > 0.6 && clusT[i] >= 150 && clusT[i] <= 160)
         {
            double diff = fabs(clusT[i] - 155.5);
            if (diff < bestDiff)
            {
               secondBestDiff = bestDiff;
               secondBestIndex = bestIndex;
               bestDiff = diff;
               bestIndex = i;
            }
            else if (diff < secondBestDiff)
            {
               secondBestDiff = diff;
               secondBestIndex = i;
            }
         }
      }

      // If the event qualifies as "intime", loop over blocks to fill the histograms.
      if (intimeEvent)
      {
         for (int ib = 0; ib < nBlocks_evt; ib++)
         {
            int cid = static_cast<int>(block_clusterID[ib]);
            if (cid >= 0 && cid < numClus && goodCluster[cid])
            {
               if (block_e[ib] > 0)
               {
                  double bTime = block_t[ib];
                  if (bTime > 100 && bTime < 200)
                  {
                     // Fill hIntime with all blocks passing the criteria.
                     hIntime->Fill(bTime);
                     // For hIntimeClusters, fill only if the block is associated with one of the two chosen clusters.
                     if (cid == bestIndex || cid == secondBestIndex)
                        hIntimeClusters->Fill(bTime);
                  }
               }
            }
         }
      }
   }
   cout << "Processed " << processedEvents << " events for run " << nrun << endl;

   // ------------------------------
   // Fit the Raw Cluster Time Histogram (background subtraction method)
   // ------------------------------
   double fitRangeMinVal = 150.0;
   double fitRangeMaxVal = 160.0;

   TF1 *bgLeftFit = new TF1("bgLeft", "pol2", 117, fitRangeMinVal);
   hClusterTime->Fit(bgLeftFit, "R+");
   TF1 *bgRightFit = new TF1("bgRight", "pol2", fitRangeMaxVal, 184);
   hClusterTime->Fit(bgRightFit, "R+");

   TH1F *hClusterTimeNoBg = (TH1F *)hClusterTime->Clone("hClusterTimeNoBg");
   for (int bin = 1; bin <= hClusterTime->GetNbinsX(); bin++)
   {
      double binCenter = hClusterTime->GetBinCenter(bin);
      double origContent = hClusterTime->GetBinContent(bin);
      double bgVal = (binCenter < fitRangeMinVal) ? bgLeftFit->Eval(binCenter) : bgRightFit->Eval(binCenter);
      double newContent = origContent - bgVal;
      hClusterTimeNoBg->SetBinContent(bin, newContent);
   }

   TF1 *gausFitRaw = new TF1("gausFitRaw", "gaus", fitRangeMinVal, fitRangeMaxVal);
   hClusterTimeNoBg->Fit(gausFitRaw, "R");
   double meanFit = gausFitRaw->GetParameter(1);
   double sigmaFit = gausFitRaw->GetParameter(2);
   cout << "Raw Cluster Time Fit (background subtracted): Mean = " << meanFit
        << " ns, Sigma = " << sigmaFit << " ns" << endl;

   // ------------------------------
   // Fit the Good Cluster Block Time Histogram with a Gaussian plus constant background.
   // ------------------------------
   double fitRangeBlockMinVal = 152.0;
   double fitRangeBlockMaxVal = 158.0;
   TF1 *gausPlusBGBlockFit = new TF1("gausPlusBGBlock", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",
                                     fitRangeBlockMinVal, fitRangeBlockMaxVal);
   gausPlusBGBlockFit->SetParameters(60000, 155, 1.0, 100);
   hGoodBlockTime->Fit(gausPlusBGBlockFit, "R");
   double ampBlock = gausPlusBGBlockFit->GetParameter(0);
   double meanBlock = gausPlusBGBlockFit->GetParameter(1);
   double sigmaBlock = gausPlusBGBlockFit->GetParameter(2);
   double constBlock = gausPlusBGBlockFit->GetParameter(3);
   cout << "Good Cluster Block Time Fit (Gaussian + BG, narrower range):" << endl
        << "  Amplitude = " << ampBlock << endl
        << "  Mean      = " << meanBlock << " ns" << endl
        << "  Sigma     = " << sigmaBlock << " ns" << endl
        << "  BG const  = " << constBlock << endl;

   // ------------------------------
   // Fit the In-Time Histogram (blocks) with a Gaussian plus constant background.
   // ------------------------------
   TF1 *gausPlusBGIntime = new TF1("gausPlusBGIntime", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", 150, 160);
   gausPlusBGIntime->SetParameters(12000, 155, 1.0, 2000);
   hIntime->Fit(gausPlusBGIntime, "R");
   double meanIntime = gausPlusBGIntime->GetParameter(1);
   double sigmaIntime = gausPlusBGIntime->GetParameter(2);
   cout << "In-Time Histogram Fit: Mean = " << meanIntime << " ns, Sigma = " << sigmaIntime << " ns" << endl;

   // ------------------------------
   // Fit the In-Time Cluster Time Histogram with a Gaussian plus constant background.
   // ------------------------------
   TF1 *gausPlusBGIntimeClust = new TF1("gausPlusBGIntimeClust","gaus", 154, 157);
   gausPlusBGIntimeClust->SetParameters(1500, 155.5, 1.0);
   hIntimeClusters->Fit(gausPlusBGIntimeClust, "R");
   double meanIntimeClust = gausPlusBGIntimeClust->GetParameter(1);
   double sigmaIntimeClust = gausPlusBGIntimeClust->GetParameter(2);
   cout << "In-Time Cluster Time Histogram Fit:" << endl
        << "  Mean  = " << meanIntimeClust << " ns" << endl
        << "  Sigma = " << sigmaIntimeClust << " ns" << endl;

   TString outDir = "/volatile/hallc/nps/jpcrafts/Plots";
   gSystem->Exec(Form("mkdir -p %s", outDir.Data()));

   TCanvas *c1_all = new TCanvas("c1_all", "All 1D Histograms", 2100, 1200);
   c1_all->Divide(2, 2);

   // Pad 1: Raw Cluster Time Histogram.
   c1_all->cd(1);
   hClusterTime->SetStats(0);
   hClusterTime->SetLineColor(kBlue);
   hClusterTime->Draw("hist");
   hClusterTimeNoBg->SetLineColorAlpha(kRed, 0.0);
   hClusterTimeNoBg->Draw("hist same");
   gausFitRaw->SetLineColor(kMagenta);
   gausFitRaw->Draw("same");
   bgLeftFit->SetLineColor(kGreen);
   bgLeftFit->SetLineStyle(2);
   bgLeftFit->Draw("same");
   bgRightFit->SetLineColor(kGreen);
   bgRightFit->SetLineStyle(2);
   bgRightFit->Draw("same");
   TPaveText *pt1 = new TPaveText(0.65, 0.75, 0.90, 0.85, "NDC");
   pt1->AddText(Form("Mean = %.2f ns", meanFit));
   pt1->AddText(Form("Sigma = %.2f ns", sigmaFit));
   pt1->SetFillColor(0);
   pt1->SetBorderSize(0);
   pt1->Draw();

   // Pad 2: Good Cluster Block Time Histogram.
   c1_all->cd(2);
   hGoodBlockTime->SetStats(1);
   hGoodBlockTime->SetLineColor(kBlue);
   hGoodBlockTime->Draw("hist");
   gausPlusBGBlockFit->SetLineColor(kMagenta);
   gausPlusBGBlockFit->Draw("same");
   TPaveText *pt2 = new TPaveText(meanBlock - 1, hGoodBlockTime->GetMaximum() * 0.75,
                                  meanBlock + 1, hGoodBlockTime->GetMaximum() * 0.85, "brNDC");
   pt2->AddText(Form("Amp = %.0f", ampBlock));
   pt2->AddText(Form("Mean = %.2f ns", meanBlock));
   pt2->AddText(Form("Sigma = %.2f ns", sigmaBlock));
   pt2->AddText(Form("BG const = %.0f", constBlock));
   pt2->SetFillColor(0);
   pt2->SetBorderSize(0);
   pt2->Draw();

   // Pad 3: In-Time Histogram (blocks).
   c1_all->cd(3);
   hIntime->SetStats(1);
   hIntime->SetLineColor(kBlue);
   hIntime->Draw("hist");
   gausPlusBGIntime->SetLineColor(kMagenta);
   gausPlusBGIntime->Draw("same");

   // Pad 4: In-Time Histogram (clusters, selected clusters only).
   c1_all->cd(4);
   hIntimeClusters->SetStats(1);
   hIntimeClusters->SetLineColor(kBlue);
   hIntimeClusters->Draw("hist");
   gausPlusBGIntimeClust->SetLineColor(kMagenta);
   gausPlusBGIntimeClust->Draw("same");

   c1_all->SaveAs(Form("%s/All1DHistograms_run%d.png", outDir.Data(), nrun));

   // (2D maps drawing is disabled.)

   // Cleanup.
   f->Close();
   delete f;

   return 0;
}
