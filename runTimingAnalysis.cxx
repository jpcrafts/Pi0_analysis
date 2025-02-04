// File: runTimingAnalysis.cpp
// Description:
//   1) Reads a ROOT file from a specified directory,
//   2) Applies HMS cuts and fills two 1D timing histograms using TTreeReader,
//      and also fills two 2D maps using the detector geometry (36 rows x 30 columns):
//         - A Block Contribution Map (TH2F): counts, for each block, the number of events in which that block contributed
//           to a "good" cluster (i.e. its block time is in [100,200] ns and its associated cluster passes the cut).
//         - An Alternative Block Deviation Map: calculated by accumulating (bTime - 155.5) values in each bin (TH2F hSumDiff)
//           and counting the number of entries (TH2F hCount). After processing, a new histogram (hAvgDiff) is computed
//           by dividing hSumDiff by hCount on a bin‐by‐bin basis. For bins with zero count, the cell is left blank (set to NaN).
//   3) Additionally, the script fits the Raw Cluster Time histogram (after background subtraction)
//      with a Gaussian (in the range 150–160 ns) and the Good Cluster Block Time histogram with a
//      Gaussian plus constant background (in [152,158] ns).
//   4) The script draws a combined canvas for the two 1D histograms and saves it as a PNG.
//   5) The Block Contribution Map is drawn on its own canvas (720×864) and saved as a PNG.
//   6) The Alternative Block Deviation Map (hAvgDiff) is drawn on its own separate canvas (720×864) using a custom
//      red-green-blue gradient. The pad is forced to have a fixed aspect ratio so that the individual bins appear square.
//      Additionally, SetStats(0) is called on both 2D maps so that their stats boxes are hidden.
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
#include "TPaveText.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TColor.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <limits>

using namespace std;

static const int ROWS = 36;
static const int COLS = 30;
static const int N_BLOCKS = ROWS * COLS;

int main(int argc, char* argv[]){
   if(argc < 2){
      cout << "Usage: " << argv[0] << " <nrun>" << endl;
      return 1;
   }
   
   // Enable full fit stats.
   gStyle->SetOptFit(1111);
   
   int nrun = atoi(argv[1]);
   TString inFileName = Form("/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%d_1_-1.root", nrun);
   TFile *f = TFile::Open(inFileName, "READ");
   if(!f || f->IsZombie()){
      cout << "Error opening file " << inFileName << endl;
      return 1;
   }
   
   // Get the TTree and disable unused branches.
   TTree *tree = (TTree*)f->Get("T");
   if(!tree){
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
   
   // Create 1D histograms (range 100–200 ns).
   TH1F *hClusterTime = new TH1F("hClusterTime", "Raw Cluster Time Histogram (HMS Cuts);Cluster Time (ns);Counts", 200, 100, 200);
   TH1F *hGoodBlockTime = new TH1F("hGoodBlockTime", "Good Cluster Block Time Histogram;Block Time (ns);Counts", 200, 100, 200);
   
   // Create 2D maps.
   TH2F *hBlockContrib = new TH2F("hBlockContrib", "Block Contribution Map;Column;Row", 
                                  COLS, 0, COLS, ROWS, 0, ROWS);
   // Alternative Block Deviation Map using two TH2F's.
   TH2F *hSumDiff = new TH2F("hSumDiff", "Sum of (t-155.5);Column;Row", 
                             COLS, 0, COLS, ROWS, 0, ROWS);
   TH2F *hCount   = new TH2F("hCount", "Count per Block;Column;Row", 
                             COLS, 0, COLS, ROWS, 0, ROWS);
   // hAvgDiff will be computed after the event loop.
   TH2F *hAvgDiff = (TH2F*)hSumDiff->Clone("hAvgDiff");
   
   // Loop over events.
   int processedEvents = 0;
   while(reader.Next()){
      // Apply HMS cuts.
      if(*edtmtdc >= 0.1 || *hdelta <= -12 || *hdelta >= 12 || *hcaltot <= 0.6 || *hcernpe <= 1.0)
         continue;
      
      processedEvents++;
      
      int numClus = static_cast<int>(*nclust);
      vector<bool> goodCluster(numClus, false);
      
      // Fill cluster time histogram and mark "good" clusters.
      for(int i=0; i<numClus; i++){
         double tVal = clusT[i];
         if(tVal > 100 && tVal < 200)
            hClusterTime->Fill(tVal);
         if(clusE[i] > 0.6 && tVal > 130 && tVal < 170)
            goodCluster[i] = true;
      }
      
      // Fill the 1D Good Cluster Block Time histogram.
      int nBlocks_evt = block_t.GetSize();
      for(int ib=0; ib<nBlocks_evt; ib++){
         int cid = static_cast<int>(block_clusterID[ib]);
         if(cid >= 0 && cid < numClus && goodCluster[cid]){
            if(block_e[ib] > 0){
               double bTime = block_t[ib];
               if(bTime > 100 && bTime < 200)
                  hGoodBlockTime->Fill(bTime);
            }
         }
      }
      
      // Update the 2D Block Contribution Map and the Alternative Block Deviation Map.
      for(int ib=0; ib < N_BLOCKS; ib++){
         int cid = static_cast<int>(block_clusterID[ib]);
         if(cid >= 0 && cid < numClus && goodCluster[cid]){
            if(block_e[ib] > 0){
               double bTime = block_t[ib];
               if(bTime > 100 && bTime < 200){
                  int row = ib / COLS;
                  int col = ib % COLS;
                  hBlockContrib->Fill(col + 0.5, row + 0.5);
                  // Accumulate the signed difference (bTime - 155.5) for this block.
                  hSumDiff->Fill(col + 0.5, row + 0.5, (bTime - 155.5));
                  hCount->Fill(col + 0.5, row + 0.5, 1);
               }
            }
         }
      }
   }
   cout << "Processed " << processedEvents << " events for run " << nrun << endl;
   
   // Compute hAvgDiff: For each bin, if hCount is zero, set bin content to NaN; else, average.
   for(int binx = 1; binx <= hAvgDiff->GetNbinsX(); binx++){
      for(int biny = 1; biny <= hAvgDiff->GetNbinsY(); biny++){
         double cnt = hCount->GetBinContent(binx, biny);
         if(cnt == 0){
            hAvgDiff->SetBinContent(binx, biny, std::numeric_limits<double>::quiet_NaN());
         } else {
            double sum = hSumDiff->GetBinContent(binx, biny);
            hAvgDiff->SetBinContent(binx, biny, sum / cnt);
         }
      }
   }
   
   // ------------------------------
   // Fit the Raw Cluster Time Histogram (background subtraction method)
   // ------------------------------
   double fitRangeMinVal = 150.0;
   double fitRangeMaxVal = 160.0;
   
   TF1 *bgLeft = new TF1("bgLeft", "pol2", 117, fitRangeMinVal);
   hClusterTime->Fit(bgLeft, "R+");
   TF1 *bgRight = new TF1("bgRight", "pol2", fitRangeMaxVal, 184);
   hClusterTime->Fit(bgRight, "R+");
   
   TH1F *hClusterTimeNoBg = (TH1F*)hClusterTime->Clone("hClusterTimeNoBg");
   for(int bin = 1; bin <= hClusterTime->GetNbinsX(); bin++){
      double binCenter = hClusterTime->GetBinCenter(bin);
      double origContent = hClusterTime->GetBinContent(bin);
      double bgVal = (binCenter < fitRangeMinVal) ? bgLeft->Eval(binCenter) : bgRight->Eval(binCenter);
      double newContent = origContent - bgVal;
      hClusterTimeNoBg->SetBinContent(bin, newContent);
   }
   
   TF1 *gausFit = new TF1("gausFit", "gaus", fitRangeMinVal, fitRangeMaxVal);
   hClusterTimeNoBg->Fit(gausFit, "R");
   double mean = gausFit->GetParameter(1);
   double sigma = gausFit->GetParameter(2);
   cout << "Raw Cluster Time Fit (background subtracted): Mean = " << mean 
        << " ns, Sigma = " << sigma << " ns" << endl;
   
   // ------------------------------
   // Fit the Good Cluster Block Time Histogram (unchanged)
   // ------------------------------
   double fitRangeBlockMin = 152.0;
   double fitRangeBlockMax = 158.0;
   TF1 *gausPlusBGBlock = new TF1("gausPlusBGBlock",
                                  "[0]*exp(-0.5*((x-[1])/[2])^2) + pol0(3)",
                                  fitRangeBlockMin, fitRangeBlockMax);
   gausPlusBGBlock->SetParameters(60000, 155, 1.0, 100);
   hGoodBlockTime->Fit(gausPlusBGBlock, "R");
   double ampBlock    = gausPlusBGBlock->GetParameter(0);
   double meanBlock   = gausPlusBGBlock->GetParameter(1);
   double sigmaBlock  = gausPlusBGBlock->GetParameter(2);
   double constBlock  = gausPlusBGBlock->GetParameter(3);
   cout << "Good Block Time Fit (Gaussian + BG, narrower range):" << endl
        << "  Amplitude = " << ampBlock << endl
        << "  Mean      = " << meanBlock << " ns" << endl
        << "  Sigma     = " << sigmaBlock << " ns" << endl
        << "  BG const  = " << constBlock << endl;
   
   // ------------------------------
   // Draw Combined 1D Histograms (unchanged)
   // ------------------------------
   TString outDir = "/volatile/hallc/nps/jpcrafts/Plots";
   gSystem->Exec(Form("mkdir -p %s", outDir.Data()));
   
   TCanvas *c1_canvas = new TCanvas("c1_canvas", "Combined Timing Histograms and Fits", 1400, 600);
   c1_canvas->Divide(2, 1);
   
   c1_canvas->cd(1);
   hClusterTime->SetStats(0);
   hClusterTime->SetLineColor(kBlue);
   hClusterTime->Draw("hist");
   hClusterTimeNoBg->SetLineColorAlpha(kRed, 0.0);
   hClusterTimeNoBg->Draw("hist same");
   gausFit->SetLineColor(kMagenta);
   gausFit->Draw("same");
   bgLeft->SetLineColor(kGreen);
   bgLeft->SetLineStyle(2);
   bgLeft->Draw("same");
   bgRight->SetLineColor(kGreen);
   bgRight->SetLineStyle(2);
   bgRight->Draw("same");
   TPaveText *pt1 = new TPaveText(0.65, 0.75, 0.90, 0.85, "NDC");
   pt1->AddText(Form("Mean = %.2f ns", mean));
   pt1->AddText(Form("Sigma = %.2f ns", sigma));
   pt1->SetFillColor(0);
   pt1->SetBorderSize(0);
   pt1->Draw();
   
   c1_canvas->cd(2);
   hGoodBlockTime->SetLineColor(kBlue);
   hGoodBlockTime->Draw("hist");
   gausPlusBGBlock->SetLineColor(kMagenta);
   gausPlusBGBlock->Draw("same");
   TPaveText *pt2 = new TPaveText(meanBlock - 1, hGoodBlockTime->GetMaximum()*0.75,
                                  meanBlock + 1, hGoodBlockTime->GetMaximum()*0.85, "brNDC");
   pt2->AddText(Form("Amp = %.0f", ampBlock));
   pt2->AddText(Form("Mean = %.2f ns", meanBlock));
   pt2->AddText(Form("Sigma = %.2f ns", sigmaBlock));
   pt2->AddText(Form("BG const = %.0f", constBlock));
   pt2->SetFillColor(0);
   pt2->SetBorderSize(0);
   pt2->Draw();
   
   c1_canvas->SaveAs(Form("%s/CombinedTimingHistograms_run%d.png", outDir.Data(), nrun));
   
   // ------------------------------
   // Draw Block Contribution Map on its own canvas.
   // ------------------------------
   TCanvas *c2_contrib = new TCanvas("c2_contrib", "Block Contribution Map", 720, 864);
   c2_contrib->SetRightMargin(0.15);
   hBlockContrib->SetContour(99);
   hBlockContrib->Draw("COLZ");
   hBlockContrib->GetXaxis()->SetTitle("Column");
   hBlockContrib->GetYaxis()->SetTitle("Row");
   gPad->SetFixedAspectRatio();
   hBlockContrib->SetStats(0);  // Hide stats for this 2D histogram.
   TPaletteAxis* pal1 = (TPaletteAxis*)hBlockContrib->GetListOfFunctions()->FindObject("palette");
   if(pal1) {
      pal1->SetX1NDC(1.0);
   }
   c2_contrib->SaveAs(Form("%s/BlockContributionMap_run%d.png", outDir.Data(), nrun));
   
   // ------------------------------
   // Draw Block Deviation Map on its own canvas.
   // ------------------------------
   TCanvas *c2_deviation = new TCanvas("c2_deviation", "Block Deviation Map", 720, 864);
   // Set up a red-green-blue gradient for the deviation map.
   const Int_t NRGBs = 3;
   Double_t stops[NRGBs] = {0.0, 0.5, 1.0};
   // For negative differences (faster than 155.5): blue; near zero: green; positive differences (slower): red.
   Double_t red[NRGBs]   = {0.0, 0.0, 1.0};
   Double_t green[NRGBs] = {0.0, 1.0, 0.0};
   Double_t blue[NRGBs]  = {1.0, 0.0, 0.0};
   Int_t nColors = 100;
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, nColors);
   gStyle->SetNumberContours(nColors);
   
   hAvgDiff->SetContour(99);
   hAvgDiff->Draw("COLZ");
   hAvgDiff->GetXaxis()->SetTitle("Column");
   hAvgDiff->GetYaxis()->SetTitle("Row");
   gPad->SetFixedAspectRatio();
   // Set the X and Y axis ranges to match the detector geometry.
   hAvgDiff->GetXaxis()->SetRangeUser(0, COLS);
   hAvgDiff->GetYaxis()->SetRangeUser(0, ROWS);
   hAvgDiff->SetStats(0);  // Hide the stats box.
   // Hide the palette if desired.
   TPaletteAxis* pal2 = (TPaletteAxis*)hAvgDiff->GetListOfFunctions()->FindObject("palette");
   if(pal2) {
      pal2->SetX1NDC(1.0);
   }
   
   c2_deviation->SaveAs(Form("%s/BlockDeviationMap_run%d.png", outDir.Data(), nrun));
   
   // Cleanup.
   f->Close();
   delete f;
   
   return 0;
}
