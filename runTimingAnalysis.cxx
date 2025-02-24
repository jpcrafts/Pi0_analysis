// File: runTimingAnalysis.cpp
// Description:
//   1) Reads a ROOT file from a specified directory,
//   2) Applies HMS cuts (including the H.gtr.th cut) and fills several 1D timing histograms using TTreeReader,
//      (2D maps are disabled in this version).
//   3) The code applies block-level timing offsets (from a CSV file) to each block's ADC-TDC diff time.
//   4) It computes new energy–weighted cluster times (newClusterTimes) from the corrected block times.
//   5) Histograms are filled as follows:
//         - hClusterTime: Filled with raw cluster times (clusT) from the file.
//         - hNewClusterTime: Filled with new energy–weighted cluster times computed from corrected block times.
//         - hGoodBlockTime, hIntime, and hIntimeClusters are filled per block using corrected block times,
//           with selection criteria based on the newClusterTimes.
//         - For hIntimeClusters, only blocks from the two best clusters (those with newClusterTimes closest to 150 ns) are used.
//   6) A canvas is arranged in a 3×2 layout as follows:
//         Top row: Pad 1 = hClusterTime, Pad 2 = hNewClusterTime, Pad 3 = hGoodBlockTime.
//         Bottom row: Pad 4 = hIntime, Pad 5 = hIntimeClusters, Pad 6 = HMS Cuts Summary.
//   7) Every 10k events, progress is printed to cout.
//   8) The HMS cut values are defined as constants so that changing them updates both the selection and the summary.
//   
// Compile with:
//   g++ -O2 -o runTimingAnalysis runTimingAnalysis.cpp `root-config --cflags --libs`
// Run with:
//   ./runTimingAnalysis <nrun> <offset_csv_file>
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
#include "TPaveText.h" // For annotation and summary text
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TPaletteAxis.h"
#include "TColor.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <limits>
#include <algorithm>
#include <map>
#include <fstream>
#include <sstream>

using namespace std;

// HMS cut constants.
const double edtmtdcCut = 0.1;
const double hdeltaLow = -8.5;
const double hdeltaHigh = 8.5;
const double hcaltotCut = 0.6;
const double hcernpeCut = 1.0;
const double gtrthCut = 0.09;

static const int ROWS = 36;
static const int COLS = 30;
static const int N_BLOCKS = ROWS * COLS;

// Helper function to load block timing offsets from a CSV file.
std::map<int, double> loadOffsets(const char* filename) {
   std::map<int, double> offsets;
   ifstream file(filename);
   if (!file.is_open()){
      cout << "Error opening CSV file " << filename << endl;
      return offsets;
   }
   string line;
   // Skip header line.
   getline(file, line);
   while(getline(file, line)) {
      if(line.empty()) continue;
      istringstream ss(line);
      string blockStr, offsetStr, qualityStr;
      if(getline(ss, blockStr, ',') && getline(ss, offsetStr, ',') && getline(ss, qualityStr, ',')) {
         try {
            int blockID = stoi(blockStr);
            double offset = stod(offsetStr);
            offsets[blockID] = offset;
         } catch(const exception &e) {
            cout << "Error parsing line: " << line << " (" << e.what() << ")" << endl;
         }
      }
   }
   file.close();
   return offsets;
}

int main(int argc, char *argv[]){
   if(argc < 3){
      cout << "Usage: " << argv[0] << " <nrun> <offset_csv_file>" << endl;
      return 1;
   }

   // Load block-level offsets.
   map<int, double> blockOffset = loadOffsets(argv[2]);

   gStyle->SetOptFit(1111);
   int nrun = atoi(argv[1]);
   TString inFileName = Form("/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%d_1_-1.root", nrun);
   TFile *f = TFile::Open(inFileName, "READ");
   if(!f || f->IsZombie()){
      cout << "Error opening file " << inFileName << endl;
      return 1;
   }

   TTree *tree = (TTree *)f->Get("T");
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
   tree->SetBranchStatus("H.gtr.th", 1);

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
   TTreeReaderValue<Double_t> gtrth(reader, "H.gtr.th");

   // Define histograms.
   TH1F *hClusterTime = new TH1F("hClusterTime", "Raw Cluster Time Histogram (HMS Cuts);Cluster Time (ns);Counts", 200, 100, 200);
   TH1F *hNewClusterTime = new TH1F("hNewClusterTime", "New Energy Weighted Cluster Time;Cluster Time (ns);Counts", 200, 100, 200);
   TH1F *hGoodBlockTime = new TH1F("hGoodBlockTime", "Good Cluster Block Time Histogram;ADC-TDC Diff Time (ns);Counts", 200, 120, 180);
   TH1F *hIntime = new TH1F("hIntime", "In-Time ADC-TDC Diff Time;ADC-TDC Diff Time (ns);Counts", 200, 120, 180);
   hIntime->SetLineColor(kBlue);
   TH1F *hIntimeClusters = new TH1F("hIntimeClusters", "In-Time Cluster Time Histogram (selected clusters);ADC-TDC Diff Time (ns);Counts", 200, 140, 160);

   // --- Loop over events ---
   int processedEvents = 0;
   while(reader.Next()){
      // HMS cuts.
      if(*edtmtdc >= edtmtdcCut || *hdelta <= hdeltaLow || *hdelta >= hdeltaHigh ||
         *hcaltot <= hcaltotCut || *hcernpe <= hcernpeCut)
         continue;
      if(fabs(*gtrth) > gtrthCut)
         continue;

      processedEvents++;
      int numClus = static_cast<int>(*nclust);

      // Fill hClusterTime using raw cluster times.
      for (int i = 0; i < numClus; i++){
         double tVal = clusT[i];
         if(tVal > 100 && tVal < 200)
            hClusterTime->Fill(tVal);
      }

      // --- Compute corrected block times once per block ---
      int nBlocks_evt = block_t.GetSize();
      vector<double> correctedBlockTimes(nBlocks_evt, -1);
      for (int ib = 0; ib < nBlocks_evt; ib++){
         double offset = (blockOffset.find(ib) != blockOffset.end()) ? blockOffset[ib] : 0.0;
         correctedBlockTimes[ib] = block_t[ib] + offset;
      }

      // --- Compute new energy-weighted cluster times using corrected block times ---
      vector<double> newClusterTimes(numClus, -1);
      for (int i = 0; i < numClus; i++){
         double sumWeightedTime = 0.0;
         double sumEnergy = 0.0;
         for (int ib = 0; ib < nBlocks_evt; ib++){
            int cid = static_cast<int>(block_clusterID[ib]);
            if(cid == i){
               double energy = block_e[ib];
               if(energy > 0){
                  sumWeightedTime += correctedBlockTimes[ib] * energy;
                  sumEnergy += energy;
               }
            }
         }
         if(sumEnergy > 0)
            newClusterTimes[i] = sumWeightedTime / sumEnergy;
      }
      // Fill hNewClusterTime histogram.
      for (int i = 0; i < numClus; i++){
         if(newClusterTimes[i] > 0)
            hNewClusterTime->Fill(newClusterTimes[i]);
      }

      // --- Use newClusterTimes for selection for the last three histograms ---
      vector<bool> goodCluster(numClus, false);
      vector<double> narrowTimes;
      narrowTimes.reserve(numClus);
      for (int i = 0; i < numClus; i++){
         if(clusE[i] > 0.6 && newClusterTimes[i] > 130 && newClusterTimes[i] < 170)
            goodCluster[i] = true;
         if(clusE[i] > 0.6 && newClusterTimes[i] >= 147 && newClusterTimes[i] <= 153)
            narrowTimes.push_back(newClusterTimes[i]);
      }
      bool intimeEvent = false;
      if(narrowTimes.size() >= 2){
         sort(narrowTimes.begin(), narrowTimes.end());
         for(size_t i = 0; i < narrowTimes.size()-1; i++){
            if((narrowTimes[i+1] - narrowTimes[i]) <= 2){
               intimeEvent = true;
               break;
            }
         }
      }
      
      // --- Determine the two best clusters based on newClusterTimes (closest to 150 ns) ---
      int bestIndex = -1, secondBestIndex = -1;
      double bestDiff = numeric_limits<double>::max();
      double secondBestDiff = numeric_limits<double>::max();
      for (int i = 0; i < numClus; i++){
         if(clusE[i] > 0.6 && newClusterTimes[i] > 0){
            double diff = fabs(newClusterTimes[i] - 150.0);
            if(diff < bestDiff){
               secondBestDiff = bestDiff;
               secondBestIndex = bestIndex;
              	bestDiff = diff;
              	bestIndex = i;
            } else if(diff < secondBestDiff){
               secondBestDiff = diff;
              	secondBestIndex = i;
            }
         }
      }
      
      // --- Fill block-level histograms using corrected block times ---
      for (int ib = 0; ib < nBlocks_evt; ib++){
         int cid = static_cast<int>(block_clusterID[ib]);
         if(cid >= 0 && cid < numClus && goodCluster[cid]){
            if(block_e[ib] > 0){
               double newBlockTime = correctedBlockTimes[ib];
               if(newBlockTime > 100 && newBlockTime < 200){
                  hGoodBlockTime->Fill(newBlockTime);
                  if(intimeEvent){
                     hIntime->Fill(newBlockTime);
                     if(cid == bestIndex || cid == secondBestIndex)
                        hIntimeClusters->Fill(newBlockTime);
                  }
               }
            }
         }
      }
      
      // Print progress every 10k events.
      if(processedEvents % 10000 == 0)
         cout << "Processed " << processedEvents << " events..." << endl;
   }
   cout << "Processed " << processedEvents << " events for run " << nrun << endl;

   // ------------------------------
   // Fitting routines
   // ------------------------------
   double fitRangeMinVal = 145.0, fitRangeMaxVal = 155.0;
   TF1 *bgLeftFit = new TF1("bgLeft", "pol2", 117, fitRangeMinVal);
   hClusterTime->Fit(bgLeftFit, "R+");
   TF1 *bgRightFit = new TF1("bgRight", "pol2", fitRangeMaxVal, 184);
   hClusterTime->Fit(bgRightFit, "R+");

   TH1F *hClusterTimeNoBg = (TH1F *)hClusterTime->Clone("hClusterTimeNoBg");
   for (int bin = 1; bin <= hClusterTime->GetNbinsX(); bin++){
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

   double fitRangeBlockMinVal = 147.0, fitRangeBlockMaxVal = 153.0;
   TF1 *gausPlusBGBlockFit = new TF1("gausPlusBGBlock", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",
                                     fitRangeBlockMinVal, fitRangeBlockMaxVal);
   gausPlusBGBlockFit->SetParameters(60000, 150, 1.0, 100);
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

   TF1 *gausPlusBGIntime = new TF1("gausPlusBGIntime", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", 145, 155);
   gausPlusBGIntime->SetParameters(12000, 150, 1.0, 2000);
   hIntime->Fit(gausPlusBGIntime, "R");
   double meanIntime = gausPlusBGIntime->GetParameter(1);
   double sigmaIntime = gausPlusBGIntime->GetParameter(2);
   cout << "In-Time Histogram Fit: Mean = " << meanIntime << " ns, Sigma = " << sigmaIntime << " ns" << endl;

   TF1 *doubleGausIntimeClust = new TF1("doubleGausIntimeClust",
         "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)", 148, 152);
   doubleGausIntimeClust->SetParameters(1500, 150, 0.5, 500, 150, 1.0);
   hIntimeClusters->Fit(doubleGausIntimeClust, "R");
   double meanIntimeClust1 = doubleGausIntimeClust->GetParameter(1);
   double sigmaIntimeClust1 = doubleGausIntimeClust->GetParameter(2);
   double meanIntimeClust2 = doubleGausIntimeClust->GetParameter(4);
   double sigmaIntimeClust2 = doubleGausIntimeClust->GetParameter(5);
   cout << "In-Time Cluster Time Histogram Double Gaussian Fit:" << endl
        << "  First Gaussian:  Mean = " << meanIntimeClust1 << " ns, Sigma = " << sigmaIntimeClust1 << " ns" << endl
        << "  Second Gaussian: Mean = " << meanIntimeClust2 << " ns, Sigma = " << sigmaIntimeClust2 << " ns" << endl;

   // --- Create a canvas with 3 columns and 2 rows in rearranged order ---
   TCanvas *c1_all = new TCanvas("c1_all", "All 1D Histograms", 2100, 1200);
   c1_all->Divide(3,2);

   // Top row:
   // Pad 1: hClusterTime (Raw Cluster Time)
   c1_all->cd(1);
   hClusterTime->SetStats(0);
   hClusterTime->SetLineColor(kBlue);
   hClusterTime->Draw("hist");
   hClusterTimeNoBg->SetLineColorAlpha(kRed, 0.0);
   hClusterTimeNoBg->Draw("hist same");
   gausFitRaw->SetLineColor(kMagenta);
   gausFitRaw->Draw("same");

   // Pad 2: hNewClusterTime (New Energy Weighted Cluster Time)
   c1_all->cd(2);
   hNewClusterTime->SetStats(0);
   hNewClusterTime->SetLineColor(kBlue);
   hNewClusterTime->Draw("hist");

   // Pad 3: hGoodBlockTime (Good Cluster Block Time)
   c1_all->cd(3);
   hGoodBlockTime->SetStats(1);
   hGoodBlockTime->SetLineColor(kBlue);
   hGoodBlockTime->Draw("hist");
   gausPlusBGBlockFit->SetLineColor(kMagenta);
   gausPlusBGBlockFit->Draw("same");

   // Bottom row:
   // Pad 4: hIntime (Raw In-Time ADC-TDC Diff Time)
   c1_all->cd(4);
   hIntime->SetStats(1);
   hIntime->SetLineColor(kBlue);
   hIntime->Draw("hist");
   gausPlusBGIntime->SetLineColor(kMagenta);
   gausPlusBGIntime->Draw("same");

   // Pad 5: hIntimeClusters (In-Time Cluster Time for Selected Clusters)
   c1_all->cd(5);
   hIntimeClusters->SetStats(1);
   hIntimeClusters->SetLineColor(kBlue);
   hIntimeClusters->Draw("hist");
   doubleGausIntimeClust->SetLineColor(kMagenta);
   doubleGausIntimeClust->Draw("same");

   // Pad 6: HMS Cuts Summary (and progress info)
   c1_all->cd(6);
   TPaveText *ptSummary = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
   ptSummary->AddText(Form("Run %d HMS Cuts Summary:", nrun));
   ptSummary->AddText(Form("T.hms.hEDTM_tdcTimeRaw < %.2f", edtmtdcCut));
   ptSummary->AddText(Form("H.gtr.dp between %.0f and %.0f", hdeltaLow, hdeltaHigh));
   ptSummary->AddText(Form("H.cal.etotnorm > %.2f", hcaltotCut));
   ptSummary->AddText(Form("H.cer.npeSum > %.1f", hcernpeCut));
   ptSummary->AddText(Form("H.gtr.th within ±%.2f", gtrthCut));
   ptSummary->AddText("");
   ptSummary->AddText(Form("Processed %d events", processedEvents));
   ptSummary->SetFillColor(0);
   ptSummary->SetBorderSize(0);
   ptSummary->Draw();

   c1_all->SaveAs(Form("/volatile/hallc/nps/jpcrafts/Plots/All1DHistograms_run%d.png", nrun));

   f->Close();
   delete f;
   return 0;
}
