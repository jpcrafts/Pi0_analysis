/*  scan_scaler_timeSeries.cxx
 *  ---------------------------------------------------------------
 *  Plots raw scaler counts & instantaneous rates vs time.
 *  (No charge‑norm, no CPU‑live‑time, no EDTM subtraction.)
 *
 *  Compile:
 *    g++ -std=c++17 -O2 -Wall scan_scaler_timeSeries.cxx \
 *        $(root-config --cflags --libs) -o scan_scaler_timeSeries
 *
 *  Run (uses default directory & filename pattern below):
 *    ./scan_scaler_timeSeries 2093
 *
 *  Run (override with your own path *and* pattern):
 *    ./scan_scaler_timeSeries 2093 "/my/dir/nps_hms_skim_%d_1_-1.root"
 */

 #include <iostream>
 #include <vector>
 #include <string>
 #include <cstdio>
 
 #include "TFile.h"
 #include "TTree.h"
 #include "TString.h"
 #include "TCanvas.h"
 #include "TGraph.h"
 #include "TMultiGraph.h"
 #include "TLegend.h"
 #include "TStyle.h"
 #include "TAxis.h"          // full definition of TAxis
 
 /*------------------------------------------------------------------*/
 struct ScalVec {
   std::string  title;
   int          color;
   std::vector<double> cnt;
   std::vector<double> rate;
 };
 
 /*------------------------------------------------------------------*/
 int main(int argc, char* argv[])
 {
   if(argc < 2){
     std::cerr << "\nUsage: " << argv[0] << " RUNNUMBER "
               << "[full_pattern_with_%d]\n\n";
     return 1;
   }
 
   const int run = std::stoi(argv[1]);
 
   /* ----------------------------------------------------------------
    *  Default *full* pattern: directory + filename with “%d” for run.
    *  You can override the whole thing as the 2nd command‑line arg.
    * ----------------------------------------------------------------*/
   const std::string defaultPattern =
       "/lustre24/expphy/cache/hallc/c-nps/analysis/pass1/replays/skim/"
       "nps_hms_skim_%d_1_-1.root";
 
   const std::string pattern =
       (argc > 2) ? argv[2] : defaultPattern;
 
   /* --- build the file name -------------------------------------- */
   char buf[512];
   std::snprintf(buf, sizeof(buf), pattern.c_str(), run);
   TString fname(buf);
 
   TFile *f = TFile::Open(fname, "READ");
   if(!f || f->IsZombie()){
     std::cerr << "Cannot open " << fname << "\n";
     return 1;
   }
 
   TTree *tscal = static_cast<TTree*>(f->Get("TSH"));
   if(!tscal){
     std::cerr << "TSH tree not found in " << fname << "!\n";
     return 1;
   }
 
   /* ---------- branch handles ---------- */
   double tMHz = 0;
   double cS1X = 0, cT1 = 0, cT4 = 0, cT3 = 0, cEDTM = 0;
 
   tscal->SetBranchAddress("H.1MHz.scalerTime", &tMHz);
   tscal->SetBranchAddress("H.S1X.scaler",      &cS1X);
   tscal->SetBranchAddress("H.hTRIG1.scaler",   &cT1);
   tscal->SetBranchAddress("H.hTRIG4.scaler",   &cT4);
   tscal->SetBranchAddress("H.hTRIG3.scaler",   &cT3);
   tscal->SetBranchAddress("H.EDTM.scaler",     &cEDTM);
 
   std::vector<double> time;
   ScalVec scal[5] = {
     {"S1X",   kBlue},
     {"TRIG1", kGreen+2},
     {"TRIG4", kRed},
     {"TRIG3", kMagenta},
     {"EDTM",  kOrange+7}
   };
 
   /* ---------- read the tree ----------- */
   const Long64_t nEnt = tscal->GetEntries();
   time.reserve(nEnt);
   for(Long64_t i = 0; i < nEnt; ++i){
     tscal->GetEntry(i);
     time.push_back(tMHz);
 
     scal[0].cnt.push_back(cS1X);
     scal[1].cnt.push_back(cT1);
     scal[2].cnt.push_back(cT4);
     scal[3].cnt.push_back(cT3);
     scal[4].cnt.push_back(cEDTM);
   }
 
   /* ---- instantaneous rates ----------- */
   for(auto &sv : scal){
     sv.rate.resize(sv.cnt.size(), 0.);
     for(size_t i = 1; i < sv.cnt.size(); ++i){
       double dN = sv.cnt[i] - sv.cnt[i-1];
       double dt = time[i] - time[i-1];
       sv.rate[i] = (dt > 0) ? dN / dt : 0.;
     }
   }
 
   gStyle->SetOptStat(0);
 
   /* ---- counts vs time canvas --------- */
   TCanvas cCnt("cCnt", "Counts vs time", 1400, 900);
   TMultiGraph mgCnt;
   TLegend legCnt(0.75, 0.7, 0.9, 0.9);
 
   for(const auto &sv : scal){
     TGraph *gr = new TGraph(time.size(), time.data(), sv.cnt.data());
     gr->SetLineColor(sv.color); gr->SetLineWidth(2);
     mgCnt.Add(gr, "L");  legCnt.AddEntry(gr, sv.title.c_str(), "l");
   }
   mgCnt.Draw("A");
   mgCnt.GetXaxis()->SetTitle("Time (s)");
   mgCnt.GetYaxis()->SetTitle("Scaler counts (cumulative)");
   legCnt.Draw();
   cCnt.SaveAs(Form("counts_vs_time_%d.pdf", run));
 
   /* ---- rates vs time canvas ---------- */
   TCanvas cRate("cRate", "Rates vs time", 1400, 900);
   TMultiGraph mgRate;
   TLegend legRate(0.75, 0.7, 0.9, 0.9);
 
   for(const auto &sv : scal){
     TGraph *gr = new TGraph(time.size(), time.data(), sv.rate.data());
     gr->SetLineColor(sv.color); gr->SetLineWidth(2);
     mgRate.Add(gr, "L");  legRate.AddEntry(gr, sv.title.c_str(), "l");
   }
   mgRate.Draw("A");
   mgRate.GetXaxis()->SetTitle("Time (s)");
   mgRate.GetYaxis()->SetTitle("Instantaneous rate (Hz)");
   legRate.Draw();
   cRate.SaveAs(Form("rates_vs_time_%d.pdf", run));
 
   /* ---- console summary --------------- */
   double runTime = (time.empty() ? 0 : time.back() - time.front());
   std::cout << "\nRun " << run << "  (Δt = " << runTime << " s)\n";
   for(const auto &sv : scal){
     double avgRate = (runTime > 0) ? (sv.cnt.back() - sv.cnt.front()) / runTime : 0.;
     std::cout << "  " << sv.title << "  avg rate = " << avgRate << " Hz\n";
   }
   std::cout << std::endl;
 
   return 0;
 }
 