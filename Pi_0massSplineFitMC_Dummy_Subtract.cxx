/*************************************************************
  Pi_0massSplineFitMC_Dummy_Subtract.cxx

  This code:
    1) Reads a real data TTree from a ROOT file, applies HPC/HMS cuts 
       to build a background histogram (in the time domain) which is spline‐fitted 
       and shifted.
    2) Builds a raw signal histogram from the time information (using the first cluster time)
       and then produces a nominal subtracted histogram (hSub = raw signal – shifted BG).
    3) Separately reads a dummy (target–wall) file; splits its events into upstream vs. downstream 
       (based on H.gtr.y), scales each piece by its effective charge and wall–thickness factor,
       and subtracts the combined dummy histogram from hSub.
    4) Both the real data and dummy histograms are then scaled by their effective BCM4A charges 
       (in µC) provided on the command line so that the final histograms are in consistent units.
    5) Next, the code gathers cluster pairs (using the full array of cluster times) for a Toy MC loop 
       and produces a 5–pad canvas.
       
  IMPORTANT: The branch "NPS.cal.clusT" is set only once (as an array) and for each event we use clusT[0] 
  (if at least one cluster exists) as the representative time for the time–domain histograms.

  Usage:
    ./Pi_0massSplineFitMC_Dummy_Subtract realData.root output.pdf dummyData.root effectiveCharge_data dummyCharge_uC

  Compile example:
    g++ -O2 -std=c++17 -o Pi_0massSplineFitMC_Dummy_Subtract \
      Pi_0massSplineFitMC_Dummy_Subtract.cxx alglib_src/*.cpp \
      `root-config --cflags --libs` -I. -Ialglib_src -I`root-config --incdir`/TMVA \
      -lTMVA -lRooFit -lRooFitCore -lgsl -lgslcblas
*************************************************************/

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TSystem.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>  // for std::stod

// ============================================================================
// 1) Normalization & smoothing functions
// ============================================================================
std::vector<double> normalizePeaks(const std::vector<double>& y, int peakFilterWindow, double minPeakFraction)
{
    int n = y.size();
    std::vector<double> yNorm = y;
    std::vector<int> candidatePeaks;
    double maxVal = *std::max_element(y.begin(), y.end());
    double minPeakHeight = minPeakFraction * maxVal;

    for (int i = 1; i < n - 1; i++){
        if (y[i] > y[i-1] && y[i] > y[i+1] && y[i] > minPeakHeight){
            candidatePeaks.push_back(i);
        }
    }
    std::vector<int> filteredPeaks;
    for (int p : candidatePeaks){
        if (filteredPeaks.empty()){
            filteredPeaks.push_back(p);
        } else {
            int last = filteredPeaks.back();
            if (p - last < peakFilterWindow){
                if (y[p] > y[last])
                    filteredPeaks.back() = p;
            } else {
                filteredPeaks.push_back(p);
            }
        }
    }
    if (filteredPeaks.empty()) return yNorm;

    double sumPeaks = 0.0;
    for (int idx : filteredPeaks)
        sumPeaks += y[idx];
    double avgPeak = sumPeaks / filteredPeaks.size();

    for (int p : filteredPeaks){
        double scale = avgPeak / y[p];
        int halfWin = peakFilterWindow / 2;
        for (int i = std::max(0, p - halfWin); i <= std::min(n - 1, p + halfWin); i++){
            double weight = 1.0 - (std::fabs(i - p) / double(peakFilterWindow + 1));
            yNorm[i] = y[i]*(1 - weight) + (y[i]*scale)*weight;
        }
    }
    return yNorm;
}

std::vector<double> normalizeTroughs(const std::vector<double>& y, int troughFilterWindow, double troughToleranceFactor)
{
    int n = y.size();
    std::vector<double> yNorm = y;
    std::vector<int> candidateTroughs;
    double globalMin = *std::min_element(y.begin(), y.end());
    double globalMax = *std::max_element(y.begin(), y.end());
    double troughTolerance = troughToleranceFactor * (globalMax - globalMin);

    for (int i = 1; i < n - 1; i++){
        if (y[i] < y[i-1] && y[i] < y[i+1] && y[i] <= globalMin + troughTolerance)
            candidateTroughs.push_back(i);
    }
    std::vector<int> filteredTroughs;
    for (int t : candidateTroughs){
        if (filteredTroughs.empty()){
            filteredTroughs.push_back(t);
        } else {
            int last = filteredTroughs.back();
            if (t - last < troughFilterWindow){
                if (y[t] < y[last])
                    filteredTroughs.back() = t;
            } else {
                filteredTroughs.push_back(t);
            }
        }
    }
    if (filteredTroughs.empty()) return yNorm;

    double sumTroughs = 0.0;
    for (int idx : filteredTroughs)
        sumTroughs += y[idx];
    double avgTrough = sumTroughs / filteredTroughs.size();

    for (int t : filteredTroughs){
        double scale = avgTrough / y[t];
        int halfWin = troughFilterWindow / 2;
        for (int i = std::max(0, t - halfWin); i <= std::min(n - 1, t + halfWin); i++){
            double weight = 1.0 - (std::fabs(i - t) / double(troughFilterWindow + 1));
            yNorm[i] = y[i]*(1 - weight) + (y[i]*scale)*weight;
        }
    }
    return yNorm;
}

std::vector<double> regionSmooth(const std::vector<double>& data, int smoothWindowPeak, int smoothWindowTrough)
{
    int n = data.size();
    std::vector<double> smoothed(n);
    double globalMax = *std::max_element(data.begin(), data.end());
    double globalMin = *std::min_element(data.begin(), data.end());
    double midValue  = (globalMax + globalMin) / 2.0;
    for (int i = 0; i < n; i++){
        int window = (data[i] > midValue) ? smoothWindowPeak : smoothWindowTrough;
        double sum = 0.0;
        int count = 0;
        for (int j = i - window/2; j <= i + window/2; j++){
            if (j >= 0 && j < n){
                sum += data[j];
                count++;
            }
        }
        smoothed[i] = (count > 0 ? sum / count : data[i]);
    }
    return smoothed;
}

// ============================================================================
// 2) Cluster-level cut: Exclude clusters near the outer one-block layer
// ============================================================================
static bool isGoodCluster(double e, double t, double x, double y)
{
    bool passBasic = (e >= 0.6 && t >= 149.0 && t <= 151.0);
    double xInnerMin = -29.16, xInnerMax = 29.16;
    double yInnerMin = -35.64, yInnerMax = 35.64;
    bool insideEdges = (x > xInnerMin && x < xInnerMax && y > yInnerMin && y < yInnerMax);
    return (passBasic && insideEdges);
}

// ============================================================================
// 3) Main: Build 5-pad output with Toy MC and a shaded error region (Pi0 Mass)
// ============================================================================
struct ClusterRec {
    int timeIndex;   // time bin index (0 to nBins-1)
    double massUnw;  // unweighted invariant mass of the cluster pair
};

int main(int argc, char* argv[])
{
    if (argc < 6){
        std::cerr << "Usage: " << argv[0]
                  << " <realData.root> <output.pdf> <dummyData.root> <realCharge_uC> <dummyCharge_uC>\n";
        return 1;
    }
    std::string inputFileName  = argv[1];
    std::string outputFileName = argv[2];
    std::string dummyFileName  = argv[3];
    double effective_charge_data = std::stod(argv[4]);  // in µC
    double dummyEffectiveCharge  = std::stod(argv[5]);   // in µC

    // ----------------- Set up HPC/HMS cuts -----------------
    double edtmtdcCut   = 0.1;
    double hdeltaLowCut = -8.5, hdeltaHighCut = 8.5;
    double hcaltotCut   = 0.6, hcernpeCut = 1.0;
    double gtrthCut     = 0.09, gtrphCut = 0.09;
    auto passHMSCuts = [&](double edtmtdc, double hdelta, double hcaltot,
                             double hcernpe, double gtrth, double gtrph) -> bool {
        if (edtmtdc >= edtmtdcCut) return false;
        if (hdelta < hdeltaLowCut || hdelta > hdeltaHighCut) return false;
        if (hcaltot <= hcaltotCut) return false;
        if (hcernpe <= hcernpeCut) return false;
        if (std::fabs(gtrth) > gtrthCut) return false;
        if (std::fabs(gtrph) > gtrphCut) return false;
        return true;
    };

    // ----------------- Define regions -----------------
    double bgLow = 113.0, bgHigh = 142.5;
    int nBins = 650;  // time bins for BG and signal regions
    double shiftVal = 28.05;
    double sigLow = 141.789, sigHigh = 171.289;

    // ----------------- Open real data file and get TTree -----------------
    TFile* inFile = TFile::Open(inputFileName.c_str(), "READ");
    if (!inFile || inFile->IsZombie()){
        std::cerr << "Error: cannot open " << inputFileName << "\n";
        return 1;
    }
    TTree* tree = dynamic_cast<TTree*>( inFile->Get("T") );
    if (!tree){
        std::cerr << "Error: TTree 'T' not found.\n";
        inFile->Close();
        return 1;
    }

    // Instead of setting "NPS.cal.clusT" as a separate double, we set up the branches
    // for cluster pairing here (as an array), and then use clusT[0] for time histograms.
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    tree->SetBranchStatus("H.gtr.dp",              1);
    tree->SetBranchStatus("H.cal.etotnorm",        1);
    tree->SetBranchStatus("H.cer.npeSum",          1);
    tree->SetBranchStatus("H.gtr.th",              1);
    tree->SetBranchStatus("H.gtr.ph",              1);
    // For cluster pairing, enable the following branches:
    tree->SetBranchStatus("NPS.cal.nclust",        1);
    tree->SetBranchStatus("NPS.cal.clusE",         1);
    tree->SetBranchStatus("NPS.cal.clusT",         1);
    tree->SetBranchStatus("NPS.cal.clusX",         1);
    tree->SetBranchStatus("NPS.cal.clusY",         1);

    double edtmtdc=0, hdelta=0, hcaltot=0, hcernpe=0, gtrth=0, gtrph=0;
    tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
    tree->SetBranchAddress("H.gtr.dp",               &hdelta);
    tree->SetBranchAddress("H.cal.etotnorm",         &hcaltot);
    tree->SetBranchAddress("H.cer.npeSum",           &hcernpe);
    tree->SetBranchAddress("H.gtr.th",               &gtrth);
    tree->SetBranchAddress("H.gtr.ph",               &gtrph);

    double nclustDouble = 0.0;
    static const int MAXC = 10000;
    double clusE[MAXC], clusT[MAXC], clusX[MAXC], clusY[MAXC];
    tree->SetBranchAddress("NPS.cal.nclust", &nclustDouble);
    tree->SetBranchAddress("NPS.cal.clusE", clusE);
    tree->SetBranchAddress("NPS.cal.clusT", clusT);
    tree->SetBranchAddress("NPS.cal.clusX", clusX);
    tree->SetBranchAddress("NPS.cal.clusY", clusY);

    // ----------------- (A) Build BG and Signal histograms (real data) -----------------
    // Use the first cluster time (clusT[0]) if at least one cluster exists.
    TH1F* hBG = new TH1F("hBG", "BG Region", nBins, bgLow, bgHigh);
    hBG->SetDirectory(nullptr);
    TH1F* hSignalRaw = new TH1F("hSignalRaw", "Signal Region Raw", nBins, sigLow, sigHigh);
    hSignalRaw->SetDirectory(nullptr);

    Long64_t nAll = tree->GetEntries();
    for (Long64_t i = 0; i < nAll; i++){
        tree->GetEntry(i);
        if (!passHMSCuts(edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph))
            continue;
        if (nclustDouble < 1) continue;
        double eventTime = clusT[0];  // use first cluster time
        if (eventTime >= bgLow && eventTime <= bgHigh)
            hBG->Fill(eventTime);
        if (eventTime >= sigLow && eventTime <= sigHigh)
            hSignalRaw->Fill(eventTime);
    }

    // ----------------- (A.1) Scale real data histograms by effective charge -----------------
    hBG->Scale(1.0 / effective_charge_data);
    hSignalRaw->Scale(1.0 / effective_charge_data);

    // ----------------- (B) Spline Fit & Smoothing for BG (real data) -----------------
    std::vector<double> xBG(nBins), yBG_vals(nBins);
    for (int b = 1; b <= nBins; b++){
        xBG[b-1] = hBG->GetBinCenter(b);
        yBG_vals[b-1] = hBG->GetBinContent(b);
    }
    TGraph* grBG = new TGraph(nBins, &xBG[0], &yBG_vals[0]);
    grBG->Sort();
    TSpline3* splineBG = new TSpline3("splineBG", grBG);
    std::vector<double> ySpline(nBins);
    for (int i = 0; i < nBins; i++){
        ySpline[i] = splineBG->Eval(xBG[i]);
    }
    int peakWindow = 10;       
    double minPeakFrac = 0.75;
    int troughWindow = 15;     
    double troughFrac = 0.05;
    int swPeak = 2, swTrough = 5;
    auto yPeaks   = normalizePeaks(ySpline, peakWindow, minPeakFrac);
    auto yTroughs = normalizeTroughs(yPeaks, troughWindow, troughFrac);
    auto yFinal   = regionSmooth(yTroughs, swPeak, swTrough);

    // ----------------- (C) Create nominal subtracted histogram for real data -----------------
    std::vector<double> xShift(nBins), yShift(nBins);
    for (int i = 0; i < nBins; i++){
        xShift[i] = xBG[i] + shiftVal;
        yShift[i] = yFinal[i];
    }
    TGraph* grBGshift = new TGraph(nBins, &xShift[0], &yShift[0]);
    TH1F* hSub = (TH1F*) hSignalRaw->Clone("hSub");
    hSub->SetDirectory(nullptr);
    hSub->SetTitle("Subtracted Signal");
    for (int b = 1; b <= nBins; b++){
        double center = hSignalRaw->GetBinCenter(b);
        double dVal   = hSignalRaw->GetBinContent(b);
        double bgv    = grBGshift->Eval(center);
        hSub->SetBinContent(b, dVal - bgv);
    }

    // ----------------- (D) Dummy Subtraction with charge normalization -----------------
    double scale_factor_dummy_upstream   = 8.467;
    double scale_factor_dummy_downstream = 4.256;

    TFile* dummyFile = TFile::Open(dummyFileName.c_str(), "READ");
    if (!dummyFile || dummyFile->IsZombie()){
        std::cerr << "Warning: cannot open dummy file " << dummyFileName << ". Skipping dummy subtraction.\n";
    } 
    else {
        TTree* dummyTree = dynamic_cast<TTree*>( dummyFile->Get("T") );
        if (!dummyTree){
            std::cerr << "Warning: TTree 'T' not found in dummy file. Skipping dummy subtraction.\n";
        } else {
            dummyTree->SetBranchStatus("*", 0);
            dummyTree->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
            dummyTree->SetBranchStatus("H.gtr.dp",              1);
            dummyTree->SetBranchStatus("H.cal.etotnorm",        1);
            dummyTree->SetBranchStatus("H.cer.npeSum",          1);
            dummyTree->SetBranchStatus("H.gtr.th",              1);
            dummyTree->SetBranchStatus("H.gtr.ph",              1);
            dummyTree->SetBranchStatus("H.gtr.y",               1);
            dummyTree->SetBranchStatus("NPS.cal.clusT",         1);

            double edtmtdc_d=0, hdelta_d=0, hcaltot_d=0, hcernpe_d=0,
                   gtrth_d=0, gtrph_d=0, gtry_d=0;
            double clusT_dummy = 0; // for time histogram from dummy; assume one cluster per event
            dummyTree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc_d);
            dummyTree->SetBranchAddress("H.gtr.dp",               &hdelta_d);
            dummyTree->SetBranchAddress("H.cal.etotnorm",         &hcaltot_d);
            dummyTree->SetBranchAddress("H.cer.npeSum",           &hcernpe_d);
            dummyTree->SetBranchAddress("H.gtr.th",               &gtrth_d);
            dummyTree->SetBranchAddress("H.gtr.ph",               &gtrph_d);
            dummyTree->SetBranchAddress("H.gtr.y",                &gtry_d);
            dummyTree->SetBranchAddress("NPS.cal.clusT",          &clusT_dummy);

            // Create histograms for dummy events: upstream and downstream,
            // in the same signal window [sigLow, sigHigh].
            TH1F* hSignalRaw_dummy_up = new TH1F("hSignalRaw_dummy_up", "Dummy Target (upstream)", nBins, sigLow, sigHigh);
            hSignalRaw_dummy_up->SetDirectory(nullptr);
            TH1F* hSignalRaw_dummy_down = new TH1F("hSignalRaw_dummy_down", "Dummy Target (downstream)", nBins, sigLow, sigHigh);
            hSignalRaw_dummy_down->SetDirectory(nullptr);

            Long64_t nDummy = dummyTree->GetEntries();
            for (Long64_t iD = 0; iD < nDummy; iD++){
                dummyTree->GetEntry(iD);
                if (!passHMSCuts(edtmtdc_d, hdelta_d, hcaltot_d, hcernpe_d, gtrth_d, gtrph_d))
                    continue;
                if (clusT_dummy < sigLow || clusT_dummy > sigHigh)
                    continue;
                if (gtry_d > 0)
                    hSignalRaw_dummy_down->Fill(clusT_dummy);
                else
                    hSignalRaw_dummy_up->Fill(clusT_dummy);
            }

            // Scale the dummy histograms by 1/(dummyEffectiveCharge * scaleFactor)
            hSignalRaw_dummy_up->Scale(1.0 / (dummyEffectiveCharge * scale_factor_dummy_upstream));
            hSignalRaw_dummy_down->Scale(1.0 / (dummyEffectiveCharge * scale_factor_dummy_downstream));

            TH1F* hSignalRaw_dummy = (TH1F*) hSignalRaw_dummy_up->Clone("hSignalRaw_dummy");
            hSignalRaw_dummy->Add(hSignalRaw_dummy_down);

            std::cout << "[Dummy Subtraction] Subtracting dummy from hSub...\n";
            hSub->Add(hSignalRaw_dummy, -1.0);
        }
        dummyFile->Close();
    }

    // ----------------- (E) Build arrays for Toy MC (unchanged) -----------------
    std::vector<double> dataVal(nBins, 0.0), bgValArray(nBins, 0.0), bgErr(nBins, 0.0);
    for (int b = 1; b <= nBins; b++){
        dataVal[b-1] = hSignalRaw->GetBinContent(b);
        double c = hSignalRaw->GetBinCenter(b);
        double bEst = grBGshift->Eval(c);
        bgValArray[b-1] = bEst;
        bgErr[b-1] = std::sqrt(std::max(bEst, 0.0));
    }

    // ----------------- (F) Final pass: Gather cluster pairs (mass domain) -----------------
    // (Now using the already-set branch addresses for cluster info.)
    // Note: For cluster pairing, we loop through each event and use the array clusT for each event.
    std::vector<ClusterRec> clusterList;
    clusterList.reserve(2000000);
    double DNPS = 407.0;
    for (Long64_t iEvt = 0; iEvt < nAll; iEvt++){
        tree->GetEntry(iEvt);
        if (!passHMSCuts(edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph))
            continue;
        int nclust = static_cast<int>(nclustDouble);
        if (nclust < 2) continue;
        std::vector<int> keep;
        keep.reserve(nclust);
        for (int cID = 0; cID < nclust; cID++){
            double tC = clusT[cID];
            if (tC < sigLow || tC > sigHigh) continue;
            if (!isGoodCluster(clusE[cID], tC, clusX[cID], clusY[cID]))
                continue;
            keep.push_back(cID);
        }
        if (keep.size() < 2) continue;
        for (size_t ic = 0; ic < keep.size(); ic++){
            for (size_t jc = ic + 1; jc < keep.size(); jc++){
                int id1 = keep[ic];
                int id2 = keep[jc];
                double E1 = clusE[id1], E2 = clusE[id2];
                double x1 = clusX[id1], y1 = clusY[id1];
                double x2 = clusX[id2], y2 = clusY[id2];
                double dx = x1 - x2, dy = y1 - y2;
                double dist = std::sqrt(dx*dx + dy*dy);
                double theta = dist / DNPS;
                double s2 = std::sin(0.5 * theta);
                double Mgg = std::sqrt(4.0 * E1 * E2 * s2 * s2);
                double tMean = 0.5 * (clusT[id1] + clusT[id2]);
                int binID = hSignalRaw->FindBin(tMean);
                if (binID < 1 || binID > nBins) continue;
                int iTime = binID - 1;
                ClusterRec rec;
                rec.timeIndex = iTime;
                rec.massUnw = Mgg;
                clusterList.push_back(rec);
            }
        }
    }
    inFile->Close();
    std::cout << "Built clusterList, size = " << clusterList.size() << "\n";

    // ----------------- (G) Toy Monte Carlo Loop (unchanged) -----------------
    int Ntoys = 200;
    std::vector<TH1F*> toyMassHists;
    toyMassHists.reserve(Ntoys);
    int nMassBins = 200;
    double mLo = 0.0, mHi = 0.3;
    TRandom3 randGen(0);
    for (int itoy = 0; itoy < Ntoys; itoy++){
        std::vector<double> fraction(nBins, 0.0);
        for (int i = 0; i < nBins; i++){
            double dev = randGen.Gaus(0.0, bgErr[i]);
            double bToy = bgValArray[i] + dev;
            if (bToy < 0) bToy = 0;
            double subV = dataVal[i] - bToy;
            double f = 0.0;
            if (dataVal[i] > 1e-9)
                f = subV / dataVal[i];
            fraction[i] = f;
        }
        TH1F* hToy = new TH1F(Form("hToy_%d", itoy),
                              "Toy Pi0 Mass", nMassBins, mLo, mHi);
        for (auto & c : clusterList){
            double w = fraction[c.timeIndex];
            hToy->Fill(c.massUnw, w);
        }
        toyMassHists.push_back(hToy);
    }

    // ----------------- (H) Compute TMC Average (hMean) and RMS (hVar) -----------------
    TH1F* hMean = (TH1F*)toyMassHists[0]->Clone("hMean");
    hMean->Reset();
    TH1F* hVar = (TH1F*)toyMassHists[0]->Clone("hVar");
    hVar->Reset();
    for (int t = 0; t < Ntoys; t++){
        TH1F* hist = toyMassHists[t];
        for (int ib = 1; ib <= nMassBins; ib++){
            double val = hist->GetBinContent(ib);
            hMean->AddBinContent(ib, val);
        }
    }
    hMean->Scale(1.0 / double(Ntoys));
    for (int t = 0; t < Ntoys; t++){
        for (int ib = 1; ib <= nMassBins; ib++){
            double val = toyMassHists[t]->GetBinContent(ib);
            double diff = val - hMean->GetBinContent(ib);
            hVar->AddBinContent(ib, diff * diff);
        }
    }
    hVar->Scale(1.0 / double(Ntoys));
    for (int ib = 1; ib <= nMassBins; ib++){
        double var = hVar->GetBinContent(ib);
        hVar->SetBinContent(ib, std::sqrt(var));  // hVar now holds the RMS
    }
    TH1F* hMeanWithErr = (TH1F*) hMean->Clone("hMeanWithErr");
    for (int ib = 1; ib <= nMassBins; ib++){
        double err = hVar->GetBinContent(ib);
        hMeanWithErr->SetBinError(ib, err);
    }

    // ----------------- (I) Build a polygon for the error region -----------------
    std::vector<double> poly_x, poly_y;
    poly_x.reserve(2 * nMassBins);
    poly_y.reserve(2 * nMassBins);
    // Lower edge: hMean
    for (int i = 1; i <= nMassBins; i++){
        poly_x.push_back(hMean->GetBinCenter(i));
        poly_y.push_back(hMean->GetBinContent(i));
    }
    // Upper edge: hMean + hVar (reverse order)
    for (int i = nMassBins; i >= 1; i--){
        double scaleFactor = 1.0;
        double up = hMean->GetBinContent(i) + scaleFactor * hVar->GetBinContent(i);
        poly_x.push_back(hMean->GetBinCenter(i));
        poly_y.push_back(up);
    }
    TGraph* grErrorPoly = new TGraph(poly_x.size(), poly_x.data(), poly_y.data());
    grErrorPoly->SetFillColorAlpha(kRed, 1.0);
    grErrorPoly->SetFillStyle(3001);

    // ----------------- (J) Create a 5-pad canvas -----------------
    TCanvas* cMerged = new TCanvas("cMerged", "Combined Plots (5 Pads)", 1200, 3000);
    cMerged->Divide(1, 5);

    // Pad 1: BG region + spline fit
    cMerged->cd(1);
    // Draw original BG points clearly
    gBG.SetMarkerStyle(20);
    gBG.SetMarkerSize(0.8);
    gBG.Draw("AP");

    // Explicitly construct and plot the final smoothed BG curve (yF)
    TGraph* grBGfinal = new TGraph(nBins, &vx[0], &yF[0]);
    grBGfinal->SetLineColor(kBlue);
    grBGfinal->SetLineWidth(2);
    grBGfinal->Draw("L SAME");
    {
        TLatex lat1;
        lat1.SetNDC();
        lat1.SetTextSize(0.04);
        lat1.DrawLatex(0.15, 0.92, "BG Region + Spline Fit");
    }

    // Pad 2: Signal region overlay (raw, shifted BG, subtracted)
    cMerged->cd(2);
    double maxSig = hSignalRaw->GetMaximum();
    double maxBG_val = 0.0;
    for (int b = 1; b <= nBins; b++){
        double diff = hSignalRaw->GetBinContent(b) - hSub->GetBinContent(b);
        if (diff > maxBG_val)
            maxBG_val = diff;
    }
    double maxVal2 = std::max(maxSig, maxBG_val);
    TH2F* frame2 = new TH2F("frame2", "Signal Region;Time (ns);Counts", 10, sigLow, sigHigh, 10, 0, 1.2 * maxVal2);
    frame2->Draw("AXIS");
    TGraph* grSignalRawG = new TGraph(hSignalRaw);
    grSignalRawG->SetLineColor(kBlack);
    grSignalRawG->SetLineWidth(2);
    grSignalRawG->Draw("L same");
    TGraph* grShiftedG = new TGraph(nBins);
    for (int i = 0; i < nBins; i++){
        double xx = hSignalRaw->GetBinCenter(i+1);
        double val = grBGshift->Eval(xx);
        grShiftedG->SetPoint(i, xx, val);
    }
    grShiftedG->SetLineColor(kRed);
    grShiftedG->SetLineWidth(2);
    grShiftedG->Draw("L same");
    TGraph* grSubG = new TGraph(hSub);
    grSubG->SetLineColor(kGreen+2);
    grSubG->SetLineWidth(2);
    grSubG->Draw("L same");
    {
        TLatex lat2;
        lat2.SetNDC();
        lat2.SetTextSize(0.04);
        lat2.DrawLatex(0.15, 0.92, "Signal vs Shifted BG vs Sub");
    }

    // Pad 3: TMC average pi0 mass histogram (hMean)
    cMerged->cd(3);
    double maxMeanCanvas = 1.2 * std::max(1.0, hMean->GetMaximum());
    TH2F* frame3 = new TH2F("frame3", "TMC #pi^{0} Mass (Mean);M_{#gamma#gamma} (GeV);Counts",
                            10, 0.0, 0.3, 10, 0, maxMeanCanvas);
    frame3->Draw("AXIS");
    hMean->SetLineColor(kBlack);
    hMean->SetLineWidth(2);
    hMean->SetMarkerStyle(0);
    hMean->Draw("L same");
    {
        TLatex lat3;
        lat3.SetNDC();
        lat3.SetTextSize(0.04);
        lat3.DrawLatex(0.15, 0.92, "Toy MC Mean Pi0 Mass");
    }

    // Pad 4: TMC RMS histogram (hVar)
    cMerged->cd(4);
    double maxVarCanvas = 1.2 * hVar->GetMaximum();
    TH2F* frame4 = new TH2F("frame4", "TMC RMS;M_{#gamma#gamma} (GeV);#sigma(Counts)",
                            10, 0.0, 0.3, 10, 0, maxVarCanvas);
    frame4->Draw("AXIS");
    hVar->SetLineColor(kRed);
    hVar->SetLineWidth(2);
    hVar->SetMarkerStyle(0);
    hVar->Draw("L same");
    {
        TLatex lat4;
        lat4.SetNDC();
        lat4.SetTextSize(0.04);
        lat4.DrawLatex(0.15, 0.92, "Toy MC RMS (BG Variation)");
    }

    // Pad 5: Final pi0 mass with a shaded error band
    cMerged->cd(5);
    double max5Canvas = 1.2 * std::max(1.0, (hMean->GetMaximum() + hVar->GetMaximum()));
    TH2F* frame5 = new TH2F("frame5", "Pi0 Mass + TMC Error Region;M_{#gamma#gamma} (GeV);Counts",
                            10, 0.0, 0.3, 10, 0, max5Canvas);
    frame5->Draw("AXIS");
    grErrorPoly->Draw("F same");
    hMean->SetLineColor(kBlack);
    hMean->SetMarkerStyle(0);
    hMean->Draw("L same");
    {
        TLatex lat5;
        lat5.SetNDC();
        lat5.SetTextSize(0.04);
        lat5.DrawLatex(0.15, 0.92, "Final Pi0 Mass w/ TMC Error Region");
    }

    cMerged->Print(outputFileName.c_str());
    std::cout << "All done. Wrote combined 5-pad canvas to " << outputFileName << "\n";
    return 0;
}
