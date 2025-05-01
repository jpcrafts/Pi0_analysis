#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TF1.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

// --------------------------------------------------------------------
// 1) Normalization & smoothing (unchanged from before).
// --------------------------------------------------------------------
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
                if (y[p] > y[last]){
                    filteredPeaks.back() = p;
                }
            } else {
                filteredPeaks.push_back(p);
            }
        }
    }
    if (filteredPeaks.empty()) return yNorm;

    double sumPeaks = 0.0;
    for (int idx : filteredPeaks){
        sumPeaks += y[idx];
    }
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
        if (y[i] < y[i-1] && y[i] < y[i+1] && y[i] <= globalMin + troughTolerance){
            candidateTroughs.push_back(i);
        }
    }
    std::vector<int> filteredTroughs;
    for (int t : candidateTroughs){
        if (filteredTroughs.empty()){
            filteredTroughs.push_back(t);
        } else {
            int last = filteredTroughs.back();
            if (t - last < troughFilterWindow){
                if (y[t] < y[last]){
                    filteredTroughs.back() = t;
                }
            } else {
                filteredTroughs.push_back(t);
            }
        }
    }
    if (filteredTroughs.empty()) return yNorm;

    double sumTroughs = 0.0;
    for (int idx : filteredTroughs){
        sumTroughs += y[idx];
    }
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
        double sum  = 0.0;
        int count = 0;
        for (int j = i - window/2; j <= i + window/2; j++){
            if (j >= 0 && j < n){
                sum += data[j];
                count++;
            }
        }
        smoothed[i] = (count>0 ? sum / count : data[i]);
    }
    return smoothed;
}

// --------------------------------------------------------------------
// 2) Final cluster cut for the pi0 mass (AFTER background subtraction).
//    Now includes boundary cuts to exclude outer 1-block layer. // NEW
// --------------------------------------------------------------------

// Suppose your total array is ±31.32 cm in x, ±37.80 cm in y
// One block ~2.16 cm => we exclude that from each edge => ±29.16, ±35.64
static bool isGoodCluster(double e, double t, double x, double y)
{
    // Basic energy/time cuts
    bool passBasic = (e >= 0.6 && t >= 149.0 && t <= 151.0);

    // One-block offset for the bounding box // NEW
    double xInnerMin = -29.16;
    double xInnerMax = +29.16;
    double yInnerMin = -35.64;
    double yInnerMax = +35.64;

    bool insideEdges = (x > xInnerMin && x < xInnerMax &&
                        y > yInnerMin && y < yInnerMax);

    return (passBasic && insideEdges);
}

int main(int argc, char* argv[])
{
    if (argc < 3){
        std::cerr << "Usage: " << argv[0] << " <input.root> <output.png>\n";
        return 1;
    }
    std::string inputFileName  = argv[1];
    std::string outputFileName = argv[2];

    // ---------------------------------------------------------------------
    // HMS-level cuts (typo HPC -> HMS).
    // ---------------------------------------------------------------------
    double edtmtdcCut   = 0.1;
    double hdeltaLowCut = -8.5, hdeltaHighCut = 8.5;
    double hcaltotCut   = 0.6;
    double hcernpeCut   = 1.0;
    double gtrthCut     = 0.09;
    double gtrphCut     = 0.09;

    // BG region
    double bgLow  = 113.0;
    double bgHigh = 142.5;
    int    nBins  = 650; // bin width ~0.0454 ns

    // Signal region (shifted region)
    double shiftVal = 28.05;
    double sigLow   = 141.789;
    double sigHigh  = 171.289;

    // ---------------------------------------------------------------------
    // Define a lambda for the HMS cuts // NEW
    // ---------------------------------------------------------------------
    auto passHMSCuts = [&](double edtmtdc, double hdelta, double hcaltot,
                           double hcernpe, double gtrth, double gtrph)
    {
        if (edtmtdc >= edtmtdcCut) return false;
        if (hdelta < hdeltaLowCut || hdelta > hdeltaHighCut) return false;
        if (hcaltot <= hcaltotCut) return false;
        if (hcernpe <= hcernpeCut) return false;
        if (std::fabs(gtrth) > gtrthCut) return false;
        if (std::fabs(gtrph) > gtrphCut) return false;
        return true;
    };

    // Open file once
    TFile* inFile = TFile::Open(inputFileName.c_str(), "READ");
    if (!inFile || inFile->IsZombie()){
        std::cerr << "Error: cannot open " << inputFileName << std::endl;
        return 1;
    }
    TTree* tree = dynamic_cast<TTree*>(inFile->Get("T"));
    if (!tree){
        std::cerr << "Error: TTree 'T' not found.\n";
        inFile->Close();
        return 1;
    }

    // Branch setup
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    tree->SetBranchStatus("H.gtr.dp",              1);
    tree->SetBranchStatus("H.cal.etotnorm",        1);
    tree->SetBranchStatus("H.cer.npeSum",          1);
    tree->SetBranchStatus("H.gtr.th",              1);
    tree->SetBranchStatus("H.gtr.ph",              1);
    tree->SetBranchStatus("NPS.cal.newEWClusT",    1);

    double edtmtdc = 0, hdelta = 0, hcaltot = 0, hcernpe = 0, gtrth = 0, gtrph = 0;
    std::vector<double>* vecT = nullptr;

    tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
    tree->SetBranchAddress("H.gtr.dp",               &hdelta);
    tree->SetBranchAddress("H.cal.etotnorm",         &hcaltot);
    tree->SetBranchAddress("H.cer.npeSum",           &hcernpe);
    tree->SetBranchAddress("H.gtr.th",               &gtrth);
    tree->SetBranchAddress("H.gtr.ph",               &gtrph);
    tree->SetBranchAddress("NPS.cal.newEWClusT",     &vecT);

    // ---------------------------------
    // A) Build BG region histogram [113,142.5] (HMS only)
    // ---------------------------------
    TH1F* hBG = new TH1F("hBG", "BG region", nBins, bgLow, bgHigh);
    hBG->SetDirectory(nullptr);

    Long64_t nAll = tree->GetEntries();
    for (Long64_t i = 0; i < nAll; i++){
        tree->GetEntry(i);

        if (!passHMSCuts(edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph)) continue;
        if (!vecT) continue;

        for (double tVal : *vecT){
            if (tVal >= bgLow && tVal <= bgHigh){
                hBG->Fill(tVal);
            }
        }
    }

    // ---------------------------------
    // B) Fit + Norm + Smooth for BG region
    // ---------------------------------
    std::vector<double> xBG(nBins), yBG(nBins);
    for (int b = 1; b <= nBins; b++){
        xBG[b-1] = hBG->GetBinCenter(b);
        yBG[b-1] = hBG->GetBinContent(b);
    }
    TGraph* grBG = new TGraph(nBins, &xBG[0], &yBG[0]);
    grBG->Sort();
    TSpline3* splineBG = new TSpline3("splineBG", grBG);

    std::vector<double> ySpline(nBins);
    for (int i = 0; i < nBins; i++){
        ySpline[i] = splineBG->Eval(xBG[i]);
    }

    int peakWindow   = 10;  double minPeakFrac = 0.75;
    int troughWindow = 15;  double troughFrac  = 0.05;
    int swPeak       = 2,   swTrough          = 5;

    auto yPeaks   = normalizePeaks(ySpline, peakWindow, minPeakFrac);
    auto yTroughs = normalizeTroughs(yPeaks,  troughWindow, troughFrac);
    auto yFinal   = regionSmooth(yTroughs,   swPeak,       swTrough);

    // ---------------------------------
    // C) Build Raw Signal histogram [141.789,171.289] (HMS only)
    // ---------------------------------
    TH1F* hSignalRaw = new TH1F("hSignalRaw", "Signal Region", nBins, sigLow, sigHigh);
    hSignalRaw->SetDirectory(nullptr);

    for (Long64_t i = 0; i < nAll; i++){
        tree->GetEntry(i);

        if (!passHMSCuts(edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph)) continue;
        if (!vecT) continue;

        for (double tVal : *vecT){
            if (tVal >= sigLow && tVal <= sigHigh){
                hSignalRaw->Fill(tVal);
            }
        }
    }

    // ---------------------------------
    // D) Shift BG + subtract in the signal region
    // ---------------------------------
    std::vector<double> xShift(nBins), yShift(nBins);
    for (int i = 0; i < nBins; i++){
        xShift[i] = xBG[i] + shiftVal;
        yShift[i] = yFinal[i];
    }
    TGraph* grBGshift = new TGraph(nBins, &xShift[0], &yShift[0]);

    TH1F* hSub = (TH1F*)hSignalRaw->Clone("hSub");
    hSub->SetDirectory(nullptr);
    hSub->SetTitle("Subtracted Signal");
    for (int b = 1; b <= nBins; b++){
        double binCenter = hSignalRaw->GetBinCenter(b);
        double dataVal   = hSignalRaw->GetBinContent(b);
        double bgVal     = grBGshift->Eval(binCenter);
        hSub->SetBinContent(b, dataVal - bgVal);
    }

    // ---------------------------------
    // E) Final pass: cluster-level cuts + compute Pi0 mass
    // ---------------------------------
    tree->SetBranchStatus("NPS.cal.nclust", 1);
    tree->SetBranchStatus("NPS.cal.clusE",  1);
    tree->SetBranchStatus("NPS.cal.clusT",  1);
    tree->SetBranchStatus("NPS.cal.clusX",  1);
    tree->SetBranchStatus("NPS.cal.clusY",  1);

    double nclustDouble = 0.0;
    static const int MAXC = 10000;
    double clusE[MAXC], clusT[MAXC], clusX[MAXC], clusY[MAXC];

    tree->SetBranchAddress("NPS.cal.nclust", &nclustDouble);
    tree->SetBranchAddress("NPS.cal.clusE", clusE);
    tree->SetBranchAddress("NPS.cal.clusT", clusT);
    tree->SetBranchAddress("NPS.cal.clusX", clusX);
    tree->SetBranchAddress("NPS.cal.clusY", clusY);

    // Pi0 mass histograms
    TH1F* hMgg = new TH1F("hMgg", "Pi0 Mass (Weighted);M_{#gamma#gamma};Counts", 200, 0.0, 0.3);
    hMgg->SetDirectory(nullptr);
    TH1F* hMggNoWeight = new TH1F("hMggNoWeight", "Pi0 Mass (Unweighted);M_{#gamma#gamma};Counts", 200, 0.0, 0.3);
    hMggNoWeight->SetDirectory(nullptr);

    double DNPS = 407.0; // some NPS distance

    for (Long64_t i = 0; i < nAll; i++){
        tree->GetEntry(i);

        // HPC/HMS-level cuts
        if (!passHMSCuts(edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph)) continue;

        int nclust = (int)nclustDouble;
        if (nclust < 2) continue;

        // Build a list of clusters passing final cuts
        struct Cinfo { int idx; double weight; };
        std::vector<Cinfo> keep;
        keep.reserve(nclust);

        for (int cID = 0; cID < nclust; cID++){
            double tC = clusT[cID];
            if (tC < sigLow || tC > sigHigh) continue;

            int binID   = hSub->FindBin(tC);
            double dataVal = hSignalRaw->GetBinContent(binID);
            if (dataVal <= 0) continue;
            double subVal  = hSub->GetBinContent(binID);
            if (subVal <= 0) continue;

            double frac = subVal / dataVal;
            if (frac > 1.0) frac = 1.0;

            // boundary check // NEW
            if (!isGoodCluster(clusE[cID], tC, clusX[cID], clusY[cID])) continue;

            keep.push_back({cID, frac});
        }
        if (keep.size() < 2) continue;

        // Now form ALL pairs among keep // CHANGED
        for (size_t iC = 0; iC < keep.size(); iC++){
            for (size_t jC = iC+1; jC < keep.size(); jC++){
                int id1 = keep[iC].idx;
                int id2 = keep[jC].idx;
                double w1 = keep[iC].weight;
                double w2 = keep[jC].weight;

                double E1 = clusE[id1];
                double E2 = clusE[id2];
                double x1 = clusX[id1];
                double y1 = clusY[id1];
                double x2 = clusX[id2];
                double y2 = clusY[id2];

                double dx = x1 - x2;
                double dy = y1 - y2;
                double dist = std::sqrt(dx*dx + dy*dy);

                double theta = dist / DNPS;
                double s2    = std::sin(0.5 * theta);
                double Mgg   = std::sqrt(4.0 * E1 * E2 * s2 * s2);

                double pairWeight = w1 * w2;
                hMgg->Fill(Mgg, pairWeight);
                hMggNoWeight->Fill(Mgg);
            }
        }
    }

    inFile->Close();

    // ---------------------------------------------------------------------
    // F) Plotting: 4 vertical pads
    // ---------------------------------------------------------------------
    TCanvas* cMerged = new TCanvas("cMerged", "Combined Plots", 1600, 2400);
    cMerged->Divide(1, 4);

    // Pad 1: BG region + final smoothed spline
    cMerged->cd(1);
    grBG->SetMarkerStyle(20);
    grBG->SetMarkerSize(0.7);
    grBG->Draw("AP");
    TGraph* grBGfinal = new TGraph(nBins);
    for (int i = 0; i < nBins; i++){
        grBGfinal->SetPoint(i, xBG[i], yFinal[i]);
    }
    grBGfinal->SetLineColor(kBlue);
    grBGfinal->SetLineWidth(2);
    grBGfinal->Draw("L same");
    {
        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextSize(0.04);
        latex1.DrawLatex(0.15, 0.92, "BG Region [113,142.5] + Spline Fit");
    }

    // Pad 2: Signal region (raw, shifted BG, subtracted)
    cMerged->cd(2);
    double maxSig = hSignalRaw->GetMaximum();
    double maxBG  = 0.0;
    for (int b = 1; b <= nBins; b++){
        double diff = hSignalRaw->GetBinContent(b) - hSub->GetBinContent(b);
        if (diff > maxBG) maxBG = diff;
    }
    double maxVal2 = std::max(maxSig, maxBG);
    TH2F* frame2 = new TH2F("frame2", "Signal Region;Time (ns);Counts", 
                            10, sigLow, sigHigh, 10, 0, 1.2*maxVal2);
    frame2->Draw("AXIS");

    TGraph* grSignalRaw = new TGraph(hSignalRaw);
    grSignalRaw->SetLineColor(kBlack);
    grSignalRaw->SetLineWidth(2);
    grSignalRaw->Draw("L same");

    TGraph* grShifted = new TGraph(nBins);
    for (int i = 0; i < nBins; i++){
        double xx  = hSignalRaw->GetBinCenter(i+1);
        double val = grBGshift->Eval(xx);
        grShifted->SetPoint(i, xx, val);
    }
    grShifted->SetLineColor(kRed);
    grShifted->SetLineWidth(2);
    grShifted->Draw("L same");

    TGraph* grSub = new TGraph(hSub);
    grSub->SetLineColor(kGreen+2);
    grSub->SetLineWidth(2);
    grSub->Draw("L same");

    {
        TLatex latex2;
        latex2.SetNDC();
        latex2.SetTextSize(0.04);
        latex2.DrawLatex(0.15, 0.92, "Signal Region [141.789,171.289] Subtraction");
    }

    // Pad 3: Unweighted π⁰ mass
    cMerged->cd(3);
    double maxMassUnw = 1.2 * std::max(1.0, hMggNoWeight->GetMaximum());
    TH2F* frame3 = new TH2F("frame3", "Pi0 Mass (Unweighted);M_{#gamma#gamma} (GeV);Counts", 
                            10, 0.0, 0.3, 10, 0, maxMassUnw);
    frame3->Draw("AXIS");
    hMggNoWeight->SetLineColor(kBlack);
    hMggNoWeight->SetLineWidth(2);
    hMggNoWeight->Draw("same");

    {
        TLatex latex3;
        latex3.SetNDC();
        latex3.SetTextSize(0.04);
        latex3.DrawLatex(0.15, 0.92, "Unweighted Pi0 Mass");
    }

    // Pad 4: Weighted π⁰ mass
    cMerged->cd(4);
    std::cout << "DEBUG: hMgg entries = " << hMgg->GetEntries()
              << ", integral = " << hMgg->Integral() << std::endl;
    double maxMass = 1.2 * std::max(1.0, hMgg->GetMaximum());
    TH2F* frame4 = new TH2F("frame4", "Pi0 Mass (Weighted);M_{#gamma#gamma} (GeV);Counts",
                            10, 0.0, 0.3, 10, 0, maxMass);
    frame4->Draw("AXIS");
    hMgg->SetLineColor(kBlack);
    hMgg->SetLineWidth(2);
    hMgg->Draw("same");
    gStyle->SetOptStat(1111);

    // Example Gaussian fit near ~135 MeV
    TF1* fitG = new TF1("fitG", "gaus", 0.122, 0.141);
    fitG->SetLineColor(kRed);
    hMgg->Fit(fitG, "R");

    {
        TLatex latex4;
        latex4.SetNDC();
        latex4.SetTextSize(0.04);
        latex4.DrawLatex(0.15, 0.92, "Weighted Pi0 Mass");
    }

    cMerged->Print(outputFileName.c_str());
    std::cout << "All done. Wrote combined 4-pad canvas to " << outputFileName << "\n";

    return 0;
}
