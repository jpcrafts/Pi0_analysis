#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_spline.h>

// We focus only on TSpline3 and TSpline5

// Standard moving-average smoothing (if needed)
std::vector<double> movingAverage(const std::vector<double>& data, int smoothWindow) {
    int n = data.size();
    std::vector<double> smoothed(n);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        int count = 0;
        for (int j = i - smoothWindow/2; j <= i + smoothWindow/2; j++) {
            if (j >= 0 && j < n) {
                sum += data[j];
                count++;
            }
        }
        smoothed[i] = sum / count;
    }
    return smoothed;
}

// Region-specific moving average smoothing:
// Uses one window size for peak regions and a different window for trough regions.
// It uses the midpoint (average of global min and max) as the threshold.
std::vector<double> regionSpecificMovingAverage(const std::vector<double>& data, int smoothWindowPeak, int smoothWindowTrough) {
    int n = data.size();
    std::vector<double> smoothed(n);
    double globalMax = *std::max_element(data.begin(), data.end());
    double globalMin = *std::min_element(data.begin(), data.end());
    double midValue = (globalMax + globalMin) / 2.0;
    for (int i = 0; i < n; i++) {
        int window = (data[i] > midValue) ? smoothWindowPeak : smoothWindowTrough;
        double sum = 0.0;
        int count = 0;
        for (int j = i - window/2; j <= i + window/2; j++) {
            if (j >= 0 && j < n) {
                sum += data[j];
                count++;
            }
        }
        smoothed[i] = sum / count;
    }
    return smoothed;
}

int main(int argc, char* argv[]){
    if(argc < 3){
        std::cerr << "Usage: " << argv[0] << " <input.root> <output.png>\n";
        return 1;
    }
    std::string inputFileName  = argv[1];
    std::string outputFileName = argv[2]; // e.g., "multisplinefit.png"

    // 1. Open the ROOT file and get TTree "T"
    TFile* inFile = TFile::Open(inputFileName.c_str(), "READ");
    if(!inFile || inFile->IsZombie()){
        std::cerr << "Error: cannot open " << inputFileName << std::endl;
        return 1;
    }
    TTree* tree = dynamic_cast<TTree*>(inFile->Get("T"));
    if(!tree){
        std::cerr << "Error: TTree 'T' not found in file " << inputFileName << std::endl;
        inFile->Close();
        return 1;
    }

    // 2. Set branch addresses for cuts and data
    double edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph;
    std::vector<double>* vecNewEWClusT = nullptr;
    tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
    tree->SetBranchAddress("H.gtr.dp",               &hdelta);
    tree->SetBranchAddress("H.cal.etotnorm",         &hcaltot);
    tree->SetBranchAddress("H.cer.npeSum",           &hcernpe);
    tree->SetBranchAddress("H.gtr.th",               &gtrth);
    tree->SetBranchAddress("H.gtr.ph",               &gtrph);
    // Data branch for fitting:
    tree->SetBranchAddress("NPS.cal.newEWClusT",     &vecNewEWClusT);

    // 3. Define region and histogram parameters
    double regionLow  = 113.0;
    double regionHigh = 142.5;
    int nBins = 650;  // one point per bin
    TH1F* hRegion = new TH1F("hRegion", "Region for Spline Fit", nBins, regionLow, regionHigh);
    hRegion->SetDirectory(nullptr);  // detach so it persists after file close

    // Example cut thresholds
    double edtmtdcCut = 0.1;
    double hdeltaLowCut = -8.5;
    double hdeltaHighCut = 8.5;
    double hcaltotCut = 0.6;
    double hcernpeCut = 1.0;
    double gtrthCut = 0.09;
    double gtrphCut = 0.09;

    // 4. Fill histogram from the TTree
    Long64_t nEntries = tree->GetEntries();
    for(Long64_t i = 0; i < nEntries; i++){
        tree->GetEntry(i);
        if(edtmtdc >= edtmtdcCut)       continue;
        if(hdelta < hdeltaLowCut)       continue;
        if(hdelta > hdeltaHighCut)      continue;
        if(hcaltot <= hcaltotCut)       continue;
        if(hcernpe <= hcernpeCut)       continue;
        if(fabs(gtrth) > gtrthCut)       continue;
        if(fabs(gtrph) > gtrphCut)       continue;
        if(!vecNewEWClusT)              continue;
        for(auto &t : *vecNewEWClusT){
            if(t >= regionLow && t <= regionHigh)
                hRegion->Fill(t);
        }
    }
    inFile->Close();

    // 5. Convert histogram to (x,y) vectors (one per bin)
    std::vector<double> xvals, yvals;
    xvals.reserve(nBins);
    yvals.reserve(nBins);
    for(int b = 1; b <= nBins; b++){
        double x = hRegion->GetBinCenter(b);
        double y = hRegion->GetBinContent(b);
        xvals.push_back(x);
        yvals.push_back(y);
    }
    double sumY = 0.0;
    for(double v : yvals) sumY += v;
    std::cout << "DEBUG: sum of y-values = " << sumY << std::endl;
    if(sumY == 0.0){
        std::cerr << "ERROR: All bins are zero.\n";
        delete hRegion;
        return 1;
    }

    // 6. Sort data and remove duplicate x-values
    std::vector<std::pair<double,double>> dataPoints;
    dataPoints.reserve(xvals.size());
    for(size_t i = 0; i < xvals.size(); i++){
        dataPoints.push_back({xvals[i], yvals[i]});
    }
    std::sort(dataPoints.begin(), dataPoints.end());
    for(size_t i = 0; i < dataPoints.size(); i++){
        xvals[i] = dataPoints[i].first;
        yvals[i] = dataPoints[i].second;
    }
    std::vector<double> filtered_x, filtered_y;
    filtered_x.reserve(xvals.size());
    filtered_y.reserve(yvals.size());
    filtered_x.push_back(xvals[0]);
    filtered_y.push_back(yvals[0]);
    for(size_t i = 1; i < xvals.size(); i++){
        if(xvals[i] != xvals[i-1]){
            filtered_x.push_back(xvals[i]);
            filtered_y.push_back(yvals[i]);
        }
    }
    if(filtered_x.size() < 2){
        std::cerr << "ERROR: Not enough unique x-values for interpolation.\n";
        delete hRegion;
        return 1;
    }

    // 7. Create evaluation grid (nEval = 650 points)
    int nEval = 650;
    double dx = (regionHigh - regionLow) / (nEval - 1);
    std::vector<double> xSpline(nEval);
    for(int i = 0; i < nEval; i++){
        xSpline[i] = regionLow + i * dx;
    }

    // 8. Create TGraph from filtered data and build TSpline3 and TSpline5
    TGraph* grRegion = new TGraph(filtered_x.size(), &filtered_x[0], &filtered_y[0]);
    grRegion->SetName("grRegion");
    grRegion->SetTitle("Data Region");
    grRegion->Sort();
    TSpline3* spline3 = new TSpline3("spline3", grRegion);
    TSpline5* spline5 = new TSpline5("spline5", grRegion);

    // 9. Evaluate the original splines on the evaluation grid
    std::vector<double> ySpline3_vals(nEval), ySpline5_vals(nEval);
    for (int i = 0; i < nEval; i++){
        double xx = xSpline[i];
        ySpline3_vals[i] = spline3->Eval(xx);
        ySpline5_vals[i] = spline5->Eval(xx);
    }

    // 10. Peak normalization parameters for peaks:
    int peakFilterWindow = 10;        // window for merging nearby peaks
    double minPeakFraction = 0.75;       // only consider peaks above 50% of max

    auto normalizePeaks = [nEval, peakFilterWindow, minPeakFraction](const std::vector<double>& y) -> std::vector<double> {
        std::vector<double> yNorm = y; // Copy original data
        std::vector<int> candidatePeaks;
        double maxVal = *std::max_element(y.begin(), y.end());
        double minPeakHeight = minPeakFraction * maxVal;
        // Detect candidate peaks: local maximum above threshold
        for (int i = 1; i < nEval - 1; i++){
            if(y[i] > y[i-1] && y[i] > y[i+1] && y[i] > minPeakHeight)
                candidatePeaks.push_back(i);
        }
        // Filter nearby peaks: within peakFilterWindow indices, keep only the highest.
        std::vector<int> filteredPeaks;
        for (int p : candidatePeaks) {
            if(filteredPeaks.empty()){
                filteredPeaks.push_back(p);
            } else {
                int last = filteredPeaks.back();
                if(p - last < peakFilterWindow){
                    if(y[p] > y[last])
                        filteredPeaks.back() = p;
                } else {
                    filteredPeaks.push_back(p);
                }
            }
        }
        if(filteredPeaks.empty()) return yNorm;
        double sumPeaks = 0.0;
        for (int idx : filteredPeaks)
            sumPeaks += y[idx];
        double avgPeak = sumPeaks / filteredPeaks.size();
        // Adjust values around each filtered peak in a window of size peakFilterWindow
        for (int p : filteredPeaks) {
            double scale = avgPeak / y[p];
            for (int i = std::max(0, p - peakFilterWindow/2); i <= std::min(nEval - 1, p + peakFilterWindow/2); i++){
                double weight = 1.0 - fabs(i - p) / double(peakFilterWindow + 1);
                yNorm[i] = y[i] * (1 - weight) + (y[i] * scale) * weight;
            }
        }
        return yNorm;
    };

    // 11. Trough normalization parameters for troughs:
    int troughFilterWindow = 15;       // window for merging nearby troughs (separate from peaks)
    double troughToleranceFactor = 0.05; // Use 5% of data range as tolerance
    auto normalizeTroughs = [nEval, troughFilterWindow, troughToleranceFactor](const std::vector<double>& y) -> std::vector<double> {
        std::vector<double> yNorm = y; // Copy original data
        std::vector<int> candidateTroughs;
        double globalMin = *std::min_element(y.begin(), y.end());
        double globalMax = *std::max_element(y.begin(), y.end());
        double troughTolerance = troughToleranceFactor * (globalMax - globalMin);
        // Consider troughs that are within troughTolerance of the global minimum.
        for (int i = 1; i < nEval - 1; i++){
            if(y[i] < y[i-1] && y[i] < y[i+1] && y[i] <= globalMin + troughTolerance)
                candidateTroughs.push_back(i);
        }
        // Filter nearby troughs: within troughFilterWindow, keep only the deepest trough.
        std::vector<int> filteredTroughs;
        for (int t : candidateTroughs) {
            if(filteredTroughs.empty()){
                filteredTroughs.push_back(t);
            } else {
                int last = filteredTroughs.back();
                if(t - last < troughFilterWindow){
                    if(y[t] < y[last])
                        filteredTroughs.back() = t;
                } else {
                    filteredTroughs.push_back(t);
                }
            }
        }
        if(filteredTroughs.empty()) return yNorm;
        double sumTroughs = 0.0;
        for (int idx : filteredTroughs)
            sumTroughs += y[idx];
        double avgTrough = sumTroughs / filteredTroughs.size();
        // Adjust values around each trough within troughFilterWindow region.
        for (int t : filteredTroughs) {
            double scale = avgTrough / y[t];
            for (int i = std::max(0, t - troughFilterWindow/2); i <= std::min(nEval - 1, t + troughFilterWindow/2); i++){
                double weight = 1.0 - fabs(i - t) / double(troughFilterWindow + 1);
                yNorm[i] = y[i] * (1 - weight) + (y[i] * scale) * weight;
            }
        }
        return yNorm;
    };

    // Apply normalizations sequentially: first peaks then troughs.
    std::vector<double> ySpline3_norm = normalizeTroughs(normalizePeaks(ySpline3_vals));
    std::vector<double> ySpline5_norm = normalizeTroughs(normalizePeaks(ySpline5_vals));

    // 12. Apply region-specific moving-average smoothing to the normalized curves.
    // We now allow separate smoothing for peak and trough regions.
    int smoothWindowPeak = 2;   // Smoothing window for peaks
    int smoothWindowTrough = 5; // Smoothing window for troughs
    auto regionSmooth = [=](const std::vector<double>& data) -> std::vector<double> {
        int n = data.size();
        std::vector<double> smoothed(n);
        double globalMax = *std::max_element(data.begin(), data.end());
        double globalMin = *std::min_element(data.begin(), data.end());
        double midValue = (globalMax + globalMin) / 2.0;
        for (int i = 0; i < n; i++){
            int window = (data[i] > midValue) ? smoothWindowPeak : smoothWindowTrough;
            double sum = 0.0;
            int count = 0;
            for (int j = i - window/2; j <= i + window/2; j++){
                if(j >= 0 && j < n){
                    sum += data[j];
                    count++;
                }
            }
            smoothed[i] = sum / count;
        }
        return smoothed;
    };
    std::vector<double> ySpline3_final = regionSmooth(ySpline3_norm);
    std::vector<double> ySpline5_final = regionSmooth(ySpline5_norm);

    // 13. Create TGraph objects for original, normalized (unsmoothed), and final (normalized+smoothed) curves.
    TGraph* grSpline3Original = new TGraph(nEval, &xSpline[0], &ySpline3_vals[0]);
    TGraph* grSpline3Normalized = new TGraph(nEval, &xSpline[0], &ySpline3_norm[0]);
    TGraph* grSpline3Final = new TGraph(nEval, &xSpline[0], &ySpline3_final[0]);

    TGraph* grSpline5Original = new TGraph(nEval, &xSpline[0], &ySpline5_vals[0]);
    TGraph* grSpline5Normalized = new TGraph(nEval, &xSpline[0], &ySpline5_norm[0]);
    TGraph* grSpline5Final = new TGraph(nEval, &xSpline[0], &ySpline5_final[0]);

    // Set line styles: original normalized as dashed, final (normalized+smoothed) as solid.
    grSpline3Original->SetLineColor(kRed);
    grSpline3Original->SetLineWidth(2);
    grSpline3Original->SetLineStyle(2);  // dashed
    grSpline3Normalized->SetLineColor(kBlue);
    grSpline3Normalized->SetLineWidth(2);
    grSpline3Normalized->SetLineStyle(2);  // dashed (optional)
    grSpline3Final->SetLineColor(kGreen);
    grSpline3Final->SetLineWidth(2);
    grSpline3Final->SetLineStyle(1);  // solid

    grSpline5Original->SetLineColor(kRed);
    grSpline5Original->SetLineWidth(2);
    grSpline5Original->SetLineStyle(2);  // dashed
    grSpline5Normalized->SetLineColor(kBlue);
    grSpline5Normalized->SetLineWidth(2);
    grSpline5Normalized->SetLineStyle(2);  // dashed (optional)
    grSpline5Final->SetLineColor(kGreen);
    grSpline5Final->SetLineWidth(2);
    grSpline5Final->SetLineStyle(1);  // solid

    // 14. Create one canvas with two vertically arranged pads for comparison
    TCanvas* cAll = new TCanvas("cAll", "Spline Comparison with Peak & Trough Normalization and Region-Specific Smoothing", 1200, 1600);
    cAll->Divide(1, 2);
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);

    // Pad 1: TSpline3 – show original, normalized (unsmoothed), and final (normalized+smoothed)
    cAll->cd(1);
    grRegion->Draw("AP");
    grSpline3Original->Draw("L same");
    grSpline3Normalized->Draw("L same");
    grSpline3Final->Draw("L same");
    {
      TLegend* leg = new TLegend(0.55, 0.7, 0.9, 0.9);
      leg->AddEntry(grRegion, "Data", "P");
      leg->AddEntry(grSpline3Original, "TSpline3 (Original)", "L");
      leg->AddEntry(grSpline3Normalized, "TSpline3 (Normalized)", "L");
      leg->AddEntry(grSpline3Final, "TSpline3 (Final)", "L");
      leg->Draw();
      latex.DrawLatex(0.2, 0.92, "TSpline3 - Cubic Spline");
    }

    // Pad 2: TSpline5 – show original, normalized (unsmoothed), and final (normalized+smoothed)
    cAll->cd(2);
    grRegion->Draw("AP");
    grSpline5Original->Draw("L same");
    grSpline5Normalized->Draw("L same");
    grSpline5Final->Draw("L same");
    {
      TLegend* leg = new TLegend(0.55, 0.7, 0.9, 0.9);
      leg->AddEntry(grRegion, "Data", "P");
      leg->AddEntry(grSpline5Original, "TSpline5 (Original)", "L");
      leg->AddEntry(grSpline5Normalized, "TSpline5 (Normalized)", "L");
      leg->AddEntry(grSpline5Final, "TSpline5 (Final)", "L");
      leg->Draw();
      latex.DrawLatex(0.2, 0.92, "TSpline5 - Akima Spline");
    }

    // 15. Save the canvas as a PNG file
    cAll->Print(outputFileName.c_str());

    return 0;
}
