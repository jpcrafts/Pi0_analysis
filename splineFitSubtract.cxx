#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
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

// -------------------------------------------------------------
// 1) HELPER FUNCTIONS for smoothing & normalization
//    (identical or nearly so to your prior code).
// -------------------------------------------------------------

// Identify & scale peak regions to a common average peak height.
std::vector<double> normalizePeaks(const std::vector<double>& y, int peakFilterWindow,
                                   double minPeakFraction)
{
    int n = y.size();
    std::vector<double> yNorm = y;
    std::vector<int> candidatePeaks;
    double maxVal = *std::max_element(y.begin(), y.end());
    double minPeakHeight = minPeakFraction * maxVal;

    // Find local maxima above minPeakHeight
    for (int i = 1; i < n - 1; i++){
        if(y[i] > y[i-1] && y[i] > y[i+1] && y[i] > minPeakHeight){
            candidatePeaks.push_back(i);
        }
    }

    // Merge close peaks: if two peaks are within peakFilterWindow,
    // keep only the larger of the two
    std::vector<int> filteredPeaks;
    for (int p : candidatePeaks){
        if(filteredPeaks.empty()){
            filteredPeaks.push_back(p);
        } else {
            int last = filteredPeaks.back();
            if(p - last < peakFilterWindow){
                if(y[p] > y[last]) {
                    filteredPeaks.back() = p;
                }
            } else {
                filteredPeaks.push_back(p);
            }
        }
    }

    if(filteredPeaks.empty()) {
        return yNorm;
    }

    // Compute average peak height of the final filtered peaks
    double sumPeaks = 0.0;
    for (int idx : filteredPeaks){
        sumPeaks += y[idx];
    }
    double avgPeak = sumPeaks / filteredPeaks.size();

    // Scale each peak region to match the average peak height
    for (int p : filteredPeaks){
        double scale = avgPeak / y[p];
        int pStart = std::max(0, p - peakFilterWindow/2);
        int pEnd   = std::min(n - 1, p + peakFilterWindow/2);
        for (int i = pStart; i <= pEnd; i++){
            double weight = 1.0 - (std::fabs(i - p) / double(peakFilterWindow + 1));
            yNorm[i] = y[i] * (1 - weight) + (y[i] * scale) * weight;
        }
    }

    return yNorm;
}

// Identify & scale trough regions to a common average trough height.
std::vector<double> normalizeTroughs(const std::vector<double>& y, int troughFilterWindow,
                                     double troughToleranceFactor)
{
    int n = y.size();
    std::vector<double> yNorm = y;
    std::vector<int> candidateTroughs;

    double globalMin = *std::min_element(y.begin(), y.end());
    double globalMax = *std::max_element(y.begin(), y.end());
    double troughTolerance = troughToleranceFactor * (globalMax - globalMin);

    // Find local minima near globalMin
    for (int i = 1; i < n - 1; i++){
        if(y[i] < y[i-1] && y[i] < y[i+1] && y[i] <= globalMin + troughTolerance){
            candidateTroughs.push_back(i);
        }
    }

    // Merge close troughs: if two are within troughFilterWindow,
    // keep only the deeper one
    std::vector<int> filteredTroughs;
    for (int t : candidateTroughs){
        if(filteredTroughs.empty()){
            filteredTroughs.push_back(t);
        } else {
            int last = filteredTroughs.back();
            if(t - last < troughFilterWindow){
                if(y[t] < y[last]) {
                    filteredTroughs.back() = t;
                }
            } else {
                filteredTroughs.push_back(t);
            }
        }
    }

    if(filteredTroughs.empty()) {
        return yNorm;
    }

    // Compute average trough
    double sumTroughs = 0.0;
    for (int idx : filteredTroughs){
        sumTroughs += y[idx];
    }
    double avgTrough = sumTroughs / filteredTroughs.size();

    // Scale each trough region to match the average trough height
    for (int t : filteredTroughs){
        double scale = avgTrough / y[t];
        int tStart = std::max(0, t - troughFilterWindow/2);
        int tEnd   = std::min(n - 1, t + troughFilterWindow/2);
        for (int i = tStart; i <= tEnd; i++){
            double weight = 1.0 - (std::fabs(i - t) / double(troughFilterWindow + 1));
            yNorm[i] = y[i] * (1 - weight) + (y[i] * scale) * weight;
        }
    }

    return yNorm;
}

// Region-specific smoothing: smaller window for "peak" values, larger window for "trough" values
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
        int    count = 0;
        for (int j = i - window/2; j <= i + window/2; j++){
            if(j >= 0 && j < n){
                sum += data[j];
                count++;
            }
        }
        smoothed[i] = (count > 0) ? (sum / count) : data[i];
    }
    return smoothed;
}

// -------------------------------------------------------------
// 2) MAIN FUNCTION
// -------------------------------------------------------------
int main(int argc, char* argv[]){
    if(argc < 3){
        std::cerr << "Usage: " << argv[0] << " <input.root> <output.png>\n";
        return 1;
    }
    std::string inputFileName  = argv[1];
    std::string outputFileName = argv[2];

    // ---------------------------------------------------------
    // A) Build background region histogram: [113, 142.5], 650 bins
    // ---------------------------------------------------------
    double bgLow  = 113.0;
    double bgHigh = 142.5;
    int    nBins  = 650; // => width ~ 29.5 ns / 650 ~ 0.0453846 ns

    TFile* inFile = TFile::Open(inputFileName.c_str(),"READ");
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

    // Example branch addresses + cuts (adjust as you like)
    double edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph;
    std::vector<double>* vecNewEWClusT = nullptr;
    tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
    tree->SetBranchAddress("H.gtr.dp",               &hdelta);
    tree->SetBranchAddress("H.cal.etotnorm",         &hcaltot);
    tree->SetBranchAddress("H.cer.npeSum",           &hcernpe);
    tree->SetBranchAddress("H.gtr.th",               &gtrth);
    tree->SetBranchAddress("H.gtr.ph",               &gtrph);
    tree->SetBranchAddress("NPS.cal.newEWClusT",     &vecNewEWClusT);

    // Example cut thresholds
    double edtmtdcCut   = 0.1;
    double hdeltaLowCut = -8.5, hdeltaHighCut = 8.5;
    double hcaltotCut   = 0.6;
    double hcernpeCut   = 1.0;
    double gtrthCut     = 0.09;
    double gtrphCut     = 0.09;

    // Create background histogram
    TH1F* hRegion = new TH1F("hRegion", "Background Region", nBins, bgLow, bgHigh);
    hRegion->SetDirectory(nullptr);

    Long64_t nEntries = tree->GetEntries();
    for(Long64_t i = 0; i < nEntries; i++){
        tree->GetEntry(i);

        // apply cuts
        if(edtmtdc >= edtmtdcCut)          continue;
        if(hdelta < hdeltaLowCut)          continue;
        if(hdelta > hdeltaHighCut)         continue;
        if(hcaltot <= hcaltotCut)          continue;
        if(hcernpe <= hcernpeCut)          continue;
        if(std::fabs(gtrth) > gtrthCut)    continue;
        if(std::fabs(gtrph) > gtrphCut)    continue;

        if(!vecNewEWClusT) continue;
        for(double t : *vecNewEWClusT){
            if(t >= bgLow && t <= bgHigh){
                hRegion->Fill(t);
            }
        }
    }
    inFile->Close();

    // Convert hRegion to vectors
    std::vector<double> xvals(nBins), yvals(nBins);
    for(int b = 1; b <= nBins; b++){
        xvals[b-1] = hRegion->GetBinCenter(b);
        yvals[b-1] = hRegion->GetBinContent(b);
    }

    // Build a TGraph & Spline
    TGraph* grRegion = new TGraph(nBins, &xvals[0], &yvals[0]);
    grRegion->SetName("grRegion");
    grRegion->SetTitle("Data Region");
    grRegion->Sort();

    TSpline3* spline3 = new TSpline3("spline3", grRegion);
    std::vector<double> ySpline(nBins);
    for(int i = 0; i < nBins; i++){
        ySpline[i] = spline3->Eval(xvals[i]);
    }

    // ---------------------------------------------------------
    // B) Normalization + smoothing on the spline
    // ---------------------------------------------------------
    // 1) Peak normalization, 2) trough normalization, 3) region-based smoothing.
    int    peakFilterWindow      = 10;
    double minPeakFraction       = 0.75;  // peaks above 75% of max
    int    troughFilterWindow    = 15;
    double troughToleranceFactor = 0.05;  // 5% of data range
    int    smoothWindowPeak      = 2;
    int    smoothWindowTrough    = 5;

    std::vector<double> yPeaksNorm   = normalizePeaks(ySpline, peakFilterWindow, minPeakFraction);
    std::vector<double> yTroughsNorm = normalizeTroughs(yPeaksNorm, troughFilterWindow, troughToleranceFactor);
    std::vector<double> yFinal       = regionSmooth(yTroughsNorm, smoothWindowPeak, smoothWindowTrough);

    // Build a TGraph for the final background
    TGraph* grBGfinal = new TGraph(nBins, &xvals[0], &yFinal[0]);
    grBGfinal->SetName("grBGfinal");
    grBGfinal->SetTitle("Final Background (after normalization + smoothing)");

    // ---------------------------------------------------------
    // C) Shift the background by +28.015 ns (pure horizontal)
    // ---------------------------------------------------------
    double shiftVal = 28.05;
    std::vector<double> xShifted(nBins), yShifted(nBins);  // same y, shifted x
    for(int i = 0; i < nBins; i++){
        xShifted[i] = xvals[i] + shiftVal; 
        yShifted[i] = yFinal[i];
    }
    TGraph* grBGshifted = new TGraph(nBins, &xShifted[0], &yShifted[0]);
    grBGshifted->SetName("grBGshifted");
    grBGshifted->SetTitle("Shifted Background (+28.015 ns)");

    // ---------------------------------------------------------
    // D) Subtraction region: [141.789, 171.289], also 650 bins
    // ---------------------------------------------------------
    double newLow  = 141.789;
    double newHigh = 171.289;
    // same 650 bins => same bin width (29.5 / 650)
    TH1F* hSignal = new TH1F("hSignal", "Signal in [141.789,171.289]", nBins, newLow, newHigh);
    hSignal->SetDirectory(nullptr);

    // Reopen file to fill from TTree for the new region
    TFile* inFile2 = TFile::Open(inputFileName.c_str(),"READ");
    if(!inFile2 || inFile2->IsZombie()){
        std::cerr << "Error: cannot reopen " << inputFileName << std::endl;
        return 1;
    }
    TTree* tree2 = dynamic_cast<TTree*>(inFile2->Get("T"));
    if(!tree2){
        std::cerr << "Error: TTree 'T' not found in file " << inputFileName << std::endl;
        inFile2->Close();
        return 1;
    }
    // Re-set branches
    tree2->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
    tree2->SetBranchAddress("H.gtr.dp",               &hdelta);
    tree2->SetBranchAddress("H.cal.etotnorm",         &hcaltot);
    tree2->SetBranchAddress("H.cer.npeSum",           &hcernpe);
    tree2->SetBranchAddress("H.gtr.th",               &gtrth);
    tree2->SetBranchAddress("H.gtr.ph",               &gtrph);
    tree2->SetBranchAddress("NPS.cal.newEWClusT",     &vecNewEWClusT);

    Long64_t nEntries2 = tree2->GetEntries();
    for(Long64_t i = 0; i < nEntries2; i++){
        tree2->GetEntry(i);
        // same cuts
        if(edtmtdc >= edtmtdcCut)          continue;
        if(hdelta < hdeltaLowCut)          continue;
        if(hdelta > hdeltaHighCut)         continue;
        if(hcaltot <= hcaltotCut)          continue;
        if(hcernpe <= hcernpeCut)          continue;
        if(std::fabs(gtrth) > gtrthCut)    continue;
        if(std::fabs(gtrph) > gtrphCut)    continue;

        if(!vecNewEWClusT) continue;
        for(double t : *vecNewEWClusT){
            if(t >= newLow && t <= newHigh){
                hSignal->Fill(t);
            }
        }
    }
    inFile2->Close();

    // Build x array & evaluate the shifted BG in the new region
    std::vector<double> xNew(nBins), yBG(nBins), ySub(nBins);
    for(int b = 1; b <= nBins; b++){
        double binCenter = hSignal->GetBinCenter(b);
        xNew[b-1]  = binCenter;
        double bgVal = grBGshifted->Eval(binCenter); // Evaluate SHIFTED background
        double sigVal= hSignal->GetBinContent(b);
        yBG[b-1]   = bgVal;
        ySub[b-1]  = sigVal - bgVal;
    }

    // Make TGraphs for plotting
    TGraph* grData     = new TGraph(hSignal);  // data in new region
    TGraph* grBGregion = new TGraph(nBins, &xNew[0], &yBG[0]);
    TGraph* grSubtract = new TGraph(nBins, &xNew[0], &ySub[0]);

    grData->SetLineColor(kBlack);
    grData->SetLineWidth(2);
    grBGregion->SetLineColor(kRed);
    grBGregion->SetLineWidth(2);
    grSubtract->SetLineColor(kGreen+3);
    grSubtract->SetLineWidth(2);

    // ---------------------------------------------------------
    // E) Draw everything
    // ---------------------------------------------------------
    TCanvas* cMerged = new TCanvas("cMerged", "Spline Fit + Shifted Subtraction", 1600, 1200);
    cMerged->Divide(1,2);

    // Left pad: show background region with final (smoothed) spline
    cMerged->cd(1);
    grRegion->SetMarkerStyle(20);
    grRegion->SetMarkerSize(0.7);
    grRegion->Draw("AP");
    grBGfinal->SetLineColor(kBlue);
    grBGfinal->SetLineWidth(2);
    grBGfinal->Draw("L same");

    TLegend* leg1 = new TLegend(0.55, 0.7, 0.9, 0.9);
    leg1->AddEntry(grRegion,  "Background Data", "P");
    leg1->AddEntry(grBGfinal, "Final Spline Fit", "L");
    leg1->Draw();

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.04);
    latex1.DrawLatex(0.15, 0.92, "Background Region [113,142.5], nBins=650");

    // Right pad: new region [141.789,171.289] - show Data, Shifted BG, Subtracted
    cMerged->cd(2);
    double maxSig = hSignal->GetMaximum();
    double maxBG  = *std::max_element(yBG.begin(), yBG.end());
    double maxVal = std::max(maxSig, maxBG);
    TH2F* frame = new TH2F("frame","",10,newLow,newHigh,10,0,1.1*maxVal);
    frame->SetXTitle("Time (ns)");
    frame->SetYTitle("Counts");
    frame->Draw("AXIS");

    grData->Draw("L same");
    grBGregion->Draw("L same");
    grSubtract->Draw("L same");

    TLegend* leg2 = new TLegend(0.55, 0.65, 0.9, 0.85);
    leg2->AddEntry(grData,     "Data (141.789-171.289)", "L");
    leg2->AddEntry(grBGregion, "Shifted BG (+28.015 ns)", "L");
    leg2->AddEntry(grSubtract, "Subtracted", "L");
    leg2->Draw();

    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.04);
    latex2.DrawLatex(0.15, 0.92, "Signal Region [141.789,171.289], nBins=650");

    // Save
    cMerged->Print(outputFileName.c_str());

    return 0;
}
