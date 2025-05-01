/*************************************************************
  dummy_debug.cpp

  This minimal ROOT-based program processes a dummy file to:

    1. Fill two raw histograms (upstream and downstream) for the signal time window.
    2. Scale them using the dummy effective charge and fixed scale factors.
    3. Combine the two histograms into a normalized dummy histogram (in counts/µC).
    4. Save output PNG files of each histogram for inspection.

  Usage:
    ./dummy_debug dummyFile.root dummyEffectiveCharge

  Example:
    ./dummy_debug VolatileROOTfiles/dummy_x58_q51_p5_merged.root 254792.874

  Compile example:
    g++ -O2 -std=c++17 -o dummy_debug dummy_debug.cpp `root-config --cflags --libs`
*************************************************************/

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib> // for std::stod

// A simple function to implement the HMS cuts.
// Adjust cut values as needed.
bool passHMSCuts(double edtmtdc, double hdelta, double hcaltot, double hcernpe, double gtrth, double gtrph) {
    const double edtmtdcCut = 0.1;
    const double hdeltaLowCut = -8.5, hdeltaHighCut = 8.5;
    const double hcaltotCut = 0.6, hcernpeCut = 1.0;
    const double gtrthCut = 0.09, gtrphCut = 0.09;
    if(edtmtdc >= edtmtdcCut) return false;
    if(hdelta < hdeltaLowCut || hdelta > hdeltaHighCut) return false;
    if(hcaltot <= hcaltotCut) return false;
    if(hcernpe <= hcernpeCut) return false;
    if(std::fabs(gtrth) > gtrthCut) return false;
    if(std::fabs(gtrph) > gtrphCut) return false;
    return true;
}

int main(int argc, char* argv[]){
    if(argc < 3){
        std::cerr << "Usage: " << argv[0] << " <dummyFile.root> <dummyEffectiveCharge>\n";
        return 1;
    }

    std::string dummyFileName = argv[1];
    double dummyEffectiveCharge = std::stod(argv[2]);  // in µC

    // Hard-coded parameters for dummy processing:
    const int nBins = 650;
    const double sigLow = 141.789, sigHigh = 171.289; // signal time window
    const double scale_factor_dummy_upstream   = 8.467;
    const double scale_factor_dummy_downstream = 4.256;

    // Open the dummy file.
    TFile* dummyFile = TFile::Open(dummyFileName.c_str(), "READ");
    if(!dummyFile || dummyFile->IsZombie()){
        std::cerr << "Error: could not open dummy file " << dummyFileName << "\n";
        return 1;
    }

    TTree* dummyTree = dynamic_cast<TTree*>(dummyFile->Get("T"));
    if(!dummyTree){
        std::cerr << "Error: TTree 'T' not found in dummy file.\n";
        dummyFile->Close();
        return 1;
    }

    // Enable only needed branches.
    dummyTree->SetBranchStatus("*", 0);
    dummyTree->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    dummyTree->SetBranchStatus("H.gtr.dp",              1);
    dummyTree->SetBranchStatus("H.cal.etotnorm",        1);
    dummyTree->SetBranchStatus("H.cer.npeSum",          1);
    dummyTree->SetBranchStatus("H.gtr.th",              1);
    dummyTree->SetBranchStatus("H.gtr.ph",              1);
    dummyTree->SetBranchStatus("H.gtr.y",               1);
    dummyTree->SetBranchStatus("NPS.cal.clusT",         1);

    // Set branch addresses.
    double edtmtdc_d = 0, hdelta_d = 0, hcaltot_d = 0, hcernpe_d = 0;
    double gtrth_d = 0, gtrph_d = 0, gtry_d = 0;
    double clusT_dummy = 0;  // assume one cluster time per event in dummy
    dummyTree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc_d);
    dummyTree->SetBranchAddress("H.gtr.dp",               &hdelta_d);
    dummyTree->SetBranchAddress("H.cal.etotnorm",         &hcaltot_d);
    dummyTree->SetBranchAddress("H.cer.npeSum",           &hcernpe_d);
    dummyTree->SetBranchAddress("H.gtr.th",               &gtrth_d);
    dummyTree->SetBranchAddress("H.gtr.ph",               &gtrph_d);
    dummyTree->SetBranchAddress("H.gtr.y",                &gtry_d);
    dummyTree->SetBranchAddress("NPS.cal.clusT",          &clusT_dummy);

    // Create raw dummy histograms.
    TH1F* hDummy_up_raw = new TH1F("hDummy_up_raw", "Dummy Upstream (raw)", nBins, sigLow, sigHigh);
    hDummy_up_raw->SetDirectory(nullptr);
    TH1F* hDummy_down_raw = new TH1F("hDummy_down_raw", "Dummy Downstream (raw)", nBins, sigLow, sigHigh);
    hDummy_down_raw->SetDirectory(nullptr);

    // Loop over dummy entries.
    Long64_t nDummy = dummyTree->GetEntries();
    for(Long64_t i = 0; i < nDummy; i++){
        dummyTree->GetEntry(i);
        if(!passHMSCuts(edtmtdc_d, hdelta_d, hcaltot_d, hcernpe_d, gtrth_d, gtrph_d))
            continue;
        // Fill only if clusT_dummy is within the signal window.
        if(clusT_dummy < sigLow || clusT_dummy > sigHigh)
            continue;
        if(gtry_d > 0)
            hDummy_down_raw->Fill(clusT_dummy);
        else
            hDummy_up_raw->Fill(clusT_dummy);
    }

    // Print out dummy histogram binning info.
    std::cout << "hDummy_up_raw: " << hDummy_up_raw->GetNbinsX() << " bins, range ["
              << hDummy_up_raw->GetXaxis()->GetXmin() << ", " << hDummy_up_raw->GetXaxis()->GetXmax() << "]\n";
    std::cout << "hDummy_down_raw: " << hDummy_down_raw->GetNbinsX() << " bins, range ["
              << hDummy_down_raw->GetXaxis()->GetXmin() << ", " << hDummy_down_raw->GetXaxis()->GetXmax() << "]\n";

    // Save the raw dummy histograms for debugging.
    TCanvas* cDummyRaw = new TCanvas("cDummyRaw", "Raw Dummy Histograms", 800, 600);
    hDummy_up_raw->SetLineColor(kBlue);
    hDummy_up_raw->Draw();
    cDummyRaw->SaveAs("dummy_up_raw.png");
    hDummy_down_raw->SetLineColor(kRed);
    hDummy_down_raw->Draw("same");
    cDummyRaw->SaveAs("dummy_updown_raw.png");

    // Normalize dummy histograms (to counts per µC).
    hDummy_up_raw->Scale(1.0 / (dummyEffectiveCharge * 8.467));
    hDummy_down_raw->Scale(1.0 / (dummyEffectiveCharge * 4.256));

    // Combine normalized dummy histograms.
    TH1F* hDummy_norm = (TH1F*) hDummy_up_raw->Clone("hDummy_norm");
    hDummy_norm->Add(hDummy_down_raw);

    // Debug print for normalized dummy.
    std::cout << "hDummy_norm: " << hDummy_norm->GetNbinsX() << " bins, range ["
              << hDummy_norm->GetXaxis()->GetXmin() << ", " << hDummy_norm->GetXaxis()->GetXmax() << "]\n";

    // Save the normalized dummy histogram.
    TCanvas* cDummyNorm = new TCanvas("cDummyNorm", "Normalized Dummy", 800, 600);
    hDummy_norm->SetLineColor(kGreen+2);
    hDummy_norm->Draw();
    cDummyNorm->SaveAs("dummy_norm.png");

    dummyFile->Close();

    // End dummy-only debug program.
    std::cout << "Dummy histograms have been processed and saved as PNGs.\n";
    gSystem->Exit(0);
}
// End of dummy_debug.cpp
// *************************************************************/
//
// This program is a minimal example to demonstrate the dummy histogram processing.
// It does not include the full functionality of the original code.
// The original code is much more complex and includes additional steps for
// processing real data, performing mass calculations, and generating plots.    