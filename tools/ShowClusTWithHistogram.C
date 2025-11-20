#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TPaveStats.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TPaveText.h>

void ShowClusTWithHistogram(Int_t nrun, Int_t firstevent, Int_t lastevent, Double_t fitRangeMin = 150.0, Double_t fitRangeMax = 160.0) {
    // Define input file path
    std::string inputFileName = TString::Format("/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%i_%i_%i.root", nrun, firstevent, lastevent).Data();
    std::cout << "Input file: " << inputFileName << std::endl;

    // Hardcoded PNG directory (update as needed)
    std::string pngDir = "/volatile/hallc/nps/jpcrafts/Plots"; // Specify your desired directory here
    gSystem->mkdir(pngDir.c_str(), true);

    // Generate the output PNG file name
    std::string pngFilePath = TString::Format("%s/output_timing_%i_%i_%i.png", pngDir.c_str(), nrun, firstevent, lastevent).Data();
    std::cout << "Output file: " << pngFilePath << std::endl;

    // Open the ROOT file
    TFile* file = TFile::Open(inputFileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << inputFileName << std::endl;
        return;
    }
    std::cout << "File opened successfully: " << inputFileName << std::endl;

    // Access the TTree
    TTree* tree = dynamic_cast<TTree*>(file->Get("T"));
    if (!tree) {
        std::cerr << "Error: Cannot find TTree in file." << std::endl;
        file->Close();
        return;
    }
    std::cout << "TTree loaded successfully." << std::endl;

    // Set up branches
    Double_t nclust;
    Double_t clusT[10000];

    tree->SetBranchAddress("NPS.cal.nclust", &nclust);
    tree->SetBranchAddress("NPS.cal.clusT", clusT);

    // Create a histogram for cluster times
    TH1F* histClusT = new TH1F("histClusT", "Cluster Times;Time [ns];Counts", 100, 100, 200);

    // Event loop
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Number of entries in the tree: " << nEntries << std::endl;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // Progress counter
        if (i % 10000 == 0) {
            std::cout << "Processed " << i << " / " << nEntries << " events." << std::endl;
        }

        // Fill histogram with cluster times for each event
        for (Int_t j = 0; j < nclust; j++) {
            histClusT->Fill(clusT[j]);
        }
    }

    std::cout << "Histogram filled with cluster times." << std::endl;

    // Fit the background using a polynomial (degree 1), excluding the region of interest
    TF1* backgroundFunc = new TF1("backgroundFunc", "pol2", 117, 184); // Define the function
    backgroundFunc->SetRange(117, fitRangeMin); // Left range
    histClusT->Fit(backgroundFunc, "R+");
    TF1* backgroundFuncRight = new TF1("backgroundFuncRight", "pol2", fitRangeMax, 184); // Right range
    histClusT->Fit(backgroundFuncRight, "R+");

    // Subtract the background
    TH1F* histClusTNoBg = (TH1F*)histClusT->Clone("histClusTNoBg");
    for (int bin = 1; bin <= histClusT->GetNbinsX(); bin++) {
        Double_t binContent = histClusT->GetBinContent(bin);
        Double_t binCenter = histClusT->GetBinCenter(bin);
        Double_t bgValue = (binCenter < fitRangeMin) ? backgroundFunc->Eval(binCenter) : backgroundFuncRight->Eval(binCenter);
        histClusTNoBg->SetBinContent(bin, binContent - bgValue);
    }

    // Fit the Gaussian to the subtracted histogram
    TF1* fitFunc = new TF1("fitFunc", "gaus", fitRangeMin, fitRangeMax);
    histClusTNoBg->Fit(fitFunc, "R");

    // Extract fit parameters
    Double_t mean = fitFunc->GetParameter(1);
    Double_t sigma = fitFunc->GetParameter(2);
    std::cout << "Fit Results (after background subtraction): Mean = " << mean
              << " ns, Sigma = " << sigma << " ns" << std::endl;

    // Draw the histogram with the fit
    TCanvas* canvas = new TCanvas("canvas", "Cluster Times", 800, 600);
    histClusT->SetLineColor(kBlue); // Original histogram
    histClusTNoBg->SetLineColor(kRed); // Background-subtracted histogram
    histClusT->Draw("hist");
    histClusTNoBg->Draw("hist same");
    backgroundFunc->SetLineColor(kGreen);
    backgroundFunc->Draw("same");
    backgroundFuncRight->SetLineColor(kGreen);
    backgroundFuncRight->Draw("same");
    fitFunc->Draw("same");

    // Add mean as a text box above the Gaussian peak
    Double_t textX1 = mean - 5;
    Double_t textX2 = mean + 5;
    Double_t textY1 = histClusTNoBg->GetMaximum() * 0.75;
    Double_t textY2 = histClusTNoBg->GetMaximum() * 0.85;
    TPaveText* paveText = new TPaveText(textX1, textY1, textX2, textY2, "br");
    paveText->AddText(TString::Format("Mean: %.3f ns", mean));
    paveText->SetFillColor(0);
    paveText->SetTextSize(0.03);
    paveText->SetBorderSize(0);
    paveText->Draw();

    // Save the histogram to a file
    canvas->SaveAs(pngFilePath.c_str());

    std::cout << "Histogram saved to: " << pngFilePath << std::endl;

    // Clean up
    delete canvas;
    delete histClusT;
    delete histClusTNoBg;
    delete fitFunc;
    delete backgroundFunc;
    delete backgroundFuncRight;
    delete paveText;
    file->Close();
    delete file;
}

// Usage: Pass the run number, first event, last event, and optionally the fit range as arguments to the script
int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <run number> <first event> <last event> [fit range min] [fit range max]" << std::endl;
        return 1;
    }

    Int_t nrun = atoi(argv[1]);
    Int_t firstevent = atoi(argv[2]);
    Int_t lastevent = atoi(argv[3]);

    Double_t fitRangeMin = (argc > 4) ? atof(argv[4]) : 150.0;
    Double_t fitRangeMax = (argc > 5) ? atof(argv[5]) : 160.0;

    ShowClusTWithHistogram(nrun, firstevent, lastevent, fitRangeMin, fitRangeMax);

    return 0;
}
