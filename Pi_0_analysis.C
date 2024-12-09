#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h> // For drawing histograms
#include <TPaveText.h> // Include for TPaveText
#include <TLegend.h>
#include <TText.h>   // For adding cut information to plots
#include <TBranch.h>
#include <TSystem.h>
#include <TString.h>
#include <TCut.h> // For handling cuts
#include <TTreeFormula.h> // For evaluating cuts
#include <iostream>
#include <map>
#include <string>
#include <memory>

// Helper function to extract the base name of a file
std::string getBaseName(const std::string& filePath) {
    size_t pos = filePath.find_last_of("/\\");
    return (pos == std::string::npos) ? filePath : filePath.substr(pos + 1);
}

void processRootFile(const std::string& inputFile, const std::string& outputDir,
                     const std::map<std::string, std::tuple<std::string, std::string, int, float, float>>& branchHistMap,
                     const TCut& cut) {
    // Open the input ROOT file
    std::unique_ptr<TFile> inputRoot(TFile::Open(inputFile.c_str(), "READ"));
    if (!inputRoot || inputRoot->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputFile << std::endl;
        return;
    }

    // Get the TTree from the input file
    TTree* tree = dynamic_cast<TTree*>(inputRoot->Get("T")); // Replace "T" with your actual TTree name
    if (!tree) {
        std::cerr << "Error: TTree not found in the input file." << std::endl;
        return;
    }

    // Define the cuts using TCut
    TCut cutFormula = cut;

    // Prepare the TTreeFormula to evaluate the cut
    TTreeFormula formula("cut_formula", cutFormula.GetTitle(), tree);
    if (formula.GetNdim() == 0) {
        std::cerr << "Error: Invalid cut formula: " << cutFormula.GetTitle() << std::endl;
        return;
    }

    // Create histograms
    std::map<std::string, std::unique_ptr<TH1F>> histograms;
    for (const auto& [branchName, histConfig] : branchHistMap) {
        const auto& [histName, histTitle, bins, xMin, xMax] = histConfig;
        histograms[branchName] = std::make_unique<TH1F>(histName.c_str(), histTitle.c_str(), bins, xMin, xMax);
    }

    // Set branch addresses
    std::map<std::string, Double_t> branchValues;
    tree->SetBranchStatus("*", 0); // Disable all branches by default
    for (const auto& [branchName, _] : branchHistMap) {
        tree->SetBranchStatus(branchName.c_str(), 1); // Enable only branches in branchHistMap
        branchValues[branchName] = 0.0f;
        tree->SetBranchAddress(branchName.c_str(), &branchValues[branchName]);
    }

    // Loop over events and fill histograms if the cut is satisfied
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Check if the cut condition is satisfied
        if (!formula.EvalInstance()) continue;

        // Fill histograms
        for (const auto& [branchName, hist] : histograms) {
            hist->Fill(branchValues[branchName]);
        }

        if (i % 10000 == 0) std::cout << "Processed " << i << " events..." << std::endl;
    }

    // Create the output directory
    gSystem->mkdir(outputDir.c_str(), true);

    // Extract the base name of the input file (without the directory or extension)
    std::string baseName = getBaseName(inputFile);
    baseName = TString(baseName).ReplaceAll(".root", "").Data(); // Remove the .root extension

    // Create the output ROOT file
    std::string outputFilePath = outputDir + "/" + baseName + ".root";
    std::unique_ptr<TFile> outputRoot(TFile::Open(outputFilePath.c_str(), "RECREATE"));
    if (!outputRoot || outputRoot->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outputFilePath << std::endl;
        return;
    }

    // Write histograms to the output ROOT file
    for (const auto& [_, hist] : histograms) {
        hist->Write();
    }

    // Hardcoded PDF directory (update as needed)
    std::string pdfDir = "/volatile/hallc/nps/jpcrafts/Plots"; // Specify your desired directory here
    gSystem->mkdir(pdfDir.c_str(), true);

    // Generate the output PDF file name
    std::string pdfFilePath = pdfDir + "/" + baseName + ".pdf";

    // Create a canvas for plotting
    TCanvas canvas("canvas", "Histograms", 800, 600);
    canvas.SetLeftMargin(0.15); // Expand canvas to the left
    canvas.SetBottomMargin(0.15); // Increase bottom margin to ensure enough space
    canvas.Print((pdfFilePath + "[").c_str()); // Start the PDF file

    // Draw each histogram with cut information and save to PDF
    for (const auto& [branchName, hist] : histograms) {
        hist->SetStats(0); // Disable the default statistics box
        hist->Draw();

        // Create and move the legend to the upper-right corner, including the number of events
        TLegend legend(0.7, 0.75, 0.9, 0.9); // Adjust position to the upper-right corner
        legend.AddEntry(hist.get(), TString::Format("#splitline{%s}{Events: %d}", branchName.c_str(), static_cast<int>(hist->GetEntries())), "l");
        legend.Draw();

        // Add cut information in a resized text box above the x-axis label
        TPaveText textBox(0.1, 0.03, 0.9, 0.08, "NDC"); // Adjusted position and size
        textBox.AddText(cut.GetTitle());
        textBox.SetTextSize(0.03); // Reduced font size
        textBox.SetFillColor(0);   // Transparent background
        textBox.SetBorderSize(0);  // Remove border
        textBox.Draw();

        // Save the histogram to the PDF
        canvas.Print(pdfFilePath.c_str());
    }

    // Close the PDF file
    canvas.Print((pdfFilePath + "]").c_str());

    std::cout << "Histograms saved to " << outputFilePath << " and " << pdfFilePath << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <nrun> <firstevent> <lastevent>" << std::endl;
        return 1;
    }

    // Convert command-line arguments
    int nrun = std::stoi(argv[1]);
    int firstevent = std::stoi(argv[2]);
    int lastevent = std::stoi(argv[3]);

    // Define input file, output directory, and branches with histogram configurations
    std::string inputFile = TString::Format("/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%i_%i_%i.root", nrun, firstevent, lastevent).Data();
    std::string outputDir = "/volatile/hallc/nps/jpcrafts/ROOTfiles";

    std::map<std::string, std::tuple<std::string, std::string, int, float, float>> branchHistMap = {
        {"H.cal.etottracknorm", {"h_etottracknorm", "HMS_Etottracknorm;E/p;nevents", 150, 0.01, 1.51}},
        {"H.gtr.dp", {"h_gtr_delta", "HMS_GoldenTrack_Delta;", 160, -8.0, 8.0}},
        // Add more branches as needed
    };

    // Define the cut condition using TCut
    TCut cut = "(H.cal.etottracknorm > 0.1 && H.cal.etottracknorm < 1.3) && (H.gtr.dp > -5.0 && H.gtr.dp < 5.0)";

    // Process the file with cuts
    processRootFile(inputFile, outputDir, branchHistMap, cut);

    return 0;
}
