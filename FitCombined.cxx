#include <iostream>
#include <string>
#include <cstdlib> // for std::stod
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <cmath>  // if needed

// A helper function that returns the "basename" without the path, 
// e.g. from "/path/to/foo.root" -> "foo.root"
static std::string getBaseName(const std::string &path)
{
    // find the last slash
    size_t pos = path.find_last_of("/\\");
    if(pos == std::string::npos) {
        // no slash found, return entire string
        return path;
    }
    // substring after last slash
    return path.substr(pos+1);
}

// A helper that strips off ".root" if present
static std::string stripRootExtension(const std::string &fname)
{
    if(fname.size() > 5) {
        // check if ends with ".root"
        if(fname.compare(fname.size()-5, 5, ".root") == 0) {
            return fname.substr(0, fname.size()-5);
        }
    }
    return fname;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <combined-file.root> [fitMin=0.12] [fitMax=0.15]\n";
        return 1;
    }

    // 1) Parse arguments
    std::string fileName = argv[1];
    double fitMin = 0.12;
    double fitMax = 0.145;
    if (argc > 2) fitMin = std::stod(argv[2]);
    if (argc > 3) fitMax = std::stod(argv[3]);

    // Directory where we want the PNG
    std::string pngDir = "/volatile/hallc/nps/jpcrafts/Plots";
    // Ensure it exists
    gSystem->mkdir(pngDir.c_str(), true);

    // 2) Open the ROOT file
    TFile* f = TFile::Open(fileName.c_str(), "READ");
    if(!f || f->IsZombie()) {
        std::cerr << "Error: cannot open file " << fileName << std::endl;
        return 1;
    }

    // 3) Retrieve histogram "hMgg"
    TH1F* hMgg = (TH1F*) f->Get("hMgg");
    if(!hMgg) {
        std::cerr << "Error: histogram 'hMgg' not found in " << fileName << std::endl;
        f->Close();
        return 1;
    }

    // 4) Setup style so we see fit parameters in stats box
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    // 5) Make a canvas
    TCanvas* c1 = new TCanvas("c1","Combine Hist Fit",800,600);
    hMgg->Draw();

    // 6) Define a Gaussian function in [fitMin..fitMax]
    TF1* fitGauss = new TF1("fitGauss","gaus", fitMin, fitMax);
    fitGauss->SetLineColor(kRed);

    // Optionally set initial guesses for amplitude, mean, sigma
    fitGauss->SetParameters(1000.0, 0.135, 0.01);

    // 7) Fit
    hMgg->Fit(fitGauss, "R");

    // 8) Construct the PNG path in /volatile/hallc/nps/jpcrafts/Plots
    //    Use the "basename" stripped of path and .root, then "_fit.png"
    std::string base = getBaseName(fileName);      // e.g. "foo.root"
    base = stripRootExtension(base);               // e.g. "foo"
    // Now append "_fit.png"
    std::string outPng = pngDir + "/" + base + "_fit.png";

    // 9) Save the canvas as a PNG in /Plots
    c1->SaveAs(outPng.c_str());

    std::cout << "Fitted histogram from " << fileName
              << " in [" << fitMin << ".." << fitMax << "]"
              << "\nPNG saved to: " << outPng << std::endl;

    // Cleanup
    f->Close();
    delete c1;
    delete fitGauss;

    return 0;
}
