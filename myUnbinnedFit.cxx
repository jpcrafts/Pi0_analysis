#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <input.root> <output.pdf>\n";
        return 1;
    }
    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];

    // --------------------------------------------------------------------------
    // 1. Open the ROOT file and retrieve TTree "T"
    // --------------------------------------------------------------------------
    TFile *inFile = TFile::Open(inputFileName.c_str(), "READ");
    if (!inFile || inFile->IsZombie())
    {
        std::cerr << "Error: Cannot open file " << inputFileName << std::endl;
        return 1;
    }
    TTree *tree = dynamic_cast<TTree *>(inFile->Get("T"));
    if (!tree)
    {
        std::cerr << "Error: TTree 'T' not found in file " << inputFileName << std::endl;
        inFile->Close();
        return 1;
    }
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Loaded TTree 'T' from " << inputFileName
              << " with " << nEntries << " entries.\n";

    // --------------------------------------------------------------------------
    // 2. Set branch addresses for cuts & cluster times
    // --------------------------------------------------------------------------
    double edtmtdc;                               // T.hms.hEDTM_tdcTimeRaw
    double hdelta;                                // H.gtr.dp
    double hcaltot;                               // H.cal.etotnorm
    double hcernpe;                               // H.cer.npeSum
    double gtrth;                                 // H.gtr.th
    double gtrph;                                 // H.gtr.ph
    double nclust_d;                              // NPS.cal.nclust
    std::vector<double> *vecNewEWClusT = nullptr; // NPS.cal.newEWClusT

    tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
    tree->SetBranchAddress("H.gtr.dp", &hdelta);
    tree->SetBranchAddress("H.cal.etotnorm", &hcaltot);
    tree->SetBranchAddress("H.cer.npeSum", &hcernpe);
    tree->SetBranchAddress("H.gtr.th", &gtrth);
    tree->SetBranchAddress("H.gtr.ph", &gtrph);
    tree->SetBranchAddress("NPS.cal.nclust", &nclust_d);
    tree->SetBranchAddress("NPS.cal.newEWClusT", &vecNewEWClusT);

    // --------------------------------------------------------------------------
    // 3. Define your cut thresholds
    // --------------------------------------------------------------------------
    const double edtmtdcCut = 0.1;
    const double hdeltaLow = -8.5;
    const double hdeltaHigh = 8.5;
    const double hcaltotCut = 0.6;
    const double hcernpeCut = 1.0;
    const double gtrthCut = 0.09;
    const double gtrphCut = 0.09;
    // const int    maxClusters= 10;   // example limit

    // --------------------------------------------------------------------------
    // 4. Create a RooDataSet for the unbinned fit:
    //    We'll define "x" from [80,220] ns to capture the region of interest.
    // --------------------------------------------------------------------------
    RooRealVar x("x", "Cluster Time (ns)", 110, 190);
    RooArgSet varSet(x);
    RooDataSet data("data", "Unbinned dataset", varSet);

    Long64_t passCount = 0;
    for (Long64_t i = 0; i < nEntries; i++)
    {
        tree->GetEntry(i);

        // Apply your posted cuts:
        if (edtmtdc >= edtmtdcCut)
            continue;
        if (hdelta <= hdeltaLow)
            continue;
        if (hdelta >= hdeltaHigh)
            continue;
        if (hcaltot <= hcaltotCut)
            continue;
        if (hcernpe <= hcernpeCut)
            continue;
        if (std::fabs(gtrth) > gtrthCut)
            continue;
        if (std::fabs(gtrph) > gtrphCut)
            continue;

        // int nclust = static_cast<int>(nclust_d);
        // if (nclust > maxClusters)          continue; // or just skip if > max

        // Now each event can have multiple cluster times in vecNewEWClusT
        if (!vecNewEWClusT)
            continue; // safety check
        for (auto &timeVal : *vecNewEWClusT)
        {
            // We only add times that fall within [80,220]
            if (timeVal < 80 || timeVal > 220)
                continue;

            x.setVal(timeVal);
            data.add(varSet);
            passCount++;
        }
    }
    inFile->Close();

    std::cout << "After cuts, stored " << passCount
              << " cluster times in the RooDataSet.\n";

    // --------------------------------------------------------------------------
    // 5. Build the Constrained Summation Model:
    //    One large Gaussian around ~150 ns, plus a sum of many small Gaussians
    //    from 100 ns to 172 ns in steps of 2 ns, all equally weighted.
    // --------------------------------------------------------------------------

    // --- The Big Peak
    RooRealVar bigMean("bigMean", "Mean of Big Peak", 150, 140, 160);
    RooRealVar bigSigma("bigSigma", "Sigma of Big Peak", 5, 0.1, 20);
    
    // If you strongly believe bigMean must be ~150 +/- 1:
    bigMean.setRange(149., 151.);
    // narrower bounds
    // or completely fix it:
    //bigMean.setVal(150.);
    //bigMean.setConstant(true);

    RooGaussian gaussBig("gaussBig", "Large Gaussian", x, bigMean, bigSigma);

    // --- The Many Small Peaks
    int nPeaks = 37; // e.g., from 100 -> 172 in 2 ns steps
    double firstMeanVal = 111;
    double spacing = 2;
    RooRealVar smallSigma("smallSigma", "Common Sigma for Small Peaks", 2, 0.1, 10);

    RooArgList smallGaussians;
    for (int i = 0; i < nPeaks; i++)
    {
        double meanVal = firstMeanVal + i * spacing; // 100, 102, 104, ...
        auto *meanVar = new RooRealVar(Form("meanSmall%d", i), "small peak mean", meanVal);
        meanVar->setConstant(true); // fix the small-peak means

        auto *g = new RooGaussian(Form("gaussSmall%d", i),
                                  "small Gauss",
                                  x, *meanVar, smallSigma);
        smallGaussians.add(*g);
    }

    // We want these 37 small Gaussians to share equal fractions among themselves.
    // That means each has fraction 1/37.
    // For 37 peaks, we define (37-1)=36 fraction parameters = 1/37,
    // and the last fraction is automatically 1 - sum(others).
    RooArgList smallFracs;
    for (int i = 0; i < nPeaks - 1; i++)
    {
        double fracVal = 1.0 / double(nPeaks); // ~0.027
        auto *fracVar = new RooRealVar(Form("fracS%d", i),
                                       "fraction for small peak",
                                       fracVal);
        fracVar->setConstant(true);
        smallFracs.add(*fracVar);
    }

    RooAddPdf smallSum("smallSum", "Sum of small peaks",
                       smallGaussians, smallFracs);

    // Now combine { gaussBig, smallSum } in a 2-component RooAddPdf
    // with a single fraction parameter for the big peak
    RooRealVar fracBig("fracBig", "Fraction of Big Peak", 0.5, 0.0, 1.0);
    RooAddPdf model("model", "Big + smallSum",
                    RooArgList(gaussBig, smallSum),
                    RooArgList(fracBig));
    // => total PDF = fracBig * gaussBig + (1 - fracBig) * smallSum

    // --------------------------------------------------------------------------
    // 6. Fit unbinned
    // --------------------------------------------------------------------------
    RooFitResult *fitRes = model.fitTo(data,
                                       RooFit::Save(true),
                                       RooFit::PrintLevel(1),
                                       RooFit::NumCPU(1));
    if (fitRes)
    {
        fitRes->Print("v");
    }
    else
    {
        std::cerr << "Fit failed or was interrupted.\n";
    }

    // --------------------------------------------------------------------------
    // 7. Plot the results
    // --------------------------------------------------------------------------
    TCanvas c("c", "Constrained Summation Fit", 800, 600);
    RooPlot *frame = x.frame();
    data.plotOn(frame, RooFit::Binning(650, 100, 200));    
    model.plotOn(frame); // draw the total PDF
    frame->SetTitle("One Large Peak + 37 Small Peaks (Constrained Summation)");
    frame->Draw();

    c.SaveAs(outputFileName.c_str());
    std::cout << "Saved plot to " << outputFileName << std::endl;

    return 0;
}
