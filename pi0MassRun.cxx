#include <iostream>
#include <iomanip>
#include <cstdlib> // for atoi/atol
#include <vector>
#include <cmath>  // for sqrt, sin, fabs
#include <cstdio> // for Form()

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TEntryList.h>
#include <TStyle.h>
#include <TF1.h>
#include <TString.h>
#include <ROOT/RLogger.hxx>    // optional, for debug logs
#include <ROOT/RDataFrame.hxx> // for more advanced parallel, optional

// for EnableImplicitMT
#include <ROOT/TSeq.hxx>
// #include <ROOT/EnableImplicitMT.hxx>

// 1) HPC/EDTM TCut:
// We want T.hms.hEDTM_tdcTimeRaw <0.1 && H.gtr.dp in [-12..12]
// plus H.cal.etotnorm>0.6, H.cer.npeSum>1.0
TString HPCcut = "T.hms.hEDTM_tdcTimeRaw<0.1 && H.gtr.dp>-12 && H.gtr.dp<12 && H.cal.etotnorm>0.6 && H.cer.npeSum>1.0";

// 2) Distance from target to NPS in cm (6.07 m)
static const double DNPS = 607.0;

// 3) Check if a cluster passes the "broad" time + energy cut (clusE >=0.6, clusT in [130..170])
static bool isGoodBroad(double e, double t)
{
    if (e < 0.6)
        return false;
    if (t < 130.0 || t > 170.0)
        return false;
    return true;
}

int main(int argc, char *argv[])
{
    // Step (4) Enable implicit multithreading (you can choose # of threads, or let ROOT decide)
    // ROOT::EnableImplicitMT();  // tries to auto-detect number of cores

    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <run-number>\n";
        return 1;
    }

    int nrun = std::atoi(argv[1]);

    // Build the input ROOT file path, e.g.: nps_hms_skim_<run>_1_-1.root
    std::string inputFileName = Form("/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%d_1_-1.root", nrun);
    std::cout << "Input file: " << inputFileName << std::endl;

    // Hardcoded PNG output directory
    std::string pngDir = "/volatile/hallc/nps/jpcrafts/Plots";
    gSystem->mkdir(pngDir.c_str(), true);

    // Output file name includes run number
    std::string outName = Form("%s/pi0_mass_run%d.png", pngDir.c_str(), nrun);
    std::cout << "Will save final histogram to: " << outName << std::endl;

    // Open the ROOT file
    TFile *f = TFile::Open(inputFileName.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "Error: Cannot open file " << inputFileName << std::endl;
        return 1;
    }
    TTree *t = (TTree *)f->Get("T");
    if (!t)
    {
        std::cerr << "Error: Cannot find TTree 'T' in " << inputFileName << std::endl;
        f->Close();
        return 1;
    }

    Long64_t nEntries = t->GetEntries();
    std::cout << "Opened TTree with " << nEntries << " entries." << std::endl;

    // --------------------------------------------------------------------
    // 1) Disable all branches, then enable only the ones we need:
    // HPC/EDTM event-level, plus cluster arrays:
    // --------------------------------------------------------------------
    t->SetBranchStatus("*", 0); // disable all
    // HPC-level
    t->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    t->SetBranchStatus("H.gtr.dp", 1);
    t->SetBranchStatus("H.cal.etotnorm", 1);
    t->SetBranchStatus("H.cer.npeSum", 1);
    // cluster-level
    t->SetBranchStatus("NPS.cal.nclust", 1);
    t->SetBranchStatus("NPS.cal.clusE", 1);
    t->SetBranchStatus("NPS.cal.clusT", 1);
    t->SetBranchStatus("NPS.cal.clusX", 1);
    t->SetBranchStatus("NPS.cal.clusY", 1);

    // now we set branch addresses only for these enabled ones
    double edtmtdc = -1, hdelta = -999, hcaltot = -999, hcernpe = -999;

    double nclustDouble = 0.0;
    double clusE[10000];
    double clusT[10000];
    double clusX[10000];
    double clusY[10000];

    t->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
    t->SetBranchAddress("H.gtr.dp", &hdelta);
    t->SetBranchAddress("H.cal.etotnorm", &hcaltot);
    t->SetBranchAddress("H.cer.npeSum", &hcernpe);

    t->SetBranchAddress("NPS.cal.nclust", &nclustDouble);
    t->SetBranchAddress("NPS.cal.clusE", &clusE);
    t->SetBranchAddress("NPS.cal.clusT", &clusT);
    t->SetBranchAddress("NPS.cal.clusX", &clusX);
    t->SetBranchAddress("NPS.cal.clusY", &clusY);

    // --------------------------------------------------------------------
    // 2) Create a TEntryList using HPC/EDTM TCut
    //    Then the main loop only iterates over passing events.
    // --------------------------------------------------------------------

    // comment or remove this line if you prefer the default behavior
    // t->SetAutoFlush(-1);

    // Make a "TString" HPCcut defined above
    // HPCcut: "T.hms.hEDTM_tdcTimeRaw<0.1 && H.gtr.dp>-12 && H.gtr.dp<12 && H.cal.etotnorm>0.6 && H.cer.npeSum>1.0"
    std::cout << "Defining HPCcut = " << HPCcut << std::endl;
    t->Draw(">>elist", HPCcut, "entrylist"); // create TEntryList for events that pass HPCcut
    TEntryList *elist = (TEntryList *)gDirectory->Get("elist");
    if (!elist)
    {
        std::cerr << "Error: TEntryList not created. Possibly no events pass the HPC cut!\n";
        f->Close();
        return 1;
    }
    t->SetEntryList(elist);
    Long64_t nSelected = elist->GetN();
    std::cout << "TEntryList selected " << nSelected << " events passing HPC/EDTM.\n";

    // create the pi0 mass histogram
    TH1F *hMgg = new TH1F("hMgg", "Two-Cluster Mass;M_{#gamma#gamma} (GeV);Counts",
                          200, 0.0, 0.3);

    // main loop (only over selected events in HPC/EDTM cut)
    for (Long64_t idx = 0; idx < nSelected; idx++)
    {
        Long64_t entryIndex = elist->GetEntry(idx); // actual TTree entry
        t->GetEntry(entryIndex);

        // cluster-level
        int nclust = (int)nclustDouble;
        if (nclust < 2)
            continue;

        // gather broad clusters
        std::vector<int> goodIdx;
        goodIdx.reserve(nclust);
        for (int cID = 0; cID < nclust; cID++)
        {
            if (isGoodBroad(clusE[cID], clusT[cID]))
            {
                goodIdx.push_back(cID);
            }
        }
        if ((int)goodIdx.size() < 2)
            continue;

        // form all pairs
        for (size_t i1 = 0; i1 < goodIdx.size(); i1++)
        {
            int cID1 = goodIdx[i1];
            double E1 = clusE[cID1];
            double x1 = clusX[cID1];
            double y1 = clusY[cID1];
            double t1 = clusT[cID1];

            for (size_t i2 = i1 + 1; i2 < goodIdx.size(); i2++)
            {
                int cID2 = goodIdx[i2];
                double E2 = clusE[cID2];
                double x2 = clusX[cID2];
                double y2 = clusY[cID2];
                double t2 = clusT[cID2];

                // time difference ±3ns
                double tdiff = std::fabs(t1 - t2);
                if (tdiff > 3.0)
                    continue;

                // angle
                double dx = x1 - x2;
                double dy = y1 - y2;
                double dist = std::sqrt(dx * dx + dy * dy);

                double theta = dist / DNPS;
                double s2 = std::sin(0.5 * theta);
                double Mgg = std::sqrt(4. * E1 * E2 * s2 * s2);

                hMgg->Fill(Mgg);
            }
        }
    } // end HPC/EDTM-selected events loop

    std::cout << "Done looping over HPC/EDTM-selected events.\n"
              << "Filling final mass histogram, now drawing & fitting.\n";

    // Optional: enable stats & fit box
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    // define a TCanvas, draw the histogram
    TCanvas *c1 = new TCanvas("c1", "Pi0 Mass", 800, 600);
    hMgg->Draw();

    // define a gaussian function in [0.12..0.15]
    TF1 *fitGauss = new TF1("fitGauss", "gaus", 0.119, 0.145);
    fitGauss->SetLineColor(kRed);
    // fitGauss->SetParameters(300.,0.135,0.01); // amplitude, mean, sigma guess

    // Fit
    hMgg->Fit(fitGauss, "R");

    // Save
    std::string outPng = Form("%s/pi0_mass_run%d.png", pngDir.c_str(), nrun);
    c1->SaveAs(outPng.c_str());

    std::cout << "Saved mass histogram with fit to: " << outPng << std::endl;

    // cleanup
    delete fitGauss;
    delete c1;
    delete hMgg;
    f->Close();
    delete f;

    return 0;
}
