#include <iostream>
#include <iomanip>
#include <cstdlib> // for atoi/atol
#include <vector>
#include <cmath>   // for sqrt, sin, fabs
#include <cstdio>  // for Form()

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

// HPC/EDTM TCut
TString HPCcut =
  "T.hms.hEDTM_tdcTimeRaw<0.1 && "
  "H.gtr.dp>-8.5 && H.gtr.dp<8.5 && "
  "H.cal.etotnorm>0.6 && "
  "H.cer.npeSum>1.0";
  //HMS Colimaotor cut xptar and yptr at abs <.1

// Distance from target to NPS in cm (6.07 m)
static const double DNPS = 607.0;

// Check if a cluster passes the "broad" time + energy cut
static bool isGoodBroad(double e, double t)
{
    if (e < 0.6) return false;
    if (t < 130.0 || t > 170.0) return false;
    return true;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run-number>\n";
        return 1;
    }

    int nrun = std::atoi(argv[1]);

    // Input file path
    std::string inputFileName =
      Form("/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%d_1_-1.root", nrun);
    std::cout << "Input file: " << inputFileName << std::endl;

    // Directory for PNGs
    std::string pngDir = "/volatile/hallc/nps/jpcrafts/Plots";
    gSystem->mkdir(pngDir.c_str(), true);

    // Directory for ROOT files
    std::string rootDir = "/volatile/hallc/nps/jpcrafts/ROOTfiles/Pi_0";
    gSystem->mkdir(rootDir.c_str(), true);

    // Output file names
    std::string outPng  = Form("%s/pi0_mass_run%d.png",  pngDir.c_str(),  nrun);
    std::string outROOT = Form("%s/pi0_mass_run%d.root", rootDir.c_str(), nrun);

    std::cout << "Will save final histogram PNG to : " << outPng  << std::endl;
    std::cout << "Will also write histogram to ROOT file: " << outROOT << std::endl;

    // Open the input ROOT file
    TFile *f = TFile::Open(inputFileName.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open file " << inputFileName << std::endl;
        return 1;
    }
    TTree *t = (TTree *)f->Get("T");
    if (!t) {
        std::cerr << "Error: Cannot find TTree 'T' in " << inputFileName << std::endl;
        f->Close();
        return 1;
    }
    Long64_t nEntries = t->GetEntries();
    std::cout << "Opened TTree with " << nEntries << " entries." << std::endl;

    // Disable all branches, enable only needed
    t->SetBranchStatus("*", 0);
    // HPC-level
    t->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    t->SetBranchStatus("H.gtr.dp",               1);
    t->SetBranchStatus("H.cal.etotnorm",         1);
    t->SetBranchStatus("H.cer.npeSum",           1);
    // cluster-level
    t->SetBranchStatus("NPS.cal.nclust",         1);
    t->SetBranchStatus("NPS.cal.clusE",          1);
    t->SetBranchStatus("NPS.cal.clusT",          1);
    t->SetBranchStatus("NPS.cal.clusX",          1);
    t->SetBranchStatus("NPS.cal.clusY",          1);

    // set addresses
    double edtmtdc = -1, hdelta = -999, hcaltot = -999, hcernpe = -999;
    double nclustDouble = 0.0;
    double clusE[10000];
    double clusT[10000];
    double clusX[10000];
    double clusY[10000];

    t->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
    t->SetBranchAddress("H.gtr.dp",              &hdelta);
    t->SetBranchAddress("H.cal.etotnorm",        &hcaltot);
    t->SetBranchAddress("H.cer.npeSum",          &hcernpe);

    t->SetBranchAddress("NPS.cal.nclust", &nclustDouble);
    t->SetBranchAddress("NPS.cal.clusE",  &clusE);
    t->SetBranchAddress("NPS.cal.clusT",  &clusT);
    t->SetBranchAddress("NPS.cal.clusX",  &clusX);
    t->SetBranchAddress("NPS.cal.clusY",  &clusY);

    // Create a TEntryList using HPC/EDTM TCut
    std::cout << "Defining HPCcut = " << HPCcut << std::endl;
    t->Draw(">>elist", HPCcut, "entrylist");
    TEntryList *elist = (TEntryList *)gDirectory->Get("elist");
    if (!elist) {
        std::cerr << "Error: TEntryList not created. Possibly no events pass HPC cut!\n";
        f->Close();
        return 1;
    }
    t->SetEntryList(elist);
    Long64_t nSelected = elist->GetN();
    std::cout << "TEntryList selected " << nSelected << " events passing HPC/EDTM.\n";

    // create the pi0 mass histogram
    TH1F *hMgg = new TH1F("hMgg", "Two-Cluster Mass;M_{#gamma#gamma} (GeV);Counts",
                          200, 0.0, 0.3);

    // main loop (only HPC/EDTM-passing events)
    for (Long64_t idx = 0; idx < nSelected; idx++) {
        Long64_t entryIndex = elist->GetEntry(idx); // actual TTree entry
        t->GetEntry(entryIndex);

        int nclust = (int)nclustDouble;
        if (nclust < 2) continue;

        // gather broad clusters
        std::vector<int> goodIdx;
        goodIdx.reserve(nclust);
        for (int cID = 0; cID < nclust; cID++) {
            if (isGoodBroad(clusE[cID], clusT[cID])) {
                goodIdx.push_back(cID);
            }
        }
        if ((int)goodIdx.size() < 2) continue;

        // form pairs
        for (size_t i1 = 0; i1 < goodIdx.size(); i1++) {
            int cID1 = goodIdx[i1];
            double E1 = clusE[cID1];
            double x1 = clusX[cID1];
            double y1 = clusY[cID1];
            double t1 = clusT[cID1];

            for (size_t i2 = i1 + 1; i2 < goodIdx.size(); i2++) {
                int cID2 = goodIdx[i2];
                double E2 = clusE[cID2];
                double x2 = clusX[cID2];
                double y2 = clusY[cID2];
                double t2 = clusT[cID2];

                // time difference ±3ns
                double tdiff = std::fabs(t1 - t2);
                if (tdiff > 3.0) continue;

                // angle
                double dx = x1 - x2;
                double dy = y1 - y2;
                double dist = std::sqrt(dx * dx + dy * dy);

                double theta = dist / DNPS;
                double s2 = std::sin(0.5 * theta);
                double Mgg = std::sqrt(4.0 * E1 * E2 * s2 * s2);

                hMgg->Fill(Mgg);
            }
        }
    }

    std::cout << "Done looping over HPC-passing events. Filled final mass histogram.\n";

    // Optionally do a quick Gaussian fit around ~0.135
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    TCanvas *c1 = new TCanvas("c1", "Pi0 Mass", 800, 600);
    hMgg->Draw();

    // define a gaussian around [0.12..0.15]
    TF1 *fitGauss = new TF1("fitGauss","gaus",0.12,0.145);
    fitGauss->SetLineColor(kRed);
    //fitGauss->SetParameters(300.,0.135,0.01);

    hMgg->Fit(fitGauss,"R");

    // Save the PNG
    c1->SaveAs(outPng.c_str());
    std::cout << "Saved mass histogram with fit to: " << outPng << std::endl;

    // Now also create an output ROOT file to store the histogram and fit
    TFile *outFile = new TFile(outROOT.c_str(),"RECREATE");
    if(outFile->IsOpen()) {
        hMgg->Write();           // store the histogram
        fitGauss->Write();       // store the fit function
        c1->Write("massCanvas"); // store the canvas (optional)
        outFile->Close();

        std::cout << "Wrote histogram & fit to file: " << outROOT << std::endl;
        delete outFile;
    } else {
        std::cerr << "Cannot open output ROOT file " << outROOT << std::endl;
    }

    // cleanup
    delete fitGauss;
    delete c1;
    delete hMgg;
    f->Close();
    delete f;

    return 0;
}
