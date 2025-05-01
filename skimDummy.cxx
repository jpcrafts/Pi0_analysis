#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>  // for gSystem->AccessPathName
#include <iostream>
#include <string>
#include <cstdio>     // for sprintf

// -----------------------------------------------------------------
// This function skims only the branches you need from an input file
// into an output file, printing a status update every 100k events.
//
bool skimSelectedBranches(const std::string &inFile, const std::string &outFile)
{
    // Open input
    TFile* fIn = TFile::Open(inFile.c_str(), "READ");
    if (!fIn || fIn->IsZombie()){
        std::cerr << "[skim] Error opening input file " << inFile << std::endl;
        return false;
    }

    // Get TTree "T"
    TTree* tIn = (TTree*) fIn->Get("T");
    if (!tIn){
        std::cerr << "[skim] Error: TTree 'T' not found in " << inFile << std::endl;
        fIn->Close();
        return false;
    }

    // Disable all branches
    tIn->SetBranchStatus("*", 0);

    // Re-enable only the needed branches (adjust as needed):
    tIn->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    tIn->SetBranchStatus("H.gtr.dp",              1);
    tIn->SetBranchStatus("H.gtr.th",              1);
    tIn->SetBranchStatus("H.gtr.ph",              1);
    tIn->SetBranchStatus("H.gtr.y",               1);
    tIn->SetBranchStatus("H.cal.etotnorm",        1);
    tIn->SetBranchStatus("H.cer.npeSum",          1);

    // If you need cluster info:
    tIn->SetBranchStatus("NPS.cal.clusT",         1);
    tIn->SetBranchStatus("NPS.cal.clusE",         1);
    tIn->SetBranchStatus("NPS.cal.clusX",         1);
    tIn->SetBranchStatus("NPS.cal.clusY",         1);
    tIn->SetBranchStatus("NPS.cal.nclust",        1);

    // Create output file
    TFile* fOut = TFile::Open(outFile.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()){
        std::cerr << "[skim] Error creating output file " << outFile << std::endl;
        fIn->Close();
        return false;
    }

    // Instead of CloneTree(-1), we do CloneTree(0) plus a manual event loop
    // so we can print status updates.
    TTree* tOut = tIn->CloneTree(0, "fast");

    // Now copy event-by-event
    Long64_t nEntries = tIn->GetEntries();
    const Long64_t chunkSize = 100000; // print progress every 100k
    for (Long64_t i = 0; i < nEntries; i++){
        if ((i > 0) && (i % chunkSize == 0)) {
            std::cout << "[skim] Processed " << i 
                      << " / " << nEntries << " entries from " << inFile << "\n";
        }
        tIn->GetEntry(i);
        tOut->Fill();
    }

    // Write & close
    fOut->Write();
    fOut->Close();
    fIn->Close();

    std::cout << "[skim] Done copying " << nEntries
              << " events from " << inFile
              << " => " << outFile << std::endl;
    return true;
}

// -----------------------------------------------------------------
// main program: usage:
//    ./skimDummyRuns <runNum> <segMax>
//
// This will loop over seg = 0..segMax, building input file names:
//    /cache/hallc/c-nps/analysis/pass2/replays/production/
//        nps_hms_coin_<runNum>_<seg>_1_-1.root
//
// and output file names:
//    /volatile/hallc/nps/jpcrafts/ROOTfiles/Pi_0/
//        nps_hms_coin_<runNum>_<seg>_skim.root
//
// If a file doesn't exist or is unreadable, it will skip and continue.
//
int main(int argc, char* argv[])
{
    if (argc < 3){
        std::cerr << "Usage: " << argv[0]
                  << " <runNum> <segMax>\n\n"
                  << "Example:\n  " << argv[0] << " 6834 5\n"
                  << " --> processes seg=0..5 for run 6834\n";
        return 1;
    }

    int runNum  = std::stoi(argv[1]);
    int segMax  = std::stoi(argv[2]);

    for (int seg = 0; seg <= segMax; seg++){
        // Build input file path
        char inName[512];
        sprintf(inName,
            "/cache/hallc/c-nps/analysis/pass2/replays/production/"
            "nps_hms_coin_%d_%d_1_-1.root",
            runNum, seg);

        // Check if file exists
        if (gSystem->AccessPathName(inName, kFileExists)){
            std::cerr << "[skim] File doesn't exist or not accessible: "
                      << inName << "  -> Skipping.\n";
            continue;
        }

        // Build output path
        char outName[512];
        sprintf(outName,
            "/volatile/hallc/nps/jpcrafts/ROOTfiles/Pi_0/"
            "nps_hms_coin_%d_%d_skim.root",
            runNum, seg);

        // Skim
        skimSelectedBranches(inName, outName);
    }

    std::cout << "[skim] Done skimming seg=0.."
              << segMax << " for run " << runNum << std::endl;

    return 0;
}
