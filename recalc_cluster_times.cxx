// recalc_cluster_times_slim.cpp
//
// This program reads a ROOT file containing NPS replay data,
// selects a subset of branches, recalculates the cluster times using block-level timing offsets,
// and writes out a new slimmed ROOT file that contains these branches plus an extra branch
// "NPS.cal.newClusT" which holds the new, energy-weighted cluster times.
//
// Usage:
//   ./recalc_cluster_times_slim <RunNumber> <offset_csv_file> <MaxEvents> <output_root_file>
//
// Use MaxEvents = -1 to process all events.
//
// Example:
//   ./recalc_cluster_times_slim 972 Offsets_4205.csv -1 4196_newClusT.root

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

using namespace std;

// Helper function to load block-level timing offsets from a CSV file.
std::map<int, double> loadOffsets(const char* filename) {
    std::map<int, double> offsets;
    ifstream file(filename);
    if (!file.is_open()){
        cout << "Error opening CSV file " << filename << endl;
        return offsets;
    }
    string line;
    // Skip header
    getline(file, line);
    while(getline(file, line)) {
        if(line.empty()) continue;
        istringstream ss(line);
        string blockStr, offsetStr, quality;
        if(getline(ss, blockStr, ',') && getline(ss, offsetStr, ',') && getline(ss, quality, ',')) {
            try {
                int blockID = stoi(blockStr);
                double offset = stod(offsetStr);
                offsets[blockID] = offset;
            } catch (std::exception &e) {
                cout << "Error parsing line: " << line << " (" << e.what() << ")" << endl;
            }
        }
    }
    file.close();
    return offsets;
}

int main(int argc, char* argv[]){
    if(argc < 5){
        cout << "Usage: " << argv[0] << " <RunNumber> <offset_csv_file> <MaxEvents> <output_root_file>" << endl;
        cout << "Use MaxEvents = -1 to process all events." << endl;
        return 1;
    }
    
    Int_t RunNumber = atoi(argv[1]);
    const char* csvFile = argv[2];
    Int_t MaxEvents = atoi(argv[3]);
    // Prepend output directory to file name.
    TString outFilePath = Form("/volatile/hallc/nps/jpcrafts/ROOTfiles/Pi_0/%s", argv[4]);
    
    // Load block-level offsets.
    std::map<int, double> blockOffset = loadOffsets(csvFile);
    if(blockOffset.empty()){
        cout << "No offsets loaded. Exiting." << endl;
        return 1;
    }
    
    // Open the input ROOT file.
    TString inFileName = Form("/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%d_1_-1.root", RunNumber);
    TFile *fin = TFile::Open(inFileName, "READ");
    if(!fin || fin->IsZombie()){
        cout << "Error opening file " << inFileName << endl;
        return 1;
    }
    
    // Get the TTree (assumed to be named "T").
    TTree *tree = (TTree*) fin->Get("T");
    if(!tree){
        cout << "Error: TTree 'T' not found in file " << inFileName << endl;
        fin->Close();
        return 1;
    }
    
    // Enable only the desired branches.
    tree->SetBranchStatus("*", 0);
    // HMS-level branches.
    tree->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    tree->SetBranchStatus("H.gtr.dp",               1);
    tree->SetBranchStatus("H.cal.etotnorm",         1);
    tree->SetBranchStatus("H.cer.npeSum",           1);
    tree->SetBranchStatus("H.gtr.th",               1);
    // Cluster-level branches.
    tree->SetBranchStatus("NPS.cal.nclust",         1);
    tree->SetBranchStatus("NPS.cal.clusT",          1);
    tree->SetBranchStatus("NPS.cal.clusE",          1);
    tree->SetBranchStatus("NPS.cal.clusX",          1);
    tree->SetBranchStatus("NPS.cal.clusY",          1);
    // And the fly branches.
    tree->SetBranchStatus("NPS.cal.fly.block_clusterID", 1);
    tree->SetBranchStatus("NPS.cal.fly.goodAdcTdcDiffTime", 1);
    tree->SetBranchStatus("NPS.cal.fly.e",          1);
    
    // We'll create the new branch "NPS.cal.newClusT" in the output file.
    
    // Set branch addresses.
    // nclust and block_clusterID are stored as Double_t.
    Double_t nclust_d;
    tree->SetBranchAddress("NPS.cal.nclust", &nclust_d);
    
    const int nBlocksMax = 1080;
    Double_t *block_clusterID_d = new Double_t[nBlocksMax];
    tree->SetBranchAddress("NPS.cal.fly.block_clusterID", block_clusterID_d);
    
    Double_t *block_t = new Double_t[nBlocksMax];
    tree->SetBranchAddress("NPS.cal.fly.goodAdcTdcDiffTime", block_t);
    
    Double_t *block_e = new Double_t[nBlocksMax];
    tree->SetBranchAddress("NPS.cal.fly.e", block_e);
    
    // Raw cluster time.
    const int maxClusters = 100;
    Double_t clusT_raw[maxClusters];
    tree->SetBranchAddress("NPS.cal.clusT", clusT_raw);
    
    // Also set addresses for clusE, clusX, clusY so they are copied.
    Double_t clusE_raw[maxClusters];
    tree->SetBranchAddress("NPS.cal.clusE", clusE_raw);
    Double_t clusX_raw[maxClusters];
    tree->SetBranchAddress("NPS.cal.clusX", clusX_raw);
    Double_t clusY_raw[maxClusters];
    tree->SetBranchAddress("NPS.cal.clusY", clusY_raw);
    
    // Set branch addresses for HMS cuts.
    Double_t edtmtdc;
    Double_t hdelta;
    Double_t hcaltot;
    Double_t hcernpe;
    Double_t gtrth;
    tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
    tree->SetBranchAddress("H.gtr.dp", &hdelta);
    tree->SetBranchAddress("H.cal.etotnorm", &hcaltot);
    tree->SetBranchAddress("H.cer.npeSum", &hcernpe);
    tree->SetBranchAddress("H.gtr.th", &gtrth);
    
    Long64_t nEntries = tree->GetEntries();
    cout << "Total entries in tree: " << nEntries << endl;
    
    // If MaxEvents is -1, process all events.
    if(MaxEvents == -1) {
        MaxEvents = nEntries;
    }
    
    // Create output ROOT file.
    TFile *fout = TFile::Open(outFilePath.Data(), "RECREATE");
    // Clone the tree structure (only selected branches) with zero entries.
    TTree *newtree = tree->CloneTree(0);
    
    // Create new branch for recalculated cluster times.
    std::vector<double> newClusT;
    TBranch *bNewClusT = newtree->Branch("NPS.cal.newClusT", &newClusT);
    
    // Loop over events.
    for(Long64_t ievt = 0; ievt < nEntries && ievt < MaxEvents; ievt++){
        tree->GetEntry(ievt);
        
        int nclust = static_cast<int>(nclust_d);
        
        newClusT.clear();
        newClusT.resize(nclust, 0.0);
        vector<double> sumEnergy(nclust, 0.0);
        
        // Loop over blocks.
        for (int ib = 0; ib < nBlocksMax; ib++){
            int cid = static_cast<int>(block_clusterID_d[ib]);
            if(cid < 0 || cid >= nclust) continue;
            double offset = 0.0;
            if(blockOffset.find(ib) != blockOffset.end())
                offset = blockOffset[ib];
            double correctedTime = block_t[ib] + offset;
            newClusT[cid] += block_e[ib] * correctedTime;
            sumEnergy[cid] += block_e[ib];
        }
        
        for (int i = 0; i < nclust; i++){
            if(sumEnergy[i] > 0)
                newClusT[i] /= sumEnergy[i];
            else
                newClusT[i] = 0;
        }
        
        // Optionally, print first 10 events.
        if(ievt < 10) {
            cout << "Event " << ievt << ":\n";
            cout << "  Raw cluster times: ";
            for (int i = 0; i < nclust; i++){
                cout << clusT_raw[i] << " ";
            }
            cout << "\n  New cluster times: ";
            for (int i = 0; i < nclust; i++){
                cout << newClusT[i] << " ";
            }
            cout << "\n";
        }
        
        newtree->Fill();
        
        if((ievt+1) % 10000 == 0)
            cout << "Processed " << (ievt+1) << " events." << endl;
    }
    
    newtree->Write("", TObject::kOverwrite);
    fout->Close();
    fin->Close();
    
    delete [] block_clusterID_d;
    delete [] block_t;
    delete [] block_e;
    
    cout << "Finished processing. Output saved to " << outFilePath.Data() << endl;
    return 0;
}
