// recalc_cluster_times_slim.cpp
//
// This program reads a ROOT file containing NPS replay data,
// selects a subset of branches, recalculates the cluster times using block-level timing offsets,
// and writes out a new slimmed ROOT file that contains these branches plus two extra branches:
//   "NPS.cal.newEWClusT" holds the energy-weighted cluster times,
//   "NPS.cal.newClusT" holds the simple (arithmetic) average cluster times.
//
// HMS-level branches include T.hms.hEDTM_tdcTimeRaw, H.gtr.dp, H.cal.etotnorm,
// H.cer.npeSum, H.gtr.th, and now H.gtr.ph.
// HMS cuts applied are:
//   T.hms.hEDTM_tdcTimeRaw < 0.1,
//   H.gtr.dp between -8.5 and 8.5,
//   H.cal.etotnorm > 0.6,
//   H.cer.npeSum > 1.0,
//   |H.gtr.th| < 0.09,
//   |H.gtr.ph| < 0.09
//
// Usage:
//   ./recalc_cluster_times_slim <RunNumber> <offset_csv_file> <MaxEvents> <output_root_file>
//
// Use MaxEvents = -1 to process all events.
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

// Modified helper function: reads only the first two columns from the CSV
// and checks that the first data line starts with block 0.
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
    bool firstDataLine = true;
    while(getline(file, line)) {
        if(line.empty()) continue;
        istringstream ss(line);
        string blockToken, offsetToken;
        
        if(!getline(ss, blockToken, ',')) continue;
        int blockID;
        try {
            blockID = stoi(blockToken);
        } catch (std::exception &e) {
            cout << "Error parsing blockID in line: " << line << " (" << e.what() << ")" << endl;
            continue;
        }
        // Check that first data line starts with block 0.
        if(firstDataLine) {
            if(blockID != 0) {
                cout << "Error: The first block ID is not 0. Found " << blockID << endl;
                return std::map<int, double>();  // Return empty map.
            }
            firstDataLine = false;
        }
        
        if(!getline(ss, offsetToken, ',')) continue;
        double offsetVal;
        try {
            offsetVal = stod(offsetToken);
        } catch (std::exception &e) {
            cout << "Error parsing offset in line: " << line << " (" << e.what() << ")" << endl;
            continue;
        }
        offsets[blockID] = offsetVal;
        // Ignore any extra columns.
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
    
    // Load block-level offsets from CSV.
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
    // Global Event Number.
    tree->SetBranchStatus("g.evnum", 1);
    // HMS-level branches.
    tree->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    tree->SetBranchStatus("H.gtr.dp",               1);
    tree->SetBranchStatus("H.cal.etotnorm",         1);
    tree->SetBranchStatus("H.cer.npeSum",           1);
    tree->SetBranchStatus("H.gtr.th",               1);
    tree->SetBranchStatus("H.gtr.ph",               1);
    // Cluster-level branches.
    tree->SetBranchStatus("NPS.cal.nclust",         1);
    tree->SetBranchStatus("NPS.cal.clusT",          1);
    tree->SetBranchStatus("NPS.cal.clusE",          1);
    tree->SetBranchStatus("NPS.cal.clusX",          1);
    tree->SetBranchStatus("NPS.cal.clusY",          1);
    // Fly-level branches.
    tree->SetBranchStatus("NPS.cal.fly.block_clusterID", 1);
    tree->SetBranchStatus("NPS.cal.fly.goodAdcTdcDiffTime", 1);
    tree->SetBranchStatus("NPS.cal.fly.e",          1);
    
    // We'll create two new branches in the output file:
    // "NPS.cal.newEWClusT" for energy-weighted cluster times,
    // "NPS.cal.newClusT" for simple average cluster times.
    
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
    Double_t edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph;
    tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
    tree->SetBranchAddress("H.gtr.dp", &hdelta);
    tree->SetBranchAddress("H.cal.etotnorm", &hcaltot);
    tree->SetBranchAddress("H.cer.npeSum", &hcernpe);
    tree->SetBranchAddress("H.gtr.th", &gtrth);
    tree->SetBranchAddress("H.gtr.ph", &gtrph);
    
    Long64_t nEntries = tree->GetEntries();
    cout << "Total entries in tree: " << nEntries << endl;
    
    // If MaxEvents is -1, process all events.
    if(MaxEvents == -1) {
        MaxEvents = nEntries;
    }
    
    // Create output ROOT file.
    TFile *fout = TFile::Open(outFilePath.Data(), "RECREATE");
    // Clone the tree structure (with selected branches) with zero entries.
    TTree *newtree = tree->CloneTree(0);
    
    // Create new branches for recalculated cluster times.
    // Energy-weighted cluster times.
    std::vector<double> newEWClusT;
    TBranch *bNewEWClusT = newtree->Branch("NPS.cal.newEWClusT", &newEWClusT);
    // Simple average cluster times.
    std::vector<double> newClusT;
    TBranch *bNewClusT = newtree->Branch("NPS.cal.newClusT", &newClusT);
    
    // Loop over events.
    for(Long64_t ievt = 0; ievt < nEntries && ievt < MaxEvents; ievt++){
        tree->GetEntry(ievt);
        
        int nclust = static_cast<int>(nclust_d);
        
        newEWClusT.clear();
        newEWClusT.resize(nclust, 0.0);
        newClusT.clear();
        newClusT.resize(nclust, 0.0);
        
        // Vectors for summing values.
        vector<double> sumEnergy(nclust, 0.0);
        vector<int> countBlocks(nclust, 0);
        
        // Loop over blocks.
        for (int ib = 0; ib < nBlocksMax; ib++){
            int cid = static_cast<int>(block_clusterID_d[ib]);
            if(cid < 0 || cid >= nclust) continue;
            double offset = 0.0;
            if(blockOffset.find(ib) != blockOffset.end())
                offset = blockOffset[ib];
            double correctedTime = block_t[ib] + offset;
            // Energy-weighted accumulation.
            newEWClusT[cid] += block_e[ib] * correctedTime;
            sumEnergy[cid] += block_e[ib];
            // Simple average accumulation.
            newClusT[cid] += correctedTime;
            countBlocks[cid]++;
        }
        
        // Finalize energy-weighted average.
        for (int i = 0; i < nclust; i++){
            if(sumEnergy[i] > 0)
                newEWClusT[i] /= sumEnergy[i];
            else
                newEWClusT[i] = 0;
        }
        // Finalize simple average.
        for (int i = 0; i < nclust; i++){
            if(countBlocks[i] > 0)
                newClusT[i] /= countBlocks[i];
            else
                newClusT[i] = 0;
        }
        
        // Optionally, print first 10 events for debugging.
        if(ievt < 10) {
            cout << "Event " << ievt << ":\n";
            cout << "  Raw cluster times: ";
            for (int i = 0; i < nclust; i++){
                cout << clusT_raw[i] << " ";
            }
            cout << "\n  Energy-weighted new cluster times: ";
            for (int i = 0; i < nclust; i++){
                cout << newEWClusT[i] << " ";
            }
            cout << "\n  Simple average new cluster times: ";
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
