#include <iostream>
#include <iomanip>
#include <cstdlib>  // for atoi/atol
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TText.h> // for drawing text

// Detector dimensions
const int ROWS = 36;  // row index: 0..35
const int COLS = 30;  // col index: 0..29
const int N_BLOCKS = ROWS * COLS; // 1080 total blocks

// ------------------------------------------------------------------------
// Print the cluster layout (text view)
void printCalorimeter(const double* block_clusterID) {
    // Iterate rows from top (35) to bottom (0)
    for (int row = ROWS - 1; row >= 0; --row) {
        for (int col = 0; col < COLS; ++col) {
            int index = row * COLS + col;
            // If clusterID == -1, show 0
            int clusterID = (block_clusterID[index] != -1)
                            ? static_cast<int>(block_clusterID[index])
                            : 0;
            std::cout << std::setw(2) << clusterID << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// ------------------------------------------------------------------------
// Print valid diff-time info (text), now including which triggers fired
// for each block in this event.
void print_adc_tdc_diff_time(const double* adctdc_difftime,
                             const double* clusterID,
                             bool trig1Fired,
                             bool trig6Fired)
{
    std::cout << "Block-by-block listing:\n";

    // We'll construct a small label showing TRIG1/TRIG6/Both/None
    std::string firedForEvent;
    if (trig1Fired && trig6Fired)      firedForEvent = "TRIG1 & TRIG6";
    else if (trig1Fired)              firedForEvent = "TRIG1";
    else if (trig6Fired)              firedForEvent = "TRIG6";
    else                              firedForEvent = "None";

    for(int i = 0; i < N_BLOCKS; i++) {
        // Skip "invalid" times, e.g., > 1.e10
        if(adctdc_difftime[i] < 1.e10) {
            std::cout << "  Block " << i
                      << " | ClusterID = " << clusterID[i]
                      << " | DiffTime = " << adctdc_difftime[i]
                      << " | Triggers = " << firedForEvent
                      << std::endl;
        }
    }
}

// ------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Check usage
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <root-file> [event-number]"
                  << std::endl;
        return 1;
    }

    // 1) Root file name from argv[1]
    const char* input_file = argv[1];

    // 2) Optional single event number, default = 1003
    Long64_t eventNum = 1003;
    if (argc > 2) {
        eventNum = std::atol(argv[2]);
    }

    // Open file and retrieve TTree
    TFile* f1 = new TFile(input_file);
    if (!f1 || f1->IsZombie()) {
        std::cerr << "Error: cannot open file " << input_file << std::endl;
        return 1;
    }
    TTree* myTree = (TTree*) f1->Get("T");
    if (!myTree) {
        std::cerr << "Error: cannot find TTree 'T' in file " 
                  << input_file << std::endl;
        return 1;
    }

    // Make directory for output PNGs
    gSystem->mkdir("blocktimes", kTRUE);

    // Branch addresses
    double TRIG6 = -1000;
    double TRIG1 = -1000;
    double edtmtdc = -1000;  // We'll use this for the EDTM cut
    double block_clusterID[N_BLOCKS];
    double adctdc_difftime[N_BLOCKS];

    // Link branches
    myTree->SetBranchAddress("T.hms.npsTRIG6_tdcTimeRaw", &TRIG6);
    myTree->SetBranchAddress("T.hms.npsTRIG1_tdcTimeRaw", &TRIG1);
    myTree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw",    &edtmtdc);
    myTree->SetBranchAddress("NPS.cal.fly.block_clusterID",&block_clusterID);
    myTree->SetBranchAddress("NPS.cal.fly.goodAdcTdcDiffTime",&adctdc_difftime);

    // Check event range
    Long64_t nEntries = myTree->GetEntries();
    if (eventNum < 0 || eventNum >= nEntries) {
        std::cerr << "Error: requested event " << eventNum
                  << " out of range [0, " << (nEntries - 1) << "]" 
                  << std::endl;
        return 1;
    }

    // Retrieve this single event
    myTree->GetEntry(eventNum);

    // --- Apply EDTM cut: skip if edtmtdc > 0.1 ---
    if (edtmtdc > 0.1) {
        std::cout << "Event " << eventNum
                  << " fails EDTM cut (edtmtdc=" << edtmtdc 
                  << " > 0.1). Exiting.\n";
        f1->Close();
        delete f1;
        return 0;
    }

    // Determine which triggers fired
    // (Example logic: TRIGx is "valid" if 0 < TRIGx < 1.e9)
    bool trig1Fired = (TRIG1 > 0 && TRIG1 < 1.e9);
    bool trig6Fired = (TRIG6 > 0 && TRIG6 < 1.e9);

    // Print event-level info
    std::cout << "\n========================================\n";
    std::cout << "Event: "  << eventNum
              << "  TRIG1: " << TRIG1
              << "  TRIG6: " << TRIG6
              << "  EDTM: "  << edtmtdc
              << std::endl;

    // Quick line about which triggers fired this event
    std::cout << "Triggers that fired this event: ";
    if (trig1Fired) std::cout << "TRIG1 ";
    if (trig6Fired) std::cout << "TRIG6 ";
    if (!trig1Fired && !trig6Fired) std::cout << "None";
    std::cout << "\n========================================\n";

    // 1) Print text representation of cluster IDs
    printCalorimeter(block_clusterID);

    // 2) Print all valid block times, labeling triggers
    print_adc_tdc_diff_time(adctdc_difftime, block_clusterID,
                            trig1Fired, trig6Fired);

    // Prepare canvas
    int width  = 800;
    int height = (int)(width * double(ROWS) / double(COLS));
    TCanvas* c1 = new TCanvas("c1","ADC-TDC Diff Time",width,height);
    gStyle->SetOptStat(0);

    // Adjust margins to reduce distortion
    c1->SetLeftMargin(0.10);
    c1->SetRightMargin(0.15);
    c1->SetTopMargin(0.10);
    c1->SetBottomMargin(0.15);

    // 3) Make a 2D histogram for this event
    TH2F* h2DiffTime = new TH2F(
        Form("h2DiffTime_event%lld", eventNum),
        Form("ADC-TDC Diff Time (Event %lld); Column; Row", eventNum),
        COLS, 0, (double)COLS,
        ROWS, 0, (double)ROWS
    );

    // Track min/max of valid diffTime to set the color scale
    double minVal = 1.e20;
    double maxVal = -1.e20;

    // Fill 2D histogram & find min/max
    for(int i = 0; i < N_BLOCKS; i++) {
        int row = i / COLS;
        int col = i % COLS;
        double diffTime = adctdc_difftime[i];

        if(diffTime < 1.e10) {
            h2DiffTime->SetBinContent(col+1, row+1, diffTime);
            if(diffTime < minVal) minVal = diffTime;
            if(diffTime > maxVal) maxVal = diffTime;
        }
    }

    // If we found any valid times, set z-axis range:
    if (minVal < maxVal) {
        h2DiffTime->GetZaxis()->SetRangeUser(minVal, maxVal);
    }

    // Draw the histogram
    c1->cd();
    h2DiffTime->Draw("COLZ");

    // Overlay text for blocks except -1
    TText txt;
    txt.SetTextSize(0.03);
    txt.SetTextColor(kRed);
    txt.SetTextAlign(22); // center text
    txt.SetNDC(false);

    for(int i = 0; i < N_BLOCKS; i++) {
        int clusterID_val = (int) block_clusterID[i];
        if (clusterID_val != -1) {
            int row = i / COLS;
            int col = i % COLS;
            double xCenter = col + 0.5;
            double yCenter = row + 0.5;
            txt.DrawText(xCenter, yCenter, Form("%d", clusterID_val));
        }
    }

    // Save the PNG in "blocktimes"
    c1->SaveAs(Form("blocktimes/adc_tdc_diff_time_event%lld.png", eventNum));

    // Cleanup
    delete h2DiffTime;
    delete c1;

    f1->Close();
    delete f1;
    return 0;
}
