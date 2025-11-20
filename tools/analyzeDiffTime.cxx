#include <iostream>
#include <iomanip>
#include <cstdlib>   // for atoi/atol
#include <vector>
#include <string>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TBox.h>   // for red box around highest-energy cluster
#include <TText.h>  // for cluster ID text
#include <cstdio>   // for Form()

// Detector geometry
static const int ROWS = 36;
static const int COLS = 30;
static const int N_BLOCKS = ROWS * COLS;

// ------------------------------------------------------------------------
// Print a 36×30 table of cluster IDs
void printCalorimeter(const double* block_clusterID) {
    std::cout << "\nCalorimeter Layout (Cluster IDs):\n";
    for (int row = ROWS - 1; row >= 0; --row) {
        for (int col = 0; col < COLS; ++col) {
            int index = row * COLS + col;
            // If clusterID == -1, we print 0
            int cid = (block_clusterID[index] != -1)
                      ? static_cast<int>(block_clusterID[index])
                      : 0;
            std::cout << std::setw(3) << cid;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// ------------------------------------------------------------------------
// Print block-by-block info (DiffTime + ClusT + triggers)
void print_adc_tdc_diff_time(const double* adctdc_difftime,
                             const double* block_clusterID,
                             bool trig1Fired,
                             bool trig6Fired,
                             const double* clusT,
                             double nclustDouble)
{
    int nclustInt = static_cast<int>(nclustDouble);

    // Trigger label
    std::string firedForEvent;
    if (trig1Fired && trig6Fired) {
        firedForEvent = "TRIG1 & TRIG6";
    } else if (trig1Fired) {
        firedForEvent = "TRIG1";
    } else if (trig6Fired) {
        firedForEvent = "TRIG6";
    } else {
        firedForEvent = "None";
    }

    std::cout << "\nBlock-by-block listing (DiffTime & ClusT):\n";

    for(int i = 0; i < N_BLOCKS; i++) {
        double dt = adctdc_difftime[i];
        if (dt >= 1.e10) continue; // skip invalid placeholders

        int cid = static_cast<int>(block_clusterID[i]);
        // Default cluster time is "N/A" if invalid
        std::string clusT_str = "N/A";
        if (cid >= 0 && cid < nclustInt) {
            double thisClusTime = clusT[cid];
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%.2f", thisClusTime);
            clusT_str = buf;
        }

        std::cout 
            << "Block "       << std::setw(4) << i
            << " | ClusterID=" << std::setw(3) << cid
            << " | DiffTime="  << std::setw(7) << std::fixed << std::setprecision(2) << dt
            << " | ClusT="     << clusT_str
            << " | Triggers="  << firedForEvent
            << std::endl;
    }
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <root-file> [startEvent]\n";
        return 1;
    }

    // Root file + optional start event
    const char* input_file = argv[1];
    Long64_t startEvent = 1003;
    if (argc > 2) {
        startEvent = std::atol(argv[2]);
    }
    Long64_t stopEvent = startEvent + 100;  // We'll do 100 events

    // Open file
    TFile* f1 = new TFile(input_file);
    if (!f1 || f1->IsZombie()) {
        std::cerr << "Error: cannot open file " << input_file << "\n";
        return 1;
    }
    TTree* t = (TTree*) f1->Get("T");
    if (!t) {
        std::cerr << "Error: cannot find TTree 'T' in " << input_file << "\n";
        return 1;
    }

    Long64_t nEntries = t->GetEntries();
    if (startEvent < 0 || startEvent >= nEntries) {
        std::cerr << "Error: startEvent=" << startEvent
                  << " out of range [0, " << (nEntries -1) << "]\n";
        return 1;
    }
    if (stopEvent > nEntries) {
        stopEvent = nEntries;
    }

    // We'll store PNGs in a subdir named e.g. blocktimes/events_1003_to_1102
    std::string subDir = Form("blocktimes/events_%lld_to_%lld",
                              startEvent, stopEvent-1);
    gSystem->mkdir(subDir.c_str(), kTRUE);

    // TTree branch variables
    double TRIG6, TRIG1, edtmtdc;
    double hdelta, hcaltot, hcernpe;
    double block_clusterID[N_BLOCKS];
    double goodAdcTdcDiffTime[N_BLOCKS];
    double nclustDouble;
    double clusE[10000];
    double clusT[10000];

    // Set branch addresses
    t->SetBranchAddress("T.hms.npsTRIG6_tdcTimeRaw", &TRIG6);
    t->SetBranchAddress("T.hms.npsTRIG1_tdcTimeRaw", &TRIG1);
    t->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw",    &edtmtdc);

    t->SetBranchAddress("H.gtr.dp",       &hdelta);
    t->SetBranchAddress("H.cal.etotnorm", &hcaltot);
    t->SetBranchAddress("H.cer.npeSum",   &hcernpe);

    t->SetBranchAddress("NPS.cal.fly.block_clusterID", &block_clusterID);
    t->SetBranchAddress("NPS.cal.fly.goodAdcTdcDiffTime", &goodAdcTdcDiffTime);

    t->SetBranchAddress("NPS.cal.nclust", &nclustDouble);
    t->SetBranchAddress("NPS.cal.clusE",  &clusE);
    t->SetBranchAddress("NPS.cal.clusT",  &clusT);

    // We'll define the cluster E cut for the "broad" classification
    // (and also for "narrow" if you want the same E cut).
    const double ENERGY_CUT = 0.6;

    // Process up to 100 events
    for(Long64_t ev = startEvent; ev < stopEvent; ev++){

        t->GetEntry(ev);

        // 1) HPC/EDTM cuts
        if (
            edtmtdc > 0.1
            || hdelta < -12.0 || hdelta > 12.0
            || hcaltot < 0.6
            || hcernpe < 1.0
        ) {
            std::cout << "Event " << ev
                      << " fails HPC/EDTM cuts -> skipping.\n";
            continue;
        }

        bool trig1Fired = (TRIG1>0 && TRIG1<1.e9);
        bool trig6Fired = (TRIG6>0 && TRIG6<1.e9);

        std::cout << "\n============================================\n";
        std::cout << "Event " << ev
                  << " / TRIG1=" << TRIG1
                  << " TRIG6="   << TRIG6
                  << " EDTM="    << edtmtdc
                  << " dp="      << hdelta
                  << " calE="    << hcaltot
                  << " cer="     << hcernpe
                  << "\nTriggers: ";
        if(trig1Fired) std::cout << "TRIG1 ";
        if(trig6Fired) std::cout << "TRIG6 ";
        if(!trig1Fired && !trig6Fired) std::cout << "None";
        std::cout << "\n============================================\n";

        // (A) Print cluster layout
        printCalorimeter(block_clusterID);

        // (B) Print block-level lines
        print_adc_tdc_diff_time(
            goodAdcTdcDiffTime,
            block_clusterID,
            trig1Fired,
            trig6Fired,
            clusT,
            nclustDouble
        );

        int nclust = (int)nclustDouble;
        // We'll mark each cluster if it meets E>=0.6 and time in [130..170] or [150..160]
        std::vector<bool> clusterInBroad(nclust,false);
        std::vector<bool> clusterInNarrow(nclust,false);

        bool haveBroad = false; // track if at least one cluster meets broad
        for(int cID=0; cID<nclust; cID++){
            double eVal= clusE[cID];
            double tVal= clusT[cID];

            // broad => eVal>=0.6, t in [130..170]
            if(eVal>=ENERGY_CUT && tVal>=130. && tVal<=170.){
                clusterInBroad[cID]=true;
                haveBroad = true;
            }
            // narrow => eVal>=0.6, t in [150..160]
            if(eVal>=ENERGY_CUT && tVal>=150. && tVal<=160.){
                clusterInNarrow[cID]=true;
            }
        }

        // (C) If no cluster meets broad (which user specifically requires), skip event
        if(!haveBroad) {
            std::cout << "Event " << ev 
                      << " has no cluster passing (E>="
                      << ENERGY_CUT << " & T in [130..170]) -> skipping.\n";
            continue;
        }

        // (D) Among clusters that pass broad or narrow, find highest-E cluster
        double ehi=-1.0;
        int nhi=-1;
        for(int cID=0; cID<nclust; cID++){
            if(clusterInBroad[cID] || clusterInNarrow[cID]){
                double eVal= clusE[cID];
                if(eVal> ehi){
                    ehi=eVal;
                    nhi=cID;
                }
            }
        }

        // Build 2 hist
        TH2F* h2Broad = new TH2F(
            Form("h2Broad_ev%lld",ev),
            Form("Broad [130..170,E>%.2f], Ev%lld;Col;Row", ENERGY_CUT, ev),
            COLS, 0,(double)COLS,
            ROWS, 0,(double)ROWS
        );
        TH2F* h2Narrow= new TH2F(
            Form("h2Narrow_ev%lld",ev),
            Form("Narrow [150..160,E>%.2f], Ev%lld;Col;Row", ENERGY_CUT, ev),
            COLS, 0,(double)COLS,
            ROWS, 0,(double)ROWS
        );

        // Fill them
        for(int iBlock=0;iBlock<N_BLOCKS;iBlock++){
            int cID=(int)block_clusterID[iBlock];
            if(cID<0||cID>=nclust) continue;
            double dt= goodAdcTdcDiffTime[iBlock];
            if(dt>=1.e10) continue;

            if(clusterInBroad[cID]){
                int row=iBlock/COLS, col=iBlock%COLS;
                h2Broad->SetBinContent(col+1,row+1,dt);
            }
            if(clusterInNarrow[cID]){
                int row=iBlock/COLS, col=iBlock%COLS;
                h2Narrow->SetBinContent(col+1,row+1,dt);
            }
        }

        // Force color ranges
        double broadMin=130., broadMax=170.;
        double narMin=150.,  narMax=160.;

        // Make a single canvas => 1200x720 => each pad 600x720 => ratio=1.2 => squares
        TCanvas* c1= new TCanvas(Form("c1_ev%lld",ev),
                                 "Broad vs Narrow",1200,720);
        gStyle->SetOptStat(0);
        c1->Divide(2,1);

        // (E) Left => broad
        c1->cd(1);
        h2Broad->Draw("COLZ");
        h2Broad->GetZaxis()->SetRangeUser(broadMin,broadMax);

        // Red box if nhi is also in broad
        if(nhi>=0 && clusterInBroad[nhi]){
            for(int iBlock=0;iBlock<N_BLOCKS;iBlock++){
                if((int)block_clusterID[iBlock]==nhi){
                    int row=iBlock/COLS, col=iBlock%COLS;
                    TBox* box=new TBox(col,row,col+1,row+1);
                    box->SetLineColor(kRed);
                    box->SetLineWidth(3);
                    box->SetFillStyle(0);
                    box->Draw("same");
                }
            }
        }
        // cluster ID => only for blocks that pass broad
        {
            TText txt;
            txt.SetTextSize(0.03);
            txt.SetTextColor(kRed);
            txt.SetTextAlign(22);
            txt.SetNDC(false);

            for(int iBlock=0;iBlock<N_BLOCKS;iBlock++){
                int cID=(int)block_clusterID[iBlock];
                if(cID>=0 && cID<nclust && clusterInBroad[cID]){
                    int row=iBlock/COLS, col=iBlock%COLS;
                    double xC= col+0.5;
                    double yC= row+0.5;
                    txt.DrawText(xC,yC,Form("%d",cID));
                }
            }
        }

        // (F) Right => narrow
        c1->cd(2);
        h2Narrow->Draw("COLZ");
        h2Narrow->GetZaxis()->SetRangeUser(narMin,narMax);

        if(nhi>=0 && clusterInNarrow[nhi]){
            for(int iBlock=0;iBlock<N_BLOCKS;iBlock++){
                if((int)block_clusterID[iBlock]==nhi){
                    int row=iBlock/COLS, col=iBlock%COLS;
                    TBox* box=new TBox(col,row,col+1,row+1);
                    box->SetLineColor(kRed);
                    box->SetLineWidth(3);
                    box->SetFillStyle(0);
                    box->Draw("same");
                }
            }
        }
        // cluster ID => only for blocks that pass narrow
        {
            TText txt;
            txt.SetTextSize(0.03);
            txt.SetTextColor(kRed);
            txt.SetTextAlign(22);
            txt.SetNDC(false);

            for(int iBlock=0;iBlock<N_BLOCKS;iBlock++){
                int cID=(int)block_clusterID[iBlock];
                if(cID>=0 && cID<nclust && clusterInNarrow[cID]){
                    int row=iBlock/COLS, col=iBlock%COLS;
                    double xC= col+0.5;
                    double yC= row+0.5;
                    txt.DrawText(xC,yC,Form("%d",cID));
                }
            }
        }

        // Save
        c1->SaveAs(Form("%s/broad_and_narrow_event%lld.png",
                        subDir.c_str(), ev));

        // cleanup
        delete h2Broad;
        delete h2Narrow;
        delete c1;
    } // end event loop

    f1->Close();
    delete f1;
    return 0;
}
