/*  count_EDTM_TRIG3.cxx
 *  ----------------------------------------------------------
 *  Count how many events in the skimmed run have
 *     – an EDTM trigger
 *     – an hTRIG3 trigger
 *     – both simultaneously
 *
 *  Compile:
 *    g++ -std=c++17 -O2 -Wall count_EDTM_TRIG3.cxx \
 *        $(root-config --cflags --libs) -o count_EDTM_TRIG3
 *
 *  Run with default filename pattern:
 *    ./count_EDTM_TRIG3 2093
 *
 *  Run with a custom full pattern (must contain %d for run):
 *    ./count_EDTM_TRIG3 2093 "/my/path/nps_hms_skim_%d_1_-1.root"
 */

 #include <iostream>
 #include <string>
 #include <cstdio>
 
 #include "TFile.h"
 #include "TTree.h"
 #include "TString.h"
 
 /*------------------------------------------------------------------*/
 int main(int argc, char* argv[])
 {
   if(argc < 2){
     std::cerr << "\nUsage: " << argv[0] << " RUNNUMBER "
               << "[full_pattern_with_%d]\n\n";
     return 1;
   }
 
   const int run = std::stoi(argv[1]);
 
   /* default directory+pattern (same as before) -------------------- */
   const std::string defaultPattern =
       "/lustre24/expphy/cache/hallc/c-nps/analysis/pass1/replays/skim/"
       "nps_hms_skim_%d_1_-1.root";
 
   const std::string pattern = (argc > 2) ? argv[2] : defaultPattern;
 
   char buf[512];
   std::snprintf(buf, sizeof(buf), pattern.c_str(), run);
   TString fname(buf);
 
   TFile *f = TFile::Open(fname, "READ");
   if(!f || f->IsZombie()){
     std::cerr << "Cannot open " << fname << "\n";
     return 1;
   }
 
   TTree *T = static_cast<TTree*>(f->Get("T"));
   if(!T){
     std::cerr << "T tree not found in " << fname << "\n";
     return 1;
   }
 
   /* ---- branches we need ---------------------------------------- */
   double t_EDTM = 0, t_TRIG3 = 0;
   T->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &t_EDTM);
   T->SetBranchAddress("T.hms.hTRIG3_tdcTimeRaw", &t_TRIG3);
 
   /* ---- counters ------------------------------------------------- */
   Long64_t N_EDTM = 0, N_T3 = 0, N_both = 0;
   const Long64_t nEnt = T->GetEntries();
 
   for(Long64_t i = 0; i < nEnt; ++i){
     T->GetEntry(i);
 
     bool hasEDTM = (t_EDTM > 0);
     bool hasT3   = (t_TRIG3 > 0);
 
     if(hasEDTM) ++N_EDTM;
     if(hasT3)   ++N_T3;
     if(hasEDTM && hasT3) ++N_both;
   }
 
   /* ---- report --------------------------------------------------- */
   std::cout << "\nRun " << run << "  (" << nEnt << " events)\n"
             << "  EDTM   triggers : " << N_EDTM << '\n'
             << "  hTRIG3 triggers : " << N_T3   << '\n'
             << "  Overlap (both)  : " << N_both << '\n' << std::endl;
 
   return 0;
 }
 