#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>

const int ROWS = 36;  // 0-30
const int COLS = 30;  // 0-35

void printCalorimeter(double* block_clusterID) {
    // Iterate over rows from top to bottom to place (0,0) at the bottom-right
    for (int row = ROWS - 1; row >= 0; --row) {
        for (int col = 0; col < COLS; ++col) {
            int index = row * COLS + col;
            int clusterID = (block_clusterID[index] != 0 && block_clusterID[index] != -1) 
                            ? static_cast<int>(block_clusterID[index]) 
                            : 0;
            std::cout << std::setw(2) << clusterID << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print_adc_tdc_diff_time(double* adctdc_difftime, double* clusterID) {
  int block_num = 0;
  for(int i = 0; i < 1080; i++) {
    if(adctdc_difftime[i] < 1.e10) {
      std::cout << i << " " << clusterID[i] << " " << adctdc_difftime[i] << std::endl;
    }
  }
  return;
}

int main(int argc, char* argv[]) {
  const char* input_file = argv[1];
  TFile* f1 = new TFile(input_file);
  TTree* t = (TTree*) f1->Get("T");
  double TRIG6=-1000;t->SetBranchAddress("T.hms.npsTRIG6_tdcTimeRaw",&TRIG6);
  double TRIG1=-1000;t->SetBranchAddress("T.hms.npsTRIG1_tdcTimeRaw",&TRIG1);
  double block_clusterID[2000]; t->SetBranchAddress("NPS.cal.fly.block_clusterID",&block_clusterID);
  double adctdc_difftime[2000]; t->SetBranchAddress("NPS.cal.fly.goodAdcTdcDiffTime",&adctdc_difftime);

  // Inside your event loop:
  for (int entry = 1003; entry < 1004; ++entry) {
    t->GetEntry(entry);
    //if (TRIG6 > 100) continue;
    std::cout << "Entry: " << entry << " TRIG1: " << TRIG1 << " TRIG6: " << TRIG6 << '\n';
    printCalorimeter(block_clusterID);
    std::cout << std::endl;
    print_adc_tdc_diff_time(adctdc_difftime, block_clusterID);
    std::cout << std::endl;
  }
  
  return 0;
}
