#include "TFile.h"
#include "TVectorD.h"
#include <vector>
#include <iostream>

//   g++ -std=c++17 -O2 ./make_bins_schema_MxCorr_Mgg.C `root-config --cflags --libs`     -o make_bins_schema_MxCorr_Mgg

void make_bins_schema_MxCorr_Mgg(const char* out="bins_MxCorr_Mgg.root") {
  TFile f(out,"RECREATE");
  f.mkdir("BinEdges"); f.cd("BinEdges");

  // EDIT edges as needed:
  std::vector<double> mx = {0.0, 0.5, 1.0, 1.5, 2.0};             // Mx_corr (GeV)
  std::vector<double> mg = {0.10, 0.12, 0.135, 0.15, 0.18, 0.20}; // mgg (GeV)

  TVectorD vmx(mx.size()); for (int i=0;i<vmx.GetNoElements();++i) vmx[i]=mx[i];
  TVectorD vmg(mg.size()); for (int i=0;i<vmg.GetNoElements();++i) vmg[i]=mg[i];

  vmx.Write("MxCorr");
  vmg.Write("Mgg");
  f.Close();
  printf("Wrote %s with BinEdges/{MxCorr,Mgg}\n", out);
}

// -------------------------------------------------------------
//  Minimal CLI driver
// -------------------------------------------------------------
int main(int argc, char** argv)
{
    const char* outfile = "bins_MxCorr_Mgg.root";

    if (argc == 2) {
        outfile = argv[1];
    } else if (argc > 2) {
        std::cerr << "Usage: " << argv[0] << " [output.root]\n";
        return 1;
    }

    make_bins_schema_MxCorr_Mgg(outfile);
    return 0;
}