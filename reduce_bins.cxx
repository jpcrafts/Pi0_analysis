#include <vector>
#include <string>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TKey.h"

// g++ -std=c++17 -O0 reduce_bins.C -o reduce_bins `root-config --cflags --libs` -lTMVA -lRooFitCore -lRooFit

#include <string>

/**
 * @brief Merge per-bin histograms from multiple ROOT files.
 *
 * This function:
 *   - Verifies a consistent schema/config hash across all inputs.
 *   - Reads the `bins` TTree from each file.
 *   - Sums weight-related branches per bin ID.
 *   - Copies the `BinEdges/` directory from the first file.
 *   - Writes a merged output ROOT file with:
 *        - bins_merged TTree
 *        - bin_schema_id (TParameter<std::string>)
 *        - config_hash   (TParameter<std::string>)
 *
 * @param outFile  Path to output ROOT file to create.
 * @param inputs   List of input ROOT file paths.
 */
void reduce_bins(const char* outFile,
                 std::vector<std::string> inputs)
{
  if (inputs.empty()) { printf("[reduce] no inputs\n"); return; }

  TFile f0(inputs[0].c_str(),"READ");
  if (!f0.IsOpen()) { printf("[reduce] cannot open %s\n", inputs[0].c_str()); return; }
  auto p_schema = (TParameter<std::string>*)f0.Get("bin_schema_id");
  auto p_config = (TParameter<std::string>*)f0.Get("config_hash");
  if (!p_schema || !p_config) { printf("[reduce] missing schema/config in %s\n", inputs[0].c_str()); return; }
  std::string schema_id = p_schema->GetVal();
  std::string config_id = p_config->GetVal();

  std::map<int,double> sumW, sumW2, sumW_data, sumW_up, sumW_dn, sumW_acc;

  auto add_file = [&](const std::string& path){
    TFile fi(path.c_str(),"READ");
    if (!fi.IsOpen()) { printf("[reduce] skip (open fail): %s\n", path.c_str()); return; }
    auto s = (TParameter<std::string>*)fi.Get("bin_schema_id");
    auto c = (TParameter<std::string>*)fi.Get("config_hash");
    if (!s || !c) { printf("[reduce] skip (no schema/config): %s\n", path.c_str()); return; }
    if (s->GetVal()!=schema_id || c->GetVal()!=config_id) {
      printf("[reduce] ABORT: schema/config mismatch in %s\n", path.c_str());
      gSystem->Exit(1);
    }
    TTree* tb = (TTree*)fi.Get("bins");
    if (!tb) { printf("[reduce] skip (no bins tree): %s\n", path.c_str()); return; }
    int bin_id=0; double sw=0, sw2=0, d=0,u=0,n=0,a=0;
    tb->SetBranchAddress("bin_id",&bin_id);
    tb->SetBranchAddress("SumW",&sw);
    tb->SetBranchAddress("SumW2",&sw2);
    tb->SetBranchAddress("SumW_data",&d);
    tb->SetBranchAddress("SumW_dummyUP",&u);
    tb->SetBranchAddress("SumW_dummyDN",&n);
    tb->SetBranchAddress("SumW_acc",&a);
    Long64_t N = tb->GetEntries();
    for (Long64_t i=0;i<N;i++) {
      tb->GetEntry(i);
      sumW[bin_id]     += sw;
      sumW2[bin_id]    += sw2;
      sumW_data[bin_id]+= d;
      sumW_up[bin_id]  += u;
      sumW_dn[bin_id]  += n;
      sumW_acc[bin_id] += a;
    }
  };

  for (auto& s : inputs) add_file(s);

  TFile fo(outFile,"RECREATE");
  TTree t("bins_merged","merged per-bin yields");
  int bin_id=0; double sw=0, sw2=0, d=0,u=0,n=0,a=0;
  t.Branch("bin_id",&bin_id,"bin_id/I");
  t.Branch("SumW",&sw,"SumW/D");
  t.Branch("SumW2",&sw2,"SumW2/D");
  t.Branch("SumW_data",&d,"SumW_data/D");
  t.Branch("SumW_dummyUP",&u,"SumW_dummyUP/D");
  t.Branch("SumW_dummyDN",&n,"SumW_dummyDN/D");
  t.Branch("SumW_acc",&a,"SumW_acc/D");

  // copy BinEdges/ from first file
  fo.mkdir("BinEdges"); fo.cd("BinEdges");
  if (auto* dbe = (TDirectory*)f0.Get("BinEdges")) {
    TIter next(dbe->GetListOfKeys()); TKey* key;
    while ((key = (TKey*)next())) {
      TObject* obj = dbe->Get(key->GetName());
      if (obj) obj->Write(key->GetName());
    }
  }
  fo.cd();

  // Change TParameter to TNamed for ROOT >= 6.36 workaround
  TNamed("bin_schema_id", schema_id.c_str()).Write();
  TNamed("config_hash",  config_id.c_str()).Write();

  for (auto& kv : sumW) {
    bin_id = kv.first;
    sw     = kv.second;
    sw2    = sumW2[bin_id];
    d      = sumW_data[bin_id];
    u      = sumW_up[bin_id];
    n      = sumW_dn[bin_id];
    a      = sumW_acc[bin_id];
    t.Fill();
  }

  t.Write();
  fo.Close();
  printf("[reduce] wrote %s with %zu bins\n", outFile, sumW.size());
}

/**
 * @brief Command-line driver for reduce_bins()
 *
 * Usage:
 *   ./reduce_bins_driver output.root input1.root input2.root ...
 */
int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " output.root input1.root [input2.root ...]\n";
        return 1;
    }

    const char* outFile = argv[1];
    std::vector<std::string> inputs;
    for (int i = 2; i < argc; ++i)
        inputs.emplace_back(argv[i]);

    reduce_bins(outFile, inputs);
    return 0;
}