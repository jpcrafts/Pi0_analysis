#include <cstdio>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <iostream>

// ROOT headers
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>

// constants
constexpr double DNPS   = 407.0;   // detector distance in cm
constexpr double egmin  = 0.6;     // GeV, min single‐photon energy
constexpr double rclmin = 15.0;    // cm, min cluster separation
constexpr double mpi0   = 0.135;   // GeV

// compute photon 4-vector from cluster
inline std::array<double,3> photon_p(double E, double x, double y) {
    double L = std::sqrt(x*x + y*y + DNPS*DNPS);
    return { E*x/L, E*y/L, E*DNPS/L };
}

int main(int argc, char** argv) {
    std::vector<std::string> procs = { "excl", "dlta", "semi" };
    std::vector<std::string> targets = { "DUM" };   // just LH2
    int ikin = 17;                                   // your chosen kinematic setting

    // create output histogram
    TH1D *h_mgg = new TH1D("h_mgg",
        "π^{0} → γγ Invariant Mass;M_{#gamma#gamma} [GeV];Counts",
        200, 0.0, 0.3);

    for (auto &t : targets) {
      for (auto &p : procs) {
        char fname[256];
        snprintf(fname, sizeof(fname),
          "/work/hallc/nps/bosted/Simcfiles/simc_%d_%s_%s.txt",
          ikin, t.c_str(), p.c_str());
        std::ifstream in(fname);
        if (!in) {
          std::cerr << ">> failed to open " << fname << "\n";
          continue;
        }
        double weight, sigcc, sigcm;
        double dpe, dphie, dthe, hztar, hdcx, hdcxp, hdcy, hdcyp;
        // (we ignore the initial summary lines)
        double normfac=1, sum3=0;
        for (;;) {
          // read until weight == -1 signals summary
          in >> weight;
          if (!in || weight == -1) {
            in >> sigcm; // grab sum3
            normfac = sigcm / sum3;
            break;
          }
          // skip the rest of header
          for (int i=0; i<28; ++i) in >> std::ws; 

          // now event‐by‐event:
          // read 2 clusters: x,y,E and p_x,p_y,p_z
          double clusx[2], clusy[2], clusE[2];
          double px[2], py[2], pz[2];
          for (int k=0;k<2;++k) {
            in >> clusx[k] >> clusy[k] >> clusE[k]
               >> px[k]    >> py[k]    >> pz[k];
          }
          // --- photon cuts ---
          if (clusE[0] <= egmin || clusE[1] <= egmin) continue;
          double dx = clusx[0]-clusx[1], dy = clusy[0]-clusy[1];
          if (std::sqrt(dx*dx+dy*dy) < rclmin) continue;

          // compute two‐photon invariant mass
          auto p1 = photon_p(clusE[0], clusx[0], clusy[0]);
          auto p2 = photon_p(clusE[1], clusx[1], clusy[1]);
          double E1 = clusE[0], E2 = clusE[1];
          double dot = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2];
          double cos12 = dot/(E1*E2);
          cos12 = std::clamp(cos12, -1.0, 1.0);
          double mgg2 = 2.0*E1*E2*(1.0 - cos12);
          if (mgg2<0) continue;
          double mgg = std::sqrt(mgg2);

          // fill histogram
          h_mgg->Fill(mgg);
          ++sum3;
        }
        in.close();
      }
    }

    // draw and save
    TCanvas *c = new TCanvas("c","π0 mass",800,600);
    h_mgg->SetLineColor(kBlue);
    h_mgg->Draw();
    c->SaveAs("pi0_mass_hist.png");

    return 0;
}
