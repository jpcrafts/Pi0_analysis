#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "TH1D.h"
#include "TCanvas.h"

// simple POD for cluster
struct Cluster { double x, y, E; };

// invariant gamma-gamma mass calculation
static double m_gg(const Cluster& a, const Cluster& b, double Ldet) {
    auto dir = [&](double x, double y) {
        double R = std::sqrt(x*x + y*y + Ldet*Ldet);
        return std::array<double,3>{ x/R, y/R, Ldet/R };
    };
    auto u = dir(a.x, a.y);
    auto v = dir(b.x, b.y);
    double cos12 = std::clamp(u[0]*v[0] + u[1]*v[1] + u[2]*v[2],
                              -1.0, +1.0);
    double m2 = 2.0 * a.E * b.E * (1.0 - cos12);
    return m2 > 0.0 ? std::sqrt(m2) : 0.0;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <simc_flat_txt>\n";
        return 1;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "Failed to open " << argv[1] << "\n";
        return 2;
    }

    // analysis cuts & constants
    const double Emin   = 0.60;   // GeV
    const double minSep = 15.0;   // cm
    const double Ldet   = 407.0;  // cm

    // histogram for π0 invariant mass
    TH1D* h_pi0 = new TH1D("h_pi0",
        "#gamma#gamma invariant mass;M_{#gamma#gamma} [GeV];Yield (/µC)",
        150, 0.0, 0.3);

    long nEv = 0, nSel = 0;

    // --- Main event loop ---
    while (true) {
        double hdr[16];
        for (int i = 0; i < 16; ++i) {
            if (!(in >> hdr[i])) goto _finish;  // EOF
        }

        double w_evt = hdr[0];
        if (w_evt < 0) {
            // trailer detected, skip rest
            break;
        }

        Cluster cl[2];
        double px, py, pz;
        for (int k = 0; k < 2; ++k) {
            if (!(in >> cl[k].x >> cl[k].y >> cl[k].E >> px >> py >> pz)) {
                std::cerr << "Error reading cluster lines!\n";
                goto _finish;
            }
        }
        ++nEv;

        double sep = std::hypot(cl[0].x - cl[1].x,
                                cl[0].y - cl[1].y);
        if (cl[0].E < Emin || cl[1].E < Emin) continue;
        if (sep          < minSep)       continue;

        double m = m_gg(cl[0], cl[1], Ldet);
        double weight = w_evt * 1e6;  // No scaling needed!

        h_pi0->Fill(m, weight);
        ++nSel;
    }

_finish:

    std::cout << "Processed " << nEv
              << " events, " << nSel
              << " passed cuts.\n";

    // --- Save output ---
    TCanvas c("c", "#pi^{0} mass", 800, 600);
    h_pi0->SetLineColor(kRed);
    h_pi0->Draw("hist");
    c.Print("pi0_invariant_mass.pdf");
    std::cout << "Wrote pi0_invariant_mass.pdf\n";

    return 0;
}
