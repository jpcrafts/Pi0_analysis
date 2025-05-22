// Simc_import_take4.cxx
// Clean implementation of SIMC missing mass analysis (LH2), guided by npa.f
// Usage: ./Simc_import_take4 <input1> [<input2> ...] -o <output.root>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>

// ROOT headers
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

// Structure for cluster info (energy, momentum)
struct Cluster {
    double x, y, E, px, py, pz;
};

// Structure for event data
// Structure for event data, matching SIMC input in npa.f
struct Event {
    // First line
    double weight, sigcc, sigcm;
    double dpe, dphie, dthe, hztar, hdcx, hdcxp, hdcy, hdcyp;
    double dpe_init, xptar_init, yptar_init, yt_orig, yrast, ztarg;
    // Two clusters
    Cluster cl[2];
};

// Function to parse a SIMC event from input stream (to be implemented)
bool read_event(std::istream& in, Event& evt) {
    // Read the first line (17 values)
    if (!(in >> evt.weight >> evt.sigcc >> evt.sigcm
            >> evt.dpe >> evt.dphie >> evt.dthe >> evt.hztar
            >> evt.hdcx >> evt.hdcxp >> evt.hdcy >> evt.hdcyp
            >> evt.dpe_init >> evt.xptar_init >> evt.yptar_init
            >> evt.yt_orig >> evt.yrast >> evt.ztarg)) {
        return false;
    }
    // Read two cluster lines (each 6 values)
    for (int kk = 0; kk < 2; ++kk) {
        if (!(in >> evt.cl[kk].x >> evt.cl[kk].y >> evt.cl[kk].E
                >> evt.cl[kk].px >> evt.cl[kk].py >> evt.cl[kk].pz)) {
            return false;
        }
    }
    // Apply Fortran scaling: dthe and dphie are stored 100x larger in the file
    evt.dthe /= 100.0;
    evt.dphie /= 100.0;
    return true;
}

// Missing mass calculation (from npa.f logic)
double missing_mass2(const Event& evt, double e0, double ep0, double c1 = 1.0, double c2 = 1.0) {
    // Calculate scattered electron energy
    double ep_e = ep0 * (1.0 + evt.dpe / 100.0);

    // Calculate electron momentum components from angles (approximate, assuming small angles)
    // dthe = xptar, dphie = yptar (in radians)
    double theta_e = evt.dthe;   // xptar, radians
    double phi_e = evt.dphie;    // yptar, radians
    double pe = ep_e;            // Assume massless electron, p = E

    double px_e = pe * theta_e;  // Approximate px (for small angles, sin(theta) ~ theta)
    double py_e = pe * phi_e;    // Approximate py
    double pz_e = std::sqrt(std::max(0.0, pe * pe - px_e * px_e - py_e * py_e)); // pz

    // Cluster momenta (no geometric corrections yet)
    double reco_px[2] = {evt.cl[0].px, evt.cl[1].px};
    double reco_py[2] = {evt.cl[0].py, evt.cl[1].py};
    double reco_pz[2] = {evt.cl[0].pz, evt.cl[1].pz};

    // Use proton mass for am (LH2 target)
    const double mp = 0.938272; // GeV/c^2

    // Fortran logic: mm2_0, mm2_1, mm2_2, correction
    auto mm2_formula = [&](double c1v, double c2v) {
        double Etot = (e0 + mp) - (ep_e + evt.cl[0].E * c1v + evt.cl[1].E * c2v);
        double Px = px_e + reco_px[0] * c1v + reco_px[1] * c2v;
        double Py = py_e + reco_py[0] * c1v + reco_py[1] * c2v;
        double Pz = pz_e + reco_pz[0] * c1v + reco_pz[1] * c2v - e0;
        return Etot * Etot - (Px * Px + Py * Py + Pz * Pz);
    };

    double mm2_0 = mm2_formula(1.0, 1.0);
    double mm2_1 = mm2_formula(1.01, 1.0);
    double mm2_2 = mm2_formula(1.0, 1.01);

    double dmm2dc1 = (mm2_1 - mm2_0) / 0.01;
    double dmm2dc2 = (mm2_2 - mm2_0) / 0.01;

    // For now, use pion mass for ampi (as in Fortran)
    const double ampi = 0.134977; // GeV/c^2 (neutral pion mass)
    // m2 = sqrt(etot^2 - (p1 + p2)^2) -- for correction, use mm2_0 as m2^2
    double m2 = std::max(0.0, mm2_0);
    double dc1 = (ampi * ampi - m2) / m2 / (1.0 + dmm2dc1 / dmm2dc2);
    double dc2 = dc1 * dmm2dc1 / dmm2dc2;
    double mm2_corr = mm2_0 + dmm2dc1 * dc1 + dmm2dc2 * dc2;

    return mm2_corr;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input1> [<input2> ...] -o <output.root>\n";
        return 1;
    }

    // Parse arguments
    std::vector<std::string> input_files;
    std::string output_file;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else {
            input_files.push_back(argv[i]);
        }
    }
    if (output_file.empty() || input_files.empty()) {
        std::cerr << "Error: Must specify input file(s) and output file.\n";
        return 1;
    }

    // Set up beam and electron energies (user should adjust as needed)
    const double e0 = 4.4;   // Beam energy in GeV (example)
    const double ep0 = 2.2;  // Central scattered electron energy in GeV (example)

    // Set up ROOT output
    TFile* fout = new TFile(output_file.c_str(), "RECREATE");
    TH1D* h_mm2 = new TH1D("h_mm2", "Missing Mass^2;mm^{2} (GeV^{2});Counts", 400, -1, 3);

    // Loop over input files
    // Directory prefix for input files
    const std::string dir = "/work/hallc/nps/bosted/Simcfiles/";

    for (const auto& fname : input_files) {
        std::ifstream fin(dir + fname);
        if (!fin) {
            std::cerr << "Error: Cannot open " << dir + fname << "\n";
            continue;
        }
        std::cout << "Opened file: " << dir + fname << std::endl;
        int total_events = 0;
        int passed_events = 0;
        Event evt;
        bool first_failed_printed = false;
        while (read_event(fin, evt)) {
            ++total_events;
            if (total_events == 1) {
                std::cout << "First event in file " << fname << ":\n";
                std::cout << "  weight=" << evt.weight
                          << ", sigcc=" << evt.sigcc
                          << ", sigcm=" << evt.sigcm
                          << ", dpe=" << evt.dpe
                          << ", dphie=" << evt.dphie
                          << ", dthe=" << evt.dthe
                          << ", hztar=" << evt.hztar << "\n";
                std::cout << "  Cluster 0: x=" << evt.cl[0].x << ", y=" << evt.cl[0].y
                          << ", E=" << evt.cl[0].E << ", px=" << evt.cl[0].px
                          << ", py=" << evt.cl[0].py << ", pz=" << evt.cl[0].pz << "\n";
                std::cout << "  Cluster 1: x=" << evt.cl[1].x << ", y=" << evt.cl[1].y
                          << ", E=" << evt.cl[1].E << ", px=" << evt.cl[1].px
                          << ", py=" << evt.cl[1].py << ", pz=" << evt.cl[1].pz << "\n";
            }
            // Event selection cuts (from NPA.f lines 1965–1972)
            // Calculate etot (sum of cluster energies)
            double etot = evt.cl[0].E + evt.cl[1].E;
            // Calculate rcl (distance between clusters in x/y)
            double rcl = std::sqrt(
                std::pow(evt.cl[0].x - evt.cl[1].x, 2) +
                std::pow(evt.cl[0].y - evt.cl[1].y, 2)
            );
            // Cut values from NPA.f lines 304-306
            const double egmin = 1.0;       // GeV
            const double rclmin = 15.0;     // cm
            const double etotmin = 2.0;     // GeV

            // Apply cuts
            bool passed_cuts =
                evt.dpe > -9.0 && evt.dpe < 11.0 &&
                std::abs(evt.dphie) < 0.025 &&
                std::abs(evt.dthe) < 0.060 &&
                std::abs(evt.hztar) < 8.0 &&
                // Cluster position cuts (okcl) and dead block logic omitted for now
                rcl > rclmin &&
                etot > etotmin &&
                evt.cl[0].E > egmin &&
                evt.cl[1].E > egmin;

            if (passed_cuts) {
                ++passed_events;
                // Calculate missing mass squared
                double mm2 = missing_mass2(evt, e0, ep0);

                // Fill histogram
                h_mm2->Fill(mm2);
            } else if (!first_failed_printed) {
                std::cout << "First event in file " << fname << " that failed selection cuts:\n";
                std::cout << "  weight=" << evt.weight
                          << ", sigcc=" << evt.sigcc
                          << ", sigcm=" << evt.sigcm
                          << ", dpe=" << evt.dpe
                          << ", dphie=" << evt.dphie
                          << ", dthe=" << evt.dthe
                          << ", hztar=" << evt.hztar << "\n";
                std::cout << "  Cluster 0: x=" << evt.cl[0].x << ", y=" << evt.cl[0].y
                          << ", E=" << evt.cl[0].E << ", px=" << evt.cl[0].px
                          << ", py=" << evt.cl[0].py << ", pz=" << evt.cl[0].pz << "\n";
                std::cout << "  Cluster 1: x=" << evt.cl[1].x << ", y=" << evt.cl[1].y
                          << ", E=" << evt.cl[1].E << ", px=" << evt.cl[1].px
                          << ", py=" << evt.cl[1].py << ", pz=" << evt.cl[1].pz << "\n";
                std::cout << "  Failed criteria:";
                if (!(evt.dpe > -9.0 && evt.dpe < 11.0)) std::cout << " dpe";
                if (!(std::abs(evt.dphie) < 0.025)) std::cout << " dphie";
                if (!(std::abs(evt.dthe) < 0.060)) std::cout << " dthe";
                if (!(std::abs(evt.hztar) < 8.0)) std::cout << " hztar";
                if (!(rcl > rclmin)) std::cout << " rcl";
                if (!(etot > etotmin)) std::cout << " etot";
                if (!(evt.cl[0].E > egmin)) std::cout << " cl[0].E";
                if (!(evt.cl[1].E > egmin)) std::cout << " cl[1].E";
                std::cout << std::endl;
                first_failed_printed = true;
            }
        }
        std::cout << "File: " << fname << " -- Total events read: " << total_events
                  << ", Events passing cuts: " << passed_events << std::endl;
        if (total_events == 0) {
            std::cerr << "Warning: No events read from file " << fname << std::endl;
        }
        if (passed_events == 0 && total_events > 0) {
            std::cerr << "Warning: No events passed selection cuts in file " << fname << std::endl;
        }
    }

    // Write and close
    fout->Write();
    fout->Close();

    std::cout << "Analysis complete. Output written to " << output_file << "\n";
    return 0;
}