// simc_mm_compare.cpp
// Usage: ./simc_mm_compare simc_16_LH2_excl.txt simc_16_LH2_delta.txt simc_16_LH2_semi.txt output_base
// Input files are looked for in /work/hallc/nps/bosted/Simcfiles/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>

// Constants (as in npa.f and SIMC)
const double M_p     = 0.938272;   // Proton mass (GeV)
const double m_e     = 0.000511;   // Electron mass (GeV)
const double E_beam  = 10.540;     // Beam energy (GeV)
const double NPS_dist = 407.0;     // NPS face distance from target (cm)
const double theta0_deg = 16.48;   // HMS central angle (deg)
const double p_central = 5.878;    // HMS central momentum (GeV)

// Directory prefix for input files
const std::string dir = "/work/hallc/nps/bosted/Simcfiles/";

struct Cluster {
    double x, y, e, t, c1, c2;
};

struct SIMCEvent {
    double weight, sigcc, sigcm, dpe, dphie, dthe, hztar, hdcx, hdcxp, hdcy, hdcyp;
    double hcer, hcal, q2, etotnorm, beta, ph_q;
    Cluster c1, c2;
};

double calc_missing_mass(const SIMCEvent& evt) {
    // Incoming electron
    TLorentzVector pe_in(0, 0, E_beam, sqrt(E_beam*E_beam + m_e*m_e));
    // Target proton at rest
    TLorentzVector p_p(0, 0, 0, M_p);

    // Outgoing e- (scattered)
    double p_e = p_central * (1.0 + evt.dpe / 100.0);
    double theta_e = (evt.dthe + theta0_deg) * TMath::DegToRad(); // central angle + delta
    double phi_e = evt.dphie * TMath::DegToRad();
    double px_e = p_e * sin(theta_e) * cos(phi_e);
    double py_e = p_e * sin(theta_e) * sin(phi_e);
    double pz_e = p_e * cos(theta_e);
    TLorentzVector pe_out(px_e, py_e, pz_e, sqrt(p_e*p_e + m_e*m_e));

    // Photons from clusters (NPS face at z = NPS_dist)
    auto make_gamma = [](const Cluster& clus) -> TLorentzVector {
        double x = clus.x;
        double y = clus.y;
        double z = NPS_dist;
        double norm = sqrt(x*x + y*y + z*z);
        double px = clus.e * x / norm;
        double py = clus.e * y / norm;
        double pz = clus.e * z / norm;
        return TLorentzVector(px, py, pz, clus.e);
    };
    TLorentzVector g1 = make_gamma(evt.c1);
    TLorentzVector g2 = make_gamma(evt.c2);

    // Missing mass calculation
    TLorentzVector mmvec = pe_in + p_p - pe_out - g1 - g2;
    double mm2 = mmvec.M2();
    return (mm2 > 0) ? sqrt(mm2) : 0.0;
}

// Reads all but last line block as events; last line is for normalization
void read_simc_file(const std::string& fname, std::vector<double>& mm_list, double& normfac) {
    std::ifstream fin(fname);
    if (!fin) {
        std::cerr << "ERROR: Could not open file: " << fname << std::endl;
        exit(1);
    }
    std::vector<std::string> all_lines;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        all_lines.push_back(line);
    }
    fin.close();
    if (all_lines.size() < 4) {
        std::cerr << "ERROR: File " << fname << " is empty or incomplete!\n";
        exit(2);
    }
    // Last line for normalization (could be three lines at end; handle robustly)
    size_t n_lines = all_lines.size();
    // Last event block is 3 lines before the last line; last line is normalization
    std::istringstream lastiss(all_lines.back());
    double w, scc, scm;
    lastiss >> w >> scc >> scm;
    normfac = scm; // Fortran uses sigcm as normfac

    // Remove the last line from processing
    n_lines--;

    // Now process all remaining lines in blocks of 3
    size_t n_events = n_lines / 3;
    int nBad = 0;
    mm_list.reserve(n_events);
    for (size_t i = 0; i < n_events; ++i) {
        std::istringstream iss1(all_lines[3*i]);
        std::istringstream iss2(all_lines[3*i+1]);
        std::istringstream iss3(all_lines[3*i+2]);
        SIMCEvent evt;
        // Parse global event variables
        iss1 >> evt.weight >> evt.sigcc >> evt.sigcm >> evt.dpe >> evt.dphie >> evt.dthe
             >> evt.hztar >> evt.hdcx >> evt.hdcxp >> evt.hdcy >> evt.hdcyp
             >> evt.hcer >> evt.hcal >> evt.q2 >> evt.etotnorm >> evt.beta >> evt.ph_q;
        // Parse cluster 1
        iss2 >> evt.c1.x >> evt.c1.y >> evt.c1.e >> evt.c1.t >> evt.c1.c1 >> evt.c1.c2;
        // Parse cluster 2
        iss3 >> evt.c2.x >> evt.c2.y >> evt.c2.e >> evt.c2.t >> evt.c2.c1 >> evt.c2.c2;

        // --- SIMC event selection logic (verbatim from npa.f) ---
        // 1. Cluster distance cut (15cm)
        double dist = std::sqrt(
            std::pow(evt.c1.x - evt.c2.x, 2) +
            std::pow(evt.c1.y - evt.c2.y, 2)
        );
        if (dist < 15.0) continue;

        // 2. Minimum photon energy cut (0.06 GeV)
        if (evt.c1.e < 0.06 || evt.c2.e < 0.06) continue;

        // 3. Minimum cluster energy cut (0.2 GeV)
        if (evt.c1.e < 0.2  || evt.c2.e < 0.2)  continue;
        // --------------------------------------------------------

        double mm = calc_missing_mass(evt);
        if (mm > 0) mm_list.push_back(mm);
    }
    std::cout << "[INFO] " << fname << ": " << n_events << " event blocks processed (" << nBad << " skipped)\n";
    std::cout << "        normfac = " << normfac << "\n";
    std::cout << "        MM sample: ";
    for (size_t i = 0; i < std::min(mm_list.size(), size_t(4)); ++i) std::cout << mm_list[i] << " ";
    std::cout << (mm_list.size() > 4 ? "..." : "") << "\n";
}

void fill_and_save_hist(const std::vector<double>& mm_list, double normfac, const std::string& hist_name, const std::string& png_name) {
    // Binning/range as before (0 to 2.5 GeV, 200 bins)
    TH1F* hMM = new TH1F(hist_name.c_str(), (hist_name + ";Missing Mass [GeV];Counts").c_str(), 200, 0, 2.5);
    for (const auto& mm : mm_list) {
        hMM->Fill(mm, 1.0);
    }
    hMM->Scale(normfac); // Apply normalization

    // Diagnostics
    std::cout << "Histogram: " << hist_name << "\n";
    std::cout << "  normfac = " << normfac << "\n";
    std::cout << "  Integral (after norm) = " << hMM->Integral() << "\n";
    std::cout << "  Max bin content = " << hMM->GetMaximum() << "\n\n";

    // Save as PNG
    TCanvas* c1 = new TCanvas(("c_" + hist_name).c_str(), hist_name.c_str(), 900, 600);
    hMM->Draw("HIST");
    c1->SaveAs(png_name.c_str());
    delete c1;
    delete hMM;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cout << "Usage: " << argv[0] << " excl.txt delta.txt semi.txt output_base\n";
        std::cout << "Input files will be searched in: " << dir << "\n";
        return 0;
    }
    std::string fname_excl  = dir + argv[1];
    std::string fname_delta = dir + argv[2];
    std::string fname_semi  = dir + argv[3];
    std::string outbase     = argv[4];

    std::vector<double> mm_excl, mm_delta, mm_semi;
    double norm_excl = 1.0, norm_delta = 1.0, norm_semi = 1.0;

    read_simc_file(fname_excl,  mm_excl,  norm_excl);
    read_simc_file(fname_delta, mm_delta, norm_delta);
    read_simc_file(fname_semi,  mm_semi,  norm_semi);

    fill_and_save_hist(mm_excl,  norm_excl,  "MM_exclusive",    outbase + "_exclusive.png");
    fill_and_save_hist(mm_delta, norm_delta, "MM_delta",        outbase + "_delta.png");
    fill_and_save_hist(mm_semi,  norm_semi,  "MM_semi",         outbase + "_semi.png");

    std::cout << "Done. Plots saved as " << outbase << "_exclusive.png, " << outbase << "_delta.png, " << outbase << "_semi.png\n";
    return 0;
}
