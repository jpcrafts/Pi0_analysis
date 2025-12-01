// Simc_import_take5.cxx
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

//
//   g++ -std=c++17 -O2     ./Simc_import_take5.cxx -lyaml-cpp     `root-config --cflags
//     --libs` -lTMVA -lRooFitCore -lRooFit    -o Simc_import_take5
//

// ==== Constants for kin23 (update as needed for each kinematic) ====
const double M_p        = 0.938272;   // Proton mass (GeV)
const double m_e        = 0.000511;   // Electron mass (GeV)
const double E_beam     = 10.539;     // Beam energy (GeV) for kin23
const double NPS_dist   = 407.0;      // NPS face distance from target (cm)
const double theta0_deg = 16.93;      // HMS central angle (deg) for kin23
const double p_central  = 5.253;      // HMS central momentum (GeV) for kin23
const double thp0sv_deg = 13.43;      // NPS central angle (deg) for kin23

//const std::string dir = "/work/hallc/nps/bosted/Simcfiles/";
const std::string dir = "/group/nps/jpcrafts/Pi_0/npa_test/SimcFiles/";

const int NBINS = 200;
const double XMIN = 0.0;
const double XMAX = 5.0;

struct Cluster {
    double x, y, e, t, c1, c2;
};

struct SIMCEvent {
    double weight, sigcc, sigcm, dpe, dphie, dthe, hztar, hdcx, hdcxp, hdcy, hdcyp;
    double hcer, hcal, q2, etotnorm, beta, ph_q;
    Cluster c1, c2;
};

double calc_missing_mass(const SIMCEvent& evt) {
    TLorentzVector pe_in(0, 0, E_beam, sqrt(E_beam*E_beam + m_e*m_e));
    TLorentzVector p_p(0, 0, 0, M_p);

    double p_e = p_central * (1.0 + evt.dpe / 100.0);
    double theta_e = (evt.dthe + theta0_deg) * TMath::DegToRad();
    double phi_e = evt.dphie * TMath::DegToRad();
    double px_e = p_e * sin(theta_e) * cos(phi_e);
    double py_e = p_e * sin(theta_e) * sin(phi_e);
    double pz_e = p_e * cos(theta_e);
    TLorentzVector pe_out(px_e, py_e, pz_e, sqrt(p_e*p_e + m_e*m_e));

    // NPS central angle in radians
    double thp0_rad = thp0sv_deg * TMath::DegToRad();

    // Function to rotate cluster vector by thp0sv about y-axis
    auto make_gamma = [thp0_rad](const Cluster& clus) -> TLorentzVector {
        double x = clus.x;
        double y = clus.y;
        double z = NPS_dist;
        // Rotate (x, y, z) by thp0sv around y-axis
        double x_rot =  x;                                      // untouched
        double y_rot =  y * cos(thp0_rad) + z * sin(thp0_rad);  // +sin
        double z_rot =  z * cos(thp0_rad) - y * sin(thp0_rad);  // −sin

        double norm = sqrt(x_rot*x_rot + y_rot*y_rot + z_rot*z_rot);
        double px = clus.e * x_rot / norm;
        double py = clus.e * y_rot / norm;
        double pz = clus.e * z_rot / norm;
        return TLorentzVector(px, py, pz, clus.e);
    };
    TLorentzVector g1 = make_gamma(evt.c1);
    TLorentzVector g2 = make_gamma(evt.c2);

    TLorentzVector mmvec = pe_in + p_p - pe_out - g1 - g2;
    double mm2 = mmvec.M2();
    return (mm2 > 0) ? sqrt(mm2) : 0.0;
}

typedef std::pair<double, double> mmw_pair;

void read_simc_file(const std::string& fname, std::vector<mmw_pair>& mm_list, double& normfac) {
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
    size_t n_lines = all_lines.size();
    std::istringstream lastiss(all_lines.back());
    double w, scc, scm;
    lastiss >> w >> scc >> scm;
    normfac = scm;

    n_lines--;

    size_t n_events = n_lines / 3;
    int nBad = 0;
    mm_list.reserve(n_events);
    for (size_t i = 0; i < n_events; ++i) {
        std::istringstream iss1(all_lines[3*i]);
        std::istringstream iss2(all_lines[3*i+1]);
        std::istringstream iss3(all_lines[3*i+2]);
        SIMCEvent evt;
        iss1 >> evt.weight >> evt.sigcc >> evt.sigcm >> evt.dpe >> evt.dphie >> evt.dthe
             >> evt.hztar >> evt.hdcx >> evt.hdcxp >> evt.hdcy >> evt.hdcyp
             >> evt.hcer >> evt.hcal >> evt.q2 >> evt.etotnorm >> evt.beta >> evt.ph_q;
        iss2 >> evt.c1.x >> evt.c1.y >> evt.c1.e >> evt.c1.t >> evt.c1.c1 >> evt.c1.c2;
        iss3 >> evt.c2.x >> evt.c2.y >> evt.c2.e >> evt.c2.t >> evt.c2.c1 >> evt.c2.c2;

        // Event selection logic from npa.f:
        double dist = std::sqrt(
            std::pow(evt.c1.x - evt.c2.x, 2) +
            std::pow(evt.c1.y - evt.c2.y, 2)
        );
        if (dist < 15.0) continue;
        if (evt.c1.e < 0.06 || evt.c2.e < 0.06) continue;
        if (evt.c1.e < 0.2  || evt.c2.e < 0.2)  continue;

        double mm = calc_missing_mass(evt);
        if (mm > 0) mm_list.emplace_back(mm, evt.weight);
    }
    std::cout << "[INFO] " << fname << ": " << n_events << " event blocks processed (" << nBad << " skipped)\n";
    std::cout << "        normfac = " << normfac << "\n";
    if (!mm_list.empty()) {
        std::cout << "        MM sample: ";
        for (size_t i = 0; i < std::min(mm_list.size(), size_t(4)); ++i)
            std::cout << mm_list[i].first << " ";
        std::cout << (mm_list.size() > 4 ? "..." : "") << "\n";
    }
}

TH1F* make_hist(const std::vector<mmw_pair>& mm_list, double normfac, const std::string& name) {
    TH1F* h = new TH1F(name.c_str(), (name + ";Missing Mass [GeV];Counts").c_str(), NBINS, XMIN, XMAX);
    for (const auto& pr : mm_list) h->Fill(pr.first, pr.second);
    h->Scale(normfac / 1000.0); // normfac per mC → per µC
    return h;
}

void save_hist_png(TH1F* h, const std::string& png_name) {
    TCanvas* c = new TCanvas(("c_" + png_name).c_str(), png_name.c_str(), 900, 600);
    h->Draw("HIST");
    c->SaveAs(png_name.c_str());
    delete c;
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

    std::vector<mmw_pair> mm_excl, mm_delta, mm_semi;
    double norm_excl = 1.0, norm_delta = 1.0, norm_semi = 1.0;

    read_simc_file(fname_excl,  mm_excl,  norm_excl);
    read_simc_file(fname_delta, mm_delta, norm_delta);
    read_simc_file(fname_semi,  mm_semi,  norm_semi);

    TH1F* hExcl  = make_hist(mm_excl,  norm_excl,  "MM_exclusive");
    TH1F* hDelta = make_hist(mm_delta, norm_delta, "MM_delta");
    TH1F* hSemi  = make_hist(mm_semi,  norm_semi,  "MM_semi");

    save_hist_png(hExcl,  outbase + "_exclusive.png");
    save_hist_png(hDelta, outbase + "_delta.png");
    save_hist_png(hSemi,  outbase + "_semi.png");

    TH1F* hTotal = new TH1F("MM_total", "Total Missing Mass;Missing Mass [GeV];Counts", NBINS, XMIN, XMAX);
    hTotal->Add(hExcl);
    hTotal->Add(hDelta);
    hTotal->Add(hSemi);

    save_hist_png(hTotal, outbase + "_total.png");

    std::cout << "Done. Plots saved as "
              << outbase << "_exclusive.png, "
              << outbase << "_delta.png, "
              << outbase << "_semi.png, "
              << outbase << "_total.png\n";

    delete hExcl; delete hDelta; delete hSemi; delete hTotal;
    return 0;
}
