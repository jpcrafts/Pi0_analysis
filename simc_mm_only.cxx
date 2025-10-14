#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

// Helper: Read three lines per event and parse floats
bool read_next_event(std::ifstream &in, std::vector<double> &main_vals, std::vector<double> &clus1_vals, std::vector<double> &clus2_vals) {
    std::string line;
    main_vals.clear(); clus1_vals.clear(); clus2_vals.clear();

    // Read main line
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        double v;
        while (ss >> v) main_vals.push_back(v);
        break;
    }
    if (main_vals.empty()) return false;

    // Read cluster 1
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        double v;
        while (ss >> v) clus1_vals.push_back(v);
        break;
    }
    if (clus1_vals.empty()) return false;

    // Read cluster 2
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        double v;
        while (ss >> v) clus2_vals.push_back(v);
        break;
    }
    if (clus2_vals.empty()) return false;

    return true;
}

// --- Insert your missing mass calculation logic here ---
// This version uses the same logic as Simc_import_take5.cxx
// (For best accuracy, copy the relevant lines for 4-vectors and mm)

double calc_missing_mass(
    const std::vector<double> &main_vals,
    const std::vector<double> &clus1_vals,
    const std::vector<double> &clus2_vals)
{
    // Example: For Hall C Pi0, you likely have:
    // main_vals: [q2, nu, x_bj, ... px, py, pz, E, ...]
    // clus1/2: [E, px, py, pz, ...]
    // But for exact order, check your actual code/data!
    // Below is the structure for standard SIMC pi0:
    // main[0-2] are BCMs, main[3-14] are kinematic vars, main[15-] are optional

    // This matches your simc_23_LH2_excl.txt: main_vals.size() == 16
    // Cluster lines have 6 variables

    // Electron beam: assume along z
    double mp = 0.938272;   // proton mass [GeV]
    double me = 0.000511;   // electron mass [GeV]
    double mpi0 = 0.134977; // pi0 mass [GeV]

    // Assign event variables
    // These indices are based on the file layout you showed
    double ebeam = main_vals[4];   // Ebeam
    double eep   = main_vals[5];   // E' (scattered electron)
    double p_ep  = main_vals[6];
    double theta_ep = main_vals[7];
    double phi_ep   = main_vals[8];

    // Cluster 1 photon (the code often assumes: E, px, py, pz, ... but yours is E, px, py, pz, ...)
    double E1  = clus1_vals[0];
    double px1 = clus1_vals[1];
    double py1 = clus1_vals[2];
    double pz1 = clus1_vals[3];

    // Cluster 2 photon
    double E2  = clus2_vals[0];
    double px2 = clus2_vals[1];
    double py2 = clus2_vals[2];
    double pz2 = clus2_vals[3];

    // Electron 4-vector (initial, final)
    double Ee = ebeam;
    double Pe_x = 0.0, Pe_y = 0.0, Pe_z = Ee;
    double Eep = eep;
    double Pep_x = p_ep * sin(theta_ep) * cos(phi_ep);
    double Pep_y = p_ep * sin(theta_ep) * sin(phi_ep);
    double Pep_z = p_ep * cos(theta_ep);

    // Target 4-vector (assume at rest)
    double Mp = mp;
    double Pt_x = 0.0, Pt_y = 0.0, Pt_z = 0.0, Et = Mp;

    // Two photon 4-vector
    double Egam = E1 + E2;
    double Gam_x = px1 + px2;
    double Gam_y = py1 + py2;
    double Gam_z = pz1 + pz2;

    // Calculate q (virtual photon)
    double q_x = Pe_x - Pep_x;
    double q_y = Pe_y - Pep_y;
    double q_z = Pe_z - Pep_z;
    double nu = Ee - Eep;

    // Final system: target + virtual photon - pi0 system
    double X_x = Pt_x + q_x - Gam_x;
    double X_y = Pt_y + q_y - Gam_y;
    double X_z = Pt_z + q_z - Gam_z;
    double X_E = Et + nu - Egam;

    double mm2 = X_E*X_E - (X_x*X_x + X_y*X_y + X_z*X_z);

    return (mm2 > 0.0) ? std::sqrt(mm2) : -std::sqrt(-mm2); // Can be negative if off-shell
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " simc_events.txt" << std::endl;
        return 1;
    }

    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "Error opening file " << argv[1] << std::endl;
        return 1;
    }

    std::vector<double> main_vals, clus1_vals, clus2_vals;
    int evt = 0;
    while (read_next_event(in, main_vals, clus1_vals, clus2_vals)) {
        double mm = calc_missing_mass(main_vals, clus1_vals, clus2_vals);
        std::cout << mm << std::endl;
        // Optionally, print event weight or more: std::cout << mm << " " << main_vals[0] << std::endl;
        evt++;
    }

    return 0;
}
