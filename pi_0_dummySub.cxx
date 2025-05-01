#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <simc_file>\n";
        return 1;
    }
    const char* simc_file = argv[1];
    std::ifstream in(simc_file);
    if (!in) {
        std::cerr << "Failed to open " << simc_file << "\n";
        return 1;
    }

    // main event loop
    while (true) {
        // 1) read the 16 header values
        double hdr[16];
        for (int i = 0; i < 16; ++i) {
            if (!(in >> hdr[i])) {
                // end-of-file or error -> finish
                return 0;
            }
        }
        // hdr[0] = weight, hdr[1] = sigcc, hdr[2] = sigcm, hdr[3] = dpe, ... hdr[15] = ztarg

        // 2) now read the two photon-cluster lines (6 values each)
        double clusX[2], clusY[2], clusE[2];
        double pX[3], pY[3], pZ[3];
        for (int kk = 0; kk < 2; ++kk) {
            double clusx, clusy, cluse;
            double px, py, pz;
            if (!(in >> clusx >> clusy >> cluse
                     >> px     >> py     >> pz)) {
                throw std::runtime_error(
                    "Malformed SIMC file: unable to read photon cluster #" + std::to_string(kk+1)
                );
            }
            clusX[kk] = clusx;
            clusY[kk] = clusy;
            clusE[kk] = cluse;
            pX[kk+1]  = px;
            pY[kk+1]  = py;
            pZ[kk+1]  = pz;
        }

        // 3) now you have:
        //    hdr[0..15]    -> the 16 SIMC event header values
        //    clusX[0..1],  clusY[0..1],  clusE[0..1]
        //    pX[1..2],     pY[1..2],     pZ[1..2]
        //    ... proceed with fiducial cuts, mass calc, histograms, etc.

        // --- your existing processing goes here ---
    }

    return 0;
}
