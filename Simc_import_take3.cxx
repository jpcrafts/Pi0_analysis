/*************************************************************
  Simc_import_take3.cxx   – 2025‑05‑04 (finalised)
*************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"

constexpr double deg = M_PI/180.0;

/* ---------- 56‑point kinematic tables copied from npa.f ---------- */
static const double e0sv[56]={10.540,10.540,10.540,10.540,10.540,10.540,10.540,10.540,
 8.457,8.457,8.457,8.457,8.457,8.457,8.457,10.538,10.540,10.540,10.540,10.540,10.540,10.540,
10.540,10.540,10.540,10.540,10.540,10.540, 6.369, 6.369, 6.369, 6.369, 6.369, 6.369, 6.369,
 6.369, 6.369, 6.369, 8.457, 8.456, 8.456, 6.370, 6.370, 6.371, 6.371, 8.456, 8.456,10.543,
10.543,10.543,10.539,10.539,10.539,10.539,10.539,10.539};
static const double ep0sv[56]={5.890,6.110,4.640,4.640,5.233,5.233,5.233,5.878,4.042,4.042,
4.042,4.726,4.726,3.803,3.803,4.637,6.667,6.667,5.253,5.253,4.637,5.878,5.878,5.038,5.038,
2.416,2.416,4.149,2.638,2.638,2.638,2.638,2.638,1.719,1.734,1.956,3.805,2.562,2.131,4.726,
2.638,2.638,1.956,1.719,2.638,4.726,4.042,3.803,6.667,5.878,5.253,5.038,4.637,4.149,2.416,
6.117};
static const double the0sv[56]={16.480,12.370,16.440,16.433,16.910,16.910,16.910,16.483,
17.010,17.010,17.010,16.750,16.750,22.920,22.920,16.440,12.480,12.480,16.930,16.930,16.440,
16.480,16.480,19.350,19.350,26.800,26.800,15.200,25.930,25.930,25.930,25.930,25.930,39.810,
25.120,28.320,22.940,24.780,23.700,16.740,25.940,25.940,28.340,39.810,25.940,16.740,17.000,
22.920,12.490,16.480,16.910,19.360,16.430,15.200,26.850,12.500};
/* ----------------------------------------------------------------- */

struct Cluster{double x,y,E,px,py,pz;};

static double m_gg(const Cluster&a,const Cluster&b,double Ldet){
    auto dir=[&](double x,double y){double R=std::sqrt(x*x+y*y+Ldet*Ldet);
        return std::array<double,3>{x/R,y/R,Ldet/R};};
    auto u=dir(a.x,a.y), v=dir(b.x,b.y);
    double c=std::clamp(u[0]*v[0]+u[1]*v[1]+u[2]*v[2],-1.0,1.0);
    double m2=2*a.E*b.E*(1-c); return m2>0? std::sqrt(m2):0.0;
}

enum class Proc{EXCL,DELTA,SEMI,UNK};
static Proc proc_from(const std::string&s){
    if(s.find("_excl")!=std::string::npos) return Proc::EXCL;
    if(s.find("_dlta")!=std::string::npos||s.find("_delta")!=std::string::npos) return Proc::DELTA;
    if(s.find("_semi")!=std::string::npos) return Proc::SEMI;
    return Proc::UNK;
}
static int kin_idx(const std::string&f){
    auto p=f.find("simc_"); if(p==std::string::npos) return 0; p+=5;
    auto q=f.find('_',p); if(q==std::string::npos) return 0;
    return std::stoi(f.substr(p,q-p));
}
static double read_sigma(const std::string&p){
    std::ifstream fin(p); std::string l,last;
    while(std::getline(fin,l)) if(!l.empty()) last=l;
    std::istringstream iss(last); double a,b,s=0; iss>>a>>b>>s; return s;
}

int main(int argc,char**argv){
    if(argc<2){std::cerr<<"Usage: "<<argv[0]<<" simc_files\n";return 1;}
    const std::string dir="/work/hallc/nps/bosted/Simcfiles/";
    const double mp=0.938272, Ldet=407.0, Emin=0.6, minSep=15, etotmin=1.5, scale=3000;

    /* histograms */
    TH1D *h_pi[4]{new TH1D("h_excl","Exclusive;M_{#gamma#gamma} [GeV];Yield (/µC)",150,0,0.3),
                  new TH1D("h_delta","#Delta;M_{#gamma#gamma} [GeV];Yield (/µC)",150,0,0.3),
                  new TH1D("h_semi","Semi‑incl.;M_{#gamma#gamma} [GeV];Yield (/µC)",150,0,0.3),
                  new TH1D("h_all" ,"Combined;M_{#gamma#gamma} [GeV];Yield (/µC)",150,0,0.3)};
    TH1D *h_mm[4]{new TH1D("mm_excl","Exclusive;M_{miss} [GeV];Yield (/µC)",150,0,2.0),
                  new TH1D("mm_delta","#Delta;M_{miss} [GeV];Yield (/µC)",150,0,2.0),
                  new TH1D("mm_semi","Semi‑incl.;M_{miss} [GeV];Yield (/µC)",150,0,2.0),
                  new TH1D("mm_all" ,"Combined;M_{miss} [GeV];Yield (/µC)",150,0,2.0)};

    for(int a=1;a<argc;++a){
        std::string path=argv[a]; if(path.find('/')==std::string::npos) path=dir+path;
        Proc P=proc_from(path); int k=kin_idx(path);
        if(k<1||k>56){std::cerr<<"kin idx fail "<<path<<"\n";continue;}
        double e0=e0sv[k-1], ep0=ep0sv[k-1], th0=the0sv[k-1]*deg;

        long N=0; {std::ifstream in(path); double w;
            while(in>>w){double h; for(int i=0;i<15;++i) in>>h;
                Cluster c; for(int i=0;i<2;++i) in>>c.x>>c.y>>c.E>>c.px>>c.py>>c.pz;
                if(!in) break; ++N;}}
        double sigma=read_sigma(path); if(sigma<=0){std::cerr<<"Bad trailer "<<path<<"\n";continue;}
        double norm=sigma/N; if(P==Proc::EXCL||P==Proc::DELTA) norm*=0.7;

        std::ifstream in(path); long kept=0;
        while(true){
            double w_evt; if(!(in>>w_evt)) break;
            double h[15]; for(int i=0;i<15;++i) in>>h[i];
            Cluster cl[2]; for(int i=0;i<2;++i) in>>cl[i].x>>cl[i].y>>cl[i].E>>cl[i].px>>cl[i].py>>cl[i].pz;
            if(!in) break;

            // (Removed mm to cm conversion; assume input is already in cm)

            double dpe=h[3];             // %
            double xptar=h[4];            // radians already
            double yptar=h[5];            // radians
            // Apply Fortran-like scaling: divide by 100 to match NPA.f logic
            xptar /= 100.0;
            yptar /= 100.0;
            double hztar = h[6];          // target z position

            // Set thp0r as a constant for your kinematic setting (example: 18 degrees)
            // TODO: Set this value appropriately for your data!
            constexpr double thp0r = 18.0 * M_PI / 180.0; // replace 18.0 with your actual angle

            // Set DTARG_CONST as appropriate for your run (example: 307.0)
            // TODO: Set this value appropriately for your data!
            constexpr double DTARG_CONST = 307.0; // or 357.0, 407.0, etc.

            // Compute dtarg per event as in NPA.f
            double dtarg = DTARG_CONST - hztar * std::cos(thp0r);

            // dpe is fractional deviation in percent; reconstruct scattered electron energy
            double ep_e = ep0 * (1.0 + dpe / 100.0);

            // Reconstruct electron kinematics (physics_angles logic)
            double Nvec=std::sqrt(1+xptar*xptar+yptar*yptar);
            double ux=(std::sin(th0)+xptar*std::cos(th0))/Nvec;
            double uy= yptar/Nvec;
            double uz=(std::cos(th0)-xptar*std::sin(th0))/Nvec;
            double px_e=ep_e*ux, py_e=ep_e*uy, pz_e=ep_e*uz;

            // Reconstruct cluster momenta from positions and energies (NPA.f lines 1940–1946)
            struct RecoMom {
                double px, py, pz;
            } reco[2];
            for(int kk=0; kk<2; ++kk) {
                double xcorr = cl[kk].x - hztar;
                double pz0 = cl[kk].E * std::sqrt(dtarg*dtarg - xcorr*xcorr - cl[kk].y*cl[kk].y) / dtarg;
                double px = -pz0 * cl[kk].y / dtarg;
                double py0 = pz0 * xcorr / dtarg;
                double py = py0 * std::cos(thp0r) + pz0 * std::sin(thp0r);
                double pz = pz0 * std::cos(thp0r) - py0 * std::sin(thp0r);
                reco[kk] = {px, py, pz};
            }

            static int debug_cut_count = 0;
            bool cutfail = false;
            if(debug_cut_count < 20) {
                std::cout << "Event " << debug_cut_count+1 << ":\n";
                std::cout << "  dpe = " << dpe << "\n";
                std::cout << "  xptar (dthe) = " << xptar << "\n";
                std::cout << "  yptar (dphie) = " << yptar << "\n";
                std::cout << "  hztar = " << hztar << "\n";
                std::cout << "  cl[0].E = " << cl[0].E << ", cl[1].E = " << cl[1].E << "\n";
                std::cout << "  etot = " << cl[0].E + cl[1].E << "\n";
                std::cout << "  sep = " << std::hypot(cl[0].x-cl[1].x,cl[0].y-cl[1].y) << "\n";
            }
            // Restore yptar and xptar cuts (after division by 100)
            if(std::abs(yptar) >= 0.025) { if(debug_cut_count < 20) { std::cout << "  FAIL: dphie (yptar)\n"; ++debug_cut_count; } continue; }
            if(std::abs(xptar) >= 0.060) { if(debug_cut_count < 20) { std::cout << "  FAIL: dthe (xptar)\n"; ++debug_cut_count; } continue; }
            if(cl[0].E<Emin||cl[1].E<Emin) { if(debug_cut_count < 20) { std::cout << "  FAIL: Emin\n"; ++debug_cut_count; } continue; }
            if(std::hypot(cl[0].x-cl[1].x,cl[0].y-cl[1].y)<minSep) { if(debug_cut_count < 20) { std::cout << "  FAIL: minSep\n"; ++debug_cut_count; } continue; }
            if(dpe <= -9.0 || dpe >= 11.0) { if(debug_cut_count < 20) { std::cout << "  FAIL: dpe\n"; ++debug_cut_count; } continue; }
            if(std::abs(yptar) >= 0.025) { if(debug_cut_count < 20) { std::cout << "  FAIL: dphie (yptar)\n"; ++debug_cut_count; } continue; }
            if(std::abs(xptar) >= 0.060) { if(debug_cut_count < 20) { std::cout << "  FAIL: dthe (xptar)\n"; ++debug_cut_count; } continue; }
            if(std::abs(hztar) >= 8.0) { if(debug_cut_count < 20) { std::cout << "  FAIL: hztar\n"; ++debug_cut_count; } continue; }
            double etot = cl[0].E + cl[1].E;
            if(etot < etotmin) { if(debug_cut_count < 20) { std::cout << "  FAIL: etotmin\n"; ++debug_cut_count; } continue; }

            double mpi=m_gg(cl[0],cl[1],Ldet);

            // Missing mass calculation as in NPA.f (mm2_0, lines 2035–2038)
            // Apply cluster energy/momentum scaling as in NPA.f (c1=1.01, c2=1.0)
            double c1 = 1.0, c2 = 1.0;
            double Etot = (e0 + mp) - (ep_e + cl[0].E * c1 + cl[1].E * c2);
            double Px = px_e + reco[0].px * c1 + reco[1].px * c2;
            double Py = py_e + reco[0].py * c1 + reco[1].py * c2;
            double Pz = pz_e + reco[0].pz * c1 + reco[1].pz * c2 - e0;
            double mm2 = Etot * Etot - (Px * Px + Py * Py + Pz * Pz);
            double mm = mm2 > 0 ? std::sqrt(mm2) : 0;

            // Debug print for the first 20 events
            static int debug_count = 0;
            if (debug_count < 20) {
                std::cout << "Event " << debug_count+1 << ":\n";
                std::cout << "  Etot = " << Etot << "\n";
                std::cout << "  Px = " << Px << "\n";
                std::cout << "  Py = " << Py << "\n";
                std::cout << "  Pz = " << Pz << "\n";
                std::cout << "  mm2 = " << mm2 << "\n";
                std::cout << "  mm = " << mm << "\n";
                std::cout << "  e0 = " << e0 << ", mp = " << mp << ", ep_e = " << ep_e
                          << ", cl[0].E = " << cl[0].E << ", cl[1].E = " << cl[1].E << "\n";
                std::cout << "  px_e = " << px_e << ", py_e = " << py_e << ", pz_e = " << pz_e << "\n";
                std::cout << "  reco[0].px = " << reco[0].px << ", py = " << reco[0].py << ", pz = " << reco[0].pz << "\n";
                std::cout << "  reco[1].px = " << reco[1].px << ", py = " << reco[1].py << ", pz = " << reco[1].pz << "\n";
                std::cout << "  dtarg = " << dtarg << ", thp0r = " << thp0r << ", hztar = " << hztar << "\n";
                ++debug_count;
            }

            // EXPLANATION:
            // - Cluster momenta are now reconstructed from positions and energies, matching NPA.f.
            // - The missing mass formula now uses these reconstructed momenta.
            // - Electron kinematics are reconstructed as before, matching the physics_angles logic.
            // - This should yield missing mass distributions much closer to NPA.f.

            double wt=w_evt*norm/scale;
            int idx=P==Proc::EXCL?0:P==Proc::DELTA?1:P==Proc::SEMI?2:3;
            if(idx<3){ h_pi[idx]->Fill(mpi,wt); h_mm[idx]->Fill(mm,wt);}
            h_pi[3]->Fill(mpi,wt); h_mm[3]->Fill(mm,wt);
            ++kept;
        }
        std::cout<<path<<" kept "<<kept<<" / "<<N<<"\n";
    }

    /* mask first missing‑mass bin in copies for plotting */
    for(auto h:{h_mm[0],h_mm[1],h_mm[2],h_mm[3]}) h->SetBinContent(1,0);

    /* multipage PDF */
    const char* pdf="pi0_simc_plots.pdf";
    TCanvas open("o","",10,10); open.Print((std::string(pdf)+"[").c_str());

    TCanvas c1("c1","π0 mass",1000,800); c1.Divide(2,2);
    c1.cd(1); h_pi[0]->SetLineColor(kBlue); h_pi[0]->Draw("hist");
    c1.cd(2); h_pi[1]->SetLineColor(kGreen); h_pi[1]->Draw("hist");
    c1.cd(3); h_pi[2]->SetLineColor(kMagenta); h_pi[2]->Draw("hist");
    c1.cd(4); h_pi[3]->SetLineColor(kRed); h_pi[3]->Draw("hist");
    c1.Update(); c1.Print(pdf);

    TCanvas c2("c2","Missing mass",1000,800); c2.Divide(2,2);
    c2.cd(1); h_mm[0]->SetLineColor(kBlue); h_mm[0]->Draw("hist");
    c2.cd(2); h_mm[1]->SetLineColor(kGreen); h_mm[1]->Draw("hist");
    c2.cd(3); h_mm[2]->SetLineColor(kMagenta); h_mm[2]->Draw("hist");
    c2.cd(4); h_mm[3]->SetLineColor(kRed); h_mm[3]->Draw("hist");
    c2.Update(); c2.Print(pdf);

    open.Print((std::string(pdf)+"]").c_str());
    std::cout<<"Wrote "<<pdf<<"\n";

    TFile fout("simc_yield_histos.root","RECREATE");
    for(auto h:h_pi) h->Write();
    for(auto h:h_mm) h->Write();
    fout.Close();
    std::cout<<"Histograms saved to simc_yield_histos.root\n";
    return 0;
}
