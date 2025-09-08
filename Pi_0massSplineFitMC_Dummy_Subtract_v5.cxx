/*************************************************************
  Pi_0massSplineFitMC_Dummy_Subtract_v5.cxx

  – arithmetic in raw counts, global 1/realCharge scale at the end
  – dummy already normalised to counts / µC, then converted back to raw
  – builds Toy-MC-weighted π0-mass and Q² distributions (8-pad canvas)

  MM-only updates previously approved:
    * Missing-mass electron 4-vector from measured H.gtr.{px,py,pz,p}
    * NPS rotation kept (−13.43 deg, cm units)
    * Missing-mass photons: pick top-2 by energy (no opening-angle cut)
    * ToyMC π0-pair loop unchanged (all pairs, 3D geometry)
    * Background/spline/dummy untouched

  NEW (2025-08-12):
    * Toggleable debug (MM_DEBUG env var)
    * Prints first N events and any event with mm^2<=0, with detailed components
    *
 Compile instructrions: g++ -O0 -g -std=c++17 -o Pi_0massSplineFitMC_Dummy_Subtract_v5 Pi_0massSplineFitMC_Dummy_Subtract_v5.cxx alglib_src/*.cpp -I. -Ialglib_src `root-config --cflags --libs` -lTMVA -lRooFit
    Usage:
        ./Pi_0massSplineFitMC_Dummy_Subtract_v3 \
            real.root out.pdf dummy.root realCharge_uC dummyCharge_uC
        example(Pass2_4205): ./Pi_0massSplineFitMC_Dummy_Subtract_v5 /cache/hallc/c-nps/analysis/pass2/replays/production/nps_hms_coin_4205_0_1_-1.root test_v5_2_4205.pdf VolatileROOTfiles/dummy_x58_q51_p5_merged.root 29446.282 254792.874
*************************************************************/

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TFitResult.h>
#include "interpolation.h" // alglib

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <limits>
#include <type_traits>

template <typename S>
inline double EvalShiftGeneric(const S& s, double t) {
    if constexpr (std::is_pointer<S>::value) {
        return s ? s->Eval(t) : 0.0;
    } else if constexpr (std::is_same_v<S, alglib::spline1dinterpolant>) {
        return alglib::spline1dcalc(s, t);
    } else {
        return s.Eval(t);
    }
}


// ───────────────────────────────── smoothing helpers
std::vector<double> normalizePeaks(const std::vector<double> &y, int win, double frac)
{
    int n = y.size();
    std::vector<double> out = y;
    double ymax = *std::max_element(y.begin(), y.end()), thr = frac * ymax;
    std::vector<int> cand;
    for (int i = 1; i < n - 1; ++i)
        if (y[i] > y[i - 1] && y[i] > y[i + 1] && y[i] > thr)
            cand.push_back(i);
    std::vector<int> peaks;
    for (int p : cand)
    {
        if (peaks.empty())
            peaks.push_back(p);
        else
        {
            int last = peaks.back();
            if (p - last < win)
            {
                if (y[p] > y[last])
                    peaks.back() = p;
            }
            else
                peaks.push_back(p);
        }
    }
    if (peaks.empty())
        return out;
    double avg = 0;
    for (int p : peaks)
        avg += y[p];
    avg /= peaks.size();
    int h = win / 2;
    for (int p : peaks)
    {
        double scl = avg / y[p];
        for (int i = std::max(0, p - h); i <= std::min(n - 1, p + h); ++i)
        {
            double w = 1.0 - std::fabs(i - p) / double(win + 1);
            out[i] = y[i] * (1 - w) + (y[i] * scl) * w;
        }
    }
    return out;
}
std::vector<double> normalizeTroughs(const std::vector<double> &y, int win, double tol)
{
    int n = y.size();
    std::vector<double> out = y;
    double ymin = *std::min_element(y.begin(), y.end()),
           ymax = *std::max_element(y.begin(), y.end());
    std::vector<int> cand;
    for (int i = 1; i < n - 1; ++i)
        if (y[i] < y[i - 1] && y[i] < y[i + 1] && y[i] <= ymin + tol * (ymax - ymin))
            cand.push_back(i);
    std::vector<int> tr;
    for (int t : cand)
    {
        if (tr.empty())
            tr.push_back(t);
        else
        {
            int last = tr.back();
            if (t - last < win)
            {
                if (y[t] < y[last])
                    tr.back() = t;
            }
            else
                tr.push_back(t);
        }
    }
    if (tr.empty())
        return out;
    double avg = 0;
    for (int t : tr)
        avg += y[t];
    avg /= tr.size();
    int h = win / 2;
    for (int t : tr)
    {
        double scl = avg / y[t];
        for (int i = std::max(0, t - h); i <= std::min(n - 1, t + h); ++i)
        {
            double w = 1.0 - std::fabs(i - t) / double(win + 1);
            out[i] = y[i] * (1 - w) + (y[i] * scl) * w;
        }
    }
    return out;
}
std::vector<double> regionSmooth(const std::vector<double> &d, int swP, int swT)
{
    int n = d.size();
    std::vector<double> s(n);
    double ymax = *std::max_element(d.begin(), d.end()),
           ymin = *std::min_element(d.begin(), d.end()),
           mid = 0.5 * (ymax + ymin);
    for (int i = 0; i < n; ++i)
    {
        int win = (d[i] > mid) ? swP : swT;
        double sum = 0;
        int cnt = 0;
        for (int j = i - win / 2; j <= i + win / 2; ++j)
            if (j >= 0 && j < n)
            {
                sum += d[j];
                ++cnt;
            }
        s[i] = (cnt ? sum / cnt : d[i]);
    }
    return s;
}

// helper: bind a branch only if it exists
template <typename T>
bool bindIfPresent(TTree *t, const char *name, T *addr)
{
    if (!t)
        return false;
    if (t->GetBranch(name))
    {
        t->SetBranchStatus(name, 1);
        t->SetBranchAddress(name, addr);
        return true;
    }
    return false;
}

// ───────────────────────────────── simple cuts
inline bool passHMSCuts(double edt, double dp, double et, double npe,
                        double th, double ph)
{
    return (edt < 0.1 && std::fabs(dp) <= 8.5 && et > 0.6 && npe > 1.0 &&
            std::fabs(th) <= 0.09 && std::fabs(ph) <= 0.09);
}
inline bool isGoodCluster(double e, double t, double x, double y)
{
    return (e >= 0.6 && t >= 149 && t <= 151 &&
            x > -29.16 && x < 29.16 && y > -35.64 && y < 35.64);
}

// ───────────────────────────────── quick Q² helper
const double theta0_deg = 16.48;                    // HMS central angle (deg)
const double theta0_rad = theta0_deg * M_PI / 180.; // radians
inline double compQ2(double E0, double Ep, double th, double ph)
{
    double th_tot = theta0_rad + th;
    double cosT = std::cos(th_tot) * std::cos(ph); // to O(θ²)
    return 2.0 * E0 * Ep * (1.0 - cosT);
}

// —— NPS rotation (left-side ⇒ negative angle) ———————————
constexpr double NPS_theta_deg = -13.43;                      // deg
constexpr double NPS_theta_rad = NPS_theta_deg * M_PI / 180.; // rad
constexpr double NPS_dist = 407.0;                            // cm

inline void rotateNPS(double xDet, double yDet,
                      double &xHall, double &yHall, double &zHall,
                      double zDet = NPS_dist)
{
    // Passive rotation about +y by NPS_theta_deg (left side ⇒ angle is negative)
    const double th = NPS_theta_rad; // = -20.15° in radians
    const double c = std::cos(th), s = std::sin(th);

    // (xDet, yDet, zDet) treated as a vector in NPS frame → hall frame
    xHall = c * xDet - s * zDet;
    yHall = yDet;
    zHall = s * xDet + c * zDet;
}

// Passive rotations (match Python): rotate axes, not the vector
inline void rotZ_passive(double x, double y, double z, double deg,
                         double &xo, double &yo, double &zo)
{
    double th = deg * M_PI / 180.0, c = std::cos(th), s = std::sin(th);
    xo = c * x + s * y;
    yo = -s * x + c * y;
    zo = z;
}

inline void rotY_passive(double x, double y, double z, double deg,
                         double &xo, double &yo, double &zo)
{
    double th = deg * M_PI / 180.0, c = std::cos(th), s = std::sin(th);
    xo = c * x - s * z;
    yo = y;
    zo = s * x + c * z;
}

// Fill ordered pairs (i != j). Produces symmetric plot like your example.
static inline void FillNpsTimePairs(
    TH2F& H, int ncl, const double* E, const double* T,
    double Emin, double tmin, double tmax)
{
    for (int i = 0; i < ncl; ++i) {
        if (E[i] <= Emin) continue;
        const double ti = T[i];
        if (ti < tmin || ti > tmax) continue;
        for (int j = 0; j < ncl; ++j) {
            if (j == i || E[j] <= Emin) continue;
            const double tj = T[j];
            if (tj < tmin || tj > tmax) continue;
            H.Fill(ti, tj);
        }
    }
}


// ──────────────────────────────────────────── main
int main(int argc, char *argv[])
{
    if (argc < 6)
    {
        std::cerr << "Usage: " << argv[0]
                  << " real.root out.pdf dummy.root realCharge_uC dummyCharge_uC\n";
        return 1;
    }
    std::string realF = argv[1];
    std::string outPDF = argv[2];
    std::string dumF = argv[3];
    double realQ = std::stod(argv[4]); // µC
    double dumQ = std::stod(argv[5]);  // µC

    // Debug toggle (env var)
    const char *dbgEnv = gSystem->Getenv("MM_DEBUG");
    const bool DEBUG_MM = (dbgEnv && std::string(dbgEnv) != "0");
    const int DEBUG_MAX = 20;
    int dbgPrinted = 0;
    auto dbgOk = [&](bool force = false)
    { return DEBUG_MM && (force || dbgPrinted < DEBUG_MAX); };
    auto dbgBump = [&]()
    { ++dbgPrinted; };

    //----------------------------------------------------------------
    // constants (current kinematic = E₀ ~10.54 GeV)
    //----------------------------------------------------------------
    const double e0_nom = 10.54350201; // GeV (beam)
    const double ep0_nom = 5.878;      // GeV (central scattered electron)
    const double bgLo = 113, bgHi = 142.5;
    const double sigLo = 141.789, sigHi = 171.289;
    const double shift_time = 28.05;
    const int nBins = 650;

    // ────────────────────────────────── 1) dummy histogram (counts / µC)
    TH1F *hDumNorm = nullptr;
    TH1F hSub_before("hSub_before", "", nBins, sigLo, sigHi);
    TH1F hSub_after("hSub_after", "", nBins, sigLo, sigHi);
    TH1F hD("hD", "", nBins, sigLo, sigHi);
    TH1F hMissMass("hMissMass", "Missing Mass;M_{miss} [GeV/c^{2}];Counts", 300, 0, 5);

// --- NPS cluster timing pairs: before vs after timing subtraction (real data only) ---
const int    NB_TIM  = 200;
const double TMIN_TS = 140.0, TMAX_TS = 160.0;   // display range

TH2F hPairs_preTS ("hPairs_preTS",
    "NPS Cluster Timing Pairs (Raw; E>0.6 GeV);Cluster i Timing [ns];Cluster j Timing [ns]",
    NB_TIM, TMIN_TS, TMAX_TS, NB_TIM, TMIN_TS, TMAX_TS);

TH2F hPairs_postTS("hPairs_postTS",
    "NPS Cluster Timing Pairs (After timing subtraction; E>0.6 GeV);Cluster i Timing [ns];Cluster j Timing [ns]",
    NB_TIM, TMIN_TS, TMAX_TS, NB_TIM, TMIN_TS, TMAX_TS);

hPairs_preTS .SetDirectory(nullptr);
hPairs_postTS.SetDirectory(nullptr);


    // --- QA opening-angle histograms (no selections; SIDIS-friendly) ---
    TH1F h_open_ang("h_open_ang",
                    "Opening angle #theta_{12};#theta_{12} [rad];Counts",
                    200, 0.0, 0.20);

    TH1F h_open_ang_resid("h_open_ang_resid",
                          "Residual #Delta#theta;#theta_{12}-#theta_{ideal} [rad];Counts",
                          200, -0.05, 0.05);

    TH2F h_theta12_vs_mgg("h_theta12_vs_mgg",
                          "#theta_{12} vs m_{#gamma#gamma};#theta_{12} [rad];m_{#gamma#gamma} [GeV]",
                          200, 0.0, 0.20, 200, 0.0, 0.30);

    TH2F h_resid_vs_asym("h_resid_vs_asym",
                         "#Delta#theta vs A;A=(E_{1}-E_{2})/(E_{1}+E_{2});#theta_{12}-#theta_{ideal} [rad]",
                         200, -1.0, 1.0, 200, -0.05, 0.05);

    {
        TFile fd(dumF.c_str(), "READ");
        if (!fd.IsZombie())
        {
            TTree *t = dynamic_cast<TTree *>(fd.Get("T"));
            if (t)
            {
                double edt = 0, dp = 0, et = 0, npe = 0, th = 0, ph = 0, gy = 0, ct = 0;
                t->SetBranchStatus("*", 0);
                t->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
                t->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edt);
                t->SetBranchStatus("H.gtr.dp", 1);
                t->SetBranchAddress("H.gtr.dp", &dp);
                t->SetBranchStatus("H.cal.etotnorm", 1);
                t->SetBranchAddress("H.cal.etotnorm", &et);
                t->SetBranchStatus("H.cer.npeSum", 1);
                t->SetBranchAddress("H.cer.npeSum", &npe);
                t->SetBranchStatus("H.gtr.th", 1);
                t->SetBranchAddress("H.gtr.th", &th);
                t->SetBranchStatus("H.gtr.ph", 1);
                t->SetBranchAddress("H.gtr.ph", &ph);
                t->SetBranchStatus("H.gtr.y", 1);
                t->SetBranchAddress("H.gtr.y", &gy);
                t->SetBranchStatus("NPS.cal.clusT", 1);
                t->SetBranchAddress("NPS.cal.clusT", &ct);

                TH1F *hUp = new TH1F("dum_up", "", nBins, sigLo, sigHi);
                TH1F *hDn = new TH1F("dum_dn", "", nBins, sigLo, sigHi);
                hUp->SetDirectory(nullptr);
                hDn->SetDirectory(nullptr);

                const Long64_t N = t->GetEntries();
                for (Long64_t i = 0; i < N; ++i)
                {
                    t->GetEntry(i);
                    if (!passHMSCuts(edt, dp, et, npe, th, ph))
                        continue;
                    if (ct < sigLo || ct > sigHi)
                        continue;
                    (gy > 0 ? hDn : hUp)->Fill(ct);
                }
                hUp->Scale(1.0 / (dumQ * 8.467));
                hDn->Scale(1.0 / (dumQ * 4.256));
                hDumNorm = static_cast<TH1F *>(hUp->Clone("hDumNorm"));
                hDumNorm->Add(hDn);
                hDumNorm->SetDirectory(nullptr);
            }
        }
    }

    // ────────────────────────────────── 2) real file (raw)
    TFile fr(realF.c_str(), "READ");
    if (fr.IsZombie())
    {
        std::cerr << "Cannot open " << realF << "\n";
        return 1;
    }
    TTree *tr = dynamic_cast<TTree *>(fr.Get("T"));
    if (!tr)
    {
        std::cerr << "No T tree in " << realF << "\n";
        return 1;
    }

    double edt = 0, dp = 0, et = 0, npe = 0, th = 0, ph = 0;
    tr->SetBranchStatus("*", 0);
    tr->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    tr->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edt);
    tr->SetBranchStatus("H.gtr.dp", 1);
    tr->SetBranchAddress("H.gtr.dp", &dp);
    tr->SetBranchStatus("H.cal.etotnorm", 1);
    tr->SetBranchAddress("H.cal.etotnorm", &et);
    tr->SetBranchStatus("H.cer.npeSum", 1);
    tr->SetBranchAddress("H.cer.npeSum", &npe);
    tr->SetBranchStatus("H.gtr.th", 1);
    tr->SetBranchAddress("H.gtr.th", &th);
    tr->SetBranchStatus("H.gtr.ph", 1);
    tr->SetBranchAddress("H.gtr.ph", &ph);

    // measured electron momentum components for MM
    double hp = 0, hpx = 0, hpy = 0, hpz = 0;
    tr->SetBranchStatus("H.gtr.p", 1);
    tr->SetBranchAddress("H.gtr.p", &hp);
    tr->SetBranchStatus("H.gtr.px", 1);
    tr->SetBranchAddress("H.gtr.px", &hpx);
    tr->SetBranchStatus("H.gtr.py", 1);
    tr->SetBranchAddress("H.gtr.py", &hpy);
    tr->SetBranchStatus("H.gtr.pz", 1);
    tr->SetBranchAddress("H.gtr.pz", &hpz);

    double ncl = 0;
    static const int MAX = 10000;
    double cE[MAX], cT[MAX], cX[MAX], cY[MAX];
    tr->SetBranchStatus("NPS.cal.nclust", 1);
    tr->SetBranchAddress("NPS.cal.nclust", &ncl);
    tr->SetBranchStatus("NPS.cal.clusE", 1);
    tr->SetBranchAddress("NPS.cal.clusE", cE);
    tr->SetBranchStatus("NPS.cal.clusT", 1);
    tr->SetBranchAddress("NPS.cal.clusT", cT);
    tr->SetBranchStatus("NPS.cal.clusX", 1);
    tr->SetBranchAddress("NPS.cal.clusX", cX);
    tr->SetBranchStatus("NPS.cal.clusY", 1);
    tr->SetBranchAddress("NPS.cal.clusY", cY);

    // New: Histogram for all clusT[0] values (full range)
    const double t0Min = 0.0;
    const double t0Max = 200.0;
    const int nRegionBins = 650;
    const double regionWidth = bgHi - bgLo;
    const double binWidth = regionWidth / nRegionBins;
    const int nAllBins = static_cast<int>((t0Max - t0Min) / binWidth + 0.5);

    TH1F hAll("hAll", "All clusT[0] (charge normalized)", nAllBins, t0Min, t0Max);
    hAll.SetDirectory(nullptr);

    Long64_t Nr = tr->GetEntries();
    for (Long64_t i = 0; i < Nr; ++i)
    {
        tr->GetEntry(i);
        if (ncl < 1)
            continue;
        double t0 = cT[0];
        hAll.Fill(t0);
    }
    hAll.Scale(1.0 / realQ);
    {
        TFile outHist("charge_normalized_raw_yield.root", "RECREATE");
        hAll.Write();
    }

    // Second pass: Apply HMS cuts and process for background/signal as before
    TH1F hBG("bg", "", nRegionBins, bgLo, bgHi);
    hBG.SetDirectory(nullptr);
    TH1F hSig("sig", "", nRegionBins, sigLo, sigHi);
    hSig.SetDirectory(nullptr);

    for (Long64_t i = 0; i < Nr; ++i)
    {
        tr->GetEntry(i);
        if (!passHMSCuts(edt, dp, et, npe, th, ph))
            continue;
        if (ncl < 1)
            continue;
        double t0 = cT[0];
        if (t0 >= bgLo && t0 <= bgHi)
            hBG.Fill(t0);
        if (t0 >= sigLo && t0 <= sigHi)
            hSig.Fill(t0);
    }
    hBG.SetBins(nBins, bgLo, bgHi);
    hSig.SetBins(nBins, sigLo, sigHi);
    hSig.SetDirectory(nullptr);


    for (Long64_t i = 0; i < Nr; ++i)
    {
        tr->GetEntry(i);
        if (!passHMSCuts(edt, dp, et, npe, th, ph))
            continue;
        if (ncl < 1)
            continue;
        double t0 = cT[0];
        if (t0 >= bgLo && t0 <= bgHi)
            hBG.Fill(t0);
        if (t0 >= sigLo && t0 <= sigHi)
            hSig.Fill(t0);
    }

    // ────────────────────────────────── 3) spline BG and subtract
    std::vector<double> vx(nBins), vy(nBins);
    for (int b = 1; b <= nBins; ++b)
    {
        vx[b - 1] = hBG.GetBinCenter(b);
        vy[b - 1] = hBG.GetBinContent(b);
    }
    TGraph gBG(nBins, &vx[0], &vy[0]);
    gBG.Sort();
    TSpline3 spl("spl", &gBG);
    std::vector<double> yS(nBins);
    for (int i = 0; i < nBins; ++i)
        yS[i] = spl.Eval(vx[i]);
    auto yP = normalizePeaks(yS, 10, 0.75);
    auto yT = normalizeTroughs(yP, 15, 0.05);
    auto yF = regionSmooth(yT, 2, 5);

    std::vector<double> vxs(nBins), vys(nBins);
    for (int i = 0; i < nBins; ++i)
    {
        vxs[i] = vx[i] + shift_time;
        vys[i] = yF[i];
    }
    TGraph gShift(nBins, &vxs[0], &vys[0]);

    TH1F hSub = *static_cast<TH1F *>(hSig.Clone("hSub"));
    hSub.SetDirectory(nullptr);
    for (int b = 1; b <= nBins; ++b)
    {
        double c = hSub.GetBinCenter(b);
        hSub.SetBinContent(b, hSig.GetBinContent(b) - gShift.Eval(c));
    }

    if (hDumNorm)
    {
        hSub_before = hSub;
        hD = *hDumNorm;
        hD.Scale(realQ); // back to raw
        hSub.Add(&hD, -1.0);
        hSub_after = hSub;
    }
    std::cout << "Integral hSub (raw) = " << hSub.Integral() << "\n";

// ===== Second pass: fill RAW and AFTER–timing–subtraction maps together =====
// Assumes: tr, Nr, passHMSCuts, ncl, cE[], cT[], sigLo/sigHi, hSig, TMIN_TS/TMAX_TS, gShift exist.

Long64_t nPre=0, nPost=0;  // tiny counters for sanity

for (Long64_t i = 0; i < Nr; ++i) {
    tr->GetEntry(i);
    if (!passHMSCuts(edt, dp, et, npe, th, ph)) continue;

    const int N = (int)std::min(ncl, (double)MAX);
    if (N < 2) continue;

    for (int a = 0; a < N; ++a) {
        if (cE[a] <= 0.6) continue;
        const double ti = cT[a];
        if (ti < TMIN_TS || ti > TMAX_TS) continue;

        for (int b = 0; b < N; ++b) {
            if (b == a || cE[b] <= 0.6) continue;
            const double tj = cT[b];
            if (tj < TMIN_TS || tj > TMAX_TS) continue;

            // --- always fill RAW (no timing subtraction)
            hPairs_preTS.Fill(ti, tj); ++nPre;

            // --- timing-subtracted weight, defined only inside [sigLo, sigHi]
            const double tbar = 0.5*(ti + tj);
            if (tbar < sigLo || tbar > sigHi) continue;

            const int    bSig = hSig.FindBin(tbar);
            const double dSig = hSig.GetBinContent(bSig);
            const double bEst = gShift.Eval(tbar);   // if gShift is a pointer, use: gShift ? gShift->Eval(tbar) : 0.0

            double wTS = (dSig > 1e-9) ? (1.0 - bEst / dSig) : 0.0;
            if (wTS < 0.0) wTS = 0.0;
            if (wTS > 1.0) wTS = 1.0;

            hPairs_postTS.Fill(ti, tj, wTS); ++nPost;
        }
    }
}

printf("[Pairs] filled RAW=%lld, POST=%lld\n", (long long)nPre, (long long)nPost);


    // ────────────────────────────────── 4) arrays for Toy MC
    std::vector<double> dataVal(nBins), bgVal(nBins), bgErr(nBins);
    for (int b = 1; b <= nBins; ++b)
    {
        dataVal[b - 1] = hSub.GetBinContent(b);
        double est = gShift.Eval(hSub.GetBinCenter(b));
        bgVal[b - 1] = est;
        bgErr[b - 1] = std::sqrt(std::max(est, 0.0));
    }

    // ────────────────────────────────── 5) build cluster list + event list
    const double DNPS = 407.0; // cm
    struct Cl
    {
        int tbin;
        double m;
    };
    struct Ev
    {
        int tbin;
        double q2;
    };

    std::vector<Cl> cls;
    cls.reserve(1'000'000);
    std::vector<Ev> evtList;
    evtList.reserve(500'000);

    const double me = 0.000511; // GeV
    const double mp = 0.938272; // GeV

    if (DEBUG_MM)
    {
        std::cout << "[MM-DEBUG] Enabled. Will print up to " << DEBUG_MAX
                  << " events and any with mm^2<=0\n";
        std::cout << "[MM-DEBUG] Branch presence: "
                  << " H.gtr.p=" << (tr->GetBranch("H.gtr.p") ? "Y" : "N")
                  << " px=" << (tr->GetBranch("H.gtr.px") ? "Y" : "N")
                  << " py=" << (tr->GetBranch("H.gtr.py") ? "Y" : "N")
                  << " pz=" << (tr->GetBranch("H.gtr.pz") ? "Y" : "N")
                  << "\n";
    }

    for (Long64_t ie = 0; ie < Nr; ++ie)
    {
        tr->GetEntry(ie);
        if (!passHMSCuts(edt, dp, et, npe, th, ph))
            continue;

        // event-by-event Q² (unchanged)
        double ep_evt = ep0_nom * (1.0 + dp / 100.0);
        double q2_evt = compQ2(e0_nom, ep_evt, th, ph);

        if (ncl < 1)
            continue;
        double tfirst = cT[0];
        if (tfirst < sigLo || tfirst > sigHi)
            continue;

        int tbin_evt = hSub.FindBin(tfirst) - 1;
        if (tbin_evt < 0 || tbin_evt >= nBins)
            continue;
        evtList.push_back({tbin_evt, q2_evt});

        // Select good clusters
        std::vector<int> keep;
        keep.reserve((int)ncl);
        for (int c = 0; c < (int)ncl; ++c)
        {
            if (cT[c] < sigLo || cT[c] > sigHi)
                continue;
            if (!isGoodCluster(cE[c], cT[c], cX[c], cY[c]))
                continue;
            keep.push_back(c);
        }
        if (keep.size() < 2)
            continue;

        // --- Missing-mass path ---
        std::nth_element(keep.begin(), keep.begin() + 1, keep.end(),
                         [&](int a, int b)
                         { return cE[a] > cE[b]; });
        int i1 = keep[0], i2 = keep[1];

        // Use slopes-built electron (hall/beam frame) for missing-mass
        /* const double me = 0.000511; // GeV
        double px_s = hp * ph;
        double py_s = hp * th;
        double pz_s = std::sqrt(std::max(0.0, hp*hp - px_s*px_s - py_s*py_s));
        double Ep   = std::sqrt(hp*hp + me*me);  // energy from |p|

        double peprime_px = px_s;
        double peprime_py = py_s;
        double peprime_pz = pz_s; */

        // --- Electron 4-vector (Python-style) ---
        // Assumes these branches are already bound:
        //   double hp  = H.gtr.p;     // magnitude
        //   double hpx = H.gtr.px;    // measured components
        //   double hpy = H.gtr.py;
        //   double hpz = H.gtr.pz;

        const double me = 0.000511; // GeV

        // Use measured components directly (no slope rebuild)
        double peprime_px = hpx;
        double peprime_py = hpy;
        double peprime_pz = hpz;

        // Energy: choose ONE of the following

        // (B) Exact match to the notebook (massless approx):
        double Ep = hp;

        // (A) Slightly more physical (tiny difference): include electron mass
        // double p_mag = std::sqrt(peprime_px*peprime_px + peprime_py*peprime_py + peprime_pz*peprime_pz);
        // double Ep     = std::sqrt(p_mag*p_mag + me*me);

        // Build the TLorentzVector (avoid constructor arg order gotchas)
        TLorentzVector peprime;
        peprime.SetPxPyPzE(peprime_px, peprime_py, peprime_pz, Ep);

        // Beam + target at rest
        double e0 = e0_nom;
        double Ein = e0 + mp;
        double Pin_x = 0.0, Pin_y = 0.0, Pin_z = e0;

        double zDet_evt = NPS_dist;

        // ---------------- Photon 1 (NPS → Hall via passive Ry(NPS_theta_deg)) ----------------
        double px1_h, py1_h, pz1_h; // hall-frame direction (unnormalized)
        rotY_passive(/*xDet=*/cX[i1], /*yDet=*/cY[i1], /*zDet=*/zDet_evt,
                     /*deg =*/NPS_theta_deg, px1_h, py1_h, pz1_h);

        double r1 = std::sqrt(px1_h * px1_h + py1_h * py1_h + pz1_h * pz1_h);
        double u1x = (r1 > 0) ? px1_h / r1 : 0.0;
        double u1y = (r1 > 0) ? py1_h / r1 : 0.0;
        double u1z = (r1 > 0) ? pz1_h / r1 : 0.0;

        double E1 = cE[i1]; // GeV
        double p1x = E1 * u1x;
        double p1y = E1 * u1y;
        double p1z = E1 * u1z;

        // ---------------- Photon 2 (same transform) -----------------------------------------
        double px2_h, py2_h, pz2_h;
        rotY_passive(/*xDet=*/cX[i2], /*yDet=*/cY[i2], /*zDet=*/zDet_evt,
                     /*deg =*/NPS_theta_deg, px2_h, py2_h, pz2_h);

        double r2 = std::sqrt(px2_h * px2_h + py2_h * py2_h + pz2_h * pz2_h);
        double u2x = (r2 > 0) ? px2_h / r2 : 0.0;
        double u2y = (r2 > 0) ? py2_h / r2 : 0.0;
        double u2z = (r2 > 0) ? pz2_h / r2 : 0.0;

        double E2 = cE[i2]; // GeV
        double p2x = E2 * u2x;
        double p2y = E2 * u2y;
        double p2z = E2 * u2z;

        // Outgoing sum
        double E_out = Ep + E1 + E2;
        double px_out = peprime_px + p1x + p2x;
        double py_out = peprime_py + p1y + p2y;
        double pz_out = peprime_pz + p1z + p2z;

        // Missing mass
        double mm2 = std::pow(Ein - E_out, 2) - std::pow(Pin_x - px_out, 2) - std::pow(Pin_y - py_out, 2) - std::pow(Pin_z - pz_out, 2);
        double mm = (mm2 > 0 && std::isfinite(mm2)) ? std::sqrt(mm2) : 0.0;
        hMissMass.Fill(mm);

        // === Open-angle QA (no cuts added) ==================================
        const double cos12_raw = u1x * u2x + u1y * u2y + u1z * u2z;
        const double cos12 = std::max(-1.0, std::min(1.0, cos12_raw));
        const double theta12 = std::acos(cos12);

        // π0-consistent "ideal" opening angle from the measured energies only
        const double mpi0 = 0.1349768; // GeV
        const double s = (E1 > 0 && E2 > 0) ? (mpi0 * mpi0) / (2.0 * E1 * E2) : 0.0;
        // guard numerical corner cases so acos argument is in [-1,1]
        const double arg = 1.0 - std::max(0.0, std::min(2.0, s));
        const double theta_ideal = std::acos(arg);

        // Residual and energy asymmetry
        const double dtheta = theta12 - theta_ideal;
        const double asym = (E1 > 0 || E2 > 0) ? (E1 - E2) / (E1 + E2) : 0.0;

        // Invariant mass (from energies + angle) for the 2D QA
        const double mgg2 = std::max(0.0, 2.0 * E1 * E2 * (1.0 - cos12));
        const double mgg = std::sqrt(mgg2);

        // Fill QA histograms
        h_open_ang.Fill(theta12);
        h_open_ang_resid.Fill(dtheta);
        h_theta12_vs_mgg.Fill(theta12, mgg);
        h_resid_vs_asym.Fill(asym, dtheta);
        // ====================================================================

        // ---- DEBUG PRINTS ----
        bool bad = (!(std::isfinite(peprime_px) && std::isfinite(peprime_py) && std::isfinite(peprime_pz) && std::isfinite(Ep))) || (r1 <= 0 || r2 <= 0) || (!std::isfinite(E1) || !std::isfinite(E2)) || (!std::isfinite(mm2));
        if (dbgOk() || bad || mm2 <= 0)
        {
            std::cout.setf(std::ios::fixed);
            std::cout.precision(6);
            std::cout << "[MM-DEBUG] ie=" << ie
                      << " t0=" << tfirst
                      << " keep=" << keep.size()
                      << " i1=" << i1 << " E1=" << E1
                      << " i2=" << i2 << " E2=" << E2
                      << "\n";

            std::cout << "  rot1(x,y,z)=" << px1_h << "," << py1_h << "," << pz1_h << " r1=" << r1
                      << " u1=" << u1x << "," << u1y << "," << u1z << "\n";
            std::cout << "  rot2(x,y,z)=" << px2_h << "," << py2_h << "," << pz2_h << " r2=" << r2
                      << " u2=" << u2x << "," << u2y << "," << u2z << "\n";
            std::cout << "  sums: Ein=" << Ein << " Eout=" << E_out
                      << " PinZ=" << Pin_z << " poutZ=" << pz_out
                      << " (px,py,pz)_out=" << px_out << "," << py_out << "," << pz_out << "\n";
            std::cout << "  mm2=" << mm2 << "  mm=" << mm
                      << ((mm2 <= 0) ? "  [mm2<=0]" : "")
                      << (bad ? "  [NaN/Inf or bad r]" : "")
                      << "\n";
            dbgBump();
        }

        // ---------- π0 candidate list for Toy-MC (UNCHANGED) ----------
        for (size_t a = 0; a < keep.size(); ++a)
            for (size_t b = a + 1; b < keep.size(); ++b)
            {
                int j1 = keep[a], j2 = keep[b];

                double xra, yra, zra;
                rotateNPS(cX[j1], cY[j1], xra, yra, zra, zDet_evt);
                double xrb, yrb, zrb;
                rotateNPS(cX[j2], cY[j2], xrb, yrb, zrb, zDet_evt);

                double rA = std::sqrt(xra * xra + yra * yra + zra * zra);
                double rB = std::sqrt(xrb * xrb + yrb * yrb + zrb * zrb);
                if (rA <= 0 || rB <= 0)
                    continue;

                double ux1 = xra / rA, uy1 = yra / rA, uz1 = zra / rA;
                double ux2 = xrb / rB, uy2 = yrb / rB, uz2 = zrb / rB;

                double cosg = ux1 * ux2 + uy1 * uy2 + uz1 * uz2;
                cosg = std::max(-1.0, std::min(1.0, cosg));

                double m = std::sqrt(std::max(0.0, 2.0 * cE[j1] * cE[j2] * (1.0 - cosg)));

                double tM = 0.5 * (cT[j1] + cT[j2]);
                int tb = hSub.FindBin(tM) - 1;
                if (tb < 0 || tb >= nBins)
                    continue;
                cls.push_back({tb, m});
            }
    }

    std::cout << "Cluster list size = " << cls.size() << "\n";
    std::cout << "Event   list size = " << evtList.size() << "\n";

    // ────────────────────────────────── 6) Toy MC (UNCHANGED)
    const int Ntoys = 200, nMB = 200;
    const double mLo = 0, mHi = 0.3;
    const int nQ2B = 120;
    const double q2Lo = 0, q2Hi = 12; // GeV²

    TRandom3 rng(0);
    std::vector<TH1F *> toyPi, toyQ2;
    toyPi.reserve(Ntoys);
    toyQ2.reserve(Ntoys);

    for (int it = 0; it < Ntoys; ++it)
    {
        std::vector<double> w(nBins, 0.0);
        for (int i = 0; i < nBins; ++i)
        {
            double dev = rng.Gaus(0.0, bgErr[i]);
            double bt = bgVal[i] + dev;
            if (bt < 0)
                bt = 0;
            double d = dataVal[i];
            if (d < 0)
                d = 0;
            double f = (d > 1e-6) ? ((d - bt) / d) : 0;
            if (f > 1)
                f = 1;
            if (f < -1)
                f = -1;
            w[i] = f;
        }
        TH1F *hM = new TH1F(Form("toyM%d", it), "", nMB, mLo, mHi);
        TH1F *hQ = new TH1F(Form("toyQ%d", it), "", nQ2B, q2Lo, q2Hi);
        for (const auto &c : cls)
            hM->Fill(c.m, w[c.tbin]);
        for (const auto &ev : evtList)
            hQ->Fill(ev.q2, w[ev.tbin]);
        toyPi.push_back(hM);
        toyQ2.push_back(hQ);
    }

    TH1F *hMean = (TH1F *)toyPi[0]->Clone("hMean");
    hMean->Reset();
    TH1F *hVar = (TH1F *)toyPi[0]->Clone("hVar");
    hVar->Reset();
    for (auto h : toyPi)
        hMean->Add(h);
    hMean->Scale(1.0 / Ntoys);
    for (auto h : toyPi)
        for (int b = 1; b <= nMB; ++b)
        {
            double diff = h->GetBinContent(b) - hMean->GetBinContent(b);
            hVar->AddBinContent(b, diff * diff);
        }
    hVar->Scale(1.0 / Ntoys);
    for (int b = 1; b <= nMB; ++b)
        hVar->SetBinContent(b, std::sqrt(hVar->GetBinContent(b)));

    TH1F *hQmean = (TH1F *)toyQ2[0]->Clone("hQmean");
    hQmean->Reset();
    TH1F *hQvar = (TH1F *)toyQ2[0]->Clone("hQvar");
    hQvar->Reset();
    for (auto h : toyQ2)
        hQmean->Add(h);
    hQmean->Scale(1.0 / Ntoys);
    for (auto h : toyQ2)
        for (int b = 1; b <= nQ2B; ++b)
        {
            double diff = h->GetBinContent(b) - hQmean->GetBinContent(b);
            hQvar->AddBinContent(b, diff * diff);
        }
    hQvar->Scale(1.0 / Ntoys);
    for (int b = 1; b <= nQ2B; ++b)
        hQvar->SetBinContent(b, std::sqrt(hQvar->GetBinContent(b)));

    double sf = 1.0 / realQ;
    hMean->Scale(sf);
    hVar->Scale(sf);
    for (int b = 1; b <= hMean->GetNbinsX(); ++b)
    {
        double var = hVar->GetBinContent(b);
        double err = (Ntoys > 1) ? std::sqrt(var / (Ntoys - 1)) : 0.0;
        hMean->SetBinError(b, err);
        TFile fMissMass("missing_mass.root", "RECREATE");
        hMissMass.Write();
    }
    hQmean->Scale(sf);
    hQvar->Scale(sf);

    auto makeBand = [&](TH1F *hC, TH1F *hE) -> TGraph *
    {
        int nb = hC->GetNbinsX();
        std::vector<double> px, py;
        px.reserve(2 * nb);
        py.reserve(2 * nb);
        for (int b = 1; b <= nb; ++b)
        {
            px.push_back(hC->GetBinCenter(b));
            py.push_back(hC->GetBinContent(b));
        }
        for (int b = nb; b >= 1; --b)
        {
            px.push_back(hC->GetBinCenter(b));
            py.push_back(hC->GetBinContent(b) + hE->GetBinContent(b));
        }
        auto *g = new TGraph((int)px.size(), px.data(), py.data());
        g->SetFillColorAlpha(kBlue, 0.35);
        return g;
    };
    TGraph *gBandM = makeBand(hMean, hVar);
    TGraph *gBandQ = makeBand(hQmean, hQvar);

    TCanvas c("c", "", 1200, 4100);
    c.Divide(1, 8);

    // Pad 1: Missing mass
    c.cd(1);
    std::cout << "hMissMass entries: " << hMissMass.GetEntries() << std::endl;

    hMissMass.SetLineColor(kRed + 1);
    hMissMass.SetLineWidth(2);
    hMissMass.SetStats(1);

    // --- Set X range first ---
    const double MM_XMIN = 0.0;
    const double MM_XMAX = 2.5;
    hMissMass.GetXaxis()->SetRangeUser(MM_XMIN, MM_XMAX);

    // --- Compute Y max only over the visible X range ---
    int ibLo = hMissMass.GetXaxis()->FindBin(MM_XMIN + 1e-9);
    int ibHi = hMissMass.GetXaxis()->FindBin(MM_XMAX - 1e-9);
    double maxVis = 0.0;
    for (int ib = ibLo; ib <= ibHi; ++ib)
    {
        double c = hMissMass.GetBinContent(ib);
        if (c > maxVis)
            maxVis = c;
    }
    if (maxVis > 0.0)
    {
        hMissMass.GetYaxis()->SetRangeUser(0.0, maxVis * 1.2);
    }

    // Draw
    hMissMass.Draw("HIST");

    // Pad 2: hAll (charge-normalized raw yield)
    c.cd(2);
    hAll.SetLineColor(kBlack);
    hAll.SetTitle("All clusT[0] (charge normalized);clusT[0];Counts/\\muC");
    hAll.GetXaxis()->SetRangeUser(100, hAll.GetXaxis()->GetXmax());
    hAll.Draw("HIST");

    // Pad 3: Background region with smoothed fit (normalized)
    c.cd(3);
    hBG.SetLineColor(kBlack);
    hBG.SetMarkerStyle(20);
    hBG.SetMarkerSize(0.8);
    hBG.SetTitle("Background region (normalized) with smoothed fit;clusT[0];Counts/\\muC");
    hBG.Draw("P");
    std::vector<double> yF_norm(nBins);
    for (int i = 0; i < nBins; ++i)
        yF_norm[i] = yF[i] / realQ;
    TGraph *grBGfinal_norm = new TGraph(nBins, &vx[0], &yF_norm[0]);
    grBGfinal_norm->SetLineColor(kBlue);
    grBGfinal_norm->SetLineWidth(1);
    grBGfinal_norm->Draw("L SAME");

    // Prepare shifted background fit (normalized)
    std::vector<double> vys_shifted_norm(nBins);
    for (int i = 0; i < nBins; ++i)
        vys_shifted_norm[i] = vys[i] / realQ;
    TGraph *grBGfit_shifted_norm = new TGraph(nBins, &vxs[0], &vys_shifted_norm[0]);

    // Pad 4: Signal vs shifted BG vs subtracted (normalized)
    c.cd(4);
    hSig.SetLineColor(kRed);
    hSig.SetTitle("Signal, Background (fit, shifted), and Subtracted (all normalized);clusT[0];Counts/\\muC");
    hSig.Draw("HIST");
    grBGfit_shifted_norm->SetLineColor(kBlue);
    grBGfit_shifted_norm->SetLineWidth(1);
    grBGfit_shifted_norm->Draw("L SAME");
    hSub.SetLineColor(kGreen);
    hSub.SetLineStyle(1);
    hSub.Draw("HIST SAME");
    {
        auto leg3 = new TLegend(0.6, 0.7, 0.88, 0.88);
        leg3->AddEntry(&hSig, "Signal (hSig)", "l");
        leg3->AddEntry(grBGfit_shifted_norm, "Background fit (shifted, norm)", "l");
        leg3->AddEntry(&hSub_after, "Subtracted (hSub_after)", "l");
        leg3->Draw();
    }

    // Pad 5: ToyMC π0 mass mean
    c.cd(5);
    hMean->SetLineColor(kBlack);
    hMean->SetTitle("ToyMC #pi^{0} mass mean;M_{#gamma#gamma} (GeV);Counts/\\muC");
    hMean->Draw("HIST");

    // Pad 6: ToyMC π0 mass mean + band (+ optional SIMC overlay if present)
    c.cd(6);
    {
        double ymaxM = 1.2 * (hMean->GetMaximum() + hVar->GetMaximum());
        TH2F frame5("f5", ";M_{#gamma#gamma} (GeV);Counts/\\muC", 10, 0, 0.30, 10, 0, std::max(ymaxM, 1e-6));
        frame5.Draw("AXIS");
        TGraph *gBandM2 = makeBand(hMean, hVar);
        gBandM2->Draw("F SAME");
        hMean->Draw("HIST SAME");

        hMean->Fit("gaus", "SAME", "", 0.122, 0.142);
        TF1 *fitFunc = hMean->GetFunction("gaus");
        if (fitFunc)
        {
            fitFunc->SetLineWidth(1);
            fitFunc->Draw("SAME");
            double chi2 = fitFunc->GetChisquare();
            int ndf = fitFunc->GetNDF();
            double mean = fitFunc->GetParameter(1);
            double meanErr = fitFunc->GetParError(1);
            double sigma = fitFunc->GetParameter(2);
            double sigmaErr = fitFunc->GetParError(2);
            double constant = fitFunc->GetParameter(0);
            double constantErr = fitFunc->GetParError(0);
            int bin1 = hMean->FindBin(0.12);
            int bin2 = hMean->FindBin(0.14);
            double hist_integral = hMean->Integral(bin1, bin2);
            TPaveText *pint = new TPaveText(0.6, 0.65, 0.88, 0.7, "NDC");
            pint->SetFillColor(0);
            pint->SetTextAlign(12);
            pint->AddText(Form("Hist Int[0.12,0.14] = %.4f", hist_integral));
            pint->Draw("SAME");
            TPaveText *stats = new TPaveText(0.6, 0.7, 0.88, 0.88, "NDC");
            stats->SetFillColor(0);
            stats->SetTextAlign(12);
            stats->AddText(Form("Entries = %.0f", hMean->GetEntries()));
            stats->AddText(Form("Mean = %.4f #pm %.4f", mean, meanErr));
            stats->AddText(Form("Sigma = %.4f #pm %.4f", sigma, sigmaErr));
            stats->AddText(Form("Const = %.4f #pm %.4f", constant, constantErr));
            stats->AddText(Form("#chi^{2}/NDF = %.2f / %d", chi2, ndf));
            stats->Draw("SAME");
        }
    }

    // Pad 7: Q² mean + error band
    c.cd(7);
    {
        double ymaxQ = 1.2 * (hQmean->GetMaximum() + hQvar->GetMaximum());
        TH2F frame6("f6", ";Q^{2} (GeV^{2});Counts/\\muC", 10, q2Lo, q2Hi, 10, 0, std::max(ymaxQ, 1e-8));
        frame6.Draw("AXIS");
        gBandQ->Draw("F SAME");
        hQmean->Draw("HIST SAME");
    }

    // Pad 8: Dummy subtraction visualization (raw)
    c.cd(8);
    hSub_before.SetLineColor(kBlack);
    hSub_before.SetTitle("Dummy Subtraction Step;clusT[0];Counts/raw");
    hSub_before.SetLineStyle(1);
    hSub_before.Draw("HIST");
    hD.SetLineColor(kOrange + 1);
    hD.SetLineStyle(2);
    hD.Draw("HIST SAME");
    hSub_after.SetLineColor(kRed);
    hSub_after.SetLineStyle(1);
    hSub_after.Draw("HIST SAME");
    {
        auto leg7 = new TLegend(0.6, 0.7, 0.88, 0.88);
        leg7->AddEntry(&hSub_before, "Before subtraction", "l");
        leg7->AddEntry(&hD, "Scaled dummy", "l");
        leg7->AddEntry(&hSub_after, "After subtraction", "l");
        leg7->Draw();
    }

    c.Print(outPDF.c_str());
    std::cout << "All done. Canvas saved to " << outPDF << "\n";

    // --- Optional: quick QA page for opening-angle diagnostics ----------
    {
        // Derive "<out>_qa.pdf" from the main PDF name
        std::string outPDFqa = outPDF;
        auto dot = outPDFqa.find_last_of('.');
        if (dot != std::string::npos)
            outPDFqa.insert(dot, "_qa");
        else
            outPDFqa += "_qa.pdf";

        TCanvas cQA("cQA", "Open-angle QA", 1200, 1200);
        cQA.Divide(2, 2);

        cQA.cd(1);
        h_open_ang.SetLineWidth(2);
        h_open_ang.Draw("HIST");
        cQA.cd(2);
        h_open_ang_resid.SetLineWidth(2);
        h_open_ang_resid.Draw("HIST");
        cQA.cd(3);
        h_theta12_vs_mgg.Draw("COLZ");
        cQA.cd(4);
        h_resid_vs_asym.Draw("COLZ");

        cQA.Print(outPDFqa.c_str());
    }
    // --------------------------------------------------------------------

    // ----- Full-page Missing Mass PDF (rebinned 100 bins, 0.5→2.5, matplotlib-like aspect) -----
    {
        // Derive "<out>_MM.pdf" from your existing outPDF name
        std::string outMM = outPDF;
        const auto dot = outMM.find_last_of('.');
        if (dot != std::string::npos)
            outMM.insert(dot, "_MM");
        else
            outMM += "_MM.pdf";

        // Aspect ratio from your reference image: 775x486 ≈ 1.595:1 (width:height)
        const int W = 1600;              // pick page width; height follows from aspect
        const double AR = 775.0 / 486.0; // ≈ 1.595
        const int H = int(W / AR + 0.5);
        TCanvas cMM("cMM", "Missing Mass (full page)", W, H);
        cMM.cd();

        // --- Rebin/view setup (keeps your current binning; ngroup=1 is fine) ---
        const double xlo = 0.5, xhi = 2.5;
        const int ngroup = 1; // your source hist is already fine-grained

        TH1F *hMM_view = (TH1F *)hMissMass.Rebin(ngroup, "hMM_view");
        hMM_view->SetDirectory(nullptr);

        // Limit display to [0.5, 2.5] GeV
        hMM_view->GetXaxis()->SetRangeUser(xlo, xhi);

        // Auto Y-range only on visible part
        double maxVis = 0.0;
        for (int b = 1; b <= hMM_view->GetNbinsX(); ++b)
        {
            double xc = hMM_view->GetBinCenter(b);
            if (xc >= xlo && xc <= xhi)
                maxVis = std::max(maxVis, (double)hMM_view->GetBinContent(b));
        }
        if (maxVis > 0.0)
            hMM_view->GetYaxis()->SetRangeUser(0.0, maxVis * 1.2);

        // Style + draw
        hMM_view->SetLineColor(kRed + 1);
        hMM_view->SetLineWidth(2);
        hMM_view->SetStats(1);
        hMM_view->GetXaxis()->SetTitle("Missing Mass m_{X}  [GeV]");
        hMM_view->GetYaxis()->SetTitle("Counts");
        hMM_view->Draw("HIST");
        hMM_view->SetStats(0);

        // ==================== Gaussian fit around proton peak ====================

        double fit_lo = 0.70;
        double fit_hi = 1.00;

        // amplitude guess = max bin in window
        int i1 = hMM_view->GetXaxis()->FindBin(fit_lo + 1e-9);
        int i2 = hMM_view->GetXaxis()->FindBin(fit_hi - 1e-9);
        double A0 = 0.0;
        for (int i = i1; i <= i2; ++i)
            A0 = std::max(A0, (double)hMM_view->GetBinContent(i));

        // define & seed
        TF1 fG("fG", "gaus", fit_lo, fit_hi); // [0]=A, [1]=mu, [2]=sigma
        fG.SetParameters(A0, 0.938, 0.020);

        // style (and make a smooth curve)
        fG.SetLineColor(kOrange + 7);
        fG.SetLineWidth(3);
        fG.SetNpx(800);

        // fit in range, store result, but don't auto-draw
        TFitResultPtr fr = hMM_view->Fit(&fG, "RQ0S"); // R=range, Q=quiet, 0=no draw, S=store

        // --- Extend the drawn Gaussian to ±N·sigma (clamped to visible axis) ---
        const double Nsig = 6.0; // 5–7 typical
        const double mu_fit = fG.GetParameter(1);
        const double sigma_fit = fG.GetParameter(2);
        double left = std::max(mu_fit - Nsig * sigma_fit, hMM_view->GetXaxis()->GetXmin());
        double right = std::min(mu_fit + Nsig * sigma_fit, hMM_view->GetXaxis()->GetXmax());

        // reuse the fitted function, just widen its draw range and redraw once
        fG.SetRange(left, right);
        fG.SetNpx(1200);
        fG.Draw("same");

        // optional: markers for m_p and fit window
        gPad->Update();
        const double yTop = gPad->GetUymax();
        auto vline = [&](double x, Color_t col, Style_t sty)
        {
            TLine *L = new TLine(x, 0.0, x, yTop);
            L->SetLineColor(col);
            L->SetLineStyle(sty);
            L->SetLineWidth(3);
            L->Draw("same");
        };
        vline(0.938272, kBlue + 2, 2);
        vline(fit_lo, kAzure + 6, 3);
        vline(fit_hi, kAzure + 6, 3);

        // numbers
        const double A = fG.GetParameter(0);
        const double mu = fG.GetParameter(1);
        const double sigma = fG.GetParameter(2);
        const double eA = fG.GetParError(0);
        const double emu = fG.GetParError(1);
        const double esig = fG.GetParError(2);
        const double chi2 = fG.GetChisquare();
        const int ndf = fG.GetNDF();
        const double yield = std::sqrt(2 * M_PI) * sigma * A; // ∫ Gaussian
        const double yerr = std::sqrt(2 * M_PI) * std::hypot(sigma * eA, A * esig);

        // annotate
        TPaveText box(0.58, 0.68, 0.88, 0.90, "NDC");
        box.SetFillColor(0);
        box.SetFillStyle(0);
        box.SetBorderSize(0);
        box.AddText(Form("Gaussian fit  [%.2f, %.2f] GeV", fit_lo, fit_hi));
        box.AddText(Form("#mu = %.6f #pm %.6f GeV", mu, emu));
        box.AddText(Form("#sigma = %.4f #pm %.4f GeV", sigma, esig));
        box.AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf));
        box.AddText(Form("Yield = %.1f #pm %.1f", yield, yerr));
        box.Draw("same");

        TLegend leg(0.12, 0.80, 0.46, 0.90);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.AddEntry(hMM_view, "Missing mass", "l");
        leg.AddEntry(&fG, "Gaussian fit", "l");
        leg.AddEntry((TObject *)0, "Markers: m_{p}, window", "");
        leg.Draw("same");

        // console
        printf("[MM simple gauss] A=%.3f±%.3f  mu=%.6f±%.6f  sigma=%.4f±%.4f  chi2/ndf=%.2f/%d  Yield=%.1f±%.1f\n",
               A, eA, mu, emu, sigma, esig, chi2, ndf, yield, yerr);

        // ========================================================================

        cMM.Print(outMM.c_str());
        delete hMM_view; // keep memory tidy
    }

// --- Draw timing-pair maps: raw vs after timing subtraction ---
std::string outPairsTS = outPDF;
if (auto dot = outPairsTS.find_last_of('.'); dot != std::string::npos) outPairsTS.insert(dot, "_NPSpairsTime");
else outPairsTS += "_NPSpairsTime.pdf";

gStyle->SetOptStat(0);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,18,0)
gStyle->SetPalette(kViridis);
#else
gStyle->SetPalette(kBird);
#endif

// Match color scales for fair visual comparison
double zmax_ts = std::max(hPairs_preTS.GetMaximum(), hPairs_postTS.GetMaximum());
if (zmax_ts > 0) { hPairs_preTS.SetMaximum(zmax_ts); hPairs_postTS.SetMaximum(zmax_ts); }

TCanvas cPairsTS("cPairsTS","NPS Cluster Timing Pairs (Raw vs Time-Sub)", 1600, 800);
cPairsTS.Divide(2,1);

cPairsTS.cd(1);
gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
hPairs_preTS.GetZaxis()->SetTitle("Counts");
hPairs_preTS.GetZaxis()->SetTitleOffset(1.2);
hPairs_preTS.Draw("COLZ");

cPairsTS.cd(2);
gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
hPairs_postTS.GetZaxis()->SetTitle("Weighted Counts (time-subtracted)");
hPairs_postTS.GetZaxis()->SetTitleOffset(1.2);
hPairs_postTS.Draw("COLZ");

cPairsTS.Print(outPairsTS.c_str());


    // -----------------------------------------------
    gSystem->Exit(0);
    return 0;
}
