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
 Compile instructrions: g++ -O0 -g -std=c++17 -o Pi_0massSplineFitMC_Dummy_Subtract_v5 Pi_0massSplineFitMC_Dummy_Subtract_v5.cxx alglib_src/*.cpp -I. -Ialglib_src `root-config --cflags --libs` -lTMVA -lRooFit -lRooFitCore
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

// --- RooFit ---
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "TROOT.h"

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

// ---- Global-scope helpers (place above any function definitions) ----
#include <algorithm> // if not already included

template <typename W>
static inline double totalWidth(const W* wins, int n) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        const double w = wins[i].hi - wins[i].lo;
        if (w > 0.0) s += w;     // ignore inverted/empty windows
    }
    return s;
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

    // --- constants used across the whole analysis (put near other constants)
const double mp = 0.938272;   // GeV  (proton mass)
const double me = 0.000511;   // GeV  (electron mass, if you need it)

// --- Time-proximity selector parameters (no TOF) ---
const double T_COIN   = 150.0;   // ns, fixed Hall-C coincidence reference
const double SIG_TBAR = 0.80;    // ns, tolerance for t̄ closeness to 150
const double SIG_DTEL = 0.60;    // ns, tolerance for Δt closeness to 0



    // ────────────────────────────────── 1) dummy histogram (counts / µC)
    TH1F *hDumNorm = nullptr;
    TH1F hSub_before("hSub_before", "", nBins, sigLo, sigHi);
    TH1F hSub_after("hSub_after", "", nBins, sigLo, sigHi);
    TH1F hD("hD", "", nBins, sigLo, sigHi);
    TH1F hMissMass("hMissMass", "Missing Mass;M_{miss} [GeV/c^{2}];Counts", 300, 0, 5);

// ===== Absolute timing windows (ns) used by timing-template classification & B_est =====
struct Win { double lo, hi; };

// Central (coincidence) window fixed at 150 ns ±1 ns
const double C_LO = 149.0;
const double C_HI = 151.0;

// Far side windows on BOTH sides; skip the adjacent bunches [147,149] and [151,153]
static const Win posWins[] = { {153.0,155.0}, {155.0,157.0}, {157.0,159.0} };
static const Win negWins[] = { {145.0,147.0}, {143.0,145.0}, {141.0,143.0} };

const int NPOS = (int)(sizeof(posWins)/sizeof(posWins[0]));
const int NNEG = (int)(sizeof(negWins)/sizeof(negWins[0]));

// >>> Add the helpers RIGHT HERE (still inside main) <<<
    auto whichPos = [&](double t)->int {
        for (int k = 0; k < NPOS; ++k)
            if (t >= posWins[k].lo && t <= posWins[k].hi) return k;
        return -1;
    };
    auto whichNeg = [&](double t)->int {
        for (int k = 0; k < NNEG; ++k)
            if (t >= negWins[k].lo && t <= negWins[k].hi) return k;
        return -1;
    };

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

// Ratio map (post / pre) and Δt vs t̄ projections
TH2F hPairs_ratio("hPairs_ratio",
  "NPS Timing Pairs Ratio (Post / Raw);Cluster i Timing [ns];Cluster j Timing [ns]",
  NB_TIM, TMIN_TS, TMAX_TS, NB_TIM, TMIN_TS, TMAX_TS);
hPairs_ratio.SetDirectory(nullptr);
hPairs_ratio.Sumw2(true); // keep errors

const int NB_DT = 200;
TH2F hDT_Tbar_pre("hDT_Tbar_pre",
  "#Delta t vs #bar{t} (Raw);#bar{t} [ns];#Delta t = t_{i}-t_{j} [ns]",
  NB_TIM, TMIN_TS, TMAX_TS, NB_DT, -10.0, 10.0);
TH2F hDT_Tbar_post("hDT_Tbar_post",
  "#Delta t vs #bar{t} (After timing subtraction);#bar{t} [ns];#Delta t = t_{i}-t_{j} [ns]",
  NB_TIM, TMIN_TS, TMAX_TS, NB_DT, -10.0, 10.0);
hDT_Tbar_pre .SetDirectory(nullptr);
hDT_Tbar_post.SetDirectory(nullptr);

// --- t1 vs t2 for the time-proximity selected pair (one per event) ---
TH2F hT12_timePick(
    "hT12_timePick",
    "t_{#gamma1} vs t_{#gamma2} (time-proximity pair);t_{#gamma1} [ns];t_{#gamma2} [ns]",
    hPairs_preTS.GetNbinsX(), hPairs_preTS.GetXaxis()->GetXmin(), hPairs_preTS.GetXaxis()->GetXmax(),
    hPairs_preTS.GetNbinsY(), hPairs_preTS.GetYaxis()->GetXmin(), hPairs_preTS.GetYaxis()->GetXmax()
);
hT12_timePick.SetDirectory(nullptr);


// --- Mass-gated Δt vs t̄ maps (π0 window vs sideband)
const double M_PI0_PDGB = 0.1349766;   // GeV (unused but handy)
double MWIN_LO = 0.115, MWIN_HI = 0.155;   // π0 mass window (adjust as you like)
double MSB_LO  = 0.180, MSB_HI  = 0.240;   // sideband window

TH2F hDT_Tbar_pre_pi0 ("hDT_Tbar_pre_pi0",
  "#Delta t vs #bar{t} (Raw, m_{#gamma#gamma} #in [0.115,0.155] GeV);#bar{t} [ns];#Delta t = t_{j}-t_{i} [ns]",
  NB_TIM, TMIN_TS, TMAX_TS, 200, -10.0, 10.0);

TH2F hDT_Tbar_post_pi0("hDT_Tbar_post_pi0",
  "#Delta t vs #bar{t} (After timing sub, m_{#gamma#gamma} #in [0.115,0.155] GeV);#bar{t} [ns];#Delta t = t_{j}-t_{i} [ns]",
  NB_TIM, TMIN_TS, TMAX_TS, 200, -10.0, 10.0);

TH2F hDT_Tbar_pre_side ("hDT_Tbar_pre_side",
  "#Delta t vs #bar{t} (Raw, sideband m_{#gamma#gamma} #in [0.18,0.24] GeV);#bar{t} [ns];#Delta t = t_{j}-t_{i} [ns]",
  NB_TIM, TMIN_TS, TMAX_TS, 200, -10.0, 10.0);

TH2F hDT_Tbar_post_side("hDT_Tbar_post_side",
  "#Delta t vs #bar{t} (After timing sub, sideband m_{#gamma#gamma} #in [0.18,0.24] GeV);#bar{t} [ns];#Delta t = t_{j}-t_{i} [ns]",
  NB_TIM, TMIN_TS, TMAX_TS, 200, -10.0, 10.0);

hDT_Tbar_pre_pi0 .SetDirectory(nullptr);
hDT_Tbar_post_pi0.SetDirectory(nullptr);
hDT_Tbar_pre_side .SetDirectory(nullptr);
hDT_Tbar_post_side.SetDirectory(nullptr);

// --- Timing-template maps in (t1, t2) for top-2 E clusters (energy+fid only) ---
TH2F hT12_all   ("hT12_all",   "t_{#gamma1} vs t_{#gamma2} (ALL cand);t_{#gamma1} [ns];t_{#gamma2} [ns]",
                 hPairs_preTS.GetNbinsX(), hPairs_preTS.GetXaxis()->GetXmin(), hPairs_preTS.GetXaxis()->GetXmax(),
                 hPairs_preTS.GetNbinsY(), hPairs_preTS.GetYaxis()->GetXmin(), hPairs_preTS.GetYaxis()->GetXmax());
TH2F hT12_CC    ("hT12_CC",    "CC (C#timesC);t_{#gamma1} [ns];t_{#gamma2} [ns]",
                 hT12_all.GetNbinsX(), hT12_all.GetXaxis()->GetXmin(), hT12_all.GetXaxis()->GetXmax(),
                 hT12_all.GetNbinsY(), hT12_all.GetYaxis()->GetXmin(), hT12_all.GetYaxis()->GetXmax());
TH2F hT12_V     ("hT12_V",     "C#timesA (vertical);t_{#gamma1} [ns];t_{#gamma2} [ns]",
                 hT12_all.GetNbinsX(), hT12_all.GetXaxis()->GetXmin(), hT12_all.GetXaxis()->GetXmax(),
                 hT12_all.GetNbinsY(), hT12_all.GetYaxis()->GetXmin(), hT12_all.GetYaxis()->GetXmax());
TH2F hT12_H     ("hT12_H",     "A#timesC (horizontal);t_{#gamma1} [ns];t_{#gamma2} [ns]",
                 hT12_all.GetNbinsX(), hT12_all.GetXaxis()->GetXmin(), hT12_all.GetXaxis()->GetXmax(),
                 hT12_all.GetNbinsY(), hT12_all.GetYaxis()->GetXmin(), hT12_all.GetYaxis()->GetXmax());
TH2F hT12_Dsame ("hT12_Dsame", "A#timesA (same side/index);t_{#gamma1} [ns];t_{#gamma2} [ns]",
                 hT12_all.GetNbinsX(), hT12_all.GetXaxis()->GetXmin(), hT12_all.GetXaxis()->GetXmax(),
                 hT12_all.GetNbinsY(), hT12_all.GetYaxis()->GetXmin(), hT12_all.GetYaxis()->GetXmax());
TH2F hT12_Dopp  ("hT12_Dopp",  "A#timesA (opposite/off-diag);t_{#gamma1} [ns];t_{#gamma2} [ns]",
                 hT12_all.GetNbinsX(), hT12_all.GetXaxis()->GetXmin(), hT12_all.GetXaxis()->GetXmax(),
                 hT12_all.GetNbinsY(), hT12_all.GetYaxis()->GetXmin(), hT12_all.GetYaxis()->GetXmax());

// --- Timing-template maps (ALL pairs, HMS + E>0.6 only) ---
TH2F hT12_allPairs   ("hT12_allPairs",   "t_{#gamma1} vs t_{#gamma2} (ALL pairs);t_{#gamma1} [ns];t_{#gamma2} [ns]",
    hPairs_preTS.GetNbinsX(), hPairs_preTS.GetXaxis()->GetXmin(), hPairs_preTS.GetXaxis()->GetXmax(),
    hPairs_preTS.GetNbinsY(), hPairs_preTS.GetYaxis()->GetXmin(), hPairs_preTS.GetYaxis()->GetXmax());
TH2F hT12_CC_pairs   ("hT12_CC_pairs",   "CC (C#timesC) — ALL pairs; t_{#gamma1} [ns];t_{#gamma2} [ns]",
    hT12_allPairs.GetNbinsX(), hT12_allPairs.GetXaxis()->GetXmin(), hT12_allPairs.GetXaxis()->GetXmax(),
    hT12_allPairs.GetNbinsY(), hT12_allPairs.GetYaxis()->GetXmin(), hT12_allPairs.GetYaxis()->GetXmax());
TH2F hT12_V_pairs    ("hT12_V_pairs",    "C#timesA (vertical) — ALL pairs; t_{#gamma1} [ns];t_{#gamma2} [ns]",
    hT12_allPairs.GetNbinsX(), hT12_allPairs.GetXaxis()->GetXmin(), hT12_allPairs.GetXaxis()->GetXmax(),
    hT12_allPairs.GetNbinsY(), hT12_allPairs.GetYaxis()->GetXmin(), hT12_allPairs.GetYaxis()->GetXmax());
TH2F hT12_H_pairs    ("hT12_H_pairs",    "A#timesC (horizontal) — ALL pairs; t_{#gamma1} [ns];t_{#gamma2} [ns]",
    hT12_allPairs.GetNbinsX(), hT12_allPairs.GetXaxis()->GetXmin(), hT12_allPairs.GetXaxis()->GetXmax(),
    hT12_allPairs.GetNbinsY(), hT12_allPairs.GetYaxis()->GetXmin(), hT12_allPairs.GetYaxis()->GetXmax());
TH2F hT12_Dsame_pairs("hT12_Dsame_pairs","A#timesA (same side/index) — ALL pairs; t_{#gamma1} [ns];t_{#gamma2} [ns]",
    hT12_allPairs.GetNbinsX(), hT12_allPairs.GetXaxis()->GetXmin(), hT12_allPairs.GetXaxis()->GetXmax(),
    hT12_allPairs.GetNbinsY(), hT12_allPairs.GetYaxis()->GetXmin(), hT12_allPairs.GetYaxis()->GetXmax());
TH2F hT12_Dopp_pairs ("hT12_Dopp_pairs", "A#timesA (opposite/off-diag) — ALL pairs; t_{#gamma1} [ns];t_{#gamma2} [ns]",
    hT12_allPairs.GetNbinsX(), hT12_allPairs.GetXaxis()->GetXmin(), hT12_allPairs.GetXaxis()->GetXmax(),
    hT12_allPairs.GetNbinsY(), hT12_allPairs.GetYaxis()->GetXmin(), hT12_allPairs.GetYaxis()->GetXmax());


TH2F hDT_Tbar_pre_good(
    "hDT_Tbar_pre_good",
    "#Delta t vs #bar{t} (GOOD clusters only, pre);#bar{t} [ns];#Delta t [ns]",
    hDT_Tbar_pre.GetNbinsX(),
    hDT_Tbar_pre.GetXaxis()->GetXmin(),
    hDT_Tbar_pre.GetXaxis()->GetXmax(),
    hDT_Tbar_pre.GetNbinsY(),
    hDT_Tbar_pre.GetYaxis()->GetXmin(),
    hDT_Tbar_pre.GetYaxis()->GetXmax());

TH2F hDT_Tbar_post_good(
    "hDT_Tbar_post_good",
    "#Delta t vs #bar{t} (GOOD clusters only, post/weighted);#bar{t} [ns];#Delta t [ns]",
    hDT_Tbar_post.GetNbinsX(),
    hDT_Tbar_post.GetXaxis()->GetXmin(),
    hDT_Tbar_post.GetXaxis()->GetXmax(),
    hDT_Tbar_post.GetNbinsY(),
    hDT_Tbar_post.GetYaxis()->GetXmin(),
    hDT_Tbar_post.GetYaxis()->GetXmax());

// ---------- MM by timing-class regions (use hMissMass binning) ----------
const int    nMM  = hMissMass.GetNbinsX();
const double mmLo = hMissMass.GetXaxis()->GetXmin();
const double mmHi = hMissMass.GetXaxis()->GetXmax();

TH1F hMM_CC        ("hMM_CC",        "MM;M_{X} (GeV);Counts",                 nMM, mmLo, mmHi);
TH1F hMM_vert      ("hMM_vert",      "MM (vertical C#timesA);M_{X} (GeV);",   nMM, mmLo, mmHi);
TH1F hMM_horiz     ("hMM_horiz",     "MM (horizontal A#timesC);M_{X} (GeV);", nMM, mmLo, mmHi);
TH1F hMM_diag      ("hMM_diag",      "MM (diagonal A#timesA same);M_{X};",    nMM, mmLo, mmHi);
TH1F hMM_pure      ("hMM_pure",      "MM (pure random A#timesA opp.);M_{X};", nMM, mmLo, mmHi);
TH1F hMM_subtracted("hMM_subtracted","MM (CC - B_{est});M_{X} (GeV);Counts",  nMM, mmLo, mmHi);

// ---------- MM by timing-class regions (ALL PAIRS; HMS + E>0.6 only) ----------
const int    nMM_pairs  = hMissMass.GetNbinsX();
const double mmLo_pairs = hMissMass.GetXaxis()->GetXmin();
const double mmHi_pairs = hMissMass.GetXaxis()->GetXmax();

TH1F hMM_CC_pairs        ("hMM_CC_pairs",        "MM (ALL pairs) CC;M_{X} (GeV);Counts",                 nMM_pairs, mmLo_pairs, mmHi_pairs);
TH1F hMM_vert_pairs      ("hMM_vert_pairs",      "MM (ALL pairs) C#timesA;M_{X} (GeV);Counts",           nMM_pairs, mmLo_pairs, mmHi_pairs);
TH1F hMM_horiz_pairs     ("hMM_horiz_pairs",     "MM (ALL pairs) A#timesC;M_{X} (GeV);Counts",           nMM_pairs, mmLo_pairs, mmHi_pairs);
TH1F hMM_diag_pairs      ("hMM_diag_pairs",      "MM (ALL pairs) A#timesA same;M_{X} (GeV);Counts",      nMM_pairs, mmLo_pairs, mmHi_pairs);
TH1F hMM_pure_pairs      ("hMM_pure_pairs",      "MM (ALL pairs) A#timesA opposite;M_{X} (GeV);Counts",  nMM_pairs, mmLo_pairs, mmHi_pairs);
TH1F hMM_subtracted_pairs("hMM_subtracted_pairs","MM (ALL pairs) CC - B_{est};M_{X} (GeV);Counts",       nMM_pairs, mmLo_pairs, mmHi_pairs);

// --- Missing mass using "time-proximity" chosen pair (one pair per event) ---
TH1F hMM_timePick("hMM_timePick",
                  "MM — time-proximity pair;M_{X} (GeV);Counts",
                  hMissMass.GetNbinsX(),
                  hMissMass.GetXaxis()->GetXmin(),
                  hMissMass.GetXaxis()->GetXmax());

// ===== Invariant mass (mgg) timing-template histograms: TOP-2 (or time-prox) =====
TH1F hMgg_CC   ("hMgg_CC",   "m_{#gamma#gamma} CC; m_{#gamma#gamma} (GeV);Counts",             160, 0.00, 0.40);
TH1F hMgg_vert ("hMgg_vert", "m_{#gamma#gamma} vertical; m_{#gamma#gamma} (GeV);Counts",       160, 0.00, 0.40);
TH1F hMgg_horiz("hMgg_horiz","m_{#gamma#gamma} horizontal; m_{#gamma#gamma} (GeV);Counts",     160, 0.00, 0.40);
TH1F hMgg_diag ("hMgg_diag", "m_{#gamma#gamma} diagonal; m_{#gamma#gamma} (GeV);Counts",       160, 0.00, 0.40);
TH1F hMgg_pure ("hMgg_pure", "m_{#gamma#gamma} pure randoms; m_{#gamma#gamma} (GeV);Counts",   160, 0.00, 0.40);
TH1F hMgg_Best ("hMgg_Best", "m_{#gamma#gamma} background estimate; m_{#gamma#gamma} (GeV);Counts", 160, 0.00, 0.40);
TH1F hMgg_sub  ("hMgg_sub",  "m_{#gamma#gamma} CC - B_{est}; m_{#gamma#gamma} (GeV);Counts",   160, 0.00, 0.40);
hMgg_CC.SetDirectory(nullptr); hMgg_vert.SetDirectory(nullptr); hMgg_horiz.SetDirectory(nullptr);
hMgg_diag.SetDirectory(nullptr); hMgg_pure.SetDirectory(nullptr); hMgg_Best.SetDirectory(nullptr);
hMgg_sub.SetDirectory(nullptr);

// ===== Invariant mass (mgg) timing-template histograms: ALL-PAIRS =====
TH1F hMgg_CC_pairs   ("hMgg_CC_pairs",   "m_{#gamma#gamma} CC (ALL pairs); m_{#gamma#gamma} (GeV);Counts",           160, 0.00, 0.40);
TH1F hMgg_vert_pairs ("hMgg_vert_pairs", "m_{#gamma#gamma} vertical (ALL pairs); m_{#gamma#gamma} (GeV);Counts",     160, 0.00, 0.40);
TH1F hMgg_horiz_pairs("hMgg_horiz_pairs","m_{#gamma#gamma} horizontal (ALL pairs); m_{#gamma#gamma} (GeV);Counts",   160, 0.00, 0.40);
TH1F hMgg_diag_pairs ("hMgg_diag_pairs", "m_{#gamma#gamma} diagonal (ALL pairs); m_{#gamma#gamma} (GeV);Counts",     160, 0.00, 0.40);
TH1F hMgg_pure_pairs ("hMgg_pure_pairs", "m_{#gamma#gamma} pure (ALL pairs); m_{#gamma#gamma} (GeV);Counts",         160, 0.00, 0.40);
TH1F hMgg_Best_pairs ("hMgg_Best_pairs", "m_{#gamma#gamma} B_{est} (ALL pairs); m_{#gamma#gamma} (GeV);Counts",      160, 0.00, 0.40);
TH1F hMgg_sub_pairs  ("hMgg_sub_pairs",  "m_{#gamma#gamma} CC - B_{est} (ALL pairs); m_{#gamma#gamma} (GeV);Counts", 160, 0.00, 0.40);
hMgg_CC_pairs.SetDirectory(nullptr); hMgg_vert_pairs.SetDirectory(nullptr); hMgg_horiz_pairs.SetDirectory(nullptr);
hMgg_diag_pairs.SetDirectory(nullptr); hMgg_pure_pairs.SetDirectory(nullptr); hMgg_Best_pairs.SetDirectory(nullptr);
hMgg_sub_pairs.SetDirectory(nullptr);

// -------------------- Missing-mass histograms for three methods --------------------
const int    NB_MM   = 200;           // adjust if you like
const double MMLO    = 0;         // GeV (match your existing MM range)
const double MMHI    = 3.5;         // GeV

TH1D* hMM_top2_sig  = new TH1D("hMM_top2_sig",  "MM (Top-2) signal window", NB_MM, MMLO, MMHI);
TH1D* hMM_top2_acc  = new TH1D("hMM_top2_acc",  "MM (Top-2) accidentals",   NB_MM, MMLO, MMHI);
TH1D* hMM_all_sig   = new TH1D("hMM_all_sig",   "MM (All pairs) signal",    NB_MM, MMLO, MMHI);
TH1D* hMM_all_acc   = new TH1D("hMM_all_acc",   "MM (All pairs) accidentals",NB_MM,MMLO,MMHI);
TH1D* hMM_temp_sig  = new TH1D("hMM_temp_sig",  "MM (Template) signal",     NB_MM, MMLO, MMHI);
TH1D* hMM_temp_acc  = new TH1D("hMM_temp_acc",  "MM (Template) accidentals",NB_MM,MMLO,MMHI);

// Proper error treatment for subtractions
hMM_top2_sig->Sumw2();  hMM_top2_acc->Sumw2();
hMM_all_sig ->Sumw2();  hMM_all_acc ->Sumw2();
hMM_temp_sig->Sumw2();  hMM_temp_acc->Sumw2();

// Consistent colors for the final overlay
hMM_top2_sig->SetLineColor(kRed+1);
hMM_all_sig ->SetLineColor(kBlue+1);
hMM_temp_sig->SetLineColor(kGreen+2);

// --- Absolute-time helpers (use the same posWins/negWins & C_LO/C_HI you defined)
// -------------------- Timing helpers (absolute times) --------------------
auto anyPos = [&](double t)->bool { return whichPos(t) >= 0; };
auto anyNeg = [&](double t)->bool { return whichNeg(t) >= 0; };

auto inCentral = [&](double ti, double tj)->bool {
    return (ti >= C_LO && ti <= C_HI) && (tj >= C_LO && tj <= C_HI);
};

// Accidentals if EITHER time is in any side window (OR-of-two selection)
auto inFar = [&](double ti, double tj)->bool {
    return anyPos(ti) || anyNeg(ti) || anyPos(tj) || anyNeg(tj);
};

// Small utilities for later (unchanged)
auto totalWidth = [](const Win* arr, int N)->double {
    double s = 0.0; for (int i=0; i<N; ++i) s += (arr[i].hi - arr[i].lo); return s;
};

// Call these with (mm, ti, tj)
auto FillMM_Top2 = [&](double mm, double ti, double tj) {
    if (inCentral(ti,tj))      hMM_top2_sig->Fill(mm);
    else if (inFar(ti,tj))     hMM_top2_acc->Fill(mm);
};
auto FillMM_All  = [&](double mm, double ti, double tj) {
    if (inCentral(ti,tj))      hMM_all_sig->Fill(mm);
    else if (inFar(ti,tj))     hMM_all_acc->Fill(mm);
};
auto FillMM_Temp = [&](double mm, double ti, double tj) {
    if (inCentral(ti,tj))      hMM_temp_sig->Fill(mm);
    else if (inFar(ti,tj))     hMM_temp_acc->Fill(mm);
};

// --- Option B: per-event sideband histos (best-in-class per event) ---
TH1 *hMM_vert_evBest  = (TH1*)hMM_vert .Clone("hMM_vert_evBest");  hMM_vert_evBest ->SetDirectory(nullptr); hMM_vert_evBest ->Reset("ICESM"); hMM_vert_evBest ->Sumw2();
TH1 *hMM_horiz_evBest = (TH1*)hMM_horiz.Clone("hMM_horiz_evBest"); hMM_horiz_evBest->SetDirectory(nullptr); hMM_horiz_evBest->Reset("ICESM"); hMM_horiz_evBest->Sumw2();
TH1 *hMM_diag_evBest  = (TH1*)hMM_diag .Clone("hMM_diag_evBest");  hMM_diag_evBest ->SetDirectory(nullptr); hMM_diag_evBest ->Reset("ICESM"); hMM_diag_evBest ->Sumw2();
TH1 *hMM_pure_evBest  = (TH1*)hMM_pure .Clone("hMM_pure_evBest");  hMM_pure_evBest ->SetDirectory(nullptr); hMM_pure_evBest ->Reset("ICESM"); hMM_pure_evBest ->Sumw2();

// --- Option B (union style): per-event histos ---
TH1 *hMM_side_evBest = (TH1*)hMM_CC.Clone("hMM_side_evBest"); // V∪H (exactly one stripe entry per event)
hMM_side_evBest->SetDirectory(nullptr); hMM_side_evBest->Reset("ICESM"); hMM_side_evBest->Sumw2();

TH1 *hMM_AA_evBest   = (TH1*)hMM_CC.Clone("hMM_AA_evBest");   // A×A corners (exactly one corner entry per event)
hMM_AA_evBest  ->SetDirectory(nullptr); hMM_AA_evBest  ->Reset("ICESM"); hMM_AA_evBest  ->Sumw2();

// 2D maps that show EXACTLY what feeds each template background
TH2F hPairs_usedA(
    "hPairs_usedA",
    "Pairs used in Option A background (V + H + A#timesA);t_{i} [ns];t_{j} [ns]",
    NB_TIM, TMIN_TS, TMAX_TS, NB_TIM, TMIN_TS, TMAX_TS
);
hPairs_usedA.SetDirectory(nullptr);
hPairs_usedA.Sumw2();

TH2F hPairs_usedB(
    "hPairs_usedB",
    "Pairs used in Option B background (1 stripe + 1 corner per event);t_{i} [ns];t_{j} [ns]",
    NB_TIM, TMIN_TS, TMAX_TS, NB_TIM, TMIN_TS, TMAX_TS
);
hPairs_usedB.SetDirectory(nullptr);
hPairs_usedB.Sumw2();

// --- Template A/B pointers (declare ONCE early in main) ---
TH1* hTplA_raw = nullptr;
TH1* hTplA_bkg = nullptr;
TH1* hTplA_sub = nullptr;

TH1* hTplB_raw = nullptr;
TH1* hTplB_bkg = nullptr;
TH1* hTplB_sub = nullptr;


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
// Assumes: tr, Nr, passHMSCuts, ncl, cE[], cT[], sigLo/sigHi, hSig,
//          TMIN_TS/TMAX_TS, MAX, gShift, and the histos exist.

Long64_t nPre = 0, nPost = 0;

for (Long64_t i = 0; i < Nr; ++i) {
    tr->GetEntry(i);
    if (!passHMSCuts(edt, dp, et, npe, th, ph)) continue;

    const int N = (int)std::min(ncl, (double)MAX);
    if (N < 2) continue;

// -- per-event "best by time" bookkeeping
int    best_a   = -1, best_b = -1;
double bestS    = 1e99;
double best_ti  = 0.0, best_tj = 0.0;

// --- Option B trackers (reset each event) ---
int    bV_a=-1, bV_b=-1; double bV_S=1e300, bV_mm=0, bV_ti=0, bV_tj=0;
int    bH_a=-1, bH_b=-1; double bH_S=1e300, bH_mm=0, bH_ti=0, bH_tj=0;
int    bD_a=-1, bD_b=-1; double bD_S=1e300, bD_mm=0, bD_ti=0, bD_tj=0;
int    bP_a=-1, bP_b=-1; double bP_S=1e300, bP_mm=0, bP_ti=0, bP_tj=0;

// scoring parameters (same spirit as your time-proximity diagnostic)
const double T0        = 150.0;        // your coincidence center (matches your debug prints ~150.xx)
const double SIG_TBAR  = 1.0;          // scale for (tbar - T0)
const double SIG_DELT  = 1.0;          // scale for Δt
auto score_time = [&](double ti, double tj)->double {
  const double tbar = 0.5*(ti+tj);
  const double dtt  = tj - ti;
  return ((tbar - T0)/SIG_TBAR)*((tbar - T0)/SIG_TBAR) + (dtt/SIG_DELT)*(dtt/SIG_DELT);
};


    for (int a = 0; a < N; ++a) {
        if (cE[a] <= 0.6) continue;
        const double ti = cT[a];
        if (ti < TMIN_TS || ti > TMAX_TS) continue;

        for (int b = 0; b < N; ++b) {
            if (b == a || cE[b] <= 0.6) continue;
            const double tj = cT[b];
            if (tj < TMIN_TS || tj > TMAX_TS) continue;

            // --- Timing coordinates once per pair
            const double tbar = 0.5 * (ti + tj);
            const double dt   = tj - ti;  // Δt = t_j - t_i

// -- time-proximity score: prefer t̄≈150 ns and Δt≈0 (no TOF)
const double z1 = (tbar - T_COIN) / SIG_TBAR;
const double z2 = (dt)            / SIG_DTEL;
const double S_time = z1*z1 + z2*z2;

if (S_time < bestS) {
    bestS   = S_time;
    best_a  = a;
    best_b  = b;
    best_ti = ti;
    best_tj = tj;
}
            // --- Always fill RAW timing maps
            hPairs_preTS.Fill(ti, tj); ++nPre;
            hDT_Tbar_pre.Fill(tbar, dt);

// --- ALL-PAIRS timing-template MM fills (HMS cuts already passed; only E>0.6 used)
// classify (ti,tj) into timing windows (absolute, both sides)
auto whichPos = [&](double t)->int { for (int k=0;k<NPOS;++k) if (t>=posWins[k].lo && t<=posWins[k].hi) return k; return -1; };
auto whichNeg = [&](double t)->int { for (int k=0;k<NNEG;++k) if (t>=negWins[k].lo && t<=negWins[k].hi) return k; return -1; };

const bool C1 = (ti>=C_LO && ti<=C_HI);
const bool C2 = (tj>=C_LO && tj<=C_HI);
const int  p1 = whichPos(ti), n1 = whichNeg(ti);
const int  p2 = whichPos(tj), n2 = whichNeg(tj);
const bool A1 = (p1>=0 || n1>=0);
const bool A2 = (p2>=0 || n2>=0);

const bool CC    = (C1 && C2);
const bool VERT  = (C1 && A2);
const bool HORIZ = (A1 && C2);
const bool sameIdx = ( (p1>=0 && p2>=0 && p1==p2) || (n1>=0 && n2>=0 && n1==n2) );
const bool AAny   = (A1 && A2);
const bool DIAG   = (AAny && sameIdx);
const bool PURE   = (AAny && !sameIdx);

// build MM for THIS ordered pair (a,b) exactly like your selected-pair path
double px1_h, py1_h, pz1_h, px2_h, py2_h, pz2_h;
rotY_passive(cX[a], cY[a], NPS_dist, NPS_theta_deg, px1_h, py1_h, pz1_h);
rotY_passive(cX[b], cY[b], NPS_dist, NPS_theta_deg, px2_h, py2_h, pz2_h);

const double r1 = std::sqrt(px1_h*px1_h + py1_h*py1_h + pz1_h*pz1_h);
const double r2 = std::sqrt(px2_h*px2_h + py2_h*py2_h + pz2_h*pz2_h);
if (r1>0.0 && r2>0.0) {
    const double u1x = px1_h/r1, u1y = py1_h/r1, u1z = pz1_h/r1;
    const double u2x = px2_h/r2, u2y = py2_h/r2, u2z = pz2_h/r2;
    const double E1 = cE[a],     E2 = cE[b];

    // electron & beam (same as elsewhere)
    const double peprime_px = hpx;
    const double peprime_py = hpy;
    const double peprime_pz = hpz;
    const double Ep         = hp;
    const double e0         = e0_nom;
    const double Ein        = e0 + mp;
    const double Pin_x = 0.0, Pin_y = 0.0, Pin_z = e0;

    // photon 4-momenta
    const double p1x = E1*u1x, p1y = E1*u1y, p1z = E1*u1z;
    const double p2x = E2*u2x, p2y = E2*u2y, p2z = E2*u2z;

    // missing mass
    const double E_out  = Ep + E1 + E2;
    const double px_out = peprime_px + p1x + p2x;
    const double py_out = peprime_py + p1y + p2y;
    const double pz_out = peprime_pz + p1z + p2z;

    const double mm2 = std::pow(Ein - E_out, 2)
                     - std::pow(Pin_x - px_out, 2)
                     - std::pow(Pin_y - py_out, 2)
                     - std::pow(Pin_z - pz_out, 2);
    if (mm2>0 && std::isfinite(mm2)) {
        const double mm = std::sqrt(mm2);

        // fill class spectra (ordered pairs, to match RAW multiplicity)
        if (CC)    hMM_CC_pairs   .Fill(mm);
        if (VERT)  hMM_vert_pairs .Fill(mm);
        if (HORIZ) hMM_horiz_pairs.Fill(mm);
        if (DIAG)  hMM_diag_pairs .Fill(mm);
        if (PURE)  hMM_pure_pairs .Fill(mm);
        FillMM_All (mm, ti, tj);
        //FillMM_Temp(mm, ti, tj);
        // --- Option B: fill per-event best sideband entries (require physical mm2>0 like elsewhere)

     // -------- Option B (union): fill at most ONE stripe and ONE corner per event --------
const bool hasV = (bV_a >= 0 && bV_mm > 0);
const bool hasH = (bH_a >= 0 && bH_mm > 0);
const bool hasD = (bD_a >= 0 && bD_mm > 0);
const bool hasP = (bP_a >= 0 && bP_mm > 0);

// --- Option B (union): choose exactly ONE stripe (V or H) and ONE corner (D or P)
// Declare outside the if-scopes so they are visible below.
bool   pickV = false;  // true ⇒ use best vertical (C×A); false ⇒ use best horizontal (A×C)
bool   pickD = false;  // true ⇒ use best diagonal (A×A same-index); false ⇒ use best off-diag (A×A opposite)
double selStripe_ti = std::numeric_limits<double>::quiet_NaN();
double selStripe_tj = std::numeric_limits<double>::quiet_NaN();
double selCorner_ti = std::numeric_limits<double>::quiet_NaN();
double selCorner_tj = std::numeric_limits<double>::quiet_NaN();

// Choose best stripe (V vs H) by the time-score; fill exactly one if any exist
if (hasV || hasH) {
  pickV = (hasV && (!hasH || bV_S <= bH_S));
  const double mmStripe = pickV ? bV_mm : bH_mm;
  hMM_side_evBest->Fill(mmStripe);
  if (pickV && hasV) { selStripe_ti = bV_ti; selStripe_tj = bV_tj; }
  else               { selStripe_ti = bH_ti; selStripe_tj = bH_tj; }
}

// Choose best corner (D vs P); fill exactly one if any exist
if (hasD || hasP) {
  pickD = (hasD && (!hasP || bD_S <= bP_S));
  const double mmCorner = pickD ? bD_mm : bP_mm;
  hMM_AA_evBest->Fill(mmCorner);
  if (pickD && hasD) { selCorner_ti = bD_ti; selCorner_tj = bD_tj; }
  else               { selCorner_ti = bP_ti; selCorner_tj = bP_tj; }
}

// Record the ACTUAL pairs chosen by Option B (one stripe + one corner per event)
if (std::isfinite(selStripe_ti) && std::isfinite(selStripe_tj)) {
  hPairs_usedB.Fill(selStripe_ti, selStripe_tj);
}
if (std::isfinite(selCorner_ti) && std::isfinite(selCorner_tj)) {
  hPairs_usedB.Fill(selCorner_ti, selCorner_tj);
}

// Reset trackers for next event (keep whatever reset block you already use)        
// -------- Option B: UPDATE per-event best sideband candidates (no filling here) --------
  // requires you already declared & reset per-event trackers earlier in the event:
  //   bV_a,bV_b,bV_S,bV_mm,bV_ti,bV_tj, etc., and the score_time(ti,tj) helper
  const double S_this = score_time(ti, tj);

  if (VERT) {
    if (S_this < bV_S) { bV_S = S_this; bV_a = a; bV_b = b; bV_mm = mm; bV_ti = ti; bV_tj = tj; }
  }
  if (HORIZ) {
    if (S_this < bH_S) { bH_S = S_this; bH_a = a; bH_b = b; bH_mm = mm; bH_ti = ti; bH_tj = tj; }
  }
  if (DIAG) {
    if (S_this < bD_S) { bD_S = S_this; bD_a = a; bD_b = b; bD_mm = mm; bD_ti = ti; bD_tj = tj; }
  }
  if (PURE) {
    if (S_this < bP_S) { bP_S = S_this; bP_a = a; bP_b = b; bP_mm = mm; bP_ti = ti; bP_tj = tj; }
  }


    }
}

// fill all-pairs timing-template maps (ordered pairs, to match your RAW map multiplicity)
hT12_allPairs.Fill(ti, tj);
if (CC)    hT12_CC_pairs   .Fill(ti, tj);
if (VERT)  hT12_V_pairs    .Fill(ti, tj);
if (HORIZ) hT12_H_pairs    .Fill(ti, tj);
if (DIAG)  hT12_Dsame_pairs.Fill(ti, tj);
if (PURE)  hT12_Dopp_pairs .Fill(ti, tj);

            // --- GOOD clusters (per-cluster timing/energy/fiducial gate)
const bool goodA = isGoodCluster(cE[a], cT[a], cX[a], cY[a]);
const bool goodB = isGoodCluster(cE[b], cT[b], cX[b], cY[b]);
if (goodA && goodB) {
    // Pre (unweighted) fill for good clusters only
    hDT_Tbar_pre_good.Fill(tbar, dt);
}


            // ------------------------------------------------------------------
            // Di-photon mass using your existing convention:
            // rotate NPS (cX,cY,NPS_dist) to Hall frame with rotY_passive(...),
            // normalize to unit vectors, then mgg = sqrt[ 2 E1 E2 (1 - cos12) ].
            // ------------------------------------------------------------------
            double pxA_h, pyA_h, pzA_h;
            double pxB_h, pyB_h, pzB_h;

            // Use NPS_dist here (avoid zDet_evt scope issues in this block)
            rotY_passive(/*xDet=*/cX[a], /*yDet=*/cY[a], /*zDet=*/NPS_dist,
                         /*deg=*/NPS_theta_deg, pxA_h, pyA_h, pzA_h);
            rotY_passive(/*xDet=*/cX[b], /*yDet=*/cY[b], /*zDet=*/NPS_dist,
                         /*deg=*/NPS_theta_deg, pxB_h, pyB_h, pzB_h);

            const double rA = std::sqrt(pxA_h*pxA_h + pyA_h*pyA_h + pzA_h*pzA_h);
            const double rB = std::sqrt(pxB_h*pxB_h + pyB_h*pyB_h + pzB_h*pzB_h);
            if (rA <= 0.0 || rB <= 0.0) continue; // guard bad geometry

            const double uAx = pxA_h / rA, uAy = pyA_h / rA, uAz = pzA_h / rA;
            const double uBx = pxB_h / rB, uBy = pyB_h / rB, uBz = pzB_h / rB;

            const double cos12_raw = uAx*uBx + uAy*uBy + uAz*uBz;
            const double cos12     = std::max(-1.0, std::min(1.0, cos12_raw));

            const double E1 = cE[a];
            const double E2 = cE[b];

            const double mgg2 = std::max(0.0, 2.0 * E1 * E2 * (1.0 - cos12));
            const double mgg  = std::sqrt(mgg2);

// --- fill mgg ALL-pairs, using the same class flags as MM pairs
if (mgg > 0.0) {
    if (CC)      hMgg_CC_pairs.Fill(mgg);
    else if (VERT)  hMgg_vert_pairs.Fill(mgg);
    else if (HORIZ) hMgg_horiz_pairs.Fill(mgg);
    else if (DIAG)  hMgg_diag_pairs.Fill(mgg);
    else if (PURE)  hMgg_pure_pairs.Fill(mgg);
}
if (bV_a>=0 && bV_mm>0) hMM_vert_evBest ->Fill(bV_mm);
if (bH_a>=0 && bH_mm>0) hMM_horiz_evBest->Fill(bH_mm);
if (bD_a>=0 && bD_mm>0) hMM_diag_evBest ->Fill(bD_mm);
if (bP_a>=0 && bP_mm>0) hMM_pure_evBest ->Fill(bP_mm);
            // --- RAW mass-gated Δt–t̄ maps (π0 window & sidebands)
            if (mgg > 0.0) {
                if (mgg >= MWIN_LO && mgg <= MWIN_HI) {
                    hDT_Tbar_pre_pi0.Fill(tbar, dt);
                } else if (mgg >= MSB_LO && mgg <= MSB_HI) {
                    hDT_Tbar_pre_side.Fill(tbar, dt);
                }
            }

            // --- Timing-subtracted maps (weighted), only inside the signal t̄ window
            if (tbar >= sigLo && tbar <= sigHi) {
                const int    bSig = hSig.FindBin(tbar);
                const double dSig = hSig.GetBinContent(bSig);

                // Background estimate at t̄ (gShift is your fitted background model)
                const double bEst = gShift.Eval(tbar);

                // Weight = (1 − BG/Sig), clamped to [0,1]
                double wTS = (dSig > 1e-9) ? (1.0 - bEst / dSig) : 0.0;
                if (wTS < 0.0) wTS = 0.0;
                if (wTS > 1.0) wTS = 1.0;

                // POST timing maps (weighted)
                hPairs_postTS.Fill(ti, tj, wTS); ++nPost;
                hDT_Tbar_post.Fill(tbar, dt, wTS);

   //  add this: post (weighted) fill for good clusters only
    if (goodA && goodB) {
        hDT_Tbar_post_good.Fill(tbar, dt, wTS);
    }


                // POST mass-gated maps (weighted)
                if (mgg > 0.0) {
                    if (mgg >= MWIN_LO && mgg <= MWIN_HI) {
                        hDT_Tbar_post_pi0.Fill(tbar, dt, wTS);
                    } else if (mgg >= MSB_LO && mgg <= MSB_HI) {
                        hDT_Tbar_post_side.Fill(tbar, dt, wTS);
                    }
                }
            } // end if tbar in [sigLo, sigHi]
        } // end for b
    } // end for a

// -- after scanning all pairs in this event: fill MM for the best-by-time pair
if (best_a >= 0 && best_b >= 0) {
   
   // -- fill the 2D time map for the chosen pair
    hT12_timePick.Fill(best_ti, best_tj);

    // build photon directions (same recipe you already use elsewhere)
    double px1_h, py1_h, pz1_h, px2_h, py2_h, pz2_h;
    rotY_passive(cX[best_a], cY[best_a], NPS_dist, NPS_theta_deg, px1_h, py1_h, pz1_h);
    rotY_passive(cX[best_b], cY[best_b], NPS_dist, NPS_theta_deg, px2_h, py2_h, pz2_h);

    const double r1 = std::sqrt(px1_h*px1_h + py1_h*py1_h + pz1_h*pz1_h);
    const double r2 = std::sqrt(px2_h*px2_h + py2_h*py2_h + pz2_h*pz2_h);
    if (r1 > 0.0 && r2 > 0.0) {
        const double u1x = px1_h/r1, u1y = py1_h/r1, u1z = pz1_h/r1;
        const double u2x = px2_h/r2, u2y = py2_h/r2, u2z = pz2_h/r2;
        const double E1 = cE[best_a], E2 = cE[best_b];

        // electron & beam (as in your main MM calc)
        const double peprime_px = hpx, peprime_py = hpy, peprime_pz = hpz;
        const double Ep = hp;
        const double e0 = e0_nom, Ein = e0 + mp;
        const double Pin_x = 0.0, Pin_y = 0.0, Pin_z = e0;

        const double p1x = E1*u1x, p1y = E1*u1y, p1z = E1*u1z;
        const double p2x = E2*u2x, p2y = E2*u2y, p2z = E2*u2z;

        const double E_out  = Ep + E1 + E2;
        const double px_out = peprime_px + p1x + p2x;
        const double py_out = peprime_py + p1y + p2y;
        const double pz_out = peprime_pz + p1z + p2z;

        const double mm2 = std::pow(Ein - E_out, 2)
                         - std::pow(Pin_x - px_out, 2)
                         - std::pow(Pin_y - py_out, 2)
                         - std::pow(Pin_z - pz_out, 2);
        if (mm2 > 0 && std::isfinite(mm2)) {
            const double mm_pair = std::sqrt(mm2);
            hMM_timePick.Fill(std::sqrt(mm2));
            FillMM_Top2(mm_pair, best_ti, best_tj);  
        } 
 // --- compute mgg for the chosen pair and fill TOP-2/time-prox mass hists
const double cos12_raw_tp = u1x*u2x + u1y*u2y + u1z*u2z;
const double cos12_tp     = std::max(-1.0, std::min(1.0, cos12_raw_tp));
const double mgg2_tp      = std::max(0.0, 2.0 * E1 * E2 * (1.0 - cos12_tp));
const double mgg_tp       = std::sqrt(mgg2_tp);

// classify the chosen pair by timing windows (reuse the same lambdas / C_LO/C_HI)
const bool C1_tp = (best_ti >= C_LO && best_ti <= C_HI);
const bool C2_tp = (best_tj >= C_LO && best_tj <= C_HI);
const int  p1_tp = whichPos(best_ti), n1_tp = whichNeg(best_ti);
const int  p2_tp = whichPos(best_tj), n2_tp = whichNeg(best_tj);
const bool A1_tp = (p1_tp>=0 || n1_tp>=0);
const bool A2_tp = (p2_tp>=0 || n2_tp>=0);
const bool CC_tp    = (C1_tp && C2_tp);
const bool VERT_tp  = (C1_tp && A2_tp);
const bool HORIZ_tp = (A1_tp && C2_tp);
const bool sameIdx_tp = ( (p1_tp>=0 && p2_tp>=0 && p1_tp==p2_tp) || (n1_tp>=0 && n2_tp>=0 && n1_tp==n2_tp) );
const bool AAny_tp  = (A1_tp && A2_tp);
const bool DIAG_tp  = (AAny_tp && sameIdx_tp);
const bool PURE_tp  = (AAny_tp && !sameIdx_tp);

if (mgg_tp > 0.0) {
    if (CC_tp)         hMgg_CC.Fill(mgg_tp);
    else if (VERT_tp)  hMgg_vert.Fill(mgg_tp);
    else if (HORIZ_tp) hMgg_horiz.Fill(mgg_tp);
    else if (DIAG_tp)  hMgg_diag.Fill(mgg_tp);
    else if (PURE_tp)  hMgg_pure.Fill(mgg_tp);
}


    }

}
} // end for i


printf("[Pairs] filled RAW=%lld, POST=%lld\n",
       (long long)nPre, (long long)nPost);

// --- Build ratio = post / pre, guarding zeros ---
hPairs_ratio.Reset();
for (int ix = 1; ix <= hPairs_preTS.GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hPairs_preTS.GetNbinsY(); ++iy) {
        const double pre  = hPairs_preTS.GetBinContent(ix, iy);
        const double post = hPairs_postTS.GetBinContent(ix, iy);
        hPairs_ratio.SetBinContent(ix, iy, (pre > 0.0 ? post/pre : 0.0));
    }
}


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

// ---- timing-template classification (absolute times, both sides; no timing pre-cuts) ----
// Central (C) = [149,151] ns. Exclude adjacent [147,149] and [151,153].
// Use far side windows on BOTH sides of 150 ns:
//   positives: centers 154, 156, 158 → [153,155], [155,157], [157,159]
//   negatives: centers 146, 144, 142 → [145,147], [143,145], [141,143]
// Keep only energy + fiducial for candidate choice here (no timing cut), so side windows populate.

// Build candidate list: energy + fiducial only (NO timing cut here)
std::vector<int> cand;
cand.reserve((int)ncl);
for (int c = 0; c < (int)ncl; ++c) {
    if (cE[c] < 0.6) continue;                                    // energy
    if (!(cX[c] > -29.16 && cX[c] < 29.16 && cY[c] > -35.64 && cY[c] < 35.64)) continue; // fiducial
    cand.push_back(c);
}
if (cand.size() >= 2) {
    // pick top-2 energy clusters (unbiased by timing)
    std::nth_element(cand.begin(), cand.begin()+1, cand.end(),
                     [&](int a, int b){ return cE[a] > cE[b]; });
    const int j1 = cand[0];
    const int j2 = cand[1];

    const double t1 = cT[j1];
    const double t2 = cT[j2];

    // Central membership
    const bool C1 = (t1 >= C_LO && t1 <= C_HI);
    const bool C2 = (t2 >= C_LO && t2 <= C_HI);

    // Side-window membership: which side (+/-) and which index (0..)
    auto whichPos = [&](double t)->int {
        for (int k = 0; k < NPOS; ++k) if (t >= posWins[k].lo && t <= posWins[k].hi) return k;
        return -1;
    };
    auto whichNeg = [&](double t)->int {
        for (int k = 0; k < NNEG; ++k) if (t >= negWins[k].lo && t <= negWins[k].hi) return k;
        return -1;
    };
    const int p1 = whichPos(t1), n1 = whichNeg(t1);
    const int p2 = whichPos(t2), n2 = whichNeg(t2);

    const bool A1 = (p1 >= 0 || n1 >= 0);
    const bool A2 = (p2 >= 0 || n2 >= 0);

    // Timing classes
    const bool CC    = (C1 && C2);
    const bool VERT  = (C1 && A2);                // C × A
    const bool HORIZ = (A1 && C2);                // A × C

    // A×A split: “same” = both in the SAME far window (same side & same index)
    const bool A_sameSideSameIdx =
        (p1 >= 0 && p2 >= 0 && p1 == p2) || (n1 >= 0 && n2 >= 0 && n1 == n2);
    const bool A_anyFar = (A1 && A2);
    const bool DIAG = (A_anyFar && A_sameSideSameIdx);
    const bool PURE = (A_anyFar && !A_sameSideSameIdx);  // all off-diagonals (opp. sides or diff. indices)

    // ---- Build missing mass for this (j1,j2) pair (same recipe you use for i1,i2) ----
    // Electron 4-vector (massless approx like your main path)
    double peprime_px = hpx;
    double peprime_py = hpy;
    double peprime_pz = hpz;
    double Ep         = hp;

    // Beam + target at rest
    double e0   = e0_nom;
    double Ein  = e0 + mp;
    double Pin_x = 0.0, Pin_y = 0.0, Pin_z = e0;

    // Photons j1, j2: rotate to Hall frame, unit vectors, 4-momenta
    double zDet_evt = NPS_dist;

    double px1_h, py1_h, pz1_h;
    rotY_passive(cX[j1], cY[j1], zDet_evt, NPS_theta_deg, px1_h, py1_h, pz1_h);
    double r1 = std::sqrt(px1_h*px1_h + py1_h*py1_h + pz1_h*pz1_h);
    double u1x = (r1>0)? px1_h/r1 : 0.0, u1y = (r1>0)? py1_h/r1 : 0.0, u1z = (r1>0)? pz1_h/r1 : 0.0;
    double E1 = cE[j1];
    double p1x = E1*u1x, p1y = E1*u1y, p1z = E1*u1z;

    double px2_h, py2_h, pz2_h;
    rotY_passive(cX[j2], cY[j2], zDet_evt, NPS_theta_deg, px2_h, py2_h, pz2_h);
    double r2 = std::sqrt(px2_h*px2_h + py2_h*py2_h + pz2_h*pz2_h);
    double u2x = (r2>0)? px2_h/r2 : 0.0, u2y = (r2>0)? py2_h/r2 : 0.0, u2z = (r2>0)? pz2_h/r2 : 0.0;
    double E2 = cE[j2];
    double p2x = E2*u2x, p2y = E2*u2y, p2z = E2*u2z;

    // Outgoing sum
    double E_out  = Ep + E1 + E2;
    double px_out = peprime_px + p1x + p2x;
    double py_out = peprime_py + p1y + p2y;
    double pz_out = peprime_pz + p1z + p2z;

    // Missing mass
    double mm2 = std::pow(Ein - E_out, 2)
               - std::pow(Pin_x - px_out, 2)
               - std::pow(Pin_y - py_out, 2)
               - std::pow(Pin_z - pz_out, 2);
    double mm = (mm2 > 0 && std::isfinite(mm2)) ? std::sqrt(mm2) : 0.0;

// Fill the overview map for all candidates (top-2 E, energy+fid only)
hT12_all.Fill(t1, t2);

// Fill per-class maps based on your flags
if (CC)    hT12_CC   .Fill(t1, t2);
if (VERT)  hT12_V    .Fill(t1, t2);
if (HORIZ) hT12_H    .Fill(t1, t2);
if (DIAG)  hT12_Dsame.Fill(t1, t2);
if (PURE)  hT12_Dopp .Fill(t1, t2);



    // Fill timing-template spectra (unweighted)
    if (CC)    hMM_CC   .Fill(mm);
    if (VERT)  hMM_vert .Fill(mm);
    if (HORIZ) hMM_horiz.Fill(mm);
    if (DIAG)  hMM_diag .Fill(mm);
    if (PURE)  hMM_pure .Fill(mm);
}


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

// --- Draw ratio page ---
std::string outPairsRatio = outPDF;
if (auto dot = outPairsRatio.find_last_of('.'); dot != std::string::npos) outPairsRatio.insert(dot, "_NPSpairsTime_ratio");
else outPairsRatio += "_NPSpairsTime_ratio.pdf";

gStyle->SetOptStat(0);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,18,0)
gStyle->SetPalette(kViridis);
#else
gStyle->SetPalette(kBird);
#endif

TCanvas cPairsRatio("cPairsRatio","NPS Timing Pairs Ratio", 900, 850);
cPairsRatio.cd();
gPad->SetRightMargin(0.16);
gPad->SetLeftMargin(0.12);
gPad->SetBottomMargin(0.12);

hPairs_ratio.GetZaxis()->SetTitle("Post / Raw");
hPairs_ratio.GetZaxis()->SetRangeUser(0.0, 1.0);  // clamp to [0,1] for interpretability
hPairs_ratio.GetZaxis()->SetTitleOffset(1.2);
hPairs_ratio.Draw("COLZ");

cPairsRatio.Print(outPairsRatio.c_str());

std::string outDT = outPDF;
if (auto dot = outDT.find_last_of('.'); dot != std::string::npos) outDT.insert(dot, "_DTtbar");
else outDT += "_DTtbar.pdf";

gStyle->SetOptStat(0);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,18,0)
gStyle->SetPalette(kViridis);
#else
gStyle->SetPalette(kBird);
#endif

TCanvas cDT("cDT","#Delta t vs #bar{t}", 1600, 800);
cDT.Divide(2,1);

cDT.cd(1);
gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
hDT_Tbar_pre.GetZaxis()->SetTitle("Counts");
hDT_Tbar_pre.GetZaxis()->SetTitleOffset(1.2);
hDT_Tbar_pre.Draw("COLZ");

cDT.cd(2);
gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
hDT_Tbar_post.GetZaxis()->SetTitle("Weighted Counts (time-subtracted)");
hDT_Tbar_post.GetZaxis()->SetTitleOffset(1.2);
hDT_Tbar_post.Draw("COLZ");

cDT.Print(outDT.c_str());


// --- π0 window page
{
  std::string outPi0 = outPDF;
  if (auto dot = outPi0.find_last_of('.'); dot != std::string::npos) outPi0.insert(dot, "_DTtbar_pi0");
  else outPi0 += "_DTtbar_pi0.pdf";

  gStyle->SetOptStat(0);
  #if ROOT_VERSION_CODE >= ROOT_VERSION(6,18,0)
    gStyle->SetPalette(kViridis);
  #else
    gStyle->SetPalette(kBird);
  #endif

  TCanvas cPi0("cPi0","#Delta t vs #bar{t} (π^{0} mass window)", 1600, 800);
  cPi0.Divide(2,1);

  cPi0.cd(1);
  gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
  hDT_Tbar_pre_pi0.GetZaxis()->SetTitle("Counts");
  hDT_Tbar_pre_pi0.GetZaxis()->SetTitleOffset(1.2);
  hDT_Tbar_pre_pi0.Draw("COLZ");

  cPi0.cd(2);
  gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
  hDT_Tbar_post_pi0.GetZaxis()->SetTitle("Weighted Counts (time-subtracted)");
  hDT_Tbar_post_pi0.GetZaxis()->SetTitleOffset(1.2);
  hDT_Tbar_post_pi0.Draw("COLZ");

  cPi0.Print(outPi0.c_str());
}

// --- sideband page
{
  std::string outSB = outPDF;
  if (auto dot = outSB.find_last_of('.'); dot != std::string::npos) outSB.insert(dot, "_DTtbar_side");
  else outSB += "_DTtbar_side.pdf";

  TCanvas cSB("cSB","#Delta t vs #bar{t} (sideband)", 1600, 800);
  cSB.Divide(2,1);

  cSB.cd(1);
  gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
  hDT_Tbar_pre_side.GetZaxis()->SetTitle("Counts");
  hDT_Tbar_pre_side.GetZaxis()->SetTitleOffset(1.2);
  hDT_Tbar_pre_side.Draw("COLZ");

  cSB.cd(2);
  gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
  hDT_Tbar_post_side.GetZaxis()->SetTitle("Weighted Counts (time-subtracted)");
  hDT_Tbar_post_side.GetZaxis()->SetTitleOffset(1.2);
  hDT_Tbar_post_side.Draw("COLZ");

  cSB.Print(outSB.c_str());
}

// ---------- Δt vs t̄ (GOOD clusters only) ----------
std::string outDTgood = outPDF;
if (auto dot = outDTgood.find_last_of('.'); dot != std::string::npos)
    outDTgood.insert(dot, "_DTtbar_good");
else
    outDTgood += "_DTtbar_good.pdf";

gStyle->SetOptStat(0);
gStyle->SetPalette(kViridis);

TCanvas cDTgood("cDTgood", "DT vs tbar (GOOD clusters only)", 1200, 600);
cDTgood.Divide(2,1);

// Left: pre
cDTgood.cd(1);
gPad->SetRightMargin(0.16);
gPad->SetLeftMargin(0.12);
gPad->SetBottomMargin(0.12);
hDT_Tbar_pre_good.GetZaxis()->SetTitle("Counts");
hDT_Tbar_pre_good.GetZaxis()->SetTitleOffset(1.2);
hDT_Tbar_pre_good.Draw("COLZ");

// Right: post (weighted)
cDTgood.cd(2);
gPad->SetRightMargin(0.16);
gPad->SetLeftMargin(0.12);
gPad->SetBottomMargin(0.12);
hDT_Tbar_post_good.GetZaxis()->SetTitle("Weighted Counts (time-subtracted)");
hDT_Tbar_post_good.GetZaxis()->SetTitleOffset(1.2);
hDT_Tbar_post_good.Draw("COLZ");

cDTgood.Print(outDTgood.c_str());

// ---------- Build B_est with area normalization (absolute windows, both sides) ----------
const double Wc = (C_HI - C_LO);

// Sum widths of all far windows (both sides)
double sumPos = 0.0, sumNeg = 0.0, sumSq = 0.0;
for (int k = 0; k < NPOS; ++k) { double w = posWins[k].hi - posWins[k].lo; sumPos += w; sumSq += w*w; }
for (int k = 0; k < NNEG; ++k) { double w = negWins[k].hi - negWins[k].lo; sumNeg += w; sumSq += w*w; }

const double sumAll = sumPos + sumNeg;

// Geometric “areas” of the timing classes mapped to the CC box
const double ACC = Wc * Wc;                       // CC area
const double AV  = Wc * sumAll;                   // vertical (C×A, both sides, all far wins)
const double AH  = Wc * sumAll;                   // horizontal (A×C, both sides)
const double AD  = sumSq;                         // A×A same (sum of squares of each far win)
const double AP  = (sumAll * sumAll) - AD;        // A×A off-diagonals (all other A×A combos)

// Coefficients to predict CC-sized contamination from each class
const double aV = (AV > 0) ? (ACC / AV) : 0.0;
const double aH = (AH > 0) ? (ACC / AH) : 0.0;
const double aD = (AD > 0) ? (ACC / AD) : 0.0;
const double aP = (AP > 0) ? (ACC / AP) : 0.0;

TH1F hMM_Best("hMM_Best","Estimated accidental background;M_{X} (GeV);Counts",
              hMM_CC.GetNbinsX(), hMM_CC.GetXaxis()->GetXmin(), hMM_CC.GetXaxis()->GetXmax());
hMM_Best.Reset();
hMM_Best.Add(&hMM_diag,   aD);    // out-of-time real π0
hMM_Best.Add(&hMM_vert,   aV);    // single-γ accidentals (C×A)
hMM_Best.Add(&hMM_horiz,  aH);    // single-γ accidentals (A×C)
hMM_Best.Add(&hMM_pure,  -aP);    // remove uncorrelated off-diagonals

hMM_subtracted.Reset();
hMM_subtracted.Add(&hMM_CC, 1.0);
hMM_subtracted.Add(&hMM_Best, -1.0);

// After hMM_CC (RAW), hMM_Best (A-method B_est), hMM_subtracted are fully built:
hTplA_raw = (TH1*)&hMM_CC;         // RAW = central (per-event)
hTplA_bkg = (TH1*)&hMM_Best;       // B_est from timing-template (A method)
hTplA_sub = (TH1*)&hMM_subtracted; // SUB = RAW - B_est


// ---------- Plot the timing-template background and subtraction ----------
std::string outMM_templ = outPDF;
if (auto dot = outMM_templ.find_last_of('.'); dot != std::string::npos)
    outMM_templ.insert(dot, "_MM_timingTemplate");
else
    outMM_templ += "_MM_timingTemplate.pdf";

TCanvas cMMtempl("cMMtempl","MM timing-template subtraction",1200,600);
cMMtempl.Divide(2,1);

// Left: CC vs B_est overlay
cMMtempl.cd(1);
gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
hMM_CC.SetLineColor(kBlack);
hMM_CC.SetLineWidth(2);
hMM_CC.SetTitle("CC (central) and B_{est}");
hMM_CC.Draw("HIST");
hMM_Best.SetLineColor(kRed+1);
hMM_Best.SetLineStyle(2);
hMM_Best.SetLineWidth(2);
hMM_Best.Draw("HIST SAME");
auto leg1 = new TLegend(0.60,0.72,0.88,0.88);
leg1->AddEntry(&hMM_CC,  "CC data", "l");
leg1->AddEntry(&hMM_Best,"B_{est} (timing template)", "l");
leg1->Draw();

// Right: subtracted CC - B_est
cMMtempl.cd(2);
gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
hMM_subtracted.SetLineColor(kBlue+1);
hMM_subtracted.SetLineWidth(2);
hMM_subtracted.SetTitle("CC - B_{est}");
hMM_subtracted.Draw("HIST");

cMMtempl.Print(outMM_templ.c_str());

// ===== Build mgg background estimate from TOP-2 classes and subtract
hMgg_Best.Reset();
hMgg_Best.Add(&hMgg_vert,  1.0/6.0);
hMgg_Best.Add(&hMgg_horiz, 1.0/6.0);
hMgg_Best.Add(&hMgg_diag,  1.0/6.0);
hMgg_Best.Add(&hMgg_pure, -1.0/30.0);

hMgg_sub.Reset();
hMgg_sub.Add(&hMgg_CC, 1.0);
hMgg_sub.Add(&hMgg_Best, -1.0);

// ===== Build mgg background estimate for ALL-PAIRS and subtract
hMgg_Best_pairs.Reset();
hMgg_Best_pairs.Add(&hMgg_vert_pairs,  1.0/6.0);
hMgg_Best_pairs.Add(&hMgg_horiz_pairs, 1.0/6.0);
hMgg_Best_pairs.Add(&hMgg_diag_pairs,  1.0/6.0);
hMgg_Best_pairs.Add(&hMgg_pure_pairs, -1.0/30.0);

hMgg_sub_pairs.Reset();
hMgg_sub_pairs.Add(&hMgg_CC_pairs, 1.0);
hMgg_sub_pairs.Add(&hMgg_Best_pairs, -1.0);


// ===== RooFit: Gaussian peak + template background (TOP-2) =====
{
    using namespace RooFit;

    // Observable with the same range as your MM hists
    RooRealVar mx("mx","M_{X} (GeV)", hMM_CC.GetXaxis()->GetXmin(), hMM_CC.GetXaxis()->GetXmax());

    // Data: CC spectrum (no subtraction)
    RooDataHist dh_data("dh_data","CC data", RooArgList(mx), &hMM_CC);

    // Background template: B_est (shape only; free overall yield)
    RooDataHist dh_b("dh_b","B_{est} template", RooArgList(mx), &hMM_Best);
    RooHistPdf  bpdf("bpdf","bkg pdf", RooArgSet(mx), dh_b);

    // Signal: a Gaussian around proton mass
    RooRealVar mean("mean","mean", 0.94, 0.80, 1.15);   // center floats; start near m_p
    RooRealVar sigma("sigma","sigma", 0.10, 0.02, 0.30); // width floats
    RooGaussian sig("sig","signal", mx, mean, sigma);

    // Extended yields
    RooRealVar Ns("Ns","Ns", hMM_CC.Integral()*0.5, 0.0, 1e9);
    RooRealVar Nb("Nb","Nb", hMM_CC.Integral()*0.5, 0.0, 1e9);

    // Full model
    RooAddPdf model("model","sig+bkg", RooArgList(sig, bpdf), RooArgList(Ns, Nb));

    // Restrict fit to a sensible window around the peak (keeps SIDIS high-mass out)
    mx.setRange("pk", 0.60, 1.40);

    // Fit
    std::unique_ptr<RooFitResult> fr(model.fitTo(dh_data, Save(true), Extended(true), Range("pk")));

    // Plot
    TCanvas cFit("cFit","MM fit (TOP-2)", 900, 700);
   RooPlot* frame = mx.frame(RooFit::Title("MM: CC data = Gaussian + B_{est} (TOP-2)"));
dh_data.plotOn(frame, RooFit::CutRange("pk"));
model.plotOn(frame, RooFit::Range("pk"));
model.plotOn(frame, RooFit::Components(bpdf), RooFit::LineStyle(kDashed), RooFit::LineColor(kBlue), RooFit::Range("pk"));
model.plotOn(frame, RooFit::Components(sig),  RooFit::LineStyle(kSolid),  RooFit::LineColor(kRed),  RooFit::Range("pk"));
    frame->GetXaxis()->SetTitle("M_{X} (GeV)");
    frame->GetYaxis()->SetTitle("Counts");
    frame->Draw();

    // Print to PDF (new page)
    std::string outFitTop2 = outPDF;
    if (auto dot = outFitTop2.find_last_of('.'); dot != std::string::npos)
        outFitTop2.insert(dot, "_MM_fit_template");
    else
        outFitTop2 += "_MM_fit_template.pdf";
    cFit.Print(outFitTop2.c_str());

    // Quick numbers to stdout
    std::cout << "[MM fit TOP-2] mean=" << mean.getVal() << " ± " << mean.getError()
              << "  sigma=" << sigma.getVal() << " ± " << sigma.getError()
              << "  Ns=" << Ns.getVal() << "  Nb=" << Nb.getVal() << "\n";
}


// ===== TOP-2 timing-template: integrated summary + f_acc(MX) plot =====
{
    auto HInt = [](const TH1F &h)->double { return h.Integral(1, h.GetNbinsX()); };

    const double I_CC  = HInt(hMM_CC);
    const double I_V   = HInt(hMM_vert);
    const double I_H   = HInt(hMM_horiz);
    const double I_D   = HInt(hMM_diag);
    const double I_P   = HInt(hMM_pure);
    const double I_B   = HInt(hMM_Best);
    const double I_Sub = HInt(hMM_subtracted);
    const double f_acc = (I_CC > 0.0) ? (I_B / I_CC) : 0.0;

    printf("[TimingTemplate TOP-2]  CC=%g  V=%g  H=%g  D=%g  P=%g  |  B_est=%g  Sub=%g  f_acc=%g\n",
           I_CC, I_V, I_H, I_D, I_P, I_B, I_Sub, f_acc);

    // Per-bin accidental fraction: f_acc(MX) = B_est / CC  (clamped to [0,1])
    TH1F hMM_frac_top2("hMM_frac_top2",
                       "Accidental fraction f_{acc}(M_{X}) (TOP-2);M_{X} (GeV);f_{acc}",
                       hMM_CC.GetNbinsX(),
                       hMM_CC.GetXaxis()->GetXmin(),
                       hMM_CC.GetXaxis()->GetXmax());
    hMM_frac_top2.SetDirectory(nullptr);
    for (int b = 1; b <= hMM_CC.GetNbinsX(); ++b) {
        const double cc = hMM_CC.GetBinContent(b);
        const double be = hMM_Best.GetBinContent(b);
        if (cc > 0.0) {
            double frac = be / cc;
            if (frac < 0.0) frac = 0.0;
            if (frac > 1.0) frac = 1.0;
            hMM_frac_top2.SetBinContent(b, frac);
        }
    }

    std::string outMM_frac = outPDF;
    if (auto dot = outMM_frac.find_last_of('.'); dot != std::string::npos)
        outMM_frac.insert(dot, "_MM_accFrac");
    else
        outMM_frac += "_MM_accFrac.pdf";

    TCanvas cFracTop2("cFracTop2", "f_acc vs M_X (TOP-2)", 900, 600);
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
    hMM_frac_top2.SetMinimum(0.0);
    hMM_frac_top2.SetMaximum(1.0);
    hMM_frac_top2.Draw("HIST");
    cFracTop2.Print(outMM_frac.c_str());
}

// ===== RooFit on m_gg (TOP-2/time-prox) : Gaussian + HistPdf(B_est) =====
{
    using namespace RooFit;

    RooRealVar mgg("mgg","m_{#gamma#gamma} (GeV)", hMgg_CC.GetXaxis()->GetXmin(), hMgg_CC.GetXaxis()->GetXmax());
    RooDataHist dh_data("dh_mgg","mgg CC data", RooArgList(mgg), &hMgg_CC);
    RooDataHist dh_b   ("dh_b_mgg","mgg B_{est}", RooArgList(mgg), &hMgg_Best);
    RooHistPdf  bpdf   ("bpdf_mgg","bkg pdf", RooArgSet(mgg), dh_b);

    RooRealVar mean("mean","mean", 0.135, 0.05, 0.30);
    RooRealVar sigma("sigma","sigma", 0.020, 0.005, 0.080);
    RooGaussian sig("sig_mgg","signal", mgg, mean, sigma);

    RooRealVar Ns("Ns","Ns", hMgg_CC.Integral()*0.5, 0.0, 1e9);
    RooRealVar Nb("Nb","Nb", hMgg_CC.Integral()*0.5, 0.0, 1e9);
    RooAddPdf model("model_mgg","sig+bkg", RooArgList(sig, bpdf), RooArgList(Ns, Nb));

    mgg.setRange("pk", 0.06, 0.24); // fit window around pi0

    std::unique_ptr<RooFitResult> fr(model.fitTo(dh_data, Save(true), Extended(true), Range("pk")));

    TCanvas cMggFit("cMggFit","mgg fit (TOP-2)", 900, 700);
    RooPlot* frg = mgg.frame(RooFit::Title("m_{#gamma#gamma}: CC = Gaussian + B_{est} (TOP-2)"));
    dh_data.plotOn(frg, RooFit::CutRange("pk"));
    model.plotOn(frg, RooFit::Range("pk"));
    model.plotOn(frg, RooFit::Components(bpdf), RooFit::LineStyle(kDashed), RooFit::LineColor(kBlue), RooFit::Range("pk"));
    model.plotOn(frg, RooFit::Components(sig),  RooFit::LineStyle(kSolid),  RooFit::LineColor(kRed),  RooFit::Range("pk"));
    frg->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    frg->GetYaxis()->SetTitle("Counts");
    frg->Draw();
// --- stats box (TOP-2) ---
auto pt = new TPaveText(0.62, 0.62, 0.98, 0.92, "NDC");
pt->SetFillColor(0);
pt->SetFillStyle(0);
pt->SetLineColor(0);
pt->SetTextAlign(12);
pt->SetTextSize(0.035);
pt->AddText(Form("#mu = %.4f #pm %.4f GeV",  mean.getVal(),  mean.getError()));
pt->AddText(Form("#sigma = %.4f #pm %.4f GeV", sigma.getVal(), sigma.getError()));
pt->AddText(Form("N_{s} = %.0f", Ns.getVal()));
pt->AddText(Form("N_{b} = %.0f", Nb.getVal()));
pt->Draw();                 // <— draw on pad (do NOT addObject)
gPad->Modified();
gPad->Update();

}

// ===== Summary printout of timing-template integrals =====

// (optional) quick counts to stdout
printf("[timing-template] CC=%g  V=%g  H=%g  D=%g  P=%g  |  Best=%g  Sub=%g\n",
       hMM_CC.Integral(), hMM_vert.Integral(), hMM_horiz.Integral(),
       hMM_diag.Integral(), hMM_pure.Integral(),
       hMM_Best.Integral(), hMM_subtracted.Integral());


       // ---------- Visualize timing-template classes in (t1, t2) ----------
std::string outTTmaps = outPDF;
if (auto dot = outTTmaps.find_last_of('.'); dot != std::string::npos)
    outTTmaps.insert(dot, "_TimingTemplate_maps");
else
    outTTmaps += "_TimingTemplate_maps.pdf";

TCanvas cTT("cTT","Timing-template maps",1600,900);
cTT.Divide(3,2); // ALL, CC, V, H, Dsame, Dopp

auto drawZ = [](){ gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); };

cTT.cd(1); drawZ(); hT12_all.Draw("COLZ");
cTT.cd(2); drawZ(); hT12_CC.Draw("COLZ");
cTT.cd(3); drawZ(); hT12_V.Draw("COLZ");
cTT.cd(4); drawZ(); hT12_H.Draw("COLZ");
cTT.cd(5); drawZ(); hT12_Dsame.Draw("COLZ");
cTT.cd(6); drawZ(); hT12_Dopp.Draw("COLZ");
cTT.Print((outTTmaps+"(").c_str()); // open multi-page

// ---------- Second page: ALL pairs with window boxes overlaid ----------
TCanvas cBoxes("cBoxes","Windows overlay",1200,1000);
gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
hT12_all.SetTitle("t_{#gamma1} vs t_{#gamma2} (ALL cand) with timing windows");
hT12_all.Draw("COLZ");

// Draw central box
TBox bC(C_LO, C_LO, C_HI, C_HI); bC.SetLineColor(kRed); bC.SetLineWidth(2); bC.SetFillStyle(0); bC.Draw("same");

// Draw vertical (C × A) and horizontal (A × C) bands, and diagonal squares
std::vector<TBox*> boxes; boxes.reserve(2*(NPOS+NNEG) + (NPOS+NNEG) + (NPOS*NNEG));
for (int k = 0; k < NPOS; ++k) {
    // vertical: C × pos
    boxes.push_back(new TBox(C_LO, posWins[k].lo, C_HI, posWins[k].hi));
    // horizontal: pos × C
    boxes.push_back(new TBox(posWins[k].lo, C_LO, posWins[k].hi, C_HI));
    // diagonal same-side square: pos × pos
    boxes.push_back(new TBox(posWins[k].lo, posWins[k].lo, posWins[k].hi, posWins[k].hi));
}
for (int k = 0; k < NNEG; ++k) {
    // vertical: C × neg
    boxes.push_back(new TBox(C_LO, negWins[k].lo, C_HI, negWins[k].hi));
    // horizontal: neg × C
    boxes.push_back(new TBox(negWins[k].lo, C_LO, negWins[k].hi, C_HI));
    // diagonal same-side square: neg × neg
    boxes.push_back(new TBox(negWins[k].lo, negWins[k].lo, negWins[k].hi, negWins[k].hi));
}
// opposite-side boxes are many; outline one example corner to illustrate
for (auto *bx : boxes) { bx->SetLineColor(kBlack); bx->SetLineStyle(3); bx->SetFillStyle(0); bx->Draw("same"); }

cBoxes.Print((outTTmaps+")").c_str()); // close multi-page

// ---------- ALL-PAIRS timing-template: area-normalized B_est and plot ----------
const double Wc_pairs = (C_HI - C_LO);

// totals over far windows (both sides)
double sumPos_pairs = 0.0, sumNeg_pairs = 0.0, sumSq_pairs = 0.0;
for (int k=0;k<NPOS;++k){ double w = posWins[k].hi - posWins[k].lo; sumPos_pairs += w; sumSq_pairs += w*w; }
for (int k=0;k<NNEG;++k){ double w = negWins[k].hi - negWins[k].lo; sumNeg_pairs += w; sumSq_pairs += w*w; }
const double sumAll_pairs = sumPos_pairs + sumNeg_pairs;

const double ACC_pairs = Wc_pairs*Wc_pairs;
const double AV_pairs  = Wc_pairs*sumAll_pairs;          // C×A (both sides)
const double AH_pairs  = Wc_pairs*sumAll_pairs;          // A×C (both sides)
const double AD_pairs  = sumSq_pairs;                    // A×A same (sum of squares)
const double AP_pairs  = (sumAll_pairs*sumAll_pairs)-AD_pairs; // A×A off-diagonals

const double aV_pairs = (AV_pairs>0)? ACC_pairs/AV_pairs : 0.0;
const double aH_pairs = (AH_pairs>0)? ACC_pairs/AH_pairs : 0.0;
const double aD_pairs = (AD_pairs>0)? ACC_pairs/AD_pairs : 0.0;
const double aP_pairs = (AP_pairs>0)? ACC_pairs/AP_pairs : 0.0;

TH1F hMM_Best_pairs("hMM_Best_pairs","Estimated accidental B_{est} (ALL pairs);M_{X} (GeV);Counts",
                    hMM_CC_pairs.GetNbinsX(), hMM_CC_pairs.GetXaxis()->GetXmin(), hMM_CC_pairs.GetXaxis()->GetXmax());
hMM_Best_pairs.Reset();
hMM_Best_pairs.Add(&hMM_diag_pairs,   aD_pairs);
hMM_Best_pairs.Add(&hMM_vert_pairs,   aV_pairs);
hMM_Best_pairs.Add(&hMM_horiz_pairs,  aH_pairs);
hMM_Best_pairs.Add(&hMM_pure_pairs,  -aP_pairs);

hMM_subtracted_pairs.Reset();
hMM_subtracted_pairs.Add(&hMM_CC_pairs, 1.0);
hMM_subtracted_pairs.Add(&hMM_Best_pairs, -1.0);

// --- plot
std::string outMM_pairs = outPDF;
if (auto dot = outMM_pairs.find_last_of('.'); dot != std::string::npos)
    outMM_pairs.insert(dot, "_MM_timingTemplate_allPairs");
else
    outMM_pairs += "_MM_timingTemplate_allPairs.pdf";

TCanvas cMMpairs("cMMpairs","MM timing-template (ALL pairs)",1200,600);
cMMpairs.Divide(2,1);

cMMpairs.cd(1);
gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
hMM_CC_pairs.SetLineColor(kBlack); hMM_CC_pairs.SetLineWidth(2);
hMM_CC_pairs.SetTitle("ALL pairs: CC and B_{est}");
hMM_CC_pairs.Draw("HIST");
hMM_Best_pairs.SetLineColor(kRed+1); hMM_Best_pairs.SetLineStyle(2); hMM_Best_pairs.SetLineWidth(2);
hMM_Best_pairs.Draw("HIST SAME");
auto legP = new TLegend(0.60,0.72,0.88,0.88);
legP->AddEntry(&hMM_CC_pairs,  "CC (ALL pairs)", "l");
legP->AddEntry(&hMM_Best_pairs,"B_{est} (ALL pairs)", "l");
legP->Draw();

cMMpairs.cd(2);
gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
hMM_subtracted_pairs.SetLineColor(kBlue+1); hMM_subtracted_pairs.SetLineWidth(2);
hMM_subtracted_pairs.SetTitle("ALL pairs: CC - B_{est}");
hMM_subtracted_pairs.Draw("HIST");

cMMpairs.Print(outMM_pairs.c_str());

// ===== RooFit: Gaussian peak + template background (ALL pairs) =====
{
    using namespace RooFit;

    RooRealVar mx("mx","M_{X} (GeV)", hMM_CC_pairs.GetXaxis()->GetXmin(), hMM_CC_pairs.GetXaxis()->GetXmax());

    RooDataHist dh_data("dh_data","CC data (ALL pairs)", RooArgList(mx), &hMM_CC_pairs);
    RooDataHist dh_b("dh_b","B_{est} template (ALL pairs)", RooArgList(mx), &hMM_Best_pairs);
    RooHistPdf  bpdf("bpdf","bkg pdf", RooArgSet(mx), dh_b);

    RooRealVar mean("mean","mean", 0.94, 0.80, 1.15);
    RooRealVar sigma("sigma","sigma", 0.10, 0.02, 0.30);
    RooGaussian sig("sig","signal", mx, mean, sigma);

    RooRealVar Ns("Ns","Ns", hMM_CC_pairs.Integral()*0.5, 0.0, 1e9);
    RooRealVar Nb("Nb","Nb", hMM_CC_pairs.Integral()*0.5, 0.0, 1e9);

    RooAddPdf model("model","sig+bkg", RooArgList(sig, bpdf), RooArgList(Ns, Nb));

    mx.setRange("pk", 0.60, 1.40);
    std::unique_ptr<RooFitResult> fr(model.fitTo(dh_data, Save(true), Extended(true), Range("pk")));

    TCanvas cFit("cFit","MM fit (ALL pairs)", 900, 700);
  RooPlot* frame = mx.frame(RooFit::Title("MM: CC data = Gaussian + B_{est} (ALL pairs)"));
dh_data.plotOn(frame, RooFit::CutRange("pk"));
model.plotOn(frame, RooFit::Range("pk"));
model.plotOn(frame, RooFit::Components(bpdf), RooFit::LineStyle(kDashed), RooFit::LineColor(kBlue), RooFit::Range("pk"));
model.plotOn(frame, RooFit::Components(sig),  RooFit::LineStyle(kSolid),  RooFit::LineColor(kRed),  RooFit::Range("pk"));
 
    frame->GetXaxis()->SetTitle("M_{X} (GeV)");
    frame->GetYaxis()->SetTitle("Counts");
    frame->Draw();

    std::string outFitPairs = outPDF;
    if (auto dot = outFitPairs.find_last_of('.'); dot != std::string::npos)
        outFitPairs.insert(dot, "_MM_fit_template_allPairs");
    else
        outFitPairs += "_MM_fit_template_allPairs.pdf";
    cFit.Print(outFitPairs.c_str());

    std::cout << "[MM fit ALL-pairs] mean=" << mean.getVal() << " ± " << mean.getError()
              << "  sigma=" << sigma.getVal() << " ± " << sigma.getError()
              << "  Ns=" << Ns.getVal() << "  Nb=" << Nb.getVal() << "\n";
}

// ===== ALL-PAIRS timing-template: integrated summary + f_acc(MX) plot =====
{
    auto HInt = [](const TH1F &h)->double { return h.Integral(1, h.GetNbinsX()); };

    const double I_CC  = HInt(hMM_CC_pairs);
    const double I_V   = HInt(hMM_vert_pairs);
    const double I_H   = HInt(hMM_horiz_pairs);
    const double I_D   = HInt(hMM_diag_pairs);
    const double I_P   = HInt(hMM_pure_pairs);
    const double I_B   = HInt(hMM_Best_pairs);
    const double I_Sub = HInt(hMM_subtracted_pairs);
    const double f_acc = (I_CC > 0.0) ? (I_B / I_CC) : 0.0;

    printf("[TimingTemplate ALL-PAIRS]  CC=%g  V=%g  H=%g  D=%g  P=%g  |  B_est=%g  Sub=%g  f_acc=%g\n",
           I_CC, I_V, I_H, I_D, I_P, I_B, I_Sub, f_acc);

    TH1F hMM_frac_pairs("hMM_frac_pairs",
                        "Accidental fraction f_{acc}(M_{X}) (ALL pairs);M_{X} (GeV);f_{acc}",
                        hMM_CC_pairs.GetNbinsX(),
                        hMM_CC_pairs.GetXaxis()->GetXmin(),
                        hMM_CC_pairs.GetXaxis()->GetXmax());
    hMM_frac_pairs.SetDirectory(nullptr);
    for (int b = 1; b <= hMM_CC_pairs.GetNbinsX(); ++b) {
        const double cc = hMM_CC_pairs.GetBinContent(b);
        const double be = hMM_Best_pairs.GetBinContent(b);
        if (cc > 0.0) {
            double frac = be / cc;
            if (frac < 0.0) frac = 0.0;
            if (frac > 1.0) frac = 1.0;
            hMM_frac_pairs.SetBinContent(b, frac);
        }
    }

    std::string outMM_frac_pairs = outPDF;
    if (auto dot = outMM_frac_pairs.find_last_of('.'); dot != std::string::npos)
        outMM_frac_pairs.insert(dot, "_MM_accFrac_allPairs");
    else
        outMM_frac_pairs += "_MM_accFrac_allPairs.pdf";

    TCanvas cFracPairs("cFracPairs", "f_acc vs M_X (ALL pairs)", 900, 600);
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
    hMM_frac_pairs.SetMinimum(0.0);
    hMM_frac_pairs.SetMaximum(1.0);
    hMM_frac_pairs.Draw("HIST");
    
    cFracPairs.Print(outMM_frac_pairs.c_str());
}

// ===== RooFit on m_gg (ALL pairs) : Gaussian + HistPdf(B_est) =====
{
    using namespace RooFit;

    RooRealVar mgg("mgg","m_{#gamma#gamma} (GeV)", hMgg_CC_pairs.GetXaxis()->GetXmin(), hMgg_CC_pairs.GetXaxis()->GetXmax());
    RooDataHist dh_data("dh_mgg_pairs","mgg CC (ALL pairs)", RooArgList(mgg), &hMgg_CC_pairs);
    RooDataHist dh_b   ("dh_b_mgg_pairs","mgg B_{est} (ALL pairs)", RooArgList(mgg), &hMgg_Best_pairs);
    RooHistPdf  bpdf   ("bpdf_mgg_pairs","bkg pdf", RooArgSet(mgg), dh_b);

    RooRealVar mean("mean","mean", 0.135, 0.05, 0.30);
    RooRealVar sigma("sigma","sigma", 0.020, 0.005, 0.080);
    RooGaussian sig("sig_mgg_pairs","signal", mgg, mean, sigma);

    RooRealVar Ns("Ns","Ns", hMgg_CC_pairs.Integral()*0.5, 0.0, 1e9);
    RooRealVar Nb("Nb","Nb", hMgg_CC_pairs.Integral()*0.5, 0.0, 1e9);
    RooAddPdf model("model_mgg_pairs","sig+bkg", RooArgList(sig, bpdf), RooArgList(Ns, Nb));

    mgg.setRange("pk", 0.06, 0.24);

    std::unique_ptr<RooFitResult> fr(model.fitTo(dh_data, Save(true), Extended(true), Range("pk")));

    TCanvas cMggFitPairs("cMggFitPairs","mgg fit (ALL pairs)", 900, 700);
    RooPlot* frg = mgg.frame(RooFit::Title("m_{#gamma#gamma}: CC = Gaussian + B_{est} (ALL pairs)"));
    dh_data.plotOn(frg, RooFit::CutRange("pk"));
    model.plotOn(frg, RooFit::Range("pk"));
    model.plotOn(frg, RooFit::Components(bpdf), RooFit::LineStyle(kDashed), RooFit::LineColor(kBlue), RooFit::Range("pk"));
    model.plotOn(frg, RooFit::Components(sig),  RooFit::LineStyle(kSolid),  RooFit::LineColor(kRed),  RooFit::Range("pk"));
    frg->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
    frg->GetYaxis()->SetTitle("Counts");
    frg->Draw();

// --- stats box (ALL pairs) ---
// roomy lower-right box (taller + inset from edges)
auto ptPairs = new TPaveText(0.64, 0.22, 0.98, 0.62, "NDC");
ptPairs->SetFillColor(0);
ptPairs->SetFillStyle(0);
ptPairs->SetLineColor(0);
ptPairs->SetBorderSize(0);
ptPairs->SetTextAlign(12);
ptPairs->SetTextSize(0.030);   // a bit smaller so 4 lines fit comfortably
ptPairs->SetMargin(0.12);      // inner left margin so text isn’t flush

ptPairs->AddText(Form("#mu = %.4f #pm %.4f GeV",  mean.getVal(),  mean.getError()));
ptPairs->AddText(Form("#sigma = %.4f #pm %.4f GeV", sigma.getVal(), sigma.getError()));
ptPairs->AddText(Form("N_{s} = %.0f", Ns.getVal()));
ptPairs->AddText(Form("N_{b} = %.0f", Nb.getVal()));

ptPairs->Draw("same");
gPad->Modified(); gPad->Update();


    std::string outFitMggPairs = outPDF;
    if (auto dot = outFitMggPairs.find_last_of('.'); dot != std::string::npos) outFitMggPairs.insert(dot, "_Mgg_fit_template_allPairs");
    else outFitMggPairs += "_Mgg_fit_template_allPairs.pdf";
    cMggFitPairs.Print(outFitMggPairs.c_str());

    std::cout << "[mgg fit ALL-pairs] mean=" << mean.getVal() << " ± " << mean.getError()
              << "  sigma=" << sigma.getVal() << " ± " << sigma.getError()
              << "  Ns=" << Ns.getVal() << "  Nb=" << Nb.getVal() << "\n";
}


// quick integrals
printf("[ALL-pairs templ] CC=%g  V=%g  H=%g  D=%g  P=%g  |  Best=%g  Sub=%g\n",
       hMM_CC_pairs.Integral(), hMM_vert_pairs.Integral(), hMM_horiz_pairs.Integral(),
       hMM_diag_pairs.Integral(), hMM_pure_pairs.Integral(),
       hMM_Best_pairs.Integral(), hMM_subtracted_pairs.Integral());


// --- ALL-pairs timing-template maps page ---
std::string outTTpairs = outPDF;
if (auto dot = outTTpairs.find_last_of('.'); dot != std::string::npos)
    outTTpairs.insert(dot, "_TimingTemplate_maps_allPairs");
else
    outTTpairs += "_TimingTemplate_maps_allPairs.pdf";

TCanvas cTTpairs("cTTpairs","Timing-template maps (ALL pairs)",1600,900);
cTTpairs.Divide(3,2);
auto padZ = [](){ gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); };

// Build visualization for Option A: what pairs actually feed B_pairs
hPairs_usedA.Reset("ICES");
hPairs_usedA.Add(&hT12_V_pairs,      1.0);   // or aV if you want coefficient-weighted look
hPairs_usedA.Add(&hT12_H_pairs,      1.0);   // or aH
hPairs_usedA.Add(&hT12_Dsame_pairs,  1.0);   // or aA
hPairs_usedA.Add(&hT12_Dopp_pairs,   1.0);   // or aA


cTTpairs.cd(1); padZ(); hT12_allPairs.Draw("COLZ");
cTTpairs.cd(2); padZ(); hT12_CC_pairs.Draw("COLZ");
cTTpairs.cd(3); padZ(); hT12_V_pairs.Draw("COLZ");
cTTpairs.cd(4); padZ(); hT12_H_pairs.Draw("COLZ");
cTTpairs.cd(5); padZ(); hT12_Dsame_pairs.Draw("COLZ");
cTTpairs.cd(6); padZ(); hT12_Dopp_pairs.Draw("COLZ");

cTTpairs.Print(outTTpairs.c_str());

// --- Compare standard MM to time-proximity selection
std::string outMM_timePick = outPDF;
if (auto dot = outMM_timePick.find_last_of('.'); dot != std::string::npos)
    outMM_timePick.insert(dot, "_MM_timeProximity");
else
    outMM_timePick += "_MM_timeProximity.pdf";

TCanvas cMMtp("cMMtp","MM: time-proximity vs reference", 900, 600);
gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);

// choose a reference to overlay; hMissMass is common in your script
hMissMass.SetLineColor(kBlack); hMissMass.SetLineWidth(2);
hMissMass.SetTitle("Missing Mass;M_{X} (GeV);Counts");
hMissMass.Draw("HIST");

hMM_timePick.SetLineColor(kBlue+1); hMM_timePick.SetLineStyle(2); hMM_timePick.SetLineWidth(2);
hMM_timePick.Draw("HIST SAME");

auto legTP = new TLegend(0.58,0.72,0.88,0.88);
legTP->AddEntry(&hMissMass,     "MM (your current selection)", "l");
legTP->AddEntry(&hMM_timePick,  "MM (time-proximity pair)",   "l");
legTP->Draw();

cMMtp.Print(outMM_timePick.c_str());

// --- Plot the time-proximity t1 vs t2 map
std::string outTtimePick = outPDF;
if (auto dot = outTtimePick.find_last_of('.'); dot != std::string::npos)
    outTtimePick.insert(dot, "_Timing_timeProximity");
else
    outTtimePick += "_Timing_timeProximity.pdf";

TCanvas cTP("cTP","Time-proximity: t1 vs t2", 900, 800);
gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);

hT12_timePick.SetTitle("t_{#gamma1} vs t_{#gamma2} (time-proximity pair)");
hT12_timePick.Draw("COLZ");

// overlay central coincidence box [149,151] × [149,151]
TBox bC_tp(C_LO, C_LO, C_HI, C_HI);
bC_tp.SetLineColor(kRed);
bC_tp.SetLineWidth(2);
bC_tp.SetFillStyle(0);
bC_tp.Draw("same");

cTP.Print(outTtimePick.c_str());

// -------------------- Build subtracted spectra and overlay --------------------
const double Wacc = totalWidth(posWins, NPOS) + totalWidth(negWins, NNEG);
const double Wsig = (C_HI - C_LO);
const double alpha = (Wacc > 0) ? (Wsig / (2.0*Wacc)) : 0.0;

printf("[overlay] ALL(sig,acc)=(%g,%g)  TOP2(sig,acc)=(%g,%g)  TEMPLATE(CC,B_est,Sub)=(%g,%g,%g)\n",
       hMM_all_sig->Integral(),  hMM_all_acc->Integral(),
       hMM_top2_sig->Integral(), hMM_top2_acc->Integral(),
       hMM_CC.Integral(),        hMM_Best.Integral(), hMM_subtracted.Integral());

// --- Build overlay histograms ---
// keep All-pairs and Top-2 as window-normalized (alpha) spectra:
TH1* hMM_top2 = (TH1*)hMM_top2_sig->Clone("hMM_top2");   hMM_top2->Add(hMM_top2_acc, -alpha);
TH1* hMM_all  = (TH1*)hMM_all_sig ->Clone("hMM_all");    hMM_all ->Add(hMM_all_acc,  -alpha);

// BUT: the Timing-template curve should come from your per-event template result (CC − B_est),
// which you already computed as `hMM_subtracted` in the “TimingTemplate TOP-2” section.
TH1* hMM_temp = (TH1*)hMM_subtracted.Clone("hMM_temp");


hMM_top2->SetLineColor(kRed+1);   hMM_top2->SetLineWidth(2);
hMM_all ->SetLineColor(kBlue+1);  hMM_all ->SetLineWidth(2);
hMM_temp->SetLineColor(kGreen+2); hMM_temp->SetLineWidth(2);

// Choose a common y-max
double ymax = std::max({ hMM_top2->GetMaximum(), hMM_all->GetMaximum(), hMM_temp->GetMaximum() });

// Canvas
// ======================= MM overlay: Top-2, All-pairs, Old-Template, + A method =======================
{
  // local prefix for output name
  TString __prefix = argv[2]; __prefix.ReplaceAll(".pdf","");

  // helper: find a TH1 by any of several known names
  auto findH1 = [](const std::initializer_list<const char*>& names)->TH1*{
    for (auto nm : names) {
      if (!nm) continue;
      TObject* o = ROOT::GetROOT() ? ROOT::GetROOT()->FindObject(nm) : nullptr;
      if (!o) continue;
      if (auto h = dynamic_cast<TH1*>(o)) return h;
    }
    return nullptr;
  };

  // your existing three (use whatever names your file actually uses)
  TH1* hTop2 = findH1({"hMM_top2_sub","hMM_top2_SUB","hMM_top2"});     // Top-2 (subtracted) 
  TH1* hAll  = findH1({"hMM_all_sub","hMM_all_SUB","hMM_all"});        // All-pairs (subtracted)
  TH1* hOld  = findH1({"hMM_temp_sub","hMM_temp_SUB","hMM_temp"});     // Old timing template (subtracted)

  // A-method (new) — already built above; use pointer if available
  TH1* hA    = hTplA_sub;

  // If any of the three legacy histos are actually your in-memory variables, you can also hard-wire:
  // if (!hTop2) hTop2 = /* pointer to your in-code Top-2 subtracted histo */;
  // if (!hAll)  hAll  = /* pointer to your in-code All-pairs subtracted histo */;
  // if (!hOld)  hOld  = /* pointer to your in-code Old-template subtracted histo */;

  // Build a common frame for axes & ymax
  double ymax = 0.0;
  auto updMax = [&](TH1* h){ if (h) ymax = std::max(ymax, h->GetMaximum()); };
  updMax(hTop2); updMax(hAll); updMax(hOld); updMax(hA);

  // If nothing found, skip gracefully
  if (ymax <= 0.0) {
    printf("[MM overlay] No histograms found to overlay; skipping page.\n");
  } else {
    TCanvas cMMOverlay("cMMOverlay","MM Overlay (methods)", 1000, 700);
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);

    // Create a blank axis frame using the x-range from whichever hist is non-null
    TH1* href =
      hTop2 ? hTop2 :
      hAll  ? hAll  :
      hOld  ? hOld  :
      hA; // must be non-null since ymax>0
    TH1* hAxis = (TH1*)href->Clone("hAxis_overlay");
    hAxis->Reset("ICESM");
    hAxis->SetLineColor(0); hAxis->SetLineWidth(0);
    hAxis->SetTitle("Missing Mass — method comparison;M_{X} (GeV);Counts");
    hAxis->SetMaximum(1.20*ymax);
    hAxis->Draw("AXIS");

    // Style and draw each (skip nulls)
    auto styleDraw = [](TH1* h, Color_t c, Style_t ls=1, int lw=2){
      if (!h) return;
      h->SetLineColor(c); h->SetLineStyle(ls); h->SetLineWidth(lw);
      h->SetMarkerStyle(0);
      h->Draw("HIST SAME");
    };

    // colors: keep old three distinct, add A in a new color
    styleDraw(hTop2, kRed+1,    1, 2);
    styleDraw(hAll,  kBlue+1,   1, 2);
    styleDraw(hOld,  kGreen+2,  2, 2);
    styleDraw(hA,    kMagenta+2,1, 2);  // A-method subtraction

    // legend
    auto leg = new TLegend(0.55,0.62,0.88,0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    if (hTop2) leg->AddEntry(hTop2, "Top-2 method (SUB)", "l");
    if (hAll)  leg->AddEntry(hAll,  "All-pairs method (SUB)", "l");
    if (hOld)  leg->AddEntry(hOld,  "Old timing template (SUB)", "l");
    if (hA)    leg->AddEntry(hA,    "A-method (CC - B_{est})", "l");
    leg->Draw();

    cMMOverlay.Print(Form("%s_MM_overlay_methods.pdf", __prefix.Data()));
    delete hAxis;
  }
}


// =====================================================================
// Missing-mass breakdowns: RAW vs BACKGROUND vs SUBTRACTED (three plots)
// Insert right after the methods overlay legend is drawn.
// =====================================================================

// local prefix derived from argv[2] (avoid relying on outer-scope vars)
TString __mm_prefix_local = argv[2];
__mm_prefix_local.ReplaceAll(".pdf","");

// ---------- Top-2 breakdown ----------
{
  // RAW = signal-window contents; BKG = alpha * accidentals; SUB = RAW - BKG
  TH1D* h_top2_raw = (TH1D*)hMM_top2_sig->Clone("h_top2_raw");   h_top2_raw->SetDirectory(nullptr); h_top2_raw->Sumw2();
  TH1D* h_top2_bkg = (TH1D*)hMM_top2_acc->Clone("h_top2_bkg");   h_top2_bkg->SetDirectory(nullptr); h_top2_bkg->Sumw2(); h_top2_bkg->Scale(alpha);
  TH1D* h_top2_sub = (TH1D*)h_top2_raw->Clone("h_top2_sub");     h_top2_sub->Add(h_top2_bkg, -1.0);

  // style
  h_top2_raw->SetLineColor(kBlack);    h_top2_raw->SetLineWidth(2);
  h_top2_bkg->SetLineColor(kOrange+7); h_top2_bkg->SetLineWidth(2);
  h_top2_sub->SetLineColor(kBlue+1);   h_top2_sub->SetLineWidth(2);

  // y-range (include negatives if any)
  auto __minbin = [](TH1* h){ double m=1e300; for(int b=1;b<=h->GetNbinsX();++b) m = std::min(m, h->GetBinContent(b)); return m; };
  const double __ymax_top2 = std::max({ h_top2_raw->GetMaximum(), h_top2_bkg->GetMaximum(), h_top2_sub->GetMaximum() });
  const double __ymin_top2 = std::min({ __minbin(h_top2_raw), __minbin(h_top2_bkg), __minbin(h_top2_sub), 0.0 });

  TCanvas cMM_Top2_breakdown("cMM_Top2_breakdown","Top-2: raw / bkg / sub",1000,750);
  cMM_Top2_breakdown.SetRightMargin(0.08); cMM_Top2_breakdown.SetTopMargin(0.08);

  h_top2_raw->SetTitle("Top-2 Missing Mass: raw vs. background vs. subtracted");
  h_top2_raw->GetXaxis()->SetTitle("M_{X} (GeV)");
  h_top2_raw->GetYaxis()->SetTitle("Counts");
  h_top2_raw->SetMaximum(1.15*std::max(1e-12, __ymax_top2));
  h_top2_raw->SetMinimum(__ymin_top2 < 0 ? 1.15*__ymin_top2 : 0.0);

  h_top2_raw->Draw("HIST");
  h_top2_bkg->Draw("HIST SAME");
  h_top2_sub->Draw("HIST SAME");

  auto leg_top2 = new TLegend(0.58,0.68,0.88,0.88);
  leg_top2->SetBorderSize(0); leg_top2->SetFillStyle(0);
  leg_top2->AddEntry(h_top2_raw, "Raw (signal window)", "l");
  leg_top2->AddEntry(h_top2_bkg, Form("Background (#alpha=%.4f)", alpha), "l");
  leg_top2->AddEntry(h_top2_sub, "Subtracted = Raw - Bkg", "l");
  leg_top2->Draw();

  cMM_Top2_breakdown.Print(Form("%s_MM_breakdown_top2.pdf", __mm_prefix_local.Data()));
}

// ---------- All-pairs breakdown ----------
{
  TH1D* h_all_raw = (TH1D*)hMM_all_sig->Clone("h_all_raw");   h_all_raw->SetDirectory(nullptr); h_all_raw->Sumw2();
  TH1D* h_all_bkg = (TH1D*)hMM_all_acc->Clone("h_all_bkg");   h_all_bkg->SetDirectory(nullptr); h_all_bkg->Sumw2(); h_all_bkg->Scale(alpha);
  TH1D* h_all_sub = (TH1D*)h_all_raw->Clone("h_all_sub");     h_all_sub->Add(h_all_bkg, -1.0);

  h_all_raw->SetLineColor(kBlack);    h_all_raw->SetLineWidth(2);
  h_all_bkg->SetLineColor(kOrange+7); h_all_bkg->SetLineWidth(2);
  h_all_sub->SetLineColor(kGreen+2);  h_all_sub->SetLineWidth(2);

  auto __minbin = [](TH1* h){ double m=1e300; for(int b=1;b<=h->GetNbinsX();++b) m = std::min(m, h->GetBinContent(b)); return m; };
  const double __ymax_all = std::max({ h_all_raw->GetMaximum(), h_all_bkg->GetMaximum(), h_all_sub->GetMaximum() });
  const double __ymin_all = std::min({ __minbin(h_all_raw), __minbin(h_all_bkg), __minbin(h_all_sub), 0.0 });

  TCanvas cMM_All_breakdown("cMM_All_breakdown","All-pairs: raw / bkg / sub",1000,750);
  cMM_All_breakdown.SetRightMargin(0.08); cMM_All_breakdown.SetTopMargin(0.08);

  h_all_raw->SetTitle("All-pairs Missing Mass: raw vs. background vs. subtracted");
  h_all_raw->GetXaxis()->SetTitle("M_{X} (GeV)");
  h_all_raw->GetYaxis()->SetTitle("Counts");
  h_all_raw->SetMaximum(1.15*std::max(1e-12, __ymax_all));
  h_all_raw->SetMinimum(__ymin_all < 0 ? 1.15*__ymin_all : 0.0);

  h_all_raw->Draw("HIST");
  h_all_bkg->Draw("HIST SAME");
  h_all_sub->Draw("HIST SAME");

  auto leg_all = new TLegend(0.58,0.68,0.88,0.88);
  leg_all->SetBorderSize(0); leg_all->SetFillStyle(0);
  leg_all->AddEntry(h_all_raw, "Raw (signal window)", "l");
  leg_all->AddEntry(h_all_bkg, Form("Background (#alpha=%.4f)", alpha), "l");
  leg_all->AddEntry(h_all_sub, "Subtracted = Raw - Bkg", "l");
  leg_all->Draw();

  cMM_All_breakdown.Print(Form("%s_MM_breakdown_allpairs.pdf", __mm_prefix_local.Data()));
}

// ---------- Timing-template breakdown (robust, with diagnostics) ----------
{
  // RAW and SUB come from your template outputs (both are TH1F *objects*, not pointers)
  TH1* h_temp_raw = (TH1*)hMM_CC.Clone("h_temp_raw");                h_temp_raw->SetDirectory(nullptr); h_temp_raw->Sumw2();
  TH1* h_temp_sub = (TH1*)hMM_subtracted.Clone("h_temp_sub");        h_temp_sub->SetDirectory(nullptr); h_temp_sub->Sumw2();

  // Background from template truth: bkg_true = RAW - SUB
  TH1* h_temp_bkg = (TH1*)h_temp_raw->Clone("h_temp_bkg");           h_temp_bkg->Add(h_temp_sub, -1.0);

  // --- Optional cross-check: compose background from components if present in this scope ---
  // Some builds name the extra component "hMM_pile" vs "hMM_pure". We handle both safely.
  TH1* h_temp_bkg_comp = nullptr;
  bool have_comp = false;
  {
    // We’ll detect presence via address-of; if not compiled in this block, the name won’t exist.
    // Comment/uncomment ONE of the following compilation paths depending on your file:

    // Path A: components exist as TH1F *objects* in this scope:
    // #define HAS_TEMPLATE_COMPONENTS

#ifdef HAS_TEMPLATE_COMPONENTS
    TH1* h_comp_sum = (TH1*)hMM_vert.Clone("h_temp_comp_sum"); h_comp_sum->SetDirectory(nullptr); h_comp_sum->Sumw2();
    h_comp_sum->Add(&hMM_horiz);
    h_comp_sum->Add(&hMM_diag);
    // If your build has "hMM_pile" use that; else if it has "hMM_pure" use that; else skip.
    // Uncomment ONE of these to match your file:
    // h_comp_sum->Add(&hMM_pile);
    // h_comp_sum->Add(&hMM_pure);

    // Scale components to match the integral of the true background (derives an effective f_acc)
    h_temp_bkg_comp = (TH1*)h_comp_sum->Clone("h_temp_bkg_comp"); h_temp_bkg_comp->SetDirectory(nullptr); h_temp_bkg_comp->Sumw2();
    const double int_true = h_temp_bkg->Integral();
    const double int_comp = h_comp_sum->Integral();
    const double facc_eff = (int_comp > 0 ? int_true / int_comp : 0.0);
    if (int_comp > 0) h_temp_bkg_comp->Scale(facc_eff);
    have_comp = (int_comp > 0);
#endif
  }

  // Style
  h_temp_raw->SetLineColor(kBlack);    h_temp_raw->SetLineWidth(2);
  h_temp_bkg->SetLineColor(kOrange+7); h_temp_bkg->SetLineStyle(7); h_temp_bkg->SetLineWidth(2);
  h_temp_sub->SetLineColor(kRed+1);    h_temp_sub->SetLineWidth(2);
  if (have_comp) { h_temp_bkg_comp->SetLineColor(kGreen+2); h_temp_bkg_comp->SetLineStyle(9); h_temp_bkg_comp->SetLineWidth(2); }

  // Axis range (include negatives if any)
  auto __minbin = [](TH1* h){ double m=1e300; for (int b=1;b<=h->GetNbinsX();++b) m = std::min(m, h->GetBinContent(b)); return m; };
  const double __ymax_tpl = std::max({
    h_temp_raw->GetMaximum(),
    h_temp_bkg->GetMaximum(),
    h_temp_sub->GetMaximum(),
    have_comp ? h_temp_bkg_comp->GetMaximum() : 0.0
  });
  const double __ymin_tpl = std::min({
    __minbin(h_temp_raw),
    __minbin(h_temp_bkg),
    __minbin(h_temp_sub),
    have_comp ? __minbin(h_temp_bkg_comp) : 0.0,
    0.0
  });

  TCanvas cMM_Temp_breakdown("cMM_Temp_breakdown","Timing-template: raw / bkg / sub",1000,750);
  cMM_Temp_breakdown.SetRightMargin(0.08); cMM_Temp_breakdown.SetTopMargin(0.08);

  h_temp_raw->SetTitle("Timing-template Missing Mass: raw vs. background vs. subtracted");
  h_temp_raw->GetXaxis()->SetTitle("M_{X} (GeV)");
  h_temp_raw->GetYaxis()->SetTitle("Counts");
  h_temp_raw->SetMaximum(1.15*std::max(1e-12, __ymax_tpl));
  h_temp_raw->SetMinimum(__ymin_tpl < 0 ? 1.15*__ymin_tpl : 0.0);

  h_temp_raw->Draw("HIST");
  h_temp_bkg->Draw("HIST SAME");
  h_temp_sub->Draw("HIST SAME");
  if (have_comp) h_temp_bkg_comp->Draw("HIST SAME");

  auto leg_tpl = new TLegend(0.56,0.64,0.88,0.88);
  leg_tpl->SetBorderSize(0); leg_tpl->SetFillStyle(0);
  leg_tpl->AddEntry(h_temp_raw, "Raw (CC)", "l");
  leg_tpl->AddEntry(h_temp_bkg, "Background (RAW - SUB)", "l");
  if (have_comp) leg_tpl->AddEntry(h_temp_bkg_comp, "Background (scaled V+H+D+X)", "l");
  leg_tpl->AddEntry(h_temp_sub, "Subtracted = Raw - Bkg", "l");
  leg_tpl->Draw();

  // Diagnostics: integrals and peak-bin fractions
  const int ibin_pk = h_temp_raw->GetMaximumBin();
  const double raw_pk = h_temp_raw->GetBinContent(ibin_pk);
  const double bkg_pk = h_temp_bkg->GetBinContent(ibin_pk);
  const double sub_pk = h_temp_sub->GetBinContent(ibin_pk);
  printf("[tpl breakdown] integrals  RAW=%g  BKG=%g  SUB=%g  | peak-bin (x=%.4f): raw=%g  bkg=%g  sub=%g  frac_bkg=%.4f\n",
         h_temp_raw->Integral(), h_temp_bkg->Integral(), h_temp_sub->Integral(),
         h_temp_raw->GetBinCenter(ibin_pk), raw_pk, bkg_pk, sub_pk,
         (raw_pk>0? bkg_pk/raw_pk : 0.0));
  if (have_comp) {
    const double bkgc_pk = h_temp_bkg_comp->GetBinContent(ibin_pk);
    printf("[tpl breakdown] components check: BKG_true_int=%g  BKG_comp_int=%g  peak_comp=%g  (should be close if same selection)\n",
           h_temp_bkg->Integral(), h_temp_bkg_comp->Integral(), bkgc_pk);
  }

  cMM_Temp_breakdown.Print(Form("%s_MM_breakdown_template.pdf", __mm_prefix_local.Data()));
}

// =====================================================================
// Timing-template (RELAXED) background: widen central exposure virtually
// without changing the selection. Produces a second breakdown figure.
// Paste this immediately after the strict template breakdown block.
// =====================================================================
{
  // Local output prefix (no reliance on outer vars)
  TString __tpl_prefix2 = argv[2]; __tpl_prefix2.ReplaceAll(".pdf","");

  // 1) Reuse RAW and SUB from your strict template products
  TH1* h_tpl_raw_strict = (TH1*)hMM_CC.Clone("h_tpl_raw_strict");            h_tpl_raw_strict->SetDirectory(nullptr); h_tpl_raw_strict->Sumw2();
  TH1* h_tpl_sub_strict = (TH1*)hMM_subtracted.Clone("h_tpl_sub_strict");    h_tpl_sub_strict->SetDirectory(nullptr); h_tpl_sub_strict->Sumw2();

  // The strict template background we already verified as: B_true = RAW - SUB
  TH1* h_tpl_bkg_strict = (TH1*)h_tpl_raw_strict->Clone("h_tpl_bkg_strict"); h_tpl_bkg_strict->Add(h_tpl_sub_strict, -1.0);

  // 2) Compute current (strict) exposure and a RELAXED exposure
  //    alpha_strict = Wsig / Wacc, with Wsig = C_HI - C_LO and Wacc = sum of sideband widths
  double __Wacc = 0.0;
  for (int i=0;i<NPOS;++i) __Wacc += (posWins[i].hi - posWins[i].lo);
  for (int i=0;i<NNEG;++i) __Wacc += (negWins[i].hi - negWins[i].lo);
  const double __Wsig_strict = (C_HI - C_LO);
  const double __alpha_strict = (__Wacc > 0 ? __Wsig_strict / __Wacc : 0.0);

  // --- KNOB: edge padding for the coincidence box (in ns), used ONLY for background scaling
  const double TPL_EDGE_PAD_NS = 1.0;  // try 0.25–1.0 ns; does not change which events are selected

  const double __Wsig_relax = std::max(0.0, __Wsig_strict + 2.0*TPL_EDGE_PAD_NS);
  const double __alpha_relax = (__Wacc > 0 ? __Wsig_relax / __Wacc : 0.0);

  // Scale factor to boost the background consistently across m:
  const double __k_relax = (__alpha_strict > 0 ? (__alpha_relax / __alpha_strict) : 1.0);

  // 3) Build RELAXED background and subtracted spectrum
  TH1* h_tpl_bkg_relax = (TH1*)h_tpl_bkg_strict->Clone("h_tpl_bkg_relax"); h_tpl_bkg_relax->Scale(__k_relax);
  TH1* h_tpl_sub_relax = (TH1*)h_tpl_raw_strict->Clone("h_tpl_sub_relax"); h_tpl_sub_relax->Add(h_tpl_bkg_relax, -1.0);

  // 4) Plot breakdown (RAW, BKG_relax, SUB_relax) with diagnostics
  // Styling
  h_tpl_raw_strict->SetLineColor(kBlack);    h_tpl_raw_strict->SetLineWidth(2);
  h_tpl_bkg_relax ->SetLineColor(kOrange+7); h_tpl_bkg_relax ->SetLineStyle(7); h_tpl_bkg_relax->SetLineWidth(2);
  h_tpl_sub_relax ->SetLineColor(kMagenta+2);h_tpl_sub_relax ->SetLineWidth(2);

  auto __minbin = [](TH1* h){ double m=1e300; for(int b=1;b<=h->GetNbinsX(); ++b) m = std::min(m, h->GetBinContent(b)); return m; };
  const double __ymax = std::max({ h_tpl_raw_strict->GetMaximum(), h_tpl_bkg_relax->GetMaximum(), h_tpl_sub_relax->GetMaximum() });
  const double __ymin = std::min({ __minbin(h_tpl_raw_strict), __minbin(h_tpl_bkg_relax), __minbin(h_tpl_sub_relax), 0.0 });

  TCanvas cMM_Tpl_breakdown_relax("cMM_Tpl_breakdown_relax","Timing-template (RELAXED): raw / bkg / sub",1000,750);
  cMM_Tpl_breakdown_relax.SetRightMargin(0.08); cMM_Tpl_breakdown_relax.SetTopMargin(0.08);

  h_tpl_raw_strict->SetTitle(Form("Timing-template (RELAXED, pad=%.2f ns): raw vs. background vs. subtracted", TPL_EDGE_PAD_NS));
  h_tpl_raw_strict->GetXaxis()->SetTitle("M_{X} (GeV)");
  h_tpl_raw_strict->GetYaxis()->SetTitle("Counts");
  h_tpl_raw_strict->SetMaximum(1.15*std::max(1e-12, __ymax));
  h_tpl_raw_strict->SetMinimum(__ymin < 0 ? 1.15*__ymin : 0.0);

  h_tpl_raw_strict->Draw("HIST");
  h_tpl_bkg_relax ->Draw("HIST SAME");
  h_tpl_sub_relax ->Draw("HIST SAME");

  auto leg_rel = new TLegend(0.54,0.66,0.88,0.88);
  leg_rel->SetBorderSize(0); leg_rel->SetFillStyle(0);
  leg_rel->AddEntry(h_tpl_raw_strict, "Raw (CC)", "l");
  leg_rel->AddEntry(h_tpl_bkg_relax,  Form("Background (alpha_relax = %.4f)", __alpha_relax), "l");
  leg_rel->AddEntry(h_tpl_sub_relax,  "Subtracted = Raw - Bkg", "l");
  leg_rel->Draw();

  // Print diagnostics so you can see the scaling clearly
  printf("[tpl relaxed] Wsig=%.3f ns  Wacc=%.3f ns  alpha_strict=%.6f  alpha_relax=%.6f  k_relax=%.4f | "
         "Int: RAW=%g  BKG_strict=%g  BKG_relax=%g  SUB_relax=%g\n",
         __Wsig_strict, __Wacc, __alpha_strict, __alpha_relax, __k_relax,
         h_tpl_raw_strict->Integral(), h_tpl_bkg_strict->Integral(),
         h_tpl_bkg_relax->Integral(),  h_tpl_sub_relax->Integral());

  cMM_Tpl_breakdown_relax.Print(Form("%s_MM_breakdown_template_relaxed.pdf", __tpl_prefix2.Data()));
}
std::cout << "preTS: entries=" << hPairs_preTS.GetEntries()
          << " integral=" << hPairs_preTS.Integral() << "\n";
std::cout << "DT_pre: entries=" << hDT_Tbar_pre.GetEntries()
          << " integral=" << hDT_Tbar_pre.Integral() << "\n";
std::cout << "[Pairs] filled RAW=" << (long long)nPre
          << " POST=" << (long long)nPost << "\n";

// =====================================================================
// Timing-template OPTION A:
//   Signal: per-event CC (hMM_CC)
//   Background: from ALL-PAIRS sidebands (proper coeffs) scaled to event space
//   B_A(m) = (Wsig/Wacc)*[VERT_pairs + HORIZ_pairs] - (Wsig/Wacc)^2*[DIAG_pairs + PURE_pairs]
//   SUB_A(m) = CC_event(m) - B_A_event(m)
// =====================================================================
{
  TString __prefixA = argv[2]; __prefixA.ReplaceAll(".pdf","");

  // --- Geometry and coefficients ---
  double Wacc = 0.0;
  for (int i=0;i<NPOS;++i) Wacc += (posWins[i].hi - posWins[i].lo);
  for (int i=0;i<NNEG;++i) Wacc += (negWins[i].hi - negWins[i].lo);
  const double Wsig = (C_HI - C_LO);
  const double ACC  = Wsig * Wsig;

  const double aV = (Wacc > 0 ? Wsig / Wacc : 0.0);          // C×A
  const double aH = (Wacc > 0 ? Wsig / Wacc : 0.0);          // A×C
  const double aA = (Wacc > 0 ? (ACC / (Wacc*Wacc)) : 0.0);  // A×A total

  // --- Build background in PAIR space, then convert to EVENT space ---
  hTplA_bkg = (TH1*)hMM_CC_pairs.Clone("hTplA_bkg");
  hTplA_bkg->SetDirectory(nullptr); hTplA_bkg->Sumw2(); hTplA_bkg->Reset("ICESM");

  hTplA_bkg->Add(&hMM_vert_pairs , +aV);
  hTplA_bkg->Add(&hMM_horiz_pairs, +aH);
  hTplA_bkg->Add(&hMM_diag_pairs , -aA);
  hTplA_bkg->Add(&hMM_pure_pairs , -aA);

  // pairs → event normalization using central counts
  const double cc_evt  = hMM_CC.Integral();        // per-event CC (RAW for template)
  const double cc_pair = hMM_CC_pairs.Integral();  // all-pairs CC
  const double k_pairs2event = (cc_pair > 0 ? cc_evt / cc_pair : 1.0);
  hTplA_bkg->Scale(k_pairs2event);

  // --- Subtract from per-event CC ---
  hTplA_raw = (TH1*)hMM_CC.Clone("hTplA_raw"); hTplA_raw->SetDirectory(nullptr); hTplA_raw->Sumw2();
  hTplA_sub = (TH1*)hTplA_raw->Clone("hTplA_sub"); hTplA_sub->Add(hTplA_bkg, -1.0);

  // --- Plot breakdown ---
  hTplA_raw->SetLineColor(kBlack);     hTplA_raw->SetLineWidth(2);
  hTplA_bkg->SetLineColor(kOrange+7);  hTplA_bkg->SetLineStyle(7); hTplA_bkg->SetLineWidth(2);
  hTplA_sub->SetLineColor(kMagenta+2); hTplA_sub->SetLineWidth(2);

  auto __minbin = [](TH1* h){ double m=1e300; for(int b=1;b<=h->GetNbinsX();++b) m=std::min(m,h->GetBinContent(b)); return m; };
  const double __ymax = std::max({ hTplA_raw->GetMaximum(), hTplA_bkg->GetMaximum(), hTplA_sub->GetMaximum() });
  const double __ymin = std::min({ __minbin(hTplA_raw), __minbin(hTplA_bkg), __minbin(hTplA_sub), 0.0 });

  TCanvas cMM_TplA("cMM_TplA","Template Option A",1000,750);
  cMM_TplA.SetRightMargin(0.08); cMM_TplA.SetTopMargin(0.08);

  hTplA_raw->SetTitle("Timing-template (Option A): per-event CC vs all-pairs background vs subtracted");
  hTplA_raw->GetXaxis()->SetTitle("M_{X} (GeV)");
  hTplA_raw->GetYaxis()->SetTitle("Counts");
  hTplA_raw->SetMaximum(1.15*std::max(1e-12,__ymax));
  hTplA_raw->SetMinimum(__ymin < 0 ? 1.15*__ymin : 0.0);

  hTplA_raw->Draw("HIST");
  hTplA_bkg->Draw("HIST SAME");
  hTplA_sub->Draw("HIST SAME");

  auto legA = new TLegend(0.50,0.66,0.88,0.88);
  legA->SetBorderSize(0); legA->SetFillStyle(0);
  legA->AddEntry(hTplA_raw, "Raw (CC, per-event chosen pair)", "l");
  legA->AddEntry(hTplA_bkg, Form("Bkg from ALL-PAIRS (aV=aH=%.4f, aA=%.4f) × k=%.3g", aV, aA, k_pairs2event), "l");
  legA->AddEntry(hTplA_sub, "Subtracted = Raw - Bkg_pairs→event", "l");
  legA->Draw();

  const double int_raw = hTplA_raw->Integral();
  const double int_bkg = hTplA_bkg->Integral();
  const double int_sub = hTplA_sub->Integral();
  printf("[TPL Option A] Wsig=%.3f Wacc=%.3f  aV=aH=%.6f  aA=%.6f  |  CC_evt=%g  CC_pairs=%g  k=%g  |  Int: RAW=%g  BKG=%g  SUB=%g  f_acc_A=%.5f\n",
         Wsig,Wacc,aV,aH,aA, cc_evt,cc_pair,k_pairs2event, int_raw,int_bkg,int_sub, (int_raw>0?int_bkg/int_raw:0.0));

  cMM_TplA.Print(Form("%s_MM_breakdown_template_optionA.pdf", __prefixA.Data()));
}

// =====================================================================
// Option A: 2D timing map (pairs that feed B_pairs) + MM(A) breakdown
//   Left  = hPairs_usedA  (V + H + A×A contributing pairs)
//   Right = hTplA_*       (RAW/BKG/SUB for Option A)
// =====================================================================
{
  TString __prefixA = argv[2]; __prefixA.ReplaceAll(".pdf","");

  auto __min3 = [](double a,double b,double c){ return std::min(a,std::min(b,c)); };
  auto __max3 = [](double a,double b,double c){ return std::max(a,std::max(b,c)); };

  double tmin = __min3(C_LO, C_HI, C_LO - 1.0);
  double tmax = __max3(C_LO, C_HI, C_HI + 1.0);
  for (int i=0;i<NPOS;++i){ tmin = std::min(tmin, posWins[i].lo); tmax = std::max(tmax, posWins[i].hi); }
  for (int i=0;i<NNEG;++i){ tmin = std::min(tmin, negWins[i].lo); tmax = std::max(tmax, negWins[i].hi); }
  const double pad = 0.05*(tmax - tmin);
  tmin -= pad; tmax += pad;

  TCanvas cA_panel("cA_panel","Timing vs MM (Option A)",1200,600);
  cA_panel.Divide(2,1);

  // ---------------- Left: DATA used by Option A background ----------------
  cA_panel.cd(1);
  gPad->SetRightMargin(0.14); gPad->SetTopMargin(0.08);
  hPairs_usedA.SetTitle("Option A: contributing pairs (V + H + A#timesA); t_{i} vs t_{j}");
  hPairs_usedA.GetXaxis()->SetTitle("Cluster 1 time (ns)");
  hPairs_usedA.GetYaxis()->SetTitle("Cluster 2 time (ns)");
  hPairs_usedA.SetContour(99);
  hPairs_usedA.Draw("COLZ");

  // Overlays (same cosmetics as before)
  TBox *bCC = new TBox(C_LO, C_LO, C_HI, C_HI);
  bCC->SetLineColor(kGreen+2); bCC->SetLineWidth(2); bCC->SetFillStyle(0); bCC->Draw("SAME");
  for(int i=0;i<NPOS;++i) { TBox *b=new TBox(C_LO, posWins[i].lo, C_HI, posWins[i].hi);
    b->SetLineColor(kOrange+7); b->SetLineStyle(2); b->SetFillStyle(0); b->Draw("SAME"); }
  for(int i=0;i<NNEG;++i) { TBox *b=new TBox(C_LO, negWins[i].lo, C_HI, negWins[i].hi);
    b->SetLineColor(kOrange+7); b->SetLineStyle(2); b->SetFillStyle(0); b->Draw("SAME"); }
  for(int i=0;i<NPOS;++i) { TBox *b=new TBox(posWins[i].lo, C_LO, posWins[i].hi, C_HI);
    b->SetLineColor(kAzure+2);  b->SetLineStyle(2); b->SetFillStyle(0); b->Draw("SAME"); }
  for(int i=0;i<NNEG;++i) { TBox *b=new TBox(negWins[i].lo, C_LO, negWins[i].hi, C_HI);
    b->SetLineColor(kAzure+2);  b->SetLineStyle(2); b->SetFillStyle(0); b->Draw("SAME"); }
  gPad->RedrawAxis();

  // ---------------- Right: MM breakdown (Option A) ----------------
  cA_panel.cd(2);
  gPad->SetRightMargin(0.08); gPad->SetTopMargin(0.08);

  TH1* hAmax = (hTplA_raw ? (TH1*)hTplA_raw->Clone("hAmax") : nullptr);
  if (hAmax) {
    hAmax->SetLineColor(0); hAmax->SetLineWidth(0);
    double ymaxA = std::max(std::max(hTplA_raw->GetMaximum(), hTplA_bkg->GetMaximum()), hTplA_sub->GetMaximum());
    hAmax->SetTitle("Option A: RAW / BKG (pairs) / SUB");
    hAmax->GetXaxis()->SetTitle("M_{X} (GeV)"); hAmax->GetYaxis()->SetTitle("Counts");
    hAmax->SetMaximum(1.15*std::max(1e-12, ymaxA));
    hAmax->Draw("AXIS");
    hTplA_raw->Draw("HIST SAME"); hTplA_bkg->Draw("HIST SAME"); hTplA_sub->Draw("HIST SAME");
    auto leg = new TLegend(0.52,0.66,0.88,0.88); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(hTplA_raw, "RAW (CC per-event)", "l");
    leg->AddEntry(hTplA_bkg, "BKG from all-pairs (Option A)", "l");
    leg->AddEntry(hTplA_sub, "SUB = RAW - BKG", "l");
    leg->Draw();
  }

  cA_panel.Print(Form("%s_TimingPlusMM_optionA.pdf", __prefixA.Data()));
}

// =====================================================================
// Timing-template OPTION B (union per-event):
//   Signal: per-event CC (hMM_CC)
//   Background built from per-event unions:
//     H_side ≡ V∪H (one stripe entry/event)
//     H_AA   ≡ A×A (one corner entry/event)
//   B_B(m) = (Wsig/Wacc)*H_side(m) - (Wsig/Wacc)^2*H_AA(m)
//   SUB_B  = CC_event - B_B
// =====================================================================
{
  TString __prefixB = argv[2]; __prefixB.ReplaceAll(".pdf","");

  // Geometry & coeffs
  double Wacc = 0.0;
  for (int i=0;i<NPOS;++i) Wacc += (posWins[i].hi - posWins[i].lo);
  for (int i=0;i<NNEG;++i) Wacc += (negWins[i].hi - negWins[i].lo);
  const double Wsig = (C_HI - C_LO);
  const double ACC  = Wsig * Wsig;

  const double a1 = (Wacc > 0 ? Wsig / Wacc : 0.0);            // stripes
  const double aA = (Wacc > 0 ? (ACC / (Wacc*Wacc)) : 0.0);    // corners

  // Background from per-event unions
  hTplB_bkg = (TH1*)hMM_CC.Clone("hTplB_bkg");
  hTplB_bkg->SetDirectory(nullptr); hTplB_bkg->Sumw2(); hTplB_bkg->Reset("ICESM");
  hTplB_bkg->Add(hMM_side_evBest, +a1);
  hTplB_bkg->Add(hMM_AA_evBest,   -aA);

  // Subtract from per-event CC
  hTplB_raw = (TH1*)hMM_CC.Clone("hTplB_raw"); hTplB_raw->SetDirectory(nullptr); hTplB_raw->Sumw2();
  hTplB_sub = (TH1*)hTplB_raw->Clone("hTplB_sub"); hTplB_sub->Add(hTplB_bkg, -1.0);

  // Plot breakdown
  hTplB_raw->SetLineColor(kBlack);    hTplB_raw->SetLineWidth(2);
  hTplB_bkg->SetLineColor(kOrange+7); hTplB_bkg->SetLineStyle(7); hTplB_bkg->SetLineWidth(2);
  hTplB_sub->SetLineColor(kAzure+2);  hTplB_sub->SetLineWidth(2);

  auto __minbin = [](TH1* h){ double m=1e300; for(int b=1;b<=h->GetNbinsX();++b) m=std::min(m,h->GetBinContent(b)); return m; };
  const double __ymax = std::max({ hTplB_raw->GetMaximum(), hTplB_bkg->GetMaximum(), hTplB_sub->GetMaximum() });
  const double __ymin = std::min({ __minbin(hTplB_raw), __minbin(hTplB_bkg), __minbin(hTplB_sub), 0.0 });

  TCanvas cMM_TplB("cMM_TplB","Template Option B (union)",1000,750);
  cMM_TplB.SetRightMargin(0.08); cMM_TplB.SetTopMargin(0.08);

  hTplB_raw->SetTitle("Timing-template (Option B): per-event CC vs union-sidebands vs subtracted");
  hTplB_raw->GetXaxis()->SetTitle("M_{X} (GeV)");
  hTplB_raw->GetYaxis()->SetTitle("Counts");
  hTplB_raw->SetMaximum(1.15*std::max(1e-12,__ymax));
  hTplB_raw->SetMinimum(__ymin < 0 ? 1.15*__ymin : 0.0);

  hTplB_raw->Draw("HIST");
  hTplB_bkg->Draw("HIST SAME");
  hTplB_sub->Draw("HIST SAME");

  auto legB = new TLegend(0.50,0.66,0.88,0.88);
  legB->SetBorderSize(0); legB->SetFillStyle(0);
  legB->AddEntry(hTplB_raw, "Raw (CC, per-event chosen pair)", "l");
  legB->AddEntry(hTplB_bkg, Form("Bkg: a1*(V∪H) - aA*(A×A), a1=%.4f, aA=%.4f", a1, aA), "l");
  legB->AddEntry(hTplB_sub, "Subtracted = Raw - Bkg_union", "l");
  legB->Draw();

  const double int_raw = hTplB_raw->Integral();
  const double int_bkg = hTplB_bkg->Integral();
  const double int_sub = hTplB_sub->Integral();
  printf("[TPL Option B union] Wsig=%.3f Wacc=%.3f  a1=%.6f  aA=%.6f  |  Int: RAW=%g  BKG=%g  SUB=%g  f_acc_B=%.5f\n",
         Wsig, Wacc, a1, aA, int_raw, int_bkg, int_sub, (int_raw>0?int_bkg/int_raw:0.0));

  cMM_TplB.Print(Form("%s_MM_breakdown_template_optionB.pdf", __prefixB.Data()));
}

// =====================================================================
// Option B: 2D timing map (per-event unions actually used) + MM(B)
//   Left  = hPairs_usedB (ONE stripe + ONE corner per event; chosen pairs)
//   Right = hTplB_*      (RAW/BKG/SUB for Option B unions)
// =====================================================================
{
  TString __prefixB = argv[2]; __prefixB.ReplaceAll(".pdf","");

  auto __min3 = [](double a,double b,double c){ return std::min(a,std::min(b,c)); };
  auto __max3 = [](double a,double b,double c){ return std::max(a,std::max(b,c)); };

  double tmin = __min3(C_LO, C_HI, C_LO - 1.0);
  double tmax = __max3(C_LO, C_HI, C_HI + 1.0);
  for (int i=0;i<NPOS;++i){ tmin = std::min(tmin, posWins[i].lo); tmax = std::max(tmax, posWins[i].hi); }
  for (int i=0;i<NNEG;++i){ tmin = std::min(tmin, negWins[i].lo); tmax = std::max(tmax, negWins[i].hi); }
  const double pad = 0.05*(tmax - tmin);
  tmin -= pad; tmax += pad;

  TCanvas cB_panel("cB_panel","Timing vs MM (Option B)",1200,600);
  cB_panel.Divide(2,1);

  // ---------------- Left: DATA used by Option B background ----------------
  cB_panel.cd(1);
  gPad->SetRightMargin(0.14); gPad->SetTopMargin(0.08);
  hPairs_usedB.SetTitle("Option B: per-event chosen pairs (V∪H + A×A); t_{i} vs t_{j}");
  hPairs_usedB.GetXaxis()->SetTitle("Cluster 1 time (ns)");
  hPairs_usedB.GetYaxis()->SetTitle("Cluster 2 time (ns)");
  hPairs_usedB.SetContour(99);
  hPairs_usedB.Draw("COLZ");

  // Region overlays (same as before)
  TBox *bCCB = new TBox(C_LO, C_LO, C_HI, C_HI);
  bCCB->SetLineColor(kGreen+2); bCCB->SetLineWidth(2); bCCB->SetFillStyle(0); bCCB->Draw("SAME");
  for(int i=0;i<NPOS;++i){ TBox *b=new TBox(C_LO, posWins[i].lo, C_HI, posWins[i].hi);
    b->SetLineColor(kOrange+7); b->SetLineStyle(2); b->SetFillStyle(0); b->Draw("SAME"); }
  for(int i=0;i<NNEG;++i){ TBox *b=new TBox(C_LO, negWins[i].lo, C_HI, negWins[i].hi);
    b->SetLineColor(kOrange+7); b->SetLineStyle(2); b->SetFillStyle(0); b->Draw("SAME"); }
  for(int i=0;i<NPOS;++i){ TBox *b=new TBox(posWins[i].lo, C_LO, posWins[i].hi, C_HI);
    b->SetLineColor(kAzure+2);  b->SetLineStyle(2); b->SetFillStyle(0); b->Draw("SAME"); }
  for(int i=0;i<NNEG;++i){ TBox *b=new TBox(negWins[i].lo, C_LO, negWins[i].hi, C_HI);
    b->SetLineColor(kAzure+2);  b->SetLineStyle(2); b->SetFillStyle(0); b->Draw("SAME"); }
  gPad->RedrawAxis();

  // ---------------- Right: MM breakdown (Option B) ----------------
  cB_panel.cd(2);
  gPad->SetRightMargin(0.08); gPad->SetTopMargin(0.08);

  TH1* hBmax = (hTplB_raw ? (TH1*)hTplB_raw->Clone("hBmax") : nullptr);
  if (hBmax) {
    hBmax->SetLineColor(0); hBmax->SetLineWidth(0);
    double ymaxB = std::max(std::max(hTplB_raw->GetMaximum(), hTplB_bkg->GetMaximum()), hTplB_sub->GetMaximum());
    hBmax->SetTitle("Option B: RAW / BKG (per-event unions) / SUB");
    hBmax->GetXaxis()->SetTitle("M_{X} (GeV)"); hBmax->GetYaxis()->SetTitle("Counts");
    hBmax->SetMaximum(1.15*std::max(1e-12, ymaxB));
    hBmax->Draw("AXIS");
    hTplB_raw->Draw("HIST SAME"); hTplB_bkg->Draw("HIST SAME"); hTplB_sub->Draw("HIST SAME");
    auto leg = new TLegend(0.50,0.66,0.88,0.88); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(hTplB_raw, "RAW (CC, per-event)", "l");
    leg->AddEntry(hTplB_bkg, "BKG from per-event unions", "l");
    leg->AddEntry(hTplB_sub, "SUB = RAW - BKG", "l");
    leg->Draw();
  }

  cB_panel.Print(Form("%s_TimingPlusMM_optionB.pdf", __prefixB.Data()));
}

// (optional) one-line sanity print
printf("[MM breakdowns] top2: raw=%g bkg=%g sub=%g | all: raw=%g bkg=%g sub=%g | tpl: raw=%g bkg=%g sub=%g\n",
       hMM_top2_sig->Integral(), alpha*hMM_top2_acc->Integral(), ((TH1*)gDirectory->Get("h_top2_sub")) ? ((TH1*)gDirectory->Get("h_top2_sub"))->Integral() : 0.0,
       hMM_all_sig ->Integral(), alpha*hMM_all_acc ->Integral(), ((TH1*)gDirectory->Get("h_all_sub"))  ? ((TH1*)gDirectory->Get("h_all_sub")) ->Integral() : 0.0,
       hMM_CC.Integral(),        ((TH1*)gDirectory->Get("h_temp_bkg")) ? ((TH1*)gDirectory->Get("h_temp_bkg"))->Integral() : 0.0,
       hMM_subtracted.Integral());

    // -----------------------------------------------
    gSystem->Exit(0);
    return 0;
}