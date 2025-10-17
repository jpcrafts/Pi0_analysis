// Pi0_end2end_dummyFirst_Amethod.cxx
// End-to-end slim analysis:
//   1) Read DATA and DUMMY trees (v5 branch names)
//   2) HMS cuts
//   3) Build ordered photon pairs; compute MM and Mgg per pair
//   4) Fill timing categories: CC, V, H, AA (for DATA and DUMMY)
//      - For DUMMY, also fill Upstream (UP) / Downstream (DN) splits by y-sign
//   5) Dummy subtraction FIRST (per category) using v5 UD factors (8.467, 4.256), charge-normalized
//   6) Option-A timing-accidental subtraction on the dummy-subtracted categories
//   7) Write final, charge-normalized, dummy-subtracted & accidental-subtracted MM and Mgg
//   8) Produce simple overlay canvases
//
// Usage:
//   g++ -O0 -g -std=c++17 -o Pi_0_v5_slim Pi_0_v5_slim.cxx alglib_src/*.cpp -I. -Ialglib_src `root-config --cflags --libs` -lTMVA -lRooFit -lRooFitCore
//   ./Pi0_end2end data.root dummy.root out.root [Qdata] [Qdummy]
//   ./Pi_0_v5_slim /cache/hallc/c-nps/analysis/pass2/replays/production/nps_hms_coin_4205_0_1_-1.root VolatileROOTfiles/dummy_x58_q51_p5_merged.root 4205_v5_slim.root 29446.282 254792.874
//
// Notes:
// - Tree name is assumed "T" (same as v5). Edit if needed.
// - HMS & NPS branches copied from v5; adjust branch names if your files differ.
// - MM kinematics uses beam+proton target and HMS e' 4-vector; NPS photon dirs from (x,y,NPS_dist) rotated by NPS theta.
// - Mgg uses a simple small-angle estimate from cluster separations; replace with your full routine if desired.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip> // for std::fixed, std::setprecision

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TBox.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH2F.h"
// Enable ALL-PAIRS Template-A with pairs→event scaling (matches full v5)
#ifndef USE_PAIR_TEMPLATE
#define USE_PAIR_TEMPLATE 1 // set to 0 to switch the multiplicity on
#endif

// ───────────────────────── constants / config ─────────────────────────
static constexpr double mp = 0.938272081;     // GeV
static constexpr double e0_nom = 10.54350201; // GeV beam energy (tune per run if needed)

static constexpr double NPS_theta_deg = -13.43; // deg
static constexpr double NPS_theta_rad = NPS_theta_deg * M_PI / 180.0;
static constexpr double NPS_dist_cm = 407.0; // cm

struct Win
{
  double lo, hi;
};
static constexpr double C_LO = 149.0; // ns
static constexpr double C_HI = 151.0; // ns
static const Win posWins[] = {{153.0, 155.0}, {155.0, 157.0}, {157.0, 159.0}};
static const Win negWins[] = {{145.0, 147.0}, {143.0, 145.0}, {141.0, 143.0}};
static const int NPOS = int(sizeof(posWins) / sizeof(posWins[0]));
static const int NNEG = int(sizeof(negWins) / sizeof(negWins[0]));

// Binning
static constexpr int nMM = 300;
static constexpr double mmLo = 0.0, mmHi = 5.0;
static constexpr int nMG = 200;
static constexpr double mgLo = 0.0, mgHi = 1.0;

// v5 timing dummy normalization factors
static constexpr double kUP_v5 = 8.467; // upstream
static constexpr double kDN_v5 = 4.256; // downstream

// ─────────────────────────── helpers ───────────────────────────

// Fiducial gate (adjust numbers to your full script if needed)
inline bool goodXY(double x, double y)
{
  return (x > -29.16 && x < 29.16 && y > -35.64 && y < 35.64);
}

inline void rotY_passive(double x, double y, double z, double deg,
                         double &xo, double &yo, double &zo)
{
  double th = deg * M_PI / 180.0, c = std::cos(th), s = std::sin(th);
  xo = c * x - s * z;
  yo = y;
  zo = s * x + c * z;
}

inline bool passHMSCuts(double edt, double dp, double et, double npe, double th, double ph)
{
  return (edt < 0.1 && std::fabs(dp) <= 8.5 && et > 0.6 && npe > 1.0 &&
          std::fabs(th) <= 0.09 && std::fabs(ph) <= 0.09);
}

inline bool passHMSCuts(double edt, double dp, double et, double npe, double th, double ph,
                        int ncl, const double *cE, const double *cX, const double *cY)
{
  // First, require the original HMS cuts
  if (!(edt < 0.1 && std::fabs(dp) <= 8.5 && et > 0.6 && npe > 1.0 &&
        std::fabs(th) <= 0.09 && std::fabs(ph) <= 0.09))
  {
    return false;
  }

  // Then, require at least TWO fiducial NPS clusters with basic energy quality
  int nFid = 0;
  for (int i = 0; i < ncl; ++i)
  {
    if (cE[i] >= 0.6 && goodXY(cX[i], cY[i]))
      ++nFid;
    if (nFid >= 2)
      return true; // early exit
  }
  return false;
}

static inline bool inWin(double t, const Win &w) { return (t >= w.lo && t <= w.hi); }
static inline bool inAny(double t, const Win *arr, int n)
{
  for (int i = 0; i < n; ++i)
    if (inWin(t, arr[i]))
      return true;
  return false;
}
static inline bool isCC(double ti, double tj) { return inWin(ti, {C_LO, C_HI}) && inWin(tj, {C_LO, C_HI}); }
static inline bool isV(double ti, double tj) { return inWin(ti, {C_LO, C_HI}) && (inAny(tj, posWins, NPOS) || inAny(tj, negWins, NNEG)); }
static inline bool isH(double ti, double tj) { return (inAny(ti, posWins, NPOS) || inAny(ti, negWins, NNEG)) && inWin(tj, {C_LO, C_HI}); }
static inline bool isAA(double ti, double tj)
{
  bool Ai = (inAny(ti, posWins, NPOS) || inAny(ti, negWins, NNEG));
  bool Aj = (inAny(tj, posWins, NPOS) || inAny(tj, negWins, NNEG));
  return (Ai && Aj);
}

static inline bool computeMM(double E1, double E2,
                             double x1, double y1, double x2, double y2, // cm at NPS plane
                             double Ep, double epx, double epy, double epz,
                             double &mm_out)
{
  // photon directions from (x,y,NPS_dist)
  double px1_h, py1_h, pz1_h, px2_h, py2_h, pz2_h;
  rotY_passive(x1, y1, NPS_dist_cm, NPS_theta_deg, px1_h, py1_h, pz1_h);
  rotY_passive(x2, y2, NPS_dist_cm, NPS_theta_deg, px2_h, py2_h, pz2_h);
  const double r1 = std::sqrt(px1_h * px1_h + py1_h * py1_h + pz1_h * pz1_h);
  const double r2 = std::sqrt(px2_h * px2_h + py2_h * py2_h + pz2_h * pz2_h);
  if (!(r1 > 0 && r2 > 0))
    return false;
  const double u1x = px1_h / r1, u1y = py1_h / r1, u1z = pz1_h / r1;
  const double u2x = px2_h / r2, u2y = py2_h / r2, u2z = pz2_h / r2;

  // beam + target-at-rest
  const double Ein = e0_nom + mp;
  const double Pinx = 0.0, Piny = 0.0, Pinz = e0_nom;

  // photons (massless)
  const double p1x = E1 * u1x, p1y = E1 * u1y, p1z = E1 * u1z;
  const double p2x = E2 * u2x, p2y = E2 * u2y, p2z = E2 * u2z;

  const double E_out = Ep + E1 + E2;
  const double px_out = epx + p1x + p2x;
  const double py_out = epy + p1y + p2y;
  const double pz_out = epz + p1z + p2z;

  const double mm2 = std::pow(Ein - E_out, 2) - std::pow(Pinx - px_out, 2) - std::pow(Piny - py_out, 2) - std::pow(Pinz - pz_out, 2);
  if (!(mm2 > 0) || !std::isfinite(mm2))
    return false;
  mm_out = std::sqrt(mm2);
  return true;
}

static inline bool computeMgg(double E1, double E2,
                              double x1, double y1, double x2, double y2,
                              double &mgg_out)
{
  // small-angle approx from separation at NPS plane
  const double dx = (x1 - x2);
  const double dy = (y1 - y2);
  const double d = std::sqrt(dx * dx + dy * dy);
  const double theta12 = std::atan2(d, NPS_dist_cm);
  const double cos12 = std::cos(theta12);
  const double m2 = 2.0 * E1 * E2 * (1.0 - cos12);
  if (!(m2 > 0) || !std::isfinite(m2))
    return false;
  mgg_out = std::sqrt(m2);
  return true;
}

static inline bool inPos(double t) { return inAny(t, posWins, NPOS); }
static inline bool inNeg(double t) { return inAny(t, negWins, NNEG); }

// Both photons in sidebands on the SAME side (pos–pos or neg–neg)
static inline int whichPosIdx(double t)
{
  for (int k = 0; k < NPOS; ++k)
    if (inWin(t, posWins[k]))
      return k;
  return -1;
}
static inline int whichNegIdx(double t)
{
  for (int k = 0; k < NNEG; ++k)
    if (inWin(t, negWins[k]))
      return k;
  return -1;
}

static inline bool isAA_diag(double ti, double tj)
{
  const int ip = whichPosIdx(ti), jp = whichPosIdx(tj);
  if (ip >= 0 && jp >= 0)
    return (ip == jp); // upper-right: same positive stripe
  const int in = whichNegIdx(ti), jn = whichNegIdx(tj);
  if (in >= 0 && jn >= 0)
    return (in == jn); // lower-left: same negative stripe
  return false;
}

// Photons in sidebands on OPPOSITE sides (pos–neg or neg–pos)
static inline bool isAA_pure(double ti, double tj)
{
  return ((inPos(ti) && inNeg(tj)) || (inNeg(ti) && inPos(tj)));
}

// ─────────────────────── histogram pack ───────────────────────
struct Pack
{
  // Data categories (raw)
  TH1F hMM_CC, hMM_V, hMM_H, hMM_AA;
  TH1F hMG_CC, hMG_V, hMG_H, hMG_AA; // Mgg categories (raw)

  // Dummy UP/DN splits (only filled for dummy)
  TH1F hMM_CC_UP, hMM_CC_DN, hMM_V_UP, hMM_V_DN, hMM_H_UP, hMM_H_DN, hMM_AA_UP, hMM_AA_DN;
  TH1F hMG_CC_UP, hMG_CC_DN, hMG_V_UP, hMG_V_DN, hMG_H_UP, hMG_H_DN, hMG_AA_UP, hMG_AA_DN;

  // Derived (after dummy then A-method)
  TH1F hMM_Best, hMM_Sub; // MM B_est and SUB
  TH1F hMG_Best, hMG_Sub; // Mgg B_est and SUB

  // QA overlays: CC after dummy (the "before subtraction" curve)
  TH1F hMM_CC_afterDummy;
  TH1F hMG_CC_afterDummy;

  // Per-event CC (best-pair) histos (for pairs→event scaling and final subtraction)
  TH1F hMM_CC_evt, hMG_CC_evt;       // per-event CC (one entry per event)
  TH1F hMM_CC_evt_UP, hMM_CC_evt_DN; // per-event CC UP/DN (dummy-only fills)
  TH1F hMG_CC_evt_UP, hMG_CC_evt_DN; // per-event Mgg CC UP/DN (dummy-only)

  // MM (pairs)
  TH1F hMM_AD;               // AA diagonal (pos–pos + neg–neg)
  TH1F hMM_AP;               // AA pure/opposite (pos–neg + neg–pos)
  TH1F hMM_AD_UP, hMM_AD_DN; // dummy splits
  TH1F hMM_AP_UP, hMM_AP_DN;

  // Mγγ (pairs)
  TH1F hMG_AD, hMG_AP;
  TH1F hMG_AD_UP, hMG_AD_DN;
  TH1F hMG_AP_UP, hMG_AP_DN;

  // 2D timing QA (pairs): t_i vs t_j
  TH2F hTT_pairs;
  TH2F hTT_pairs_UP, hTT_pairs_DN; // for optional dummy subtraction

  Pack(const char *tag)
      : hMM_CC(TString::Format("hMM_CC_%s", tag), "MM CC;M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMM_V(TString::Format("hMM_V_%s", tag), "MM V;M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMM_H(TString::Format("hMM_H_%s", tag), "MM H;M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMM_AA(TString::Format("hMM_AA_%s", tag), "MM AA;M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMG_CC(TString::Format("hMG_CC_%s", tag), "M_{#gamma#gamma} CC;M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi), hMG_V(TString::Format("hMG_V_%s", tag), "M_{#gamma#gamma} V;M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi), hMG_H(TString::Format("hMG_H_%s", tag), "M_{#gamma#gamma} H;M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi), hMG_AA(TString::Format("hMG_AA_%s", tag), "M_{#gamma#gamma} AA;M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi), hMM_CC_afterDummy(TString::Format("hMM_CC_afterDummy_%s", tag),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "MM CC (after dummy);M_{X} (GeV);Counts", nMM, mmLo, mmHi),
        hMG_CC_afterDummy(TString::Format("hMG_CC_afterDummy_%s", tag),
                          "M_{#gamma#gamma} CC (after dummy);M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi)

        ,
        hMM_CC_evt(TString::Format("hMM_CC_evt_%s", tag), "MM CC (per-event);M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMG_CC_evt(TString::Format("hMG_CC_evt_%s", tag), "M_{#gamma#gamma} CC (per-event);M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi), hMM_CC_evt_UP(TString::Format("hMM_CC_evt_UP_%s", tag), "MM CC evt UP", nMM, mmLo, mmHi), hMM_CC_evt_DN(TString::Format("hMM_CC_evt_DN_%s", tag), "MM CC evt DN", nMM, mmLo, mmHi), hMG_CC_evt_UP(TString::Format("hMG_CC_evt_UP_%s", tag), "Mgg CC evt UP", nMG, mgLo, mgHi), hMG_CC_evt_DN(TString::Format("hMG_CC_evt_DN_%s", tag), "Mgg CC evt DN", nMG, mgLo, mgHi), hMM_CC_UP(TString::Format("hMM_CC_UP_%s", tag), "MM CC UP", nMM, mmLo, mmHi), hMM_CC_DN(TString::Format("hMM_CC_DN_%s", tag), "MM CC DN", nMM, mmLo, mmHi), hMM_V_UP(TString::Format("hMM_V_UP_%s", tag), "MM V UP", nMM, mmLo, mmHi), hMM_V_DN(TString::Format("hMM_V_DN_%s", tag), "MM V DN", nMM, mmLo, mmHi), hMM_H_UP(TString::Format("hMM_H_UP_%s", tag), "MM H UP", nMM, mmLo, mmHi), hMM_H_DN(TString::Format("hMM_H_DN_%s", tag), "MM H DN", nMM, mmLo, mmHi), hMM_AA_UP(TString::Format("hMM_AA_UP_%s", tag), "MM AA UP", nMM, mmLo, mmHi), hMM_AA_DN(TString::Format("hMM_AA_DN_%s", tag), "MM AA DN", nMM, mmLo, mmHi)

        ,
        hMG_CC_UP(TString::Format("hMG_CC_UP_%s", tag), "Mgg CC UP", nMG, mgLo, mgHi), hMG_CC_DN(TString::Format("hMG_CC_DN_%s", tag), "Mgg CC DN", nMG, mgLo, mgHi), hMG_V_UP(TString::Format("hMG_V_UP_%s", tag), "Mgg V UP", nMG, mgLo, mgHi), hMG_V_DN(TString::Format("hMG_V_DN_%s", tag), "Mgg V DN", nMG, mgLo, mgHi), hMG_H_UP(TString::Format("hMG_H_UP_%s", tag), "Mgg H UP", nMG, mgLo, mgHi), hMG_H_DN(TString::Format("hMG_H_DN_%s", tag), "Mgg H DN", nMG, mgLo, mgHi), hMG_AA_UP(TString::Format("hMG_AA_UP_%s", tag), "Mgg AA UP", nMG, mgLo, mgHi), hMG_AA_DN(TString::Format("hMG_AA_DN_%s", tag), "Mgg AA DN", nMG, mgLo, mgHi)

        ,
        hMM_Best(TString::Format("hMM_Best_%s", tag), "MM B_{est};M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMM_Sub(TString::Format("hMM_Sub_%s", tag), "MM SUB;M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMG_Best(TString::Format("hMG_Best_%s", tag), "Mgg B_{est};M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi), hMG_Sub(TString::Format("hMG_Sub_%s", tag), "Mgg SUB;M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi)

        ,
        hMM_AD(Form("hMM_AD_%s", tag), "MM AA (diag);M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMM_AP(Form("hMM_AP_%s", tag), "MM AA (pure);M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMM_AD_UP(Form("hMM_AD_UP_%s", tag), "MM AA (diag) UP;M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMM_AD_DN(Form("hMM_AD_DN_%s", tag), "MM AA (diag) DN;M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMM_AP_UP(Form("hMM_AP_UP_%s", tag), "MM AA (pure) UP;M_{X} (GeV);Counts", nMM, mmLo, mmHi), hMM_AP_DN(Form("hMM_AP_DN_%s", tag), "MM AA (pure) DN;M_{X} (GeV);Counts", nMM, mmLo, mmHi)

        ,
        hMG_AD(Form("hMG_AD_%s", tag), "m_{#gamma#gamma} AA (diag);m_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi), hMG_AP(Form("hMG_AP_%s", tag), "m_{#gamma#gamma} AA (pure);m_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi), hMG_AD_UP(Form("hMG_AD_UP_%s", tag), "m_{#gamma#gamma} AA (diag) UP;...", nMG, mgLo, mgHi), hMG_AD_DN(Form("hMG_AD_DN_%s", tag), "m_{#gamma#gamma} AA (diag) DN;...", nMG, mgLo, mgHi), hMG_AP_UP(Form("hMG_AP_UP_%s", tag), "m_{#gamma#gamma} AA (pure) UP;...", nMG, mgLo, mgHi), hMG_AP_DN(Form("hMG_AP_DN_%s", tag), "m_{#gamma#gamma} AA (pure) DN;...", nMG, mgLo, mgHi), hTT_pairs(Form("hTT_pairs_%s", tag), "t_{i} vs t_{j};t_{i} (ns);t_{j} (ns)", 200, 140, 160, 200, 140, 160), hTT_pairs_UP(Form("hTT_pairs_UP_%s", tag), "t_{i} vs t_{j} (UP)", 200, 140, 160, 200, 140, 160), hTT_pairs_DN(Form("hTT_pairs_DN_%s", tag), "t_{i} vs t_{j} (DN)", 200, 140, 160, 200, 140, 160)
  {
    auto s2 = [&](TH1F &h)
    { h.Sumw2(); };
    s2(hMM_CC);
    s2(hMM_V);
    s2(hMM_H);
    s2(hMM_AA);
    s2(hMG_CC);
    s2(hMG_V);
    s2(hMG_H);
    s2(hMG_AA);
    s2(hMM_CC_UP);
    s2(hMM_CC_DN);
    s2(hMM_V_UP);
    s2(hMM_V_DN);
    s2(hMM_H_UP);
    s2(hMM_H_DN);
    s2(hMM_AA_UP);
    s2(hMM_AA_DN);
    s2(hMG_CC_UP);
    s2(hMG_CC_DN);
    s2(hMG_V_UP);
    s2(hMG_V_DN);
    s2(hMG_H_UP);
    s2(hMG_H_DN);
    s2(hMG_AA_UP);
    s2(hMG_AA_DN);
    s2(hMM_Best);
    s2(hMM_Sub);
    s2(hMG_Best);
    s2(hMG_Sub);
  }
};

// ─────────────────────── file processing ───────────────────────
static void fillFromFile(const TString &inF, const char *tag, Pack &O)
{
  TFile f(inF, "READ");
  if (f.IsZombie())
  {
    std::cerr << "Cannot open " << inF << "\n";
    return;
  }
  TTree *tr = dynamic_cast<TTree *>(f.Get("T"));
  if (!tr)
  {
    std::cerr << "No T tree in " << inF << "\n";
    return;
  }

  // HMS branches (v5 names)
  double edt = 0, dp = 0, et = 0, npe = 0, th = 0, ph = 0;
  double hp = 0, hpx = 0, hpy = 0, hpz = 0;
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
  tr->SetBranchStatus("H.gtr.p", 1);
  tr->SetBranchAddress("H.gtr.p", &hp);
  tr->SetBranchStatus("H.gtr.px", 1);
  tr->SetBranchAddress("H.gtr.px", &hpx);
  tr->SetBranchStatus("H.gtr.py", 1);
  tr->SetBranchAddress("H.gtr.py", &hpy);
  tr->SetBranchStatus("H.gtr.pz", 1);
  tr->SetBranchAddress("H.gtr.pz", &hpz);

  // NPS clusters (arrays)
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

  const bool isDummy = (std::string(tag) == "dummy");

  const Long64_t N = tr->GetEntries();
  for (Long64_t ie = 0; ie < N; ++ie)
  {
    // Track best CC pair per event (closest to 150 ns)
    const double t0 = 0.5 * (C_LO + C_HI);
    bool haveBestCC = false;
    double bestScore = 1e99;
    double bestMM = 0.0;
    double bestMG = -1.0;
    double bestYavg = 0.0;

    tr->GetEntry(ie);
    if (!passHMSCuts(edt, dp, et, npe, th, ph, (int)ncl, cE, cX, cY))
      continue;

    for (int a = 0; a < (int)ncl; ++a)
    {
      // loose cluster quality similar to v5 "good cluster" (t window here is not a cut for categories)
      if (cE[a] < 0.6)
        continue;
      for (int b = 0; b < (int)ncl; ++b)
      {
        if (b == a)
          continue;
        if (cE[b] < 0.6)
          continue;
        const double ti = cT[a], tj = cT[b];

        double mm = 0.0;
        if (!computeMM(cE[a], cE[b], cX[a], cY[a], cX[b], cY[b], hp, hpx, hpy, hpz, mm))
          continue;
        double mg = 0.0;
        (void)computeMgg(cE[a], cE[b], cX[a], cY[a], cX[b], cY[b], mg); // ok if false

        const double yavg = 0.5 * (cY[a] + cY[b]); // sign for UP/DN split

        O.hTT_pairs.Fill(ti, tj);
        if (isDummy)
        {
          (yavg >= 0 ? O.hTT_pairs_DN : O.hTT_pairs_UP).Fill(ti, tj);
        }

        if (isCC(ti, tj))
        {
          O.hMM_CC.Fill(mm);
          if (mg > 0)
            O.hMG_CC.Fill(mg);
          if (isDummy)
          {
            (yavg >= 0 ? O.hMM_CC_DN : O.hMM_CC_UP).Fill(mm);
            if (mg > 0)
              (yavg >= 0 ? O.hMG_CC_DN : O.hMG_CC_UP).Fill(mg);
          }

          // update best-per-event CC candidate
          const double score = std::fabs(ti - t0) + std::fabs(tj - t0);
          if (score < bestScore)
          {
            haveBestCC = true;
            bestScore = score;
            bestMM = mm;
            bestMG = mg;
            bestYavg = yavg;
          }
        }
        else if (isV(ti, tj))
        {
          O.hMM_V.Fill(mm);
          if (mg > 0)
            O.hMG_V.Fill(mg);
          if (isDummy)
          {
            (yavg >= 0 ? O.hMM_V_DN : O.hMM_V_UP).Fill(mm);
            if (mg > 0)
              (yavg >= 0 ? O.hMG_V_DN : O.hMG_V_UP).Fill(mg);
          }
        }
        else if (isH(ti, tj))
        {
          O.hMM_H.Fill(mm);
          if (mg > 0)
            O.hMG_H.Fill(mg);
          if (isDummy)
          {
            (yavg >= 0 ? O.hMM_H_DN : O.hMM_H_UP).Fill(mm);
            if (mg > 0)
              (yavg >= 0 ? O.hMG_H_DN : O.hMG_H_UP).Fill(mg);
          }
        }
        else if (isAA_diag(ti, tj))
        {
          O.hMM_AD.Fill(mm);
          if (mg > 0)
            O.hMG_AD.Fill(mg);
          if (isDummy)
          {
            (yavg >= 0 ? O.hMM_AD_DN : O.hMM_AD_UP).Fill(mm);
            if (mg > 0)
              (yavg >= 0 ? O.hMG_AD_DN : O.hMG_AD_UP).Fill(mg);
          }
        }
        else if (isAA_pure(ti, tj))
        {
          O.hMM_AP.Fill(mm);
          if (mg > 0)
            O.hMG_AP.Fill(mg);
          if (isDummy)
          {
            (yavg >= 0 ? O.hMM_AP_DN : O.hMM_AP_UP).Fill(mm);
            if (mg > 0)
              (yavg >= 0 ? O.hMG_AP_DN : O.hMG_AP_UP).Fill(mg);
          }
        }
      }
    }

    // per-event CC fill using selected best pair
    if (haveBestCC)
    {
      O.hMM_CC_evt.Fill(bestMM);
      if (bestMG > 0)
        O.hMG_CC_evt.Fill(bestMG);
      if (isDummy)
      {
        if (bestYavg >= 0)
        {
          O.hMM_CC_evt_DN.Fill(bestMM);
          if (bestMG > 0)
            O.hMG_CC_evt_DN.Fill(bestMG);
        }
        else
        {
          O.hMM_CC_evt_UP.Fill(bestMM);
          if (bestMG > 0)
            O.hMG_CC_evt_UP.Fill(bestMG);
        }
      }
    }
  }
}

// ───────────── dummy-first then A-method (for both MM and Mgg) ─────────────
static void doDummyThenA(const Pack &D, const Pack &M, double Qdata, double Qdum,
                         double kUP, double kDN, Pack &OUT)
{
  auto norm_dummy = [&](const TH1F &hUP, const TH1F &hDN) -> TH1F
  {
    TH1F h = hUP;
    h.Reset("ICESM");
    if (Qdum > 0)
    {
      TH1F hup = hUP;
      hup.Scale(1.0 / (Qdum * kUP));
      TH1F hdn = hDN;
      hdn.Scale(1.0 / (Qdum * kDN));
      h.Add(&hup, 1.0);
      h.Add(&hdn, 1.0); // per µC
      h.Scale(Qdata);   // to data exposure
    }
    return h;
  };

  // MM categories after dummy
  TH1F hCC_mm = D.hMM_CC;
  TH1F hCCm = norm_dummy(M.hMM_CC_UP, M.hMM_CC_DN);
  hCC_mm.Add(&hCCm, -1.0);
  TH1F hV_mm = D.hMM_V;
  TH1F hVm = norm_dummy(M.hMM_V_UP, M.hMM_V_DN);
  hV_mm.Add(&hVm, -1.0);
  TH1F hH_mm = D.hMM_H;
  TH1F hHm = norm_dummy(M.hMM_H_UP, M.hMM_H_DN);
  hH_mm.Add(&hHm, -1.0);
  TH1F hAD_mm = D.hMM_AD;
  {
    TH1F tmp = norm_dummy(M.hMM_AD_UP, M.hMM_AD_DN);
    hAD_mm.Add(&tmp, -1.0);
  }
  TH1F hAP_mm = D.hMM_AP;
  {
    TH1F tmp = norm_dummy(M.hMM_AP_UP, M.hMM_AP_DN);
    hAP_mm.Add(&tmp, -1.0);
  }
  // Rebuild combined AA for QA output/writes:
  TH1F hAA_mm = hAD_mm;
  hAA_mm.Add(&hAP_mm, +1.0);

  // Mgg categories after dummy
  TH1F hCC_mg = D.hMG_CC;
  TH1F hCCg = norm_dummy(M.hMG_CC_UP, M.hMG_CC_DN);
  hCC_mg.Add(&hCCg, -1.0);
  TH1F hV_mg = D.hMG_V;
  TH1F hVg = norm_dummy(M.hMG_V_UP, M.hMG_V_DN);
  hV_mg.Add(&hVg, -1.0);
  TH1F hH_mg = D.hMG_H;
  TH1F hHg = norm_dummy(M.hMG_H_UP, M.hMG_H_DN);
  hH_mg.Add(&hHg, -1.0);
  TH1F hAD_mg = D.hMG_AD;
  {
    TH1F tmp = norm_dummy(M.hMG_AD_UP, M.hMG_AD_DN);
    hAD_mg.Add(&tmp, -1.0);
  }
  TH1F hAP_mg = D.hMG_AP;
  {
    TH1F tmp = norm_dummy(M.hMG_AP_UP, M.hMG_AP_DN);
    hAP_mg.Add(&tmp, -1.0);
  }
  // Rebuild combined AA for QA output/writes:
  TH1F hAA_mg = hAD_mg;
  hAA_mg.Add(&hAP_mg, +1.0);

  // NEW: area-aware coefficient for AD (diag-squares) + AP (full off-diagonal)
  double Wpos = 0.0, Wneg = 0.0, sumSqPos = 0.0, sumSqNeg = 0.0;
  for (int i = 0; i < NPOS; ++i)
  {
    const double w = posWins[i].hi - posWins[i].lo;
    Wpos += w;
    sumSqPos += w * w;
  }
  for (int i = 0; i < NNEG; ++i)
  {
    const double w = negWins[i].hi - negWins[i].lo;
    Wneg += w;
    sumSqNeg += w * w;
  }
  // Areas in (ti,tj) plane:
  const double A_AD = sumSqPos + sumSqNeg;                 // diag “same-stripe” squares (LL and UR only)
  const double A_AP = 2.0 * Wpos * Wneg;                   // full off-diagonal corners (UL + LR)
                                                           // Central (CC) width and area
  const double Wsig = (C_LO < C_HI) ? (C_HI - C_LO) : 0.0; // ~2 ns
  const double ACC = Wsig * Wsig;                          // CC box area

  // Stripe scale factors for Template A
  // V and H stripes each cover: Wsig × (Wpos + Wneg)  → aV = aH = ACC / [Wsig*(Wpos+Wneg)] = Wsig/(Wpos+Wneg)
  const double aV = (Wpos + Wneg) > 0.0 ? (Wsig / (Wpos + Wneg)) : 0.0;
  const double aH = aV;

  // AA is split into: AD (diag narrow squares only) and AP (opposite-side big corners).
  // We are *using only* AD + AP (not the same-side off-diagonal squares), so normalize to (A_AD + A_AP).
  const double aA = (A_AD + A_AP) > 0.0 ? (ACC / (A_AD + A_AP)) : 0.0;

  // Use the same coefficient on both pieces so that
  //   aAD*AD + aAP*AP  ≡  aA*(AD+AP)  (identical to classic AA when recombined).
  const double aA_split = aA;

  // (optional) debug print
  std::cout << std::fixed << std::setprecision(3)
            << "[geom] Wsig=" << Wsig
            << "  Wpos=" << Wpos << "  Wneg=" << Wneg
            << "  A_AD=" << A_AD << "  A_AP=" << A_AP
            << "  aV=" << aV << "  aH=" << aH
            << "  aA_total=" << aA << "  aA_split=" << aA_split << "\n";

#if USE_PAIR_TEMPLATE
  // ── Build B_est in ALL-PAIRS space and rescale to per-event CC (Template A) ──
  TH1F hMM_Best_pairs("hMM_Best_pairs", "MM B_{est} (ALL pairs);M_{X} (GeV);Counts", nMM, mmLo, mmHi);
  hMM_Best_pairs.Reset();
  hMM_Best_pairs.Add(&hV_mm, +aV);
  hMM_Best_pairs.Add(&hH_mm, +aH);
  hMM_Best_pairs.Add(&hAD_mm, -aA_split);
  hMM_Best_pairs.Add(&hAP_mm, -aA_split);

  TH1F hMG_Best_pairs("hMG_Best_pairs", "Mgg B_{est} (ALL pairs);M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi);
  hMG_Best_pairs.Reset();
  hMG_Best_pairs.Add(&hV_mg, +aV);
  hMG_Best_pairs.Add(&hH_mg, +aH);
  hMG_Best_pairs.Add(&hAD_mg, -aA_split);
  hMG_Best_pairs.Add(&hAP_mg, -aA_split);

  // Per-event CC after dummy (for scaling target)
  TH1F hCC_mm_evt = D.hMM_CC_evt;
  TH1F hCCm_evt = norm_dummy(M.hMM_CC_evt_UP, M.hMM_CC_evt_DN);
  hCC_mm_evt.Add(&hCCm_evt, -1.0);
  TH1F hCC_mg_evt = D.hMG_CC_evt;
  TH1F hCCg_evt = norm_dummy(M.hMG_CC_evt_UP, M.hMG_CC_evt_DN);
  hCC_mg_evt.Add(&hCCg_evt, -1.0);

  std::cerr << "[doDummyThenA] per-event CC built: "
            << "MM=" << hCC_mm_evt.Integral()
            << "  Mgg=" << hCC_mg_evt.Integral() << "\n";

  OUT.hMM_CC_afterDummy.Reset("ICESM");
  OUT.hMM_CC_afterDummy.Add(&hCC_mm_evt, 1.0);

  OUT.hMG_CC_afterDummy.Reset("ICESM");
  OUT.hMG_CC_afterDummy.Add(&hCC_mg_evt, 1.0);

  std::cerr << "[doDummyThenA] overlays seeded from per-event CC\n";

  // pairs→event scales (after dummy)
  const double cc_evt_mm = hCC_mm_evt.Integral();
  const double cc_pair_mm = hCC_mm.Integral();
  const double k_pairs2event_mm = (cc_pair_mm > 0 ? cc_evt_mm / cc_pair_mm : 1.0);

  const double cc_evt_mg = hCC_mg_evt.Integral();
  const double cc_pair_mg = hCC_mg.Integral();
  const double k_pairs2event_mg = (cc_pair_mg > 0 ? cc_evt_mg / cc_pair_mg : 1.0);

  // Scale to event space and subtract
  OUT.hMM_Best = hMM_Best_pairs;
  OUT.hMM_Best.Scale(k_pairs2event_mm);
  OUT.hMM_Sub = hCC_mm_evt;
  OUT.hMM_Sub.Add(&OUT.hMM_Best, -1.0);

  OUT.hMG_Best = hMG_Best_pairs;
  OUT.hMG_Best.Scale(k_pairs2event_mg);
  OUT.hMG_Sub = hCC_mg_evt;
  OUT.hMG_Sub.Add(&OUT.hMG_Best, -1.0);

  std::cerr << "[doDummyThenA] B_est/Sub done: "
            << "MM_Best=" << OUT.hMM_Best.Integral()
            << "  MM_Sub=" << OUT.hMM_Sub.Integral()
            << "  Mgg_Best=" << OUT.hMG_Best.Integral()
            << "  Mgg_Sub=" << OUT.hMG_Sub.Integral() << "\n";

#else
  // MM A-method on after-dummy categories (event-space analogue)
  OUT.hMM_Best.Reset();
  OUT.hMM_Best.Add(&hV_mm, +aV);
  OUT.hMM_Best.Add(&hH_mm, +aH);
  OUT.hMM_Best.Add(&hAD_mm, -aA_split);
  OUT.hMM_Best.Add(&hAP_mm, -aA_split);

  OUT.hMM_Sub.Reset();
  OUT.hMM_Sub.Add(&hCC_mm, 1.0);
  OUT.hMM_Sub.Add(&OUT.hMM_Best, -1.0);

  // Mgg A-method on after-dummy categories
  OUT.hMG_Best.Reset();
  OUT.hMG_Best.Add(&hV_mg, +aV);
  OUT.hMG_Best.Add(&hH_mg, +aH);
  OUT.hMG_Best.Add(&hAD_mg, -aA_split);
  OUT.hMG_Best.Add(&hAP_mg, -aA_split);

  OUT.hMG_Sub.Reset();
  OUT.hMG_Sub.Add(&hCC_mg, 1.0);
  OUT.hMG_Sub.Add(&OUT.hMG_Best, -1.0);

  OUT.hMM_CC_afterDummy.Write();
  OUT.hMG_CC_afterDummy.Write();

#endif
  // Optionally write after-dummy categories for QA
  hCC_mm.SetName("hMM_CC_afterDummy");
  hV_mm.SetName("hMM_V_afterDummy");
  hH_mm.SetName("hMM_H_afterDummy");
  hAA_mm.SetName("hMM_AA_afterDummy");
  hCC_mg.SetName("hMG_CC_afterDummy");
  hV_mg.SetName("hMG_V_afterDummy");
  hH_mg.SetName("hMG_H_afterDummy");
  hAA_mg.SetName("hMG_AA_afterDummy");
  hCC_mm.Write();
  hV_mm.Write();
  hH_mm.Write();
  hAA_mm.Write();
  hCC_mg.Write();
  hV_mg.Write();
  hH_mg.Write();
  hAA_mg.Write();
}

// ─────────────────────────────────── main ───────────────────────────────────
int main(int argc, char **argv)
{
  gROOT->SetBatch(kTRUE);
  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0] << " <data.root> <dummy.root> <out.root> [Qdata] [Qdummy]\n";
    return 1;
  }
  TString dataF = argv[1];
  TString dummyF = argv[2];
  TString outF = argv[3];

  double Qdata = 1.0, Qdum = 1.0;
  if (argc >= 5)
    Qdata = atof(argv[4]);
  if (argc >= 6)
    Qdum = atof(argv[5]);

  std::cout << "[Inputs] data=" << dataF << "  dummy=" << dummyF << "  out=" << outF << "\n";
  std::cout << "[Norm]   Q_data=" << Qdata << "  Q_dummy=" << Qdum
            << "  (kUP=" << kUP_v5 << ", kDN=" << kDN_v5 << ")\n";

  // Fill from files
  Pack D("data");
  fillFromFile(dataF, "data", D);
  Pack M("dummy");
  fillFromFile(dummyF, "dummy", M);

  // Open output FIRST so any internal Write() calls (e.g. in doDummyThenA) work
  TFile fout(outF, "RECREATE");
  fout.cd();

  // Do dummy-first then Option-A for both MM and Mgg
  Pack OUT("final");
  doDummyThenA(D, M, Qdata, Qdum, kUP_v5, kDN_v5, OUT);

  // Raw categories
  fout.cd();
  D.hMM_CC.Write();
  D.hMM_V.Write();
  D.hMM_H.Write();
  D.hMM_AA.Write();
  D.hMG_CC.Write();
  D.hMG_V.Write();
  D.hMG_H.Write();
  D.hMG_AA.Write();
  M.hMM_CC.Write();
  M.hMM_V.Write();
  M.hMM_H.Write();
  M.hMM_AA.Write();
  M.hMG_CC.Write();
  M.hMG_V.Write();
  M.hMG_H.Write();
  M.hMG_AA.Write();

  // Dummy splits
  M.hMM_CC_UP.Write();
  M.hMM_CC_DN.Write();
  M.hMM_V_UP.Write();
  M.hMM_V_DN.Write();
  M.hMM_H_UP.Write();
  M.hMM_H_DN.Write();
  M.hMM_AA_UP.Write();
  M.hMM_AA_DN.Write();
  M.hMG_CC_UP.Write();
  M.hMG_CC_DN.Write();
  M.hMG_V_UP.Write();
  M.hMG_V_DN.Write();
  M.hMG_H_UP.Write();
  M.hMG_H_DN.Write();
  M.hMG_AA_UP.Write();
  M.hMG_AA_DN.Write();

  // Final background & subtracted spectra
  OUT.hMM_Best.Write();
  OUT.hMM_Sub.Write();
  OUT.hMG_Best.Write();
  OUT.hMG_Sub.Write();

  // Simple overlays
  {
    double ymax = 0.0;
    ymax = std::max(ymax, OUT.hMM_Sub.GetMaximum());
    ymax = std::max(ymax, OUT.hMM_Best.GetMaximum());
    ymax = std::max(ymax, OUT.hMM_CC_afterDummy.GetMaximum());
    if (ymax <= 0)
      ymax = 1.0;

    TCanvas c("cMM", "MM (dummy-first, Option A)", 1000, 700);
    TH1F axis("axis", "Missing Mass;M_{X} (GeV);Counts", nMM, mmLo, mmHi);
    axis.SetMaximum(1.15 * ymax);
    axis.Draw("AXIS");

    OUT.hMM_CC_afterDummy.SetMarkerStyle(20); // NEW
    OUT.hMM_CC_afterDummy.SetMarkerSize(0.8); // NEW
    OUT.hMM_CC_afterDummy.SetLineWidth(2);    // NEW
    OUT.hMM_CC_afterDummy.Draw("E1 SAME");    // NEW
    OUT.hMM_Best.SetLineColor(kOrange + 7);
    OUT.hMM_Best.SetLineWidth(2);
    OUT.hMM_Best.Draw("HIST SAME");
    OUT.hMM_Sub.SetLineColor(kGreen + 2);
    OUT.hMM_Sub.SetLineWidth(3);
    OUT.hMM_Sub.Draw("HIST SAME");
    auto leg = new TLegend(0.58, 0.68, 0.90, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(&OUT.hMM_CC_afterDummy, "MM CC (after dummy)", "lep"); // NEW
    leg->AddEntry(&OUT.hMM_Best, "MM B_{est} (after dummy)", "l");
    leg->AddEntry(&OUT.hMM_Sub, "MM Final SUB", "l");
    leg->Draw();
    c.Write("MM_overlay_final");
  }
  {
    double ymax = 0.0;
    ymax = std::max(ymax, OUT.hMG_Sub.GetMaximum());
    ymax = std::max(ymax, OUT.hMG_Best.GetMaximum());
    ymax = std::max(ymax, OUT.hMG_CC_afterDummy.GetMaximum()); // NEW
    if (ymax <= 0)
      ymax = 1.0;

    TCanvas c("cMG", "Mgg (dummy-first, Option A)", 1000, 700);
    TH1F axis("axis", "M_{#gamma#gamma};M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi);
    axis.SetMaximum(1.15 * ymax);
    axis.Draw("AXIS");
    OUT.hMG_CC_afterDummy.SetMarkerStyle(20); // NEW
    OUT.hMG_CC_afterDummy.SetMarkerSize(0.8); // NEW
    OUT.hMG_CC_afterDummy.SetLineWidth(2);    // NEW
    OUT.hMG_CC_afterDummy.Draw("E1 SAME");    // NEW
    OUT.hMG_Best.SetLineColor(kOrange + 7);
    OUT.hMG_Best.SetLineWidth(2);
    OUT.hMG_Best.Draw("HIST SAME");
    OUT.hMG_Sub.SetLineColor(kGreen + 2);
    OUT.hMG_Sub.SetLineWidth(3);
    OUT.hMG_Sub.Draw("HIST SAME");
    auto leg = new TLegend(0.58, 0.68, 0.90, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(&OUT.hMG_CC_afterDummy, "Mgg CC (after dummy)", "lep"); // NEW
    leg->AddEntry(&OUT.hMG_Best, "Mgg B_{est} (after dummy)", "l");
    leg->AddEntry(&OUT.hMG_Sub, "Mgg Final SUB", "l");
    leg->Draw();
    std::cerr << "[doDummyThenA] wrote MM_overlay_final\n";

    c.Write("Mgg_overlay_final");
  }

  // ───────────── 2D timing map with region boxes (DATA pairs; CC not subtracted) ─────────────
  {
    // Use DATA pairs for the density map (un-subtracted CC as requested)
    TH2F hTT = D.hTT_pairs; // if you prefer after-dummy, see commented block below

    // Optional: do after-dummy in 2D (uncomment to enable)
    // TH2F hTTdum = M.hTT_pairs_UP; hTTdum.Scale(kUP);
    // { TH2F tmp = M.hTT_pairs_DN; tmp.Scale(kDN); hTTdum.Add(&tmp); }
    // if (Qdum > 0) hTTdum.Scale(Qdata / Qdum);
    // hTT.Add(&hTTdum, -1.0);

    TCanvas cTT("cTT", "Timing map with regions", 1000, 900);
    gPad->SetRightMargin(0.12);
    gPad->SetLogz();
    hTT.SetTitle("Pairs timing map; t_{i} (ns); t_{j} (ns)");
    hTT.Draw("COLZ");

    // Helpers to collect total pos/neg spans for AP “big corners”
    double pos_lo = +1e9, pos_hi = -1e9, neg_lo = +1e9, neg_hi = -1e9;
    for (int i = 0; i < NPOS; ++i)
    {
      pos_lo = std::min(pos_lo, posWins[i].lo);
      pos_hi = std::max(pos_hi, posWins[i].hi);
    }
    for (int i = 0; i < NNEG; ++i)
    {
      neg_lo = std::min(neg_lo, negWins[i].lo);
      neg_hi = std::max(neg_hi, negWins[i].hi);
    }

    // Style helper
    auto makeBox = [](double x1, double y1, double x2, double y2, Color_t col, int lw = 3) -> TBox *
    {
      auto b = new TBox(x1, y1, x2, y2);
      b->SetFillStyle(0);
      b->SetLineColor(col);
      b->SetLineWidth(lw);
      return b;
    };

    // CC box
    auto bCC = makeBox(C_LO, C_LO, C_HI, C_HI, kBlack, 4);
    bCC->Draw("same");

    // Vertical (V): ti in CC, tj in each sideband stripe
    std::vector<TBox *> vBoxes;
    for (int i = 0; i < NPOS; ++i)
    {
      vBoxes.push_back(makeBox(C_LO, posWins[i].lo, C_HI, posWins[i].hi, kBlue));
    }
    for (int i = 0; i < NNEG; ++i)
    {
      vBoxes.push_back(makeBox(C_LO, negWins[i].lo, C_HI, negWins[i].hi, kBlue));
    }
    for (auto *b : vBoxes)
      b->Draw("same");

    // Horizontal (H): tj in CC, ti in each sideband stripe
    std::vector<TBox *> hBoxes;
    for (int i = 0; i < NPOS; ++i)
    {
      hBoxes.push_back(makeBox(posWins[i].lo, C_LO, posWins[i].hi, C_HI, kGreen + 2));
    }
    for (int i = 0; i < NNEG; ++i)
    {
      hBoxes.push_back(makeBox(negWins[i].lo, C_LO, negWins[i].hi, C_HI, kGreen + 2));
    }
    for (auto *b : hBoxes)
      b->Draw("same");

    // AD (“diagonal squares”): same stripe on same side (LL and UR), EXACT squares only
    std::vector<TBox *> adBoxes;
    for (int i = 0; i < NPOS; ++i)
    {
      adBoxes.push_back(makeBox(posWins[i].lo, posWins[i].lo, posWins[i].hi, posWins[i].hi, kMagenta + 1));
    }
    for (int i = 0; i < NNEG; ++i)
    {
      adBoxes.push_back(makeBox(negWins[i].lo, negWins[i].lo, negWins[i].hi, negWins[i].hi, kMagenta + 1));
    }
    for (auto *b : adBoxes)
      b->Draw("same");

    // AP (“pure A”): big off-diagonal corners (UL and LR)
    auto bAP_UL = makeBox(neg_lo, pos_lo, neg_hi, pos_hi, kOrange + 1);
    auto bAP_LR = makeBox(pos_lo, neg_lo, pos_hi, neg_hi, kOrange + 1);
    bAP_UL->Draw("same");
    bAP_LR->Draw("same");

    // Legend
    auto leg = new TLegend(0.13, 0.78, 0.45, 0.93);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(bCC, "CC window", "l");
    leg->AddEntry(vBoxes.front(), "Vertical (V)", "l");
    leg->AddEntry(hBoxes.front(), "Horizontal (H)", "l");
    leg->AddEntry(adBoxes.front(), "Diag A (AD)", "l");
    leg->AddEntry(bAP_UL, "Pure A (AP)", "l");
    leg->Draw();

    cTT.Write("TT_regions_overlay");
  }
  fout.Close();
  std::cout << "Wrote " << outF << "\n";
  return 0;
}
