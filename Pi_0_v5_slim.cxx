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

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TDirectory.h"

// ───────────────────────── constants / config ─────────────────────────
static constexpr double mp = 0.938272081;      // GeV
static constexpr double e0_nom = 10.54350201;  // GeV beam energy (tune per run if needed)

static constexpr double NPS_theta_deg = -13.43;              // deg
static constexpr double NPS_theta_rad = NPS_theta_deg * M_PI / 180.0;
static constexpr double NPS_dist_cm   = 407.0;               // cm

struct Win { double lo, hi; };
static constexpr double C_LO = 149.0;   // ns
static constexpr double C_HI = 151.0;   // ns
static const Win posWins[] = { {153.0,155.0}, {155.0,157.0}, {157.0,159.0} };
static const Win negWins[] = { {145.0,147.0}, {143.0,145.0}, {141.0,143.0} };
static const int NPOS = int(sizeof(posWins)/sizeof(posWins[0]));
static const int NNEG = int(sizeof(negWins)/sizeof(negWins[0]));

// Binning
static constexpr int    nMM  = 300; static constexpr double mmLo=0.0, mmHi=5.0;
static constexpr int    nMG  = 200; static constexpr double mgLo=0.0, mgHi=1.0;

// v5 timing dummy normalization factors
static constexpr double kUP_v5 = 8.467; // upstream
static constexpr double kDN_v5 = 4.256; // downstream

// ─────────────────────────── helpers ───────────────────────────
inline void rotY_passive(double x, double y, double z, double deg,
                         double &xo, double &yo, double &zo)
{
  double th = deg * M_PI / 180.0, c = std::cos(th), s = std::sin(th);
  xo =  c*x - s*z;
  yo =  y;
  zo =  s*x + c*z;
}

inline bool passHMSCuts(double edt, double dp, double et, double npe, double th, double ph)
{
  return (edt < 0.1 && std::fabs(dp) <= 8.5 && et > 0.6 && npe > 1.0 &&
          std::fabs(th) <= 0.09 && std::fabs(ph) <= 0.09);
}

static inline bool inWin(double t, const Win &w){ return (t >= w.lo && t <= w.hi); }
static inline bool inAny(double t, const Win *arr, int n){ for(int i=0;i<n;++i) if (inWin(t,arr[i])) return true; return false; }
static inline bool isCC (double ti,double tj){ return inWin(ti,{C_LO,C_HI}) && inWin(tj,{C_LO,C_HI}); }
static inline bool isV  (double ti,double tj){ return inWin(ti,{C_LO,C_HI}) && (inAny(tj,posWins,NPOS)||inAny(tj,negWins,NNEG)); }
static inline bool isH  (double ti,double tj){ return (inAny(ti,posWins,NPOS)||inAny(ti,negWins,NNEG)) && inWin(tj,{C_LO,C_HI}); }
static inline bool isAA (double ti,double tj){ bool Ai=(inAny(ti,posWins,NPOS)||inAny(ti,negWins,NNEG)); bool Aj=(inAny(tj,posWins,NPOS)||inAny(tj,negWins,NNEG)); return (Ai&&Aj); }

static inline bool computeMM(double E1, double E2,
                             double x1, double y1, double x2, double y2, // cm at NPS plane
                             double Ep, double epx, double epy, double epz,
                             double &mm_out)
{
  // photon directions from (x,y,NPS_dist)
  double px1_h, py1_h, pz1_h, px2_h, py2_h, pz2_h;
  rotY_passive(x1, y1, NPS_dist_cm, NPS_theta_deg, px1_h, py1_h, pz1_h);
  rotY_passive(x2, y2, NPS_dist_cm, NPS_theta_deg, px2_h, py2_h, pz2_h);
  const double r1 = std::sqrt(px1_h*px1_h + py1_h*py1_h + pz1_h*pz1_h);
  const double r2 = std::sqrt(px2_h*px2_h + py2_h*py2_h + pz2_h*pz2_h);
  if (!(r1>0 && r2>0)) return false;
  const double u1x = px1_h/r1, u1y = py1_h/r1, u1z = pz1_h/r1;
  const double u2x = px2_h/r2, u2y = py2_h/r2, u2z = pz2_h/r2;

  // beam + target-at-rest
  const double Ein = e0_nom + mp;
  const double Pinx=0.0, Piny=0.0, Pinz=e0_nom;

  // photons (massless)
  const double p1x = E1*u1x, p1y = E1*u1y, p1z = E1*u1z;
  const double p2x = E2*u2x, p2y = E2*u2y, p2z = E2*u2z;

  const double E_out  = Ep + E1 + E2;
  const double px_out = epx + p1x + p2x;
  const double py_out = epy + p1y + p2y;
  const double pz_out = epz + p1z + p2z;

  const double mm2 = std::pow(Ein - E_out, 2)
                   - std::pow(Pinx - px_out, 2)
                   - std::pow(Piny - py_out, 2)
                   - std::pow(Pinz - pz_out, 2);
  if (!(mm2>0) || !std::isfinite(mm2)) return false;
  mm_out = std::sqrt(mm2);
  return true;
}

static inline bool computeMgg(double E1,double E2,
                              double x1,double y1,double x2,double y2,
                              double &mgg_out)
{
  // small-angle approx from separation at NPS plane
  const double dx = (x1 - x2);
  const double dy = (y1 - y2);
  const double d  = std::sqrt(dx*dx + dy*dy);
  const double theta12 = std::atan2(d, NPS_dist_cm);
  const double cos12 = std::cos(theta12);
  const double m2 = 2.0 * E1 * E2 * (1.0 - cos12);
  if (!(m2>0) || !std::isfinite(m2)) return false;
  mgg_out = std::sqrt(m2);
  return true;
}

// ─────────────────────── histogram pack ───────────────────────
struct Pack {
  // Data categories (raw)
  TH1F hMM_CC, hMM_V, hMM_H, hMM_AA;
  TH1F hMG_CC, hMG_V, hMG_H, hMG_AA; // Mgg categories (raw)

  // Dummy UP/DN splits (only filled for dummy)
  TH1F hMM_CC_UP, hMM_CC_DN, hMM_V_UP, hMM_V_DN, hMM_H_UP, hMM_H_DN, hMM_AA_UP, hMM_AA_DN;
  TH1F hMG_CC_UP, hMG_CC_DN, hMG_V_UP, hMG_V_DN, hMG_H_UP, hMG_H_DN, hMG_AA_UP, hMG_AA_DN;

  // Derived (after dummy then A-method)
  TH1F hMM_Best, hMM_Sub;  // MM B_est and SUB
  TH1F hMG_Best, hMG_Sub;  // Mgg B_est and SUB

  Pack(const char* tag)
  : hMM_CC (TString::Format("hMM_CC_%s",tag),  "MM CC;M_{X} (GeV);Counts", nMM, mmLo, mmHi)
  , hMM_V  (TString::Format("hMM_V_%s",tag),   "MM V;M_{X} (GeV);Counts",   nMM, mmLo, mmHi)
  , hMM_H  (TString::Format("hMM_H_%s",tag),   "MM H;M_{X} (GeV);Counts",   nMM, mmLo, mmHi)
  , hMM_AA (TString::Format("hMM_AA_%s",tag),  "MM AA;M_{X} (GeV);Counts",  nMM, mmLo, mmHi)
  , hMG_CC (TString::Format("hMG_CC_%s",tag),  "M_{#gamma#gamma} CC;M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi)
  , hMG_V  (TString::Format("hMG_V_%s",tag),   "M_{#gamma#gamma} V;M_{#gamma#gamma} (GeV);Counts",   nMG, mgLo, mgHi)
  , hMG_H  (TString::Format("hMG_H_%s",tag),   "M_{#gamma#gamma} H;M_{#gamma#gamma} (GeV);Counts",   nMG, mgLo, mgHi)
  , hMG_AA (TString::Format("hMG_AA_%s",tag),  "M_{#gamma#gamma} AA;M_{#gamma#gamma} (GeV);Counts",  nMG, mgLo, mgHi)

  , hMM_CC_UP (TString::Format("hMM_CC_UP_%s",tag),  "MM CC UP", nMM, mmLo, mmHi)
  , hMM_CC_DN (TString::Format("hMM_CC_DN_%s",tag),  "MM CC DN", nMM, mmLo, mmHi)
  , hMM_V_UP  (TString::Format("hMM_V_UP_%s",tag),   "MM V UP",  nMM, mmLo, mmHi)
  , hMM_V_DN  (TString::Format("hMM_V_DN_%s",tag),   "MM V DN",  nMM, mmLo, mmHi)
  , hMM_H_UP  (TString::Format("hMM_H_UP_%s",tag),   "MM H UP",  nMM, mmLo, mmHi)
  , hMM_H_DN  (TString::Format("hMM_H_DN_%s",tag),   "MM H DN",  nMM, mmLo, mmHi)
  , hMM_AA_UP (TString::Format("hMM_AA_UP_%s",tag),  "MM AA UP", nMM, mmLo, mmHi)
  , hMM_AA_DN (TString::Format("hMM_AA_DN_%s",tag),  "MM AA DN", nMM, mmLo, mmHi)

  , hMG_CC_UP (TString::Format("hMG_CC_UP_%s",tag),  "Mgg CC UP", nMG, mgLo, mgHi)
  , hMG_CC_DN (TString::Format("hMG_CC_DN_%s",tag),  "Mgg CC DN", nMG, mgLo, mgHi)
  , hMG_V_UP  (TString::Format("hMG_V_UP_%s",tag),   "Mgg V UP",  nMG, mgLo, mgHi)
  , hMG_V_DN  (TString::Format("hMG_V_DN_%s",tag),   "Mgg V DN",  nMG, mgLo, mgHi)
  , hMG_H_UP  (TString::Format("hMG_H_UP_%s",tag),   "Mgg H UP",  nMG, mgLo, mgHi)
  , hMG_H_DN  (TString::Format("hMG_H_DN_%s",tag),   "Mgg H DN",  nMG, mgLo, mgHi)
  , hMG_AA_UP (TString::Format("hMG_AA_UP_%s",tag),  "Mgg AA UP", nMG, mgLo, mgHi)
  , hMG_AA_DN (TString::Format("hMG_AA_DN_%s",tag),  "Mgg AA DN", nMG, mgLo, mgHi)

  , hMM_Best(TString::Format("hMM_Best_%s",tag), "MM B_{est};M_{X} (GeV);Counts", nMM, mmLo, mmHi)
  , hMM_Sub (TString::Format("hMM_Sub_%s",tag),  "MM SUB;M_{X} (GeV);Counts",      nMM, mmLo, mmHi)
  , hMG_Best(TString::Format("hMG_Best_%s",tag), "Mgg B_{est};M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi)
  , hMG_Sub (TString::Format("hMG_Sub_%s",tag),  "Mgg SUB;M_{#gamma#gamma} (GeV);Counts",      nMG, mgLo, mgHi)
  {
    auto s2=[&](TH1F& h){ h.Sumw2(); };
    s2(hMM_CC); s2(hMM_V); s2(hMM_H); s2(hMM_AA);
    s2(hMG_CC); s2(hMG_V); s2(hMG_H); s2(hMG_AA);
    s2(hMM_CC_UP); s2(hMM_CC_DN); s2(hMM_V_UP); s2(hMM_V_DN);
    s2(hMM_H_UP);  s2(hMM_H_DN);  s2(hMM_AA_UP); s2(hMM_AA_DN);
    s2(hMG_CC_UP); s2(hMG_CC_DN); s2(hMG_V_UP); s2(hMG_V_DN);
    s2(hMG_H_UP);  s2(hMG_H_DN);  s2(hMG_AA_UP); s2(hMG_AA_DN);
    s2(hMM_Best); s2(hMM_Sub); s2(hMG_Best); s2(hMG_Sub);
  }
};

// ─────────────────────── file processing ───────────────────────
static void fillFromFile(const TString& inF, const char* tag, Pack& O)
{
  TFile f(inF, "READ");
  if (f.IsZombie()) { std::cerr << "Cannot open " << inF << "\n"; return; }
  TTree *tr = dynamic_cast<TTree*>(f.Get("T"));
  if (!tr)        { std::cerr << "No T tree in " << inF << "\n"; return; }

  // HMS branches (v5 names)
  double edt=0, dp=0, et=0, npe=0, th=0, ph=0;
  double hp=0, hpx=0, hpy=0, hpz=0;
  tr->SetBranchStatus("*", 0);
  tr->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1); tr->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edt);
  tr->SetBranchStatus("H.gtr.dp", 1);              tr->SetBranchAddress("H.gtr.dp", &dp);
  tr->SetBranchStatus("H.cal.etotnorm", 1);        tr->SetBranchAddress("H.cal.etotnorm", &et);
  tr->SetBranchStatus("H.cer.npeSum", 1);          tr->SetBranchAddress("H.cer.npeSum", &npe);
  tr->SetBranchStatus("H.gtr.th", 1);              tr->SetBranchAddress("H.gtr.th", &th);
  tr->SetBranchStatus("H.gtr.ph", 1);              tr->SetBranchAddress("H.gtr.ph", &ph);
  tr->SetBranchStatus("H.gtr.p", 1);               tr->SetBranchAddress("H.gtr.p",  &hp);
  tr->SetBranchStatus("H.gtr.px", 1);              tr->SetBranchAddress("H.gtr.px", &hpx);
  tr->SetBranchStatus("H.gtr.py", 1);              tr->SetBranchAddress("H.gtr.py", &hpy);
  tr->SetBranchStatus("H.gtr.pz", 1);              tr->SetBranchAddress("H.gtr.pz", &hpz);

  // NPS clusters (arrays)
  double ncl = 0; static const int MAX=10000; double cE[MAX], cT[MAX], cX[MAX], cY[MAX];
  tr->SetBranchStatus("NPS.cal.nclust", 1); tr->SetBranchAddress("NPS.cal.nclust", &ncl);
  tr->SetBranchStatus("NPS.cal.clusE", 1);  tr->SetBranchAddress("NPS.cal.clusE", cE);
  tr->SetBranchStatus("NPS.cal.clusT", 1);  tr->SetBranchAddress("NPS.cal.clusT", cT);
  tr->SetBranchStatus("NPS.cal.clusX", 1);  tr->SetBranchAddress("NPS.cal.clusX", cX);
  tr->SetBranchStatus("NPS.cal.clusY", 1);  tr->SetBranchAddress("NPS.cal.clusY", cY);

  const bool isDummy = (std::string(tag) == "dummy");

  const Long64_t N = tr->GetEntries();
  for (Long64_t ie=0; ie<N; ++ie){
    tr->GetEntry(ie);
    if (!passHMSCuts(edt,dp,et,npe,th,ph)) continue;
    if (ncl < 2) continue;

    for (int a=0; a<(int)ncl; ++a){
      // loose cluster quality similar to v5 "good cluster" (t window here is not a cut for categories)
      if (cE[a] < 0.6) continue;
      for (int b=0; b<(int)ncl; ++b){ if (b==a) continue; if (cE[b] < 0.6) continue;
        const double ti = cT[a], tj = cT[b];

        double mm=0.0; if (!computeMM(cE[a],cE[b], cX[a],cY[a], cX[b],cY[b], hp,hpx,hpy,hpz, mm)) continue;
        double mg=0.0; (void)computeMgg(cE[a],cE[b], cX[a],cY[a], cX[b],cY[b], mg); // ok if false

        const double yavg = 0.5*(cY[a] + cY[b]); // sign for UP/DN split

        if (isCC(ti,tj)) {
          O.hMM_CC.Fill(mm); if (mg>0) O.hMG_CC.Fill(mg);
          if (isDummy) { (yavg>=0? O.hMM_CC_DN : O.hMM_CC_UP).Fill(mm);
                         if (mg>0) (yavg>=0? O.hMG_CC_DN : O.hMG_CC_UP).Fill(mg); }
        } else if (isV(ti,tj)) {
          O.hMM_V.Fill(mm);  if (mg>0) O.hMG_V.Fill(mg);
          if (isDummy) { (yavg>=0? O.hMM_V_DN : O.hMM_V_UP).Fill(mm);
                         if (mg>0) (yavg>=0? O.hMG_V_DN : O.hMG_V_UP).Fill(mg); }
        } else if (isH(ti,tj)) {
          O.hMM_H.Fill(mm);  if (mg>0) O.hMG_H.Fill(mg);
          if (isDummy) { (yavg>=0? O.hMM_H_DN : O.hMM_H_UP).Fill(mm);
                         if (mg>0) (yavg>=0? O.hMG_H_DN : O.hMG_H_UP).Fill(mg); }
        } else if (isAA(ti,tj)) {
          O.hMM_AA.Fill(mm); if (mg>0) O.hMG_AA.Fill(mg);
          if (isDummy) { (yavg>=0? O.hMM_AA_DN : O.hMM_AA_UP).Fill(mm);
                         if (mg>0) (yavg>=0? O.hMG_AA_DN : O.hMG_AA_UP).Fill(mg); }
        }
      }
    }
  }
}

// ───────────── dummy-first then A-method (for both MM and Mgg) ─────────────
static void doDummyThenA(const Pack& D, const Pack& M, double Qdata, double Qdum,
                         double kUP, double kDN, Pack& OUT)
{
  auto norm_dummy = [&](const TH1F& hUP, const TH1F& hDN)->TH1F{
    TH1F h = hUP; h.Reset("ICESM");
    if (Qdum>0) {
      TH1F hup = hUP; hup.Scale(1.0 / (Qdum * kUP));
      TH1F hdn = hDN; hdn.Scale(1.0 / (Qdum * kDN));
      h.Add(&hup, 1.0); h.Add(&hdn, 1.0); // per µC
      h.Scale(Qdata);                     // to data exposure
    }
    return h;
  };

  // MM categories after dummy
  TH1F hCC_mm = D.hMM_CC; TH1F hCCm = norm_dummy(M.hMM_CC_UP, M.hMM_CC_DN); hCC_mm.Add(&hCCm, -1.0);
  TH1F hV_mm  = D.hMM_V;  TH1F hVm  = norm_dummy(M.hMM_V_UP,  M.hMM_V_DN);  hV_mm .Add(&hVm,  -1.0);
  TH1F hH_mm  = D.hMM_H;  TH1F hHm  = norm_dummy(M.hMM_H_UP,  M.hMM_H_DN);  hH_mm .Add(&hHm,  -1.0);
  TH1F hAA_mm = D.hMM_AA; TH1F hAAm = norm_dummy(M.hMM_AA_UP, M.hMM_AA_DN); hAA_mm.Add(&hAAm, -1.0);

  
  // Mgg categories after dummy
  TH1F hCC_mg = D.hMG_CC; TH1F hCCg = norm_dummy(M.hMG_CC_UP, M.hMG_CC_DN); hCC_mg.Add(&hCCg, -1.0);
  TH1F hV_mg  = D.hMG_V;  TH1F hVg  = norm_dummy(M.hMG_V_UP,  M.hMG_V_DN);  hV_mg .Add(&hVg,  -1.0);
  TH1F hH_mg  = D.hMG_H;  TH1F hHg  = norm_dummy(M.hMG_H_UP,  M.hMG_H_DN);  hH_mg .Add(&hHg,  -1.0);
  TH1F hAA_mg = D.hMG_AA; TH1F hAAg = norm_dummy(M.hMG_AA_UP, M.hMG_AA_DN); hAA_mg.Add(&hAAg, -1.0);

  // Coefficients (geometric Option A)
  double Wacc=0.0; for(int i=0;i<NPOS;++i) Wacc += (posWins[i].hi-posWins[i].lo);
                   for(int i=0;i<NNEG;++i) Wacc += (negWins[i].hi-negWins[i].lo);
  const double Wsig = (C_HI - C_LO);
  const double ACC  = Wsig * Wsig;
  const double aV = (Wacc>0 ? Wsig / Wacc        : 0.0);
  const double aH = (Wacc>0 ? Wsig / Wacc        : 0.0);
  const double aA = (Wacc>0 ? ACC  / (Wacc*Wacc) : 0.0);

  // MM A-method on after-dummy categories
  OUT.hMM_Best.Reset();
  OUT.hMM_Best.Add(&hV_mm,  +aV);
  OUT.hMM_Best.Add(&hH_mm,  +aH);
  OUT.hMM_Best.Add(&hAA_mm, -aA);

  OUT.hMM_Sub.Reset();
  OUT.hMM_Sub.Add(&hCC_mm, 1.0);
  OUT.hMM_Sub.Add(&OUT.hMM_Best, -1.0);

  // Mgg A-method on after-dummy categories
  OUT.hMG_Best.Reset();
  OUT.hMG_Best.Add(&hV_mg,  +aV);
  OUT.hMG_Best.Add(&hH_mg,  +aH);
  OUT.hMG_Best.Add(&hAA_mg, -aA);

  OUT.hMG_Sub.Reset();
  OUT.hMG_Sub.Add(&hCC_mg, 1.0);
  OUT.hMG_Sub.Add(&OUT.hMG_Best, -1.0);

  // Optionally write after-dummy categories for QA
  hCC_mm.SetName("hMM_CC_afterDummy"); hV_mm.SetName("hMM_V_afterDummy"); hH_mm.SetName("hMM_H_afterDummy"); hAA_mm.SetName("hMM_AA_afterDummy");
  hCC_mg.SetName("hMG_CC_afterDummy"); hV_mg.SetName("hMG_V_afterDummy"); hH_mg.SetName("hMG_H_afterDummy"); hAA_mg.SetName("hMG_AA_afterDummy");
  hCC_mm.Write(); hV_mm.Write(); hH_mm.Write(); hAA_mm.Write();
  hCC_mg.Write(); hV_mg.Write(); hH_mg.Write(); hAA_mg.Write();
}

// ─────────────────────────────────── main ───────────────────────────────────
int main(int argc, char** argv)
{
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <data.root> <dummy.root> <out.root> [Qdata] [Qdummy]\n";
    return 1;
  }
  TString dataF  = argv[1];
  TString dummyF = argv[2];
  TString outF   = argv[3];

  double Qdata = 1.0, Qdum = 1.0;
  if (argc >= 5) Qdata = atof(argv[4]);
  if (argc >= 6) Qdum  = atof(argv[5]);

  std::cout << "[Inputs] data=" << dataF << "  dummy=" << dummyF << "  out=" << outF << "\n";
  std::cout << "[Norm]   Q_data=" << Qdata << "  Q_dummy=" << Qdum
            << "  (kUP=" << kUP_v5 << ", kDN=" << kDN_v5 << ")\n";

  // Fill from files
  Pack D("data");  fillFromFile(dataF,  "data",  D);
  Pack M("dummy"); fillFromFile(dummyF, "dummy", M);

  // Write outputs
  TFile fout(outF, "RECREATE");

// Open output FIRST so any internal Write() calls (e.g. in doDummyThenA) work
   fout.cd();

  // Do dummy-first then Option-A for both MM and Mgg
  Pack OUT("final");
  doDummyThenA(D, M, Qdata, Qdum, kUP_v5, kDN_v5, OUT);

  // Raw categories
  fout.cd();
  D.hMM_CC.Write(); D.hMM_V.Write(); D.hMM_H.Write(); D.hMM_AA.Write();
  D.hMG_CC.Write(); D.hMG_V.Write(); D.hMG_H.Write(); D.hMG_AA.Write();
  M.hMM_CC.Write(); M.hMM_V.Write(); M.hMM_H.Write(); M.hMM_AA.Write();
  M.hMG_CC.Write(); M.hMG_V.Write(); M.hMG_H.Write(); M.hMG_AA.Write();

  // Dummy splits
  M.hMM_CC_UP.Write(); M.hMM_CC_DN.Write(); M.hMM_V_UP.Write(); M.hMM_V_DN.Write();
  M.hMM_H_UP.Write();  M.hMM_H_DN.Write();  M.hMM_AA_UP.Write(); M.hMM_AA_DN.Write();
  M.hMG_CC_UP.Write(); M.hMG_CC_DN.Write(); M.hMG_V_UP.Write(); M.hMG_V_DN.Write();
  M.hMG_H_UP.Write();  M.hMG_H_DN.Write();  M.hMG_AA_UP.Write(); M.hMG_AA_DN.Write();

  // Final background & subtracted spectra
  OUT.hMM_Best.Write(); OUT.hMM_Sub.Write();
  OUT.hMG_Best.Write(); OUT.hMG_Sub.Write();

  // Simple overlays
  {
    double ymax = 0.0;
    ymax = std::max(ymax, OUT.hMM_Sub.GetMaximum());
    ymax = std::max(ymax, OUT.hMM_Best.GetMaximum());
    if (ymax <= 0) ymax = 1.0;

    TCanvas c("cMM","MM (dummy-first, Option A)", 1000, 700);
    TH1F axis("axis","Missing Mass;M_{X} (GeV);Counts", nMM, mmLo, mmHi); axis.SetMaximum(1.15*ymax); axis.Draw("AXIS");
    OUT.hMM_Best.SetLineColor(kOrange+7); OUT.hMM_Best.SetLineWidth(2); OUT.hMM_Best.Draw("HIST SAME");
    OUT.hMM_Sub .SetLineColor(kGreen+2);  OUT.hMM_Sub .SetLineWidth(3); OUT.hMM_Sub .Draw("HIST SAME");
    auto leg = new TLegend(0.58,0.68,0.90,0.90); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(&OUT.hMM_Best, "MM B_{est} (after dummy)", "l");
    leg->AddEntry(&OUT.hMM_Sub,  "MM Final SUB", "l");
    leg->Draw();
    c.Write("MM_overlay_final");
  }
  {
    double ymax = 0.0;
    ymax = std::max(ymax, OUT.hMG_Sub.GetMaximum());
    ymax = std::max(ymax, OUT.hMG_Best.GetMaximum());
    if (ymax <= 0) ymax = 1.0;

    TCanvas c("cMG","Mgg (dummy-first, Option A)", 1000, 700);
    TH1F axis("axis","M_{#gamma#gamma};M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi); axis.SetMaximum(1.15*ymax); axis.Draw("AXIS");
    OUT.hMG_Best.SetLineColor(kOrange+7); OUT.hMG_Best.SetLineWidth(2); OUT.hMG_Best.Draw("HIST SAME");
    OUT.hMG_Sub .SetLineColor(kGreen+2);  OUT.hMG_Sub .SetLineWidth(3); OUT.hMG_Sub .Draw("HIST SAME");
    auto leg = new TLegend(0.58,0.68,0.90,0.90); leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(&OUT.hMG_Best, "Mgg B_{est} (after dummy)", "l");
    leg->AddEntry(&OUT.hMG_Sub,  "Mgg Final SUB", "l");
    leg->Draw();
    c.Write("Mgg_overlay_final");
  }

  fout.Close();
  std::cout << "Wrote " << outF << "\n";

  return 0;
}
