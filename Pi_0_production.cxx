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
//   g++ -O0 -g -std=c++17 -o Pi_0_v5_slim_pairsA Pi_0_v5_slim_pairsA.cxx alglib_src/*.cpp -I. -Ialglib_src `root-config --cflags --libs` -lTMVA -lRooFitCore -lRooFit
//   ./Pi0_end2end data.root dummy.root out.root [Qdata] [Qdummy]
//   ./Pi_0_v5_slim /cache/hallc/c-nps/analysis/pass2/replays/production/nps_hms_coin_4205_0_1_-1.root VolatileROOTfiles/dummy_x58_q51_p5_merged.root 4205_v5_slim.root 29446.282 254792.874
//
// Notes:
// - Tree name is assumed "T" (same as v5). Edit if needed.
// - HMS & NPS branches copied from v5; adjust branch names if your files differ.
// - MM kinematics uses beam+proton target and HMS e' 4-vector; NPS photon dirs from (x,y,NPS_dist) rotated by NPS theta.
// - Mgg uses a simple small-angle estimate from cluster separations; replace with your full routine if desired.

#include <regex>
#include <memory>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip> // for std::fixed, std::setprecision
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooPlot.h>
#include <RooBernstein.h>
#include <yaml-cpp/yaml.h> // for YAML::Node / YAML::LoadFile

#ifdef HAVE_YAML_CPP
#include <yaml-cpp/yaml.h>
#endif

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
#include "TParameter.h"
#include "TVectorD.h"
#include "TUUID.h"
#include "TH2D.h"
#include "TLeaf.h"        // for TLeaf::GetName(), GetTypeName()
#include "RooFitResult.h" // for owning/return type of fitTo(...)

// ======= [ADD] structs for normalization & mgg-window =======
struct RunNorm
{
    double Q_data = 0.0;  // from CLI
    double Q_dummy = 0.0; // from CLI
    double kUP = 1.0;     // from config
    double kDN = 1.0;     // from config
} gNorm;

struct MGW
{
    bool ok = false;
    double mu = 0.0;
    double sigma = 0.0;
    double lo = 0.0;
    double hi = 0.0;
    double nsig = 3.0; // default: 2σ window

    // Added to support YAML-configurable fit:
    int order = 2;                     // Bernstein(k) order
    int rebin = 3;                     // histogram rebin factor
    std::string signal_mode = "gauss"; // "gauss" | "doubleG" | "cb"

    // --- Added for S+B RooFit diagnostics and event-level signal weights ---
    double nS_total = 0.0; // fitted total signal yield
    double nB_total = 0.0; // fitted total background yield
    double S_window = 0.0; // signal within ±nsig window
    double B_window = 0.0; // background within ±nsig window
    std::vector<double> Sfrac_bins; // per-bin S/(S+B) fractions
} g_MGW;

// ======= [ADD] simple PDF bin integral helper =======
static double pdfYieldInBin(const RooAbsPdf &pdf, RooRealVar &x,
                            double xlo, double xhi, double Ntot)
{
    const double oldLo = x.getMin();
    const double oldHi = x.getMax();
    x.setRange("bin", xlo, xhi);
    const std::unique_ptr<RooAbsReal> I(pdf.createIntegral(x, RooFit::NormSet(x), RooFit::Range("bin")));
    const double frac = I->getVal();
    x.setRange(oldLo, oldHi);
    return Ntot * frac;
}

// Enable ALL-PAIRS Template-A with pairs→event scaling (matches full v5)
#ifndef USE_PAIR_TEMPLATE
#define USE_PAIR_TEMPLATE 1 // set to 0 to switch the multiplicity on
#endif

// Δm_eff = (m_inv - m_pi0) after taper  vs  corr_fac actually used
static TH2D hDeltaM_vs_CorrFac("hDeltaM_vs_CorrFac",
                               "tapered #Delta m vs corr_{fac};corr_{fac} [GeV];tapered #Delta m [GeV]",
                               300, -60.0, 60.0,  // corr_fac axis
                               200, -0.20, 0.20); // tapered delta_m axis

// Mx_corr vs Mx_raw (both in GeV)
static TH2D hMxCorr_vs_MxRaw("hMxCorr_vs_MxRaw",
                             "M_{x}^{corr} vs M_{x}^{raw};M_{x}^{raw} [GeV];M_{x}^{corr} [GeV]",
                             300, 0.0, 3.0,
                             300, 0.0, 3.0);

// ───────────────────────── constants / config ─────────────────────────
static constexpr double mp = 0.938272081;     // GeV
static constexpr double e0_nom = 10.54350201; // GeV beam energy (tune per run if needed)
static constexpr double m_pi0 = 0.1349768;    // GeV (PDG π0 mass)

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
static constexpr double mgLo = 0.05, mgHi = 0.2;

// v5 timing dummy normalization factors
static double kUP_v5 = 8.467; // upstream
static double kDN_v5 = 4.256; // downstream

// ---- mγγ window (placed early so MR_compute_config_hash() can see it) ----
struct MggWindow
{
    double mu = 0.135, sigma = 0.006, lo = 0.129, hi = 0.141, nsig = 2.0;
    std::string signal_mode = "gauss";
    int order = 2, rebin = 1;
    bool ok = false;
};

// ===================== Map→Reduce support (schema + accumulator) =====================
struct BinSchema
{
    std::vector<std::string> axisNames;     // e.g., {"MxCorr","Mgg"}
    std::vector<std::vector<double>> edges; // edges per axis
    std::string id;                         // stable hash
    bool valid() const { return !axisNames.empty() && axisNames.size() == edges.size(); }
};

// FNV-1a hash for edge arrays
static std::string MR_hashEdges(const std::vector<std::vector<double>> &E)
{
    uint64_t h = 1469598103934665603ull;
    auto mix = [&](const void *p, size_t n)
    {
        const unsigned char *s = (const unsigned char *)p;
        for (size_t i = 0; i < n; i++)
        {
            h ^= s[i];
            h *= 1099511628211ull;
        }
    };
    for (auto &v : E)
    {
        if (!v.empty())
            mix(v.data(), v.size() * sizeof(double));
    }
    char buf[32];
    snprintf(buf, sizeof(buf), "%016llx", (unsigned long long)h);
    return std::string(buf);
}

// Load schema from a ROOT file at BinEdges/{AxisName}
static BinSchema MR_loadBinSchema(const std::string &path)
{
    BinSchema s;
    if (path.empty())
        return s;
    TFile f(path.c_str(), "READ");
    if (!f.IsOpen())
    {
        std::cerr << "[binschema] cannot open " << path << "\n";
        return s;
    }
    TDirectory *d = (TDirectory *)f.Get("BinEdges");
    if (!d)
    {
        std::cerr << "[binschema] no BinEdges/ in " << path << "\n";
        return s;
    }
    f.cd("BinEdges");

    // NOTE: start with MxCorr and Mgg; extend later with Q2,xB,z,pT,phi,tprime
    std::vector<std::string> want = {"MxCorr", "Mgg"};
    for (auto &nm : want)
    {
        if (auto *v = (TVectorD *)gDirectory->Get(nm.c_str()))
        {
            s.axisNames.push_back(nm);
            s.edges.emplace_back(v->GetNoElements());
            for (int i = 0; i < v->GetNoElements(); ++i)
                s.edges.back()[i] = (*v)[i];
        }
    }
    f.cd();
    if (!s.valid())
    {
        std::cerr << "[binschema] invalid/empty " << path << "\n";
        return s;
    }
    s.id = MR_hashEdges(s.edges);
    std::cerr << "[binschema] axes=";
    for (size_t i = 0; i < s.axisNames.size(); ++i)
        std::cerr << s.axisNames[i] << (i + 1 < s.axisNames.size() ? ", " : "");
    std::cerr << "  id=" << s.id << "\n";
    return s;
};

// Our 5 categories (match your code): CC, V, H, AD (=AA_diag), AP (=AA_pure)
struct MR_BinRow
{
    // tallies before algebra (counts)
    double CC_data = 0, V_data = 0, H_data = 0, AD_data = 0, AP_data = 0;
    double CC_up = 0, V_up = 0, H_up = 0, AD_up = 0, AP_up = 0;
    double CC_dn = 0, V_dn = 0, H_dn = 0, AD_dn = 0, AP_dn = 0;
    // final algebra results
    double Sum = 0, Sum2 = 0;
};

struct MR_Accumulator
{
    BinSchema schema;
    std::vector<MR_BinRow> bins;

    void init(const BinSchema &s)
    {
        schema = s;
        int nTot = 1;
        for (auto &e : s.edges)
            nTot *= std::max(0, (int)e.size() - 1);
        bins.assign(std::max(0, nTot), MR_BinRow{});
    }

    // locate bin for up to 2 axes (MxCorr, Mgg); extend if/when you add more
    int locate(double ax0, double ax1) const
    {
        if (!schema.valid())
            return -1;
        int ia = -1, ib = -1;
        for (size_t i = 0; i < schema.axisNames.size(); ++i)
        {
            auto &nm = schema.axisNames[i];
            auto &e = schema.edges[i];
            if (nm == "MxCorr")
            {
                ia = (int)(std::upper_bound(e.begin(), e.end(), ax0) - e.begin()) - 1;
                if (ia < 0 || ia + 1 >= (int)e.size())
                    return -1;
            }
            else if (nm == "Mgg")
            {
                ib = (int)(std::upper_bound(e.begin(), e.end(), ax1) - e.begin()) - 1;
                if (ib < 0 || ib + 1 >= (int)e.size())
                    return -1;
            }
        }
        if (schema.axisNames.size() == 1)
            return ia;
        if (schema.axisNames.size() == 2)
        {
            int nx = (int)schema.edges[0].size() - 1;
            return (ia >= 0 && ib >= 0) ? ib * nx + ia : -1;
        }
        return -1;
    }
};

// globals for Map→Reduce
static std::string MR_run_uid = "";
static int MR_chunk_id = 0;
static std::string MR_schema_path = "";
static bool MR_emit_micro = false; // default OFF
static MR_Accumulator MR_acc;
static BinSchema MR_schema; // global

// DATA event number (bound to input tree branch "evnum")
static ULong64_t g_evnum = 0;
bool g_ev_is_int = false, g_ev_is_dbl = false;
Int_t g_ev_i = 0;
Long64_t g_ev_l = 0;
Double_t g_ev_d = 0.0;

// ---- MR helpers (forward decls) ----
static std::string MR_compute_config_hash();
static inline void MR_extract_axes(double mm_corr, double mgg,
                                   double &ax_mx, double &ax_mg);

// ─────────────────────────── helpers ───────────────────────────
// --- Per-entry record bound to the DATA tree ---
struct GEntry
{
    ULong64_t evnum = 0; // event number (matches ROOT leaf "evnum"; use ULong64_t even if leaf is I)
    double xB = 0, Q2 = 0, z = 0, pT = 0, t = 0, phi = 0;
    double Mx_raw = 0, Mx_corr = 0;
    double mgg = 0; // invariant mass per event/pair, as you compute/fill it
} g;

// === Build ±nsig window from a Gaussian ⊕ Bernstein(k) fit on the final mgg (dummy+acc subtracted) ===
static void BuildMggWindow_Bernstein(TH1 *hFinal_mgg, int bern_order, int rebin_factor)
{
    if (!hFinal_mgg)
    {
        std::cerr << "[auto-fit] null hFinal_mgg\n";
        g_MGW.ok = false;
        return;
    }

    // work on a local clone; keep original binning for later lookups
    TH1 *hFit = (TH1 *)hFinal_mgg->Clone("hFinal_mgg_fit");
    hFit->SetDirectory(nullptr);
    hFit->Sumw2();
    if (rebin_factor > 1)
    {
        TH1 *tmp = hFit->Rebin(rebin_factor, "hFinal_mgg_fit_reb");
        hFit = dynamic_cast<TH1 *>(tmp);
        if (!hFit)
        {
            std::cerr << "[auto-fit] rebin failed\n";
            g_MGW.ok = false;
            return;
        }
        hFit->SetDirectory(nullptr);
        hFit->Sumw2();
    }

    // fit range (use your study’s if you expose them; these are safe defaults)
    double mmin = 0.08, mmax = 0.18;

    RooRealVar m("m", "m_{#gamma#gamma} [GeV]", mmin, mmax);
    m.setBins(std::max(20, hFit->GetNbinsX()));
    RooDataHist dh("dh", "dh", RooArgList(m), hFit);

    // Signal: Gaussian, seeded from any prior g_MGW if present
    RooRealVar mu("mu", "mu", (g_MGW.mu > 0 ? g_MGW.mu : 0.135), 0.120, 0.150);
    RooRealVar sg("sigma", "sigma", (g_MGW.sigma > 0 ? g_MGW.sigma : 0.010), 0.003, 0.020);
    RooGaussian gaus("gaus", "gaus", m, mu, sg);

    // Background: Bernstein(k)
    const int K = std::max(1, bern_order);
    RooArgList coeffs;
    std::vector<std::unique_ptr<RooRealVar>> keep;
    keep.reserve(K + 1);
    for (int i = 0; i <= K; i++)
    {
        auto v = std::make_unique<RooRealVar>(Form("c%d", i), Form("c%d", i), 0.1, 0.0, 1e6);
        coeffs.add(*v);
        keep.emplace_back(std::move(v));
    }
    RooBernstein bkg("bkg", "bernstein", m, coeffs);

    // Extended S+B
    const double N = hFit->Integral();
    RooRealVar nS("nS", "signal", 0.6 * N, 0.0, 10.0 * N + 1.0);
    RooRealVar nB("nB", "bkg", 0.4 * N, 0.0, 10.0 * N + 1.0);
    RooAddPdf model("model", "S+B", RooArgList(gaus, bkg), RooArgList(nS, nB));
    model.fitTo(dh, RooFit::Extended(true), RooFit::Save(true),
                RooFit::PrintLevel(-1), RooFit::Warnings(false), RooFit::SumW2Error(true));

    // commit μ, σ, window from your configured nsig
    g_MGW.mu = mu.getVal();
    g_MGW.sigma = sg.getVal();
    g_MGW.lo = g_MGW.mu - g_MGW.nsig * g_MGW.sigma;
    g_MGW.hi = g_MGW.mu + g_MGW.nsig * g_MGW.sigma;
    g_MGW.ok = true;

    // commit μ, σ, window from your configured nsig
    g_MGW.mu = mu.getVal();
    g_MGW.sigma = sg.getVal();
    g_MGW.lo = g_MGW.mu - g_MGW.nsig * g_MGW.sigma;
    g_MGW.hi = g_MGW.mu + g_MGW.nsig * g_MGW.sigma;
    g_MGW.ok = true;

    // --- NEW: per-bin S and B yields and S/(S+B) fractions ---
    const int nbins = hFinal_mgg->GetNbinsX();

    g_MGW.Sfrac_bins.assign(nbins + 1, 0.0); // indices 1..nbins
    g_MGW.nS_total = nS.getVal();
    g_MGW.nB_total = nB.getVal();
    g_MGW.S_window = 0.0;
    g_MGW.B_window = 0.0;

    for (int ib = 1; ib <= nbins; ++ib)
    {
        const double xlo = hFinal_mgg->GetXaxis()->GetBinLowEdge(ib);
        const double xhi = hFinal_mgg->GetXaxis()->GetBinUpEdge(ib);

        const double Sj = pdfYieldInBin(gaus, m, xlo, xhi, g_MGW.nS_total);
        const double Bj = pdfYieldInBin(bkg,  m, xlo, xhi, g_MGW.nB_total);
        const double denom = Sj + Bj;

        const double fracS = (denom > 0.0) ? (Sj / denom) : 0.0;
        g_MGW.Sfrac_bins[ib] = fracS;

        // optional, if you want total S/B in the ±nsig window
        if (xhi >= g_MGW.lo && xlo <= g_MGW.hi)
        {
            g_MGW.S_window += Sj;
            g_MGW.B_window += Bj;
        }
    }


    // provenance (optional)
    TTree *tCut = new TTree("MggCut", "MggCut");
    double mu_out = g_MGW.mu, sg_out = g_MGW.sigma, lo_out = g_MGW.lo, hi_out = g_MGW.hi, ns = g_MGW.nsig;
    int ord_out = K, reb_out = rebin_factor;
    tCut->Branch("mu", &mu_out, "mu/D");
    tCut->Branch("sigma", &sg_out, "sigma/D");
    tCut->Branch("lo", &lo_out, "lo/D");
    tCut->Branch("hi", &hi_out, "hi/D");
    tCut->Branch("nsig", &ns, "nsig/D");
    tCut->Branch("order", &ord_out, "order/I");
    tCut->Branch("rebin", &reb_out, "rebin/I");
    tCut->Fill();
    tCut->Write();

    std::cout << "[auto-fit] mu=" << g_MGW.mu << " sigma=" << g_MGW.sigma
              << " nsig=" << g_MGW.nsig << " window=[" << g_MGW.lo << "," << g_MGW.hi << "]\n";
}

// ---- MR helpers (definitions) ----
static std::string MR_compute_config_hash()
{
    std::ostringstream ss;
    ss << std::setprecision(8)
       << "CWIN:" << C_LO << "," << C_HI
       << ";SIDES:" << NPOS << "," << NNEG
       << ";UD:" << kUP_v5 << "," << kDN_v5
       // NOTE: use g_MGW.mu etc. (no parentheses)
       << ";g_MGW:" << g_MGW.mu << "," << g_MGW.sigma << "," << g_MGW.nsig
       << ";MM:" << nMM << "," << mmLo << "," << mmHi
       << ";MG:" << nMG << "," << mgLo << "," << mgHi;

    auto s = ss.str();
    uint64_t h = 1469598103934665603ull; // FNV-1a
    for (unsigned char c : s)
    {
        h ^= c;
        h *= 1099511628211ull;
    }
    char buf[32];
    snprintf(buf, sizeof(buf), "%016llx", (unsigned long long)h);
    return std::string(buf);
}

static inline void MR_extract_axes(double mm_corr, double mgg,
                                   double &ax_mx, double &ax_mg)
{
    ax_mx = mm_corr;
    ax_mg = mgg;
}

static double parse_json_double(const std::string &s, const char *key, double defval)
{
    std::string pat = std::string("\"") + key + "\"";
    size_t p = s.find(pat);
    if (p == std::string::npos)
        return defval;
    p = s.find(':', p);
    if (p == std::string::npos)
        return defval;
    while (p < s.size() && (s[p] == ':' || s[p] == ' '))
        ++p;
    char *endp = nullptr;
    double v = std::strtod(s.c_str() + p, &endp);
    return (endp == s.c_str() + p) ? defval : v;
}

static std::string parse_json_string(const std::string &s, const char *key, const char *defval)
{
    std::string pat = std::string("\"") + key + "\"";
    size_t p = s.find(pat);
    if (p == std::string::npos)
        return defval;
    p = s.find(':', p);
    if (p == std::string::npos)
        return defval;
    p = s.find('"', p);
    if (p == std::string::npos)
        return defval;
    size_t q = s.find('"', p + 1);
    if (q == std::string::npos)
        return defval;
    return s.substr(p + 1, q - (p + 1));
}

static MggWindow load_mgg_window_from_json(const std::string &json_path, double nsig = 3.0)
{
    MggWindow w;
    w.nsig = nsig;
    std::ifstream in(json_path);
    if (!in)
        return w;
    std::stringstream buf;
    buf << in.rdbuf();
    const std::string s = buf.str();
    std::string sm = parse_json_string(s, "signal_mode", "gauss");
    double mu = parse_json_double(s, "mu", 0.135);
    double sigma = parse_json_double(s, sm == "doubleG" ? "sigma_eff" : (sm == "cb" ? "sigma_cb" : "sigma"), 0.006);
    int order = (int)parse_json_double(s, "order", 2);
    int rebin = (int)parse_json_double(s, "rebin", 1);
    w.signal_mode = sm;
    w.mu = mu;
    w.sigma = sigma;
    w.order = order;
    w.rebin = rebin;
    w.lo = mu - nsig * sigma;
    w.hi = mu + nsig * sigma;
    w.ok = true;
    return w;
}

struct MggWeights
{
    TH1 *hWsig = nullptr; // per-bin w_sig
    TH1 *hUsed = nullptr; // rebinned data used in fit (binning reference)
    bool ok() const { return hWsig && hUsed; }
};

static MggWeights load_mgg_weights(const std::string &path = "fit_mgg_roofit_out.root")
{
    MggWeights W;
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie())
        return W;
    W.hWsig = dynamic_cast<TH1 *>(f->Get("hW_sig"));
    W.hUsed = dynamic_cast<TH1 *>(f->Get("h_mgg_used"));
    if (W.hWsig)
        W.hWsig->SetDirectory(nullptr);
    if (W.hUsed)
        W.hUsed->SetDirectory(nullptr);
    return W;
}

// Globals so fillFromFile() can see them

static MggWeights g_MGWts;

// Optional QA + skim (created in main, filled in fillFromFile)
static TH1 *g_hMG_pass = nullptr;
static TTree *g_tSkim = nullptr;

// Skim branches (keep trivial; you can extend with 4-vectors later)
static ULong64_t g_skim_eventnum = 0; // event number from g.evnum / evnum / eventnum
static double g_skim_mgg = 0.0;
static double g_skim_wsig = 1.0;
static int g_skim_pass_pi0 = 0;

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

// Avnish method with taper/cap to prevent high-Mx smearing.
// Returns corrected missing mass  mm_corr  [GeV]; falls back to sqrt(mm2) on guard failure.
//
// Conventions: matches your computeMM (nu = Ebeam - Ep_used).
static inline double MxCorr_Avnish_pair(double mm2,
                                        double Ebeam,
                                        double Ep_used, // pass hp here
                                        // electron 3-momentum components [GeV]
                                        double px_e, double py_e, double pz_e,
                                        // photon energies [GeV] and NPS-plane coords [cm]
                                        double E1, double E2,
                                        double x1, double y1, double x2, double y2)
{
    // ---- Tunables to control tails ----
    // Smoothly taper the effective (m_inv - m_pi0) using tanh with a 30 MeV scale.
    // Values much larger than DELTA_M_TAPER MeV are strongly down-weighted.
    static constexpr double DELTA_M_TAPER = 0.030; // GeV
    // Hard cap on corr_fac magnitude to prevent explosive rotations in edge cases.
    // corr_fac has units of GeV (so corr_fac * delta_m has GeV^2).
    static constexpr double CORR_FAC_CAP = 40.0; // GeV (adjust 25–50 as needed)
    // Keep your working choices:
    static constexpr bool USE_Q_EQ_KPRIME_MINUS_K = false;
    static constexpr int ROT_TERM_SIGN = -1; // you found +1 looks better

    // --- Photon unit vectors (same geometry as computeMM) ---
    double px1_h, py1_h, pz1_h, px2_h, py2_h, pz2_h;
    rotY_passive(x1, y1, NPS_dist_cm, NPS_theta_deg, px1_h, py1_h, pz1_h);
    rotY_passive(x2, y2, NPS_dist_cm, NPS_theta_deg, px2_h, py2_h, pz2_h);

    const double r1 = std::sqrt(px1_h * px1_h + py1_h * py1_h + pz1_h * pz1_h);
    const double r2 = std::sqrt(px2_h * px2_h + py2_h * py2_h + pz2_h * pz2_h);
    if (!(r1 > 0.0 && r2 > 0.0))
        return (mm2 > 0.0 ? std::sqrt(mm2) : 0.0);

    const double u1x = px1_h / r1, u1y = py1_h / r1, u1z = pz1_h / r1;
    const double u2x = px2_h / r2, u2y = py2_h / r2, u2z = pz2_h / r2;

    // π0 3-momentum and invariant mass (massless photons)
    const double p1x = E1 * u1x, p1y = E1 * u1y, p1z = E1 * u1z;
    const double p2x = E2 * u2x, p2y = E2 * u2y, p2z = E2 * u2z;
    const double pi_px = p1x + p2x;
    const double pi_py = p1y + p2y;
    const double pi_pz = p1z + p2z;
    const double pi_p = std::sqrt(pi_px * pi_px + pi_py * pi_py + pi_pz * pi_pz);

    const double e12 = E1 + E2;
    const double m_inv2 = e12 * e12 - pi_p * pi_p;
    if (!(m_inv2 > 0.0 && std::isfinite(m_inv2)))
        return (mm2 > 0.0 ? std::sqrt(mm2) : 0.0);
    const double m_inv = std::sqrt(m_inv2);

    // --- Virtual photon 3-vector ---
    const double qx = (USE_Q_EQ_KPRIME_MINUS_K ? px_e : -px_e);
    const double qy = (USE_Q_EQ_KPRIME_MINUS_K ? py_e : -py_e);
    const double qz = (USE_Q_EQ_KPRIME_MINUS_K ? (pz_e - Ebeam) : (Ebeam - pz_e));
    const double q_mag = std::sqrt(qx * qx + qy * qy + qz * qz);
    if (!(q_mag > 0.0 && std::isfinite(q_mag)))
        return (mm2 > 0.0 ? std::sqrt(mm2) : 0.0);

    // cos(theta_{π,γ*})
    const double upx = (pi_p > 0 ? pi_px / pi_p : 0.0);
    const double upy = (pi_p > 0 ? pi_py / pi_p : 0.0);
    const double upz = (pi_p > 0 ? pi_pz / pi_p : 0.0);
    double costh = (qx * upx + qy * upy + qz * upz) / q_mag;
    if (costh > 1.0)
        costh = 1.0;
    if (costh < -1.0)
        costh = -1.0;

    // Use the same nu as your Mx builder (Ep_used = hp)
    const double nu = Ebeam - Ep_used;

    // Raw correction factor
    const double corr_fac_raw =
        (2.0 / m_inv) * (m_inv2 - e12 * (nu + mp) + q_mag * costh * pi_p);

    // Cap the correction factor
    double corr_fac = corr_fac_raw;
    if (corr_fac > CORR_FAC_CAP)
        corr_fac = CORR_FAC_CAP;
    if (corr_fac < -CORR_FAC_CAP)
        corr_fac = -CORR_FAC_CAP;

    // Smoothly taper delta_m so outliers don't over-rotate:
    // delta_eff = Δ * tanh( (m_inv - m_pi0) / Δ )
    const double delta_m = (m_inv - m_pi0);
    const double delta_eff = DELTA_M_TAPER * std::tanh(delta_m / DELTA_M_TAPER);

    // Apply rotation with taper & cap
    double mm2_new = mm2 + ROT_TERM_SIGN * corr_fac * delta_eff;

    /// ---- Diagnostics fill (safe, side-effect-free) ----
    {
        // 1) Δm_eff vs corr_fac actually used
        //    Reuse the values you ALREADY computed for the correction path:
        //    - m_inv (computed earlier in this function)
        //    - delta_eff (your tapered Δm actually applied)
        //    - corr_fac  (already capped just above, before the rotation)
        const double m_pi0 = 0.1349768; // GeV
        const double delta_m_diag = (m_inv - m_pi0);

        // Use the EXACT tapered value you applied in the rotation:
        const double delta_m_eff = delta_eff; // (already = DELTA_M_TAPER * tanh(...))

        // Mirror rotation ordering: cap first (already done), then apply sign:
        double corr_fac_eff = corr_fac; // capped corr_fac you just used
        corr_fac_eff *= ROT_TERM_SIGN;  // apply the same sign you used

        hDeltaM_vs_CorrFac.Fill(corr_fac_eff, delta_m_eff);

        // 2) Mx_corr vs Mx_raw
        const double Mx_raw = (mm2 > 0.0 ? std::sqrt(mm2) : 0.0);
        const double Mx_corr = (mm2_new > 0.0 ? std::sqrt(mm2_new) : 0.0);
        hMxCorr_vs_MxRaw.Fill(Mx_raw, Mx_corr);
    }
    // ---- end diagnostics ----

    if (!std::isfinite(mm2_new) || mm2_new <= 0.0)
        return 0.0;
    return std::sqrt(mm2_new);
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

    // 2D decorrelation QA (pairs): Mx^2 vs mgg, before/after rotation
    TH2F h2Mx2_vs_mgg_uncorr; // global (all timing categories)
    TH2F h2Mx2_vs_mgg_corr;
    TH2F h2Mx2_vs_mgg_CC_uncorr; // CC-only
    TH2F h2Mx2_vs_mgg_CC_corr;

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

          ,
          // X: Mx^2 [GeV^2], Y: mgg [GeV]
          h2Mx2_vs_mgg_uncorr(
              Form("h2Mx2_vs_mgg_uncorr_%s", tag),
              "M_{x}^{2} vs m_{#gamma#gamma} (uncorr);M_{x}^{2} [GeV^{2}];m_{#gamma#gamma} [GeV]",
              120, 0.50, 2.50, // X = Mx^2
              100, 0.1, 0.16   // Y = mgg (your existing mg window)
              ),
          h2Mx2_vs_mgg_corr(
              Form("h2Mx2_vs_mgg_corr_%s", tag),
              "M_{x}^{2} vs m_{#gamma#gamma} (corr);M_{x}^{2} [GeV^{2}];m_{#gamma#gamma} [GeV]",
              120, 0.50, 2.50,
              100, 0.1, 0.16),
          h2Mx2_vs_mgg_CC_uncorr(
              Form("h2Mx2_vs_mgg_CC_uncorr_%s", tag),
              "CC: M_{x}^{2} vs m_{#gamma#gamma} (uncorr);M_{x}^{2} [GeV^{2}];m_{#gamma#gamma} [GeV]",
              120, 0.50, 2.50,
              100, 0.1, 0.16),
          h2Mx2_vs_mgg_CC_corr(
              Form("h2Mx2_vs_mgg_CC_corr_%s", tag),
              "CC: M_{x}^{2} vs m_{#gamma#gamma} (corr);M_{x}^{2} [GeV^{2}];m_{#gamma#gamma} [GeV]",
              120, 0.50, 2.50,
              100, 0.1, 0.16)

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

        auto s2_2d = [&](TH2F &h)
        { h.Sumw2(); };
        s2_2d(hTT_pairs);
        s2_2d(hTT_pairs_UP);
        s2_2d(hTT_pairs_DN);
        s2_2d(h2Mx2_vs_mgg_uncorr);
        s2_2d(h2Mx2_vs_mgg_corr);
        s2_2d(h2Mx2_vs_mgg_CC_uncorr);
        s2_2d(h2Mx2_vs_mgg_CC_corr);
    }
};

// ─────────────────────── file processing ───────────────────────
static void fillFromFile(const TString &inF, const char *tag, Pack &O)
{
    const bool isDummy = (tag && std::strcmp(tag, "dummy") == 0);
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

    // ---- Bind event number (handles g.evnum / evnum / eventnum) ----
    // NOTE: detect the leaf's type *before* calling SetBranchAddress to avoid ROOT type errors.
    static bool ev_is_int = false;
    static bool ev_is_dbl = false;
    static Int_t ev_i = 0;      // Int_t / UInt_t
    static Long64_t ev_l = 0;   // Long64_t / ULong64_t
    static Double_t ev_d = 0.0; // Double_t

    if (!isDummy)
    {
        TLeaf *lf = tr->GetLeaf("g.evnum");
        if (!lf)
            lf = tr->GetLeaf("evnum");
        if (!lf)
            lf = tr->GetLeaf("eventnum");

        if (!lf)
        {
            std::cerr << "[bind] ERROR: no evnum leaf (tried g.evnum/evnum/eventnum)\n";
        }
        else
        {
            const char *lname = lf->GetName();
            const char *tname = lf->GetTypeName(); // "Int_t","UInt_t","Long64_t","ULong64_t","Double_t", ...
            tr->SetBranchStatus(lname, 1);

            if (!strcmp(tname, "Int_t") || !strcmp(tname, "UInt_t"))
            {
                ev_is_int = true;
                ev_is_dbl = false;
                tr->SetBranchAddress(lname, &ev_i);
            }
            else if (!strcmp(tname, "Long64_t") || !strcmp(tname, "ULong64_t"))
            {
                ev_is_int = false;
                ev_is_dbl = false;
                tr->SetBranchAddress(lname, &ev_l);
            }
            else if (!strcmp(tname, "Double_t"))
            {
                ev_is_int = false;
                ev_is_dbl = true;
                tr->SetBranchAddress(lname, &ev_d);
            }
            else
            {
                std::cerr << "[bind] ERROR: unsupported evnum type: " << tname << "\n";
            }
            std::cout << "[bind] event number bound: " << lname << " (" << tname << ")\n";
        }
    }
    else
    {
        // Dummy file: skip binding evnum (branch often absent)
        ev_is_int = ev_is_dbl = false;
        std::cerr << "[bind] dummy: skipping evnum bind\n";
    }

    // ===== Map→Reduce: per-category tally =====
    // cat: 0=CC, 1=V, 2=H, 3=AD (AA_diag), 4=AP (AA_pure)
    auto MR_tally_cat = [&](int bin_id, bool isData, int updn, int cat)
    {
        if (!MR_schema.valid())
            return;
        if (bin_id < 0 || bin_id >= (int)MR_acc.bins.size())
            return;
        MR_BinRow &r = MR_acc.bins[bin_id];
        auto bump = [&](double &x)
        { x += 1.0; };
        if (isData)
        {
            if (cat == 0)
                bump(r.CC_data);
            else if (cat == 1)
                bump(r.V_data);
            else if (cat == 2)
                bump(r.H_data);
            else if (cat == 3)
                bump(r.AD_data);
            else
                bump(r.AP_data);
        }
        else
        {
            if (updn > 0)
            { // +1 = DN; -1 = UP (matches your yavg sign)
                if (cat == 0)
                    bump(r.CC_dn);
                else if (cat == 1)
                    bump(r.V_dn);
                else if (cat == 2)
                    bump(r.H_dn);
                else if (cat == 3)
                    bump(r.AD_dn);
                else
                    bump(r.AP_dn);
            }
            else
            {
                if (cat == 0)
                    bump(r.CC_up);
                else if (cat == 1)
                    bump(r.V_up);
                else if (cat == 2)
                    bump(r.H_up);
                else if (cat == 3)
                    bump(r.AD_up);
                else
                    bump(r.AP_up);
            }
        }
    };

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
                                                                                // Avnish correction → corrected Mx for this pair
                double mm_corr = mm;
                if (mg > 0.0)
                {
                    // Ee ≈ hp (ultrarel.), electron components hpx/hpy/hpz are set from branches
                    mm_corr = MxCorr_Avnish_pair(mm * mm, e0_nom, hp, hpx, hpy, hpz,
                                                 cE[a], cE[b], cX[a], cY[a], cX[b], cY[b]);
                }

                // [ADD] global before/after decorrelation views
                if (mg > 0.0 && mm > 0.0)
                    O.h2Mx2_vs_mgg_uncorr.Fill(mm * mm, mg);
                if (mg > 0.0 && mm_corr > 0.0)
                    O.h2Mx2_vs_mgg_corr.Fill(mm_corr * mm_corr, mg);

                const double yavg = 0.5 * (cY[a] + cY[b]); // sign for UP/DN split

                O.hTT_pairs.Fill(ti, tj);
                if (isDummy)
                {
                    (yavg >= 0 ? O.hTT_pairs_DN : O.hTT_pairs_UP).Fill(ti, tj);
                }

                if (isCC(ti, tj))
                {
                    O.hMM_CC.Fill(mm_corr);
                    // [ADD] CC-only before/after decorrelation views
                    if (mg > 0.0 && mm > 0.0)
                        O.h2Mx2_vs_mgg_CC_uncorr.Fill(mm * mm, mg);
                    if (mg > 0.0 && mm_corr > 0.0)
                        O.h2Mx2_vs_mgg_CC_corr.Fill(mm_corr * mm_corr, mg);

                    if (mg > 0)
                        O.hMG_CC.Fill(mg);
                    // ---- Map→Reduce tally (CC) ----
                    if (MR_schema.valid())
                    {
                        double ax_mx = 0, ax_mg = 0;
                        MR_extract_axes(mm_corr, mg, ax_mx, ax_mg);
                        int bid = MR_acc.locate(ax_mx, ax_mg);
                        int updn = (isDummy ? (yavg >= 0 ? +1 : -1) : 0);
                        MR_tally_cat(bid, /*isData=*/!isDummy, /*updn*/ updn, /*cat=*/0);
                    }
                    // --------------------------------

                    if (isDummy)
                    {
                        (yavg >= 0 ? O.hMM_CC_DN : O.hMM_CC_UP).Fill(mm_corr);
                        if (mg > 0)
                            (yavg >= 0 ? O.hMG_CC_DN : O.hMG_CC_UP).Fill(mg);
                    }

                    // update best-per-event CC candidate
                    const double score = std::fabs(ti - t0) + std::fabs(tj - t0);
                    if (score < bestScore)
                    {
                        haveBestCC = true;
                        bestScore = score;
                        bestMM = mm_corr;
                        bestMG = mg;
                        bestYavg = yavg;
                    }
                }
                else if (isV(ti, tj))
                {
                    O.hMM_V.Fill(mm_corr);
                    if (mg > 0)
                        O.hMG_V.Fill(mg);
                    // ---- Map→Reduce tally (V) ----
                    if (MR_schema.valid())
                    {
                        double ax_mx = 0, ax_mg = 0;
                        MR_extract_axes(mm_corr, mg, ax_mx, ax_mg);
                        int bid = MR_acc.locate(ax_mx, ax_mg);
                        int updn = (isDummy ? (yavg >= 0 ? +1 : -1) : 0);
                        MR_tally_cat(bid, /*isData=*/!isDummy, /*updn*/ updn, /*cat=*/1);
                    }
                    // --------------------------------

                    if (isDummy)
                    {
                        (yavg >= 0 ? O.hMM_V_DN : O.hMM_V_UP).Fill(mm_corr);
                        if (mg > 0)
                            (yavg >= 0 ? O.hMG_V_DN : O.hMG_V_UP).Fill(mg);
                    }
                }
                else if (isH(ti, tj))
                {
                    O.hMM_H.Fill(mm_corr);
                    if (mg > 0)
                        O.hMG_H.Fill(mg);
                    // ---- Map→Reduce tally (H) ----
                    if (MR_schema.valid())
                    {
                        double ax_mx = 0, ax_mg = 0;
                        MR_extract_axes(mm_corr, mg, ax_mx, ax_mg);
                        int bid = MR_acc.locate(ax_mx, ax_mg);
                        int updn = (isDummy ? (yavg >= 0 ? +1 : -1) : 0);
                        MR_tally_cat(bid, /*isData=*/!isDummy, /*updn*/ updn, /*cat=*/2);
                    }
                    // --------------------------------

                    if (isDummy)
                    {
                        (yavg >= 0 ? O.hMM_H_DN : O.hMM_H_UP).Fill(mm_corr);
                        if (mg > 0)
                            (yavg >= 0 ? O.hMG_H_DN : O.hMG_H_UP).Fill(mg);
                    }
                }
                else if (isAA_diag(ti, tj))
                {
                    O.hMM_AD.Fill(mm_corr);
                    if (mg > 0)
                        O.hMG_AD.Fill(mg);
                    // ---- Map→Reduce tally (AD) ----
                    if (MR_schema.valid())
                    {
                        double ax_mx = 0, ax_mg = 0;
                        MR_extract_axes(mm_corr, mg, ax_mx, ax_mg);
                        int bid = MR_acc.locate(ax_mx, ax_mg);
                        int updn = (isDummy ? (yavg >= 0 ? +1 : -1) : 0);
                        MR_tally_cat(bid, /*isData=*/!isDummy, /*updn*/ updn, /*cat=*/3);
                    }
                    // --------------------------------

                    if (isDummy)
                    {
                        (yavg >= 0 ? O.hMM_AD_DN : O.hMM_AD_UP).Fill(mm_corr);
                        if (mg > 0)
                            (yavg >= 0 ? O.hMG_AD_DN : O.hMG_AD_UP).Fill(mg);
                    }
                }
                else if (isAA_pure(ti, tj))
                {
                    O.hMM_AP.Fill(mm_corr);
                    if (mg > 0)
                        O.hMG_AP.Fill(mg);
                    // ---- Map→Reduce tally (AP) ----
                    if (MR_schema.valid())
                    {
                        double ax_mx = 0, ax_mg = 0;
                        MR_extract_axes(mm_corr, mg, ax_mx, ax_mg);
                        int bid = MR_acc.locate(ax_mx, ax_mg);
                        int updn = (isDummy ? (yavg >= 0 ? +1 : -1) : 0);
                        MR_tally_cat(bid, /*isData=*/!isDummy, /*updn*/ updn, /*cat=*/4);
                    }
                    // --------------------------------

                    if (isDummy)
                    {
                        (yavg >= 0 ? O.hMM_AP_DN : O.hMM_AP_UP).Fill(mm_corr);
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

            // --- pi0 2σ selection + per-bin background weights (global sidecars) ---
            if (bestMG > 0 && g_MGW.ok)
            {
                const bool pass_mgg_window = (bestMG >= g_MGW.lo && bestMG <= g_MGW.hi);
                if (pass_mgg_window)
                {
                    double w_sig = 1.0;
                    if (g_MGWts.ok())
                    {
                        int ibin = g_MGWts.hUsed->GetXaxis()->FindBin(bestMG);
                        if (ibin < 1)
                            ibin = 1;
                        if (ibin > g_MGWts.hUsed->GetNbinsX())
                            ibin = g_MGWts.hUsed->GetNbinsX();
                        w_sig = g_MGWts.hWsig->GetBinContent(ibin);
                        // If you prefer to forbid negatives, uncomment:
                        // if (w_sig < 0.0) w_sig = 0.0;
                        // if (w_sig > 1.0) w_sig = 1.0;
                    }
                    if (g_hMG_pass)
                        g_hMG_pass->Fill(bestMG, w_sig);
                    if (g_tSkim)
                    {
                        ULong64_t evnum_for_skim =
                            ev_is_dbl ? (ULong64_t)llround(ev_d) : (ev_is_int ? (ULong64_t)ev_i : (ULong64_t)ev_l);

                        g_skim_eventnum = evnum_for_skim;
                        g_skim_mgg = bestMG;
                        g_skim_wsig = w_sig; // store actual signal weight
                        g_skim_pass_pi0 = 1;
                        g_tSkim->Fill();
                    }
                }
            }
        }

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

    // --- Split AA coefficients (AD, AP) derived separately from geometry ---
    const double aAD = (A_AD > 0) ? (ACC / A_AD) : 0.0; // e.g. 1/6 for 2 ns geometry
    const double aAP = (A_AP > 0) ? (ACC / A_AP) : 0.0; // e.g. 1/18 for 2 ns geometry
    const double aV_half = 0.5 * aV;
    const double aH_half = 0.5 * aH;

    std::cout << "[geom α] "
              << "aV=" << aV << " aH=" << aH
              << " aAD=" << aAD << " aAP=" << aAP
              << "  (aV/2=" << aV_half << " aH/2=" << aH_half << ")"
              << std::endl;

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
    hMM_Best_pairs.Add(&hV_mm, +0.5 * aV);
    hMM_Best_pairs.Add(&hH_mm, +0.5 * aH);
    hMM_Best_pairs.Add(&hAD_mm, +aAD);
    hMM_Best_pairs.Add(&hAP_mm, -aAP);

    TH1F hMG_Best_pairs("hMG_Best_pairs", "Mgg B_{est} (ALL pairs);M_{#gamma#gamma} (GeV);Counts", nMG, mgLo, mgHi);
    hMG_Best_pairs.Reset();
    hMG_Best_pairs.Add(&hV_mg, +0.5 * aV);
    hMG_Best_pairs.Add(&hH_mg, +0.5 * aH);
    hMG_Best_pairs.Add(&hAD_mg, +aAD);
    hMG_Best_pairs.Add(&hAP_mg, -aAP);

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
    OUT.hMM_Best.Add(&hV_mm, +0.5 * aV);
    OUT.hMM_Best.Add(&hH_mm, +0.5 * aH);
    OUT.hMM_Best.Add(&hAD_mm, +aAD);
    OUT.hMM_Best.Add(&hAP_mm, -aAP);

    OUT.hMM_Sub.Reset();
    OUT.hMM_Sub.Add(&hCC_mm, 1.0);
    OUT.hMM_Sub.Add(&OUT.hMM_Best, -1.0);

    // Mgg A-method on after-dummy categories
    OUT.hMG_Best.Reset();
    OUT.hMG_Best.Add(&hV_mg, +0.5 * aV);
    OUT.hMG_Best.Add(&hH_mg, +0.5 * aH);
    OUT.hMG_Best.Add(&hAD_mg, +aAD);
    OUT.hMG_Best.Add(&hAP_mg, -aAP);

    OUT.hMG_Sub.Reset();
    OUT.hMG_Sub.Add(&hCC_mg, 1.0);
    OUT.hMG_Sub.Add(&OUT.hMG_Best, -1.0);

#endif

 // >>> build ±nsig window from the final (data−acc−dummy) mγγ <<<
    BuildMggWindow_Bernstein(&OUT.hMG_Sub, g_MGW.order, g_MGW.rebin);

    // Optional: quick print so you see what window was derived
    if (g_MGW.ok)
    {
        std::cout << "[mgg-window] mu=" << g_MGW.mu
                  << "  sigma=" << g_MGW.sigma
                  << "  nsig=" << g_MGW.nsig
                  << "  window=[" << g_MGW.lo << "," << g_MGW.hi << "]\n";
    }
    else
    {
        std::cerr << "[mgg-window] window not available (fit failed)\n";
    }

    OUT.hMM_CC_afterDummy.Write();
    OUT.hMG_CC_afterDummy.Write();

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

    // --- Open output FIRST so anything that writes has a real file target ---
    TFile fout(outF, "RECREATE");
    fout.cd();

    double Qdata = 1.0, Qdum = 1.0;
    if (argc >= 5)
        Qdata = atof(argv[4]);
    if (argc >= 6)
        Qdum = atof(argv[5]);

    // ===== Map→Reduce CLI (non-breaking) =====
    // Accept extra flags anywhere after the positional args
    for (int i = 1; i < argc; i++)
    {

        // --- config file (YAML) ---
        if (!strcmp(argv[i], "--cfg") && i + 1 < argc)
        {
            const std::string cfgPath = argv[++i];
            try
            {
                YAML::Node cfg = YAML::LoadFile(cfgPath);

                // mgg-fit settings (defaults already in g_MGW)
                if (cfg["mgg_fit"])
                {
                    auto m = cfg["mgg_fit"];
                    if (m["mode"])
                        g_MGW.signal_mode = m["mode"].as<std::string>(); // "gauss"|"doubleG"|"cb"
                    if (m["order"])
                        g_MGW.order = m["order"].as<int>(); // Bernstein(k)
                    if (m["rebin"])
                        g_MGW.rebin = m["rebin"].as<int>();
                    if (m["nsig"])
                        g_MGW.nsig = m["nsig"].as<double>(); // 2 or 3
                }

                // timing dummy normalization
                if (cfg["timing"] && cfg["timing"]["kUP"] && cfg["timing"]["kDN"])
                {
                    // If you keep kUP_v5/kDN_v5 as constexpr now, switch them to globals to override here.
                    // e.g., make them `static double kUP_v5 = 8.467, kDN_v5 = 4.256;` at file-scope,
                    // then:
                    kUP_v5 = cfg["timing"]["kUP"].as<double>();
                    kDN_v5 = cfg["timing"]["kDN"].as<double>();
                }

                // (optional) invariant-mass fit window
                if (cfg["mgg_fit"] && (cfg["mgg_fit"]["mmin"] || cfg["mgg_fit"]["mmax"]))
                {
                    // if you store mmin/mmax in your code, capture them from YAML here
                    // (define two globals or locals that your fit block will read)
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "[cfg] FAILED to load " << cfgPath << " : " << e.what() << "\n";
            }
            continue;
        }

        if (!strcmp(argv[i], "--bin-schema") && i + 1 < argc)
        {
            MR_schema_path = argv[++i];
            continue;
        }
        if (!strcmp(argv[i], "--map-run-uid") && i + 1 < argc)
        {
            MR_run_uid = argv[++i];
            continue;
        }
        if (!strcmp(argv[i], "--chunk-id") && i + 1 < argc)
        {
            MR_chunk_id = atoi(argv[++i]);
            continue;
        }
        if (!strcmp(argv[i], "--emit-microntuples") && i + 1 < argc)
        {
            MR_emit_micro = (std::string(argv[++i]) == "on");
            continue;
        }
    }
    // after parsing args and before the event loop
    const std::string MR_CONFIG_HASH = MR_compute_config_hash();

    if (!MR_schema_path.empty())
    {
        MR_schema = MR_loadBinSchema(MR_schema_path);
        if (MR_schema.valid())
            MR_acc.init(MR_schema);
    }

    if (MR_run_uid.empty())
    {
        MR_run_uid = TUUID().AsString();
    }
    // ===== End Map→Reduce CLI =====

    std::cout << "[Inputs] data=" << dataF << "  dummy=" << dummyF << "  out=" << outF << "\n";
    std::cout << "[Norm]   Q_data=" << Qdata << "  Q_dummy=" << Qdum
              << "  (kUP=" << kUP_v5 << ", kDN=" << kDN_v5 << ")\n";

    // Fill from files

    // --- mgg sidecars (2σ window and weights) ---
    const std::string mgg_json_sidecar = "mgg_roofit_sb.json";
    const std::string mgg_root_sidecar = "fit_mgg_roofit_out.root";

    // --- mgg window from JSON sidecar (preserve YAML-set nsig) ---
    {
        auto tmp = load_mgg_window_from_json(mgg_json_sidecar, g_MGW.nsig); // returns MggWindow (not MGW)
        if (tmp.ok)
        {
            // copy common fields
            g_MGW.ok = true;
            g_MGW.mu = tmp.mu;
            g_MGW.sigma = tmp.sigma;

            // prefer explicit lo/hi from sidecar if present/non-degenerate; otherwise compute from mu, sigma, nsig
            bool have_lohi = (tmp.hi > tmp.lo) && std::isfinite(tmp.lo) && std::isfinite(tmp.hi);
            if (have_lohi)
            {
                g_MGW.lo = tmp.lo;
                g_MGW.hi = tmp.hi;
            }
            else
            {
                g_MGW.lo = g_MGW.mu - g_MGW.nsig * g_MGW.sigma;
                g_MGW.hi = g_MGW.mu + g_MGW.nsig * g_MGW.sigma;
            }
        }
        else
        {
            std::cerr << "[mgg-cut] WARN: cannot read " << mgg_json_sidecar
                      << " — fallback mu=0.135 sigma=0.006 (±" << g_MGW.nsig << "σ)\n";
            g_MGW.ok = true;
            g_MGW.mu = 0.135;
            g_MGW.sigma = 0.006;
            g_MGW.lo = g_MGW.mu - g_MGW.nsig * g_MGW.sigma;
            g_MGW.hi = g_MGW.mu + g_MGW.nsig * g_MGW.sigma;
        }

        // status line
        std::cout << std::fixed << std::setprecision(6)
                  << "[mgg-cut] mu=" << g_MGW.mu << "  sigma=" << g_MGW.sigma
                  << "  window=[" << g_MGW.lo << "," << g_MGW.hi << "]"
                  << "  (mode=" << g_MGW.signal_mode
                  << ", ord=" << g_MGW.order
                  << ", rebin=" << g_MGW.rebin
                  << ", nsig=" << g_MGW.nsig << ")\n";
    }

    // weights (unchanged)
    g_MGWts = load_mgg_weights(mgg_root_sidecar);
    if (!g_MGWts.ok())
    {
        std::cerr << "[mgg-weights] WARN: no hW_sig/h_mgg_used in " << mgg_root_sidecar
                  << " — using w_sig=1 for passes.\n";
    }

    // Book a “passed ±nsig (weighted)” QA hist with your global mγγ binning
    if (!g_hMG_pass)
    {
        g_hMG_pass = new TH1D(
            "hMG_CC_evt_pass_nsig", // name: no hard-coded “2sigma”
            Form("m_{#gamma#gamma} passed #pm%.1f#sigma (weighted);M_{#gamma#gamma} (GeV);Counts",
                 g_MGW.nsig), // title reflects current nsig
            nMG, mgLo, mgHi);
        g_hMG_pass->SetDirectory(nullptr); // memory-resident is OK; you already call g_hMG_pass->Write() later
    }

    // --- CREATE & ATTACH SKIM TREE HERE (before any filling) ---
    fout.cd();
    // --- CREATE SKIM TREE (pi0 2sigma) IN MEMORY ---
    g_tSkim = new TTree("Skim", "Skim (pi0 2sigma)");
    g_tSkim->SetDirectory(nullptr); // keep it in memory; we'll rebuild a final Skim at the end
    g_tSkim->SetAutoSave(0);        // no autosave needed
    g_tSkim->SetAutoFlush(0);

    // branch backing variables (init once)
    g_skim_eventnum = 0;
    g_skim_mgg = 0.0;
    g_skim_wsig = 1.0;
    g_skim_pass_pi0 = 0;

    // define branches ONCE
    g_tSkim->Branch("eventnum", &g_skim_eventnum, "eventnum/L");
    g_tSkim->Branch("mgg", &g_skim_mgg, "mgg/D");
    g_tSkim->Branch("w_sig", &g_skim_wsig, "w_sig/D");
    g_tSkim->Branch("pass_pi0", &g_skim_pass_pi0, "pass_pi0/I");

    // ---- Now run event loops (data and dummy) ----
    Pack D("data");
    fillFromFile(dataF, "data", D);
    Pack M("dummy");
    fillFromFile(dummyF, "dummy", M);

    // Do dummy-first then Option-A for both MM and Mgg
    fout.cd();
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

    // 2D decorrelation QA (Mx^2 vs mgg) — write D and M (these are the ones we filled)
    D.h2Mx2_vs_mgg_uncorr.Write();
    D.h2Mx2_vs_mgg_corr.Write();
    D.h2Mx2_vs_mgg_CC_uncorr.Write();
    D.h2Mx2_vs_mgg_CC_corr.Write();

    M.h2Mx2_vs_mgg_uncorr.Write();
    M.h2Mx2_vs_mgg_corr.Write();
    M.h2Mx2_vs_mgg_CC_uncorr.Write();
    M.h2Mx2_vs_mgg_CC_corr.Write();

    // Final background & subtracted spectra
    OUT.hMM_Best.Write();
    OUT.hMM_Sub.Write();
    OUT.hMG_Best.Write();
    OUT.hMG_Sub.Write();

    // --- Persist what cut/weights we used this run
    {
        // Locals pulled from globals/paths already set earlier
        double mu_v = g_MGW.mu;
        double sigma_v = g_MGW.sigma;
        double lo_v = g_MGW.lo;
        double hi_v = g_MGW.hi;
        double nsig_v = g_MGW.nsig;
        int order_v = g_MGW.order;
        int rebin_v = g_MGW.rebin;
        std::string mode_v = g_MGW.signal_mode;
        std::string json = mgg_json_sidecar;
        std::string weights = mgg_root_sidecar;

        TTree *tMggCut = new TTree("MggCut", "MggCut");
        tMggCut->Branch("mu", &mu_v);
        tMggCut->Branch("sigma", &sigma_v);
        tMggCut->Branch("lo", &lo_v);
        tMggCut->Branch("hi", &hi_v);
        tMggCut->Branch("nsig", &nsig_v);
        tMggCut->Branch("order", &order_v);
        tMggCut->Branch("rebin", &rebin_v);
        tMggCut->Branch("mode", &mode_v);
        tMggCut->Branch("json_path", &json);
        tMggCut->Branch("weights_path", &weights);
        tMggCut->Fill();
        tMggCut->Write();
    }

    // Write the new QA histogram and skim (if they were created)
    if (g_hMG_pass)
        g_hMG_pass->Write();
    if (g_tSkim)
        g_tSkim->Write();

    // ================== Map→Reduce: write per-run per-bin table ==================
    if (MR_schema.valid())
    {
        // Recompute the same geometry scalars used in doDummyThenA() so algebra matches
        const double Wsig = (C_LO < C_HI) ? (C_HI - C_LO) : 0.0;
        double Wpos = 0.0, Wneg = 0.0;
        for (int i = 0; i < NPOS; ++i)
            Wpos += (posWins[i].hi - posWins[i].lo);
        for (int i = 0; i < NNEG; ++i)
            Wneg += (negWins[i].hi - negWins[i].lo);

        const double ACC = Wsig; // your code uses this to normalize stripes
        const double aV = (Wpos + Wneg) > 0.0 ? (Wsig / (Wpos + Wneg)) : 0.0;
        const double aH = aV;

        // A_AD (diagonals) and A_AP (anti-diagonals) geometry weights (match doDummyThenA)
        auto sumSquares = [&](const Win *w, int n)
        {
            double s = 0.0;
            for (int i = 0; i < n; i++)
                s += (w[i].hi - w[i].lo) * (w[i].hi - w[i].lo);
            return s;
        };
        const double sumSqPos = sumSquares(posWins, NPOS);
        const double sumSqNeg = sumSquares(negWins, NNEG);
        const double A_AD = sumSqPos + sumSqNeg; // same-stripe squares (LL + UR)
        const double A_AP = 2.0 * Wpos * Wneg;   // off-diagonal corners (UL + LR)

        const double aAD = (A_AD > 0) ? (ACC / A_AD) : 0.0;
        const double aAP = (A_AP > 0) ? (ACC / A_AP) : 0.0;

        // Build TTREE
        fout.cd();
        TTree tb("bins", "per-run per-bin yields (dummy-first, Option-A)");

        int bin_id = -1;
        int iMx = -1, iMg = -1;
        double edgeMx_lo = 0, edgeMx_hi = 0, edgeMg_lo = 0, edgeMg_hi = 0;
        double SumW = 0, SumW2 = 0;
        double SumW_data = 0, SumW_dummyUP = 0, SumW_dummyDN = 0, SumW_acc = 0;

        tb.Branch("bin_id", &bin_id, "bin_id/I");
        if (MR_schema.axisNames.size() >= 1)
        {
            tb.Branch("iMx", &iMx, "iMx/I");
            tb.Branch("edgeMx_lo", &edgeMx_lo, "edgeMx_lo/D");
            tb.Branch("edgeMx_hi", &edgeMx_hi, "edgeMx_hi/D");
        }
        if (MR_schema.axisNames.size() >= 2)
        {
            tb.Branch("iMg", &iMg, "iMg/I");
            tb.Branch("edgeMg_lo", &edgeMg_lo, "edgeMg_lo/D");
            tb.Branch("edgeMg_hi", &edgeMg_hi, "edgeMg_hi/D");
        }
        tb.Branch("SumW", &SumW, "SumW/D");
        tb.Branch("SumW2", &SumW2, "SumW2/D");
        tb.Branch("SumW_data", &SumW_data, "SumW_data/D");
        tb.Branch("SumW_dummyUP", &SumW_dummyUP, "SumW_dummyUP/D");
        tb.Branch("SumW_dummyDN", &SumW_dummyDN, "SumW_dummyDN/D");
        tb.Branch("SumW_acc", &SumW_acc, "SumW_acc/D");

        // Store BinEdges/ once for convenience
        TDirectory *dBE = fout.GetDirectory("BinEdges");
        if (!dBE)
            dBE = fout.mkdir("BinEdges");
        dBE->cd();
        for (size_t ax = 0; ax < MR_schema.axisNames.size(); ++ax)
        {
            TVectorD v((int)MR_schema.edges[ax].size());
            for (int i = 0; i < v.GetNoElements(); ++i)
                v[i] = MR_schema.edges[ax][i];
            v.Write(MR_schema.axisNames[ax].c_str(), TObject::kOverwrite);
        }
        fout.cd();

        // Algebra per bin (dummy first → A-method), using same factors as in doDummyThenA()
        int nx = (MR_schema.edges[0].size() > 0) ? (int)MR_schema.edges[0].size() - 1 : 1;

        for (size_t b = 0; b < MR_acc.bins.size(); ++b)
        {
            const MR_BinRow &r = MR_acc.bins[b];

            // 1) Dummy normalization per category, scaled from dummy→data exposure
            auto norm_dummy = [&](double up, double dn) -> double
            {
                if (Qdum <= 0)
                    return 0.0;
                double cup = 1.0 / (Qdum * kUP_v5);
                double cdn = 1.0 / (Qdum * kDN_v5);
                double per_uC = up * cup + dn * cdn; // per unit charge
                return per_uC * Qdata;               // scale to data charge seen
            };

            double CC_afterDummy = r.CC_data - norm_dummy(r.CC_up, r.CC_dn);
            double V_afterDummy = r.V_data - norm_dummy(r.V_up, r.V_dn);
            double H_afterDummy = r.H_data - norm_dummy(r.H_up, r.H_dn);
            double AD_afterDummy = r.AD_data - norm_dummy(r.AD_up, r.AD_dn);
            double AP_afterDummy = r.AP_data - norm_dummy(r.AP_up, r.AP_dn);

            // 2) A-method Best (exactly mirrors your doDummyThenA combination)
            double Best = 0.5 * aV * V_afterDummy + 0.5 * aH * H_afterDummy + aAD * AD_afterDummy - aAP * AP_afterDummy;

            double FinalSUB = CC_afterDummy - Best;

            // 3) Variances (Poisson for counts; include scale factors^2 for dummy)
            auto var_dummy = [&](double up, double dn) -> double
            {
                if (Qdum <= 0)
                    return 0.0;
                double cup = 1.0 / (Qdum * kUP_v5);
                double cdn = 1.0 / (Qdum * kDN_v5);
                double var_per_uC = up * cup * cup + dn * cdn * cdn;
                return var_per_uC * Qdata * Qdata;
            };

            double var_CC_afterDummy = r.CC_data + var_dummy(r.CC_up, r.CC_dn);
            double var_V_afterDummy = r.V_data + var_dummy(r.V_up, r.V_dn);
            double var_H_afterDummy = r.H_data + var_dummy(r.H_up, r.H_dn);
            double var_AD_afterDummy = r.AD_data + var_dummy(r.AD_up, r.AD_dn);
            double var_AP_afterDummy = r.AP_data + var_dummy(r.AP_up, r.AP_dn);

            double var_Best = (0.5 * aV) * (0.5 * aV) * var_V_afterDummy + (0.5 * aH) * (0.5 * aH) * var_H_afterDummy + (aAD) * (aAD)*var_AD_afterDummy + (-aAP) * (-aAP) * var_AP_afterDummy; // subtractive term

            double var_Final = var_CC_afterDummy + var_Best; // CC - Best ⇒ add variances

            // Fill row
            bin_id = (int)b;
            if (MR_schema.axisNames.size() == 1)
            {
                iMx = bin_id;
                edgeMx_lo = MR_schema.edges[0][iMx];
                edgeMx_hi = MR_schema.edges[0][iMx + 1];
            }
            else if (MR_schema.axisNames.size() == 2)
            {
                iMx = (int)(b % nx);
                iMg = (int)(b / nx);
                edgeMx_lo = MR_schema.edges[0][iMx];
                edgeMx_hi = MR_schema.edges[0][iMx + 1];
                edgeMg_lo = MR_schema.edges[1][iMg];
                edgeMg_hi = MR_schema.edges[1][iMg + 1];
            }

            SumW = FinalSUB;
            SumW2 = std::max(0.0, var_Final);

            // Audits (compact)
            SumW_data = r.CC_data;
            SumW_dummyUP = r.CC_up;
            SumW_dummyDN = r.CC_dn;
            SumW_acc = Best;

            tb.Fill();
        }

        // Metadata for reducer guardrails
        TNamed("bin_schema_id", MR_schema.id.c_str()).Write();
        TNamed("config_hash", MR_CONFIG_HASH.c_str()).Write();
        TNamed("run_uid", MR_run_uid.c_str()).Write();
        TParameter<int>("chunk_id", MR_chunk_id).Write("chunk_id");
        TParameter<double>("Q_data", Qdata).Write("Q_data");
        TParameter<double>("Q_dummy", Qdum).Write("Q_dummy");

        tb.Write("", TObject::kOverwrite);
        std::cerr << "[bins] wrote per-run table with " << MR_acc.bins.size() << " bins\n";
    }
    // =============================================================================

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
        std::cerr << "[doDummyThenA] wrote Mgg_overlay_final\n";

        c.Write("Mgg_overlay_final");

        // Persist the final Mgg B_est and SUB histograms for QA / downstream checks
        fout.cd();
        OUT.hMG_Best.Write(); // name is like hMG_Best_final
        OUT.hMG_Sub.Write();  // name is like hMG_Sub_final

        // --- write diagnostics added earlier at file scope ---
        fout.cd();
        hDeltaM_vs_CorrFac.Write();
        hMxCorr_vs_MxRaw.Write();
        if (g_hMG_pass)
            g_hMG_pass->Write(); // QA hist (optional)

        // --- Build final Skim tree with per-event physics weights ---
        if (g_tSkim && g_hMG_pass)
        {
            // hSub: final Mgg after dummy + Option-A (already charge-normalized)
            TH1 *hSub = &OUT.hMG_Sub;
            TH1 *hPass = g_hMG_pass; // mgg of passing events, weighted by w_sig

            const int nbins = hSub->GetNbinsX();
            std::vector<double> s_dummyA(nbins + 1, 0.0); // indices 1..nbins

            // Per-bin scaling so that sum_i w_phys_i in bin j = hSub(j)
            for (int ib = 1; ib <= nbins; ++ib)
            {
                const double num = hSub->GetBinContent(ib);
                const double den = hPass->GetBinContent(ib);
                s_dummyA[ib] = (den > 0.0) ? (num / den) : 0.0;
            }

            // Input skim: eventnum, mgg, w_sig, pass_pi0
            Long64_t in_evnum = 0;
            double in_mgg = 0.0;
            double in_wsig = 0.0;
            int in_pass = 0;

            g_tSkim->SetBranchStatus("*", 0);
            g_tSkim->SetBranchStatus("eventnum", 1);
            g_tSkim->SetBranchStatus("mgg", 1);
            g_tSkim->SetBranchStatus("w_sig", 1);
            g_tSkim->SetBranchStatus("pass_pi0", 1);

            g_tSkim->SetBranchAddress("eventnum", &in_evnum);
            g_tSkim->SetBranchAddress("mgg", &in_mgg);
            g_tSkim->SetBranchAddress("w_sig", &in_wsig);
            g_tSkim->SetBranchAddress("pass_pi0", &in_pass);

            // Output skim written to the file: adds w_phys
            TTree *tSkimOut = new TTree("Skim", "Skim (pi0 2sigma, dummy+OptionA physics weights)");
            tSkimOut->SetDirectory(&fout);

            Long64_t out_evnum = 0;
            double out_mgg = 0.0;
            double out_wsig = 0.0;
            double out_wphys = 0.0;
            int out_pass = 0;
            double out_wsignal = 0.0;

            tSkimOut->Branch("eventnum", &out_evnum, "eventnum/l");
            tSkimOut->Branch("mgg", &out_mgg, "mgg/D");
            tSkimOut->Branch("w_sig", &out_wsig, "w_sig/D");
            tSkimOut->Branch("w_phys", &out_wphys, "w_phys/D");
            tSkimOut->Branch("pass_pi0", &out_pass, "pass_pi0/I");
            tSkimOut->Branch("w_signal", &out_wsignal, "w_signal/D");

            const Long64_t nentries = g_tSkim->GetEntries();
            for (Long64_t i = 0; i < nentries; ++i)
            {
                g_tSkim->GetEntry(i);

                out_evnum = in_evnum;
                out_mgg = in_mgg;
                out_wsig = in_wsig;
                out_pass = in_pass;

                // Map mgg to Mgg bin in OUT.hMG_Sub
                int ib = hSub->GetXaxis()->FindBin(in_mgg);
                if (ib < 1)
                    ib = 1;
                if (ib > nbins)
                    ib = nbins;

                const double sDA = s_dummyA[ib];
                out_wphys = out_wsig * sDA; // per-event physics weight: dummy+Option-A applied

                const double fracS = (ib < (int)g_MGW.Sfrac_bins.size())
                                         ? g_MGW.Sfrac_bins[ib]
                                         : 0.0;

                out_wsignal = out_wphys * fracS; // final S+B–weighted signal-only weight

                tSkimOut->Fill();
            }

            tSkimOut->Write("", TObject::kOverwrite);
        }
        else if (g_tSkim)
        {
            // Fallback: if QA hist is missing, at least persist the raw skim
            g_tSkim->SetDirectory(&fout);
            g_tSkim->Write();
        }
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
