// fit_mgg_minimal.C
// Minimal background-subtract + Gaussian fit on mgg_final_sub_data.
//
// Usage (batch):
//   root -l -b -q 'fit_mgg_minimal.C("your_output.root")'
//
// Optional args:
//   root -l -b -q 'fit_mgg_minimal.C("your_output.root",
//                                     "mgg_overlay_final/mgg_final_sub_data",
//                                     0.10, 0.17,   // fit range
//                                     0.100,0.115,  // left sideband
//                                     0.155,0.170,  // right sideband
//                                     2,            // background order (2=quadratic)
//                                     true)'        // save QA ROOT/PNGs
//
// Output (stdout):
//   One-line summary with mu, sigma, yield, chi2/ndf, negatives-in-fit, etc.
//
// If saveQA==true, writes 'mgg_fit_minimal_out.root' with:
//   hData : original input hist (clone)
//   grSB  : sideband points used in the bkg fit
//   fBkg  : fitted background (TF1)
//   hSub  : data - bkg (in [fit_min, fit_max])
//   fGaus : Gaussian fit to residual
//   plus two PNGs: mgg_bkgfit.png, mgg_sub_gauss.png
//
//   g++ -std=c++17 fit_mgg_minimal_sb.C     -lyaml-cpp \
//    `root-config --cflags --libs` -lTMVA -lRooFitCore -lRooFit -lRIO \
//    -o fit_mgg_minimal
//

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>

#include <yaml-cpp/yaml.h>

static TH1* fetchHist(TFile* f, const std::string& path) {
  if (!f) return nullptr;
  TObject* o = f->Get(path.c_str());
  if (o) {
    TH1* h = dynamic_cast<TH1*>(o);
    if (h) return h;
  }
  // try splitting "dir/hist"
  auto pos = path.rfind('/');
  if (pos != std::string::npos) {
    std::string dname = path.substr(0, pos);
    std::string hname = path.substr(pos+1);
    TDirectory* d = f->GetDirectory(dname.c_str());
    if (d) {
      TObject* o2 = d->Get(hname.c_str());
      if (o2) return dynamic_cast<TH1*>(o2);
    }
  }
  return nullptr;
}

void fit_mgg_minimal(const char* infile,
                     const char* histpath="hMG_CC_evt_data",
                     double fit_min=0.10, double fit_max=0.17,
                     double sbL1=0.100, double sbL2=0.115,
                     double sbR1=0.155, double sbR2=0.170,
                     int bkg_order=2,
                     bool saveQA=true)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // Open file and fetch histogram
  std::unique_ptr<TFile> f(TFile::Open(infile, "READ"));
  if (!f || f->IsZombie()) {
    std::cerr << "[ERR] Cannot open file: " << infile << "\n";
    return;
  }
  TH1* hIn = fetchHist(f.get(), histpath);
  if (!hIn) {
    std::cerr << "[ERR] Cannot find histogram: " << histpath << "\n";
    return;
  }

  // Clone as TH1D, ensure Sumw2
  TH1D* hData = dynamic_cast<TH1D*>(hIn->Clone("hData"));
  if (!hData) {
    hData = new TH1D("hData", hIn->GetTitle(),
                     hIn->GetNbinsX(),
                     hIn->GetXaxis()->GetXmin(),
                     hIn->GetXaxis()->GetXmax());
    for (int i=1;i<=hIn->GetNbinsX();++i) {
      hData->SetBinContent(i, hIn->GetBinContent(i));
      hData->SetBinError  (i, hIn->GetBinError(i));
    }
  }
  hData->SetDirectory(nullptr);
  hData->Sumw2();
  hData->GetXaxis()->SetTitle("m_{#gamma#gamma}  [GeV]");
  hData->GetYaxis()->SetTitle("Counts");

  // Sideband graph (points with errors)
  std::vector<double> vx, vy, vex, vey;
  const int nb = hData->GetNbinsX();
  for (int i=1;i<=nb;++i) {
    double x  = hData->GetXaxis()->GetBinCenter(i);
    double y  = hData->GetBinContent(i);
    double ey = hData->GetBinError(i);
    bool inLeft  = (x>=sbL1 && x<=sbL2);
    bool inRight = (x>=sbR1 && x<=sbR2);
    if (inLeft || inRight) {
      vx.push_back(x); vy.push_back(y);
      vex.push_back(0.0); vey.push_back((ey>0 ? ey : TMath::Sqrt(std::max(0.0,y))));
    }
  }
  if (vx.size() < (size_t)(bkg_order+1)) {
    std::cerr << "[ERR] Not enough sideband points for polynomial order " << bkg_order
              << " (points=" << vx.size() << ")\n";
    return;
  }
  auto grSB = std::make_unique<TGraphErrors>((int)vx.size(), vx.data(), vy.data(), vex.data(), vey.data());
  grSB->SetName("grSB");
  grSB->SetTitle("Sideband points for background fit");

  // Background TF1
  std::unique_ptr<TF1> fBkg;
  if (bkg_order==1)      fBkg.reset(new TF1("fBkg","pol1", fit_min, fit_max));
  else if (bkg_order==2) fBkg.reset(new TF1("fBkg","pol2", fit_min, fit_max));
  else if (bkg_order==3) fBkg.reset(new TF1("fBkg","pol3", fit_min, fit_max));
  else {
    std::cerr << "[WARN] Unsupported bkg_order=" << bkg_order << ", forcing 2.\n";
    fBkg.reset(new TF1("fBkg","pol2", fit_min, fit_max));
  }

  // Fit sidebands
  grSB->Fit(fBkg.get(), "Q0"); // quiet, no draw
  // Build background histogram evaluated at bin centers
  auto hBkg = std::unique_ptr<TH1D>((TH1D*)hData->Clone("hBkg"));
  hBkg->Reset("ICESM");
  hBkg->Sumw2();
  for (int i=1;i<=nb;++i) {
    double x  = hData->GetXaxis()->GetBinCenter(i);
    double bv = fBkg->Eval(x);
    hBkg->SetBinContent(i, bv);
    // Error on background model is not trivial; keep zero here (fit captures uncertainty implicitly via χ²)
    hBkg->SetBinError(i, 0.0);
  }

  // Subtract background to get residual S_hat
  auto hSub = std::unique_ptr<TH1D>((TH1D*)hData->Clone("hSub"));
  hSub->Add(hBkg.get(), -1.0);
  hSub->Sumw2();
  hSub->GetYaxis()->SetTitle("Counts (subtracted)");

  // Count negatives within fit range
  int nfitbins=0, nneg=0;
  for (int i=1;i<=nb;++i) {
    double x = hSub->GetXaxis()->GetBinCenter(i);
    if (x>=fit_min && x<=fit_max) {
      ++nfitbins;
      if (hSub->GetBinContent(i) <= 0) ++nneg;
    }
  }

  // Gaussian fit in [fit_min, fit_max] (binned chi2)
  TF1 fGaus("fGaus","gaus", fit_min, fit_max);
  // Crude seeds from sub-moment in window
  double mean_est = 0.135;
  double rms_est  = 0.010;
  // Use data around the peak to refine seed
  double wsum=0, xsum=0, x2sum=0;
  for (int i=1;i<=nb;++i) {
    double x = hSub->GetXaxis()->GetBinCenter(i);
    if (x<fit_min || x>fit_max) continue;
    double y = std::max(0.0, (double)hSub->GetBinContent(i)); // ignore negatives for seeding
    wsum += y; xsum += x*y; x2sum += x*x*y;
  }
  if (wsum>0) {
    mean_est = xsum/wsum;
    double var = x2sum/wsum - mean_est*mean_est;
    if (var>1e-6) rms_est = TMath::Sqrt(var);
  }
  double amp_est = hSub->GetBinContent(hSub->GetXaxis()->FindBin(mean_est));
  fGaus.SetParameters(std::max(1.0, amp_est), mean_est, std::min(std::max(rms_est, 0.003), 0.03));
  fGaus.SetParNames("A","mu","sigma");

  int fitStatus = hSub->Fit(&fGaus, "RQ0"); // R: use range, Q: quiet, 0: no draw
  double A     = fGaus.GetParameter(0);
  double mu    = fGaus.GetParameter(1);
  double sigma = fGaus.GetParameter(2);
  double Aerr     = fGaus.GetParError(0);
  double muerr    = fGaus.GetParError(1);
  double sigmaerr = fGaus.GetParError(2);
  double chi2     = fGaus.GetChisquare();
  int    ndf      = fGaus.GetNDF();

  // Yield inside fit window (analytic integral of Gaussian from a to b)
  auto gauss_int = [](double A, double mu, double sig, double a, double b){
    const double rt2 = TMath::Sqrt2();
    double z1 = (a - mu)/(sig*rt2);
    double z2 = (b - mu)/(sig*rt2);
    return A * sig * TMath::Sqrt(2*TMath::Pi()) * 0.5 * (TMath::Erf(z2) - TMath::Erf(z1));
  };
  double yield_fit = gauss_int(A, mu, sigma, fit_min, fit_max);

  // One-line summary
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "[mgg-fit] range=[" << fit_min << "," << fit_max << "]  "
            << "SB(L=[" << sbL1 << "," << sbL2 << "], R=[" << sbR1 << "," << sbR2 << "])  "
            << "bkgOrder=" << bkg_order << "  nfitbins=" << nfitbins << "  neg_in_fit=" << nneg
            << "  mu=" << mu << "±" << muerr
            << "  sigma=" << sigma << "±" << sigmaerr
            << "  A=" << A << "±" << Aerr
            << "  chi2/ndf=" << (ndf>0 ? chi2/ndf : -1)
            << "  yield[" << fit_min << "," << fit_max << "]=" << yield_fit
            << "  status=" << fitStatus
            << "\n";

  if (!saveQA) return;

  // QA plots
  TCanvas c1("c1","mgg background fit",900,700);
  c1.SetMargin(0.12,0.04,0.12,0.05);
  hData->SetLineWidth(2);
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(0.9);
  hData->Draw("E1");
  fBkg->SetLineColor(kBlue+1);
  fBkg->SetLineWidth(3);
  fBkg->Draw("SAME");
  grSB->SetMarkerStyle(24);
  grSB->SetMarkerColor(kBlue+3);
  grSB->Draw("P SAME");
  {
    TLegend leg(0.55,0.73,0.88,0.88);
    leg.AddEntry(hData, "Data (selected pairs)", "lep");
    leg.AddEntry(fBkg.get(), "Sideband bkg fit (poly)", "l");
    leg.AddEntry(grSB.get(), "Sideband points", "p");
    leg.Draw();
  }
  c1.SaveAs("mgg_bkgfit.png");

  TCanvas c2("c2","subtracted & gaussian",900,700);
  c2.SetMargin(0.12,0.04,0.12,0.05);
  hSub->SetLineWidth(2);
  hSub->SetMarkerStyle(20);
  hSub->SetMarkerSize(0.9);
  hSub->Draw("E1");
  fGaus.SetLineColor(kRed+1);
  fGaus.SetLineWidth(3);
  fGaus.Draw("SAME");
  {
    TLegend leg(0.55,0.73,0.88,0.88);
    leg.AddEntry(hSub.get(), "Data - bkg (residual)", "lep");
    leg.AddEntry(&fGaus, "Gaussian fit", "l");
    leg.Draw();
  }
  c2.SaveAs("mgg_sub_gauss.png");

  // Save QA objects
  std::unique_ptr<TFile> fout(TFile::Open("mgg_fit_minimal_out.root","RECREATE"));
  if (fout && !fout->IsZombie()) {
    hData->Write("hData");
    grSB->Write("grSB");
    fBkg->Write("fBkg");
    hSub->Write("hSub");
    fGaus.Write("fGaus");
    fout->Write();
    fout->Close();
  }
}

// -------------------------------------------------------------
//  YAML-CPP driver
// -------------------------------------------------------------
int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] 
                  << " config.yaml\n";
        return 1;
    }

    std::string yamlfile = argv[1];

    YAML::Node cfg;
    try {
        cfg = YAML::LoadFile(yamlfile);
    } catch (const YAML::Exception& e) {
        std::cerr << "Error reading YAML config: " 
                  << e.what() << "\n";
        return 1;
    }

    // Required
    std::string infile  = cfg["infile"].as<std::string>();
    std::string histpath = cfg["histpath"].as<std::string>();

    // Fit range
    double fit_min = cfg["fit_range"]["min"].as<double>();
    double fit_max = cfg["fit_range"]["max"].as<double>();

    // Sidebands
    double sbL1 = cfg["sidebands"]["left"]["min"].as<double>();
    double sbL2 = cfg["sidebands"]["left"]["max"].as<double>();
    double sbR1 = cfg["sidebands"]["right"]["min"].as<double>();
    double sbR2 = cfg["sidebands"]["right"]["max"].as<double>();

    // Background polynomial order
    int bkg_order = cfg["background_order"].as<int>();

    // Save QA plots?
    bool saveQA = cfg["save_QA"].as<bool>();

    // ---------------------------------------------------------
    //  Run the user function (no changes needed)
    // ---------------------------------------------------------
    fit_mgg_minimal(infile.c_str(),
                    histpath.c_str(),
                    fit_min, fit_max,
                    sbL1, sbL2,
                    sbR1, sbR2,
                    bkg_order,
                    saveQA);

    return 0;
}