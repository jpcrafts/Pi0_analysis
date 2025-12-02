// fit_mgg_roofit_sb.C
// Unsubtracted S+B fit on mgg with Bernstein background (order ord_min..ord_max).
// Signal options: "gauss" (default), "doubleG" (shared μ), or "cb" (CrystalBall).
//
// Usage (batch):
//   root -l -b -q 'fit_mgg_roofit_sb.C("4205_v5_slim_pairsA_modified_mgg.root",
//                                      "hMG_CC_evt_data", 0.08, 0.18,
//                                      true, 2, 2, "cb")'
//
// Outputs (if saveQA):
//   fit_mgg_roofit_out.root : RooFitResult, RooPlot, chosen model snapshot
//   mgg_roofit_sb.png       : plot with data/model/components
//   mgg_roofit_sb.json      : μ, σ (or σ_eff / CB σ), Ns/Nb, order, AICc, χ²/ndf, etc.
//
//   g++ -std=c++17 -O2     fit_mgg_roofit_sb.C     -lyaml-cpp     `root-config --cflags --libs` -lTMVA -lRooFitCore -lRooFit    -o fit_mgg_roofit_sb
//

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <TString.h>
#include <TLegend.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TLine.h>

#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooCurve.h>
#include <RooBernstein.h>
#include <RooCBShape.h> // Crystal Ball

#include <memory>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include <yaml-cpp/yaml.h>

using namespace RooFit;

static TH1 *fetchHist(TFile *f, const std::string &path)
{
    if (!f)
        return nullptr;
    if (auto *o = f->Get(path.c_str()))
        return dynamic_cast<TH1 *>(o);
    auto pos = path.rfind('/');
    if (pos != std::string::npos)
    {
        std::string dname = path.substr(0, pos);
        std::string hname = path.substr(pos + 1);
        if (auto *d = f->GetDirectory(dname.c_str()))
        {
            if (auto *o2 = d->Get(hname.c_str()))
                return dynamic_cast<TH1 *>(o2);
        }
    }
    return nullptr;
}

void fit_mgg_roofit_sb(const char *infile,
                       const char *histpath = "hMG_CC_evt_data",
                       double mmin = 0.08,
                       double mmax = 0.18,
                       bool saveQA = true,
                       int ord_min = 1,
                       int ord_max = 3,
                       const char *signal_mode = "gauss",
                       int rebin = 1) // <--- NEW

{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    std::unique_ptr<TFile> fin(TFile::Open(infile, "READ"));
    if (!fin || fin->IsZombie())
    {
        printf("[ERR] Cannot open %s\n", infile);
        return;
    }
    TH1 *hIn = fetchHist(fin.get(), histpath);
    if (!hIn)
    {
        printf("[ERR] Cannot find hist: %s\n", histpath);
        fin->ls();
        return;
    }

    TH1D *h = dynamic_cast<TH1D *>(hIn->Clone("h_mgg_in"));
    if (!h)
    {
        h = new TH1D("h_mgg_in", hIn->GetTitle(), hIn->GetNbinsX(),
                     hIn->GetXaxis()->GetXmin(), hIn->GetXaxis()->GetXmax());
        for (int i = 1; i <= hIn->GetNbinsX(); ++i)
        {
            h->SetBinContent(i, hIn->GetBinContent(i));
            h->SetBinError(i, hIn->GetBinError(i));
        }
    }
    h->SetDirectory(nullptr);
    h->Sumw2();

    // Optional rebin before building RooDataHist
    if (rebin > 1)
    {
        TH1 *hTmp = h->Rebin(rebin, "h_mgg_in_rebinned"); // returns TH1*
        TH1D *hReb = dynamic_cast<TH1D *>(hTmp);
        if (!hReb)
        {
            // If ROOT returned a non-TH1D, rebuild as TH1D with same binning
            hReb = new TH1D("h_mgg_in_rebinned", hTmp->GetTitle(),
                            hTmp->GetNbinsX(),
                            hTmp->GetXaxis()->GetXmin(),
                            hTmp->GetXaxis()->GetXmax());
            hReb->Sumw2();
            for (int i = 1; i <= hTmp->GetNbinsX(); ++i)
            {
                hReb->SetBinContent(i, hTmp->GetBinContent(i));
                hReb->SetBinError(i, hTmp->GetBinError(i));
            }
        }
        h = hReb; // take ownership of the rebinned histogram
        h->SetDirectory(nullptr);
        h->Sumw2();
    }

    RooRealVar m("m", "m_{#gamma#gamma} [GeV]", mmin, mmax);

    int firstBin = h->GetXaxis()->FindBin(mmin + 1e-6);
    int lastBin = h->GetXaxis()->FindBin(mmax - 1e-6);
    int nbinsInRange = lastBin - firstBin + 1;
    if (nbinsInRange < 20)
        nbinsInRange = std::max(20, (int)((mmax - mmin) / 0.001));
    m.setBins(nbinsInRange);

    std::unique_ptr<TH1D> hRange(new TH1D("h_mgg_used", h->GetTitle(), nbinsInRange, mmin, mmax));
    hRange->Sumw2();
    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
        double x = h->GetXaxis()->GetBinCenter(i);
        if (x < mmin || x > mmax)
            continue;
        int j = hRange->FindBin(x);
        hRange->SetBinContent(j, h->GetBinContent(i));
        hRange->SetBinError(j, h->GetBinError(i));
    }

    RooDataHist data("data", "mgg data", RooArgList(m), hRange.get());
    int N_nonzero = 0;
    for (int i = 1; i <= hRange->GetNbinsX(); ++i)
        if (hRange->GetBinContent(i) > 0)
            ++N_nonzero;
    int N_for_AICc = (N_nonzero > 0 ? N_nonzero : hRange->GetNbinsX());

    struct FitPick
    {
        int order = -1, status = -999;
        double aicc = 1e300, chi2ndf = 1e300, nll = 1e300;
        double mu = 0, muerr = 0;
        double sig = 0, sigerr = 0;                                                 // σ for Gauss/CB, or σ_eff for doubleG
        double sig1 = 0, sig1err = 0, sig2 = 0, sig2err = 0, frac = 0, fracerr = 0; // doubleG details
        double alpha = 0, alphaerr = 0, n = 0, nerr = 0;                            // CB details
        double Ns = 0, Nserr = 0, Nb = 0, Nberr = 0;
        RooFitResult *res = nullptr;
        RooAddPdf *model = nullptr;
        RooPlot *frame = nullptr;
        std::string smode = "gauss";
    } best;

    if (ord_min < 1)
        ord_min = 1;
    if (ord_max < ord_min)
        ord_max = ord_min;
    std::string smode_in(signal_mode);

    for (int k = ord_min; k <= ord_max; ++k)
    {
        // --- Signal
        RooRealVar mu("mu", "mu", 0.135, 0.12, 0.15);

        std::unique_ptr<RooAbsPdf> sig_pdf;
        std::unique_ptr<RooRealVar> sigma, sigma1, sigma2, frac, alpha, n;
        std::unique_ptr<RooGaussian> g1, g2;
        std::unique_ptr<RooCBShape> cb;

        if (smode_in == "doubleG")
        {
            sigma1.reset(new RooRealVar("sigma1", "sigma1", 0.006, 0.003, 0.020));
            sigma2.reset(new RooRealVar("sigma2", "sigma2", 0.012, 0.004, 0.040));
            frac.reset(new RooRealVar("f", "f", 0.60, 0.00, 1.00));
            g1.reset(new RooGaussian("g1", "G1", m, mu, *sigma1));
            g2.reset(new RooGaussian("g2", "G2", m, mu, *sigma2));
            sig_pdf.reset(new RooAddPdf("sig", "f*G1+(1-f)*G2", RooArgList(*g1, *g2), RooArgList(*frac)));
        }
        else if (smode_in == "cb")
        {
            sigma.reset(new RooRealVar("sigma", "sigma", 0.010, 0.003, 0.030));
            alpha.reset(new RooRealVar("alpha", "alpha", 1.5, 0.1, 5.0)); // tail
            n.reset(new RooRealVar("n", "n", 5.0, 1.1, 20.0));            // power
            cb.reset(new RooCBShape("sig", "CrystalBall", m, mu, *sigma, *alpha, *n));
            sig_pdf.reset((RooAbsPdf *)cb->clone("sig"));
        }
        else
        { // "gauss"
            sigma.reset(new RooRealVar("sigma", "sigma", 0.010, 0.003, 0.020));
            g1.reset(new RooGaussian("g1", "Gaussian", m, mu, *sigma));
            sig_pdf.reset((RooAbsPdf *)g1->clone("sig"));
        }

        // --- Background: Bernstein(k)
        RooArgList coeffs;
        std::vector<std::unique_ptr<RooRealVar>> cstore;
        cstore.reserve(k + 1);
        for (int i = 0; i <= k; i++)
        {
            TString nm;
            nm.Form("c%d", i);
            cstore.emplace_back(new RooRealVar(nm, nm, 0.1, 0.0, 1e6));
            coeffs.add(*cstore.back());
        }
        RooBernstein bkg("bkg", "Bernstein background", m, coeffs);

        // --- Extended yields
        double sumY = hRange->Integral();
        RooRealVar Ns("Ns", "signal yield", 0.5 * sumY, 0.0, 10.0 * sumY + 1.0);
        RooRealVar Nb("Nb", "bkg yield", 0.5 * sumY, 0.0, 10.0 * sumY + 1.0);

        RooAddPdf model("model", "sig+bkg", RooArgList(*sig_pdf, bkg), RooArgList(Ns, Nb));

        std::unique_ptr<RooFitResult> res(model.fitTo(data, Extended(true), Save(true),
                                                      PrintLevel(-1), Warnings(false),
                                                      SumW2Error(true)));

        int status = (res ? res->status() : 999);
        double nll = (res ? res->minNll() : 1e300);

        // Parameter count K for AICc
        // base: Ns, Nb (2) + mu (1) + background (k+1)
        // gauss adds sigma (1); doubleG adds sigma1, sigma2, frac (3); cb adds sigma, alpha, n (3)
        int add_sig = (smode_in == "doubleG") ? 3 : ((smode_in == "cb") ? 3 : 1);
        int K = 2 + 1 + (k + 1) + add_sig;

        double AIC = 2.0 * K + 2.0 * nll;
        double AICc = AIC;
        if (N_for_AICc > (K + 1))
            AICc += (2.0 * K * (K + 1)) / ((double)N_for_AICc - K - 1.0);

        // Plot & chi2/ndf
        std::unique_ptr<RooPlot> fr(m.frame(Bins(nbinsInRange), Title("m_{#gamma#gamma} fit")));
        data.plotOn(fr.get(), DataError(RooAbsData::SumW2));
        model.plotOn(fr.get(), Components("bkg"), LineStyle(kDashed), LineColor(kBlue + 1));
        model.plotOn(fr.get(), Components("sig"), LineStyle(kDotted), LineColor(kRed + 1));
        model.plotOn(fr.get());
        double chi2 = fr->chiSquare();

        // Extract params
        double mu_v = mu.getVal(), mu_e = mu.getError();
        double Ns_v = Ns.getVal(), Ns_e = Ns.getError();
        double Nb_v = Nb.getVal(), Nb_e = Nb.getError();

        double sig_report = 0, sig_err = 0;
        double s1_v = 0, s1_e = 0, s2_v = 0, s2_e = 0, f_v = 0, f_e = 0, a_v = 0, a_e = 0, n_v = 0, n_e = 0;

        bool sane = (status == 0 || status == 1) && (mu_v >= 0.12 && mu_v <= 0.15) && (Ns_v > 0);

        if (smode_in == "doubleG")
        {
            s1_v = sigma1->getVal();
            s1_e = sigma1->getError();
            s2_v = sigma2->getVal();
            s2_e = sigma2->getError();
            f_v = frac->getVal();
            f_e = frac->getError();
            sig_report = std::sqrt(std::max(0.0, f_v * s1_v * s1_v + (1.0 - f_v) * s2_v * s2_v)); // σ_eff
            sig_err = 0.5 * (s1_e + s2_e);
            sane = sane && (s1_v >= 0.003 && s1_v <= 0.040) && (s2_v >= 0.003 && s2_v <= 0.080) &&
                   (sig_report >= 0.003 && sig_report <= 0.030);
        }
        else if (smode_in == "cb")
        {
            sig_report = sigma->getVal();
            sig_err = sigma->getError();
            a_v = alpha->getVal();
            a_e = alpha->getError();
            n_v = n->getVal();
            n_e = n->getError();
            sane = sane && (sig_report >= 0.003 && sig_report <= 0.030) &&
                   (a_v > 0.05 && a_v < 8.0) && (n_v > 1.0 && n_v < 50.0);
        }
        else
        { // gauss
            sig_report = sigma->getVal();
            sig_err = sigma->getError();
            sane = sane && (sig_report >= 0.003 && sig_report <= 0.020);
        }

        bool better = (best.order < 0) || (AICc < best.aicc);
        if (better)
        {
            if (best.res)
                delete best.res;
            if (best.model)
                delete best.model;
            if (best.frame)
                delete best.frame;

            best.order = k;
            best.aicc = AICc;
            best.chi2ndf = chi2;
            best.nll = nll;
            best.mu = mu_v;
            best.muerr = mu_e;
            best.sig = sig_report;
            best.sigerr = sig_err;
            best.Ns = Ns_v;
            best.Nserr = Ns_e;
            best.Nb = Nb_v;
            best.Nberr = Nb_e;
            best.status = status;
            best.smode = smode_in;

            if (smode_in == "doubleG")
            {
                best.sig1 = s1_v;
                best.sig1err = s1_e;
                best.sig2 = s2_v;
                best.sig2err = s2_e;
                best.frac = f_v;
                best.fracerr = f_e;
            }
            else if (smode_in == "cb")
            {
                best.alpha = a_v;
                best.alphaerr = a_e;
                best.n = n_v;
                best.nerr = n_e;
            }

            best.res = (RooFitResult *)res.release();
            res.reset();
            best.model = (RooAddPdf *)model.cloneTree();
            best.frame = (RooPlot *)fr.release();
        }
    }

    std::cout << std::fixed << std::setprecision(6);
    if (best.smode == "doubleG")
    {
        std::cout << "[mgg-sb-fit] range=[" << mmin << "," << mmax << "]  "
                  << "order=" << best.order
                  << "  mu=" << best.mu << "±" << best.muerr
                  << "  sigma_eff=" << best.sig << "±" << best.sigerr
                  << "  Ns=" << best.Ns << "±" << best.Nserr
                  << "  Nb=" << best.Nb << "±" << best.Nberr
                  << "  chi2/ndf=" << best.chi2ndf
                  << "  AICc=" << best.aicc
                  << "  status=" << best.status
                  << "  N(bins used)=" << N_for_AICc
                  << "\n";
    }
    else if (best.smode == "cb")
    {
        std::cout << "[mgg-sb-fit] range=[" << mmin << "," << mmax << "]  "
                  << "order=" << best.order
                  << "  mu=" << best.mu << "±" << best.muerr
                  << "  sigma(CB)=" << best.sig << "±" << best.sigerr
                  << "  alpha=" << best.alpha << "±" << best.alphaerr
                  << "  n=" << best.n << "±" << best.nerr
                  << "  Ns=" << best.Ns << "±" << best.Nserr
                  << "  Nb=" << best.Nb << "±" << best.Nberr
                  << "  chi2/ndf=" << best.chi2ndf
                  << "  AICc=" << best.aicc
                  << "  status=" << best.status
                  << "  N(bins used)=" << N_for_AICc
                  << "\n";
    }
    else
    {
        std::cout << "[mgg-sb-fit] range=[" << mmin << "," << mmax << "]  "
                  << "order=" << best.order
                  << "  mu=" << best.mu << "±" << best.muerr
                  << "  sigma=" << best.sig << "±" << best.sigerr
                  << "  Ns=" << best.Ns << "±" << best.Nserr
                  << "  Nb=" << best.Nb << "±" << best.Nberr
                  << "  chi2/ndf=" << best.chi2ndf
                  << "  AICc=" << best.aicc
                  << "  status=" << best.status
                  << "  N(bins used)=" << N_for_AICc
                  << "\n";
    }

    if (!saveQA)
    {
        delete best.res;
        delete best.model;
        delete best.frame;
        return;
    }

    TCanvas c("c", "mgg S+B (RooFit)", 900, 700);
    c.SetMargin(0.12, 0.04, 0.12, 0.06);
    RooPlot *fr = m.frame(Bins(nbinsInRange), Title("m_{#gamma#gamma} S+B fit"));
    data.plotOn(fr, DataError(RooAbsData::SumW2));
    best.model->plotOn(fr, Components("bkg"), LineStyle(kDashed), LineColor(kBlue + 1));
    best.model->plotOn(fr, Components("sig"), LineStyle(kDotted), LineColor(kRed + 1));
    best.model->plotOn(fr);
    fr->GetXaxis()->SetTitle("m_{#gamma#gamma}  [GeV]");
    fr->GetYaxis()->SetTitle("Counts per bin");
    fr->Draw();
    TLegend leg(0.58, 0.72, 0.88, 0.88);
    leg.AddEntry((TObject *)nullptr, Form("order=%d", best.order), "");
    leg.AddEntry((TObject *)nullptr, Form("#mu=%.6f #pm %.6f", best.mu, best.muerr), "");
    if (best.smode == "doubleG")
    {
        leg.AddEntry((TObject *)nullptr, Form("#sigma_{eff}=%.6f", best.sig), "");
    }
    else if (best.smode == "cb")
    {
        leg.AddEntry((TObject *)nullptr, Form("#sigma_{CB}=%.6f, #alpha=%.2f, n=%.1f", best.sig, best.alpha, best.n), "");
    }
    else
    {
        leg.AddEntry((TObject *)nullptr, Form("#sigma=%.6f #pm %.6f", best.sig, best.sigerr), "");
    }
    leg.AddEntry((TObject *)nullptr, Form("#chi^{2}/ndf=%.3f", best.chi2ndf), "");
    leg.Draw();
    // --- Shaded sigma bands on the plot (auto-updated from fit)
    double mu = best.mu;
    double sig = best.sig; // for doubleG we store sigma_eff here; for CB this is CB sigma
    double ylo = 0.0;
    double yhi = 1.10 * std::max(1.0, hRange->GetMaximum()); // simple, robust y-extent

    // 3σ (widest) at the back, then 2.5σ, then 2σ on top
    TBox *b3 = new TBox(mu - 3.0 * sig, ylo, mu + 3.0 * sig, yhi);
    TBox *b25 = new TBox(mu - 2.5 * sig, ylo, mu + 2.5 * sig, yhi);
    TBox *b2 = new TBox(mu - 2.0 * sig, ylo, mu + 2.0 * sig, yhi);

    // Soft translucent fills so data/model remain visible
    b3->SetFillColorAlpha(kAzure - 9, 0.12);
    b3->SetLineColor(0);
    b25->SetFillColorAlpha(kTeal - 7, 0.18);
    b25->SetLineColor(0);
    b2->SetFillColorAlpha(kOrange - 3, 0.25);
    b2->SetLineColor(0);

    // Draw after the frame so they overlay behind the legend
    b3->Draw("same");
    b25->Draw("same");
    b2->Draw("same");

    // Vertical line at μ
    TLine *lmu = new TLine(mu, ylo, mu, yhi);
    lmu->SetLineColor(kGray + 2);
    lmu->SetLineStyle(2);
    lmu->SetLineWidth(2);
    lmu->Draw("same");

    // Extend the existing legend with band labels
    leg.AddEntry(b2, "#mu #pm 2#sigma", "f");
    leg.AddEntry(b25, "#mu #pm 2.5#sigma", "f");
    leg.AddEntry(b3, "#mu #pm 3#sigma", "f");
    leg.AddEntry(lmu, "#mu", "l");

    c.SaveAs("mgg_roofit_sb.png");

    // ---- Build per-bin background expectation and weights (on h_mgg_used)
    TH1D *hUsed = (TH1D *)best.frame ? (TH1D *)gDirectory->Get("h_mgg_used") : nullptr; // we stored it later; use our local pointer instead:
    hUsed = (TH1D *)gDirectory->Get("h_mgg_used");                                      // in case ROOT registered it
    if (!hUsed)
        hUsed = hRange.get(); // fallback to local

    const int nbw = hRange->GetNbinsX();
    TH1D *hBkgExp = (TH1D *)hRange->Clone("hBkg_expect");
    hBkgExp->Reset("ICESM");
    TH1D *hSigSub = (TH1D *)hRange->Clone("hSig_subtracted");
    hSigSub->Reset("ICESM");
    TH1D *hWsig = (TH1D *)hRange->Clone("hW_sig");
    hWsig->Reset("ICESM");
    hBkgExp->SetTitle("Fitted background expectation per bin");
    hSigSub->SetTitle("Data - fitted background per bin");
    hWsig->SetTitle("Per-bin signal weight w_sig = (D-B)/D");

    // We need the background PDF from best.model
    RooAbsPdf *model_final = best.model; // RooAddPdf(sig + bkg)
    RooAbsPdf *bkg_pdf = (RooAbsPdf *)model_final->getComponents()->find("bkg");
    RooRealVar *Nb_var = (RooRealVar *)model_final->getVariables()->find("Nb");
    double Nb_val = Nb_var ? Nb_var->getVal() : best.Nb;

    for (int i = 1; i <= nbw; ++i)
    {
        double xlo = hRange->GetXaxis()->GetBinLowEdge(i);
        double xhi = hRange->GetXaxis()->GetBinUpEdge(i);
        // define a named range for this bin
        TString rname;
        rname.Form("bin_%d", i);
        m.setRange(rname, xlo, xhi);
        // integral of background PDF over this bin range (normalized over [mmin,mmax])
        std::unique_ptr<RooAbsReal> I(bkg_pdf->createIntegral(m, NormSet(m), Range(rname)));
        double prob = I->getVal(); // fraction in this bin
        double Bi = Nb_val * prob; // expected background counts in bin
        double Di = hRange->GetBinContent(i);
        double Si = Di - Bi;
        double wi = (Di > 0.0 ? Si / Di : 0.0);

        hBkgExp->SetBinContent(i, Bi);
        hSigSub->SetBinContent(i, Si);
        hWsig->SetBinContent(i, wi);
    }
    // (optional) copy bin errors from data into S for quick visual χ check
    for (int i = 1; i <= nbw; ++i)
        hSigSub->SetBinError(i, hRange->GetBinError(i));

    std::unique_ptr<TFile> fout(TFile::Open("fit_mgg_roofit_out.root", "RECREATE"));
    if (fout && !fout->IsZombie())
    {
        fout->WriteTObject(best.res, "fitResult", "Overwrite");
        fout->WriteTObject(fr, "frame", "Overwrite");
        fout->WriteTObject(best.model, "model", "Overwrite");
        fout->WriteTObject(hRange.get(), "h_mgg_used", "Overwrite");
        fout->WriteTObject(hBkgExp, "hBkg_expect", "Overwrite");
        fout->WriteTObject(hSigSub, "hSig_sub", "Overwrite");
        fout->WriteTObject(hWsig, "hW_sig", "Overwrite");
        fout->Write();
        fout->Close();
    }

    std::ofstream js("mgg_roofit_sb.json");
    if (js)
    {
        js << std::setprecision(8);
        js << "{\n";
        js << "  \"input_file\": \"" << infile << "\",\n";
        js << "  \"hist_path\": \"" << histpath << "\",\n";
        js << "  \"range\": [" << mmin << ", " << mmax << "],\n";
        js << "  \"order\": " << best.order << ",\n";
        js << "  \"signal_mode\": \"" << best.smode << "\",\n";
        js << "  \"mu\": " << best.mu << ", \"mu_err\": " << best.muerr << ",\n";
        if (best.smode == "doubleG")
        {
            js << "  \"sigma_eff\": " << best.sig << ", \"sigma_eff_err\": " << best.sigerr << ",\n";
            js << "  \"Ns\": " << best.Ns << ", \"Ns_err\": " << best.Nserr << ",\n";
        }
        else if (best.smode == "cb")
        {
            js << "  \"sigma_cb\": " << best.sig << ", \"sigma_cb_err\": " << best.sigerr << ",\n";
            js << "  \"alpha\": " << best.alpha << ", \"alpha_err\": " << best.alphaerr << ",\n";
            js << "  \"n\": " << best.n << ", \"n_err\": " << best.nerr << ",\n";
            js << "  \"Ns\": " << best.Ns << ", \"Ns_err\": " << best.Nserr << ",\n";
        }
        else
        {
            js << "  \"sigma\": " << best.sig << ", \"sigma_err\": " << best.sigerr << ",\n";
            js << "  \"Ns\": " << best.Ns << ", \"Ns_err\": " << best.Nserr << ",\n";
        }
        js << "  \"Nb\": " << best.Nb << ", \"Nb_err\": " << best.Nberr << ",\n";
        js << "  \"AICc\": " << best.aicc << ",\n";
        js << "  \"chi2_ndf\": " << best.chi2ndf << ",\n";
        js << "  \"status\": " << best.status << ",\n";
        js << "  \"bins_used\": " << N_for_AICc << "\n";
        js << "  \"rebin\": " << rebin << ",\n";
        js << "}\n";
        js.close();
    }

    delete best.res;
    delete best.model;
    delete best.frame;
}

// -------------------------------------------------------------
// YAML DRIVER
// -------------------------------------------------------------
int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] 
                  << " config.yaml\n";
        return 1;
    }

    YAML::Node cfg;
    try {
        cfg = YAML::LoadFile(argv[1]);
    } catch (const YAML::Exception &e) {
        std::cerr << "Failed to load YAML: " << e.what() << "\n";
        return 1;
    }

    // Required parameters
    std::string infile  = cfg["infile"].as<std::string>();
    std::string histpath = cfg["histpath"].as<std::string>();

    // Fit range (maps naturally to mmin, mmax)
    double mmin = cfg["fit_range"]["min"].as<double>();
    double mmax = cfg["fit_range"]["max"].as<double>();

    // Optional / defaulted
    int bkg_order = cfg["background_order"].as<int>();   // use for ord_min/ord_max
    bool saveQA = cfg["save_QA"].as<bool>();

    // Additional options
    std::string signal_mode = "gauss";
    if (cfg["signal_mode"])
        signal_mode = cfg["signal_mode"].as<std::string>();

    int rebin = 1;
    if (cfg["rebin"])
        rebin = cfg["rebin"].as<int>();

    // ---------------------------------------------------------
    // Call your analysis function
    // ---------------------------------------------------------
    fit_mgg_roofit_sb(
        infile.c_str(),
        histpath.c_str(),
        mmin,
        mmax,
        saveQA,
        bkg_order,  // ord_min
        bkg_order,  // ord_max (same unless you want a range)
        signal_mode.c_str(),
        rebin
    );

    return 0;
}