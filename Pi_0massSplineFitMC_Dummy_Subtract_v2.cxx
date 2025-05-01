/*************************************************************
  Pi_0massSplineFitMC_Dummy_Subtract_v2.cpp
  – arithmetic in raw counts, global 1/realCharge scale at the end
  – dummy already normalised to counts / µC, then converted back to raw
  – builds Toy-MC-weighted π0-mass and Q² distributions (6-pad canvas)

  Usage:
    ./Pi_0massSplineFitMC_Dummy_Subtract_v2 \
        real.root out.pdf dummy.root realCharge_uC dummyCharge_uC
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

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

// ───────────────────────────────── smoothing helpers
std::vector<double> normalizePeaks(const std::vector<double>& y,int win,double frac){
    int n=y.size(); std::vector<double> out=y;
    double ymax=*std::max_element(y.begin(),y.end()),thr=frac*ymax;
    std::vector<int> cand;
    for(int i=1;i<n-1;++i)
        if(y[i]>y[i-1]&&y[i]>y[i+1]&&y[i]>thr) cand.push_back(i);
    std::vector<int> peaks;
    for(int p:cand){
        if(peaks.empty()) peaks.push_back(p);
        else{
            int last=peaks.back();
            if(p-last<win){ if(y[p]>y[last]) peaks.back()=p; }
            else peaks.push_back(p);
        }}
    if(peaks.empty()) return out;
    double avg=0; for(int p:peaks) avg+=y[p]; avg/=peaks.size();
    int h=win/2;
    for(int p:peaks){
        double scl=avg/y[p];
        for(int i=std::max(0,p-h);i<=std::min(n-1,p+h);++i){
            double w=1.0-std::fabs(i-p)/double(win+1);
            out[i]=y[i]*(1-w)+(y[i]*scl)*w;
        }}
    return out;
}
std::vector<double> normalizeTroughs(const std::vector<double>& y,int win,double tol){
    int n=y.size(); std::vector<double> out=y;
    double ymin=*std::min_element(y.begin(),y.end()),
           ymax=*std::max_element(y.begin(),y.end());
    std::vector<int> cand;
    for(int i=1;i<n-1;++i)
        if(y[i]<y[i-1]&&y[i]<y[i+1]&&y[i]<=ymin+tol*(ymax-ymin)) cand.push_back(i);
    std::vector<int> tr;
    for(int t:cand){
        if(tr.empty()) tr.push_back(t);
        else{
            int last=tr.back();
            if(t-last<win){ if(y[t]<y[last]) tr.back()=t; }
            else tr.push_back(t);
        }}
    if(tr.empty()) return out;
    double avg=0; for(int t:tr) avg+=y[t]; avg/=tr.size();
    int h=win/2;
    for(int t:tr){
        double scl=avg/y[t];
        for(int i=std::max(0,t-h);i<=std::min(n-1,t+h);++i){
            double w=1.0-std::fabs(i-t)/double(win+1);
            out[i]=y[i]*(1-w)+(y[i]*scl)*w;
        }}
    return out;
}
std::vector<double> regionSmooth(const std::vector<double>& d,int swP,int swT){
    int n=d.size(); std::vector<double> s(n);
    double ymax=*std::max_element(d.begin(),d.end()),
           ymin=*std::min_element(d.begin(),d.end()),
           mid =0.5*(ymax+ymin);
    for(int i=0;i<n;++i){
        int win=(d[i]>mid)?swP:swT; double sum=0; int cnt=0;
        for(int j=i-win/2;j<=i+win/2;++j)
            if(j>=0&&j<n){ sum+=d[j]; ++cnt; }
        s[i]=(cnt?sum/cnt:d[i]);
    }
    return s;
}

// ───────────────────────────────── simple cuts
inline bool passHMSCuts(double edt,double dp,double et,double npe,
                        double th ,double ph){
    return (edt<0.1 && std::fabs(dp)<=8.5 && et>0.6 && npe>1.0 &&
            std::fabs(th)<=0.09 && std::fabs(ph)<=0.09);
}
inline bool isGoodCluster(double e,double t,double x,double y){
    return (e>=0.6 && t>=149&&t<=151 &&
            x>-29.16&&x<29.16 && y>-35.64&&y<35.64);
}

// ───────────────────────────────── quick Q² helper
// --- constants for current kinematics --------------------
const double theta0_deg = 16.44;                  // HMS central angle (deg)
const double theta0_rad = theta0_deg * M_PI/180.; // radians
// ---------------------------------------------------------
inline double compQ2(double E0,double Ep,double th,double ph)
{
    // track slopes are small: th≈δy/p , ph≈δx/p
    // total polar deflection: add HMS setting angle to vertical slope
    double th_tot = theta0_rad + th;
    double cosT   = std::cos(th_tot) * std::cos(ph);   // to O(θ²)
    return 2.0*E0*Ep*(1.0 - cosT);
}

// ──────────────────────────────────────────── main
int main(int argc,char* argv[])
{
    if(argc<6){
        std::cerr<<"Usage: "<<argv[0]
                 <<" real.root out.pdf dummy.root realCharge_uC dummyCharge_uC\n";
        return 1;
    }
    std::string realF  = argv[1];
    std::string outPDF = argv[2];
    std::string dumF   = argv[3];
    double      realQ  = std::stod(argv[4]);   // µC
    double      dumQ   = std::stod(argv[5]);   // µC

    //----------------------------------------------------------------
    // constants (current kinematic = E₀ 10.538 GeV, ep₀ 4.637 GeV)
    //----------------------------------------------------------------
    const double e0_nom   = 10.538;   // GeV (beam)
    const double ep0_nom  =  4.637;   // GeV (central scattered electron)
    const double bgLo=113, bgHi=142.5;
    const double sigLo=141.789, sigHi=171.289;
    const double shift=28.05;
    const int    nBins=650;

    //────────────────────────────────── 1) dummy histogram (counts / µC)
    TH1F *hDumNorm=nullptr;
    {
        TFile fd(dumF.c_str(),"READ");
        if(!fd.IsZombie()){
            TTree *t = dynamic_cast<TTree*>( fd.Get("T") );
            if(t){
                // raw variables
                double edt=0,dp=0,et=0,npe=0,th=0,ph=0,gy=0,ct=0;
                t->SetBranchStatus("*",0);
                t->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw",1);
                t->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw",&edt);
                t->SetBranchStatus("H.gtr.dp",1);
                t->SetBranchAddress("H.gtr.dp",&dp);
                t->SetBranchStatus("H.cal.etotnorm",1);
                t->SetBranchAddress("H.cal.etotnorm",&et);
                t->SetBranchStatus("H.cer.npeSum",1);
                t->SetBranchAddress("H.cer.npeSum",&npe);
                t->SetBranchStatus("H.gtr.th",1);
                t->SetBranchAddress("H.gtr.th",&th);
                t->SetBranchStatus("H.gtr.ph",1);
                t->SetBranchAddress("H.gtr.ph",&ph);
                t->SetBranchStatus("H.gtr.y",1);
                t->SetBranchAddress("H.gtr.y",&gy);
                t->SetBranchStatus("NPS.cal.clusT",1);
                t->SetBranchAddress("NPS.cal.clusT",&ct);

                TH1F *hUp=new TH1F("dum_up","",nBins,sigLo,sigHi);
                TH1F *hDn=new TH1F("dum_dn","",nBins,sigLo,sigHi);
                hUp->SetDirectory(nullptr); hDn->SetDirectory(nullptr);

                const Long64_t N=t->GetEntries();
                for(Long64_t i=0;i<N;++i){
                    t->GetEntry(i);
                    if(!passHMSCuts(edt,dp,et,npe,th,ph)) continue;
                    if(ct<sigLo||ct>sigHi) continue;
                    (gy>0? hDn:hUp)->Fill(ct);
                }
                hUp->Scale( 1.0/(dumQ*8.467) );
                hDn->Scale( 1.0/(dumQ*4.256) );
                hDumNorm = static_cast<TH1F*>( hUp->Clone("hDumNorm") );
                hDumNorm->Add(hDn);
                hDumNorm->SetDirectory(nullptr);
            }
        }
    }

    //────────────────────────────────── 2) real file (raw)
    TFile fr(realF.c_str(),"READ");
    if(fr.IsZombie()){ std::cerr<<"Cannot open "<<realF<<"\n"; return 1; }
    TTree *tr = dynamic_cast<TTree*>( fr.Get("T") );
    if(!tr){ std::cerr<<"No T tree in "<<realF<<"\n"; return 1; }

    double edt=0,dp=0,et=0,npe=0,th=0,ph=0;
    tr->SetBranchStatus("*",0);
    tr->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw",1);
    tr->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw",&edt);
    tr->SetBranchStatus("H.gtr.dp",1);
    tr->SetBranchAddress("H.gtr.dp",&dp);
    tr->SetBranchStatus("H.cal.etotnorm",1);
    tr->SetBranchAddress("H.cal.etotnorm",&et);
    tr->SetBranchStatus("H.cer.npeSum",1);
    tr->SetBranchAddress("H.cer.npeSum",&npe);
    tr->SetBranchStatus("H.gtr.th",1);
    tr->SetBranchAddress("H.gtr.th",&th);
    tr->SetBranchStatus("H.gtr.ph",1);
    tr->SetBranchAddress("H.gtr.ph",&ph);

    double ncl=0; static const int MAX=10000;
    double cE[MAX],cT[MAX],cX[MAX],cY[MAX];
    tr->SetBranchStatus("NPS.cal.nclust",1);
    tr->SetBranchAddress("NPS.cal.nclust",&ncl);
    tr->SetBranchStatus("NPS.cal.clusE",1);
    tr->SetBranchAddress("NPS.cal.clusE",cE);
    tr->SetBranchStatus("NPS.cal.clusT",1);
    tr->SetBranchAddress("NPS.cal.clusT",cT);
    tr->SetBranchStatus("NPS.cal.clusX",1);
    tr->SetBranchAddress("NPS.cal.clusX",cX);
    tr->SetBranchStatus("NPS.cal.clusY",1);
    tr->SetBranchAddress("NPS.cal.clusY",cY);

    TH1F hBG ("bg","",nBins,bgLo,bgHi);  hBG .SetDirectory(nullptr);
    TH1F hSig("sig","",nBins,sigLo,sigHi);hSig.SetDirectory(nullptr);

    Long64_t Nr = tr->GetEntries();
    for(Long64_t i=0;i<Nr;++i){
        tr->GetEntry(i);
        if(!passHMSCuts(edt,dp,et,npe,th,ph)) continue;
        if(ncl<1) continue;
        double t0=cT[0];
        if(t0>=bgLo&&t0<=bgHi ) hBG .Fill(t0);
        if(t0>=sigLo&&t0<=sigHi) hSig.Fill(t0);
    }

    //────────────────────────────────── 3) spline BG and subtract
    std::vector<double> vx(nBins),vy(nBins);
    for(int b=1;b<=nBins;++b){
        vx[b-1]=hBG.GetBinCenter(b);
        vy[b-1]=hBG.GetBinContent(b);
    }
    TGraph gBG(nBins,&vx[0],&vy[0]); gBG.Sort();
    TSpline3 spl("spl",&gBG);
    std::vector<double> yS(nBins);
    for(int i=0;i<nBins;++i) yS[i]=spl.Eval(vx[i]);
    auto yP=normalizePeaks  (yS,10,0.75);
    auto yT=normalizeTroughs(yP,15,0.05);
    auto yF=regionSmooth    (yT,2,5);

    std::vector<double> vxs(nBins),vys(nBins);
    for(int i=0;i<nBins;++i){ vxs[i]=vx[i]+shift; vys[i]=yF[i]; }
    TGraph gShift(nBins,&vxs[0],&vys[0]);

    TH1F hSub = *static_cast<TH1F*>(hSig.Clone("hSub"));
    hSub.SetDirectory(nullptr);
    for(int b=1;b<=nBins;++b){
        double c=hSub.GetBinCenter(b);
        hSub.SetBinContent(b, hSig.GetBinContent(b)-gShift.Eval(c));
    }

    if(hDumNorm){
        TH1F hD=*hDumNorm; hD.Scale(realQ); // convert back to raw
        hSub.Add(&hD,-1.0);
    }
    std::cout<<"Integral hSub (raw) = "<<hSub.Integral()<<"\n";

    //────────────────────────────────── 4) arrays for Toy MC
    std::vector<double> dataVal(nBins),bgVal(nBins),bgErr(nBins);
    for(int b=1;b<=nBins;++b){
        dataVal[b-1]=hSub.GetBinContent(b);
        double est=gShift.Eval(hSub.GetBinCenter(b));
        bgVal[b-1]=est;
        bgErr[b-1]=std::sqrt(std::max(est,0.0));
    }

    //────────────────────────────────── 5) build cluster list + event list
    struct Cl  {int tbin; double m;};
    struct Ev  {int tbin; double q2;};
    std::vector<Cl> cls; cls.reserve(1000000);
    std::vector<Ev> evtList; evtList.reserve(500000);

    const double DNPS=407.0; // distance for cluster-pair angle
    for(Long64_t ie=0;ie<Nr;++ie){
        tr->GetEntry(ie);
        if(!passHMSCuts(edt,dp,et,npe,th,ph)) continue;

        // event-by-event Q² using measured slopes
        double ep_evt = ep0_nom*(1.0+dp/100.0);
        double q2_evt = compQ2(e0_nom,ep_evt,th,ph);

        if(ncl<1) continue;
        double tfirst = cT[0];
        if(tfirst<sigLo||tfirst>sigHi) continue;
        int tbin_evt = hSub.FindBin(tfirst)-1;
        if(tbin_evt<0||tbin_evt>=nBins) continue;
        evtList.push_back({tbin_evt,q2_evt});

        int N=(int)ncl; if(N<2) continue;
        std::vector<int> keep; keep.reserve(N);
        for(int c=0;c<N;++c){
            if(cT[c]<sigLo||cT[c]>sigHi) continue;
            if(!isGoodCluster(cE[c],cT[c],cX[c],cY[c])) continue;
            keep.push_back(c);
        }
        if(keep.size()<2) continue;

        for(size_t a=0;a<keep.size();++a)
            for(size_t b=a+1;b<keep.size();++b){
                int i1=keep[a], i2=keep[b];
                double dx=cX[i1]-cX[i2], dy=cY[i1]-cY[i2];
                double theta=std::sqrt(dx*dx+dy*dy)/DNPS;
                double s2   = std::sin(0.5*theta);
                double m    = std::sqrt(4.0*cE[i1]*cE[i2]*s2*s2);
                double tM   = 0.5*(cT[i1]+cT[i2]);
                int tb = hSub.FindBin(tM)-1;
                if(tb<0||tb>=nBins) continue;
                cls.push_back({tb,m});
            }
    }
    std::cout<<"Cluster list size = "<<cls.size()<<"\n";
    std::cout<<"Event   list size = "<<evtList.size()<<"\n";

    //────────────────────────────────── 6) Toy MC
    const int Ntoys=200, nMB=200; const double mLo=0, mHi=0.3;
    const int nQ2B=120; const double q2Lo=0, q2Hi=12;        // GeV²

    TRandom3 rng(0);
    std::vector<TH1F*> toyPi, toyQ2;
    toyPi.reserve(Ntoys); toyQ2.reserve(Ntoys);

    for(int it=0;it<Ntoys;++it){
        std::vector<double> w(nBins,0.0);
        for(int i=0;i<nBins;++i){
            double dev=rng.Gaus(0.0,bgErr[i]);
            double bt = bgVal[i]+dev; if(bt<0) bt=0;
            double d  = dataVal[i];   if(d<0) d=0;
            double f  = (d>1e-6)?((d-bt)/d):0;
            if(f> 1) f= 1;
            if(f<-1) f=-1;
            w[i]=f;
        }
        TH1F *hM=new TH1F(Form("toyM%d",it),"",nMB ,mLo ,mHi );
        TH1F *hQ=new TH1F(Form("toyQ%d",it),"",nQ2B,q2Lo,q2Hi);
        for(const auto& c:cls)      hM->Fill(c.m , w[c.tbin]);
        for(const auto& ev:evtList) hQ->Fill(ev.q2, w[ev.tbin]);
        toyPi.push_back(hM); toyQ2.push_back(hQ);
    }

    // mean & rms for π0 mass
    TH1F* hMean=(TH1F*)toyPi[0]->Clone("hMean"); hMean->Reset();
    TH1F* hVar =(TH1F*)toyPi[0]->Clone("hVar" ); hVar ->Reset();
    for(auto h:toyPi) hMean->Add(h);
    hMean->Scale(1.0/Ntoys);
    for(auto h:toyPi)
        for(int b=1;b<=nMB;++b){
            double diff=h->GetBinContent(b)-hMean->GetBinContent(b);
            hVar->AddBinContent(b,diff*diff);}
    hVar->Scale(1.0/Ntoys);
    for(int b=1;b<=nMB;++b)
        hVar->SetBinContent(b,std::sqrt(hVar->GetBinContent(b)));

    // mean & rms for Q²
    TH1F* hQmean=(TH1F*)toyQ2[0]->Clone("hQmean"); hQmean->Reset();
    TH1F* hQvar =(TH1F*)toyQ2[0]->Clone("hQvar" ); hQvar ->Reset();
    for(auto h:toyQ2) hQmean->Add(h);
    hQmean->Scale(1.0/Ntoys);
    for(auto h:toyQ2)
        for(int b=1;b<=nQ2B;++b){
            double diff=h->GetBinContent(b)-hQmean->GetBinContent(b);
            hQvar->AddBinContent(b,diff*diff);}
    hQvar->Scale(1.0/Ntoys);
    for(int b=1;b<=nQ2B;++b)
        hQvar->SetBinContent(b,std::sqrt(hQvar->GetBinContent(b)));

    // global scale → counts/µC
    double sf=1.0/realQ;
    hMean ->Scale(sf);  hVar ->Scale(sf);
    hQmean->Scale(sf);  hQvar->Scale(sf);

    // error-band helpers
    auto makeBand=[&](TH1F* hC,TH1F* hE)->TGraph*{
        int nb=hC->GetNbinsX();
        std::vector<double> px,py; px.reserve(2*nb); py.reserve(2*nb);
        for(int b=1;b<=nb;++b){
            px.push_back(hC->GetBinCenter(b));
            py.push_back(hC->GetBinContent(b));
        }
        for(int b=nb;b>=1;--b){
            px.push_back(hC->GetBinCenter(b));
            py.push_back(hC->GetBinContent(b)+hE->GetBinContent(b));
        }
        auto *g=new TGraph(px.size(),px.data(),py.data());
        g->SetFillColorAlpha(kBlue,0.35);
        return g;
    };
    TGraph *gBandM = makeBand(hMean ,hVar );
    TGraph *gBandQ = makeBand(hQmean,hQvar);

    //────────────────────────────────── 7) Canvas (6 pads)
    TCanvas c("c","",1200,3600); c.Divide(1,6);

    // pad 1 : BG points + spline
    c.cd(1); 
    // Draw original BG points clearly
    gBG.SetMarkerStyle(20);
    gBG.SetMarkerSize(0.8);
    gBG.Draw("AP");

    // Explicitly construct and plot the final smoothed BG curve (yF)
    TGraph* grBGfinal = new TGraph(nBins, &vx[0], &yF[0]);
    grBGfinal->SetLineColor(kBlue);
    grBGfinal->SetLineWidth(2);
    grBGfinal->Draw("L SAME");
    
    // Now explicitly overlay the spline curve:
    spl.SetLineColor(kRed);
    spl.SetLineWidth(2);
    spl.Draw("C SAME");  // clearly overlay spline fit (red curve)

    c.cd(2); hSig.Draw("HIST");
             hSub.SetLineColor(kGreen+2); hSub.Draw("HIST SAME");

    c.cd(3); hMean->Draw("HIST");

    c.cd(4); hVar ->SetLineColor(kRed); hVar->Draw("HIST");

    // pad 5 : π0 mass + band
    c.cd(5);
    double ymaxM = 1.2*(hMean->GetMaximum()+hVar->GetMaximum());
    if(ymaxM<1e-6) ymaxM=1e-6;
    TH2F frame5("f5",";M_{#gamma#gamma} (GeV);Counts/µC",
                10,0,0.30,10,0,ymaxM);
    frame5.Draw("AXIS");
    gBandM->Draw("F SAME");
    hMean ->Draw("HIST SAME");

    // pad 6 : Q² + band
    c.cd(6);
    double ymaxQ = 1.2*(hQmean->GetMaximum()+hQvar->GetMaximum());
    if(ymaxQ<1e-8) ymaxQ=1e-8;
    TH2F frame6("f6",";Q^{2} (GeV^{2});Counts/µC",
                10,q2Lo,q2Hi,10,0,ymaxQ);
    frame6.Draw("AXIS");
    gBandQ->Draw("F SAME");
    hQmean->Draw("HIST SAME");

    c.Print(outPDF.c_str());
    std::cout<<"All done. Canvas saved to "<<outPDF<<"\n";

    gSystem->Exit(0);
    return 0;
}


