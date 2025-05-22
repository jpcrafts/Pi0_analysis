    /*************************************************************
     Pi_0massSplineFitMC_Dummy_Subtract_v3.cpp
    – arithmetic in raw counts, global 1/realCharge scale at the end
    – dummy already normalised to counts / µC, then converted back to raw
    – builds Toy-MC-weighted π0-mass and Q² distributions (6-pad canvas)

    Usage:
        ./Pi_0massSplineFitMC_Dummy_Subtract_v3 \
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
    #include <TLegend.h>
    #include <TPaveText.h>
    #include <TF1.h>

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
    const double theta0_deg = 16.48;                  // HMS central angle (deg)
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
        const double e0_nom   = 10.54350201;   // GeV (beam)
        const double ep0_nom  =  5.878;   // GeV (central scattered electron)
        const double bgLo=113, bgHi=142.5;
        const double sigLo=141.789, sigHi=171.289;
        const double shift=28.05;
        const int    nBins=650;

        //────────────────────────────────── 1) dummy histogram (counts / µC)
        TH1F *hDumNorm=nullptr;
        // Declare histograms for Pad 7 visualization (matching hSub binning)
        TH1F hSub_before("hSub_before", "", nBins, sigLo, sigHi);
        TH1F hSub_after("hSub_after", "", nBins, sigLo, sigHi);
        TH1F hD("hD", "", nBins, sigLo, sigHi);
        // Histogram for missing mass (200 bins from 0 to 5 GeV/c^2 for debug)
        TH1F hMissMass("hMissMass", "Missing Mass;M_{miss} [GeV/c^{2}];Counts", 200, 0, 5);
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

        // New: Histogram for all clusT[0] values (full range)
        const double t0Min = 0.0; // Adjust as needed for full clusT[0] range
        const double t0Max = 200.0; // Adjust as needed for full clusT[0] range
        // Set bin width so that the background region (bgHi-bgLo) has 650 bins
        const int nRegionBins = 650;
        const double regionWidth = bgHi - bgLo; // or use sigHi-sigLo if preferred
        const double binWidth = regionWidth / nRegionBins;
        // Extend hAll to cover the full t0 range with this bin width
        // const double t0Min = 0.0; // Removed redeclaration
        // const double t0Max = 200.0; // Removed redeclaration
        const int nAllBins = static_cast<int>((t0Max - t0Min) / binWidth + 0.5);

        // First pass: Fill hAll with all clusT[0] values (no HMS cuts, only ncl>0)
        TH1F hAll("hAll", "All clusT[0] (charge normalized)", nAllBins, t0Min, t0Max);
        hAll.SetDirectory(nullptr);

        Long64_t Nr = tr->GetEntries();
        for(Long64_t i=0;i<Nr;++i){
            tr->GetEntry(i);
            if(ncl<1) continue;
            double t0 = cT[0];
            hAll.Fill(t0); // Fill all t0 values, no cuts
        }

        // Normalize hAll by realQ to get charge-normalized raw yield
        hAll.Scale(1.0 / realQ);

        // Output hAll as the charge-normalized raw yield histogram
        TFile outHist("charge_normalized_raw_yield.root", "RECREATE");
        hAll.Write();
        outHist.Close();

        // Second pass: Apply HMS cuts and process for background/signal as before
        // Create hBG and hSig as 650-bin histograms over their respective regions
        TH1F hBG("bg", "", nRegionBins, bgLo, bgHi);  hBG.SetDirectory(nullptr);
        TH1F hSig("sig", "", nRegionBins, sigLo, sigHi); hSig.SetDirectory(nullptr);

        for(Long64_t i=0;i<Nr;++i){
            tr->GetEntry(i);
            if(!passHMSCuts(edt,dp,et,npe,th,ph)) continue;
            if(ncl<1) continue;
            double t0 = cT[0];
            // Fill hBG and hSig by mapping t0 to the correct bin in each region
            if(t0>=bgLo&&t0<=bgHi ) hBG.Fill(t0);
            if(t0>=sigLo&&t0<=sigHi) hSig.Fill(t0);
        }
        hBG.SetBins(nBins, bgLo, bgHi);  // Reuse existing hBG, reset bins
        hSig.SetBins(nBins, sigLo, sigHi); hSig.SetDirectory(nullptr); // Reuse existing hSig, reset bins

        for(Long64_t i=0;i<Nr;++i){
            tr->GetEntry(i);
            if(!passHMSCuts(edt,dp,et,npe,th,ph)) continue;
            if(ncl<1) continue;
            double t0 = cT[0];
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
            hSub_before = hSub; // Save hSub before dummy subtraction
            hD = *hDumNorm; hD.Scale(realQ); // convert back to raw
            hSub.Add(&hD,-1.0);
            hSub_after = hSub; // Save hSub after dummy subtraction
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
        
            // --- Missing mass calculation for events with at least two good clusters ---
            // Use the first two clusters in 'keep'
            int i1 = keep[0], i2 = keep[1];
        
            // Constants
            const double am = 0.938; // Proton mass [GeV/c^2]
            // Beam electron four-vector: (E, px, py, pz)
            double e0 = e0_nom;
            double pe_px = 0, pe_py = 0, pe_pz = e0_nom;
        
            // Scattered electron
            double ep = ep0_nom * (1.0 + dp / 100.0); // Scattered electron energy
            // HMS angles: th (vertical), ph (horizontal), both in radians
            // Assume small angles: px = p*sin(ph), py = p*sin(th), pz = p*cos(th)*cos(ph)
            double peprime_p = std::sqrt(ep*ep - 0.000511*0.000511); // Neglect electron mass
            double peprime_px = peprime_p * ph;
            double peprime_py = peprime_p * th;
            double peprime_pz = peprime_p * (1.0 - 0.5*(th*th + ph*ph)); // cos(th)*cos(ph) ≈ 1 - (th^2+ph^2)/2
        
            // Photon 1
            double E1 = cE[i1];
            double x1 = cX[i1], y1 = cY[i1];
            double r1 = std::sqrt(x1*x1 + y1*y1 + DNPS*DNPS);
            double p1x = E1 * x1 / r1;
            double p1y = E1 * y1 / r1;
            double p1z = E1 * DNPS / r1;
        
            // Photon 2
            double E2 = cE[i2];
            double x2 = cX[i2], y2 = cY[i2];
            double r2 = std::sqrt(x2*x2 + y2*y2 + DNPS*DNPS);
            double p2x = E2 * x2 / r2;
            double p2y = E2 * y2 / r2;
            double p2z = E2 * DNPS / r2;
        
            // Total outgoing energy and momentum
            double E_out = ep + E1 + E2;
            double px_out = peprime_px + p1x + p2x;
            double py_out = peprime_py + p1y + p2y;
            double pz_out = peprime_pz + p1z + p2z;
        
            // Initial state: beam electron + proton at rest
            double E_in = e0 + am;
            double px_in = 0.0;
            double py_in = 0.0;
            double pz_in = e0;
        
            // Missing mass squared
            double mm2 = std::pow(E_in - E_out, 2)
                       - std::pow(px_in - px_out, 2)
                       - std::pow(py_in - py_out, 2)
                       - std::pow(pz_in - pz_out, 2);
        
            double mm = (mm2 > 0) ? std::sqrt(mm2) : 0.0;
            hMissMass.Fill(mm);

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
// Set bin errors of hMean using variance from hVar (across ToyMC samples)
        for (int b = 1; b <= hMean->GetNbinsX(); ++b) {
            double var = hVar->GetBinContent(b);
            double err = (Ntoys > 1) ? std::sqrt(var / (Ntoys - 1)) : 0.0;
            hMean->SetBinError(b, err);
            // Save missing mass histogram to file
            TFile fMissMass("missing_mass.root", "RECREATE");
            hMissMass.Write();
            fMissMass.Close();
        }
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

        //────────────────────────────────── 7) Output: 6 Pads (restored ToyMC/missing mass/Q²)

        // Normalize all relevant histograms and background fit by realQ
        hBG.Scale(1.0 / realQ);
        hSig.Scale(1.0 / realQ);
        hSub.Scale(1.0 / realQ);
        std::vector<double> yF_norm(nBins);
        for (int i = 0; i < nBins; ++i) yF_norm[i] = yF[i] / realQ;

        // ToyMC/missing mass/Q² histograms are already normalized below (sf = 1.0/realQ)

        TCanvas c("c", "", 1200, 4100); // Increased height for 8 pads
        c.Divide(1, 8);
        
        // Pad 1: hAll (charge-normalized raw yield)
        c.cd(1);
        // Pad 1: Missing mass histogram (now at the bottom)
        std::cout << "hMissMass entries: " << hMissMass.GetEntries() << std::endl;
        std::cout << "hMissMass x-axis: [" << hMissMass.GetXaxis()->GetXmin() << ", " << hMissMass.GetXaxis()->GetXmax() << "]" << std::endl;
        std::cout << "First bin content: " << hMissMass.GetBinContent(1) << std::endl;
        std::cout << "Last bin content: " << hMissMass.GetBinContent(hMissMass.GetNbinsX()) << std::endl;
        std::cout << "Underflow: " << hMissMass.GetBinContent(0) << std::endl;
        std::cout << "Overflow: " << hMissMass.GetBinContent(hMissMass.GetNbinsX()+1) << std::endl;
        if (hMissMass.GetEntries() == 0) {
            // Test: fill with dummy value to check pad
            hMissMass.Fill(0.5);
            std::cout << "Filled hMissMass with test value 0.5" << std::endl;
        }
        hMissMass.SetLineColor(kRed+1);
        hMissMass.SetLineWidth(2);
        hMissMass.SetStats(1);
        double maxMM = hMissMass.GetMaximum();
        if (maxMM > 0) {
            hMissMass.SetAxisRange(0, maxMM * 1.2, "Y");
        }
        hMissMass.Draw("HIST");
        gPad->SetLogy(0); // Linear y-axis for missing mass, adjust as needed
        gPad->Update();

        // Pad 2: hAll (charge-normalized raw yield)
        c.cd(2);
        hAll.SetLineColor(kBlack);
        hAll.SetTitle("All clusT[0] (charge normalized);clusT[0];Counts/\\muC");
        hAll.GetXaxis()->SetRangeUser(100, hAll.GetXaxis()->GetXmax());
        hAll.Draw("HIST");

        // ... (existing pad drawing code for pads 2-7, now shifted up by 1) ...
        
            // Pad 8: Missing mass histogram
            // (Pad 8 is now used for the last plot, not missing mass. This block is removed.)

        // Pad 2: hBG (background region, normalized) + smoothed fit (yF, normalized)
        c.cd(3);
        hBG.SetLineColor(kBlack);
        hBG.SetMarkerStyle(20);
        hBG.SetMarkerSize(0.8);
        hBG.SetTitle("Background region (normalized) with smoothed fit;clusT[0];Counts/\\muC");
        hBG.Draw("P");
        TGraph* grBGfinal_norm = new TGraph(nBins, &vx[0], &yF_norm[0]);
        grBGfinal_norm->SetLineColor(kBlue);
        grBGfinal_norm->SetLineWidth(1);
        grBGfinal_norm->Draw("L SAME");

        // Prepare shifted background fit (normalized) for Pad 3
        std::vector<double> vys_shifted_norm(nBins);
        for (int i = 0; i < nBins; ++i) vys_shifted_norm[i] = vys[i] / realQ;
        TGraph* grBGfit_shifted_norm = new TGraph(nBins, &vxs[0], &vys_shifted_norm[0]);
    
        // Pad 3: hSig (signal, normalized), shifted background from spline (normalized), and hSub (subtracted, normalized)
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
        auto leg3 = new TLegend(0.6,0.7,0.88,0.88);
        leg3->AddEntry(&hSig,"Signal (hSig)","l");
        leg3->AddEntry(grBGfit_shifted_norm,"Background fit (shifted, norm)","l");
        leg3->AddEntry(&hSub_after,"Subtracted (hSub_after)","l");
        leg3->Draw();
        // Use shifted and normalized background fit
        grBGfit_shifted_norm->SetLineColor(kBlue);

        // Pad 7: Visualize dummy subtraction step
        c.cd(8);
        hSub_before.SetLineColor(kBlack);
        hSub_before.SetTitle("Dummy Subtraction Step;clusT[0];Counts/raw");
        hSub_before.SetLineStyle(1);
        hSub_before.Draw("HIST");
        hD.SetLineColor(kOrange+1);
        hD.SetLineStyle(2);
        hD.Draw("HIST SAME");
        hSub_after.SetLineColor(kRed);
        hSub_after.SetLineStyle(1);
        hSub_after.Draw("HIST SAME");
        auto leg7 = new TLegend(0.6,0.7,0.88,0.88);
        leg7->AddEntry(&hSub_before,"Before subtraction","l");
        leg7->AddEntry(&hD,"Scaled dummy","l");
        leg7->AddEntry(&hSub_after,"After subtraction","l");
        leg7->Draw();

        // Add legend for clarity
        auto leg = new TLegend(0.6, 0.7, 0.88, 0.88);
        leg->AddEntry(&hSig, "Signal (hSig)", "l");
        leg->AddEntry(grBGfit_shifted_norm, "Background (fit, shifted)", "l");
        leg->AddEntry(&hSub, "Subtracted (hSub)", "l");
        leg->Draw();

        // Pad 4: π0 mass mean (ToyMC)
        c.cd(5);
        hMean->SetLineColor(kBlack);
        hMean->SetTitle("ToyMC π^{0} mass mean;M_{#gamma#gamma} (GeV);Counts/\\muC");
        hMean->Draw("HIST");

        // Pad 5: π0 mass mean + error band (ToyMC)
        c.cd(6);
        double ymaxM = 1.2 * (hMean->GetMaximum() + hVar->GetMaximum());
        double simcMax = 0.0;
        
        // Overlay simulation histogram from simc_yield_histos.root
        TFile* fSimc = TFile::Open("simc_yield_histos.root");
        TH1* hSimcToDraw = nullptr;
        if (fSimc && !fSimc->IsZombie()) {
            TH1* hSimc = dynamic_cast<TH1*>(fSimc->Get("h_all"));
            if (hSimc) {
                // Print binning and range info for debugging
                std::cout << "hMean: bins=" << hMean->GetNbinsX()
                          << ", xmin=" << hMean->GetXaxis()->GetXmin()
                          << ", xmax=" << hMean->GetXaxis()->GetXmax() << std::endl;
                std::cout << "h_all: bins=" << hSimc->GetNbinsX()
                          << ", xmin=" << hSimc->GetXaxis()->GetXmin()
                          << ", xmax=" << hSimc->GetXaxis()->GetXmax() << std::endl;
                std::cout << "h_all: integral=" << hSimc->Integral() << std::endl;
        
                // Set axis range to match hMean/frame5
                hSimc->GetXaxis()->SetRangeUser(0, 0.30);
        
                // If binning doesn't match, rebin to match hMean
                if (hSimc->GetNbinsX() != hMean->GetNbinsX() ||
                    hSimc->GetXaxis()->GetXmin() != hMean->GetXaxis()->GetXmin() ||
                    hSimc->GetXaxis()->GetXmax() != hMean->GetXaxis()->GetXmax()) {
                    std::cerr << "Warning: h_all binning does not match hMean. Attempting to rebin." << std::endl;
                    int nBins = hMean->GetNbinsX();
                    double xMin = hMean->GetXaxis()->GetXmin();
                    double xMax = hMean->GetXaxis()->GetXmax();
                    TH1* hSimcRebinned = (TH1*)hSimc->Rebin(nBins / hSimc->GetNbinsX(), "hSimcRebinned");
                    hSimcRebinned->GetXaxis()->SetRangeUser(xMin, xMax);
                    hSimcToDraw = hSimcRebinned;
                } else {
                    hSimcToDraw = hSimc;
                }
                if (hSimcToDraw) {
                    hSimcToDraw->SetLineColor(kBlue);
                    hSimcToDraw->SetLineWidth(1);
                    hSimcToDraw->SetLineStyle(1); // solid
                    simcMax = hSimcToDraw->GetMaximum();
                }
            } else {
                std::cerr << "Could not find histogram 'h_all' in simc_yield_histos.root" << std::endl;
            }
        } else {
            std::cerr << "Could not open simc_yield_histos.root" << std::endl;
        }
        
        // Increase y-axis to accommodate both hMean/hVar and hSimc
        if (simcMax > 0.0) {
            double meanMax = hMean->GetMaximum() + hVar->GetMaximum();
            ymaxM = 1.2 * std::max(meanMax, simcMax);
        }
        if (ymaxM < 1e-6) ymaxM = 1e-6;
        
        TH2F frame5("f5", ";M_{#gamma#gamma} (GeV);Counts/\\muC", 10, 0, 0.30, 10, 0, ymaxM);
        frame5.Draw("AXIS");
        gBandM->Draw("F SAME");
        hMean->Draw("HIST SAME");
        if (hSimcToDraw) hSimcToDraw->Draw("HIST SAME");
        
        hMean->Fit("gaus", "SAME", "", 0.122, 0.142);
        TF1* fitFunc = hMean->GetFunction("gaus");
        if (fitFunc) fitFunc->SetLineWidth(1);
        if (fitFunc) fitFunc->Draw("SAME");
        // Custom stat box for fit and histogram
        if (fitFunc) {
            double chi2 = fitFunc->GetChisquare();
            int ndf = fitFunc->GetNDF();
            double mean = fitFunc->GetParameter(1);
// Calculate and display the integral of the histogram hMean from 0.12 to 0.14
int bin1 = hMean->FindBin(0.12);
int bin2 = hMean->FindBin(0.14);
double hist_integral = hMean->Integral(bin1, bin2);
TPaveText* pint = new TPaveText(0.6, 0.65, 0.88, 0.7, "NDC");
pint->SetFillColor(0);
pint->SetTextAlign(12);
pint->AddText(Form("Hist Int[0.12,0.14] = %.4f", hist_integral));
pint->Draw("SAME");
            double meanErr = fitFunc->GetParError(1);
            double sigma = fitFunc->GetParameter(2);
            double sigmaErr = fitFunc->GetParError(2);
            double constant = fitFunc->GetParameter(0);
            double constantErr = fitFunc->GetParError(0);

            TPaveText* stats = new TPaveText(0.6, 0.7, 0.88, 0.88, "NDC");
            stats->SetFillColor(0);
            stats->SetTextAlign(12);
            stats->AddText(Form("Entries = %.0f", hMean->GetEntries()));
            stats->AddText(Form("Mean = %.4f #pm %.4f", mean, meanErr));
            stats->AddText(Form("Sigma = %.4f #pm %.4f", sigma, sigmaErr));
            stats->AddText(Form("Const = %.4f #pm %.4f", constant, constantErr));
            stats->AddText(Form("#chi^{2}/NDF = %.2f / %d", chi2, ndf));
            stats->Draw("SAME");
        }
        gPad->Update();

        // Pad 6: Q² mean + error band (ToyMC)
        c.cd(7);
        double ymaxQ = 1.2 * (hQmean->GetMaximum() + hQvar->GetMaximum());
        if (ymaxQ < 1e-8) ymaxQ = 1e-8;
        TH2F frame6("f6", ";Q^{2} (GeV^{2});Counts/\\muC", 10, q2Lo, q2Hi, 10, 0, ymaxQ);
        frame6.Draw("AXIS");
        gBandQ->Draw("F SAME");
        hQmean->Draw("HIST SAME");

        c.Print(outPDF.c_str());
        std::cout << "All done. Canvas saved to " << outPDF << "\n";

        gSystem->Exit(0);
        return 0;
    }


