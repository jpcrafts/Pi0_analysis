#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooFormulaVar.h>
#include <RooPlot.h>
#include <RooFitResult.h>

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>

int main(int argc, char* argv[]){
  if(argc < 3){
    std::cerr << "Usage: " << argv[0] << " <input.root> <output.pdf>\n";
    return 1;
  }
  
  std::string inputFileName  = argv[1];
  std::string outputFileName = argv[2];
  
  // -------------------------------
  // 1. Open the ROOT file and get TTree "T"
  // -------------------------------
  TFile* inFile = TFile::Open(inputFileName.c_str(),"READ");
  if(!inFile || inFile->IsZombie()){
    std::cerr << "Error: cannot open " << inputFileName << std::endl;
    return 1;
  }
  TTree* tree = dynamic_cast<TTree*>(inFile->Get("T"));
  if(!tree){
    std::cerr << "Error: TTree 'T' not found in file " << inputFileName << std::endl;
    inFile->Close();
    return 1;
  }
  Long64_t nEntries = tree->GetEntries();
  std::cout << "Loaded TTree 'T' with " << nEntries << " entries.\n";
  
  // ---------------------------------------------------
  // 2. Set branch addresses (including HMS cuts variables)
  // ---------------------------------------------------
  double edtmtdc, hdelta, hcaltot, hcernpe, gtrth, gtrph, nclust_d;
  std::vector<double>* vecNewEWClusT = nullptr; // cluster times
  
  tree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
  tree->SetBranchAddress("H.gtr.dp",               &hdelta);
  tree->SetBranchAddress("H.cal.etotnorm",         &hcaltot);
  tree->SetBranchAddress("H.cer.npeSum",           &hcernpe);
  tree->SetBranchAddress("H.gtr.th",               &gtrth);
  tree->SetBranchAddress("H.gtr.ph",               &gtrph);
  tree->SetBranchAddress("NPS.cal.nclust",         &nclust_d);
  tree->SetBranchAddress("NPS.cal.newEWClusT",     &vecNewEWClusT);
  
  // -------------------------------
  // 3. Define HMS and other cut thresholds
  // -------------------------------
  const double edtmtdcCut = 0.1;
  const double hdeltaLow  = -8.5;
  const double hdeltaHigh =  8.5;
  const double hcaltotCut = 0.6;
  const double hcernpeCut = 1.0;
  const double gtrthCut   = 0.09;
  const double gtrphCut   = 0.09;
  
  // -------------------------------
  // 4. Create a RooDataSet for x in [111,144] ns.
  // This range should cover all 17 peaks.
  // -------------------------------
  RooRealVar x("x","Cluster Time (ns)", 111, 144.5);
  RooArgSet varSet(x);
  RooDataSet data("data","Data for first 17 peaks", varSet);
  
  for(Long64_t i = 0; i < nEntries; i++){
    tree->GetEntry(i);
    // Apply HMS and event-level cuts:
    if(edtmtdc >= edtmtdcCut)         continue;
    if(hdelta  <= hdeltaLow)          continue;
    if(hdelta  >= hdeltaHigh)         continue;
    if(hcaltot <= hcaltotCut)         continue;
    if(hcernpe <= hcernpeCut)         continue;
    if(std::fabs(gtrth) > gtrthCut)   continue;
    if(std::fabs(gtrph) > gtrphCut)   continue;
    
    if(!vecNewEWClusT) continue;
    for(auto& t : *vecNewEWClusT){
      if(t >= 111 && t <= 144.5){
        x.setVal(t);
        data.add(varSet);
      }
    }
  }
  inFile->Close();
  std::cout << "After cuts, dataset has " << data.numEntries() 
            << " entries in [111,144] ns.\n";
  
  // -----------------------------------------------------
  // 5. Build the signal model: Sum of 17 Gaussians for the first 17 peaks.
  // Each peak's mean is defined by:
  //   mean_i = firstMean + i * spacing   for i = 0,1,...,16.
  // We initialize firstMean near 111.3 ns and spacing near 2 ns.
  // The common width (smallSigma) is constrained around 0.5 ns.
  // -----------------------------------------------------
  int nPeaks = 17;
  RooRealVar firstMean("firstMean", "Mean of first peak", 111.78, 111.7, 111.8);
  RooRealVar spacing("spacing", "Spacing between peaks", 2.0, 1.9, 2.1);
  RooRealVar smallSigma("smallSigma", "Common width", 0.5, 0.45, 0.55);
  
  RooArgList gaussList;
  RooArgList fracList;
  for (int i = 0; i < nPeaks; i++){
    TString formula = Form("@0 + %d*@1", i);
    RooFormulaVar* mean_i = new RooFormulaVar(Form("mean%d", i),
                                               "Peak mean",
                                               formula,
                                               RooArgList(firstMean, spacing));
    RooGaussian* gauss = new RooGaussian(Form("gauss%d", i),
                                         Form("Peak %d", i),
                                         x,
                                         *mean_i,
                                         smallSigma);
    gaussList.add(*gauss);
    // For i < nPeaks-1, define fraction parameters fixed to 1/nPeaks.
    if(i < nPeaks - 1){
      RooRealVar* fracVar = new RooRealVar(Form("frac%d", i),
                                           "Fraction for peak",
                                           1.0/double(nPeaks));
      fracVar->setConstant(true);
      fracList.add(*fracVar);
    }
  }
  RooAddPdf signal("signal", "Sum of 17 peaks", gaussList, fracList);
  
  // -----------------------------------------------------
  // 6. Fit the signal model to the data.
  // -----------------------------------------------------
  RooFitResult* fitRes = signal.fitTo(data, RooFit::Save(true), RooFit::PrintLevel(1));
  if(fitRes) fitRes->Print("v");
  else std::cerr << "Fit failed or was interrupted.\n";
  
  // Compute the derived peak means:
  std::cout << "Fitted firstMean = " << firstMean.getVal() << " ns\n";
  std::cout << "Fitted spacing = " << spacing.getVal() << " ns\n";
  for (int i = 0; i < nPeaks; i++){
    double peakMean = firstMean.getVal() + i * spacing.getVal();
    std::cout << "Peak" << i << " mean = " << peakMean << " ns\n";
  }
  std::cout << "Fitted sigma = " << smallSigma.getVal() << " ns\n";
  
  // -----------------------------------------------------
  // 7. Plot the fit result with high resolution.
  // -----------------------------------------------------
  TCanvas c("c", "Seventeen-Peak Fit", 800, 600);
  RooPlot* frame = x.frame();
  // Use 2000 bins for high resolution.
  data.plotOn(frame, RooFit::Binning(2000));
  signal.plotOn(frame);
  frame->SetTitle("Fit of 17 Peaks in [111,144] ns");
  frame->Draw();
  c.SaveAs(outputFileName.c_str());
  
  std::cout << "Saved plot to " << outputFileName << std::endl;
  return 0;
}
