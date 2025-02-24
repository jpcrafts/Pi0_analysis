#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>

int main(int argc, char* argv[]){
  if(argc < 3){
    std::cerr << "Usage: " << argv[0] << " <input.root> <output.pdf>\n";
    return 1;
  }
  
  std::string inputFileName  = argv[1];
  std::string outputFileName = argv[2];
  
  // Open the ROOT file and get TTree "T"
  TFile* inFile = TFile::Open(inputFileName.c_str(), "READ");
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
  
  // Set branch addresses (including HMS cuts)
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
  
  // Define HMS and other cut thresholds
  const double edtmtdcCut = 0.1;
  const double hdeltaLow  = -8.5;
  const double hdeltaHigh =  8.5;
  const double hcaltotCut = 0.6;
  const double hcernpeCut = 1.0;
  const double gtrthCut   = 0.09;
  const double gtrphCut   = 0.09;
  
  // Create a RooDataSet for x in [111,144] ns covering all 17 peaks.
  RooRealVar x("x", "Cluster Time (ns)", 111, 144);
  RooArgSet varSet(x);
  RooDataSet data("data", "Data for 17 peaks", varSet);
  
  for(Long64_t i = 0; i < nEntries; i++){
    tree->GetEntry(i);
    // Apply HMS cuts:
    if(edtmtdc >= edtmtdcCut)         continue;
    if(hdelta  <= hdeltaLow)          continue;
    if(hdelta  >= hdeltaHigh)         continue;
    if(hcaltot <= hcaltotCut)         continue;
    if(hcernpe <= hcernpeCut)         continue;
    if(std::fabs(gtrth) > gtrthCut)   continue;
    if(std::fabs(gtrph) > gtrphCut)   continue;
    
    if(!vecNewEWClusT) continue;
    for(auto& t : *vecNewEWClusT){
      if(t >= 111 && t <= 144){
        x.setVal(t);
        data.add(varSet);
      }
    }
  }
  inFile->Close();
  std::cout << "After cuts, dataset has " << data.numEntries() 
            << " entries in [111,144] ns.\n";
  
  // Global estimates from previous 17-peak fit
  double globalFirstMean = 111.773;  // example value
  double globalSpacing   = 2.00104;   // example value
  double globalSigma     = 0.530157;  // example value
  
  // Loop over each peak index (0 to 16) and fit a single Gaussian.
  std::ofstream outfile("peak_fit_results.txt");
  if (!outfile) {
    std::cerr << "Error opening output file for writing.\n";
    return 1;
  }
  
  outfile << "Individual Peak Fits:\n";
  for (int i = 0; i < 17; i++){
    double predictedMean = globalFirstMean + i * globalSpacing;
    // Define a narrow subrange around the predicted mean:
    double lowRange = predictedMean - 0.8;
    double highRange = predictedMean + 0.8;
    RooRealVar x_peak("x_peak", "Cluster Time (ns)", lowRange, highRange);
    
    // Reduce the main dataset to this subrange:
    TString cut = Form("x >= %f && x <= %f", lowRange, highRange);
    RooAbsData* data_peak = data.reduce(cut);
    outfile << "Peak " << i << " using subrange [" << lowRange << ", " << highRange 
            << "] ns, entries = " << data_peak->numEntries() << "\n";
    
    // Define a Gaussian for this peak.
    RooRealVar mean("mean", "Peak Mean", predictedMean, predictedMean - 0.1, predictedMean + 0.1);
    RooRealVar sigma("sigma", "Peak Width", globalSigma, globalSigma * 0.8, globalSigma * 1.2);
    RooGaussian gauss("gauss", "Gaussian Peak", x_peak, mean, sigma);
    
    RooFitResult* fitRes = gauss.fitTo(*data_peak, RooFit::Save(true), RooFit::PrintLevel(0));
    if(fitRes) {
      outfile << "Peak " << i << " fitted mean = " << mean.getVal() 
              << " ns, sigma = " << sigma.getVal() << " ns\n";
    } else {
      outfile << "Peak " << i << " fit failed.\n";
    }
  }
  outfile.close();
  std::cout << "Peak fit results saved to peak_fit_results.txt\n";
  
  // Optionally, produce an overall high-resolution plot of the data:
  TCanvas c("c", "Overall Data", 800, 600);
  RooPlot* frame = x.frame();
  data.plotOn(frame, RooFit::Binning(2000));
  frame->SetTitle("Overall Data in [111,144] ns");
  frame->Draw();
  TCanvas c2("c2", "Overall Data", 800, 600);
  c2.SaveAs("overall_data.pdf");
  std::cout << "Saved overall data plot to overall_data.pdf\n";
  
  return 0;
}
