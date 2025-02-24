#include <TCanvas.h>
#include <RooRealVar.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooPlot.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <fstream>
#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[]){
  // This code builds an analytic background model from fixed parameters 
  // derived from a previous 17-peak fit, using the central 15 peaks (indices 1 through 15).
  // It then evaluates the background model on a fine grid and writes the (x, f(x))
  // values to a text file, and produces a high-resolution plot.
  
  if(argc < 2){
    std::cerr << "Usage: " << argv[0] << " <output.pdf>\n";
    return 1;
  }
  std::string outputFileName = argv[1];
  
  // Define the time variable x over the range covering the central 15 peaks.
  // For example, if the central 15 peaks lie approximately between 113 ns and 142 ns:
  RooRealVar x("x", "Time (ns)", 113, 142);
  
  // Set fixed parameters from the previous 17-peak fit.
  RooRealVar firstMean("firstMean", "Mean of first full peak", 111.773, 111.773, 111.773);
  firstMean.setConstant(true);
  
  RooRealVar spacing("spacing", "Spacing between peaks", 2.00104, 2.00104, 2.00104);
  spacing.setConstant(true);
  
  RooRealVar smallSigma("smallSigma", "Common width", 0.530157, 0.530157, 0.530157);
  smallSigma.setConstant(true);
  
  // Build the background model using the central 15 peaks (peak indices 1 to 15).
  int nBgPeaks = 15;
  RooArgList gaussList;
  RooArgList fracList;
  
  for (int i = 1; i <= nBgPeaks; i++){
    // The mean for the i-th peak is given by: firstMean + i*spacing.
    TString formula = Form("@0 + %d*@1", i);
    RooFormulaVar* mean_i = new RooFormulaVar(Form("mean%d", i),
                                               "Background Peak Mean",
                                               formula,
                                               RooArgList(firstMean, spacing));
    RooGaussian* gauss = new RooGaussian(Form("gauss%d", i),
                                         "Background Peak",
                                         x,
                                         *mean_i,
                                         smallSigma);
    gaussList.add(*gauss);
    // For peaks 1 to 14, provide fraction parameters fixed to 1/15.
    if(i < nBgPeaks){
      RooRealVar* fracVar = new RooRealVar(Form("frac%d", i),
                                           "Fraction for peak",
                                           1.0/double(nBgPeaks));
      fracVar->setConstant(true);
      fracList.add(*fracVar);
    }
  }
  
  RooAddPdf bgModel("bgModel", "Background Model (central 15 peaks)", gaussList, fracList);
  
  // Fix the coefficient normalization using the observable x.
  bgModel.fixCoefNormalization(RooArgSet(x));
  
  // Evaluate the background model over a fine grid.
  int nPoints = 100;
  double x_min = 113;
  double x_max = 142;
  double dx = (x_max - x_min) / (nPoints - 1);
  
  std::ofstream outfile("background_model_values.txt");
  if(!outfile){
    std::cerr << "Error opening output file for writing.\n";
    return 1;
  }
  outfile << "x (ns)\tbgModel(x)\n";
  for (int i = 0; i < nPoints; i++){
    double x_val = x_min + i * dx;
    x.setVal(x_val);
    // Use getValV() workaround:
    double val = bgModel.getVal();
    outfile << x_val << "\t" << val << "\n";
  }
  outfile.close();
  std::cout << "Background model values saved to background_model_values.txt\n";
  
  // Plot the background model.
  TCanvas c("c", "Background Model", 800, 600);
  RooPlot* frame = x.frame();
  // Instead of SetNpx(), simply plot with the default resolution.
  bgModel.plotOn(frame);
  frame->SetTitle("Background Model from Central 15 Peaks");
  frame->Draw();
  c.SaveAs(outputFileName.c_str());
  
  std::cout << "Background model saved to " << outputFileName << std::endl;
  return 0;
}
