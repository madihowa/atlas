#include "convert_csv_ttree.C"
#include "ClusterTreeCat.C"
#include <TStyle.h>
#include <TCanvas.h>
#include "PlotHisto.C"

void runTStage2(){


  convert("results.csv", "results.root");

  PlotHisto(gSystem->pwd());
}
