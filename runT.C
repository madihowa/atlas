#include "convert_csv_ttree.C"
#include "ClusterTreeCat.C"
#include <TStyle.h>
#include <TCanvas.h>
#include "PlotHisto.C"

void runT(){
  const char * name_out = "results.root";
  if (!gSystem->IsFileInIncludePath(name_out)) convert("results.csv", name_out);
  ClusterTree cluster;
  // topo.InitHist()
  cluster.Loop();
  PlotHisto(gSystem->pwd());
}
