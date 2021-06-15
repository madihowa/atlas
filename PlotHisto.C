#include "PlotHisto.h"



void PlotHisto(const char *dir = "Results_2020-07-17_run_11"){

  gROOT->SetBatch(kTRUE);

  gSystem->cd(dir);

  std::cout << gSystem->pwd() << endl;
  TFile *f = new TFile("results.root", "UPDATE");
  TTree *t1, *t2, *t3;

  Bounds *b1 = new Bounds;



  f->GetObject("EM_tree", t2);
  f->GetObject("Had_tree", t1);
  f->GetObject("ClusterTree", t3);

  auto h_EM = Plot_performance(t2, 1);
  gPad->Print("PerformanceEM.png");
  auto h_Had = Plot_performance(t1, 2);
  gPad->Print("PerformanceHad.png");
  auto h_all_old = Plot_performance_Old(t3, 0);
  auto h_all = Plot_performance(t3, 0);

  // std::cout << "Number of Hadronic entries: " << t1->GetEntriesFast() << endl;
  std::cout << "Number of Electromagetic entries: " << t2->GetEntriesFast() << endl;
  std::cout << "Number of entries: " << t3->GetEntriesFast() << endl;
  std::cout << "Number of entires in ratio branch: " << t3->GetBranch("Ratio")->GetEntries() << endl;
  //
  //
  //
  //

  // return;

  auto g_Had = GraphMean(h_Had);
  auto g_EM = GraphMean(h_EM);
  auto g_all = GraphMean(h_all);
  auto g_all_old = GraphMean(h_all_old);


  //
  // g_Had->SetXTitle("True Energy [GeV]");
  // g_Had->SetYTitle("Median of Response");
  //

  Int_t color_Had = 2;
  // g_Had->SetLineColor(color_Had);

  Int_t color_EM = 3;
  g_EM->SetLineColor(color_EM);

  Int_t color_all_old = 3;
  g_all_old->SetLineColor(color_all_old);


  Int_t color_all = 4;
  g_all->SetLineColor(color_all);

  gStyle->SetOptStat(0);
  AtlasStyle();
  SetAtlasStyle();
  TCanvas * c_Mean = new TCanvas("c_Ratio", "Means");
  TMultiGraph *mg = new TMultiGraph();
  c_Mean->SetLogx();
  // mg->Add(g_Had);
  // mg->Add(g_EM);
  mg->Add(g_all);
  mg->Add(g_all_old);

  mg->Draw("AC");
  mg->GetXaxis()->SetTitle("Energy [GeV]");
  mg->GetYaxis()->SetTitle("#frac{Calibrated Energy}{Cluster Energy}");
  // g_EM->Draw("same AC");
  // g_all->Draw("same AC");
  // myLineText(.60, .90, .02, color_Had, "Hadronic Showers");
  myLineText(.60, .85, .02, color_all_old, "LCW Calibration");
  myLineText(.60, .80, .02, color_all, "RNN all inputs");


  gPad->Print("Mean.png");


  b1->maxX = t3->GetMaximum("cluster_ENG_CALIB_TOT");
  b1->maxX = 10;
  b1->minX = t3->GetMinimum("cluster_ENG_CALIB_TOT");

  b1->maxY = t3->GetMaximum("clusterE");
  b1->maxY = 10;
  b1->minY = t3->GetMinimum("clusterE");;


  std::cout << b1->maxX << b1->minX << b1->maxY << b1->minY << endl;
  TH2D * temp = Plot_TH2D(t3->GetLeaf("cluster_ENG_CALIB_TOT"), t3->GetLeaf("clusterE"), "True Energy [GeV]", "Cluster Energy [GeV]", b1, 1, 1, 100);
  gPad->Print("EngvsEng.png");

  //
  // Plot_Energy(t3);
  TH2D *h1 = Plot_TH2D(t3->GetLeaf("cluster_ENG_CALIB_TOT"), t3->GetLeaf("CalibratedE"), "True Energy [GeV]", "Cluster Energy [GeV]", b1, 1, 1, 100);
  gPad->Print("AAAA.png");
  TCanvas * c_IQR = new TCanvas("c_Ratio", "Means");
  c_IQR->SetLogx();


  auto g_IQR = IQR(h_all);
  TMultiGraph *mg_2 = new TMultiGraph();
  mg_2->Add(g_IQR);
  mg_2->Draw("AC");
  mg_2->GetXaxis()->SetTitle("Energy [GeV]");
  mg_2->GetYaxis()->SetTitle("#frac{IQR}{mean}");

  gPad->Print("IQR.png");
  // Plot_Calib(t1->GetLeaf("cluster_ENG_CALIB_TOT"), t1->GetLeaf("CalibratedE"));
  // gPad->Print("plot1.png");
  //
  //
  // Plot_Calib(t2->GetLeaf("cluster_ENG_CALIB_TOT"), t2->GetLeaf("CalibratedE"));
  // gPad->Print("plot.png");
  //

}
