#ifdef __CLING__
// these are not headers - do not treat them as such - needed for ROOT6
#include "AtlasLabels.C"
#include "AtlasUtils.C"
#endif
#include "AtlasStyle.h"
#include "AtlasStyle.C"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// #include <TMath.h>

class Bounds
{
public:
  Double_t maxX;
  Double_t maxY;
  Double_t minX;
  Double_t minY;
  const char *TitleX;
  const char *TitleY;



};
class TitleAxis{
public:
  const char *TitleX;
  const char *TitleY;
  const char *AxisTitle;
  const char *SaveName;
  Bool_t LogY = 0;
  Bool_t LogX = 0;
};

TH2D * Plot_performance(TTree *t, Int_t switch_name);
void Plot_Calib(TLeaf *First, TLeaf *Second);
// void Plot_TH2D(TLeaf *First, TLeaf *Second, const char *titleX, const char *titleY, Bounds *bound, Bool_t LogX = "False", Bool_t LogY = "False", Int_t n_bins = 100);
Double_t sum(TLeaf *l);
void Plot_TH1D(TLeaf *l1, TLeaf *l2, Bounds *b);
TH2D * TH2Subset(TH2D *h,Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2);
TGraph * GraphMean(TH2D* h);
void Plot_Energy(TTree *t);
TitleAxis * GetTit(const char* branch_name);
TH2D * Plot_performance_Old(TTree *t, Int_t switch_name);




TH2D * Plot_performance(TTree *t, Int_t switch_name){

  // TFile *f1 = new TFile("FinalResults.root");



  gStyle->SetOptStat(0);
  AtlasStyle();
  SetAtlasStyle();
  gStyle->SetPadRightMargin(0.16);

  Double_t cluster_ENG_CALIB_TOT, CalibratedE, Ratio;
  // TLeaf *l_performance = new TLeaf();
  const char *list = "performance";
  TBranch *b_Ratio = t->Branch("Ratio", &Ratio, "Ratio/D");
  // TBranch *b_Ratio = new TBranch ();
  t->SetBranchAddress("Ratio", &Ratio, &b_Ratio);
  // TLeaf *l_true = t->GetLeaf("cluster_ENG_CALIB_TOT");
  TLeaf *l_true = t->GetLeaf("cluster_ENG_CALIB_TOT");
  TLeaf *l_calib = t->GetLeaf("CalibratedE");




  Int_t n1 = l_true->GetBranch()->GetEntries();
  Int_t n2 = l_calib->GetBranch()->GetEntries();

  if (n1 != n2){
    std:cout << "Error something majorly wrong. different number of entries in different brnaches" << endl;
    }

  Int_t n_entries = n1;

  Double_t MinX = .05;
  Double_t MaxX = t->GetMaximum("cluster_ENG_CALIB_TOT");
  // MaxX = 1;
  Double_t MaxY = 2;
  Double_t MinY = 0;
  Int_t n_bins = 100;
  Double_t LogWidth [n_bins + 1];
  Double_t yWidth [n_bins + 1];
  auto c_Ratio = new TCanvas("c_Ratio", "performance");
  c_Ratio->SetLogx();
  c_Ratio->SetLogz();
  // c_Ratio->SetLogy();
  // std:cout << "Min Y; " << Min
  for (Int_t i = 0; i <= n_bins; i++){
    LogWidth[i] = pow(10, log10(MinX) + (log10(MaxX) - log10(MinX))/double(n_bins)*i);
    // yWidth[i] = pow(10, log10(MinY) + (log10(MaxY) - log10(MinY))/double(n_bins)*i);
    yWidth[i] = MinY + i*(MaxY - MinY)/n_bins;
    // yWidth[i] = MinY + i*(MaxY - MinY)/n_bins;
  }
  auto h_ratio = new TH2D("h_ratio", "performance", n_bins, LogWidth, n_bins, yWidth);


  for (Int_t i = 0; i < n_entries; i++){
    t->GetEntry(i);
    cluster_ENG_CALIB_TOT = l_true->GetValue(0);
    CalibratedE = l_calib->GetValue(0);
    Ratio = CalibratedE/cluster_ENG_CALIB_TOT;
    b_Ratio->Fill();
    // std::cout << "Filling Ratio with: " <<b_Ratio->GetLeaf("Ratio")->GetValue(0)<< endl;
    h_ratio->Fill(cluster_ENG_CALIB_TOT, Ratio);

  }
  h_ratio->SetXTitle("Cluster Energy [GeV]");
  h_ratio->SetYTitle("#frac{Calibrated Energy}{Cluster Energy}");
  gStyle->SetPalette(112);
  h_ratio->Draw("COLZ");

  TGraph * g_mean = GraphMean(h_ratio);
  g_mean->SetLineColor(2);
  g_mean->Draw("same");

  // h_ratio->GetYaxis()->SetLimits(0, 2);
  if (switch_name == 0){
    myText(.55, .85, 2, "#pi^{0} and #pi^{#pm} clusters");
    c_Ratio->Print("performance.png");}
  if (switch_name == 1){
    myText(.55, .85, 2, "#pi^{0} clusters");
    c_Ratio->Print("EMperformance.png");}
  if (switch_name == 2){
    myText(.55, .85, 2, "#pi^{#pm} clusters");
    c_Ratio->Print("HADperformance.png");}

    // h_ratio->GetYaxis()->SetLimits(0, t->GetMaximum("Ratio"));

    t->Write();

  return h_ratio;

}

TH2D * Plot_performance_Old(TTree *t, Int_t switch_name){

  // TFile *f1 = new TFile("FinalResults.root");



  gStyle->SetOptStat(0);
  AtlasStyle();
  SetAtlasStyle();
  gStyle->SetPadRightMargin(0.16);

  Double_t cluster_ENG_CALIB_TOT, clusterECalib, Ratio_old;
  // TLeaf *l_performance = new TLeaf();
  const char *list = "performance";
  TBranch *b_Ratio_Old = t->Branch("Ratio_old", &Ratio_old, "Ratio_old/D");
  t->SetBranchAddress("Ratio_old", &Ratio_old, &b_Ratio_Old);
  TLeaf *l_true = t->GetLeaf("cluster_ENG_CALIB_TOT");
  TLeaf *l_calib = t->GetLeaf("clusterECalib");




  Int_t n1 = l_true->GetBranch()->GetEntries();
  Int_t n2 = l_calib->GetBranch()->GetEntries();

  if (n1 != n2){
    std:cout << "Error something majorly wrong. different number of entries in different brnaches" << endl;
    }

  Int_t n_entries = n1;

  Double_t MinX = .05;
  Double_t MaxX = t->GetMaximum("cluster_ENG_CALIB_TOT");
  // MaxX = 1;
  Double_t MaxY = 2;
  Double_t MinY = 0;
  Int_t n_bins = 100;
  Double_t LogWidth [n_bins + 1];
  Double_t yWidth [n_bins + 1];
  auto c_Ratio = new TCanvas("c_Ratio", "performance");
  c_Ratio->SetLogx();
  c_Ratio->SetLogz();
  // c_Ratio->SetLogy();
  // std:cout << "Min Y; " << Min
  for (Int_t i = 0; i <= n_bins; i++){
    LogWidth[i] = pow(10, log10(MinX) + (log10(MaxX) - log10(MinX))/double(n_bins)*i);
    // yWidth[i] = pow(10, log10(MinY) + (log10(MaxY) - log10(MinY))/double(n_bins)*i);
    yWidth[i] = MinY + i*(MaxY - MinY)/n_bins;
    // yWidth[i] = MinY + i*(MaxY - MinY)/n_bins;
  }
  auto h_ratio_old = new TH2D("h_ratio_old", "performance", n_bins, LogWidth, n_bins, yWidth);


  for (Int_t i = 0; i < n_entries; i++){
    t->GetEntry(i);
    cluster_ENG_CALIB_TOT = l_true->GetValue(0);
    clusterECalib = l_calib->GetValue(0);
    Ratio_old = clusterECalib/cluster_ENG_CALIB_TOT;
    b_Ratio_Old->Fill();
    h_ratio_old->Fill(cluster_ENG_CALIB_TOT, Ratio_old);

  }
  h_ratio_old->SetXTitle("Cluster Energy [GeV]");
  h_ratio_old->SetYTitle("#frac{LCW Calibrated Energy}{Cluster Energy}");
  gStyle->SetPalette(112);
  h_ratio_old->Draw("COLZ");

  TGraph * g_mean = GraphMean(h_ratio_old);
  g_mean->SetLineColor(2);
  g_mean->Draw("same");

  if (switch_name == 0){
    myText(.55, .85, 1, "#pi^{0} and #pi^{#pm} clusters");
    c_Ratio->Print("LCWperformance.png");}
  if (switch_name == 1){
    myText(.55, .85, 1, "#pi^{0} clusters");
    c_Ratio->Print("LCWEMperformance.png");}
  if (switch_name == 2){
    myText(.55, .85, 1, "#pi^{#pm} clusters");
    c_Ratio->Print("LCWHADperformance.png");}





  return h_ratio_old;
}

void Plot_TH1D(TLeaf *l1, TLeaf *l2, Bounds *b){

  gStyle->SetOptStat(0);
  AtlasStyle();
  SetAtlasStyle();

  auto c1 = new TCanvas("c1", "Fit");
  auto h1 = new TH1D("h1", "Title", 1001, b->minX, b->maxX);
  auto h2 = new TH1D("h2", "Title?", 1001, b->minX, b->maxX);

  Int_t nentries_1 = l1->GetBranch()->GetEntries();
  Int_t nentries_2 = l2->GetBranch()->GetEntries();


  TTree *t1 = l1->GetBranch()->GetTree();
  TTree *t2 = l2->GetBranch()->GetTree();

  TLeaf *Pro1 = t1->GetBranch("EM_Pro")->GetLeaf("EM_Pro");
  TLeaf *Pro2 = t2->GetBranch("EM_Pro")->GetLeaf("EM_Pro");

  Double_t sum_probabilities1 = sum(Pro1);
  Double_t sum_probabilities2 = sum(Pro2);

  sum_probabilities1 = nentries_1 - sum_probabilities1;
  sum_probabilities2 = nentries_2 - sum_probabilities2;

  Double_t temp_scale;

  for (Int_t i =0; i < nentries_1; i++){
    l1->GetBranch()->GetEntry(i);
    t1->GetEntry(i);
    // std::cout << Pro1->GetValue() << endl;
    h1->Fill(l1->GetValue(),(1.0-Pro1->GetValue())/sum_probabilities1);
  }



  for (Int_t i =0; i < nentries_1; i++){
    l2->GetBranch()->GetEntry(i);
    t2->GetEntry(i);
    h2->Fill(l2->GetValue(), (1.0-Pro2->GetValue())/sum_probabilities2);
  }


  h1->SetXTitle("Energy");
  h1->SetYTitle("");

  Int_t color_h1 = 4, color_h2 = 8;

  // h1->SetFillColor(3);
  // h1->SetMarkerStyle(0);
  // h1->SetLineStyle(1);

  h1->SetLineColor(color_h1);
  h2->SetLineColor(color_h2);


  h2->Draw("HIST");
  h1->Draw("HIST same");
  myText(       0.57, 0.90, 1, "Shower Probability");
  myLineText(    0.55, 0.80, 0.05, color_h1, "H1");
  myLineText(    0.55, 0.70, 0.05, color_h2, "H2");
  // myMarkerText( 0.55, 0.75, 1, 20, "Data 2009",1.3);

  return;


}

void Plot_Calib(TLeaf *First, TLeaf *Second){
  gStyle->SetOptStat(0);
  AtlasStyle();
  SetAtlasStyle();

  Double_t xMax = First->GetMaximum();

  if (xMax < Second->GetMaximum()) xMax = Second->GetMaximum();
  auto c1 = new TCanvas("c1", "Fit");
  auto h1 = new TH1D("h1", "EM Probability", 101, 0, xMax);
  auto h2 = new TH1D("h2", "Title?", 101, 0, xMax);

  Int_t nentries_1 = First->GetBranch()->GetEntries();
  for (Int_t i =0; i < nentries_1; i++){
    First->GetBranch()->GetEntry(i);
    h1->Fill(First->GetValue(), 1.0/(Double_t)nentries_1);
  }
  Int_t nentries_2 = First->GetBranch()->GetEntries();
  for (Int_t i =0; i < nentries_2; i++){
    Second->GetBranch()->GetEntry(i);
    h2->Fill(Second->GetValue(),1.0/(Double_t)nentries_2);
  }


  h1->SetXTitle("Energy");
  h1->SetYTitle("Percentage of clusters");
  h2->SetXTitle("Energy");
  h2->SetYTitle("Percentage of clusters");


  Int_t color_h1 = 4, color_h2 = 8;

  // h1->SetFillColor(3);
  // h1->SetMarkerStyle(0);
  // h1->SetLineStyle(1);

  h1->SetLineColor(color_h1);
  h2->SetLineColor(color_h2);


  h1->Draw("HIST same");
  h2->Draw("HIST same");
  // myText(       0.57, 0.90, 1, "Shower Probability");
  myText(       0.5, 0.87, 1, "Hadronic Clusters (#pi_{#pm})");
  myLineText(    0.55, 0.80, 0.05, color_h1, "Calibrated Energy");
  myLineText(    0.55, 0.72, 0.05, color_h2, "True Energy");
  return;

}

TH2D * Plot_TH2D(TLeaf *First, TLeaf *Second, const char *titleX, const char *titleY, Bounds *bound, Bool_t LogX = 0, Bool_t LogY = 0, Int_t n_bins = 100){
  gStyle->SetOptStat(0);
  AtlasStyle();
  SetAtlasStyle();
  gStyle->SetPadRightMargin(0.16);


  auto c2 = new TCanvas("c2", "Fit");

  Double_t xWidth[n_bins+1];
  Double_t yWidth[n_bins+1];

  Double_t xData[n_bins], yData[n_bins];

  if (LogX){c2->SetLogx(); std::cout << "Setting Logx\n"; }
  if (LogY){c2->SetLogy();}

  for (Int_t i = 0; i <= n_bins; i++){
    // yWidth[i] = pow(10, log10(MinY) + (log10(MaxY) - log10(MinY))/double(n_bins)*i);
    yWidth[i] = bound->minY + i*(bound->maxY - bound->minY)/n_bins;
    xWidth[i] = bound->minX + i*(bound->maxX - bound->minX)/n_bins;
    // std::cout << yWidth[i] <<endl;
    if (LogX) {xWidth[i] = pow(10, log10(bound->minX) + (log10(bound->maxX) - log10(bound->minX))/double(n_bins)*i);}
    if (LogY) {yWidth[i] = pow(10, log10(bound->minY) + (log10(bound->maxY) - log10(bound->minY))/double(n_bins)*i);}
  }
  // auto h_Mean = new TH1D("h_Mean", "", n_bins, xWidth);
  auto h3 = new TH2D("h3", "EM Probability", n_bins, xWidth, n_bins, yWidth);
  h3->Sumw2();


  Int_t nentries = First->GetBranch()->GetEntries();
  for (Int_t i =0; i < nentries; i++){
    First->GetBranch()->GetEntry(i);
    Second->GetBranch()->GetEntry(i);
    h3->Fill(First->GetValue(),Second->GetValue());
  }

  // for (Int_t i = 0; i < n_bins; i++){
  //   yData[i] = (TH2Subset(h3, i, i+1, i, i+1)->GetMean(1));
  //   xData[i] = xWidth[i];
  // }
  //
  // auto Array_weights = h3->GetSumw2();

  // auto h_Mean = h3->ProjectionY("projection", 0, n_bins);

  // TGraph *g1 = new TGraph(n_bins, xData,  yData);

  h3->SetXTitle(titleX);
  h3->SetYTitle(titleY);
  h3->SetLineColor(38);
  h3->Draw("COLZ");
  // h_Mean->Draw("same");
  // g1->Draw("same");
  // myLineText(    0.55, 0.80, 0.05, 38, "title");
  return h3;

}

TGraph * GraphMean(TH2D* h){
  Double_t ColumnSum = 0, ColumnWeightedSum = 0;
  auto n_BinsX = h->GetNbinsX();
  auto n_BinsY = h->GetNbinsY();
  auto Xaxis = h->GetXaxis();
  auto Yaxis = h->GetYaxis();
  Double_t ArrX[n_BinsX];
  Double_t ArrY[n_BinsX];
  Int_t n_NAN = 1;
  for(Int_t i = 1; i <= n_BinsX; i++){
    ColumnSum = 0;
    ColumnWeightedSum = 0;
    for(Int_t j = 0; j < n_BinsY; j++){
      ColumnWeightedSum += Yaxis->GetBinCenter(j) * h->GetBinContent(i, j);
      ColumnSum += h->GetBinContent(i, j);
      if (ColumnWeightedSum > 10000){std::cout << i << "\t " << j << "\t " << Yaxis->GetBinCenter(j) << "\t " << h->GetBinContent(i, j) << endl;}
    }
    // std::cout << i << "\t" << ColumnWeightedSum << "\t" << ColumnSum << endl;
    if (ColumnSum != 0) ArrY[i] = (Double_t)ColumnWeightedSum / (Double_t)ColumnSum;
    else {
      ArrY[i] = 0;
      n_NAN = i+1;
    }
    if (ArrY[i] > 2) n_NAN += 1;
    ArrX[i] = Xaxis->GetBinCenter(i);
    std::cout << ArrY[i] << '\t' << ArrX[i] << endl;
  }
  TGraph * g = new TGraph(n_BinsX-n_NAN, &ArrX[n_NAN], &ArrY[n_NAN]);

  g->SetLineColor(2);

  return g;




}


Double_t sum(TLeaf *l){

  Int_t nentries = l->GetBranch()->GetEntries();
  Double_t sum = 0;
  for (Int_t i = 0; i < nentries; i++){
    l->GetBranch()->GetEntry(i);
    sum += l->GetValue();
    // if (l->GetValue() != 1) std::cout << l->GetBranch()->GetEntries() << endl;
    // std::cout << "iteration: " << i << endl;
  }
  // std::cout << "The sum is: "<<sum;
  return sum;
}
TH2D * TH2Subset(TH2D *h,Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2){
   Double_t x1 = h->GetXaxis()->GetBinLowEdge(binx1);
   Double_t x2 = h->GetXaxis()->GetBinLowEdge(binx2)+h->GetXaxis()->GetBinWidth(binx2);
   Double_t y1 = h->GetYaxis()->GetBinLowEdge(biny1);
   Double_t y2 = h->GetYaxis()->GetBinLowEdge(biny2)+h->GetYaxis()->GetBinWidth(biny2);
   Int_t    n1 = binx1+binx2+1;
   Int_t    n2 = biny1+biny2+1;

   Int_t nbinsx = h->GetNbinsX();
   Int_t nbinsy = h->GetNbinsY();

   TH2D *hs = new TH2D(h->GetName(),
                       h->GetTitle(),
                       n1, x1, x2,
                       n2, y1, y2);

   Double_t content, x, y;

   for (Int_t i=1; i<=nbinsx; i++) {
      for (Int_t j=1; j<=nbinsx; j++) {
         content = h->GetBinContent(i,j);
         x = h->GetXaxis()->GetBinCenter(i);
         y = h->GetYaxis()->GetBinCenter(j);
         hs->Fill(x, y, content);
      }
   }

   return hs;
}

void Plot_Energy(TTree *t){

  Bounds *b = new Bounds;

  Int_t n_Branch = t->GetListOfLeaves()->GetEntriesFast();



  // TBranch *b_Energy = t->GetLeaf("cluster_ENG_CALIB_TOT");

  TitleAxis *tit = new TitleAxis;
  b->minX = t->GetMinimum("cluster_ENG_CALIB_TOT");
  b->maxX = t->GetMaximum("cluster_ENG_CALIB_TOT");
  tit->TitleX = "Energy [Gev]";

  TString PNG_name;
  TString str1 = (TString)".png";
  TString str2 = (TString)"Energy_";

  for (Int_t i = 0; i < n_Branch - 1; i++){
    // std::cout << t->GetListOfLeaves()->At(i)->GetName() << "\t";

    tit = GetTit(t->GetListOfLeaves()->At(i)->GetName());
    if (!strcmp(t->GetListOfLeaves()->At(i)->GetName(), "Ratio_old")) continue;

    std::cout << t->GetListOfLeaves()->At(i)->GetName() << endl;
    std::cout << "Entries found: " << t->GetBranch(t->GetListOfLeaves()->At(i)->GetName())->GetEntries() << endl;
    tit->TitleX = "True cluster energy [GeV]";
    tit->LogX = 1;
    b->minY = t->GetMinimum(t->GetListOfLeaves()->At(i)->GetName());
    b->maxY = t->GetMaximum(t->GetListOfLeaves()->At(i)->GetName());
    if (b->minY <= 0) b->minY += .00000000001;
    if (b->minY >= b->maxX) b->maxX = b->minY + .01;
    std::cout <<"min Y value: " << b->minY << endl;

    //Strange conditions I want to impose.
    if (!strcmp(t->GetListOfLeaves()->At(i)->GetName(), "cluster_nCells_tot")) b->maxY = 10;


    TH2D * temp = Plot_TH2D(t->GetLeaf("cluster_ENG_CALIB_TOT"), t->GetLeaf(t->GetListOfLeaves()->At(i)->GetName()), tit->TitleX, tit->TitleY, b, tit->LogX, tit->LogY);
    PNG_name = str2 + (TString)tit->SaveName + str1;
    gPad->Print(PNG_name);
    gPad->Clear();
  }



}

TitleAxis * GetTit(const char* branch_name){
  const char *tit = "unkonw";
  const char *name = "NameHEre";
  TitleAxis *names = new TitleAxis;

    if (!strcmp(branch_name, "runNumber")){
      tit = "Run Number";
      name = "runNumber";
      }
    else if (!strcmp(branch_name, "eventNumber")){
      tit = "Event Number";
      name = "eventNumber";
      }
    else if (!strcmp(branch_name, "truthE")){
      tit = "True pion energy [GeV]";
      name = "PionEnergy";
      }
    else if (!strcmp(branch_name, "truthPt")){
      tit = "True pion transverse momentum [GeV]";
      name = "PionPt";
      }
    else if (!strcmp(branch_name, "truthEta")){
      tit = "True pion pseudo-rapidity";
      name = "PionEta";
      }
    else if (!strcmp(branch_name, "truthPhi")){
      tit = "True pion azimuth angle";
      name = "PionPhi";
      }
    else if (!strcmp(branch_name, "truthPDG")){
      tit = "Particle code";
      name = "PionCode";
      }
    else if (!strcmp(branch_name, "nCluster")){
      tit = "Number of clusters/pion";
      name = "nClusters";
      }
    else if (!strcmp(branch_name, "clusterIndex")){
      tit = "Index of cluster for pion";
      name = "IndexPion";
      }
    else if (!strcmp(branch_name, "cluster_nCells")){
      tit = "Number of cells in cluster with E > 0";
      name = "nCellsNonZero";
      }
    else if (!strcmp(branch_name, "cluster_nCells_tot")){
      tit = "Total number of cells in cluster";
      name = "nCells";
      }
    else if (!strcmp(branch_name, "clusterIndex_1")){
      tit = "Index of cluster for pion";
      name = "IndexPion";
      }
    else if (!strcmp(branch_name, "clusterECalib")){
      tit = "Calibrated cluster energy [GeV]";
      name = "CalibEnergy";
      }
    else if (!strcmp(branch_name, "clusterPtCalib")){
      tit = "Calibrated cluster transverse momentum [GeV]";
      name = "CalibPt";
      }
    else if (!strcmp(branch_name, "clusterEtaCalib")){
      tit = "Cluster pseudorapidity after LCW calibration";
      name = "CalibEta";
      }
    else if (!strcmp(branch_name, "clusterPhiCalib")){
      tit = "Cluster azimuth after LCW calibration";
      name = "CalibPhi";
      }
    else if (!strcmp(branch_name, "cluster_sumCellECalib")){
      tit = "Sum of calibrated cell energies for E > 0 [GeV]";
      name = "CalibCellENonZero";
      }
    else if (!strcmp(branch_name, "clusterE")){
      tit = "Cluster EM scale energy [GeV]";
      name = "EMEnergyScale";
      names->LogY = 1;
      }
    else if (!strcmp(branch_name, "clusterPt")){
      tit = "Cluster EM scale transverse momentum [GeV]";
      name = "ClusterPt";
      }
    else if (!strcmp(branch_name, "clusterEta")){
      tit = "Cluster pseudorapidity";
      name = "ClusterEta";
      }
    else if (!strcmp(branch_name, "clusterPhi")){
      tit = "Cluster azimuth";
      name = "ClusterPhi";
      }
    else if (!strcmp(branch_name, "cluster_sumCellE")){
      tit = "Sum of EM scale energies for cells with E > 0";
      name = "sumCellE";
      }
    else if (!strcmp(branch_name, "cluster_EM_PROBABILITY")){
      tit = "Probability for cluster to be EM";
      name = "EMPro";
      }
    else if (!strcmp(branch_name, "cluster_HAD_WEIGHT")){
      tit = "Effective cluster signal weight after hcw";
      name = "HADWeight";
      }
    else if (!strcmp(branch_name, "cluster_OOC_WEIGHT")){
      tit = "Effecitive weight after ooc correction";
      name = "OOCWeight";
      }
    else if (!strcmp(branch_name, "cluster_DM_WEIGHT")){
      tit = "Effective cluster signal weight after dm";
      name = "DMWeight";
      }
    else if (!strcmp(branch_name, "cluster_ENG_CALIB_TOT")){
      tit = "True cluster energy [GeV]";
      name = "EngTrue";
      }
    else if (!strcmp(branch_name, "cluster_ENG_CALIB_OUT_T")){
      tit = "Energy depostited outside cluster but associated with it [GeV]";
      name = "EngOut";
      }
    else if (!strcmp(branch_name, "cluster_ENG_CALIB_DEAD_TOT")){
      tit = "Energy deposited in dead materiral";
      name = "EngDead";
      }
    else if (!strcmp(branch_name, "cluster_CENTER_MAG")){
      tit = "Distance From vertex [mm]";
      name = "CenterMag";
      }
    else if (!strcmp(branch_name, "cluster_FIRST_ENG_DENS")){
      tit = "EnergyDensity";
      name = "EngDen";
      }
    else if (!strcmp(branch_name, "cluster_FIRST_PHI")){
      tit = "Energy weighted azimuthal distance [rad]";
      name = "EngWeighPhi";
      }
    else if (!strcmp(branch_name, "cluster_FIRST_ETA")){
      tit = "Energy weighted pseudorapidity";
      name = "EngWeighEta";
      }
    else if (!strcmp(branch_name, "cluster_SECOND_R")){
      tit = "Second radial distnace of cells [mm#^2]";
      name = "SecR";
      }
    else if (!strcmp(branch_name, "cluster_SECOND_LAMBDA")){
      tit = "Second moment of distance of cells";
      name = "SecLambda";
      }
    else if (!strcmp(branch_name, "cluster_DELTA_PHI")){
      tit = "azimuthal distance of shower axis to nominal axis [rad]";
      name = "DeltaPhi [rad]";
      }
    else if (!strcmp(branch_name, "cluster_DELTA_THETA")){
      tit = "Polar distance between shower and nominal axis";
      name = "DelatTheta [rad]";
      }
    else if (!strcmp(branch_name, "cluster_DELTA_ALPHA")){
      tit = "Angular distance between shower and nominal axis";
      name = "DeltaAlpha [rad]";
      }
    else if (!strcmp(branch_name, "cluster_CENTER_X")){
      tit = "Cluster center of gavity x-coordinate";
      name = "CenterX [mm]";
      }
    else if (!strcmp(branch_name, "cluster_CENTER_Y")){
      tit = "Cluster center of gavity y-coordinate";
      name = "CenterY [mm]";
      }
    else if (!strcmp(branch_name, "cluster_CENTER_Z")){
      tit = "Cluster center of gavity z-coordinate";
      name = "CenterZ [mm]";
      }
    else if (!strcmp(branch_name, "cluster_CENTER_LAMBDA")){
      tit = "Distance from calorimeter front face [mm]";
      name = "CenterLambda";
      }
    else if (!strcmp(branch_name, "cluster_LATERAL")){
      tit = "Measure of lateral energy dispersion";
      name = "Lateral";
      names->LogY = 0;
      }
    else if (!strcmp(branch_name, "cluster_LONGITUDINAL")){
      tit = "Measure of longitudinal energy dispersion";
      name = "Longitudinal";
      names->LogY = 0;
      }
    else if (!strcmp(branch_name, "cluster_ENG_FRAC_EM")){
      tit = "Fraction of energy in EM calorimeter";
      name = "EngFrac";
      }
    else if (!strcmp(branch_name, "cluster_ENG_FRAC_MAX")){
      tit = "Most energitc cell signal fraction in cluster";
      name = "EngFracMax";
      }
    else if (!strcmp(branch_name, "cluster_ENG_FRAC_CORE")){
      tit = "Energy fraction in core of cluster";
      name = "EngFracCore";
      }
    else if (!strcmp(branch_name, "cluster_SECOND_ENG_DENS")){
      tit = "Energy weighted second moment of cell density [(#frac{GeV}{mm#^{3}})#^{2}]";
      name = "SecEngDensity";
      }
    else if (!strcmp(branch_name, "cluster_ISOLATION")){
      tit = "Measure for cluster isolation";
      name = "Isolation";
      names->LogY = 0;
      }
    else if (!strcmp(branch_name, "cluster_ENG_BAD_CELLS")){
      tit = "Energy in bad cells [GeV]";
      name = "EBadCells";
      }
    else if (!strcmp(branch_name, "cluster_N_BAD_CELLS")){
      tit = "Number of Bad Cells";
      name = "nBadCells";
      }
    else if (!strcmp(branch_name, "cluster_N_BAD_CELLS_CORR")){
      tit = "Bad number of cells correction";
      name = "nBadCellsCorr";
      }
    else if (!strcmp(branch_name, "cluster_BAD_CELLS_CORR_E")){
      tit = "Bad cells energy correction";
      name = "EBadcellsCorr";
      }
    else if(!strcmp(branch_name, "cluster_BADLARQ_FRAC")){
      tit = "Bad larq fraction";
      name = "BADLARQFrac";
      }
    else if (!strcmp(branch_name, "cluster_ENG_POS")){
      tit = "Sum of postive energy in cluster [GeV]";
      name = "EngPos";
      }
    else if (!strcmp(branch_name, "cluster_SIGNIFICANCE")){
      tit = "Cluster Significance";
      name = "Significance";
      names->LogY = 1;
      }
    else if (!strcmp(branch_name, "cluster_CELL_SIGNIFICANCE")){
      tit = "Max cell significance in cluster";
      names->LogY = 1;
      name = "MaxSignificance";
      }
    else if (!strcmp(branch_name, "cluster_CELL_SIG_SAMPLING")){
      tit = "Samping of ID at highest cell significance";
      name = "SigSampling";
      }
    else if (!strcmp(branch_name, "cluster_AVG_LAR_Q")){
      tit = "Average LAR Q";
      name = "AVGLARQ";
      }
    else if (!strcmp(branch_name, "cluster_AVG_TILE_Q")){
      tit = "Average Tile Q";
      name = "AVGTILEQ";
      }
    else if (!strcmp(branch_name, "cluster_ENG_BAD_HV_CELLS")){
      tit = "Energy of bad cells?";
      name = "EngBadHVCells";
      }
    else if (!strcmp(branch_name, "cluster_N_BAD_HV_CELLS")){
      tit = "Number of bad cells in HV?";
      name = "nBadHVCells";
      }
    else if (!strcmp(branch_name, "cluster_PTD")){
      tit = "Cluster longitudinal fragmentation";
      name = "PTD";
      }
    else if (!strcmp(branch_name, "cluster_MASS")){
      tit = "Cluster mass using only E>0";
      name = "Mass";
      }
    else if (!strcmp(branch_name, "EM_Shower")){
      tit = "EM shower";
      name = "EMShower";
      }
    else if (!strcmp(branch_name, "EM_Pro")){
      tit = "Probability of shower being EM";
      name = "EMProNetwork";
      }
    else if (!strcmp(branch_name, "CalibratedE")){
      tit = "Calibrated energy by network [GeV]";
      name = "NetworkCalibE";
      }
    else if (!strcmp(branch_name, "Delta_E")){
      tit = "Difference between calibrated energy and true energy [GeV]";
      name = "DeltaE";
      }
    else if (!strcmp(branch_name, "Delta_Calib_E")){
      tit = "Don;t know difference?";
      name = "DeltaDifferenceCalib";
      }
    else if (!strcmp(branch_name, "Delta_E")){
      tit = "Difference between calib and actual";
      name = "DeltaEA";
      }
    else if (!strcmp(branch_name, "Delta_Calib_E")){
      tit = "Repeat Value";
      name = "RepeatDelta";
      }
    else if (!strcmp(branch_name, "Ratio")){
      tit = "#E^{Calib}_{Cluster}/#E^{True}_{Cluster}";
      name = "Ratio";
      std::cout << "Running Branch Ratio." <<endl;
      }
    else if (!strcmp(branch_name, "Ratio_old")){
      tit = "#E^{LCW}_{Cluster}/#E^{True}_{Cluster}";
      name = "LCWRatio";
      std::cout << "Running Branch Ratio." <<endl;
      }

  names->TitleY = tit;
  names->SaveName = name;
  return names;
}



TGraph * IQR(TH2D * h){
  Double_t ColumnSum = 0, ColumnWeightedSum = 0;
  auto n_BinsX = h->GetNbinsX();
  auto n_BinsY = h->GetNbinsY();
  auto Xaxis = h->GetXaxis();
  auto Yaxis = h->GetYaxis();
  Double_t ArrX[n_BinsX];
  Double_t ArrY[n_BinsX];
  TH2D *subTH2D;
  Double_t x[n_BinsY];
  Double_t y[n_BinsY];
  Double_t xq[3] = {.25, .5, .75}, yq[3];
  Int_t bin;
  Double_t Median, q1, q3;
  Int_t Median_bin;
  Int_t cnt_t = 0, cnt_c;
  // string const &str;
  // Int_t n_NAN = 0;

  Int_t n_NAN = 1;
  for(Int_t i = 1; i <= n_BinsX; i++){
    cnt_t = 0;
    cnt_c = 0;
    for(Int_t j = 1; j <= n_BinsY; j++){
      bin = h->GetBin(i, j);
      x[j] = Yaxis->GetBinCenter(j);
      y[j] = h->GetBinContent(bin);
      cnt_t += h->GetBinContent(bin);
      // std::cout << "Bin: " << bin << "\tx: " << x[j] << "\ty: " << y[j] << endl;
    }
    q1 = 0;
    q3 = 0;
    Median = -1;
    for (Int_t j = 1; j <= n_BinsY; j++){
      cnt_c += y[j];
      if (cnt_c >= cnt_t/4 && !q1) {q1 = x[j];}
      if (cnt_c >= cnt_t/2 && Median == -1) {
        Median = TMath::Median(n_BinsY, &x[0], &y[0]);
        std::cout << "Median: " << Median << "\tQ2: " << x[j] << endl;
      }
      if (cnt_c >= 3*cnt_t/4 && !q3) {q3 = x[j];}

    }




    // std::stringstream s;
    // // TH2D * TH2Subset(TH2D *h,Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2){
    // subTH2D = TH2Subset(h, i, i+1, 0, n_BinsY);
    ArrX[i] = Xaxis->GetBinCenter(i);
    //
    // (h->ProjectionY("Name", i, i+1))->GetQuantiles(3, xq, yq);
    // h->ProjectionY("Name", i, i+1)->Draw("H");
    // if (i == 10) gPad->Print("name10.png");
    // else if (i == 20) gPad->Print("name20.png");
    // else if (i == 30) gPad->Print("name30.png");
    // else if (i == 40) gPad->Print("name40.png");
    // else if (i == 50) gPad->Print("name50.png");
    // else if (i == 60) gPad->Print("name60.png");
    // else if (i == 70) gPad->Print("name70.png");
    // else if (i == 80) gPad->Print("name80.png");
    // else if (i == 90) gPad->Print("name90.png");
    ArrY[i] = (q3 - q1) / TMath::Median(n_BinsY, &x[0], &y[0]);// - h->GetBinContent(h->FindBin(yq[2]))); // / TMath::Median(n_BinsY, &x[0], &y[0]);
    std::cout << "X Values: " << ArrX[i] << "\tY Values: " << ArrY[i] << "\tQuantile 1: " << q1<<"\tQuantile 3: " << q3<< "\tMedian: " << TMath::Median(n_BinsY, &x[0], &y[0]) << endl;
    if (Median == 0 || Xaxis->GetBinCenter(i) < .1) n_NAN++;
  }



  TGraph * g = new TGraph(n_BinsX-n_NAN, &ArrX[n_NAN], &ArrY[n_NAN]);

  g->SetLineColor(2);

  return g;




}
