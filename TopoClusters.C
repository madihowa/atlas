#define TopoClusters_cxx
#include "TopoClusters.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void TopoClusters::Loop()
{


  TFile *f = new TFile("results.root", "UPDATE");

//   In a ROOT session, you can do:
//      root> .L TopoClusters.C
//      root> TopoClusters t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results.root");
   // if (!f || !f->IsOpen()) {
   //    f = new TFile("results.root", "RECREATE");
   // }

   TTree *t_temp;
   f->GetObject("ClusterTree", t_temp);
   TTree *t_EM = new TTree("t_EM", "EM Tree");
   TTree *t_Had = new TTree("t_Had", "Hardonic shower Tree");

   Fill_Correlation();
   // Calculate_Correlation();
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   Init_EM(t_EM);
   Init_HAD(t_Had);
   TLeaf *L = fChain->GetLeaf("EM_Shower");
   Int_t val;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      val = L->GetValue(0);
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (ientry % 131072 == 0){
        printf("%lli\n", jentry);
      }
      if (val == 0){
        t_Had->Fill();
      }
      else if (val == 1){
        t_EM->Fill();
      }
      // else printf("Error Shower not EM or Had");
      // if (Cut(ientry) < 0) continue;
      std::cout << "Shower: " << val << endl;
   }
for (Int_t i =0; i < fChain->GetListOfBranches()->GetEntriesFast(); i++){
  printf("%s\t", titles[i]);
}
std:cout << endl;
for (Int_t j = 0; j < 20; j++){
  for (Int_t i =0; i < fChain->GetListOfBranches()->GetEntriesFast(); i++){
    std::cout << fChain->GetLeaf(titles[i])->GetValue(0) << "\t";
  }
  std::cout << endl;
}
   // His->Fill(x);
   printf("The total size is: %lld\n", nbytes);
   printf("The current working directory is: %s\n", gSystem->pwd());
   printf("Total Number of entries: %lli\n\n", nentries);
   //
   // Calculate_Correlation();
   // Corr->Print();

   // for (Int_t i = 0; i < 5; i++){
   //   for (Int_t j = 0; j < 5; j++){
   //     std::cout << (*Corr)[i][j] << "\t";
   //     Corr_Histogram->Fill(i, j, (*Corr)[i][j]);
   //   }
   //   printf("\n");
   // }
   //
   // // TFile *g = new TFile("Final_results.root", "RECREATE");
   // Histograms_List->Write("Histograms", TObject::kSingleKey);
   // TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
   // c1->Write("Here");
   // Corr_Histogram->Draw("colz");
   // // fChain->Write();
   // // g->Close();
   // t_Had->Print();
   // t_Had->AutoSave();
   // t_EM->Print();
   // t_EM->AutoSave();
   f->Write();
}
