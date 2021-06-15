#define ClusterTree_cxx
#include "ClusterTreeCat.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ClusterTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ClusterTree.C
//      root> ClusterTree t
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

   TFile *f = new TFile("results.root", "RECREATE");
   // TTree *t;
   // f->GetObject("ClusterTree", t);

   std::cout << "Current number of Entries: " << fChain->GetEntriesFast() << endl;

   TTree *EM_tree = new TTree("EM_tree", "The EM tree");
   TTree *Had_tree = new TTree("Had_tree", "The Hadronic tree");
   TTree *Orignal_tree = fChain->CloneTree();
   std::cout << "Tree Cloned" << endl;
   // TTree *tree;
   std::cout << "Current number of Entries: " << fChain->GetEntriesFast() << endl;
   // f->GetObject("ClusterTree", tree);
   Init_EM(EM_tree);
   Init_Had(Had_tree);
   std::cout << "Current number of Entries: " << fChain->GetEntriesFast() << endl;

   //
   // TBranch * Had_b_clusterECalib = Had_tree->Branch("clusterECalib", &clusterECalib, "clusterECalib/D");  //!
   // TBranch * Had_b_clusterPtCalib = Had_tree->Branch("clusterPtCalib", &clusterPtCalib, "clusterPtCalib/D");   //!
   // TBranch * Had_b_cluster_CENTER_LAMBDA = Had_tree->Branch("cluster_CENTER_LAMBDA", &cluster_CENTER_LAMBDA, "cluster_CENTER_LAMBDA/D");   //!
   // TBranch * Had_b_clusterEtaCalib = Had_tree->Branch("clusterEtaCalib", &clusterEtaCalib, "clusterEtaCalib/D");   //!
   // TBranch * Had_b_clusterPhiCalib = Had_tree->Branch("clusterPhiCalib", &clusterPhiCalib, "clusterPhiCalib/D");  //!
   // TBranch * Had_b_cluster_PTD = Had_tree->Branch("cluster_PTD", &cluster_PTD, "cluster_PTD/D");    //!
   // TBranch * Had_b_EM_Shower = Had_tree->Branch("EM_Shower", &EM_Shower, "EM_Shower/I");   //!
   // TBranch * Had_b_EM_Pro = Had_tree->Branch("EM_Pro", &EM_Pro, "EM_Pro/D");
   //


   std::cout << "Current number of Entries: " << fChain->GetEntriesFast() << endl;

   TLeaf* Leaf_Shower = fChain->GetLeaf("EM_Shower");

   if (fChain == 0) return;
   Long64_t r_EM = 0, r_Had = 0;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t n_EM = 0, n_Had = 0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Delta_E = cluster_ENG_CALIB_TOT - CalibratedE;
      Delta_Calib_E = cluster_ENG_CALIB_TOT - clusterECalib;
      // std::cout << "Current Value of the EM Shower: " <<EM_Shower << endl;
      // std::cout << "Current Value of the leaf EM Shower " << Leaf_Shower->GetValue(0) << endl;
      Orignal_tree->Fill();
      std::cout << "EM shower is: "<< EM_Shower << endl;
      if (EM_Shower >= .8){
        EM_tree->Fill();
        n_EM++;
        // std::cout << (EM_Pro) << endl;
        if (EM_Pro > .5) {r_EM++;}
      }
      else if (EM_Shower == 0)
      {
        Had_tree->Fill();
        n_Had++;
        if (EM_Pro < .5) {r_Had++;}
      }
      else std::cout << "Something is still wrong" << endl;
      // if (Cut(ientry) < 0) continue;

   }
   // f = new TFile("results.root", "RECREATE");
   // f->Write();
   EM_tree->Write();
   printf("Number of entries: %lld\nNumber of entries in hadronic tree: %lld\nNumber of entries in Electromagetic tree: %lld\n",nentries, n_Had, n_EM);
   std::ofstream ReportFile("../NetworkHistory.txt", ios::app);
   ReportFile << (Double_t)r_EM / (Double_t)n_EM << "\t|" << (Double_t)r_Had / (Double_t)n_Had << endl;
   std::cout << (r_EM) << endl;

   f->Close();
}
