//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 15 16:47:22 2020 by ROOT version 6.19/01
// from TTree TopoClusters/TopoClusters
// found on file: results.root
//////////////////////////////////////////////////////////

#ifndef TopoClusters_h
#define TopoClusters_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMatrix.h>
// Header file for the classes stored in the TTree if any.

class TopoClusters {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runNumber;
   Int_t           eventNumber;
   Float_t         truthE;
   Float_t         truthPt;
   Float_t         truthEta;
   Float_t         truthPhi;
   Int_t           truthPDG;
   Int_t           nCluster;
   Int_t           clusterIndex;
   Int_t           cluster_nCells;
   Int_t           cluster_nCells_tot;
   Float_t         clusterECalib;
   Float_t         clusterPtCalib;
   Float_t         clusterEtaCalib;
   Float_t         clusterPhiCalib;
   Float_t         cluster_sumCellECalib;
   Float_t         clusterE;
   Float_t         clusterPt;
   Float_t         clusterEta;
   Float_t         clusterPhi;
   Float_t         cluster_sumCellE;
   Float_t         cluster_EM_PROBABILITY;
   Float_t         cluster_HAD_WEIGHT;
   Float_t         cluster_OOC_WEIGHT;
   Float_t         cluster_DM_WEIGHT;
   Float_t         cluster_ENG_CALIB_TOT;
   Float_t         cluster_ENG_CALIB_OUT_T;
   Float_t         cluster_ENG_CALIB_DEAD_TOT;
   Float_t         cluster_CENTER_MAG;
   Float_t         cluster_FIRST_ENG_DENS;
   Float_t         cluster_FIRST_PHI;
   Float_t         cluster_FIRST_ETA;
   Float_t         cluster_SECOND_R;
   Float_t         cluster_SECOND_LAMBDA;
   Float_t         cluster_DELTA_PHI;
   Float_t         cluster_DELTA_THETA;
   Float_t         cluster_DELTA_ALPHA;
   Float_t         cluster_CENTER_X;
   Float_t         cluster_CENTER_Y;
   Float_t         cluster_CENTER_Z;
   Float_t         cluster_CENTER_LAMBDA;
   Float_t         cluster_LATERAL;
   Float_t         cluster_LONGITUDINAL;
   Float_t         cluster_ENG_FRAC_EM;
   Float_t         cluster_ENG_FRAC_MAX;
   Float_t         cluster_ENG_FRAC_CORE;
   Float_t         cluster_SECOND_ENG_DENS;
   Float_t         cluster_ISOLATION;
   Float_t         cluster_ENG_BAD_CELLS;
   Float_t         cluster_N_BAD_CELLS;
   Float_t         cluster_N_BAD_CELLS_CORR;
   Float_t         cluster_BAD_CELLS_CORR_E;
   Float_t         cluster_BADLARQ_FRAC;
   Float_t         cluster_ENG_POS;
   Float_t         cluster_SIGNIFICANCE;
   Float_t         cluster_CELL_SIGNIFICANCE;
   Float_t         cluster_CELL_SIG_SAMPLING;
   Float_t         cluster_AVG_LAR_Q;
   Float_t         cluster_AVG_TILE_Q;
   Float_t         cluster_ENG_BAD_HV_CELLS;
   Float_t         cluster_N_BAD_HV_CELLS;
   Float_t         cluster_PTD;
   Float_t         cluster_MASS;
   Float_t         EM_Pro;
   Int_t           EM_Shower; //1 is em shower 0 is hadronic

   Double_t        TrueE;
   Double_t        CalibratedE;
   TH1D            *His;
   TH2D            *Corr_Histogram;
   TList           *Histograms_List;
   TMatrix         *Corr;
   TMatrix         Temp_Matrix;
   TMatrix         Sum_Matrix;
   TMatrix         Product_Matrix;
   Int_t           n_branches;
   const char      *titles[124];
   // List of branches




      TBranch        *EM_b_runNumber;   //!
      TBranch        *EM_b_eventNumber;   //!
      TBranch        *EM_b_truthE;   //!
      TBranch        *EM_b_truthPt;   //!
      TBranch        *EM_b_truthEta;   //!
      TBranch        *EM_b_truthPhi;   //!
      TBranch        *EM_b_truthPDG;   //!
      TBranch        *EM_b_nCluster;   //!
      TBranch        *EM_b_cluster_nCells;   //!
      TBranch        *EM_b_cluster_nCells_tot;   //!
      TBranch        *EM_b_clusterIndex;   //!
      TBranch        *EM_b_clusterECalib;   //!
      TBranch        *EM_b_clusterPtCalib;   //!
      TBranch        *EM_b_clusterEtaCalib;   //!
      TBranch        *EM_b_clusterPhiCalib;      // t_Had->Print();
      // t_Had->AutoSave();
      // t_EM->Print();
      // t_EM->AutoSave();//!
      TBranch        *EM_b_cluster_sumCellECAlib;   //!
      TBranch        *EM_b_clusterE;   //!
      TBranch        *EM_b_clusterPt;   //!
      TBranch        *EM_b_clusterEta;   //!
      TBranch        *EM_b_clusterPhi;   //!
      TBranch        *EM_b_cluster_sumCellE;   //!
      TBranch        *EM_b_cluster_EM_PROBABILITY;   //!
      TBranch        *EM_b_cluster_HAD_WEIGHT;   //!
      TBranch        *EM_b_cluster_OOC_WEIGHT;   //!
      TBranch        *EM_b_cluster_DM_WEIGHT;   //!
      TBranch        *EM_b_cluster_ENG_CALIB_TOT;   //!
      TBranch        *EM_b_cluster_ENG_CALIB_OUT_T;   //!
      TBranch        *EM_b_cluster_ENG_CALIB_DEAD_TOT;   //!
      TBranch        *EM_b_cluster_CENTER_MAG;   //!
      TBranch        *EM_b_cluster_FIRST_ENG_DENS;   //!
      TBranch        *EM_b_cluster_FIRST_PHI;   //!
      TBranch        *EM_b_cluster_FIRST_ETA;   //!
      TBranch        *EM_b_cluster_SECOND_R;   //!
      TBranch        *EM_b_cluster_SECOND_LAMBDA;   //!
      TBranch        *EM_b_cluster_DELTA_PHI;   //!
      TBranch        *EM_b_cluster_DELTA_THETA;   //!
      TBranch        *EM_b_cluster_DELTA_ALPHA;   //!
      TBranch        *EM_b_cluster_CENTER_X;   //!
      TBranch        *EM_b_cluster_CENTER_Y;   //!
      TBranch        *EM_b_cluster_CENTER_Z;   //!
      TBranch        *EM_b_cluster_CENTER_LAMBDA;   //!
      TBranch        *EM_b_cluster_LATERAL;   //!
      TBranch        *EM_b_cluster_LONGITUDINAL;   //!
      TBranch        *EM_b_cluster_ENG_FRAC_EM;   //!
      TBranch        *EM_b_cluster_ENG_FRAC_MAX;   //!
      TBranch        *EM_b_cluster_ENG_FRAC_CORE;   //!
      TBranch        *EM_b_cluster_SECOND_ENG_DENS;   //!
      TBranch        *EM_b_cluster_ISOLATION;   //!
      TBranch        *EM_b_cluster_ENG_BAD_CELLS;   //!
      TBranch        *EM_b_cluster_N_BAD_CELLS;   //!
      TBranch        *EM_b_cluster_N_BAD_CELLS_CORR;   //!
      TBranch        *EM_b_cluster_BAD_CELLS_CORR_E;   //!
      TBranch        *EM_b_cluster_BADLARQ_FRAC;   //!
      TBranch        *EM_b_cluster_ENG_POS;   //!
      TBranch        *EM_b_cluster_SIGNIFICANCE;   //!
      TBranch        *EM_b_cluster_CELL_SIGNIFICANCE;   //!
      TBranch        *EM_b_cluster_CELL_SIG_SAMPLING;   //!
      TBranch        *EM_b_cluster_AVG_LAR_Q;   //!
      TBranch        *EM_b_cluster_AVG_TILE_Q;   //!
      TBranch        *EM_b_cluster_ENG_BAD_HV_CELLS;   //!
      TBranch        *EM_b_cluster_N_BAD_HV_CELLS;   //!
      TBranch        *EM_b_cluster_PTD;   //!
      TBranch        *EM_b_cluster_MASS;   //!
      TBranch        *EM_b_TrueE;   //!
      TBranch        *EM_b_CalibratedE;   //!
      // TBranch        *EM_b_His; //!
      // TBranch        *EM_b_Histograms_List; //!
      TBranch        *EM_b_Corr;
      TBranch        *EM_b_EM_Pro;
      TBranch        *EM_b_EM_Shower;




      TBranch        *Had_b_runNumber;   //!
      TBranch        *Had_b_eventNumber;   //!
      TBranch        *Had_b_truthE;   //!
      TBranch        *Had_b_truthPt;   //!
      TBranch        *Had_b_truthEta;   //!
      TBranch        *Had_b_truthPhi;   //!
      TBranch        *Had_b_truthPDG;   //!
      TBranch        *Had_b_nCluster;   //!
      TBranch        *Had_b_cluster_nCells;   //!
      TBranch        *Had_b_cluster_nCells_tot;   //!
      TBranch        *Had_b_clusterIndex;   //!
      TBranch        *Had_b_clusterECalib;   //!
      TBranch        *Had_b_clusterPtCalib;   //!
      TBranch        *Had_b_clusterEtaCalib;   //!
      TBranch        *Had_b_clusterPhiCalib;      // t_Had->Print();
      // t_Had->AutoSave();
      // t_EM->Print();
      // t_EM->AutoSave();//!
      TBranch        *Had_b_cluster_sumCellECAlib;   //!
      TBranch        *Had_b_clusterE;   //!
      TBranch        *Had_b_clusterPt;   //!
      TBranch        *Had_b_clusterEta;   //!
      TBranch        *Had_b_clusterPhi;   //!
      TBranch        *Had_b_cluster_sumCellE;   //!
      TBranch        *Had_b_cluster_EM_PROBABILITY;   //!
      TBranch        *Had_b_cluster_HAD_WEIGHT;   //!
      TBranch        *Had_b_cluster_OOC_WEIGHT;   //!
      TBranch        *Had_b_cluster_DM_WEIGHT;   //!
      TBranch        *Had_b_cluster_ENG_CALIB_TOT;   //!
      TBranch        *Had_b_cluster_ENG_CALIB_OUT_T;   //!
      TBranch        *Had_b_cluster_ENG_CALIB_DEAD_TOT;   //!
      TBranch        *Had_b_cluster_CENTER_MAG;   //!
      TBranch        *Had_b_cluster_FIRST_ENG_DENS;   //!
      TBranch        *Had_b_cluster_FIRST_PHI;   //!
      TBranch        *Had_b_cluster_FIRST_ETA;   //!
      TBranch        *Had_b_cluster_SECOND_R;   //!
      TBranch        *Had_b_cluster_SECOND_LAMBDA;   //!
      TBranch        *Had_b_cluster_DELTA_PHI;   //!
      TBranch        *Had_b_cluster_DELTA_THETA;   //!
      TBranch        *Had_b_cluster_DELTA_ALPHA;   //!
      TBranch        *Had_b_cluster_CENTER_X;   //!
      TBranch        *Had_b_cluster_CENTER_Y;   //!
      TBranch        *Had_b_cluster_CENTER_Z;   //!
      TBranch        *Had_b_cluster_CENTER_LAMBDA;   //!
      TBranch        *Had_b_cluster_LATERAL;   //!
      TBranch        *Had_b_cluster_LONGITUDINAL;   //!
      TBranch        *Had_b_cluster_ENG_FRAC_EM;   //!
      TBranch        *Had_b_cluster_ENG_FRAC_MAX;   //!
      TBranch        *Had_b_cluster_ENG_FRAC_CORE;   //!
      TBranch        *Had_b_cluster_SECOND_ENG_DENS;   //!
      TBranch        *Had_b_cluster_ISOLATION;   //!
      TBranch        *Had_b_cluster_ENG_BAD_CELLS;   //!
      TBranch        *Had_b_cluster_N_BAD_CELLS;   //!
      TBranch        *Had_b_cluster_N_BAD_CELLS_CORR;   //!
      TBranch        *Had_b_cluster_BAD_CELLS_CORR_E;   //!
      TBranch        *Had_b_cluster_BADLARQ_FRAC;   //!
      TBranch        *Had_b_cluster_ENG_POS;   //!
      TBranch        *Had_b_cluster_SIGNIFICANCE;   //!
      TBranch        *Had_b_cluster_CELL_SIGNIFICANCE;   //!
      TBranch        *Had_b_cluster_CELL_SIG_SAMPLING;   //!
      TBranch        *Had_b_cluster_AVG_LAR_Q;   //!
      TBranch        *Had_b_cluster_AVG_TILE_Q;   //!
      TBranch        *Had_b_cluster_ENG_BAD_HV_CELLS;   //!
      TBranch        *Had_b_cluster_N_BAD_HV_CELLS;   //!
      TBranch        *Had_b_cluster_PTD;   //!
      TBranch        *Had_b_cluster_MASS;   //!
      TBranch        *Had_b_TrueE;   //!
      TBranch        *Had_b_CalibratedE;   //!
      // TBranch        *Had_b_His; //!
      // TBranch        *Had_b_Histograms_List; //!
      TBranch        *Had_b_Corr;
      TBranch        *Had_b_EM_Pro;
      TBranch        *Had_b_EM_Shower;


   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_truthE;   //!
   TBranch        *b_truthPt;   //!
   TBranch        *b_truthEta;   //!
   TBranch        *b_truthPhi;   //!
   TBranch        *b_truthPDG;   //!
   TBranch        *b_nCluster;   //!
   TBranch        *b_cluster_nCells;   //!
   TBranch        *b_cluster_nCells_tot;   //!
   TBranch        *b_clusterIndex;   //!
   TBranch        *b_clusterECalib;   //!
   TBranch        *b_clusterPtCalib;   //!
   TBranch        *b_clusterEtaCalib;   //!
   TBranch        *b_clusterPhiCalib;      // t_Had->Print();
   // t_Had->AutoSave();
   // t_EM->Print();
   // t_EM->AutoSave();//!
   TBranch        *b_cluster_sumCellECAlib;   //!
   TBranch        *b_clusterE;   //!
   TBranch        *b_clusterPt;   //!
   TBranch        *b_clusterEta;   //!
   TBranch        *b_clusterPhi;   //!
   TBranch        *b_cluster_sumCellE;   //!
   TBranch        *b_cluster_EM_PROBABILITY;   //!
   TBranch        *b_cluster_HAD_WEIGHT;   //!
   TBranch        *b_cluster_OOC_WEIGHT;   //!
   TBranch        *b_cluster_DM_WEIGHT;   //!
   TBranch        *b_cluster_ENG_CALIB_TOT;   //!
   TBranch        *b_cluster_ENG_CALIB_OUT_T;   //!
   TBranch        *b_cluster_ENG_CALIB_DEAD_TOT;   //!
   TBranch        *b_cluster_CENTER_MAG;   //!
   TBranch        *b_cluster_FIRST_ENG_DENS;   //!
   TBranch        *b_cluster_FIRST_PHI;   //!
   TBranch        *b_cluster_FIRST_ETA;   //!
   TBranch        *b_cluster_SECOND_R;   //!
   TBranch        *b_cluster_SECOND_LAMBDA;   //!
   TBranch        *b_cluster_DELTA_PHI;   //!
   TBranch        *b_cluster_DELTA_THETA;   //!
   TBranch        *b_cluster_DELTA_ALPHA;   //!
   TBranch        *b_cluster_CENTER_X;   //!
   TBranch        *b_cluster_CENTER_Y;   //!
   TBranch        *b_cluster_CENTER_Z;   //!
   TBranch        *b_cluster_CENTER_LAMBDA;   //!
   TBranch        *b_cluster_LATERAL;   //!
   TBranch        *b_cluster_LONGITUDINAL;   //!
   TBranch        *b_cluster_ENG_FRAC_EM;   //!
   TBranch        *b_cluster_ENG_FRAC_MAX;   //!
   TBranch        *b_cluster_ENG_FRAC_CORE;   //!
   TBranch        *b_cluster_SECOND_ENG_DENS;   //!
   TBranch        *b_cluster_ISOLATION;   //!
   TBranch        *b_cluster_ENG_BAD_CELLS;   //!
   TBranch        *b_cluster_N_BAD_CELLS;   //!
   TBranch        *b_cluster_N_BAD_CELLS_CORR;   //!
   TBranch        *b_cluster_BAD_CELLS_CORR_E;   //!
   TBranch        *b_cluster_BADLARQ_FRAC;   //!
   TBranch        *b_cluster_ENG_POS;   //!
   TBranch        *b_cluster_SIGNIFICANCE;   //!
   TBranch        *b_cluster_CELL_SIGNIFICANCE;   //!
   TBranch        *b_cluster_CELL_SIG_SAMPLING;   //!
   TBranch        *b_cluster_AVG_LAR_Q;   //!
   TBranch        *b_cluster_AVG_TILE_Q;   //!
   TBranch        *b_cluster_ENG_BAD_HV_CELLS;   //!
   TBranch        *b_cluster_N_BAD_HV_CELLS;   //!
   TBranch        *b_cluster_PTD;   //!
   TBranch        *b_cluster_MASS;   //!
   TBranch        *b_TrueE;   //!
   TBranch        *b_CalibratedE;   //!
   // TBranch        *b_His; //!
   // TBranch        *b_Histograms_List; //!
   TBranch        *b_Corr;
   TBranch        *b_EM_Pro;
   TBranch        *b_EM_Shower;

   TopoClusters(TTree *tree=0);
   virtual ~TopoClusters();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Init_EM(TTree *tree);
   virtual void     Init_HAD(TTree *tree);
   virtual void     InitHist(TTree *tree);
   virtual void     Init_Matrix();
   virtual void     Calculate_Correlation();
   virtual void     Fill_Correlation();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TopoClusters_cxx
TopoClusters::TopoClusters(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("results.root");
      }
      f->GetObject("ClusterTree",tree);
    }

   Init(tree);
   InitHist(tree);
   std::cout << tree << endl;
   n_branches = tree->GetListOfBranches()->GetEntriesFast();
   Init_Matrix();
   for(Int_t i = 0; i < n_branches; i++){
     titles[i] = tree->GetListOfLeaves()->operator[](i)->GetName();
   }
}
TopoClusters::~TopoClusters()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TopoClusters::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TopoClusters::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TopoClusters::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("truthE", &truthE, &b_truthE);
   fChain->SetBranchAddress("truthPt", &truthPt, &b_truthPt);
   fChain->SetBranchAddress("truthEta", &truthEta, &b_truthEta);
   fChain->SetBranchAddress("truthPhi", &truthPhi, &b_truthPhi);
   fChain->SetBranchAddress("truthPDG", &truthPDG, &b_truthPDG);
   fChain->SetBranchAddress("nCluster", &nCluster, &b_nCluster);
   fChain->SetBranchAddress("clusterIndex", &clusterIndex, &b_clusterIndex);
   fChain->SetBranchAddress("cluster_nCells", &cluster_nCells, &b_cluster_nCells);
   fChain->SetBranchAddress("cluster_nCells_tot", &cluster_nCells_tot, &b_cluster_nCells_tot);
//    fChain->SetBranchAddress("clusterIndex", &clusterIndex, &b_clusterIndex);
   fChain->SetBranchAddress("clusterECalib", &clusterECalib, &b_clusterECalib);
   fChain->SetBranchAddress("clusterPtCalib", &clusterPtCalib, &b_clusterPtCalib);
   fChain->SetBranchAddress("clusterEtaCalib", &clusterEtaCalib, &b_clusterEtaCalib);
   fChain->SetBranchAddress("clusterPhiCalib", &clusterPhiCalib, &b_clusterPhiCalib);
   fChain->SetBranchAddress("cluster_sumCellECalib", &cluster_sumCellECalib, &b_cluster_sumCellECAlib);
   fChain->SetBranchAddress("clusterE", &clusterE, &b_clusterE);
   fChain->SetBranchAddress("clusterPt", &clusterPt, &b_clusterPt);
   fChain->SetBranchAddress("clusterEta", &clusterEta, &b_clusterEta);
   fChain->SetBranchAddress("clusterPhi", &clusterPhi, &b_clusterPhi);
   fChain->SetBranchAddress("cluster_sumCellE", &cluster_sumCellE, &b_cluster_sumCellE);
   fChain->SetBranchAddress("cluster_EM_PROBABILITY", &cluster_EM_PROBABILITY, &b_cluster_EM_PROBABILITY);
   fChain->SetBranchAddress("cluster_HAD_WEIGHT", &cluster_HAD_WEIGHT, &b_cluster_HAD_WEIGHT);
   fChain->SetBranchAddress("cluster_OOC_WEIGHT", &cluster_OOC_WEIGHT, &b_cluster_OOC_WEIGHT);
   fChain->SetBranchAddress("cluster_DM_WEIGHT", &cluster_DM_WEIGHT, &b_cluster_DM_WEIGHT);
   fChain->SetBranchAddress("cluster_ENG_CALIB_TOT", &cluster_ENG_CALIB_TOT, &b_cluster_ENG_CALIB_TOT);
   fChain->SetBranchAddress("cluster_ENG_CALIB_OUT_T", &cluster_ENG_CALIB_OUT_T, &b_cluster_ENG_CALIB_OUT_T);
   fChain->SetBranchAddress("cluster_ENG_CALIB_DEAD_TOT", &cluster_ENG_CALIB_DEAD_TOT, &b_cluster_ENG_CALIB_DEAD_TOT);
   fChain->SetBranchAddress("cluster_CENTER_MAG", &cluster_CENTER_MAG, &b_cluster_CENTER_MAG);
   fChain->SetBranchAddress("cluster_FIRST_ENG_DENS", &cluster_FIRST_ENG_DENS, &b_cluster_FIRST_ENG_DENS);
   fChain->SetBranchAddress("cluster_FIRST_PHI", &cluster_FIRST_PHI, &b_cluster_FIRST_PHI);
   fChain->SetBranchAddress("cluster_FIRST_ETA", &cluster_FIRST_ETA, &b_cluster_FIRST_ETA);
   fChain->SetBranchAddress("cluster_SECOND_R", &cluster_SECOND_R, &b_cluster_SECOND_R);
   fChain->SetBranchAddress("cluster_SECOND_LAMBDA", &cluster_SECOND_LAMBDA, &b_cluster_SECOND_LAMBDA);
   fChain->SetBranchAddress("cluster_DELTA_PHI", &cluster_DELTA_PHI, &b_cluster_DELTA_PHI);
   fChain->SetBranchAddress("cluster_DELTA_THETA", &cluster_DELTA_THETA, &b_cluster_DELTA_THETA);
   fChain->SetBranchAddress("cluster_DELTA_ALPHA", &cluster_DELTA_ALPHA, &b_cluster_DELTA_ALPHA);
   fChain->SetBranchAddress("cluster_CENTER_X", &cluster_CENTER_X, &b_cluster_CENTER_X);
   fChain->SetBranchAddress("cluster_CENTER_Y", &cluster_CENTER_Y, &b_cluster_CENTER_Y);
   fChain->SetBranchAddress("cluster_CENTER_Z", &cluster_CENTER_Z, &b_cluster_CENTER_Z);
   fChain->SetBranchAddress("cluster_CENTER_LAMBDA", &cluster_CENTER_LAMBDA, &b_cluster_CENTER_LAMBDA);
   fChain->SetBranchAddress("cluster_LATERAL", &cluster_LATERAL, &b_cluster_LATERAL);
   fChain->SetBranchAddress("cluster_LONGITUDINAL", &cluster_LONGITUDINAL, &b_cluster_LONGITUDINAL);
   fChain->SetBranchAddress("cluster_ENG_FRAC_EM", &cluster_ENG_FRAC_EM, &b_cluster_ENG_FRAC_EM);
   fChain->SetBranchAddress("cluster_ENG_FRAC_MAX", &cluster_ENG_FRAC_MAX, &b_cluster_ENG_FRAC_MAX);
   fChain->SetBranchAddress("cluster_ENG_FRAC_CORE", &cluster_ENG_FRAC_CORE, &b_cluster_ENG_FRAC_CORE);
   fChain->SetBranchAddress("cluster_SECOND_ENG_DENS", &cluster_SECOND_ENG_DENS, &b_cluster_SECOND_ENG_DENS);
   fChain->SetBranchAddress("cluster_ISOLATION", &cluster_ISOLATION, &b_cluster_ISOLATION);
   fChain->SetBranchAddress("cluster_ENG_BAD_CELLS", &cluster_ENG_BAD_CELLS, &b_cluster_ENG_BAD_CELLS);
   fChain->SetBranchAddress("cluster_N_BAD_CELLS", &cluster_N_BAD_CELLS, &b_cluster_N_BAD_CELLS);
   fChain->SetBranchAddress("cluster_N_BAD_CELLS_CORR", &cluster_N_BAD_CELLS_CORR, &b_cluster_N_BAD_CELLS_CORR);
   fChain->SetBranchAddress("cluster_BAD_CELLS_CORR_E", &cluster_BAD_CELLS_CORR_E, &b_cluster_BAD_CELLS_CORR_E);
   fChain->SetBranchAddress("cluster_BADLARQ_FRAC", &cluster_BADLARQ_FRAC, &b_cluster_BADLARQ_FRAC);
   fChain->SetBranchAddress("cluster_ENG_POS", &cluster_ENG_POS, &b_cluster_ENG_POS);
   fChain->SetBranchAddress("cluster_SIGNIFICANCE", &cluster_SIGNIFICANCE, &b_cluster_SIGNIFICANCE);
   fChain->SetBranchAddress("cluster_CELL_SIGNIFICANCE", &cluster_CELL_SIGNIFICANCE, &b_cluster_CELL_SIGNIFICANCE);
   fChain->SetBranchAddress("cluster_CELL_SIG_SAMPLING", &cluster_CELL_SIG_SAMPLING, &b_cluster_CELL_SIG_SAMPLING);
   fChain->SetBranchAddress("cluster_AVG_LAR_Q", &cluster_AVG_LAR_Q, &b_cluster_AVG_LAR_Q);
   fChain->SetBranchAddress("cluster_AVG_TILE_Q", &cluster_AVG_TILE_Q, &b_cluster_AVG_TILE_Q);
   fChain->SetBranchAddress("cluster_ENG_BAD_HV_CELLS", &cluster_ENG_BAD_HV_CELLS, &b_cluster_ENG_BAD_HV_CELLS);
   fChain->SetBranchAddress("cluster_N_BAD_HV_CELLS", &cluster_N_BAD_HV_CELLS, &b_cluster_N_BAD_HV_CELLS);
   fChain->SetBranchAddress("cluster_PTD", &cluster_PTD, &b_cluster_PTD);
   fChain->SetBranchAddress("cluster_MASS", &cluster_MASS, &b_cluster_MASS);


   fChain->SetBranchAddress("EM_Shower", &EM_Shower, &b_EM_Shower);
   fChain->SetBranchAddress("EM_Pro", &EM_Pro, &b_EM_Pro);

   fChain->SetBranchAddress("TrueE", &TrueE, &b_TrueE);
   fChain->SetBranchAddress("CalibratedE", &CalibratedE, &b_CalibratedE);
   fChain->SetBranchAddress("Corr", &Corr, &b_Corr);
   Notify();
   // InitHist(tree);
}
void TopoClusters::Init_EM(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fCurrent = -1;
   tree->SetMakeClass(1);

   tree->SetBranchAddress("runNumber", &runNumber, &EM_b_runNumber);
   tree->SetBranchAddress("eventNumber", &eventNumber, &EM_b_eventNumber);
   tree->SetBranchAddress("truthE", &truthE, &EM_b_truthE);
   tree->SetBranchAddress("truthPt", &truthPt, &EM_b_truthPt);
   tree->SetBranchAddress("truthEta", &truthEta, &EM_b_truthEta);
   tree->SetBranchAddress("truthPhi", &truthPhi, &EM_b_truthPhi);
   tree->SetBranchAddress("truthPDG", &truthPDG, &EM_b_truthPDG);
   tree->SetBranchAddress("nCluster", &nCluster, &EM_b_nCluster);
   tree->SetBranchAddress("clusterIndex", &clusterIndex, &EM_b_clusterIndex);
   tree->SetBranchAddress("cluster_nCells", &cluster_nCells, &EM_b_cluster_nCells);
   tree->SetBranchAddress("cluster_nCells_tot", &cluster_nCells_tot, &EM_b_cluster_nCells_tot);
//    tree->SetBranchAddress("clusterIndex", &clusterIndex, &EM_b_clusterIndex);
   tree->SetBranchAddress("clusterECalib", &clusterECalib, &EM_b_clusterECalib);
   tree->SetBranchAddress("clusterPtCalib", &clusterPtCalib, &EM_b_clusterPtCalib);
   tree->SetBranchAddress("clusterEtaCalib", &clusterEtaCalib, &EM_b_clusterEtaCalib);
   tree->SetBranchAddress("clusterPhiCalib", &clusterPhiCalib, &EM_b_clusterPhiCalib);
   tree->SetBranchAddress("cluster_sumCellECalib", &cluster_sumCellECalib, &EM_b_cluster_sumCellECAlib);
   tree->SetBranchAddress("clusterE", &clusterE, &EM_b_clusterE);
   tree->SetBranchAddress("clusterPt", &clusterPt, &EM_b_clusterPt);
   tree->SetBranchAddress("clusterEta", &clusterEta, &EM_b_clusterEta);
   tree->SetBranchAddress("clusterPhi", &clusterPhi, &EM_b_clusterPhi);
   tree->SetBranchAddress("cluster_sumCellE", &cluster_sumCellE, &EM_b_cluster_sumCellE);
   tree->SetBranchAddress("cluster_EM_PROBABILITY", &cluster_EM_PROBABILITY, &EM_b_cluster_EM_PROBABILITY);
   tree->SetBranchAddress("cluster_HAD_WEIGHT", &cluster_HAD_WEIGHT, &EM_b_cluster_HAD_WEIGHT);
   tree->SetBranchAddress("cluster_OOC_WEIGHT", &cluster_OOC_WEIGHT, &EM_b_cluster_OOC_WEIGHT);
   tree->SetBranchAddress("cluster_DM_WEIGHT", &cluster_DM_WEIGHT, &EM_b_cluster_DM_WEIGHT);
   tree->SetBranchAddress("cluster_ENG_CALIB_TOT", &cluster_ENG_CALIB_TOT, &EM_b_cluster_ENG_CALIB_TOT);
   tree->SetBranchAddress("cluster_ENG_CALIB_OUT_T", &cluster_ENG_CALIB_OUT_T, &EM_b_cluster_ENG_CALIB_OUT_T);
   tree->SetBranchAddress("cluster_ENG_CALIB_DEAD_TOT", &cluster_ENG_CALIB_DEAD_TOT, &EM_b_cluster_ENG_CALIB_DEAD_TOT);
   tree->SetBranchAddress("cluster_CENTER_MAG", &cluster_CENTER_MAG, &EM_b_cluster_CENTER_MAG);
   tree->SetBranchAddress("cluster_FIRST_ENG_DENS", &cluster_FIRST_ENG_DENS, &EM_b_cluster_FIRST_ENG_DENS);
   tree->SetBranchAddress("cluster_FIRST_PHI", &cluster_FIRST_PHI, &EM_b_cluster_FIRST_PHI);
   tree->SetBranchAddress("cluster_FIRST_ETA", &cluster_FIRST_ETA, &EM_b_cluster_FIRST_ETA);
   tree->SetBranchAddress("cluster_SECOND_R", &cluster_SECOND_R, &EM_b_cluster_SECOND_R);
   tree->SetBranchAddress("cluster_SECOND_LAMBDA", &cluster_SECOND_LAMBDA, &EM_b_cluster_SECOND_LAMBDA);
   tree->SetBranchAddress("cluster_DELTA_PHI", &cluster_DELTA_PHI, &EM_b_cluster_DELTA_PHI);
   tree->SetBranchAddress("cluster_DELTA_THETA", &cluster_DELTA_THETA, &EM_b_cluster_DELTA_THETA);
   tree->SetBranchAddress("cluster_DELTA_ALPHA", &cluster_DELTA_ALPHA, &EM_b_cluster_DELTA_ALPHA);
   tree->SetBranchAddress("cluster_CENTER_X", &cluster_CENTER_X, &EM_b_cluster_CENTER_X);
   tree->SetBranchAddress("cluster_CENTER_Y", &cluster_CENTER_Y, &EM_b_cluster_CENTER_Y);
   tree->SetBranchAddress("cluster_CENTER_Z", &cluster_CENTER_Z, &EM_b_cluster_CENTER_Z);
   tree->SetBranchAddress("cluster_CENTER_LAMBDA", &cluster_CENTER_LAMBDA, &EM_b_cluster_CENTER_LAMBDA);
   tree->SetBranchAddress("cluster_LATERAL", &cluster_LATERAL, &EM_b_cluster_LATERAL);
   tree->SetBranchAddress("cluster_LONGITUDINAL", &cluster_LONGITUDINAL, &EM_b_cluster_LONGITUDINAL);
   tree->SetBranchAddress("cluster_ENG_FRAC_EM", &cluster_ENG_FRAC_EM, &EM_b_cluster_ENG_FRAC_EM);
   tree->SetBranchAddress("cluster_ENG_FRAC_MAX", &cluster_ENG_FRAC_MAX, &EM_b_cluster_ENG_FRAC_MAX);
   tree->SetBranchAddress("cluster_ENG_FRAC_CORE", &cluster_ENG_FRAC_CORE, &EM_b_cluster_ENG_FRAC_CORE);
   tree->SetBranchAddress("cluster_SECOND_ENG_DENS", &cluster_SECOND_ENG_DENS, &EM_b_cluster_SECOND_ENG_DENS);
   tree->SetBranchAddress("cluster_ISOLATION", &cluster_ISOLATION, &EM_b_cluster_ISOLATION);
   tree->SetBranchAddress("cluster_ENG_BAD_CELLS", &cluster_ENG_BAD_CELLS, &EM_b_cluster_ENG_BAD_CELLS);
   tree->SetBranchAddress("cluster_N_BAD_CELLS", &cluster_N_BAD_CELLS, &EM_b_cluster_N_BAD_CELLS);
   tree->SetBranchAddress("cluster_N_BAD_CELLS_CORR", &cluster_N_BAD_CELLS_CORR, &EM_b_cluster_N_BAD_CELLS_CORR);
   tree->SetBranchAddress("cluster_BAD_CELLS_CORR_E", &cluster_BAD_CELLS_CORR_E, &EM_b_cluster_BAD_CELLS_CORR_E);
   tree->SetBranchAddress("cluster_BADLARQ_FRAC", &cluster_BADLARQ_FRAC, &EM_b_cluster_BADLARQ_FRAC);
   tree->SetBranchAddress("cluster_ENG_POS", &cluster_ENG_POS, &EM_b_cluster_ENG_POS);
   tree->SetBranchAddress("cluster_SIGNIFICANCE", &cluster_SIGNIFICANCE, &EM_b_cluster_SIGNIFICANCE);
   tree->SetBranchAddress("cluster_CELL_SIGNIFICANCE", &cluster_CELL_SIGNIFICANCE, &EM_b_cluster_CELL_SIGNIFICANCE);
   tree->SetBranchAddress("cluster_CELL_SIG_SAMPLING", &cluster_CELL_SIG_SAMPLING, &EM_b_cluster_CELL_SIG_SAMPLING);
   tree->SetBranchAddress("cluster_AVG_LAR_Q", &cluster_AVG_LAR_Q, &EM_b_cluster_AVG_LAR_Q);
   tree->SetBranchAddress("cluster_AVG_TILE_Q", &cluster_AVG_TILE_Q, &EM_b_cluster_AVG_TILE_Q);
   tree->SetBranchAddress("cluster_ENG_BAD_HV_CELLS", &cluster_ENG_BAD_HV_CELLS, &EM_b_cluster_ENG_BAD_HV_CELLS);
   tree->SetBranchAddress("cluster_N_BAD_HV_CELLS", &cluster_N_BAD_HV_CELLS, &EM_b_cluster_N_BAD_HV_CELLS);
   tree->SetBranchAddress("cluster_PTD", &cluster_PTD, &EM_b_cluster_PTD);
   tree->SetBranchAddress("cluster_MASS", &cluster_MASS, &EM_b_cluster_MASS);


   tree->SetBranchAddress("EM_Shower", &EM_Shower, &EM_b_EM_Shower);
   tree->SetBranchAddress("EM_Pro", &EM_Pro, &EM_b_EM_Pro);

   tree->SetBranchAddress("TrueE", &TrueE, &EM_b_TrueE);
   tree->SetBranchAddress("CalibratedE", &CalibratedE, &EM_b_CalibratedE);
   tree->SetBranchAddress("Corr", &Corr, &EM_b_Corr);
   Notify();
   // InitHist(tree);
}

void TopoClusters::Init_HAD(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fCurrent = -1;
   tree->SetMakeClass(1);

   tree->SetBranchAddress("runNumber", &runNumber, &Had_b_runNumber);
   tree->SetBranchAddress("eventNumber", &eventNumber, &Had_b_eventNumber);
   tree->SetBranchAddress("truthE", &truthE, &Had_b_truthE);
   tree->SetBranchAddress("truthPt", &truthPt, &Had_b_truthPt);
   tree->SetBranchAddress("truthEta", &truthEta, &Had_b_truthEta);
   tree->SetBranchAddress("truthPhi", &truthPhi, &Had_b_truthPhi);
   tree->SetBranchAddress("truthPDG", &truthPDG, &Had_b_truthPDG);
   tree->SetBranchAddress("nCluster", &nCluster, &Had_b_nCluster);
   tree->SetBranchAddress("clusterIndex", &clusterIndex, &Had_b_clusterIndex);
   tree->SetBranchAddress("cluster_nCells", &cluster_nCells, &Had_b_cluster_nCells);
   tree->SetBranchAddress("cluster_nCells_tot", &cluster_nCells_tot, &Had_b_cluster_nCells_tot);
//    tree->SetBranchAddress("clusterIndex", &clusterIndex, &Had_b_clusterIndex);
   tree->SetBranchAddress("clusterECalib", &clusterECalib, &Had_b_clusterECalib);
   tree->SetBranchAddress("clusterPtCalib", &clusterPtCalib, &Had_b_clusterPtCalib);
   tree->SetBranchAddress("clusterEtaCalib", &clusterEtaCalib, &Had_b_clusterEtaCalib);
   tree->SetBranchAddress("clusterPhiCalib", &clusterPhiCalib, &Had_b_clusterPhiCalib);
   tree->SetBranchAddress("cluster_sumCellECalib", &cluster_sumCellECalib, &Had_b_cluster_sumCellECAlib);
   tree->SetBranchAddress("clusterE", &clusterE, &Had_b_clusterE);
   tree->SetBranchAddress("clusterPt", &clusterPt, &Had_b_clusterPt);
   tree->SetBranchAddress("clusterEta", &clusterEta, &Had_b_clusterEta);
   tree->SetBranchAddress("clusterPhi", &clusterPhi, &Had_b_clusterPhi);
   tree->SetBranchAddress("cluster_sumCellE", &cluster_sumCellE, &Had_b_cluster_sumCellE);
   tree->SetBranchAddress("cluster_EM_PROBABILITY", &cluster_EM_PROBABILITY, &Had_b_cluster_EM_PROBABILITY);
   tree->SetBranchAddress("cluster_HAD_WEIGHT", &cluster_HAD_WEIGHT, &Had_b_cluster_HAD_WEIGHT);
   tree->SetBranchAddress("cluster_OOC_WEIGHT", &cluster_OOC_WEIGHT, &Had_b_cluster_OOC_WEIGHT);
   tree->SetBranchAddress("cluster_DM_WEIGHT", &cluster_DM_WEIGHT, &Had_b_cluster_DM_WEIGHT);
   tree->SetBranchAddress("cluster_ENG_CALIB_TOT", &cluster_ENG_CALIB_TOT, &Had_b_cluster_ENG_CALIB_TOT);
   tree->SetBranchAddress("cluster_ENG_CALIB_OUT_T", &cluster_ENG_CALIB_OUT_T, &Had_b_cluster_ENG_CALIB_OUT_T);
   tree->SetBranchAddress("cluster_ENG_CALIB_DEAD_TOT", &cluster_ENG_CALIB_DEAD_TOT, &Had_b_cluster_ENG_CALIB_DEAD_TOT);
   tree->SetBranchAddress("cluster_CENTER_MAG", &cluster_CENTER_MAG, &Had_b_cluster_CENTER_MAG);
   tree->SetBranchAddress("cluster_FIRST_ENG_DENS", &cluster_FIRST_ENG_DENS, &Had_b_cluster_FIRST_ENG_DENS);
   tree->SetBranchAddress("cluster_FIRST_PHI", &cluster_FIRST_PHI, &Had_b_cluster_FIRST_PHI);
   tree->SetBranchAddress("cluster_FIRST_ETA", &cluster_FIRST_ETA, &Had_b_cluster_FIRST_ETA);
   tree->SetBranchAddress("cluster_SECOND_R", &cluster_SECOND_R, &Had_b_cluster_SECOND_R);
   tree->SetBranchAddress("cluster_SECOND_LAMBDA", &cluster_SECOND_LAMBDA, &Had_b_cluster_SECOND_LAMBDA);
   tree->SetBranchAddress("cluster_DELTA_PHI", &cluster_DELTA_PHI, &Had_b_cluster_DELTA_PHI);
   tree->SetBranchAddress("cluster_DELTA_THETA", &cluster_DELTA_THETA, &Had_b_cluster_DELTA_THETA);
   tree->SetBranchAddress("cluster_DELTA_ALPHA", &cluster_DELTA_ALPHA, &Had_b_cluster_DELTA_ALPHA);
   tree->SetBranchAddress("cluster_CENTER_X", &cluster_CENTER_X, &Had_b_cluster_CENTER_X);
   tree->SetBranchAddress("cluster_CENTER_Y", &cluster_CENTER_Y, &Had_b_cluster_CENTER_Y);
   tree->SetBranchAddress("cluster_CENTER_Z", &cluster_CENTER_Z, &Had_b_cluster_CENTER_Z);
   tree->SetBranchAddress("cluster_CENTER_LAMBDA", &cluster_CENTER_LAMBDA, &Had_b_cluster_CENTER_LAMBDA);
   tree->SetBranchAddress("cluster_LATERAL", &cluster_LATERAL, &Had_b_cluster_LATERAL);
   tree->SetBranchAddress("cluster_LONGITUDINAL", &cluster_LONGITUDINAL, &Had_b_cluster_LONGITUDINAL);
   tree->SetBranchAddress("cluster_ENG_FRAC_EM", &cluster_ENG_FRAC_EM, &Had_b_cluster_ENG_FRAC_EM);
   tree->SetBranchAddress("cluster_ENG_FRAC_MAX", &cluster_ENG_FRAC_MAX, &Had_b_cluster_ENG_FRAC_MAX);
   tree->SetBranchAddress("cluster_ENG_FRAC_CORE", &cluster_ENG_FRAC_CORE, &Had_b_cluster_ENG_FRAC_CORE);
   tree->SetBranchAddress("cluster_SECOND_ENG_DENS", &cluster_SECOND_ENG_DENS, &Had_b_cluster_SECOND_ENG_DENS);
   tree->SetBranchAddress("cluster_ISOLATION", &cluster_ISOLATION, &Had_b_cluster_ISOLATION);
   tree->SetBranchAddress("cluster_ENG_BAD_CELLS", &cluster_ENG_BAD_CELLS, &Had_b_cluster_ENG_BAD_CELLS);
   tree->SetBranchAddress("cluster_N_BAD_CELLS", &cluster_N_BAD_CELLS, &Had_b_cluster_N_BAD_CELLS);
   tree->SetBranchAddress("cluster_N_BAD_CELLS_CORR", &cluster_N_BAD_CELLS_CORR, &Had_b_cluster_N_BAD_CELLS_CORR);
   tree->SetBranchAddress("cluster_BAD_CELLS_CORR_E", &cluster_BAD_CELLS_CORR_E, &Had_b_cluster_BAD_CELLS_CORR_E);
   tree->SetBranchAddress("cluster_BADLARQ_FRAC", &cluster_BADLARQ_FRAC, &Had_b_cluster_BADLARQ_FRAC);
   tree->SetBranchAddress("cluster_ENG_POS", &cluster_ENG_POS, &Had_b_cluster_ENG_POS);
   tree->SetBranchAddress("cluster_SIGNIFICANCE", &cluster_SIGNIFICANCE, &Had_b_cluster_SIGNIFICANCE);
   tree->SetBranchAddress("cluster_CELL_SIGNIFICANCE", &cluster_CELL_SIGNIFICANCE, &Had_b_cluster_CELL_SIGNIFICANCE);
   tree->SetBranchAddress("cluster_CELL_SIG_SAMPLING", &cluster_CELL_SIG_SAMPLING, &Had_b_cluster_CELL_SIG_SAMPLING);
   tree->SetBranchAddress("cluster_AVG_LAR_Q", &cluster_AVG_LAR_Q, &Had_b_cluster_AVG_LAR_Q);
   tree->SetBranchAddress("cluster_AVG_TILE_Q", &cluster_AVG_TILE_Q, &Had_b_cluster_AVG_TILE_Q);
   tree->SetBranchAddress("cluster_ENG_BAD_HV_CELLS", &cluster_ENG_BAD_HV_CELLS, &Had_b_cluster_ENG_BAD_HV_CELLS);
   tree->SetBranchAddress("cluster_N_BAD_HV_CELLS", &cluster_N_BAD_HV_CELLS, &Had_b_cluster_N_BAD_HV_CELLS);
   tree->SetBranchAddress("cluster_PTD", &cluster_PTD, &Had_b_cluster_PTD);
   tree->SetBranchAddress("cluster_MASS", &cluster_MASS, &Had_b_cluster_MASS);


   tree->SetBranchAddress("EM_Shower", &EM_Shower, &Had_b_EM_Shower);
   tree->SetBranchAddress("EM_Pro", &EM_Pro, &Had_b_EM_Pro);

   tree->SetBranchAddress("TrueE", &TrueE, &Had_b_TrueE);
   tree->SetBranchAddress("CalibratedE", &CalibratedE, &Had_b_CalibratedE);
   tree->SetBranchAddress("Corr", &Corr, &Had_b_Corr);
   Notify();
   // InitHist(tree);
}

void TopoClusters::InitHist(TTree *tree){

    // Make list of Histograms
    Histograms_List = new TList();
    if (!tree) return;
    fChain = tree;
    fChain->SetMakeClass(1);

    //Create Each Histogram and add it to the list.
    His = new TH1D("His", "X Distrubtion", 111, 0, 110); Histograms_List->Add(His);
    Corr_Histogram = new TH2D("Corr_Histogram", "Correlation fo the varibles", 5, -.5, 4.5, 5, -.5, 4.5); Histograms_List->Add(Corr_Histogram);
    // b_His = fChain->Branch("b_His", "TH1D", &His);
    // b_Histograms_List = fChain->Branch("Histograms_List", "TList", &Histograms_List);
    return;

}

void TopoClusters::Init_Matrix(){

  Product_Matrix.ResizeTo(n_branches, n_branches);
  Temp_Matrix.ResizeTo(n_branches, n_branches);
  Sum_Matrix.ResizeTo(n_branches,n_branches);
  // Corr->ResizeTo(5, 5);

}

void TopoClusters::Calculate_Correlation(){
  if (fChain == 0) return;
  // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results.root");
  // if (!f || !f->IsOpen()) {
  //    f = new TFile("results.root", "RECREATE");
  // }

        // printf("Interation: %lld\n", jentry);
     // if (Cut(ientry) < 0) continue;
  Long64_t nentries = fChain->GetEntriesFast();
  printf("Total Number of entries: %lli\n\n", nentries);
  Product_Matrix.Print();
  for (Int_t i = 0; i < 5; i++){
    for (Int_t j = 0; j < 5; j++){
        printf("Current Index: (%i, %i)\n", i, j);
        Temp_Matrix[i][j] = (nentries*(Product_Matrix[i][j])-Sum_Matrix[i][i]*Sum_Matrix[j][j])
        /sqrt((nentries*Product_Matrix[i][i]-Sum_Matrix[i][i]*Sum_Matrix[i][i])
        *(nentries*Product_Matrix[j][j]-Sum_Matrix[j][j]*Sum_Matrix[j][j]));
        // printf("\n");
    }
  }
  Corr = &Temp_Matrix;
  Corr->ResizeTo(Temp_Matrix);
  Corr->Print();
  return;
}

void TopoClusters::Fill_Correlation(){

  TTree *t = fChain;

  TLeaf *temp_leaf_j, *temp_leaf_i;
  Double_t i_val, j_val;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    for (Int_t i = 0; i < n_branches; i++){
      temp_leaf_i = t->FindLeaf(titles[i]);
      i_val = temp_leaf_i->GetValue(0);
      Sum_Matrix[i][i] += i_val;
      // printf("%f\n", i_val);
        for(Int_t j = i; j < n_branches; j++){
                temp_leaf_j = t->FindLeaf(titles[j]);
                j_val = temp_leaf_j->GetValue(0);
                Product_Matrix[i][j] += j_val * i_val;
          }
      }

      // printf("i value: %d\tj value: %d\n", i, j);
    }
    // for(Int_t i = 0; i<n_branches; i++){
    //   Sum_Matrix[i][i] = Sum_Matrix[i][i]/n_branches;
    // }
  return;
}


Bool_t TopoClusters::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TopoClusters::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TopoClusters::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TopoClusters_cxx
