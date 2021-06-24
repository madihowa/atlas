//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 26 15:59:48 2020 by ROOT version 6.19/01
// from TTree ClusterTree/ClusterTree
// found on file: results.root
//////////////////////////////////////////////////////////

#ifndef ClusterTree_h
#define ClusterTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ClusterTree {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent; //!current Tree number in a TChain

        // Fixed size dimensions of array or collections stored in the TTree if any.

        // Declaration of leaf types and branches
        Long64_t        runNumber;
        Double_t        eventNumber;
        Double_t        truthE;
        Double_t        truthPt;
        Double_t        truthEta;
        Double_t        truthPhi;
        Long64_t        truthPDG;
        Long64_t        nCluster;
        Long64_t        clusterIndex;
        Long64_t        cluster_nCells;
        Long64_t        cluster_nCells_tot;
        Long64_t        clusterIndex_1;
        Double_t        clusterECalib;
        Double_t        clusterPtCalib;
        Double_t        clusterEtaCalib;
        Double_t        clusterPhiCalib;
        Double_t        cluster_sumCellECalib;
        Double_t        clusterE;
        Double_t        clusterPt;
        Double_t        clusterEta;
        Double_t        clusterPhi;
        Double_t        cluster_sumCellE;
        Double_t        cluster_EM_PROBABILITY;
        Double_t        cluster_HAD_WEIGHT;
        Double_t        cluster_OOC_WEIGHT;
        Double_t        cluster_DM_WEIGHT;
        Double_t        cluster_ENG_CALIB_TOT;
        Double_t        cluster_ENG_CALIB_OUT_T;
        Double_t        cluster_ENG_CALIB_DEAD_TOT;
        Double_t        cluster_CENTER_MAG;
        Double_t        cluster_FIRST_ENG_DENS;
        Double_t        cluster_FIRST_PHI;
        Double_t        cluster_FIRST_ETA;
        Double_t        cluster_SECOND_R;
        Double_t        cluster_SECOND_LAMBDA;
        Double_t        cluster_DELTA_PHI;
        Double_t        cluster_DELTA_THETA;
        Double_t        cluster_DELTA_ALPHA;
        Double_t        cluster_CENTER_X;
        Double_t        cluster_CENTER_Y;
        Double_t        cluster_CENTER_Z;
        Double_t        cluster_CENTER_LAMBDA;
        Double_t        cluster_LATERAL;
        Double_t        cluster_LONGITUDINAL;
        Double_t        cluster_ENG_FRAC_EM;
        Double_t        cluster_ENG_FRAC_MAX;
        Double_t        cluster_ENG_FRAC_CORE;
        Double_t        cluster_SECOND_ENG_DENS;
        Double_t        cluster_ISOLATION;
        Double_t        cluster_ENG_BAD_CELLS;
        Long64_t        cluster_N_BAD_CELLS;
        Long64_t        cluster_N_BAD_CELLS_CORR;
        Double_t        cluster_BAD_CELLS_CORR_E;
        Long64_t        cluster_BADLARQ_FRAC;
        Double_t        cluster_ENG_POS;
        Double_t        cluster_SIGNIFICANCE;
        Double_t        cluster_CELL_SIGNIFICANCE;
        Long64_t        cluster_CELL_SIG_SAMPLING;
        Double_t        cluster_AVG_LAR_Q;
        Double_t        cluster_AVG_TILE_Q;
        Double_t        cluster_ENG_BAD_HV_CELLS;
        Long64_t        cluster_N_BAD_HV_CELLS;
        Double_t        cluster_PTD;
        Double_t        cluster_MASS;
        Long64_t        EM_Shower;
        Double_t        EM_Pro;
        Double_t        CalibratedE;
        Double_t        Delta_E;
        Double_t        Delta_Calib_E;


        // List of branches
        TBranch        *b_runNumber;   //!
        TBranch        *b_eventNumber;   //!
        TBranch        *b_truthE;   //!
        TBranch        *b_truthPt;   //!
        TBranch        *b_truthEta;   //!
        TBranch        *b_truthPhi;   //!
        TBranch        *b_truthPDG;   //!
        TBranch        *b_nCluster;   //!
        TBranch        *b_clusterIndex;   //!
        TBranch        *b_cluster_nCells;   //!
        TBranch        *b_cluster_nCells_tot;   //!
        TBranch        *b_clusterIndex_1;   //!
        TBranch        *b_clusterECalib;   //!
        TBranch        *b_clusterPtCalib;   //!
        TBranch        *b_clusterEtaCalib;   //!
        TBranch        *b_clusterPhiCalib;   //!
        TBranch        *b_cluster_sumCellECalib;   //!
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
        TBranch        *b_EM_Shower;   //!
        TBranch        *b_EM_Pro;   //!
        TBranch        *b_CalibratedE;   //!
        TBranch        *b_Delta_E;   //!
        TBranch        *b_Delta_Calib_E;   //!


        // List of branches
        TBranch        *EM_b_runNumber;   //!
        TBranch        *EM_b_eventNumber;   //!
        TBranch        *EM_b_truthE;   //!
        TBranch        *EM_b_truthPt;   //!
        TBranch        *EM_b_truthEta;   //!
        TBranch        *EM_b_truthPhi;   //!
        TBranch        *EM_b_truthPDG;   //!
        TBranch        *EM_b_nCluster;   //!
        TBranch        *EM_b_clusterIndex;   //!
        TBranch        *EM_b_cluster_nCells;   //!
        TBranch        *EM_b_cluster_nCells_tot;   //!
        TBranch        *EM_b_clusterIndex_1;   //!
        TBranch        *EM_b_clusterECalib;   //!
        TBranch        *EM_b_clusterPtCalib;   //!
        TBranch        *EM_b_clusterEtaCalib;   //!
        TBranch        *EM_b_clusterPhiCalib;   //!
        TBranch        *EM_b_cluster_sumCellECalib;   //!
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
        TBranch        *EM_b_EM_Shower;   //!
        TBranch        *EM_b_EM_Pro;   //!
        TBranch        *EM_b_CalibratedE;   //!
        TBranch        *EM_b_Delta_E;   //!
        TBranch        *EM_b_Delta_Calib_E;   //!


        // List of branches
        TBranch        *Had_b_runNumber;   //!
        TBranch        *Had_b_eventNumber;   //!
        TBranch        *Had_b_truthE;   //!
        TBranch        *Had_b_truthPt;   //!
        TBranch        *Had_b_truthEta;   //!
        TBranch        *Had_b_truthPhi;   //!
        TBranch        *Had_b_truthPDG;   //!
        TBranch        *Had_b_nCluster;   //!
        TBranch        *Had_b_clusterIndex;   //!
        TBranch        *Had_b_cluster_nCells;   //!
        TBranch        *Had_b_cluster_nCells_tot;   //!
        TBranch        *Had_b_clusterIndex_1;   //!
        TBranch        *Had_b_clusterECalib;   //!
        TBranch        *Had_b_clusterPtCalib;   //!
        TBranch        *Had_b_clusterEtaCalib;   //!
        TBranch        *Had_b_clusterPhiCalib;   //!
        TBranch        *Had_b_cluster_sumCellECalib;   //!
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
        TBranch        *Had_b_EM_Shower;   //!
        TBranch        *Had_b_EM_Pro;   //!
        TBranch        *Had_b_CalibratedE;   //!
        TBranch        *Had_b_Delta_E;   //!
        TBranch        *Had_b_Delta_Calib_E;   //!

        ClusterTree(TTree *tree=0);
        virtual ~ClusterTree();
        virtual Int_t    Cut(Long64_t entry);
        virtual Int_t    GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void     Init(TTree *tree);
        virtual void     Init_EM(TTree *tree);
        virtual void     Init_Had(TTree *tree);
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ClusterTree_cxx
ClusterTree::ClusterTree(TTree *tree) : fChain(0)
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

    TTree *Had_ptr = new TTree("Had_tree", "Hadron shower Tree");
    TTree *EM_ptr = new TTree("EM_tree", "EM shower Tree");
    // EM_tree = EM_ptr;
    // Had_tree = Had_ptr;
    Init(tree);
    // Init_Had(Had_tree);
    // Init_EM(EM_tree);
    fChain=tree;
    std::cout << "1 Current number of Entries: " << fChain->GetEntriesFast() << endl;

}

ClusterTree::~ClusterTree()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t ClusterTree::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t ClusterTree::LoadTree(Long64_t entry)
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

void ClusterTree::Init(TTree *tree)
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

    // fChain->SetBranchAddress("clusterECalib", &clusterECalib, &b_clusterECalib);
    // fChain->SetBranchAddress("clusterPtCalib", &clusterPtCalib, &b_clusterPtCalib);
    // fChain->SetBranchAddress("clusterEtaCalib", &clusterEtaCalib, &b_clusterEtaCalib);
    // fChain->SetBranchAddress("clusterPhiCalib", &clusterPhiCalib, &b_clusterPhiCalib);

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
    fChain->SetBranchAddress("clusterIndex_1", &clusterIndex_1, &b_clusterIndex_1);
    fChain->SetBranchAddress("clusterECalib", &clusterECalib, &b_clusterECalib);
    fChain->SetBranchAddress("clusterPtCalib", &clusterPtCalib, &b_clusterPtCalib);
    fChain->SetBranchAddress("clusterEtaCalib", &clusterEtaCalib, &b_clusterEtaCalib);
    fChain->SetBranchAddress("clusterPhiCalib", &clusterPhiCalib, &b_clusterPhiCalib);
    fChain->SetBranchAddress("cluster_sumCellECalib", &cluster_sumCellECalib, &b_cluster_sumCellECalib);
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
    fChain->SetBranchAddress("CalibratedE", &CalibratedE, &b_CalibratedE);
    b_Delta_E = fChain->Branch("Delta_E", &Delta_E, "Delta_E/D");
    b_Delta_Calib_E = fChain->Branch("Delta_Calib_E", &Delta_Calib_E, "Delta_Calib_E/D");
    fChain->SetBranchAddress("Delta_E", &Delta_E, &b_Delta_E);
    fChain->SetBranchAddress("Delta_Calib_E", &Delta_Calib_E, &b_Delta_Calib_E);

    Notify();
}

void ClusterTree::Init_EM(TTree *tree)
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

    // EM_b_clusterECalib = tree->Branch("clusterECalib", &clusterECalib, "clusterECalib/D");  //!
    // EM_b_clusterPtCalib = tree->Branch("clusterPtCalib", &clusterPtCalib, "clusterPtCalib/D");   //!
    // EM_b_clusterEtaCalib = tree->Branch("clusterEtaCalib", &clusterEtaCalib, "clusterEtaCalib/D");   //!
    // EM_b_clusterPhiCalib = tree->Branch("clusterPhiCalib", &clusterPhiCalib, "clusterPhiCalib/D");  //!

    // EM_b_clusterE = tree->Branch("clusterE", &clusterE, "clusterE/D");  //!
    // EM_b_clusterPt = tree->Branch("clusterPt", &clusterPt, "clusterPt/D");   //!
    // EM_b_clusterEta = tree->Branch("clusterEta", &clusterEta, "clusterEta/D");   //!
    // EM_b_clusterPhi = tree->Branch("clusterPhi", &clusterPhi, "clusterPhi/D");  //!
    // EM_b_cluster_PTD = tree->Branch("cluster_PTD", &cluster_PTD, "cluster_PTD/D");    //!
    //
    //
    // EM_b_cluster_CENTER_LAMBDA = tree->Branch("cluster_CENTER_LAMBDA", &cluster_CENTER_LAMBDA, "cluster_CENTER_LAMBDA/D");   //!
    // EM_b_cluster_FIRST_ENG_DENS = tree->Branch("cluster_FIRST_ENG_DENS", &cluster_FIRST_ENG_DENS, "cluster_FIRST_ENG_DENS/D");  //!
    //
    // // EM_b_cluster_EM_PROBABILITY = tree->Branch("cluster_EM_PROBABILITY", &cluster_EM_PROBABILITY, "cluster_EM_PROBABILITY/D");  //!
    // EM_b_EM_Shower = tree->Branch("EM_Shower", &EM_Shower, "EM_Shower/I");   //!
    // EM_b_CalibratedE = tree->Branch("CalibratedE", &EM_Shower, "CalibratedE/I");   //!
    // EM_b_TruthE = tree->Branch("TruthE", &EM_Shower, "TruthE/I");   //!
    // EM_b_EM_Pro = tree->Branch("EM_Pro", &EM_Pro, "EM_Pro/D");
    // EM_b_Delta_E = tree->Branch("Delta_E", &Delta_E, "Delta_E/D");
    //

    EM_b_runNumber = tree->Branch("runNumber", &runNumber, "runNumber/I");
    EM_b_eventNumber = tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
    EM_b_truthE = tree->Branch("truthE", &truthE, "truthE/D");
    EM_b_truthPt = tree->Branch("truthPt", &truthPt, "truthPt/D");
    EM_b_truthEta = tree->Branch("truthEta", &truthEta, "truthEta/D");
    EM_b_truthPhi = tree->Branch("truthPhi", &truthPhi, "truthPhi/D");
    EM_b_truthPDG = tree->Branch("truthPDG", &truthPDG, "truthPDG/I");
    EM_b_nCluster = tree->Branch("nCluster", &nCluster, "nCluster/I");
    EM_b_clusterIndex = tree->Branch("clusterIndex", &clusterIndex, "clusterIndex/I");
    EM_b_cluster_nCells = tree->Branch("cluster_nCells", &cluster_nCells, "cluster_nCells/I");
    EM_b_cluster_nCells_tot = tree->Branch("cluster_nCells_tot", &cluster_nCells_tot, "cluster_nCells_tot/I");
    EM_b_clusterIndex_1 = tree->Branch("clusterIndex_1", &clusterIndex_1, "clusterIndex_1/I");
    EM_b_clusterECalib = tree->Branch("clusterECalib", &clusterECalib, "clusterECalib/D");
    EM_b_clusterPtCalib = tree->Branch("clusterPtCalib", &clusterPtCalib, "clusterPtCalib/D");
    EM_b_clusterEtaCalib = tree->Branch("clusterEtaCalib", &clusterEtaCalib, "clusterEtaCalib/D");
    EM_b_clusterPhiCalib = tree->Branch("clusterPhiCalib", &clusterPhiCalib, "clusterPhiCalib/D");
    EM_b_cluster_sumCellECalib = tree->Branch("cluster_sumCellECalib", &cluster_sumCellECalib, "cluster_sumCellECalib/D");
    EM_b_clusterE = tree->Branch("clusterE", &clusterE, "clusterE/D");
    EM_b_clusterPt = tree->Branch("clusterPt", &clusterPt, "clusterPt/D");
    EM_b_clusterEta = tree->Branch("clusterEta", &clusterEta, "clusterEta/D");
    EM_b_clusterPhi = tree->Branch("clusterPhi", &clusterPhi, "clusterPhi/D");
    EM_b_cluster_sumCellE = tree->Branch("cluster_sumCellE", &cluster_sumCellE, "cluster_sumCellE/D");
    EM_b_cluster_EM_PROBABILITY = tree->Branch("cluster_EM_PROBABILITY", &cluster_EM_PROBABILITY, "cluster_EM_PROBABILITY/D");
    EM_b_cluster_HAD_WEIGHT = tree->Branch("Blacluster_HAD_WEIGHTnk", &cluster_HAD_WEIGHT, "cluster_HAD_WEIGHT/D");
    EM_b_cluster_OOC_WEIGHT = tree->Branch("cluster_OOC_WEIGHT", &cluster_OOC_WEIGHT, "cluster_OOC_WEIGHT/D");
    EM_b_cluster_DM_WEIGHT = tree->Branch("cluster_DM_WEIGHT", &cluster_DM_WEIGHT, "cluster_DM_WEIGHT/D");
    EM_b_cluster_ENG_CALIB_TOT = tree->Branch("cluster_ENG_CALIB_TOT", &cluster_ENG_CALIB_TOT, "cluster_ENG_CALIB_TOT/D");
    EM_b_cluster_ENG_CALIB_OUT_T = tree->Branch("cluster_ENG_CALIB_OUT_T", &cluster_ENG_CALIB_OUT_T, "cluster_ENG_CALIB_OUT_T/D");
    EM_b_cluster_ENG_CALIB_DEAD_TOT = tree->Branch("cluster_ENG_CALIB_DEAD_TOT", &cluster_ENG_CALIB_DEAD_TOT, "cluster_ENG_CALIB_DEAD_TOT/D");
    EM_b_cluster_CENTER_MAG = tree->Branch("cluster_CENTER_MAG", &cluster_CENTER_MAG, "cluster_CENTER_MAG/D");
    EM_b_cluster_FIRST_ENG_DENS = tree->Branch("cluster_FIRST_ENG_DENS", &cluster_FIRST_ENG_DENS, "cluster_FIRST_ENG_DENS/D");
    EM_b_cluster_FIRST_PHI = tree->Branch("cluster_FIRST_PHI", &cluster_FIRST_PHI, "cluster_FIRST_PHI/D");
    EM_b_cluster_FIRST_ETA = tree->Branch("cluster_FIRST_ETA", &cluster_FIRST_ETA, "cluster_FIRST_ETA/D");
    EM_b_cluster_SECOND_R = tree->Branch("cluster_SECOND_R", &cluster_SECOND_R, "cluster_SECOND_R/D");
    EM_b_cluster_SECOND_LAMBDA = tree->Branch("cluster_SECOND_LAMBDA", &cluster_SECOND_LAMBDA, "cluster_SECOND_LAMBDA/D");
    EM_b_cluster_DELTA_PHI = tree->Branch("cluster_DELTA_PHI", &cluster_DELTA_PHI, "cluster_DELTA_PHI/D");
    EM_b_cluster_DELTA_THETA = tree->Branch("cluster_DELTA_THETA", &cluster_DELTA_THETA, "cluster_DELTA_THETA/D");
    EM_b_cluster_DELTA_ALPHA = tree->Branch("cluster_DELTA_ALPHA", &cluster_DELTA_ALPHA, "cluster_DELTA_ALPHA/D");
    EM_b_cluster_CENTER_X = tree->Branch("cluster_CENTER_X", &cluster_CENTER_X, "cluster_CENTER_X/D");
    EM_b_cluster_CENTER_Y = tree->Branch("cluster_CENTER_Y", &cluster_CENTER_Y, "cluster_CENTER_Y/D");
    EM_b_cluster_CENTER_Z = tree->Branch("cluster_CENTER_Z", &cluster_CENTER_Z, "cluster_CENTER_Z/D");
    EM_b_cluster_CENTER_LAMBDA = tree->Branch("cluster_CENTER_LAMBDA", &cluster_CENTER_LAMBDA, "cluster_CENTER_LAMBDA/D");
    EM_b_cluster_LATERAL = tree->Branch("cluster_LATERAL", &cluster_LATERAL, "cluster_LATERAL/D");
    EM_b_cluster_LONGITUDINAL = tree->Branch("cluster_LONGITUDINAL", &cluster_LONGITUDINAL, "cluster_LONGITUDINAL/D");
    EM_b_cluster_ENG_FRAC_EM = tree->Branch("cluster_ENG_FRAC_EM", &cluster_ENG_FRAC_EM, "cluster_ENG_FRAC_EM/D");
    EM_b_cluster_ENG_FRAC_MAX = tree->Branch("cluster_ENG_FRAC_MAX", &cluster_ENG_FRAC_MAX, "cluster_ENG_FRAC_MAX/D");
    EM_b_cluster_ENG_FRAC_CORE = tree->Branch("cluster_ENG_FRAC_CORE", &cluster_ENG_FRAC_CORE, "cluster_ENG_FRAC_CORE/D");
    EM_b_cluster_SECOND_ENG_DENS = tree->Branch("cluster_SECOND_ENG_DENS", &cluster_SECOND_ENG_DENS, "cluster_SECOND_ENG_DENS/D");
    EM_b_cluster_ISOLATION = tree->Branch("cluster_ISOLATION", &cluster_ISOLATION, "cluster_ISOLATION/D");
    EM_b_cluster_ENG_BAD_CELLS = tree->Branch("cluster_ENG_BAD_CELLS", &cluster_ENG_BAD_CELLS, "cluster_ENG_BAD_CELLS/D");
    EM_b_cluster_N_BAD_CELLS = tree->Branch("cluster_N_BAD_CELLS", &cluster_N_BAD_CELLS, "cluster_N_BAD_CELLS/D");
    EM_b_cluster_N_BAD_CELLS_CORR = tree->Branch("cluster_N_BAD_CELLS_CORR", &cluster_N_BAD_CELLS_CORR, "cluster_N_BAD_CELLS_CORR/D");
    EM_b_cluster_BAD_CELLS_CORR_E = tree->Branch("cluster_BAD_CELLS_CORR_E", &cluster_BAD_CELLS_CORR_E, "cluster_BAD_CELLS_CORR_E/D");
    EM_b_cluster_BADLARQ_FRAC = tree->Branch("cluster_BADLARQ_FRAC", &cluster_BADLARQ_FRAC, "cluster_BADLARQ_FRAC/D");
    EM_b_cluster_ENG_POS = tree->Branch("cluster_ENG_POS", &cluster_ENG_POS, "cluster_ENG_POS/D");
    EM_b_cluster_SIGNIFICANCE = tree->Branch("cluster_SIGNIFICANCE", &cluster_SIGNIFICANCE, "cluster_SIGNIFICANCE/D");
    EM_b_cluster_CELL_SIGNIFICANCE = tree->Branch("cluster_CELL_SIGNIFICANCE", &cluster_CELL_SIGNIFICANCE, "cluster_CELL_SIGNIFICANCE/D");
    EM_b_cluster_CELL_SIG_SAMPLING = tree->Branch("cluster_CELL_SIG_SAMPLING", &cluster_CELL_SIG_SAMPLING, "cluster_CELL_SIG_SAMPLING/D");
    EM_b_cluster_AVG_LAR_Q = tree->Branch("cluster_AVG_LAR_Q", &cluster_AVG_LAR_Q, "cluster_AVG_LAR_Q/D");
    EM_b_cluster_AVG_TILE_Q = tree->Branch("cluster_AVG_TILE_Q", &cluster_AVG_TILE_Q, "cluster_AVG_TILE_Q/D");
    EM_b_cluster_ENG_BAD_HV_CELLS = tree->Branch("cluster_ENG_BAD_HV_CELLS", &cluster_ENG_BAD_HV_CELLS, "cluster_ENG_BAD_HV_CELLS/D");
    EM_b_cluster_N_BAD_HV_CELLS = tree->Branch("cluster_N_BAD_HV_CELLS", &cluster_N_BAD_HV_CELLS, "cluster_N_BAD_HV_CELLS/D");
    EM_b_cluster_PTD = tree->Branch("cluster_PTD", &cluster_PTD, "cluster_PTD/D");
    EM_b_cluster_MASS = tree->Branch("cluster_MASS", &cluster_MASS, "cluster_MASS/D");
    EM_b_EM_Shower = tree->Branch("EM_Shower", &EM_Shower, "EM_Shower/D");
    EM_b_EM_Pro = tree->Branch("EM_Pro", &EM_Pro, "EM_Pro/D");
    EM_b_CalibratedE = tree->Branch("CalibratedE", &CalibratedE, "CalibratedE/D");
    EM_b_Delta_E = tree->Branch("Delta_E", &Delta_E, "Delta_E/D");
    EM_b_Delta_Calib_E = tree->Branch("Delta_Calib_E", &Delta_Calib_E, "Delta_Calib_E/D");







    Notify();
    return;
}

void ClusterTree::Init_Had(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // Init() will be called many times when running on PROOF
    // code, but the routine can be extended by the user if needed.
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;


    // Had_b_clusterECalib = tree->Branch("clusterECalib", &clusterECalib, "clusterECalib/D");  //!
    // Had_b_clusterPtCalib = tree->Branch("clusterPtCalib", &clusterPtCalib, "clusterPtCalib/D");   //!
    // Had_b_clusterEtaCalib = tree->Branch("clusterEtaCalib", &clusterEtaCalib, "clusterEtaCalib/D");   //!
    // Had_b_clusterPhiCalib = tree->Branch("clusterPhiCalib", &clusterPhiCalib, "clusterPhiCalib/D");  //!

    // Had_b_clusterE = tree->Branch("clusterE", &clusterE, "clusterE/D");  //!
    // Had_b_clusterPt = tree->Branch("clusterPt", &clusterPt, "clusterPt/D");   //!
    // Had_b_clusterEta = tree->Branch("clusterEta", &clusterEta, "clusterEta/D");   //!
    // Had_b_clusterPhi = tree->Branch("clusterPhi", &clusterPhi, "clusterPhi/D");  //!
    // Had_b_cluster_PTD = tree->Branch("cluster_PTD", &cluster_PTD, "cluster_PTD/D");    //!
    //
    //
    // Had_b_cluster_CENTER_LAMBDA = tree->Branch("cluster_CENTER_LAMBDA", &cluster_CENTER_LAMBDA, "cluster_CENTER_LAMBDA/D");   //!
    // Had_b_cluster_FIRST_ENG_DENS = tree->Branch("cluster_FIRST_ENG_DENS", &cluster_FIRST_ENG_DENS, "cluster_FIRST_ENG_DENS/D");  //!
    //
    // // Had_b_cluster_EM_PROBABILITY = tree->Branch("cluster_EM_PROBABILITY", &cluster_EM_PROBABILITY, "cluster_EM_PROBABILITY/D");  //!
    // Had_b_EM_Shower = tree->Branch("EM_Shower", &EM_Shower, "EM_Shower/I");   //!
    // Had_b_EM_Pro = tree->Branch("EM_Pro", &EM_Pro, "EM_Pro/D");
    // Had_b_Delta_E = tree->Branch("Delta_E", &Delta_E, "Delta_E/D");
    // Had_b_CalibratedE = tree->Branch("CalibratedE", &EM_Shower, "CalibratedE/I");   //!
    // Had_b_TruthE = tree->Branch("TruthE", &EM_Shower, "TruthE/I");   //!

    Had_b_runNumber = tree->Branch("runNumber", &runNumber, "runNumber/I");
    Had_b_eventNumber = tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
    Had_b_truthE = tree->Branch("truthE", &truthE, "truthE/D");
    Had_b_truthPt = tree->Branch("truthPt", &truthPt, "truthPt/D");
    Had_b_truthEta = tree->Branch("truthEta", &truthEta, "truthEta/D");
    Had_b_truthPhi = tree->Branch("truthPhi", &truthPhi, "truthPhi/D");
    Had_b_truthPDG = tree->Branch("truthPDG", &truthPDG, "truthPDG/I");
    Had_b_nCluster = tree->Branch("nCluster", &nCluster, "nCluster/I");
    Had_b_clusterIndex = tree->Branch("clusterIndex", &clusterIndex, "clusterIndex/I");
    Had_b_cluster_nCells = tree->Branch("cluster_nCells", &cluster_nCells, "cluster_nCells/I");
    Had_b_cluster_nCells_tot = tree->Branch("cluster_nCells_tot", &cluster_nCells_tot, "cluster_nCells_tot/I");
    Had_b_clusterIndex_1 = tree->Branch("clusterIndex_1", &clusterIndex_1, "clusterIndex_1/I");
    Had_b_clusterECalib = tree->Branch("clusterECalib", &clusterECalib, "clusterECalib/D");
    Had_b_clusterPtCalib = tree->Branch("clusterPtCalib", &clusterPtCalib, "clusterPtCalib/D");
    Had_b_clusterEtaCalib = tree->Branch("clusterEtaCalib", &clusterEtaCalib, "clusterEtaCalib/D");
    Had_b_clusterPhiCalib = tree->Branch("clusterPhiCalib", &clusterPhiCalib, "clusterPhiCalib/D");
    Had_b_cluster_sumCellECalib = tree->Branch("cluster_sumCellECalib", &cluster_sumCellECalib, "cluster_sumCellECalib/D");
    Had_b_clusterE = tree->Branch("clusterE", &clusterE, "clusterE/D");
    Had_b_clusterPt = tree->Branch("clusterPt", &clusterPt, "clusterPt/D");
    Had_b_clusterEta = tree->Branch("clusterEta", &clusterEta, "clusterEta/D");
    Had_b_clusterPhi = tree->Branch("clusterPhi", &clusterPhi, "clusterPhi/D");
    Had_b_cluster_sumCellE = tree->Branch("cluster_sumCellE", &cluster_sumCellE, "cluster_sumCellE/D");
    Had_b_cluster_EM_PROBABILITY = tree->Branch("cluster_EM_PROBABILITY", &cluster_EM_PROBABILITY, "cluster_EM_PROBABILITY/D");
    Had_b_cluster_HAD_WEIGHT = tree->Branch("Blacluster_HAD_WEIGHTnk", &cluster_HAD_WEIGHT, "cluster_HAD_WEIGHT/D");
    Had_b_cluster_OOC_WEIGHT = tree->Branch("cluster_OOC_WEIGHT", &cluster_OOC_WEIGHT, "cluster_OOC_WEIGHT/D");
    Had_b_cluster_DM_WEIGHT = tree->Branch("cluster_DM_WEIGHT", &cluster_DM_WEIGHT, "cluster_DM_WEIGHT/D");
    Had_b_cluster_ENG_CALIB_TOT = tree->Branch("cluster_ENG_CALIB_TOT", &cluster_ENG_CALIB_TOT, "cluster_ENG_CALIB_TOT/D");
    Had_b_cluster_ENG_CALIB_OUT_T = tree->Branch("cluster_ENG_CALIB_OUT_T", &cluster_ENG_CALIB_OUT_T, "cluster_ENG_CALIB_OUT_T/D");
    Had_b_cluster_ENG_CALIB_DEAD_TOT = tree->Branch("cluster_ENG_CALIB_DEAD_TOT", &cluster_ENG_CALIB_DEAD_TOT, "cluster_ENG_CALIB_DEAD_TOT/D");
    Had_b_cluster_CENTER_MAG = tree->Branch("cluster_CENTER_MAG", &cluster_CENTER_MAG, "cluster_CENTER_MAG/D");
    Had_b_cluster_FIRST_ENG_DENS = tree->Branch("cluster_FIRST_ENG_DENS", &cluster_FIRST_ENG_DENS, "cluster_FIRST_ENG_DENS/D");
    Had_b_cluster_FIRST_PHI = tree->Branch("cluster_FIRST_PHI", &cluster_FIRST_PHI, "cluster_FIRST_PHI/D");
    Had_b_cluster_FIRST_ETA = tree->Branch("cluster_FIRST_ETA", &cluster_FIRST_ETA, "cluster_FIRST_ETA/D");
    Had_b_cluster_SECOND_R = tree->Branch("cluster_SECOND_R", &cluster_SECOND_R, "cluster_SECOND_R/D");
    Had_b_cluster_SECOND_LAMBDA = tree->Branch("cluster_SECOND_LAMBDA", &cluster_SECOND_LAMBDA, "cluster_SECOND_LAMBDA/D");
    Had_b_cluster_DELTA_PHI = tree->Branch("cluster_DELTA_PHI", &cluster_DELTA_PHI, "cluster_DELTA_PHI/D");
    Had_b_cluster_DELTA_THETA = tree->Branch("cluster_DELTA_THETA", &cluster_DELTA_THETA, "cluster_DELTA_THETA/D");
    Had_b_cluster_DELTA_ALPHA = tree->Branch("cluster_DELTA_ALPHA", &cluster_DELTA_ALPHA, "cluster_DELTA_ALPHA/D");
    Had_b_cluster_CENTER_X = tree->Branch("cluster_CENTER_X", &cluster_CENTER_X, "cluster_CENTER_X/D");
    Had_b_cluster_CENTER_Y = tree->Branch("cluster_CENTER_Y", &cluster_CENTER_Y, "cluster_CENTER_Y/D");
    Had_b_cluster_CENTER_Z = tree->Branch("cluster_CENTER_Z", &cluster_CENTER_Z, "cluster_CENTER_Z/D");
    Had_b_cluster_CENTER_LAMBDA = tree->Branch("cluster_CENTER_LAMBDA", &cluster_CENTER_LAMBDA, "cluster_CENTER_LAMBDA/D");
    Had_b_cluster_LATERAL = tree->Branch("cluster_LATERAL", &cluster_LATERAL, "cluster_LATERAL/D");
    Had_b_cluster_LONGITUDINAL = tree->Branch("cluster_LONGITUDINAL", &cluster_LONGITUDINAL, "cluster_LONGITUDINAL/D");
    Had_b_cluster_ENG_FRAC_EM = tree->Branch("cluster_ENG_FRAC_EM", &cluster_ENG_FRAC_EM, "cluster_ENG_FRAC_EM/D");
    Had_b_cluster_ENG_FRAC_MAX = tree->Branch("cluster_ENG_FRAC_MAX", &cluster_ENG_FRAC_MAX, "cluster_ENG_FRAC_MAX/D");
    Had_b_cluster_ENG_FRAC_CORE = tree->Branch("cluster_ENG_FRAC_CORE", &cluster_ENG_FRAC_CORE, "cluster_ENG_FRAC_CORE/D");
    Had_b_cluster_SECOND_ENG_DENS = tree->Branch("cluster_SECOND_ENG_DENS", &cluster_SECOND_ENG_DENS, "cluster_SECOND_ENG_DENS/D");
    Had_b_cluster_ISOLATION = tree->Branch("cluster_ISOLATION", &cluster_ISOLATION, "cluster_ISOLATION/D");
    Had_b_cluster_ENG_BAD_CELLS = tree->Branch("cluster_ENG_BAD_CELLS", &cluster_ENG_BAD_CELLS, "cluster_ENG_BAD_CELLS/D");
    Had_b_cluster_N_BAD_CELLS = tree->Branch("cluster_N_BAD_CELLS", &cluster_N_BAD_CELLS, "cluster_N_BAD_CELLS/D");
    Had_b_cluster_N_BAD_CELLS_CORR = tree->Branch("cluster_N_BAD_CELLS_CORR", &cluster_N_BAD_CELLS_CORR, "cluster_N_BAD_CELLS_CORR/D");
    Had_b_cluster_BAD_CELLS_CORR_E = tree->Branch("cluster_BAD_CELLS_CORR_E", &cluster_BAD_CELLS_CORR_E, "cluster_BAD_CELLS_CORR_E/D");
    Had_b_cluster_BADLARQ_FRAC = tree->Branch("cluster_BADLARQ_FRAC", &cluster_BADLARQ_FRAC, "cluster_BADLARQ_FRAC/D");
    Had_b_cluster_ENG_POS = tree->Branch("cluster_ENG_POS", &cluster_ENG_POS, "cluster_ENG_POS/D");
    Had_b_cluster_SIGNIFICANCE = tree->Branch("cluster_SIGNIFICANCE", &cluster_SIGNIFICANCE, "cluster_SIGNIFICANCE/D");
    Had_b_cluster_CELL_SIGNIFICANCE = tree->Branch("cluster_CELL_SIGNIFICANCE", &cluster_CELL_SIGNIFICANCE, "cluster_CELL_SIGNIFICANCE/D");
    Had_b_cluster_CELL_SIG_SAMPLING = tree->Branch("cluster_CELL_SIG_SAMPLING", &cluster_CELL_SIG_SAMPLING, "cluster_CELL_SIG_SAMPLING/D");
    Had_b_cluster_AVG_LAR_Q = tree->Branch("cluster_AVG_LAR_Q", &cluster_AVG_LAR_Q, "cluster_AVG_LAR_Q/D");
    Had_b_cluster_AVG_TILE_Q = tree->Branch("cluster_AVG_TILE_Q", &cluster_AVG_TILE_Q, "cluster_AVG_TILE_Q/D");
    Had_b_cluster_ENG_BAD_HV_CELLS = tree->Branch("cluster_ENG_BAD_HV_CELLS", &cluster_ENG_BAD_HV_CELLS, "cluster_ENG_BAD_HV_CELLS/D");
    Had_b_cluster_N_BAD_HV_CELLS = tree->Branch("cluster_N_BAD_HV_CELLS", &cluster_N_BAD_HV_CELLS, "cluster_N_BAD_HV_CELLS/D");
    Had_b_cluster_PTD = tree->Branch("cluster_PTD", &cluster_PTD, "cluster_PTD/D");
    Had_b_cluster_MASS = tree->Branch("cluster_MASS", &cluster_MASS, "cluster_MASS/D");
    Had_b_EM_Shower = tree->Branch("EM_Shower", &EM_Shower, "EM_Shower/D");
    Had_b_EM_Pro = tree->Branch("EM_Pro", &EM_Pro, "EM_Pro/D");
    Had_b_CalibratedE = tree->Branch("CalibratedE", &CalibratedE, "CalibratedE/D");
    Had_b_Delta_E = tree->Branch("Delta_E", &Delta_E, "Delta_E/D");
    Had_b_Delta_Calib_E = tree->Branch("Delta_Calib_E", &Delta_Calib_E, "Delta_Calib_E/D");





    Notify();
    return;
}

Bool_t ClusterTree::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void ClusterTree::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t ClusterTree::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef ClusterTree_cxx
