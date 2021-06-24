#define ClusterTree_cxx
#include "ClusterTree.h"
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
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    TLeaf* Temp;
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        Temp = fChain->FindLeaf(fChain->GetListOfBranches()->operator[](1)->GetName());
std:cout <<fChain->GetListOfBranches()->operator[](1)->GetName() << ": " << Temp->GetValue(0) << " " << eventNumber <<  endl;
    }



    TTree *t = fChain;
    ofstream csv_file;
    csv_file.open("ML_input.csv");
    char buffer[1024];

    Int_t n_branches = t->GetListOfBranches()->GetEntriesFast();

    const char * titles[n_branches];

    for(Int_t i = 0; i < n_branches; i++){
        titles[i] = t->GetListOfLeaves()->operator[](i)->GetName();
        csv_file << t->GetListOfBranches()->operator[](i)->GetName();
        if (i < n_branches -1){
            csv_file << ", ";
        }
    }
    csv_file << endl;

    Int_t n_entries = fChain->GetEntriesFast();

    Double_t error_amount = 0.0;
    TLeaf *temp_leaf;

    for(Int_t i = 0; i < n_entries; i++){
        Long64_t ientry = LoadTree(i);
        nb = fChain->GetEntry(i);
        for(Int_t j = 0; j < n_branches; j++){
            temp_leaf = t->FindLeaf(titles[j]);
            csv_file << temp_leaf->GetValue(0);
            // printf("Current Title:\n%s\n", titles[j]);
            cout << temp_leaf << endl;
            if(j < n_branches - 2){
                csv_file << ", ";
            }
            else{
                error_amount += cluster_MASS - temp_leaf->GetValue(0);
            }
        }
        csv_file << endl;
    }
    printf("%f\n", error_amount);
}

void ClusterTree::CombineFiles(){





    TFile * plus = new TFile("topo-cluster.pi+.root");
    TFile * minus = new TFile("topo-cluster.pi-.root");
    TFile * naught = new TFile("topo-cluster.pi0.root");

    TTree *Pi_plus = (TTree*)plus->Get("ClusterTree");
    TTree *Pi_minus = (TTree*)minus->Get("ClusterTree");
    TTree *Pi_naught = (TTree*)naught->Get("ClusterTree");


    TTree *test_tree;
    TTree *train_tree;


    TFile *f_train = (TFile*)gROOT->GetListOfFiles()->FindObject("topo-train.root");
    if (!f_train || !f_train->IsOpen()) {
        f_train = new TFile("topo-train.root");
    }

    TFile *f_test = (TFile*)gROOT->GetListOfFiles()->FindObject("topo-test.root");
    if (!f_test || !f_test->IsOpen()) {
        f_test = new TFile("topo-test.root");
    }


    f_test->GetObject("ClusterTree",test_tree);
    f_train->GetObject("ClusterTree",train_tree);

    Int_t n_plus_entries = (Int_t)Pi_plus->GetEntriesFast();
    Int_t n_minus_entries = (Int_t)Pi_minus->GetEntriesFast();
    Int_t n_naught_entries = (Int_t)Pi_naught->GetEntriesFast();

    Int_t n_max = n_plus;
    if (n_minus_entries > n_max) n_max = n_minus_entries;
    if (n_naught_entries > n_max) n_max = n_naught_entries;



    for (Int_t i = 0; i < (int)(n_max / 2) ; i++){
        Pi_plus->GetEntry(i);
        if (i < n_entries * .9) {
            train_tree->Fill();
        }
        else{
            test_tree->Fill();
        }
    }
    for (Int_t i = 0; i < (int)(n_max / 2); i++){
        Pi_minus->GetEntry(i);
        if (i < n_entries * .9) {
            train_tree->Fill();
        }
        else{
            test_tree->Fill();
        }
    }
    for (Int_t i = 0; i < n_max; i++){
        Pi_naught->GetEntry(i);
        if (i < n_entries * .9) {
            train_tree->Fill();
        }
        else{
            test_tree->Fill();
        }
    }
    train_tree->Print();

}
