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
   Int_t           id;
   Double_t        x;
   Double_t        y;
   Double_t        z;
   Double_t        Target;
   Double_t        Predicted;
   TH1D            *His;
   TH2D            *Corr_Histogram;
   TList           *Histograms_List;
   TMatrix         *Corr;
   TMatrix         Temp_Matrix;
   TMatrix         Sum_Matrix;
   TMatrix         Product_Matrix;
   // List of branches
   TBranch        *b_id;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_Target;   //!
   TBranch        *b_Predicted;   //!
   // TBranch        *b_His; //!
   // TBranch        *b_Histograms_List; //!
   TBranch        *b_Corr;

   TopoClusters(TTree *tree=0);
   virtual ~TopoClusters();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
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
      f->GetObject("TopoClusters",tree);

   }
   Init(tree);
   InitHist(tree);
   Init_Matrix();
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

   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("Target", &Target, &b_Target);
   fChain->SetBranchAddress("Predicted", &Predicted, &b_Predicted);
   fChain->SetBranchAddress("Corr", &Corr, &b_Corr);
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

  Product_Matrix.ResizeTo(5, 5);
  Temp_Matrix.ResizeTo(5, 5);
  Sum_Matrix.ResizeTo(5,5);
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
    Product_Matrix[0][0] += x*x;
    Product_Matrix[0][1] += x*y;
    Product_Matrix[0][2] += x*z;
    Product_Matrix[0][3] += x*Target;
    Product_Matrix[0][4] += x*Predicted;

    Product_Matrix[1][0] += y*x;
    Product_Matrix[1][1] += y*y;
    Product_Matrix[1][2] += y*z;
    Product_Matrix[1][3] += y*Target;
    Product_Matrix[1][4] += y*Predicted;

    Product_Matrix[2][0] += z*x;
    Product_Matrix[2][1] += z*y;
    Product_Matrix[2][2] += z*z;
    Product_Matrix[2][3] += z*Target;
    Product_Matrix[2][4] += z*Predicted;

    Product_Matrix[3][0] += Target*x;
    Product_Matrix[3][1] += Target*y;
    Product_Matrix[3][2] += Target*z;
    Product_Matrix[3][3] += Target*Target;
    Product_Matrix[3][4] += Target*Predicted;

    Product_Matrix[4][0] += Predicted*x;
    Product_Matrix[4][1] += Predicted*y;
    Product_Matrix[4][2] += Predicted*z;
    Product_Matrix[4][3] += Predicted*Target;
    Product_Matrix[4][4] += Predicted*Predicted;

    Sum_Matrix[0][0] += x;
    Sum_Matrix[1][1] += y;
    Sum_Matrix[2][2] += z;
    Sum_Matrix[3][3] += Target;
    Sum_Matrix[4][4] += Predicted;
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
