

#include "TTree.h"
#include "TFile.h"


#include <string>
#include <fstream>
#include <iostream>
#include <cstdio>


void convert_ntuple_csv(const std::string& fileName="NONE",const std::string& outFile="NONE", const std::string& outDir=""){

    if ( fileName == "NONE" || outFile == "NONE" ) { return; }
    //
    TFile *f = new TFile("Ran_Tuple.root","READ");
    if ( !f ) {
        printf("[corrvar] cannot open file name \042%s\042\n",fileName.c_str());
        return;
    }
    ofstream csv_file;
    csv_file.open(outDir + "/ML_input.csv");



    char buffer[1024];

    TTree *t1 = (TTree*)f->Get("ClusterTree");

    Float_t x, y, z, w;
    t1->SetBranchAddress("x", &x);
    t1->SetBranchAddress("y", &y);
    t1->SetBranchAddress("z", &z);
    t1->SetBranchAddress("w", &w);

    csv_file << "x, y, z, w" << endl;

    Int_t nentries = (Int_t)t1->GetEntries();
    for (int i = 0; i<nentries; i++){
        t1->GetEntry(i);
        // printf("%f, %f, %f, %f\n", x, y, z, w);
        csv_file << x << ", " <<  y << ", " << z << ", " << w << endl;
    }
}
