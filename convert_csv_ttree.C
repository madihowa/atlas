// using namespace ROOT;
#include "TTree.h"
#include "TFile.h"
// #include "RDF.h"
#include <string>
#include <fstream>

#include <cstdio>

void convert(const std::string& fileName="NONE",const std::string& outFile="NONE") {
    if ( fileName == "NONE" || outFile == "NONE" ) { return; }
    //
    // std::ifstream ifs;
    // ifs.open(fileName.c_str());
    // if ( !ifs.is_open() ) {
    //   printf("[corrvar] cannot open file name \042%s\042\n",fileName.c_str());
    //   return;
    // }


    ROOT::RDataFrame rdf = ROOT::RDF::MakeCsvDataFrame(fileName);

    rdf.Snapshot("ClusterTree", outFile);


    //
    // //
    // char buffer[1024];
    // ifs.getline(buffer,1024);
    // printf("[corrvar] read header \042%s\042\n",buffer);
    // //
    // TFile* f = new TFile(outFile.c_str(),"RECREATE");
    // //
    // TTree* t = new TTree("TopoClusters","TopoClusters");
    // // int iflg(0); double pt(0.); double y(0.); double m(0.); double w(0.); double ptd(0.); int nc(0); double ptf(0.); double ncs(0.);
    // // t->Branch("iflg", &iflg, "iflg/I");
    // // t->Branch("pt",   &pt,   "pt/D"  );
    // // t->Branch("y",    &y,    "y/D"   );
    // // t->Branch("m",    &m,    "m/D"   );
    // // t->Branch("w",    &w,    "w/D"   );
    // // t->Branch("ptd",  &ptd,  "ptd/D" );
    // // t->Branch("nc",   &nc,   "nc/I"  );
    // // t->Branch("ptf",  &ptf,  "ptf/D" );
    // // t->Branch("ncs",  &ncs,  "ncs/D" );
    // //
    // printf("\n[corrvar] converting data from file \042%s\042 to tree <%s> in file \042%s\042\n",fileName.c_str(),t->GetName(),outFile.c_str());
    // // printf("[corrvar] converted %i entries\n\n",(int)t->ReadStream(ifs,"BCha/D:QCha/D:Flag/I:PT/D:y/D:m/D:w/D:ptD/D:Cons/D",','));
    // printf("[corrvar] converted %i entries\n\n",(int)t->ReadStream(ifs,"runNumber/I:eventNumber/I:truthPt/D:truthEta/D:truthPhi/D:truthPDG/D:nCluster/I:clusterIndex/I:cluster_nCells/I:cluster_nCells_tot/I:clusterIndex.1/D:clusterECalib/D:clusterPtCalib/D:clusterEtaCalib/D:clusterPhiCalib/D:cluster_sumCellECalib/D:clusterE/D:clusterPt/D:clusterEta/D:clusterPhi/D:cluster_sumCellE/D:cluster_EM_PROBABILITY/D:cluster_HAD_WEIGHT/D:cluster_OOC_WEIGHT/D:cluster_DM_WEIGHT/D:cluster_ENG_CALIB_TOT/D:cluster_ENG_CALIB_OUT_T/D:cluster_ENG_CALIB_DEAD_TOT/D:cluster_CENTER_MAG/D:cluster_FIRST_ENG_DENS/D:cluster_FIRST_PHI/D:cluster_FIRST_ETA/D:cluster_SECOND_R/D:cluster_SECOND_LAMBDA/D:cluster_DELTA_PHI/D:cluster_DELTA_THETA/D:cluster_DELTA_ALPHA/D:cluster_CENTER_X/D:cluster_CENTER_Y/D:cluster_CENTER_Z/D:cluster_CENTER_LAMBDA/D:cluster_LATERAL/D:cluster_LONGITUDINAL/D:cluster_ENG_FRAC_EM/D:cluster_ENG_FRAC_MAX/D:cluster_ENG_FRAC_CORE/D:cluster_SECOND_ENG_DENS/D:cluster_ISOLATION/D:cluster_ENG_BAD_CELLS/D:cluster_N_BAD_CELLS/D:cluster_N_BAD_CELLS_CORR/D:cluster_BAD_CELLS_CORR_E/D:cluster_BADLARQ_FRAC/D:cluster_ENG_POS/D:cluster_SIGNIFICANCE/D:cluster_CELL_SIGNIFICANCE/D:cluster_CELL_SIG_SAMPLING/D:cluster_AVG_LAR_Q/D:cluster_AVG_TILE_Q/D:cluster_ENG_BAD_HV_CELLS/D:cluster_N_BAD_HV_CELLS/D:cluster_PTD/D:cluster_MASS,truthE,Predicted Value",','));
    // //
    // t->Write();
    // f->Close();
}
