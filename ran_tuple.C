


void ran_tuple(){

  TFile* f = new TFile("Ran_Tuple.root","RECREATE");

  TNtuple *tuple = new TNtuple("tuple", "Randomly created tuples", "x:y:z:w");
  Double_t x, y, z, w;

  int number_tuples = 10000;
  TRandom3 * gen_random = new TRandom3();

  for (int i = 0; i < number_tuples; i++){
    x = 100 * gen_random->Rndm();
    y = 100 * gen_random->Rndm();
    z = 100 * gen_random->Rndm();
    w = x * gen_random->Gaus(1, .1) + 3 * y * gen_random->Gaus(1, .1) + 10 * z * gen_random->Gaus(1, .1);
    tuple->Fill(x, y, z, w);
  }
  tuple->Write();
  f->Close();
}
