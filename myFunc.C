void myFunc()
{

    TRandom3 rng(10);

    //Create output file
    auto file = TFile::Open("myfile.root", "RECREATE");

    //Define a TF1 function with two parameters
    auto f1 = new TF1("fa","[0]*TMath::BreitWigner(x, [1], [2])",-5,5);
    f1->SetParameters(100,0,1);
    f1->Draw();

    auto h1 = new TH1F("h1","test",100,-5,5);
    h1->FillRandom("fa",2000, &rng);

    h1->Fit("fa","R");

    h1->Write();
    file->Close();

}