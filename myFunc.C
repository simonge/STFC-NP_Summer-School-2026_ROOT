void myFunc()
{

    //Create output file
    auto file = TFile::Open("myfile.root", "RECREATE");

    //Define a TF1 function with two parameters
    auto f1 = new TF1("fa","gaus",-3,3);
    f1->SetParameters(0.2,1.3,1);
    f1->Draw();

    auto h1 = new TH1F("h1","test",100,-3,3);
    h1->FillRandom("fa",2000);

    f1->SetParameters(0.1,1.0,11);
    h1->Fit("fa");

    h1->Write();
    file->Close();

}