#include "TRandom.h"
#include "TStopwatch.h"

void ToyGenerator( size_t nEntries = 100000,
               const char* filename = "AsymEvents.root" )
{

    TStopwatch timer;
    timer.Start();

    TRandom rng(0);

    // Create output file and tree
    auto file = TFile::Open(filename, "RECREATE");
    TTree tree("events", "Simulated ep events");

    //Create branches for theta, phi, missing mass, photon energy and time
    double theta, phi, missingMass, photonEnergy, time;
    tree.Branch("theta", &theta);
    tree.Branch("phi", &phi);
    tree.Branch("missingMass", &missingMass);
    tree.Branch("photonEnergy", &photonEnergy);
    tree.Branch("time", &time);

    // Max photon energy from beam
    double Emax = 1500;

    // Pi0 photoproduction threshold off a proton
    double Ethr = 144.7;

    //Function for the threshold Breit-Wigner distribution  
    TF1 *fXS = new TF1("fXS",
        "[0] * sqrt(x - [1]) * TMath::BreitWigner(x, [2], [3])",
        Ethr, 1500);

    fXS->SetParameters(1.0, Ethr, 1232.0, 120.0); // A, threshold, mean, width

    //Phi function
    TF1 *fPhi = new TF1("fPhi", "[0] * (1 + [1] * cos(2 * x))", 0, 2 * 3.14159);
    fPhi->SetParameters(1.0, 0.5); // A, asymmetry parameter

    // Loop over events
    for (size_t i = 0; i < nEntries; ++i) {
        // Generate signal event
        missingMass = rng.Gaus(938, 20); // missing mass around proton mass
        double truePhotonEnergy = fXS->GetRandom()-938;
        photonEnergy = truePhotonEnergy; // photon energy with exponential distribution

        double costheta = rng.Uniform(-1, 1); // costheta in [-1, 1]
        theta = TMath::ACos(costheta); // theta in [0, pi]
        // Phi asummetry based on the generated theta value and energy of the photon
        fPhi->SetParameter(1, (1 - costheta*costheta)*photonEnergy/Emax); // asymmetry parameter depends on theta and photon energy
        phi = fPhi->GetRandom();
        time = rng.Gaus(0, 2); // time in [0, 100] ns

        // Fill the tree with the generated values
        tree.Fill();

        // Generate out of time background events with a poisson distribution
        int nBackground = rng.Poisson(100); // average 20 background events per signal event
        for (int j = 0; j < nBackground; ++j) {
            photonEnergy = rng.Exp(200); // different exponential distribution
            // if (photonEnergy > Emax) j--; // reject events with photon energy above Emax
            // continue; // skip to next iteration if photon energy is above Emax
            missingMass = missingMass+photonEnergy-truePhotonEnergy; // missing mass is the sum of the signal and background photon energies
            time = rng.Uniform(-100, 100);

            tree.Fill();
        }

    }

    // Write the tree to the file
    file->cd();
    tree.Write();
    file->Close();

    std::cout << "Generated " << nEntries << " events in " << timer.RealTime() << " seconds." << std::endl;
    std::cout << "Output file: " << filename << std::endl;

}



