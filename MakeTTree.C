#include <TFile.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TStopwatch.h>

#include <Math/Vector4D.h>
#include <Math/Vector3D.h>

#include <cmath>
#include <iostream>
#include <vector>

void MakeTTree(size_t nEntries = 100000,
               const char* filename = "TTreeFile.root")
{
   TStopwatch timer;
   timer.Start();

   //-------------------------------------------------------------
   // Output file and tree
   //-------------------------------------------------------------

   auto file = TFile::Open(filename, "RECREATE");

   TTree tree("events", "Simulated ep events");

   //-------------------------------------------------------------
   // Branch variables
   //-------------------------------------------------------------

   int eventID;

   ROOT::Math::PxPyPzEVector beamP4;
   ROOT::Math::PxPyPzEVector targetP4;

   int nParticles;

   std::vector<int> particlePDG;
   std::vector<int> particleCharge;

   std::vector<ROOT::Math::PxPyPzEVector> particleP4;

   ROOT::Math::XYZVector vertex;

   //-------------------------------------------------------------
   // Branches
   //-------------------------------------------------------------

   tree.Branch("eventID", &eventID);

   tree.Branch("beamP4", &beamP4);
   tree.Branch("targetP4", &targetP4);

   tree.Branch("nParticles", &nParticles);

   tree.Branch("particlePDG", &particlePDG);
   tree.Branch("particleCharge", &particleCharge);
   tree.Branch("particleP4", &particleP4);

   tree.Branch("vertex", &vertex);

   //-------------------------------------------------------------
   // PDG information
   //-------------------------------------------------------------

   auto db = TDatabasePDG::Instance();

   auto eMinus = db->GetParticle(11);
   auto piPlus = db->GetParticle(211);
   auto proton = db->GetParticle(2212);

   const double mElectron = eMinus->Mass();
   const double mPion     = piPlus->Mass();
   const double mProton   = proton->Mass();

   //-------------------------------------------------------------
   // Random generator
   //-------------------------------------------------------------

   TRandom3 rng(42);

   //-------------------------------------------------------------
   // Event loop
   //-------------------------------------------------------------

   for (size_t i = 0; i < nEntries; ++i)
   {
      eventID = i;

      //----------------------------------------------------------
      // Beam electron
      //----------------------------------------------------------

      constexpr double Ebeam = 10.0;

      double pzBeam =
         std::sqrt(Ebeam * Ebeam -
                   mElectron * mElectron);

      beamP4 =
         ROOT::Math::PxPyPzEVector(
            0.0,
            0.0,
            pzBeam,
            Ebeam);

      //----------------------------------------------------------
      // Target proton at rest
      //----------------------------------------------------------

      targetP4 =
         ROOT::Math::PxPyPzEVector(
            0.0,
            0.0,
            0.0,
            mProton);

      //----------------------------------------------------------
      // Vertex
      //----------------------------------------------------------

      double vx = rng.Gaus(0.0, 0.01);
      double vy = rng.Gaus(0.0, 0.01);
      double vz = rng.Gaus(0.0, 0.5);

      vertex =
         ROOT::Math::XYZVector(vx, vy, vz);

      //----------------------------------------------------------
      // Multiplicity correlated with vertex-z
      //----------------------------------------------------------

      double lambda =
         4.5 * std::exp(-std::abs(vz) / 1.5);

      if (lambda < 0.5)
         lambda = 0.5;

      nParticles = rng.Poisson(lambda);

      if (nParticles < 2)
         nParticles = 2;

      particlePDG.clear();
      particleCharge.clear();
      particleP4.clear();

      particlePDG.resize(nParticles);
      particleCharge.resize(nParticles);
      particleP4.resize(nParticles);

      //----------------------------------------------------------
      // Multiplicity-energy correlation
      //----------------------------------------------------------

      double momentumScale =
         1.0 + 0.10 * nParticles;

      double sumPx = 0.0;
      double sumPy = 0.0;

      //----------------------------------------------------------
      // Particle loop
      //----------------------------------------------------------

      for (int j = 0; j < nParticles; ++j)
      {
         int pdg;
         int charge;
         double mass;
         double meanPz;

         //-------------------------------------------------------
         // Frequently produce charge-balanced pairs
         //-------------------------------------------------------

         if (j > 0 && (j % 2 == 1) && rng.Rndm() < 0.70)
         {
            pdg = -particlePDG[j - 1];

            if (std::abs(pdg) == 11)
            {
               mass   = mElectron;
               charge = (pdg > 0) ? -1 : 1;
               meanPz = 1.0;
            }
            else if (std::abs(pdg) == 211)
            {
               mass   = mPion;
               charge = (pdg > 0) ? 1 : -1;
               meanPz = 2.0;
            }
            else
            {
               mass   = mProton;
               charge = 1;
               meanPz = 3.0;
            }
         }
         else
         {
            double r = rng.Rndm();

            if (r < 0.40)
            {
               pdg = (rng.Rndm() < 0.5) ? 11 : -11;

               mass   = mElectron;
               charge = (pdg > 0) ? -1 : 1;
               meanPz = 1.0;
            }
            else if (r < 0.85)
            {
               pdg = (rng.Rndm() < 0.5) ? 211 : -211;

               mass   = mPion;
               charge = (pdg > 0) ? 1 : -1;
               meanPz = 2.0;
            }
            else
            {
               pdg = 2212;

               mass   = mProton;
               charge = 1;
               meanPz = 3.0;
            }
         }

         //-------------------------------------------------------
         // Approximate transverse momentum conservation
         //-------------------------------------------------------

         double px;
         double py;

         if (j == nParticles - 1)
         {
            px = -sumPx + rng.Gaus(0.0, 0.05);
            py = -sumPy + rng.Gaus(0.0, 0.05);
         }
         else
         {
            px = rng.Gaus(
               0.0,
               0.30 * momentumScale);

            py = rng.Gaus(
               0.0,
               0.30 * momentumScale);

            sumPx += px;
            sumPy += py;
         }

         //-------------------------------------------------------
         // Particle-type dependent momentum spectra
         //-------------------------------------------------------

         double pz =
            rng.Gaus(meanPz * momentumScale,
                     0.6 * momentumScale);

         if (pz < 0.05)
            pz = 0.05;

         //-------------------------------------------------------
         // Detector smearing
         //-------------------------------------------------------

         px *= rng.Gaus(1.0, 0.02);
         py *= rng.Gaus(1.0, 0.02);
         pz *= rng.Gaus(1.0, 0.02);

         double E =
            std::sqrt(px * px +
                      py * py +
                      pz * pz +
                      mass * mass);

         particlePDG[j]    = pdg;
         particleCharge[j] = charge;

         particleP4[j] =
            ROOT::Math::PxPyPzEVector(
               px,
               py,
               pz,
               E);
      }

      tree.Fill();
   }

   //-------------------------------------------------------------
   // Write file
   //-------------------------------------------------------------

   file->cd();
   tree.Write();
   file->Close();

   timer.Stop();

   std::cout
      << "Time to fill TTree: "
      << timer.RealTime()
      << " seconds\n";

   std::cout
      << "Created "
      << filename
      << std::endl;
}