#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TStopwatch.h>

#include <cmath>
#include <iostream>
#include <vector>

void MakeRNTuple(size_t nEntries = 100000,
                 const char *filename = "RNTupleFile.root")
{
   TStopwatch timer;
   timer.Start();

   auto model = ROOT::RNTupleModel::Create();

   auto eventID = model->MakeField<int>("eventID");

   // Beam and target
   auto beamP4 =
      model->MakeField<ROOT::Math::PxPyPzEVector>("beamP4");

   auto targetP4 =
      model->MakeField<ROOT::Math::PxPyPzEVector>("targetP4");

   // Event information
   auto nParticles =
      model->MakeField<int>("nParticles");

   auto particlePDG =
      model->MakeField<std::vector<int>>("particlePDG");

   auto particleQ =
      model->MakeField<std::vector<int>>("particleCharge");

   auto particleP4 =
      model->MakeField<
         std::vector<ROOT::Math::PxPyPzEVector>
      >("particleP4");

   auto vertex =
      model->MakeField<ROOT::Math::XYZVector>("vertex");

   auto writer =
      ROOT::RNTupleWriter::Recreate(std::move(model),
                              "events",
                              filename);

   auto db = TDatabasePDG::Instance();

   auto eMinus = db->GetParticle(11);
   auto ePlus  = db->GetParticle(-11);

   auto piPlus  = db->GetParticle(211);
   auto piMinus = db->GetParticle(-211);

   auto proton  = db->GetParticle(2212);

   const double mElectron = eMinus->Mass();
   const double mPion     = piPlus->Mass();
   const double mProton   = proton->Mass();

   TRandom3 rng(42);

   for (size_t i = 0; i < nEntries; ++i)
   {
      *eventID = i;

      //---------------------------------------------------------
      // Beam and target
      //---------------------------------------------------------

      const double Ebeam = 10.0;

      double pzBeam =
         std::sqrt(Ebeam * Ebeam -
                   mElectron * mElectron);

      *beamP4 =
         ROOT::Math::PxPyPzEVector(
            0.0, 0.0, pzBeam, Ebeam);

      *targetP4 =
         ROOT::Math::PxPyPzEVector(
            0.0, 0.0, 0.0, mProton);

      //---------------------------------------------------------
      // Vertex
      //---------------------------------------------------------

      double vx = rng.Gaus(0.0, 0.01);
      double vy = rng.Gaus(0.0, 0.01);
      double vz = rng.Gaus(0.0, 0.5);

      *vertex =
         ROOT::Math::XYZVector(vx, vy, vz);

      //---------------------------------------------------------
      // Multiplicity depends on vertex-z
      //---------------------------------------------------------

      double lambda =
         4.5 * std::exp(-std::abs(vz)/1.5);

      if (lambda < 0.5)
         lambda = 0.5;

      int n = rng.Poisson(lambda);

      if (n < 2)
         n = 2;

      *nParticles = n;

      particlePDG->resize(n);
      particleQ->resize(n);
      particleP4->resize(n);

      //---------------------------------------------------------
      // Multiplicity-energy correlation
      //---------------------------------------------------------

      double momentumScale =
         1.0 + 0.10 * n;

      //---------------------------------------------------------
      // Store momentum sums for balancing
      //---------------------------------------------------------

      double sumPx = 0.0;
      double sumPy = 0.0;

      //---------------------------------------------------------
      // Generate particles
      //---------------------------------------------------------

      for (int j = 0; j < n; ++j)
      {
         int pdg;
         int charge;
         double mass;
         double meanPz;

         //------------------------------------------------------
         // Force charge-balanced pairs frequently
         //------------------------------------------------------

         if (j > 0 && (j % 2 == 1) && rng.Rndm() < 0.70)
         {
            pdg = -(*particlePDG)[j-1];

            if (std::abs(pdg) == 11)
            {
               mass = mElectron;
               charge = pdg > 0 ? -1 : +1;
               meanPz = 1.0;
            }
            else if (std::abs(pdg) == 211)
            {
               mass = mPion;
               charge = pdg > 0 ? +1 : -1;
               meanPz = 2.0;
            }
            else
            {
               mass = mProton;
               charge = +1;
               meanPz = 3.0;
            }
         }
         else
         {
            double r = rng.Rndm();

            if (r < 0.40)
            {
               pdg = rng.Rndm() < 0.5 ? 11 : -11;
               mass = mElectron;
               charge = pdg > 0 ? -1 : +1;
               meanPz = 1.0;
            }
            else if (r < 0.85)
            {
               pdg = rng.Rndm() < 0.5 ? 211 : -211;
               mass = mPion;
               charge = pdg > 0 ? +1 : -1;
               meanPz = 2.0;
            }
            else
            {
               pdg = 2212;
               mass = mProton;
               charge = +1;
               meanPz = 3.0;
            }
         }

         //------------------------------------------------------
         // Approximate momentum conservation
         //------------------------------------------------------

         double px;
         double py;

         if (j == n - 1)
         {
            px = -sumPx + rng.Gaus(0.0, 0.05);
            py = -sumPy + rng.Gaus(0.0, 0.05);
         }
         else
         {
            px = rng.Gaus(0.0, 0.30*momentumScale);
            py = rng.Gaus(0.0, 0.30*momentumScale);

            sumPx += px;
            sumPy += py;
         }

         //------------------------------------------------------
         // Species-dependent pz spectra
         //------------------------------------------------------

         double pz =
            rng.Gaus(meanPz * momentumScale,
                     0.6 * momentumScale);

         if (pz < 0.05)
            pz = 0.05;

         //------------------------------------------------------
         // Detector-like smearing
         //------------------------------------------------------

         px *= rng.Gaus(1.0, 0.02);
         py *= rng.Gaus(1.0, 0.02);
         pz *= rng.Gaus(1.0, 0.02);

         double E =
            std::sqrt(px*px +
                      py*py +
                      pz*pz +
                      mass*mass);

         (*particlePDG)[j] = pdg;
         (*particleQ)[j]   = charge;

         (*particleP4)[j] =
            ROOT::Math::PxPyPzEVector(
               px, py, pz, E);
      }

      writer->Fill();
   }

   timer.Stop();

   std::cout
      << "Time to fill RNTuple: "
      << timer.RealTime()
      << " s\n";

   std::cout
      << "Created "
      << filename
      << std::endl;
}