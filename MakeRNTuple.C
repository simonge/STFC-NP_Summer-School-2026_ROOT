#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RField.hxx>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <TRandom3.h>
#include <TStopwatch.h>

void MakeRNTuple(size_t nEntries = 100000, const char* filename = "RNTupleFile.root")
{
   TStopwatch timer;
   timer.Start();

   auto model = ROOT::RNTupleModel::Create();

   auto eventID     = model->MakeField<int>("eventID");

   // Beam and target
   auto beamP4   = model->MakeField<ROOT::Math::PxPyPzEVector>("beamP4");
   auto targetP4 = model->MakeField<ROOT::Math::PxPyPzEVector>("targetP4");

   auto nParticles  = model->MakeField<int>("nParticles");
   auto particlePDG = model->MakeField<std::vector<int>>("particlePDG");
   auto particleQ   = model->MakeField<std::vector<int>>("particleCharge");
   auto particleP4  = model->MakeField<std::vector<ROOT::Math::PxPyPzEVector>>("particleP4");
   auto vertex      = model->MakeField<ROOT::Math::XYZVector>("vertex");

   auto writer = ROOT::RNTupleWriter::Recreate(std::move(model),
					       "events",
					       filename);
                      
   // PDG database
    auto db = TDatabasePDG::Instance();

    auto pElectron = db->GetParticle(11);
    auto pPion     = db->GetParticle(211);
    auto pProton   = db->GetParticle(2212);

    float mElectron = pElectron->Mass();
    float mPion     = pPion->Mass();
    float mProton   = pProton->Mass();
                      
   
   TRandom3 rng(42);

   for (size_t i = 0; i < nEntries; i++) {

      *eventID = i;

      // Beam electron (10 GeV)
      float Ebeam = 10.0;
      float pzbeam = std::sqrt(Ebeam*Ebeam - mElectron*mElectron);
      *beamP4 = ROOT::Math::PxPyPzEVector(0, 0, pzbeam, Ebeam);

      // Target proton at rest
      *targetP4 = ROOT::Math::PxPyPzEVector(0, 0, 0, mProton);


      // Vertex (simple Gaussian)
      *vertex = ROOT::Math::XYZVector(rng.Gaus(0, 0.01),
                                       rng.Gaus(0, 0.01),
                     rng.Gaus(0, 0.5));
      // Reconstructed particles
      int n = rng.Poisson(3);
      *nParticles = n;

      particlePDG->resize(n);
      particleQ->resize(n);
      particleP4->resize(n);

      for (int j = 0; j < n; j++) {

         // Choose particle type
         int pdg;
         float mass;
         int charge;

         float r = rng.Rndm();
         if (r < 0.5) {
               pdg = 11;
               mass = mElectron;
               charge = pElectron->Charge() / 3.0;
         } else if (r < 0.8) {
               pdg = 211;
               mass = mPion;
               charge = pPion->Charge() / 3.0;
         } else {
               pdg = 2212;
               mass = mProton;
               charge = pProton->Charge() / 3.0;
         }

         (*particlePDG)[j] = pdg;
         (*particleQ)[j]   = charge;

         // Simple momentum model
         float px = rng.Gaus(0, 0.3)/n;
         float py = rng.Gaus(0, 0.3)/n;
         float pz = rng.Gaus(2.0, 1.0)/n;
         float E  = std::sqrt(px*px + py*py + pz*pz + mass*mass);

         (*particleP4)[j] = ROOT::Math::PxPyPzEVector(px, py, pz, E);
      }

      writer->Fill();
   }

   timer.Stop();
   std::cout << "Time to fill RNTuple: " << timer.RealTime() << " seconds" << std::endl;


   std::cout << "Created " << filename << std::endl;
}
