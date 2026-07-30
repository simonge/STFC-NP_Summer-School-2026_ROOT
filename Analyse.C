#include <ROOT/RDataFrame.hxx>

#include <TCanvas.h>

void Analyse( const char* filename = "RNTupleFile.root", const char* outFilename = "AnalysisOutput.root", 
              const char* plotFilename = "plots.pdf")
{
   ROOT::EnableImplicitMT();

   ROOT::RDataFrame df(
      "events",
      filename);

   
   auto selected =
     df.Filter("nParticles > 3")
       .Define("vertexR", "sqrt(vertex.X()*vertex.X() + vertex.Y()*vertex.Y()+ vertex.Z()*vertex.Z())")
         .Filter("vertexR < 1.0")
         .Define("chargeSum", "std::accumulate(particleCharge.begin(), particleCharge.end(), 0)");
   
   std::cout
     << "Events passing cuts = "
     << *selected.Count()
     << std::endl;

   auto hCharge =
     selected.Histo1D(
         {"hCharge",
          "Charge Sum;#Sigma Q;Events",
          22,
          -11.0,
          11.0},
         "chargeSum");

   auto hVertexR =
     selected.Histo1D(
         {"hVertexR",
          "Vertex Radius;#sqrt{x^{2}+y^{2}+z^{2}} [cm];Events",
          100,
          0.0,
          1.0},
         "vertexR");

   TCanvas c("c","c",1200,500);

   c.Divide(2,1);

   c.cd(1);
   hCharge->Draw();

   c.cd(2);
   hVertexR->Draw();

   c.SaveAs(plotFilename);

   //Create snapshot of the filtered data frame
   ROOT::RDF::RSnapshotOptions opts;
   opts.fOutputFormat = ROOT::RDF::ESnapshotOutputFormat::kRNTuple;
   selected.Snapshot("events", outFilename,
                     {"vertexR", "chargeSum", "nParticles"},opts);

   //Save histograms to same file
   TFile outFile(outFilename, "UPDATE");
   hCharge->Write();
   hVertexR->Write();
   outFile.Close();

   // Save the computational graph to the dot file
   ROOT::RDF::SaveGraph(hCharge, "computational_graph.dot");

}
