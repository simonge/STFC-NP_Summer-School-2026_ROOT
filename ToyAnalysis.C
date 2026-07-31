#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>

void ToyAnalysis(const char* filename = "AsymEvents.root",
                 const char* outfilename = "AsymAnalysis_RDF_bins_MM.root")
{
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("events", filename);

    // 1. Timing-based event weights (prompt–random subtraction)
    double sigW   = 1.0;
    double scale  = 10.0 / 160.0; // signal window / background window

    auto dfw = df.Define(
        "weight",
        [sigW, scale](double t) {
            if (t > -5.0 && t < 5.0)   return sigW;      // prompt
            if (t < -20.0 || t > 20.0) return -scale;    // random
            return 0.0;                                  // unused
        },
        {"time"}
    ).Filter("weight != 0.0");

    // 2. Define binning for energy and theta
    std::vector<double> E_bins = {400,600,1000,1500};
    std::vector<double> T_bins = {0,1.0,2.0,3.2};

    auto dfb = dfw
        .Define("Ebin", [E_bins](double E){
            for (int i=0;i<(int)E_bins.size()-1;++i)
                if (E>=E_bins[i] && E<E_bins[i+1]) return i;
            return -1;
        }, {"photonEnergy"})
        .Define("Tbin", [T_bins](double th){
            for (int i=0;i<(int)T_bins.size()-1;++i)
                if (th>=T_bins[i] && th<T_bins[i+1]) return i;
            return -1;
        }, {"theta"})
        .Filter("Ebin>=0 && Tbin>=0");

    int nE = (int)E_bins.size()-1;
    int nT = (int)T_bins.size()-1;

    // 3. Global missing-mass (prompt–random subtracted)
    auto hMM_global =
        dfw.Histo1D({"hMM_global",
                     "Missing mass (prompt-random subtracted);M_{miss} (MeV);Weighted counts",
                     200, 900, 980},
                    "missingMass", "weight");

    // 4. Per-bin missing-mass histograms
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> hMM(
        nE, std::vector<ROOT::RDF::RResultPtr<TH1D>>(nT));

    for (int iE=0;iE<nE;++iE) {
        for (int iT=0;iT<nT;++iT) {
            hMM[iE][iT] =
                dfb.Filter(Form("Ebin==%d && Tbin==%d", iE, iT))
                   .Histo1D(
                       {Form("hMM_E%d_T%d",iE,iT),
                        Form("Missing mass Ebin=%d Tbin=%d;M_{miss} (MeV);Weighted counts",iE,iT),
                        200, 900, 980},
                       "missingMass", "weight"
                   );
        }
    }

    // 5. Per-bin phi histograms (weighted)
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH1D>>> hphi(
        nE, std::vector<ROOT::RDF::RResultPtr<TH1D>>(nT));

    for (int iE=0;iE<nE;++iE) {
        for (int iT=0;iT<nT;++iT) {
            hphi[iE][iT] =
                dfb.Filter(Form("Ebin==%d && Tbin==%d", iE, iT))
                   .Histo1D(
                       {Form("hphi_E%d_T%d",iE,iT),
                        Form("phi Ebin=%d Tbin=%d;#phi;Weighted counts",iE,iT),
                        60, 0, 2*TMath::Pi()},
                       "phi", "weight"
                   );
        }
    }

    // 6. Fit N(1 + A cos(2phi)) in each bin
    TF1 fPhi("fPhi","[0]*(1 + [1]*cos(2*x))",0,2*TMath::Pi());

    std::vector<double> Avals;
    std::vector<double> Aerrs;
    std::vector<double> Ecenters;
    std::vector<double> Tcenters;

    for (int iE=0;iE<nE;++iE) {
        for (int iT=0;iT<nT;++iT) {

            TH1D* h = hphi[iE][iT].GetPtr();
            if (!h || h->GetEntries() < 40) continue;

            fPhi.SetParameters(h->GetMaximum(), 0.1);
            h->Fit(&fPhi, "Q");

            Avals.push_back(fPhi.GetParameter(1));
            Aerrs.push_back(fPhi.GetParError(1));
            Ecenters.push_back(0.5*(E_bins[iE] + E_bins[iE+1]));
            Tcenters.push_back(0.5*(T_bins[iT] + T_bins[iT+1]));
        }
    }

    // 7. Build graphs with errors
    int N = (int)Avals.size();

    TGraphErrors gA_E(N);
    gA_E.SetName("gA_E");
    gA_E.SetTitle("Asymmetry vs Energy;E_{#gamma} (MeV);A");

    TGraphErrors gA_T(N);
    gA_T.SetName("gA_T");
    gA_T.SetTitle("Asymmetry vs Theta;#theta (rad);A");

    for (int i=0;i<N;++i) {
        gA_E.SetPoint(i, Ecenters[i], Avals[i]);
        gA_E.SetPointError(i, 0.0, Aerrs[i]);

        gA_T.SetPoint(i, Tcenters[i], Avals[i]);
        gA_T.SetPointError(i, 0.0, Aerrs[i]);
    }

    // 8. Plot key results
    TCanvas cMM("cMM","Missing mass (prompt-random subtracted)",800,600);
    hMM_global->Draw();

    TCanvas cAE("cAE","Asymmetry vs Energy",800,600);
    gA_E.Draw("AP");

    TCanvas cAT("cAT","Asymmetry vs Theta",800,600);
    gA_T.Draw("AP");

    // 9. Save everything
    TFile fout(outfilename,"RECREATE");

    hMM_global->Write();

    for (int iE=0;iE<nE;++iE)
        for (int iT=0;iT<nT;++iT)
            hMM[iE][iT]->Write();

    for (int iE=0;iE<nE;++iE)
        for (int iT=0;iT<nT;++iT)
            hphi[iE][iT]->Write();

    gA_E.Write();
    gA_T.Write();
    cMM.Write();
    cAE.Write();
    cAT.Write();

    fout.Close();
}
