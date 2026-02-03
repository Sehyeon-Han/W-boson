void RecoilHist_MC(const char* DYMC = "/data2/kplee/Lecture/CMSOpenData/MC2016/DY_M50_aMCNLO/*.root",
                   const char* outFile = "Recoil_ZMC.root")
{
    TChain chain("Events");
    chain.Add(DYMC);

    TTreeReader reader(&chain);

    TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
    TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
    TTreeReaderArray<int> Muon_charge(reader, "Muon_charge");
    TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
    TTreeReaderValue<float> MET_pt(reader, "MET_pt");
    TTreeReaderValue<float> MET_phi(reader, "MET_phi");
    TTreeReaderValue<float> genWeight(reader, "genWeight");
    TTreeReaderArray<float> Muon_mass(reader, "Muon_mass");

    TH2D *h1 = new TH2D("h1", "U_{1} vs p_{T}^{Z};p_{T}^{Z}[GeV];U_{1}[GeV]", 20, 0, 100, 20, -80, 0);
    TH2D *h2 = new TH2D("h2", "U_{2} vs p_{T}^{Z};p_{T}^{Z}[GeV];U_{2}[GeV}", 20, 0, 100, 20, -80, 80);
    
    auto Muon = [&](int i)
    {
        TLorentzVector mu;
        mu.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
        return mu;
    };

    while(reader.Next())
    {   
        if(!(*HLT_IsoMu24)) continue;

        vector<int> MuIdx;
        
        const int n = Muon_pt.GetSize();
        for(int i = 0; i < n; ++i)
        {
            if(Muon_pt[i] > 25.0 && fabs(Muon_eta[i]) < 2.4 && Muon_iso[i] < 0.12 && Muon_tightId[i])
            {
                MuIdx.push_back(i);
            }
        }
        if(MuIdx.size() < 2) continue;

        int iMu1 = MuIdx[0];
        int iMu2 = MuIdx[1];

        if(Muon_charge[iMu1] * Muon_charge[iMu2] >= 0) continue;
        
        TLorentzVector mu1 = Muon(iMu1);
        TLorentzVector mu2 = Muon(iMu2);
        TLorentzVector Z = mu1 + mu2;
        
        double mass = Z.M();
        if (mass < 60.0 || mass > 120.0) continue;
        
        double ptZ = Z.Pt();
        double phiZ = Z.Phi();

        double metPx = (*MET_pt) * cos(*MET_phi);
        double metPy = (*MET_pt) * sin(*MET_phi);

        double Ux = -(metPx + mu1.Px() + mu2.Px());
        double Uy = -(metPy + mu1.Py() + mu2.Py());

        double U1 = Ux * cos(phiZ) + Uy * sin(phiZ);
        double U2 = -Ux * sin(phiZ) + Uy * cos(phiZ);

        double weight = *genWeight;

        h1->Fill(ptZ, U1, weight);
        h2->Fill(ptZ, U2, weight);
    }

    TFile *fout = TFile::Open(outFile, "RECREATE");
    h1->Write();
    h2->Write();
    fout->Close();

    printf(">>> RecoilHist done. Output = %s\n", outFile);
}
