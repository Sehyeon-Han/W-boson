void Recoil_Z_MC(const char* inFile = "/data2/kplee/Lecture/CMSOpenData/MC2016/DY_M50_aMCNLO/*.root",
                 const char* outFile = "recoil_Z_MC.root",
                 double ptMin = 0,
                 double ptMax = 100,
                 double dPt = 2)
{
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);
    
    TChain chain("Events");
    chain.Add(inFile);

    TTreeReader reader(&chain);

    TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<float> Muon_mass(reader, "Muon_mass");
    TTreeReaderArray<int> Muon_charge(reader, "Muon_charge");
    TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
    TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
    TTreeReaderValue<float> MET_pt(reader, "MET_pt");
    TTreeReaderValue<float> MET_phi(reader, "MET_phi");
    TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
    TTreeReaderValue<float> genWeight(reader, "genWeight");

    int nBins = (ptMax - ptMin)/dPt;

    vector<TH1D*> hU1(nBins), hU2(nBins);

    for(int i = 0; i < nBins; i++)
    {
        double low = ptMin + i * dPt;
        double high = low + dPt;
        hU1[i] = new TH1D(Form("hU1_%d",i), Form("U1 %.0f-%.0f;U1 [GeV]", low, high), 240, -120, 120);
        hU2[i] = new TH1D(Form("hU2_%d",i), Form("U2 %.0f-%.0f;U2 [GeV]", low, high), 240, -120, 120);
        hU1[i]->Sumw2();
        hU2[i]->Sumw2();
    }

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
        if(ptZ<ptMin || ptZ>=ptMax) continue;
        int ibin = (ptZ-ptMin)/dPt;
        double phiZ = Z.Phi();

        double metPx = (*MET_pt) * cos(*MET_phi);
        double metPy = (*MET_pt) * sin(*MET_phi);

        double Ux = -(metPx + mu1.Px() + mu2.Px());
        double Uy = -(metPy + mu1.Py() + mu2.Py());

        double U1 = Ux * cos(phiZ) + Uy * sin(phiZ);
        double U2 = -Ux * sin(phiZ) + Uy * cos(phiZ);
        
        double weight = *genWeight;

        hU1[ibin]->Fill(U1, weight);
        hU2[ibin]->Fill(U2, weight);
    }

    vector<double> x(nBins), ex(nBins, 0.0);
    vector<double> mean1(nBins, 0.0), emean1(nBins, 0.0);
    vector<double> mean2(nBins, 0.0), emean2(nBins, 0.0);
    vector<double> sig1(nBins, 0.0), esig1(nBins, 0.0);
    vector<double> sig2(nBins, 0.0), esig2(nBins, 0.0);
    vector<TF1*> fitU1(nBins, nullptr);
    vector<TF1*> fitU2(nBins, nullptr);

    for(int i = 0; i < nBins; i++)
    {
        x[i] = ptMin + (i + 0.5) * dPt;
        
        if(hU1[i]->GetEntries() > 50)
        {
            double m = hU1[i]->GetMean();
            double r = hU1[i]->GetRMS();
            double lo = hU1[i]->GetXaxis()->GetXmin();
            double hi = hU1[i]->GetXaxis()->GetXmax();

            fitU1[i] = new TF1(Form("fitU1_%d", i), "gaus", lo, hi);
            fitU1[i]->SetParameters(hU1[i]->GetMaximum(), m, r);
            hU1[i]->Fit(fitU1[i], "RQ0");

            mean1[i]  = fitU1[i]->GetParameter(1);
            emean1[i] = fitU1[i]->GetParError(1);

            sig1[i] = fitU1[i]->GetParameter(2);
            esig1[i] = fitU1[i]->GetParError(2);
        }
        
        if(hU2[i]->GetEntries() > 50)
        {
            double m = hU2[i]->GetMean();
            double r = hU2[i]->GetRMS();
            double lo = hU2[i]->GetXaxis()->GetXmin();
            double hi = hU2[i]->GetXaxis()->GetXmax();
            
            fitU2[i] = new TF1(Form("fitU2_%d", i), "gaus", lo, hi);
            fitU2[i]->SetParameters(hU2[i]->GetMaximum(), m, r);
            hU2[i]->Fit(fitU2[i], "RQ0");

            mean2[i] = fitU2[i]->GetParameter(1);
            emean2[i] = fitU2[i]->GetParError(1);

            sig2[i] = fitU2[i]->GetParameter(2);
            esig2[i] = fitU2[i]->GetParError(2);
        }
    }
    
    auto gMeanU1 = new TGraphErrors(nBins, &x[0], &mean1[0], &ex[0], &emean1[0]);
    gMeanU1->SetName("g_mean_U1");

    auto gMeanU2 = new TGraphErrors(nBins, &x[0], &mean2[0], &ex[0], &emean2[0]);
    gMeanU2->SetName("g_mean_U2");

    auto gSigU1 = new TGraphErrors(nBins, &x[0], &sig1[0], &ex[0], &esig1[0]);
    gSigU1->SetName("g_sigma_U1");

    auto gSigU2 = new TGraphErrors(nBins, &x[0], &sig2[0], &ex[0], &esig2[0]);
    gSigU2->SetName("g_sigma_U2");

    TFile *fout = new TFile(outFile, "RECREATE");

    for(int i = 0; i < nBins; i++)
    {
        hU1[i]->Write();
        hU2[i]->Write();
        if(fitU1[i]) fitU1[i]->Write();
        if(fitU2[i]) fitU2[i]->Write();
    }

    gMeanU1->Write();
    gMeanU2->Write();
    gSigU1->Write();
    gSigU2->Write();

    fout->Close();
}

void DrawBinFit(const char* rootFile = "recoil_Z_MC.root",
                int ibin = 0,
                const char* which = "U1")
{
    TFile *f = TFile::Open(rootFile);

    TH1D* h = (TH1D*)f->Get(Form("h%s_%d", which, ibin));
    TF1* g = (TF1*) f->Get(Form("fit%s_%d", which, ibin));

    TCanvas *c = new TCanvas(Form("c_%s_%d", which, ibin), "Fit Gaussian", 800, 600);
    h->Draw("E");
    if(g) g->Draw("same");
    c->Update();
}

void DrawRecoilGraphs(const char* rootFile = "recoil_Z_MC.root")
{
    TFile *f = TFile::Open(rootFile, "UPDATE");

    auto gMeanU1 = (TGraphErrors*)f->Get("g_mean_U1");
    auto gMeanU2 = (TGraphErrors*)f->Get("g_mean_U2");
    auto gSigU1 = (TGraphErrors*)f->Get("g_sigma_U1");
    auto gSigU2 = (TGraphErrors*)f->Get("g_sigma_U2");

    if(gMeanU1)
    {
        TCanvas *c1 = new TCanvas("c_meanU1","<U1> vs pT^Z", 800, 600);
        gMeanU1->SetTitle(";Z p_{T} [GeV];#LT U_{1} #GT [GeV]");
        gMeanU1->SetMarkerStyle(kOpenCircle);
        gMeanU1->SetMarkerSize(0.9);
        gMeanU1->SetLineWidth(1);
        gMeanU1->Fit("pol1");
        gMeanU1->Draw("AP");

        f->cd();
        TF1* f1 = gMeanU1->GetFunction("pol1");
        f1->SetName("f_meanU1");
        f1->Write();
    }

    if(gMeanU2)
    {
        TCanvas *c2 = new TCanvas("c_meanU2","<U2> vs pT^Z", 800, 600);
        gMeanU2->SetTitle(";Z p_{T} [GeV];#LT U_{2} #GT [GeV]");
        gMeanU2->SetMarkerStyle(kOpenCircle);
        gMeanU2->SetMarkerSize(0.9);
        gMeanU2->SetLineWidth(1);
        gMeanU2->Fit("pol1");
        gMeanU2->Draw("AP");

        f->cd();
        TF1* f2 = gMeanU2->GetFunction("pol1");
        f2->SetName("f_meanU2");
        f2->Write();
    }

    if(gSigU1)
    {
        TCanvas *c3 = new TCanvas("c_sigU1","Sigma U1", 800, 600);
        gSigU1->SetTitle(";Z p_{T} [GeV];#sigma(U_{1}) [GeV]");
        gSigU1->SetMarkerStyle(kOpenCircle);
        gSigU1->SetMarkerSize(0.9);
        gSigU1->SetLineWidth(1);
        gSigU1->Fit("pol2");
        gSigU1->Draw("AP");

        f->cd();
        TF1* f3 = gSigU1->GetFunction("pol2");
        f3->SetName("f_sigU1");
        f3->Write();
    }

    if(gSigU2)
    {
        TCanvas *c4 = new TCanvas("c_sigU2","Sigma U2", 800, 600);
        gSigU2->SetTitle(";Z p_{T} [GeV];#sigma(U_{2}) [GeV]");
        gSigU2->SetMarkerStyle(kOpenCircle);
        gSigU2->SetMarkerSize(0.9);
        gSigU2->SetLineWidth(1);
        gSigU2->Fit("pol2");
        gSigU2->Draw("AP");

        f->cd();
        TF1* f4 = gSigU2->GetFunction("pol2");
        f4->SetName("f_sigU2");
        f4->Write();
    }

    f->Close();
}
