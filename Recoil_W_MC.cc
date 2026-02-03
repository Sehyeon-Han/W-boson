void Recoil_W_MC(const char* inFile  = "/data2/kplee/Lecture/CMSOpenData/MC2016/WJetsToLNu_aMCNLO/*.root",
                 const char* outFile = "recoil_W_MC.root",
                 double ptMin = 0,
                 double ptMax = 100,
                 double dPt   = 2)
{
    gStyle->SetOptFit(1110);
    gStyle->SetOptStat(1);

    TChain chain("Events");
    chain.Add(inFile);

    TTreeReader reader(&chain);

    TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
    TTreeReaderArray<bool>  Muon_tightId(reader, "Muon_tightId");
    TTreeReaderValue<float> MET_pt(reader, "MET_pt");
    TTreeReaderValue<float> MET_phi(reader, "MET_phi");
    TTreeReaderValue<bool>  HLT_IsoMu24(reader, "HLT_IsoMu24");
    TTreeReaderArray<int>   GenPart_pdgId(reader, "GenPart_pdgId");
    TTreeReaderArray<float> GenPart_pt(reader,   "GenPart_pt");
    TTreeReaderArray<float> GenPart_phi(reader,  "GenPart_phi");
    TTreeReaderValue<float> genWeight(reader, "genWeight");
    TTreeReaderArray<int> GenPart_statusFlags(reader, "GenPart_statusFlags");

    const int nBins = (int)((ptMax - ptMin) / dPt);

    std::vector<TH1D*> hU1(nBins, nullptr), hU2(nBins, nullptr);

    for(int i = 0; i < nBins; ++i)
    {
        double low  = ptMin + i * dPt;
        double high = low + dPt;
        hU1[i] = new TH1D(Form("hU1_%d", i), Form("U1 %.0f-%.0f;U1 [GeV];Events", low, high), 240, -120, 120);
        hU2[i] = new TH1D(Form("hU2_%d", i), Form("U2 %.0f-%.0f;U2 [GeV];Events", low, high), 240, -120, 120);
        hU1[i]->Sumw2();
        hU2[i]->Sumw2();
    }

    std::vector<TF1*> fitU1(nBins, nullptr), fitU2(nBins, nullptr);
    std::vector<double> x(nBins, 0.0);
    std::vector<double> mean1(nBins, 0.0), emean1(nBins, 0.0), sig1(nBins, 0.0), esig1(nBins, 0.0);
    std::vector<double> mean2(nBins, 0.0), emean2(nBins, 0.0), sig2(nBins, 0.0), esig2(nBins, 0.0);

    while(reader.Next())
    {
        if(!(*HLT_IsoMu24)) continue;

        int selCount = 0;
        int selIdx = -1;
        const int nMu = Muon_pt.GetSize();

        for(int i = 0; i < nMu; ++i)
        {
            if(Muon_pt[i] > 25.0 && std::fabs(Muon_eta[i]) < 2.4 && Muon_iso[i] < 0.12 && Muon_tightId[i])
            {
                ++selCount;
                selIdx = i;
                if(selCount > 1) break;
            }
        }
        if(!(selCount == 1 && selIdx >= 0)) continue;

        int wIdx = -1;
        float wMaxPt = -1.f;
        const int nGen = GenPart_pdgId.GetSize();
        for(int ig = 0; ig < nGen; ++ig)
        {
            if(abs(GenPart_pdgId[ig]) != 24) continue;
            if(!(abs(GenPart_statusFlags[ig]) & 1)) continue;

            if(GenPart_pt[ig] > wMaxPt)
            {
                wMaxPt = GenPart_pt[ig];
                wIdx = ig;
            }
        }
        if(wIdx < 0) continue;

        const double ptW_gen  = GenPart_pt[wIdx];
        const double phiW_gen = GenPart_phi[wIdx];

        if(ptW_gen < ptMin || ptW_gen >= ptMax) continue;
        const int ibin = (int)((ptW_gen - ptMin) / dPt);
        if(ibin < 0 || ibin >= nBins) continue;

        const double muPt  = Muon_pt[selIdx];
        const double muPhi = Muon_phi[selIdx];
        const double muPx  = muPt * std::cos(muPhi);
        const double muPy  = muPt * std::sin(muPhi);

        const double metPx = (*MET_pt) * std::cos(*MET_phi);
        const double metPy = (*MET_pt) * std::sin(*MET_phi);

        const double Ux = -(metPx + muPx);
        const double Uy = -(metPy + muPy);

        const double px = std::cos(phiW_gen);
        const double py = std::sin(phiW_gen);
        const double nx = -std::sin(phiW_gen);
        const double ny =  std::cos(phiW_gen);

        const double U1 = Ux * px + Uy * py;
        const double U2 = Ux * nx + Uy * ny;

        hU1[ibin]->Fill(U1, *genWeight);
        hU2[ibin]->Fill(U2, *genWeight);
    }

    for(int i = 0; i < nBins; ++i)
    {
        x[i] = ptMin + (i + 0.5) * dPt;

        if(hU1[i]->GetEffectiveEntries() > 50)
        {
            double m = hU1[i]->GetMean();
            double r = hU1[i]->GetRMS();
            if(r > 1e-6)
            {
                double lo = hU1[i]->GetXaxis()->GetXmin();
                double hi = hU1[i]->GetXaxis()->GetXmax();
                fitU1[i] = new TF1(Form("fitU1_%d", i), "gaus", lo, hi);
                fitU1[i]->SetParameters(hU1[i]->GetMaximum(), m, r);
                hU1[i]->Fit(fitU1[i], "RQ0");
                mean1[i]  = fitU1[i]->GetParameter(1);
                emean1[i] = fitU1[i]->GetParError(1);
                sig1[i]   = std::fabs(fitU1[i]->GetParameter(2));
                esig1[i]  = fitU1[i]->GetParError(2);
            }
        }

        if(hU2[i]->GetEffectiveEntries() > 50)
        {
            double m = hU2[i]->GetMean();
            double r = hU2[i]->GetRMS();
            if(r > 1e-6)
            {
                double lo = hU2[i]->GetXaxis()->GetXmin();
                double hi = hU2[i]->GetXaxis()->GetXmax();
                fitU2[i] = new TF1(Form("fitU2_%d", i), "gaus", lo, hi);
                fitU2[i]->SetParameters(hU2[i]->GetMaximum(), m, r);
                hU2[i]->Fit(fitU2[i], "RQ0");
                mean2[i]  = fitU2[i]->GetParameter(1);
                emean2[i] = fitU2[i]->GetParError(1);
                sig2[i]   = std::fabs(fitU2[i]->GetParameter(2));
                esig2[i]  = fitU2[i]->GetParError(2);
            }
        }
    }

    std::vector<double> xx1, exx1, my1, emy1, sy1, esy1;
    std::vector<double> xx2, exx2, my2, emy2, sy2, esy2;

    for(int i = 0; i < nBins; ++i)
    {
        if(fitU1[i])
        {
            xx1.push_back(x[i]); exx1.push_back(0.0);
            my1.push_back(mean1[i]); emy1.push_back(emean1[i]);
            sy1.push_back(sig1[i]);  esy1.push_back(esig1[i]);
        }
        if(fitU2[i])
        {
            xx2.push_back(x[i]); exx2.push_back(0.0);
            my2.push_back(mean2[i]); emy2.push_back(emean2[i]);
            sy2.push_back(sig2[i]);  esy2.push_back(esig2[i]);
        }
    }

    auto gMeanU1 = new TGraphErrors((int)xx1.size(), xx1.data(), my1.data(), exx1.data(), emy1.data());
    auto gSigU1  = new TGraphErrors((int)xx1.size(), xx1.data(), sy1.data(), exx1.data(), esy1.data());
    auto gMeanU2 = new TGraphErrors((int)xx2.size(), xx2.data(), my2.data(), exx2.data(), emy2.data());
    auto gSigU2  = new TGraphErrors((int)xx2.size(), xx2.data(), sy2.data(), exx2.data(), esy2.data());

    gMeanU1->SetName("g_mean_U1");
    gMeanU1->SetTitle(";W p_{T} [GeV];#LT U_{1} #GT [GeV]");
    gSigU1->SetName("g_sigma_U1");
    gSigU1->SetTitle(";W p_{T} [GeV];#sigma(U_{1}) [GeV]");
    gMeanU2->SetName("g_mean_U2");
    gMeanU2->SetTitle(";W p_{T} [GeV];#LT U_{2} #GT [GeV]");
    gSigU2->SetName("g_sigma_U2");
    gSigU2->SetTitle(";W p_{T} [GeV];#sigma(U_{2}) [GeV]");

    TF1 f_meanU1("f_meanU1", "pol1", ptMin, ptMax);
    TF1 f_meanU2("f_meanU2", "pol1", ptMin, ptMax);
    TF1 f_sigU1 ("f_sigU1",  "pol2", ptMin, ptMax);
    TF1 f_sigU2 ("f_sigU2",  "pol2", ptMin, ptMax);

    if(gMeanU1->GetN() >= 3) gMeanU1->Fit(&f_meanU1, "Q0R");
    if(gMeanU2->GetN() >= 3) gMeanU2->Fit(&f_meanU2, "Q0R");
    if(gSigU1->GetN()  >= 4) gSigU1 ->Fit(&f_sigU1,  "Q0R");
    if(gSigU2->GetN()  >= 4) gSigU2 ->Fit(&f_sigU2,  "Q0R");

    gMeanU1->SetMarkerStyle(kOpenCircle);
    gMeanU1->SetMarkerSize(0.9);
    gMeanU1->SetLineWidth(1);

    gMeanU2->SetMarkerStyle(kOpenCircle);
    gMeanU2->SetMarkerSize(0.9);
    gMeanU2->SetLineWidth(1);

    gSigU1->SetMarkerStyle(kOpenCircle);
    gSigU1->SetMarkerSize(0.9);
    gSigU1->SetLineWidth(1);

    gSigU2->SetMarkerStyle(kOpenCircle);
    gSigU2->SetMarkerSize(0.9);
    gSigU2->SetLineWidth(1);

    f_meanU1.SetLineColor(kRed);
    f_meanU1.SetLineWidth(2);
    f_meanU2.SetLineColor(kRed);
    f_meanU2.SetLineWidth(2);
    f_sigU1.SetLineColor(kRed);
    f_sigU1.SetLineWidth(2);
    f_sigU2.SetLineColor(kRed);
    f_sigU2.SetLineWidth(2);

    TCanvas *cMean = new TCanvas("cMean", "Mean(U1), Mean(U2)", 1200, 450);
    cMean->Divide(2, 1);

    cMean->cd(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);
    gMeanU1->Draw("AP");
    gMeanU1->Fit(&f_meanU1, "R");
    gPad->Update();

    cMean->cd(2);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);
    gMeanU2->Draw("AP");
    gMeanU2->Fit(&f_meanU2, "R");  
    gPad->Update();

    cMean->SaveAs("MeanU1U2_fit.png");

    TCanvas *cSig = new TCanvas("cSig", "Sigma(U1), Sigma(U2)", 1200, 450);
    cSig->Divide(2, 1);

    cSig->cd(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);
    gSigU1->Draw("AP");
    gSigU1->Fit(&f_sigU1, "R");  
    gPad->Update();

    cSig->cd(2);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);
    gSigU2->Draw("AP");
    gSigU2->Fit(&f_sigU2, "R");
    gPad->Update();

    cSig->SaveAs("SigmaU1U2_fit.png");

    TFile fout(outFile, "RECREATE");
    for(int i = 0; i < nBins; ++i)
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

    f_meanU1.Write("f_meanU1");
    f_meanU2.Write("f_meanU2");
    f_sigU1.Write("f_sigU1");
    f_sigU2.Write("f_sigU2");

    fout.Close();
}
