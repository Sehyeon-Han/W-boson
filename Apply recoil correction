void Recoil_correction(const char* WMC  = "/data2/kplee/Lecture/CMSOpenData/MC2016/WJetsToLNu_aMCNLO/*.root",
                       const char* fitZMC = "recoil_Z_MC.root",
                       const char* fitZData = "recoil_Z_data.root",
                       const char* fitWMC = "recoil_W_MC.root",
                       const char* outFile = "WMC_recoil_correction.root")
{
    gStyle->SetOptFit(1110);
    gStyle->SetOptStat(1);

    TFile fzmc(fitZMC, "READ");
    TFile fzdata(fitZData, "READ");
    TFile fwmc(fitWMC, "READ");

    auto zmc_meanU1 = (TF1*)fzmc.Get("f_meanU1");
    auto zmc_meanU2 = (TF1*)fzmc.Get("f_meanU2");
    auto zmc_sigU1 = (TF1*)fzmc.Get("f_sigU1");
    auto zmc_sigU2 = (TF1*)fzmc.Get("f_sigU2");

    auto zdata_meanU1 = (TF1*)fzdata.Get("f_meanU1");
    auto zdata_meanU2 = (TF1*)fzdata.Get("f_meanU2");
    auto zdata_sigU1 = (TF1*)fzdata.Get("f_sigU1");
    auto zdata_sigU2 = (TF1*)fzdata.Get("f_sigU2");

    auto wmc_meanU1 = (TF1*)fwmc.Get("f_meanU1");
    auto wmc_meanU2 = (TF1*)fwmc.Get("f_meanU2");
    auto wmc_sigU1 = (TF1*)fwmc.Get("f_sigU1");
    auto wmc_sigU2 = (TF1*)fwmc.Get("f_sigU2");

    TChain chain("Events");
    chain.Add(WMC);

    TTreeReader reader(&chain);

    TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
    TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
    TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
    TTreeReaderValue<float> MET_pt(reader, "MET_pt");
    TTreeReaderValue<float> MET_phi(reader, "MET_phi");
    TTreeReaderArray<int> GenPart_pdgId(reader, "GenPart_pdgId");
    TTreeReaderArray<float> GenPart_pt(reader, "GenPart_pt");
    TTreeReaderArray<float> GenPart_phi(reader, "GenPart_phi");
    TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
    TTreeReaderValue<float> genWeight(reader, "genWeight");
    TTreeReaderArray<int> GenPart_statusFlags(reader, "GenPart_statusFlags");
    TTreeReaderValue<UInt_t> run(reader, "run");
    TTreeReaderValue<UInt_t> luminosityBlock(reader, "luminosityBlock");
    TTreeReaderValue<ULong64_t> event(reader, "event");

    TFile fout(outFile, "RECREATE");
    TTree out("Events", "W MC with recoil correction");

    float U1_old = 0;
    float U2_old = 0;
    float U1_new = 0;
    float U2_new = 0;
    float MET_corr_pt = 0;
    float MET_corr_phi = 0;
    float ptW_gen = 0;
    float mT_old = 0;
    float mT_new = 0;
    UInt_t out_run = 0;
    UInt_t out_lumi = 0;
    ULong64_t out_event = 0;
    float out_genWeight = 0;

    out.Branch("U1_old", &U1_old, "U1_old/F");
    out.Branch("U2_old", &U2_old, "U2_old/F");
    out.Branch("U1_new", &U1_new, "U1_new/F");
    out.Branch("U2_new", &U2_new, "U2_new/F");
    out.Branch("MET_corr_pt", &MET_corr_pt, "MET_corr_pt/F");
    out.Branch("MET_corr_phi", &MET_corr_phi, "MET_corr_phi/F");
    out.Branch("ptW_gen", &ptW_gen, "ptW_gen/F");
    out.Branch("mT_old", &mT_old, "mT_old/F");
    out.Branch("mT_new", &mT_new, "mT_new/F");
    out.Branch("run", &out_run, "run/i");
    out.Branch("luminosityBlock", &out_lumi, "luminosityBlock/i");
    out.Branch("event", &out_event, "event/l");
    out.Branch("genWeight", &out_genWeight, "genWeight/F");

    const double ptMin = 0.0;
    const double ptMax = 100;
    const double dPt = 2.0;
    const int nBin = (ptMax - ptMin) / dPt;

    vector<TH1D*> hU1(nBin, nullptr);
    vector<TH1D*> hU2(nBin, nullptr);

    for(int i = 0; i < nBin; ++i)
    {
        double lo = ptMin + i * dPt;
        double hi = lo + dPt;

        hU1[i] = new TH1D(Form("hU1_corr_bin%d", i), Form("U1 corr (%.0f-%.0f GeV);U1;Events", lo, hi), 200, -200, 200);
        hU1[i]->Sumw2();
        hU2[i] = new TH1D(Form("hU2_corr_bin%d", i), Form("U2 corr (%.0f-%.0f GeV);U2;Events", lo, hi), 200, -200, 200);
        hU2[i]->Sumw2();
    }
    
    TH1::AddDirectory(kFALSE);

    TH1D* h_mT_old = new TH1D("h_mT_old", "m_{T} (old MET);m_{T} [GeV];Events", 200, 0, 200);
    TH1D* h_mT_new = new TH1D("h_mT_new", "m_{T} (recoil-corrected MET);m_{T} [GeV];Events", 200, 0, 200);
    TH1D* h_METcorr = new TH1D("h_METcorr", "MET corr;MET [GeV];Events", 50, 0, 200);
    TH1D* h_METtrue = new TH1D("h_METtrue", "MET true (MC);MET [GeV];Events", 50, 0, 200);

    h_mT_old->Sumw2();
    h_mT_new->Sumw2();
    h_METcorr->Sumw2();
    h_METtrue->Sumw2();

    while(reader.Next())
    {
        int selCount = 0;
        int selIdx = -1;

        const int n = Muon_pt.GetSize();
        for(int i = 0; i < n; ++i)
        {
            if(Muon_pt[i] > 25.0 && fabs(Muon_eta[i]) < 2.4 && Muon_iso[i] < 0.12 && Muon_tightId[i])
            {
                ++selCount;
                selIdx = i;

                if(selCount > 1)
                    break;
            }
        }

        if(!(selCount == 1 && selIdx >= 0 && *HLT_IsoMu24)) continue;
        const int imu = selIdx;
        
        int wIdx = -1;
        float wMaxPt = -1.0f;

        const int nGen = GenPart_pdgId.GetSize();
        for(int i = 0; i < nGen; ++i)
        {
            if(abs(GenPart_pdgId[i]) != 24) continue;
            if(!(abs(GenPart_statusFlags[i]) & 1)) continue;

            if(GenPart_pt[i] > wMaxPt)
            {
                wMaxPt = GenPart_pt[i];
                wIdx = i;
            }
        }

        if(wIdx < 0) continue;
        const int iW = wIdx;

        ptW_gen = GenPart_pt[iW];
        float phiW = GenPart_phi[iW];

        float lep_px = Muon_pt[imu] * cos(Muon_phi[imu]);
        float lep_py = Muon_pt[imu] * sin(Muon_phi[imu]);

        float met_px = (*MET_pt) * cos(*MET_phi);
        float met_py = (*MET_pt) * sin(*MET_phi);

        float Ux = -(met_px + lep_px);
        float Uy = -(met_py + lep_py);
        
        float px = cos(phiW);
        float py = sin(phiW);
        float nx = -sin(phiW);
        float ny = cos(phiW);

        U1_old = Ux * px + Uy * py;
        U2_old = Ux * nx + Uy * ny;

        float zmc_m1 = zmc_meanU1->Eval(ptW_gen);
        float zmc_m2 = zmc_meanU2->Eval(ptW_gen);
        float zmc_s1 = zmc_sigU1->Eval(ptW_gen);
        float zmc_s2 = zmc_sigU2->Eval(ptW_gen);

        float zd_m1 = zdata_meanU1->Eval(ptW_gen);
        float zd_m2 = zdata_meanU2->Eval(ptW_gen);
        float zd_s1 = zdata_sigU1->Eval(ptW_gen);
        float zd_s2 = zdata_sigU2->Eval(ptW_gen);

        float k_meanU1 = zd_m1 / zmc_m1;
        float k_meanU2 = zd_m2 / zmc_m2;
        float k_sigU1 = zd_s1 / zmc_s1;
        float k_sigU2 = zd_s2 / zmc_s2;

        float meanU1_corr = k_meanU1 * wmc_meanU1->Eval(ptW_gen);
        float meanU2_corr = k_meanU2 * wmc_meanU2->Eval(ptW_gen);
        float sigU1_corr = k_sigU1 * wmc_sigU1->Eval(ptW_gen);
        float sigU2_corr = k_sigU2 * wmc_sigU2->Eval(ptW_gen);

        const float wmc_m1 = wmc_meanU1->Eval(ptW_gen);
        const float wmc_m2 = wmc_meanU2->Eval(ptW_gen);
        const float wmc_s1 = wmc_sigU1 ->Eval(ptW_gen);
        const float wmc_s2 = wmc_sigU2 ->Eval(ptW_gen);

        float r1 = 1.0;
        float r2 = 1.0;
        r1 = sigU1_corr / wmc_s1;
        r2 = sigU2_corr / wmc_s2;

        U1_new = meanU1_corr + (U1_old - wmc_m1) * r1;
        U2_new = meanU2_corr + (U2_old - wmc_m2) * r2;

        int i = (ptW_gen - ptMin) / dPt;

        if(i >= 0 && i < nBin)
        {
            hU1[i]->Fill(U1_new, *genWeight);
            hU2[i]->Fill(U2_new, *genWeight);
        }

        float u_corr_px = (U1_new) * px + (U2_new) * nx;
        float u_corr_py = (U1_new) * py + (U2_new) * ny;

        float metc_px = -(u_corr_px + lep_px);
        float metc_py = -(u_corr_py + lep_py);

        MET_corr_pt = sqrt(metc_px * metc_px + metc_py * metc_py);
        MET_corr_phi = atan2(metc_py, metc_px);

        const float dphi_old = TVector2::Phi_mpi_pi(Muon_phi[imu] - (*MET_phi));
        mT_old = sqrt(2 * Muon_pt[imu] * (*MET_pt) * (1 - cos(dphi_old)));

        const float dphi_new = TVector2::Phi_mpi_pi(Muon_phi[imu] - MET_corr_phi);
        mT_new = sqrt(2 * Muon_pt[imu] * MET_corr_pt * (1 - cos(dphi_new)));

        h_mT_old->Fill(mT_old, *genWeight);
        h_mT_new->Fill(mT_new, *genWeight);
        h_METtrue->Fill(*MET_pt, *genWeight);
        h_METcorr->Fill(MET_corr_pt, *genWeight);
        
        out_run = *run;
        out_lumi = *luminosityBlock;
        out_event = *event;
        out_genWeight = *genWeight;

        out.Fill();
    }

    vector<double> x(nBin), ex(nBin, 0.0);
    vector<double> m1(nBin), em1(nBin), s1(nBin), es1(nBin);
    vector<double> m2(nBin), em2(nBin), s2(nBin), es2(nBin);

    for(int i = 0; i < nBin; ++i)
    {
        double lo = ptMin + i * dPt;
        double hi = lo + dPt;
        x[i] = 0.5 * (lo + hi);

        if (hU1[i]->GetEntries() < 200 || hU2[i]->GetEntries() < 200)
        {
            m1[i] = m2[i] = 0;
            s1[i] = s2[i] = 0;
            em1[i] = em2[i] = 0;
            es1[i] = es2[i] = 0;
            continue;
        }

        TF1 g1("g1", "gaus", hU1[i]->GetMean() - 2 * hU1[i]->GetRMS(), hU1[i]->GetMean() + 2 * hU1[i]->GetRMS());
        hU1[i]->Fit(&g1, "Q0R");

        m1[i] = g1.GetParameter(1);
        s1[i] = fabs(g1.GetParameter(2));
        em1[i] = g1.GetParError(1);
        es1[i] = g1.GetParError(2);

        TF1 g2("g2", "gaus", hU2[i]->GetMean() - 2 * hU2[i]->GetRMS(), hU2[i]->GetMean() + 2 * hU2[i]->GetRMS());
        hU2[i]->Fit(&g2, "Q0R");

        m2[i] = g2.GetParameter(1);
        s2[i] = fabs(g2.GetParameter(2));
        em2[i] = g2.GetParError(1);
        es2[i] = g2.GetParError(2);
    }

    TGraphErrors gr_mU1(nBin, x.data(), m1.data(), ex.data(), em1.data());
    TGraphErrors gr_mU2(nBin, x.data(), m2.data(), ex.data(), em2.data());
    TGraphErrors gr_sU1(nBin, x.data(), s1.data(), ex.data(), es1.data());
    TGraphErrors gr_sU2(nBin, x.data(), s2.data(), ex.data(), es2.data());

    gr_mU1.SetName("gr_meanU1_corr");
    gr_mU1.SetTitle(";W p_{T} [GeV];#LT U_{1} #GT [GeV]");
    gr_mU2.SetName("gr_meanU2_corr");
    gr_mU2.SetTitle(";W p_{T} [GeV];#LT U_{2} #GT [GeV]");
    gr_sU1.SetName("gr_sigU1_corr");
    gr_sU1.SetTitle(";W p_{T} [GeV];#sigma(U_{1}) [GeV]");
    gr_sU2.SetName("gr_sigU2_corr");    
    gr_sU2.SetTitle(";W p_{T} [GeV];#sigma(U_{2}) [GeV]");

    gr_mU1.SetMarkerStyle(kOpenCircle);
    gr_mU1.SetMarkerSize(0.9);
    gr_mU1.SetLineWidth(1);

    gr_mU2.SetMarkerStyle(kOpenCircle);
    gr_mU2.SetMarkerSize(0.9);
    gr_mU2.SetLineWidth(1);

    gr_sU1.SetMarkerStyle(kOpenCircle);
    gr_sU1.SetMarkerSize(0.9);
    gr_sU1.SetLineWidth(1);

    gr_sU2.SetMarkerStyle(kOpenCircle);
    gr_sU2.SetMarkerSize(0.9);
    gr_sU2.SetLineWidth(1);

    TF1 f_meanU1("f_meanU1","pol1", ptMin, ptMax);
    TF1 f_meanU2("f_meanU2","pol1", ptMin, ptMax);
    TF1 f_sigU1 ("f_sigU1" ,"pol2", ptMin, ptMax);
    TF1 f_sigU2 ("f_sigU2" ,"pol2", ptMin, ptMax);      

    gr_mU1.Fit(&f_meanU1,"Q0R");
    gr_mU2.Fit(&f_meanU2,"Q0R");
    gr_sU1.Fit(&f_sigU1 ,"Q0R");
    gr_sU2.Fit(&f_sigU2 ,"Q0R");

    TH1D* h_frac = (TH1D*)h_METcorr->Clone("h_MET_fracDiff");
    h_frac->Add(h_METtrue, -1.0);     // pred - true
    h_frac->Divide(h_METtrue);        // (pred-true)/true
    h_frac->SetTitle("Frac Difference;MET [GeV];(Pred-True)/True");

    TCanvas* cMET = new TCanvas("c_MET_closure", "MET closure", 900, 800);
    TPad* p1 = new TPad("p1", "", 0, 0.20, 1, 1.00);
    TPad* p2 = new TPad("p2", "", 0, 0.00, 1, 0.20);

    p1->SetBottomMargin(0.12);
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.30);
    p1->Draw();
    p2->Draw();

    p1->cd();
    h_METtrue->SetLineWidth(2);
    h_METtrue->SetLineColor(kRed);
    h_METtrue->Draw("HIST");

    h_METcorr->SetLineWidth(2);
    h_METcorr->SetLineColor(kBlue);
    h_METcorr->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.58, 0.75, 0.88, 0.88);
    leg->AddEntry(h_METtrue, "Uncorrected", "l");
    leg->AddEntry(h_METcorr, "Corrected", "l");
    leg->SetBorderSize(0);
    leg->Draw();

    p2->cd();
    h_frac->SetTitle("");
    h_frac->SetMarkerStyle(20);
    h_frac->SetMarkerSize(0.8);
    h_frac->GetYaxis()->SetTitle("(Corr - Uncorr) / Uncorr");
    h_frac->GetYaxis()->SetTitleSize(0.10);
    h_frac->GetYaxis()->SetTitleOffset(0.45);
    h_frac->GetYaxis()->SetLabelSize(0.10);
    h_frac->GetXaxis()->SetTitleSize(0.12);
    h_frac->GetXaxis()->SetLabelSize(0.10);
    h_frac->GetYaxis()->SetNdivisions(505);
    h_frac->SetMinimum(-0.5);
    h_frac->SetMaximum(0.5);
    h_frac->Draw("E1");

    h_METtrue->SetStats(0);
    h_METcorr->SetStats(0);
    h_frac->SetStats(0);

    TLine* l0 = new TLine(h_frac->GetXaxis()->GetXmin(), 0.0, h_frac->GetXaxis()->GetXmax(), 0.0);

    l0->SetLineStyle(2);
    l0->Draw("SAME");

    fout.cd();

    out.Write();

    gr_mU1.Write();
    gr_mU2.Write();
    gr_sU1.Write();
    gr_sU2.Write();
    f_meanU1.Write("f_meanU1");
    f_meanU2.Write("f_meanU2");
    f_sigU1.Write("f_sigU1");
    f_sigU2.Write("f_sigU2");

    h_mT_old->Write();
    h_mT_new->Write();
    h_METcorr->Write();
    h_METtrue->Write();

    cMET->Modified();
    cMET->Update();
    cMET->Write();
    cMET->SaveAs("MET_closure.png");

    TCanvas *cMean = new TCanvas("cMean", "Mean(U1), Mean(U2)", 1200, 450);
    cMean->Divide(2, 1);

    cMean->cd(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);
    gr_mU1.Draw("AP");
    gr_mU1.Fit(&f_meanU1, "R");
    gPad->Update();

    cMean->cd(2);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);
    gr_mU2.Draw("AP");
    gr_mU2.Fit(&f_meanU2, "R");  
    gPad->Update();

    cMean->SaveAs("MeanU1U2_fit.png");

    TCanvas *cSig = new TCanvas("cSig", "Sigma(U1), Sigma(U2)", 1200, 450);
    cSig->Divide(2, 1);

    cSig->cd(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);
    gr_sU1.Draw("AP");
    gr_sU1.Fit(&f_sigU1, "R");  
    gPad->Update();

    cSig->cd(2);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.12);
    gr_sU2.Draw("AP");
    gr_sU2.Fit(&f_sigU2, "R");
    gPad->Update();

    fout.Write();
    fout.Close();
}
