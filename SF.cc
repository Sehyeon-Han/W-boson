void Muon_pt(const char* data = "/data2/kplee/Lecture/CMSOpenData/Data/SingleMuon/Run2016H/*.root",
             const char* WMC = "/data2/kplee/Lecture/CMSOpenData/MC2016/WJetsToLNu_aMCNLO/*.root",
             const char* DYMC = "/data2/kplee/Lecture/CMSOpenData/MC2016/DY_M50_aMCNLO/*.root",
             const char* TTMC = "/data2/kplee/Lecture/CMSOpenData/MC2016/TTTo2L2Nu_Powheg/*.root",
             double lumi_pb = 8740.1,
             double xsec_W = 61526.7,
             double xsec_DY = 6019.95,
             double xsec_TT = 88.51)
{
    TH1::SetDefaultSumw2();

    gStyle->SetOptStat(0);

    TChain chain_data("Events");
    chain_data.Add(data);
    TChain chain_WMC("Events");
    chain_WMC.Add(WMC);
    TChain chain_DYMC("Events");
    chain_DYMC.Add(DYMC);
    TChain chain_TTMC("Events");
    chain_TTMC.Add(TTMC);

    TH1D *h_data = new TH1D("h_data", "Muon_p_{T}; p_{T} [GeV]; Events", 50, 0, 200);
    TH1D *h_WMC = new TH1D("h_WMC", "Muon_p_{T}; p_{T} [GeV]; Events", 50, 0, 200);
    TH1D *h_DYMC = new TH1D("h_DYMC", "Muon_p_{T}; p_{T} [GeV]; Events", 50, 0, 200);
    TH1D *h_TTMC = new TH1D("h_TTMC", "Muon_p_{T}; p_{T} [GeV]; Events", 50, 0, 200);

    TFile *f_SF_ID = TFile::Open("EffSF_Muon_ID.root");
    TFile *f_SF_ISO = TFile::Open("EffSF_Muon_ISO.root");
    TFile *f_SF_TRIG = TFile::Open("EffSF_Muon_TRIG.root");

    TH2D* h2D_SF_ID = (TH2D*)f_SF_ID->Get("NUM_HighPtID_DEN_TrackerMuons_abseta_pt")->Clone();
    TH2D* h2D_SF_ISO = (TH2D*)f_SF_ISO->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt")->Clone();
    TH2D* h2D_eff_data_TRIG = (TH2D*)f_SF_TRIG->Get("NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt_efficiencyData")->Clone();
    TH2D* h2D_eff_mc_TRIG = (TH2D*)f_SF_TRIG->Get("NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt_efficiencyMC")->Clone();

    TH2D* h2D_SF_TRIG = (TH2D*)h2D_eff_data_TRIG->Clone("h2D_SF_TRIG");
    h2D_SF_TRIG->Divide(h2D_eff_mc_TRIG);

    {
        TTreeReader reader(&chain_data);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
        
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

            if(selCount == 1 && selIdx >= 0 && *HLT_IsoMu24)
            {
                h_data->Fill(Muon_pt[selIdx]);
            }
        }
    }

    double sumw_WMC = 0.0;
    {
        TTreeReader reader(&chain_WMC);
        
        TTreeReaderValue<float> genWeight(reader, "genWeight");

        while(reader.Next())
        {
            sumw_WMC += (double)(*genWeight);
        }
    }

    double sumw_DYMC = 0.0;
    {
        TTreeReader reader(&chain_DYMC);
        
        TTreeReaderValue<float> genWeight(reader, "genWeight");

        while(reader.Next())
        {
            sumw_DYMC += (double)(*genWeight);
        }
    }

    double sumw_TTMC = 0.0;
    {
        TTreeReader reader(&chain_TTMC);
        
        TTreeReaderValue<float> genWeight(reader, "genWeight");

        while(reader.Next())
        {
            sumw_TTMC += (double)(*genWeight);
        }
    }

    const double NFWMC = (xsec_W * lumi_pb) / sumw_WMC;
    const double NFDYMC = (xsec_DY * lumi_pb) / sumw_DYMC;
    const double NFTTMC = (xsec_TT * lumi_pb) / sumw_TTMC;

    {
        TTreeReader reader(&chain_WMC);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
        TTreeReaderValue<float> genWeight(reader, "genWeight");

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

            if(selCount == 1 && selIdx >=0 && *HLT_IsoMu24)
            {   
                double pt = Muon_pt[selIdx];
                double aeta = fabs(Muon_eta[selIdx]);

                int binx_id = h2D_SF_ID->GetXaxis()->FindBin(aeta);
                int biny_id = h2D_SF_ID->GetYaxis()->FindBin(pt);

                int binx_iso = h2D_SF_ISO->GetXaxis()->FindBin(aeta);
                int biny_iso = h2D_SF_ISO->GetYaxis()->FindBin(pt);

                int binx_trig = h2D_SF_TRIG->GetXaxis()->FindBin(aeta);
                int biny_trig = h2D_SF_TRIG->GetYaxis()->FindBin(pt);

                double sf_id = h2D_SF_ID->GetBinContent(binx_id, biny_id);
                double sf_iso = h2D_SF_ISO->GetBinContent(binx_iso, biny_iso);
                double sf_trig = h2D_SF_TRIG->GetBinContent(binx_trig, biny_trig);

                double sf = sf_id * sf_iso * sf_trig;

                double w = (double)(*genWeight) * NFWMC * sf;
                h_WMC->Fill(Muon_pt[selIdx], w);
            }
        }
    }

    {
        TTreeReader reader(&chain_DYMC);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
        TTreeReaderValue<float> genWeight(reader, "genWeight");

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

            if(selCount == 1 && selIdx >=0 && *HLT_IsoMu24)
            {
                double pt = Muon_pt[selIdx];
                double aeta = fabs(Muon_eta[selIdx]);

                int binx_id = h2D_SF_ID->GetXaxis()->FindBin(aeta);
                int biny_id = h2D_SF_ID->GetYaxis()->FindBin(pt);

                int binx_iso = h2D_SF_ISO->GetXaxis()->FindBin(aeta);
                int biny_iso = h2D_SF_ISO->GetYaxis()->FindBin(pt);

                int binx_trig = h2D_SF_TRIG->GetXaxis()->FindBin(aeta);
                int biny_trig = h2D_SF_TRIG->GetYaxis()->FindBin(pt);

                double sf_id = h2D_SF_ID->GetBinContent(binx_id, biny_id);
                double sf_iso = h2D_SF_ISO->GetBinContent(binx_iso, biny_iso);
                double sf_trig = h2D_SF_TRIG->GetBinContent(binx_trig, biny_trig);

                double sf = sf_id * sf_iso * sf_trig;

                double w = (double)(*genWeight) * NFDYMC * sf;
                h_DYMC->Fill(Muon_pt[selIdx], w);  
            }
        }
    }

    {
        TTreeReader reader(&chain_TTMC);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
        TTreeReaderValue<float> genWeight(reader, "genWeight");

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

            if(selCount == 1 && selIdx >=0 && *HLT_IsoMu24)
            {
                double pt = Muon_pt[selIdx];
                double aeta = fabs(Muon_eta[selIdx]);

                int binx_id = h2D_SF_ID->GetXaxis()->FindBin(aeta);
                int biny_id = h2D_SF_ID->GetYaxis()->FindBin(pt);

                int binx_iso = h2D_SF_ISO->GetXaxis()->FindBin(aeta);
                int biny_iso = h2D_SF_ISO->GetYaxis()->FindBin(pt);

                int binx_trig = h2D_SF_TRIG->GetXaxis()->FindBin(aeta);
                int biny_trig = h2D_SF_TRIG->GetYaxis()->FindBin(pt);

                double sf_id = h2D_SF_ID->GetBinContent(binx_id, biny_id);
                double sf_iso = h2D_SF_ISO->GetBinContent(binx_iso, biny_iso);
                double sf_trig = h2D_SF_TRIG->GetBinContent(binx_trig, biny_trig);

                double sf = sf_id * sf_iso * sf_trig;

                double w = (double)(*genWeight) * NFTTMC * sf;
                h_TTMC->Fill(Muon_pt[selIdx], w);
            }
        }
    }

    TH1D* hMCsum = (TH1D*)h_WMC->Clone("hMCsum");
    hMCsum->Add(h_DYMC);
    hMCsum->Add(h_TTMC);

    const double ymax = max(hMCsum->GetMaximum(), h_data->GetMaximum());
    TCanvas *c = new TCanvas("c", "Muon p_{T}", 900, 800);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.18, 1, 1);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    THStack *stack = new THStack("hstack", "Muon p_{T};p_{T} [GeV];Events");
    h_TTMC->SetFillColor(kGreen);
    h_DYMC->SetFillColor(kRed);
    h_WMC->SetFillColor(kOrange);

    stack->Add(h_TTMC);
    stack->Add(h_DYMC);
    stack->Add(h_WMC);

    stack->Draw("HIST");
    stack->SetMinimum(0.0);
    stack->SetMaximum(1.25 * ymax);

    h_data->SetLineColor(kBlack);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1.0);
    h_data->Draw("SAME");

    stack->GetXaxis()->SetLabelSize(0);
    stack->GetXaxis()->SetTitleSize(0);
    stack->GetXaxis()->SetTickLength(0);

    auto leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->AddEntry(h_data, "Data", "lep");
    leg->AddEntry(h_WMC,   "W#rightarrow#mu#nu", "f");
    leg->AddEntry(h_DYMC,  "DY (Z/#gamma*#rightarrow ll)", "f");
    leg->AddEntry(h_TTMC,  "t#bar{t}", "f");
    leg->Draw();

    pad1->Update();

    gStyle->SetOptStat(0);

    c->cd();

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.18);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.20);
    pad2->Draw();
    pad2->cd();

    TH1D *h_ratio = (TH1D*)h_data->Clone("h_ratio");

    h_ratio->SetTitle("");
    h_ratio->Divide(hMCsum);

    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerSize(0.8);
    h_ratio->SetLineColor(kBlack);

    h_ratio->GetYaxis()->SetTitle("Data / MC");
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.10);
    h_ratio->GetYaxis()->SetLabelSize(0.08);

    h_ratio->GetXaxis()->SetTitle("mT [GeV]");
    h_ratio->GetXaxis()->SetTitleSize(0.10);
    h_ratio->GetXaxis()->SetLabelSize(0.08);

    h_ratio->SetMinimum(0.0);
    h_ratio->SetMaximum(2.0);
    h_ratio->Draw();

    TLine *line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw("SAME");

    pad2->Update();

    c->Update();

}
