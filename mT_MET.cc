void mT_MET(const char* data = "/data2/kplee/Lecture/CMSOpenData/Data/SingleMuon/Run2016H/*.root",
            const char* WMC = "/data2/kplee/Lecture/CMSOpenData/MC2016/WJetsToLNu_aMCNLO/*.root",
            const char* DYMC = "/data2/kplee/Lecture/CMSOpenData/MC2016/DY_M50_aMCNLO/*.root",
            const char* TTMC = "/data2/kplee/Lecture/CMSOpenData/MC2016/TTTo2L2Nu_Powheg/*.root",
             
            const char* WrecoilFile = "WMC_recoil_correction.root",
            const char* hWmTName = "h_mT_new",
            const char* hWMETName = "h_METcorr",
             
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
    TChain chain_DY("Events");
    chain_DY.Add(DYMC);
    TChain chain_TT("Events");
    chain_TT.Add(TTMC);

    TH1D *h_data_mT = new TH1D("h_data_mT", "Transverse Mass;m_{T} [GeV];Events", 50, 0, 200);
    TH1D *h_DY_mT = new TH1D("h_DY_mT", "Transverse Mass;m_{T} [GeV];Events", 50, 0, 200);
    TH1D *h_TT_mT = new TH1D("h_TT_mT", "Transverse Mass;m_{T} [GeV];Events", 50, 0, 200);

    TH1D *h_data_MET = new TH1D("h_data_MET", "MET;MET [GeV];Events", 50, 0, 200);
    TH1D *h_DY_MET = new TH1D("h_DY_MET", "MET;MET [GeV];Events", 50, 0, 200);
    TH1D *h_TT_MET = new TH1D("h_TT_MET", "MET;MET [GeV];Events", 50, 0, 200);

    {
        TTreeReader reader(&chain_data);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
        TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
        TTreeReaderValue<float> MET_pt(reader, "MET_pt");
        TTreeReaderValue<float> MET_phi(reader, "MET_phi");

        while(reader.Next())
        {
            if(!(*HLT_IsoMu24)) continue;

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

            if(!(selCount == 1 && selIdx >= 0)) continue;

            float dphi = TVector2::Phi_mpi_pi(Muon_phi[selIdx] - (*MET_phi));
            float mT = sqrt(2 * Muon_pt[selIdx] * (*MET_pt) * (1 - cos(dphi)));
            h_data_mT->Fill(mT);
        }
    }

    {
        TTreeReader reader(&chain_data);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
        TTreeReaderValue<float> MET_pt(reader, "MET_pt");

        while(reader.Next())
        {
            if(!(*HLT_IsoMu24)) continue;

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

            if(!(selCount == 1 && selIdx >= 0)) continue;

            h_data_MET->Fill(*MET_pt);
        }
    }

    {
        TTreeReader reader(&chain_DY);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
        TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
        TTreeReaderValue<float> genWeight(reader, "genWeight");
        TTreeReaderValue<float> MET_pt(reader, "MET_pt");
        TTreeReaderValue<float> MET_phi(reader, "MET_phi");

        while(reader.Next())
        {
            const double w = (double)(*genWeight);

            if(!(*HLT_IsoMu24)) continue;

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

            if(!(selCount == 1 && selIdx >= 0)) continue;

            float dphi = TVector2::Phi_mpi_pi(Muon_phi[selIdx] - (*MET_phi));
            float mT = sqrt(2 * Muon_pt[selIdx] * (*MET_pt) * (1 - cos(dphi)));
            h_DY_mT->Fill(mT, w);
        }
    }
    
    {
        TTreeReader reader(&chain_TT);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_phi(reader, "Muon_phi");
        TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<float> genWeight(reader, "genWeight");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
        TTreeReaderValue<float> MET_pt(reader, "MET_pt");
        TTreeReaderValue<float> MET_phi(reader, "MET_phi");

        while(reader.Next())
        {
            const double w = (double)(*genWeight);
            if(!(*HLT_IsoMu24)) continue;

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

            if(!(selCount == 1 && selIdx >= 0)) continue;

            float dphi = TVector2::Phi_mpi_pi(Muon_phi[selIdx] - (*MET_phi));
            float mT = sqrt(2 * Muon_pt[selIdx] * (*MET_pt) * (1 - cos(dphi)));
            h_TT_mT->Fill(mT, w);
        }
    }

    {
        TTreeReader reader(&chain_DY);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
        TTreeReaderValue<float> genWeight(reader, "genWeight");
        TTreeReaderValue<float> MET_pt(reader, "MET_pt");

        while(reader.Next())
        {
            const double w = (double)(*genWeight);

            if(!(*HLT_IsoMu24))continue;

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

            if(!(selCount == 1 && selIdx >= 0)) continue;

            h_DY_MET->Fill(*MET_pt, w);
        }
    }

    {
        TTreeReader reader(&chain_TT);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
        TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");
        TTreeReaderValue<float> genWeight(reader, "genWeight");
        TTreeReaderValue<float> MET_pt(reader, "MET_pt");

        while(reader.Next())
        {
            const double w = (double)(*genWeight);

            if(!(*HLT_IsoMu24)) continue;

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

            if(!(selCount == 1 && selIdx >= 0)) continue;

            h_TT_MET->Fill(*MET_pt, w);
        }
    }

    TFile fW(WrecoilFile, "READ");

    TH1D *hWmT_in = (TH1D*)fW.Get(hWmTName);
    TH1D *hWMET_in = (TH1D*)fW.Get(hWMETName);

    TH1D *h_W_mT = (TH1D*)hWmT_in->Clone("h_W_mT_recoil");
    TH1D *h_W_MET = (TH1D*)hWMET_in->Clone("h_W_MET_recoil");

    h_W_mT ->SetDirectory(0);
    h_W_MET->SetDirectory(0);
    fW.Close();

    double sW = 0.0, sDY = 0.0, sTT = 0.0;

    {
        TTreeReader reader(&chain_WMC);

        TTreeReaderValue<float> genWeight(reader, "genWeight");

        while(reader.Next())
        sW += (double)(*genWeight);
    }

    {
        TTreeReader reader(&chain_DY);

        TTreeReaderValue<float> genWeight(reader, "genWeight");

        while(reader.Next())
        sDY += (double)(*genWeight);
    }

    {
        TTreeReader reader(&chain_TT);

        TTreeReaderValue<float> genWeight(reader, "genWeight");

        while(reader.Next())
        sTT += (double)(*genWeight);
    }

    const double kW  = (sW  != 0.0) ? (lumi_pb * xsec_W ) / sW  : 0.0;
    const double kDY = (sDY != 0.0) ? (lumi_pb * xsec_DY) / sDY : 0.0;
    const double kTT = (sTT != 0.0) ? (lumi_pb * xsec_TT) / sTT : 0.0;

    h_W_mT ->Scale(kW);
    h_W_MET->Scale(kW);
    h_DY_mT ->Scale(kDY);
    h_TT_mT ->Scale(kTT);
    h_DY_MET->Scale(kDY);
    h_TT_MET->Scale(kTT);

    {
        TH1D *hMCsum = (TH1D*)h_W_mT->Clone("c_mT_MCsum");
        hMCsum->Add(h_DY_mT);
        hMCsum->Add(h_TT_mT);

        const double ymax = max(hMCsum->GetMaximum(), h_data_mT->GetMaximum());

        TCanvas *c = new TCanvas("c_mT", "c_mT", 700, 700);

        TPad *pad1 = new TPad("c_mT_pad1", "pad1", 0, 0.22, 1, 1);
        pad1->SetBottomMargin(0.02);
        pad1->Draw();
        pad1->cd();

        THStack *stack = new THStack("c_mT_stack", "");
        h_TT_mT->SetFillColor(kGreen);
        h_TT_mT->SetLineColor(kBlack);
        h_DY_mT->SetFillColor(kRed);
        h_DY_mT->SetLineColor(kBlack);
        h_W_mT ->SetFillColor(kOrange);
        h_W_mT ->SetLineColor(kBlack);

        stack->Add(h_TT_mT);
        stack->Add(h_DY_mT);
        stack->Add(h_W_mT);

        stack->Draw("HIST");
        stack->SetMinimum(0.0);
        stack->SetMaximum(1.25 * ymax);
        stack->GetXaxis()->SetLimits(0, 200);
        stack->GetXaxis()->SetLabelSize(0);
        stack->GetXaxis()->SetTitleSize(0);
        stack->GetYaxis()->SetTitle("Events");

        h_data_mT->SetMarkerStyle(20);
        h_data_mT->SetMarkerSize(1.0);
        h_data_mT->SetLineColor(kBlack);
        h_data_mT->Draw("E1 SAME");

        TLegend *leg = new TLegend(0.60, 0.64, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_data_mT, "Data", "lep");
        leg->AddEntry(h_W_mT,  "W#rightarrow#mu#nu (recoil)", "f");
        leg->AddEntry(h_DY_mT, "DY (Z/#gamma*#rightarrow ll)", "f");
        leg->AddEntry(h_TT_mT, "t#bar{t}", "f");
        leg->Draw();

        c->cd();
        TPad *pad2 = new TPad("c_mT_pad2", "pad2", 0, 0.00, 1, 0.22);
        pad2->SetTopMargin(0.03);
        pad2->SetBottomMargin(0.30);
        pad2->Draw();
        pad2->cd();

        TH1D *h_ratio = (TH1D*)h_data_mT->Clone("c_mT_ratio");
        h_ratio->Divide(hMCsum);
        h_ratio->SetTitle("");
        h_ratio->SetMarkerStyle(20);
        h_ratio->SetMarkerSize(0.8);
        h_ratio->SetLineColor(kBlack);

        h_ratio->GetYaxis()->SetTitle("Data / MC");
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.12);
        h_ratio->GetYaxis()->SetLabelSize(0.10);
        h_ratio->GetYaxis()->SetTitleOffset(0.45);

        h_ratio->GetXaxis()->SetTitle("m_{T} [GeV]");
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->GetXaxis()->SetLabelSize(0.10);

        h_ratio->SetMinimum(0.0);
        h_ratio->SetMaximum(2.0);
        h_ratio->Draw("E1");

        TLine *line = new TLine(0, 1.0, 200, 1.0);
        line->SetLineStyle(2);
        line->SetLineWidth(1);
        line->Draw("SAME");

        c->Update();
    }

    {
        TH1D *hMCsum = (TH1D*)h_W_MET->Clone("c_MET_MCsum");
        hMCsum->Add(h_DY_MET);
        hMCsum->Add(h_TT_MET);

        const double ymax = max(hMCsum->GetMaximum(), h_data_MET->GetMaximum());

        TCanvas *c = new TCanvas("c_MET", "c_MET", 700, 700);

        TPad *pad1 = new TPad("c_MET_pad1", "pad1", 0, 0.22, 1, 1);
        pad1->SetBottomMargin(0.02);
        pad1->Draw();
        pad1->cd();

        THStack *stack = new THStack("c_MET_stack", "");
        h_TT_MET->SetFillColor(kGreen);
        h_TT_MET->SetLineColor(kBlack);
        h_DY_MET->SetFillColor(kRed);
        h_DY_MET->SetLineColor(kBlack);
        h_W_MET ->SetFillColor(kOrange);
        h_W_MET ->SetLineColor(kBlack);

        stack->Add(h_TT_MET);
        stack->Add(h_DY_MET);
        stack->Add(h_W_MET);

        stack->Draw("HIST");
        stack->SetMinimum(0.0);
        stack->SetMaximum(1.25 * ymax);
        stack->GetXaxis()->SetLimits(0, 200);
        stack->GetXaxis()->SetLabelSize(0);
        stack->GetXaxis()->SetTitleSize(0);
        stack->GetYaxis()->SetTitle("Events");

        h_data_MET->SetMarkerStyle(20);
        h_data_MET->SetMarkerSize(1.0);
        h_data_MET->SetLineColor(kBlack);
        h_data_MET->Draw("E1 SAME");

        TLegend *leg = new TLegend(0.60, 0.64, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_data_MET, "Data", "lep");
        leg->AddEntry(h_W_MET,  "W#rightarrow#mu#nu (recoil)", "f");
        leg->AddEntry(h_DY_MET, "DY (Z/#gamma*#rightarrow ll)", "f");
        leg->AddEntry(h_TT_MET, "t#bar{t}", "f");
        leg->Draw();

        c->cd();
        TPad *pad2 = new TPad("c_MET_pad2", "pad2", 0, 0.00, 1, 0.22);
        pad2->SetTopMargin(0.03);
        pad2->SetBottomMargin(0.30);
        pad2->Draw();
        pad2->cd();

        TH1D *h_ratio = (TH1D*)h_data_MET->Clone("c_MET_ratio");
        h_ratio->Divide(hMCsum);
        h_ratio->SetTitle("");
        h_ratio->SetMarkerStyle(20);
        h_ratio->SetMarkerSize(0.8);
        h_ratio->SetLineColor(kBlack);

        h_ratio->GetYaxis()->SetTitle("Data / MC");
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.12);
        h_ratio->GetYaxis()->SetLabelSize(0.10);
        h_ratio->GetYaxis()->SetTitleOffset(0.45);

        h_ratio->GetXaxis()->SetTitle("MET [GeV]");
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->GetXaxis()->SetLabelSize(0.10);

        h_ratio->SetMinimum(0.0);
        h_ratio->SetMaximum(2.0);
        h_ratio->Draw("E1");

        TLine *line = new TLine(0, 1.0, 200, 1.0);
        line->SetLineStyle(2);
        line->SetLineWidth(1);
        line->Draw("SAME");

        c->Update();
    }

    std::cout << "sumW(raw)  = " << sW  << "\n";
    std::cout << "sumDY(raw) = " << sDY << "\n";
    std::cout << "sumTT(raw) = " << sTT << "\n";
    std::cout << "kW  = " << kW  << "\n";
    std::cout << "kDY = " << kDY << "\n";
    std::cout << "kTT = " << kTT << "\n";

    std::cout << "Data mT integral = " << h_data_mT->Integral() << "\n";
    std::cout << "MC   mT integral = " << (h_W_mT->Integral() + h_DY_mT->Integral() + h_TT_mT->Integral()) << "\n";

    std::cout << "Data MET integral = " << h_data_MET->Integral() << "\n";
    std::cout << "MC   MET integral = " << (h_W_MET->Integral() + h_DY_MET->Integral() + h_TT_MET->Integral()) << "\n";
}
