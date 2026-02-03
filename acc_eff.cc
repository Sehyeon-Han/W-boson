void acc_eff(const char* WMC  = "/data2/kplee/Lecture/CMSOpenData/MC2016/WJetsToLNu_aMCNLO/*.root")
{
    TH1::SetDefaultSumw2();

    TChain chain("Events");
    chain.Add(WMC);
    
    TTreeReader reader(&chain);

    TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
    TTreeReaderArray<float> Muon_eta(reader, "Muon_eta");
    TTreeReaderArray<float> Muon_iso(reader, "Muon_pfRelIso04_all");
    TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
    TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");

    TTreeReaderValue<float> genWeight(reader,"genWeight");
    TTreeReaderArray<float> GenPart_pt(reader, "GenPart_pt");
    TTreeReaderArray<float> GenPart_eta(reader, "GenPart_eta");
    TTreeReaderArray<int> GenPart_pdgId(reader, "GenPart_pdgId");
    TTreeReaderArray<int> GenPart_status(reader, "GenPart_status");
    TTreeReaderArray<int> GenPart_statusFlags(reader, "GenPart_statusFlags");
    
    double sumWeight_tot = 0;
    double sumWeight_acc = 0;
    double sumWeight_eff = 0;
    
    while(reader.Next())
    {
        bool pass_total = false;
        bool pass_acc = false;
        bool pass_eff = false;

        const int nacc = GenPart_pt.GetSize();
        int selIdx = -1;

        for(int i_lep = 0; i_lep < nacc; ++i_lep)
        {
            if(abs(GenPart_pdgId[i_lep]) != 13) continue;
            if(GenPart_status[i_lep] != 1) continue;
            if(!(abs(GenPart_statusFlags[i_lep]) & 256)) continue;

            pass_total = true;
            selIdx = i_lep;

            if(GenPart_pt[i_lep] > 25.0 && fabs(GenPart_eta[i_lep]) < 2.4)
                pass_acc = true;
                break;
        }

        if(pass_acc)
        {
            if(*HLT_IsoMu24)
            {
                const int neff = Muon_pt.GetSize();
                for(int i = 0; i < neff; ++i)
                {
                    if(Muon_pt[i] > 25.0 && fabs(Muon_eta[i]) < 2.4 && Muon_iso[i] < 0.12 && Muon_tightId[i])
                    {
                        pass_eff = true;
                        break;
                    }
                }
            }
        }

        if(pass_total)
            sumWeight_tot += *genWeight;
        if(pass_acc)
            sumWeight_acc += *genWeight;
        if(pass_eff)
            sumWeight_eff += *genWeight;
    }

    double acceptance = (sumWeight_acc / sumWeight_tot);
    double efficiency = (sumWeight_eff / sumWeight_acc);
        
    cout << "Acceptance = " << acceptance << endl;
    cout << "Efficiency = " << efficiency << endl;
}
