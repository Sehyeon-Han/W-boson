void N_obs(const char* data = "/data2/kplee/Lecture/CMSOpenData/Data/SingleMuon/Run2016H/*.root",
           const char* DYMC = "/data2/kplee/Lecture/CMSOpenData/MC2016/DY_M50_aMCNLO/*.root",
           const char* TTMC = "/data2/kplee/Lecture/CMSOpenData/MC2016/TTTo2L2Nu_Powheg/*.root",
           double lumi_pb = 8740.1,
           double xsec_DY = 6019.95,
           double xsec_TT = 88.51)
{
    TH1::SetDefaultSumw2();

    TChain chain_data("Events");
    chain_data.Add(data);
    TChain chain_DYMC("Events");
    chain_DYMC.Add(DYMC);
    TChain chain_TTMC("Events");
    chain_TTMC.Add(TTMC);
    
    auto yield = [](TChain &chain)
    {
        struct Result
        {
            double sumWeight_tot;
            double sumWeight_sel;
        } out = {0.0, 0.0};
        
        TTreeReader reader(&chain);

        TTreeReaderArray<float> Muon_pt(reader, "Muon_pt");
        TTreeReaderArray<float> Muon_eta(reader,"Muon_eta");
        TTreeReaderArray<float> Muon_iso(reader,"Muon_pfRelIso04_all");
        TTreeReaderArray<bool> Muon_tightId(reader, "Muon_tightId");
        TTreeReaderValue<bool> HLT_IsoMu24(reader, "HLT_IsoMu24");

        const bool hasGenWeight = (chain.GetListOfBranches()->FindObject("genWeight")) != nullptr;
        TTreeReaderValue<float> *genWeight = nullptr;
        if(hasGenWeight)
            genWeight = new TTreeReaderValue<float> (reader, "genWeight");

        while(reader.Next())
        {   
            double w = 1.0;
            if(genWeight)
            w = **genWeight;
            out.sumWeight_tot += w;

            int selCount = 0;
            int selIdx = -1;

            const int n = Muon_pt.GetSize();
            for(int i = 0; i < n; ++i)
            {
                if(Muon_pt[i] > 25.0 && fabs(Muon_eta[i]) < 2.1 && Muon_iso[i] < 0.12 && Muon_tightId[i])
                {
                    ++selCount;
                    selIdx = i;

                    if(selCount > 1)
                    break;
                }
            }

            if(selCount == 1 && selIdx >= 0 && *HLT_IsoMu24)
            {
                out.sumWeight_sel += w;
            }
        }
        return out;
    };

    auto y_data = yield(chain_data);
    auto y_DYMC = yield(chain_DYMC);
    auto y_TTMC = yield(chain_TTMC);

    double kDY = lumi_pb * xsec_DY / y_DYMC.sumWeight_tot;
    double N_DYMC = kDY * y_DYMC.sumWeight_sel;

    double kTT = lumi_pb * xsec_TT / y_TTMC.sumWeight_tot;
    double N_TTMC = kTT * y_TTMC.sumWeight_sel;

    double N_data = y_data.sumWeight_sel;
    double N_obs = N_data - (N_DYMC + N_TTMC);
    
    cout << fixed << setprecision(6);
    cout << "Number of data = " << N_data << endl;
    cout << "Number of DYMC = " << N_DYMC << endl;
    cout << "Number of TTMC = " << N_TTMC << endl;
    cout << "Number of Wobs = " << N_obs << endl;
}
