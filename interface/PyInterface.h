/*! Definition of c++ methods used in python code.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#pragma once

enum class LegType { e = 1, mu = 2, tau = 4, jet = 8 };

class TriggerMatchProvider {
public:
    struct MatchDescriptor {
        ULong64_t match_mask{0};
        std::vector<UInt_t> filter_hashes;
        int min_run{-1}, max_run{-1};
        float hltObj_pt{-1}, l1Tau_pt{-1};
        int l1Tau_hwIso{-1};

        MatchDescriptor() {}
        MatchDescriptor(ULong64_t _match_mask, const std::vector<UInt_t>& _filter_hashes, int _min_run, int _max_run,
                        float _hltObj_pt, float _l1Tau_pt, float _l1Tau_hwIso) :
                match_mask(_match_mask), filter_hashes(_filter_hashes), min_run(_min_run), max_run(_max_run),
                hltObj_pt(_hltObj_pt), l1Tau_pt(_l1Tau_pt), l1Tau_hwIso(_l1Tau_hwIso)
        {
        }
    };

    void Add(int channel_id, const MatchDescriptor& desc)
    {
        channel_matches[channel_id].push_back(desc);
    }

    bool Pass(int channel_id, UInt_t run, float tau_eta, float tau_phi, ULong64_t hlt_accept, float deltaRThr,
              const ROOT::VecOps::RVec<UInt_t>& hltObj_types, const ROOT::VecOps::RVec<float>& hltObj_pt,
              const ROOT::VecOps::RVec<float>& hltObj_eta, const ROOT::VecOps::RVec<float>& hltObj_phi,
              const ROOT::VecOps::RVec<ULong64_t>& hltObj_hasPathName, const ROOT::VecOps::RVec<UInt_t>& filter_hltObj,
              const ROOT::VecOps::RVec<UInt_t>& filter_hash) const
    {
        const auto desc_iter = channel_matches.find(channel_id);
        if(desc_iter != channel_matches.end()) {
            const float deltaRThr2 = std::pow(deltaRThr, 2);
            for(const MatchDescriptor& match_desc : desc_iter->second) {
                //if((hlt_accept & match_desc.match_mask) == 0) continue;
                if(match_desc.min_run >= 0 && run < match_desc.min_run) continue;
                if(match_desc.max_run >= 0 && run >= match_desc.max_run) continue;
                //if(match_desc.l1Tau_pt >= 0 && l1Tau_pt <= match_desc.l1Tau_pt) continue;
                //if(match_desc.l1Tau_hwIso >= 0 && l1Tau_hwIso <= match_desc.l1Tau_hwIso) continue;

                for(size_t n = 0; n < hltObj_pt.size(); ++n) {
                    //if((hltObj_types.at(n) & static_cast<UInt_t>(LegType::tau)) == 0) continue;
                    //if((hltObj_hasPathName.at(n) & match_desc.match_mask) == 0) continue;
                    if(match_desc.hltObj_pt >= 0 && hltObj_pt.at(n) <= match_desc.hltObj_pt) continue;
                    const float deta = tau_eta - hltObj_eta.at(n);
                    const float dphi = ROOT::Math::VectorUtil::Phi_mpi_pi(tau_phi - hltObj_phi.at(n));
                    const float deltaR2 = std::pow(deta, 2) + std::pow(dphi, 2);
                    if(deltaR2 >= deltaRThr2) continue;
                    //if(PassFilters(match_desc.filter_hashes, n, filter_hltObj, filter_hash)) return true;
                    // test for botao
                    return true;
                    // end test
                }
            }
        }
        return false;
    }

    static TriggerMatchProvider& Initialize()
    {
        default_provider.reset(new TriggerMatchProvider());
        return GetDefault();
    }

    static TriggerMatchProvider& GetDefault()
    {
        if(!default_provider)
            throw std::runtime_error("Default TriggerMatchProvider is not initialized.");
        return *default_provider;
    }

private:
    static std::unique_ptr<TriggerMatchProvider> default_provider;

    bool PassFilters(const std::vector<UInt_t>& filter_hashes, size_t hltObj_index,
                     const ROOT::VecOps::RVec<UInt_t>& filter_hltObj,
                     const ROOT::VecOps::RVec<UInt_t>& filter_hash) const
    {
        for(UInt_t filter_ref : filter_hashes) {
            bool filter_found = false;
            for(size_t n = 0; n < filter_hltObj.size() && !filter_found; ++n) {
                filter_found = filter_hltObj.at(n) == hltObj_index && filter_hash.at(n) == filter_ref;
            }
            if(!filter_found) return false;
        }
        return true;
    }

private:
    std::map<int, std::vector<MatchDescriptor>> channel_matches;
};

class PileUpWeightProvider {
public:
    PileUpWeightProvider(const TH1D& data_pu_orig, const TH1D& mc_pu_orig)
    {
        TH1D data_pu(data_pu_orig);
        data_pu.Scale(1. / data_pu.Integral());
        TH1D mc_pu(mc_pu_orig);
        mc_pu.Scale(1. / mc_pu.Integral());
        ratio.reset(new TH1D(data_pu));
        ratio->Divide(&mc_pu);
    }

    float GetWeight(int npu) const
    {
        int bin = ratio->FindBin(npu);
        if(bin < 1 || bin > ratio->GetNbinsX())
            return 0;
        return ratio->GetBinContent(bin);
    }

    static void Initialize(const TH1D& data_pu, const TH1D& mc_pu)
    {
        default_provider.reset(new PileUpWeightProvider(data_pu, mc_pu));
    }

    static const PileUpWeightProvider& GetDefault()
    {
        if(!default_provider)
            throw std::runtime_error("Default PileUpWeightProvider is not initialized.");
        return *default_provider;
    }

private:
    static std::unique_ptr<PileUpWeightProvider> default_provider;

private:
    std::unique_ptr<TH1D> ratio;
};

std::unique_ptr<TriggerMatchProvider> TriggerMatchProvider::default_provider;
std::unique_ptr<PileUpWeightProvider> PileUpWeightProvider::default_provider;
