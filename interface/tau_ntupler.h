#include "ROOT/RVec.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <nlohmann/json.hpp>
 
using namespace ROOT;
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using str = const std::string &;

using cRVecF = const ROOT::RVecF &;
using cRVecI = const ROOT::RVecI &;
using cRVecC = const ROOT::RVecC &;
using cRVecU = const ROOT::RVecU &;
using cRVecB = const ROOT::RVecB &;

float deltaR(float eta_1, float eta_2, float phi_1, float phi_2){
   const float deta = eta_1 - eta_2;
   const float dphi = ROOT::Math::VectorUtil::Phi_mpi_pi(phi_1 - phi_2);
   const float dRsq = std::pow(deta,2) + std::pow(dphi,2);

   return sqrt(dRsq);
}

auto jsonFilterlambda(uint run, uint luminosity) {
  std::ifstream i("/eos/user/b/boguo/botao/CMSSW_10_6_29/src/PhysicsTools/NanoAODTools/TAU-Trigger-NANO/Collisions24_13p6TeV_378981_381384_DCSOnly_TkPx.txt");
  nlohmann::json golden_json;
  i >> golden_json;
  bool matched = false;
  // check if the run exists
  if (golden_json.find(std::to_string(run)) != golden_json.end()) {
    // std::cout << "run : " << run << std::endl;
    // std::cout << "luminosity : " << luminosity << std::endl;
    for (auto &luminosityrange : golden_json[std::to_string(run)]) {
      if (luminosity >= luminosityrange[0] &&
        luminosity <= luminosityrange[1]) {
        // std::cout << "luminosity : " << luminosity << std::endl;
        matched = true;
        break;
      }
    }
  }
  return matched;
};

auto jsonFilterlambda2023(uint run, uint luminosity) {
  std::ifstream i("/eos/user/b/boguo/botao/CMSSW_10_6_29/src/PhysicsTools/NanoAODTools/TAU-Trigger-NANO/Cert_Collisions2023_366442_370790_GoldenJSON.txt");
  nlohmann::json golden_json;
  i >> golden_json;
  bool matched = false;
  // check if the run exists
  if (golden_json.find(std::to_string(run)) != golden_json.end()) {
  // now loop over all luminosity blocks and check if the event is
  // valid
    for (auto &luminosityrange : golden_json[std::to_string(run)]) {
      if (luminosity >= luminosityrange[0] &&
        luminosity <= luminosityrange[1]) {
        matched = true;
        break;
      }
    }
  }
  return matched;
};

// Muon_pfRelIso04_all
// muon_pt > 24, muon_eta < 2.1, muon_iso < 0.1, medium_id
int SelectMuon(cRVecF muon_pt, cRVecF muon_eta, cRVecF muon_phi, cRVecF muon_iso, cRVecB muon_mediumId)
{
  int idx = -1;
  float min_muon_iso = 0.1;
  for(auto i = 0; i < muon_pt.size(); i++) {
    if (muon_pt[i] > 24 && abs(muon_eta[i]) < 2.1 && muon_iso[i] < 0.1 && muon_mediumId[i] == true) {
      // pick the index of the muon that has min isolation
      if (muon_iso[i] < min_muon_iso) {
        idx = i;
        min_muon_iso = muon_iso[i];
      }
    }
  }
  return idx;
}
// match the trig obj and muon with dR < 0.5
// cRVecU trig_id, cRVecF trig_eta, cRVecF trig_phi, match trigobj dR < 0.5
bool Muon_match(cRVecU trig_id, cRVecF trig_eta, cRVecF trig_phi, float sig_muon_eta, float sig_muon_phi)
{
  for(auto i = 0; i < trig_eta.size(); i++){
    float dR = deltaR(trig_eta[i], sig_muon_eta, trig_phi[i], sig_muon_phi);
    if (dR < 0.5) {
      if (trig_id[i] == 13) {
        return true;
      }
    }
  }
  return false;
}

int SelectTau(cRVecF tau_pt, cRVecF tau_eta, cRVecC tau_idDeepTau2018v2p5VSe, cRVecC tau_idDeepTau2018v2p5VSmu, cRVecC tau_idDeepTau2018v2p5VSjet)
{
  int idx = -1;
  float max_tau_pt = -1.0;
  for(auto i =0; i < tau_pt.size(); i++) {
    if (tau_pt[i] > 20 && abs(tau_eta[i]) < 2.3) {
      if (tau_idDeepTau2018v2p5VSe[i] >= 2 && tau_idDeepTau2018v2p5VSmu[i] >= 4 && tau_idDeepTau2018v2p5VSjet[i] >=5) {
        // pick the index of the tau that has max pt
        if (tau_pt[i] > max_tau_pt) {
          idx = i;
          max_tau_pt = tau_pt[i];
        }
      }
    }
  }
  return idx;
}

bool VetoEle(cRVecF ele_pt, cRVecF ele_eta, cRVecF ele_mvaIso)
{
  for(auto i = 0; i < ele_pt.size(); i++) {
    if (ele_pt[i] > 10 && abs(ele_eta[i]) < 2.5 && ele_mvaIso[i] > 0.5) {
      return false;
    }
  }
  return true;
}

bool Tau_match(cRVecU trig_id, cRVecF trig_eta, cRVecF trig_phi, float leading_tau_eta, float leading_tau_phi)
{
  for(auto i = 0; i < trig_eta.size(); i++){
    float dR = deltaR(trig_eta[i], leading_tau_eta, trig_phi[i], leading_tau_phi);
    if (dR < 0.5) {
      if (trig_id[i] == 15) {
        return true;
      }
    }
  }
  return false;
}

// return math.sqrt( 2.0 * lepton_v4.Pt() * met.Pt() * (1.0 - math.cos(delta_phi)) )
float MT(float pt_1, float pt_2, float phi_1, float phi_2) {
  const float dphi = ROOT::Math::VectorUtil::Phi_mpi_pi(phi_1 - phi_2);
  const float mt = sqrt(2.0 * pt_1 * pt_2 * (1.0 - std::cos(dphi)));
  return mt;
}

float VisMass(float muon_pt, float muon_eta, float muon_phi, float muon_mass, float tau_pt, float tau_eta, float tau_phi, float tau_mass) {
  const ROOT::Math::PtEtaPhiMVector muon_v4(muon_pt, muon_eta, muon_phi, muon_mass);
  const ROOT::Math::PtEtaPhiMVector tau_v4(tau_pt, tau_eta, tau_phi, tau_mass);
  float vis_mass = (muon_v4 + tau_v4).M();
  return vis_mass;
}

int FindLeadingTau(cRVecF tau_pt)
{
    int idx = -1; // Return -1 if no good lepton is found.
    double leading_pt = 0.0;
    for(auto i = 0; i < tau_pt.size(); i++) {
        if (tau_pt[i] > leading_pt) {
            leading_pt = tau_pt[i];
            idx = i;
        }
    }
    return idx;
}

int FindSubLeadingTau(cRVecF tau_pt, int leading_idx)
{
    int idx = -1; // Return -1 if no sub leading lepton is found.
    double subleading_pt = 0.0;
    for(auto i = 0; i < tau_pt.size(); i++) {
        if (i == leading_idx)
          continue;
        if (tau_pt[i] > subleading_pt) {
            subleading_pt = tau_pt[i];
            idx = i;
        }
    }
    return idx;
}

bool SelectVBFJets(cRVecF jet_pt, cRVecF jet_eta, cRVecF jet_phi, cRVecF jet_mass, float tau_eta,float tau_phi)
{
  int jet_num = 0;
  for(auto i =0; i < jet_pt.size(); i++) {
    float dR = deltaR(jet_eta[i],tau_eta,jet_phi[i],tau_phi);
    if (dR > 0.5) {
      if (jet_pt[i] > 45 && abs(jet_eta[i]) < 4.7) {
        jet_num +=1;
      }
    }
  }
  if (jet_num >= 2) return true;
  else return false;
}
bool SelectJet(cRVecF jet_pt, cRVecF jet_eta, cRVecF jet_phi, cRVecF jet_mass, float tau_eta,float tau_phi)
{
  int jet_num = 0;
  for(auto i =0; i < jet_pt.size(); i++) {
    float dR = deltaR(jet_eta[i],tau_eta,jet_phi[i],tau_phi);
    if (dR > 0.5) {
      if (jet_pt[i] > 60 && abs(jet_eta[i]) < 4.7) {
        jet_num +=1;
      }
    }
  }
  if (jet_num >= 1) return true;
  else return false;
}

// 0 => Loose, 1 => Medium, 2 => Tight, 3 => DeepTau no spec WP, 4 => PNet no specified WP, 
// 5 => ChargedIso, 6 => Dxy, 7 => e-tau inside filter, 8 => mu-tau inside filter, 
// 9 => Single Tau, 10 => VBF DiTau, 11 => di-tau, 12 => e-tau, 13 => mu-tau, 14 => di-tau + PFJet, 
// 15 => e-tau displaced, 16 => mu-tau displaced, 17 => di-tau displaced, 18 => Monitoring, 
// 19 => VBF SingleTau Monitoring, 20 => DiTau+Jet Monitoring, 21 => Monitoring muTau displaced, 
// 22 => OneProng, 23 => DiTau Monitoring, 24 => OverlapFilter, 25 => VBF DiTau monitoring, 
// 26 => SingleTau Monitoring, 27 => MatchL1HLT, 28 => HPS, 29 => single PF-tau inside filter, 30 => VBF SingleTau for Tau

// ditau pnet
// HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3
// HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1 (+Tight)
// PNet monitoring bit 1, 4, 23
bool PassDiTauPNet(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi, int wp){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 4 => PNet no specified WP, 23 => DiTau Monitoring
      if((trig_bits[i] & (1<<wp)) != 0 && (trig_bits[i] & (1<<4)) != 0 && (trig_bits[i] & (1<<23)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 30 && trig_l1pt[i] > 34 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}

// ditau deeptau
// HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1
// HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1
// DeepTau monitoring bit 3, 23
bool PassDiTauDeepTau(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 23 => DiTau Monitoring
      if((trig_bits[i] & (1<<3)) != 0 && (trig_bits[i] & (1<<23)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 35 && trig_l1pt[i] > 34 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}

// mutau pnet
// HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1 (+Medium, Tight)
// HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1
// HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1
// PNet monitoring bit 0, 4, 13
// PNet Medium 1, 4, 13
// PNet Tight 2, 4, 13
// take int wp as input, should in [0,1,2]
bool PassMuTauPNet(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi, int wp){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 0 => Loose, 4 => PNet no specified WP, 13 => mu-tau
      if((trig_bits[i] & (1<<wp)) != 0 && (trig_bits[i] & (1<<4)) != 0 && (trig_bits[i] & (1<<13)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 27 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// mutau deeptau
// HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1
// DeepTau monitoring bit 1, 3, 13
bool PassMuTauDeepTau(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 13 => mu-tau
      if((trig_bits[i] & (1<<3)) != 0 && (trig_bits[i] & (1<<13)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 27 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// ETau PNet
bool PassETauPNet(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi, int wp){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 0 => Loose, 4 => PNet no specified WP, 27 => MatchL1HLT
      if((trig_bits[i] & (1<<wp)) != 0 && (trig_bits[i] & (1<<4)) != 0 && (trig_bits[i] & (1<<27)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 30 && trig_l1iso[i] > 0 && trig_l1pt[i] > 26 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// ETau PNet check
bool PassETauPNet_nobit4(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi, int wp){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 0 => Loose, 4 => PNet no specified WP, 27 => MatchL1HLT
      if((trig_bits[i] & (1<<wp)) != 0 && (trig_bits[i] & (1<<27)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 30 && trig_l1iso[i] > 0 && trig_l1pt[i] > 26 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// ETau PNet check
bool PassETauPNet_nobit27(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi, int wp){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 0 => Loose, 4 => PNet no specified WP, 27 => MatchL1HLT
      if((trig_bits[i] & (1<<wp)) != 0 && (trig_bits[i] & (1<<4)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 30 && trig_l1iso[i] > 0 && trig_l1pt[i] > 26 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// ETau PNet check (using mutau bit and extra etau selection)
bool PassETauPNet_withMutaubit(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi, int wp){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 0 => Loose, 4 => PNet no specified WP, 13 => mu-tau
      if((trig_bits[i] & (1<<wp)) != 0 && (trig_bits[i] & (1<<4)) != 0 && (trig_bits[i] & (1<<13)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 30 && trig_l1iso[i] > 0 && trig_l1pt[i] > 26 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}

// ETau DeepTau
bool PassETauDeepTau(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 27 => MatchL1HLT
      if((trig_bits[i] & (1<<3)) != 0 && (trig_bits[i] & (1<<27)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 30 && trig_l1iso[i] > 0 && trig_l1pt[i] > 26 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}
bool PassETauDeepTau_withMutaubit(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 27 => MatchL1HLT
      if((trig_bits[i] & (1<<3)) != 0 && (trig_bits[i] & (1<<13)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 30 && trig_l1iso[i] > 0 && trig_l1pt[i] > 26 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// DiTaujet PNet
bool PassDiTaujetPNet(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 4 => PNet no specified WP, 20 => DiTau+Jet Monitoring
      if((trig_bits[i] & (1<<4)) != 0 && (trig_bits[i] & (1<<20)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 26 && trig_l1iso[i] > 0 && trig_l1pt[i] > 26 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// DiTaujet DeepTau
bool PassDiTaujetDeepTau(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 20 => DiTau+Jet Monitoring
      if((trig_bits[i] & (1<<3)) != 0 && (trig_bits[i] & (1<<20)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 30 && trig_l1iso[i] > 0 && trig_l1pt[i] > 26 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}

// singletau PNet
bool PassSingleTauPNet(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi, int wp){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 0 => Loose, 4 => PNet no specified WP, 26 => SingleTau Monitoring
      if((trig_bits[i] & (1<<wp)) != 0 && (trig_bits[i] & (1<<4)) != 0 && (trig_bits[i] & (1<<26)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 130 && trig_l1pt[i] > 130 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}
/// singletau PNet no filter bit
bool PassSingleTauPNet_nofilterbit(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi, int wp){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 0 => Loose, 4 => PNet no specified WP, 26 => SingleTau Monitoring
      if(1){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 130 && trig_l1pt[i] > 130 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// singletau DeepTau
bool PassSingleTauDeepTau(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 26 => SingleTau Monitoring
      if((trig_bits[i] & (1<<3)) != 0 && (trig_bits[i] & (1<<26)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 180 && trig_l1pt[i] > 130 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// singletau DeepTau no filter bit
bool PassSingleTauDeepTau_nofilterbit(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 26 => SingleTau Monitoring
      if(1){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 180 && trig_l1pt[i] > 130 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}

// VBF singletau PNet
bool PassVBFSingleTauPNet(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 4 => PNet no specified WP, 30 => VBF SingleTau for Tau
      if((trig_bits[i] & (1<<4)) != 0 && (trig_bits[i] & (1<<30)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 45 && trig_l1iso[i] > 0 && trig_l1pt[i] > 45 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// VBF singletau DeepTau
bool PassVBFSingleTauDeepTau(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 30 => VBF SingleTau for Tau
      if((trig_bits[i] & (1<<3)) != 0 && (trig_bits[i] & (1<<30)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 45 && trig_l1iso[i] > 0 && trig_l1pt[i] > 45 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}

// vbfditau PNet
bool PassVBFDiTauPNet(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 4 => PNet no specified WP, 25 => VBF DiTau monitoring
      if((trig_bits[i] & (1<<4)) != 0 && (trig_bits[i] & (1<<25)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 20 && abs(trig_eta[i]) < 2.2 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// vbfditau DeepTau
bool PassVBFDiTauDeepTau(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 25 => VBF DiTau monitoring
      if((trig_bits[i] & (1<<3)) != 0 && (trig_bits[i] & (1<<25)) != 0){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 20 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}


/////////
/////////
/////////
// VBF singletau PNet check
bool PassVBFSingleTauPNet_nofilterbit(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 4 => PNet no specified WP, 30 => VBF SingleTau for Tau
      if(1){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 45 && trig_l1iso[i] > 0 && trig_l1pt[i] > 45 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// VBF singletau DeepTau
bool PassVBFSingleTauDeepTau_nofilterbit(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecI trig_l1iso, cRVecF trig_l1pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 30 => VBF SingleTau for Tau
      if(1){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 45 && trig_l1iso[i] > 0 && trig_l1pt[i] > 45 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}

bool PassMuTauPNet_nofilterbit(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi, int wp){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 0 => Loose, 4 => PNet no specified WP, 13 => mu-tau
      if(1){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 27 && abs(trig_eta[i]) < 2.3 ) {
          return true;
        }
      }
    }
  }
  return false;
}
// mutau deeptau
// HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1
// DeepTau monitoring bit 1, 3, 13
bool PassMuTauDeepTau_nofilterbit(cRVecU trig_id, cRVecI trig_bits, cRVecF trig_pt, cRVecF trig_eta, cRVecF trig_phi, float tau_pt, float tau_eta, float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(auto i=0; i < trig_pt.size(); i++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[i],trig_eta[i],trig_phi[i],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 3 => DeepTau no specified WP, 13 => mu-tau
      if(1){ 
        if ( trig_id[i] == 15 && trig_pt[i] > 27 && abs(trig_eta[i]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}


////////
// HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1
bool PassMuTauTrig2023DeepTau(uint ntrig,cRVecU trig_id,cRVecI trig_bits,cRVecF trig_pt,cRVecF trig_eta,cRVecF trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & (1<<7)) != 0 && (trig_bits[it] & (1<<26)) != 0 && trig_id[it] == 15){ 
        if ( trig_id[it] == 15 && trig_pt[it] > 27 && abs(trig_eta[it]) < 2.1 ) {
          return true;
        }
      }
    }
  }
  return false;
}
