#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
 
using namespace ROOT;
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;

using Vec_t = const ROOT::RVec<float>&;
using Vec_i = const ROOT::RVec<int>&;


float getFloatValue(Vec_t vec, Int_t index) {
  if (index >= 0) return vec[index];
  else return -999.;
}

int getIntValue(Vec_i vec, Int_t index) {
  if (index >= 0) return vec[index];
  else return -999;
}

float ZMass(TLorentzVector tau_p4,TLorentzVector muon_p4) {
  if (tau_p4.Pt() <= 0 || muon_p4.Pt() <= 0)
    return -999.;
  return (tau_p4 + muon_p4).M();
}

float CalcMT(TLorentzVector lep_p4,float met_pt,float met_phi){
  if (lep_p4.Pt() <= 0)
    return -999.;
  TLorentzVector met_p4;
  met_p4.SetPtEtaPhiE(met_pt,0,met_phi,0);
  const double delta_phi = ROOT::Math::VectorUtil::DeltaPhi(lep_p4, met_p4);
  return std::sqrt( 2.0 * lep_p4.Pt() * met_p4.Pt() * ( 1.0 - std::cos(delta_phi) ) );

}

float deltaR(float eta_1, float eta_2, float phi_1, float phi_2){
   const float deta = eta_1 - eta_2;
   const float dphi = ROOT::Math::VectorUtil::Phi_mpi_pi(phi_1 - phi_2);
   const float dRsq = std::pow(deta,2) + std::pow(dphi,2);

   return sqrt(dRsq);
}

int PassTagFilter(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float muon_pt,float muon_eta,float muon_phi){
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),muon_eta,trig.Phi(),muon_phi);
    if (dR < 0.5){
      if((trig_bits[it] & 8) != 0 && trig_id[it] == 13){
        return it;
      }
    }
  }
  return -1;
}

int MuonIndex(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,UInt_t nMu,Vec_t pt_1, Vec_t eta_1, Vec_t phi_1, Vec_t mass_1,Vec_t pfIso){
  int mu_index = -1;
  if(nMu > 0){
    for(int imu = nMu - 1; imu >= 0; imu--){
      const ROOT::Math::PtEtaPhiMVector muon(pt_1[imu], eta_1[imu], phi_1[imu], mass_1[imu]);
      float mu_iso = pfIso[imu]/muon.Pt();
      if(muon.Pt() < 24)continue;
      if(std::fabs(muon.Eta()) > 2.1)continue;
      //if(mu_iso > 0.1)continue;
      if(PassTagFilter(ntrig,trig_id,trig_bits,trig_pt,trig_eta,trig_phi,muon.Pt(),muon.Eta(),muon.Phi()) < 0) continue;
      mu_index = imu;
    }
  }
  return mu_index;
}

int MuonIndexFull(UInt_t nTau, UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,UInt_t nMu,Vec_t pt_1, Vec_t eta_1, Vec_t phi_1, Vec_t mass_1, Vec_t pfIso) {
  if (nMu == 1 && nTau >= 1) {
    return MuonIndex(ntrig, trig_id, trig_bits, trig_pt, trig_eta, trig_phi, nMu, pt_1, eta_1, phi_1, mass_1, pfIso);
  }
  return -999;
}

int TauIndex(UInt_t ntau, Vec_t pt_1, Vec_t eta_1, Vec_t phi_1, Vec_t mass_1, Vec_t dz_1, TLorentzVector muon_p4,Vec_t Tau_rawIsodR03){
  int tau_index = -1;
  float tau_iso = -9999;
  if(ntau > 0 && muon_p4.Pt() > 0.){
    // for (int itau = 0; itau <= ntau; itau++){
      for (int itau = ntau - 1; itau >= 0; itau--){
      const ROOT::Math::PtEtaPhiMVector tau(pt_1[itau], eta_1[itau], phi_1[itau], mass_1[itau]);
      if(dz_1[itau] > 0.2)continue;
      if(deltaR(tau.Eta(), muon_p4.Eta(), tau.Phi(), muon_p4.Phi()) < 0.5) continue;
      if(tau.Pt() < 18) continue;
      if(std::fabs(tau.Eta()) > 2.1) continue; 
      if(Tau_rawIsodR03[itau] > tau_iso){
        tau_iso = Tau_rawIsodR03[itau];
        tau_index = itau;
      }
    }  
  }
  return tau_index;
}

// int JetIndex(UInt_t njet, Vec_t pt_1, Vec_t eta_1, Vec_t phi_1, Vec_t mass_1, Vec_t pu_id, Vec_t jet_id, TLorentzVector muon_p4, TLorentzVector tau_p4){
int JetIndex(UInt_t njet, Vec_t pt_1, Vec_t eta_1, Vec_t phi_1, Vec_t mass_1, Vec_t jet_id, TLorentzVector muon_p4, TLorentzVector tau_p4) {
  int jet_index = -1;
  if(njet > 0 && muon_p4.Pt() > 0 && tau_p4.Pt() > 0) {
    for (int ijet = njet - 1; ijet >= 0; ijet--){
      const ROOT::Math::PtEtaPhiMVector jet(pt_1[ijet], eta_1[ijet], phi_1[ijet], mass_1[ijet]);
      // if ((pu_id[ijet] < 4 && pt_1[ijet] <= 50) || jet_id[ijet] < 2) continue;
      if (jet_id[ijet] < 2) continue;
      if (deltaR(jet.Eta(), muon_p4.Eta(), jet.Phi(), muon_p4.Phi()) < 0.5) continue;
      if (deltaR(jet.Eta(), tau_p4.Eta(), jet.Phi(), tau_p4.Phi()) < 0.5) continue;
      if (jet.Pt() < 18) continue;
      jet_index = ijet;
    }  
  }
  return jet_index;
}

int PassDiTauFilter(UInt_t ntrig, Vec_i trig_id, Vec_i trig_bits, Vec_t trig_pt, Vec_t trig_eta, Vec_t trig_phi, float tau_pt, float tau_eta, float tau_phi) {
  if (tau_pt <= 0)
    return -1;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (1){ //dR < 0.5
      if((trig_bits[it] & 512) != 0 && (trig_bits[it] & 1024) != 0){ 
        return it;
      }
    }
  }
  return -1;
}

int PassDiTauJetFilter(UInt_t ntrig, Vec_i trig_id, Vec_i trig_bits, Vec_t trig_pt, Vec_t trig_eta, Vec_t trig_phi, float jet_pt, float jet_eta, float jet_phi){
  if (jet_pt <= 0)
    return -1;
  for(int it = 0; it < ntrig; it++) {
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it], trig_eta[it], trig_phi[it], 0);
    float dR = deltaR(trig.Eta(), jet_eta, trig.Phi(), jet_phi);
    if (dR < 0.5) {
      if((trig_bits[it] & 2097152) != 0 && trig_id[it] == 1) {
        if(trig_id[it] == 1) {
          return it;
        }
      }
    }
  }
  return -1;
}

int PassMuTauFilter(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return -1;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (1){ //dR < 0.5
      if((trig_bits[it] & 512) != 0){ 
          return it;
      }
    }
  }
  return -1;
}

int PassElTauFilter(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,Vec_t trig_l1pt,Vec_i trig_l1iso,float tau_pt,float tau_eta,float tau_phi){
   
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (1){ //dR < 0.5
      if((trig_bits[it] & 512) != 0){ 
        if(trig_l1pt[it] > 26 && trig_l1iso[it] > 0)
          return it;
      }
    }
  }
  return -1;
}

float LeadingTauPT(Vec_t tau_pt,Vec_t tau_eta,Vec_t tau_phi,Vec_t tau_m,int index){
  float pt = 0;
  if(index >= 0){
    const ROOT::Math::PtEtaPhiMVector tau(tau_pt[index],tau_eta[index],tau_phi[index],tau_m[index]);
    pt = tau.Pt();
  }
  return pt;
}

float LeadingTauEta(Vec_t tau_pt,Vec_t tau_eta,Vec_t tau_phi,Vec_t tau_m,int index){
  float eta = -999;
  if(index >= 0){
    const ROOT::Math::PtEtaPhiMVector tau(tau_pt[index],tau_eta[index],tau_phi[index],tau_m[index]);
    eta = tau.Eta();
  }
  return eta;
}
TLorentzVector Obj_p4(int index,Vec_t pt_1, Vec_t eta_1, Vec_t phi_1, Vec_t mass_1){
  TLorentzVector vec_p4(0, 0, 0, 0);
  if(index >= 0){
    const ROOT::Math::PtEtaPhiMVector p4(pt_1[index], eta_1[index], phi_1[index], mass_1[index]);
    vec_p4.SetPtEtaPhiM(p4.Pt(), p4.Eta(), p4.Phi(), p4.M());
  }
  return vec_p4;
}

bool PassBtagVeto(TLorentzVector muon_p4, TLorentzVector tau_p4,UInt_t njet, Vec_t jet_pt, Vec_t jet_eta, Vec_t jet_phi, Vec_t jet_m, Vec_t Jet_btagCSVV2){
  if(njet > 0){
    for(int j = 0; j <= njet; j++){
      const ROOT::Math::PtEtaPhiMVector jet(jet_pt[j], jet_eta[j], jet_phi[j], jet_m[j]);
      if(deltaR(muon_p4.Eta(),jet.Eta(),muon_p4.Phi(),jet.Phi()) > 0.5 && deltaR(tau_p4.Eta(),jet.Eta(),tau_p4.Phi(),jet.Phi()) > 0.5 && jet.Pt() > 20 && std::fabs(jet.Eta()) < 2.4 && Jet_btagCSVV2[j] > 0.0494){
	return true;
      }
    }
  }
  return false;
}
float TauL1_PT(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_l1pt,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  float taul1_pt = 0;
   for(int it=0; it < ntrig; it++){
     const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
     float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
     if (dR < 0.5){
       if( trig_id[it] == 15){ 
	 taul1_pt = trig_l1pt[it];
      }
   }
  }
  return taul1_pt;
}


// add by botao
bool PassMuTauTrig(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & 512) != 0 && trig_id[it] == 15 && trig_pt[it] > 27){ 
          return true;
      }
    }
  }
  return false;
}
// (trig_bits[it] & 512) != 0 && trig_id[it] == 15 && trig_pt[it] > 27

bool PassElTauTrig(UInt_t ntrig,Vec_t trig_l1pt,Vec_i trig_l1iso,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & 512) != 0 && trig_id[it] == 15){ 
          if(trig_l1pt[it] > 26 && trig_l1iso[it] > 0 && trig_pt[it] > 30)
            return true;
      }
    }
  }
  return false;
}
// TrigObj_id==15 && (TrigObj_filterBits&256)!=0

bool PassDiTauTrig(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & 512) != 0 && (trig_bits[it] & 1024) != 0){ 
          return true;
      }
    }
  }
  return false;
}
// (trig_bits[it] & 64) != 0 && trig_id[it] == 15 && ( ((trig_bits[it] & 4) != 0 && (trig_bits[it] & 16) != 0) || (trig_pt[it] > 40 && ( ((trig_bits[it] & 2) != 0 && (trig_bits[it] & 16) != 0) || (trig_bits[it] & 4) != 0 ) ) )

bool PassDiTauTrigMC(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & 512) != 0 && (trig_bits[it] & 1024) != 0){ 
          return true;
      }
    }
  }
  return false;
}
// (trig_bits[it] & 64) != 0 && trig_id[it] == 15

/// 2016 
bool PassMuTauTrig2016(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & 32) != 0 && trig_id[it] == 15 && (trig_bits[it] & 1) != 0){ 
          return true;
      }
    }
  }
  return false;
}
bool PassElTauTrig2016(UInt_t ntrig,Vec_t trig_l1pt,Vec_i trig_l1iso,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi,int run){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & 32) != 0 && trig_id[it] == 15 && (trig_bits[it] & 1) != 0){ 
          // if(trig_l1pt[it] > 26 && trig_l1iso[it] > 0 && trig_pt[it] > 30)
          if (run<276215) {
            return true;
          } else if (run >= 276215 && run < 278270) {
            return true;
          } else if (run >= 278270) {
            if(trig_l1pt[it] > 26 && trig_id[it] == 15 && trig_pt[it] > 30) {
              return true;
            }
          }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2016(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & 2) != 0 && trig_id[it] == 15){ 
          return true;
      }
    }
  }
  return false;
}

// nano v12
// 0 => LooseChargedIso, 1 => MediumChargedIso, 2 => TightChargedIso, 3 => DeepTau, 4 => TightID OOSC photons, 5 => HPS, 
// 6 => charged iso di-tau, 7 => deeptau di-tau, 8 => e-tau, 9 => mu-tau, 10 => single-tau/tau+MET, 11 => run 2 VBF+ditau, 
// 12 => run 3 VBF+ditau, 13 => run 3 double PF jets + ditau, 14 => di-tau + PFJet, 15 => Displaced Tau, 16 => Monitoring, 
// 17 => regional paths, 18 => L1 seeded paths, 19 => 1 prong tau paths for Tau

// 7 -> deeptau di-tau, 128
// 8 -> e-tau, 256
// 9 -> mu-tau, 512
// 10 -> single-tau, 1024
// 12 -> run 3 VBF+ditau, 4096
// 14 -> di-tau + PFJet, 16384

// new sample
// 0 => Loose, 1 => Medium, 2 => Tight, 3 => DeepTau no spec WP, 4 => ChargedIso, 5 => HPS, 
// 6 => e-tau inside filter, 7 => mu-tau inside filter, 8 => single-tau inside filter, 9 => VBF matching, 10 => di-tau, 11 => e-tau, 12 => mu-tau, 
// 13 => di-tau + PFJet, 14 => e-tau displaced, 15 => mu-tau displaced, 16 => di-tau displaced, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, 
// 19 => 'Monitoring di-tau + PFJet, 20 => 'Monitoring muTau displaced, 21 => OneProng, 22 => DiJetCorr, 23 => OverlapFilter, 24 => Dxy, 25 => MatchL1HLT, 
// 26 => MatchL1HLT, 27 => VBF + DoubleTau Monitoring, 28 => For matching to monitoring trigger for 20 GeV tau leg of VBF triggers, 29 => single PF-tau inside filter for Tau

// 10, 11, 12 -> ditau 1024,  etau and mutau 4096
// 6, 7, 8 -> 64, 128, 256
// 13 ditau PFJet 8192, 19 'Monitoring di-tau + PFJet,
// 28 For matching to monitoring trigger for 20 GeV tau leg of VBF triggers ? 268435456
// 29 -> single PF-tau inside filter for Tau 536870912

bool PassMuTauTrig2022(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & (1<<7)) != 0 && (trig_bits[it] & (1<<26)) != 0 && trig_id[it] == 15){ 
          return true;
      }
    }
  }
  return false;
}
bool PassEleTauTrig2022(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & (1<<7)) != 0 && (trig_bits[it] & (1<<26)) != 0 && trig_id[it] == 15 && trig_pt[it] > 30 && trig_l1iso[it] > 0 && trig_l1pt[it] > 26){ 
          return true;
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<17)) != 0 && (trig_bits[it] & (1<<18)) == 0 && trig_id[it] == 15 && trig_pt[it] > 35){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        }
      }
    }
  }
  return false;
}
/// checking about ditau path
bool PassDiTauTrig2022_withptiso_nobitcut(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if(trig_id[it] == 15 && trig_pt[it] > 35){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_withptiso_bit1(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && trig_id[it] == 15 && trig_pt[it] > 35){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_withptiso_bit1_bit17(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<17)) != 0 && trig_id[it] == 15 && trig_pt[it] > 35){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_withptiso_bit1_bit17_0bit18(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<17)) != 0 && (trig_bits[it] & (1<<18)) == 0 && trig_id[it] == 15 && trig_pt[it] > 35){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        }
      }
    }
  }
  return false;
}

bool PassDiTauTrig2022_nobit17(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<18)) == 0 && trig_id[it] == 15 && trig_pt[it] > 35){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_bit1_bit17_0bit18_TrigPt35(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<17)) != 0 && (trig_bits[it] & (1<<18)) == 0 && trig_id[it] == 15 && trig_pt[it] > 35){ 
        // if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        // }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_bit1_bit17_0bit18_TrigPt35_l1Iso(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<17)) != 0 && (trig_bits[it] & (1<<18)) == 0 && trig_id[it] == 15 && trig_pt[it] > 35){ 
        if ( (trig_l1iso[it] > 0)) {
          return true;
        }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_bit1_bit17_0bit18_TrigPt35_l1Iso_ORl1pt70(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<17)) != 0 && (trig_bits[it] & (1<<18)) == 0 && trig_id[it] == 15 && trig_pt[it] > 35){ 
        if ( (trig_l1iso[it] > 0) || trig_l1pt[it] > 70 ) {
          return true;
        }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_bit1_bit17_0bit18_TrigPt35_l1Iso_ANDl1pt32(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<17)) != 0 && (trig_bits[it] & (1<<18)) == 0 && trig_id[it] == 15 && trig_pt[it] > 35){ 
        if ( (trig_l1iso[it] > 0) && trig_l1pt[it] > 32 ) {
          return true;
        }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_bit1(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && trig_id[it] == 15){ 
        // if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        // }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_bit1_bit17(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<17)) != 0 && trig_id[it] == 15){ 
        // if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        // }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_bit17(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<17)) != 0 && trig_id[it] == 15){ 
        // if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        // }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_0bit18(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<18)) == 0 && trig_id[it] == 15){ 
        // if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        // }
      }
    }
  }
  return false;
}
bool PassDiTauTrig2022_bit1_bit17_0bit18(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<17)) != 0 && (trig_bits[it] & (1<<18)) == 0 && trig_id[it] == 15){ 
        // if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) || trig_l1pt[it] > 70 ) {
          return true;
        // }
      }
    }
  }
  return false;
}
bool PassDiTauJetTrig2022(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<17)) != 0 && (trig_bits[it] & (1<<18)) == 0 && trig_id[it] == 15){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 32) )  
          return true;
      }
    }
  }
  return false;
}
bool PassSingleTauTrig2022(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & 536870912) != 0 && trig_id[it] == 15){ 
          return true;
      }
    }
  }
  return false;
}
bool PassVBFDiTauTrig2022(UInt_t ntrig,Vec_i trig_id,Vec_i trig_bits,Vec_t trig_pt,Vec_t trig_eta,Vec_t trig_phi,float tau_pt,float tau_eta,float tau_phi){
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if((trig_bits[it] & 268435456) != 0 && trig_id[it] == 15){ 
          return true;
      }
    }
  }
  return false;
}