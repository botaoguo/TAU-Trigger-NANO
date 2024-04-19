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
using cRVecF = const ROOT::RVecF &;
using cRVecI = const ROOT::RVecI &;
using cRVecC = const ROOT::RVecC &;
using cRVecU = const ROOT::RVecU &;

float deltaR(float eta_1, float eta_2, float phi_1, float phi_2){
   const float deta = eta_1 - eta_2;
   const float dphi = ROOT::Math::VectorUtil::Phi_mpi_pi(phi_1 - phi_2);
   const float dRsq = std::pow(deta,2) + std::pow(dphi,2);

   return sqrt(dRsq);
}

bool SelectTau(cRVecF tau_pt, cRVecF tau_eta, cRVecC tau_idDeepTau2018v2p5VSe, cRVecC tau_idDeepTau2018v2p5VSmu, cRVecC tau_idDeepTau2018v2p5VSjet)
{
  int tau_num = 0;
  for(auto i =0; i < tau_pt.size(); i++) {
    if (tau_pt[i] > 20 && abs(tau_eta[i]) < 2.3) {
      if (tau_idDeepTau2018v2p5VSe[i] >= 2 && tau_idDeepTau2018v2p5VSmu[i] >= 1 && tau_idDeepTau2018v2p5VSjet[i] >=5) {
        tau_num +=1;
      }
    }
  }
  if (tau_num >= 2) return true;
  else return false;
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
