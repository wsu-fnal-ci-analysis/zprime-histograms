//==============================================================
//          Analysis code for Z' boson to e+ e- analysis       =  
//          In this code we select the HEEP events             = 
//          To run over MINIAOD MC with fixed trigger          = 
//                  Author:  Sherif Elgammal                   = 
//                                                             = 
//                       01/01/2016                            = 
//==============================================================
#ifndef ZprimeMuMuPatMiniAodNewData_h
#define ZprimeMuMuPatMiniAodNewData_h
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include <math.h>
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include <time.h>
#include <iostream>
#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <sstream>
#include "TVector3.h"
// Header file for the classes stored in the TTree if any.
#include <vector>
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class ZprimeEleElePatMiniAodNewData {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event_runNo;
   Int_t           event_evtNo;
   Int_t           event_lumi;
   Int_t           event_bunch;
   vector<int>     *iGenJet;
   vector<int>     *idGenJet;
   vector<int>     *statusGenJet;
   vector<int>     *chargeGenJet;
   vector<float>   *ptGenJet;
   vector<float>   *etaGenJet;
   vector<float>   *phiGenJet;
   vector<int>     *iGen;
   vector<int>     *idGen;
   vector<int>     *statusGen;
   vector<float>   *ptGen;
   vector<float>   *etaGen;
   vector<float>   *phiGen;
   vector<int>     *chargeGen;
   vector<float>   *EnergyGen;
   vector<float>   *pxGen;
   vector<float>   *pyGen;
   vector<float>   *pzGen;
   vector<int>     *Mu_nbMuon;
   vector<bool>    *Mu_isTightMuon;
   vector<bool>    *Mu_isLooseMuon;
   vector<bool>    *Mu_isGlobalMuon;
   vector<bool>    *Mu_isPF;
   vector<bool>    *Mu_isTrackerMuon;
   vector<float>   *Mu_et;
   vector<float>   *Mu_en;
   vector<float>   *Mu_pt;
   vector<float>   *Mu_eta;
   vector<float>   *Mu_phi;
   vector<float>   *Mu_charge;
   vector<float>   *Mu_ptTunePMuonBestTrack;
   vector<float>   *Mu_pxTunePMuonBestTrack;
   vector<float>   *Mu_pyTunePMuonBestTrack;
   vector<float>   *Mu_pzTunePMuonBestTrack;
   vector<float>   *Mu_pTunePMuonBestTrack;
   vector<float>   *Mu_etaTunePMuonBestTrack;
   vector<float>   *Mu_phiTunePMuonBestTrack;
   vector<float>   *Mu_thetaTunePMuonBestTrack;
   vector<float>   *Mu_chargeTunePMuonBestTrack;
   vector<float>   *Mu_dPToverPTTunePMuonBestTrack;
   vector<float>   *Mu_absdxyTunePMuonBestTrack;
   vector<float>   *Mu_absdzTunePMuonBestTrack;
   vector<float>   *Mu_ptInnerTrack;
   vector<float>   *Mu_pxInnerTrack;
   vector<float>   *Mu_pyInnerTrack;
   vector<float>   *Mu_pzInnerTrack;
   vector<float>   *Mu_pInnerTrack;
   vector<float>   *Mu_etaInnerTrack;
   vector<float>   *Mu_phiInnerTrack;
   vector<float>   *Mu_thetaInnerTrack;
   vector<float>   *Mu_chargeInnerTrack;
   vector<float>   *Mu_dPToverPTInnerTrack;
   vector<float>   *Mu_absdxyInnerTrack;
   vector<float>   *Mu_absdzInnerTrack;
   vector<float>   *Mu_normalizedChi2;
   vector<float>   *Mu_absdxy;
   vector<float>   *Mu_absdz;
   vector<float>   *Mu_vtxMass;
   vector<float>   *Mu_vtxNormChi2;
   vector<float>   *Mu_vtxMassLept;
   vector<int>     *Mu_numberOfMatchedStations;
   vector<int>     *Mu_numberOfValidPixelHits;
   vector<int>     *Mu_numberOfValidMuonHits;
   vector<int>     *Mu_numberOftrackerLayersWithMeasurement;
   vector<float>   *Mu_emIso;
   vector<float>   *Mu_hadIso;
   vector<float>   *Mu_trackiso;
   vector<float>   *Mu_pfSumChargedHadronPt;
   vector<float>   *Mu_pfSumNeutralHadronEt;
   vector<float>   *Mu_PFSumPhotonEt;
   vector<float>   *Mu_pfSumPUPt;
   vector<int>     *Mu_nbofpv;
   vector<float>   *Mu_patDeltaBeta;
   vector<int>     *HLT_nb;
   vector<string>  *HLT_name;
   vector<bool>    *HLT_isaccept;
   vector<int>     *HLTObj_nbObj;
   vector<float>   *HLTObj_pt;
   vector<float>   *HLTObj_eta;
   vector<float>   *HLTObj_phi;
   vector<string>  *HLTObj_collection;
   Double_t        GenMet_pt;
   Double_t        PFMet_et_cor;
   Double_t        PFMet_pt_cor;
   Double_t        PFMet_phi_cor;
   Double_t        PFMet_en_cor;
   Double_t        PFMet_px_cor;
   Double_t        PFMet_py_cor;
   Double_t        PFMet_pz_cor;
   Double_t        PFMet_sumEt_cor;
   Double_t        CaloMet_pt;
   Double_t        CaloMet_phi;
   Double_t        CaloMet_sumEt;
   Double_t        PFMet_shiftedPt_JetEnUp;
   Double_t        PFMet_shiftedPt_JetEnDown;
   vector<int>     *jet_nb;
   vector<float>   *jet_charge;
   vector<float>   *jet_et;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_en;
   vector<float>   *jet_theta;
   vector<float>   *jet_beta;
   vector<float>   *jet_pileup_mva_disc;
   vector<int>     *Nb_bDiscriminators;
   vector<float>   *jet_btag_pt;
   vector<float>   *jet_btag_eta;
   vector<float>   *jet_btag_phi;
   vector<int>     *jet_btag_flavor;
   vector<float>   *jet_btag_pfCSVv2IVF_discriminator;
   vector<int>     *Nb_taus;
   vector<float>   *Tau_pt;
   vector<float>   *Tau_eta;
   vector<float>   *Tau_phi;
   vector<int>     *Tau_id;
   vector<float>   *Tau_LooseCombinedIsolationDeltaBetaCorr3Hits;
   Int_t           pfphoton_size;
   vector<float>   *pfphoton_pt;
   vector<float>   *pfphoton_eta;
   vector<float>   *pfphoton_phi;
   vector<float>   *pfphoton_theta;
   Int_t           num_PU_vertices;
   Int_t           PU_BunchCrossing;
   Int_t           num_PU_gen_vertices;
   Float_t         Rho;
   vector<float>   *MC_weighting;
   vector<int>     *Ele_nbElectrons;
   vector<bool>    *Ele_isEcalDrivenSeed;
   vector<bool>    *Ele_isPassConversionVeto;
   vector<int>     *Ele_charge;
   vector<int>     *Ele_nbOfMissingHits;
   vector<int>     *Ele_nbVtx;
   vector<float>   *Ele_Et;
   vector<float>   *Ele_EtFromCaloEn;
   vector<float>   *Ele_pt;
   vector<float>   *Ele_thetaSC;
   vector<float>   *Ele_etaSC;
   vector<float>   *Ele_phiSC;
   vector<float>   *Ele_energySC;
   vector<float>   *Ele_preshowerEnergySC;
   vector<float>   *Ele_thetaTrack;
   vector<float>   *Ele_etaTrack;
   vector<float>   *Ele_phiTrack;
   vector<float>   *Ele_hadronicOverEm;
   vector<float>   *Ele_deltaEtaInSeedCluster;
   vector<float>   *Ele_deltaPhiInSeedCluster;
   vector<float>   *Ele_deltaEtaInSC;
   vector<float>   *Ele_deltaPhiInSC;
   vector<float>   *Ele_sigmaIetaIeta;
   vector<float>   *Ele_e2x5Max;
   vector<float>   *Ele_e1x5;
   vector<float>   *Ele_frac15;
   vector<float>   *Ele_frac51;
   vector<float>   *Ele_e5x5;
   vector<float>   *Ele3x3;
   vector<float>   *Ele_e2x5MaxOver5x5;
   vector<float>   *Ele_e1x5Over5x5;
   vector<float>   *Ele_sigmaIetaIetaFull5x5;
   vector<float>   *Ele_e2x5MaxFull5x5;
   vector<float>   *Ele_e1x5Full5x5;
   vector<float>   *Ele_e5x5Full5x5;
   vector<float>   *Ele_e2x5MaxOver5x5Full5x5;
   vector<float>   *Ele_e1x5Over5x5Full5x5;
   vector<float>   *Ele_e2x5Right;
   vector<float>   *Ele_e2x5Left;
   vector<float>   *Ele_e2x5Top;
   vector<float>   *Ele_e2x5Bottom;
   vector<float>   *Ele_eMax;
   vector<float>   *Ele_eRight;
   vector<float>   *Ele_eLeft;
   vector<float>   *Ele_eTop;
   vector<float>   *Ele_eBottom;
   vector<float>   *Ele_dxy;
   vector<float>   *Ele_dz;
   vector<float>   *Ele_rhoIso;
   vector<float>   *Ele_fbrem;
   vector<float>   *Ele_EoverP;
   vector<float>   *Ele_Xposition;
   vector<float>   *Ele_Yposition;
   vector<float>   *Ele_EcalPlusHcald1iso;
   vector<float>   *Ele_dr03EcalRecHitSumEt;
   vector<float>   *Ele_dr03HcalDepth1TowerSumEt;
   vector<float>   *Ele_dr03HcalDepth1TowerSumEtBc;
   vector<float>   *Ele_hcalDepth1OverEcal;
   vector<float>   *Ele_hcalDepth2OverEcal;
   vector<float>   *Ele_dr03HcalDepth2TowerSumEt;
   vector<float>   *Ele_hcalDepth2TowerSumEtNoVeto;
   vector<float>   *Ele_hcalDepth1TowerSumEtNoVeto;
   vector<float>   *Ele_pfSumPhotonEt;
   vector<float>   *Ele_pfSumChargedHadronPt;
   vector<float>   *Ele_pfSumNeutralHadronEt;
   vector<float>   *Ele_pfSumPUPt;
   vector<float>   *Ele_pfDeltaBeta;
   vector<int>     *Ele_rawId;
   vector<float>   *Ele_x;
   vector<float>   *Ele_y;
   vector<float>   *Ele_z;
   vector<float>   *Ele_zTrackPositionAtVtx;
   vector<int>     *Ele_ieta;
   vector<float>   *Ele_phiWidth;
   vector<float>   *Ele_etaWidth;
   vector<float>   *Ele_dr03TkSumPt;
   vector<float>   *Ele_dr03TkSumPt_corrected;

   // List of branches
   TBranch        *b_event_runNo;   //!
   TBranch        *b_event_evtNo;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_event_bunch;   //!
   TBranch        *b_iGenJet;   //!
   TBranch        *b_idGenJet;   //!
   TBranch        *b_statusGenJet;   //!
   TBranch        *b_chargeGenJet;   //!
   TBranch        *b_ptGenJet;   //!
   TBranch        *b_etaGenJet;   //!
   TBranch        *b_phiGenJet;   //!
   TBranch        *b_iGen;   //!
   TBranch        *b_idGen;   //!
   TBranch        *b_statusGen;   //!
   TBranch        *b_ptGen;   //!
   TBranch        *b_etaGen;   //!
   TBranch        *b_phiGen;   //!
   TBranch        *b_chargeGen;   //!
   TBranch        *b_EnergyGen;   //!
   TBranch        *b_pxGen;   //!
   TBranch        *b_pyGen;   //!
   TBranch        *b_pzGen;   //!
   TBranch        *b_Mu_nbMuon;   //!
   TBranch        *b_Mu_isTightMuon;   //!
   TBranch        *b_Mu_isLooseMuon;   //!
   TBranch        *b_Mu_isGlobalMuon;   //!
   TBranch        *b_Mu_isPF;   //!
   TBranch        *b_Mu_isTrackerMuon;   //!
   TBranch        *b_Mu_et;   //!
   TBranch        *b_Mu_en;   //!
   TBranch        *b_Mu_pt;   //!
   TBranch        *b_Mu_eta;   //!
   TBranch        *b_Mu_phi;   //!
   TBranch        *b_Mu_charge;   //!
   TBranch        *b_Mu_ptTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pxTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pyTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pzTunePMuonBestTrack;   //!
   TBranch        *b_Mu_pTunePMuonBestTrack;   //!
   TBranch        *b_Mu_etaTunePMuonBestTrack;   //!
   TBranch        *b_Mu_phiTunePMuonBestTrack;   //!
   TBranch        *b_Mu_thetaTunePMuonBestTrack;   //!
   TBranch        *b_Mu_chargeTunePMuonBestTrack;   //!
   TBranch        *b_Mu_dPToverPTTunePMuonBestTrack;   //!
   TBranch        *b_Mu_absdxyTunePMuonBestTrack;   //!
   TBranch        *b_Mu_absdzTunePMuonBestTrack;   //!
   TBranch        *b_Mu_ptInnerTrack;   //!
   TBranch        *b_Mu_pxInnerTrack;   //!
   TBranch        *b_Mu_pyInnerTrack;   //!
   TBranch        *b_Mu_pzInnerTrack;   //!
   TBranch        *b_Mu_pInnerTrack;   //!
   TBranch        *b_Mu_etaInnerTrack;   //!
   TBranch        *b_Mu_phiInnerTrack;   //!
   TBranch        *b_Mu_thetaInnerTrack;   //!
   TBranch        *b_Mu_chargeInnerTrack;   //!
   TBranch        *b_Mu_dPToverPTInnerTrack;   //!
   TBranch        *b_Mu_absdxyInnerTrack;   //!
   TBranch        *b_Mu_absdzInnerTrack;   //!
   TBranch        *b_Mu_normalizedChi2;   //!
   TBranch        *b_Mu_absdxy;   //!
   TBranch        *b_Mu_absdz;   //!
   TBranch        *b_Mu_vtxMass;   //!
   TBranch        *b_Mu_vtxNormChi2;   //!
   TBranch        *b_Mu_vtxMassLept;   //!
   TBranch        *b_Mu_numberOfMatchedStations;   //!
   TBranch        *b_Mu_numberOfValidPixelHits;   //!
   TBranch        *b_Mu_numberOfValidMuonHits;   //!
   TBranch        *b_Mu_numberOftrackerLayersWithMeasurement;   //!
   TBranch        *b_Mu_emIso;   //!
   TBranch        *b_Mu_hadIso;   //!
   TBranch        *b_Mu_trackiso;   //!
   TBranch        *b_Mu_pfSumChargedHadronPt;   //!
   TBranch        *b_Mu_pfSumNeutralHadronEt;   //!
   TBranch        *b_Mu_PFSumPhotonEt;   //!
   TBranch        *b_Mu_pfSumPUPt;   //!
   TBranch        *b_Mu_nbofpv;   //!
   TBranch        *b_Mu_patDeltaBeta;   //!
   TBranch        *b_HLT_nb;   //!
   TBranch        *b_HLT_name;   //!
   TBranch        *b_HLT_isaccept;   //!
   TBranch        *b_HLTObj_nbObj;   //!
   TBranch        *b_HLTObj_pt;   //!
   TBranch        *b_HLTObj_eta;   //!
   TBranch        *b_HLTObj_phi;   //!
   TBranch        *b_HLTObj_collection;   //!
   TBranch        *b_GenMet_pt;   //!
   TBranch        *b_PFMet_et_cor;   //!
   TBranch        *b_PFMet_pt_cor;   //!
   TBranch        *b_PFMet_phi_cor;   //!
   TBranch        *b_PFMet_en_cor;   //!
   TBranch        *b_PFMet_px_cor;   //!
   TBranch        *b_PFMet_py_cor;   //!
   TBranch        *b_PFMet_pz_cor;   //!
   TBranch        *b_PFMet_sumEt_cor;   //!
   TBranch        *b_CaloMet_pt;   //!
   TBranch        *b_CaloMet_phi;   //!
   TBranch        *b_CaloMet_sumEt;   //!
   TBranch        *b_PFMet_shiftedPt_JetEnUp;   //!
   TBranch        *b_PFMet_shiftedPt_JetEnDown;   //!
   TBranch        *b_jet_nb;   //!
   TBranch        *b_jet_charge;   //!
   TBranch        *b_jet_et;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_en;   //!
   TBranch        *b_jet_theta;   //!
   TBranch        *b_jet_beta;   //!
   TBranch        *b_jet_pileup_mva_disc;   //!
   TBranch        *b_Nb_bDiscriminators;   //!
   TBranch        *b_jet_btag_pt;   //!
   TBranch        *b_jet_btag_eta;   //!
   TBranch        *b_jet_btag_phi;   //!
   TBranch        *b_jet_btag_flavor;   //!
   TBranch        *b_jet_btag_pfCSVv2IVF_discriminator;   //!
   TBranch        *b_Nb_taus;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_id;   //!
   TBranch        *b_Tau_LooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_pfphoton_size;   //!
   TBranch        *b_pfphoton_pt;   //!
   TBranch        *b_pfphoton_eta;   //!
   TBranch        *b_pfphoton_phi;   //!
   TBranch        *b_pfphoton_theta;   //!
   TBranch        *b_num_PU_vertices;   //!
   TBranch        *b_PU_BunchCrossing;   //!
   TBranch        *b_num_PU_gen_vertices;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_MC_weighting;   //!
   TBranch        *b_Ele_nbElectrons;   //!
   TBranch        *b_Ele_isEcalDrivenSeed;   //!
   TBranch        *b_Ele_isPassConversionVeto;   //!
   TBranch        *b_Ele_charge;   //!
   TBranch        *b_Ele_nbOfMissingHits;   //!
   TBranch        *b_Ele_nbVtx;   //!
   TBranch        *b_Ele_Et;   //!
   TBranch        *b_Ele_EtFromCaloEn;   //!
   TBranch        *b_Ele_pt;   //!
   TBranch        *b_Ele_thetaSC;   //!
   TBranch        *b_Ele_etaSC;   //!
   TBranch        *b_Ele_phiSC;   //!
   TBranch        *b_Ele_energySC;   //!
   TBranch        *b_Ele_preshowerEnergySC;   //!
   TBranch        *b_Ele_thetaTrack;   //!
   TBranch        *b_Ele_etaTrack;   //!
   TBranch        *b_Ele_phiTrack;   //!
   TBranch        *b_Ele_hadronicOverEm;   //!
   TBranch        *b_Ele_deltaEtaInSeedCluster;   //!
   TBranch        *b_Ele_deltaPhiInSeedCluster;   //!
   TBranch        *b_Ele_deltaEtaInSC;   //!
   TBranch        *b_Ele_deltaPhiInSC;   //!
   TBranch        *b_Ele_sigmaIetaIeta;   //!
   TBranch        *b_Ele_e2x5Max;   //!
   TBranch        *b_Ele_e1x5;   //!
   TBranch        *b_Ele_frac15;   //!
   TBranch        *b_Ele_frac51;   //!
   TBranch        *b_Ele_e5x5;   //!
   TBranch        *b_Ele3x3;   //!
   TBranch        *b_Ele_e2x5MaxOver5x5;   //!
   TBranch        *b_Ele_e1x5Over5x5;   //!
   TBranch        *b_Ele_sigmaIetaIetaFull5x5;   //!
   TBranch        *b_Ele_e2x5MaxFull5x5;   //!
   TBranch        *b_Ele_e1x5Full5x5;   //!
   TBranch        *b_Ele_e5x5Full5x5;   //!
   TBranch        *b_Ele_e2x5MaxOver5x5Full5x5;   //!
   TBranch        *b_Ele_e1x5Over5x5Full5x5;   //!
   TBranch        *b_Ele_e2x5Right;   //!
   TBranch        *b_Ele_e2x5Left;   //!
   TBranch        *b_Ele_e2x5Top;   //!
   TBranch        *b_Ele_e2x5Bottom;   //!
   TBranch        *b_Ele_eMax;   //!
   TBranch        *b_Ele_eRight;   //!
   TBranch        *b_Ele_eLeft;   //!
   TBranch        *b_Ele_eTop;   //!
   TBranch        *b_Ele_eBottom;   //!
   TBranch        *b_Ele_dxy;   //!
   TBranch        *b_Ele_dz;   //!
   TBranch        *b_Ele_rhoIso;   //!
   TBranch        *b_Ele_fbrem;   //!
   TBranch        *b_Ele_EoverP;   //!
   TBranch        *b_Ele_Xposition;   //!
   TBranch        *b_Ele_Yposition;   //!
   TBranch        *b_Ele_EcalPlusHcald1iso;   //!
   TBranch        *b_Ele_dr03EcalRecHitSumEt;   //!
   TBranch        *b_Ele_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Ele_dr03HcalDepth1TowerSumEtBc;   //!
   TBranch        *b_Ele_hcalDepth1OverEcal;   //!
   TBranch        *b_Ele_hcalDepth2OverEcal;   //!
   TBranch        *b_Ele_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_Ele_hcalDepth2TowerSumEtNoVeto;   //!
   TBranch        *b_Ele_hcalDepth1TowerSumEtNoVeto;   //!
   TBranch        *b_Ele_pfSumPhotonEt;   //!
   TBranch        *b_Ele_pfSumChargedHadronPt;   //!
   TBranch        *b_Ele_pfSumNeutralHadronEt;   //!
   TBranch        *b_Ele_pfSumPUPt;   //!
   TBranch        *b_Ele_pfDeltaBeta;   //!
   TBranch        *b_Ele_rawId;   //!
   TBranch        *b_Ele_x;   //!
   TBranch        *b_Ele_y;   //!
   TBranch        *b_Ele_z;   //!
   TBranch        *b_Ele_zTrackPositionAtVtx;   //!
   TBranch        *b_Ele_ieta;   //!
   TBranch        *b_Ele_phiWidth;   //!
   TBranch        *b_Ele_etaWidth;   //!
   TBranch        *b_Ele_dr03TkSumPt;   //!
   TBranch        *b_Ele_dr03TkSumPt_corrected;   //!

   ZprimeEleElePatMiniAodNewData(TTree *tree=0);
   virtual ~ZprimeEleElePatMiniAodNewData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   bool SelectFirstEle(float &ETele1,float &Enele1,float &Etaele1,
		       float &Phiele1,int &ChargeEle1,float &EtaSCele1,
		       float &PhiSCele1,unsigned &FlagEle1);
   bool SelectSecondEle(int ChargeEle1,unsigned FlagEle1,float ETele1,float Etaele1,float Phiele1,
			float &ETele2,float &Enele2,float &Etaele2,float &Phiele2,int &ChargeEle2,
			float &EtaSCele2,float &PhiSCele2);
   float Mass(float Pt1,float Eta1,float Phi1,float En1,
              float Pt2,float Eta2,float Phi2,float En2);
   void PlotRecoInfo(float MassEle,float etaEle1,float etaEle2);
   bool SelectFirstGenEle(float &ETEle1,float &PhiSCEle1,
                         float &EtaSCEle1,float &EnEle1,
                         int &IDele1,int &Statele1,
                         unsigned &GenFlag1);
   bool SelectSecondGenEle(unsigned GenFlag1,float ETEle1,float &ETEle2,float &PhiSCEle2,
                          float &EtaSCEle2,float &EnEle2,int &IDele2,int &Statele2);
   bool GenRecoMatchEle(float RecoEta1,float RecoPhi1);
   bool RecoHLTEleMatching(float RecoEta,float RecoPhi);
   bool isPassHLT();
   void CosThetaCollinSoper(float Et1,float Eta1,float Phi1,float En1,
			    float Et2,float Eta2,float Phi2,float En2,
			    float ChargeEle1,float RecoMass);
   float delR(float eta1,float phi1,float eta2,float phi2);
   //================================================================================
   float HLT_pt,HLT_eta,HLT_phi;
   int weight;
   float Etele1,Enele1,EtaTrakele1,PhiTrakele1,EtaSCele1,PhiSCele1;
   int Chargeele1;
   float Etele2,Enele2,EtaTrakele2,PhiTrakele2,EtaSCele2,PhiSCele2;
   int Chargeele2;
   float ptEffCut;
   float PtDYTRecMu1,PtDYTRecMu2,PtRecTunePMu1,PtRecTunePMu2,PtRecMuBestTrack1,PtRecMuBestTrack2;
   float RecoHLTMatchingDeltaRcut,deltaRcut,minMassCut,maxMassCut;
   float DiEleMass;
   float mPtGen1,mPhiGen1,mEtaGen1,mEnGen1;
   unsigned mGenFlag1;
   float mPtGen2,mPhiGen2,mEtaGen2,mEnGen2;
   int ChargeRecMu1,ChargeRecMu2;
   unsigned flagmu1;
   unsigned flag1;
   float PtRecTunePMuBestTrack1,EnRecMu1,EtaRecMu1,PhiRecMu1;
   float PtRecTunePMuBestTrack2,EnRecMu2,EtaRecMu2,PhiRecMu2;
   float pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,dxyRecMu1;
   float pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2;
   float genET1,genPhi1,genEta1,genEn1;
   int genID1,genStat1;
   float genET2,genPhi2,genEta2,genEn2;
   int genID2,genStat2;
   float MassGen,RecoMass;
   int NbGen,NbReco;
   int nbTP,nbTT,nbTF;
   float TagProbeEtaCut;
   float Eff;
   float MassCutMin,MassCutMax;
   float MassResolution;
   float EtaCut;
   int NbFireHLT; 
   ofstream output_txt;    
   TH1F* h1_CosAngleCollinSoperCorrect60Mass120_;
   TH1F* h1_CosAngleCollinSoperCorrect120Mass300_;
   TH1F* h1_CosAngleCollinSoperCorrect300Mass700_;
   TH1F* h1_CosAngleCollinSoperCorrect700Mass3000_;
   TH1F* h1_CosAngleCollinSoperCorrect4900Mass5100_;
   TH1F* h1_absCosAngleCollinSoperCorrect4500Mass5500_;
   TH1F *h1_ZprimeRecomassBinWidth_;
   TH1F* h1_ZprimeRecomassBinWidthBB_;
   TH1F* h1_ZprimeRecomassBinWidthBE_;
   TH1F* h1_ZprimeRecomassBinWidthEE_;
   TH1F* h1_ZprimeRecomass60to120BE_;
   TH1F* h1_ZprimeRecomass60to120EE_;
   TH1F* h1_ZprimeRecomass60to120BB_;
   TH1F* h1_ZprimeRecomass60to120_;
};

#endif

#ifdef ZprimeEleElePatMiniAodNewData_cxx
ZprimeEleElePatMiniAodNewData::ZprimeEleElePatMiniAodNewData(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MC.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MC.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

ZprimeEleElePatMiniAodNewData::~ZprimeEleElePatMiniAodNewData()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ZprimeEleElePatMiniAodNewData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZprimeEleElePatMiniAodNewData::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ZprimeEleElePatMiniAodNewData::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   iGenJet = 0;
   idGenJet = 0;
   statusGenJet = 0;
   chargeGenJet = 0;
   ptGenJet = 0;
   etaGenJet = 0;
   phiGenJet = 0;
   iGen = 0;
   idGen = 0;
   statusGen = 0;
   ptGen = 0;
   etaGen = 0;
   phiGen = 0;
   chargeGen = 0;
   EnergyGen = 0;
   pxGen = 0;
   pyGen = 0;
   pzGen = 0;
   Mu_nbMuon = 0;
   Mu_isTightMuon = 0;
   Mu_isLooseMuon = 0;
   Mu_isGlobalMuon = 0;
   Mu_isPF = 0;
   Mu_isTrackerMuon = 0;
   Mu_et = 0;
   Mu_en = 0;
   Mu_pt = 0;
   Mu_eta = 0;
   Mu_phi = 0;
   Mu_charge = 0;
   Mu_ptTunePMuonBestTrack = 0;
   Mu_pxTunePMuonBestTrack = 0;
   Mu_pyTunePMuonBestTrack = 0;
   Mu_pzTunePMuonBestTrack = 0;
   Mu_pTunePMuonBestTrack = 0;
   Mu_etaTunePMuonBestTrack = 0;
   Mu_phiTunePMuonBestTrack = 0;
   Mu_thetaTunePMuonBestTrack = 0;
   Mu_chargeTunePMuonBestTrack = 0;
   Mu_dPToverPTTunePMuonBestTrack = 0;
   Mu_absdxyTunePMuonBestTrack = 0;
   Mu_absdzTunePMuonBestTrack = 0;
   Mu_ptInnerTrack = 0;
   Mu_pxInnerTrack = 0;
   Mu_pyInnerTrack = 0;
   Mu_pzInnerTrack = 0;
   Mu_pInnerTrack = 0;
   Mu_etaInnerTrack = 0;
   Mu_phiInnerTrack = 0;
   Mu_thetaInnerTrack = 0;
   Mu_chargeInnerTrack = 0;
   Mu_dPToverPTInnerTrack = 0;
   Mu_absdxyInnerTrack = 0;
   Mu_absdzInnerTrack = 0;
   Mu_normalizedChi2 = 0;
   Mu_absdxy = 0;
   Mu_absdz = 0;
   Mu_vtxMass = 0;
   Mu_vtxNormChi2 = 0;
   Mu_vtxMassLept = 0;
   Mu_numberOfMatchedStations = 0;
   Mu_numberOfValidPixelHits = 0;
   Mu_numberOfValidMuonHits = 0;
   Mu_numberOftrackerLayersWithMeasurement = 0;
   Mu_emIso = 0;
   Mu_hadIso = 0;
   Mu_trackiso = 0;
   Mu_pfSumChargedHadronPt = 0;
   Mu_pfSumNeutralHadronEt = 0;
   Mu_PFSumPhotonEt = 0;
   Mu_pfSumPUPt = 0;
   Mu_nbofpv = 0;
   Mu_patDeltaBeta = 0;
   HLT_nb = 0;
   HLT_name = 0;
   HLT_isaccept = 0;
   HLTObj_nbObj = 0;
   HLTObj_pt = 0;
   HLTObj_eta = 0;
   HLTObj_phi = 0;
   HLTObj_collection = 0;
   jet_nb = 0;
   jet_charge = 0;
   jet_et = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_en = 0;
   jet_theta = 0;
   jet_beta = 0;
   jet_pileup_mva_disc = 0;
   Nb_bDiscriminators = 0;
   jet_btag_pt = 0;
   jet_btag_eta = 0;
   jet_btag_phi = 0;
   jet_btag_flavor = 0;
   jet_btag_pfCSVv2IVF_discriminator = 0;
   Nb_taus = 0;
   Tau_pt = 0;
   Tau_eta = 0;
   Tau_phi = 0;
   Tau_id = 0;
   Tau_LooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   pfphoton_pt = 0;
   pfphoton_eta = 0;
   pfphoton_phi = 0;
   pfphoton_theta = 0;
   MC_weighting = 0;
   Ele_nbElectrons = 0;
   Ele_isEcalDrivenSeed = 0;
   Ele_isPassConversionVeto = 0;
   Ele_charge = 0;
   Ele_nbOfMissingHits = 0;
   Ele_nbVtx = 0;
   Ele_Et = 0;
   Ele_EtFromCaloEn = 0;
   Ele_pt = 0;
   Ele_thetaSC = 0;
   Ele_etaSC = 0;
   Ele_phiSC = 0;
   Ele_energySC = 0;
   Ele_preshowerEnergySC = 0;
   Ele_thetaTrack = 0;
   Ele_etaTrack = 0;
   Ele_phiTrack = 0;
   Ele_hadronicOverEm = 0;
   Ele_deltaEtaInSeedCluster = 0;
   Ele_deltaPhiInSeedCluster = 0;
   Ele_deltaEtaInSC = 0;
   Ele_deltaPhiInSC = 0;
   Ele_sigmaIetaIeta = 0;
   Ele_e2x5Max = 0;
   Ele_e1x5 = 0;
   Ele_frac15 = 0;
   Ele_frac51 = 0;
   Ele_e5x5 = 0;
   Ele3x3 = 0;
   Ele_e2x5MaxOver5x5 = 0;
   Ele_e1x5Over5x5 = 0;
   Ele_sigmaIetaIetaFull5x5 = 0;
   Ele_e2x5MaxFull5x5 = 0;
   Ele_e1x5Full5x5 = 0;
   Ele_e5x5Full5x5 = 0;
   Ele_e2x5MaxOver5x5Full5x5 = 0;
   Ele_e1x5Over5x5Full5x5 = 0;
   Ele_e2x5Right = 0;
   Ele_e2x5Left = 0;
   Ele_e2x5Top = 0;
   Ele_e2x5Bottom = 0;
   Ele_eMax = 0;
   Ele_eRight = 0;
   Ele_eLeft = 0;
   Ele_eTop = 0;
   Ele_eBottom = 0;
   Ele_dxy = 0;
   Ele_dz = 0;
   Ele_rhoIso = 0;
   Ele_fbrem = 0;
   Ele_EoverP = 0;
   Ele_Xposition = 0;
   Ele_Yposition = 0;
   Ele_EcalPlusHcald1iso = 0;
   Ele_dr03EcalRecHitSumEt = 0;
   Ele_dr03HcalDepth1TowerSumEt = 0;
   Ele_dr03HcalDepth1TowerSumEtBc = 0;
   Ele_hcalDepth1OverEcal = 0;
   Ele_hcalDepth2OverEcal = 0;
   Ele_dr03HcalDepth2TowerSumEt = 0;
   Ele_hcalDepth2TowerSumEtNoVeto = 0;
   Ele_hcalDepth1TowerSumEtNoVeto = 0;
   Ele_pfSumPhotonEt = 0;
   Ele_pfSumChargedHadronPt = 0;
   Ele_pfSumNeutralHadronEt = 0;
   Ele_pfSumPUPt = 0;
   Ele_pfDeltaBeta = 0;
   Ele_rawId = 0;
   Ele_x = 0;
   Ele_y = 0;
   Ele_z = 0;
   Ele_zTrackPositionAtVtx = 0;
   Ele_ieta = 0;
   Ele_phiWidth = 0;
   Ele_etaWidth = 0;
   Ele_dr03TkSumPt = 0;
   Ele_dr03TkSumPt_corrected = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_runNo", &event_runNo, &b_event_runNo);
   fChain->SetBranchAddress("event_evtNo", &event_evtNo, &b_event_evtNo);
   fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
   fChain->SetBranchAddress("event_bunch", &event_bunch, &b_event_bunch);
   fChain->SetBranchAddress("iGenJet", &iGenJet, &b_iGenJet);
   fChain->SetBranchAddress("idGenJet", &idGenJet, &b_idGenJet);
   fChain->SetBranchAddress("statusGenJet", &statusGenJet, &b_statusGenJet);
   fChain->SetBranchAddress("chargeGenJet", &chargeGenJet, &b_chargeGenJet);
   fChain->SetBranchAddress("ptGenJet", &ptGenJet, &b_ptGenJet);
   fChain->SetBranchAddress("etaGenJet", &etaGenJet, &b_etaGenJet);
   fChain->SetBranchAddress("phiGenJet", &phiGenJet, &b_phiGenJet);
   fChain->SetBranchAddress("iGen", &iGen, &b_iGen);
   fChain->SetBranchAddress("idGen", &idGen, &b_idGen);
   fChain->SetBranchAddress("statusGen", &statusGen, &b_statusGen);
   fChain->SetBranchAddress("ptGen", &ptGen, &b_ptGen);
   fChain->SetBranchAddress("etaGen", &etaGen, &b_etaGen);
   fChain->SetBranchAddress("phiGen", &phiGen, &b_phiGen);
   fChain->SetBranchAddress("chargeGen", &chargeGen, &b_chargeGen);
   fChain->SetBranchAddress("EnergyGen", &EnergyGen, &b_EnergyGen);
   fChain->SetBranchAddress("pxGen", &pxGen, &b_pxGen);
   fChain->SetBranchAddress("pyGen", &pyGen, &b_pyGen);
   fChain->SetBranchAddress("pzGen", &pzGen, &b_pzGen);
   fChain->SetBranchAddress("Mu_nbMuon", &Mu_nbMuon, &b_Mu_nbMuon);
   fChain->SetBranchAddress("Mu_isTightMuon", &Mu_isTightMuon, &b_Mu_isTightMuon);
   fChain->SetBranchAddress("Mu_isLooseMuon", &Mu_isLooseMuon, &b_Mu_isLooseMuon);
   fChain->SetBranchAddress("Mu_isGlobalMuon", &Mu_isGlobalMuon, &b_Mu_isGlobalMuon);
   fChain->SetBranchAddress("Mu_isPF", &Mu_isPF, &b_Mu_isPF);
   fChain->SetBranchAddress("Mu_isTrackerMuon", &Mu_isTrackerMuon, &b_Mu_isTrackerMuon);
   fChain->SetBranchAddress("Mu_et", &Mu_et, &b_Mu_et);
   fChain->SetBranchAddress("Mu_en", &Mu_en, &b_Mu_en);
   fChain->SetBranchAddress("Mu_pt", &Mu_pt, &b_Mu_pt);
   fChain->SetBranchAddress("Mu_eta", &Mu_eta, &b_Mu_eta);
   fChain->SetBranchAddress("Mu_phi", &Mu_phi, &b_Mu_phi);
   fChain->SetBranchAddress("Mu_charge", &Mu_charge, &b_Mu_charge);
   fChain->SetBranchAddress("Mu_ptTunePMuonBestTrack", &Mu_ptTunePMuonBestTrack, &b_Mu_ptTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pxTunePMuonBestTrack", &Mu_pxTunePMuonBestTrack, &b_Mu_pxTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pyTunePMuonBestTrack", &Mu_pyTunePMuonBestTrack, &b_Mu_pyTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pzTunePMuonBestTrack", &Mu_pzTunePMuonBestTrack, &b_Mu_pzTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_pTunePMuonBestTrack", &Mu_pTunePMuonBestTrack, &b_Mu_pTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_etaTunePMuonBestTrack", &Mu_etaTunePMuonBestTrack, &b_Mu_etaTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_phiTunePMuonBestTrack", &Mu_phiTunePMuonBestTrack, &b_Mu_phiTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_thetaTunePMuonBestTrack", &Mu_thetaTunePMuonBestTrack, &b_Mu_thetaTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_chargeTunePMuonBestTrack", &Mu_chargeTunePMuonBestTrack, &b_Mu_chargeTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_dPToverPTTunePMuonBestTrack", &Mu_dPToverPTTunePMuonBestTrack, &b_Mu_dPToverPTTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_absdxyTunePMuonBestTrack", &Mu_absdxyTunePMuonBestTrack, &b_Mu_absdxyTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_absdzTunePMuonBestTrack", &Mu_absdzTunePMuonBestTrack, &b_Mu_absdzTunePMuonBestTrack);
   fChain->SetBranchAddress("Mu_ptInnerTrack", &Mu_ptInnerTrack, &b_Mu_ptInnerTrack);
   fChain->SetBranchAddress("Mu_pxInnerTrack", &Mu_pxInnerTrack, &b_Mu_pxInnerTrack);
   fChain->SetBranchAddress("Mu_pyInnerTrack", &Mu_pyInnerTrack, &b_Mu_pyInnerTrack);
   fChain->SetBranchAddress("Mu_pzInnerTrack", &Mu_pzInnerTrack, &b_Mu_pzInnerTrack);
   fChain->SetBranchAddress("Mu_pInnerTrack", &Mu_pInnerTrack, &b_Mu_pInnerTrack);
   fChain->SetBranchAddress("Mu_etaInnerTrack", &Mu_etaInnerTrack, &b_Mu_etaInnerTrack);
   fChain->SetBranchAddress("Mu_phiInnerTrack", &Mu_phiInnerTrack, &b_Mu_phiInnerTrack);
   fChain->SetBranchAddress("Mu_thetaInnerTrack", &Mu_thetaInnerTrack, &b_Mu_thetaInnerTrack);
   fChain->SetBranchAddress("Mu_chargeInnerTrack", &Mu_chargeInnerTrack, &b_Mu_chargeInnerTrack);
   fChain->SetBranchAddress("Mu_dPToverPTInnerTrack", &Mu_dPToverPTInnerTrack, &b_Mu_dPToverPTInnerTrack);
   fChain->SetBranchAddress("Mu_absdxyInnerTrack", &Mu_absdxyInnerTrack, &b_Mu_absdxyInnerTrack);
   fChain->SetBranchAddress("Mu_absdzInnerTrack", &Mu_absdzInnerTrack, &b_Mu_absdzInnerTrack);
   fChain->SetBranchAddress("Mu_normalizedChi2", &Mu_normalizedChi2, &b_Mu_normalizedChi2);
   fChain->SetBranchAddress("Mu_absdxy", &Mu_absdxy, &b_Mu_absdxy);
   fChain->SetBranchAddress("Mu_absdz", &Mu_absdz, &b_Mu_absdz);
   fChain->SetBranchAddress("Mu_vtxMass", &Mu_vtxMass, &b_Mu_vtxMass);
   fChain->SetBranchAddress("Mu_vtxNormChi2", &Mu_vtxNormChi2, &b_Mu_vtxNormChi2);
   fChain->SetBranchAddress("Mu_vtxMassLept", &Mu_vtxMassLept, &b_Mu_vtxMassLept);
   fChain->SetBranchAddress("Mu_numberOfMatchedStations", &Mu_numberOfMatchedStations, &b_Mu_numberOfMatchedStations);
   fChain->SetBranchAddress("Mu_numberOfValidPixelHits", &Mu_numberOfValidPixelHits, &b_Mu_numberOfValidPixelHits);
   fChain->SetBranchAddress("Mu_numberOfValidMuonHits", &Mu_numberOfValidMuonHits, &b_Mu_numberOfValidMuonHits);
   fChain->SetBranchAddress("Mu_numberOftrackerLayersWithMeasurement", &Mu_numberOftrackerLayersWithMeasurement, &b_Mu_numberOftrackerLayersWithMeasurement);
   fChain->SetBranchAddress("Mu_emIso", &Mu_emIso, &b_Mu_emIso);
   fChain->SetBranchAddress("Mu_hadIso", &Mu_hadIso, &b_Mu_hadIso);
   fChain->SetBranchAddress("Mu_trackiso", &Mu_trackiso, &b_Mu_trackiso);
   fChain->SetBranchAddress("Mu_pfSumChargedHadronPt", &Mu_pfSumChargedHadronPt, &b_Mu_pfSumChargedHadronPt);
   fChain->SetBranchAddress("Mu_pfSumNeutralHadronEt", &Mu_pfSumNeutralHadronEt, &b_Mu_pfSumNeutralHadronEt);
   fChain->SetBranchAddress("Mu_PFSumPhotonEt", &Mu_PFSumPhotonEt, &b_Mu_PFSumPhotonEt);
   fChain->SetBranchAddress("Mu_pfSumPUPt", &Mu_pfSumPUPt, &b_Mu_pfSumPUPt);
   fChain->SetBranchAddress("Mu_nbofpv", &Mu_nbofpv, &b_Mu_nbofpv);
   fChain->SetBranchAddress("Mu_patDeltaBeta", &Mu_patDeltaBeta, &b_Mu_patDeltaBeta);
   fChain->SetBranchAddress("HLT_nb", &HLT_nb, &b_HLT_nb);
   fChain->SetBranchAddress("HLT_name", &HLT_name, &b_HLT_name);
   fChain->SetBranchAddress("HLT_isaccept", &HLT_isaccept, &b_HLT_isaccept);
   fChain->SetBranchAddress("HLTObj_nbObj", &HLTObj_nbObj, &b_HLTObj_nbObj);
   fChain->SetBranchAddress("HLTObj_pt", &HLTObj_pt, &b_HLTObj_pt);
   fChain->SetBranchAddress("HLTObj_eta", &HLTObj_eta, &b_HLTObj_eta);
   fChain->SetBranchAddress("HLTObj_phi", &HLTObj_phi, &b_HLTObj_phi);
   fChain->SetBranchAddress("HLTObj_collection", &HLTObj_collection, &b_HLTObj_collection);
   fChain->SetBranchAddress("GenMet_pt", &GenMet_pt, &b_GenMet_pt);
   fChain->SetBranchAddress("PFMet_et_cor", &PFMet_et_cor, &b_PFMet_et_cor);
   fChain->SetBranchAddress("PFMet_pt_cor", &PFMet_pt_cor, &b_PFMet_pt_cor);
   fChain->SetBranchAddress("PFMet_phi_cor", &PFMet_phi_cor, &b_PFMet_phi_cor);
   fChain->SetBranchAddress("PFMet_en_cor", &PFMet_en_cor, &b_PFMet_en_cor);
   fChain->SetBranchAddress("PFMet_px_cor", &PFMet_px_cor, &b_PFMet_px_cor);
   fChain->SetBranchAddress("PFMet_py_cor", &PFMet_py_cor, &b_PFMet_py_cor);
   fChain->SetBranchAddress("PFMet_pz_cor", &PFMet_pz_cor, &b_PFMet_pz_cor);
   fChain->SetBranchAddress("PFMet_sumEt_cor", &PFMet_sumEt_cor, &b_PFMet_sumEt_cor);
   fChain->SetBranchAddress("CaloMet_pt", &CaloMet_pt, &b_CaloMet_pt);
   fChain->SetBranchAddress("CaloMet_phi", &CaloMet_phi, &b_CaloMet_phi);
   fChain->SetBranchAddress("CaloMet_sumEt", &CaloMet_sumEt, &b_CaloMet_sumEt);
   fChain->SetBranchAddress("PFMet_shiftedPt_JetEnUp", &PFMet_shiftedPt_JetEnUp, &b_PFMet_shiftedPt_JetEnUp);
   fChain->SetBranchAddress("PFMet_shiftedPt_JetEnDown", &PFMet_shiftedPt_JetEnDown, &b_PFMet_shiftedPt_JetEnDown);
   fChain->SetBranchAddress("jet_nb", &jet_nb, &b_jet_nb);
   fChain->SetBranchAddress("jet_charge", &jet_charge, &b_jet_charge);
   fChain->SetBranchAddress("jet_et", &jet_et, &b_jet_et);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_en", &jet_en, &b_jet_en);
   fChain->SetBranchAddress("jet_theta", &jet_theta, &b_jet_theta);
   fChain->SetBranchAddress("jet_beta", &jet_beta, &b_jet_beta);
   fChain->SetBranchAddress("jet_pileup_mva_disc", &jet_pileup_mva_disc, &b_jet_pileup_mva_disc);
   fChain->SetBranchAddress("Nb_bDiscriminators", &Nb_bDiscriminators, &b_Nb_bDiscriminators);
   fChain->SetBranchAddress("jet_btag_pt", &jet_btag_pt, &b_jet_btag_pt);
   fChain->SetBranchAddress("jet_btag_eta", &jet_btag_eta, &b_jet_btag_eta);
   fChain->SetBranchAddress("jet_btag_phi", &jet_btag_phi, &b_jet_btag_phi);
   fChain->SetBranchAddress("jet_btag_flavor", &jet_btag_flavor, &b_jet_btag_flavor);
   fChain->SetBranchAddress("jet_btag_pfCSVv2IVF_discriminator", &jet_btag_pfCSVv2IVF_discriminator, &b_jet_btag_pfCSVv2IVF_discriminator);
   fChain->SetBranchAddress("Nb_taus", &Nb_taus, &b_Nb_taus);
   fChain->SetBranchAddress("Tau_pt", &Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_eta", &Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_phi", &Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_id", &Tau_id, &b_Tau_id);
   fChain->SetBranchAddress("Tau_LooseCombinedIsolationDeltaBetaCorr3Hits", &Tau_LooseCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_LooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("pfphoton_size", &pfphoton_size, &b_pfphoton_size);
   fChain->SetBranchAddress("pfphoton_pt", &pfphoton_pt, &b_pfphoton_pt);
   fChain->SetBranchAddress("pfphoton_eta", &pfphoton_eta, &b_pfphoton_eta);
   fChain->SetBranchAddress("pfphoton_phi", &pfphoton_phi, &b_pfphoton_phi);
   fChain->SetBranchAddress("pfphoton_theta", &pfphoton_theta, &b_pfphoton_theta);
   fChain->SetBranchAddress("num_PU_vertices", &num_PU_vertices, &b_num_PU_vertices);
   fChain->SetBranchAddress("PU_BunchCrossing", &PU_BunchCrossing, &b_PU_BunchCrossing);
   fChain->SetBranchAddress("num_PU_gen_vertices", &num_PU_gen_vertices, &b_num_PU_gen_vertices);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("MC_weighting", &MC_weighting, &b_MC_weighting);
   fChain->SetBranchAddress("Ele_nbElectrons", &Ele_nbElectrons, &b_Ele_nbElectrons);
   fChain->SetBranchAddress("Ele_isEcalDrivenSeed", &Ele_isEcalDrivenSeed, &b_Ele_isEcalDrivenSeed);
   fChain->SetBranchAddress("Ele_isPassConversionVeto", &Ele_isPassConversionVeto, &b_Ele_isPassConversionVeto);
   fChain->SetBranchAddress("Ele_charge", &Ele_charge, &b_Ele_charge);
   fChain->SetBranchAddress("Ele_nbOfMissingHits", &Ele_nbOfMissingHits, &b_Ele_nbOfMissingHits);
   fChain->SetBranchAddress("Ele_nbVtx", &Ele_nbVtx, &b_Ele_nbVtx);
   fChain->SetBranchAddress("Ele_Et", &Ele_Et, &b_Ele_Et);
   fChain->SetBranchAddress("Ele_EtFromCaloEn", &Ele_EtFromCaloEn, &b_Ele_EtFromCaloEn);
   fChain->SetBranchAddress("Ele_pt", &Ele_pt, &b_Ele_pt);
   fChain->SetBranchAddress("Ele_thetaSC", &Ele_thetaSC, &b_Ele_thetaSC);
   fChain->SetBranchAddress("Ele_etaSC", &Ele_etaSC, &b_Ele_etaSC);
   fChain->SetBranchAddress("Ele_phiSC", &Ele_phiSC, &b_Ele_phiSC);
   fChain->SetBranchAddress("Ele_energySC", &Ele_energySC, &b_Ele_energySC);
   fChain->SetBranchAddress("Ele_preshowerEnergySC", &Ele_preshowerEnergySC, &b_Ele_preshowerEnergySC);
   fChain->SetBranchAddress("Ele_thetaTrack", &Ele_thetaTrack, &b_Ele_thetaTrack);
   fChain->SetBranchAddress("Ele_etaTrack", &Ele_etaTrack, &b_Ele_etaTrack);
   fChain->SetBranchAddress("Ele_phiTrack", &Ele_phiTrack, &b_Ele_phiTrack);
   fChain->SetBranchAddress("Ele_hadronicOverEm", &Ele_hadronicOverEm, &b_Ele_hadronicOverEm);
   fChain->SetBranchAddress("Ele_deltaEtaInSeedCluster", &Ele_deltaEtaInSeedCluster, &b_Ele_deltaEtaInSeedCluster);
   fChain->SetBranchAddress("Ele_deltaPhiInSeedCluster", &Ele_deltaPhiInSeedCluster, &b_Ele_deltaPhiInSeedCluster);
   fChain->SetBranchAddress("Ele_deltaEtaInSC", &Ele_deltaEtaInSC, &b_Ele_deltaEtaInSC);
   fChain->SetBranchAddress("Ele_deltaPhiInSC", &Ele_deltaPhiInSC, &b_Ele_deltaPhiInSC);
   fChain->SetBranchAddress("Ele_sigmaIetaIeta", &Ele_sigmaIetaIeta, &b_Ele_sigmaIetaIeta);
   fChain->SetBranchAddress("Ele_e2x5Max", &Ele_e2x5Max, &b_Ele_e2x5Max);
   fChain->SetBranchAddress("Ele_e1x5", &Ele_e1x5, &b_Ele_e1x5);
   fChain->SetBranchAddress("Ele_frac15", &Ele_frac15, &b_Ele_frac15);
   fChain->SetBranchAddress("Ele_frac51", &Ele_frac51, &b_Ele_frac51);
   fChain->SetBranchAddress("Ele_e5x5", &Ele_e5x5, &b_Ele_e5x5);
   fChain->SetBranchAddress("Ele3x3", &Ele3x3, &b_Ele3x3);
   fChain->SetBranchAddress("Ele_e2x5MaxOver5x5", &Ele_e2x5MaxOver5x5, &b_Ele_e2x5MaxOver5x5);
   fChain->SetBranchAddress("Ele_e1x5Over5x5", &Ele_e1x5Over5x5, &b_Ele_e1x5Over5x5);
   fChain->SetBranchAddress("Ele_sigmaIetaIetaFull5x5", &Ele_sigmaIetaIetaFull5x5, &b_Ele_sigmaIetaIetaFull5x5);
   fChain->SetBranchAddress("Ele_e2x5MaxFull5x5", &Ele_e2x5MaxFull5x5, &b_Ele_e2x5MaxFull5x5);
   fChain->SetBranchAddress("Ele_e1x5Full5x5", &Ele_e1x5Full5x5, &b_Ele_e1x5Full5x5);
   fChain->SetBranchAddress("Ele_e5x5Full5x5", &Ele_e5x5Full5x5, &b_Ele_e5x5Full5x5);
   fChain->SetBranchAddress("Ele_e2x5MaxOver5x5Full5x5", &Ele_e2x5MaxOver5x5Full5x5, &b_Ele_e2x5MaxOver5x5Full5x5);
   fChain->SetBranchAddress("Ele_e1x5Over5x5Full5x5", &Ele_e1x5Over5x5Full5x5, &b_Ele_e1x5Over5x5Full5x5);
   fChain->SetBranchAddress("Ele_e2x5Right", &Ele_e2x5Right, &b_Ele_e2x5Right);
   fChain->SetBranchAddress("Ele_e2x5Left", &Ele_e2x5Left, &b_Ele_e2x5Left);
   fChain->SetBranchAddress("Ele_e2x5Top", &Ele_e2x5Top, &b_Ele_e2x5Top);
   fChain->SetBranchAddress("Ele_e2x5Bottom", &Ele_e2x5Bottom, &b_Ele_e2x5Bottom);
   fChain->SetBranchAddress("Ele_eMax", &Ele_eMax, &b_Ele_eMax);
   fChain->SetBranchAddress("Ele_eRight", &Ele_eRight, &b_Ele_eRight);
   fChain->SetBranchAddress("Ele_eLeft", &Ele_eLeft, &b_Ele_eLeft);
   fChain->SetBranchAddress("Ele_eTop", &Ele_eTop, &b_Ele_eTop);
   fChain->SetBranchAddress("Ele_eBottom", &Ele_eBottom, &b_Ele_eBottom);
   fChain->SetBranchAddress("Ele_dxy", &Ele_dxy, &b_Ele_dxy);
   fChain->SetBranchAddress("Ele_dz", &Ele_dz, &b_Ele_dz);
   fChain->SetBranchAddress("Ele_rhoIso", &Ele_rhoIso, &b_Ele_rhoIso);
   fChain->SetBranchAddress("Ele_fbrem", &Ele_fbrem, &b_Ele_fbrem);
   fChain->SetBranchAddress("Ele_EoverP", &Ele_EoverP, &b_Ele_EoverP);
   fChain->SetBranchAddress("Ele_Xposition", &Ele_Xposition, &b_Ele_Xposition);
   fChain->SetBranchAddress("Ele_Yposition", &Ele_Yposition, &b_Ele_Yposition);
   fChain->SetBranchAddress("Ele_EcalPlusHcald1iso", &Ele_EcalPlusHcald1iso, &b_Ele_EcalPlusHcald1iso);
   fChain->SetBranchAddress("Ele_dr03EcalRecHitSumEt", &Ele_dr03EcalRecHitSumEt, &b_Ele_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Ele_dr03HcalDepth1TowerSumEt", &Ele_dr03HcalDepth1TowerSumEt, &b_Ele_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Ele_dr03HcalDepth1TowerSumEtBc", &Ele_dr03HcalDepth1TowerSumEtBc, &b_Ele_dr03HcalDepth1TowerSumEtBc);
   fChain->SetBranchAddress("Ele_hcalDepth1OverEcal", &Ele_hcalDepth1OverEcal, &b_Ele_hcalDepth1OverEcal);
   fChain->SetBranchAddress("Ele_hcalDepth2OverEcal", &Ele_hcalDepth2OverEcal, &b_Ele_hcalDepth2OverEcal);
   fChain->SetBranchAddress("Ele_dr03HcalDepth2TowerSumEt", &Ele_dr03HcalDepth2TowerSumEt, &b_Ele_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("Ele_hcalDepth2TowerSumEtNoVeto", &Ele_hcalDepth2TowerSumEtNoVeto, &b_Ele_hcalDepth2TowerSumEtNoVeto);
   fChain->SetBranchAddress("Ele_hcalDepth1TowerSumEtNoVeto", &Ele_hcalDepth1TowerSumEtNoVeto, &b_Ele_hcalDepth1TowerSumEtNoVeto);
   fChain->SetBranchAddress("Ele_pfSumPhotonEt", &Ele_pfSumPhotonEt, &b_Ele_pfSumPhotonEt);
   fChain->SetBranchAddress("Ele_pfSumChargedHadronPt", &Ele_pfSumChargedHadronPt, &b_Ele_pfSumChargedHadronPt);
   fChain->SetBranchAddress("Ele_pfSumNeutralHadronEt", &Ele_pfSumNeutralHadronEt, &b_Ele_pfSumNeutralHadronEt);
   fChain->SetBranchAddress("Ele_pfSumPUPt", &Ele_pfSumPUPt, &b_Ele_pfSumPUPt);
   fChain->SetBranchAddress("Ele_pfDeltaBeta", &Ele_pfDeltaBeta, &b_Ele_pfDeltaBeta);
   fChain->SetBranchAddress("Ele_rawId", &Ele_rawId, &b_Ele_rawId);
   fChain->SetBranchAddress("Ele_x", &Ele_x, &b_Ele_x);
   fChain->SetBranchAddress("Ele_y", &Ele_y, &b_Ele_y);
   fChain->SetBranchAddress("Ele_z", &Ele_z, &b_Ele_z);
   fChain->SetBranchAddress("Ele_zTrackPositionAtVtx", &Ele_zTrackPositionAtVtx, &b_Ele_zTrackPositionAtVtx);
   fChain->SetBranchAddress("Ele_ieta", &Ele_ieta, &b_Ele_ieta);
   fChain->SetBranchAddress("Ele_phiWidth", &Ele_phiWidth, &b_Ele_phiWidth);
   fChain->SetBranchAddress("Ele_etaWidth", &Ele_etaWidth, &b_Ele_etaWidth);
   fChain->SetBranchAddress("Ele_dr03TkSumPt", &Ele_dr03TkSumPt, &b_Ele_dr03TkSumPt);
   fChain->SetBranchAddress("Ele_dr03TkSumPt_corrected", &Ele_dr03TkSumPt_corrected, &b_Ele_dr03TkSumPt_corrected);
   Notify();
}

Bool_t ZprimeEleElePatMiniAodNewData::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZprimeEleElePatMiniAodNewData::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ZprimeEleElePatMiniAodNewData::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ZprimeEleElePatMiniAodNewData_cxx
