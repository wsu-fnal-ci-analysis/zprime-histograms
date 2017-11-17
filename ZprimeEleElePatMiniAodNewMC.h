//==============================================================
//          Analysis code for Z' boson to e+ e- analysis       =
//          In this code we select the HEEP events             =
//          To run over MINIAOD MC with fixed trigger          =
//                  Author:  Sherif Elgammal                   =
//                                                             =
//                       11/03/2016                            =
//==============================================================
#ifndef ZprimeEleElePatMiniAodNewMC_h
#define ZprimeEleElePatMiniAodNewMC_h

#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include <iostream>
#include <memory>
#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <sstream>
#include "TVector3.h"
// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class ZprimeEleElePatMiniAodNewMC {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  UInt_t          event_runNo;
  UInt_t          event_evtNo;
  UInt_t          event_lumi;
  UInt_t          event_bunch;

  // special for CI smaples
  bool            isCISample;
  Double_t        xsWeight;
  Bool_t          passPreFSRMInvCut;
  Bool_t          passMInvCut;

  std::vector<int>     *Ele_nbElectrons;
  std::vector<bool>    *Ele_isEcalDrivenSeed;
  std::vector<bool>    *Ele_isPassConversionVeto;
  std::vector<int>     *Ele_charge;
  std::vector<int>     *Ele_nbOfMissingHits;
  std::vector<int>     *Ele_nbVtx;
  std::vector<float>   *Ele_Et;
  std::vector<float>   *Ele_EtFromCaloEn;
  std::vector<float>   *Ele_pt;
  std::vector<float>   *Ele_thetaSC;
  std::vector<float>   *Ele_etaSC;
  std::vector<float>   *Ele_phiSC;
  std::vector<float>   *Ele_energySC;
  std::vector<float>   *Ele_preshowerEnergySC;
  std::vector<float>   *Ele_thetaTrack;
  std::vector<float>   *Ele_etaTrack;
  std::vector<float>   *Ele_phiTrack;
  std::vector<float>   *Ele_hadronicOverEm;
  std::vector<float>   *Ele_deltaEtaInSeedCluster;
  std::vector<float>   *Ele_deltaPhiInSeedCluster;
  std::vector<float>   *Ele_deltaEtaInSC;
  std::vector<float>   *Ele_deltaPhiInSC;
  std::vector<float>   *Ele_sigmaIetaIeta;
  std::vector<float>   *Ele_e2x5Max;
  std::vector<float>   *Ele_e1x5;
  std::vector<float>   *Ele_frac15;
  std::vector<float>   *Ele_frac51;
  std::vector<float>   *Ele_e5x5;
  std::vector<float>   *Ele3x3;
  std::vector<float>   *Ele_e2x5MaxOver5x5;
  std::vector<float>   *Ele_e1x5Over5x5;
  std::vector<float>   *Ele_sigmaIetaIetaFull5x5;
  std::vector<float>   *Ele_e2x5MaxFull5x5;
  std::vector<float>   *Ele_e1x5Full5x5;
  std::vector<float>   *Ele_e5x5Full5x5;
  std::vector<float>   *Ele_e2x5MaxOver5x5Full5x5;
  std::vector<float>   *Ele_e1x5Over5x5Full5x5;
  std::vector<float>   *Ele_e2x5Right;
  std::vector<float>   *Ele_e2x5Left;
  std::vector<float>   *Ele_e2x5Top;
  std::vector<float>   *Ele_e2x5Bottom;
  std::vector<float>   *Ele_eMax;
  std::vector<float>   *Ele_eRight;
  std::vector<float>   *Ele_eLeft;
  std::vector<float>   *Ele_eTop;
  std::vector<float>   *Ele_eBottom;
  std::vector<float>   *Ele_dxy;
  std::vector<float>   *Ele_dz;
  std::vector<float>   *Ele_rhoIso;
  std::vector<float>   *Ele_fbrem;
  std::vector<float>   *Ele_EoverP;
  std::vector<float>   *Ele_Xposition;
  std::vector<float>   *Ele_Yposition;
  std::vector<float>   *Ele_EcalPlusHcald1iso;
  std::vector<float>   *Ele_dr03EcalRecHitSumEt;
  std::vector<float>   *Ele_dr03HcalDepth1TowerSumEt;
  std::vector<float>   *Ele_dr03HcalDepth1TowerSumEtBc;
  std::vector<float>   *Ele_hcalDepth1OverEcal;
  std::vector<float>   *Ele_hcalDepth2OverEcal;
  std::vector<float>   *Ele_dr03HcalDepth2TowerSumEt;
  std::vector<float>   *Ele_hcalDepth2TowerSumEtNoVeto;
  std::vector<float>   *Ele_hcalDepth1TowerSumEtNoVeto;
  std::vector<float>   *Ele_pfSumPhotonEt;
  std::vector<float>   *Ele_pfSumChargedHadronPt;
  std::vector<float>   *Ele_pfSumNeutralHadronEt;
  std::vector<float>   *Ele_pfSumPUPt;
  std::vector<float>   *Ele_pfDeltaBeta;
  std::vector<int>     *Ele_rawId;
  std::vector<float>   *Ele_x;
  std::vector<float>   *Ele_y;
  std::vector<float>   *Ele_z;
  std::vector<float>   *Ele_zTrackPositionAtVtx;
  std::vector<int>     *Ele_ieta;
  std::vector<float>   *Ele_phiWidth;
  std::vector<float>   *Ele_etaWidth;
  std::vector<float>   *Ele_dr03TkSumPt;
  std::vector<float>   *Ele_dr03TkSumPt_corrected;
  std::vector<float>   *Ele_nrSatCrys;
  std::vector<bool>    *Ele_isPassHeepID;
  std::vector<int>     *HLT_nb;
  std::vector<std::string>  *HLT_name;
  std::vector<bool>    *HLT_isaccept;
  std::vector<int>     *HLTObj_nbObj;
  std::vector<float>   *HLTObj_pt;
  std::vector<float>   *HLTObj_eta;
  std::vector<float>   *HLTObj_phi;
  std::vector<std::string>  *HLTObj_collection;
  std::vector<int>     *iGenJet;
  std::vector<int>     *idGenJet;
  std::vector<int>     *statusGenJet;
  std::vector<int>     *chargeGenJet;
  std::vector<float>   *ptGenJet;
  std::vector<float>   *etaGenJet;
  std::vector<float>   *phiGenJet;
  std::vector<int>     *iGen;
  std::vector<int>     *idGen;
  std::vector<int>     *statusGen;
  std::vector<float>   *ptGen;
  std::vector<float>   *etaGen;
  std::vector<float>   *phiGen;
  std::vector<int>     *chargeGen;
  std::vector<float>   *EnergyGen;
  std::vector<float>   *pxGen;
  std::vector<float>   *pyGen;
  std::vector<float>   *pzGen;
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
  std::vector<int>     *jet_nb;
  std::vector<float>   *jet_charge;
  std::vector<float>   *jet_et;
  std::vector<float>   *jet_pt;
  std::vector<float>   *jet_eta;
  std::vector<float>   *jet_phi;
  std::vector<float>   *jet_en;
  std::vector<float>   *jet_theta;
  std::vector<float>   *jet_beta;
  std::vector<float>   *jet_pileup_mva_disc;
  std::vector<int>     *Nb_bDiscriminators;
  std::vector<float>   *jet_btag_pt;
  std::vector<float>   *jet_btag_eta;
  std::vector<float>   *jet_btag_phi;
  std::vector<int>     *jet_btag_flavor;
  std::vector<float>   *jet_btag_pfCSVv2IVF_discriminator;
  std::vector<int>     *Nb_taus;
  std::vector<float>   *Tau_pt;
  std::vector<float>   *Tau_eta;
  std::vector<float>   *Tau_phi;
  std::vector<int>     *Tau_id;
  std::vector<float>   *Tau_LooseCombinedIsolationDeltaBetaCorr3Hits;
  Int_t           pfphoton_size;
  std::vector<float>   *pfphoton_pt;
  std::vector<float>   *pfphoton_eta;
  std::vector<float>   *pfphoton_phi;
  std::vector<float>   *pfphoton_theta;
  Int_t           num_PU_vertices;
  Int_t           PU_BunchCrossing;
  Int_t           num_PU_gen_vertices;
  Float_t         Rho;
  std::vector<float>   *MC_weighting;

  // List of branches
  TBranch        *b_event_runNo;   //!
  TBranch        *b_event_evtNo;   //!
  TBranch        *b_event_lumi;   //!
  TBranch        *b_event_bunch;   //!

  // special for CI smaples
  TBranch        *b_xsWeight;
  TBranch        *b_passPreFSRMInvCut;
  TBranch        *b_passMInvCut;

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
  TBranch        *b_Ele_nrSatCrys;   //!
  TBranch        *b_Ele_isPassHeepID;   //!
  TBranch        *b_HLT_nb;   //!
  TBranch        *b_HLT_name;   //!
  TBranch        *b_HLT_isaccept;   //!
  TBranch        *b_HLTObj_nbObj;   //!
  TBranch        *b_HLTObj_pt;   //!
  TBranch        *b_HLTObj_eta;   //!
  TBranch        *b_HLTObj_phi;   //!
  TBranch        *b_HLTObj_collection;   //!
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

  ZprimeEleElePatMiniAodNewMC(Char_t namechar_[300], TTree *tree=0, Double_t weight_=1.,
                              std::string DATA_type_="DATA", std::string MC_type_="MC");
  virtual ~ZprimeEleElePatMiniAodNewMC();
  Double_t m_weight;
  Char_t name[300];
  std::string DATA_type,MC_type;
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(bool debug=false);
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  double MassCorrection(float M);
  bool SelectFirstEle(float &ETele1,float &Enele1,float &Etaele1,
                      float &Phiele1,int &ChargeEle1,float &EtaSCele1,
                      float &PhiSCele1,unsigned &FlagEle1,
                      float &genEle1Pt, float &genEle1Eta, float &genEle1Phi, float &genEle1En);
  bool SelectSecondEle(int ChargeEle1,unsigned FlagEle1,float ETele1,float Etaele1,float Phiele1,
                       float &ETele2,float &Enele2,float &Etaele2,float &Phiele2,int &ChargeEle2,
                       float &EtaSCele2,float &PhiSCele2,
                       float &genEle2Pt, float &genEle2Eta, float &genEle2Phi, float &genEle2En);
  float Mass(float Pt1,float Eta1,float Phi1,float En1,
             float Pt2,float Eta2,float Phi2,float En2);
  float smearedMass(float Eta1,float Eta2,float vtxMass,float genMass, float &scaleUnc);
  void PlotRecoInfo(float MassEle,float etaEle1,float etaEle2);
  bool SelectFirstGenEle(float &ETEle1,float &PhiSCEle1,
                         float &EtaSCEle1,float &EnEle1,
                         int &IDele1,int &Statele1,
                         unsigned &GenFlag1);
  bool SelectSecondGenEle(unsigned GenFlag1,float ETEle1,float &ETEle2,float &PhiSCEle2,
                          float &EtaSCEle2,float &EnEle2,int &IDele2,int &Statele2);
  float GenMass(float Pt1,float Eta1,float Phi1,float En1,
                float Pt2,float Eta2,float Phi2,float En2);
  bool GenRecoMatchEle(float RecoEta1,float RecoPhi1,
                       float &ptGele,float &etaGele,float &phiGele,float &enGele);
  bool RecoHLTEleMatching(float RecoEta,float RecoPhi);
  bool isPassHLT();
  float CosThetaCollinSoper(float Et1,float Eta1,float Phi1,float En1,
			    float Et2,float Eta2,float Phi2,float En2,
			    float ChargeEle1,float RecoMass);
  float delR(float eta1,float phi1,float eta2,float phi2);

  //================================================================================
  float HLT_pt,HLT_eta,HLT_phi;
  /* int m_weight; BUG BUG BUG*/
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
  unsigned flagel1;
  unsigned flag1;
  float PtRecTunePMuBestTrack1,EnRecMu1,EtaRecMu1,PhiRecMu1;
  float PtRecTunePMuBestTrack2,EnRecMu2,EtaRecMu2,PhiRecMu2;
  float pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,dxyRecMu1;
  float pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2;
  float m_genET1,m_genPhi1,m_genEta1,m_genEn1;
  int m_genID1,m_genStat1;
  float m_genET2,m_genPhi2,m_genEta2,m_genEn2;
  int m_genID2,m_genStat2;
  float m_genMass,m_recoMass,m_recoMassSmeared,m_recoMassCorr,m_scaleUnc,m_csAngle;
  int NbGen,NbReco;
  int nbTP,nbTT,nbTF;
  float TagProbeEtaCut;
  float Eff;
  float MassCutMin,MassCutMax;
  float MassResolution;
  float EtaCut;
  int NbFireHLT;
  std::ofstream output_txt;

  std::shared_ptr<TH1D> h1_CosAngleCollinSoperCorrect60Mass120_;
  std::shared_ptr<TH1D> h1_CosAngleCollinSoperCorrect120Mass300_;
  std::shared_ptr<TH1D> h1_CosAngleCollinSoperCorrect300Mass700_;
  std::shared_ptr<TH1D> h1_CosAngleCollinSoperCorrect700Mass3000_;
  std::shared_ptr<TH1D> h1_CosAngleCollinSoperCorrect4900Mass5100_;
  std::shared_ptr<TH1D> h1_absCosAngleCollinSoperCorrect4500Mass5500_;

  /* std::array<std::array<TH1D*,4> 3> h1_SmearedMassBinned_; */
  /* std::array<std::array<TH1D*,4> 3> h1_MassBinned_       ; */
  /* std::array<std::array<TH1D*,4> 3> h1_MassUpBinned_     ; */
  /* std::array<std::array<TH1D*,4> 3> h1_MassDownBinned_   ; */
  std::shared_ptr<TH2D> h2_CSSmearedMassBinned_;
  std::shared_ptr<TH2D> h2_CSMassBinned_       ;
  std::shared_ptr<TH2D> h2_CSMassUpBinned_     ;
  std::shared_ptr<TH2D> h2_CSMassDownBinned_   ;
  std::shared_ptr<TH2D> h2_CSPosSmearedMassBinned_;
  std::shared_ptr<TH2D> h2_CSPosMassBinned_       ;
  std::shared_ptr<TH2D> h2_CSPosMassUpBinned_     ;
  std::shared_ptr<TH2D> h2_CSPosMassDownBinned_   ;
  std::shared_ptr<TH2D> h2_CSNegSmearedMassBinned_;
  std::shared_ptr<TH2D> h2_CSNegMassBinned_       ;
  std::shared_ptr<TH2D> h2_CSNegMassUpBinned_     ;
  std::shared_ptr<TH2D> h2_CSNegMassDownBinned_   ;

  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidth_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthBB_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthBE_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthEE_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass60to120BE_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass60to120EE_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass60to120BB_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass60to120_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass_;
};

#endif

#ifdef ZprimeEleElePatMiniAodNewMC_cxx
ZprimeEleElePatMiniAodNewMC::ZprimeEleElePatMiniAodNewMC(Char_t namechar_[300], TTree *tree, Double_t weight_,
                                                         std::string DATA_type_, std::string MC_type_) : fChain(0)
{
  sprintf(name,"%s",namechar_);
  m_weight = weight_;
  std::cout << "Name is= " << name << " and weight is=" << m_weight << std::endl;
  DATA_type = DATA_type_;
  MC_type = MC_type_;

  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("");
    if (!f || !f->IsOpen()) {
      f = new TFile("../ZToEE_NNPDF30_13TeV-powheg_M_3500_4500_Moriond17_tree.root");
    }
    f->GetObject("tree",tree);

  }
  Init(tree);
}

ZprimeEleElePatMiniAodNewMC::~ZprimeEleElePatMiniAodNewMC()
{
  if (!fChain) return;
  std::cout << std::hex << fChain << std::endl;
  std::cout << std::hex << fChain->GetCurrentFile() << std::endl;
  if (fChain->GetCurrentFile())
    delete fChain->GetCurrentFile();
}

Int_t ZprimeEleElePatMiniAodNewMC::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t ZprimeEleElePatMiniAodNewMC::LoadTree(Long64_t entry)
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

void ZprimeEleElePatMiniAodNewMC::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
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
  Ele_nrSatCrys = 0;
  Ele_isPassHeepID = 0;
  HLT_nb = 0;
  HLT_name = 0;
  HLT_isaccept = 0;
  HLTObj_nbObj = 0;
  HLTObj_pt = 0;
  HLTObj_eta = 0;
  HLTObj_phi = 0;
  HLTObj_collection = 0;
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
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("event_runNo", &event_runNo, &b_event_runNo);
  fChain->SetBranchAddress("event_evtNo", &event_evtNo, &b_event_evtNo);
  fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
  fChain->SetBranchAddress("event_bunch", &event_bunch, &b_event_bunch);

  // special for CI smaples (will fail if not present... need to robustify
  isCISample = false;
  if (std::string(name).find("CI") != std::string::npos) {
    isCISample = true;
    fChain->SetBranchAddress("xsWeight",          &xsWeight,          &b_xsWeight);
    fChain->SetBranchAddress("passPreFSRMInvCut", &passPreFSRMInvCut, &b_passPreFSRMInvCut);
    fChain->SetBranchAddress("passMInvCut",       &passMInvCut,       &b_passMInvCut);
  }

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
  fChain->SetBranchAddress("Ele_nrSatCrys", &Ele_nrSatCrys, &b_Ele_nrSatCrys);
  fChain->SetBranchAddress("Ele_isPassHeepID", &Ele_isPassHeepID, &b_Ele_isPassHeepID);
  fChain->SetBranchAddress("HLT_nb", &HLT_nb, &b_HLT_nb);
  fChain->SetBranchAddress("HLT_name", &HLT_name, &b_HLT_name);
  fChain->SetBranchAddress("HLT_isaccept", &HLT_isaccept, &b_HLT_isaccept);
  fChain->SetBranchAddress("HLTObj_nbObj", &HLTObj_nbObj, &b_HLTObj_nbObj);
  fChain->SetBranchAddress("HLTObj_pt", &HLTObj_pt, &b_HLTObj_pt);
  fChain->SetBranchAddress("HLTObj_eta", &HLTObj_eta, &b_HLTObj_eta);
  fChain->SetBranchAddress("HLTObj_phi", &HLTObj_phi, &b_HLTObj_phi);
  fChain->SetBranchAddress("HLTObj_collection", &HLTObj_collection, &b_HLTObj_collection);
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
  Notify();
}

Bool_t ZprimeEleElePatMiniAodNewMC::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void ZprimeEleElePatMiniAodNewMC::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t ZprimeEleElePatMiniAodNewMC::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef ZprimeEleElePatMiniAodNewMC_cxx
