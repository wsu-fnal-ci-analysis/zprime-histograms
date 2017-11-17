//==============================================================
//          Analysis code for Z' boson to Mu Mu analysis       =
//          In this code we select the high pt di-muons events =
//          To run over MINIAOD MC with fixed trigger          =
//                  Author:  Sherif Elgammal                   =
//                                                             =
//                       13/05/2017                            =
//==============================================================
#ifndef ZprimeMuMuPatMiniAodNewMC_h
#define ZprimeMuMuPatMiniAodNewMC_h

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

class ZprimeMuMuPatMiniAodNewMC {
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
  std::vector<int>     *Mu_nbMuon;
  std::vector<bool>    *Mu_passOldMatchedStationsCut;
  std::vector<bool>    *Mu_passNewMatchedStationsCut;
  std::vector<bool>    *Mu_isTightMuon;
  std::vector<bool>    *Mu_isLooseMuon;
  //   std::vector<bool>    *Mu_isMuonsCleaned;
  std::vector<bool>    *Mu_isGlobalMuon;
  std::vector<bool>    *Mu_isPF;
  std::vector<bool>    *Mu_isTrackerMuon;
  std::vector<float>   *Mu_et;
  std::vector<float>   *Mu_en;
  std::vector<float>   *Mu_pt;
  std::vector<float>   *Mu_eta;
  std::vector<float>   *Mu_phi;
  std::vector<float>   *Mu_charge;
  std::vector<float>   *Mu_ptTunePMuonBestTrack;
  std::vector<float>   *Mu_pxTunePMuonBestTrack;
  std::vector<float>   *Mu_pyTunePMuonBestTrack;
  std::vector<float>   *Mu_pzTunePMuonBestTrack;
  std::vector<float>   *Mu_pTunePMuonBestTrack;
  std::vector<float>   *Mu_etaTunePMuonBestTrack;
  std::vector<float>   *Mu_phiTunePMuonBestTrack;
  std::vector<float>   *Mu_thetaTunePMuonBestTrack;
  std::vector<float>   *Mu_chargeTunePMuonBestTrack;
  std::vector<float>   *Mu_dPToverPTTunePMuonBestTrack;
  std::vector<float>   *Mu_absdxyTunePMuonBestTrack;
  std::vector<float>   *Mu_absdzTunePMuonBestTrack;
  std::vector<float>   *Mu_ptInnerTrack;
  std::vector<float>   *Mu_pxInnerTrack;
  std::vector<float>   *Mu_pyInnerTrack;
  std::vector<float>   *Mu_pzInnerTrack;
  std::vector<float>   *Mu_pInnerTrack;
  std::vector<float>   *Mu_etaInnerTrack;
  std::vector<float>   *Mu_phiInnerTrack;
  std::vector<float>   *Mu_thetaInnerTrack;
  std::vector<float>   *Mu_chargeInnerTrack;
  std::vector<float>   *Mu_dPToverPTInnerTrack;
  std::vector<float>   *Mu_absdxyInnerTrack;
  std::vector<float>   *Mu_absdzInnerTrack;
  std::vector<float>   *Mu_normalizedChi2;
  std::vector<float>   *Mu_absdxy;
  std::vector<float>   *Mu_absdz;
  std::vector<float>   *Mu_vtxMass;
  std::vector<float>   *Mu_vtxNormChi2;
  std::vector<float>   *Mu_vtxMassLept;
  std::vector<int>     *Mu_numberOfMatchedStations;
  std::vector<int>     *Mu_numberOfValidPixelHits;
  std::vector<int>     *Mu_numberOfValidMuonHits;
  std::vector<int>     *Mu_numberOftrackerLayersWithMeasurement;
  std::vector<float>   *Mu_emIso;
  std::vector<float>   *Mu_hadIso;
  std::vector<float>   *Mu_trackiso;
  std::vector<float>   *Mu_pfSumChargedHadronPt;
  std::vector<float>   *Mu_pfSumNeutralHadronEt;
  std::vector<float>   *Mu_PFSumPhotonEt;
  std::vector<float>   *Mu_pfSumPUPt;
  std::vector<int>     *Mu_nbofpv;
  std::vector<float>   *Mu_patDeltaBeta;
  std::vector<int>     *Mu_stationMask;
  std::vector<int>     *Mu_numberOfMatchedRPCLayers;
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
  TBranch        *b_Mu_nbMuon;   //!
  TBranch        *b_Mu_passOldMatchedStationsCut;   //!
  TBranch        *b_Mu_passNewMatchedStationsCut;   //!
  TBranch        *b_Mu_isTightMuon;   //!
  TBranch        *b_Mu_isLooseMuon;   //!
  //   TBranch        *b_Mu_isMuonsCleaned;   //!
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
  TBranch        *b_Mu_stationMask;   //!
  TBranch        *b_Mu_numberOfMatchedRPCLayers;   //!
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

  ZprimeMuMuPatMiniAodNewMC(Char_t namechar_[300], TTree *tree=0, Double_t weight_=1.,
                            std::string DATA_type_="DATA", std::string MC_type_="MC");
  virtual ~ZprimeMuMuPatMiniAodNewMC();
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
  void initMemberVariables();
  bool DiPFJet(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2);
  bool DiPFJetCut(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2);
  void PrintEventInformation(unsigned int runNumber, unsigned int lumiNumber, unsigned int eventNumber,
                             float vtxChi2, float vtxMass, float CosmicRejection);
  bool SelectFirstMuon(float &pTmuon1,float &Enmuon1,float &Etamuon1,
                       float &Phimuon1,int &ChargeMu1,unsigned &FlagMu1,
                       float &pxmuon1,float &pymuon1,float &pzmuon1,
                       float &pmuon1,float &dxymuon1,float &pTmuon1tuneP,
                       float &pTmuonBestTrack1,
                       float &genMu1Pt, float &genMu1Eta, float &genMu1Phi, float &genMu1En);
  bool SelectSecondMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float Etamuon1,float Phimuon1,
                        float &pTmuon2,float &Enmuon2,
                        float &Etamuon2,float &Phimuon2,int &ChargeMu2,float &pxmuon2,
                        float &pymuon2,float &pzmuon2,float &pmuon2,float &dxymuon2,
                        float &pTmuon2tuneP,float &pTmuonBestTrack2,
                        float &genMu2Pt, float &genMu2Eta, float &genMu2Phi, float &genMu2En);
  float Mass(float Pt1,float Eta1,float Phi1,float En1,
             float Pt2,float Eta2,float Phi2,float En2);
  float smearedMass(float Eta1,float Eta2,float vtxMass,float genMass, float &scaleUnc);
  void PlotRecoInfo(float CosmicMuonRejec,float vertexMassMu,float MassGenerated,
                    float PtTunePMuBestTrack,float PtTunePMu,float PtMuBestTrack,
                    float PtGenerated,float etaMu1,
                    float PtTunePMuBestTrack2,float PtTunePMu2,float PtMuBestTrack2,
                    float PtGenerated2,float etaMu2, float bosonPt);
  void PickThehighestMass(float &vtxHighestMass,float &vtxHighestChi2,int EvtNb);
  double ThreeDangle(float pxMu1,float pyMu1,float pzMu1,float pMu1,
                     float pxMu2,float pyMu2,float pzMu2,float pMu2);
  bool SelectFirstGenMu(float &ETMu1,float &PhiMu1,
                        float &EtaMu1,float &EnMu1,
                        int &IDele1,int &Statele1,
                        unsigned &GenFlag1);
  bool SelectSecondGenMu(unsigned GenFlag1,float ETMu1,float &ETMu2,float &PhiMu2,
                         float &EtaMu2,float &EnMu2,int &IDele2,int &Statele2);
  float GenMass(float Pt1,float Eta1,float Phi1,float En1,
                float Pt2,float Eta2,float Phi2,float En2);
  /*
    bool GenRecoMatchMu1(float RecoEta1,float RecoPhi1,
    float &ETMu1,float &PhiMu1,
    float &EtaMu1,float &EnMu1,
    unsigned &GenFlag1);
  */
  bool GenRecoMatchMu(float RecoEta1,float RecoPhi1,
                      float &ptGmu,float &etaGmu,float &phiGmu,float &enGmu);
  /*
    bool GenRecoMatchMu2(unsigned GenFlag1,float RecoEta2,float RecoPhi2,
    float &ETMu2,float &PhiMu2,
    float &EtaMu2,float &EnMu2);
  */
  void PlotGenInfo(float ZprimeGenMass,float EtaGenMu1,float EtaGenMu2,float PtGenMu1,
                   float PtGenMu2,float EnGenMu1,float EnGenMu2);

  bool RecoHLTMuonMatching(float RecoEta,float RecoPhi);
  //bool RecoHLTMuonMatching(float hltEta,float hltPhi,float RecoEta,float RecoPhi);
  bool isPassHLT();
  void DrawBTaggingDiscriminator(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2);
  bool BTaggingDiscriminator(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2);
  bool BTaggingDiscriminator2(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2);
  bool BTaggingDiscriminator3(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2);
  void Boson(float Pt1,float Eta1,float Phi1,float En1,
             float Pt2,float Eta2,float Phi2,float En2,
             float ChargeEle1,float MetEt,float MetPx,
             float MetPy,float MetPz,float MetEn, float &bosonPt);
  void MuonPassingID();
  void PlotPterror();
  void PlotNbTrackLayers();
  void PlotNBValidPixelHits();
  void PlotNbValidMuonHits();
  void PlotNbMatchedStations();
  void PlotTrackiso();
  void PlotAbsDxy();
  void PlotPtTuneP();
  void plotAllHighPtMuonsID();
  void MuonPassingNewID();
  void MuonPassingTightID();
  float CosThetaCollinSoper(float Et1,float Eta1,float Phi1,float En1,
                            float Et2,float Eta2,float Phi2,float En2,
                            float ChargeEle1,float RecoMass);
  float delR(float eta1,float phi1,float eta2,float phi2);
  double MassCorrection(float M, float bosonPt, float eta1, float eta2);
  void DrawDiJetMassBB();
  void DrawDiJetMassBE();
  void DrawDiJetMassEE();
  void DrawWJetsMassBB();
  void DrawWJetsMassBE1();
  void DrawWJetsMassBE2();
  void DrawWJetsMassEE();
  float FRweight(float eta, float pt);

  //=========Do any of these ever get properly initialized?=========================
  float FR_Ptcut;
  float parEB1,parEB2,parEB3,parEB4,parEB5,parEB6,parEB7,parEB8,parEB9,parEB10;
  float parEE1,parEE2,parEE3,parEE4,parEE5,parEE6,parEE7,parEE8,parEE9,parEE10;
  float parEB11,parEB22,parEB33;
  float parEE11,parEE22,parEE33;
  float HLT_pt,HLT_eta,HLT_phi;
  //int m_weight;
  float ptEffCut;
  float PtDYTRecMu1,PtDYTRecMu2,PtRecTunePMu1,PtRecTunePMu2,PtRecMuBestTrack1,PtRecMuBestTrack2;
  float RecoHLTMatchingDeltaRcut,deltaRcut,minMassCut,maxMassCut;
  float m_vtxChi2Mu,m_vtxMassMu,m_vtxMassSmearedMu,m_scaleUnc,m_csAngle;
  float m_ptGen1,m_phiGen1,m_etaGen1,m_enGen1;
  unsigned m_genFlag1;
  float m_ptGen2,m_phiGen2,m_etaGen2,m_enGen2;
  int ChargeRecMu1,ChargeRecMu2;
  unsigned flagmu1;
  unsigned flag1;
  float PtRecTunePMuBestTrack1,EnRecMu1,EtaRecMu1,PhiRecMu1;
  float PtRecTunePMuBestTrack2,EnRecMu2,EtaRecMu2,PhiRecMu2;
  float pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,dxyRecMu1;
  float pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2;

  float m_genET1,m_genPhi1,m_genEta1,m_genEn1;
  int m_genID1,m_genStat1;
  float m_genET2,m_genPhi2,m_genEta2,m_genEn2;
  int m_genID2,m_genStat2;
  float m_genMass, m_recoMass, m_recoMassCorr, m_bosonPt;  // seems not used...
  int m_nbGen,m_nbReco;
  int nbTP,nbTT,nbTF;
  float TagProbeEtaCut;
  float Eff;
  float MassCutMin,MassCutMax;
  float MassResolution;
  float EtaCut;

  int m_nbFireHLT;
  std::ofstream output_txt;

  // HISTOGRAMS
  /* std::shared_ptr<TH1D> h1_ZprimeRecomassBeforeTrigger_; */
  std::shared_ptr<TH1D> h1_ZprimeRecomass_;
  std::shared_ptr<TH1D> h1_ZprimeRecomasslogscale_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass60to120_;
  /* std::shared_ptr<TH1D> h1_ZprimeRecomassAbove400GeV_; */
  /* std::shared_ptr<TH1D> h1_ZprimeRecomassAbove1000GeV_; */
  std::shared_ptr<TH1D> h1_3Dangle_;
  /* std::shared_ptr<TH1D> h1_DxyDiff_; */
  std::shared_ptr<TH1D> h1_ZprimeGenmass_;
  std::shared_ptr<TH1D> h1_ZprimeGenEta1_;
  std::shared_ptr<TH1D> h1_ZprimeGenEta2_;
  std::shared_ptr<TH1D> h1_ZprimeGenPt1_;
  std::shared_ptr<TH1D> h1_ZprimeGenPt2_;
  std::shared_ptr<TH1D> h1_ZprimeGenEn1_;
  std::shared_ptr<TH1D> h1_ZprimeGenEn2_;
  /* std::shared_ptr<TH1D> h1_MassRecoGenDif_; */
  std::shared_ptr<TH1D> h1_dPToverPT_;
  /* std::shared_ptr<TH1D> h1_normalizedChi2_; */
  std::shared_ptr<TH1D> h1_numberOftrackerLayersWithMeasurement_;
  std::shared_ptr<TH1D> h1_numberOfValidPixelHits_;
  std::shared_ptr<TH1D> h1_numberOfValidMuonHits_;
  std::shared_ptr<TH1D> h1_numberOfMatchedStations_;
  std::shared_ptr<TH1D> h1_trackiso_;
  std::shared_ptr<TH1D> h1_absdxy_;
  std::shared_ptr<TH1D> h1_PtEffpterror_;
  std::shared_ptr<TH1D> h1_PtEffptnumberOftrackerLayers_;
  std::shared_ptr<TH1D> h1_PtEffptnumberOfPixelHits_;
  std::shared_ptr<TH1D> h1_PtEffptnumberOfMuonHits_;
  std::shared_ptr<TH1D> h1_PtEffptnumberOfMatchedStations_;
  std::shared_ptr<TH1D> h1_PtEffptTrackIso_;
  std::shared_ptr<TH1D> h1_PtEffptabsdsy_;
  /* std::shared_ptr<TH1D> h1_PtEffpfSumChargedHadron_; */
  /* std::shared_ptr<TH1D> h1_PtEffpfSumNeutralHadron_; */
  /* std::shared_ptr<TH1D> h1_PtEffpfPhotonIso_; */
  /* std::shared_ptr<TH1D> h1_FracpfSumChargedHadron_; */
  /* std::shared_ptr<TH1D> h1_FracpfSumNeutralHadron_; */
  /* std::shared_ptr<TH1D> h1_FracpfPhotonIso_; */
  std::shared_ptr<TH1D> h1_EtaEffpterror_;
  std::shared_ptr<TH1D> h1_EtaEffptnumberOftrackerLayers_;
  std::shared_ptr<TH1D> h1_EtaEffptnumberOfPixelHits_;
  std::shared_ptr<TH1D> h1_EtaEffptnumberOfMuonHits_;
  std::shared_ptr<TH1D> h1_EtaEffptnumberOfMatchedStations_;
  std::shared_ptr<TH1D> h1_EtaEffptTrackIso_;
  std::shared_ptr<TH1D> h1_EtaEffptabsdsy_;
  std::shared_ptr<TH1D> h1_nbPVID_;
  std::shared_ptr<TH1D> h1_PtID_;
  std::shared_ptr<TH1D> h1_EtaID_;
  std::shared_ptr<TH1D> h1_nbPVNewID_;
  std::shared_ptr<TH1D> h1_PtNewID_;
  std::shared_ptr<TH1D> h1_EtaNewID_;
  std::shared_ptr<TH1D> h1_nbPVTightID_;
  std::shared_ptr<TH1D> h1_PtTightID_;
  std::shared_ptr<TH1D> h1_EtaTightID_;

  std::shared_ptr<TH1D> h1_PtResolutionMBT_;
  std::shared_ptr<TH1D> h1_PtResolutionTunePMBT_;
  std::shared_ptr<TH1D> h1_PtResolutiontuneP_;
  std::shared_ptr<TH1D> h1_MassResultionEBEB1_;
  std::shared_ptr<TH1D> h1_MassResultionEBEB2_;
  std::shared_ptr<TH1D> h1_MassResultionEBEB3_;
  std::shared_ptr<TH1D> h1_MassResultionEBEB4_;
  std::shared_ptr<TH1D> h1_MassResultionEBEB5_;
  std::shared_ptr<TH1D> h1_MassResultionEBEB6_;
  std::shared_ptr<TH1D> h1_MassResultionEBEB7_;
  /* std::shared_ptr<TH1D> h1_MassResultionEBEB8_; */
  /* std::shared_ptr<TH1D> h1_MassResultionEBEB9_; */
  /* std::shared_ptr<TH1D> h1_MassResultionEBEB10_; */
  std::shared_ptr<TH1D> h1_MassGenInAccep_;
  std::shared_ptr<TH1D> h1_MassRecoInAccep_;
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

  /* std::shared_ptr<TH1D> h1_ptHistoBefor_; */
  /* std::shared_ptr<TH1D> h1_ptHistoPassingVtxChi2Mu_; */
  /* std::shared_ptr<TH1D> h1_ptHistoPassingCosmicRejec_; */
  /* std::shared_ptr<TH1D> h1_ptHistoPassingHLT_; */
  /* std::shared_ptr<TH1D> h1_etaHistoBefor_; */
  /* std::shared_ptr<TH1D> h1_etaHistoPassingVtxChi2Mu_; */
  /* std::shared_ptr<TH1D> h1_etaHistoPassingCosmicRejec_; */
  /* std::shared_ptr<TH1D> h1_etaHistoPassingHLT_; */
  /* std::shared_ptr<TH1D> h1_3DangleHisto1_; */
  /* std::shared_ptr<TH1D> h1_3DangleHisto2_; */
  /* std::shared_ptr<TH1D> h1_Fail3DangleHistoMass_; */
  /* std::shared_ptr<TH1D> h1_Fail3DangleHistoPhi_; */
  /* std::shared_ptr<TH1D> h1_ptHistoFRDum_; */
  /* std::shared_ptr<TH1D> h1_ptHistoFRNum_; */
  std::shared_ptr<TH1D> h1_ZprimeRecomass50_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBB_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassEE_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBE_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass20_;
  /* std::shared_ptr<TH2F> h2_ZprimeRecomassNewbin_; */
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidth_;

  /* std::shared_ptr<TH1D> h1_NbFireHLT; */
  std::shared_ptr<TH1D> h1_ptAfterTrigger_;
  std::shared_ptr<TH1D> h1_ptBeforeTrigger_;
  std::shared_ptr<TH1D> h1_PtTuneP_;
  std::shared_ptr<TH1D> h1_PFMetCorr_;
  std::shared_ptr<TH1D> h1_CaloMet_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthBB_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthBEpos_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthBEnev_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthEE_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass60to120BEpos_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass60to120BEnev_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass60to120EE_;
  std::shared_ptr<TH1D> h1_ZprimeRecomass60to120BB_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthAfterBtaging_;
  /* std::shared_ptr<TH1D> h1_jetBTag_; */
  /* std::shared_ptr<TH1D> h1_jetBTagB_; */
  /* std::shared_ptr<TH1D> h1_jetBTagC_; */
  /* std::shared_ptr<TH1D> h1_jetBTagUDSG_; */
  std::shared_ptr<TH1D> h1_jetBTagStep1_;
  std::shared_ptr<TH1D> h1_jetBTagStep2_;
  std::shared_ptr<TH1D> h1_jetBTagStep3_;
  std::shared_ptr<TH1D> h1_nbBTagStep1_;
  std::shared_ptr<TH1D> h1_nbBTagStep2_;
  std::shared_ptr<TH1D> h1_nbBTagStep3_;
  std::shared_ptr<TH1D> h1_BosPt_;
  std::shared_ptr<TH1D> h1_BosPhi_;
  std::shared_ptr<TH1D> h1_DeltaPhi_;
  std::shared_ptr<TH1D> h1_DeltaPtoverPt_;
  std::shared_ptr<TH1D> h1_Mt_;
  std::shared_ptr<TH1D> h1_MissingEt_;
  std::shared_ptr<TH1D> h1_BTagMassMuMu_;
  std::shared_ptr<TH1D> h1_BTagMassMuMu1GeVbin_;
  std::shared_ptr<TH1D> h1_BTagMassMuMuBinWidth_;
  std::shared_ptr<TH1D> h1_MassMuMuDijetBinWidth_;
  std::shared_ptr<TH1D> h1_MassMuMuDijet1GeVbin_;
  std::shared_ptr<TH1D> h1_MassMuMuDijetBinWidthMET_;
  std::shared_ptr<TH1D> h1_MassMuMuBinWidthMET_;
  std::shared_ptr<TH1D> h1_MassMuMuDijet1GeVbinMET_;
  std::shared_ptr<TH1D> h1_MassMuMu1GeVbinMET_;
  std::shared_ptr<TH1D> h1_MassMuMuDijetBinWidthMET100_;
  std::shared_ptr<TH1D> h1_MassMuMuDijet1GeVbinMET100_;
  std::shared_ptr<TH1D> h1_NbPFjetsAll_;
  std::shared_ptr<TH1D> h1_NbPFjets2_;
  std::shared_ptr<TH1D> h1_ptPFjetsAll_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthAllBE_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthAllEE_;
  std::shared_ptr<TH1D> h1_ZprimeRecomassBinWidthAll_;
  std::shared_ptr<TH1D> h1_DijetBinWidthBB_;
  std::shared_ptr<TH1D> h1_DijetBinWidthBE_;
  std::shared_ptr<TH1D> h1_DijetBinWidthEE_;
  std::shared_ptr<TH1D> h1_DijetBinWidthBBBE_;
  std::shared_ptr<TH1D> h1_DijetEta1_;
  std::shared_ptr<TH1D> h1_DijetEta2_;
  std::shared_ptr<TH1D> h1_WjetsBinWidthBB_;
  std::shared_ptr<TH1D> h1_WjetsBinWidthBE_;
  std::shared_ptr<TH1D> h1_WjetsBinWidthEE_;
  std::shared_ptr<TH1D> h1_WjetsBinWidthBBBE_;
  std::shared_ptr<TH1D> h1_Dijet1GeVBB_;
  std::shared_ptr<TH1D> h1_Dijet1GeVBEEE_;
  std::shared_ptr<TH1D> h1_Dijet1GeVEE_;
  std::shared_ptr<TH1D> h1_Dijet1GeVBBBEEE_;
  std::shared_ptr<TH1D> h1_Wjets1GeVBB_;
  std::shared_ptr<TH1D> h1_Wjets1GeVBEEE_;
  std::shared_ptr<TH1D> h1_Wjets1GeVEE_;
  std::shared_ptr<TH1D> h1_Wjets1GeVBBBEEE_;
  std::shared_ptr<TH1D> h1_Dijet20GeVBB_;
  std::shared_ptr<TH1D> h1_Dijet20GeVBEEE_;
  std::shared_ptr<TH1D> h1_Dijet20GeVBBBEEE_;
  std::shared_ptr<TH1D> h1_Wjets20GeVBB_;
  std::shared_ptr<TH1D> h1_Wjets20GeVBEEE_;
  std::shared_ptr<TH1D> h1_Wjets20GeVBBBEEE_;
};

#endif

#ifdef ZprimeMuMuPatMiniAodNewMC_cxx
ZprimeMuMuPatMiniAodNewMC::ZprimeMuMuPatMiniAodNewMC(Char_t namechar_[300], TTree *tree, Double_t weight_,
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
      f = new TFile("/data/DATA/temp_elgammal/sherif807/CMSSW_8_0_7/src/2016MCAfterICHEP/Keep_Moriond17_reMINIAOD_GiovaniFilter_MC_muon/SingleMuon_MCB_tree.root");
    }
    f->GetObject("tree",tree);

  }
  Init(tree);
  initMemberVariables();
}

ZprimeMuMuPatMiniAodNewMC::~ZprimeMuMuPatMiniAodNewMC()
{
  if (!fChain) return;
  std::cout << std::hex << fChain << std::endl;
  std::cout << std::hex << fChain->GetCurrentFile() << std::endl;
  if (fChain->GetCurrentFile())
    delete fChain->GetCurrentFile();
}

Int_t ZprimeMuMuPatMiniAodNewMC::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t ZprimeMuMuPatMiniAodNewMC::LoadTree(Long64_t entry)
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

void ZprimeMuMuPatMiniAodNewMC::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  /* // special for CI smaples */
  /* xsWeight          = -1; */
  /* passPreFSRMInvCut = false; */
  /* passMInvCut       = false; */

  HLT_nb       = 0;
  HLT_name     = 0;
  HLT_isaccept = 0;
  HLTObj_nbObj = 0;
  HLTObj_pt    = 0;
  HLTObj_eta   = 0;
  HLTObj_phi   = 0;
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
  Mu_nbMuon = 0;
  Mu_passOldMatchedStationsCut = 0;
  Mu_passNewMatchedStationsCut = 0;
  Mu_isTightMuon = 0;
  Mu_isLooseMuon = 0;
  //   Mu_isMuonsCleaned = 0;
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
  Mu_stationMask = 0;
  Mu_numberOfMatchedRPCLayers = 0;
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
  fChain   = tree;
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

  fChain->SetBranchAddress("HLT_nb",       &HLT_nb,       &b_HLT_nb);
  fChain->SetBranchAddress("HLT_name",     &HLT_name,     &b_HLT_name);
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
  fChain->SetBranchAddress("Mu_nbMuon", &Mu_nbMuon, &b_Mu_nbMuon);
  fChain->SetBranchAddress("Mu_passOldMatchedStationsCut", &Mu_passOldMatchedStationsCut, &b_Mu_passOldMatchedStationsCut);
  fChain->SetBranchAddress("Mu_passNewMatchedStationsCut", &Mu_passNewMatchedStationsCut, &b_Mu_passNewMatchedStationsCut);
  fChain->SetBranchAddress("Mu_isTightMuon", &Mu_isTightMuon, &b_Mu_isTightMuon);
  fChain->SetBranchAddress("Mu_isLooseMuon", &Mu_isLooseMuon, &b_Mu_isLooseMuon);
  //   fChain->SetBranchAddress("Mu_isMuonsCleaned", &Mu_isMuonsCleaned, &b_Mu_isMuonsCleaned);
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
  fChain->SetBranchAddress("Mu_stationMask", &Mu_stationMask, &b_Mu_stationMask);
  fChain->SetBranchAddress("Mu_numberOfMatchedRPCLayers", &Mu_numberOfMatchedRPCLayers, &b_Mu_numberOfMatchedRPCLayers);
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

Bool_t ZprimeMuMuPatMiniAodNewMC::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void ZprimeMuMuPatMiniAodNewMC::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t ZprimeMuMuPatMiniAodNewMC::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef ZprimeMuMuPatMiniAodNewMC_cxx
