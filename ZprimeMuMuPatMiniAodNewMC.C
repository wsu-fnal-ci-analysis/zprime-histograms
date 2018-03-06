//==============================================================
//          Analysis code for Z' boson to Mu Mu analysis       =
//          In this code we select the high pt di-muons events =
//          To run over MINIAOD MC with fixed trigger          =
//                  Author:  Sherif Elgammal                   =
//                                                             =
//                       13/05/2017                            =
//==============================================================
#define ZprimeMuMuPatMiniAodNewMC_cxx
#include "ZprimeMuMuPatMiniAodNewMC.h"
#include <math.h>
#include <array>
#include <iomanip>
#include "TStopwatch.h"
#include "TRandom3.h"
#include <time.h>
#include <algorithm>

// using namespace std;
#define PI 3.14159265358979
#define MUON_MASS 0.1056583
#define ELEC_MASS 0.000511

bool myfunction (int i,int j) { return (i<j); }
bool picklargemass (float lhs,float rhs) { return (lhs > rhs); }
TString inputfile;
float newweight = 1.;
float pu_weight=1.;

void ZprimeMuMuPatMiniAodNewMC::initMemberVariables()
{
  rand = std::make_shared<TRandom3>();

  m_nbGen        = 0;
  m_nbReco       = 0;
  MassCutMin     = 0.0;
  MassCutMax     = 2000.0;
  EtaCut         = 2.4;
  MassResolution = 0.10;
  deltaRcut      = 0.15;
  RecoHLTMatchingDeltaRcut = 0.20;
  minMassCut  = 50.0;
  maxMassCut  = 4000.0;
  ptEffCut    = 3000.0;
  FR_Ptcut    = 53.0; //53.0;

  parEB1 = -99999.;
  parEB2 = -99999.;
  parEB3 = -99999.;
  parEB4 = -99999.;
  parEB5 = -99999.;
  parEB6 = -99999.;
  parEB7 = -99999.;
  parEB8 = -99999.;
  parEB9 = -99999.;
  parEB10 = -99999.;

  parEE1 = -99999.;
  parEE2 = -99999.;
  parEE3 = -99999.;
  parEE4 = -99999.;
  parEE5 = -99999.;
  parEE6 = -99999.;
  parEE7 = -99999.;
  parEE8 = -99999.;
  parEE9 = -99999.;
  parEE10 = -99999.;

  parEB11 = -99999.;
  parEB22 = -99999.;
  parEB33 = -99999.;
  parEE11 = -99999.;
  parEE22 = -99999.;
  parEE33 = -99999.;

  HLT_pt  = -99999.;
  HLT_eta = -99999.;
  HLT_phi = -99999.;

  PtDYTRecMu1       = -99999.;
  PtDYTRecMu2       = -99999.;
  PtRecTunePMu1     = -99999.;
  PtRecTunePMu2     = -99999.;
  PtRecMuBestTrack1 = -99999.;
  PtRecMuBestTrack2 = -99999.;

  m_vtxChi2Mu = -99999.;
  m_vtxMassMu = -99999.;
  m_vtxMassSmearedMu = -99999.;
  m_vtxMassScaledMu = -99999.;
  m_scaleUnc  = -99999.;
  m_csAngle   = -99999.;

  m_ptGen1   = -99999.;
  m_phiGen1  = -99999.;
  m_etaGen1  = -99999.;
  m_enGen1   = -99999.;
  m_genFlag1 = -1;

  m_ptGen2  = -99999.;
  m_phiGen2 = -99999.;
  m_etaGen2 = -99999.;
  m_enGen2  = -99999.;

  ChargeRecMu1 = 0;
  ChargeRecMu2 = 0;

  flagmu1 = -1;
  flag1   = -1;

  PtRecTunePMuBestTrack1 = -99999.;
  EnRecMu1  = -99999.;
  EtaRecMu1 = -99999.;
  PhiRecMu1 = -99999.;

  PtRecTunePMuBestTrack2 = -99999.;
  EnRecMu2  = -99999.;
  EtaRecMu2 = -99999.;
  PhiRecMu2 = -99999.;

  pxRecMu1 = -99999.;
  pyRecMu1 = -99999.;
  pzRecMu1 = -99999.;
  pRecMu1  = -99999.;
  dxyRecMu1 = -99999.;

  pxRecMu2 = -99999.;
  pyRecMu2 = -99999.;
  pzRecMu2 = -99999.;
  pRecMu2  = -99999.;
  dxyRecMu2 = -99999.;

  m_genET1  = -99999.;
  m_genPhi1 = -99999.;
  m_genEta1 = -99999.;
  m_genEn1  = -99999.;;

  m_genID1   = -1;
  m_genStat1 = -1;

  m_genET2 = -99999.;
  m_genPhi2 = -99999.;
  m_genEta2 = -99999.;
  m_genEn2 = -99999.;

  m_genID2 = -1;
  m_genStat2 = -1;

  m_genMass = -99999.;
  // seems not used...
  m_recoMass = -99999.;

  nbTP = 0;
  nbTT = 0;
  nbTF = 0;

  TagProbeEtaCut = -99999.;

  Eff = -99999.;

  m_nbFireHLT = 0;
}

void ZprimeMuMuPatMiniAodNewMC::Loop(bool debug)
{
  time_t start,end;
  double dif;
  time (&start);
  FILE * pFile;
  pFile = fopen ("myfile.txt","w");
  //values needed for AccepXeff study
  // m_nbGen = 0;
  // m_nbReco= 0;
  int binMass   = 10000; //100; //100; //100; //100; //100; //100; //100; //100; //100; //10000;
  float minMass = 0.0;
  float maxMass = 10000.0;
  // MassCutMin = 0.0;
  // MassCutMax = 2000.0;
  // EtaCut = 2.4;
  // MassResolution = 0.10;
  // deltaRcut = 0.15;
  // RecoHLTMatchingDeltaRcut = 0.20;
  // minMassCut = 50.0;
  // maxMassCut = 4000.0;
  int ptBins = 40;
  float ptMin = 0.0;
  float ptMax = 400.0;
  // ptEffCut = 3000.0;
  // FR_Ptcut = 53.0; //53.0;

  TFile *f = new TFile("muon_systematics.root","OPEN");
  m_muon_scale_ratio_hist = static_cast<TH2D *>(f->Get("h2_muo_scale"));
  // m_weight = 1.;  // this is dumb, definitely don't want to reset the weight...
  if (DATA_type=="2015" || DATA_type=="2016" || DATA_type=="2017")
    m_weight = 1.;
  std::shared_ptr<TFile> output = std::make_shared<TFile>("ZprimeToMuMu_13TeV.root","recreate");
  // Enable Sumw2 for histograms as we'll be normalizing them
  TH1::SetDefaultSumw2(true);

  //==================================================================================
  //                                                                                 =
  //             Dijet histograms                                                    =
  //                                                                                 =
  //==================================================================================
  h1_DijetEta1_       = std::make_shared<TH1D>("DijetEta1","",100,0.0,3.0);
  h1_DijetEta2_       = std::make_shared<TH1D>("DijetEta2","",100,0.0,3.0);
  // Build the histo with constant log bin width
  const int NMBINS2 = 27;
  const double MMIN2 = 60., MMAX2 = 2000.;
  double logMbins2[NMBINS2+1];
  float binNormNr2=0.;
  for (int ibin = 0; ibin <= NMBINS2; ibin++) {
    logMbins2[ibin] = exp(log(MMIN2) + (log(MMAX2)-log(MMIN2))*ibin/NMBINS2);
    if (debug)
      std::cout << logMbins2[ibin] << std::endl;
  }
  h1_BTagMassMuMuBinWidth_        = std::make_shared<TH1D>("BTagMassMuMuBinWidth","BTagMassMuMuBinWidth",NMBINS2, logMbins2);
  h1_MassMuMuDijetBinWidth_       = std::make_shared<TH1D>("MassMuMuDijetBinWidth","MassMuMuDijetBinWidth",NMBINS2, logMbins2);
  h1_MassMuMuDijetBinWidthMET_    = std::make_shared<TH1D>("MassMuMuDijetBinWidthMET","MassMuMuDijetBinWidthMET",NMBINS2, logMbins2);
  h1_MassMuMuBinWidthMET_         = std::make_shared<TH1D>("MassMuMuBinWidthMET","MassMuMuBinWidthMET",NMBINS2, logMbins2);
  h1_MassMuMuDijet1GeVbin_        = std::make_shared<TH1D>("MassMuMuDijet1GeVbin","",3000,0.0,3000.0);
  h1_MassMuMuDijet1GeVbinMET_     = std::make_shared<TH1D>("MassMuMuDijet1GeVbinMET","",3000,0.0,3000.0);
  h1_MassMuMu1GeVbinMET_          = std::make_shared<TH1D>("MassMuMu1GeVbinMET","",3000,0.0,3000.0);
  h1_MassMuMuDijetBinWidthMET100_ = std::make_shared<TH1D>("MassMuMuDijetBinWidthMET100","MassMuMuDijetBinWidthMET100",NMBINS2, logMbins2);
  h1_MassMuMuDijet1GeVbinMET100_  = std::make_shared<TH1D>("MassMuMuDijet1GeVbinMET100","",3000,0.0,3000.0);

  //==================================================================================
  //                                                                                 =
  //             Start the histograms for CollinSoper CMF                            =
  //                                                                                 =
  //==================================================================================
  m_nbFireHLT = 0;
  int   NbBins = 10;
  float MinBin = -1.0;
  float MaxBin =  1.0;
  h1_ptPFjetsAll_        = std::make_shared<TH1D>("ptPFjetsAll","",100,0.0,2000.0);
  h1_NbPFjetsAll_        = std::make_shared<TH1D>("NbPFjetsAll","",5,0.0,5.0);
  h1_NbPFjets2_          = std::make_shared<TH1D>("NbPFjets2","",5,0.0,5.0);
  h1_BTagMassMuMu_       = std::make_shared<TH1D>("BTagMassMuMu","",100,0.0,3000.0);
  h1_BTagMassMuMu1GeVbin_= std::make_shared<TH1D>("BTagMassMuMu1GeVbin","",3000,0.0,3000.0);
  h1_nbBTagStep1_        = std::make_shared<TH1D>("nbBTagStep1","",5,0.0,5.0);
  h1_jetBTagStep1_       = std::make_shared<TH1D>("jetBTagStep1","",50,0.0,1.0);
  h1_nbBTagStep2_        = std::make_shared<TH1D>("nbBTagStep2","",5,0.0,5.0);
  h1_jetBTagStep2_       = std::make_shared<TH1D>("jetBTagStep2","",50,0.0,1.0);
  h1_nbBTagStep3_        = std::make_shared<TH1D>("nbBTagStep3","",5,0.0,5.0);
  h1_jetBTagStep3_       = std::make_shared<TH1D>("jetBTagStep3","",50,0.0,1.0);
  h1_MissingEt_          = std::make_shared<TH1D>("MissingEt","",100,0.0,2000.0);
  h1_Mt_                 = std::make_shared<TH1D>("Mt","",100,0.0,2000.0);
  h1_DeltaPtoverPt_      = std::make_shared<TH1D>("DeltaPtoverPt","",100,0.0,1.0);
  h1_DeltaPhi_           = std::make_shared<TH1D>("DeltaPhi","",50,-4.0,4.0);
  h1_BosPt_              = std::make_shared<TH1D>("BosPt","",40,0.0,200.0);
  h1_BosPhi_             = std::make_shared<TH1D>("BosPhi","",50,-4.0,4.0);
  // NOT FILLED h1_jetBTagB_           = std::make_shared<TH1D>("jetBTagB","",50,0.0,1.0);
  // NOT FILLED h1_jetBTagC_           = std::make_shared<TH1D>("jetBTagC","",50,0.0,1.0);
  // NOT FILLED h1_jetBTagUDSG_        = std::make_shared<TH1D>("jetBTagUDSG","",50,0.0,1.0);
  h1_PFMetCorr_          = std::make_shared<TH1D>("PFMetCorr","",100,0.0,2000.0);
  h1_CaloMet_            = std::make_shared<TH1D>("CaloMet","",100,0.0,2000.0);
  // NOT FILLED h1_NbFireHLT           = std::make_shared<TH1D>("NbFireHLT", "NbFireHLT", 30 , 0. , 30. );
  // NOT FILLED h1_ptBeforeTrigger_    = std::make_shared<TH1D>("ptBeforeTrigger","",50,0.0,2000.0);
  // NOT FILLED h1_ptAfterTrigger_     = std::make_shared<TH1D>("ptAfterTrigger","",50,0.0,2000.0);
  h1_CosAngleCollinSoperCorrect60Mass120_    = std::make_shared<TH1D>("CosAngleCollinSoperCorrect60Mass120",
								      "",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect120Mass300_   = std::make_shared<TH1D>("CosAngleCollinSoperCorrect120Mass300",
								      "",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect300Mass700_   = std::make_shared<TH1D>("CosAngleCollinSoperCorrect300Mass700",
								      "",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect700Mass3000_  = std::make_shared<TH1D>("CosAngleCollinSoperCorrect700Mass3000",
								      "",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect4900Mass5100_ = std::make_shared<TH1D>("CosAngleCollinSoperCorrect4900Mass5100",
								      "",NbBins,MinBin,MaxBin);
  h1_absCosAngleCollinSoperCorrect4500Mass5500_ = std::make_shared<TH1D>("absCosAngleCollinSoperCorrect4500Mass5500",
									 "",5,0.0,1.0);

  double etaBins[] = {-2.4,-1.2,0.,1.2,2.4};
  std::array<std::string,9> etaBinLabels{"All","BB","BE","EE","BE+","BE-","E+E-","E-E-","E+E+"};
  std::array<std::string,3> csBinLabels{"","CS < 0;","CS > 0;"};
  double massBins[] = {0., 200., 400., 500., 700., 1100., 1900., 3500., 5000.};
  int nEtaBins  = sizeof(etaBins)/sizeof(etaBins[0]);
  int nMassBins = sizeof(massBins)/sizeof(massBins[0]);

  std::array<std::string,4> etaBin{"BB","BE","EE","Inc"};
  std::array<std::string,3> csBin{"Pos","Neg","Inc"};
  // int eb = 0;
  // for (auto e: etaBin) {
  //   int cb = 0
  //   for (auto c: csBin) {
  //     h1_SmearedMassBinned_[cb][eb] = std::make_shared<TH1D>(("CS"+c+"SmearedMass"+e,"", 20000,0,20000);
  //     h1_MassBinned_[cb][eb]        = std::make_shared<TH1D>(("CS"+c+"Mass"+e       ,"", 20000,0,20000);
  //     h1_MassUpBinned_[cb][eb]      = std::make_shared<TH1D>(("CS"+c+"MassUp"+e     ,"", 20000,0,20000);
  //     h1_MassDownBinned_[cb][eb]    = std::make_shared<TH1D>(("CS"+c+"MassDown"+e   ,"", 20000,0,20000);
  //     ++cb;
  //   }
  //   ++eb;
  // }
  h2_CSSmearedMassBinned_ = std::make_shared<TH2D>("CSSmearedMassBinned","", 20000,0.,20000., 30,-0.5,29.5);
  h2_CSMassBinned_        = std::make_shared<TH2D>("CSMassBinned"       ,"", 20000,0.,20000., 30,-0.5,29.5);
  h2_CSMassMuIDBinned_    = std::make_shared<TH2D>("CSMassMuIDBinned"   ,"", 20000,0.,20000., 30,-0.5,29.5);
  h2_CSMassUpBinned_      = std::make_shared<TH2D>("CSMassUpBinned"     ,"", 20000,0.,20000., 30,-0.5,29.5);
  h2_CSMassDownBinned_    = std::make_shared<TH2D>("CSMassDownBinned"   ,"", 20000,0.,20000., 30,-0.5,29.5);
  // h2_CSPosSmearedMassBinned_ = std::make_shared<TH2D>("CSPosSmearedMassBinned","", 20000,0.,20000., 30,-0.5,29.5);
  // h2_CSPosMassBinned_        = std::make_shared<TH2D>("CSPosMassBinned"       ,"", 20000,0.,20000., 30,-0.5,29.5);
  // h2_CSPosMassUpBinned_      = std::make_shared<TH2D>("CSPosMassUpBinned"     ,"", 20000,0.,20000., 30,-0.5,29.5);
  // h2_CSPosMassDownBinned_    = std::make_shared<TH2D>("CSPosMassDownBinned"   ,"", 20000,0.,20000., 30,-0.5,29.5);
  // h2_CSNegSmearedMassBinned_ = std::make_shared<TH2D>("CSNegSmearedMassBinned","", 20000,0.,20000., 30,-0.5,29.5);
  // h2_CSNegMassBinned_        = std::make_shared<TH2D>("CSNegMassBinned"       ,"", 20000,0.,20000., 30,-0.5,29.5);
  // h2_CSNegMassUpBinned_      = std::make_shared<TH2D>("CSNegMassUpBinned"     ,"", 20000,0.,20000., 30,-0.5,29.5);
  // h2_CSNegMassDownBinned_    = std::make_shared<TH2D>("CSNegMassDownBinned"   ,"", 20000,0.,20000., 30,-0.5,29.5);
  for (int eb = 0; eb < etaBinLabels.size(); ++eb) {
    for (int cb = 0; cb < csBinLabels.size(); ++cb) {
      h2_CSSmearedMassBinned_->GetYaxis()->SetBinLabel((eb*csBinLabels.size())+cb+1, (csBinLabels[cb]+etaBinLabels[eb]).c_str());
      h2_CSMassBinned_       ->GetYaxis()->SetBinLabel((eb*csBinLabels.size())+cb+1, (csBinLabels[cb]+etaBinLabels[eb]).c_str());
      h2_CSMassMuIDBinned_   ->GetYaxis()->SetBinLabel((eb*csBinLabels.size())+cb+1, (csBinLabels[cb]+etaBinLabels[eb]).c_str());
      h2_CSMassUpBinned_     ->GetYaxis()->SetBinLabel((eb*csBinLabels.size())+cb+1, (csBinLabels[cb]+etaBinLabels[eb]).c_str());
      h2_CSMassDownBinned_   ->GetYaxis()->SetBinLabel((eb*csBinLabels.size())+cb+1, (csBinLabels[cb]+etaBinLabels[eb]).c_str());
      // h2_CSPosSmearedMassBinned_->GetYaxis()->SetBinLabel(eb+1+cb, (csBinLabels[cb]+", "+etaBinLabels[eb]).c_str());
      // h2_CSPosMassBinned_       ->GetYaxis()->SetBinLabel(eb+1+cb, (csBinLabels[cb]+", "+etaBinLabels[eb]).c_str());
      // h2_CSPosMassUpBinned_     ->GetYaxis()->SetBinLabel(eb+1+cb, (csBinLabels[cb]+", "+etaBinLabels[eb]).c_str());
      // h2_CSPosMassDownBinned_   ->GetYaxis()->SetBinLabel(eb+1+cb, (csBinLabels[cb]+", "+etaBinLabels[eb]).c_str());
      // h2_CSNegSmearedMassBinned_->GetYaxis()->SetBinLabel(eb+1+cb, (csBinLabels[cb]+", "+etaBinLabels[eb]).c_str());
      // h2_CSNegMassBinned_       ->GetYaxis()->SetBinLabel(eb+1+cb, (csBinLabels[cb]+", "+etaBinLabels[eb]).c_str());
      // h2_CSNegMassUpBinned_     ->GetYaxis()->SetBinLabel(eb+1+cb, (csBinLabels[cb]+", "+etaBinLabels[eb]).c_str());
      // h2_CSNegMassDownBinned_   ->GetYaxis()->SetBinLabel(eb+1+cb, (csBinLabels[cb]+", "+etaBinLabels[eb]).c_str());
    }
  }
  //==================================================================================
  // NOT USED h1_ptHistoBefor_               = std::make_shared<TH1D>("ptHistoBefor","",ptBins,ptMin,ptMax);
  // NOT USED h1_ptHistoPassingVtxChi2Mu_    = std::make_shared<TH1D>("ptHistoPassingVtxChi2Mu","",ptBins,ptMin,ptMax);
  // NOT USED h1_ptHistoPassingCosmicRejec_  = std::make_shared<TH1D>("ptHistoPassingCosmicRejec","",ptBins,ptMin,ptMax);
  // NOT USED h1_ptHistoPassingHLT_          = std::make_shared<TH1D>("ptHistoPassingHLT","",ptBins,ptMin,ptMax);
  // NOT USED h1_etaHistoBefor_              = std::make_shared<TH1D>("etaHistoBefor","",30,0.0,3.0);
  // NOT USED h1_etaHistoPassingVtxChi2Mu_   = std::make_shared<TH1D>("etaHistoPassingVtxChi2Mu","",30,0.0,3.0);
  // NOT USED h1_etaHistoPassingCosmicRejec_ = std::make_shared<TH1D>("etaHistoPassingCosmicRejec","",30,0.0,3.0);
  // NOT USED h1_etaHistoPassingHLT_         = std::make_shared<TH1D>("etaHistoPassingHLT","",30,0.0,3.0);
  //==================================================================================
  //                                                                                 =
  //            Start a histograms for Mass resolution                               =
  //                                                                                 =
  //==================================================================================
  h1_MassResultionEBEB1_      = std::make_shared<TH1D>("MassResultionEBEB1","",100,-0.5,0.5);
  h1_MassResultionEBEB2_      = std::make_shared<TH1D>("MassResultionEBEB2","",100,-0.5,0.5);
  h1_MassResultionEBEB3_      = std::make_shared<TH1D>("MassResultionEBEB3","",100,-0.5,0.5);
  h1_MassResultionEBEB4_      = std::make_shared<TH1D>("MassResultionEBEB4","",100,-0.5,0.5);
  h1_MassResultionEBEB5_      = std::make_shared<TH1D>("MassResultionEBEB5","",100,-0.5,0.5);
  h1_MassResultionEBEB6_      = std::make_shared<TH1D>("MassResultionEBEB6","",100,-0.5,0.5);
  h1_MassResultionEBEB7_      = std::make_shared<TH1D>("MassResultionEBEB7","",100,-0.5,0.5);
  // NOT USED h1_MassResultionEBEB8_      = std::make_shared<TH1D>("MassResultionEBEB8","",100,-0.5,0.5);
  // NOT USED h1_MassResultionEBEB9_      = std::make_shared<TH1D>("MassResultionEBEB9","",100,-0.5,0.5);
  // NOT USED h1_MassResultionEBEB10_     = std::make_shared<TH1D>("MassResultionEBEB10","",100,-0.5,0.5);
  //==================================================================================
  //                                                                                 =
  //            Start the histograms for the mass of Z                               =
  //                                                                                 =
  //==================================================================================
  h1_ZprimeRecomasslogscale_     = std::make_shared<TH1D>("ZprimeRecomasslogscale","",42,log10(60.0),log10(1200.0));
  // NOT USED h2_ZprimeRecomassNewbin_       = std::make_shared<TH2D>("ZprimeRecomassNewbin","ZprimeRecomassNewbin",80,50.,1000.,500,0.001,1000.);
  h1_MassGenInAccep_   = std::make_shared<TH1D>("MassGenInAccep","",58,0.0,5000.0);
  h1_MassRecoInAccep_  = std::make_shared<TH1D>("MassRecoInAccep","",58,0.0,5000.0);
  h1_ZprimeRecomass_   = std::make_shared<TH1D>("ZprimeRecomass","",binMass,minMass,maxMass);
  h1_ZprimeRecomass20_ = std::make_shared<TH1D>("ZprimeRecomass20","",300,0.0,6000.0);
  h1_ZprimeRecomassBB_ = std::make_shared<TH1D>("ZprimeRecomassBB","",binMass,minMass,maxMass);
  h1_ZprimeRecomassEE_ = std::make_shared<TH1D>("ZprimeRecomassEE","",binMass,minMass,maxMass);
  h1_ZprimeRecomassBE_ = std::make_shared<TH1D>("ZprimeRecomassBE","",binMass,minMass,maxMass);
  h1_ZprimeRecomass50_ = std::make_shared<TH1D>("ZprimeRecomass50","",120,0.0,6000.0);
  // h1_ZprimeRecomass60to120_      = std::make_shared<TH1D>("ZprimeRecomass60to120","",binMass,minMass,maxMass);
  // NOT USED h1_ZprimeRecomassAbove1000GeV_ = std::make_shared<TH1D>("ZprimeRecomassAbove1000GeV","",binMass,minMass,maxMass);
  h1_ZprimeGenmass_ = std::make_shared<TH1D>("ZprimeGenmass","",binMass,minMass,maxMass);
  h1_ZprimeGenEta1_ = std::make_shared<TH1D>("ZprimeGenEta1","",100,-8.0,8.0);
  h1_ZprimeGenEta2_ = std::make_shared<TH1D>("ZprimeGenEta2","",100,-8.0,8.0);
  h1_ZprimeGenPt1_  = std::make_shared<TH1D>("ZprimeGenPt1","",100,0.0,2000.0);
  h1_ZprimeGenPt2_  = std::make_shared<TH1D>("ZprimeGenPt2","",100,0.0,2000.0);
  h1_ZprimeGenEn1_  = std::make_shared<TH1D>("ZprimeGenEn1","",100,0.0,2000.0);
  h1_ZprimeGenEn2_  = std::make_shared<TH1D>("ZprimeGenEn2","",100,0.0,2000.0);
  h1_3Dangle_       = std::make_shared<TH1D>("3Dangle","",100,-2.0,2.0);
  // NOT USED h1_DxyDiff_                    = std::make_shared<TH1D>("DxyDiff","",100,10.0,10.0);
  // NOT USED h1_MassRecoGenDif_             = std::make_shared<TH1D>("MassRecoGenDif","",100,-0.5,0.5);
  h1_PtResolutionTunePMBT_       = std::make_shared<TH1D>("PtResolutionTunePMBT","",100,-0.5,0.5);
  h1_PtResolutiontuneP_          = std::make_shared<TH1D>("PtResolutiontuneP","",100,-0.5,0.5);
  h1_PtResolutionMBT_            = std::make_shared<TH1D>("PtResolutionMBT","",100,-0.5,0.5);
  //==================================================================================
  //                                                                                 =
  //                 Start the histograms for N-1 dist                               =
  //                                                                                 =
  //==================================================================================
  h1_dPToverPT_                            = std::make_shared<TH1D>("dPToverPT","",100,0.0,0.5);
  // NOT USED h1_normalizedChi2_                       = std::make_shared<TH1D>("normalizedChi2","",100,0.0,20.0);
  h1_numberOftrackerLayersWithMeasurement_ = std::make_shared<TH1D>("numberOftrackerLayersWithMeasurement","",20,0.0,20.0);
  h1_numberOfValidPixelHits_               = std::make_shared<TH1D>("numberOfValidPixelHits","",10,0.0,10.0);
  h1_numberOfValidMuonHits_                = std::make_shared<TH1D>("numberOfValidMuonHits","",60,0.0,60.0);
  h1_numberOfMatchedStations_              = std::make_shared<TH1D>("numberOfMatchedStations","",10,0.0,10.0);
  h1_trackiso_                             = std::make_shared<TH1D>("trackiso","",50,0.0,0.3);
  h1_absdxy_                               = std::make_shared<TH1D>("absdxy","",100,0.0,0.3);
  h1_PtEffpterror_                   = std::make_shared<TH1D>("PtEffpterror","",ptBins,ptMin,ptMax);
  h1_PtEffptnumberOftrackerLayers_   = std::make_shared<TH1D>("PtEffptnumberOftrackerLayers","",ptBins,ptMin,ptMax);
  h1_PtEffptnumberOfPixelHits_       = std::make_shared<TH1D>("PtEffptnumberOfPixelHits","",ptBins,ptMin,ptMax);
  h1_PtEffptnumberOfMuonHits_        = std::make_shared<TH1D>("PtEffptnumberOfMuonHits","",ptBins,ptMin,ptMax);
  h1_PtEffptnumberOfMatchedStations_ = std::make_shared<TH1D>("PtEffptnumberOfMatchedStations","",ptBins,ptMin,ptMax);
  h1_PtEffptTrackIso_                = std::make_shared<TH1D>("PtEffptTrackIso","",ptBins,ptMin,ptMax);
  h1_PtEffptabsdsy_                  = std::make_shared<TH1D>("PtEffptabsdsy","",ptBins,ptMin,ptMax);
  // NOT USED h1_PtEffpfSumChargedHadron_        = std::make_shared<TH1D>("PtEffpfSumChargedHadron","",ptBins,ptMin,ptMax);
  // NOT USED h1_PtEffpfSumNeutralHadron_        = std::make_shared<TH1D>("PtEffpfSumNeutralHadron","",ptBins,ptMin,ptMax);
  // NOT USED h1_PtEffpfPhotonIso_               = std::make_shared<TH1D>("PtEffpfPhotonIso","",ptBins,ptMin,ptMax);
  h1_EtaEffpterror_                  = std::make_shared<TH1D>("EtaEffpterror","",30,0.0,3.0);
  h1_EtaEffptnumberOftrackerLayers_  = std::make_shared<TH1D>("EtaEffptnumberOftrackerLayers","",30,0.0,3.0);
  h1_EtaEffptnumberOfPixelHits_      = std::make_shared<TH1D>("EtaEffptnumberOfPixelHits","",30,0.0,3.0);
  h1_EtaEffptnumberOfMuonHits_       = std::make_shared<TH1D>("EtaEffptnumberOfMuonHits","",30,0.0,3.0);
  h1_EtaEffptnumberOfMatchedStations_= std::make_shared<TH1D>("EtaEffptnumberOfMatchedStations","",30,0.0,3.0);
  h1_EtaEffptTrackIso_               = std::make_shared<TH1D>("EtaEffptTrackIso","",30,0.0,3.0);
  h1_EtaEffptabsdsy_                 = std::make_shared<TH1D>("EtaEffptabsdsy","",30,0.0,3.0);
  h1_nbPVID_               = std::make_shared<TH1D>("nbPVID","",50,0.0,50.0);
  h1_PtID_                 = std::make_shared<TH1D>("PtID","",ptBins,ptMin,ptMax);
  h1_EtaID_                = std::make_shared<TH1D>("EtaID","",30,0.0,3.0);
  h1_nbPVNewID_            = std::make_shared<TH1D>("nbPVNewID","",50,0.0,50.0);
  h1_PtNewID_              = std::make_shared<TH1D>("PtNewID","",ptBins,ptMin,ptMax);
  h1_EtaNewID_             = std::make_shared<TH1D>("EtaNewID","",30,0.0,3.0);
  h1_nbPVTightID_          = std::make_shared<TH1D>("nbPVTightID","",50,0.0,50.0);
  h1_PtTightID_            = std::make_shared<TH1D>("PtTightID","",ptBins,ptMin,ptMax);
  h1_EtaTightID_           = std::make_shared<TH1D>("EtaTightID","",30,0.0,3.0);
  // h1_3DangleHisto1_        = std::make_shared<TH1D>("3DangleHisto1","",50,0.00001,0.1);
  // h1_3DangleHisto2_        = std::make_shared<TH1D>("3DangleHisto2","",100,0.00001,10.0);  //100,0.1,10.0);

  // NOT USED h1_3DangleHisto1_        = std::make_shared<TH1D>("3DangleHisto1","",1000,0.00001,1.0);
  // NOT USED h1_3DangleHisto2_        = std::make_shared<TH1D>("3DangleHisto2","",1000,0.00001,10.0);
  // NOT USED h1_Fail3DangleHistoMass_ = std::make_shared<TH1D>("Fail3DangleHistoMass","",100,100.0,12000.0);
  // NOT USED h1_Fail3DangleHistoPhi_  = std::make_shared<TH1D>("Fail3DangleHistoPhi","",100,-4.0,4.0);
  //-----------------------------------------------------------------------------------
  // NOT USED h1_ptHistoFRDum_ = std::make_shared<TH1D>("ptHistoFRDum","",ptBins,ptMin,ptMax);
  // NOT USED h1_ptHistoFRNum_ = std::make_shared<TH1D>("ptHistoFRNum","",ptBins,ptMin,ptMax);
  h1_PtTuneP_      = std::make_shared<TH1D>("PtTuneP","",200,0.0,10000.0);
  //----------------------------------------------------------------------------------------------

  // Build the histo with constant log bin width
  const int NMBINS = 100;
  const double MMIN = 60., MMAX = 6000.;
  double logMbins[NMBINS+1];
  float binNormNr=0.;

  for (int ibin = 0; ibin <= NMBINS; ibin++) {
    logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
    //std::cout << logMbins[ibin] << std::endl;
  }
  h1_ZprimeRecomassBinWidthAll_     = std::make_shared<TH1D>("ZprimeRecomassBinWidthAll",
							     "ZprimeRecomassBinWidthAll",NMBINS, logMbins);
  h1_ZprimeRecomassBinWidth_        = std::make_shared<TH1D>("ZprimeRecomassBinWidth",
							     "ZprimeRecomassBinWidth",NMBINS, logMbins);
  h1_ZprimeRecomassBinWidthAllBE_   = std::make_shared<TH1D>("ZprimeRecomassBinWidthAllBE",
							     "ZprimeRecomassBinWidthAllBE",NMBINS, logMbins);
  h1_ZprimeRecomassBinWidthAllEE_   = std::make_shared<TH1D>("ZprimeRecomassBinWidthAllEE",
							     "ZprimeRecomassBinWidthAllEE",NMBINS, logMbins);
  h1_ZprimeRecomassBinWidthEE_      = std::make_shared<TH1D>("ZprimeRecomassBinWidthEE",   "",NMBINS, logMbins);
  h1_ZprimeRecomassBinWidthBB_      = std::make_shared<TH1D>("ZprimeRecomassBinWidthBB",   "",NMBINS, logMbins);
  h1_ZprimeRecomassBinWidthBEpos_   = std::make_shared<TH1D>("ZprimeRecomassBinWidthBEpos","",NMBINS, logMbins);
  h1_ZprimeRecomassBinWidthBEnev_   = std::make_shared<TH1D>("ZprimeRecomassBinWidthBEnev","",NMBINS, logMbins);
  h1_ZprimeRecomass60to120BEpos_    = std::make_shared<TH1D>("ZprimeRecomass60to120BEpos", "",60,60.0,120.0);
  h1_ZprimeRecomass60to120BEnev_    = std::make_shared<TH1D>("ZprimeRecomass60to120BEnev", "",60,60.0,120.0);
  h1_ZprimeRecomass60to120EE_       = std::make_shared<TH1D>("ZprimeRecomass60to120EE",    "",60,60.0,120.0);
  h1_ZprimeRecomass60to120BB_       = std::make_shared<TH1D>("ZprimeRecomass60to120BB",    "",60,60.0,120.0);
  h1_ZprimeRecomass60to120_         = std::make_shared<TH1D>("ZprimeRecomass60to120",      "",60,60.0,120.0);
  h1_ZprimeRecomassBinWidthAfterBtaging_ = std::make_shared<TH1D>("ZprimeRecomassBinWidthAfterBtaging",
								  "ZprimeRecomassBinWidthAfterBtaging",NMBINS, logMbins);
  h1_DijetBinWidthBB_   = std::make_shared<TH1D>("DijetBinWidthBB",  "",NMBINS, logMbins);
  h1_DijetBinWidthBE_   = std::make_shared<TH1D>("DijetBinWidthBE",  "",NMBINS, logMbins);
  h1_DijetBinWidthEE_   = std::make_shared<TH1D>("DijetBinWidthEE",  "",NMBINS, logMbins);
  h1_DijetBinWidthBBBE_ = std::make_shared<TH1D>("DijetBinWidthBBBE","",NMBINS, logMbins);
  h1_WjetsBinWidthBB_   = std::make_shared<TH1D>("WjetsBinWidthBB",  "",NMBINS, logMbins);
  h1_WjetsBinWidthBE_   = std::make_shared<TH1D>("WjetsBinWidthBE",  "",NMBINS, logMbins);
  h1_WjetsBinWidthEE_   = std::make_shared<TH1D>("WjetsBinWidthEE",  "",NMBINS, logMbins);
  h1_WjetsBinWidthBBBE_ = std::make_shared<TH1D>("WjetsBinWidthBBBE","",NMBINS, logMbins);

  h1_Dijet1GeVBB_     = std::make_shared<TH1D>("Dijet1GeVBB",    "",3000,0.0,3000.0);
  h1_Dijet1GeVBEEE_   = std::make_shared<TH1D>("Dijet1GeVBEEE",  "",3000,0.0,3000.0);
  h1_Dijet1GeVEE_     = std::make_shared<TH1D>("Dijet1GeVEE",    "",3000,0.0,3000.0);
  h1_Dijet1GeVBBBEEE_ = std::make_shared<TH1D>("Dijet1GeVBBBEEE","",3000,0.0,3000.0);
  h1_Wjets1GeVBB_     = std::make_shared<TH1D>("Wjets1GeVBB",    "",3000,0.0,3000.0);
  h1_Wjets1GeVBEEE_   = std::make_shared<TH1D>("Wjets1GeVBEEE",  "",3000,0.0,3000.0);
  h1_Wjets1GeVEE_     = std::make_shared<TH1D>("Wjets1GeVEE",    "",3000,0.0,3000.0);
  h1_Wjets1GeVBBBEEE_ = std::make_shared<TH1D>("Wjets1GeVBBBEEE","",3000,0.0,3000.0);

  h1_Dijet20GeVBB_     = std::make_shared<TH1D>("Dijet20GeVBB",    "",300,0.0,6000.0);
  h1_Dijet20GeVBEEE_   = std::make_shared<TH1D>("Dijet20GeVBEEE",  "",300,0.0,6000.0);
  h1_Dijet20GeVBBBEEE_ = std::make_shared<TH1D>("Dijet20GeVBBBEEE","",300,0.0,6000.0);
  h1_Wjets20GeVBB_     = std::make_shared<TH1D>("Wjet20GeVBB",     "",300,0.0,6000.0);
  h1_Wjets20GeVBEEE_   = std::make_shared<TH1D>("Wjets20GeVBEEE",  "",300,0.0,6000.0);
  h1_Wjets20GeVBBBEEE_ = std::make_shared<TH1D>("Wjets20GeVBBBEEE","",300,0.0,6000.0);

   // Pileup reweighting 2016 data vs Spring16 MC in 80x
  TFile *filePU= TFile::Open("puWeightsMoriond17_v2.root");
  TH1F *puweight = (TH1F*)filePU->Get("weights");

  // Book txt file for candidate events
  Char_t txtOUT[500];
  sprintf(txtOUT,"ZprimeToMuMu_13TeV_cand.txt");
  output_txt.open(txtOUT);
  output_txt << "CANDIDATES Events:" << std::endl;
  Char_t outform[20000];
  sprintf (outform,"run: lumi: event: dil_mass: pTmu1: pTmu2: Etamu1: Etamu2:");
  output_txt  << outform << std::endl;

  TString inputfile=name;
  inputfile=name;
  std::cout << "Name of the input file is= " << inputfile.Data() << std::endl;
  std::cout << "Weight of the sample is= " << m_weight << std::endl;

  //==================================================================================
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  if (debug)
    nentries = 1000;

  // Timing information
  int decade  = 0;
  int century = 0;
  TStopwatch tsw;
  int tenpcount = 1;
  int onepcount = 1;
  int tenthpcount = 1;
  std::streamsize defpres = std::cout.precision();

  int wwto2l2nu_input(0),wwto2l2nu_fail_gen_mass(0),wwto2l2nu_fail_reco_mass(0);
  int ttto2l2nu_input(0),ttto2l2nu_fail_gen_mass(0),ttto2l2nu_fail_reco_mass(0);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //Timing information
    if (jentry==0) {
      tsw.Start();
      std::cout << "." << std::flush;
    }
    if ((jentry*10)/nentries == tenpcount ) {
      tsw.Stop();
      Double_t time = tsw.RealTime();
      tsw.Start(kFALSE);
      Double_t finTime(0.);
      Double_t frac = (jentry*1.0)/(nentries*1.0);
      if (frac>0)
	finTime = time / frac - time;
      Double_t finMin = finTime / 60.;
      std::cout << tenpcount*10 << "% done.  "
		<< "t = "  << std::setprecision(4) << std::setw(7) << time
		<< " projected finish =" << std::setw(7) << std::setprecision(4) << finTime << "s"
		<< " (" << std::setw(4) << std::setprecision(2) << finMin << " min).   "
		<< std::setprecision(defpres) << std::resetiosflags << std::endl;
      std::cout << std::flush;
      tenpcount++;
      // } else if ( (jentry*100)/nentries == onepcount ) {
      //   std::cout << ".";
      //   std::cout << std::flush;
      //   onepcount++;
      // }
    } else if ( (jentry*1000)/nentries == tenthpcount ) {
      std::cout << ".";
      std::cout << std::flush;
      tenthpcount++;
    }
    ///

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    newweight = m_weight;

    // Pileup Reweighting
    Int_t binx = puweight->GetXaxis()->FindBin(num_PU_vertices);
    pu_weight=double(puweight->GetBinContent(binx));
    // Changing the weight for pileup
    newweight = m_weight*pu_weight;
    if (debug)
      std::cout << "Starting weight + pileup, old weight and new = " << m_weight << " * " << newweight << std::endl;

    if (isCISample) {
      // have to choose which cut to use
      // if (!passMInvCut)
      // if (!passPreFSRMInvCut) {
      // if (!passST1MInvCut) {
      // if (!passST23MInvCut) {
      if (!passHSMInvCut) {
	if (debug)
	  std::cout << "failed CI gen cut" << std::endl;
	continue;
      }
    }

    // if (inputfile.Contains("amcatnlo") && !(inputfile.Contains("DYJetsToTauTau"))) {
    //   std::cout << "Reweighting sample of amcatnlo with weight= " << MC_weighting->at(0)<< std::endl;
    //   m_weight = m_weight*MC_weighting->at(0);
    // }
    // if (event_runNo==251249) continue;
    /*std<<cout << "=======> jentry = " << jentry
                << "=======> Evt = "    << event_evtNo
                << "=======> Run = "    << event_runNo
                << "=======> Lumi = "   << event_lumi
                << "=======> bunch = "  << event_bunch
	        << std::endl;
    */
    //==================================================================
    //                                                                 =
    //               Calling methods to get 2FR estimate               =
    //                                                                 =
    //==================================================================
    DrawDiJetMassBB();
    DrawDiJetMassBE();
    DrawDiJetMassEE();
    //==================================================================
    //                                                                 =
    //               Calling methods to get 1FR estimate               =
    //                                                                 =
    //==================================================================
    DrawWJetsMassBB();
    DrawWJetsMassBE1();
    DrawWJetsMassBE2();
    DrawWJetsMassEE();
    //=========================================================
    //                                                        =
    // Calling methods to get events with 2 muons passing ID  =
    //                                                        =
    //=========================================================
    bool firstMuFinal  = SelectFirstMuon(PtRecTunePMuBestTrack1,EnRecMu1,EtaRecMu1,PhiRecMu1,ChargeRecMu1,flagmu1,
					 pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,dxyRecMu1,PtRecTunePMu1,
					 PtRecMuBestTrack1,
					 m_genET1, m_genEta1, m_genPhi1, m_genEn1);

    bool secondMuFinal = SelectSecondMuon(ChargeRecMu1,flagmu1,PtRecTunePMuBestTrack1,EtaRecMu1,PhiRecMu1,
					  PtRecTunePMuBestTrack2,EnRecMu2,
					  EtaRecMu2,PhiRecMu2,ChargeRecMu2,pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2,
					  PtRecTunePMu2,PtRecMuBestTrack2,
					  m_genET2, m_genEta2, m_genPhi2, m_genEn2);

    PickThehighestMass(m_vtxMassMu,m_vtxChi2Mu,event_evtNo);

    m_genMass = GenMass(m_genET1, m_genPhi1, m_genEta1, m_genEn1,
			m_genET2, m_genPhi2, m_genEta2, m_genEn2);
    m_vtxMassSmearedMu = smearedMass(EtaRecMu1, PhiRecMu1, PtRecTunePMuBestTrack1, EtaRecMu2, PhiRecMu2, PtRecTunePMuBestTrack2, m_vtxMassMu);
    m_vtxMassScaledMu  = scaledMass(EtaRecMu1, PhiRecMu1, PtRecTunePMuBestTrack1, ChargeRecMu1,  EtaRecMu2, PhiRecMu2, PtRecTunePMuBestTrack2, ChargeRecMu2, m_vtxMassMu);

    double CosmicRejec = ThreeDangle(pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,
				     pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2);
    //=========================================================
    //        call the method for N-1 plots                   =
    //                                                        =
    //=========================================================
    //std::cout << "firstMu= " << firstMuFinal << " " << "secondMu= " << secondMuFinal << std::endl;
    if (firstMuFinal == 0 || secondMuFinal == 0) continue;

    //std::cout << "Vertex mass mu= " << m_vtxMassMu << std::endl;
    //if (m_vtxMassMu<60 || m_vtxMassMu>1200) continue;
    if (m_vtxMassMu<60) continue;

    //=========================================================
    //        start doing matching between reco & HLT         =
    //                                                        =
    //=========================================================
    bool fireHLT2 = isPassHLT();

    if (fireHLT2 == 0) continue;

    bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(EtaRecMu1,PhiRecMu1);
    bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(EtaRecMu2,PhiRecMu2);

    if (RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
      plotAllHighPtMuonsID();
      //PrintEventInformation(256843,465,665539990,m_vtxChi2Mu,m_vtxMassMu,CosmicRejec);
      if (m_vtxChi2Mu<20.0 && CosmicRejec>-0.9998) {
        DrawBTaggingDiscriminator(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
        if (PFMet_et_cor > 50.0) {
          h1_PFMetCorr_->Fill(PFMet_et_cor,m_weight);
          h1_CaloMet_->Fill(CaloMet_pt,m_weight);
          h1_MassMuMuBinWidthMET_->Fill(m_vtxMassMu,m_weight);
          h1_MassMuMu1GeVbinMET_->Fill(m_vtxMassMu,m_weight);
        }
        bool passDijet = DiPFJet(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
        if (passDijet==1) {
          h1_MassMuMuDijetBinWidth_->Fill(m_vtxMassMu,m_weight);
          h1_MassMuMuDijet1GeVbin_->Fill(m_vtxMassMu,m_weight);
        }

        bool passDijetcuts = DiPFJetCut(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
        if (passDijetcuts==1 && PFMet_et_cor > 50.0) {
          h1_MassMuMuDijetBinWidthMET_->Fill(m_vtxMassMu,m_weight);
          h1_MassMuMuDijet1GeVbinMET_->Fill(m_vtxMassMu,m_weight);
        }
        if (passDijetcuts==1 && PFMet_et_cor > 100.0) {
          h1_MassMuMuDijetBinWidthMET100_->Fill(m_vtxMassMu,m_weight);
          h1_MassMuMuDijet1GeVbinMET100_->Fill(m_vtxMassMu,m_weight);
        }
        bool passBTaggingDiscriminator2 = BTaggingDiscriminator2(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
        if (passBTaggingDiscriminator2==1) {
          h1_BTagMassMuMu_->Fill(m_vtxMassMu,m_weight);
        }
        bool passBTaggingDiscriminator3 = BTaggingDiscriminator3(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
        if (passBTaggingDiscriminator3==1) {
          h1_BTagMassMuMu_->Fill(m_vtxMassMu,m_weight);
        }

        // Special inclusively binned samples, used only in the low mass region?
        //  - Scaling is off...
        // WW sample
        // need to keep a running track of how many events are rejcted, and adjust the weight appropriately
        if (inputfile.Contains("WWTo2L2Nu_13TeV")) {
          wwto2l2nu_input++;
          if (m_genMass > 600.) {
            if (debug)
              std::cout << "Reweighting sample of WWTo2L2Nu with weight=0: gen("
                        << m_genMass << ") reco("
                        << m_vtxMassMu << ")" << std::endl;
            // newweight = 0;
            wwto2l2nu_fail_gen_mass++;
          }
          if (m_vtxMassMu > 600.) {  // WHY CUT ON THE RECO MASS???
            if (debug)
              std::cout << "Reweighting sample of WWTo2L2Nu with weight=0: gen("
                        << m_genMass << ") reco("
                        << m_vtxMassMu << ")" << std::endl;
            newweight = 0;
            wwto2l2nu_fail_reco_mass++;
          }
        } else if (inputfile.Contains("WWTo2L2Nu_Mll")) {
          wwto2l2nu_input++;
          if (m_genMass < 600.) {
            if (debug)
              std::cout << "Reweighting sample of WWTo2L2Nu_Mll with weight=0: gen("
                        << m_genMass << ") reco("
                        << m_vtxMassMu << ")" << std::endl;
            // newweight = 0;
            wwto2l2nu_fail_gen_mass++;
          }
          if (m_vtxMassMu < 600.) {  // WHY CUT ON THE RECO MASS???
            if (debug)
              std::cout << "Reweighting sample of WWTo2L2Nu_Mll with weight=0: gen("
                        << m_genMass << ") reco("
                        << m_vtxMassMu << ")" << std::endl;
            newweight = 0;
            wwto2l2nu_fail_reco_mass++;
          }
        } else if (inputfile.Contains("TTTo2L2Nu_Tune")) {
          // TTTo2L2Nu sample
	  if (debug)
	    std::cout << "Checking reweighting of inclusive TTTo2L2Nu sample: gen("
		      << m_genMass << ") reco("
		      << m_vtxMassMu << ")" << std::endl;
          ttto2l2nu_input++;
          if (m_genMass > 500.) {
            if (debug)
	      std::cout << "Reweighting inclusive TTTo2L2Nu sample with weight=0: gen("
                        << m_genMass << ") reco("
                        << m_vtxMassMu << ")" << std::endl;
            // newweight = 0;
            ttto2l2nu_fail_gen_mass++;
          }
          if (m_vtxMassMu > 500.) {  // WHY CUT ON THE RECO MASS???
            if (debug)
	      std::cout << "Reweighting inclusive TTTo2L2Nu sample with weight=0: gen("
                        << m_genMass << ") reco("
                        << m_vtxMassMu << ")" << std::endl;
            newweight = 0;
            ttto2l2nu_fail_reco_mass++;
          }
	} else if (inputfile.Contains("TTTo2L2Nu_M") || inputfile.Contains("TTToLL_MLL_")) {
	  if (debug)
	    std::cout << "Checking reweighting of mass binned TTTo2L2Nu sample: gen("
		      << m_genMass << ") reco("
		      << m_vtxMassMu << ")" << std::endl;
          ttto2l2nu_input++;
          if (m_genMass > 500.) {
            if (debug)
	      std::cout << "Reweighting mass binned TTTo2L2Nu sample with weight=0: gen("
                        << m_genMass << ") reco("
                        << m_vtxMassMu << ")" << std::endl;
            // newweight = 0;
            ttto2l2nu_fail_gen_mass++;
          }
          if (m_vtxMassMu<500.) {
            if (debug)
	      std::cout << "Reweighting mass binned TTTo2L2Nu sample with weight=0: gen("
                        << m_genMass << ") reco("
                        << m_vtxMassMu << ")" << std::endl;
            newweight = 0;
            ttto2l2nu_fail_reco_mass++;
          }
        }


        Boson(pxRecMu1,pyRecMu1,pzRecMu1,EnRecMu1,pxRecMu2,pyRecMu2,pzRecMu2,EnRecMu2,
              ChargeRecMu1,PFMet_et_cor,PFMet_px_cor,PFMet_py_cor,PFMet_pz_cor,PFMet_en_cor, m_bosonPt);
        PlotRecoInfo(CosmicRejec,m_vtxMassMu,m_genMass,PtRecTunePMuBestTrack1,PtRecTunePMu1,PtRecMuBestTrack1,m_ptGen1,EtaRecMu1, pRecMu1,
                     PtRecTunePMuBestTrack2,PtRecTunePMu2,PtRecMuBestTrack2,m_ptGen2,EtaRecMu2, pRecMu2, m_bosonPt, inputfile);
        PlotGenInfo(m_genMass,m_genEta1,m_genEta2,m_genET1,m_genET2,m_genEn1,m_genEn2);
        m_csAngle = CosThetaCollinSoper(PtRecTunePMuBestTrack1,EtaRecMu1,PhiRecMu1,EnRecMu1,
                                        PtRecTunePMuBestTrack2,EtaRecMu2,PhiRecMu2,EnRecMu2,
                                        ChargeRecMu1,m_vtxMassMu);
        bool passBTaggingDiscriminator = BTaggingDiscriminator(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
        if (passBTaggingDiscriminator==1) {
          h1_BTagMassMuMuBinWidth_->Fill(m_vtxMassMu,m_weight);
          h1_BTagMassMuMu1GeVbin_->Fill(m_vtxMassMu,m_weight);
        }
        if (passBTaggingDiscriminator==0) {
          h1_ZprimeRecomassBinWidthAfterBtaging_->Fill(m_vtxMassMu,m_weight);
        }
      }
    }
  }

  std::cout << "100% done!" << std::endl;
  std::cout << std::flush;

  if (inputfile.Contains("TTToLL") || inputfile.Contains("TTTo2L2Nu"))
    std::cout << "===Low mass TTTo2L2Nu sample info===="   << std::endl
	      << "Total:     " << ttto2l2nu_input          << std::endl
	      << "Fail GEN:  " << ttto2l2nu_fail_gen_mass  << std::endl
	      << "Fail RECO: " << ttto2l2nu_fail_reco_mass << std::endl;
  if (inputfile.Contains("WWTo2L2Nu"))
    std::cout << "===Low mass WWTo2L2Nu sample info===="   << std::endl
	      << "Total:     " << wwto2l2nu_input          << std::endl
	      << "Fail GEN:  " << wwto2l2nu_fail_gen_mass  << std::endl
	      << "Fail RECO: " << wwto2l2nu_fail_reco_mass << std::endl;

  ///===================================================
  fclose (pFile);
  //==================================================================================
  //                                                                                 =
  //                         writing Histograms to a file                            =
  //                                                                                 =
  //==================================================================================
  output->cd();
  output->Write();
  output->Close();
  output_txt.close();
  //========================================================================
  time (&end);
  dif = difftime (end,start);
  printf ("It took you %.2lf minutes to run your program.\n", (dif/60.0) );
}


//==================================================================================
//                                                                                 =
//                                    Start methods                                =
//                                                                                 =
//==================================================================================

//----------------------------------------------------
//                                                   -
//       Part for Gen & Reco Matching                -
//                                                   -
//----------------------------------------------------
//
//
//

float ZprimeMuMuPatMiniAodNewMC::GetScaleBias(float eta, float phi, float pt, float charge)
{
  double shift = m_muon_scale_ratio_hist->GetBinContent(
							m_muon_scale_ratio_hist->FindBin(eta, phi));
  double uncertainty = m_muon_scale_ratio_hist->GetBinError(
							    m_muon_scale_ratio_hist->FindBin(eta, phi));


  //rand = new TRandom3();
  // Central value correction + gaussian uncertainty
  if (std::abs(eta) < 1.2) {
    shift = 0.0;
    uncertainty = 0.025;
  }
  double ratio   = 1. + (rand->Gaus(shift, uncertainty) * pt / 1000./ charge);
  return 1./ratio;
}

float ZprimeMuMuPatMiniAodNewMC::GetResultionUncert(float pt,float eta)
{
  //rand = new TRandom3();
  double smearing = 0.0;
  if (pt < 200.) {
    smearing = 0.003;
  } else if (pt < 500.) {
    smearing = 0.005;
  } else {
    smearing = 0.01;
  }
  // Double smearing for muons in the endcaps
  if (std::abs(eta) > 1.2) smearing *= 2;
  // Return Gaussian smearing
  double ratio = rand->Gaus(1, smearing);
  return ratio;
}

TLorentzVector ZprimeMuMuPatMiniAodNewMC::GetShiftedMuon(float px, float py, float pz, float E, float ratio)
{
  // TLorentzVector * result = new TLorentzVector();
  std::shared_ptr<TLorentzVector> result = std::make_shared<TLorentzVector>();
  result->SetPxPyPzE(ratio*px,
		     ratio*py,
		     ratio*pz,
		     ratio*E);
  return *result;
}


float ZprimeMuMuPatMiniAodNewMC::delR(float eta1,float phi1,float eta2,float phi2)
{
  float mpi=M_PI;
  float dp=std::abs(phi1-phi2);
  if (dp>mpi) dp-=float(2*mpi);
  return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
}

bool ZprimeMuMuPatMiniAodNewMC::GenRecoMatchMu(float RecoEta1,float RecoPhi1,
					       float &ptGmu,float &etaGmu,float &phiGmu,float &enGmu)
{
  int NbHighPtmu = 0;
  unsigned iflag = -10;
  ptGmu  = -10.;
  etaGmu = -10.;
  phiGmu = -10.;
  enGmu  = -10.;
  for (unsigned i=0; i<iGen->size(); i++) {
    if (fabs(idGen->at(i)) != 13) continue;
    if (statusGen->at(i) != 1)    continue;
    float deltaR1   = delR(RecoEta1,RecoPhi1,etaGen->at(i),phiGen->at(i));
    if (fabs(deltaR1)>deltaRcut)  continue;
    iflag  = i;
    NbHighPtmu++;
    ptGmu  = ptGen->at(i);
    etaGmu = etaGen->at(i);
    phiGmu = phiGen->at(i);
    enGmu  = EnergyGen->at(i);
  }
  if (NbHighPtmu > 0) {
    return true;
  } else {
    return false;
  }
}

//============================ Method to select first high pt muon ========================
bool ZprimeMuMuPatMiniAodNewMC::SelectFirstMuon(float &pTmuon1,float &Enmuon1,float &Etamuon1,
						float &Phimuon1,int &ChargeMu1,unsigned &FlagMu1,
						float &pxmuon1,float &pymuon1,float &pzmuon1,
						float &pmuon1,float &dxymuon1,float &pTmuon1tuneP,
						float &pTmuonBestTrack1,
						float &genMu1Pt, float &genMu1Eta, float &genMu1Phi, float &genMu1En)
{
  int NbHighPtmu = 0;
  unsigned iflag = -10;
  float highestpt=-999.;

  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    //    if (Mu_isMuonsCleaned->at(i) != Mu_isPF->at(i)) continue;
    if (Mu_isTrackerMuon->at(i) == 1 &&
       	Mu_isGlobalMuon->at(i) == 1 &&
	fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
	Mu_ptTunePMuonBestTrack->at(i) > 53.0 &&
	Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
	(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
	Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
	Mu_numberOfValidPixelHits->at(i) > 0 &&
	Mu_numberOfValidMuonHits->at(i) > 0 &&
	Mu_passNewMatchedStationsCut->at(i) == 1 &&
	Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {

      if (Mu_ptTunePMuonBestTrack->at(i)>highestpt) {
	bool GenRecoMatch1 = GenRecoMatchMu(Mu_etaTunePMuonBestTrack->at(i),Mu_phiTunePMuonBestTrack->at(i),
					    genMu1Pt, genMu1Eta, genMu1Phi, genMu1En);
	//if (GenRecoMatch1 == 0) continue;
	highestpt=Mu_ptTunePMuonBestTrack->at(i);
	iflag  = i;
	NbHighPtmu++;
      }
    } else {
      continue;
    }
  }

  if (NbHighPtmu > 0) {
    FlagMu1             = iflag;
    pTmuon1             = Mu_ptTunePMuonBestTrack->at(iflag);
    Enmuon1             = Mu_en->at(iflag);
    Etamuon1            = Mu_etaTunePMuonBestTrack->at(iflag);
    Phimuon1            = Mu_phiTunePMuonBestTrack->at(iflag);
    ChargeMu1           = Mu_chargeTunePMuonBestTrack->at(iflag);
    pxmuon1             = Mu_pxTunePMuonBestTrack->at(iflag);
    pymuon1             = Mu_pyTunePMuonBestTrack->at(iflag);
    pzmuon1             = Mu_pzTunePMuonBestTrack->at(iflag);
    pmuon1              = Mu_pTunePMuonBestTrack->at(iflag);
    dxymuon1            = Mu_absdxyTunePMuonBestTrack->at(iflag);
    pTmuon1tuneP        = Mu_ptTunePMuonBestTrack->at(iflag);
    pTmuonBestTrack1    = Mu_ptTunePMuonBestTrack->at(iflag);
    //std::cout << "First Muon ChargeMu1= " << ChargeMu1 << " Pt1= " << pTmuonBestTrack1 << std::endl;
    return true;
  } else {
    return false;
  }
}

//============================ Method to select second high pt muon ========================
bool ZprimeMuMuPatMiniAodNewMC::SelectSecondMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float Etamuon1,float Phimuon1,
						 float &pTmuon2,float &Enmuon2,
						 float &Etamuon2,float &Phimuon2,int &ChargeMu2,float &pxmuon2,
						 float &pymuon2,float &pzmuon2,float &pmuon2,float &dxymuon2,
						 float &pTmuon2tuneP,float &pTmuonBestTrack2,
						 float &genMu2Pt, float &genMu2Eta, float &genMu2Phi, float &genMu2En)
{
  int NbHighPtmu = 0;
  unsigned iflag = -10;
  float highestpt=-999.;
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (i == FlagMu1) continue;
    if (Mu_ptTunePMuonBestTrack->at(i) == pTmuon1) continue;
    if (Mu_etaTunePMuonBestTrack->at(i) == Etamuon1) continue;
    if (Mu_phiTunePMuonBestTrack->at(i) == Phimuon1) continue;
    //std::cout << "Charge=" << ChargeMu1 << " " << Mu_chargeTunePMuonBestTrack->at(i) << " pT= " <<  Mu_ptTunePMuonBestTrack->at(i) <<  std::endl;
    if (ChargeMu1*Mu_chargeTunePMuonBestTrack->at(i)>0) continue;
    //if (ChargeMu1*Mu_chargeTunePMuonBestTrack->at(i)<0) continue;
    //    if (Mu_isMuonsCleaned->at(i) != Mu_isPF->at(i)) continue;
    if (Mu_isTrackerMuon->at(i) == 1 &&
	Mu_isGlobalMuon->at(i) == 1 &&
	fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
	Mu_ptTunePMuonBestTrack->at(i) > 53.0 &&
	Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
	(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
	Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
	Mu_numberOfValidPixelHits->at(i) > 0 &&
	Mu_numberOfValidMuonHits->at(i) > 0 &&
	Mu_passNewMatchedStationsCut->at(i) == 1 &&
	Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      if (Mu_ptTunePMuonBestTrack->at(i)>highestpt) {
	bool GenRecoMatch2 = GenRecoMatchMu(Mu_etaTunePMuonBestTrack->at(i),Mu_phiTunePMuonBestTrack->at(i),
					    genMu2Pt, genMu2Eta, genMu2Phi, genMu2En);
	//if (GenRecoMatch2 == 0) continue;
	highestpt=Mu_ptTunePMuonBestTrack->at(i);
	//std::cout << "Highest PT second lepton has pt= " << highestpt << std::endl;
	iflag  = i;
	NbHighPtmu++;
      }
    } else {
      continue;
    }
  }
  if (NbHighPtmu > 0 ) {
    pTmuon2          = Mu_ptTunePMuonBestTrack->at(iflag);
    Enmuon2          = Mu_en->at(iflag);
    Etamuon2         = Mu_etaTunePMuonBestTrack->at(iflag);
    Phimuon2         = Mu_phiTunePMuonBestTrack->at(iflag);
    ChargeMu2        = Mu_chargeTunePMuonBestTrack->at(iflag);
    pxmuon2          = Mu_pxTunePMuonBestTrack->at(iflag);
    pymuon2          = Mu_pyTunePMuonBestTrack->at(iflag);
    pzmuon2          = Mu_pzTunePMuonBestTrack->at(iflag);
    pmuon2           = Mu_pTunePMuonBestTrack->at(iflag);
    dxymuon2         = Mu_absdxyTunePMuonBestTrack->at(iflag);
    pTmuon2tuneP     = Mu_ptTunePMuonBestTrack->at(iflag);
    pTmuonBestTrack2 = Mu_ptTunePMuonBestTrack->at(iflag);
    //std::cout << "Second ChargeMu2= " << ChargeMu2 << " Pt2= " << pTmuonBestTrack2 << std::endl;
    return true;
  } else {
    return false;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotRecoInfo(float CosmicMuonRejec, float vertexMassMu,float MassGenerated,
					     float PtTunePMuBestTrack,float PtTunePMu,float PtMuBestTrack,
					     float PtGenerated, float etaMu1, float pMu1,
					     float PtTunePMuBestTrack2,float PtTunePMu2,float PtMuBestTrack2,
					     float PtGenerated2,float etaMu2, float pMu2, float bosonPt, TString name)
{
  //----------------------------------------------------------
  if (vertexMassMu>900.0) {
    output_txt << event_runNo
	       << "   "      << event_lumi
	       << "       "  << event_evtNo
	       << "        " << vertexMassMu
	       << "        " << PtTunePMuBestTrack
	       << "        " << PtTunePMuBestTrack2
	       << "        " << etaMu1
	       << "        " << etaMu2
	       << std::endl;

  }

  // only for DY POWHEG??
  float weight10 = 1.;
  if (name.Contains("NNPDF30"))
    weight10 = MassCorrection(MassGenerated, bosonPt, etaMu1, etaMu2);

  newweight = newweight*weight10;
  m_recoMassCorr = vertexMassMu*weight10;

  //float weight2 = std::min(1.01696-7.73522E-5*vertexMassMu+6.69239E-9*vertexMassMu*vertexMassMu,1);
  //----------------------------------------------------------
  h1_ZprimeRecomassBinWidth_->Fill(vertexMassMu,newweight);
  h1_ZprimeRecomasslogscale_->Fill(log10(vertexMassMu),newweight);
  h1_ZprimeRecomass_->Fill(vertexMassMu,newweight);
  h1_MassRecoInAccep_->Fill(MassGenerated,newweight);

  float SF1 = 1.;
  float SF2 = 1.;

  if (fabs(etaMu1) <= 1.6 && pMu1 > 100) SF1 = (0.994 - 4.08e-6 * pMu1)/(0.994 - 4.08e-6 * 100);
  else if (fabs(etaMu1) > 1.6 && pMu1 > 200)  SF1 = ((0.9784 - 4.73e-5 * pMu1)/(0.9908 - 1.26e-5 * pMu1)) / ((0.9784 - 4.73e-5 * 200)/(0.9908 - 1.26e-5 * 200)) ;
  if (fabs(etaMu2) <= 1.6 && pMu2 > 100) SF2 = (0.994 - 4.08e-6 * pMu2)/(0.994 - 4.08e-6 * 100);
  else if (fabs(etaMu2) > 1.6 && pMu2 > 200)  SF2 = ((0.9784 - 4.73e-5 * pMu2)/(0.9908 - 1.26e-5 * pMu2)) / ((0.9784 - 4.73e-5 * 200)/(0.9908 - 1.26e-5 * 200) ) ;




  h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        0.,newweight);
  h2_CSMassBinned_       ->Fill(m_vtxMassMu,               0.,newweight);
  h2_CSMassMuIDBinned_   ->Fill(m_vtxMassMu,               0.,newweight*SF1*SF2);
  h2_CSMassUpBinned_     ->Fill(m_vtxMassScaledMu,         0.,newweight);
  h2_CSMassDownBinned_   ->Fill(m_vtxMassScaledMu,         0.,newweight);
  if (m_csAngle > 0) {
    h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        2.,newweight);
    h2_CSMassBinned_       ->Fill(m_vtxMassMu,               2.,newweight);
    h2_CSMassMuIDBinned_   ->Fill(m_vtxMassMu,               2.,newweight*SF1*SF2);
    h2_CSMassUpBinned_     ->Fill(m_vtxMassScaledMu,               2.,newweight);
    h2_CSMassDownBinned_   ->Fill(m_vtxMassScaledMu,               2.,newweight);
  } else {
    h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        1.,newweight);
    h2_CSMassBinned_       ->Fill(m_vtxMassMu,               1.,newweight);
    h2_CSMassMuIDBinned_   ->Fill(m_vtxMassMu,               1.,newweight*SF1*SF2);
    h2_CSMassUpBinned_     ->Fill(m_vtxMassScaledMu,         1.,newweight);
    h2_CSMassDownBinned_   ->Fill(m_vtxMassScaledMu,         1.,newweight);
  }

  /*
  std::cout << h1_ZprimeRecomass_->GetEntries()
	    << ", " << h2_CSMassBinned_->GetEntries();
  for (int i = 0; i < 25; ++i)
    std::cout << ", " << h2_CSMassBinned_->ProjectionX("tmp",i,i)->GetEntries();
  std::cout << std::endl;
  if (fabs(h2_CSMassBinned_->ProjectionX("tmp",1,1)->GetEntries()-h1_ZprimeRecomass_->GetEntries()) > 0.01)
    std::cout << "Mass values different"
	      << "  h2_CSMassBinned_->ProjectionX(\"tmp\",1,1)->GetEntries(): "
	      << h2_CSMassBinned_->ProjectionX("tmp",1,1)->GetEntries()
	      << "  h1_ZprimeRecomass_->GetEntries(): " << h1_ZprimeRecomass_->GetEntries()
	      << std::endl;
  if (fabs(m_weight-newweight) > 0)
    std::cout << "Weights different"
	      << "  weight: " << m_weight
	      << "  newweight: " << newweight
	      << std::endl;
  if (fabs(m_vtxMassMu-vertexMassMu) > 0)
    std::cout << "Mass values different"
	      << "  m_vtxMassMu: " << m_vtxMassMu
	      << "  vertexMassMu: " << vertexMassMu
	      << std::endl;
  */

  int priEtaBin = -1;
  int secEtaBin = -1;

  if (fabs(etaMu1) < 1.2 && fabs(etaMu2) < 1.2) {  //BB
    h1_ZprimeRecomassBinWidthBB_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomass60to120BB_->Fill(vertexMassMu,newweight);
    // SHOULD BE FILLED FOR ALL EVENTS???
    h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
    priEtaBin = 1;
  } else if ((fabs(etaMu1) < 1.2 && (fabs(etaMu2) > 1.2 && fabs(etaMu2) < 2.4)) ||
	     (fabs(etaMu2) < 1.2 && (fabs(etaMu1) > 1.2 && fabs(etaMu1) < 2.4))) {  //BE
    h1_ZprimeRecomassBinWidthAllBE_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
    priEtaBin = 2;

    //==================Sub-categories============================================
    if ((fabs(etaMu1) < 1.2 && (etaMu2 > 1.2 && etaMu2 < 2.4)) ||
	(fabs(etaMu2) < 1.2 && (etaMu1 > 1.2 && etaMu1 < 2.4))) {  //BE+
      h1_ZprimeRecomassBinWidthBEpos_->Fill(vertexMassMu,newweight);
      h1_ZprimeRecomass60to120BEpos_->Fill(vertexMassMu,newweight);
      // h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
      // h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
      secEtaBin = 4;
    } else if ((fabs(etaMu1) < 1.2 && (etaMu2 > -2.4 && etaMu2 < -1.2)) ||
	       (fabs(etaMu2) < 1.2 && (etaMu1 > -2.4 && etaMu1 < -1.2))) {  //BE-
      h1_ZprimeRecomassBinWidthBEnev_->Fill(vertexMassMu,newweight);
      h1_ZprimeRecomass60to120BEnev_->Fill(vertexMassMu,newweight);
      // WHY FILLED FOR BE+ BUT NOT BE-??
      // h1_ZprimeRecomass_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
      // h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
      // h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
      secEtaBin = 5;
    }
  } else if ((fabs(etaMu1) > 1.2 && fabs(etaMu1) < 2.4) &&
	     (fabs(etaMu2) > 1.2 && fabs(etaMu2) < 2.4)) {  //EE
    h1_ZprimeRecomassBinWidthAllEE_->Fill(vertexMassMu,newweight);
    priEtaBin = 3;
  //==================Sub-categories============================================
    if (((etaMu1 > 1.2 && etaMu1 < 2.4) && (etaMu2 > -2.4 && etaMu2 < -1.2)) ||
	((etaMu2 > 1.2 && etaMu2 < 2.4) && (etaMu1 > -2.4 && etaMu1 < -1.2))) {  //E+E-
      h1_ZprimeRecomassBinWidthEE_->Fill(vertexMassMu,newweight);
      h1_ZprimeRecomass60to120EE_->Fill(vertexMassMu,newweight);
      //h1_ZprimeRecomass_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
      //h1_ZprimeRecomassBinWidth_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
      h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
      secEtaBin = 6;
    // currently ignored cases?
    } else if ((etaMu1 > -2.4 && etaMu1 < -1.2) && (etaMu2 > -2.4 && etaMu2 < -1.2)) {  //E-E-
      secEtaBin = 7;
    } else if ((etaMu1 > 1.2 && etaMu1 < 2.4) && (etaMu2 > 1.2 && etaMu2 < 2.4)) {  //E+E+
      secEtaBin = 8;
    }
  }

  /*
  //==================Sub-categories============================================
  //BE+
  if (fabs(etaMu1) < 1.2 && (etaMu2 > 1.2 && etaMu2 < 2.4)) {
    h1_ZprimeRecomassBinWidthBEpos_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomass60to120BEpos_->Fill(vertexMassMu,newweight);
    //h1_ZprimeRecomass_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
    h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
    secEtaBin = 4;
  }

  if (fabs(etaMu2) < 1.2 && (etaMu1 > 1.2 && etaMu1 < 2.4)) {
    h1_ZprimeRecomassBinWidthBEpos_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomassBinWidthBEpos_->Fill(vertexMassMu,newweight);
    //h1_ZprimeRecomass_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
    h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
    secEtaBin = 4;
    }

  //BE-
  if (fabs(etaMu1) < 1.2 && (etaMu2 > -2.4 && etaMu2 < -1.2)) {
    h1_ZprimeRecomassBinWidthBEnev_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomass60to120BEnev_->Fill(vertexMassMu,newweight);
    //h1_ZprimeRecomass_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
    //h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
    secEtaBin = 5;
  }

  if (fabs(etaMu2) < 1.2 && (etaMu1 > -2.4 && etaMu1 < -1.2)) {
    h1_ZprimeRecomassBinWidthBEnev_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomass60to120BEnev_->Fill(vertexMassMu,newweight);
    //h1_ZprimeRecomass_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
    //h1_ZprimeRecomassBinWidth_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
    h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
    secEtaBin = 5;
    }

  //E+E-
  if ((etaMu1 > 1.2 && etaMu1 < 2.4) && (etaMu2 > -2.4 && etaMu2 < -1.2)) {
    h1_ZprimeRecomassBinWidthEE_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomass60to120EE_->Fill(vertexMassMu,newweight);
    //h1_ZprimeRecomass_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
    //h1_ZprimeRecomassBinWidth_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
    h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
    secEtaBin = 6;
  }
  if ((etaMu2 > 1.2 && etaMu2 < 2.4) && (etaMu1 > -2.4 && etaMu1 < -1.2)) {
    h1_ZprimeRecomassBinWidthEE_->Fill(vertexMassMu,newweight);
    h1_ZprimeRecomass60to120EE_->Fill(vertexMassMu,newweight);
    //h1_ZprimeRecomass_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
    //h1_ZprimeRecomassBinWidth_->Fill(vertexMassMu,newweight); // DUPLICATE!!!!!
    h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
    secEtaBin = 6;
  }
  // currently ignored cases
  //E-E-
  if ((etaMu1 > -2.4 && etaMu1 < -1.2) &&
      (etaMu2 > -2.4 && etaMu2 < -1.2)) {
    secEtaBin = 7;
  }
  //E+E+
  if ((etaMu1 < 2.4 && etaMu1 > 1.2) &&
      (etaMu2 < 2.4 && etaMu2 > 1.2)) {
    secEtaBin = 8;
  }
  */

  /*
  std::cout << "priEtaBin=" << priEtaBin
	    << ", secEtaBin=" << secEtaBin
	    << ", etaMu1=" << etaMu1
	    << ", etaMu2=" << etaMu2
	    << std::endl;
  std::cout << "primary bin (" << priEtaBin << "*3)+0(" << (priEtaBin*3)+0 << ")" << std::endl;
  */
  h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (priEtaBin*3)+0,newweight);
  h2_CSMassBinned_       ->Fill(m_vtxMassMu,               (priEtaBin*3)+0,newweight);
  h2_CSMassMuIDBinned_   ->Fill(m_vtxMassMu,               (priEtaBin*3)+0,newweight*SF1*SF2);
  h2_CSMassUpBinned_     ->Fill(m_vtxMassScaledMu,         (priEtaBin*3)+0,newweight);
  h2_CSMassDownBinned_   ->Fill(m_vtxMassScaledMu          ,(priEtaBin*3)+0,newweight);
  if (secEtaBin > 0) {
    // std::cout << "secondary bin (" << secEtaBin << "*3)+0(" << (secEtaBin*3)+0 << ")" << std::endl;
    h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (secEtaBin*3)+0,newweight);
    h2_CSMassBinned_       ->Fill(m_vtxMassMu,               (secEtaBin*3)+0,newweight);
    h2_CSMassMuIDBinned_   ->Fill(m_vtxMassMu,               (secEtaBin*3)+0,newweight*SF1*SF2);
    h2_CSMassUpBinned_     ->Fill(m_vtxMassScaledMu,         (secEtaBin*3)+0,newweight);
    h2_CSMassDownBinned_   ->Fill(m_vtxMassScaledMu,         (secEtaBin*3)+0,newweight);
  }
  if (m_csAngle > 0) {
    // std::cout << "primary bin (" << priEtaBin << "*3)+2(" << (priEtaBin*3)+2 << ")" << std::endl;
    h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (priEtaBin*3)+2,newweight);
    h2_CSMassBinned_       ->Fill(m_vtxMassMu,               (priEtaBin*3)+2,newweight);
    h2_CSMassMuIDBinned_   ->Fill(m_vtxMassMu,               (priEtaBin*3)+2,newweight*SF1*SF2);
    h2_CSMassUpBinned_     ->Fill(m_vtxMassScaledMu,         (priEtaBin*3)+2,newweight);
    h2_CSMassDownBinned_   ->Fill(m_vtxMassScaledMu,         (priEtaBin*3)+2,newweight);
    if (secEtaBin > 0) {
      // std::cout << "secondary bin (" << secEtaBin << "*3)+2(" << (secEtaBin*3)+2 << ")" << std::endl;
      h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (secEtaBin*3)+2,newweight);
      h2_CSMassBinned_       ->Fill(m_vtxMassMu,               (secEtaBin*3)+2,newweight);
      h2_CSMassMuIDBinned_   ->Fill(m_vtxMassMu,               (secEtaBin*3)+2,newweight*SF1*SF2);
      h2_CSMassUpBinned_     ->Fill(m_vtxMassScaledMu,         (secEtaBin*3)+2,newweight);
      h2_CSMassDownBinned_   ->Fill(m_vtxMassScaledMu,         (secEtaBin*3)+2,newweight);
    }
  } else {
    // std::cout << "primary bin (" << priEtaBin << "*3)+1(" << (priEtaBin*3)+1 << ")" << std::endl;
    h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (priEtaBin*3)+1,newweight);
    h2_CSMassBinned_       ->Fill(m_vtxMassMu,               (priEtaBin*3)+1,newweight);
    h2_CSMassMuIDBinned_   ->Fill(m_vtxMassMu,               (priEtaBin*3)+1,newweight*SF1*SF2);
    h2_CSMassUpBinned_     ->Fill(m_vtxMassScaledMu,         (priEtaBin*3)+1,newweight);
    h2_CSMassDownBinned_   ->Fill(m_vtxMassScaledMu,         (priEtaBin*3)+1,newweight);
    if (secEtaBin > 0) {
      // std::cout << "secondary bin (" << secEtaBin << "*3)+1(" << (secEtaBin*3)+1 << ")" << std::endl;
      h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (secEtaBin*3)+1,newweight);
      h2_CSMassBinned_       ->Fill(m_vtxMassMu,               (secEtaBin*3)+1,newweight);
      h2_CSMassMuIDBinned_   ->Fill(m_vtxMassMu,               (secEtaBin*3)+1,newweight*SF1*SF2);
      h2_CSMassUpBinned_     ->Fill(m_vtxMassScaledMu,         (secEtaBin*3)+1,newweight);
      h2_CSMassDownBinned_   ->Fill(m_vtxMassScaledMu,         (secEtaBin*3)+1,newweight);
    }
  }

  h1_ZprimeRecomass50_->Fill(vertexMassMu,newweight);
  h1_ZprimeRecomass20_->Fill(vertexMassMu,newweight);
  if (fabs(etaMu1) < 1.2 && fabs(etaMu2) < 1.2) {
    h1_ZprimeRecomassBB_->Fill(vertexMassMu,newweight);
  } else if (fabs(etaMu1) > 1.2 && fabs(etaMu2) > 1.2) {
    h1_ZprimeRecomassEE_->Fill(vertexMassMu,newweight);
  } else if ((fabs(etaMu1) < 1.2 && fabs(etaMu2) > 1.2) ||
	     (fabs(etaMu1) > 1.2 && fabs(etaMu2) < 1.2)) {
    h1_ZprimeRecomassBE_->Fill(vertexMassMu,newweight);
  }

  /*if (vertexMassMu>60 && vertexMassMu<120)
    h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);*/
  h1_3Dangle_->Fill(CosmicMuonRejec,newweight);

  // part for Pt resolution
  h1_PtResolutionTunePMBT_->Fill((PtTunePMuBestTrack-PtGenerated)/PtGenerated,newweight);
  h1_PtResolutiontuneP_->Fill((PtTunePMu-PtGenerated)/PtGenerated,newweight);
  h1_PtResolutionMBT_->Fill((PtMuBestTrack-PtGenerated)/PtGenerated,newweight);

  //part for mass resolution
  if (vertexMassMu > 0.0 && vertexMassMu < 250.0 ) {
    h1_MassResultionEBEB1_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);
  } else if (vertexMassMu > 250 && vertexMassMu < 750.0 ) {
    h1_MassResultionEBEB2_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);
  } else if (vertexMassMu > 750 && vertexMassMu < 1250.0 ) {
    h1_MassResultionEBEB3_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);
  } else if (vertexMassMu > 1250 && vertexMassMu < 1750.0 ) {
    h1_MassResultionEBEB4_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);
  } else if (vertexMassMu > 1750 && vertexMassMu < 2250.0 ) {
    h1_MassResultionEBEB5_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);
  } else if (vertexMassMu > 2000 && vertexMassMu < 4000.0 ) {
    h1_MassResultionEBEB6_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);
  } else if (vertexMassMu > 4000 && vertexMassMu < 6000.0) {
    h1_MassResultionEBEB7_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);
  }
}

//===================== Method to calculate the mass ========================
float ZprimeMuMuPatMiniAodNewMC::Mass(float Pt1,float Eta1,float Phi1,float En1,
				      float Pt2,float Eta2,float Phi2,float En2)
{
  float MuMuMass = 0.0;
  TLorentzVector Mu1;
  TLorentzVector Mu2;
  Mu1.SetPtEtaPhiM(Pt1,Eta1,Phi1,En1);
  Mu2.SetPtEtaPhiM(Pt2,Eta2,Phi2,En2);
  MuMuMass = sqrt((Mu1+Mu2).M2());
  return MuMuMass;
}

//===================== Method to calculate the smeared mass ========================
float ZprimeMuMuPatMiniAodNewMC::smearedMass(float Eta1, float Phi1, float Pt1, float Eta2, float Phi2, float Pt2,
					     float vtxMass)
{

  double ratio1 = GetResultionUncert(Pt1,Eta1);
  double ratio2 = GetResultionUncert(Pt2,Eta2);

  // TLorentzVector * muon1 = new TLorentzVector();
  std::shared_ptr<TLorentzVector> muon1 = std::make_shared<TLorentzVector>();
  muon1->SetPtEtaPhiM(Pt1,Eta1,Phi1,MUON_MASS);
  TLorentzVector shiftedMuon1 = GetShiftedMuon(muon1->Px(),muon1->Py(),muon1->Pz(),muon1->E(), ratio1);

  // TLorentzVector * muon2 = new TLorentzVector();
  std::shared_ptr<TLorentzVector> muon2 = std::make_shared<TLorentzVector>();
  muon2->SetPtEtaPhiM(Pt2,Eta2,Phi2,MUON_MASS);
  TLorentzVector shiftedMuon2 = GetShiftedMuon(muon2->Px(),muon2->Py(),muon2->Pz(),muon2->E(), ratio2);

  double mass = (*muon1 + *muon2).M();
  double massShifted = (shiftedMuon1 + shiftedMuon2).M();

  return vtxMass*massShifted/mass;

}
float ZprimeMuMuPatMiniAodNewMC::scaledMass(float Eta1, float Phi1, float Pt1, float Charge1, float Eta2, float Phi2, float Pt2, float Charge2,
					     float vtxMass)
{

  double ratio1 = GetScaleBias(Eta1,Phi1, Pt1, Charge1);
  double ratio2 = GetScaleBias(Eta2,Phi2, Pt2, Charge2);

  // TLorentzVector * muon1 = new TLorentzVector();
  std::shared_ptr<TLorentzVector> muon1 = std::make_shared<TLorentzVector>();
  muon1->SetPtEtaPhiM(Pt1,Eta1,Phi1,MUON_MASS);
  TLorentzVector shiftedMuon1 = GetShiftedMuon(muon1->Px(),muon1->Py(),muon1->Pz(),muon1->E(), ratio1);

  // TLorentzVector * muon2 = new TLorentzVector();
  std::shared_ptr<TLorentzVector> muon2 = std::make_shared<TLorentzVector>();
  muon2->SetPtEtaPhiM(Pt2,Eta2,Phi2,MUON_MASS);
  TLorentzVector shiftedMuon2 = GetShiftedMuon(muon2->Px(),muon2->Py(),muon2->Pz(),muon2->E(), ratio2);

  double mass = (*muon1 + *muon2).M();
  double massShifted = (shiftedMuon1 + shiftedMuon2).M();

  return vtxMass*massShifted/mass;

}

void ZprimeMuMuPatMiniAodNewMC::PickThehighestMass(float &vtxHighestMass,float &vtxHighestChi2,int EvtNb)
{
  float Massinv  = -10.0;
  unsigned iflag = -10;
  vtxHighestMass = -10.0;
  vtxHighestChi2 = 1000.0;
  int NbMu = 0;
  int Nb = 0;
  int countlept=0;
  for (unsigned i=0; i<Mu_vtxMass->size(); i++) {
    Nb++;
    countlept=2*i;
    //std::cout << "vtx mass" << Mu_vtxMass->at(i) << " Chi2= " << Mu_vtxNormChi2->at(i)<< std::endl;
    //std::cout << "vtx Mass lepton= " << Mu_vtxMassLept->at(countlept) << " " <<  Mu_vtxMassLept->at(countlept+1)<< std::endl;
    float chargepair=0;
    for (unsigned j=0; j<Mu_nbMuon->size(); j++) {
      if (Mu_ptTunePMuonBestTrack->at(j)==Mu_vtxMassLept->at(countlept)) chargepair=Mu_chargeTunePMuonBestTrack->at(j);
      if (Mu_ptTunePMuonBestTrack->at(j)==Mu_vtxMassLept->at(countlept+1)) chargepair=chargepair*Mu_chargeTunePMuonBestTrack->at(j);
    }

    //std::cout << "Chargepair for vtxmass= " <<  Mu_vtxMass->at(i) << "   is " << chargepair << std::endl;
    if (chargepair!=-1) continue;
    //std::cout << "Chargepair for vtxmass= " <<  Mu_vtxMass->at(i) << "   is " << chargepair << std::endl;


    int leptmatchBest=0;
    if (Mu_vtxMassLept->at(countlept)==PtRecTunePMu1 && Mu_vtxMassLept->at(countlept+1)==PtRecTunePMu2)
      leptmatchBest=2;
    if (Mu_vtxMassLept->at(countlept)==PtRecTunePMu2 && Mu_vtxMassLept->at(countlept+1)==PtRecTunePMu1)
      leptmatchBest=2;

    if (leptmatchBest!=2)
      continue;
    //std::cout << "Chargepair surviving matching for vtxmass= " <<  Mu_vtxMass->at(i) << "   with " << leptmatchBest << " leptons" << std::endl;

    if (Mu_vtxNormChi2->at(i)> 20)
      continue;
    if (Mu_vtxMass->at(i)>Massinv) {
      Massinv = Mu_vtxMass->at(i);
      iflag  = i;
      NbMu++;
    }
  }
  if ( NbMu > 0 ) {
    vtxHighestMass = Mu_vtxMass->at(iflag);
    vtxHighestChi2 = Mu_vtxNormChi2->at(iflag);
  }
}

double ZprimeMuMuPatMiniAodNewMC::ThreeDangle(float pxMu1,float pyMu1,float pzMu1,float pMu1,
					      float pxMu2,float pyMu2,float pzMu2,float pMu2)
{
  TVector3 Mu1(pxMu1,pyMu1,pzMu1);
  TVector3 Mu2(pxMu2,pyMu2,pzMu2);
  double cos_angle = Mu1.Dot(Mu2) / pMu1 / pMu2;
  return cos_angle;
}

//----------------------------------------------------
//                                                   -
//       Part for Gen & Reco Matching                -
//                                                   -
//----------------------------------------------------
//========================== Method to select firt Gen Mu =======================
bool ZprimeMuMuPatMiniAodNewMC::SelectFirstGenMu(float &ETMu1,float &PhiMu1,
						 float &EtaMu1,float &EnMu1,
						 int &IDmu1,int &Statele1,
						 unsigned &GenFlag1)
{
  int NbHighPtmu = 0;
  int iflag = -10;
  ETMu1 = 0.0;
  for (unsigned i=0; i<iGen->size(); i++) {
    if (fabs(idGen->at(i)) != 13 ) continue;
    if (statusGen->at(i) != 1 )  continue;
    if (ptGen->at(i) > ETMu1) {
      ETMu1 = ptGen->at(i);
      iflag = i;
      NbHighPtmu++;
    } else {
      continue;
    }
  }
  if (NbHighPtmu>0) {
    GenFlag1 = iflag;
    ETMu1    = ptGen->at(iflag);
    PhiMu1   = phiGen->at(iflag);
    EtaMu1   = etaGen->at(iflag);
    EnMu1    = EnergyGen->at(iflag);
    IDmu1    = idGen->at(iflag);
    Statele1 = statusGen->at(iflag);
    return true;
  } else {
    return false;
  }
}
//============================ Method to select second Gen Mu ========================
bool ZprimeMuMuPatMiniAodNewMC::SelectSecondGenMu(unsigned GenFlag1, float ETMu1,
						  float &ETMu2, float &PhiMu2, float &EtaMu2,
						  float &EnMu2,int &IDmu2,int &Statele2)
{
  int NbHighPtmu = 0;
  int iflag = -10;
  ETMu2 = 0.0;
  for (unsigned i=0; i<iGen->size(); i++) {
    if (fabs(idGen->at(i)) != 13 ) continue;
    if (statusGen->at(i) != 1 )  continue;
    if (i == GenFlag1) continue;
    if (fabs(ptGen->at(i) - ETMu1) <0.00001 ) continue;
    if (ptGen->at(i) > ETMu2) {
      ETMu2 = ptGen->at(i);
      iflag = i;
      NbHighPtmu++;
    } else {
      continue;
    }
  }
  if (NbHighPtmu>0) {
    ETMu2    = ptGen->at(iflag);
    PhiMu2   = phiGen->at(iflag);
    EtaMu2   = etaGen->at(iflag);
    EnMu2    = EnergyGen->at(iflag);
    IDmu2    = idGen->at(iflag);
    Statele2 = statusGen->at(iflag);
    return true;
  } else {
    return false;
  }
}

//============================ Method to compute gen level invariant mass ========================
float ZprimeMuMuPatMiniAodNewMC::GenMass(float ETMu1, float PhiMu1, float EtaMu1,float EnMu1,
					 float ETMu2, float PhiMu2, float EtaMu2,float EnMu2)
{
  TLorentzVector mu1, mu2;
  // mu1.SetPtEtaPhiE(ETMu1,EtaMu1,PhiMu1,EnMu1);
  // mu2.SetPtEtaPhiE(ETMu2,EtaMu2,PhiMu2,EnMu2);
  mu1.SetPtEtaPhiM(ETMu1,EtaMu1,PhiMu1,MUON_MASS);
  mu2.SetPtEtaPhiM(ETMu2,EtaMu2,PhiMu2,MUON_MASS);
  return (mu1+mu2).M();
}

//============================ Method to plot Gen Mu ========================
void ZprimeMuMuPatMiniAodNewMC::PlotGenInfo(float ZprimeGenMass,float EtaGenMu1,float EtaGenMu2,float PtGenMu1,
					    float PtGenMu2,float EnGenMu1,float EnGenMu2)
{
  h1_MassGenInAccep_->Fill(ZprimeGenMass,newweight);
  h1_ZprimeGenmass_->Fill(ZprimeGenMass,newweight);
  h1_ZprimeGenEta1_->Fill(EtaGenMu1,newweight);
  h1_ZprimeGenEta2_->Fill(EtaGenMu2,newweight);
  h1_ZprimeGenPt1_->Fill(PtGenMu1,newweight);
  h1_ZprimeGenPt2_->Fill(PtGenMu2,newweight);
  h1_ZprimeGenEn1_->Fill(EnGenMu1,newweight);
  h1_ZprimeGenEn2_->Fill(EnGenMu2,newweight);
}
//----------------------------------------------------
//                                                   -
//                           N-1 eff.                -
//                                                   -
//----------------------------------------------------
void ZprimeMuMuPatMiniAodNewMC::MuonPassingID()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      h1_PtID_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaID_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
      h1_nbPVID_->Fill(Mu_nbofpv->at(i),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotPterror()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1) {
      h1_dPToverPT_->Fill(Mu_dPToverPTTunePMuonBestTrack->at(i),newweight );
      h1_PtEffpterror_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffpterror_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotNbTrackLayers()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       //Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      h1_numberOftrackerLayersWithMeasurement_->Fill( Mu_numberOftrackerLayersWithMeasurement->at(i),newweight );
      h1_PtEffptnumberOftrackerLayers_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptnumberOftrackerLayers_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotNBValidPixelHits()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       //Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      h1_numberOfValidPixelHits_->Fill( Mu_numberOfValidPixelHits->at(i),newweight );
      h1_PtEffptnumberOfPixelHits_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptnumberOfPixelHits_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotNbValidMuonHits()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       //Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      h1_numberOfValidMuonHits_->Fill( Mu_numberOfValidMuonHits->at(i),newweight );
      h1_PtEffptnumberOfMuonHits_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptnumberOfMuonHits_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotNbMatchedStations()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       //Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      h1_numberOfMatchedStations_->Fill( Mu_numberOfMatchedStations->at(i),newweight );
      h1_PtEffptnumberOfMatchedStations_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptnumberOfMatchedStations_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotTrackiso()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       //(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      h1_trackiso_->Fill( Mu_trackiso->at(i)/Mu_ptTunePMuonBestTrack->at(i),newweight );
      h1_PtEffptTrackIso_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptTrackIso_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotAbsDxy()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       //Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      h1_absdxy_->Fill( Mu_absdxyTunePMuonBestTrack->at(i) ,newweight);
      h1_PtEffptabsdsy_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptabsdsy_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotPtTuneP()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      h1_PtTuneP_->Fill( Mu_ptTunePMuonBestTrack->at(i) ,newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::plotAllHighPtMuonsID()
{
  MuonPassingID();
  PlotPterror();
  PlotNbTrackLayers();
  PlotNBValidPixelHits();
  PlotNbValidMuonHits();
  PlotNbMatchedStations();
  PlotTrackiso();
  PlotAbsDxy();
  MuonPassingNewID();
  MuonPassingTightID();
  PlotPtTuneP();
}

void ZprimeMuMuPatMiniAodNewMC::MuonPassingNewID()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.02 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      h1_PtNewID_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaNewID_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
      h1_nbPVNewID_->Fill(Mu_nbofpv->at(i),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::MuonPassingTightID()
{
  for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
    if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.01 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ) {
      h1_PtTightID_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaTightID_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
      h1_nbPVTightID_->Fill(Mu_nbofpv->at(i),newweight);
    }
    else continue;
  }
}

float ZprimeMuMuPatMiniAodNewMC::CosThetaCollinSoper(float Et1,float Eta1,float Phi1,float En1,
						     float Et2,float Eta2,float Phi2,float En2,
						     float ChargeEle1, float RecoMass) {

  TLorentzVector Ele;
  TLorentzVector Elebar;
  if (ChargeEle1<0) {
    Ele.SetPtEtaPhiE(Et1,Eta1,Phi1,En1);
    Elebar.SetPtEtaPhiE(Et2,Eta2,Phi2,En2);
  }
  if (ChargeEle1>0) {
    Ele.SetPtEtaPhiE(Et2,Eta2,Phi2,En2);
    Elebar.SetPtEtaPhiE(Et1,Eta1,Phi1,En1);
  }
  TLorentzVector Q(Ele+Elebar);

  //************************************************************************
  //
  // 1) cos(theta) = 2 Q^-1 (Q^2+Qt^2)^-1 (Ele^+ Elebar^- - Ele^- Elebar^+)
  //
  //
  //************************************************************************
  double Eleplus  = 1.0/sqrt(2.0) * (Ele.E() + Ele.Z());
  double Eleminus = 1.0/sqrt(2.0) * (Ele.E() - Ele.Z());

  double Elebarplus  = 1.0/sqrt(2.0) * (Elebar.E() + Elebar.Z());
  double Elebarminus = 1.0/sqrt(2.0) * (Elebar.E() - Elebar.Z());

  double costheta = 2.0 / Q.Mag() / sqrt(pow(Q.Mag(), 2) + pow(Q.Pt(), 2)) *
    (Eleplus * Elebarminus - Eleminus * Elebarplus);

  if (Q.Pz()<0.0) costheta = -costheta;

  if (RecoMass > 60.0 && RecoMass < 120.0) {
    h1_CosAngleCollinSoperCorrect60Mass120_->Fill(costheta,newweight);
  }

  if (RecoMass > 120.0 && RecoMass < 300.0) {
    h1_CosAngleCollinSoperCorrect120Mass300_->Fill(costheta,newweight);
  }

  if (RecoMass > 300.0 && RecoMass < 700.0) {
    h1_CosAngleCollinSoperCorrect300Mass700_->Fill(costheta,newweight);
  }

  if (RecoMass > 700.0 && RecoMass < 3000.0) {
    h1_CosAngleCollinSoperCorrect700Mass3000_->Fill(costheta,newweight);
  }

  if (RecoMass > 4500.0 && RecoMass < 6000.0) {
    h1_CosAngleCollinSoperCorrect4900Mass5100_->Fill(costheta,newweight);
    h1_absCosAngleCollinSoperCorrect4500Mass5500_->Fill(fabs(costheta),newweight);
  }

  /************************************************************************
   *
   * 2) tanphi = (Q^2 + Qt^2)^1/2 / Q (Dt dot R unit) /(Dt dot Qt unit)
   *
   ************************************************************************/
  TLorentzVector Pbeam(0.0, 0.0,  4000., 4000.); // beam momentum in lab frame
  TLorentzVector Ptarget(0.0, 0.0, -4000., 4000.); // beam momentum in lab frame
  TLorentzVector D(Ele-Elebar);
  // unit vector on R direction
  TVector3 R = Pbeam.Vect().Cross(Q.Vect());
  TVector3 Runit = R.Unit();
  // unit vector on Qt
  TVector3 Qt = Q.Vect(); Qt.SetZ(0);
  TVector3 Qtunit = Qt.Unit();
  TVector3 Dt = D.Vect(); Dt.SetZ(0);
  double tanphi = sqrt(pow(Q.Mag(), 2) + pow(Q.Perp(), 2)) / Q.Mag() * Dt.Dot(Runit) / Dt.Dot(Qtunit);
  if (Q.Pz()<0.0) tanphi = -tanphi;
  //h1_TanPhiCollinSoperCorrect_->Fill(tanphi,newweight);
  /************************************************************************
   *
   * 3) sin2(theta) = Q^-2 Dt^2 - Q^-2 (Q^2 + Qt^2)^-1 * (Dt dot Qt)^2
   *
   ************************************************************************/
  //double dt_qt = D.X()*Q.X() + D.Y()*Q.Y();
  //double sin2theta = pow(D.Pt()/Q.Mag(), 2)
  //- 1.0/pow(Q.Mag(), 2)/(pow(Q.Mag(), 2) + pow(Q.Pt(), 2))*pow(dt_qt, 2);
  //h1_Sin2AngleCollinSoperCorrect_->Fill(sin2theta,newweight);

  return costheta;
}


void ZprimeMuMuPatMiniAodNewMC::PrintEventInformation(unsigned int runNumber, unsigned int lumiNumber, unsigned int eventNumber,
						      float vtxChi2, float vtxMass, float CosmicRejection)
{
  if (event_runNo == runNumber && event_lumi == lumiNumber && event_evtNo == eventNumber) {
    output_txt << event_runNo
	       << "        " << event_lumi
	       << "        " << event_evtNo
	       << "        " << vtxChi2
	       << "        " << vtxMass << std::endl;
    for (unsigned i=0; i<Mu_nbMuon->size(); i++) {
      /*if (fabs(Mu_etaTunePMuonBestTrack->at(i)) < EtaCut ) {
	std::cout<<"[1] eta="<<Mu_etaTunePMuonBestTrack->at(i) << std::endl;
	}*/
      std::cout<<"[0] phi="<<Mu_phiTunePMuonBestTrack->at(i) << std::endl;
      std::cout<<"[1] eta="<<Mu_etaTunePMuonBestTrack->at(i) << std::endl;
      if (Mu_isGlobalMuon->at(i) == 1) {
	std::cout<<"[2] isGlobal="<<Mu_isGlobalMuon->at(i) << std::endl;
      }
      if (Mu_ptTunePMuonBestTrack->at(i) > 53.0) {
	std::cout<<"[3] ptcocktail="<<Mu_ptTunePMuonBestTrack->at(i) << std::endl;
      }
      if (Mu_absdxyTunePMuonBestTrack->at(i) < 0.2) {
	std::cout<<"[4] absdxy="<<Mu_absdxyTunePMuonBestTrack->at(i) << std::endl;
      }
      if (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i) < 0.10) {
	std::cout<<"[5] trackiso="<<Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i) << std::endl;
      }
      if (Mu_numberOftrackerLayersWithMeasurement->at(i) > 5) {
	std::cout<<"[6] nbTrackerLayer="<<Mu_numberOftrackerLayersWithMeasurement->at(i) << std::endl;
      }
      /*if (Mu_numberOfValidPixelHits->at(i) > 0) {
	std::cout<<"[7] nbPixelHits="<<Mu_numberOfValidPixelHits->at(i) << std::endl;
	}*/
      std::cout<<"[7] nbPixelHits (global tk) ="<<Mu_numberOfValidPixelHits->at(i) << std::endl;
      /*std::cout<<"[7bar] nbPixelHits (inner tk) ="<<Mu_innerTK_numberOfValidPixelHits->at(i) << std::endl;*/
      if (Mu_numberOfValidMuonHits->at(i) > 0) {
	std::cout<<"[8] nbMuonHits="<<Mu_numberOfValidMuonHits->at(i) << std::endl;
      }
      if (Mu_numberOfMatchedStations->at(i) > 1) {
	std::cout<<"[9] nbStation="<<Mu_numberOfMatchedStations->at(i) << std::endl;
      }
      if (Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3) {
	std::cout<<"[10] DeltaPterror="<<Mu_dPToverPTTunePMuonBestTrack->at(i) << std::endl;
      }
      std::cout<<"[11] Charge="<<Mu_chargeTunePMuonBestTrack->at(i) << std::endl;
    }
    std::cout<<"[000] vtxMassMu="        << vtxMass            << std::endl;
    std::cout<<"[000] vtxMassSmearedMu=" << m_vtxMassSmearedMu << std::endl;
    std::cout<<"[000] vtxChi2Mu="        << vtxChi2            << std::endl;
    std::cout<<"[000] CosAngle="         << CosmicRejection    << std::endl;
  }
}

//----------------------------------------------------
//                                                   -
//       Part for HLT & Reco Matching                -
//                                                   -
//----------------------------------------------------
bool ZprimeMuMuPatMiniAodNewMC::isPassHLT()
{
  int nbMatch = 0;
  for (unsigned i=0; i<HLT_nb->size(); i++) {
    if ((HLT_name->at(i) == "HLT_Mu50_v1" ||
         HLT_name->at(i) == "HLT_Mu50_v2" ||
         HLT_name->at(i) == "HLT_Mu50_v3" ||
         HLT_name->at(i) == "HLT_Mu50_v4" ||
         HLT_name->at(i) == "HLT_Mu50_v5" ||
         HLT_name->at(i) == "HLT_Mu50_v6" ||
         HLT_name->at(i) == "HLT_Mu50_v7" ||
         HLT_name->at(i) == "HLT_Mu50_v8" ||
         HLT_name->at(i) == "HLT_Mu50_v9" ||
         HLT_name->at(i) == "HLT_Mu50_v10"||
	 HLT_name->at(i) == "HLT_TkMu50_v1" ||
         HLT_name->at(i) == "HLT_TkMu50_v2" ||
         HLT_name->at(i) == "HLT_TkMu50_v3" ||
         HLT_name->at(i) == "HLT_TkMu50_v4" ||
         HLT_name->at(i) == "HLT_TkMu50_v5" ||
         HLT_name->at(i) == "HLT_TkMu50_v6" ||
         HLT_name->at(i) == "HLT_TkMu50_v7" ||
         HLT_name->at(i) == "HLT_TkMu50_v8" ||
         HLT_name->at(i) == "HLT_TkMu50_v9" ||
         HLT_name->at(i) == "HLT_TkMu50_v10"
	 ) && HLT_isaccept->at(i) == 1 ) {
      //std::cout<<"triggerName = "<<triggerName << std::endl;
      nbMatch++;
    }
  }
  if (nbMatch>0) {
    return true;
  } else return false;
}

bool ZprimeMuMuPatMiniAodNewMC::RecoHLTMuonMatching(float RecoEta,float RecoPhi) {
  int nbMatch = 0;
  float deltaR   = -10000.0;
  for (unsigned i=0; i<HLTObj_nbObj->size(); i++) {
    //std::cout<<"[before]triggerName"<<HLTObj_collection->at(i)  << std::endl;
    if (HLTObj_collection->at(i) == "HLT_Mu50_v1" ||
	HLTObj_collection->at(i) == "HLT_Mu50_v2" ||
	HLTObj_collection->at(i) == "HLT_Mu50_v3" ||
        HLTObj_collection->at(i) == "HLT_Mu50_v4" ||
	HLTObj_collection->at(i) == "HLT_Mu50_v5" ||
        HLTObj_collection->at(i) == "HLT_Mu50_v6" ||
        HLTObj_collection->at(i) == "HLT_Mu50_v7" ||
        HLTObj_collection->at(i) == "HLT_Mu50_v8" ||
        HLTObj_collection->at(i) == "HLT_Mu50_v9" ||
        HLTObj_collection->at(i) == "HLT_Mu50_v10"||
	HLTObj_collection->at(i) == "HLT_TkMu50_v1" ||
	HLTObj_collection->at(i) == "HLT_TkMu50_v2" ||
	HLTObj_collection->at(i) == "HLT_TkMu50_v3" ||
        HLTObj_collection->at(i) == "HLT_TkMu50_v4" ||
	HLTObj_collection->at(i) == "HLT_TkMu50_v5" ||
        HLTObj_collection->at(i) == "HLT_TkMu50_v6" ||
        HLTObj_collection->at(i) == "HLT_TkMu50_v7" ||
        HLTObj_collection->at(i) == "HLT_TkMu50_v8" ||
        HLTObj_collection->at(i) == "HLT_TkMu50_v9" ||
        HLTObj_collection->at(i) == "HLT_TkMu50_v10"
	) {
      //std::cout<<"[after]triggerName"<<HLTObj_collection->at(i)  << std::endl;
      deltaR   = delR(HLTObj_eta->at(i),HLTObj_phi->at(i),RecoEta,RecoPhi);
      //printf ("HLT_Eta = %f  HLT_Phi = %f recoEta = %f recoPhi = %f DelR_trigger = %f\n",HLTObj_eta->at(i),HLTObj_phi->at(i),RecoEta,RecoPhi,deltaR);
      if (fabs(deltaR)>RecoHLTMatchingDeltaRcut) continue;
      nbMatch++;
    }
  }
  if (nbMatch>0) return true;
  else return false;
}



//if (flavor==5 ) b jet
//if (flavor==4 ) c jets
//else light-flavor jet
bool ZprimeMuMuPatMiniAodNewMC::BTaggingDiscriminator(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2) {
  int nbBTag  = 0;
  int nbMatch = 0;
  for (unsigned i=0; i<Nb_bDiscriminators->size(); i++) {
    if (jet_btag_pt->at(i) < 35.0 || fabs(jet_btag_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if (deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if (deltaR2 < 0.5) continue;
    nbBTag++;
    if (nbBTag>1 && nbBTag<3 && jet_btag_pfCSVv2IVF_discriminator->at(i)>0.5426) {
      nbMatch++;
    }
    else continue;
  }
  if (nbMatch>0) return true;
  else return false;
}

bool ZprimeMuMuPatMiniAodNewMC::BTaggingDiscriminator2(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2) {
  int nbBTag  = 0;
  int nbMatch = 0;
  for (unsigned i=0; i<Nb_bDiscriminators->size(); i++) {
    if (jet_btag_pt->at(i) < 35.0 || fabs(jet_btag_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if (deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if (deltaR2 < 0.5) continue;
    nbBTag++;
    if (nbBTag>1 && nbBTag<3 && jet_btag_pfCSVv2IVF_discriminator->at(i)>0.800 &&jet_btag_pfCSVv2IVF_discriminator->at(i)<0.935) {
      nbMatch++;
    }
    else continue;
  }
  if (nbMatch>0) return true;
  else return false;
}

bool ZprimeMuMuPatMiniAodNewMC::BTaggingDiscriminator3(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2) {
  int nbBTag  = 0;
  int nbMatch = 0;
  for (unsigned i=0; i<Nb_bDiscriminators->size(); i++) {
    if (jet_btag_pt->at(i) < 35.0 || fabs(jet_btag_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if (deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if (deltaR2 < 0.5) continue;
    nbBTag++;
    if (nbBTag>1 && nbBTag<3 && jet_btag_pfCSVv2IVF_discriminator->at(i)>0.935) {
      nbMatch++;
    }
    else continue;
  }
  if (nbMatch>0) return true;
  else return false;
}


void ZprimeMuMuPatMiniAodNewMC::DrawBTaggingDiscriminator(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2) {
  int nbBTag = 0;
  for (unsigned i=0; i<Nb_bDiscriminators->size(); i++) {
    if (jet_btag_pt->at(i) < 35.0 || fabs(jet_btag_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if (deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if (deltaR2 < 0.5) continue;
    nbBTag++;
    h1_nbBTagStep1_->Fill(nbBTag,newweight);
    h1_jetBTagStep1_->Fill(jet_btag_pfCSVv2IVF_discriminator->at(i));
    if (nbBTag>1 && nbBTag<3) {
      h1_nbBTagStep2_->Fill(nbBTag,newweight);
      h1_jetBTagStep2_->Fill(jet_btag_pfCSVv2IVF_discriminator->at(i));
    }
    if (nbBTag>1 && nbBTag<3 && jet_btag_pfCSVv2IVF_discriminator->at(i)>0.5426) {
      h1_nbBTagStep3_->Fill(nbBTag,newweight);
      h1_jetBTagStep3_->Fill(jet_btag_pfCSVv2IVF_discriminator->at(i));
    }
  }
}

void ZprimeMuMuPatMiniAodNewMC::Boson(float Px1,float Py1,float Pz1,float En1,
				      float Px2,float Py2,float Pz2,float En2,
				      float ChargeMu1,float MetEt,float MetPx,
				      float MetPy,float MetPz,float MetEn, float &bosonPt) {

  TLorentzVector Mu;
  TLorentzVector Mubar;
  TLorentzVector MissingParticle;
  if (ChargeMu1<0) {
    Mu.SetPxPyPzE(Px1,Py1,Pz1,En1);
    Mubar.SetPxPyPzE(Px2,Py2,Pz2,En2);
  }
  if (ChargeMu1>0) {
    Mu.SetPxPyPzE(Px2,Py2,Pz2,En2);
    Mubar.SetPxPyPzE(Px1,Py1,Pz1,En1);
  }
  TLorentzVector Boson(Mu+Mubar);
  MissingParticle.SetPxPyPzE(MetPx,MetPy,MetPz,MetEn);
  float a = Boson.Angle(MissingParticle.Vect()); // get angle between v1 and v2
  //BosPt  = Boson.Pt();
  //BosPhi = Boson.Phi();
  //  if (Boson.Pt()>60 && MetEt>0 && a>2.8 && (fabs(MetEt-Boson.Pt())/Boson.Pt())<0.4)
  if (Boson.Pt()>60) {
    h1_BosPt_->Fill(Boson.Pt());
    h1_BosPhi_->Fill(Boson.Phi());
    h1_DeltaPhi_->Fill(a,newweight); //Fill(MetPhi-Boson.Phi());
    h1_DeltaPtoverPt_->Fill(fabs(MetEt-Boson.Pt())/Boson.Pt());
    float Mt = sqrt(2.0*Boson.Pt()*MetEt*(1.0-cos(a)));
    h1_Mt_->Fill(Mt,newweight);
    h1_MissingEt_->Fill(MetEt,newweight);
  }
  bosonPt = Boson.Pt();
}

bool ZprimeMuMuPatMiniAodNewMC::DiPFJet(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2) {
  int nbBTag  = 0;
  int nbMatch = 0;
  for (unsigned i=0; i<jet_nb->size(); i++) {
    if (jet_pt->at(i) < 35.0 || fabs(jet_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_eta->at(i),jet_phi->at(i));
    if (deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_eta->at(i),jet_phi->at(i));
    if (deltaR2 < 0.5) continue;
    nbBTag++;
    h1_NbPFjetsAll_->Fill(nbBTag,newweight);
    if (nbBTag>1 && nbBTag<3) {
      h1_NbPFjets2_->Fill(nbBTag,newweight);
      h1_ptPFjetsAll_->Fill(jet_pt->at(i),newweight);
      nbMatch++;
    }
    else continue;
  }
  if (nbMatch>0) return true;
  else return false;
}

bool ZprimeMuMuPatMiniAodNewMC::DiPFJetCut(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2) {
  int nbBTag  = 0;
  int nbMatch = 0;
  for (unsigned i=0; i<jet_nb->size(); i++) {
    if (jet_pt->at(i) < 35.0 || fabs(jet_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_eta->at(i),jet_phi->at(i));
    if (deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_eta->at(i),jet_phi->at(i));
    if (deltaR2 < 0.5) continue;
    nbBTag++;
    if (nbBTag>1 && nbBTag<3) {
      nbMatch++;
    }
    else continue;
  }
  if (nbMatch>0) return true;
  else return false;
}

double ZprimeMuMuPatMiniAodNewMC::MassCorrection(float M, float pT, float Eta1, float Eta2)
{
  float a = 0.;
  float b = 0.;
  float c = 0.;
  float d = 0.;
  float mAdj = M;

  if (pT > 30 && pT < 170) {
    mAdj = M - 170;
    if (fabs(Eta1) <= 1.2 && fabs(Eta2) <= 1.2) { //sigma BB
      a =  1.003;
      b = -0.0002904;
      c = -3.281e-06;
      d =  5.258e-09;
    } else if ((fabs(Eta1) < 1.2 && (fabs(Eta2) > 1.2 && fabs(Eta2) < 2.4)) ||
	       (fabs(Eta2) < 1.2 && (fabs(Eta1) > 1.2 && fabs(Eta1) < 2.4))) {  //BE
      a =  1.012;
      b = -0.001607;
      c =  8.796e-07;
      d =  1.401e-06;
    } else if ((fabs(Eta1) > 1.2 && fabs(Eta1) < 2.4) &&
	       (fabs(Eta2) > 1.2 && fabs(Eta2) < 2.4)) {  //EE
      a =  1.012;
      b = -0.001607;
      c =  8.796e-07;
      d =  1.401e-06;
    } else { // other?
      a =  1.012;
      b = -0.001607;
      c =  8.796e-07;
      d =  1.401e-06;
    }
  } else { // what about pT < 30?
    mAdj = M - 400;
    if (fabs(Eta1) <= 1.2 && fabs(Eta2) <= 1.2) { //sigma BB
      a =  1.036;
      b = -0.0001441;
      c =  5.058e-08;
      d = -7.581e-12;
    } else if ((fabs(Eta1) < 1.2 && (fabs(Eta2) > 1.2 && fabs(Eta2) < 2.4)) ||
	       (fabs(Eta2) < 1.2 && (fabs(Eta1) > 1.2 && fabs(Eta1) < 2.4))) {  //BE
      a =  1.052;
      b = -0.0001471;
      c =  5.903e-08;
      d = -9.037e-12;
    } else if ((fabs(Eta1) > 1.2 && fabs(Eta1) < 2.4) &&
	       (fabs(Eta2) > 1.2 && fabs(Eta2) < 2.4)) {  //EE
      a =  1.052;
      b = -0.0001471;
      c =  5.903e-08;
      d = -9.037e-12;
    } else { // other?
      a =  1.052;
      b = -0.0001471;
      c =  5.903e-08;
      d = -9.037e-12;
    }
  }

  double function = d*pow(mAdj,3) + c*pow(mAdj,2) + b*pow(mAdj,1) + a;
  return function;
}



//============================ Method to di-jet (Barrel-Barrel) ========================
void ZprimeMuMuPatMiniAodNewMC::DrawDiJetMassBB()
{
  float invmass = -10;
  for (unsigned jet1=0; jet1<Mu_nbMuon->size(); jet1++) {
    if (Mu_isGlobalMuon->at(jet1) == 1 &&
	//	Mu_isMuonsCleaned->at(jet1) ==  Mu_isPF->at(jet1) &&
	Mu_ptTunePMuonBestTrack->at(jet1) > FR_Ptcut &&
        Mu_absdxyTunePMuonBestTrack->at(jet1) < 0.2 &&
        (Mu_trackiso->at(jet1)/Mu_ptInnerTrack->at(jet1)) < 0.10  &&
        Mu_numberOftrackerLayersWithMeasurement->at(jet1) > 5 &&
        Mu_numberOfValidPixelHits->at(jet1) > 0 &&
        Mu_numberOfValidMuonHits->at(jet1) > 0 &&
	Mu_passNewMatchedStationsCut->at(jet1) == 1 &&
        Mu_dPToverPTTunePMuonBestTrack->at(jet1) < 0.3 ) continue; //to get rid of real muons
    if ((Mu_isGlobalMuon->at(jet1) == 0 || Mu_isTrackerMuon->at(jet1) == 0) ||
	Mu_ptTunePMuonBestTrack->at(jet1) < FR_Ptcut ||
        Mu_absdxyTunePMuonBestTrack->at(jet1) > 0.2 ||
	Mu_absdzTunePMuonBestTrack->at(jet1) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet1) < 5  ||
        Mu_numberOfValidPixelHits->at(jet1) < 0  ) continue;
    for (unsigned jet2=0; jet2<Mu_nbMuon->size(); jet2++) {
      if (jet2 == jet1) continue;
      if (Mu_chargeTunePMuonBestTrack->at(jet1)*Mu_chargeTunePMuonBestTrack->at(jet2) > 0 ) continue; //OS
      if (Mu_isGlobalMuon->at(jet2) == 1 &&
	  //Mu_isMuonsCleaned->at(jet2) ==  Mu_isPF->at(jet2) &&
	  Mu_ptTunePMuonBestTrack->at(jet2) > FR_Ptcut &&
	  Mu_absdxyTunePMuonBestTrack->at(jet2) < 0.2 &&
	  (Mu_trackiso->at(jet2)/Mu_ptInnerTrack->at(jet2)) < 0.10  &&
	  Mu_numberOftrackerLayersWithMeasurement->at(jet2) > 5 &&
	  Mu_numberOfValidPixelHits->at(jet2) > 0 &&
	  Mu_numberOfValidMuonHits->at(jet2) > 0 &&
	  Mu_passNewMatchedStationsCut->at(jet2) == 1 &&
	  Mu_dPToverPTTunePMuonBestTrack->at(jet2) < 0.3 ) continue; //to get rid of real muons
      if ((Mu_isGlobalMuon->at(jet2) == 0 || Mu_isTrackerMuon->at(jet2) == 0) ||
	  Mu_ptTunePMuonBestTrack->at(jet2) < FR_Ptcut ||
          Mu_absdxyTunePMuonBestTrack->at(jet2) > 0.2 ||
	  Mu_absdzTunePMuonBestTrack->at(jet2) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet2) < 5  ||
	  Mu_numberOfValidPixelHits->at(jet2) < 0  ) continue;
      if (fabs(Mu_etaTunePMuonBestTrack->at(jet1)) > 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet2)) > 1.2) continue;
      float mEl = 0.105658371500;
      float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(jet1),Mu_etaTunePMuonBestTrack->at(jet1),
			     Mu_phiTunePMuonBestTrack->at(jet1),mEl,
			     Mu_ptTunePMuonBestTrack->at(jet2),Mu_etaTunePMuonBestTrack->at(jet2),
			     Mu_phiTunePMuonBestTrack->at(jet2),mEl);
      //pick highest mass dijet
      if (MassDiJet > 60.0) {
	invmass = MassDiJet;
	bool fireHLT1 = isPassHLT();
	if (fireHLT1 == 0) continue;
	bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet1),
							     Mu_phiTunePMuonBestTrack->at(jet1));
	bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet2),
							     Mu_phiTunePMuonBestTrack->at(jet2));
	if (RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
	  float weight1 = FRweight(Mu_etaTunePMuonBestTrack->at(jet1),Mu_ptTunePMuonBestTrack->at(jet1));
	  float weight2 = FRweight(Mu_etaTunePMuonBestTrack->at(jet2),Mu_ptTunePMuonBestTrack->at(jet2));
	  h1_DijetBinWidthBB_->Fill(invmass,(weight1*weight2));
	  h1_DijetBinWidthBBBE_->Fill(invmass,(weight1*weight2));
          h1_Dijet20GeVBB_->Fill(invmass,(weight1*weight2));
	  h1_Dijet20GeVBBBEEE_->Fill(invmass,(weight1*weight2));
	  h1_Dijet1GeVBB_->Fill(invmass,(weight1*weight2));
	  h1_Dijet1GeVBBBEEE_->Fill(invmass,(weight1*weight2));
	}
      }
    }
  }
}

//============================ Method to di-jet (Barrel-Endcaps) ========================
void ZprimeMuMuPatMiniAodNewMC::DrawDiJetMassBE()
{
  float invmass = -10;
  for (unsigned jet1=0; jet1<Mu_nbMuon->size(); jet1++) {
    if (Mu_isGlobalMuon->at(jet1) == 1 &&
	//	Mu_isMuonsCleaned->at(jet1) ==  Mu_isPF->at(jet1) &&
        Mu_ptTunePMuonBestTrack->at(jet1) > FR_Ptcut &&
        Mu_absdxyTunePMuonBestTrack->at(jet1) < 0.2 &&
        (Mu_trackiso->at(jet1)/Mu_ptInnerTrack->at(jet1)) < 0.10  &&
        Mu_numberOftrackerLayersWithMeasurement->at(jet1) > 5 &&
        Mu_numberOfValidPixelHits->at(jet1) > 0 &&
        Mu_numberOfValidMuonHits->at(jet1) > 0 &&
	Mu_passNewMatchedStationsCut->at(jet1) == 1 &&
        Mu_dPToverPTTunePMuonBestTrack->at(jet1) < 0.3 ) continue; //to get rid of real muons
    if ((Mu_isGlobalMuon->at(jet1) == 0 || Mu_isTrackerMuon->at(jet1) == 0) ||
	fabs(Mu_etaTunePMuonBestTrack->at(jet1)) > 2.4 ||
	Mu_ptTunePMuonBestTrack->at(jet1) < FR_Ptcut ||
        Mu_absdxyTunePMuonBestTrack->at(jet1) > 0.2 ||
	Mu_absdzTunePMuonBestTrack->at(jet1) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet1) < 5  ||
        Mu_numberOfValidPixelHits->at(jet1) < 0  ) continue;
    for (unsigned jet2=0; jet2<Mu_nbMuon->size(); jet2++) {
      if (jet2 == jet1) continue;
      if (Mu_chargeTunePMuonBestTrack->at(jet1)*Mu_chargeTunePMuonBestTrack->at(jet2) > 0 ) continue; //OS
      if (Mu_isGlobalMuon->at(jet2) == 1 &&
	  //	  Mu_isMuonsCleaned->at(jet2) ==  Mu_isPF->at(jet2) &&
	  Mu_ptTunePMuonBestTrack->at(jet2) > FR_Ptcut &&
	  Mu_absdxyTunePMuonBestTrack->at(jet2) < 0.2 &&
	  (Mu_trackiso->at(jet2)/Mu_ptInnerTrack->at(jet2)) < 0.10  &&
	  Mu_numberOftrackerLayersWithMeasurement->at(jet2) > 5 &&
	  Mu_numberOfValidPixelHits->at(jet2) > 0 &&
	  Mu_numberOfValidMuonHits->at(jet2) > 0 &&
	  Mu_passNewMatchedStationsCut->at(jet2) == 1 &&
	  Mu_dPToverPTTunePMuonBestTrack->at(jet2) < 0.3 ) continue; //to get rid of real muons
      if ((Mu_isGlobalMuon->at(jet2) == 0 || Mu_isTrackerMuon->at(jet2) == 0) ||
	  fabs(Mu_etaTunePMuonBestTrack->at(jet2)) > 2.4 ||
	  Mu_ptTunePMuonBestTrack->at(jet2) < FR_Ptcut ||
          Mu_absdxyTunePMuonBestTrack->at(jet2) > 0.2 ||
	  Mu_absdzTunePMuonBestTrack->at(jet2) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet2) < 5  ||
	  Mu_numberOfValidPixelHits->at(jet2) < 0  ) continue;
      if (fabs(Mu_etaTunePMuonBestTrack->at(jet1)) > 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet2)) < 1.2) continue;
      float mEl = 0.105658371500;
      float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(jet1),Mu_etaTunePMuonBestTrack->at(jet1),
			     Mu_phiTunePMuonBestTrack->at(jet1),mEl,
			     Mu_ptTunePMuonBestTrack->at(jet2),Mu_etaTunePMuonBestTrack->at(jet2),
			     Mu_phiTunePMuonBestTrack->at(jet2),mEl);
      //pick highest mass dijet
      if (MassDiJet > 60.0) {
	invmass = MassDiJet;
	bool fireHLT1 = isPassHLT();
	if (fireHLT1 == 0) continue;
	bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet1),
							     Mu_phiTunePMuonBestTrack->at(jet1));
	bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet2),
							     Mu_phiTunePMuonBestTrack->at(jet2));
	if (RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
	  float weight1 = FRweight(Mu_etaTunePMuonBestTrack->at(jet1),Mu_ptTunePMuonBestTrack->at(jet1));
	  float weight2 = FRweight(Mu_etaTunePMuonBestTrack->at(jet2),Mu_ptTunePMuonBestTrack->at(jet2));
	  h1_DijetBinWidthBE_->Fill(invmass,(weight1*weight2));
	  h1_DijetBinWidthBBBE_->Fill(invmass,(weight1*weight2));
          h1_DijetEta1_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(jet1)));
	  h1_DijetEta2_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(jet2)));
          h1_Dijet20GeVBEEE_->Fill(invmass,(weight1*weight2));
	  h1_Dijet20GeVBBBEEE_->Fill(invmass,(weight1*weight2));
	  h1_Dijet1GeVBEEE_->Fill(invmass,(weight1*weight2));
	  h1_Dijet1GeVBBBEEE_->Fill(invmass,(weight1*weight2));
	}
      }
    }
  }
}
//============================ Method to di-jet (Endcaps-Endcaps) ========================
void ZprimeMuMuPatMiniAodNewMC::DrawDiJetMassEE()
{
  float invmass = -10;
  for (unsigned jet1=0; jet1<Mu_nbMuon->size(); jet1++) {
    if (Mu_isGlobalMuon->at(jet1) == 1 &&
	//	Mu_isMuonsCleaned->at(jet1) ==  Mu_isPF->at(jet1) &&
        Mu_ptTunePMuonBestTrack->at(jet1) > FR_Ptcut &&
        Mu_absdxyTunePMuonBestTrack->at(jet1) < 0.2 &&
        (Mu_trackiso->at(jet1)/Mu_ptInnerTrack->at(jet1)) < 0.10  &&
        Mu_numberOftrackerLayersWithMeasurement->at(jet1) > 5 &&
        Mu_numberOfValidPixelHits->at(jet1) > 0 &&
        Mu_numberOfValidMuonHits->at(jet1) > 0 &&
	Mu_passNewMatchedStationsCut->at(jet1) == 1 &&
        Mu_dPToverPTTunePMuonBestTrack->at(jet1) < 0.3 ) continue; //to get rid of real muons
    if ((Mu_isGlobalMuon->at(jet1) == 0 || Mu_isTrackerMuon->at(jet1) == 0) ||
	fabs(Mu_etaTunePMuonBestTrack->at(jet1)) > 2.4 ||
	Mu_ptTunePMuonBestTrack->at(jet1) < FR_Ptcut ||
        Mu_absdxyTunePMuonBestTrack->at(jet1) > 0.2 ||
	Mu_absdzTunePMuonBestTrack->at(jet1) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet1) < 5  ||
        Mu_numberOfValidPixelHits->at(jet1) < 0  ) continue;
    for (unsigned jet2=0; jet2<Mu_nbMuon->size(); jet2++) {
      if (jet2 == jet1) continue;
      if (Mu_chargeTunePMuonBestTrack->at(jet1)*Mu_chargeTunePMuonBestTrack->at(jet2) > 0 ) continue; //OS
      if (Mu_isGlobalMuon->at(jet2) == 1 &&
	  // Mu_isMuonsCleaned->at(jet2) ==  Mu_isPF->at(jet2) &&
	  Mu_ptTunePMuonBestTrack->at(jet2) > FR_Ptcut &&
	  Mu_absdxyTunePMuonBestTrack->at(jet2) < 0.2 &&
	  (Mu_trackiso->at(jet2)/Mu_ptInnerTrack->at(jet2)) < 0.10  &&
	  Mu_numberOftrackerLayersWithMeasurement->at(jet2) > 5 &&
	  Mu_numberOfValidPixelHits->at(jet2) > 0 &&
	  Mu_numberOfValidMuonHits->at(jet2) > 0 &&
	  Mu_passNewMatchedStationsCut->at(jet2) == 1 &&
	  Mu_dPToverPTTunePMuonBestTrack->at(jet2) < 0.3 ) continue; //to get rid of real muons
      if ((Mu_isGlobalMuon->at(jet2) == 0 || Mu_isTrackerMuon->at(jet2) == 0) ||
	  fabs(Mu_etaTunePMuonBestTrack->at(jet2)) > 2.4 ||
	  Mu_ptTunePMuonBestTrack->at(jet2) < FR_Ptcut ||
          Mu_absdxyTunePMuonBestTrack->at(jet2) > 0.2 ||
	  Mu_absdzTunePMuonBestTrack->at(jet2) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet2) < 5  ||
	  Mu_numberOfValidPixelHits->at(jet2) < 0  ) continue;
      if (fabs(Mu_etaTunePMuonBestTrack->at(jet1)) < 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet2)) < 1.2) continue;
      float mEl = 0.105658371500;
      float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(jet1),Mu_etaTunePMuonBestTrack->at(jet1),
			     Mu_phiTunePMuonBestTrack->at(jet1),mEl,
			     Mu_ptTunePMuonBestTrack->at(jet2),Mu_etaTunePMuonBestTrack->at(jet2),
			     Mu_phiTunePMuonBestTrack->at(jet2),mEl);
      //pick highest mass dijet
      if (MassDiJet > 60.0) {
	invmass = MassDiJet;
	bool fireHLT1 = isPassHLT();
	if (fireHLT1 == 0) continue;
	bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet1),
							     Mu_phiTunePMuonBestTrack->at(jet1));
	bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet2),
							     Mu_phiTunePMuonBestTrack->at(jet2));
	if (RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
	  float weight1 = FRweight(Mu_etaTunePMuonBestTrack->at(jet1),Mu_ptTunePMuonBestTrack->at(jet1));
	  float weight2 = FRweight(Mu_etaTunePMuonBestTrack->at(jet2),Mu_ptTunePMuonBestTrack->at(jet2));
	  h1_DijetBinWidthEE_->Fill(invmass,(weight1*weight2));
	  h1_DijetBinWidthBBBE_->Fill(invmass,(weight1*weight2));
          h1_Dijet20GeVBEEE_->Fill(invmass,(weight1*weight2));
	  h1_Dijet20GeVBBBEEE_->Fill(invmass,(weight1*weight2));
	  h1_Dijet1GeVEE_->Fill(invmass,(weight1*weight2));
	  h1_Dijet1GeVBEEE_->Fill(invmass,(weight1*weight2));
	  h1_Dijet1GeVBBBEEE_->Fill(invmass,(weight1*weight2));
        }
      }
    }
  }
}

//============================ Method to W-jets (Barrel-Barrel) ========================
void ZprimeMuMuPatMiniAodNewMC::DrawWJetsMassBB()
{
  float invmass = -10;
  for (unsigned muon=0; muon<Mu_nbMuon->size(); muon++) {
    if (Mu_isGlobalMuon->at(muon) == 1 &&
	//        Mu_isMuonsCleaned->at(muon) ==  Mu_isPF->at(muon) &&
	fabs(Mu_etaTunePMuonBestTrack->at(muon)) < 2.4 &&
        Mu_ptTunePMuonBestTrack->at(muon) > FR_Ptcut &&
        Mu_absdxyTunePMuonBestTrack->at(muon) < 0.2 &&
        (Mu_trackiso->at(muon)/Mu_ptInnerTrack->at(muon)) < 0.10  &&
        Mu_numberOftrackerLayersWithMeasurement->at(muon) > 5 &&
        Mu_numberOfValidPixelHits->at(muon) > 0 &&
        Mu_numberOfValidMuonHits->at(muon) > 0 &&
	Mu_passNewMatchedStationsCut->at(muon) == 1 &&
        Mu_dPToverPTTunePMuonBestTrack->at(muon) < 0.3 ) {
      for (unsigned jet=0; jet<Mu_nbMuon->size(); jet++) {
	if (jet == muon) continue;
	if (Mu_chargeTunePMuonBestTrack->at(muon)*Mu_chargeTunePMuonBestTrack->at(jet) > 0 ) continue; //OS
	if (Mu_isGlobalMuon->at(jet) == 1 &&
	    //            Mu_isMuonsCleaned->at(jet) ==  Mu_isPF->at(jet) &&
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) < 2.4 &&
	    Mu_ptTunePMuonBestTrack->at(jet) > FR_Ptcut &&
	    Mu_absdxyTunePMuonBestTrack->at(jet) < 0.2 &&
	    (Mu_trackiso->at(jet)/Mu_ptInnerTrack->at(jet)) < 0.10  &&
	    Mu_numberOftrackerLayersWithMeasurement->at(jet) > 5 &&
	    Mu_numberOfValidPixelHits->at(jet) > 0 &&
	    Mu_numberOfValidMuonHits->at(jet) > 0 &&
	    Mu_passNewMatchedStationsCut->at(jet) == 1 &&
	    Mu_dPToverPTTunePMuonBestTrack->at(jet) < 0.3 ) continue; //to get rid of real muons
	if ((Mu_isGlobalMuon->at(jet) == 0 || Mu_isTrackerMuon->at(jet) == 0) ||
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 2.4 ||
	    Mu_ptTunePMuonBestTrack->at(jet) < FR_Ptcut ||
	    Mu_absdxyTunePMuonBestTrack->at(jet) > 0.2 ||
	    Mu_absdzTunePMuonBestTrack->at(jet) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet) < 5  ||
	    Mu_numberOfValidPixelHits->at(jet) < 0  ) continue;
	float deltaR1 = delR(Mu_etaTunePMuonBestTrack->at(muon),Mu_phiTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(jet),Mu_phiTunePMuonBestTrack->at(jet));
        if (deltaR1 < 0.5) continue;
	if (fabs(Mu_etaTunePMuonBestTrack->at(muon)) > 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 1.2) continue;
	float mEl = 0.105658371500;
	float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(muon),
			       Mu_phiTunePMuonBestTrack->at(muon),mEl,
			       Mu_ptTunePMuonBestTrack->at(jet),Mu_etaTunePMuonBestTrack->at(jet),
			       Mu_phiTunePMuonBestTrack->at(jet),mEl);
	//pick highest mass dijet
	if (MassDiJet > 60.0) {
	  invmass = MassDiJet;
	  bool fireHLT1 = isPassHLT();
	  if (fireHLT1 == 0) continue;
	  bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(muon),
							       Mu_phiTunePMuonBestTrack->at(muon));
	  bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet),
							       Mu_phiTunePMuonBestTrack->at(jet));
	  if (RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
	    float weight = FRweight(Mu_etaTunePMuonBestTrack->at(jet),Mu_ptTunePMuonBestTrack->at(jet));
	    h1_WjetsBinWidthBB_->Fill(invmass,newweight);
	    h1_WjetsBinWidthBBBE_->Fill(invmass,newweight);
            h1_Wjets20GeVBB_->Fill(invmass,newweight);
	    h1_Wjets20GeVBBBEEE_->Fill(invmass,newweight);
	    h1_Wjets1GeVBB_->Fill(invmass,newweight);
            h1_Wjets1GeVBBBEEE_->Fill(invmass,newweight);
	  }
	}
      }
    }
  }
}

//============================ Method to W-jets (Barrel-Endcaps) ========================
void ZprimeMuMuPatMiniAodNewMC::DrawWJetsMassBE1()
{
  float invmass = -10;
  for (unsigned muon=0; muon<Mu_nbMuon->size(); muon++) {
    if (Mu_isGlobalMuon->at(muon) == 1 &&
	//        Mu_isMuonsCleaned->at(muon) ==  Mu_isPF->at(muon) &&
	fabs(Mu_etaTunePMuonBestTrack->at(muon)) < 2.4 &&
        Mu_ptTunePMuonBestTrack->at(muon) > FR_Ptcut &&
        Mu_absdxyTunePMuonBestTrack->at(muon) < 0.2 &&
        (Mu_trackiso->at(muon)/Mu_ptInnerTrack->at(muon)) < 0.10  &&
        Mu_numberOftrackerLayersWithMeasurement->at(muon) > 5 &&
        Mu_numberOfValidPixelHits->at(muon) > 0 &&
        Mu_numberOfValidMuonHits->at(muon) > 0 &&
	Mu_passNewMatchedStationsCut->at(muon) == 1 &&
        Mu_dPToverPTTunePMuonBestTrack->at(muon) < 0.3 ) {
      for (unsigned jet=0; jet<Mu_nbMuon->size(); jet++) {
	if (jet == muon) continue;
	if (Mu_chargeTunePMuonBestTrack->at(muon)*Mu_chargeTunePMuonBestTrack->at(jet) > 0 ) continue; //OS
	if (Mu_isGlobalMuon->at(jet) == 1 &&
	    //            Mu_isMuonsCleaned->at(jet) ==  Mu_isPF->at(jet) &&
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) < 2.4 &&
	    Mu_ptTunePMuonBestTrack->at(jet) > FR_Ptcut &&
	    Mu_absdxyTunePMuonBestTrack->at(jet) < 0.2 &&
	    (Mu_trackiso->at(jet)/Mu_ptInnerTrack->at(jet)) < 0.10  &&
	    Mu_numberOftrackerLayersWithMeasurement->at(jet) > 5 &&
	    Mu_numberOfValidPixelHits->at(jet) > 0 &&
	    Mu_numberOfValidMuonHits->at(jet) > 0 &&
	    Mu_passNewMatchedStationsCut->at(jet) == 1 &&
	    Mu_dPToverPTTunePMuonBestTrack->at(jet) < 0.3 ) continue; //to get rid of real muons
	if ((Mu_isGlobalMuon->at(jet) == 0 || Mu_isTrackerMuon->at(jet) == 0) ||
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 2.4 ||
	    Mu_ptTunePMuonBestTrack->at(jet) < FR_Ptcut ||
	    Mu_absdxyTunePMuonBestTrack->at(jet) > 0.2 ||
	    Mu_absdzTunePMuonBestTrack->at(jet) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet) < 5  ||
	    Mu_numberOfValidPixelHits->at(jet) < 0  ) continue;
	float deltaR1 = delR(Mu_etaTunePMuonBestTrack->at(muon),Mu_phiTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(jet),Mu_phiTunePMuonBestTrack->at(jet));
        if (deltaR1 < 0.5) continue;
	if (fabs(Mu_etaTunePMuonBestTrack->at(muon)) > 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet)) < 1.2) continue;
	float mEl = 0.105658371500;
	float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(muon),
			       Mu_phiTunePMuonBestTrack->at(muon),mEl,
			       Mu_ptTunePMuonBestTrack->at(jet),Mu_etaTunePMuonBestTrack->at(jet),
			       Mu_phiTunePMuonBestTrack->at(jet),mEl);
	//pick highest mass dijet
	if (MassDiJet > 60.0) {
	  invmass = MassDiJet;
	  bool fireHLT1 = isPassHLT();
	  if (fireHLT1 == 0) continue;
	  bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(muon),
							       Mu_phiTunePMuonBestTrack->at(muon));
	  bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet),
							       Mu_phiTunePMuonBestTrack->at(jet));
	  if (RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
	    float weight = FRweight(Mu_etaTunePMuonBestTrack->at(jet),Mu_ptTunePMuonBestTrack->at(jet));
	    h1_WjetsBinWidthBE_->Fill(invmass,newweight);
	    h1_WjetsBinWidthBBBE_->Fill(invmass,newweight);
            h1_Wjets20GeVBEEE_->Fill(invmass,newweight);
	    h1_Wjets20GeVBBBEEE_->Fill(invmass,newweight);
	    h1_Wjets1GeVBEEE_->Fill(invmass,newweight);
            h1_Wjets1GeVBBBEEE_->Fill(invmass,newweight);
	  }
	}
      }
    }
  }
}

void ZprimeMuMuPatMiniAodNewMC::DrawWJetsMassBE2()
{
  float invmass = -10;
  for (unsigned muon=0; muon<Mu_nbMuon->size(); muon++) {
    if (Mu_isGlobalMuon->at(muon) == 1 &&
	//	Mu_isMuonsCleaned->at(muon) ==  Mu_isPF->at(muon) &&
        fabs(Mu_etaTunePMuonBestTrack->at(muon)) < 2.4 &&
        Mu_ptTunePMuonBestTrack->at(muon) > FR_Ptcut &&
        Mu_absdxyTunePMuonBestTrack->at(muon) < 0.2 &&
        (Mu_trackiso->at(muon)/Mu_ptInnerTrack->at(muon)) < 0.10  &&
        Mu_numberOftrackerLayersWithMeasurement->at(muon) > 5 &&
        Mu_numberOfValidPixelHits->at(muon) > 0 &&
        Mu_numberOfValidMuonHits->at(muon) > 0 &&
	Mu_passNewMatchedStationsCut->at(muon) == 1 &&
        Mu_dPToverPTTunePMuonBestTrack->at(muon) < 0.3 ) {
      for (unsigned jet=0; jet<Mu_nbMuon->size(); jet++) {
	if (jet == muon) continue;
	if (Mu_chargeTunePMuonBestTrack->at(muon)*Mu_chargeTunePMuonBestTrack->at(jet) > 0 ) continue; //OS
	if (Mu_isGlobalMuon->at(jet) == 1 &&
	    //            Mu_isMuonsCleaned->at(jet) ==  Mu_isPF->at(jet) &&
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) < 2.4 &&
	    Mu_ptTunePMuonBestTrack->at(jet) > FR_Ptcut &&
	    Mu_absdxyTunePMuonBestTrack->at(jet) < 0.2 &&
	    (Mu_trackiso->at(jet)/Mu_ptInnerTrack->at(jet)) < 0.10  &&
	    Mu_numberOftrackerLayersWithMeasurement->at(jet) > 5 &&
	    Mu_numberOfValidPixelHits->at(jet) > 0 &&
	    Mu_numberOfValidMuonHits->at(jet) > 0 &&
	    Mu_passNewMatchedStationsCut->at(jet) == 1 &&
	    Mu_dPToverPTTunePMuonBestTrack->at(jet) < 0.3 ) continue; //to get rid of real muons
	float deltaR1 = delR(Mu_etaTunePMuonBestTrack->at(muon),Mu_phiTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(jet),Mu_phiTunePMuonBestTrack->at(jet));
        if (deltaR1 < 0.5) continue;
	if ((Mu_isGlobalMuon->at(jet) == 0 || Mu_isTrackerMuon->at(jet) == 0) ||
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 2.4 ||
	    Mu_ptTunePMuonBestTrack->at(jet) < FR_Ptcut ||
	    Mu_absdxyTunePMuonBestTrack->at(jet) > 0.2 ||
	    Mu_absdzTunePMuonBestTrack->at(jet) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet) < 5  ||
	    Mu_numberOfValidPixelHits->at(jet) < 0  ) continue;
	if (fabs(Mu_etaTunePMuonBestTrack->at(muon)) < 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 1.2) continue;
	float mEl = 0.105658371500;
	float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(muon),
			       Mu_phiTunePMuonBestTrack->at(muon),mEl,
			       Mu_ptTunePMuonBestTrack->at(jet),Mu_etaTunePMuonBestTrack->at(jet),
			       Mu_phiTunePMuonBestTrack->at(jet),mEl);
	//pick highest mass dijet
	if (MassDiJet > 60.0) {
	  invmass = MassDiJet;
	  bool fireHLT1 = isPassHLT();
	  if (fireHLT1 == 0) continue;
	  bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(muon),
							       Mu_phiTunePMuonBestTrack->at(muon));
	  bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet),
							       Mu_phiTunePMuonBestTrack->at(jet));
	  if (RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
	    float weight = FRweight(Mu_etaTunePMuonBestTrack->at(jet),Mu_ptTunePMuonBestTrack->at(jet));
	    h1_WjetsBinWidthBE_->Fill(invmass,newweight);
	    h1_WjetsBinWidthBBBE_->Fill(invmass,newweight);
            h1_Wjets20GeVBEEE_->Fill(invmass,newweight);
	    h1_Wjets20GeVBBBEEE_->Fill(invmass,newweight);
	    h1_Wjets1GeVBEEE_->Fill(invmass,newweight);
            h1_Wjets1GeVBBBEEE_->Fill(invmass,newweight);
	  }
	}
      }
    }
  }
}
//============================ Method to W-jets (Endcap-Endcap) ========================
void ZprimeMuMuPatMiniAodNewMC::DrawWJetsMassEE()
{
  float invmass = -10;
  for (unsigned muon=0; muon<Mu_nbMuon->size(); muon++) {
    if (Mu_isGlobalMuon->at(muon) == 1 &&
	//        Mu_isMuonsCleaned->at(muon) ==  Mu_isPF->at(muon) &&
	fabs(Mu_etaTunePMuonBestTrack->at(muon)) < 2.4 &&
        Mu_ptTunePMuonBestTrack->at(muon) > FR_Ptcut &&
        Mu_absdxyTunePMuonBestTrack->at(muon) < 0.2 &&
        (Mu_trackiso->at(muon)/Mu_ptInnerTrack->at(muon)) < 0.10  &&
        Mu_numberOftrackerLayersWithMeasurement->at(muon) > 5 &&
        Mu_numberOfValidPixelHits->at(muon) > 0 &&
        Mu_numberOfValidMuonHits->at(muon) > 0 &&
	Mu_passNewMatchedStationsCut->at(muon) == 1 &&
        Mu_dPToverPTTunePMuonBestTrack->at(muon) < 0.3 ) {
      for (unsigned jet=0; jet<Mu_nbMuon->size(); jet++) {
	if (jet == muon) continue;
	if (Mu_chargeTunePMuonBestTrack->at(muon)*Mu_chargeTunePMuonBestTrack->at(jet) > 0 ) continue; //OS
	if (Mu_isGlobalMuon->at(jet) == 1 &&
	    //            Mu_isMuonsCleaned->at(jet) ==  Mu_isPF->at(jet) &&
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) < 2.4 &&
	    Mu_ptTunePMuonBestTrack->at(jet) > FR_Ptcut &&
	    Mu_absdxyTunePMuonBestTrack->at(jet) < 0.2 &&
	    (Mu_trackiso->at(jet)/Mu_ptInnerTrack->at(jet)) < 0.10  &&
	    Mu_numberOftrackerLayersWithMeasurement->at(jet) > 5 &&
	    Mu_numberOfValidPixelHits->at(jet) > 0 &&
	    Mu_numberOfValidMuonHits->at(jet) > 0 &&
	    Mu_passNewMatchedStationsCut->at(jet) == 1 &&
	    Mu_dPToverPTTunePMuonBestTrack->at(jet) < 0.3 ) continue; //to get rid of real muons
	if ((Mu_isGlobalMuon->at(jet) == 0 || Mu_isTrackerMuon->at(jet) == 0) ||
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 2.4 ||
	    Mu_ptTunePMuonBestTrack->at(jet) < FR_Ptcut ||
	    Mu_absdxyTunePMuonBestTrack->at(jet) > 0.2 ||
	    Mu_absdzTunePMuonBestTrack->at(jet) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet) < 5  ||
	    Mu_numberOfValidPixelHits->at(jet) < 0  ) continue;
	float deltaR1 = delR(Mu_etaTunePMuonBestTrack->at(muon),Mu_phiTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(jet),Mu_phiTunePMuonBestTrack->at(jet));
        if (deltaR1 < 0.5) continue;
	if (fabs(Mu_etaTunePMuonBestTrack->at(muon)) < 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet)) < 1.2) continue;
	float mEl = 0.105658371500;
	float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(muon),
			       Mu_phiTunePMuonBestTrack->at(muon),mEl,
			       Mu_ptTunePMuonBestTrack->at(jet),Mu_etaTunePMuonBestTrack->at(jet),
			       Mu_phiTunePMuonBestTrack->at(jet),mEl);
	//pick highest mass dijet
	if (MassDiJet > 60.0) {
	  invmass = MassDiJet;
	  bool fireHLT1 = isPassHLT();
	  if (fireHLT1 == 0) continue;
	  bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(muon),
							       Mu_phiTunePMuonBestTrack->at(muon));
	  bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet),
							       Mu_phiTunePMuonBestTrack->at(jet));
	  if (RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
	    float weight = FRweight(Mu_etaTunePMuonBestTrack->at(jet),Mu_ptTunePMuonBestTrack->at(jet));
	    h1_WjetsBinWidthEE_->Fill(invmass,newweight);
	    h1_WjetsBinWidthBE_->Fill(invmass,newweight);
            h1_Wjets20GeVBEEE_->Fill(invmass,newweight);
	    h1_Wjets20GeVBBBEEE_->Fill(invmass,newweight);
	    h1_Wjets1GeVEE_->Fill(invmass,newweight);
	    h1_Wjets1GeVBEEE_->Fill(invmass,newweight);
            h1_Wjets1GeVBBBEEE_->Fill(invmass,newweight);
	  }
	}
      }
    }
  }
}

/*
//================= Method to select first high pt muon ==============
float ZprimeMuMuPatMiniAodNewMC::FRweight(float eta, float pt) {
  float FR = 0.0;
  parEB1 = 4.39759e-01;
  parEB2 = -6.62025e+02;
  parEB3 = 1.61874e+03;
  parEB4 = 0.1733;
  parEE1 = 6.62323e-01;
  parEE2 = -4.90668e+02;
  parEE3 = 7.96427e+02;
  parEE4 = 0.3570;
  if (fabs(eta)<1.2 && pt>50.0 && pt<700.0) {
    FR = parEB1 + parEB2 / (parEB3 + pt );
  } else if (fabs(eta)<1.2 && pt>700.0) {
    FR = parEB4;
  } else if (fabs(eta)>1.2 && pt>50.0 && pt<700.0) {
    FR = parEE1 + parEE2 / (parEE3 + pt );
  } else if (fabs(eta)>1.2 && fabs(eta)<2.4 && pt>700.0) {
    FR = parEE4;
  } else {
    std::cout<<"out of FR range" << std::endl;
  }
  return (FR/(1-FR));
}
*/
/*
//================= Method to select first high pt muon ==============
float ZprimeMuMuPatMiniAodNewMC::FRweight(float eta, float pt) {
  float FR = 0.0;
  parEB1 = 1.39489e-01;
  parEB2 = -1.97664e-03;
  parEB3 = 1.08015e-05;
  parEB4 = -1.23318e-08;
  parEB5 = 4.81540e-01;
  parEE1 = 2.12117e-01;
  parEE2 = -3.08913e-03;
  parEE3 = 2.00294e-05;
  parEE4 = -2.51248e-08;
  parEE5 = 2.40380e-01;
  if (fabs(eta)<1.2 && pt>FR_Ptcut && pt<=600.0) {
    FR = parEB1 + parEB2*pow(pt,1) + parEB3*pow(pt,2) + parEB4*pow(pt,3);
  } else if (fabs(eta)<1.2 && pt>600.0) {
    FR = parEB5;
  } else if (fabs(eta)>1.2  && fabs(eta)<2.4 && pt>FR_Ptcut && pt<=600.0) {
    FR = parEE1 + parEE2*pow(pt,1) + parEE3*pow(pt,2) + parEE4*pow(pt,3);
  } else if (fabs(eta)>1.2 && fabs(eta)<2.4 && pt>600.0) {
    FR = parEE5;
  } else {
    std::cout<<"out of FR range" << std::endl;
  }
  return (FR/(1-FR));
}
*/


 //================= Method to select first high pt muon ==============

float ZprimeMuMuPatMiniAodNewMC::FRweight(float eta, float pt) {
  /* without Giovanni filter
     float FR = 0.0;
     parEB1 = 1.11018e-01;
     parEB2 = -1.21961e-03;
     parEB3 = 5.27578e-06;
     parEB4 = 2.40160e+00;
     parEB5 = -5.90161e+03;
     parEB6 = 2.28745e+03;
     parEB7 = 4.04533e-01;
     parEE1 = 1.37269e-01;
     parEE2 =-1.08720e-03;
     parEE3 =  5.62699e-06;
     parEE4 = 2.33319e-01;*/

  /* with Giovanni filter*/
  float FR = 0.0;
  parEB1 = 1.11040e-01;
  parEB2 = -1.21997e-03;
  parEB3 = 5.27748e-06;
  parEB4 = 2.38002e+00;
  parEB5 = -5.72651e+03;
  parEB6 = 2.23475e+03;
  parEB7 = 4.12519e-01;
  parEE1 = 1.37352e-01;
  parEE2 = -1.08994e-03;
  parEE3 = 5.64918e-06;
  parEE4 = 2.38683e-01;

  if (fabs(eta)<1.2 && pt>FR_Ptcut && pt<=200.0) {
      FR = parEB1 + parEB2*pow(pt,1) + parEB3*pow(pt,2);
    } else if (fabs(eta)<1.2 && pt>200 && pt<=800.0) {
      FR = parEB4 + parEB5 / (parEB6 + pt );
    } else if (fabs(eta)<1.2 && pt>800.0) {
      FR = parEB7;
    } else if (fabs(eta)>1.2  && fabs(eta)<2.4 && pt>FR_Ptcut && pt<=250.0) {
      FR = parEE1 + parEE2*pow(pt,1) + parEE3*pow(pt,2);
    } else if (fabs(eta)>1.2 && fabs(eta)<2.4 && pt>250.0) {
      FR = parEE4;
    } else {
    std::cout<<"out of FR range" << std::endl;
  }
  return (FR/(1-FR));
}

/*
float ZprimeMuMuPatMiniAodNewMC::FRweight(float eta, float pt)
{
  float FR = 0.0;
  parEB1  = 1.09268e-01;
  parEB2  = -1.17369e-03;
  parEB3  = 5.01282e-06;
  parEB4  = -1.94865e+00;
  parEB5  = 6.59464e+03;
  parEB6  = 9.59248e+02;
  parEB7  = -1.17298e+06;
  parEB8  = 1.09544e+05;
  parEB9  = 1.45555e+08;
  parEB10 = 2.73227e+07;

  parEE1  =  1.37269e-01;
  parEE2  = -1.08720e-03;
  parEE3  = 5.62699e-06;
  parEE4  = -1.72358e-01;
  parEE5  = 3.52325e+02;
  parEE6  = 2.03179e+01;
  parEE7  = 1.25204e+04;
  parEE8  =  6.45460e+03;
  parEE9  = -2.15325e+07;
  parEE10 = 2.63281e+06;

  if (fabs(eta)<1.2 && pt>FR_Ptcut && pt<=250.0) {
    FR = parEB1 + parEB2*pow(pt,1) + parEB3*pow(pt,2);
  } else if (fabs(eta)<1.2 && pt>250) {
    FR = parEB4+(parEB5/(pt+parEB6))+(parEB7/(pt*pt+parEB8))+(parEB9/(pt*pt*pt+parEB10));
  } else if (fabs(eta)>1.2  && fabs(eta)<2.4 && pt>FR_Ptcut && pt<=250.0) {
    FR = parEE1 + parEE2*pow(pt,1) + parEE3*pow(pt,2);
  } else if (fabs(eta)>1.2 && fabs(eta)<2.4 && pt>250.0) {
    FR = parEE4+(parEE5/(pt+parEE6))+(parEE7/(pt*pt+parEE8))+(parEE9/(pt*pt*pt+parEE10));
  } else {
    std::cout<<"out of FR range" << std::endl;
  }
  return (FR/(1-FR));
}
*/
