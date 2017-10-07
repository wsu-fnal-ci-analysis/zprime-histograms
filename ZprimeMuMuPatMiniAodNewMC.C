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
#include <algorithm>
#include <fstream>
#include <sstream>
using namespace std;
#define PI 3.14159265
bool myfunction (int i,int j) { return (i<j); }
bool picklargemass (float lhs,float rhs) { return (lhs > rhs); }
TString inputfile;
float newweight=1.;

void ZprimeMuMuPatMiniAodNewMC::Loop()
{
time_t start,end;
double dif;
time (&start);
FILE * pFile;
   pFile = fopen ("myfile.txt","w");
   //values needed for AccepXeff study
   NbGen = 0;
   NbReco= 0;
   int binMass   = 10000; //100; //100; //100; //100; //100; //100; //100; //100; //100; //10000;
   float minMass = 0.0;
   float maxMass = 10000.0;
   MassCutMin = 0.0;
   MassCutMax = 2000.0;
   EtaCut = 2.4;
   MassResolution = 0.10;
   deltaRcut = 0.15;
   RecoHLTMatchingDeltaRcut = 0.20;
   minMassCut = 50.0;
   maxMassCut = 4000.0;
   int ptBins = 40;
   float ptMin = 0.0;
   float ptMax = 400.0;
   ptEffCut = 3000.0;
   FR_Ptcut = 53.0; //53.0;
   double muon_mass = 0.1056583;
   weight=1.;
   if( DATA_type=="2016") weight=1.;
   TFile *output = new TFile("CMSSW803-Analyse_ZprimeToMuMu_13TeV.root","recreate");
   //==================================================================================
   //                                                                                 =
   //             Start the histograms for CollinSoper CMF                            =
   //                                                                                 =
   //==================================================================================
   h1_DijetEta1_       = new TH1F("DijetEta1","",100,0.0,3.0);
   h1_DijetEta2_       = new TH1F("DijetEta2","",100,0.0,3.0);
   // Build the histo with constant log bin width
   const int NMBINS2 = 27;
   const double MMIN2 = 60., MMAX2 = 2000.;
   double logMbins2[NMBINS2+1];
   float binNormNr2=0.;
   for (int ibin = 0; ibin <= NMBINS2; ibin++) {
   	logMbins2[ibin] = exp(log(MMIN2) + (log(MMAX2)-log(MMIN2))*ibin/NMBINS2);
   	//cout << logMbins2[ibin] << endl;
   }
   h1_BTagMassMuMuBinWidth_        = new TH1F("BTagMassMuMuBinWidth","BTagMassMuMuBinWidth",NMBINS2, logMbins2);
   h1_MassMuMuDijetBinWidth_       = new TH1F("MassMuMuDijetBinWidth","MassMuMuDijetBinWidth",NMBINS2, logMbins2);
   h1_MassMuMuDijetBinWidthMET_    = new TH1F("MassMuMuDijetBinWidthMET","MassMuMuDijetBinWidthMET",NMBINS2, logMbins2);
   h1_MassMuMuBinWidthMET_         = new TH1F("MassMuMuBinWidthMET","MassMuMuBinWidthMET",NMBINS2, logMbins2);
   h1_MassMuMuDijet1GeVbin_        = new TH1F("MassMuMuDijet1GeVbin","",3000,0.0,3000.0);
   h1_MassMuMuDijet1GeVbinMET_     = new TH1F("MassMuMuDijet1GeVbinMET","",3000,0.0,3000.0);
   h1_MassMuMu1GeVbinMET_          = new TH1F("MassMuMu1GeVbinMET","",3000,0.0,3000.0);
   h1_MassMuMuDijetBinWidthMET100_ = new TH1F("MassMuMuDijetBinWidthMET100","MassMuMuDijetBinWidthMET100",NMBINS2, logMbins2);
   h1_MassMuMuDijet1GeVbinMET100_  = new TH1F("MassMuMuDijet1GeVbinMET100","",3000,0.0,3000.0);

   NbFireHLT = 0;
   int NbBins   = 10;
   float MinBin = -1.0;
   float MaxBin =  1.0;
   h1_ptPFjetsAll_ = new TH1F("ptPFjetsAll","",100,0.0,2000.0);
   h1_NbPFjetsAll_ = new TH1F("NbPFjetsAll","",5,0.0,5.0);
   h1_NbPFjets2_   = new TH1F("NbPFjets2","",5,0.0,5.0);
   h1_BTagMassMuMu_       = new TH1F("BTagMassMuMu","",100,0.0,3000.0);
   h1_BTagMassMuMu1GeVbin_= new TH1F("BTagMassMuMu1GeVbin","",3000,0.0,3000.0);
   h1_nbBTagStep1_   = new TH1F("nbBTagStep1","",5,0.0,5.0);
   h1_jetBTagStep1_  = new TH1F("jetBTagStep1","",50,0.0,1.0);
   h1_nbBTagStep2_   = new TH1F("nbBTagStep2","",5,0.0,5.0);
   h1_jetBTagStep2_  = new TH1F("jetBTagStep2","",50,0.0,1.0);
   h1_nbBTagStep3_   = new TH1F("nbBTagStep3","",5,0.0,5.0);
   h1_jetBTagStep3_  = new TH1F("jetBTagStep3","",50,0.0,1.0);
   h1_MissingEt_     = new TH1F("MissingEt","",100,0.0,2000.0);
   h1_Mt_            = new TH1F("Mt","",100,0.0,2000.0);
   h1_DeltaPtoverPt_ = new TH1F("DeltaPtoverPt","",100,0.0,1.0);
   h1_DeltaPhi_      = new TH1F("DeltaPhi","",50,-4.0,4.0);
   h1_BosPt_         = new TH1F("BosPt","",40,0.0,200.0);
   h1_BosPhi_        = new TH1F("BosPhi","",50,-4.0,4.0);
   h1_jetBTagB_      = new TH1F("jetBTagB","",50,0.0,1.0);
   h1_jetBTagC_      = new TH1F("jetBTagC","",50,0.0,1.0);
   h1_jetBTagUDSG_   = new TH1F("jetBTagUDSG","",50,0.0,1.0);
   h1_PFMetCorr_ = new TH1F("PFMetCorr","",100,0.0,2000.0);
   h1_CaloMet_   = new TH1F("CaloMet","",100,0.0,2000.0);
   h_NbFireHLT   = new TH1F("NbFireHLT", "NbFireHLT", 30 , 0. , 30. );
   h1_ptBeforeTrigger_  = new TH1F("ptBeforeTrigger","",50,0.0,2000.0);
   h1_ptAfterTrigger_   = new TH1F("ptAfterTrigger","",50,0.0,2000.0);
   h1_CosAngleCollinSoperCorrect60Mass120_    = new TH1F("CosAngleCollinSoperCorrect60Mass120","",NbBins,MinBin,MaxBin);
   h1_CosAngleCollinSoperCorrect120Mass300_   = new TH1F("CosAngleCollinSoperCorrect120Mass300","",NbBins,MinBin,MaxBin);
   h1_CosAngleCollinSoperCorrect300Mass700_   = new TH1F("CosAngleCollinSoperCorrect300Mass700","",NbBins,MinBin,MaxBin);
   h1_CosAngleCollinSoperCorrect700Mass3000_  = new TH1F("CosAngleCollinSoperCorrect700Mass3000","",NbBins,MinBin,MaxBin);
   h1_CosAngleCollinSoperCorrect4900Mass5100_ = new TH1F("CosAngleCollinSoperCorrect4900Mass5100","",NbBins,MinBin,MaxBin);
   h1_absCosAngleCollinSoperCorrect4500Mass5500_ = new TH1F("absCosAngleCollinSoperCorrect4500Mass5500","",5,0.0,1.0);
   //==================================================================================
   h1_ptHistoBefor_              = new TH1F("ptHistoBefor","",ptBins,ptMin,ptMax);
   h1_ptHistoPassingVtxChi2Mu_   = new TH1F("ptHistoPassingVtxChi2Mu","",ptBins,ptMin,ptMax);
   h1_ptHistoPassingCosmicRejec_ = new TH1F("ptHistoPassingCosmicRejec","",ptBins,ptMin,ptMax);
   h1_ptHistoPassingHLT_         = new TH1F("ptHistoPassingHLT","",ptBins,ptMin,ptMax);
   h1_etaHistoBefor_              = new TH1F("etaHistoBefor","",30,0.0,3.0);
   h1_etaHistoPassingVtxChi2Mu_   = new TH1F("etaHistoPassingVtxChi2Mu","",30,0.0,3.0);
   h1_etaHistoPassingCosmicRejec_ = new TH1F("etaHistoPassingCosmicRejec","",30,0.0,3.0);
   h1_etaHistoPassingHLT_         = new TH1F("etaHistoPassingHLT","",30,0.0,3.0);
   //==================================================================================
   //                                                                                 =
   //            Start a histograms for Mass resolution                               =
   //                                                                                 =
   //==================================================================================
   h1_MassResultionEBEB1_      = new TH1F("MassResultionEBEB1","",100,-0.5,0.5);
   h1_MassResultionEBEB2_      = new TH1F("MassResultionEBEB2","",100,-0.5,0.5);
   h1_MassResultionEBEB3_      = new TH1F("MassResultionEBEB3","",100,-0.5,0.5);
   h1_MassResultionEBEB4_      = new TH1F("MassResultionEBEB4","",100,-0.5,0.5);
   h1_MassResultionEBEB5_      = new TH1F("MassResultionEBEB5","",100,-0.5,0.5);
   h1_MassResultionEBEB6_      = new TH1F("MassResultionEBEB6","",100,-0.5,0.5);
   h1_MassResultionEBEB7_      = new TH1F("MassResultionEBEB7","",100,-0.5,0.5);
   h1_MassResultionEBEB8_      = new TH1F("MassResultionEBEB8","",100,-0.5,0.5);
   h1_MassResultionEBEB9_      = new TH1F("MassResultionEBEB9","",100,-0.5,0.5);
   h1_MassResultionEBEB10_     = new TH1F("MassResultionEBEB10","",100,-0.5,0.5);
   //==================================================================================
   //                                                                                 =
   //            Start the histograms for the mass of Z                               =
   //                                                                                 =
   //==================================================================================
   h1_ZprimeRecomasslogscale_     = new TH1F("ZprimeRecomasslogscale","",42,log10(60.0),log10(1200.0));
   h2_ZprimeRecomassNewbin_       = new TH2F("ZprimeRecomassNewbin","ZprimeRecomassNewbin",80,50.,1000.,500,0.001,1000.);
   h1_MassGenInAccep_             = new TH1F("MassGenInAccep","",58,0.0,5000.0);
   h1_MassRecoInAccep_            = new TH1F("MassRecoInAccep","",58,0.0,5000.0);
   h1_ZprimeRecomass_             = new TH1F("ZprimeRecomass","",binMass,minMass,maxMass);
   h1_ZprimeRecomass20_           = new TH1F("ZprimeRecomass20","",300,0.0,6000.0);
   h1_ZprimeRecomassBB_           = new TH1F("ZprimeRecomassBB","",binMass,minMass,maxMass);
   h1_ZprimeRecomassEE_           = new TH1F("ZprimeRecomassEE","",binMass,minMass,maxMass);
   h1_ZprimeRecomassBE_           = new TH1F("ZprimeRecomassBE","",binMass,minMass,maxMass);
   h1_ZprimeRecomass50_           = new TH1F("ZprimeRecomass50","",120,0.0,6000.0);
   //h1_ZprimeRecomass60to120_    = new TH1F("ZprimeRecomass60to120","",binMass,minMass,maxMass);
   h1_ZprimeRecomassAbove1000GeV_ = new TH1F("ZprimeRecomassAbove1000GeV","",binMass,minMass,maxMass);
   h1_ZprimeGenmass_              = new TH1F("ZprimeGenmass","",binMass,minMass,maxMass);
   h1_ZprimeGenEta1_              = new TH1F("ZprimeGenEta1","",100,-8.0,8.0);
   h1_ZprimeGenEta2_              = new TH1F("ZprimeGenEta2","",100,-8.0,8.0);
   h1_ZprimeGenPt1_               = new TH1F("ZprimeGenPt1","",100,0.0,2000.0);
   h1_ZprimeGenPt2_               = new TH1F("ZprimeGenPt2","",100,0.0,2000.0);
   h1_ZprimeGenEn1_               = new TH1F("ZprimeGenEn1","",100,0.0,2000.0);
   h1_ZprimeGenEn2_               = new TH1F("ZprimeGenEn2","",100,0.0,2000.0);
   h1_3Dangle_                    = new TH1F("3Dangle","",100,-2.0,2.0);
   h1_DxyDiff_                    = new TH1F("DxyDiff","",100,10.0,10.0);
   h1_MassRecoGenDif_             = new TH1F("MassRecoGenDif","",100,-0.5,0.5);
   h1_PtResolutionTunePMBT_       =  new TH1F("PtResolutionTunePMBT","",100,-0.5,0.5);
   h1_PtResolutiontuneP_          =  new TH1F("PtResolutiontuneP","",100,-0.5,0.5);
   h1_PtResolutionMBT_            =  new TH1F("PtResolutionMBT","",100,-0.5,0.5);
   //==================================================================================
   //                                                                                 =
   //                 Start the histograms for N-1 dist                               =
   //                                                                                 =
   //==================================================================================
   h1_dPToverPT_                            = new TH1F("dPToverPT","",100,0.0,0.5);
   h1_normalizedChi2_                       = new TH1F("normalizedChi2","",100,0.0,20.0);
   h1_numberOftrackerLayersWithMeasurement_ = new TH1F("numberOftrackerLayersWithMeasurement","",20,0.0,20.0);
   h1_numberOfValidPixelHits_               = new TH1F("numberOfValidPixelHits","",10,0.0,10.0);
   h1_numberOfValidMuonHits_                = new TH1F("numberOfValidMuonHits","",60,0.0,60.0);
   h1_numberOfMatchedStations_              = new TH1F("numberOfMatchedStations","",10,0.0,10.0);
   h1_trackiso_                             = new TH1F("trackiso","",50,0.0,0.3);
   h1_absdxy_                               = new TH1F("absdxy","",100,0.0,0.3);
   h1_PtEffpterror_                   = new TH1F("PtEffpterror","",ptBins,ptMin,ptMax);
   h1_PtEffptnumberOftrackerLayers_   = new TH1F("PtEffptnumberOftrackerLayers","",ptBins,ptMin,ptMax);
   h1_PtEffptnumberOfPixelHits_       = new TH1F("PtEffptnumberOfPixelHits","",ptBins,ptMin,ptMax);
   h1_PtEffptnumberOfMuonHits_        = new TH1F("PtEffptnumberOfMuonHits","",ptBins,ptMin,ptMax);
   h1_PtEffptnumberOfMatchedStations_ = new TH1F("PtEffptnumberOfMatchedStations","",ptBins,ptMin,ptMax);
   h1_PtEffptTrackIso_                = new TH1F("PtEffptTrackIso","",ptBins,ptMin,ptMax);
   h1_PtEffptabsdsy_                  = new TH1F("PtEffptabsdsy","",ptBins,ptMin,ptMax);
   h1_PtEffpfSumChargedHadron_        = new TH1F("PtEffpfSumChargedHadron","",ptBins,ptMin,ptMax);
   h1_PtEffpfSumNeutralHadron_        = new TH1F("PtEffpfSumNeutralHadron","",ptBins,ptMin,ptMax);
   h1_PtEffpfPhotonIso_               = new TH1F("PtEffpfPhotonIso","",ptBins,ptMin,ptMax);
   h1_EtaEffpterror_                  = new TH1F("EtaEffpterror","",30,0.0,3.0);
   h1_EtaEffptnumberOftrackerLayers_  = new TH1F("EtaEffptnumberOftrackerLayers","",30,0.0,3.0);
   h1_EtaEffptnumberOfPixelHits_      = new TH1F("EtaEffptnumberOfPixelHits","",30,0.0,3.0);
   h1_EtaEffptnumberOfMuonHits_       = new TH1F("EtaEffptnumberOfMuonHits","",30,0.0,3.0);
   h1_EtaEffptnumberOfMatchedStations_= new TH1F("EtaEffptnumberOfMatchedStations","",30,0.0,3.0);
   h1_EtaEffptTrackIso_               = new TH1F("EtaEffptTrackIso","",30,0.0,3.0);
   h1_EtaEffptabsdsy_                 = new TH1F("EtaEffptabsdsy","",30,0.0,3.0);
   h1_nbPVID_               = new TH1F("nbPVID","",50,0.0,50.0);
   h1_PtID_                 = new TH1F("PtID","",ptBins,ptMin,ptMax);
   h1_EtaID_                = new TH1F("EtaID","",30,0.0,3.0);
   h1_nbPVNewID_            = new TH1F("nbPVNewID","",50,0.0,50.0);
   h1_PtNewID_              = new TH1F("PtNewID","",ptBins,ptMin,ptMax);
   h1_EtaNewID_             = new TH1F("EtaNewID","",30,0.0,3.0);
   h1_nbPVTightID_          = new TH1F("nbPVTightID","",50,0.0,50.0);
   h1_PtTightID_            = new TH1F("PtTightID","",ptBins,ptMin,ptMax);
   h1_EtaTightID_           = new TH1F("EtaTightID","",30,0.0,3.0);
   //h1_3DangleHisto1_        = new TH1F("3DangleHisto1","",50,0.00001,0.1);
   //h1_3DangleHisto2_        = new TH1F("3DangleHisto2","",100,0.00001,10.0);  //100,0.1,10.0);
   //
   h1_3DangleHisto1_        = new TH1F("3DangleHisto1","",1000,0.00001,1.0);
   h1_3DangleHisto2_        = new TH1F("3DangleHisto2","",1000,0.00001,10.0);
   h1_Fail3DangleHistoMass_ = new TH1F("Fail3DangleHistoMass","",100,100.0,12000.0);
   h1_Fail3DangleHistoPhi_  = new TH1F("Fail3DangleHistoPhi","",100,-4.0,4.0);
   //-----------------------------------------------------------------------------------
   h_ptHistoFRDum_ = new TH1F("ptHistoFRDum","",ptBins,ptMin,ptMax);
   h_ptHistoFRNum_ = new TH1F("ptHistoFRNum","",ptBins,ptMin,ptMax);
   h1_PtTuneP_     = new TH1F("PtTuneP","",200,0.0,10000.0);
   //----------------------------------------------------------------------------------------------

   // Build the histo with constant log bin width
   const int NMBINS = 100;
   const double MMIN = 60., MMAX = 6000.;
   double logMbins[NMBINS+1];
   float binNormNr=0.;

   for (int ibin = 0; ibin <= NMBINS; ibin++) {
     logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
     //cout << logMbins[ibin] << endl;
   }
   h1_ZprimeRecomassBinWidthAll_     = new TH1F("ZprimeRecomassBinWidthAll","ZprimeRecomassBinWidthAll",NMBINS, logMbins);
   h1_ZprimeRecomassBinWidth_        = new TH1F("ZprimeRecomassBinWidth","ZprimeRecomassBinWidth",NMBINS, logMbins);
   h1_ZprimeRecomassBinWidthAllBE_   = new TH1F("ZprimeRecomassBinWidthAllBE","ZprimeRecomassBinWidthAllBE",NMBINS, logMbins);
   h1_ZprimeRecomassBinWidthAllEE_   = new TH1F("ZprimeRecomassBinWidthAllEE","ZprimeRecomassBinWidthAllEE",NMBINS, logMbins);
   h1_ZprimeRecomassBinWidthEE_      = new TH1F("ZprimeRecomassBinWidthEE","",NMBINS, logMbins);
   h1_ZprimeRecomassBinWidthBB_      = new TH1F("ZprimeRecomassBinWidthBB","",NMBINS, logMbins);
   h1_ZprimeRecomassBinWidthBEpos_   = new TH1F("ZprimeRecomassBinWidthBEpos","",NMBINS, logMbins);
   h1_ZprimeRecomassBinWidthBEnev_   = new TH1F("ZprimeRecomassBinWidthBEnev","",NMBINS, logMbins);
   h1_ZprimeRecomass60to120BEpos_    = new TH1F("ZprimeRecomass60to120BEpos","",60,60.0,120.0);
   h1_ZprimeRecomass60to120BEnev_    = new TH1F("ZprimeRecomass60to120BEnev","",60,60.0,120.0);
   h1_ZprimeRecomass60to120EE_       = new TH1F("ZprimeRecomass60to120EE","",60,60.0,120.0);
   h1_ZprimeRecomass60to120BB_       = new TH1F("ZprimeRecomass60to120BB","",60,60.0,120.0);
   h1_ZprimeRecomass60to120_         = new TH1F("ZprimeRecomass60to120","",60,60.0,120.0);
   h1_ZprimeRecomassBinWidthAfterBtaging_ = new TH1F("ZprimeRecomassBinWidthAfterBtaging","ZprimeRecomassBinWidthAfterBtaging",NMBINS, logMbins);
   h1_DijetBinWidthBB_         = new TH1F("DijetBinWidthBB","",NMBINS, logMbins);
   h1_DijetBinWidthBE_         = new TH1F("DijetBinWidthBE","",NMBINS, logMbins);
   h1_DijetBinWidthEE_         = new TH1F("DijetBinWidthEE","",NMBINS, logMbins);
   h1_DijetBinWidthBBBE_       = new TH1F("DijetBinWidthBBBE","",NMBINS, logMbins);
   h1_WjetsBinWidthBB_         = new TH1F("WjetsBinWidthBB","",NMBINS, logMbins);
   h1_WjetsBinWidthBE_         = new TH1F("WjetsBinWidthBE","",NMBINS, logMbins);
   h1_WjetsBinWidthEE_         = new TH1F("WjetsBinWidthEE","",NMBINS, logMbins);
   h1_WjetsBinWidthBBBE_       = new TH1F("WjetsBinWidthBBBE","",NMBINS, logMbins);

   h1_Dijet1GeVBB_     = new TH1F("Dijet1GeVBB","",3000,0.0,3000.0);
   h1_Dijet1GeVBEEE_   = new TH1F("Dijet1GeVBEEE","",3000,0.0,3000.0);
   h1_Dijet1GeVEE_     = new TH1F("Dijet1GeVEE","",3000,0.0,3000.0);
   h1_Dijet1GeVBBBEEE_ = new TH1F("Dijet1GeVBBBEEE","",3000,0.0,3000.0);
   h1_Wjets1GeVBB_     = new TH1F("Wjets1GeVBB","",3000,0.0,3000.0);
   h1_Wjets1GeVBEEE_   = new TH1F("Wjets1GeVBEEE","",3000,0.0,3000.0);
   h1_Wjets1GeVEE_     = new TH1F("Wjets1GeVEE","",3000,0.0,3000.0);
   h1_Wjets1GeVBBBEEE_ = new TH1F("Wjets1GeVBBBEEE","",3000,0.0,3000.0);

   h1_Dijet20GeVBB_     = new TH1F("Dijet20GeVBB","",300,0.0,6000.0);
   h1_Dijet20GeVBEEE_   = new TH1F("Dijet20GeVBEEE","",300,0.0,6000.0);
   h1_Dijet20GeVBBBEEE_ = new TH1F("Dijet20GeVBBBEEE","",300,0.0,6000.0);
   h1_Wjets20GeVBB_     = new TH1F("Wjet20GeVBB","",300,0.0,6000.0);
   h1_Wjets20GeVBEEE_   = new TH1F("Wjets20GeVBEEE","",300,0.0,6000.0);
   h1_Wjets20GeVBBBEEE_ = new TH1F("Wjets20GeVBBBEEE","",300,0.0,6000.0);

   // Book txt file for candidate events
   Char_t txtOUT[500];
   sprintf(txtOUT,"CMSSW745-Analyse_ZprimeToMuMu_13TeV_cand.txt");
   output_txt.open(txtOUT);
   output_txt << "CANDIDATES Events:" << endl;
   Char_t outform[20000];
   sprintf (outform,"run: lumi: event: dil_mass: pTmu1: pTmu2: Etamu1: Etamu2:");
   output_txt  << outform << endl;

   TString inputfile=name;
   inputfile=name;
   cout << "Name of the input file is= " << inputfile.Data() << endl;

   //==================================================================================
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     newweight=weight;

     //if( inputfile.Contains("amcatnlo") && !(inputfile.Contains("DYJetsToTauTau"))) {
       //cout << "Reweighting sample of amcatnlo with weight= " << MC_weighting->at(0)<< endl;
       //weight=weight*MC_weighting->at(0);
     //}
     //if (event_runNo==251249) continue;
     /*cout<<"=======> jentry = "<<jentry<<
       "=======> Evt = "<<event_evtNo<<
       "=======> Run = "<<event_runNo<<
       "=======> Lumi = "<<event_lumi<<
       "=======> bunch = "<<event_bunch<<endl; */
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
					  PtRecMuBestTrack1);

     bool secondMuFinal = SelectSecondMuon(ChargeRecMu1,flagmu1,PtRecTunePMuBestTrack1,EtaRecMu1,PhiRecMu1,
					   PtRecTunePMuBestTrack2,EnRecMu2,
					   EtaRecMu2,PhiRecMu2,ChargeRecMu2,pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2,dxyRecMu2,
					   PtRecTunePMu2,PtRecMuBestTrack2);

     PickThehighestMass(vtxMassMu,vtxChi2Mu,event_evtNo);
     double CosmicRejec = ThreeDangle(pxRecMu1,pyRecMu1,pzRecMu1,pRecMu1,
				      pxRecMu2,pyRecMu2,pzRecMu2,pRecMu2);
     //=========================================================
     //        call the method for N-1 plots                   =
     //                                                        =
     //=========================================================
     //cout << "firstMu= " << firstMuFinal << " " << "secondMu= " << secondMuFinal << endl;
     if(firstMuFinal == 0 || secondMuFinal == 0) continue;
     //cout << "Vertex mass mu= " << vtxMassMu << endl;
     //if(vtxMassMu<60 || vtxMassMu>1200) continue;
     if(vtxMassMu<60) continue;
     //=========================================================
     //        start doing matching between reco & HLT         =
     //                                                        =
     //=========================================================
     bool fireHLT2 = isPassHLT();
     if(fireHLT2 == 0) continue;
     bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(EtaRecMu1,PhiRecMu1);
     bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(EtaRecMu2,PhiRecMu2);
     if(RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1)
       {
	 plotAllHighPtMuonsID();
	 //PrintEventInformation(256843,465,665539990,vtxChi2Mu,vtxMassMu,CosmicRejec);
	 if(vtxChi2Mu<20.0 && CosmicRejec>-0.9998)
	   {
	     DrawBTaggingDiscriminator(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
	     if(PFMet_et_cor > 50.0) {
	       h1_PFMetCorr_->Fill(PFMet_et_cor);
	       h1_CaloMet_->Fill(CaloMet_pt);
	       h1_MassMuMuBinWidthMET_->Fill(vtxMassMu);
	       h1_MassMuMu1GeVbinMET_->Fill(vtxMassMu);
	     }
	     bool passDijet = DiPFJet(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
	     if(passDijet==1) {
	       h1_MassMuMuDijetBinWidth_->Fill(vtxMassMu);
	       h1_MassMuMuDijet1GeVbin_->Fill(vtxMassMu);
	     }

	     bool passDijetcuts = DiPFJetCut(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
	     if(passDijetcuts==1 && PFMet_et_cor > 50.0) {
	       h1_MassMuMuDijetBinWidthMET_->Fill(vtxMassMu);
	       h1_MassMuMuDijet1GeVbinMET_->Fill(vtxMassMu);
	     }
             if(passDijetcuts==1 && PFMet_et_cor > 100.0) {
               h1_MassMuMuDijetBinWidthMET100_->Fill(vtxMassMu);
               h1_MassMuMuDijet1GeVbinMET100_->Fill(vtxMassMu);
             }
             bool passBTaggingDiscriminator2 = BTaggingDiscriminator2(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
	     if(passBTaggingDiscriminator2==1) {
	       h1_BTagMassMuMu_->Fill(vtxMassMu);
	     }
	     bool passBTaggingDiscriminator3 = BTaggingDiscriminator3(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
	     if(passBTaggingDiscriminator3==1) {
	       h1_BTagMassMuMu_->Fill(vtxMassMu);
             }

             // WW sample
             //if( datasetName.Contains("WWTo2L2Nu")) {
             if ( inputfile.Contains("WWTo2L2Nu")){
               //cout << "Reweighting sample of WWTo2L2Nu with weight=0" << endl;
               //if (MassGen>600.) newweight=0;
               if (vtxMassMu>600.) {
		 cout << "Reweighting sample of WWTo2L2Nu with weight=0" << endl;
		 newweight=0;
	       }
             }

             // TTTo2L2Nu sample
             //if( datasetName.Contains("TTTo2L2Nu")) {
             if ( inputfile.Contains("TTTo2L2Nu")){
               //cout << "Reweighting sample of TTTo2L2Nu with weight=0" << endl;
               //if (MassGen>500.) newweight=0;
               if (vtxMassMu>500.) {
		 cout << "Reweighting sample of TTTo2L2Nu with weight=0" << endl;
		 newweight=0;
	       }
             }

             // CITo2Mu_Lam10TeV_LLConM300 sample
             if ( inputfile.Contains("CITo2Mu_Lam10TeV_LLConM300")){
               //cout << "Reweighting sample of CITo2Mu_Lam10TeV_LLConM300 with weight=0" << endl;
               if (vtxMassMu<300.) {
		 cout << "Reweighting sample of CITo2Mu_Lam10TeV_LLConM300 with weight=0" << endl;
		 newweight=0;
	       }
             }

             // CITo2Mu_Lam10TeV_LLConM800 sample

             //if ( inputfile.Contains("CITo2Mu_Lam10TeV_LLConM800")){
             //  cout << "Reweighting sample of CITo2Mu_Lam10TeV_LLConM800 with weight=0" << endl;
             //  if (vtxMassMu<800.) weight=0;
             //}

	     PlotRecoInfo(CosmicRejec,vtxMassMu,MassGen,PtRecTunePMuBestTrack1,PtRecTunePMu1,PtRecMuBestTrack1,mPtGen1,EtaRecMu1,
			  PtRecTunePMuBestTrack2,PtRecTunePMu2,PtRecMuBestTrack2,mPtGen2,EtaRecMu2);
	     CosThetaCollinSoper(PtRecTunePMuBestTrack1,EtaRecMu1,PhiRecMu1,EnRecMu1,
				 PtRecTunePMuBestTrack2,EtaRecMu2,PhiRecMu2,EnRecMu2,
				 ChargeRecMu1,vtxMassMu);
	     Boson(pxRecMu1,pyRecMu1,pzRecMu1,EnRecMu1,pxRecMu2,pyRecMu2,pzRecMu2,EnRecMu2,
		   ChargeRecMu1,PFMet_et_cor,PFMet_px_cor,PFMet_py_cor,PFMet_pz_cor,PFMet_en_cor);
	     bool passBTaggingDiscriminator = BTaggingDiscriminator(EtaRecMu1,PhiRecMu1,EtaRecMu2,PhiRecMu2);
	     if(passBTaggingDiscriminator==1)
	       {
		 h1_BTagMassMuMuBinWidth_->Fill(vtxMassMu);
		 h1_BTagMassMuMu1GeVbin_->Fill(vtxMassMu);
	       }
	     if(passBTaggingDiscriminator==0)
	       {
		 h1_ZprimeRecomassBinWidthAfterBtaging_->Fill(vtxMassMu);
	       }
	   }
       }
   }
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
 float ZprimeMuMuPatMiniAodNewMC::delR(float eta1,float phi1,float eta2,float phi2){
   float mpi=3.14;
   float dp=std::abs(phi1-phi2);
   if (dp>mpi) dp-=float(2*mpi);
   return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
 }

 bool ZprimeMuMuPatMiniAodNewMC::GenRecoMatchMu(float RecoEta1,float RecoPhi1){
   int NbHEEPele = 0;
   unsigned iflag = -10;
   for(unsigned i=0; i<iGen->size(); i++){
     float deltaR1   = delR(RecoEta1,RecoPhi1,etaGen->at(i),phiGen->at(i));
     if( fabs(idGen->at(i)) != 13 ) continue;
     if( statusGen->at(i) != 1 )  continue;
     if(fabs(deltaR1)>deltaRcut) continue;
     iflag  = i;
     NbHEEPele ++;
   }
   if(NbHEEPele > 0) {
     return true;
   }
   else return false;
 }
//============================ Method to select first high pt muon ========================
bool ZprimeMuMuPatMiniAodNewMC::SelectFirstMuon(float &pTmuon1,float &Enmuon1,float &Etamuon1,
						  float &Phimuon1,int &ChargeMu1,unsigned &FlagMu1,
						  float &pxmuon1,float &pymuon1,float &pzmuon1,
						  float &pmuon1,float &dxymuon1,float &pTmuon1tuneP,
						  float &pTmuonBestTrack1){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  float highestpt=-999.;

  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    //    if(Mu_isMuonsCleaned->at(i) != Mu_isPF->at(i)) continue;
    if( Mu_isTrackerMuon->at(i) == 1 &&
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
	//bool GenRecoMatch1 = GenRecoMatchMu(Mu_etaTunePMuonBestTrack->at(i),Mu_phiTunePMuonBestTrack->at(i));
	//if(GenRecoMatch1 == 0) continue;
	highestpt=Mu_ptTunePMuonBestTrack->at(i);
	iflag  = i;
	NbHEEPele ++;
      }
    }
    else continue;
  }

  if( NbHEEPele > 0 ){
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
    //cout << "First Muon ChargeMu1= " << ChargeMu1 << " Pt1= " << pTmuonBestTrack1 << endl;
    return true;
  }
  else return false;
}
//============================ Method to select second high pt muon ========================
bool ZprimeMuMuPatMiniAodNewMC::SelectSecondMuon(int ChargeMu1,unsigned FlagMu1,float pTmuon1,float Etamuon1,float Phimuon1,
						   float &pTmuon2,float &Enmuon2,
						   float &Etamuon2,float &Phimuon2,int &ChargeMu2,float &pxmuon2,
						   float &pymuon2,float &pzmuon2,float &pmuon2,float &dxymuon2,
						   float &pTmuon2tuneP,float &pTmuonBestTrack2)
{
  int NbHEEPele = 0;
  unsigned iflag = -10;
  float highestpt=-999.;
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(i == FlagMu1) continue;
    if(Mu_ptTunePMuonBestTrack->at(i) == pTmuon1) continue;
    if(Mu_etaTunePMuonBestTrack->at(i) == Etamuon1) continue;
    if(Mu_phiTunePMuonBestTrack->at(i) == Phimuon1) continue;
    //cout << "Charge=" << ChargeMu1 << " " << Mu_chargeTunePMuonBestTrack->at(i) << " pT= " <<  Mu_ptTunePMuonBestTrack->at(i) <<  endl;
    if(ChargeMu1*Mu_chargeTunePMuonBestTrack->at(i)>0) continue;
    //if(ChargeMu1*Mu_chargeTunePMuonBestTrack->at(i)<0) continue;
    //    if(Mu_isMuonsCleaned->at(i) != Mu_isPF->at(i)) continue;
    if( Mu_isTrackerMuon->at(i) == 1 &&
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
      if(Mu_ptTunePMuonBestTrack->at(i)>highestpt) {
	//bool GenRecoMatch2 = GenRecoMatchMu(Mu_etaTunePMuonBestTrack->at(i),Mu_phiTunePMuonBestTrack->at(i));
	//if(GenRecoMatch2 == 0) continue;
	highestpt=Mu_ptTunePMuonBestTrack->at(i);
	//cout << "Highest PT second lepton has pt= " << highestpt << endl;
	iflag  = i;
	NbHEEPele ++;
      }
    }
    else continue;
  }
  if( NbHEEPele > 0 ){
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
    //cout << "Second ChargeMu2= " << ChargeMu2 << " Pt2= " << pTmuonBestTrack2 << endl;
    return true;
  }
  else return false;
}
void ZprimeMuMuPatMiniAodNewMC::PlotRecoInfo(float CosmicMuonRejec, float vertexMassMu,float MassGenerated,
				 float PtTunePMuBestTrack,float PtTunePMu,float PtMuBestTrack,
				 float PtGenerated, float etaMu1,
				 float PtTunePMuBestTrack2,float PtTunePMu2,float PtMuBestTrack2,
				 float PtGenerated2,float etaMu2){
  //----------------------------------------------------------
  if (vertexMassMu>900.0) {
  output_txt << event_runNo
             << "   " << event_lumi
             << "       " << event_evtNo
             << "        " << vertexMassMu
             << "        " << PtTunePMuBestTrack
             << "        " << PtTunePMuBestTrack2
             << "        " << etaMu1
             << "        " << etaMu2
             << endl;

   }



  float weight10 = MassCorrection(vertexMassMu);
  //float weight2 = std::min(1.01696-7.73522E-5*vertexMassMu+6.69239E-9*vertexMassMu*vertexMassMu,1);
  //----------------------------------------------------------
  h1_ZprimeRecomassBinWidth_->Fill(vertexMassMu,newweight);
  h1_ZprimeRecomasslogscale_->Fill(log10(vertexMassMu));
  h1_ZprimeRecomass_->Fill(vertexMassMu,newweight);
  h1_MassRecoInAccep_->Fill(MassGenerated,newweight);
  //if (!(inputfile.Contains("WW") && vertexMassMu>2000.) ) {
    //BB
    if(fabs(etaMu1)<1.2 && fabs(etaMu2)<1.2)
      {
	h1_ZprimeRecomassBinWidthBB_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120BB_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
        h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
      }
    //BE
    if( fabs(etaMu1)<1.2 && (fabs(etaMu2) > 1.2 && fabs(etaMu2) < 2.4) )
      {
	h1_ZprimeRecomassBinWidthAllBE_->Fill(vertexMassMu,newweight);
      }
    if( fabs(etaMu2)<1.2 && (fabs(etaMu1) > 1.2 && fabs(etaMu1) < 2.4) )
      {
	h1_ZprimeRecomassBinWidthAllBE_->Fill(vertexMassMu,newweight);
      }
    //EE
    if( (fabs(etaMu1) > 1.2 && fabs(etaMu1) < 2.4) && (fabs(etaMu2) > 1.2 && fabs(etaMu2) < 2.4) )
      {
	h1_ZprimeRecomassBinWidthAllEE_->Fill(vertexMassMu,newweight);
      }
    //============================================================================
    //BE+
    if( fabs(etaMu1)<1.2 && (etaMu2 > 1.2 && etaMu2 < 2.4) )
      {
	h1_ZprimeRecomassBinWidthBEpos_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120BEpos_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
      }

    if( fabs(etaMu2)<1.2 && (etaMu1 > 1.2 && etaMu1 < 2.4) )
      {
	h1_ZprimeRecomassBinWidthBEpos_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomassBinWidthBEpos_->Fill(vertexMassMu,newweight);
	//	h1_ZprimeRecomass_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
      }

    //BE-
    if( fabs(etaMu1)<1.2 && (etaMu2 > -2.4 && etaMu2 < -1.2) )
      {
	h1_ZprimeRecomassBinWidthBEnev_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120BEnev_->Fill(vertexMassMu,newweight);
	//h1_ZprimeRecomass_->Fill(vertexMassMu,newweight);
	//h1_ZprimeRecomassBinWidthAll_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
      }

    if( fabs(etaMu2)<1.2 && (etaMu1 > -2.4 && etaMu1 < -1.2) )
      {
	h1_ZprimeRecomassBinWidthBEnev_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120BEnev_->Fill(vertexMassMu,newweight);
	//h1_ZprimeRecomass_->Fill(vertexMassMu,newweight);
	//h1_ZprimeRecomassBinWidth_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
      }

    //E+E-
    if( (etaMu1 > 1.2 && etaMu1 < 2.4) && (etaMu2 > -2.4 && etaMu2 < -1.2) )
      {
	h1_ZprimeRecomassBinWidthEE_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120EE_->Fill(vertexMassMu,newweight);
	//h1_ZprimeRecomass_->Fill(vertexMassMu,newweight);
	//h1_ZprimeRecomassBinWidth_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
      }
    if( (etaMu2 > 1.2 && etaMu2 < 2.4) && (etaMu1 > -2.4 && etaMu1 < -1.2) )
      {
	h1_ZprimeRecomassBinWidthEE_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120EE_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass_->Fill(vertexMassMu,newweight);
	//h1_ZprimeRecomassBinWidth_->Fill(vertexMassMu,newweight);
	h1_ZprimeRecomass60to120_->Fill(vertexMassMu,newweight);
      }
    //}
  h1_ZprimeRecomass50_->Fill(vertexMassMu);
  h1_ZprimeRecomass20_->Fill(vertexMassMu);
  if(fabs(etaMu1)<1.2 && fabs(etaMu2)<1.2){h1_ZprimeRecomassBB_->Fill(vertexMassMu);}
  if(fabs(etaMu1)>1.2 && fabs(etaMu2)>1.2){h1_ZprimeRecomassEE_->Fill(vertexMassMu);}
  if(fabs(etaMu1)<1.2 && fabs(etaMu2)>1.2){h1_ZprimeRecomassBE_->Fill(vertexMassMu);}
  if(fabs(etaMu1)>1.2 && fabs(etaMu2)<1.2){h1_ZprimeRecomassBE_->Fill(vertexMassMu);}
  //if(vertexMassMu>60 && vertexMassMu<120) h1_ZprimeRecomass60to120_->Fill(vertexMassMu);
  h1_3Dangle_->Fill(CosmicMuonRejec,newweight);
  //part for Pt resolution
  h1_PtResolutionTunePMBT_->Fill((PtTunePMuBestTrack-PtGenerated)/PtGenerated,newweight);
  h1_PtResolutiontuneP_->Fill((PtTunePMu-PtGenerated)/PtGenerated,newweight);
  h1_PtResolutionMBT_->Fill((PtMuBestTrack-PtGenerated)/PtGenerated,newweight);
  //part for mass resolution
  if( vertexMassMu > 0.0 && vertexMassMu < 250.0 ){h1_MassResultionEBEB1_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);}
  if( vertexMassMu > 250 && vertexMassMu < 750.0 ){h1_MassResultionEBEB2_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);}
  if( vertexMassMu > 750 && vertexMassMu < 1250.0 ){h1_MassResultionEBEB3_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);}
  if( vertexMassMu > 1250 && vertexMassMu < 1750.0 ){h1_MassResultionEBEB4_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);}
  if( vertexMassMu > 1750 && vertexMassMu < 2250.0 ){h1_MassResultionEBEB5_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);}
  if( vertexMassMu > 2000 && vertexMassMu < 4000.0 ){h1_MassResultionEBEB6_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);}
  if( vertexMassMu > 4000 && vertexMassMu < 6000.0 ){h1_MassResultionEBEB7_->Fill((vertexMassMu-MassGenerated)/MassGenerated,newweight);}
}
//===================== Methode to calculate the mass ========================
float ZprimeMuMuPatMiniAodNewMC::Mass(float Pt1,float Eta1,float Phi1,float En1,
              		                float Pt2,float Eta2,float Phi2,float En2){
  float MuMuMass = 0.0;
  TLorentzVector Mu1;
  TLorentzVector Mu2;
  Mu1.SetPtEtaPhiM(Pt1,Eta1,Phi1,En1);
  Mu2.SetPtEtaPhiM(Pt2,Eta2,Phi2,En2);
  MuMuMass = (Mu1 + Mu2).M();
  return MuMuMass;
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
  for(unsigned i=0; i<Mu_vtxMass->size(); i++)
    {
      Nb++;
      countlept=2*i;
      //cout << "vtx mass" << Mu_vtxMass->at(i) << " Chi2= " << Mu_vtxNormChi2->at(i)<< endl;
      //cout << "vtx Mass lepton= " << Mu_vtxMassLept->at(countlept) << " " <<  Mu_vtxMassLept->at(countlept+1)<< endl;
      float chargepair=0;
      for(unsigned j=0; j<Mu_nbMuon->size(); j++){
	if (Mu_ptTunePMuonBestTrack->at(j)==Mu_vtxMassLept->at(countlept)) chargepair=Mu_chargeTunePMuonBestTrack->at(j);
	if (Mu_ptTunePMuonBestTrack->at(j)==Mu_vtxMassLept->at(countlept+1)) chargepair=chargepair*Mu_chargeTunePMuonBestTrack->at(j);
      }

      //cout << "Chargepair for vtxmass= " <<  Mu_vtxMass->at(i) << "   is " << chargepair << endl;
      if (chargepair!=-1) continue;
      //cout << "Chargepair for vtxmass= " <<  Mu_vtxMass->at(i) << "   is " << chargepair << endl;


      int leptmatchBest=0;
      if (Mu_vtxMassLept->at(countlept)==PtRecTunePMu1 && Mu_vtxMassLept->at(countlept+1)==PtRecTunePMu2) leptmatchBest=2;
      if (Mu_vtxMassLept->at(countlept)==PtRecTunePMu2 && Mu_vtxMassLept->at(countlept+1)==PtRecTunePMu1) leptmatchBest=2;

      if (leptmatchBest!=2) continue;
      //cout << "Chargepair surviving matching for vtxmass= " <<  Mu_vtxMass->at(i) << "   with " << leptmatchBest << " leptons" << endl;

      if(Mu_vtxNormChi2->at(i)> 20) continue;
      if(Mu_vtxMass->at(i)>Massinv){
	Massinv = Mu_vtxMass->at(i);
	iflag  = i;
	NbMu++;
      }
    }
  if( NbMu > 0 ){
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
bool ZprimeMuMuPatMiniAodNewMC::SelectFirstGenMu(float &ETMu1,float &PhiSCMu1,
					    float &EtaSCMu1,float &EnMu1,
					    int &IDele1,int &Statele1,
					    unsigned &GenFlag1){
  int NbHEEPele = 0;
  int iflag = -10;
  ETMu1 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    if( fabs(idGen->at(i)) != 13 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if( ptGen->at(i) > ETMu1) {
      ETMu1 = ptGen->at(i);
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if(NbHEEPele>0) {
    GenFlag1       = iflag;
    ETMu1          = ptGen->at(iflag);
    PhiSCMu1       = phiGen->at(iflag);
    EtaSCMu1       = etaGen->at(iflag);
    EnMu1          = EnergyGen->at(iflag);
    IDele1         = idGen->at(iflag);
    Statele1       = statusGen->at(iflag);
    return true;
  }
  else return false;
}
//============================ Method to select second Gen Mu ========================
bool ZprimeMuMuPatMiniAodNewMC::SelectSecondGenMu(unsigned GenFlag1,float ETMu1,float &ETMu2,float &PhiSCMu2,
					     float &EtaSCMu2,float &EnMu2,int &IDele2,int &Statele2){
  int NbHEEPele = 0;
  int iflag = -10;
  ETMu2 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    if( fabs(idGen->at(i)) != 13 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if(i == GenFlag1) continue;
    if( fabs(ptGen->at(i) - ETMu1) <0.00001 ) continue;
    if( ptGen->at(i) > ETMu2) {
      ETMu2 = ptGen->at(i);
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if(NbHEEPele>0) {
    ETMu2      = ptGen->at(iflag);
    PhiSCMu2   = phiGen->at(iflag);
    EtaSCMu2   = etaGen->at(iflag);
    EnMu2      = EnergyGen->at(iflag);
    IDele2     = idGen->at(iflag);
    Statele2   = statusGen->at(iflag);
    return true;
  }
  else return false;
}


//============================ Method to plot Gen Mu ========================
void ZprimeMuMuPatMiniAodNewMC::PlotGenInfo(float ZprimeGenMass,float EtaGenMu1,float EtaGenMu2,float PtGenMu1,
				       float PtGenMu2,float EnGenMu1,float EnGenMu2){
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
void ZprimeMuMuPatMiniAodNewMC::MuonPassingID(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_PtID_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaID_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
      h1_nbPVID_->Fill(Mu_nbofpv->at(i),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotPterror(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1){
      h1_dPToverPT_->Fill(Mu_dPToverPTTunePMuonBestTrack->at(i),newweight );
      h1_PtEffpterror_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffpterror_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotNbTrackLayers(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       //Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_numberOftrackerLayersWithMeasurement_->Fill( Mu_numberOftrackerLayersWithMeasurement->at(i),newweight );
      h1_PtEffptnumberOftrackerLayers_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptnumberOftrackerLayers_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}


void ZprimeMuMuPatMiniAodNewMC::PlotNBValidPixelHits(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       //Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_numberOfValidPixelHits_->Fill( Mu_numberOfValidPixelHits->at(i),newweight );
      h1_PtEffptnumberOfPixelHits_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptnumberOfPixelHits_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotNbValidMuonHits(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       //Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_numberOfValidMuonHits_->Fill( Mu_numberOfValidMuonHits->at(i),newweight );
      h1_PtEffptnumberOfMuonHits_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptnumberOfMuonHits_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}


void ZprimeMuMuPatMiniAodNewMC::PlotNbMatchedStations(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       //Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_numberOfMatchedStations_->Fill( Mu_numberOfMatchedStations->at(i),newweight );
      h1_PtEffptnumberOfMatchedStations_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptnumberOfMatchedStations_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}


void ZprimeMuMuPatMiniAodNewMC::PlotTrackiso(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       //(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_trackiso_->Fill( Mu_trackiso->at(i)/Mu_ptTunePMuonBestTrack->at(i),newweight );
      h1_PtEffptTrackIso_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptTrackIso_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}


void ZprimeMuMuPatMiniAodNewMC::PlotAbsDxy(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       //Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_absdxy_->Fill( Mu_absdxyTunePMuonBestTrack->at(i) ,newweight);
      h1_PtEffptabsdsy_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaEffptabsdsy_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
    }
    else continue;
  }
}

void ZprimeMuMuPatMiniAodNewMC::PlotPtTuneP(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.2 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_PtTuneP_->Fill( Mu_ptTunePMuonBestTrack->at(i) ,newweight);
    }
    else continue;
  }
}



void ZprimeMuMuPatMiniAodNewMC::plotAllHighPtMuonsID(){
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

void ZprimeMuMuPatMiniAodNewMC::MuonPassingNewID(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.02 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_PtNewID_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaNewID_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
      h1_nbPVNewID_->Fill(Mu_nbofpv->at(i),newweight);
    }
    else continue;
  }
}



void ZprimeMuMuPatMiniAodNewMC::MuonPassingTightID(){
  for(unsigned i=0; i<Mu_nbMuon->size(); i++){
    if(fabs(Mu_etaTunePMuonBestTrack->at(i)) < 2.4 &&
       Mu_isGlobalMuon->at(i) == 1 &&
       (Mu_ptTunePMuonBestTrack->at(i) > 53.0 && Mu_ptTunePMuonBestTrack->at(i) < ptEffCut) &&
       Mu_absdxyTunePMuonBestTrack->at(i) < 0.01 &&
       (Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)) < 0.10  &&
       Mu_numberOftrackerLayersWithMeasurement->at(i) > 5 &&
       Mu_numberOfValidPixelHits->at(i) > 0 &&
       Mu_numberOfValidMuonHits->at(i) > 0 &&
       Mu_numberOfMatchedStations->at(i) > 1 &&
       Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3 ){
      h1_PtTightID_->Fill(Mu_ptTunePMuonBestTrack->at(i),newweight);
      h1_EtaTightID_->Fill(fabs(Mu_etaTunePMuonBestTrack->at(i)),newweight);
      h1_nbPVTightID_->Fill(Mu_nbofpv->at(i),newweight);
    }
    else continue;
  }
}


void ZprimeMuMuPatMiniAodNewMC::CosThetaCollinSoper(float Et1,float Eta1,float Phi1,float En1,
						  float Et2,float Eta2,float Phi2,float En2,
						  float ChargeEle1,float RecoMass){

  TLorentzVector Ele;
  TLorentzVector Elebar;
  if(ChargeEle1<0) {
    Ele.SetPtEtaPhiE(Et1,Eta1,Phi1,En1);
    Elebar.SetPtEtaPhiE(Et2,Eta2,Phi2,En2);
  }
  if(ChargeEle1>0) {
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

  if( RecoMass > 60.0 && RecoMass < 120.0 )
    {
      h1_CosAngleCollinSoperCorrect60Mass120_->Fill(costheta,newweight);
    }


  if( RecoMass > 120.0 && RecoMass < 300.0 )
    {
      h1_CosAngleCollinSoperCorrect120Mass300_->Fill(costheta,newweight);
    }

  if( RecoMass > 300.0 && RecoMass < 700.0 )
    {
      h1_CosAngleCollinSoperCorrect300Mass700_->Fill(costheta,newweight);
    }

  if( RecoMass > 700.0 && RecoMass < 3000.0 )
    {
      h1_CosAngleCollinSoperCorrect700Mass3000_->Fill(costheta,newweight);
    }

  if( RecoMass > 4500.0 && RecoMass < 6000.0 )
    {
      h1_CosAngleCollinSoperCorrect4900Mass5100_->Fill(costheta,newweight);
      h1_absCosAngleCollinSoperCorrect4500Mass5500_->Fill(fabs(costheta),newweight);
    }

  /************************************************************************
   *
   * 2) tanphi = (Q^2 + Qt^2)^1/2 / Q (Dt dot R unit) /(Dt dot Qt unit)
   *
   ************************************************************************/
  TLorentzVector   Pbeam(0.0, 0.0,  4000., 4000.); // beam momentum in lab frame
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
}


void ZprimeMuMuPatMiniAodNewMC::PrintEventInformation(unsigned int runNumber, unsigned int lumiNumber, unsigned int eventNumber,
					  float vtxChi2, float vtxMass, float CosmicRejection)
{
  if(event_runNo == runNumber && event_lumi == lumiNumber && event_evtNo == eventNumber)
    {
      output_txt << event_runNo
                 << "        " << event_lumi
                 << "        " << event_evtNo
                 << "        " << vtxChi2
                 << "        " << vtxMass << endl;
      for(unsigned i=0; i<Mu_nbMuon->size(); i++){
        //if( fabs(Mu_etaTunePMuonBestTrack->at(i)) < EtaCut ){cout<<"[1] eta="<<Mu_etaTunePMuonBestTrack->at(i)<<endl;}
        cout<<"[0] phi="<<Mu_phiTunePMuonBestTrack->at(i)<<endl;
        cout<<"[1] eta="<<Mu_etaTunePMuonBestTrack->at(i)<<endl;
        if(Mu_isGlobalMuon->at(i) == 1) {cout<<"[2] isGlobal="<<Mu_isGlobalMuon->at(i)<<endl;}
        if(Mu_ptTunePMuonBestTrack->at(i) > 53.0) {cout<<"[3] ptcocktail="<<Mu_ptTunePMuonBestTrack->at(i)<<endl;}
        if(Mu_absdxyTunePMuonBestTrack->at(i) < 0.2) {cout<<"[4] absdxy="<<Mu_absdxyTunePMuonBestTrack->at(i)<<endl;}
        if(Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i) < 0.10) {cout<<"[5] trackiso="<<Mu_trackiso->at(i)/Mu_ptInnerTrack->at(i)<<endl;}
        if(Mu_numberOftrackerLayersWithMeasurement->at(i) > 5) {cout<<"[6] nbTrackerLayer="<<Mu_numberOftrackerLayersWithMeasurement->at(i)<<endl;}
        //if(Mu_numberOfValidPixelHits->at(i) > 0) {cout<<"[7] nbPixelHits="<<Mu_numberOfValidPixelHits->at(i)<<endl;}
        cout<<"[7] nbPixelHits (global tk) ="<<Mu_numberOfValidPixelHits->at(i)<<endl;
	//cout<<"[7bar] nbPixelHits (inner tk) ="<<Mu_innerTK_numberOfValidPixelHits->at(i)<<endl;
	if(Mu_numberOfValidMuonHits->at(i) > 0) {cout<<"[8] nbMuonHits="<<Mu_numberOfValidMuonHits->at(i)<<endl;}
        if(Mu_numberOfMatchedStations->at(i) > 1) {cout<<"[9] nbStation="<<Mu_numberOfMatchedStations->at(i)<<endl;}
        if(Mu_dPToverPTTunePMuonBestTrack->at(i) < 0.3) {cout<<"[10] DeltaPterror="<<Mu_dPToverPTTunePMuonBestTrack->at(i)<<endl;}
        cout<<"[11] Charge="<<Mu_chargeTunePMuonBestTrack->at(i)<<endl;
      }
      cout<<"[000] vtxMassMu="<<vtxMass<<endl;
      cout<<"[000] vtxChi2Mu="<<vtxChi2<<endl;
      cout<<"[000] CosAngle="<<CosmicRejection<<endl;
    }
}

//----------------------------------------------------
//                                                   -
//       Part for HLT & Reco Matching                -
//                                                   -
//----------------------------------------------------
bool ZprimeMuMuPatMiniAodNewMC::isPassHLT(){
  int nbMatch = 0;
  for(unsigned i=0; i<HLT_nb->size(); i++){
    if( (HLT_name->at(i) == "HLT_Mu50_v1" ||
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
      //cout<<"triggerName = "<<triggerName<<endl;
      nbMatch++;
    }
  }
  if(nbMatch>0) {
    return true;
  }
  else return false;
}

bool ZprimeMuMuPatMiniAodNewMC::RecoHLTMuonMatching(float RecoEta,float RecoPhi){
  int nbMatch = 0;
  float deltaR   = -10000.0;
  for(unsigned i=0; i<HLTObj_nbObj->size(); i++){
    //cout<<"[before]triggerName"<<HLTObj_collection->at(i) <<endl;
    if( HLTObj_collection->at(i) == "HLT_Mu50_v1" ||
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
	){
      //cout<<"[after]triggerName"<<HLTObj_collection->at(i) <<endl;
	deltaR   = delR(HLTObj_eta->at(i),HLTObj_phi->at(i),RecoEta,RecoPhi);
        //printf ("HLT_Eta = %f  HLT_Phi = %f recoEta = %f recoPhi = %f DelR_trigger = %f\n",HLTObj_eta->at(i),HLTObj_phi->at(i),RecoEta,RecoPhi,deltaR);
      if(fabs(deltaR)>RecoHLTMatchingDeltaRcut) continue;
      nbMatch++;
    }
  }
  if(nbMatch>0) return true;
  else return false;
}



//if( flavor==5 ) b jet
//if( flavor==4 ) c jets
//else light-flavor jet
bool ZprimeMuMuPatMiniAodNewMC::BTaggingDiscriminator(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2){
  int nbBTag  = 0;
  int nbMatch = 0;
  for(unsigned i=0; i<Nb_bDiscriminators->size(); i++){
    if( jet_btag_pt->at(i) < 35.0 || fabs(jet_btag_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if(deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if(deltaR2 < 0.5) continue;
    nbBTag++;
    if(nbBTag>1 && nbBTag<3 && jet_btag_pfCSVv2IVF_discriminator->at(i)>0.5426) {
      nbMatch++;
    }
    else continue;
  }
  if(nbMatch>0) return true;
  else return false;
}

bool ZprimeMuMuPatMiniAodNewMC::BTaggingDiscriminator2(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2){
  int nbBTag  = 0;
  int nbMatch = 0;
  for(unsigned i=0; i<Nb_bDiscriminators->size(); i++){
    if( jet_btag_pt->at(i) < 35.0 || fabs(jet_btag_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if(deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if(deltaR2 < 0.5) continue;
    nbBTag++;
    if(nbBTag>1 && nbBTag<3 && jet_btag_pfCSVv2IVF_discriminator->at(i)>0.800 &&jet_btag_pfCSVv2IVF_discriminator->at(i)<0.935) {
      nbMatch++;
    }
    else continue;
  }
  if(nbMatch>0) return true;
  else return false;
}

bool ZprimeMuMuPatMiniAodNewMC::BTaggingDiscriminator3(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2){
  int nbBTag  = 0;
  int nbMatch = 0;
  for(unsigned i=0; i<Nb_bDiscriminators->size(); i++){
    if( jet_btag_pt->at(i) < 35.0 || fabs(jet_btag_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if(deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if(deltaR2 < 0.5) continue;
    nbBTag++;
    if(nbBTag>1 && nbBTag<3 && jet_btag_pfCSVv2IVF_discriminator->at(i)>0.935) {
      nbMatch++;
    }
    else continue;
  }
  if(nbMatch>0) return true;
  else return false;
}


void ZprimeMuMuPatMiniAodNewMC::DrawBTaggingDiscriminator(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2){
  int nbBTag = 0;
  for(unsigned i=0; i<Nb_bDiscriminators->size(); i++){
    if( jet_btag_pt->at(i) < 35.0 || fabs(jet_btag_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if(deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_btag_eta->at(i),jet_btag_phi->at(i));
    if(deltaR2 < 0.5) continue;
    nbBTag++;
    h1_nbBTagStep1_->Fill(nbBTag);
    h1_jetBTagStep1_->Fill(jet_btag_pfCSVv2IVF_discriminator->at(i));
    if(nbBTag>1 && nbBTag<3)
      {
	h1_nbBTagStep2_->Fill(nbBTag);
	h1_jetBTagStep2_->Fill(jet_btag_pfCSVv2IVF_discriminator->at(i));
      }
    if(nbBTag>1 && nbBTag<3 && jet_btag_pfCSVv2IVF_discriminator->at(i)>0.5426)
      {
	h1_nbBTagStep3_->Fill(nbBTag);
	h1_jetBTagStep3_->Fill(jet_btag_pfCSVv2IVF_discriminator->at(i));
      }
  }
}

void ZprimeMuMuPatMiniAodNewMC::Boson(float Px1,float Py1,float Pz1,float En1,
					float Px2,float Py2,float Pz2,float En2,
					float ChargeEle1,float MetEt,float MetPx,
					float MetPy,float MetPz,float MetEn){

  TLorentzVector Ele;
  TLorentzVector Elebar;
  TLorentzVector MissingParticle;
  if(ChargeEle1<0) {
    Ele.SetPxPyPzE(Px1,Py1,Pz1,En1);
    Elebar.SetPxPyPzE(Px2,Py2,Pz2,En2);
  }
  if(ChargeEle1>0) {
    Ele.SetPxPyPzE(Px2,Py2,Pz2,En2);
    Elebar.SetPxPyPzE(Px1,Py1,Pz1,En1);
  }
  TLorentzVector Bosson(Ele+Elebar);
  MissingParticle.SetPxPyPzE(MetPx,MetPy,MetPz,MetEn);
  float a = Bosson.Angle(MissingParticle.Vect()); // get angle between v1 and v2
  //BosPt  = Bosson.Pt();
  //BosPhi = Bosson.Phi();
  //  if(Bosson.Pt()>60 && MetEt>0 && a>2.8 && (fabs(MetEt-Bosson.Pt())/Bosson.Pt())<0.4)
  if(Bosson.Pt()>60)
    {
      h1_BosPt_->Fill(Bosson.Pt());
      h1_BosPhi_->Fill(Bosson.Phi());
      h1_DeltaPhi_->Fill(a); //Fill(MetPhi-Bosson.Phi());
      h1_DeltaPtoverPt_->Fill(fabs(MetEt-Bosson.Pt())/Bosson.Pt());
      float Mt = sqrt(2.0*Bosson.Pt()*MetEt*(1.0-cos(a)));
      h1_Mt_->Fill(Mt);
      h1_MissingEt_->Fill(MetEt);
    }
}

bool ZprimeMuMuPatMiniAodNewMC::DiPFJet(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2){
  int nbBTag  = 0;
  int nbMatch = 0;
  for(unsigned i=0; i<jet_nb->size(); i++){
    if( jet_pt->at(i) < 35.0 || fabs(jet_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_eta->at(i),jet_phi->at(i));
    if(deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_eta->at(i),jet_phi->at(i));
    if(deltaR2 < 0.5) continue;
    nbBTag++;
    h1_NbPFjetsAll_->Fill(nbBTag);
    if(nbBTag>1 && nbBTag<3) {
      h1_NbPFjets2_->Fill(nbBTag);
      h1_ptPFjetsAll_->Fill(jet_pt->at(i));
      nbMatch++;
    }
    else continue;
  }
  if(nbMatch>0) return true;
  else return false;
}

bool ZprimeMuMuPatMiniAodNewMC::DiPFJetCut(float MuonEta1,float MuonPhi1,float MuonEta2,float MuonPhi2){
  int nbBTag  = 0;
  int nbMatch = 0;
  for(unsigned i=0; i<jet_nb->size(); i++){
    if( jet_pt->at(i) < 35.0 || fabs(jet_eta->at(i)) > 2.5 ) continue;
    float deltaR1 = delR(MuonEta1,MuonPhi1,jet_eta->at(i),jet_phi->at(i));
    if(deltaR1 < 0.5) continue;
    float deltaR2 = delR(MuonEta2,MuonPhi2,jet_eta->at(i),jet_phi->at(i));
    if(deltaR2 < 0.5) continue;
    nbBTag++;
    if(nbBTag>1 && nbBTag<3) {
      nbMatch++;
    }
    else continue;
  }
  if(nbMatch>0) return true;
  else return false;
}

double ZprimeMuMuPatMiniAodNewMC::MassCorrection(float M)
{
  float a = 1.06780e+00;
  float b = -1.20666e-04;
  float c = 3.22646e-08;
  float d = -3.94886e-12;
  double function = d*pow(M,3) + c*pow(M,2) + b*pow(M,1) + a;
  return function;
}



//============================ Method to di-jet (Barrel-Barrel) ========================
void ZprimeMuMuPatMiniAodNewMC::DrawDiJetMassBB(){
  float invmass = -10;
  for(unsigned jet1=0; jet1<Mu_nbMuon->size(); jet1++){
    if( Mu_isGlobalMuon->at(jet1) == 1 &&
	//	Mu_isMuonsCleaned->at(jet1) ==  Mu_isPF->at(jet1) &&
	Mu_ptTunePMuonBestTrack->at(jet1) > FR_Ptcut &&
        Mu_absdxyTunePMuonBestTrack->at(jet1) < 0.2 &&
        (Mu_trackiso->at(jet1)/Mu_ptInnerTrack->at(jet1)) < 0.10  &&
        Mu_numberOftrackerLayersWithMeasurement->at(jet1) > 5 &&
        Mu_numberOfValidPixelHits->at(jet1) > 0 &&
        Mu_numberOfValidMuonHits->at(jet1) > 0 &&
	Mu_passNewMatchedStationsCut->at(jet1) == 1 &&
        Mu_dPToverPTTunePMuonBestTrack->at(jet1) < 0.3 ) continue; //to get rid of real muons
    if( (Mu_isGlobalMuon->at(jet1) == 0 || Mu_isTrackerMuon->at(jet1) == 0) ||
	Mu_ptTunePMuonBestTrack->at(jet1) < FR_Ptcut ||
        Mu_absdxyTunePMuonBestTrack->at(jet1) > 0.2 ||
	Mu_absdzTunePMuonBestTrack->at(jet1) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet1) < 5  ||
        Mu_numberOfValidPixelHits->at(jet1) < 0  ) continue;
    for(unsigned jet2=0; jet2<Mu_nbMuon->size(); jet2++){
      if(jet2 == jet1) continue;
      if( Mu_chargeTunePMuonBestTrack->at(jet1)*Mu_chargeTunePMuonBestTrack->at(jet2) > 0 ) continue; //OS
      if( Mu_isGlobalMuon->at(jet2) == 1 &&
	  //Mu_isMuonsCleaned->at(jet2) ==  Mu_isPF->at(jet2) &&
	  Mu_ptTunePMuonBestTrack->at(jet2) > FR_Ptcut &&
	  Mu_absdxyTunePMuonBestTrack->at(jet2) < 0.2 &&
	  (Mu_trackiso->at(jet2)/Mu_ptInnerTrack->at(jet2)) < 0.10  &&
	  Mu_numberOftrackerLayersWithMeasurement->at(jet2) > 5 &&
	  Mu_numberOfValidPixelHits->at(jet2) > 0 &&
	  Mu_numberOfValidMuonHits->at(jet2) > 0 &&
	  Mu_passNewMatchedStationsCut->at(jet2) == 1 &&
	  Mu_dPToverPTTunePMuonBestTrack->at(jet2) < 0.3 ) continue; //to get rid of real muons
      if( (Mu_isGlobalMuon->at(jet2) == 0 || Mu_isTrackerMuon->at(jet2) == 0) ||
	  Mu_ptTunePMuonBestTrack->at(jet2) < FR_Ptcut ||
          Mu_absdxyTunePMuonBestTrack->at(jet2) > 0.2 ||
	  Mu_absdzTunePMuonBestTrack->at(jet2) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet2) < 5  ||
	  Mu_numberOfValidPixelHits->at(jet2) < 0  ) continue;
      if(fabs(Mu_etaTunePMuonBestTrack->at(jet1)) > 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet2)) > 1.2) continue;
      float mEl = 0.105658371500;
      float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(jet1),Mu_etaTunePMuonBestTrack->at(jet1),
			     Mu_phiTunePMuonBestTrack->at(jet1),mEl,
			     Mu_ptTunePMuonBestTrack->at(jet2),Mu_etaTunePMuonBestTrack->at(jet2),
			     Mu_phiTunePMuonBestTrack->at(jet2),mEl);
      //pick highest mass dijet
      if( MassDiJet > 60.0) {
	invmass = MassDiJet;
	bool fireHLT1 = isPassHLT();
	if(fireHLT1 == 0) continue;
	bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet1),
							     Mu_phiTunePMuonBestTrack->at(jet1));
	bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet2),
							     Mu_phiTunePMuonBestTrack->at(jet2));
	if(RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
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
void ZprimeMuMuPatMiniAodNewMC::DrawDiJetMassBE(){
  float invmass = -10;
  for(unsigned jet1=0; jet1<Mu_nbMuon->size(); jet1++){
    if( Mu_isGlobalMuon->at(jet1) == 1 &&
	//	Mu_isMuonsCleaned->at(jet1) ==  Mu_isPF->at(jet1) &&
        Mu_ptTunePMuonBestTrack->at(jet1) > FR_Ptcut &&
        Mu_absdxyTunePMuonBestTrack->at(jet1) < 0.2 &&
        (Mu_trackiso->at(jet1)/Mu_ptInnerTrack->at(jet1)) < 0.10  &&
        Mu_numberOftrackerLayersWithMeasurement->at(jet1) > 5 &&
        Mu_numberOfValidPixelHits->at(jet1) > 0 &&
        Mu_numberOfValidMuonHits->at(jet1) > 0 &&
	Mu_passNewMatchedStationsCut->at(jet1) == 1 &&
        Mu_dPToverPTTunePMuonBestTrack->at(jet1) < 0.3 ) continue; //to get rid of real muons
    if( (Mu_isGlobalMuon->at(jet1) == 0 || Mu_isTrackerMuon->at(jet1) == 0) ||
	fabs(Mu_etaTunePMuonBestTrack->at(jet1)) > 2.4 ||
	Mu_ptTunePMuonBestTrack->at(jet1) < FR_Ptcut ||
        Mu_absdxyTunePMuonBestTrack->at(jet1) > 0.2 ||
	Mu_absdzTunePMuonBestTrack->at(jet1) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet1) < 5  ||
        Mu_numberOfValidPixelHits->at(jet1) < 0  ) continue;
    for(unsigned jet2=0; jet2<Mu_nbMuon->size(); jet2++){
      if(jet2 == jet1) continue;
      if( Mu_chargeTunePMuonBestTrack->at(jet1)*Mu_chargeTunePMuonBestTrack->at(jet2) > 0 ) continue; //OS
      if( Mu_isGlobalMuon->at(jet2) == 1 &&
	  //	  Mu_isMuonsCleaned->at(jet2) ==  Mu_isPF->at(jet2) &&
	  Mu_ptTunePMuonBestTrack->at(jet2) > FR_Ptcut &&
	  Mu_absdxyTunePMuonBestTrack->at(jet2) < 0.2 &&
	  (Mu_trackiso->at(jet2)/Mu_ptInnerTrack->at(jet2)) < 0.10  &&
	  Mu_numberOftrackerLayersWithMeasurement->at(jet2) > 5 &&
	  Mu_numberOfValidPixelHits->at(jet2) > 0 &&
	  Mu_numberOfValidMuonHits->at(jet2) > 0 &&
	  Mu_passNewMatchedStationsCut->at(jet2) == 1 &&
	  Mu_dPToverPTTunePMuonBestTrack->at(jet2) < 0.3 ) continue; //to get rid of real muons
      if( (Mu_isGlobalMuon->at(jet2) == 0 || Mu_isTrackerMuon->at(jet2) == 0) ||
	  fabs(Mu_etaTunePMuonBestTrack->at(jet2)) > 2.4 ||
	  Mu_ptTunePMuonBestTrack->at(jet2) < FR_Ptcut ||
          Mu_absdxyTunePMuonBestTrack->at(jet2) > 0.2 ||
	  Mu_absdzTunePMuonBestTrack->at(jet2) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet2) < 5  ||
	  Mu_numberOfValidPixelHits->at(jet2) < 0  ) continue;
      if(fabs(Mu_etaTunePMuonBestTrack->at(jet1)) > 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet2)) < 1.2) continue;
      float mEl = 0.105658371500;
      float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(jet1),Mu_etaTunePMuonBestTrack->at(jet1),
			     Mu_phiTunePMuonBestTrack->at(jet1),mEl,
			     Mu_ptTunePMuonBestTrack->at(jet2),Mu_etaTunePMuonBestTrack->at(jet2),
			     Mu_phiTunePMuonBestTrack->at(jet2),mEl);
      //pick highest mass dijet
      if( MassDiJet > 60.0) {
	invmass = MassDiJet;
	bool fireHLT1 = isPassHLT();
	if(fireHLT1 == 0) continue;
	bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet1),
							     Mu_phiTunePMuonBestTrack->at(jet1));
	bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet2),
							     Mu_phiTunePMuonBestTrack->at(jet2));
	if(RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
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
void ZprimeMuMuPatMiniAodNewMC::DrawDiJetMassEE(){
  float invmass = -10;
  for(unsigned jet1=0; jet1<Mu_nbMuon->size(); jet1++){
    if( Mu_isGlobalMuon->at(jet1) == 1 &&
	//	Mu_isMuonsCleaned->at(jet1) ==  Mu_isPF->at(jet1) &&
        Mu_ptTunePMuonBestTrack->at(jet1) > FR_Ptcut &&
        Mu_absdxyTunePMuonBestTrack->at(jet1) < 0.2 &&
        (Mu_trackiso->at(jet1)/Mu_ptInnerTrack->at(jet1)) < 0.10  &&
        Mu_numberOftrackerLayersWithMeasurement->at(jet1) > 5 &&
        Mu_numberOfValidPixelHits->at(jet1) > 0 &&
        Mu_numberOfValidMuonHits->at(jet1) > 0 &&
	Mu_passNewMatchedStationsCut->at(jet1) == 1 &&
        Mu_dPToverPTTunePMuonBestTrack->at(jet1) < 0.3 ) continue; //to get rid of real muons
    if( (Mu_isGlobalMuon->at(jet1) == 0 || Mu_isTrackerMuon->at(jet1) == 0) ||
	fabs(Mu_etaTunePMuonBestTrack->at(jet1)) > 2.4 ||
	Mu_ptTunePMuonBestTrack->at(jet1) < FR_Ptcut ||
        Mu_absdxyTunePMuonBestTrack->at(jet1) > 0.2 ||
	Mu_absdzTunePMuonBestTrack->at(jet1) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet1) < 5  ||
        Mu_numberOfValidPixelHits->at(jet1) < 0  ) continue;
    for(unsigned jet2=0; jet2<Mu_nbMuon->size(); jet2++){
      if(jet2 == jet1) continue;
      if( Mu_chargeTunePMuonBestTrack->at(jet1)*Mu_chargeTunePMuonBestTrack->at(jet2) > 0 ) continue; //OS
      if( Mu_isGlobalMuon->at(jet2) == 1 &&
	  // Mu_isMuonsCleaned->at(jet2) ==  Mu_isPF->at(jet2) &&
	  Mu_ptTunePMuonBestTrack->at(jet2) > FR_Ptcut &&
	  Mu_absdxyTunePMuonBestTrack->at(jet2) < 0.2 &&
	  (Mu_trackiso->at(jet2)/Mu_ptInnerTrack->at(jet2)) < 0.10  &&
	  Mu_numberOftrackerLayersWithMeasurement->at(jet2) > 5 &&
	  Mu_numberOfValidPixelHits->at(jet2) > 0 &&
	  Mu_numberOfValidMuonHits->at(jet2) > 0 &&
	  Mu_passNewMatchedStationsCut->at(jet2) == 1 &&
	  Mu_dPToverPTTunePMuonBestTrack->at(jet2) < 0.3 ) continue; //to get rid of real muons
      if( (Mu_isGlobalMuon->at(jet2) == 0 || Mu_isTrackerMuon->at(jet2) == 0) ||
	  fabs(Mu_etaTunePMuonBestTrack->at(jet2)) > 2.4 ||
	  Mu_ptTunePMuonBestTrack->at(jet2) < FR_Ptcut ||
          Mu_absdxyTunePMuonBestTrack->at(jet2) > 0.2 ||
	  Mu_absdzTunePMuonBestTrack->at(jet2) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet2) < 5  ||
	  Mu_numberOfValidPixelHits->at(jet2) < 0  ) continue;
      if(fabs(Mu_etaTunePMuonBestTrack->at(jet1)) < 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet2)) < 1.2) continue;
      float mEl = 0.105658371500;
      float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(jet1),Mu_etaTunePMuonBestTrack->at(jet1),
			     Mu_phiTunePMuonBestTrack->at(jet1),mEl,
			     Mu_ptTunePMuonBestTrack->at(jet2),Mu_etaTunePMuonBestTrack->at(jet2),
			     Mu_phiTunePMuonBestTrack->at(jet2),mEl);
      //pick highest mass dijet
      if( MassDiJet > 60.0) {
	invmass = MassDiJet;
	bool fireHLT1 = isPassHLT();
	if(fireHLT1 == 0) continue;
	bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet1),
							     Mu_phiTunePMuonBestTrack->at(jet1));
	bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet2),
							     Mu_phiTunePMuonBestTrack->at(jet2));
	if(RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
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
void ZprimeMuMuPatMiniAodNewMC::DrawWJetsMassBB(){
  float invmass = -10;
  for(unsigned muon=0; muon<Mu_nbMuon->size(); muon++){
    if( Mu_isGlobalMuon->at(muon) == 1 &&
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
      for(unsigned jet=0; jet<Mu_nbMuon->size(); jet++){
	if(jet == muon) continue;
	if( Mu_chargeTunePMuonBestTrack->at(muon)*Mu_chargeTunePMuonBestTrack->at(jet) > 0 ) continue; //OS
	if( Mu_isGlobalMuon->at(jet) == 1 &&
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
	if( (Mu_isGlobalMuon->at(jet) == 0 || Mu_isTrackerMuon->at(jet) == 0) ||
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 2.4 ||
	    Mu_ptTunePMuonBestTrack->at(jet) < FR_Ptcut ||
	    Mu_absdxyTunePMuonBestTrack->at(jet) > 0.2 ||
	    Mu_absdzTunePMuonBestTrack->at(jet) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet) < 5  ||
	    Mu_numberOfValidPixelHits->at(jet) < 0  ) continue;
	float deltaR1 = delR(Mu_etaTunePMuonBestTrack->at(muon),Mu_phiTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(jet),Mu_phiTunePMuonBestTrack->at(jet));
        if(deltaR1 < 0.5) continue;
	if(fabs(Mu_etaTunePMuonBestTrack->at(muon)) > 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 1.2) continue;
	float mEl = 0.105658371500;
	float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(muon),
			       Mu_phiTunePMuonBestTrack->at(muon),mEl,
			       Mu_ptTunePMuonBestTrack->at(jet),Mu_etaTunePMuonBestTrack->at(jet),
			       Mu_phiTunePMuonBestTrack->at(jet),mEl);
	//pick highest mass dijet
	if( MassDiJet > 60.0) {
	  invmass = MassDiJet;
	  bool fireHLT1 = isPassHLT();
	  if(fireHLT1 == 0) continue;
	  bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(muon),
							       Mu_phiTunePMuonBestTrack->at(muon));
	  bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet),
							       Mu_phiTunePMuonBestTrack->at(jet));
	  if(RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
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
void ZprimeMuMuPatMiniAodNewMC::DrawWJetsMassBE1(){
  float invmass = -10;
  for(unsigned muon=0; muon<Mu_nbMuon->size(); muon++){
    if( Mu_isGlobalMuon->at(muon) == 1 &&
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
      for(unsigned jet=0; jet<Mu_nbMuon->size(); jet++){
	if(jet == muon) continue;
	if( Mu_chargeTunePMuonBestTrack->at(muon)*Mu_chargeTunePMuonBestTrack->at(jet) > 0 ) continue; //OS
	if( Mu_isGlobalMuon->at(jet) == 1 &&
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
	if( (Mu_isGlobalMuon->at(jet) == 0 || Mu_isTrackerMuon->at(jet) == 0) ||
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 2.4 ||
	    Mu_ptTunePMuonBestTrack->at(jet) < FR_Ptcut ||
	    Mu_absdxyTunePMuonBestTrack->at(jet) > 0.2 ||
	    Mu_absdzTunePMuonBestTrack->at(jet) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet) < 5  ||
	    Mu_numberOfValidPixelHits->at(jet) < 0  ) continue;
	float deltaR1 = delR(Mu_etaTunePMuonBestTrack->at(muon),Mu_phiTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(jet),Mu_phiTunePMuonBestTrack->at(jet));
        if(deltaR1 < 0.5) continue;
	if(fabs(Mu_etaTunePMuonBestTrack->at(muon)) > 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet)) < 1.2) continue;
	float mEl = 0.105658371500;
	float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(muon),
			       Mu_phiTunePMuonBestTrack->at(muon),mEl,
			       Mu_ptTunePMuonBestTrack->at(jet),Mu_etaTunePMuonBestTrack->at(jet),
			       Mu_phiTunePMuonBestTrack->at(jet),mEl);
	//pick highest mass dijet
	if( MassDiJet > 60.0) {
	  invmass = MassDiJet;
	  bool fireHLT1 = isPassHLT();
	  if(fireHLT1 == 0) continue;
	  bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(muon),
							       Mu_phiTunePMuonBestTrack->at(muon));
	  bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet),
							       Mu_phiTunePMuonBestTrack->at(jet));
	  if(RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
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

void ZprimeMuMuPatMiniAodNewMC::DrawWJetsMassBE2(){
  float invmass = -10;
  for(unsigned muon=0; muon<Mu_nbMuon->size(); muon++){
    if( Mu_isGlobalMuon->at(muon) == 1 &&
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
      for(unsigned jet=0; jet<Mu_nbMuon->size(); jet++){
	if(jet == muon) continue;
	if( Mu_chargeTunePMuonBestTrack->at(muon)*Mu_chargeTunePMuonBestTrack->at(jet) > 0 ) continue; //OS
	if( Mu_isGlobalMuon->at(jet) == 1 &&
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
        if(deltaR1 < 0.5) continue;
	if( (Mu_isGlobalMuon->at(jet) == 0 || Mu_isTrackerMuon->at(jet) == 0) ||
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 2.4 ||
	    Mu_ptTunePMuonBestTrack->at(jet) < FR_Ptcut ||
	    Mu_absdxyTunePMuonBestTrack->at(jet) > 0.2 ||
	    Mu_absdzTunePMuonBestTrack->at(jet) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet) < 5  ||
	    Mu_numberOfValidPixelHits->at(jet) < 0  ) continue;
	if(fabs(Mu_etaTunePMuonBestTrack->at(muon)) < 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 1.2) continue;
	float mEl = 0.105658371500;
	float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(muon),
			       Mu_phiTunePMuonBestTrack->at(muon),mEl,
			       Mu_ptTunePMuonBestTrack->at(jet),Mu_etaTunePMuonBestTrack->at(jet),
			       Mu_phiTunePMuonBestTrack->at(jet),mEl);
	//pick highest mass dijet
	if( MassDiJet > 60.0) {
	  invmass = MassDiJet;
	  bool fireHLT1 = isPassHLT();
	  if(fireHLT1 == 0) continue;
	  bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(muon),
							       Mu_phiTunePMuonBestTrack->at(muon));
	  bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet),
							       Mu_phiTunePMuonBestTrack->at(jet));
	  if(RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
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
void ZprimeMuMuPatMiniAodNewMC::DrawWJetsMassEE(){
  float invmass = -10;
  for(unsigned muon=0; muon<Mu_nbMuon->size(); muon++){
    if( Mu_isGlobalMuon->at(muon) == 1 &&
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
      for(unsigned jet=0; jet<Mu_nbMuon->size(); jet++){
	if(jet == muon) continue;
	if( Mu_chargeTunePMuonBestTrack->at(muon)*Mu_chargeTunePMuonBestTrack->at(jet) > 0 ) continue; //OS
	if( Mu_isGlobalMuon->at(jet) == 1 &&
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
	if( (Mu_isGlobalMuon->at(jet) == 0 || Mu_isTrackerMuon->at(jet) == 0) ||
	    fabs(Mu_etaTunePMuonBestTrack->at(jet)) > 2.4 ||
	    Mu_ptTunePMuonBestTrack->at(jet) < FR_Ptcut ||
	    Mu_absdxyTunePMuonBestTrack->at(jet) > 0.2 ||
	    Mu_absdzTunePMuonBestTrack->at(jet) > 1.0 || Mu_numberOftrackerLayersWithMeasurement->at(jet) < 5  ||
	    Mu_numberOfValidPixelHits->at(jet) < 0  ) continue;
	float deltaR1 = delR(Mu_etaTunePMuonBestTrack->at(muon),Mu_phiTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(jet),Mu_phiTunePMuonBestTrack->at(jet));
        if(deltaR1 < 0.5) continue;
	if(fabs(Mu_etaTunePMuonBestTrack->at(muon)) < 1.2 || fabs(Mu_etaTunePMuonBestTrack->at(jet)) < 1.2) continue;
	float mEl = 0.105658371500;
	float MassDiJet = Mass(Mu_ptTunePMuonBestTrack->at(muon),Mu_etaTunePMuonBestTrack->at(muon),
			       Mu_phiTunePMuonBestTrack->at(muon),mEl,
			       Mu_ptTunePMuonBestTrack->at(jet),Mu_etaTunePMuonBestTrack->at(jet),
			       Mu_phiTunePMuonBestTrack->at(jet),mEl);
	//pick highest mass dijet
	if( MassDiJet > 60.0) {
	  invmass = MassDiJet;
	  bool fireHLT1 = isPassHLT();
	  if(fireHLT1 == 0) continue;
	  bool RecoMuon1MatchingWithHLT1 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(muon),
							       Mu_phiTunePMuonBestTrack->at(muon));
	  bool RecoMuon2MatchingWithHLT2 = RecoHLTMuonMatching(Mu_etaTunePMuonBestTrack->at(jet),
							       Mu_phiTunePMuonBestTrack->at(jet));
	  if(RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1) {
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
float ZprimeMuMuPatMiniAodNewMC::FRweight(float eta, float pt){
  float FR = 0.0;
  parEB1 = 4.39759e-01;
  parEB2 = -6.62025e+02;
  parEB3 = 1.61874e+03;
  parEB4 = 0.1733;
  parEE1 = 6.62323e-01;
  parEE2 = -4.90668e+02;
  parEE3 = 7.96427e+02;
  parEE4 = 0.3570;
  if(fabs(eta)<1.2 && pt>50.0 && pt<700.0)
    {
      FR = parEB1 + parEB2 / (parEB3 + pt );
    }
  else if(fabs(eta)<1.2 && pt>700.0)
    {
      FR = parEB4;
    }
  else if(fabs(eta)>1.2 && pt>50.0 && pt<700.0)
    {
      FR = parEE1 + parEE2 / (parEE3 + pt );
    }
  else if(fabs(eta)>1.2 && fabs(eta)<2.4 && pt>700.0)
    {
      FR = parEE4;
    }
  else {
    cout<<"out of FR range"<<endl;
  }
  return (FR/(1-FR));
}
*/
/*
//================= Method to select first high pt muon ==============
float ZprimeMuMuPatMiniAodNewMC::FRweight(float eta, float pt){
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
  if(fabs(eta)<1.2 && pt>FR_Ptcut && pt<=600.0)
    {
      FR = parEB1 + parEB2*pow(pt,1) + parEB3*pow(pt,2) + parEB4*pow(pt,3);
    }
  else if(fabs(eta)<1.2 && pt>600.0)
    {
      FR = parEB5;
    }
  else if(fabs(eta)>1.2  && fabs(eta)<2.4 && pt>FR_Ptcut && pt<=600.0)
    {
      FR = parEE1 + parEE2*pow(pt,1) + parEE3*pow(pt,2) + parEE4*pow(pt,3);
    }
  else if(fabs(eta)>1.2 && fabs(eta)<2.4 && pt>600.0)
    {
      FR = parEE5;
    }
  else {
    cout<<"out of FR range"<<endl;
  }
  return (FR/(1-FR));
}
*/


//================= Method to select first high pt muon ==============

float ZprimeMuMuPatMiniAodNewMC::FRweight(float eta, float pt){
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

  if(fabs(eta)<1.2 && pt>FR_Ptcut && pt<=200.0)
    {
      FR = parEB1 + parEB2*pow(pt,1) + parEB3*pow(pt,2);
    }
  else if(fabs(eta)<1.2 && pt>200 && pt<=800.0)
    {
      FR = parEB4 + parEB5 / (parEB6 + pt );
    }
  else if(fabs(eta)<1.2 && pt>800.0)
    {
      FR = parEB7;
    }
  else if(fabs(eta)>1.2  && fabs(eta)<2.4 && pt>FR_Ptcut && pt<=250.0)
    {
      FR = parEE1 + parEE2*pow(pt,1) + parEE3*pow(pt,2);
    }
  else if(fabs(eta)>1.2 && fabs(eta)<2.4 && pt>250.0)
    {
      FR = parEE4;
    }
  else {
    cout<<"out of FR range"<<endl;
  }
  return (FR/(1-FR));
}

/*
float ZprimeMuMuPatMiniAodNewMC::FRweight(float eta, float pt){
  float FR = 0.0;
  parEB1  = 1.09268e-01;
  parEB2  = -1.17369e-03;
  parEB3  = 5.01282e-06;
  parEB4  = -1.94865e+00;
  parEB5  = 6.59464e+03 ;
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

  if(fabs(eta)<1.2 && pt>FR_Ptcut && pt<=250.0)
    {
      FR = parEB1 + parEB2*pow(pt,1) + parEB3*pow(pt,2);
    }
  else if(fabs(eta)<1.2 && pt>250)
    {
      FR = parEB4+(parEB5/(pt+parEB6))+(parEB7/(pt*pt+parEB8))+(parEB9/(pt*pt*pt+parEB10));
    }
  else if(fabs(eta)>1.2  && fabs(eta)<2.4 && pt>FR_Ptcut && pt<=250.0)
    {
      FR = parEE1 + parEE2*pow(pt,1) + parEE3*pow(pt,2);
    }
  else if(fabs(eta)>1.2 && fabs(eta)<2.4 && pt>250.0)
    {
      FR = parEE4+(parEE5/(pt+parEE6))+(parEE7/(pt*pt+parEE8))+(parEE9/(pt*pt*pt+parEE10));
    }
  else {
    cout<<"out of FR range"<<endl;
  }
  return (FR/(1-FR));
}
*/
