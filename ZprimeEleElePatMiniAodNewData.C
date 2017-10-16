//==============================================================
//          Analysis code for Z' boson to e+ e- analysis       =  
//          In this code we select the HEEP events             = 
//          To run over MINIAOD MC with fixed trigger          = 
//                  Author:  Sherif Elgammal                   = 
//                                                             = 
//                       01/01/2016                            = 
//==============================================================
#define ZprimeEleElePatMiniAodNewData_cxx
#include "ZprimeEleElePatMiniAodNewData.h"
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

void ZprimeEleElePatMiniAodNewData::Loop()
{
time_t start,end;
double dif;
time (&start);
FILE * pFile;
   pFile = fopen ("myfile.txt","w");
   //values needed for AccepXeff study
   NbGen = 0;
   NbReco= 0;
   int binMass   = 100; //100; //100; //100; //100; //100; //100; //100; //100; //100; //10000;
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
   double muon_mass = 0.1056583;
   weight=1.;
   //if( DATA_type=="2015") weight=1.;
   TFile *output = new TFile("Data_B_rereco_json_analysis.root","recreate");
   // Enable Sumw2 for histograms as we'll be normalizing them
   // TH1::SetDefaultSumw2();
   //================================================================================== 
   //                                                                                 =
   //             Start the histograms for CollinSoper CMF                            =
   //                                                                                 =
   //==================================================================================
   NbFireHLT = 0;
   int NbBins   = 10;
   float MinBin = -1.0;
   float MaxBin =  1.0;
   h1_CosAngleCollinSoperCorrect60Mass120_    = new TH1F("CosAngleCollinSoperCorrect60Mass120","",NbBins,MinBin,MaxBin);
   h1_CosAngleCollinSoperCorrect120Mass300_   = new TH1F("CosAngleCollinSoperCorrect120Mass300","",NbBins,MinBin,MaxBin);
   h1_CosAngleCollinSoperCorrect300Mass700_   = new TH1F("CosAngleCollinSoperCorrect300Mass700","",NbBins,MinBin,MaxBin);
   h1_CosAngleCollinSoperCorrect700Mass3000_  = new TH1F("CosAngleCollinSoperCorrect700Mass3000","",NbBins,MinBin,MaxBin);
   h1_CosAngleCollinSoperCorrect4900Mass5100_ = new TH1F("CosAngleCollinSoperCorrect4900Mass5100","",NbBins,MinBin,MaxBin); 
   h1_absCosAngleCollinSoperCorrect4500Mass5500_ = new TH1F("absCosAngleCollinSoperCorrect4500Mass5500","",5,0.0,1.0);
   //==================================================================================
   //================================================================================== 
   //                                                                                 =
   //            Start the histograms for the mass of Z                               =  
   //                                                                                 =
   //==================================================================================
   // Build the histo with constant log bin width
   const int NMBINS = 100;
   const double MMIN = 60., MMAX = 6000.;
   double logMbins[NMBINS+1];
   float binNormNr=0.;
   for (int ibin = 0; ibin <= NMBINS; ibin++) {
     logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
     cout << logMbins[ibin] << endl;
   }
   h1_ZprimeRecomassBinWidth_        = new TH1F("ZprimeRecomassBinWidth","ZprimeRecomassBinWidth",NMBINS, logMbins);
   h1_ZprimeRecomassBinWidthEE_      = new TH1F("ZprimeRecomassBinWidthEE","",NMBINS, logMbins);
   h1_ZprimeRecomassBinWidthBB_      = new TH1F("ZprimeRecomassBinWidthBB","",NMBINS, logMbins);
   h1_ZprimeRecomassBinWidthBE_      = new TH1F("ZprimeRecomassBinWidthBE","",NMBINS, logMbins);
   h1_ZprimeRecomass60to120EE_       = new TH1F("ZprimeRecomass60to120EE","",60,60.0,120.0);
   h1_ZprimeRecomass60to120BB_       = new TH1F("ZprimeRecomass60to120BB","",60,60.0,120.0);
   h1_ZprimeRecomass60to120BE_       = new TH1F("ZprimeRecomass60to120BE","",60,60.0,120.0);
   h1_ZprimeRecomass60to120_         = new TH1F("ZprimeRecomass60to120","",60,60.0,120.0);
   // Book txt file for candidate events
   Char_t txtOUT[500];
   sprintf(txtOUT,"ZprimeToMuMu_13TeV_cand.txt");
   output_txt.open(txtOUT);
   output_txt << "CANDIDATES Events:" << endl;
   Char_t outform[20000];
   sprintf (outform,"run: lumi: event: dil_mass: pTmu1: pTmu2: Etamu1: Etamu2:");
   output_txt  << outform << endl; 
   //==================================================================================  
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if( inputfile.Contains("amcatnlo") && !(inputfile.Contains("DYJetsToTauTau"))) {      
       //cout << "Reweighting sample of amcatnlo with weight= " << MC_weighting->at(0)<< endl;
       //weight=weight*MC_weighting->at(0);
     }

     cout<<"=======> jentry = "<<jentry<< 
       "=======> Evt = "<<event_evtNo<< 
       "=======> Run = "<<event_runNo<< 
       "=======> Lumi = "<<event_lumi<< 
       "=======> bunch = "<<event_bunch<<endl; 

     //=====================================================
     //                                                    =
     //  Calling methods to get the generated information  =
     //                                                    =
     //=====================================================
     /* bool firstGenMu  = SelectFirstGenMu(genET1,genPhi1,genEta1,genEn1,genID1,genStat1,flag1);
     bool secondGenMu = SelectSecondGenMu(flag1,genET1,genET2,genPhi2,genEta2,genEn2,genID2,genStat2);
     h1_ptBeforeTrigger_->Fill(genET1);
     bool fireHLT = isPassHLT();
     if(fireHLT == 1) { 
       NbFireHLT++;
       h_NbFireHLT->Fill(NbFireHLT);
       h1_ptAfterTrigger_->Fill(genET1);
     }
     if(firstGenMu == 0 || secondGenMu == 0) continue;
     MassGen = Mass(genET1,genEta1,genPhi1,genEn1,genET2,genEta2,genPhi2,genEn2);
     PlotGenInfo(MassGen,genEta1,genEta2,genET1,genET2,genEn1,genEn2); */
     //=========================================================
     //                                                        =
     // Calling methods to get events with 2 muons passing ID  =   
     //                                                        =
     //=========================================================
     bool firstMuFinal  = SelectFirstEle(Etele1,Enele1,EtaTrakele1,PhiTrakele1,Chargeele1,
					 EtaSCele1,PhiSCele1,flagmu1);
     bool secondMuFinal = SelectSecondEle(Chargeele1,flagmu1,Etele1,EtaTrakele1,PhiTrakele1,
					  Etele2,Enele2,EtaTrakele2,PhiTrakele2,Chargeele2,
					  EtaSCele2,PhiSCele2);
     //=========================================================
     //        call the method for N-1 plots                   =
     //                                                        =
     //=========================================================
     //cout << "firstMu= " << firstMuFinal << " " << "secondMu= " << secondMuFinal << endl;
     if(firstMuFinal == 0 || secondMuFinal == 0) continue;
     DiEleMass = Mass(Etele1,EtaTrakele1,PhiTrakele1,Enele1,
		      Etele2,EtaTrakele2,PhiTrakele2,Enele2);
     if(DiEleMass<60) continue;
     //=========================================================
     //        start doing matching between reco & HLT         =
     //                                                        =
     //=========================================================
     /*bool fireHLT2 = isPassHLT();
     if(fireHLT2 == 0) continue; 
     bool RecoMuon1MatchingWithHLT1 = RecoHLTEleMatching(EtaSCele1,PhiSCele1);
     bool RecoMuon2MatchingWithHLT2 = RecoHLTEleMatching(EtaSCele2,PhiSCele2);
     if(RecoMuon1MatchingWithHLT1==1 || RecoMuon2MatchingWithHLT2==1)
     {*/
     PlotRecoInfo(DiEleMass,EtaSCele1,EtaSCele2);
     CosThetaCollinSoper(Etele1,EtaSCele1,PhiSCele1,Enele1,
			 Etele2,EtaSCele2,PhiSCele2,Enele2,
			 Chargeele1,DiEleMass);
     //}
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
float ZprimeEleElePatMiniAodNewData::delR(float eta1,float phi1,float eta2,float phi2){
  float mpi=3.14;
  float dp=std::abs(phi1-phi2);
  if (dp>mpi) dp-=float(2*mpi);
  return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
}

bool ZprimeEleElePatMiniAodNewData::GenRecoMatchEle(float RecoEta1,float RecoPhi1){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  for(unsigned i=0; i<iGen->size(); i++){
    float deltaR1   = delR(RecoEta1,RecoPhi1,etaGen->at(i),phiGen->at(i));
    if( fabs(idGen->at(i)) != 11 ) continue;
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

//============================ Method to select first high pt ele ========================
bool ZprimeEleElePatMiniAodNewData::SelectFirstEle(float &ETele1,float &Enele1,float &Etaele1,
						   float &Phiele1,int &ChargeEle1,float &EtaSCele1,
						   float &PhiSCele1,unsigned &FlagEle1){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  float highestpt=-999.;
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    float ET = Ele_EtFromCaloEn->at(i);
    //Barrel
    if ( ET > 35 && fabs(Ele_etaSC->at(i)) < 1.4442 &&
	 Ele_isEcalDrivenSeed->at(i) &&
	 fabs(Ele_deltaEtaInSeedCluster->at(i)) < 0.004 &&
	 fabs(Ele_deltaPhiInSeedCluster->at(i)) < 0.06 &&
	 Ele_hadronicOverEm->at(i) < (0.05 + 1.0/ Ele_energySC->at(i)) &&
	 (Ele_e1x5Over5x5Full5x5->at(i) > 0.83 || Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94) &&
	 Ele_nbOfMissingHits->at(i) < 2 &&
	 fabs(Ele_dxy->at(i)) < 0.02 &&
	 Ele_dr03EcalRecHitSumEt->at(i) + Ele_dr03HcalDepth1TowerSumEt->at(i) < 2 + 0.03 * ET + 0.28 * Ele_rhoIso->at(i) &&
	 Ele_dr03TkSumPt_corrected->at(i) < 5.0 ){
      if (ET>highestpt) {
	//bool GenRecoMatch1 = GenRecoMatchEle(Ele_etaSC->at(i),Ele_phiSC->at(i));
	//if(GenRecoMatch1 == 0) continue;
	highestpt=ET;
	iflag  = i;
	NbHEEPele ++;
      }
    }
    //endcap
    if ( ET > 35 && (fabs(Ele_etaSC->at(i)) > 1.566 && (abs(Ele_etaSC->at(i)) < 2.5) )&&
	 Ele_isEcalDrivenSeed->at(i) &&
	 fabs(Ele_deltaEtaInSeedCluster->at(i)) < 0.006 &&
	 fabs(Ele_deltaPhiInSeedCluster->at(i)) < 0.06 &&
	 Ele_hadronicOverEm->at(i) < (0.05 + 5.0/ Ele_energySC->at(i)) &&
	 Ele_sigmaIetaIetaFull5x5->at(i) <0.03 &&
	 Ele_nbOfMissingHits->at(i) < 2 &&
	 fabs(Ele_dxy->at(i)) < 0.05 &&
	 (( ET < 50 && Ele_dr03EcalRecHitSumEt->at(i) + Ele_dr03HcalDepth1TowerSumEt->at(i) < 2.5 + 0.28 * Ele_rhoIso->at(i)) ||
	  ( ET > 50 && Ele_dr03EcalRecHitSumEt->at(i) + Ele_dr03HcalDepth1TowerSumEt->at(i) < 2.5 + 0.03 * (ET-50) + 0.28 * Ele_rhoIso->at(i))) &&
	 Ele_dr03TkSumPt_corrected->at(i) < 5.0 ){
      if (ET>highestpt) {
	//bool GenRecoMatch1 = GenRecoMatchEle(Ele_etaSC->at(i),Ele_phiSC->at(i));
	//if(GenRecoMatch1 == 0) continue;
	highestpt=ET;
	iflag  = i;
	NbHEEPele ++;
      }
    }
    else continue;
  }
  if( NbHEEPele > 0 ){
    FlagEle1           = iflag;
    ETele1             = Ele_EtFromCaloEn->at(iflag);
    Enele1             = Ele_energySC->at(iflag);
    Etaele1            = Ele_etaTrack->at(iflag);
    Phiele1            = Ele_phiTrack->at(iflag);
    ChargeEle1         = Ele_charge->at(iflag);
    EtaSCele1          = Ele_etaSC->at(iflag);
    PhiSCele1          = Ele_phiSC->at(iflag);
    return true;  
  }
  else return false;
}
//============================ Method to select second high pt ele ========================
bool ZprimeEleElePatMiniAodNewData::SelectSecondEle(int ChargeEle1,unsigned FlagEle1,float ETele1,float Etaele1,float Phiele1,
						    float &ETele2,float &Enele2,float &Etaele2,float &Phiele2,int &ChargeEle2,
						    float &EtaSCele2,float &PhiSCele2){
  int NbHEEPele = 0;
  unsigned iflag = -10;
  float highestpt=-999.;
  for(unsigned i=0; i<Ele_nbElectrons->size(); i++){
    if(i == FlagEle1) continue;
    if(Ele_EtFromCaloEn->at(i) == ETele1) continue;
    if(Ele_etaSC->at(i) == Etaele1) continue;
    if(Ele_phiSC->at(i) == Phiele1) continue;
    float ET = Ele_EtFromCaloEn->at(i);
    //Barrel
    if ( ET > 35 && fabs(Ele_etaSC->at(i)) < 1.4442 &&
	 Ele_isEcalDrivenSeed->at(i) &&
	 fabs(Ele_deltaEtaInSeedCluster->at(i)) < 0.004 &&
	 fabs(Ele_deltaPhiInSeedCluster->at(i)) < 0.06 &&
	 Ele_hadronicOverEm->at(i) < (0.05 + 1.0/ Ele_energySC->at(i)) &&
	 (Ele_e1x5Over5x5Full5x5->at(i) > 0.83 || Ele_e2x5MaxOver5x5Full5x5->at(i) > 0.94) &&
	 Ele_nbOfMissingHits->at(i) < 2 &&
	 fabs(Ele_dxy->at(i)) < 0.02 &&
	 Ele_dr03EcalRecHitSumEt->at(i) + Ele_dr03HcalDepth1TowerSumEt->at(i) < 2 + 0.03 * ET + 0.28 * Ele_rhoIso->at(i) &&
	 Ele_dr03TkSumPt_corrected->at(i) < 5.0 ){
      if (ET>highestpt) {
	//bool GenRecoMatch1 = GenRecoMatchEle(Ele_etaSC->at(i),Ele_phiSC->at(i));
	//if(GenRecoMatch1 == 0) continue;
	highestpt=ET;
	iflag  = i;
	NbHEEPele ++;
      }
    }
    //endcap
    if ( ET > 35 && (fabs(Ele_etaSC->at(i)) > 1.566 && (abs(Ele_etaSC->at(i)) < 2.5) )&&
	 Ele_isEcalDrivenSeed->at(i) &&
	 fabs(Ele_deltaEtaInSeedCluster->at(i)) < 0.006 &&
	 fabs(Ele_deltaPhiInSeedCluster->at(i)) < 0.06 &&
	 Ele_hadronicOverEm->at(i) < (0.05 + 5.0/ Ele_energySC->at(i)) &&
	 Ele_sigmaIetaIetaFull5x5->at(i) <0.03 &&
	 Ele_nbOfMissingHits->at(i) < 2 &&
	 fabs(Ele_dxy->at(i)) < 0.05 &&
	 (( ET < 50 && Ele_dr03EcalRecHitSumEt->at(i) + Ele_dr03HcalDepth1TowerSumEt->at(i) < 2.5 + 0.28 * Ele_rhoIso->at(i)) ||
	  ( ET > 50 && Ele_dr03EcalRecHitSumEt->at(i) + Ele_dr03HcalDepth1TowerSumEt->at(i) < 2.5 + 0.03 * (ET-50) + 0.28 * Ele_rhoIso->at(i))) &&
	 Ele_dr03TkSumPt_corrected->at(i) < 5.0 ){
      if (ET>highestpt) {
	//bool GenRecoMatch1 = GenRecoMatchEle(Ele_etaSC->at(i),Ele_phiSC->at(i));
	//if(GenRecoMatch1 == 0) continue;
	highestpt=ET;
	iflag  = i;
	NbHEEPele ++;
      }
    }     
    else continue;
  }
  if( NbHEEPele > 0 ){
    ETele2             = Ele_EtFromCaloEn->at(iflag);
    Enele2             = Ele_energySC->at(iflag);
    Etaele2            = Ele_etaTrack->at(iflag);
    Phiele2            = Ele_phiTrack->at(iflag);
    ChargeEle2         = Ele_charge->at(iflag);
    EtaSCele2          = Ele_etaSC->at(iflag);
    PhiSCele2          = Ele_phiSC->at(iflag);
    return true;  
  }
  else return false;
}
void ZprimeEleElePatMiniAodNewData::PlotRecoInfo(float MassEle,float etaEle1,float etaEle2){
  //----------------------------------------------------------
  if (MassEle>1800.0) {
  output_txt << event_runNo
             << "   " << event_lumi
             << "       " << event_evtNo
             << "        " << MassEle
    //<< "        " << PtTunePEleBestTrack
    //<< "        " << PtTunePEleBestTrack2
             << "        " << etaEle1
             << "        " << etaEle2
             << endl;
  
   }
  //float weight2 = std::min(1.01696-7.73522E-5*MassEle+6.69239E-9*MassEle*MassEle,1);
  //----------------------------------------------------------
  if (!(inputfile.Contains("WW") && MassEle>2000.) ) {
    //BB
    if(fabs(etaEle1)<1.4442 && fabs(etaEle2)<1.4442)
      {
	h1_ZprimeRecomassBinWidthBB_->Fill(MassEle,weight);
	h1_ZprimeRecomass60to120BB_->Fill(MassEle,weight);
	h1_ZprimeRecomassBinWidth_->Fill(MassEle,weight);
	h1_ZprimeRecomass60to120_->Fill(MassEle,weight);
      }
      
    //BE
    if( fabs(etaEle1)<1.4442 && (etaEle2 > 1.566 && etaEle2 < 2.5) )
      {
	h1_ZprimeRecomassBinWidthBE_->Fill(MassEle,weight);
	h1_ZprimeRecomass60to120BE_->Fill(MassEle,weight);
	h1_ZprimeRecomassBinWidth_->Fill(MassEle,weight);
	h1_ZprimeRecomass60to120_->Fill(MassEle,weight);
      }
    if( fabs(etaEle2)<1.4442 && (etaEle1 > 1.566 && etaEle1 < 2.5) )
      {
	h1_ZprimeRecomassBinWidthBE_->Fill(MassEle,weight);
	h1_ZprimeRecomass60to120BE_->Fill(MassEle,weight);
	h1_ZprimeRecomassBinWidth_->Fill(MassEle,weight);
	h1_ZprimeRecomass60to120_->Fill(MassEle,weight);
      }
    
    //EE
    if( (etaEle1 > 1.566 && etaEle1 < 2.5) && (etaEle2 > 1.566 && etaEle2 < 2.5) )
      {
	h1_ZprimeRecomassBinWidthEE_->Fill(MassEle,weight);
	h1_ZprimeRecomass60to120EE_->Fill(MassEle,weight);
      }

  }
  //h1_ZprimeRecomass50_->Fill(MassEle);
  //h1_ZprimeRecomass20_->Fill(MassEle);
}
//===================== Methode to calculate the mass ========================
float ZprimeEleElePatMiniAodNewData::Mass(float Pt1,float Eta1,float Phi1,float En1,
					  float Pt2,float Eta2,float Phi2,float En2){
  float MuMuMass = 0.0;
  TLorentzVector Mu1;
  TLorentzVector Mu2;
  Mu1.SetPtEtaPhiE(Pt1,Eta1,Phi1,En1);
  Mu2.SetPtEtaPhiE(Pt2,Eta2,Phi2,En2);
  MuMuMass = (Mu1 + Mu2).M();
  return MuMuMass;
}


//----------------------------------------------------
//                                                   -
//       Part for Gen & Reco Matching                -
//                                                   -  
//----------------------------------------------------
//========================== Method to select firt Gen Ele =======================
bool ZprimeEleElePatMiniAodNewData::SelectFirstGenEle(float &ETEle1,float &PhiSCEle1,
						      float &EtaSCEle1,float &EnEle1,
						      int &IDele1,int &Statele1,
						      unsigned &GenFlag1){
  int NbHEEPele = 0;
  int iflag = -10;
  ETEle1 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    if( fabs(idGen->at(i)) != 11 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if( ptGen->at(i) > ETEle1) {
      ETEle1 = ptGen->at(i);
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if(NbHEEPele>0) {
    GenFlag1       = iflag;
    ETEle1          = ptGen->at(iflag);
    PhiSCEle1       = phiGen->at(iflag);
    EtaSCEle1       = etaGen->at(iflag);
    EnEle1          = EnergyGen->at(iflag);
    IDele1         = idGen->at(iflag);
    Statele1       = statusGen->at(iflag);
    return true;
  }         
  else return false;
}
//============================ Method to select second Gen Ele ========================
bool ZprimeEleElePatMiniAodNewData::SelectSecondGenEle(unsigned GenFlag1,float ETEle1,float &ETEle2,float &PhiSCEle2,
					     float &EtaSCEle2,float &EnEle2,int &IDele2,int &Statele2){
  int NbHEEPele = 0;
  int iflag = -10;
  ETEle2 = 0.0;
  for(unsigned i=0; i<iGen->size(); i++){
    if( fabs(idGen->at(i)) != 11 ) continue;
    if( statusGen->at(i) != 1 )  continue;
    if(i == GenFlag1) continue;
    if( fabs(ptGen->at(i) - ETEle1) <0.00001 ) continue;
    if( ptGen->at(i) > ETEle2) {
      ETEle2 = ptGen->at(i);
      iflag  = i;
      NbHEEPele ++;
    }
    else continue;
  }
  if(NbHEEPele>0) {
    ETEle2      = ptGen->at(iflag);
    PhiSCEle2   = phiGen->at(iflag);
    EtaSCEle2   = etaGen->at(iflag);
    EnEle2      = EnergyGen->at(iflag);
    IDele2     = idGen->at(iflag);
    Statele2   = statusGen->at(iflag); 
    return true;
  }
  else return false;
}


void ZprimeEleElePatMiniAodNewData::CosThetaCollinSoper(float Et1,float Eta1,float Phi1,float En1,
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
      h1_CosAngleCollinSoperCorrect60Mass120_->Fill(costheta,weight);
    }
  
  
  if( RecoMass > 120.0 && RecoMass < 300.0 ) 
    {
      h1_CosAngleCollinSoperCorrect120Mass300_->Fill(costheta,weight);
    }
  
  if( RecoMass > 300.0 && RecoMass < 700.0 ) 
    {
      h1_CosAngleCollinSoperCorrect300Mass700_->Fill(costheta,weight);
    }
  
  if( RecoMass > 700.0 && RecoMass < 3000.0 ) 
    {
      h1_CosAngleCollinSoperCorrect700Mass3000_->Fill(costheta,weight);
    }
  
  if( RecoMass > 4500.0 && RecoMass < 6000.0 ) 
    {
      h1_CosAngleCollinSoperCorrect4900Mass5100_->Fill(costheta,weight);
      h1_absCosAngleCollinSoperCorrect4500Mass5500_->Fill(fabs(costheta),weight);
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
  //h1_TanPhiCollinSoperCorrect_->Fill(tanphi,weight);
  /************************************************************************
   *
   * 3) sin2(theta) = Q^-2 Dt^2 - Q^-2 (Q^2 + Qt^2)^-1 * (Dt dot Qt)^2
   *
   ************************************************************************/
  //double dt_qt = D.X()*Q.X() + D.Y()*Q.Y();
  //double sin2theta = pow(D.Pt()/Q.Mag(), 2)
  //- 1.0/pow(Q.Mag(), 2)/(pow(Q.Mag(), 2) + pow(Q.Pt(), 2))*pow(dt_qt, 2);
  //h1_Sin2AngleCollinSoperCorrect_->Fill(sin2theta,weight);
}
//----------------------------------------------------
//                                                   -
//       Part for HLT & Reco Matching                -
//                                                   -  
//----------------------------------------------------
bool ZprimeEleElePatMiniAodNewData::isPassHLT(){ 
  int nbMatch = 0;
  for(unsigned i=0; i<HLT_nb->size(); i++){
    if( (HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v1" || 
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v2" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v3" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v4" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v5" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v6" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v7" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v8" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v9" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v10") && HLT_isaccept->at(i) == 1 ) {
      //cout<<"triggerName = "<<triggerName<<endl;
      nbMatch++;
    }
  }
  if(nbMatch>0) {
    return true;
  }
  else return false;
}

bool ZprimeEleElePatMiniAodNewData::RecoHLTEleMatching(float RecoEta,float RecoPhi){
  int nbMatch = 0;
  float deltaR   = -10000.0;
  for(unsigned i=0; i<HLTObj_nbObj->size(); i++){
    //cout<<"[before]triggerName"<<HLTObj_collection->at(i) <<endl;
    if( HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v1" || 
	HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v2" ||
	HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v3" || 
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v4" ||
	HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v5" || 
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v6" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v7" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v8" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v9" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v10" ){
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


