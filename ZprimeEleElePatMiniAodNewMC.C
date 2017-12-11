//==============================================================
//          Analysis code for Z' boson to e+ e- analysis       =
//          In this code we select the HEEP events             =
//          To run over MINIAOD MC with fixed trigger          =
//                  Author:  Sherif Elgammal                   =
//                                                             =
//                       11/03/2016                            =
//==============================================================
#define ZprimeEleElePatMiniAodNewMC_cxx
#include "ZprimeEleElePatMiniAodNewMC.h"
#include <math.h>
#include <array>
#include <iomanip>
#include "TStopwatch.h"
#include "TRandom3.h"
#include <time.h>
#include <algorithm>

// using namespace std;
#define PI 3.14159265
#define MUON_MASS 0.1056583
#define ELEC_MASS 0.000511

bool myfunction (int i,int j) { return (i < j); }
bool picklargemass (float lhs,float rhs) { return (lhs > rhs); }
TString inputfile;
float newweight = 1.;

void ZprimeEleElePatMiniAodNewMC::Loop(bool debug)
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
  deltaRcut = 0.10;
  RecoHLTMatchingDeltaRcut = 0.10;
  minMassCut = 50.0;
  maxMassCut = 4000.0;
  int ptBins = 40;
  float ptMin = 0.0;
  float ptMax = 400.0;
  ptEffCut = 3000.0;

  // m_weight = 1.;  // this is dumb, definitely don't want to reset the weight...
  if (DATA_type=="2015" || DATA_type=="2016" || DATA_type=="2017")
    m_weight = 1.;
  std::shared_ptr<TFile> output = std::make_shared<TFile>("ZprimeToEleEle_13TeV.root","recreate");
  // Enable Sumw2 for histograms as we'll be normalizing them
  TH1::SetDefaultSumw2();

  //==================================================================================
  //                                                                                 =
  //             Start the histograms for CollinSoper CMF                            =
  //                                                                                 =
  //==================================================================================
  NbFireHLT = 0;
  int NbBins   = 10;
  float MinBin = -1.0;
  float MaxBin =  1.0;
  // h1_ZprimeRecomass_                         = std::make_shared<TH1D>("ZprimeRecomass","",6000,0.0,6000.0);
  h1_ZprimeRecomass_                         = std::make_shared<TH1D>("ZprimeRecomass","",binMass,minMass,maxMass);
  h1_CosAngleCollinSoperCorrect60Mass120_    = std::make_shared<TH1D>("CosAngleCollinSoperCorrect60Mass120","",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect120Mass300_   = std::make_shared<TH1D>("CosAngleCollinSoperCorrect120Mass300","",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect300Mass700_   = std::make_shared<TH1D>("CosAngleCollinSoperCorrect300Mass700","",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect700Mass3000_  = std::make_shared<TH1D>("CosAngleCollinSoperCorrect700Mass3000","",NbBins,MinBin,MaxBin);
  h1_CosAngleCollinSoperCorrect4900Mass5100_ = std::make_shared<TH1D>("CosAngleCollinSoperCorrect4900Mass5100","",NbBins,MinBin,MaxBin);
  h1_absCosAngleCollinSoperCorrect4500Mass5500_ = std::make_shared<TH1D>("absCosAngleCollinSoperCorrect4500Mass5500","",5,0.0,1.0);

  double etaBins[] = {-2.5,-1.566,-1.4442,0.,1.4442,1.566,2.5};
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
    if (debug)
      std::cout << logMbins[ibin] << std::endl;
  }
  h1_ZprimeRecomassBinWidth_   = std::make_shared<TH1D>("ZprimeRecomassBinWidth",  "ZprimeRecomassBinWidth",NMBINS, logMbins);
  h1_ZprimeRecomassBinWidthEE_ = std::make_shared<TH1D>("ZprimeRecomassBinWidthEE","",NMBINS, logMbins);
  h1_ZprimeRecomassBinWidthBB_ = std::make_shared<TH1D>("ZprimeRecomassBinWidthBB","",NMBINS, logMbins);
  h1_ZprimeRecomassBinWidthBE_ = std::make_shared<TH1D>("ZprimeRecomassBinWidthBE","",NMBINS, logMbins);
  h1_ZprimeRecomass60to120EE_  = std::make_shared<TH1D>("ZprimeRecomass60to120EE", "",60,60.0,120.0);
  h1_ZprimeRecomass60to120BB_  = std::make_shared<TH1D>("ZprimeRecomass60to120BB", "",60,60.0,120.0);
  h1_ZprimeRecomass60to120BE_  = std::make_shared<TH1D>("ZprimeRecomass60to120BE", "",60,60.0,120.0);
  h1_ZprimeRecomass60to120_    = std::make_shared<TH1D>("ZprimeRecomass60to120",   "",60,60.0,120.0);
  // Book txt file for candidate events
  Char_t txtOUT[500];
  sprintf(txtOUT,"ZprimeToEleEle_13TeV_cand.txt");
  output_txt.open(txtOUT);
  output_txt << "CANDIDATES Events:" << std::endl;
  Char_t outform[20000];
  sprintf(outform,"run: lumi: event: dil_mass: pTele1: pTele2: Etaele1: Etaele2:");
  output_txt << outform << std::endl;

  TString inputfile=name;
  inputfile=name;
  std::cout << "Name of the input file is= " << inputfile.Data() << std::endl;
  std::cout << "Weight of the sample is= " << m_weight << std::endl;

  //==================================================================================
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
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
  for (Long64_t jentry=0; jentry < nentries;jentry++) {
    //Timing information
    if (jentry==0) {
      tsw.Start();
      std::cout << "." << std::flush;
    }
    if ((jentry*10)/nentries == tenpcount) {
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
      // } else if ( (jentry*100)/nentries == onepcount) {
      //   std::cout << ".";
      //   std::cout << std::flush;
      //   onepcount++;
      // }
    } else if ( (jentry*1000)/nentries == tenthpcount) {
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

    if (isCISample) {
      // have to choose which cut to use
      // if (!passMInvCut)
      // 	 continue;
      if (!passPreFSRMInvCut) {
	if (debug)
	  std::cout << "failed CI gen cut" << std::endl;
	continue;
      }
    }
    /*std::cout <<"=======> jentry = "<<jentry<<
      "=======> Evt = "<<event_evtNo<<
      "=======> Run = "<<event_runNo<<
      "=======> Lumi = "<<event_lumi<<
      "=======> bunch = "<<event_bunch<< std::endl; */
    //============================================================
    //                                                           =
    // Calling methods to get events with 2 electrons passing ID =
    //                                                           =
    //============================================================
    bool firstEleFinal  = SelectFirstEle(Etele1,Enele1,EtaTrakele1,PhiTrakele1,Chargeele1,
                                         EtaSCele1,PhiSCele1,flagel1,
                                         m_genET1,m_genPhi1,m_genEta1,m_genEn1);
    bool secondEleFinal = SelectSecondEle(Chargeele1,flagel1,Etele1,EtaTrakele1,PhiTrakele1,
                                          Etele2,Enele2,EtaTrakele2,PhiTrakele2,Chargeele2,
                                          EtaSCele2,PhiSCele2,
                                          m_genET2,m_genPhi2,m_genEta2,m_genEn2);
    //=========================================================
    //        call the method for N-1 plots                   =
    //                                                        =
    //=========================================================
    if (firstEleFinal == 0 || secondEleFinal == 0) {
      if (debug)
	std::cout << "failed el selection" << std::endl;
      continue;
    }

    m_recoMass = Mass(Etele1,EtaTrakele1,PhiTrakele1,ELEC_MASS,
		      Etele2,EtaTrakele2,PhiTrakele2,ELEC_MASS);

    if (m_recoMass < 60) {
      if (debug)
	std::cout << "failed m_recoMass" << std::endl;
      continue;
    }

    m_genMass = GenMass(m_genET1, m_genEta1, m_genPhi1, m_genEn1,
			m_genET2, m_genEta2, m_genPhi2, m_genEn2);

    m_recoMassSmeared = smearedMass(EtaTrakele1, EtaTrakele2, m_recoMass, m_genMass, m_scaleUnc);

    //=========================================================
    //        start doing matching between reco & HLT         =
    //                                                        =
    //=========================================================
    bool fireHLT2 = isPassHLT();
    if (fireHLT2 == 0) {
      if (debug)
	std::cout << "failed HLT" << std::endl;
      continue;
    }

    bool RecoEle1MatchingWithHLT1 = RecoHLTEleMatching(EtaSCele1,PhiSCele1);
    bool RecoEle2MatchingWithHLT2 = RecoHLTEleMatching(EtaSCele2,PhiSCele2);

    if (RecoEle1MatchingWithHLT1==1 && RecoEle2MatchingWithHLT2==1) {
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
		      << m_recoMass << ")" << std::endl;
	  newweight = 0;
	  wwto2l2nu_fail_gen_mass++;
	}
	if (m_recoMass > 600.) {  // WHY CUT ON THE RECO MASS???
	  if (debug)
	    std::cout << "Reweighting sample of WWTo2L2Nu with weight=0: gen("
		      << m_genMass << ") reco("
		      << m_recoMass << ")" << std::endl;
	  // newweight = 0;
	  wwto2l2nu_fail_reco_mass++;
	}
      } else if (inputfile.Contains("WWTo2L2Nu_Mll")) {
	wwto2l2nu_input++;
	if (m_genMass < 600.) {
	  if (debug)
	    std::cout << "Reweighting sample of WWTo2L2Nu_Mll with weight=0: gen("
		      << m_genMass << ") reco("
		      << m_recoMass << ")" << std::endl;
	  newweight = 0;
	  wwto2l2nu_fail_gen_mass++;
	}
	if (m_recoMass < 600.) {  // WHY CUT ON THE RECO MASS???
	  if (debug)
	    std::cout << "Reweighting sample of WWTo2L2Nu_Mll with weight=0: gen("
		      << m_genMass << ") reco("
		      << m_recoMass << ")" << std::endl;
	  // newweight = 0;
	  wwto2l2nu_fail_reco_mass++;
	}
      } else if (inputfile.Contains("TTTo2L2Nu_Tune")) {
	// Inclusive TTTo2L2Nu sample
	if (debug)
	  std::cout << "Checking reweighting of inclusive TTTo2L2Nu sample: gen("
		    << m_genMass << ") reco("
		    << m_recoMass << ")" << std::endl;
	ttto2l2nu_input++;
	if (m_genMass > 500.) {
	  if (debug)
	    std::cout << "Reweighting inclusive TTTo2L2Nu sample with weight=0: gen("
		      << m_genMass << ") reco("
		      << m_recoMass << ")" << std::endl;
	  newweight = 0;
	  ttto2l2nu_fail_gen_mass++;
	}
	if (m_recoMass > 500.) {  // WHY CUT ON THE RECO MASS???
	  if (debug)
	    std::cout << "Reweighting inclusive TTTo2L2Nu sample with weight=0: gen("
		      << m_genMass << ") reco("
		      << m_recoMass << ")" << std::endl;
	  // newweight = 0;
	  ttto2l2nu_fail_reco_mass++;
	}
      } else if (inputfile.Contains("TTTo2L2Nu_M") || inputfile.Contains("TTToLL_MLL_")) {
	if (debug)
	  std::cout << "Checking reweighting of mass binned TTTo2L2Nu sample: gen("
		    << m_genMass << ") reco("
		    << m_recoMass << ")" << std::endl;
	ttto2l2nu_input++;
	if (m_genMass > 500.) {
	  if (debug)
	    std::cout << "Reweighting mass binned TTTo2L2Nu sample with weight=0: gen("
		      << m_genMass << ") reco("
		      << m_recoMass << ")" << std::endl;
	  newweight = 0;
	  ttto2l2nu_fail_gen_mass++;
	}
	if (m_recoMass<500.) {
	  if (debug)
	    std::cout << "Reweighting mass binned TTTo2L2Nu sample with weight=0: gen("
		      << m_genMass << ") reco("
		      << m_recoMass << ")" << std::endl;
	  // newweight = 0;
	  ttto2l2nu_fail_reco_mass++;
	}
      }

      m_csAngle = CosThetaCollinSoper(Etele1,EtaSCele1,PhiSCele1,Enele1,
				      Etele2,EtaSCele2,PhiSCele2,Enele2,
				      Chargeele1,m_recoMass);
      PlotRecoInfo(m_recoMass,EtaSCele1,EtaSCele2);
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
    std::cout << "===Low mass WWTo2L2Nu sample info===="  << std::endl
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
  printf("It took you %.2lf minutes to run your program.\n", (dif/60.0));
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
float ZprimeEleElePatMiniAodNewMC::delR(float eta1,float phi1,float eta2,float phi2)
{
  float mpi=3.14;
  float dp=std::abs(phi1-phi2);
  if (dp>mpi) dp-=float(2*mpi);
  return sqrt((eta1-eta2)*(eta1-eta2) + dp*dp);
}

bool ZprimeEleElePatMiniAodNewMC::GenRecoMatchEle(float RecoEta1,float RecoPhi1,
                                                  float &ptGele,float &etaGele,float &phiGele,float &enGele)
{
  int NbHEEPele = 0;
  unsigned iflag = -10;
  ptGele  = -10.;
  etaGele = -10.;
  phiGele = -10.;
  enGele  = -10.;
  for (unsigned i=0; i<iGen->size(); i++) {
    if (fabs(idGen->at(i)) != 11) continue;
    if (statusGen->at(i) != 1)    continue;
    float deltaR1   = delR(RecoEta1,RecoPhi1,etaGen->at(i),phiGen->at(i));
    if (fabs(deltaR1)>deltaRcut)  continue;
    iflag  = i;
    NbHEEPele++;
    ptGele  = ptGen->at(i);
    etaGele = etaGen->at(i);
    phiGele = phiGen->at(i);
    enGele  = EnergyGen->at(i);
  }
  if (NbHEEPele > 0) {
    return true;
  } else {
    return false;
  }
}

//============================ Method to select first high pt ele ========================
bool ZprimeEleElePatMiniAodNewMC::SelectFirstEle(float &ETele1,float &Enele1,float &Etaele1,
                                                 float &Phiele1,int &ChargeEle1,float &EtaSCele1,
                                                 float &PhiSCele1,unsigned &FlagEle1,
                                                 float &genEle1Pt, float &genEle1Eta, float &genEle1Phi, float &genEle1En)
{
  int NbHEEPele = 0;
  unsigned iflag = -10;
  float highestpt=-999.;
  float ET = 0.0;
  for (unsigned i=0; i<Ele_nbElectrons->size(); i++) {
    ET = Ele_EtFromCaloEn->at(i);
    //Barrel
    if (ET > 35 && fabs(Ele_etaSC->at(i)) < 1.4442 && Ele_isPassHeepID->at(i)==1) {
      if (ET>highestpt) {
	bool GenRecoMatch1 = GenRecoMatchEle(Ele_etaSC->at(i),Ele_phiSC->at(i),
                                             genEle1Pt, genEle1Eta, genEle1Phi, genEle1En);
	if (GenRecoMatch1 == 0) continue;
	highestpt=ET;
	iflag  = i;
	NbHEEPele++;
      }
    } else if (ET > 35 && fabs(Ele_etaSC->at(i)) > 1.566 && fabs(Ele_etaSC->at(i)) < 2.5 && Ele_isPassHeepID->at(i)==1) { //endcap
      if (ET>highestpt) {
	bool GenRecoMatch2 = GenRecoMatchEle(Ele_etaSC->at(i),Ele_phiSC->at(i),
                                             genEle1Pt, genEle1Eta, genEle1Phi, genEle1En);
	if (GenRecoMatch2 == 0) continue;
	highestpt=ET;
	iflag  = i;
	NbHEEPele++;
      }
    } else {
      continue;
    }
  }

  if (NbHEEPele > 0) {
    FlagEle1           = iflag;
    ETele1             = Ele_EtFromCaloEn->at(iflag);
    Enele1             = Ele_energySC->at(iflag);
    Etaele1            = Ele_etaTrack->at(iflag);
    Phiele1            = Ele_phiTrack->at(iflag);
    ChargeEle1         = Ele_charge->at(iflag);
    EtaSCele1          = Ele_etaSC->at(iflag);
    PhiSCele1          = Ele_phiSC->at(iflag);
    return true;
  } else {
    return false;
  }
}
//============================ Method to select second high pt ele ========================
bool ZprimeEleElePatMiniAodNewMC::SelectSecondEle(int ChargeEle1,unsigned FlagEle1,float ETele1,float Etaele1,float Phiele1,
						  float &ETele2,float &Enele2,float &Etaele2,float &Phiele2,int &ChargeEle2,
						  float &EtaSCele2,float &PhiSCele2,
                                                  float &genEle1Pt, float &genEle1Eta, float &genEle1Phi, float &genEle1En)
{
  int NbHEEPele = 0;
  unsigned iflag = -10;
  float highestpt=-999.;
  float ET = 0.0;
  for (unsigned i=0; i<Ele_nbElectrons->size(); i++) {
    if (i == FlagEle1) continue;
    if (Ele_EtFromCaloEn->at(i) == ETele1) continue;
    if (Ele_etaSC->at(i) == Etaele1) continue;
    if (Ele_phiSC->at(i) == Phiele1) continue;
    ET = Ele_EtFromCaloEn->at(i);

    if (ET > 35 && fabs(Ele_etaSC->at(i)) < 1.4442 && Ele_isPassHeepID->at(i)==1) {  //Barrel
      if (ET>highestpt) {
	bool GenRecoMatch1 = GenRecoMatchEle(Ele_etaSC->at(i),Ele_phiSC->at(i),
                                             genEle1Pt, genEle1Eta, genEle1Phi, genEle1En);
	if (GenRecoMatch1 == 0) continue;
	highestpt=ET;
	iflag  = i;
	NbHEEPele++;
      }
    } else if (ET > 35 && fabs(Ele_etaSC->at(i)) > 1.566 && fabs(Ele_etaSC->at(i)) < 2.5 && Ele_isPassHeepID->at(i)==1) {  //endcap
      if (ET>highestpt) {
	bool GenRecoMatch2 = GenRecoMatchEle(Ele_etaSC->at(i),Ele_phiSC->at(i),
                                             genEle1Pt, genEle1Eta, genEle1Phi, genEle1En);
	if (GenRecoMatch2 == 0) continue;
	highestpt=ET;
	iflag  = i;
	NbHEEPele++;
      }
    } else {
      continue;
    }
  }

  if (NbHEEPele > 0) {
    ETele2             = Ele_EtFromCaloEn->at(iflag);
    Enele2             = Ele_energySC->at(iflag);
    Etaele2            = Ele_etaTrack->at(iflag);
    Phiele2            = Ele_phiTrack->at(iflag);
    ChargeEle2         = Ele_charge->at(iflag);
    EtaSCele2          = Ele_etaSC->at(iflag);
    PhiSCele2          = Ele_phiSC->at(iflag);
    return true;
  } else {
    return false;
  }
}
void ZprimeEleElePatMiniAodNewMC::PlotRecoInfo(float MassEle,float etaEle1,float etaEle2)
{
  //----------------------------------------------------------
  if (fabs(etaEle1) < 1.4442 && fabs(etaEle2) < 1.4442) {
    if (MassEle>900.0) {
      output_txt << event_runNo
                 << "   " << event_lumi
                 << "       " << event_evtNo
                 << "        " << MassEle
                 << "        " << etaEle1
                 << "        " << etaEle2
                 << std::endl;

    }
  }

  if (fabs(etaEle1) < 1.4442 && (fabs(etaEle2) > 1.566 && fabs(etaEle2) < 2.5)) {
    if (MassEle>900.0) {
      output_txt << event_runNo
                 << "   " << event_lumi
                 << "       " << event_evtNo
                 << "        " << MassEle
                 << "        " << etaEle1
                 << "        " << etaEle2
                 << std::endl;

    }
  }
  if (fabs(etaEle2) < 1.4442 && (fabs(etaEle1) > 1.566 && fabs(etaEle1) < 2.5)) {
    if (MassEle>900.0) {
      output_txt << event_runNo
                 << "   " << event_lumi
                 << "       " << event_evtNo
                 << "        " << MassEle
                 << "        " << etaEle1
                 << "        " << etaEle2
                 << std::endl;

    }
  }

  // only for DY POWHEG??
  float weight10 = 1.;
  if (inputfile.Contains("NNPDF30"))
    weight10 = MassCorrection(MassEle);

  newweight = newweight*weight10;
  m_recoMassCorr = m_recoMass*weight10;

  //----------------------------------------------------------
  if (!(inputfile.Contains("WW") && MassEle>2000.)) {
    // h2_CSSmearedMassBinned_->Fill(m_smearedMass,        0.,newweight);
    h2_CSMassBinned_       ->Fill(MassEle,               0.,newweight);
    h2_CSMassUpBinned_     ->Fill(MassEle*(1+m_scaleUnc),0.,newweight);
    h2_CSMassDownBinned_   ->Fill(MassEle*(1-m_scaleUnc),0.,newweight);
    if (m_csAngle > 0) {
      // h2_CSSmearedMassBinned_->Fill(m_smearedMass,        2.,newweight);
      h2_CSMassBinned_       ->Fill(MassEle,               2.,newweight);
      h2_CSMassUpBinned_     ->Fill(MassEle*(1+m_scaleUnc),2.,newweight);
      h2_CSMassDownBinned_   ->Fill(MassEle*(1-m_scaleUnc),2.,newweight);
    } else {
      // h2_CSSmearedMassBinned_->Fill(m_smearedMass,        1.,newweight);
      h2_CSMassBinned_       ->Fill(MassEle,               1.,newweight);
      h2_CSMassUpBinned_     ->Fill(MassEle*(1+m_scaleUnc),1.,newweight);
      h2_CSMassDownBinned_   ->Fill(MassEle*(1-m_scaleUnc),1.,newweight);
    }

    int priEtaBin = -1;
    int secEtaBin = -1;

    if (fabs(etaEle1) < 1.4442 && fabs(etaEle2) < 1.4442) {  //BB
      h1_ZprimeRecomassBinWidthBB_->Fill(MassEle,newweight);
      h1_ZprimeRecomass60to120BB_->Fill(MassEle,newweight);
      h1_ZprimeRecomassBinWidth_->Fill(MassEle,newweight);
      h1_ZprimeRecomass60to120_->Fill(MassEle,newweight);
      h1_ZprimeRecomass_->Fill(MassEle,newweight);
      priEtaBin = 1;
    } else if ((fabs(etaEle1) < 1.4442 && (fabs(etaEle2) > 1.566 && fabs(etaEle2) < 2.5)) ||
               (fabs(etaEle2) < 1.4442 && (fabs(etaEle1) > 1.566 && fabs(etaEle1) < 2.5))) {  //BE
      h1_ZprimeRecomassBinWidthBE_->Fill(MassEle,newweight);
      h1_ZprimeRecomass60to120BE_->Fill(MassEle,newweight);
      h1_ZprimeRecomassBinWidth_->Fill(MassEle,newweight);
      h1_ZprimeRecomass60to120_->Fill(MassEle,newweight);
      h1_ZprimeRecomass_->Fill(MassEle,newweight);
      priEtaBin = 2;
    } else if ((fabs(etaEle1) > 1.566 && fabs(etaEle1) < 2.5) &&
               (fabs(etaEle2) > 1.566 && fabs(etaEle2) < 2.5)) {  //EE
      h1_ZprimeRecomassBinWidthEE_->Fill(MassEle,newweight);
      h1_ZprimeRecomass60to120EE_->Fill(MassEle,newweight);
      h1_ZprimeRecomass_->Fill(MassEle,newweight);
      priEtaBin = 3;
    }
    // h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (priEtaBin*3)+0,newweight);
    h2_CSMassBinned_       ->Fill(MassEle,               (priEtaBin*3)+0,newweight);
    h2_CSMassUpBinned_     ->Fill(MassEle*(1+m_scaleUnc),(priEtaBin*3)+0,newweight);
    h2_CSMassDownBinned_   ->Fill(MassEle*(1-m_scaleUnc),(priEtaBin*3)+0,newweight);
    if (secEtaBin > 0) {
      // std::cout << "secondary bin (" << secEtaBin << "*3)+0(" << (secEtaBin*3)+0 << ")" << std::endl;
      // h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (secEtaBin*3)+0,newweight);
      h2_CSMassBinned_       ->Fill(MassEle,               (secEtaBin*3)+0,newweight);
      h2_CSMassUpBinned_     ->Fill(MassEle*(1+m_scaleUnc),(secEtaBin*3)+0,newweight);
      h2_CSMassDownBinned_   ->Fill(MassEle*(1-m_scaleUnc),(secEtaBin*3)+0,newweight);
    }
    if (m_csAngle > 0) {
      // std::cout << "primary bin (" << priEtaBin << "*3)+2(" << (priEtaBin*3)+2 << ")" << std::endl;
      // h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (priEtaBin*3)+2,newweight);
      h2_CSMassBinned_       ->Fill(MassEle,               (priEtaBin*3)+2,newweight);
      h2_CSMassUpBinned_     ->Fill(MassEle*(1+m_scaleUnc),(priEtaBin*3)+2,newweight);
      h2_CSMassDownBinned_   ->Fill(MassEle*(1-m_scaleUnc),(priEtaBin*3)+2,newweight);
      if (secEtaBin > 0) {
	// std::cout << "secondary bin (" << secEtaBin << "*3)+2(" << (secEtaBin*3)+2 << ")" << std::endl;
	// h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (secEtaBin*3)+2,newweight);
	h2_CSMassBinned_       ->Fill(MassEle,               (secEtaBin*3)+2,newweight);
	h2_CSMassUpBinned_     ->Fill(MassEle*(1+m_scaleUnc),(secEtaBin*3)+2,newweight);
	h2_CSMassDownBinned_   ->Fill(MassEle*(1-m_scaleUnc),(secEtaBin*3)+2,newweight);
      }
    } else {
      // std::cout << "primary bin (" << priEtaBin << "*3)+1(" << (priEtaBin*3)+1 << ")" << std::endl;
      // h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (priEtaBin*3)+1,newweight);
      h2_CSMassBinned_       ->Fill(MassEle,               (priEtaBin*3)+1,newweight);
      h2_CSMassUpBinned_     ->Fill(MassEle*(1+m_scaleUnc),(priEtaBin*3)+1,newweight);
      h2_CSMassDownBinned_   ->Fill(MassEle*(1-m_scaleUnc),(priEtaBin*3)+1,newweight);
      if (secEtaBin > 0) {
	// std::cout << "secondary bin (" << secEtaBin << "*3)+1(" << (secEtaBin*3)+1 << ")" << std::endl;
	// h2_CSSmearedMassBinned_->Fill(m_vtxMassSmearedMu,        (secEtaBin*3)+1,newweight);
	h2_CSMassBinned_       ->Fill(MassEle,               (secEtaBin*3)+1,newweight);
	h2_CSMassUpBinned_     ->Fill(MassEle*(1+m_scaleUnc),(secEtaBin*3)+1,newweight);
	h2_CSMassDownBinned_   ->Fill(MassEle*(1-m_scaleUnc),(secEtaBin*3)+1,newweight);
      }
    }
  }
}

//===================== Methode to calculate the mass ========================
float ZprimeEleElePatMiniAodNewMC::Mass(float Pt1,float Eta1,float Phi1,float En1,
                                        float Pt2,float Eta2,float Phi2,float En2)
{
  float EleEleMass = 0.0;
  TLorentzVector Ele1;
  TLorentzVector Ele2;
  /*Ele1.SetPtEtaPhiE(Pt1,Eta1,Phi1,En1);
    Ele2.SetPtEtaPhiE(Pt2,Eta2,Phi2,En2);*/
  Ele1.SetPtEtaPhiM(Pt1,Eta1,Phi1,En1);
  Ele2.SetPtEtaPhiM(Pt2,Eta2,Phi2,En2);
  EleEleMass = (Ele1+Ele2).M();
  return EleEleMass;
}



//===================== Method to calculate the smeared mass ========================
float ZprimeEleElePatMiniAodNewMC::smearedMass(float Eta1, float Eta2,
					       float recoMass, float genMass, float &scaleUnc)
{
  double a = 0.;
  double b = 0.;
  double c = 0.;
  double d = 0.;
  double e = 0.;

  //sigma BB
  if (fabs(Eta1) <= 1.4442 && fabs(Eta2) <= 1.4442) {
    scaleUnc = 0.01;
    a =  0.00701;
    b =  3.32e-05;
    c = -1.29e-08;
    d =  2.73e-12;
    e = -2.05e-16;
  } else if ((fabs(Eta1) < 1.4442 && (fabs(Eta2) > 1.566 && fabs(Eta2) < 2.5)) ||
	     (fabs(Eta2) < 1.4442 && (fabs(Eta1) > 1.566 && fabs(Eta1) < 2.5))) {  //BE
    scaleUnc = 0.02;
    a =  0.0124;
    b =  3.75e-05;
    c = -1.52e-08;
    d =  3.44e-12;
    e = -2.85e-16;
  } else if ((fabs(Eta1) > 1.566 && fabs(Eta1) < 2.5) &&
	     (fabs(Eta2) > 1.566 && fabs(Eta2) < 2.5)) {  //EE
    scaleUnc = 0.02; //?
    a =  0.0124;
    b =  3.75e-05;
    c = -1.52e-08;
    d =  3.44e-12;
    e = -2.85e-16;
  } else { // other?
    scaleUnc = 0.03; //?
    a =  0.0124;
    b =  3.75e-05;
    c = -1.52e-08;
    d =  3.44e-12;
    e = -2.85e-16;
  }
  double res = a + b*genMass + c*genMass*genMass + d*pow(genMass,3) + e*pow(genMass,4);

  double extraSmear = res*0.15;

  std::shared_ptr<TRandom3> rand = std::make_shared<TRandom3>(0); // why not a different random seed?
  return recoMass*rand->Gaus(1,extraSmear);
}

//----------------------------------------------------
//                                                   -
//       Part for Gen & Reco Matching                -
//                                                   -
//----------------------------------------------------
//========================== Method to select firt Gen Ele =======================
bool ZprimeEleElePatMiniAodNewMC::SelectFirstGenEle(float &ETEle1,float &PhiSCEle1,
                                                    float &EtaSCEle1,float &EnEle1,
                                                    int &IDele1,int &Statele1,
                                                    unsigned &GenFlag1)
{
  int NbHEEPele = 0;
  int iflag = -10;
  ETEle1 = 0.0;
  for (unsigned i=0; i<iGen->size(); i++) {
    if (fabs(idGen->at(i)) != 11) continue;
    if (statusGen->at(i) != 1)  continue;
    if (ptGen->at(i) > ETEle1) {
      ETEle1 = ptGen->at(i);
      iflag  = i;
      NbHEEPele++;
    }
    else continue;
  }
  if (NbHEEPele>0) {
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
bool ZprimeEleElePatMiniAodNewMC::SelectSecondGenEle(unsigned GenFlag1,float ETEle1,float &ETEle2,float &PhiSCEle2,
                                                     float &EtaSCEle2,float &EnEle2,int &IDele2,int &Statele2)
{
  int NbHEEPele = 0;
  int iflag = -10;
  ETEle2 = 0.0;
  for (unsigned i=0; i<iGen->size(); i++) {
    if (fabs(idGen->at(i)) != 11) continue;
    if (statusGen->at(i) != 1)  continue;
    if (i == GenFlag1) continue;
    if (fabs(ptGen->at(i) - ETEle1) <0.00001) continue;
    if (ptGen->at(i) > ETEle2) {
      ETEle2 = ptGen->at(i);
      iflag  = i;
      NbHEEPele++;
    }
    else continue;
  }
  if (NbHEEPele>0) {
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

//============================ Method to compute gen level invariant mass ========================
float ZprimeEleElePatMiniAodNewMC::GenMass(float ETEle1, float PhiEle1, float EtaEle1,float EnEle1,
					   float ETEle2, float PhiEle2, float EtaEle2,float EnEle2)
{
  TLorentzVector ele1, ele2;
  // ele1.SetPtEtaPhiE(ETEle1,EtaEle1,PhiEle1,EnEle1);
  // ele2.SetPtEtaPhiE(ETEle2,EtaEle2,PhiEle2,EnEle2);
  ele1.SetPtEtaPhiM(ETEle1,EtaEle1,PhiEle1,ELEC_MASS);
  ele2.SetPtEtaPhiM(ETEle2,EtaEle2,PhiEle2,ELEC_MASS);

  return (ele1+ele2).M();
}

float ZprimeEleElePatMiniAodNewMC::CosThetaCollinSoper(float Et1,float Eta1,float Phi1,float En1,
						       float Et2,float Eta2,float Phi2,float En2,
						       float ChargeEle1,float RecoMass)
{

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

  return costheta;
}
//----------------------------------------------------
//                                                   -
//       Part for HLT & Reco Matching                -
//                                                   -
//----------------------------------------------------
bool ZprimeEleElePatMiniAodNewMC::isPassHLT()
{
  int nbMatch = 0;
  for (unsigned i=0; i<HLT_nb->size(); i++) {
    if ((HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v2" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v4" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v5" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v6" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v8" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v9" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v10"||
	 HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v1" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v2" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v3" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v4" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v5" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v6" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v7" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v8" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v9" ||
         HLT_name->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v10") && HLT_isaccept->at(i) == 1) {
      //std::cout <<"triggerName = "<<triggerName<< std::endl;
      nbMatch++;
    }
  }
  if (nbMatch>0) {
    return true;
  } else {
    return false;
  }
}

bool ZprimeEleElePatMiniAodNewMC::RecoHLTEleMatching(float RecoEta,float RecoPhi)
{
  int nbMatch = 0;
  float deltaR   = -10000.0;
  for (unsigned i=0; i<HLTObj_nbObj->size(); i++) {
    //std::cout <<"[before]triggerName"<<HLTObj_collection->at(i) << std::endl;
    if (HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v2" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v4" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v5" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v6" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v7" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v8" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v9" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v10" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v1" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v2" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v3" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v4" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v5" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v6" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v7" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v8" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v9" ||
        HLTObj_collection->at(i) == "HLT_DoubleEle33_CaloIdL_MW_v10") {
      //std::cout <<"[after]triggerName"<<HLTObj_collection->at(i) << std::endl;
      deltaR   = delR(HLTObj_eta->at(i),HLTObj_phi->at(i),RecoEta,RecoPhi);
      //printf("HLT_Eta = %f  HLT_Phi = %f recoEta = %f recoPhi = %f DelR_trigger = %f\n",HLTObj_eta->at(i),HLTObj_phi->at(i),RecoEta,RecoPhi,deltaR);
      if (fabs(deltaR)>RecoHLTMatchingDeltaRcut)
	continue;
      nbMatch++;
    }
  }
  if (nbMatch>0) {
    return true;
  } else {
    return false;
  }
}


double ZprimeEleElePatMiniAodNewMC::MassCorrection(float M)
{
  float a =  1.06780e+00;
  float b = -1.20666e-04;
  float c =  3.22646e-08;
  float d = -3.94886e-12;
  double function = d*pow(M,3) + c*pow(M,2) + b*pow(M,1) + a;
  return function;
}
