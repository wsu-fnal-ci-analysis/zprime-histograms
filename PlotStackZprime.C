#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "ZZStyle.C"
#include "TFile.h"
#include "TColor.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"
#include "TMath.h"
#include "TSystem.h"
#include <libgen.h>

// using namespace std;

// Usage:
// .L PlotStackZprime.C+
// PlotStackZprime()

class PlotStackZprime{

public:
  PlotStackZprime(std::string const&, std::string const&, std::string const&, int, int, bool, bool, bool, bool, bool);
  void plotm4l(std::string);
  void setSamplesNames4l();
  void printnumbers(char*, TH1F*);
  /*void createdatacards(
		       float Higgsm, float channel, float energy,
		       float masslow, float masshigh,
		       float ggH, float qqH, float WH, float ZH, float ttH,
		       float bkg_qqzz, float bkg_ggzz, float bkg_zjets
		       );*/
  //void getMassWindow(float Higgsm);
  TH1F *DrawOverflow(TH1F *h);

private:
  std::vector<string> Vdatasetnamebkg,Vdatasetnamesig,Vdatasetnamedata,Vdatasetnamebkgdata;
  std::vector<string> Vlabelbkg,Vlabelsig,Vlabeldata,Vlabelbkgdata;
  std::vector<float> Vxsectionbkg,Vxsectionsig,Vxsectiondata,Vxsectionbkgdata;
  std::vector<Color_t> Vcolorbkg, Vcolorsig/*, Vcolordata*/;

  double Nbins;
  double Xmin;
  double Xmax;
  int nRebin,nRebinZ_X;
  int histErrType;
  double Ymax;
  bool m_setLogY,m_setLogX;
  bool m_useDYJets,m_useDYJetsFromData,m_useDiJetsFromFakeRateFromData,m_useWJetsFromFakeRateFromMC;
  bool m_useDiJetsWJetsFromFakeRateFromData;
  std::string histosdir;
  std::string inputfile;
  std::string ciSampleType;
  std::string whichchannel,whichchanneltex,whichenergy,whichsample,whichinvmasstex;
  char histotitle[500];
  ofstream outputyields;
  float errorZZ,errorH125,errorH126,errorH200,errorH350,errorH500,errorH800;
  TSystem LoadLib;
public:
  //float Higgsm, channel, energy, masslow, masshigh;
  //float ggH, qqH, WH, ZH, ttH, bkg_qqzz, bkg_ggzz, bkg_zjets;

};

PlotStackZprime::PlotStackZprime(std::string const& histolabel,
				 std::string const& sampletype="Lam16TeVConLL",
				 std::string const& leptonflavour="MuMu",
				 int histerrors=1, // 0 for kNormal, 1 for kPoisson, 2 for kPoisson2
				 int rebins=20,
				 bool dyjets=true,
				 bool dyjetsdata=false,
				 bool dijetsfrdata=false,
				 bool wjetsfrmc=false,
				 bool diwjetsfrdata=false)
{
  //TSystem LoadLib;
  //LoadLib.Load("/cmshome/nicola/slc6/MonoHiggs/Analysis13TeV/CMSSW_7_2_0/lib/slc6_amd64_gcc481/libHiggsHiggs_CS_and_Width.so");
  //getMassWindow(500.);

  if (leptonflavour.rfind("MuMu") != std::string::npos)
    inputfile="filelist_zprime_SingleMuon_2016_Spring16_25ns_AN.txt";
  else if (leptonflavour.rfind("EE") != std::string::npos)
    inputfile="filelist_zprime_DoubleEG_2016_Spring16_25ns_AN.txt";
  else
    return;

  ciSampleType = sampletype;
  histErrType  = histerrors;

  setSamplesNames4l();
  std::cout << "\t Analysing samples for " << whichchannel << " analysis" << std::endl;


  //WARNING: depending on histolabel, modify the declaration and the settings of hframe below
  //also choose a sensible value for nRebin

  //std::string histolabel = "hPUvertices";    // numPU
  //std::string histolabel = "hPUvertices_ReWeighted";    // numPY reweighted

  //std::string histolabel = "ZprimeRecomass";
  //std::string histolabel = "ZprimeRecomassBinWidth";

  //std::string histolabel = "dPToverPT";
  //std::string histolabel = "numberOfValidMuonHits";
  //std::string histolabel = "numberOfValidPixelHits";
  //std::string histolabel = "numberOfMatchedStations";
  //std::string histolabel = "numberOftrackerLayersWithMeasurement";
  //std::string histolabel = "trackiso";
  //std::string histolabel = "absdxy";
  //std::string histolabel = "PtID";

  //std::string histolabel = "CosAngleCollinSoperCorrect60Mass120";
  //std::string histolabel = "CosAngleCollinSoperCorrect120Mass300";
  //std::string histolabel = "CosAngleCollinSoperCorrect300Mass700";
  //std::string histolabel = "CosAngleCollinSoperCorrect700Mass3000";
  //std::string histolabel = "CosAngleCollinSoperCorrect4900Mass5100";

  m_setLogY = true;
  m_setLogX = false;

  m_useDYJets                          = dyjets;
  m_useDYJetsFromData                  = dyjetsdata;
  m_useDiJetsFromFakeRateFromData      = dijetsfrdata;
  m_useWJetsFromFakeRateFromMC         = wjetsfrmc;
  m_useDiJetsWJetsFromFakeRateFromData = diwjetsfrdata;

  nRebin = rebins;
  std::cout << "Histogram label is= " << histolabel << std::endl;

  // Final yields
  system("mkdir plots");

  Char_t yieldsOUT[500];
  sprintf(yieldsOUT,"plots/yields_%s_%s.txt",whichchannel.c_str(),whichenergy.c_str());
  if (histolabel.find("hM4l_8") < 10) {
    std::cout << "Opening a file for final numbers= " << yieldsOUT << std::endl;
    outputyields.open(yieldsOUT);
  }

  // Execute the analysis
  errorZZ=0.,errorH125=0.,errorH126=0.,errorH200=0.,errorH350=0.,errorH500=0.,errorH800=0.;
  plotm4l(histolabel);

  // close file for final yields
  if (histolabel.find("hM4l_8") < 10) outputyields.close();

}

void PlotStackZprime::plotm4l(std::string histlabel)
{
  TStyle * style = getStyle("ZZ");
  //style->SetMarkerSize(0.8);
  style->cd();
  style->SetNdivisions(508, "X");
  style->SetNdivisions(508, "Y");
  style->SetMarkerSize(0.8);
  //style->SetOptStat(111111);
//   double ytitleOffset = 1.36;
//   double xtitleOffset = 1.18;
//   double labelSize = 0.05;
//   double titleSize = 0.05;
//   double lineWidth = 2;

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->cd();
  c1->SetTicks(1,1);

  //  ratioPad->cd();

  // TString lumist="5.2 fb^{-1}";
//   TPaveText *ll = new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC");
//   ll->SetTextSize(0.03);
//   ll->SetTextFont(42);
//   ll->SetFillColor(0);
//   ll->SetBorderSize(0);
//   ll->SetMargin(0.01);
//   ll->SetTextAlign(12); // align left
//   TString text = "CMS Preliminary 2011";
//   ll->AddText(0.01,0.5,text);
//   text = "#sqrt{s} = 7 TeV  L = ";
//   text = text + lumist;
//   //  ll->SetTextAlign(32); // align right
//   ll->AddText(0.65, 0.6, text);


  TPaveText *ll = new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC");
  ll->SetTextSize(0.027);
  ll->SetTextFont(42);
  ll->SetFillColor(0);
  ll->SetBorderSize(0);
  ll->SetMargin(0.01);
  ll->SetTextAlign(12); // align left
  TString text = "CMS Preliminary";
  //TString text = "CMS";
  ll->AddText(0.01,0.5,text);
  std::cout << "Energy= " << whichenergy << std::endl;
  if (whichenergy.find("RunI") < 100) {
    text = "#sqrt{s} = 7 TeV, L = 5.05 fb^{-1} ; #sqrt{s} = 8 TeV, L = 19.71 fb^{-1}" ;
    ll->AddText(0.3, 0.6, text);
  }
  else if (whichenergy.find("7TeV") < 100) {
    text = "#sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" ;
    ll->AddText(0.65, 0.6, text);
  }
  else if (whichenergy.find("8TeV") < 100) {
    text = "#sqrt{s} = 8 TeV, L = 19.71 fb^{-1}" ;
    ll->AddText(0.65, 0.6, text);
  }
  else if (whichenergy.find("13TeV") < 100) {
    //text = "#sqrt{s} = 13 TeV, L = 50.852 pb^{-1}" ;
    //text = "#sqrt{s} = 13 TeV, L = 77.346 pb^{-1}" ;
    text = "#sqrt{s} = 13 TeV, L = 36.3 fb^{-1}" ;
    ll->AddText(0.65, 0.6, text);
  }
  //ll->Draw();



  TLegend *leg0 = new TLegend(0.6,0.40,0.8,0.90,NULL,"brNDC");
  leg0->SetTextSize(0.020);
  leg0->SetLineColor(0);
  leg0->SetLineWidth(1);
  leg0->SetFillColor(kWhite);
  leg0->SetBorderSize(0);

  TLegend* legend = new TLegend( 0.5, 0.74, 0.9, 0.92);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.020);


//   TLegend *leg1 = new TLegend(0.25,0.80,0.4,0.9,NULL,"brNDC");
//   leg1->SetTextSize(0.030);
//   leg1->SetLineColor(0);
//   leg1->SetLineWidth(1);
//   leg1->SetFillColor(kWhite);
//   leg1->SetBorderSize(0);

//   Color_t redBgColor = kGreen-5;
  Color_t ZZBgColor = kAzure-9;
//   Color_t h200Color = kOrange;
//   Color_t h350Color = kRed+1;
//   Color_t h400Color = kRed;




  if (m_setLogY)
    c1->SetLogy(1);
  else
    c1->SetLogy(0);

  if (m_setLogX)
    c1->SetLogx(1);
  else
    c1->SetLogx(0);

  std::cout << "Vdatasetnamedata " <<  Vdatasetnamedata.size() << std::endl;

  TFile *ff = NULL;

  if (Vdatasetnamedata.size() > 0)
    ff = TFile::Open(Vdatasetnamedata.at(0).c_str());  // just a random file, to determine the binning
  else if (Vdatasetnamebkg.size() > 0)
    ff = TFile::Open(Vdatasetnamebkg.at(0).c_str());
  else if (Vdatasetnamesig.size() > 0)
    ff = TFile::Open(Vdatasetnamesig.at(0).c_str());


  TH1F *hfourlepbestmass_4l_afterSel_new = (TH1F*)ff->Get(histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new"*/);
  //hfourlepbestmass_4l_afterSel_new->SetBins(9940,60.,10000.);

  Nbins = hfourlepbestmass_4l_afterSel_new->GetNbinsX() / nRebin;
  Xmin  = hfourlepbestmass_4l_afterSel_new->GetXaxis()->GetXmin();
  Xmax  = hfourlepbestmass_4l_afterSel_new->GetXaxis()->GetXmax() ;
  Ymax  = hfourlepbestmass_4l_afterSel_new->GetBinContent(hfourlepbestmass_4l_afterSel_new->GetMaximumBin()) * 580.;

  std::cout << "Ymax = " << Ymax << std::endl;

  TH2F *hframe = NULL, *hframe2 = NULL;

  hframe= new TH2F("hframe","hframe",80,70.,4000.,500,0.0005,150000.);// di-mu  mass nrebin=10  GeV
  hframe2= new TH2F("hframe2","hframe2",80,70.,4000.,500, 0.5, 1.5);// di-mu mass


  if (histlabel.find("ZprimeRecomass") < 20 && nRebin == 20) {
    hframe= new TH2F("hframe","hframe",100,70.,3000.,500,0.002,5000000.);// di-mu  mass nrebin=10  GeV
    hframe2= new TH2F("hframe2","hframe2",100,70.,3000.,1000, 0., 2.);// di-mu mass
  }

  if (histlabel.find("ZprimeRecomass") < 20 && nRebin == 5) {
    hframe= new TH2F("hframe","hframe",100,70.,3000.,500,0.0005,250000.);// di-mu  mass nrebin=5  GeV
    hframe2= new TH2F("hframe2","hframe2",100,70.,3000.,1000, 0.5, 1.5);// di-mu mass
  }

  if (histlabel.find("ZprimeRecomass") < 20 && nRebin == 10) {
    hframe= new TH2F("hframe","hframe",100,70.,3000.,1000,0.001,250000.);// mZ mumu
    hframe2= new TH2F("hframe2","hframe2",100, 70., 3000., 1000, 0.5, 1.5);// mZ mumu
  }

  if (histlabel.find("dPToverPT") < 10) {
    hframe= new TH2F("hframe","hframe",100,0.0,0.5,500,0.001,30000050.);// dPToverPT in control region
    hframe2= new TH2F("hframe2","hframe2",100,0.0,0.5, 1000, 0.5, 2.);// dPToverPT in control region
  }

  if (histlabel.find("numberOfValidMuonHits") < 10) {
    hframe= new TH2F("hframe","hframe",100,0.,60.,500,0.0001,100000000050.);// numberOfValidMuonHits in control region
    hframe2= new TH2F("hframe2","hframe2",100,0.,60., 1000, 0.5, 2.);// numberOfValidMuonHits in control region
  }

  if (histlabel.find("numberOfValidPixelHits") < 10) {
    hframe= new TH2F("hframe","hframe",100,0.,10.,500,0.001,30000050.);// numberOfValidPixelHits in control region
    hframe2= new TH2F("hframe2","hframe2",100,0., 10., 1000, 0.5, 2.);// numberOfValidPixelHits in control region
  }

  if (histlabel.find("numberOfMatchedStations") < 10) {
    hframe= new TH2F("hframe","hframe",100,0.,10.,500,0.0001,100000000050.);// numberOfMatchedStations in control region
    hframe2= new TH2F("hframe2","hframe2",100,0., 10., 1000, 0.5, 2.);// numberOfMatchedStations in control region
  }

  if (histlabel.find("numberOftrackerLayersWithMeasurement") < 10) {
    hframe= new TH2F("hframe","hframe",100,0.,20.,500,0.001,30000050.);// numberOftrackerLayersWithMeasurement in control region
    hframe2= new TH2F("hframe2","hframe2",100,0., 20., 1000, 0.5, 2.);// numberOftrackerLayersWithMeasurement in control region
  }

  if (histlabel.find("trackiso") < 10) {
    hframe= new TH2F("hframe","hframe",100,0.,0.3,500,0.001,100000000050.);// trackiso in control region
    hframe2= new TH2F("hframe2","hframe2",100,0., 0.3, 1000, 0.5, 2.);// trackiso in control region
  }

  if (histlabel.find("absdxy") < 10) {
    hframe= new TH2F("hframe","hframe",100,0.,0.3,500,0.001,100000000050.);// trackiso in control region
    hframe2= new TH2F("hframe2","hframe2",100,0., 0.3, 1000, 0.5, 2.);// trackiso in control region
  }

  if (histlabel.find("PtID") < 10) {
    hframe= new TH2F("hframe","hframe",100,50.,400.,500,0.001,100000000050.);// pT  in control region
    hframe2= new TH2F("hframe2","hframe2",100,50.,400., 1000, 0.5, 2.);// pT in control region
  }

  if (histlabel.find("CosAngleCollinSoperCorrect60Mass120") < 10) {
    hframe= new TH2F("hframe","hframe",100,-1.,1.,500,0.001,100000000050.);// costheeta CS
    hframe2= new TH2F("hframe2","hframe2",100,-1.,1., 1000, 0.5, 2.);// costheta CS
  }

  if (histlabel.find("CosAngleCollinSoperCorrect300Mass700") < 10) {
    hframe= new TH2F("hframe","hframe",100,-1.,1.,500,0.001,1000.);// costheeta CS
    hframe2= new TH2F("hframe2","hframe2",100,-1.,1., 1000, 0.5, 2.);// costheta CS
  }

  if (nRebin == 5) hframe->SetYTitle("Events/5 GeV");
  if (nRebin == 1) hframe->SetYTitle("Events/1 GeV");
  if (nRebin == 10) hframe->SetYTitle("Events/10 GeV");
  if (nRebin == 20) hframe->SetYTitle("Events/20 GeV");

  //hframe->SetXTitle("M_{Z2} [GeV]");
  if (histlabel.find("ZprimeRecomass") < 20) hframe->SetXTitle("m_{" + TString(whichinvmasstex) + "} [GeV]");
  if (histlabel.find("Den_Pt") < 10 || histlabel.find("Num_Pt") < 10) hframe->SetXTitle("p_{T} [GeV]");

  if (histlabel.find("CosAngleCollinSoperCorrect60Mass120") < 10) {
    if (nRebin == 10) hframe->SetYTitle("Events/ bin=10");
    hframe->SetXTitle("cos(#theta), 60<m_{" + TString(whichinvmasstex) + "}<120 GeV");
  }

  if (histlabel.find("CosAngleCollinSoperCorrect300Mass700") < 10) {
    if (nRebin == 10) hframe->SetYTitle("Events/ bin=10");
    hframe->SetXTitle("cos(#theta), 300<m_{" + TString(whichinvmasstex) + "}<700 GeV");
  }

  hframe->GetYaxis()->SetLabelSize(0.03);
  hframe->GetXaxis()->SetLabelSize(0.03);

  hframe->GetYaxis()->SetTitleSize(0.03);
  hframe->GetXaxis()->SetTitleSize(0.03);

  hframe->GetXaxis()->SetLabelOffset(0.007);
  hframe->GetXaxis()->SetTitleOffset(1.5);
  hframe->GetYaxis()->SetLabelOffset(0.007);
  hframe->GetYaxis()->SetTitleOffset(2.);

  hframe->Draw();

  TH1F *htotaldata = new TH1F("htotaldata", "htotaldata", Nbins, Xmin, Xmax);
  //TH1F *htotaldata = new TH1F("htotaldata", "htotaldata", 51, 95., 605.);
  //TH1F *htotaldata = new TH1F("htotaldata", "htotaldata", 50, 99.5, 599.5);
  htotaldata->SetMarkerColor(1);
  htotaldata->SetMarkerStyle(20);
  THStack *htotalbkg     = new THStack("bkgstack","");
  THStack *htotalsig     = new THStack("sigstack","");
  THStack *htotal        = new THStack("Nicola","");
  TH1F *htotalHisto      = new TH1F("htotalHisto", "htotalHisto", Nbins, Xmin, Xmax);
  TH1F *htotalHistoRatio = new TH1F("htotalHistoRatio", "htotalHistoRatio", Nbins, Xmin, Xmax);

  // data

  for (unsigned int datasetIdData=0; datasetIdData<Vdatasetnamedata.size(); datasetIdData++) {
    char dataset[328];
    sprintf(dataset,"%s",Vdatasetnamedata.at(datasetIdData).c_str());
    std::cout << "Root-ple= " << dataset << std::endl;
    // std::cout << "Counter=" << datasetIdData << " Root-ple=" << dataset << " Label=" << Vlabelbkg.at(datasetIdData) <<endl;

    TFile *f1 = TFile::Open(dataset);
    hfourlepbestmass_4l_afterSel_new = (TH1F*)f1->Get(histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new"*/);
    //hfourlepbestmass_4l_afterSel_new->SetBins(3000,0.,3000.);
    //hfourlepbestmass_4l_afterSel_new->SetBins(9940,60.,10000.);

    TH1 *hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebin, histlabel.c_str()/*"hfourlepbestmass_4l_afterSel_new_new"*/);
    hfourlepbestmass_4l_afterSel_new_new->SetMarkerColor(1);
    hfourlepbestmass_4l_afterSel_new_new->SetMarkerStyle(20);
    hfourlepbestmass_4l_afterSel_new_new->SetMarkerSize(0.95);
    // hfourlepbestmass_4l_afterSel_new_new->Draw("EPsame");

    //leg0->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabeldata.at(datasetIdData).c_str(), "P");

    //std:cout << "Nbins=" << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX() << std::endl;
    //std:cout << "htotaldata nBins = " << htotaldata->GetNbinsX() << ", hfourlepbestmass_4l_afterSel_new_new nBins = " << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX() << std::endl;
    //std:cout << "htotaldata lowestX = " << htotaldata->GetXaxis()->GetXmin() <<  ", htotaldata highestX = " << htotaldata->GetXaxis()->GetXmax() << ", hfourlepbestmass_4l_afterSel_new_new lowestX = " << hfourlepbestmass_4l_afterSel_new_new->GetXaxis()->GetXmin() << ", hfourlepbestmass_4l_afterSel_new_new highestX = " << hfourlepbestmass_4l_afterSel_new_new->GetXaxis()->GetXmax() << std::endl;

    htotaldata->Add(hfourlepbestmass_4l_afterSel_new_new);
    std::cout << "(l400)Label= " << Vlabeldata.at(datasetIdData)
	      << "  Entries= "   << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
	      << "  Integral= "  << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) << std::endl;
    if (datasetIdData == (Vdatasetnamedata.size()-1)) {
      leg0->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"Data 2015", "P");
      //legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"ExpressPhysics 2015B - Run 251027-251721", "P");
      //legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"ExpressPhysics 2015B - Run 251027-251883", "P");
      //legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"SingleMuon 2015B - Run 246908-251883", "P");
      //legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"SingleMuon 2015B+2015C+2015D", "P");
      legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,"Data", "P");
      //htotaldata->Draw("EPsame");
    }
  }

  // Set Errors as in http://www-cdf.fnal.gov/physics/statistics/notes/pois_eb.txt
  //Float_t x[51],y[51],exl[51],exh[51],eyl[51],eyh[51];

  int* arraysize = new int[1];
  arraysize[0] = htotaldata->GetNbinsX();
  //std::cout << "arraysize = " << arraysize[0] << std::endl;
  Float_t x[arraysize[0]],y[arraysize[0]],exl[arraysize[0]],exh[arraysize[0]],eyl[arraysize[0]],eyh[arraysize[0]];
  delete [] arraysize;

  float totaldataentries=0.,totaldataentries100=0.;

  if (histErrType == 0)
    htotaldata->SetBinErrorOption(TH1::kNormal);
  else if (histErrType == 1)
    htotaldata->SetBinErrorOption(TH1::kPoisson);
  else if (histErrType == 2)
    htotaldata->SetBinErrorOption(TH1::kPoisson2);

  for (int nbins=1;nbins<=htotaldata->GetNbinsX(); nbins++) {
    // std::cout << "BinCenter=" << htotaldata->GetBinCenter(nbins) << " BinContent=" << htotaldata->GetBinContent(nbins) << " BinErrorContent=" << htotaldata->GetBinError(nbins) << std::endl;
    x[nbins-1]=htotaldata->GetBinCenter(nbins);
    y[nbins-1]=htotaldata->GetBinContent(nbins);
    exl[nbins-1]=0.;
    exh[nbins-1]=0.;
    totaldataentries=totaldataentries+htotaldata->GetBinContent(nbins);
    if (htotaldata->GetBinCenter(nbins) > 100. && htotaldata->GetBinCenter(nbins) < 800.)
      totaldataentries100=totaldataentries100+htotaldata->GetBinContent(nbins);
    if (htotaldata->GetBinContent(nbins) > 0) {
        eyh[nbins-1]=+0.5 + sqrt(htotaldata->GetBinContent(nbins)+0.25);
        eyl[nbins-1]=-0.5 + sqrt(htotaldata->GetBinContent(nbins)+0.25);
    } else {
      x[nbins-1]   = 0.;
      eyl[nbins-1] = 0.;
      eyh[nbins-1] = 0.;
    }
  }

  std::cout << "Total data= " << totaldataentries << std::endl;
  if (histlabel.find("hM4l_8") < 10) outputyields << "Data " << totaldataentries << " +/- 0" << std::endl;
  std::cout << "Total data in the range m4l>100 and  m4l<800= " << totaldataentries100 << std::endl;

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(Nbins,x,y,exl,exh,eyl,eyh);
  gr->SetMarkerColor(1);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.95);
  //

  //Dijets by FR from data
  TH1 *hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new_new = NULL;
  if (m_useDiJetsFromFakeRateFromData == true) {
    std::cout << "Estimating the di-jet contribution from data with the fake rate method" << std::endl;
    if (histlabel.find("ZprimeRecomass") < 20) {
      TFile *fDiJets = TFile::Open("Data-total-jets.root");
      TH1F *hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData = (TH1F*)fDiJets->Get("TotalJets");
      hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData->SetBins(10000,0.,10000.);
      //hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData->SetBins(9940,60.,10000.);

      hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new_new = hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData->Rebin(nRebin,"hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData");

      hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new_new->SetLineColor(kYellow);
      hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new_new->SetFillColor(kYellow);
      hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new_new->SetMarkerStyle(24);
      legend->AddEntry(hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new_new,"di-jets (FR from data)", "F");

    }
  }

  //Wjets by FR from MC
  TH1 *hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC_new_new = NULL;
  if (m_useWJetsFromFakeRateFromMC == true) {
    std::cout << "Estimating the W+jets contribution from MC with the fake rate method" << std::endl;
    if (histlabel.find("ZprimeRecomass") < 20) {
      TFile *fWJets = TFile::Open("plots/FR-estimate-Dijets-Wjets-MC-OS-1GeVBin-2800pb.root");
      TH1F *hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC = (TH1F*)fWJets->Get("WjetsHisto");
      hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC->SetBins(10000,0.,10000.);
      //hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC->SetBins(9940,60.,10000.);

      hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC_new_new = hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC->Rebin(nRebin,"hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC");
      hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC_new_new->SetLineColor(kYellow+2);
      hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC_new_new->SetFillColor(kYellow+2);
      hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC_new_new->SetMarkerStyle(24);
      legend->AddEntry(hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC_new_new,"W-Jets (FR from MC)", "F");

    }
  }

  //Dijets+WJets by FR from data
  TH1 *hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new = NULL;
  if (m_useDiJetsWJetsFromFakeRateFromData == true) {
    std::cout << "Estimating the di-jet and W+jets contribution from data with the fake rate method" << std::endl;
    if (histlabel.find("ZprimeRecomass") < 20) {
      char dataset[328];
      sprintf(dataset,"%s",Vdatasetnamebkgdata.at(0).c_str()); // just one file
      //TFile *fDiJetsWJets = TFile::Open("Data-total-jets.root");
      TFile *fDiJetsWJets = TFile::Open(dataset);
      TH1F *hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData = (TH1F*)fDiJetsWJets->Get("TotalJets"); //binned with 1 GeV bins

      // const int NMBINS = 5940;
      // const double MMIN = 60., MMAX = 6000.;

      // double linMbins[NMBINS+1];

      // for (int ibin = 0; ibin <= NMBINS; ibin++) {
      // 	linMbins[ibin] = MMIN + ibin;
      // 	cout << "Bin= " << ibin+1 << " Value= " << linMbins[ibin] << endl;
      // }

      //hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new = hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData->Rebin(99,"hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new",linMbins);
      //hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData->SetBins(10000,0.,10000.);
      // //hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData->SetBins(9940,60.,10000.);

      hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new = hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData->Rebin(nRebin,"TotalJets");

      hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new->Scale(0.9714); // to be understood

      hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new->SetLineColor(kYellow);
      hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new->SetFillColor(kYellow);
      hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new->SetMarkerStyle(24);

      htotal->Add(hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new);
      htotalHisto->Add(hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new);
      legend->AddEntry(hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new,"Jets", "F");

    }
  }


  // Z+X from data
  if (m_useDYJets == false && m_useDYJetsFromData == true) {
    std::cout << "Estimating the DY+jets contribution from data" << std::endl;
    if (histlabel.find("hM4l_7") < 10 || histlabel.find("hM4l_8") < 10) {
      for (unsigned int datasetIdbkgData=0; datasetIdbkgData<Vdatasetnamebkgdata.size(); datasetIdbkgData++) {
	char dataset[328];
	sprintf(dataset,"%s",Vdatasetnamebkgdata.at(datasetIdbkgData).c_str());
	cout << "Root-ple= " << dataset << std::endl;
	TFile *f3 = TFile::Open(dataset);
	TH1F *hfourlepbestmass_4l_afterSel_orig = new TH1F("hfourlepbestmass_4l_afterSel_orig", "Mass of four leptons after fullselection", 2460, 0.,1230.);
        hfourlepbestmass_4l_afterSel_orig = (TH1F*)f3->Get("h_3P1F_2P2F");
	//hfourlepbestmass_4l_afterSel_orig->SetBins(3000,0.,3000.);
	//hfourlepbestmass_4l_afterSel_orig->SetBins(9940,60.,10000.);

        hfourlepbestmass_4l_afterSel_new = new TH1F("hfourlepbestmass_4l_afterSel_new", "Mass of four leptons after fullselection", 2400, 4.5,1204.5);

        int mbins=1;
        for (int nbins=1;nbins<=hfourlepbestmass_4l_afterSel_orig->GetNbinsX(); nbins++) {
          if (hfourlepbestmass_4l_afterSel_orig->GetBinCenter(nbins) > 4.5 && hfourlepbestmass_4l_afterSel_orig->GetBinCenter(nbins) < 1204.5) {
            hfourlepbestmass_4l_afterSel_new->SetBinContent(mbins,double(hfourlepbestmass_4l_afterSel_orig->GetBinContent(nbins)));
            mbins++;
          }
        }

	TH1 *hfourlepbestmass_4l_afterSel_new_new;

	nRebinZ_X=nRebin*2;
	hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebinZ_X, "h_3P1F_2P2F");
	hfourlepbestmass_4l_afterSel_new_new->SetLineColor(kCyan-2);
	hfourlepbestmass_4l_afterSel_new_new->SetFillColor(kCyan-2);
	hfourlepbestmass_4l_afterSel_new_new->SetMarkerStyle(24);
	hfourlepbestmass_4l_afterSel_new_new->SetLineWidth(1);

	char temp[328];
        sprintf(temp,"%s",histosdir.c_str());

	if (histlabel.find("hM4l_7") < 10 && Vdatasetnamebkgdata.at(datasetIdbkgData).find("m4l_gt_70") < 85) {
	  std::cout << "Adding Z+X for m4l > 70. GeV" << std::endl;
	  std::cout << "N bins Z+X= " << hfourlepbestmass_4l_afterSel_new_new->GetNbinsX() << std::endl;
	  std::cout << "Z+X entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) << std::endl;
	  htotal->Add(hfourlepbestmass_4l_afterSel_new_new);
	  htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_new);
	  if (Vdatasetnamebkgdata.at(datasetIdbkgData).find(temp) <200 &&
	      (Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichenergy) < 200 || Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichsample) < 200) &&
	      Vdatasetnamebkgdata.at(datasetIdbkgData).find("2mu2e") > 85) {
	    std::cout << "Adding legend for Z+X" << std::endl;
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkgdata.at(datasetIdbkgData).c_str(), "F"); }
	} else if (Vdatasetnamebkgdata.at(datasetIdbkgData).find(temp) <200 &&
		   (Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichenergy) < 200 || Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichsample) < 200) &&
		   histlabel.find("hM4l_8") < 10 && !(Vdatasetnamebkgdata.at(datasetIdbkgData).find("m4l_gt_70") < 85)) {
	  std::cout << "Adding Z+X for m4l > 100. GeV" << std::endl;
	  std::cout << "Z+X entries= " << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) << std::endl;
	  htotal->Add(hfourlepbestmass_4l_afterSel_new_new);
	  htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_new);
	  if (Vdatasetnamebkgdata.at(datasetIdbkgData).find(temp) <200 &&
	      (Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichenergy) < 200 || Vdatasetnamebkgdata.at(datasetIdbkgData).find(whichsample) < 200) &&
	      Vdatasetnamebkgdata.at(datasetIdbkgData).find("2mu2e") > 85) {
	    std::cout << "Adding legend for Z+X" << std::endl;
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkgdata.at(datasetIdbkgData).c_str(), "F");}
	}
      }
    }
  }
  //std:cout << "Total Z+X is= " << htotal->GetHistogram()->GetEntries() << std::endl;
  std::cout << "Total Z+X is= " << htotalHisto->Integral(0,-1)  << std::endl;
  //outputyields << "Z+X "   << htotal->GetHistogram()->GetEntries() << " +/- 0"<< std::endl;
  if (histlabel.find("hM4l_8") < 10)
    outputyields << "Z+X " << htotalHisto->Integral(0,-1) << " +/- 0" << std::endl;

  // Background

  TH1F *hfourlepbestmass_4l_afterSel_new_qqZZ    = new TH1F("hfourlepbestmass_4l_afterSel_new_qqZZ", "hfourlepbestmass_4l_afterSel_new_qqZZ", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_ggZZ    = new TH1F("hfourlepbestmass_4l_afterSel_new_ggZZ", "hfourlepbestmass_4l_afterSel_new_ggZZ", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_ZZ      = new TH1F("hfourlepbestmass_4l_afterSel_new_ZZ",   "hfourlepbestmass_4l_afterSel_new_ZZ",   Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_qcdDEM= new TH1F("hfourlepbestmass_4l_afterSel_new_qcdDEM", "hfourlepbestmass_4l_afterSel_new_qcdDEM", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_qcdMu = new TH1F("hfourlepbestmass_4l_afterSel_new_qcdMu",  "hfourlepbestmass_4l_afterSel_new_qcdMu",  Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_qcdBC = new TH1F("hfourlepbestmass_4l_afterSel_new_qcdBC",  "hfourlepbestmass_4l_afterSel_new_qcdBC",  Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_qcd   = new TH1F("hfourlepbestmass_4l_afterSel_new_qcd",    "hfourlepbestmass_4l_afterSel_new_qcd",    Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_singlet = new TH1F("hfourlepbestmass_4l_afterSel_new_singlet", "hfourlepbestmass_4l_afterSel_new_singlet", Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_DY      = new TH1F("hfourlepbestmass_4l_afterSel_new_DY",      "hfourlepbestmass_4l_afterSel_new_DY",      Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_DYlight = new TH1F("hfourlepbestmass_4l_afterSel_new_DYlight", "hfourlepbestmass_4l_afterSel_new_DYlight", Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_DYbb     = new TH1F("hfourlepbestmass_4l_afterSel_new_DYbb",     "hfourlepbestmass_4l_afterSel_new_DYbb",     Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_DYcc     = new TH1F("hfourlepbestmass_4l_afterSel_new_DYcc",     "hfourlepbestmass_4l_afterSel_new_DYcc",     Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_DYtautau = new TH1F("hfourlepbestmass_4l_afterSel_new_DYtautau", "hfourlepbestmass_4l_afterSel_new_DYtautau", Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_WW      = new TH1F("hfourlepbestmass_4l_afterSel_new_WW",      "hfourlepbestmass_4l_afterSel_new_WW",      Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_WZ      = new TH1F("hfourlepbestmass_4l_afterSel_new_WZ",      "hfourlepbestmass_4l_afterSel_new_WZ",      Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_diBoson = new TH1F("hfourlepbestmass_4l_afterSel_new_diBoson", "hfourlepbestmass_4l_afterSel_new_diBoson", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_TT      = new TH1F("hfourlepbestmass_4l_afterSel_new_TT",      "hfourlepbestmass_4l_afterSel_new_TT",      Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_Tlike   = new TH1F("hfourlepbestmass_4l_afterSel_new_Tlike",   "hfourlepbestmass_4l_afterSel_new_Tlike",   Nbins, Xmin, Xmax);

  TH1F *hfourlepbestmass_4l_afterSel_new_Wj = new TH1F("hfourlepbestmass_4l_afterSel_new_Wj", "hfourlepbestmass_4l_afterSel_new_Wj",Nbins, Xmin, Xmax);

  // Higgs as background
  TH1F *hfourlepbestmass_4l_afterSel_new_ggH  = new TH1F("hfourlepbestmass_4l_afterSel_new_ggH",  "hfourlepbestmass_4l_afterSel_new_ggH",  Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_VH   = new TH1F("hfourlepbestmass_4l_afterSel_new_VH",   "hfourlepbestmass_4l_afterSel_new_VH",   Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_VBFH = new TH1F("hfourlepbestmass_4l_afterSel_new_VBFH", "hfourlepbestmass_4l_afterSel_new_VBFH", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_afterSel_new_ttH  = new TH1F("hfourlepbestmass_4l_afterSel_new_ttH",  "hfourlepbestmass_4l_afterSel_new_ttH",  Nbins, Xmin, Xmax);

  for (int datasetId=Vdatasetnamebkg.size()-1; datasetId >=0; datasetId--) {
    char dataset[328];
    sprintf(dataset,"%s",Vdatasetnamebkg.at(datasetId).c_str());
    //std:cout << "Root-ple= " << dataset << "N entries= " <<  hfourlepbestmass_4l_afterSel_new->GetEntries() << std::endl;
    std::cout << "Counter=" << datasetId << " Root-ple=" << dataset << " Label=" << Vlabelbkg.at(datasetId) <<endl;

    std::string datasetnamebkg = "";
    datasetnamebkg = Vdatasetnamebkg.at(datasetId);

    TFile *f2 = TFile::Open(dataset);
    hfourlepbestmass_4l_afterSel_new = (TH1F*)f2->Get(histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new"*/);
    //hfourlepbestmass_4l_afterSel_new->SetBins(3000,0.,3000.);
    //hfourlepbestmass_4l_afterSel_new->SetBins(2940,60.,3000.);

    TH1 *hfourlepbestmass_4l_afterSel_new_new;

    if (datasetnamebkg.find("WZ")              < 200 ||
	datasetnamebkg.find("WW")              < 200 ||
	datasetnamebkg.find("ZToEE")           < 200 ||
	datasetnamebkg.find("ZToMuMu")         < 200 ||
	// datasetnamebkg.find("DYJetsToLL")    < 200 ||
	// datasetnamebkg.find("DYtoMuMu")      < 200 ||
	datasetnamebkg.find("DYlightJetsToLL") < 200 ||
	datasetnamebkg.find("DYbbJetsToLL")    < 200 ||
	datasetnamebkg.find("DYccJetsToLL")    < 200 ||
	// datasetnamebkg.find("DYToEE")         < 200 ||
	// datasetnamebkg.find("DYToMuMu")       < 200 ||
	// datasetnamebkg.find("TTTo2L2Nu")      < 200 ||
	// datasetnamebkg.find("TTJets")         < 200 ||
	// datasetnamebkg.find("TT_Tune")        < 200 ||
	// datasetnamebkg.find("TTbar")          < 200 ||
	datasetnamebkg.find("TTTo2L2Nu")       < 200 ||
	datasetnamebkg.find("TTToLL_MLL")      < 200 ||
	datasetnamebkg.find("WJetsToLNu")      < 200 ||
	datasetnamebkg.find("M-125")           < 200
	) {
      hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebin, histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
      hfourlepbestmass_4l_afterSel_new_new->SetLineColor(Vcolorbkg.at(datasetId)/*datasetId+2*/);
      hfourlepbestmass_4l_afterSel_new_new->SetFillColor(Vcolorbkg.at(datasetId)/*datasetId+2*/);
      hfourlepbestmass_4l_afterSel_new_new->SetMarkerStyle(24);
      hfourlepbestmass_4l_afterSel_new_new->SetLineWidth(1);

      std::cout << "(l644)Label= " << Vlabelbkg.at(datasetId)
		<< "  Entries= "   << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
		<< "  Integral= "  << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1);

      if (hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0)
	std::cout << "  Error= " << sqrt(hfourlepbestmass_4l_afterSel_new_new->GetEntries())*hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1)/hfourlepbestmass_4l_afterSel_new_new->GetEntries();
      std::cout << std::endl;

      // Higgs as background
      if (datasetnamebkg.find("GluGluToHToZZ") < 200) {
        std::cout << "ggH" << std::endl;
	hfourlepbestmass_4l_afterSel_new_ggH->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_ggH->SetMarkerColor(kOrange-3);
	hfourlepbestmass_4l_afterSel_new_ggH->SetFillColor(kOrange-3);

	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if (datasetnamebkg.find(temp) < 200 &&
	   (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	   hfourlepbestmass_4l_afterSel_new_ggH->GetEntries() > 0.)
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
	//hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
      }
      if (datasetnamebkg.find("TTbarH") < 200) {
        std::cout << "ttH" << std::endl;
	hfourlepbestmass_4l_afterSel_new_ttH->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_ttH->SetMarkerColor(kOrange);
	hfourlepbestmass_4l_afterSel_new_ttH->SetFillColor(kOrange);

	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if (datasetnamebkg.find(temp) < 200 &&
	    (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	    hfourlepbestmass_4l_afterSel_new_ttH->GetEntries() > 0.)
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
        //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
      }
      if (datasetnamebkg.find("VBF") < 200) {
        std::cout << "VBF" << std::endl;
	hfourlepbestmass_4l_afterSel_new_VBFH->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_VBFH->SetMarkerColor(kOrange-2);
	hfourlepbestmass_4l_afterSel_new_VBFH->SetFillColor(kOrange-2);

	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if (datasetnamebkg.find(temp) < 200 &&
	    (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	    hfourlepbestmass_4l_afterSel_new_VBFH->GetEntries() > 0.)
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
	//hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
      }

      if (datasetnamebkg.find("WH_ZH") < 200) {
        std::cout << "WH_ZH" << std::endl;
	hfourlepbestmass_4l_afterSel_new_VH->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_VH->SetMarkerColor(kOrange-1);
	hfourlepbestmass_4l_afterSel_new_VH->SetFillColor(kOrange-1);

	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if (datasetnamebkg.find(temp) < 200 &&
	    (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	    hfourlepbestmass_4l_afterSel_new_VH->GetEntries() > 0.)
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
	//hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
       }

      // DYJetsToLL check normalization
      //if (datasetnamebkg.find("DYJetsToLL") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0) {

      // DYJetsToTauTau
      if (datasetnamebkg.find("DYJetsToTauTau") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0) {
        hfourlepbestmass_4l_afterSel_new_DYtautau->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_DYtautau->SetMarkerColor(kAzure+4);
        hfourlepbestmass_4l_afterSel_new_DYtautau->SetFillColor(kAzure+4);
	hfourlepbestmass_4l_afterSel_new_DYtautau->SetLineColor(kAzure-9);

	hfourlepbestmass_4l_afterSel_new_DY->Add(hfourlepbestmass_4l_afterSel_new_new);
        hfourlepbestmass_4l_afterSel_new_DY->SetMarkerColor(kAzure-9);
        hfourlepbestmass_4l_afterSel_new_DY->SetFillColor(kAzure-9);
	hfourlepbestmass_4l_afterSel_new_DY->SetLineColor(kAzure-9);

        char temp[328];
        sprintf(temp,"%s",histosdir.c_str());
        if (datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 && datasetnamebkg.find("DYJetsToTauTau") < 200)) {
          std::cout << "DY->tautau= " << hfourlepbestmass_4l_afterSel_new_DYtautau->Integral(0,-1) << std::endl;
          /*
	  if (m_useDYJets == true)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
	  */
        }
      }

      if (datasetnamebkg.find("ZToMuMu") < 200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0) {
	//if (datasetnamebkg.find("DYtoMuMu") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0) {
	//hfourlepbestmass_4l_afterSel_new_new->Scale(double(6532812.*9153492../36277961.)*double(hfourlepbestmass_4l_afterSel_new_new->GetEntries()/12145114./hfourlepbestmass_4l_afterSel_new_new->GetEntries()));
	//hfourlepbestmass_4l_afterSel_new_new->Scale(double(4710.*3048.*11974371./80767910.)*double(hfourlepbestmass_4l_afterSel_new_new->GetEntries()/11974371./hfourlepbestmass_4l_afterSel_new_new->GetEntries()));
	// hfourlepbestmass_4l_afterSel_new_new->Scale(double(4710.*3048.*12138430./36257961.)*double(hfourlepbestmass_4l_afterSel_new_new->GetEntries()/12138430./hfourlepbestmass_4l_afterSel_new_new->GetEntries()));
	hfourlepbestmass_4l_afterSel_new_DY->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DY->SetMarkerColor(kAzure-9);
	hfourlepbestmass_4l_afterSel_new_DY->SetFillColor(kAzure-9);
	hfourlepbestmass_4l_afterSel_new_DY->SetLineColor(kAzure-9);

	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	//	if (datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && (datasetnamebkg.find("DYJetsToLL_M-50_TuneZ2Star") < 200 || datasetnamebkg.find("DYJetsToLL_M-50") < 200)) {
	//if (datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 && datasetnamebkg.find("CMSSW8012_MC_reHLT_DYtoMuMu120to200") < 200)) {
	if (datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 && datasetnamebkg.find("ZToMuMu_NNPDF30_13TeV-powheg_M_120_200") < 200)) {
          std::cout << "DY= " << hfourlepbestmass_4l_afterSel_new_DY->Integral(0,-1) << std::endl;
	  if (m_useDYJets == true)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
	}
	//hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
      }

      if (datasetnamebkg.find("ZToEE") < 200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0) {
	//if (datasetnamebkg.find("DYtoEE") <200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0) {
	//hfourlepbestmass_4l_afterSel_new_new->Scale(double(6532812.*9153492../36277961.)*double(hfourlepbestmass_4l_afterSel_new_new->GetEntries()/12145114./hfourlepbestmass_4l_afterSel_new_new->GetEntries()));
	//hfourlepbestmass_4l_afterSel_new_new->Scale(double(4710.*3048.*11974371./80767910.)*double(hfourlepbestmass_4l_afterSel_new_new->GetEntries()/11974371./hfourlepbestmass_4l_afterSel_new_new->GetEntries()));
	// hfourlepbestmass_4l_afterSel_new_new->Scale(double(4710.*3048.*12138430./36257961.)*double(hfourlepbestmass_4l_afterSel_new_new->GetEntries()/12138430./hfourlepbestmass_4l_afterSel_new_new->GetEntries()));
	hfourlepbestmass_4l_afterSel_new_DY->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DY->SetMarkerColor(kAzure-9);
	hfourlepbestmass_4l_afterSel_new_DY->SetFillColor(kAzure-9);
	hfourlepbestmass_4l_afterSel_new_DY->SetLineColor(kAzure-9);

	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	//	if (datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && (datasetnamebkg.find("DYJetsToLL_M-50_TuneZ2Star") < 200 || datasetnamebkg.find("DYJetsToLL_M-50") < 200)) {
	//if (datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 && datasetnamebkg.find("CMSSW8012_MC_reHLT_DYtoEE120to200") < 200)) {
	if (datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 && datasetnamebkg.find("ZToEE_NNPDF30_13TeV-powheg_M_120_200") < 200)) {
          std::cout << "DY= " << hfourlepbestmass_4l_afterSel_new_DY->Integral(0,-1) << std::endl;
	  if (m_useDYJets == true)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
	}
	//hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
      }

      // DYlightJetsToLL check normalization
      if (datasetnamebkg.find("DYlightJetsToLL") <200 &&
	  hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0) {
	hfourlepbestmass_4l_afterSel_new_DYlight->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DYlight->SetMarkerColor(kAzure+6);
	//hfourlepbestmass_4l_afterSel_new_DYlight->SetLineColor(kAzure+6);
	//hfourlepbestmass_4l_afterSel_new_DYlight->SetLineWidth(2);
	hfourlepbestmass_4l_afterSel_new_DYlight->SetFillColor(kAzure+6);

	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if (datasetnamebkg.find(temp) <200 && (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) && datasetnamebkg.find("DYlightJetsToLL_TuneZ2_M-50") < 200) {
          std::cout << "DYlight= " << hfourlepbestmass_4l_afterSel_new_DYlight->Integral(0,-1) << std::endl;
	  if (m_useDYJets == false)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
        }
	//hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
      }

      // DYbb
      if (datasetnamebkg.find("DYbbJetsToLL") <200 &&
	  hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0) {
	hfourlepbestmass_4l_afterSel_new_DYbb->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DYbb->SetMarkerColor(kAzure+2);
	//	hfourlepbestmass_4l_afterSel_new_DYbb->SetLineColor(kAzure+2);
	//	hfourlepbestmass_4l_afterSel_new_DYbb->SetLineWidth(2);
	hfourlepbestmass_4l_afterSel_new_DYbb->SetFillColor(kAzure+2);

	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if (datasetnamebkg.find(temp) <200 &&
	    (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	    datasetnamebkg.find("DYbbJetsToLL_TuneZ2_M-50") < 200) {
          std::cout << "DYbb= " << hfourlepbestmass_4l_afterSel_new_DYbb->Integral(0,-1) << std::endl;
	  if (m_useDYJets == false)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
        }
      }

      //DYCC
      if (datasetnamebkg.find("DYccJetsToLL") < 200 && hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0) {

	hfourlepbestmass_4l_afterSel_new_DYcc->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_DYcc->SetMarkerColor(kRed+0);
	// hfourlepbestmass_4l_afterSel_new_DYcc->SetLineColor(kRed+0);
	// hfourlepbestmass_4l_afterSel_new_DYcc->SetLineWidth(2);
       	hfourlepbestmass_4l_afterSel_new_DYcc->SetFillColor(kRed+0);
	char temp[328];
	sprintf(temp,"%s",histosdir.c_str());
	if (datasetnamebkg.find(temp) <200 &&
	    (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	    datasetnamebkg.find("DYccJetsToLL_M-50_TuneZ2Star") < 200) {
          std::cout << "DYcc= " << hfourlepbestmass_4l_afterSel_new_DYcc->Integral(0,-1) << std::endl;
	  if (m_useDYJets == false)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
        }
      }

      if (m_useDYJetsFromData == false) {
	// WW
	if (datasetnamebkg.find("WW") < 200) {
	  hfourlepbestmass_4l_afterSel_new_WW->Add(hfourlepbestmass_4l_afterSel_new_new);
	  hfourlepbestmass_4l_afterSel_new_WW->SetMarkerColor(kCyan+3);
	  hfourlepbestmass_4l_afterSel_new_WW->SetFillColor(kCyan+3);

	  hfourlepbestmass_4l_afterSel_new_diBoson->Add(hfourlepbestmass_4l_afterSel_new_new);
	  hfourlepbestmass_4l_afterSel_new_diBoson->SetMarkerColor(kGreen+2);
          hfourlepbestmass_4l_afterSel_new_diBoson->SetFillColor(kGreen+2);
	  hfourlepbestmass_4l_afterSel_new_diBoson->SetLineColor(kGreen+2);


	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  /*if (datasetnamebkg.find(temp) < 200 &&
	    (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	    hfourlepbestmass_4l_afterSel_new_WW->GetEntries() > 0.)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");*/

	  if (datasetnamebkg.find("WWTo2L2Nu_13TeV") < 200 &&
	      datasetnamebkg.find(temp) < 200 &&
	      (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	      hfourlepbestmass_4l_afterSel_new_diBoson->GetEntries() > 0.)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_diBoson,"WW, WZ, ZZ", "F");

	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
	}

	// WZ
	if (datasetnamebkg.find("WZ") < 200) {
	  hfourlepbestmass_4l_afterSel_new_WZ->Add(hfourlepbestmass_4l_afterSel_new_new);
	  hfourlepbestmass_4l_afterSel_new_WZ->SetMarkerColor(kCyan-2);
	  hfourlepbestmass_4l_afterSel_new_WZ->SetFillColor(kCyan-2);

	  hfourlepbestmass_4l_afterSel_new_diBoson->Add(hfourlepbestmass_4l_afterSel_new_new);
	  hfourlepbestmass_4l_afterSel_new_diBoson->SetMarkerColor(kGreen+2);
          hfourlepbestmass_4l_afterSel_new_diBoson->SetFillColor(kGreen+2);
	  hfourlepbestmass_4l_afterSel_new_diBoson->SetLineColor(kGreen+2);

	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  /*if (datasetnamebkg.find(temp) < 200 &&
	    (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	    hfourlepbestmass_4l_afterSel_new_WZ->GetEntries() > 0.)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");*/
	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
	}

	// TTT
	//if (datasetnamebkg.find("TTJets") < 200) {
	if (datasetnamebkg.find("TTToLL_MLL") < 200 || datasetnamebkg.find("TTTo2L2Nu") < 200) {
	  // this dataset is improperly weighted due to the way someone set weights to 0...
	  hfourlepbestmass_4l_afterSel_new_TT->Add(hfourlepbestmass_4l_afterSel_new_new);
	  hfourlepbestmass_4l_afterSel_new_TT->SetMarkerColor(kTeal-6);
	  hfourlepbestmass_4l_afterSel_new_TT->SetFillColor(kTeal-6);

	  hfourlepbestmass_4l_afterSel_new_Tlike->Add(hfourlepbestmass_4l_afterSel_new_new);
	  hfourlepbestmass_4l_afterSel_new_Tlike->SetFillColor(kRed);
	  hfourlepbestmass_4l_afterSel_new_Tlike->SetLineColor(kRed);


	  std::cout << "TT+jets= " << hfourlepbestmass_4l_afterSel_new_TT->GetEntries() << std::endl;
	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  /*if (datasetnamebkg.find(temp) < 200 &&
	    (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	    hfourlepbestmass_4l_afterSel_new_TT->GetEntries() > 0.)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");*/
	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
	  if (datasetnamebkg.find("TTTo2L2Nu_Tune") < 200  &&
	      datasetnamebkg.find(temp) < 200 &&
	      (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	      hfourlepbestmass_4l_afterSel_new_Tlike->GetEntries() > 0.)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_Tlike,"t#bar{t}, tW, #bar{t}W", "F");
	}

	// W+jets
	if (datasetnamebkg.find("_WJ") < 200 && m_useWJetsFromFakeRateFromMC == false) {
	  std::cout << "Wjets" << std::endl;
	  hfourlepbestmass_4l_afterSel_new_Wj->Add(hfourlepbestmass_4l_afterSel_new_new);
	  hfourlepbestmass_4l_afterSel_new_Wj->SetMarkerColor(kSpring);
	  hfourlepbestmass_4l_afterSel_new_Wj->SetFillColor(kSpring);

	  char temp[328];
	  sprintf(temp,"%s",histosdir.c_str());
	  if (datasetnamebkg.find(temp) < 200 &&
	      (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	      hfourlepbestmass_4l_afterSel_new_Wj->GetEntries() > 0.)
	    legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
	  //hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
	}
      }
    }
    //
    else if (datasetnamebkg.find("GluGluToZZTo") < 200) {
      std::cout << "Adding sample ggZZ" << std::endl;
      hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
      if (hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0.)
	errorZZ=errorZZ+pow(sqrt(hfourlepbestmass_4l_afterSel_new_new->GetEntries())*hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1)/hfourlepbestmass_4l_afterSel_new_new->GetEntries(),2);
      //std:cout << sqrt(errorZZ) << std::endl;
      hfourlepbestmass_4l_afterSel_new_ggZZ->Add(hfourlepbestmass_4l_afterSel_new_new);
      //hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001);
      hfourlepbestmass_4l_afterSel_new_ggZZ->SetLineColor(1);
      hfourlepbestmass_4l_afterSel_new_ggZZ->SetFillColor(kPink+5);
      hfourlepbestmass_4l_afterSel_new_ggZZ->SetLineWidth(1);

      hfourlepbestmass_4l_afterSel_new_ZZ->Add(hfourlepbestmass_4l_afterSel_new_new);
      //hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001);
      hfourlepbestmass_4l_afterSel_new_ZZ->SetLineColor(1);
      hfourlepbestmass_4l_afterSel_new_ZZ->SetFillColor(ZZBgColor);
      hfourlepbestmass_4l_afterSel_new_ZZ->SetLineWidth(1);

      hfourlepbestmass_4l_afterSel_new_diBoson->Add(hfourlepbestmass_4l_afterSel_new_new);
      hfourlepbestmass_4l_afterSel_new_diBoson->SetMarkerColor(kGreen+2);
      hfourlepbestmass_4l_afterSel_new_diBoson->SetFillColor(kGreen+2);
      hfourlepbestmass_4l_afterSel_new_diBoson->SetLineColor(kGreen+2);

      std::cout << "(l895)Label= " << Vlabelbkg.at(datasetId)
		<< "  Entries= "   << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
		<< "  Integral= "  << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1);

      if (hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0.)
	std::cout << "  Error= " << sqrt(hfourlepbestmass_4l_afterSel_new_new->GetEntries())*hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1)/hfourlepbestmass_4l_afterSel_new_new->GetEntries();
      std::cout <<endl;
      //legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");

      if (datasetnamebkg.find("GluGluToZZTo2L2L") < 200 &&
	  (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200)) {
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_ggZZ,"ggZZ", "F");
	std::cout << "(l907)Label= ggZZ      Integral= " << hfourlepbestmass_4l_afterSel_new_ggZZ->Integral(0,-1) << std::endl;
	std::cout << "(l907)Label= Total ZZ  Integral= " << hfourlepbestmass_4l_afterSel_new_ZZ->Integral(0,-1)   << std::endl;
      }
    }
    //else if ( datasetnamebkg.find("_ZZTo") < 200) {
    else if ( datasetnamebkg.find("_ZZ_") < 200) {
      std::cout << "Adding sample qqZZ" << std::endl;
      hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
      if (hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0.)
	errorZZ=errorZZ+pow(sqrt(hfourlepbestmass_4l_afterSel_new_new->GetEntries())*hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1)/hfourlepbestmass_4l_afterSel_new_new->GetEntries(),2);
      //std:cout << sqrt(errorZZ) << std::endl;
      hfourlepbestmass_4l_afterSel_new_qqZZ->Add(hfourlepbestmass_4l_afterSel_new_new);
      //hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001);
      hfourlepbestmass_4l_afterSel_new_qqZZ->SetLineColor(1);
      hfourlepbestmass_4l_afterSel_new_qqZZ->SetFillColor(kPink+5);
      hfourlepbestmass_4l_afterSel_new_qqZZ->SetLineWidth(1);

      hfourlepbestmass_4l_afterSel_new_ZZ->Add(hfourlepbestmass_4l_afterSel_new_new);
      //hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001);
      hfourlepbestmass_4l_afterSel_new_ZZ->SetLineColor(1);
      hfourlepbestmass_4l_afterSel_new_ZZ->SetFillColor(ZZBgColor);
      hfourlepbestmass_4l_afterSel_new_ZZ->SetLineWidth(1);

      hfourlepbestmass_4l_afterSel_new_diBoson->Add(hfourlepbestmass_4l_afterSel_new_new);
      hfourlepbestmass_4l_afterSel_new_diBoson->SetMarkerColor(kGreen);
      hfourlepbestmass_4l_afterSel_new_diBoson->SetFillColor(kGreen);
      hfourlepbestmass_4l_afterSel_new_diBoson->SetLineColor(kGreen+2);
      hfourlepbestmass_4l_afterSel_new_diBoson->SetLineWidth(1);

      char temp[328];
      if (whichsample.find("8TeV") < 200)
	sprintf(temp,"%s/output_ZZTo2e2mu_%s",histosdir.c_str(),whichsample.c_str());
      else if (whichsample.find("7TeV") < 200)
	sprintf(temp,"%s/output_ZZTo2e2mu_mll4_%s",histosdir.c_str(),whichsample.c_str());
      /*else if (whichsample.find("13TeV") < 200)
	sprintf(temp,"%s/output_ZZTo4L_%s",histosdir.c_str(),whichsample.c_str());*/
      else if (whichsample.find("13TeV") < 200)
	sprintf(temp,"%s/output_ZZ_%s",histosdir.c_str(),whichsample.c_str());

      if (datasetnamebkg.find(temp) < 200) {
	//std:cout << "Ciao" << std::endl;
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_qqZZ,"ZZ", "F");
	if (errorZZ>0.) errorZZ=sqrt(errorZZ);
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_ZZ,"Z#gamma^{*}, ZZ", "F");
	std::cout << "(l590)Label= qqZZ      Integral= " << hfourlepbestmass_4l_afterSel_new_qqZZ->Integral(0,-1) << std::endl;
        std::cout << "(l590)Label= qqZZ+ggZZ Integral= " << hfourlepbestmass_4l_afterSel_new_ZZ->Integral(0,-1)
		  << "  Error= " << errorZZ << std::endl;
	if (histlabel.find("hM4l_8") < 10)
	  outputyields << "ZZ " << hfourlepbestmass_4l_afterSel_new_ZZ->Integral(0,-1) << " +/- " << errorZZ << std::endl;

	int bin100=0,bin800=0;
	bool b_bin100=0,b_bin800=0;
	for (int nbins=1;nbins<=hfourlepbestmass_4l_afterSel_new_ZZ->GetNbinsX(); nbins++) {
	  if (hfourlepbestmass_4l_afterSel_new_ZZ->GetBinCenter(nbins) >=100. && b_bin100 == false) {
	    b_bin100=true;
	    bin100=nbins;
	  }
	  if (hfourlepbestmass_4l_afterSel_new_ZZ->GetBinCenter(nbins) >=800. && b_bin800 == false) {
	    b_bin800=true;
	    bin800=nbins;
	  }
	}
	std::cout << "Label= qqZZ+ggZZ (mll>100 and mll<800)  bins=" << bin100 << "-" << bin800
		  << " Integral= " << hfourlepbestmass_4l_afterSel_new_ZZ->Integral(bin100,bin800) << std::endl;
      }

      std::cout << "(l972)Label= " << Vlabelbkg.at(datasetId)
		<< "  Entries= "   << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
		<< "  Integral= "  << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1);
      if (hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0.)
	std::cout << "  Error= " << sqrt(hfourlepbestmass_4l_afterSel_new_new->GetEntries())*hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1)/hfourlepbestmass_4l_afterSel_new_new->GetEntries();
      std::cout << std::endl;

      //legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
    }
    else if (m_useDYJetsFromData == false &&  m_useDiJetsFromFakeRateFromData == false) {
      if (datasetnamebkg.find("_MuPt5Enriched") < 200) {
	hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
	hfourlepbestmass_4l_afterSel_new_qcdMu->Add(hfourlepbestmass_4l_afterSel_new_new);
	//hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001);
	hfourlepbestmass_4l_afterSel_new_qcdMu->SetLineColor(1);
	hfourlepbestmass_4l_afterSel_new_qcdMu->SetFillColor(kTeal-8);
	hfourlepbestmass_4l_afterSel_new_qcdMu->SetLineWidth(1);

	char temp[328];
	sprintf(temp,"%s/output_QCD_Pt-15to20_MuPt5Enriched",histosdir.c_str());

	if (datasetnamebkg.find(temp) < 200) { // provided that this is the last single-top sample
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_qcdMu,"QCD MuPt5", "F");
	  std::cout << "Label= QCD MuPt5    Integral= " << hfourlepbestmass_4l_afterSel_new_qcdMu->Integral(0,-1) << std::endl;
	}
	std::cout << "(l997)Label= " << Vlabelbkg.at(datasetId)
		  << "  Entries= "   << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
		  << "  Integral= "  << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) << std::endl;
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
      }
      else if (datasetnamebkg.find("_doubleEMEnriched") < 200) {
	hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
	hfourlepbestmass_4l_afterSel_new_qcdDEM->Add(hfourlepbestmass_4l_afterSel_new_new);
	//hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001);
	hfourlepbestmass_4l_afterSel_new_qcdDEM->SetLineColor(1);
	hfourlepbestmass_4l_afterSel_new_qcdDEM->SetFillColor(kTeal+8);
	hfourlepbestmass_4l_afterSel_new_qcdDEM->SetLineWidth(1);

	char temp[328];
	sprintf(temp,"%s/output_QCD_Pt-80_doubleEMEnriched",histosdir.c_str());

	if (datasetnamebkg.find(temp) < 80) {
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_qcdDEM,"QCD doubleEM", "F");
	  std::cout << "Label= QCD double EM    Integral= " << hfourlepbestmass_4l_afterSel_new_qcdMu->Integral(0,-1) <<endl;
	}
	std::cout << "(l1017)Label= " << Vlabelbkg.at(datasetId)
		  << "  Entries= "    << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
		  << "  Integral= "   << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) << std::endl;
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
      }
      else if (datasetnamebkg.find("_BCtoE") < 200) {
	hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
	hfourlepbestmass_4l_afterSel_new_qcdBC->Add(hfourlepbestmass_4l_afterSel_new_new);
	//hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001);
	hfourlepbestmass_4l_afterSel_new_qcdBC->SetLineColor(1);
	hfourlepbestmass_4l_afterSel_new_qcdBC->SetFillColor(kTeal-2);
	hfourlepbestmass_4l_afterSel_new_qcdBC->SetLineWidth(1);

	char temp[328];
	sprintf(temp,"%s/output_QCD_Pt-20to30_BCtoE",histosdir.c_str());

	if (datasetnamebkg.find(temp) < 200) { // provided that this is the last single-top sample
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_qcdBC,"QCD BCtoE", "F");
	  std::cout << "Label= QCD BCtoE    Integral= " << hfourlepbestmass_4l_afterSel_new_qcdBC->Integral(0,-1) << std::endl;
	}
	std::cout << "(l1037)Label= " << Vlabelbkg.at(datasetId)
		  << "  Entries= "    << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
		  << "  Integral= "   << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) << std::endl;
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
      }
      else if (datasetnamebkg.find("QCD_Pt") < 200) {

	hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
	hfourlepbestmass_4l_afterSel_new_qcd->Add(hfourlepbestmass_4l_afterSel_new_new);
	//hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(1001);
	hfourlepbestmass_4l_afterSel_new_qcd->SetLineColor(kTeal-2);
	hfourlepbestmass_4l_afterSel_new_qcd->SetFillColor(kTeal-2);
	hfourlepbestmass_4l_afterSel_new_qcd->SetLineWidth(1);

	char temp[328];
	sprintf(temp,"%s/output_QCD_Pt_1000to1400",histosdir.c_str());

	cout << "alpha" << temp << datasetnamebkg.find(temp) << std::endl;
	//sprintf(temp,"%s",histosdir.c_str());
	/*if (datasetnamebkg.find(temp) < 200 &&
	  (datasetnamebkg.find(whichenergy) < 200 || datasetnamebkg.find(whichsample) < 200) &&
	  hfourlepbestmass_4l_afterSel_new_Wj->GetEntries() > 0.)
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");*/

	if (datasetnamebkg.find(temp) < 200) { // provided that this is the last single-top sample
	  leg0->AddEntry(hfourlepbestmass_4l_afterSel_new_qcd,"QCD", "F");
	  legend->AddEntry(hfourlepbestmass_4l_afterSel_new_qcd,"QCD", "F");
	  std::cout << "Label= QCD Integral= " << hfourlepbestmass_4l_afterSel_new_qcd->Integral(0,-1) << std::endl;
	}
	std::cout << "(l1063)Label= " << Vlabelbkg.at(datasetId)
		  << "  Entries= "    << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
		  << "  Integral= "   << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) << std::endl;
	//legend->AddEntry(hfourlepbestmass_4l_afterSel_new_new,Vlabelbkg.at(datasetId).c_str(), "F");
      }

      // single top
      //else if (datasetId>=12 && datasetId<=14) {
      else if (datasetnamebkg.find("ST_") < 200) {
	hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebin,histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new_new"*/);
	hfourlepbestmass_4l_afterSel_new_singlet->Add(hfourlepbestmass_4l_afterSel_new_new);
	// hfourlepbestmass_4l_afterSel_new_new->SetMarkerColor(datasetId+4);
	// hfourlepbestmass_4l_afterSel_new_new->SetMarkerStyle(26);
	// hfourlepbestmass_4l_afterSel_new_new->SetLineWidth(2);
	// hfourlepbestmass_4l_afterSel_new_new->Draw("sameP");
	hfourlepbestmass_4l_afterSel_new_singlet->SetLineColor(Vcolorbkg.at(datasetId)/*datasetId-12+5*/);
	hfourlepbestmass_4l_afterSel_new_singlet->SetFillColor(kViolet);
	// hfourlepbestmass_4l_afterSel_new_singlet->SetFillStyle(3004);
	hfourlepbestmass_4l_afterSel_new_singlet->SetLineWidth(1);

	hfourlepbestmass_4l_afterSel_new_Tlike->Add(hfourlepbestmass_4l_afterSel_new_new);
	hfourlepbestmass_4l_afterSel_new_Tlike->SetFillColor(kRed);
	hfourlepbestmass_4l_afterSel_new_Tlike->SetLineColor(kRed);

	char temp[328];
	sprintf(temp,"%s/output_ST_t-channel_antitop_",histosdir.c_str());

	if (datasetnamebkg.find(temp) < 200 && datasetnamebkg.find("ST_t-channel_antitop") < 200) { // provided that this is the last single-top sample
	  //hfourlepbestmass_4l_afterSel_new_singlet->Draw("sameP");
	  leg0->AddEntry(hfourlepbestmass_4l_afterSel_new_singlet,"Single Top", "F");
	  //legend->AddEntry(hfourlepbestmass_4l_afterSel_new_singlet,"Single Top", "F");
	  std::cout << "Label= single t     Integral= " << hfourlepbestmass_4l_afterSel_new_singlet->Integral(0,-1) << std::endl;
	}
	std::cout << "(l1096)Label= " << Vlabelbkg.at(datasetId)
		  << "  Entries= "    << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
		  << "  Integral= "   << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1) << std::endl;
      }
    }

    char tempp[328];
    if (whichsample.find("8TeV") < 200)
      sprintf(tempp,"%s/output_ZZTo2e2mu_%s",histosdir.c_str(),whichsample.c_str());
    else if (whichsample.find("7TeV") < 200)
      sprintf(tempp,"%s/output_ZZTo2e2mu_mll4_%s",histosdir.c_str(),whichsample.c_str());
    /*else if (whichsample.find("13TeV") < 200)
      sprintf(tempp,"%s/output_ZZTo4L_%s",histosdir.c_str(),whichsample.c_str());*/
    else if (whichsample.find("13TeV") < 200)
      sprintf(tempp,"%s/output_ZZ_%s",histosdir.c_str(),whichsample.c_str());

    //std:cout << "tempp is " << tempp << std::endl;

    if (datasetnamebkg.find(tempp) < 300) {
      std::cout << "Stacking ZZ " << tempp << std::endl;
      //htotal->Add(hfourlepbestmass_4l_afterSel_new_ZZ);
      //htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_ZZ);
    }
    else {

       char temppp[328];
       sprintf(temppp,"%s",histosdir.c_str());

       // Higgs bkg
       if (datasetnamebkg.find(temppp) < 200 &&
	   (datasetnamebkg.find("output_GluGluToHToZZ") < 200 ||
	    (datasetnamebkg.find("output_GluGluToHToZZ") < 200 && datasetnamebkg.find(whichenergy.c_str()) < 200)
	   )
	  ) {
	 std::cout << "Adding ggH" << std::endl;
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_ggH);
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_ggH);
       }
       if (datasetnamebkg.find(temppp) < 200 &&
	   (datasetnamebkg.find("output_WH_ZH") < 200 ||
	    (datasetnamebkg.find("output_WH_ZH") < 200 && datasetnamebkg.find(whichenergy.c_str()) < 200)
	   )
	  ) {
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_VH);
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_VH);
       }
       if (datasetnamebkg.find(temppp) < 200 &&
	   (datasetnamebkg.find("output_TTbarH") < 200 ||
	    (datasetnamebkg.find("output_TTbarH") < 200 && datasetnamebkg.find(whichenergy.c_str()) < 200)
	   )
	  ) {
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_ttH);
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_ttH);
       }
       if (datasetnamebkg.find(temppp) < 200 &&
	   (datasetnamebkg.find("output_VBF") < 200 ||
	    (datasetnamebkg.find("output_VBF") < 200 && datasetnamebkg.find(whichenergy.c_str()) < 200)
	   )
	  ) {
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_VBFH);
	 htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_VBFH);
       }

       // other backgrounds
       if (m_useDYJets == true) {
	 if (datasetnamebkg.find(temppp) < 200 &&
	     (
	      // datasetnamebkg.find("output_DYJetsToLL_M-50_TuneZ2Star") < 200 ||
	      //(datasetnamebkg.find("output_DYJetsToLL_M-50") < 200 && datasetnamebkg.find(whichenergy.c_str()) < 200) ||
	      datasetnamebkg.find("output_ZToMuMu_NNPDF30_13TeV-powheg_M_120_200") <200
	      //datasetnamebkg.find("output_CMSSW8012_MC_reHLT_DYtoMuMu120to200_13TeV") < 200
	     )
	    ) {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_DY);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_DY);
	 }

	 //if (datasetnamebkg.find(temppp) < 200 &&
	 //   (datasetnamebkg.find("output_DYJetToTauTau") < 200)
         //  ) {
         //  htotal->Add(hfourlepbestmass_4l_afterSel_new_DYtautau);
         //  htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_DYtautau);
         //}

       }
       else if (m_useDYJets == false) {
	 if (datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_DYlightJetsToLL_TuneZ2_M-50") < 200) {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_DYlight);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_DYlight);
	 }
	 if (datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_DYbbJetsToLL_TuneZ2_M-50") < 200) {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_DYbb);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_DYbb);
	 }
	 if (datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_DYccJetsToLL_TuneZ2_M-50") < 200) {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_DYcc);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_DYcc);
	 }
       }


       if (datasetnamebkg.find(temppp) < 200 &&
	   (datasetnamebkg.find("output_WWTo2L2Nu_13TeV") < 200 ||
	    (datasetnamebkg.find("output_WWTo2L2Nu_13TeV") < 200 && datasetnamebkg.find(whichenergy.c_str()) < 200)
	   )
	  ) {
	 std::cout << "Adding WW or diBoson" << std::endl;
	 //htotal->Add(hfourlepbestmass_4l_afterSel_new_WW);
	 //htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_WW);
         htotal->Add(hfourlepbestmass_4l_afterSel_new_diBoson);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_diBoson);
       }
       if (datasetnamebkg.find(temppp) < 200 &&
	   (datasetnamebkg.find("output_WZ") < 200 ||
	    (datasetnamebkg.find("output_WZ") < 200 && datasetnamebkg.find(whichenergy.c_str()) < 200)
	   )
	  ) {
	 std::cout << "Adding WZ" << std::endl;
	 //htotal->Add(hfourlepbestmass_4l_afterSel_new_WZ);
	 //htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_WZ);
       }
       if (datasetnamebkg.find(temppp) < 200 &&
	   (datasetnamebkg.find("output_TTTo2L2Nu") < 200 ||
	    (datasetnamebkg.find("output_TTTo2L2Nu") < 200 && datasetnamebkg.find(whichenergy.c_str()) < 200)
	   )
	  ) {
	 std::cout << "Adding TT" << std::endl;
	 //htotal->Add(hfourlepbestmass_4l_afterSel_new_TT);
	 //htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_TT);
	 htotal->Add(hfourlepbestmass_4l_afterSel_new_Tlike);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_Tlike);
       }

       //if (datasetnamebkg.find("GluGluToZZTo2L2L") < 200) // provided that this is the last single-top sample
       //htotal->Add(hfourlepbestmass_4l_afterSel_new_ggZZ);

       //if (datasetnamebkg.find(temppp) < 200 datasetnamebkg.find("output_ZZTo2e2mu_8TeV") < 200)
       //	htotal->Add(hfourlepbestmass_4l_afterSel_new_ZZ);
       ////htotal->Add(hfourlepbestmass_4l_afterSel_new_qqZZ);

       if (m_useDiJetsFromFakeRateFromData == true) {
	 std::cout << "Adding Dijets with FR method from data" << std::endl;
         htotal->Add(hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new_new);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new_new);
	 m_useDiJetsFromFakeRateFromData = false;
       } else if (m_useDiJetsFromFakeRateFromData == false) {
	 if (datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_QCD_Pt-15to20_MuPt5Enriched") < 200) {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_qcdMu);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_qcdMu);
	 }

	 if (datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_QCD_Pt-40_doubleEMEnriched") < 200) {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_qcdDEM);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_qcdDEM);
	 }

	 if (datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_QCD_Pt_20to30_BCtoE") < 200) {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_qcdBC);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_qcdBC);
	 }

	 if (datasetnamebkg.find(temppp) < 200 && datasetnamebkg.find("output_QCD_Pt_1000to1400") < 200) {
	   htotal->Add(hfourlepbestmass_4l_afterSel_new_qcd);
	   htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_qcd);
	 }
       }

       if (m_useWJetsFromFakeRateFromMC == false && datasetnamebkg.find(temppp) < 200 &&
	   (datasetnamebkg.find("output_WJ") < 200 ||
	    (datasetnamebkg.find("output_WJ") < 200 &&
	     datasetnamebkg.find(whichenergy.c_str()) < 200)
	   )
	  ) {
	 std::cout << "Stacking W+jets" << std::endl;
         htotal->Add(hfourlepbestmass_4l_afterSel_new_Wj);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_Wj);
       }
       if (m_useWJetsFromFakeRateFromMC == true) {
         htotal->Add(hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC_new_new);
         htotalHisto->Add(hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC_new_new);
         m_useWJetsFromFakeRateFromMC = false;
       }

       if (datasetnamebkg.find(temppp) < 200 &&
	   datasetnamebkg.find("output_ST_") < 200 &&
	   datasetnamebkg.find("t-channel_antitop") < 200) {
	 std::cout << "Adding Single Top" << std::endl;
	 //htotal->Add(hfourlepbestmass_4l_afterSel_new_singlet);
	 //htotalHisto->Add(hfourlepbestmass_4l_afterSel_new_singlet);
       }
    }
     // htotal->Add(hfourlepbestmass_4l_afterSel_new_new);
  }

  // Building the ratio
  htotalHistoRatio->Divide(htotaldata,htotalHisto,1.,1.);
  /*
  for (int nbins=1;nbins<=htotaldata->GetNbinsX(); nbins++) {
    //std:cout << "Total: BinCenter=" << htotalHisto->GetBinCenter(nbins) << " BinContent=" << htotalHisto->GetBinContent(nbins) << " BinErrorContent=" << htotalHisto->GetBinError(nbins) << std::endl;
    if (htotalHisto->GetBinContent(nbins) > 0.) {
      htotalHistoRatio->SetBinContent(nbins,double(htotaldata->GetBinContent(nbins)/htotalHisto->GetBinContent(nbins)));
      //htotalHistoRatio->SetBinError(nbins,double(sqrt(htotaldata->GetBinContent(nbins))/htotalHisto->GetBinContent(nbins)));
      htotalHistoRatio->SetBinError(nbins,double(sqrt(
		 (1./(htotalHisto->GetBinContent(nbins)*htotalHisto->GetBinContent(nbins)))*htotaldata->GetBinContent(nbins) +
		 (htotaldata->GetBinContent(nbins)*htotaldata->GetBinContent(nbins)/pow(htotalHisto->GetBinContent(nbins),4))
		 *htotalHisto->GetBinContent(nbins)
		)));

    }
  }
  */


  // Signal  Z'prime and Contact Interaction
  TH1F *hfourlepbestmass_4l_afterSel_new_Zprime5000 = new TH1F("hfourlepbestmass_4l_afterSel_new_Zprime5000", "hfourlepbestmass_4l_afterSel_new_Zprime5000", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_CIToLL_M300  = new TH1F("hfourlepbestmass_4l_CIToLL_M300",  "hfourlepbestmass_4l_CIToLL_M300",  Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_CIToLL_M800  = new TH1F("hfourlepbestmass_4l_CIToLL_M800",  "hfourlepbestmass_4l_CIToLL_M800",  Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_CIToLL_M1300 = new TH1F("hfourlepbestmass_4l_CIToLL_M1300", "hfourlepbestmass_4l_CIToLL_M1300", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_CIToLL_M2000 = new TH1F("hfourlepbestmass_4l_CIToLL_M2000", "hfourlepbestmass_4l_CIToLL_M2000", Nbins, Xmin, Xmax);
  TH1F *hfourlepbestmass_4l_CIToLL       = new TH1F("hfourlepbestmass_4l_CIToLL",        "hfourlepbestmass_4l_CIToLL",      Nbins, Xmin, Xmax);

  std::string ciLFlav, ciLVal, ciIntf, ciHeli;
  for (int datasetIdSig=Vdatasetnamesig.size()-1; datasetIdSig >=0; datasetIdSig--) {

    char dataset[500];
    sprintf(dataset,"%s",Vdatasetnamesig.at(datasetIdSig).c_str());
    std::cout << "Counter=" << datasetIdSig << " Root-ple=" << dataset << " Label=" << Vlabelsig.at(datasetIdSig) <<endl;

    std::string datasetnamesig = "";
    datasetnamesig = Vdatasetnamesig.at(datasetIdSig);

    TFile *f0 = TFile::Open(dataset);
    hfourlepbestmass_4l_afterSel_new = (TH1F*)f0->Get(histlabel.c_str() /*"hfourlepbestmass_4l_afterSel_new"*/);
    //hfourlepbestmass_4l_afterSel_new->SetBins(3000,0.,3000.);
    //hfourlepbestmass_4l_afterSel_new->SetBins(2940,60.,3000.);

    TH1 *hfourlepbestmass_4l_afterSel_new_new = hfourlepbestmass_4l_afterSel_new->Rebin(nRebin, histlabel.c_str()/*"hfourlepbestmass_4l_afterSel_new_new"*/);

    //hfourlepbestmass_4l_afterSel_new_new->SetLineColor();
    //hfourlepbestmass_4l_afterSel_new_new->SetFillColor(0);
    //hfourlepbestmass_4l_afterSel_new_new->SetFillStyle(3244);
    hfourlepbestmass_4l_afterSel_new_new->SetMarkerSize(1.5);
    //hfourlepbestmass_4l_afterSel_new_new->Draw("same");

    std::cout << "(l1336)Label= " << Vlabelsig.at(datasetIdSig)
	      << "  Entries= "    << hfourlepbestmass_4l_afterSel_new_new->GetEntries()
	      << "  Integral= "   << hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1);
    ciLFlav = Vlabelsig.at(datasetIdSig).substr(Vlabelsig.at(datasetIdSig).find("TeV, m_{")+9,
						Vlabelsig.at(datasetIdSig).find("}=")-(Vlabelsig.at(datasetIdSig).find("TeV, m_{")+9));
    ciLVal = Vlabelsig.at(datasetIdSig).substr(Vlabelsig.at(datasetIdSig).find("#lambda")+7,
					       Vlabelsig.at(datasetIdSig).find(" TeV")-(Vlabelsig.at(datasetIdSig).find("#lambda")+7));
    ciIntf = Vlabelsig.at(datasetIdSig).substr(Vlabelsig.at(datasetIdSig).find("mode = ")+7,
					       Vlabelsig.at(datasetIdSig).find(", #eta")-(Vlabelsig.at(datasetIdSig).find("mode = ")+7));
    ciHeli = Vlabelsig.at(datasetIdSig).substr(Vlabelsig.at(datasetIdSig).find("#eta")+4,
					       2+(Vlabelsig.at(datasetIdSig).find("#eta")+4));
    if (hfourlepbestmass_4l_afterSel_new_new->GetEntries() > 0)
      std::cout << "  Error= " << sqrt(hfourlepbestmass_4l_afterSel_new_new->GetEntries())*hfourlepbestmass_4l_afterSel_new_new->Integral(0,-1)/hfourlepbestmass_4l_afterSel_new_new->GetEntries();
    std::cout << std::endl;

    double ci_factor = 1.0;
    if (datasetnamesig.find("100kTeV") < 200)
      ci_factor = -1.0;

    if ((datasetnamesig.find("CITo2Mu") < 200) || (datasetnamesig.find("CITo2E") < 200))
      hfourlepbestmass_4l_afterSel_new_new->Scale(1.3); // apply the k-factor

    if ((datasetnamesig.find("CITo2Mu_M300") < 200) || (datasetnamesig.find("CITo2E_M300") < 200)) {
      hfourlepbestmass_4l_CIToLL_M300->Add(hfourlepbestmass_4l_afterSel_new_new, ci_factor);
      hfourlepbestmass_4l_CIToLL->Add(hfourlepbestmass_4l_afterSel_new_new, ci_factor);
    }
    if ((datasetnamesig.find("CITo2Mu_M800") < 200) || (datasetnamesig.find("CITo2E_M800") < 200)) {
      hfourlepbestmass_4l_CIToLL_M800->Add(hfourlepbestmass_4l_afterSel_new_new, ci_factor);
      hfourlepbestmass_4l_CIToLL->Add(hfourlepbestmass_4l_afterSel_new_new, ci_factor);
    }
    if ((datasetnamesig.find("CITo2Mu_M1300") < 200) || (datasetnamesig.find("CITo2E_M1300") < 200)) {
      hfourlepbestmass_4l_CIToLL_M1300->Add(hfourlepbestmass_4l_afterSel_new_new, ci_factor);
      hfourlepbestmass_4l_CIToLL->Add(hfourlepbestmass_4l_afterSel_new_new, ci_factor);
    }
    if ((datasetnamesig.find("CITo2Mu_M2000") < 200) || (datasetnamesig.find("CITo2E_M2000") < 200)) {
      hfourlepbestmass_4l_CIToLL_M2000->Add(hfourlepbestmass_4l_afterSel_new_new, ci_factor);
      hfourlepbestmass_4l_CIToLL->Add(hfourlepbestmass_4l_afterSel_new_new, ci_factor);
    }

    if (datasetnamesig.find("Zprime") < 200 && datasetnamesig.find("5000") < 200)
      hfourlepbestmass_4l_afterSel_new_Zprime5000->Add(hfourlepbestmass_4l_afterSel_new_new);
  }

  // Zprime signal
  Double_t intErr = 0.;
  std::cout << "Zprime Signal expected at m_Z'=5000 GeV is " << hfourlepbestmass_4l_afterSel_new_Zprime5000->Integral(0,-1) << " +/- " << errorH125 << std::endl;
  if (histlabel.find("ZprimeRecomass") < 20)
    outputyields << "m_Z'=5 TeV " << hfourlepbestmass_4l_afterSel_new_Zprime5000->Integral(0,-1) << " +/- " << errorH125 << std::endl;
  hfourlepbestmass_4l_afterSel_new_Zprime5000->SetMarkerSize(0.95);
  hfourlepbestmass_4l_afterSel_new_Zprime5000->SetMarkerColor(kRed-4);
  hfourlepbestmass_4l_afterSel_new_Zprime5000->SetLineColor(kRed-4);

  // Contact Interaction
  intErr = 0.;
  std::cout << "CI Signal expected at m_X=300 GeV is " << hfourlepbestmass_4l_CIToLL_M300->IntegralAndError(0,-1,intErr)
	    << " +/- " << sqrt(hfourlepbestmass_4l_CIToLL_M300->GetEntries())*hfourlepbestmass_4l_CIToLL_M300->Integral(0,-1)/hfourlepbestmass_4l_CIToLL_M300->GetEntries()
	    << " +/- " << intErr // errors correct ???
	    << std::endl;
  if ((histlabel.find("CITo2Mu_M300") < 20) || (histlabel.find("CITo2E_M300") < 20))
    outputyields << "CI:m_X=300 GeV " << hfourlepbestmass_4l_CIToLL_M300->Integral(0,-1) << std::endl;

  intErr = 0.;
  std::cout << "CI Signal expected at m_X=800 GeV is " << hfourlepbestmass_4l_CIToLL_M800->IntegralAndError(0,-1,intErr)
	    << " +/- " << sqrt(hfourlepbestmass_4l_CIToLL_M800->GetEntries())*hfourlepbestmass_4l_CIToLL_M800->Integral(0,-1)/hfourlepbestmass_4l_CIToLL_M800->GetEntries()
	    << " +/- " << intErr // errors correct ???
	    << std::endl;
  if ((histlabel.find("CITo2Mu_M800") < 20) || (histlabel.find("CITo2E_M800") < 20))
    outputyields << "CI: m_X=800 GeV" << hfourlepbestmass_4l_CIToLL_M800->Integral(0,-1) << std::endl;

  intErr = 0.;
  std::cout << "CI Signal expected at m_X=1300 GeV is " << hfourlepbestmass_4l_CIToLL_M1300->IntegralAndError(0,-1,intErr)
	    << " +/- " << sqrt(hfourlepbestmass_4l_CIToLL_M1300->GetEntries())*hfourlepbestmass_4l_CIToLL_M1300->Integral(0,-1)/hfourlepbestmass_4l_CIToLL_M1300->GetEntries()
	    << " +/- " << intErr // errors correct ???
	    << std::endl;
  if ((histlabel.find("CITo2Mu_M1300") < 20) || (histlabel.find("CITo2E_M1300") < 20))
    outputyields << "CI: m_X=1300 GeV" << hfourlepbestmass_4l_CIToLL_M1300->Integral(0,-1) << std::endl;

  intErr = 0.;
  std::cout << "CI Signal expected at m_X=2000 GeV is " << hfourlepbestmass_4l_CIToLL_M2000->IntegralAndError(0,-1,intErr)
	    << " +/- " << sqrt(hfourlepbestmass_4l_CIToLL_M2000->GetEntries())*hfourlepbestmass_4l_CIToLL_M2000->Integral(0,-1)/hfourlepbestmass_4l_CIToLL_M2000->GetEntries()
	    << " +/- " << intErr // errors correct ???
	    << std::endl;
  if ((histlabel.find("CITo2Mu_M2000") < 20) || (histlabel.find("CITo2E_M2000") < 20))
    outputyields << "CI: m_X=2000 GeV" << hfourlepbestmass_4l_CIToLL_M2000->Integral(0,-1) << std::endl;

  hfourlepbestmass_4l_CIToLL_M300->SetMarkerSize(0.95);
  hfourlepbestmass_4l_CIToLL_M300->SetMarkerColor(kRed-3);
  hfourlepbestmass_4l_CIToLL_M300->SetLineColor(kRed-3);
  hfourlepbestmass_4l_CIToLL_M300->SetLineWidth(2);

  hfourlepbestmass_4l_CIToLL_M800->SetMarkerSize(0.95);
  hfourlepbestmass_4l_CIToLL_M800->SetMarkerColor(kRed-4);
  hfourlepbestmass_4l_CIToLL_M800->SetLineColor(kRed-4);
  hfourlepbestmass_4l_CIToLL_M800->SetLineWidth(2);

  hfourlepbestmass_4l_CIToLL_M1300->SetMarkerSize(0.95);
  hfourlepbestmass_4l_CIToLL_M1300->SetMarkerColor(kRed-7);
  hfourlepbestmass_4l_CIToLL_M1300->SetLineColor(kRed-7);
  hfourlepbestmass_4l_CIToLL_M1300->SetLineWidth(2);

  hfourlepbestmass_4l_CIToLL_M2000->SetMarkerSize(0.95);
  hfourlepbestmass_4l_CIToLL_M2000->SetMarkerColor(kRed-9);
  hfourlepbestmass_4l_CIToLL_M2000->SetLineColor(kRed-9);
  hfourlepbestmass_4l_CIToLL_M2000->SetLineWidth(2);

  intErr = 0.;
  std::cout << "Total CI Signal expected is " << hfourlepbestmass_4l_CIToLL->IntegralAndError(0,-1,intErr)
	    << " +/- " << sqrt(hfourlepbestmass_4l_CIToLL->GetEntries())*hfourlepbestmass_4l_CIToLL->Integral(0,-1)/hfourlepbestmass_4l_CIToLL->GetEntries() // errors correct ???
	    << " +/- " << intErr // errors correct ???
	    << std::endl;
  hfourlepbestmass_4l_CIToLL->SetMarkerSize(0.95);
  hfourlepbestmass_4l_CIToLL->SetMarkerColor(kRed-7);
  hfourlepbestmass_4l_CIToLL->SetLineColor(kRed-7);
  // hfourlepbestmass_4l_CIToLL->SetFillColor(kRed-7);
  hfourlepbestmass_4l_CIToLL->SetLineWidth(2);

  //htotal->Draw("hist same");
  //htotaldata->Draw("EPsame");
  //gr->Draw("EPsame");

  // Zprime
  //hfourlepbestmass_4l_afterSel_new_Zprime5000->Draw("same");
  //htotal->Add(hfourlepbestmass_4l_afterSel_new_Zprime5000);

  // Contact Interaction
  // hfourlepbestmass_4l_CIToLL_M300->Draw("same");
  // hfourlepbestmass_4l_CIToLL_M800->Draw("same");
  // hfourlepbestmass_4l_CIToLL_M1300->Draw("same");
  // hfourlepbestmass_4l_CIToLL_M2000->Draw("same");

  // htotal->Add(hfourlepbestmass_4l_CIToLL_M300);
  // htotal->Add(hfourlepbestmass_4l_CIToLL_M800);
  // htotal->Add(hfourlepbestmass_4l_CIToLL_M1300);
  // htotal->Add(hfourlepbestmass_4l_CIToLL_M2000);
  htotal->Add(hfourlepbestmass_4l_CIToLL);
  // why not get these from the Vlabelsig?? why are we doing things in the worst possible way????
  // Vlabelsig.at(datasetIdSig)
  std::cout << "Histlabel " << histlabel << std::endl;
  // legend->AddEntry(hfourlepbestmass_4l_CIToLL,TString("CITo2Mu #Lambda="+ciLVal+" TeV, #eta="+ciHeli+", "+ciIntf), "F");
  legend->AddEntry(hfourlepbestmass_4l_CIToLL,TString("CI#rightarrow"+ciLFlav+" #Lambda="+ciLVal+" TeV, #eta="+ciHeli+", "+ciIntf), "F");

  //hfourlepbestmass_4l_CIToLL_M800->Draw("same");
  //htotal->Add(hfourlepbestmass_4l_CIToLL_M800);
  //legend->AddEntry(hfourlepbestmass_4l_CIToLL_M800,"#Lambda = 10 TeV (m_{#mu^{+} #mu^{-}} = 800 GeV)", "L");

  if (histlabel.find("ZprimeRecomass") < 20) {
    std::cout << "Plotting di-lepton mass" << std::endl;
    //htotal->Add(hfourlepbestmass_4l_afterSel_new_Zprime5000);
    //legend->AddEntry(hfourlepbestmass_4l_afterSel_new_Zprime5000,"m_{Z'}=5 TeV", "L");
  }


  // Zoom
  // htotal->SetMinimum(0);
  //   htotal->SetMaximum(12.4); // 3 GeV bin
  //   htotal->Draw();
  //   htotal->GetXaxis()->SetRange(23,59); // 3 GeV bin
  //   //htotal->GetYaxis()->SetRangeUser(0.,16.);
  //   htotal->GetYaxis()->SetTitle("Events / 3 GeV"); // 3 GeV bin
  //   sprintf(histotitle,"m_{%s} [GeV]",whichchanneltex.c_str());
  //   htotal->GetXaxis()->SetTitle(histotitle);
  //   htotal->GetXaxis()->SetLabelSize(0.045);
  //   htotal->GetXaxis()->SetTitleSize(0.05);
  //   htotal->GetXaxis()->SetTitleOffset(1.15);
  //   htotal->GetXaxis()->SetTitleFont(42);
  //   htotal->GetYaxis()->SetLabelSize(0.045);
  //   htotal->GetYaxis()->SetTitleSize(0.05);
  //   htotal->GetYaxis()->SetTitleOffset(1.15);
  //   htotal->GetYaxis()->SetTitleFont(42);
  //   gr->Draw("EPsame");

  htotal->Draw("hist same");
  gr->Draw("E0P0same");

  //leg0->Draw("same");
  //leg1->Draw("same");

  legend->Draw("same");
  ll->Draw("same");

  std::string saveaspdfzoom = m_setLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom_log.pdf" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom.pdf";
  std::string saveaspngzoom = m_setLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom_log.png" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom.png";
  std::string saveasepszoom = m_setLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom_log.eps" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom.eps";
  std::string saveasrootzoom = m_setLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom_log.root" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_zoom.root";

  //   c1->SaveAs(saveaspdfzoom.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.pdf"*/);
  //   c1->SaveAs(saveaspngzoom.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.png"*/);
  //   c1->SaveAs(saveasepszoom.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.eps"*/);
  //   c1->SaveAs(saveasrootzoom.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.eps"*/);

  gPad->RedrawAxis();

  //  c1->Update();
  double canvasratio = 0.3;
  c1->SetBottomMargin(canvasratio + (1-canvasratio)*c1->GetBottomMargin()-canvasratio*c1->GetTopMargin());
  //std:cout << "Canvas= " << canvasratio + (1-canvasratio)*c1->GetBottomMargin()-canvasratio*c1->GetTopMargin() << std::endl;

  // Ratio: data / total bkg
  canvasratio = 0.16;
  TPad *ratioPad = new TPad("BottomPad","",0,0,1,1);
  ratioPad->SetTopMargin((1-canvasratio) - (1-canvasratio)*ratioPad->GetBottomMargin()+canvasratio*ratioPad->GetTopMargin());
  ratioPad->SetFillStyle(4000);
  ratioPad->SetFillColor(4000);
  ratioPad->SetFrameFillColor(4000);
  ratioPad->SetFrameFillStyle(4000);
  ratioPad->SetFrameBorderMode(0);
  ratioPad->SetTicks(1,1);
  ratioPad->SetGrid(1,1);
  if (m_setLogX)
    ratioPad->SetLogx();
  ratioPad->Draw();
  ratioPad->cd();

  //TH2F *hframe2= new TH2F("hframe2","hframe2",6000, 0., 2.2, 1000, 0.5, 2.);// iso

  hframe2->GetYaxis()->SetLabelSize(0.010);
  hframe2->GetXaxis()->SetLabelSize(0.010);
  hframe2->GetYaxis()->SetTitleSize(0.03);
  hframe2->GetXaxis()->SetTitleSize(0.03);
  //  hframe2->GetYaxis()->SetTitleSize(0.047);
  hframe2->SetYTitle("Data/MC");
  //  hframe2->GetYaxis()->SetRangeUser(-10,10);
  hframe2->GetYaxis()->SetNdivisions(503);
  //hframe2->GetXaxis()->SetTitleOffset(1.25);
  hframe2->Draw("");

  htotalHistoRatio->SetMarkerStyle(20);
  htotalHistoRatio->SetMarkerSize(0.95);
  htotalHistoRatio->SetMarkerColor(kBlack);
  htotalHistoRatio->Draw("Psame");

  c1->Update();

  std::string saveaspdfratio = m_setLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.pdf" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.pdf";
  std::string saveaspngratio = m_setLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.png" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.png";
  std::string saveasepsratio = m_setLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.eps" : "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.eps";
  std::string saveasrootratio= m_setLogY ? "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio_log.root": "plots/h_"+histlabel+"_"+whichchannel+"_"+whichenergy+"_ratio.root";
  std::cout << saveasrootratio.c_str() << std::endl;

  c1->SaveAs(saveaspdfratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.pdf"*/);
  c1->SaveAs(saveaspngratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.png"*/);
  c1->SaveAs(saveasepsratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.eps"*/);
  c1->SaveAs(saveasrootratio.c_str()/*"plots/hfourlepbestmass_4l_afterSel_new_m4l.root"*/);

  char htotal_root[300];
  if (nRebin == 20) {
    sprintf(htotal_root,"plots/htotal_root_%s.root",histlabel.c_str());
    TFile *file1 = new TFile(htotal_root, "RECREATE");
    file1->cd();
    htotaldata->Write();
    htotalHisto->Write();
    htotalHistoRatio->Write();
    hfourlepbestmass_4l_afterSel_new_DY->Write();
    hfourlepbestmass_4l_afterSel_new_diBoson->Write();
    hfourlepbestmass_4l_afterSel_new_Tlike->Write();
    if (m_useDiJetsFromFakeRateFromData)
      hfourlepbestmass_4l_afterSel_DiJetsFromFakeRateFromData_new_new->Write();
    if (m_useWJetsFromFakeRateFromMC)
      hfourlepbestmass_4l_afterSel_WJetsFromFakeRateFromMC_new_new->Write();
    if (m_useDiJetsWJetsFromFakeRateFromData)
      hfourlepbestmass_4l_afterSel_DiJetsWJetsFromFakeRateFromData_new_new->Write();
    file1->Write();
    file1->Close();
  }


}

TH1F* PlotStackZprime::DrawOverflow(TH1F *h)
{
  // This function paint the histogram h with an extra bin for overflows
  UInt_t nx    = h->GetNbinsX()+1;
  Double_t *xbins= new Double_t[nx+1];
  for (UInt_t i=0;i<nx;i++)
    xbins[i]=h->GetBinLowEdge(i+1);
  xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
  char *tempName= new char[strlen(h->GetName())+10];
  sprintf(tempName,"%swtOverFlow",h->GetName());
  // Book a temporary histogram having ab extra bin for overflows
  TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
  // Reset the axis labels
  htmp->SetXTitle(h->GetXaxis()->GetTitle());
  htmp->SetYTitle(h->GetYaxis()->GetTitle());
  // Fill the new hitogram including the extra bin for overflows
  for (UInt_t i=1; i<=nx; i++)
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
  // Fill the underflows
  htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());
  // FillStyle and color
  htmp->SetFillStyle(h->GetFillStyle());
  htmp->SetFillColor(h->GetFillColor());
  return htmp;
}

void PlotStackZprime::setSamplesNames4l()
{

  std::ifstream infile;
  infile.open(inputfile.c_str());

  if (inputfile.find("2011") < 100) {
    whichenergy="7TeV";
    whichsample="7TeV";
  } else if (inputfile.find("2012") < 100) {
    whichenergy="8TeV";
    whichsample="8TeV";
  } else if (inputfile.find("RunI") < 100) {
    whichenergy="RunI";
    whichsample="8TeV";
  } else if (inputfile.find("Spring16") < 100) {
    whichenergy="13TeV";
    whichsample="13TeV";
  }

  std::cout << "Doing plot for " << whichenergy.c_str() << "  and " << whichsample.c_str() << std::endl;

  if (inputfile.find("4l") < 25) {
    std::cout << "Plotting 4e+4mu+2e2mu combined in 4lepton plots" << std::endl;
    whichchannel="4l";
    whichchanneltex="4l";
    whichinvmasstex="4l";
    histosdir="histos4mu";
  } else if (inputfile.find("4mu") < 25) {
    std::cout << "Plotting 4mu" << std::endl;
    whichchannel="4mu";
    whichchanneltex="4#mu";
    whichinvmasstex="4#mu";
    histosdir="histos4mu";
  } else if (inputfile.find("4e") < 25) {
    std::cout << "Plotting 4e" << std::endl;
    whichchannel="4e";
    whichchanneltex="4e";
    whichinvmasstex="4e";
    histosdir="histos4e";
  } else if (inputfile.find("2e2mu") < 25) {
    std::cout << "Plotting 2e2mu" << std::endl;
    whichchannel="2e2mu";
    whichchanneltex="2e2#mu";
    whichinvmasstex="2e2#mu";
    histosdir="histos2e2mu";
  } else if (inputfile.find("zprime_SingleMuon") < 25) {
    std::cout << "Plotting Z->mumu" << std::endl;
    whichchannel="2mu";
    whichchanneltex="#mu#mu";
    whichinvmasstex="#mu^{+}#mu^{-}";
    histosdir="histos/histosZprimeMuMu";
  } else if (inputfile.find("zprime_DoubleEG") < 25) {
    std::cout << "Plotting Z->ee" << std::endl;
    whichchannel="2e";
    whichchanneltex="ee";
    whichinvmasstex="e^{+}e^{-}";
    histosdir="histos/histosZprimeEleEle";
  }

  Vdatasetnamedata.clear();
  Vdatasetnamebkg.clear();
  Vdatasetnamesig.clear();

  std::string inputfilename;

  while(std::getline(infile,inputfilename)) {
    // DATA
    if (inputfilename.find("#") < 1)
      continue;

    if (inputfilename.rfind("CITo2") != std::string::npos) {
      if ((inputfilename.rfind("M2000") != std::string::npos) && 
	  (ciSampleType.rfind("ConLL") == std::string::npos))
	continue;

      if (inputfilename.rfind("Lam100k") == std::string::npos) {
	// modify based on input requirement
	// lambda
	size_t charLp = ciSampleType.find("Lam");
	size_t charTp = ciSampleType.find("TeV"); //first instance
	std::string lVal = ciSampleType.substr(charLp,charTp-(charLp));
	std::string iVal = ciSampleType.substr(charTp+3,3);
	std::string hVal = ciSampleType.substr(charTp+6,2);

	if (inputfilename.find("Lam16") != std::string::npos) {
	  inputfilename.replace(inputfilename.find("Lam16"),std::string("Lam16").size(),lVal);
	}
	// interference
	if (inputfilename.find("Con") != std::string::npos) {
	  inputfilename.replace(inputfilename.find("Con"),std::string("Con").size(),iVal);
	}
	// helicity
	if (inputfilename.find("LL") != std::string::npos) {
	  inputfilename.replace(inputfilename.find("LL"),std::string("LL").size(),hVal);
	}
      }
    }

    std::cout << "Reading " << inputfilename.c_str() << std::endl;

    if (inputfilename.find("_SingleElectron") < 53) {
      Vdatasetnamedata.push_back(inputfilename);
      Vlabeldata.push_back("2016");
      Vxsectiondata.push_back(1.); //pb
    } else if (inputfilename.find("_SingleMuon") < 200) { // as many times as it occurs in the input file
      Vdatasetnamedata.push_back(inputfilename);
      Vlabeldata.push_back("2016");
      Vxsectiondata.push_back(1.); //pb
    } else if (inputfilename.find("_DoubleEG") < 200) {
      Vdatasetnamedata.push_back(inputfilename);
      Vlabeldata.push_back("2016");
      Vxsectiondata.push_back(1.); //pb
    } else if (inputfilename.find("Z+X") < 85) { //Z+X from data
      Vdatasetnamebkgdata.push_back(inputfilename);
      Vlabelbkgdata.push_back("Z+X");
      Vxsectionbkgdata.push_back(1.); //pb
    } else if (inputfilename.find("Data-total-jets") < 85) { //multijet from data
      Vdatasetnamebkgdata.push_back(inputfilename);
      Vlabelbkgdata.push_back("Jets");
      Vxsectionbkgdata.push_back(1.); //pb
    }

    // SIGNAL
    // need **way** more configurable and robust switching here...
    else if (inputfilename.find("ZprimeToMuMu_M-5000") < 200) {
      Vdatasetnamesig.push_back(inputfilename);
      Vlabelsig.push_back("m_{Z'}= 5 TeV/c^{2}");
      Vxsectionsig.push_back(1.); //pb
      Vcolorsig.push_back(kOrange-3);
    }
    // CI Sample looks like :CITo2{lepton}_M{mass}_CUETP8M1_Lam{lambda}TeV{interference}{helicity}_13TeV_Pythia8_Corrected-v4_ntuple
    // but also, don't want to have separate smples plotted...
    else if (inputfilename.find("CITo") < 200) {
      size_t charSp = inputfilename.find("CITo2");
      size_t charMp = inputfilename.find("_M");
      size_t charPp = inputfilename.find("_CUET");
      size_t charLp = inputfilename.find("_Lam");
      size_t charTp = inputfilename.find("TeV"); //first instance
      std::string lVal = inputfilename.substr(charLp+4,charTp-(charLp+4));
      std::string mVal = inputfilename.substr(charMp+2,(charPp)-(charMp+2));
      std::string iVal = inputfilename.substr(charTp+3,3);
      std::string hVal = inputfilename.substr(charTp+6,2);
      std::string lFlav = inputfilename.substr(charSp+6,charMp);

      Vdatasetnamesig.push_back(inputfilename);
      Vlabelsig.push_back("#lambda="+lVal+" TeV, m_{" + (lFlav == "EE" ? "e^{+}e^{-}" : "mu^{+}#mu^{-}")
			  + "}= "+mVal+" GeV, mode = "+iVal+(iVal == "Des" ? "tructive" :  "structive")+", #eta="+hVal);
      Vxsectionsig.push_back(1.); //pb
      Vcolorsig.push_back(kOrange-3);
    } else if (inputfilename.find("CITo2Mu_Lam10TeV_LLConM300") < 200) {
      Vdatasetnamesig.push_back(inputfilename);
      Vlabelsig.push_back("#lambda=10 TeV, m_{X}= 300 GeV");
      Vxsectionsig.push_back(1.); //pb
      Vcolorsig.push_back(kOrange-3);
    } else if (inputfilename.find("CITo2Mu_Lam10TeV_LLConM800") < 200) {
      Vdatasetnamesig.push_back(inputfilename);
      Vlabelsig.push_back("#lambda=10 TeV, m_{X}= 800 GeV");
      Vxsectionsig.push_back(1.); //pb
      Vcolorsig.push_back(kOrange-3);
    }

    // background higgs
    else if (inputfilename.find("GluGluToHToZZTo4L") < 100 ||
	    inputfilename.find("SMHiggsToZZTo4L") < 100   ||
	    inputfilename.find("VBF_HToZZTo4L") < 100     ||
	    inputfilename.find("VBF_ToHToZZTo4L") < 100   ||
	    inputfilename.find("TTbarH_HToZZTo4L") < 100  ||
	    inputfilename.find("WH_HToZZTo4L") < 100  ||
	    inputfilename.find("ZH_HToZZTo4L") < 100  ||
	    inputfilename.find("WH_ZH_HToZZ_4LFilter") < 100) {  // provided that signal samples contain 'GluGluHToZZTo4L'
      //Vdatasetnamebkg.push_back(inputfilename);


      if ( (inputfilename.find("M-125") < 200 && (inputfilename.find("GluGluToHToZZTo4L") < 200))) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("ggH, m_{H}=125 GeV/c^{2}");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kOrange-3);
	//std:cout << "ggH" << std::endl;
      }

      if ( (inputfilename.find("M-125") < 200 && (inputfilename.find("VBF") < 200))) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("VBF H, m_{H}=125 GeV/c^{2}");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kOrange-2);
	//std:cout << "VBF H" << std::endl;
      }

      if ( (inputfilename.find("M-125") < 200 && (inputfilename.find("WH_ZH") < 200))) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("VH, m_{H}=125 GeV/c^{2}");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kOrange-1);
	//std:cout << "VH" << std::endl;
      }

      if ( (inputfilename.find("M-125") < 200 && (inputfilename.find("TTbarH") < 200))) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("ttH, m_{H}=125 GeV/c^{2}");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kOrange);
	//std:cout << "ttH" << std::endl;
      }


      if ( (inputfilename.find("M-126") < 200 && !(inputfilename.find("GluGluToHToZZTo4L") < 200))) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("ggH, m_{H}=126 GeV/c^{2}");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kOrange-3);
      }

    }

    // BACKGROUND from other sources
    else{  // provided that everything that is neither data nor signal is background
      //Vdatasetnamebkg.push_back(inputfilename);

       // qqZZ to 4l
      if (inputfilename.find("_ZZ_TuneCUETP8M1_13TeV-pythia8") < 200 || inputfilename.find("_ZZ_TuneCUETP8M1_13TeV-pythia8") < 200 || inputfilename.find("_ZZ_") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("ZZ");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kPink+5);
      }

      if (inputfilename.find("_ZZTo4L_Tune4C_13TeV-powheg-pythia8") < 200 || inputfilename.find("_ZZTo4L_13TeV_powheg_pythia8") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("ZZ->4l");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kPink+5);
      }

      // qqZZ to 4mu
      if (inputfilename.find("_ZZTo4mu") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Z#gamma^{*},ZZ");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kPink+5);
      }

      // qqZZ to 4e
      if (inputfilename.find("_ZZTo4e") < 200) {
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("Z#gamma^{*],ZZ");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kPink+5);
      }

      // qqZZ to 4tau
      if (inputfilename.find("_ZZTo4tau") < 200) {
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("Z\\#gamma^{*],ZZ");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kPink+5);
      }

      // qqZZ to 2e2mu
      if (inputfilename.find("_ZZTo2e2mu") < 200) {
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("Z#gamma^{*],ZZ");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kPink+5);
      }

      // qqZZ to 2mu2tau
      if (inputfilename.find("_ZZTo2mu2tau") < 200) {
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("Z#gamma^{*],ZZ");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kPink+5);
      }

      // qqZZ to 2e2tau
      if (inputfilename.find("_ZZTo2e2tau") < 200) {
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("ZZ");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kPink+5);
      }



      // ggZZ 2L2L'
      if (inputfilename.find("GluGluToZZTo2L2L") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("GluGluToZZTo2L2L_8TeV-gg2zz-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kPink+5);
      }

      // ggZZ 4L
      if (inputfilename.find("GluGluToZZTo4L") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("GluGluToZZTo4L_8TeV-gg2zz-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kPink+5);
      }

      // WZ
      if (inputfilename.find("WZ") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("WZ");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kCyan-2);
      }

      // WWT2L2Nu
      if (inputfilename.find("WW") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("WW");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kCyan+3);
      }

      // DYJetsToLL_TuneZ2_M-50_8TeV-madgraph-tauola
      if (inputfilename.find("DYJetsToLL_TuneZ2_M-50") < 200 ||
	  inputfilename.find("DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Z+jets, m_{ll}>50 GeV");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kAzure+2);
      }

      // DYJetsToLL_M-10To50
      if (inputfilename.find("DYJetsToLL_M-10To50") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Z+jets");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kAzure+2);
      }


      // DYlighJetsToLL_TuneZ2_M-50_8TeV-madgraph-tauola
      if (inputfilename.find("DYlightJetsToLL_TuneZ2_M-50") < 200) {//DYlightJetsToLL_TuneZ2_M-50
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Zlight");
	Vxsectionbkg.push_back(1.);
	Vcolorbkg.push_back(kAzure+6);
      }

      // DYlighJetsToLL_TuneZ2_M-10To50_8TeV-madgraph-tauola
      if (inputfilename.find("DYlightJetsToLL_TuneZ2_M-10To50") < 200) {//DYlightJetsToLL_TuneZ2_M-10To50
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Zlight");
	Vxsectionbkg.push_back(1.);
	Vcolorbkg.push_back(kAzure+6);
      }


      // DYbbJetsToLL_TuneZ2_M-50_8TeV-madgraph-tauola
      if (inputfilename.find("DYbbJetsToLL_TuneZ2_M-50") < 200) {//DYbbJetsToLL_TuneZ2_M-50
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Zbb");
	Vxsectionbkg.push_back(1.);
	Vcolorbkg.push_back(kAzure+2);
      }

      // DYbbJetsToLL_TuneZ2_M-10To50
      if (inputfilename.find("DYbbJetsToLL_TuneZ2_M-10To50") < 200) {//DYbbJetsToLL_TuneZ2_M-10To50
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Zbb");
	Vxsectionbkg.push_back(1.);
	Vcolorbkg.push_back(kAzure+2);
      }

      // DYccJetsToLL_TuneZ2_M-50_8TeV-madgraph-tauola
      if (inputfilename.find("DYccJetsToLL_TuneZ2_M-50") < 200) {//DYccJetsToLL_TuneZ2_M-50
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Zcc");
	Vxsectionbkg.push_back(1.);
	Vcolorbkg.push_back(kRed+0);
      }

      // DYccJetsToLL_TuneZ2_M-10To50
      if (inputfilename.find("DYccJetsToLL_TuneZ2_M-10To50") < 200) {//DYccJetsToLL_TuneZ2_M-10To50
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Zcc");
	Vxsectionbkg.push_back(1.);
	Vcolorbkg.push_back(kRed+0);
      }

      // // DYToEE_M-10To20_TuneZ2_8TeV-pythia6
//       if (inputfilename.find("DYToEE")!=-1) {
// 	Vdatasetnamebkg.push_back(inputfilename);
// 	Vlabelbkg.push_back("DYToEE_M-10To20_TuneZ2_8TeV-pythia6");
// 	Vxsectionbkg.push_back(1.); //pb
// 	Vcolorbkg.push_back(kAzure+9);
//       }

//       // DYToMuMu_M-20_TuneZ2_8TeV-pythia6
//       if (inputfilename.find("DYToMuMu")!=-1) {
// 	Vdatasetnamebkg.push_back(inputfilename);
// 	Vlabelbkg.push_back("DYToMuMu_M-20_TuneZ2_8TeV-pythia");
// 	Vxsectionbkg.push_back(1.); //pb
// 	Vcolorbkg.push_back(kAzure-7);
//       }

      // DYToMuMu_M-20_TuneZ2_8TeV-pythia6
      if (inputfilename.find("DYtoMuMu") < 200) {
 	Vdatasetnamebkg.push_back(inputfilename);
 	Vlabelbkg.push_back("#gamma^{*} /Z#rightarrow #mu^{+} #mu^{-}");
 	Vxsectionbkg.push_back(1.); //pb
 	Vcolorbkg.push_back(kAzure-9);
      }

      // ZToMuMu_NNPDF30_13TeV-powheg_M
      if (inputfilename.find("ZToMuMu_NNPDF30_13TeV-powheg_M") < 200) {//ZToMuMu_NNPDF30_13TeV-powheg_M
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Z+jets");
	Vxsectionbkg.push_back(1.);
	Vcolorbkg.push_back(kAzure-9);
      }
      // ZToEE_NNPDF30_13TeV-powheg_M
      if (inputfilename.find("ZToEE_NNPDF30_13TeV-powheg_M") < 200) {//ZToEE_NNPDF30_13TeV-powheg_M
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("Z+jets");
	Vxsectionbkg.push_back(1.);
	Vcolorbkg.push_back(kAzure-9);
      }

      // QCD_Pt-30to40_doubleEMEnriched_TuneZ2_8TeV-pythia6
      if (inputfilename.find("30to40_doubleEMEnriched_TuneZ2_8TeV-pythia6") <52) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-30to40_doubleEMEnriched_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal+8);
      }
      // QCD_Pt-40_doubleEMEnriched_TuneZ2_8TeV-pythia6
      if (inputfilename.find("-40_doubleEMEnriched_TuneZ2_8TeV-pythia6") <52) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-40_doubleEMEnriched_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal+8);
      }


      // QCD_Pt-15to20_MuPt5Enriched_TuneZ2_8TeV-pythia6
      if (inputfilename.find("QCD_Pt-15to20_MuPt5Enriched_TuneZ2_8TeV-pythia6") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-15to20_MuPt5Enriched_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-8);
      }

      // QCD_Pt-20to30_MuPt5Enriched_TuneZ2_8TeV-pythia6
      if (inputfilename.find("QCD_Pt-20to30_MuPt5Enriched_TuneZ2_8TeV-pythia6") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-20to30_MuPt5Enriched_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-8);
      }

      // QCD_Pt-30to50_MuPt5Enriched_TuneZ2_8TeV-pythia6
      if (inputfilename.find("QCD_Pt-30to50_MuPt5Enriched_TuneZ2_8TeV-pythia6") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-30to50_MuPt5Enriched_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-8);
      }

      // QCD_Pt-50to80_MuPt5Enriched_TuneZ2_8TeV-pythia6
      if (inputfilename.find("QCD_Pt-50to80_MuPt5Enriched_TuneZ2_8TeV-pythia6") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-50to80_MuPt5Enriched_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-8);
      }

      // QCD_Pt-80to120_MuPt5Enriched_TuneZ2_8TeV-pythia6
      if (inputfilename.find("QCD_Pt-80to120_MuPt5Enriched_TuneZ2_8TeV-pythia6") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-80to120_MuPt5Enriched_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-8);
      }

      // QCD_Pt-120to150_MuPt5Enriched_TuneZ2_8TeV-pythia6
      if (inputfilename.find("QCD_Pt-120to150_MuPt5Enriched_TuneZ2_8TeV-pythia6") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-120to150_MuPt5Enriched_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-8);
      }

      // QCD_Pt-150_MuPt5Enriched_TuneZ2_8TeV-pythia6
      if (inputfilename.find("QCD_Pt-150_MuPt5Enriched_TuneZ2_8TeV-pythia6") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-150_MuPt5Enriched_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-8);
      }

      // QCD_Pt-20to30_BCtoE_TuneZ2_8TeV-pythia6
      if (inputfilename.find("QCD_Pt-20to30_BCtoE_TuneZ2_8TeV-pythia6") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-20to30_BCtoE_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      // QCD_Pt-30to80_BCtoE_TuneZ2_8TeV-pythia6
      if (inputfilename.find("QCD_Pt-30to80_BCtoE_TuneZ2_8TeV-pythia6") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt-30to80_BCtoE_TuneZ2_8TeV-pythia6");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }


      // NEW QCD Pythia 8
      if (inputfilename.find("QCD_Pt_1000to1400") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_10to15") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_120to170") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_1400to1800") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_15to30") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_170to300") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_1800to2400") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_2400to3200") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_300to470") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_30to50") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_3200toInf") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_470to600") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_50to80") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_5to10") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_600to800") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_800to1000") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }

      if (inputfilename.find("QCD_Pt_80to120") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-2);
      }


      // TTbar
      if (inputfilename.find("TTbar") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("t#bar{t} + jets");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-6);
      }
      if (inputfilename.find("TTJets") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("t#bar{t} + jets");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-6);
      }
      // TTTo2L2Nu2B_8TeV-powheg-pythia6
      if (inputfilename.find("TTTo2L2Nu") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("t#bar{t}");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kTeal-6);
      }

      // TTToLL_MLL
      if (inputfilename.find("TTToLL_MLL") < 200) {
        Vdatasetnamebkg.push_back(inputfilename);
        Vlabelbkg.push_back("t#bar{t}");
        Vxsectionbkg.push_back(1.); //pb
        Vcolorbkg.push_back(kTeal-6);
      }

      // Single Top
      // T_TuneZ2_s-channel_13TeV-madgraph
      if (inputfilename.find("ST_s-channel_top") < 200 ||
	  inputfilename.find("T_TuneZ2_s-channel") < 200 ||
	  inputfilename.find("singletop") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("T_TuneZ2_s-channel_8TeV-powheg-tauola");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      }

      // T_TuneZ2_t-channel_13TeV-madgraph
      if (inputfilename.find("ST_t-channel_top") < 200 ||
	  inputfilename.find("T_TuneZ2_t-channel") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      }

      // T_TuneZ2_tW-channel_8TeV-madgraph
      if (inputfilename.find("tW_top") < 200 ||
	  inputfilename.find("T_TuneZ2_tW-channel") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      }

      // Tbar_TuneZ2_s-channel_8TeV-madgraph
      if (inputfilename.find("ST_s-channel_antitop") < 200 ||
	  inputfilename.find("Tbar_TuneZ2_s-channel") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("ST_s-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      }

      // Tbar_TuneZ2_t-channel_8TeV-madgraph
      if (inputfilename.find("ST_t-channel_antitop") < 200 ||
	  inputfilename.find("Tbar_TuneZ2_t-channel") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      }

      // Tbar_TuneZ2_tW-channel-DR_8TeV-madgraph
      if (inputfilename.find("ST_tW_antitop") < 200 ||
	  inputfilename.find("tbarW") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kViolet);
      }

      // W+jets
      // WJetsToLNu_TuneZ2_8TeV-madgraph-tauola
      if (inputfilename.find("WJetsToLNu") < 200) {
	Vdatasetnamebkg.push_back(inputfilename);
	Vlabelbkg.push_back("W+Jets");
	Vxsectionbkg.push_back(1.); //pb
	Vcolorbkg.push_back(kSpring);
      }
    }
  }
  infile.close();
}

void PlotStackZprime::printnumbers(char name[400], TH1F *h)
{
  int nSkim_MC=h->GetBinContent(1);
  int nZLQpt_MC=h->GetBinContent(2);
  int nZLQptLI_MC=h->GetBinContent(3);
  int n4LQptLI_MC=h->GetBinContent(4);
  int n4lAfterPresel_MC=h->GetBinContent(5);
  int nIso_MC=h->GetBinContent(6);
  int nIsopt_MC=h->GetBinContent(7);
  int nIP_MC=h->GetBinContent(8);
  int nKin_MC=h->GetBinContent(9);

  std::cout << "Sample= "  <<  name << std::endl;
  std::cout << "Skim MC= " <<  nSkim_MC
	    << " \n ZLQpt_MC="         << nZLQpt_MC
	    << " \n ZLQptLI_MC="       << nZLQptLI_MC
	    << " \n 4LQptLI_MC="       << n4LQptLI_MC
	    << " \n BestCand_MC="      << n4lAfterPresel_MC
	    << " \n Iso_MC="           << nIso_MC
	    << " \n Iso(pt)_MC="       << nIsopt_MC
	    << " \n IP_MC="            << nIP_MC
	    << " \n Kin_MC="           << nKin_MC
	    << std::endl;
}

/*
void PlotStackZprime::createdatacards(float Higgsm, float channel, float energy, float masslow, float masshigh, float ggH, float qqH, float WH, float ZH, float ttH, float bkg_qqzz, float bkg_ggzz, float bkg_zjets)
{


  Char_t txtOUT[500];
  //sprintf(txtOUT,"datacards/%s/hzz4l_%s_txt.txt",Higgsm,datasetName.Data());
  std::cout << "Opening a datacard file " << txtOUT << std::endl;
  ofstream output_txt;
  output_txt.open(txtOUT);

  output_txt << "imax 1" << std::endl;
  output_txt << "jmax 7" << std::endl;
  output_txt << "kmax *" << std::endl;
  output_txt <<"------------" << std::endl;
  output_txt <<"shapes * * hzz4l_4muS_8TeV.input.root w:$PROCESS" << std::endl;
  output_txt <<"------------" << std::endl;
  output_txt <<"bin a1" << std::endl;
  output_txt <<"observation 16" << std::endl;
  output_txt <<"------------" << std::endl;
  output_txt <<"## mass window [105.0,140.0]" << std::endl;
  output_txt <<"bin a1 a1 a1 a1 a1 a1 a1 a1" << std::endl;
  output_txt <<"process ggH qqH WH ZH ttH bkg2d_qqzz bkg2d_ggzz bkg2d_zjets" << std::endl;
  output_txt <<"process -4 -3 -2 -1 0 1 2 3" << std::endl;
  output_txt <<"rate 1.0000 1.0000 1.0000 1.0000 1.0000 7.6204 0.1543 1.1796" << std::endl;
  output_txt <<"------------" << std::endl;
  output_txt <<"lumi_8TeV lnN 1.026 1.026 1.026 1.026 1.026 1.026 1.026 -" << std::endl;
  output_txt <<"pdf_gg lnN 1.0720 - - - 1.0780 - 1.0708 -" << std::endl;
  output_txt <<"pdf_qqbar lnN - 1.0270 1.0350 1.0350 - 1.0341 - -" << std::endl;
  output_txt <<"pdf_hzz4l_accept lnN 1.02 1.02 1.02 1.02 1.02 - - -" << std::endl;
  output_txt <<"QCDscale_ggH lnN 1.0750 - - - - - - -" << std::endl;
  output_txt <<"QCDscale_qqH lnN - 1.0020 - - - - - -" << std::endl;
  output_txt <<"QCDscale_VH lnN - - 1.0040 1.0155 - - - -" << std::endl;
  output_txt <<"QCDscale_ttH lnN - - - - 1.0655 - - -" << std::endl;
  output_txt <<"QCDscale_ggVV lnN - - - - - - 1.2431 -" << std::endl;
  output_txt <<"QCDscale_VV lnN - - - - - 1.0284 - -" << std::endl;
  output_txt <<"BRhiggs_hzz4l lnN 1.02 1.02 1.02 1.02 1.02 - - -" << std::endl;
  output_txt <<"CMS_eff_m lnN 1.043 1.043 1.043 1.043 1.043 1.043 1.043 -" << std::endl;
  output_txt <<"CMS_hzz4mu_Zjets lnN - - - - - - - 0.6/1.4" << std::endl;
  output_txt <<"CMS_zz4l_bkgMELA param 0  1  [-3,3]" << std::endl;
  output_txt <<"CMS_zz4l_mean_m_sig param 0.0 1.0" << std::endl;
  output_txt <<"## CMS_zz4l_mean_m_sig = 0.001" << std::endl;
  output_txt <<"CMS_zz4l_sigma_m_sig param 0.0 0.2" << std::endl;
  output_txt <<"CMS_zz4l_n_sig_1_8 param 0.0 0.01" << std::endl;
  output_txt <<"interf_ggH param 0 1 [-1,1]" << std::endl;

  output_txt.close();
}
*/

/*
void PlotStackZprime::getMassWindow(float Higgsm)
{

  std::string fileLoc = "/cmshome/nicola/tmp/test/Paper/last/last/CMSSW_5_3_9/src/Higgs/Higgs_CS_and_Width/txtFiles";
  HiggsCSandWidth myCSW = HiggsCSandWidth(fileLoc);

  double widthHVal =  myCSW.HiggsWidth(0,Higgsm);
  double windowVal = max( widthHVal, 1.0);
  double lowside = 100.0;
  double highside = 1000.0;

  if (Higgsm >= 275) {
    lowside = 180.0;
    highside = 650.0;
  }
  if (Higgsm >= 350) {
    lowside = 200.0;
    highside = 900.0;
  }
  if (Higgsm >= 500) {
     lowside = 250.0;
     highside = 1000.0;
  }
  if (Higgsm >= 700) {
    lowside = 350.0;
    highside = 1400.0;
  }

   float low_M = max( (Higgsm - 20.*windowVal), lowside);
   float high_M = min( (Higgsm + 15.*windowVal), highside);

   std::cout << "Mass window: " << low_M << " " << high_M << std::endl;
}
*/

//  LocalWords:  TH1F hfourlepbestmass
