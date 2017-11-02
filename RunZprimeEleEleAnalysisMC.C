#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <string>
#include <memory>
#include <iostream>
#include <TSystem.h>
#include <TH2.h>
#include "TChain.h"
#include <stdlib.h>
#include <TDCacheFile.h>
#include "ZprimeEleElePatMiniAodNewMC.h"
#endif

int main(int argc, char ** argv)
{
  std::cout << "This is " << argv[1] << std::endl;
  std::string sampletype=argv[1];
  std::ifstream fdata;
  int nlines;
  std::string dataconf="";
  std::string mcconf="";

  if (sampletype.find("sig") < 10) {
    fdata.open(argv[2]);
    nlines = atoi(argv[3]);
    mcconf="Spring16";
  } else if (sampletype.find("bkg") < 10) {
    fdata.open(argv[4]);
    nlines = atoi(argv[5]);
    mcconf="Spring16";
  } else if (sampletype.find("data") < 10) {
    fdata.open(argv[6]);
    nlines = atoi(argv[7]);
    dataconf="2016";
  }

  std::string samples[nlines];
  float ninput[nlines];
  float nhlt[nlines];
  float nskim[nlines];
  float xsection[nlines];

  for (int i = 0; i < nlines; ++i) {
    fdata >> samples[i] >> ninput[i] >> nhlt[i] >> nskim[i] >> xsection[i];
    std::cout << "Sample=" << samples[i] << " Ninput=" << ninput[i]
              << " NHLT="  << nhlt[i]    << " NSkim="  << nskim[i]
              << " Xsection(pb)=" << xsection[i] << std::endl;
  }

  //
  float lumifb=0.;
  if (mcconf.find("Spring16") < 5)
    lumifb=36.3 ; // 2016

  std::string site=argv[8];
  std::cout << "Site is " << site.c_str() << " MC conf.= " << mcconf.c_str() << " data conf.= " << dataconf.c_str() << std::endl;

  // Run on data
  for (int i = 0; i < nlines; ++i) {
    std::string name;
    if (mcconf.find("Spring16") < 10)
      name = samples[i]+".root";
    if (dataconf.find("2016") < 10)
      name = samples[i]+".root";

    TString dirInput;
    if (site.find("CERN") < 5) {
      dirInput="/castor/cern.ch/user/n/ndefilip/Paper/MCFall11";    // to run at CERN
    } else if (site.find("DESY") < 5) {
      if (dataconf.find("2016") < 50) {
        dirInput="/nfs/dust2/cms/group/DAS2016/ZprimeDiLepton/Data2016_ZprimeEE_13TeV_merged_HLT"; //to run at DESY
      }
      if (mcconf.find("Spring16") < 50) {
        dirInput="/nfs/dust2/cms/group/DAS2016/ZprimeDiLepton//Spring16_ZprimeEE_13TeV_merged";
        if (name.find("reHLT_DYtoEE") < 100) dirInput="/nfs/dust2/cms/group/DAS2016/ZprimeDiLepton/Spring16_ZprimeEE_13TeV_merged_HLT";
      }
    } else if (site.find("FNAL") < 5 && mcconf.find("Spring15_combined") < 5) {
      dirInput="root://cmseos.fnal.gov///store/user/cmsdas/2016/LONG_EXERCISES/ZprimeDiLeptons/Spring15_25ns_merged";
    } else if (site.find("FNAL") < 5 && dataconf.find("2015") < 5) {
      dirInput="root://cmseos.fnal.gov///store/user/cmsdas/2016/LONG_EXERCISES/ZprimeDiLeptons/Data2015_ZprimeEE_13TeV_merged";
    } else if (mcconf.find("Spring16") < 50) {
      dirInput="root://cmseos.fnal.gov///store/group/lpcci2dileptons/ZprimeDiLeptonsAnalysis2017/MonteCarlo_Moriond";
      if (name.find("reHLT_DYtoEE") < 100) dirInput="/lustre/cms/store/user/defilip/ZprimeAnalysis/Spring16_ZprimeEE_13TeV_merged_HLT";
      if (name.find("CMSSW803_MC_DYtoEE") < 100) dirInput="/lustre/cms/store/user/defilip/ZprimeAnalysis/Spring16_merged";
      //  dirInput="/lustre/cms/store/user/aqamesh/ZprimeAnalysis_QCD_Tree";
      //dirInput="/lustre/cms/store/user/selgammal/ZprimeEE/MCs";
    } else if (dataconf.find("2016") < 50) {
      //dirInput="/lustre/cms/store/user/selgammal/ZprimeEE/2016Data/";
      //dirInput="/lustre/cms/store/user/defilip/ZprimeAnalysis/Data2016_ZprimeEE_13TeV_merged_HLT";
      dirInput="root://cmseos.fnal.gov///store/group/lpcci2dileptons/ZprimeDiLeptonsAnalysis2017/Keep_Moriond17_reMINIAOD_Data_ele_final";
    }

    TString File = name;
    Char_t namechar[300];
    sprintf(namechar,"%s/%s",dirInput.Data(),File.Data());
    float weight= -999.;
    if (mcconf.find("Spring16") < 50) {
      if (name.find("reHLT_DYtoEE") < 100 ) {
	weight=0.96*lumifb*(xsection[i]*1000.*nskim[i]/ninput[i])/nskim[i];
      } else if (name.find("ZToEE") < 50) {
	weight=0.9714*lumifb*(xsection[i]*1000.*nskim[i]/ninput[i])/nskim[i];
      } else {
	weight=lumifb*(xsection[i]*1000.*nskim[i]/ninput[i])/nskim[i];
      }
    } else if (dataconf.find("2016") < 50) {
      weight=1.;
    }

    std::cout << "weight is " << weight << std::endl;

    // TFile *file3 = TFile::Open(namechar);
    // std::shared_ptr<TFile> file3 = std::make_shared<TFile>((TFile*)TFile::Open(namechar));
    std::shared_ptr<TFile> file3(TFile::Open(namechar));
    std::cout << "Read file (" << std::hex << file3 << std::dec << ") with name: " << namechar << std::endl;

    TTree *tree3 = (TTree*)file3->Get("tree");
    // std::shared_ptr<TTree> tree3 = std::make_shared<TTree>((TTree*)file3->Get("tree"));
    // std::shared_ptr<TTree> tree3((TTree*)file3->Get("tree"));
    std::cout << "Read file with name: " << namechar << " " << tree3->GetEntries() << std::endl;

    ZprimeEleElePatMiniAodNewMC b(namechar,tree3,weight,dataconf,mcconf);
    // ZprimeEleElePatMiniAodNewMC b(namechar,tree3.get(),weight,dataconf,mcconf);
    b.Loop(false);
    // tree3 = nullptr;
    // delete tree3;  // didn't create with new, no need to delete
    file3->Close();
    // file3 = nullptr;
  }
  return 0;
}
