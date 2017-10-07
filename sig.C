
using namespace RooFit;
using namespace RooStats;

void makeModel(RooWorkspace& w) {

   TFile *_file0 = TFile::Open("plots/htotal_root_ZprimeRecomass.root");
   TH1F *Histo = (TH1F*)_file0->Get("htotaldata");
   RooRealVar invm("invm","invm",200.,4000.);    
   RooDataHist* data = new RooDataHist("data","data",invm,Import(*Histo)) ;
   

//   TTree* tree = new TTree("simple","data from ascii file");
//   Long64_t nlines = tree->ReadFile("list_mll_200_2016.txt","x1:x2:x3:invm:x5:x6");
//   Long64_t nlines = tree->ReadFile("a.txt","x1:x2:x3:invm:x5:x6");
//   printf(" found %lld pointsn",nlines);
//   tree->Write();
//   tree->GetEntries();

   RooRealVar mass("mass","mass", 300., 200., 1600.);
   RooRealVar nsig("nsig","Number of signal events", 0., 5000.);
   RooRealVar nbkg("nbkg","Number of background events", 0., 300000.);
   w.import(mass);
   w.import(nsig);
   w.import(nbkg);

//   RooRealVar invm("invm","Invariant mass", 200., 4000.);
//   RooDataSet* data = new RooDataSet("data", "Data", invm, RooFit::Import(*tree));

   data->Print("v");
   w.import(invm);
   w.import(*data);
 
   w.factory("expr::sigma('invm*(0.01292 + 0.00001835 * invm - 0.0000000002733 * invm*invm)',invm)");
   w.factory("expr::width('0.03*invm',invm)");
 
   w.factory("CEXPR::bkgpdf('exp(24.9327 - 2.39287e-03*invm + 3.19926e-07*invm*invm - 3.38799e-11*invm*invm*invm)*pow(invm,-3.3634)',invm)");
   w.factory("Voigtian::sigpdf(invm,mass,width,sigma)");

   w.factory("SUM::model(nbkg*bkgpdf, nsig*sigpdf)");

   RooAbsPdf* sigpdf = w.pdf("sigpdf");
   RooAbsPdf* bkgpdf = w.pdf("bkgpdf");
   RooAbsPdf* model  = w.pdf("model");

   RooStats::ModelConfig* mc = new ModelConfig("mc",&w);
   mc->SetPdf(*w.pdf("model"));
   mc->SetParametersOfInterest(*w.var("nsig"));
   mc->SetObservables(*w.var("invm"));
   w.defineSet("nuisParams","nbkg");

   mc->SetNuisanceParameters(*w.set("nuisParams"));
   w.var("mass")->setConstant(true);

   w.import(*mc);

   w.Print("tree");

   w.writeToFile("MyModel_workspace.root");
  
   TCanvas* c1 = new TCanvas("c1","Control Plots", 900, 700);
   RooPlot* plot = w.var("invm")->frame();
   w.data("data")->plotOn(plot);
   w.pdf("model")->plotOn(plot);
   w.pdf("model")->plotOn(plot, Components("bkgpdf"),LineStyle(kDashed));
   w.pdf("model")->plotOn(plot, Components("sigpdf"),LineColor(kRed));
   plot->Draw();

   return;

}


void significance(RooWorkspace& w ) {

  ModelConfig* mc = (ModelConfig*)w.obj("mc");
  RooDataSet* data = (RooDataSet*)w.data("data");
  //data->Print();

  // define the S+B snapshot (this is used for computing the expected significance)
  ModelConfig* sbModel = mc->Clone();
  sbModel->SetName("S+B Model");
  RooRealVar* poi = (RooRealVar*) sbModel->GetParametersOfInterest()->first();
  poi->setVal(50);
  sbModel->SetSnapshot(*poi);

  ModelConfig * bModel = (ModelConfig*) sbModel->Clone();
  bModel->SetName("B model");
  poi->setVal(0);
  bModel->SetSnapshot(*poi);

  vector<double> masses;
  vector<double> p0values;
  vector<double> p0valuesExpected;
  vector<double> sigvalues;

  double massMin = 200;
  double massMax = 2500;
  int nbins = 100;

  // loop on the mass values 
  for ( double mass=massMin; mass<=massMax; mass += (massMax-massMin)/nbins ) {

    w.var("mass")->setVal( mass );

    // create the AsymptoticCalculator from data,alt model, null model
    AsymptoticCalculator * ac = new AsymptoticCalculator(*data, *sbModel, *bModel);
    ac->SetOneSidedDiscovery(true);  // for one-side discovery test                                      
    AsymptoticCalculator::SetPrintLevel(-1);

    // run the calculator
    HypoTestResult* asymCalcResult = ac->GetHypoTest();
    asymCalcResult->Print();
     
    double pvalue = asymCalcResult->NullPValue();
    double sigvalue = asymCalcResult->Significance();
    double expectedP0 = AsymptoticCalculator::GetExpectedPValues(asymCalcResult->NullPValue(),asymCalcResult->AlternatePValue(), 0, false);

    masses.push_back(mass);
    p0values.push_back(pvalue);
    p0valuesExpected.push_back(expectedP0);
    sigvalues.push_back(sigvalue);
    std::cout << "** Mass = " << mass << " p0-value = " << expectedP0 << " p-value = " << pvalue << " significance = " << sigvalue << std::endl;

  }

  TGraph* graph1  = new TGraph(masses.size(),&masses[0],&p0values[0]);
  TGraph* graph2  = new TGraph(masses.size(),&masses[0],&p0valuesExpected[0]);
  TGraph* graph3  = new TGraph(masses.size(),&masses[0],&sigvalues[0]);

  TCanvas* c2 = new TCanvas("c2","Significance", 900, 700);
  c2->Divide(1,2);
  c2->cd(1);
  graph1->SetMarkerStyle(10);
  //graph1->Draw("APC");
  graph1->Draw("AC");
  graph2->SetLineStyle(2);
  graph2->Draw("C");
  graph1->GetXaxis()->SetTitle("Mass [GeV]");
  graph1->GetYaxis()->SetTitle("p0 value");
  graph1->SetTitle("P-value vs Mass");
  graph1->SetMinimum(graph2->GetMinimum());
  graph1->SetLineColor(kBlue);
  graph2->SetLineColor(kRed);
  gPad->SetLogy(true);

  c2->cd(2);
  graph3->SetMarkerStyle(10);
  graph3->Draw("AC");
  graph3->SetLineStyle(1);
  graph3->SetLineColor(kRed);
  graph3->GetXaxis()->SetTitle("Mass [GeV]");
  graph3->GetYaxis()->SetTitle("Significance");
  graph3->SetTitle("Significance vs Mass");
  gPad->SetLogy(false);

  c2->SaveAs("significance.pdf");
  c2->SaveAs("significance.png");
}


void sig() {

  RooWorkspace w("w");

  // create model
  makeModel(w);

  // calculate significance
  significance(w);

  return;

}
