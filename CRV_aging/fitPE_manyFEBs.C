#include <string>
#include <fstream>
#include <TMath.h>
#include <stdio.h>

struct fitInfo {
  float maxvalx, gauspar, gauserr, integral, chi2, rchi2;
};

typedef struct fitInfo Struct;

int nFEB = 4;
float NC = 64 * nFEB;

Double_t lowerEdge = 0.;
Double_t upperEdge = 3000.;
Int_t nbins = 30;
TH1F* chi2hist = new TH1F("chi2hist","chi2hist",500,0.,500.);
TH1F* gauserrhist = new TH1F("gauserrhist","gauserrhist",500,0.,500.);
TH1F* rchi2hist = new TH1F("rchi2hist","rchi2hist",500,0.,50.);
TH1F* nentriesHist = new TH1F("nentriesHist","integral of histo",256,0.,256.);
TH1F* nPEsHist = new TH1F("nPEsHist","Corrected PE Yield per Channel",256,0.,256.);
TH2F* PE_run[256];
TH2F* PE_corrected[256];
TH2F* gaus_run[256];
TH2F* gaus_corrected[256];
const int npars = 6;
double par[npars];

Struct fitFunc(TH1F *h,  int ithHist){
  Struct fitPars;

  float mpv = h->GetBinCenter(h->GetMaximumBin());
  TF1 *g1 = new TF1("g1", "gaus", mpv-20, mpv+20);
  g1->SetLineColor(kBlue);
  //TF1 *g2 = new TF1("g2", "landau", mpv+20, 120);
  //TF1 *total = new TF1("total", "gaus(0)+landau(3)", mpv-20, mpv+60);
  //TF1 *total = new TF1("total", "gaus(0)+landau(3)+pol2(6)", mpv-20, mpv+60);
  //TF1 *total = new TF1("total", "gaus(0)+pol2(3)", mpv-20, mpv+60);
  TF1 *total = new TF1("total", "gaus(0)+pol2(3)", mpv-30, mpv+30);

  h->Fit(g1, "R+");
  //h->Fit(g2, "R+");

  g1->GetParameters(&par[0]);
  //g2->GetParameters(&par[3]);

  for(int i=0; i<3; i++){
    total->SetParLimits(i, 0., 1e4);
  }
  for(int i=3; i<npars; i++){
    total->SetParLimits(i, -1.e4, 1e4);
  }

 total->SetParameters(par);

  //h->Smooth(1);
  //  h->Fit(total, "R");
  Int_t fitStatus = h->Fit(total,"R+");
  std::cout << "ithHist, number of entries = " << ithHist << " " << h->Integral() << std::endl;
 std:cout << "ithHist convergence = " << ithHist << " " << fitStatus << std::endl;
  Double_t number = h->Integral();
  if (number < 1.) {number = 0.01;}
  std::cout << "integral = " << number << std::endl;
  nentriesHist->Fill(static_cast<float>(ithHist)+0.1,number);
  //need the +1 see https://root.cern.ch/doc/master/classTH1.html on bin def'ns.  Bin 0 is underflows, Bin 1 is first bin
  nentriesHist->SetBinError(ithHist+1,0.01);
  total->GetParameters(&par[0]);
  if (number < 2000.){
    fitPars.gauspar = 0.;
    fitPars.gauserr = 0.;
    fitPars.integral = number;
    fitPars.maxvalx = 0.;
    fitPars.chi2 = 0.;
    fitPars.rchi2 = 0.;
  }
  if (number > 2000.){
    float p1 = total->GetParameter(1);
    float p1_err = total->GetParError(1);
    fitPars.gauspar = p1;
    fitPars.gauserr = p1_err;
    fitPars.integral = number;
    gauserrhist->Fill(p1_err);

    // extract chi2
    TF1* fit = h->GetFunction("total");
    Float_t chi2 = fit->GetChisquare();
    Int_t   ndf  = fit->GetNDF();
    Float_t rchi = chi2/ndf;
    fitPars.chi2 = chi2;
    fitPars.rchi2 = rchi;
    if (chi2 > 200){std::cout << "ithhist, chi2 = " << ithHist << " "  << chi2 << std::endl;}
    chi2hist->Fill(chi2);
    rchi2hist->Fill(rchi);
    // now total is the final fit, get peak value of total
    total->Draw();
    g1->Draw("sames");
    int maxval  = total->GetMaximum(20,100);
    float maxvalx = total->GetMaximumX(20,100);
    fitPars.maxvalx = maxvalx;
    std::cout << "ithhist = " << ithHist << " " << "maxval = " << maxval << ", maxvalx = " << maxvalx << std::endl;
    ++ithHist;
  }
  return fitPars;
}

float tempCorrection(float PE, float temp){
  float PE_corr;
  float ref_temp = 25.0;

  PE_corr = PE + 0.7506 * (temp - ref_temp);
  
  return PE_corr;
}

void fitPE_manyFEBs(){

  ifstream input("tryMe_0.txt");
  string oneline;
  //TFile *file = TFile::Open("/pnfs/mu2e/tape/phy-rec/rec/mu2e/CRV_wideband_cosmics/crvaging-001/root/e3/a2/rec.mu2e.CRV_wideband_cosmics.crvaging-001.001009_000.root");
  //TFile *file = TFile::Open("/pnfs/mu2e/tape/phy-rec/rec/mu2e/CRV_wideband_cosmics/crvaging-001/root/91/5f/rec.mu2e.CRV_wideband_cosmics.crvaging-001.000119_000.root");

  while(std::getline(input, oneline)){
    const char* cfile = oneline.c_str();
    TFile *file = TFile::Open(cfile);
    std::cout << "file " << cfile << " is opening..." << std::endl;

    gStyle->SetOptStat(1111111);

    TIter keyList(file->GetListOfKeys());
    TKey* key;
    TString key_name;
    while((key = (TKey*)keyList())){
      key_name = key->GetName();
      std::cout << key_name << std::endl;
    }
    TTree *tree = file->Get<TTree>(key_name);
    TString filename;
    filename = file->GetName();
    std::cout << "filename: " << filename << " length: " << filename.Length() << std::endl;
    int run_num;

    TString short_runstr = filename(filename.Length()-15, filename.Length()-9);
    run_num = short_runstr.Atoi();
    int subrun_num;
    TString short_subrunstr = filename(filename.Length()-8, filename.Length()-5);
    subrun_num = short_subrunstr.Atoi();
    std::cout << "run number: " << run_num << ", subrun number: " << subrun_num << std::endl;

    //Create output file. If it already exists, recreate it.
    TFile *ofile = new TFile(Form("correctedPlots_posgaus_run%i_%i.root", run_num, subrun_num),"RECREATE");

    TCanvas *page[512];
    //page->Divide(4,4);
    TH1F *hPEs[256];
    TH1F* htemp[256];
    TH1F* hG_PEerr[256];
    TH1F* hG_PEcor_err[256];

    // Create a tree with branches (variables)
    TTree *newtree = new TTree("newtree","Data holder tree");
    Float_t max_xval;
    Float_t mean_temp;
    newtree->Branch("max_xval",&max_xval,"max_xval/F");
    newtree->Branch("mean_temp",&mean_temp,"mean_temp/F");

    string name;
    const char* cname;
    string tname;
    const char* tcname;
    Struct thisFit;
    float this_xval = -1.;
    float this_gaus = -1.;
    float this_gauserr = -1.;
    float this_chi2 = -1.;
    float this_rchi2 = -1.;
    float this_int  = -1.;
    float corr_xval = -1.;
    float corr_gaus = -1.;
    int ithHist = 0;
    for (int iFEB = 0; iFEB < 4; iFEB++){
      for(int i=0; i<64; i++){
	hPEs[i] = new TH1F(Form("hist_%d",(i + iFEB*64)),Form("hist_%d",(i + iFEB*64)), 100, 10, 120);
	name = "can" + to_string(i + iFEB*64);
	cname = name.c_str();
	tree->Draw(Form("PEs[%d][%d]>>%s",iFEB, i, hPEs[i]->GetName()) , "", "", 1000000, 0 ); // plot only 1M events
	hPEs[i]->Draw();
	hPEs[i]->GetXaxis()->SetTitle("Number of Photoelectrons");
	hPEs[i]->GetYaxis()->SetTitle("Events"); 
	float PE_mean = hPEs[i]->GetMean();
	std::cout << "PE mean for channel " << i << " in FEB " << iFEB << " is " << PE_mean << std::endl;
	if(PE_mean > 0.) {
	  thisFit = fitFunc(hPEs[i],ithHist);
	  this_xval = thisFit.maxvalx;
	  this_gaus = thisFit.gauspar;
	  this_gauserr = thisFit.gauserr;
	  this_int  = thisFit.integral;
	  this_chi2 = thisFit.chi2;
	  this_rchi2 = thisFit.rchi2;
	}
	for (int ithPar = 0; ithPar < npars; ++ithPar)
	  {
	    std::cout << "par " << i << " " << par[ithPar] << std::endl;
	  }
	if(PE_mean <= 0.) {
	  this_xval = 0.;
	  this_gaus = 0.;
	  this_gauserr = 0.;
	  this_chi2 = 0.;
	  this_rchi2 = 0.;
	  this_int  = thisFit.integral;
	  std::cout << "PE_mean < 0 " << ithHist << " " <<std::endl; 
	  nentriesHist->Fill(static_cast<float>(ithHist)+0.1,1.);
	  nentriesHist->SetBinError(ithHist+1,0.1);
	}
	hPEs[i]->Write();
	++ithHist;
	max_xval = this_xval;
	//std::cout << "max value from return = "<< max_xval << std::endl;

	htemp[i] = new TH1F(Form("temp_%d",i+(iFEB*64)),Form("temp_%d",i+(iFEB*64)), 100, 10, 120);
	tname = "tcan" + to_string(i + iFEB*64);
	tcname = tname.c_str();
	tree->Draw(Form("temperature[%d][%d]>>%s",iFEB, i, htemp[i]->GetName()) , "", "", 1000000, 0 ); // plot only 1M events
	htemp[i]->Draw();
	float temp = htemp[i]->GetMean();
	mean_temp = temp;
	htemp[i]->Write();

	//if(this_xval != 0.){
	  ofile->cd();
	  newtree->Fill();
	  PE_run[i] = new TH2F(Form("PE_vs_run_%d",(i + iFEB*64)), Form("PE_vs_run_%d",(i + iFEB*64)), 1100, 0., 1100., 200, 30., 50.);
	  PE_run[i]->SetMarkerStyle(20);
	  PE_run[i]->GetXaxis()->SetTitle("Run Number");
	  PE_run[i]->GetYaxis()->SetTitle("PE Yield (peak value from fit)");      
	  PE_run[i]->Fill(run_num, this_xval);
	  PE_run[i]->Write();

	  gaus_run[i] = new TH2F(Form("gaus_vs_run_%d",(i + iFEB*64)), Form("gaus_vs_run_%d",(i + iFEB*64)), 1100, 0., 1100., 200, 25., 50.);
	  gaus_run[i]->SetMarkerStyle(20);
	  gaus_run[i]->GetXaxis()->SetTitle("Run Number");
	  gaus_run[i]->GetYaxis()->SetTitle("PE Yield (Gaussian mean from fit)");   
	  gaus_run[i]->Fill(run_num, this_gaus);
	  gaus_run[i]->Write();

	  hG_PEerr[i] = new TH1F(Form("gaus_vs_run_err_%d",(i + iFEB*64)), Form("gaus_vs_run_err_%d",(i + iFEB*64)), 1100, 0., 1100.);
	  hG_PEerr[i]->SetMarkerStyle(20);
	  hG_PEerr[i]->GetXaxis()->SetTitle("Run Number");
	  hG_PEerr[i]->GetYaxis()->SetTitle("PE Yield (Gaussian mean from fit) with chi2 as error");   
	  hG_PEerr[i]->Sumw2();
	  hG_PEerr[i]->SetBinContent(run_num, this_gaus);
	  hG_PEerr[i]->SetBinError(run_num, this_chi2);
	  hG_PEerr[i]->Draw("E1");
	  hG_PEerr[i]->Write();

	  if(PE_mean > 0.) {
	    corr_xval = tempCorrection(this_xval, temp);
	    corr_gaus = tempCorrection(this_gaus, temp);
	  }
	  if(PE_mean <= 0.) {
	    corr_xval = 0.;
	    corr_gaus = 0.;
	  }

	  if((corr_gaus < 35.) || (corr_gaus > 65.)){
	    std::cout << "extreme PE channel: channel " << i << ", corrected Gaus PE of " << corr_gaus << " with temp " << temp << " and input " << this_gaus << std::endl;
	  }

	  PE_corrected[i] = new TH2F(Form("PE_corrected_%d",(i + iFEB*64)), Form("PE_corrected_%d",(i + iFEB*64)), 1100, 0., 1100., 200, 30., 50.);
	  PE_corrected[i]->SetMarkerStyle(20);
	  PE_corrected[i]->GetXaxis()->SetTitle("Run Number");
	  PE_corrected[i]->GetYaxis()->SetTitle("Corrected PE Yield (peak value from fit)");      
	  PE_corrected[i]->Fill(run_num, corr_xval);
	  PE_corrected[i]->Write();

	  gaus_corrected[i] = new TH2F(Form("gaus_corrected_%d",(i + iFEB*64)), Form("gaus_corrected_%d",(i + iFEB*64)), 1100, 0., 1100., 200, 25., 50.);
	  gaus_corrected[i]->SetMarkerStyle(20);
	  gaus_corrected[i]->GetXaxis()->SetTitle("Run Number");
	  gaus_corrected[i]->GetYaxis()->SetTitle("Corrected PE Yield (Gaussian mean from fit)");      
	  gaus_corrected[i]->Fill(run_num, corr_gaus);
	  gaus_corrected[i]->Write();

	  hG_PEcor_err[i] = new TH1F(Form("gaus_vs_run_cor_err_%d",(i + iFEB*64)), Form("gaus_vs_run_cor_err_%d",(i + iFEB*64)), 1100, 0., 1100.);
	  hG_PEcor_err[i]->SetMarkerStyle(20);
	  hG_PEcor_err[i]->GetXaxis()->SetTitle("Run Number");
	  hG_PEcor_err[i]->GetYaxis()->SetTitle("Corrected PE Yield (Gaussian mean from fit)");   
	  hG_PEcor_err[i]->Sumw2();
	  hG_PEcor_err[i]->SetBinContent(run_num, corr_gaus);
	  hG_PEcor_err[i]->SetBinError(run_num, this_gauserr);
	  hG_PEcor_err[i]->Draw("E1");
	  hG_PEcor_err[i]->Write();

	  nPEsHist->SetBinContent(i + iFEB*64 + 1, corr_gaus);
	  //nPEsHist->SetBinError(i + iFEB*64 + 1, this_gauserr);
	  //}
      }
    }
  
    nentriesHist->Write();
    nPEsHist->Write();
    chi2hist->Write();
    rchi2hist->Write();
    gauserrhist->Write();
    newtree->Write();

    nentriesHist->Reset("ICESM");
    nPEsHist->Reset("ICESM");
    chi2hist->Reset("ICESM");
    rchi2hist->Reset("ICESM");
    gauserrhist->Reset("ICESM");
  }
}
