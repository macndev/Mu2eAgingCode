#include <string>
#include <fstream>
#include <TMath.h>
#include <stdio.h>

struct fitInfo {
  float maxvalx, gauspar, gauserr, integral, chi2, rchi2;
};

typedef struct fitInfo Struct;

const int cnFEB = 6;
const int cnChan = 64;
const int cnLayer = 32;

Double_t lowerEdge = 0.;
Double_t upperEdge = 3000.;
Int_t nbins = 30;
TH1F* chi2hist = new TH1F("chi2hist","chi2hist",500,0.,500.);
TH1F* gauserrhist = new TH1F("gauserrhist","gauserrhist",500,0.,500.);
TH1F* rchi2hist = new TH1F("rchi2hist","rchi2hist",500,0.,50.);
TH1F* nentriesHist = new TH1F("nentriesHist","integral of histo",384,0.,384.);
TH1F* nPEsHist = new TH1F("nPEsHist","Corrected PE Yield per Channel",384,0.,384.);
TH2F* PE_run[cnFEB][cnChan];
TH2F* gaus_run[cnFEB][cnChan];
TH2F* CPE_run[cnFEB][cnChan];
TH2F* Cgaus_run[cnFEB][cnChan];
const int npars = 6;
double par[npars];

Struct fitFunc(TH1F *h,  int ithHist){
  Struct fitPars;

  float mpv = h->GetBinCenter(h->GetMaximumBin());
  TF1 *g1 = new TF1("g1", "gaus", mpv-20, mpv+20);
  g1->SetLineColor(kBlue);
  //TF1 *total = new TF1("total", "gaus(0)+landau(3)", mpv-20, mpv+60);
  //TF1 *total = new TF1("total", "gaus(0)+landau(3)+pol2(6)", mpv-20, mpv+60);
  //TF1 *total = new TF1("total", "gaus(0)+pol2(3)", mpv-20, mpv+60);
  TF1 *total = new TF1("total", "gaus(0)+pol2(3)", mpv-30, mpv+30);

  h->Fit(g1, "R+");

  g1->GetParameters(&par[0]);

  for(int i=0; i<3; i++){
    total->SetParLimits(i, 0., 1e4);
  }
  for(int i=3; i<npars; i++){
    total->SetParLimits(i, -1.e4, 1e4);
  }

 total->SetParameters(par);

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
    //g1->Draw("sames");
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

void aging_fitPE_v4(){

  ifstream input("module168_files_v3_all.txt");
  string oneline;
  //TFile *file = TFile::Open("/pnfs/mu2e/tape/phy-rec/rec/mu2e/CRV_wideband_cosmics/crvaging-001/root/e3/a2/rec.mu2e.CRV_wideband_cosmics.crvaging-001.001009_000.root");

  std::map < std::string, std::vector<Int_t> > badbias_map;

  int line_count = 0;
  int run_num;
  int subrun_num;
  int run_pos;
  int und_pos;
  int run_len;
  int sr_len;
  std::string run_s;
  std::string subrun_s;
  std::string run_subrun;
  int badspill;
  std::string badspill_s;
  std::vector<Int_t> badspills_thisrun; // initialize vector of ints for a particular sr

  // get bad spills from looking at bias voltages
  ifstream bad_spills("bad_bias_spills_168_v3.txt");
  string biasline;
  while(std::getline(bad_spills, biasline)){
    const char* linetext = biasline.c_str();
    std::string slinetext = linetext;
    int length = slinetext.size();

    if(length > 10){ // this is a header line for the sr to bad spill lines, get run_sr
      if(line_count != 0) {
	std::cout << "insert these: run/sr" << run_subrun << " bad spills ";
	for(Int_t i : badspills_thisrun){
	  std::cout << i;
	}
	std::cout << "\n";
	badbias_map.insert({run_subrun, badspills_thisrun}); // insert entry for each header before next one
	std::cout << "in the insert loop" << std::endl;
      }
      run_pos = slinetext.find("run");
      und_pos = slinetext.find("_",run_pos);
      run_len = und_pos - (run_pos+4);
      sr_len = (length - 1) - und_pos;
      run_s = slinetext.substr(run_pos+4,run_len);
      subrun_s = slinetext.substr(und_pos+1, sr_len-1);
      run_num = std::stoi(run_s);
      subrun_num = std::stoi(subrun_s);
      run_subrun = std::to_string(run_num) + "_" + std::to_string(subrun_num);
      std::cout << "run_subrun string: " << run_subrun << std::endl;
      badspills_thisrun.clear(); // make empty vector of ints for bad spills after each header
      }

    if(length <= 10){
      badspill_s = slinetext;
      badspill = std::stoi(badspill_s);
      badspills_thisrun.push_back(badspill);
      std::cout << "bad entry number: " << badspill << std::endl;
    }
    line_count++;
    std::cout << line_count << std::endl;
  }
  badbias_map.insert({run_subrun, badspills_thisrun});

  // print bad bias map to check it out
  std::cout << "bad bias map: " << std::endl;
  for(const auto& pair : badbias_map){
    std::cout << "{" << pair.first << ": ";
    for(Int_t f : pair.second){
      std::cout << f << ", ";
    }
    std::cout << "}\n";
  }

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
    TTree *tree = file->Get<TTree>("run");
    TString filename;
    filename = file->GetName();
    std::cout << "filename: " << filename << " length: " << filename.Length() << std::endl;
    int frun_num;
    std::string filerun_subrun;

    int nEvents;
    int nMuon;
    int nTriggered3;
    int nTriggered4;
    int nLayers_with_PEs;
    float nPE_per_layer_cutoff = 25.0;
    nEvents = tree->GetEntries();
    std::cout << "number of events in this file: " << nEvents << std::endl;

    TString short_runstr = filename(filename.Length()-15, filename.Length()-9);
    frun_num = short_runstr.Atoi();
    int fsubrun_num;
    TString short_subrunstr = filename(filename.Length()-8, filename.Length()-5);
    fsubrun_num = short_subrunstr.Atoi();
    filerun_subrun = std::to_string(frun_num) + "_" + std::to_string(fsubrun_num);
    std::cout << "run number: " << frun_num << ", subrun number: " << fsubrun_num << " , r_s: " << filerun_subrun << std::endl;

    //Create output file. If it already exists, recreate it.
    TFile *ofile = new TFile(Form("aging_fitPE_v3_mod168_run%i_%i.root", frun_num, fsubrun_num),"RECREATE");

    //std::cout << "const int constants: cnFEB = " << cnFEB << ", cnChan = " << cnChan << std::endl;
    //Declare histogram objects
    TH1F *hPEs[cnFEB][cnChan];
    TH1F *hCPEs[cnFEB][cnChan];
    TH1F *htemp[cnFEB][cnChan];
    TH1F *hG_PEcor_err[cnFEB][cnChan];
    TH1F *hG_PEun_err[cnFEB][cnChan];
    TH1F *hPEchi[cnFEB][cnChan];
    TH1F *hCPEchi[cnFEB][cnChan];
    TH1F *hPEfullchi[cnFEB][cnChan];

    // Create a tree with branches (variables)
    TTree *newtree = new TTree("newtree","Data holder tree");
    Float_t max_xval;
    Float_t mean_temp;
    newtree->Branch("max_xval",&max_xval,"max_xval/F");
    newtree->Branch("mean_temp",&mean_temp,"mean_temp/F");

    Struct thisFit;
    float PEs[cnFEB][cnChan], temperature[cnFEB][cnChan];
    float this_xval = -1.;
    float this_gaus = -1.;
    float this_gauserr = -1.;
    float this_chi2 = -1.;
    float this_rchi2 = -1.;
    float this_int  = -1.;
    float corr_xval = -1.;
    float corr_gaus = -1.;
    float nPE_in_layer = 0;
    float uPE_val = -1.;
    float cPE_val = -1.;

    int ithHist = 0;

    // declare histograms, do this only once
    for (int iFEB = 0; iFEB < cnFEB; iFEB++){
      for (int iChan = 0; iChan < cnChan; iChan++){
	hPEs[iFEB][iChan] = new TH1F(Form("hist_%d",(iChan + iFEB*cnChan)),Form("hist_%d; Number of Photoelectrons; Events",(iChan + iFEB*cnChan)), 100, 10, 120);
	hCPEs[iFEB][iChan] = new TH1F(Form("Chist_%d",(iChan + iFEB*cnChan)),Form("Chist_%d; Number of Photoelectrons; Events",(iChan + iFEB*cnChan)), 100, 10, 120);
	htemp[iFEB][iChan] = new TH1F(Form("temp_%d",iChan+(iFEB*cnChan)),Form("temp_%d; Temperature; Events",iChan+(iFEB*cnChan)), 100, 10, 120);
	PE_run[iFEB][iChan] = new TH2F(Form("PE_vs_run_%d",(iChan + iFEB*cnChan)), Form("PE_vs_run_%d; Run Number; Uncorrected PE Yield, max X value",(iChan + iFEB*cnChan)), 200, 1000., 1200., 200, 30., 70.);
	gaus_run[iFEB][iChan] = new TH2F(Form("gaus_vs_run_%d",(iChan + iFEB*cnChan)), Form("gaus_vs_run_%d; Run Number; Uncorrected PE Yield, Gaussian mean",(iChan + iFEB*cnChan)), 200, 1000., 1200., 200, 30., 70.);
	hG_PEun_err[iFEB][iChan] = new TH1F(Form("gaus_vs_run_un_err_%d",(iChan + iFEB*cnChan)), Form("gaus_vs_run_un_err_%d; Run Number; Uncorrected PE Yield, Gaussian mean",(iChan + iFEB*cnChan)), 200, 1000., 1200.);
	CPE_run[iFEB][iChan] = new TH2F(Form("CPE_vs_run_%d",(iChan + iFEB*cnChan)), Form("CPE_vs_run_%d; Run Number; Corrected PE Yield, max X value",(iChan + iFEB*cnChan)), 200, 1000., 1200., 200, 30., 70.);
	Cgaus_run[iFEB][iChan] = new TH2F(Form("Cgaus_vs_run_%d",(iChan + iFEB*cnChan)), Form("Cgaus_vs_run_%d; Run Number; Corrected PE Yield, Gaussian mean",(iChan + iFEB*cnChan)), 200, 1000., 1200., 200, 30., 70.);
	hG_PEcor_err[iFEB][iChan] = new TH1F(Form("gaus_vs_run_cor_err_%d",(iChan + iFEB*cnChan)), Form("gaus_vs_run_cor_err_%d; Run Number; Corrected PE Yield, Gaussian mean",(iChan + iFEB*cnChan)), 200, 1000., 1200.);
	hPEchi[iFEB][iChan] = new TH1F(Form("PEchi_%d",(iChan + iFEB*cnChan)), Form("PEchi_%d; Run Number; Reduced Chi2 of Uncorrected Fit",(iChan + iFEB*cnChan)), 100, 0., 10.);
	hCPEchi[iFEB][iChan] = new TH1F(Form("CPEchi_%d",(iChan + iFEB*cnChan)), Form("CPEchi_%d; Run Number; Reduced Chi2 of Corrected Fit",(iChan + iFEB*cnChan)), 100, 0., 10.);
	hPEfullchi[iFEB][iChan] = new TH1F(Form("PEfullchi_%d",(iChan + iFEB*cnChan)), Form("PEfullchi_%d; Run Number; Chi2 of Uncorrected Fit",(iChan + iFEB*cnChan)), 500, 0., 500.);
      }
    }

    // get entries to skip from badbias_map
    std::vector<Int_t> badentries_toskip;
    /*std::string this_string;
    for(const auto& pair : map){
      this_string = pair.first;
      if(this_string == filerun_subrun){
	badentries_toskip = pair.second;
      }
    }*/

    auto this_mapentry = badbias_map.find(filerun_subrun);
    badentries_toskip = this_mapentry->second;

    tree->SetBranchAddress("PEs", &PEs);
    tree->SetBranchAddress("temperature", &temperature);
    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
      tree->GetEntry(iEvent);
      bool badEntry = false;
      if(std::find(badentries_toskip.begin(), badentries_toskip.end(), iEvent) != badentries_toskip.end()){
	badEntry = true; //skip bad entries
	std::cout << "skipping entry " << iEvent << " from file " << filerun_subrun << std::endl;
      }
      nLayers_with_PEs = 0;
      int nFlag = 1;
      for (int iFEB = 0; iFEB < 1 && !badEntry; iFEB++) {
	for (int iLayer = 0; iLayer < 2; iLayer++){
	  nPE_in_layer = 0;
	  for (int iChan = 0; iChan < cnLayer; iChan++) {
	    uPE_val = PEs[iFEB][iChan+iLayer*cnLayer];

	    //get temperature by event (by spill)
	    mean_temp = temperature[iFEB][iChan+iLayer*cnLayer];
	    cPE_val = tempCorrection(uPE_val, mean_temp);
	    nPE_in_layer += cPE_val;

	    //std::cout << "we got here: temp = " << mean_temp << " PE vals, un = " << uPE_val << ", cor = " << cPE_val << std::endl;
	    //std::cout << "with indices iChan = " << iChan << ", iLayer = " << iLayer << ", cnLayer = " << cnLayer << std::endl;

	    hPEs[iFEB][iChan+iLayer*cnLayer]->Fill(uPE_val);
	    hCPEs[iFEB][iChan+iLayer*cnLayer]->Fill(cPE_val);
	    htemp[iFEB][iChan+iLayer*cnLayer]->Fill(mean_temp);
	  }
	  if (nPE_in_layer <= nPE_per_layer_cutoff){
	    nFlag *= 0;
	  }
	}
      }
      if (nFlag == 1){
	// do stuff
	nMuon += 1;
	for (int iFEB = 1; iFEB < cnFEB && !badEntry; iFEB++){
	  for (int iLayer = 0; iLayer < 2; iLayer++){
	    nPE_in_layer = 0;
	    for (int iChan = 0; iChan < cnLayer; iChan++) {
	      uPE_val = PEs[iFEB][iChan+iLayer*cnLayer];

	      mean_temp = temperature[iFEB][iChan+iLayer*cnLayer];
	      cPE_val = tempCorrection(uPE_val, mean_temp);
	      nPE_in_layer += cPE_val;

	      hPEs[iFEB][iChan+iLayer*cnLayer]->Fill(uPE_val);
	      hCPEs[iFEB][iChan+iLayer*cnLayer]->Fill(cPE_val);
	      htemp[iFEB][iChan+iLayer*cnLayer]->Fill(mean_temp);
	    }
	    if (nPE_in_layer <= nPE_per_layer_cutoff){
	      nLayers_with_PEs += 1;
	    }
	  }
	}
	if (nLayers_with_PEs == 3){
	  nTriggered3 += 1;
	}
	if (nLayers_with_PEs == 4){
	  nTriggered4 += 1;
	}
      }
      else {
	for (int iFEB = 1; iFEB < cnFEB && !badEntry; iFEB++){
	  for (int iLayer = 0; iLayer < 2; iLayer++){
	    nPE_in_layer = 0;
	    for (int iChan = 0; iChan < cnLayer; iChan++) {
	      uPE_val = PEs[iFEB][iChan+iLayer*cnLayer];
	      mean_temp = temperature[iFEB][iChan+iLayer*cnLayer];
	      cPE_val = tempCorrection(uPE_val, mean_temp);
	      hPEs[iFEB][iChan+iLayer*cnLayer]->Fill(uPE_val);
	      hCPEs[iFEB][iChan+iLayer*cnLayer]->Fill(cPE_val);
	      htemp[iFEB][iChan+iLayer*cnLayer]->Fill(mean_temp);
	    }
	  }
	}
      }
    }

    // draw histograms, do this only once
    for (int iFEB = 0; iFEB < cnFEB; iFEB++){
      for (int iChan = 0; iChan < cnChan; iChan++){
	hPEs[iFEB][iChan]->Draw();
	hCPEs[iFEB][iChan]->Draw();
	float PE_mean = hCPEs[iFEB][iChan]->GetMean();
	if(PE_mean > 0.) {
	  thisFit = fitFunc(hCPEs[iFEB][iChan],ithHist); // corrected fit
	  this_xval = thisFit.maxvalx;
	  this_gaus = thisFit.gauspar;
	  this_gauserr = thisFit.gauserr;
	  this_int  = thisFit.integral;
	  this_chi2 = thisFit.chi2;
	  this_rchi2 = thisFit.rchi2;
	  nPEsHist->SetBinContent(ithHist + 1, this_gaus);
	  
	  CPE_run[iFEB][iChan]->SetMarkerStyle(20);
	  CPE_run[iFEB][iChan]->Fill(frun_num, this_xval);
	  Cgaus_run[iFEB][iChan]->SetMarkerStyle(20);
	  Cgaus_run[iFEB][iChan]->Fill(frun_num, this_gaus);
	  hG_PEcor_err[iFEB][iChan]->SetMarkerStyle(20);
	  hG_PEcor_err[iFEB][iChan]->Sumw2();
	  hG_PEcor_err[iFEB][iChan]->SetBinContent(frun_num-999, this_gaus);
	  hG_PEcor_err[iFEB][iChan]->SetBinError(frun_num-999, this_gauserr);
	  hCPEchi[iFEB][iChan]->Fill(this_rchi2);
	  //std::cout << "filling cor error hist with: " << frun_num << ", " << this_gaus << ", " << this_gauserr << std::endl;
	  
	  thisFit = fitFunc(hPEs[iFEB][iChan],ithHist); // uncorrected fit
	  this_xval = thisFit.maxvalx;
	  this_gaus = thisFit.gauspar;
	  this_gauserr = thisFit.gauserr;
	  this_int  = thisFit.integral;
	  this_chi2 = thisFit.chi2;
	  this_rchi2 = thisFit.rchi2;
	  
	  PE_run[iFEB][iChan]->SetMarkerStyle(20);
	  PE_run[iFEB][iChan]->Fill(frun_num, this_xval);
	  gaus_run[iFEB][iChan]->SetMarkerStyle(20);
	  gaus_run[iFEB][iChan]->Fill(frun_num, this_gaus);
	  hG_PEun_err[iFEB][iChan]->SetMarkerStyle(20);
	  hG_PEun_err[iFEB][iChan]->Sumw2();
	  hG_PEun_err[iFEB][iChan]->SetBinContent(frun_num-999, this_gaus);
	  hG_PEun_err[iFEB][iChan]->SetBinError(frun_num-999, this_gauserr);
	  hPEchi[iFEB][iChan]->Fill(this_rchi2);
	  hPEfullchi[iFEB][iChan]->Fill(this_chi2);
	  //std::cout << "filling un error hist with: " << frun_num << ", " << this_gaus << ", " << this_gauserr << std::endl;
	  
	  ++ithHist;
	}
	htemp[iFEB][iChan]->Draw();
	hPEs[iFEB][iChan]->Write();
	hCPEs[iFEB][iChan]->Write();
	PE_run[iFEB][iChan]->Write();
	gaus_run[iFEB][iChan]->Write();
	hG_PEun_err[iFEB][iChan]->Draw("E1");
	hG_PEun_err[iFEB][iChan]->Write();
	CPE_run[iFEB][iChan]->Write();
	Cgaus_run[iFEB][iChan]->Write();
	hG_PEcor_err[iFEB][iChan]->Draw("E1");
	hG_PEcor_err[iFEB][iChan]->Write();
	hPEchi[iFEB][iChan]->Write();
	hCPEchi[iFEB][iChan]->Write();
	hPEfullchi[iFEB][iChan]->Write();
      }
    }
  
    std::cout << "writing histograms for " << frun_num << "_" << fsubrun_num << std::endl;
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

    ofile->Close();
    file->Close();
  }
}
