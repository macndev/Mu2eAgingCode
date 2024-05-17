#include <string>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double runningAverage(int lower, int upper, double array[]){
  double avg = 0;
  int range = upper - lower + 1;
  for(int i = lower; i <= upper; i++){
    avg += array[i];
    //if(lower > 980) std::cout << "array [i] in runningAvg function: " << array[i] << std::endl;
  }
  avg /= double(range);
  
  //if(lower > 980) std::cout << "my running average building blocks: range = " << range << " and avg = " << avg << std::endl;

  return avg;
}

const char* short_timechar;

void plotNewTemps(){

  ifstream input("notetempfiles.txt");
  string oneline;

  while(std::getline(input, oneline)){
    const char* cfile = oneline.c_str();
    TFile *file = TFile::Open(cfile);
    std::cout << "file " << cfile << " is opening..." << std::endl;

    TString filename;
    TString runstr;
    TString subrunstr;
    filename = file->GetName();
    int run_num;
    TString short_runstr = filename(filename.Length()-15, filename.Length()-9);
    run_num = short_runstr.Atoi();
    int subrun_num;
    TString short_subrunstr = filename(filename.Length()-8, filename.Length()-5);
    subrun_num = short_subrunstr.Atoi();
    std::cout << "working on file: " << run_num << "_" << subrun_num << std::endl;
    // pad Tstrings with zeros
    if(subrun_num < 10){
      subrunstr = "00" + std::to_string(subrun_num);
    }
    if(subrun_num >= 10){
      subrunstr = "0" + std::to_string(subrun_num);
    }
    if(run_num < 10){
      runstr = "00000" + std::to_string(run_num);
    }
    if((run_num < 100) && (run_num >= 10)){
      runstr = "0000" + std::to_string(run_num);
    }
    if((run_num < 1000) && (run_num >= 100)){
      runstr = "000" + std::to_string(run_num);
    }
    if((run_num < 10000) && (run_num >= 1000)){
      runstr = "00" + std::to_string(run_num);
    }
    //std::cout << "strings: " << short_runstr << "_" << short_subrunstr << std::endl;
    const char* runchar = runstr.Data();
    const char* subrunchar = subrunstr.Data();

    //Create output file. If it already exists, recreate it.
    TFile *ofile = new TFile(Form("output_FEBCMBtemp_v3_vs_spill_run%s_%s.root", runchar, subrunchar),"RECREATE");

    gStyle->SetOptStat(1111111);

    TTree *tree = file->Get<TTree>("otree");
    
    const int nFEB  = 6;
    const int nChan = 64;

    //get info from input files
    int spillNumber;
    float FEBtemp[nFEB];
    float CMBtemp[nFEB*nChan];

    tree->SetBranchAddress("spillNumber", &spillNumber);
    tree->SetBranchAddress("FEBtemp", &FEBtemp);
    tree->SetBranchAddress("CMBtemp", &CMBtemp);

    ofile->cd();

    //make containers for everything
    
    //declare TGraph objects, one for each channel or FEB
    TGraph *ctemps_vs_spill[nFEB*nChan];
    TGraph *ctemps_vs_spill_scaled[nFEB*nChan];
    TGraph *ctemps_vs_spill_running[nFEB*nChan];
    TGraph *ctemps_vs_spill_residual[nFEB*nChan];
    TGraph *ftemps_vs_spill[nFEB];
    TGraph *ftemps_vs_spill_running[nFEB];
    TGraph *ftemps_vs_spill_residual[nFEB];

    TGraph *ctemps_vs_hours[nFEB*nChan];
    TGraph *ctemps_vs_hours_running[nFEB*nChan];
    TGraph *ftemps_vs_hours[nFEB];
    TGraph *ftemps_vs_hours_running[nFEB];

    int nEntries = tree->GetEntries();
    double spill_num[10000];
    double hours[10000];
    double FEB_temp[10000];
    double CMB_temp[10000];
    double CMB_temp_scaled[10000];
    double FEB_temp_running[10000];
    double CMB_temp_running[10000];
    double FEB_temp_residual[10000];
    double CMB_temp_residual[10000];

    int maxchan = nFEB*nChan;

    // get the timestamp of the run
    TNamed *timestamp = (TNamed*)file->Get("start_time");
    const char *timechar = timestamp->GetTitle();
    std::cout << "start time is: " << timechar << std::endl;
    std::string timestring = (std::string)timechar;
    std::string short_timestring(timestring.substr(0,25));
    //std::cout << "string is: " << timestring << std::endl;
    //std::cout << " shortened to " << short_timestring << std::endl;
    short_timechar = short_timestring.c_str();
    std::cout << " shortened to " << short_timechar << std::endl;

    for(int iFEB = 0; iFEB < nFEB; iFEB++){
      for(int iEntry = 0; iEntry < nEntries; iEntry++){
	tree->GetEntry(iEntry);

	spill_num[iEntry] = spillNumber;
	hours[iEntry] = spillNumber * 0.05; // conversion to hours, 1 spill = 3 minutes = 0.05 hours
      
	FEB_temp[iEntry] = FEBtemp[iFEB];

	//std::cout << "FEB temp of FEB " << iFEB << " is " << FEBtemp[iFEB] <<std::endl;
      }
      ftemps_vs_spill[iFEB] = new TGraph(nEntries, spill_num, FEB_temp);
      ftemps_vs_spill[iFEB]->SetTitle(Form("FEBtemp_%i", iFEB));
      ftemps_vs_spill[iFEB]->SetName(Form("FEBtemp_%i", iFEB));
      ftemps_vs_spill[iFEB]->GetXaxis()->SetTitle("Spill Number");
      ftemps_vs_spill[iFEB]->GetYaxis()->SetTitle("FEB Temperature (C)");
      ftemps_vs_spill[iFEB]->SetMarkerColor(kRed);
      ftemps_vs_spill[iFEB]->SetLineColor(kRed);
      ftemps_vs_spill[iFEB]->Write();

      ftemps_vs_hours[iFEB] = new TGraph(nEntries, hours, FEB_temp);
      ftemps_vs_hours[iFEB]->SetTitle(Form("FEBtemp_hours_%i", iFEB));
      ftemps_vs_hours[iFEB]->SetName(Form("FEBtemp_hours_%i", iFEB));
      ftemps_vs_hours[iFEB]->GetXaxis()->SetTitle("Hours Since Start of Subrun");
      ftemps_vs_hours[iFEB]->GetYaxis()->SetTitle("FEB Temperature (C)");
      ftemps_vs_hours[iFEB]->SetMarkerColor(kRed);
      ftemps_vs_hours[iFEB]->SetLineColor(kRed);
      ftemps_vs_hours[iFEB]->Write();

      // do running averages to make the plots smoother
      for(int iEntry = 0; iEntry < nEntries; iEntry++){
	int lowerEnt = iEntry - 3;
	int upperEnt = iEntry + 3;

	if(lowerEnt < 0) lowerEnt = 0; // need to always start with first index
	if(upperEnt >= nEntries) upperEnt = nEntries - 1; // nEntries is 1 index past the end, stop at end

	FEB_temp_running[iEntry] = runningAverage(lowerEnt, upperEnt, FEB_temp);
	FEB_temp_residual[iEntry] = FEB_temp[iEntry] - FEB_temp_running[iEntry];
      }

      //now make running average and residual plots
      ftemps_vs_spill_running[iFEB] = new TGraph(nEntries, spill_num, FEB_temp_running);
      ftemps_vs_spill_running[iFEB]->SetTitle(Form("FEBtemp_running_%i", iFEB));
      ftemps_vs_spill_running[iFEB]->SetName(Form("FEBtemp_running_%i", iFEB));
      ftemps_vs_spill_running[iFEB]->GetXaxis()->SetTitle("Spill Number");
      ftemps_vs_spill_running[iFEB]->GetYaxis()->SetTitle("FEB Temperature (C)");
      ftemps_vs_spill_running[iFEB]->SetMarkerColor(kRed);
      ftemps_vs_spill_running[iFEB]->SetLineColor(kRed);
      ftemps_vs_spill_running[iFEB]->Write();

      ftemps_vs_hours_running[iFEB] = new TGraph(nEntries, hours, FEB_temp_running);
      ftemps_vs_hours_running[iFEB]->SetTitle(Form("FEBtemp_running_hours_%i", iFEB));
      ftemps_vs_hours_running[iFEB]->SetName(Form("FEBtemp_running_hours_%i", iFEB));
      ftemps_vs_hours_running[iFEB]->GetXaxis()->SetTitle("Hours Since Start of Subrun");
      ftemps_vs_hours_running[iFEB]->GetYaxis()->SetTitle("FEB Temperature (C)");
      ftemps_vs_hours_running[iFEB]->SetMarkerColor(kRed);
      ftemps_vs_hours_running[iFEB]->SetLineColor(kRed);
      ftemps_vs_hours_running[iFEB]->Write();

      ftemps_vs_spill_residual[iFEB] = new TGraph(nEntries, spill_num, FEB_temp_residual);
      ftemps_vs_spill_residual[iFEB]->SetTitle(Form("FEBtemp_residual_%i", iFEB));
      ftemps_vs_spill_residual[iFEB]->SetName(Form("FEBtemp_residual_%i", iFEB));
      ftemps_vs_spill_residual[iFEB]->GetXaxis()->SetTitle("Spill Number");
      ftemps_vs_spill_residual[iFEB]->GetYaxis()->SetTitle("FEB Temperature (C)");
      ftemps_vs_spill_residual[iFEB]->SetMarkerColor(kRed);
      ftemps_vs_spill_residual[iFEB]->SetLineColor(kRed);
      ftemps_vs_spill_residual[iFEB]->Write();
    }

    for(int chan = 0; chan < maxchan; chan++){

      for(int iEntry = 0; iEntry < nEntries; iEntry++){
	tree->GetEntry(iEntry);
      
	spill_num[iEntry] = spillNumber;

	CMB_temp[iEntry] = CMBtemp[chan];
	CMB_temp_scaled[iEntry] = CMBtemp[chan]*2.5;
      }
      ctemps_vs_spill[chan] = new TGraph(nEntries, spill_num, CMB_temp);
      ctemps_vs_spill[chan]->SetTitle(Form("CMBtemp_%i", chan));
      ctemps_vs_spill[chan]->SetName(Form("CMBtemp_%i", chan));
      ctemps_vs_spill[chan]->GetXaxis()->SetTitle("Spill Number");
      ctemps_vs_spill[chan]->GetYaxis()->SetTitle("CMB Temperature (C)");
      ctemps_vs_spill[chan]->SetMarkerColor(kBlue);
      ctemps_vs_spill[chan]->SetLineColor(kBlue);
      ctemps_vs_spill[chan]->Write();

      ctemps_vs_hours[chan] = new TGraph(nEntries, hours, CMB_temp);
      ctemps_vs_hours[chan]->SetTitle(Form("CMBtemp_hours_%i", chan));
      ctemps_vs_hours[chan]->SetName(Form("CMBtemp_hours_%i", chan));
      ctemps_vs_hours[chan]->GetXaxis()->SetTitle("Hours Since Start of Subrun");
      ctemps_vs_hours[chan]->GetYaxis()->SetTitle("CMB Temperature (C)");
      ctemps_vs_hours[chan]->SetMarkerColor(kBlue);
      ctemps_vs_hours[chan]->SetLineColor(kBlue);
      ctemps_vs_hours[chan]->Write();

      ctemps_vs_spill_scaled[chan] = new TGraph(nEntries, spill_num, CMB_temp_scaled);
      ctemps_vs_spill_scaled[chan]->SetTitle(Form("CMBtempscaled_%i", chan));
      ctemps_vs_spill_scaled[chan]->SetName(Form("CMBtempscaled_%i", chan));
      ctemps_vs_spill_scaled[chan]->GetXaxis()->SetTitle("Spill Number");
      ctemps_vs_spill_scaled[chan]->GetYaxis()->SetTitle("CMB Temperature * 2.5 (C)");
      ctemps_vs_spill_scaled[chan]->SetMarkerColor(kBlue);
      ctemps_vs_spill_scaled[chan]->SetLineColor(kBlue);
      ctemps_vs_spill_scaled[chan]->Write();

      // do running averages to make the plots smoother
      for(int iEntry = 0; iEntry < nEntries; iEntry++){
	int lowerEnt = iEntry - 3;
	int upperEnt = iEntry + 3;

	if(lowerEnt < 0) lowerEnt = 0; // need to always start with first index
	if(upperEnt >= nEntries) upperEnt = nEntries - 1; // nEntries is 1 index past the end, stop at end

	CMB_temp_running[iEntry] = runningAverage(lowerEnt, upperEnt, CMB_temp);
	CMB_temp_residual[iEntry] = CMB_temp[iEntry] - CMB_temp_running[iEntry];
      }
      // now make running average and residual plots
      ctemps_vs_spill_running[chan] = new TGraph(nEntries, spill_num, CMB_temp_running);
      ctemps_vs_spill_running[chan]->SetTitle(Form("CMBtemp_running_%i", chan));
      ctemps_vs_spill_running[chan]->SetName(Form("CMBtemp_running_%i", chan));
      ctemps_vs_spill_running[chan]->GetXaxis()->SetTitle("Spill Number");
      ctemps_vs_spill_running[chan]->GetYaxis()->SetTitle("CMB Temperature (C)");
      ctemps_vs_spill_running[chan]->SetMarkerColor(kBlue);
      ctemps_vs_spill_running[chan]->SetLineColor(kBlue);
      ctemps_vs_spill_running[chan]->Write();

      ctemps_vs_hours_running[chan] = new TGraph(nEntries, hours, CMB_temp_running);
      ctemps_vs_hours_running[chan]->SetTitle(Form("CMBtemp_running_hours_%i", chan));
      ctemps_vs_hours_running[chan]->SetName(Form("CMBtemp_running_hours_%i", chan));
      ctemps_vs_hours_running[chan]->GetXaxis()->SetTitle("Hours Since Start of Subrun");
      ctemps_vs_hours_running[chan]->GetYaxis()->SetTitle("CMB Temperature (C)");
      ctemps_vs_hours_running[chan]->SetMarkerColor(kBlue);
      ctemps_vs_hours_running[chan]->SetLineColor(kBlue);
      ctemps_vs_hours_running[chan]->Write();

      ctemps_vs_spill_residual[chan] = new TGraph(nEntries, spill_num, CMB_temp_residual);
      ctemps_vs_spill_residual[chan]->SetTitle(Form("CMBtemp_residual_%i", chan));
      ctemps_vs_spill_residual[chan]->SetName(Form("CMBtemp_residual_%i", chan));
      ctemps_vs_spill_residual[chan]->GetXaxis()->SetTitle("Spill Number");
      ctemps_vs_spill_residual[chan]->GetYaxis()->SetTitle("CMB Temperature (C)");
      ctemps_vs_spill_residual[chan]->SetMarkerColor(kBlue);
      ctemps_vs_spill_residual[chan]->SetLineColor(kBlue);
      ctemps_vs_spill_residual[chan]->Write();
    }
    
    // now, overlay the plots so we can look for trends in both
    TMultiGraph* multi[nFEB*nChan];
    TMultiGraph* multi_running[nFEB*nChan];
    TCanvas* comp_temp[nFEB*nChan];
    TCanvas* comp_temp_hours[nFEB*nChan];
    TCanvas* comp_ravg_FEB[nFEB];
    TCanvas* comp_ravg_CMB[nFEB*nChan];
    TCanvas* FEB_0and1;
    TCanvas* FEB_1and2;

    TPaveText *pt = new TPaveText(.8,.85,.95,.95, "NDC");
    pt->AddText(short_timechar);

    for(int chan = 0; chan < maxchan; chan++){
      comp_temp[chan] = new TCanvas(Form("comp_temp_%i", chan), Form("Temperature Comparison for Channel %i", chan), 900, 600);
      comp_temp[chan]->Divide(1,2);
      comp_temp_hours[chan] = new TCanvas(Form("comp_temp_hours_%i", chan), Form("Temperature Comparison for Channel %i", chan), 900, 600);
      comp_temp_hours[chan]->Divide(1,2);
      if(chan < 4){
	comp_ravg_FEB[chan] = new TCanvas(Form("comp_ravg_FEB_%i", chan), Form("Running Avg Residual for FEB %i", chan), 900, 600);
	comp_ravg_FEB[chan]->Divide(1,3);
      }
      comp_ravg_CMB[chan] = new TCanvas(Form("comp_ravg_CMB_%i", chan), Form("Running Avg Residual for CMB %i", chan), 900, 600);
      comp_ravg_CMB[chan]->Divide(1,3);
      multi[chan] = new TMultiGraph(Form("overlay_%i", chan), Form("Temperature Overlay for Channel %i", chan));
      //multi[chan]->GetXaxis()->SetTitle("Spill Number");
      //multi[chan]->GetYaxis()->SetTitle("Temperature (C); Red = FEB, Blue = CMB*2.5");
      multi_running[chan] = new TMultiGraph(Form("overlay_runningavg_%i", chan), Form("Temperature Overlay for Channel %i with Running Avg", chan));
      //multi_running[chan]->GetXaxis()->SetTitle("Spill Number");
      //multi_running[chan]->GetYaxis()->SetTitle("Running Avg Temperature (C); Red = FEB, Blue = CMB*2.5");
      
      // draw correct FEB temp graph, groups of 64 channels (FEB temp is higher, draw first)
      if(chan < nChan){
	comp_temp[chan]->cd(1);
	ftemps_vs_spill_running[0]->Draw();
	comp_temp_hours[chan]->cd(1);
	ftemps_vs_hours_running[0]->Draw();
	//pt->SetFillStyle(0);
	pt->Draw();
	if(chan < 4){
	  comp_ravg_FEB[chan]->cd(1);
	  ftemps_vs_spill[0]->Draw();
	  comp_ravg_FEB[chan]->cd(2);
	  ftemps_vs_spill_running[0]->Draw();
	  comp_ravg_FEB[chan]->cd(3);
	  ftemps_vs_spill_residual[0]->Draw();
	}
	ofile->cd();
	multi[chan]->Add(ftemps_vs_spill[0]);
	multi_running[chan]->Add(ftemps_vs_spill_running[0]);
      }
      if((chan >= nChan) && (chan < nChan*2)){
	comp_temp[chan]->cd(1);
	ftemps_vs_spill_running[1]->Draw();
	comp_temp_hours[chan]->cd(1);
	ftemps_vs_hours_running[1]->Draw();
	if(chan < 4){
	  comp_ravg_FEB[chan]->cd(1);
	  ftemps_vs_spill[1]->Draw();
	  comp_ravg_FEB[chan]->cd(2);
	  ftemps_vs_spill_running[1]->Draw();
	  comp_ravg_FEB[chan]->cd(3);
	  ftemps_vs_spill_residual[1]->Draw();
	}
	ofile->cd();
	multi[chan]->Add(ftemps_vs_spill[1]);
	multi_running[chan]->Add(ftemps_vs_spill_running[1]);
      }
      if((chan >= nChan*2) && (chan < nChan*3)){
	comp_temp[chan]->cd(1);
	ftemps_vs_spill_running[2]->Draw();
	comp_temp_hours[chan]->cd(1);
	ftemps_vs_hours_running[2]->Draw();
	if(chan < 4){
	  comp_ravg_FEB[chan]->cd(1);
	  ftemps_vs_spill[2]->Draw();
	  comp_ravg_FEB[chan]->cd(2);
	  ftemps_vs_spill_running[2]->Draw();
	  comp_ravg_FEB[chan]->cd(3);
	  ftemps_vs_spill_residual[2]->Draw();
	}
	ofile->cd();
	multi[chan]->Add(ftemps_vs_spill[2]);
	multi_running[chan]->Add(ftemps_vs_spill_running[2]);
      }
      if((chan >= nChan*3) && (chan < nChan*4)){
	comp_temp[chan]->cd(1);
	ftemps_vs_spill_running[3]->Draw();
	comp_temp_hours[chan]->cd(1);
	ftemps_vs_hours_running[3]->Draw();
	if(chan < 4){
	  comp_ravg_FEB[chan]->cd(1);
	  ftemps_vs_spill[3]->Draw();
	  comp_ravg_FEB[chan]->cd(2);
	  ftemps_vs_spill_running[3]->Draw();
	  comp_ravg_FEB[chan]->cd(3);
	  ftemps_vs_spill_residual[3]->Draw();
	}
	ofile->cd();
	multi[chan]->Add(ftemps_vs_spill[3]);
	multi_running[chan]->Add(ftemps_vs_spill_running[3]);
      }

      comp_temp[chan]->cd(2);
      ctemps_vs_spill_running[chan]->Draw();
      comp_temp[chan]->Write();
      comp_temp[chan]->Close();
      comp_temp_hours[chan]->cd(2);
      ctemps_vs_hours_running[chan]->Draw();
      comp_temp_hours[chan]->Write();
      comp_temp_hours[chan]->Close();
      if(chan < 4){ 
	comp_ravg_FEB[chan]->Write();
	comp_ravg_FEB[chan]->Close();
      }
      comp_ravg_CMB[chan]->cd(1);
      ctemps_vs_spill[chan]->Draw();
      comp_ravg_CMB[chan]->cd(2);
      ctemps_vs_spill_running[chan]->Draw();
      comp_ravg_CMB[chan]->cd(3);
      ctemps_vs_spill_residual[chan]->Draw();
      comp_ravg_CMB[chan]->Write();
      comp_ravg_CMB[chan]->Close();
      ofile->cd();
      multi[chan]->Add(ctemps_vs_spill_scaled[chan]);
      multi_running[chan]->Add(ctemps_vs_spill_running[chan]);
      multi[chan]->Draw("AP");
      multi_running[chan]->Draw("AP");
      multi[chan]->Write();
      multi_running[chan]->Write();
    }

    FEB_0and1 = new TCanvas("FEB_0and1", "FEB_0and1", 900, 600);
    ftemps_vs_hours_running[0]->Draw();
    ftemps_vs_hours_running[1]->SetMarkerColor(kBlue);
    ftemps_vs_hours_running[1]->SetLineColor(kBlue);
    ftemps_vs_hours_running[1]->Draw("same");
    FEB_0and1->Update();
    ofile->cd();
    FEB_0and1->Write();
    FEB_0and1->Close();

    FEB_1and2 = new TCanvas("FEB_1and2", "FEB_1and2", 900, 600);
    ftemps_vs_hours_running[1]->Draw();
    ftemps_vs_hours_running[2]->SetMarkerColor(kBlue);
    ftemps_vs_hours_running[2]->SetLineColor(kBlue);
    ftemps_vs_hours_running[2]->Draw("same");
    FEB_1and2->Update();
    ofile->cd();
    FEB_1and2->Write();
    FEB_1and2->Close();

    ofile->Close();
    file->Close();
  }
}
