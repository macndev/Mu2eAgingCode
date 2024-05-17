#include <string>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

float getTime(int run){
  float frac = -1;

  if(run == 66) frac = 0.1699;
  if(run == 94) frac = 0.5945;
  if(run == 119) frac = 0.7781;
  if(run == 1009) frac = 0.8027;
  if(run == 1010) frac = 0.8055;
  if(run == 1022) frac = 0.8630;
  if(run == 1027) frac = 0.9534;
  if(run == 1029) frac = 0.9534;
  if(run == 1030) frac = 0.9562;
  if(run == 1031) frac = 0.9726;
  if(run == 1032) frac = 0.9753;
  if(run == 1033) frac = 0.9781;
  if(run == 1034) frac = 0.9863;
  if(run == 1035) frac = 1.0110;
  if(run == 1036) frac = 1.0164;
  if(run == 1037) frac = 1.0740;
  if(run == 1038) frac = 1.0849;
  if(run == 1039) frac = 1.0904;
  if(run == 1040) frac = 1.1041;
  if(run == 1041) frac = 1.1123;
  if(run == 1042) frac = 1.1233;
  if(run == 1043) frac = 1.1315;

  return frac;
}

void PEs_vs_TimeFrac_datesplit(){

 ifstream input("correctedPEfiles_before317.txt");
 //ifstream input("first.txt");
 string oneline;

 //Create output file. If it already exists, recreate it.
 TFile *ofile = new TFile("PEs_vs_time_plots_beforeMarch_label.root","RECREATE");

 int iFEB  = 2;
 int iChan = 64;
 int nhists = iFEB * iChan;

 TH2F *hPE_un[128];
 TH2F *hPE_cor[128];
 TH2F *hPE_un_frac[128];
 TH2F *hPE_cor_frac[128];
 TH2F *hPE_G_un[128];
 TH2F *hPE_G_cor[128];
 TH2F *hPE_G_un_frac[128];
 TH2F *hPE_G_cor_frac[128];
 TGraphErrors* tgPE_G_un_err[128];
 TGraphErrors* tgPE_G_cor_err[128];
 TGraphErrors* tgPE_G_un_err_frac[128];
 TGraphErrors* tgPE_G_cor_err_frac[128];
 TGraphErrors* tgPE_G_cor_avg;
 TGraphErrors* tgPE_G_cor_avg_frac;
 TH1F *hPE_chi[128];
 TH1F *htotal_chi;
 TH1F *hPE_slopes;
 TH1F *hPE_slope_frac;
 TH1F *hslopes1D;
 TH1F *hFEB0slopes;
 TH1F *hFEB1slopes;
 TCanvas* page[128][4];
 //TCanvas* avg[2];
 TCanvas* slopeslabel[2];

 gStyle->SetOptStat(1111111);

 //avg[0] = new TCanvas("avg_vs_run", "avg_vs_run", 900, 600);
 //avg[1] = new TCanvas("avg_vs_frac", "avg_vs_frac", 900, 600);
 slopeslabel[0] = new TCanvas("slopes_labeled", "slopes_labeled", 900, 600);
 slopeslabel[1] = new TCanvas("slopesfrac_labeled", "slopesfrac_labeled", 900, 600);

 TPaveText *pt = new TPaveText(.55,.85,.75,.95, "NDC");
 pt->AddText("July 8, 2021 to Feb 26, 2022");

 //initialize empty hists only once, set style
 for(int ihist = 0; ihist < nhists; ihist++){
   string sname;
   const char* cname;
   int ncansets = 4;
   for(int ican = 0; ican < ncansets; ican++) {
     sname = "can" + std::to_string(ihist) + "_" + std::to_string(ican);
     cname = sname.c_str();
     page[ihist][ican] = new TCanvas(cname, cname, 900, 600);
   }

   hPE_un[ihist] = new TH2F(Form("hPE_un_%i",ihist), "Uncorrected PE vs Run", 1050, 0., 1050., 200, 30., 50.);
   hPE_un[ihist]->GetXaxis()->SetTitle("Run Number");
   hPE_un[ihist]->GetYaxis()->SetTitle("PE Yield (peak value from fit)"); 

   hPE_cor[ihist] = new TH2F(Form("hPE_cor_%i",ihist), "Corrected PE vs Run", 1050, 0., 1050., 200, 30., 50.);
   hPE_cor[ihist]->GetXaxis()->SetTitle("Run Number");
   hPE_cor[ihist]->GetYaxis()->SetTitle("PE Yield (corrected peak value from fit)"); 

   hPE_un_frac[ihist] = new TH2F(Form("hPE_un_frac_%i",ihist), "Uncorrected PE vs Time Fraction", 1000, 0., 1.2, 200, 30., 50.);
   hPE_un_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 5/7/2021)");
   hPE_un_frac[ihist]->GetYaxis()->SetTitle("PE Yield (peak value from fit)"); 

   hPE_cor_frac[ihist] = new TH2F(Form("hPE_cor_frac_%i",ihist), "Corrected PE vs Time Fraction", 1000, 0., 1.2, 200, 30., 50.);
   hPE_cor_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 5/7/2021)");
   hPE_cor_frac[ihist]->GetYaxis()->SetTitle("PE Yield (corrected peak value from fit)"); 

   hPE_G_un[ihist] = new TH2F(Form("hPE_un_%i",ihist), "Uncorrected PE vs Run", 1050, 0., 1050., 200, 30., 50.);
   hPE_G_un[ihist]->GetXaxis()->SetTitle("Run Number");
   hPE_G_un[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean from fit)"); 

   hPE_G_cor[ihist] = new TH2F(Form("hPE_G_cor_%i",ihist), "Corrected PE vs Run", 1050, 0., 1050., 200, 30., 50.);
   hPE_G_cor[ihist]->GetXaxis()->SetTitle("Run Number");
   hPE_G_cor[ihist]->GetYaxis()->SetTitle("PE Yield (corrected Gaussian mean value from fit)"); 

   hPE_G_un_frac[ihist] = new TH2F(Form("hPE_G_un_frac_%i",ihist), "Uncorrected PE vs Time Fraction", 1000, 0., 1.1, 200, 30., 50.);
   hPE_G_un_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 5/7/2021)");
   hPE_G_un_frac[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean from fit)"); 

   hPE_G_cor_frac[ihist] = new TH2F(Form("hPE_G_cor_frac_%i",ihist), "Corrected PE vs Time Fraction", 1000, 0., 1.1, 200, 30., 50.);
   hPE_G_cor_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 5/7/2021)");
   hPE_G_cor_frac[ihist]->GetYaxis()->SetTitle("PE Yield (corrected Gaussian mean value from fit)");

   hPE_chi[ihist] = new TH1F(Form("hPE_chi_%i",ihist), "Chi2 from PE Fits", 500, 0., 500.);
   hPE_chi[ihist]->GetXaxis()->SetTitle("Chi2 from Fit");
   hPE_chi[ihist]->GetYaxis()->SetTitle("Counts");

   htotal_chi = new TH1F(Form("htotal_chi_%i",ihist), "Chi2 from PE Fits", 500, 0., 500.);
   htotal_chi->GetXaxis()->SetTitle("Chi2 from Fit");
   htotal_chi->GetYaxis()->SetTitle("Counts");


   hPE_slopes = new TH1F("hPE_slopes", "Aging Slope vs Channel Number", 128, 0., 128.);
   hPE_slopes->GetXaxis()->SetTitle("Channel Number");
   hPE_slopes->GetYaxis()->SetTitle("Slope from Linear Aging Fit");
   hPE_slopes->GetYaxis()->SetRangeUser(-20.,20.);
   //-8 to 0

   hPE_slope_frac = new TH1F("hPE_slope_frac", "Slope in Percentage of PE vs Channel Number", 128, 0., 128.);
   hPE_slope_frac->GetXaxis()->SetTitle("Channel Number");
   hPE_slope_frac->GetYaxis()->SetTitle("Slope/PE Yield Percentage");
   hPE_slope_frac->GetYaxis()->SetRangeUser(-50.,50.);
   //-20 to 0

   hslopes1D = new TH1F("hslopes1D", "Aging Slopes", 100, -20., 20.);
   hslopes1D->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
   hslopes1D->GetYaxis()->SetTitle("Counts");
   //-20 to 10

   hFEB0slopes = new TH1F("hFEB0slopes", "Aging Slopes", 100, -20., 20.);
   hFEB0slopes->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
   hFEB0slopes->GetYaxis()->SetTitle("Counts");

   hFEB1slopes = new TH1F("hFEB1slopes", "Aging Slopes", 100, -20., 10.);
   hFEB1slopes->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
   hFEB1slopes->GetYaxis()->SetTitle("Counts");
 }

 Float_t aPE_G_un[128][100];
 Float_t aPE_G_cor[128][100];
 Float_t aPE_G_err[128][100];
 Float_t aPE_G_chi[128][100];
 Float_t aempty_xerr[100];
 Float_t arun_bins[100];
 Float_t atimefrac_bins[100];
 Float_t aslopes[128];
 Float_t aslope_frac[128];

 int ifile = 0;
 while(std::getline(input, oneline)){
   const char* cfile = oneline.c_str();
   TFile *file = TFile::Open(cfile);

   TString filename;
   filename = file->GetName();
   std::string sfilename;
   sfilename = (std::string) filename.Data();
   std::cout << "filename string: " << sfilename << " length: " << sfilename.length() << std::endl;
   int run_num;
   int subrun_num;
   int run_pos;
   int und_pos;
   int dot_pos;
   int run_len;
   int sr_len;
   std::string run_s;
   std::string subrun_s;
   run_pos = sfilename.find("run");
   und_pos = sfilename.find("_",run_pos);
   dot_pos = sfilename.find(".",und_pos);
   run_len = und_pos - (run_pos+3);
   run_s = sfilename.substr(run_pos+3,run_len);
   sr_len = dot_pos - (und_pos+1);
   subrun_s = sfilename.substr(und_pos+1, sr_len);
   run_num = std::stoi(run_s);
   Float_t frun_num = (Float_t) run_num;
   subrun_num = std::stoi(subrun_s);
   std::cout << "run number: " << run_s << ", run length: " << run_len << "; subrun number: " << subrun_s << ", subrun length: " << sr_len << std::endl;

   //fill x bins for both versions of plots, fill empty x errors
   Float_t ftimefrac = getTime(run_num);
   arun_bins[ifile] = frun_num;
   atimefrac_bins[ifile] = ftimefrac;
   Float_t fzero = 0.;
   aempty_xerr[ifile] = fzero;

   //declare temp histograms
   TH2F *htemp_u = new TH2F("htemp_u", "temporary_un", 1050, 0., 1050., 200, 30., 50.);
   TH2F *htemp_ug = new TH2F("htemp_ug", "temporary_un_g", 1050, 0., 1050., 200, 30., 50.);
   TH2F *htemp_c = new TH2F("htemp_c", "temporary_cor", 1050, 0., 1050., 200, 30., 50.);
   TH2F *htemp_cg = new TH2F("htemp_cg", "temporary_cor_g", 1050, 0., 1050., 200, 30., 50.);
   TH1F *htemp_cge = new TH1F("htemp_cge", "temporary_cor_gerr", 1050, 0., 1050.);
   TH1F *htemp_chi = new TH1F("htemp_chi", "temporary_chi2", 1050, 0., 1050.);

   //get all of the interesting quantities from input files
   for(int ihist = 0; ihist < nhists; ihist++){
     //TH2F *htemp_u = new TH2F("htemp_u", "temporary_un", 1050, 0., 1050., 200, 30., 50.);
     htemp_u = (TH2F*)file->Get(Form("PE_vs_run_%i",ihist));
     float run_num = htemp_u->GetMean(1);
     float PE_un   = htemp_u->GetMean(2);
     int  irun_num = (int) run_num;
     htemp_u->Delete();

     //TH2F *htemp_ug = new TH2F("htemp_ug", "temporary_un_g", 1050, 0., 1050., 200, 30., 50.);
     htemp_ug = (TH2F*)file->Get(Form("gaus_vs_run_%i",ihist));
     Float_t PE_un_g = htemp_ug->GetMean(2);
     aPE_G_un[ihist][ifile] = PE_un_g;
     htemp_ug->Delete();
     

     //TH2F *htemp_c = new TH2F("htemp_c", "temporary_cor", 1050, 0., 1050., 200, 30., 50.);
     htemp_c = (TH2F*)file->Get(Form("PE_corrected_%i",ihist));
     float PE_cor = htemp_c->GetMean(2);
     htemp_c->Delete();

     //TH2F *htemp_cg = new TH2F("htemp_cg", "temporary_cor_g", 1050, 0., 1050., 200, 30., 50.);
     htemp_cg = (TH2F*)file->Get(Form("gaus_corrected_%i",ihist));
     Float_t PE_cor_g = htemp_cg->GetMean(2);
     if(PE_cor_g < 30.){
       std::cout << "low PE yield in this channel (channel " << ihist << " in run " << run_num << " subrun " << subrun_num << "), PE yield is: " << PE_cor_g << std::endl;
     }
     aPE_G_cor[ihist][ifile] = PE_cor_g;
     htemp_cg->Delete();

     //TH1F *htemp_cge = new TH1F("htemp_cge", "temporary_cor_gerr", 1050, 0., 1050.);
     htemp_cge = (TH1F*)file->Get(Form("gaus_vs_run_cor_err_%i",ihist));
     Float_t PE_g_err = htemp_cge->GetBinError(irun_num);
     // if(PE_g_err > 100.){
     //   std::cout << "error for this channel (channel " << ihist << " in run " << run_num << " subrun " << subrun_num << ")  is: " << PE_g_err << std::endl;
     // }
     // if(PE_un == 0.){
     //   std::cout << "zero PE channel (channel " << ihist << " in run " << run_num << " subrun " << subrun_num << ")  zero: " << PE_un << " with error: " << PE_g_err << std::endl;
     // }
     aPE_G_err[ihist][ifile] = PE_g_err;
     htemp_cge->Delete();

     htemp_chi = (TH1F*)file->Get(Form("gaus_vs_run_err_%i",ihist));
     Float_t PE_chi = htemp_chi->GetBinError(irun_num);
     aPE_G_chi[ihist][ifile] = PE_chi;
     htemp_chi->Delete();

     float timefrac = getTime(run_num);

     //fill histograms
     if((PE_un != 0.) && (PE_g_err < 100.)){
       ofile->cd();
       hPE_un[ihist]->SetMarkerStyle(20);
       hPE_cor[ihist]->SetMarkerStyle(20);
       hPE_un_frac[ihist]->SetMarkerStyle(20);
       hPE_cor_frac[ihist]->SetMarkerStyle(20);
       hPE_un[ihist]->Fill(run_num,PE_un);
       hPE_cor[ihist]->Fill(run_num,PE_cor);
       hPE_un_frac[ihist]->Fill(timefrac,PE_un);
       hPE_cor_frac[ihist]->Fill(timefrac,PE_cor);
       hPE_G_un[ihist]->SetMarkerStyle(20);
       hPE_G_cor[ihist]->SetMarkerStyle(20);
       hPE_G_un_frac[ihist]->SetMarkerStyle(20);
       hPE_G_cor_frac[ihist]->SetMarkerStyle(20);
       hPE_G_un[ihist]->Fill(run_num,PE_un_g);
       hPE_G_cor[ihist]->Fill(run_num,PE_cor_g);
       hPE_G_un_frac[ihist]->Fill(timefrac,PE_un_g);
       hPE_G_cor_frac[ihist]->Fill(timefrac,PE_cor_g);
       hPE_chi[ihist]->Fill(PE_chi);
       htotal_chi->Fill(PE_chi);
     }
   }
   file->Close(); 
   ifile++;
 }

 Float_t PE_G_cor_avg[ifile];
 Float_t PE_G_avg_err[ifile];
 int num_full[ifile];
 //write histograms only once
 for(int ihist = 0; ihist < nhists; ihist++){
   ofile->cd();
   hPE_un[ihist]->Write();
   hPE_cor[ihist]->Write();
   hPE_un_frac[ihist]->Write();
   hPE_cor_frac[ihist]->Write();
   hPE_G_un[ihist]->Write();
   hPE_G_cor[ihist]->Write();
   hPE_G_un_frac[ihist]->Write();
   hPE_G_cor_frac[ihist]->Write();
   hPE_chi[ihist]->Write();

   Float_t thisPE_G_un[ifile];
   Float_t thisPE_G_cor[ifile];
   Float_t thisPE_G_err[ifile];
   for(int i = 0; i < ifile; i++){
     if((aPE_G_un[ihist][i] == 0.) || (aPE_G_err[ihist][i] > 100.) || (aPE_G_chi[ihist][i] >= 90.47)){
       thisPE_G_un[i] = 0.;
       thisPE_G_cor[i] = 0.;
       thisPE_G_err[i] = 0.;
     }
     if((aPE_G_un[ihist][i] != 0.) && (aPE_G_err[ihist][i] < 100.) && (aPE_G_chi[ihist][i] < 90.47)){
       thisPE_G_un[i] = aPE_G_un[ihist][i];
       thisPE_G_cor[i] = aPE_G_cor[ihist][i];
       thisPE_G_err[i] = aPE_G_err[ihist][i];
       PE_G_cor_avg[i] += aPE_G_cor[ihist][i];
       PE_G_avg_err[i] += aPE_G_err[ihist][i];
       num_full[i]     += 1;
     }
   }

   std::cout << "ifile = " << ifile << std::endl;
   /*page[ihist][0]->Close();
   page[ihist][1]->Close();
   page[ihist][2]->Close();
   page[ihist][3]->Close();*/
   std::cout << "ihist = " << ihist << std::endl;
   page[ihist][0]->cd();
   tgPE_G_un_err[ihist] = new TGraphErrors(ifile,arun_bins,thisPE_G_un,aempty_xerr,thisPE_G_err);
   tgPE_G_un_err[ihist]->SetTitle(Form("Uncorrected PE Yield vs Run (Channel %i)",ihist));
   tgPE_G_un_err[ihist]->GetXaxis()->SetTitle("Run Number");
   tgPE_G_un_err[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean)"); 
   tgPE_G_un_err[ihist]->GetYaxis()->SetRangeUser(20.,50.);
   tgPE_G_un_err[ihist]->SetMarkerStyle(20);
   tgPE_G_un_err[ihist]->Draw("AP");
   page[ihist][0]->Update();
   ofile->cd();
   page[ihist][0]->Write();
   page[ihist][0]->Close();
   page[ihist][1]->cd();
   tgPE_G_cor_err[ihist] = new TGraphErrors(ifile,arun_bins,thisPE_G_cor,aempty_xerr,thisPE_G_err);
   tgPE_G_cor_err[ihist]->SetTitle(Form("Corrected PE Yield vs Run (Gaus mean, channel %i)", ihist));
   tgPE_G_cor_err[ihist]->GetXaxis()->SetTitle("Run Number");
   tgPE_G_cor_err[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean)"); 
   tgPE_G_cor_err[ihist]->GetYaxis()->SetRangeUser(20.,50.);
   tgPE_G_cor_err[ihist]->SetMarkerStyle(20);
   tgPE_G_cor_err[ihist]->Draw("AP");
   page[ihist][1]->Update();
   ofile->cd();
   gStyle->SetOptFit(1);
   page[ihist][1]->Write();
   page[ihist][1]->Close();
   page[ihist][2]->cd();
   tgPE_G_un_err_frac[ihist] = new TGraphErrors(ifile,atimefrac_bins,thisPE_G_un,aempty_xerr,thisPE_G_err);
   tgPE_G_un_err_frac[ihist]->SetTitle(Form("Uncorrected PE Yield vs Time Frac (Gaus mean, channel %i)",ihist));
   tgPE_G_un_err_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 5/7/2021)");
   tgPE_G_un_err_frac[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean)"); 
   tgPE_G_un_err_frac[ihist]->GetYaxis()->SetRangeUser(20.,50.);
   tgPE_G_un_err_frac[ihist]->SetMarkerStyle(20);
   tgPE_G_un_err_frac[ihist]->Draw("AP");
   page[ihist][2]->Update();
   ofile->cd();
   page[ihist][2]->Write();
   page[ihist][2]->Close();
   page[ihist][3]->cd();
   tgPE_G_cor_err_frac[ihist] = new TGraphErrors(ifile,atimefrac_bins,thisPE_G_cor,aempty_xerr,thisPE_G_err);
   tgPE_G_cor_err_frac[ihist]->SetTitle(Form("Corrected PE Yield vs Time Frac (Gaus mean, channel %i)",ihist));
   tgPE_G_cor_err_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 5/7/2021)");
   tgPE_G_cor_err_frac[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean)"); 
   //tgPE_G_cor_err_frac[ihist]->GetYaxis()->SetRangeUser(20.,50.);
   tgPE_G_cor_err_frac[ihist]->SetMarkerStyle(20);
   TF1 *newfit = new TF1("newfit","pol1",0.16,0.81);
   tgPE_G_cor_err_frac[ihist]->Fit("newfit","R");
   //TF1 *fit = (TF1*)tgPE_G_cor_err_frac[ihist]->GetListOfFunctions()->FindObject("pol1");
   float fitp0 = newfit->GetParameter(0);
   float fitp1 = newfit->GetParameter(1);
   aslopes[ihist] = fitp1;
   std::cout << "slope for this channel is: " << aslopes[ihist] << std::endl;
   float lowerbound = fitp0 + (fitp1*0.16);
   float upperbound = fitp0 + (fitp1*0.81);
   std::cout << "bounds for this channel are: " << lowerbound << ", " << upperbound << std::endl;
   tgPE_G_cor_err_frac[ihist]->GetYaxis()->SetRangeUser(lowerbound-2, upperbound+2);

   int nonzerobins = 0;
   Float_t avgPE = 0.;
   for(int i = 0; i < ifile; i++){
     if(thisPE_G_cor[i] != 0){
       avgPE += thisPE_G_cor[i];
       nonzerobins += 1;
       if((thisPE_G_cor[i] > 1.) && (thisPE_G_cor[i] < 32.)){std::cout << "Low Nonzero PE Yield, channel " << ihist << " PE = " << thisPE_G_cor[i] << " from file " << i << std::endl;}
       if(thisPE_G_cor[i] > 45.){std::cout << "High PE Yield, channel " << ihist << " PE = " << thisPE_G_cor[i] << " from file " << i << std::endl;}
     }
   }
   avgPE /= nonzerobins;


   if(avgPE > 0.){aslope_frac[ihist] = (aslopes[ihist]/avgPE)*100;}
   if(avgPE == 0.){aslope_frac[ihist] = -50.;}
   if(isnan(avgPE) > 0.){
     aslope_frac[ihist] = -50.;
     aslopes[ihist] = -50.;
   }
   std::cout << "slope fraction for this channel is: " << aslope_frac[ihist] << " = " << aslopes[ihist] << " / " << avgPE << std::endl;
     std::cout << "WAS : " << aslopes[ihist] << " / " << tgPE_G_cor_err_frac[ihist]->GetMean(2) << std::endl;
   if((aslope_frac[ihist] > 10.) || (aslope_frac[ihist] < -20.)){
     std::cout << "BAD CHANNEL: " << ihist << " slope is " << aslope_frac[ihist] << std::endl;
   }
   tgPE_G_cor_err_frac[ihist]->Draw("AP");
   page[ihist][3]->Update();
   ofile->cd();
   page[ihist][3]->Write();
   page[ihist][3]->Close();

   hPE_slopes->SetBinContent(ihist+1,aslopes[ihist]);
   hPE_slope_frac->SetBinContent(ihist+1,aslope_frac[ihist]);
   hslopes1D->Fill(aslope_frac[ihist]);
   if(ihist < 64.) {hFEB0slopes->Fill(aslope_frac[ihist]);}
   else {hFEB1slopes->Fill(aslope_frac[ihist]);}
 } 
 for(int i = 0; i < ifile; i++){
   PE_G_cor_avg[i] /= num_full[i];
   PE_G_avg_err[i] /= num_full[i];
 }

 /*avg[0]->cd();
 tgPE_G_cor_avg = new TGraphErrors(ifile,arun_bins,PE_G_cor_avg,aempty_xerr,PE_G_avg_err);
 tgPE_G_cor_avg->SetTitle("Avg Corrected PE Yield vs Run Number (Gaus mean)");
 tgPE_G_cor_avg->GetXaxis()->SetTitle("Run Number");
 tgPE_G_cor_avg->GetYaxis()->SetTitle("PE Yield (Gaussian mean)"); 
 tgPE_G_cor_avg->GetYaxis()->SetRangeUser(20.,50.);
 tgPE_G_cor_avg->SetMarkerStyle(20);
 tgPE_G_cor_avg->Fit("pol1");
 tgPE_G_cor_avg->Draw("AP");
 avg[0]->Update();
 ofile->cd();
 avg[0]->Write();
 avg[0]->Close();

 avg[1]->cd();
 tgPE_G_cor_avg_frac = new TGraphErrors(ifile,atimefrac_bins,PE_G_cor_avg,aempty_xerr,PE_G_avg_err);
 tgPE_G_cor_avg_frac->SetTitle("Avg Corrected PE Yield vs Time Frac (Gaus mean)");
 tgPE_G_cor_avg_frac->GetXaxis()->SetTitle("Time Fraction (years since 5/7/2021)");
 tgPE_G_cor_avg_frac->GetYaxis()->SetTitle("PE Yield (Gaussian mean)"); 
 tgPE_G_cor_avg_frac->GetYaxis()->SetRangeUser(20.,50.);
 tgPE_G_cor_avg_frac->SetMarkerStyle(20);
 tgPE_G_cor_avg_frac->Fit("pol1");
 tgPE_G_cor_avg_frac->Draw("AP");
 avg[1]->Update();
 ofile->cd();
 avg[1]->Write();
 avg[1]->Close();*/

 slopeslabel[0]->cd();
 hPE_slopes->Draw();
 pt->Draw("same");
 slopeslabel[0]->Update();
 ofile->cd();
 slopeslabel[0]->Write();
 slopeslabel[0]->Close();

 slopeslabel[1]->cd();
 hPE_slope_frac->Draw();
 pt->Draw("same");
 slopeslabel[1]->Update();
 ofile->cd();
 slopeslabel[1]->Write();
 slopeslabel[1]->Close();

 htotal_chi->Write();
 hPE_slopes->Write();
 hPE_slope_frac->Write();
 hslopes1D->Write();
 hFEB0slopes->Write();
 hFEB1slopes->Write();
}

