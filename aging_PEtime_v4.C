#include <string>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const int cnFEB = 6;
const int cnchan = 64;
const int ctotalchan = cnFEB * cnchan;

float getTime(int run, int subrun){
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

  float timefrac_50hours = 0.0057;

  frac += (timefrac_50hours * subrun);

  return frac;
}

float getTime_168(int run, int subrun){
  float frac = -1;

  if(run == 1053) frac = 0.000;
  if(run == 1054) frac = 0.030;
  if(run == 1059) frac = 0.038;
  if(run == 1066) frac = 0.058;
  if(run == 1079) frac = 0.066;
  if(run == 1080) frac = 0.068;
  if(run == 1081) frac = 0.071;
  if(run == 1082) frac = 0.074;
  if(run == 1084) frac = 0.077;
  if(run == 1088) frac = 0.085;
  if(run == 1089) frac = 0.090;
  if(run == 1090) frac = 0.093;
  if(run == 1091) frac = 0.096;
  if(run == 1116) frac = 0.115;
  if(run == 1124) frac = 0.129;
  if(run == 1133) frac = 0.132;
  if(run == 1134) frac = 0.134;
  if(run == 1137) frac = 0.148;
  if(run == 1138) frac = 0.151;
  if(run == 1141) frac = 0.164;
  if(run == 1142) frac = 0.167;
  if(run == 1146) frac = 0.173;
  if(run == 1148) frac = 0.184;
  if(run == 1149) frac = 0.192;
  if(run == 1150) frac = 0.200;
  if(run == 1152) frac = 0.203;
  if(run == 1154) frac = 0.205;
  if(run == 1159) frac = 0.208;
  if(run == 1160) frac = 0.216;
  if(run == 1161) frac = 0.219;
  if(run == 1167) frac = 0.230;
  if(run == 1168) frac = 0.238;
  if(run == 1169) frac = 0.244;
  if(run == 1170) frac = 0.249;
  if(run == 1171) frac = 0.255;
  if(run == 1172) frac = 0.258;
  if(run == 1173) frac = 0.260;
  if(run == 1176) frac = 0.263;
  if(run == 1177) frac = 0.263;
  if(run == 1178) frac = 0.268;
  if(run == 1181) frac = 0.268;
  if(run == 1182) frac = 0.274;
  if(run == 1208) frac = 0.279;
  if(run == 1209) frac = 0.282;
  if(run == 1212) frac = 0.282;
  if(run == 1217) frac = 0.296;

  float timefrac_50hours = 0.0057;

  frac += (timefrac_50hours * subrun);

  return frac;
}

void aging_PEtime_v4(){

 ifstream input("aging_fitPE_files_168.txt");
 //ifstream input("first.txt");
 string oneline;

 //Create output file. If it already exists, recreate it.
 TFile *ofile = new TFile("PEtime_plots_module168_const_re.root","RECREATE");

 TH2F *hPE_un[ctotalchan];
 TH2F *hPE_cor[ctotalchan];
 TH2F *hPE_un_frac[ctotalchan];
 TH2F *hPE_cor_frac[ctotalchan];
 TH2F *hPE_G_un[ctotalchan];
 TH2F *hPE_G_cor[ctotalchan];
 TH2F *hPE_G_un_frac[ctotalchan];
 TH2F *hPE_G_cor_frac[ctotalchan];
 TGraphErrors* tgPE_G_un_err[ctotalchan];
 TGraphErrors* tgPE_G_cor_err[ctotalchan];
 TGraphErrors* tgPE_G_un_err_frac[ctotalchan];
 TGraphErrors* tgPE_G_cor_err_frac[ctotalchan];
 TGraphErrors* tgPE_G_cor_avg;
 TGraphErrors* tgPE_G_cor_avg_frac;
 TH1F *hPE_chi[ctotalchan];
 TH1F *htotal_chi;
 TH1F *hCPE_chi[ctotalchan];
 TH1F *htotal_Cchi;
 TH1F *hPE_slopes;
 TH1F *hPEhist_slopes;
 TH1F *hPEhist_slopes_un;
 TH1F *hPE_slope_frac;
 TH1F *hslopes1D;
 TH1F *hFEB0slopes;
 TH1F *hFEB1slopes;
 TH1F *hFEB2slopes;
 TH1F *hFEB3slopes;
 TH1F *hFEB4slopes;
 TH1F *hFEB5slopes;
 TH1F *hPE_slopes_un;
 TH1F *hPE_slope_frac_un;
 TH1F *hslopes1D_un;
 TH1F *hFEB0slopes_un;
 TH1F *hFEB1slopes_un;
 TH1F *hFEB2slopes_un;
 TH1F *hFEB3slopes_un;
 TH1F *hFEB4slopes_un;
 TH1F *hFEB5slopes_un;
 TCanvas* page[ctotalchan][4];
 //TCanvas* avg[2];
 TCanvas* slopeslabel[2];
 TCanvas* slopeslabel_un[2];

 gStyle->SetOptStat(1111111);

 //avg[0] = new TCanvas("avg_vs_run", "avg_vs_run", 900, 600);
 //avg[1] = new TCanvas("avg_vs_frac", "avg_vs_frac", 900, 600);
 slopeslabel[0] = new TCanvas("slopes_labeled", "slopes_labeled", 900, 600);
 slopeslabel[1] = new TCanvas("slopesfrac_labeled", "slopesfrac_labeled", 900, 600);
 slopeslabel_un[0] = new TCanvas("slopes_labeled_un", "slopes_labeled_un", 900, 600);
 slopeslabel_un[1] = new TCanvas("slopesfrac_labeled_un", "slopesfrac_labeled_un", 900, 600);

 TPaveText *pt = new TPaveText(.55,.85,.75,.95, "NDC");
 pt->AddText("July 1, 2022 to Oct 17, 2022");

 htotal_chi = new TH1F("htotal_chi", "Reduced Chi2 from Uncorrected PE Fits", 500, 0., 100.);
 htotal_chi->GetXaxis()->SetTitle("Chi2 from Fit");
 htotal_chi->GetYaxis()->SetTitle("Counts");

 htotal_Cchi = new TH1F("htotal_Cchi", "Reduced Chi2 from Corrected PE Fits", 500, 0., 100.);
 htotal_Cchi->GetXaxis()->SetTitle("Chi2 from Fit");
 htotal_Cchi->GetYaxis()->SetTitle("Counts");

 hPE_slopes = new TH1F("hPE_slopes", "Aging Slope vs Channel Number", 384, 0., 384.);
 hPE_slopes->GetXaxis()->SetTitle("Channel Number");
 hPE_slopes->GetYaxis()->SetTitle("Slope from Linear Aging Fit");
 hPE_slopes->GetYaxis()->SetRangeUser(-20.,20.);
 //-8 to 0

 hPEhist_slopes = new TH1F("hPEhist_slopes", "Hist Aging Slope vs Channel Number", 384, 0., 384.);
 hPEhist_slopes->GetXaxis()->SetTitle("Channel Number");
 hPEhist_slopes->GetYaxis()->SetTitle("Slope from Linear Aging Fit, Cor Histograms");
 hPEhist_slopes->GetYaxis()->SetRangeUser(-20.,20.);

 hPEhist_slopes_un = new TH1F("hPEhist_slopes_un", "Hist Aging Slope vs Channel Number", 384, 0., 384.);
 hPEhist_slopes_un->GetXaxis()->SetTitle("Channel Number");
 hPEhist_slopes_un->GetYaxis()->SetTitle("Slope from Linear Aging Fit, Uncor Histograms");
 hPEhist_slopes_un->GetYaxis()->SetRangeUser(-20.,20.);

 hPE_slope_frac = new TH1F("hPE_slope_frac", "Slope in Percentage of PE vs Channel Number", 384, 0., 384.);
 hPE_slope_frac->GetXaxis()->SetTitle("Channel Number");
 hPE_slope_frac->GetYaxis()->SetTitle("Slope/PE Yield Percentage");
 hPE_slope_frac->GetYaxis()->SetRangeUser(-50.,50.);
 //-20 to 0

 hslopes1D = new TH1F("hslopes1D", "Aging Slopes", 100, -25., 25.);
 hslopes1D->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hslopes1D->GetYaxis()->SetTitle("Counts");
 //-20 to 10

 hFEB0slopes = new TH1F("hFEB0slopes", "Aging Slopes", 100, -25., 25.);
 hFEB0slopes->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB0slopes->GetYaxis()->SetTitle("Counts");

 hFEB1slopes = new TH1F("hFEB1slopes", "Aging Slopes", 100, -25., 25.);
 hFEB1slopes->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB1slopes->GetYaxis()->SetTitle("Counts");

 hFEB2slopes = new TH1F("hFEB2slopes", "Aging Slopes", 100, -25., 25.);
 hFEB2slopes->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB2slopes->GetYaxis()->SetTitle("Counts");

 hFEB3slopes = new TH1F("hFEB3slopes", "Aging Slopes", 100, -25., 25.);
 hFEB3slopes->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB3slopes->GetYaxis()->SetTitle("Counts");

 hFEB4slopes = new TH1F("hFEB4slopes", "Aging Slopes", 100, -25., 25.);
 hFEB4slopes->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB4slopes->GetYaxis()->SetTitle("Counts");

 hFEB5slopes = new TH1F("hFEB5slopes", "Aging Slopes", 100, -25., 25.);
 hFEB5slopes->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB5slopes->GetYaxis()->SetTitle("Counts");

 hPE_slopes_un = new TH1F("hPE_slopes_un", "Uncorrected Aging Slope vs Channel Number", 384, 0., 384.);
 hPE_slopes_un->GetXaxis()->SetTitle("Channel Number");
 hPE_slopes_un->GetYaxis()->SetTitle("Slope from Linear Aging Fit");
 hPE_slopes_un->GetYaxis()->SetRangeUser(-20.,20.);
 //-8 to 0

 hPE_slope_frac_un = new TH1F("hPE_slope_frac_un", "Uncorrected Slope in Percentage of PE vs Channel Number", 384, 0., 384.);
 hPE_slope_frac_un->GetXaxis()->SetTitle("Channel Number");
 hPE_slope_frac_un->GetYaxis()->SetTitle("Slope/PE Yield Percentage");
 hPE_slope_frac_un->GetYaxis()->SetRangeUser(-50.,50.);
 //-20 to 0

 hslopes1D_un = new TH1F("hslopes1D_un", "Uncorrected Aging Slopes", 100, -25., 25.);
 hslopes1D_un->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hslopes1D_un->GetYaxis()->SetTitle("Counts");
 //-20 to 10

 hFEB0slopes_un = new TH1F("hFEB0slopes_un", "Uncorrected Aging Slopes", 100, -25., 25.);
 hFEB0slopes_un->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB0slopes_un->GetYaxis()->SetTitle("Counts");

 hFEB1slopes_un = new TH1F("hFEB1slopes_un", "Uncorrected Aging Slopes", 100, -25., 25.);
 hFEB1slopes_un->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB1slopes_un->GetYaxis()->SetTitle("Counts");

 hFEB2slopes_un = new TH1F("hFEB2slopes_un", "Uncorrected Aging Slopes", 100, -25., 25.);
 hFEB2slopes_un->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB2slopes_un->GetYaxis()->SetTitle("Counts");

 hFEB3slopes_un = new TH1F("hFEB3slopes_un", "Uncorrected Aging Slopes", 100, -25., 25.);
 hFEB3slopes_un->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB3slopes_un->GetYaxis()->SetTitle("Counts");

 hFEB4slopes_un = new TH1F("hFEB4slopes_un", "Uncorrected Aging Slopes", 100, -25., 25.);
 hFEB4slopes_un->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB4slopes_un->GetYaxis()->SetTitle("Counts");

 hFEB5slopes_un = new TH1F("hFEB5slopes_un", "Uncorrected Aging Slopes", 100, -25., 25.);
 hFEB5slopes_un->GetXaxis()->SetTitle("Slope (Percent of Channel Yield)");
 hFEB5slopes_un->GetYaxis()->SetTitle("Counts");

 //initialize empty hists only once, set style
 for(int ihist = 0; ihist < ctotalchan; ihist++){
   string sname;
   const char* cname;
   int ncansets = 4;
   for(int ican = 0; ican < ncansets; ican++) {
     sname = "can" + std::to_string(ihist) + "_" + std::to_string(ican);
     cname = sname.c_str();
     page[ihist][ican] = new TCanvas(cname, cname, 900, 600);
   }

   hPE_un[ihist] = new TH2F(Form("hPE_un_%i",ihist), "Uncorrected PE vs Run", 1250, 0., 1250., 200, 30., 70.);
   hPE_un[ihist]->GetXaxis()->SetTitle("Run Number");
   hPE_un[ihist]->GetYaxis()->SetTitle("PE Yield (peak value from fit)");
   hPE_un[ihist]->GetXaxis()->SetRangeUser(1050., 1250.);

   hPE_cor[ihist] = new TH2F(Form("hPE_cor_%i",ihist), "Corrected PE vs Run", 1250, 0., 1250., 200, 30., 70.);
   hPE_cor[ihist]->GetXaxis()->SetTitle("Run Number");
   hPE_cor[ihist]->GetYaxis()->SetTitle("PE Yield (corrected peak value from fit)");
   hPE_cor[ihist]->GetXaxis()->SetRangeUser(1050., 1250.); 

   hPE_un_frac[ihist] = new TH2F(Form("hPE_un_frac_%i",ihist), "Uncorrected PE vs Time Fraction", 1000, 0., 0.3, 200, 30., 70.);
   hPE_un_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 7/1/2022)");
   hPE_un_frac[ihist]->GetYaxis()->SetTitle("PE Yield (peak value from fit)"); 

   hPE_cor_frac[ihist] = new TH2F(Form("hPE_cor_frac_%i",ihist), "Corrected PE vs Time Fraction", 1000, 0., 0.3, 200, 30., 70.);
   hPE_cor_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 7/1/2022)");
   hPE_cor_frac[ihist]->GetYaxis()->SetTitle("PE Yield (corrected peak value from fit)"); 

   hPE_G_un[ihist] = new TH2F(Form("hPE_G_un_%i",ihist), "Uncorrected PE vs Run", 1250, 0., 1250., 200, 30., 70.);
   hPE_G_un[ihist]->GetXaxis()->SetTitle("Run Number");
   hPE_G_un[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean from fit)"); 
   hPE_G_un[ihist]->GetXaxis()->SetRangeUser(1050., 1250.);

   hPE_G_cor[ihist] = new TH2F(Form("hPE_G_cor_%i",ihist), "Corrected PE vs Run", 1250, 0., 1250., 200, 30., 70.);
   hPE_G_cor[ihist]->GetXaxis()->SetTitle("Run Number");
   hPE_G_cor[ihist]->GetYaxis()->SetTitle("PE Yield (corrected Gaussian mean value from fit)"); 
   hPE_G_cor[ihist]->GetXaxis()->SetRangeUser(1050., 1250.);

   hPE_G_un_frac[ihist] = new TH2F(Form("hPE_G_un_frac_%i",ihist), "Uncorrected PE vs Time Fraction", 1000, 0., 0.3, 200, 30., 70.);
   hPE_G_un_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 7/1/2022)");
   hPE_G_un_frac[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean from fit)"); 

   hPE_G_cor_frac[ihist] = new TH2F(Form("hPE_G_cor_frac_%i",ihist), "Corrected PE vs Time Fraction", 1000, 0., 0.3, 200, 30., 70.);
   hPE_G_cor_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 7/1/2022)");
   hPE_G_cor_frac[ihist]->GetYaxis()->SetTitle("PE Yield (corrected Gaussian mean value from fit)");

   hPE_chi[ihist] = new TH1F(Form("hPE_chi_%i",ihist), "Reduced Chi2 from Uncorrected PE Fits", 500, 0., 500.);
   hPE_chi[ihist]->GetXaxis()->SetTitle("Chi2 from Fit");
   hPE_chi[ihist]->GetYaxis()->SetTitle("Counts");

   hCPE_chi[ihist] = new TH1F(Form("hCPE_chi_%i",ihist), "Reduced Chi2 from Corrected PE Fits", 500, 0., 500.);
   hCPE_chi[ihist]->GetXaxis()->SetTitle("Chi2 from Fit");
   hCPE_chi[ihist]->GetYaxis()->SetTitle("Counts");
 }

 Float_t aPE_G_un[ctotalchan][100];
 Float_t aPE_G_un_err[ctotalchan][100];
 Float_t aPE_G_cor[ctotalchan][100];
 Float_t aPE_G_cor_err[ctotalchan][100];
 Float_t aPE_G_uchi[ctotalchan][100];
 Float_t aPE_G_cchi[ctotalchan][100];
 Float_t aPE_G_fullchi[ctotalchan][100];
 Float_t aempty_xerr[100];
 Float_t arun_bins[100];
 Float_t atimefrac_bins[100];
 Float_t aslopes[ctotalchan];
 Float_t aslopes_un[ctotalchan];
 Float_t aslope_frac[ctotalchan];
 Float_t aslope_frac_un[ctotalchan];
 Float_t ftimefrac = -1.;

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
    if(run_num < 1045.){
      ftimefrac = getTime(run_num, subrun_num);
    }
    if(run_num > 1045.){
      ftimefrac = getTime_168(run_num, subrun_num);
    }
   arun_bins[ifile] = frun_num;
   atimefrac_bins[ifile] = ftimefrac;
   Float_t fzero = 0.;
   aempty_xerr[ifile] = fzero;

   //declare temp histograms
   TH2F *htemp_u = new TH2F("htemp_u", "temporary_un", 1250, 0., 1250., 200, 30., 50.);
   TH2F *htemp_ug = new TH2F("htemp_ug", "temporary_un_g", 1250, 0., 1250., 200, 30., 50.);
   TH1F *htemp_uge = new TH1F("htemp_uge", "temporary_un_gerr", 1250, 0., 1250.);
   TH2F *htemp_c = new TH2F("htemp_c", "temporary_cor", 1250, 0., 1250., 200, 30., 50.);
   TH2F *htemp_cg = new TH2F("htemp_cg", "temporary_cor_g", 1250, 0., 1250., 200, 30., 50.);
   TH1F *htemp_cge = new TH1F("htemp_cge", "temporary_cor_gerr", 1250, 0., 1250.);
   TH1F *htemp_uchi = new TH1F("htemp_uchi", "temporary_uchi2", 100, 0., 10.);
   TH1F *htemp_cchi = new TH1F("htemp_cchi", "temporary_cchi2", 100, 0., 10.);
   TH1F *htemp_fullchi = new TH1F("htemp_fullchi", "temporary_fullchi2", 500, 0., 500.);

   //get all of the interesting quantities from input files
   for(int ihist = 0; ihist < ctotalchan; ihist++){
     //TH2F *htemp_u = new TH2F("htemp_u", "temporary_un", 1250, 0., 1250., 200, 30., 50.);
     htemp_u = (TH2F*)file->Get(Form("PE_vs_run_%i",ihist));
     float run_num = htemp_u->GetMean(1);
     float PE_un   = htemp_u->GetMean(2);
     int  irun_num = (int) run_num;
     htemp_u->Delete();

     //TH2F *htemp_ug = new TH2F("htemp_ug", "temporary_un_g", 1250, 0., 1250., 200, 30., 50.);
     htemp_ug = (TH2F*)file->Get(Form("gaus_vs_run_%i",ihist));
     Float_t PE_un_g = htemp_ug->GetMean(2);
     aPE_G_un[ihist][ifile] = PE_un_g;
     htemp_ug->Delete();
     
     htemp_uge = (TH1F*)file->Get(Form("gaus_vs_run_un_err_%i",ihist));
     Float_t PE_g_un_err = htemp_uge->GetBinError(irun_num);
     aPE_G_un_err[ihist][ifile] = PE_g_un_err;
     htemp_uge->Delete();

     //TH2F *htemp_c = new TH2F("htemp_c", "temporary_cor", 1250, 0., 1250., 200, 30., 50.);
     htemp_c = (TH2F*)file->Get(Form("CPE_vs_run_%i",ihist));
     float PE_cor = htemp_c->GetMean(2);
     htemp_c->Delete();

     //TH2F *htemp_cg = new TH2F("htemp_cg", "temporary_cor_g", 1250, 0., 1250., 200, 30., 50.);
     htemp_cg = (TH2F*)file->Get(Form("Cgaus_vs_run_%i",ihist));
     Float_t PE_cor_g = htemp_cg->GetMean(2);
     if(PE_cor_g < 30.){
       std::cout << "low PE yield in this channel (channel " << ihist << " in run " << irun_num << " subrun " << subrun_num << "), PE yield is: " << PE_cor_g << std::endl;
     }
     aPE_G_cor[ihist][ifile] = PE_cor_g;
     htemp_cg->Delete();

     //TH1F *htemp_cge = new TH1F("htemp_cge", "temporary_cor_gerr", 1250, 0., 1250.);
     htemp_cge = (TH1F*)file->Get(Form("gaus_vs_run_cor_err_%i",ihist));
     Float_t PE_g_cor_err = htemp_cge->GetBinError(irun_num);
     // if(PE_g_err > 100.){
     //   std::cout << "error for this channel (channel " << ihist << " in run " << run_num << " subrun " << subrun_num << ")  is: " << PE_g_err << std::endl;
     // }
     // if(PE_un == 0.){
     //   std::cout << "zero PE channel (channel " << ihist << " in run " << run_num << " subrun " << subrun_num << ")  zero: " << PE_un << " with error: " << PE_g_err << std::endl;
     // }
     aPE_G_cor_err[ihist][ifile] = PE_g_cor_err;
     htemp_cge->Delete();

     htemp_uchi = (TH1F*)file->Get(Form("PEchi_%i",ihist));
     Float_t PE_uchi = htemp_uchi->GetMean();
     aPE_G_uchi[ihist][ifile] = PE_uchi;
     htemp_uchi->Delete();

     htemp_cchi = (TH1F*)file->Get(Form("CPEchi_%i",ihist));
     Float_t PE_cchi = htemp_cchi->GetMean();
     aPE_G_cchi[ihist][ifile] = PE_cchi;
     htemp_cchi->Delete();

     htemp_fullchi = (TH1F*)file->Get(Form("PEfullchi_%i",ihist));
     Float_t PE_fullchi = htemp_fullchi->GetMean();
     aPE_G_fullchi[ihist][ifile] = PE_fullchi;
     htemp_fullchi->Delete();

     //fill histograms
     if((PE_un != 0.) && (PE_fullchi < 100.)){
       ofile->cd();
       hPE_un[ihist]->SetMarkerStyle(20);
       hPE_cor[ihist]->SetMarkerStyle(20);
       hPE_un_frac[ihist]->SetMarkerStyle(20);
       hPE_cor_frac[ihist]->SetMarkerStyle(20);
       hPE_un[ihist]->Fill(run_num,PE_un);
       hPE_cor[ihist]->Fill(run_num,PE_cor);
       hPE_un_frac[ihist]->Fill(ftimefrac,PE_un);
       hPE_cor_frac[ihist]->Fill(ftimefrac,PE_cor);
       hPE_G_un[ihist]->SetMarkerStyle(20);
       hPE_G_cor[ihist]->SetMarkerStyle(20);
       hPE_G_un_frac[ihist]->SetMarkerStyle(20);
       hPE_G_cor_frac[ihist]->SetMarkerStyle(20);
       hPE_G_un[ihist]->Fill(run_num,PE_un_g);
       hPE_G_cor[ihist]->Fill(run_num,PE_cor_g);
       hPE_G_un_frac[ihist]->Fill(ftimefrac,PE_un_g);
       hPE_G_cor_frac[ihist]->Fill(ftimefrac,PE_cor_g);
       hPE_chi[ihist]->Fill(PE_uchi);
       htotal_chi->Fill(PE_uchi);
       hCPE_chi[ihist]->Fill(PE_cchi);
       htotal_Cchi->Fill(PE_cchi);
     }
   }
   file->Close(); 
   ifile++;
 }

 Float_t PE_G_cor_avg[ifile];
 Float_t PE_G_avg_err[ifile];
 int num_full[ifile];
 //write histograms only once
 for(int ihist = 0; ihist < ctotalchan; ihist++){
   ofile->cd();
   hPE_un[ihist]->Write();
   hPE_cor[ihist]->Write();
   hPE_un_frac[ihist]->Write();
   hPE_cor_frac[ihist]->Write();
   hPE_G_un[ihist]->Write();
   hPE_G_cor[ihist]->Write();
   hPE_G_un_frac[ihist]->Fit("pol1");
   float slope_un_noerr = -100.;
   TF1 *fit_hun = (TF1*)hPE_G_un_frac[ihist]->GetListOfFunctions()->FindObject("pol1");
   if(hPE_G_un_frac[ihist]->GetMean() > 0.) slope_un_noerr = fit_hun->GetParameter(1);
   hPEhist_slopes_un->SetBinContent(ihist+1,slope_un_noerr);
   hPE_G_un_frac[ihist]->Write();
   hPE_G_cor_frac[ihist]->Fit("pol1");
   float slope_cor_noerr = -100.;
   TF1 *fit_hcor = (TF1*)hPE_G_cor_frac[ihist]->GetListOfFunctions()->FindObject("pol1");
   if(hPE_G_cor_frac[ihist]->GetMean() > 0.) slope_cor_noerr = fit_hcor->GetParameter(1);
   hPEhist_slopes->SetBinContent(ihist+1,slope_cor_noerr);
   hPE_G_cor_frac[ihist]->Write();
   hPE_chi[ihist]->Write();
   hCPE_chi[ihist]->Write();

   Float_t thisPE_G_un[ifile];
   Float_t thisPE_G_cor[ifile];
   Float_t thisPE_G_un_err[ifile];
   Float_t thisPE_G_cor_err[ifile];
   for(int i = 0; i < ifile; i++){
     if((aPE_G_un[ihist][i] == 0.) || (aPE_G_un_err[ihist][i] > 100.) || (aPE_G_fullchi[ihist][i] >= 90.47)){
       thisPE_G_un[i] = 0.;
       thisPE_G_cor[i] = 0.;
       thisPE_G_un_err[i] = 0.;
       thisPE_G_cor_err[i] = 0.;
     }
     if((aPE_G_un[ihist][i] != 0.) && (aPE_G_un_err[ihist][i] < 100.) && (aPE_G_fullchi[ihist][i] < 90.47)){
       thisPE_G_un[i] = aPE_G_un[ihist][i];
       thisPE_G_cor[i] = aPE_G_cor[ihist][i];
       thisPE_G_un_err[i] = aPE_G_un_err[ihist][i];
       thisPE_G_cor_err[i] = aPE_G_cor_err[ihist][i];
       //PE_G_cor_avg[i] += aPE_G_cor[ihist][i];
       //PE_G_avg_err[i] += aPE_G_err[ihist][i];
       num_full[i]     += 1;
     }
   }

   std::cout << "uncorrected PE array has entries: ";
   for(int i = 0; i < ifile; i++){
     std::cout << thisPE_G_un[i] << ", ";
   }
   std::cout << " end of this array" << std::endl;

   std::cout << "corrected PE array has entries: ";
   for(int i = 0; i < ifile; i++){
     std::cout << thisPE_G_cor[i] << ", ";
   }
   std::cout << " end of this array" << std::endl;

   std::cout << "uncorrected PE error array has entries: ";
   for(int i = 0; i < ifile; i++){
     std::cout << thisPE_G_un_err[i] << ", ";
   }
   std::cout << " end of this array" << std::endl;

   std::cout << "corrected PE error array has entries: ";
   for(int i = 0; i < ifile; i++){
     std::cout << thisPE_G_cor_err[i] << ", ";
   }
   std::cout << " end of this array" << std::endl;

   /* // remove zeroes from thisPE arrays, need more intermediate arrays
   Float_t thisPE_G_un0[ifile];
   Float_t thisPE_G_cor0[ifile];
   Float_t thisPE_G_un_err0[ifile];
   Float_t thisPE_G_cor_err0[ifile];
   Float_t arun_bins0[ifile];
   Float_t aempty_xerr0[ifile];
   int ifile0 = 0;
   for(int i = 0; i < ifile; i++){
     if(thisPE_G_un[i] != 0){
       thisPE_G_un0[ifile0] = thisPE_G_un[i];
       thisPE_G_cor0[ifile0] = thisPE_G_cor[i];
       thisPE_G_un_err0[ifile0] = thisPE_G_un_err[i];
       thisPE_G_cor_err0[ifile0] = thisPE_G_cor_err[i];
       ifile0 += 1;
     }
     }*/

   std::cout << "ifile = " << ifile << std::endl;
   /*page[ihist][0]->Close();
   page[ihist][1]->Close();
   page[ihist][2]->Close();
   page[ihist][3]->Close();*/
   std::cout << "ihist = " << ihist << std::endl;
   page[ihist][0]->cd();
   tgPE_G_un_err[ihist] = new TGraphErrors(ifile,arun_bins,thisPE_G_un,aempty_xerr,thisPE_G_un_err);
   tgPE_G_un_err[ihist]->SetTitle(Form("Uncorrected PE Yield vs Run (Channel %i)",ihist));
   tgPE_G_un_err[ihist]->GetXaxis()->SetTitle("Run Number");
   tgPE_G_un_err[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean)"); 
   tgPE_G_un_err[ihist]->GetYaxis()->SetRangeUser(30.,70.);
   tgPE_G_un_err[ihist]->SetMarkerStyle(20);
   tgPE_G_un_err[ihist]->Draw("AP");
   page[ihist][0]->Update();
   ofile->cd();
   page[ihist][0]->Write();
   page[ihist][0]->Close();
   page[ihist][1]->cd();
   tgPE_G_cor_err[ihist] = new TGraphErrors(ifile,arun_bins,thisPE_G_cor,aempty_xerr,thisPE_G_cor_err);
   tgPE_G_cor_err[ihist]->SetTitle(Form("Corrected PE Yield vs Run (Gaus mean, channel %i)", ihist));
   tgPE_G_cor_err[ihist]->GetXaxis()->SetTitle("Run Number");
   tgPE_G_cor_err[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean)"); 
   tgPE_G_cor_err[ihist]->GetYaxis()->SetRangeUser(30.,70.);
   tgPE_G_cor_err[ihist]->SetMarkerStyle(20);
   tgPE_G_cor_err[ihist]->Draw("AP");
   page[ihist][1]->Update();
   ofile->cd();
   gStyle->SetOptFit(1);
   page[ihist][1]->Write();
   page[ihist][1]->Close();
   page[ihist][2]->cd();
   tgPE_G_un_err_frac[ihist] = new TGraphErrors(ifile,atimefrac_bins,thisPE_G_un,aempty_xerr,thisPE_G_un_err);
   tgPE_G_un_err_frac[ihist]->SetTitle(Form("Uncorrected PE Yield vs Time Frac (Gaus mean, channel %i)",ihist));
   tgPE_G_un_err_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 7/1/2022)");
   tgPE_G_un_err_frac[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean)"); 
   //tgPE_G_un_err_frac[ihist]->GetYaxis()->SetRangeUser(20.,50.);
   tgPE_G_un_err_frac[ihist]->SetMarkerStyle(20);

   TF1 *ufunc = new TF1("ufunc", "[0]*x+[1]", -0.1, 0.31);
   ufunc->SetParameter(0, -0.5);
   ufunc->SetParameter(1, 50.);
   ufunc->SetParName(0, "m");
   ufunc->SetParName(1, "b");
   ufunc->SetParLimits(0, -10., 5.);
   ufunc->SetParLimits(1, 30., 80.);
   tgPE_G_un_err_frac[ihist]->Fit(ufunc, "R");
   float ufitp0 = ufunc->GetParameter(0);
   float ufitp1 = ufunc->GetParameter(1);
   aslopes_un[ihist] = ufitp0;
   std::cout << "slope for this channel is: " << aslopes_un[ihist] << std::endl;

   TF1* ucfunc = new TF1("ucfunc", "[0]", -0.1, 0.31);
   ucfunc->SetParameter(0, 50.);
   ucfunc->SetParLimits(0, 30., 80.);
   ucfunc->SetLineColor(kBlue);
   tgPE_G_un_err_frac[ihist]->Fit(ucfunc, "R");

   tgPE_G_un_err_frac[ihist]->GetYaxis()->SetRangeUser(30., 70.);

   tgPE_G_un_err_frac[ihist]->Draw("AP");
   page[ihist][2]->Update();
   ofile->cd();
   page[ihist][2]->Write();
   page[ihist][2]->Close();
   ufunc->Delete();
   ucfunc->Delete();
   page[ihist][3]->cd();
   tgPE_G_cor_err_frac[ihist] = new TGraphErrors(ifile,atimefrac_bins,thisPE_G_cor,aempty_xerr,thisPE_G_cor_err);
   tgPE_G_cor_err_frac[ihist]->SetTitle(Form("Corrected PE Yield vs Time Frac (Gaus mean, channel %i)",ihist));
   tgPE_G_cor_err_frac[ihist]->GetXaxis()->SetTitle("Time Fraction (years since 7/1/2022)");
   tgPE_G_cor_err_frac[ihist]->GetYaxis()->SetTitle("PE Yield (Gaussian mean)"); 
   //tgPE_G_cor_err_frac[ihist]->GetYaxis()->SetRangeUser(20.,50.);
   tgPE_G_cor_err_frac[ihist]->SetMarkerStyle(20);

   /* TF1 *newfit = new TF1("newfit","pol1");
   tgPE_G_cor_err_frac[ihist]->Fit(newfit);
   //TF1 *fit = (TF1*)tgPE_G_cor_err_frac[ihist]->GetListOfFunctions()->FindObject("pol1");
   float fitp0 = newfit->GetParameter(0);
   float fitp1 = newfit->GetParameter(1);
   aslopes[ihist] = fitp1;
   std::cout << "slope for this channel is: " << aslopes[ihist] << std::endl;
   float lowerbound = fitp0 + (fitp1*0.16);
   float upperbound = fitp0 + (fitp1*0.81);
   std::cout << "bounds for this channel are: " << lowerbound << ", " << upperbound << std::endl;
   //tgPE_G_cor_err_frac[ihist]->GetYaxis()->SetRangeUser(lowerbound-2, upperbound+2); */

   TF1 *func = new TF1("func", "[0]*x+[1]", -0.1, 0.31);
   func->SetParameter(0, -0.5);
   func->SetParameter(1, 50.);
   func->SetParName(0, "m");
   func->SetParName(1, "b");
   func->SetParLimits(0, -10., 5.);
   func->SetParLimits(1, 30., 80.);
   //TFitResultPtr r = tgPE_G_cor_err_frac[ihist].Fit(func,"S");  // I don't understand this object
   tgPE_G_cor_err_frac[ihist]->Fit(func, "R");
   float fitp0 = func->GetParameter(0);
   float fitp1 = func->GetParameter(1);
   aslopes[ihist] = fitp0;
   std::cout << "slope for this channel is: " << aslopes[ihist] << std::endl;

   TF1* cfunc = new TF1("cfunc", "[0]", -0.1, 0.31);
   cfunc->SetParameter(0, 50.);
   cfunc->SetParLimits(0, 30., 80.);
   cfunc->SetLineColor(kBlue);
   tgPE_G_cor_err_frac[ihist]->Fit(cfunc, "R");

   tgPE_G_cor_err_frac[ihist]->GetYaxis()->SetRangeUser(30., 70.);

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


   if(avgPE > 0.){
     aslope_frac[ihist] = (aslopes[ihist]/avgPE)*100;
     aslope_frac_un[ihist] = (aslopes_un[ihist]/avgPE)*100;
   }
   if(avgPE == 0.){
     aslope_frac[ihist] = -50.;
     aslope_frac_un[ihist] = -50.;
   }
   if(isnan(avgPE) > 0.){
     aslope_frac[ihist] = -50.;
     aslopes[ihist] = -50.;
     aslope_frac_un[ihist] = -50.;
     aslopes_un[ihist] = -50.;
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
   func->Delete();
   cfunc->Delete();

   hPE_slopes->SetBinContent(ihist+1,aslopes[ihist]);
   hPE_slope_frac->SetBinContent(ihist+1,aslope_frac[ihist]);
   hPE_slopes_un->SetBinContent(ihist+1,aslopes_un[ihist]);
   hPE_slope_frac_un->SetBinContent(ihist+1,aslope_frac_un[ihist]);
   hslopes1D->Fill(aslope_frac[ihist]);
   hslopes1D_un->Fill(aslope_frac_un[ihist]);
   if(ihist < 64.) {
     hFEB0slopes->Fill(aslope_frac[ihist]);
     hFEB0slopes_un->Fill(aslope_frac_un[ihist]);
   }
   if((ihist >= 64.) && (ihist < 128.)) {
     hFEB1slopes->Fill(aslope_frac[ihist]);
     hFEB1slopes_un->Fill(aslope_frac_un[ihist]);
   }
   if((ihist >= 128.) && (ihist < 192.)) {
     hFEB2slopes->Fill(aslope_frac[ihist]);
     hFEB2slopes_un->Fill(aslope_frac_un[ihist]);
   }
   if((ihist >= 192.) && (ihist < 256.)) {
     hFEB3slopes->Fill(aslope_frac[ihist]);
     hFEB3slopes_un->Fill(aslope_frac_un[ihist]);
   }
   if((ihist >= 256.) && (ihist < 320.)) {
     hFEB4slopes->Fill(aslope_frac[ihist]);
     hFEB4slopes_un->Fill(aslope_frac_un[ihist]);
   }
   if((ihist >= 320.) && (ihist < 384.)) {
     hFEB5slopes->Fill(aslope_frac[ihist]);
     hFEB5slopes_un->Fill(aslope_frac_un[ihist]);
   }
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
 tgPE_G_cor_avg_frac->GetXaxis()->SetTitle("Time Fraction (years since 7/1/2022)");
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

 slopeslabel_un[0]->cd();
 hPE_slopes_un->Draw();
 pt->Draw("same");
 slopeslabel_un[0]->Update();
 ofile->cd();
 slopeslabel_un[0]->Write();
 slopeslabel_un[0]->Close();

 slopeslabel_un[1]->cd();
 hPE_slope_frac_un->Draw();
 pt->Draw("same");
 slopeslabel_un[1]->Update();
 ofile->cd();
 slopeslabel_un[1]->Write();
 slopeslabel_un[1]->Close();

 htotal_chi->Write();
 htotal_Cchi->Write();
 hPE_slopes->Write();
 hPEhist_slopes->Write();
 hPEhist_slopes_un->Write();
 hPE_slope_frac->Write();
 hslopes1D->Write();
 hFEB0slopes->Write();
 hFEB1slopes->Write();
 hFEB2slopes->Write();
 hFEB3slopes->Write();
 hFEB4slopes->Write();
 hFEB5slopes->Write();
 hPE_slopes_un->Write();
 hPE_slope_frac_un->Write();
 hslopes1D_un->Write();
 hFEB0slopes_un->Write();
 hFEB1slopes_un->Write();
 hFEB2slopes_un->Write();
 hFEB3slopes_un->Write();
 hFEB4slopes_un->Write();
 hFEB5slopes_un->Write();
}

