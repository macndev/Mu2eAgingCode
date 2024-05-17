#include <string>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
  if(run == 1242) frac = 1.474;
  if(run == 1243) frac = 1.477;
  if(run == 1244) frac = 1.479;
  if(run == 1245) frac = 1.485;
  if(run == 1246) frac = 1.490;
  if(run == 1247) frac = 1.493;
  if(run == 1248) frac = 1.496;
  if(run == 1249) frac = 1.504;
  if(run == 1251) frac = 1.507;
  if(run == 1262) frac = 1.512;
  if(run == 1263) frac = 1.523;

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
  if(run == 1218) frac = 0.307;
  if(run == 1219) frac = 0.310;
  if(run == 1220) frac = 0.315;
  if(run == 1240) frac = 0.318;
  if(run == 1241) frac = 0.321;

  float timefrac_50hours = 0.0057;

  frac += (timefrac_50hours * subrun);

  return frac;
}

void avgNewTemps(){

  ifstream input("FEBCMBtemp_files_v3_168_foravg.txt");
  string oneline;

  //Create output file. If it already exists, recreate it.
  TFile *ofile = new TFile("output_average_FEBCMBtemp_v3_168.root","RECREATE");

  // make containers to hold average graphs

  const int nFEB  = 6;
  const int nChan = 64;
  int maxchan = nFEB * nChan;

  TGraph *ctemps_vs_time[nFEB*nChan];
  TGraph *ftemps_vs_time[nFEB];

  // make arrays for TGraphs
  Float_t timefrac_bins[100];
  Float_t avg_cmb_temps[nFEB*nChan][100];
  Float_t avg_feb_temps[nFEB][100];
  Float_t ftimefrac = -1.;

  int ifile = 0;
  while(std::getline(input, oneline)){
    const char* cfile = oneline.c_str();
    TFile *file = TFile::Open(cfile);
    std::cout << "file " << cfile << " is opening..." << std::endl;

    TString filename;
    filename = file->GetName();
    int run_num;
    TString short_runstr = filename(filename.Length()-15, filename.Length()-9);
    run_num = short_runstr.Atoi();
    int subrun_num;
    TString short_subrunstr = filename(filename.Length()-8, filename.Length()-5);
    subrun_num = short_subrunstr.Atoi();
    std::cout << "working on file: " << run_num << "_" << subrun_num << std::endl;

   // fill x bins with time fractions
    if(run_num < 1045.){
      ftimefrac = getTime(run_num, subrun_num);
    }
    if(run_num > 1045.){
      ftimefrac = getTime_168(run_num, subrun_num);
    }
    timefrac_bins[ifile] = ftimefrac;

   // get average temps from plots in each file

   for(int iFEB = 0; iFEB < nFEB; iFEB++){
     TString idstring = std::to_string(iFEB);
     TString febstring = "FEBtemp_" + idstring;
     //std::cout << "febstring is: " << febstring << std::endl;
     TGraph *FEB = (TGraph*)file->Get(febstring);
     avg_feb_temps[iFEB][ifile] = FEB->GetMean(2);
   }

   for(int iCMB = 0; iCMB < maxchan; iCMB++){
     TString idstring = std::to_string(iCMB);
     TString cmbstring = "CMBtemp_" + idstring;
     //std::cout << "cmbstring is: " << cmbstring << std::endl;
     TGraph *CMB = (TGraph*)file->Get(cmbstring);
     avg_cmb_temps[iCMB][ifile] = CMB->GetMean(2);
   }
   file->Close();
   ifile++;
  }

  ofile->cd();
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(11111);

  std::cout << "average FEB0 temp array: ";
  for(int i = 0; i < 80; i++){
    std::cout << avg_feb_temps[0][i] << ", ";
  }
  std::cout << " end" << std::endl;
  std::cout << "average CMB0 temp array: ";
  for(int i = 0; i < 80; i++){
    std::cout << avg_cmb_temps[0][i] << ", ";
  }
  std::cout << " end" << std::endl;

  for(int iFEB = 0; iFEB < nFEB; iFEB++){
    ftemps_vs_time[iFEB] = new TGraph(ifile, timefrac_bins, avg_feb_temps[iFEB]);
    ftemps_vs_time[iFEB]->SetTitle(Form("avgFEBtemp_%i", iFEB));
    ftemps_vs_time[iFEB]->SetName(Form("avgFEBtemp_%i", iFEB));
    ftemps_vs_time[iFEB]->GetXaxis()->SetTitle("Time Fraction");
    ftemps_vs_time[iFEB]->GetYaxis()->SetTitle("Avg FEB Temperature (C)");
    ftemps_vs_time[iFEB]->SetMarkerColor(kRed);
    ftemps_vs_time[iFEB]->SetMarkerStyle(20);
    ftemps_vs_time[iFEB]->SetLineColor(kRed);
    ftemps_vs_time[iFEB]->Draw("AP");
    ftemps_vs_time[iFEB]->Write();
  }
  for(int iCMB = 0; iCMB < maxchan; iCMB++){
    ctemps_vs_time[iCMB] = new TGraph(ifile, timefrac_bins, avg_cmb_temps[iCMB]);
    ctemps_vs_time[iCMB]->SetTitle(Form("avgCMBtemp_%i", iCMB));
    ctemps_vs_time[iCMB]->SetName(Form("avgCMBtemp_%i", iCMB));
    ctemps_vs_time[iCMB]->GetXaxis()->SetTitle("Time Fraction");
    ctemps_vs_time[iCMB]->GetYaxis()->SetTitle("Avg CMB Temperature (C)");
    ctemps_vs_time[iCMB]->SetMarkerColor(kBlue);
    ctemps_vs_time[iCMB]->SetMarkerStyle(20);
    ctemps_vs_time[iCMB]->SetLineColor(kBlue);
    ctemps_vs_time[iCMB]->Draw("AP");
    ctemps_vs_time[iCMB]->Write();
  }
  ofile->Close();
}
