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

void avgNewBias(){

  ifstream input("new127_storeflagtree_files.txt");
  string oneline;

  //Create output file. If it already exists, recreate it.
  TFile *ofile = new TFile("output_average_badfrac127_FEBbias.root","RECREATE");

  // make containers to hold average graphs
  const int nFEB = 6;
  
  TGraph *tgavg_bias_feb0[8];
  TGraph *tgavg_bias_feb1[8];
  TGraph *tgavg_bias_feb2[8];
  TGraph *tgavg_bias_feb3[8];
  TGraph *tgavg_bias_feb4[8];
  TGraph *tgavg_bias_feb5[8];
  
  TGraph *bad_fraction;
  std::vector<float> bad_fracs;
  std::vector<float> run_nums;

  //make arrays for TGraphs
  Float_t timefrac_bins[100];
  Float_t avg_bias_feb0[8][100];
  Float_t avg_bias_feb1[8][100];
  Float_t avg_bias_feb2[8][100];
  Float_t avg_bias_feb3[8][100];
  Float_t avg_bias_feb4[8][100];
  Float_t avg_bias_feb5[8][100];
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
    if((run_num < 1045.) || (run_num > 1241.)){
      ftimefrac = getTime(run_num, subrun_num);
    }
    if((run_num > 1045.) && (run_num < 1241.)){
      ftimefrac = getTime_168(run_num, subrun_num);
    }
   timefrac_bins[ifile] = ftimefrac;

   // get avg bias voltages from plots in each file

   for(int iBus = 0; iBus < 8; iBus++){
     TString busstring = std::to_string(iBus);
     TString plotstring0 = "biasBus0_" + busstring;
     TString plotstring1 = "biasBus1_" + busstring;
     TString plotstring2 = "biasBus2_" + busstring;
     TString plotstring3 = "biasBus3_" + busstring;
     TString plotstring4 = "biasBus4_" + busstring;
     TString plotstring5 = "biasBus5_" + busstring;
     TGraph *bias0 = (TGraph*)file->Get(plotstring0);
     TGraph *bias1 = (TGraph*)file->Get(plotstring1);
     TGraph *bias2 = (TGraph*)file->Get(plotstring2);
     TGraph *bias3 = (TGraph*)file->Get(plotstring3);
     TGraph *bias4 = (TGraph*)file->Get(plotstring4);
     TGraph *bias5 = (TGraph*)file->Get(plotstring5);
     avg_bias_feb0[iBus][ifile] = bias0->GetMean(2);
     avg_bias_feb1[iBus][ifile] = bias1->GetMean(2);
     avg_bias_feb2[iBus][ifile] = bias2->GetMean(2);
     avg_bias_feb3[iBus][ifile] = bias3->GetMean(2);
     avg_bias_feb4[iBus][ifile] = bias4->GetMean(2);
     avg_bias_feb5[iBus][ifile] = bias5->GetMean(2);
   }

   // get badtree with bad fraction stored
   TTree *obadtree = file->Get<TTree>("obadtree");
   float obad_frac;
   obadtree->SetBranchAddress("obad_frac",&obad_frac);
   obadtree->GetEntry(0);
   float run_float = (float)run_num;
   run_nums.push_back(run_float);
   bad_fracs.push_back(obad_frac);

   file->Close();
   ifile++;
  }

  ofile->cd();
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(11111);

  for(int iBus = 0; iBus < 8; iBus++){
    tgavg_bias_feb0[iBus] = new TGraph(ifile, timefrac_bins, avg_bias_feb0[iBus]);    
    tgavg_bias_feb0[iBus]->SetTitle(Form("avg_bias_feb0_%i", iBus));
    tgavg_bias_feb0[iBus]->SetName(Form("avg_bias_feb0_%i", iBus));
    tgavg_bias_feb0[iBus]->GetXaxis()->SetTitle("Time Fraction");
    tgavg_bias_feb0[iBus]->GetYaxis()->SetTitle("Avg Bias Voltage (V)");
    tgavg_bias_feb0[iBus]->SetMarkerStyle(20);
    tgavg_bias_feb0[iBus]->SetMarkerColor(kBlue);
    tgavg_bias_feb0[iBus]->Draw("AP");
    tgavg_bias_feb0[iBus]->Write();

    tgavg_bias_feb1[iBus] = new TGraph(ifile, timefrac_bins, avg_bias_feb1[iBus]);    
    tgavg_bias_feb1[iBus]->SetTitle(Form("avg_bias_feb1_%i", iBus));
    tgavg_bias_feb1[iBus]->SetName(Form("avg_bias_feb1_%i", iBus));
    tgavg_bias_feb1[iBus]->GetXaxis()->SetTitle("Time Fraction");
    tgavg_bias_feb1[iBus]->GetYaxis()->SetTitle("Avg Bias Voltage (V)");
    tgavg_bias_feb1[iBus]->SetMarkerStyle(20);
    tgavg_bias_feb1[iBus]->SetMarkerColor(kBlue);
    tgavg_bias_feb1[iBus]->Draw("AP");
    tgavg_bias_feb1[iBus]->Write();

    tgavg_bias_feb2[iBus] = new TGraph(ifile, timefrac_bins, avg_bias_feb2[iBus]);    
    tgavg_bias_feb2[iBus]->SetTitle(Form("avg_bias_feb2_%i", iBus));
    tgavg_bias_feb2[iBus]->SetName(Form("avg_bias_feb2_%i", iBus));
    tgavg_bias_feb2[iBus]->GetXaxis()->SetTitle("Time Fraction");
    tgavg_bias_feb2[iBus]->GetYaxis()->SetTitle("Avg Bias Voltage (V)");
    tgavg_bias_feb2[iBus]->SetMarkerStyle(20);
    tgavg_bias_feb2[iBus]->SetMarkerColor(kBlue);
    tgavg_bias_feb2[iBus]->Draw("AP");
    tgavg_bias_feb2[iBus]->Write();

    tgavg_bias_feb3[iBus] = new TGraph(ifile, timefrac_bins, avg_bias_feb3[iBus]);    
    tgavg_bias_feb3[iBus]->SetTitle(Form("avg_bias_feb3_%i", iBus));
    tgavg_bias_feb3[iBus]->SetName(Form("avg_bias_feb3_%i", iBus));
    tgavg_bias_feb3[iBus]->GetXaxis()->SetTitle("Time Fraction");
    tgavg_bias_feb3[iBus]->GetYaxis()->SetTitle("Avg Bias Voltage (V)");
    tgavg_bias_feb3[iBus]->SetMarkerStyle(20);
    tgavg_bias_feb3[iBus]->SetMarkerColor(kBlue);
    tgavg_bias_feb3[iBus]->Draw("AP");
    tgavg_bias_feb3[iBus]->Write();

    tgavg_bias_feb4[iBus] = new TGraph(ifile, timefrac_bins, avg_bias_feb4[iBus]);    
    tgavg_bias_feb4[iBus]->SetTitle(Form("avg_bias_feb4_%i", iBus));
    tgavg_bias_feb4[iBus]->SetName(Form("avg_bias_feb4_%i", iBus));
    tgavg_bias_feb4[iBus]->GetXaxis()->SetTitle("Time Fraction");
    tgavg_bias_feb4[iBus]->GetYaxis()->SetTitle("Avg Bias Voltage (V)");
    tgavg_bias_feb4[iBus]->SetMarkerStyle(20);
    tgavg_bias_feb4[iBus]->SetMarkerColor(kBlue);
    tgavg_bias_feb4[iBus]->Draw("AP");
    tgavg_bias_feb4[iBus]->Write();

    tgavg_bias_feb5[iBus] = new TGraph(ifile, timefrac_bins, avg_bias_feb5[iBus]);    
    tgavg_bias_feb5[iBus]->SetTitle(Form("avg_bias_feb5_%i", iBus));
    tgavg_bias_feb5[iBus]->SetName(Form("avg_bias_feb5_%i", iBus));
    tgavg_bias_feb5[iBus]->GetXaxis()->SetTitle("Time Fraction");
    tgavg_bias_feb5[iBus]->GetYaxis()->SetTitle("Avg Bias Voltage (V)");
    tgavg_bias_feb5[iBus]->SetMarkerStyle(20);
    tgavg_bias_feb5[iBus]->SetMarkerColor(kBlue);
    tgavg_bias_feb5[iBus]->Draw("AP");
    tgavg_bias_feb5[iBus]->Write(); 
  }

  // make plot of bad events per file
  int nruns = run_nums.size();

  float run_nums_a[nruns];
  float bad_fracs_a[nruns];
  std::copy(run_nums.begin(), run_nums.end(), run_nums_a);
  std::copy(bad_fracs.begin(), bad_fracs.end(), bad_fracs_a);

  bad_fraction = new TGraph(nruns, run_nums_a, bad_fracs_a);
  bad_fraction->SetTitle("Fraction of Bad Events per Run");
  bad_fraction->SetName("bad_fraction");
  bad_fraction->GetXaxis()->SetTitle("Run Number");
  bad_fraction->GetYaxis()->SetTitle("Fraction of Bad Events");
  bad_fraction->SetMarkerStyle(20);
  bad_fraction->Draw("AP");
  bad_fraction->Write();

  ofile->Close();
}
