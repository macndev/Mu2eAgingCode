// this file takes the reflectance data histograms as an input and will do some
// fancier analysis like making averages for a particular sample, taking 
// differences between two samples, and plotting data over time

TH1F* difference_hist(TH1F* hist1, TH1F* hist2, const char* name){
  int   maxnbins      = 43;
  float diff_bin      = 0.;
  int   nbins         = 0;

  TH1F* diff_hist = new TH1F(name, name, 43, 360., 780.);
  diff_hist->GetXaxis()->SetTitle("Wavelength (nm)");
  diff_hist->GetYaxis()->SetTitle("Reflectance (R) (%)");
  
  for(int bin = 1; bin <= maxnbins;bin++){
    float val_1 = hist1->GetBinContent(bin);
    float val_2 = hist2->GetBinContent(bin);
    if((val_1 != 0) && (val_2 != 0)){
      diff_bin = val_1 - val_2;
    }
    else {
      diff_bin = 0;
    }
    diff_hist->SetBinContent(bin,diff_bin);
  }
  return diff_hist;
}

TH1F* average_hist(std::vector<TH1F*> hists, std::vector<int> indices, const char* name){
  int   maxnbins     = 43;
  float avg_bins[43] = {0};

  TH1F* avg_hist  = new TH1F(name, name, 43, 360., 780.);
  avg_hist->GetXaxis()->SetTitle("Wavelength (nm)");
  avg_hist->GetYaxis()->SetTitle("Reflectance (R) (%)");

  for(int i = 0;i < indices.size();i++){
    int this_ind = indices[i];
    int nbins    = hists[this_ind]->GetNbinsX();
    for(int bin = 1;bin <= nbins;bin++){
      avg_bins[bin] += hists[this_ind]->GetBinContent(bin);
    }
  }
  for(int ib = 1;ib <= maxnbins;ib++){
    avg_bins[ib] /= indices.size();
    avg_hist->SetBinContent(ib,avg_bins[ib]);
  }
  return avg_hist;
  
}

float stdev_error(std::vector<TH1F*> hists, std::vector<int> indices, int nm){
  std::vector<float> these_points;
  int bin;

  if(nm == 380) bin = 3;
  if(nm == 390) bin = 4;
  if(nm == 400) bin = 5;
  if(nm == 410) bin = 6;
  if(nm == 420) bin = 7;
  if(nm == 430) bin = 8;
  if(nm == 440) bin = 9;
  if(nm == 450) bin = 10;
  if(nm == 460) bin = 11;
  if(nm == 470) bin = 12;
  if(nm == 480) bin = 13;
  if(nm == 490) bin = 14;
  if(nm == 500) bin = 15;

  for(int ind = 0;ind < indices.size();ind++){
    int hist = indices[ind];
    these_points.push_back(hists[hist]->GetBinContent(bin)); // get relevant wavelength bin and store 
  }
  float stdev = TMath::StdDev(these_points.begin(), these_points.end());
  return stdev;
}

void reflectance_analysis_aging(const char* filename){

  TFile* file = TFile::Open(Form("%s", filename));
  if (file->IsZombie()) {
    std::cout << "Problem opening file " << filename << std::endl;
    return;
  }

  //Create output file. If it already exists, recreate it.
  TFile *ofile = new TFile("output_TiO2_analysis_aging_1010.root","RECREATE");

  // define average and difference histograms for each sample
  TH1F* USPN085939_thin_mackenzie_avg  = new TH1F("USPN085939_thin_mackenzie_avg", "USPN085939_thin_mackenzie_avg",43,360.,780.);
  TH1F* USPN085939_thick_mackenzie_avg = new TH1F("USPN085939_thick_mackenzie_avg", "USPN085939_thick_mackenzie_avg",43,360.,780.);
  TH1F* USPN085939_thin_brian_avg      = new TH1F("USPN085939_thin_brian_avg", "USPN085939_thin_brian_avg",43,360.,780.);
  TH1F* USPN085939_thick_brian_avg     = new TH1F("USPN085939_thick_brian_avg", "USPN085939_thick_brian_avg",43,360.,780.);
  TH1F* USPN085939_thin_diff           = new TH1F("USPN085939_thin_diff_b-m", "USPN085939_thin_diff_b-m",43,360.,780.);
  TH1F* USPN085939_thick_diff          = new TH1F("USPN085939_thick_diff_b-m", "USPN085939_thick_diff_b-m",43,360.,780.);

  TH1F* USPN024188_thin_mackenzie_avg  = new TH1F("USPN024188_thin_mackenzie_avg", "USPN024188_thin_mackenzie_avg",43,360.,780.);
  TH1F* USPN024188_thick_mackenzie_avg = new TH1F("USPN024188_thick_mackenzie_avg", "USPN024188_thick_mackenzie_avg",43,360.,780.);
  TH1F* USPN024188_thin_brian_avg      = new TH1F("USPN024188_thin_brian_avg", "USPN024188_thin_brian_avg",43,360.,780.);
  TH1F* USPN024188_thick_brian_avg     = new TH1F("USPN024188_thick_brian_avg", "USPN024188_thick_brian_avg",43,360.,780.);
  TH1F* USPN024188_thin_diff           = new TH1F("USPN024188_thin_diff_b-m", "USPN024188_thin_diff_b-m",43,360.,780.);
  TH1F* USPN024188_thick_diff          = new TH1F("USPN024188_thick_diff_b-m", "USPN024188_thick_diff_b-m",43,360.,780.);

  TH1F* USPN035597_thin_mackenzie_avg  = new TH1F("USPN035597_thin_mackenzie_avg", "USPN035597_thin_mackenzie_avg",43,360.,780.);
  TH1F* USPN035597_thick_mackenzie_avg = new TH1F("USPN035597_thick_mackenzie_avg", "USPN035597_thick_mackenzie_avg",43,360.,780.);
  TH1F* USPN035597_thin_brian_avg      = new TH1F("USPN035597_thin_brian_avg", "USPN035597_thin_brian_avg",43,360.,780.);
  TH1F* USPN035597_thick_brian_avg     = new TH1F("USPN035597_thick_brian_avg", "USPN035597_thick_brian_avg",43,360.,780.);
  TH1F* USPN035597_thin_diff           = new TH1F("USPN035597_thin_diff_b-m", "USPN035597_thin_diff_b-m",43,360.,780.);
  TH1F* USPN035597_thick_diff          = new TH1F("USPN035597_thick_diff_b-m", "USPN035597_thick_diff_b-m",43,360.,780.);

  TH1F* USPN100423_thin_mackenzie_avg  = new TH1F("USPN100423_thin_mackenzie_avg", "USPN100423_thin_mackenzie_avg",43,360.,780.);
  TH1F* USPN100423_thick_mackenzie_avg = new TH1F("USPN100423_thick_mackenzie_avg", "USPN100423_thick_mackenzie_avg",43,360.,780.);
  TH1F* USPN100423_thin_brian_avg      = new TH1F("USPN100423_thin_brian_avg", "USPN100423_thin_brian_avg",43,360.,780.);
  TH1F* USPN100423_thick_brian_avg     = new TH1F("USPN100423_thick_brian_avg", "USPN100423_thick_brian_avg",43,360.,780.);
  TH1F* USPN100423_thin_diff           = new TH1F("USPN100423_thin_diff_b-m", "USPN100423_thin_diff_b-m",43,360.,780.);
  TH1F* USPN100423_thick_diff          = new TH1F("USPN100423_thick_diff_b-m", "USPN100423_thick_diff_b-m",43,360.,780.);

  TH1F* USPN015794_thin_mackenzie_avg  = new TH1F("USPN015794_thin_mackenzie_avg", "USPN015794_thin_mackenzie_avg",43,360.,780.);
  TH1F* USPN015794_thick_mackenzie_avg = new TH1F("USPN015794_thick_mackenzie_avg", "USPN015794_thick_mackenzie_avg",43,360.,780.);
  TH1F* USPN015794_thin_brian_avg      = new TH1F("USPN015794_thin_brian_avg", "USPN015794_thin_brian_avg",43,360.,780.);
  TH1F* USPN015794_thick_brian_avg     = new TH1F("USPN015794_thick_brian_avg", "USPN015794_thick_brian_avg",43,360.,780.);
  TH1F* USPN015794_thin_diff           = new TH1F("USPN015794_thin_diff_b-m", "USPN015794_thin_diff_b-m",43,360.,780.);
  TH1F* USPN015794_thick_diff          = new TH1F("USPN015794_thick_diff_b-m", "USPN015794_thick_diff_b-m",43,360.,780.);

  TH1F* USPN014804_thin_mackenzie_avg  = new TH1F("USPN014804_thin_mackenzie_avg", "USPN014804_thin_mackenzie_avg",43,360.,780.);
  TH1F* USPN014804_thick_mackenzie_avg = new TH1F("USPN014804_thick_mackenzie_avg", "USPN014804_thick_mackenzie_avg",43,360.,780.);

  TH1F* KE08508_thin_mackenzie_avg  = new TH1F("KE08508_thin_mackenzie_avg", "KE08508_thin_mackenzie_avg",43,360.,780.);
  TH1F* KE08508_thick_mackenzie_avg = new TH1F("KE08508_thick_mackenzie_avg", "KE08508_thick_mackenzie_avg",43,360.,780.);

  TH1F* KE05173_thin_mackenzie_avg  = new TH1F("KE05173_thin_mackenzie_avg", "KE05173_thin_mackenzie_avg",43,360.,780.);
  TH1F* KE05173_mid_mackenzie_avg   = new TH1F("KE05173_mid_mackenzie_avg", "KE05173_mid_mackenzie_avg",43,360.,780.);
  TH1F* KE05173_thick_mackenzie_avg = new TH1F("KE05173_thick_mackenzie_avg", "KE05173_thick_mackenzie_avg",43,360.,780.);

  TH1F* WHC26311A_thin_mackenzie_avg  = new TH1F("WHC26311A_thin_mackenzie_avg", "WHC26311A_thin_mackenzie_avg",43,360.,780.);
  TH1F* WHC26311A_mid_mackenzie_avg   = new TH1F("WHC26311A_mid_mackenzie_avg", "WHC26311A_mid_mackenzie_avg",43,360.,780.);
  TH1F* WHC26311A_thick_mackenzie_avg = new TH1F("WHC26311A_thick_mackenzie_avg", "WHC26311A_thick_mackenzie_avg",43,360.,780.);

  TH1F* KE19161_thin_brian_avg      = new TH1F("KE19161_thin_brian_avg", "KE19161_thin_brian_avg",43,360.,780.);
  TH1F* KE19161_thick_brian_avg     = new TH1F("KE19161_thick_brian_avg", "KE19161_thick_brian_avg",43,360.,780.);
  TH1F* KE19161_thin_mackenzie_avg  = new TH1F("KE19161_thin_mackenzie_avg", "KE19161_thin_mackenzie_avg",43,360.,780.);
  TH1F* KE19161_thick_mackenzie_avg = new TH1F("KE19161_thick_mackenzie_avg", "KE19161_thick_mackenzie_avg",43,360.,780.);
  TH1F* KE19161_thin_diff           = new TH1F("KE19161_thin_diff_b-m", "KE19161_thin_diff_b-m",43,360.,780.);
  TH1F* KE19161_thick_diff          = new TH1F("KE19161_thick_diff_b-m", "KE19161_thick_diff_b-m",43,360.,780.);

  TH1F* KE12291_thin_brian_avg      = new TH1F("KE12291_thin_brian_avg", "KE12291_thin_brian_avg",43,360.,780.);
  TH1F* KE12291_thick_brian_avg     = new TH1F("KE12291_thick_brian_avg", "KE12291_thick_brian_avg",43,360.,780.);
  TH1F* KE12291_thin_mackenzie_avg  = new TH1F("KE12291_thin_mackenzie_avg", "KE12291_thin_mackenzie_avg",43,360.,780.);
  TH1F* KE12291_mid_mackenzie_avg   = new TH1F("KE12291_mid_mackenzie_avg", "KE12291_mid_mackenzie_avg",43,360.,780.);
  TH1F* KE12291_thick_mackenzie_avg = new TH1F("KE12291_thick_mackenzie_avg", "KE12291_thick_mackenzie_avg",43,360.,780.);
  TH1F* KE12291_thin_diff           = new TH1F("KE12291_thin_diff_b-m", "KE12291_thin_diff_b-m",43,360.,780.);
  TH1F* KE12291_thick_diff          = new TH1F("KE12291_thick_diff_b-m", "KE12291_thick_diff_b-m",43,360.,780.);

  TH1F* KE23167_thin_brian_avg      = new TH1F("KE23167_thin_brian_avg", "KE23167_thin_brian_avg",43,360.,780.);
  TH1F* KE23167_thick_brian_avg     = new TH1F("KE23167_thick_brian_avg", "KE23167_thick_brian_avg",43,360.,780.);
  TH1F* KE23167_thin_mackenzie_avg  = new TH1F("KE23167_thin_mackenzie_avg", "KE23167_thin_mackenzie_avg",43,360.,780.);
  TH1F* KE23167_thick_mackenzie_avg = new TH1F("KE23167_thick_mackenzie_avg", "KE23167_thick_mackenzie_avg",43,360.,780.);
  TH1F* KE23167_thin_diff           = new TH1F("KE23167_thin_diff_b-m", "KE23167_thin_diff_b-m",43,360.,780.);
  TH1F* KE23167_thick_diff          = new TH1F("KE23167_thick_diff_b-m", "KE23167_thick_diff_b-m",43,360.,780.);

  TH1F* JE05501_thin_mackenzie_avg  = new TH1F("JE05501_thin_mackenzie_avg", "JE05501_thin_mackenzie_avg",43,360.,780.);
  TH1F* JE05501_thick_mackenzie_avg = new TH1F("JE05501_thick_mackenzie_avg", "JE05501_thick_mackenzie_avg",43,360.,780.);

  TH1F* KE05130_thin_mackenzie_avg  = new TH1F("KE05130_thin_mackenzie_avg", "KE05130_thin_mackenzie_avg",43,360.,780.);
  TH1F* KE05130_thick_mackenzie_avg = new TH1F("KE05130_thick_mackenzie_avg", "KE05130_thick_mackenzie_avg",43,360.,780.);

  TH1F* IE08144_thin_mackenzie_avg  = new TH1F("IE08144_thin_mackenzie_avg", "IE08144_thin_mackenzie_avg",43,360.,780.);
  TH1F* IE08144_thick_mackenzie_avg = new TH1F("IE08144_thick_mackenzie_avg", "IE08144_thick_mackenzie_avg",43,360.,780.);

  TH1F* IE08146_thin_mackenzie_avg  = new TH1F("IE08146_thin_mackenzie_avg", "IE08146_thin_mackenzie_avg",43,360.,780.);
  TH1F* IE08146_thick_mackenzie_avg = new TH1F("IE08146_thick_mackenzie_avg", "IE08146_thick_mackenzie_avg",43,360.,780.);

  TH1F* IE00108_thin_mackenzie_avg  = new TH1F("IE00108_thin_mackenzie_avg", "IE00108_thin_mackenzie_avg",43,360.,780.);
  TH1F* IE00108_mid_mackenzie_avg   = new TH1F("IE00108_mid_mackenzie_avg", "IE00108_mid_mackenzie_avg",43,360.,780.);
  TH1F* IE00108_thick_mackenzie_avg = new TH1F("IE00108_thick_mackenzie_avg", "IE00108_thick_mackenzie_avg",43,360.,780.);

  TH1F* IE00109_thin_mackenzie_avg  = new TH1F("IE00109_thin_mackenzie_avg", "IE00109_thin_mackenzie_avg",43,360.,780.);
  TH1F* IE00109_mid_mackenzie_avg   = new TH1F("IE00109_mid_mackenzie_avg", "IE00109_mid_mackenzie_avg",43,360.,780.);
  TH1F* IE00109_thick_mackenzie_avg = new TH1F("IE00109_thick_mackenzie_avg", "IE00109_thick_mackenzie_avg",43,360.,780.);

  TH1F* USPN115101_thin_mackenzie_317_avg  = new TH1F("USPN115101_thin_mackenzie_317_avg", "USPN115101_thin_mackenzie_317_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_317_avg = new TH1F("USPN115101_thick_mackenzie_317_avg", "USPN115101_thick_mackenzie_317_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_324_avg  = new TH1F("USPN115101_thin_mackenzie_324_avg", "USPN115101_thin_mackenzie_324_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_324_avg = new TH1F("USPN115101_thick_mackenzie_324_avg", "USPN115101_thick_mackenzie_324_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_407_avg  = new TH1F("USPN115101_thin_mackenzie_407_avg", "USPN115101_thin_mackenzie_407_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_407_avg = new TH1F("USPN115101_thick_mackenzie_407_avg", "USPN115101_thick_mackenzie_407_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_414_avg  = new TH1F("USPN115101_thin_mackenzie_414_avg", "USPN115101_thin_mackenzie_414_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_414_avg = new TH1F("USPN115101_thick_mackenzie_414_avg", "USPN115101_thick_mackenzie_414_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_421_avg  = new TH1F("USPN115101_thin_mackenzie_421_avg", "USPN115101_thin_mackenzie_421_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_421_avg = new TH1F("USPN115101_thick_mackenzie_421_avg", "USPN115101_thick_mackenzie_421_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_428_avg  = new TH1F("USPN115101_thin_mackenzie_428_avg", "USPN115101_thin_mackenzie_428_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_428_avg = new TH1F("USPN115101_thick_mackenzie_428_avg", "USPN115101_thick_mackenzie_428_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_505_avg  = new TH1F("USPN115101_thin_mackenzie_505_avg", "USPN115101_thin_mackenzie_505_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_505_avg = new TH1F("USPN115101_thick_mackenzie_505_avg", "USPN115101_thick_mackenzie_505_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_512_avg  = new TH1F("USPN115101_thin_mackenzie_512_avg", "USPN115101_thin_mackenzie_512_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_512_avg = new TH1F("USPN115101_thick_mackenzie_512_avg", "USPN115101_thick_mackenzie_512_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_526_avg  = new TH1F("USPN115101_thin_mackenzie_526_avg", "USPN115101_thin_mackenzie_526_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_526_avg = new TH1F("USPN115101_thick_mackenzie_526_avg", "USPN115101_thick_mackenzie_526_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_602_avg  = new TH1F("USPN115101_thin_mackenzie_602_avg", "USPN115101_thin_mackenzie_602_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_602_avg = new TH1F("USPN115101_thick_mackenzie_602_avg", "USPN115101_thick_mackenzie_602_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_616_avg  = new TH1F("USPN115101_thin_mackenzie_616_avg", "USPN115101_thin_mackenzie_616_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_616_avg = new TH1F("USPN115101_thick_mackenzie_616_avg", "USPN115101_thick_mackenzie_616_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_714_avg  = new TH1F("USPN115101_thin_mackenzie_714_avg", "USPN115101_thin_mackenzie_714_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_714_avg = new TH1F("USPN115101_thick_mackenzie_714_avg", "USPN115101_thick_mackenzie_714_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_728_avg  = new TH1F("USPN115101_thin_mackenzie_728_avg", "USPN115101_thin_mackenzie_728_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_728_avg = new TH1F("USPN115101_thick_mackenzie_728_avg", "USPN115101_thick_mackenzie_728_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_811_avg  = new TH1F("USPN115101_thin_mackenzie_811_avg", "USPN115101_thin_mackenzie_811_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_811_avg = new TH1F("USPN115101_thick_mackenzie_811_avg", "USPN115101_thick_mackenzie_811_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_825_avg  = new TH1F("USPN115101_thin_mackenzie_825_avg", "USPN115101_thin_mackenzie_825_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_825_avg = new TH1F("USPN115101_thick_mackenzie_825_avg", "USPN115101_thick_mackenzie_825_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_908_avg  = new TH1F("USPN115101_thin_mackenzie_908_avg", "USPN115101_thin_mackenzie_908_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_908_avg = new TH1F("USPN115101_thick_mackenzie_908_avg", "USPN115101_thick_mackenzie_908_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_920_avg  = new TH1F("USPN115101_thin_mackenzie_920_avg", "USPN115101_thin_mackenzie_920_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_920_avg = new TH1F("USPN115101_thick_mackenzie_920_avg", "USPN115101_thick_mackenzie_920_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_1004_avg  = new TH1F("USPN115101_thin_mackenzie_1004_avg", "USPN115101_thin_mackenzie_1004_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_1004_avg = new TH1F("USPN115101_thick_mackenzie_1004_avg", "USPN115101_thick_mackenzie_1004_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_1018_avg  = new TH1F("USPN115101_thin_mackenzie_1018_avg", "USPN115101_thin_mackenzie_1018_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_1018_avg = new TH1F("USPN115101_thick_mackenzie_1018_avg", "USPN115101_thick_mackenzie_1018_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_1101_avg  = new TH1F("USPN115101_thin_mackenzie_1101_avg", "USPN115101_thin_mackenzie_1101_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_1101_avg = new TH1F("USPN115101_thick_mackenzie_1101_avg", "USPN115101_thick_mackenzie_1101_avg",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_1115_avg  = new TH1F("USPN115101_thin_mackenzie_1115_avg", "USPN115101_thin_mackenzie_1115_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_1115_avg = new TH1F("USPN115101_thick_mackenzie_1115_avg", "USPN115101_thick_mackenzie_1115_avg",43,360.,780.);
  TH1F* USPN115101_mackenzie_1115_diff = new TH1F("USPN115101_mackenzie_1115_diff", "USPN115101_mackenzie_1115_diff",43,360.,780.);
  TH1F* USPN115101_thin_mackenzie_1129_avg  = new TH1F("USPN115101_thin_mackenzie_1129_avg", "USPN115101_thin_mackenzie_1129_avg",43,360.,780.);
  TH1F* USPN115101_thick_mackenzie_1129_avg = new TH1F("USPN115101_thick_mackenzie_1129_avg", "USPN115101_thick_mackenzie_1129_avg",43,360.,780.);

  TH1F* USPN115101_thin_brian_020123_avg  = new TH1F("USPN115101_thin_brian_020123_avg", "USPN115101_thin_brian_020123_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_020123_avg = new TH1F("USPN115101_thick_brian_020123_avg", "USPN115101_thick_brian_020123_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_020923_avg  = new TH1F("USPN115101_thin_brian_020923_avg", "USPN115101_thin_brian_020923_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_020923_avg = new TH1F("USPN115101_thick_brian_020923_avg", "USPN115101_thick_brian_020923_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_021723_avg  = new TH1F("USPN115101_thin_brian_021723_avg", "USPN115101_thin_brian_021723_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_021723_avg = new TH1F("USPN115101_thick_brian_021723_avg", "USPN115101_thick_brian_021723_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_030923_avg  = new TH1F("USPN115101_thin_brian_030923_avg", "USPN115101_thin_brian_030923_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_030923_avg = new TH1F("USPN115101_thick_brian_030923_avg", "USPN115101_thick_brian_030923_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_031623_avg  = new TH1F("USPN115101_thin_brian_031623_avg", "USPN115101_thin_brian_031623_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_031623_avg = new TH1F("USPN115101_thick_brian_031623_avg", "USPN115101_thick_brian_031623_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_032323_avg  = new TH1F("USPN115101_thin_brian_032323_avg", "USPN115101_thin_brian_032323_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_032323_avg = new TH1F("USPN115101_thick_brian_032323_avg", "USPN115101_thick_brian_032323_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_040623_avg  = new TH1F("USPN115101_thin_brian_040623_avg", "USPN115101_thin_brian_040623_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_040623_avg = new TH1F("USPN115101_thick_brian_040623_avg", "USPN115101_thick_brian_040623_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_041323_avg  = new TH1F("USPN115101_thin_brian_041323_avg", "USPN115101_thin_brian_041323_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_041323_avg = new TH1F("USPN115101_thick_brian_041323_avg", "USPN115101_thick_brian_041323_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_042123_avg  = new TH1F("USPN115101_thin_brian_042123_avg", "USPN115101_thin_brian_042123_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_042123_avg = new TH1F("USPN115101_thick_brian_042123_avg", "USPN115101_thick_brian_042123_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_052423_avg  = new TH1F("USPN115101_thin_brian_052423_avg", "USPN115101_thin_brian_052423_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_052423_avg = new TH1F("USPN115101_thick_brian_052423_avg", "USPN115101_thick_brian_052423_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_060823_avg  = new TH1F("USPN115101_thin_brian_060823_avg", "USPN115101_thin_brian_060823_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_060823_avg = new TH1F("USPN115101_thick_brian_060823_avg", "USPN115101_thick_brian_060823_avg",43,360.,780.);
  TH1F* USPN115101_thin_brian_062323_avg  = new TH1F("USPN115101_thin_brian_062323_avg", "USPN115101_thin_brian_062323_avg",43,360.,780.);
  TH1F* USPN115101_thick_brian_062323_avg = new TH1F("USPN115101_thick_brian_062323_avg", "USPN115101_thick_brian_062323_avg",43,360.,780.);


  // get list of keys in old file
  TIter keyList(file->GetListOfKeys());
  TKey* key;
  std::vector<TString> keyNames;
  std::vector<string> key_strs;
  std::vector<TH1F*> input_hists;
  while((key = (TKey*)keyList())){
    // std::cout << key->GetName() << " " << key->GetClassName() << std::endl;
    TString key_name = key->GetName();
    keyNames.push_back(key_name);
    input_hists.push_back(new TH1F(key_name, key_name,43,360.,780.));
  }

  // get histograms from input file
  int nkeys = keyNames.size();
  for(int i = 0;i < nkeys;i++){
    //std::cout << keyNames[i] << std::endl;
    input_hists[i] = (TH1F*)file->Get(keyNames[i]);
    std::string key_string(keyNames[i].Data());
    key_strs.push_back(key_string);
  }

  // determine groups of files that are interesting...
  // here, you want to collect vectors of indices into the key vectors to do operations on a particular group of hists
  std::vector<int> avg_USPN085939_brian_thin;
  std::vector<int> avg_USPN085939_brian_thick;
  std::vector<int> avg_USPN085939_mackenzie_thin;
  std::vector<int> avg_USPN085939_mackenzie_thick;

  std::vector<int> avg_USPN024188_brian_thin;
  std::vector<int> avg_USPN024188_brian_thick;
  std::vector<int> avg_USPN024188_mackenzie_thin;
  std::vector<int> avg_USPN024188_mackenzie_thick;

  std::vector<int> avg_USPN035597_brian_thin;
  std::vector<int> avg_USPN035597_brian_thick;
  std::vector<int> avg_USPN035597_mackenzie_thin;
  std::vector<int> avg_USPN035597_mackenzie_thick;

  std::vector<int> avg_USPN100423_brian_thin;
  std::vector<int> avg_USPN100423_brian_thick;
  std::vector<int> avg_USPN100423_mackenzie_thin;
  std::vector<int> avg_USPN100423_mackenzie_thick;

  std::vector<int> avg_USPN015794_brian_thin;
  std::vector<int> avg_USPN015794_brian_thick;
  std::vector<int> avg_USPN015794_mackenzie_thin;
  std::vector<int> avg_USPN015794_mackenzie_thick;

  std::vector<int> avg_USPN014804_mackenzie_thin;
  std::vector<int> avg_USPN014804_mackenzie_thick;

  std::vector<int> avg_KE08508_mackenzie_thin;
  std::vector<int> avg_KE08508_mackenzie_thick;

  std::vector<int> avg_KE05173_mackenzie_thin;
  std::vector<int> avg_KE05173_mackenzie_mid;
  std::vector<int> avg_KE05173_mackenzie_thick;

  std::vector<int> avg_WHC26311A_mackenzie_thin;
  std::vector<int> avg_WHC26311A_mackenzie_mid;
  std::vector<int> avg_WHC26311A_mackenzie_thick;

  std::vector<int> avg_KE19161_mackenzie_thin;
  std::vector<int> avg_KE19161_mackenzie_thick;
  std::vector<int> avg_KE19161_brian_thin;
  std::vector<int> avg_KE19161_brian_thick;

  std::vector<int> avg_KE12291_mackenzie_thin;
  std::vector<int> avg_KE12291_mackenzie_mid;
  std::vector<int> avg_KE12291_mackenzie_thick;
  std::vector<int> avg_KE12291_brian_thin;
  std::vector<int> avg_KE12291_brian_thick;

  std::vector<int> avg_KE23167_mackenzie_thin;
  std::vector<int> avg_KE23167_mackenzie_thick;
  std::vector<int> avg_KE23167_brian_thin;
  std::vector<int> avg_KE23167_brian_thick;

  std::vector<int> avg_JE05501_mackenzie_thin;
  std::vector<int> avg_JE05501_mackenzie_thick;

  std::vector<int> avg_KE05130_mackenzie_thin;
  std::vector<int> avg_KE05130_mackenzie_thick;

  std::vector<int> avg_IE08144_mackenzie_thin;
  std::vector<int> avg_IE08144_mackenzie_thick;

  std::vector<int> avg_IE08146_mackenzie_thin;
  std::vector<int> avg_IE08146_mackenzie_thick;

  std::vector<int> avg_IE00108_mackenzie_thin;
  std::vector<int> avg_IE00108_mackenzie_mid;
  std::vector<int> avg_IE00108_mackenzie_thick;

  std::vector<int> avg_IE00109_mackenzie_thin;
  std::vector<int> avg_IE00109_mackenzie_mid;
  std::vector<int> avg_IE00109_mackenzie_thick;

  std::vector<int> avg_USPN115101_mackenzie_thin_317;
  std::vector<int> avg_USPN115101_mackenzie_thick_317;
  std::vector<int> avg_USPN115101_mackenzie_thin_324;
  std::vector<int> avg_USPN115101_mackenzie_thick_324;
  std::vector<int> avg_USPN115101_mackenzie_thin_407;
  std::vector<int> avg_USPN115101_mackenzie_thick_407;
  std::vector<int> avg_USPN115101_mackenzie_thin_414;
  std::vector<int> avg_USPN115101_mackenzie_thick_414;
  std::vector<int> avg_USPN115101_mackenzie_thin_421;
  std::vector<int> avg_USPN115101_mackenzie_thick_421;
  std::vector<int> avg_USPN115101_mackenzie_thin_428;
  std::vector<int> avg_USPN115101_mackenzie_thick_428;
  std::vector<int> avg_USPN115101_mackenzie_thin_505;
  std::vector<int> avg_USPN115101_mackenzie_thick_505;  
  std::vector<int> avg_USPN115101_mackenzie_thin_512;
  std::vector<int> avg_USPN115101_mackenzie_thick_512;
  std::vector<int> avg_USPN115101_mackenzie_thin_526;
  std::vector<int> avg_USPN115101_mackenzie_thick_526;
  std::vector<int> avg_USPN115101_mackenzie_thin_602;
  std::vector<int> avg_USPN115101_mackenzie_thick_602;
  std::vector<int> avg_USPN115101_mackenzie_thin_616;
  std::vector<int> avg_USPN115101_mackenzie_thick_616;
  std::vector<int> avg_USPN115101_mackenzie_thin_714;
  std::vector<int> avg_USPN115101_mackenzie_thick_714;
  std::vector<int> avg_USPN115101_mackenzie_thin_728;
  std::vector<int> avg_USPN115101_mackenzie_thick_728;
  std::vector<int> avg_USPN115101_mackenzie_thin_811;
  std::vector<int> avg_USPN115101_mackenzie_thick_811;
  std::vector<int> avg_USPN115101_mackenzie_thin_825;
  std::vector<int> avg_USPN115101_mackenzie_thick_825;
  std::vector<int> avg_USPN115101_mackenzie_thin_908;
  std::vector<int> avg_USPN115101_mackenzie_thick_908;
  std::vector<int> avg_USPN115101_mackenzie_thin_920;
  std::vector<int> avg_USPN115101_mackenzie_thick_920;
  std::vector<int> avg_USPN115101_mackenzie_thin_1004;
  std::vector<int> avg_USPN115101_mackenzie_thick_1004;
  std::vector<int> avg_USPN115101_mackenzie_thin_1018;
  std::vector<int> avg_USPN115101_mackenzie_thick_1018;
  std::vector<int> avg_USPN115101_mackenzie_thin_1101;
  std::vector<int> avg_USPN115101_mackenzie_thick_1101;
  std::vector<int> avg_USPN115101_mackenzie_thin_1115;
  std::vector<int> avg_USPN115101_mackenzie_thick_1115;
  std::vector<int> avg_USPN115101_mackenzie_thin_1129;
  std::vector<int> avg_USPN115101_mackenzie_thick_1129;

  std::vector<int> avg_USPN115101_brian_thin_020123;
  std::vector<int> avg_USPN115101_brian_thick_020123;
  std::vector<int> avg_USPN115101_brian_thin_020923;
  std::vector<int> avg_USPN115101_brian_thick_020923;
  std::vector<int> avg_USPN115101_brian_thin_021723;
  std::vector<int> avg_USPN115101_brian_thick_021723;
  std::vector<int> avg_USPN115101_brian_thin_030923;
  std::vector<int> avg_USPN115101_brian_thick_030923;
  std::vector<int> avg_USPN115101_brian_thin_031623;
  std::vector<int> avg_USPN115101_brian_thick_031623;
  std::vector<int> avg_USPN115101_brian_thin_032323;
  std::vector<int> avg_USPN115101_brian_thick_032323;
  std::vector<int> avg_USPN115101_brian_thin_040623;
  std::vector<int> avg_USPN115101_brian_thick_040623;
  std::vector<int> avg_USPN115101_brian_thin_041323;
  std::vector<int> avg_USPN115101_brian_thick_041323;
  std::vector<int> avg_USPN115101_brian_thin_042123;
  std::vector<int> avg_USPN115101_brian_thick_042123;
  std::vector<int> avg_USPN115101_brian_thin_052423;
  std::vector<int> avg_USPN115101_brian_thick_052423;
  std::vector<int> avg_USPN115101_brian_thin_060823;
  std::vector<int> avg_USPN115101_brian_thick_060823;
  std::vector<int> avg_USPN115101_brian_thin_062323;
  std::vector<int> avg_USPN115101_brian_thick_062323;

  std::vector<int> stab_USPN115101_mackenzie_thin_920;
  std::vector<int> stab_USPN115101_mackenzie_thick_920;

  
  // here is where we actually do the grouping and push indices into the above vectors
  int index = 0;
  for(std::string& s : key_strs){
    if(s.find("USPN085939_mackenzie_thin") != std::string::npos){
      //std::cout << s << std::endl;
      avg_USPN085939_mackenzie_thin.push_back(index);
      //std::cout << "index = " << index << " val check: " << keyNames[index].Data() << std::endl; 
    }
    if(s.find("USPN085939_mackenzie_thick") != std::string::npos){
      avg_USPN085939_mackenzie_thick.push_back(index);
    }    
    if(s.find("USPN085939_brian_thin") != std::string::npos){
      avg_USPN085939_brian_thin.push_back(index);
    }
    if(s.find("USPN085939_brian_thick") != std::string::npos){
      avg_USPN085939_brian_thick.push_back(index);
    }

    if(s.find("USPN024188_mackenzie_thick") != std::string::npos){
      avg_USPN024188_mackenzie_thick.push_back(index);
    }
    if(s.find("USPN024188_mackenzie_thin") != std::string::npos){
      avg_USPN024188_mackenzie_thin.push_back(index);
    }
    if(s.find("USPN024188_brian_thick") != std::string::npos){
      avg_USPN024188_brian_thick.push_back(index);
    }
    if(s.find("USPN024188_brian_thin") != std::string::npos){
      avg_USPN024188_brian_thin.push_back(index);
    }

    if(s.find("USPN035597_mackenzie_thick") != std::string::npos){
      avg_USPN035597_mackenzie_thick.push_back(index);
    }
    if(s.find("USPN035597_mackenzie_thin") != std::string::npos){
      avg_USPN035597_mackenzie_thin.push_back(index);
    }
    if(s.find("USPN035597_brian_thick") != std::string::npos){
      avg_USPN035597_brian_thick.push_back(index);
    }
    if(s.find("USPN035597_brian_thin") != std::string::npos){
      avg_USPN035597_brian_thin.push_back(index);
    }

    if(s.find("USPN100423_mackenzie_thick") != std::string::npos){
      avg_USPN100423_mackenzie_thick.push_back(index);
    }
    if(s.find("USPN100423_mackenzie_thin") != std::string::npos){
      avg_USPN100423_mackenzie_thin.push_back(index);
    }
    if(s.find("USPN100423_brian_thick") != std::string::npos){
      avg_USPN100423_brian_thick.push_back(index);
    }
    if(s.find("USPN100423_brian_thin") != std::string::npos){
      avg_USPN100423_brian_thin.push_back(index);
    }

    if(s.find("USPN015794_mackenzie_thick") != std::string::npos){
      avg_USPN015794_mackenzie_thick.push_back(index);
    }
    if(s.find("USPN015794_mackenzie_thin") != std::string::npos){
      avg_USPN015794_mackenzie_thin.push_back(index);
    }
    if(s.find("USPN015794_brian_thick") != std::string::npos){
      avg_USPN015794_brian_thick.push_back(index);
    }
    if(s.find("USPN015794_brian_thin") != std::string::npos){
      avg_USPN015794_brian_thin.push_back(index);
    }

    if(s.find("USPN014804_mackenzie_thick") != std::string::npos){
      avg_USPN014804_mackenzie_thick.push_back(index);
    }
    if(s.find("USPN014804_mackenzie_thin") != std::string::npos){
      avg_USPN014804_mackenzie_thin.push_back(index);
    }

    if(s.find("KE08508_mackenzie_thick") != std::string::npos){
      avg_KE08508_mackenzie_thick.push_back(index);
    }
    if(s.find("KE08508_mackenzie_thin") != std::string::npos){
      avg_KE08508_mackenzie_thin.push_back(index);
    }

    if(s.find("KE05173_mackenzie_thick") != std::string::npos){
      avg_KE05173_mackenzie_thick.push_back(index);
    }
    if(s.find("KE05173_mackenzie_mid") != std::string::npos){
      avg_KE05173_mackenzie_mid.push_back(index);
    }
    if(s.find("KE05173_mackenzie_thin") != std::string::npos){
      avg_KE05173_mackenzie_thin.push_back(index);
    }

    if(s.find("WHC26311A_mackenzie_thick") != std::string::npos){
      avg_WHC26311A_mackenzie_thick.push_back(index);
    }
    if(s.find("WHC26311A_mackenzie_mid") != std::string::npos){
      avg_WHC26311A_mackenzie_mid.push_back(index);
    }
    if(s.find("WHC26311A_mackenzie_thin") != std::string::npos){
      avg_WHC26311A_mackenzie_thin.push_back(index);
    }

    if(s.find("KE19161_mackenzie_thick") != std::string::npos){
      avg_KE19161_mackenzie_thick.push_back(index);
    }
    if(s.find("KE19161_mackenzie_thin") != std::string::npos){
      avg_KE19161_mackenzie_thin.push_back(index);
    }
    if(s.find("KE19161_brian_thick") != std::string::npos){
      avg_KE19161_brian_thick.push_back(index);
    }
    if(s.find("KE19161_brian_thin") != std::string::npos){
      avg_KE19161_brian_thin.push_back(index);
    }

    if(s.find("KE12291_mackenzie_thick") != std::string::npos){
      avg_KE12291_mackenzie_thick.push_back(index);
    }
    if(s.find("KE12291_mackenzie_mid") != std::string::npos){
      avg_KE12291_mackenzie_mid.push_back(index);
    }
    if(s.find("KE12291_mackenzie_thin") != std::string::npos){
      avg_KE12291_mackenzie_thin.push_back(index);
    }
    if(s.find("KE12291_brian_thick") != std::string::npos){
      avg_KE12291_brian_thick.push_back(index);
    }
    if(s.find("KE12291_brian_thin") != std::string::npos){
      avg_KE12291_brian_thin.push_back(index);
    }

    if(s.find("KE23167_mackenzie_thick") != std::string::npos){
      avg_KE23167_mackenzie_thick.push_back(index);
    }
    if(s.find("KE23167_mackenzie_thin") != std::string::npos){
      avg_KE23167_mackenzie_thin.push_back(index);
    }
    if(s.find("KE23167_brian_thick") != std::string::npos){
      avg_KE23167_brian_thick.push_back(index);
    }
    if(s.find("KE23167_brian_thin") != std::string::npos){
      avg_KE23167_brian_thin.push_back(index);
    }

    if(s.find("JE05501_mackenzie_thick") != std::string::npos){
      avg_JE05501_mackenzie_thick.push_back(index);
    }
    if(s.find("JE05501_mackenzie_thin") != std::string::npos){
      avg_JE05501_mackenzie_thin.push_back(index);
    }

    if(s.find("KE05130_mackenzie_thick") != std::string::npos){
      avg_KE05130_mackenzie_thick.push_back(index);
    }
    if(s.find("KE05130_mackenzie_thin") != std::string::npos){
      avg_KE05130_mackenzie_thin.push_back(index);
    }

    if(s.find("IE08144_mackenzie_thick") != std::string::npos){
      avg_IE08144_mackenzie_thick.push_back(index);
    }
    if(s.find("IE08144_mackenzie_thin") != std::string::npos){
      avg_IE08144_mackenzie_thin.push_back(index);
    }

    if(s.find("IE08146_mackenzie_thick") != std::string::npos){
      avg_IE08146_mackenzie_thick.push_back(index);
    }
    if(s.find("IE08146_mackenzie_thin") != std::string::npos){
      avg_IE08146_mackenzie_thin.push_back(index);
    }

    if(s.find("IE00108_mackenzie_thick") != std::string::npos){
      avg_IE00108_mackenzie_thick.push_back(index);
    }
    if(s.find("IE00108_mackenzie_mid") != std::string::npos){
      avg_IE00108_mackenzie_mid.push_back(index);
    }
    if(s.find("IE00108_mackenzie_thin") != std::string::npos){
      avg_IE00108_mackenzie_thin.push_back(index);
    }

    if(s.find("IE00109_mackenzie_thick") != std::string::npos){
      avg_IE00109_mackenzie_thick.push_back(index);
    }
    if(s.find("IE00109_mackenzie_mid") != std::string::npos){
      avg_IE00109_mackenzie_mid.push_back(index);
    }
    if(s.find("IE00109_mackenzie_thin") != std::string::npos){
      avg_IE00109_mackenzie_thin.push_back(index);
    }


    // new sample USPN115101
    if(s.find("USPN115101_mackenzie_thick_317") != std::string::npos){
      avg_USPN115101_mackenzie_thick_317.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_317") != std::string::npos){
      avg_USPN115101_mackenzie_thin_317.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_324") != std::string::npos){
      avg_USPN115101_mackenzie_thick_324.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_324") != std::string::npos){
      avg_USPN115101_mackenzie_thin_324.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_407") != std::string::npos){
      avg_USPN115101_mackenzie_thick_407.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_407") != std::string::npos){
      avg_USPN115101_mackenzie_thin_407.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_414") != std::string::npos){
      avg_USPN115101_mackenzie_thick_414.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_414") != std::string::npos){
      avg_USPN115101_mackenzie_thin_414.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_421") != std::string::npos){
      avg_USPN115101_mackenzie_thick_421.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_421") != std::string::npos){
      avg_USPN115101_mackenzie_thin_421.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_428") != std::string::npos){
      avg_USPN115101_mackenzie_thick_428.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_428") != std::string::npos){
      avg_USPN115101_mackenzie_thin_428.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_505") != std::string::npos){
      avg_USPN115101_mackenzie_thick_505.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_505") != std::string::npos){
      avg_USPN115101_mackenzie_thin_505.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_512") != std::string::npos){
      avg_USPN115101_mackenzie_thick_512.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_512") != std::string::npos){
      avg_USPN115101_mackenzie_thin_512.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_526") != std::string::npos){
      avg_USPN115101_mackenzie_thick_526.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_526") != std::string::npos){
      avg_USPN115101_mackenzie_thin_526.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_602") != std::string::npos){
      avg_USPN115101_mackenzie_thick_602.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_602") != std::string::npos){
      avg_USPN115101_mackenzie_thin_602.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_616") != std::string::npos){
      avg_USPN115101_mackenzie_thick_616.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_616") != std::string::npos){
      avg_USPN115101_mackenzie_thin_616.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_714") != std::string::npos){
      avg_USPN115101_mackenzie_thick_714.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_714") != std::string::npos){
      avg_USPN115101_mackenzie_thin_714.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_728") != std::string::npos){
      avg_USPN115101_mackenzie_thick_728.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_728") != std::string::npos){
      avg_USPN115101_mackenzie_thin_728.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_811") != std::string::npos){
      avg_USPN115101_mackenzie_thick_811.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_811") != std::string::npos){
      avg_USPN115101_mackenzie_thin_811.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_825") != std::string::npos){
      avg_USPN115101_mackenzie_thick_825.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_825") != std::string::npos){
      avg_USPN115101_mackenzie_thin_825.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_908") != std::string::npos){
      avg_USPN115101_mackenzie_thick_908.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_908") != std::string::npos){
      avg_USPN115101_mackenzie_thin_908.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_920") != std::string::npos){
      avg_USPN115101_mackenzie_thick_920.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_920") != std::string::npos){
      avg_USPN115101_mackenzie_thin_920.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_1004") != std::string::npos){
      avg_USPN115101_mackenzie_thick_1004.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_1004") != std::string::npos){
      avg_USPN115101_mackenzie_thin_1004.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_1018") != std::string::npos){
      avg_USPN115101_mackenzie_thick_1018.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_1018") != std::string::npos){
      avg_USPN115101_mackenzie_thin_1018.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_1101") != std::string::npos){
      avg_USPN115101_mackenzie_thick_1101.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_1101") != std::string::npos){
      avg_USPN115101_mackenzie_thin_1101.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_1115") != std::string::npos){
      avg_USPN115101_mackenzie_thick_1115.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_1115") != std::string::npos){
      avg_USPN115101_mackenzie_thin_1115.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thick_1129") != std::string::npos){
      avg_USPN115101_mackenzie_thick_1129.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_1129") != std::string::npos){
      avg_USPN115101_mackenzie_thin_1129.push_back(index);
    }

    if(s.find("USPN115101_brian_thick_020123") != std::string::npos){
      avg_USPN115101_brian_thick_020123.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_020123") != std::string::npos){
      avg_USPN115101_brian_thin_020123.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_020923") != std::string::npos){
      avg_USPN115101_brian_thick_020923.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_020923") != std::string::npos){
      avg_USPN115101_brian_thin_020923.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_021723") != std::string::npos){
      avg_USPN115101_brian_thick_021723.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_021723") != std::string::npos){
      avg_USPN115101_brian_thin_021723.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_030923") != std::string::npos){
      avg_USPN115101_brian_thick_030923.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_030923") != std::string::npos){
      avg_USPN115101_brian_thin_030923.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_031623") != std::string::npos){
      avg_USPN115101_brian_thick_031623.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_031623") != std::string::npos){
      avg_USPN115101_brian_thin_031623.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_032323") != std::string::npos){
      avg_USPN115101_brian_thick_032323.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_032323") != std::string::npos){
      avg_USPN115101_brian_thin_032323.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_040623") != std::string::npos){
      avg_USPN115101_brian_thick_040623.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_040623") != std::string::npos){
      avg_USPN115101_brian_thin_040623.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_041323") != std::string::npos){
      avg_USPN115101_brian_thick_041323.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_041323") != std::string::npos){
      avg_USPN115101_brian_thin_041323.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_042123") != std::string::npos){
      avg_USPN115101_brian_thick_042123.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_042123") != std::string::npos){
      avg_USPN115101_brian_thin_042123.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_052423") != std::string::npos){
      avg_USPN115101_brian_thick_052423.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_052423") != std::string::npos){
      avg_USPN115101_brian_thin_052423.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_060823") != std::string::npos){
      avg_USPN115101_brian_thick_060823.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_060823") != std::string::npos){
      avg_USPN115101_brian_thin_060823.push_back(index);
    }
    if(s.find("USPN115101_brian_thick_062323") != std::string::npos){
      avg_USPN115101_brian_thick_062323.push_back(index);
    }
    if(s.find("USPN115101_brian_thin_062323") != std::string::npos){
      avg_USPN115101_brian_thin_062323.push_back(index);
    }

    if(s.find("USPN115101_mackenzie_thick_920_3") != std::string::npos){
      stab_USPN115101_mackenzie_thick_920.push_back(index);
    }
    if(s.find("USPN115101_mackenzie_thin_920_3") != std::string::npos){
      stab_USPN115101_mackenzie_thin_920.push_back(index);
    }

    ++index;
  }

  // now, do fanicier analysis here using groups from above. Averages, differences, etc.

  // average histogram function
  USPN085939_thin_mackenzie_avg = average_hist(input_hists, avg_USPN085939_mackenzie_thin, "avg_USPN085939_mackenzie_thin");
  USPN085939_thin_brian_avg = average_hist(input_hists, avg_USPN085939_brian_thin, "avg_USPN085939_brian_thin");
  USPN085939_thick_mackenzie_avg = average_hist(input_hists, avg_USPN085939_mackenzie_thick, "avg_USPN085939_mackenzie_thick");
  USPN085939_thick_brian_avg = average_hist(input_hists, avg_USPN085939_brian_thick, "avg_USPN085939_brian_thick");

  USPN024188_thin_mackenzie_avg = average_hist(input_hists, avg_USPN024188_mackenzie_thin, "avg_USPN024188_mackenzie_thin");
  USPN024188_thin_brian_avg = average_hist(input_hists, avg_USPN024188_brian_thin, "avg_USPN024188_brian_thin");
  USPN024188_thick_mackenzie_avg = average_hist(input_hists, avg_USPN024188_mackenzie_thick, "avg_USPN024188_mackenzie_thick");
  USPN024188_thick_brian_avg = average_hist(input_hists, avg_USPN024188_brian_thick, "avg_USPN024188_brian_thick");

  USPN035597_thin_mackenzie_avg = average_hist(input_hists, avg_USPN035597_mackenzie_thin, "avg_USPN035597_mackenzie_thin");
  USPN035597_thin_brian_avg = average_hist(input_hists, avg_USPN035597_brian_thin, "avg_USPN035597_brian_thin");
  USPN035597_thick_mackenzie_avg = average_hist(input_hists, avg_USPN035597_mackenzie_thick, "avg_USPN035597_mackenzie_thick");
  USPN035597_thick_brian_avg = average_hist(input_hists, avg_USPN035597_brian_thick, "avg_USPN035597_brian_thick");

  USPN100423_thin_mackenzie_avg = average_hist(input_hists, avg_USPN100423_mackenzie_thin, "avg_USPN100423_mackenzie_thin");
  USPN100423_thin_brian_avg = average_hist(input_hists, avg_USPN100423_brian_thin, "avg_USPN100423_brian_thin");
  USPN100423_thick_mackenzie_avg = average_hist(input_hists, avg_USPN100423_mackenzie_thick, "avg_USPN100423_mackenzie_thick");
  USPN100423_thick_brian_avg = average_hist(input_hists, avg_USPN100423_brian_thick, "avg_USPN100423_brian_thick");

  USPN015794_thin_mackenzie_avg = average_hist(input_hists, avg_USPN015794_mackenzie_thin, "avg_USPN015794_mackenzie_thin");
  USPN015794_thin_brian_avg = average_hist(input_hists, avg_USPN015794_brian_thin, "avg_USPN015794_brian_thin");
  USPN015794_thick_mackenzie_avg = average_hist(input_hists, avg_USPN015794_mackenzie_thick, "avg_USPN015794_mackenzie_thick");
  USPN015794_thick_brian_avg = average_hist(input_hists, avg_USPN015794_brian_thick, "avg_USPN015794_brian_thick");

  USPN014804_thin_mackenzie_avg = average_hist(input_hists, avg_USPN014804_mackenzie_thin, "avg_USPN014804_mackenzie_thin");
  USPN014804_thick_mackenzie_avg = average_hist(input_hists, avg_USPN014804_mackenzie_thick, "avg_USPN014804_mackenzie_thick");

  KE08508_thin_mackenzie_avg = average_hist(input_hists, avg_KE08508_mackenzie_thin, "avg_KE08508_mackenzie_thin");
  KE08508_thick_mackenzie_avg = average_hist(input_hists, avg_KE08508_mackenzie_thick, "avg_KE08508_mackenzie_thick");

  KE05173_thin_mackenzie_avg = average_hist(input_hists, avg_KE05173_mackenzie_thin, "avg_KE05173_mackenzie_thin");
  KE05173_mid_mackenzie_avg = average_hist(input_hists, avg_KE05173_mackenzie_mid, "avg_KE05173_mackenzie_mid");
  KE05173_thick_mackenzie_avg = average_hist(input_hists, avg_KE05173_mackenzie_thick, "avg_KE05173_mackenzie_thick");

  WHC26311A_thin_mackenzie_avg = average_hist(input_hists, avg_WHC26311A_mackenzie_thin, "avg_WHC26311A_mackenzie_thin");
  WHC26311A_mid_mackenzie_avg = average_hist(input_hists, avg_WHC26311A_mackenzie_mid, "avg_WHC26311A_mackenzie_mid");
  WHC26311A_thick_mackenzie_avg = average_hist(input_hists, avg_WHC26311A_mackenzie_thick, "avg_WHC26311A_mackenzie_thick");

  KE19161_thin_brian_avg = average_hist(input_hists, avg_KE19161_brian_thin, "avg_KE19161_brian_thin");
  KE19161_thick_brian_avg = average_hist(input_hists, avg_KE19161_brian_thick, "avg_KE19161_brian_thick");
  KE19161_thin_mackenzie_avg = average_hist(input_hists, avg_KE19161_mackenzie_thin, "avg_KE19161_mackenzie_thin");
  KE19161_thick_mackenzie_avg = average_hist(input_hists, avg_KE19161_mackenzie_thick, "avg_KE19161_mackenzie_thick");

  KE12291_thin_brian_avg = average_hist(input_hists, avg_KE12291_brian_thin, "avg_KE12291_brian_thin");
  KE12291_thick_brian_avg = average_hist(input_hists, avg_KE12291_brian_thick, "avg_KE12291_brian_thick");
  KE12291_thin_mackenzie_avg = average_hist(input_hists, avg_KE12291_mackenzie_thin, "avg_KE12291_mackenzie_thin");
  KE12291_mid_mackenzie_avg = average_hist(input_hists, avg_KE12291_mackenzie_mid, "avg_KE12291_mackenzie_mid");
  KE12291_thick_mackenzie_avg = average_hist(input_hists, avg_KE12291_mackenzie_thick, "avg_KE12291_mackenzie_thick");

  KE23167_thin_brian_avg = average_hist(input_hists, avg_KE23167_brian_thin, "avg_KE23167_brian_thin");
  KE23167_thick_brian_avg = average_hist(input_hists, avg_KE23167_brian_thick, "avg_KE23167_brian_thick");
  KE23167_thin_mackenzie_avg = average_hist(input_hists, avg_KE23167_mackenzie_thin, "avg_KE23167_mackenzie_thin");
  KE23167_thick_mackenzie_avg = average_hist(input_hists, avg_KE23167_mackenzie_thick, "avg_KE23167_mackenzie_thick");

  JE05501_thin_mackenzie_avg = average_hist(input_hists, avg_JE05501_mackenzie_thin, "avg_JE05501_mackenzie_thin");
  JE05501_thick_mackenzie_avg = average_hist(input_hists, avg_JE05501_mackenzie_thick, "avg_JE05501_mackenzie_thick");

  KE05130_thin_mackenzie_avg = average_hist(input_hists, avg_KE05130_mackenzie_thin, "avg_KE05130_mackenzie_thin");
  KE05130_thick_mackenzie_avg = average_hist(input_hists, avg_KE05130_mackenzie_thick, "avg_KE05130_mackenzie_thick");

  IE08144_thin_mackenzie_avg = average_hist(input_hists, avg_IE08144_mackenzie_thin, "avg_IE08144_mackenzie_thin");
  IE08144_thick_mackenzie_avg = average_hist(input_hists, avg_IE08144_mackenzie_thick, "avg_IE08144_mackenzie_thick");

  IE08146_thin_mackenzie_avg = average_hist(input_hists, avg_IE08146_mackenzie_thin, "avg_IE08146_mackenzie_thin");
  IE08146_thick_mackenzie_avg = average_hist(input_hists, avg_IE08146_mackenzie_thick, "avg_IE08146_mackenzie_thick");

  IE00108_thin_mackenzie_avg = average_hist(input_hists, avg_IE00108_mackenzie_thin, "avg_IE00108_mackenzie_thin");
  IE00108_mid_mackenzie_avg = average_hist(input_hists, avg_IE00108_mackenzie_mid, "avg_IE00108_mackenzie_mid");
  IE00108_thick_mackenzie_avg = average_hist(input_hists, avg_IE00108_mackenzie_thick, "avg_IE00108_mackenzie_thick");

  IE00109_thin_mackenzie_avg = average_hist(input_hists, avg_IE00109_mackenzie_thin, "avg_IE00109_mackenzie_thin");
  IE00109_mid_mackenzie_avg = average_hist(input_hists, avg_IE00109_mackenzie_mid, "avg_IE00109_mackenzie_mid");
  IE00109_thick_mackenzie_avg = average_hist(input_hists, avg_IE00109_mackenzie_thick, "avg_IE00109_mackenzie_thick");

  USPN115101_thin_mackenzie_317_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_317, "avg_USPN115101_mackenzie_thin_317");
  USPN115101_thick_mackenzie_317_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_317, "avg_USPN115101_mackenzie_thick_317");
  USPN115101_thin_mackenzie_324_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_324, "avg_USPN115101_mackenzie_thin_324");
  USPN115101_thick_mackenzie_324_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_324, "avg_USPN115101_mackenzie_thick_324");
  USPN115101_thin_mackenzie_407_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_407, "avg_USPN115101_mackenzie_thin_407");
  USPN115101_thick_mackenzie_407_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_407, "avg_USPN115101_mackenzie_thick_407");
  USPN115101_thin_mackenzie_414_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_414, "avg_USPN115101_mackenzie_thin_414");
  USPN115101_thick_mackenzie_414_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_414, "avg_USPN115101_mackenzie_thick_414");
  USPN115101_thin_mackenzie_421_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_421, "avg_USPN115101_mackenzie_thin_421");
  USPN115101_thick_mackenzie_421_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_421, "avg_USPN115101_mackenzie_thick_421");
  USPN115101_thin_mackenzie_428_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_428, "avg_USPN115101_mackenzie_thin_428");
  USPN115101_thick_mackenzie_428_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_428, "avg_USPN115101_mackenzie_thick_428");
  USPN115101_thin_mackenzie_505_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_505, "avg_USPN115101_mackenzie_thin_505");
  USPN115101_thick_mackenzie_505_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_505, "avg_USPN115101_mackenzie_thick_505");
  USPN115101_thin_mackenzie_512_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_512, "avg_USPN115101_mackenzie_thin_512");
  USPN115101_thick_mackenzie_512_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_512, "avg_USPN115101_mackenzie_thick_512");
  USPN115101_thin_mackenzie_526_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_526, "avg_USPN115101_mackenzie_thin_526");
  USPN115101_thick_mackenzie_526_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_526, "avg_USPN115101_mackenzie_thick_526");
  USPN115101_thin_mackenzie_602_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_602, "avg_USPN115101_mackenzie_thin_602");
  USPN115101_thick_mackenzie_602_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_602, "avg_USPN115101_mackenzie_thick_602");
  USPN115101_thin_mackenzie_616_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_616, "avg_USPN115101_mackenzie_thin_616");
  USPN115101_thick_mackenzie_616_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_616, "avg_USPN115101_mackenzie_thick_616");
  USPN115101_thin_mackenzie_714_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_714, "avg_USPN115101_mackenzie_thin_714");
  USPN115101_thick_mackenzie_714_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_714, "avg_USPN115101_mackenzie_thick_714");
  USPN115101_thin_mackenzie_728_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_728, "avg_USPN115101_mackenzie_thin_728");
  USPN115101_thick_mackenzie_728_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_728, "avg_USPN115101_mackenzie_thick_728");
  USPN115101_thin_mackenzie_811_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_811, "avg_USPN115101_mackenzie_thin_811");
  USPN115101_thick_mackenzie_811_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_811, "avg_USPN115101_mackenzie_thick_811");
  USPN115101_thin_mackenzie_825_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_825, "avg_USPN115101_mackenzie_thin_825");
  USPN115101_thick_mackenzie_825_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_825, "avg_USPN115101_mackenzie_thick_825");
  USPN115101_thin_mackenzie_908_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_908, "avg_USPN115101_mackenzie_thin_908");
  USPN115101_thick_mackenzie_908_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_908, "avg_USPN115101_mackenzie_thick_908");
  USPN115101_thin_mackenzie_920_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_920, "avg_USPN115101_mackenzie_thin_920");
  USPN115101_thick_mackenzie_920_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_920, "avg_USPN115101_mackenzie_thick_920");
  USPN115101_thin_mackenzie_1004_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_1004, "avg_USPN115101_mackenzie_thin_1004");
  USPN115101_thick_mackenzie_1004_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_1004, "avg_USPN115101_mackenzie_thick_1004");
  USPN115101_thin_mackenzie_1018_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_1018, "avg_USPN115101_mackenzie_thin_1018");
  USPN115101_thick_mackenzie_1018_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_1018, "avg_USPN115101_mackenzie_thick_1018");
  USPN115101_thin_mackenzie_1101_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_1101, "avg_USPN115101_mackenzie_thin_1101");
  USPN115101_thick_mackenzie_1101_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_1101, "avg_USPN115101_mackenzie_thick_1101");
  USPN115101_thin_mackenzie_1115_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_1115, "avg_USPN115101_mackenzie_thin_1115");
  USPN115101_thick_mackenzie_1115_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_1115, "avg_USPN115101_mackenzie_thick_1115");
  USPN115101_thin_mackenzie_1129_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thin_1129, "avg_USPN115101_mackenzie_thin_1129");
  USPN115101_thick_mackenzie_1129_avg = average_hist(input_hists, avg_USPN115101_mackenzie_thick_1129, "avg_USPN115101_mackenzie_thick_1129");

  USPN115101_thin_brian_020123_avg = average_hist(input_hists, avg_USPN115101_brian_thin_020123, "avg_USPN115101_brian_thin_020123");
  USPN115101_thick_brian_020123_avg = average_hist(input_hists, avg_USPN115101_brian_thick_020123, "avg_USPN115101_brian_thick_020123");
  USPN115101_thin_brian_020923_avg = average_hist(input_hists, avg_USPN115101_brian_thin_020923, "avg_USPN115101_brian_thin_020923");
  USPN115101_thick_brian_020923_avg = average_hist(input_hists, avg_USPN115101_brian_thick_020923, "avg_USPN115101_brian_thick_020923");
  USPN115101_thin_brian_021723_avg = average_hist(input_hists, avg_USPN115101_brian_thin_021723, "avg_USPN115101_brian_thin_021723");
  USPN115101_thick_brian_021723_avg = average_hist(input_hists, avg_USPN115101_brian_thick_021723, "avg_USPN115101_brian_thick_021723");
  USPN115101_thin_brian_030923_avg = average_hist(input_hists, avg_USPN115101_brian_thin_030923, "avg_USPN115101_brian_thin_030923");
  USPN115101_thick_brian_030923_avg = average_hist(input_hists, avg_USPN115101_brian_thick_030923, "avg_USPN115101_brian_thick_030923");
  USPN115101_thin_brian_031623_avg = average_hist(input_hists, avg_USPN115101_brian_thin_031623, "avg_USPN115101_brian_thin_031623");
  USPN115101_thick_brian_031623_avg = average_hist(input_hists, avg_USPN115101_brian_thick_031623, "avg_USPN115101_brian_thick_031623");
  USPN115101_thin_brian_032323_avg = average_hist(input_hists, avg_USPN115101_brian_thin_032323, "avg_USPN115101_brian_thin_032323");
  USPN115101_thick_brian_032323_avg = average_hist(input_hists, avg_USPN115101_brian_thick_032323, "avg_USPN115101_brian_thick_032323");
  USPN115101_thin_brian_040623_avg = average_hist(input_hists, avg_USPN115101_brian_thin_040623, "avg_USPN115101_brian_thin_040623");
  USPN115101_thick_brian_040623_avg = average_hist(input_hists, avg_USPN115101_brian_thick_040623, "avg_USPN115101_brian_thick_040623");
  USPN115101_thin_brian_041323_avg = average_hist(input_hists, avg_USPN115101_brian_thin_041323, "avg_USPN115101_brian_thin_041323");
  USPN115101_thick_brian_041323_avg = average_hist(input_hists, avg_USPN115101_brian_thick_041323, "avg_USPN115101_brian_thick_041323");
  USPN115101_thin_brian_042123_avg = average_hist(input_hists, avg_USPN115101_brian_thin_042123, "avg_USPN115101_brian_thin_042123");
  USPN115101_thick_brian_042123_avg = average_hist(input_hists, avg_USPN115101_brian_thick_042123, "avg_USPN115101_brian_thick_042123");
  USPN115101_thin_brian_052423_avg = average_hist(input_hists, avg_USPN115101_brian_thin_052423, "avg_USPN115101_brian_thin_052423");
  USPN115101_thick_brian_052423_avg = average_hist(input_hists, avg_USPN115101_brian_thick_052423, "avg_USPN115101_brian_thick_052423");
  USPN115101_thin_brian_060823_avg = average_hist(input_hists, avg_USPN115101_brian_thin_060823, "avg_USPN115101_brian_thin_060823");
  USPN115101_thick_brian_060823_avg = average_hist(input_hists, avg_USPN115101_brian_thick_060823, "avg_USPN115101_brian_thick_060823");
  USPN115101_thin_brian_062323_avg = average_hist(input_hists, avg_USPN115101_brian_thin_062323, "avg_USPN115101_brian_thin_062323");
  USPN115101_thick_brian_062323_avg = average_hist(input_hists, avg_USPN115101_brian_thick_062323, "avg_USPN115101_brian_thick_062323");


  // difference histogram function
  USPN085939_thin_diff = difference_hist(USPN085939_thin_brian_avg, USPN085939_thin_mackenzie_avg, "diff_USPN085939_thin_b-m");
  USPN085939_thick_diff = difference_hist(USPN085939_thick_brian_avg, USPN085939_thick_mackenzie_avg, "diff_USPN085939_thick_b-m");

  USPN024188_thin_diff = difference_hist(USPN024188_thin_brian_avg, USPN024188_thin_mackenzie_avg, "diff_USPN024188_thin_b-m");
  USPN024188_thick_diff = difference_hist(USPN024188_thick_brian_avg, USPN024188_thick_mackenzie_avg, "diff_USPN024188_thick_b-m");

  USPN035597_thin_diff = difference_hist(USPN035597_thin_brian_avg, USPN035597_thin_mackenzie_avg, "diff_USPN035597_thin_b-m");
  USPN035597_thick_diff = difference_hist(USPN035597_thick_brian_avg, USPN035597_thick_mackenzie_avg, "diff_USPN035597_thick_b-m");

  USPN100423_thin_diff = difference_hist(USPN100423_thin_brian_avg, USPN100423_thin_mackenzie_avg, "diff_USPN100423_thin_b-m");
  USPN100423_thick_diff = difference_hist(USPN100423_thick_brian_avg, USPN100423_thick_mackenzie_avg, "diff_USPN100423_thick_b-m");

  USPN015794_thin_diff = difference_hist(USPN015794_thin_brian_avg, USPN015794_thin_mackenzie_avg, "diff_USPN015794_thin_b-m");
  USPN015794_thick_diff = difference_hist(USPN015794_thick_brian_avg, USPN015794_thick_mackenzie_avg, "diff_USPN015794_thick_b-m");

  KE19161_thin_diff = difference_hist(KE19161_thin_brian_avg, KE19161_thin_mackenzie_avg, "diff_KE19161_thin_b-m");
  KE19161_thick_diff = difference_hist(KE19161_thick_brian_avg, KE19161_thick_mackenzie_avg, "diff_KE19161_thick_b-m");

  KE12291_thin_diff = difference_hist(KE12291_thin_brian_avg, KE12291_thin_mackenzie_avg, "diff_KE12291_thin_b-m");
  KE12291_thick_diff = difference_hist(KE12291_thick_brian_avg, KE12291_thick_mackenzie_avg, "diff_KE12291_thick_b-m");

  KE23167_thin_diff = difference_hist(KE23167_thin_brian_avg, KE23167_thin_mackenzie_avg, "diff_KE23167_thin_b-m");
  KE23167_thick_diff = difference_hist(KE23167_thick_brian_avg, KE23167_thick_mackenzie_avg, "diff_KE23167_thick_b-m");

  USPN115101_mackenzie_1115_diff = difference_hist(USPN115101_thick_mackenzie_1115_avg, USPN115101_thin_mackenzie_1115_avg, "diff_USPN115101_m_thick-thin");

  // Make sure ROOT knows we want to write to the output file
  ofile->cd();

  // write the histograms to file
  USPN085939_thin_mackenzie_avg->Write();
  USPN085939_thin_brian_avg->Write();
  USPN085939_thin_diff->Write();
  USPN085939_thick_mackenzie_avg->Write();
  USPN085939_thick_brian_avg->Write();
  USPN085939_thick_diff->Write();

  USPN024188_thin_mackenzie_avg->Write();
  USPN024188_thin_brian_avg->Write();
  USPN024188_thin_diff->Write();
  USPN024188_thick_mackenzie_avg->Write();
  USPN024188_thick_brian_avg->Write();
  USPN024188_thick_diff->Write();

  USPN035597_thin_mackenzie_avg->Write();
  USPN035597_thin_brian_avg->Write();
  USPN035597_thin_diff->Write();
  USPN035597_thick_mackenzie_avg->Write();
  USPN035597_thick_brian_avg->Write();
  USPN035597_thick_diff->Write();

  USPN100423_thin_mackenzie_avg->Write();
  USPN100423_thin_brian_avg->Write();
  USPN100423_thin_diff->Write();
  USPN100423_thick_mackenzie_avg->Write();
  USPN100423_thick_brian_avg->Write();
  USPN100423_thick_diff->Write();

  USPN015794_thin_mackenzie_avg->Write();
  USPN015794_thin_brian_avg->Write();
  USPN015794_thin_diff->Write();
  USPN015794_thick_mackenzie_avg->Write();
  USPN015794_thick_brian_avg->Write();
  USPN015794_thick_diff->Write();

  USPN014804_thin_mackenzie_avg->Write();
  USPN014804_thick_mackenzie_avg->Write();

  KE08508_thin_mackenzie_avg->Write();
  KE08508_thick_mackenzie_avg->Write();

  KE05173_thin_mackenzie_avg->Write();
  KE05173_mid_mackenzie_avg->Write();
  KE05173_thick_mackenzie_avg->Write();

  WHC26311A_thin_mackenzie_avg->Write();
  WHC26311A_mid_mackenzie_avg->Write();
  WHC26311A_thick_mackenzie_avg->Write();

  KE19161_thin_mackenzie_avg->Write();
  KE19161_thin_brian_avg->Write();
  KE19161_thin_diff->Write();
  KE19161_thick_mackenzie_avg->Write();
  KE19161_thick_brian_avg->Write();
  KE19161_thick_diff->Write();

  KE12291_thin_mackenzie_avg->Write();
  KE12291_thin_brian_avg->Write();
  KE12291_thin_diff->Write();
  KE12291_mid_mackenzie_avg->Write();
  KE12291_thick_mackenzie_avg->Write();
  KE12291_thick_brian_avg->Write();
  KE12291_thick_diff->Write();

  KE23167_thin_mackenzie_avg->Write();
  KE23167_thin_brian_avg->Write();
  KE23167_thin_diff->Write();
  KE23167_thick_mackenzie_avg->Write();
  KE23167_thick_brian_avg->Write();
  KE23167_thick_diff->Write();

  JE05501_thin_mackenzie_avg->Write();
  JE05501_thick_mackenzie_avg->Write();

  KE05130_thin_mackenzie_avg->Write();
  KE05130_thick_mackenzie_avg->Write();

  IE08144_thin_mackenzie_avg->Write();
  IE08144_thick_mackenzie_avg->Write();

  IE08146_thin_mackenzie_avg->Write();
  IE08146_thick_mackenzie_avg->Write();

  IE00108_thin_mackenzie_avg->Write();
  IE00108_mid_mackenzie_avg->Write();
  IE00108_thick_mackenzie_avg->Write();

  IE00109_thin_mackenzie_avg->Write();
  IE00109_mid_mackenzie_avg->Write();
  IE00109_thick_mackenzie_avg->Write();

  USPN115101_thin_mackenzie_317_avg->Write();
  USPN115101_thick_mackenzie_317_avg->Write();
  USPN115101_thin_mackenzie_324_avg->Write();
  USPN115101_thick_mackenzie_324_avg->Write();
  USPN115101_thin_mackenzie_407_avg->Write();
  USPN115101_thick_mackenzie_407_avg->Write();
  USPN115101_thin_mackenzie_414_avg->Write();
  USPN115101_thick_mackenzie_414_avg->Write();
  USPN115101_thin_mackenzie_421_avg->Write();
  USPN115101_thick_mackenzie_421_avg->Write();
  USPN115101_thin_mackenzie_428_avg->Write();
  USPN115101_thick_mackenzie_428_avg->Write();
  USPN115101_thin_mackenzie_505_avg->Write();
  USPN115101_thick_mackenzie_505_avg->Write();
  USPN115101_thin_mackenzie_512_avg->Write();
  USPN115101_thick_mackenzie_512_avg->Write();
  USPN115101_thin_mackenzie_526_avg->Write();
  USPN115101_thick_mackenzie_526_avg->Write();
  USPN115101_thin_mackenzie_602_avg->Write();
  USPN115101_thick_mackenzie_602_avg->Write();
  USPN115101_thin_mackenzie_616_avg->Write();
  USPN115101_thick_mackenzie_616_avg->Write();
  USPN115101_thin_mackenzie_714_avg->Write();
  USPN115101_thick_mackenzie_714_avg->Write();
  USPN115101_thin_mackenzie_728_avg->Write();
  USPN115101_thick_mackenzie_728_avg->Write();
  USPN115101_thin_mackenzie_811_avg->Write();
  USPN115101_thick_mackenzie_811_avg->Write();
  USPN115101_thin_mackenzie_825_avg->Write();
  USPN115101_thick_mackenzie_825_avg->Write();
  USPN115101_thin_mackenzie_908_avg->Write();
  USPN115101_thick_mackenzie_908_avg->Write();
  USPN115101_thin_mackenzie_920_avg->Write();
  USPN115101_thick_mackenzie_920_avg->Write();
  USPN115101_thin_mackenzie_1004_avg->Write();
  USPN115101_thick_mackenzie_1004_avg->Write();
  USPN115101_thin_mackenzie_1018_avg->Write();
  USPN115101_thick_mackenzie_1018_avg->Write();
  USPN115101_thin_mackenzie_1101_avg->Write();
  USPN115101_thick_mackenzie_1101_avg->Write();
  USPN115101_thin_mackenzie_1115_avg->Write();
  USPN115101_thick_mackenzie_1115_avg->Write();
  USPN115101_mackenzie_1115_diff->Write();
  USPN115101_thin_mackenzie_1129_avg->Write();
  USPN115101_thick_mackenzie_1129_avg->Write();

  USPN115101_thin_brian_020123_avg->Write();
  USPN115101_thick_brian_020123_avg->Write();
  USPN115101_thin_brian_020923_avg->Write();
  USPN115101_thick_brian_020923_avg->Write();
  USPN115101_thin_brian_021723_avg->Write();
  USPN115101_thick_brian_021723_avg->Write();
  USPN115101_thin_brian_030923_avg->Write();
  USPN115101_thick_brian_030923_avg->Write();
  USPN115101_thin_brian_031623_avg->Write();
  USPN115101_thick_brian_031623_avg->Write();
  USPN115101_thin_brian_032323_avg->Write();
  USPN115101_thick_brian_032323_avg->Write();
  USPN115101_thin_brian_040623_avg->Write();
  USPN115101_thick_brian_040623_avg->Write();
  USPN115101_thin_brian_041323_avg->Write();
  USPN115101_thick_brian_041323_avg->Write();
  USPN115101_thin_brian_042123_avg->Write();
  USPN115101_thick_brian_042123_avg->Write();
  USPN115101_thin_brian_052423_avg->Write();
  USPN115101_thick_brian_052423_avg->Write();
  USPN115101_thin_brian_060823_avg->Write();
  USPN115101_thick_brian_060823_avg->Write();
  USPN115101_thin_brian_062323_avg->Write();
  USPN115101_thick_brian_062323_avg->Write();

  // make collections of average data for wavelengths of interest, for aging
  std::vector<float> aging_380_115101_thin;
  std::vector<float> aging_380_115101_thick;
  std::vector<float> aging_380_115101_thin_err;
  std::vector<float> aging_380_115101_thick_err;
  std::vector<float> aging_390_115101_thin;
  std::vector<float> aging_390_115101_thick;
  std::vector<float> aging_390_115101_thin_err;
  std::vector<float> aging_390_115101_thick_err;
  std::vector<float> aging_400_115101_thin;
  std::vector<float> aging_400_115101_thick;
  std::vector<float> aging_400_115101_thin_err;
  std::vector<float> aging_400_115101_thick_err;
  std::vector<float> aging_410_115101_thin;
  std::vector<float> aging_410_115101_thick;
  std::vector<float> aging_410_115101_thin_err;
  std::vector<float> aging_410_115101_thick_err;
  std::vector<float> aging_420_115101_thin;
  std::vector<float> aging_420_115101_thick;
  std::vector<float> aging_420_115101_thin_err;
  std::vector<float> aging_420_115101_thick_err;
  std::vector<float> aging_430_115101_thin;
  std::vector<float> aging_430_115101_thick;
  std::vector<float> aging_430_115101_thin_err;
  std::vector<float> aging_430_115101_thick_err;
  std::vector<float> aging_440_115101_thin;
  std::vector<float> aging_440_115101_thick;
  std::vector<float> aging_440_115101_thin_err;
  std::vector<float> aging_440_115101_thick_err;
  std::vector<float> aging_450_115101_thin;
  std::vector<float> aging_450_115101_thick;
  std::vector<float> aging_450_115101_thin_err;
  std::vector<float> aging_450_115101_thick_err;

  std::vector<float> stab_115101_thin;
  std::vector<float> stab_115101_thick;

  // make the x-axis points array, time fraction of years
  const int npoints = 34;
  float aging_timefrac_xbins[npoints] = {0.0, 0.019, 0.058, 0.077, 0.096, 0.115, 0.134, 0.153, 0.192, 0.211, 0.249, 0.326, 0.364, 0.403, 0.441, 0.479, 0.512, 0.551, 0.589, 0.627, 0.666, 0.704, 0.879, 0.901, 0.923, 0.978, 0.997, 1.016, 1.055, 1.074, 1.096, 1.186, 1.227, 1.268};

  float aging_empty_xbins_err[npoints];
  std::fill(aging_empty_xbins_err, aging_empty_xbins_err+npoints, 0);
  //std::fill(aging_400_115101_thin_err, aging_400_115101_thin_err+15, 0);
  //std::fill(aging_400_115101_thick_err, aging_400_115101_thick_err+15, 0);

  std::cout << "bin number for 500 wavelength = " << USPN024188_thick_mackenzie_avg->FindBin(500.) << std::endl;
  std::cout << "bin number for 400 wavelength = " << USPN024188_thick_mackenzie_avg->FindBin(400.) << std::endl;

  // wavelength = 380 bin content
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_317_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_317_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_324_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_324_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_407_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_407_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_414_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_414_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_421_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_421_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_428_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_428_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_505_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_505_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_512_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_512_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_526_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_526_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_602_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_602_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_616_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_616_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_714_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_714_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_728_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_728_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_811_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_811_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_825_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_825_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_908_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_908_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_920_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_920_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_1004_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_1004_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_1018_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_1018_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_1101_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_1101_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_1115_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_1115_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_mackenzie_1129_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_mackenzie_1129_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_020123_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_020123_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_020923_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_020923_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_021723_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_021723_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_030923_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_030923_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_031623_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_031623_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_032323_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_032323_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_040623_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_040623_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_041323_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_041323_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_042123_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_042123_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_052423_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_052423_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_060823_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_060823_avg->GetBinContent(3));
  aging_380_115101_thin.push_back(USPN115101_thin_brian_062323_avg->GetBinContent(3));
  aging_380_115101_thick.push_back(USPN115101_thick_brian_062323_avg->GetBinContent(3));
  // wavelength = 390 bin content
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_317_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_317_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_324_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_324_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_407_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_407_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_414_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_414_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_421_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_421_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_428_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_428_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_505_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_505_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_512_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_512_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_526_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_526_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_602_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_602_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_616_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_616_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_714_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_714_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_728_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_728_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_811_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_811_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_825_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_825_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_908_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_908_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_920_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_920_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_1004_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_1004_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_1018_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_1018_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_1101_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_1101_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_1115_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_1115_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_mackenzie_1129_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_mackenzie_1129_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_020123_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_020123_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_020923_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_020923_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_021723_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_021723_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_030923_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_030923_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_031623_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_031623_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_032323_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_032323_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_040623_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_040623_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_041323_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_041323_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_042123_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_042123_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_052423_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_052423_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_060823_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_060823_avg->GetBinContent(4));
  aging_390_115101_thin.push_back(USPN115101_thin_brian_062323_avg->GetBinContent(4));
  aging_390_115101_thick.push_back(USPN115101_thick_brian_062323_avg->GetBinContent(4));
  // wavelength = 400 bin content
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_317_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_317_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_324_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_324_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_407_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_407_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_414_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_414_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_421_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_421_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_428_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_428_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_505_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_505_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_512_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_512_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_526_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_526_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_602_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_602_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_616_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_616_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_714_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_714_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_728_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_728_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_811_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_811_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_825_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_825_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_908_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_908_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_920_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_920_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_1004_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_1004_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_1018_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_1018_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_1101_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_1101_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_1115_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_1115_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_mackenzie_1129_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_mackenzie_1129_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_020123_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_020123_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_020923_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_020923_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_021723_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_021723_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_030923_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_030923_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_031623_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_031623_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_032323_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_032323_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_040623_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_040623_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_041323_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_041323_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_042123_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_042123_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_052423_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_052423_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_060823_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_060823_avg->GetBinContent(5));
  aging_400_115101_thin.push_back(USPN115101_thin_brian_062323_avg->GetBinContent(5));
  aging_400_115101_thick.push_back(USPN115101_thick_brian_062323_avg->GetBinContent(5));
  // wavelength = 410 bin content
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_317_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_317_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_324_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_324_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_407_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_407_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_414_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_414_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_421_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_421_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_428_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_428_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_505_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_505_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_512_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_512_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_526_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_526_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_602_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_602_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_616_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_616_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_714_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_714_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_728_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_728_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_811_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_811_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_825_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_825_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_908_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_908_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_920_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_920_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_1004_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_1004_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_1018_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_1018_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_1101_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_1101_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_1115_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_1115_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_mackenzie_1129_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_mackenzie_1129_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_020123_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_020123_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_020923_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_020923_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_021723_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_021723_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_030923_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_030923_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_031623_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_031623_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_032323_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_032323_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_040623_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_040623_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_041323_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_041323_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_042123_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_042123_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_052423_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_052423_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_060823_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_060823_avg->GetBinContent(6));
  aging_410_115101_thin.push_back(USPN115101_thin_brian_062323_avg->GetBinContent(6));
  aging_410_115101_thick.push_back(USPN115101_thick_brian_062323_avg->GetBinContent(6));
  // wavelength = 420 bin content
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_317_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_317_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_324_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_324_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_407_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_407_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_414_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_414_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_421_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_421_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_428_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_428_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_505_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_505_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_512_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_512_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_526_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_526_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_602_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_602_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_616_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_616_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_714_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_714_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_728_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_728_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_811_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_811_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_825_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_825_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_908_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_908_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_920_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_920_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_1004_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_1004_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_1018_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_1018_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_1101_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_1101_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_1115_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_1115_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_mackenzie_1129_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_mackenzie_1129_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_020123_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_020123_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_020923_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_020923_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_021723_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_021723_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_030923_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_030923_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_031623_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_031623_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_032323_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_032323_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_040623_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_040623_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_041323_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_041323_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_042123_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_042123_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_052423_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_052423_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_060823_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_060823_avg->GetBinContent(7));
  aging_420_115101_thin.push_back(USPN115101_thin_brian_062323_avg->GetBinContent(7));
  aging_420_115101_thick.push_back(USPN115101_thick_brian_062323_avg->GetBinContent(7));
  // wavelength = 430 bin content
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_317_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_317_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_324_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_324_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_407_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_407_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_414_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_414_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_421_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_421_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_428_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_428_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_505_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_505_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_512_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_512_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_526_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_526_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_602_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_602_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_616_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_616_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_714_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_714_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_728_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_728_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_811_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_811_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_825_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_825_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_908_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_908_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_920_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_920_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_1004_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_1004_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_1018_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_1018_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_1101_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_1101_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_1115_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_1115_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_mackenzie_1129_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_mackenzie_1129_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_020123_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_020123_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_020923_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_020923_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_021723_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_021723_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_030923_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_030923_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_031623_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_031623_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_032323_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_032323_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_040623_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_040623_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_041323_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_041323_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_042123_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_042123_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_052423_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_052423_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_060823_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_060823_avg->GetBinContent(8));
  aging_430_115101_thin.push_back(USPN115101_thin_brian_062323_avg->GetBinContent(8));
  aging_430_115101_thick.push_back(USPN115101_thick_brian_062323_avg->GetBinContent(8));
  // wavelength = 440 bin content
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_317_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_317_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_324_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_324_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_407_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_407_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_414_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_414_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_421_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_421_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_428_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_428_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_505_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_505_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_512_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_512_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_526_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_526_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_602_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_602_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_616_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_616_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_714_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_714_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_728_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_728_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_811_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_811_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_825_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_825_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_908_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_908_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_920_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_920_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_1004_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_1004_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_1018_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_1018_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_1101_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_1101_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_1115_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_1115_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_mackenzie_1129_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_mackenzie_1129_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_020123_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_020123_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_020923_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_020923_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_021723_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_021723_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_030923_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_030923_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_031623_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_031623_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_032323_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_032323_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_040623_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_040623_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_041323_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_041323_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_042123_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_042123_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_052423_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_052423_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_060823_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_060823_avg->GetBinContent(9));
  aging_440_115101_thin.push_back(USPN115101_thin_brian_062323_avg->GetBinContent(9));
  aging_440_115101_thick.push_back(USPN115101_thick_brian_062323_avg->GetBinContent(9));
  // wavelength = 450 bin content
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_317_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_317_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_324_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_324_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_407_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_407_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_414_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_414_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_421_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_421_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_428_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_428_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_505_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_505_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_512_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_512_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_526_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_526_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_602_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_602_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_616_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_616_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_714_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_714_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_728_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_728_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_811_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_811_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_825_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_825_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_908_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_908_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_920_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_920_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_1004_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_1004_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_1018_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_1018_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_1101_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_1101_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_1115_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_1115_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_mackenzie_1129_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_mackenzie_1129_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_020123_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_020123_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_020923_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_020923_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_021723_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_021723_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_030923_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_030923_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_031623_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_031623_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_032323_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_032323_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_040623_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_040623_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_041323_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_041323_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_042123_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_042123_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_052423_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_052423_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_060823_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_060823_avg->GetBinContent(10));
  aging_450_115101_thin.push_back(USPN115101_thin_brian_062323_avg->GetBinContent(10));
  aging_450_115101_thick.push_back(USPN115101_thick_brian_062323_avg->GetBinContent(10));

  // do errors here, standard deviations

  // wavelength = 380nm
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_317, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_324, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_407, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_414, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_421, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_428, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_505, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_512, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_526, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_602, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_616, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_714, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_728, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_811, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_825, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_908, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_920, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1004, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1018, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1101, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1115, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1129, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020123, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020923, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_021723, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_030923, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_031623, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_032323, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_040623, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_041323, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_042123, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_052423, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_060823, 380));
  aging_380_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_062323, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_317, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_324, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_407, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_414, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_421, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_428, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_505, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_512, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_526, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_602, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_616, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_714, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_728, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_811, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_825, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_908, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_920, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1004, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1018, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1101, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1115, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1129, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020123, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020923, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_021723, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_030923, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_031623, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_032323, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_040623, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_041323, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_042123, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_052423, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_060823, 380));
  aging_380_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_062323, 380));

  // wavelength = 390nm
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_317, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_324, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_407, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_414, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_421, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_428, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_505, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_512, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_526, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_602, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_616, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_714, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_728, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_811, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_825, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_908, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_920, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1004, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1018, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1101, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1115, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1129, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020123, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020923, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_021723, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_030923, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_031623, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_032323, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_040623, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_041323, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_042123, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_052423, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_060823, 390));
  aging_390_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_062323, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_317, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_324, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_407, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_414, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_421, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_428, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_505, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_512, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_526, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_602, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_616, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_714, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_728, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_811, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_825, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_908, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_920, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1004, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1018, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1101, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1115, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1129, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020123, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020923, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_021723, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_030923, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_031623, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_032323, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_040623, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_041323, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_042123, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_052423, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_060823, 390));
  aging_390_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_062323, 390));

  // wavelength = 400nm
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_317, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_324, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_407, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_414, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_421, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_428, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_505, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_512, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_526, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_602, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_616, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_714, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_728, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_811, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_825, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_908, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_920, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1004, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1018, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1101, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1115, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1129, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020123, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020923, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_021723, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_030923, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_031623, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_032323, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_040623, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_041323, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_042123, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_052423, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_060823, 400));
  aging_400_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_062323, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_317, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_324, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_407, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_414, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_421, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_428, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_505, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_512, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_526, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_602, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_616, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_714, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_728, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_811, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_825, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_908, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_920, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1004, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1018, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1101, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1115, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1129, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020123, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020923, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_021723, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_030923, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_031623, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_032323, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_040623, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_041323, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_042123, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_052423, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_060823, 400));
  aging_400_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_062323, 400));

  // wavelength = 410nm
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_317, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_324, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_407, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_414, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_421, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_428, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_505, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_512, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_526, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_602, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_616, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_714, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_728, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_811, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_825, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_908, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_920, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1004, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1018, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1101, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1115, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1129, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020123, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020923, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_021723, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_030923, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_031623, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_032323, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_040623, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_041323, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_042123, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_052423, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_060823, 410));
  aging_410_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_062323, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_317, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_324, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_407, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_414, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_421, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_428, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_505, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_512, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_526, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_602, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_616, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_714, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_728, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_811, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_825, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_908, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_920, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1004, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1018, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1101, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1115, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1129, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020123, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020923, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_021723, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_030923, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_031623, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_032323, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_040623, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_041323, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_042123, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_052423, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_060823, 410));
  aging_410_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_062323, 410));

  // wavelength = 420nm
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_317, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_324, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_407, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_414, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_421, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_428, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_505, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_512, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_526, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_602, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_616, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_714, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_728, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_811, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_825, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_908, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_920, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1004, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1018, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1101, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1115, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1129, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020123, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020923, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_021723, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_030923, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_031623, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_032323, 420));  
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_040623, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_041323, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_042123, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_052423, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_060823, 420));
  aging_420_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_062323, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_317, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_324, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_407, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_414, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_421, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_428, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_505, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_512, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_526, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_602, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_616, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_714, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_728, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_811, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_825, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_908, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_920, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1004, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1018, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1101, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1115, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1129, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020123, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020923, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_021723, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_030923, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_031623, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_032323, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_040623, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_041323, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_042123, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_052423, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_060823, 420));
  aging_420_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_062323, 420));

  // wavelength = 430nm
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_317, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_324, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_407, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_414, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_421, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_428, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_505, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_512, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_526, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_602, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_616, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_714, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_728, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_811, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_825, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_908, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_920, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1004, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1018, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1101, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1115, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1129, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020123, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020923, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_021723, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_030923, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_031623, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_032323, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_040623, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_041323, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_042123, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_052423, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_060823, 430));
  aging_430_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_062323, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_317, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_324, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_407, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_414, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_421, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_428, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_505, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_512, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_526, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_602, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_616, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_714, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_728, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_811, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_825, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_908, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_920, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1004, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1018, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1101, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1115, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1129, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020123, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020923, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_021723, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_030923, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_031623, 430));  
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_032323, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_040623, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_041323, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_042123, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_052423, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_060823, 430));
  aging_430_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_062323, 430));

  // wavelength = 440nm
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_317, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_324, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_407, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_414, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_421, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_428, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_505, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_512, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_526, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_602, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_616, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_714, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_728, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_811, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_825, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_908, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_920, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1004, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1018, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1101, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1115, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1129, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020123, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020923, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_021723, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_030923, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_031623, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_032323, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_040623, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_041323, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_042123, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_052423, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_060823, 440));
  aging_440_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_062323, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_317, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_324, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_407, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_414, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_421, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_428, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_505, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_512, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_526, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_602, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_616, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_714, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_728, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_811, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_825, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_908, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_920, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1004, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1018, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1101, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1115, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1129, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020123, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020923, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_021723, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_030923, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_031623, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_032323, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_040623, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_041323, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_042123, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_052423, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_060823, 440));
  aging_440_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_062323, 440));

  // wavelength = 450nm
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_317, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_324, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_407, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_414, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_421, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_428, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_505, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_512, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_526, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_602, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_616, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_714, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_728, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_811, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_825, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_908, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_920, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1004, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1018, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1101, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1115, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thin_1129, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020123, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_020923, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_021723, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_030923, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_031623, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_032323, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_040623, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_041323, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_042123, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_052423, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_060823, 450));
  aging_450_115101_thin_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thin_062323, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_317, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_324, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_407, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_414, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_421, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_428, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_505, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_512, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_526, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_602, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_616, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_714, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_728, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_811, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_825, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_908, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_920, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1004, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1018, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1101, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1115, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_mackenzie_thick_1129, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020123, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_020923, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_021723, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_030923, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_031623, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_032323, 450));  
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_040623, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_041323, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_042123, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_052423, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_060823, 450));
  aging_450_115101_thick_err.push_back(stdev_error(input_hists, avg_USPN115101_brian_thick_062323, 450));

  stab_115101_thin.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thin_920, 380));
  stab_115101_thin.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thin_920, 390));
  stab_115101_thin.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thin_920, 400));
  stab_115101_thin.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thin_920, 410));
  stab_115101_thin.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thin_920, 420));
  stab_115101_thin.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thin_920, 430));
  stab_115101_thin.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thin_920, 440));
  stab_115101_thin.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thin_920, 450));

  stab_115101_thick.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thick_920, 380));
  stab_115101_thick.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thick_920, 390));
  stab_115101_thick.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thick_920, 400));
  stab_115101_thick.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thick_920, 410));
  stab_115101_thick.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thick_920, 420));
  stab_115101_thick.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thick_920, 430));
  stab_115101_thick.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thick_920, 440));
  stab_115101_thick.push_back(stdev_error(input_hists, stab_USPN115101_mackenzie_thick_920, 450));

  // make TCanvases and TGraphErrors objects to hold aging plots
  gStyle->SetOptFit(11111);

  TCanvas* page[16];
  TCanvas* page_ce[4];
  TCanvas* stability[2];
  TGraphErrors* tg_115101_thin_frac_380;
  TGraphErrors* tg_115101_thick_frac_380;
  TGraphErrors* tg_115101_thin_frac_380_ce;
  TGraphErrors* tg_115101_thick_frac_380_ce;
  TGraphErrors* tg_115101_thin_frac_390;
  TGraphErrors* tg_115101_thick_frac_390;
  TGraphErrors* tg_115101_thin_frac_400;
  TGraphErrors* tg_115101_thick_frac_400;
  TGraphErrors* tg_115101_thin_frac_410;
  TGraphErrors* tg_115101_thick_frac_410;
  TGraphErrors* tg_115101_thin_frac_410_ce;
  TGraphErrors* tg_115101_thick_frac_410_ce;
  TGraphErrors* tg_115101_thin_frac_420;
  TGraphErrors* tg_115101_thick_frac_420;
  TGraphErrors* tg_115101_thin_frac_430;
  TGraphErrors* tg_115101_thick_frac_430;
  TGraphErrors* tg_115101_thin_frac_440;
  TGraphErrors* tg_115101_thick_frac_440;
  TGraphErrors* tg_115101_thin_frac_450;
  TGraphErrors* tg_115101_thick_frac_450;
  TGraph* tg_stab_115101_thin;
  TGraph* tg_stab_115101_thick;
  TH1F* hslopes1D_thin;
  TH1F* hslopes1D_thick;

  page[12] = new TCanvas("thin_frac_380","thin_frac_380",900,600);
  page[13] = new TCanvas("thick_frac_380","thick_frac_380",900,600);
  page[14] = new TCanvas("thin_frac_390","thin_frac_390",900,600);
  page[15] = new TCanvas("thick_frac_390","thick_frac_390",900,600);
  page[0] = new TCanvas("thin_frac_400","thin_frac_400",900,600);
  page[1] = new TCanvas("thick_frac_400","thick_frac_400",900,600);
  page[2] = new TCanvas("thin_frac_410","thin_frac_410",900,600);
  page[3] = new TCanvas("thick_frac_410","thick_frac_410",900,600);
  page[4] = new TCanvas("thin_frac_420","thin_frac_420",900,600);
  page[5] = new TCanvas("thick_frac_420","thick_frac_420",900,600);
  page[6] = new TCanvas("thin_frac_430","thin_frac_430",900,600);
  page[7] = new TCanvas("thick_frac_430","thick_frac_430",900,600);
  page[8] = new TCanvas("thin_frac_440","thin_frac_440",900,600);
  page[9] = new TCanvas("thick_frac_440","thick_frac_440",900,600);
  page[10] = new TCanvas("thin_frac_450","thin_frac_450",900,600);
  page[11] = new TCanvas("thick_frac_450","thick_frac_450",900,600);

  page_ce[0] = new TCanvas("thin_frac_380_ce","thin_frac_380_ce",900,600);
  page_ce[1] = new TCanvas("thick_frac_380_ce","thick_frac_380_ce",900,600);
  page_ce[2] = new TCanvas("thin_frac_410_ce","thin_frac_410_ce",900,600);
  page_ce[3] = new TCanvas("thick_frac_410_ce","thick_frac_410_ce",900,600);

  stability[0] = new TCanvas("thin_stability", "thin_stability",900,600);
  stability[1] = new TCanvas("thick_stability", "thick_stability",900,600);

  hslopes1D_thin = new TH1F("hslopes1D_thin", "hslopes1D_thin", 100, -1., 1.);
  hslopes1D_thin->SetTitle("USPN115101 Aging Slopes, Thin-side");
  hslopes1D_thin->GetXaxis()->SetTitle("Thin-side Slopes (Percent of R Lost per Year)");
  hslopes1D_thin->GetYaxis()->SetTitle("Counts");

  hslopes1D_thick = new TH1F("hslopes1D_thick", "hslopes1D_thick", 100, -1., 1.);
  hslopes1D_thick->SetTitle("USPN115101 Aging Slopes, Thick-side");
  hslopes1D_thick->GetXaxis()->SetTitle("Thick-side Slopes (Percent of R Lost per Year)");
  hslopes1D_thick->GetYaxis()->SetTitle("Counts");

  // need arrays for TGraphErrors

  // wavelength = 380nm
  float aging_380_115101_thin_a[npoints];
  float aging_380_115101_thick_a[npoints];
  float aging_380_115101_thin_err_a[npoints];
  float aging_380_115101_thick_err_a[npoints];
  std::copy(aging_380_115101_thin.begin(), aging_380_115101_thin.end(), aging_380_115101_thin_a);
  std::copy(aging_380_115101_thick.begin(), aging_380_115101_thick.end(), aging_380_115101_thick_a);
  std::copy(aging_380_115101_thin_err.begin(), aging_380_115101_thin_err.end(), aging_380_115101_thin_err_a);
  std::copy(aging_380_115101_thick_err.begin(), aging_380_115101_thick_err.end(), aging_380_115101_thick_err_a);
  
  //constant error array
  float aging_115101_err_a_ce[npoints];
  std::fill(aging_115101_err_a_ce, aging_115101_err_a_ce+npoints, 0.1);

  // wavelength = 390nm
  float aging_390_115101_thin_a[npoints];
  float aging_390_115101_thick_a[npoints];
  float aging_390_115101_thin_err_a[npoints];
  float aging_390_115101_thick_err_a[npoints];
  std::copy(aging_390_115101_thin.begin(), aging_390_115101_thin.end(), aging_390_115101_thin_a);
  std::copy(aging_390_115101_thick.begin(), aging_390_115101_thick.end(), aging_390_115101_thick_a);
  std::copy(aging_390_115101_thin_err.begin(), aging_390_115101_thin_err.end(), aging_390_115101_thin_err_a);
  std::copy(aging_390_115101_thick_err.begin(), aging_390_115101_thick_err.end(), aging_390_115101_thick_err_a);

  // wavelength = 400nm
  float aging_400_115101_thin_a[npoints];
  float aging_400_115101_thick_a[npoints];
  float aging_400_115101_thin_err_a[npoints];
  float aging_400_115101_thick_err_a[npoints];
  std::copy(aging_400_115101_thin.begin(), aging_400_115101_thin.end(), aging_400_115101_thin_a);
  std::copy(aging_400_115101_thick.begin(), aging_400_115101_thick.end(), aging_400_115101_thick_a);
  std::copy(aging_400_115101_thin_err.begin(), aging_400_115101_thin_err.end(), aging_400_115101_thin_err_a);
  std::copy(aging_400_115101_thick_err.begin(), aging_400_115101_thick_err.end(), aging_400_115101_thick_err_a);

  // wavelength = 410nm
  float aging_410_115101_thin_a[npoints];
  float aging_410_115101_thick_a[npoints];
  float aging_410_115101_thin_err_a[npoints];
  float aging_410_115101_thick_err_a[npoints];
  std::copy(aging_410_115101_thin.begin(), aging_410_115101_thin.end(), aging_410_115101_thin_a);
  std::copy(aging_410_115101_thick.begin(), aging_410_115101_thick.end(), aging_410_115101_thick_a);
  std::copy(aging_410_115101_thin_err.begin(), aging_410_115101_thin_err.end(), aging_410_115101_thin_err_a);
  std::copy(aging_410_115101_thick_err.begin(), aging_410_115101_thick_err.end(), aging_410_115101_thick_err_a);

  // wavelength = 420nm
  float aging_420_115101_thin_a[npoints];
  float aging_420_115101_thick_a[npoints];
  float aging_420_115101_thin_err_a[npoints];
  float aging_420_115101_thick_err_a[npoints];
  std::copy(aging_420_115101_thin.begin(), aging_420_115101_thin.end(), aging_420_115101_thin_a);
  std::copy(aging_420_115101_thick.begin(), aging_420_115101_thick.end(), aging_420_115101_thick_a);
  std::copy(aging_420_115101_thin_err.begin(), aging_420_115101_thin_err.end(), aging_420_115101_thin_err_a);
  std::copy(aging_420_115101_thick_err.begin(), aging_420_115101_thick_err.end(), aging_420_115101_thick_err_a);

  // wavelength = 430nm
  float aging_430_115101_thin_a[npoints];
  float aging_430_115101_thick_a[npoints];
  float aging_430_115101_thin_err_a[npoints];
  float aging_430_115101_thick_err_a[npoints];
  std::copy(aging_430_115101_thin.begin(), aging_430_115101_thin.end(), aging_430_115101_thin_a);
  std::copy(aging_430_115101_thick.begin(), aging_430_115101_thick.end(), aging_430_115101_thick_a);
  std::copy(aging_430_115101_thin_err.begin(), aging_430_115101_thin_err.end(), aging_430_115101_thin_err_a);
  std::copy(aging_430_115101_thick_err.begin(), aging_430_115101_thick_err.end(), aging_430_115101_thick_err_a);

  // wavelength = 440nm
  float aging_440_115101_thin_a[npoints];
  float aging_440_115101_thick_a[npoints];
  float aging_440_115101_thin_err_a[npoints];
  float aging_440_115101_thick_err_a[npoints];
  std::copy(aging_440_115101_thin.begin(), aging_440_115101_thin.end(), aging_440_115101_thin_a);
  std::copy(aging_440_115101_thick.begin(), aging_440_115101_thick.end(), aging_440_115101_thick_a);
  std::copy(aging_440_115101_thin_err.begin(), aging_440_115101_thin_err.end(), aging_440_115101_thin_err_a);
  std::copy(aging_440_115101_thick_err.begin(), aging_440_115101_thick_err.end(), aging_440_115101_thick_err_a);

  // wavelength = 450nm
  float aging_450_115101_thin_a[npoints];
  float aging_450_115101_thick_a[npoints];
  float aging_450_115101_thin_err_a[npoints];
  float aging_450_115101_thick_err_a[npoints];
  std::copy(aging_450_115101_thin.begin(), aging_450_115101_thin.end(), aging_450_115101_thin_a);
  std::copy(aging_450_115101_thick.begin(), aging_450_115101_thick.end(), aging_450_115101_thick_a);
  std::copy(aging_450_115101_thin_err.begin(), aging_450_115101_thin_err.end(), aging_450_115101_thin_err_a);
  std::copy(aging_450_115101_thick_err.begin(), aging_450_115101_thick_err.end(), aging_450_115101_thick_err_a);

  float stab_115101_thin_a[8];
  float stab_115101_thick_a[8];
  std::copy(stab_115101_thin.begin(), stab_115101_thin.end(), stab_115101_thin_a);
  std::copy(stab_115101_thick.begin(), stab_115101_thick.end(), stab_115101_thick_a);

  std::cout << "thin aging errors: ";
  for(int i = 0;i < npoints;i++){
    if(i < 17){
      std::cout << aging_400_115101_thin_err_a[i] << ", ";
    }
    if(i == 17){
      std::cout << aging_400_115101_thin_err_a[i] << std::endl;
    }
  }
  std::vector<float> thin_nm_slopes;
  std::vector<float> thick_nm_slopes;
  std::vector<float> thin_nm_slopes_err;
  std::vector<float> thick_nm_slopes_err;

  page[12]->cd();
  tg_115101_thin_frac_380 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_380_115101_thin_a, aging_empty_xbins_err, aging_380_115101_thin_err_a);
  tg_115101_thin_frac_380->SetTitle("Average Thin-side R at 380nm of 115101 vs Time");
  tg_115101_thin_frac_380->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thin_frac_380->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thin_frac_380->SetMarkerStyle(20);
  tg_115101_thin_frac_380->Fit("pol1");
  tg_115101_thin_frac_380->Draw("AP");
  TF1 *fit12 = (TF1*)tg_115101_thin_frac_380->GetListOfFunctions()->FindObject("pol1");
  thin_nm_slopes.push_back(fit12->GetParameter(1));
  thin_nm_slopes_err.push_back(fit12->GetParError(1));
  page[12]->Update();
  ofile->cd();
  page[12]->Write();
  page[12]->Close();

  page[13]->cd();
  tg_115101_thick_frac_380 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_380_115101_thick_a, aging_empty_xbins_err, aging_380_115101_thick_err_a);
  tg_115101_thick_frac_380->SetTitle("Average Thick-side R at 380nm of 115101 vs Time");
  tg_115101_thick_frac_380->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thick_frac_380->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thick_frac_380->SetMarkerStyle(20);
  tg_115101_thick_frac_380->Fit("pol1");
  tg_115101_thick_frac_380->Draw("AP");
  TF1 *fit13 = (TF1*)tg_115101_thick_frac_380->GetListOfFunctions()->FindObject("pol1");
  thick_nm_slopes.push_back(fit13->GetParameter(1));
  thick_nm_slopes_err.push_back(fit13->GetParError(1));
  page[13]->Update();
  ofile->cd();
  page[13]->Write();
  page[13]->Close();

  page_ce[0]->cd();
  tg_115101_thin_frac_380_ce = new TGraphErrors(npoints, aging_timefrac_xbins, aging_380_115101_thin_a, aging_empty_xbins_err, aging_115101_err_a_ce);
  tg_115101_thin_frac_380_ce->SetTitle("Average Thin-side R at 380nm of 115101 vs Time");
  tg_115101_thin_frac_380_ce->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thin_frac_380_ce->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thin_frac_380_ce->SetMarkerStyle(20);
  tg_115101_thin_frac_380_ce->Fit("pol1");
  tg_115101_thin_frac_380_ce->Draw("AP");
  page_ce[0]->Update();
  ofile->cd();
  page_ce[0]->Write();
  page_ce[0]->Close();

  page_ce[1]->cd();
  tg_115101_thick_frac_380_ce = new TGraphErrors(npoints, aging_timefrac_xbins, aging_380_115101_thick_a, aging_empty_xbins_err, aging_115101_err_a_ce);
  tg_115101_thick_frac_380_ce->SetTitle("Average Thick-side R at 380nm of 115101 vs Time");
  tg_115101_thick_frac_380_ce->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thick_frac_380_ce->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thick_frac_380_ce->SetMarkerStyle(20);
  tg_115101_thick_frac_380_ce->Fit("pol1");
  tg_115101_thick_frac_380_ce->Draw("AP");
  page_ce[1]->Update();
  ofile->cd();
  page_ce[1]->Write();
  page_ce[1]->Close();

  page[14]->cd();
  tg_115101_thin_frac_390 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_390_115101_thin_a, aging_empty_xbins_err, aging_390_115101_thin_err_a);
  tg_115101_thin_frac_390->SetTitle("Average Thin-side R at 390nm of 115101 vs Time");
  tg_115101_thin_frac_390->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thin_frac_390->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thin_frac_390->SetMarkerStyle(20);
  tg_115101_thin_frac_390->Fit("pol1");
  tg_115101_thin_frac_390->Draw("AP");
  TF1 *fit14 = (TF1*)tg_115101_thin_frac_390->GetListOfFunctions()->FindObject("pol1");
  thin_nm_slopes.push_back(fit14->GetParameter(1));
  thin_nm_slopes_err.push_back(fit14->GetParError(1));
  page[14]->Update();
  ofile->cd();
  page[14]->Write();
  page[14]->Close();

  page[15]->cd();
  tg_115101_thick_frac_390 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_390_115101_thick_a, aging_empty_xbins_err, aging_390_115101_thick_err_a);
  tg_115101_thick_frac_390->SetTitle("Average Thick-side R at 390nm of 115101 vs Time");
  tg_115101_thick_frac_390->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thick_frac_390->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thick_frac_390->SetMarkerStyle(20);
  tg_115101_thick_frac_390->Fit("pol1");
  tg_115101_thick_frac_390->Draw("AP");
  TF1 *fit15 = (TF1*)tg_115101_thick_frac_390->GetListOfFunctions()->FindObject("pol1");
  thick_nm_slopes.push_back(fit15->GetParameter(1));
  thick_nm_slopes_err.push_back(fit15->GetParError(1));
  page[15]->Update();
  ofile->cd();
  page[15]->Write();
  page[15]->Close();

  page[0]->cd();
  tg_115101_thin_frac_400 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_400_115101_thin_a, aging_empty_xbins_err, aging_400_115101_thin_err_a);
  tg_115101_thin_frac_400->SetTitle("Average Thin-side R at 400nm of 115101 vs Time");
  tg_115101_thin_frac_400->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thin_frac_400->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thin_frac_400->SetMarkerStyle(20);
  tg_115101_thin_frac_400->Fit("pol1");
  tg_115101_thin_frac_400->Draw("AP");
  TF1 *fit0 = (TF1*)tg_115101_thin_frac_400->GetListOfFunctions()->FindObject("pol1");
  thin_nm_slopes.push_back(fit0->GetParameter(1));
  thin_nm_slopes_err.push_back(fit0->GetParError(1));
  page[0]->Update();
  ofile->cd();
  page[0]->Write();
  page[0]->Close();

  page[1]->cd();
  tg_115101_thick_frac_400 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_400_115101_thick_a, aging_empty_xbins_err, aging_400_115101_thick_err_a);
  tg_115101_thick_frac_400->SetTitle("Average Thick-side R at 400nm of 115101 vs Time");
  tg_115101_thick_frac_400->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thick_frac_400->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thick_frac_400->SetMarkerStyle(20);
  tg_115101_thick_frac_400->Fit("pol1");
  tg_115101_thick_frac_400->Draw("AP");
  TF1 *fit1 = (TF1*)tg_115101_thick_frac_400->GetListOfFunctions()->FindObject("pol1");
  thick_nm_slopes.push_back(fit1->GetParameter(1));
  thick_nm_slopes_err.push_back(fit1->GetParError(1));
  page[1]->Update();
  ofile->cd();
  page[1]->Write();
  page[1]->Close();

  page[2]->cd();
  tg_115101_thin_frac_410 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_410_115101_thin_a, aging_empty_xbins_err, aging_410_115101_thin_err_a);
  tg_115101_thin_frac_410->SetTitle("Average Thin-side R at 410nm of 115101 vs Time");
  tg_115101_thin_frac_410->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thin_frac_410->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thin_frac_410->SetMarkerStyle(20);
  tg_115101_thin_frac_410->Fit("pol1");
  tg_115101_thin_frac_410->Draw("AP");
  TF1 *fit2 = (TF1*)tg_115101_thin_frac_410->GetListOfFunctions()->FindObject("pol1");
  thin_nm_slopes.push_back(fit2->GetParameter(1));
  thin_nm_slopes_err.push_back(fit2->GetParError(1));
  page[2]->Update();
  ofile->cd();
  page[2]->Write();
  page[2]->Close();

  page[3]->cd();
  tg_115101_thick_frac_410 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_410_115101_thick_a, aging_empty_xbins_err, aging_410_115101_thick_err_a);
  tg_115101_thick_frac_410->SetTitle("Average Thick-side R at 410nm of 115101 vs Time");
  tg_115101_thick_frac_410->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thick_frac_410->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thick_frac_410->SetMarkerStyle(20);
  tg_115101_thick_frac_410->Fit("pol1");
  tg_115101_thick_frac_410->Draw("AP");
  TF1 *fit3 = (TF1*)tg_115101_thick_frac_410->GetListOfFunctions()->FindObject("pol1");
  thick_nm_slopes.push_back(fit3->GetParameter(1));
  thick_nm_slopes_err.push_back(fit3->GetParError(1));
  page[3]->Update();
  ofile->cd();
  page[3]->Write();
  page[3]->Close();

  page_ce[2]->cd();
  tg_115101_thin_frac_410_ce = new TGraphErrors(npoints, aging_timefrac_xbins, aging_410_115101_thin_a, aging_empty_xbins_err, aging_115101_err_a_ce);
  tg_115101_thin_frac_410_ce->SetTitle("Average Thin-side R at 410nm of 115101 vs Time");
  tg_115101_thin_frac_410_ce->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thin_frac_410_ce->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thin_frac_410_ce->SetMarkerStyle(20);
  tg_115101_thin_frac_410_ce->Fit("pol1");
  tg_115101_thin_frac_410_ce->Draw("AP");
  page_ce[2]->Update();
  ofile->cd();
  page_ce[2]->Write();
  page_ce[2]->Close();

  page_ce[3]->cd();
  tg_115101_thick_frac_410_ce = new TGraphErrors(npoints, aging_timefrac_xbins, aging_410_115101_thick_a, aging_empty_xbins_err, aging_115101_err_a_ce);
  tg_115101_thick_frac_410_ce->SetTitle("Average Thick-side R at 410nm of 115101 vs Time");
  tg_115101_thick_frac_410_ce->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thick_frac_410_ce->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thick_frac_410_ce->SetMarkerStyle(20);
  tg_115101_thick_frac_410_ce->Fit("pol1");
  tg_115101_thick_frac_410_ce->Draw("AP");
  page_ce[3]->Update();
  ofile->cd();
  page_ce[3]->Write();
  page_ce[3]->Close();

  page[4]->cd();
  tg_115101_thin_frac_420 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_420_115101_thin_a, aging_empty_xbins_err, aging_420_115101_thin_err_a);
  tg_115101_thin_frac_420->SetTitle("Average Thin-side R at 420nm of 115101 vs Time");
  tg_115101_thin_frac_420->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thin_frac_420->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thin_frac_420->SetMarkerStyle(20);
  tg_115101_thin_frac_420->Fit("pol1");
  tg_115101_thin_frac_420->Draw("AP");
  TF1 *fit4 = (TF1*)tg_115101_thin_frac_420->GetListOfFunctions()->FindObject("pol1");
  thin_nm_slopes.push_back(fit4->GetParameter(1));
  thin_nm_slopes_err.push_back(fit4->GetParError(1));
  page[4]->Update();
  ofile->cd();
  page[4]->Write();
  page[4]->Close();

  page[5]->cd();
  tg_115101_thick_frac_420 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_420_115101_thick_a, aging_empty_xbins_err, aging_420_115101_thick_err_a);
  tg_115101_thick_frac_420->SetTitle("Average Thick-side R at 420nm of 115101 vs Time");
  tg_115101_thick_frac_420->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thick_frac_420->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thick_frac_420->SetMarkerStyle(20);
  tg_115101_thick_frac_420->Fit("pol1");
  tg_115101_thick_frac_420->Draw("AP");
  TF1 *fit5 = (TF1*)tg_115101_thick_frac_420->GetListOfFunctions()->FindObject("pol1");
  thick_nm_slopes.push_back(fit5->GetParameter(1));
  thick_nm_slopes_err.push_back(fit5->GetParError(1));
  page[5]->Update();
  ofile->cd();
  page[5]->Write();
  page[5]->Close();

  page[6]->cd();
  tg_115101_thin_frac_430 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_430_115101_thin_a, aging_empty_xbins_err, aging_430_115101_thin_err_a);
  tg_115101_thin_frac_430->SetTitle("Average Thin-side R at 430nm of 115101 vs Time");
  tg_115101_thin_frac_430->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thin_frac_430->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thin_frac_430->SetMarkerStyle(20);
  tg_115101_thin_frac_430->Fit("pol1");
  tg_115101_thin_frac_430->Draw("AP");
  TF1 *fit6 = (TF1*)tg_115101_thin_frac_430->GetListOfFunctions()->FindObject("pol1");
  thin_nm_slopes.push_back(fit6->GetParameter(1));
  thin_nm_slopes_err.push_back(fit6->GetParError(1));
  page[6]->Update();
  ofile->cd();
  page[6]->Write();
  page[6]->Close();

  page[7]->cd();
  tg_115101_thick_frac_430 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_430_115101_thick_a, aging_empty_xbins_err, aging_430_115101_thick_err_a);
  tg_115101_thick_frac_430->SetTitle("Average Thick-side R at 430nm of 115101 vs Time");
  tg_115101_thick_frac_430->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thick_frac_430->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thick_frac_430->SetMarkerStyle(20);
  tg_115101_thick_frac_430->Fit("pol1");
  tg_115101_thick_frac_430->Draw("AP");
  TF1 *fit7 = (TF1*)tg_115101_thick_frac_430->GetListOfFunctions()->FindObject("pol1");
  thick_nm_slopes.push_back(fit7->GetParameter(1));
  thick_nm_slopes_err.push_back(fit7->GetParError(1));
  page[7]->Update();
  ofile->cd();
  page[7]->Write();
  page[7]->Close();

  page[8]->cd();
  tg_115101_thin_frac_440 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_440_115101_thin_a, aging_empty_xbins_err, aging_440_115101_thin_err_a);
  tg_115101_thin_frac_440->SetTitle("Average Thin-side R at 440nm of 115101 vs Time");
  tg_115101_thin_frac_440->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thin_frac_440->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thin_frac_440->SetMarkerStyle(20);
  tg_115101_thin_frac_440->Fit("pol1");
  tg_115101_thin_frac_440->Draw("AP");
  TF1 *fit8 = (TF1*)tg_115101_thin_frac_440->GetListOfFunctions()->FindObject("pol1");
  thin_nm_slopes.push_back(fit8->GetParameter(1));
  thin_nm_slopes_err.push_back(fit8->GetParError(1));
  page[8]->Update();
  ofile->cd();
  page[8]->Write();
  page[8]->Close();

  page[9]->cd();
  tg_115101_thick_frac_440 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_440_115101_thick_a, aging_empty_xbins_err, aging_440_115101_thick_err_a);
  tg_115101_thick_frac_440->SetTitle("Average Thick-side R at 440nm of 115101 vs Time");
  tg_115101_thick_frac_440->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thick_frac_440->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thick_frac_440->SetMarkerStyle(20);
  tg_115101_thick_frac_440->Fit("pol1");
  tg_115101_thick_frac_440->Draw("AP");
  TF1 *fit9 = (TF1*)tg_115101_thick_frac_440->GetListOfFunctions()->FindObject("pol1");
  thick_nm_slopes.push_back(fit9->GetParameter(1));
  thick_nm_slopes_err.push_back(fit9->GetParError(1));
  page[9]->Update();
  ofile->cd();
  page[9]->Write();
  page[9]->Close();

  page[10]->cd();
  tg_115101_thin_frac_450 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_450_115101_thin_a, aging_empty_xbins_err, aging_450_115101_thin_err_a);
  tg_115101_thin_frac_450->SetTitle("Average Thin-side R at 450nm of 115101 vs Time");
  tg_115101_thin_frac_450->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thin_frac_450->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thin_frac_450->SetMarkerStyle(20);
  tg_115101_thin_frac_450->Fit("pol1");
  tg_115101_thin_frac_450->Draw("AP");
  TF1 *fit10 = (TF1*)tg_115101_thin_frac_450->GetListOfFunctions()->FindObject("pol1");
  thin_nm_slopes.push_back(fit10->GetParameter(1));
  thin_nm_slopes_err.push_back(fit10->GetParError(1));
  page[10]->Update();
  ofile->cd();
  page[10]->Write();
  page[10]->Close();

  page[11]->cd();
  tg_115101_thick_frac_450 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_450_115101_thick_a, aging_empty_xbins_err, aging_450_115101_thick_err_a);
  tg_115101_thick_frac_450->SetTitle("Average Thick-side R at 450nm of 115101 vs Time");
  tg_115101_thick_frac_450->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  tg_115101_thick_frac_450->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_115101_thick_frac_450->SetMarkerStyle(20);
  tg_115101_thick_frac_450->Fit("pol1");
  tg_115101_thick_frac_450->Draw("AP");
  TF1 *fit11 = (TF1*)tg_115101_thick_frac_450->GetListOfFunctions()->FindObject("pol1");
  thick_nm_slopes.push_back(fit11->GetParameter(1));
  thick_nm_slopes_err.push_back(fit11->GetParError(1));
  page[11]->Update();
  ofile->cd();
  page[11]->Write();
  page[11]->Close();

  // slopes as a function of wavelength

  TCanvas* slopes = new TCanvas("aging_slope_vs_wavelength","aging_slope_vs_wavelength", 900, 600);
  TCanvas* slopes_err = new TCanvas("aging_slope_vs_wavelength_err","aging_slope_vs_wavelength_err", 900, 600);
  TCanvas* slopes_err2 = new TCanvas("aging_slope_vs_wavelength_err2","aging_slope_vs_wavelength_err2", 900, 600);

  float thin_nm_slopes_a[8];
  float thick_nm_slopes_a[8];
  std::copy(thin_nm_slopes.begin(), thin_nm_slopes.end(), thin_nm_slopes_a);
  std::copy(thick_nm_slopes.begin(), thick_nm_slopes.end(), thick_nm_slopes_a);
  float thin_nm_slopes_err_a[8];
  float thick_nm_slopes_err_a[8];
  std::copy(thin_nm_slopes_err.begin(), thin_nm_slopes_err.end(), thin_nm_slopes_err_a);
  std::copy(thick_nm_slopes_err.begin(), thick_nm_slopes_err.end(), thick_nm_slopes_err_a);
  float wavelengths[8] = { 380., 390., 400., 410., 420., 430., 440., 450. };
  float empty_wl_err[8] = { 0., 0., 0., 0., 0., 0., 0., 0. };

  for(int i = 0; i < 8; i++){
    hslopes1D_thin->Fill(thin_nm_slopes_a[i]);
    hslopes1D_thick->Fill(thick_nm_slopes_a[i]);
  }
  hslopes1D_thin->Write();
  hslopes1D_thick->Write();

  slopes->cd();
  TGraphErrors* tg_thin_nm_slopes = new TGraphErrors(8, wavelengths, thin_nm_slopes_a);
  tg_thin_nm_slopes->SetMarkerStyle(24);
  TGraphErrors* tg_thick_nm_slopes = new TGraphErrors(8, wavelengths, thick_nm_slopes_a);
  tg_thick_nm_slopes->SetTitle("Aging Slopes vs Wavelength, 115101 Coupons");
  tg_thick_nm_slopes->GetXaxis()->SetTitle("Wavelength (nm)");
  tg_thick_nm_slopes->GetYaxis()->SetTitle("Aging Slope (R% Loss/Year)");
  tg_thick_nm_slopes->GetYaxis()->SetRangeUser(-1.2,0.8);
  tg_thick_nm_slopes->SetMarkerStyle(20);
  tg_thick_nm_slopes->Draw("AP");
  tg_thin_nm_slopes->Draw("P");
  slopes->Update();
  ofile->cd();
  slopes->Write();
  slopes->Close();

  slopes_err->cd();
  TGraphErrors* tg_thin_nm_slopes_err = new TGraphErrors(8, wavelengths,thin_nm_slopes_a, empty_wl_err, thin_nm_slopes_err_a);
  tg_thin_nm_slopes_err->SetMarkerStyle(24);
  TGraphErrors* tg_thick_nm_slopes_err = new TGraphErrors(8, wavelengths, thick_nm_slopes_a, empty_wl_err, thick_nm_slopes_err_a);
  tg_thick_nm_slopes_err->SetTitle("Aging Slopes vs Wavelength, USPN115101 Coupons");
  tg_thick_nm_slopes_err->GetXaxis()->SetTitle("Wavelength (nm)");
  tg_thick_nm_slopes_err->GetYaxis()->SetTitle("Aging Slope (R% Loss/Year)");
  tg_thick_nm_slopes_err->GetYaxis()->SetRangeUser(-0.8,0.2);
  tg_thick_nm_slopes_err->SetMarkerStyle(20);
  //tg_thick_nm_slopes_err->Fit("pol0");
  tg_thick_nm_slopes_err->Draw("AP");
  //tg_thin_nm_slopes_err->Fit("pol0");
  //tg_thin_nm_slopes_err->GetFunction("pol0")->SetLineColor(kBlue);
  tg_thin_nm_slopes_err->Draw("P");

  TLegend* slope_leg = new TLegend(0.12, 0.13, 0.4, 0.2);
  slope_leg->AddEntry(tg_thick_nm_slopes_err,"Thick Side Slopes", "lep");
  slope_leg->AddEntry(tg_thin_nm_slopes_err,"Thin Side Slopes", "lep");
  slope_leg->Draw();

  slopes_err->Update();
  ofile->cd();
  slopes_err->Write();
  slopes_err->Close();

  slopes_err2->cd(); 
  tg_thin_nm_slopes_err->SetTitle("Aging Slopes vs Wavelength, USPN115101 Coupons");
  tg_thin_nm_slopes_err->GetXaxis()->SetTitle("Wavelength (nm)");
  tg_thin_nm_slopes_err->GetYaxis()->SetTitle("Aging Slope (R% Loss/Year)");
  tg_thin_nm_slopes_err->GetYaxis()->SetRangeUser(-0.8,0.2);
  tg_thin_nm_slopes_err->SetMarkerStyle(24); 
  //tg_thin_nm_slopes_err->Fit("pol0");
  tg_thin_nm_slopes_err->Draw("AP");
  slopes_err2->Update();
  ofile->cd();
  slopes_err2->Write();
  slopes_err2->Close();

  stability[0]->cd();
  tg_stab_115101_thin = new TGraph(8, wavelengths, stab_115101_thin_a);
  tg_stab_115101_thin->SetMarkerStyle(20);
  tg_stab_115101_thin->SetTitle("Standard Deviation of Repeated Thin-Side Measurements");
  tg_stab_115101_thin->GetXaxis()->SetTitle("Wavelength (nm)");
  tg_stab_115101_thin->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  tg_stab_115101_thin->Draw("AP");
  stability[0]->Update();
  ofile->cd();
  stability[0]->Write();
  stability[0]->Close();

  stability[1]->cd();
  tg_stab_115101_thick = new TGraph(8, wavelengths, stab_115101_thick_a);
  tg_stab_115101_thick->SetMarkerStyle(20);
  tg_stab_115101_thick->SetTitle("Standard Deviation of Repeated Thick-Side Measurements");
  tg_stab_115101_thick->GetXaxis()->SetTitle("Wavelength (nm)");
  tg_stab_115101_thick->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  tg_stab_115101_thick->Draw("AP");
  stability[1]->Update();
  ofile->cd();
  stability[1]->Write();
  stability[1]->Close();

  // Close the output file
  ofile->Close();
}
