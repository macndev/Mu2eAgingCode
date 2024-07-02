// this file takes the extrusion data histograms as an input and will do some
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

void extrusion_analysis(const char* filename){

  TFile* file = TFile::Open(Form("%s", filename));
  if (file->IsZombie()) {
    std::cout << "Problem opening file " << filename << std::endl;
    return;
  }

  //Create output file. If it already exists, recreate it.
  TFile *ofile = new TFile("output_extrusion_analysis_aging_0814_no0524.root","RECREATE");

  // define average and difference histograms for each sample
  TH1F* ext_thin_poly_906_avg  = new TH1F("ext_thick_poly_906_avg", "ext_thick_poly_906_avg",43,360.,780.);
  TH1F* ext_thin_clad_906_avg  = new TH1F("ext_thick_clad_906_avg", "ext_thick_clad_906_avg",43,360.,780.);

  TH1F* ext_thin_poly_913_avg  = new TH1F("ext_thick_poly_913_avg", "ext_thick_poly_913_avg",43,360.,780.);
  TH1F* ext_thin_clad_913_avg  = new TH1F("ext_thick_clad_913_avg", "ext_thick_clad_913_avg",43,360.,780.);

  TH1F* ext_thin_poly_920_avg  = new TH1F("ext_thick_poly_920_avg", "ext_thick_poly_920_avg",43,360.,780.);
  TH1F* ext_thin_clad_920_avg  = new TH1F("ext_thick_clad_920_avg", "ext_thick_clad_920_avg",43,360.,780.);

  TH1F* ext_thin_poly_927_avg  = new TH1F("ext_thick_poly_927_avg", "ext_thick_poly_927_avg",43,360.,780.);
  TH1F* ext_thin_clad_927_avg  = new TH1F("ext_thick_clad_927_avg", "ext_thick_clad_927_avg",43,360.,780.);

  TH1F* ext_thin_poly_1004_avg  = new TH1F("ext_thick_poly_1004_avg", "ext_thick_poly_1004_avg",43,360.,780.);
  TH1F* ext_thin_clad_1004_avg  = new TH1F("ext_thick_clad_1004_avg", "ext_thick_clad_1004_avg",43,360.,780.);

  TH1F* ext_thin_poly_1011_avg  = new TH1F("ext_thick_poly_1011_avg", "ext_thick_poly_1011_avg",43,360.,780.);
  TH1F* ext_thin_clad_1011_avg  = new TH1F("ext_thick_clad_1011_avg", "ext_thick_clad_1011_avg",43,360.,780.);

  TH1F* ext_thin_poly_1018_avg  = new TH1F("ext_thick_poly_1018_avg", "ext_thick_poly_1018_avg",43,360.,780.);
  TH1F* ext_thin_clad_1018_avg  = new TH1F("ext_thick_clad_1018_avg", "ext_thick_clad_1018_avg",43,360.,780.);

  TH1F* ext_thin_poly_1025_avg  = new TH1F("ext_thick_poly_1025_avg", "ext_thick_poly_1025_avg",43,360.,780.);
  TH1F* ext_thin_clad_1025_avg  = new TH1F("ext_thick_clad_1025_avg", "ext_thick_clad_1025_avg",43,360.,780.);

  TH1F* ext_thin_poly_1101_avg  = new TH1F("ext_thick_poly_1101_avg", "ext_thick_poly_1101_avg",43,360.,780.);
  TH1F* ext_thin_clad_1101_avg  = new TH1F("ext_thick_clad_1101_avg", "ext_thick_clad_1101_avg",43,360.,780.);

  TH1F* ext_thin_poly_1108_avg  = new TH1F("ext_thick_poly_1108_avg", "ext_thick_poly_1108_avg",43,360.,780.);
  TH1F* ext_thin_clad_1108_avg  = new TH1F("ext_thick_clad_1108_avg", "ext_thick_clad_1108_avg",43,360.,780.);

  TH1F* ext_thin_poly_1115_avg  = new TH1F("ext_thick_poly_1115_avg", "ext_thick_poly_1115_avg",43,360.,780.);
  TH1F* ext_thin_clad_1115_avg  = new TH1F("ext_thick_clad_1115_avg", "ext_thick_clad_1115_avg",43,360.,780.);

  TH1F* ext_thin_poly_1122_avg  = new TH1F("ext_thick_poly_1122_avg", "ext_thick_poly_1122_avg",43,360.,780.);
  TH1F* ext_thin_clad_1122_avg  = new TH1F("ext_thick_clad_1122_avg", "ext_thick_clad_1122_avg",43,360.,780.);

  TH1F* ext_thin_poly_1129_avg  = new TH1F("ext_thick_poly_1129_avg", "ext_thick_poly_1129_avg",43,360.,780.);
  TH1F* ext_thin_clad_1129_avg  = new TH1F("ext_thick_clad_1129_avg", "ext_thick_clad_1129_avg",43,360.,780.);

  TH1F* ext_thin_poly_1206_avg  = new TH1F("ext_thick_poly_1206_avg", "ext_thick_poly_1206_avg",43,360.,780.);
  TH1F* ext_thin_clad_1206_avg  = new TH1F("ext_thick_clad_1206_avg", "ext_thick_clad_1206_avg",43,360.,780.);

  TH1F* ext_thin_poly_020123_avg  = new TH1F("ext_thick_poly_020123_avg", "ext_thick_poly_020123_avg",43,360.,780.);
  TH1F* ext_thin_clad_020123_avg  = new TH1F("ext_thick_clad_020123_avg", "ext_thick_clad_020123_avg",43,360.,780.);

  TH1F* ext_thin_poly_020923_avg  = new TH1F("ext_thick_poly_020923_avg", "ext_thick_poly_020923_avg",43,360.,780.);
  TH1F* ext_thin_clad_020923_avg  = new TH1F("ext_thick_clad_020923_avg", "ext_thick_clad_020923_avg",43,360.,780.);

  TH1F* ext_thin_poly_021723_avg  = new TH1F("ext_thick_poly_021723_avg", "ext_thick_poly_021723_avg",43,360.,780.);
  TH1F* ext_thin_clad_021723_avg  = new TH1F("ext_thick_clad_021723_avg", "ext_thick_clad_021723_avg",43,360.,780.);

  TH1F* ext_thin_poly_030923_avg  = new TH1F("ext_thick_poly_030923_avg", "ext_thick_poly_030923_avg",43,360.,780.);
  TH1F* ext_thin_clad_030923_avg  = new TH1F("ext_thick_clad_030923_avg", "ext_thick_clad_030923_avg",43,360.,780.);

  TH1F* ext_thin_poly_031623_avg  = new TH1F("ext_thick_poly_031623_avg", "ext_thick_poly_031623_avg",43,360.,780.);
  TH1F* ext_thin_clad_031623_avg  = new TH1F("ext_thick_clad_031623_avg", "ext_thick_clad_031623_avg",43,360.,780.);

  TH1F* ext_thin_poly_032323_avg  = new TH1F("ext_thick_poly_032323_avg", "ext_thick_poly_032323_avg",43,360.,780.);
  TH1F* ext_thin_clad_032323_avg  = new TH1F("ext_thick_clad_032323_avg", "ext_thick_clad_032323_avg",43,360.,780.);

  TH1F* ext_thin_poly_040623_avg  = new TH1F("ext_thick_poly_040623_avg", "ext_thick_poly_040623_avg",43,360.,780.);
  TH1F* ext_thin_clad_040623_avg  = new TH1F("ext_thick_clad_040623_avg", "ext_thick_clad_040623_avg",43,360.,780.);

  TH1F* ext_thin_poly_041323_avg  = new TH1F("ext_thick_poly_041323_avg", "ext_thick_poly_041323_avg",43,360.,780.);
  TH1F* ext_thin_clad_041323_avg  = new TH1F("ext_thick_clad_041323_avg", "ext_thick_clad_041323_avg",43,360.,780.);

  TH1F* ext_thin_poly_042123_avg  = new TH1F("ext_thick_poly_042123_avg", "ext_thick_poly_042123_avg",43,360.,780.);
  TH1F* ext_thin_clad_042123_avg  = new TH1F("ext_thick_clad_042123_avg", "ext_thick_clad_042123_avg",43,360.,780.);

  TH1F* ext_thin_clad_052423_avg  = new TH1F("ext_thick_clad_052423_avg", "ext_thick_clad_052423_avg",43,360.,780.);
  TH1F* ext_thin_clad_060823_avg  = new TH1F("ext_thick_clad_060823_avg", "ext_thick_clad_060823_avg",43,360.,780.);
  TH1F* ext_thin_clad_062323_avg  = new TH1F("ext_thick_clad_062323_avg", "ext_thick_clad_062323_avg",43,360.,780.);

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
  std::vector<int> avg_ext_thin_poly_906;
  std::vector<int> avg_ext_thin_clad_906;
  std::vector<int> avg_ext_thin_poly_913;
  std::vector<int> avg_ext_thin_clad_913;
  std::vector<int> avg_ext_thin_poly_920;
  std::vector<int> avg_ext_thin_clad_920;
  std::vector<int> avg_ext_thin_poly_927;
  std::vector<int> avg_ext_thin_clad_927;
  std::vector<int> avg_ext_thin_poly_1004;
  std::vector<int> avg_ext_thin_clad_1004;
  std::vector<int> avg_ext_thin_poly_1011;
  std::vector<int> avg_ext_thin_clad_1011;
  std::vector<int> avg_ext_thin_poly_1018;
  std::vector<int> avg_ext_thin_clad_1018;
  std::vector<int> avg_ext_thin_poly_1025;
  std::vector<int> avg_ext_thin_clad_1025;
  std::vector<int> avg_ext_thin_poly_1101;
  std::vector<int> avg_ext_thin_clad_1101;
  std::vector<int> avg_ext_thin_poly_1108;
  std::vector<int> avg_ext_thin_clad_1108;
  std::vector<int> avg_ext_thin_poly_1115;
  std::vector<int> avg_ext_thin_clad_1115;
  std::vector<int> avg_ext_thin_poly_1122;
  std::vector<int> avg_ext_thin_clad_1122;
  std::vector<int> avg_ext_thin_poly_1129;
  std::vector<int> avg_ext_thin_clad_1129;
  std::vector<int> avg_ext_thin_poly_1206;
  std::vector<int> avg_ext_thin_clad_1206;
  std::vector<int> avg_ext_thin_poly_020123;
  std::vector<int> avg_ext_thin_clad_020123;
  std::vector<int> avg_ext_thin_poly_020923;
  std::vector<int> avg_ext_thin_clad_020923;
  std::vector<int> avg_ext_thin_poly_021723;
  std::vector<int> avg_ext_thin_clad_021723;
  std::vector<int> avg_ext_thin_poly_030923;
  std::vector<int> avg_ext_thin_clad_030923;
  std::vector<int> avg_ext_thin_poly_031623;
  std::vector<int> avg_ext_thin_clad_031623;
  std::vector<int> avg_ext_thin_poly_032323;
  std::vector<int> avg_ext_thin_clad_032323;
  std::vector<int> avg_ext_thin_poly_040623;
  std::vector<int> avg_ext_thin_clad_040623;
  std::vector<int> avg_ext_thin_poly_041323;
  std::vector<int> avg_ext_thin_clad_041323;
  std::vector<int> avg_ext_thin_poly_042123;
  std::vector<int> avg_ext_thin_clad_042123;
  std::vector<int> avg_ext_thin_clad_052423;
  std::vector<int> avg_ext_thin_clad_060823;
  std::vector<int> avg_ext_thin_clad_062323;

  std::vector<int> stab_ext_thin_poly_920;
  std::vector<int> stab_ext_thin_clad_920;

  // here is where we actually do the grouping and push indices into the above vectors
  int index = 0;
  for(std::string& s : key_strs){
    if(s.find("ext_thin_poly_906") != std::string::npos){
      avg_ext_thin_poly_906.push_back(index);
    }
    if(s.find("ext_thin_clad_906") != std::string::npos){
      avg_ext_thin_clad_906.push_back(index);
    }

    if(s.find("ext_thin_poly_913") != std::string::npos){
      avg_ext_thin_poly_913.push_back(index);
    }
    if(s.find("ext_thin_clad_913") != std::string::npos){
      avg_ext_thin_clad_913.push_back(index);
    }

    if(s.find("ext_thin_poly_920") != std::string::npos){
      avg_ext_thin_poly_920.push_back(index);
    }
    if(s.find("ext_thin_clad_920") != std::string::npos){
      avg_ext_thin_clad_920.push_back(index);
    }

    if(s.find("ext_thin_poly_927") != std::string::npos){
      avg_ext_thin_poly_927.push_back(index);
    }
    if(s.find("ext_thin_clad_927") != std::string::npos){
      avg_ext_thin_clad_927.push_back(index);
    }

    if(s.find("ext_thin_poly_1004") != std::string::npos){
      avg_ext_thin_poly_1004.push_back(index);
    }
    if(s.find("ext_thin_clad_1004") != std::string::npos){
      avg_ext_thin_clad_1004.push_back(index);
    }

    if(s.find("ext_thin_poly_1011") != std::string::npos){
      avg_ext_thin_poly_1011.push_back(index);
    }
    if(s.find("ext_thin_clad_1011") != std::string::npos){
      avg_ext_thin_clad_1011.push_back(index);
    }

    if(s.find("ext_thin_poly_1018") != std::string::npos){
      avg_ext_thin_poly_1018.push_back(index);
    }
    if(s.find("ext_thin_clad_1018") != std::string::npos){
      avg_ext_thin_clad_1018.push_back(index);
    }

    if(s.find("ext_thin_poly_1025") != std::string::npos){
      avg_ext_thin_poly_1025.push_back(index);
    }
    if(s.find("ext_thin_clad_1025") != std::string::npos){
      avg_ext_thin_clad_1025.push_back(index);
    }

    if(s.find("ext_thin_poly_1101") != std::string::npos){
      avg_ext_thin_poly_1101.push_back(index);
    }
    if(s.find("ext_thin_clad_1101") != std::string::npos){
      avg_ext_thin_clad_1101.push_back(index);
    }

    if(s.find("ext_thin_poly_1108") != std::string::npos){
      avg_ext_thin_poly_1108.push_back(index);
    }
    if(s.find("ext_thin_clad_1108") != std::string::npos){
      avg_ext_thin_clad_1108.push_back(index);
    }

    if(s.find("ext_thin_poly_1115") != std::string::npos){
      avg_ext_thin_poly_1115.push_back(index);
    }
    if(s.find("ext_thin_clad_1115") != std::string::npos){
      avg_ext_thin_clad_1115.push_back(index);
    }

    if(s.find("ext_thin_poly_1122") != std::string::npos){
      avg_ext_thin_poly_1122.push_back(index);
    }
    if(s.find("ext_thin_clad_1122") != std::string::npos){
      avg_ext_thin_clad_1122.push_back(index);
    }

    if(s.find("ext_thin_poly_1129") != std::string::npos){
      avg_ext_thin_poly_1129.push_back(index);
    }
    if(s.find("ext_thin_clad_1129") != std::string::npos){
      avg_ext_thin_clad_1129.push_back(index);
    }

    if(s.find("ext_thin_poly_1206") != std::string::npos){
      avg_ext_thin_poly_1206.push_back(index);
    }
    if(s.find("ext_thin_clad_1206") != std::string::npos){
      avg_ext_thin_clad_1206.push_back(index);
    }

    if(s.find("ext_thin_poly_020123") != std::string::npos){
      avg_ext_thin_poly_020123.push_back(index);
    }
    if(s.find("ext_thin_clad_020123") != std::string::npos){
      avg_ext_thin_clad_020123.push_back(index);
    }

    if(s.find("ext_thin_poly_020923") != std::string::npos){
      avg_ext_thin_poly_020923.push_back(index);
    }
    if(s.find("ext_thin_clad_020923") != std::string::npos){
      avg_ext_thin_clad_020923.push_back(index);
    }

    if(s.find("ext_thin_poly_021723") != std::string::npos){
      avg_ext_thin_poly_021723.push_back(index);
    }
    if(s.find("ext_thin_clad_021723") != std::string::npos){
      avg_ext_thin_clad_021723.push_back(index);
    }

    if(s.find("ext_thin_poly_030923") != std::string::npos){
      avg_ext_thin_poly_030923.push_back(index);
    }
    if(s.find("ext_thin_clad_030923") != std::string::npos){
      avg_ext_thin_clad_030923.push_back(index);
    }

    if(s.find("ext_thin_poly_031623") != std::string::npos){
      avg_ext_thin_poly_031623.push_back(index);
    }
    if(s.find("ext_thin_clad_031623") != std::string::npos){
      avg_ext_thin_clad_031623.push_back(index);
    }

    if(s.find("ext_thin_poly_032323") != std::string::npos){
      avg_ext_thin_poly_032323.push_back(index);
    }
    if(s.find("ext_thin_clad_032323") != std::string::npos){
      avg_ext_thin_clad_032323.push_back(index);
    }

    if(s.find("ext_thin_poly_040623") != std::string::npos){
      avg_ext_thin_poly_040623.push_back(index);
    }
    if(s.find("ext_thin_clad_040623") != std::string::npos){
      avg_ext_thin_clad_040623.push_back(index);
    }

    if(s.find("ext_thin_poly_041323") != std::string::npos){
      avg_ext_thin_poly_041323.push_back(index);
    }
    if(s.find("ext_thin_clad_041323") != std::string::npos){
      avg_ext_thin_clad_041323.push_back(index);
    }

    if(s.find("ext_thin_poly_042123") != std::string::npos){
      avg_ext_thin_poly_042123.push_back(index);
    }
    if(s.find("ext_thin_clad_042123") != std::string::npos){
      avg_ext_thin_clad_042123.push_back(index);
    }

    if(s.find("ext_thin_clad_052423") != std::string::npos){
      avg_ext_thin_clad_052423.push_back(index);
    }
    if(s.find("ext_thin_clad_060823") != std::string::npos){
      avg_ext_thin_clad_060823.push_back(index);
    }
    if(s.find("ext_thin_clad_062323") != std::string::npos){
      avg_ext_thin_clad_062323.push_back(index);
    }


    if(s.find("ext_thin_poly_920_4") != std::string::npos){
      stab_ext_thin_poly_920.push_back(index);
    }
    if(s.find("ext_thin_clad_920_4") != std::string::npos){
      stab_ext_thin_clad_920.push_back(index);
    }

    index++;
  }

  // now, do fanicier analysis here using groups from above. Averages, differences, etc.

  // average histogram function
  ext_thin_poly_906_avg = average_hist(input_hists, avg_ext_thin_poly_906, "avg_ext_thin_poly_906");
  ext_thin_clad_906_avg = average_hist(input_hists, avg_ext_thin_clad_906, "avg_ext_thin_clad_906");

  ext_thin_poly_913_avg = average_hist(input_hists, avg_ext_thin_poly_913, "avg_ext_thin_poly_913");
  ext_thin_clad_913_avg = average_hist(input_hists, avg_ext_thin_clad_913, "avg_ext_thin_clad_913");

  ext_thin_poly_920_avg = average_hist(input_hists, avg_ext_thin_poly_920, "avg_ext_thin_poly_920");
  ext_thin_clad_920_avg = average_hist(input_hists, avg_ext_thin_clad_920, "avg_ext_thin_clad_920");

  ext_thin_poly_927_avg = average_hist(input_hists, avg_ext_thin_poly_927, "avg_ext_thin_poly_927");
  ext_thin_clad_927_avg = average_hist(input_hists, avg_ext_thin_clad_927, "avg_ext_thin_clad_927");

  ext_thin_poly_1004_avg = average_hist(input_hists, avg_ext_thin_poly_1004, "avg_ext_thin_poly_1004");
  ext_thin_clad_1004_avg = average_hist(input_hists, avg_ext_thin_clad_1004, "avg_ext_thin_clad_1004");

  ext_thin_poly_1011_avg = average_hist(input_hists, avg_ext_thin_poly_1011, "avg_ext_thin_poly_1011");
  ext_thin_clad_1011_avg = average_hist(input_hists, avg_ext_thin_clad_1011, "avg_ext_thin_clad_1011");

  ext_thin_poly_1018_avg = average_hist(input_hists, avg_ext_thin_poly_1018, "avg_ext_thin_poly_1018");
  ext_thin_clad_1018_avg = average_hist(input_hists, avg_ext_thin_clad_1018, "avg_ext_thin_clad_1018");

  ext_thin_poly_1025_avg = average_hist(input_hists, avg_ext_thin_poly_1025, "avg_ext_thin_poly_1025");
  ext_thin_clad_1025_avg = average_hist(input_hists, avg_ext_thin_clad_1025, "avg_ext_thin_clad_1025");

  ext_thin_poly_1101_avg = average_hist(input_hists, avg_ext_thin_poly_1101, "avg_ext_thin_poly_1101");
  ext_thin_clad_1101_avg = average_hist(input_hists, avg_ext_thin_clad_1101, "avg_ext_thin_clad_1101");

  ext_thin_poly_1108_avg = average_hist(input_hists, avg_ext_thin_poly_1108, "avg_ext_thin_poly_1108");
  ext_thin_clad_1108_avg = average_hist(input_hists, avg_ext_thin_clad_1108, "avg_ext_thin_clad_1108");

  ext_thin_poly_1115_avg = average_hist(input_hists, avg_ext_thin_poly_1115, "avg_ext_thin_poly_1115");
  ext_thin_clad_1115_avg = average_hist(input_hists, avg_ext_thin_clad_1115, "avg_ext_thin_clad_1115");

  ext_thin_poly_1122_avg = average_hist(input_hists, avg_ext_thin_poly_1122, "avg_ext_thin_poly_1122");
  ext_thin_clad_1122_avg = average_hist(input_hists, avg_ext_thin_clad_1122, "avg_ext_thin_clad_1122");

  ext_thin_poly_1129_avg = average_hist(input_hists, avg_ext_thin_poly_1129, "avg_ext_thin_poly_1129");
  ext_thin_clad_1129_avg = average_hist(input_hists, avg_ext_thin_clad_1129, "avg_ext_thin_clad_1129");

  ext_thin_poly_1206_avg = average_hist(input_hists, avg_ext_thin_poly_1206, "avg_ext_thin_poly_1206");
  ext_thin_clad_1206_avg = average_hist(input_hists, avg_ext_thin_clad_1206, "avg_ext_thin_clad_1206");

  ext_thin_poly_020123_avg = average_hist(input_hists, avg_ext_thin_poly_020123, "avg_ext_thin_poly_020123");
  ext_thin_clad_020123_avg = average_hist(input_hists, avg_ext_thin_clad_020123, "avg_ext_thin_clad_020123");

  ext_thin_poly_020923_avg = average_hist(input_hists, avg_ext_thin_poly_020923, "avg_ext_thin_poly_020923");
  ext_thin_clad_020923_avg = average_hist(input_hists, avg_ext_thin_clad_020923, "avg_ext_thin_clad_020923");

  ext_thin_poly_021723_avg = average_hist(input_hists, avg_ext_thin_poly_021723, "avg_ext_thin_poly_021723");
  ext_thin_clad_021723_avg = average_hist(input_hists, avg_ext_thin_clad_021723, "avg_ext_thin_clad_021723");

  ext_thin_poly_030923_avg = average_hist(input_hists, avg_ext_thin_poly_030923, "avg_ext_thin_poly_030923");
  ext_thin_clad_030923_avg = average_hist(input_hists, avg_ext_thin_clad_030923, "avg_ext_thin_clad_030923");

  ext_thin_poly_031623_avg = average_hist(input_hists, avg_ext_thin_poly_031623, "avg_ext_thin_poly_031623");
  ext_thin_clad_031623_avg = average_hist(input_hists, avg_ext_thin_clad_031623, "avg_ext_thin_clad_031623");

  ext_thin_poly_032323_avg = average_hist(input_hists, avg_ext_thin_poly_032323, "avg_ext_thin_poly_032323");
  ext_thin_clad_032323_avg = average_hist(input_hists, avg_ext_thin_clad_032323, "avg_ext_thin_clad_032323");

  ext_thin_poly_040623_avg = average_hist(input_hists, avg_ext_thin_poly_040623, "avg_ext_thin_poly_040623");
  ext_thin_clad_040623_avg = average_hist(input_hists, avg_ext_thin_clad_040623, "avg_ext_thin_clad_040623");

  ext_thin_poly_041323_avg = average_hist(input_hists, avg_ext_thin_poly_041323, "avg_ext_thin_poly_041323");
  ext_thin_clad_041323_avg = average_hist(input_hists, avg_ext_thin_clad_041323, "avg_ext_thin_clad_041323");

  ext_thin_poly_042123_avg = average_hist(input_hists, avg_ext_thin_poly_042123, "avg_ext_thin_poly_042123");
  ext_thin_clad_042123_avg = average_hist(input_hists, avg_ext_thin_clad_042123, "avg_ext_thin_clad_042123");

  ext_thin_clad_052423_avg = average_hist(input_hists, avg_ext_thin_clad_052423, "avg_ext_thin_clad_052423");
  ext_thin_clad_060823_avg = average_hist(input_hists, avg_ext_thin_clad_060823, "avg_ext_thin_clad_060823");
  ext_thin_clad_062323_avg = average_hist(input_hists, avg_ext_thin_clad_062323, "avg_ext_thin_clad_062323");

  // Make sure ROOT knows we want to write to the output file
  ofile->cd();

  // write the histograms to file
  ext_thin_poly_906_avg->Write();
  ext_thin_clad_906_avg->Write();
  ext_thin_poly_913_avg->Write();
  ext_thin_clad_913_avg->Write();
  ext_thin_poly_920_avg->Write();
  ext_thin_clad_920_avg->Write();
  ext_thin_poly_927_avg->Write();
  ext_thin_clad_927_avg->Write();
  ext_thin_poly_1004_avg->Write();
  ext_thin_clad_1004_avg->Write();
  ext_thin_poly_1011_avg->Write();
  ext_thin_clad_1011_avg->Write();
  ext_thin_poly_1018_avg->Write();
  ext_thin_clad_1018_avg->Write();
  ext_thin_poly_1025_avg->Write();
  ext_thin_clad_1025_avg->Write();
  ext_thin_poly_1101_avg->Write();
  ext_thin_clad_1101_avg->Write();
  ext_thin_poly_1108_avg->Write();
  ext_thin_clad_1108_avg->Write();
  ext_thin_poly_1115_avg->Write();
  ext_thin_clad_1115_avg->Write();
  ext_thin_poly_1122_avg->Write();
  ext_thin_clad_1122_avg->Write();
  ext_thin_poly_1129_avg->Write();
  ext_thin_clad_1129_avg->Write();
  ext_thin_poly_1206_avg->Write();
  ext_thin_clad_1206_avg->Write();
  ext_thin_poly_020123_avg->Write();
  ext_thin_clad_020123_avg->Write();
  ext_thin_poly_020923_avg->Write();
  ext_thin_clad_020923_avg->Write();  
  ext_thin_poly_021723_avg->Write();
  ext_thin_clad_021723_avg->Write();
  ext_thin_poly_030923_avg->Write();
  ext_thin_clad_030923_avg->Write();
  ext_thin_poly_031623_avg->Write();
  ext_thin_clad_031623_avg->Write();
  ext_thin_poly_032323_avg->Write();
  ext_thin_clad_032323_avg->Write();
  ext_thin_poly_040623_avg->Write();
  ext_thin_clad_040623_avg->Write();
  ext_thin_poly_041323_avg->Write();
  ext_thin_clad_041323_avg->Write();
  ext_thin_poly_042123_avg->Write();
  ext_thin_clad_042123_avg->Write();
  ext_thin_clad_052423_avg->Write();
  ext_thin_clad_060823_avg->Write();
  ext_thin_clad_062323_avg->Write();

  // make collections of average data for wavelengths of interest, for aging
  std::vector<float> aging_380_ext_thin_poly;
  std::vector<float> aging_380_ext_thin_clad;
  std::vector<float> aging_380_ext_thin_poly_err;
  std::vector<float> aging_380_ext_thin_clad_err;
  std::vector<float> aging_390_ext_thin_poly;
  std::vector<float> aging_390_ext_thin_clad;
  std::vector<float> aging_390_ext_thin_poly_err;
  std::vector<float> aging_390_ext_thin_clad_err;
  std::vector<float> aging_400_ext_thin_poly;
  std::vector<float> aging_400_ext_thin_clad;
  std::vector<float> aging_400_ext_thin_poly_err;
  std::vector<float> aging_400_ext_thin_clad_err;
  std::vector<float> aging_410_ext_thin_poly;
  std::vector<float> aging_410_ext_thin_clad;
  std::vector<float> aging_410_ext_thin_poly_err;
  std::vector<float> aging_410_ext_thin_clad_err;
  std::vector<float> aging_420_ext_thin_poly;
  std::vector<float> aging_420_ext_thin_clad;
  std::vector<float> aging_420_ext_thin_poly_err;
  std::vector<float> aging_420_ext_thin_clad_err;
  std::vector<float> aging_430_ext_thin_poly;
  std::vector<float> aging_430_ext_thin_clad;
  std::vector<float> aging_430_ext_thin_poly_err;
  std::vector<float> aging_430_ext_thin_clad_err;
  std::vector<float> aging_440_ext_thin_poly;
  std::vector<float> aging_440_ext_thin_clad;
  std::vector<float> aging_440_ext_thin_poly_err;
  std::vector<float> aging_440_ext_thin_clad_err;
  std::vector<float> aging_450_ext_thin_poly;
  std::vector<float> aging_450_ext_thin_clad;
  std::vector<float> aging_450_ext_thin_poly_err;
  std::vector<float> aging_450_ext_thin_clad_err;
  std::vector<float> stab_ext_thin_poly;
  std::vector<float> stab_ext_thin_clad;

  // make the x-axis points array, time fraction of years
  const int npoints = 25;
  float aging_timefrac_xbins[npoints] = {0.0, 0.019, 0.038, 0.058, 0.077, 0.096, 0.115, 0.134, 0.153, 0.173, 0.192, 0.211, 0.230, 0.249, 0.405, 0.427, 0.449, 0.504, 0.523, 0.542, 0.581, 0.600, 0.622, 0.753, 0.795};

  float aging_empty_xbins_err[npoints];
  std::fill(aging_empty_xbins_err, aging_empty_xbins_err+npoints, 0);

  // wavelength = 380 bin content
  aging_380_ext_thin_poly.push_back(ext_thin_poly_906_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_906_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_913_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_913_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_920_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_920_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_927_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_927_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_1004_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_1004_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_1011_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_1011_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_1018_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_1018_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_1025_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_1025_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_1101_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_1101_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_1108_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_1108_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_1115_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_1115_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_1122_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_1122_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_1129_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_1129_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_1206_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_1206_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_020123_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_020123_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_020923_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_020923_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_021723_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_021723_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_030923_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_030923_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_031623_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_031623_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_032323_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_032323_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_040623_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_040623_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_041323_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_041323_avg->GetBinContent(3));
  aging_380_ext_thin_poly.push_back(ext_thin_poly_042123_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_042123_avg->GetBinContent(3));
  //aging_380_ext_thin_clad.push_back(ext_thin_clad_052423_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_060823_avg->GetBinContent(3));
  aging_380_ext_thin_clad.push_back(ext_thin_clad_062323_avg->GetBinContent(3));
  // wavelength = 390 bin content
  aging_390_ext_thin_poly.push_back(ext_thin_poly_906_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_906_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_913_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_913_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_920_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_920_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_927_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_927_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_1004_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_1004_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_1011_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_1011_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_1018_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_1018_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_1025_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_1025_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_1101_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_1101_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_1108_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_1108_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_1115_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_1115_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_1122_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_1122_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_1129_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_1129_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_1206_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_1206_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_020123_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_020123_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_020923_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_020923_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_021723_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_021723_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_030923_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_030923_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_031623_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_031623_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_032323_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_032323_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_040623_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_040623_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_041323_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_041323_avg->GetBinContent(4));
  aging_390_ext_thin_poly.push_back(ext_thin_poly_042123_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_042123_avg->GetBinContent(4));
  //aging_390_ext_thin_clad.push_back(ext_thin_clad_052423_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_060823_avg->GetBinContent(4));
  aging_390_ext_thin_clad.push_back(ext_thin_clad_062323_avg->GetBinContent(4));
  // wavelength = 400 bin content
  aging_400_ext_thin_poly.push_back(ext_thin_poly_906_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_906_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_913_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_913_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_920_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_920_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_927_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_927_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_1004_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_1004_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_1011_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_1011_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_1018_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_1018_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_1025_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_1025_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_1101_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_1101_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_1108_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_1108_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_1115_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_1115_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_1122_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_1122_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_1129_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_1129_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_1206_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_1206_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_020123_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_020123_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_020923_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_020923_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_021723_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_021723_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_030923_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_030923_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_031623_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_031623_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_032323_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_032323_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_040623_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_040623_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_041323_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_041323_avg->GetBinContent(5));
  aging_400_ext_thin_poly.push_back(ext_thin_poly_042123_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_042123_avg->GetBinContent(5));
  //aging_400_ext_thin_clad.push_back(ext_thin_clad_052423_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_060823_avg->GetBinContent(5));
  aging_400_ext_thin_clad.push_back(ext_thin_clad_062323_avg->GetBinContent(5));
  // wavelength = 410 bin content
  aging_410_ext_thin_poly.push_back(ext_thin_poly_906_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_906_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_913_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_913_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_920_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_920_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_927_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_927_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_1004_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_1004_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_1011_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_1011_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_1018_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_1018_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_1025_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_1025_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_1101_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_1101_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_1108_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_1108_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_1115_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_1115_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_1122_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_1122_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_1129_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_1129_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_1206_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_1206_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_020123_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_020123_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_020923_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_020923_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_021723_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_021723_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_030923_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_030923_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_031623_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_031623_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_032323_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_032323_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_040623_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_040623_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_041323_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_041323_avg->GetBinContent(6));
  aging_410_ext_thin_poly.push_back(ext_thin_poly_042123_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_042123_avg->GetBinContent(6));
  //aging_410_ext_thin_clad.push_back(ext_thin_clad_052423_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_060823_avg->GetBinContent(6));
  aging_410_ext_thin_clad.push_back(ext_thin_clad_062323_avg->GetBinContent(6));
  // wavelength = 420 bin content
  aging_420_ext_thin_poly.push_back(ext_thin_poly_906_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_906_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_913_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_913_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_920_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_920_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_927_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_927_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_1004_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_1004_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_1011_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_1011_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_1018_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_1018_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_1025_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_1025_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_1101_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_1101_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_1108_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_1108_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_1115_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_1115_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_1122_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_1122_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_1129_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_1129_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_1206_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_1206_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_020123_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_020123_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_020923_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_020923_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_021723_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_021723_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_030923_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_030923_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_031623_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_031623_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_032323_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_032323_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_040623_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_040623_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_041323_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_041323_avg->GetBinContent(7));
  aging_420_ext_thin_poly.push_back(ext_thin_poly_042123_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_042123_avg->GetBinContent(7));
  //aging_420_ext_thin_clad.push_back(ext_thin_clad_052423_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_060823_avg->GetBinContent(7));
  aging_420_ext_thin_clad.push_back(ext_thin_clad_062323_avg->GetBinContent(7));
  // wavelength = 430 bin content
  aging_430_ext_thin_poly.push_back(ext_thin_poly_906_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_906_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_913_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_913_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_920_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_920_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_927_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_927_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_1004_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_1004_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_1011_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_1011_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_1018_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_1018_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_1025_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_1025_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_1101_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_1101_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_1108_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_1108_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_1115_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_1115_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_1122_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_1122_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_1129_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_1129_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_1206_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_1206_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_020123_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_020123_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_020923_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_020923_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_021723_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_021723_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_030923_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_030923_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_031623_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_031623_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_032323_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_032323_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_040623_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_040623_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_041323_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_041323_avg->GetBinContent(8));
  aging_430_ext_thin_poly.push_back(ext_thin_poly_042123_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_042123_avg->GetBinContent(8));
  //aging_430_ext_thin_clad.push_back(ext_thin_clad_052423_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_060823_avg->GetBinContent(8));
  aging_430_ext_thin_clad.push_back(ext_thin_clad_062323_avg->GetBinContent(8));
  // wavelength = 440 bin content
  aging_440_ext_thin_poly.push_back(ext_thin_poly_906_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_906_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_913_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_913_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_920_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_920_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_927_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_927_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_1004_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_1004_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_1011_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_1011_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_1018_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_1018_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_1025_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_1025_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_1101_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_1101_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_1108_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_1108_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_1115_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_1115_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_1122_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_1122_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_1129_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_1129_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_1206_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_1206_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_020123_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_020123_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_020923_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_020923_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_021723_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_021723_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_030923_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_030923_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_031623_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_031623_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_032323_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_032323_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_040623_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_040623_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_041323_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_041323_avg->GetBinContent(9));
  aging_440_ext_thin_poly.push_back(ext_thin_poly_042123_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_042123_avg->GetBinContent(9));
  //aging_440_ext_thin_clad.push_back(ext_thin_clad_052423_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_060823_avg->GetBinContent(9));
  aging_440_ext_thin_clad.push_back(ext_thin_clad_062323_avg->GetBinContent(9));
  // wavelength = 450 bin content
  aging_450_ext_thin_poly.push_back(ext_thin_poly_906_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_906_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_913_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_913_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_920_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_920_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_927_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_927_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_1004_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_1004_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_1011_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_1011_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_1018_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_1018_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_1025_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_1025_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_1101_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_1101_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_1108_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_1108_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_1115_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_1115_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_1122_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_1122_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_1129_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_1129_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_1206_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_1206_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_020123_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_020123_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_020923_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_020923_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_021723_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_021723_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_030923_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_030923_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_031623_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_031623_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_032323_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_032323_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_040623_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_040623_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_041323_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_041323_avg->GetBinContent(10));
  aging_450_ext_thin_poly.push_back(ext_thin_poly_042123_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_042123_avg->GetBinContent(10));
  //aging_450_ext_thin_clad.push_back(ext_thin_clad_052423_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_060823_avg->GetBinContent(10));
  aging_450_ext_thin_clad.push_back(ext_thin_clad_062323_avg->GetBinContent(10));

  // do errors here, standard deviations

  // wavelength = 380nm
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_906, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_906, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_913, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_913, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_920, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_920, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_927, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_927, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1004, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1004, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1011, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1011, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1018, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1018, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1025, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1025, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1101, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1101, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1108, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1108, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1115, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1115, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1122, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1122, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1129, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1129, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1206, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1206, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020123, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020123, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020923, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020923, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_021723, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_021723, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_030923, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_030923, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_031623, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_031623, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_032323, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_032323, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_040623, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_040623, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_041323, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_041323, 380));
  aging_380_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_042123, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_042123, 380));
  //aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_052423, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_060823, 380));
  aging_380_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_062323, 380));
  // wavelength = 390nm
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_906, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_906, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_913, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_913, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_920, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_920, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_927, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_927, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1004, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1004, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1011, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1011, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1018, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1018, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1025, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1025, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1101, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1101, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1108, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1108, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1115, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1115, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1122, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1122, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1129, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1129, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1206, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1206, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020123, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020123, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020923, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020923, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_021723, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_021723, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_030923, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_030923, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_031623, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_031623, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_032323, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_032323, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_040623, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_040623, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_041323, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_041323, 390));
  aging_390_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_042123, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_042123, 390));
  //aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_052423, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_060823, 390));
  aging_390_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_062323, 390));
  // wavelength = 400nm
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_906, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_906, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_913, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_913, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_920, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_920, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_927, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_927, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1004, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1004, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1011, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1011, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1018, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1018, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1025, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1025, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1101, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1101, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1108, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1108, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1115, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1115, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1122, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1122, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1129, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1129, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1206, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1206, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020123, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020123, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020923, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020923, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_021723, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_021723, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_030923, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_030923, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_031623, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_031623, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_032323, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_032323, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_040623, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_040623, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_041323, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_041323, 400));
  aging_400_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_042123, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_042123, 400));
  //aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_052423, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_060823, 400));
  aging_400_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_062323, 400));
  // wavelength = 410nm
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_906, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_906, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_913, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_913, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_920, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_920, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_927, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_927, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1004, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1004, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1011, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1011, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1018, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1018, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1025, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1025, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1101, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1101, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1108, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1108, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1115, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1115, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1122, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1122, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1129, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1129, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1206, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1206, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020123, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020123, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020923, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020923, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_021723, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_021723, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_030923, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_030923, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_031623, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_031623, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_032323, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_032323, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_040623, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_040623, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_041323, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_041323, 410));
  aging_410_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_042123, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_042123, 410));
  //aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_052423, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_060823, 410));
  aging_410_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_062323, 410));
  // wavelength = 420nm
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_906, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_906, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_913, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_913, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_920, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_920, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_927, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_927, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1004, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1004, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1011, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1011, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1018, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1018, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1025, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1025, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1101, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1101, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1108, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1108, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1115, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1115, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1122, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1122, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1129, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1129, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1206, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1206, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020123, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020123, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020923, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020923, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_021723, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_021723, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_030923, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_030923, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_031623, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_031623, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_032323, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_032323, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_040623, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_040623, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_041323, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_041323, 420));
  aging_420_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_042123, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_042123, 420));
  //aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_052423, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_060823, 420));
  aging_420_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_062323, 420));
  // wavelength = 430nm
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_906, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_906, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_913, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_913, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_920, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_920, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_927, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_927, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1004, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1004, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1011, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1011, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1018, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1018, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1025, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1025, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1101, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1101, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1108, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1108, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1115, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1115, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1122, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1122, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1129, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1129, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1206, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1206, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020123, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020123, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020923, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020923, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_021723, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_021723, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_030923, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_030923, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_031623, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_031623, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_032323, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_032323, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_040623, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_040623, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_041323, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_041323, 430));
  aging_430_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_042123, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_042123, 430));
  //aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_052423, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_060823, 430));
  aging_430_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_062323, 430));
  // wavelength = 440nm
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_906, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_906, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_913, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_913, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_920, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_920, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_927, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_927, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1004, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1004, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1011, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1011, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1018, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1018, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1025, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1025, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1101, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1101, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1108, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1108, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1115, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1115, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1122, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1122, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1129, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1129, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1206, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1206, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020123, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020123, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020923, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020923, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_021723, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_021723, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_030923, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_030923, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_031623, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_031623, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_032323, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_032323, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_040623, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_040623, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_041323, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_041323, 440));
  aging_440_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_042123, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_042123, 440));
  //aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_052423, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_060823, 440));
  aging_440_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_062323, 440));
  // wavelength = 450nm
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_906, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_906, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_913, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_913, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_920, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_920, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_927, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_927, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1004, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1004, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1011, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1011, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1018, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1018, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1025, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1025, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1101, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1101, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1108, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1108, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1115, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1115, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1122, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1122, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1129, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1129, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_1206, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_1206, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020123, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020123, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_020923, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_020923, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_021723, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_021723, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_030923, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_030923, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_031623, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_031623, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_032323, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_032323, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_040623, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_040623, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_041323, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_041323, 450));
  aging_450_ext_thin_poly_err.push_back(stdev_error(input_hists, avg_ext_thin_poly_042123, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_042123, 450));
  //aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_052423, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_060823, 450));
  aging_450_ext_thin_clad_err.push_back(stdev_error(input_hists, avg_ext_thin_clad_062323, 450));

  // stability standard deviations
  stab_ext_thin_poly.push_back(stdev_error(input_hists, stab_ext_thin_poly_920, 380));
  stab_ext_thin_poly.push_back(stdev_error(input_hists, stab_ext_thin_poly_920, 390));
  stab_ext_thin_poly.push_back(stdev_error(input_hists, stab_ext_thin_poly_920, 400));
  stab_ext_thin_poly.push_back(stdev_error(input_hists, stab_ext_thin_poly_920, 410));
  stab_ext_thin_poly.push_back(stdev_error(input_hists, stab_ext_thin_poly_920, 420));
  stab_ext_thin_poly.push_back(stdev_error(input_hists, stab_ext_thin_poly_920, 430));
  stab_ext_thin_poly.push_back(stdev_error(input_hists, stab_ext_thin_poly_920, 440));
  stab_ext_thin_poly.push_back(stdev_error(input_hists, stab_ext_thin_poly_920, 450));
  stab_ext_thin_clad.push_back(stdev_error(input_hists, stab_ext_thin_clad_920, 380));
  stab_ext_thin_clad.push_back(stdev_error(input_hists, stab_ext_thin_clad_920, 390));
  stab_ext_thin_clad.push_back(stdev_error(input_hists, stab_ext_thin_clad_920, 400));
  stab_ext_thin_clad.push_back(stdev_error(input_hists, stab_ext_thin_clad_920, 410));
  stab_ext_thin_clad.push_back(stdev_error(input_hists, stab_ext_thin_clad_920, 420));
  stab_ext_thin_clad.push_back(stdev_error(input_hists, stab_ext_thin_clad_920, 430));
  stab_ext_thin_clad.push_back(stdev_error(input_hists, stab_ext_thin_clad_920, 440));
  stab_ext_thin_clad.push_back(stdev_error(input_hists, stab_ext_thin_clad_920, 450));

  // make TCanvases and TGraphErrors objects to hold aging plots
  gStyle->SetOptFit(11111);

  TCanvas* page[16];
  TCanvas* stability[2];
  TGraphErrors* tg_ext_thin_poly_frac_380;
  TGraphErrors* tg_ext_thin_clad_frac_380;
  TGraphErrors* tg_ext_thin_poly_frac_390;
  TGraphErrors* tg_ext_thin_clad_frac_390;
  TGraphErrors* tg_ext_thin_poly_frac_400;
  TGraphErrors* tg_ext_thin_clad_frac_400;
  TGraphErrors* tg_ext_thin_poly_frac_410;
  TGraphErrors* tg_ext_thin_clad_frac_410;
  TGraphErrors* tg_ext_thin_poly_frac_420;
  TGraphErrors* tg_ext_thin_clad_frac_420;
  TGraphErrors* tg_ext_thin_poly_frac_430;
  TGraphErrors* tg_ext_thin_clad_frac_430;
  TGraphErrors* tg_ext_thin_poly_frac_440;
  TGraphErrors* tg_ext_thin_clad_frac_440;
  TGraphErrors* tg_ext_thin_poly_frac_450;
  TGraphErrors* tg_ext_thin_clad_frac_450;
  TGraph* tg_poly_stab;
  TGraph* tg_clad_stab;

  page[12] = new TCanvas("poly_frac_380","poly_frac_380",900,600);
  page[13] = new TCanvas("clad_frac_380","clad_frac_380",900,600);
  page[14] = new TCanvas("poly_frac_390","poly_frac_390",900,600);
  page[15] = new TCanvas("clad_frac_390","clad_frac_390",900,600);
  page[0] = new TCanvas("poly_frac_400","poly_frac_400",900,600);
  page[1] = new TCanvas("clad_frac_400","clad_frac_400",900,600);
  page[2] = new TCanvas("poly_frac_410","poly_frac_410",900,600);
  page[3] = new TCanvas("clad_frac_410","clad_frac_410",900,600);
  page[4] = new TCanvas("poly_frac_420","poly_frac_420",900,600);
  page[5] = new TCanvas("clad_frac_420","clad_frac_420",900,600);
  page[6] = new TCanvas("poly_frac_430","poly_frac_430",900,600);
  page[7] = new TCanvas("clad_frac_430","clad_frac_430",900,600);
  page[8] = new TCanvas("poly_frac_440","poly_frac_440",900,600);
  page[9] = new TCanvas("clad_frac_440","clad_frac_440",900,600);
  page[10] = new TCanvas("poly_frac_450","poly_frac_450",900,600);
  page[11] = new TCanvas("clad_frac_450","clad_frac_450",900,600);

  stability[0] = new TCanvas("poly_stability","poly_stability",900,600);
  stability[1] = new TCanvas("clad_stability","clad_stability",900,600);

  // need arrays for TGraphErrors

  // wavelength = 380nm
  float aging_380_ext_thin_poly_a[npoints];
  float aging_380_ext_thin_clad_a[npoints];
  float aging_380_ext_thin_poly_err_a[npoints];
  float aging_380_ext_thin_clad_err_a[npoints];
  std::copy(aging_380_ext_thin_poly.begin(), aging_380_ext_thin_poly.end(), aging_380_ext_thin_poly_a);
  std::copy(aging_380_ext_thin_clad.begin(), aging_380_ext_thin_clad.end(), aging_380_ext_thin_clad_a);
  std::copy(aging_380_ext_thin_poly_err.begin(), aging_380_ext_thin_poly_err.end(), aging_380_ext_thin_poly_err_a);
  std::copy(aging_380_ext_thin_clad_err.begin(), aging_380_ext_thin_clad_err.end(), aging_380_ext_thin_clad_err_a);
  // wavelength = 390nm
  float aging_390_ext_thin_poly_a[npoints];
  float aging_390_ext_thin_clad_a[npoints];
  float aging_390_ext_thin_poly_err_a[npoints];
  float aging_390_ext_thin_clad_err_a[npoints];
  std::copy(aging_390_ext_thin_poly.begin(), aging_390_ext_thin_poly.end(), aging_390_ext_thin_poly_a);
  std::copy(aging_390_ext_thin_clad.begin(), aging_390_ext_thin_clad.end(), aging_390_ext_thin_clad_a);
  std::copy(aging_390_ext_thin_poly_err.begin(), aging_390_ext_thin_poly_err.end(), aging_390_ext_thin_poly_err_a);
  std::copy(aging_390_ext_thin_clad_err.begin(), aging_390_ext_thin_clad_err.end(), aging_390_ext_thin_clad_err_a);
  // wavelength = 400nm
  float aging_400_ext_thin_poly_a[npoints];
  float aging_400_ext_thin_clad_a[npoints];
  float aging_400_ext_thin_poly_err_a[npoints];
  float aging_400_ext_thin_clad_err_a[npoints];
  std::copy(aging_400_ext_thin_poly.begin(), aging_400_ext_thin_poly.end(), aging_400_ext_thin_poly_a);
  std::copy(aging_400_ext_thin_clad.begin(), aging_400_ext_thin_clad.end(), aging_400_ext_thin_clad_a);
  std::copy(aging_400_ext_thin_poly_err.begin(), aging_400_ext_thin_poly_err.end(), aging_400_ext_thin_poly_err_a);
  std::copy(aging_400_ext_thin_clad_err.begin(), aging_400_ext_thin_clad_err.end(), aging_400_ext_thin_clad_err_a);
  // wavelength = 410nm
  float aging_410_ext_thin_poly_a[npoints];
  float aging_410_ext_thin_clad_a[npoints];
  float aging_410_ext_thin_poly_err_a[npoints];
  float aging_410_ext_thin_clad_err_a[npoints];
  std::copy(aging_410_ext_thin_poly.begin(), aging_410_ext_thin_poly.end(), aging_410_ext_thin_poly_a);
  std::copy(aging_410_ext_thin_clad.begin(), aging_410_ext_thin_clad.end(), aging_410_ext_thin_clad_a);
  std::copy(aging_410_ext_thin_poly_err.begin(), aging_410_ext_thin_poly_err.end(), aging_410_ext_thin_poly_err_a);
  std::copy(aging_410_ext_thin_clad_err.begin(), aging_410_ext_thin_clad_err.end(), aging_410_ext_thin_clad_err_a);
  // wavelength = 420nm
  float aging_420_ext_thin_poly_a[npoints];
  float aging_420_ext_thin_clad_a[npoints];
  float aging_420_ext_thin_poly_err_a[npoints];
  float aging_420_ext_thin_clad_err_a[npoints];
  std::copy(aging_420_ext_thin_poly.begin(), aging_420_ext_thin_poly.end(), aging_420_ext_thin_poly_a);
  std::copy(aging_420_ext_thin_clad.begin(), aging_420_ext_thin_clad.end(), aging_420_ext_thin_clad_a);
  std::copy(aging_420_ext_thin_poly_err.begin(), aging_420_ext_thin_poly_err.end(), aging_420_ext_thin_poly_err_a);
  std::copy(aging_420_ext_thin_clad_err.begin(), aging_420_ext_thin_clad_err.end(), aging_420_ext_thin_clad_err_a);
  // wavelength = 430nm
  float aging_430_ext_thin_poly_a[npoints];
  float aging_430_ext_thin_clad_a[npoints];
  float aging_430_ext_thin_poly_err_a[npoints];
  float aging_430_ext_thin_clad_err_a[npoints];
  std::copy(aging_430_ext_thin_poly.begin(), aging_430_ext_thin_poly.end(), aging_430_ext_thin_poly_a);
  std::copy(aging_430_ext_thin_clad.begin(), aging_430_ext_thin_clad.end(), aging_430_ext_thin_clad_a);
  std::copy(aging_430_ext_thin_poly_err.begin(), aging_430_ext_thin_poly_err.end(), aging_430_ext_thin_poly_err_a);
  std::copy(aging_430_ext_thin_clad_err.begin(), aging_430_ext_thin_clad_err.end(), aging_430_ext_thin_clad_err_a);
  // wavelength = 440nm
  float aging_440_ext_thin_poly_a[npoints];
  float aging_440_ext_thin_clad_a[npoints];
  float aging_440_ext_thin_poly_err_a[npoints];
  float aging_440_ext_thin_clad_err_a[npoints];
  std::copy(aging_440_ext_thin_poly.begin(), aging_440_ext_thin_poly.end(), aging_440_ext_thin_poly_a);
  std::copy(aging_440_ext_thin_clad.begin(), aging_440_ext_thin_clad.end(), aging_440_ext_thin_clad_a);
  std::copy(aging_440_ext_thin_poly_err.begin(), aging_440_ext_thin_poly_err.end(), aging_440_ext_thin_poly_err_a);
  std::copy(aging_440_ext_thin_clad_err.begin(), aging_440_ext_thin_clad_err.end(), aging_440_ext_thin_clad_err_a);
  // wavelength = 450nm
  float aging_450_ext_thin_poly_a[npoints];
  float aging_450_ext_thin_clad_a[npoints];
  float aging_450_ext_thin_poly_err_a[npoints];
  float aging_450_ext_thin_clad_err_a[npoints];
  std::copy(aging_450_ext_thin_poly.begin(), aging_450_ext_thin_poly.end(), aging_450_ext_thin_poly_a);
  std::copy(aging_450_ext_thin_clad.begin(), aging_450_ext_thin_clad.end(), aging_450_ext_thin_clad_a);
  std::copy(aging_450_ext_thin_poly_err.begin(), aging_450_ext_thin_poly_err.end(), aging_450_ext_thin_poly_err_a);
  std::copy(aging_450_ext_thin_clad_err.begin(), aging_450_ext_thin_clad_err.end(), aging_450_ext_thin_clad_err_a);

  float stab_ext_thin_poly_a[8];
  float stab_ext_thin_clad_a[8];
  std::copy(stab_ext_thin_poly.begin(), stab_ext_thin_poly.end(), stab_ext_thin_poly_a);
  std::copy(stab_ext_thin_clad.begin(), stab_ext_thin_clad.end(), stab_ext_thin_clad_a);

  std::vector<float> poly_nm_slopes;
  std::vector<float> clad_nm_slopes;
  std::vector<float> poly_nm_slopes_err;
  std::vector<float> clad_nm_slopes_err;

  page[12]->cd();
  tg_ext_thin_poly_frac_380 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_380_ext_thin_poly_a, aging_empty_xbins_err, aging_380_ext_thin_poly_err_a);
  tg_ext_thin_poly_frac_380->SetTitle("Average Poly-side R at 380nm of New Ext vs Time");
  tg_ext_thin_poly_frac_380->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_poly_frac_380->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_poly_frac_380->SetMarkerStyle(20);
  tg_ext_thin_poly_frac_380->Fit("pol1");
  tg_ext_thin_poly_frac_380->Draw("AP");
  TF1 *fit12 = (TF1*)tg_ext_thin_poly_frac_380->GetListOfFunctions()->FindObject("pol1");
  poly_nm_slopes.push_back(fit12->GetParameter(1));
  poly_nm_slopes_err.push_back(fit12->GetParError(1));
  page[12]->Update();
  ofile->cd();
  page[12]->Write();
  page[12]->Close();

  page[13]->cd();
  tg_ext_thin_clad_frac_380 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_380_ext_thin_clad_a, aging_empty_xbins_err, aging_380_ext_thin_clad_err_a);
  tg_ext_thin_clad_frac_380->SetTitle("Average Clad-side R at 380nm of New Ext vs Time");
  tg_ext_thin_clad_frac_380->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_clad_frac_380->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_clad_frac_380->SetMarkerStyle(20);
  tg_ext_thin_clad_frac_380->Fit("pol1");
  tg_ext_thin_clad_frac_380->Draw("AP");
  TF1 *fit13 = (TF1*)tg_ext_thin_clad_frac_380->GetListOfFunctions()->FindObject("pol1");
  clad_nm_slopes.push_back(fit13->GetParameter(1));
  clad_nm_slopes_err.push_back(fit13->GetParError(1));
  page[13]->Update();
  ofile->cd();
  page[13]->Write();
  page[13]->Close();

  page[14]->cd();
  tg_ext_thin_poly_frac_390 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_390_ext_thin_poly_a, aging_empty_xbins_err, aging_390_ext_thin_poly_err_a);
  tg_ext_thin_poly_frac_390->SetTitle("Average Poly-side R at 390nm of New Ext vs Time");
  tg_ext_thin_poly_frac_390->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_poly_frac_390->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_poly_frac_390->SetMarkerStyle(20);
  tg_ext_thin_poly_frac_390->Fit("pol1");
  tg_ext_thin_poly_frac_390->Draw("AP");
  TF1 *fit14 = (TF1*)tg_ext_thin_poly_frac_390->GetListOfFunctions()->FindObject("pol1");
  poly_nm_slopes.push_back(fit14->GetParameter(1));
  poly_nm_slopes_err.push_back(fit14->GetParError(1));
  page[14]->Update();
  ofile->cd();
  page[14]->Write();
  page[14]->Close();

  page[15]->cd();
  tg_ext_thin_clad_frac_390 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_390_ext_thin_clad_a, aging_empty_xbins_err, aging_390_ext_thin_clad_err_a);
  tg_ext_thin_clad_frac_390->SetTitle("Average Clad-side R at 390nm of New Ext vs Time");
  tg_ext_thin_clad_frac_390->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_clad_frac_390->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_clad_frac_390->SetMarkerStyle(20);
  tg_ext_thin_clad_frac_390->Fit("pol1");
  tg_ext_thin_clad_frac_390->Draw("AP");
  TF1 *fit15 = (TF1*)tg_ext_thin_clad_frac_390->GetListOfFunctions()->FindObject("pol1");
  clad_nm_slopes.push_back(fit15->GetParameter(1));
  clad_nm_slopes_err.push_back(fit15->GetParError(1));
  page[15]->Update();
  ofile->cd();
  page[15]->Write();
  page[15]->Close();

  page[0]->cd();
  tg_ext_thin_poly_frac_400 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_400_ext_thin_poly_a, aging_empty_xbins_err, aging_400_ext_thin_poly_err_a);
  tg_ext_thin_poly_frac_400->SetTitle("Average Poly-side R at 400nm of New Ext vs Time");
  tg_ext_thin_poly_frac_400->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_poly_frac_400->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_poly_frac_400->SetMarkerStyle(20);
  tg_ext_thin_poly_frac_400->Fit("pol1");
  tg_ext_thin_poly_frac_400->Draw("AP");
  TF1 *fit0 = (TF1*)tg_ext_thin_poly_frac_400->GetListOfFunctions()->FindObject("pol1");
  poly_nm_slopes.push_back(fit0->GetParameter(1));
  poly_nm_slopes_err.push_back(fit0->GetParError(1));
  page[0]->Update();
  ofile->cd();
  page[0]->Write();
  page[0]->Close();

  page[1]->cd();
  tg_ext_thin_clad_frac_400 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_400_ext_thin_clad_a, aging_empty_xbins_err, aging_400_ext_thin_clad_err_a);
  tg_ext_thin_clad_frac_400->SetTitle("Average Clad-side R at 400nm of New Ext vs Time");
  tg_ext_thin_clad_frac_400->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_clad_frac_400->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_clad_frac_400->SetMarkerStyle(20);
  tg_ext_thin_clad_frac_400->Fit("pol1");
  tg_ext_thin_clad_frac_400->Draw("AP");
  TF1 *fit1 = (TF1*)tg_ext_thin_clad_frac_400->GetListOfFunctions()->FindObject("pol1");
  clad_nm_slopes.push_back(fit1->GetParameter(1));
  clad_nm_slopes_err.push_back(fit1->GetParError(1));
  page[1]->Update();
  ofile->cd();
  page[1]->Write();
  page[1]->Close();

  page[2]->cd();
  tg_ext_thin_poly_frac_410 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_410_ext_thin_poly_a, aging_empty_xbins_err, aging_410_ext_thin_poly_err_a);
  tg_ext_thin_poly_frac_410->SetTitle("Average Poly-side R at 410nm of New Ext vs Time");
  tg_ext_thin_poly_frac_410->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_poly_frac_410->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_poly_frac_410->SetMarkerStyle(20);
  tg_ext_thin_poly_frac_410->Fit("pol1");
  tg_ext_thin_poly_frac_410->Draw("AP");
  TF1 *fit2 = (TF1*)tg_ext_thin_poly_frac_410->GetListOfFunctions()->FindObject("pol1");
  poly_nm_slopes.push_back(fit2->GetParameter(1));
  poly_nm_slopes_err.push_back(fit2->GetParError(1));
  page[2]->Update();
  ofile->cd();
  page[2]->Write();
  page[2]->Close();

  page[3]->cd();
  tg_ext_thin_clad_frac_410 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_410_ext_thin_clad_a, aging_empty_xbins_err, aging_410_ext_thin_clad_err_a);
  tg_ext_thin_clad_frac_410->SetTitle("Average Clad-side R at 410nm of New Ext vs Time");
  tg_ext_thin_clad_frac_410->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_clad_frac_410->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_clad_frac_410->SetMarkerStyle(20);
  tg_ext_thin_clad_frac_410->Fit("pol1");
  tg_ext_thin_clad_frac_410->Draw("AP");
  TF1 *fit3 = (TF1*)tg_ext_thin_clad_frac_410->GetListOfFunctions()->FindObject("pol1");
  clad_nm_slopes.push_back(fit3->GetParameter(1));
  clad_nm_slopes_err.push_back(fit3->GetParError(1));
  page[3]->Update();
  ofile->cd();
  page[3]->Write();
  page[3]->Close();

  page[4]->cd();
  tg_ext_thin_poly_frac_420 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_420_ext_thin_poly_a, aging_empty_xbins_err, aging_420_ext_thin_poly_err_a);
  tg_ext_thin_poly_frac_420->SetTitle("Average Poly-side R at 420nm of New Ext vs Time");
  tg_ext_thin_poly_frac_420->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_poly_frac_420->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_poly_frac_420->SetMarkerStyle(20);
  tg_ext_thin_poly_frac_420->Fit("pol1");
  tg_ext_thin_poly_frac_420->Draw("AP");
  TF1 *fit4 = (TF1*)tg_ext_thin_poly_frac_420->GetListOfFunctions()->FindObject("pol1");
  poly_nm_slopes.push_back(fit4->GetParameter(1));
  poly_nm_slopes_err.push_back(fit4->GetParError(1));
  page[4]->Update();
  ofile->cd();
  page[4]->Write();
  page[4]->Close();

  page[5]->cd();
  tg_ext_thin_clad_frac_420 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_420_ext_thin_clad_a, aging_empty_xbins_err, aging_420_ext_thin_clad_err_a);
  tg_ext_thin_clad_frac_420->SetTitle("Average Clad-side R at 420nm of New Ext vs Time");
  tg_ext_thin_clad_frac_420->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_clad_frac_420->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_clad_frac_420->SetMarkerStyle(20);
  tg_ext_thin_clad_frac_420->Fit("pol1");
  tg_ext_thin_clad_frac_420->Draw("AP");
  TF1 *fit5 = (TF1*)tg_ext_thin_clad_frac_420->GetListOfFunctions()->FindObject("pol1");
  clad_nm_slopes.push_back(fit5->GetParameter(1));
  clad_nm_slopes_err.push_back(fit5->GetParError(1));
  page[5]->Update();
  ofile->cd();
  page[5]->Write();
  page[5]->Close();

  page[6]->cd();
  tg_ext_thin_poly_frac_430 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_430_ext_thin_poly_a, aging_empty_xbins_err, aging_430_ext_thin_poly_err_a);
  tg_ext_thin_poly_frac_430->SetTitle("Average Poly-side R at 430nm of New Ext vs Time");
  tg_ext_thin_poly_frac_430->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_poly_frac_430->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_poly_frac_430->SetMarkerStyle(20);
  tg_ext_thin_poly_frac_430->Fit("pol1");
  tg_ext_thin_poly_frac_430->Draw("AP");
  TF1 *fit6 = (TF1*)tg_ext_thin_poly_frac_430->GetListOfFunctions()->FindObject("pol1");
  poly_nm_slopes.push_back(fit6->GetParameter(1));
  poly_nm_slopes_err.push_back(fit6->GetParError(1));
  page[6]->Update();
  ofile->cd();
  page[6]->Write();
  page[6]->Close();

  page[7]->cd();
  tg_ext_thin_clad_frac_430 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_430_ext_thin_clad_a, aging_empty_xbins_err, aging_430_ext_thin_clad_err_a);
  tg_ext_thin_clad_frac_430->SetTitle("Average Clad-side R at 430nm of New Ext vs Time");
  tg_ext_thin_clad_frac_430->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_clad_frac_430->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_clad_frac_430->SetMarkerStyle(20);
  tg_ext_thin_clad_frac_430->Fit("pol1");
  tg_ext_thin_clad_frac_430->Draw("AP");
  TF1 *fit7 = (TF1*)tg_ext_thin_clad_frac_430->GetListOfFunctions()->FindObject("pol1");
  clad_nm_slopes.push_back(fit7->GetParameter(1));
  clad_nm_slopes_err.push_back(fit7->GetParError(1));
  page[7]->Update();
  ofile->cd();
  page[7]->Write();
  page[7]->Close();

  page[8]->cd();
  tg_ext_thin_poly_frac_440 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_440_ext_thin_poly_a, aging_empty_xbins_err, aging_440_ext_thin_poly_err_a);
  tg_ext_thin_poly_frac_440->SetTitle("Average Poly-side R at 440nm of New Ext vs Time");
  tg_ext_thin_poly_frac_440->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_poly_frac_440->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_poly_frac_440->SetMarkerStyle(20);
  tg_ext_thin_poly_frac_440->Fit("pol1");
  tg_ext_thin_poly_frac_440->Draw("AP");
  TF1 *fit8 = (TF1*)tg_ext_thin_poly_frac_440->GetListOfFunctions()->FindObject("pol1");
  poly_nm_slopes.push_back(fit8->GetParameter(1));
  poly_nm_slopes_err.push_back(fit8->GetParError(1));
  page[8]->Update();
  ofile->cd();
  page[8]->Write();
  page[8]->Close();

  page[9]->cd();
  tg_ext_thin_clad_frac_440 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_440_ext_thin_clad_a, aging_empty_xbins_err, aging_440_ext_thin_clad_err_a);
  tg_ext_thin_clad_frac_440->SetTitle("Average Clad-side R at 440nm of New Ext vs Time");
  tg_ext_thin_clad_frac_440->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_clad_frac_440->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_clad_frac_440->SetMarkerStyle(20);
  tg_ext_thin_clad_frac_440->Fit("pol1");
  tg_ext_thin_clad_frac_440->Draw("AP");
  TF1 *fit9 = (TF1*)tg_ext_thin_clad_frac_440->GetListOfFunctions()->FindObject("pol1");
  clad_nm_slopes.push_back(fit9->GetParameter(1));
  clad_nm_slopes_err.push_back(fit9->GetParError(1));
  page[9]->Update();
  ofile->cd();
  page[9]->Write();
  page[9]->Close();

  page[10]->cd();
  tg_ext_thin_poly_frac_450 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_450_ext_thin_poly_a, aging_empty_xbins_err, aging_450_ext_thin_poly_err_a);
  tg_ext_thin_poly_frac_450->SetTitle("Average Poly-side R at 450nm of New Ext vs Time");
  tg_ext_thin_poly_frac_450->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_poly_frac_450->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_poly_frac_450->SetMarkerStyle(20);
  tg_ext_thin_poly_frac_450->Fit("pol1");
  tg_ext_thin_poly_frac_450->Draw("AP");
  TF1 *fit10 = (TF1*)tg_ext_thin_poly_frac_450->GetListOfFunctions()->FindObject("pol1");
  poly_nm_slopes.push_back(fit10->GetParameter(1));
  poly_nm_slopes_err.push_back(fit10->GetParError(1));
  page[10]->Update();
  ofile->cd();
  page[10]->Write();
  page[10]->Close();

  page[11]->cd();
  tg_ext_thin_clad_frac_450 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_450_ext_thin_clad_a, aging_empty_xbins_err, aging_450_ext_thin_clad_err_a);
  tg_ext_thin_clad_frac_450->SetTitle("Average Clad-side R at 450nm of New Ext vs Time");
  tg_ext_thin_clad_frac_450->GetXaxis()->SetTitle("Time Fraction (% of Years Since 9/6/22)");
  tg_ext_thin_clad_frac_450->GetYaxis()->SetTitle("Reflectance (R) (%)");
  tg_ext_thin_clad_frac_450->SetMarkerStyle(20);
  tg_ext_thin_clad_frac_450->Fit("pol1");
  tg_ext_thin_clad_frac_450->Draw("AP");
  TF1 *fit11 = (TF1*)tg_ext_thin_clad_frac_450->GetListOfFunctions()->FindObject("pol1");
  clad_nm_slopes.push_back(fit11->GetParameter(1));
  clad_nm_slopes_err.push_back(fit11->GetParError(1));
  page[11]->Update();
  ofile->cd();
  page[11]->Write();
  page[11]->Close();

 // slopes as a function of wavelength

  TCanvas* slopes = new TCanvas("aging_slope_vs_wavelength","aging_slope_vs_wavelength", 900, 600);

  TH1F* hslopes1D_poly;
  TH1F* hslopes1D_clad;
  hslopes1D_poly = new TH1F("hslopes1D_poly", "hslopes1D_poly", 100, -3., 3.);
  hslopes1D_poly->GetXaxis()->SetTitle("Poly-side Slopes (Percent of R Lost per Year)");
  hslopes1D_poly->GetYaxis()->SetTitle("Counts");

  hslopes1D_clad = new TH1F("hslopes1D_clad", "hslopes1D_clad", 100, -3., 3.);
  hslopes1D_clad->SetTitle("Extrusion Sample Aging Slopes, Cladding side");
  hslopes1D_clad->GetXaxis()->SetTitle("Clad-side Slopes (R% Lost per Year)");
  hslopes1D_clad->GetYaxis()->SetTitle("Counts");

  Float_t poly_nm_slopes_a[8];
  Float_t clad_nm_slopes_a[8];
  std::copy(poly_nm_slopes.begin(), poly_nm_slopes.end(), poly_nm_slopes_a);
  std::copy(clad_nm_slopes.begin(), clad_nm_slopes.end(), clad_nm_slopes_a);
  Float_t poly_nm_slopes_err_a[8];
  Float_t clad_nm_slopes_err_a[8];
  std::copy(poly_nm_slopes_err.begin(), poly_nm_slopes_err.end(), poly_nm_slopes_err_a);
  std::copy(clad_nm_slopes_err.begin(), clad_nm_slopes_err.end(), clad_nm_slopes_err_a);
  Float_t wavelengths[8] = { 380., 390., 400., 410., 420., 430., 440., 450. };
  Float_t wavelength_err[8] = {0};

  for(int i = 0; i < 8; i++){
    hslopes1D_poly->Fill(poly_nm_slopes_a[i]);
    hslopes1D_clad->Fill(clad_nm_slopes_a[i]);
  }
  hslopes1D_poly->Write();
  hslopes1D_clad->Write();

  slopes->cd();
  TGraphErrors* tg_poly_nm_slopes = new TGraphErrors(8, wavelengths, poly_nm_slopes_a, wavelength_err, poly_nm_slopes_err_a);
  tg_poly_nm_slopes->SetMarkerStyle(24);
  TGraphErrors* tg_clad_nm_slopes = new TGraphErrors(8, wavelengths, clad_nm_slopes_a, wavelength_err, clad_nm_slopes_err_a);
  tg_clad_nm_slopes->SetTitle("Aging Slopes vs Wavelength, New Extrusions");
  tg_clad_nm_slopes->GetXaxis()->SetTitle("Wavelength (nm)");
  tg_clad_nm_slopes->GetYaxis()->SetTitle("Aging Slope (R% Loss/Year)");
  tg_clad_nm_slopes->GetYaxis()->SetRangeUser(-2,2);
  tg_clad_nm_slopes->SetMarkerStyle(20);
  //tg_clad_nm_slopes->SetMarkerColor(4);
  //tg_clad_nm_slopes->SetLineColor(4);
  tg_poly_nm_slopes->SetMarkerColor(2);
  tg_poly_nm_slopes->SetLineColor(2);
  tg_clad_nm_slopes->Draw("AP");
  //tg_poly_nm_slopes->Draw("P");

  TLegend* slope_leg = new TLegend(0.12, 0.13, 0.4, 0.2);
  slope_leg->AddEntry(tg_clad_nm_slopes,"Cladding Side Slopes", "lep");
  slope_leg->AddEntry(tg_poly_nm_slopes,"Polymer Side Slopes", "lep");
  //slope_leg->Draw();

  slopes->Update();
  ofile->cd();
  slopes->Write();
  slopes->Close();

  stability[0]->cd();
  tg_poly_stab = new TGraph(8, wavelengths, stab_ext_thin_poly_a);
  tg_poly_stab->SetMarkerStyle(20);
  tg_poly_stab->SetTitle("Standard Deviation of Repeated Polymer-Side Measurements");
  tg_poly_stab->GetXaxis()->SetTitle("Wavelength (nm)");
  tg_poly_stab->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  tg_poly_stab->Draw("AP");
  stability[0]->Update();
  ofile->cd();
  stability[0]->Write();
  stability[0]->Close();

  stability[1]->cd();
  tg_clad_stab = new TGraph(8, wavelengths, stab_ext_thin_clad_a);
  tg_clad_stab->SetMarkerStyle(20);
  tg_clad_stab->SetTitle("Standard Deviation of Repeated Cladding-Side Measurements");
  tg_clad_stab->GetXaxis()->SetTitle("Wavelength (nm)");
  tg_clad_stab->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  tg_clad_stab->Draw("AP");
  stability[1]->Update();
  ofile->cd();
  stability[1]->Write();
  stability[1]->Close();

  // Close the output file
  ofile->Close();
}
