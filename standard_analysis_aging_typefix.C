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
  if(nm == 600) bin = 25;

  for(int ind = 0;ind < indices.size();ind++){
    int hist = indices[ind];
    these_points.push_back(hists[hist]->GetBinContent(bin)); // get relevant wavelength bin and store 
    //if(nm == 600) std::cout << "600nm points are: " << hists[hist]->GetBinContent(bin) << std::endl;
  }
  //if(nm == 600) std::cout << "\n" << std::endl;
  float stdev = TMath::StdDev(these_points.begin(), these_points.end());
  return stdev;
}


std::vector<int> slicing(std::vector<int>& vec, int x, int y){

  auto start = vec.begin() + x;
  auto end = vec.begin() + y + 1;

  std::vector<int> result(y - x + 1);

  std::copy(start, end, result.begin());

  return result;
}


void standard_analysis_aging_typefix(const char* filename){

  TFile* file = TFile::Open(Form("%s", filename));
  if (file->IsZombie()) {
    std::cout << "Problem opening file " << filename << std::endl;
    return;
  }

  //Create output file. If it already exists, recreate it.
  TFile *ofile = new TFile("output_TiO2_standards_analysis_aging_typefix_0803_no0524_goreerr.root","RECREATE");

  // define histograms
  TH1F* brian_92316_avg  = new TH1F("brian_92316_avg", "brian_92316_avg",43,360.,780.);
  TH1F* brian_111516_avg  = new TH1F("brian_111516_avg", "brian_111516_avg",43,360.,780.);
  TH1F* brian_11618_avg  = new TH1F("brian_11618_avg", "brian_11618_avg",43,360.,780.);
  TH1F* brian_40318_avg  = new TH1F("brian_40318_avg", "brian_40318_avg",43,360.,780.);
  TH1F* brian_110420_avg  = new TH1F("brian_110420_avg", "brian_110420_avg",43,360.,780.);

  TH1F* mackenzie_920_avg  = new TH1F("mackenzie_920_avg", "mackenzie_920_avg",43,360.,780.);
  TH1F* mackenzie_1011_bef_avg  = new TH1F("mackenzie_1011_bef_avg", "mackenzie_1011_bef_avg",43,360.,780.);
  TH1F* mackenzie_1011_aft_avg  = new TH1F("mackenzie_1011_aft_avg", "mackenzie_1011_aft_avg",43,360.,780.);
  TH1F* mackenzie_1018_bef_avg  = new TH1F("mackenzie_1018_bef_avg", "mackenzie_1018_bef_avg",43,360.,780.);
  TH1F* mackenzie_1018_aft_avg  = new TH1F("mackenzie_1018_aft_avg", "mackenzie_1018_aft_avg",43,360.,780.);
  TH1F* mackenzie_1025_bef_avg  = new TH1F("mackenzie_1025_bef_avg", "mackenzie_1025_bef_avg",43,360.,780.);
  TH1F* mackenzie_1025_aft_avg  = new TH1F("mackenzie_1025_aft_avg", "mackenzie_1025_aft_avg",43,360.,780.);
  TH1F* mackenzie_1101_bef_avg  = new TH1F("mackenzie_1101_bef_avg", "mackenzie_1101_bef_avg",43,360.,780.);
  TH1F* mackenzie_1101_aft_avg  = new TH1F("mackenzie_1101_aft_avg", "mackenzie_1101_aft_avg",43,360.,780.);
  TH1F* mackenzie_1103_bef_avg  = new TH1F("mackenzie_1103_bef_avg", "mackenzie_1103_bef_avg",43,360.,780.);
  TH1F* mackenzie_1103_aft_avg  = new TH1F("mackenzie_1103_aft_avg", "mackenzie_1103_aft_avg",43,360.,780.);
  TH1F* mackenzie_1108_bef_avg  = new TH1F("mackenzie_1108_bef_avg", "mackenzie_1108_bef_avg",43,360.,780.);
  TH1F* mackenzie_1108_aft_avg  = new TH1F("mackenzie_1108_aft_avg", "mackenzie_1108_aft_avg",43,360.,780.);
  TH1F* mackenzie_1115_bef_avg  = new TH1F("mackenzie_1115_bef_avg", "mackenzie_1115_bef_avg",43,360.,780.);
  TH1F* mackenzie_1115_aft_avg  = new TH1F("mackenzie_1115_aft_avg", "mackenzie_1115_aft_avg",43,360.,780.);
  TH1F* mackenzie_1122_bef_avg  = new TH1F("mackenzie_1122_bef_avg", "mackenzie_1122_bef_avg",43,360.,780.);
  TH1F* mackenzie_1122_aft_avg  = new TH1F("mackenzie_1122_aft_avg", "mackenzie_1122_aft_avg",43,360.,780.);
  TH1F* mackenzie_1129_bef_avg  = new TH1F("mackenzie_1129_bef_avg", "mackenzie_1129_bef_avg",43,360.,780.);
  TH1F* mackenzie_1129_aft_avg  = new TH1F("mackenzie_1129_aft_avg", "mackenzie_1129_aft_avg",43,360.,780.);
  TH1F* mackenzie_1206_bef_avg  = new TH1F("mackenzie_1206_bef_avg", "mackenzie_1206_bef_avg",43,360.,780.);
  TH1F* mackenzie_1206_aft_avg  = new TH1F("mackenzie_1206_aft_avg", "mackenzie_1206_aft_avg",43,360.,780.);
  TH1F* brian_020123_bef_avg  = new TH1F("brian_020123_bef_avg", "brian_020123_bef_avg",43,360.,780.);
  TH1F* brian_020123_aft_avg  = new TH1F("brian_020123_aft_avg", "brian_020123_aft_avg",43,360.,780.);
  TH1F* brian_020923_bef_avg  = new TH1F("brian_020923_bef_avg", "brian_020923_bef_avg",43,360.,780.);
  TH1F* brian_020923_aft_avg  = new TH1F("brian_020923_aft_avg", "brian_020923_aft_avg",43,360.,780.);
  TH1F* brian_021723_bef_avg  = new TH1F("brian_021723_bef_avg", "brian_021723_bef_avg",43,360.,780.);
  TH1F* brian_021723_aft_avg  = new TH1F("brian_021723_aft_avg", "brian_021723_aft_avg",43,360.,780.);
  TH1F* brian_030923_bef_avg  = new TH1F("brian_030923_bef_avg", "brian_030923_bef_avg",43,360.,780.);
  TH1F* brian_030923_aft_avg  = new TH1F("brian_030923_aft_avg", "brian_030923_aft_avg",43,360.,780.);
  TH1F* brian_031623_bef_avg  = new TH1F("brian_031623_bef_avg", "brian_031623_bef_avg",43,360.,780.);
  TH1F* brian_031623_aft_avg  = new TH1F("brian_031623_aft_avg", "brian_031623_aft_avg",43,360.,780.);
  TH1F* brian_032323_bef_avg  = new TH1F("brian_032323_bef_avg", "brian_032323_bef_avg",43,360.,780.);
  TH1F* brian_032323_aft_avg  = new TH1F("brian_032323_aft_avg", "brian_032323_aft_avg",43,360.,780.);
  TH1F* brian_040623_bef_avg  = new TH1F("brian_040623_bef_avg", "brian_040623_bef_avg",43,360.,780.);
  TH1F* brian_040623_aft_avg  = new TH1F("brian_040623_aft_avg", "brian_040623_aft_avg",43,360.,780.);
  TH1F* brian_041323_bef_avg  = new TH1F("brian_041323_bef_avg", "brian_041323_bef_avg",43,360.,780.);
  TH1F* brian_041323_aft_avg  = new TH1F("brian_041323_aft_avg", "brian_041323_aft_avg",43,360.,780.);
  TH1F* brian_042123_bef_avg  = new TH1F("brian_042123_bef_avg", "brian_042123_bef_avg",43,360.,780.);
  TH1F* brian_042123_aft_avg  = new TH1F("brian_042123_aft_avg", "brian_042123_aft_avg",43,360.,780.);
  TH1F* brian_052423_bef_avg  = new TH1F("brian_052423_bef_avg", "brian_052423_bef_avg",43,360.,780.);
  TH1F* brian_052423_aft_avg  = new TH1F("brian_052423_aft_avg", "brian_052423_aft_avg",43,360.,780.);
  TH1F* brian_060823_bef_avg  = new TH1F("brian_060823_bef_avg", "brian_060823_bef_avg",43,360.,780.);
  TH1F* brian_060823_aft_avg  = new TH1F("brian_060823_aft_avg", "brian_060823_aft_avg",43,360.,780.);
  TH1F* brian_062323_bef_avg  = new TH1F("brian_062323_bef_avg", "brian_062323_bef_avg",43,360.,780.);
  TH1F* brian_062323_aft_avg  = new TH1F("brian_062323_aft_avg", "brian_062323_aft_avg",43,360.,780.);

  TH1F* mackenzie_920_gore_avg       = new TH1F("mackenzie_920_gore_avg", "mackenzie_920_gore_avg",43,360.,780.);
  TH1F* mackenzie_1018_gore_avg      = new TH1F("mackenzie_1018_gore_avg", "mackenzie_1018_gore_avg",43,360.,780.);
  TH1F* mackenzie_1025_gore_bef_avg  = new TH1F("mackenzie_1025_gore_bef_avg", "mackenzie_1025_gore_bef_avg",43,360.,780.);
  TH1F* mackenzie_1025_gore_aft_avg  = new TH1F("mackenzie_1025_gore_aft_avg", "mackenzie_1025_gore_aft_avg",43,360.,780.);
  TH1F* mackenzie_1101_gore_bef_avg  = new TH1F("mackenzie_1101_gore_bef_avg", "mackenzie_1101_gore_bef_avg",43,360.,780.);
  TH1F* mackenzie_1101_gore_aft_avg  = new TH1F("mackenzie_1101_gore_aft_avg", "mackenzie_1101_gore_aft_avg",43,360.,780.);
  TH1F* mackenzie_1108_gore_avg      = new TH1F("mackenzie_1108_gore_avg", "mackenzie_1108_gore_avg",43,360.,780.);
  TH1F* mackenzie_1115_gore_bef_avg  = new TH1F("mackenzie_1115_gore_bef_avg", "mackenzie_1115_gore_bef_avg",43,360.,780.);
  TH1F* mackenzie_1115_gore_aft_avg  = new TH1F("mackenzie_1115_gore_aft_avg", "mackenzie_1115_gore_aft_avg",43,360.,780.);
  TH1F* mackenzie_1122_gore_bef_avg  = new TH1F("mackenzie_1122_gore_bef_avg", "mackenzie_1122_gore_bef_avg",43,360.,780.);
  TH1F* mackenzie_1122_gore_aft_avg  = new TH1F("mackenzie_1122_gore_aft_avg", "mackenzie_1122_gore_aft_avg",43,360.,780.);
  TH1F* mackenzie_1129_gore_bef_avg  = new TH1F("mackenzie_1129_gore_bef_avg", "mackenzie_1129_gore_bef_avg",43,360.,780.);
  TH1F* mackenzie_1129_gore_aft_avg  = new TH1F("mackenzie_1129_gore_aft_avg", "mackenzie_1129_gore_aft_avg",43,360.,780.);
  TH1F* mackenzie_1206_gore_bef_avg  = new TH1F("mackenzie_1206_gore_bef_avg", "mackenzie_1206_gore_bef_avg",43,360.,780.);
  TH1F* mackenzie_1206_gore_aft_avg  = new TH1F("mackenzie_1206_gore_aft_avg", "mackenzie_1206_gore_aft_avg",43,360.,780.);
  TH1F* brian_020123_gore_bef_avg  = new TH1F("brian_020123_gore_bef_avg", "brian_020123_gore_bef_avg",43,360.,780.);
  TH1F* brian_020123_gore_aft_avg  = new TH1F("brian_020123_gore_aft_avg", "brian_020123_gore_aft_avg",43,360.,780.);
  TH1F* brian_020923_gore_bef_avg  = new TH1F("brian_020923_gore_bef_avg", "brian_020923_gore_bef_avg",43,360.,780.);
  TH1F* brian_020923_gore_aft_avg  = new TH1F("brian_020923_gore_aft_avg", "brian_020923_gore_aft_avg",43,360.,780.);
  TH1F* brian_021723_gore_bef_avg  = new TH1F("brian_021723_gore_bef_avg", "brian_021723_gore_bef_avg",43,360.,780.);
  TH1F* brian_021723_gore_aft_avg  = new TH1F("brian_021723_gore_aft_avg", "brian_021723_gore_aft_avg",43,360.,780.);
  TH1F* brian_030923_gore_bef_avg  = new TH1F("brian_030923_gore_bef_avg", "brian_030923_gore_bef_avg",43,360.,780.);
  TH1F* brian_030923_gore_aft_avg  = new TH1F("brian_030923_gore_aft_avg", "brian_030923_gore_aft_avg",43,360.,780.);
  TH1F* brian_031623_gore_bef_avg  = new TH1F("brian_031623_gore_bef_avg", "brian_031623_gore_bef_avg",43,360.,780.);
  TH1F* brian_031623_gore_aft_avg  = new TH1F("brian_031623_gore_aft_avg", "brian_031623_gore_aft_avg",43,360.,780.);
  TH1F* brian_032323_gore_bef_avg  = new TH1F("brian_032323_gore_bef_avg", "brian_032323_gore_bef_avg",43,360.,780.);
  TH1F* brian_032323_gore_aft_avg  = new TH1F("brian_032323_gore_aft_avg", "brian_032323_gore_aft_avg",43,360.,780.);
  TH1F* brian_040623_gore_bef_avg  = new TH1F("brian_040623_gore_bef_avg", "brian_040623_gore_bef_avg",43,360.,780.);
  TH1F* brian_040623_gore_aft_avg  = new TH1F("brian_040623_gore_aft_avg", "brian_040623_gore_aft_avg",43,360.,780.);
  TH1F* brian_041323_gore_bef_avg  = new TH1F("brian_041323_gore_bef_avg", "brian_041323_gore_bef_avg",43,360.,780.);
  TH1F* brian_041323_gore_aft_avg  = new TH1F("brian_041323_gore_aft_avg", "brian_041323_gore_aft_avg",43,360.,780.);
  TH1F* brian_042123_gore_bef_avg  = new TH1F("brian_042123_gore_bef_avg", "brian_042123_gore_bef_avg",43,360.,780.);
  TH1F* brian_042123_gore_aft_avg  = new TH1F("brian_042123_gore_aft_avg", "brian_042123_gore_aft_avg",43,360.,780.);
  TH1F* brian_052423_gore_bef_avg  = new TH1F("brian_052423_gore_bef_avg", "brian_052423_gore_bef_avg",43,360.,780.);
  TH1F* brian_052423_gore_aft_avg  = new TH1F("brian_052423_gore_aft_avg", "brian_052423_gore_aft_avg",43,360.,780.);
  TH1F* brian_060823_gore_bef_avg  = new TH1F("brian_060823_gore_bef_avg", "brian_060823_gore_bef_avg",43,360.,780.);
  TH1F* brian_060823_gore_aft_avg  = new TH1F("brian_060823_gore_aft_avg", "brian_060823_gore_aft_avg",43,360.,780.);
  TH1F* brian_062323_gore_bef_avg  = new TH1F("brian_062323_gore_bef_avg", "brian_062323_gore_bef_avg",43,360.,780.);
  TH1F* brian_062323_gore_aft_avg  = new TH1F("brian_062323_gore_aft_avg", "brian_062323_gore_aft_avg",43,360.,780.);

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
  std::vector<int> avg_brian_92316;
  std::vector<int> avg_brian_111516;
  std::vector<int> avg_brian_40318;
  std::vector<int> avg_brian_11618;
  std::vector<int> avg_brian_110420;

  std::vector<int> mackenzie_nova_920;
  std::vector<int> mackenzie_nova_1011;
  std::vector<int> mackenzie_nova_1018;
  std::vector<int> mackenzie_nova_1025;
  std::vector<int> mackenzie_nova_1101;
  std::vector<int> mackenzie_nova_1103;
  std::vector<int> mackenzie_nova_1108;
  std::vector<int> mackenzie_nova_1115;
  std::vector<int> mackenzie_nova_1122;
  std::vector<int> mackenzie_nova_1129;
  std::vector<int> mackenzie_nova_1206;
  std::vector<int> brian_nova_020123;
  std::vector<int> brian_nova_020923;
  std::vector<int> brian_nova_021723;
  std::vector<int> brian_nova_030923;
  std::vector<int> brian_nova_031623;
  std::vector<int> brian_nova_032323;
  std::vector<int> brian_nova_040623;
  std::vector<int> brian_nova_041323;
  std::vector<int> brian_nova_042123;
  std::vector<int> brian_nova_052423;
  std::vector<int> brian_nova_060823;
  std::vector<int> brian_nova_062323;

  std::vector<int> mackenzie_goreDRP_920;
  std::vector<int> mackenzie_goreDRP_1018;
  std::vector<int> mackenzie_goreDRP_1025;
  std::vector<int> mackenzie_goreDRP_1101;
  std::vector<int> mackenzie_goreDRP_1108;
  std::vector<int> mackenzie_goreDRP_1115;
  std::vector<int> mackenzie_goreDRP_1122;
  std::vector<int> mackenzie_goreDRP_1129;
  std::vector<int> mackenzie_goreDRP_1206;
  std::vector<int> brian_goreDRP_020123;
  std::vector<int> brian_goreDRP_020923;
  std::vector<int> brian_goreDRP_021723;
  std::vector<int> brian_goreDRP_030923;
  std::vector<int> brian_goreDRP_031623;
  std::vector<int> brian_goreDRP_032323;
  std::vector<int> brian_goreDRP_040623;
  std::vector<int> brian_goreDRP_041323;
  std::vector<int> brian_goreDRP_042123;
  std::vector<int> brian_goreDRP_052423;
  std::vector<int> brian_goreDRP_060823;
  std::vector<int> brian_goreDRP_062323;

  std::vector<int> ordered_aging_indices;
  int mackenzie_317;
  int mackenzie_324;
  int mackenzie_407;
  int mackenzie_414;
  int mackenzie_421;
  int mackenzie_428;
  int mackenzie_505;
  int mackenzie_512;
  int mackenzie_526;
  int mackenzie_602;
  int mackenzie_616;
  int mackenzie_714;
  int mackenzie_728;
  int mackenzie_811;
  int mackenzie_825;
  int mackenzie_908;
  int mackenzie_913;
  int mackenzie_920;
  int mackenzie_927;
  int mackenzie_1004;
  int mackenzie_goreDRP_831;
  int mackenzie_goreDRP_901;
  int mackenzie_goreDRP_927;
  int mackenzie_goreDRP_1011;

  // here is where we actually do the grouping and push indices into the above vectors
  int index = 0;
  for(std::string& s : key_strs){
    if(s.find("brian_92316") != std::string::npos){
      avg_brian_92316.push_back(index);
    }
    if(s.find("brian_111516") != std::string::npos){
      avg_brian_111516.push_back(index);
    }
    if(s.find("brian_40318") != std::string::npos){
      avg_brian_40318.push_back(index);
    }
    if(s.find("brian_11618") != std::string::npos){
      avg_brian_11618.push_back(index);
    }
    if(s.find("brian_110420") != std::string::npos){
      avg_brian_110420.push_back(index);
    }

    if(s.find("mackenzie_317") != std::string::npos){
      mackenzie_317 = index;
    }
    if(s.find("mackenzie_324") != std::string::npos){
      mackenzie_324 = index;
    }
    if(s.find("mackenzie_407") != std::string::npos){
      mackenzie_407 = index;
    }
    if(s.find("mackenzie_414") != std::string::npos){
      mackenzie_414 = index;
    }
    if(s.find("mackenzie_421") != std::string::npos){
      mackenzie_421 = index;
    }
    if(s.find("mackenzie_428") != std::string::npos){
      mackenzie_428 = index;
    }
    if(s.find("mackenzie_505") != std::string::npos){
      mackenzie_505 = index;
    }
    if(s.find("mackenzie_512") != std::string::npos){
      mackenzie_512 = index;
    }
    if(s.find("mackenzie_526") != std::string::npos){
      mackenzie_526 = index;
    }
    if(s.find("mackenzie_602") != std::string::npos){
      mackenzie_602 = index;
    }
    if(s.find("mackenzie_616") != std::string::npos){
      mackenzie_616 = index;
    }
    if(s.find("mackenzie_714") != std::string::npos){
      mackenzie_714 = index;
    }
    if(s.find("mackenzie_728") != std::string::npos){
      mackenzie_728 = index;
    }
    if(s.find("mackenzie_811") != std::string::npos){
      mackenzie_811 = index;
    }
    if(s.find("mackenzie_825") != std::string::npos){
      mackenzie_825 = index;
    }
    if(s.find("mackenzie_908") != std::string::npos){
      mackenzie_908 = index;
    }
    if(s.find("mackenzie_913") != std::string::npos){
      mackenzie_913 = index;
    }
    if(s.find("mackenzie_920_0") != std::string::npos){
      mackenzie_920 = index;
    }
    if(s.find("mackenzie_927") != std::string::npos){
      mackenzie_927 = index;
    }
    if(s.find("mackenzie_1004") != std::string::npos){
      mackenzie_1004 = index;
    }

    if(s.find("mackenzie_920") != std::string::npos){
      mackenzie_nova_920.push_back(index);
    }
    if(s.find("mackenzie_1011") != std::string::npos){
      mackenzie_nova_1011.push_back(index);
    }
    if(s.find("mackenzie_1018") != std::string::npos){
      mackenzie_nova_1018.push_back(index);
    }
    if(s.find("mackenzie_1025") != std::string::npos){
      mackenzie_nova_1025.push_back(index);
    }
    if(s.find("mackenzie_1101") != std::string::npos){
      mackenzie_nova_1101.push_back(index);
    }
    if(s.find("mackenzie_1103") != std::string::npos){
      mackenzie_nova_1103.push_back(index);
    }
    if(s.find("mackenzie_1108") != std::string::npos){
      mackenzie_nova_1108.push_back(index);
    }
    if(s.find("mackenzie_1115") != std::string::npos){
      mackenzie_nova_1115.push_back(index);
    }
    if(s.find("mackenzie_1122") != std::string::npos){
      mackenzie_nova_1122.push_back(index);
    }
    if(s.find("mackenzie_1129") != std::string::npos){
      mackenzie_nova_1129.push_back(index);
    }
    if(s.find("mackenzie_1206") != std::string::npos){
      mackenzie_nova_1206.push_back(index);
    }
    if(s.find("brian_020123") != std::string::npos){
      brian_nova_020123.push_back(index);
    }
    if(s.find("brian_020923") != std::string::npos){
      brian_nova_020923.push_back(index);
    }
    if(s.find("brian_021723") != std::string::npos){
      brian_nova_021723.push_back(index);
    }
    if(s.find("brian_030923") != std::string::npos){
      brian_nova_030923.push_back(index);
    }
    if(s.find("brian_031623") != std::string::npos){
      brian_nova_031623.push_back(index);
    }
    if(s.find("brian_032323") != std::string::npos){
      brian_nova_032323.push_back(index);
    }
    if(s.find("brian_040623") != std::string::npos){
      brian_nova_040623.push_back(index);
    }
    if(s.find("brian_041323") != std::string::npos){
      brian_nova_041323.push_back(index);
    }
    if(s.find("brian_042123") != std::string::npos){
      brian_nova_042123.push_back(index);
    }
    if(s.find("brian_052423") != std::string::npos){
      brian_nova_052423.push_back(index);
    }
    if(s.find("brian_060823") != std::string::npos){
      brian_nova_060823.push_back(index);
    }
    if(s.find("brian_062323") != std::string::npos){
      brian_nova_062323.push_back(index);
    }


    if(s.find("goreDRP_831") != std::string::npos){
      mackenzie_goreDRP_831 = index;
    }
    if(s.find("goreDRP_901") != std::string::npos){
      mackenzie_goreDRP_901 = index;
    }
    if(s.find("goreDRP_920") != std::string::npos){
      mackenzie_goreDRP_920.push_back(index);
    }
    if(s.find("goreDRP_927") != std::string::npos){
      mackenzie_goreDRP_927 = index;
    }
    if(s.find("goreDRP_1011") != std::string::npos){
      mackenzie_goreDRP_1011 = index;
    }
    if(s.find("goreDRP_1018") != std::string::npos){
      mackenzie_goreDRP_1018.push_back(index);
    }
    if(s.find("goreDRP_1025") != std::string::npos){
      mackenzie_goreDRP_1025.push_back(index);
    }
    if(s.find("goreDRP_1101") != std::string::npos){
      mackenzie_goreDRP_1101.push_back(index);
    }
    if(s.find("goreDRP_1108") != std::string::npos){
      mackenzie_goreDRP_1108.push_back(index);
    }
    if(s.find("goreDRP_1115") != std::string::npos){
      mackenzie_goreDRP_1115.push_back(index);
    }
    if(s.find("goreDRP_1122") != std::string::npos){
      mackenzie_goreDRP_1122.push_back(index);
    }
    if(s.find("goreDRP_1129") != std::string::npos){
      mackenzie_goreDRP_1129.push_back(index);
    }
    if(s.find("goreDRP_1206") != std::string::npos){
      mackenzie_goreDRP_1206.push_back(index);
    }
    if(s.find("goreDRP_020123") != std::string::npos){
      brian_goreDRP_020123.push_back(index);
    }
    if(s.find("goreDRP_020923") != std::string::npos){
      brian_goreDRP_020923.push_back(index);
    }
    if(s.find("goreDRP_021723") != std::string::npos){
      brian_goreDRP_021723.push_back(index);
    }
    if(s.find("goreDRP_030923") != std::string::npos){
      brian_goreDRP_030923.push_back(index);
    }
    if(s.find("goreDRP_031623") != std::string::npos){
      brian_goreDRP_031623.push_back(index);
    }
    if(s.find("goreDRP_032323") != std::string::npos){
      brian_goreDRP_032323.push_back(index);
    }
    if(s.find("goreDRP_040623") != std::string::npos){
      brian_goreDRP_040623.push_back(index);
    }
    if(s.find("goreDRP_041323") != std::string::npos){
      brian_goreDRP_041323.push_back(index);
    }
    if(s.find("goreDRP_042123") != std::string::npos){
      brian_goreDRP_042123.push_back(index);
    }
    if(s.find("goreDRP_052423") != std::string::npos){
      brian_goreDRP_052423.push_back(index);
    }
    if(s.find("goreDRP_060823") != std::string::npos){
      brian_goreDRP_060823.push_back(index);
    }
    if(s.find("goreDRP_062323") != std::string::npos){
      brian_goreDRP_062323.push_back(index);
    }

    ++index;
  }

  // make a vector of indices into input hists in order for aging
  // this list only has 1 histogram per date, deal with averages for >1
  ordered_aging_indices.push_back(mackenzie_317);
  ordered_aging_indices.push_back(mackenzie_324);
  ordered_aging_indices.push_back(mackenzie_407);
  ordered_aging_indices.push_back(mackenzie_414);
  ordered_aging_indices.push_back(mackenzie_421);
  ordered_aging_indices.push_back(mackenzie_428);
  ordered_aging_indices.push_back(mackenzie_505);
  ordered_aging_indices.push_back(mackenzie_512);
  ordered_aging_indices.push_back(mackenzie_526);
  ordered_aging_indices.push_back(mackenzie_602);
  ordered_aging_indices.push_back(mackenzie_616);
  ordered_aging_indices.push_back(mackenzie_714);
  ordered_aging_indices.push_back(mackenzie_728);
  ordered_aging_indices.push_back(mackenzie_811);
  ordered_aging_indices.push_back(mackenzie_825);
  ordered_aging_indices.push_back(mackenzie_908);
  ordered_aging_indices.push_back(mackenzie_913);
  ordered_aging_indices.push_back(mackenzie_920);
  ordered_aging_indices.push_back(mackenzie_927);
  ordered_aging_indices.push_back(mackenzie_1004);

  // now, do fanicier analysis here using groups from above. Averages, differences, etc.

  // average histogram function
  brian_92316_avg = average_hist(input_hists, avg_brian_92316, "avg_brian_92316");
  brian_111516_avg = average_hist(input_hists, avg_brian_111516, "avg_brian_111516");
  brian_11618_avg = average_hist(input_hists, avg_brian_11618, "avg_brian_11618");
  brian_40318_avg = average_hist(input_hists, avg_brian_40318, "avg_brian_40318");
  brian_110420_avg = average_hist(input_hists, avg_brian_110420, "avg_brian_110420");

  mackenzie_920_avg = average_hist(input_hists, mackenzie_nova_920, "avg_mackenzie_nova_920");
  mackenzie_1011_bef_avg = average_hist(input_hists, slicing(mackenzie_nova_1011,0,2), "avg_mackenzie_nova_1011_bef");
  mackenzie_1011_aft_avg = average_hist(input_hists, slicing(mackenzie_nova_1011,3,5), "avg_mackenzie_nova_1011_aft");
  mackenzie_1018_bef_avg = average_hist(input_hists, slicing(mackenzie_nova_1018,0,2), "avg_mackenzie_nova_1018_bef");
  mackenzie_1018_aft_avg = average_hist(input_hists, slicing(mackenzie_nova_1018,3,5), "avg_mackenzie_nova_1018_aft");
  mackenzie_1025_bef_avg = average_hist(input_hists, slicing(mackenzie_nova_1025,0,2), "avg_mackenzie_nova_1025_bef");
  mackenzie_1025_aft_avg = average_hist(input_hists, slicing(mackenzie_nova_1025,3,5), "avg_mackenzie_nova_1025_aft");
  mackenzie_1101_bef_avg = average_hist(input_hists, slicing(mackenzie_nova_1101,0,2), "avg_mackenzie_nova_1101_bef");
  mackenzie_1101_aft_avg = average_hist(input_hists, slicing(mackenzie_nova_1101,3,5), "avg_mackenzie_nova_1101_aft");
  mackenzie_1103_bef_avg = average_hist(input_hists, slicing(mackenzie_nova_1103,0,2), "avg_mackenzie_nova_1103_bef");
  mackenzie_1103_aft_avg = average_hist(input_hists, slicing(mackenzie_nova_1103,3,5), "avg_mackenzie_nova_1103_aft");
  mackenzie_1108_bef_avg = average_hist(input_hists, slicing(mackenzie_nova_1108,0,2), "avg_mackenzie_nova_1108_bef");
  mackenzie_1108_aft_avg = average_hist(input_hists, slicing(mackenzie_nova_1108,3,5), "avg_mackenzie_nova_1108_aft");
  mackenzie_1115_bef_avg = average_hist(input_hists, slicing(mackenzie_nova_1115,0,2), "avg_mackenzie_nova_1115_bef");
  mackenzie_1115_aft_avg = average_hist(input_hists, slicing(mackenzie_nova_1115,3,5), "avg_mackenzie_nova_1115_aft");
  mackenzie_1122_bef_avg = average_hist(input_hists, slicing(mackenzie_nova_1122,0,2), "avg_mackenzie_nova_1122_bef");
  mackenzie_1122_aft_avg = average_hist(input_hists, slicing(mackenzie_nova_1122,3,5), "avg_mackenzie_nova_1122_aft");
  mackenzie_1129_bef_avg = average_hist(input_hists, slicing(mackenzie_nova_1129,0,2), "avg_mackenzie_nova_1129_bef");
  mackenzie_1129_aft_avg = average_hist(input_hists, slicing(mackenzie_nova_1129,3,5), "avg_mackenzie_nova_1129_aft");
  mackenzie_1206_bef_avg = average_hist(input_hists, slicing(mackenzie_nova_1206,0,2), "avg_mackenzie_nova_1206_bef");
  mackenzie_1206_aft_avg = average_hist(input_hists, slicing(mackenzie_nova_1206,3,5), "avg_mackenzie_nova_1206_aft");
  brian_020123_bef_avg = average_hist(input_hists, slicing(brian_nova_020123,0,2), "avg_brian_nova_020123_bef");
  brian_020123_aft_avg = average_hist(input_hists, slicing(brian_nova_020123,3,5), "avg_brian_nova_020123_aft");
  brian_020923_bef_avg = average_hist(input_hists, slicing(brian_nova_020923,0,2), "avg_brian_nova_020923_bef");
  brian_020923_aft_avg = average_hist(input_hists, slicing(brian_nova_020923,3,5), "avg_brian_nova_020923_aft");
  brian_021723_bef_avg = average_hist(input_hists, slicing(brian_nova_021723,0,2), "avg_brian_nova_021723_bef");
  brian_021723_aft_avg = average_hist(input_hists, slicing(brian_nova_021723,3,5), "avg_brian_nova_021723_aft");
  brian_030923_bef_avg = average_hist(input_hists, slicing(brian_nova_030923,0,2), "avg_brian_nova_030923_bef");
  brian_030923_aft_avg = average_hist(input_hists, slicing(brian_nova_030923,3,5), "avg_brian_nova_030923_aft");
  brian_031623_bef_avg = average_hist(input_hists, slicing(brian_nova_031623,0,2), "avg_brian_nova_031623_bef");
  brian_031623_aft_avg = average_hist(input_hists, slicing(brian_nova_031623,3,5), "avg_brian_nova_031623_aft");
  brian_032323_bef_avg = average_hist(input_hists, slicing(brian_nova_032323,0,2), "avg_brian_nova_032323_bef");
  brian_032323_aft_avg = average_hist(input_hists, slicing(brian_nova_032323,3,5), "avg_brian_nova_032323_aft");
  brian_040623_bef_avg = average_hist(input_hists, slicing(brian_nova_040623,0,2), "avg_brian_nova_040623_bef");
  brian_040623_aft_avg = average_hist(input_hists, slicing(brian_nova_040623,3,5), "avg_brian_nova_040623_aft");
  brian_041323_bef_avg = average_hist(input_hists, slicing(brian_nova_041323,0,2), "avg_brian_nova_041323_bef");
  brian_041323_aft_avg = average_hist(input_hists, slicing(brian_nova_041323,3,5), "avg_brian_nova_041323_aft");
  brian_042123_bef_avg = average_hist(input_hists, slicing(brian_nova_042123,0,2), "avg_brian_nova_042123_bef");
  brian_042123_aft_avg = average_hist(input_hists, slicing(brian_nova_042123,3,5), "avg_brian_nova_042123_aft");
  brian_052423_bef_avg = average_hist(input_hists, slicing(brian_nova_052423,0,2), "avg_brian_nova_052423_bef");
  brian_052423_aft_avg = average_hist(input_hists, slicing(brian_nova_052423,3,5), "avg_brian_nova_052423_aft");
  brian_060823_bef_avg = average_hist(input_hists, slicing(brian_nova_060823,0,2), "avg_brian_nova_060823_bef");
  brian_060823_aft_avg = average_hist(input_hists, slicing(brian_nova_060823,3,5), "avg_brian_nova_060823_aft");
  brian_062323_bef_avg = average_hist(input_hists, slicing(brian_nova_062323,0,2), "avg_brian_nova_062323_bef");
  brian_062323_aft_avg = average_hist(input_hists, slicing(brian_nova_062323,3,5), "avg_brian_nova_062323_aft");

  mackenzie_920_gore_avg      = average_hist(input_hists, mackenzie_goreDRP_920, "avg_mackenzie_goreDRP_920");
  mackenzie_1018_gore_avg     = average_hist(input_hists, mackenzie_goreDRP_1018, "avg_mackenzie_goreDRP_1018");
  mackenzie_1025_gore_bef_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1025,0,2), "avg_mackenzie_goreDRP_1025_bef");
  mackenzie_1025_gore_aft_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1025,3,5), "avg_mackenzie_goreDRP_1025_aft");
  mackenzie_1101_gore_bef_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1101,0,2), "avg_mackenzie_goreDRP_1101_bef");
  mackenzie_1101_gore_aft_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1101,3,5), "avg_mackenzie_goreDRP_1101_aft");
  mackenzie_1108_gore_avg     = average_hist(input_hists, mackenzie_goreDRP_1108, "avg_mackenzie_goreDRP_1108");
  mackenzie_1115_gore_bef_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1115,0,2), "avg_mackenzie_goreDRP_1115_bef");
  mackenzie_1115_gore_aft_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1115,3,5), "avg_mackenzie_goreDRP_1115_aft");
  mackenzie_1122_gore_bef_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1122,0,2), "avg_mackenzie_goreDRP_1122_bef");
  mackenzie_1122_gore_aft_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1122,3,5), "avg_mackenzie_goreDRP_1122_aft");
  mackenzie_1129_gore_bef_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1129,0,2), "avg_mackenzie_goreDRP_1129_bef");
  mackenzie_1129_gore_aft_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1129,3,5), "avg_mackenzie_goreDRP_1129_aft");
  mackenzie_1206_gore_bef_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1206,0,2), "avg_mackenzie_goreDRP_1206_bef");
  mackenzie_1206_gore_aft_avg = average_hist(input_hists, slicing(mackenzie_goreDRP_1206,3,5), "avg_mackenzie_goreDRP_1206_aft");
  brian_020123_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_020123,0,2), "avg_brian_goreDRP_020123_bef");
  brian_020123_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_020123,3,5), "avg_brian_goreDRP_020123_aft");
  brian_020923_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_020923,0,2), "avg_brian_goreDRP_020923_bef");
  brian_020923_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_020923,3,5), "avg_brian_goreDRP_020923_aft");
  brian_021723_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_021723,0,2), "avg_brian_goreDRP_021723_bef");
  brian_021723_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_021723,3,5), "avg_brian_goreDRP_021723_aft");
  brian_030923_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_030923,0,2), "avg_brian_goreDRP_030923_bef");
  brian_030923_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_030923,3,5), "avg_brian_goreDRP_030923_aft");
  brian_031623_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_031623,0,2), "avg_brian_goreDRP_031623_bef");
  brian_031623_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_031623,3,5), "avg_brian_goreDRP_031623_aft");
  brian_032323_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_032323,0,2), "avg_brian_goreDRP_032323_bef");
  brian_032323_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_032323,3,5), "avg_brian_goreDRP_032323_aft");
  brian_040623_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_040623,0,2), "avg_brian_goreDRP_040623_bef");
  brian_040623_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_040623,3,5), "avg_brian_goreDRP_040623_aft");
  brian_041323_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_041323,0,2), "avg_brian_goreDRP_041323_bef");
  brian_041323_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_041323,3,5), "avg_brian_goreDRP_041323_aft");
  brian_042123_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_042123,0,2), "avg_brian_goreDRP_042123_bef");
  brian_042123_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_042123,3,5), "avg_brian_goreDRP_042123_aft");
  brian_052423_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_052423,0,2), "avg_brian_goreDRP_052423_bef");
  brian_052423_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_052423,3,5), "avg_brian_goreDRP_052423_aft");
  brian_060823_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_060823,0,2), "avg_brian_goreDRP_060823_bef");
  brian_060823_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_060823,3,5), "avg_brian_goreDRP_060823_aft");
  brian_062323_gore_bef_avg = average_hist(input_hists, slicing(brian_goreDRP_062323,0,2), "avg_brian_goreDRP_062323_bef");
  brian_062323_gore_aft_avg = average_hist(input_hists, slicing(brian_goreDRP_062323,3,5), "avg_brian_goreDRP_062323_aft");
  

  // difference histogram function

  // Make sure ROOT knows we want to write to the output file
  ofile->cd();

  // write the histograms to file
  brian_92316_avg->Write();
  brian_111516_avg->Write();
  brian_11618_avg->Write();
  brian_40318_avg->Write();
  brian_110420_avg->Write();
  mackenzie_920_avg->Write();
  mackenzie_1011_bef_avg->Write();
  mackenzie_1011_aft_avg->Write();
  mackenzie_1018_bef_avg->Write();
  mackenzie_1018_aft_avg->Write();
  mackenzie_1025_bef_avg->Write();
  mackenzie_1025_aft_avg->Write();
  mackenzie_1101_bef_avg->Write();
  mackenzie_1101_aft_avg->Write();
  mackenzie_1103_bef_avg->Write();
  mackenzie_1103_aft_avg->Write();
  mackenzie_1108_bef_avg->Write();
  mackenzie_1108_aft_avg->Write();
  mackenzie_1115_bef_avg->Write();
  mackenzie_1115_aft_avg->Write();
  mackenzie_1122_bef_avg->Write();
  mackenzie_1122_aft_avg->Write();
  mackenzie_1129_bef_avg->Write();
  mackenzie_1129_aft_avg->Write();
  mackenzie_1206_bef_avg->Write();
  mackenzie_1206_aft_avg->Write();
  brian_020123_bef_avg->Write();
  brian_020123_aft_avg->Write();
  brian_020923_bef_avg->Write();
  brian_020923_aft_avg->Write();
  brian_021723_bef_avg->Write();
  brian_021723_aft_avg->Write();
  brian_030923_bef_avg->Write();
  brian_030923_aft_avg->Write();
  brian_031623_bef_avg->Write();
  brian_031623_aft_avg->Write();
  brian_032323_bef_avg->Write();
  brian_032323_aft_avg->Write();
  brian_040623_bef_avg->Write();
  brian_040623_aft_avg->Write();
  brian_041323_bef_avg->Write();
  brian_041323_aft_avg->Write();
  brian_042123_bef_avg->Write();
  brian_042123_aft_avg->Write();
  brian_052423_bef_avg->Write();
  brian_052423_aft_avg->Write();
  brian_060823_bef_avg->Write();
  brian_060823_aft_avg->Write();
  brian_062323_bef_avg->Write();
  brian_062323_aft_avg->Write();

  mackenzie_920_gore_avg->Write();
  mackenzie_1018_gore_avg->Write();
  mackenzie_1025_gore_bef_avg->Write();
  mackenzie_1025_gore_aft_avg->Write();
  mackenzie_1101_gore_bef_avg->Write();
  mackenzie_1101_gore_aft_avg->Write();
  mackenzie_1108_gore_avg->Write();
  mackenzie_1115_gore_bef_avg->Write();
  mackenzie_1115_gore_aft_avg->Write();
  mackenzie_1122_gore_bef_avg->Write();
  mackenzie_1122_gore_aft_avg->Write();
  mackenzie_1129_gore_bef_avg->Write();
  mackenzie_1129_gore_aft_avg->Write();
  mackenzie_1206_gore_bef_avg->Write();
  mackenzie_1206_gore_aft_avg->Write();
  brian_020123_gore_bef_avg->Write();
  brian_020123_gore_aft_avg->Write();
  brian_020923_gore_bef_avg->Write();
  brian_020923_gore_aft_avg->Write();
  brian_021723_gore_bef_avg->Write();
  brian_021723_gore_aft_avg->Write();
  brian_030923_gore_bef_avg->Write();
  brian_030923_gore_aft_avg->Write();
  brian_031623_gore_bef_avg->Write();
  brian_031623_gore_aft_avg->Write();
  brian_032323_gore_bef_avg->Write();
  brian_032323_gore_aft_avg->Write();
  brian_040623_gore_bef_avg->Write();
  brian_040623_gore_aft_avg->Write();
  brian_041323_gore_bef_avg->Write();
  brian_041323_gore_aft_avg->Write();
  brian_042123_gore_bef_avg->Write();
  brian_042123_gore_aft_avg->Write();
  brian_052423_gore_bef_avg->Write();
  brian_052423_gore_aft_avg->Write();
  brian_060823_gore_bef_avg->Write();
  brian_060823_gore_aft_avg->Write();
  brian_062323_gore_bef_avg->Write();
  brian_062323_gore_aft_avg->Write();


  // make collections of data at specific wavelengths for aging plots
  std::vector<Float_t> aging_380;
  std::vector<Float_t> aging_390;
  std::vector<Float_t> aging_400;
  std::vector<Float_t> aging_410;
  std::vector<Float_t> aging_420;
  std::vector<Float_t> aging_430;
  std::vector<Float_t> aging_440;
  std::vector<Float_t> aging_450;
  std::vector<Float_t> aging_600;
  std::vector<Float_t> aging_600_bef;
  std::vector<Float_t> aging_600_aft;

  std::vector<Float_t> aging_380_gore;
  std::vector<Float_t> aging_390_gore;
  std::vector<Float_t> aging_400_gore;
  std::vector<Float_t> aging_410_gore;
  std::vector<Float_t> aging_420_gore;
  std::vector<Float_t> aging_430_gore;
  std::vector<Float_t> aging_440_gore;
  std::vector<Float_t> aging_450_gore;
  std::vector<Float_t> aging_600_gore;

  // make the x-axis points array, time fraction of years
  const int npoints = 41;
  Float_t aging_timefrac_xbins[npoints] = {0.0, 0.019, 0.058, 0.077, 0.096, 0.115, 0.134, 0.153, 0.192, 0.211, 0.249, 0.326, 0.364, 0.403, 0.441, 0.479, 0.493, 0.512, 0.532, 0.551, 0.570, 0.589, 0.608, 0.627, 0.633, 0.647, 0.666, 0.685, 0.704, 0.723, 0.879, 0.901, 0.923, 0.978, 0.997, 1.016, 1.055, 1.074, 1.096, 1.227, 1.268};

  // for single histogram points early in time...
  for(int ind : ordered_aging_indices){
    // add wavelength = 380nm content
    aging_380.push_back(input_hists[ind]->GetBinContent(3));
    std::cout << "printing 380nm points: " << input_hists[ind]->GetBinContent(3) << std::endl;

    // add wavelength = 390nm content
    aging_390.push_back(input_hists[ind]->GetBinContent(4));
    //std::cout << "printing point at 390nm: " << input_hists[ind]->GetBinContent(4) << std::endl;

    // add wavelength = 400nm content
    aging_400.push_back(input_hists[ind]->GetBinContent(5));
    //std::cout << "printing 400nm points: " << input_hists[ind]->GetBinContent(5) << std::endl;

    // add wavelength = 410nm content
    aging_410.push_back(input_hists[ind]->GetBinContent(6));
    //std::cout << "printing 410nm points: " << input_hists[ind]->GetBinContent(6) << std::endl;

    // add wavelength = 420nm content
    aging_420.push_back(input_hists[ind]->GetBinContent(7));
    //std::cout << "printing 420nm points: " << input_hists[ind]->GetBinContent(7) << std::endl;
    
    // add wavelength = 430nm content
    aging_430.push_back(input_hists[ind]->GetBinContent(8));
    //std::cout << "printing 430nm points: " << input_hists[ind]->GetBinContent(8) << std::endl;

    // add wavelength = 440nm content
    aging_440.push_back(input_hists[ind]->GetBinContent(9));
    //std::cout << "printing 440nm points: " << input_hists[ind]->GetBinContent(9) << std::endl;

    // add wavelength = 450nm content
    aging_450.push_back(input_hists[ind]->GetBinContent(10));
    //std::cout << "printing 450nm points: " << input_hists[ind]->GetBinContent(10) << std::endl;

    // add wavelength = 600nm content
    aging_600.push_back(input_hists[ind]->GetBinContent(25));
  }

  //add content from average hists to end of lists, nova
  aging_380.push_back(mackenzie_1011_bef_avg->GetBinContent(3));
  aging_380.push_back(mackenzie_1018_bef_avg->GetBinContent(3));
  aging_380.push_back(mackenzie_1025_bef_avg->GetBinContent(3));
  aging_380.push_back(mackenzie_1101_bef_avg->GetBinContent(3));
  aging_380.push_back(mackenzie_1103_bef_avg->GetBinContent(3));
  aging_380.push_back(mackenzie_1108_bef_avg->GetBinContent(3));
  aging_380.push_back(mackenzie_1115_bef_avg->GetBinContent(3));
  aging_380.push_back(mackenzie_1122_bef_avg->GetBinContent(3));
  aging_380.push_back(mackenzie_1129_bef_avg->GetBinContent(3));
  aging_380.push_back(mackenzie_1206_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_020123_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_020923_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_021723_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_030923_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_031623_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_032323_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_040623_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_041323_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_042123_bef_avg->GetBinContent(3));
  //aging_380.push_back(brian_052423_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_060823_bef_avg->GetBinContent(3));
  aging_380.push_back(brian_062323_bef_avg->GetBinContent(3));

  aging_390.push_back(mackenzie_1011_bef_avg->GetBinContent(4));
  aging_390.push_back(mackenzie_1018_bef_avg->GetBinContent(4));
  aging_390.push_back(mackenzie_1025_bef_avg->GetBinContent(4));
  aging_390.push_back(mackenzie_1101_bef_avg->GetBinContent(4));
  aging_390.push_back(mackenzie_1103_bef_avg->GetBinContent(4));
  aging_390.push_back(mackenzie_1108_bef_avg->GetBinContent(4));
  aging_390.push_back(mackenzie_1115_bef_avg->GetBinContent(4));
  aging_390.push_back(mackenzie_1122_bef_avg->GetBinContent(4));
  aging_390.push_back(mackenzie_1129_bef_avg->GetBinContent(4));
  aging_390.push_back(mackenzie_1206_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_020123_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_020923_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_021723_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_030923_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_031623_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_032323_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_040623_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_041323_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_042123_bef_avg->GetBinContent(4));
  //aging_390.push_back(brian_052423_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_060823_bef_avg->GetBinContent(4));
  aging_390.push_back(brian_062323_bef_avg->GetBinContent(4));
  /*std::cout << "printing 390nm vector after all push backs: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_390[i] << ", ";
  }
  std::cout << "\n";*/

  aging_400.push_back(mackenzie_1011_bef_avg->GetBinContent(5));
  aging_400.push_back(mackenzie_1018_bef_avg->GetBinContent(5));
  aging_400.push_back(mackenzie_1025_bef_avg->GetBinContent(5));
  aging_400.push_back(mackenzie_1101_bef_avg->GetBinContent(5));
  aging_400.push_back(mackenzie_1103_bef_avg->GetBinContent(5));
  aging_400.push_back(mackenzie_1108_bef_avg->GetBinContent(5));
  aging_400.push_back(mackenzie_1115_bef_avg->GetBinContent(5));
  aging_400.push_back(mackenzie_1122_bef_avg->GetBinContent(5));
  aging_400.push_back(mackenzie_1129_bef_avg->GetBinContent(5));
  aging_400.push_back(mackenzie_1206_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_020123_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_020923_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_021723_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_030923_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_031623_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_032323_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_040623_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_041323_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_042123_bef_avg->GetBinContent(5));
  //aging_400.push_back(brian_052423_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_060823_bef_avg->GetBinContent(5));
  aging_400.push_back(brian_062323_bef_avg->GetBinContent(5));
  /*std::cout << "printing 400nm vector after all push backs: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_400[i] << ", ";
  }
  std::cout << "\n";*/

  aging_410.push_back(mackenzie_1011_bef_avg->GetBinContent(6));
  aging_410.push_back(mackenzie_1018_bef_avg->GetBinContent(6));
  aging_410.push_back(mackenzie_1025_bef_avg->GetBinContent(6));
  aging_410.push_back(mackenzie_1101_bef_avg->GetBinContent(6));
  aging_410.push_back(mackenzie_1103_bef_avg->GetBinContent(6));
  aging_410.push_back(mackenzie_1108_bef_avg->GetBinContent(6));
  aging_410.push_back(mackenzie_1115_bef_avg->GetBinContent(6));
  aging_410.push_back(mackenzie_1122_bef_avg->GetBinContent(6));
  aging_410.push_back(mackenzie_1129_bef_avg->GetBinContent(6));
  aging_410.push_back(mackenzie_1206_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_020123_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_020923_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_021723_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_030923_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_031623_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_032323_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_040623_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_041323_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_042123_bef_avg->GetBinContent(6));
  //aging_410.push_back(brian_052423_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_060823_bef_avg->GetBinContent(6));
  aging_410.push_back(brian_062323_bef_avg->GetBinContent(6));
  /*std::cout << "printing 410nm vector after all push backs: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_410[i] << ", ";
  }
  std::cout << "\n";*/

  aging_420.push_back(mackenzie_1011_bef_avg->GetBinContent(7));
  aging_420.push_back(mackenzie_1018_bef_avg->GetBinContent(7));
  aging_420.push_back(mackenzie_1025_bef_avg->GetBinContent(7));
  aging_420.push_back(mackenzie_1101_bef_avg->GetBinContent(7));
  aging_420.push_back(mackenzie_1103_bef_avg->GetBinContent(7));
  aging_420.push_back(mackenzie_1108_bef_avg->GetBinContent(7));
  aging_420.push_back(mackenzie_1115_bef_avg->GetBinContent(7));
  aging_420.push_back(mackenzie_1122_bef_avg->GetBinContent(7));
  aging_420.push_back(mackenzie_1129_bef_avg->GetBinContent(7));
  aging_420.push_back(mackenzie_1206_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_020123_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_020923_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_021723_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_030923_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_031623_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_032323_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_040623_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_041323_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_042123_bef_avg->GetBinContent(7));
  //aging_420.push_back(brian_052423_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_060823_bef_avg->GetBinContent(7));
  aging_420.push_back(brian_062323_bef_avg->GetBinContent(7));
  /*std::cout << "printing 420nm vector after all push backs: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_420[i] << ", ";
  }
  std::cout << "\n";*/

  aging_430.push_back(mackenzie_1011_bef_avg->GetBinContent(8));
  aging_430.push_back(mackenzie_1018_bef_avg->GetBinContent(8));
  aging_430.push_back(mackenzie_1025_bef_avg->GetBinContent(8));
  aging_430.push_back(mackenzie_1101_bef_avg->GetBinContent(8));
  aging_430.push_back(mackenzie_1103_bef_avg->GetBinContent(8));
  aging_430.push_back(mackenzie_1108_bef_avg->GetBinContent(8));
  aging_430.push_back(mackenzie_1115_bef_avg->GetBinContent(8));
  aging_430.push_back(mackenzie_1122_bef_avg->GetBinContent(8));
  aging_430.push_back(mackenzie_1129_bef_avg->GetBinContent(8));
  aging_430.push_back(mackenzie_1206_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_020123_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_020923_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_021723_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_030923_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_031623_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_032323_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_040623_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_041323_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_042123_bef_avg->GetBinContent(8));
  //aging_430.push_back(brian_052423_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_060823_bef_avg->GetBinContent(8));
  aging_430.push_back(brian_062323_bef_avg->GetBinContent(8));
  /*std::cout << "printing 430nm vector after all push backs: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_430[i] << ", ";
  }
  std::cout << "\n";*/

  aging_440.push_back(mackenzie_1011_bef_avg->GetBinContent(9));
  aging_440.push_back(mackenzie_1018_bef_avg->GetBinContent(9));
  aging_440.push_back(mackenzie_1025_bef_avg->GetBinContent(9));
  aging_440.push_back(mackenzie_1101_bef_avg->GetBinContent(9));
  aging_440.push_back(mackenzie_1103_bef_avg->GetBinContent(9));
  aging_440.push_back(mackenzie_1108_bef_avg->GetBinContent(9));
  aging_440.push_back(mackenzie_1115_bef_avg->GetBinContent(9));
  aging_440.push_back(mackenzie_1122_bef_avg->GetBinContent(9));
  aging_440.push_back(mackenzie_1129_bef_avg->GetBinContent(9));
  aging_440.push_back(mackenzie_1206_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_020123_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_020923_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_021723_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_030923_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_031623_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_032323_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_040623_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_041323_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_042123_bef_avg->GetBinContent(9));
  //aging_440.push_back(brian_052423_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_060823_bef_avg->GetBinContent(9));
  aging_440.push_back(brian_062323_bef_avg->GetBinContent(9));

  aging_450.push_back(mackenzie_1011_bef_avg->GetBinContent(10));
  aging_450.push_back(mackenzie_1018_bef_avg->GetBinContent(10));
  aging_450.push_back(mackenzie_1025_bef_avg->GetBinContent(10));
  aging_450.push_back(mackenzie_1101_bef_avg->GetBinContent(10));
  aging_450.push_back(mackenzie_1103_bef_avg->GetBinContent(10));
  aging_450.push_back(mackenzie_1108_bef_avg->GetBinContent(10));
  aging_450.push_back(mackenzie_1115_bef_avg->GetBinContent(10));
  aging_450.push_back(mackenzie_1122_bef_avg->GetBinContent(10));
  aging_450.push_back(mackenzie_1129_bef_avg->GetBinContent(10));
  aging_450.push_back(mackenzie_1206_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_020123_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_020923_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_021723_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_030923_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_031623_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_032323_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_040623_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_041323_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_042123_bef_avg->GetBinContent(10));
  //aging_450.push_back(brian_052423_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_060823_bef_avg->GetBinContent(10));
  aging_450.push_back(brian_062323_bef_avg->GetBinContent(10));

  aging_600.push_back(mackenzie_1011_bef_avg->GetBinContent(25));
  aging_600.push_back(mackenzie_1018_bef_avg->GetBinContent(25));
  aging_600.push_back(mackenzie_1025_bef_avg->GetBinContent(25));
  aging_600.push_back(mackenzie_1101_bef_avg->GetBinContent(25));
  aging_600.push_back(mackenzie_1103_bef_avg->GetBinContent(25));
  aging_600.push_back(mackenzie_1108_bef_avg->GetBinContent(25));
  aging_600.push_back(mackenzie_1115_bef_avg->GetBinContent(25));
  aging_600.push_back(mackenzie_1122_bef_avg->GetBinContent(25));
  aging_600.push_back(mackenzie_1129_bef_avg->GetBinContent(25));
  aging_600.push_back(mackenzie_1206_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_020123_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_020923_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_021723_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_030923_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_031623_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_032323_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_040623_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_041323_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_042123_bef_avg->GetBinContent(25));
  //aging_600.push_back(brian_052423_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_060823_bef_avg->GetBinContent(25));
  aging_600.push_back(brian_062323_bef_avg->GetBinContent(25));

  aging_600_bef.push_back(mackenzie_1011_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(mackenzie_1018_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(mackenzie_1025_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(mackenzie_1101_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(mackenzie_1103_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(mackenzie_1108_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(mackenzie_1115_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(mackenzie_1122_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(mackenzie_1129_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(mackenzie_1206_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_020123_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_020923_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_021723_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_030923_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_031623_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_032323_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_040623_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_041323_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_042123_bef_avg->GetBinContent(25));
  //aging_600_bef.push_back(brian_052423_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_060823_bef_avg->GetBinContent(25));
  aging_600_bef.push_back(brian_062323_bef_avg->GetBinContent(25));

  aging_600_aft.push_back(mackenzie_1011_aft_avg->GetBinContent(25));
  aging_600_aft.push_back(mackenzie_1018_aft_avg->GetBinContent(25));
  aging_600_aft.push_back(mackenzie_1025_aft_avg->GetBinContent(25));
  aging_600_aft.push_back(mackenzie_1101_aft_avg->GetBinContent(25));
  aging_600_aft.push_back(mackenzie_1103_aft_avg->GetBinContent(25));
  aging_600_aft.push_back(mackenzie_1108_aft_avg->GetBinContent(25));
  aging_600_aft.push_back(mackenzie_1115_aft_avg->GetBinContent(25));
  aging_600_aft.push_back(mackenzie_1122_aft_avg->GetBinContent(25));
  aging_600_aft.push_back(mackenzie_1129_aft_avg->GetBinContent(25));
  aging_600_aft.push_back(mackenzie_1206_aft_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_020123_bef_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_020923_bef_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_021723_bef_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_030923_bef_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_031623_bef_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_032323_bef_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_040623_bef_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_041323_bef_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_042123_bef_avg->GetBinContent(25));
  //aging_600_aft.push_back(brian_052423_bef_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_060823_bef_avg->GetBinContent(25));
  aging_600_aft.push_back(brian_062323_bef_avg->GetBinContent(25));


  // now get Gore-DRP values

  aging_380_gore.push_back(input_hists[mackenzie_goreDRP_831]->GetBinContent(3));
  aging_380_gore.push_back(input_hists[mackenzie_goreDRP_901]->GetBinContent(3));
  aging_380_gore.push_back(mackenzie_920_gore_avg->GetBinContent(3));
  aging_380_gore.push_back(input_hists[mackenzie_goreDRP_927]->GetBinContent(3));
  aging_380_gore.push_back(input_hists[mackenzie_goreDRP_1011]->GetBinContent(3));
  aging_380_gore.push_back(mackenzie_1018_gore_avg->GetBinContent(3));
  aging_380_gore.push_back(mackenzie_1025_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(mackenzie_1101_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(mackenzie_1108_gore_avg->GetBinContent(3));
  aging_380_gore.push_back(mackenzie_1115_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(mackenzie_1122_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(mackenzie_1129_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(mackenzie_1206_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(brian_020123_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(brian_020923_gore_bef_avg->GetBinContent(3)); 
  aging_380_gore.push_back(brian_021723_gore_bef_avg->GetBinContent(3)); 
  aging_380_gore.push_back(brian_030923_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(brian_031623_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(brian_032323_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(brian_040623_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(brian_041323_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(brian_042123_gore_bef_avg->GetBinContent(3));
  //aging_380_gore.push_back(brian_052423_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(brian_060823_gore_bef_avg->GetBinContent(3));
  aging_380_gore.push_back(brian_062323_gore_bef_avg->GetBinContent(3));

  aging_390_gore.push_back(input_hists[mackenzie_goreDRP_831]->GetBinContent(4));
  aging_390_gore.push_back(input_hists[mackenzie_goreDRP_901]->GetBinContent(4));
  aging_390_gore.push_back(mackenzie_920_gore_avg->GetBinContent(4));
  aging_390_gore.push_back(input_hists[mackenzie_goreDRP_927]->GetBinContent(4));
  aging_390_gore.push_back(input_hists[mackenzie_goreDRP_1011]->GetBinContent(4));
  aging_390_gore.push_back(mackenzie_1018_gore_avg->GetBinContent(4));
  aging_390_gore.push_back(mackenzie_1025_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(mackenzie_1101_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(mackenzie_1108_gore_avg->GetBinContent(4));
  aging_390_gore.push_back(mackenzie_1115_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(mackenzie_1122_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(mackenzie_1129_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(mackenzie_1206_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(brian_020123_gore_bef_avg->GetBinContent(4)); 
  aging_390_gore.push_back(brian_020923_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(brian_021723_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(brian_030923_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(brian_031623_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(brian_032323_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(brian_040623_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(brian_041323_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(brian_042123_gore_bef_avg->GetBinContent(4));
  //aging_390_gore.push_back(brian_052423_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(brian_060823_gore_bef_avg->GetBinContent(4));
  aging_390_gore.push_back(brian_062323_gore_bef_avg->GetBinContent(4));

  aging_400_gore.push_back(input_hists[mackenzie_goreDRP_831]->GetBinContent(5));
  aging_400_gore.push_back(input_hists[mackenzie_goreDRP_901]->GetBinContent(5));
  aging_400_gore.push_back(mackenzie_920_gore_avg->GetBinContent(5));
  aging_400_gore.push_back(input_hists[mackenzie_goreDRP_927]->GetBinContent(5));
  aging_400_gore.push_back(input_hists[mackenzie_goreDRP_1011]->GetBinContent(5));
  aging_400_gore.push_back(mackenzie_1018_gore_avg->GetBinContent(5));
  aging_400_gore.push_back(mackenzie_1025_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(mackenzie_1101_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(mackenzie_1108_gore_avg->GetBinContent(5));
  aging_400_gore.push_back(mackenzie_1115_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(mackenzie_1122_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(mackenzie_1129_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(mackenzie_1206_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(brian_020123_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(brian_020923_gore_bef_avg->GetBinContent(5)); 
  aging_400_gore.push_back(brian_021723_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(brian_030923_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(brian_031623_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(brian_032323_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(brian_040623_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(brian_041323_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(brian_042123_gore_bef_avg->GetBinContent(5));
  //aging_400_gore.push_back(brian_052423_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(brian_060823_gore_bef_avg->GetBinContent(5));
  aging_400_gore.push_back(brian_062323_gore_bef_avg->GetBinContent(5));

  aging_410_gore.push_back(input_hists[mackenzie_goreDRP_831]->GetBinContent(6));
  aging_410_gore.push_back(input_hists[mackenzie_goreDRP_901]->GetBinContent(6));
  aging_410_gore.push_back(mackenzie_920_gore_avg->GetBinContent(6));
  aging_410_gore.push_back(input_hists[mackenzie_goreDRP_927]->GetBinContent(6));
  aging_410_gore.push_back(input_hists[mackenzie_goreDRP_1011]->GetBinContent(6));
  aging_410_gore.push_back(mackenzie_1018_gore_avg->GetBinContent(6));
  aging_410_gore.push_back(mackenzie_1025_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(mackenzie_1101_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(mackenzie_1108_gore_avg->GetBinContent(6));
  aging_410_gore.push_back(mackenzie_1115_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(mackenzie_1122_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(mackenzie_1129_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(mackenzie_1206_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_020123_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_020923_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_021723_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_030923_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_031623_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_032323_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_040623_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_041323_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_042123_gore_bef_avg->GetBinContent(6));
  //aging_410_gore.push_back(brian_052423_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_060823_gore_bef_avg->GetBinContent(6));
  aging_410_gore.push_back(brian_062323_gore_bef_avg->GetBinContent(6));

  aging_420_gore.push_back(input_hists[mackenzie_goreDRP_831]->GetBinContent(7));
  aging_420_gore.push_back(input_hists[mackenzie_goreDRP_901]->GetBinContent(7));
  aging_420_gore.push_back(mackenzie_920_gore_avg->GetBinContent(7));
  aging_420_gore.push_back(input_hists[mackenzie_goreDRP_927]->GetBinContent(7));
  aging_420_gore.push_back(input_hists[mackenzie_goreDRP_1011]->GetBinContent(7));
  aging_420_gore.push_back(mackenzie_1018_gore_avg->GetBinContent(7));
  aging_420_gore.push_back(mackenzie_1025_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(mackenzie_1101_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(mackenzie_1108_gore_avg->GetBinContent(7));
  aging_420_gore.push_back(mackenzie_1115_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(mackenzie_1122_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(mackenzie_1129_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(mackenzie_1206_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_020123_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_020923_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_021723_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_030923_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_031623_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_032323_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_040623_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_041323_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_042123_gore_bef_avg->GetBinContent(7));
  //aging_420_gore.push_back(brian_052423_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_060823_gore_bef_avg->GetBinContent(7));
  aging_420_gore.push_back(brian_062323_gore_bef_avg->GetBinContent(7));

  aging_430_gore.push_back(input_hists[mackenzie_goreDRP_831]->GetBinContent(8));
  aging_430_gore.push_back(input_hists[mackenzie_goreDRP_901]->GetBinContent(8));
  aging_430_gore.push_back(mackenzie_920_gore_avg->GetBinContent(8));
  aging_430_gore.push_back(input_hists[mackenzie_goreDRP_927]->GetBinContent(8));
  aging_430_gore.push_back(input_hists[mackenzie_goreDRP_1011]->GetBinContent(8));
  aging_430_gore.push_back(mackenzie_1018_gore_avg->GetBinContent(8));
  aging_430_gore.push_back(mackenzie_1025_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(mackenzie_1101_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(mackenzie_1108_gore_avg->GetBinContent(8));
  aging_430_gore.push_back(mackenzie_1115_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(mackenzie_1122_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(mackenzie_1129_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(mackenzie_1206_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_020123_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_020923_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_021723_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_030923_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_031623_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_032323_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_040623_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_041323_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_042123_gore_bef_avg->GetBinContent(8));
  //aging_430_gore.push_back(brian_052423_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_060823_gore_bef_avg->GetBinContent(8));
  aging_430_gore.push_back(brian_062323_gore_bef_avg->GetBinContent(8));

  aging_440_gore.push_back(input_hists[mackenzie_goreDRP_831]->GetBinContent(9));
  aging_440_gore.push_back(input_hists[mackenzie_goreDRP_901]->GetBinContent(9));
  aging_440_gore.push_back(mackenzie_920_gore_avg->GetBinContent(9));
  aging_440_gore.push_back(input_hists[mackenzie_goreDRP_927]->GetBinContent(9));
  aging_440_gore.push_back(input_hists[mackenzie_goreDRP_1011]->GetBinContent(9));
  aging_440_gore.push_back(mackenzie_1018_gore_avg->GetBinContent(9));
  aging_440_gore.push_back(mackenzie_1025_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(mackenzie_1101_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(mackenzie_1108_gore_avg->GetBinContent(9));
  aging_440_gore.push_back(mackenzie_1115_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(mackenzie_1122_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(mackenzie_1129_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(mackenzie_1206_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_020123_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_020923_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_021723_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_030923_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_031623_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_032323_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_040623_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_041323_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_042123_gore_bef_avg->GetBinContent(9));
  //aging_440_gore.push_back(brian_052423_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_060823_gore_bef_avg->GetBinContent(9));
  aging_440_gore.push_back(brian_062323_gore_bef_avg->GetBinContent(9));

  aging_450_gore.push_back(input_hists[mackenzie_goreDRP_831]->GetBinContent(10));
  aging_450_gore.push_back(input_hists[mackenzie_goreDRP_901]->GetBinContent(10));
  aging_450_gore.push_back(mackenzie_920_gore_avg->GetBinContent(10));
  aging_450_gore.push_back(input_hists[mackenzie_goreDRP_927]->GetBinContent(10));
  aging_450_gore.push_back(input_hists[mackenzie_goreDRP_1011]->GetBinContent(10));
  aging_450_gore.push_back(mackenzie_1018_gore_avg->GetBinContent(10));
  aging_450_gore.push_back(mackenzie_1025_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(mackenzie_1101_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(mackenzie_1108_gore_avg->GetBinContent(10));
  aging_450_gore.push_back(mackenzie_1115_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(mackenzie_1122_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(mackenzie_1129_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(mackenzie_1206_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_020123_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_020923_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_021723_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_030923_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_031623_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_032323_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_040623_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_041323_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_042123_gore_bef_avg->GetBinContent(10));
  //aging_450_gore.push_back(brian_052423_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_060823_gore_bef_avg->GetBinContent(10));
  aging_450_gore.push_back(brian_062323_gore_bef_avg->GetBinContent(10));

  aging_600_gore.push_back(input_hists[mackenzie_goreDRP_831]->GetBinContent(25));
  aging_600_gore.push_back(input_hists[mackenzie_goreDRP_901]->GetBinContent(25));
  aging_600_gore.push_back(mackenzie_920_gore_avg->GetBinContent(25));
  aging_600_gore.push_back(input_hists[mackenzie_goreDRP_927]->GetBinContent(25));
  aging_600_gore.push_back(input_hists[mackenzie_goreDRP_1011]->GetBinContent(25));
  aging_600_gore.push_back(mackenzie_1018_gore_avg->GetBinContent(25));
  aging_600_gore.push_back(mackenzie_1025_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(mackenzie_1101_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(mackenzie_1108_gore_avg->GetBinContent(25));
  aging_600_gore.push_back(mackenzie_1115_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(mackenzie_1122_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(mackenzie_1129_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(mackenzie_1206_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_020123_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_020923_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_021723_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_030923_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_031623_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_032323_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_040623_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_041323_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_042123_gore_bef_avg->GetBinContent(25));
  //aging_600_gore.push_back(brian_052423_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_060823_gore_bef_avg->GetBinContent(25));
  aging_600_gore.push_back(brian_062323_gore_bef_avg->GetBinContent(25));
  

  // get errors, only from 9/20 data and put into array for TGraph
  Float_t nova_err[9];
  nova_err[0] = stdev_error(input_hists, mackenzie_nova_920, 380);
  nova_err[1] = stdev_error(input_hists, mackenzie_nova_920, 390);
  nova_err[2] = stdev_error(input_hists, mackenzie_nova_920, 400);
  nova_err[3] = stdev_error(input_hists, mackenzie_nova_920, 410);
  nova_err[4] = stdev_error(input_hists, mackenzie_nova_920, 420);
  nova_err[5] = stdev_error(input_hists, mackenzie_nova_920, 430);
  nova_err[6] = stdev_error(input_hists, mackenzie_nova_920, 440);
  nova_err[7] = stdev_error(input_hists, mackenzie_nova_920, 450);
  nova_err[8] = stdev_error(input_hists, mackenzie_nova_920, 600);

  // make TGraph stuff for Gore-DRP
  Float_t gore_err[9];
  gore_err[0] = stdev_error(input_hists, mackenzie_goreDRP_920, 380);
  gore_err[1] = stdev_error(input_hists, mackenzie_goreDRP_920, 390);
  gore_err[2] = stdev_error(input_hists, mackenzie_goreDRP_920, 400);
  gore_err[3] = stdev_error(input_hists, mackenzie_goreDRP_920, 410);
  gore_err[4] = stdev_error(input_hists, mackenzie_goreDRP_920, 420);
  gore_err[5] = stdev_error(input_hists, mackenzie_goreDRP_920, 430);
  gore_err[6] = stdev_error(input_hists, mackenzie_goreDRP_920, 440);
  gore_err[7] = stdev_error(input_hists, mackenzie_goreDRP_920, 450);
  gore_err[8] = stdev_error(input_hists, mackenzie_goreDRP_920, 600);

  // plot 600nm points to monitor stability, gather collections
  // want these plots to have ~6 points on each wavelength... with before/after clear
  // TMultiGraph (?) with one graph for each standard set

  const int nmultpoints = 21;
  Float_t befaft_time[nmultpoints] = {0.570, 0.589, 0.608, 0.627, 0.633, 0.647, 0.666, 0.685, 0.704, 0.723, 0.879, 0.901, 0.923, 0.978, 0.997, 1.016, 1.055, 1.074, 1.096, 1.227, 1.268}; // for 600 stab and later plots
  Float_t befaft_err[nmultpoints] = {0};

  const int nmultgorepoints = 18;
  Float_t befaft_gore_time[nmultgorepoints] = {0.000, 0.019, 0.058, 0.077, 0.096, 0.115, 0.271, 0.293, 0.315, 0.370, 0.389, 0.408, 0.447, 0.466, 0.488, 0.578, 0.619, 0.660};
  Float_t befaft_gore_err[nmultgorepoints] = {0};

  const int ngorepoints = 24;
  Float_t gore_time[ngorepoints] = {0.000, 0.003, 0.055, 0.074, 0.112, 0.132, 0.151, 0.170, 0.189, 0.208, 0.227, 0.247, 0.266, 0.422, 0.444, 0.466, 0.521, 0.540, 0.559, 0.597, 0.616, 0.638, 0.770, 0.811};

  // collect points in vector, then turn them into arrays, before and after different now, bin 25 = 600nm
  std::vector<Float_t> stab_600[6];
  for(int i = 0; i < 6; i++){
    stab_600[i].push_back(input_hists[mackenzie_nova_1011[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[mackenzie_nova_1018[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[mackenzie_nova_1025[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[mackenzie_nova_1101[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[mackenzie_nova_1103[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[mackenzie_nova_1108[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[mackenzie_nova_1115[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[mackenzie_nova_1122[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[mackenzie_nova_1129[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[mackenzie_nova_1206[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_020123[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_020923[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_021723[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_030923[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_031623[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_032323[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_040623[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_041323[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_042123[i]]->GetBinContent(25));
    //stab_600[i].push_back(input_hists[brian_nova_052423[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_060823[i]]->GetBinContent(25));
    stab_600[i].push_back(input_hists[brian_nova_062323[i]]->GetBinContent(25));
  }

  Float_t stab_600_a[6][nmultpoints];
  for(int i = 0; i < 6; i++){
    std::copy(stab_600[i].begin(), stab_600[i].end(), stab_600_a[i]);
  }

  TCanvas* stability_600;
  stability_600 = new TCanvas("stability_600", "stability_600",900,600);
  TGraph* stab_nova_600[6];
  for(int i = 0; i < 6; i++){
    stab_nova_600[i] = new TGraph(nmultpoints, befaft_time, stab_600_a[i]);
  }
  stability_600->cd();
  for(int i = 0; i < 3; i++) stab_nova_600[i]->SetMarkerStyle(24);
  for(int i = 3; i < 6; i++) stab_nova_600[i]->SetMarkerStyle(20);
  stab_nova_600[0]->SetTitle("Repeated Nova Measurements at 600nm");
  stab_nova_600[0]->GetXaxis()->SetTitle("Time, (Fraction of Years)");
  stab_nova_600[0]->GetYaxis()->SetTitle("R Values at 600nm");
  stab_nova_600[0]->GetYaxis()->SetRangeUser(93.,96.);
  stab_nova_600[0]->Draw("AP");
  for(int i = 1; i < 6; i++) stab_nova_600[i]->Draw("P");
  stability_600->Update();
  stability_600->Write();
  stability_600->Close();
  
  
  gStyle->SetOptFit(11111);

  // now use these to make aging plots
  TCanvas* page[8];
  TCanvas* stability[2];
  TCanvas* beforeafter[9];
  TCanvas* beforeafter_gore[9];
  TCanvas* page600[3];
  TGraphErrors* nova_380;
  TGraphErrors* nova_390;
  TGraphErrors* nova_400;
  TGraphErrors* nova_410;
  TGraphErrors* nova_420;
  TGraphErrors* nova_430;
  TGraphErrors* nova_440;
  TGraphErrors* nova_450;
  TGraph* nova_bef_380;
  TGraph* nova_aft_380;
  TGraph* nova_bef_390;
  TGraph* nova_aft_390;
  TGraph* nova_bef_400;
  TGraph* nova_aft_400;
  TGraph* nova_bef_410;
  TGraph* nova_aft_410;
  TGraph* nova_bef_420;
  TGraph* nova_aft_420;
  TGraph* nova_bef_430;
  TGraph* nova_aft_430;
  TGraph* nova_bef_440;
  TGraph* nova_aft_440;
  TGraph* nova_bef_450;
  TGraph* nova_aft_450;
  TGraph* nova_stability;
  TGraph* gore_stability;
  TGraphErrors* nova_600;
  TGraph* nova_bef_600;
  TGraph* nova_aft_600;
  TGraph* nova_600_bef;
  TGraph* nova_600_aft;

  TCanvas* gorepage[9];
  TGraphErrors* gore_380;
  TGraphErrors* gore_390;
  TGraphErrors* gore_400;
  TGraphErrors* gore_410;
  TGraphErrors* gore_420;
  TGraphErrors* gore_430;
  TGraphErrors* gore_440;
  TGraphErrors* gore_450;
  TGraphErrors* gore_600;
  TGraph* gore_bef_380;
  TGraph* gore_aft_380;
  TGraph* gore_bef_390;
  TGraph* gore_aft_390;
  TGraph* gore_bef_400;
  TGraph* gore_aft_400;
  TGraph* gore_bef_410;
  TGraph* gore_aft_410;
  TGraph* gore_bef_420;
  TGraph* gore_aft_420;
  TGraph* gore_bef_430;
  TGraph* gore_aft_430;
  TGraph* gore_bef_440;
  TGraph* gore_aft_440;
  TGraph* gore_bef_450;
  TGraph* gore_aft_450;
  TGraph* gore_bef_600;
  TGraph* gore_aft_600;

  page[0] = new TCanvas("nova_380","nova_380",900,600);
  page[1] = new TCanvas("nova_390","nova_390",900,600);
  page[2] = new TCanvas("nova_400","nova_400",900,600);
  page[3] = new TCanvas("nova_410","nova_410",900,600);
  page[4] = new TCanvas("nova_420","nova_420",900,600);
  page[5] = new TCanvas("nova_430","nova_430",900,600);
  page[6] = new TCanvas("nova_440","nova_440",900,600);
  page[7] = new TCanvas("nova_450","nova_450",900,600);
  page600[0] = new TCanvas("nova_600","nova_600",900,600);
  page600[1] = new TCanvas("nova_600_bef","nova_600_bef",900,600);
  page600[2] = new TCanvas("nova_600_aft","nova_600_aft",900,600);
  stability[0] = new TCanvas("nova_stability", "nova_stability",900,600);
  stability[1] = new TCanvas("gore_stability", "gore_stability",900,600);

  beforeafter[0] = new TCanvas("380_beforeafter", "380_beforeafter",900,600);
  beforeafter[1] = new TCanvas("390_beforeafter", "390_beforeafter",900,600);
  beforeafter[2] = new TCanvas("400_beforeafter", "400_beforeafter",900,600);
  beforeafter[3] = new TCanvas("410_beforeafter", "410_beforeafter",900,600);
  beforeafter[4] = new TCanvas("420_beforeafter", "420_beforeafter",900,600);
  beforeafter[5] = new TCanvas("430_beforeafter", "430_beforeafter",900,600);
  beforeafter[6] = new TCanvas("440_beforeafter", "440_beforeafter",900,600);
  beforeafter[7] = new TCanvas("450_beforeafter", "450_beforeafter",900,600);
  beforeafter[8] = new TCanvas("600_beforeafter", "600_beforeafter",900,600);

  beforeafter_gore[0] = new TCanvas("380_beforeafter_gore", "380_beforeafter_gore",900,600);
  beforeafter_gore[1] = new TCanvas("390_beforeafter_gore", "390_beforeafter_gore",900,600);
  beforeafter_gore[2] = new TCanvas("400_beforeafter_gore", "400_beforeafter_gore",900,600);
  beforeafter_gore[3] = new TCanvas("410_beforeafter_gore", "410_beforeafter_gore",900,600);
  beforeafter_gore[4] = new TCanvas("420_beforeafter_gore", "420_beforeafter_gore",900,600);
  beforeafter_gore[5] = new TCanvas("430_beforeafter_gore", "430_beforeafter_gore",900,600);
  beforeafter_gore[6] = new TCanvas("440_beforeafter_gore", "440_beforeafter_gore",900,600);
  beforeafter_gore[7] = new TCanvas("450_beforeafter_gore", "450_beforeafter_gore",900,600);
  beforeafter_gore[8] = new TCanvas("600_beforeafter_gore", "600_beforeafter_gore",900,600);

  gorepage[0] = new TCanvas("gore_380","gore_380",900,600);
  gorepage[1] = new TCanvas("gore_390","gore_390",900,600);
  gorepage[2] = new TCanvas("gore_400","gore_400",900,600);
  gorepage[3] = new TCanvas("gore_410","gore_410",900,600);
  gorepage[4] = new TCanvas("gore_420","gore_420",900,600);
  gorepage[5] = new TCanvas("gore_430","gore_430",900,600);
  gorepage[6] = new TCanvas("gore_440","gore_440",900,600);
  gorepage[7] = new TCanvas("gore_450","gore_450",900,600);
  gorepage[8] = new TCanvas("gore_600","gore_600",900,600);

  gStyle->SetOptFit(11111);

  // need arrays for TGraphs
  // first 18 points all get 920 error, rest have individual errors

  // wavelength = 380nm
  Float_t empty_x_err[npoints];
  std::fill(empty_x_err, empty_x_err+npoints, 0);

  Float_t aging_380_a[npoints];
  std::copy(aging_380.begin(), aging_380.end(), aging_380_a);
  Float_t aging_380_err[npoints];
  Float_t nova_380_bef_err[nmultpoints];
  Float_t nova_380_aft_err[nmultpoints];
  std::fill(aging_380_err, aging_380_err+20, nova_err[0]);
  aging_380_err[20] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 380);
  aging_380_err[21] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 380);
  aging_380_err[22] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 380);
  aging_380_err[23] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 380);
  aging_380_err[24] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 380);
  aging_380_err[25] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 380);
  aging_380_err[26] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 380);
  aging_380_err[27] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 380);
  aging_380_err[28] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 380);
  aging_380_err[29] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 380);
  aging_380_err[30] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 380);
  aging_380_err[31] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 380);
  aging_380_err[32] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 380);
  aging_380_err[33] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 380);
  aging_380_err[34] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 380);
  aging_380_err[35] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 380);
  aging_380_err[36] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 380);
  aging_380_err[37] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 380);
  aging_380_err[38] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 380);
  //aging_380_err[39] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 380);
  aging_380_err[39] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 380);
  aging_380_err[40] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 380);
  nova_380_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 380);
  nova_380_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 380);
  nova_380_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 380);
  nova_380_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 380);
  nova_380_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 380);
  nova_380_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 380);
  nova_380_bef_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 380);
  nova_380_bef_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 380);
  nova_380_bef_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 380);
  nova_380_bef_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 380);
  nova_380_bef_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 380);
  nova_380_bef_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 380);
  nova_380_bef_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 380);
  nova_380_bef_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 380);
  nova_380_bef_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 380);
  nova_380_bef_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 380);
  nova_380_bef_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 380);
  nova_380_bef_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 380);
  nova_380_bef_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 380);
  //nova_380_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 380);
  nova_380_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 380);
  nova_380_bef_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 380);
  nova_380_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,3,5), 380);
  nova_380_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,3,5), 380);
  nova_380_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,3,5), 380);
  nova_380_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,3,5), 380);
  nova_380_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,3,5), 380);
  nova_380_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,3,5), 380);
  nova_380_aft_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,3,5), 380);
  nova_380_aft_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,3,5), 380);
  nova_380_aft_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,3,5), 380);
  nova_380_aft_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,3,5), 380);
  nova_380_aft_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,3,5), 380);
  nova_380_aft_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,3,5), 380);
  nova_380_aft_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,3,5), 380);
  nova_380_aft_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,3,5), 380);
  nova_380_aft_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,3,5), 380);
  nova_380_aft_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,3,5), 380);
  nova_380_aft_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,3,5), 380);
  nova_380_aft_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,3,5), 380);
  nova_380_aft_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,3,5), 380);
  //nova_380_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,3,5), 380);
  nova_380_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,3,5), 380);
  nova_380_aft_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,3,5), 380);

  // wavelength = 390nm
  Float_t aging_390_a[npoints];
  std::copy(aging_390.begin(), aging_390.end(), aging_390_a);
  /*std::cout << "printing 390nm array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_390_a[i] << ", ";
  }
  std::cout << "\n";*/
  Float_t aging_390_err[npoints];
  Float_t nova_390_bef_err[nmultpoints];
  Float_t nova_390_aft_err[nmultpoints];
  std::fill(aging_390_err, aging_390_err+20, nova_err[1]);
  aging_390_err[20] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 390);
  aging_390_err[21] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 390);
  aging_390_err[22] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 390);
  aging_390_err[23] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 390);
  aging_390_err[24] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 390);
  aging_390_err[25] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 390);
  aging_390_err[26] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 390);
  aging_390_err[27] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 390);
  aging_390_err[28] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 390);
  aging_390_err[29] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 390);
  aging_390_err[30] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 390);
  aging_390_err[31] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 390);
  aging_390_err[32] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 390);
  aging_390_err[33] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 390);
  aging_390_err[34] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 390);
  aging_390_err[35] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 390);
  aging_390_err[36] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 390);
  aging_390_err[37] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 390);
  aging_390_err[38] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 390);
  //aging_390_err[39] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 390);
  aging_390_err[39] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 390);
  aging_390_err[40] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 390);
  nova_390_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 390);
  nova_390_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 390);
  nova_390_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 390);
  nova_390_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 390);
  nova_390_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 390);
  nova_390_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 390);
  nova_390_bef_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 390);
  nova_390_bef_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 390);
  nova_390_bef_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 390);
  nova_390_bef_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 390);
  nova_390_bef_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 390);
  nova_390_bef_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 390);
  nova_390_bef_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 390);
  nova_390_bef_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 390);  
  nova_390_bef_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 390);
  nova_390_bef_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 390);
  nova_390_bef_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 390);
  nova_390_bef_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 390);
  nova_390_bef_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 390);
  //nova_390_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 390);
  nova_390_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 390);
  nova_390_bef_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 390);
  nova_390_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,3,5), 390);
  nova_390_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,3,5), 390);
  nova_390_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,3,5), 390);
  nova_390_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,3,5), 390);
  nova_390_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,3,5), 390);
  nova_390_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,3,5), 390);
  nova_390_aft_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,3,5), 390);
  nova_390_aft_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,3,5), 390);
  nova_390_aft_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,3,5), 390);
  nova_390_aft_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,3,5), 390);
  nova_390_aft_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,3,5), 390);
  nova_390_aft_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,3,5), 390);
  nova_390_aft_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,3,5), 390);
  nova_390_aft_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,3,5), 390);
  nova_390_aft_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,3,5), 390);
  nova_390_aft_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,3,5), 390);
  nova_390_aft_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,3,5), 390);
  nova_390_aft_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,3,5), 390);
  nova_390_aft_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,3,5), 390);
  //nova_390_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,3,5), 390);
  nova_390_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,3,5), 390);
  nova_390_aft_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,3,5), 390);

  // wavelength = 400nm
  Float_t aging_400_a[npoints];
  std::copy(aging_400.begin(), aging_400.end(), aging_400_a);
  Float_t aging_400_err[npoints];
  Float_t nova_400_bef_err[nmultpoints];
  Float_t nova_400_aft_err[nmultpoints];
  std::fill(aging_400_err, aging_400_err+20, nova_err[2]);
  aging_400_err[20] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 400);
  aging_400_err[21] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 400);
  aging_400_err[22] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 400);
  aging_400_err[23] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 400);
  aging_400_err[24] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 400);
  aging_400_err[25] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 400);
  aging_400_err[26] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 400);
  aging_400_err[27] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 400);
  aging_400_err[28] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 400);
  aging_400_err[29] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 400);
  aging_400_err[30] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 400);
  aging_400_err[31] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 400);
  aging_400_err[32] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 400);
  aging_400_err[33] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 400);
  aging_400_err[34] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 400);
  aging_400_err[35] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 400);
  aging_400_err[36] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 400);  
  aging_400_err[37] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 400);
  aging_400_err[38] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 400);
  //aging_400_err[39] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 400);
  aging_400_err[39] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 400);
  aging_400_err[40] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 400);
  nova_400_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 400);
  nova_400_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 400);
  nova_400_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 400);
  nova_400_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 400);
  nova_400_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 400);
  nova_400_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 400);
  nova_400_bef_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 400);
  nova_400_bef_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 400);
  nova_400_bef_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 400);
  nova_400_bef_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 400);
  nova_400_bef_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 400);
  nova_400_bef_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 400);
  nova_400_bef_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 400);
  nova_400_bef_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 400);
  nova_400_bef_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 400); 
  nova_400_bef_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 400);
  nova_400_bef_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 400);
  nova_400_bef_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 400);
  nova_400_bef_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 400);
  //nova_400_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 400);
  nova_400_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 400);
  nova_400_bef_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 400);
  nova_400_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,3,5), 400);
  nova_400_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,3,5), 400);
  nova_400_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,3,5), 400);
  nova_400_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,3,5), 400);
  nova_400_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,3,5), 400);
  nova_400_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,3,5), 400);
  nova_400_aft_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,3,5), 400);
  nova_400_aft_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,3,5), 400);
  nova_400_aft_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,3,5), 400);
  nova_400_aft_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,3,5), 400);
  nova_400_aft_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,3,5), 400);
  nova_400_aft_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,3,5), 400);
  nova_400_aft_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,3,5), 400);
  nova_400_aft_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,3,5), 400);
  nova_400_aft_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,3,5), 400);
  nova_400_aft_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,3,5), 400);
  nova_400_aft_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,3,5), 400);
  nova_400_aft_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,3,5), 400);
  nova_400_aft_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,3,5), 400);
  //nova_400_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,3,5), 400);
  nova_400_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,3,5), 400);
  nova_400_aft_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,3,5), 400);

  // wavelength = 410nm
  /*std::cout << "printing 410nm vector before copy to array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_410[i] << ", ";
  }
  std::cout << "\n";*/
  Float_t aging_410_a[npoints];
  std::copy(aging_410.begin(), aging_410.end(), aging_410_a);
  /*std::cout << "printing 410nm array after copy from vector: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_410_a[i] << ", ";
  }
  std::cout << "\n";*/
  Float_t aging_410_err[npoints];
  Float_t nova_410_bef_err[nmultpoints];
  Float_t nova_410_aft_err[nmultpoints];
  std::fill(aging_410_err, aging_410_err+20, nova_err[3]);
  aging_410_err[20] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 410);
  aging_410_err[21] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 410);
  aging_410_err[22] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 410);
  aging_410_err[23] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 410);
  aging_410_err[24] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 410);
  aging_410_err[25] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 410);
  aging_410_err[26] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 410);
  aging_410_err[27] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 410);
  aging_410_err[28] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 410);
  aging_410_err[29] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 410);
  aging_410_err[30] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 410);
  aging_410_err[31] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 410);
  aging_410_err[32] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 410);
  aging_410_err[33] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 410);
  aging_410_err[34] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 410);
  aging_410_err[35] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 410);
  aging_410_err[36] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 410);
  aging_410_err[37] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 410);
  aging_410_err[38] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 410);
  //aging_410_err[39] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 410);
  aging_410_err[39] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 410);
  aging_410_err[40] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 410);
  nova_410_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 410);
  nova_410_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 410);
  nova_410_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 410);
  nova_410_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 410);
  nova_410_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 410);
  nova_410_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 410);
  nova_410_bef_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 410);
  nova_410_bef_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 410);
  nova_410_bef_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 410);
  nova_410_bef_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 410);
  nova_410_bef_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 410);
  nova_410_bef_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 410);
  nova_410_bef_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 410);
  nova_410_bef_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 410);
  nova_410_bef_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 410);
  nova_410_bef_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 410);
  nova_410_bef_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 410);
  nova_410_bef_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 410);
  nova_410_bef_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 410);
  //nova_410_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 410);
  nova_410_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 410);
  nova_410_bef_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 410);
  nova_410_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,3,5), 410);
  nova_410_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,3,5), 410);
  nova_410_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,3,5), 410);
  nova_410_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,3,5), 410);
  nova_410_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,3,5), 410);
  nova_410_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,3,5), 410);
  nova_410_aft_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,3,5), 410);
  nova_410_aft_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,3,5), 410);
  nova_410_aft_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,3,5), 410);
  nova_410_aft_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,3,5), 410);
  nova_410_aft_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,3,5), 410);
  nova_410_aft_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,3,5), 410);
  nova_410_aft_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,3,5), 410);
  nova_410_aft_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,3,5), 410);
  nova_410_aft_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,3,5), 410);
  nova_410_aft_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,3,5), 410);
  nova_410_aft_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,3,5), 410);
  nova_410_aft_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,3,5), 410);
  nova_410_aft_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,3,5), 410);
  //nova_410_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,3,5), 410);
  nova_410_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,3,5), 410);
  nova_410_aft_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,3,5), 410);

  // wavelength = 420nm
  Float_t aging_420_a[npoints];
  std::copy(aging_420.begin(), aging_420.end(), aging_420_a);
  Float_t aging_420_err[npoints];
  Float_t nova_420_bef_err[nmultpoints];
  Float_t nova_420_aft_err[nmultpoints];
  std::fill(aging_420_err, aging_420_err+20, nova_err[4]);
  aging_420_err[20] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 420);
  aging_420_err[21] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 420);
  aging_420_err[22] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 420);
  aging_420_err[23] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 420);
  aging_420_err[24] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 420);
  aging_420_err[25] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 420);
  aging_420_err[26] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 420);
  aging_420_err[27] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 420);
  aging_420_err[28] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 420);
  aging_420_err[29] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 420);
  aging_420_err[30] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 420);
  aging_420_err[31] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 420);
  aging_420_err[32] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 420);
  aging_420_err[33] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 420);
  aging_420_err[34] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 420);
  aging_420_err[35] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 420);
  aging_420_err[36] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 420);
  aging_420_err[37] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 420);
  aging_420_err[38] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 420);
  //aging_420_err[39] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 420);
  aging_420_err[39] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 420);
  aging_420_err[40] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 420);
  nova_420_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 420);
  nova_420_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 420);
  nova_420_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 420);
  nova_420_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 420);
  nova_420_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 420);
  nova_420_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 420);
  nova_420_bef_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 420);
  nova_420_bef_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 420);
  nova_420_bef_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 420);
  nova_420_bef_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 420);
  nova_420_bef_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 420);
  nova_420_bef_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 420);
  nova_420_bef_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 420);
  nova_420_bef_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 420);
  nova_420_bef_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 420);
  nova_420_bef_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 420);
  nova_420_bef_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 420);
  nova_420_bef_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 420);
  nova_420_bef_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 420);
  //nova_420_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 420);
  nova_420_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 420);
  nova_420_bef_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 420);
  nova_420_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,3,5), 420);
  nova_420_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,3,5), 420);
  nova_420_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,3,5), 420);
  nova_420_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,3,5), 420);
  nova_420_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,3,5), 420);
  nova_420_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,3,5), 420);
  nova_420_aft_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,3,5), 420);
  nova_420_aft_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,3,5), 420);
  nova_420_aft_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,3,5), 420);
  nova_420_aft_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,3,5), 420);
  nova_420_aft_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,3,5), 420);
  nova_420_aft_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,3,5), 420);
  nova_420_aft_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,3,5), 420);
  nova_420_aft_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,3,5), 420);
  nova_420_aft_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,3,5), 420);
  nova_420_aft_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,3,5), 420);
  nova_420_aft_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,3,5), 420);
  nova_420_aft_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,3,5), 420);
  nova_420_aft_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,3,5), 420);
  //nova_420_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,3,5), 420);
  nova_420_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,3,5), 420);
  nova_420_aft_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,3,5), 420);

  // wavelength = 430nm
  Float_t aging_430_a[npoints];
  std::copy(aging_430.begin(), aging_430.end(), aging_430_a);
  Float_t aging_430_err[npoints];
  Float_t nova_430_bef_err[nmultpoints];
  Float_t nova_430_aft_err[nmultpoints];
  std::fill(aging_430_err, aging_430_err+20, nova_err[5]);
  aging_430_err[20] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 430);
  aging_430_err[21] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 430);
  aging_430_err[22] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 430);
  aging_430_err[23] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 430);
  aging_430_err[24] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 430);
  aging_430_err[25] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 430);
  aging_430_err[26] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 430);
  aging_430_err[27] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 430);
  aging_430_err[28] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 430);
  aging_430_err[29] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 430);
  aging_430_err[30] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 430);
  aging_430_err[31] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 430);
  aging_430_err[32] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 430);
  aging_430_err[33] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 430);
  aging_430_err[34] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 430);
  aging_430_err[35] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 430);
  aging_430_err[36] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 430);
  aging_430_err[37] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 430);
  aging_430_err[38] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 430);
  //aging_430_err[39] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 430);
  aging_430_err[39] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 430);
  aging_430_err[40] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 430);
  nova_430_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 430);
  nova_430_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 430);
  nova_430_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 430);
  nova_430_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 430);
  nova_430_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 430);
  nova_430_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 430);
  nova_430_bef_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 430);
  nova_430_bef_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 430);
  nova_430_bef_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 430);
  nova_430_bef_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 430);
  nova_430_bef_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 430);
  nova_430_bef_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 430);
  nova_430_bef_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 430);
  nova_430_bef_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 430);
  nova_430_bef_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 430);
  nova_430_bef_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 430);
  nova_430_bef_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 430);
  nova_430_bef_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 430);
  nova_430_bef_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 430);
  //nova_430_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 430);
  nova_430_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 430);
  nova_430_bef_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 430);
  nova_430_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,3,5), 430);
  nova_430_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,3,5), 430);
  nova_430_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,3,5), 430);
  nova_430_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,3,5), 430);
  nova_430_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,3,5), 430);
  nova_430_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,3,5), 430);
  nova_430_aft_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,3,5), 430);
  nova_430_aft_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,3,5), 430);
  nova_430_aft_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,3,5), 430);
  nova_430_aft_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,3,5), 430);
  nova_430_aft_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,3,5), 430);
  nova_430_aft_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,3,5), 430);
  nova_430_aft_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,3,5), 430);
  nova_430_aft_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,3,5), 430);
  nova_430_aft_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,3,5), 430);
  nova_430_aft_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,3,5), 430);
  nova_430_aft_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,3,5), 430);
  nova_430_aft_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,3,5), 430);
  nova_430_aft_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,3,5), 430);
  //nova_430_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,3,5), 430);
  nova_430_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,3,5), 430);
  nova_430_aft_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,3,5), 430);

  // wavelength = 440nm
  Float_t aging_440_a[npoints];
  std::copy(aging_440.begin(), aging_440.end(), aging_440_a);
  Float_t aging_440_err[npoints];
  Float_t nova_440_bef_err[nmultpoints];
  Float_t nova_440_aft_err[nmultpoints];
  std::fill(aging_440_err, aging_440_err+20, nova_err[6]);
  aging_440_err[20] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 440);
  aging_440_err[21] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 440);
  aging_440_err[22] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 440);
  aging_440_err[23] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 440);
  aging_440_err[24] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 440);
  aging_440_err[25] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 440);
  aging_440_err[26] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 440);
  aging_440_err[27] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 440);
  aging_440_err[28] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 440);
  aging_440_err[29] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 440);
  aging_440_err[30] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 440);
  aging_440_err[31] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 440);
  aging_440_err[32] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 440);
  aging_440_err[33] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 440);
  aging_440_err[34] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 440);
  aging_440_err[35] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 440);
  aging_440_err[36] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 440);
  aging_440_err[37] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 440);
  aging_440_err[38] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 440);
  //aging_440_err[39] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 440);
  aging_440_err[39] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 440);
  aging_440_err[40] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 440);
  nova_440_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 440);
  nova_440_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 440);
  nova_440_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 440);
  nova_440_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 440);
  nova_440_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 440);
  nova_440_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 440);
  nova_440_bef_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 440);
  nova_440_bef_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 440);
  nova_440_bef_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 440);
  nova_440_bef_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 440);
  nova_440_bef_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 440);
  nova_440_bef_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 440);
  nova_440_bef_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 440);
  nova_440_bef_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 440);
  nova_440_bef_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 440);
  nova_440_bef_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 440);
  nova_440_bef_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 440);
  nova_440_bef_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 440);
  nova_440_bef_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 440);
  //nova_440_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 440);
  nova_440_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 440);
  nova_440_bef_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 440);
  nova_440_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,3,5), 440);
  nova_440_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,3,5), 440);
  nova_440_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,3,5), 440);
  nova_440_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,3,5), 440);
  nova_440_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,3,5), 440);
  nova_440_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,3,5), 440);
  nova_440_aft_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,3,5), 440);
  nova_440_aft_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,3,5), 440);
  nova_440_aft_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,3,5), 440);
  nova_440_aft_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,3,5), 440);
  nova_440_aft_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,3,5), 440);
  nova_440_aft_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,3,5), 440);
  nova_440_aft_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,3,5), 440);
  nova_440_aft_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,3,5), 440);
  nova_440_aft_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,3,5), 440);
  nova_440_aft_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,3,5), 440);
  nova_440_aft_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,3,5), 440);
  nova_440_aft_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,3,5), 440);
  nova_440_aft_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,3,5), 440);
  //nova_440_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,3,5), 440);
  nova_440_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,3,5), 440);
  nova_440_aft_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,3,5), 440);

  // wavelength = 450nm
  Float_t aging_450_a[npoints];
  std::copy(aging_450.begin(), aging_450.end(), aging_450_a);
  Float_t aging_450_err[npoints];
  Float_t nova_450_bef_err[nmultpoints];
  Float_t nova_450_aft_err[nmultpoints];
  std::fill(aging_450_err, aging_450_err+20, nova_err[7]);
  aging_450_err[20] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 450);
  aging_450_err[21] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 450);
  aging_450_err[22] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 450);
  aging_450_err[23] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 450);
  aging_450_err[24] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 450);
  aging_450_err[25] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 450);
  aging_450_err[26] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 450);
  aging_450_err[27] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 450);
  aging_450_err[28] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 450);
  aging_450_err[29] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 450);
  aging_450_err[30] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 450);
  aging_450_err[31] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 450);
  aging_450_err[32] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 450);
  aging_450_err[33] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 450);
  aging_450_err[34] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 450);
  aging_450_err[35] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 450);
  aging_450_err[36] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 450);
  aging_450_err[37] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 450);
  aging_450_err[38] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 450);
  //aging_450_err[39] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 450);
  aging_450_err[39] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 450);
  aging_450_err[40] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 450);
  nova_450_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 450);
  nova_450_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 450);
  nova_450_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 450);
  nova_450_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 450);
  nova_450_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 450);
  nova_450_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 450);
  nova_450_bef_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 450);
  nova_450_bef_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 450);
  nova_450_bef_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 450);
  nova_450_bef_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 450);
  nova_450_bef_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 450);
  nova_450_bef_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 450);
  nova_450_bef_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 450);
  nova_450_bef_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 450);
  nova_450_bef_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 450);
  nova_450_bef_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 450);
  nova_450_bef_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 450);
  nova_450_bef_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 450);
  nova_450_bef_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 450);
  //nova_450_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 450);
  nova_450_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 450);
  nova_450_bef_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 450);
  nova_450_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,3,5), 450);
  nova_450_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,3,5), 450);
  nova_450_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,3,5), 450);
  nova_450_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,3,5), 450);
  nova_450_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,3,5), 450);
  nova_450_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,3,5), 450);
  nova_450_aft_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,3,5), 450);
  nova_450_aft_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,3,5), 450);
  nova_450_aft_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,3,5), 450);
  nova_450_aft_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,3,5), 450);
  nova_450_aft_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,3,5), 450);
  nova_450_aft_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,3,5), 450);
  nova_450_aft_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,3,5), 450);
  nova_450_aft_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,3,5), 450);
  nova_450_aft_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,3,5), 450);
  nova_450_aft_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,3,5), 450);
  nova_450_aft_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,3,5), 450);
  nova_450_aft_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,3,5), 450);
  nova_450_aft_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,3,5), 450);
  //nova_450_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,3,5), 450);
  nova_450_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,3,5), 450);
  nova_450_aft_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,3,5), 450);

  // wavelength = 600nm
  Float_t aging_600_a[npoints];
  Float_t aging_600_bef_a[nmultpoints];
  Float_t aging_600_aft_a[nmultpoints];
  std::copy(aging_600.begin(), aging_600.end(), aging_600_a);
  std::copy(aging_600_bef.begin(), aging_600_bef.end(), aging_600_bef_a);
  std::copy(aging_600_aft.begin(), aging_600_aft.end(), aging_600_aft_a);
  Float_t aging_600_err[npoints];
  Float_t nova_600_bef_err[nmultpoints];
  Float_t nova_600_aft_err[nmultpoints];
  std::fill(aging_600_err, aging_600_err+20, nova_err[8]);
  aging_600_err[20] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 600);
  aging_600_err[21] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 600);
  aging_600_err[22] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 600);
  aging_600_err[23] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 600);
  aging_600_err[24] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 600);
  aging_600_err[25] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 600);
  aging_600_err[26] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 600);
  aging_600_err[27] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 600);
  aging_600_err[28] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 600);
  aging_600_err[29] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 600);
  aging_600_err[30] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 600);
  aging_600_err[31] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 600);
  aging_600_err[32] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 600);
  aging_600_err[33] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 600);
  aging_600_err[34] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 600);
  aging_600_err[35] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 600);
  aging_600_err[36] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 600);
  aging_600_err[37] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 600);
  aging_600_err[38] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 600);
  //aging_600_err[39] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 600);
  aging_600_err[39] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 600);
  aging_600_err[40] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 600);
  nova_600_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,0,2), 600);
  nova_600_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,0,2), 600);
  nova_600_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,0,2), 600);
  nova_600_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,0,2), 600);
  nova_600_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,0,2), 600);
  nova_600_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,0,2), 600);
  nova_600_bef_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,0,2), 600);
  nova_600_bef_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,0,2), 600);
  nova_600_bef_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,0,2), 600);
  nova_600_bef_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,0,2), 600);
  nova_600_bef_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,0,2), 600);
  nova_600_bef_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,0,2), 600);
  nova_600_bef_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,0,2), 600);
  nova_600_bef_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,0,2), 600);
  nova_600_bef_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,0,2), 600);
  nova_600_bef_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,0,2), 600);
  nova_600_bef_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,0,2), 600);
  nova_600_bef_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,0,2), 600);
  nova_600_bef_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,0,2), 600);
  //nova_600_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,0,2), 600);
  nova_600_bef_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,0,2), 600);
  nova_600_bef_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,0,2), 600);
  nova_600_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_nova_1011,3,5), 600);
  nova_600_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_nova_1018,3,5), 600);
  nova_600_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_nova_1025,3,5), 600);
  nova_600_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_nova_1101,3,5), 600);
  nova_600_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_nova_1103,3,5), 600);
  nova_600_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_nova_1108,3,5), 600);
  nova_600_aft_err[6] = stdev_error(input_hists, slicing(mackenzie_nova_1115,3,5), 600);
  nova_600_aft_err[7] = stdev_error(input_hists, slicing(mackenzie_nova_1122,3,5), 600);
  nova_600_aft_err[8] = stdev_error(input_hists, slicing(mackenzie_nova_1129,3,5), 600);
  nova_600_aft_err[9] = stdev_error(input_hists, slicing(mackenzie_nova_1206,3,5), 600);
  nova_600_aft_err[10] = stdev_error(input_hists, slicing(brian_nova_020123,3,5), 600);
  nova_600_aft_err[11] = stdev_error(input_hists, slicing(brian_nova_020923,3,5), 600);
  nova_600_aft_err[12] = stdev_error(input_hists, slicing(brian_nova_021723,3,5), 600);
  nova_600_aft_err[13] = stdev_error(input_hists, slicing(brian_nova_030923,3,5), 600);
  nova_600_aft_err[14] = stdev_error(input_hists, slicing(brian_nova_031623,3,5), 600);
  nova_600_aft_err[15] = stdev_error(input_hists, slicing(brian_nova_032323,3,5), 600);
  nova_600_aft_err[16] = stdev_error(input_hists, slicing(brian_nova_040623,3,5), 600);
  nova_600_aft_err[17] = stdev_error(input_hists, slicing(brian_nova_041323,3,5), 600);
  nova_600_aft_err[18] = stdev_error(input_hists, slicing(brian_nova_042123,3,5), 600);
  //nova_600_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_052423,3,5), 600);
  nova_600_aft_err[19] = stdev_error(input_hists, slicing(brian_nova_060823,3,5), 600);
  nova_600_aft_err[20] = stdev_error(input_hists, slicing(brian_nova_062323,3,5), 600);
  /*std::cout << "600nm errors: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_600_err[i] << ", ";
  }
  std::cout << "\n";*/

  // now goreDRP errors
  Float_t empty_x_err_gore[ngorepoints];
  std::fill(empty_x_err_gore, empty_x_err_gore+ngorepoints, 0);

  Float_t aging_380_gore_a[ngorepoints];
  std::copy(aging_380_gore.begin(), aging_380_gore.end(), aging_380_gore_a);
  Float_t aging_380_gore_err[ngorepoints];
  std::fill(aging_380_gore_err, aging_380_gore_err+5, gore_err[0]);
  aging_380_gore_err[5] = stdev_error(input_hists, mackenzie_goreDRP_1018, 380);
  std::cout << "6th error point, 380nm = "  << aging_380_gore_err[5] << std::endl;
  aging_380_gore_err[6] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 380);
  aging_380_gore_err[7] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 380);
  aging_380_gore_err[8] = stdev_error(input_hists, mackenzie_goreDRP_1108, 380);
  aging_380_gore_err[9] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 380);
  aging_380_gore_err[10] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 380);
  aging_380_gore_err[11] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 380);
  aging_380_gore_err[12] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 380);
  aging_380_gore_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 380);
  aging_380_gore_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 380);
  aging_380_gore_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 380);
  aging_380_gore_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 380);
  aging_380_gore_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 380);
  aging_380_gore_err[18] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 380);
  aging_380_gore_err[19] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 380);
  aging_380_gore_err[20] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 380);
  aging_380_gore_err[21] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 380);
  //aging_380_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 380);
  aging_380_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 380);
  aging_380_gore_err[23] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 380);
  Float_t aging_380_gore_bef_err[nmultgorepoints];
  Float_t aging_380_gore_aft_err[nmultgorepoints];
  aging_380_gore_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 380);
  aging_380_gore_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 380);
  aging_380_gore_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 380);
  aging_380_gore_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 380);
  aging_380_gore_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 380);
  aging_380_gore_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 380);
  aging_380_gore_bef_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 380);
  aging_380_gore_bef_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 380);
  aging_380_gore_bef_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 380);
  aging_380_gore_bef_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 380);
  aging_380_gore_bef_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 380);
  aging_380_gore_bef_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 380);
  aging_380_gore_bef_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 380);
  aging_380_gore_bef_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 380);
  aging_380_gore_bef_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 380);
  aging_380_gore_bef_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 380);
  aging_380_gore_bef_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 380);
  aging_380_gore_bef_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 380);
  aging_380_gore_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,3,5), 380);
  aging_380_gore_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,3,5), 380);
  aging_380_gore_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,3,5), 380);
  aging_380_gore_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,3,5), 380);
  aging_380_gore_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,3,5), 380);
  aging_380_gore_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,3,5), 380);
  aging_380_gore_aft_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,3,5), 380);
  aging_380_gore_aft_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,3,5), 380);
  aging_380_gore_aft_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,3,5), 380);
  aging_380_gore_aft_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,3,5), 380);
  aging_380_gore_aft_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,3,5), 380);
  aging_380_gore_aft_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,3,5), 380);
  aging_380_gore_aft_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,3,5), 380);
  aging_380_gore_aft_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,3,5), 380);
  aging_380_gore_aft_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,3,5), 380);
  aging_380_gore_aft_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,3,5), 380);
  aging_380_gore_aft_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,3,5), 380);
  aging_380_gore_aft_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,3,5), 380);

  Float_t aging_390_gore_a[ngorepoints];
  std::copy(aging_390_gore.begin(), aging_390_gore.end(), aging_390_gore_a);
  Float_t aging_390_gore_err[ngorepoints];
  std::fill(aging_390_gore_err, aging_390_gore_err+5, gore_err[1]);
  aging_390_gore_err[5] = stdev_error(input_hists, mackenzie_goreDRP_1018, 390);
  aging_390_gore_err[6] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 390);
  aging_390_gore_err[7] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 390);
  aging_390_gore_err[8] = stdev_error(input_hists, mackenzie_goreDRP_1108, 390);
  aging_390_gore_err[9] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 390);
  aging_390_gore_err[10] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 390);
  aging_390_gore_err[11] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 390);
  aging_390_gore_err[12] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 390);
  aging_390_gore_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 390);
  aging_390_gore_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 390);
  aging_390_gore_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 390);
  aging_390_gore_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 390);
  aging_390_gore_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 390);
  aging_390_gore_err[18] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 390);
  aging_390_gore_err[19] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 390);
  aging_390_gore_err[20] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 390);
  aging_390_gore_err[21] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 390);
  //aging_390_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 390);
  aging_390_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 390);
  aging_390_gore_err[23] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 390);
  Float_t aging_390_gore_bef_err[nmultgorepoints];
  Float_t aging_390_gore_aft_err[nmultgorepoints];
  aging_390_gore_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 390);
  aging_390_gore_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 390);
  aging_390_gore_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 390);
  aging_390_gore_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 390);
  aging_390_gore_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 390);
  aging_390_gore_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 390);
  aging_390_gore_bef_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 390);
  aging_390_gore_bef_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 390);
  aging_390_gore_bef_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 390);
  aging_390_gore_bef_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 390);
  aging_390_gore_bef_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 390);
  aging_390_gore_bef_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 390);
  aging_390_gore_bef_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 390);
  aging_390_gore_bef_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 390);
  aging_390_gore_bef_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 390);
  aging_390_gore_bef_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 390);
  aging_390_gore_bef_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 390);
  aging_390_gore_bef_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 390);
  aging_390_gore_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,3,5), 390);
  aging_390_gore_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,3,5), 390);
  aging_390_gore_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,3,5), 390);
  aging_390_gore_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,3,5), 390);
  aging_390_gore_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,3,5), 390);
  aging_390_gore_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,3,5), 390);
  aging_390_gore_aft_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,3,5), 390);
  aging_390_gore_aft_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,3,5), 390);
  aging_390_gore_aft_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,3,5), 390);
  aging_390_gore_aft_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,3,5), 390);
  aging_390_gore_aft_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,3,5), 390);
  aging_390_gore_aft_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,3,5), 390);
  aging_390_gore_aft_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,3,5), 390);
  aging_390_gore_aft_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,3,5), 390);
  aging_390_gore_aft_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,3,5), 390);
  aging_390_gore_aft_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,3,5), 390);
  aging_390_gore_aft_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,3,5), 390);
  aging_390_gore_aft_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,3,5), 390);

  Float_t aging_400_gore_a[ngorepoints];
  std::copy(aging_400_gore.begin(), aging_400_gore.end(), aging_400_gore_a);
  Float_t aging_400_gore_err[ngorepoints];
  std::fill(aging_400_gore_err, aging_400_gore_err+5, gore_err[2]);
  aging_400_gore_err[5] = stdev_error(input_hists, mackenzie_goreDRP_1018, 400);
  aging_400_gore_err[6] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 400);
  aging_400_gore_err[7] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 400);
  aging_400_gore_err[8] = stdev_error(input_hists, mackenzie_goreDRP_1108, 400);
  aging_400_gore_err[9] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 400);
  aging_400_gore_err[10] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 400);
  aging_400_gore_err[11] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 400);
  aging_400_gore_err[12] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 400);
  aging_400_gore_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 400);
  aging_400_gore_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 400);
  aging_400_gore_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 400);
  aging_400_gore_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 400);
  aging_400_gore_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 400);
  aging_400_gore_err[18] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 400);
  aging_400_gore_err[19] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 400);
  aging_400_gore_err[20] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 400);
  aging_400_gore_err[21] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 400);
  //aging_400_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 400);
  aging_400_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 400);
  aging_400_gore_err[23] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 400);
  Float_t aging_400_gore_bef_err[nmultgorepoints];
  Float_t aging_400_gore_aft_err[nmultgorepoints];
  aging_400_gore_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 400);
  aging_400_gore_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 400);
  aging_400_gore_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 400);
  aging_400_gore_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 400);
  aging_400_gore_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 400);
  aging_400_gore_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 400);
  aging_400_gore_bef_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 400);
  aging_400_gore_bef_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 400);
  aging_400_gore_bef_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 400);
  aging_400_gore_bef_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 400);
  aging_400_gore_bef_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 400);
  aging_400_gore_bef_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 400);
  aging_400_gore_bef_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 400);
  aging_400_gore_bef_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 400);
  aging_400_gore_bef_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 400);
  aging_400_gore_bef_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 400);
  aging_400_gore_bef_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 400);
  aging_400_gore_bef_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 400);
  aging_400_gore_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,3,5), 400);
  aging_400_gore_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,3,5), 400);
  aging_400_gore_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,3,5), 400);
  aging_400_gore_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,3,5), 400);
  aging_400_gore_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,3,5), 400);
  aging_400_gore_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,3,5), 400);
  aging_400_gore_aft_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,3,5), 400);
  aging_400_gore_aft_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,3,5), 400);
  aging_400_gore_aft_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,3,5), 400);
  aging_400_gore_aft_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,3,5), 400);
  aging_400_gore_aft_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,3,5), 400);
  aging_400_gore_aft_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,3,5), 400);
  aging_400_gore_aft_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,3,5), 400);
  aging_400_gore_aft_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,3,5), 400);
  aging_400_gore_aft_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,3,5), 400);
  aging_400_gore_aft_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,3,5), 400);
  aging_400_gore_aft_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,3,5), 400);
  aging_400_gore_aft_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,3,5), 400);

  Float_t aging_410_gore_a[ngorepoints];
  std::copy(aging_410_gore.begin(), aging_410_gore.end(), aging_410_gore_a);
  Float_t aging_410_gore_err[ngorepoints];
  std::fill(aging_410_gore_err, aging_410_gore_err+5, gore_err[3]);
  aging_410_gore_err[5] = stdev_error(input_hists, mackenzie_goreDRP_1018, 410);
  aging_410_gore_err[6] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 410);
  aging_410_gore_err[7] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 410);
  aging_410_gore_err[8] = stdev_error(input_hists, mackenzie_goreDRP_1108, 410);
  aging_410_gore_err[9] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 410);
  aging_410_gore_err[10] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 410);
  aging_410_gore_err[11] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 410);
  aging_410_gore_err[12] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 410);
  aging_410_gore_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 410);
  aging_410_gore_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 410);
  aging_410_gore_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 410);
  aging_410_gore_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 410);
  aging_410_gore_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 410);
  aging_410_gore_err[18] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 410);
  aging_410_gore_err[19] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 410);
  aging_410_gore_err[20] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 410);
  aging_410_gore_err[21] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 410);
  //aging_410_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 410);
  aging_410_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 410);
  aging_410_gore_err[23] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 410);
  Float_t aging_410_gore_bef_err[nmultgorepoints];
  Float_t aging_410_gore_aft_err[nmultgorepoints];
  aging_410_gore_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 410);
  aging_410_gore_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 410);
  aging_410_gore_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 410);
  aging_410_gore_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 410);
  aging_410_gore_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 410);
  aging_410_gore_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 410);
  aging_410_gore_bef_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 410);
  aging_410_gore_bef_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 410);
  aging_410_gore_bef_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 410);
  aging_410_gore_bef_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 410);
  aging_410_gore_bef_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 410);
  aging_410_gore_bef_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 410);
  aging_410_gore_bef_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 410);
  aging_410_gore_bef_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 410);
  aging_410_gore_bef_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 410);
  aging_410_gore_bef_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 410);
  aging_410_gore_bef_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 410);
  aging_410_gore_bef_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 410);
  aging_410_gore_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,3,5), 410);
  aging_410_gore_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,3,5), 410);
  aging_410_gore_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,3,5), 410);
  aging_410_gore_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,3,5), 410);
  aging_410_gore_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,3,5), 410);
  aging_410_gore_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,3,5), 410);
  aging_410_gore_aft_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,3,5), 410);
  aging_410_gore_aft_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,3,5), 410);
  aging_410_gore_aft_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,3,5), 410);
  aging_410_gore_aft_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,3,5), 410);
  aging_410_gore_aft_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,3,5), 410);
  aging_410_gore_aft_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,3,5), 410);
  aging_410_gore_aft_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,3,5), 410);
  aging_410_gore_aft_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,3,5), 410);
  aging_410_gore_aft_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,3,5), 410);
  aging_410_gore_aft_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,3,5), 410);
  aging_410_gore_aft_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,3,5), 410);
  aging_410_gore_aft_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,3,5), 410);

  Float_t aging_420_gore_a[ngorepoints];
  std::copy(aging_420_gore.begin(), aging_420_gore.end(), aging_420_gore_a);
  Float_t aging_420_gore_err[ngorepoints];
  std::fill(aging_420_gore_err, aging_420_gore_err+5, gore_err[4]);
  aging_420_gore_err[5] = stdev_error(input_hists, mackenzie_goreDRP_1018, 420);
  aging_420_gore_err[6] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 420);
  aging_420_gore_err[7] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 420);
  aging_420_gore_err[8] = stdev_error(input_hists, mackenzie_goreDRP_1108, 420);
  aging_420_gore_err[9] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 420);
  aging_420_gore_err[10] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 420);
  aging_420_gore_err[11] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 420);
  aging_420_gore_err[12] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 420);
  aging_420_gore_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 420);
  aging_420_gore_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 420);
  aging_420_gore_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 420);
  aging_420_gore_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 420);
  aging_420_gore_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 420);
  aging_420_gore_err[18] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 420);
  aging_420_gore_err[19] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 420);
  aging_420_gore_err[20] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 420);
  aging_420_gore_err[21] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 420);
  //aging_420_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 420);
  aging_420_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 420);
  aging_420_gore_err[23] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 420);
  Float_t aging_420_gore_bef_err[nmultgorepoints];
  Float_t aging_420_gore_aft_err[nmultgorepoints];
  aging_420_gore_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 420);
  aging_420_gore_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 420);
  aging_420_gore_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 420);
  aging_420_gore_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 420);
  aging_420_gore_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 420);
  aging_420_gore_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 420);
  aging_420_gore_bef_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 420);
  aging_420_gore_bef_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 420);
  aging_420_gore_bef_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 420);
  aging_420_gore_bef_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 420);
  aging_420_gore_bef_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 420);
  aging_420_gore_bef_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 420);
  aging_420_gore_bef_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 420);
  aging_420_gore_bef_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 420);
  aging_420_gore_bef_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 420);
  aging_420_gore_bef_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 420);
  aging_420_gore_bef_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 420);
  aging_420_gore_bef_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 420);
  aging_420_gore_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,3,5), 420);
  aging_420_gore_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,3,5), 420);
  aging_420_gore_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,3,5), 420);
  aging_420_gore_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,3,5), 420);
  aging_420_gore_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,3,5), 420);
  aging_420_gore_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,3,5), 420);
  aging_420_gore_aft_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,3,5), 420);
  aging_420_gore_aft_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,3,5), 420);
  aging_420_gore_aft_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,3,5), 420);
  aging_420_gore_aft_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,3,5), 420);
  aging_420_gore_aft_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,3,5), 420);
  aging_420_gore_aft_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,3,5), 420);
  aging_420_gore_aft_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,3,5), 420);
  aging_420_gore_aft_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,3,5), 420);
  aging_420_gore_aft_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,3,5), 420);
  aging_420_gore_aft_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,3,5), 420);
  aging_420_gore_aft_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,3,5), 420);
  aging_420_gore_aft_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,3,5), 420);

  Float_t aging_430_gore_a[ngorepoints];
  std::copy(aging_430_gore.begin(), aging_430_gore.end(), aging_430_gore_a);
  Float_t aging_430_gore_err[ngorepoints];
  std::fill(aging_430_gore_err, aging_430_gore_err+5, gore_err[5]);
  aging_430_gore_err[5] = stdev_error(input_hists, mackenzie_goreDRP_1018, 430);
  aging_430_gore_err[6] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 430);
  aging_430_gore_err[7] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 430);
  aging_430_gore_err[8] = stdev_error(input_hists, mackenzie_goreDRP_1108, 430);
  aging_430_gore_err[9] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 430);
  aging_430_gore_err[10] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 430);
  aging_430_gore_err[11] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 430);
  aging_430_gore_err[12] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 430);
  aging_430_gore_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 430);
  aging_430_gore_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 430);
  aging_430_gore_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 430);
  aging_430_gore_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 430);
  aging_430_gore_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 430);
  aging_430_gore_err[18] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 430);
  aging_430_gore_err[19] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 430);
  aging_430_gore_err[20] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 430);
  aging_430_gore_err[21] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 430);
  //aging_430_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 430);
  aging_430_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 430);
  aging_430_gore_err[23] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 430);
  Float_t aging_430_gore_bef_err[nmultgorepoints];
  Float_t aging_430_gore_aft_err[nmultgorepoints];
  aging_430_gore_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 430);
  aging_430_gore_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 430);
  aging_430_gore_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 430);
  aging_430_gore_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 430);
  aging_430_gore_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 430);
  aging_430_gore_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 430);
  aging_430_gore_bef_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 430);
  aging_430_gore_bef_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 430);
  aging_430_gore_bef_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 430);
  aging_430_gore_bef_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 430);
  aging_430_gore_bef_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 430);
  aging_430_gore_bef_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 430);
  aging_430_gore_bef_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 430);
  aging_430_gore_bef_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 430);
  aging_430_gore_bef_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 430);
  aging_430_gore_bef_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 430);
  aging_430_gore_bef_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 430);
  aging_430_gore_bef_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 430);
  aging_430_gore_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,3,5), 430);
  aging_430_gore_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,3,5), 430);
  aging_430_gore_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,3,5), 430);
  aging_430_gore_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,3,5), 430);
  aging_430_gore_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,3,5), 430);
  aging_430_gore_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,3,5), 430);
  aging_430_gore_aft_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,3,5), 430);
  aging_430_gore_aft_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,3,5), 430);
  aging_430_gore_aft_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,3,5), 430);
  aging_430_gore_aft_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,3,5), 430);
  aging_430_gore_aft_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,3,5), 430);
  aging_430_gore_aft_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,3,5), 430);
  aging_430_gore_aft_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,3,5), 430);
  aging_430_gore_aft_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,3,5), 430);
  aging_430_gore_aft_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,3,5), 430);
  aging_430_gore_aft_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,3,5), 430);
  aging_430_gore_aft_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,3,5), 430);
  aging_430_gore_aft_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,3,5), 430);

  Float_t aging_440_gore_a[ngorepoints];
  std::copy(aging_440_gore.begin(), aging_440_gore.end(), aging_440_gore_a);
  Float_t aging_440_gore_err[ngorepoints];
  std::fill(aging_440_gore_err, aging_440_gore_err+5, gore_err[6]);
  aging_440_gore_err[5] = stdev_error(input_hists, mackenzie_goreDRP_1018, 440);
  aging_440_gore_err[6] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 440);
  aging_440_gore_err[7] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 440);
  std::cout << "Gore error pt 8, 440nm = " << aging_440_gore_err[7] << std::endl;
  aging_440_gore_err[8] = stdev_error(input_hists, mackenzie_goreDRP_1108, 440);
  aging_440_gore_err[9] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 440);
  aging_440_gore_err[10] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 440);
  aging_440_gore_err[11] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 440);
  aging_440_gore_err[12] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 440);
  aging_440_gore_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 440);
  aging_440_gore_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 440);
  aging_440_gore_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 440);
  aging_440_gore_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 440);
  aging_440_gore_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 440);
  aging_440_gore_err[18] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 440);
  aging_440_gore_err[19] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 440);
  aging_440_gore_err[20] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 440);
  aging_440_gore_err[21] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 440);
  //aging_440_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 440);
  aging_440_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 440);
  aging_440_gore_err[23] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 440);
  Float_t aging_440_gore_bef_err[nmultgorepoints];
  Float_t aging_440_gore_aft_err[nmultgorepoints];
  aging_440_gore_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 440);
  aging_440_gore_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 440);
  aging_440_gore_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 440);
  aging_440_gore_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 440);
  aging_440_gore_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 440);
  aging_440_gore_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 440);
  aging_440_gore_bef_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 440);
  aging_440_gore_bef_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 440);
  aging_440_gore_bef_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 440);
  aging_440_gore_bef_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 440);
  aging_440_gore_bef_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 440);
  aging_440_gore_bef_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 440);
  aging_440_gore_bef_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 440);
  aging_440_gore_bef_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 440);
  aging_440_gore_bef_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 440);
  aging_440_gore_bef_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 440);
  aging_440_gore_bef_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 440);
  aging_440_gore_bef_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 440);
  aging_440_gore_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,3,5), 440);
  aging_440_gore_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,3,5), 440);
  aging_440_gore_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,3,5), 440);
  aging_440_gore_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,3,5), 440);
  aging_440_gore_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,3,5), 440);
  aging_440_gore_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,3,5), 440);
  aging_440_gore_aft_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,3,5), 440);
  aging_440_gore_aft_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,3,5), 440);
  aging_440_gore_aft_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,3,5), 440);
  aging_440_gore_aft_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,3,5), 440);
  aging_440_gore_aft_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,3,5), 440);
  aging_440_gore_aft_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,3,5), 440);
  aging_440_gore_aft_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,3,5), 440);
  aging_440_gore_aft_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,3,5), 440);
  aging_440_gore_aft_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,3,5), 440);
  aging_440_gore_aft_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,3,5), 440);
  aging_440_gore_aft_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,3,5), 440);
  aging_440_gore_aft_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,3,5), 440);

  Float_t aging_450_gore_a[ngorepoints];
  std::copy(aging_450_gore.begin(), aging_450_gore.end(), aging_450_gore_a);
  Float_t aging_450_gore_err[ngorepoints];
  std::fill(aging_450_gore_err, aging_450_gore_err+5, gore_err[7]);
  aging_450_gore_err[5] = stdev_error(input_hists, mackenzie_goreDRP_1018, 450);
  aging_450_gore_err[6] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 450);
  aging_450_gore_err[7] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 450);
  aging_450_gore_err[8] = stdev_error(input_hists, mackenzie_goreDRP_1108, 450);
  aging_450_gore_err[9] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 450);
  aging_450_gore_err[10] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 450);
  aging_450_gore_err[11] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 450);
  aging_450_gore_err[12] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 450);
  aging_450_gore_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 450);
  aging_450_gore_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 450);
  aging_450_gore_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 450);
  aging_450_gore_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 450);
  aging_450_gore_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 450);
  aging_450_gore_err[18] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 450);
  aging_450_gore_err[19] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 450);
  aging_450_gore_err[20] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 450);
  aging_450_gore_err[21] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 450);
  //aging_450_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 450);
  aging_450_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 450);
  aging_450_gore_err[23] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 450);
  Float_t aging_450_gore_bef_err[nmultgorepoints];
  Float_t aging_450_gore_aft_err[nmultgorepoints];
  aging_450_gore_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 450);
  aging_450_gore_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 450);
  aging_450_gore_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 450);
  aging_450_gore_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 450);
  aging_450_gore_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 450);
  aging_450_gore_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 450);
  aging_450_gore_bef_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 450);
  aging_450_gore_bef_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 450);
  aging_450_gore_bef_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 450);
  aging_450_gore_bef_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 450);
  aging_450_gore_bef_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 450);
  aging_450_gore_bef_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 450);
  aging_450_gore_bef_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 450);
  aging_450_gore_bef_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 450);
  aging_450_gore_bef_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 450);
  aging_450_gore_bef_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 450);
  aging_450_gore_bef_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 450);
  aging_450_gore_bef_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 450);
  aging_450_gore_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,3,5), 450);
  aging_450_gore_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,3,5), 450);
  aging_450_gore_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,3,5), 450);
  aging_450_gore_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,3,5), 450);
  aging_450_gore_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,3,5), 450);
  aging_450_gore_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,3,5), 450);
  aging_450_gore_aft_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,3,5), 450);
  aging_450_gore_aft_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,3,5), 450);
  aging_450_gore_aft_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,3,5), 450);
  aging_450_gore_aft_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,3,5), 450);
  aging_450_gore_aft_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,3,5), 450);
  aging_450_gore_aft_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,3,5), 450);
  aging_450_gore_aft_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,3,5), 450);
  aging_450_gore_aft_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,3,5), 450);
  aging_450_gore_aft_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,3,5), 450);
  aging_450_gore_aft_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,3,5), 450);
  aging_450_gore_aft_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,3,5), 450);
  aging_450_gore_aft_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,3,5), 450);

  Float_t aging_600_gore_a[ngorepoints];
  std::copy(aging_600_gore.begin(), aging_600_gore.end(), aging_600_gore_a);
  Float_t aging_600_gore_err[ngorepoints];
  std::fill(aging_600_gore_err, aging_600_gore_err+5, gore_err[8]);
  aging_600_gore_err[5] = stdev_error(input_hists, mackenzie_goreDRP_1018, 600);
  aging_600_gore_err[6] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 600);
  aging_600_gore_err[7] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 600);
  aging_600_gore_err[8] = stdev_error(input_hists, mackenzie_goreDRP_1108, 600);
  aging_600_gore_err[9] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 600);
  aging_600_gore_err[10] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 600);
  aging_600_gore_err[11] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 600);
  aging_600_gore_err[12] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 600);
  aging_600_gore_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 600);
  aging_600_gore_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 600);
  aging_600_gore_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 600);
  aging_600_gore_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 600);
  aging_600_gore_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 600);
  aging_600_gore_err[18] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 600);
  aging_600_gore_err[19] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 600);
  aging_600_gore_err[20] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 600);
  aging_600_gore_err[21] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 600);
  //aging_600_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 600);
  aging_600_gore_err[22] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 600);
  aging_600_gore_err[23] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 600);
  Float_t aging_600_gore_bef_err[nmultgorepoints];
  Float_t aging_600_gore_aft_err[nmultgorepoints];
  aging_600_gore_bef_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,0,2), 600);
  aging_600_gore_bef_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,0,2), 600);
  aging_600_gore_bef_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,0,2), 600);
  aging_600_gore_bef_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,0,2), 600);
  aging_600_gore_bef_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,0,2), 600);
  aging_600_gore_bef_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,0,2), 600);
  aging_600_gore_bef_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,0,2), 600);
  aging_600_gore_bef_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,0,2), 600);
  aging_600_gore_bef_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,0,2), 600);
  aging_600_gore_bef_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,0,2), 600);
  aging_600_gore_bef_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,0,2), 600);
  aging_600_gore_bef_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,0,2), 600);
  aging_600_gore_bef_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,0,2), 600);
  aging_600_gore_bef_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,0,2), 600);
  aging_600_gore_bef_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,0,2), 600);
  aging_600_gore_bef_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,0,2), 600);
  aging_600_gore_bef_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,0,2), 600);
  aging_600_gore_bef_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,0,2), 600);
  aging_600_gore_aft_err[0] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1025,3,5), 600);
  aging_600_gore_aft_err[1] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1101,3,5), 600);
  aging_600_gore_aft_err[2] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1115,3,5), 600);
  aging_600_gore_aft_err[3] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1122,3,5), 600);
  aging_600_gore_aft_err[4] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1129,3,5), 600);
  aging_600_gore_aft_err[5] = stdev_error(input_hists, slicing(mackenzie_goreDRP_1206,3,5), 600);
  aging_600_gore_aft_err[6] = stdev_error(input_hists, slicing(brian_goreDRP_020123,3,5), 600);
  aging_600_gore_aft_err[7] = stdev_error(input_hists, slicing(brian_goreDRP_020923,3,5), 600);
  aging_600_gore_aft_err[8] = stdev_error(input_hists, slicing(brian_goreDRP_021723,3,5), 600);
  aging_600_gore_aft_err[9] = stdev_error(input_hists, slicing(brian_goreDRP_030923,3,5), 600);
  aging_600_gore_aft_err[10] = stdev_error(input_hists, slicing(brian_goreDRP_031623,3,5), 600);
  aging_600_gore_aft_err[11] = stdev_error(input_hists, slicing(brian_goreDRP_032323,3,5), 600);
  aging_600_gore_aft_err[12] = stdev_error(input_hists, slicing(brian_goreDRP_040623,3,5), 600);
  aging_600_gore_aft_err[13] = stdev_error(input_hists, slicing(brian_goreDRP_041323,3,5), 600);
  aging_600_gore_aft_err[14] = stdev_error(input_hists, slicing(brian_goreDRP_042123,3,5), 600);
  aging_600_gore_aft_err[15] = stdev_error(input_hists, slicing(brian_goreDRP_052423,3,5), 600);
  aging_600_gore_aft_err[16] = stdev_error(input_hists, slicing(brian_goreDRP_060823,3,5), 600);
  aging_600_gore_aft_err[17] = stdev_error(input_hists, slicing(brian_goreDRP_062323,3,5), 600);

  std::cout << "We are getting all goreDRP errors." << std::endl;

  std::vector<Float_t> aging_slopes;
  std::vector<Float_t> aging_slopes_err;
  std::vector<Float_t> aging_slopes_gore;
  std::vector<Float_t> aging_slopes_gore_err;

  page[0]->cd();
  nova_380 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_380_a, empty_x_err, aging_380_err);
  nova_380->SetName("tg_nova_380");
  /*std::cout << "printing 380nm array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_380_a[i] << ", ";
  }
  std::cout << "\n";*/
  nova_380->SetTitle("R at 380nm of NOvA Standard vs Time");
  nova_380->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_380->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_380->SetMarkerStyle(20);
  nova_380->Fit("pol1");
  nova_380->Draw("AP");
  TF1 *fit0 = (TF1*)nova_380->GetListOfFunctions()->FindObject("pol1");
  aging_slopes.push_back(fit0->GetParameter(1));
  aging_slopes_err.push_back(fit0->GetParError(1));
  std::cout << "380nm slope error = " << fit0->GetParError(1) << std::endl;
  page[0]->Update();
  ofile->cd();
  page[0]->Write();
  page[0]->Close();

  for(auto i : aging_slopes_err) std::cout << "errors = " << i << std::endl;

  page[1]->cd();
  nova_390 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_390_a, empty_x_err, aging_390_err);
  nova_390->SetName("tg_nova_390");
  /*std::cout << "printing 390nm array again, after TG plotting line: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_390_a[i] << ", ";
  }
  std::cout << "\n";*/
  nova_390->SetTitle("R at 390nm of NOvA Standard vs Time");
  nova_390->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_390->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_390->SetMarkerStyle(20);
  nova_390->Fit("pol1");
  nova_390->Draw("AP");
  TF1 *fit1 = (TF1*)nova_390->GetListOfFunctions()->FindObject("pol1");
  aging_slopes.push_back(fit1->GetParameter(1));
  aging_slopes_err.push_back(fit1->GetParError(1));
  page[1]->Update();
  ofile->cd();
  page[1]->Write();
  page[1]->Close();

  page[2]->cd();
  nova_400 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_400_a, empty_x_err, aging_400_err);
  nova_400->SetName("tg_nova_400");
  /*std::cout << "printing 400nm array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_400_a[i] << ", ";
  }
  std::cout << "\n";*/
  nova_400->SetTitle("R at 400nm of NOvA Standard vs Time");
  nova_400->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_400->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_400->SetMarkerStyle(20);
  nova_400->Fit("pol1");
  nova_400->Draw("AP");
  TF1 *fit2 = (TF1*)nova_400->GetListOfFunctions()->FindObject("pol1");
  aging_slopes.push_back(fit2->GetParameter(1));
  aging_slopes_err.push_back(fit2->GetParError(1));
  page[2]->Update();
  ofile->cd();
  page[2]->Write();
  page[2]->Close();

  page[3]->cd();
  nova_410 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_410_a, empty_x_err, aging_410_err);
  nova_410->SetName("tg_nova_410");
  /*std::cout << "printing 410nm array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_410_a[i] << ", ";
  }
  std::cout << "\n";*/
  nova_410->SetTitle("R at 410nm of NOvA Standard vs Time");
  nova_410->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_410->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_410->SetMarkerStyle(20);
  nova_410->Fit("pol1");
  nova_410->Draw("AP");
  TF1 *fit3 = (TF1*)nova_410->GetListOfFunctions()->FindObject("pol1");
  aging_slopes.push_back(fit3->GetParameter(1));
  aging_slopes_err.push_back(fit3->GetParError(1));
  page[3]->Update();
  ofile->cd();
  page[3]->Write();
  page[3]->Close();

  page[4]->cd();
  nova_420 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_420_a, empty_x_err, aging_420_err);
  nova_420->SetName("tg_nova_420");
  /*std::cout << "printing 420nm array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_420_a[i] << ", ";
  }
  std::cout << "\n";*/
  nova_420->SetTitle("R at 420nm of NOvA Standard vs Time");
  nova_420->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_420->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_420->SetMarkerStyle(20);
  nova_420->Fit("pol1");
  nova_420->Draw("AP");
  TF1 *fit4 = (TF1*)nova_420->GetListOfFunctions()->FindObject("pol1");
  aging_slopes.push_back(fit4->GetParameter(1));
  aging_slopes_err.push_back(fit4->GetParError(1));
  page[4]->Update();
  ofile->cd();
  page[4]->Write();
  page[4]->Close();

  page[5]->cd();
  nova_430 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_430_a, empty_x_err, aging_430_err);
  nova_430->SetName("tg_nova_430");
  /*std::cout << "printing 430nm array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_430_a[i] << ", ";
  }
  std::cout << "\n";*/
  nova_430->SetTitle("R at 430nm of NOvA Standard vs Time");
  nova_430->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_430->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_430->SetMarkerStyle(20);
  nova_430->Fit("pol1");
  nova_430->Draw("AP");
  TF1 *fit5 = (TF1*)nova_430->GetListOfFunctions()->FindObject("pol1");
  aging_slopes.push_back(fit5->GetParameter(1));
  aging_slopes_err.push_back(fit5->GetParError(1));
  page[5]->Update();
  ofile->cd();
  page[5]->Write();
  page[5]->Close();

  page[6]->cd();
  nova_440 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_440_a, empty_x_err, aging_440_err);
  nova_440->SetName("tg_nova_440");
  /*std::cout << "printing 440nm array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_440_a[i] << ", ";
  }
  std::cout << "\n";*/
  nova_440->SetTitle("R at 440nm of NOvA Standard vs Time");
  nova_440->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_440->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_440->SetMarkerStyle(20);
  nova_440->Fit("pol1");
  nova_440->Draw("AP");
  TF1 *fit6 = (TF1*)nova_440->GetListOfFunctions()->FindObject("pol1");
  aging_slopes.push_back(fit6->GetParameter(1));
  aging_slopes_err.push_back(fit6->GetParError(1));
  page[6]->Update();
  ofile->cd();
  page[6]->Write();
  page[6]->Close();

  page[7]->cd();
  nova_450 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_450_a, empty_x_err, aging_450_err);
  nova_450->SetName("tg_nova_450");
  /*std::cout << "printing 450nm array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_450_a[i] << ", ";
  }
  std::cout << "\n";*/
  nova_450->SetTitle("R at 450nm of NOvA Standard vs Time");
  nova_450->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_450->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_450->SetMarkerStyle(20);
  nova_450->Fit("pol1");
  nova_450->Draw("AP");
  TF1 *fit7 = (TF1*)nova_450->GetListOfFunctions()->FindObject("pol1");
  aging_slopes.push_back(fit7->GetParameter(1));
  aging_slopes_err.push_back(fit7->GetParError(1));
  page[7]->Update();
  ofile->cd();
  page[7]->Write();
  page[7]->Close();

  page600[0]->cd();
  nova_600 = new TGraphErrors(npoints, aging_timefrac_xbins, aging_600_a, empty_x_err, aging_600_err);
  nova_600->SetName("tg_nova_600");
  std::cout << "printing 600nm array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_600_a[i] << ", ";
  }
  std::cout << "\n";
  std::cout << "printing 600nm error array: ";
  for(int i = 0; i < npoints; i++){
    std::cout << aging_600_err[i] << ", ";
  }
  std::cout << "\n";
  nova_600->SetTitle("R at 600nm of NOvA Standard vs Time");
  nova_600->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_600->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_600->SetMarkerStyle(20);
  nova_600->Fit("pol1");
  nova_600->Draw("AP");
  page600[0]->Update();
  ofile->cd();
  page600[0]->Write();
  page600[0]->Close();

  // before points
  page600[1]->cd();
  nova_600_bef = new TGraph(nmultpoints, befaft_time, aging_600_bef_a);
  nova_600_bef->SetName("tg_nova_600_bef");
  nova_600_bef->SetTitle("R at 600nm of NOvA Standard vs Time");
  nova_600_bef->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_600_bef->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_600_bef->SetMarkerStyle(20);
  nova_600_bef->Fit("pol1");
  nova_600_bef->Draw("AP");
  page600[1]->Update();
  ofile->cd();
  page600[1]->Write();
  page600[1]->Close();

  // after points
  page600[2]->cd();
  nova_600_aft = new TGraph(nmultpoints, befaft_time, aging_600_aft_a);
  nova_600_aft->SetName("tg_nova_600_aft");
  nova_600_aft->SetTitle("R at 600nm of NOvA Standard vs Time");
  nova_600_aft->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  nova_600_aft->GetYaxis()->SetTitle("Reflectance (R) (%)");
  nova_600_aft->SetMarkerStyle(20);
  nova_600_aft->Fit("pol1");
  nova_600_aft->Draw("AP");
  page600[2]->Update();
  ofile->cd();
  page600[2]->Write();
  page600[2]->Close();

  for(auto i : aging_slopes_err) std::cout << "errors after all plots = " << i << std::endl;

  // now for goreDRP aging plots

  gorepage[0]->cd();
  gore_380 = new TGraphErrors(ngorepoints, gore_time, aging_380_gore_a, empty_x_err_gore, aging_380_gore_err);
  gore_380->SetName("tg_gore_380");
  gore_380->SetTitle("R at 380nm of GoreDRP Standard vs Time");
  gore_380->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  gore_380->GetYaxis()->SetTitle("Reflectance (R) (%)");
  gore_380->SetMarkerStyle(20);
  gore_380->Fit("pol1");
  gore_380->Draw("AP");
  TF1 *fitg0 = (TF1*)gore_380->GetListOfFunctions()->FindObject("pol1");
  aging_slopes_gore.push_back(fitg0->GetParameter(1));
  aging_slopes_gore_err.push_back(fitg0->GetParError(1));
  gorepage[0]->Update();
  ofile->cd();
  gorepage[0]->Write();
  gorepage[0]->Close();

  gorepage[1]->cd();
  gore_390 = new TGraphErrors(ngorepoints, gore_time, aging_390_gore_a, empty_x_err_gore, aging_390_gore_err);
  gore_390->SetName("tg_gore_390");
  gore_390->SetTitle("R at 390nm of GoreDRP Standard vs Time");
  gore_390->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  gore_390->GetYaxis()->SetTitle("Reflectance (R) (%)");
  gore_390->SetMarkerStyle(20);
  gore_390->Fit("pol1");
  gore_390->Draw("AP");
  TF1 *fitg1 = (TF1*)gore_390->GetListOfFunctions()->FindObject("pol1");
  aging_slopes_gore.push_back(fitg1->GetParameter(1));
  aging_slopes_gore_err.push_back(fitg1->GetParError(1));
  gorepage[1]->Update();
  ofile->cd();
  gorepage[1]->Write();
  gorepage[1]->Close();

  gorepage[2]->cd();
  gore_400 = new TGraphErrors(ngorepoints, gore_time, aging_400_gore_a, empty_x_err_gore, aging_400_gore_err);
  gore_400->SetName("tg_gore_400");
  gore_400->SetTitle("R at 400nm of GoreDRP Standard vs Time");
  gore_400->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  gore_400->GetYaxis()->SetTitle("Reflectance (R) (%)");
  gore_400->SetMarkerStyle(20);
  gore_400->Fit("pol1");
  gore_400->Draw("AP");
  TF1 *fitg2 = (TF1*)gore_400->GetListOfFunctions()->FindObject("pol1");
  aging_slopes_gore.push_back(fitg2->GetParameter(1));
  aging_slopes_gore_err.push_back(fitg2->GetParError(1));
  gorepage[2]->Update();
  ofile->cd();
  gorepage[2]->Write();
  gorepage[2]->Close();

  gorepage[3]->cd();
  gore_410 = new TGraphErrors(ngorepoints, gore_time, aging_410_gore_a, empty_x_err_gore, aging_410_gore_err);
  gore_410->SetName("tg_gore_410");
  gore_410->SetTitle("R at 410nm of GoreDRP Standard vs Time");
  gore_410->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  gore_410->GetYaxis()->SetTitle("Reflectance (R) (%)");
  gore_410->SetMarkerStyle(20);
  gore_410->Fit("pol1");
  gore_410->Draw("AP");
  TF1 *fitg3 = (TF1*)gore_410->GetListOfFunctions()->FindObject("pol1");
  aging_slopes_gore.push_back(fitg3->GetParameter(1));
  aging_slopes_gore_err.push_back(fitg3->GetParError(1));
  gorepage[3]->Update();
  ofile->cd();
  gorepage[3]->Write();
  gorepage[3]->Close();

  gorepage[4]->cd();
  gore_420 = new TGraphErrors(ngorepoints, gore_time, aging_420_gore_a, empty_x_err_gore, aging_420_gore_err);
  gore_420->SetName("tg_gore_420");
  gore_420->SetTitle("R at 420nm of GoreDRP Standard vs Time");
  gore_420->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  gore_420->GetYaxis()->SetTitle("Reflectance (R) (%)");
  gore_420->SetMarkerStyle(20);
  gore_420->Fit("pol1");
  gore_420->Draw("AP");
  TF1 *fitg4 = (TF1*)gore_420->GetListOfFunctions()->FindObject("pol1");
  aging_slopes_gore.push_back(fitg4->GetParameter(1));
  aging_slopes_gore_err.push_back(fitg4->GetParError(1));
  gorepage[4]->Update();
  ofile->cd();
  gorepage[4]->Write();
  gorepage[4]->Close();

  gorepage[5]->cd();
  gore_430 = new TGraphErrors(ngorepoints, gore_time, aging_430_gore_a, empty_x_err_gore, aging_430_gore_err);
  gore_430->SetName("tg_gore_430");
  gore_430->SetTitle("R at 430nm of GoreDRP Standard vs Time");
  gore_430->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  gore_430->GetYaxis()->SetTitle("Reflectance (R) (%)");
  gore_430->SetMarkerStyle(20);
  gore_430->Fit("pol1");
  gore_430->Draw("AP");
  TF1 *fitg5 = (TF1*)gore_430->GetListOfFunctions()->FindObject("pol1");
  aging_slopes_gore.push_back(fitg5->GetParameter(1));
  aging_slopes_gore_err.push_back(fitg5->GetParError(1));
  gorepage[5]->Update();
  ofile->cd();
  gorepage[5]->Write();
  gorepage[5]->Close();

  gorepage[6]->cd();
  gore_440 = new TGraphErrors(ngorepoints, gore_time, aging_440_gore_a, empty_x_err_gore, aging_440_gore_err);
  gore_440->SetName("tg_gore_440");
  gore_440->SetTitle("R at 440nm of GoreDRP Standard vs Time");
  gore_440->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  gore_440->GetYaxis()->SetTitle("Reflectance (R) (%)");
  gore_440->SetMarkerStyle(20);
  gore_440->Fit("pol1");
  gore_440->Draw("AP");
  TF1 *fitg6 = (TF1*)gore_440->GetListOfFunctions()->FindObject("pol1");
  aging_slopes_gore.push_back(fitg6->GetParameter(1));
  aging_slopes_gore_err.push_back(fitg6->GetParError(1));
  gorepage[6]->Update();
  ofile->cd();
  gorepage[6]->Write();
  gorepage[6]->Close();

  gorepage[7]->cd();
  gore_450 = new TGraphErrors(ngorepoints, gore_time, aging_450_gore_a, empty_x_err_gore, aging_450_gore_err);
  gore_450->SetName("tg_gore_450");
  gore_450->SetTitle("R at 450nm of GoreDRP Standard vs Time");
  gore_450->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  gore_450->GetYaxis()->SetTitle("Reflectance (R) (%)");
  gore_450->SetMarkerStyle(20);
  gore_450->Fit("pol1");
  gore_450->Draw("AP");
  TF1 *fitg7 = (TF1*)gore_450->GetListOfFunctions()->FindObject("pol1");
  aging_slopes_gore.push_back(fitg7->GetParameter(1));
  aging_slopes_gore_err.push_back(fitg7->GetParError(1));
  gorepage[7]->Update();
  ofile->cd();
  gorepage[7]->Write();
  gorepage[7]->Close();

  gorepage[8]->cd();
  gore_600 = new TGraphErrors(ngorepoints, gore_time, aging_600_gore_a, empty_x_err_gore, aging_600_gore_err);
  gore_600->SetName("tg_gore_600");
  gore_600->SetTitle("R at 600nm of GoreDRP Standard vs Time");
  gore_600->GetXaxis()->SetTitle("Time Fraction (% of Years Since 3/17/22)");
  gore_600->GetYaxis()->SetTitle("Reflectance (R) (%)");
  gore_600->SetMarkerStyle(20);
  gore_600->Fit("pol1");
  gore_600->Draw("AP");
  TF1 *fitg8 = (TF1*)gore_600->GetListOfFunctions()->FindObject("pol1");
  aging_slopes_gore.push_back(fitg8->GetParameter(1));
  gorepage[8]->Update();
  ofile->cd();
  gorepage[8]->Write();
  gorepage[8]->Close();

  // slopes as a function of wavelength

  TCanvas* slopes = new TCanvas("nova_aging_slope_vs_wavelength","nova_aging_slope_vs_wavelength", 900, 600);
  TCanvas* gore_slopes = new TCanvas("gore_aging_slope_vs_wavelength","gore_aging_slope_vs_wavelength", 900, 600);

  TH1F* hslopes1D_nova;
  TH1F* hslopes1D_gore;
  hslopes1D_nova = new TH1F("hslopes1D_nova", "hslopes1D_nova", 100, -1., 1.);
  hslopes1D_nova->GetXaxis()->SetTitle("Nova N-27-09-NC Slopes (Percent of R Lost per Year)");
  hslopes1D_nova->GetYaxis()->SetTitle("Counts");

  hslopes1D_gore = new TH1F("hslopes1D_gore", "hslopes1D_gore", 100, -1., 1.);
  hslopes1D_gore->GetXaxis()->SetTitle("GoreDRP Slopes (Percent of R Lost per Year)");
  hslopes1D_gore->GetYaxis()->SetTitle("Counts");

  Float_t aging_slopes_a[8];
  std::copy(aging_slopes.begin(), aging_slopes.end(), aging_slopes_a);
  Float_t aging_slopes_err_a[8];
  std::copy(aging_slopes_err.begin(), aging_slopes_err.end(), aging_slopes_err_a);
  Float_t aging_slopes_gore_a[8];
  std::copy(aging_slopes_gore.begin(), aging_slopes_gore.end(), aging_slopes_gore_a);
  Float_t aging_slopes_gore_err_a[8];
  std::copy(aging_slopes_gore_err.begin(), aging_slopes_gore_err.end(), aging_slopes_gore_err_a);
  Float_t wavelengths[8] = { 380., 390., 400., 410., 420., 430., 440., 450. };
  Float_t wavelength_err[8] = {0};

  for(int i; i < 8; i++) std::cout << "error array = " << aging_slopes_err_a[i] << std::endl;

  for(int i = 0; i < 8; i++){
    hslopes1D_nova->Fill(aging_slopes_a[i]);
    hslopes1D_gore->Fill(aging_slopes_gore_a[i]);
  }
  hslopes1D_nova->Write();
  hslopes1D_gore->Write();

  slopes->cd();
  TGraphErrors* tgaging_slopes = new TGraphErrors(8, wavelengths, aging_slopes_a, wavelength_err, aging_slopes_err_a);
  tgaging_slopes->SetName("tg_aging_slopes");
  tgaging_slopes->SetMarkerStyle(20);
  tgaging_slopes->SetTitle("Aging Slopes vs Wavelength, NOvA Standard N-27-09-NC");
  tgaging_slopes->GetXaxis()->SetTitle("Wavelength (nm)");
  tgaging_slopes->GetYaxis()->SetTitle("Aging Slope (R% Loss/Year)");
  tgaging_slopes->Fit("pol0");
  tgaging_slopes->Draw("AP");
  slopes->Update();
  ofile->cd();
  slopes->Write();
  slopes->Close();

  gore_slopes->cd();
  TGraphErrors* tgaging_slopes_gore = new TGraphErrors(8, wavelengths, aging_slopes_gore_a, wavelength_err, aging_slopes_gore_err_a);
  tgaging_slopes_gore->SetName("tg_aging_slopes_gore");
  tgaging_slopes_gore->SetMarkerStyle(20);
  tgaging_slopes_gore->SetTitle("Aging Slopes vs Wavelength, GoreDRP Standard");
  tgaging_slopes_gore->GetXaxis()->SetTitle("Wavelength (nm)");
  tgaging_slopes_gore->GetYaxis()->SetTitle("Aging Slope (R% Loss/Year)");
  tgaging_slopes_gore->Fit("pol0");
  tgaging_slopes_gore->Draw("AP");
  gore_slopes->Update();
  ofile->cd();
  gore_slopes->Write();
  gore_slopes->Close();

  stability[0]->cd();
  nova_stability = new TGraph(8, wavelengths, nova_err);
  nova_stability->SetName("tg_nova_stability");
  nova_stability->SetMarkerStyle(20);
  nova_stability->SetTitle("Standard Deviation of Repeated Nova Measurements");
  nova_stability->GetXaxis()->SetTitle("Wavelength (nm)");
  nova_stability->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  nova_stability->Draw("AP");
  stability[0]->Update();
  ofile->cd();
  stability[0]->Write();
  stability[0]->Close();

  stability[1]->cd();
  gore_stability = new TGraph(8, wavelengths, gore_err);
  gore_stability->SetName("tg_gore_stability");
  gore_stability->SetMarkerStyle(20);
  gore_stability->SetTitle("Standard Deviation of Repeated GoreDRP Measurements");
  gore_stability->GetXaxis()->SetTitle("Wavelength (nm)");
  gore_stability->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  gore_stability->Draw("AP");
  stability[1]->Update();
  ofile->cd();
  stability[1]->Write();
  stability[1]->Close();

  beforeafter[0]->cd();
  nova_bef_380 = new TGraph(nmultpoints, befaft_time, nova_380_bef_err);
  nova_bef_380->SetName("tg_nova_bef_380");
  nova_bef_380->SetMarkerStyle(24);
  nova_bef_380->SetTitle("Standard Deviation of Repeated Nova Measurements");
  nova_bef_380->GetXaxis()->SetTitle("Time, 380nm (Fraction of Years)");
  nova_bef_380->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  nova_bef_380->Fit("pol0");
  nova_bef_380->Draw("AP");
  nova_aft_380 = new TGraph(nmultpoints, befaft_time, nova_380_aft_err);
  nova_aft_380->SetName("tg_nova_aft_380");
  nova_aft_380->SetMarkerStyle(20);
  nova_aft_380->Fit("pol0");
  nova_aft_380->GetFunction("pol0")->SetLineColor(kBlue);
  nova_aft_380->Draw("P");
  beforeafter[0]->Update();
  ofile->cd();
  beforeafter[0]->Write();
  beforeafter[0]->Close();

  beforeafter[1]->cd();
  nova_bef_390 = new TGraph(nmultpoints, befaft_time, nova_390_bef_err);
  nova_bef_390->SetName("tg_nova_bef_390");
  nova_bef_390->SetMarkerStyle(24);
  nova_bef_390->SetTitle("Standard Deviation of Repeated Nova Measurements");
  nova_bef_390->GetXaxis()->SetTitle("Time, 390nm (Fraction of Years)");
  nova_bef_390->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  nova_bef_390->Fit("pol0");
  nova_bef_390->Draw("AP");
  nova_aft_390 = new TGraph(nmultpoints, befaft_time, nova_390_aft_err);
  nova_aft_390->SetName("tg_nova_aft_390");
  nova_aft_390->SetMarkerStyle(20);
  nova_aft_390->Fit("pol0");
  nova_aft_390->GetFunction("pol0")->SetLineColor(kBlue);
  nova_aft_390->Draw("P");
  beforeafter[1]->Update();
  ofile->cd();
  beforeafter[1]->Write();
  beforeafter[1]->Close();

  beforeafter[2]->cd();
  nova_bef_400 = new TGraph(nmultpoints, befaft_time, nova_400_bef_err);
  nova_bef_400->SetName("tg_nova_bef_400");
  nova_bef_400->SetMarkerStyle(24);
  nova_bef_400->SetTitle("Standard Deviation of Repeated Nova Measurements");
  nova_bef_400->GetXaxis()->SetTitle("Time, 400nm (Fraction of Years)");
  nova_bef_400->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  nova_bef_400->Fit("pol0");
  nova_bef_400->Draw("AP");
  nova_aft_400 = new TGraph(nmultpoints, befaft_time, nova_400_aft_err);
  nova_aft_400->SetName("tg_nova_aft_400");
  nova_aft_400->SetMarkerStyle(20);
  nova_aft_400->Fit("pol0");
  nova_aft_400->GetFunction("pol0")->SetLineColor(kBlue);
  nova_aft_400->Draw("P");
  beforeafter[2]->Update();
  ofile->cd();
  beforeafter[2]->Write();
  beforeafter[2]->Close();

  beforeafter[3]->cd();
  nova_bef_410 = new TGraph(nmultpoints, befaft_time, nova_410_bef_err);
  nova_bef_410->SetName("tg_nova_bef_410");
  nova_bef_410->SetMarkerStyle(24);
  nova_bef_410->SetTitle("Standard Deviation of Repeated Nova Measurements");
  nova_bef_410->GetXaxis()->SetTitle("Time, 410nm (Fraction of Years)");
  nova_bef_410->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  nova_bef_410->Fit("pol0");
  nova_bef_410->Draw("AP");
  nova_aft_410 = new TGraph(nmultpoints, befaft_time, nova_410_aft_err);
  nova_aft_410->SetName("tg_nova_aft_410");
  nova_aft_410->SetMarkerStyle(20);
  nova_aft_410->Fit("pol0");
  nova_aft_410->GetFunction("pol0")->SetLineColor(kBlue);
  nova_aft_410->Draw("P");
  beforeafter[3]->Update();
  ofile->cd();
  beforeafter[3]->Write();
  beforeafter[3]->Close();

  beforeafter[4]->cd();
  nova_bef_420 = new TGraph(nmultpoints, befaft_time, nova_420_bef_err);
  nova_bef_420->SetName("tg_nova_bef_420");
  nova_bef_420->SetMarkerStyle(24);
  nova_bef_420->SetTitle("Standard Deviation of Repeated Nova Measurements");
  nova_bef_420->GetXaxis()->SetTitle("Time, 420nm (Fraction of Years)");
  nova_bef_420->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  nova_bef_420->Fit("pol0");
  nova_bef_420->Draw("AP");
  nova_aft_420 = new TGraph(nmultpoints, befaft_time, nova_420_aft_err);
  nova_aft_420->SetName("tg_nova_aft_420");
  nova_aft_420->SetMarkerStyle(20);
  nova_aft_420->Fit("pol0");
  nova_aft_420->GetFunction("pol0")->SetLineColor(kBlue);
  nova_aft_420->Draw("P");
  beforeafter[4]->Update();
  ofile->cd();
  beforeafter[4]->Write();
  beforeafter[4]->Close();

  beforeafter[5]->cd();
  nova_bef_430 = new TGraph(nmultpoints, befaft_time, nova_430_bef_err);
  nova_bef_430->SetName("tg_nova_bef_430");
  nova_bef_430->SetMarkerStyle(24);
  nova_bef_430->SetTitle("Standard Deviation of Repeated Nova Measurements");
  nova_bef_430->GetXaxis()->SetTitle("Time, 430nm (Fraction of Years)");
  nova_bef_430->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  nova_bef_430->Fit("pol0");
  nova_bef_430->Draw("AP");
  nova_aft_430 = new TGraph(nmultpoints, befaft_time, nova_430_aft_err);
  nova_aft_430->SetName("tg_nova_aft_430");
  nova_aft_430->SetMarkerStyle(20);
  nova_aft_430->Fit("pol0");
  nova_aft_430->GetFunction("pol0")->SetLineColor(kBlue);
  nova_aft_430->Draw("P");
  beforeafter[5]->Update();
  ofile->cd();
  beforeafter[5]->Write();
  beforeafter[5]->Close();

  beforeafter[6]->cd();
  nova_bef_440 = new TGraph(nmultpoints, befaft_time, nova_440_bef_err);
  nova_bef_440->SetName("tg_nova_bef_440");
  nova_bef_440->SetMarkerStyle(24);
  nova_bef_440->SetTitle("Standard Deviation of Repeated Nova Measurements");
  nova_bef_440->GetXaxis()->SetTitle("Time, 440nm (Fraction of Years)");
  nova_bef_440->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  nova_bef_440->Fit("pol0");
  nova_bef_440->Draw("AP");
  nova_aft_440 = new TGraph(nmultpoints, befaft_time, nova_440_aft_err);
  nova_aft_440->SetName("tg_nova_aft_440");
  nova_aft_440->SetMarkerStyle(20);
  nova_aft_440->Fit("pol0");
  nova_aft_440->GetFunction("pol0")->SetLineColor(kBlue);
  nova_aft_440->Draw("P");
  beforeafter[6]->Update();
  ofile->cd();
  beforeafter[6]->Write();
  beforeafter[6]->Close();

  beforeafter[7]->cd();
  nova_bef_450 = new TGraph(nmultpoints, befaft_time, nova_450_bef_err);
  nova_bef_450->SetName("tg_nova_bef_450");
  nova_bef_450->SetMarkerStyle(24);
  nova_bef_450->SetTitle("Standard Deviation of Repeated Nova Measurements");
  nova_bef_450->GetXaxis()->SetTitle("Time, 450nm (Fraction of Years)");
  nova_bef_450->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  nova_bef_450->Fit("pol0");
  nova_bef_450->Draw("AP");
  nova_aft_450 = new TGraph(nmultpoints, befaft_time, nova_450_aft_err);
  nova_aft_450->SetName("tg_nova_aft_450");
  nova_aft_450->SetMarkerStyle(20);
  nova_aft_450->Fit("pol0");
  nova_aft_450->GetFunction("pol0")->SetLineColor(kBlue);
  nova_aft_450->Draw("P");
  beforeafter[7]->Update();
  ofile->cd();
  beforeafter[7]->Write();
  beforeafter[7]->Close();

  beforeafter[8]->cd();
  nova_bef_600 = new TGraph(nmultpoints, befaft_time, nova_600_bef_err);
  nova_bef_600->SetName("tg_nova_bef_600");
  nova_bef_600->SetMarkerStyle(24);
  nova_bef_600->SetTitle("Standard Deviation of Repeated Nova Measurements");
  nova_bef_600->GetXaxis()->SetTitle("Time, 600nm (Fraction of Years)");
  nova_bef_600->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  nova_bef_600->Fit("pol0");
  nova_bef_600->Draw("AP");
  nova_aft_600 = new TGraph(nmultpoints, befaft_time, nova_600_aft_err);
  nova_aft_600->SetName("tg_nova_aft_600");
  nova_aft_600->SetMarkerStyle(20);
  nova_aft_600->Fit("pol0");
  nova_aft_600->GetFunction("pol0")->SetLineColor(kBlue);
  nova_aft_600->Draw("P");
  beforeafter[8]->Update();
  ofile->cd();
  beforeafter[8]->Write();
  beforeafter[8]->Close();

  beforeafter_gore[0]->cd();
  gore_bef_380 = new TGraph(nmultgorepoints, befaft_gore_time, aging_380_gore_bef_err);
  gore_bef_380->SetName("tg_gore_bef_380");
  gore_bef_380->SetMarkerStyle(24);
  gore_bef_380->SetTitle("Standard Deviation of Repeated GoreDRP Measurements");
  gore_bef_380->GetXaxis()->SetTitle("Time, 380nm (Fraction of Years)");
  gore_bef_380->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  gore_bef_380->Fit("pol0");
  gore_bef_380->Draw("AP");
  gore_aft_380 = new TGraph(nmultgorepoints, befaft_gore_time, aging_380_gore_aft_err);
  gore_aft_380->SetName("tg_gore_aft_380");
  gore_aft_380->SetMarkerStyle(20);
  gore_aft_380->Fit("pol0");
  gore_aft_380->GetFunction("pol0")->SetLineColor(kBlue);
  gore_aft_380->Draw("P");
  beforeafter_gore[0]->Update();
  ofile->cd();
  beforeafter_gore[0]->Write();
  beforeafter_gore[0]->Close();

  beforeafter_gore[1]->cd();
  gore_bef_390 = new TGraph(nmultgorepoints, befaft_gore_time, aging_390_gore_bef_err);
  gore_bef_390->SetName("tg_gore_bef_390");
  gore_bef_390->SetMarkerStyle(24);
  gore_bef_390->SetTitle("Standard Deviation of Repeated GoreDRP Measurements");
  gore_bef_390->GetXaxis()->SetTitle("Time, 390nm (Fraction of Years)");
  gore_bef_390->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  gore_bef_390->Fit("pol0");
  gore_bef_390->Draw("AP");
  gore_aft_390 = new TGraph(nmultgorepoints, befaft_gore_time, aging_390_gore_aft_err);
  gore_aft_390->SetName("tg_gore_aft_390");
  gore_aft_390->SetMarkerStyle(20);
  gore_aft_390->Fit("pol0");
  gore_aft_390->GetFunction("pol0")->SetLineColor(kBlue);
  gore_aft_390->Draw("P");
  beforeafter_gore[1]->Update();
  ofile->cd();
  beforeafter_gore[1]->Write();
  beforeafter_gore[1]->Close();

  beforeafter_gore[2]->cd();
  gore_bef_400 = new TGraph(nmultgorepoints, befaft_gore_time, aging_400_gore_bef_err);
  gore_bef_400->SetName("tg_gore_bef_400");
  gore_bef_400->SetMarkerStyle(24);
  gore_bef_400->SetTitle("Standard Deviation of Repeated GoreDRP Measurements");
  gore_bef_400->GetXaxis()->SetTitle("Time, 400nm (Fraction of Years)");
  gore_bef_400->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  gore_bef_400->Fit("pol0");
  gore_bef_400->Draw("AP");
  gore_aft_400 = new TGraph(nmultgorepoints, befaft_gore_time, aging_400_gore_aft_err);
  gore_aft_400->SetName("tg_gore_aft_400");
  gore_aft_400->SetMarkerStyle(20);
  gore_aft_400->Fit("pol0");
  gore_aft_400->GetFunction("pol0")->SetLineColor(kBlue);
  gore_aft_400->Draw("P");
  beforeafter_gore[2]->Update();
  ofile->cd();
  beforeafter_gore[2]->Write();
  beforeafter_gore[2]->Close();

  beforeafter_gore[3]->cd();
  gore_bef_410 = new TGraph(nmultgorepoints, befaft_gore_time, aging_410_gore_bef_err);
  gore_bef_410->SetName("tg_gore_bef_410");
  gore_bef_410->SetMarkerStyle(24);
  gore_bef_410->SetTitle("Standard Deviation of Repeated GoreDRP Measurements");
  gore_bef_410->GetXaxis()->SetTitle("Time, 410nm (Fraction of Years)");
  gore_bef_410->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  gore_bef_410->Fit("pol0");
  gore_bef_410->Draw("AP");
  gore_aft_410 = new TGraph(nmultgorepoints, befaft_gore_time, aging_410_gore_aft_err);
  gore_aft_410->SetName("tg_gore_aft_410");
  gore_aft_410->SetMarkerStyle(20);
  gore_aft_410->Fit("pol0");
  gore_aft_410->GetFunction("pol0")->SetLineColor(kBlue);
  gore_aft_410->Draw("P");
  beforeafter_gore[3]->Update();
  ofile->cd();
  beforeafter_gore[3]->Write();
  beforeafter_gore[3]->Close();

  beforeafter_gore[4]->cd();
  gore_bef_420 = new TGraph(nmultgorepoints, befaft_gore_time, aging_420_gore_bef_err);
  gore_bef_420->SetName("tg_gore_bef_420");
  gore_bef_420->SetMarkerStyle(24);
  gore_bef_420->SetTitle("Standard Deviation of Repeated GoreDRP Measurements");
  gore_bef_420->GetXaxis()->SetTitle("Time, 420nm (Fraction of Years)");
  gore_bef_420->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  gore_bef_420->Fit("pol0");
  gore_bef_420->Draw("AP");
  gore_aft_420 = new TGraph(nmultgorepoints, befaft_gore_time, aging_420_gore_aft_err);
  gore_aft_420->SetName("tg_gore_aft_420");
  gore_aft_420->SetMarkerStyle(20);
  gore_aft_420->Fit("pol0");
  gore_aft_420->GetFunction("pol0")->SetLineColor(kBlue);
  gore_aft_420->Draw("P");
  beforeafter_gore[4]->Update();
  ofile->cd();
  beforeafter_gore[4]->Write();
  beforeafter_gore[4]->Close();

  beforeafter_gore[5]->cd();
  gore_bef_430 = new TGraph(nmultgorepoints, befaft_gore_time, aging_430_gore_bef_err);
  gore_bef_430->SetName("tg_gore_bef_430");
  gore_bef_430->SetMarkerStyle(24);
  gore_bef_430->SetTitle("Standard Deviation of Repeated GoreDRP Measurements");
  gore_bef_430->GetXaxis()->SetTitle("Time, 430nm (Fraction of Years)");
  gore_bef_430->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  gore_bef_430->Fit("pol0");
  gore_bef_430->Draw("AP");
  gore_aft_430 = new TGraph(nmultgorepoints, befaft_gore_time, aging_430_gore_aft_err);
  gore_aft_430->SetName("tg_gore_aft_430");
  gore_aft_430->SetMarkerStyle(20);
  gore_aft_430->Fit("pol0");
  gore_aft_430->GetFunction("pol0")->SetLineColor(kBlue);
  gore_aft_430->Draw("P");
  beforeafter_gore[5]->Update();
  ofile->cd();
  beforeafter_gore[5]->Write();
  beforeafter_gore[5]->Close();

  beforeafter_gore[6]->cd();
  gore_bef_440 = new TGraph(nmultgorepoints, befaft_gore_time, aging_440_gore_bef_err);
  gore_bef_440->SetName("tg_gore_bef_440");
  gore_bef_440->SetMarkerStyle(24);
  gore_bef_440->SetTitle("Standard Deviation of Repeated GoreDRP Measurements");
  gore_bef_440->GetXaxis()->SetTitle("Time, 440nm (Fraction of Years)");
  gore_bef_440->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  gore_bef_440->Fit("pol0");
  gore_bef_440->Draw("AP");
  gore_aft_440 = new TGraph(nmultgorepoints, befaft_gore_time, aging_440_gore_aft_err);
  gore_aft_440->SetName("tg_gore_aft_440");
  gore_aft_440->SetMarkerStyle(20);
  gore_aft_440->Fit("pol0");
  gore_aft_440->GetFunction("pol0")->SetLineColor(kBlue);
  gore_aft_440->Draw("P");
  beforeafter_gore[6]->Update();
  ofile->cd();
  beforeafter_gore[6]->Write();
  beforeafter_gore[6]->Close();

  beforeafter_gore[7]->cd();
  gore_bef_450 = new TGraph(nmultgorepoints, befaft_gore_time, aging_450_gore_bef_err);
  gore_bef_450->SetName("tg_gore_bef_450");
  gore_bef_450->SetMarkerStyle(24);
  gore_bef_450->SetTitle("Standard Deviation of Repeated GoreDRP Measurements");
  gore_bef_450->GetXaxis()->SetTitle("Time, 450nm (Fraction of Years)");
  gore_bef_450->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  gore_bef_450->Fit("pol0");
  gore_bef_450->Draw("AP");
  gore_aft_450 = new TGraph(nmultgorepoints, befaft_gore_time, aging_450_gore_aft_err);
  gore_aft_450->SetName("tg_gore_aft_450");
  gore_aft_450->SetMarkerStyle(20);
  gore_aft_450->Fit("pol0");
  gore_aft_450->GetFunction("pol0")->SetLineColor(kBlue);
  gore_aft_450->Draw("P");
  beforeafter_gore[7]->Update();
  ofile->cd();
  beforeafter_gore[7]->Write();
  beforeafter_gore[7]->Close();

  beforeafter_gore[8]->cd();
  gore_bef_600 = new TGraph(nmultgorepoints, befaft_gore_time, aging_600_gore_bef_err);
  gore_bef_600->SetName("tg_gore_bef_600");
  gore_bef_600->SetMarkerStyle(24);
  gore_bef_600->SetTitle("Standard Deviation of Repeated GoreDRP Measurements");
  gore_bef_600->GetXaxis()->SetTitle("Time, 600nm (Fraction of Years)");
  gore_bef_600->GetYaxis()->SetTitle("Standard Deviation of 3 Repeated Measurements");
  gore_bef_600->Fit("pol0");
  gore_bef_600->Draw("AP");
  gore_aft_600 = new TGraph(nmultgorepoints, befaft_gore_time, aging_600_gore_aft_err);
  gore_aft_600->SetName("tg_gore_aft_600");
  gore_aft_600->SetMarkerStyle(20);
  gore_aft_600->Fit("pol0");
  gore_aft_600->GetFunction("pol0")->SetLineColor(kBlue);
  gore_aft_600->Draw("P");
  beforeafter_gore[8]->Update();
  ofile->cd();
  beforeafter_gore[8]->Write();
  beforeafter_gore[8]->Close();

  // Close the output file and input file
  ofile->Close();
  file->Close();
}
