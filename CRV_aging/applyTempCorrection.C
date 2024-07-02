float tempCorrection(float PE, float temp){
  float PE_corr;
  float ref_temp = 25.0;

  PE_corr = PE + 0.7506 * (temp - ref_temp);
  
  return PE_corr;
}

void applyTempCorrection(){

  TFile *file = TFile::Open("/pnfs/mu2e/tape/phy-rec/rec/mu2e/CRV_wideband_cosmics/crvaging-001/root/ce/35/rec.mu2e.CRV_wideband_cosmics.crvaging-001.001029_000.root");

  TIter keyList(file->GetListOfKeys());
  TKey* key;
  TString key_name;
  while((key = (TKey*)keyList())){
    key_name = key->GetName();
    std::cout << key_name << std::endl;
  }
  TTree *tree = file->Get<TTree>(key_name);
  int run_num;
  std::cout << "key length = " << key_name.Length() << std::endl;
  TString short_str = key_name(4, key_name.Length());
  run_num = short_str.Atoi();
  std::cout << "run number: " << run_num << std::endl;

  //Create output file. If it already exists, recreate it.
  TFile *ofile = new TFile("output_tempcorrected_1029.root","RECREATE");

  //From minimizing chi2 between 3 channels we measured at test beam:
  // PE_Tcorr = -0.7506 * temp + 62.29
  //Correct everything to 25 degrees... low PE @ 30 should become higher @ 25

  float PE_corrected;
  float PE_0;
  float temp_0;

  float_t PE_tree;
  float_t temp_tree;
  tree->SetBranchAddress("PE_tree",&PE_tree);
  tree->SetBranchAddress("temp_tree",&temp_tree);

  TCanvas *page[256];
  TH1F *hPEs[128];
  TH1F *htemp[128];

  string name;
  const char* cname;
  string tname;
  const char* tcname;
  int ithHist = 0;
  for (int iFEB = 0; iFEB < 2; iFEB++){
    for(int i=0; i<64; i++){
      hPEs[i] = new TH1F(Form("hist_%d",i+(iFEB*64)),Form("hist_%d",i+(iFEB*64)), 100, 10, 120);
      name = "can" + to_string(i + iFEB*64);
      cname = name.c_str();
      page[i] = new TCanvas(cname, cname, 900, 600);
      page[i]->cd();
      tree->Draw(Form("PEs[%d][%d]>>%s",iFEB, i, hPEs[i]->GetName()) , "", "", 1000000, 0 ); // plot only 1M events
      hPEs[i]->Draw();
      ++ithHist;

      htemp[i] = new TH1F(Form("temp_%d",i+(iFEB*64)),Form("temp_%d",i+(iFEB*64)), 100, 10, 120);
      tname = "tcan" + to_string(i + iFEB*64);
      tcname = tname.c_str();
      page[i+128] = new TCanvas(tcname, tcname, 900, 600);
      page[i+128]->cd();
      tree->Draw(Form("temperature[%d][%d]>>%s",iFEB, i, htemp[i]->GetName()) , "", "", 1000000, 0 ); // plot only 1M events
      htemp[i]->Draw();

      ofile->cd();
      page[i]->Write();
      page[i]->Close();
      page[i+128]->Write();
      page[i+128]->Close();
      
      int PE_ent = hPEs[i]->GetEntries();
      int temp_ent = htemp[i]->GetEntries();

      std::cout << "entries in PE tree = " << PE_ent << ", entries in temp tree = " << temp_ent << std::endl;

      if(i == 43){
	for(int i_ent = 0; i_ent < PE_ent; i_ent++){
	  tree->GetEntry(i);
	  std::cout << "PEs: " << PE_tree << std::endl;
	  std::cout << "temp : " << temp_tree << std::endl;
	  PE_corrected = tempCorrection(PE_0, temp_0);

	  std::cout << "for tree 43, temp corrections = " << PE_corrected << std::endl;
	}
      }
    }
  }
  

}
