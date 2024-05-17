void PE_stability(int conf_num, int run_num){

  TFile *file = TFile::Open(Form("/pnfs/mu2e/scratch/users/macndev/wideband/newreco/rec.mu2e.CRV_wideband_cosmics.crvaging-00%d.00%d_000.root", conf_num, run_num));

  TIter keyList(file->GetListOfKeys());
  TKey* key;
  TString key_name;
  while((key = (TKey*)keyList())){
    key_name = key->GetName();
    std::cout << key_name << std::endl;
  }
  TTree *tree = file->Get<TTree>(key_name);
  // int run_num;
  // std::cout << "key length = " << key_name.Length() << std::endl;
  // TString short_str = key_name(4, key_name.Length());
  // run_num = short_str.Atoi();
  // std::cout << "run number: " << run_num << std::endl;

  //Create output file. If it already exists, recreate it.
  TFile *ofile = new TFile(Form("PE_stability_plots_%d.root", run_num),"RECREATE");

  TCanvas *page[256];
  TH2F* hPE_vs_spill[256];

  string name;
  const char* cname;
  for (int iFEB = 0; iFEB < 4; iFEB++){
    for(int i=0; i<64; i++){
      //hPE_vs_spill[i] = new TH2F(Form("hist_%d",i+(iFEB*64)),Form("hist_%d",i+(iFEB*64)), 500, 0., 500., 120, 0., 120.);
      name = "can" + to_string(i + iFEB*64);
      cname = name.c_str();
      page[i] = new TCanvas(cname, cname, 900, 600);
      page[i]->cd();
      tree->Draw(Form("PEs[%d][%d]:spillNumber>>htemp2",iFEB, i) , Form("PEs[%d][%d]>10",iFEB, i), "prof", 1000000, 0 ); // plot only 1M events
      hPE_vs_spill[i] = (TH2F*)gPad->GetPrimitive("htemp2"); // 2D
      hPE_vs_spill[i]->SetNameTitle(Form("hist_%d",i+(iFEB*64)),Form("hist_%d",i+(iFEB*64)));

      tree->Draw(Form("PEs[%d][%d]:spillNumber",iFEB,i), Form("PEs[%d][%d]>10",iFEB, i), "prof", 1000000, 0 ); // plot only 1M events
      tree->Draw(Form("temperature[%d][%d]:spillNumber",iFEB, i) , Form("PEs[%d][%d]>10",iFEB, i), "same prof", 1000000, 0 ); // plot only 1M events
      
      float avg_PE = hPE_vs_spill[i]->GetMean(2);
      std::cout << "Mean Y = " << avg_PE << std::endl;
      ofile->cd();
      hPE_vs_spill[i]->Write();
      page[i]->Write();
      page[i]->Close();
    }
  }
  
    //run1030->Draw("PEs[0][1]:spillNumber", "PEs[0][1]>10", "prof");
    //run1030->Draw("temperature[0][1]:spillNumber", "PEs[0][1]>10", "same prof");

  //final version = remove 1M event limit in Draw lines

}
