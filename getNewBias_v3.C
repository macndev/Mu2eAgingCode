void getNewBias_v3(){

  ifstream input("new127files.txt");
  string oneline;

  ofstream output("bad_bias_spills_new127.txt");

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
    std::cout << "padded strings are: " << runstr << "_" << subrunstr << std::endl;
    const char* runchar = runstr.Data();
    const char* subrunchar = subrunstr.Data();

    //Create output file. If it already exists, recreate it.
    TFile *ofile = new TFile(Form("output_FEBbias_cuts_v3_run%s_%s.root", runchar, subrunchar),"RECREATE");

    gStyle->SetOptStat(1111111);

    TTree *runtree = file->Get<TTree>("run");
    TTree *spilltree = file->Get<TTree>("spills");
    TTree *rstree = file->Get<TTree>("runSummary");

    const int nFEB  = 6;
    const int nChan = 64;

    int spill_num; 
    Long64_t timestamp;
    int spill_boardStatus[nFEB][22];

    spilltree->SetBranchAddress("spill_boardStatus", spill_boardStatus);
    spilltree->SetBranchAddress("spill_num", &spill_num);
    rstree->SetBranchAddress("timestamp", &timestamp);

    ofile->cd();
    TTree *otree = new TTree("otree","New Run Format Data");

    Int_t spillNumber;
    Float_t biasBus[nFEB][8];

    otree->Branch("spillNumber",&spillNumber,"spillNumber/I");
    otree->Branch("biasBus",&biasBus,"biasBus[6][8]/F");

    int this_spill = 0;
    int old_spill = 0;
    int last_spill_bad = 0;
    int nSpills = spilltree->GetEntries();
    std::vector<int> badEntries;

    for(int iSpill = 0; iSpill < nSpills; iSpill++){
      spilltree->GetEntry(iSpill);
      this_spill = spill_num;
      if(iSpill == 0){
	rstree->GetEntry(0);
	Long64_t seconds_since_1970 = timestamp;
	Long64_t seconds_in_5hrs = 5 * 60 * 60;
	seconds_since_1970 += seconds_in_5hrs;
	TTimeStamp time(seconds_since_1970);
	std::string timestring;
	timestring = time.AsString("l");
	//std::cout << "raw timestamp is: " << seconds_since_1970 << " from timestamp: " << timestamp << std::endl;
	std::cout << "timestring is: " << timestring << std::endl;
	TNamed* timename = new TNamed("start_time", timestring);
	timename->Write();

	old_spill = spill_num;   // set old_spill to first spill number on first pass through loop
      }
    
      if((iSpill == 0) || (last_spill_bad == 1)){  // special cases
	spillNumber = spill_num;
	for(int iFEB = 0; iFEB < nFEB; iFEB++){
	  for(int iBus = 0; iBus < 8; iBus++){
	    int iBus_entry = iBus + 11;
	    if((spill_boardStatus[iFEB][iBus_entry]/50.) > 0.05){
	      biasBus[iFEB][iBus] = spill_boardStatus[iFEB][iBus_entry]/50.;
	    }
	  }
	}
      } 
      if((iSpill > 0) && (this_spill == (old_spill + 1))){
	spillNumber = spill_num;
	for(int iFEB = 0; iFEB < nFEB; iFEB++){
	  for(int iBus = 0; iBus < 8; iBus++){
	    int iBus_entry = iBus + 11;
	    if((spill_boardStatus[iFEB][iBus_entry]/50.) > 0.05){
	      biasBus[iFEB][iBus] = spill_boardStatus[iFEB][iBus_entry]/50.;
	    }
	  }
	}
      }
      else if ((iSpill != 0) && (last_spill_bad != 1)) {
	last_spill_bad = 1;           // say this spill was bad for next comparison
	badEntries.push_back(iSpill);
	std::cout << "bad entry number: " << iSpill << std::endl;
      } 
      else if ((iSpill == 0) || (last_spill_bad == 1)) {
	last_spill_bad = 0; // we don't see consecutive bad entries
      }

      ofile->cd();
      otree->Fill();

      if((iSpill % 100) == 0) std::cout << "through spill number " << iSpill << std::endl; // progress print line

      old_spill = this_spill; // increment old_spill for next comparison
    }
    otree->Write();

    ofile->Write();
    ofile->Close();
    file->Close();

    output << "bad entries for run " << run_num << "_" << subrun_num << ": \n";
    for(int i : badEntries) output << i << "\n";
  }
  output.close();
}
