void getNewTemps_v3(){

  ifstream input("notetemps.txt");
  string oneline;

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
    TFile *ofile = new TFile(Form("output_FEBCMBtemp_v3_run%s_%s.root", runchar, subrunchar),"RECREATE");

    gStyle->SetOptStat(1111111);

    TTree *runtree = file->Get<TTree>("run");
    TTree *spilltree = file->Get<TTree>("spills");
    TTree *rstree = file->Get<TTree>("runSummary");

    const int nFEB  = 6;
    const int nChan = 64;

    int spill_nevents, spill_neventsActual, spill_num; 
    Long64_t timestamp;
    int spill_boardStatus[nFEB][22];
    float temperature[nFEB][nChan];

    spilltree->SetBranchAddress("spill_nevents", &spill_nevents);
    spilltree->SetBranchAddress("spill_neventsActual", &spill_neventsActual);
    spilltree->SetBranchAddress("spill_boardStatus", spill_boardStatus);
    spilltree->SetBranchAddress("spill_num", &spill_num);
    runtree->SetBranchAddress("temperature", temperature);
    rstree->SetBranchAddress("timestamp", &timestamp);

    ofile->cd();
    TTree *otree = new TTree("otree","New Run Format Data");

    Int_t   spillNumber;
    Float_t CMBtemp[nFEB*nChan];
    Float_t FEBtemp[nFEB];
          
    otree->Branch("spillNumber",&spillNumber,"spillNumber/I");
    otree->Branch("CMBtemp",&CMBtemp,"CMBtemp[384]/F");
    otree->Branch("FEBtemp",&FEBtemp,"FEBtemp[6]/F");

    int index = 0;
    int nSpills = spilltree->GetEntries();

    for(int iSpill = 0; iSpill < nSpills; iSpill++){
      spilltree->GetEntry(iSpill);
      spillNumber = spill_num;
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
      }

      //std::cout << "neventsActual = " << spill_neventsActual << ", index = " << index << std::endl;
      
      for(int iFEB = 0; iFEB < nFEB; iFEB++){
	FEBtemp[iFEB] = spill_boardStatus[iFEB][2]/100.;
      }

      if(spill_neventsActual == 0){
	for(int iFEB = 0; iFEB < nFEB; iFEB++){
	  for(int iChan = 0; iChan < nChan; iChan++){
	    CMBtemp[iChan+iFEB*nChan] = -999.;
	  }
	}
	continue;
      }
      if(spill_neventsActual > 0){
	runtree->GetEntry(index);
	//std::cout << "got entry from run tree with index " << index << std::endl;

	for(int iFEB = 0; iFEB < nFEB; iFEB++){
	  for(int iChan = 0; iChan < nChan; iChan++){
	    CMBtemp[iChan+iFEB*nChan] = temperature[iFEB][iChan];
	  }
	}
	index += spill_neventsActual;
      }
      ofile->cd();
      otree->Fill();

      if((iSpill % 100) == 0) std::cout << "through spill number " << iSpill << std::endl;
    }
    otree->Write();

    ofile->Write();
    ofile->Close();
    file->Close();
  }
}
