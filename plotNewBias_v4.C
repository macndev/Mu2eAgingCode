const char* short_timechar;

void plotNewBias_v4(){

  ifstream input("notetempfiles_bias.txt");
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
    //std::cout << "strings: " << short_runstr << "_" << short_subrunstr << std::endl;
    const char* runchar = runstr.Data();
    const char* subrunchar = subrunstr.Data();

    //Create output file. If it already exists, recreate it.
    TFile *ofile = new TFile(Form("output_FEBbias_storeflagtree_vs_spill_run%s_%s.root", runchar, subrunchar),"RECREATE");

    gStyle->SetOptStat(1111111);

    TTree *tree = file->Get<TTree>("otree");
    TTree *badtree = file->Get<TTree>("badtree");
    
    const int nFEB  = 6;

    //get info from input files
    int spillNumber;
    float biasBus[nFEB][8];
    float bad_frac;

    tree->SetBranchAddress("spillNumber", &spillNumber);
    tree->SetBranchAddress("biasBus", &biasBus);
    badtree->SetBranchAddress("bad_frac", &bad_frac);

    ofile->cd();
    TTree *obadtree = new TTree("obadtree","Fraction of Bad Events per SR");

    Float_t obad_frac;
    obadtree->Branch("obad_frac", &obad_frac, "obad_frac/F");

    TGraph* tgbiasBus[nFEB][8];
    TGraph* tgbiasBus_hours[nFEB][8];

    int nEntries = tree->GetEntries();
    double spill_num[10000];
    double hours[10000];
    double biasBus_feb0_0[10000];
    double biasBus_feb0_1[10000];
    double biasBus_feb0_2[10000];
    double biasBus_feb0_3[10000];
    double biasBus_feb0_4[10000];
    double biasBus_feb0_5[10000];
    double biasBus_feb0_6[10000];
    double biasBus_feb0_7[10000];
    double biasBus_feb1_0[10000];
    double biasBus_feb1_1[10000];
    double biasBus_feb1_2[10000];
    double biasBus_feb1_3[10000];
    double biasBus_feb1_4[10000];
    double biasBus_feb1_5[10000];
    double biasBus_feb1_6[10000];
    double biasBus_feb1_7[10000];
    double biasBus_feb2_0[10000];
    double biasBus_feb2_1[10000];
    double biasBus_feb2_2[10000];
    double biasBus_feb2_3[10000];
    double biasBus_feb2_4[10000];
    double biasBus_feb2_5[10000];
    double biasBus_feb2_6[10000];
    double biasBus_feb2_7[10000];
    double biasBus_feb3_0[10000];
    double biasBus_feb3_1[10000];
    double biasBus_feb3_2[10000];
    double biasBus_feb3_3[10000];
    double biasBus_feb3_4[10000];
    double biasBus_feb3_5[10000];
    double biasBus_feb3_6[10000];
    double biasBus_feb3_7[10000];
    double biasBus_feb4_0[10000];
    double biasBus_feb4_1[10000];
    double biasBus_feb4_2[10000];
    double biasBus_feb4_3[10000];
    double biasBus_feb4_4[10000];
    double biasBus_feb4_5[10000];
    double biasBus_feb4_6[10000];
    double biasBus_feb4_7[10000];
    double biasBus_feb5_0[10000];
    double biasBus_feb5_1[10000];
    double biasBus_feb5_2[10000];
    double biasBus_feb5_3[10000];
    double biasBus_feb5_4[10000];
    double biasBus_feb5_5[10000];
    double biasBus_feb5_6[10000];
    double biasBus_feb5_7[10000];

    // get the timestamp of the run
    TNamed *timestamp = (TNamed*)file->Get("start_time");
    const char *timechar = timestamp->GetTitle();
    std::cout << "start time is: " << timechar << std::endl;
    std::string timestring = (std::string)timechar;
    std::string short_timestring(timestring.substr(0,25));
    //std::cout << "string is: " << timestring << std::endl;
    //std::cout << " shortened to " << short_timestring << std::endl;
    short_timechar = short_timestring.c_str();
    std::cout << " shortened to " << short_timechar << std::endl;

    for(int iFEB = 0; iFEB < nFEB; iFEB++){
      for(int iEntry = 0; iEntry < nEntries; iEntry++){
	tree->GetEntry(iEntry);

	spill_num[iEntry] = spillNumber;
	hours[iEntry] = spillNumber * 0.05; // conversion to hours, 1 spill = 3 minutes = 0.05 hours
      
	if(iFEB == 0){
	  biasBus_feb0_0[iEntry] = biasBus[iFEB][0];
	  biasBus_feb0_1[iEntry] = biasBus[iFEB][1];
	  biasBus_feb0_2[iEntry] = biasBus[iFEB][2];
	  biasBus_feb0_3[iEntry] = biasBus[iFEB][3];
	  biasBus_feb0_4[iEntry] = biasBus[iFEB][4];
	  biasBus_feb0_5[iEntry] = biasBus[iFEB][5];
	  biasBus_feb0_6[iEntry] = biasBus[iFEB][6];
	  biasBus_feb0_7[iEntry] = biasBus[iFEB][7];
	}
	if(iFEB == 1){
	  biasBus_feb1_0[iEntry] = biasBus[iFEB][0];
	  biasBus_feb1_1[iEntry] = biasBus[iFEB][1];
	  biasBus_feb1_2[iEntry] = biasBus[iFEB][2];
	  biasBus_feb1_3[iEntry] = biasBus[iFEB][3];
	  biasBus_feb1_4[iEntry] = biasBus[iFEB][4];
	  biasBus_feb1_5[iEntry] = biasBus[iFEB][5];
	  biasBus_feb1_6[iEntry] = biasBus[iFEB][6];
	  biasBus_feb1_7[iEntry] = biasBus[iFEB][7];
	}
	if(iFEB == 2){
	  biasBus_feb2_0[iEntry] = biasBus[iFEB][0];
	  biasBus_feb2_1[iEntry] = biasBus[iFEB][1];
	  biasBus_feb2_2[iEntry] = biasBus[iFEB][2];
	  biasBus_feb2_3[iEntry] = biasBus[iFEB][3];
	  biasBus_feb2_4[iEntry] = biasBus[iFEB][4];
	  biasBus_feb2_5[iEntry] = biasBus[iFEB][5];
	  biasBus_feb2_6[iEntry] = biasBus[iFEB][6];
	  biasBus_feb2_7[iEntry] = biasBus[iFEB][7];
	}
	if(iFEB == 3){
	  biasBus_feb3_0[iEntry] = biasBus[iFEB][0];
	  biasBus_feb3_1[iEntry] = biasBus[iFEB][1];
	  biasBus_feb3_2[iEntry] = biasBus[iFEB][2];
	  biasBus_feb3_3[iEntry] = biasBus[iFEB][3];
	  biasBus_feb3_4[iEntry] = biasBus[iFEB][4];
	  biasBus_feb3_5[iEntry] = biasBus[iFEB][5];
	  biasBus_feb3_6[iEntry] = biasBus[iFEB][6];
	  biasBus_feb3_7[iEntry] = biasBus[iFEB][7];
	}
	if(iFEB == 4){
	  biasBus_feb4_0[iEntry] = biasBus[iFEB][0];
	  biasBus_feb4_1[iEntry] = biasBus[iFEB][1];
	  biasBus_feb4_2[iEntry] = biasBus[iFEB][2];
	  biasBus_feb4_3[iEntry] = biasBus[iFEB][3];
	  biasBus_feb4_4[iEntry] = biasBus[iFEB][4];
	  biasBus_feb4_5[iEntry] = biasBus[iFEB][5];
	  biasBus_feb4_6[iEntry] = biasBus[iFEB][6];
	  biasBus_feb4_7[iEntry] = biasBus[iFEB][7];
	}
	if(iFEB == 5){
	  biasBus_feb5_0[iEntry] = biasBus[iFEB][0];
	  biasBus_feb5_1[iEntry] = biasBus[iFEB][1];
	  biasBus_feb5_2[iEntry] = biasBus[iFEB][2];
	  biasBus_feb5_3[iEntry] = biasBus[iFEB][3];
	  biasBus_feb5_4[iEntry] = biasBus[iFEB][4];
	  biasBus_feb5_5[iEntry] = biasBus[iFEB][5];
	  biasBus_feb5_6[iEntry] = biasBus[iFEB][6];
	  biasBus_feb5_7[iEntry] = biasBus[iFEB][7];
	}
      }
    }
    for(int i = 0; i < 8; i++){
      tgbiasBus[0][0] = new TGraph(nEntries, spill_num, biasBus_feb0_0);
      tgbiasBus[0][1] = new TGraph(nEntries, spill_num, biasBus_feb0_1);
      tgbiasBus[0][2] = new TGraph(nEntries, spill_num, biasBus_feb0_2);
      tgbiasBus[0][3] = new TGraph(nEntries, spill_num, biasBus_feb0_3);
      tgbiasBus[0][4] = new TGraph(nEntries, spill_num, biasBus_feb0_4);
      tgbiasBus[0][5] = new TGraph(nEntries, spill_num, biasBus_feb0_5);
      tgbiasBus[0][6] = new TGraph(nEntries, spill_num, biasBus_feb0_6);
      tgbiasBus[0][7] = new TGraph(nEntries, spill_num, biasBus_feb0_7);
      tgbiasBus[0][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 0",i));
      tgbiasBus[0][i]->SetName(Form("biasBus0_%i",i));
      tgbiasBus[0][i]->GetXaxis()->SetTitle("Spill Number");
      tgbiasBus[0][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus[0][i]->Write();
      tgbiasBus_hours[0][0] = new TGraph(nEntries, hours, biasBus_feb0_0);
      tgbiasBus_hours[0][1] = new TGraph(nEntries, hours, biasBus_feb0_1);
      tgbiasBus_hours[0][2] = new TGraph(nEntries, hours, biasBus_feb0_2);
      tgbiasBus_hours[0][3] = new TGraph(nEntries, hours, biasBus_feb0_3);
      tgbiasBus_hours[0][4] = new TGraph(nEntries, hours, biasBus_feb0_4);
      tgbiasBus_hours[0][5] = new TGraph(nEntries, hours, biasBus_feb0_5);
      tgbiasBus_hours[0][6] = new TGraph(nEntries, hours, biasBus_feb0_6);
      tgbiasBus_hours[0][7] = new TGraph(nEntries, hours, biasBus_feb0_7);
      tgbiasBus_hours[0][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 0",i));
      tgbiasBus_hours[0][i]->SetName(Form("biasBus_hours0_%i",i));
      tgbiasBus_hours[0][i]->GetXaxis()->SetTitle("Hours Since Start of Run");
      tgbiasBus_hours[0][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus_hours[0][i]->Write();

      tgbiasBus[1][0] = new TGraph(nEntries, spill_num, biasBus_feb1_0);
      tgbiasBus[1][1] = new TGraph(nEntries, spill_num, biasBus_feb1_1);
      tgbiasBus[1][2] = new TGraph(nEntries, spill_num, biasBus_feb1_2);
      tgbiasBus[1][3] = new TGraph(nEntries, spill_num, biasBus_feb1_3);
      tgbiasBus[1][4] = new TGraph(nEntries, spill_num, biasBus_feb1_4);
      tgbiasBus[1][5] = new TGraph(nEntries, spill_num, biasBus_feb1_5);
      tgbiasBus[1][6] = new TGraph(nEntries, spill_num, biasBus_feb1_6);
      tgbiasBus[1][7] = new TGraph(nEntries, spill_num, biasBus_feb1_7);
      tgbiasBus[1][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 1",i));
      tgbiasBus[1][i]->SetName(Form("biasBus1_%i",i));
      tgbiasBus[1][i]->GetXaxis()->SetTitle("Spill Number");
      tgbiasBus[1][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus[1][i]->Write();
      tgbiasBus_hours[1][0] = new TGraph(nEntries, hours, biasBus_feb1_0);
      tgbiasBus_hours[1][1] = new TGraph(nEntries, hours, biasBus_feb1_1);
      tgbiasBus_hours[1][2] = new TGraph(nEntries, hours, biasBus_feb1_2);
      tgbiasBus_hours[1][3] = new TGraph(nEntries, hours, biasBus_feb1_3);
      tgbiasBus_hours[1][4] = new TGraph(nEntries, hours, biasBus_feb1_4);
      tgbiasBus_hours[1][5] = new TGraph(nEntries, hours, biasBus_feb1_5);
      tgbiasBus_hours[1][6] = new TGraph(nEntries, hours, biasBus_feb1_6);
      tgbiasBus_hours[1][7] = new TGraph(nEntries, hours, biasBus_feb1_7);
      tgbiasBus_hours[1][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 1",i));
      tgbiasBus_hours[1][i]->SetName(Form("biasBus_hours1_%i",i));
      tgbiasBus_hours[1][i]->GetXaxis()->SetTitle("Hours Since Start of Run");
      tgbiasBus_hours[1][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus_hours[1][i]->Write();

      tgbiasBus[2][0] = new TGraph(nEntries, spill_num, biasBus_feb2_0);
      tgbiasBus[2][1] = new TGraph(nEntries, spill_num, biasBus_feb2_1);
      tgbiasBus[2][2] = new TGraph(nEntries, spill_num, biasBus_feb2_2);
      tgbiasBus[2][3] = new TGraph(nEntries, spill_num, biasBus_feb2_3);
      tgbiasBus[2][4] = new TGraph(nEntries, spill_num, biasBus_feb2_4);
      tgbiasBus[2][5] = new TGraph(nEntries, spill_num, biasBus_feb2_5);
      tgbiasBus[2][6] = new TGraph(nEntries, spill_num, biasBus_feb2_6);
      tgbiasBus[2][7] = new TGraph(nEntries, spill_num, biasBus_feb2_7);
      tgbiasBus[2][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 2",i));
      tgbiasBus[2][i]->SetName(Form("biasBus2_%i",i));
      tgbiasBus[2][i]->GetXaxis()->SetTitle("Spill Number");
      tgbiasBus[2][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus[2][i]->Write();
      tgbiasBus_hours[2][0] = new TGraph(nEntries, hours, biasBus_feb2_0);
      tgbiasBus_hours[2][1] = new TGraph(nEntries, hours, biasBus_feb2_1);
      tgbiasBus_hours[2][2] = new TGraph(nEntries, hours, biasBus_feb2_2);
      tgbiasBus_hours[2][3] = new TGraph(nEntries, hours, biasBus_feb2_3);
      tgbiasBus_hours[2][4] = new TGraph(nEntries, hours, biasBus_feb2_4);
      tgbiasBus_hours[2][5] = new TGraph(nEntries, hours, biasBus_feb2_5);
      tgbiasBus_hours[2][6] = new TGraph(nEntries, hours, biasBus_feb2_6);
      tgbiasBus_hours[2][7] = new TGraph(nEntries, hours, biasBus_feb2_7);
      tgbiasBus_hours[2][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 2",i));
      tgbiasBus_hours[2][i]->SetName(Form("biasBus_hours2_%i",i));
      tgbiasBus_hours[2][i]->GetXaxis()->SetTitle("Hours Since Start of Run");
      tgbiasBus_hours[2][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus_hours[2][i]->Write();

      tgbiasBus[3][0] = new TGraph(nEntries, spill_num, biasBus_feb3_0);
      tgbiasBus[3][1] = new TGraph(nEntries, spill_num, biasBus_feb3_1);
      tgbiasBus[3][2] = new TGraph(nEntries, spill_num, biasBus_feb3_2);
      tgbiasBus[3][3] = new TGraph(nEntries, spill_num, biasBus_feb3_3);
      tgbiasBus[3][4] = new TGraph(nEntries, spill_num, biasBus_feb3_4);
      tgbiasBus[3][5] = new TGraph(nEntries, spill_num, biasBus_feb3_5);
      tgbiasBus[3][6] = new TGraph(nEntries, spill_num, biasBus_feb3_6);
      tgbiasBus[3][7] = new TGraph(nEntries, spill_num, biasBus_feb3_7);
      tgbiasBus[3][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 3",i));
      tgbiasBus[3][i]->SetName(Form("biasBus3_%i",i));
      tgbiasBus[3][i]->GetXaxis()->SetTitle("Spill Number");
      tgbiasBus[3][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus[3][i]->Write();
      tgbiasBus_hours[3][0] = new TGraph(nEntries, hours, biasBus_feb3_0);
      tgbiasBus_hours[3][1] = new TGraph(nEntries, hours, biasBus_feb3_1);
      tgbiasBus_hours[3][2] = new TGraph(nEntries, hours, biasBus_feb3_2);
      tgbiasBus_hours[3][3] = new TGraph(nEntries, hours, biasBus_feb3_3);
      tgbiasBus_hours[3][4] = new TGraph(nEntries, hours, biasBus_feb3_4);
      tgbiasBus_hours[3][5] = new TGraph(nEntries, hours, biasBus_feb3_5);
      tgbiasBus_hours[3][6] = new TGraph(nEntries, hours, biasBus_feb3_6);
      tgbiasBus_hours[3][7] = new TGraph(nEntries, hours, biasBus_feb3_7);
      tgbiasBus_hours[3][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 3",i));
      tgbiasBus_hours[3][i]->SetName(Form("biasBus_hours3_%i",i));
      tgbiasBus_hours[3][i]->GetXaxis()->SetTitle("Hours Since Start of Run");
      tgbiasBus_hours[3][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus_hours[3][i]->Write();

      tgbiasBus[4][0] = new TGraph(nEntries, spill_num, biasBus_feb4_0);
      tgbiasBus[4][1] = new TGraph(nEntries, spill_num, biasBus_feb4_1);
      tgbiasBus[4][2] = new TGraph(nEntries, spill_num, biasBus_feb4_2);
      tgbiasBus[4][3] = new TGraph(nEntries, spill_num, biasBus_feb4_3);
      tgbiasBus[4][4] = new TGraph(nEntries, spill_num, biasBus_feb4_4);
      tgbiasBus[4][5] = new TGraph(nEntries, spill_num, biasBus_feb4_5);
      tgbiasBus[4][6] = new TGraph(nEntries, spill_num, biasBus_feb4_6);
      tgbiasBus[4][7] = new TGraph(nEntries, spill_num, biasBus_feb4_7);
      tgbiasBus[4][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 4",i));
      tgbiasBus[4][i]->SetName(Form("biasBus4_%i",i));
      tgbiasBus[4][i]->GetXaxis()->SetTitle("Spill Number");
      tgbiasBus[4][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus[4][i]->Write();
      tgbiasBus_hours[4][0] = new TGraph(nEntries, hours, biasBus_feb4_0);
      tgbiasBus_hours[4][1] = new TGraph(nEntries, hours, biasBus_feb4_1);
      tgbiasBus_hours[4][2] = new TGraph(nEntries, hours, biasBus_feb4_2);
      tgbiasBus_hours[4][3] = new TGraph(nEntries, hours, biasBus_feb4_3);
      tgbiasBus_hours[4][4] = new TGraph(nEntries, hours, biasBus_feb4_4);
      tgbiasBus_hours[4][5] = new TGraph(nEntries, hours, biasBus_feb4_5);
      tgbiasBus_hours[4][6] = new TGraph(nEntries, hours, biasBus_feb4_6);
      tgbiasBus_hours[4][7] = new TGraph(nEntries, hours, biasBus_feb4_7);
      tgbiasBus_hours[4][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 4",i));
      tgbiasBus_hours[4][i]->SetName(Form("biasBus_hours4_%i",i));
      tgbiasBus_hours[4][i]->GetXaxis()->SetTitle("Hours Since Start of Run");
      tgbiasBus_hours[4][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus_hours[4][i]->Write();

      tgbiasBus[5][0] = new TGraph(nEntries, spill_num, biasBus_feb5_0);
      tgbiasBus[5][1] = new TGraph(nEntries, spill_num, biasBus_feb5_1);
      tgbiasBus[5][2] = new TGraph(nEntries, spill_num, biasBus_feb5_2);
      tgbiasBus[5][3] = new TGraph(nEntries, spill_num, biasBus_feb5_3);
      tgbiasBus[5][4] = new TGraph(nEntries, spill_num, biasBus_feb5_4);
      tgbiasBus[5][5] = new TGraph(nEntries, spill_num, biasBus_feb5_5);
      tgbiasBus[5][6] = new TGraph(nEntries, spill_num, biasBus_feb5_6);
      tgbiasBus[5][7] = new TGraph(nEntries, spill_num, biasBus_feb5_7);
      tgbiasBus[5][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 5",i));
      tgbiasBus[5][i]->SetName(Form("biasBus5_%i",i));
      tgbiasBus[5][i]->GetXaxis()->SetTitle("Spill Number");
      tgbiasBus[5][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus[5][i]->Write();
      tgbiasBus_hours[5][0] = new TGraph(nEntries, hours, biasBus_feb5_0);
      tgbiasBus_hours[5][1] = new TGraph(nEntries, hours, biasBus_feb5_1);
      tgbiasBus_hours[5][2] = new TGraph(nEntries, hours, biasBus_feb5_2);
      tgbiasBus_hours[5][3] = new TGraph(nEntries, hours, biasBus_feb5_3);
      tgbiasBus_hours[5][4] = new TGraph(nEntries, hours, biasBus_feb5_4);
      tgbiasBus_hours[5][5] = new TGraph(nEntries, hours, biasBus_feb5_5);
      tgbiasBus_hours[5][6] = new TGraph(nEntries, hours, biasBus_feb5_6);
      tgbiasBus_hours[5][7] = new TGraph(nEntries, hours, biasBus_feb5_7);
      tgbiasBus_hours[5][i]->SetTitle(Form("Bias Voltage for Bus %i on FEB 5",i));
      tgbiasBus_hours[5][i]->SetName(Form("biasBus_hours5_%i",i));
      tgbiasBus_hours[5][i]->GetXaxis()->SetTitle("Hours Since Start of Run");
      tgbiasBus_hours[5][i]->GetYaxis()->SetTitle("Bias Voltage (V)");
      tgbiasBus_hours[5][i]->Write();
    }					    
    badtree->GetEntry(0);
    obad_frac = bad_frac;
    obadtree->Fill();
    obadtree->Write();

    ofile->Close();
    file->Close();
  }
}
