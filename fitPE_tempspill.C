
#include <TMath.h>
TH1F* chi2hist = new TH1F("chi2hist","chi2hist",500,0.,500.);
TH2F* PE_temp[16];
float fitFunc(TH1F *h,  int ithHist){
  float mpv = h->GetBinCenter(h->GetMaximumBin());
  TF1 *g1 = new TF1("g1", "gaus", mpv-20, mpv+20);
  TF1 *g2 = new TF1("g2", "landau", mpv+20, 120);
  TF1 *total = new TF1("total", "gaus(0)+landau(3)+pol2(6)", mpv-20, mpv+60);

  h->Fit(g1, "R");
  h->Fit(g2, "R+");

  const int npars = 9;
  double par[npars];
  g1->GetParameters(&par[0]);
  g2->GetParameters(&par[3]);

  for(int i=0; i<npars; i++)
  total->SetParLimits(i, 0, 1e4);

  total->SetParameters(par);

  //h->Smooth(1);
  //  h->Fit(total, "R");
   Int_t fitStatus = h->Fit(total,"R");
 std:cout << "ithHist convergence = " << ithHist << " " << fitStatus << std::endl;
  total->GetParameters(&par[0]);
  // extract chi2
  TF1* fit = h->GetFunction("total");
  Float_t chi2 = fit->GetChisquare();
  if (chi2 > 200){std::cout << "ithhist, chi2 = " << ithHist << " "  << chi2 << std::endl;}
  chi2hist->Fill(chi2);
  // now total is the final fit, get peak value of total
  int maxval  = total->GetMaximum(20,100);
  float maxvalx = total->GetMaximumX(20,100);
  std::cout << "ithhist = " << ithHist << " " << "maxval = " << maxval << ", maxvalx = " << maxvalx << std::endl;
  ++ithHist;
  return maxvalx;
}


void fitPE_tempspill(){

  TFile *file = TFile::Open("/pnfs/mu2e/scratch/users/ehrlich/testbeam_temperature/crvreco/crv.reco.testbeam2022.run4778032.root");

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
  TFile *ofile = new TFile("output_fitPE_tempspill4778032.root","RECREATE");

  TCanvas *page[32];
  //page->Divide(4,4);
  TH1F *hPEs[16];
  TH1F *htemp[16];

  // Create a tree with branches (variables)
  TTree *newtree = new TTree("newtree","Data holder tree");
  Float_t max_xval;
  newtree->Branch("max_xval",&max_xval,"max_xval/F");

  string name;
  const char* cname;
  string tname;
  const char* tcname;
  float this_xval = -1.;
  int ithHist = 0;
  for (int iFEB = 0; iFEB < 1; iFEB++){
    for(int i=0; i<16; i++){
      hPEs[i] = new TH1F(Form("hist_%d",i+(iFEB*64)),Form("hist_%d",i+(iFEB*64)), 100, 10, 120);
      name = "can" + to_string(i + iFEB*64);
      cname = name.c_str();
      page[i] = new TCanvas(cname, cname, 900, 600);
      page[i]->cd();
      tree->Draw(Form("PEs[%d][%d]>>%s",iFEB, i, hPEs[i]->GetName()) , "", "", 1000000, 0 ); // plot only 1M events
      hPEs[i]->Draw();
      if(hPEs[i]->Integral()>1000)
	this_xval = fitFunc(hPEs[i],ithHist);
      ++ithHist;
      max_xval = this_xval;
      std::cout << "max value from return = "<< max_xval << std::endl;

      htemp[i] = new TH1F(Form("temp_%d",i+(iFEB*64)),Form("temp_%d",i+(iFEB*64)), 100, 10, 120);
      tname = "tcan" + to_string(i + iFEB*64);
      tcname = tname.c_str();
      page[i+16] = new TCanvas(tcname, tcname, 900, 600);
      page[i+16]->cd();
      tree->Draw(Form("temperature[%d][%d]>>%s",iFEB, i, htemp[i]->GetName()) , "", "", 1000000, 0 ); // plot only 1M events
      htemp[i]->Draw();
      float temp = htemp[i]->GetMean();
      if(temp > 43){
	std::cout << "temperature = " << temp << " , this_xval = " << this_xval << std::endl;
      }
      ofile->cd();
      newtree->Fill();
      page[i]->Write();
      page[i]->Close();
      page[i+16]->Write();
      page[i+16]->Close();
      PE_temp[i] = new TH2F(Form("PE_vs_temp_%d",(i + iFEB*64)), Form("PE_vs_temp_%d",(i + iFEB*64)), 150, 20., 50., 200, 30., 50.);
      PE_temp[i]->SetMarkerStyle(20);
      PE_temp[i]->GetXaxis()->SetTitle("Temperature");
      PE_temp[i]->GetYaxis()->SetTitle("PE Yield (peak value from fit)");      
      PE_temp[i]->Fill(temp, this_xval);
      PE_temp[i]->Write();
    }
  }
  chi2hist->Write();
  newtree->Write();
}
