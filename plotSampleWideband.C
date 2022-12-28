{
  const int CHANNEL_PER_FEB=64;
//  const int NUMBER_OF_FEBS=4;
  const int NUMBER_OF_FEBS=2;
  const int NUMBER_OF_SAMPLES=128-1;

  const std::string runNumber="000105_007";
  const int sampleStart=0;
  const int sampleEnd=126;
  const int eventNumberStart=0;
  const int eventNumberEnd=2;
  const int channelStart=30;
  const int channelEnd=40;
  const int FEB=0;

  TFile *file = TFile::Open(Form("/pnfs/mu2e/scratch/outstage/ehrlich/wideband7/crvparsed/ntd.mu2e.CRV_wideband_cosmics.crvaging-001.%s.root",runNumber.c_str()));
  TTree *tree = (TTree*)file->FindObjectAny("run");

  short adc[NUMBER_OF_FEBS][CHANNEL_PER_FEB][NUMBER_OF_SAMPLES];
  tree->SetBranchAddress("runtree_adc", adc);

  for(int eventNumber=eventNumberStart; eventNumber<=eventNumberEnd; eventNumber++)
  {
    if(eventNumber>=tree->GetEntries()) continue;
    tree->GetEntry(eventNumber);

    std::string title=Form("run%s_event%i",runNumber.c_str(),eventNumber);
    TCanvas *c = new TCanvas(title.c_str(),title.c_str(),800,600);
    TMultiGraph *mg = new TMultiGraph();
    TGraph *g[NUMBER_OF_FEBS*CHANNEL_PER_FEB];
    for(int channel=channelStart; channel<=channelEnd; channel++)
    {
      int nSamples=sampleEnd-sampleStart+1;
      int index=FEB*CHANNEL_PER_FEB+channel;
      g[index]=new TGraph(nSamples);
      for(int sample=sampleStart; sample<=sampleEnd; sample++)
      {
        double t=sample*12.55;
        double v=adc[FEB][channel][sample];
        g[index]->SetPoint(sample-sampleStart, t, v);
      }

      g[index]->SetTitle(title.c_str());
      g[index]->SetMarkerStyle(20);
      g[index]->SetMarkerColor(channel-channelStart+1);
      g[index]->GetHistogram()->GetXaxis()->SetTitle("Time [ns]");
      g[index]->GetHistogram()->GetYaxis()->SetTitle("ADC");
      g[index]->GetHistogram()->GetXaxis()->SetTitleOffset(1.1);
      g[index]->GetHistogram()->GetYaxis()->SetTitleOffset(1.1);
      mg->Add(g[index]);
    }
    mg->SetTitle((title+";t[ns]").c_str());
    mg->Draw("AP");
    stringstream filename;
    filename<<"Run"<<runNumber<<"_Event"<<eventNumber<<".png";
//    c->SaveAs(filename.str().c_str());
  }

}
