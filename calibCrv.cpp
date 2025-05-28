# include <cstdio>
# include <cmath>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <string>
# include <cassert>
# include <vector>
# include <algorithm>
# include <experimental/filesystem>

# include "TCanvas.h"
# include "TROOT.h"
# include "TGraphErrors.h"
# include "TGraph.h"
# include "TH1.h"
# include "TH1F.h"
# include "TF1.h"
# include "TLegend.h"
# include "TLatex.h"
# include "TStyle.h"
# include "TApplication.h"
# include "TMultiGraph.h"
# include "TMath.h"
# include "TTree.h"
# include "TNtuple.h"
# include "TFile.h"
# include "TLine.h"
# include "TPaveText.h"
# include "TFitResult.h"
# include "TMarker.h"
# include "TSystem.h"

const double CRV_TDC_RATE    = 159.324e6;  // Hz
const double RATE=(CRV_TDC_RATE/2.0)/1.0e9; // GHZ
const int BOARD_STATUS_REGISTERS=22;
const int FPGA_BLOCK_REGISTERS=38;
const int FPGA_BLOCKS=4;

struct TemperatureCorrections
{
  double calibOvervoltageChange{125.5};   //ADC*ns/V
  double calibTemperatureChangeAFE{-1.46};  //ADC*ns/K  additional calib change due to temperature change in AFE
  double overvoltageTemperatureChangeCMB{-0.0554};  //V/K
  double overvoltageTemperatureChangeFEB{0.00409};  //V/K
  double referenceTemperatureCMB{20.0};   //degC
  double referenceTemperatureFEB{40.0};   //degC
};

class CrvEvent
{
  public:
  CrvEvent(const int numberOfFebs, const int channelsPerFeb, const int numberOfSamples, TTree *tree, const int preSignalSamples,
           const double noiseThreshold, const double minCalibMaxBin, const TemperatureCorrections &temperatureCorrections);
  void     FillPedestalHistograms(int entry);
  void     CalculatePedestal(double pedestalCorrection);
  void     FillCalibrationHistograms(int entry);
  void     CalculateCalibrationConstants(const std::string &pdfFileName, int nPEpeaksToFit);
  void     StorePedestalAndCalibrationConstants(std::ofstream &calibFile);
  void     Summarize(const std::string &pdfFileName);
  void     StoreHistFile(const std::string &histFileName, const std::string &runNumber);

  private:
  CrvEvent();

  int _numberOfPreSignalSamples;
  int _numberOfFebs;
  int _channelsPerFeb;
  int _numberOfSamples;

  TTree *_tree;

  double _noiseThreshold;
  double _minCalibMaxBin;

  const TemperatureCorrections _tc;

  int     _spillIndex;
  int    *_lastSpillIndex;
  int     _eventNumber;
  int    *_tdcSinceSpill;
  double *_timeSinceSpill;
  short  *_adc;
  float  *_temperature;
  int    *_boardStatus;
  int    *_FPGABlocks;

  std::vector<TCanvas*> _canvas;
  std::vector<bool>     _draw1PE;
  std::vector<bool>     _draw2PE;
  std::vector<TH1F*>    _pedestalHist;
  std::vector<TGraph*>  _temperatureHist;
  std::vector<TGraph*>  _temperatureFEBHist;
  std::vector<double>   _pedestals;
  std::vector<bool>     _deadChannels;
  std::vector<bool>     _noPedestal;
  std::vector<int>      _nPreSignalRegions;
  std::vector<int>      _nNoiseHits;

  //need two instances: for non-temperature corrected and temperature corrected calibrations
  struct calibStruct
  {
    TH1F*   _calibrationHist;
    bool    _noCalibration;
    size_t  _nPEpeaks;
    double  _calibrationConstant;
    double  _noiseRate;
    double  _xtalkProbability;
  };
  std::vector<calibStruct> _calibVector[2];

};
CrvEvent::CrvEvent(const int numberOfFebs, const int channelsPerFeb, const int numberOfSamples, TTree *tree, const int preSignalSamples,
                   const double noiseThreshold, const double minCalibMaxBin, const TemperatureCorrections &temperatureCorrections) :
                   _numberOfPreSignalSamples(preSignalSamples), _numberOfFebs(numberOfFebs), _channelsPerFeb(channelsPerFeb), _numberOfSamples(numberOfSamples),
                   _tree(tree), _noiseThreshold(noiseThreshold), _minCalibMaxBin(minCalibMaxBin), _tc(temperatureCorrections)
{
  _lastSpillIndex = new int[_numberOfFebs*_channelsPerFeb];
  _adc            = new short[_numberOfFebs*_channelsPerFeb*_numberOfSamples];
//  _timeSinceSpill = new double[_numberOfFebs];  //OLD
  _timeSinceSpill = new double[_numberOfFebs*_channelsPerFeb];
  _temperature    = new float[_numberOfFebs*_channelsPerFeb];
  _boardStatus    = new int[_numberOfFebs*BOARD_STATUS_REGISTERS];
  _FPGABlocks     = new int[_numberOfFebs*FPGA_BLOCKS*FPGA_BLOCK_REGISTERS];

  tree->SetBranchAddress("runtree_spill_index", &_spillIndex);
  tree->SetBranchAddress("runtree_adc", _adc);
  tree->SetBranchAddress("runtree_time_since_spill", _timeSinceSpill);
  tree->SetBranchAddress("runtree_temperature", _temperature);
  tree->SetBranchAddress("runtree_boardStatus", _boardStatus);
  tree->SetBranchAddress("runtree_FPGABlocks", _FPGABlocks);

  _canvas.resize(_numberOfFebs*_channelsPerFeb);
  _draw1PE.resize(_numberOfFebs*_channelsPerFeb);
  _draw2PE.resize(_numberOfFebs*_channelsPerFeb);
  _pedestalHist.resize(_numberOfFebs*_channelsPerFeb);
  _temperatureHist.resize(_numberOfFebs*_channelsPerFeb);
  _temperatureFEBHist.resize(_numberOfFebs*_channelsPerFeb);
  _pedestals.resize(_numberOfFebs*_channelsPerFeb);
  _deadChannels.resize(_numberOfFebs*_channelsPerFeb);
  _noPedestal.resize(_numberOfFebs*_channelsPerFeb);
  _nPreSignalRegions.resize(_numberOfFebs*_channelsPerFeb);
  _nNoiseHits.resize(_numberOfFebs*_channelsPerFeb);
  _calibVector[0].resize(_numberOfFebs*_channelsPerFeb);
  _calibVector[1].resize(_numberOfFebs*_channelsPerFeb);

  if(_tree->GetEntries()>0) _tree->GetEntry(0);

  for(int i=0; i<_numberOfFebs; i++)
  {
    for(int j=0; j<_channelsPerFeb; j++)
    {
      int index=i*_channelsPerFeb+j;  //used for _variable[i][j]

      _lastSpillIndex[index]=-1;
      _draw1PE[index]=true;
      _draw2PE[index]=true;
      _canvas[index]=new TCanvas();
      _canvas[index]->Divide(2,2);
      _canvas[index]->cd(1);
      gPad->Divide(2,1);  //pedestal and sample pulses
      gPad->cd(2);
      gPad->SetLogy(0);
      gPad->Divide(1,2);  //1PE and 2PE sample
      _canvas[index]->cd(2);
      gPad->Divide(2,1);  //temperatures and fits
      gPad->cd(1);
      gPad->Divide(1,2);  //temperature and FEB temperature
      _canvas[index]->cd(2);
      gPad->cd(2);
      gPad->SetLogy(0);
      gPad->Divide(1,2);  //calib fit and temperature corrected calib fit

      _pedestalHist[index]=new TH1F(Form("FEB%i Chan%i ped", i, j), Form("FEB%i Chan%i Pedestal; ADC; Counts", i, j), 1001, -50.05, 50.05);

      _temperatureHist[index]=new TGraph();
      _temperatureHist[index]->SetMarkerStyle(20);
      _temperatureHist[index]->SetMarkerSize(0.5);
      _temperatureHist[index]->SetMarkerColor(kBlack);
      _temperatureHist[index]->GetYaxis()->SetLabelSize(0.06);
      _temperatureHist[index]->SetNameTitle(Form("hT_%i",index),Form("Temperature FEB%i Chan%i;Spill;Temperature [deg C]", i, j));

      _temperatureFEBHist[index]=new TGraph();
      _temperatureFEBHist[index]->SetMarkerStyle(20);
      _temperatureFEBHist[index]->SetMarkerSize(0.5);
      _temperatureFEBHist[index]->SetMarkerColor(kBlack);
      _temperatureFEBHist[index]->GetYaxis()->SetLabelSize(0.06);
      _temperatureFEBHist[index]->SetNameTitle(Form("hTFEB_%i",index),Form("FEB Temperature FEB%i Chan%i;Spill;Temperature [deg C]", i, j));

      _calibVector[0][index]._calibrationHist=new TH1F(Form("FEB%i Chan%i cal0", i, j), Form("FEB%i Chan%i Calibration; Pulse area [ADC*ns]; Counts", i, j), 300, 0, 3000);
      _calibVector[1][index]._calibrationHist=new TH1F(Form("FEB%i Chan%i cal1", i, j), Form("FEB%i Chan%i Temp. corrected calibration; Pulse area [ADC*ns]; Counts", i, j), 300, 0, 3000);
      for(int k=0; k<2; ++k)
      {
        _calibVector[k][index]._calibrationHist->GetXaxis()->SetTitleSize(0.04);
        _calibVector[k][index]._calibrationHist->GetXaxis()->SetLabelSize(0.04);
        _calibVector[k][index]._calibrationHist->GetYaxis()->SetTitleSize(0.04);
        _calibVector[k][index]._calibrationHist->GetYaxis()->SetLabelSize(0.04);
      }
    }
  }
}
void CrvEvent::FillPedestalHistograms(int entry)
{
  _tree->GetEntry(entry);

if(entry%1000==0) std::cout<<"P "<<entry<<std::endl;

  for(int i=0; i<_numberOfFebs; i++)
  {
//    if(std::isnan(_timeSinceSpill[i])) continue;  //for missing FEBs  //OLD
    for(int j=0; j<_channelsPerFeb; j++)
    {
      int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
      if(std::isnan(_timeSinceSpill[index])) continue;  //missing FEB/channel in raw data

      const short *data = &(_adc[index*_numberOfSamples]);

      if(data[0]==0 && data[1]==0 && data[3]==0) continue; //FIXME temporary check for bad events
                                                           //where other channels work, so that timSinceSpill wasn't marked as NAN

      //divide prespill region into three parts to remove portions of the waveform that may contain pulses
      int numberOfRegions=3;
      for(int i=0; i<numberOfRegions; i++)
      {
        double average=0;
        double minADC=NAN;
        double maxADC=NAN;
        for(int j=0; j<_numberOfPreSignalSamples/numberOfRegions; j++)
        {
          short ADC=data[i*_numberOfPreSignalSamples/numberOfRegions+j];
          average+=ADC;
          if(ADC<minADC || std::isnan(minADC)) minADC=ADC;
          if(ADC>maxADC || std::isnan(maxADC)) maxADC=ADC;
        }
        average/=_numberOfPreSignalSamples/numberOfRegions;
        if(maxADC-average>2.5) continue;
        if(minADC-average<-2.5) continue;

        for(int j=0; j<_numberOfPreSignalSamples/numberOfRegions; j++)
        {
          short ADC=data[i*_numberOfPreSignalSamples/numberOfRegions+j];
          _pedestalHist[index]->Fill(ADC);
        }
      }
    }
  }
}
void CrvEvent::CalculatePedestal(double pedestalCorrection)
{
  for(int i=0; i<_numberOfFebs; i++)
  {
    for(int j=0; j<_channelsPerFeb; j++)
    {
      int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
      _canvas[index]->cd(1);
      gPad->cd(1);
      _pedestalHist[index]->SetLineColor(kBlack);
      _pedestalHist[index]->DrawClone();

      _pedestals[index]=0;
      if(_pedestalHist[index]->GetEntries()<200) continue;
      int n=_pedestalHist[index]->GetNbinsX();
      double overflow=_pedestalHist[index]->GetBinContent(0)+_pedestalHist[index]->GetBinContent(n+1);
      if(overflow/((double)_pedestalHist[index]->GetEntries())>0.1) continue;

      int maxbinPedestal = _pedestalHist[index]->GetMaximumBin();
      double peakPedestal = _pedestalHist[index]->GetBinCenter(maxbinPedestal);
      TF1 funcPedestal("f0", "gaus",peakPedestal-4,peakPedestal+4);
      funcPedestal.SetLineWidth(1);
      funcPedestal.SetLineColor(kRed);
      funcPedestal.DrawClone("SAME");
      _pedestalHist[index]->Fit(&funcPedestal, "QR");
      _pedestals[index] = funcPedestal.GetParameter(1);
      _pedestals[index]+= pedestalCorrection;

      double noiseStdDev=_pedestalHist[index]->GetStdDev();

      TPaveText t1(.15, .7, .50, .8, "NDC");
      t1.SetFillColor(0);
      t1.AddText(Form("Pedestal = %6.2lf", _pedestals[index]));
      t1.AddText(Form("Noise StdDev = %6.2lf", noiseStdDev));
      t1.SetTextAlign(12);
      t1.DrawClone("SAME");
    }
  }
}
void CrvEvent::FillCalibrationHistograms(int entry)
{
  _tree->GetEntry(entry);

if(entry%1000==0) std::cout<<"C "<<entry<<std::endl;

  for(int i=0; i<_numberOfFebs; i++)
  {
//    if(std::isnan(_timeSinceSpill[i])) continue;  //for missing FEBs  //OLD
    for(int j=0; j<_channelsPerFeb; j++)
    {
      int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
      if(std::isnan(_timeSinceSpill[index])) continue;  //missing FEB/channel in raw data

      const short *data = &(_adc[index*_numberOfSamples]);

      if(data[0]==0 && data[1]==0 && data[3]==0) continue; //FIXME temporary check for bad events
                                                           //where other channels work, so that timSinceSpill wasn't marked as NAN

      const double &pedestal = _pedestals[index];
      _canvas[index]->cd(1);
      gPad->cd(2); //pad where sample pulses are drawn

      int nBins = _numberOfPreSignalSamples;  //pre-signal region

      //remove the pedestal and find the local maxima
      std::vector<double> waveform;
      std::vector<std::pair<int,bool> > peaks;
      //double sum=0;
      for(int bin=0; bin<nBins; bin++)
      {
        waveform.push_back(data[bin]-pedestal);
        if(bin>1 && bin<nBins-3)  //don't search for peaks too close to the sample start or end
        {
          if(data[bin-1]<data[bin] && data[bin]>data[bin+1] && data[bin]-pedestal>_noiseThreshold) peaks.emplace_back(bin,false);
          if(data[bin-1]<data[bin] && data[bin]==data[bin+1] && data[bin+1]>data[bin+2] && data[bin]-pedestal>_noiseThreshold) peaks.emplace_back(bin,true);
        }
      //sum+=fabs(data[bin]-pedestal);
      }

      //don't use noisy events  //doesn't seem to do much
      //if(sum/nBins>1.0) return;  //FIXME

      _nPreSignalRegions[index]++;
      _nNoiseHits[index]+=peaks.size();

      //fit all peaks
      for(size_t iPeak=0; iPeak<peaks.size(); ++iPeak)
      {
        //select a range of up to 4 points before and after the maximum point
        //-find up to 5 points before and after the maximum point for which the waveform is stricly decreasing
        //-remove 1 point on each side. this removes potentially "bad points" belonging to a second pulse (i.e. in double pulses)
        int maxBin=peaks[iPeak].first;
        int startBin=maxBin;
        int endBin=maxBin;
        for(int bin=maxBin-1; bin>=0 && bin>=maxBin-5; bin--)
        {
          if(waveform[bin]<=waveform[bin+1]) startBin=bin;
          else break;
        }
        for(int bin=maxBin+1; bin<nBins && bin<=maxBin+5; bin++)
        {
          if(waveform[bin]<=waveform[bin-1]) endBin=bin;
          else break;
        }
        if(maxBin-startBin>1) startBin++;
        if(endBin-maxBin>1) endBin--;

        //fill the graph
        double binWidth=1.0/RATE;
        TGraph g;
        for(int bin=startBin; bin<=endBin; bin++)
        {
          double t=bin*binWidth;
          double v=waveform[bin];
          g.SetPoint(g.GetN(), t, v);
        }

        //set the fit function
        TF1 f("peakfinder","[0]*(TMath::Exp(-(x-[1])/[2]-TMath::Exp(-(x-[1])/[2])))");
        f.SetParameter(0, waveform[maxBin]*TMath::E());
        f.SetParameter(1, maxBin*binWidth);
        f.SetParameter(2, 12.6);
        if(peaks[iPeak].second) f.SetParameter(1, (maxBin+0.5)*binWidth);

        //do the fit
        TFitResultPtr fr = g.Fit(&f,(_draw1PE[index]||_draw2PE[index])?"QS":"NQS");

        if(!fr->IsValid()) continue;
        if(fr->Parameter(0)<=0 || fr->Parameter(2)<4) continue; //probably misreconstructed
        if(fr->Parameter(2)>25) continue;  //probably not a noise hit
        if(fabs(fr->Parameter(1)-maxBin*binWidth)>30) continue;
        if(fr->Parameter(0)/(waveform[maxBin]*TMath::E())>1.5) continue;

        double pulseArea = fr->Parameter(0)*fr->Parameter(2);
        _calibVector[0][index]._calibrationHist->Fill(pulseArea);

        //temperature correction of noise pulse area
        double temperatureFEB=-1000;
        if(_boardStatus[i*BOARD_STATUS_REGISTERS]!=-1) //i-th FEB was read for this spill
        {
          //temperature of i-th FEB
          temperatureFEB=_boardStatus[i*BOARD_STATUS_REGISTERS+2]*0.01;  //TODO: document seems to indicate a factor of 10.0
        }

        if(_temperature[index]>-300 && temperatureFEB>-300) //temperature of -1000 means no temperature found
        {
          //overvoltage difference for actual CMB and FEB temperatures w.r.t. reference CMB and FEB temperatures
          //deltaOvervoltage = overvoltageTemperatureChangeCMB*(TCMB-TrefCMB) + overvoltageTemperatureChangeFEB*(TFEB-TrefFEB)
          float deltaOvervoltage = _tc.overvoltageTemperatureChangeCMB*(_temperature[index]-_tc.referenceTemperatureCMB)
                                 + _tc.overvoltageTemperatureChangeFEB*(temperatureFEB-_tc.referenceTemperatureFEB);

          //calibConst(TCMB,TFEB) = calibConst(TrefCMB,TrefFEB) + calibOvervoltageChange*deltaOvervoltage(TCMB,TFEB) + calibTempChangeAFE*(TFEB-TrefFEB)
          double pulseAreaAtTrefs = pulseArea - _tc.calibOvervoltageChange*deltaOvervoltage
                                              - _tc.calibTemperatureChangeAFE*(temperatureFEB-_tc.referenceTemperatureFEB);
          _calibVector[1][index]._calibrationHist->Fill(pulseAreaAtTrefs);
        }

        if(_draw1PE[index] && waveform[maxBin]<_noiseThreshold*3.5)
        {
          TVirtualPad *p=gPad;
          gPad->cd(1);
          g.SetTitle("Example of a 1PE dark noise pulse; Time [ns]; ADC   ");
          g.SetMarkerStyle(20);
          g.SetMarkerColor(kBlack);
          g.GetHistogram()->GetXaxis()->SetTitleSize(0.05);
          g.GetHistogram()->GetXaxis()->SetLabelSize(0.05);
          g.GetHistogram()->GetYaxis()->SetTitleSize(0.05);
          g.GetHistogram()->GetYaxis()->SetLabelSize(0.05);
          g.GetHistogram()->GetYaxis()->SetTitleOffset(-0.6);
          g.DrawClone("AP");
          f.SetLineColor(kRed);
          f.DrawCopy("same");
          _draw1PE[index]=false;
          p->cd();
        }
        if(_draw2PE[index] && waveform[maxBin]>_noiseThreshold*5.0 && waveform[maxBin]<_noiseThreshold*7.0)
        {
          TVirtualPad *p=gPad;
          gPad->cd(2);
          g.SetTitle("Example of a 2PE dark noise pulse; Time [ns]; ADC   ");
          g.SetMarkerStyle(20);
          g.SetMarkerColor(kBlack);
          g.GetHistogram()->GetXaxis()->SetTitleSize(0.05);
          g.GetHistogram()->GetXaxis()->SetLabelSize(0.05);
          g.GetHistogram()->GetYaxis()->SetTitleSize(0.05);
          g.GetHistogram()->GetYaxis()->SetLabelSize(0.05);
          g.GetHistogram()->GetYaxis()->SetTitleOffset(-0.6);
          g.DrawClone("AP");
          f.SetLineColor(kRed);
          f.DrawCopy("same");
          _draw2PE[index]=false;
          p->cd();
        }
      }//peaks

      //fill temperature plot
      if(_temperature[index]>-300 && _lastSpillIndex[index]!=_spillIndex)  //temperature of -1000 means no temperature found
      {
        _temperatureHist[index]->SetPoint(_temperatureHist[index]->GetN(),_spillIndex,_temperature[index]);
        double temperatureFEB=_boardStatus[i*BOARD_STATUS_REGISTERS+2]*0.01;  //TODO: document seems to indicate a factor of 10.0
        _temperatureFEBHist[index]->SetPoint(_temperatureFEBHist[index]->GetN(),_spillIndex,temperatureFEB);
        _lastSpillIndex[index]=_spillIndex;
      }

    }//channel
  }//feb
}
void CrvEvent::CalculateCalibrationConstants(const std::string &pdfFileName, int nPEpeaksToFit)
{
//std::cout<<"REMOVED 2ND PE PEAK"<<std::endl;
  for(int i=0; i<_numberOfFebs; i++)
  {
    for(int j=0; j<_channelsPerFeb; j++)
    {
      int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
      double avgTemperature=-1000;
      double avgTemperatureFEB=-1000;

      if(_temperatureHist[index]->GetN()>0)
      {
        _canvas[index]->cd(2);
        gPad->cd(1);
        gPad->cd(1);
        _temperatureHist[index]->Draw("AP");
        avgTemperature=_temperatureHist[index]->GetMean(2);
      }

      if(_temperatureFEBHist[index]->GetN()>0)
      {
        _canvas[index]->cd(2);
        gPad->cd(1);
        gPad->cd(2);
        _temperatureFEBHist[index]->Draw("AP");
        avgTemperatureFEB=_temperatureFEBHist[index]->GetMean(2);
      }

      double minCalibMaxBinAtTref = _minCalibMaxBin;
      if(avgTemperature>-300 && avgTemperatureFEB>-300) //temperature of -1000 means no temperature found
      {
        float deltaOvervoltage = _tc.overvoltageTemperatureChangeCMB*(avgTemperature-_tc.referenceTemperatureCMB)
                               + _tc.overvoltageTemperatureChangeFEB*(avgTemperatureFEB-_tc.referenceTemperatureFEB);
        minCalibMaxBinAtTref = _minCalibMaxBin - _tc.calibOvervoltageChange*deltaOvervoltage
                                               - _tc.calibTemperatureChangeAFE*(avgTemperatureFEB-_tc.referenceTemperatureFEB);
      }

      _deadChannels[index]=false;
      _noPedestal[index]=false;

      if(std::isnan(_pedestals[index]))
      {
        _noPedestal[index]=true;
        _canvas[index]->Print(pdfFileName.c_str(), "pdf");
        continue;
      }

      if(_calibVector[0][index]._calibrationHist->GetEntries()<200)
      {
        _deadChannels[index]=true;
        _canvas[index]->Print(pdfFileName.c_str(), "pdf");
        continue;
      }

      for(int k=0; k<2; ++k) //k=0: no temperature correction, k=1: temperature correction
      {
        calibStruct &CS =_calibVector[k][index];
        CS._noCalibration=true;
        CS._nPEpeaks=0;
        CS._noiseRate=0;
        CS._xtalkProbability=0;
        CS._calibrationConstant=0;

        if(k==1 && _calibVector[k][index]._calibrationHist->GetEntries()<200)
        {
          CS._noCalibration=true;
          _canvas[index]->Print(pdfFileName.c_str(), "pdf");
          continue;
        }

        //Calibration plot
        _canvas[index]->cd(k+3);
        gPad->SetLogy();
        CS._calibrationHist->SetLineColor(kBlue);
        CS._calibrationHist->DrawClone("");

        int maxbin = 0;
        double maxbinContent = 0;
        for(int bin=1; bin<CS._calibrationHist->GetNbinsX(); bin++)
        {
          if(CS._calibrationHist->GetBinCenter(bin)<(k==0?_minCalibMaxBin:minCalibMaxBinAtTref)) continue;   //find 1PE maximum only between 250 and 1500
          if(CS._calibrationHist->GetBinCenter(bin)>1500) break;     //1500
          double binContent = CS._calibrationHist->GetBinContent(bin);
          if(binContent>maxbinContent)
          {
            maxbin=bin;
            maxbinContent=binContent;
          }
        }

        double peak1PE = CS._calibrationHist->GetBinCenter(maxbin);
        TF1 func1("f1", "gaus",peak1PE*0.8,peak1PE*1.2);
        func1.SetParameter(1,peak1PE);
        func1.SetLineWidth(1);
        func1.SetLineColor(kRed);
        CS._calibrationHist->Fit(&func1, "QR");
        peak1PE = func1.GetParameter(1);

        std::vector<double> peaks;
        if(peak1PE>200.0 && peak1PE<1800.0)
        {
          peaks.push_back(peak1PE);
          func1.DrawClone("same");
          for(int iPeak=2; iPeak<=nPEpeaksToFit; ++iPeak)
          {
            double peakTmp = iPeak*peak1PE;
            double rangeStart = peakTmp-0.35*peak1PE;
            double rangeEnd = peakTmp+0.35*peak1PE;
            if(rangeEnd>CS._calibrationHist->GetXaxis()->GetXmax()) break;
            if(CS._calibrationHist->Integral(CS._calibrationHist->FindFixBin(rangeStart),CS._calibrationHist->FindFixBin(rangeEnd))<50) break;
            TF1 func2(Form("f%i",iPeak), "gaus",rangeStart,rangeEnd);
            func2.SetParameter(1,peakTmp);
            func2.SetLineWidth(1);
            func2.SetLineColor(kRed);
            if((int)CS._calibrationHist->Fit(&func2, "QR")!=0) break;
            peakTmp = func2.GetParameter(1);
            if(peakTmp/peak1PE>iPeak*0.95 && peakTmp/peak1PE<iPeak*1.05)
            {
              peaks.push_back(peakTmp);
              func2.DrawClone("same");
            }
            else break;
          }
        }
        CS._nPEpeaks=peaks.size();

        TPaveText t2(.50, .65, .89, .89, "NDC");
        t2.SetFillColor(0);
        for(size_t iPeak=0; iPeak<peaks.size(); ++iPeak)
        {
          t2.AddText(Form("%luPE = %4.0lf ADC*ns", iPeak+1,peaks[iPeak]));
        }

        if(k==1)
        {
          t2.AddText("Temperature corrected values");
          t2.AddText(Form("reference temp %.1f deg C (CMB), %.1f deg C (FEB)",_tc.referenceTemperatureCMB, _tc.referenceTemperatureFEB));
          t2.AddText(Form("overvoltage %.5f V/K (CMB), %.5f V/K (FEB)",_tc.overvoltageTemperatureChangeCMB, _tc.overvoltageTemperatureChangeFEB));
          t2.AddText(Form("calib %.1f ADC*ns/V, calibAFE %.2f ADC*ns/K (FEB)", _tc.calibOvervoltageChange, _tc.calibTemperatureChangeAFE));
        }

        if(peaks.empty())
        {
          t2.DrawClone("SAME");
          if(k==1) _canvas[index]->Print(pdfFileName.c_str(), "pdf");
          continue;
        }

        //cross-talk and noise
        double events1PEandMore=0;
        double events2PEandMore=0;
        for(int bin=1; bin<CS._calibrationHist->GetNbinsX(); bin++)
        {
          double area = CS._calibrationHist->GetBinCenter(bin);
          double binContent = CS._calibrationHist->GetBinContent(bin);
          if(area>0.5*peak1PE) events1PEandMore+=binContent;
          if(area>1.5*peak1PE) events2PEandMore+=binContent;
        }

        CS._noiseRate=events1PEandMore/(_nPreSignalRegions[index]*_numberOfPreSignalSamples/(RATE*1e9));
        CS._xtalkProbability=events2PEandMore/events1PEandMore;
        t2.AddText(Form("Noise rate = %4.2lf MHz", CS._noiseRate/1.0e6));
        t2.AddText(Form("Xtalk prob. = %4.2lf", CS._xtalkProbability));
        t2.DrawClone("SAME");

        //Fit plot to determine calibration constants
        _canvas[index]->cd(2);  //pad for temperature and calib fits
        gPad->cd(2);  //pad for both calib fits
        gPad->cd(k+1);  //non-temperature corrected and temperature corrected calib fits
        gPad->SetLogy(0);

        TGraph *graph=new TGraph();
        graph->SetPoint(0,0,0);
        for(size_t iPeak=0; iPeak<peaks.size(); ++iPeak) graph->SetPoint(iPeak+1,iPeak+1,peaks[iPeak]);
        if(k==0) graph->SetTitle(Form("FEB%i Ch%i Calib. Fit; n PE peak; Area [ADC*ns]",i,j));
        else graph->SetTitle(Form("FEB%i Ch%i Temp. corrected calib fit; n PE peak; Area [ADC*ns]",i,j));
        graph->GetXaxis()->SetRangeUser(0, peaks.size()+0.5);
        graph->GetYaxis()->SetRangeUser(0, peaks[0]*(peaks.size()+0.5));
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kBlack);
        graph->Draw("AP");

        TF1 funcFit("calibration","[0]*x", -0.5, peaks.size()+0.5);
        funcFit.SetLineColor(kRed);
        graph->Fit(&funcFit, "QR");
        funcFit.DrawClone("SAME");
        CS._calibrationConstant = funcFit.GetParameter(0);

        TPaveText t3(.15, .65, .7, .88, "NDC");
        t3.SetLineColorAlpha(0,0);
        t3.SetFillStyle(0);
        t3.AddText("Calibration constant");
        t3.AddText(Form("%.0lf ADC*ns/PE", CS._calibrationConstant));
        t3.SetTextAlign(12);
        t3.DrawClone("SAME");

        if(k==1) _canvas[index]->Print(pdfFileName.c_str(), "pdf");
      }
    }
  }
}
void CrvEvent::StorePedestalAndCalibrationConstants(std::ofstream &calibFile)
{
  //new format / temperature corrected
  calibFile<<"calib v2  ";
  calibFile<<"reference temp: "<<_tc.referenceTemperatureCMB<<" deg C (CMB), "<<_tc.referenceTemperatureFEB<<" deg C (FEB)   ";
  calibFile<<"overvoltage: "<<_tc.overvoltageTemperatureChangeCMB<<" V/K (CMB), "<<_tc.overvoltageTemperatureChangeFEB<<" V/K (FEB)   ";
  calibFile<<"calib: "<<_tc.calibOvervoltageChange<<" ADC*ns/V, calibAFE: "<<_tc.calibTemperatureChangeAFE<<" ADC*ns/K (FEB)"<<std::endl;
  calibFile<<"  FEB  Channel  Pedestal  Calib     CalibT   Noise  Xtalk"<<std::endl;

  std::cout<<"reference temp: "<<_tc.referenceTemperatureCMB<<" deg C (CMB), "<<_tc.referenceTemperatureFEB<<" deg C (FEB)   ";
  std::cout<<"overvoltage: "<<_tc.overvoltageTemperatureChangeCMB<<" V/K (CMB), "<<_tc.overvoltageTemperatureChangeFEB<<" V/K (FEB)   ";
  std::cout<<"calib: "<<_tc.calibOvervoltageChange<<" ADC*ns/V, calibAFE: "<<_tc.calibTemperatureChangeAFE<<" ADC*ns/K (FEB)"<<std::endl;
  std::cout<<"         FEB  Channel  Pedestal  Calib     CalibT    Noise  Xtalk"<<std::endl;

  for(int i=0; i<_numberOfFebs; i++)
  {
    for(int j=0; j<_channelsPerFeb; j++)
    {
      int index=i*_channelsPerFeb+j;  //used for _variable[i][j]

      calibFile<<std::setw(4)<<i<<"  "<<std::setw(4)<<j<<"     ";
      calibFile<<std::setw(8)<<_pedestals[index]<<"  ";
      calibFile<<std::setw(8)<<_calibVector[0][index]._calibrationConstant<<"  ";
      calibFile<<std::setw(8)<<_calibVector[1][index]._calibrationConstant<<"  ";
      calibFile<<std::setw(8)<<_calibVector[1][index]._noiseRate/1.0e6<<"  ";
      calibFile<<std::setw(8)<<_calibVector[1][index]._xtalkProbability<<std::endl;

      std::cout<<"FEB/ch "<<std::setw(4)<<i<<"/"<<std::setw(4)<<j<<"     ";
      std::cout<<"pedestal "<<std::setw(8)<<_pedestals[index]<<"  ";
      std::cout<<"calibConst "<<std::setw(8)<<_calibVector[0][index]._calibrationConstant<<"  ";
      std::cout<<"calibConstT "<<std::setw(8)<<_calibVector[1][index]._calibrationConstant<<"  ";
      std::cout<<"noise "<<std::setw(8)<<_calibVector[1][index]._noiseRate/1.0e6<<"MHz"<<"  ";
      std::cout<<"xtalk "<<std::setw(8)<<_calibVector[1][index]._xtalkProbability<<std::endl;
    }
  }
}
std::string CreateSequenceString(const std::vector<int> channels)
{
  std::string sequence;
  for(size_t i=0; i<channels.size(); i++)
  {
    if(i==0) sequence.append(std::to_string(channels[i]));
    else
    {
      if(channels[i]==channels[i-1]+1) //channel is in a sequence
      {
        if(i==channels.size()-1) sequence.append("-"+std::to_string(channels[i]));  //last channel
        else
        {
          if(channels[i]!=channels[i+1]-1) sequence.append("-"+std::to_string(channels[i])); //last channel in sequence
        }
      }
      else sequence.append(","+std::to_string(channels[i]));
    }
  }
  return sequence;
}
void CrvEvent::Summarize(const std::string &pdfFileName)
{
  std::ifstream settingsFile;
  settingsFile.open("config.txt");
  if(!settingsFile.is_open()) {std::cerr<<"Could not open config.txt."<<std::endl; exit(1);}

  double minCalibConst, maxCalibConst, minNoiseRate;
  std::string settingsKey, settingsValue;
  while(settingsFile>>settingsKey>>settingsValue)
  {
    if(settingsKey=="calibMinConst")       minCalibConst=atof(settingsValue.c_str());
    if(settingsKey=="calibMaxConst")       maxCalibConst=atof(settingsValue.c_str());
    if(settingsKey=="calibMinNoise")       minNoiseRate=atof(settingsValue.c_str());
  }
  settingsFile.close();

  TCanvas summaryCanvas;
  summaryCanvas.Divide(2,4);

  for(int k=0; k<2; ++k)
  {
    summaryCanvas.cd(1+4*k);
    TH1F *hPedestal = new TH1F(Form("hPedestal%i",k),"Pedestal;ADC;Number of Channels",100,-50,50);
    for(int i=0; i<_numberOfFebs*_channelsPerFeb; i++)
    {
      if(!_noPedestal[i]) hPedestal->Fill(_pedestals[i]);
    }
    hPedestal->Draw();

    summaryCanvas.cd(2+4*k);
    TH1F *hCalibrationConstant = new TH1F(Form("hCalibrationConstant%i",k),"Calibration Constant;ADC*ns/PE;Number of Channels",100,0,1000);
    if(k==1) hCalibrationConstant->SetTitle("Temperature corrected Calibration Constant;ADC*ns/PE;Number of Channels");
    for(int i=0; i<_numberOfFebs*_channelsPerFeb; i++)
    {
      if(!_noPedestal[i] && !_deadChannels[i] && !_calibVector[k][i]._noCalibration) hCalibrationConstant->Fill(_calibVector[k][i]._calibrationConstant);
    }
    hCalibrationConstant->Draw();

    summaryCanvas.cd(3+4*k);
    TPaveText *t1=new TPaveText(.0, .75, 1., 1., "NDC");
    t1->SetFillColor(kWhite);
    t1->SetTextColor(kBlack);
    t1->SetTextAlign(12);
    TText *t1Header = t1->AddText("Channels which don't have enough noise hits, probably dead channels");
    t1Header->SetTextColor(kRed);
    for(int i=0; i<_numberOfFebs; i++)
    {
      std::vector<int> channelNumbers;
      for(int j=0; j<_channelsPerFeb; j++)
      {
        int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
        if(_deadChannels[index]) channelNumbers.push_back(j);
      }
      std::string sequence=CreateSequenceString(channelNumbers);
      t1->AddText(Form("FEB %i: %s", i, sequence.c_str()));
    }
    t1->Draw();

    TPaveText *t2=new TPaveText(.0, .45, 1., .7, "NDC");
    t2->SetFillColor(kWhite);
    t2->SetTextColor(kBlack);
    t2->SetTextAlign(12);
    TText *t2Header = t2->AddText("Channels for which no pedestal was found");
    t2Header->SetTextColor(kRed);
    for(int i=0; i<_numberOfFebs; i++)
    {
      std::vector<int> channelNumbers;
      for(int j=0; j<_channelsPerFeb; j++)
      {
        int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
        if(_noPedestal[index]) channelNumbers.push_back(j);
      }
      std::string sequence=CreateSequenceString(channelNumbers);
      t2->AddText(Form("FEB %i: %s", i, sequence.c_str()));
    }
    t2->Draw();

    summaryCanvas.cd(4+4*k);
    TPaveText *t3=new TPaveText(.0, .75, 1., 1., "NDC");
    t3->SetFillColor(kWhite);
    t3->SetTextColor(kBlack);
    t3->SetTextAlign(12);
    TText *t3Header = t3->AddText("Channels for which no calibration constant was found");
    t3Header->SetTextColor(kRed);
    for(int i=0; i<_numberOfFebs; i++)
    {
      std::vector<int> channelNumbers;
      for(int j=0; j<_channelsPerFeb; j++)
      {
        int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
        if(_calibVector[k][index]._noCalibration) channelNumbers.push_back(j);
      }
      std::string sequence=CreateSequenceString(channelNumbers);
      t3->AddText(Form("FEB %i: %s", i, sequence.c_str()));
    }
    t3->Draw();

    TPaveText *t4=new TPaveText(.0, .45, 1., .7, "NDC");
    t4->SetFillColor(kWhite);
    t4->SetTextColor(kBlack);
    t4->SetTextAlign(12);
    TText *t4Header = t4->AddText("Channels for which the 2nd peak in the calibration was not found");
    t4Header->SetTextColor(kRed);
    for(int i=0; i<_numberOfFebs; i++)
    {
      std::vector<int> channelNumbers;
      for(int j=0; j<_channelsPerFeb; j++)
      {
        int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
        if(_calibVector[k][index]._nPEpeaks<2) channelNumbers.push_back(j);
      }
      std::string sequence=CreateSequenceString(channelNumbers);
      t4->AddText(Form("FEB %i: %s", i, sequence.c_str()));
    }
    t4->Draw();

    TPaveText *t5=new TPaveText(.0, .15, 1., .4, "NDC");
    t5->SetFillColor(kWhite);
    t5->SetTextColor(kBlack);
    t5->SetTextAlign(12);
    TText *t5Header = t5->AddText(Form("Channels with calibration constants outside of %.0f and %.0f",minCalibConst,maxCalibConst));
    t5Header->SetTextColor(kRed);
    for(int i=0; i<_numberOfFebs; i++)
    {
      std::vector<int> channelNumbers;
      for(int j=0; j<_channelsPerFeb; j++)
      {
        int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
        if(_calibVector[k][index]._calibrationConstant<minCalibConst || _calibVector[k][index]._calibrationConstant>maxCalibConst) channelNumbers.push_back(j);
      }
      std::string sequence=CreateSequenceString(channelNumbers);
      t5->AddText(Form("FEB %i: %s", i, sequence.c_str()));
    }
    t5->Draw();
  }

  summaryCanvas.Print(pdfFileName.c_str(), "pdf");

  TCanvas summaryCanvas2;
  summaryCanvas2.Divide(2,4);

  for(int k=0; k<2; ++k)
  {
    summaryCanvas2.cd(1+4*k);
    TH1F *hNoise=new TH1F(Form("hNoise%i",k),"Noise Rate;Noise Rate [MHz];Number of Channels",100,0,1.0);
    if(k==1) hNoise->SetTitle("Temperature corrected Noise Rate;Noise Rate [MHz];Number of Channels");
    for(int i=0; i<_numberOfFebs*_channelsPerFeb; i++)
    {
      if(!_noPedestal[i] && !_deadChannels[i] && !_calibVector[k][i]._noCalibration) hNoise->Fill(_calibVector[k][i]._noiseRate/1.0e6);
    }
    hNoise->Draw();

    summaryCanvas2.cd(2+4*k);
    TH1F *hCrossTalk=new TH1F(Form("hCrossTalk%i",k),"XTalk Probability;XTalk Probability;Number of Channels",100,0,0.25);
    if(k==1) hCrossTalk->SetTitle("Temperature corrected XTalk Probability;XTalk Probability;Number of Channels");
    for(int i=0; i<_numberOfFebs*_channelsPerFeb; i++)
    {
      if(!_noPedestal[i] && !_deadChannels[i] && !_calibVector[k][i]._noCalibration) hCrossTalk->Fill(_calibVector[k][i]._xtalkProbability);
    }
    hCrossTalk->Draw();

    summaryCanvas2.cd(3+4*k);
    TPaveText *t6=new TPaveText(.0, .75, 1., 1., "NDC");
    t6->SetFillColor(kWhite);
    t6->SetTextColor(kBlack);
    t6->SetTextAlign(12);
    TText *t6Header = t6->AddText("Channels which don't have enough noise hits, probably dead channels");
    t6Header->SetTextColor(kRed);
    for(int i=0; i<_numberOfFebs; i++)
    {
      std::vector<int> channelNumbers;
      for(int j=0; j<_channelsPerFeb; j++)
      {
        int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
        if(_deadChannels[index]) channelNumbers.push_back(j);
      }
      std::string sequence=CreateSequenceString(channelNumbers);
      t6->AddText(Form("FEB %i: %s", i, sequence.c_str()));
    }
    t6->Draw();

    TPaveText *t7=new TPaveText(.0, .45, 1., .7, "NDC");
    t7->SetFillColor(kWhite);
    t7->SetTextColor(kBlack);
    t7->SetTextAlign(12);
    TText *t7Header = t7->AddText(Form("Channels which have a noise rate above %0.1f MHz",minNoiseRate/1.0e6));
    t7Header->SetTextColor(kRed);
    for(int i=0; i<_numberOfFebs; i++)
    {
      std::vector<int> channelNumbers;
      for(int j=0; j<_channelsPerFeb; j++)
      {
        int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
        if(_calibVector[k][index]._noiseRate>minNoiseRate) channelNumbers.push_back(j);
      }
      std::string sequence=CreateSequenceString(channelNumbers);
      t7->AddText(Form("FEB %i: %s", i, sequence.c_str()));
    }
    t7->Draw();
  }

  summaryCanvas2.Print(pdfFileName.c_str(), "pdf");
}

void CrvEvent::StoreHistFile(const std::string &histFileName, const std::string &runNumber)
{
  TFile histFile(histFileName.c_str(), "RECREATE");
  if(!histFile.IsOpen()) {std::cerr<<"Could not open hist file for run "<<runNumber<<std::endl; exit(1);}

  for(int i=0; i<_numberOfFebs; ++i)
  {
    for(int j=0; j<_channelsPerFeb; ++j)
    {
      int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
      _pedestalHist[index]->Write();
      _calibVector[0][index]._calibrationHist->Write();
      _calibVector[1][index]._calibrationHist->Write();  //temperature corrected
    }
  }

  histFile.Close();
}

void process(const std::string &runNumber, const std::string &inFileName, const std::string &calibFileName, const std::string &pdfFileName, const std::string &histFileName,
             int nPEpeaksToFit, int preSignalSamples, double pedestalCorrection, double noiseThreshold, double minCalibMaxBin, const TemperatureCorrections &tc)
{
  TFile file(inFileName.c_str(), "READ");
  if(!file.IsOpen()) {std::cerr<<"Could not read CRV file for run "<<runNumber<<std::endl; exit(1);}

  std::ofstream calibFile;
  calibFile.open(calibFileName.c_str(),std::ios_base::trunc);

  gSystem->Unlink(pdfFileName.c_str());
  TCanvas c0;
  c0.Print(Form("%s[", pdfFileName.c_str()), "pdf");

  TTree *tree = (TTree*)file.Get("run");
  TTree *treeSpills = (TTree*)file.Get("spills");
  //try older tree names
  if(tree==NULL) tree = (TTree*)file.Get(Form("run%04i",atoi(runNumber.c_str())));
  if(treeSpills==NULL) treeSpills = (TTree*)file.Get(Form("run%04i_spills", atoi(runNumber.c_str())));
  if(tree==NULL || treeSpills==NULL) {std::cerr<<"Could not find tree or spill tree"<<std::endl; exit(1);}

  int numberOfFebs;
  int channelsPerFeb;
  int numberOfSamples;
  treeSpills->SetBranchAddress("spill_number_of_febs", &numberOfFebs);
  treeSpills->SetBranchAddress("spill_channels_per_feb", &channelsPerFeb);
  treeSpills->SetBranchAddress("spill_number_of_samples", &numberOfSamples);
  treeSpills->GetEntry(0);  //to read the channelsPerFeb and numberOfSamples

  CrvEvent event(numberOfFebs, channelsPerFeb, numberOfSamples, tree, preSignalSamples, noiseThreshold, minCalibMaxBin, tc);

  int nEvents = tree->GetEntries();
//std::cout<<"USING A WRONG NUMBER OF EVENTS"<<std::endl;
//if(nEvents>100000) nEvents=2000;   //FIXME
  for(int i=0; i<nEvents; i++) event.FillPedestalHistograms(i);
  event.CalculatePedestal(pedestalCorrection);
  for(int i=0; i<nEvents; i++) event.FillCalibrationHistograms(i);

  event.CalculateCalibrationConstants(pdfFileName, nPEpeaksToFit);
  event.StorePedestalAndCalibrationConstants(calibFile);
  event.Summarize(pdfFileName);
  event.StoreHistFile(histFileName, runNumber);

  c0.Print(Form("%s]", pdfFileName.c_str()), "pdf");

  calibFile.close();
  file.Close();
}

void makeFileNames(const std::string &runNumber, std::string &inFileName, std::string &calibFileName, std::string &pdfFileName, std::string &histFileName)
{
  std::ifstream dirFile;
  dirFile.open("config.txt");
  if(!dirFile.is_open()) {std::cerr<<"Could not open config.txt."<<std::endl; exit(1);}

  std::string inDirName, outDirName;
  std::string dirType, dir;
  while(dirFile>>dirType>>dir)
  {
    if(dirType=="crvparsed") inDirName=dir;
    if(dirType=="crvcalib") outDirName=dir;
  }
  dirFile.close();

  bool found=false;
  for(const auto& dirEntry : std::experimental::filesystem::directory_iterator(inDirName))
  {
    const std::string s = dirEntry.path().filename().string();
    const std::string s0 = "crv.parsed.";
    const std::string s1 = ".run"+runNumber+".root";
    if(s.compare(0,s0.length(),s0)!=0) continue;
    if(s.compare(s.length()-s1.length(),s1.length(),s1)!=0) continue;

    found=true;
    inFileName = dirEntry.path().string();
    calibFileName = outDirName+"crv.calib."+dirEntry.path().stem().string().substr(s0.length())+".txt";
    pdfFileName = outDirName+"log.crv.calib."+dirEntry.path().stem().string().substr(s0.length())+".pdf";
    histFileName = outDirName+"crv.calib."+dirEntry.path().stem().string().substr(s0.length())+".root";
    break;
  }

  if(!found)
  {
//try Ray's version of file names
    for(const auto& dirEntry : std::experimental::filesystem::directory_iterator(inDirName))
    {
      const std::string s = dirEntry.path().filename().string();
      const std::string s0 = "ntd.mu2e.";
      const std::string s1 = "."+runNumber+".root";
      if(s.compare(0,s0.length(),s0)!=0) continue;
      if(s.compare(s.length()-s1.length(),s1.length(),s1)!=0) continue;

      found=true;
      inFileName = dirEntry.path().string();
      calibFileName = outDirName+"cal.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".txt";
      pdfFileName = outDirName+"cal.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".pdf";
      histFileName = outDirName+"cal.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".root";
      break;
    }
  }


  if(!found) {std::cerr<<"Could not open input file for run "<<runNumber<<"."<<std::endl; exit(1);}
}

void printHelp()
{
  std::cout<<"Usa as"<<std::endl;
  std::cout<<"calibCrv -h                   Prints this help."<<std::endl;
  std::cout<<"calibCrv RUNNUMBER [OPTIONS]  Calibrates a run."<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Options:"<<std::endl;
  std::cout<<"-p ###           number of PE peaks fitted"<<std::endl;
  std::cout<<"                 default: 1"<<std::endl;
  std::cout<<"-c ###           length of presignal region for the calibration"<<std::endl;
  std::cout<<"                 default: 60"<<std::endl;
  std::cout<<"-a ###           correction to the pedestal"<<std::endl;
  std::cout<<"                 default: -0.13"<<std::endl;
  std::cout<<"--noiseThreshold ###              peaks below this threshold (above pedestal) are ignored"<<std::endl;
  std::cout<<"                                  default: 5.0 ADC (use 10.0 ADC for high gain)"<<std::endl;
  std::cout<<"--minCalibMaxBin ###              minimum bin in calib histogram where the maximum bin for the 1PE peak is searched"<<std::endl;
  std::cout<<"                                  default: 250.0 ADC*ns (use 500.0 ADC*ns for gain setting 0x000)"<<std::endl;
  std::cout<<"--calibOvervoltageChange ###      change of calibration constant per overvoltage change"<<std::endl;
  std::cout<<"                                  default: 125.5 ADC*ns/V (use 412.9 ADC*ns/V for gain setting 0x000)"<<std::endl;
  std::cout<<"--calibTempChangeAFE ###          additional change of calibration constant due to temperature change in AFE"<<std::endl;
  std::cout<<"                                  default: -1.46 ADC*ns/K (use -4.80 ADC*ns/K for gain setting 0x000)"<<std::endl;
  std::cout<<"--overvoltageTempChangeCMB ###    overvoltage change due to temperature change at the CMB"<<std::endl;
  std::cout<<"                                  default: -0.0554 V/K"<<std::endl;
  std::cout<<"--overvoltageTempChangeFEB ###    overvoltage change due to temperature change at the FEB"<<std::endl;
  std::cout<<"                                  default: +0.00409 V/K"<<std::endl;
  std::cout<<"--referenceTempCMB ###            CMB temperature to which the overvoltage is adjusted to"<<std::endl;
  std::cout<<"                                  default: 20.0 degC"<<std::endl;
  std::cout<<"--referenceTempFEB ###            FEB temperature to which the overvoltage is adjusted to"<<std::endl;
  std::cout<<"                                  default: 40.0 degC"<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Note: There needs to be a file config.txt at the current location,"<<std::endl;
  std::cout<<"which lists the location of the raw files, parsed files, etc."<<std::endl;
  exit(0);
}

int main(int argc, char **argv)
{
  if(argc==1) printHelp();
  else if(strcmp(argv[1],"-h")==0) printHelp();

  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::string runNumber=argv[1];
  std::string inFileName;
  std::string calibFileName;
  std::string pdfFileName;
  std::string histFileName;
  makeFileNames(runNumber, inFileName, calibFileName, pdfFileName, histFileName);

  int    nPEpeaksToFit      =   1;
  int    preSignalSamples   =  60;
  double pedestalCorrection =  -0.13;
  double noiseThreshold     =   5.0;  //ADC  (use 10.0 for high gain)
  double minCalibMaxBin     = 250.0;  //ADC*ns  (use 500.0 for high gain)
  TemperatureCorrections tc;  //default values are set in the definition of this struct
  for(int i=2; i<argc-1; i++)
  {
    if(strcmp(argv[i],"-p")==0) nPEpeaksToFit=atoi(argv[i+1]);
    if(strcmp(argv[i],"-c")==0) preSignalSamples=atoi(argv[i+1]);
    if(strcmp(argv[i],"-a")==0) pedestalCorrection=atof(argv[i+1]);
    if(strcmp(argv[i],"--noiseThreshold")==0)            noiseThreshold=atof(argv[i+1]);
    if(strcmp(argv[i],"--minCalibMaxBin")==0)            minCalibMaxBin=atof(argv[i+1]);
    if(strcmp(argv[i],"--calibOvervoltageChange")==0)    tc.calibOvervoltageChange=atof(argv[i+1]);
    if(strcmp(argv[i],"--calibTempChangeAFE")==0)        tc.calibTemperatureChangeAFE=atof(argv[i+1]);
    if(strcmp(argv[i],"--overvoltageTempChangeCMB")==0)  tc.overvoltageTemperatureChangeCMB=atof(argv[i+1]);
    if(strcmp(argv[i],"--overvoltageTempChangeFEB")==0)  tc.overvoltageTemperatureChangeFEB=atof(argv[i+1]);
    if(strcmp(argv[i],"--referenceTempCMB")==0)          tc.referenceTemperatureCMB=atof(argv[i+1]);
    if(strcmp(argv[i],"--referenceTempFEB")==0)          tc.referenceTemperatureFEB=atof(argv[i+1]);
  }

  process(runNumber, inFileName, calibFileName, pdfFileName, histFileName, nPEpeaksToFit, preSignalSamples,
          pedestalCorrection, noiseThreshold, minCalibMaxBin, tc);

  return 0;
}
