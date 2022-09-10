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
# include "TLegendEntry.h"
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

const float CRV_TDC_RATE    = 159.324e6;  // Hz
const float RATE=(CRV_TDC_RATE/2.0)/1.0e9; // GHZ
const int BOARD_STATUS_REGISTERS=22;
const float DEFAULT_BETA=16.0;

class Calibration
{
  public:
  Calibration(const std::string &calibFileName, const int numberOfFebs, const int channelsPerFeb);
  bool  IsNewFormat() const {return _newFormat;}
  const std::vector<float> &GetPedestals() const {return _pedestal;}
  const std::vector<float> &GetCalibrationFactors() const {return _calibrationFactor;}
  const std::vector<float> &GetCalibrationFactorsTemperatureCorrected() const {return _calibrationFactorTemperatureCorrected;}

  private:
  int                 _numberOfFebs;
  int                 _channelsPerFeb;
  bool                _newFormat;  //temperature corrected
  std::vector<float>  _pedestal;
  std::vector<float>  _calibrationFactor;
  std::vector<float>  _calibrationFactorTemperatureCorrected;

  Calibration();
};

Calibration::Calibration(const std::string &calibFileName, const int numberOfFebs, const int channelsPerFeb) : 
                         _numberOfFebs(numberOfFebs), _channelsPerFeb(channelsPerFeb), _newFormat(false)
{
  std::ifstream calibFile;
  calibFile.open(calibFileName.c_str());
  if(!calibFile.is_open()) {std::cerr<<"Could not open calib file."<<std::endl; exit(1);}

  _pedestal.resize(_numberOfFebs*_channelsPerFeb);
  _calibrationFactor.resize(_numberOfFebs*_channelsPerFeb);
  _calibrationFactorTemperatureCorrected.resize(_numberOfFebs*_channelsPerFeb);
  std::string tmp;

  getline(calibFile,tmp);
  if(tmp.find("calib v2")==0)
  {
    getline(calibFile,tmp); //table header
    _newFormat=true;
    int i, j;
    float pedestalTmp, calibrationFactorTmp, calibrationFactorTemperatureCorrectedTmp;
    while(calibFile >> i >> j >> pedestalTmp >> calibrationFactorTmp >> calibrationFactorTemperatureCorrectedTmp)
    {
      _pedestal.at(i*_channelsPerFeb+j)=pedestalTmp;
      _calibrationFactor.at(i*_channelsPerFeb+j)=calibrationFactorTmp;
      _calibrationFactorTemperatureCorrected.at(i*_channelsPerFeb+j)=calibrationFactorTemperatureCorrectedTmp;
    }
    calibFile.close();
    return;
  }

  //try old format
  calibFile.seekg(0);
  for (int i=0; i<numberOfFebs; i++)
  {
    for (int j=0; j<_channelsPerFeb; j++)
    {
      calibFile>>tmp;
      _pedestal[i*_channelsPerFeb+j]=atof(tmp.c_str());
    }
  }
  for (int i=0; i<numberOfFebs; i++)
  {
    for (int j=0; j<_channelsPerFeb; j++) 
    {
      calibFile>>tmp;
      _calibrationFactor[i*_channelsPerFeb+j]=atof(tmp.c_str());
      _calibrationFactorTemperatureCorrected[i*_channelsPerFeb+j]=0.0;
    }
  }
  calibFile.close();
}

class CrvRecoEvent
{
  CrvRecoEvent();

  public:
  CrvRecoEvent(int signalRegionStart, int signalRegionEnd) : _f("peakfitter","[0]*(TMath::Exp(-(x-[1])/[2]-TMath::Exp(-(x-[1])/[2])))"), _signalRegionStart(signalRegionStart), _signalRegionEnd(signalRegionEnd) {Init();} 
  void Init();
  bool FailedFit(TFitResultPtr fr);
  void PeakFitter(const int* data, int numberOfSamples, float pedestal, float calibrationFactor, bool &draw);

  TF1    _f;
  int    _signalRegionStart;
  int    _signalRegionEnd;

  int    _fitStatus[2]; //0: no pulse, 1: everything ok, 2: fit failed
  float  _PEs[2];
  float  _pulseHeight[2];
  float  _beta[2];
  float  _time[2];
  float  _LEtime[2];
  int    _recoStartBin[2], _recoEndBin[2];
};
void CrvRecoEvent::Init()
{
  for(int i=0; i<2; ++i)
  {
    _fitStatus[i]=0; _PEs[i]=0; _pulseHeight[i]=NAN; _beta[i]=NAN; _time[i]=NAN; _LEtime[i]=NAN; _recoStartBin[i]=-1; _recoEndBin[i]=-1;
  }
}
bool CrvRecoEvent::FailedFit(TFitResultPtr fr)
{
  if(fr!=0) return true;
  if(!fr->IsValid()) return true;

  const double tolerance=0.01; //TODO: Try to ask the minimizer, if the parameter is at the limit
  for(int i=0; i<=2; ++i)
  {
    double v=fr->Parameter(i);
    double lower, upper;
    fr->ParameterBounds(i,lower,upper);
    if((v-lower)/(upper-lower)<tolerance) return true;
    if((upper-v)/(upper-lower)<tolerance) return true;
  }
  return false;
}
void CrvRecoEvent::PeakFitter(const int* data, int numberOfSamples, float pedestal, float calibrationFactor, bool &draw)
{
  if(std::isnan(calibrationFactor) || calibrationFactor==0) return;

  if(data[0]==0 && data[1]==0 && data[2]==0 && data[3]==0) return; //FIXME temporary check for bad events
                                                                   //where other channels work, so that timSinceSpill wasn't marked as NAN

  //remove the pedestal and find the maxima in the signal region
  std::vector<float> waveform;
  std::vector<std::pair<float, std::pair<int,int> > > peaks;    //std::pair<ADC, std::pair<peakbinsStart,peakbinsEnd> >
                                                                //for peaks where two more more consecutive ADC values are equal
                                                                //(e.g. if the max ADC value is reached)
                                                                //the beginning and end of these peak bins is indicated by 
                                                                //peakbinsStart and peakbinsEnd.
                                                                //the actual peak is probably in between this interval
                                                                //(relevant for the seed value of the fit)
  int peakBinsStart=0;
  int peakBinsEnd=0;
  for(int bin=_signalRegionStart; bin<=std::min(_signalRegionEnd,numberOfSamples); bin++) 
  {
    waveform.push_back(data[bin]-pedestal);
    if(bin<=_signalRegionStart) continue;

    if(data[bin-1]<data[bin])  //rising edge
    {
      peakBinsStart=bin;
      peakBinsEnd=bin;
    }
    if(data[bin-1]==data[bin])  //potentially at a peak with consecutive ADC values which are equal
    {
      peakBinsEnd=bin;
    }
    if(data[bin-1]>data[bin])  //falling edge
    {
      if(peakBinsStart>0)   //found a peak
      {
        if(data[peakBinsStart]-pedestal>6) //ignores fluctuations of the baseline
          peaks.emplace_back(data[peakBinsStart]-pedestal, std::make_pair(peakBinsStart,peakBinsEnd));
      }
      peakBinsStart=0;  //so that the loop has to wait for the next rising edge
    }
  }

  //ignore events without pulses
  if(peaks.size()==0) return;

//only look at the largest peak (and if available the peak after that for a potential reflected pulse)
  size_t iPeak=std::max_element(peaks.begin(),peaks.end())-peaks.begin();
  for(int i=0; i<2 && iPeak<peaks.size(); ++i, ++iPeak)
  {
    peakBinsStart = peaks[iPeak].second.first-_signalRegionStart;
    peakBinsEnd   = peaks[iPeak].second.second-_signalRegionStart;
    float averagePeakBin = 0.5*(peakBinsStart+peakBinsEnd);

    //select a range of up to 4 points before and after the peak
    //-find up to 5 points before and after the peak for which the waveform is stricly decreasing
    //-remove 1 point on each side. this removes potentially "bad points" belonging to a second pulse (i.e. in double pulses)

    int nBins = waveform.size();
    _recoStartBin[i]=peakBinsStart-1;
    _recoEndBin[i]=peakBinsEnd+1;
    for(int bin=peakBinsStart-1; bin>=0 && bin>=peakBinsStart-5; bin--)
    {
      if(waveform[bin]<=waveform[bin+1]) _recoStartBin[i]=bin;
      else break;
    }
    for(int bin=peakBinsEnd+1; bin<=nBins && bin<=peakBinsEnd+5; bin++)
    {
      if(waveform[bin]<=waveform[bin-1]) _recoEndBin[i]=bin;
      else break;
    }
    if(peakBinsStart-_recoStartBin[i]>1) _recoStartBin[i]++;
    if(_recoEndBin[i]-peakBinsEnd>1) _recoEndBin[i]--;

    //fill the graph
    float binWidth=1.0/RATE;
    TGraph g;
    for(int bin=_recoStartBin[i]; bin<=_recoEndBin[i]; bin++) 
    {
      float t=bin*binWidth;
      float v=waveform[bin];
      g.SetPoint(g.GetN(), t, v);
    }

    //set the fit function
    _f.SetParameter(0, waveform[peakBinsStart]*TMath::E());
    _f.SetParameter(1, averagePeakBin*binWidth);
    _f.SetParameter(2, DEFAULT_BETA);
    _f.SetParLimits(0, waveform[peakBinsStart]*TMath::E()*0.7,waveform[peakBinsStart]*TMath::E()*1.5);
    _f.SetParLimits(1, averagePeakBin*binWidth-15.0,averagePeakBin*binWidth+15.0);
    _f.SetParLimits(2, 5.0, 40.0);

    //do the fit
    TFitResultPtr fr = g.Fit(&_f,(draw && i==0)?"QS":"NQS");
    bool invalidFit=FailedFit(fr);
    if(invalidFit)
    {
      _PEs[i]         = waveform[peakBinsStart]*TMath::E()*DEFAULT_BETA/calibrationFactor; //using maximum ADC value of this pulse and a typical value of beta
      _pulseHeight[i] = waveform[peakBinsStart];
      _time[i]        = averagePeakBin*binWidth;
      _beta[i]        = DEFAULT_BETA; 
      _LEtime[i]      = _time[i]-0.985*DEFAULT_BETA;   //time-0.985*beta for 50% pulse height
      _fitStatus[i]   = 2;
      draw            = false;
    }
    else
    {
      _PEs[i]         = fr->Parameter(0)*fr->Parameter(2) / calibrationFactor;
      _pulseHeight[i] = fr->Parameter(0)/TMath::E();
      _time[i]        = fr->Parameter(1);
      _beta[i]        = fr->Parameter(2);
      _LEtime[i]      = _time[i]-0.985*_beta[i];   //at 50% of pulse height
      _fitStatus[i]   = 1;
    }

    _time[i]+=_signalRegionStart*binWidth;
    _LEtime[i]+=_signalRegionStart*binWidth;

    if(_PEs[i]<20 || i>0) draw=false;
    if(draw)
    {
      g.SetTitle(Form("Example of a %0.0f PE pulse; Time [ns]; ADC   ",_PEs[i]));
      g.SetMarkerStyle(20);
      g.SetMarkerColor(kBlack);
      g.GetHistogram()->GetXaxis()->SetTitleSize(0.05);
      g.GetHistogram()->GetXaxis()->SetLabelSize(0.05);
      g.GetHistogram()->GetYaxis()->SetTitleSize(0.05);
      g.GetHistogram()->GetYaxis()->SetLabelSize(0.05);
      g.GetHistogram()->GetYaxis()->SetTitleOffset(-0.6);
      g.DrawClone("AP");
      _f.SetLineColor(kRed);
      _f.DrawCopy("same");
    }
  }
  return;
}

class CrvEvent
{
  public:
  CrvEvent(const std::string &runNumber, const int numberOfFebs, const int channelsPerFeb, const int numberOfSamples, 
           TTree *tree, TTree *recoTree);
  void     Reconstruct(int entry, const Calibration &calib);
  TCanvas *GetCanvas(int feb, int channel) {return _canvas[feb*_channelsPerFeb+channel];}
  TH1F    *GetHistPEs(int i, int feb, int channel) 
           {
             return (i==0?_histPEs[feb*_channelsPerFeb+channel]:_histPEsTemperatureCorrected[feb*_channelsPerFeb+channel]);
           }
  TGraph  *GetHistTemperatures(int feb, int channel) {return _histTemperatures[feb*_channelsPerFeb+channel];}
  float    GetReferenceTemperature() {return _referenceTemperature;}
  float    GetPETemperatureIntercept() {return _PETemperatureIntercept;}
  float    GetCalibTemperatureIntercept() {return _calibTemperatureIntercept;}

  private:
  CrvEvent();

  int   _signalRegionStart;
  int   _signalRegionEnd;
  float _referenceTemperature;
  float _PETemperatureIntercept;
  float _calibTemperatureIntercept;

  std::string _runNumber;
  int _numberOfFebs;
  int _channelsPerFeb;
  int _numberOfSamples;
  
  TTree *_tree;
  TTree *_recoTree;

  int     _spillNumber;
  int    *_lastSpillNumber;
  int     _eventNumber;
  //int    *_tdcSinceSpill;  //OLD
  Long64_t *_tdcSinceSpill;
  double *_timeSinceSpill;
  int    *_adc;
  float  *_temperature;

  int    *_fitStatus;
  float  *_PEs;
  float  *_PEsTemperatureCorrected;
  float  *_pulseHeight;
  float  *_beta;
  float  *_time;
  float  *_LEtime;
  int    *_recoStartBin;
  int    *_recoEndBin;
  float  *_pedestal;

  int    *_fitStatusReflectedPulse;
  float  *_PEsReflectedPulse;
  float  *_PEsTemperatureCorrectedReflectedPulse;
  float  *_pulseHeightReflectedPulse;
  float  *_betaReflectedPulse;
  float  *_timeReflectedPulse;
  float  *_LEtimeReflectedPulse;
  int    *_recoStartBinReflectedPulse;
  int    *_recoEndBinReflectedPulse;

  std::vector<TCanvas*> _canvas;
  std::vector<int>      _plot;
  std::vector<TH1F*>    _histPEs, _histPEsTemperatureCorrected;
  std::vector<TGraph*>  _histTemperatures;
};
CrvEvent::CrvEvent(const std::string &runNumber, const int numberOfFebs, const int channelsPerFeb, const int numberOfSamples,
                   TTree *tree, TTree *recoTree) :
                   _runNumber(runNumber), _numberOfFebs(numberOfFebs), _channelsPerFeb(channelsPerFeb), _numberOfSamples(numberOfSamples), 
                   _tree(tree), _recoTree(recoTree)
{
  std::ifstream configFile;
  configFile.open("config.txt");
  if(!configFile.is_open()) {std::cerr<<"Could not open config.txt."<<std::endl; exit(1);}

  _PETemperatureIntercept=-1.0;
  _calibTemperatureIntercept=-1.0;
  _referenceTemperature=25.0;
  std::string configKey, configValue;
  while(configFile>>configKey>>configValue)
  {
    if(configKey=="signalRegionStart") _signalRegionStart=atoi(configValue.c_str());
    if(configKey=="signalRegionEnd")   _signalRegionEnd=atoi(configValue.c_str());
    if(configKey=="referenceTemperature")   _referenceTemperature=atof(configValue.c_str());
    if(configKey=="PETemperatureIntercept") _PETemperatureIntercept=atof(configValue.c_str());
    if(configKey=="calibTemperatureIntercept") _calibTemperatureIntercept=atof(configValue.c_str());
  }
  if(_PETemperatureIntercept==-1.0) {std::cerr<<"PETemperatureIntercept missing in config.txt"<<std::endl; exit(1);}
  if(_calibTemperatureIntercept==-1.0) {std::cerr<<"calibTemperatureIntercept missing in config.txt"<<std::endl; exit(1);}

  _lastSpillNumber= new int[_numberOfFebs*_channelsPerFeb];
//  _tdcSinceSpill  = new int[_numberOfFebs];  //OLD
//  _timeSinceSpill = new double[_numberOfFebs];  //OLD
  _tdcSinceSpill  = new Long64_t[_numberOfFebs*_channelsPerFeb];
  _timeSinceSpill = new double[_numberOfFebs*channelsPerFeb];
  _adc            = new int[_numberOfFebs*_channelsPerFeb*_numberOfSamples];
  _temperature    = new float[_numberOfFebs*_channelsPerFeb];

  tree->SetBranchAddress("runtree_spill_num", &_spillNumber);
  tree->SetBranchAddress("runtree_event_num",&_eventNumber);
  tree->SetBranchAddress("runtree_tdc_since_spill", _tdcSinceSpill);
  tree->SetBranchAddress("runtree_time_since_spill", _timeSinceSpill);
  tree->SetBranchAddress("runtree_adc", _adc);
  tree->SetBranchAddress("runtree_temperature", _temperature);

  _fitStatus               = new int[_numberOfFebs*_channelsPerFeb];
  _PEs                     = new float[_numberOfFebs*_channelsPerFeb];
  _PEsTemperatureCorrected = new float[_numberOfFebs*_channelsPerFeb];
  _pulseHeight             = new float[_numberOfFebs*_channelsPerFeb];
  _beta                    = new float[_numberOfFebs*_channelsPerFeb];
  _time                    = new float[_numberOfFebs*_channelsPerFeb];
  _LEtime                  = new float[_numberOfFebs*_channelsPerFeb];
  _recoStartBin            = new int[_numberOfFebs*_channelsPerFeb];
  _recoEndBin              = new int[_numberOfFebs*_channelsPerFeb];
  _pedestal                = new float[_numberOfFebs*_channelsPerFeb];

  _fitStatusReflectedPulse               = new int[_numberOfFebs*_channelsPerFeb];
  _PEsReflectedPulse                     = new float[_numberOfFebs*_channelsPerFeb];
  _PEsTemperatureCorrectedReflectedPulse = new float[_numberOfFebs*_channelsPerFeb];
  _pulseHeightReflectedPulse             = new float[_numberOfFebs*_channelsPerFeb];
  _betaReflectedPulse                    = new float[_numberOfFebs*_channelsPerFeb];
  _timeReflectedPulse                    = new float[_numberOfFebs*_channelsPerFeb];
  _LEtimeReflectedPulse                  = new float[_numberOfFebs*_channelsPerFeb];
  _recoStartBinReflectedPulse            = new int[_numberOfFebs*_channelsPerFeb];
  _recoEndBinReflectedPulse              = new int[_numberOfFebs*_channelsPerFeb];

  recoTree->Branch("spillNumber", &_spillNumber, "spillNumber/I");
  recoTree->Branch("eventNumber", &_eventNumber, "eventNumber/I");
  recoTree->Branch("tdcSinceSpill", _tdcSinceSpill, Form("tdcSinceSpill[%i][%i]/L",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("timeSinceSpill", _timeSinceSpill, Form("timeSinceSpill[%i][%i]/D",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("fitStatus", _fitStatus, Form("fitStatus[%i][%i]/I",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("PEs", _PEs, Form("PEs[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("PEsTemperatureCorrected", _PEsTemperatureCorrected, Form("PEsTemperatureCorrected[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("temperature", _temperature, Form("temperature[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("pulseHeight", _pulseHeight, Form("pulseHeight[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("beta", _beta, Form("beta[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("time", _time, Form("time[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("LEtime", _LEtime, Form("LEtime[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("adc", _adc, Form("adc[%i][%i][%i]/I",_numberOfFebs,_channelsPerFeb,_numberOfSamples));
  recoTree->Branch("recoStartBin", _recoStartBin, Form("recoStartBin[%i][%i]/I",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("recoEndBin", _recoEndBin, Form("recoEndBin[%i][%i]/I",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("pedestal", _pedestal, Form("pedestal[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("fitStatusReflectedPulse", _fitStatusReflectedPulse, Form("fitStatusReflectedPulse[%i][%i]/I",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("PEsReflectedPulse", _PEsReflectedPulse, Form("PEsReflectedPulse[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("PEsTemperatureCorrectedReflectedPulse", _PEsTemperatureCorrectedReflectedPulse, Form("PEsTemperatureCorrectedReflectedPulse[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("pulseHeightReflectedPulse", _pulseHeightReflectedPulse, Form("pulseHeightReflectedPulse[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("betaReflectedPulse", _betaReflectedPulse, Form("betaReflectedPulse[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("timeReflectedPulse", _timeReflectedPulse, Form("timeReflectedPulse[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("LEtimeReflectedPulse", _LEtimeReflectedPulse, Form("LEtimeReflectedPulse[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("recoStartBinReflectedPulse", _recoStartBinReflectedPulse, Form("recoStartBinReflectedPulse[%i][%i]/I",_numberOfFebs,_channelsPerFeb));
  recoTree->Branch("recoEndBinReflectedPulse", _recoEndBinReflectedPulse, Form("recoEndBinReflectedPulse[%i][%i]/I",_numberOfFebs,_channelsPerFeb));

  _canvas.resize(_numberOfFebs*_channelsPerFeb);
  _plot.resize(_numberOfFebs*_channelsPerFeb);
  _histPEs.resize(_numberOfFebs*_channelsPerFeb);
  _histPEsTemperatureCorrected.resize(_numberOfFebs*_channelsPerFeb);
  _histTemperatures.resize(_numberOfFebs*_channelsPerFeb);
  for(int i=0; i<_numberOfFebs; i++)
  {
    for(int j=0; j<_channelsPerFeb; j++)
    {
      int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
      _lastSpillNumber[index]=-1;
      _plot[index]=4;
      _canvas[index]=new TCanvas();
      _canvas[index]->Divide(2,2);
      for(int k=0; k<4; ++k)
      {
        _canvas[index]->cd(k+1);
//        gPad->SetPad(0,0,0.5,1);
      }
      _canvas[index]->cd(3);
      gPad->Divide(2,2);

      _histPEs[index]=new TH1F(Form("h0_%i",index),Form("PE Distribution Run %s FEB %i Channel %i;PE;Counts",_runNumber.c_str(),i,j),150,0,150);
      _histPEsTemperatureCorrected[index]=new TH1F(Form("h1_%i",index),Form("PE Distribution Run %s FEB %i Channel %i;PE;Counts",_runNumber.c_str(),i,j),150,0,150);
      _histTemperatures[index]=new TGraph();

      _histPEs[index]->SetLineColor(kBlack);
      _histPEsTemperatureCorrected[index]->SetLineColor(kBlack);
      _histTemperatures[index]->SetMarkerStyle(20);
      _histTemperatures[index]->SetMarkerSize(0.5);
      _histTemperatures[index]->SetMarkerColor(kBlack);
      _histTemperatures[index]->SetNameTitle(Form("hT_%i",index),Form("Temperature Run %s FEB %i Channel %i;Spill;Temperature [deg C]",_runNumber.c_str(),i,j));
    }
  }
}
void CrvEvent::Reconstruct(int entry, const Calibration &calib)
{
   CrvRecoEvent reco(_signalRegionStart,_signalRegionEnd);

  _tree->GetEntry(entry);

if(entry%1000==0) std::cout<<"R "<<entry<<std::endl;

  for(int i=0; i<_numberOfFebs; i++)
  {
    for(int j=0; j<_channelsPerFeb; j++)
    {
      reco.Init();

      int index=i*_channelsPerFeb+j;  //used for _variable[i][j]

      bool draw = (_plot[index]>0);
      if(draw)
      {
        _canvas[index]->cd(3);
        gPad->cd(_plot[index]);
      }

      float pedestal = calib.GetPedestals().at(index);
      float calibrationFactor = calib.GetCalibrationFactors().at(index);
      float calibrationFactorTemperatureCorrected = calib.GetCalibrationFactorsTemperatureCorrected().at(index);

      if(!isnan(_timeSinceSpill[index]))  //missing FEB/channel in raw data
        reco.PeakFitter(&(_adc[index*_numberOfSamples]), _numberOfSamples, pedestal, calibrationFactor, draw);

      //main pulse
      _fitStatus[index]               = reco._fitStatus[0];
      _PEs[index]                     = reco._PEs[0];
      _PEsTemperatureCorrected[index] = -1;

      //PEref = PE * (calibT0-Tref)/(calibT0-T) * (PET0-Tref)/(PET0-T)  //calib correction and PE correction
      float temperatureCorrection=(_PETemperatureIntercept-_referenceTemperature)/(_PETemperatureIntercept-_temperature[index]);
      if(calib.IsNewFormat() && calibrationFactorTemperatureCorrected!=0)
      {
        //old calib format doesn't have temperature correction
        float calibTemperatureCorrection=(_calibTemperatureIntercept-_referenceTemperature)/(_calibTemperatureIntercept-_temperature[index]);
        //use temperature-corrected calibration constant instead of the non-temperature-corrected calibration constant
        temperatureCorrection*=calibTemperatureCorrection*calibrationFactor/calibrationFactorTemperatureCorrected;
      }
      if(_temperature[index]!=0) _PEsTemperatureCorrected[index] = reco._PEs[0]*temperatureCorrection;

      _pulseHeight[index]             = reco._pulseHeight[0];
      _beta[index]                    = reco._beta[0];
      _time[index]                    = reco._time[0];
      _LEtime[index]                  = reco._LEtime[0];
      _recoStartBin[index]            = reco._recoStartBin[0]+_signalRegionStart;
      _recoEndBin[index]              = reco._recoEndBin[0]+_signalRegionStart;
      _pedestal[index]                = pedestal;

      if(draw && reco._fitStatus[0]==1) _plot[index]--;

      //reflected pulse
      _fitStatusReflectedPulse[index]               = reco._fitStatus[1];
      _PEsReflectedPulse[index]                     = reco._PEs[1];
      _PEsTemperatureCorrectedReflectedPulse[index] = -1;
      if(_temperature[index]!=0) _PEsTemperatureCorrectedReflectedPulse[index] = reco._PEs[1]*temperatureCorrection;
      _pulseHeightReflectedPulse[index]             = reco._pulseHeight[1];
      _betaReflectedPulse[index]                    = reco._beta[1];
      _timeReflectedPulse[index]                    = reco._time[1];
      _LEtimeReflectedPulse[index]                  = reco._LEtime[1];
      _recoStartBinReflectedPulse[index]            = reco._recoStartBin[1]+_signalRegionStart;
      _recoEndBinReflectedPulse[index]              = reco._recoEndBin[1]+_signalRegionStart;

      //fill histograms
      if(reco._fitStatus[0]==1)
      {
        _histPEs[index]->Fill(reco._PEs[0]);
        _histPEsTemperatureCorrected[index]->Fill(_PEsTemperatureCorrected[index]);
      }
      if(_temperature[index]!=0 && _lastSpillNumber[index]!=_spillNumber)
      {
        _histTemperatures[index]->SetPoint(_histTemperatures[index]->GetN(),_histTemperatures[index]->GetN(),_temperature[index]);
        _lastSpillNumber[index]=_spillNumber;
      }
    }
  }

  _recoTree->Fill();
}

double LandauGaussFunction(double *x, double *par) 
{
    //From $ROOTSYS/tutorials/fit/langaus.C
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation),
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.

    // Numeric constants
    Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    Double_t mpshift  = -0.22278298;       // Landau maximum location

    // Control constants
    Double_t np = 100.0;      // number of convolution steps
    Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

    // Variables
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;

    // MP shift correction
    mpc = par[1] - mpshift * par[0];

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    step = (xupp-xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) 
    {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);

      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
}
void LandauGauss(TH1F &h, float &mpv, float &fwhm, float &signals)
{
    float maxX=0;
    float maxValue=0;
    for(int i=15; i<=h.GetNbinsX(); i++)
    {
      if(h.GetBinContent(i)>maxValue) {maxValue=h.GetBinContent(i); maxX=h.GetBinCenter(i);}
    }
//    if(maxValue<50) return;
    if(maxValue<20) return;

    //Fit range
    Double_t fitRangeStart=0.7*maxX;
    Double_t fitRangeEnd  =1.5*maxX;
    if(maxValue<50)
    {
      fitRangeStart=0.5*maxX;
      fitRangeEnd  =1.7*maxX;
    }
    if(fitRangeStart<20.0) fitRangeStart=20.0;

    //Parameters 
    Double_t startValues[4], parLimitsLow[4], parLimitsHigh[4];
    //Most probable value
    startValues[1]=maxX;
    parLimitsLow[1]=fitRangeStart;
    parLimitsHigh[1]=fitRangeEnd;
    //Area
    startValues[2]=h.Integral(h.FindBin(fitRangeStart),h.FindBin(fitRangeEnd));
    parLimitsLow[2]=0.01*startValues[2];
    parLimitsHigh[2]=100*startValues[2];
    //Other parameters
    startValues[0]=4.0;   startValues[3]=8.0;
    parLimitsLow[0]=1.0;  parLimitsLow[3]=1.0;
    parLimitsHigh[0]=20.0; parLimitsHigh[3]=20.0;

    TF1 fit("LandauGauss",LandauGaussFunction,fitRangeStart,fitRangeEnd,4);
    fit.SetParameters(startValues);
    fit.SetLineColor(kRed);
    fit.SetParNames("Width","MP","Area","GSigma");
    for(int i=0; i<4; i++) fit.SetParLimits(i, parLimitsLow[i], parLimitsHigh[i]);
    h.Fit(&fit,"QR");
    fit.Draw("same");

    mpv = fit.GetMaximumX();
    float halfMaximum = fit.Eval(mpv)/2.0;
    float leftX = fit.GetX(halfMaximum,0.0,mpv);
    float rightX = fit.GetX(halfMaximum,mpv,10.0*mpv);
    fwhm = rightX-leftX;
    signals = h.Integral(15,150);
}

void BoardRegisters(TTree *treeSpills, std::ofstream &txtFile, const int numberOfFebs, int *febID, int &nSpillsActual, int *nFebSpillsActual, 
                    float *febTemperaturesAvg, float *supplyMonitorsAvg, float *biasVoltagesAvg, int *pipeline, int *samples, time_t &timestamp)
{
  if(!treeSpills->GetBranch("spill_boardStatus")) return;  //older file: board status was not stored

  int   nEventsExpected;
  int   nEventsActual;
  bool  spillStored;
  struct tm  timestampStruct;
  int  *boardStatus = new int[numberOfFebs*BOARD_STATUS_REGISTERS];
  treeSpills->SetBranchAddress("spill_nevents", &nEventsExpected);
  treeSpills->SetBranchAddress("spill_neventsActual", &nEventsActual);
  treeSpills->SetBranchAddress("spill_stored", &spillStored);
  treeSpills->SetBranchAddress("spill_boardStatus", boardStatus);
  treeSpills->SetBranchAddress("spill_timestamp_sec", &timestampStruct.tm_sec);
  treeSpills->SetBranchAddress("spill_timestamp_min", &timestampStruct.tm_min);
  treeSpills->SetBranchAddress("spill_timestamp_hour", &timestampStruct.tm_hour);
  treeSpills->SetBranchAddress("spill_timestamp_mday", &timestampStruct.tm_mday);
  treeSpills->SetBranchAddress("spill_timestamp_mon", &timestampStruct.tm_mon);
  treeSpills->SetBranchAddress("spill_timestamp_year", &timestampStruct.tm_year);
  treeSpills->SetBranchAddress("spill_timestamp_wday", &timestampStruct.tm_wday);
  treeSpills->SetBranchAddress("spill_timestamp_yday", &timestampStruct.tm_yday);
  treeSpills->SetBranchAddress("spill_timestamp_isdst", &timestampStruct.tm_isdst);

  treeSpills->GetEntry(0);
  txtFile<<"timestamp: "<<asctime(&timestampStruct);
  timestamp=mktime(&timestampStruct);

  int nEventsExpectedTotal=0;   //of the spills that were stored
  int nEventsActualTotal=0;
  int nSpillsExpected=treeSpills->GetEntries();
  nSpillsActual=0;

  bool  foundFeb[numberOfFebs]={false};
  for(int feb=0; feb<numberOfFebs; ++feb)
  {
    febID[feb]=0; 
    febTemperaturesAvg[feb]=0;
    pipeline[feb]=0;
    samples[feb]=0;
    nFebSpillsActual[feb]=0;
  }
  for(int i=0; i<numberOfFebs*8; ++i)
  {
    supplyMonitorsAvg[i]=0;
    biasVoltagesAvg[i]=0;
  }
  float febTemperatures[numberOfFebs]={0};
  float supplyMonitors[numberOfFebs][8]={0};
  float biasVoltages[numberOfFebs][8]={0};
  for(int iSpill=0; iSpill<nSpillsExpected; ++iSpill)
  {
    treeSpills->GetEntry(iSpill);

    if(spillStored)
    {
      nEventsExpectedTotal+=nEventsExpected;
      nEventsActualTotal+=nEventsActual;
      ++nSpillsActual;
    }

    for(int feb=0; feb<numberOfFebs; ++feb)
    {
      if(boardStatus[feb*BOARD_STATUS_REGISTERS]==-1) continue; //FEB was not read for this spill

      if(!foundFeb[feb])
      {
        foundFeb[feb]=true;
        febID[feb]=boardStatus[feb*BOARD_STATUS_REGISTERS];
        pipeline[feb]=boardStatus[feb*BOARD_STATUS_REGISTERS+20];
        samples[feb]=boardStatus[feb*BOARD_STATUS_REGISTERS+21];
      }

      ++nFebSpillsActual[feb];
      febTemperatures[feb]+=boardStatus[feb*BOARD_STATUS_REGISTERS+2]*0.01;  //TODO: document seems to indicate a factor of 10.0
      supplyMonitors[feb][0]+=boardStatus[feb*BOARD_STATUS_REGISTERS+3+0]*0.001;
      supplyMonitors[feb][1]+=boardStatus[feb*BOARD_STATUS_REGISTERS+3+1]*0.001;
      supplyMonitors[feb][2]+=boardStatus[feb*BOARD_STATUS_REGISTERS+3+2]*0.002;
      supplyMonitors[feb][3]+=boardStatus[feb*BOARD_STATUS_REGISTERS+3+3]*0.004;
      supplyMonitors[feb][4]+=boardStatus[feb*BOARD_STATUS_REGISTERS+3+4]*0.001;
      supplyMonitors[feb][5]+=boardStatus[feb*BOARD_STATUS_REGISTERS+3+5]*0.002;
      supplyMonitors[feb][6]+=boardStatus[feb*BOARD_STATUS_REGISTERS+3+6]*0.006;
      supplyMonitors[feb][7]+=boardStatus[feb*BOARD_STATUS_REGISTERS+3+7]*0.001;
      for(int i=0; i<8; ++i) biasVoltages[feb][i]+=boardStatus[feb*BOARD_STATUS_REGISTERS+11+i]*0.02;
    }
  }


  txtFile<<"Expected spills: "<<nSpillsExpected<<std::endl;
  txtFile<<"Stored spills:   "<<nSpillsActual<<std::endl;
  txtFile<<"Expected number of events (of all stored spills): "<<nEventsExpectedTotal<<std::endl;
  txtFile<<"Actual number of events                         : "<<nEventsActualTotal<<std::endl;

  txtFile<<std::fixed;
  txtFile<<std::setprecision(2);
  txtFile<<"FEB  ID   spills  FEBtemp  15Vmon  10Vmon  5Vmon   -5Vmon  3.3Vmon 2.5Vmon 1.8Vmon  1.2Vmon"<<std::endl;
  for(int feb=0; feb<numberOfFebs; ++feb)
  {
    if(nFebSpillsActual[feb]>0)
    {
      febTemperaturesAvg[feb]=febTemperatures[feb]/nFebSpillsActual[feb];
      txtFile<<feb<<"    "<<febID[feb]<<"   "<<nFebSpillsActual[feb]<<"     "<<febTemperaturesAvg[feb]<<"    ";
      for(int i=0; i<8; ++i) {supplyMonitorsAvg[feb*8+i]=supplyMonitors[feb][i]/nFebSpillsActual[feb]; txtFile<<supplyMonitorsAvg[feb*8+i]<<"    ";}
      txtFile<<std::endl;
    }
    else
    {
      txtFile<<feb<<"    no spills"<<std::endl;
    }
  }

  txtFile<<"FEB  bias0  bias1  bias2  bias3  bias4  bias5  bias6  bias7  pipeline  samples"<<std::endl;
  for(int feb=0; feb<numberOfFebs; ++feb)
  {
    if(nFebSpillsActual[feb]>0)
    {
      txtFile<<feb<<"    ";
      for(int i=0; i<8; ++i) {biasVoltagesAvg[feb*8+i]=biasVoltages[feb][i]/nFebSpillsActual[feb]; txtFile<<biasVoltagesAvg[feb*8+i]<<"  ";}
      txtFile<<pipeline[feb]<<"         "<<samples[feb]<<"  ";
      txtFile<<std::endl;
    }
    else
    {
      txtFile<<feb<<"    no spills"<<std::endl;
    }
  }

  txtFile<<std::endl;
}

void StorePEyields(const std::string &txtFileName, const int numberOfFebs, const int channelsPerFeb,
                   const std::vector<float> mpvs[2], const std::vector<float> fwhms[2], const std::vector<float> signals[2],
                   const std::vector<float> &meanTemperatures, TTree *treeSpills, int nEventsActual, const Calibration &calib)
{
  std::ifstream settingsFile;
  settingsFile.open("config.txt");
  if(!settingsFile.is_open()) {std::cerr<<"Could not open config.txt."<<std::endl; exit(1);}

  float referenceTemperature, PETemperatureIntercept, calibTemperatureIntercept;
  std::string settingsKey, settingsValue;
  while(settingsFile>>settingsKey>>settingsValue)
  {
    if(settingsKey=="referenceTemperature")   referenceTemperature=atof(settingsValue.c_str());
    if(settingsKey=="PETemperatureIntercept") PETemperatureIntercept=atof(settingsValue.c_str());
    if(settingsKey=="calibTemperatureIntercept") calibTemperatureIntercept=atof(settingsValue.c_str());
  }
  settingsFile.close();

  time_t timestamp;  //=long
  int   *febID = new int[numberOfFebs];
  int    nSpillsActual;
  int   *nFebSpillsActual = new int[numberOfFebs];
  float *febTemperaturesAvg = new float[numberOfFebs];
  float *supplyMonitorsAvg = new float[numberOfFebs*8];
  float *biasVoltagesAvg = new float[numberOfFebs*8];
  int   *pipeline = new int[numberOfFebs];
  int   *samples = new int[numberOfFebs];
  float *mpvsSummary = const_cast<float*>(mpvs[0].data());
  float *mpvsTSummary = const_cast<float*>(mpvs[1].data());
  float *fwhmsSummary = const_cast<float*>(fwhms[0].data());
  float *fwhmsTSummary = const_cast<float*>(fwhms[1].data());
  float *signalsSummary = const_cast<float*>(signals[0].data());
  float *signalsTSummary = const_cast<float*>(signals[1].data());
  float *meanTemperaturesSummary = const_cast<float*>(meanTemperatures.data());
  float *pedestals = const_cast<float*>(calib.GetPedestals().data());
  float *calibConstants = const_cast<float*>(calib.GetCalibrationFactors().data());
  float *calibConstantsT = const_cast<float*>(calib.GetCalibrationFactorsTemperatureCorrected().data());
  TTree *recoTreeSummary = new TTree("runSummary","runSummary");
  recoTreeSummary->Branch("timestamp", &timestamp, "timestamp/L");
  recoTreeSummary->Branch("febID", febID, Form("febID[%i]/I",numberOfFebs));
  recoTreeSummary->Branch("spillsRecorded", &nSpillsActual, "spillsRecorded/I");
  recoTreeSummary->Branch("eventsRecorded", &nEventsActual, "eventsRecorded/I");
  recoTreeSummary->Branch("febSpills", nFebSpillsActual, Form("febSpills[%i]/I",numberOfFebs));
  recoTreeSummary->Branch("febTemperaturesAvg", febTemperaturesAvg, Form("febTemperaturesAvg[%i]/F",numberOfFebs));
  recoTreeSummary->Branch("supplyMonitorsAvg", supplyMonitorsAvg, Form("supplyMonitorsAvg[%i][%i]/F",numberOfFebs,8));
  recoTreeSummary->Branch("biasVoltagesAvg", biasVoltagesAvg, Form("biasVoltagesAvg[%i][%i]/F",numberOfFebs,8));
  recoTreeSummary->Branch("pipeline", pipeline, Form("pipeline[%i]/I",numberOfFebs));
  recoTreeSummary->Branch("samples", samples, Form("samples[%i]/I",numberOfFebs));
  recoTreeSummary->Branch("PEs", mpvsSummary, Form("PEs[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("PEsTemperatureCorrected", mpvsTSummary, Form("PEsTemperatureCorrected[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("FWHMs", fwhmsSummary, Form("FWHMs[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("FWHMsTemperatureCorrected", fwhmsTSummary, Form("FWHMsTemperatureCorrected[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("signals", signalsSummary, Form("signals[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("signalsTemperatureCorrected", signalsTSummary, Form("signalsTemperatureCorrected[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("meanTemperatures", meanTemperaturesSummary, Form("meanTemperatures[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("pedestals", pedestals, Form("pedestals[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("calibConstants", calibConstants, Form("calibConstants[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("calibConstantsTemperatureCorrected", calibConstantsT, Form("calibConstantsTemperatureCorrected[%i][%i]/F",numberOfFebs,channelsPerFeb));

  std::ofstream txtFile;
  txtFile.open(txtFileName.c_str());
 
  timestamp=0;
  nSpillsActual=0; 
  BoardRegisters(treeSpills, txtFile, numberOfFebs, febID, nSpillsActual, nFebSpillsActual, febTemperaturesAvg, supplyMonitorsAvg, biasVoltagesAvg, pipeline, samples, timestamp);

  txtFile<<"referenceTemperature: "<<referenceTemperature<<" deg C  ";
  txtFile<<"calibTemperatureIntercept: "<<calibTemperatureIntercept<<" deg C  ";
  txtFile<<"PETemperatureIntercept: "<<PETemperatureIntercept<<" deg C"<<std::endl;
  txtFile<<"  FEB  Channel     PE     PE Temp Corr    FWHM                   Signals              meanTemp    pedestal  calibConst  calibConst Temp Corr"<<std::endl;

  std::cout<<"referenceTemperature: "<<referenceTemperature<<" deg C  ";
  std::cout<<"calibTemperatureIntercept: "<<calibTemperatureIntercept<<" deg C  ";
  std::cout<<"PETemperatureIntercept: "<<PETemperatureIntercept<<" deg C"<<std::endl;
  std::cout<<"         FEB  Channel     PE     PE Temp Corr    FWHM                   Signals              meanTemp    pedestal  calibConst  calibConst Temp Corr"<<std::endl;

  for(int i=0; i<numberOfFebs; i++)
  {
    for(int j=0; j<channelsPerFeb; j++)
    {
      int index=i*channelsPerFeb+j;  //used for _variable[i][j]

      txtFile<<std::fixed;
      txtFile<<std::setprecision(1);
      txtFile<<std::setw(4)<<i<<"  "<<std::setw(4)<<j<<"     ";
      txtFile<<std::setw(8)<<mpvs[0][index]<<"  ";
      txtFile<<std::setw(8)<<mpvs[1][index]<<"     ";
      txtFile<<std::setw(8)<<fwhms[0][index]<<"  ";
      txtFile<<std::setw(8)<<fwhms[1][index]<<"     ";
      txtFile<<std::setprecision(0);
      txtFile<<std::setw(8)<<signals[0][index]<<"  ";
      txtFile<<std::setw(8)<<signals[1][index]<<"     ";
      txtFile<<std::setprecision(1);
      txtFile<<std::setw(8)<<meanTemperatures[index]<<"     ";
      txtFile<<std::setw(8)<<calib.GetPedestals()[index]<<"  ";
      txtFile<<std::setw(8)<<calib.GetCalibrationFactors()[index]<<"  ";
      txtFile<<std::setw(8)<<calib.GetCalibrationFactorsTemperatureCorrected()[index]<<std::endl;

      std::cout<<std::fixed;
      std::cout<<std::setprecision(1);
      std::cout<<"FEB/ch "<<std::setw(4)<<i<<"  "<<std::setw(4)<<j<<"     ";
      std::cout<<std::setw(8)<<mpvs[0][index]<<"  ";
      std::cout<<std::setw(8)<<mpvs[1][index]<<"     ";
      std::cout<<std::setw(8)<<fwhms[0][index]<<"  ";
      std::cout<<std::setw(8)<<fwhms[1][index]<<"     ";
      std::cout<<std::setprecision(0);
      std::cout<<std::setw(8)<<signals[0][index]<<"  ";
      std::cout<<std::setw(8)<<signals[1][index]<<"     ";
      std::cout<<std::setprecision(1);
      std::cout<<std::setw(8)<<meanTemperatures[index]<<"     ";
      std::cout<<std::setw(8)<<calib.GetPedestals()[index]<<"  ";
      std::cout<<std::setw(8)<<calib.GetCalibrationFactors()[index]<<"  ";
      std::cout<<std::setw(8)<<calib.GetCalibrationFactorsTemperatureCorrected()[index]<<std::endl;
    }
  }
  txtFile.close();
  recoTreeSummary->Fill();
  recoTreeSummary->Write("", TObject::kOverwrite);
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
void Summarize(const std::string &pdfFileName, const std::string &txtFileName, const int numberOfFebs, const int channelsPerFeb, 
               const std::vector<float> mpvs[2], const std::vector<float> fwhms[2])
{
  std::ifstream settingsFile;
  settingsFile.open("config.txt");
  if(!settingsFile.is_open()) {std::cerr<<"Could not open config.txt."<<std::endl; exit(1);}

  float minPE, maxPE, minPEfarSide, maxPEfarSide;
  float referenceTemperature, PETemperatureIntercept, calibTemperatureIntercept;
  std::string settingsKey, settingsValue;
  while(settingsFile>>settingsKey>>settingsValue)
  {
    if(settingsKey=="recoMinPE")        minPE=atof(settingsValue.c_str());
    if(settingsKey=="recoMaxPE")        maxPE=atof(settingsValue.c_str());
    if(settingsKey=="recoMinPEfarSide") minPEfarSide=atof(settingsValue.c_str());
    if(settingsKey=="recoMaxPEfarSide") maxPEfarSide=atof(settingsValue.c_str());
    if(settingsKey=="referenceTemperature")   referenceTemperature=atof(settingsValue.c_str());
    if(settingsKey=="PETemperatureIntercept") PETemperatureIntercept=atof(settingsValue.c_str());
    if(settingsKey=="calibTemperatureIntercept") calibTemperatureIntercept=atof(settingsValue.c_str());
  }
  settingsFile.close();

  TCanvas summaryCanvas;
  summaryCanvas.Divide(2,3);

  TH1F *hMpvNear[2];
  TH1F *hMpvFar[2];
  TH1F *hFwhmNear[2];
  TH1F *hFwhmFar[2];
  TPaveText *tMpv[2];
  float meanNear[2], meanFar[2];
  for(int j=0; j<2; ++j)
  {
    summaryCanvas.cd(1+j*2);
    hMpvNear[j]=new TH1F(Form("hMPVnear%i",j),"MPV;PE;Number of Channels",100,0,100);
    hMpvFar[j]=new TH1F(Form("hMPVfar%i",j),"MPV;PE;Number of Channels",100,0,100);
    for(size_t i=0; i<(size_t)2*channelsPerFeb && i<mpvs[j].size(); i++)
    {
      if(mpvs[j][i]>0.0) hMpvNear[j]->Fill(mpvs[j][i]);
    }
    for(size_t i=2*channelsPerFeb; i<mpvs[j].size(); i++)
    {
      if(mpvs[j][i]>0.0) hMpvFar[j]->Fill(mpvs[j][i]);
    }
    TF1 fFar(Form("fitFar%i",j),"gaus(0)");
    fFar.SetLineColor(kRed);
    fFar.SetLineWidth(1);
    TF1 fNear(Form("fitNear%i",j),"gaus(0)");
    fNear.SetLineColor(kRed);
    fNear.SetLineWidth(1);
    hMpvNear[j]->SetMaximum(std::max(hMpvNear[j]->GetMaximum(),hMpvFar[j]->GetMaximum())*1.1);
    hMpvNear[j]->SetLineColor(kBlue);
    hMpvNear[j]->Draw();
    meanNear[j] = 0;
    if(hMpvNear[j]->GetEntries()>0)
    {
      TFitResultPtr fitNear=hMpvNear[j]->Fit(&fNear,"QS");
      if(fitNear.Get()!=0)
      {
        meanNear[j]=fitNear->Parameter(1);
        fNear.Draw("same");
      }
    }
    hMpvFar[j]->SetLineColor(kGreen);
    hMpvFar[j]->Draw("same");
    meanFar[j] = 0;
    if(hMpvFar[j]->GetEntries()>0)
    {
      TFitResultPtr fitFar=hMpvFar[j]->Fit(&fFar,"QS");
      if(fitFar.Get()!=0)
      {
        meanFar[j]=fitFar->Parameter(1);
        fFar.Draw("same");
      }
    }
    tMpv[j]=new TPaveText(0.5, 0.4, 0.85, 0.65, "NDC");
    tMpv[j]->SetFillColor(kWhite);
    tMpv[j]->SetTextAlign(12);
    tMpv[j]->AddText(Form("Near side,  mpv of fit %.1f PEs",meanNear[j]))->SetTextColor(kBlue);
    tMpv[j]->AddText(Form("accepted range: %.0f ... %.0f PEs",minPE,maxPE))->SetTextColor(kBlue);
    tMpv[j]->AddText(Form("Far side,  mpv of fit %.1f PEs",meanFar[j]))->SetTextColor(kGreen);
    tMpv[j]->AddText(Form("accepted range: %.0f ... %.0f PEs",minPEfarSide,maxPEfarSide))->SetTextColor(kGreen);
    tMpv[j]->Draw("same");

    summaryCanvas.cd(2+j*2);
    hFwhmNear[j]=new TH1F(Form("hFWHMnear%i",j),"FWHM;FWHM;Number of Channels",100,0,100);
    hFwhmFar[j]=new TH1F(Form("hFWHMfar%i",j),"FWHM;FWHM;Number of Channels",100,0,100);
    for(size_t i=0; i<(size_t)2*channelsPerFeb && i<fwhms[j].size(); i++)
    {
      if(fwhms[j][i]>0.0) hFwhmNear[j]->Fill(fwhms[j][i]);
    }
    for(size_t i=2*channelsPerFeb; i<fwhms[j].size(); i++)
    {
      if(fwhms[j][i]>0.0) hFwhmFar[j]->Fill(fwhms[j][i]);
    }
    hFwhmNear[j]->SetMaximum(std::max(hFwhmNear[j]->GetMaximum(),hFwhmFar[j]->GetMaximum())*1.1);
    hFwhmNear[j]->SetLineColor(kBlue);
    hFwhmNear[j]->Draw();
    hFwhmFar[j]->SetLineColor(kGreen);
    hFwhmFar[j]->Draw("same");
  }

  TPaveText tT(0.5, 0.65, 0.85, 0.85, "NDC");
  tT.SetFillColor(kWhite);
  tT.SetTextColor(kBlack);
  tT.SetTextAlign(12);
  tT.AddText("Temperature corrected values");
  tT.AddText(Form("refT %.1f deg C, calibT0 %.1f deg C, PET0 %.1f deg C",referenceTemperature,calibTemperatureIntercept, PETemperatureIntercept));
  summaryCanvas.cd(3);
  tT.Draw("same");
  summaryCanvas.cd(4);
  tT.Draw("same");

  summaryCanvas.cd(5);
  TPaveText t1(.0, .6, 1., 1., "NDC");
  t1.SetFillColor(kWhite);
  t1.SetTextColor(kBlack);
  t1.SetTextAlign(12);
  TText *t1Header = t1.AddText("Channels which seem to be dead");
  t1Header->SetTextColor(kRed);
  for(int i=0; i<numberOfFebs; i++)
  {
    std::vector<int> channelNumbers;
    for(int j=0; j<channelsPerFeb; j++)
    {
      int index=i*channelsPerFeb+j;  //used for _variable[i][j]
      if(mpvs[0][index]==0.0) channelNumbers.push_back(j);
    }
    std::string sequence=CreateSequenceString(channelNumbers);
    t1.AddText(Form("FEB %i: %s", i, sequence.c_str()));
  }
  t1.Draw();

  summaryCanvas.cd(6);
  TPaveText t2(.0, .6, 1., 1., "NDC");
  t2.SetFillColor(kWhite);
  t2.SetTextColor(kBlack);
  t2.SetTextAlign(12);
  TText *t2Header = t2.AddText("Channels which have a PE yield outside of the accepted range");
  t2Header->SetTextColor(kRed);
  for(int i=0; i<numberOfFebs; i++)
  {
    std::vector<int> channelNumbers;
    for(int j=0; j<channelsPerFeb; j++)
    {
      int index=i*channelsPerFeb+j;  //used for _variable[i][j]
      if(i<2)
      {
        if(mpvs[0][index]<minPE || mpvs[0][index]>maxPE) channelNumbers.push_back(j);
      }
      else
      {
        if(mpvs[0][index]<minPEfarSide || mpvs[0][index]>maxPEfarSide) channelNumbers.push_back(j);
      }
    }
    std::string sequence=CreateSequenceString(channelNumbers);
    t2.AddText(Form("FEB %i: %s", i, sequence.c_str()));
  }
  t2.Draw();

  summaryCanvas.Print(pdfFileName.c_str(), "pdf");

  for(int j=0; j<2; ++j)
  {
    delete hMpvNear[j];
    delete hMpvFar[j];
    delete hFwhmNear[j];
    delete hFwhmFar[j];
  }

  std::ofstream txtFile;
  txtFile.open(txtFileName.c_str(),std::ios_base::app);
  txtFile<<"Mean near side    "<<std::setw(8)<<meanNear[0]<<"  "<<std::setw(8)<<meanNear[1]<<std::endl;
  txtFile<<"Mean far side    "<<std::setw(8)<<meanFar[0]<<"  "<<std::setw(8)<<meanFar[1]<<std::endl;
  txtFile.close();

  std::cout<<"Mean near side    "<<std::setw(8)<<meanNear[0]<<"  "<<std::setw(8)<<meanNear[1]<<std::endl;
  std::cout<<"Mean far side    "<<std::setw(8)<<meanFar[0]<<"  "<<std::setw(8)<<meanFar[1]<<std::endl;
}

void process(const std::string &runNumber, const std::string &inFileName, const std::string &calibFileName, const std::string &recoFileName, const std::string &pdfFileName, const std::string &txtFileName)
{
  TFile file(inFileName.c_str(), "READ");
  if(!file.IsOpen()) {std::cerr<<"Could not read CRV file for run "<<runNumber<<std::endl; exit(1);}

  TTree *tree = (TTree*)file.Get("run");
  TTree *treeSpills = (TTree*)file.Get("spills");
  //try older tree names
  if(tree==NULL) tree = (TTree*)file.Get(Form("run%04i",atoi(runNumber.c_str())));
  if(treeSpills==NULL) treeSpills = (TTree*)file.Get(Form("run%04i_spills", atoi(runNumber.c_str())));
  if(tree==NULL || treeSpills==NULL) {std::cerr<<"Could not find tree or spill tree"<<std::endl; exit(1);}

  TFile recoFile(recoFileName.c_str(), "RECREATE");
  if(!recoFile.IsOpen()) {std::cerr<<"Could not create reco file for run "<<runNumber<<std::endl; exit(1);}

  TTree *recoTree = new TTree("run","run");
  TTree *recoTreeSpill = treeSpills->CloneTree();
	
  int numberOfFebs;
  int channelsPerFeb;
  int numberOfSamples;
  treeSpills->SetBranchAddress("spill_number_of_febs", &numberOfFebs);
  treeSpills->SetBranchAddress("spill_channels_per_feb", &channelsPerFeb);
  treeSpills->SetBranchAddress("spill_number_of_samples", &numberOfSamples);
  treeSpills->GetEntry(0);  //to read the numberOfFebs, channelsPerFeb, and numberOfSamples

  Calibration calib(calibFileName, numberOfFebs, channelsPerFeb); 
  CrvEvent event(runNumber, numberOfFebs, channelsPerFeb, numberOfSamples, tree, recoTree);

  int nEvents = tree->GetEntries();
//std::cout<<"USING A WRONG NUMBER OF EVENTS"<<std::endl;
//if(nEvents>10000) nEvents=10000;

  TCanvas c0;
  c0.Print(Form("%s[", pdfFileName.c_str()), "pdf");

  for(int i=0; i<nEvents; i++) event.Reconstruct(i, calib);

  std::vector<float> mpvs[2];
  std::vector<float> fwhms[2];
  std::vector<float> signals[2];
  std::vector<float> meanTemperatures;
  for(int feb=0; feb<numberOfFebs; feb++)
  {
    for(int channel=0; channel<channelsPerFeb; channel++)
    {
      TH1F      *h[2];
      TPaveText *t[2];
      for(int i=0; i<2; ++i)
      {
        event.GetCanvas(feb, channel)->cd(i+1);
        h[i]=event.GetHistPEs(i,feb,channel);
//        h[i]->Draw();

        float mpv=0;
        float fwhm=0;
        float nsignals=0;
        LandauGauss(*h[i], mpv, fwhm, nsignals);
        mpvs[i].push_back(mpv);
        fwhms[i].push_back(fwhm);
        signals[i].push_back(nsignals);
        h[i]->GetXaxis()->SetRange(15,150);
        h[i]->Draw();

        t[i]=new TPaveText(0.65,0.65,0.85,0.75,"NDC");
        t[i]->SetFillColor(kWhite);
        t[i]->SetTextColor(kRed);
        t[i]->SetTextAlign(12);
        t[i]->AddText(Form("MPV: %.1f PEs",mpv));
        t[i]->AddText(Form("FWHM: %.1f PEs",fwhm));
        t[i]->AddText(Form("Signals: %.0f",nsignals));
        t[i]->Draw("same");
      }    

      event.GetCanvas(feb, channel)->cd(2);
      TPaveText tt(0.45,0.75,0.85,0.85,"NDC");
      tt.SetFillColor(kWhite);
      tt.SetTextColor(kBlack);
      tt.SetTextAlign(12);
      tt.AddText("Temperature corrected PE values");
      tt.AddText(Form("refT %.1f deg C, calibT0 %.1f deg C, PET0 %.1f deg C",event.GetReferenceTemperature(),event.GetCalibTemperatureIntercept(),event.GetPETemperatureIntercept()));
      tt.Draw("same");

      TGraph *g=event.GetHistTemperatures(feb, channel);
      if(g->GetN()>0)
      {
        event.GetCanvas(feb, channel)->cd(4);
        g->Draw("AP");
        meanTemperatures.push_back(g->GetMean(2));
      }
      else meanTemperatures.push_back(0);

      event.GetCanvas(feb,channel)->Update();
      event.GetCanvas(feb,channel)->Print(pdfFileName.c_str(), "pdf");
    }
  }

  StorePEyields(txtFileName, numberOfFebs, channelsPerFeb, mpvs, fwhms, signals, meanTemperatures, treeSpills, nEvents, calib);

  Summarize(pdfFileName, txtFileName, numberOfFebs, channelsPerFeb, mpvs, fwhms);

  c0.Print(Form("%s]", pdfFileName.c_str()), "pdf");

  recoTree->Write("", TObject::kOverwrite);
  recoTreeSpill->Write("", TObject::kOverwrite);

  recoFile.Close();
  file.Close();
}

void makeFileNames(const std::string &runNumber, std::string &inFileName, std::string &calibFileName, std::string &recoFileName, std::string &pdfFileName, std::string &txtFileName)
{
  std::ifstream dirFile;
  dirFile.open("config.txt");
  if(!dirFile.is_open()) {std::cerr<<"Could not open config.txt."<<std::endl; exit(1);}

  std::string crvDirName, calibDirName, recoDirName;
  std::string dirType, dir;
  while(dirFile>>dirType>>dir)
  {
    if(dirType=="crvparsed") crvDirName=dir;
    if(dirType=="crvcalib")  calibDirName=dir;
    if(dirType=="crvreco")   recoDirName=dir;
  }
  dirFile.close();

  bool found=false;
  for(const auto& dirEntry : std::experimental::filesystem::directory_iterator(crvDirName)) 
  {
    const std::string s = dirEntry.path().filename().string();
    const std::string s0 = "crv.parsed.";
    const std::string s1 = ".run"+runNumber+".root";
    if(s.compare(0,s0.length(),s0)!=0) continue;
    if(s.compare(s.length()-s1.length(),s1.length(),s1)!=0) continue;

    found=true;
    inFileName = dirEntry.path().string();
    calibFileName = calibDirName+"crv.calib."+dirEntry.path().stem().string().substr(s0.length())+".txt";
    recoFileName = recoDirName+"crv.reco."+dirEntry.path().stem().string().substr(s0.length())+".root";
    pdfFileName = recoDirName+"log.crv.reco."+dirEntry.path().stem().string().substr(s0.length())+".pdf";
    txtFileName = recoDirName+"log.crv.reco."+dirEntry.path().stem().string().substr(s0.length())+".txt";
    break;
  }

  if(!found)
  {
//try Ray's version of file names
    for(const auto& dirEntry : std::experimental::filesystem::directory_iterator(crvDirName)) 
    {
      const std::string s = dirEntry.path().filename().string();
      const std::string s0 = "ntd.mu2e.";
      const std::string s1 = "."+runNumber+".root";
      if(s.compare(0,s0.length(),s0)!=0) continue;
      if(s.compare(s.length()-s1.length(),s1.length(),s1)!=0) continue;

      found=true;
      inFileName = dirEntry.path().string();
      calibFileName = calibDirName+"cal.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".txt";
      recoFileName = recoDirName+"rec.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".root";
      pdfFileName = recoDirName+"rec.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".pdf";
      txtFileName = recoDirName+"rec.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".txt";
      break;
    }
  }

  if(!found) {std::cerr<<"Could not open input file for run "<<runNumber<<"."<<std::endl; exit(1);}
}

void printHelp()
{
  std::cout<<"Use as"<<std::endl;
  std::cout<<"recoCrv -h             Prints this help."<<std::endl;
  std::cout<<"recoCrv RUNNUMBER      Reconstructs a run."<<std::endl;
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
  std::string recoFileName;
  std::string pdfFileName;
  std::string txtFileName;
  makeFileNames(runNumber, inFileName, calibFileName, recoFileName, pdfFileName, txtFileName);

  process(runNumber, inFileName, calibFileName, recoFileName, pdfFileName, txtFileName);

//#pragma message "USING WRONG CALIB CONSTANTS"
//std::cout<<"USING WRONG CALIB CONSTANTS"<<std::endl;
//#pragma message "USING WRONG SIGNAL REGION"
//std::cout<<"USING WRONG SIGNAL REGION"<<std::endl;
  return 0;
}
