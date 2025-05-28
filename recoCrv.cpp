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
# include "TKey.h"
# include "TLine.h"
# include "TPaveText.h"
# include "TFitResult.h"
# include "TSystem.h"

const float CRV_TDC_RATE    = 159.324e6;  // Hz
const float RATE=(CRV_TDC_RATE/2.0)/1.0e9; // GHZ
const int BOARD_STATUS_REGISTERS=22;
const int FPGA_BLOCK_REGISTERS=38;
const int FPGA_BLOCKS=4;
const float DEFAULT_BETA=16.0;
const float COUNTER_WIDTH=51.34; //mm

struct TemperatureCorrections
{
  double PEOvervoltageChange{0.229}; //1/V
  double calibOvervoltageChange{125.5};   //ADC*ns/V
  double calibTemperatureChangeAFE{-1.46};  //ADC*ns/K  additional calib change due to temperature change in AFE
  double overvoltageTemperatureChangeCMB{-0.0554};  //V/K
  double overvoltageTemperatureChangeFEB{0.00409};  //V/K
  double referenceTemperatureCMB{20.0};   //degC
  double referenceTemperatureFEB{40.0};   //degC
};

struct ChannelStruct
{
  int _side;
  int _sector;
  float _x,_y;
  int _ignoreForFit;
  ChannelStruct() : _side(0), _sector(0), _x(0), _y(0), _ignoreForFit(0) {}
  ChannelStruct(int side, int sector, float x, float y, int ignoreForFit) : _side(side), _sector(sector), _x(x), _y(y), _ignoreForFit(ignoreForFit) {}
};
typedef std::map<std::pair<int,std::pair<int,int> >, ChannelStruct>  ChannelMapType;   //feb,(channel1,channel2) --> channelStruct

class Calibration
{
  public:
  Calibration(const std::string &calibFileName, const int numberOfFebs, const int channelsPerFeb);
  bool  IsNewFormat() const {return _newFormat;}
  const std::vector<float> &GetPedestals() const {return _pedestal;}
  const std::vector<float> &GetCalibrationFactors() const {return _calibrationFactor;}
  const std::vector<float> &GetCalibrationFactorsTemperatureCorrected() const {return _calibrationFactorTemperatureCorrected;}
  const std::vector<float> &GetNoiseRate() const {return _noise;}
  const std::vector<float> &GetXtalkProbability() const {return _xtalk;}

  private:
  int                 _numberOfFebs;
  int                 _channelsPerFeb;
  bool                _newFormat;  //temperature corrected
  std::vector<float>  _pedestal;
  std::vector<float>  _calibrationFactor;
  std::vector<float>  _calibrationFactorTemperatureCorrected;
  std::vector<float>  _noise;
  std::vector<float>  _xtalk;

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
  _noise.resize(_numberOfFebs*_channelsPerFeb);
  _xtalk.resize(_numberOfFebs*_channelsPerFeb);
  std::string tmp;

  getline(calibFile,tmp);
  if(tmp.find("calib v2")==0)
  {
    getline(calibFile,tmp); //table header
    _newFormat=true;
    int i, j;
    float pedestalTmp, calibrationFactorTmp, calibrationFactorTemperatureCorrectedTmp, noiseTmp, xtalkTmp;
    while(calibFile >> i >> j >> pedestalTmp >> calibrationFactorTmp >> calibrationFactorTemperatureCorrectedTmp >> noiseTmp >> xtalkTmp)
    {
      _pedestal.at(i*_channelsPerFeb+j)=pedestalTmp;
      _calibrationFactor.at(i*_channelsPerFeb+j)=calibrationFactorTmp;
      _calibrationFactorTemperatureCorrected.at(i*_channelsPerFeb+j)=calibrationFactorTemperatureCorrectedTmp;
      _noise.at(i*_channelsPerFeb+j)=noiseTmp;
      _xtalk.at(i*_channelsPerFeb+j)=xtalkTmp;
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
  void PeakFitter(const short* data, int numberOfSamples, float pedestal, float calibrationFactor, bool &draw);

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
void CrvRecoEvent::PeakFitter(const short* data, int numberOfSamples, float pedestal, float calibrationFactor, bool &draw)
{
  if(std::isnan(calibrationFactor) || calibrationFactor==0) return;

//  if(data[0]==0 && data[1]==0 && data[2]==0 && data[3]==0) return; //FIXME temporary check for bad events
//                                                                   //where other channels work, so that timSinceSpill wasn't marked as NAN

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
           TTree *tree, TTree *recoTree, const TemperatureCorrections &temperatureCorrections, float PEcut,
           const std::string &channelMapFile);
  void     Reconstruct(int entry, const Calibration &calib);
  TCanvas *GetCanvas(int feb, int channel) {return _canvas[feb*_channelsPerFeb+channel];}
  TH1F    *GetHistPEs(int i, int feb, int channel)
           {
             return (i==0?_histPEs[feb*_channelsPerFeb+channel]:_histPEsTemperatureCorrected[feb*_channelsPerFeb+channel]);
           }
  TGraph  *GetHistTemperatures(int feb, int channel) {return _histTemperatures[feb*_channelsPerFeb+channel];}
  TGraph  *GetHistTemperaturesFEB(int feb, int channel) {return _histTemperaturesFEB[feb*_channelsPerFeb+channel];}
  void     ReadChannelMap(const std::string &channelMapFile);
  void     TrackFit();
  int      GetMaxedOutEvents(int feb, int channel) {return _maxedOut[feb*_channelsPerFeb+channel];}

  private:
  CrvEvent();

  int   _signalRegionStart;
  int   _signalRegionEnd;

  std::string _runNumber;
  int _run;
  int _subrun;
  int _numberOfFebs;
  int _channelsPerFeb;
  int _numberOfSamples;

  TTree *_tree;
  TTree *_recoTree;

  const TemperatureCorrections _tc;

  int     _spillIndex;
  int     _spillNumber;
  int    *_lastSpillIndex;
  int     _eventNumber;
  //int    *_tdcSinceSpill;  //OLD
  Long64_t *_tdcSinceSpill;
  double *_timeSinceSpill;
  short  *_adc;
  float  *_temperature;
  int    *_boardStatus;
  int    *_FPGABlocks;
  Long64_t  _timestamp;  //time_t

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
  int    *_maxedOut;

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
  std::vector<TGraph*>  _histTemperaturesFEB;

  //for track fits (arrays are sorted by sectors; entry 0 is used for all sectors)
  ChannelMapType _channelMap;
  int    _numberOfSectors;
  float  _PEcut;
  float *_trackSlope;      //using slope=dx/dy to avoid inf for vertical tracks
  float *_trackIntercept;  //x value, where y=0
  float *_trackChi2;
  int   *_trackPoints;
  float *_trackPEs;
};
CrvEvent::CrvEvent(const std::string &runNumber, const int numberOfFebs, const int channelsPerFeb, const int numberOfSamples,
                   TTree *tree, TTree *recoTree, const TemperatureCorrections &temperatureCorrections, float PEcut,
                   const std::string &channelMapFile) :
                   _runNumber(runNumber), _numberOfFebs(numberOfFebs), _channelsPerFeb(channelsPerFeb), _numberOfSamples(numberOfSamples),
                   _tree(tree), _recoTree(recoTree), _tc(temperatureCorrections), _PEcut(PEcut)
{
  std::ifstream configFile;
  configFile.open("config.txt");
  if(!configFile.is_open()) {std::cerr<<"Could not open config.txt."<<std::endl; exit(1);}

  std::string configKey, configValue;
  while(configFile>>configKey>>configValue)
  {
    if(configKey=="signalRegionStart") _signalRegionStart=atoi(configValue.c_str());
    if(configKey=="signalRegionEnd")   _signalRegionEnd=atoi(configValue.c_str());
  }
  configFile.close();

  _lastSpillIndex = new int[_numberOfFebs*_channelsPerFeb];
//  _tdcSinceSpill  = new int[_numberOfFebs];  //OLD
//  _timeSinceSpill = new double[_numberOfFebs];  //OLD
  _tdcSinceSpill  = new Long64_t[_numberOfFebs*_channelsPerFeb];
  _timeSinceSpill = new double[_numberOfFebs*channelsPerFeb];
  _adc            = new short[_numberOfFebs*_channelsPerFeb*_numberOfSamples];
  _temperature    = new float[_numberOfFebs*_channelsPerFeb];
  _boardStatus    = new int[_numberOfFebs*BOARD_STATUS_REGISTERS];
  _FPGABlocks     = new int[_numberOfFebs*FPGA_BLOCKS*FPGA_BLOCK_REGISTERS];

  tree->SetBranchAddress("runtree_run_num", &_run);
  tree->SetBranchAddress("runtree_subrun_num", &_subrun);
  tree->SetBranchAddress("runtree_spill_index", &_spillIndex);
  tree->SetBranchAddress("runtree_spill_num", &_spillNumber);
  tree->SetBranchAddress("runtree_event_num",&_eventNumber);
  tree->SetBranchAddress("runtree_tdc_since_spill", _tdcSinceSpill);
  tree->SetBranchAddress("runtree_time_since_spill", _timeSinceSpill);
  tree->SetBranchAddress("runtree_adc", _adc);
  tree->SetBranchAddress("runtree_temperature", _temperature);
  tree->SetBranchAddress("runtree_boardStatus", _boardStatus);
  tree->SetBranchAddress("runtree_FPGABlocks", _FPGABlocks);
  tree->SetBranchAddress("runtree_spillTimestamp", &_timestamp);

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
  _maxedOut                = new int[_numberOfFebs*_channelsPerFeb];

  _fitStatusReflectedPulse               = new int[_numberOfFebs*_channelsPerFeb];
  _PEsReflectedPulse                     = new float[_numberOfFebs*_channelsPerFeb];
  _PEsTemperatureCorrectedReflectedPulse = new float[_numberOfFebs*_channelsPerFeb];
  _pulseHeightReflectedPulse             = new float[_numberOfFebs*_channelsPerFeb];
  _betaReflectedPulse                    = new float[_numberOfFebs*_channelsPerFeb];
  _timeReflectedPulse                    = new float[_numberOfFebs*_channelsPerFeb];
  _LEtimeReflectedPulse                  = new float[_numberOfFebs*_channelsPerFeb];
  _recoStartBinReflectedPulse            = new int[_numberOfFebs*_channelsPerFeb];
  _recoEndBinReflectedPulse              = new int[_numberOfFebs*_channelsPerFeb];

  _numberOfSectors=0;
  if(channelMapFile!="") ReadChannelMap(channelMapFile);
  _trackSlope     = new float[_numberOfSectors+1];
  _trackIntercept = new float[_numberOfSectors+1];
  _trackChi2      = new float[_numberOfSectors+1];
  _trackPoints    = new int[_numberOfSectors+1];
  _trackPEs       = new float[_numberOfSectors+1];

  recoTree->Branch("runNumber", &_run, "runNumber/I");
  recoTree->Branch("subrunNumber", &_subrun, "subrunNumber/I");
  recoTree->Branch("spillIndex", &_spillIndex, "spillIndex/I");
  recoTree->Branch("spillNumber", &_spillNumber, "spillNumber/I");
  recoTree->Branch("boardStatus", _boardStatus, Form("boardStatus[%i][%i]/I",_numberOfFebs,BOARD_STATUS_REGISTERS));
  recoTree->Branch("FPGABlocks", _FPGABlocks, Form("FPGABlocks[%i][%i][%i]/I",_numberOfFebs,FPGA_BLOCKS,FPGA_BLOCK_REGISTERS));
  recoTree->Branch("spillTimestamp", &_timestamp, "spillTimestamp/L");
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
  recoTree->Branch("adc", _adc, Form("adc[%i][%i][%i]/S",_numberOfFebs,_channelsPerFeb,_numberOfSamples));
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
  recoTree->Branch("trackSlope", _trackSlope, Form("trackSlope[%i]/F",_numberOfSectors+1));
  recoTree->Branch("trackIntercept", _trackIntercept, Form("trackIntercept[%i]/F",_numberOfSectors+1));
  recoTree->Branch("trackChi2", _trackChi2, Form("trackChi2[%i]/F",_numberOfSectors+1));
  recoTree->Branch("trackPoints", _trackPoints, Form("trackPoints[%i]/I",_numberOfSectors+1));
  recoTree->Branch("trackPEs", _trackPEs, Form("trackPEs[%i]/F",_numberOfSectors+1));

  _canvas.resize(_numberOfFebs*_channelsPerFeb);
  _plot.resize(_numberOfFebs*_channelsPerFeb);
  _histPEs.resize(_numberOfFebs*_channelsPerFeb);
  _histPEsTemperatureCorrected.resize(_numberOfFebs*_channelsPerFeb);
  _histTemperatures.resize(_numberOfFebs*_channelsPerFeb);
  _histTemperaturesFEB.resize(_numberOfFebs*_channelsPerFeb);
  for(int i=0; i<_numberOfFebs; i++)
  {
    for(int j=0; j<_channelsPerFeb; j++)
    {
      int index=i*_channelsPerFeb+j;  //used for _variable[i][j]
      _maxedOut[index]=0;
      _lastSpillIndex[index]=-1;
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
      _canvas[index]->cd(4);
      gPad->Divide(1,2);

      _histPEs[index]=new TH1F(Form("PEs_%i_%i",i,j),Form("PE Distribution Run %s FEB %i Channel %i;PE;Counts",_runNumber.c_str(),i,j),75,0,150);
      _histPEsTemperatureCorrected[index]=new TH1F(Form("PEsTempCorrected_%i_%i",i,j),Form("PE Distribution (temp. corrected) Run %s FEB %i Channel %i;PE;Counts",_runNumber.c_str(),i,j),75,0,150);
      _histTemperatures[index]=new TGraph();
      _histTemperaturesFEB[index]=new TGraph();

      _histPEs[index]->SetLineColor(kBlack);
      _histPEsTemperatureCorrected[index]->SetLineColor(kBlack);
      _histTemperatures[index]->SetMarkerStyle(20);
      _histTemperatures[index]->SetMarkerSize(0.5);
      _histTemperatures[index]->SetMarkerColor(kBlack);
      _histTemperatures[index]->SetNameTitle(Form("hT_%i",index),Form("Temperature Run %s FEB %i Channel %i;Spill;Temperature [deg C]",_runNumber.c_str(),i,j));
      _histTemperaturesFEB[index]->SetMarkerStyle(20);
      _histTemperaturesFEB[index]->SetMarkerSize(0.5);
      _histTemperaturesFEB[index]->SetMarkerColor(kBlack);
      _histTemperaturesFEB[index]->SetNameTitle(Form("hTFEB_%i",index),Form("FEB Temperature Run %s FEB %i Channel %i;Spill;Temperature [deg C]",_runNumber.c_str(),i,j));
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

      if(!std::isnan(_timeSinceSpill[index]))  //missing FEB/channel in raw data
        reco.PeakFitter(&(_adc[index*_numberOfSamples]), _numberOfSamples, pedestal, calibrationFactor, draw);

      //main pulse
      _fitStatus[index]               = reco._fitStatus[0];
      _PEs[index]                     = reco._PEs[0];
      _PEsTemperatureCorrected[index] = -1;

      double temperatureFEB=-1000;
      if(_boardStatus[i*BOARD_STATUS_REGISTERS]!=-1) //i-th FEB was read for this spill
      {
        //temperature of i-th FEB
        temperatureFEB=_boardStatus[i*BOARD_STATUS_REGISTERS+2]*0.01;  //TODO: document seems to indicate a factor of 10.0
      }

      if(fabs(_temperature[index])<100 && fabs(temperatureFEB)<100) //temperature of -1000 means no temperature found
      {
        //overvoltage difference for actual CMB and FEB temperatures w.r.t. reference CMB and FEB temperatures
        //deltaOvervoltage = overvoltageTemperatureChangeCMB*(TCMB-TrefCMB) + overvoltageTemperatureChangeFEB*(TFEB-TrefFEB)
        float deltaOvervoltage = _tc.overvoltageTemperatureChangeCMB*(_temperature[index]-_tc.referenceTemperatureCMB)
                               + _tc.overvoltageTemperatureChangeFEB*(temperatureFEB-_tc.referenceTemperatureFEB);

        _PEsTemperatureCorrected[index] = reco._PEs[0];

        if(calib.IsNewFormat() && calibrationFactorTemperatureCorrected!=0) //The old calib format doesn't have temperature-corrected calib constants.
        {
          //The stored temperature corrected calibration constants were adjusted to CMB/FEB reference temperatures TrefCMB/TrefFEB of e.g. 20 degC / 40 degC.
          //The calibration constant needs to be adjusted to the temperature T at which the signal pulse happened.
          //calibConst(TCMB,TFEB) = calibConst(TrefCMB,TrefFEB) + calibOvervoltageChange*deltaOvervoltage(TCMB,TFEB) + calibTempChangeAFE*(TFEB-TrefFEB)
          float calibrationFactorAtT = calibrationFactorTemperatureCorrected + _tc.calibOvervoltageChange*deltaOvervoltage
                                                                             + _tc.calibTemperatureChangeAFE*(temperatureFEB-_tc.referenceTemperatureFEB);
          //The PeakFitter used the un-corrected calibration factor to calculate the PEs
          //This needs to be undone first to get the pulse area
          //pulseArea = PEsFromPeakFitter * uncorrectedCalibrationFactor
          //PEsWithCorrectCalibConstant = pulseArea / calibrationFactorAtT
          _PEsTemperatureCorrected[index]*=calibrationFactor/calibrationFactorAtT;
        }

        //The PEs that were measured at a particular temperature need to be adjusted to the reference temperature of e.g. 20 degC
        //PE(TCMB,TFEB)/PE(TrefCMB,TrefFEB) = 1.0 + PEOvervoltageChange*deltaOvervoltage(TCMB,TFEB)
        _PEsTemperatureCorrected[index]/=1.0 + _tc.PEOvervoltageChange*deltaOvervoltage;
      }

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
      if(_temperature[index]>-300 && temperatureFEB>-300) //temperature of -1000 means no temperature found
      {
        float deltaOvervoltage = _tc.overvoltageTemperatureChangeCMB*(_temperature[index]-_tc.referenceTemperatureCMB)
                               + _tc.overvoltageTemperatureChangeFEB*(temperatureFEB-_tc.referenceTemperatureFEB);
        _PEsTemperatureCorrectedReflectedPulse[index] = reco._PEs[1];
        if(calib.IsNewFormat() && calibrationFactorTemperatureCorrected!=0)
        {
          float calibrationFactorAtT = calibrationFactorTemperatureCorrected + _tc.calibOvervoltageChange*deltaOvervoltage
                                                                             + _tc.calibTemperatureChangeAFE*(temperatureFEB-_tc.referenceTemperatureFEB);
          _PEsTemperatureCorrectedReflectedPulse[index]*=calibrationFactor/calibrationFactorAtT;
        }
        _PEsTemperatureCorrectedReflectedPulse[index]/=1.0 + _tc.PEOvervoltageChange*deltaOvervoltage;
      }
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
      if(_temperature[index]>-300 && _lastSpillIndex[index]!=_spillIndex) //temperature of -1000 means no temperature found
      {
        _histTemperatures[index]->SetPoint(_histTemperatures[index]->GetN(),_spillIndex,_temperature[index]);
        double temperatureFEB=_boardStatus[i*BOARD_STATUS_REGISTERS+2]*0.01;  //TODO: document seems to indicate a factor of 10.0
        _histTemperaturesFEB[index]->SetPoint(_histTemperaturesFEB[index]->GetN(),_spillIndex,temperatureFEB);
        _lastSpillIndex[index]=_spillIndex;
      }

      //check, if ADCs have maxed-out
      for(int iSample=0; iSample<_numberOfSamples; ++iSample)
      {
        if(_adc[index*_numberOfSamples+iSample]==2047) {++_maxedOut[index]; break;}
      }
    }
  }

  //Fit track, if channel map is provided
  if(!_channelMap.empty())
  {
    TrackFit();
  }

  _recoTree->Fill();
}

void CrvEvent::TrackFit()
{
  std::vector<float> sumX;
  std::vector<float> sumY;
  std::vector<float> sumXY;
  std::vector<float> sumYY;
  sumX.resize(_numberOfSectors+1);
  sumY.resize(_numberOfSectors+1);
  sumXY.resize(_numberOfSectors+1);
  sumYY.resize(_numberOfSectors+1);

  for(int sector=0; sector<=_numberOfSectors; ++sector)
  {
    sumX[sector]=0;
    sumY[sector]=0;
    sumXY[sector]=0;
    sumYY[sector]=0;
    _trackSlope[sector]=0;
    _trackIntercept[sector]=0;
    _trackChi2[sector]=-1;
    _trackPoints[sector]=0;
    _trackPEs[sector]=0;
  }

  //loop through the channel map
  for(ChannelMapType::iterator channelIter=_channelMap.begin(); channelIter!=_channelMap.end(); ++channelIter)
  {
    int feb=channelIter->first.first;
    int channel1=channelIter->first.second.first;
    int channel2=channelIter->first.second.second;
    int sector=channelIter->second._sector;
    int ignoreForFit=channelIter->second._ignoreForFit;
    float x=channelIter->second._x;
    float y=channelIter->second._y;

    if(feb<0 || feb>=_numberOfFebs) continue;  //feb not in event tree
    if(channel1<0 || channel2<0 || channel1>=_channelsPerFeb || channel2>=_channelsPerFeb) continue;  //channels not in event tree

    int index1=feb*_channelsPerFeb+channel1;  //used for _variable[i][j]
    int index2=feb*_channelsPerFeb+channel2;  //used for _variable[i][j]

    float PE1     = _PEsTemperatureCorrected[index1];
    float PE2     = _PEsTemperatureCorrected[index2];
    if(PE1<=0 || _fitStatus[index1]==0) PE1=0;
    if(PE2<=0 || _fitStatus[index2]==0) PE2=0;

    //collect information for the fit
    //uses hit information from both sides
    float PE  = PE1+PE2;
    if(PE<_PEcut) continue;

    if(ignoreForFit==0)
    {
      sumX[0] +=x*PE;
      sumY[0] +=y*PE;
      sumXY[0]+=x*y*PE;
      sumYY[0]+=y*y*PE;
      _trackPEs[0]+=PE;
      ++_trackPoints[0];
    }

    if(sector!=0)
    {
      sumX[sector] +=x*PE;
      sumY[sector] +=y*PE;
      sumXY[sector]+=x*y*PE;
      sumYY[sector]+=y*y*PE;
      _trackPEs[sector]+=PE;
      ++_trackPoints[sector];
    }
  }

//do the fit
  for(int sector=0; sector<=_numberOfSectors; ++sector)
  {
    if(_trackPEs[sector]>=2*_PEcut && _trackPoints[sector]>1)
    {
      if(_trackPEs[sector]*sumYY[sector]-sumY[sector]*sumY[sector]!=0)
      {
        _trackSlope[sector]=(_trackPEs[sector]*sumXY[sector]-sumX[sector]*sumY[sector])/(_trackPEs[sector]*sumYY[sector]-sumY[sector]*sumY[sector]);
        _trackIntercept[sector]=(sumX[sector]-_trackSlope[sector]*sumY[sector])/_trackPEs[sector];

        //find chi2
        _trackChi2[sector]=0;
        for(ChannelMapType::iterator channelIter=_channelMap.begin(); channelIter!=_channelMap.end(); ++channelIter)
        {
          int feb=channelIter->first.first;
          int channel1=channelIter->first.second.first;
          int channel2=channelIter->first.second.second;
          int thisSector=channelIter->second._sector;
          int ignoreForFit=channelIter->second._ignoreForFit;
          if(sector!=thisSector && sector!=0) continue;
          if(sector==0 && ignoreForFit!=0) continue;

          float x=channelIter->second._x;
          float y=channelIter->second._y;

          if(feb<0 || feb>=_numberOfFebs) continue;  //feb not in event tree
          if(channel1<0 || channel2<0 || channel1>=_channelsPerFeb || channel2>=_channelsPerFeb) continue;  //channels not in event tree

          int index1=feb*_channelsPerFeb+channel1;  //used for _variable[i][j]
          int index2=feb*_channelsPerFeb+channel2;  //used for _variable[i][j]

          float PE1     = _PEsTemperatureCorrected[index1];
          float PE2     = _PEsTemperatureCorrected[index2];
          if(PE1<=0 || _fitStatus[index1]==0) PE1=0;
          if(PE2<=0 || _fitStatus[index2]==0) PE2=0;

          float PE  = PE1+PE2;
          if(PE<_PEcut) continue;

          float xFit = _trackSlope[sector]*y + _trackIntercept[sector];
          _trackChi2[sector]+=(xFit-x)*(xFit-x)*PE;  //PE-weighted chi2
        }
        _trackChi2[sector]/=(COUNTER_WIDTH*COUNTER_WIDTH/12.0)*(_trackPEs[sector]/_trackPoints[sector]);
      }
    }
  }
}

void CrvEvent::ReadChannelMap(const std::string &channelMapFile)
{
  std::ifstream file(channelMapFile);
  if(!file.is_open()) {std::cout<<"Channel map file could not be opened."<<std::endl; return;}

  std::string header;
  getline(file,header);
  bool hasSectors=false;
  if(header.find("sector")!=std::string::npos) hasSectors=true;
  bool hasIgnoreForFit=false;
  if(header.find("ignoreForFit")!=std::string::npos) hasIgnoreForFit=true;

  int febA, channelA1, channelA2;
  int febB, channelB1, channelB2;
  float x, y;
  int sector=0;
  int ignoreForFit=0;
  _numberOfSectors=0;

  while(file >> febA >> channelA1 >> channelA2 >> febB >> channelB1 >> channelB2 >> x >> y)
  {
    if(hasSectors) file >> sector;
    if(sector>_numberOfSectors) _numberOfSectors=sector;
    if(hasIgnoreForFit) file >> ignoreForFit;

    std::pair<int,int> channelPairA(channelA1,channelA2);
    std::pair<int,int> channelPairB(channelB1,channelB2);
    std::pair<int,std::pair<int,int> > counterA(febA,channelPairA);
    std::pair<int,std::pair<int,int> > counterB(febB,channelPairB);
    std::pair<float,float> counterPos(x,y);
    _channelMap[counterA]=ChannelStruct(0,sector,x,y,ignoreForFit);
    _channelMap[counterB]=ChannelStruct(1,sector,x,y,ignoreForFit);
  }

  file.close();
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
void LandauGauss(TH1F &h, float &mpv, float &fwhm, float &signals, float &chi2, float &error)
{
    std::multimap<float,float> bins;  //binContent,binCenter
    for(int i=8; i<=h.GetNbinsX(); i++) bins.emplace(h.GetBinContent(i),h.GetBinCenter(i));  //ordered from smallest to largest bin entries
    if(bins.size()<4) return;
    if(bins.rbegin()->first<20) return;  //low statistics
    int nBins=0;
    float binSum=0;
    for(auto bin=bins.rbegin(); bin!=bins.rend(); ++bin)
    {
      nBins++;
      binSum+=bin->second;
      if(nBins==4) break;
    }
    float maxX=binSum/4;
    float fitRangeStart=0.7*maxX;  //0.6 @ 24
    float fitRangeEnd  =2.0*maxX;
    if(fitRangeStart<15.0) fitRangeStart=15.0;

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
    startValues[0]=5.0;   startValues[3]=10.0;
    parLimitsLow[0]=2.0;  parLimitsLow[3]=2.0;
    parLimitsHigh[0]=15.0; parLimitsHigh[3]=20.0; //7 and 15 @ 21  //6 and 13 @ 23

    TF1 fit("LandauGauss",LandauGaussFunction,fitRangeStart,fitRangeEnd,4);
    fit.SetParameters(startValues);
    fit.SetLineColor(kRed);
    fit.SetParNames("Width","MP","Area","GSigma");
    for(int i=0; i<4; i++) fit.SetParLimits(i, parLimitsLow[i], parLimitsHigh[i]);
    TFitResultPtr fr = h.Fit(&fit,"LQRS");
    fit.Draw("same");

    mpv = fit.GetMaximumX();  //fit.GetParameter(1) does not give the correct result, probably because it is not a pure Landau function
    chi2 = (fr->Ndf()>0?fr->Chi2()/fr->Ndf():NAN);
    if(mpv==fitRangeStart) {mpv=0; return;}
    float halfMaximum = fit.Eval(mpv)/2.0;
    float leftX = fit.GetX(halfMaximum,0.0,mpv);
    float rightX = fit.GetX(halfMaximum,mpv,10.0*mpv);
    fwhm = rightX-leftX;
    signals = fit.Integral(0,150,1e-3)/h.GetBinWidth(1);  //need to divide by bin width.
                                                     //if the bin width is 2 and one has e.g. 20 events for 50PEs and 20 events for 51PEs,
                                                     //the combined bin of x=50/51 gets 40 entries and the integral assumes that there are 40 entries for x=50 and x=51.
    error = fit.GetParError(1);
}
double PoissonFunction(double *x, double *par)
{
    //par[1] is the single parameter of the original TMath::Poisson
    //par[0] is the vertial scaling factor
    //par[2] is the horizontal scaling factor
    if(par[0]<0 || par[1]<0 || par[2]<0) return TMath::QuietNaN();
    if(x[0]<0) return 0;
    if(x[0]==0.0) return TMath::Exp(-par[1]);

    const double scaledX=x[0]*par[2];
    return par[0]*TMath::Exp(scaledX * log(par[1]) - TMath::LnGamma(scaledX + 1.) - par[1]);
}
void Poisson(TH1F &h, float &mpv, float &fwhm, float &signals, float &chi2, float &error)
{
    float maxX=0;
    float maxValue=0;
    for(int i=0; i<=h.GetNbinsX(); i++)
    {
      if(h.GetBinContent(i)>maxValue) {maxValue=h.GetBinContent(i); maxX=h.GetBinCenter(i);}
    }
    if(maxValue<20) return;

    //Fit range
    Double_t fitRangeStart=0.7*maxX;
    Double_t fitRangeEnd  =1.3*maxX;

    TF1 fit("Poisson",PoissonFunction,fitRangeStart,fitRangeEnd,3);
    fit.SetParameter(0,maxValue/TMath::Poisson(maxX,maxX));  //"height" of the function
    fit.SetParameter(1,maxX);   //expected value
    fit.SetParameter(2,1.0);    //horizontal (x) scaling factor
    fit.SetLineColor(kRed);
    fit.SetParNames("Const","mean");
    TFitResultPtr fr = h.Fit(&fit,"QRS");
    fit.Draw("same");

    mpv = fit.GetMaximumX();
    chi2 = (fr->Ndf()>0?fr->Chi2()/fr->Ndf():NAN);
    if(mpv==fitRangeStart) {mpv=0; return;}
    float halfMaximum = fit.Eval(mpv)/2.0;
    float leftX = fit.GetX(halfMaximum,0.0,mpv);
    float rightX = fit.GetX(halfMaximum,mpv,10.0*mpv);
    fwhm = rightX-leftX;
    signals = fit.Integral(0,150,1e-3)/h.GetBinWidth(1); //explanation see Landau-Gauss
    error = fit.GetParError(1);
    if(fit.GetParameter(2)!=0) error/=fit.GetParameter(2);
}

void BoardRegisters(TTree *treeSpills, std::ofstream &txtFile, const int numberOfFebs, int *febID, int &nSpillsActual, int *nFebSpillsActual,
                    float *febTemperaturesAvg, float *supplyMonitorsAvg, float *biasVoltagesAvg, int *pipeline, int *samples, Long64_t &timestamp) //time_t
{
  if(!treeSpills->GetBranch("spill_boardStatus")) return;  //older file: board status was not stored

  Long64_t  timestampTmp;  //time_t
  int   nEventsExpected;
  int   nEventsActual;
  bool  spillStored;
  bool  timestampFound=false;
  int  *boardStatus = new int[numberOfFebs*BOARD_STATUS_REGISTERS];
  treeSpills->SetBranchAddress("spill_nevents", &nEventsExpected);
  treeSpills->SetBranchAddress("spill_neventsActual", &nEventsActual);
  treeSpills->SetBranchAddress("spill_stored", &spillStored);
  treeSpills->SetBranchAddress("spill_boardStatus", boardStatus);
  treeSpills->SetBranchAddress("spill_timestamp", &timestampTmp);

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
    if(timestampTmp!=0 && !timestampFound)  //the first spill may not have a time stamp
    {
      timestampFound=true;
      timestamp=timestampTmp;  //to get the time stamp, that gets overwritten at the next GetEntry calls
      long timestampL=timestamp;  //long required below
      txtFile<<"timestamp: "<<ctime(&timestampL);
    }

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
  txtFile<<"FEB  ID   spills  FEBtemp 1.2Vmon  1.8Vmon  5Vmon  10Vmon  2.5Vmon -5Vmon  15Vmon  3.3Vmon"<<std::endl;
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

void StorePEyields(const std::string &runNumber, const std::string &txtFileName, const int numberOfFebs, const int channelsPerFeb,
                   const std::vector<float> mpvs[2], const std::vector<float> fwhms[2], const std::vector<float> signals[2],
                   const std::vector<float> chi2s[2], const std::vector<float> errors[2],
                   const std::vector<float> &meanTemperatures, const std::vector<float> &stddevTemperatures, const std::vector<float> &maxedOutFraction,
                   TTree *treeSpills, int nEventsActual, const Calibration &calib, const TemperatureCorrections &tc)
{
  int   run=0;
  int   subrun=0;
  Long64_t  timestamp; //time_t
  int   *febID = new int[numberOfFebs];
  int    nSpillsTotal=treeSpills->GetEntries();
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
  float *chi2sSummary = const_cast<float*>(chi2s[0].data());
  float *chi2sTSummary = const_cast<float*>(chi2s[1].data());
  float *errorsSummary = const_cast<float*>(errors[0].data());
  float *errorsTSummary = const_cast<float*>(errors[1].data());
  float *meanTemperaturesSummary = const_cast<float*>(meanTemperatures.data());
  float *stddevTemperaturesSummary = const_cast<float*>(stddevTemperatures.data());
  float *maxedOutFractionSummary = const_cast<float*>(maxedOutFraction.data());
  float *pedestals = const_cast<float*>(calib.GetPedestals().data());
  float *calibConstants = const_cast<float*>(calib.GetCalibrationFactors().data());
  float *calibConstantsT = const_cast<float*>(calib.GetCalibrationFactorsTemperatureCorrected().data());
  float *noiseRate = const_cast<float*>(calib.GetNoiseRate().data());
  float *xtalkProbability = const_cast<float*>(calib.GetXtalkProbability().data());
  TTree *recoTreeSummary = new TTree("runSummary","runSummary");
  recoTreeSummary->Branch("runNumber", &run, "runNumber/I");
  recoTreeSummary->Branch("subrunNumber", &subrun, "subrunNumber/I");
  recoTreeSummary->Branch("timestamp", &timestamp, "timestamp/L");
  recoTreeSummary->Branch("febID", febID, Form("febID[%i]/I",numberOfFebs));
  recoTreeSummary->Branch("spillsTotal", &nSpillsTotal, "spillsTotal/I");
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
  recoTreeSummary->Branch("chi2s", chi2sSummary, Form("chi2s[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("chi2sTemperatureCorrected", chi2sTSummary, Form("chi2sTemperatureCorrected[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("errors", errorsSummary, Form("errors[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("errorsTemperatureCorrected", errorsTSummary, Form("errorsTemperatureCorrected[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("meanTemperatures", meanTemperaturesSummary, Form("meanTemperatures[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("stddevTemperatures", stddevTemperaturesSummary, Form("stddevTemperatures[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("maxedOutFraction", maxedOutFractionSummary, Form("maxedOutFraction[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("pedestals", pedestals, Form("pedestals[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("calibConstants", calibConstants, Form("calibConstants[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("calibConstantsTemperatureCorrected", calibConstantsT, Form("calibConstantsTemperatureCorrected[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("noiseRate", noiseRate, Form("noiseRate[%i][%i]/F",numberOfFebs,channelsPerFeb));
  recoTreeSummary->Branch("xtalkProbability", xtalkProbability, Form("xtalkProbability[%i][%i]/F",numberOfFebs,channelsPerFeb));

  size_t underscorePos = runNumber.rfind('_');
  if(underscorePos==std::string::npos) run=atoi(runNumber.c_str());
  else
  {
    run=atoi(runNumber.substr(0,underscorePos).c_str());
    if(underscorePos+1<runNumber.size()) subrun=atoi(runNumber.substr(underscorePos+1).c_str());
  }

  std::ofstream txtFile;
  txtFile.open(txtFileName.c_str(),std::ios_base::trunc);

  timestamp=0;
  nSpillsActual=0;
  BoardRegisters(treeSpills, txtFile, numberOfFebs, febID, nSpillsActual, nFebSpillsActual, febTemperaturesAvg, supplyMonitorsAvg, biasVoltagesAvg, pipeline, samples, timestamp);

  txtFile<<"reference temp: "<<tc.referenceTemperatureCMB<<" deg C (CMB), "<<tc.referenceTemperatureFEB<<" deg C (FEB)   ";
  txtFile<<"overvoltage: "<<tc.overvoltageTemperatureChangeCMB<<" V/K (CMB), "<<tc.overvoltageTemperatureChangeFEB<<" V/K (FEB)   ";
  txtFile<<"calib: "<<tc.calibOvervoltageChange<<" ADC*ns/V, calibAFE: "<<tc.calibTemperatureChangeAFE<<" ADC*ns/K (FEB)   ";
  txtFile<<"PEs: "<<tc.PEOvervoltageChange<<" 1/V"<<std::endl;
  txtFile<<"  FEB  Channel     PE     PE Temp Corr    FWHM                   Signals              meanTemp    pedestal  calibConst  calibConst Temp Corr maxedOut  noise  xtalk"<<std::endl;

  std::cout<<"reference temp: "<<tc.referenceTemperatureCMB<<" deg C (CMB), "<<tc.referenceTemperatureFEB<<" deg C (FEB)   ";
  std::cout<<"overvoltage: "<<tc.overvoltageTemperatureChangeCMB<<" V/K (CMB), "<<tc.overvoltageTemperatureChangeFEB<<" V/K (FEB)   ";
  std::cout<<"calib: "<<tc.calibOvervoltageChange<<" ADC*ns/V, calibAFE: "<<tc.calibTemperatureChangeAFE<<" ADC*ns/K (FEB)   ";
  std::cout<<"PEs: "<<tc.PEOvervoltageChange<<" 1/V"<<std::endl;
  std::cout<<"         FEB  Channel     PE     PE Temp Corr    FWHM                   Signals              meanTemp    pedestal  calibConst  calibConst Temp Corr maxedOut  noise  xtalk"<<std::endl;

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
      txtFile<<std::setw(8)<<calib.GetCalibrationFactorsTemperatureCorrected()[index]<<"  ";
      txtFile<<std::setprecision(3);
      txtFile<<std::setw(8)<<maxedOutFraction[index]<<"  ";
      txtFile<<std::setw(8)<<calib.GetNoiseRate()[index]<<"  ";
      txtFile<<std::setw(8)<<calib.GetXtalkProbability()[index]<<std::endl;

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
      std::cout<<std::setw(8)<<calib.GetCalibrationFactorsTemperatureCorrected()[index]<<"  ";
      std::cout<<std::setprecision(3);
      std::cout<<std::setw(8)<<maxedOutFraction[index]<<"  ";
      std::cout<<std::setw(8)<<calib.GetNoiseRate()[index]<<"  ";
      std::cout<<std::setw(8)<<calib.GetXtalkProbability()[index]<<std::endl;
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
               const std::vector<float> mpvs[2], const std::vector<float> fwhms[2], const TemperatureCorrections &tc)
{
  std::ifstream settingsFile;
  settingsFile.open("config.txt");
  if(!settingsFile.is_open()) {std::cerr<<"Could not open config.txt."<<std::endl; exit(1);}

  float minPE, maxPE, minPEfarSide, maxPEfarSide;
  std::string settingsKey, settingsValue;
  while(settingsFile>>settingsKey>>settingsValue)
  {
    if(settingsKey=="recoMinPE")        minPE=atof(settingsValue.c_str());
    if(settingsKey=="recoMaxPE")        maxPE=atof(settingsValue.c_str());
    if(settingsKey=="recoMinPEfarSide") minPEfarSide=atof(settingsValue.c_str());
    if(settingsKey=="recoMaxPEfarSide") maxPEfarSide=atof(settingsValue.c_str());
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
  tT.AddText(Form("reference temp %.1f deg C (CMB), %.1f deg C (FEB)",tc.referenceTemperatureCMB, tc.referenceTemperatureFEB));
  tT.AddText(Form("overvoltage %.5f V/K (CMB), %.5f V/K (FEB)",tc.overvoltageTemperatureChangeCMB, tc.overvoltageTemperatureChangeFEB));
  tT.AddText(Form("calib %.1f ADC*ns/V, calibAFE %.2f ADC*ns/K (FEB), PEs %.4f 1/V", tc.calibOvervoltageChange, tc.calibTemperatureChangeAFE, tc.PEOvervoltageChange));

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

void fillDqmFile(const std::string &dqmFileName, TTree *treeMetaData, TTree *recoTree, TTree *recoTreeSpill, TTree *recoTreeSummary)
{
  int         run;
  int         subrun;
  std::string *configuration = new std::string;
  treeMetaData->SetBranchAddress("runNumber", &run);
  treeMetaData->SetBranchAddress("subrunNumber", &subrun);
  treeMetaData->SetBranchAddress("configuration", &configuration);
  treeMetaData->GetEntry(0); //only one entry per file

  int       nEventsExpected;
  int       nEventsActual;
  Long64_t  timestampTmp=0;  //time_t
  recoTreeSpill->SetBranchAddress("spill_nevents", &nEventsExpected);
  recoTreeSpill->SetBranchAddress("spill_neventsActual", &nEventsActual);
  recoTreeSpill->SetBranchAddress("spill_timestamp", &timestampTmp);

  long totalEventsExpected=0;
  long totalEventsActual=0;
  Long64_t firstTimestamp=0;
  Long64_t lastTimestamp=0;
  for(int i=0; i<recoTreeSpill->GetEntries(); ++i)
  {
    recoTreeSpill->GetEntry(i);
    if(firstTimestamp==0) firstTimestamp=timestampTmp;  //the first spill may not have a time stamp
    if(timestampTmp!=0) lastTimestamp=timestampTmp;  //the last spill may not have a time stamp

    totalEventsExpected+=nEventsExpected;
    totalEventsActual+=nEventsActual;
  }

  //dqm data in config file
  std::ifstream configFile;
  configFile.open("config.txt");
  if(!configFile.is_open()) {std::cerr<<"Could not open config.txt."<<std::endl; exit(1);}

  std::map<std::string,float> metaDataFValues;
  std::map<std::string,int> metaDataIValues;
  std::string configKey, configValue;
  while(configFile>>configKey>>configValue)
  {
    if(configKey.compare(0,9,"metaDataF")==0) metaDataFValues[configKey.substr(10)]=atof(configValue.c_str());
    if(configKey.compare(0,9,"metaDataI")==0) metaDataIValues[configKey.substr(10)]=atoi(configValue.c_str());
  }
  configFile.close();

  TFile dqmFile(dqmFileName.c_str(), "RECREATE");

  //meta data in dqm file
  TTree dqmTreeMetaData("metaData","metaData");
  dqmTreeMetaData.Branch("runNumber", &run);
  dqmTreeMetaData.Branch("subrunNumber", &subrun);
  dqmTreeMetaData.Branch("configuration", configuration);
  dqmTreeMetaData.Branch("firstSpillTime", &firstTimestamp);
  dqmTreeMetaData.Branch("lastSpillTime", &lastTimestamp);
  dqmTreeMetaData.Branch("eventsExpected", &totalEventsExpected);
  dqmTreeMetaData.Branch("eventsActual", &totalEventsActual);
  for(auto i=metaDataFValues.begin(); i!=metaDataFValues.end(); ++i)
  {
    dqmTreeMetaData.Branch(i->first.c_str(), &i->second);
  }
  for(auto i=metaDataIValues.begin(); i!=metaDataIValues.end(); ++i)
  {
    dqmTreeMetaData.Branch(i->first.c_str(), &i->second);
  }

  TH1F *hdqmPEs = new TH1F("PEs","PEs;PE yield [PEs];count",150,0,150);
  TH1F *hdqmPEsTemperatureCorrected = new TH1F("PEsTemperatureCorrected","PEsTemperatureCorrected;PE yield [PEs];count",150,0,150);
  TH1F *hdqmTime = new TH1F("time","time;time [ns];count",100,700,1700);
  TH1F *hdqmMpvPEs = new TH1F("mpvPEs","mpvPEs;PE yield [PEs];count",150,0,150);
  TH1F *hdqmMpvPEsTemperatureCorrected = new TH1F("mpvPEsTemperatureCorrected","mpvPEsTemperatureCorrected;PE yield [PEs];count",150,0,150);
  TH1F *hdqmPedestals = new TH1F("pedestals","pedestals;pedestal [ADC];count",100,-50,50);
  TH1F *hdqmCalibConstants = new TH1F("calibConstants","calibConstants;calib const [ADC*ns/PEs];count",100,200,700);
  TH1F *hdqmCalibConstantsTemperatureCorrected = new TH1F("calibConstantsTemperatureCorrected","calibConstantsTemperatureCorrected;calib const [ADC*ns/PEs];count",100,200,700);
  TH1F *hdqmMaxedOutFraction = new TH1F("maxedOutFraction","maxedOutFraction;fraction;count",100,0,0.01);
  TH1F *hdqmNoiseRate = new TH1F("noiseRate","noiseRate;rate [MHz];count",100,0,0.4);
  TH1F *hdqmXtalkProbability = new TH1F("xtalkProbability","xtalkProbability;probability;count",100,0,0.2);
  TH1F *hdqmMeanTemperatures = new TH1F("meanTemperatures","meanTemperatures;temp [deg C];count",40,0,40);
  TH1F *hdqmFebTemperaturesAvg = new TH1F("febTemperaturesAvg","febTemperaturesAvg;temp [deg C];count",60,10,70);
  TH1F *hdqmBiasVoltagesAvg = new TH1F("biasVoltagesAvg","biasVoltagesAvg;bias [V];count",100,50,60);
  TH1I *hdqmFebID = new TH1I("febID","febID;febID;count",100,0,100);
  recoTree->Draw("PEs>>+PEs");
  recoTree->Draw("PEsTemperatureCorrected>>+PEsTemperatureCorrected");
  recoTree->Draw("time>>+time");

  //fit time peak
  int maxBin=hdqmTime->GetMaximumBin();
  float maxTime=hdqmTime->GetBinCenter(maxBin);
  hdqmTime->Fit("gaus","S","",maxTime-100.0,maxTime+100.0);

  recoTreeSummary->Draw("PEs>>+mpvPEs");
  recoTreeSummary->Draw("PEsTemperatureCorrected>>+mpvPEsTemperatureCorrected");
  recoTreeSummary->Draw("pedestals>>+pedestals");
  recoTreeSummary->Draw("calibConstants>>+calibConstants");
  recoTreeSummary->Draw("calibConstantsTemperatureCorrected>>+calibConstantsTemperatureCorrected");
  recoTreeSummary->Draw("maxedOutFraction>>+maxedOutFraction");
  recoTreeSummary->Draw("noiseRate>>+noiseRate");
  recoTreeSummary->Draw("xtalkProbability>>+xtalkProbability");
  recoTreeSummary->Draw("meanTemperatures>>+meanTemperatures");
  recoTreeSummary->Draw("febTemperaturesAvg>>+febTemperaturesAvg");
  recoTreeSummary->Draw("biasVoltagesAvg>>+biasVoltagesAvg");
  recoTreeSummary->Draw("febID>>+febID");
  hdqmPEs->Write();
  hdqmPEsTemperatureCorrected->Write();
  hdqmTime->Write();
  hdqmMpvPEs->Write();
  hdqmMpvPEsTemperatureCorrected->Write();
  hdqmPedestals->Write();
  hdqmCalibConstants->Write();
  hdqmCalibConstantsTemperatureCorrected->Write();
  hdqmMaxedOutFraction->Write();
  hdqmNoiseRate->Write();
  hdqmXtalkProbability->Write();
  hdqmMeanTemperatures->Write();
  hdqmFebTemperaturesAvg->Write();
  hdqmBiasVoltagesAvg->Write();
  hdqmFebID->Write();

  dqmTreeMetaData.Fill();
  dqmTreeMetaData.Write();

  dqmFile.Close();
}

void process(const std::string &runNumber, const std::string &inFileName, const std::string &calibFileName, const std::string &recoFileName, const std::string &recoFileName2,
             const std::string &pdfFileName, const std::string &txtFileName, const std::string &dqmFileName,
             bool usePoisson, const TemperatureCorrections &tc, const std::string &channelMapFile, float PEcut)
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

  recoFile.mkdir("plots");
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
  CrvEvent event(runNumber, numberOfFebs, channelsPerFeb, numberOfSamples, tree, recoTree, tc, PEcut, channelMapFile);

  int nEvents = tree->GetEntries();
//std::cout<<"USING A WRONG NUMBER OF EVENTS"<<std::endl;
//if(nEvents>100) nEvents=100;

  gSystem->Unlink(pdfFileName.c_str());
  TCanvas c0;
  c0.Print(Form("%s[", pdfFileName.c_str()), "pdf");

  for(int i=0; i<nEvents; i++) event.Reconstruct(i, calib);

  std::vector<float> mpvs[2];
  std::vector<float> fwhms[2];
  std::vector<float> signals[2];
  std::vector<float> chi2s[2];
  std::vector<float> errors[2];
  std::vector<float> meanTemperatures;
  std::vector<float> stddevTemperatures;
  std::vector<float> maxedOutFraction;
  for(int feb=0; feb<numberOfFebs; feb++)
  {
    for(int channel=0; channel<channelsPerFeb; channel++)
    {
      maxedOutFraction.push_back(((float)event.GetMaxedOutEvents(feb,channel))/nEvents);

      TH1F      *h[2];
      TPaveText *t[2];
      for(int i=0; i<2; ++i)
      {
        event.GetCanvas(feb, channel)->cd(i+1);
//        gPad->SetLogy(1);
        h[i]=event.GetHistPEs(i,feb,channel);

        float mpv=0;
        float fwhm=0;
        float nsignals=0;
        float chi2=0;
        float error=0;
        if(!usePoisson) LandauGauss(*h[i], mpv, fwhm, nsignals, chi2, error);
        else Poisson(*h[i], mpv, fwhm, nsignals, chi2, error);
        mpvs[i].push_back(mpv);
        fwhms[i].push_back(fwhm);
        signals[i].push_back(nsignals);
        chi2s[i].push_back(chi2);
        errors[i].push_back(error);
        if(!usePoisson) h[i]->GetXaxis()->SetRangeUser(10,150);
        h[i]->Draw();
        TDirectory *tempDirectory = gDirectory;
        gDirectory->cd("plots");
        h[i]->Write();  //write to root file
        tempDirectory->cd();

        t[i]=new TPaveText(0.65,0.58,0.88,0.72,"NDC");
        t[i]->SetFillColor(kWhite);
        t[i]->SetTextColor(kRed);
        t[i]->SetTextAlign(12);
        t[i]->AddText(Form("MPV: %.1f PEs",mpv));
        t[i]->AddText(Form("FWHM: %.1f PEs",fwhm));
        t[i]->AddText(Form("Signals: %.0f",nsignals));
        t[i]->AddText(Form("Chi2/Ndf: %.3f",chi2));
        t[i]->AddText(Form("maxed out: %.3f",maxedOutFraction.back()));
        t[i]->Draw("same");
      }

      event.GetCanvas(feb, channel)->cd(2);
      TPaveText tt(0.48,0.72,0.88,0.88,"NDC");
      tt.SetFillColor(kWhite);
      tt.SetTextColor(kBlack);
      tt.SetTextAlign(12);
      tt.AddText("Temperature corrected PE values");
      tt.AddText(Form("reference temp %.1f deg C (CMB), %.1f deg C (FEB)",tc.referenceTemperatureCMB, tc.referenceTemperatureFEB));
      tt.AddText(Form("overvoltage %.4f V/K (CMB), %.4f V/K (FEB)",tc.overvoltageTemperatureChangeCMB, tc.overvoltageTemperatureChangeFEB));
      tt.AddText(Form("calib %.1f ADC*ns/V, calibAFE %.2f ADC*ns/K (FEB), PEs %.4f 1/V", tc.calibOvervoltageChange, tc.calibTemperatureChangeAFE, tc.PEOvervoltageChange));
      tt.Draw("same");

      TGraph *g=event.GetHistTemperatures(feb, channel);
      if(g->GetN()>0)
      {
        event.GetCanvas(feb, channel)->cd(4);
        gPad->cd(1);
        g->Draw("AP");
        meanTemperatures.push_back(g->GetMean(2));
        stddevTemperatures.push_back(g->GetRMS(2));
      }
      else
      {
        meanTemperatures.push_back(0);
        stddevTemperatures.push_back(0);
      }

      TGraph *g2=event.GetHistTemperaturesFEB(feb, channel);
      if(g2->GetN()>0)
      {
        event.GetCanvas(feb, channel)->cd(4);
        gPad->cd(2);
        g2->Draw("AP");
      }

      event.GetCanvas(feb,channel)->Update();
      event.GetCanvas(feb,channel)->Print(pdfFileName.c_str(), "pdf");
    }
  }

  StorePEyields(runNumber, txtFileName, numberOfFebs, channelsPerFeb, mpvs, fwhms, signals, chi2s, errors, meanTemperatures, stddevTemperatures, maxedOutFraction, treeSpills, nEvents, calib, tc);

  Summarize(pdfFileName, txtFileName, numberOfFebs, channelsPerFeb, mpvs, fwhms, tc);

  c0.Print(Form("%s]", pdfFileName.c_str()), "pdf");

  recoTree->Write("", TObject::kOverwrite);
  recoTreeSpill->Write("", TObject::kOverwrite);

  //store a copy of the reco file without the adc
  recoTree->SetBranchStatus("adc",0);
  gDirectory->cd("plots");
  TDirectory *plotDirectory = gDirectory;

  TFile recoFile2(recoFileName2.c_str(), "RECREATE");
  TTree *recoTree2 = recoTree->CloneTree();
  recoTree2->Write("", TObject::kOverwrite);
  TTree *recoTreeSpill2 = recoTreeSpill->CloneTree();
  recoTreeSpill2->Write("", TObject::kOverwrite);

  TTree *recoTreeSummary = (TTree*)recoFile.Get("runSummary");
  TTree *recoTreeSummary2 = recoTreeSummary->CloneTree();
  recoTreeSummary2->Write("", TObject::kOverwrite);

  TDirectory *plotDirectory2 = recoFile2.mkdir("plots");
  plotDirectory2->cd();
  size_t nHistograms = plotDirectory->GetListOfKeys()->GetSize();
  for(size_t i=0; i<nHistograms; ++i)
  {
    TObject *hist = ((TKey*)plotDirectory->GetListOfKeys()->At(i))->ReadObj();
    hist->Write();
    delete hist;
  }

  //meta data
  TTree *treeMetaData = (TTree*)file.Get("metaData");
  fillDqmFile(dqmFileName, treeMetaData, recoTree, recoTreeSpill, recoTreeSummary);

  recoFile.Close();
  recoFile2.Close();
  file.Close();
}

void makeFileNames(const std::string &runNumber, std::string &inFileName, std::string &calibFileName, std::string &recoFileName, std::string &recoFileName2, std::string &pdfFileName, std::string &txtFileName, std::string &dqmFileName)
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
    recoFileName2 = recoDirName+"crv.reco2."+dirEntry.path().stem().string().substr(s0.length())+".root";
    pdfFileName = recoDirName+"log.crv.reco."+dirEntry.path().stem().string().substr(s0.length())+".pdf";
    txtFileName = recoDirName+"log.crv.reco."+dirEntry.path().stem().string().substr(s0.length())+".txt";
    dqmFileName = recoDirName+"crv.dqm."+dirEntry.path().stem().string().substr(s0.length())+".root";
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
      recoFileName2 = recoDirName+"rec2.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".root";
      pdfFileName = recoDirName+"rec.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".pdf";
      txtFileName = recoDirName+"rec.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".txt";
      dqmFileName = recoDirName+"dqm.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".root";
      break;
    }
  }

  if(!found) {std::cerr<<"Could not open input file for run "<<runNumber<<"."<<std::endl; exit(1);}
}

void printHelp()
{
  std::cout<<"Use as"<<std::endl;
  std::cout<<"recoCrv -h                        Prints this help."<<std::endl;
  std::cout<<"recoCrv RUNNUMBER [-p] [options] Reconstructs a run."<<std::endl;
  std::cout<<"-p                                Uses Poisson fit for the PE yield distribution, instead of a convoluted Gauss-Landau fit."<<std::endl;
  std::cout<<"                                  (used for LED light sources)"<<std::endl;
  std::cout<<"--PEOvervoltageChange ###         change of PE yield per overvoltage change (percentage change of the PE yield per V)"<<std::endl;
  std::cout<<"                                  default: 0.229 1/V"<<std::endl;
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
  std::cout<<"--channelMap ###                  file name of the channel map to be used to fit tracks"<<std::endl;
  std::cout<<"                                  default: no file (=no track fits)"<<std::endl;
  std::cout<<"--PEcut ###                       mininum number of PEs required to be used for the track fit"<<std::endl;
  std::cout<<"                                  default: 5 PEs"<<std::endl;
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
  std::string channelMapFile;
  bool usePoisson = false;
  TemperatureCorrections tc;  //default values are set in the definition of this struct
  float PEcut = 5;
  for(int i=2; i<argc; i++)  //loop through options that don't require a second argument
  {
    if(strcmp(argv[i],"-p")==0) usePoisson=true;
  }
  for(int i=2; i<argc-1; i++)  //loop through options that require a second argument
  {
    if(strcmp(argv[i],"--PEOvervoltageChange")==0)       tc.PEOvervoltageChange=atof(argv[i+1]);
    if(strcmp(argv[i],"--calibOvervoltageChange")==0)    tc.calibOvervoltageChange=atof(argv[i+1]);
    if(strcmp(argv[i],"--calibTempChangeAFE")==0)        tc.calibTemperatureChangeAFE=atof(argv[i+1]);
    if(strcmp(argv[i],"--overvoltageTempChangeCMB")==0)  tc.overvoltageTemperatureChangeCMB=atof(argv[i+1]);
    if(strcmp(argv[i],"--overvoltageTempChangeFEB")==0)  tc.overvoltageTemperatureChangeFEB=atof(argv[i+1]);
    if(strcmp(argv[i],"--referenceTempCMB")==0)          tc.referenceTemperatureCMB=atof(argv[i+1]);
    if(strcmp(argv[i],"--referenceTempFEB")==0)          tc.referenceTemperatureFEB=atof(argv[i+1]);
    if(strcmp(argv[i],"--channelMap")==0)                channelMapFile=argv[i+1];
    if(strcmp(argv[i],"--PEcut")==0)                     PEcut=atof(argv[i+1]);
  }
  std::string inFileName;
  std::string calibFileName;
  std::string recoFileName;
  std::string recoFileName2;
  std::string pdfFileName;
  std::string txtFileName;
  std::string dqmFileName;
  makeFileNames(runNumber, inFileName, calibFileName, recoFileName, recoFileName2, pdfFileName, txtFileName, dqmFileName);

  process(runNumber, inFileName, calibFileName, recoFileName, recoFileName2, pdfFileName, txtFileName, dqmFileName, usePoisson, tc, channelMapFile, PEcut);

  return 0;
}
