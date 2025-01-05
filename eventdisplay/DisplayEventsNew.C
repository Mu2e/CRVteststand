#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGResourcePool.h>
#include <TBox.h>
#include <vector>
#include <map>
#include <iostream>
#include <cstdlib>
#include <CLHEP/Vector/ThreeVector.h>

class CrvCounter : public TWbox
{
  public:
  CrvCounter(int side, int sector, int ignoreForFit, int feb, int channel1, int channel2, TText *textPEs,
             TPad *padAdc, TPad *padReco1, TPad *padReco2, double x, double y, double xDisplay, double yDisplay) :
             TWbox(xDisplay-25.65, yDisplay-9.9, xDisplay+25.65, yDisplay+9.9, 0,(feb>=0?3:1),1),
             _side(side), _sector(sector), _ignoreForFit(ignoreForFit), _feb(feb), _channel1(channel1), _channel2(channel2), _PE1(0), _PE2(0),
             _x(x), _y(y), _textPEs(textPEs), _padAdc(padAdc)
  {
    _padReco[0]=padReco1;
    _padReco[1]=padReco2;
    _fitStatus[0]=0;
    _fitStatus[1]=0;
  }
  void GetChannels(int &feb, int &channel1, int &channel2)
  {
    feb=_feb; channel1=_channel1; channel2=_channel2;
  }
  void GetPos(float &x, float &y)
  {
    x=_x; y=_y;
  }
  int GetSector()
  {
    return _sector;
  }
  int IgnoreForFit()
  {
    return _ignoreForFit;
  }
  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py)
  {
    TWbox::ExecuteEvent(event, px, py);
    if(event==kMouseMotion)
    {
      _padReco[0]->GetCanvas()->cd();
      _padReco[0]->GetCanvas()->Modified();
      _padReco[0]->GetCanvas()->Update();
      _padReco[1]->GetCanvas()->cd();
      _padReco[1]->GetCanvas()->Modified();
      _padReco[1]->GetCanvas()->Update();
    }
    if(event==kButton1Down)
    {
      //remove old things
      _padAdc->cd();
      _padAdc->Clear();
      for(int i=0; i<2; i++)
      {
        _padReco[i]->cd();
        _padReco[i]->Clear();
      }

      //PEs
      if(_feb>=0) _textPEs->SetTitle(Form("FEB %i   channels %i/%i   PEs %i/%i",_feb,_channel1,_channel2,_PE1,_PE2));
      else _textPEs->SetTitle("no FEB connected to this counter");

      if(_feb<0 || _adc[0].empty() || _adc[1].empty())
      {
        gPad->GetCanvas()->Modified();
        gPad->GetCanvas()->Update();
        return;
      }

      //waveforms
      _padAdc->cd();
      TMultiGraph *mg = new TMultiGraph(_side==0?"adcGraphA":"adcGraphB",";Time [ns];ADC");
      for(int i=0; i<2; i++)
      {
        TGraph *g = new TGraph(_adc[i].size());
        for(int sample=0; sample<_adc[i].size(); ++sample)
        {
          double t=sample*12.55;
          double v=_adc[i][sample];
          g->SetPoint(sample, t, v);
        }
        g->SetMarkerStyle(20);
        g->SetMarkerColor(i+2);
        mg->Add(g);
      }
      mg->Draw("AP");
      mg->GetYaxis()->SetTitleOffset(0.7);

      TText *channel1 = new TText(0.2,0.7,Form("channel %i",_channel1));
      TText *channel2 = new TText(0.2,0.8,Form("channel %i",_channel2));
      channel1->SetNDC(true);
      channel2->SetNDC(true);
      channel1->SetTextColor(2);
      channel2->SetTextColor(3);
      channel1->Draw();
      channel2->Draw();

      //reco fits
      int firstRecoBin=std::min(_recoStartBin[0],_recoStartBin[1]);
      int lastRecoBin=std::max(_recoEndBin[0],_recoEndBin[1]);
      unsigned int n = lastRecoBin-firstRecoBin+1;
      for(int i=0; i<2; i++)
      {
        _padReco[i]->cd();

        if(_fitStatus[i]!=1) continue;

//        unsigned int n = _recoEndBin[i]-_recoStartBin[i]+1;
        if(n==0) continue;
        double *t = new double[n];
        double *v = new double[n];
/*
        for(unsigned int j=0; j<n; j++)
        {
          t[j]=(_recoStartBin[i]+j)*12.55;
          v[j]=_adc[i][_recoStartBin[i]+j];
        }
*/
        for(unsigned int j=0; j<n; j++)
        {
          t[j]=(firstRecoBin+j)*12.55;
          v[j]=_adc[i][firstRecoBin+j];
        }
        TGraph *g = new TGraph(n,t,v);
        g->SetTitle(";time [ns];ADC");
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.5);
        g->SetMarkerColor(i+2);
        g->Draw("AP");
        g->GetXaxis()->SetLabelSize(0.07);
        g->GetYaxis()->SetLabelSize(0.07);
        g->GetXaxis()->SetNdivisions(505);

        delete[] t;
        delete[] v;

        double tF1=_recoStartBin[i]*12.55;
        double tF2=_recoEndBin[i]*12.55;
        int nF=(tF2-tF1)/1.0 + 1;
        double *tF = new double[nF];
        double *vF = new double[nF];
        for(int jF=0; jF<nF; jF++)
        {
          double p0 = _pulseHeight[i]*2.718;
          double p1 = _time[i];
          double p2 = _beta[i];
          tF[jF] = tF1 + jF*1.0;
          vF[jF] = p0*TMath::Exp(-(tF[jF]-p1)/p2-TMath::Exp(-(tF[jF]-p1)/p2));
          vF[jF]+= _pedestal[i];
          if(std::isnan(vF[jF])) nF=0;
        }
        if(nF>0)
        {
          TGraph *gF=new TGraph();
          gF->SetLineWidth(2);
          gF->SetLineColor(kBlack);
          gF->DrawGraph(nF,tF,vF,"sameL");
        }

        delete[] tF;
        delete[] vF;
      }

      gPad->GetCanvas()->Modified();
      gPad->GetCanvas()->Update();
    }
  }
  void SetPEs(int PE1, int PE2)
  {
    _PE1=PE1; _PE2=PE2;
  }
//  void SetADCs(int adc1[], int adc2[])
  void SetADCs(short adc1[], short adc2[])
  {
    _adc[0].clear();
    _adc[1].clear();
    for(int i=0; i<127; i++)
    {
      _adc[0].push_back(adc1[i]);
      _adc[1].push_back(adc2[i]);
    }
  }
  void SetRecoValues(int fitStatus1, int fitStatus2,
                     int recoStartBin1, int recoStartBin2,
                     int recoEndBin1, int recoEndBin2,
                     double time1, double time2,
                     double beta1, double beta2,
                     double pulseHeight1, double pulseHeight2,
                     double pedestal1, double pedestal2)
  {
    _fitStatus[0]=fitStatus1;
    _recoStartBin[0]=recoStartBin1;
    _recoEndBin[0]=recoEndBin1;
    _time[0]=time1;
    _beta[0]=beta1;
    _pulseHeight[0]=pulseHeight1;
    _pedestal[0]=pedestal1;

    _fitStatus[1]=fitStatus2;
    _recoStartBin[1]=recoStartBin2;
    _recoEndBin[1]=recoEndBin2;
    _time[1]=time2;
    _beta[1]=beta2;
    _pulseHeight[1]=pulseHeight2;
    _pedestal[1]=pedestal2;
  }

  private:
  int _side;
  int _sector;
  int _ignoreForFit;
  int _feb;
  int _channel1, _channel2;
  int _PE1, _PE2;
  float _x, _y;
  std::vector<int> _adc[2];
  int              _fitStatus[2];
  int              _recoStartBin[2], _recoEndBin[2];
  double           _time[2], _beta[2], _pulseHeight[2], _pedestal[2];
  TText *_textPEs;
  TPad  *_padAdc;
  TPad  *_padReco[2];
};

class EventWindow : public TGMainFrame
{
  public:
  EventWindow(const TGWindow *p, UInt_t w, UInt_t h, const std::string &runFileName, const std::string &channelMapFileName);
  virtual ~EventWindow();
  void DoGoToEvent();
  void DoPrevEvent();
  void DoNextEvent();
  void DoSave();
  void DoExit();
  bool DoEvent();
  void DrawFront(int side);
  void ReadChannelMap();

  private:
  std::string         _runFileName;
  std::string         _runNumber;
  TFile*              _file;
  TTree*              _tree;
  int                 _nEntries;
  int                 _currentEntry;
  int                 _currentSpill;
  int                 _currentEvent;
  int                 _numberOfFebs;
  int                 _channelsPerFeb;
  int                 _numberOfSamples;

  std::map<int,float>  _minX;
  std::map<int,float>  _maxX;
  std::map<int,float>  _minY;
  std::map<int,float>  _maxY;
  float  _centerX;

  TPad*                     _pads[2];
  TPad*                     _padsAdc[2];
  TPad*                     _padsReco1[2];
  TPad*                     _padsReco2[2];
  std::vector<CrvCounter*>  _boxes[2];
  std::map<int,TLine*>      _fit[2];  //sorted by sector (entry 0 is used for all sectors)
  TText*                    _textPEs[2];
  TText*                    _textTotalPEs[2];

  TGNumberEntryField* _gotoSpill;
  TGNumberEntryField* _gotoEvent;
  TGLabel*            _textRunNumber;
  TGLabel*            _textSpillNumber;
  TGLabel*            _textEventNumber;
  TGCheckButton*      _fitButton;

  struct channelStruct
  {
    int _side;
    int _sector;
    float _x,_y;
    int _ignoreForFit;  //for the fit using all sectors
    channelStruct() : _side(0), _sector(0), _x(0), _y(0), _ignoreForFit(0) {}
    channelStruct(int side, int sector, float x, float y, int ignoreForFit) : _side(side), _sector(sector), _x(x), _y(y), _ignoreForFit(ignoreForFit) {}
  };

  std::string                                                  _channelMapFileName;
  std::map<std::pair<int,std::pair<int,int> >, channelStruct>  _channelMap;   //feb,(channel1,channel2) --> channelStruct

  ClassDef(EventWindow, 0)
};

void EventWindow::DoGoToEvent()
{
  int newSpill = _gotoSpill->GetNumber();
  int newEvent = _gotoEvent->GetNumber();

  //find tree enty of newSpill/newEvent
  int spillNumber;
  int eventNumber;
  _tree->SetBranchStatus("*",0);
  _tree->SetBranchStatus("spillNumber",1);
  _tree->SetBranchStatus("eventNumber",1);
  _tree->SetBranchAddress("spillNumber", &spillNumber);
  _tree->SetBranchAddress("eventNumber", &eventNumber);
  for(int i=0; i<_nEntries; ++i)
  {
    _tree->GetEvent(i);
    if(newSpill==spillNumber && newEvent==eventNumber)
    {
      _currentEntry=i;
      break;
    }
  }
  _tree->SetBranchStatus("*",1);

  DoEvent();  //resets the spill/event number, if spill/event wasn't found
}

void EventWindow::DoPrevEvent()
{
  while(1) //while loop and bool return value of DoEvent may become useful, if an event filter gets implemented
  {
    _currentEntry--;
    if(_currentEntry<0) _currentEntry=_nEntries-1;
    if(DoEvent()) break;
  }
}

void EventWindow::DoNextEvent()
{
  while(1) //while loop and bool return value of DoEvent may become useful, if an event filter gets implemented
  {
    _currentEntry++;
    if(_currentEntry>=_nEntries) _currentEntry=0;
    if(DoEvent()) break;
  }
}

void EventWindow::DoSave()
{
  SaveAs(Form("Run%s_Spill%i_Event%i.png",_runNumber.c_str(),_currentSpill,_currentEvent));
}

void EventWindow::ReadChannelMap()
{
  std::ifstream file(_channelMapFileName);
  if(!file.is_open()) {std::cerr<<"Channel map file could not be opened."<<std::endl; exit(1);}

  string header;
  getline(file,header);
  bool hasSectors=false;
  if(header.find("sector")!=string::npos) hasSectors=true;
  bool hasIgnoreForFit=false;
  if(header.find("ignoreForFit")!=string::npos) hasIgnoreForFit=true;

  int febA, channelA1, channelA2;
  int febB, channelB1, channelB2;
  float x, y;
  int sector=0;
  int ignoreForFit=0;

  while(file >> febA >> channelA1 >> channelA2 >> febB >> channelB1 >> channelB2 >> x >> y)
  {
    if(hasSectors) file >> sector;
    if(hasIgnoreForFit) file >> ignoreForFit;
    std::pair<int,int> channelPairA(channelA1,channelA2);
    std::pair<int,int> channelPairB(channelB1,channelB2);
    std::pair<int,std::pair<int,int> > counterA(febA,channelPairA);
    std::pair<int,std::pair<int,int> > counterB(febB,channelPairB);

    _channelMap[counterA]=channelStruct(0,sector,x,y,ignoreForFit);
    _channelMap[counterB]=channelStruct(1,sector,x,y,ignoreForFit);

    //corners of display
    if(_minX.find(-1)==_minX.end())
    {
      _minX[-1]=x;
      _maxX[-1]=x;
      _minY[-1]=y;
      _maxY[-1]=y;
    }

    if(x<_minX[-1]) _minX[-1]=x;
    if(x>_maxX[-1]) _maxX[-1]=x;
    if(y<_minY[-1]) _minY[-1]=y;
    if(y>_maxY[-1]) _maxY[-1]=y;

    //area of fit through all modules
    if(ignoreForFit==0)
    {
      if(_minX.find(0)==_minX.end())
      {
        _minX[0]=x;
        _maxX[0]=x;
        _minY[0]=y;
        _maxY[0]=y;
      }

      if(x<_minX[0]) _minX[0]=x;
      if(x>_maxX[0]) _maxX[0]=x;
      if(y<_minY[0]) _minY[0]=y;
      if(y>_maxY[0]) _maxY[0]=y;
    }

    //area of fit through individual sectors
    if(_minX.find(sector)==_minX.end())
    {
      _minX[sector]=x;
      _maxX[sector]=x;
      _minY[sector]=y;
      _maxY[sector]=y;
    }

    if(x<_minX[sector]) _minX[sector]=x;
    if(x>_maxX[sector]) _maxX[sector]=x;
    if(y<_minY[sector]) _minY[sector]=y;
    if(y>_maxY[sector]) _maxY[sector]=y;
  }

  file.close();
}

bool EventWindow::DoEvent()
{
  int     spillNumber;
  int     eventNumber;
  float   PEs[_numberOfFebs][_channelsPerFeb];
  int     fitStatus[_numberOfFebs][_channelsPerFeb];
//  int     adc[_numberOfFebs][_channelsPerFeb][_numberOfSamples];
  short   adc[_numberOfFebs][_channelsPerFeb][_numberOfSamples];
  float   time[_numberOfFebs][_channelsPerFeb];
  float   beta[_numberOfFebs][_channelsPerFeb];
  float   pulseHeight[_numberOfFebs][_channelsPerFeb];
  int     recoStartBin[_numberOfFebs][_channelsPerFeb];
  int     recoEndBin[_numberOfFebs][_channelsPerFeb];
  float   pedestal[_numberOfFebs][_channelsPerFeb];

  _tree->SetBranchAddress("spillNumber", &spillNumber);
  _tree->SetBranchAddress("eventNumber", &eventNumber);
  _tree->SetBranchAddress("PEs", PEs);
  _tree->SetBranchAddress("PEs", PEs);
  _tree->SetBranchAddress("fitStatus", fitStatus);
  _tree->SetBranchAddress("adc", adc);
  _tree->SetBranchAddress("time", time);
  _tree->SetBranchAddress("beta", beta);
  _tree->SetBranchAddress("pulseHeight", pulseHeight);
  _tree->SetBranchAddress("recoStartBin", recoStartBin);
  _tree->SetBranchAddress("recoEndBin", recoEndBin);
  _tree->SetBranchAddress("pedestal", pedestal);
  _tree->GetEvent(_currentEntry);

  _textSpillNumber->ChangeText(Form("Spill %05i",spillNumber));
  _textEventNumber->ChangeText(Form("Event %05i",eventNumber));
  _gotoSpill->SetNumber(spillNumber);
  _gotoEvent->SetNumber(eventNumber);
  _textPEs[0]->SetTitle("");
  _textPEs[1]->SetTitle("");

  float totalPEs[2]={0,0};

  //track fit variables, sorted by sector (entry of 0 uses all sectors)
  std::map<int,float> sumX;
  std::map<int,float> sumY;
  std::map<int,float> sumXY;
  std::map<int,float> sumYY;
  std::map<int,float> sumPEs;
  std::map<int,int>   nPoints;

  sumX[0]=0;
  sumY[0]=0;
  sumXY[0]=0;
  sumYY[0]=0;
  sumPEs[0]=0;
  nPoints[0]=0;

  for(int side=0; side<2; ++side)
  {
    for(auto boxIter=_boxes[side].begin(); boxIter!=_boxes[side].end(); ++boxIter)
    {
      CrvCounter *counter=*boxIter;

      int feb, channel1, channel2;
      counter->GetChannels(feb,channel1,channel2);
      if(feb<0 || feb>=_numberOfFebs) continue;  //feb not in event tree
      if(channel1<0 || channel2<0 || channel1>=_channelsPerFeb || channel2>=_channelsPerFeb) continue;  //channels not in event tree

      float PE1     = PEs[feb][channel1];
      float PE2     = PEs[feb][channel2];
      if(PE1<=0 || fitStatus[feb][channel1]==0) PE1=0;
      if(PE2<=0 || fitStatus[feb][channel2]==0) PE2=0;

      float PE  = PE1+PE2;
      int color = 100.0*(PE/100.0)+2000.0;
      if(PE<=0) color=0;
      if(color>=2100) color=2099;

      //collect information for the fit
      //uses hit information from both sides
      if(PE>=5)
      {
        float x,y;
        counter->GetPos(x,y);

        if(counter->IgnoreForFit()==0)
        {
          sumX[0] +=x*PE;
          sumY[0] +=y*PE;
          sumXY[0]+=x*y*PE;
          sumYY[0]+=y*y*PE;
          sumPEs[0]+=PE;
          ++nPoints[0];
        }

        int sector=counter->GetSector();
        if(sector!=0)
        {
          if(nPoints.find(sector)==nPoints.end())
          {
            sumX[sector]=0;
            sumY[sector]=0;
            sumXY[sector]=0;
            sumYY[sector]=0;
            sumPEs[sector]=0;
            nPoints[sector]=0;
          }
          sumX[sector] +=x*PE;
          sumY[sector] +=y*PE;
          sumXY[sector]+=x*y*PE;
          sumYY[sector]+=y*y*PE;
          sumPEs[sector]+=PE;
          ++nPoints[sector];
        }
      }

      _pads[side]->cd();
      counter->SetFillColor(color);
      counter->SetPEs(PE1,PE2);
      counter->SetADCs(adc[feb][channel1],adc[feb][channel2]);
      counter->SetRecoValues(fitStatus[feb][channel1], fitStatus[feb][channel2],
                             recoStartBin[feb][channel1],recoStartBin[feb][channel2],
                             recoEndBin[feb][channel1],recoEndBin[feb][channel2],
                             time[feb][channel1], time[feb][channel2],
                             beta[feb][channel1], beta[feb][channel2],
                             pulseHeight[feb][channel1], pulseHeight[feb][channel2],
                             pedestal[feb][channel1], pedestal[feb][channel2]);

      totalPEs[side]+=PE;
    }
  }

  _textTotalPEs[0]->SetTitle(Form("Side A - Total number of PEs: %.0f",totalPEs[0]));
  _textTotalPEs[1]->SetTitle(Form("Side B - Total number of PEs: %.0f",totalPEs[1]));

//ADC part
  for(int i=0; i<2; i++)
  {
    _padsAdc[i]->Clear();
    _padsReco1[i]->Clear();
    _padsReco2[i]->Clear();
  }

//fit part
  for(int side=0; side<2; ++side)
  {
    _pads[side]->cd();
    //delete old fit lines
    for(auto oldFitIter=_fit[side].begin(); oldFitIter!=_fit[side].end(); ++oldFitIter) {delete oldFitIter->second;}
    _fit[side].clear();
  }
  for(auto sectorIter=nPoints.begin(); sectorIter!=nPoints.end(); ++sectorIter)
  {
    int sector=sectorIter->first;
    if(_fitButton->IsOn() && sumPEs[sector]>=10 && nPoints[sector]>1)
    {
      if(sumPEs[sector]*sumYY[sector]-sumY[sector]*sumY[sector]!=0)
      {
        float slope=(sumPEs[sector]*sumXY[sector]-sumX[sector]*sumY[sector])/(sumPEs[sector]*sumYY[sector]-sumY[sector]*sumY[sector]);
        float intercept=(sumX[sector]-slope*sumY[sector])/sumPEs[sector];
        float x1=intercept+slope*(_minY[sector]-(sector==0?40:20));
        float x2=intercept+slope*(_maxY[sector]+(sector==0?40:20));
        for(int side=0; side<2; ++side)
        {
          if(side==1)
          {
            x1=_centerX+(_centerX-x1);
            x2=_centerX+(_centerX-x2);
          }
          _pads[side]->cd();
          _fit[side][sector] = new TLine(x1,_minY[sector]-(sector==0?40:20),x2,_maxY[sector]+(sector==0?40:20));
          _fit[side][sector]->SetLineColor(sector==0?3:6);
          _fit[side][sector]->SetLineWidth(3);
          _fit[side][sector]->Draw();
        }
      }
    }
  }

  for(int side=1; side>=0; --side)
  {
    _pads[side]->cd();
    _pads[side]->Modified();
    _pads[side]->Update();
    _pads[side]->GetCanvas()->cd();
    _pads[side]->GetCanvas()->Modified();
    _pads[side]->GetCanvas()->Update();
  }
  _pads[0]->SetName("padA");
  _pads[1]->SetName("padB");
  _pads[0]->Modified();
  _pads[1]->Modified();
  _pads[0]->GetCanvas()->Update();
  _pads[1]->GetCanvas()->Update();

  return true;
}

void EventWindow::DoExit()
{
  gApplication->Terminate(0);
}

void EventWindow::DrawFront(int side)
{
  _pads[side]=new TPad(side==0?"padA":"padB","",0,0.4,1,1);
  _pads[side]->Draw();

  _padsAdc[side]=new TPad(side==0?"padAdcA":"padAdcB","",0,0,0.7,0.4);
  _padsAdc[side]->Draw();

  _padsReco1[side]=new TPad(side==0?"padReco1A":"padReco1B","",0.7,0,1,0.2);
  _padsReco1[side]->Draw();

  _padsReco2[side]=new TPad(side==0?"padReco2A":"padReco2B","",0.7,0.2,1,0.4);
  _padsReco2[side]->Draw();

  _pads[side]->cd();

  TGraph *g = new TGraph();
  g->SetPoint(0,_minX[-1]-20,_minY[-1]-20);
  g->SetPoint(1,_maxX[-1]+20,_maxY[-1]+20);
  g->GetXaxis()->SetTitle("x [mm]");
  g->GetYaxis()->SetTitle("y [mm]");
  g->Draw(side==0?"A P":"A P RX");

  double minXpad=g->GetXaxis()->GetXmin();
  double maxXpad=g->GetXaxis()->GetXmax();
  double minYpad=g->GetYaxis()->GetXmin();
  double maxYpad=g->GetYaxis()->GetXmax();
        _centerX=0.5*(minXpad+maxXpad);
  double centerY=0.5*(minYpad+maxYpad);
  double diffX=maxXpad-minXpad;
  double diffY=maxYpad-minYpad;

  _textPEs[side] = new TText(_centerX,centerY-diffY*0.58,"");
  _textPEs[side]->SetTextAlign(22);
  _textPEs[side]->SetTextColor(4);
  _textPEs[side]->SetTextSize(0.05);
  _textPEs[side]->Draw();

  _textTotalPEs[side] = new TText(_centerX,centerY+diffY*0.56,side==0?"Side A - Total number of PEs: 0":"Side B - Total number of PEs: 0");
  _textTotalPEs[side]->SetTextAlign(22);
  _textTotalPEs[side]->SetTextColor(1);
  _textTotalPEs[side]->SetTextSize(0.05);
  _textTotalPEs[side]->Draw();

  for(auto channelIter=_channelMap.begin(); channelIter!=_channelMap.end(); ++channelIter)
  {
    if(channelIter->second._side!=side) continue; //drawing the other side at the moment
    int sector = channelIter->second._sector;
    int ignoreForFit = channelIter->second._ignoreForFit;
    double x = channelIter->second._x;
    double y = channelIter->second._y;

    double xDisplay=(side==1?_centerX+(_centerX-x):x);
    double yDisplay=y;

    int feb=channelIter->first.first;
    int channel1=channelIter->first.second.first;
    int channel2=channelIter->first.second.second;

    CrvCounter *box = new CrvCounter(side, sector, ignoreForFit, feb, channel1, channel2, _textPEs[side],
                                     _padsAdc[side], _padsReco1[side], _padsReco2[side], x, y, xDisplay, yDisplay);
    box->Draw();
    _boxes[side].push_back(box);

    for(int side=0; side<2; ++side) _fit[side][sector] = NULL;
  }
}

EventWindow::EventWindow(const TGWindow *p, UInt_t w, UInt_t h, const std::string &runFileName, const string &channelMapFileName) :
                         TGMainFrame(p, w, h), _runFileName(runFileName), _channelMapFileName(channelMapFileName)
{
   size_t dotPos = _runFileName.rfind('.');
   if(dotPos==std::string::npos) {std::cerr<<"The run file name has an unknown pattern."<<std::endl; return;}
   std::string runNumberTmp = _runFileName.substr(0,dotPos);
   dotPos = runNumberTmp.rfind('.');
   if(dotPos+1>=runNumberTmp.length()) {std::cerr<<"The run file name has an unknown pattern."<<std::endl; return;}
   _runNumber = runNumberTmp.substr(dotPos+1);

   _file = new TFile(_runFileName.c_str());
   _tree = (TTree*)_file->Get("run");
   if(_tree==NULL) {std::cerr<<"Could not find event tree."<<std::endl; return;}
   _nEntries = _tree->GetEntries();
   if(_nEntries<1)  {std::cerr<<"No events found!"<<std::endl; return; }

   //to read the numberOfFebs, channelsPerFeb, and numberOfSamples
   TTree *treeSpills = (TTree*)_file->Get("spills");
   treeSpills->SetBranchAddress("spill_number_of_febs", &_numberOfFebs);
   treeSpills->SetBranchAddress("spill_channels_per_feb", &_channelsPerFeb);
   treeSpills->SetBranchAddress("spill_number_of_samples", &_numberOfSamples);
   treeSpills->GetEntry(0);

   ReadChannelMap();

   for(int i=0; i<100; i++)
   {
     float r,g,b;
     new TColor(i+2000, (100.0-i)/100.0, (100.0-i)/100.0, 1.0);
   }

   TGHorizontalFrame *buttonFrame = new TGHorizontalFrame(this, 1600, 200);
   AddFrame(buttonFrame);

   //const TGFont *font1= gClient->GetFont("-adobe-helvetica-medium-r-*-*-25-*-*-*-*-*-iso8859-1");
   const TGFont *font1= gClient->GetFont("-adobe-helvetica-medium-r-*-*-20-*-*-*-*-*-iso8859-1");
   FontStruct_t ft1 = font1->GetFontStruct();


   TGLabel      *gotoTextSpill = new TGLabel(buttonFrame, "Go to spill");
   _gotoSpill = new TGNumberEntryField(buttonFrame,-1,0,TGNumberFormat::kNESInteger);
   TGLabel      *gotoTextEvent = new TGLabel(buttonFrame, "event");
   _gotoEvent = new TGNumberEntryField(buttonFrame,-1,0,TGNumberFormat::kNESInteger);
   TGTextButton *gotoButton = new TGTextButton(buttonFrame, "&Go");
   TGTextButton *prevButton = new TGTextButton(buttonFrame, "&Prev Event");
   TGTextButton *nextButton = new TGTextButton(buttonFrame, "&Next Event");
   _fitButton = new TGCheckButton(buttonFrame, "fit tracks");
   TGTextButton *saveButton = new TGTextButton(buttonFrame, "&Save");
   TGTextButton *exitButton = new TGTextButton(buttonFrame, "&Exit");

   gotoTextSpill->SetTextFont(ft1);
   gotoTextEvent->SetTextFont(ft1);
   _gotoSpill->SetFont(ft1);
   _gotoEvent->SetFont(ft1);
   gotoButton->SetFont(ft1);
   prevButton->SetFont(ft1);
   nextButton->SetFont(ft1);
   _fitButton->SetFont(ft1);
   saveButton->SetFont(ft1);
   exitButton->SetFont(ft1);

   gotoTextSpill->SetMargins(5,5,5,5);
   gotoTextEvent->SetMargins(5,5,5,5);
   gotoButton->SetMargins(5,5,5,5);
   prevButton->SetMargins(5,5,5,5);
   nextButton->SetMargins(5,5,5,5);
   saveButton->SetMargins(5,5,5,5);
   exitButton->SetMargins(5,5,5,5);

   gotoButton->Connect("Clicked()", "EventWindow", this, "DoGoToEvent()");
   prevButton->Connect("Clicked()", "EventWindow", this, "DoPrevEvent()");
   nextButton->Connect("Clicked()", "EventWindow", this, "DoNextEvent()");
   _fitButton->Connect("Clicked()", "EventWindow", this, "DoEvent()");  //redraws the same event with/without the fit
   saveButton->Connect("Clicked()", "EventWindow", this, "DoSave()");
   exitButton->Connect("Pressed()", "EventWindow", this, "DoExit()");

   buttonFrame->AddFrame(gotoTextSpill, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 5, 5));
   buttonFrame->AddFrame(_gotoSpill, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 5, 5));
   buttonFrame->AddFrame(gotoTextEvent, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 5, 5));
   buttonFrame->AddFrame(_gotoEvent, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 5, 5));
   buttonFrame->AddFrame(gotoButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 5, 5));
   buttonFrame->AddFrame(prevButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 5, 5));
   buttonFrame->AddFrame(nextButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 5, 5));
   buttonFrame->AddFrame(_fitButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 5, 5));
   buttonFrame->AddFrame(saveButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 5, 5));
   buttonFrame->AddFrame(exitButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 5, 5, 5, 5));

//   const TGFont *font2 = gClient->GetFont("-adobe-courier-bold-r-normal--34-240-100-100-m-200-iso8859-1");
   const TGFont *font2 = gClient->GetFont("-adobe-courier-bold-r-normal--28-240-100-100-m-200-iso8859-1");
   FontStruct_t ft2 = font2->GetFontStruct();
   _textRunNumber = new TGLabel(buttonFrame, Form("Run %s",_runNumber.c_str()));
   _textRunNumber->SetTextColor(gROOT->GetColor(6));
   _textRunNumber->SetTextFont(ft2);
   _textSpillNumber = new TGLabel(buttonFrame, "Spill -----");
   _textSpillNumber->SetTextColor(gROOT->GetColor(8));
   _textSpillNumber->SetTextFont(ft2);
   _textEventNumber = new TGLabel(buttonFrame, "Event -----");
   _textEventNumber->SetTextColor(gROOT->GetColor(8));
   _textEventNumber->SetTextFont(ft2);
   buttonFrame->AddFrame(_textRunNumber, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 10, 10, 5, 5));
   buttonFrame->AddFrame(_textSpillNumber, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 10, 10, 5, 5));
   buttonFrame->AddFrame(_textEventNumber, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 10, 10, 5, 5));

   TGHorizontalFrame *contentFrame = new TGHorizontalFrame(this, 1600, 800);
   AddFrame(contentFrame);

   TRootEmbeddedCanvas *canvasFrontA = new TRootEmbeddedCanvas("CanvasFrontA",contentFrame,800,800);
   TRootEmbeddedCanvas *canvasFrontB = new TRootEmbeddedCanvas("CanvasFrontB",contentFrame,800,800);
   contentFrame->AddFrame(canvasFrontA);
   contentFrame->AddFrame(canvasFrontB);

   canvasFrontA->GetCanvas()->cd();
   DrawFront(0);

   canvasFrontB->GetCanvas()->cd();
   DrawFront(1);

   _currentEntry = 0;
   _currentSpill = 0;
   _currentEvent = 0;
   DoEvent();
   _pads[0]->Modified();
   _pads[1]->Modified();
   _pads[0]->GetCanvas()->Update();
   _pads[1]->GetCanvas()->Update();

   // Set a name to the main frame
   SetWindowName("CRV Teststand Event Display");
   MapSubwindows();

   // Initialize the layout algorithm via Resize()
   Resize(GetDefaultSize());

   // Map main frame
   MapWindow();
}

EventWindow::~EventWindow()
{
   // Clean up main frame...
   Cleanup();
}

void DisplayEventsNew(const std::string &runFileName, const std::string &channelMapFileName)
{
   // Popup the GUI...
   new EventWindow(gClient->GetRoot(), 1600, 1600, runFileName, channelMapFileName);
}
