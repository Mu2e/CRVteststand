# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <string>
# include <vector>
# include <map>
# include <experimental/filesystem>

# include "TTree.h"
# include "TFile.h"

const int CHANNEL_PER_FPGA=16;
const double CRV_TDC_RATE = 159.324e6;  // Hz
const int BOARD_STATUS_REGISTERS=22;

struct Event
{
  bool                                            _badEvent;
//  std::map<int,int>                               _tdcSinceSpill;  //the map key is the FEB //OLD
  std::map<std::pair<int,int> ,long>              _tdcSinceSpill;  //TDC for every FEB/channel pair
  std::map<std::pair<int,int>, std::vector<int> > _adc;            //ADC samples for every FEB/channel pair
  std::map<std::pair<int,int>, float>             _temperature;    //CMB temperature

  Event() : _badEvent(false) {}
};

class EventTree
{
  public:
  EventTree(const std::string &runNumber, const std::string &inFileName, const std::string &outFileName,
            const int numberOfFebs, const int channelsPerFeb, const int numberOfSamples);
  ~EventTree();
  void Finish();
  void Clear();
  void ReadSpill();
  void FillSpill();
  bool AtEOF() {return _inFile.eof();}
  int  GetSpillNumber() {return _spillNumber;}

  private:
  std::string    _runNumber;
  std::string    _inFileName, _outFileName;
  std::ifstream  _inFile;
  TFile         *_file;
  TTree         *_tree, *_treeSpills; 

  bool   _binary;

  //temp event data
  std::map<int,Event> _spill;    //map key is the event number
  bool   _missingFebs;           //missing FEBs or incomplete stabs (while stab marker was present)
  bool   _spillStored;           //spill was stored in the event tree
  size_t _nEventsActual;         //number of events that were stored in the event tree

  //tree event data
  int    _numberOfFebs;
  int    _channelsPerFeb;
  int    _numberOfSamples;
  size_t _nEvents;  //from spill header
  int    _spillNumber;
  int    _eventNumber;
  int   *_boardStatus;  //need one board status for each FEB

  //time stamp stored only for the first spill of the first subrun
  struct tm  _timestamp;

//  int    *_tdcSinceSpill;  //OLD
  long   *_tdcSinceSpill;
  double *_timeSinceSpill;
  int    *_adc;
  float  *_temperature;

  EventTree();
  void PrepareTree();
  bool ParseData(std::vector<unsigned> &buffer, size_t n, bool expectEOL=true, bool hexdec=false);
  bool ReadSpillHeader();
  bool ReadEventHeaderStart();
  bool ReadEventHeader(int &eventNumber, long &tdcSinceSpill, size_t &dataSize, bool missingBytes);
  //bool ReadEvent(Event &theEvent, int feb, size_t dataSize, bool checkMissingBytes, bool &missingBytes);  //OLD
  bool ReadEvent(Event &theEvent, long tdcSinceSpill, const std::vector<float> &temperatures, int feb, size_t dataSize, bool checkMissingBytes, bool &missingBytes);
  bool ReadStab(std::vector<float> &temperatures, int boardStatus[BOARD_STATUS_REGISTERS]);
};

EventTree::EventTree(const std::string &runNumber, const std::string &inFileName, const std::string &outFileName,
                     const int numberOfFebs, const int channelsPerFeb, const int numberOfSamples) : 
                     _runNumber(runNumber), _inFileName(inFileName), _outFileName(outFileName), _binary(false), 
                     _numberOfFebs(numberOfFebs), _channelsPerFeb(channelsPerFeb), _numberOfSamples(numberOfSamples)
{
  _boardStatus    = new int[_numberOfFebs*BOARD_STATUS_REGISTERS];
//  _tdcSinceSpill  = new int[_numberOfFebs];   //OLD
//  _timeSinceSpill = new double[_numberOfFebs];  //OLD
  _tdcSinceSpill  = new long[_numberOfFebs*_channelsPerFeb];
  _timeSinceSpill = new double[_numberOfFebs*_channelsPerFeb];
  _adc            = new int[_numberOfFebs*_channelsPerFeb*_numberOfSamples];
  _temperature    = new float[_numberOfFebs*_channelsPerFeb];

  _file = new TFile(_outFileName.c_str(), "recreate");
  _tree = new TTree("run","run");
  _treeSpills = new TTree("spills","spills");
  
  _inFile.open(_inFileName.c_str());
  if(!_inFile.good())
  {
    std::cout<<"Couldn't open input file "<<_inFileName<<"."<<std::endl;
    exit(0);
  }

  Clear();
  PrepareTree();
}

void EventTree::Clear()
{
  _missingFebs=false;
  _spillStored=true;
  _nEventsActual=0;
  for(int i=0; i<_numberOfFebs*BOARD_STATUS_REGISTERS; i++)
  {

    _boardStatus[i]=-1;
  }
  for(int i=0; i<_numberOfFebs*_channelsPerFeb; i++)
  {
    _tdcSinceSpill[i]=0;
    _timeSinceSpill[i]=0;
    _temperature[i]=0;
  }
  for(int i=0; i<_numberOfFebs*_channelsPerFeb*_numberOfSamples; i++) _adc[i]=0;
}

EventTree::~EventTree()
{
  delete [] _boardStatus;
  delete [] _tdcSinceSpill;
  delete [] _timeSinceSpill;
  delete [] _adc;
  delete [] _temperature;
}

void EventTree::PrepareTree()
{
  _tree->Branch("runtree_spill_num", &_spillNumber, "runtree_spill_num/I");
  _tree->Branch("runtree_event_num", &_eventNumber, "runtree_event_num/I");
//  _tree->Branch("runtree_tdc_since_spill", _tdcSinceSpill, Form("runtree_tdc_since_spill[%i]/I",_numberOfFebs)); //OLD
//  _tree->Branch("runtree_time_since_spill", _timeSinceSpill, Form("runtree_time_since_spill[%i]/D",_numberOfFebs)); //OLD
  _tree->Branch("runtree_tdc_since_spill", _tdcSinceSpill, Form("runtree_tdc_since_spill[%i][%i]/L",_numberOfFebs,_channelsPerFeb));
  _tree->Branch("runtree_time_since_spill", _timeSinceSpill, Form("runtree_time_since_spill[%i][%i]/D",_numberOfFebs,_channelsPerFeb));
  _tree->Branch("runtree_adc", _adc, Form("runtree_adc[%i][%i][%i]/I",_numberOfFebs,_channelsPerFeb,_numberOfSamples));
  _tree->Branch("runtree_temperature", _temperature, Form("runtree_temperature[%i][%i]/F",_numberOfFebs,_channelsPerFeb));

  _treeSpills->Branch("spill_num", &_spillNumber, "spill_num/I");
  _treeSpills->Branch("spill_nevents", &_nEvents, "spill_nevents/I");
  _treeSpills->Branch("spill_neventsActual", &_nEventsActual, "spill_neventsActual/I");
  _treeSpills->Branch("spill_stored", &_spillStored, "spill_stored/O");
  _treeSpills->Branch("spill_number_of_febs", &_numberOfFebs, "spill_number_of_febs/I");
  _treeSpills->Branch("spill_channels_per_feb", &_channelsPerFeb, "spill_channels_per_feb/I");
  _treeSpills->Branch("spill_number_of_samples", &_numberOfSamples, "spill_number_of_samples/I");
  _treeSpills->Branch("spill_boardStatus", _boardStatus, Form("spill_boardStatus[%i][%i]/I",_numberOfFebs,BOARD_STATUS_REGISTERS));
  _treeSpills->Branch("spill_timestamp_sec", &_timestamp.tm_sec);
  _treeSpills->Branch("spill_timestamp_min", &_timestamp.tm_min);
  _treeSpills->Branch("spill_timestamp_hour", &_timestamp.tm_hour);
  _treeSpills->Branch("spill_timestamp_mday", &_timestamp.tm_mday);
  _treeSpills->Branch("spill_timestamp_mon", &_timestamp.tm_mon);
  _treeSpills->Branch("spill_timestamp_year", &_timestamp.tm_year);
  _treeSpills->Branch("spill_timestamp_wday", &_timestamp.tm_wday);
  _treeSpills->Branch("spill_timestamp_yday", &_timestamp.tm_yday);
  _treeSpills->Branch("spill_timestamp_isdst", &_timestamp.tm_isdst);
}

void EventTree::Finish()
{
  _tree->Write("", TObject::kOverwrite);
  _treeSpills->Write("", TObject::kOverwrite);
  _file->Close();
  _inFile.close();
}

bool EventTree::ParseData(std::vector<unsigned> &buffer, size_t n, bool expectEOL, bool hexdec)
{
  if(!_binary)
  {
    for(size_t i=0; i<n; i++)
    {
      if(hexdec) _inFile >> std::hex >> buffer[i];
      else _inFile >> buffer[i];
      if(_inFile.fail())
      {
        _inFile.clear();
        return false;
      }
    }
  }
  else
  {
    size_t additionalBytes=0;
    if(expectEOL) additionalBytes=3*n/16;    //3 additional bytes for line breaks every 16 bytes
    if(n<16 && expectEOL) additionalBytes=3; //3 additional bytes for line break at event header 
                                             //(which is only 14 bytes due to the additional "8 8" at the beginning)
    unsigned char binaryLine[n+additionalBytes];
    _inFile.read(reinterpret_cast<char*>(binaryLine),n+additionalBytes);
    if(_inFile.fail())
    {
      _inFile.clear();
      return false;
    }
    for(size_t i=0, iBuffer=0; i<n+additionalBytes && iBuffer<buffer.size(); i++, iBuffer++) 
    {
      buffer[iBuffer]=binaryLine[i];
      if(expectEOL && iBuffer%16==15) i+=3;
    }
  }
  return true;
}

bool EventTree::ReadSpillHeader()
{
  std::vector<unsigned> buffer(16);
  if(!ParseData(buffer, 16)) return false;

  _spillNumber=buffer.at(9)+256*buffer.at(8);
  _nEvents=buffer.at(7)+256*buffer.at(6)+256*256*buffer.at(5)+256*256*256*buffer.at(4);
//std::cout<<"Spill#, number of events "<<_spillNumber<<"  "<<_nEvents<<std::endl;

  return true;
}

bool EventTree::ReadEventHeaderStart()  //needs to be skipped, when the previous event had two missing bytes,
                                        //which was noted, because the "8 8 " mark showed at the end of the previous event
{
  std::vector<unsigned> buffer(2);
  ParseData(buffer, 2, false);
  if(buffer[0]==8 && buffer[1]==8) return true;  //the event needs to start with "8 8 "
//  if(buffer[0]==15 && buffer[1]==232) return true;  //the event needs to start with "15 232 " for the special case of the flash gate
  else return false;
}

bool EventTree::ReadEventHeader(int &eventNumber, long &tdcSinceSpill, size_t &dataSize, bool missingBytes)
{
  if(!missingBytes)
  {
    if(!ReadEventHeaderStart()) return false;
  }

  std::vector<unsigned> buffer(14);
  if(!ParseData(buffer, 14)) return false;

  eventNumber=buffer.at(7)+256*buffer.at(6)+256*256*buffer.at(5)+256*256*256*buffer.at(4);
  tdcSinceSpill=buffer.at(3)+256*buffer.at(2)+256*256*buffer.at(1)+256*256*256*buffer.at(0);
  dataSize=buffer.at(9)+256*buffer.at(8);
//std::cout<<"Event# "<<eventNumber<<"   tdcSinceSpill "<<tdcSinceSpill<<"   dataSize "<<dataSize<<std::endl;

  return true;
}

//bool EventTree::ReadEvent(Event &theEvent, int feb, size_t dataSize, bool checkMissingBytes, bool &missingBytes) //OLD
bool EventTree::ReadEvent(Event &theEvent, long tdcSinceSpill, const std::vector<float> &temperatures, int feb, size_t dataSize, bool checkMissingBytes, bool &missingBytes)
{
  //get the data block for this event and channel
  std::vector<unsigned> buffer(dataSize*2);
  if(!ParseData(buffer, dataSize*2))
  {
    theEvent._badEvent=true;
std::cout<<"Couldn't parse data"<<std::endl;
    return false;
  }

  int check128=buffer.at(0);
  if(check128!=128) theEvent._badEvent=true; //every channel in each event needs to start with 128
                                             //will check at the end whether it's a "two missing bytes" problem,
                                             //but need to continue reading the file
  if(check128!=128) std::cout<<"The 128 signature was missing!"<<std::endl;

  int channel=buffer.at(1);
  if(channel>=_channelsPerFeb || channel<0)
  {
    theEvent._badEvent=true;
std::cout<<"unknown channel "<<channel<<std::endl;
    return false;
  }
//std::cout<<"reading event for FEB, ch "<<feb<<"  "<<channel<<std::endl;
if(theEvent._badEvent) std::cout<<"bad Event at FEB "<<feb<<" ch "<<channel<<std::endl;

  theEvent._tdcSinceSpill[std::pair<int,int>(feb,channel)]=tdcSinceSpill;  //TDCs may differ between different FPGAs
  std::vector<int> &adcSamples = theEvent._adc[std::pair<int,int>(feb,channel)];  //get the ADC samples for this FEB/channel pair
  for(size_t i=1; i<dataSize; i++)
  {
    int adcSample=buffer.at(i*2+1)+256*buffer.at(i*2);
    if(adcSample>2048) adcSample-=4096;
    adcSamples.push_back(adcSample);
//std::cout<<"ADCsample "<<i<<"  "<<adcSample<<std::endl;
  } 
  theEvent._temperature[std::pair<int,int>(feb,channel)]=temperatures[channel];

  if(checkMissingBytes && buffer.at(dataSize*2-2)==8 && buffer.at(dataSize*2-1)==8)   //at last channel of an FPGA
  {
std::cout<<"missingBytes FEB/ch "<<feb<<"  "<<channel<<std::endl;
    theEvent._badEvent=true;  //this event is bad, but the next events can still be read
    missingBytes=true;        //this event is short by two bytes. 
                              //the last two bytes which were read (" 8 8 "), are actually the first two bytes
                              //of the event header of the next event.
                              //therefore, when the next event header is read, the first two bytes need to be skipped,
                              //i.e. read only 14 instead of 16 bytes.
  }

  return true;
}

bool EventTree::ReadStab(std::vector<float> &temperatures, int boardStatus[BOARD_STATUS_REGISTERS])
{
  std::string line;
  std::vector<unsigned> buffer(16);

  //line with a marker and two board block lines
  //new FEBs: marker shows up in a separate line
  //older FEBs: marker shows up in first board block line
  std::streampos startOfStab = _inFile.tellg();
  if(!getline(_inFile, line)) return false;
  if(line.size()!=3) //no separate marker: old FEB
  {
    _inFile.seekg(startOfStab);  //need to go back one line
    _inFile.seekg(3,std::ios_base::cur);  //skip the marker
  }

  //1st board block line
  bool binaryTmp=_binary;
  _binary=false;
  if(!ParseData(buffer, 16, true, true)) {_binary=binaryTmp; return false;}
  for(int i=0; i<16; ++i) boardStatus[i]=buffer[i];

  //2nd board block line
  if(!ParseData(buffer, 6, true, true)) {_binary=binaryTmp; return false;}
  _binary=binaryTmp;
  for(int i=0; i<6; ++i) boardStatus[i+16]=buffer[i];

  if(!getline(_inFile, line)) return false;  //seems to be needed to get to the end of the line

  //FPGA blocks
//  int nFPGA = _channelsPerFeb/16;  //FIXME: make sure that nChannels are always multiples of 16 or switch to nFPGAs for user parameter
  int nFPGA = 4;   //FIXME: the DAQ program always writes the data of all 4 FPGAs
  for(int iFPGA=0; iFPGA<nFPGA; ++iFPGA)
  {
    if(!getline(_inFile, line)) return false;  //not used

    bool binaryTmp=_binary;
    _binary=false;
    if(!ParseData(buffer, 16, true, true)) {_binary=binaryTmp; return false;}
    _binary=binaryTmp;
    for(int iCMB=0; iCMB<4; ++iCMB)
    {
      for(int i=0; i<4; ++i) temperatures[16*iFPGA+4*iCMB+i]=buffer[8+iCMB]*0.062;
    }

    if(!getline(_inFile, line)) return false;  //for EOL at 2nd line
    if(!getline(_inFile, line)) return false;  //not used
  }

  return true;
}

void EventTree::ReadSpill()
{
  _spill.clear();
  _timestamp = {0,0,0,0,0,0,0,0,0};

std::cout<<"start new spill"<<std::endl;
  int feb=-1;

  while(1)
  {
    std::streampos startOfSpill = _inFile.tellg();
    std::string line;
    if(!getline(_inFile, line)) return;

    if(line.find("START OF RUN")!=std::string::npos)
    {
      std::string timestampString=line.substr(19);
      if(strptime(timestampString.c_str(), "%m/%d/%Y %I:%M:%S %p", &_timestamp)!=NULL)
      {
        std::cout<<"Found time stamp "<<timestampString<<std::endl;
      }
    }
  
    if(line.find("--Begin of spill")>1) continue;  //not at the spill, yet (at 0: text file, at 1: binary file)
    if(line.at(0)==18) _binary=true;

    std::streampos potential2ndBeginOfSpill = _inFile.tellg();
    if(!getline(_inFile, line)) return;   //end of file
    if(line.find("--Begin of spill")<=1)        //Found a subsequent "Begin of spill" marker.
    {                                           //This can happen, if there was a problem with
      _inFile.seekg(potential2ndBeginOfSpill);  //the FEB from the previous "Begin of spill" marker.
      _missingFebs=true;
      std::cout<<"missing data for FEB. skipping this FEB."<<std::endl;
      continue;
    }

    //Check whether there is stab data in the file
    std::vector<float> CMBtemperatures;
//    CMBtemperatures.resize(_channelsPerFeb);
    CMBtemperatures.resize(64);  //FIXME: stabs always contain the temperatures of all 4 FPGAs (i.e. all 64 channels)
    int boardStatusTmp[BOARD_STATUS_REGISTERS];
    if(line.find("stab")<2)
    {
      if(!ReadStab(CMBtemperatures,boardStatusTmp))
      {
        std::cout<<"incomplete stabs. skipping this FEB."<<std::endl;
        _missingFebs=true;
        continue;
      }
      if(!getline(_inFile, line)) return;
    }

    //FEB numbers start at 0 in new files and start at 1 in older files
    std::string searchString="--** SOURCE = FEB";
    size_t searchStringPos=line.find(searchString);
    if(searchStringPos==std::string::npos) continue;

    int newfeb = std::stoi(line.substr(searchStringPos+searchString.size()));
    if(newfeb<=feb)
    {
      //seems to have found a new spill (and not just a new FEB of the same spill). go back in the file to the beginning of this spill
      _inFile.seekg(startOfSpill);
      return;
    }
    if(newfeb>=_numberOfFebs)
    {
      std::cout<<"Found FEB number "<<feb<<" which is higher than expected."<<std::endl;
      continue;
    }
    feb=newfeb;
    std::cout<<"Reading feb "<<feb<<std::endl;

    for(int iRegister=0; iRegister<BOARD_STATUS_REGISTERS; ++iRegister)
    {
      _boardStatus[feb*BOARD_STATUS_REGISTERS+iRegister]=boardStatusTmp[iRegister];
    }

    //in the spill for the current FEB. reading the events for all channels of this FEB now
    if(!ReadSpillHeader())
    {
      //couldn't read the spill header.
      std::cout<<"Couldn't read spill header. Moving to the next FEB."<<std::endl;
      continue;
    }

    bool missingBytes=false;
    while(1)
    {
      int    eventNumber;
      long   tdcSinceSpill;
      size_t dataSize;
      if(!ReadEventHeader(eventNumber,tdcSinceSpill,dataSize,missingBytes))
      {
//std::cout<<"Couldn't read event header. Usually at end of spill for current FEB. Also happens for corrupted events. Moving to the next FEB."<<std::endl;
        break;
        //couldn't read the event header. 
        //-probably end of spill for the current FEB.
        //-could also be a corrupted event (but not due to missing bytes)
        //leaving the loop and go to the next FEB.
      }

      missingBytes=false;
      Event &theEvent = _spill[eventNumber]; 
//      theEvent._tdcSinceSpill[feb]=tdcSinceSpill;  //OLD This is now done in ReadEvent

      for(int channelInFpga=0; channelInFpga<CHANNEL_PER_FPGA; channelInFpga++)
      {
//std::cout<<"EVENT  "<<feb<<"  "<<eventNumber<<"  "<<channelInFpga<<"  "<<dataSize<<"  "<<tdcSinceSpill<<"         "<<_spill.size()<<std::endl;
        bool checkMissingBytes=(channelInFpga==CHANNEL_PER_FPGA-1);
//        if(!ReadEvent(theEvent, feb, dataSize, checkMissingBytes, missingBytes))  //OLD
        if(!ReadEvent(theEvent, tdcSinceSpill, CMBtemperatures, feb, dataSize, checkMissingBytes, missingBytes))
        {
          std::cout<<"Event "<<eventNumber<<" at FEB "<<feb<<" cannot be read. ";
          std::cout<<"Skipping event. May not be able to recover this spill."<<std::endl;
          break; //It couldn't be parsed due to non-numerical data. 
                 //Won't be able to parse the next event header either, 
                 //since the ifstream >> operator won't get past the non-numerical data.
                 //Therefore leaving the while loop completely, and go to the next FEB
        }
        else
        {
          if(missingBytes)
          {
            std::cout<<"Event "<<eventNumber<<" at FEB "<<feb<<" is damaged due to missing two bytes. ";
            std::cout<<"Skipping event."<<std::endl;
            break;
          }
          else if(theEvent._badEvent && checkMissingBytes) //It checked for the "missing bytes situation" at the end of the 16 channels 
                                                           //(by searching for the " 8 8 " string), but there were no missing bytes
                                                           //(it would have gone into the missingBytes branch above).
                                                           //So the bad event flag must have been caused by something else.
                                                           //The spill is now considered unrecoverable.
                                                           //It will fail at the next event header, since the " 8 8 " will most likely not be there.
          {
            std::cout<<"Event "<<eventNumber<<" at FEB "<<feb<<" cannot be read. ";
            std::cout<<"Skipping event. May not be able to recover this spill."<<std::endl;
            break;
          }
        }
      }
    }//end while(1) / loop over all FPGAs of the FEB for this spill
  }
}

void EventTree::FillSpill()
{
  if(_spill.empty())
  {
    std::cout<<"Spill "<<_spillNumber<<" is empty."<<std::endl;
    return;
  }

  if(_missingFebs)
  {
    std::cout<<"Spill "<<_spillNumber<<" has a missing FEB or incomplete stabs."<<std::endl;
    std::cout<<"Events of this spill are not stored, but spill information gets added to spill tree."<<std::endl;
    _spillStored=false;
    _treeSpills->Fill();
    return;
  }

  if(_spill.size()!=_nEvents)
  {
    std::cout<<"Expected "<<_nEvents<<" events in spill "<<_spillNumber<<". Found "<<_spill.size()<<" events."<<std::endl;
  }

  std::map<int,Event>::const_iterator event;
  for(event=_spill.begin(); event!=_spill.end(); event++)
  {
    _eventNumber = event->first;

    const Event &theEvent = event->second;
    if(theEvent._badEvent) continue;   //don't store this event, even if it was bad in only one FEB

    for(int feb=0; feb<_numberOfFebs; feb++)
    {
/*  //OLD
      if(theEvent._tdcSinceSpill.find(feb)==theEvent._tdcSinceSpill.end())
      {
//        std::cout<<"Missing FEB "<<feb<<" at event number "<<_eventNumber<<". "<<std::endl;
        _timeSinceSpill[feb]=NAN;  //will be used as identifier for missing FEBs in an event
        continue;
      }

      _tdcSinceSpill[feb] = theEvent._tdcSinceSpill.at(feb);
      _timeSinceSpill[feb] = _tdcSinceSpill[feb] / CRV_TDC_RATE;
*/

      for(int channel=0; channel<_channelsPerFeb; channel++)
      {
        int indexTDC=feb*_channelsPerFeb+channel;
        std::map<std::pair<int,int>, long>::const_iterator iterTDC=theEvent._tdcSinceSpill.find(std::pair<int,int>(feb,channel));
        if(iterTDC!=theEvent._tdcSinceSpill.end()) 
        {
          _tdcSinceSpill[indexTDC]=iterTDC->second; //_tdcSinceSpill[feb][channel];
          _timeSinceSpill[indexTDC]=iterTDC->second / (CRV_TDC_RATE/1e9); //_timeSinceSpill[feb][channel];
        }
        else
        {
//          std::cout<<"Missing FEB "<<feb<<" at event number "<<_eventNumber<<". "<<std::endl;
          //NAN will be used as identifier for missing FEBs/channels
          _timeSinceSpill[indexTDC]=NAN; //_timeSinceSpill[feb][channel];
          continue;   //no need to fill ADCs below; they are set to zero in the constructor
        }

        std::map<std::pair<int,int>, std::vector<int> >::const_iterator iterADC=theEvent._adc.find(std::pair<int,int>(feb,channel));
        if(iterADC!=theEvent._adc.end())
        {
          int nSamples=static_cast<int>(iterADC->second.size());
          for(int sample=0; sample<_numberOfSamples && sample<nSamples; sample++)
          {
            int index=feb*_channelsPerFeb*_numberOfSamples+channel*_numberOfSamples+sample;
            _adc[index]=iterADC->second.at(sample); //_adc[feb][channel][sample];
          }
        }

        int indexTemperature=feb*_channelsPerFeb+channel;
        std::map<std::pair<int,int>, float>::const_iterator iterTemperature=theEvent._temperature.find(std::pair<int,int>(feb,channel));
        if(iterTemperature!=theEvent._temperature.end()) 
        {
          _temperature[indexTemperature]=iterTemperature->second; //_temperature[feb][channel];
        }
      }//channels
    }//febs

    _tree->Fill();
    ++_nEventsActual;
  }//events

  _treeSpills->Fill();

  std::cout<<"Found "<<_nEventsActual<<" good events out of "<<_nEvents<<" total events in spill "<<_spillNumber<<"."<<std::endl;
}

void makeFileNames(const std::string &runNumber, std::string &inFileName, std::string &outFileName)
{
  std::ifstream dirFile;
  dirFile.open("config.txt");
  if(!dirFile.is_open()) {std::cerr<<"Could not open config.txt."<<std::endl; exit(1);}

  std::string inDirName, outDirName;
  std::string dirType, dir;
  while(dirFile>>dirType>>dir)
  {
    if(dirType=="crvraw") inDirName=dir;
    if(dirType=="crvparsed") outDirName=dir;
  }
  dirFile.close();

  bool found=false;
  for(const auto& dirEntry : std::experimental::filesystem::directory_iterator(inDirName)) 
  {
    const std::string s = dirEntry.path().filename().string();
    const std::string s0 = "crv.raw.";
    const std::string s1 = ".run"+runNumber+".txt";
    if(s.compare(0,s0.length(),s0)!=0) continue;
    if(s.compare(s.length()-s1.length(),s1.length(),s1)!=0) continue;

    found=true;
    inFileName = dirEntry.path().string();
    outFileName = outDirName+"crv.parsed."+dirEntry.path().stem().string().substr(s0.length())+".root";
    break;
  }

  if(!found)
  {
//try Ray's version of file names
    for(const auto& dirEntry : std::experimental::filesystem::directory_iterator(inDirName)) 
    {
      const std::string s = dirEntry.path().filename().string();
      const std::string s0 = "raw.mu2e.";
      const std::string s1 = "."+runNumber+".dat";
      if(s.compare(0,s0.length(),s0)!=0) continue;
      if(s.compare(s.length()-s1.length(),s1.length(),s1)!=0) continue;

      found=true;
      inFileName = dirEntry.path().string();
      outFileName = outDirName+"ntd.mu2e."+dirEntry.path().stem().string().substr(s0.length())+".root";
      break;
    }
  }

  if(!found) {std::cerr<<"Could not open input file for run "<<runNumber<<"."<<std::endl; exit(1);}
}

void printHelp()
{
  std::cout<<"Use as"<<std::endl;
  std::cout<<"parseCrv -h                      Prints this help."<<std::endl;
  std::cout<<"parseCrv RUNNUMBER [options]     Parses a raw file."<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Options:"<<std::endl;
  std::cout<<"-f   number of FEBs"<<std::endl;
  std::cout<<"     default: 4"<<std::endl;
  std::cout<<"-c   number of channels per FEB"<<std::endl;
  std::cout<<"     default: 64"<<std::endl;
  std::cout<<"-s   number of samples"<<std::endl;
  std::cout<<"     default: 127"<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Note: There needs to be a file config.txt at the current location,"<<std::endl;
  std::cout<<"which lists the location of the raw files, parsed files, etc."<<std::endl;
  exit(0);
}

int main(int argc, char **argv)
{
  if(argc==1) printHelp();
  else if(strcmp(argv[1],"-h")==0) printHelp();

  std::string runNumber=argv[1];
  std::string inFileName;
  std::string outFileName;
  makeFileNames(runNumber, inFileName, outFileName);

  int nFebs=4;
  int nChannels=64;
  int nSamples=127;
  for(int i=2; i<argc-1; i++)
  {
    if(strcmp(argv[i],"-f")==0) nFebs=atoi(argv[i+1]);
    if(strcmp(argv[i],"-c")==0) nChannels=atoi(argv[i+1]);
    if(strcmp(argv[i],"-s")==0) nSamples=atoi(argv[i+1]);
  }

  EventTree eventTree(runNumber, inFileName, outFileName, nFebs, nChannels, nSamples);

  while(1)
  {
    eventTree.Clear(); //sets all variables going to the tree to 0
    eventTree.ReadSpill();
    eventTree.FillSpill();
    if(eventTree.AtEOF()) break;
//    if(eventTree.GetSpillNumber()>4000) break;
  }
  eventTree.Finish();

  return 0;
}
