# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <string>
# include <vector>
# include <map>
# include <set>
# include <experimental/filesystem>

# include "TMap.h"
# include "TObjString.h"
# include "TTree.h"
# include "TFile.h"

const int CHANNEL_PER_FEB=64;
const int CHANNEL_PER_FPGA=16;
const double CRV_TDC_RATE = 159.324e6;  // Hz
const int BOARD_STATUS_REGISTERS=22;
const int FPGA_BLOCK_REGISTERS=38;
const int FPGA_BLOCKS=4;

struct Event
{
  bool                                              _badEvent;
  std::map<std::pair<int,int> ,long>                _tdcSinceSpill;  //TDC for every FEB/channel pair
  std::map<std::pair<int,int>, std::vector<short> > _adc;            //ADC samples for every FEB/channel pair

  Event() : _badEvent(false) {}
};

class EventTree
{
  public:
  EventTree(const std::string &runNumber, const std::string &inFileName, const std::string &outFileName,
            const int channelsPerFeb, const int numberOfSamples, const std::string &configuration);
  ~EventTree();
  void Finish();
  void Clear();
  void ReadSpill();
  void FillSpill();
  bool AtEOF() {return _inFile.eof();}

  private:
  std::string    _runNumber;
  std::string    _inFileName, _outFileName;
  std::ifstream  _inFile;
  TFile         *_file;
  TTree         *_tree, *_treeSpills, *_treeMetaData;

  bool   _binary;

  //temp event data
  std::map<int,Event> _spill;    //map key is the event number
  bool   _missingFebs;           //missing FEBs or incomplete stabs (while stab marker was present)
  bool   _spillStored;           //spill was stored in the event tree
  size_t _nEventsActual;         //number of events that were stored in the event tree
  bool   _oldFileVersion;        //older files start with FEB 1

  //tree event data
  int    _numberOfFebs;
  int    _channelsPerFeb;
  int    _numberOfSamples;
  size_t _nEvents;  //from spill header
  int    _nFPGAs;   //from word count in spill header
  int    _run;
  int    _subrun;
  int    _spillIndex; //increased for every spill
  int    _spillNumber; //from FEB spill header
  int    _eventNumber;
  int   *_boardStatus;  //need one board status for each FEB
  int   *_FPGABlocks;   //need 4 FPGA blocks for each FEB regardless of how many FPGAs have data

  //older files: time stamp stored only for the first spill of the first subrun
  struct tm  _timestampStruct;
  Long64_t   _timestamp;   //time_t

  long   *_tdcSinceSpill;
  double *_timeSinceSpill;
  short  *_adc;
  float  *_biasVoltage; //for each channel
  float  *_temperature; //CMB temperature for each channel

  std::string _configuration;

  EventTree();
  void PrepareTree();
  bool ParseData(std::vector<unsigned> &buffer, size_t n, bool expectEOL=true, bool hexdec=false);
  bool ReadSpillHeader();
  bool ReadEventHeaderStart();
  bool ReadEventHeader(int &eventNumber, long &tdcSinceSpill, size_t &dataSize, bool missingBytes);
  bool ReadEvent(Event &theEvent, long tdcSinceSpill, int feb, size_t dataSize, bool checkMissingBytes, bool &missingBytes);
  bool ReadStab(float biasVoltageTmp[CHANNEL_PER_FEB], float temperatureTmp[CHANNEL_PER_FEB], int boardStatus[BOARD_STATUS_REGISTERS], int FPGABlocks[FPGA_BLOCKS*FPGA_BLOCK_REGISTERS]);
};

EventTree::EventTree(const std::string &runNumber, const std::string &inFileName, const std::string &outFileName,
                     const int channelsPerFeb, const int numberOfSamples, const std::string &configuration) :
                     _runNumber(runNumber), _inFileName(inFileName), _outFileName(outFileName), _binary(false),
                     _channelsPerFeb(channelsPerFeb), _numberOfSamples(numberOfSamples), _spillIndex(0), _configuration(configuration)
{
  _file = new TFile(_outFileName.c_str(), "recreate");
  _tree = new TTree("run","run");
  _treeSpills = new TTree("spills","spills");
  _treeMetaData = new TTree("metaData","metaData");

  _inFile.open(_inFileName.c_str());
  if(!_inFile.good())
  {
    std::cout<<"Couldn't open input file "<<_inFileName<<"."<<std::endl;
    exit(0);
  }

  //find number of FEBs in file
  std::cout<<"Determining the number of FEBs. This may take a while."<<std::endl;
  std::string line;
  std::string searchString="--** SOURCE = FEB";
  std::set<int> febs;
  while(getline(_inFile, line))
  {
    size_t searchStringPos=line.find(searchString);
    if(searchStringPos==std::string::npos) continue;
    int feb = std::stoi(line.substr(searchStringPos+searchString.size()));
    febs.insert(feb);
  }
  _inFile.clear();
  _inFile.seekg(0,std::ios_base::beg);

  //FEB numbers start at 0 in new files and start at 1 in older files
  _oldFileVersion=false;
  if(*febs.begin()==1) _oldFileVersion=true;
  _numberOfFebs=febs.size();
  std::cout<<"Found "<<_numberOfFebs<<" FEBs"<<std::endl;


  _boardStatus    = new int[_numberOfFebs*BOARD_STATUS_REGISTERS];
  _FPGABlocks     = new int[_numberOfFebs*FPGA_BLOCKS*FPGA_BLOCK_REGISTERS];
  _tdcSinceSpill  = new long[_numberOfFebs*_channelsPerFeb];
  _timeSinceSpill = new double[_numberOfFebs*_channelsPerFeb];
  _adc            = new short[_numberOfFebs*_channelsPerFeb*_numberOfSamples];
  _biasVoltage    = new float[_numberOfFebs*_channelsPerFeb];
  _temperature    = new float[_numberOfFebs*_channelsPerFeb];

  size_t underscorePos = runNumber.rfind('_');
  if(underscorePos==std::string::npos) _run=atoi(runNumber.c_str());
  else
  {
    _run=atoi(runNumber.substr(0,underscorePos).c_str());
    if(underscorePos+1<runNumber.size()) _subrun=atoi(runNumber.substr(underscorePos+1).c_str());
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
  for(int i=0; i<_numberOfFebs*FPGA_BLOCKS*FPGA_BLOCK_REGISTERS; i++)
  {
    _FPGABlocks[i]=-1;
  }
  for(int i=0; i<_numberOfFebs*_channelsPerFeb; i++)
  {
    _tdcSinceSpill[i]=0;
    _timeSinceSpill[i]=0;
    _biasVoltage[i]=0;
    _temperature[i]=-1000;
  }
  for(int i=0; i<_numberOfFebs*_channelsPerFeb*_numberOfSamples; i++) _adc[i]=0;
}

EventTree::~EventTree()
{
  delete [] _boardStatus;
  delete [] _FPGABlocks;
  delete [] _tdcSinceSpill;
  delete [] _timeSinceSpill;
  delete [] _adc;
  delete [] _biasVoltage;
  delete [] _temperature;
}

void EventTree::PrepareTree()
{
  _tree->Branch("runtree_run_num", &_run, "runtree_run_num/I");
  _tree->Branch("runtree_subrun_num", &_subrun, "runtree_subrun_num/I");
  _tree->Branch("runtree_spill_index", &_spillIndex, "runtree_spill_index/I");
  _tree->Branch("runtree_spill_num", &_spillNumber, "runtree_spill_num/I");
  _tree->Branch("runtree_event_num", &_eventNumber, "runtree_event_num/I");
  _tree->Branch("runtree_tdc_since_spill", _tdcSinceSpill, Form("runtree_tdc_since_spill[%i][%i]/L",_numberOfFebs,_channelsPerFeb));
  _tree->Branch("runtree_time_since_spill", _timeSinceSpill, Form("runtree_time_since_spill[%i][%i]/D",_numberOfFebs,_channelsPerFeb));
  _tree->Branch("runtree_adc", _adc, Form("runtree_adc[%i][%i][%i]/S",_numberOfFebs,_channelsPerFeb,_numberOfSamples));
  _tree->Branch("runtree_biasVoltage", _biasVoltage, Form("runtree_biasVoltage[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  _tree->Branch("runtree_temperature", _temperature, Form("runtree_temperature[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  _tree->Branch("runtree_boardStatus", _boardStatus, Form("runtree_boardStatus[%i][%i]/I",_numberOfFebs,BOARD_STATUS_REGISTERS));
  _tree->Branch("runtree_FPGABlocks", _FPGABlocks, Form("runtree_FPGABlocks[%i][%i][%i]/I",_numberOfFebs,FPGA_BLOCKS,FPGA_BLOCK_REGISTERS));
  _tree->Branch("runtree_spillTimestamp", &_timestamp,"runtree_spillTimestamp/L");

  _treeSpills->Branch("runNumber", &_run, "runNumber/I");
  _treeSpills->Branch("subrunNumber", &_subrun, "subrunNumber/I");
  _treeSpills->Branch("spill_index", &_spillIndex, "spill_index/I");
  _treeSpills->Branch("spill_num", &_spillNumber, "spill_num/I");
  _treeSpills->Branch("spill_nevents", &_nEvents, "spill_nevents/I");
  _treeSpills->Branch("spill_neventsActual", &_nEventsActual, "spill_neventsActual/I");
  _treeSpills->Branch("spill_stored", &_spillStored, "spill_stored/O");
  _treeSpills->Branch("spill_number_of_febs", &_numberOfFebs, "spill_number_of_febs/I");
  _treeSpills->Branch("spill_channels_per_feb", &_channelsPerFeb, "spill_channels_per_feb/I");
  _treeSpills->Branch("spill_number_of_samples", &_numberOfSamples, "spill_number_of_samples/I");
  _treeSpills->Branch("spill_biasVoltage", _biasVoltage, Form("spill_biasVoltage[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  _treeSpills->Branch("spill_temperature", _temperature, Form("spill_temperature[%i][%i]/F",_numberOfFebs,_channelsPerFeb));
  _treeSpills->Branch("spill_boardStatus", _boardStatus, Form("spill_boardStatus[%i][%i]/I",_numberOfFebs,BOARD_STATUS_REGISTERS));
  _treeSpills->Branch("spill_FPGABlocks", _FPGABlocks, Form("spill_FPGABlocks[%i][%i][%i]/I",_numberOfFebs,FPGA_BLOCKS,FPGA_BLOCK_REGISTERS));
  _treeSpills->Branch("spill_timestamp", &_timestamp,"spill_timestamp/L");
  _treeSpills->Branch("spill_timestamp_sec", &_timestampStruct.tm_sec);
  _treeSpills->Branch("spill_timestamp_min", &_timestampStruct.tm_min);
  _treeSpills->Branch("spill_timestamp_hour", &_timestampStruct.tm_hour);
  _treeSpills->Branch("spill_timestamp_mday", &_timestampStruct.tm_mday);
  _treeSpills->Branch("spill_timestamp_mon", &_timestampStruct.tm_mon);
  _treeSpills->Branch("spill_timestamp_year", &_timestampStruct.tm_year);
  _treeSpills->Branch("spill_timestamp_wday", &_timestampStruct.tm_wday);
  _treeSpills->Branch("spill_timestamp_yday", &_timestampStruct.tm_yday);
  _treeSpills->Branch("spill_timestamp_isdst", &_timestampStruct.tm_isdst);

  _treeMetaData->Branch("runNumber", &_run, "runNumber/I");
  _treeMetaData->Branch("subrunNumber", &_subrun, "subrunNumber/I");
  _treeMetaData->Branch("configuration", &_configuration);
  _treeMetaData->Fill();  //only one entry per file
}

void EventTree::Finish()
{
  _tree->Write("", TObject::kOverwrite);
  _treeSpills->Write("", TObject::kOverwrite);
  _treeMetaData->Write("", TObject::kOverwrite);
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
  int wordCount=buffer.at(3)+256*buffer.at(2)+256*256*buffer.at(1)+256*256*256*buffer.at(0);
  if(_nEvents>0)
  {
    if((wordCount-8)%((16*128+8)*_nEvents)==0) _nFPGAs=(wordCount-8)/((16*128+8)*_nEvents);
    else
    {
      std::cout<<"spill header word count "<<wordCount<<" doesn't match the number of events "<<_nEvents<<std::endl;
      _nFPGAs=4;
    }
  }
  else
  {
    std::cout<<"spill header indicates 0 events"<<std::endl;
    _nFPGAs=4;
  }
//std::cout<<"Spill# "<<_spillNumber<<"    number of events "<<_nEvents<<"    word count "<<wordCount<<"   nFPGAs "<<_nFPGAs<<std::endl;

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

bool EventTree::ReadEvent(Event &theEvent, long tdcSinceSpill, int feb, size_t dataSize, bool checkMissingBytes, bool &missingBytes)
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
  std::vector<short> &adcSamples = theEvent._adc[std::pair<int,int>(feb,channel)];  //get the ADC samples for this FEB/channel pair
  for(size_t i=1; i<dataSize; i++)
  {
    short adcSample=buffer.at(i*2+1)+256*buffer.at(i*2);   //buffer stays within unsigned char for ADC values
    if(adcSample>2048) adcSample-=4096;
    adcSamples.push_back(adcSample);
//std::cout<<"ADCsample "<<i<<"  "<<adcSample<<std::endl;
  }

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

bool EventTree::ReadStab(float biasVoltageTmp[CHANNEL_PER_FEB], float temperatureTmp[CHANNEL_PER_FEB], int boardStatus[BOARD_STATUS_REGISTERS], int FPGABlocks[FPGA_BLOCKS*FPGA_BLOCK_REGISTERS])
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
  for(int i=0; i<6; ++i) boardStatus[i+16]=buffer[i];
  _binary=binaryTmp;

  for(int iAFE=0; iAFE<8; ++iAFE)
  {
     for(int i=0; i<8; ++i) biasVoltageTmp[iAFE*8+i]=boardStatus[11+iAFE]*0.02;
  }

  if(!getline(_inFile, line)) return false;  //seems to be needed to get to the end of the line

  //FPGA blocks
  //the stab command always returns the data of all 4 FPGAs, even if some FPGAs are turned off
  for(int iFPGA=0; iFPGA<FPGA_BLOCKS; ++iFPGA)
  {
    binaryTmp=_binary;
    _binary=false;

    //1st line
    if(!ParseData(buffer, 16, true, true)) {_binary=binaryTmp; return false;}
    for(int i=0; i<16; ++i) FPGABlocks[iFPGA*FPGA_BLOCK_REGISTERS+i]=buffer[i];

    //2nd line
    if(!ParseData(buffer, 16, true, true)) {_binary=binaryTmp; return false;}
    for(int i=0; i<16; ++i) FPGABlocks[iFPGA*FPGA_BLOCK_REGISTERS+i+16]=buffer[i];
    for(int iCMB=0; iCMB<4; ++iCMB)
    {
      for(int i=0; i<4; ++i)
      {
        temperatureTmp[16*iFPGA+4*iCMB+i]=buffer[8+iCMB]/16.0;
        if((buffer[8+iCMB]>>11) & 1) temperatureTmp[16*iFPGA+4*iCMB+i]-=4096.0;
      }
    }

    //3rd line
    if(!ParseData(buffer, 6, true, true)) {_binary=binaryTmp; return false;}
    for(int i=0; i<6; ++i) FPGABlocks[iFPGA*FPGA_BLOCK_REGISTERS+i+32]=buffer[i];

    _binary=binaryTmp;

    if(!getline(_inFile, line)) return false;  //for EOL at 3rd line
  }

  return true;
}

void EventTree::ReadSpill()
{
  _spill.clear();
  _timestampStruct = {0,0,0,0,0,0,0,0,-1};
  _timestamp = 0;

  std::cout<<"start new spill"<<std::endl;
  ++_spillIndex;
  int feb=-1;

  while(1)
  {
    std::streampos startOfSpill = _inFile.tellg();
    std::string line;
    if(!getline(_inFile, line)) return;

    if(line.find("START OF RUN")!=std::string::npos)
    {
      std::string timestampString=line.substr(19);
      if(strptime(timestampString.c_str(), "%m/%d/%Y %I:%M:%S %p", &_timestampStruct)!=NULL)
      {
        _timestamp=mktime(&_timestampStruct);
        std::cout<<"Found time stamp "<<timestampString<<std::endl;
      }
    }

    if(line.find("timestamp")!=std::string::npos)
    {
      if(feb!=-1) //haven't started a new spill, yet
      {
        _inFile.seekg(startOfSpill);
        return;
      }
      std::string timestampString=line.substr(10);
      if(strptime(timestampString.c_str(), "%m/%d/%Y %I:%M:%S %p", &_timestampStruct)!=NULL)
      {
        _timestamp=mktime(&_timestampStruct);
        std::cout<<"Found time stamp "<<timestampString<<std::endl;
      }
    }

    if(line.find("--Begin of spill")>1) continue;  //not at the spill, yet (at 0: text file, at 1: binary file)
    if(line.at(0)==18) _binary=true;

    std::streampos potential2ndBeginOfSpill = _inFile.tellg();
    if(!getline(_inFile, line)) return;   //end of file
    if(line.find("--Begin of spill")<=1 || line.find("timestamp")!=std::string::npos)        //Found a subsequent "Begin of spill" marker.
    {                                           //This can happen, if there was a problem with
      _inFile.seekg(potential2ndBeginOfSpill);  //the FEB from the previous "Begin of spill" marker.
      _missingFebs=true;
      std::cout<<"missing data for FEB. skipping this FEB."<<std::endl;
      continue;
    }

    //Check whether there is stab data in the file
    //We don't know yet, what FEB number it is. That's why we have to use temporary variables
    //stabs always contain the temperatures and bias voltages of all 4 FPGAs (i.e. all 64 channels)
    float biasVoltageTmp[CHANNEL_PER_FEB];
    float temperatureTmp[CHANNEL_PER_FEB];
    int   boardStatusTmp[BOARD_STATUS_REGISTERS];
    int   FPGABlocksTmp[FPGA_BLOCKS*FPGA_BLOCK_REGISTERS];
    bool  foundStab=false;
    if(line.find("stab")<2)
    {
      if(!ReadStab(biasVoltageTmp,temperatureTmp,boardStatusTmp,FPGABlocksTmp))
      {
        std::cout<<"incomplete stabs. skipping this FEB."<<std::endl;
        _missingFebs=true;
        continue;
      }
      if(!getline(_inFile, line)) return;
      foundStab=true;
    }

    //FEB numbers start at 0 in new files and start at 1 in older files
    std::string searchString="--** SOURCE = FEB";
    size_t searchStringPos=line.find(searchString);
    if(searchStringPos==std::string::npos) continue;

    int newfeb = std::stoi(line.substr(searchStringPos+searchString.size()));
    if(_oldFileVersion) --newfeb;
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

    if(foundStab)
    {
      for(int iChannel=0; iChannel<_channelsPerFeb; ++iChannel)
      {
        int indexChannel=feb*_channelsPerFeb+iChannel;
        _temperature[indexChannel]=temperatureTmp[iChannel];
        _biasVoltage[indexChannel]=biasVoltageTmp[iChannel];
      }
      for(int iRegister=0; iRegister<BOARD_STATUS_REGISTERS; ++iRegister)
      {
        _boardStatus[feb*BOARD_STATUS_REGISTERS+iRegister]=boardStatusTmp[iRegister];
      }
      for(int iRegister=0; iRegister<FPGA_BLOCKS*FPGA_BLOCK_REGISTERS; ++iRegister)
      {
        _FPGABlocks[feb*FPGA_BLOCKS*FPGA_BLOCK_REGISTERS+iRegister]=FPGABlocksTmp[iRegister];
      }
    }

    //in the spill for the current FEB. reading the events for all channels of this FEB now
    if(!ReadSpillHeader())
    {
      //couldn't read the spill header.
      std::cout<<"Couldn't read spill header. Moving to the next FEB."<<std::endl;
      _missingFebs=true;
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

      for(int channelInFpga=0; channelInFpga<CHANNEL_PER_FPGA; channelInFpga++)
      {
//std::cout<<"EVENT  "<<feb<<"  "<<eventNumber<<"  "<<channelInFpga<<"  "<<dataSize<<"  "<<tdcSinceSpill<<"         "<<_spill.size()<<std::endl;
        bool checkMissingBytes=(channelInFpga==CHANNEL_PER_FPGA-1);
        if(!ReadEvent(theEvent, tdcSinceSpill, feb, dataSize, checkMissingBytes, missingBytes))
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

  //check whether all channels are present
  size_t eventsWithMissingChannels[_numberOfFebs]={};
  std::map<int,Event>::iterator event;
  for(event=_spill.begin(); event!=_spill.end(); event++)
  {
    Event &theEvent = event->second;
    for(int feb=0; feb<_numberOfFebs; feb++)
    {
      for(int channel=0; channel<16*_nFPGAs; channel++)
//      for(int channel=0; channel<_channelsPerFeb; channel++)
      {
        if(theEvent._tdcSinceSpill.find(std::pair<int,int>(feb,channel))==theEvent._tdcSinceSpill.end())
        {
          theEvent._badEvent=true;
          ++eventsWithMissingChannels[feb];
          break;  //to avoid double counting this event
        }
      }
    }
  }

  for(int feb=0; feb<_numberOfFebs; feb++)
  {
    if(eventsWithMissingChannels[feb]>=_spill.size())
    {
      std::cout<<"All events of spill "<<_spillNumber<<" have missing channels at FEB "<<feb<<"."<<std::endl;
      _boardStatus[feb*BOARD_STATUS_REGISTERS]=-1;  //this marks the FEB as bad for this spill
      _spillStored=false;
    }
  }
  if(!_spillStored)
  {
    _treeSpills->Fill();
    std::cout<<"This spill will not be stored, but spill information gets added to spill tree."<<std::endl;
    return;
  }

  for(int feb=0; feb<_numberOfFebs; feb++)
  {
    if(eventsWithMissingChannels[feb]>0)
    {
      std::cout<<eventsWithMissingChannels<<" events have missing channels at FEB "<<feb<<". These events are not stored."<<std::endl;
    }
  }

  for(event=_spill.begin(); event!=_spill.end(); event++)
  {
    _eventNumber = event->first;

    const Event &theEvent = event->second;
    if(theEvent._badEvent) continue;   //don't store this event, even if it was bad in only one FEB

    for(int feb=0; feb<_numberOfFebs; feb++)
    {
      for(int channel=0; channel<16*_nFPGAs; channel++)
//      for(int channel=0; channel<_channelsPerFeb; channel++)
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
          //this should actually not happen, because spills with missing channels are filtered out above
          std::cout<<"Channel "<<channel<<" at "<<feb<<" is missing for event number "<<_eventNumber<<". "<<std::endl;
          //NAN will be used as identifier for missing FEBs/channels
          _timeSinceSpill[indexTDC]=NAN; //_timeSinceSpill[feb][channel];
          continue;   //no need to fill ADCs below; they are set to zero in the constructor
        }

        std::map<std::pair<int,int>, std::vector<short> >::const_iterator iterADC=theEvent._adc.find(std::pair<int,int>(feb,channel));
        if(iterADC!=theEvent._adc.end())
        {
          int nSamples=static_cast<int>(iterADC->second.size());
          for(int sample=0; sample<_numberOfSamples && sample<nSamples; sample++)
          {
            int index=feb*_channelsPerFeb*_numberOfSamples+channel*_numberOfSamples+sample;
            _adc[index]=iterADC->second.at(sample); //_adc[feb][channel][sample];
          }
        }
      }//channels
    }//febs

    _tree->Fill();
    ++_nEventsActual;
  }//events

  _treeSpills->Fill();

  std::cout<<"Found "<<_nEventsActual<<" good events out of "<<_nEvents<<" total events in spill "<<_spillNumber<<"."<<std::endl;
}

void makeFileNames(const std::string &runNumber, std::string &inFileName, std::string &outFileName, std::string &configuration)
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

    //4th field of file name
    std::stringstream ss(s);
    std::string field;
    for(int i=0; i<4 && getline(ss,field,','); ++i) {if(i==3) configuration=field;}

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

      //4th field of file name
      std::stringstream ss(s);
      std::string field;
      for(int i=0; i<4 && getline(ss,field,'.'); ++i) {if(i==3) configuration=field;}

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
  std::string configuration;
  makeFileNames(runNumber, inFileName, outFileName, configuration);

  int nChannels=64;
  int nSamples=127;
  for(int i=2; i<argc-1; i++)
  {
    if(strcmp(argv[i],"-f")==0) std::cout<<"The -f argument is not needed anymore. The parser finds the number of FEBs"<<std::endl;
    if(strcmp(argv[i],"-c")==0) nChannels=atoi(argv[i+1]);
    if(strcmp(argv[i],"-s")==0) nSamples=atoi(argv[i+1]);
  }

  EventTree eventTree(runNumber, inFileName, outFileName, nChannels, nSamples, configuration);

  while(1)
  {
    eventTree.Clear(); //sets all variables going to the tree to 0
    eventTree.ReadSpill();
    eventTree.FillSpill();
    if(eventTree.AtEOF()) break;
  }
  eventTree.Finish();

  return 0;
}
