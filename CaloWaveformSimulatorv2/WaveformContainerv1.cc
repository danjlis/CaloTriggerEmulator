
#include "WaveformContainerv1.h"
#include <phool/PHObject.h>

#include <TClonesArray.h>
#include <cmath>
#include <iostream>
#include <cassert>

WaveformContainerv1::WaveformContainerv1()
{
  _det = WaveformContainerv1::DETECTOR::ALL;

  Reset();

}

WaveformContainerv1::~WaveformContainerv1()
{
  Reset();
}

//______________________________________
void WaveformContainerv1::Reset()
{
  while (_waveforms.begin() != _waveforms.end())
    {
      delete _waveforms.begin()->second;
      _waveforms.erase(_waveforms.begin());
    } 

  _packet_event_number.clear();
  _packet_event_number.clear();
  _fem_clock_number.clear();
  _fem_clock_number.clear();
}

//______________________________________
void WaveformContainerv1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a WaveformContainer object for " << _det << std::endl;
  out << "Size : " << _waveforms.size() << std::endl;
  ConstRange range = getWaveforms();
  for (auto iter = range.first; iter != range.second ; ++iter)
    {
      out << "Wave "<<(*iter).first<<": ";
      for ( auto i = (*iter).second->begin(); i != (*iter).second->end(); ++i) out << (*i) << " ";
      out <<" "<<std::endl;
    }
}


void WaveformContainerv1::set_Detector(std::string &det_name)
{
  if (strcmp(det_name.c_str(),"ALL") == 0)
    {
      _det = WaveformContainer::DETECTOR::ALL;
    }
 if (strcmp(det_name.c_str(),"CEMC") == 0)
    {
      _det = WaveformContainer::DETECTOR::CEMC;
    }
  if (strcmp(det_name.c_str(),"HCALIN") == 0)
    {
      _det = WaveformContainer::DETECTOR::HCALIN;
    }
  if (strcmp(det_name.c_str(),"HCALOUT") == 0)
    {
      _det = WaveformContainer::DETECTOR::HCALOUT;
    }
  if (strcmp(det_name.c_str(),"SEPD") == 0)
    {
      _det = WaveformContainer::DETECTOR::SEPD;
    }
  if (strcmp(det_name.c_str(),"MBD") == 0)
    {
      _det = WaveformContainer::DETECTOR::MBD;
    }
  else
    {
      _det = WaveformContainer::DETECTOR::DETECTOR_INVALID;
    }
  return;
}


int WaveformContainerv1::isValid() const 
{
  return (!_waveforms.empty());
}

std::vector<int>* WaveformContainerv1::get_waveform_at_channel(int index)
{
  if (!_waveforms[index]) return nullptr;
  
  return _waveforms[index];
}

std::vector<int>* WaveformContainerv1::get_waveform_at_key(int key)
{

  if (!_waveforms[key]) return nullptr;
  
  return _waveforms[key];
}


WaveformContainerv1::ConstRange WaveformContainerv1::getWaveforms() const
{
  return make_pair(_waveforms.begin(), _waveforms.end());
}

WaveformContainerv1::Range WaveformContainerv1::getWaveforms()
{
  return make_pair(_waveforms.begin(), _waveforms.end());
}


void WaveformContainerv1::set_waveform_at_channel(int index, std::vector<int>* wave)
{
  _waveforms[index] = wave;
}

void WaveformContainerv1::set_waveform_at_key(int key, std::vector<int>* wave)
{
  int index = key; 
  _waveforms[index] = wave;
}
