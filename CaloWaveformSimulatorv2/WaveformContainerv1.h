#ifndef __WAVEFORMCONTAINERV1_H
#define __WAVEFORMCONTAINERV1_H

#include <string>
#include <ostream>
#include <iostream>
#include <phool/PHObject.h>
#include "WaveformContainer.h"
#include <map>
#include <utility>

//#include <TClonesArray.h>
//#include <calobase/RawTowerDefs.h> 
///

class WaveformContainerv1 : public WaveformContainer
{
 public:

  typedef std::map<unsigned int, unsigned int> MapInfo;
  typedef MapInfo::const_iterator ConstIterInfo;
  typedef MapInfo::iterator IterInfo;
  typedef std::pair<IterInfo, IterInfo> RangeInfo;
  typedef std::pair<ConstIterInfo, ConstIterInfo> ConstRangeInfo;

  ///
  WaveformContainerv1();
  ///
  ~WaveformContainerv1();

  /// Clear Event from memory
  virtual void Reset() override;
  virtual void identify(std::ostream& os = std::cout) const override;
  virtual int isValid() const override;

  virtual void set_Detector(DETECTOR detector) { _det = detector; }  
  virtual void set_Detector(std::string &det_name);

  virtual std::vector<int>* get_waveform_at_channel(int /* index */ ) override;
  virtual std::vector<int>* get_waveform_at_key(int /* index */ ) ;

  virtual void set_waveform_at_channel(int /* index */, std::vector<int>* ) override;
  virtual void set_waveform_at_key(int /* index */, std::vector<int>* );

  ConstRange getWaveforms(void) const;
  Range getWaveforms(void);

  virtual void add_packet_event(unsigned int key, unsigned int evt) { _packet_event_number[key] = evt;}
  virtual void add_packet_clock(unsigned int key, unsigned int clk) { _packet_clock_number[key] = clk;}
  virtual void add_fem_event(unsigned int key, unsigned int evt) { _fem_event_number[key] = evt;}
  virtual void add_fem_clock(unsigned int key, unsigned int clk) { _fem_clock_number[key] = clk;}
  // To do encode key and decode key

  virtual RangeInfo get_packet_events()
  {
    return std::make_pair(_packet_event_number.begin(), _packet_event_number.end());
  }
  virtual RangeInfo get_packet_clocks()
  {
    return std::make_pair(_packet_clock_number.begin(), _packet_clock_number.end());
  }
  virtual RangeInfo get_fem_events()
  {
    return std::make_pair(_fem_event_number.begin(), _fem_event_number.end());
  }
  virtual RangeInfo get_fem_clocks()
  {
    return std::make_pair(_fem_clock_number.begin(), _fem_clock_number.end());
  }
  virtual ConstRangeInfo get_packet_events() const
  {
    return std::make_pair(_packet_event_number.begin(), _packet_event_number.end());
  }
  virtual ConstRangeInfo get_packet_clocks() const
  {
    return std::make_pair(_packet_clock_number.begin(), _packet_clock_number.end());
  }
  virtual ConstRangeInfo get_fem_events() const
  {
    return std::make_pair(_fem_event_number.begin(), _fem_event_number.end());
  }
  virtual ConstRangeInfo get_fem_clocks() const
  {
    return std::make_pair(_fem_clock_number.begin(), _fem_clock_number.end());
  }

  virtual size_t size() override { return _waveforms.size();}

 protected:
  
  Map _waveforms;
  DETECTOR _det;
  MapInfo _packet_event_number;
  MapInfo _packet_clock_number;
  MapInfo _fem_event_number;
  MapInfo _fem_clock_number;

 private: // so the ClassDef does not show up with doc++
  ClassDefOverride(WaveformContainerv1,1);
};

#endif
