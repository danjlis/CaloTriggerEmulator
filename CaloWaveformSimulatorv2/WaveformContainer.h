#ifndef __WAVEFORMCONTAINER_H
#define __WAVEFORMCONTAINER_H

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <ostream>
#include <phool/PHObject.h>
#include <map>
#include <climits>
///
class WaveformContainer : public PHObject
{
 public:
  typedef std::map<unsigned int, std::vector<int>*> Map;
  typedef Map::const_iterator ConstIter;
  typedef Map::iterator Iter;
  typedef std::pair<Iter, Iter> Range;
  typedef std::pair<ConstIter, ConstIter> ConstRange;

  enum DETECTOR
  {
    ALL = 0,
    EMCAL = 1,
    HCALOUT = 2,
    HCALIN = 3,
    SEPD = 4,
    MBD = 5,
    DETECTOR_INVALID = 9999
  };

  WaveformContainer();
  virtual ~WaveformContainer();
  
  /// Clear Event from memory
  virtual void Reset() override;
  virtual void identify(std::ostream& os = std::cout);
  virtual int isValid();

  virtual std::vector<int>* get_waveform_at_eta_phi(int /* eta */, int /* phi */) { return nullptr; }
  virtual std::vector<int>* get_waveform_at_channel(int /* index */ ) { return nullptr; }

  virtual void set_waveform_at_eta_phi(int /* eta */ , int /* phi */, std::vector<int>*);
  virtual void set_waveform_at_channel(int /* index */ , std::vector<int>*);

  virtual size_t size() {return 0;}

 private:
  
  DETECTOR _det;
  Map _waveforms;

 private: // so the ClassDef does not show up with doc++
  ClassDefOverride(WaveformContainer,1);
};

#endif
