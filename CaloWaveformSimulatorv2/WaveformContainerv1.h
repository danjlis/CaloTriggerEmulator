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

  virtual void set_waveform_at_channel(int /* index */, std::vector<int>* ) override;

  ConstRange getWaveforms(void) const;
  Range getWaveforms(void);

  virtual size_t size() override { return _waveforms.size();}

 protected:
  
  Map _waveforms;
  DETECTOR _det;

 private: // so the ClassDef does not show up with doc++
  ClassDefOverride(WaveformContainerv1,1);
};

#endif
