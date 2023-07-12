#ifndef __LL1OUTV2_H
#define __LL1OUTV2_H

#include <string>
#include <ostream>
#include <phool/PHObject.h>
#include "LL1Outv1.h"
#include "LL1Out.h"
#include "TriggerDefs.h"
#include "TriggerPrimitive.h"
#include "TriggerPrimitiveContainerv1.h"
#include <vector>


using namespace std;
///
class LL1Outv2 : public LL1Outv1
{
 public:
  
  ///
  LL1Outv2();

  LL1Outv2(std::string triggertype, std::string ll1type);
  ///
  virtual ~LL1Outv2() override;

  /// Clear Event from memory
  virtual void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
  */
  virtual void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  virtual int isValid() const override;

  // Get Trigger Type
  virtual std::string getLL1Type() const {return _ll1_type;}

  // Set Trigger Type  
  virtual void setLL1Type(std::string &ll1type) { _ll1_type = ll1type;}
  virtual void setTriggerType(std::string &triggertype) { _trigger_type = triggertype;}

  virtual TriggerDefs::TriggerKey getTriggerKey() const {return _trigger_key;}
  virtual void setTriggerKey(TriggerDefs::TriggerKey key) {_trigger_key = key;}

  virtual TriggerPrimitiveContainerv1* GetTriggerPrimitiveContainer() {return _trigger_primitives;}
  virtual std::vector<unsigned int>* GetTriggerBits() {return _trigger_bits;}
  
 protected:

  virtual void Init() override;
  
 private:
  
  std::string _ll1_type;
  std::string _trigger_type;

  TriggerDefs::TriggerKey _trigger_key = TriggerDefs::TRIGGERKEYMAX;

  int idx;

  vector<unsigned int> *_trigger_bits;
  vector<vector<unsigned int>> _trigger_words;
  TriggerPrimitiveContainerv1 *_trigger_primitives;
  unsigned int _thresholds[10];

 private: // so the ClassDef does not show up with doc++
  ClassDefOverride(LL1Outv2,1);
};

#endif
