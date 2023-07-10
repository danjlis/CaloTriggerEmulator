
#include "LL1ReturnCodes.h"
#include "LL1Outv2.h"
#include "TriggerDefs.h"
#include "TriggerPrimitive.h"
#include "TriggerPrimitiveContainerv1.h"

#include <cmath>
#include <iostream>

ClassImp(LL1Outv2)

LL1Outv2::LL1Outv2()
{
  _ll1_type = "NONE";
  _trigger_type = "NONE";

  setTriggerKey(TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(_trigger_type)));


  Init();
}


LL1Outv2::LL1Outv2(std::string triggertype, std::string ll1type)
{

  setTriggerType(triggertype);
  setTriggerKey(TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(triggertype)));

  setLL1Type(ll1type);

  Init();
}

LL1Outv2::~LL1Outv2()
{
  Reset();
}

//______________________________________
void LL1Outv2::Init()
{
  _trigger_primitives = new TriggerPrimitiveContainerv1();
  _trigger_primitives->setTriggerType(TriggerDefs::GetTriggerId(_trigger_type));

}


//______________________________________
void LL1Outv2::Reset()
{
  _trigger_bits.clear();
  _trigger_words.clear();
  _trigger_primitives->Reset();
}

//______________________________________
void LL1Outv2::identify(std::ostream& out) const
{
  out << __FILE__ << __FUNCTION__ <<" LL1Out: Triggertype = " << _trigger_type << " LL1 type = " << _ll1_type << std::endl;

  out << __FILE__ << __FUNCTION__ <<" Size of trigger bits       : "<< _trigger_bits.size();
  out << __FILE__ << __FUNCTION__ <<" Size of trigger words      : "<< _trigger_words.size();
  out << __FILE__ << __FUNCTION__ <<" Trigger Primitive Container: "<< _trigger_primitives->size();
  _trigger_primitives->identify();

}

int LL1Outv2::isValid() const
{
  return _trigger_primitives->isValid();
}

