
#include "TriggerPrimitive.h"
#include "TriggerDefs.h"
#include <cmath>
#include <iostream>

ClassImp(TriggerPrimitive)

TriggerPrimitive::TriggerPrimitive()
{

}

TriggerPrimitive::~TriggerPrimitive()
{

}

//______________________________________
void TriggerPrimitive::Reset()
{
  _sums.clear();
}

void TriggerPrimitive::add_sum(TriggerDefs::TriggerSumKey key, std::vector<unsigned int> *sum)
{
  _sums[key] = sum;
}
std::vector<unsigned int>*  TriggerPrimitive::get_sum_at_key(TriggerDefs::TriggerSumKey key)
{
  if (!_sums[key]) return nullptr;

  return _sums[key];
}


//______________________________________
void TriggerPrimitive::identify(std::ostream& out)
{
  out << "identify yourself: I am a TriggerPrimitive object" << std::endl;

}

int TriggerPrimitive::isValid()
{
  return (!_sums.empty());
}
