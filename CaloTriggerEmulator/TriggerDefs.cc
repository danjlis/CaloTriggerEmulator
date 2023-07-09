#include "TriggerDefs.h"
#include <bitset>

uint32_t
TriggerDefs::getTriggerKey(const TriggerDefs::TriggerId triggerId)
{
  uint32_t tmp = triggerId << kBitShiftTriggerId;
  return tmp;
}

uint32_t
TriggerDefs::getTriggerPrimKey(const TriggerDefs::TriggerId triggerId, const TriggerDefs::DetectorId detectorId, const TriggerDefs::PrimitiveId primitiveId, const uint16_t primlocid)
{
  uint32_t tmp = triggerId;
  uint32_t key = triggerId << kBitShiftTriggerId;
  tmp = detectorId;
  key |= tmp << kBitShiftDetectorId;
  tmp = primitiveId;
  key |= tmp << kBitShiftPrimitiveId;
  tmp = primlocid;
  key |= tmp << kBitShiftPrimitiveLocId;
  return key;
}

uint32_t
TriggerDefs::getTriggerSumKey(const TriggerDefs::TriggerId triggerId, const TriggerDefs::DetectorId detectorId, const TriggerDefs::PrimitiveId primitiveId, const uint16_t primlocid, const uint16_t sumlocid)
{
  uint32_t tmp = triggerId;
  uint32_t key = triggerId << kBitShiftTriggerId;
  tmp = detectorId;
  key |= tmp << kBitShiftDetectorId;
  tmp = primitiveId;
  key |= tmp << kBitShiftPrimitiveId;
  tmp = primlocid;
  key |= tmp << kBitShiftPrimitiveLocId;
  tmp = sumlocid;
  key |= tmp << kBitShiftSumLocId;
  return key;
}

uint32_t
TriggerDefs::getTriggerId_from_TriggerKey(const TriggerDefs::TriggerKey triggerkey)
{
  uint32_t tmp = (triggerkey >> kBitShiftTriggerId);
  return tmp;
}

uint32_t
TriggerDefs::getTriggerId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t tmp =  (triggerprimkey >> kBitShiftTriggerId);
  return tmp;
}

uint32_t 
TriggerDefs::getTriggerId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t tmp =  (triggersumkey >> kBitShiftTriggerId);
  return tmp;
}

uint32_t
TriggerDefs::getDetectorId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t tmp =  (triggerprimkey >> kBitShiftDetectorId) & 0xf;
  return tmp;
}

uint32_t 
TriggerDefs::getDetectorId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t tmp =  (triggersumkey >> kBitShiftDetectorId) & 0xf;
  return tmp;
}

uint32_t 
TriggerDefs::getPrimitiveId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t tmp =  (triggerprimkey >> kBitShiftPrimitiveId) & 0xf;
  return tmp;
}

uint32_t 
TriggerDefs::getPrimitiveId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t tmp =  (triggersumkey >> kBitShiftPrimitiveId) & 0xf;
  return tmp;
}

uint16_t 
TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint16_t tmp =  (triggerprimkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  return tmp;
}

uint16_t 
TriggerDefs::getPrimitiveLocId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint16_t tmp =  (triggersumkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  return tmp;
}

uint16_t 
TriggerDefs::getPrimitivePhiId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t detId= getDetectorId_from_TriggerPrimKey(triggerprimkey);

  uint16_t tmp = (triggerprimkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
      return tmp/3;
      break;
    case TriggerDefs::DetectorId::hcaloutDId :
      return tmp/3;
      break;
    case TriggerDefs::DetectorId::emcalDId :
      return tmp/12;
      break;
    default :
      return UINT16_MAX;
      break;
    }
  return UINT16_MAX;
}

uint16_t 
TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t detId= getDetectorId_from_TriggerSumKey(triggersumkey);

  uint16_t tmp = (triggersumkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
      return tmp/3;
      break;
    case TriggerDefs::DetectorId::hcaloutDId :
      return tmp/3;
      break;
    case TriggerDefs::DetectorId::emcalDId :
      return tmp/12;
      break;
    default :
      return UINT16_MAX;
      break;
    }
  return UINT16_MAX;
}

uint16_t 
TriggerDefs::getPrimitiveEtaId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t detId= getDetectorId_from_TriggerPrimKey(triggerprimkey);

  uint16_t tmp = (triggerprimkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
      return tmp%3;
      break;
    case TriggerDefs::DetectorId::hcaloutDId :
      return tmp%3;
      break;
    case TriggerDefs::DetectorId::emcalDId :
      return tmp%12;
      break;
    default :
      return UINT16_MAX;
      break;
    }
  return UINT16_MAX;
}

uint16_t 
TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t detId= getDetectorId_from_TriggerSumKey(triggersumkey);

  uint16_t tmp = (triggersumkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
      return tmp%3;
      break;
    case TriggerDefs::DetectorId::hcaloutDId :
      return tmp%3;
      break;
    case TriggerDefs::DetectorId::emcalDId :
      return tmp%12;
      break;
    default :
      return UINT16_MAX;
      break;
    }
  return UINT16_MAX;
}

uint16_t 
TriggerDefs::getSumLocId(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint16_t tmp = (triggersumkey >> kBitShiftSumLocId) & 0xf;
  return tmp;
}
uint16_t 
TriggerDefs::getSumPhiId(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t detId= getDetectorId_from_TriggerSumKey(triggersumkey);

  uint16_t tmp = TriggerDefs::getSumLocId(triggersumkey);
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
      return tmp/4;
      break;
    case TriggerDefs::DetectorId::hcaloutDId :
      return tmp/4;
      break;
    case TriggerDefs::DetectorId::emcalDId :
      return tmp/4;
      break;
    default :
      return UINT16_MAX;
      break;
    }
  return UINT16_MAX;
}
uint16_t 
TriggerDefs::getSumEtaId(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t detId= getDetectorId_from_TriggerSumKey(triggersumkey);

  uint16_t tmp = TriggerDefs::getSumLocId(triggersumkey);
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
      return tmp%4;
      break;
    case TriggerDefs::DetectorId::hcaloutDId :
      return tmp%4;
      break;
    case TriggerDefs::DetectorId::emcalDId :
      return tmp%4;
      break;
    default :
      return UINT16_MAX;
      break;
    }
  return UINT16_MAX;
}

