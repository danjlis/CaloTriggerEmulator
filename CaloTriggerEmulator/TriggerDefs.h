#ifndef __TRIGGERDEFS_H__
#define __TRIGGERDEFS_H__

#include <cstdint>
#include <iostream>

namespace TriggerDefs
{

  typedef uint32_t TriggerKey;
  typedef uint32_t TriggerPrimKey;
  typedef uint32_t TriggerSumKey;
  
  static TriggerKey TRIGGERKEYMAX __attribute__((unused)) = UINT32_MAX;
  static TriggerPrimKey TRIGGERPRIMKEYMAX __attribute__((unused)) = UINT32_MAX;
  static TriggerPrimKey TRIGGERSUMKEYMAX __attribute__((unused)) = UINT32_MAX;

  static const unsigned int kBitShiftTriggerId __attribute__((unused)) = 24;
  static const unsigned int kBitShiftDetectorId __attribute__((unused)) = 20;
  static const unsigned int kBitShiftPrimitiveId __attribute__((unused)) = 16;

  enum TriggerId
  {
    mbdTId = 0,
    jetTId = 1,
    pairTId = 2,
    cosmicTId = 3,
  };

  enum DetectorId
  {
    mbdDId = 0,
    hcalinDId = 1,
    hcaloutDId = 2,
    emcalDId = 3,
  };

  enum PrimitiveId
  {
    mbdPId = 0,
    calPId = 1,
  };

  static const unsigned int kBitShiftPrimitiveLocId __attribute__((unused)) = 4;
  static const unsigned int kBitShiftSumLocId __attribute__((unused)) = 0;

  uint32_t getTriggerKey(const TriggerDefs::TriggerId triggerId);
  uint32_t getTriggerPrimKey(const TriggerDefs::TriggerId triggerId, const TriggerDefs::DetectorId detectorId, const TriggerDefs::PrimitiveId primitiveId, const uint16_t primlocid);
  uint32_t getTriggerSumKey(const TriggerDefs::TriggerId triggerId, const TriggerDefs::DetectorId detectorId, const TriggerDefs::PrimitiveId primitiveId, const uint16_t primlocid, const uint16_t sumlocid);

  uint32_t getTriggerId_from_TriggerKey(const TriggerDefs::TriggerKey triggerkey);
  uint32_t getTriggerId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey);
  uint32_t getTriggerId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey);

  uint32_t getDetectorId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey);
  uint32_t getDetectorId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey);

  uint32_t getPrimitiveId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey);
  uint32_t getPrimitiveId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey);

  uint16_t getPrimitiveLocId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey);
  uint16_t getPrimitiveLocId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey);

  uint16_t getPrimitivePhiId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey);
  uint16_t getPrimitivePhiId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey);

  uint16_t getPrimitiveEtaId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey);
  uint16_t getPrimitiveEtaId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey);

  uint16_t getSumLocId(const TriggerDefs::TriggerSumKey triggersumkey);

  uint16_t getSumPhiId(const TriggerDefs::TriggerSumKey triggersumkey);

  uint16_t getSumEtaId(const TriggerDefs::TriggerSumKey triggersumkey);

};

#endif
