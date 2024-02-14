#ifndef _MBDTRIGGEREMULATOR_H__
#define _MBDTRIGGEREMULATOR_H__
#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>
#include <calowaveformsim/WaveformContainerv1.h>

#include <cstdint>

#include <fun4all/SubsysReco.h>
#include <TTree.h>
#include "LL1Outv1.h"

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;

class MBDTriggerEmulator : public SubsysReco
{
 public:
  //! constructor
  MBDTriggerEmulator(const std::string &name = "MBDTriggerEmulator");

  //! destructor
  virtual ~MBDTriggerEmulator();

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! full initialization
  int InitRun(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  //! end of run method
  int End(PHCompositeNode *) override;

  //! reset variables
  int ResetEvent(PHCompositeNode *) override;

  //! Get Nodes
  void GetNodes(PHCompositeNode *);

  //! Create Nodes
  void CreateNodes(PHCompositeNode *);

  //! MakePrimitives
  int process_waveforms();

  //! MakeTriggerOutput
  int process_primitives();

  //! Set TriggerType
  void TriggerType(const std::string &name) { _trigger = name; }

  void SetThresholds(int nhit1, int nhit2, int timediff1, int timediff2, int timediff3);

  void Verbosity(const int verbosity) { _verbose = verbosity; }

  void SetNSamples(const int ns) { m_nsamples = ns;}

 protected:

  //!Trigger Type
  std::string _trigger;
  
  //! Waveform conatiner
  WaveformContainerv1 *_waveforms_mbd;

  //! LL1 Out
  LL1Outv1 *_ll1_mbd;

  //! Lookup tables
  unsigned int m_l1_adc_table[1024];
  unsigned int m_l1_slewing_table[4096];

  //! Trigger primitives
  unsigned int m_trig_charge[4][8];
  unsigned int m_trig_nhit[4];
  unsigned int m_trig_time[4][4];

  //! Trigger ouputs
  unsigned int m_out_tsum[2];
  unsigned int m_out_nhit[2];
  unsigned int m_out_tavg[2];
  unsigned int m_out_trem[2];
  unsigned int m_out_vtx_sub;
  unsigned int m_out_vtx_add;

  unsigned int m_trigger_word;
  //! Thresholds
  unsigned int m_nhit1, m_nhit2, m_timediff1, m_timediff2, m_timediff3;
  //! Waveforms
  int m_waveforms_mbd[256][31];
  int m_peak_sub_ped[256][25];
  
  //! Verbosity
  int _verbose;
  int _nevent;

  int m_isdata;
  int m_nsamples;
  int m_trig_sub_delay;
  int m_trig_sample_phase;
  int m_packet_low, m_packet_high;

  Event *_event;

};

#endif
