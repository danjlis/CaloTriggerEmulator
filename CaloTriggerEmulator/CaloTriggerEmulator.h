#ifndef _CALOTRIGGEREMULATOR_H__
#define _CALOTRIGGEREMULATOR_H__
#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>
#include <calowaveformsim/WaveformContainerv1.h>

#include <cstdint>

#include <fun4all/SubsysReco.h>
#include <TTree.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TH2.h>
#include "LL1Outv1.h"

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TProfile;
class TEfficiency;
class TH2D;

class CaloTriggerEmulator : public SubsysReco
{
 public:
  //! constructor
  CaloTriggerEmulator(const std::string &name = "CaloTriggerEmulator", const std::string &fname = "MyNtuple.root");

  //! destructor
  virtual ~CaloTriggerEmulator();

  //! full initialization
  int Init(PHCompositeNode *);

  //! full initialization
  int InitRun(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  //! reset variables
  void reset_vars();

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

  void SetThresholds();

  void Verbosity(const int verbosity) { _verbose = verbosity; }

 protected:
  std::string outfilename;
  Fun4AllHistoManager *hm;
  TFile *outfile;
  TTree *_tree;

  //!Trigger Type
  std::string _trigger;
  
  //! Waveform conatiner
  WaveformContainerv1 *_waveforms_hcal;

  //! LL1 Out
  LL1Outv1 *_ll1_hcal;
  std::vector<std::vector<unsigned int>> m_trigger_primitives;

  TProfile *avg_primitive[24];
  TH2D *peak_primitive[24];
  TH2D *primitives[24];
  TEfficiency *trigger_fire_map[24];
  TEfficiency *full_fire_map;
  
  //! Lookup tables
  unsigned int m_l1_adc_table[1024];

  //! Trigger primitives
  unsigned int m_trig_sums[24][16];

  unsigned int m_trigger_word;
  //! Thresholds
  unsigned int m_nhit1, m_nhit2, m_timediff1, m_timediff2, m_timediff3;
  //! Waveforms
  int m_waveforms_hcal[24*64][31];
  int m_peak_sub_ped[24*64][25];
  
  //! Verbosity
  int _verbose;
  int _nevent;

  int m_isdata;
  int m_nsamples = 16;
  int m_packet_low, m_packet_high;

  Event *_event;

};

#endif
