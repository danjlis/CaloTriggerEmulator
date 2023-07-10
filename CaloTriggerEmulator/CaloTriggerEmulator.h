#ifndef _CALOTRIGGEREMULATOR_H__
#define _CALOTRIGGEREMULATOR_H__

#include <calowaveformsim/WaveformContainerv1.h>
#include "TriggerPrimitive.h"
#include "TriggerPrimitiveContainerv1.h"
#include <cstdint>

#include <fun4all/SubsysReco.h>
#include <TTree.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TH2.h>
#include "LL1Outv2.h"

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
  void setTriggerType(const std::string &name) { _trigger = name; }

  void Verbosity(const int verbosity) { _verbose = verbosity; }

  void identify();
 protected:
  std::string outfilename;
  std::string _ll1_nodename;
  std::string _waveform_nodename;

  Fun4AllHistoManager *hm;
  TFile *outfile;
  TTree *_tree;

  //!Trigger Type
  std::string _trigger;
  
  //! Waveform conatiner
  WaveformContainerv1 *_waveforms;

  //! LL1 Out
  LL1Outv2 *_ll1out;
  TriggerPrimitiveContainerv1 *_primitives;
  TriggerPrimitive *_primitive;
  std::vector<unsigned int> *_sum;

  TProfile *avg_primitive;
  TH2D *peak_primitive;
  TH2D *primitives;
  TEfficiency *trigger_fire_map;

  std::vector<TProfile*> v_avg_primitive;
  std::vector<TH2D*> v_peak_primitive;
  std::vector<TH2D*> v_primitives;
  std::vector<TEfficiency*> v_trigger_fire_map;
  
  //! Lookup tables
  unsigned int m_l1_adc_table[1024];

  std::map<int, std::vector<int>*> m_peak_sub_ped;
  
  //! Verbosity
  int _verbose;
  int _nevent;

  int _n_sums;
  int _n_primitives;
  int _m_trig_sub_delay;

  int m_isdata;
  int m_nsamples = 31;

  std::map<std::string, int> _n_prim_map;
};

#endif
