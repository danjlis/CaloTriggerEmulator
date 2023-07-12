#ifndef _CALOTRIGGEREMULATOR_H__
#define _CALOTRIGGEREMULATOR_H__

#include <calowaveformsim/WaveformContainerv1.h>
#include "TriggerPrimitive.h"
#include "TriggerPrimitiveContainerv1.h"
#include "TriggerDefs.h"
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

  int process_trigger();

  //! Set TriggerType
  void setTriggerType(const std::string &name);
  bool CheckFiberMasks(TriggerDefs::TriggerPrimKey key);
  bool CheckChannelMasks(TriggerDefs::TriggerSumKey key);
  void Verbosity(const int verbosity) { _verbose = verbosity; }

  void identify();
 protected:
  std::string outfilename;
  std::string _ll1_nodename;
  std::string _waveform_nodename;

  Fun4AllHistoManager *hm;
  TFile *outfile;

  //!Trigger Type
  std::string _trigger;

  TriggerDefs::TriggerId _triggerid = TriggerDefs::TriggerId::noneTId;

  bool _do_hcalin;
  bool _do_hcalout;
  bool _do_cemc;
  bool _do_mbd;

  //! Waveform conatiner
  WaveformContainerv1 *_waveforms_hcalin;
  WaveformContainerv1 *_waveforms_hcalout;
  WaveformContainerv1 *_waveforms_cemc;
  WaveformContainerv1 *_waveforms_mbd;

  //! LL1 Out
  LL1Outv2 *_ll1out;
  TriggerPrimitiveContainerv1 *_primitives;

  LL1Outv2 *_ll1out_hcalin;
  TriggerPrimitiveContainerv1 *_primitives_hcalin;

  LL1Outv2 *_ll1out_hcalout;
  TriggerPrimitiveContainerv1 *_primitives_hcalout;

  LL1Outv2 *_ll1out_cemc;
  TriggerPrimitiveContainerv1 *_primitives_cemc;

  TriggerPrimitive *_primitive;
  std::vector<unsigned int> *_sum;
  std::vector<unsigned int> *_bits;
  TProfile *avg_primitive;
  TH2D *peak_primitive;
  TH2D *primitives;
  TH2D *trigger_fire_map;
  TH2D *h2_line_up;
  TH1D *h_nhit;
  TH1D *h_mbd_time;
  TH1D *h_mbd_charge;

  std::vector<TProfile*> v_avg_primitive_cemc;
  std::vector<TH2D*> v_peak_primitive_cemc;
  std::vector<TH2D*> v_primitives_cemc;
  std::vector<TH2D*> v_trigger_fire_map_cemc;

  std::vector<TProfile*> v_avg_primitive_hcalin;
  std::vector<TH2D*> v_peak_primitive_hcalin;
  std::vector<TH2D*> v_primitives_hcalin;
  std::vector<TH2D*> v_trigger_fire_map_hcalin;

  std::vector<TProfile*> v_avg_primitive_hcalout;
  std::vector<TH2D*> v_peak_primitive_hcalout;
  std::vector<TH2D*> v_primitives_hcalout;
  std::vector<TH2D*> v_trigger_fire_map_hcalout;

  std::vector<TH1D*> v_nhit;
  std::vector<TH2D*> v_line_up;
  std::map<std::string, TH1D*> v_mbd_charge;
  std::map<std::string, TH1D*> v_mbd_time;


  //! Lookup tables
  unsigned int m_l1_adc_table[1024];
  unsigned int m_l1_adc_table_time[1024];
  unsigned int m_l1_slewing_table[4096];

  //! Trigger primitives
  unsigned int m_trig_charge[8];
  unsigned int m_trig_nhit;
  unsigned int m_trig_time[4];


  //! Trigger primitives
  unsigned int m2_trig_charge[4][8];
  unsigned int m2_trig_nhit[4];
  unsigned int m2_trig_time[4][4];

  std::vector<std::vector<unsigned int>*> _sum_mbd;

  //! Trigger ouputs
  std::vector<std::vector<unsigned int>*> _word_mbd;
  unsigned int m_out_tsum[2];
  unsigned int m_out_tavg[2];
  unsigned int m_out_trem[2];
  unsigned int m_out_nhit[2];
  unsigned int m_out_vtx_sub;
  unsigned int m_out_vtx_add;

  unsigned int m_nhit1, m_nhit2, m_timediff1, m_timediff2, m_timediff3;


  std::map<int, std::vector<int>*> m_peak_sub_ped_cemc;
  std::map<int, std::vector<int>*> m_peak_sub_ped_mbd;
  std::map<int, std::vector<int>*> m_peak_sub_ped_hcalin;
  std::map<int, std::vector<int>*> m_peak_sub_ped_hcalout;
  //! Verbosity.
  int _verbose;
  int _nevent;
  int _npassed;
  int _n_sums;
  int _n_primitives;
  int _m_trig_sub_delay;
  unsigned int _m_threshold;
  int m_isdata;
  int m_nsamples = 31;
  int _idx;

  std::vector<unsigned int> _masks_fiber;
  std::vector<unsigned int> _masks_channel;
  std::map<TriggerDefs::DetectorId, int> _m_prim_map;
  std::map<TriggerDefs::TriggerId, std::vector<std::string>> _m_det_map;
};

#endif
