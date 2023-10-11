#ifndef MBDEMULATORTREEMAKER_H
#define MBDEMULATORTREEMAKER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

#include <calowaveformsim/WaveformContainerv1.h>
#include <calowaveformsim/LL1Outv2.h>
#include <calowaveformsim/TriggerPrimitive.h>
#include <calowaveformsim/TriggerPrimitiveContainerv1.h>
#include <calowaveformsim/TriggerDefs.h>
#include "TTree.h"
#include "TFile.h"


class PHCompositeNode;
class TriggerPrimitive;
class TriggerPrimitiveContainerv1;
class LL1Outv2;

class MBDEmulatorTreeMaker : public SubsysReco
{
 public:

  MBDEmulatorTreeMaker(const std::string &name = "MBDEmulatorTreeMaker", const std::string &outfilename = "trees_mbd_trig.root", const std::string &nodename = "LL1OUT_MBD_RAW");

  virtual ~MBDEmulatorTreeMaker();

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;
  
  void GetNodes (PHCompositeNode *topNode);

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void SetVerbosity(int verbo) ;
 private:
  void reset_tree_vars();

  int _verbosity;

  TFile *_f;
  TTree *_tree_mbd;
  std::string _foutname;
  std::string _nodename;
  int _i_event;

  //  LL1Outv2 *_ll1_mbd_raw;

  LL1Outv2 *_ll1_mbd;

  TriggerPrimitive *_trigger_primitive;
  TriggerPrimitiveContainerv1 *_trigger_primitives;

  std::vector<unsigned int> *_trigger_bits;

  unsigned int b_clock_mbd;
  unsigned int b_event_mbd;
  std::vector<unsigned int> b_trigger_bits_mbd;
  std::vector<unsigned int> b_trigger_words_mbd[8];
  std::vector<unsigned int> b_trigger_charge_mbd[4][8];
  std::vector<unsigned int> b_trigger_nhit_mbd[4];
  std::vector<unsigned int> b_trigger_time_mbd[4][4];

  
  /* unsigned int b_clock_ll1; */
  /* unsigned int b_event_ll1; */
  /* std::vector<unsigned int> b_trigger_bits_ll1; */
  /* std::vector<unsigned int> b_trigger_words_ll1[8]; */
  /* std::vector<unsigned int> b_trigger_charge_ll1[4][8]; */
  /* std::vector<unsigned int> b_trigger_nhit_ll1[4]; */
  /* std::vector<unsigned int> b_trigger_time_ll1[4][4]; */
  
  
};

#endif 
