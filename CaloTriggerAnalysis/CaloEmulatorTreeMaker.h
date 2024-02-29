#ifndef CALOEMULATORTREEMAKER_H
#define CALOEMULATORTREEMAKER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

#include <calowaveformsim/LL1Outv2.h>
#include <calowaveformsim/TriggerPrimitive.h>
#include <calowaveformsim/TriggerPrimitiveContainerv1.h>
#include <calowaveformsim/TriggerDefs.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>

#include "TTree.h"
#include "TFile.h"


class PHCompositeNode;
class TriggerPrimitive;
class TriggerPrimitiveContainerv1;
class LL1Outv2;
class TowerInfoContainer;
class TowerInfo;
class PHG4TruthInfoContainer;
class CaloEmulatorTreeMaker : public SubsysReco
{
 public:

  CaloEmulatorTreeMaker(const std::string &name = "CaloEmulatorTreeMaker", const std::string &outfilename = "trees_calo_trig.root");

  virtual ~CaloEmulatorTreeMaker();

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;
  
  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  void SetTriggerID(TriggerDefs::TriggerId tid) { _triggerid = tid; }
  void SetTriggerType(const std::string trigger);

 private:

  TFile *_f;
  TTree *_tree;
  std::string _foutname;
  TriggerDefs::TriggerId _triggerid = TriggerDefs::TriggerId::noneTId;
  
  LL1Outv2 *_ll1out{nullptr};

  TriggerPrimitive *_trigger_primitive{nullptr};
  TriggerPrimitiveContainerv1 *_trigger_primitive_container{nullptr};

  TowerInfo *_tower;
  TowerInfoContainer* _towers;

  int b_pass;
  std::vector<int> b_ll1_8x8_nonoverlappingsum = {};
  std::vector<int> b_ll1_8x8_nonoverlappingsum_eta = {};
  std::vector<int> b_ll1_8x8_nonoverlappingsum_phi = {};
  std::vector<int> b_ll1_8x8_nonoverlappingsum_sample = {};

  std::vector<int> b_ll1_4x4_overlappingsum = {};
  std::vector<int> b_ll1_4x4_overlappingsum_eta = {};
  std::vector<int> b_ll1_4x4_overlappingsum_phi = {};
  std::vector<int> b_ll1_4x4_overlappingsum_sample = {};

  std::vector<int> b_ll1_2x2_emcal = {};
  std::vector<int> b_ll1_2x2_emcal_eta = {};
  std::vector<int> b_ll1_2x2_emcal_phi = {};
  std::vector<int> b_ll1_2x2_emcal_sample = {};

  std::vector<float> b_hcalin_energy_sim = {};
  std::vector<float> b_hcalin_etabin_sim = {};
  std::vector<float> b_hcalin_phibin_sim = {};

  std::vector<float> b_hcalout_energy_sim = {};
  std::vector<float> b_hcalout_etabin_sim = {};
  std::vector<float> b_hcalout_phibin_sim = {};

  std::vector<float> b_emcal_energy_sim = {};
  std::vector<float> b_emcal_etabin_sim = {};
  std::vector<float> b_emcal_phibin_sim = {};

  std::vector<float> b_hcalin_energy_calib = {};
  std::vector<float> b_hcalin_etabin_calib = {};
  std::vector<float> b_hcalin_phibin_calib = {};

  std::vector<float> b_hcalout_energy_calib = {};
  std::vector<float> b_hcalout_etabin_calib = {};
  std::vector<float> b_hcalout_phibin_calib = {};

  std::vector<float> b_emcal_energy_calib = {};
  std::vector<float> b_emcal_etabin_calib = {};
  std::vector<float> b_emcal_phibin_calib = {};

  std::vector<float> b_hcalin_energy_wave = {};
  std::vector<float> b_hcalin_etabin_wave = {};
  std::vector<float> b_hcalin_phibin_wave = {};

  std::vector<float> b_hcalout_energy_wave = {};
  std::vector<float> b_hcalout_etabin_wave = {};
  std::vector<float> b_hcalout_phibin_wave = {};

  std::vector<float> b_emcal_energy_wave = {};
  std::vector<float> b_emcal_etabin_wave = {};
  std::vector<float> b_emcal_phibin_wave = {};

  std::vector<float> b_hcalin_energy_raw = {};
  std::vector<float> b_hcalin_etabin_raw = {};
  std::vector<float> b_hcalin_phibin_raw = {};

  std::vector<float> b_hcalout_energy_raw = {};
  std::vector<float> b_hcalout_etabin_raw = {};
  std::vector<float> b_hcalout_phibin_raw = {};

  std::vector<float> b_emcal_energy_raw = {};
  std::vector<float> b_emcal_etabin_raw = {};
  std::vector<float> b_emcal_phibin_raw = {};

  int b_truth_particle_n = 0;
  std::vector<int> b_truth_particle_pid = {};
  std::vector<float> b_truth_particle_pt = {};
  std::vector<float>  b_truth_particle_eta = {};
  std::vector<float>  b_truth_particle_phi = {};
 
};

#endif 
