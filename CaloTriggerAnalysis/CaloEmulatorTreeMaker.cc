#include "CaloEmulatorTreeMaker.h"
#include <vector>
#include <TMath.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
// G4Cells includes

#include <iostream>

#include <map>

//____________________________________________________________________________..
CaloEmulatorTreeMaker::CaloEmulatorTreeMaker(const std::string &name, const std::string &outfilename):
  SubsysReco(name)
  
{
  _foutname = outfilename;  
}

//____________________________________________________________________________..
CaloEmulatorTreeMaker::~CaloEmulatorTreeMaker()
{

}
//____________________________________________________________________________..
void CaloEmulatorTreeMaker::SetTriggerType(const std::string trigger)
{
  SetTriggerID(TriggerDefs::GetTriggerId(trigger));
}
//____________________________________________________________________________..
int CaloEmulatorTreeMaker::Init(PHCompositeNode *topNode)
{

  _f = new TFile( _foutname.c_str(), "RECREATE");

  std::cout << " making a file = " <<  _foutname.c_str() << " , _f = " << _f << std::endl;
  
  _tree = new TTree("ttree","a persevering date tree");
  if (_triggerid == TriggerDefs::TriggerId::jetTId)
    {  
      _tree->Branch("ll1_8x8_nonovelappingsum", &b_ll1_8x8_nonoverlappingsum);
      _tree->Branch("ll1_8x8_nonovelappingsum_eta", &b_ll1_8x8_nonoverlappingsum_eta);
      _tree->Branch("ll1_8x8_nonovelappingsum_phi", &b_ll1_8x8_nonoverlappingsum_phi);
      _tree->Branch("ll1_8x8_nonovelappingsum_sample", &b_ll1_8x8_nonoverlappingsum_sample);
    }
  else if (_triggerid == TriggerDefs::TriggerId::photonTId)
    {
      _tree->Branch("ll1_4x4_overlappingsum", &b_ll1_4x4_overlappingsum);
      _tree->Branch("ll1_4x4_overlappingsum_eta", &b_ll1_4x4_overlappingsum_eta);
      _tree->Branch("ll1_4x4_overlappingsum_phi", &b_ll1_4x4_overlappingsum_phi);
      _tree->Branch("ll1_4x4_overlappingsum_sample", &b_ll1_4x4_overlappingsum_sample);
    }
  _tree->Branch("ll1_2x2_emcal", &b_ll1_2x2_emcal);
  _tree->Branch("ll1_2x2_emcal_eta", &b_ll1_2x2_emcal_eta);
  _tree->Branch("ll1_2x2_emcal_phi", &b_ll1_2x2_emcal_phi);
  _tree->Branch("ll1_2x2_emcal_sample", &b_ll1_2x2_emcal_sample);

  _tree->Branch("hcalin_energy_sim",&b_hcalin_energy_sim);
  _tree->Branch("hcalin_phibin_sim",&b_hcalin_phibin_sim);
  _tree->Branch("hcalin_etabin_sim",&b_hcalin_etabin_sim);
  _tree->Branch("hcalout_energy_sim",&b_hcalout_energy_sim);
  _tree->Branch("hcalout_phibin_sim",&b_hcalout_phibin_sim);
  _tree->Branch("hcalout_etabin_sim",&b_hcalout_etabin_sim);
  _tree->Branch("emcal_energy_sim",&b_emcal_energy_sim);
  _tree->Branch("emcal_phibin_sim",&b_emcal_phibin_sim);
  _tree->Branch("emcal_etabin_sim",&b_emcal_etabin_sim);

  _tree->Branch("hcalin_energy_raw",&b_hcalin_energy_raw);
  _tree->Branch("hcalin_phibin_raw",&b_hcalin_phibin_raw);
  _tree->Branch("hcalin_etabin_raw",&b_hcalin_etabin_raw);
  _tree->Branch("hcalout_energy_raw",&b_hcalout_energy_raw);
  _tree->Branch("hcalout_phibin_raw",&b_hcalout_phibin_raw);
  _tree->Branch("hcalout_etabin_raw",&b_hcalout_etabin_raw);
  _tree->Branch("emcal_energy_raw",&b_emcal_energy_raw);
  _tree->Branch("emcal_phibin_raw",&b_emcal_phibin_raw);
  _tree->Branch("emcal_etabin_raw",&b_emcal_etabin_raw);

  _tree->Branch("hcalin_energy_calib",&b_hcalin_energy_calib);
  _tree->Branch("hcalin_phibin_calib",&b_hcalin_phibin_calib);
  _tree->Branch("hcalin_etabin_calib",&b_hcalin_etabin_calib);
  _tree->Branch("hcalout_energy_calib",&b_hcalout_energy_calib);
  _tree->Branch("hcalout_phibin_calib",&b_hcalout_phibin_calib);
  _tree->Branch("hcalout_etabin_calib",&b_hcalout_etabin_calib);
  _tree->Branch("emcal_energy_calib",&b_emcal_energy_calib);
  _tree->Branch("emcal_phibin_calib",&b_emcal_phibin_calib);
  _tree->Branch("emcal_etabin_calib",&b_emcal_etabin_calib);

  _tree->Branch("hcalin_energy_wave",&b_hcalin_energy_wave);
  _tree->Branch("hcalin_phibin_wave",&b_hcalin_phibin_wave);
  _tree->Branch("hcalin_etabin_wave",&b_hcalin_etabin_wave);
  _tree->Branch("hcalout_energy_wave",&b_hcalout_energy_wave);
  _tree->Branch("hcalout_phibin_wave",&b_hcalout_phibin_wave);
  _tree->Branch("hcalout_etabin_wave",&b_hcalout_etabin_wave);
  _tree->Branch("emcal_energy_wave",&b_emcal_energy_wave);
  _tree->Branch("emcal_phibin_wave",&b_emcal_phibin_wave);
  _tree->Branch("emcal_etabin_wave",&b_emcal_etabin_wave);

  _tree->Branch("_truth_particle_n", &b_truth_particle_n);
  _tree->Branch("_truth_particle_pid",&b_truth_particle_pid);
  _tree->Branch("_truth_particle_pt",&b_truth_particle_pt);
  _tree->Branch("_truth_particle_eta",&b_truth_particle_eta);
  _tree->Branch("_truth_particle_phi",&b_truth_particle_phi);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloEmulatorTreeMaker::InitRun(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloEmulatorTreeMaker::process_event(PHCompositeNode *topNode)
{

  if (Verbosity()) std::cout << __FILE__ << " "<< __LINE__<<" "<<std::endl;

  PHG4TruthInfoContainer *truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!truth_info)
    {

      std::cout << "no truth info" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }


  map<int, PHG4Particle*>::const_iterator particle_iter;

  PHG4TruthInfoContainer::ConstRange primary_range =
    truth_info->GetPrimaryParticleRange();

  for (PHG4TruthInfoContainer::ConstIterator particle_iter = primary_range.first;
       particle_iter != primary_range.second; ++particle_iter)
    {
      PHG4Particle *particle = particle_iter->second;
      b_truth_particle_pid.push_back(particle->get_pid());

      float eta = 0.5*log((particle->get_e()+particle->get_pz())/
  			  (particle->get_e()-particle->get_pz()));
      float pt = sqrt(TMath::Power(particle->get_px(),2)+TMath::Power(particle->get_py(), 2));
      float phi = atan2(particle->get_py(), particle->get_px());

      b_truth_particle_pt.push_back(pt);
      b_truth_particle_eta.push_back(eta);
      b_truth_particle_phi.push_back(phi);
      b_truth_particle_n++;
    }
  if (_triggerid == TriggerDefs::TriggerId::jetTId)
    {
      _ll1out = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_JET");
    }
  else if (_triggerid == TriggerDefs::TriggerId::photonTId)
    {
      _ll1out = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_PHOTON");
    }
  else 
    {
      _ll1out = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_NULL");
    }

  if (!_ll1out) 
    {
      std::cout << "No LL1Out Calo node... " << std::endl;
      exit(1);
    }

  _trigger_primitive_container = _ll1out->GetTriggerPrimitiveContainer();

  if (!_trigger_primitive_container)
    {
      cout << "no jet primitive_container" <<endl;
      exit(1);
    }

  TriggerPrimitiveContainerv1::Range range = _trigger_primitive_container->getTriggerPrimitives();
  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
    {
        
      _trigger_primitive = (*iter).second;
      uint16_t primid = TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey((*iter).first);
      TriggerPrimitive::Range srange = _trigger_primitive->getSums();
      for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
	{
	  uint16_t sumid = TriggerDefs::getSumLocId((*siter).first);
	  int peak = 0;
	  int i = 0, sample = -99;

	  for (auto isp = (*siter).second->begin(); isp != (*siter).second->end(); ++isp, i++)
	    {
	      if (peak < (int)(*isp) ) 
		{
		  peak = (int)(*isp);
		  sample = i;
		}
	    }
		  
	  if (_triggerid == TriggerDefs::TriggerId::jetTId)
	    {
	      b_ll1_8x8_nonoverlappingsum.push_back(peak);
	      
	      b_ll1_8x8_nonoverlappingsum_eta.push_back(sumid%12);
	      b_ll1_8x8_nonoverlappingsum_phi.push_back(sumid/2 + primid*2);
	      b_ll1_8x8_nonoverlappingsum_sample.push_back(sample);
	    }
	  else if (_triggerid == TriggerDefs::TriggerId::photonTId)
	    {
	      if (peak > 0)
		{
		  b_ll1_4x4_overlappingsum.push_back(peak);
	      
		  b_ll1_4x4_overlappingsum_eta.push_back(sumid%12);
		  b_ll1_4x4_overlappingsum_phi.push_back(sumid/2 + primid*2);
		  b_ll1_4x4_overlappingsum_sample.push_back(sample);
		}
	    }
	}
    }
  _ll1out = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_EMCAL");

  if (!_ll1out) 
    {
      std::cout << "No LL1Out Calo node... " << std::endl;
      exit(1);
    }

  _trigger_primitive_container = _ll1out->GetTriggerPrimitiveContainer();

  if (!_trigger_primitive_container)
    {
      cout << "no jet primitive_container" <<endl;
      exit(1);
    }

  range = _trigger_primitive_container->getTriggerPrimitives();
  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
    {
        
      _trigger_primitive = (*iter).second;
      uint16_t primid = TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey((*iter).first);
      TriggerPrimitive::Range srange = _trigger_primitive->getSums();
      for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
	{
	  uint16_t sumid = TriggerDefs::getSumLocId((*siter).first);
	  int peak = 0;
	  int i = 0, sample = -99;

	  for (auto isp = (*siter).second->begin(); isp != (*siter).second->end(); ++isp, i++)
	    {
	      if (peak < (int)(*isp) ) 
		{
		  peak = (int)(*isp);
		  sample = i;
		}
	    }
		  
	  b_ll1_2x2_emcal.push_back(peak);
	  b_ll1_2x2_emcal_eta.push_back(4*(primid%12) + sumid%4);
	  b_ll1_2x2_emcal_phi.push_back(sumid/4 + 4*(primid/12));
	  b_ll1_2x2_emcal_sample.push_back(sample);
	}
    }
    
  int size;  

  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_SIM_HCALIN");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_hcalin_energy_sim.push_back(energy);
      b_hcalin_etabin_sim.push_back(ieta);
      b_hcalin_phibin_sim.push_back(iphi);
    }

  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_SIM_HCALOUT");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_hcalout_energy_sim.push_back(energy);
      b_hcalout_etabin_sim.push_back(ieta);
      b_hcalout_phibin_sim.push_back(iphi);
    }
  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_SIM_CEMC");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_emcal_energy_sim.push_back(energy);
      b_emcal_etabin_sim.push_back(ieta);
      b_emcal_phibin_sim.push_back(iphi);
    }

  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_HCALIN");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);

      int16_t ped = _tower->get_waveform_value(0);
      int16_t peak = ped;
      for (int i = 0; i < 12; i ++)
	{
	  if (peak < _tower->get_waveform_value(i) ) peak = _tower->get_waveform_value(i);
	}
      float energy = static_cast<float> (peak - ped);
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_hcalin_energy_raw.push_back(energy);
      b_hcalin_etabin_raw.push_back(ieta);
      b_hcalin_phibin_raw.push_back(iphi);
    }
  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_HCALOUT");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);

      int16_t ped = _tower->get_waveform_value(0);
      int16_t peak = ped;
      for (int i = 0; i < 12; i ++)
	{
	  if (peak < _tower->get_waveform_value(i) ) peak = _tower->get_waveform_value(i);
	}
      float energy = static_cast<float> (peak - ped);
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_hcalout_energy_raw.push_back(energy);
      b_hcalout_etabin_raw.push_back(ieta);
      b_hcalout_phibin_raw.push_back(iphi);
    }
  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_CEMC");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);

      int16_t ped = _tower->get_waveform_value(0);
      int16_t peak = ped;
      for (int i = 0; i < 12; i ++)
	{
	  if (peak < _tower->get_waveform_value(i) ) peak = _tower->get_waveform_value(i);
	}
      float energy = static_cast<float> (peak - ped);
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_emcal_energy_raw.push_back(energy);
      b_emcal_etabin_raw.push_back(ieta);
      b_emcal_phibin_raw.push_back(iphi);
    }

  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_hcalin_energy_calib.push_back(energy);
      b_hcalin_etabin_calib.push_back(ieta);
      b_hcalin_phibin_calib.push_back(iphi);
    }
  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_hcalout_energy_calib.push_back(energy);
      b_hcalout_etabin_calib.push_back(ieta);
      b_hcalout_phibin_calib.push_back(iphi);
    }
  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_emcal_energy_calib.push_back(energy);
      b_emcal_etabin_calib.push_back(ieta);
      b_emcal_phibin_calib.push_back(iphi);
    }
  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERSWAVEFORM_CALIB_HCALIN");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_hcalin_energy_wave.push_back(energy);
      b_hcalin_etabin_wave.push_back(ieta);
      b_hcalin_phibin_wave.push_back(iphi);
    }
  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERSWAVEFORM_CALIB_HCALOUT");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_hcalout_energy_wave.push_back(energy);
      b_hcalout_etabin_wave.push_back(ieta);
      b_hcalout_phibin_wave.push_back(iphi);
    }
  if (Verbosity()) std::cout << __FUNCTION__<<__LINE__<<std::endl;
  _towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERSWAVEFORM_CALIB_CEMC");

  assert(_towers);

  size = _towers->size(); //online towers should be the same!
  for (int channel = 0; channel < size;channel++)
    {
      _tower = _towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      unsigned int towerkey = _towers->encode_key(channel);
      int ieta = _towers->getTowerEtaBin(towerkey);
      int iphi = _towers->getTowerPhiBin(towerkey);
	
      b_emcal_energy_wave.push_back(energy);
      b_emcal_etabin_wave.push_back(ieta);
      b_emcal_phibin_wave.push_back(iphi);
    }


  _tree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloEmulatorTreeMaker::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "CaloEmulatorTreeMaker::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }

  b_truth_particle_n = 0;
  b_truth_particle_pid.clear();
  b_truth_particle_pt.clear();
  b_truth_particle_eta.clear();
  b_truth_particle_phi.clear();


  b_hcalin_energy_sim.clear();
  b_hcalin_phibin_sim.clear();
  b_hcalin_etabin_sim.clear();
  b_hcalout_energy_sim.clear();
  b_hcalout_phibin_sim.clear();
  b_hcalout_etabin_sim.clear();
  b_emcal_energy_sim.clear();
  b_emcal_phibin_sim.clear();
  b_emcal_etabin_sim.clear();

  b_hcalin_energy_calib.clear();
  b_hcalin_phibin_calib.clear();
  b_hcalin_etabin_calib.clear();
  b_hcalout_energy_calib.clear();
  b_hcalout_phibin_calib.clear();
  b_hcalout_etabin_calib.clear();
  b_emcal_energy_calib.clear();
  b_emcal_phibin_calib.clear();
  b_emcal_etabin_calib.clear();

  b_hcalin_energy_wave.clear();
  b_hcalin_phibin_wave.clear();
  b_hcalin_etabin_wave.clear();
  b_hcalout_energy_wave.clear();
  b_hcalout_phibin_wave.clear();
  b_hcalout_etabin_wave.clear();
  b_emcal_energy_wave.clear();
  b_emcal_phibin_wave.clear();
  b_emcal_etabin_wave.clear();


  b_hcalin_energy_raw.clear();
  b_hcalin_phibin_raw.clear();
  b_hcalin_etabin_raw.clear();
  b_hcalout_energy_raw.clear();
  b_hcalout_phibin_raw.clear();
  b_hcalout_etabin_raw.clear();
  b_emcal_energy_raw.clear();
  b_emcal_phibin_raw.clear();
  b_emcal_etabin_raw.clear();

  b_ll1_8x8_nonoverlappingsum.clear();
  b_ll1_8x8_nonoverlappingsum_eta.clear();
  b_ll1_8x8_nonoverlappingsum_phi.clear();
  b_ll1_8x8_nonoverlappingsum_sample.clear();

  b_ll1_4x4_overlappingsum.clear();
  b_ll1_4x4_overlappingsum_eta.clear();
  b_ll1_4x4_overlappingsum_phi.clear();
  b_ll1_4x4_overlappingsum_sample.clear();

  b_ll1_2x2_emcal.clear();
  b_ll1_2x2_emcal_eta.clear();
  b_ll1_2x2_emcal_phi.clear();
  b_ll1_2x2_emcal_sample.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloEmulatorTreeMaker::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "CaloEmulatorTreeMaker::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloEmulatorTreeMaker::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "CaloEmulatorTreeMaker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  _f->Write();
  _f->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

