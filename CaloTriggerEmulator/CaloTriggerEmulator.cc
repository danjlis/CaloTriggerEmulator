#include "CaloTriggerEmulator.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <TFile.h>
#include <cassert>
#include <sstream>
#include <string>
#include <cstdint>
#include <bitset>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include "LL1Defs.h"
#include "LL1Outv2.h"

using namespace std;

CaloTriggerEmulator::CaloTriggerEmulator(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , outfilename(filename)
  , hm(nullptr)
  , outfile(nullptr)
{

  _verbose = 0;
  _trigger = "NONE";
  _nevent = 0;
  _npassed = 0;
  m_nsamples = 31;

  for (unsigned int i = 0; i < 1024; i++)
    {
      m_l1_adc_table[i] = (i) & 0x3ff;
    }

  _ll1out = 0;
  _waveforms_cemc = 0;
  _waveforms_hcalin = 0;
  _waveforms_hcalout = 0;

  _primitives = 0;
  _primitives_cemc = 0;
  _primitives_hcalin = 0;
  _primitives_hcalout = 0;

  _primitive = 0;
  _n_primitives = 0;
  _n_sums = 16;
  _m_trig_sub_delay = 6;
  _m_threshold = 10;
  _m_det_map[TriggerDefs::TriggerId::noneTId] = {};
  _m_det_map[TriggerDefs::TriggerId::jetTId] = {"CEMC", "HCALOUT", "HCALIN"};
  _m_det_map[TriggerDefs::TriggerId::mbdTId] = {"MBD"};
  _m_det_map[TriggerDefs::TriggerId::cosmicTId] = {"HCALIN", "HCALOUT"};
  _m_det_map[TriggerDefs::TriggerId::pairTId] = {"CEMC"};

  _m_prim_map[TriggerDefs::DetectorId::noneDId] = 0;
  _m_prim_map[TriggerDefs::DetectorId::cemcDId] = 384;
  _m_prim_map[TriggerDefs::DetectorId::hcalinDId] = 24;
  _m_prim_map[TriggerDefs::DetectorId::hcaloutDId] = 24;
  _m_prim_map[TriggerDefs::DetectorId::mbdDId] = 4;

  _do_cemc = false;
  _do_hcalin = false;
  _do_hcalout = false;

  _masks_channel = {69337149, 70385703};
  _masks_fiber = {69337136, 70385696};
}

CaloTriggerEmulator::~CaloTriggerEmulator()
{
  delete hm;
}

bool CaloTriggerEmulator::CheckChannelMasks(TriggerDefs::TriggerSumKey key)
{
  for (auto it = _masks_channel.begin(); it != _masks_channel.end(); ++it)
    {
      if (key == (*it)) return true;
    }
  return false;
}
bool CaloTriggerEmulator::CheckFiberMasks(TriggerDefs::TriggerPrimKey key)
{
  for (auto it = _masks_fiber.begin(); it != _masks_fiber.end(); ++it)
    {
      if (key == (*it)) return true;
    }
  return false;
}

void CaloTriggerEmulator::setTriggerType(const std::string &name)
{
  _trigger = name;
  _triggerid = TriggerDefs::GetTriggerId(_trigger);
  std::cout << "Setting Trigger type: "<<_trigger<< " (" <<_triggerid<<")"<<std::endl;  
}
int CaloTriggerEmulator::Init(PHCompositeNode* topNode)
{
  return 0;
}

int CaloTriggerEmulator::InitRun(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FUNCTION__ << std::endl;



  // Get number of primitives to construct;

  if (_triggerid == TriggerDefs::TriggerId::jetTId)
    {
      if (_verbose) std::cout << "Using Jet Trigger."<<std::endl;
      _do_cemc = true;
      _do_hcalin = true;
      _do_hcalout = true;
    }
  else if (_triggerid == TriggerDefs::TriggerId::mbdTId)
    {
      if (_verbose) std::cout << "Using MBD Trigger."<<std::endl;
      _do_cemc = false;
      _do_hcalin = false;
      _do_hcalout = false;      
    }
  else if (_triggerid == TriggerDefs::TriggerId::cosmicTId)
    {
      if (_verbose) std::cout << "Using Cosmic Trigger."<<std::endl;
      _do_cemc = false;
      _do_hcalin = true;
      _do_hcalout = true;
    }
  else if (_triggerid == TriggerDefs::TriggerId::pairTId)
    {
      if (_verbose) std::cout << "Using Pair Trigger."<<std::endl;
      _do_cemc = true;
      _do_hcalin = false;
      _do_hcalout = false;
    }
  else
    {
      cout << __FUNCTION__ << " : No trigger selected "<<endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }


  _ll1_nodename = "LL1OUT_" + _trigger;

  // Files build

  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  if (_verbose) std::cout << __FUNCTION__ << std::endl;


  _tree = new TTree("ttree"," a persevering date tree");

  for ( auto iter = _m_det_map[_triggerid].begin() ; iter != _m_det_map[_triggerid].end() ; ++iter)
    { 
      for (int i = 0; i < _m_prim_map[TriggerDefs::GetDetectorId(*iter)]; i++)
	{
	  peak_primitive = new TH2D(Form("peak_primitive_%s_%d", (*iter).c_str(), i), ";primitive;peak;counts", 16, -0.5, 15.5, m_nsamples - _m_trig_sub_delay, -0.5, m_nsamples - 1 - _m_trig_sub_delay);
	  avg_primitive = new TProfile(Form("avg_primitive_%s_%d",(*iter).c_str(),  i), ";primitive;avg", 16, -0.5, 15.5);
	  primitives = new TH2D(Form("primitives_%s_%d",(*iter).c_str(),  i), ";primitives;", 16, -0.5, 15.5, 64, 0, 256);
	  trigger_fire_map = new TH2D(Form("trigger_fire_map_%s_%d",(*iter).c_str(),  i), ";ch;ch", 4, -0.5, 3.5, 4, -0.5, 3.5);

	  if (strcmp((*iter).c_str(), "CEMC") ==0)
	    {
	      v_peak_primitive_cemc.push_back(peak_primitive);
	      v_avg_primitive_cemc.push_back(avg_primitive);
	      v_primitives_cemc.push_back(primitives);
	      v_trigger_fire_map_cemc.push_back(trigger_fire_map);
	    }
	  if (strcmp((*iter).c_str(), "HCALIN") ==0)
	    {
	      v_peak_primitive_hcalin.push_back(peak_primitive);
	      v_avg_primitive_hcalin.push_back(avg_primitive);
	      v_primitives_hcalin.push_back(primitives);
	      v_trigger_fire_map_hcalin.push_back(trigger_fire_map);
	    }
	  if (strcmp((*iter).c_str(), "HCALOUT") ==0)
	    {
	      v_peak_primitive_hcalout.push_back(peak_primitive);
	      v_avg_primitive_hcalout.push_back(avg_primitive);
	      v_primitives_hcalout.push_back(primitives);
	      v_trigger_fire_map_hcalout.push_back(trigger_fire_map);
	    }
	  hm->registerHisto(peak_primitive);
	  hm->registerHisto(avg_primitive);
	  hm->registerHisto(primitives);
	  hm->registerHisto(trigger_fire_map);
	}
    }

  CreateNodes(topNode);

  return 0;
}
int CaloTriggerEmulator::process_event(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FUNCTION__ << ": event " <<_nevent<<std::endl;

  GetNodes(topNode);

  if (process_waveforms() == Fun4AllReturnCodes::ABORTEVENT) return Fun4AllReturnCodes::ABORTEVENT;
  
  process_primitives();

  if (process_trigger() == Fun4AllReturnCodes::ABORTRUN) return Fun4AllReturnCodes::ABORTRUN;
  
  _nevent++;

  if (_verbose) identify();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTriggerEmulator::reset_vars()
{

  _waveforms_cemc->Reset();
  _waveforms_hcalin->Reset();
  _waveforms_hcalout->Reset();

  _primitives_cemc->Reset();
  _primitives_hcalin->Reset();
  _primitives_hcalout->Reset();

  _primitive->Reset();


  while (m_peak_sub_ped_cemc.begin() != m_peak_sub_ped_cemc.end())
    {
      delete m_peak_sub_ped_cemc.begin()->second;
      m_peak_sub_ped_cemc.erase(m_peak_sub_ped_cemc.begin());
    }

  while (m_peak_sub_ped_hcalin.begin() != m_peak_sub_ped_hcalin.end())
    {
      delete m_peak_sub_ped_hcalin.begin()->second;
      m_peak_sub_ped_hcalin.erase(m_peak_sub_ped_hcalin.begin());
    }

  while (m_peak_sub_ped_hcalout.begin() != m_peak_sub_ped_hcalout.end())
    {
      delete m_peak_sub_ped_hcalout.begin()->second;
      m_peak_sub_ped_hcalout.erase(m_peak_sub_ped_hcalout.begin());
    }

}

int CaloTriggerEmulator::process_waveforms()
{
  // Get range of waveforms
  
  if (_do_cemc)
    {  
      WaveformContainerv1::Range begin_end = _waveforms_cemc->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;
      
      
      int peak_sub_ped;
      std::vector<int> *v_peak_sub_ped;
      std::vector<int> wave;
      int ij = 0;
      for (; iwave != begin_end.second; ++iwave)
	{	  
	  wave = *(iwave->second);
	  v_peak_sub_ped = new std::vector<int>();
	  peak_sub_ped = 0;
	  
	  for (int i = 0; i < m_nsamples - _m_trig_sub_delay;i++)
	    {
	      peak_sub_ped = static_cast<int>(wave.at(_m_trig_sub_delay + i)) - static_cast<int>(wave.at(i));
	      if (peak_sub_ped < 0) peak_sub_ped = 0;
	      v_peak_sub_ped->push_back(peak_sub_ped);
	    }
	  m_peak_sub_ped_cemc[ij] = v_peak_sub_ped;
	  ij++;
	}
      if (_verbose) std::cout << "Processed waves: "<<ij <<std::endl;
      if (!ij) return Fun4AllReturnCodes::ABORTEVENT;
    }
  if (_do_hcalin)
    {
      WaveformContainerv1::Range begin_end = _waveforms_hcalin->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;
      
      
      int peak_sub_ped;
      std::vector<int> *v_peak_sub_ped;
      std::vector<int> wave;
      int ij = 0;
      for (; iwave != begin_end.second; ++iwave)
	{	  
	  wave = *(iwave->second);
	  v_peak_sub_ped = new std::vector<int>();
	  peak_sub_ped = 0;
	  
	  for (int i = 0; i < m_nsamples - _m_trig_sub_delay;i++)
	    {
	      peak_sub_ped = static_cast<int>(wave.at(_m_trig_sub_delay + i)) - static_cast<int>(wave.at(i));
	      if (peak_sub_ped < 0) peak_sub_ped = 0;
	      v_peak_sub_ped->push_back(peak_sub_ped);
	    }
	  m_peak_sub_ped_hcalin[ij] = v_peak_sub_ped;
	  ij++;
	}
      if (_verbose) std::cout << "Processed waves: "<<ij <<std::endl;
      if (!ij) return Fun4AllReturnCodes::ABORTEVENT;      
    }

  if (_do_hcalout)
    {
      WaveformContainerv1::Range begin_end = _waveforms_hcalout->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;
      
      
      int peak_sub_ped;
      std::vector<int> *v_peak_sub_ped;
      std::vector<int> wave;
      int ij = 0;
      for (; iwave != begin_end.second; ++iwave)
	{	  
	  wave = *(iwave->second);
	  v_peak_sub_ped = new std::vector<int>();
	  peak_sub_ped = 0;
	  
	  for (int i = 0; i < m_nsamples - _m_trig_sub_delay;i++)
	    {
	      peak_sub_ped = static_cast<int>(wave.at(_m_trig_sub_delay + i)) - static_cast<int>(wave.at(i));
	      if (peak_sub_ped < 0) peak_sub_ped = 0;
	      v_peak_sub_ped->push_back(peak_sub_ped);
	    }
	  m_peak_sub_ped_hcalout[ij] = v_peak_sub_ped;
	  ij++;
	}
      if (_verbose) std::cout << "Processed waves: "<<ij <<std::endl;
      if (!ij) return Fun4AllReturnCodes::ABORTEVENT;      
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
  
int CaloTriggerEmulator::process_primitives()
{
  
  if (_verbose) std::cout << __FILE__<<"::"<<__FUNCTION__ <<"::"<<__LINE__<< std::endl;
  
  unsigned int ip, j;
  int id_peak;
  unsigned int peak;  
  int i;
  bool mask;
  if (_do_cemc)
    {
      if (_verbose) std::cout << __FUNCTION__<<" "<<__LINE__<<"Gathering CEMC primitives."<<std::endl;
      ip = 0;
      _n_primitives = _m_prim_map[TriggerDefs::DetectorId::cemcDId];
      for (i = 0; i < _n_primitives; i++, ip++)
	{
	  unsigned int tmp;  
	  TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("CEMC"), TriggerDefs::GetPrimitiveId("CEMC"), ip);
	  _primitive = new TriggerPrimitive(primkey);
	  unsigned int sum;
	  mask = CheckFiberMasks(primkey);

	  for (int isum = 0; isum < _n_sums; isum++)
	    {
	      id_peak = -1;
	      peak = 0;
	      TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("CEMC"), TriggerDefs::GetPrimitiveId("CEMC"), ip, isum);
	      _sum = new std::vector<unsigned int>();
	      mask |= CheckChannelMasks(sumkey);
	      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
		{
		  sum = 0;
		  if (!mask)
		    {
		      for (j = 0; j < 4;j++)
			{
			  tmp = m_l1_adc_table[m_peak_sub_ped_cemc[64*ip + isum*4 + j]->at(is) >> 4];
			  sum += (tmp & 0x3ff);
			}
		      sum = (sum & 0x3ff) >> 2;
		      if (peak < sum) 
			{
			  peak = sum;
			  id_peak = is;
			}
		      _sum->push_back(sum);
		    }
		  
		  if (peak > _m_threshold)
		    {
		      v_avg_primitive_cemc.at(i)->Fill(isum, peak);
		      v_primitives_cemc.at(i)->Fill(isum, peak);
		      
		      v_trigger_fire_map_cemc.at(i)->Fill(isum%4, isum/4);
		      v_peak_primitive_cemc.at(i)->Fill(isum, id_peak);
		    }
		}
	      _primitive->add_sum(sumkey, _sum);
	      
	    }
	  
	  _primitives_cemc->add_primitive(primkey, _primitive);
	  if (_verbose) std::cout << "Total primitives in cemc: "<<_primitive->size()<<std::endl;	  	  
	}  

    }
  if (_do_hcalout)
    {
      if (_verbose) std::cout << __FUNCTION__<<" "<<__LINE__<<"Gathering HCALOU primitives."<<std::endl;
      ip = 0;
      _n_primitives = _m_prim_map[TriggerDefs::DetectorId::hcaloutDId];
      for (i = 0; i < _n_primitives; i++, ip++)
	{
	  unsigned int tmp;  
	  TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("HCALOUT"), TriggerDefs::GetPrimitiveId("HCALOUT"), ip);
	  _primitive = new TriggerPrimitive(primkey);
	  unsigned int sum;
	  mask = CheckFiberMasks(primkey);
	  for (int isum = 0; isum < _n_sums; isum++)
	    {
	      id_peak = -1;
	      peak = 0;
	      TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("HCALOUT"), TriggerDefs::GetPrimitiveId("HCALOUT"), ip, isum);
	      _sum = new std::vector<unsigned int>();
	      mask |= CheckChannelMasks(sumkey);
	      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
		{
		  sum = 0;
		  if (!mask)
		    {
		      for (j = 0; j < 4;j++)
			{
			  tmp = m_l1_adc_table[m_peak_sub_ped_hcalout[64*ip + isum*4 + j]->at(is) >> 4];
			  sum += (tmp & 0x3ff);
			}
		      sum = (sum & 0x3ff) >> 2;
		      if (peak < sum) 
			{
			  peak = sum;
			  id_peak = is;
			}
		      _sum->push_back(sum);
		    }
		  if (peak > _m_threshold)
		    {
		      v_peak_primitive_hcalout.at(i)->Fill(isum, id_peak);
		      v_avg_primitive_hcalout.at(i)->Fill(isum, peak);
		      v_primitives_hcalout.at(i)->Fill(isum, peak);
		      v_trigger_fire_map_hcalout.at(i)->Fill(isum%4, isum/4);
		    }
		}
	      _primitive->add_sum(sumkey, _sum);
	      
	    }
	  
	  _primitives_hcalout->add_primitive(primkey, _primitive);
	  if (_verbose) std::cout << "Total primitives in hcalout: "<<_primitive->size()<<std::endl;	  
	}  

    }
  if (_do_hcalin)
    {
      if (_verbose) std::cout << __FUNCTION__<<" "<<__LINE__<<"Gathering HCALIN primitives."<<std::endl;
      ip = 0;
      _n_primitives = _m_prim_map[TriggerDefs::DetectorId::hcalinDId];
      for (i = 0; i < _n_primitives; i++, ip++)
	{
	  unsigned int tmp;  
	  TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("HCALIN"), TriggerDefs::GetPrimitiveId("HCALIN"), ip);
	  _primitive = new TriggerPrimitive(primkey);
	  unsigned int sum;
	  mask = CheckFiberMasks(primkey);
	  for (int isum = 0; isum < _n_sums; isum++)
	    {
	      id_peak = -1;
	      peak = 0;
	      TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("HCALIN"), TriggerDefs::GetPrimitiveId("HCALIN"), ip, isum);
	      _sum = new std::vector<unsigned int>();
	      mask |= CheckChannelMasks(sumkey);
	      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
		{
		  sum = 0;
		  if (!mask)
		    {
		      for (j = 0; j < 4;j++)
			{
			  tmp = m_l1_adc_table[m_peak_sub_ped_hcalin[64*ip + isum*4 + j]->at(is) >> 4];
			  sum += (tmp & 0x3ff);
			}
		      sum = (sum & 0x3ff) >> 2;
		      if (peak < sum) 
			{
			  peak = sum;
			  id_peak = is;
			}
		    }
		  _sum->push_back(sum);
		}
	      if (peak > _m_threshold)
		{
		  v_peak_primitive_hcalin.at(i)->Fill(isum, id_peak);
		  v_avg_primitive_hcalin.at(i)->Fill(isum, peak);
		  v_primitives_hcalin.at(i)->Fill(isum, peak);

		  v_trigger_fire_map_hcalin.at(i)->Fill(isum%4, isum/4);
		}
	      _primitive->add_sum(sumkey, _sum);
	      
	    }
	  if (_verbose) std::cout << "Total primitives in hcalin: "<<_primitive->size()<<std::endl;	  
	  _primitives_hcalin->add_primitive(primkey, _primitive);
	  
	}  
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTriggerEmulator::process_trigger()
{
  std::vector<unsigned int> bits;

  for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
    {
      bits.push_back(0);
    }


  if (_triggerid == TriggerDefs::TriggerId::cosmicTId)
    {

      if (_verbose) 
	{
	  std::cout <<__FUNCTION__<<" "<<__LINE__<<" processing COSMIC trigger , bits before: "<< _bits->size();
	}
      if (!_primitives_hcalout || !_primitives_hcalin)
	{
	  std::cout << "There is no primitive container" << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      TriggerPrimitiveContainerv1::Range range;      
      range = _primitives_hcalout->getTriggerPrimitives();
      if (_verbose) std::cout << __FUNCTION__<<" "<<__LINE__<<" hcalout primitives size: "<<_primitives_hcalout->size()<<std::endl; 
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  TriggerDefs::TriggerPrimKey key = (*iter).first;
	  if (CheckFiberMasks(key)) {
	    if (_verbose) std::cout << "masked: "<<key<<std::endl;
	    continue;
	  }

	  _primitive  = (*iter).second;
	  TriggerPrimitive::Range sumrange = _primitive->getSums(); 
	  if (_verbose) std::cout << __FUNCTION__<<" "<<__LINE__<<" key: "<<key<<" size: "<<_primitive->size()<<std::endl; 
	  for (TriggerPrimitive::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
	    {
	      TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;
	      int i = 0;
	      if (CheckChannelMasks(sumkey)) continue;
	      if (_verbose) cout <<" sum " << sumkey << " size "<<(*iter_sum).second->size()<<std::endl;
	      for (auto it_s = (*iter_sum).second->begin(); it_s != (*iter_sum).second->end(); ++it_s)
		{
		  
		  if ((*it_s) >= _m_threshold) {
		    bits.at(i) |= 1;
		  }
		  i++;
		}
	    }
	  
	  
	}

      range = _primitives_hcalin->getTriggerPrimitives();
      if (_verbose) std::cout << __FUNCTION__<<" "<<__LINE__<<" hcalout primitives size: "<<_primitives_hcalin->size()<<std::endl; 
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  TriggerDefs::TriggerPrimKey key = (*iter).first;
	  if (CheckFiberMasks(key)) {
	    if (_verbose) std::cout << "masked: "<<key<<std::endl;
	    continue;
	  }
	  _primitive = (*iter).second;
	  if (_verbose) std::cout << __FUNCTION__<<" "<<__LINE__<<" key: "<<key<<" size: "<<_primitive->size()<<std::endl; 
	  TriggerPrimitive::Range sumrange = _primitive->getSums(); 
	  for (TriggerPrimitive::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
	    {
	      TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;
	      int i = 0;
	      if (CheckChannelMasks(sumkey)) {
		if (_verbose) std::cout << "masked: "<<sumkey<<std::endl;
		continue;
	      }
	      if (_verbose) cout <<" sum " << sumkey << " size "<<(*iter_sum).second->size()<<std::endl;
	      for (auto it_s = (*iter_sum).second->begin(); it_s != (*iter_sum).second->end(); ++it_s)
		{

		  if ((*it_s) >= _m_threshold){
		    bits.at(i) |= 1;
		  }
		  i++;
		}
	    }
	  
	  
	}
      _bits->clear();
      int pass = 0;
      if (_verbose) std::cout << "bits after: ";
      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
	{
	  if (_verbose)std::cout <<" "<<bits.at(is);
	  _bits->push_back(bits.at(is));
	  if (bits.at(is) == 1) pass = 1;
	}
      _npassed += pass;
      if (_verbose)std::cout <<" "<<std::endl;
    }
  else 
    {
      std::cout << "Trigger "<<_trigger<< " not implemented"<<std::endl; 
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTriggerEmulator::GetNodes(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FUNCTION__ << std::endl;
  _ll1out = findNode::getClass<LL1Outv2>(topNode, _ll1_nodename);

  if (!_ll1out) 
    {
      std::cout << "No LL1Out found... " << std::endl;
      exit(1);
    }

  _primitives = _ll1out->GetTriggerPrimitiveContainer();
  _bits = _ll1out->GetTriggerBits();

  if (_do_hcalout)
    { 
      _waveforms_hcalout = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALOUT");

      if (!_waveforms_hcalout) 
	{
	  std::cout << "No HCALOUT Waveforms found... " << std::endl;
	  exit(1);
	}
      _ll1out_hcalout = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_HCALOUT");

      if (!_ll1out_hcalout) 
	{
	  std::cout << "No HCAL Waveforms found... " << std::endl;
	  exit(1);
	}
      _primitives_hcalout = _ll1out_hcalout->GetTriggerPrimitiveContainer();
    }

  if (_do_hcalin)
    { 
      _waveforms_hcalin = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALIN");

      if (!_waveforms_hcalin) 
	{
	  std::cout << "No HCAL Waveforms found... " << std::endl;
	  exit(1);
	}

      _ll1out_hcalin = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_HCALIN");

      if (!_ll1out_hcalin) 
	{
	  std::cout << "No HCAL Waveforms found... " << std::endl;
	  exit(1);
	}
      _primitives_hcalin = _ll1out_hcalin->GetTriggerPrimitiveContainer();
    }
  if (_do_cemc)
    { 
      _waveforms_cemc = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_CEMC");

      if (!_waveforms_cemc) 
	{
	  std::cout << "No HCAL Waveforms found... " << std::endl;
	  exit(1);
	}

      _ll1out_cemc = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_CEMC");

      if (!_ll1out_cemc) 
	{
	  std::cout << "No HCAL Waveforms found... " << std::endl;
	  exit(1);
	}

      _primitives_cemc = _ll1out_cemc->GetTriggerPrimitiveContainer();
    }

  return;

}
void CaloTriggerEmulator::CreateNodes(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FUNCTION__ << std::endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    }

  PHCompositeNode *ll1Node = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "LL1"));
  if (!ll1Node)
    {
      ll1Node = new PHCompositeNode("LL1");
      dstNode->addNode(ll1Node);
    }

  LL1Outv2 *ll1out = findNode::getClass<LL1Outv2>(ll1Node, _ll1_nodename);
  if (!ll1out)
    {
      ll1out = new LL1Outv2(_trigger, "NONE");
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, _ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }

  if (_do_cemc)
    {
      std::string ll1_nodename = "LL1OUT_CEMC";
      LL1Outv2 *ll1out_d = findNode::getClass<LL1Outv2>(ll1Node, ll1_nodename);
      if (!ll1out_d)
	{
	  ll1out_d = new LL1Outv2(_trigger, "CEMC");
	  PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d, ll1_nodename, "PHObject");
	  ll1Node->addNode(LL1OutNode);
	}

    }
  if (_do_hcalout)
    {
      std::string ll1_nodename = "LL1OUT_HCALOUT";
      LL1Outv2 *ll1out_d = findNode::getClass<LL1Outv2>(ll1Node, ll1_nodename);
      if (!ll1out_d)
	{
	  ll1out_d = new LL1Outv2(_trigger, "HCALOUT");
	  PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d, ll1_nodename, "PHObject");
	  ll1Node->addNode(LL1OutNode);
	}

    }
  if (_do_hcalin)
    {
      std::string ll1_nodename = "LL1OUT_HCALIN";
      LL1Outv2 *ll1out_d = findNode::getClass<LL1Outv2>(ll1Node, ll1_nodename);
      if (!ll1out_d)
	{
	  ll1out_d = new LL1Outv2(_trigger, "HCALIN");
	  PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d, ll1_nodename, "PHObject");
	  ll1Node->addNode(LL1OutNode);
	}

    }
}

int CaloTriggerEmulator::End(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << "Processed " << _nevent << " events. " << std::endl;
  std::cout << "------------------------" <<std::endl;
  std::cout << "Total passed: " <<_npassed<<"/"<<_nevent <<std::endl;
  std::cout << "------------------------" <<std::endl;

  outfile->cd();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");

  return 0;
}

void CaloTriggerEmulator::identify()
{
  std::cout <<  " CaloTriggerEmulator: "<< _trigger << std::endl;
  std::cout <<  " LL1Out: "<< std::endl;
  if(_ll1out) _ll1out->identify();
  if(_ll1out_cemc) _ll1out_cemc->identify();
  if(_ll1out_hcalin) _ll1out_hcalin->identify();
  if(_ll1out_hcalout) _ll1out_hcalout->identify();
  std::cout << " Processed " << _nevent << std::endl;
}
