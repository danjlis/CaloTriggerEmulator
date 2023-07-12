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
  _idx = 12;
  for (unsigned int i = 0; i < 1024; i++)
    {
      m_l1_adc_table[i] = (i) & 0x3ff;
      m_l1_adc_table_time[i] = (i) & 0x3ff;
    }
  for (unsigned int i = 0; i < 4096; i++)
    {
      m_l1_slewing_table[i] = (i) & 0x1ff;
    }

  _ll1out = 0;
  _waveforms_cemc = 0;
  _waveforms_hcalin = 0;
  _waveforms_hcalout = 0;
  _waveforms_mbd = 0;

  _primitives = 0;
  _primitives_cemc = 0;
  _primitives_hcalin = 0;
  _primitives_hcalout = 0;

  _primitive = 0;
  _n_primitives = 0;
  _n_sums = 16;
  _m_trig_sub_delay = 6;
  _m_threshold = 10;

  m_nhit1 = 2;
  m_nhit2 = 10;
  m_timediff1 = 10;
  m_timediff2 = 20;
  m_timediff3 = 30;

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
  _do_mbd = false;
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
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  if (_verbose) std::cout << __FUNCTION__ << std::endl;

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
      _do_mbd = false;
    }
  else if (_triggerid == TriggerDefs::TriggerId::mbdTId)
    {
      if (_verbose) std::cout << "Using MBD Trigger."<<std::endl;
      _do_cemc = false;
      _do_hcalin = false;
      _do_hcalout = false;      
      _do_mbd = true;
    }
  else if (_triggerid == TriggerDefs::TriggerId::cosmicTId)
    {
      if (_verbose) std::cout << "Using Cosmic Trigger."<<std::endl;
      _do_cemc = false;
      _do_hcalin = true;
      _do_hcalout = true;
      _do_mbd = false;
    }
  else if (_triggerid == TriggerDefs::TriggerId::pairTId)
    {
      if (_verbose) std::cout << "Using Pair Trigger."<<std::endl;
      _do_cemc = true;
      _do_hcalin = false;
      _do_hcalout = false;
      _do_mbd = false;
    }
  else
    {
      cout << __FUNCTION__ << " : No trigger selected "<<endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }


  _ll1_nodename = "LL1OUT_" + _trigger;

  // Files build


  for ( auto iter = _m_det_map[_triggerid].begin() ; iter != _m_det_map[_triggerid].end() ; ++iter)
    { 
      for (int i = 0; i < _m_prim_map[TriggerDefs::GetDetectorId(*iter)]; i++)
	{
	  if (strcmp((*iter).c_str(), "MBD") == 0)
	    {

	      h_nhit = new TH1D(Form("h_nhit_%d",i),"", 33, -0.5, 32.5);
	      h2_line_up = new TH2D(Form("h_line_up_%d",i),"", m_nsamples, -0.5, static_cast<float>(m_nsamples) - 0.5, 60, -0.5, 59.5);

	      hm->registerHisto(h_nhit);
	      hm->registerHisto(h2_line_up);

	      v_nhit.push_back(h_nhit);
	      v_line_up.push_back(h2_line_up);
	      for (int j = 0; j < 8; j++)
		{

		  h_mbd_charge = new TH1D(Form("h_mbd_charge_%d_%d", j, i), "", 257, 0, 257);
		  v_mbd_charge[Form("h_mbd_charge_%d_%d", j, i)] = h_mbd_charge;
		  hm->registerHisto(h_mbd_charge);
		}	      
	      for (int j = 0; j < 4; j++)
		{
		  h_mbd_time = new TH1D(Form("h_mbd_time_%d_%d", j, i), "", 256, 0, 4096);
		  v_mbd_time[Form("h_mbd_time_%d_%d", j, i)] = h_mbd_time;
		  hm->registerHisto(h_mbd_time);
		}

	      continue;
	    }

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
  _waveforms_mbd->Reset();

  _primitives_cemc->Reset();
  _primitives_hcalin->Reset();
  _primitives_hcalout->Reset();
  _primitives->Reset();

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

  while (m_peak_sub_ped_mbd.begin() != m_peak_sub_ped_mbd.end())
    {
      delete m_peak_sub_ped_mbd.begin()->second;
      m_peak_sub_ped_mbd.erase(m_peak_sub_ped_mbd.begin());
    }

  if (!_do_mbd) return;

  while (_sum_mbd.begin() != _sum_mbd.end())
    {
      delete *(_sum_mbd.begin());
      _sum_mbd.erase(_sum_mbd.begin());
    }

  while (_word_mbd.begin() != _word_mbd.end())
    {
      delete *(_word_mbd.begin());
      _word_mbd.erase(_word_mbd.begin());
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
  if (_do_mbd)
    {  
      WaveformContainerv1::Range begin_end = _waveforms_mbd->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;

      int peak_sub_ped;
      int peak_id = -1;
      float peak_height = 0;
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
	      if (wave.at(i) > peak_height)
		{
		  peak_height = wave.at(i);
		  peak_id = i;
		}
	      peak_sub_ped = static_cast<int>(wave.at(_m_trig_sub_delay + i)) - static_cast<int>(wave.at(i));
	      if (peak_sub_ped < 0) peak_sub_ped = 0;
	      v_peak_sub_ped->push_back(peak_sub_ped);
	    }
	  m_peak_sub_ped_mbd[ij] = v_peak_sub_ped;

	  v_line_up.at(ij/64)->Fill(peak_id, ij%64, 1);
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
  if (_do_mbd)
    {
      if (_verbose) std::cout << __FUNCTION__<<" "<<__LINE__<<"Gathering MBD primitives."<<std::endl;

      ip = 0;
      _n_primitives = _m_prim_map[TriggerDefs::DetectorId::mbdDId];

      for (i = 0; i < _n_primitives; i++, ip++)
	{

	  TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("MBD"), TriggerDefs::GetPrimitiveId("MBD"), ip);
	  _primitive = new TriggerPrimitive(primkey);
	  mask = CheckFiberMasks(primkey);
	  std::vector<unsigned int> *sum_mbd = nullptr;;

	  _sum_mbd.clear();
	  for (int j = 0 ; j < 13; j++)
	    {
	      sum_mbd = new std::vector<unsigned int>();
	      _sum_mbd.push_back(sum_mbd);
	    }
	      
	  for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++) 
	    {

	      for (int j = 0; j < 8; j++)
		{
		  m_trig_charge[j] = 0;
		}
	      m_trig_nhit = 0;
	      for (int j = 0; j < 4; j++)
		{
		  m_trig_time[j] = 0;
		}

	      unsigned int tmp, tmp2;  
	      unsigned int qadd[32];
	      for (int isec = 0; isec < 4; isec++)
		{
		  for (j = 0; j < 8;j++)
		    {
		      tmp = m_l1_adc_table[m_peak_sub_ped_mbd[ i*64 + 8 + isec*16 + j ]->at(is) >> 4];
	
		      qadd[isec*8+j] = (tmp & 0x380) >> 7;
		      
		      m_trig_charge[isec*2 + j/4] += tmp & 0x7ff;
		      
		    }
		}
	      for (int isec = 0; isec < 4; isec++)
		{
      
		  for (j = 0; j < 8;j++)
		    {
		      tmp = m_l1_adc_table[m_peak_sub_ped_mbd[ i*64 + isec*16 + j ]->at(is) >> 4];
		      
		      m_trig_nhit += (tmp & 0x200) >> 9;
		      
		      tmp2 = m_l1_slewing_table[(qadd[isec*8+j] << 9) + (tmp & 0x01ff)];
		      
		      m_trig_time[isec] += tmp2;
		      
		    }
		}
	      for (int j = 0; j < 8; j++)
		{
		  _sum_mbd[j]->push_back(m_trig_charge[j]);
		}
	      _sum_mbd[8]->push_back(m_trig_nhit);

	      for (int j = 0; j < 4; j++)
		{
		  _sum_mbd[9+j]->push_back(m_trig_time[j]);
		}

	      if (is == _idx)
		{

		  for (int j = 0; j < 8; j++)
		    {
		      v_mbd_charge[Form("h_mbd_charge_%d_%d", j, ip)]->Fill(m_trig_charge[j]);
		    }

		  v_nhit[ip]->Fill(m_trig_nhit);
		  		  
		  for (int j = 0; j < 4; j++)
		    {
		      v_mbd_time[Form("h_mbd_time_%d_%d", j, ip)]->Fill(m_trig_time[j]);
		    }
		}
	    }
	  for (int j = 0; j < 13; j++)
	    {
	      TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("MBD"), TriggerDefs::GetPrimitiveId("MBD"), ip, j);
	      _primitive->add_sum(sumkey, _sum_mbd[j]);
	    }
	  _primitives->add_primitive(primkey, _primitive);
	}

      if (_verbose) std::cout << "Total primitives in mbd: "<<_primitives->size()<<std::endl; 

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
  else if (_triggerid == TriggerDefs::TriggerId::mbdTId)
    {

      if (_verbose) 
	{
	  std::cout <<__FUNCTION__<<" "<<__LINE__<<" processing MBD trigger , bits before: "<< _bits->size();
	}
      if (!_primitives)
	{
	  std::cout << "There is no primitive container" << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}




      TriggerPrimitiveContainerv1::Range range;      
      TriggerPrimitive::Range sumrange;
      int ip, isum, i, j;

      range = _primitives->getTriggerPrimitives();
      
      if (_verbose) std::cout << __FUNCTION__<<" "<<__LINE__<<" mbd primitives size: "<<_primitives->size()<<std::endl; 

      std::vector<unsigned int> *word_mbd = nullptr;

      _word_mbd.clear();
      for (int j = 0 ; j < 8; j++)
	{
	  word_mbd = new std::vector<unsigned int>();
	  _word_mbd.push_back(word_mbd);
	}


      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
	{
	  ip = 0;
	  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter, ip++)
	    {
	      _primitive  = (*iter).second;
	      sumrange = _primitive->getSums(); 
	      isum = 0;
	      for (TriggerPrimitive::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum, isum++)
		{
		  if (isum < 8)
		    {
		      m2_trig_charge[ip][isum] = ((*iter_sum).second)->at(is);
		    }
		  else if (isum == 8)  m2_trig_nhit[ip] = ((*iter_sum).second)->at(is);
		  else 
		    {
		      m2_trig_time[ip][isum - 9] = ((*iter_sum).second)->at(is);
		    }
		}
	    }

	  if (_verbose && is == 11) {
	    for (int q = 0; q < 8; q++)
	      {

		std::cout <<"Q"<<std::dec<<q<<": "; 
		for (int ipp = 0; ipp < 4; ipp++) std::cout << std::hex<<m2_trig_charge[ipp][q]<<" ";
		std::cout <<" "<<std::endl;
	      }
	    std::cout <<"NH: "; 
	    for (int ipp = 0; ipp < 4; ipp++) std::cout << std::hex<<m2_trig_nhit[ipp]<<" ";
	    std::cout <<" "<<std::endl;

	    for (int q = 0; q < 4; q++)
	      {

		std::cout <<"T"<<std::dec<<q<<": "; 
		for (int ipp = 0; ipp < 4; ipp++) std::cout << std::hex<<m2_trig_time[ipp][q]<<" ";
		std::cout <<" "<<std::endl;
	      }

	  }

	  m_out_tsum[0] = 0;
	  m_out_tsum[1] = 0;
	  m_out_nhit[0] = 0;
	  m_out_nhit[1] = 0;
	  m_out_tavg[0] = 0;
	  m_out_tavg[1] = 0;
	  m_out_trem[0] = 0;
	  m_out_trem[1] = 0;
	  m_out_vtx_sub = 0;
	  m_out_vtx_add = 0;

	  for (i = 0; i < 2; i++)
	    {
	      for (j = 0 ; j < 4 ; j++)
		{
		  m_out_tsum[0] += m2_trig_time[i][j];
		  m_out_tsum[1] += m2_trig_time[i+2][j];
		}
	      m_out_nhit[0] += m2_trig_nhit[i];
	      m_out_nhit[1] += m2_trig_nhit[i+2];
	    }

	  if (m_out_nhit[0] != 0) 
	    {
	      m_out_tavg[0] = m_out_tsum[0]/m_out_nhit[0];
	      m_out_trem[0] = m_out_tsum[0]%m_out_nhit[0]; 
	    }
	  if (m_out_nhit[1] != 0) 
	    {
	      m_out_tavg[1] = m_out_tsum[1]/m_out_nhit[1];
	      m_out_trem[1] = m_out_tsum[1]%m_out_nhit[1]; 
	    }

	  unsigned int max = m_out_tavg[0];
	  unsigned int min = m_out_tavg[1];
	  if (min > max) 
	    {
	      max = m_out_tavg[1];
	      min = m_out_tavg[0];
	    }
      
	  m_out_vtx_sub = (max - min) & 0x1ff;
	  m_out_vtx_add = (m_out_tavg[0] + m_out_tavg[1]) & 0x3ff;

	  _word_mbd[0]->push_back(m_out_tavg[0]);
	  _word_mbd[1]->push_back(m_out_tavg[1]);
	  _word_mbd[2]->push_back(m_out_nhit[0]);
	  _word_mbd[3]->push_back(m_out_nhit[1]);
	  _word_mbd[4]->push_back(m_out_trem[0]);
	  _word_mbd[5]->push_back(m_out_trem[1]);
	  _word_mbd[6]->push_back(m_out_vtx_sub);
	  _word_mbd[7]->push_back(m_out_vtx_add);

	  if (m_out_nhit[0] >= m_nhit1) bits.at(is) ^= 1 << 0;
	  if (m_out_nhit[1] >= m_nhit1) bits.at(is) ^= 1 << 1;
	  if (m_out_nhit[0] >= m_nhit2) bits.at(is) ^= 1 << 2;
	  if (m_out_nhit[1] >= m_nhit2) bits.at(is) ^= 1 << 3;
      
	  if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff1) bits.at(is) ^= 1 << 4;
	  if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff2) bits.at(is) ^= 1 << 5;
	  if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff3) bits.at(is) ^= 1 << 6;
	  if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff1) bits.at(is) ^= 1 << 7;
	  if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff2) bits.at(is) ^= 1 << 8;
	  if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff3) bits.at(is) ^= 1 << 9;

	  if (_verbose )   {
	    std::cout << "Trigger Word : "<<std::bitset<16>(bits.at(is)) << std::dec<<std::endl;  
	  }
	}

      _bits->clear();
      if (_verbose) std::cout << "bits after: ";
      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
	{
	  if (_verbose)std::cout <<" "<<bits.at(is);
	  _bits->push_back(bits.at(is));
	}
      if (_verbose)std::cout <<" "<<std::endl;
      
      for (int is = 0; is < 8; is++)
	{
	  _ll1out->add_word(is, _word_mbd[is]);
	}
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
  if (_do_mbd)
    { 
      _waveforms_mbd = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_MBD");

      if (!_waveforms_mbd) 
	{
	  std::cout << "No HCAL Waveforms found... " << std::endl;
	  exit(1);
	}

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
