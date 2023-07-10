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
  m_nsamples = 31;

  for (unsigned int i = 0; i < 1024; i++)
    {
      m_l1_adc_table[i] = (i) & 0x3ff;
    }

  _ll1out = 0;
  _waveforms = 0;
  _primitives = 0;
  _primitive = 0;
  _n_primitives = 0;
  _n_sums = 16;
  _m_trig_sub_delay = 6;


  _n_prim_map["NONE"] = 0;
  _n_prim_map["EMCAL"] = 384;
  _n_prim_map["HCALIN"] = 24;
  _n_prim_map["HCALOUT"] = 24;

}

CaloTriggerEmulator::~CaloTriggerEmulator()
{
  delete hm;
}

int CaloTriggerEmulator::Init(PHCompositeNode* topNode)
{


  // Get number of primitives to construct;

  if (!_n_prim_map[_trigger])
    {
      cerr << __FUNCTION__ << " : No trigger selected "<<endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  _n_primitives = _n_prim_map[_trigger];
  _ll1_nodename = "LL1Out_" + _trigger;
  _waveform_nodename = "WAVEFORMS_" + _trigger; 

  v_primitives.reserve(_n_primitives);
  v_avg_primitive.reserve(_n_primitives);
  v_peak_primitive.reserve(_n_primitives);
  v_trigger_fire_map.reserve(_n_primitives);
  
  // Files build

  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  if (_verbose) std::cout << __FUNCTION__ << std::endl;


  _tree = new TTree("ttree"," a persevering date tree");

  
  for (int i = 0; i < _n_primitives; i++)
    {
      peak_primitive = new TH2D(Form("peak_primitive_%d", i), ";primitive;peak;counts", 16, -0.5, 15.5, m_nsamples - _m_trig_sub_delay, -0.5, m_nsamples - 1 - _m_trig_sub_delay);
      avg_primitive = new TProfile(Form("avg_primitive_%d", i), ";primitive;avg", 16, -0.5, 15.5);
      primitives = new TH2D(Form("primitives_%d", i), ";primitives;", 16, -0.5, 15.5, 64, 0, 256);
      trigger_fire_map = new TEfficiency(Form("trigger_fire_map_%d", i), ";ch;ch", 4, -0.5, 3.5, 4, -0.5, 3.5);

      v_peak_primitive.push_back(peak_primitive);
      v_avg_primitive.push_back(avg_primitive);
      v_primitives.push_back(primitives);
      v_trigger_fire_map.push_back(trigger_fire_map);

      hm->registerHisto(peak_primitive);
      hm->registerHisto(avg_primitive);
      hm->registerHisto(primitives);
      hm->registerHisto(trigger_fire_map);
    }

  return 0;
}

int CaloTriggerEmulator::InitRun(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FUNCTION__ << std::endl;
  CreateNodes(topNode);

  return 0;
}
int CaloTriggerEmulator::process_event(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FUNCTION__ << ": event " <<_nevent<<std::endl;

  GetNodes(topNode);

  if (process_waveforms() == Fun4AllReturnCodes::ABORTEVENT) return Fun4AllReturnCodes::ABORTEVENT;
  
  process_primitives();

    
  _nevent++;

  if (_verbose) identify();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTriggerEmulator::reset_vars()
{

  _waveforms->Reset();

  _primitives->Reset();

  _primitive->Reset();

  while (m_peak_sub_ped.begin() != m_peak_sub_ped.end())
    {
      delete m_peak_sub_ped.begin()->second;
      m_peak_sub_ped.erase(m_peak_sub_ped.begin());
    }
}

int CaloTriggerEmulator::process_waveforms()
{

  // Get range of waveforms
  WaveformContainerv1::Range begin_end = _waveforms->getWaveforms();
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
      m_peak_sub_ped[ij] = v_peak_sub_ped;
      ij++;
    }
  if (_verbose) std::cout << "Processed waves: "<<ij <<std::endl;
  if (!ij) return Fun4AllReturnCodes::ABORTEVENT;
  return Fun4AllReturnCodes::EVENT_OK;
}
int CaloTriggerEmulator::process_primitives()
{

  if (_verbose) std::cout << __FILE__<<"::"<<__FUNCTION__ <<"::"<<__LINE__<< std::endl;

  unsigned int ip, j;
  int id_peak;
  unsigned int peak;  
  int i;
  ip = 0;

  for (i = 0; i < _n_primitives; i++, ip++)
    {
      unsigned int tmp;  
      TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId("NONE"), TriggerDefs::GetDetectorId(_trigger), TriggerDefs::GetPrimitiveId(_trigger), ip);
      _primitive = new TriggerPrimitive(primkey);
      unsigned int sum;
      for (int isum = 0; isum < _n_sums; isum++)
	{
	  id_peak = -1;
	  peak = 0;
	  TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId("NONE"), TriggerDefs::GetDetectorId(_trigger), TriggerDefs::GetPrimitiveId(_trigger), ip, isum);
	  _sum = new std::vector<unsigned int>();
	  for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
	    {
	      sum = 0;
	      for (j = 0; j < 4;j++)
		{
		  tmp = m_l1_adc_table[m_peak_sub_ped[64*ip + isum*4 + j]->at(is) >> 4];
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
	  v_peak_primitive.at(i)->Fill(isum, id_peak);
	  v_avg_primitive.at(i)->Fill(isum, peak);
	  v_primitives.at(i)->Fill(isum, peak);
	  v_trigger_fire_map.at(i)->Fill(peak > 1, isum%4, isum/4);

	  _primitive->add_sum(sumkey, _sum);

	}

      _primitives->add_primitive(primkey, _primitive);

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

  _waveforms = findNode::getClass<WaveformContainerv1>(topNode, _waveform_nodename);

  if (!_waveforms) 
    {
      std::cout << "No HCAL Waveforms found... " << std::endl;
      exit(1);
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
      ll1out = new LL1Outv2("NONE", _trigger);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, _ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
}

int CaloTriggerEmulator::End(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << "Processed " << _nevent << " events. " << std::endl;

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
  std::cout <<  " Waveforms: "<< std::endl;
  if (_waveforms) _waveforms->identify();
  std::cout << " Processed " << _nevent << std::endl;
}
