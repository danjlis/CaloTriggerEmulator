#include "MBDTriggerEmulator.h"

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
#include "LL1Outv1.h"
#include <Event/Event.h>
#include <Event/packet.h>


using namespace std;

MBDTriggerEmulator::MBDTriggerEmulator(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , outfilename(filename)
  , hm(nullptr)
  , outfile(nullptr)
{

  _verbose = 0;
  _trigger = "NONE";
  _nevent = 0;

  m_nhit1 = 2;
  m_nhit2 = 5;
  m_timediff1 = 10;
  m_timediff2 = 20; 
  m_timediff3 = 30;

  for (unsigned int i = 0; i < 1024; i++)
    {
      m_l1_adc_table[i] = (i) & 0x3ff;
    }

  for (unsigned int i = 0; i < 4096; i++)
    {
      m_l1_slewing_table[i] = (i) & 0x1ff;
    }
}

MBDTriggerEmulator::~MBDTriggerEmulator()
{
  delete hm;
}

int MBDTriggerEmulator::Init(PHCompositeNode* topNode)
{
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  if (_verbose) std::cout << __FUNCTION__ << std::endl;
  return 0;
}

int MBDTriggerEmulator::InitRun(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FUNCTION__ << std::endl;
  CreateNodes(topNode);
  return 0;
}
int MBDTriggerEmulator::process_event(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FUNCTION__ << ": event " <<_nevent<<std::endl;

  GetNodes(topNode);

  reset_vars();

  process_waveforms();
  
  process_primitives();
    
  _nevent++;

  if (_verbose) std::cout << __FUNCTION__ << "end event"<<std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void MBDTriggerEmulator::reset_vars()
{
  for (int i = 0; i < 256; i++)
    {
      for (int j = 0; j < 25; j++)
	{
	  m_peak_sub_ped[i][j] = 0; 
	}      
      for (int j = 0; j < 31; j++)
	{
	  m_waveforms_mbd[i][j] = 0;
	}
    }
  for (int i = 0; i < 4; i++)
    {
      for (int j = 0; j < 8; j++)
	{
	  m_trig_charge[i][j] = 0;
	}
      m_trig_nhit[i] = 0;
      for (int j = 0; j < 4; j++)
	{
	  m_trig_time[i][j] = 0;
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

  m_trigger_word = 0;

  
}

int MBDTriggerEmulator::process_waveforms()
{

  if (_verbose) std::cout << __FILE__<<"::"<<__FUNCTION__ <<"::"<<__LINE__<< std::endl;

  WaveformContainerv1::Range begin_end = _waveforms_mbd->getWaveforms();
  WaveformContainerv1::Iter iwave = begin_end.first;
  int j = 0;

  for (; iwave != begin_end.second; ++iwave)
    {
      std::vector<int> wave = *(iwave->second);
      for (int i = 0; i < m_nsamples; i++)
	{
	  m_waveforms_mbd[iwave->first][i] = static_cast<int>(wave.at(i));
	  if (m_waveforms_mbd[iwave->first][i] > (0x3fff)) m_waveforms_mbd[iwave->first][i] = 0x3fff;
	}
      for (int i = 0; i < m_nsamples - 6;i++)
	{
	  m_peak_sub_ped[iwave->first][i] = m_waveforms_mbd[iwave->first][5 + i] - m_waveforms_mbd[iwave->first][i];
	  if (m_peak_sub_ped[iwave->first][i] < 0) m_peak_sub_ped[iwave->first][i] = 0;
	  if (_verbose > 3)
	    {
	      if (j%8 == 0) 
		{
		  std::cout << " " << std::endl; 
		  std::cout << " " << std::endl; 
		}
	      
	      if (j%64 == 0)
		{
		  std::cout <<" ADC: "<<(j/64)<<std::endl;
		}
	      
	      j++;
	      
	      std::cout << std::hex << m_peak_sub_ped[iwave->first] << " ";
	    }
	}
      if (_verbose > 3)  std::cout <<" "<<std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}
int MBDTriggerEmulator::process_primitives()
{

  if (_verbose) std::cout << __FILE__<<"::"<<__FUNCTION__ <<"::"<<__LINE__<< std::endl;

  unsigned int i, j, is;

  for (is = 0; is < 25; is++)
    {
      for (i = 0; i < 4; i++)
	{
	  unsigned int tmp, tmp2;  
	  unsigned int qadd[32];
	  if (_verbose > 3) std::cout << "i: "<<i<<": "<<std::endl;
	  
	  for (int isec = 0; isec < 4; isec++)
	    {
	      for (j = 0; j < 8;j++)
		{
		  tmp = m_l1_adc_table[m_peak_sub_ped[ i*64 + 8 + isec*16 + j ][is] >> 4];
		  
		  qadd[j] = (tmp & 0x380) >> 7;
		  
		  m_trig_charge[i][j/4] += tmp & 0x7ff;
		  
		}
	    }
	  m_trig_nhit[i] = 0;
	  
	  for (j = 0; j < 8;j++)
	    {
	      tmp = m_l1_adc_table[m_peak_sub_ped[ i*64 + j ][is] >> 4];
	      
	      m_trig_nhit[i] += (tmp & 0x200) >> 9;
	      
	      tmp2 = m_l1_slewing_table[(qadd[j] << 9) + (tmp & 0x01ff)];
	      
	      m_trig_time[i][j/8] += tmp2;
	      if (_verbose > 3) std::cout << i*64 + j <<" : "<<(m_peak_sub_ped[ i*64 + j ][is] >> 4)<<" : "<< tmp <<" : "<< tmp2 <<" : "<<((tmp&0x200) >> 9)<< std::endl;
	      
	    }
	  
	}
      
      if (_verbose) std::cout << __FILE__<<"::"<<__FUNCTION__ <<"::"<<__LINE__<< std::endl;
      
      for (i = 0; i < 2; i++)
	{
	  for (j = 0 ; j < 4 ; j++)
	    {
	      m_out_tsum[0] += m_trig_time[i][j];
	      m_out_tsum[1] += m_trig_time[i+2][j];
	    }
	  m_out_nhit[0] += m_trig_nhit[i];
	  m_out_nhit[1] += m_trig_nhit[i+2];
	}
      for (i = 0; i < 2; i++)
	{
	  m_out_tavg[i] = 0;
	  m_out_trem[i] = 0;
	  if (m_out_nhit[i] == 0) continue;
	  m_out_tavg[i] = m_out_tsum[i]/m_out_nhit[i];
	  m_out_trem[i] = m_out_tsum[i]%m_out_nhit[i];
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
      
      if (_verbose) std::cout << __FILE__<<"::"<<__FUNCTION__ <<"::"<<__LINE__<< std::endl;
      // Construct Trigger Word
      if (_verbose)
	{ 
	  std::cout << "Trigger Primitives:" <<std::endl;
	  std::cout << "  ADC  | Q | NH | T |"<<std::endl;
	  for (i = 0; i < 4; i++)
	    {
	      std::cout << std::dec <<"  "<<i<<" | ";
	      for (j = 0; j < 8;j++) std::cout << std::hex << m_trig_charge[i][j] <<" ";
	      std::cout <<" | "<<m_trig_nhit[i] <<" | ";
	      for (j = 0; j < 4;j++) std::cout << std::hex << m_trig_time[i][j] <<" ";
	      std::cout << " | "<< std::endl;
	    }
	  std::cout << "Trigger Output: "<<std::endl;
	  std::cout << std::hex << "North : NHIT ="<<m_out_nhit[0]<< "    South: NHIT = "<<m_out_nhit[1]<<std::endl;
	  std::cout << std::hex << "        TSUM ="<<m_out_tsum[0]<< "           TSUM = "<<m_out_tsum[1]<<std::endl;
	  std::cout << std::hex << "        TAVG ="<<m_out_tavg[0]<< "           TAVG = "<<m_out_tavg[1]<<std::endl;
	  std::cout << std::hex << "        TREM ="<<m_out_trem[0]<< "           TREM = "<<m_out_trem[1]<<std::endl;
	  std::cout << std::hex << " "<<std::endl;
	  std::cout << std::hex << "        VTX_SUB ="<<m_out_vtx_sub<<std::endl;
	  std::cout << std::hex << "        VTX_ADD ="<<m_out_vtx_add<<std::endl;
	}
      
      
      std::vector<unsigned int> trigger_words;
      trigger_words.push_back(m_out_nhit[0]);
      trigger_words.push_back(m_out_nhit[1]);
      trigger_words.push_back(m_out_tsum[0]);
      trigger_words.push_back(m_out_tsum[1]);
      trigger_words.push_back(m_out_tavg[0]);
      trigger_words.push_back(m_out_tavg[1]);
      trigger_words.push_back(m_out_trem[0]);
      trigger_words.push_back(m_out_trem[1]);
      trigger_words.push_back(m_out_vtx_sub);
      trigger_words.push_back(m_out_vtx_add);
      
      _ll1_mbd->AddTriggerWords(trigger_words);
      
      m_trigger_word = 0;
      if (m_out_nhit[0] >= m_nhit1) m_trigger_word ^= 1 << 0;
      if (m_out_nhit[1] >= m_nhit1) m_trigger_word ^= 1 << 1;
      if (m_out_nhit[0] >= m_nhit2) m_trigger_word ^= 1 << 2;
      if (m_out_nhit[1] >= m_nhit2) m_trigger_word ^= 1 << 3;
      
      if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff1) m_trigger_word ^= 1 << 4;
      if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff2) m_trigger_word ^= 1 << 5;
      if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff3) m_trigger_word ^= 1 << 6;
      if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff1) m_trigger_word ^= 1 << 7;
      if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff2) m_trigger_word ^= 1 << 8;
      if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff3) m_trigger_word ^= 1 << 9;

      if (_verbose )   std::cout << "Trigger Word : "<<std::bitset<16>(m_trigger_word) << std::dec<<std::endl;  
      _ll1_mbd->AddTriggerBits(m_trigger_word);
    }  
  return Fun4AllReturnCodes::EVENT_OK;
}
void MBDTriggerEmulator::GetNodes(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FUNCTION__ << std::endl;
  _ll1_mbd = findNode::getClass<LL1Outv1>(topNode, "LL1Out_MBD");

  if (!_ll1_mbd) 
    {
      std::cout << "No LL1Out for MBD found... " << std::endl;
      exit(1);
    }

  _waveforms_mbd = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_BBC");

  if (!_waveforms_mbd) 
    {
      std::cout << "No BBC Waveforms found... " << std::endl;
      exit(1);
    }

  return;

}
void MBDTriggerEmulator::CreateNodes(PHCompositeNode* topNode)
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

  LL1Outv1 *ll1out = findNode::getClass<LL1Outv1>(ll1Node, "LL1Out_MBD");
  if (!ll1out)
    {
      ll1out = new LL1Outv1();
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, "LL1Out_MBD", "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
}

int MBDTriggerEmulator::End(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FUNCTION__ << std::endl;

  if (_verbose) std::cout << "Processed " << _nevent << " events. " << std::endl;

  outfile->cd();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");
  return 0;
}
