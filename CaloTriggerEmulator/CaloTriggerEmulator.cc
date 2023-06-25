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
#include "LL1Outv1.h"
#include <Event/Event.h>
#include <Event/packet.h>


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
  m_nhit1 = 2;
  m_nhit2 = 5;

  for (unsigned int i = 0; i < 1024; i++)
    {
      m_l1_adc_table[i] = (i) & 0x3ff;
    }

}

CaloTriggerEmulator::~CaloTriggerEmulator()
{
  delete hm;
}

int CaloTriggerEmulator::Init(PHCompositeNode* topNode)
{
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  if (_verbose) std::cout << __FUNCTION__ << std::endl;


  _tree = new TTree("ttree"," a persevering date tree");
  
  //  _tree->Branch("trigger_primitives",m_trigger_primitives);
  
  for (int i = 0; i < 24; i++)
    {
      peak_primitive[i] = new TH2D(Form("peak_primitive_%d", i), ";primitive;peak;counts", 16, -0.5, 15.5, m_nsamples - 6, -0.5, m_nsamples - 5);
      avg_primitive[i] = new TProfile(Form("avg_primitive_%d", i), ";primitive;avg", 16, -0.5, 15.5);
      primitives[i] = new TH2D(Form("primitives_%d", i), ";primitives;", 16, -0.5, 15.5, 64, 0, 256);
      trigger_fire_map[i] = new TEfficiency(Form("trigger_fire_map_%d", i), ";ch;ch", 4, -0.5, 3.5, 4, -0.5, 3.5);

      hm->registerHisto(peak_primitive[i]);
      hm->registerHisto(avg_primitive[i]);
      hm->registerHisto(primitives[i]);
      hm->registerHisto(trigger_fire_map[i]);
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

  reset_vars();

  process_waveforms();
  
  process_primitives();
    
  _nevent++;

  if (_verbose) std::cout << __FUNCTION__ << "end event"<<std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTriggerEmulator::reset_vars()
{
  for (int i = 0; i < 64*24; i++)
    {
      for (int j = 0; j < m_nsamples  - 6; j++)
	{
	  m_peak_sub_ped[i][j] = 0; 
	}      
      for (int j = 0; j < m_nsamples; j++)
	{
	  m_waveforms_hcal[i][j] = 0;
	}
    }
  for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < m_nsamples - 6; j++)
	{
	  m_trig_sums[i][j] = 0; 
	}      
    }
}

int CaloTriggerEmulator::process_waveforms()
{

  if (_verbose) std::cout << __FILE__<<"::"<<__FUNCTION__ <<"::"<<__LINE__<< std::endl;

  WaveformContainerv1::Range begin_end = _waveforms_hcal->getWaveforms();
  WaveformContainerv1::Iter iwave = begin_end.first;
  int j = 0;
  for (; iwave != begin_end.second; ++iwave)
    {
      std::vector<int> wave = *(iwave->second);
      for (int i = 0; i < m_nsamples; i++)
	{
	  m_waveforms_hcal[iwave->first][i] = static_cast<int>(wave.at(i));
	  if (m_waveforms_hcal[iwave->first][i] > (0x3fff)) m_waveforms_hcal[iwave->first][i] = 0x3fff;
	}

      for (int i = 0; i < m_nsamples - 6;i++)
	{
	  m_peak_sub_ped[iwave->first][i] = m_waveforms_hcal[iwave->first][5 + i] - m_waveforms_hcal[iwave->first][i];
	  if (m_peak_sub_ped[iwave->first][i] < 0) m_peak_sub_ped[iwave->first][i] = 0;
	}
      if (_verbose > 3) 
	{
	  std::cout << "ADC "<< j/64<<" channel "<< j%64 <<" : ";

	  for (int i = 0; i < m_nsamples; i++)
	    {
	      std::cout << m_waveforms_hcal[iwave->first][i] << " ";
	    }
	  if (_verbose > 3) std::cout <<"."<< std::endl;
	  if (_verbose > 3) std::cout << "ADC "<< j/64<<" pedpeak "<< j%64 <<" : ";
	  for (int i = 0; i < m_nsamples; i++)
	    {
	      if (_verbose > 3) std::cout << m_peak_sub_ped[iwave->first][i] << " ";       
	    }	      
	  if (_verbose > 3) std::cout <<"."<< std::endl;
      
	}

      j++;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}
int CaloTriggerEmulator::process_primitives()
{

  if (_verbose) std::cout << __FILE__<<"::"<<__FUNCTION__ <<"::"<<__LINE__<< std::endl;

  unsigned int i, j;
  m_trigger_primitives.clear();
  
  std::vector<unsigned int> trigger_prims;
  for (i = 0; i < 24; i++)
    {
      
      unsigned int tmp;  
      if (_verbose > 3) std::cout << "i: "<<i<<": "<<std::endl;
	  
      for (int isum = 0; isum < 16; isum++)
	{
	  int id_peak = -1;
	  unsigned int peak = 0;
	  for (int is = 0; is < m_nsamples - 6; is++)
	    {
	      m_trig_sums[i][isum] = 0;
	      for (j = 0; j < 4;j++)
		{
		  tmp = m_l1_adc_table[m_peak_sub_ped[64*i + isum*4 + j][is] >> 4];
		      
		  m_trig_sums[i][isum] += (tmp & 0x3ff);
		      
		}
	      m_trig_sums[i][isum] = (m_trig_sums[i][isum] & 0x3ff) >> 2;
	      if (peak < m_trig_sums[i][isum]) 
		{
		  peak = m_trig_sums[i][isum];
		  id_peak = is;
		}
	      trigger_prims.push_back(m_trig_sums[i][isum]);
	    }
	  peak_primitive[i]->Fill(isum, id_peak);
	  avg_primitive[i]->Fill(isum, peak);
	  trigger_fire_map[i]->Fill(peak > 1, isum%4, isum/4);
	  primitives[i]->Fill(isum, peak);
	  if (_verbose > 3) std::cout<< "." << std::endl;
	  m_trigger_primitives.push_back(trigger_prims);
	}
      if (_verbose) std::cout << __FILE__<<"::"<<__FUNCTION__ <<"::"<<__LINE__<< std::endl;      
    }  
  
  _tree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}
void CaloTriggerEmulator::GetNodes(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FUNCTION__ << std::endl;
  _ll1_hcal = findNode::getClass<LL1Outv1>(topNode, "LL1Out_HCALOUT");

  if (!_ll1_hcal) 
    {
      std::cout << "No LL1Out for Calo found... " << std::endl;
      exit(1);
    }

  _waveforms_hcal = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALOUT");

  if (!_waveforms_hcal) 
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

  LL1Outv1 *ll1out = findNode::getClass<LL1Outv1>(ll1Node, "LL1Out_HCALOUT");
  if (!ll1out)
    {
      ll1out = new LL1Outv1();
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, "LL1Out_HCALOUT", "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
}

int CaloTriggerEmulator::End(PHCompositeNode* topNode)
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
