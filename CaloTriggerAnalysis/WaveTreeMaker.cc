#include "WaveTreeMaker.h"
#include <vector>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
// G4Cells includes
#include <calowaveformsim/ADCDefs.h>
#include <TLorentzVector.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos


#include <iostream>

#include <map>

//____________________________________________________________________________..
WaveTreeMaker::WaveTreeMaker(const std::string &name, const std::string &outfilename):
  SubsysReco(name)
{
  _foutname = outfilename;  
  _verbosity = 0;
}

//____________________________________________________________________________..
WaveTreeMaker::~WaveTreeMaker()
{

}

//____________________________________________________________________________..
int WaveTreeMaker::Init(PHCompositeNode *topNode)
{

  _f = new TFile( _foutname.c_str(), "RECREATE");

  std::cout << " making a file = " <<  _foutname.c_str() << " , _f = " << _f << std::endl;
  
  _tree = new TTree("ttree","a persevering date tree");

  _tree->Branch("evtnr", &m_evtnr,"evtnr/i");
  _tree->Branch("clk", &m_clk,"clk/i");
  _tree->Branch("femevtnr", &m_femevtnr,"femevtnr/i");
  _tree->Branch("femclk", &m_femclk,"femclk/i");
  _tree->Branch("waveform_emcal",&m_waveforms_cemc,"waveform_emcal[24576][31]/I");
  _tree->Branch("waveform_ihcal",&m_waveforms_hcalin,"waveform_ihcal[1536][31]/I");
  _tree->Branch("waveform_ohcal",&m_waveforms_hcalout,"waveform_ohcal[1536][31]/I");
  _tree->Branch("waveform_bbc",&m_waveforms_bbc,"waveform_bbc[256][31]/I");
  _tree->Branch("waveform_zdc",&m_waveforms_zdc,"waveform_zdc[16][31]/I");
  _i_event = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int WaveTreeMaker::InitRun(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

void WaveTreeMaker::SetVerbosity(int verbo){
  _verbosity = verbo;
  return;
}

void WaveTreeMaker::reset_tree_vars()
{
  m_evtnr = 0;
  m_clk = 0;
  m_femevtnr = 0;
  m_femclk = 0;


  for (int i = 0; i < 24576; i++)
    {

      for (int j = 0; j < 31; j++) m_waveforms_cemc[i][j] = 0;
    }
  for (int i = 0; i < 1536; i++)
    {
      for (int j = 0; j < 31; j++) {
	m_waveforms_hcalout[i][j] = 0;
	m_waveforms_hcalin[i][j] = 0;
      }
    }
  for (int i = 0; i < 128; i++)
    {
      for (int j = 0; j < 31; j++) m_waveforms_bbc[i][j] = 0;
    }
  for (int i = 0; i < 16; i++)
    {
      for (int j = 0; j < 31; j++) m_waveforms_zdc[i][j] = 0;
    }


  
  return;
}
void WaveTreeMaker::GetNodes(PHCompositeNode* topNode)
{

  if (_verbosity) std::cout << __FUNCTION__ << std::endl;

  _waveforms_cemc = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_CEMC");


  _waveforms_hcalout = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALOUT");


  _waveforms_hcalin = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALIN");


  _waveforms_bbc = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_MBD");

  _waveforms_zdc = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_ZDC");


  return;

}


void WaveTreeMaker::process_waveforms()
{
  if (_waveforms_cemc)
    {
      WaveformContainerv1::Range begin_end = _waveforms_cemc->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;
      for (; iwave != begin_end.second; ++iwave)
	{
	  std::vector<int> wave = *(iwave->second);
	  for (int i = 0; i < 31; i++)
	    {
	      m_waveforms_cemc[iwave->first][i] = wave.at(i);
	    }
	}
    }
  if (_waveforms_hcalout)
    {
      WaveformContainerv1::Range begin_end = _waveforms_hcalout->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;
      for (; iwave != begin_end.second; ++iwave)
	{
	  std::vector<int> wave = *(iwave->second);
	  for (int i = 0; i < 31; i++)
	    {
	      m_waveforms_hcalout[iwave->first][i] = wave.at(i);
	    }
	}
    }
  if (_waveforms_hcalin)
    {
      WaveformContainerv1::Range begin_end = _waveforms_hcalin->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;
      for (; iwave != begin_end.second; ++iwave)
	{
	  std::vector<int> wave = *(iwave->second);
	  for (int i = 0; i < 31; i++)
	    {
	      m_waveforms_hcalin[iwave->first][i] = wave.at(i);
	    }
	}
    }
  if (_waveforms_bbc)
    {
      WaveformContainerv1::Range begin_end = _waveforms_bbc->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;
      int ich = 0;
      for (; iwave != begin_end.second; ++iwave)
	{
	 
	  std::vector<int> wave = *(iwave->second);
	  
	  for (int i = 0; i < 31; i++)
	    {
	      m_waveforms_bbc[ich][i] = wave.at(i);
	    }
	  ich++;
	}
    }

  if (_waveforms_zdc)
    {

      int ich = 0;

      WaveformContainerv1::RangeInfo begin_end_info;
      WaveformContainerv1::IterInfo iinfo;

      begin_end_info = _waveforms_zdc->get_packet_events();
      iinfo = begin_end_info.first;
      m_evtnr = iinfo->second;
      begin_end_info = _waveforms_zdc->get_packet_clocks();
      iinfo = begin_end_info.first;
      m_clk = iinfo->second;
      begin_end_info = _waveforms_zdc->get_fem_events();
      iinfo = begin_end_info.first;
      m_femevtnr = iinfo->second;
      begin_end_info = _waveforms_zdc->get_fem_clocks();
      iinfo = begin_end_info.first;
      m_femclk = iinfo->second;


      WaveformContainerv1::Range begin_end = _waveforms_zdc->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;


      for (; iwave != begin_end.second; ++iwave)
	{

	  std::vector<int> wave = *(iwave->second);
	  
	  for (int i = 0; i < 31; i++)
	    {
	      m_waveforms_zdc[ich][i] = wave.at(i);
	    }
	  ich++;
	}
    }

}

int WaveTreeMaker::process_event(PHCompositeNode *topNode)
{
  _i_event++;
  GetNodes(topNode);
  reset_tree_vars();
  process_waveforms();

  _tree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}


//int TrigTreeMaker::CreateNode(PHCompositeNode *topNode)
//{

//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
int WaveTreeMaker::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "WaveTreeMaker::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int WaveTreeMaker::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "WaveTreeMaker::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int WaveTreeMaker::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "WaveTreeMaker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  std::cout<<"Total events: "<<_i_event<<std::endl;
  _f->Write();
  _f->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int WaveTreeMaker::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "WaveTreeMaker::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void WaveTreeMaker::Print(const std::string &what) const
{
  std::cout << "WaveTreeMaker::Print(const std::string &what) const Printing info for " << what << std::endl;
}
