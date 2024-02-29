#include "CaloWaveFormToy.h"
#include "WaveformContainerv1.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>


#include <TFile.h>
#include <TNtuple.h>
#include <TMath.h>
#include <cassert>
#include <sstream>
#include <string>
#include <TF1.h>


using namespace std;

CaloWaveFormToy::CaloWaveFormToy(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , detector("ALL")
{
  _type = 0;
  _occupancy = 1.;
}

CaloWaveFormToy::~CaloWaveFormToy()
{
}

int CaloWaveFormToy::Init(PHCompositeNode*)
{
  rnd = new TRandom3(0);

  //----------------------------------------------------------------------------------------------------
  //Read in the template file, this currently points to a tim local area file, 
  //but a copy of this file is in the git repository.
  //----------------------------------------------------------------------------------------------------


  for (int i = 0 ; i < 24576;i++)
    {
      for (int j = 0; j < _nsamples; j++) m_waveform_emcal[i].push_back(0.);
    }

  for (int i = 0 ; i < 1536;i++)
    {
      for (int j = 0; j < _nsamples; j++) m_waveform_ihcal[i].push_back(0.);
      for (int j = 0; j < _nsamples; j++) m_waveform_ohcal[i].push_back(0.);
    }

  for (int i = 0 ; i < 256;i++)
    {
      for (int j = 0; j < _nsamples; j++) m_waveform_mbd[i].push_back(0.);
    }


  return 0;
}


int CaloWaveFormToy::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity()) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;
  CreateNodes(topNode);
  return 0;
}


void CaloWaveFormToy::CreateNodes(PHCompositeNode* topNode)
{

  if (Verbosity()) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    }
  PHNodeIterator dstIter(dstNode);

  // Create nodes for EMCAL
  if (IsDetector("MBD"))
    {
      PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "MBD"));
      if (!detNode)
	{
	  std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
	  detNode = new PHCompositeNode("EMCAL");
	  dstNode->addNode(detNode);
	}
      
      WaveformContainerv1 *waveforms = findNode::getClass<WaveformContainerv1>(detNode, "WAVEFORMS_MBD");
      if (!waveforms)
	{
	  waveforms = new WaveformContainerv1();
	  PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(waveforms, "WAVEFORMS_MBD", "PHObject");
	  detNode->addNode(waveformcontainerNode);
	}
    }


  // Create nodes for EMCAL
  if (IsDetector("EMCAL"))
    {
      PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "EMCAL"));
      if (!detNode)
	{
	  std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
	  detNode = new PHCompositeNode("EMCAL");
	  dstNode->addNode(detNode);
	}
      
      WaveformContainerv1 *waveforms = findNode::getClass<WaveformContainerv1>(detNode, "WAVEFORMS_EMCAL");
      if (!waveforms)
	{
	  waveforms = new WaveformContainerv1();
	  PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(waveforms, "WAVEFORMS_EMCAL", "PHObject");
	  detNode->addNode(waveformcontainerNode);
	}
    }

  // Create nodes for HCALIN
  if (IsDetector("HCALIN"))
    {
      PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "HCALIN"));
      if (!detNode)
	{
	  std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
	  detNode = new PHCompositeNode("HCALIN");
	  dstNode->addNode(detNode);
	}

      WaveformContainerv1 *waveforms = findNode::getClass<WaveformContainerv1>(detNode, "WAVEFORMS_HCALIN");
      if (!waveforms)
	{
	  waveforms = new WaveformContainerv1();
	  PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(waveforms, "WAVEFORMS_HCALIN", "PHObject");
	  detNode->addNode(waveformcontainerNode);
	}
    }
  if (IsDetector("HCALOUT"))
    {
      PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "HCALOUT"));
      if (!detNode)
	{
	  std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
	  detNode = new PHCompositeNode("HCALOUT");
	  dstNode->addNode(detNode);
	}

      WaveformContainerv1 *waveforms = findNode::getClass<WaveformContainerv1>(detNode, "WAVEFORMS_HCALOUT");
      if (!waveforms)
	{
	  waveforms = new WaveformContainerv1();
	  PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(waveforms, "WAVEFORMS_HCALOUT", "PHObject");
	  detNode->addNode(waveformcontainerNode);
	}

    }
  return;
}


int CaloWaveFormToy::process_event(PHCompositeNode* topNode)
{
  
  if (Verbosity()) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  if (IsDetector("EMCAL"))
    {  
      if (Verbosity()) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      waveforms_emcal = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_EMCAL");
      if (!waveforms_emcal)
	{
	  std::cout << "Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}

    }
  if (IsDetector("HCALIN"))
    {
      if (Verbosity()) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      waveforms_hcalin = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALIN");
      if (!waveforms_hcalin)
	{
	  std::cout << "Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}
      
    }
  if (IsDetector("HCALOUT"))
    {
      if (Verbosity()) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      waveforms_hcalout = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALOUT");
      if (!waveforms_hcalout)
	{
	  std::cout << "Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}      
      
    }

  if (IsDetector("MBD"))
    {
      if (Verbosity()) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      waveforms_mbd = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_MBD");
      if (!waveforms_mbd)
	{
	  std::cout << "Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}      
    }

  for (int i = 0 ; i < 24576;i++)
    {
      m_waveform_emcal[i].clear();
      for (int j = 0; j < _nsamples; j++) m_waveform_emcal[i].push_back(0.);
    }

  for (int i = 0 ; i < 1536;i++)
    {
      m_waveform_ihcal[i].clear();
      m_waveform_ohcal[i].clear();
      for (int j = 0; j < _nsamples; j++) m_waveform_ihcal[i].push_back(0.);
      for (int j = 0; j < _nsamples; j++) m_waveform_ohcal[i].push_back(0.);
    }

  for (int i = 0 ; i < 256;i++)
    {
      m_waveform_mbd[i].clear();
      for (int j = 0; j < _nsamples; j++) m_waveform_mbd[i].push_back(0.);
    }

  
  //---------------------------------------------------------
  // do EMCAL
  //---------------------------------------------------------


  if (IsDetector("EMCAL"))
    {
      
      for (int i = 0 ; i < 24576; i++)
	{
	  
	  for (int j = 0; j < _nsamples;j++) // 16 is the number of time samples
	    {
	      m_waveform_emcal[i][j] = (j >= 4 && j < 8 ? (0x2f2e & 0x3fff) : 0);
	    }
	}
    }

  //--------------
  // do IHCAL
  //--------------
  if (IsDetector("HCALIN"))
    {
      for (int i = 0 ; i < 1536; i++)
	{
	  for (int j = 0; j < _nsamples;j++) // 16 is the number of time samples
	    {
	      m_waveform_ihcal[i][j] = (j >= 4 && j < 8 ? (0x2f2e & 0x3fff) : 0);
	    }
	}
    }
  //--------------
  // do OHCAL
  //--------------

  if (IsDetector("HCALOUT"))
    {
      for (int i = 0 ; i < 1536; i++)
	{
	  for (int j = 0; j < _nsamples;j++) // 16 is the number of time samples
	    {
	      m_waveform_ohcal[i][j] = (j >= 4 && j < 8 ? (0x2f2e & 0x3fff) : 0);
	    }
	}
    }

  
  if (IsDetector("MBD"))
    {
      for (int i = 0 ; i < 256; i++)
	{
	  for (int j = 0; j < _nsamples;j++) // 16 is the number of time samples
	    {
	      m_waveform_mbd[i][j] = (j >= 4 && j < 8 ? (0x2f2e & 0x3fff) : 0);
	    }
	}
    }  

  if (IsDetector("EMCAL"))
    {
      for (int i = 0; i < 24576;i++)
	{
	  std::vector<int> *wave = new vector<int>;
	  wave->clear();
      
	  for (int k = 0; k < _nsamples;k++)
	    {
	      wave->push_back(static_cast<int>(m_waveform_emcal[i][k]));
	    }
	  waveforms_emcal->set_waveform_at_channel(i, wave);
	}
    }
  //---------------------------
  // do noise for ihcal:
  //---------------------------
  if (IsDetector("HCALIN"))
    {  
      for (int i = 0; i < 1536;i++)
	{
	  std::vector<int> *wave = new vector<int>;
	  wave->clear();
	  for (int k = 0; k < 16;k++)
	    {
	      wave->push_back(static_cast<int>(m_waveform_ihcal[i][k]));
	    }
	  waveforms_hcalin->set_waveform_at_channel(i, wave);
	}
    }
  //---------------------------
  // do noise for ohcal:
  //---------------------------
  if (IsDetector("HCALOUT"))
    {
      for (int i = 0; i < 1536;i++)
	{
	  std::vector<int> *wave = new vector<int>;
	  wave->clear();      
      
	  for (int k = 0; k < 16;k++)
	    {
	      wave->push_back(static_cast<int>(m_waveform_ohcal[i][k]));
	    }
	  waveforms_hcalout->set_waveform_at_channel(i, wave);
	}
    }

  if (IsDetector("MBD"))
    {
      for (int i = 0; i < 256;i++)
	{
	  std::vector<int> *wave = new vector<int>;
	  wave->clear();
      
	  for (int k = 0; k < 16;k++)
	    {
	      wave->push_back(static_cast<int>(m_waveform_mbd[i][k]));
	    }
	  waveforms_mbd->set_waveform_at_channel(i, wave);
	}
    }


  return Fun4AllReturnCodes::EVENT_OK;
}



int CaloWaveFormToy::End(PHCompositeNode* topNode)
{
  return 0;
}

bool CaloWaveFormToy::IsDetector(const std::string &det)
{
  if (strcmp("ALL", detector.c_str()) == 0) return true;
  if (strcmp(det.c_str(), detector.c_str()) == 0) return true;
  return false;
}
