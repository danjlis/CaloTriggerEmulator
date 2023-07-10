
#include "CaloPacketGetter.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/packet.h>

#include <iostream>  // for operator<<, endl, basic...
#include <memory>    // for allocator_traits<>::val...
#include <vector>    // for vector

//____________________________________________________________________________..
CaloPacketGetter::CaloPacketGetter(const std::string &name, const std::string &detector)
  : SubsysReco(name)
  , m_WaveformContainer(nullptr)
  , m_packet_low(INT_MIN)
  , m_packet_high(INT_MIN)
  , m_nsamples(16)
  , m_nchannels(192)
  , m_isdata(true)
{
  m_detectors["CEMC"] = DetectorSystem::CEMC;
  m_detectors["HCALIN"] = DetectorSystem::HCALIN;
  m_detectors["HCALOUT"] = DetectorSystem::HCALOUT;
  m_detectors["SEPD"] = DetectorSystem::SEPD;

  m_detector_type = m_detectors[detector];
}

//____________________________________________________________________________..
CaloPacketGetter::~CaloPacketGetter()
{

}

//____________________________________________________________________________..
int CaloPacketGetter::InitRun(PHCompositeNode *topNode)
{

  std::cout << m_detector_type <<std::endl;
  switch (m_detector_type) {
    case DetectorSystem::CEMC:
      std::cout <<"CEMC Packet Getter"<<std::endl;
      m_packet_low = 6001;
      m_packet_high = 6128;
      m_nchannels = 192;

      break;
    case DetectorSystem::HCALIN:
      std::cout <<"HCALIN Packet Getter"<<std::endl;
      m_packet_low = 7001;
      m_packet_high = 7008;
      m_nchannels = 192;

      break;
    case DetectorSystem::HCALOUT:
      std::cout <<"HCALOUT Packet Getter"<<std::endl;
      m_packet_low = 8001;
      m_packet_high = 8008;
      m_nchannels = 192;

      break;
    case DetectorSystem::SEPD:
      std::cout <<"SEPD Packet Getter"<<std::endl;
      m_packet_low = 9001;
      m_packet_high = 9006;
      m_nchannels = 256;

      break;

    }
  
  CreateNodeTree(topNode);
  topNode->print();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloPacketGetter::process_event(PHCompositeNode *topNode)
{
  //  std::cout << m_detector_type <<std::endl;
  switch(m_detector_type){
  case DetectorSystem::CEMC:
    {
      
      m_WaveformContainer = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_CEMC");
      if (!m_WaveformContainer)
	{
	  std::cout << "CEMC Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}
      break;
    }
  case DetectorSystem::HCALIN:
    {
      
      m_WaveformContainer = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALIN");
      if (!m_WaveformContainer)
	{
	  std::cout << "HCALIN Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}
      break;
    }
  case DetectorSystem::HCALOUT:
    {
      
      m_WaveformContainer = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALOUT");
      if (!m_WaveformContainer)
	{
	  std::cout << "HCALOUT Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}
      break;
    }

  case DetectorSystem::SEPD:
    {
      
      m_WaveformContainer = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_SEPD");
      if (!m_WaveformContainer)
	{
	  std::cout << "Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}
      break;
    }
  }
  unsigned int key;
  std::vector<unsigned int> keys;
  std::vector<std::vector<int>> waveforms;
  if (m_isdata)
  {
    Event *_event = findNode::getClass<Event>(topNode, "PRDF");
    if (_event == nullptr)
    {
      std::cout << "CaloUnpackPRDF::Process_Event - Event not found" << std::endl;
      return -1;
    }
    if (_event->getEvtType() >= 8)  /// special event where we do not read out the calorimeters
    {
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
    for (int pid = m_packet_low; pid <= m_packet_high; pid++)
    {
      Packet *packet = _event->getPacket(pid);
      if (packet)
	{

	  int nchannels = packet->iValue(0, "CHANNELS");
	  if (nchannels > m_nchannels) // packet is corrupted and reports too many channels
	    {
	      return Fun4AllReturnCodes::DISCARDEVENT;
	    }
	  for (int channel = 0; channel < nchannels; channel++)
	    {
	      //	      std::cout << "Channel "<<channel<< " : ";
	      std::vector<int> waveform;
	      key = 0;
	      key |= (pid << 16);
	      key |= channel;
	      keys.push_back(key);
	      waveform.reserve(m_nsamples);
	      for (int samp = 0; samp < m_nsamples; samp++)
		{

		  waveform.push_back(packet->iValue(samp, channel));
		  //  std::cout << " "<< waveform[samp];
		}
	      //	      std::cout << " "<<std::endl;
	      waveforms.push_back(waveform);
	      waveform.clear();
	    }
	  delete packet;
	}
      else // if the packet is missing treat constitutent channels as zero suppressed 
	{
	  for (int channel = 0; channel < m_nchannels; channel++)
	    {
	      std::vector<int> waveform;
	      key = 0;
	      key |= (pid << 16);
	      key |= channel;
	      keys.push_back(key);
	      waveform.reserve(m_nsamples);
	      for (int samp = 0; samp < m_nsamples; samp++)
		{

		  waveform.push_back(1500);

		}
	      waveforms.push_back(waveform);
	      waveform.clear();
	    }

	}
    }
  }
  else  // placeholder for adding simulation
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  
  int iwave = 0;
  for (auto it = waveforms.begin(); it != waveforms.end(); it++)
    {
      std::vector<int> *wave = new std::vector<int>();
      wave->clear();

      //      std::cout << "After "<<iwave<< " : ";
	    
         
      for (int k = 0; k < m_nsamples;k++)
	{
	  wave->push_back(waveforms[iwave][k]);

	  //  std::cout << " "<< wave->at(k);
	}

      //      std::cout << " "<<std::endl;
	
      m_WaveformContainer->set_waveform_at_key(keys[iwave], wave);
      iwave++;
    }

  waveforms.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloPacketGetter::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator nodeItr(topNode);
  // DST node
  PHCompositeNode *dst_node = dynamic_cast<PHCompositeNode *>(
      nodeItr.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    std::cout << "PHComposite node created: DST" << std::endl;
    dst_node = new PHCompositeNode("DST");
    topNode->addNode(dst_node);
  }
  // towers
  PHNodeIterator dstIter(dst_node);
  switch(m_detector_type)
    {
      
    case DetectorSystem::CEMC:
      {
	PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "CEMC"));
	if (!detNode)
	  {
	    std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
	    detNode = new PHCompositeNode("CEMC");
	    dst_node->addNode(detNode);
	  }
	
	WaveformContainerv1 *waveforms = findNode::getClass<WaveformContainerv1>(detNode, "WAVEFORMS_CEMC");
	if (!waveforms)
	  {
	    waveforms = new WaveformContainerv1();
	    PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(waveforms, "WAVEFORMS_CEMC", "PHObject");
	    detNode->addNode(waveformcontainerNode);
	  }
	
	break;
      }
    case DetectorSystem::HCALIN:
      {
	PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "HCALIN"));
	if (!detNode)
	  {
	    std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
	    detNode = new PHCompositeNode("HCALIN");
	    dst_node->addNode(detNode);
	  }
	
	WaveformContainerv1 *waveforms = findNode::getClass<WaveformContainerv1>(detNode, "WAVEFORMS_HCALIN");
	if (!waveforms)
	  {
	    waveforms = new WaveformContainerv1();
	    PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(waveforms, "WAVEFORMS_HCALIN", "PHObject");
	    detNode->addNode(waveformcontainerNode);
	  }
	
	break;
      }
    case DetectorSystem::HCALOUT:
      {
	PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "HCALOUT"));
	if (!detNode)
	  {
	    std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
	    detNode = new PHCompositeNode("HCALOUT");
	    dst_node->addNode(detNode);
	  }
	
	WaveformContainerv1 *waveforms = findNode::getClass<WaveformContainerv1>(detNode, "WAVEFORMS_HCALOUT");
	if (!waveforms)
	  {
	    waveforms = new WaveformContainerv1();
	    PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(waveforms, "WAVEFORMS_HCALOUT", "PHObject");
	    detNode->addNode(waveformcontainerNode);
	  }
	
	break;
      }
    case DetectorSystem::SEPD:
      {
	PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SEPD"));
	if (!detNode)
	  {
	    std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
	    detNode = new PHCompositeNode("SEPD");
	    dst_node->addNode(detNode);
	  }
	
	WaveformContainerv1 *waveforms = findNode::getClass<WaveformContainerv1>(detNode, "WAVEFORMS_SEPD");
	if (!waveforms)
	  {
	    waveforms = new WaveformContainerv1();
	    PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(waveforms, "WAVEFORMS_SEPD", "PHObject");
	    detNode->addNode(waveformcontainerNode);
	  }
	
	break;
      }
    }
}
