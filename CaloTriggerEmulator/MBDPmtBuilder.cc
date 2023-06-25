#include "MBDPmtBuilder.h"

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
MBDPmtBuilder::MBDPmtBuilder(const std::string &name)
  : SubsysReco(name)
  , m_dettype(MBDPmtBuilder::CEMC)
  , m_WaveformContainer(nullptr)
  , m_packet_low(INT_MIN)
  , m_packet_high(INT_MIN)
  , m_nsamples(16)
  , m_nchannels(128)
  , m_isdata(true)
{

}

//____________________________________________________________________________..
MBDPmtBuilder::~MBDPmtBuilder()
{

}

//____________________________________________________________________________..
int MBDPmtBuilder::InitRun(PHCompositeNode *topNode)
{


  m_packet_low = 1001;
  m_packet_high = 1002;
  m_nchannels = 128;
  
 
  CreateNodeTree(topNode);
  topNode->print();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MBDPmtBuilder::process_event(PHCompositeNode *topNode)
{
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
	      std::vector<int> waveform;
	      waveform.reserve(m_nsamples);
	      for (int samp = 0; samp < m_nsamples; samp++)
		{
		  waveform.push_back(packet->iValue(samp, channel));
		}
	      waveforms.push_back(waveform);
	      waveform.clear();
	    }
	  delete packet;
	}
      else // if the packet is missing treat constitutent channels as zero suppressed 
	{
	  return Fun4AllReturnCodes::DISCARDEVENT;
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
      wave->reserve(m_nsamples);
      for (int k = 0; k < m_nsamples;k++)
	{
	  wave->push_back(static_cast<int>(waveforms[iwave][k]));
	}

      m_WaveformContainer->set_waveform_at_channel(iwave, wave);
      iwave++;
    }

  waveforms.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

void MBDPmtBuilder::CreateNodeTree(PHCompositeNode *topNode)
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

  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "BBC"));
  if (!detNode)
    {
      std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
      detNode = new PHCompositeNode("BBC");
      dst_node->addNode(detNode);
    }

  WaveformContainerv1 *waveforms = findNode::getClass<WaveformContainerv1>(detNode, "WAVEFORMS_BBC");
  if (!waveforms)
    {
      waveforms = new WaveformContainerv1();
      PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(waveforms, "WAVEFORMS_BBC", "PHObject");
      detNode->addNode(waveformcontainerNode);
    }

}
