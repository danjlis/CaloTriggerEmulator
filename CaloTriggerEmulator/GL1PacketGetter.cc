#include "GL1PacketGetter.h"

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
GL1PacketGetter::GL1PacketGetter(const std::string &name, const std::string &outfilename)
  : SubsysReco(name)
  , m_packet_low(INT_MIN)
  , m_packet_high(INT_MIN)
{

  _foutname = outfilename;  
  _verbose = 0;
}

//____________________________________________________________________________..
GL1PacketGetter::~GL1PacketGetter()
{

}

//____________________________________________________________________________..
int GL1PacketGetter::Init(PHCompositeNode *topNode)
{

  _f = new TFile( _foutname.c_str(), "RECREATE");

  std::cout << " making a file = " <<  _foutname.c_str() << " , _f = " << _f << std::endl;

  _t = new TTree( "ttree","a persevering date tree");
  
  _t->Branch("packet", &b_packet);
  _t->Branch("bco", &b_bco);
  for (int i = 0; i < 64; i++)
    {
      _t->Branch(Form("raw_%d", i), &b_raw[i]);
      _t->Branch(Form("live_%d", i), &b_live[i]);
      _t->Branch(Form("scaled_%d", i), &b_scaled[i]);
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int GL1PacketGetter::InitRun(PHCompositeNode *topNode)
{

  m_packet_low = 14001;
  m_packet_high = 14001;
  m_ntriggers = 64;
  topNode->print();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int GL1PacketGetter::process_event(PHCompositeNode *topNode)
{

  b_packet = 0;
  b_bco = 0;


  for (int i = 0 ; i < m_ntriggers;i++)
    {

      b_raw[i] = 0;
      b_live[i] = 0;
      b_scaled[i] = 0;
    }

  Event *_event = findNode::getClass<Event>(topNode, "PRDF");
  if (_event == nullptr)
    {
      std::cout << "GL1UnpackPRDF::Process_Event - Event not found" << std::endl;
      return -1;
    }
  for (int pid = m_packet_low; pid <= m_packet_high; pid++)
    {
      Packet *packet = _event->getPacket(pid);
      if (packet)
	{
	  b_packet = packet->iValue(0);
	  b_bco = packet->lValue(0, "BCO");

	  for (int i = 0 ; i < m_ntriggers;i++)
	    {

	      b_raw[i] = packet->lValue(i, 0);
	      b_live[i] = packet->lValue(i, 1);
	      b_scaled[i] = packet->lValue(i, 2);

	    }
	  
	  delete packet;
	}
      else // if the packet is missing treat constitutent channels as zero suppressed 
	{
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
    }

  _t->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}


int GL1PacketGetter::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "GL1PacketGetter::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }

  _f->Write();
  _f->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}
