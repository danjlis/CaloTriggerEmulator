#include "MBDEmulatorTreeMaker.h"
#include <vector>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
// G4Cells includes

#include <iostream>

#include <map>

//____________________________________________________________________________..
MBDEmulatorTreeMaker::MBDEmulatorTreeMaker(const std::string &name, const std::string &outfilename, const std::string &nodename):
  SubsysReco(name)
  
{
  _nodename = nodename;
  _foutname = outfilename;  
  _verbosity = 0;
}

//____________________________________________________________________________..
MBDEmulatorTreeMaker::~MBDEmulatorTreeMaker()
{

}

//____________________________________________________________________________..
int MBDEmulatorTreeMaker::Init(PHCompositeNode *topNode)
{

  _f = new TFile( _foutname.c_str(), "RECREATE");

  std::cout << " making a file = " <<  _foutname.c_str() << " , _f = " << _f << std::endl;
  
  _tree_mbd = new TTree("ttree_mbd","a persevering date tree");
  //  _tree_ll1 = new TTree("ttree_ll1","a persevering local date tree");

  _tree_mbd->Branch("clock_mbd", &b_clock_mbd, "clock_mbd/I");
  //  _tree_ll1->Branch("clock_ll1", &b_clock_ll1, "clock_ll1/I");
  _tree_mbd->Branch("event_mbd", &b_event_mbd, "event_mbd/I");
  // _tree_ll1->Branch("event_ll1", &b_event_ll1, "event_ll1/I");

  _tree_mbd->Branch("trigger_bits_mbd",&b_trigger_bits_mbd);
  // _tree_ll1->Branch("trigger_bits_ll1",b_trigger_bits_ll1);
  
  for (int i = 0; i < 8; i++)
    {
      _tree_mbd->Branch(Form("trigger_words_mbd_%d", i),&b_trigger_words_mbd[i]);
      // _tree_ll1->Branch(Form("trigger_words_ll1_%d", i),&b_trigger_words_ll1[i]);
    }

  for (int i = 0; i < 8; i++)
    {
      for (int j = 0; j < 4; j++)
	{
	  _tree_mbd->Branch(Form("trigger_charge_mbd_%d", i),&b_trigger_charge_mbd[j][i]);
	  // _tree_ll1->Branch(Form("trigger_charge_ll1_%d", i),&b_trigger_charge_ll1[j][i]);
	}
    }

  for (int i = 0; i < 4; i++)
    {
      _tree_mbd->Branch(Form("trigger_nhit_mbd_%d", i),&b_trigger_nhit_mbd[i]);
      // _tree_ll1->Branch(Form("trigger_nhit_ll1_%d", i),&b_trigger_nhit_ll1[i]);

    }

  for (int i = 0; i < 4; i++)
    {
      for (int j = 0; j < 4; j++)
	{
	  _tree_mbd->Branch(Form("trigger_time_mbd_%d", i),&b_trigger_time_mbd[j][i]);
	  // _tree_ll1->Branch(Form("trigger_time_ll1_%d", i),&b_trigger_time_ll1[j][i]);
	}
    }

  _i_event = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MBDEmulatorTreeMaker::InitRun(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

void MBDEmulatorTreeMaker::SetVerbosity(int verbo){
  _verbosity = verbo;
  return;
}

void MBDEmulatorTreeMaker::reset_tree_vars()
{
  
  b_clock_mbd = 0;
  b_event_mbd = 0;
  b_trigger_bits_mbd.clear();

  //  b_clock_ll1 = 0;
  //b_event_ll1 = 0;
  //  b_trigger_bits_ll1.clear();

  for (int i = 0; i < 8; i++)
    {
      b_trigger_words_mbd[i].clear();
      //      b_trigger_words_ll1[i].clear();
      for (int j = 0; j < 4; j++)
	{
	  b_trigger_charge_mbd[j][i].clear();
	  //  b_trigger_charge_ll1[j][i].clear();
	}
    }

  for (int i = 0; i < 4; i++)
    {
      b_trigger_nhit_mbd[i].clear();
      //      b_trigger_nhit_ll1[i].clear();
      for (int j = 0; j < 4; j++)
	{
	  b_trigger_time_mbd[j][i].clear();
	  //	  b_trigger_timee_ll1[j][i].clear();
	}
    }

  return;
}

int MBDEmulatorTreeMaker::process_event(PHCompositeNode *topNode)
{
  _i_event++;

  GetNodes(topNode);
  reset_tree_vars();

  b_clock_mbd = _ll1_mbd->get_clock_number();
  b_event_mbd = _ll1_mbd->get_event_number();

  _trigger_primitives = _ll1_mbd->GetTriggerPrimitiveContainer();
  _trigger_bits = _ll1_mbd->GetTriggerBits();
  LL1Outv2::Range word_range = _ll1_mbd->getTriggerWords();
  int i = 0;
  int j;
  for (auto iter = _trigger_bits->begin(); iter < _trigger_bits->end(); i++, ++iter)
    {
      b_trigger_bits_mbd.push_back((*iter));
    } 
  i = 0;
  for (auto iter = word_range.first ; iter != word_range.second; i++, ++iter)
    {
      for (auto witer = (*iter).second->begin(); witer != (*iter).second->end(); ++witer)
	{
	  b_trigger_words_mbd[i].push_back(*witer);
	}
    }

  i = 0;
  TriggerPrimitiveContainerv1::Range range = _trigger_primitives->getTriggerPrimitives();
  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter, i++)
    {
      _trigger_primitive = (*iter).second;
      TriggerPrimitive::Range srange = _trigger_primitive->getSums();
      j = 0;
      for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter, j++)
	{
	  std::vector<unsigned int> *sum = (*siter).second;
	  for (auto viter = sum->begin(); viter != sum->end(); ++viter)
	    {
	      if (j < 8)
		{
		  b_trigger_charge_mbd[i][j].push_back((*viter));
		}
	      else if (j ==8)
		{
		  b_trigger_nhit_mbd[i].push_back((*viter));
		}
	      else if ( j < 13 )
		{
		  b_trigger_time_mbd[i][j - 9].push_back((*viter));
		}
	      else break;
	     
	    }
	}
    }

  _tree_mbd->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}



void MBDEmulatorTreeMaker::GetNodes(PHCompositeNode* topNode)
{

  _ll1_mbd = findNode::getClass<LL1Outv2>(topNode, _nodename);

  if (!_ll1_mbd) 
    {
      std::cout << "No LL1Out MBD node... " << _nodename << std::endl;
      exit(1);
    }

}

int MBDEmulatorTreeMaker::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "MBDEmulatorTreeMaker::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MBDEmulatorTreeMaker::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "MBDEmulatorTreeMaker::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MBDEmulatorTreeMaker::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "MBDEmulatorTreeMaker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  std::cout<<"Total events: "<<_i_event<<std::endl;
  _f->Write();
  _f->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MBDEmulatorTreeMaker::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "MBDEmulatorTreeMaker::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

