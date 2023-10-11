#include "HCALEmulatorTreeMaker.h"
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
HCALEmulatorTreeMaker::HCALEmulatorTreeMaker(const std::string &name, const std::string &outfilename, const std::string &nodename):
  SubsysReco(name)
  
{
  useCaloTowerBuilder = false;
  _nodename = nodename;
  _foutname = outfilename;  
  _verbosity = 0;
}

//____________________________________________________________________________..
HCALEmulatorTreeMaker::~HCALEmulatorTreeMaker()
{

}

//____________________________________________________________________________..
int HCALEmulatorTreeMaker::Init(PHCompositeNode *topNode)
{

  _f = new TFile( _foutname.c_str(), "RECREATE");

  std::cout << " making a file = " <<  _foutname.c_str() << " , _f = " << _f << std::endl;
  
  _tree_hcal = new TTree("ttree_hcal","a persevering date tree");
  //  _tree_ll1 = new TTree("ttree_ll1","a persevering local date tree");

  _tree_hcal->Branch("clock_hcalout", &b_clock_hcalout, "clock_hcalout/I");
  _tree_hcal->Branch("event_hcalout", &b_event_hcalout, "event_hcalout/I");

  _tree_hcal->Branch("clock_hcalin", &b_clock_hcalin, "clock_hcalin/I");
  _tree_hcal->Branch("event_hcalin", &b_event_hcalin, "event_hcalin/I");

  _tree_hcal->Branch("trigger_bits_hcal",&b_trigger_bits_hcal);
  
  for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < 16; j++)
	{
	  _tree_hcal->Branch(Form("trigger_sum_hcalin_%d", i),&b_trigger_sum_hcalin[i][j]);
	  _tree_hcal->Branch(Form("trigger_sum_hcalout_%d", i),&b_trigger_sum_hcalout[i][j]);
	}
    }

  _i_event = 0;

  if (useCaloTowerBuilder)
    {
      _tree_hcal->Branch("hcalin_energy",&b_hcalin_energy);
      _tree_hcal->Branch("hcalin_phibin",&b_hcalin_phibin);
      _tree_hcal->Branch("hcalin_etabin",&b_hcalin_etabin);
      _tree_hcal->Branch("hcalout_energy",&b_hcalout_energy);
      _tree_hcal->Branch("hcalout_phibin",&b_hcalout_phibin);
      _tree_hcal->Branch("hcalout_etabin",&b_hcalout_etabin);

    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HCALEmulatorTreeMaker::InitRun(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

void HCALEmulatorTreeMaker::SetVerbosity(int verbo){
  _verbosity = verbo;
  return;
}

void HCALEmulatorTreeMaker::reset_tree_vars()
{
  
  b_clock_hcalin = 0;
  b_event_hcalin = 0;
  b_clock_hcalout = 0;
  b_event_hcalout = 0;

  b_trigger_bits_hcal.clear();
  b_hcalin_energy.clear();
  b_hcalin_phibin.clear();
  b_hcalin_etabin.clear();
  b_hcalout_energy.clear();
  b_hcalout_phibin.clear();
  b_hcalout_etabin.clear();

  for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < 16; j++)
	{
	  b_trigger_sum_hcalin[i][j].clear();
	  b_trigger_sum_hcalout[i][j].clear();
	}
    }

  return;
}

int HCALEmulatorTreeMaker::process_event(PHCompositeNode *topNode)
{

  if (_verbosity) std::cout << __FILE__ << " "<< __LINE__<<" "<<std::endl;
  int i;
  int j;

  _i_event++;


  reset_tree_vars();

  _ll1_hcal = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_HCALIN");

  if (!_ll1_hcal) 
    {
      std::cout << "No LL1Out HCAL node... " << _nodename << std::endl;
      exit(1);
    }

  if (_verbosity) std::cout << __FILE__ << " "<< __LINE__<<" "<<std::endl;
  b_clock_hcalin = _ll1_hcal->get_clock_number();
  b_event_hcalin = _ll1_hcal->get_event_number();

  _trigger_primitives = _ll1_hcal->GetTriggerPrimitiveContainer();

  if (_verbosity) std::cout << __FILE__ << " "<< __LINE__<<" "<<std::endl;
  if (!_trigger_primitives)
    {
      cout << "no hcalin primitives" <<endl;
      exit(1);
    }
  i = 0;
  {
  TriggerPrimitiveContainerv1::Range range = _trigger_primitives->getTriggerPrimitives();
  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter, i++)
    {
      if (i >= 24)
	{
	  std::cout << "__LINE__ out of bounds: prim "<<i<<" sum"<<j<<endl;
	  exit(1);
	}
      
      _trigger_primitive = (*iter).second;
      TriggerPrimitive::Range srange = _trigger_primitive->getSums();
      j = 0;
      for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter, j++)
	{
	  if (j >= 16)
	    {
	      std::cout << "__LINE__ out of bounds: prim "<<i<<" sum"<<j<<endl;
	    }
	  std::vector<unsigned int> *sum = (*siter).second;
	  for (auto viter = sum->begin(); viter != sum->end(); ++viter)
	    {
		b_trigger_sum_hcalin[i][j].push_back((*viter));
	    }
	}
    }
  }
  if (_verbosity) std::cout << __FILE__ << " "<< __LINE__<<" "<<std::endl;
  _ll1_hcal = findNode::getClass<LL1Outv2>(topNode, "LL1OUT_HCALOUT");

  if (!_ll1_hcal) 
    {
      std::cout << "No LL1Out HCAL node... " << _nodename << std::endl;
      exit(1);
    }

  b_clock_hcalout = _ll1_hcal->get_clock_number();
  b_event_hcalout = _ll1_hcal->get_event_number();

  _trigger_primitives = _ll1_hcal->GetTriggerPrimitiveContainer();

  if (!_trigger_primitives)
    {
      cout << "no hcalout primitives" <<endl;
      exit(1);
    }
  i = 0;
  {
  TriggerPrimitiveContainerv1::Range range = _trigger_primitives->getTriggerPrimitives();
  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter, i++)
    {
      if (i >= 24)
	{
	  std::cout << "__LINE__ out of bounds: prim "<<i<<" sum"<<j<<endl;
	  exit(1);
	}

      _trigger_primitive = (*iter).second;

      TriggerPrimitive::Range srange = _trigger_primitive->getSums();
      j = 0;
      for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter, j++)
	{
	  std::vector<unsigned int> *sum = (*siter).second;
	  for (auto viter = sum->begin(); viter != sum->end(); ++viter)
	    {
		b_trigger_sum_hcalout[i][j].push_back((*viter));
	    }
	}
    }

  }
  if (_verbosity) std::cout << __FILE__ << " "<< __LINE__<<" "<<std::endl;
  _ll1_hcal = findNode::getClass<LL1Outv2>(topNode, _nodename);

  if (!_ll1_hcal) 
    {
      std::cout << "No LL1Out HCAL node... " << _nodename << std::endl;
      exit(1);
    }
  if (_verbosity) std::cout << __FILE__ << " "<< __LINE__<<" "<<std::endl;
  _trigger_bits = _ll1_hcal->GetTriggerBits();
  if (!_trigger_bits)
    {
      std::cout <<" no trigger bits..." <<std::endl;
      exit(1);
    }
  i = 0;
  for (auto iter = _trigger_bits->begin(); iter < _trigger_bits->end(); i++, ++iter)
    {
      b_trigger_bits_hcal.push_back((*iter));
    } 

  if (useCaloTowerBuilder)
    {
      _towers = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALIN");
      if (!_towers)
	{
	  cout << "No HCALIN towers..." <<endl;
	  exit(1);	  
	}
      int size;

      size = _towers->size(); //online towers should be the same!
      for (int channel = 0; channel < size;channel++)
	{
	  _tower = _towers->get_tower_at_channel(channel);
	  float energy = _tower->get_energy();
	  unsigned int towerkey = _towers->encode_key(channel);
	  int ieta = _towers->getTowerEtaBin(towerkey);
	  int iphi = _towers->getTowerPhiBin(towerkey);
	
	  b_hcalin_energy.push_back(energy);
	  b_hcalin_etabin.push_back(ieta);
	  b_hcalin_phibin.push_back(iphi);
      }

      _towers = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_HCALOUT");
      if (!_towers)
	{
	  cout << "No HCALOUT towers..." <<endl;
	  exit(1);	  
	}

      size = _towers->size(); //online towers should be the same!
      for (int channel = 0; channel < size;channel++)
	{
	  _tower = _towers->get_tower_at_channel(channel);
	  float energy = _tower->get_energy();
	  unsigned int towerkey = _towers->encode_key(channel);
	  int ieta = _towers->getTowerEtaBin(towerkey);
	  int iphi = _towers->getTowerPhiBin(towerkey);
	
	  b_hcalout_energy.push_back(energy);
	  b_hcalout_etabin.push_back(ieta);
	  b_hcalout_phibin.push_back(iphi);
      }

    }
  if (_verbosity) std::cout << __FILE__ << " "<< __LINE__<<" "<<std::endl;

  _tree_hcal->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}



void HCALEmulatorTreeMaker::GetNodes(PHCompositeNode* topNode)
{


}

int HCALEmulatorTreeMaker::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "HCALEmulatorTreeMaker::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HCALEmulatorTreeMaker::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "HCALEmulatorTreeMaker::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HCALEmulatorTreeMaker::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "HCALEmulatorTreeMaker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  std::cout<<"Total events: "<<_i_event<<std::endl;
  _f->Write();
  _f->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HCALEmulatorTreeMaker::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "HCALEmulatorTreeMaker::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

