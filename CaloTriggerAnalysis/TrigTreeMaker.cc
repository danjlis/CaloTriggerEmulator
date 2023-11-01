#include "TrigTreeMaker.h"
#include <vector>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitDefs.h>  // for hit_idbits
#include <g4detectors/PHG4CylinderGeom_Spacalv3.h>
// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>  // for genkey, keytype

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

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <iostream>

#include <map>

//____________________________________________________________________________..
TrigTreeMaker::TrigTreeMaker(const std::string &name, const std::string &outfilename):
  SubsysReco(name)
{
  _foutname = outfilename;  
  _verbosity = 0;
}

//____________________________________________________________________________..
TrigTreeMaker::~TrigTreeMaker()
{

}

//____________________________________________________________________________..
int TrigTreeMaker::Init(PHCompositeNode *topNode)
{

  _f = new TFile( _foutname.c_str(), "RECREATE");

  std::cout << " making a file = " <<  _foutname.c_str() << " , _f = " << _f << std::endl;
  
  _tree = new TTree("ttree","a persevering date tree");

  _tree->Branch("primpt",&m_primpt);
  _tree->Branch("primeta",&m_primeta);
  _tree->Branch("primphi",&m_primphi);
  _tree->Branch("waveform_emcal",&m_waveforms_cemc,"waveform_emcal[24576][16]/I");
  _tree->Branch("waveform_ihcal",&m_waveforms_hcalin,"waveform_ihcal[1536][16]/I");
  _tree->Branch("waveform_ohcal",&m_waveforms_hcalout,"waveform_ohcal[1536][16]/I");
  _tree->Branch("waveform_bbc",&m_waveforms_bbc,"waveform_bbc[256][16]/I");
  _tree->Branch("lightyield_emcal",&m_lightyield_cemc,"lightyield_emcal[24576]/F");
  _tree->Branch("edep_emcal",&m_edep_cemc,"edep_emcal[24576]/F");
  _tree->Branch("lightyield_hcalin",&m_lightyield_hcalin,"lightyield_hcalin[24576]/F");
  _tree->Branch("edep_hcalin",&m_edep_hcalin,"edep_hcalin[24576]/F");
  _tree->Branch("lightyield_hcalout",&m_lightyield_hcalout,"lightyield_hcalout[1536]/F");
  _tree->Branch("edep_hcalout",&m_edep_hcalout,"edep_hcalout[1536]/F");

  _tree->Branch("npmt_bbc",&m_npmt,"npmt_bbc/I");
  _tree->Branch("adc_bbc",&m_adc_bbc,"adc_bbc[128]/F");
  _tree->Branch("tdc0_bbc",&m_tdc0_bbc,"tdc0_bbc[128]/F");
  _tree->Branch("tdc1_bbc",&m_tdc1_bbc,"tdc1_bbc[128]/F");
  _tree->Branch("bbc_vtx_z", &m_bbc_vtx_z,"bbc_vtx_z/F");
  _tree->Branch("bbc_vtx_t0", &m_bbc_vtx_t0,"bbc_vtx_t0/F");
  _tree->Branch("bbc_vtx_n", &m_bbc_vtx_n,"bbc_vtx_n/I");

  _tree->Branch("mbd_trigger_bits",&m_mbd_trigger_bits,"mbd_trigger_bits/I");
  _tree->Branch("mbd_trigger_words",&m_mbd_trigger_words,"mbd_trigger_words[13]/I");

  _i_event = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrigTreeMaker::InitRun(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

void TrigTreeMaker::SetVerbosity(int verbo){
  _verbosity = verbo;
  return;
}

void TrigTreeMaker::reset_tree_vars()
{
  m_primpt.clear();
  m_primeta.clear();
  m_primphi.clear();
  m_npmt = 0;
  for (int i = 0; i < 24576; i++)
    {
      m_lightyield_cemc[i] = 0.;
      m_edep_cemc[i] = 0.;
      for (int j = 0; j < 16; j++) m_waveforms_cemc[i][j] = 0;
    }
  for (int i = 0; i < 1536; i++)
    {
      m_lightyield_hcalin[i] = 0.;
      m_edep_hcalin[i] = 0.;
      m_lightyield_hcalout[i] = 0.;
      m_edep_hcalout[i] = 0.;

      for (int j = 0; j < 16; j++) {
	m_waveforms_hcalout[i][j] = 0;
	m_waveforms_hcalin[i][j] = 0;
      }
    }
  for (int i = 0; i < 128; i++)
    {
      m_adc_bbc[i] = 0.;
      m_tdc0_bbc[i] = 0.;
      m_tdc1_bbc[i] = 0.;
      for (int j = 0; j < 16; j++) m_waveforms_bbc[i][j] = 0;
    }

  m_mbd_trigger_bits = 0;
  for (int i = 0; i < 13; i++)
    {
      m_mbd_trigger_words[i] = 0;
    }

  m_bbc_vtx_z = 0.0;
  m_bbc_vtx_t0 = 0.0;
  m_bbc_vtx_n = 0;

  return;
}
void TrigTreeMaker::GetNodes(PHCompositeNode* topNode)
{

  if (_verbosity) std::cout << __FUNCTION__ << std::endl;

  _waveforms_cemc = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_CEMC");

  if (!_waveforms_cemc) 
    {
      std::cout << "No CEMC Waveforms found... " << std::endl;
      exit(1);
    }

  _waveforms_hcalout = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALOUT");

  if (!_waveforms_hcalout) 
    {
      std::cout << "No HCALOUT Waveforms found... " << std::endl;
      exit(1);
    }

  _waveforms_hcalin = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALIN");

  if (!_waveforms_hcalin) 
    {
      std::cout << "No HCALIN Waveforms found... " << std::endl;
      exit(1);
    }

  _waveforms_bbc = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_MBD");

  if (!_waveforms_bbc) 
    {
      std::cout << "No BBC Waveforms found... " << std::endl;
      exit(1);
    }


  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truthinfo) 
    {
      std::cout << "No Truthinfo Waveforms found... " << std::endl;
      exit(1);
    }

  _hits_cemc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_CEMC");
  if (!_hits_cemc) 
    {
      std::cout << "No CEMC HITS Waveforms found... " << std::endl;
      exit(1);
    }

  _layergeo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_CEMC");
  if (!_layergeo) 
    {
      std::cout << "No Layergeometry Waveforms found... " << std::endl;
      exit(1);
    }
  _seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_CEMC");
  if (!_seggeo) 
    {
      std::cout << " No CellGeom found... " << std::endl;
      exit(1);
    }

  _hits_hcalin = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_HCALIN");
  if (!_hits_hcalin) 
    {
      std::cout << "No HCALIN HITS Waveforms found... " << std::endl;
      exit(1);
    }

  _hits_hcalout = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_HCALOUT");
  if (!_hits_hcalout) 
    {
      std::cout << "No HCALOUT HITS Waveforms found... " << std::endl;
      exit(1);
    }

  _bbc_pmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (!_bbc_pmts) 
    {
      std::cout << "MBD Pmt Container not found" << std::endl;
      exit(1);
    }

  _ll1_mbd = findNode::getClass<LL1Outv1>(topNode, "LL1Out_MBD");

  if (!_ll1_mbd) 
    {
      std::cout << "No LL1Out MBD node... " << std::endl;
      exit(1);
    }

  _bbcvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  if (!_bbcvertexmap)
    {
      std::cout << "No MBD Vertex node... " << std::endl;
      exit(1);
    }
  return;

}


void TrigTreeMaker::process_waveforms()
{
  if (_waveforms_cemc)
    {
      WaveformContainerv1::Range begin_end = _waveforms_cemc->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;
      for (; iwave != begin_end.second; ++iwave)
	{
	  std::vector<int> wave = *(iwave->second);
	  for (int i = 0; i < 16; i++)
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
	  for (int i = 0; i < 16; i++)
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
	  for (int i = 0; i < 16; i++)
	    {
	      m_waveforms_hcalin[iwave->first][i] = wave.at(i);
	    }
	}
    }
  if (_waveforms_bbc)
    {
      WaveformContainerv1::Range begin_end = _waveforms_bbc->getWaveforms();
      WaveformContainerv1::Iter iwave = begin_end.first;
      for (; iwave != begin_end.second; ++iwave)
	{
	  std::vector<int> wave = *(iwave->second);
	  for (int i = 0; i < 16; i++)
	    {
	      m_waveforms_bbc[iwave->first][i] = wave.at(i);
	    }
	}
    }

}

void TrigTreeMaker::process_truth()
{

  PHG4TruthInfoContainer::Range range = _truthinfo->GetPrimaryParticleRange();

  if (_verbosity) std::cout<<"Truth Particles: (pt/eta/phi) - > (pid/track id/vtx id/parent id/primary id/ barcode)"<<std::endl;

  for (auto iter = range.first; iter != range.second; ++iter) {

    PHG4Particle* g4particle = iter->second;
    
    if ( _truthinfo->isEmbeded (g4particle->get_track_id ()) != 1) continue;                                                                                                                     

    TLorentzVector t;
    t.SetPxPyPzE (g4particle->get_px (), g4particle->get_py (), g4particle->get_pz (), g4particle->get_e ());
    float ipt = t.Pt();
    float ieta = t.Eta();
    float iphi = t.Phi();
    if (_verbosity) std::cout<<"( "<< ipt<<" / "<< ieta << " / "<< iphi <<" ) "<<std::endl;

    m_primpt.push_back(ipt);
    m_primeta.push_back(ieta);
    m_primphi.push_back(iphi);
    
  }
  
  //


  if (_verbosity) std::cout <<  __FILE__ << __FUNCTION__<<__LINE__ << std::endl;
  if (_hits_cemc)
  {
    if (_verbosity) std::cout <<  __FILE__ << __FUNCTION__<<__LINE__ << std::endl;

    //-----------------------------------------------------------------------
    //Loop over G4Hits to build a waveform simulation
    //-----------------------------------------------------------------------
    const PHG4CylinderGeom *layergeom_raw = _layergeo->GetFirstLayerGeom();
    assert(layergeom_raw);
    if (_verbosity) std::cout <<  __FILE__ << __FUNCTION__<<__LINE__ << std::endl;
    const PHG4CylinderGeom_Spacalv3 *layergeom =
      dynamic_cast<const PHG4CylinderGeom_Spacalv3 *>(layergeom_raw);
    assert(layergeom);
    if (_verbosity) std::cout <<  __FILE__ << __FUNCTION__<<__LINE__ << std::endl;

    PHG4CylinderCellGeom *geo_raw = _seggeo->GetFirstLayerCellGeom();
    if (_verbosity) std::cout <<  __FILE__ << __FUNCTION__<<__LINE__ << std::endl;
    PHG4CylinderCellGeom_Spacalv1 *geo = dynamic_cast<PHG4CylinderCellGeom_Spacalv1 *>(geo_raw);
    if (_verbosity) std::cout <<  __FILE__ << __FUNCTION__<<__LINE__ << std::endl;
    PHG4HitContainer::ConstRange hit_range = _hits_cemc->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
      {
	//-----------------------------------------------------------------
	//Extract position information for each G4hit 
	//information for each G4hit in the calorimeter
	//-----------------------------------------------------------------

	int scint_id = hit_iter->second->get_scint_id();
	PHG4CylinderGeom_Spacalv3::scint_id_coder decoder(scint_id);
	std::pair<int, int> tower_z_phi_ID = layergeom->get_tower_z_phi_ID(decoder.tower_ID, decoder.sector_ID);
	const int &tower_ID_z = tower_z_phi_ID.first;
	const int &tower_ID_phi = tower_z_phi_ID.second;

	PHG4CylinderGeom_Spacalv3::tower_map_t::const_iterator it_tower =
	  layergeom->get_sector_tower_map().find(decoder.tower_ID);
	assert(it_tower != layergeom->get_sector_tower_map().end());
      
	const int etabin_cell = geo->get_etabin_block(tower_ID_z);  // block eta bin
	const int sub_tower_ID_x = it_tower->second.get_sub_tower_ID_x(decoder.fiber_ID);
	const int sub_tower_ID_y = it_tower->second.get_sub_tower_ID_y(decoder.fiber_ID);
	unsigned short etabinshort = etabin_cell * layergeom->get_n_subtower_eta() + sub_tower_ID_y;
	unsigned short phibin_cell = tower_ID_phi * layergeom->get_n_subtower_phi() + sub_tower_ID_x;
    
	//----------------------------------------------------------------------------------------------------
	//Extract light yield from g4hit and correct for light collection efficiency
	//----------------------------------------------------------------------------------------------------
	double light_yield = hit_iter->second->get_light_yield();
	double edep = hit_iter->second->get_edep();
	//-------------------------------------------------------------------------
	int etabin = etabinshort;
	int phibin = phibin_cell;

	//------------------------------------------------------------------------
	// Map Calo Tower to channel number
	//------------------------------------------------------------------------
	int ADC = (phibin/8)*12 + (etabin/8);
	int oADC = (etabin<48 ? ADC - 6*(ADC/12):192 + (5 - (ADC%6)) + (ADC/12)*6);
	int channelnumber = (oADC < 64 ? ADCDefs::emcadc[etabin%8][phibin%8] : ADCDefs::emcadc[7 - etabin%8][7 - phibin%8]);
	int towernumber = oADC*64 + channelnumber;
      
	m_lightyield_cemc[towernumber] = light_yield;
	m_edep_cemc[towernumber] = edep;
      
      }
  }  
  

  if (_hits_hcalin)
    {
      //-----------------------------------------------------------------------
      //Loop over G4Hits to build a waveform simulation
      //-----------------------------------------------------------------------
      PHG4HitContainer::ConstRange hit_range = _hits_hcalin->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
	{
	  //-----------------------------------------------------------------
	  //Extract position information for each G4hit 
	  //information for each G4hit in the calorimeter
	  //-----------------------------------------------------------------
	  short icolumn = hit_iter->second->get_scint_id();
	  int introw = (hit_iter->second->get_hit_id() >> PHG4HitDefs::hit_idbits);

	  if (introw >= ROWDIM || introw < 0)
	    {
	      std::cout << __PRETTY_FUNCTION__ << " row " << introw
			<< " exceed array size: " << ROWDIM
			<< " adjust ROWDIM and recompile" << std::endl;
	      exit(1);
	    }
	  int towerphi = introw/4;
	  int towereta = icolumn;

	  //----------------------------------------------------------------------------------------------------
	  //Extract light yield from g4hit and correct for light collection efficiency
	  //----------------------------------------------------------------------------------------------------
	  double light_yield = hit_iter->second->get_light_yield();  //raw_light_yield has no MEPHI maps applied, light_yield aoppplies the maps change at some point     
	  double edep = hit_iter->second->get_edep();

	  int etabin = towereta;
	  int phibin =towerphi;
	  //------------------------------------------------------------------------
	  // Map Calo Tower to channel number
	  //------------------------------------------------------------------------
	  int ADC = (phibin/8)*3 + (etabin/8);
	  int channelnumber = ADCDefs::hcaladc[etabin%8][(phibin%8)%2] + 16*((phibin%8)/2);
	  int towernumber = ADC*64 + channelnumber;
      
	  m_lightyield_hcalin[towernumber] = light_yield;
	  m_edep_hcalin[towernumber] = edep;
      
	}
    }

  if (_hits_hcalout)
    {
      //-----------------------------------------------------------------------
      //Loop over G4Hits to build a waveform simulation
      //-----------------------------------------------------------------------
      PHG4HitContainer::ConstRange hit_range = _hits_hcalout->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
	{
	  //-----------------------------------------------------------------
	  //Extract position information for each G4hit 
	  //information for each G4hit in the calorimeter
	  //-----------------------------------------------------------------
	  short icolumn = hit_iter->second->get_scint_id();
	  int introw = (hit_iter->second->get_hit_id() >> PHG4HitDefs::hit_idbits);


	  if (introw >= ROWDIM || introw < 0)
	    {
	      std::cout << __PRETTY_FUNCTION__ << " row " << introw
			<< " exceed array size: " << ROWDIM
			<< " adjust ROWDIM and recompile" << std::endl;
	      exit(1);
	    }


	  int towerphi = introw/5;
	  int towereta = icolumn;

	  //----------------------------------------------------------------------------------------------------
	  //Extract light yield from g4hit and correct for light collection efficiency
	  //----------------------------------------------------------------------------------------------------

	  double light_yield = hit_iter->second->get_light_yield();  //raw_light_yield has no MEPHI maps applied, light_yield aoppplies the maps change at some point
	  double edep = hit_iter->second->get_edep();
	  //-------------------------------------------------------------------------
	  //Map the G4hits to the corresponding CEMC tower
	  //-------------------------------------------------------------------------
      
	  int etabin = towereta;
	  int phibin = towerphi;
	  //------------------------------------------------------------------------
	  // Map Calo Tower to channel number
	  //------------------------------------------------------------------------
	  int ADC = (phibin/8)*3 + (etabin/8);
	  int channelnumber = ADCDefs::hcaladc[etabin%8][(phibin%8)%2] + 16*((phibin%8)/2);
	  int towernumber = ADC*64 + channelnumber;

	  m_lightyield_hcalout[towernumber] = light_yield;
	  m_edep_hcalout[towernumber] = edep;
	}
    } 

  if (_bbc_pmts)
    {
      //---------------------------------------------------------------------
      // Get PMT charges
      //---------------------------------------------------------------------

      int ich;
      m_npmt = _bbc_pmts->get_npmt();

      for (ich = 0; ich < m_npmt ; ich++)
        {
          short ipmt = _bbc_pmts->get_pmt( ich);
          float adc  = _bbc_pmts->get_adc( ich);
          float tdc0 = _bbc_pmts->get_tdc0(ich);
          float tdc1 = _bbc_pmts->get_tdc1(ich);

	  m_adc_bbc[ipmt] = adc;
	  m_tdc0_bbc[ipmt] = tdc0;
	  m_tdc1_bbc[ipmt] = tdc1;
	}
    }

  if (_ll1_mbd)
    {
      m_mbd_trigger_bits = _ll1_mbd->GetTriggerBits(11);
      for (int i = 0 ; i < 13; i++)
	{
	  m_mbd_trigger_words[i] = _ll1_mbd->GetTriggerWord(i,11);
	}
    }

  if (_bbcvertexmap)
    {
      m_bbc_vtx_n = 0;
      for (auto it = _bbcvertexmap->begin(); it != _bbcvertexmap->end(); ++it)
	{
	  BbcVertex *v = (*it).second;
	  m_bbc_vtx_z = v->get_z();
	  m_bbc_vtx_t0 = v->get_t();
	  m_bbcn_vtx_n++;
	}
    }
  return;
}

int TrigTreeMaker::process_event(PHCompositeNode *topNode)
{
  _i_event++;
  if (_i_event%((_verbosity > 0)? 1: 100)== 0) {
    std::cout<<"------------------------------------------------"<<std::endl;
    std::cout<<"Event "<<_i_event<<std::endl;
    std::cout<<"------------------------------------------------"<<std::endl;
  }
  GetNodes(topNode);
  reset_tree_vars();
  process_truth();
  process_waveforms();

  _tree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}


//int TrigTreeMaker::CreateNode(PHCompositeNode *topNode)
//{

//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
int TrigTreeMaker::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "TrigTreeMaker::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrigTreeMaker::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "TrigTreeMaker::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrigTreeMaker::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "TrigTreeMaker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  std::cout<<"Total events: "<<_i_event<<std::endl;
  _f->Write();
  _f->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrigTreeMaker::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "TrigTreeMaker::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void TrigTreeMaker::Print(const std::string &what) const
{
  std::cout << "TrigTreeMaker::Print(const std::string &what) const Printing info for " << what << std::endl;
}
