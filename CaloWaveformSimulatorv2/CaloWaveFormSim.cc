#include "CaloWaveFormSim.h"
#include "WaveformContainerv1.h"
#include "ADCDefs.h"
// G4Hits includes
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitDefs.h>  // for hit_idbits

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>  // for genkey, keytype

// Tower includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerDefs.h>
// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>


#include "g4detectors/PHG4CylinderGeomContainer.h"
#include "g4detectors/PHG4CylinderGeom_Spacalv1.h"  // for PHG4CylinderGeom_Spaca...
#include "g4detectors/PHG4CylinderGeom_Spacalv3.h"
#include "g4detectors/PHG4CylinderCellGeomContainer.h"
#include "g4detectors/PHG4CylinderCellGeom_Spacalv1.h"
#include "g4detectors/PHG4FullProjSpacalCellReco.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>


#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtHitV1.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TMath.h>
#include <cassert>
#include <sstream>
#include <string>
#include <TF1.h>
#include <phool/onnxlib.h>

class TProfile;


TProfile *CaloWaveFormSim::h_template_emcal;
TProfile *CaloWaveFormSim::h_template_ohcal;
TProfile *CaloWaveFormSim::h_template_ihcal;
TProfile *CaloWaveFormSim::h_template_mbd;


using namespace std;

CaloWaveFormSim::CaloWaveFormSim(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , detector("ALL")
  , _verbose(0) 
{
  _gain_opts["LOW"] = GAIN::LOW;
  _gain_opts["HIGH"] = GAIN::HIGH;

  _raw_towers_ohcal = new TowerInfoContainerv1(TowerInfoContainer::DETECTOR::HCAL);
}

CaloWaveFormSim::~CaloWaveFormSim()
{
}

int CaloWaveFormSim::Init(PHCompositeNode*)
{
  rnd = new TRandom3(0);

  const char *noiseenv = getenv("CALOWAVEFORMSIM_NOISE");
  std::string noisefilename;
  if (noiseenv == nullptr)
    {
      const char *offline_main = getenv("MYINSTALL");
      assert(offline_main);
      noisefilename = offline_main;
      noisefilename += "/share/calowaveformsim/";
      
	switch(_noiselevel)
	  {
	  case 0 :
	    noisefilename += "/no_raddmgnoise.csv";
	    noise = new TTree("noise_norad", "tree");
	    break;
	  case 1:
	    noisefilename += "/low_raddmgnoise.csv";
	    noise = new TTree("noise_lowrad", "tree");
	    break;
	  case 2:
	    noisefilename += "/medium_raddmgnoise.csv";
	    noise = new TTree("noise_midrad", "tree");
	    break;
	  default:
	    noisefilename += "/no_raddmgnoise.csv";
	    noise = new TTree("noise_norad", "tree");
	    break;
	  } 
    }
  else 
    {
      noisefilename = noiseenv;
    }

  noise->ReadFile(noisefilename.c_str(), "a1:a2:a3:a4:a5:a6:a7:a8:a9:a10:a11:a12:a13:a14:a15:a16:a17:a18:a19:a20:a21:a22:a23:a24:a25:a26:a27:a28:a29:a30:a31");
	


  for (int i = 0; i < 31;i++)
    {
      noise->SetBranchAddress(Form("a%d",i+1),&noise_val[i]);
    }
  //----------------------------------------------------------------------------------------------------
  //Read in the template file, this currently points to a tim local area file, 
  //but a copy of this file is in the git repository.
  //----------------------------------------------------------------------------------------------------


  const char *emcalenv = getenv("CALOWAVEFORMSIM_TEMPLATE_CEMC");
  std::string emcalfilename;
  if (noiseenv == nullptr)
    {
      const char *offline_main = getenv("MYINSTALL");
      assert(offline_main);
      emcalfilename = offline_main;
      emcalfilename += "/share/calowaveformsim/testbeam_emcal_template.root";     
    }
  else 
    {
      emcalfilename = emcalenv;

    }
  const char *ihcalenv = getenv("CALOWAVEFORMSIM_TEMPLATE_IHCAL");
  std::string ihcalfilename;
  if (noiseenv == nullptr)
    {
      const char *offline_main = getenv("MYINSTALL");
      assert(offline_main);
      ihcalfilename = offline_main;
      ihcalfilename += "/share/calowaveformsim/testbeam_ihcal_template.root";     
    }
  else 
    {
      ihcalfilename = ihcalenv;

    }
  const char *ohcalenv = getenv("CALOWAVEFORMSIM_TEMPLATE_OHCAL");
  std::string ohcalfilename;
  if (noiseenv == nullptr)
    {
      const char *offline_main = getenv("MYINSTALL");
      assert(offline_main);
      ohcalfilename = offline_main;
      ohcalfilename += "/share/calowaveformsim/testbeam_ohcal_template.root";     
    }
  else 
    {
      ohcalfilename = ohcalenv;

    }
  const char *mbdenv = getenv("CALOWAVEFORMSIM_TEMPLATE_BBC");
  std::string mbdfilename;
  if (noiseenv == nullptr)
    {
      const char *offline_main = getenv("MYINSTALL");
      assert(offline_main);
      mbdfilename = offline_main;
      mbdfilename += "/share/calowaveformsim/laser_mbd_template.root";     
    }
  else 
    {
      mbdfilename = mbdenv;
    }


  TFile* fin1 = TFile::Open(emcalfilename.c_str());
  assert(fin1);
  assert(fin1->IsOpen());
  h_template_emcal = static_cast<TProfile*>(fin1->Get("waveform_template"));

  TFile* fin2 = TFile::Open(ihcalfilename.c_str());
  assert(fin2);
  assert(fin2->IsOpen());
  h_template_ihcal = static_cast<TProfile*>(fin2->Get("waveform_template"));

  TFile* fin3 = TFile::Open(ohcalfilename.c_str());
  assert(fin3);
  assert(fin3->IsOpen());
  h_template_ohcal = static_cast<TProfile*>(fin3->Get("waveform_template"));

  TFile* fin4 = TFile::Open(mbdfilename.c_str());
  assert(fin4);
  assert(fin4->IsOpen());
  h_template_mbd = static_cast<TProfile*>(fin4->Get("waveform_template"));

  light_collection_model.load_data_file(string(getenv("CALIBRATIONROOT")) + string("/CEMC/LightCollection/Prototype3Module.xml"),
					"data_grid_light_guide_efficiency", "data_grid_fiber_trans");


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

double CaloWaveFormSim::template_function_mbd(double *x, double *par)
{ 
  return par[0]*h_template_mbd->Interpolate(x[0]-par[1])+par[2];
}

double CaloWaveFormSim::template_function_emcal(double *x, double *par)
{ 
  return par[0]*h_template_emcal->Interpolate(x[0]-par[1])+par[2];
}

double CaloWaveFormSim::template_function_ihcal(double *x, double *par)
{ 
  return par[0]*h_template_ihcal->Interpolate(x[0]-par[1])+par[2];
}

double CaloWaveFormSim::template_function_ohcal(double *x, double *par)
{ 
  return par[0]*h_template_ohcal->Interpolate(x[0]-par[1])+par[2];
}

int CaloWaveFormSim::InitRun(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;
  CreateNodes(topNode);
  return 0;
}


void CaloWaveFormSim::CreateNodes(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
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
	  detNode = new PHCompositeNode("MBD");
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
      PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "CEMC"));
      if (!detNode)
	{
	  std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
	  detNode = new PHCompositeNode("CEMC");
	  dstNode->addNode(detNode);
	}
      
      WaveformContainerv1 *waveforms = findNode::getClass<WaveformContainerv1>(detNode, "WAVEFORMS_CEMC");
      if (!waveforms)
	{
	  waveforms = new WaveformContainerv1();
	  PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(waveforms, "WAVEFORMS_CEMC", "PHObject");
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


int CaloWaveFormSim::process_event(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
  process_g4hits(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloWaveFormSim::process_g4hits(PHCompositeNode* topNode)
{
  
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  if (IsDetector("EMCAL"))
    {  
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      waveforms_emcal = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_CEMC");
      if (!waveforms_emcal)
	{
	  std::cout << "Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}

      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      _layergeo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_CEMC");
      if (!_layergeo)
	{
	  std::cout << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - Could not locate sim geometry node "
		    <<  std::endl;
	  exit(1);
	}

      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      _seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_CEMC");
      if (!_seggeo)
	{
	  std::cout << "Seggeo node not found "<<std::endl;
	  exit(1);
	}

      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      _hits_emcal = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_CEMC");
      if (!_hits_emcal)
	{
	  std::cout << "Hits CEMC not found " << std::endl;
	  exit(1);
	}
    }
  if (IsDetector("HCALIN"))
    {
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      waveforms_hcalin = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALIN");
      if (!waveforms_hcalin)
	{
	  std::cout << "Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      _hits_ihcal = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_HCALIN");
      if (!_hits_ihcal)
	{
	  std::cout << "IHCAL hits not found - Fatal Error" << std::endl;
	  exit(1);
	}
      
    }
  if (IsDetector("HCALOUT"))
    {
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      waveforms_hcalout = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_HCALOUT");
      if (!waveforms_hcalout)
	{
	  std::cout << "Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}

      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      // _hits_ohcal = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_HCALOUT");
      // if (!_hits_ohcal)
      // 	{
      // 	  std::cout << "OHCal hits not found - Fatal Error" << std::endl;
      // 	  exit(1);
      // 	}
      _slats_ohcal = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_HCALOUT");
      if (!_slats_ohcal)
      	{
      	  std::cout << "OHCal hits not found - Fatal Error" << std::endl;
	  exit(1);
      	}
      
      
    }

  if (IsDetector("MBD"))
    {
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      waveforms_mbd = findNode::getClass<WaveformContainerv1>(topNode, "WAVEFORMS_MBD");
      if (!waveforms_mbd)
	{
	  std::cout << "Waveforms not found - Fatal Error" << std::endl;
	  exit(1);
	}

      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

      _mbdpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
      if (!_mbdpmts)
	{
	  std::cout << "mbd pmts not found - Fatal Error" << std::endl;
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
  //Load in the template function as a TF1
  //for use in waveform generation
  //---------------------------------------------------------

  TF1 *f_fit_mbd = new TF1("f_fit_mbd",template_function_mbd,0,31,3);
  f_fit_mbd->SetParameters(1,0,0);

  TF1 *f_fit_emcal = new TF1("f_fit_emcal",template_function_emcal,0,31,3);
  f_fit_emcal->SetParameters(1,0,0);

  TF1 *f_fit_ihcal = new TF1("f_fit_ihcal",template_function_ihcal,0,31,3);
  f_fit_ihcal->SetParameters(1,0,0);

  TF1 *f_fit_ohcal = new TF1("f_fit_ohcal",template_function_ohcal,0,31,3);
  f_fit_ohcal->SetParameters(1,0,0);

  //-----------------------------------------------------
  //Set the timeing in of the prompt 
  //signal peak to be 4 time samples into
  //the waveform
  //------------------------------------------------------
  float _shiftval_mbd = 4-f_fit_mbd->GetMaximumX();
  f_fit_mbd->SetParameters(1,_shiftval_mbd,0);

  float _shiftval_emcal = 4-f_fit_emcal->GetMaximumX();
  f_fit_emcal->SetParameters(1,_shiftval_emcal,0);

  float _shiftval_ihcal = 4-f_fit_ihcal->GetMaximumX();
  f_fit_ihcal->SetParameters(1,_shiftval_ihcal,0);

  float _shiftval_ohcal = 4-f_fit_ohcal->GetMaximumX();
  f_fit_ohcal->SetParameters(1,_shiftval_ohcal,0);


  if (IsDetector("EMCAL"))
    {
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      //-----------------------------------------------------------------------
      //Loop over G4Hits to build a waveform simulation
      //-----------------------------------------------------------------------

      const PHG4CylinderGeom *layergeom_raw = _layergeo->GetFirstLayerGeom();
      assert(layergeom_raw);
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      const PHG4CylinderGeom_Spacalv3 *layergeom = dynamic_cast<const PHG4CylinderGeom_Spacalv3 *>(layergeom_raw);
      assert(layergeom);
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

      PHG4CylinderCellGeom *geo_raw = _seggeo->GetFirstLayerCellGeom();
      PHG4CylinderCellGeom_Spacalv1 *geo = dynamic_cast<PHG4CylinderCellGeom_Spacalv1 *>(geo_raw);
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      PHG4HitContainer::ConstRange hit_range = _hits_emcal->getHits();

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
	  

	  {
	    const double z = 0.5 * (hit_iter->second->get_local_z(0) + hit_iter->second->get_local_z(1));
	    assert(not std::isnan(z));
	    light_yield *= light_collection_model.get_fiber_transmission(z);
	  }
	  {
	    const double x = it_tower->second.get_position_fraction_x_in_sub_tower(decoder.fiber_ID);
	    const double y = it_tower->second.get_position_fraction_y_in_sub_tower(decoder.fiber_ID);
	    light_yield *= light_collection_model.get_light_guide_efficiency(x, y);
	  }

	  //-------------------------------------------------------------------------
	  //Map the G4hits to the corresponding CEMC tower
	  //-------------------------------------------------------------------------
	  int etabin = etabinshort;
	  int phibin = phibin_cell;

	  //------------------------------------------------------------------------
	  // Map Calo Tower to channel number
	  //------------------------------------------------------------------------
	  int ADC = (phibin/8)*12 + (etabin/8);
	  int oADC = (ADC/12)*12 + (etabin<48 ? 5 - (ADC%12): (ADC%12));
	  int channelnumber = (oADC < 64 ? ADCDefs::emcadc[etabin%8][phibin%8] : ADCDefs::emcadc[7 - etabin%8][7 - phibin%8]);
	  int towernumber = oADC*64 + channelnumber;

	  //---------------------------------------------------------------------
	  //Convert the G4hit into a waveform contribution
	  //---------------------------------------------------------------------
	  // float t0 = 0.5*(hit_iter->second->get_t(0)+hit_iter->second->get_t(1)) / 16.66667;   //Average of g4hit time downscaled by 16.667 ns/time sample 
	  float t0 = (hit_iter->second->get_t(0)) / 16.66667;   //Place waveform at the starting time of the G4hit, avoids issues caused by excessively long lived g4hits
	  float tmax = 16.667*_nsamples;
	  float tmin = -20;
	  f_fit_emcal->SetParameters(light_yield*26000,_shiftval_emcal+t0,0);            //Set the waveform template to match the expected signal from such a hit

	  //-------------------------------------------------------------------------------------------------------------
	  //For each tower add the new waveform contribution to the total waveform
	  //-------------------------------------------------------------------------------------------------------------
	  if (hit_iter->second->get_edep()*26000 > 1 && hit_iter->second->get_t(1) >= tmin && hit_iter->second->get_t(0) <= tmax)
	    {
	      std::cout <<"here: "<< towernumber << ": " ;
	      for (int j = 0; j < _nsamples;j++) // 16 is the number of time samples
		{
		  std::cout << " "<<f_fit_emcal->Eval(j);
		  m_waveform_emcal[towernumber][j] += f_fit_emcal->Eval(j);
		}
	      std::cout << " " <<std::endl;
	    }
	}
    }



  //--------------
  // do IHCAL
  //--------------
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  if (IsDetector("HCALIN"))
    {
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      //-----------------------------------------------------------------------
      //Loop over G4Hits to build a waveform simulation
      //-----------------------------------------------------------------------
      PHG4HitContainer::ConstRange hit_range = _hits_ihcal->getHits();
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
	  //-------------------------------------------------------------------------
	  //Map the G4hits to the corresponding CEMC tower
	  //-------------------------------------------------------------------------
      
	  int etabin = towereta;
	  int phibin =towerphi;
	  //------------------------------------------------------------------------
	  // Map Calo Tower to channel number
	  //------------------------------------------------------------------------
	  int ADC = (phibin/8)*3 + (etabin/8);
	  int channelnumber = ADCDefs::hcaladc[etabin%8][(phibin%8)%2] + 16*((phibin%8)/2);
	  int towernumber = ADC*64 + channelnumber;

     
	  //---------------------------------------------------------------------
	  //Convert the G4hit into a waveform contribution
	  //---------------------------------------------------------------------
	  // float t0 = 0.5*(hit_iter->second->get_t(0)+hit_iter->second->get_t(1)) / 16.66667;   //Average of g4hit time downscaled by 16.667 ns/time sample 
	  float t0 = (hit_iter->second->get_t(0)) / 16.66667;   //Place waveform at the starting time of the G4hit, avoids issues caused by excessively long lived g4hits
	  float tmax = 16.667*_nsamples;
	  float tmin = -20;
	  f_fit_ihcal->SetParameters(light_yield*2600,_shiftval_ihcal+t0,0);            //Set the waveform template to match the expected signal from such a hit
	  //-------------------------------------------------------------------------------------------------------------
	  //For each tower add the new waveform contribution to the total waveform
	  //-------------------------------------------------------------------------------------------------------------
	  if (hit_iter->second->get_edep()*2600 > 1 &&hit_iter->second->get_t(1) >= tmin && hit_iter->second->get_t(0) <= tmax )
	    {
	      for (int j = 0; j < _nsamples;j++) // 16 is the number of time samples
		{
		  m_waveform_ihcal[towernumber][j] +=f_fit_ihcal->Eval(j);
		}
	    }
	}
    }
  


  //--------------
  // do OHCAL
  //--------------

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  if (IsDetector("HCALOUT"))
    {
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
      _raw_towers_ohcal->Reset();
      //-----------------------------------------------------------------------
      //Loop over G4Hits to build a waveform simulation
      //-----------------------------------------------------------------------
      PHG4CellContainer::ConstRange cell_range = _slats_ohcal->getCells();
      for (PHG4CellContainer::ConstIterator cell_iter = cell_range.first; cell_iter != cell_range.second; cell_iter++)
	{
	  PHG4Cell *cell = cell_iter->second;
	  short twrrow = PHG4CellDefs::ScintillatorSlatBinning::get_row(cell->get_cellid());
	  
	  double light_yield = cell->get_light_yield();
      
	  TowerInfo *towerinfo;
	  
	  unsigned int etabin = PHG4CellDefs::ScintillatorSlatBinning::get_column(cell->get_cellid());
	  unsigned int phibin = twrrow;
	  
	  unsigned int towerkey = (etabin << 16U) + phibin;

	  towerinfo = _raw_towers_ohcal->get_tower_at_key(towerkey);
	  if (!towerinfo)
	    {
	      std::cout << __PRETTY_FUNCTION__ << ": missing towerkey = " << towerkey << " in m_TowerInfoContainer!";
	      exit(1);
	    }
	  else
	    {
	      towerinfo->set_energy(towerinfo->get_energy() + light_yield);
	    }

	}


      for (auto tower_iter = 0 ; tower_iter < static_cast<int>(_raw_towers_ohcal->size()); tower_iter++)
	{
	  int key = _raw_towers_ohcal->encode_hcal(tower_iter);
	  TowerInfo *towerinfo = _raw_towers_ohcal->get_tower_at_key(key);
	  unsigned int phibin = _raw_towers_ohcal->getTowerPhiBin(key);
	  unsigned int etabin = _raw_towers_ohcal->getTowerEtaBin(key);
	  int ADC = (phibin/8)*3 + (etabin/8);
	  int channelnumber = ADCDefs::hcaladc[etabin%8][(phibin%8)%2] + 16*((phibin%8)/2);
	  int towernumber = ADC*64 + channelnumber;
      
	  //---------------------------------------------------------------------
	  //Convert the G4hit into a waveform contribution
	  //---------------------------------------------------------------------
	  // float t0 = 0.5*(hit_iter->second->get_t(0)+hit_iter->second->get_t(1)) / 16.66667;   //Average of g4hit time downscaled by 16.667 ns/time sample 

	  float adc_signal = towerinfo->get_energy() * _hcalout_lightyield_to_ADC;
	  float t0 = 0;
	  float tmax =16.667*_nsamples ;
	  float tmin = -20;
	  f_fit_ohcal->SetParameters(adc_signal,_shiftval_ohcal+t0,0);            //Set the waveform template to match the expected signal from such a hit

	  //-------------------------------------------------------------------------------------------------------------
	  //For each tower add the new waveform contribution to the total waveform
	  //-------------------------------------------------------------------------------------------------------------
      
	  if (towerinfo->get_energy()*5000 > 1 && towerinfo->get_time() >= tmin && towerinfo->get_time() <= tmax)
	    {
	      for (int j = 0; j < _nsamples;j++) // _nsamples is the number of time samples
		{
		  m_waveform_ohcal[towernumber][j] +=f_fit_ohcal->Eval(j);
		}
	    }
	}
    }
  
  if (IsDetector("MBD"))
    {
      //---------------------------------------------------------------------
      // Get PMT charges
      //---------------------------------------------------------------------
      
      int ich;
      int npmt = _mbdpmts->get_npmt();

      for (ich = 0; ich < npmt ; ich++)
	{
	  MbdPmtHit *tmp_pmt = _mbdpmts->get_pmt( ich);
	  int ipmt = ich;
	  int ipmtq = 32 + ((ipmt/32)*64) + (ipmt%32);
	  int ipmtt = ((ipmt/32)*64) + (ipmt%32);

	  float adc  = tmp_pmt->get_q()*70;

	  float tdc0 = tmp_pmt->get_time();
	  float tdc0_dig = 13000.* (1 - (1./22.5)*(tdc0 - 2.5));
	  float t0 = (tdc0) / 16.66667;   //Place waveform at the starting time of the G4hit, avoids issues caused by excessively long lived g4hits

	  f_fit_mbd->SetParameters(adc,_shiftval_mbd+t0,0);            //Set the wavefor
	  
	  for (int j = 0; j < _nsamples;j++) // 16 is the number of time samples
	    {
	      m_waveform_mbd[ipmtq][j] += f_fit_mbd->Eval(j);
	    }
	    
	  f_fit_mbd->SetParameters(tdc0_dig,_shiftval_mbd+t0,0);            //Set the wavefor

	  for (int j = 0; j < _nsamples;j++) // _nsamples is the number of time samples
	    {
	      m_waveform_mbd[ipmtt][j] += f_fit_mbd->Eval(j);
	    }

	}
    }

  
  //----------------------------------------------------------------------------------------
  //For each tower loop over add a noise waveform 
  //from cosmic data, gain is too high but is corrected for
  //----------------------------------------------------------------------------------------

  //-----------------------------
  // do noise for EMCal:
  //-----------------------------
  if (IsDetector("EMCAL"))
    {
      for (int i = 0; i < 24576;i++)
	{
	  int noise_waveform  = (int)rnd->Uniform(0,1990);
	  noise->GetEntry(noise_waveform);
      
	  std::vector<int> *wave = new vector<int>;
	  wave->clear();

	  for (int k = 0; k < _nsamples;k++)
	    {
	      m_waveform_emcal[i][k] = m_waveform_emcal[i][k]+(noise_val[k]-1500)/16.0+1500;
	      wave->push_back(static_cast<int>(m_waveform_emcal[i][k]));
	      if (m_waveform_emcal[i][k] > 1700.) cout << "HIT " << i << std::endl;
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
	  int noise_waveform  = (int)rnd->Uniform(0,1990);
	  std::vector<int> *wave = new vector<int>;
	  wave->clear();
	  noise->GetEntry(noise_waveform);     
      
	  for (int k = 0; k < _nsamples;k++)
	    {
	      m_waveform_ihcal[i][k] = m_waveform_ihcal[i][k]+(noise_val[k]-1500)/16.0+1500;
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
	  int noise_waveform  = (int)rnd->Uniform(0,1990);
	  noise->GetEntry(noise_waveform);
	  std::vector<int> *wave = new vector<int>;
	  wave->clear();      
      
	  for (int k = 0; k < _nsamples;k++)
	    {
	      m_waveform_ohcal[i][k] = m_waveform_ohcal[i][k]+(noise_val[k]-1500)/16.0+1500;
	      wave->push_back(static_cast<int>(m_waveform_ohcal[i][k]));
	    }
	  waveforms_hcalout->set_waveform_at_channel(i, wave);
	}
    }

  if (IsDetector("MBD"))
    {
      for (int i = 0; i < 256;i++)
	{
	  int noise_waveform  = (int)rnd->Uniform(0,1990);
	  noise->GetEntry(noise_waveform);
      
	  std::vector<int> *wave = new vector<int>;
	  wave->clear();
      
	  for (int k = 0; k < _nsamples;k++)
	    {
	      m_waveform_mbd[i][k] = m_waveform_mbd[i][k]+(noise_val[k]-1500)/16.0+1500;
	      wave->push_back(static_cast<int>(m_waveform_mbd[i][k]));
	    }
	  waveforms_mbd->set_waveform_at_channel(i, wave);
	}
    }


  return Fun4AllReturnCodes::EVENT_OK;
}



int CaloWaveFormSim::End(PHCompositeNode* topNode)
{
  return 0;
}

bool CaloWaveFormSim::IsDetector(const std::string &det)
{
  if (strcmp("ALL", detector.c_str()) == 0) return true;
  if (strcmp(det.c_str(), detector.c_str()) == 0) return true;
  return false;
}
