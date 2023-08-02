#include "CaloWaveFormSim.h"

// G4Hits includes

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <phool/getClass.h>
#include <TProfile.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TMath.h>
#include <cassert>
#include <sstream>
#include <string>
#include <TF1.h>
#include <phool/onnxlib.h>

using namespace std;
TProfile* h_template = nullptr;
TProfile* h_template_ihcal = nullptr;
TProfile* h_template_ohcal = nullptr;


#define ROWDIM 320
#define COLUMNDIM 27



double CaloWaveFormSim::template_function(double *x, double *par)
{
  Double_t v1 = par[0]*h_template->Interpolate(x[0]-par[1])+par[2];
  return v1;
}

double CaloWaveFormSim::template_function_ihcal(double *x, double *par)
{
  Double_t v1 = par[0]*h_template_ihcal->Interpolate(x[0]-par[1])+par[2];
  return v1;
}

double CaloWaveFormSim::template_function_ohcal(double *x, double *par)
{
  Double_t v1 = par[0]*h_template_ohcal->Interpolate(x[0]-par[1])+par[2];
  return v1;
}



CaloWaveFormSim::CaloWaveFormSim(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , detector("CEMC")
  , outfilename(filename)
  , hm(nullptr)
  , outfile(nullptr)
  , g4hitntuple(nullptr)
{
  _verbose = 0;
}

CaloWaveFormSim::~CaloWaveFormSim()
{
  delete hm;
  delete g4hitntuple;
}


int CaloWaveFormSim::InitRun(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;
  CreateNodes(topNode);
  return 0;
}

int CaloWaveFormSim::Init(PHCompositeNode*)
{
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;

  rnd = new TRandom3(0);
  //----------------------------------------------------------------------------------------------------
  //Read in the noise file, this currently points to a tim local area file, 
  //but a copy of this file is in the git repository.
  //----------------------------------------------------------------------------------------------------
  noise_midrad = new TTree("noise_midrad", "tree");
  // noise->ReadFile("/gpfs/mnt/gpfs02/sphenix/user/trinn/fitting_algorithm_playing/slowneutronsignals/CaloAna/noise.csv", "a1:a2:a3:a4:a5:a6:a7:a8:a9:a10:a11:a12:a13:a14:a15:a16:a17:a18:a19:a20:a21:a22:a23:a24:a25:a26:a27:a28:a29:a30:a31");
  noise_midrad->ReadFile("/gpfs/mnt/gpfs02/sphenix/user/trinn/sPHENIX_emcal_cosmics_sector0/noise_waveforms/medium_raddmgnoise.csv", "a1:a2:a3:a4:a5:a6:a7:a8:a9:a10:a11:a12:a13:a14:a15:a16:a17:a18:a19:a20:a21:a22:a23:a24:a25:a26:a27:a28:a29:a30:a31");

  noise_lowrad = new TTree("noise_lowrad", "tree");
  noise_lowrad->ReadFile("/gpfs/mnt/gpfs02/sphenix/user/trinn/sPHENIX_emcal_cosmics_sector0/noise_waveforms/low_raddmgnoise.csv", "a1:a2:a3:a4:a5:a6:a7:a8:a9:a10:a11:a12:a13:a14:a15:a16:a17:a18:a19:a20:a21:a22:a23:a24:a25:a26:a27:a28:a29:a30:a31");


  noise_norad = new TTree("noise_norad", "tree");
  noise_norad->ReadFile("/gpfs/mnt/gpfs02/sphenix/user/trinn/sPHENIX_emcal_cosmics_sector0/noise_waveforms/no_raddmgnoise.csv", "a1:a2:a3:a4:a5:a6:a7:a8:a9:a10:a11:a12:a13:a14:a15:a16:a17:a18:a19:a20:a21:a22:a23:a24:a25:a26:a27:a28:a29:a30:a31");

  for (int i = 0; i < 31;i++)
    {
      noise_midrad->SetBranchAddress(Form("a%d",i+1),&noise_val_midrad[i]);
      noise_lowrad->SetBranchAddress(Form("a%d",i+1),&noise_val_lowrad[i]);
      noise_norad->SetBranchAddress(Form("a%d",i+1),&noise_val_norad[i]);
    }
  
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  g4hitntuple = new TTree("tree", "tree");

  g4hitntuple->Branch("primpt",&m_primpt);
  g4hitntuple->Branch("primeta",&m_primeta);
  g4hitntuple->Branch("primphi",&m_primphi);
  g4hitntuple->Branch("tedep_emcal",&m_tedep,"tedep_emcal[24576]/F");
  g4hitntuple->Branch("tedep_ihcal",&m_tedep_ihcal,"tedep_ihcal[1536]/F");
  g4hitntuple->Branch("tedep_ohcal",&m_tedep_ohcal,"tedep_ohcal[1536]/F");
  g4hitntuple->Branch("extractedadc_emcal",& m_extractedadc,"extractedadc_emcal[24576]/F");
  g4hitntuple->Branch("extractedtime_emcal",& m_extractedtime,"extractedtime_emcal[24576]/F");
  g4hitntuple->Branch("extractedadc_ihcal",& m_extractedadc_ihcal,"extractedadc_ihcal[1536]/F");
  g4hitntuple->Branch("extractedtime_ihcal",& m_extractedtime_ihcal,"extractedtime_ihcal[1536]/F");
  g4hitntuple->Branch("extractedadc_ohcal",& m_extractedadc_ohcal,"extractedadc_ohcal[1536]/F");
  g4hitntuple->Branch("extractedtime_ohcal",& m_extractedtime_ohcal,"extractedtime_ohcal[1536]/F");
  g4hitntuple->Branch("toweradc_emcal",& m_toweradc,"toweradc_emcal[24576]/F");
  g4hitntuple->Branch("toweradc_ihcal",& m_toweradc_ihcal,"toweradc_ihcal[1536]/F");
  g4hitntuple->Branch("toweradc_ohcal",& m_toweradc_ohcal,"toweradc_ohcal[1536]/F");
  g4hitntuple->Branch("waveform_emcal",& m_waveform,"waveform_emcal[24576][16]/I");
  g4hitntuple->Branch("waveform_ihcal",& m_waveform_ihcal,"waveform_ihcal[1536][16]/I");
  g4hitntuple->Branch("waveform_ohcal",& m_waveform_ohcal,"waveform_ohcal[1536][16]/I");

  // g4hitntuple->Branch("npeaks_ihcal",& m_npeaks_ihcal,"npeaks_ihcall[1536]/I");
  // g4hitntuple->Branch("npeaks_ohcal",& m_npeaks_ohcal,"npeaks_ohcal[1536]/I");

  //----------------------------------------------------------------------------------------------------
  //Read in the template file, this currently points to a tim local area file, 
  //but a copy of this file is in the git repository.
  //----------------------------------------------------------------------------------------------------
  // std::string template_input_file = "/gpfs/mnt/gpfs02/sphenix/user/trinn/fitting_algorithm_playing/prdfcode/prototype/offline/packages/Prototype4/templates.root";
  std::string cemc_template_input_file = "/gpfs/mnt/gpfs02/sphenix/user/trinn/fitting_algorithm_playing/waveform_simulation/calibrations/WaveformProcessing/templates/testbeam_cemc_template.root";
  std::string ihcal_template_input_file = "/gpfs/mnt/gpfs02/sphenix/user/trinn/fitting_algorithm_playing/waveform_simulation/calibrations/WaveformProcessing/templates/testbeam_ihcal_template.root";
  std::string ohcal_template_input_file = "/gpfs/mnt/gpfs02/sphenix/user/trinn/fitting_algorithm_playing/waveform_simulation/calibrations/WaveformProcessing/templates/testbeam_ohcal_template.root";

  TFile* fin1 = TFile::Open(cemc_template_input_file.c_str());
  assert(fin1);
  assert(fin1->IsOpen());
  h_template = static_cast<TProfile*>(fin1->Get("waveform_template"));

  TFile* fin2 = TFile::Open(ihcal_template_input_file.c_str());
  assert(fin2);
  assert(fin2->IsOpen());

  h_template_ihcal = static_cast<TProfile*>(fin2->Get("waveform_template"));

  TFile* fin3 = TFile::Open(ohcal_template_input_file.c_str());
  assert(fin3);
  assert(fin3->IsOpen());

  h_template_ohcal = static_cast<TProfile*>(fin3->Get("waveform_template"));


  WaveformProcessing = new CaloWaveformProcessing();
  WaveformProcessing->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  WaveformProcessing->set_template_file("testbeam_cemc_template.root");
  WaveformProcessing->initialize_processing();



  WaveformProcessing_ihcal = new CaloWaveformProcessing();
  WaveformProcessing_ihcal->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  WaveformProcessing_ihcal->set_template_file("testbeam_ihcal_template.root");
  WaveformProcessing_ihcal->initialize_processing();


  WaveformProcessing_ohcal = new CaloWaveformProcessing();
  WaveformProcessing_ohcal->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  WaveformProcessing_ohcal->set_template_file("testbeam_ohcal_template.root");
  WaveformProcessing_ohcal->initialize_processing();


  light_collection_model.load_data_file(string(getenv("CALIBRATIONROOT")) + string("/CEMC/LightCollection/Prototype3Module.xml"),
								   "data_grid_light_guide_efficiency", "data_grid_fiber_trans");

  return 0;
}

int CaloWaveFormSim::process_event(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;
  GetNodes(topNode);
  process_g4hits(topNode);
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
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
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", detector));
  if (!detNode)
    {
      std::cout << PHWHERE << "Detector Electronics Node missing, making one"<<std::endl;
      detNode = new PHCompositeNode(detector);
      dstNode->addNode(detNode);
    }
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
  //-* contains final physics products (nmips, t0, etc)
  std::string nodename = "WAVEFORMS_" + detector;
  _waveforms= findNode::getClass<WaveformContainerv1>(detNode, nodename.c_str());
  if (!_waveforms)
    {
      _waveforms = new WaveformContainerv1();
      PHIODataNode<PHObject> *waveformcontainerNode = new PHIODataNode<PHObject>(_waveforms, nodename.c_str(), "PHObject");
      detNode->addNode(waveformcontainerNode);
    }

}

void CaloWaveFormSim::GetNodes(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
  ostringstream nodename;
  nodename.str("");
  nodename << "G4HIT_" << detector;

  _layergeo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_CEMC");

  if (!_layergeo)
    {
      std::cout << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - Could not locate sim geometry node "
		<<  std::endl;
      exit(1);
    }
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  _seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_CEMC");
  _hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  _hits_ihcal = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_HCALIN");
  _hits_ohcal = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_HCALOUT");
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
}

int CaloWaveFormSim::process_g4hits(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  const PHG4CylinderGeom *layergeom_raw = _layergeo->GetFirstLayerGeom();
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
  assert(layergeom_raw);

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  const PHG4CylinderGeom_Spacalv3 *layergeom =
      dynamic_cast<const PHG4CylinderGeom_Spacalv3 *>(layergeom_raw);
  assert(layergeom);
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  PHG4CylinderCellGeom *geo_raw = _seggeo->GetFirstLayerCellGeom();
  PHG4CylinderCellGeom_Spacalv1 *geo = dynamic_cast<PHG4CylinderCellGeom_Spacalv1 *>(geo_raw);
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  float waveform[24576][16];
  float waveform_ihcal[1536][16];
  float waveform_ohcal[1536][16];
  float tedep[24576];
  float tedep_ihcal[1536];
  float tedep_ohcal[1536];
  for (int i = 0; i < 24576;i++)
    {
      for (int j = 0; j < 16;j++)
	{
	  m_waveform[i][j] = 0.0;
	  waveform[i][j] = 0.0;
	}     
      m_extractedadc[i] = 0;
      m_extractedtime[i] = 0;
      m_toweradc[i] = 0;
      m_ndep[i] = 0.0;
      m_tedep[i] = 0.0;
      tedep[i] = 0.0;
    }

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  for (int i = 0; i < 1536;i++)
    {
      tedep_ihcal[i] = 0.0;  
      m_extractedadc_ihcal[i] = 0;
      m_extractedtime_ihcal[i] = 0;
      tedep_ohcal[i] = 0.0;  
      m_extractedadc_ohcal[i] = 0;
      m_extractedtime_ohcal[i] = 0; 
      m_toweradc_ihcal[i] = 0.0;
      m_toweradc_ohcal[i] = 0.0;
      for (int j = 0; j < 16;j++)
	{
	  m_waveform_ihcal[i][j] = 0;
	  waveform_ihcal[i][j] = 0;
	  m_waveform_ohcal[i][j] = 0;
	  waveform_ohcal[i][j] = 0;
	}
    }
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  //---------------------------------------------------------
  //Load in the template function as a TF1
  //for use in waveform generation
  //---------------------------------------------------------

  TF1* f_fit = new TF1("f_fit",template_function,0,31,3);
  f_fit->SetParameters(1,0,0);


  TF1* f_fit_ihcal = new TF1("f_fit_ihcal",template_function_ihcal,0,31,3);
  f_fit_ihcal->SetParameters(1,0,0);

  TF1* f_fit_ohcal = new TF1("f_fit_ohcal",template_function_ohcal,0,31,3);
  f_fit_ohcal->SetParameters(1,0,0);

  //-----------------------------------------------------
  //Set the timeing in of the prompt 
  //signal peak to be 4 time samples into
  //the waveform
  //------------------------------------------------------
  float _shiftval = 4-f_fit->GetMaximumX();
  f_fit->SetParameters(1,_shiftval,0);

  float _shiftval_ihcal = 4-f_fit_ihcal->GetMaximumX();
  f_fit_ihcal->SetParameters(1,_shiftval_ihcal,0);

  float _shiftval_ohcal = 4-f_fit_ohcal->GetMaximumX();
  f_fit_ohcal->SetParameters(1,_shiftval_ohcal,0);


  //-----------------------------------------------------------------------------------------
  //Loop over truth primary particles and record their kinematics
  //-----------------------------------------------------------------------------------------
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  if (_hits)
  {
    //-----------------------------------------------------------------------
    //Loop over G4Hits to build a waveform simulation
    //-----------------------------------------------------------------------
    PHG4HitContainer::ConstRange hit_range = _hits->getHits();
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
      if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

      //-------------------------------------------------------------------------
      //Map the G4hits to the corresponding CEMC tower
      //-------------------------------------------------------------------------
      int etabin = etabinshort;
      int phibin = phibin_cell;
      int towernumber = etabin + 96*phibin;
      //---------------------------------------------------------------------
      //Convert the G4hit into a waveform contribution
      //---------------------------------------------------------------------
      // float t0 = 0.5*(hit_iter->second->get_t(0)+hit_iter->second->get_t(1)) / 16.66667;   //Average of g4hit time downscaled by 16.667 ns/time sample 
      float t0 = (hit_iter->second->get_t(0)) / 16.66667;   //Place waveform at the starting time of the G4hit, avoids issues caused by excessively long lived g4hits
      float tmax = 16.667*16;
      float tmin = -20;
      f_fit->SetParameters(light_yield*26000,_shiftval+t0,0);            //Set the waveform template to match the expected signal from such a hit
      if (hit_iter->second->get_t(1) >= tmin && hit_iter->second->get_t(0) <= tmax) {
	tedep[towernumber] += light_yield*26000;    // add g4hit adc deposition to the total deposition  
      }
      //-------------------------------------------------------------------------------------------------------------
      //For each tower add the new waveform contribution to the total waveform
      //-------------------------------------------------------------------------------------------------------------
      if (hit_iter->second->get_edep()*26000 > 1 && hit_iter->second->get_t(1) >= tmin && hit_iter->second->get_t(0) <= tmax)
	{
	  for (int j = 0; j < 16;j++) // 16 is the number of time samples
	    {
	      waveform[towernumber][j] += f_fit->Eval(j);
	    
	    }
	  m_ndep[towernumber] +=1;
	}
      m_phibin.push_back(phibin);
      m_etabin.push_back(etabin);
    }
  }
  
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;


  //--------------
  // do IHCAL
  //--------------

  if (_hits_ihcal)
  {
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
      int towernumber = etabin + 24*phibin;
     
      //---------------------------------------------------------------------
      //Convert the G4hit into a waveform contribution
      //---------------------------------------------------------------------
      // float t0 = 0.5*(hit_iter->second->get_t(0)+hit_iter->second->get_t(1)) / 16.66667;   //Average of g4hit time downscaled by 16.667 ns/time sample 
      float t0 = (hit_iter->second->get_t(0)) / 16.66667;   //Place waveform at the starting time of the G4hit, avoids issues caused by excessively long lived g4hits
      float tmax = 16.667*16;
      float tmin = -20;
      f_fit_ihcal->SetParameters(light_yield*2600,_shiftval_ihcal+t0,0);            //Set the waveform template to match the expected signal from such a hit
      if (hit_iter->second->get_t(1) >= tmin && hit_iter->second->get_t(0) <= tmax) {
	tedep_ihcal[towernumber] += light_yield*2600;    // add g4hit adc deposition to the total deposition
      }
      //-------------------------------------------------------------------------------------------------------------
      //For each tower add the new waveform contribution to the total waveform
       //-------------------------------------------------------------------------------------------------------------
     if (hit_iter->second->get_edep()*2600 > 1 &&hit_iter->second->get_t(1) >= tmin && hit_iter->second->get_t(0) <= tmax )
	{
	  for (int j = 0; j < 16;j++) // 16 is the number of time samples
	    {
	      waveform_ihcal[towernumber][j] +=f_fit_ihcal->Eval(j);
	    }
	}
    }
  }


  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  //--------------
  // do OHCAL
  //--------------

  if (_hits_ohcal)
  {
    //-----------------------------------------------------------------------
    //Loop over G4Hits to build a waveform simulation
    //-----------------------------------------------------------------------
    PHG4HitContainer::ConstRange hit_range = _hits_ohcal->getHits();
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
      //-------------------------------------------------------------------------
      //Map the G4hits to the corresponding CEMC tower
      //-------------------------------------------------------------------------
      
      int etabin = towereta;
      int phibin =towerphi;
      int towernumber = etabin + 24*phibin;
      
      //---------------------------------------------------------------------
      //Convert the G4hit into a waveform contribution
      //---------------------------------------------------------------------
      // float t0 = 0.5*(hit_iter->second->get_t(0)+hit_iter->second->get_t(1)) / 16.66667;   //Average of g4hit time downscaled by 16.667 ns/time sample 
      float t0 = (hit_iter->second->get_t(0)) / 16.66667;   //Place waveform at the starting time of the G4hit, avoids issues caused by excessively long lived g4hits
      float tmax =16.667*16 ;
      float tmin = -20;
      f_fit_ohcal->SetParameters(light_yield*5000,_shiftval_ohcal+t0,0);            //Set the waveform template to match the expected signal from such a hit
      if (hit_iter->second->get_t(0) < tmax && hit_iter->second->get_t(1) > tmin) {
	  {
	    tedep_ohcal[towernumber] += light_yield*5000;    // add g4hit adc deposition to the total deposition 
	  }
      }
      //-------------------------------------------------------------------------------------------------------------
      //For each tower add the new waveform contribution to the total waveform
      //-------------------------------------------------------------------------------------------------------------

     if (hit_iter->second->get_edep()*5000 > 1 && hit_iter->second->get_t(1) >= tmin && hit_iter->second->get_t(0) <= tmax)
	{
	  for (int j = 0; j < 16;j++) // 16 is the number of time samples
	    {
	      waveform_ohcal[towernumber][j] +=f_fit_ohcal->Eval(j);
	    }
	}
    }
  }


  //----------------------------------------------------------------------------------------
  //For each tower loop over add a noise waveform 
  //from cosmic data, gain is too high but is corrected for
  //----------------------------------------------------------------------------------------
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  //-----------------------------
  // do noise for EMCal:
  //-----------------------------
  for (int i = 0; i < 24576;i++)
    {
      int noise_waveform  = (int)rnd->Uniform(0,1990);
      noise_midrad->GetEntry(noise_waveform);
      m_tedep[i] = tedep[i];
      
      std::vector<int> v_waveform;
      for (int k = 0; k < 16;k++)
	{
	  m_waveform[i][k] = waveform[i][k]+(noise_val_midrad[k]-1500)/16.0+1500;
	  v_waveform.push_back(m_waveform[i][k]);
	}
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

      _waveforms->set_waveform_at_eta_phi(i%96, i/96, &v_waveform);
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

    }
  //---------------------------
  // do noise for ihcal:
  //---------------------------
  for (int i = 0; i < 1536;i++)
    {
      int noise_waveform  = (int)rnd->Uniform(0,1990);
      noise_lowrad->GetEntry(noise_waveform);     
      m_tedep_ihcal[i] = tedep_ihcal[i];
      for (int k = 0; k < 16;k++)
	{
	   m_waveform_ihcal[i][k] = waveform_ihcal[i][k]+(noise_val_lowrad[k]-1500)/16.0+1500;
	  // m_waveform_ihcal[i][k] = waveform_ihcal[i][k]+1500;
	}
    }
  //---------------------------
  // do noise for ohcal:
  //---------------------------
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  for (int i = 0; i < 1536;i++)
    {
      int noise_waveform  = (int)rnd->Uniform(0,1990);
      noise_norad->GetEntry(noise_waveform);
      
      m_tedep_ohcal[i] = tedep_ohcal[i];
      for (int k = 0; k < 16;k++)
	{
	  m_waveform_ohcal[i][k] = waveform_ohcal[i][k]+(noise_val_norad[k]-1500)/16.0+1500;
	  // m_waveform_ohcal[i][k] = waveform_ohcal[i][k]+1500;
	}
    }
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
  

  //------------------------------
  //Clear vector content
  //------------------------------
 
  m_primid.clear();
  m_primtrkid.clear();
  m_primpt.clear();
  m_primeta.clear();
  m_primphi.clear();
  m_etabin.clear();
  m_phibin.clear();

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}



int CaloWaveFormSim::End(PHCompositeNode* topNode)
{
  outfile->cd();
  g4hitntuple->Write();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");
  return 0;
}
