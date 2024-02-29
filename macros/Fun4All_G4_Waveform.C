
#include <GlobalVariables.C>
#include <CaloTriggerEmulator.h>
#include <CaloEmulatorTreeMaker.h>
#include <G4_CEmc_Spacal.C>
#include <G4_HcalIn_ref.C>
#include <G4_HcalOut_ref.C>
#include <G4_Input.C>
#include <G4_Production.C>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <caloreco/CaloGeomMapping.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloWaveformProcessing.h>

#include <calowaveformsim/CaloWaveformSim.h>

#include <phool/recoConsts.h>
#include <DLUtility.h>

#include <frog/FROG.h>

R__LOAD_LIBRARY(libFROG.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libdlutility.so)
R__LOAD_LIBRARY(libffamodules.so)
void Fun4All_G4_Waveform(
			 const int nEvents = 10,
			   const int keep_it = -1,
			 const string &filename_1 = "/sphenix/user/dlis/Projects/CaloTriggerEmulator/macros/dst_truth.list",
			 const string &filename_2 = "/sphenix/user/dlis/Projects/CaloTriggerEmulator/macros/g4hits.list",
			 const string &filename_3 = "/sphenix/user/dlis/Projects/CaloTriggerEmulator/macros/dst_calo_cluster.list",
    
			 const string &outputFile = "DST_CALO_WAVEFORM.root",
			 const string &outdir = ".",
			   const string &cdbtag = "MDC2")

{

  std::ostringstream ostr;
  if (keep_it == -1)
    {
      ostr << std::setw(4) << std::setfill('0') << 0;
    }
  else {
    ostr << std::setw(4) << std::setfill('0') << keep_it;
  }  
  
  const char* env_p = std::getenv("CALO_TRIGGER_EMULATOR");

  if(!env_p)
    {
      std::cout << "no CALO_TRIGGER_EMULATOR set."<<endl;
      return;
    }

  const char* env_r = std::getenv("PP_RUNTYPE");

  if(!env_r)
    {
      std::cout << "no PP_RUNTYPE set"<<endl;
      return;
    }

  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << stoi(env_r);

  Fun4AllServer *se = Fun4AllServer::instance();

  recoConsts *rc = recoConsts::instance();

  //===============
  // conditions DB flags
  //===============
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", cdbtag);
  rc->set_uint64Flag("TIMESTAMP", CDB::timestamp);
  CDBInterface::instance()->Verbosity(1);


  // Fun4AllInputManager *in_1 = new Fun4AllDstInputManager("DSTin1");
  // if (keep_it == -1) in_1->AddListFile(filename_1,1);
  // else AddFileFromList(keep_it, filename_1, in_1);

  // se->registerInputManager( in_1 );
  
  Fun4AllInputManager *in_2 = new Fun4AllDstInputManager("DSTin2");
  if (keep_it == -1) in_2->AddListFile(filename_2,1);
  else AddFileFromList(keep_it, filename_2, in_2);

  se->registerInputManager( in_2 );

  Fun4AllInputManager *in_3 = new Fun4AllDstInputManager("DSTin3");
  if (keep_it == -1) in_3->AddListFile(filename_3,1);
  else AddFileFromList(keep_it, filename_3, in_3);
  se->registerInputManager( in_3 );

  // Input options
  //===============
  // verbosity setting (applies to all input managers)

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  // set up production relatedstuff
  Enable::PRODUCTION = true;

  //======================
  // Write the DST
  //======================

  Enable::DSTOUT = true;
  Enable::DSTOUT_COMPRESS = false;
  DstOut::OutputDir = outdir;
  DstOut::OutputFile = outputFile;

  //======================
  // What to run
  //======================

  

  CaloWaveformSim *caloWaveformSim = new CaloWaveformSim();
  caloWaveformSim->set_detector_type(CaloTowerDefs::HCALOUT);
  caloWaveformSim->set_detector("HCALOUT");
  caloWaveformSim->set_nsamples(12);
  caloWaveformSim->set_noise_type(CaloWaveformSim::NOISE_NONE);
  se->registerSubsystem(caloWaveformSim);
  

  caloWaveformSim = new CaloWaveformSim();
  caloWaveformSim->set_detector_type(CaloTowerDefs::HCALIN);
  caloWaveformSim->set_detector("HCALIN");
  caloWaveformSim->set_nsamples(12);
  caloWaveformSim->set_noise_type(CaloWaveformSim::NOISE_NONE);
  se->registerSubsystem(caloWaveformSim);


 

  caloWaveformSim = new CaloWaveformSim();
  caloWaveformSim->set_detector_type(CaloTowerDefs::CEMC);
  caloWaveformSim->set_detector("CEMC");
  caloWaveformSim->set_nsamples(12);
  caloWaveformSim->set_calibName("cemc_pi0_twrSlope_v1_default");
  
  caloWaveformSim->set_noise_type(CaloWaveformSim::NOISE_NONE);
  
  caloWaveformSim->get_light_collection_model().load_data_file(
							       string(getenv("CALIBRATIONROOT")) +
							       string("/CEMC/LightCollection/Prototype3Module.xml"),
							       "data_grid_light_guide_efficiency", "data_grid_fiber_trans");
  
  se->registerSubsystem(caloWaveformSim);
          
  
 
  
  

  CaloTowerBuilder *ca2 = new CaloTowerBuilder();
  ca2->set_detector_type(CaloTowerDefs::HCALOUT);
  ca2->set_nsamples(12);
  ca2->set_dataflag(false);
  ca2->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  ca2->set_builder_type(CaloTowerDefs::kWaveformTowerv2);
  se->registerSubsystem(ca2);

  ca2 = new CaloTowerBuilder();
  ca2->set_detector_type(CaloTowerDefs::HCALIN);
  ca2->set_nsamples(12);
  ca2->set_dataflag(false);
  ca2->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  ca2->set_builder_type(CaloTowerDefs::kWaveformTowerv2);
  se->registerSubsystem(ca2);

  ca2 = new CaloTowerBuilder();
  ca2->set_detector_type(CaloTowerDefs::CEMC);
  ca2->set_nsamples(12);
  ca2->set_dataflag(false);
  ca2->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  ca2->set_builder_type(CaloTowerDefs::kWaveformTowerv2);
  se->registerSubsystem(ca2);

  

  // tower calib
  
  CaloTowerCalib *calib = new CaloTowerCalib();
  calib->set_detector_type(CaloTowerDefs::HCALOUT);
  calib->set_outputNodePrefix("TOWERSWAVEFORM_CALIB_");
  se->registerSubsystem(calib);

  calib = new CaloTowerCalib();
  calib->set_detector_type(CaloTowerDefs::HCALIN);
  calib->set_outputNodePrefix("TOWERSWAVEFORM_CALIB_");
  se->registerSubsystem(calib);

  calib = new CaloTowerCalib();
  calib->set_detector_type(CaloTowerDefs::CEMC);
  calib->set_outputNodePrefix("TOWERSWAVEFORM_CALIB_");
  se->registerSubsystem(calib);
  

  if (Enable::DSTOUT)
    {
      char *dir = new char[100];
      if (keep_it == -1)
	{
	  sprintf(dir, "macros");
	}
      else
	{
	  sprintf(dir, "outputs/waveform");
	}

      string FullOutFile = Form("%s/%s/DST_WAVEFORM_%s_%s.root", env_p, dir, rstr.str().c_str(), ostr.str().c_str());

      Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
      out->AddNode("Sync");
      out->AddNode("EventHeader");
      // Inner Hcal
      out->AddNode("WAVEFORM_HCALIN");
      out->AddNode("TOWERSWAVEFORM_CALIB_HCALIN");
      out->AddNode("TOWERINFO_SIM_HCALIN");
      out->AddNode("TOWERINFO_CALIB_HCALIN");

      // Outer Hcal
      out->AddNode("WAVEFORM_HCALOUT");
      out->AddNode("TOWERSWAVEFORM_CALIB_HCALOUT");
      out->AddNode("TOWERINFO_SIM_HCALOUT");
      out->AddNode("TOWERINFO_CALIB_HCALOUT");

      // CEmc
      out->AddNode("WAVEFORM_CEMC");
      out->AddNode("TOWERSWAVEFORM_CALIB_CEMC");
      out->AddNode("TOWERINFO_SIM_CEMC");
      out->AddNode("TOWERINFO_CALIB_CEMC");

      se->registerOutputManager(out);
    }

  //-----------------
  // Event processing
  //-----------------
  // if we use a negative number of events we go back to the command line here
  if (nEvents < 0)
    {
      return 0;
    }
  se->run(nEvents);

  //-----
  // Exit
  //-----

  CDBInterface::instance()->Print(); // print used DB files
  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;
  delete se;
  if (Enable::PRODUCTION)
    {
      Production_MoveOutput();
    }

  gSystem->Exit(0);
}


