
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
#include <CaloTriggerEmulator.h>
#include <TriggerDefs.h>
#include <CaloEmulatorTreeMaker.h>
#include <calowaveformsim/CaloWaveformSim.h>

#include <phool/recoConsts.h>
#include <DLUtility.h>

#include <frog/FROG.h>

R__LOAD_LIBRARY(libFROG.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libcalotriggeremulator.so)
R__LOAD_LIBRARY(libcaloemulatortreemaker.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libdlutility.so)
R__LOAD_LIBRARY(libffamodules.so)
void Fun4All_G4_CaloTrigSim(
			 const int nEvents = 10,
			 const int keep_it = -1,
			 const string &dirname = "/sphenix/user/dlis/Projects/CaloTriggerEmulator/outputs/waveform/",
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

  const char* env_c = std::getenv("CALOTRIGGERTYPE");

  if(!env_c)
    {
      std::cout << "no CALOTRIGGERTYPE set"<<endl;
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

  string filename_1 = "/sphenix/user/dlis/Projects/CaloTriggerEmulator/macros/g4hits.list";
  Fun4AllInputManager *in_1 = new Fun4AllDstInputManager("DSTin1");
  if (keep_it == -1) in_1->AddListFile(filename_1,1);
  else AddFileFromList(keep_it, filename_1, in_1);

  se->registerInputManager( in_1 );
  
  Fun4AllInputManager *in_2 = new Fun4AllDstInputManager("DSTin2");
  in_2->AddFile(Form("%s/DST_WAVEFORM_%s_%s.root",dirname.c_str(), rstr.str().c_str(), ostr.str().c_str()));
  se->registerInputManager( in_2 );
  
  // Input options
  //===============
  // verbosity setting (applies to all input managers)
  
  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);
  
  CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR",Form("%s/outputs/calotrigsim/hist_%s_%s_%s.root", env_p, rstr.str().c_str(), env_c, ostr.str().c_str()));
  te->setTriggerType(env_c);
  te->setEmcalScale(0.005);
  te->setNSamples(12);
  //  te->Verbosity(3);
  se->registerSubsystem(te);


  
  CaloEmulatorTreeMaker *ce = new CaloEmulatorTreeMaker("CALOEMULATORTREEMAKER", Form("%s/outputs/calotrigsim/tree_%s_%s_%s.root", env_p, rstr.str().c_str(), env_c, ostr.str().c_str()));

  ce->SetTriggerType(env_c);
  se->registerSubsystem(ce);
    
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

  gSystem->Exit(0);
}


