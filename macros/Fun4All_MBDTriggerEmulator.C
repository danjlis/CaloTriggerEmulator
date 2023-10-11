#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <CaloWaveFormSim.h>
#include <CaloTriggerEmulator.h>
#include <CaloPacketGetter.h>
#include <calowaveformsim/MBDEmulatorTreeMaker.h>

R__LOAD_LIBRARY(libmbdemulatortreemaker.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalowaveformsim.so)
R__LOAD_LIBRARY(libcalotriggeremulator.so)
R__LOAD_LIBRARY(libcalopacketgetter.so)
#endif

  void Fun4All_MBDTriggerEmulator(const int runnumber)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libcalotriggeremulator");
  gSystem->Load("libcalopacketgetter");
  gSystem->Load("libtmbdemulatortreemaker");

  const char *outfile1 = Form("mbdemu_ZDC_000%d.root", runnumber);
  const char *outfile2 = Form("mbdemuhist_ZDC_000%d.root", runnumber);

  std::string fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/mbd/beam/beam_seb18-000%d-0000.prdf", runnumber);;


  Fun4AllServer *se = Fun4AllServer::instance();

  CaloPacketGetter *ca = new CaloPacketGetter("CALOPACKETGETTER_MBD","MBD");
  ca->set_nsamples(31);
  se->registerSubsystem(ca);
  
  CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR_MBD",outfile2);
  te->Verbosity(0);
  te->setTriggerType("MBD");
  se->registerSubsystem(te);

  MBDEmulatorTreeMaker *tt1 = new MBDEmulatorTreeMaker("MBDEMULATORTREEMAKER_MBD",outfile1, "LL1OUT_MBD");
  se->registerSubsystem(tt1);

  
  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

// Fun4All

  se->run(100000);
  se->End();
}
