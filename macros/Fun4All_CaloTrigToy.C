#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <CaloWaveFormToy.h>
#include <CaloTriggerEmulator.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalowaveformsim.so)
R__LOAD_LIBRARY(libcalotriggeremulator.so)
#endif

void Fun4All_CaloTrigToy(const char *outfile = "trees_calo.root", const char *outfile2 = "trees2_calo.root")
{
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libcalotriggeremulator");

  Fun4AllServer *se = Fun4AllServer::instance();

  CaloWaveFormToy *ca = new CaloWaveFormToy("CALOWAVEFORMTOY",outfile);
  ca->Detector("EMCAL");
  ca->Verbosity(2);
  ca->setNSamples(16);
  se->registerSubsystem(ca);

  
  CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR",outfile2);
  te->setTriggerType("JET");
  te->setNSamples(16);
  te->Verbosity(1);
  se->registerSubsystem(te);
  

  se->run(3);
  se->End();
}
