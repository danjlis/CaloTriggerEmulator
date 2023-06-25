#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <CaloWaveFormSim.h>
#include <CaloTriggerEmulator.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalowaveformsim.so)
R__LOAD_LIBRARY(libcalotriggeremulator.so)
#endif

void Fun4All_CaloTrigSim(const char *outfile = "trees.root", const char *outfile2 = "trees2.root")
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libcalotriggeremulator");

  Fun4AllServer *se = Fun4AllServer::instance();

  CaloWaveFormSim *ca = new CaloWaveFormSim("CALOWAVEFORMSIM",outfile);
  ca->Detector("HCALOUT");
  se->registerSubsystem(ca);
  
  CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR",outfile2);
  te->TriggerType("OHCAL");
  te->Verbosity(5);
  se->registerSubsystem(te);
  
  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  in->AddListFile("dst_truth.list");

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("in2");
  in2->AddListFile("dst_calo_cluster.list");

  Fun4AllInputManager *in3 = new Fun4AllDstInputManager("in3");
  in3->AddListFile("g4hits.list");

// Fun4All
  se->registerInputManager(in3);
  se->registerInputManager(in2);
  se->registerInputManager(in);
  se->run(2);
  se->End();
}
