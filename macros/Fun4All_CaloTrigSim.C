#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <g4waveform/CaloWaveformSim.h>
#include <CaloTriggerEmulator.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libcalotriggeremulator.so)
#endif

void Fun4All_CaloTrigSim(const char *outfile = "trees.root", const char *outfile2 = "trees2.root")
{
  gSystem->Load("libg4dst");
  gSystem->Load("libCaloWaveformSim");
  gSystem->Load("libcalotriggeremulator");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(5);

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

  CaloWaveformSim *ca = new CaloWaveformSim("CALOWAVEFORMSIM");
  ca->set_nsamples(16);
  ca->set_detector("EMCAL");
  se->registerSubsystem(ca);
  
  // CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR",outfile2);
  // te->setTriggerType("JET");
  // te->setNSamples(16);
  // te->Verbosity(2);
  // se->registerSubsystem(te);
  
  se->run(10);
  se->End();
}
