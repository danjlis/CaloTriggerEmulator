#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <CaloWaveFormSim.h>
#include <CaloTriggerEmulator.h>
#include <CaloPacketGetter.h>

R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalowaveformsim.so)
R__LOAD_LIBRARY(libcalotriggeremulator.so)
R__LOAD_LIBRARY(libcalopacketgetter.so)
#endif

  void Fun4All_CaloTrigEmulator(const std::string &fname1 = "/sphenix/user/dlis/Projects/raw/beam_Hcal-00020096-0000.prdf", const char *outfile = "trees.root", const char *outfile2 = "trees2.root", const char *outfile3 = "trees3.root")
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libcalotriggeremulator");
  gSystem->Load("libcalopacketgetter");

  Fun4AllServer *se = Fun4AllServer::instance();

  CaloPacketGetter *ca = new CaloPacketGetter("CALOPACKETGETTER_HCALOUT","HCALOUT");
  ca->set_nsamples(31);
  se->registerSubsystem(ca);

  CaloPacketGetter *ca1 = new CaloPacketGetter("CALOPACKETGETTER_HCALIN","HCALIN");
  ca1->set_nsamples(31);
  se->registerSubsystem(ca1);
  
  CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR_COSMIC",outfile2);
  te->Verbosity(2);
  te->setTriggerType("COSMIC");
  se->registerSubsystem(te);
  
  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

// Fun4All

  se->run(3);
  se->End();
}
