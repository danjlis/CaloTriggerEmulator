#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloWaveformProcessing.h>

#include <CaloWaveFormSim.h>
#include <CaloTriggerEmulator.h>
#include <CaloPacketGetter.h>
#include <HCALEmulatorTreeMaker.h>
#include <CaloAna.h>

R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalowaveformsim.so)
R__LOAD_LIBRARY(libcalotriggeremulator.so)
R__LOAD_LIBRARY(libcalopacketgetter.so)
R__LOAD_LIBRARY(libcaloana.so)
R__LOAD_LIBRARY(libhcalemulatortreemaker.so)

#endif

  void Fun4All_CosmicTriggerEmulator(const std::string &fname1 = "/sphenix/user/dlis/Projects/raw/cosmics-00025811-0000.prdf", const char *outfile = "trees_25811.root", const char *outfile2 = "trees2_25811.root", const char *outfile3 = "trees3_25811.root")
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libcaloana");
  gSystem->Load("libcalotriggeremulator");
  gSystem->Load("libcalopacketgetter");
  gSystem->Load("libhcalemulatortreemaker");

  Fun4AllServer *se = Fun4AllServer::instance();

  // CaloPacketGetter *ca = new CaloPacketGetter("CALOPACKETGETTER_HCALOUT","HCALOUT");
  // ca->set_nsamples(12);
  // se->registerSubsystem(ca);

  // CaloPacketGetter *ca1 = new CaloPacketGetter("CALOPACKETGETTER_HCALIN","HCALIN");
  // ca1->set_nsamples(12);
  // se->registerSubsystem(ca1);

  CaloTowerBuilder *ctb1 = new CaloTowerBuilder();
  ctb1->set_detector_type(CaloTowerBuilder::HCALIN);
  ctb1->set_nsamples(12);
  ctb1->set_processing_type(CaloWaveformProcessing::FAST);
  se->registerSubsystem(ctb1);

  CaloTowerBuilder *ctb2 = new CaloTowerBuilder();
  ctb2->set_detector_type(CaloTowerBuilder::HCALOUT);
  ctb2->set_nsamples(12);
  ctb2->set_processing_type(CaloWaveformProcessing::FAST);
  se->registerSubsystem(ctb2);

  CaloAna* caloana = new CaloAna("ana","hcalout_25811.root");
  se->registerSubsystem(caloana);

  // CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR_COSMIC",outfile2);
  // te->Verbosity(1);
  // te->setTriggerType("COSMIC_COIN");
  // te->setNSamples(12);
  // te->setThreshold(10);
  // se->registerSubsystem(te);
  
  // HCALEmulatorTreeMaker *tt1 = new HCALEmulatorTreeMaker("HCALEMULATORTREEMAKER","trig_tree_hcal.root", "LL1OUT_COSMIC_COIN");
  // tt1->UseCaloTowerBuilder(true);
  // se->registerSubsystem(tt1);

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

// Fun4All

  se->run();
  se->End();
}
