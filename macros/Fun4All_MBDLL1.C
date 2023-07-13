#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <CaloTriggerEmulator.h>
#include <CaloWaveFormSim.h>
#include <LL1PacketGetter.h>
#include <CaloPacketGetter.h>
#include <calowaveformsim/MBDEmulatorTreeMaker.h>

R__LOAD_LIBRARY(libmbdemulatortreemaker.so)
R__LOAD_LIBRARY(libcalotriggeremulator.so)
R__LOAD_LIBRARY(libcalowaveformsim.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libll1packetgetter.so)
R__LOAD_LIBRARY(libcalopacketgetter.so)
#endif

void Fun4All_MBDLL1(const std::string &fname1 = "/sphenix/user/dlis/Projects/raw/mbdll1_00014232-0000.prdf", const char *outfile = "trees.root", const char *outfile2 = "trees2.root", const char *outfile3 = "trees3.root")
{
  gSystem->Load("libtmbdemulatortreemaker");
  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libcalotriggeremulator");
  gSystem->Load("libll1packetgetter");
  gSystem->Load("libcalopacketgetter");
  
  Fun4AllServer *se = Fun4AllServer::instance();
  
  CaloPacketGetter *ca = new CaloPacketGetter("CALOPACKETGETTER_MBD","MBD");
  ca->set_nsamples(31);
  se->registerSubsystem(ca);
  
  LL1PacketGetter *la = new LL1PacketGetter("LL1PACKETGETTER_MBD","MBD", "MBD");
  la->set_nsamples(20);
  la->SetVerbosity(0);
  se->registerSubsystem(la);
  
  CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR_MBD",outfile2);
  te->Verbosity(0);
  te->setTriggerType("MBD");
  se->registerSubsystem(te);
  
  MBDEmulatorTreeMaker *tt1 = new MBDEmulatorTreeMaker("MBDEMULATORTREEMAKER_MBD","trig_tree_mbd.root", "LL1OUT_MBD");
  se->registerSubsystem(tt1);

  MBDEmulatorTreeMaker *tt2 = new MBDEmulatorTreeMaker("MBDEMULATORTREEMAKER_LL1","trig_tree_ll1.root", "LL1OUT_RAW_MBD");
  se->registerSubsystem(tt2);

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

// Fun4All

  se->run(3);
  se->End();
}
