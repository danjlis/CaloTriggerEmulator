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


R__LOAD_LIBRARY(libcalotriggeremulator.so)
R__LOAD_LIBRARY(libcalowaveformsim.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libll1packetgetter.so)
R__LOAD_LIBRARY(libcalopacketgetter.so)
#endif

void Fun4All_MBDLL1(const int runnumber)
{ 

  gSystem->Load("libtmbdemulatortreemaker");
  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libcalotriggeremulator");
  gSystem->Load("libll1packetgetter");
  gSystem->Load("libcalopacketgetter");

  const char *outfile1 = Form("mbd_000%d.root", runnumber);
  const char *outfile2 = Form("mbdemu_000%d.root", runnumber);
  const char *outfile3 = Form("ll1_000%d.root", runnumber);

  std::string fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/LL1/beam/beam_LL1-000%d-0000.prdf", runnumber);;
  Fun4AllServer *se = Fun4AllServer::instance();
  
  
  LL1PacketGetter *la = new LL1PacketGetter("LL1PACKETGETTER_MBD","MBD", "MBD");
  la->set_nsamples(20);
  la->SetVerbosity(0);
  se->registerSubsystem(la);

  MBDEmulatorTreeMaker *tt2 = new MBDEmulatorTreeMaker("MBDEMULATORTREEMAKER_LL1",outfile3, "LL1OUT_RAW_MBD");
  se->registerSubsystem(tt2);

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

// Fun4All

  se->run(100000);
  se->End();
}
