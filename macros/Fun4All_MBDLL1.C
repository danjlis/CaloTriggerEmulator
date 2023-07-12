#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <CaloWaveFormSim.h>
#include <LL1PacketGetter.h>

R__LOAD_LIBRARY(libcalowaveformsim.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libll1packetgetter.so)
#endif

void Fun4All_MBDLL1(const std::string &fname1 = "/sphenix/lustre01/sphnxpro/commissioning/LL1/beam/beam_LL1-00014232-0000.prdf", const char *outfile = "trees.root", const char *outfile2 = "trees2.root", const char *outfile3 = "trees3.root")
{
  gSystem->Load("libg4dst");
  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libcalotriggeremulator");
  gSystem->Load("libll1packetgetter");
  
  Fun4AllServer *se = Fun4AllServer::instance();
  
  LL1PacketGetter *ca = new LL1PacketGetter("LL1PACKETGETTER_MBD","MBD", "MBD");
  ca->set_nsamples(20);
  ca->SetVerbosity(1);
  se->registerSubsystem(ca);
    
  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

// Fun4All

  se->run(3);
  se->End();
}
