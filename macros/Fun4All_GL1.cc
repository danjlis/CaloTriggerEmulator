#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <GL1PacketGetter.h>

R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)

R__LOAD_LIBRARY(libgl1packetgetter.so)
#endif

  void Fun4All_GL1(const std::string &fname1 = "/sphenix/lustre01/sphnxpro/commissioning/GL1/beam/GL1_beam_gl1daq-00021372-0000.prdf", const char *outfile = "gl1_21372.root", const char *outfile2 = "trees2.root", const char *outfile3 = "trees3.root")
{
  gSystem->Load("libg4dst");
  gSystem->Load("libgl1packetgetter");

  Fun4AllServer *se = Fun4AllServer::instance();

  GL1PacketGetter *ca = new GL1PacketGetter("GL1PACKETGETTER",outfile);
  se->registerSubsystem(ca);

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

// Fun4All

  se->run(100000);
  se->End();
}
