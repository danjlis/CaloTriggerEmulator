#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>

#include <calowaveformsim/CaloWaveFormSim.h>
#include <g4bbc/BbcDigitization.h>
#include <bbc/BbcReconstruction.h>
#include <calowaveformsim/TrigTreeMaker.h>
#include <MBDTriggerEmulator.h>
#include <DLUtility.h>

R__LOAD_LIBRARY(libtrigtreemaker.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalowaveformsim.so)
R__LOAD_LIBRARY(libdlutility.so)
R__LOAD_LIBRARY(libmbdtriggeremulator.so)
R__LOAD_LIBRARY(libg4bbc.so)

#endif

void Fun4All_MBDTrigSim( const std::string &config_file = "inputs/MBD_config.config")
{

  TEnv *config_p = new TEnv(config_file.c_str());
  const int nEvents = config_p->GetValue("NEVENTS",1);
  const string outfile1 = config_p->GetValue("OUTFILENAME1","testy.root");
  const string outfile2 = config_p->GetValue("OUTFILENAME2","testy2.root");
  const string outfile3 = config_p->GetValue("OUTFILENAME3","testy3.root");
  const string outputDir = config_p->GetValue("OUTPUTDIR", "./outputs/");
  const int verbosity = config_p->GetValue("VERBOSITY", 0);
  const int keep_it = config_p->GetValue("KEEPIT", -1);

  std::cout << "Starting my fun4all script." <<endl;
  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libg4bbc");
  gSystem->Load("libtrigtreemaker");
  gSystem->Load("libdlutility");
  gSystem->Load("libmbdtriggeremulator");

  Fun4AllServer *se = Fun4AllServer::instance();

  auto bbcdigi = new BbcDigitization();
  bbcdigi->Verbosity(verbosity);
  se->registerSubsystem(bbcdigi);

  auto bbcreco = new BbcReconstruction();
  bbcreco->Verbosity(verbosity);
  se->registerSubsystem(bbcreco);

  const string outpath1 = outputDir + outfile1;
  const string outpath2 = outputDir + outfile2;
  const string outpath3 = outputDir + outfile3;

  CaloWaveFormSim *ca = new CaloWaveFormSim("CALOWAVEFORMSIM",outpath1);
  ca->Detector("ALL");
  se->registerSubsystem(ca);

  MBDTriggerEmulator *te = new MBDTriggerEmulator("MBDTRIGGEREMULATOR",outpath2);
  te->TriggerType("MBD");
  te->Verbosity(0);
  se->registerSubsystem(te);

  TrigTreeMaker *tt = new TrigTreeMaker("TRIGTREEMAKER", outpath3);
  se->registerSubsystem(tt);  

  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  if (keep_it == -1) in->AddListFile("dst_truth.list");
  else AddFileFromList(keep_it, "dst_truth.list", in);

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("in2");
  if (keep_it == -1) in2->AddListFile("dst_calo_cluster.list");
  else AddFileFromList(keep_it, "dst_calo_cluster.list", in2);

  Fun4AllInputManager *in3 = new Fun4AllDstInputManager("in3");
  if (keep_it == -1) in3->AddListFile("g4hits.list");
  else AddFileFromList(keep_it, "g4hits.list", in3);

// Fun4All
  se->registerInputManager(in3);
  se->registerInputManager(in2);
  se->registerInputManager(in);
  se->run(nEvents);
  se->End();
}
