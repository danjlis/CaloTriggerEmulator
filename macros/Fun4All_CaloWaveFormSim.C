#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <calowaveformsim/CaloWaveformSim.h>
#include <DLUtility.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libtrigtreemaker.so)
R__LOAD_LIBRARY(libdlutility.so)

#endif

void Fun4All_CaloWaveFormSim( const std::string &config_file = "waveform_config.config", const char *outfile = "trees.root")
{

  TEnv *config_p = new TEnv(config_file.c_str());
  const int nEvents = config_p->GetValue("NEVENTS",1);
  const string outputFile = config_p->GetValue("OUTFILENAME","testy.root");
  const string outputDir = config_p->GetValue("OUTPUTDIR", "./output/");
  const int verbo = config_p->GetValue("VERBOSITY", 0);
  const int keep_it = config_p->GetValue("KEEPIT", -1);

  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libtrigtreemaker");

  Fun4AllServer *se = Fun4AllServer::instance();

  CaloWaveFormSim *ca = new CaloWaveFormSim("CALOWAVEFORMSIM",outfile);
  ca->Detector("ALL");
  se->registerSubsystem(ca);

  TrigTreeMaker *tt = new TrigTreeMaker("TRIGTREEMAKER", "trigtree.root");
  se->registerSubsystem(tt);


  const std::string filename_1 = "g4hits.list";
  const std::string filename_2 = "dst_calo_cluster.list";
  const std::string filename_3 = "dst_truth.list";

  Fun4AllInputManager *in1 = new Fun4AllDstInputManager("in1");
  if (keep_it == -1) in1->AddListFile(filename_1,1);
  else AddFileFromList(keep_it, filename_1, in1);  
  
  //  in1->AddListFile("dst_truth.list");
  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("in2");
  if (keep_it == -1) in2->AddListFile(filename_2,1);
  else AddFileFromList(keep_it, filename_2, in2);

  //in2->AddListFile("dst_calo_cluster.list");

  Fun4AllInputManager *in3 = new Fun4AllDstInputManager("in3");
  if (keep_it == -1) in3->AddListFile(filename_3,1);
  else AddFileFromList(keep_it, filename_3, in3);
  //in3->AddListFile("g4hits.list");

// Fun4All
  se->registerInputManager(in1);
  se->registerInputManager(in3);
  se->registerInputManager(in2);
  se->run( nEvents );
  se->End();
  std::cout << "ALL Done "<<std::endl;
  delete se;

  return 0;
}
