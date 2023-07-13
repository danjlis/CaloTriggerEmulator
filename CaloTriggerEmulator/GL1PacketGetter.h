// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GL1PACKETGETTER_H
#define GL1PACKETGETTER_H

#include <fun4all/SubsysReco.h>
#include <climits>
#include <string>
#include <TTree.h>
#include <TFile.h>

class PHCompositeNode;
class TFile;
class TTree;

class GL1PacketGetter : public SubsysReco
{
 public:

  explicit GL1PacketGetter(const std::string &name = "GL1PacketGetter", const std::string &outfilename = "gl1.root");
  ~GL1PacketGetter() override;
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void SetVerbosity(int v){_verbose = v;}
  int End(PHCompositeNode *topNode) override;
 private:
  std::string _foutname;

  TFile *_f;
  TTree *_t;

  unsigned  long long int b_packet;
  unsigned  long long int b_bco;
  
  unsigned  long long int b_raw[64];
  unsigned  long long int b_live[64];
  unsigned  long long int b_scaled[64];

  int m_packet_low;
  int m_packet_high;

  int m_ntriggers;
  int _verbose;

};

#endif  // GL1TOWERBUILDER_H
