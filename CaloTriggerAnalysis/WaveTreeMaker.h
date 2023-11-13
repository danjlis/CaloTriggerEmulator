#ifndef TRIGTREEMAKER_H
#define TRIGTREEMAKER_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <string>
#include <vector>

#include <calowaveformsim/WaveformContainerv1.h>

#include "TTree.h"
#include "TFile.h"


class PHCompositeNode;

class WaveTreeMaker : public SubsysReco
{
 public:

  WaveTreeMaker(const std::string &name = "WaveTreeMaker", const std::string &outfilename = "trees_trig.root");

  virtual ~WaveTreeMaker();

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;
  
  void GetNodes (PHCompositeNode *topNode);

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void SetVerbosity(int verbo) ;
 private:
  void reset_tree_vars();

  void process_waveforms();
 
  int _verbosity;

  TFile *_f;
  TTree *_tree;
  std::string _foutname;

  int _i_event;
  WaveformContainerv1 *_waveforms_cemc;
  WaveformContainerv1 *_waveforms_hcalin;
  WaveformContainerv1 *_waveforms_hcalout;
  WaveformContainerv1 *_waveforms_bbc;
  WaveformContainerv1 *_waveforms_zdc;


  unsigned int m_evtnr;
  unsigned int m_clk;
  unsigned int m_femevtnr;
  unsigned int m_femclk;

  int m_waveforms_cemc[24576][31];
  int m_waveforms_hcalin[1536][31];
  int m_waveforms_hcalout[1536][31];
  int m_waveforms_bbc[256][31];
  int m_waveforms_zdc[16][31];

  int ROWDIM = 320;
  int COLUMNDIM = 27;

};

#endif 
