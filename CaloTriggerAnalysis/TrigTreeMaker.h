#ifndef TRIGTREEMAKER_H
#define TRIGTREEMAKER_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <string>
#include <vector>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdPmtContainer.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include "g4detectors/PHG4CylinderGeomContainer.h"
#include "g4detectors/PHG4CylinderGeom_Spacalv1.h"  // for PHG4CylinderGeom_Spaca...
#include "g4detectors/PHG4CylinderGeom_Spacalv3.h"
#include "g4detectors/PHG4CylinderCellGeomContainer.h"
#include "g4detectors/PHG4CylinderCellGeom_Spacalv1.h"

#include <g4main/PHG4TruthInfoContainer.h>
#include <calowaveformsim/WaveformContainerv1.h>
#include <calowaveformsim/LL1Outv1.h>
#include "TTree.h"
#include "TFile.h"


class PHCompositeNode;

class TrigTreeMaker : public SubsysReco
{
 public:

  TrigTreeMaker(const std::string &name = "TrigTreeMaker", const std::string &outfilename = "trees_trig.root");

  virtual ~TrigTreeMaker();

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
 
  void process_truth();

  int _verbosity;

  TFile *_f;
  TTree *_tree;
  std::string _foutname;

  int _i_event;
  WaveformContainerv1 *_waveforms_cemc;
  WaveformContainerv1 *_waveforms_hcalin;
  WaveformContainerv1 *_waveforms_hcalout;
  WaveformContainerv1 *_waveforms_bbc;
  LL1Outv1 *_ll1_mbd;
  MbdVertexMap *_bbcvertexmap;

  PHG4TruthInfoContainer* _truthinfo;
  PHG4HitContainer* _hits_cemc;
  PHG4CylinderGeomContainer* _layergeo;
  PHG4CylinderCellGeomContainer *_seggeo; 
  PHG4HitContainer* _hits_hcalin;
  PHG4HitContainer* _hits_hcalout;
  MbdPmtContainer* _bbc_pmts;

  int m_waveforms_cemc[24576][16];
  int m_waveforms_hcalin[1536][16];
  int m_waveforms_hcalout[1536][16];
  int m_waveforms_bbc[256][16];

  int m_mbd_trigger_bits;
  int m_mbd_trigger_words[13];

  float m_bbc_vtx_z, m_bbc_vtx_t0;
  int m_bbc_vtx_n;
  std::vector<float> m_primpt;
  std::vector<float> m_primeta;
  std::vector<float> m_primphi;
  
  float m_lightyield_cemc[24576];
  float m_edep_cemc[24576];

  float m_lightyield_hcalout[1536];
  float m_edep_hcalout[1536];

  float m_lightyield_hcalin[1536];
  float m_edep_hcalin[1536];

  int m_npmt;
  float m_adc_bbc[128];
  float m_tdc0_bbc[128];
  float m_tdc1_bbc[128];

  int ROWDIM = 320;
  int COLUMNDIM = 27;

};

#endif 
