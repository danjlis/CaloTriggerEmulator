#ifndef CALOWAVEFORMSIM_H__
#define CALOWAVEFORMSIM_H__
#include "g4detectors/PHG4CylinderCellGeomContainer.h"
#include <bbc/BbcPmtContainer.h>
#include <bbc/BbcPmtHit.h>
#include <fun4all/SubsysReco.h>
#include <TTree.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4FullProjSpacalCellReco.h>
#include "WaveformContainerv1.h"
#include <g4main/PHG4HitContainer.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TProfile;
class TFile;
class TNtuple;
class TTree;

class CaloWaveFormSim : public SubsysReco
{
 public:
  //! constructor
  CaloWaveFormSim(const std::string &name = "CaloWaveFormSim", const std::string &fname = "MyNtuple.root");

  //! destructor
  virtual ~CaloWaveFormSim();


  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
  void CreateNodes(PHCompositeNode*);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  int process_g4hits(PHCompositeNode *);
  int process_g4cells(PHCompositeNode *);
  int process_towers(PHCompositeNode *);
  int process_clusters(PHCompositeNode *);

  void Detector(const std::string &name) { detector = name; }
  void SetNsamples(int nsamples) { _nsamples = nsamples; }
  void SetNoise(int noiselevel) { _noiselevel = noiselevel; }
  void SetGain(std::string gain){ _gain = _gain_opts[gain];}

 private:
  std::string detector;
  std::vector<float> m_waveform_cemc[24576];
  std::vector<float> m_waveform_ihcal[1536];
  std::vector<float> m_waveform_ohcal[1536];
  std::vector<float> m_waveform_bbc[256];
  
  static TProfile* h_template_bbc;
  static TProfile* h_template_emcal;
  static TProfile* h_template_ihcal;
  static TProfile* h_template_ohcal;


  LightCollectionModel light_collection_model;

  int _noiselevel = 0;
  TTree* noise;

  float noise_val[31];

  TRandom3* rnd;
  int _verbose = 1;
  int _nsamples = 16;

  enum GAIN {
    LOW = 0,
    HIGH = 1
  };

  int _gain= GAIN::LOW;



  float _hcalout_lightyield_to_ADC = 0.001;
  
  std::map<std::string, int>  _gain_opts;

  int ROWDIM = 320;
  int COLUMNDIM = 27;
 
  static double template_function_bbc(double *x, double *par);
  static double template_function_cemc(double *x, double *par);
  static double template_function_ihcal(double *x, double *par);
  static double template_function_ohcal(double *x, double *par);

  bool IsDetector(const std::string &det);

  WaveformContainerv1 *waveforms_cemc;
  WaveformContainerv1 *waveforms_hcalout;
  WaveformContainerv1 *waveforms_hcalin;
  WaveformContainerv1 *waveforms_bbc;
  PHG4CylinderGeomContainer *_layergeo;
  PHG4CylinderCellGeomContainer *_seggeo;
  PHG4HitContainer *_hits_cemc;
  PHG4HitContainer *_hits_ihcal;
  PHG4CellContainer *_slats_ohcal;
  TowerInfoContainerv1 *_raw_towers_ohcal;

  BbcPmtContainer *_bbcpmts;
};

#endif
