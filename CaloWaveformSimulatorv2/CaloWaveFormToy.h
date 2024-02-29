#ifndef CALOWAVEFORMTOY_H__
#define CALOWAVEFORMTOY_H__
#include <fun4all/SubsysReco.h>
#include <TTree.h>
#include <TProfile.h>
#include <TRandom3.h>
#include "WaveformContainerv1.h"

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TProfile;
class TFile;
class TNtuple;
class TTree;

class CaloWaveFormToy : public SubsysReco
{
 public:
  //! constructor
  CaloWaveFormToy(const std::string &name = "CaloWaveFormToy", const std::string &fname = "MyNtuple.root");

  //! destructor
  virtual ~CaloWaveFormToy();


  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
  void CreateNodes(PHCompositeNode*);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  void Detector(const std::string &name) { detector = name; }
  void setNSamples(int nsamples) { _nsamples = nsamples; }

  void setToyType(int type) { _type = type; }
  void setOccupancy(float occ){ _occupancy = occ; }

 private:
  std::string detector;
  std::vector<float> m_waveform_emcal[24576];
  std::vector<float> m_waveform_ihcal[1536];
  std::vector<float> m_waveform_ohcal[1536];
  std::vector<float> m_waveform_mbd[256];
  
  static TProfile* h_template_mbd;
  static TProfile* h_template_emcal;
  static TProfile* h_template_ihcal;
  static TProfile* h_template_ohcal;


  TRandom3* rnd;
  int _verbose = 1;
  int _nsamples = 16;

  enum TOY {
    SQUARE = 0,
    PULSE = 1
  };


  int _type = 0;
  float _occupancy = 1.0;
  int ROWDIM = 320;

  int COLUMNDIM = 27;
 
  static double template_function_mbd(double *x, double *par);
  static double template_function_emcal(double *x, double *par);
  static double template_function_ihcal(double *x, double *par);
  static double template_function_ohcal(double *x, double *par);

  bool IsDetector(const std::string &det);

  WaveformContainerv1 *waveforms_emcal;
  WaveformContainerv1 *waveforms_hcalout;
  WaveformContainerv1 *waveforms_hcalin;
  WaveformContainerv1 *waveforms_mbd;

};

#endif
