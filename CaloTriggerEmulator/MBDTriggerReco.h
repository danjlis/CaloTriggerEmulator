// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MBDTRIGGERRECO_H
#define MBDTRIGGERRECO_H

#include <fun4all/SubsysReco.h>
#include <calowaveformsim/WaveformContainerv1.h>
#include <climits>
#include <string>


class PHCompositeNode;
class WaveformContainerv1;

class MBDTriggerReco : public SubsysReco
{
 public:
  explicit MBDTriggerReco(const std::string &name = "MBDTriggerReco");
  ~MBDTriggerReco() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void CreateNodeTree(PHCompositeNode *topNode);

  enum DetectorSystem
  {
    CEMC = 0,
    HCALIN = 1,
    HCALOUT = 2,
    EPD = 3
  };

  void set_nsamples(int _nsamples)
  {
    m_nsamples = _nsamples;
    return;
  }
  void set_dataflag(bool flag)
  {
    m_isdata = flag;
    return;
  }

 private:

  LL1Outv1 *m_LL1Out;
  std::string m_detector; 
  int m_packet_low;
  int m_packet_high;
  int m_nsamples;
  int m_nchannels;
  int m_nadc;
  bool m_isdata;
};

#endif  // CALOTOWERBUILDER_H
