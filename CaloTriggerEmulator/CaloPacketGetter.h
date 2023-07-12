// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOPACKETGETTER_H
#define CALOPACKETGETTER_H

#include <fun4all/SubsysReco.h>
#include <calowaveformsim/WaveformContainerv1.h>
#include <climits>
#include <string>


class PHCompositeNode;
class WaveformContainerv1;

class CaloPacketGetter : public SubsysReco
{
 public:
  explicit CaloPacketGetter(const std::string &name = "CaloPacketGetter",const std::string &detector = "HCALOUT");
  ~CaloPacketGetter() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void CreateNodeTree(PHCompositeNode *topNode);

  
  enum DetectorSystem
  {
    CEMC = 0,
    HCALIN = 1,
    HCALOUT = 2,
    SEPD = 3,
    MBD = 4
  };

  int m_detector_type;
  std::map<std::string, int> m_detectors;
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

  WaveformContainerv1 *m_WaveformContainer;  //! Calo info
  int m_packet_low;
  int m_packet_high;
  int m_nsamples;
  int m_nchannels;
  bool m_isdata;
};

#endif  // CALOTOWERBUILDER_H
