// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MBDPACKETGETTER_H
#define MBDPACKETGETTER_H

#include <fun4all/SubsysReco.h>
#include <calowaveformsim/WaveformContainerv1.h>
#include <climits>
#include <string>


class PHCompositeNode;
class WaveformContainerv1;

class MBDPacketGetter : public SubsysReco
{
 public:
  explicit MBDPacketGetter(const std::string &name = "MBDPacketGetter");
  ~MBDPacketGetter() override;

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

  WaveformContainerv1 *m_WaveformContainer;  //! Calo info
  int m_packet_low;
  int m_packet_high;
  int m_nsamples;
  int m_nchannels;
  bool m_isdata;
};

#endif  // CALOTOWERBUILDER_H
