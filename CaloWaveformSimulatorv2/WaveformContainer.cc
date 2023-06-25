
#include "WaveformReturnCodes.h"
#include "WaveformContainer.h"

#include <cmath>
#include <iostream>

ClassImp(WaveformContainer)

WaveformContainer::WaveformContainer()
{
}
WaveformContainer::~WaveformContainer()
{

}

//______________________________________
void WaveformContainer::Reset()
{
}


void WaveformContainer::set_waveform_at_eta_phi(int eta, int phi, std::vector<int>*)
{

}

void WaveformContainer::set_waveform_at_channel(int index , std::vector<int>*)
{

}

//______________________________________
void WaveformContainer::identify(std::ostream& out)
{
  out << "identify yourself: I am a WaveformContainer object" << std::endl;

}

int WaveformContainer::isValid()
{
  return (!_waveforms.empty());
}
