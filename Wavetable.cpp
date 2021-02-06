#include "Wavetable.h"

template class WavetableOscillator<double>;

WavetableOscillator<double>::WavetableOscillator(double startPhase, double startFreq)
  : iplug::IOscillator<double>(startPhase* mTableSizeM1, startFreq), mPrevFreq(static_cast<int>(startFreq))
{
  SetWavetable("Hydrogen");
}

WavetableOscillator<double>::WavetableOscillator(std::vector<std::string> tablePaths, double startPhase, double startFreq)
  : iplug::IOscillator<double>(startPhase* mTableSizeM1, startFreq), mPrevFreq(static_cast<int>(startFreq))
{
  LoadTableBank(tablePaths);
  SetWavetable(0);
}

void WavetableOscillator<double>::LoadTableBank(std::vector<std::string> tablePaths)
{
  WavetableOscillator<double>::mWavetables.clear();
  for (auto path : tablePaths)
    WavetableOscillator<double>::mWavetables.push_back(path);
}

void WavetableOscillator<double>::SetWavetable(int tableIdx)
{
  assert(WavetableOscillator<double>::mWavetables.size() > tableIdx);
  SetWavetable((WavetableOscillator<double>::mWavetables.at(tableIdx)));
}

void WavetableOscillator<double>::SetWavetable(std::string path)
{
  //WavFile wav(path);
  //SetWavetable(wav);
  WtFile wt(path);
  SetWavetable(wt);
}

void WavetableOscillator<double>::SetWavetable(WavFile& wav)
{
  if (wav.Success())
  {
    std::unique_lock<std::mutex> lock(mWtMutex);
    mWtReady = false;
    if (mWT != nullptr)
      delete mWT;
    mWT = new Wavetable<double>(wav, 2);
    mWtReady = true;
  }
}