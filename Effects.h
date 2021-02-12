#pragma once

#include "SignalProcessing.h"
#include "Filter.h"



template <typename T>
class SaturationEQ
{
public:
  SaturationEQ(double sampleRate=41000.) : mLowShelf(sampleRate, 0.05, 0., 1.5)
  {
  }

  void SetSampleRate(double sampleRate)
  {
    mLowShelf.SetSampleRate(sampleRate);
  }

  inline void SetLevel(T lvl)
  {
    mGain = 1. + lvl;
    mLowShelf.SetGain(lvl * 10.);
  }

  inline T Process(T s)
  {
    return std::tanh(mLowShelf.Process(s) * mGain);
  }

private:
  T mGain;
  ShelvingFilter<T, true> mLowShelf;
};
