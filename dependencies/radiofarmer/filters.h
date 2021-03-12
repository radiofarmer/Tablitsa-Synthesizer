#pragma once

#include "radiofarmer_config.h"

#include <cmath>

BEGIN_DSP_NAMESPACE

enum EOnePoleModes
{
  kLowpass1P,
  kHighpass1P,
  kBandpass1P,
  kNotch1P,
  kAllpass1P
};

/* Base class for one-pole topology-preserving transform filters */
class OnePoleTPTFilter
{
protected:
  static constexpr sample_t pi = 3.141592653589793;
  static constexpr sample_t piOver2 = 1.57079632679;
public:
  OnePoleTPTFilter(sample_t sampleRate=DEFAULT_SRATE, sample_t fc=0.25) : mSampleRate(sampleRate), mFc(fc)
  {
    CalculateCoefficients();
  }

  void SetSampleRate(const sample_t sampleRate)
  {
    mSampleRate = sampleRate;
    CalculateCoefficients();
  }

  void SetCutoff(const sample_t fc)
  {
    mFc = fc;
    CalculateCoefficients();
  }

  /* Cutoff prewarping, to be specialized for different filter types */
  const sample_t CutoffAdj() const
  {
    return mFc * piOver2;
  }

  void CalculateCoefficients()
  {
    const sample_t g = CutoffAdj();
    mG = g / (g + 1.);
  }

  void CalculateMatrix();

  sample_v __vectorcall Process4Lowpass(sample_t* x);
  sample_v __vectorcall Process4Lowpass(sample_v x);

protected:
  sample_t mSampleRate;
  sample_t mFc;
  sample_t mG;
  sample_t mZ{ 0. };
  sample_t mCoefMat[5][4];
};

class LowpassOnePole : public OnePoleTPTFilter
{
public:
  LowpassOnePole(sample_t sampleRate = DEFAULT_SRATE, sample_t fc = 0.25) :
    OnePoleTPTFilter(sampleRate, fc) {}

  sample_t CutoffAdj();

  sample_t Process(sample_t x);

  void Process4(sample_t* x);
  sample_v __vectorcall Process4(sample_v x);
};

class HighpassOnePole : public OnePoleTPTFilter
{
public:
  HighpassOnePole(sample_t sampleRate = DEFAULT_SRATE, sample_t fc = 0.25, sample_t fc_max=1.) :
    OnePoleTPTFilter(sampleRate, fc), mFcMax(fc_max), mFcMaxMapped(std::tan(fc_max / 2.)) {}

  sample_t CutoffAdj();

  sample_t Process(sample_t x);
  void Process4(sample_t* x);

protected:
  sample_t mFcMax;
  const sample_t mFcMaxMapped;
};

class PeakOnePole : public OnePoleTPTFilter
{
public:
  PeakOnePole(sample_t sampleRate = DEFAULT_SRATE, sample_t fc = 0.25) :
    OnePoleTPTFilter(sampleRate, fc) {}

  void SetPeakGain(sample_t g)
  {
    mPeakGain = g;
  }
  sample_t CutoffAdj();

  sample_t Process(sample_t x);
  void Process4(sample_t* x);
  sample_v __vectorcall Process4(sample_v x_v);

private:
  sample_t mPeakGain{ 1. };
};

END_DSP_NAMESPACE