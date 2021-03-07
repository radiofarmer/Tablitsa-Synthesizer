#pragma once

#include "radiofarmer_config.h"
#include "delayline.h"
#include "biquad.h"

BEGIN_DSP_NAMESPACE

/* First-Order allpass filter */
class Allpass1
{
public:
  Allpass1(sample_t sampleRate = DEFAULT_SRATE, sample_t delayMS=10., sample_t fb=0.5);

  void SetSampleRate(sample_t sampleRate);
  void SetDelayMS(sample_t delayMS);
  void SetFeedbackGain(sample_t fb);

  sample_t Process(const sample_t s);

protected:
  sample_t mSampleRate;
  sample_t mDelayMS;
  sample_t mFB;
  int mDelay;

  DelayLine<16384> mZ;
};

class AbsorbantAllpass : public Allpass1
{
public:
  AbsorbantAllpass(sample_t sampleRate = DEFAULT_SRATE, sample_t lpFreq=0.02, sample_t lpGain=0., sample_t delayMS = 10., sample_t fb = 0.5);

  void SetLPGain(sample_t gain);

  void SetLPF(sample_t dbGain, sample_t freq, sample_t bandwidth)
  {
    mLpf.CalculateCoefficients(dbGain, freq, bandwidth);
  }

  void SetLPF()
  {
    SetLPF(mLpGain, mLpFreq, mBandwidth);
  }

  sample_t Process(const sample_t s);

private:
  sample_t mLpFreq;
  sample_t mLpGain;
  sample_t mBandwidth{ 2. };
  Biquad mLpf;
};

END_DSP_NAMESPACE