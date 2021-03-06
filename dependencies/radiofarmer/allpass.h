#pragma once

#include "radiofarmer_config.h"
#include "delayline.h"

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

END_DSP_NAMESPACE