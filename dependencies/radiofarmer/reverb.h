#pragma once

#include "radiofarmer_config.h"
#include "allpass.h"

#include <assert.h>

BEGIN_DSP_NAMESPACE

template<int NM=3, int NS=2>
class ReverbCascade
{
public:
  ReverbCascade(sample_t sampleRate=DEFAULT_SRATE, sample_t maxDelayMS=50., sample_t minDelayMS = 10., sample_t maxFeedback=0.75, sample_t minFeedback=0.6, sample_t gain=0.5) :
    mSampleRate(sampleRate),
    mMaxDelay(maxDelayMS),
    mMinDelay(minDelayMS),
    mMaxFeedback(maxFeedback),
    mMinFeedback(minFeedback),
    mGain(gain)
  {
    SetSampleRate(mSampleRate);
    AdjustAllFilters();
  }

  void SetSampleRate(sample_t sampleRate)
  {
    mSampleRate = sampleRate;
    for (int i{ 0 }; i < NM; ++i)
    {
      mMonoFilters[i].SetSampleRate(mSampleRate);
    }
    for (int i{ 0 }; i < NS; ++i)
    {
      mStereoFilters[i][0].SetSampleRate(mSampleRate);
      mStereoFilters[i][1].SetSampleRate(mSampleRate);
    }
  }

  void SetDelay(sample_t maxDelay, sample_t minDelay = 10., const bool adjust=true)
  {
    assert(minDelay > 1. && maxDelay > minDelay && "Delay value out of range!");
    mMaxDelay = maxDelay;
    mMinDelay = minDelay;
    if (adjust)
      AdjustAllFilters();
  }

  void SetFeedback(sample_t maxFeedback, sample_t minFeedback = 0.6, const bool adjust=true)
  {
    assert(maxFeedback >= minFeedback && maxFeedback < 1. && minFeedback >= 0. && "Feedback value out of range!");
    mMaxFeedback = maxFeedback;
    mMinFeedback = minFeedback;
    if (adjust)
      AdjustAllFilters();
  }

  void SetGain(sample_t gain)
  {
    mGain = gain;
  }

  void AdjustAllFilters()
  {
    const int n = NM + NS;
    sample_t delays[n];
    sample_t fb[n];

    const sample_t delayRange = mMaxDelay - mMinDelay;
    const sample_t fbRange = mMaxFeedback - mMinFeedback;

    delays[0] = mMaxDelay;
    fb[0] = mMaxFeedback;

    for (int i{ 1 }; i < n; ++i)
    {
      delays[i] = delays[i-1] * 0.78;
      fb[i] = mMaxFeedback - static_cast<sample_t>(std::rand() % 100) / 1000.;
    }

    // Mono filters
    for (int i{ 0 }; i < NM; ++i)
    {
      mMonoFilters[i].SetDelayMS(delays[i]);
      mMonoFilters[i].SetFeedbackGain(fb[i]);
    }
    for (int i{ 0 }; i < NS; ++i)
    {
      mStereoFilters[i][0].SetDelayMS(delays[i + NM]);
      mStereoFilters[i][1].SetDelayMS(delays[i + NM]);

      mStereoFilters[i][0].SetFeedbackGain(fb[i + NM]);
      mStereoFilters[i][1].SetFeedbackGain(fb[i + NM]);
    }
  }

  inline sample_t Process(sample_t s)
  {
    StereoSample<> s_stereo{ s, s };
    ProcessStereo(s_stereo);
    return 0.5 * (s_stereo.l + s_stereo.r);
  }

  inline void ProcessStereo(StereoSample<>& s)
  {
    sample_t sum1 = mGain * (s.l + s.r);

    // Process mono filters
    sample_t monoLine{ sum1 };
    for (int i{ 0 }; i < NM; ++i)
    {
      monoLine = mMonoFilters[i].Process(monoLine);
    }
    // Process StereoFilters
    sample_t l_in, r_in;
    l_in = r_in = monoLine;
    for (int i{ 0 }; i < NS; ++i)
    {
      l_in = mStereoFilters[i][0].Process(l_in);
      r_in = mStereoFilters[i][1].Process(r_in);
    }

    s.l += l_in;
    s.r += r_in;
  }

protected:
  sample_t mSampleRate;
  sample_t mMaxDelay;
  sample_t mMinDelay;
  sample_t mMaxFeedback;
  sample_t mMinFeedback;
  sample_t mGain;

  Allpass1 mMonoFilters[NM];
  Allpass1 mStereoFilters[NS][2];
};

END_DSP_NAMESPACE