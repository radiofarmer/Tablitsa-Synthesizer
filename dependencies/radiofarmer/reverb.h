#pragma once

#include "radiofarmer_config.h"
#include "allpass.h"

#include <assert.h>
#include <cmath>

BEGIN_DSP_NAMESPACE

template<int NM=3, int NS=2>
class CascadeReverb
{
public:
  CascadeReverb(sample_t sampleRate=DEFAULT_SRATE, sample_t maxDelayMS=50., sample_t minDelayMS = 10., sample_t maxFeedback=0.75, sample_t minFeedback=0.6, sample_t gain=0.5) :
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

class UDFNReverb
{
  static constexpr int MaxDelaySamples{ 16384 };
  static constexpr int NEarlyDelays{ 5 };
  static constexpr int NLateDelays{ 1 };
  static constexpr int NLateAPFilters{ 5 };
  static constexpr int LateDelayPos{ 4 }; // The late delay comes between the fourth and fifth AP filters

  static constexpr sample_t DefaultEarlyDelayTimes[6]{
    23., 27., 29., 33., 35., 37.
  };
  static constexpr sample_t DefaultLateDelayTimes[7]{
    51., 55., 59., 61., 63., 67., 71.
  };

  static constexpr sample_t Primes[]{ 1., 2., 3., 5., 7., 11., 13., 17. };
  static constexpr int NPrimes{ sizeof(Primes) / 8 };

  static constexpr sample_t MinEarlyDelay{ 15. };
  static constexpr sample_t MaxEarlyDelay{ 100. };
  static constexpr sample_t MinLateDelay{ 45. };
  static constexpr sample_t MaxLateDelay{ 6000. };

  static inline sample_t PseudoPrimeScalar(const sample_t scale, const int idx)
  {
    return Primes[idx] + scale * (Primes[idx + 1] - Primes[idx]);
  }

public:
  UDFNReverb(sample_t sampleRate, sample_t lpFreq = 0.2);

  // Processing Parameters
  void SetSampleRate(sample_t sampleRate);
  void SetDelays();
  void SetAAPFilters();
  void ProcessStereo(StereoSample<sample_t>& s);

  // Reverb Quality Parameters
  void SetDiffusion(const sample_t diff);
  void SetDecayTime(const sample_t tNorm);
  void SetDamping(const sample_t damp);
  void SetHFDecay(const sample_t freqNorm);
  void SetColor(const sample_t freqNorm);
  void SetEarlyReflectionsLevel(const sample_t erLevel) { mERMix = erLevel; }
  void SetMixLevel(const sample_t mix) { mMix = mix; }

  void CalcNorm();

protected:
  sample_t mSampleRate;

  // Early reflections delay lines
  Biquad mEarlyLPF[2];
  DelayLine<MaxDelaySamples> mEarlyDelays[2][5];
  AbsorbantAllpass mEarlyAP[2];

  // Late reverb delay lines and filters (one per channel)
  AbsorbantAllpass mLateAP[2][5];
  DelayLine<MaxDelaySamples> mLateDelays[2];
  Biquad mLateLPF[2];

  // Coefficients
  sample_t mLPCutoff{ 0.1 };
  sample_t mG[2]{ 0.5, 0.5 };
  sample_t m_t[NLateAPFilters + 1]{ 1. };
  sample_t m_a{ 0.9 };
  sample_t m_g{ 0.6 };
  sample_t mNorm{ m_g * m_g + (1 - m_g * m_g) * ((m_a * m_a) / (1 - m_a * m_a * m_g * m_g)) };
  sample_t mERMix = 1.;
  sample_t mMix = 1.;

  sample_t mEarlyDelayTimes[6]{
    23., 27., 29., 33., 35., 37.
  };
  sample_t mLateDelayTimes[7]{
    51., 55., 59., 61., 63., 67., 71.
  };
};

END_DSP_NAMESPACE