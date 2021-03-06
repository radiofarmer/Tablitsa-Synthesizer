#include "allpass.h"

#include <algorithm>

using namespace radiofarmer;

Allpass1::Allpass1(sample_t sampleRate, sample_t delayMS, sample_t fb) : mSampleRate(sampleRate), mDelayMS(delayMS), mFB(fb)
{
  SetSampleRate(sampleRate);
}

void Allpass1::SetSampleRate(sample_t sampleRate)
{
  mSampleRate = sampleRate;
  SetDelayMS(mDelayMS);
}

void Allpass1::SetDelayMS(sample_t delayMS)
{
  mDelayMS = delayMS;
  mDelay = static_cast<int>(mDelayMS * mSampleRate / 1000.);
}

void Allpass1::SetFeedbackGain(sample_t fb)
{
  mFB = std::clamp(fb, 0., 1.);
}

sample_t Allpass1::Process(const sample_t s)
{
  sample_t sum = s + mZ[mDelay] * mFB;
  sample_t out = mZ[mDelay] - mFB * sum;
  mZ.push(sum);
  return out;
}