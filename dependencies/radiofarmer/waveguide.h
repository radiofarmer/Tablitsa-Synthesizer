#pragma once

#include "radiofarmer_config.h"
#include <functional>
#include <cmath>

BEGIN_DSP_NAMESPACE

static inline sample_t fast_tanh(sample_t x)
{
  return x / (1. + std::abs(2 * x));
}

/* Basic allpass filters */
struct NonlinearAllpass
{
  sample_t a;
  sample_t g;
  sample_t z{0.};

  NonlinearAllpass() : a(1.), g(1.) {}

  inline sample_t s_func(sample_t x)
  {
    return std::sin(3.14159 * a * fast_tanh(g * x));
  }

  inline sample_t c_func(sample_t x)
  {
    return std::cos(3.14159 * a * fast_tanh(g * x));
  }
};

template<int N>
class AllpassLadder
{
public:
  AllpassLadder(sample_t sampleRate=DEFAULT_SRATE) : mSampleRate(sampleRate)
  {}

  void SetSampleRate(sample_t sampleRate)
  {
    mSampleRate = sampleRate;
  }

  void SetCoefs(sample_t a, sample_t g)
  {
    for (int i{ 0 }; i < N; ++i)
    {
      mLadder[i].a = a;
      mLadder[i].g = g;
    }
  }

  sample_t Process(sample_t s_in)
  {
    NonlinearAllpass* stage = mLadder;
    sample_t output[N];
    for (int i{ 0 }; i < N - 1; ++i)
    {
      const sample_t s_adj = stage->s_func(s_in);
      const sample_t c_adj = stage->c_func(s_in);
      output[i] = stage->z * c_adj + s_in * s_adj;
      s_in = s_in * c_adj - stage->z * s_adj;
      stage++;
    }

    // Last stage
    const sample_t s_adj = stage->s_func(s_in);
    const sample_t c_adj = stage->c_func(s_in);
    output[N-1] = stage->z * c_adj + s_in * s_adj;
    stage->z = s_in * c_adj - stage->z * s_adj;

    for (int i{ N - 1 }; i > 0; --i)
    {
      (--stage)->z = output[i];
    }

    return output[0];
  }

protected:
  sample_t mSampleRate;
  NonlinearAllpass mLadder[N];
};

END_DSP_NAMESPACE