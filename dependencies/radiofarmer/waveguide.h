#pragma once

#include "radiofarmer_config.h"
#include "filters.h"

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
    mLP.SetSampleRate(sampleRate);
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
    sample_t ipt_lp = mLP.Process(s_in);
    NonlinearAllpass* stage = mLadder;
    sample_t output[N];

    for (int i{ 0 }; i < N - 1; ++i)
    {
      const sample_t s_adj = stage->s_func(ipt_lp);
      const sample_t c_adj = stage->c_func(ipt_lp);
      output[i] = stage->z * c_adj + s_in * s_adj;
      s_in = s_in * c_adj - stage->z * s_adj;
      stage++;
    }

    // Last stage
    const sample_t s_adj = stage->s_func(ipt_lp);
    const sample_t c_adj = stage->c_func(ipt_lp);
    output[N-1] = stage->z * c_adj + s_in * s_adj;
    stage->z = s_in * c_adj - stage->z * s_adj;

    for (int i{ N - 1 }; i > 0; --i)
    {
      (--stage)->z = output[i];
    }

    return output[0];
  }

  sample_v Process4(sample_v x)
  {
    sample_v ipt_lp = mLP.Process4(x);
    sample_t outputs[4];

    for (int s{ 0 }; s < 4; ++s)
    {
      NonlinearAllpass* stage = mLadder;
      sample_t stage_outputs[N];
      sample_t s_in = x[s];

      for (int i{ 0 }; i < N - 1; ++i)
      {
        const sample_t s_adj = stage->s_func(ipt_lp[s]);
        const sample_t c_adj = stage->c_func(ipt_lp[s]);
        stage_outputs[i] = stage->z * c_adj + s_in * s_adj;
        s_in = s_in * c_adj - stage->z * s_adj;
        stage++;
      }

      // Last stage
      const sample_t s_adj = stage->s_func(ipt_lp[s]);
      const sample_t c_adj = stage->c_func(ipt_lp[s]);
      stage_outputs[N - 1] = stage->z * c_adj + s_in * s_adj;
      stage->z = s_in * c_adj - stage->z * s_adj;

      for (int i{ N - 1 }; i > 0; --i)
      {
        (--stage)->z = stage_outputs[i];
      }

      outputs[s] = stage_outputs[0];
    }

    return sample_v().load(outputs);
  }

protected:
  sample_t mSampleRate;
  NonlinearAllpass mLadder[N];
  LowpassOnePole mLP{ DEFAULT_SRATE, 0.1 };
};

END_DSP_NAMESPACE