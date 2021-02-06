#pragma once

#define FILTER_TYPE_LIST "VSF 24dB/Oct", "Moog Ladder", "Comb"

template<typename T>
class Filter
{
  Filter(double sampleRate = 44100, double cutoff=0.5, double q=0., bool cutoffIsNormalized=true) : mSampleRate(sampleRate), mFc(cutoff), mQ(q)
  {
    if (!cutoffIsNormalized)
    {
      mFc /= mSampleRate;
    }
  }

  void UpdateSampleRate(double sampleRate)
  {
    mSampleRate = sampleRate;
  }

  virtual void SetCutoff(double cutoffNorm) = 0;

  virtual void SetQ(double q) = 0;

  virtual T Evaluate(T s) = 0;

protected:
  double mFc;
  double mQ;
  double mSampleRate;
};

template<typename T>
class VSF24dB : Filter
{
  VSF() : Filter()
  {}

  T Evaluate(T s) override
  {
    return s;
  }
};