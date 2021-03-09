#include "filters.h"

using namespace radiofarmer;

sample_t LowpassOnePole::CutoffAdj()
{
  const sample_t k = std::tan(mFc * piOver2);
  return k;
}

sample_t LowpassOnePole::Process(sample_t x)
{
  const sample_t sum = mG * (x - mZ);
  const sample_t out = mZ + sum;
  mZ = out + sum;
  return out;
}

/* Highpass */

sample_t HighpassOnePole::CutoffAdj()
{
  const sample_t k1 = std::tan(mFc * piOver2);
  const sample_t k2 = mFcMaxMapped + mFcMaxMapped * (k1 * k1 + 1.);
  return std::min(k1, k2);
}

sample_t HighpassOnePole::Process(sample_t x)
{
  const sample_t sum = mG * (x - mZ);
  const sample_t out = mZ + sum;
  mZ = out + sum;
  return x - out;
}

/* Peak */

sample_t PeakOnePole::CutoffAdj()
{
  const sample_t k = std::tan(mFc * piOver2);
  return k;
}

sample_t PeakOnePole::Process(sample_t x)
{
  const sample_t sum = mG * (x - mZ);
  const sample_t lp = mZ + sum;
  mZ = lp + sum;
  const sample_t hp = x - lp;
  return x + mPeakGain *  2. * (x - lp);
}