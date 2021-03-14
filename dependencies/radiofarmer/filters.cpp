#include "filters.h"

using namespace radiofarmer;

sample_v __vectorcall OnePoleTPTFilter::Process4Lowpass(sample_t* x)
{
  CalculateMatrix();
  const sample_v z1{ mZ };

  const sample_v x0{ x[0] };
  const sample_v x1{ x[1] };
  const sample_v x2{ x[2] };
  const sample_v x3{ x[3] };

  sample_v y{ 0. };
  y += sample_v().load(mCoefMat[0]) * x3;
  y += sample_v().load(mCoefMat[1]) * x2;
  y += sample_v().load(mCoefMat[2]) * x1;
  y += sample_v().load(mCoefMat[3]) * x0;
  y += sample_v().load(mCoefMat[4]) * z1;

  mZ = y[3];

  return y;
}


sample_v __vectorcall OnePoleTPTFilter::Process4Lowpass(sample_v x)
{
  CalculateMatrix();
  const sample_v z1{ mZ };

  const sample_v x0{ x[0] };
  const sample_v x1{ x[1] };
  const sample_v x2{ x[2] };
  const sample_v x3{ x[3] };

  sample_v y{ 0. };
  y += sample_v().load(mCoefMat[0]) * x3;
  y += sample_v().load(mCoefMat[1]) * x2;
  y += sample_v().load(mCoefMat[2]) * x1;
  y += sample_v().load(mCoefMat[3]) * x0;
  y += sample_v().load(mCoefMat[4]) * z1;

  mZ = y[3];

  return y;
}

/* Lowpass */

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

void LowpassOnePole::Process4(sample_t* x)
{
  sample_v y{ OnePoleTPTFilter::Process4Lowpass(x) };
  y.store(x);
}

sample_v __vectorcall LowpassOnePole::Process4(sample_v x)
{
  return OnePoleTPTFilter::Process4Lowpass(x);
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
  return x + mPeakGain *  2. * (x - lp);
}

void PeakOnePole::Process4(sample_t* x)
{
  sample_v y_v{ OnePoleTPTFilter::Process4Lowpass(x) };
  const sample_v x_v{ sample_v().load(x) };
  y_v = x_v + mPeakGain * 2. * (x_v - y_v);
  y_v.store(x);
}

sample_v __vectorcall PeakOnePole::Process4(sample_v x_v)
{
  sample_v y_v{ OnePoleTPTFilter::Process4Lowpass(x_v) };
  y_v = x_v + mPeakGain * 2. * (x_v - y_v);
  return y_v;
}