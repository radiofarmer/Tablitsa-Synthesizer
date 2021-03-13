#include "filters.h"

using namespace radiofarmer;

void OnePoleTPTFilter::CalculateMatrix()
{
  CalculateCoefficients();
  // y[n] = y[n-1] + g * (x[n] - y[n-1])
  // y[n] = g * x[n] + (1 - g) * y[n-1]
  /*
  *         x[n+3]  x[n+2]  x[n+1]  x[n]  y[n+2]  y[n+1]  y[n]  y[n-1]
  y[n]    | 0       0       0       g     0       0       0     1-g
  y[n+1]  | 0       0       g       0     0       0       1-g   0
  y[n+2]  | 0       g       0       0     0       1-g     0     0
  y[n+3]  | g       0       0       0     1-g     0       0     0
  */
  const sample_t oneMinusG = 1. - mG;
  sample_t coefs[4][5]{
    {0.,    0.,   0.,   mG,   oneMinusG},
    {0.,    0.,   mG,   0.,   0.},
    {0.,    mG,   0.,   0.,   0.},
    {mG,    0.,   0.,   0.,   0.}
  };
  // Store first row (transposed)
  for (int ii{ 0 }; ii < 5; ++ii)
    mCoefMat[ii][0] = coefs[0][ii];
  // Calculate and fill rows 1-3
  for (int i{ 1 }; i < 4; ++i)
  {
    for (int j{ 0 }; j < 5; ++j)
    {
      coefs[i][j] += oneMinusG * coefs[i - 1][j];
    }

    // Store transposed row
    for (int ii{ 0 }; ii < 5; ++ii)
      mCoefMat[ii][i] = coefs[i][ii];
  }
}

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