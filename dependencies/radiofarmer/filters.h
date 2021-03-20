#pragma once

#include "radiofarmer_config.h"

#include <cmath>

BEGIN_DSP_NAMESPACE

static constexpr sample_t pi = 3.1415926;

struct Integrator
{
  sample_t g;
  sample_t gp1_recip;
  sample_t z{ 0. };

  Integrator(sample_t cutoffNorm = 0.)
  {
    Recalculate(cutoffNorm);
  }

  inline void Recalculate(sample_t cutoffNorm)
  {
    g = cutoffNorm * pi * 0.5;
    Recalculate();
  }

  inline void Recalculate()
  {
    gp1_recip = 1. / (1 + g);
  }

  inline sample_t Process(const sample_t x)
  {
    const sample_t out = (g * x + z) * gp1_recip;
    z = out - x;
    return out;
  }
};

enum EOnePoleModes
{
  kLowpass1P,
  kHighpass1P,
  kBandpass1P,
  kNotch1P,
  kAllpass1P
};

/* Base class for one-pole topology-preserving transform filters */
class OnePoleTPTFilter
{
protected:
  static constexpr sample_t pi = 3.141592653589793;
  static constexpr sample_t piOver2 = 1.57079632679;
public:
  OnePoleTPTFilter(sample_t sampleRate=DEFAULT_SRATE, sample_t fc=0.25) : mSampleRate(sampleRate), mFc(fc)
  {
    CalculateCoefficients();
  }

  void SetSampleRate(const sample_t sampleRate)
  {
    mSampleRate = sampleRate;
    CalculateCoefficients();
  }

  void SetCutoff(const sample_t fc)
  {
    mFc = fc;
    CalculateCoefficients();
  }

  /* Cutoff prewarping, to be specialized for different filter types */
  const sample_t CutoffAdj() const
  {
    return mFc * piOver2;
  }

  void CalculateCoefficients()
  {
    const sample_t g = CutoffAdj();
    mG = g / (g + 1.);
  }


  inline void CalculateMatrix()
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
#pragma clang loop unroll(full)
    for (int ii{ 0 }; ii < 5; ++ii)
      mCoefMat[ii][0] = coefs[0][ii];
    // Calculate and fill rows 1-3
#pragma clang loop unroll(full)
    for (int i{ 1 }; i < 4; ++i)
    {
#pragma clang loop unroll(full)
      for (int j{ 0 }; j < 5; ++j)
      {
        coefs[i][j] += oneMinusG * coefs[i - 1][j];
      }

      // Store transposed row
#pragma clang loop unroll(full)
      for (int ii{ 0 }; ii < 5; ++ii)
        mCoefMat[ii][i] = coefs[i][ii];
    }
  }

  sample_v __vectorcall Process4Lowpass(sample_t* x);
  sample_v __vectorcall Process4Lowpass(sample_v x);

protected:
  sample_t mSampleRate;
  sample_t mFc;
  sample_t mG;
  sample_t mZ{ 0. };
  sample_t mCoefMat[5][4];
};

class LowpassOnePole : public OnePoleTPTFilter
{
public:
  LowpassOnePole(sample_t sampleRate = DEFAULT_SRATE, sample_t fc = 0.25) :
    OnePoleTPTFilter(sampleRate, fc) {}

  sample_t CutoffAdj();

  sample_t Process(sample_t x);

  void Process4(sample_t* x);
  sample_v __vectorcall Process4(sample_v x);
};

class HighpassOnePole : public OnePoleTPTFilter
{
public:
  HighpassOnePole(sample_t sampleRate = DEFAULT_SRATE, sample_t fc = 0.25, sample_t fc_max=1.) :
    OnePoleTPTFilter(sampleRate, fc), mFcMax(fc_max), mFcMaxMapped(std::tan(fc_max / 2.)) {}

  sample_t CutoffAdj();

  sample_t Process(sample_t x);
  void Process4(sample_t* x);
  sample_v __vectorcall Process4(sample_v x);

protected:
  sample_t mFcMax;
  const sample_t mFcMaxMapped;
};

class PeakOnePole : public OnePoleTPTFilter
{
public:
  PeakOnePole(sample_t sampleRate = DEFAULT_SRATE, sample_t fc = 0.25) :
    OnePoleTPTFilter(sampleRate, fc) {}

  void SetPeakGain(sample_t g)
  {
    mPeakGain = g;
  }
  sample_t CutoffAdj();

  sample_t Process(sample_t x);
  void Process4(sample_t* x);
  sample_v __vectorcall Process4(sample_v x_v);

private:
  sample_t mPeakGain{ 1. };
};

class TwoPoleTPTFilter
{
private:
  static constexpr sample_t piOver2{ 1.57079632679 };

  inline void Process(const sample_t x)
  {
    m_hp = m_denom * (x - m_g1 * m_z1 - m_z2);
    sample_t u1 = m_g * m_hp;

    m_bp = u1 + m_z1;
    sample_t u2 = m_g * m_bp;

    m_lp = u2 + m_z2;

    // Set delays
    m_z1 = m_bp + u1;
    m_z2 = m_lp + u2;
  }

public:
  TwoPoleTPTFilter(sample_t sampleRate = DEFAULT_SRATE, sample_t cutoffNorm = 0., sample_t res = 1.) :
    mSampleRate(sampleRate), mFc(cutoffNorm), mQ(res)
  {}

  inline void SetSampleRate(const sample_t sampleRate)
  {
    mFc *= mSampleRate / sampleRate;
    mSampleRate = sampleRate;
    CalculateCoefficients();
  }

  inline void SetCutoff(const sample_t cutoffNorm)
  {
    mFc = cutoffNorm;
    CalculateCoefficients();
  }

  inline void SetResonance(const sample_t q)
  {
    mQ = 1. - q;
    CalculateCoefficients();
  }

  inline void SetCutoffAndResonance(const sample_t cutoffNorm, const sample_t q)
  {
    mFc = piOver2 * cutoffNorm;
    mQ = 1. - q;
    CalculateCoefficients();
  }

  inline void CalculateCoefficients()
  {
    m_g = mFc * 0.5 * pi;
    m_g1 = mQ * 2. + m_g;
    m_denom = 1. / (1. + 2 * mQ * m_g + m_g * m_g);
  }

  inline sample_t ProcessHP(const sample_t x)
  {
    Process(x);
    return m_hp;
  }
  inline sample_t ProcessBP(const sample_t x)
  {
    Process(x);
    return m_bp * mQ * 2.;
  }
  inline sample_t ProcessLP(const sample_t x)
  {
    Process(x);
    return m_lp;
  }
  inline sample_t ProcessAP(const sample_t x)
  {
    return x - 2. * ProcessBP(x);
  }

protected:
  sample_t mSampleRate;
  sample_t mFc;
  sample_t mQ;

  sample_t m_g{ 0. };
  sample_t m_g1{ 0. };
  sample_t m_denom{ 1. };

  sample_t m_hp{ 0. };
  sample_t m_bp{ 0. };
  sample_t m_lp{ 0. };

  sample_t m_z1{ 0. };
  sample_t m_z2{ 0. };
};

END_DSP_NAMESPACE