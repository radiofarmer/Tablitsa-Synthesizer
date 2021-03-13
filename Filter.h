#pragma once

#include "radiofarmerDSP.h"

#include "VectorFunctions.h"

#include <vector>
#include <cmath>
#include <functional>

using namespace radiofarmer;

#define FILTER_TYPE_LIST "None", "VSF", "Moog Ladder", "Comb"
#define FILTER_MODE_LIST_VSF  "Lowpass", "Highpass", "Bandpass", "Allpass"
#define FILTER_MODE_LIST_MOOG  "Lowpass I", "Highpass I", "Bandpass I", "Lowpass II", "Highpass II", "Bandpass II"
#define FILTER_MODE_LIST_COMB  "N/A"

#define COMB_MAX_DELAY 1024


enum EFilters
{
  kNoFilter=0,
  kVSF,
  kMoog,
  kComb,
  kNumFilters,
};

/* General function for processing a biquad filter in direct form II*/
template<typename T>
inline T ProcessBiquadII(T s, T* a, T* b, DelayLine<2>& z)
{
  T sum = s - a[0] * z[0] - a[1] * z[1];
  T out = sum * b[0] + z[0] * b[1] + z[1] * b[2];
  z.push(sum);
  return out;
}

template<typename T, class V=void>
class Filter
{
public:
  Filter(double sampleRate=44100., double cutoff = 0.5, double q = 0., bool cutoffIsNormalized = true) : mSampleRate(sampleRate), mFc(cutoff), mQ(q)
  {
    if (!cutoffIsNormalized)
    {
      mFc /= mSampleRate;
    };
  }

  virtual void Reset() {}

  virtual void SetSampleRate(double sampleRate)
  {
    SetCutoff(mFc * mSampleRate / sampleRate);
    mSampleRate = sampleRate;
  }

  virtual inline void SetCutoff(double cutoffNorm)
  {
    mFc = cutoffNorm;
  }

  virtual inline void SetQ(double q)
  {
    mQ = q;
  }

  virtual inline void SetDrive(double drive) { mDrive = drive; }

  virtual inline void CalculateCoefficients() {}

  virtual void SetMode(int mode)
  {
    mMode = mode;
    Reset();
  }

  virtual T Process(double s) { return s;  }
  virtual Vec4d __vectorcall Process_Vector(Vec4d s) { return s; }

protected:
  const double pi{ 3.14159265359 };
  const double twoPi{ 6.28318530718 };
  T mSampleRate;
  T mFc;
  T mQ;
  T mDrive{ 0. };
  int mMode{ 0 };
};

template<typename T>
class NullFilter : public Filter<T>
{
public:
  NullFilter(double sampleRate = 44100.) : Filter<T>(sampleRate) {}

  void Reset() {}
  
  inline T Process(double s)
  {
    return s;
  }
};

template<typename T>
class SVF2 : public Filter<T>
{
  enum EModes
  {
    kLowpass = 0,
    kHighpass,
    kBandpass,
    kAllpass,
    kNumModes
  };

public:
  typedef T (SVF2::* SVF2ProcessFunc)(T s);

  SVF2(double sampleRate=44100., double cutoff = 0.5, double q = 0., bool cutoffIsNormalized = true) : Filter<T>(sampleRate, cutoff, q, cutoffIsNormalized)
  {
    SetCutoff(mFc);
    SetQ(mQ);
    SetMode(kLowpass);
  }

  void SetCutoff(double cutoffNorm) override
  {
    mFc = std::max(std::min(cutoffNorm, 0.99), 0.001);
  }

  void SetQ(double q) override
  {
#ifdef STATE_SPACE_FILTER
    mQ = q;
#else
    if (mMode == kBandpass)
      mQ = mFc / (mMaxBandwidth * (1. - q));
    else
      mQ = 1. / std::max(1. - q, 0.001);
#endif
  }

  void SetDrive(double drive) override
  {
    mDrive = drive * 2.;
  }

  void Reset()
  {
    mZ.reset();
  }

  inline void CalculateCoefficients()
  {
    const double twoPiFc{ twoPi * mFc };
    const double minusPiFcOverQ{ -twoPiFc / (2. * mQ) };

    mA = 2. * std::cos(twoPiFc) * std::exp(minusPiFcOverQ);
    mB = std::exp(2 * minusPiFcOverQ);

    if (mMode == kBandpass)
      mC = std::cos(pi * mFc) * std::sqrt(1 - mA + mB);
    else if (mMode == kAllpass)
      mC = 1.;
    else
      mC = 1 - mA + mB;
  }

  inline T ProcessLowpass(T s)
  {
    T s_out = mA * mZ[0] - mB * mZ[1] + mC * s;
    mZ.push(s_out);
    return s_out;
  }

  inline T ProcessHighpass(T s)
  {
    T sum = mA * mZ[0] - mB * mZ[1] + s; // Do not multiply s by C!
    T s_out = sum - 2. * mZ[0] + mZ[1];
    mZ.push(sum);
    return s_out;
  }

  inline T ProcessBandpass(T s)
  {
    T sum = s * mC + mA * mZ[0] - mB * mZ[1];
    T s_out = sum - mZ[0];
    mZ.push(sum);
    return s_out;
  }

  inline T ProcessAllpass(T s)
  {
    T sum = mB * (mA * mZ[0] - mB * mZ[1] + s);
    T s_out = sum - mA * mZ[0] + mZ[1];
    mZ.push(sum);
    return s_out;
  }

  inline T Process(double s)
  {
    CalculateCoefficients();
    T overdrive = SoftClip<T, 5>(s * (1. + mDrive));
    s += mDrive * (overdrive - s);
    return (this->*mProcessFunctions[mMode])(s);
  }

private:
  T mA{ 0. };
  T mB{ 0. };
  T mC{ 1. };

  const double mMaxBandwidth{ 0.05 };
  DelayLine<2> mZ;
  SVF2ProcessFunc mProcessFunctions[kNumFilters]{ &SVF2::ProcessLowpass, &SVF2::ProcessHighpass, &SVF2::ProcessBandpass, &SVF2::ProcessAllpass };
};

template<typename T>
class MoogLadder : public Filter<T>
{
  enum EModes
  {
    kLowpass12db=0,
    kHighpass12db,
    kBandpass12db,
    kLowpass24db,
    kHighpass24db,
    kBandpass24db
  };

#pragma mark Moog One-Pole -
public:
  template<typename T>
  class MoogOnePole : public Filter<T>
  {
  public:
    MoogOnePole(double sampleRate = 44100., double cutoff = 0.5, double q = 0., bool cutoffIsNormalized = true) : Filter<T>(sampleRate, cutoff, q, cutoffIsNormalized)
    {}

    void Reset()
    {
      mFF = 0.;
      mFB = 0.;
    }

    inline double Process(double s)
    {
      const double sum = s * mCoeff1 + mFF * mCoeff2 - mFB;
      const double out = sum * mG + mFB;
      mFB = out;
      mFF = s;
      return out;
    }

    double mG{};

  private:
    const double mCoeff1{ 1. / 1.3 };
    const double mCoeff2{ 0.3 / 1.3 };
    Vec4d mCoeffsV{ 1. / 1.3, 0.3 / 1.3, 0., 0. };
    double mFF{ 0 };
    double mFB{ 0 };
  };

public:
  MoogLadder(double sampleRate = 44100., double cutoff = 0.5, double q = 0., bool cutoffIsNormalized = true, double drive = 0.) :
    Filter<T>(sampleRate, cutoff, q, cutoffIsNormalized), mGComp(drive)
  {
    for (int i{ 0 }; i < 4; ++i)
    {
      mLadder[i] = new MoogOnePole<T>(sampleRate, cutoff, q, cutoffIsNormalized);
    }
  }

  inline void SetCutoff(double cutoffNorm) 
  {
    mFc = cutoffNorm * pi;
  }

  inline void SetDrive(double drive) override
  {
    mGComp = drive;
  }

  void SetMode(int mode)
  {
    switch (mode)
    {
    case kLowpass12db:
      mA = mB = mD = mE = 0.;
      mC = 1.;
      break;
    case kHighpass12db:
      mA = mC = 1.;
      mB = -2.;
      mD = mE = 0.;
      break;
    case kBandpass12db:
      mA = mD = mE = 0.;
      mB = 2.;
      mC = -2.;
      break;
    case kLowpass24db:
      mA = mB = mC = mD = 0.;
      mE = 1.;
    case kHighpass24db:
      mA = mE = 1.;
      mB = mD = -4.;
      mC = 6.;
    case kBandpass24db:
      mA = mB = 0.;
      mC = mE = 4.;
      mD = -8.;
    default:
      break;
    }
  }

  void Reset()
  {
    mZ = 0.;
  }

  inline void CalculateCoefficients()
  {
    const double fcSquared{ mFc * mFc };
    const double fcCubed{ fcSquared * mFc };
    const double g{ 0.9892 * mFc - 0.4342 * fcSquared + 0.1381 * fcCubed - 0.0202 * fcCubed * mFc };

    // Set g of each one-pole filter
    for (int i{ 0 }; i < 4; ++i)
    {
      mLadder[i]->mG = g;
    }
    mGRes = mQ * (1.0029 + 0.0526 * mFc - 0.0926 * fcSquared + 0.0218 * fcCubed );
  }

  T Process(double s)
  {
    CalculateCoefficients();
    const double sumFF{ mZ - s * mGComp };
    const double sumFB{ s - 4. * mGRes * sumFF };

    const double nonlin{ std::tanh(sumFB) };
    const double step1{ mLadder[0]->Process(nonlin) };
    const double step2{ mLadder[1]->Process(step1) };
    const double step3{ mLadder[2]->Process(step2) };
    mZ = mLadder[3]->Process(step3);

    return mA * nonlin + mB * step1 + mC * step2 + mD * step3 + mE * mZ;
  }

private:
  MoogOnePole<T>* mLadder[4];
  double mG{};
  double mGRes{};
  double mGComp;
  double mA{ 0 };
  double mB{ 0 };
  double mC{ 1 };
  double mD{ 0 };
  double mE{ 0 };
  double mZ{ 0 };
};

template<typename T>
class CombFilter : public Filter<T>
{
public:
  CombFilter(double sampleRate = 44100., double feedforward = 0.5, double feedback = 0., bool cutoffIsNormalized = true, double delayLength = 8.) :
    Filter<T>(sampleRate, feedforward, feedback, cutoffIsNormalized), mDelayLength(static_cast<int>(delayLength))
  {
    SetDelay(mDelayLength);
  }

  inline void SetCutoff(T cutoffNorm) override
  {
    mFc = cutoffNorm * 2;
  }

  inline void SetDrive(T delayLength) override
  {
    SetDelay(std::min(static_cast<int>(delayLength * mMaxDelay), mMaxDelay));
  }

  inline void SetDelay(int samples)
  {
    mDelayLength = samples;
  }

  inline T Process(T s) override
  {
    T out = s + mFF * mDelayIn[mDelayLength] - mFB * mDelayOut[mDelayLength / 2];
    mDelayIn.push(s);
    mDelayOut.push(out);
    return out;
  }

private:
  const int mMaxDelay{ COMB_MAX_DELAY };
  int mDelayLength;
  T& mFF{ mFc }; // Feedforward used as alias for cutoff
  T& mFB{ mQ }; // Feedback used as alias for resonance
  DelayLine<COMB_MAX_DELAY> mDelayIn;
  DelayLine<COMB_MAX_DELAY> mDelayOut;
};

template<typename T, class V=void, bool LowShelf=true>
class ShelvingFilter : public Filter<T>
{
public:
  ShelvingFilter(double sampleRate=41000., double cutoff=0.1, double gain=0., bool cutoffIsNormalized = true) :
    Filter<T>(mSampleRate, cutoff, 0., cutoffIsNormalized), mGain(gain)
  {
    if (LowShelf)
      mCoefficientsFunc = [this](double& k, double& v0) { return LowShelfCoefficients(k, v0); };
    else
      mCoefficientsFunc = [this](double& k, double& v0) { return HighShelfCoefficients(k, v0); };

    CalculateCoefficients();
  }

  inline void SetGain(T gain)
  {
    mGain = gain;
    CalculateCoefficients();
  }

  inline void SetCutoff(double cutoffNorm) override
  {
    mFc = std::min(cutoffNorm, 0.4);
  }

  inline void ShiftCutoff(double shiftNorm)
  {
    SetCutoff(std::max(0.001, mFc + shiftNorm));
  }

  inline void CalculateCoefficients() override
  {
    double k = std::tan(Filter<T>::pi * mFc);
    double v0 = std::pow(10., mGain / 20.);

    mCoefficientsFunc(k, v0);
  }

  inline void LowShelfCoefficients(double& k, double& v0)
  {
    // Boost
    mA = (k - 1.) / (k + 1.);
    mH0d2 = (v0 - 1.) / 2.;
    mSign = 1.;
  }

  inline void HighShelfCoefficients(double& k, double& v0)
  {
    mA = (k - 1.) / (k + 1.);
    mH0d2 = (v0 - 1.) / 2.;
    mSign = -1.;
  }

  inline T Process(T s) override
  {
    T sum = mA * s + mX_prev - mA * mY_prev;
    T out = mH0d2 * (s + mSign * sum) + s;
    mX_prev = s;
    mY_prev = sum;
    return out;
  }

private:
  double mGain; // Gain in decibels
  double mA{ 0. };
  double mH0d2{ -1. };
  double mSign{ 1. };
  T mX_prev{ 0. };
  T mY_prev{ 0. };
  std::function<void(double&, double&)> mCoefficientsFunc;
};