#pragma once

#include "SignalProcessing.h"

#include "vectorclass.h"

#include <vector>
#include <cmath>
#include <functional>


#define FILTER_TYPE_LIST "None", "VSF", "Moog Ladder", "Comb"
#define FILTER_MODE_LIST_VSF  "Lowpass", "Highpass", "Bandpass", "Allpass"
#define FILTER_MODE_LIST_MOOG  "Lowpass", "Highpass", "Bandpass"
#define FILTER_MODE_LIST_COMB  "N/A"

#define COMB_MAX_DELAY 512


enum EFilters
{
  kNoFilter=0,
  kVSF,
  kMoog,
  kComb,
  kNumFilters,
};


template<typename T>
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

  virtual void SetCutoff(double cutoffNorm)
  {
    mFc = cutoffNorm;
  }

  virtual void SetQ(double q)
  {
    mQ = q;
  }

  virtual void SetDrive(double drive) {}

  virtual void CalculateCoefficients() {}

  virtual void SetMode(int mode)
  {
    mMode = mode;
    Reset();
  }

  virtual T Process(double s) { return s;  }

protected:
  const double pi{ 3.14159265359 };
  const double twoPi{ 6.28318530718 };
  T mSampleRate;
  T mFc;
  T mQ;
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

  inline void SetCutoff(double cutoffNorm)
  {
    mFc = std::max(std::min(cutoffNorm, 0.99), 0.001);
  }

  inline void SetQ(double q)
  {
    mQ = 1. / std::max(1 - q, 0.001);
  }

  void Reset()
  {
    mZ[0] = 0;
    mZ[1] = 0;
  }

  inline void CalculateCoefficients()
  {
    if (mMode == kBandpass)
      mQAdj = mFc / (mQ * mBandwidth);
    else
      mQAdj = mQ;

    const double twoPiFc{ twoPi * mFc };
    const double minusPiFcOverQ{ -twoPiFc / (2. * mQAdj) };

    mA = 2. * std::cos(twoPiFc) * std::exp(minusPiFcOverQ);
    mB = std::exp(2 * minusPiFcOverQ);

  }
#ifdef VECTOR_FLT
  inline Vec4d ProcessLowPass(Vec4d s)
  {
    Vec4d x;
  }
#else
  inline T ProcessLowpass(double s)
  {
    T s_out = mA * mZ[0] - mB * mZ[1] + (1 - mA + mB) * s;
    mZ[1] = mZ[0];
    mZ[0] = s_out;
    return s_out;
  }
#endif

  inline T ProcessHighpass(T s)
  {
    T sum = s * (1 - mB + mA) + mA * mZ[0] - mB * mZ[1];
    T s_out = sum - 2 * mZ[0] + mZ[1];
    mZ[1] = mZ[0];
    mZ[0] = sum;
    return s_out;
  }

  inline T ProcessBandpass(T s)
  {
    T sum = s * std::cos(pi * mFc) * std::sqrt(1 - mA + mB) + mA * mZ[0] - mB * mZ[1];
    T s_out = sum - 2 * mZ[0] + mZ[1];
    mZ[1] = mZ[0];
    mZ[0] = sum;
    return s_out;
  }

  inline T ProcessAllpass(T s)
  {
    T sum = mB * (mA * mZ[0] - mB * mZ[1] + s);
    T s_out = sum - mA * mZ[0] + mZ[1];
    mZ[1] = mZ[0];
    mZ[0] = sum;
    return s_out;
  }

  inline T Process(double s)
  {
    CalculateCoefficients();
    return (this->*mProcessFunctions[mMode])(s);
  }

private:
//  double mFc;
//  double mQ;
  T mA{ 0. };
  T mB{ 0. };
  T mBandwidth{ 0.25 };
  T mQAdj{ 0. };
  T mZ[2]{};
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

  inline void SetDrive(double drive) 
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
      mA = mD = mE = 0;
      mB = 2;
      mC = -2;
      break;
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

template <typename T>
class ChebyshevBL
{
public:
  template <typename T>
  class Biquad
  {
  public:
    Biquad(std::initializer_list<double> coefs) : mCoefs(coefs)
    {
      mB = Vec4d(mCoefs[0], mCoefs[1], mCoefs[2], 0.);
      mA = Vec4d(1., -1. * mCoefs[4], -1. * mCoefs[5], 0.);
    }

    inline T Process(T s)
    {
      // Direct Form 1
      T sum = mCoefs[0] * s + mCoefs[1] * mX[0] + mCoefs[2] * mX[1] - mCoefs[4] * mY[0] - mCoefs[5] * mY[1];
      mX.push(s);
      mY.push(sum);
      return sum;
    }

    inline T Process2(T s)
    {
      // Direct Form 2 - summed last two terms in each expression are the values stored in mAddends
      T sum = s - mCoefs[4] * mZ[0] - mCoefs[5] * mZ[1]; // SUM1
      T out = mCoefs[0] * sum + mCoefs[1] * mZ[0] + mCoefs[2] * mZ[1]; // SUM2
      mZ.push(sum);
      return out;
    }

    inline T ProcessPrecomp(T s)
    {
      // Direct Form 2
      T sum = s - mAddends[0]; // SUM1 
      T out = mCoefs[0] * sum + mAddends[1]; // SUM2
      mZ.push(sum);
      return out;
    }

    inline T Process2V(T s)
    {
      mDelay = blend4<4, 0, 1, 2>(mDelay, Vec4d(s, 0., 0., 0.)); // {s, z-1, z-2, z-3}
      Vec4d sum{ horizontal_add(mDelay * mA), 0., 0., 0. }; // {SUM1, 0., 0., 0.}
      mDelay = blend4<4, 1, 2, 3>(mDelay, sum); // {SUM1, z-1, z-2, z-3}
      return horizontal_add(mDelay * mB);
    }

    void SetAddendPtr(double* ptr)
    {
      mAddends = ptr;
    }

  protected:
    std::vector<double> mCoefs;
    DelayLine mZ{ 2 };
    Vec4d mB;
    Vec4d mA;
    Vec4d mDelay{ 0. };
    double* mAddends;

    friend class ChebyshevBL<T>;
  };

public:
  ChebyshevBL()
  {
    double z1[12];
    double z2[12];

    for (int i{0}; i < 6; ++i)
    {
      // Provide each second-order filter with a pointer to the addend terms
      mSOS[i].SetAddendPtr(&mAddendPtrs[2 * i]);
      // Store the coefficients in the parent filter object
      z1[2 * i] = mSOS[i].mCoefs[4];
      z1[2 * i + 1] = mSOS[i].mCoefs[1];
      z2[2 * i] = mSOS[i].mCoefs[5];
      z2[2 * i + 1] = mSOS[i].mCoefs[2];
    }

    mZ1_Coefs1.load(z1);
    mZ1_Coefs2.load(z1 + 8);
    mZ2_Coefs1.load(z2);
    mZ2_Coefs2.load(z2 + 8);
  }

  inline void Precompute()
  {
    // Tap 1
    Vec8d z1_1( mSOS[0].mZ[0], mSOS[0].mZ[0], mSOS[1].mZ[0], mSOS[1].mZ[0], mSOS[2].mZ[0], mSOS[2].mZ[0], mSOS[3].mZ[0], mSOS[3].mZ[0] );
    Vec4d z1_2( mSOS[4].mZ[0], mSOS[4].mZ[0], mSOS[5].mZ[0], mSOS[5].mZ[0] );
    // Tap 2
    Vec8d z2_1( mSOS[0].mZ[1], mSOS[0].mZ[1], mSOS[1].mZ[1], mSOS[1].mZ[1], mSOS[2].mZ[1], mSOS[2].mZ[1], mSOS[3].mZ[1], mSOS[3].mZ[1] );
    Vec4d z2_2( mSOS[4].mZ[1], mSOS[4].mZ[1], mSOS[5].mZ[1], mSOS[5].mZ[1] );

    // First 8 terms
    Vec8d addends1 = (z1_1 * mZ1_Coefs1) + (z2_1 * mZ2_Coefs1);
    addends1.store(mAddendPtrs);
    // Last 4 terms
    Vec4d addends2 = (z1_2 * mZ1_Coefs2) + (z2_2 * mZ2_Coefs2);
    addends2.store(mAddendPtrs + 8);
  }

  inline T Process(T s)
  {
    T sum = s;
    Precompute();
    for (int i{0}; i < 6; ++i)
      sum = mSOS[i].ProcessPrecomp(sum);
    return sum;
  }

  inline T ProcessAndDownsample(T* s)
  {
    Process(s[0]); // Process and throw away sample
    return Process(s[1]); // Process and return sample
  }

private:
  Biquad<T> mSOS[6]{ { 1.0194978892532574e-05, 2.0389957785065147e-05, 1.0194978892532574e-05, 1.0, -1.6759049412919762, 0.7386853939104704 },
    { 1.0, 2.0, 1.0, 1.0, -1.4134265459338102, 0.7727803769380925 },
    { 1.0, 2.0, 1.0, 1.0, -1.0347899844575996, 0.8248787327348859 },
    { 1.0, 2.0, 1.0, 1.0, -0.6862027551047557, 0.8794409054482396 },
    { 1.0, 2.0, 1.0, 1.0, -0.4401516515832734, 0.9299883719435977 },
    { 1.0, 2.0, 1.0, 1.0, -0.317836352796945, 0.9768802444563486 } };
  double mAddendPtrs[12]{};
  std::vector<DelayLine*> mDelayPtrs;
  Vec8d mZ1_Coefs1;
  Vec8d mZ2_Coefs1;
  Vec4d mZ1_Coefs2;
  Vec4d mZ2_Coefs2;
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

  inline virtual void SetCutoff(T cutoffNorm)
  {
    mFc = cutoffNorm * 2;
  }

  inline virtual void SetDrive(T delayLength)
  {
    SetDelay(std::min(static_cast<int>(delayLength), mMaxDelay));
  }

  inline void SetDelay(int samples)
  {
    mDelayLength = samples;
  }

  inline T Process(T s)
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
  DelayLine mDelayIn{ COMB_MAX_DELAY + 1 };
  DelayLine mDelayOut{ COMB_MAX_DELAY + 1 };
};

template<typename T, bool LowShelf=true>
class ShelvingFilter : public Filter<T>
{
public:
  ShelvingFilter(double sampleRate=41000., double cutoff=0.5, double resonance=0., double gain=0., bool cutoffIsNormalized = true) :
    Filter<T>(mSampleRate, cutoff, resonance, cutoffIsNormalized), mGain(gain)
  {
    CalculateCoefficients();

    // this should probably switch the coefficient calculation function rather than the process one
    if (LowShelf)
      mProcessFunction = [this](T s) { return ProcessLowShelf(s); };
    else
      mProcessFunction = [this](T s) { return ProcessHighShelf(s); };
  }

  inline void SetGain(T gain)
  {
    mGain = gain;
    CalculateCoefficients();
  }

  inline void CalculateCoefficients() override
  {
    double k = std::tan(Filter<T>::pi * mFc);
    double v0 = std::pow(10., mGain / 20.);
    double slope = 1. / (1. - std::min(0.99, mQ));

    double kSqr = k * k;
    double vkSqr = v0 * kSqr;
    double denom = 1 + slope * k + kSqr;

    // Biquad Coefficients (low shelf)
    mB[0] = 1. + (std::sqrt(v0) * slope * k + vkSqr) / denom;
    mB[1] = 2. * (vkSqr - 1.) / denom;
    mB[2] = (1. - std::sqrt(v0) * slope * k + vkSqr) / denom;
    mA[0] = 2. * (kSqr - 1.) / denom;
    mA[1] = (1 - slope * k * kSqr) / denom;
  }

  inline T Process(T s) override
  {
    return mProcessFunction(s);
  }

  inline T ProcessLowShelf(T s)
  {
    T sum = s - mA[0] * mZ[0] - mA[1] * mZ[1];
    T out = sum * mB[0] + mZ[0] * mB[1] + mZ[1] * mB[2];
    mZ.push(sum);
    return out;
  }

  inline T ProcessHighShelf(T s)
  {
    return s;
  }

private:
  double mGain;
  double mA[2]{ 0. };
  double mB[3]{ 0. };
  DelayLine mZ{ 2 };
  std::function<T(T)> mProcessFunction;
};