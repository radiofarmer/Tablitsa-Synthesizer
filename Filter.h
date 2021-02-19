#pragma once

#include "SignalProcessing.h"

#include "vectorclass.h"

#include <vector>
#include <cmath>
#include <functional>


#define FILTER_TYPE_LIST "None", "VSF", "Moog Ladder", "Comb"
#define FILTER_MODE_LIST_VSF  "Lowpass", "Highpass", "Bandpass", "Allpass"
#define FILTER_MODE_LIST_MOOG  "Lowpass I", "Highpass I", "Bandpass I", "Lowpass II", "Highpass II", "Bandpass II"
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

/* General function for processing a biquad filter in direct form II*/
template<typename T>
inline T ProcessBiquadII(T s, T* a, T* b, DelayLine& z)
{
  T sum = s - a[0] * z[0] - a[1] * z[1];
  T out = sum * b[0] + z[0] * b[1] + z[1] * b[2];
  z.push(sum);
  return out;
}

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
    T overdrive = SoftClip(s, 1. + mDrive * 2.);
    s += mDrive * (overdrive - s);
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

    inline T ProcessDFII(T s)
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
      mMatrixDelay[1] = mMatrixDelay[0];
      mMatrixDelay[0] = sum;
      return out;
    }

    /* Process four samples simultaneously */
    inline T ProcessPrecomp(T* s)
    {
      Vec4d x = Vec4d().load(s);
    }

    /* Process one sample and return a vector the sample and the delay line */
    inline Vec4d __vectorcall ProcessPrecomp_Vector(T s)
    {
      T out = ProcessPrecomp(s);
      mDelayVector = blend4<0, 1, 4, V_DC>(Vec4d(out, sum), mDelayVector);
      return mDelayVector
    }

    void SetAddendPtr(double* ptr)
    {
      mAddends = ptr;
    }

    void SetDelayPtr(double* ptr)
    {
      mMatrixDelay = ptr;
    }

  protected:
    std::vector<double> mCoefs;
    DelayLine mZ{ 2 };
    Vec4d mB;
    Vec4d mA;
    Vec4d mDelayVector = Vec4d(0.);

    double* mAddends = nullptr;
    double* mMatrixDelay = nullptr;

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
      mSOS[i].SetDelayPtr(mDelayPtrs[i]);
      // Store the coefficients in this parent filter object
      z1[2 * i] = mSOS[i].mCoefs[4];
      z1[2 * i + 1] = mSOS[i].mCoefs[1];
      z2[2 * i] = mSOS[i].mCoefs[5];
      z2[2 * i + 1] = mSOS[i].mCoefs[2];
    }

    mZ1_Coefs1.load(z1);
    mZ1_Coefs2.load(z1 + 4);
    mZ1_Coefs3.load(z1 + 8);

    mZ2_Coefs1.load(z2);
    mZ2_Coefs2.load(z2 + 4);
    mZ2_Coefs3.load(z2 + 8);
  }

  inline void Precompute()
  {
    // Get values from the 6x2 delay matrix
    Vec4d flt01_taps = Vec4d().load(mSOS[0].mMatrixDelay);
    Vec4d flt23_taps = Vec4d().load(mSOS[2].mMatrixDelay);
    Vec4d flt45_taps = Vec4d().load(mSOS[4].mMatrixDelay);
    // Alternatively:
//  Vec4d z1_0 = blend<0, 0, 4, 4>(mSOS[0].mDelayVector, mSOS[1].mDelayVector);
    // etc...

    // Tap 1
    Vec4d z1_0 = permute4<0, 0, 2, 2>(flt01_taps);
    Vec4d z1_1 = permute4<0, 0, 2, 2>(flt23_taps);
    Vec4d z1_2 = permute4<0, 0, 2, 2>(flt45_taps);
    // Tap 2
    Vec4d z2_0 = permute4<1, 1, 3, 3>(flt01_taps);
    Vec4d z2_1 = permute4<1, 1, 3, 3>(flt23_taps);
    Vec4d z2_2 = permute4<1, 1, 3, 3>(flt45_taps);

    // First 4 terms
    Vec4d addends1 = (z1_0 * mZ1_Coefs1) + (z2_0 * mZ2_Coefs1);
    addends1.store(mAddendPtrs);
    // Next 4 terms
    Vec4d addends2 = (z1_1 * mZ1_Coefs2) + (z2_1 * mZ2_Coefs2);
    addends2.store(mAddendPtrs + 4);
    // Last 4 terms
    Vec4d addends3 = (z1_2 * mZ1_Coefs3) + (z2_2 * mZ2_Coefs3);
    addends3.store(mAddendPtrs + 8);
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

  /*__vectorcall inline Vec4d ProcessAndDownsample_Vector(Vec4d s)
  {

  }*/

private:
  Biquad<T> mSOS[6]{ { 1.0194978892532574e-05, 2.0389957785065147e-05, 1.0194978892532574e-05, 1.0, -1.6759049412919762, 0.7386853939104704 },
    { 1.0, 2.0, 1.0, 1.0, -1.4134265459338102, 0.7727803769380925 },
    { 1.0, 2.0, 1.0, 1.0, -1.0347899844575996, 0.8248787327348859 },
    { 1.0, 2.0, 1.0, 1.0, -0.6862027551047557, 0.8794409054482396 },
    { 1.0, 2.0, 1.0, 1.0, -0.4401516515832734, 0.9299883719435977 },
    { 1.0, 2.0, 1.0, 1.0, -0.317836352796945, 0.9768802444563486 } };
  double mAddendPtrs[12]{ 0. };
  double mDelayPtrs[6][2]{ 0. };

  Vec4d mZ1_Coefs1, mZ1_Coefs2, mZ1_Coefs3;
  Vec4d mZ2_Coefs1, mZ2_Coefs2, mZ2_Coefs3;
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
  DelayLine mDelayIn{ COMB_MAX_DELAY + 1 };
  DelayLine mDelayOut{ COMB_MAX_DELAY + 1 };
};

template<typename T, bool LowShelf=true>
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