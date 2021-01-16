#pragma once

#include "IPlugParameter.h"
#include "LFO.h"

#include <vector>

enum EModulators
{
  kInitial=0,
  kEnv1,
  kEnv2,
  kAmpEnv,
  kLFO1,
  kLFO2,
  kSequencer,
  kVelocity,
  kKeytrack,
  kTriggerRandom,
  kNumMods
};

/* THIS IS A BODGE. Make a proper map function class/struct at some point. */
#define MOD_MAP_LINEAR [](double value){ return value; }

struct MapFunction
{
  MapFunction() {}
  virtual MapFunction* Copy() = 0;
  virtual double Map(double valueNorm) = 0;
};

struct MapFunctionLinear : public MapFunction
{
  MapFunctionLinear() {}
  virtual MapFunction* Copy() { return new MapFunctionLinear(*this); }
  virtual double Map(double valueNorm) { return valueNorm; };
};

struct MapFunctionPower : public MapFunction
{
  MapFunctionPower(double exponent=3.) : mExp(exponent) {}
  virtual MapFunction* Copy() { return new MapFunctionPower(*this); }
  virtual double Map(double valueNorm) { return std::pow(valueNorm, mExp); }
  double mExp;
};

class ParameterModulator
{
public:
  ParameterModulator() {}

  ParameterModulator(double min, double max, MapFunction& map=MapFunctionLinear()) : mMin(min), mMax(max)
  {
    mRange = mMax - mMin;
    mMapFunction = std::unique_ptr<MapFunction>(map.Copy());
  }

  void SetMinMax(double min, double max)
  {
    mMin = min;
    mMax = max;
    mRange = mMax - mMin;
  }

  void SetValue(int idx, double value)
  {
    if (idx == kInitial)
      mInitialValue = value;
    else
    {
      // Scale value according to mapping function
      //...
      switch (idx) {
      case kEnv1:
        mEnv1Depth = value * mRange;
        break;
      case kEnv2:
        mEnv2Depth = value * mRange;
        break;
      case kAmpEnv:
        mAmpEnvDepth = value * mRange;
        break;
      case kLFO1:
        mLFO1Depth = value * mRange;
        break;
      case kLFO2:
        mLFO2Depth = value * mRange;
        break;
      default:
        break;
      }
    }
  }

  inline double GetValue(double env1 = 0., double env2 = 0., double ampEnv = 0., double lfo1 = 0., double lfo2 = 0., double sequencer = 0.)
  {
    return std::max(std::min(mInitialValue +
      env1 * mEnv1Depth +
      env2 * mEnv2Depth +
      ampEnv * mAmpEnvDepth +
      lfo1 * mLFO1Depth +
      lfo2 * mLFO2Depth, mMax), mMin);
  }

  inline double AddModulation(double initVal)
  {
    return std::max(std::min(initVal + mModValues[kEnv1] * mEnv1Depth + mModValues[kEnv2] * mEnv2Depth + mModValues[kAmpEnv] * mAmpEnvDepth +
      mModValues[kLFO1] * mLFO1Depth + mModValues[kLFO2] * mLFO2Depth + mModValues[kSequencer] * mSequencerDepth, mMax), mMin);
  }

  inline double AddModulation()
  {
    return std::max(std::min(mInitialValue + mModValues[kEnv1] * mEnv1Depth + mModValues[kEnv2] * mEnv2Depth + mModValues[kAmpEnv] * mAmpEnvDepth +
      mModValues[kLFO1] * mLFO1Depth + mModValues[kLFO2] * mLFO2Depth, mMax), mMin);
  }

  double operator[](int idx)
  {
    switch (idx) {
    case kInitial:
      return mInitialValue;
    case kEnv1:
      return mEnv1Depth;
    default:
      break;
    }
  }

  static inline void SetModValues(double* modPtr)
  {
    memcpy(ParameterModulator::mModValues + 1, modPtr, kNumMods * 8);
  }

private:
  static inline double mModValues[EModulators::kNumMods];

  double mInitialValue{ 0 };
  double mEnv1Depth{ 0 };
  double mEnv2Depth{ 0 };
  double mAmpEnvDepth{ 0 };
  double mLFO1Depth{ 0 };
  double mLFO2Depth{ 0 };
  double mSequencerDepth{ 0 };

  double mMin{ 0. };
  double mMax{ 1. };
  double mRange{ 1. };
  std::unique_ptr<MapFunction> mMapFunction;
};

BEGIN_IPLUG_NAMESPACE

/*
A fast LFO using lookup tables
*/
template<typename T = double>
class FastLFO : public LFO<T>
{
  union tabfudge
  {
    double d;
    int i[2];
  } ALIGNED(8);

public:
  FastLFO(T initialValue=0.) : LFO<T>(), mLastOutput(initialValue)
  {
    WriteLUTs();
  }

  static void WriteLUTs()
  {
    // Populate lookup tables
    for (int i{ 0 }; i < mTableSize; ++i)
    {
      double x = static_cast<double>(i) / static_cast<double>(mTableSize);
      FastLFO<T>::mLUT[LFO<T>::EShape::kTriangle][i] = (2. * (1. - std::abs((WrapPhase(x + 0.25) * 2.) - 1.))) - 1.;
      FastLFO<T>::mLUT[LFO<T>::EShape::kSquare][i] = std::copysign(1., x - 0.5);;
      FastLFO<T>::mLUT[LFO<T>::EShape::kRampUp][i] = (x * 2.) - 1.;
      FastLFO<T>::mLUT[LFO<T>::EShape::kRampDown][i] = ((1. - x) * 2.) - 1.;
      FastLFO<T>::mLUT[LFO<T>::EShape::kSine][i] = std::sin(x * 6.283185307179586);
    }
  }

  inline T Lookup(double phaseNorm)
  {
    phaseNorm += (double)UNITBIT32;

    union tabfudge tf;
    tf.d = UNITBIT32;
    const int normhipart = tf.i[HIOFFSET];

    tf.d = phaseNorm;
    const T* addr = FastLFO<T>::mLUT[mShape] + (tf.i[HIOFFSET] & mTableSizeM1);
    tf.i[HIOFFSET] = normhipart;
    const double frac = tf.d - UNITBIT32;
    const T f1 = addr[0];
    const T f2 = addr[1];
    return f1 + frac * (f2 - f1);
  }

  inline T Process(double freqHz) override
  {
    IOscillator<T>::SetFreqCPS(freqHz);

    return Process();
  }

  inline T Process(double qnPos = 0., bool transportIsRunning = false, double tempo = 120.)
  {
    T oneOverQNScalar = 1. / LFO<T>::mQNScalar;
    T phase = IOscillator<T>::mPhase;

    if (mRateMode == ERateMode::kBPM && !transportIsRunning)
      IOscillator<T>::SetFreqCPS(tempo / 60.);

    double samplesPerBeat = IOscillator<T>::mSampleRate * (60.0 / (tempo == 0.0 ? 1.0 : tempo)); // samples per beat

    T phaseIncr = IOscillator<T>::mPhaseIncr;

    double sampleAccurateQnPos = qnPos + (phase / samplesPerBeat);

    if (mRateMode == ERateMode::kBPM)
    {
      phaseIncr *= LFO<T>::mQNScalar;
      if (transportIsRunning)
        phase = std::fmod(sampleAccurateQnPos, oneOverQNScalar) * LFO<T>::mQNScalar;
    }

    T s_out = Lookup(phase);
    IOscillator<T>::mPhase = WrapPhase(IOscillator<T>::mPhase + phaseIncr * mTableSize, 0., static_cast<double>(mTableSize));
    return s_out;
  }

private:
  static inline T WrapPhase(T x, T lo = 0., T hi = 1.)
  {
    while (x >= hi)
      x -= hi;
    while (x < lo)
      x += hi - lo;
    return x;
  };

private:
  static inline constexpr int mTableSize{ 1024 };
  static inline constexpr int mTableSizeM1{ 1023 };
  T mLastOutput;
  static inline T mLUT[LFO<T>::EShape::kNumShapes][mTableSize];
//  LFO<T>::EShape mShape{ LFO<T>::EShape::kTriangle };
};

END_IPLUG_NAMESPACE

template <typename T>
class Sequencer
{
public:
  Sequencer(const int nSteps) : mNumSteps(nSteps)
  {
    mValues = new T[mNumSteps]{ 0. };
  }



private:
  int mNumSteps;
  T* mValues;
};