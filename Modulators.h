#pragma once

#include "IPlugParameter.h"
#include "LFO.h"

#include <vector>

enum EModulators
{
  kEnv1=0,
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

BEGIN_IPLUG_NAMESPACE
/*
Shaping object for modulation depth sliders that provides finer control for parameters
with large ranges, such as filter cutoff frequencies.
*/
struct ModShapePowCurve : public IParam::ShapePowCurve
{
  ModShapePowCurve(double shape) : IParam::ShapePowCurve(shape) {}

  Shape* Clone() const override { return new ModShapePowCurve(*this); }

  double NormalizedToValue(double value, const IParam& param) const
  {
    double sign = value >= 0.5 ? 1. : -1.;
    value = std::abs(0.5 - value) * 2.;
    return std::pow(value, mShape) * sign;
  }

  double ValueToNormalized(double value, const IParam& param) const
  {
    double sign = value >= 0. ? 1. : -1.;
    value = std::abs(value);
    return 0.5 + std::pow(value, 1.0 / mShape) / 2. * sign;
  }

};

END_IPLUG_NAMESPACE

class ParameterModulator
{
public:
  ParameterModulator(double min, double max, bool exponential=false) : mMin(min), mMax(max), mIsExponential(exponential)
  {
    if (mIsExponential)
    {
      mMin = std::max(mMin, 1e-6);
      mRange = std::log(mMax / mMin);
    }
    else
      mRange = mMax - mMin;
  }

  void SetMinMax(double min, double max)
  {
    mMin = min;
    mMax = max;
    if (mIsExponential)
    {
      mRange = std::log(mMax / std::max(mMin, 1.e-6));
    }
    else
      mRange = mMax - mMin;
  }

  virtual void SetValue(int idx, double value)
  {
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


  inline double GetValue(double env1 = 0., double env2 = 0., double ampEnv = 0., double lfo1 = 0., double lfo2 = 0., double sequencer = 0.)
  {
    return std::max(std::min(mInitialValue +
      env1 * mEnv1Depth +
      env2 * mEnv2Depth +
      ampEnv * mAmpEnvDepth +
      lfo1 * mLFO1Depth +
      lfo2 * mLFO2Depth, mMax), mMin);
  }

  /*
  Write a buffer of modulation buffers. Can be called directly by a ParameterModulator, or indirectly by a ParameterModulatorList.
  @param inputs The initival values (parameter values) for the parameter
  @param outputs A buffer (e.g. WDL_TypedBuf) to write the values to
  @param nFrames The length of the block
  */
  inline void ProcessBlock(double* inputs, double* outputs, int nFrames)
  {
    for (auto i{ 0 }; i < nFrames; ++i)
    {
      outputs[i] = AddModulation(inputs[i]);
    }
  }

  /* TODO: This may be better implemented with virtual functions (see ParameterModulationExp implementation below) */
  inline double AddModulation(double initVal)
  {
    if (mIsExponential)
      return std::max(std::min(initVal * std::exp(mModValues[kEnv1] * mEnv1Depth + mModValues[kEnv2] * mEnv2Depth + mModValues[kAmpEnv] * mAmpEnvDepth +
        mModValues[kLFO1] * mLFO1Depth + mModValues[kLFO2] * mLFO2Depth + mModValues[kSequencer] * mSequencerDepth), mMax), mMin);
    else
      return std::max(std::min(initVal + mModValues[kEnv1] * mEnv1Depth + mModValues[kEnv2] * mEnv2Depth + mModValues[kAmpEnv] * mAmpEnvDepth +
        mModValues[kLFO1] * mLFO1Depth + mModValues[kLFO2] * mLFO2Depth + mModValues[kSequencer] * mSequencerDepth, mMax), mMin);
  }

  inline double AddModulationExp(double initVal)
  {
    return std::max(std::min(initVal * std::exp(mModValues[kEnv1] * mEnv1Depth + mModValues[kEnv2] * mEnv2Depth + mModValues[kAmpEnv] * mAmpEnvDepth +
      mModValues[kLFO1] * mLFO1Depth + mModValues[kLFO2] * mLFO2Depth + mModValues[kSequencer] * mSequencerDepth), mMax), mMin);
  }

  double operator[](int idx)
  {
    switch (idx) {
    case kEnv1:
      return mEnv1Depth;
    case kEnv2:
      return mEnv2Depth;
    case kAmpEnv:
      return mAmpEnvDepth;
    case kLFO1:
      return mLFO1Depth;
    case kLFO2:
      return mLFO2Depth;
    case kSequencer:
      return mSequencerDepth;
    default:
      break;
    }
  }

  static inline void SetModValues(double* modPtr)
  {
    memcpy(ParameterModulator::mModValues + 1, modPtr, kNumMods * 8);
  }

protected:
  static inline double mModValues[EModulators::kNumMods];

  double mInitialValue{ 0. };
  double mEnv1Depth{ 0. };
  double mEnv2Depth{ 0. };
  double mAmpEnvDepth{ 0. };
  double mLFO1Depth{ 0. };
  double mLFO2Depth{ 0. };
  double mSequencerDepth{ 0. };
  double mVelocityDepth{ 0. };
  double mKeytrackDepth{ 0. };
  double mRandomDepth{ 0. };

  double mMin{ 0. };
  double mMax{ 1. };
  double mRange{ 1. };
  const bool mIsExponential{ false };
};

class ParameterModulatorExp : public ParameterModulator
{
public:
  ParameterModulatorExp(double min, double max) : ParameterModulator(std::min(min, 1e-6), max)
  {
    mRange = std::log(mMax / mMin);
  }

  inline double AddModulation(double initVal)
  {
    return std::max(std::min(initVal * std::exp(mModValues[kEnv1] * mEnv1Depth + mModValues[kEnv2] * mEnv2Depth + mModValues[kAmpEnv] * mAmpEnvDepth +
      mModValues[kLFO1] * mLFO1Depth + mModValues[kLFO2] * mLFO2Depth + mModValues[kSequencer] * mSequencerDepth), mMax), mMin);
  }
};

template<typename T, int NParams=1>
class ParameterModulatorList
{
public:
  ParameterModulatorList() {}

  /** Write a buffer for each of several ParameterModulator objects.
  @param inputs Pointers to the buffers of initial (unmodulated) parameter values, e.g. from the inputs buffer in the block-processing function that calls this one.s
  @param outputs Pointers to the output buffers, of dimensions (NParams, nFrames)
  @params nFrames Length of the buffer, in samples
  */
  void ProcessBlock(T** inputs, T** ouputs, int nFrames)
  {
    for (auto i{ 0 }; i < NParams; ++i)
    {
      mModulations[i]->ProcessBlock(inputs[i], outputs[i], nFrames);
    }
  }

private:
  ParameterModulator* mModulations[NParams];
} WDL_FIXALIGN;

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