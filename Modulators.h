#pragma once

#include "VectorFunctions.h"

#include "IPlugParameter.h"
#include "LFO.h"
#include "ADSREnvelope.h"

#include <cstdarg>

template<int DynMods = 6, int StatMods = 3>
class ParameterModulator
{
public:
  ParameterModulator(double min, double max, const char* name = "", bool exponential = false) : mMin(min), mMax(max), mName(name), mIsExponential(exponential)
  {
    if (mIsExponential)
    {
      mMin = std::max(mMin, 1e-6);
      mRange = std::log(mMax / mMin);
    }
    else
      mRange = mMax - mMin;

    mMinV = Vec4d(mMin);
    mMaxV = Vec4d(mMax);
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

    mMinV = Vec4d(mMin);
    mMaxV = Vec4d(mMax);
  }

  virtual void SetValue(int idx, double value)
  {
    mModDepths[idx] = value * mRange;
  }

  static inline void SetModValues(double* modPtr)
  {
    memcpy(ParameterModulator::mModValues, modPtr, DynMods * sizeof(double));
  }

  virtual void SetInitialValue(const double value)
  {
    mInitial = value;
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

  inline double AddModulation()
  {
    return AddModulation(mInitial);
  }

  /* TODO: This may be better implemented with virtual functions (see ParameterModulationExp implementation below) */
  inline double AddModulation(double initVal)
  {
    double modulation{ 0. };
    for (auto i{ 0 }; i < DynMods; ++i)
      modulation += mModDepths[i] * ParameterModulator::mModValues[i];
    modulation += mStaticModulation;
    if (!mIsExponential)
      return std::max(std::min(initVal + modulation, mMax), mMin);
    else
      return std::max(std::min(initVal * std::exp(modulation), mMax), mMin);
  }

  inline void SetStaticModulation(double* staticMods)
  {
    mStaticModulation = 0.;
    for (auto i{ 0 }; i < StatMods; ++i)
      mStaticModulation += mModDepths[i + DynMods] * staticMods[i];
  }

  inline double ClipToRange(double value)
  {
    return std::max(std::min(value, mMax), mMin);
  }

  double operator[](int idx)
  {
    return mModDepths[idx];
  }

  const double* Depths()
  {
    return mModDepths;
  }


protected:
  static inline double mModValues[DynMods + StatMods]{ 0. };
  double mStaticModulation{};

  double mInitial{ 0. };
  double mModDepths[DynMods + StatMods]{ 0. };

  double mMin{ 0. };
  double mMax{ 1. };
  double mRange{ 1. };
  std::string mName;
  const bool mIsExponential{ false };

  Vec4d mMinV = Vec4d(0.);
  Vec4d mMaxV = Vec4d(0.);
};

template<typename T>
class GenericModulator
{
public:
  GenericModulator() {}

  virtual inline T Process() = 0;

  virtual inline void SetParams(T* params) = 0;
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

template<typename T>
class Envelope : public GenericModulator<T>, public ADSREnvelope<T>
{
public:
  Envelope(const char* name = "", std::function<void()> resetFunc = nullptr, bool sustainEnabled = true) :
    ADSREnvelope<T>(name, resetFunc, sustainEnabled), GenericModulator<T>()
  {}

  inline void SetParams(T* params) override
  {
    SetSustain(params[0]);
  }

  inline void SetSustain(T susLvl)
  {
    mSustainLevel = susLvl;
  }

  inline T Process() override
  {
    return ADSREnvelope<T>::Process(mSustainLevel);
  }

protected:
  T mSustainLevel{ 0.5 };
};

/*
A fast LFO using lookup tables
*/
template<typename T = double>
class FastLFO : public LFO<T>, public GenericModulator<T>
{
  union tabfudge
  {
    double d;
    int i[2];
  } ALIGNED(8);

public:
  FastLFO(T initialValue=0.) : LFO<T>(), GenericModulator<T>(), mLastOutput(initialValue)
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

  static inline void SetTempoAndBeat(double qnPos = 0., bool transportIsRunning = false, double tempo = 120.)
  {
    mQNPos = qnPos;
    mTransportIsRunning = transportIsRunning;
    mTempo = tempo;
  }

  const double GetPhase()
  {
    return IOscillator<T>::mPhase;
  }

  void SetPhase(double newPhase)
  {
    IOscillator<T>::mPhase = newPhase;
  }


  virtual inline T Process()
  {
    return ProcessSynced(mQNPos, mTransportIsRunning, mTempo);
  }

  inline T Process(double freqHz) override
  {
    IOscillator<T>::SetFreqCPS(freqHz);

    return ProcessSynced(mQNPos, mTransportIsRunning, mTempo);
  }

  inline T ProcessSynced(double qnPos = 0., bool transportIsRunning = false, double tempo = 120.)
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
    mLastOutput = s_out;
    return s_out * mLevelScalar;
  }

  T GetLastOuput() const
  {
    return mLastOutput;
  }

  static inline T WrapPhase(T x, T lo = 0., T hi = 1.)
  {
    while (x >= hi)
      x -= hi;
    while (x < lo)
      x += hi - lo;
    return x;
  };

  inline void SetParams(T* params) override
  {
    SetScalar(params[1]);
  }

protected:
  static inline constexpr int mTableSize{ 1024 };
  static inline constexpr int mTableSizeM1{ 1023 };
  T mLastOutput;

  static inline T mLUT[LFO<T>::EShape::kNumShapes][mTableSize];

  static inline double mStaticPhase{ 0 };
  static inline double mTempo{ 120. };
  static inline bool mTransportIsRunning{ false };
  static inline double mQNPos{ 0. };
};

template <typename T, int NSteps=16>
class Sequencer : public FastLFO<T>
{
  enum ERateMode
  {
    kBPM,
    kHz
  };

public:
  Sequencer(T* stepValues) : FastLFO<T>()
  {
    mStepValues = stepValues;
  }

  void SetLength(const int numSteps)
  {
    mLength = numSteps;
  }

  void SetQNScalarFromDivision(int division) 
  {
    mQNScalar = GetQNScalar(static_cast<ETempoDivision>(Clip(division, 0, (int)kNumDivisions))) / NSteps; // Tempo-synced rate indicates the length of one step
  }

  inline void SetFreqCPS(double freqHz)
  {
    FastLFO<T>::SetFreqCPS(freqHz);
    CalculateGlideSamples();
  }

  /* Set the glide time, as a fraction of the time from one step to the next. (i.e. the fraction of the step-time required to reach the next step's value)
  */
  void SetGlide(double glideNorm)
  {
    mGlide = glideNorm;
    CalculateGlideSamples();
  }

  inline void CalculateGlideSamples()
  {
    mSamplesPerStep = 1. / mPhaseIncr / static_cast<double>(NSteps); // (cycles/sample) ^ -1 / (steps/cycle) = samples step^-1
    double samplesToTarget = mGlide * mSamplesPerStep;
    mGlidePerSample = std::min(1. / samplesToTarget, 1.);
  }

  inline T Process() override 
  {
    return ProcessSynced(FastLFO<T>::mQNPos, FastLFO<T>::mTransportIsRunning, FastLFO<T>::mTempo);
  }

  inline T GetStepValueWithGlide()
  {
    T prevValue = mPrevValue;
    T targValue = GetStepValue(mStepPos);
    T stepPhase = mPhase * NSteps - std::floor(mPhase * NSteps);
    return prevValue + std::min(mGlidePerSample * stepPhase * mSamplesPerStep, 1.) * (targValue - prevValue);
  }

  inline T ProcessSynced(double qnPos = 0., bool transportIsRunning = false, double tempo = 120.)
  {
    T oneOverQNScalar = 1. / LFO<T>::mQNScalar;
    T phase = IOscillator<T>::mPhase;
    bool tempoSync = FastLFO<T>::mRateMode == FastLFO<T>::ERateMode::kBPM;

    if (tempoSync && !transportIsRunning)
      IOscillator<T>::SetFreqCPS(tempo / 60.);

    double samplesPerBeat = IOscillator<T>::mSampleRate * (60.0 / (tempo == 0.0 ? 1.0 : tempo)); // samples per beat

    T phaseIncr = IOscillator<T>::mPhaseIncr;

    if (tempoSync)
    {
      double sampleAccurateQnPos = qnPos + (phase / samplesPerBeat);
      phaseIncr *= LFO<T>::mQNScalar;
      if (transportIsRunning)
        phase = std::fmod(sampleAccurateQnPos, oneOverQNScalar) * LFO<T>::mQNScalar;
    }

    T tempValue = GetStepValue(mStepPos);
    mStepPos = static_cast<int>(phase * NSteps) % mLength;
    T s_out = GetStepValueWithGlide();

    // When the step changes, save the previous step value (for glide calculations)
    if (mStepPos != mPrevStepPos)
    {
      mPrevValue = tempValue;
      mPrevStepPos = mStepPos;
    }

    // Increment phase
    IOscillator<T>::mPhase = FastLFO<T>::WrapPhase(IOscillator<T>::mPhase + phaseIncr, 0., 1.);
    return s_out * mLevelScalar;
  }

  inline T GetStepValue(int pos)
  {
    return mStepValues[pos];
  }

  int GetCurrentStep()
  {
    return mStepPos;
  }

protected:
  T* mStepValues;
  T mPrevValue{ 0. };
  int mStepPos{ 0 };
  int mPrevStepPos{ NSteps };
  int mPrevStep{ -1 };
  int mLength{ NSteps };
  T mGlide{ 0. };
  T mGlidePerSample{ 0. };
  T mSamplesPerStep;
};

template<typename T, class M>
class GlobalModulator : public Sequencer<T>
{
public:
  GlobalModulator() : GlobalModulator(nullptr) {}
  GlobalModulator(T* stepValues) : Sequencer<T>(stepValues) {}
  GlobalModulator(Sequencer<T>& sequencer) : Sequencer<T>(sequencer.mStepValues) {}

  inline T Process() override
  {
    return mBuffer.Get()[(mBlockPos++) % mBlockLength];
  }

  void Resize(int blockSize)
  {
    mBuffer.Resize(blockSize);
  }

  void FillBuffer(int nFrames)
  {
    mBlockPos = 0;
    mBlockLength = nFrames;
    for (int i{ 0 }; i < nFrames; ++i)
    {
      mBuffer.Get()[i] = M::Process();
    }
  }

  T GetLastOutput() const
  {
    return M::mLastOutput;
  }

private:
  WDL_TypedBuf<T> mBuffer;
  int mBlockPos{ 0 };
  int mBlockLength{ 0 };
};

END_IPLUG_NAMESPACE