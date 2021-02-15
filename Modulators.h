#pragma once

#include "VectorFunctions.h"

#include "IPlugParameter.h"
#include "LFO.h"
#include "ADSREnvelope.h"


template<typename T>
class GenericModulator
{
public:
  GenericModulator() {}

  virtual inline T Process() = 0;

  virtual inline void SetParams(T* params) = 0;
};

struct ModMetronome
{
  ModMetronome() {}

  void Set(double qnPos, double tempo, bool transportIsRunning)
  {
    mQNPos = qnPos;
    mTempo = tempo;
    mTransportIsRunning = transportIsRunning;
  }

  double mQNPos{ 0. };
  double mTempo{ iplug::DEFAULT_TEMPO };
  bool mTransportIsRunning;
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
  static constexpr T DECAY_FLOOR{ 0.00001 };
  static constexpr T HALF_LIVES_TO_FLOOR{ 16.6096404744 }; // Half-lives elapsed until magnitude is less than the floor value
  static constexpr T HALF_LIFE_SCALAR{ 1. / HALF_LIVES_TO_FLOOR };

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

  void SetStageTime(int stage, T timeMS)
  {
    switch (stage)
    {
    case kAttack:
      mAttackIncr = ADSREnvelope<T>::CalcIncrFromTimeLinear(Clip(timeMS, MIN_ENV_TIME_MS, MAX_ENV_TIME_MS), mSampleRate);
      break;
    case kDecay:
      mDecayIncr = CalcIncrFromTimeExp(Clip(timeMS, MIN_ENV_TIME_MS, MAX_ENV_TIME_MS), mSampleRate);
      mDecayIncrLinear = CalcIncrFromTimeLinear(Clip(timeMS, MIN_ENV_TIME_MS, MAX_ENV_TIME_MS), mSampleRate);
      break;
    case kRelease:
      mReleaseIncr = CalcIncrFromTimeExp(Clip(timeMS, MIN_ENV_TIME_MS, MAX_ENV_TIME_MS), mSampleRate);
      mReleaseIncrLinear = CalcIncrFromTimeLinear(Clip(timeMS, MIN_ENV_TIME_MS, MAX_ENV_TIME_MS), mSampleRate);
      break;
    default:
      //error
      break;
    }
  }

  void SetStageCurve(int stage, T expAmt)
  {
    switch (stage)
    {
    case kDecay:
      mDecCurve = Clip(expAmt, 0., 1.);
      break;
    case kRelease:
      mRelCurve = Clip(expAmt, 0., 1.);
      break;
    default:
      //error
      break;
    }
  }

  /** Process the envelope, returning the value according to the current envelope stage
  * @param sustainLevel Since the sustain level could be changed during processing, it is supplied as an argument, so that it can be smoothed extenally if nessecary, to avoid discontinuities */
  inline T Process(T sustainLevel = 0.)
  {
    T result = 0.;

    switch (mStage)
    {
    case kIdle:
      result = mEnvValue;
      break;
    case kAttack:
      mEnvValue += (mAttackIncr * mScalar);
      if (mEnvValue > ENV_VALUE_HIGH || mAttackIncr == 0.)
      {
        mStage = kDecay;
        mEnvValue = 1.;
      }
      result = mEnvValue;
      break;
    case kDecay:
      mEnvValue -= ((mDecayIncrLinear + mDecCurve * (mDecayIncr * mEnvValue - mDecayIncrLinear)) * mScalar);
      result = (mEnvValue * (1. - sustainLevel)) + sustainLevel;
      if (mEnvValue < ENV_VALUE_LOW)
      {
        if (mSustainEnabled)
        {
          mStage = kSustain;
          mEnvValue = 1.;
          result = sustainLevel;
        }
        else
          Release();
      }
      break;
    case kSustain:
      result = sustainLevel;
      break;
    case kRelease:
      mEnvValue -= ((mReleaseIncrLinear + mRelCurve * (mReleaseIncr * mEnvValue - mReleaseIncrLinear)) * mScalar);
      if (mEnvValue < ENV_VALUE_LOW || mReleaseIncr == 0.)
      {
        mStage = kIdle;
        mEnvValue = 0.;

        if (mEndReleaseFunc)
          mEndReleaseFunc();
      }
      result = mEnvValue * mReleaseLevel;
      break;
    case kReleasedToRetrigger:
      mEnvValue -= mRetriggerReleaseIncr;
      if (mEnvValue < ENV_VALUE_LOW)
      {
        mStage = kAttack;
        mLevel = mNewStartLevel;
        mEnvValue = 0.;
        mPrevResult = 0.;
        mReleaseLevel = 0.;

        if (mResetFunc)
          mResetFunc();
      }
      result = mEnvValue * mReleaseLevel;
      break;
    case kReleasedToEndEarly:
      mEnvValue -= mEarlyReleaseIncr;
      if (mEnvValue < ENV_VALUE_LOW)
      {
        mStage = kIdle;
        mLevel = 0.;
        mEnvValue = 0.;
        mPrevResult = 0.;
        mReleaseLevel = 0.;
        if (mEndReleaseFunc)
          mEndReleaseFunc();
      }
      result = mEnvValue * mReleaseLevel;
      break;
    default:
      result = mEnvValue;
      break;
    }

    mPrevResult = result;
    mPrevOutput = (result * mLevel);
    return mPrevOutput;
  }

  inline T Process() override
  {
    return Process(mSustainLevel);
  }

protected:
  inline T CalcIncrFromTimeExp(T timeMS, T sr) const
  {
    T r;

    if (timeMS <= 0.0) return 0.;
    else
    {
      r = M_LN2 / (timeMS / 1000. * sr * HALF_LIFE_SCALAR);
      if (!(r < 1.0)) r = 1.0;

      return r;
    }
  }

protected:
  T mSustainLevel{ 0.5 };
  T mDecCurve{ 1. }; // Interpolation between linear (0.) and exponential (1.) curves
  T mRelCurve{ 1. };
  T mDecayIncrLinear;
  T mReleaseIncrLinear;
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

  FastLFO(T initialValue = 0.) : FastLFO<T>(nullptr, initialValue)
  {
  }

  FastLFO<T>(ModMetronome * metronome, T initialValue = 0.) :
    mMetronome(metronome), mLastOutput(initialValue),
    LFO<T>(), GenericModulator<T>()
  {
    WriteLUTs();

    // If a global metronome was not provided,  create a local metronome
    if (mMetronome == nullptr)
      mMetronome = new ModMetronome;
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
    return ProcessSynced(mMetronome->mQNPos, mMetronome->mTransportIsRunning, mMetronome->mTempo);
  }

  inline T Process(double freqHz) override
  {
    IOscillator<T>::SetFreqCPS(freqHz);

    return Process();
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
    SetFreqCPS(params[0]);
    SetScalar(params[1]);
  }

protected:

  static inline constexpr int mTableSize{ 1024 };
  static inline constexpr int mTableSizeM1{ 1023 };
  T mLastOutput;

  static inline T mLUT[LFO<T>::EShape::kNumShapes][mTableSize];

  ModMetronome* mMetronome{ nullptr };
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
  Sequencer(ModMetronome* metronome, T* stepValues) : FastLFO<T>(metronome), mStepValues(stepValues) {}

  Sequencer(T* stepValues) : FastLFO<T>(), mStepValues(stepValues)
  {
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
    return ProcessSynced(mMetronome->mQNPos, mMetronome->mTransportIsRunning, mMetronome->mTempo);
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
  GlobalModulator(ModMetronome* metronome) : Sequencer<T>(metronome, nullptr) {}
  GlobalModulator(ModMetronome* metronome, T* stepValues) : Sequencer<T>(metronome, stepValues) {}
  GlobalModulator(Sequencer<T>& sequencer) : Sequencer<T>(sequencer.mMetronome, sequencer.mStepValues) {}

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