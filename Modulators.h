#pragma once

#include "IPlugParameter.h"
#include "LFO.h"
#include "VectorFunctions.h"

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

template<int DynMods=6, int StatMods=3>
class ParameterModulator
{
public:
  ParameterModulator(double min, double max, const char* name="", bool exponential=false) : mMin(min), mMax(max), mName(name), mIsExponential(exponential)
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

  static inline void SetModValues(double* modPtr)
  {
    memcpy(ParameterModulator::mModValues, modPtr, DynMods * sizeof(double));
  }

protected:
  static inline double mModValues[EModulators::kNumMods]{ 0. };
  double mStaticModulation{};

  double mModDepths[DynMods + StatMods]{ 0. };

  double mMin{ 0. };
  double mMax{ 1. };
  double mRange{ 1. };
  WDL_String mName;
  const bool mIsExponential{ false };

  Vec4d mMinV = Vec4d(0.);
  Vec4d mMaxV = Vec4d(0.);
};

template<typename T, int NParams = 1, int DynamicMods=6, int StaticMods=3, int VectorSize=4>
class ModulatedParameterList
{
public:
  ModulatedParameterList(std::initializer_list<ParameterModulator<DynamicMods, StaticMods>*> params) : mParams(params), mParamBlocks(NParams / 4)
  {}

  inline void SetStaticModulation(T* staticMods)
  {
    for (auto& p : mParams)
      p->SetStaticModulation(staticMods);
  }

  /** Write a buffer for each of several ParameterModulator objects.
  @param inputs_params Pointers to the buffers of initial (unmodulated) parameter values, e.g. from the inputs buffer in the block-processing function that calls this one. (NParams, nFrames)
  @param inputs_mods Pointers to the modulation values for the block (kNumMods, nFrames)
  @param outputs Pointers to the output buffers, of dimensions (NParams, nFrames)
  @params nFrames Length of the buffer, in samples
  */
  void ProcessBlock(T** param_inputs, T** mod_inputs, T** outputs, int nFrames, const int offset=0)
  {
    for (auto i{ 0 }; i < nFrames; ++i)
    {
      T modDepths[NParams]{ 0. };
      for (auto m{ 0 }; m < DynamicMods; ++m)
        modDepths[m] = mod_inputs[m][i];
      ParameterModulator<DynamicMods, StaticMods>::SetModValues(modDepths);
      for (auto p{ offset }; p < NParams; ++p)
      {
        outputs[p][i] = mParams[p]->AddModulation(param_inputs[p][i]);
      }
    }
  }

  /* Vector processing of each block of parameter values.

  TODO: Pad the modulation matrix with extra rows of zeros so that its rows are always divisible by the vector length. (See ModulatorList)
  TODO: Support for exponentially-scaled moduation
  */
  void ProcessBlockVec4d(T** param_inputs, T** mod_inputs, T** outputs, int nFrames)
  {
    for (auto i{ 0 }; i < nFrames; i+=4)
    {
      for (auto p{ 0 }; p < NParams; ++p)
      {
        /* For calculating per-sample modulation for single samples in parallel
        Vec4d modDepths;
        modDepths.load(mParams[p]->Depths()); // Modulation depths for the current parameters (constant for entire sample block)
        Vec4d modVals(mod_inputs[0][i], mod_inputs[1][i], mod_inputs[2][i], mod_inputs[3][i]); // Values of the modulators themselves for the current samples
        outputs[p][i] = mParams[p]->ClipToRange(horizontal_add(modDepths * modVals) + param_inputs[p][i]);*/

        // For calculating per-sample modulation for multiple samples in parallel
        Vec4d modSum(0.);
        Vec4d initVals;
        initVals.load(param_inputs[p]);
        for (auto m{ 0 }; m < 4; ++m)
        {
          Vec4d modSamples;
          modSamples.load(mod_inputs[m][i]);
          modSum += modSamples * mParams[p][m];
        }
        Vec4d outputV = Params[p]->ClipToRange(initVals + modSum);
        outputV.store(&outputs[p][i]);

      }
    }
    // Process the last few parameters, which don't fit exactly into the vector length
    ProcessBlock(param_inputs, mod_inputs, outputs, nFrames, NParams - (NParams % 4));
  }

  ParameterModulator<DynamicMods, StaticMods>& operator[](int idx)
  {
    return *(mParams[idx]);
  }

private:
  std::vector<ParameterModulator<DynamicMods, StaticMods>*> mParams;
  const int mParamBlocks;
  T mStaticScalar{ 0. };
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

  inline T Process(double freqHz) override
  {
    IOscillator<T>::SetFreqCPS(freqHz);

    return Process(mQNPos, mTransportIsRunning, mTempo);
  }

  inline T DoProcess() override
  {
    return Process(mQNPos, mTransportIsRunning, mTempo);
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
    return s_out * mLevelScalar;
  }

  static inline T WrapPhase(T x, T lo = 0., T hi = 1.)
  {
    while (x >= hi)
      x -= hi;
    while (x < lo)
      x += hi - lo;
    return x;
  };

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
class Sequencer final : public FastLFO<T>
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

  inline T DoProcess() override
  {
    return Process(FastLFO<T>::mQNPos, mTransportIsRunning, mTempo);
  }

  inline T Process(double qnPos = 0., bool transportIsRunning = false, double tempo = 120.) override
  {
    T oneOverQNScalar = 1. / LFO<T>::mQNScalar;
    T phase = IOscillator<T>::mPhase;

    if (FastLFO<T>::mRateMode == FastLFO<T>::ERateMode::kBPM && !transportIsRunning)
      IOscillator<T>::SetFreqCPS(tempo / 60.);

    double samplesPerBeat = IOscillator<T>::mSampleRate * (60.0 / (tempo == 0.0 ? 1.0 : tempo)); // samples per beat

    T phaseIncr = IOscillator<T>::mPhaseIncr;

    double sampleAccurateQnPos = qnPos + (phase / samplesPerBeat);

    if (FastLFO<T>::mRateMode == FastLFO<T>::ERateMode::kBPM)
    {
      phaseIncr *= LFO<T>::mQNScalar;
      if (transportIsRunning)
        phase = std::fmod(sampleAccurateQnPos, oneOverQNScalar) * LFO<T>::mQNScalar;
    }

    mStepPos = static_cast<int>(phase * NSteps) % mLength;
    T s_out = GetStepValue(mStepPos);
    IOscillator<T>::mPhase = FastLFO<T>::WrapPhase(IOscillator<T>::mPhase + phaseIncr, 0., static_cast<double>(mTableSize));
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

private:
  T* mStepValues;
  int mStepPos{ 0 };
  int mPrevStep{ -1 };
  int mLength{ NSteps };
};

/*
A container class that holds all per-voice modulators.
*/
template<typename T, template<typename> class EnvType=ADSREnvelope, template<typename> class LFOType=FastLFO, template<typename> class SequencerType=Sequencer>
class ModulatorList
{
public:
  ModulatorList() {}

  ModulatorList(EnvType<T>** envPtrs, LFOType<T>** lfoPtrs) : mEnvPtrs(envPtrs), mLFOPtrs(lfoPtrs)
  {}

  void AddModulator(EnvType<T>* env)
  {
    mEnvPtrs.push_back(env);
    mNumEnvs += 1;
  }

  void ReplaceModulator(EnvType<T>* env, size_t idx)
  {
    mEnvPtrs[idx] = env;
  }
  
  void AddModulator(FastLFO<T>* lfo)
  {
    mLFOPtrs.push_back(lfo);
    mNumLFOs += 1;
  }

  void AddModulator(Sequencer<T>* seq)
  {
    mSequencerPtrs.push_back(seq);
    mNumSeqs += 1;
  }

  /* Write modulation values to a buffer */
  void ProcessBlock(T** inputs, int nFrames)
  {
    for (auto i{ 0 }; i < nFrames; ++i)
    {
      for(auto e{0}; e < mEnvPtrs.size(); ++e)
        mModulators.GetList()[e][i] = mEnvPtrs[e]->Process(inputs[e][i]);
      for (auto l{ 0 }; l < mLFOPtrs.size(); ++l)
        mModulators.GetList()[l + mEnvPtrs.size()][i] = mLFOPtrs[l]->DoProcess();
      for (auto s{ 0 }; s < mSequencerPtrs.size(); ++s)
        mModulators.GetList()[s + mEnvPtrs.size() + mLFOPtrs.size()][i] = mSequencerPtrs[s]->DoProcess();
    }
  }

  inline T** GetList()
  {
    return mModulators.GetList();
  }

  void EmptyAndResize(int newLength, int numMods)
  {
    mNumMods = numMods;

    mModulatorRamps.Resize(newLength * mNumMods);
    mModulators.Empty();

    for (auto i = 0; i < mNumMods; i++)
    {
      mModulators.Add(mModulatorRamps.Get() + static_cast<size_t>(newLength) * i);
    }
  }

private:
  std::vector<EnvType<T>*> mEnvPtrs;
  std::vector<LFOType<T>*> mLFOPtrs;
  std::vector<SequencerType<T>*> mSequencerPtrs;
  int mNumMods{ 0 };
  int mNumEnvs{ 0 };
  int mNumLFOs{ 0 };
  int mNumSeqs{ 0 };

  WDL_PtrList<T> mModulators;
  WDL_TypedBuf<T> mModulatorRamps;
};

END_IPLUG_NAMESPACE