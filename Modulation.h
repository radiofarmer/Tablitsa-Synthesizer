#pragma once

#include "Modulators.h"

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

  inline double AddModulation_Vector(double initVal)
  {
    double modulation{ 0. };
    Vec4d modDepths;
    Vec4d modValues;
    constexpr int vectorsize = 4;
    auto i{ 0 };
    for (; i < (DynMods & vectorsize); ++i)
    {
      modDepths.load(mModDepths);
      modValues.load(ParameterModulator::mModValues);
      modulation += horizontal_add(modDepths * modValues);
    }
    modDepths.load_partial(DynMods - i, mModDepths + i);
    modDepths.cutoff(DynMods - i);
    modValues.load_partial(DynMods - i, ParameterModulator::mModValues + i);
    // modValues.cutoff(DynMods - i);
    modulation += horizontal_add(modDepths * modValues);
    modulation += mStaticModulation;
    // Clip and return
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

template<typename T, int NParams = 1, int DynamicMods = 6, int StaticMods = 3, int VectorSize = 4>
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
  void ProcessBlock(T** param_inputs, T** mod_inputs, T** outputs, int nFrames, const int offset = 0)
  {
    for (auto i{ 0 }; i < nFrames; ++i)
    {
      T modDepths[NParams]{ 0. };
      for (auto m{ 0 }; m < DynamicMods; ++m)
        modDepths[m] = mod_inputs[m][i];
      ParameterModulator<DynamicMods, StaticMods>::SetModValues(modDepths);
      for (auto p{ offset }; p < NParams; ++p)
      {
#if _DEBUG
        outputs[p][i] = mParams[p]->AddModulation(param_inputs[p][i]);
#else
        outputs[p][i] = mParams[p]->AddModulation_Vector(param_inputs[p][i]);
#endif
      }
    }
  }

  /* Vector processing of each block of parameter values.

  TODO: Pad the modulation matrix with extra rows of zeros so that its rows are always divisible by the vector length. (See ModulatorList)
  TODO: Support for exponentially-scaled moduation
  */
  void ProcessBlockVec4d(T** param_inputs, T** mod_inputs, T** outputs, int nFrames)
  {
    for (auto i{ 0 }; i < nFrames; i += 4)
    {
      for (auto p{ 0 }; p < NParams; ++p)
      {
        // For calculating per-sample modulation for single samples in parallel
        Vec4d modDepths;
        modDepths.load(mParams[p]->Depths()); // Modulation depths for the current parameters (constant for entire sample block)
        Vec4d modVals(mod_inputs[0][i], mod_inputs[1][i], mod_inputs[2][i], mod_inputs[3][i]); // Values of the modulators themselves for the current samples
        outputs[p][i] = mParams[p]->ClipToRange(horizontal_add(modDepths * modVals) + param_inputs[p][i]);

        // For calculating per-sample modulation for multiple samples in parallel
        /*Vec4d modSum(0.);
        Vec4d initVals;
        initVals.load(param_inputs[p]);
        for (auto m{ 0 }; m < 4; ++m)
        {
          Vec4d modSamples;
          modSamples.load(&mod_inputs[m][i]);
          modSum = mul_add(modSamples, mParams.at(p)[m], modSum);
        }
        Vec4d outputV = Params[p]->ClipToRange(initVals + modSum);
        outputV.store(&outputs[p][i]);*/

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


/*
A container class that holds all per-voice modulators.
*/
template<typename T, class EnvType, class LFOType, class SequencerType, int DynMods, int StatMods>
class ModulatorList
{
public:
  ModulatorList() {}

  ModulatorList(EnvType** envPtrs, LFOType** lfoPtrs) : mEnvPtrs(envPtrs), mLFOPtrs(lfoPtrs)
  {}

  void AddModulator(EnvType* env)
  {
    /* TODO: Ideally, the `mModPtrs` vector should be reconstructed from scratch after the addition of a new modulator,
    * so that the order is guranteed to be Envelopes-LFOs-Sequencers. */

    mModPtrs.push_back(env);

    mEnvPtrs.push_back(env);
    mNumEnvs += 1;
    mNumMods += 1;
    ResizeHistoryBuffer();
  }

  void AddModulator(EnvType* env, std::vector<ParameterModulator<DynMods, StatMods>*>& metaParams)
  {
    AddModulator(env);
    mMetaParams.push_back(metaParams);
  }

  void AddModulator(LFOType* lfo)
  {
    mModPtrs.push_back(lfo);

    mLFOPtrs.push_back(lfo);
    mNumLFOs += 1;
    mNumMods += 1;
    ResizeHistoryBuffer();
  }

  void AddModulator(LFOType* lfo, std::vector<ParameterModulator<DynMods, StatMods>*>& metaParams)
  {
    AddModulator(lfo);
    mMetaParams.push_back(metaParams);
  }

  void AddModulator(SequencerType* seq)
  {
    mModPtrs.push_back(seq);

    mSequencerPtrs.push_back(seq);
    mNumSeqs += 1;
    mNumMods += 1;
    ResizeHistoryBuffer();
  }

  void AddModulator(SequencerType* seq, std::vector<ParameterModulator<DynMods, StatMods>*>& metaParams)
  {
    AddModulator(seq);
    mMetaParams.push_back(metaParams);
  }

  void ReplaceModulator(EnvType* env, size_t idx)
  {
    mModPtrs[idx] = env;
    mEnvPtrs[idx] = env;
  }

  void ReplaceModulator(LFOType* lfo, size_t idx)
  {
    mModPtrs[idx + mNumEnvs] = lfo;
    mLFOPtrs[idx] = lfo;
  }

  void ReplaceModulator(SequencerType* seq, size_t idx)
  {
    mModPtrs[idx + mNumEnvs + mNumLFOs] = seq;
    mSequencerPtrs[idx] = seq;
  }

  /* Write modulation values to a buffer */
  void ProcessBlock(T** inputs, int nFrames)
  {
    T** modList = mModulators.GetList();
    for (auto i{ 0 }; i < nFrames; ++i)
    {
      for (auto e{ 0 }; e < mEnvPtrs.size(); ++e)
        modList[e][i] = mEnvPtrs[e]->Process(inputs[e][i]);
      for (auto l{ 0 }; l < mLFOPtrs.size(); ++l)
        modList[l + mEnvPtrs.size()][i] = mLFOPtrs[l]->Process();
      for (auto s{ 0 }; s < mSequencerPtrs.size(); ++s)
        modList[s + mEnvPtrs.size() + mLFOPtrs.size()][i] = mSequencerPtrs[s]->Process();
    }
  }

  /* Write meta-modulated modulation values to a buffer */
  void MetaProcessBlock(T** inputs, int nFrames)
  {
    T metaModParams[2]{ 0. };
    T** modList = mModulators.GetList();
    for (auto i{ 0 }; i < nFrames; ++i)
    {
      ParameterModulator<DynMods, StatMods>::SetModValues(mPrevValues);
      for (auto m{ 0 }; m < mModPtrs.size(); ++m)
      {
        for (auto mm{ 0 }; mm < mMetaParams[m].size(); ++mm)
        {
          metaModParams[mm] = mMetaParams[m][mm]->AddModulation();
        }
        mModPtrs[m]->SetParams(metaModParams);
        mPrevValues[m] = modList[m][i] = mModPtrs[m]->Process();
      }
    }
  }

  inline T** GetList()
  {
    return mModulators.GetList();
  }

  void EmptyAndResize(int newLength)
  {
    EmptyAndResize(newLength, mNumMods);
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

  void ResizeHistoryBuffer()
  {
    if (mPrevValues)
      delete[] mPrevValues;
    mPrevValues = new T[mNumMods];
  }

private:
  std::vector<EnvType*> mEnvPtrs;
  std::vector<LFOType*> mLFOPtrs;
  std::vector<SequencerType*> mSequencerPtrs;
  std::vector<GenericModulator<T>*> mModPtrs;
  std::vector<std::vector<ParameterModulator<DynMods, StatMods>*>> mMetaParams;
  int mNumMods{ 0 };
  int mNumEnvs{ 0 };
  int mNumLFOs{ 0 };
  int mNumSeqs{ 0 };
  T* mPrevValues{ nullptr };

  WDL_PtrList<T> mModulators;
  WDL_TypedBuf<T> mModulatorRamps;
};
