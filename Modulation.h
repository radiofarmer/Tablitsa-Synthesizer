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
    //memcpy(ParameterModulator::mModValues, modPtr, DynMods * sizeof(double));
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

  inline double AddModulation(double* modValues)
  {
    return AddModulation(mInitial, modValues);
  }

  inline double AddModulation_Vector(double* modValues)
  {
    return AddModulation_Vector(mInitial, modValues);
  }

  /* TODO: This may be better implemented with virtual functions (see ParameterModulationExp implementation below) */
  inline double AddModulation(double initVal, double* modValues) const
  {
    double modulation{ 0. };
    for (auto i{ 0 }; i < DynMods; ++i)
      modulation += mModDepths[i] * modValues[i];
    modulation += mStaticModulation;
    if (!mIsExponential)
      return std::max(std::min(initVal + modulation, mMax), mMin);
    else
      return std::max(std::min(initVal * std::exp(modulation), mMax), mMin);
  }

  inline double AddModulation_Vector(double initVal, double* modValues) const
  {
    constexpr int vectorsize = 4;

    double modulation{ 0. };
    Vec4d modDepths;
    Vec4d modValuesV;
    auto i{ 0 };
    for (; i < (DynMods & vectorsize); i+=4)
    {
      modDepths.load(mModDepths);
      modValuesV.load(modValues);
      modulation += horizontal_add(modDepths * modValuesV);
    }
    modDepths.load_partial(DynMods - i, mModDepths + i);
    modDepths.cutoff(DynMods - i);
    modValuesV.load_partial(DynMods - i, modValues + i);
    // modValues.cutoff(DynMods - i);
    modulation += horizontal_add(modDepths * modValuesV);
    modulation += mStaticModulation;
    // Clip and return
    if (!mIsExponential)
      return std::max(std::min(initVal + modulation, mMax), mMin);
    else
      return std::max(std::min(initVal * std::exp(modulation), mMax), mMin);
  }

  inline Vec4d __vectorcall AddModulation_Vector(Vec4d initVals, Vec4d* modVals)
  {
    Vec4d modulation = Vec4d(0.);
    for (auto i{ 0 }; i < DynMods; ++i)
      modulation += modVals[i] * mModDepths[i];
    modulation += mStaticModulation;
    if (!mIsExponential)
      return max(min(initVals + modulation, mMaxV), mMinV);
    else
      return max(min(initVals * exp(modulation), mMaxV), mMinV);
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
  double mModValues[DynMods + StatMods]{ 0. };
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
  inline void ProcessBlock(T** param_inputs, T** mod_inputs, T** outputs, int nFrames, const int offset = 0)
  {
    for (auto p{ 0 }; p < NParams; ++p)
    {
      for (auto i{ offset }; i < nFrames; ++i)
      {
        // Get the values of the modulators for the current frame
#ifndef VECTOR
        T modDepths[DynamicMods]{ 0. };
        for (auto m{ 0 }; m < DynamicMods; ++m)
          modDepths[m] = mod_inputs[m][i];
        outputs[p][i] = mParams[p]->AddModulation(param_inputs[p][i], modDepths);
#else
        T modDepths[DynamicMods]{ 0. };
        for (auto m{ 0 }; m < DynamicMods; ++m)
          modDepths[m] = mod_inputs[m][i];
        outputs[p][i] = mParams[p]->AddModulation_Vector(param_inputs[p][i], modDepths);
#endif
      }
    }
  }

  /* Vector processing of each block of parameter values.
  NB: The Clang compiler completely optimizes out the above `ProcessBlock` function, so implementing this function to provide further speed is probably unnecessary.
  */
  inline void ProcessBlockVec4d(T** param_inputs, T** mod_inputs, T** outputs, int nFrames)
  {
    constexpr int vectorsize = 4;

    for (auto p{ 0 }; p < NParams; ++p)
    {
      for (auto i{ 0 }; i < (nFrames & (-vectorsize)); i += vectorsize)
      {
        Vec4d modValues[DynamicMods];
        for (auto m{ 0 }; m < DynamicMods; ++m)
        {
          modValues[m] = Vec4d().load(&mod_inputs[m][i]);
        }
        Vec4d output = mParams[p]->AddModulation_Vector(Vec4d().load(&param_inputs[p][i]), modValues);
        output.store(&outputs[p][i]);
      }
    }
    // Process frames which don't fit into the vector size
    ProcessBlock(param_inputs, mod_inputs, outputs, nFrames, nFrames - (nFrames % 4));
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

  /* Write SAMPLE-ACCURATE meta-modulated modulation values to a buffer */
  inline void MetaProcessBlock_Accurate(int nFrames)
  {
    T metaModParams[2]{ 0. };
    T** ModList = mModulators.GetList(); // Pointer to first position of mModulators

    for (auto i{ 0 }; i < nFrames; ++i)
    {
      for (auto m{ 0 }; m < mModPtrs.size(); ++m)
      {
        for (auto mm{ 0 }; mm < mMetaParams[m].size(); ++mm)
        {
#if _DEBUG
          metaModParams[mm] = mMetaParams[m][mm]->AddModulation(mPrevValues);
#else
          // Possible place for optimization: Process multiple modulators simulataneously
          metaModParams[mm] = mMetaParams[m][mm]->AddModulation_Vector(mPrevValues);
#endif
        }
        mModPtrs[m]->SetParams(metaModParams);
        mPrevValues[m] = ModList[m][i] = mModPtrs[m]->Process();
      }
    }
  }

  /** Write meta-modulated modulation values to a buffer, using only the first calculated value for each
  modulator as a basis for modulation. */
  inline void MetaProcessBlock_Fast(T** inputs, int nFrames)
  {
    T metaModParams[2]{ 0. };
    T** modList = mModulators.GetList();

    // Loop through the first value of all the modulators
    for (auto m{ 0 }; m < mModPtrs.size(); ++m)
    {
      // Loop through the parameters that the current modulator requires and set their modulated values
      for (auto mm{ 0 }; mm < mMetaParams[m].size(); ++mm)
      {
#if _DEBUG
        metaModParams[mm] = mMetaParams[m][mm]->AddModulation(mPrevValues);
#else
        // Possible place for optimization: Process multiple modulators simulataneously
        metaModParams[mm] = mMetaParams[m][mm]->AddModulation_Vector(mPrevValues);
#endif
      }
      mModPtrs[m]->SetParams(metaModParams);
      mPrevValues[m] = modList[m][0] = mModPtrs[m]->Process();
    }
    for (auto m{ 0 }; m < mModPtrs.size(); ++m)
    {
      for (auto i{ 1 }; i < nFrames; ++i)
      {
        modList[m][i] = mModPtrs[m]->Process();
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
