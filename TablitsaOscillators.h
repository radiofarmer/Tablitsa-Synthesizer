#pragma once

#include "VectorFunctions.h"
#include "Oscillator.h"
#include "Wavetable.h"

#include <vectorclass.h>

//#if _DEBUG
#if 1
#define OVERSAMPLING 2
#define VECTOR_SIZE 4
#define OUTPUT_SIZE VECTOR_SIZE / OVERSAMPLING
#define VECTOR
#define RECURSION 1
#else
#define OVERSAMPLING 1
#define VECTOR_SIZE 1
#define OUTPUT_SIZE VECTOR_SIZE
#endif

BEGIN_IPLUG_NAMESPACE
template<typename T>
class VectorOscillator final : public FastSinOscillator<T>
{
  union tabfudge
  {
    double d;
    int i[2];
  } ALIGNED(8);

public:
  VectorOscillator(T startPhase=0.) : FastSinOscillator<T>(startPhase) {}
  
  inline Vec4d __vectorcall ProcessBlock4_Vector()
  {
    double phase = IOscillator<T>::mPhase + (double)UNITBIT32;
    const double phaseIncr = IOscillator<T>::mPhaseIncr * FastSinOscillator<T>::tableSize;
    Vec4d vPhase = phase + phaseIncr * mIncrVector;

    // Read 4 doubles as 8 integers (signed ints are used, but the AND operations later make this irrelevant)
    Vec8i viPhase = reinterpret_i(vPhase);
    // Upper 32 bits of 3*2^19 in ______ indices, 0xFFFF in _____: i.e. 0xFFFF, 0x18, 0xFFFF, ...
    Vec8i normhipart = blend8<HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V, 8>(reinterpret_i(Vec4d((double)UNITBIT32)), Vec8i(0xFFFFFFFF));
    // Mask the 8-item vector of 32-bit ints with one less than the table size, pad the upper bits (lower indices) with zeros, and reinterpret as a 4-item vector of 64-bit ints
    Vec8i offsets32 = viPhase & tableSizeM1;
    Vec4q offsets = reinterpret_i(permute8<HIOFFSET_V, -1, HIOFFSET_V + 2, -1, HIOFFSET_V + 4, -1, HIOFFSET_V + 6, -1>(offsets32));
    // Force the double to wrap. (AND the upper bits with the upper 32 bits of UNITBIT32)
    viPhase &= normhipart;
    Vec4d frac = reinterpret_d(viPhase) - (double)UNITBIT32; // get the fractional part
    const Vec4d f1 = lookup<tableSize>(offsets, mLUT);
    const Vec4d f2 = lookup<tableSize>(offsets, mLUT + 1);
    Vec4d output = mul_add((f2 - f1), frac, f1);

    // Restore mPhase
    tabfudge tf;
    phase += phaseIncr * 4.;
    tf.d = UNITBIT32 * tableSize;
    const int normhipart2 = tf.i[HIOFFSET];
    tf.d = phase + (UNITBIT32 * tableSize - UNITBIT32); // Remove the offset we introduced at the start of UNITBIT32.
    tf.i[HIOFFSET] = normhipart2;
    IOscillator<T>::mPhase = tf.d - UNITBIT32 * tableSize;
    return output;
  }

  inline Vec4d __vectorcall Process_Vector()
  {
//    double output[4];
//    FastSinOscillator<T>::ProcessBlock(output, 4);
//    return Vec4d().load_a(output);
    return ProcessBlock4_Vector();
  }

private:
  const Vec4d mIncrVector = Vec4d(0., 1., 2., 3.);
} ALIGNED(8);

END_IPLUG_NAMESPACE

template <typename T>
class WavetableOscillator final : public iplug::IOscillator<T>
{
  union tabfudge
  {
    double d;
    int i[2];
  } ALIGNED(8);

public:
  WavetableOscillator(const int id, const char* tableName, const double startPhase = 0., const double startFreq = 1.) :
    mID(id), IOscillator<T>(startPhase), mPrevFreq(static_cast<int>(startFreq))
  {
    WtFile table(tableName);
    LoadNewTable(table, mID);
    SetWavetable(WavetableOscillator<T>::LoadedTables[mID]);
    mWtReady = true;
  }

  WavetableOscillator(const int id, const WtFile& table, double startPhase = 0., double startFreq = 1.)
    : mID(id), IOscillator<T>(startPhase, startFreq), mPrevFreq(static_cast<int>(startFreq))
  {
    LoadNewTable(table, mID);
    SetWavetable(WavetableOscillator<T>::LoadedTables[mID]);
    mWtReady = true;
  }

  void SetSampleRate(double sampleRate)
  {
    mSampleRate = sampleRate * mProcessOS;
    mPhaseModulator.SetSampleRate(mSampleRate);
    mRingModulator.SetSampleRate(mSampleRate);
  }

  /* Load a new wavetable as a static variable */
  void LoadNewTable(WtFile& wt, int idx)
  {
    if (wt.Success())
    {
      std::unique_lock<std::mutex> lock(mMasterMutex);
      mWtReady = false;
      delete LoadedTables[idx];
      LoadedTables[idx] = new Wavetable<T>(wt);
    }
  }

  void SetWavetable(int idx)
  {
    SetWavetable(LoadedTables[idx]);
  }

  void SetWavetable(Wavetable<T>* tab)
  {
    std::unique_lock<std::mutex> lock(mWtMutex);
    mWtReady = false;
    if (tab != nullptr) // TODO: Check for nan's in the wavetable. They appear to break everything (including after new tables are loaded)
      mWT = tab;
    mPhaseIncrFactor = (1. / (mWT->mCyclesPerLevel * mProcessOS));
    mCyclesPerLevelRecip = 1. / mWT->mCyclesPerLevel;
  }

  // Chooses the proper mipmap for the current note frequency (Hz) and sample rate
  inline void SetMipmapLevel()
  {
    const double samplesPerCycle{ 1. / IOscillator<T>::mPhaseIncr };

    // Select by Index
    const double tableFact{ std::log2(mWT->GetMaxSize() / (samplesPerCycle * mTableOS)) };
    mTableInterp = tableFact - std::floor(tableFact);
    mWtIdx = std::max(static_cast<int>(std::floor(tableFact)), 0);
    SetMipmapLevel_ByIndex(mWtIdx);
  }

  inline void SetMipmapLevel_ByIndex(const int idx)
  {
    std::unique_lock<std::mutex> lock(mWtMutex);
    mCV.wait(lock, [this] { return mWtReady; });

    mWtPositionNorm = 1 - std::modf((1. - mWtPositionAbs) * (mWT->mNumTables - 1), &mWtOffset);
    int tableOffset = std::min(static_cast<int>(mWtOffset), mWT->mNumTables - 2);
    mLUTLo[0] = mWT->GetMipmapLevel_ByIndex(tableOffset, idx, mTableSize);
    mLUTLo[1] = mWT->GetMipmapLevel_ByIndex(tableOffset + 1, idx, mTableSize);
    mLUTHi[0] = mWT->GetMipmapLevel_ByIndex(tableOffset, idx + 1, mNextTableSize);
    mLUTHi[1] = mWT->GetMipmapLevel_ByIndex(tableOffset + 1, idx + 1, mNextTableSize);
    mTableSizeM1 = mTableSize - 1;
    mNextTableSizeM1 = mNextTableSize - 1;
  }

  // Chooses the proper mipmap for a particular frequency (Hz)
  inline void SetMipmapLevel(int size)
  {
    std::unique_lock<std::mutex> lock(mWtMutex);
    mCV.wait(lock, [this] { return TableLoaded; });

    assert(mWT != nullptr);
    size = std::min(size, mWT->GetMaxSize());

    mWtPositionNorm = 1. - std::modf((1. - mWtPositionAbs) * (mWT->mNumTables - 1), &mWtOffset);
    int tableOffset = static_cast<int>(mWtOffset);
    mLUTLo[0] = mWT->GetMipmapLevel_BySize(tableOffset, size, mTableSize);
    mLUTLo[1] = mWT->GetMipmapLevel_BySize(tableOffset + 1, size, mNextTableSize);
    mTableSizeM1 = mTableSize - 1;
    mNextTableSizeM1 = mNextTableSize - 1;
    // Higher-level (smaller) mipmap locations
    mLUTHi[0] = mLUTLo[0] + mTableSize;
    mLUTHi[1] = mLUTLo[1] + mTableSize;
  }

  //todo rewrite this
  inline T Process()
  {
    T output = 0.;
    ProcessBlock(&output, 1);

    return output;
  }

  inline void AdjustWavetable(double freqCPS)
  {
    IOscillator<T>::SetFreqCPS(freqCPS);

    if (mPrevFreq != static_cast<int>(freqCPS))
    {
      SetMipmapLevel();
      mPrevFreq = static_cast<int>(freqCPS);
    }
    else
    {
      SetMipmapLevel_ByIndex(mWtIdx);
    }
  }

  inline std::array<T, OUTPUT_SIZE> ProcessMultiple(double freqCPS)
  {
    AdjustWavetable(freqCPS);

    std::array<T, OUTPUT_SIZE> output{ 0. };
#if VECTOR_SIZE == 4
    Vec4d osc_out = ProcessOversamplingVec4();
    osc_out.store_partial(2, &output[0]);
#else
    ProcessOversampling(output, mProcessOS);
#endif
    return output;
  }

  inline T Process(double freqHz)
  {
    return 0.;
  }

  inline void ProcessOversampling(std::array<T, OUTPUT_SIZE>& pOutput, int nFrames)
  {

    double phase; // integer phase
    double frac, frac2; // fractional phases

    // tableOffset -= std::max(floor(tableOffset - 1e-7), 0.);
#if OVERSAMPLING > 1
    // Temporary array to be downsampled
    T* oversampled = new T[nFrames]{ 0. };
#endif
    for (auto s = 0; s < nFrames; ++s)
    {
      double phaseShift = WrapPhase(SamplePhaseShift(IOscillator<T>::mPhase) + PhaseMod()); // Adjusted phase (double between zero and one)
      frac = modf(phaseShift * mTableSize, &phase); // Split phase into integer and fractional parts
      const int offset{ static_cast<int>(phase) };
      // Higher-frequency wavetable offset
      frac2 = modf(phaseShift * mNextTableSize, &phase);
      const int halfOffset{ static_cast<int>(phase) };
      // Get two samples from each waveform
      const T* addr[4] = {
        mLUTLo[0] + offset,
        mLUTLo[1] + (offset & mTableSizeM1),
        mLUTHi[0] + halfOffset,
        mLUTHi[1] + (halfOffset & mNextTableSizeM1)
      }; // Obtain the integer portion
      // Read from wavetable
      const double sampleWtPosition{ mWtPositionNorm };
      const double sampleWtPositionInv{ 1 - sampleWtPosition };
      const T f1 = addr[0][0] * sampleWtPosition + addr[1][0] * sampleWtPositionInv;
      const T f2 = addr[0][1] * sampleWtPosition + addr[1][1] * sampleWtPositionInv;
      const T f3 = addr[2][0] * sampleWtPosition + addr[3][0] * sampleWtPositionInv;
      const T f4 = addr[2][1] * sampleWtPosition + addr[3][1] * sampleWtPositionInv;
      
      T lowerTable = (f1 + frac * (f2 - f1));
#if OVERSAMPLING > 1
      oversampled[s] = lowerTable + mTableInterp * ((f3 + frac2 * (f4 - f3)) - lowerTable);
      oversampled[s] += mRM * mRingModAmt * (RingMod() * oversampled[s] - oversampled[s]);
#else
      const T output = lowerTable + mTableInterp * ((f3 + frac2 * (f4 - f3)) - lowerTable);
      pOutput[s] = output + mRM * mRingModAmt * (RingMod() * output - output);
      mLastOutput = pOutput[s];
#endif
      // Increment Phase
      IOscillator<T>::mPhase += mPhaseIncr * mPhaseIncrFactor;
      IOscillator<T>::mPhase -= std::floor(IOscillator<T>::mPhase);
    }
#if OVERSAMPLING > 1
    for (auto s = 0; s < nFrames / mProcessOS; ++s)
    {
      mLastOutput = pOutput[s] = mAAFilter.ProcessAndDownsample(oversampled + (s * mProcessOS));
    }
    delete[] oversampled;
#endif

  }

#define VECTOR_TEST
  inline Vec4d __vectorcall ProcessOversamplingVec4()
  {
    double tableOffset{ mWtPositionAbs * (mWT->mNumTables - 1) };
    tableOffset -= std::max(floor(tableOffset - 0.0001), 0.);
#ifdef VECTOR_TEST
    double phaseNorm = SamplePhaseShift(mPhase / mTableSize); // for phase shift

    double phase = phaseNorm * mTableSize + (double)UNITBIT32;
    const double phaseIncr = mPhaseIncr * mPhaseIncrFactor * mProcessOS * mTableSize;
    const double phaseIncr2 = mPhaseIncr * mPhaseIncrFactor * mProcessOS * mNextTableSize;
    Vec4d vPhase = phase + phaseIncr * mIncrVec;

    // Read 4 doubles as 8 integers (signed ints are used, but the AND operations later make this irrelevant)
    Vec8i viPhase = reinterpret_i(vPhase);
    // Upper 32 bits of 3*2^19 in ______ indices, 0xFFFF in _____: i.e. 0xFFFF, 0x18, 0xFFFF, ...
    Vec8i normhipart = blend8<8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V>(reinterpret_i(Vec4d((double)UNITBIT32)), Vec8i(0));//Vec8i(0xFFFFFFFF));
    // Mask the 8-item vector of 32-bit ints with one less than the table size, pad the upper bits (lower indices) with zeros, and reinterpret as a 4-item vector of 64-bit ints
    Vec8i offsets32 = viPhase & mTableSizeM1;
    Vec4q offsets = reinterpret_i(permute8<HIOFFSET_V, -1, HIOFFSET_V + 2, -1, HIOFFSET_V + 4, -1, HIOFFSET_V + 6, -1>(offsets32));
    offsets = (offsets + to_int64_in_range(PhaseMod(), mTableSize)) & mTableSizeM1;
    Vec4q offsets1 = (offsets + 1); // &mTableSizeM1;
    // Force the double to wrap. (AND the upper bits with the upper 32 bits of UNITBIT32)
    viPhase &= normhipart;
    Vec4d phaseFrac = reinterpret_d(viPhase) - (double)UNITBIT32; // get the fractional part

    // Larger/lower table
    constexpr int maxTableSize = 16384 * 12;
    const Vec4d tb0s0lo = lookup<maxTableSize>(offsets, mLUTLo[0]);
    const Vec4d tb0s1lo = lookup<maxTableSize>(offsets1, mLUTLo[0]);
    const Vec4d tb1s1lo = lookup<maxTableSize>(offsets1, mLUTLo[1]);
    const Vec4d tb1s0lo = lookup<maxTableSize>(offsets, mLUTLo[1]);
    const Vec4d tb0lo = mul_add((tb0s1lo - tb0s0lo), phaseFrac, tb0s0lo);
    const Vec4d tb1lo = mul_add((tb1s1lo - tb1s0lo), phaseFrac, tb1s0lo);

    // Recalculate offsets for the smaller table
    vPhase = phaseNorm * mNextTableSize + (double)UNITBIT32 + phaseIncr2 * mIncrVec;

    // Read 4 doubles as 8 integers (signed ints are used, but the AND operations later make this irrelevant)
    viPhase = reinterpret_i(vPhase);
    normhipart = blend8<8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V>(reinterpret_i(Vec4d((double)UNITBIT32)), Vec8i(0));
    // Mask the 8-item vector of 32-bit ints with one less than the table size, pad the upper bits (lower indices) with zeros, and reinterpret as a 4-item vector of 64-bit ints
    offsets32 = viPhase & mNextTableSizeM1;
    offsets = reinterpret_i(permute8<HIOFFSET_V, -1, HIOFFSET_V + 2, -1, HIOFFSET_V + 4, -1, HIOFFSET_V + 6, -1>(offsets32));
    offsets = (offsets + to_int64_in_range(PhaseMod(), mTableSize)) & mNextTableSizeM1;
    offsets1 = (offsets + 1); // & mNextTableSizeM1;
    // Force the double to wrap. (AND the upper bits with the upper 32 bits of UNITBIT32)
    viPhase &= normhipart;
    phaseFrac = reinterpret_d(viPhase) - (double)UNITBIT32; // get the fractional part

    // Smaller/higher table
    const Vec4d tb0s0hi = lookup<maxTableSize>(offsets, mLUTHi[0]);
    const Vec4d tb0s1hi = lookup<maxTableSize>(offsets1, mLUTHi[0]);
    const Vec4d tb1s0hi = lookup<maxTableSize>(offsets, mLUTHi[1]);
    const Vec4d tb1s1hi = lookup<maxTableSize>(offsets1, mLUTHi[1]);
    // Interpolate samples
    const Vec4d tb0hi = mul_add(tb0s1hi - tb0s0hi, phaseFrac, tb0s0hi);
    const Vec4d tb1hi = mul_add(tb1s1hi - tb1s0hi, phaseFrac, tb1s0hi);

    // Interpolate mipmap levels
    Vec4d tb0 = mul_add(tb0hi - tb0lo, mTableInterp, tb0lo);
    Vec4d tb1 = mul_add(tb1hi - tb1lo, mTableInterp, tb1lo);

    // Restore mPhase
    tabfudge tf;
    phase += phaseIncr * (double)VECTOR_SIZE;
    tf.d = UNITBIT32 * mTableSize;
    const int normhipart2 = tf.i[HIOFFSET];
    tf.d = phase + (UNITBIT32 * mTableSize - UNITBIT32); // Remove the offset we introduced at the start of UNITBIT32.
    tf.i[HIOFFSET] = normhipart2; // Reset the upper 32 bits
    IOscillator<T>::mPhase = tf.d - UNITBIT32 * mTableSize;

    // Mix wavtables
    Vec4d mixed = mul_add(tb1 - tb0, 1 - tableOffset, tb0);
    Vec4d ringMod = RingMod();
    mixed = mul_add(ringMod - 1., mRM * mRingModAmt * mixed, mixed);
#else
    const double phaseIncr = mPhaseIncr * mPhaseIncrFactor * mProcessOS;
    Vec4d phaseMod = PhaseMod(); //Vec4d(PhaseMod(), PhaseMod(), PhaseMod(), PhaseMod());
    Vec4d phase = phaseMod + SamplePhaseShift(IOscillator<T>::mPhase);
    Vec4d phaseDouble = mul_add(mIncrVec, phaseIncr, phase) * mTableSize; // Next four phase positions in samples, including fractional position
    Vec4q phaseInt = truncatei(phaseDouble); // Indices of the (left) samples to read
    Vec4d phaseFrac = phaseDouble - to_double(phaseInt);
    phaseInt &= mTableSizeM1;
    Vec4q phaseInt1 = (phaseInt + 1) & mTableSizeM1;

    // Read from wavetables (lower/larger table)
    Vec4d tb0s0lo = lookup<16384 * 12>(phaseInt, mLUTLo[0]);
    Vec4d tb1s0lo = lookup<16384 * 12>(phaseInt, mLUTLo[1]);
    Vec4d tb0s1lo = lookup<16384 * 12>(phaseInt1, mLUTLo[0]);
    Vec4d tb1s1lo = lookup<16384 * 12>(phaseInt1, mLUTLo[1]);
    // Interpolate samples
    Vec4d tb0lo = mul_add((tb0s1lo - tb0s0lo), phaseFrac, tb0s0lo);
    Vec4d tb1lo = mul_add((tb1s1lo - tb1s0lo), phaseFrac, tb1s0lo);

    // Calculate indices for the higher-frequency table
    phaseDouble = mul_add(mIncrVec, phaseIncr, phase) * mNextTableSize;
    phaseInt = truncatei(phaseDouble); // Poor performance if not AVX512DQ
    phaseFrac = phaseDouble - to_double(phaseInt); // can also use truncate(phaseDouble) to get a floored double
    phaseInt &= mNextTableSizeM1;
    phaseInt1 = (phaseInt + 1) & mNextTableSizeM1;

    // Read from wavetables (higher/smaller table)
    // `lookup<n>` function requires Vec4q for Vec4d or Vec8q for Vec8d
    Vec4d tb0s0hi = lookup<16384 * 12>(phaseInt, mLUTHi[0]);
    Vec4d tb0s1hi = lookup<16384 * 12>(phaseInt1, mLUTHi[0]);
    Vec4d tb1s0hi = lookup<16384 * 12>(phaseInt, mLUTHi[1]);
    Vec4d tb1s1hi = lookup<16384 * 12>(phaseInt1, mLUTHi[1]);
    // Interpolate samples
    Vec4d tb0hi = mul_add(tb0s1hi - tb0s0hi, phaseFrac, tb0s0hi);
    Vec4d tb1hi = mul_add(tb1s1hi - tb1s0hi, phaseFrac, tb1s0hi);

    // Interpolate mipmap levels
    Vec4d tb0 = mul_add(tb0hi - tb0lo, mTableInterp, tb0lo);
    Vec4d tb1 = mul_add(tb1hi - tb1lo, mTableInterp, tb1lo);

    // Mix wavetables
    Vec4d ringMod = RingMod();
    Vec4d mixed = mul_add(tb1 - tb0, 1 - tableOffset, tb0);
    mixed = mul_add(ringMod - 1., mRM * mRingModAmt * mixed, mixed);

    IOscillator<T>::mPhase += phaseIncr * (double)VECTOR_SIZE;
    IOscillator<T>::mPhase -= floor(IOscillator<T>::mPhase);
#endif

#if !RECURSION
    T oversampled[VECTOR_SIZE];
    mixed.store(oversampled);
    for (auto s = 0; s < VECTOR_SIZE / mProcessOS; ++s)
    {
      pOutput[s] = mAAFilter.ProcessAndDownsample(oversampled + (s * mProcessOS));
    }
#else
    // Using recursive function
    return mAAFilter.ProcessAndDownsample_Recursive(mixed);
#endif

  }

  /* Returns a new phase as a power function of the current phase. Also sets mPhaseInCycle to the unshifted value. */
  inline double SamplePhaseShift(double phase)
  {
    double cycle;
    mPhaseInCycle = modf(phase * mWT->mCyclesPerLevel, &cycle);
    return (cycle + std::pow(mPhaseInCycle, 1. + (mWtBend >= 0. ? mWtBend : mWtBend / 2.))) * mCyclesPerLevelRecip;
  }

  inline double WrapPhase(double phase)
  {
    while (phase < 0.)
      phase++;
    while (phase >= 1.)
      phase--;
    return phase;
  }

  /* Returns an adjusted phase increment based on the current (cycle-normalized) phase. */
#ifndef VECTOR
  inline double PhaseMod()
  {
    return mPM * mPhaseModAmt * mPhaseModulator.Process() * mCyclesPerLevelRecip;
  }

  inline double RingMod()
  {
    return mRingModulator.Process();
  }
#else
  inline Vec4d __vectorcall RingMod()
  {
    return mRingModulator.Process_Vector();
  }

  // TODO: use enable_if to choose between vector lengths
  inline Vec4d __vectorcall PhaseMod()
  {
    return mPM * mPhaseModAmt * mPhaseModulator.Process_Vector() * mCyclesPerLevelRecip;
  }
#endif
  inline double* GetWtPosition()
  {
    return &mWtPositionAbs;
  }

  inline double* GetWtBend()
  {
    return &mWtBend;
  }

  inline void SetWtPosition(double wtPos)
  {
    mWtPositionAbs = std::min(wtPos, 0.999);
  }

  inline void SetWtBend(double wtBend)
  {
    mWtBend = wtBend;
  }

  inline double SampleWavetablePosition(double phase)
  {
    return mWtPositionAbs;
  }

  inline double GetSampleRate()
  {
    return mSampleRate;
  }

  inline void SetPhaseModulation(bool on)
  {
    if (on)
      mPM = 1.;
    else
      mPM = 0.;
  }
  inline void SetPhaseModulation(double amt, double freqCPS)
  {
    mPhaseModAmt = amt;
    mPhaseModulator.SetFreqCPS(freqCPS);
  }

  inline void SetRingModulation(bool on)
  {
    if (on)
      mRM = 1.;
    else
      mRM = 0.;
  }
  inline void SetRingModulation(double amt, double freqCPS)
  {
    mRingModAmt = amt;
    mRingModulator.SetFreqCPS(freqCPS);
  }

  void ReloadLUT()
  {
    mPrevFreq = -1;
  }

  void NotifyLoaded()
  {
    mWtReady = true;
    mCV.notify_all();
  }

  T mLastOutput = 0.;

private:
#if OVERSAMPLING == 2
  static constexpr int mProcessOS{ 2 }; // Sample processing oversampling level (number of samples processed per sample output)
#else
  static constexpr int mProcessOS{ 1 };
#endif
  ChebyshevBL<T> mAAFilter;

  // TODO: Order elements mindful of cache access:

  // Oscillator ID
  int mID;

  // Lookup Table Parameters
  int mTableSize = WT_SIZE; // Default: 2^9
  int mTableSizeM1 = WT_SIZE - 1; // Default: 2^9 -1
  int mNextTableSize = WT_SIZE / 2;
  int mNextTableSizeM1 = WT_SIZE / 2 - 1;
  int mTableOS{ 8 }; // Wavetable oversampling level (ratio of table size to maximum samples per cycle read from the table)
  double mTableInterp{ 1 };
  double mPhaseIncrFactor{ 1. }; // Reciprocal of the product of the number of cycles per wavetable level and the processing oversampling level
  double mPhaseInCycle{ 0. }; // The fractional phase (between zero and one) within a single cycle, independent of the number of cycles per wavetable level.
  double mCyclesPerLevelRecip{ 1. };

  // Wavetabe Data
  T* mLUTLo[2]{}; // Pointer for lower-frequency (higher-sample) lookup table
  T* mLUTHi[2]{};
  Wavetable<T>* mWT{ nullptr };

  // Wavetable Timbre Parameters
  double mWtPositionAbs{ 1. }; // Absolute position in full set of wavetables. Range: [0., NTables - 1]
  double mWtPositionNorm{ 1. }; // Normalized position between two wavetables. Range [0., 1.].
  double mWtOffset{ 0. };
  double mWtSpacing{ 1. };
  double mWtBend{ 0 };
  int mPrevFreq;
  int mWtIdx{ 0 }; // Mipmap level

  // Thread-related members for wavetable updates
  static inline std::mutex mWtMutex; // Mutex used when swapping out the current wavetable in each oscillator object
  static inline std::condition_variable mCV;
  bool TableLoaded{ false };
  static inline std::mutex mMasterMutex; // Static mutex used when loading a new wavetable from a file
  bool mWtReady{ false };

  // Vectors
#if VECTOR_SIZE == 4
  const Vec4d mIncrVec{ 0., 1., 2., 3. };
#else
  const Vec8d mIncrVec{ 0., 1., 2., 3., 4., 5., 6., 7. };
#endif

  // Modulation
  double mPM{ 0. };
  double mPhaseModAmt{ 0. };
  double mPhaseModFreq;
  double mRM{ 0. };
  double mRingModAmt{ 0. };
  double mRingModFreq;

  static inline constexpr double twoPi{ 6.28318530718 };

public:
  Wavetable<T>* LoadedTables[2]{ nullptr, nullptr };
#ifndef VECTOR
  iplug::FastSinOscillator<T> mPhaseModulator;
  iplug::FastSinOscillator<T> mRingModulator{ 0.5 }; // Offset start phase by half a cycle
#else
  iplug::VectorOscillator<T> mPhaseModulator;
  iplug::VectorOscillator<T> mRingModulator{ 0.5 }; // Offset start phase by half a cycle
#endif
};