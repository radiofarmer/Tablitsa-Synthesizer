#pragma once

#include "Wavetable.h"

#include <vectorclass.h>

#if !_DEBUG
#define OVERSAMPLING 2
#define VECTOR_SIZE 4
#define OUTPUT_SIZE VECTOR_SIZE / OVERSAMPLING
#define VECTOR
#else
#define OVERSAMPLING 1
#define VECTOR_SIZE 1
#define OUTPUT_SIZE VECTOR_SIZE
#endif

template<typename T>
class VectorOscillator final : public iplug::IOscillator<T>
{
public:
  VectorOscillator(T startPhase) : IOscillator<T>(startPhase) {}

  void ProcessBlockVector(T* pOutput, int nFrames)
  {
    double phase = IOscillator<T>::mPhase + (double)UNITBIT32;
    const Vec4d phaseIncr = IOscillator<T>::mPhaseIncr * mIncrVector * tableSize;

    union tabfudge tf;
    tf.d = UNITBIT32;
    const int normhipart = tf.i[HIOFFSET];

    for (auto s = 0; s < nFrames; s++)
    {
      tf.d = phase;
      phase += phaseIncr;
      const T* addr = mLUT + (tf.i[HIOFFSET] & tableSizeM1); // Obtain the integer portion
      tf.i[HIOFFSET] = normhipart; // Force the double to wrap.
      const double frac = tf.d - UNITBIT32;
      const T f1 = addr[0];
      const T f2 = addr[1];
      mLastOutput = pOutput[s] = T(f1 + frac * (f2 - f1));
    }

    // Restore mPhase
    tf.d = UNITBIT32 * tableSize;
    const int normhipart2 = tf.i[HIOFFSET];
    tf.d = phase + (UNITBIT32 * tableSize - UNITBIT32); // Remove the offset we introduced at the start of UNITBIT32.
    tf.i[HIOFFSET] = normhipart2;
    IOscillator<T>::mPhase = tf.d - UNITBIT32 * tableSize;
  }

private:
  const Vec4d mIncrVector = Vec4d(0., 1., 2., 3.);
};

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
    WavetableOscillator<T>::LoadNewTable(table, mID);
    WavetableOscillator<T>::SetWavetable(WavetableOscillator<T>::LoadedTables[mID]);
  }

  WavetableOscillator(const int id, const WtFile& table, double startPhase = 0., double startFreq = 1.)
    : mID(id), IOscillator<T>(startPhase, startFreq), mPrevFreq(static_cast<int>(startFreq))
  {
    WavetableOscillator<T>::LoadNewTable(table, mID);
    WavetableOscillator<T>::SetWavetable(WavetableOscillator<T>::LoadedTables[mID]);
  }

  void SetSampleRate(double sampleRate)
  {
    mSampleRate = sampleRate * mProcessOS;
    mPhaseModulator.SetSampleRate(mSampleRate);
    mRingModulator.SetSampleRate(mSampleRate);
  }

  /* Load a new wavetable as a static variable */
  static void LoadNewTable(WtFile& wt, int idx)
  {
    if (wt.Success())
    {
      std::unique_lock<std::mutex> lock(mMasterMutex);
      mWtReady[idx] = false;
      delete LoadedTables[idx];
      LoadedTables[idx] = new Wavetable<T>(wt);
    }
  }

  void SetWavetable(Wavetable<T>* tab)
  {
    std::unique_lock<std::mutex> lock(mWtMutex);
    mWtReady[mID] = false;
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
    mCV.wait(lock, [this] { return mWtReady[mID]; });

    int tableOffset = static_cast<int>((1 - mWtPosition) * (mWT->mNumTables - 1.0001));
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

    int tableOffset = static_cast<int>((1 - mWtPosition) * (mWT->mNumTables - 1.0001));
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
#if VECTOR_SIZE == 8
    ProcessOversamplingVec<Vec8d, Vec8q>(output);
#elif VECTOR_SIZE == 4
    ProcessOversamplingVec4(output);
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

    // Get the normalized offset of the current wavetable block
    double tableOffset{ mWtPosition * (mWT->mNumTables - 1) };
    tableOffset -= std::max(floor(tableOffset - 0.0001), 0.);
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
      const double sampleWtPosition{ tableOffset };
      const double sampleWtPositionInv{ 1 - sampleWtPosition };
      const T f1 = addr[0][0] * sampleWtPosition + addr[1][0] * sampleWtPositionInv;
      const T f2 = addr[0][1] * sampleWtPosition + addr[1][1] * sampleWtPositionInv;
      const T f3 = addr[2][0] * sampleWtPosition + addr[3][0] * sampleWtPositionInv;
      const T f4 = addr[2][1] * sampleWtPosition + addr[3][1] * sampleWtPositionInv;
      // Send output
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

  inline void ProcessOversamplingVec4(std::array<T, OUTPUT_SIZE>& pOutput)
  {
    double tableOffset{ mWtPosition * (mWT->mNumTables - 1) };
    tableOffset -= std::max(floor(tableOffset - 0.0001), 0.);

    const double phaseIncr = mPhaseIncr * mPhaseIncrFactor * mProcessOS;
    Vec4d phaseMod = Vec4d(PhaseMod(), PhaseMod(), PhaseMod(), PhaseMod());
    Vec4d phase = phaseMod + SamplePhaseShift(IOscillator<T>::mPhase);
    Vec4d phaseDouble = mul_add(mIncrVec, phaseIncr, phase) * mTableSize; // Next four phase positions in samples, including fractional position
    Vec4q phaseInt = truncatei(phaseDouble); // Indices of the (left) samples to read
    Vec4d phaseFrac = phaseDouble - to_double(phaseInt);
    phaseInt &= mTableSizeM1;
    Vec4q phaseInt1 = (phaseInt + 1) & mTableSizeM1;

    // Read from wavetables (lower/larger table)
    Vec4d tb0s0lo = lookup<16384 * 12>(phaseInt, mLUTLo[0]);
    Vec4d tb0s1lo = lookup<16384 * 12>(phaseInt1, mLUTLo[0]);
    Vec4d tb1s0lo = lookup<16384 * 12>(phaseInt, mLUTLo[1]);
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
    Vec4d ringMod(RingMod());
    Vec4d mixed = mul_add(tb1 - tb0, 1 - tableOffset, tb0);
    mixed = mul_add(ringMod - 1., mRM * mRingModAmt * mixed, mixed);

    IOscillator<T>::mPhase += phaseIncr * (double)VECTOR_SIZE;
    IOscillator<T>::mPhase -= floor(IOscillator<T>::mPhase);

    T oversampled[VECTOR_SIZE];
    mixed.store(oversampled);

//    mLastOutput = pOutput[s] = mAAFilter.ProcessAndDownsample_Vector(mixed, mProcessOS);
    for (auto s = 0; s < VECTOR_SIZE / mProcessOS; ++s)
    {
      mLastOutput = pOutput[s] = mAAFilter.ProcessAndDownsample(oversampled + (s * mProcessOS));
    }

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
  inline double PhaseMod()
  {
    return mPM * mPhaseModAmt * mPhaseModulator.Process() * mCyclesPerLevelRecip;
  }

  // TODO: use enable_if to choose between vector lengths
  inline Vec4d PhaseModVec()
  {
    return Vec4d(0.);
  }

  inline double RingMod()
  {
    return mRingModulator.Process();
  }

  inline double* GetWtPosition()
  {
    return &mWtPosition;
  }

  inline double* GetWtBend()
  {
    return &mWtBend;
  }

  inline void SetWtPosition(double wtPos)
  {
    mWtPosition = wtPos;
  }

  inline void SetWtBend(double wtBend)
  {
    mWtBend = wtBend;
  }

  inline double SampleWavetablePosition(double phase)
  {
    return mWtPosition;
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

  static void NotifyLoaded(int oscIdx)
  {
    mWtReady[oscIdx] = true;
    mCV.notify_all();
  }

  T mLastOutput = 0.;

private:
#if OVERSAMPLING == 2
  inline static int mProcessOS{ 2 }; // Sample processing oversampling level (number of samples processed per sample output)
#else
  inline static int mProcessOS{ 1 };
#endif
  ChebyshevBL<T> mAAFilter;

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
  double mWtPosition{ 0 };
  double mWtSpacing{ 1. };
  double mWtBend{ 0 };
  int mPrevFreq;
  int mWtIdx{ 0 };

  // Thread-related members for wavetable updates
  static inline std::mutex mWtMutex; // Mutex used when swapping out the current wavetable in each oscillator object
  static inline std::condition_variable mCV;
  static inline bool TableLoaded{ false };
  static inline std::mutex mMasterMutex; // Static mutex used when loading a new wavetable from a file
  static inline bool mWtReady[2]{ false, false };

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
  static inline Wavetable<T>* LoadedTables[2]{ nullptr, nullptr };
  iplug::FastSinOscillator<T> mPhaseModulator;
  iplug::FastSinOscillator<T> mRingModulator{ 0.5 }; // Offset start phase by half a cycle
};