#pragma once

#define POLYPHASE_FILTER

#include "VectorFunctions.h"
#ifdef POLYPHASE_FILTER
#include "FastOversampler.h"
#else
#include "OversamplingFilter.h"
#endif
#include "Wavetable.h"

#include "Oscillator.h"

#include <vectorclass.h>

#ifdef VECTOR
#define OVERSAMPLING 2
#define VECTOR_SIZE 4
#define OUTPUT_SIZE VECTOR_SIZE / 2
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
  VectorOscillator(T startPhase=0.) : FastSinOscillator<T>(startPhase)
  {
  }
  
  inline Vec4d __vectorcall ProcessBlock4_Vector()
  {
    double phase = IOscillator<T>::mPhase + (double)UNITBIT32;
    const double phaseIncr = IOscillator<T>::mPhaseIncr * FastSinOscillator<T>::tableSize;
    Vec4d vPhase = phase + mPhaseOffset + phaseIncr * Vec4d(0., 1., 2., 3.);
    // Read 4 doubles as 8 integers (signed ints are used, but the AND operations later make this irrelevant)
    Vec8i viPhase = reinterpret_i(vPhase);
    // Upper 32 bits of 3*2^19 in odd indices, 0xFFFF in even ones: i.e. 0xFFFF, 0x18, 0xFFFF, ...
    Vec8i normhipart = blend8<8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V>(reinterpret_i(Vec4d((double)UNITBIT32)), Vec8i(0xFFFF));
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
    mPrevPhase = mPhase;
    tabfudge tf;
    phase += phaseIncr * 4.;
    tf.d = UNITBIT32 * tableSize;
    const int normhipart2 = tf.i[HIOFFSET];
    tf.d = phase + (UNITBIT32 * tableSize - UNITBIT32); // Remove the offset we introduced at the start of UNITBIT32.
    tf.i[HIOFFSET] = normhipart2;
    IOscillator<T>::mPhase = tf.d - UNITBIT32 * tableSize;
    return output;
  }

  const bool SyncSignal()
  {
    return mPrevPhase > mPhase;
  }

  void SetPhaseOffset(double offset)
  {
    mPhaseOffset = offset * tableSize;
  }

  inline Vec4d __vectorcall Process_Vector()
  {
    return ProcessBlock4_Vector();
  }

private:
  double mPhaseOffset{ 0. };
  double mPrevPhase{ 1. };
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

  static constexpr T MaxNote = 12543.8539514; // Highest Midi note frequency (440 * 2^((127-69)/2)), for setting the maximum formant shift

public:
  WavetableOscillator(const int id, const char* tableName, const double startPhase = 0., const double startFreq = 1.) :
    mID(id), IOscillator<T>(startPhase), mPrevFreq(static_cast<int>(startFreq))
  {
    WtFile table(tableName);
    LoadNewTable(table);
    SetWavetable(mLoadedTable);
    mWtReady = true;
  }

  WavetableOscillator(const int id, const WtFile& table, double startPhase = 0., double startFreq = 1.)
    : mID(id), IOscillator<T>(startPhase, startFreq), mPrevFreq(static_cast<int>(startFreq))
  {
    LoadNewTable(table);
    SetWavetable(mLoadedTable);
    mWtReady = true;
  }

  void SetSampleRate(double sampleRate)
  {
    mSampleRate = sampleRate * mProcessOS;
//    mNyquist = sampleRate / 2.;
    mPhaseModulator.SetSampleRate(mSampleRate);
    mRingModulator.SetSampleRate(mSampleRate);
  }

  /* Load a new wavetable as a static variable */
  void LoadNewTable(WtFile& wt)
  {
    if (wt.Success())
    {
      std::unique_lock<std::mutex> lock(mMasterMutex);
      mWtReady = false;
      delete mLoadedTable;
      mLoadedTable = new Wavetable<T>(wt);
    }
  }

  void SetWavetable(Wavetable<T>* tab)
  {
    std::unique_lock<std::mutex> lock(mWtMutex);
    mWtReady = false;
    if (tab != nullptr) // TODO: Check for nan's in the wavetable. They appear to break everything (including after new tables are loaded)
      mWT = tab;
    mCyclesPerLevelRecip = 1. / mWT->mCyclesPerLevel; // TODO: Try-catch block here in case wavetable cannot be found
    mPhaseIncrFactor = mCyclesPerLevelRecip / mProcessOS;
    mCycle = 0; // Must reset the cycle to avoid reading garbage memory when switching to wavetables with fewer cycles than the previous one
  }

  // Chooses the proper mipmap for the current note frequency (Hz) and sample rate
  inline void SetMipmapLevel()
  {
    const double samplesPerCycle{ 1. / IOscillator<T>::mPhaseIncr };

    // Select by Index
    const double tableFact{ std::log2(mWT->GetMaxSize() / (samplesPerCycle * mTableOS)) }; // Factors of two by which the largest mipmap (at index 0) is larger than the required mipmap
    mTableInterp = tableFact - std::floor(tableFact);
    mMipmapIdx = std::max(static_cast<int>(std::floor(tableFact)), 0); // Nearest integer lower than the value calculated above
    SetMipmapLevel_ByIndex(mMipmapIdx);
  }

  inline void SetMipmapLevel_ByIndex(const int idx)
  {
    std::unique_lock<std::mutex> lock(mWtMutex);
    mCV.wait(lock, [this] { return mWtReady; });

    const int prevTableSize{ mTableSize };

    mWtPositionNorm = 1 - std::modf((1. - mWtPositionAbs) * (mWT->mNumTables - 1), &mWtOffset);
    int tableOffset = std::min(static_cast<int>(mWtOffset), mWT->mNumTables - 2);
    mLUTLo[0] = mWT->GetMipmapLevel_ByIndex(tableOffset, idx, mTableSize);
    mLUTLo[1] = mWT->GetMipmapLevel_ByIndex(static_cast<size_t>(tableOffset) + 1, idx, mTableSize);
    mLUTHi[0] = mWT->GetMipmapLevel_ByIndex(tableOffset, idx + 1, mNextTableSize);
    mLUTHi[1] = mWT->GetMipmapLevel_ByIndex(static_cast<size_t>(tableOffset) + 1, idx + 1, mNextTableSize);
    mTableSize /= mWT->mCyclesPerLevel;
    mTableSizeM1 = mTableSize - 1;
    mNextTableSize /= mWT->mCyclesPerLevel;
    mNextTableSizeM1 = mNextTableSize - 1;

    // Adjust phase to new range
    double phaseAdj{ static_cast<double>(mTableSize) / prevTableSize };
    mPhase = mPhase * phaseAdj;
    mPrevPhase = mPrevPhase * phaseAdj;
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
    mMaxFormant = 0.25 * (MaxNote / freqCPS);

    // Set Formant
    const double freq_adj{ mFormantOn ? freqCPS * mFormant : freqCPS };
    mFormantModulator.SetFreqCPS(freqCPS);
    IOscillator<T>::SetFreqCPS(freq_adj);

    if (std::abs(mPrevFreq - freq_adj) > 0.1)
    {
      SetMipmapLevel();
    }
    else
    {
      SetMipmapLevel_ByIndex(mMipmapIdx);
    }
    mPrevFreq = freq_adj;
  }
#ifdef VECTOR
  inline Vec4d __vectorcall ProcessMultiple(double freqCPS)
  {
    AdjustWavetable(freqCPS);
#ifdef POLYPHASE_FILTER
#if OVERSAMPLING == 4
    T temp_in[16];
    T temp_out[4];
    const Vec4d output1 = ProcessOversamplingVec4();
    const Vec4d output2 = ProcessOversamplingVec4();
    const Vec4d output3 = ProcessOversamplingVec4();
    const Vec4d output4 = ProcessOversamplingVec4();
    output1.store(temp_in); output2.store(temp_in + 4);
    output3.store(temp_in + 8); output4.store(temp_in + 12);
    mAAFilter.mDownsampler4x.process_block(temp_out, temp_in, 8);
    mAAFilter.mDownsampler2x.process_block(temp_out, temp_in, 4);
#elif OVERSAMPLING == 2
    // Generate 8 samples, then resample 4
    T temp_in[8];
    T temp_out[4];
    const Vec4d output1 = ProcessOversamplingVec4();
    const Vec4d output2 = ProcessOversamplingVec4();
    output1.store(temp_in);
    output2.store(temp_in + 4);
    mAAFilter.mDownsampler2x.process_block(temp_out, temp_in, 4);
#endif
    Vec4d output{ Vec4d().load(temp_out) };
#else // Chebyshev Anti-Aliasing Filter:
    const Vec4d output1 = ProcessOversamplingVec4();
    const Vec4d output2 = ProcessOversamplingVec4();
    Vec4d output = mAAFilter.ProcessAndDownsample(output1, output2);
#endif // !Chebyshev filter
    return output;
  }
#else // Non-vectorized Processing:
  inline std::array<T, OUTPUT_SIZE> ProcessMultiple(double freqCPS)
  {
    AdjustWavetable(freqCPS);
    std::array<T, OUTPUT_SIZE> output{ 0. };
    ProcessOversampling(output, mProcessOS);
    return output;
  }
#endif // !non-vectorized processing

  inline T Process(double freqHz)
  {
    return 0.;
  }

  inline void ProcessOversampling(std::array<T, OUTPUT_SIZE>& pOutput, int nFrames)
  {
    std::unique_lock<std::mutex> lock(mWtMutex);

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
      const int offset{ static_cast<int>(phase) + mCycle * mTableSize };
      // Higher-frequency wavetable offset
      frac2 = modf(phaseShift * mNextTableSize, &phase);
      const int halfOffset{ static_cast<int>(phase) + mCycle * mNextTableSize };
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
      mCycle = (mCycle + (mPhase < mPrevPhase)) % mWT->mCyclesPerLevel;
      mPrevPhase = mPhase;
    }
#if OVERSAMPLING > 1
    for (auto s = 0; s < nFrames / mProcessOS; ++s)
    {
      mLastOutput = pOutput[s] = mAAFilter.ProcessAndDownsample(oversampled + (s * mProcessOS));
    }
    delete[] oversampled;
#endif

  }

  inline Vec4d __vectorcall ProcessOversamplingVec4()
  {
    double tableOffset{ mWtPositionAbs * (mWT->mNumTables - 1) };
    tableOffset -= std::max(floor(tableOffset - 0.0001), 0.);

    double phaseNormShifted = SamplePhaseShift(mPhase / mTableSize); // for phase shift
    double phaseUnshifted = mPhase + (double)UNITBIT32;
    double phase = phaseNormShifted * mTableSize + (double)UNITBIT32;
    const double phaseIncr = mPhaseIncr * mTableSize;
    const double phaseIncr2 = mPhaseIncr * mNextTableSize;
    Vec4d vPhase = phase + phaseIncr * mIncrVec;

    // Read 4 doubles as 8 integers (signed ints are used, but the AND operations later make this irrelevant)
    Vec8i viPhase = reinterpret_i(vPhase);
    // Get an int64 vector of 0's and 1's representing whether any phase values belong to the next cycle.
#if INSTRSET < 8
    Vec4q cycleIncr = select(reinterpret_i(permute8<HIOFFSET_V, -1, HIOFFSET_V + 2, -1, HIOFFSET_V + 4, -1, HIOFFSET_V + 6, -1>(viPhase & (mTableSize * 2 - 1))) > Vec4q(mTableSize), Vec4q(1), Vec4q(0));
#else
    Vec4q cycleIncr = select(reinterpret_i(permute8<HIOFFSET_V, -1, HIOFFSET_V + 2, -1, HIOFFSET_V + 4, -1, HIOFFSET_V + 6, -1>(viPhase & (mTableSize * 2 - 1))) > mTableSize, Vec4q(1), Vec4q(0));
#endif
    // Upper 32 bits of 3*2^19 in upper indices, 0xFFFF (32 bits of ones) in lower: i.e. 0xFFFF, 0x18, 0xFFFF, ...
    Vec8i normhipart = blend8<8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V>(reinterpret_i(Vec4d((double)UNITBIT32)), Vec8i(0xFFFF));
    // Mask the 8-item vector of 32-bit ints with one less than the table size, pad the upper bits (lower indices) with zeros,
    // and reinterpret as a 4-item vector of 64-bit ints
    Vec8i offsets32 = viPhase & mTableSizeM1;
    Vec4q offsets = reinterpret_i(permute8<HIOFFSET_V, -1, HIOFFSET_V + 2, -1, HIOFFSET_V + 4, -1, HIOFFSET_V + 6, -1>(offsets32));
    Vec4q phaseMod = to_int64_in_range(PhaseMod(), mTableSize);
    offsets = (offsets + phaseMod) & mTableSizeM1;
    Vec4q offsets1 = (offsets + 1) & mTableSizeM1;
    // Force the double to wrap. (AND the upper bits with the upper 32 bits of UNITBIT32)
    viPhase &= normhipart;
    Vec4d phaseFrac = reinterpret_d(viPhase) - (double)UNITBIT32; // get the fractional part

    // Larger/lower table
    constexpr int maxTableSize = 16384 * 12;
    Vec4q cyclePart = (mCycle + cycleIncr) * mTableSize;
    cyclePart -= select(cyclePart >= mWT->mCyclesPerLevel, Vec4q(mWT->mCyclesPerLevel), Vec4q(0)); // if mCycle is at a maximum and cycleIncr has nonzero components
    const Vec4d tb0s0lo = lookup<maxTableSize>(offsets + cyclePart, mLUTLo[0]);
    const Vec4d tb0s1lo = lookup<maxTableSize>(offsets1 + cyclePart, mLUTLo[0]);
    const Vec4d tb1s1lo = lookup<maxTableSize>(offsets1 + cyclePart, mLUTLo[1]);
    const Vec4d tb1s0lo = lookup<maxTableSize>(offsets + cyclePart, mLUTLo[1]);
    const Vec4d tb0lo = mul_add((tb0s1lo - tb0s0lo), phaseFrac, tb0s0lo);
    const Vec4d tb1lo = mul_add((tb1s1lo - tb1s0lo), phaseFrac, tb1s0lo);

    // Recalculate offsets for the smaller table
    vPhase = phaseNormShifted * mNextTableSize + (double)UNITBIT32 + phaseIncr2 * mIncrVec;

    // Read 4 doubles as 8 integers (signed ints are used, but the AND operations later make this irrelevant)
    viPhase = reinterpret_i(vPhase);
    offsets32 = viPhase & mNextTableSizeM1;
    offsets = reinterpret_i(permute8<HIOFFSET_V, -1, HIOFFSET_V + 2, -1, HIOFFSET_V + 4, -1, HIOFFSET_V + 6, -1>(offsets32));
    // Add the phase mod, this time just bit-shifting the previous value, and AND-out bits of or above the MSB of the table size
    offsets = (offsets + (phaseMod >> (mNextTableSize != mTableSize))) & mNextTableSizeM1;
//    offsets = ((offsets >> (mNextTableSize != mTableSize)) + to_int64_in_range(PhaseMod(), mNextTableSize)) & mNextTableSizeM1;
    offsets1 = (offsets + 1) & mNextTableSizeM1;
    // Force the double to wrap. (AND the upper bits with the upper 32 bits of UNITBIT32)
    viPhase &= normhipart;
    phaseFrac = reinterpret_d(viPhase) - (double)UNITBIT32; // get the fractional part
    
    // Smaller/higher table
    Vec4q cyclePartHi = (mCycle + cycleIncr) * mNextTableSize;
    cyclePartHi -= select(cyclePartHi >= mWT->mCyclesPerLevel, Vec4q(mWT->mCyclesPerLevel), Vec4q(0)); // if mCycle is at a maximum and cycleIncr has nonzero components
    const Vec4d tb0s0hi = lookup<maxTableSize>(offsets + cyclePartHi, mLUTHi[0]);
    const Vec4d tb0s1hi = lookup<maxTableSize>(offsets1 + cyclePartHi, mLUTHi[0]);
    const Vec4d tb1s1hi = lookup<maxTableSize>(offsets1 + cyclePartHi, mLUTHi[1]);
    const Vec4d tb1s0hi = lookup<maxTableSize>(offsets + cyclePartHi, mLUTHi[1]);
    // Interpolate samples
    const Vec4d tb0hi = mul_add(tb0s1hi - tb0s0hi, phaseFrac, tb0s0hi);
    const Vec4d tb1hi = mul_add(tb1s1hi - tb1s0hi, phaseFrac, tb1s0hi);

    // Interpolate mipmap levels
    Vec4d tb0 = mul_add(tb0hi - tb0lo, mTableInterp, tb0lo);
    Vec4d tb1 = mul_add(tb1hi - tb1lo, mTableInterp, tb1lo);
    
    // Restore mPhase
    tabfudge tf;
    mPrevPhase = mPhase;
    phase = phaseUnshifted;
    phase += phaseIncr * (double)VECTOR_SIZE;
    tf.d = UNITBIT32 * mTableSize;
    const int normhipart2 = tf.i[HIOFFSET];
    tf.d = phase + (UNITBIT32 * mTableSize - UNITBIT32); // Remove the offset we introduced at the start of UNITBIT32.
    tf.i[HIOFFSET] = normhipart2; // Reset the upper 32 bits
    if (!mFormantOn || !mFormantModulator.SyncSignal())
      IOscillator<T>::mPhase = (tf.d - UNITBIT32 * mTableSize);
    else
      Reset();
    mCycle = (mCycle + (mPhase < mPrevPhase)) % mWT->mCyclesPerLevel;

    // Mix wavtables and add ring and formant modulation
    Vec4d mixed = mul_add(tb1 - tb0, 1 - tableOffset, tb0);
    if (mFormantOn)
    {
      Vec4d formantMod = mFormantModulator.Process_Vector();
      mixed *= formantMod;
    }
    mixed = mul_add(mixed * (RingMod() - 1.), mRM * mRingModAmt, mixed);

    return mixed;
  }

  /* Returns a new phase as a power function of the current phase in the wavetable (i.e. always between 0. and 1.,
  regardless of the number of fundamental cycles in the table). */
  inline double SamplePhaseShift(double phase)
  {
    return std::pow(phase, 1. + (mWtBend >= 0. ? mWtBend : mWtBend * 0.5));
  }

  inline double WrapPhase(double phase)
  {
    while (phase < 0.)
      phase++;
    while (phase >= 1.)
      phase--;
    return phase;
  }

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

  // Phase Skew
  inline void SetWtBend(const double wtBend)
  {
    mWtBend = wtBend;
  }

  // Formant: Accepts a value between zero and one
  inline void SetFormant(const double fmtNorm, const bool formantOn = false)
  {
    mFormant = std::min(1. / (1. - fmtNorm), mMaxFormant);
    mFormantRecip = 1. / mFormant;
    mFormantOn = formantOn;
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

  void Reset()
  {
    IOscillator<T>::Reset();
    mPhaseModulator.SetPhase(0.);
    mRingModulator.SetPhase(0.);
    mFormantModulator.SetPhase(0.);
    mCycle = 0;
  }

  T mLastOutput = 0.;

private:
#if OVERSAMPLING == 2
  static constexpr int mProcessOS{ 2 }; // Sample processing oversampling level (number of samples processed per sample output)
#else
  static constexpr int mProcessOS{ 1 };
#endif
//  double mNyquist{ 20000. };

  // TODO: Order elements mindful of cache access:

  // Oscillator ID
  int mID;

  // Lookup Table Parameters
  int mTableSize = WT_SIZE; // Default: 2^9
  int mTableSizeM1 = WT_SIZE - 1; // Default: 2^9 -1
  int mNextTableSize = WT_SIZE / 2;
  int mNextTableSizeM1 = WT_SIZE / 2 - 1;
  int mCycle{ 0 };
  int mTableOS{ 8 }; // Wavetable oversampling level (ratio of table size to maximum samples per cycle read from the table)
  unsigned int mMipmapIdx{ 0 }; // Mipmap level
  double mPrevPhase{ 0. };
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
  double mFormant{ 1. };
  double mFormantRecip{ 1. };
  double mMaxFormant{ 1. }; // Set according to current pitch
  double mPrevFreq;

  // Thread-related members for wavetable updates
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
  bool mFormantOn{ false };

  static inline constexpr double twoPi{ 6.28318530718 };

public:
  Wavetable<T>* mLoadedTable{ nullptr };
#ifndef VECTOR
  iplug::FastSinOscillator<T> mPhaseModulator;
  iplug::FastSinOscillator<T> mRingModulator{ 0.5 }; // Offset start phase by half a cycle
#else
  iplug::VectorOscillator<T> mPhaseModulator;
  iplug::VectorOscillator<T> mRingModulator;
  iplug::VectorOscillator<T> mFormantModulator;
#endif

#ifdef POLYPHASE_FILTER
  FastOversampler<T> mAAFilter;
#else
  ChebyshevBL<T> mAAFilter;
#endif

  static inline std::mutex mWtMutex; // Mutex used when swapping out the current wavetable in each oscillator object
};