#pragma once

#include "IPlugPlatform.h"

//#define FFT

#define WT_DIR "C:\\Users\\bramc\\VST-Development\\iPlug2\\Examples\\Tablitsa\\resources\\data\\wavetables\\"
#define WT_SIZE 1024
#define WAV16_MAX 32767
#define WT_MAX_DEFAULT 16384
#define WT_MIN_DEFAULT 32


#include <cstdint>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <mutex>
#include "vectormath_exp.h"
#include "fft.h"

#include "Filter.h"
#include "Oscillator.h"
#include "Modulators.h"

#ifndef _DEBUG
  #define OVERSAMPLING
  #define VECTOR_SIZE 4
  #define OUTPUT_SIZE 2
  #define VECTOR
#else
  #define VECTOR_SIZE 1
  #define OUTPUT_SIZE 1
#endif

#if MULTITHREAD_TEST
#include <mutex>
#endif

double* fft_lp_filter(WDL_FFT_REAL* samples, const int length, int max_bin, int decimation=1)
{
  WDL_fft_init();
  WDL_real_fft(samples, length, 0);
  const int* order = WDL_fft_permute_tab(length / 2);
  for (int i{ 0 }; i < length; ++i)
  {
    samples[i] /= (WDL_FFT_REAL)length; // Scale all bins
    if (i <= length / 2)
    {
      WDL_FFT_COMPLEX* bin = (WDL_FFT_COMPLEX*)samples + order[i]; // Get pointer to current bin
      if (i > max_bin)
      {
        bin->re = 0.;
        //bin->im = 0.;
      }
    }
  }
  // Inverse FFT
  WDL_real_fft(samples, length, 1);
  for (int i = 0; i < length / decimation; ++i)
  {
    samples[i] = samples[i * decimation] * (WDL_FFT_REAL)0.5;
  }

  return (double*)samples;
}

enum endian
{
  little,
  big
};

struct WavHeader
{
  uint32_t chunkID;
  uint32_t chunkSize;
  uint32_t format;
  uint32_t subchunk1ID;
  uint32_t subchunk1Size;
  uint16_t audioFormat;
  uint16_t numChannels;
  uint32_t sampleRate;
  uint32_t byteRate;
  uint16_t blockAlign;
  uint16_t bitsPerSample;
  uint32_t subchunk2ID;
  uint32_t subchunk2Size;
};

class WtFile
{
  struct WtHeader
  {
    uint32_t numSamples;
    uint32_t maxSamples;
    uint32_t minSamples;
    uint32_t numWaveforms;
    uint32_t cyclesPerLevel;
    uint32_t oversampling;
    //uint32_t numLevels;
  };

public:
  WtFile(std::string fname) :
    mPath(WT_DIR + fname + ".wt")
  {
    std::fstream f;
    int headerSize{ sizeof(mHeader) }; // Size of WAV header
    constexpr int byteInc = 8;

    FILE* wt = fopen(mPath.c_str(), "rb");
    if (wt != nullptr)
    {
      fread(&mHeader, headerSize, 1, wt);
      mNumSamples = mHeader.numSamples;
      
      /*for (int i{ 0 }; i < mHeader.numLevels; ++i)
      {
        int harmonics;
        fread(&harmonics, 4, 1, wt);
        mLevels.push_back(harmonics);
      }*/
      
      // Temporary array for samples
      double_t* samplesRaw{ new double_t[static_cast<size_t>(mNumSamples)] };

      // Allocate necessary amount of space
      mData = new double[static_cast<size_t>(mNumSamples)];

      for (int b = 0; b < mNumSamples; ++b)
      {
        fread(&samplesRaw[b], byteInc, 1, wt);
        mData[b] = samplesRaw[b];
      }

      delete[] samplesRaw;
      fclose(wt);
      mSuccess = true;
    }
    else
      mSuccess = false;

    f.close();
  }

  const int NumCycles()
  {
    return mHeader.cyclesPerLevel;
  }

  const int MaxLevel()
  {
    return mHeader.maxSamples;
  }

  const int MinLevel()
  {
    return mHeader.minSamples;
  }

  const int Oversampling()
  {
    return mHeader.oversampling;
  }

  const int NumSamples()
  {
    return mNumSamples;
  }

  const int NumWaveforms()
  {
    return mHeader.numWaveforms;
  }

  double GetSample(int idx)
  {
    return mData[idx];
  }

  double* GetSampleOffset(int offset = 0)
  {
    return mData + offset;
  }

  bool Success()
  {
    return mSuccess;
  }

  ~WtFile()
  {
    if (Success())
      delete[] mData;
  }

private:
  std::string mPath{};
  WtHeader mHeader{};
  int mNumSamples{ 0 };
  double* mData; // Samples
  bool mSuccess{ false };
  std::vector<int> mLevels;
};

class WavFile
{
  inline static constexpr int read_fail{ -1 };

public:
  WavFile(std::string fname) :
    mPath(WT_DIR + fname + ".wav")
  {
    std::fstream f;
    int headerSize{ sizeof(mHeader) }; // Size of WAV header

    FILE* wavFile = fopen(mPath.c_str(), "rb");
    if (wavFile != nullptr)
    {
      fread(&mHeader, headerSize, 1, wavFile);
      const int byteInc{ mHeader.bitsPerSample / 8 };
      mNumSamples = mHeader.subchunk2Size / byteInc;

      // Temporary array for samples
      int16_t* samplesRaw{ new int16_t[static_cast<size_t>(mNumSamples)] };
      // Seek to start
      // const int offset{ 16384 + 8192 + 4096 + 2048 + 1024 };
      // fseek(wavFile, offset, SEEK_CUR);

      // Allocate necessary amount of space
      mData = new double[static_cast<size_t>(mNumSamples)];

      for (int b = 0; b < mNumSamples; ++b)
      {
        fread(&samplesRaw[b], byteInc, 1, wavFile);
        mData[b] = static_cast<double>(samplesRaw[b]) / WAV16_MAX;
      }
      delete[] samplesRaw;
      fclose(wavFile);
      mSuccess = 1;
    }
    else
      mSuccess = read_fail;

    f.close();
  }

  double getSample(int idx)
  {
    return mData[idx];
  }

  double* getSampleOffset(int offset=0)
  {
    return mData + offset;
  }

  int numSamples()
  {
    return mNumSamples;
  }

  const bool Success()
  {
    if (mSuccess != read_fail)
      return true;
    else
      return false;
  }

  ~WavFile()
  {
    if (Success())
      delete[] mData;
  }

private:
  std::string mPath{}; // Path of wav file
  WavHeader mHeader{}; // Standard WAV header data
  int mNumSamples{ 0 };
  int mEndianness[14]{big, little, big, big, little, little, little, little, little, little, little, big, little, little}; // Endianness of each header component
  double* mData; // Samples
  int mSuccess{ read_fail };
};

template <typename T>
class Mipmap
{
public:
  Mipmap(T* values, int size, int maxLevel, int minLevel, int cyclesPerLevel, int factor = 2) :
    mSize(size), mMaxLevel(maxLevel), mMinLevel(minLevel), mCyclesPerLevel(cyclesPerLevel), mFactor(factor)
  {
    // This ought to be the rate-limiting step in loading a new wavetable into the oscillator
    mValues = new T[size]{};
    memcpy(mValues, values, size * sizeof(values[0]));

    // Calculate Sizes of all levels
    int pos{ 0 };
    int nextLevelSize{ mMaxLevel * mCyclesPerLevel };
    while (pos < mSize)
    {
      // Increment total length and index of the current level
      mLevelIndices.push_back(pos);
      mLevelSizes.push_back(nextLevelSize);

      // Calculate the next level size and index
      pos += nextLevelSize;
      if (nextLevelSize > mMinLevel * mCyclesPerLevel)
        nextLevelSize /= 2;
    }
  }

  Mipmap<T>(const Mipmap<T>& source) : Mipmap(source.mValues, source.mSize, source.mMaxLevel, source.mMinLevel, source.mCyclesPerLevel, source.mFactor)
  {
    //mValues = source.mValues;
  }

  T* GetLevel(int length, int& tableSize)
  {
    int curOffset{ 0 };
    int curLevel{ mMaxLevel };
    size_t levelIdx{ 0 };
    // Select the mipmap level as if all wavetables had one cycle and no size floor
    while (curLevel > length)
    {
      curOffset = mLevelIndices[levelIdx];
      curLevel /= 2;
      levelIdx++;
    }
    tableSize = mLevelSizes[levelIdx];
    return mValues + curOffset;
  }

  T* GetLevelAt(int idx, int& tableSize)
  {
    idx = std::min(idx, (int)mLevelIndices.size() - 2);
    tableSize = mLevelSizes[idx];
    return mValues + mLevelIndices[idx];
  }

  const int Size()
  {
    return mSize;
  }

  ~Mipmap()
  {
    if (mValues != nullptr)
    {
      delete[] mValues;
      mValues = nullptr;
    }
  }

protected:
  const int mMaxLevel;
  const int mMinLevel;
  const int mSize;
  const int mCyclesPerLevel;
  const int mFactor;
  std::vector<int> mLevelSizes;
  std::vector<int> mLevelIndices;
  T* mValues{ nullptr };
};

template <typename T>
class AutoMipmap
{
public:
  AutoMipmap(T* values, int start_size, int min_size, int num_levels, int cyclesPerLevel=1, int oversampling=8) :
    mStartSize(start_size), mMinSize(min_size), mNumLevels(num_levels), mCyclesPerLevel(cyclesPerLevel), mTableOS(oversampling)
  {
    mNumSamples = CalculateLength();
    mValues = new T[mNumSamples];
    memcpy(mValues, values, mStartSize * sizeof(mValues[0]));

    // Generate the rest of the mipmap
    int levelNyquist{ start_size / (4 * mTableOS) };
    T* buf = new T[mLevelSizes[0]]; // Create a temporary buffer for generating new levels
    memcpy(buf, values, mLevelSizes[0] * sizeof(mValues[0])); // Populate the temporary buffer with the full-harmonic wavetable

    /*
    For each specified mipmap level, band-limit the temporary buffer to the level Nyquist frequency, which
    is decreased by a factor of two between each level. If the last-processed level is larger than the minimum
    specified level size, decimate it by a factor of two before copying it to the permanent mValues buffer.
    */
    for (size_t i{0}; i < num_levels - 1; ++i)
    {
      // Band-limit the mipmap level
      buf = fft_lp_filter(buf, mLevelSizes[i], levelNyquist, mLevelSizes[i+1] > mMinSize ? 2 : 1);
      // Copy the filtered mipmap level to the mValues buffer
      memcpy(&(mValues[mLevelIndices[i + 1]]), buf, mLevelSizes[i + 1] * sizeof(mValues[0]));
      levelNyquist /= 2;
    }

    delete[] buf;
  }

  AutoMipmap<T>(const AutoMipmap<T>& source) :
    mValues(source.mValues), mStartSize(source.mStartSize), mMinSize(source.mMinSize), mNumLevels(source.mNumLevels), mCyclesPerLevel(source.mCyclesPerLevel), mTableOS(source.mTableOS)
  {
    mLevelSizes = source.mLevelSizes;
    mLevelIndices = source.mLevelIndices;
  }

  const int CalculateLength()
  {
    int level{ 0 };
    int nextLevelSize{ mStartSize * mCyclesPerLevel };
    int samples{ 0 };
    while (level < mNumLevels)
    {
      // Increment total length and index of the current level
      mLevelIndices.push_back(samples);
      mLevelSizes.push_back(nextLevelSize);

      // Calculate the next level size and index
      samples += nextLevelSize;
      if (nextLevelSize > mMinSize * mCyclesPerLevel)
        nextLevelSize /= 2;
      level++;
    }
    return samples;
  }

  T* GetLevel(int length, int& tableSize)
  {
    int curOffset{ 0 };
    int curLevel{ mMaxLevel };
    size_t levelIdx{ 0 };
    // Select the mipmap level as if all wavetables had one cycle and no size floor
    while (curLevel > length)
    {
      curOffset = mLevelIndices[levelIdx];
      curLevel /= 2;
      levelIdx++;
    }
    tableSize = mLevelSizes[levelIdx];
    return mValues + curOffset;
  }

  T* GetLevelAt(int idx, int& tableSize)
  {
    idx = std::min(idx, (int)mLevelIndices.size() - 2);
    tableSize = mLevelSizes[idx];
    return mValues + mLevelIndices[idx];
  }

private:
  const int mNumLevels;
  const int mStartSize;
  const int mMinSize;
  const int mCyclesPerLevel;
  int mNumSamples;
  int mTableOS;
  std::vector<int> mLevelIndices;
  std::vector<int> mLevelSizes;
  std::vector<int> mLevelNyquist;
  T* mValues{ nullptr };
};

template <typename T>
class Wavetable
{
public:
  Wavetable(WavFile& wav, int nTimbres = 2, int maxSize = WT_MAX_DEFAULT, int minSize = WT_MIN_DEFAULT, int cyclesPerLevel = 1, int oversampling = 8) :
    mNumTables(nTimbres), mMipmapMaxSize(maxSize), mMipmapMinSize(minSize), mCyclesPerLevel(cyclesPerLevel)
  {
    mWtSize = wav.numSamples() / nTimbres;
    // Load mipmaps
    for (int i{ 0 }; i < mNumTables; ++i)
    {
      Mipmap<T> mm(wav.getSampleOffset(i * mWtSize), mWtSize, mMipmapMaxSize, mMipmapMinSize, mCyclesPerLevel, oversampling);
      mWaveforms.push_back(mm);
    }
  }

  /*
  Load from WT format file
  */
  Wavetable(WtFile& wt)
  {
    mMipmapMaxSize = wt.MaxLevel();
    mMipmapMinSize = wt.MinLevel();
    mCyclesPerLevel = wt.NumCycles();
    mNumTables = wt.NumWaveforms();
    mWtSize = wt.NumSamples() / mNumTables;

    // Load mipmaps
    for (int i{ 0 }; i < mNumTables; ++i)
    {
#ifdef FFT
      AutoMipmap<T> mm{ wt.GetSampleOffset(i * mWtSize), mMipmapMaxSize * mCyclesPerLevel, 1024, 8, mCyclesPerLevel, wt.Oversampling()};
#else
      Mipmap<T> mm(wt.GetSampleOffset(i * mWtSize), mWtSize, mMipmapMaxSize, mMipmapMinSize, mCyclesPerLevel);
#endif
      mWaveforms.push_back(mm);
    }
  }

  inline T* GetMipmapLevel_ByIndex(std::size_t tableIdx, int idx, int& tableSize)
  {
    return mWaveforms.at(tableIdx).GetLevelAt(idx, tableSize);
  }

  inline T* GetMipmapLevel_BySize(std::size_t tableIdx, int length, int& tableSize)
  {
    return mWaveforms.at(tableIdx).GetLevel(length, tableSize);
  }

  T GetSample(std::size_t timbre, std::size_t length, std::size_t sample)
  {
    return (this->GetMipmapLevel_BySize(timbre, length))[sample];
  }

  int GetMaxSize()
  {
    return mMipmapMaxSize;
  }

  int GetMinSize()
  {
    return mMipmapMinSize;
  }

  ~Wavetable()
  {
  }

private:
#ifdef FFT
  std::vector<AutoMipmap<T>> mWaveforms;
#else
  std::vector<Mipmap<T>> mWaveforms;
#endif

public:
  int mNumTables;
  int mMipmapMaxSize;
  int mMipmapMinSize;
  int mCyclesPerLevel;
  int mWtSize;
};

template <typename T>
class WavetableOscillator : public iplug::IOscillator<T>
{
  union tabfudge
  {
    double d;
    int i[2];
  } ALIGNED(8);

public:
  WavetableOscillator(double startPhase = 0., double startFreq = 1.)
    : IOscillator<T>(startPhase* mTableSizeM1, startFreq), mPrevFreq(static_cast<int>(startFreq))
  {
    SetWavetable("Hydrogen");
  }

  void SetSampleRate(double sampleRate) override
  {
    mSampleRate = sampleRate * mProcessOS;
  }

  /* Load a new wavetable as a static variable */
  static void LoadNewTable(WtFile& wt, int idx)
  {
    if (wt.Success())
    {
      std::unique_lock<std::mutex> lock(mWtMutex);
      TableLoaded = false;
      delete LoadedTables[idx];
      LoadedTables[idx] = new Wavetable<T>(wt);
      //TableLoaded = true;
    }
  }

  /*
  Select wavetable with a numerical index
  */
  void SetWavetable(int tableIdx)
  {
    assert(WavetableOscillator<T>::mWavetables.size() > tableIdx);
    SetWavetable((WavetableOscillator<T>::mWavetables.at(tableIdx)));
  }

  /*
  Set wavetable by filepath
  */
  void SetWavetable(std::string path)
  {
    //WavFile wav(path);
    //SetWavetable(wav);
    WtFile wt(path);
    SetWavetable(wt);
  }

  /*
  Load a wavetable from a standard RIFF WAV file
  */
  void SetWavetable(WavFile& wav)
  {
    if (wav.Success())
    {
      std::unique_lock<std::mutex> lock(mWtMutex);
      if (mWT != nullptr)
        delete mWT;
      mWT = new Wavetable<T>(wav, 2);
    }
  }

  /*
  Load a wavetable from the custom .wt format
  */
  void SetWavetable(WtFile& wt)
  {
    if (wt.Success())
    {
      std::unique_lock<std::mutex> lock(mWtMutex);
      Wavetable<T>* newTable = new Wavetable<T>(wt);
      mWtReady = false;
      if (mWT != nullptr)
        delete mWT;
      mWT = newTable; // This step takes about as long as generating the new table in the first place. Copy function for Wavetable probably needs to be overridden and optimized
      mWtReady = true;
    }
  }

  void SetWavetable(Wavetable<T>* tab)
  {
    std::unique_lock<std::mutex> lock(mWtMutex);
    mWtReady = false;
    mWT = tab;
    mPhaseIncrFactor = (1. / (mWT->mCyclesPerLevel * mProcessOS));
    mWtReady = true;
  }

  // Chooses the proper mipmap for the current note frequency (Hz) and sample rate
  inline void SetMipmapLevel()
  {
    const double samplesPerCycle{ 1. / IOscillator<T>::mPhaseIncr };

    // Select by Index
    double tableFact{ std::log2(mWT->GetMaxSize() / (samplesPerCycle * mTableOS)) };
    mTableInterp = tableFact - std::floor(tableFact);
    mWtIdx = std::max(static_cast<int>(std::floor(tableFact)), 0);
    SetMipmapLevel_ByIndex(mWtIdx);
  }

  inline void SetMipmapLevel_ByIndex(int idx)
  {
    std::unique_lock<std::mutex> lock(mWtMutex);
    mCV.wait(lock, [this] { return mWtReady && TableLoaded; });

    assert(mWT != nullptr);

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

  inline std::array<T, OUTPUT_SIZE> ProcessMultiple(double freqCPS) override
  {
    AdjustWavetable(freqCPS);

#ifdef VECTOR
    std::array<T, OUTPUT_SIZE> output{ 0. };
    ProcessOversamplingVec4(output);
#else
    std::array<T, VECTOR_SIZE> output{ 0. };
    ProcessOversampling(output, mProcessOS);
#endif
    return output;
  }

  inline void ProcessOversampling(std::array<T, VECTOR_SIZE>& pOutput, int nFrames)
  {

    double phase; // integer phase
    double frac, frac2; // fractional phases

    // Get the normalized offset of the current wavetable block
    double tableOffset{ mWtPosition * (mWT->mNumTables - 1) };
    tableOffset -= std::max(floor(tableOffset - 0.0001), 0.);
#ifdef OVERSAMPLING
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
//      const double frac = mTableSize * phaseNorm - floor(mTableSize * phaseNorm);
      // Read from wavetable
      const double sampleWtPosition{ tableOffset  };
      const double sampleWtPositionInv{ 1 - sampleWtPosition };
      const T f1 = addr[0][0] * sampleWtPosition + addr[1][0] * sampleWtPositionInv;
      const T f2 = addr[0][1] * sampleWtPosition + addr[1][1] * sampleWtPositionInv;
      const T f3 = addr[2][0] * sampleWtPosition + addr[3][0] * sampleWtPositionInv;
      const T f4 = addr[2][1] * sampleWtPosition + addr[3][1] * sampleWtPositionInv;
      // Send output
      T lowerTable = (f1 + frac * (f2 - f1));
#ifdef OVERSAMPLING
      oversampled[s] = lowerTable + mTableInterp * ((f3 + frac2 * (f4 - f3)) - lowerTable);
#else
      const T output = lowerTable + mTableInterp * ((f3 + frac2 * (f4 - f3)) - lowerTable);
      mLastOutput = pOutput[s] = output + mRM * mRingModAmt * (RingMod() * output - output);
#endif
      // Increment Phase
      IOscillator<T>::mPhase += mPhaseIncr * mPhaseIncrFactor;
      IOscillator<T>::mPhase -= std::floor(IOscillator<T>::mPhase);
    }
#ifdef OVERSAMPLING
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

    const double phaseIncr = IOscillator<T>::mPhaseIncr * mPhaseIncrFactor;
    double phase = SamplePhaseShift(IOscillator<T>::mPhase);
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
    phaseInt = truncatei(phaseDouble);
    phaseFrac = phaseDouble - to_double(phaseInt);
    phaseInt &= mNextTableSizeM1;
    phaseInt1 = (phaseInt + 1) & mNextTableSizeM1;

    // Read from wavetables (higher/smaller table)
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
    Vec4d mixed = mul_add(tb1 - tb0, 1 - tableOffset, tb0);

    IOscillator<T>::mPhase += phaseIncr * (double)VECTOR_SIZE;
    IOscillator<T>::mPhase -= floor(IOscillator<T>::mPhase);

    T oversampled[4];
    mixed.store(oversampled);

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
    return (cycle + std::pow(mPhaseInCycle, 1. + (mWtBend >= 0. ? mWtBend : mWtBend / 2.))) / mWT->mCyclesPerLevel;
  }

  inline double WrapPhase(double phase)
  {
    while (phase < 0)
      phase++;
    while (phase > 1)
      phase--;
    return phase;
  }

  /* Returns an adjusted phase increment based on the current (cycle-normalized) phase. */
  inline double PhaseMod()
  {
    return mPM * mPhaseModAmt * mPhaseModulator.Lookup(mPhaseInCycle * mPhaseModFreq);
  }

  inline double RingMod()
  {
    return mRingModulator.Lookup(mPhaseInCycle * mRingModFreq);
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
  inline void SetPhaseModulation(double amt, double freqOffset)
  {
    mPhaseModAmt = amt;
    mPhaseModFreq = freqOffset * twoPi;
  }

  inline void SetRingModulation(bool on)
  {
    if (on)
      mRM = 1.;
    else
      mRM = 0.;
  }
  inline void SetRingModulation(double amt, double freqOffset)
  {
    mRingModAmt = amt;
    mRingModFreq = freqOffset * twoPi;
  }

  void ReloadLUT()
  {
    mPrevFreq = -1;
  }

  static void NotifyLoaded(bool isLoaded = true)
  {
    TableLoaded = true;
    mCV.notify_all();
  }

  T mLastOutput = 0.;

private:
#ifdef OVERSAMPLING
  inline static int mProcessOS{ 2 }; // Sample processing oversampling level (number of samples processed per sample output)
#else
  inline static int mProcessOS{ 1 };
#endif
  ChebyshevBL<T> mAAFilter;

  // Lookup Table Parameters
  int mTableSize = WT_SIZE; // Default: 2^9
  int mTableSizeM1 = WT_SIZE - 1; // Default: 2^9 -1
  int mNextTableSize = WT_SIZE / 2;
  int mNextTableSizeM1 = WT_SIZE / 2 - 1;
  int mTableOS{ 8 }; // Wavetable oversampling level (ratio of table size to maximum samples per cycle read from the table)
  double mTableInterp{ 1 };
  double mPhaseIncrFactor{ 1. }; // Reciprocal of the product of the number of cycles per wavetable level and the processing oversampling level
  double mPhaseInCycle{ 0. }; // The fractional phase (between zero and one) within a single cycle, independent of the number of cycles per wavetable level.

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
  static inline std::mutex mWtMutex;
  static inline std::condition_variable mCV;
  static inline bool TableLoaded{ false };
  bool mWtReady{ false };

  // Vectors
  const Vec4d mIncrVec{ 0., 1., 2., 3. };

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
  static inline iplug::FastSinOscillator<T> mPhaseModulator;
  static inline iplug::FastSinOscillator<T> mRingModulator;
};