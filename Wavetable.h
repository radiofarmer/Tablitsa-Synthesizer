#pragma once

#include "IPlugPlatform.h"
#include <ShlObj.h>
#include <Shlwapi.h>
#include <tchar.h>
#include <string>


//#define FFT
#define FFT_MAX_SIZE 32768

#if _DEBUG && !VST3_API
#define WT_DIR "..\\resources\\data\\wavetables\\"
#else
#define USE_APPDATA_PATH
#define WT_DIR "\\Tablitsa\\wavetables\\"
#if VST3_API
#include <locale>
#include <codecvt>
#endif
#endif

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


void fft_lp_filter(WDL_FFT_REAL* samples, const int length, int max_bin, int decimation=1)
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
        bin->im = 0.;
      }
    }
  }
  // Inverse FFT
  WDL_real_fft(samples, length, 1);
  for (int i = 0; i < length / decimation; ++i)
  {
    samples[i] = samples[i * decimation] * (WDL_FFT_REAL)0.5;
  }

  //return (double*)samples;
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

// TODO: To make this class portable, add the wavetable directory as a template variable
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
#ifndef USE_APPDATA_PATH
  WtFile(std::string fname) : mPath(WT_DIR + fname + ".wt")
#else
  WtFile(std::string fname)
#endif
  {
    std::fstream f;
    int headerSize{ sizeof(mHeader) }; // Size of WT header
    constexpr int byteInc = 8;
    
#ifdef USE_APPDATA_PATH
    TCHAR szPath[MAX_PATH];
    if (SUCCEEDED(SHGetFolderPath(NULL, CSIDL_APPDATA, NULL, 0, szPath)))
    {
      PathAppend(szPath, _T(WT_DIR));
    }
#if VST3_API
    // In the standalone app builds, `TCHAR` is type `char`. In VST3 builds, it's `wchar_t`, requiring a different conversion method.
    std::wstring wpath(szPath);
    std::string path = std::wstring_convert<std::codecvt_utf8<wchar_t>>().to_bytes(wpath);
#else
    std::string path = szPath;
#endif
    mPath = path + fname + ".wt";
#endif
    FILE* wt = fopen(mPath.c_str(), "rb");
    if (wt != nullptr)
    {
      fread(&mHeader, headerSize, 1, wt);
      mNumSamples = mHeader.numSamples;
      
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
 // static inline const std::string UserDir{ std::getenv("USER") };
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
    mValues = new T[size + 1]{};
    memcpy(mValues, values, size * sizeof(values[0]));
    // Set last entry equal to first (to avoid having to wrap integers twice during interpolated lookup functions)
    mValues[size] = mValues[0];

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

  T* GetLevel(const int length, int& tableSize)
  {
    const int curOffset{ 0 };
    const int curLevel{ mMaxLevel };
    const size_t levelIdx{ 0 };
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

    // Then generate the rest of the mipmap
    int levelNyquist{ (mStartSize / mCyclesPerLevel) / (4 * mTableOS) };
    // Before proceeding with band-limiting and downsampling check for further waveforms with more than 2^15 samples
    int currLevel = 0;
    do {
      memcpy(mValues + mLevelIndices[currLevel], values, mLevelSizes[currLevel] * sizeof(T)); // Copy the first wavetable and any further wavetables with more than 2^15 samples
      levelNyquist /= 2;
      currLevel++;
    } while (std::pow(2., std::ceil(std::log2(mLevelSizes[currLevel]))) > FFT_MAX_SIZE);

    T* buf = new T[static_cast<int>(FFT_MAX_SIZE)]{ 0. }; // Create a temporary FFT working buffer for generating new levels
    memcpy(buf, values + mLevelIndices[currLevel], mLevelSizes[currLevel] * sizeof(T)); // Populate the temporary FFT working buffer with the full-harmonic wavetable

    /* For each specified mipmap level, band-limit the temporary buffer to the level Nyquist frequency, which
    is decreased by a factor of two between each level. If the last-processed level is larger than the minimum
    specified level size, decimate it by a factor of two before copying it to the permanent mValues buffer. */
    for (auto i{currLevel}; i < num_levels - 1; ++i)
    {
      // Make sure the current level is within the size range for the FFT function
      assert (std::log2(mLevelSizes[i]) <= 15 && "Mipmap level exceeds maximum FFT length");
      // Get a suitable power of two for the FFT
      int powOfTwo = static_cast<int>(std::pow(2, static_cast<int>(std::ceil(std::log2(mLevelSizes[i])))));
      
      // Band-limit the mipmap level
      fft_lp_filter(buf, FFT_MAX_SIZE, levelNyquist, FFT_MAX_SIZE / mLevelSizes[i]);
      // Copy the filtered mipmap level to the mValues buffer
      memcpy(mValues + mLevelIndices[i], buf, mLevelSizes[i] * sizeof(T));

      // Reset buffer
      delete[] buf;
      T* buf = new T[static_cast<int>(FFT_MAX_SIZE)]{ 0. };
      memcpy(buf, values + mLevelIndices[currLevel], mLevelSizes[currLevel] * sizeof(T));

      // Adjust Nyquist frequency
      levelNyquist /= 2;
    }

    delete[] buf;
  }

  AutoMipmap<T>(const AutoMipmap<T>& source) :
    mValues(source.mValues), mStartSize(source.mStartSize), mMinSize(source.mMinSize), mNumLevels(source.mNumLevels),
    mCyclesPerLevel(source.mCyclesPerLevel), mTableOS(source.mTableOS), mNumSamples(source.mNumSamples)
  {
    mLevelSizes = source.mLevelSizes;
    mLevelIndices = source.mLevelIndices;
  }

  const int CalculateLength()
  {
    int level{ 0 };
    int nextLevelSize{ mStartSize };
    int samples{ 0 };
    while (level < mNumLevels)
    {
      // Increment total length and index of the current level
      mLevelIndices.push_back(samples);
      mLevelSizes.push_back(nextLevelSize);

      // Calculate the next level size and index
      samples += nextLevelSize;
      if (nextLevelSize > mMinSize) // Note: `mMinSize` has already been multipled by `mCyclesPerLevel`
        nextLevelSize /= 2;
      level++;
     }
    return samples;
  }

  T* GetLevel(const int length, int& tableSize)
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

  T* GetLevelAt(const int idx, int& tableSize)
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

    // Load mipmaps for each waveform (table position)
    for (int i{ 0 }; i < mNumTables; ++i)
    {
#ifdef FFT
      AutoMipmap<T> mm{ wt.GetSampleOffset(i * mWtSize), mMipmapMaxSize * mCyclesPerLevel, 1024 * mCyclesPerLevel, 8, mCyclesPerLevel, wt.Oversampling()};
#else
      Mipmap<T> mm(wt.GetSampleOffset(i * mWtSize), mWtSize, mMipmapMaxSize, mMipmapMinSize, mCyclesPerLevel);
#endif
      mWaveforms.push_back(mm);
    }
  }

  inline T* GetMipmapLevel_ByIndex(const std::size_t tableIdx, const int idx, int& tableSize)
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