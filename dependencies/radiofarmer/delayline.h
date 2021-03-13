#pragma once

#include "radiofarmer_config.h"

BEGIN_DSP_NAMESPACE

template<int MaxLength>
class DelayLine
{
  static constexpr int MaxLengthM1 = MaxLength - 1;

public:
  DelayLine() {}

  /* Clear the delay line and reset the read/write positions to the start. */
  void reset()
  {
    std::fill_n(mBuffer, MaxLength, (sample_t)0);
  }

  void DelayLine::SetDelay(const int d)
  {
    mLength = std::min(d, MaxLength);
  }

  /* Get a pointer to the start of the delay line */
  const sample_t* GetPointer() {return mBuffer;}

  const sample_t* GetReadPtr() { return &(mBuffer[mRead]); }

  /* Add a new sample to the first (most recent) position on the delay line, and adjust the read/write positions accordingly. */
  inline void push(const sample_t s)
  {
    mRead = mWrite;
    mBuffer[mWrite++] = s;
    if (mWrite >= mLength)
      mWrite = 0;
  }

  /* The the value in the delay line `offset` samples ago. (same as `operator[]`) */
  inline sample_t at(const int offset = 0)
  {
    int readPoint{ mRead - offset };
    if (readPoint < 0)
      readPoint += mLength;
    return mBuffer[readPoint];
  }

  /* The the value in the delay line `offset` samples ago. (same as `at`) */
  inline sample_t operator[](const int idx)
  {
    return at(idx);
  }

  const inline sample_t last()
  {
    return at(mLength - 1);
  }

  inline void push(const sample_v& s)
  {
    mRead = mWrite;
    scatter((sample_vi(0, 1, 2, 3) + mWrite) & MaxLengthM1, MaxLength, s, mBuffer);
    mWrite += 4;
  }

  inline typename sample_v __vectorcall v_at(const int offset)
  {
    static_assert((MaxLength & (MaxLength - 1)) == 0); // Delay Line length must be a power of 2 in order to use vector functions
    return lookup<MaxLength>((sample_vi(0, 1, 2, 3) + offset) & MaxLengthM1, mBuffer);
  }

private:
  sample_t mBuffer[MaxLength]{ 0. };
  int mLength{ MaxLength };
  int mRead{ 0 };
  int mWrite{ 0 };
};

END_DSP_NAMESPACE