#pragma once

#include <cmath>

template<typename T>
T SoftClip(T s, T gain = (T)1);

template<typename T, int MaxLength>
class DelayLine
{
public:
  DelayLine() {}

  /* Clear the delay line and reset the read/write positions to the start. */
  void DelayLine::reset()
  {
    std::fill_n(mBuffer, MaxLength, (T)0);
  }

  void DelayLine::SetDelay(const int d)
  {
    mLength = std::min(d, MaxDelay);
  }

  /* Get a pointer to the start of the delay line */
  const T* GetPointer() {return mBuffer;}

  const T* GetReadPtr() { return &(mBuffer[mRead]); }

  /* Add a new sample to the first (most recent) position on the delay line, and adjust the read/write positions accordingly. */
  inline void push(const T s)
  {
    mRead = mWrite;
    mBuffer[mWrite++] = s;
    if (mWrite == mLength)
      mWrite = 0;
  }

  /* The the value in the delay line `offset` samples ago. (same as `operator[]`) */
  inline T at(const int offset = 0)
  {
    int readPoint{ mRead - offset };
    if (readPoint < 0)
      readPoint += mLength;
    return mBuffer[readPoint];
  }
  /* The the value in the delay line `offset` samples ago. (same as `at`) */
  inline T operator[](const int idx)
  {
    return at(idx);
  }

private:
  T mBuffer[MaxLength]{ 0. };
  int mLength{ MaxLength };
  int mRead{ 0 };
  int mWrite{ 0 };
};

