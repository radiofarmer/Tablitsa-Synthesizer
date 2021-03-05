#pragma once

#include "VectorFunctions.h"

#include <cmath>
#include <assert.h>

template<typename T>
extern inline T SoftClip(T s, T gain = (T)1);

template<typename T1, typename T2>
extern inline T1 __vectorcall SoftClip(const T1& s, T2 gain = (T2)1);

template<typename T, int MaxLength>
class DelayLine
{
  static constexpr int MaxLengthM1 = MaxLength - 1;

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

  template<class Vd>
  inline void push(const Vd& s)
  {
    typedef typename deduce_vector_from<Vd>::int_vec Vi;
    mRead = mWrite;
    scatter((Vi(0, 1, 2, 3) + mWrite) & MaxLengthM1, MaxLength, s, mBuffer);
    mWrite += deduce_vector_from<Vd>::v_size;
  }

  inline typename deduce_vector_from<T>::double_vec4 __vectorcall v_at(const int offset)
  {
    static_assert((MaxLength & (MaxLength - 1)) == 0); // Delay Line length must be a power of 2 in order to use vector functions
    typedef typename deduce_vector_from<typename deduce_vector_from<T>::double_vec4>::int_vec Vi;
    return lookup<MaxLength>((Vi(0, 1, 2, 3) + offset) & MaxLengthM1, mBuffer );
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

