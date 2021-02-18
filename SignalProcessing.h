#pragma once

#include <cmath>

double SoftClip(double s, double gain = 3.);

class DelayLine
{
public:
  DelayLine(const int length);

  ~DelayLine()
  {
    delete[] mBuffer;
  }

  /* Clear the delay line and reset the read/write positions to the start. */
  void reset();

  void SetDelay(const int d);

  /* Get a pointer to the start of the delay line */
  const double* GetPointer() {return mBuffer;}

  const double* GetReadPtr() { return &(mBuffer[mRead]); }

  /* Add a new sample to the first (most recent) position on the delay line, and adjust the read/write positions accordingly. */
  inline void push(const double s)
  {
    mRead = mWrite;
    mBuffer[mWrite++] = s;
    if (mWrite == mLength)
      mWrite = 0;
  }

  /* The the value in the delay line `offset` samples ago. (same as `operator[]`) */
  inline double at(const int offset = 0)
  {
    int readPoint{ mRead - offset };
    if (readPoint < 0)
      readPoint += mLength;
    return mBuffer[readPoint];
  }
  /* The the value in the delay line `offset` samples ago. (same as `at`) */
  inline double operator[](const int idx)
  {
    return at(idx);
  }

private:
  int mLength;
  double* mBuffer;
  int mRead{ 0 };
  int mWrite{ 0 };
};

