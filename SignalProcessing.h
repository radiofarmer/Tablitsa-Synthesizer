#pragma once

//#include <assert.h>

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

  inline void push(const double s)
  {
    mRead = mWrite;
    mBuffer[mWrite++] = s;
    if (mWrite == mLength)
      mWrite = 0;
  }

  inline double at(const int offset = 0);

  const double* GetPointer()
  {
    return mBuffer;
  }

  double operator[](const int idx);

private:
  int mLength;
  double* mBuffer;
  int mRead{ 0 };
  int mWrite{ 0 };
};

