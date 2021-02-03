#pragma once

class DelayLine
{
public:
  DelayLine(int length) : mLength(length)
  {
    assert (mLength > 1);
    mBuffer = new double[mLength]{};
  }

  ~DelayLine()
  {
    delete[] mBuffer;
  }

  void reset()
  {
    delete[] mBuffer;
    mBuffer = new double[mLength]{};
  }

  void SetDelay(int d)
  {
    mLength = d;
  }

  inline void push(double s)
  {
    mRead = mWrite;
    mBuffer[mWrite++] = s;
    if (mWrite == mLength)
      mWrite = 0;
  }

  inline double at(int offset = 0)
  {
    int readPoint{ mRead - offset };
    if (readPoint < 0)
      readPoint += mLength;
    return mBuffer[readPoint];
  }

  const double* GetPointer()
  {
    return mBuffer;
  }

  inline double operator[](int idx)
  {
    return at(idx);
  }

private:
  int mLength;
  double* mBuffer;
  int mRead{ 0 };
  int mWrite{ 0 };
};
