#pragma once

template<typename T>
class DelayLine
{
public:
  DelayLine(int length) : mLength(length)
  {
    assert (mLength > 1);
    mBuffer = new T[mLength]{};
  }

  ~DelayLine()
  {
    delete[] mBuffer;
  }

  void reset()
  {
    delete[] mBuffer;
    mBuffer = new T[mLength]{};
  }

  void SetDelay(int d)
  {
    mLength = d;
  }

  inline void push(T s)
  {
    mRead = mWrite;
    mBuffer[mWrite++] = s;
    if (mWrite == mLength)
      mWrite = 0;
  }

  inline T at(int offset = 0)
  {
    int readPoint{ mRead - offset };
    if (readPoint < 0)
      readPoint += mLength;
    return mBuffer[readPoint];
  }

  const T* GetPointer()
  {
    return mBuffer;
  }

  inline T operator[](int idx)
  {
    return at(idx);
  }

private:
  int mLength;
  T* mBuffer;
  int mRead{ 0 };
  int mWrite{ 0 };
};
