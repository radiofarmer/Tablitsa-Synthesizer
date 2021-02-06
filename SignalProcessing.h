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

template<typename T>
class Effect
{
public:
  Effect(double sampleRate) : mSampleRate(sampleRate) {}

  virtual T Process(T s) = 0;

  void SetSampleRate(double sampleRate)
  {
    mSampleRate = sampleRate;
  }

protected:
  double mSampleRate;
};

template<typename T>
class DelayEffect final : public Effect<T>
{
public:
  DelayEffect(double sampleRate, double maxDelayMS = 5000.) :
    Effect(sampleRate),
    mMaxDelayMS(maxDelayMS),
    mMaxDelay(static_cast<int>(mMaxDelayMS / 1000. * mSampleRate)),
    mDelayL(mMaxDelay),
    mDelayR(mMaxDelay),
    mDelayLTime(mMaxDelay / 2),
    mDelayRTime(mMaxDelay / 2)
  {

  }

  T Process(T s) override
  {
    mDelayL.push(s);
    mDelayR.push(s);
    return mDelayL[mDelayLTime] * mDelayLGain + mDelayR[mDelayRTime] * mDelayRGain;
  }

private:
  const double mMaxDelayMS;
  int mMaxDelay;
  DelayLine mDelayL;
  DelayLine mDelayR;

  T mDelayLGain{ 0.5 };
  T mDelayRGain{ 0.5 };
  int mDelayLTime;
  int mDelayRTime;
};