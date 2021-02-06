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
    mDelayL(mMaxDelay + 1),
    mDelayR(mMaxDelay + 1),
    mDelayLTime(mMaxDelay / 2),
    mDelayRTime(mMaxDelay / 2)
  {

  }

  void SetDelay(T timeMS, int channel)
  {
    if (channel == 0)
    {
      mDelayLTimeMS = timeMS;
      mDelayLTime = static_cast<int>(timeMS / 1000. * mSampleRate);
    }
    else if (channel == 1)
    {
      mDelayRTimeMS = timeMS;
      mDelayRTime = static_cast<int>(timeMS / 1000. * mSampleRate);
    }
  }

  void SetFeedback(T fb)
  {
    mFeedback = fb;
  }

  void SetGain(T gain)
  {
    mDelayLGain = gain;
    mDelayRGain = gain;
  }

  T Process(T s) override
  {
    T left_out = mDelayL[mDelayLTime];
    T right_out = mDelayR[mDelayRTime];
    mDelayL.push(s + left_out * mFeedback);
    mDelayR.push(s + right_out * mFeedback);
    return left_out * mDelayLGain + right_out * mDelayRGain;
  }

  T* ProcessStereo(T sl, T sr)
  {
    T left_out = mDelayL[mDelayLTime];
    T right_out = mDelayR[mDelayRTime];
    mDelayL.push(sl + left_out * mFeedback);
    mDelayR.push(sr + right_out * mFeedback);
    T output[2]{ left_out * mDelayLGain, right_out * mDelayRGain };
    return output;
  }

private:
  const double mMaxDelayMS;
  int mMaxDelay;
  DelayLine mDelayL;
  DelayLine mDelayR;

  T mDelayLGain{ 0.5 };
  T mDelayRGain{ 0.5 };
  double mDelayLTimeMS;
  double mDelayRTimeMS;
  int mDelayLTime;
  int mDelayRTime;
  T mFeedback{ 0. };
};