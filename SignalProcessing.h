#pragma once

#include "Filter.h"

class DelayLine
{
public:
  DelayLine(const int length) : mLength(length)
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

  inline void push(const double s)
  {
    mRead = mWrite;
    mBuffer[mWrite++] = s;
    if (mWrite == mLength)
      mWrite = 0;
  }

  inline double at(const int offset = 0)
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

template<typename T>
class Effect
{
public:
  Effect(double sampleRate) : mSampleRate(sampleRate) {}

  virtual T Process(T s) = 0;

  virtual void SetSampleRate(T sampleRate)
  {
    mSampleRate = sampleRate;
  }

protected:
  double mSampleRate;
};

#define DELAY_TEMPODIV_VALIST "1/64", "1/32", "1/16T", "1/16", "1/16D", "1/8T", "1/8", "1/8D", "1/4", "1/4D", "1/2", "1/1"

template<typename T>
class DelayEffect final : public Effect<T>
{
  enum EChannels
  {
    kLeft=0,
    kRight,
    kMono
  };

public:
  DelayEffect(double sampleRate, double maxDelayMS = 5000.) :
    Effect<T>(sampleRate),
    mMaxDelayMS(maxDelayMS),
    mMaxDelay(static_cast<int>((mMaxDelayMS / 1000. + 1.) * mSampleRate)),
    mDelayL(mMaxDelay),
    mDelayR(mMaxDelay),
    mDelayLTime(mMaxDelay / 2),
    mDelayRTime(mMaxDelay / 2),
    mDelayLTimeMS(mMaxDelayMS / 2),
    mDelayRTimeMS(mMaxDelayMS / 2),
    mDelayLBeats(1.),
    mDelayRBeats(1.)
  {

  }

  void SetDelayMS(T timeMS, int channel)
  {
    if (channel == kLeft || channel == kMono)
      mDelayLTimeMS = timeMS;
    if (channel == kRight || channel == kMono)
      mDelayRTimeMS = timeMS;
    CalculateDelaySamples();
  }

  void SetDelayTempo(T beatFraction, int channel, T tempo=120.)
  {
    if (channel == kLeft || channel == kMono)
      mDelayLBeats = beatFraction;
    if (channel == kRight || channel == kMono)
      mDelayRBeats = beatFraction;
    mBPM = tempo;
    CalculateDelaySamples();
  }

  void SetTempoSync(bool sync)
  {
    mTempoSync = sync;
    CalculateDelaySamples();
  }

  void SetSampleRate(T sampleRate) override
  {
    Effect<T>::SetSampleRate(sampleRate);
    CalculateDelaySamples();
  }

  void CalculateDelaySamples()
  {
    if (mTempoSync)
    {
      mDelayLTime = static_cast<int>(mDelayLBeats * mBPM / 60. * mSampleRate);
      mDelayRTime = static_cast<int>(mDelayRBeats * mBPM / 60. * mSampleRate);
    }
    else
    {
      mDelayLTime = static_cast<int>(mDelayLTimeMS / 1000. * mSampleRate);
      mDelayRTime = static_cast<int>(mDelayRTimeMS / 1000. * mSampleRate);
    }
  }

  void SetFeedback(const T fb)
  {
    mFeedback = fb;
  }

  void SetGain(const T gain)
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

  T* ProcessStereo(const T sl, const T sr)
  {
    const T left_out = mDelayL[mDelayLTime];
    const T right_out = mDelayR[mDelayRTime];
    mDelayL.push(sl + left_out * mFeedback);
    mDelayR.push(sr + right_out * mFeedback);
    T output[2]{ left_out * mDelayLGain, right_out * mDelayRGain };
    return output;
  }

  T* ProcessStereo(T inputs[2])
  {
    const T left_out = mDelayL[mDelayLTime];
    const T right_out = mDelayR[mDelayRTime];
    mDelayL.push(inputs[0] + left_out * mFeedback);
    mDelayR.push(inputs[1] + right_out * mFeedback);
    inputs[0] = left_out * mDelayLGain;
    inputs[1] = right_out* mDelayRGain;
    return inputs;
  }

private:
  const double mMaxDelayMS;
  int mMaxDelay;
  DelayLine mDelayL;
  DelayLine mDelayR;

  T mDelayLGain{ 0.5 };
  T mDelayRGain{ 0.5 };
  // Millisecond delay times
  T mDelayLTimeMS;
  T mDelayRTimeMS;
  // Tempo-sync delay times
  T mDelayLBeats;
  T mDelayRBeats;
  // Sample Delay Times
  int mDelayLTime;
  int mDelayRTime;
  // Delay time mode
  bool mTempoSync{ false };
  T mBPM{ 120. };
  T mFeedback{ 0. };
};
