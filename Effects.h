#pragma once

#include "SignalProcessing.h"
#include "Filter.h"
#include "Modulators.h"

template<typename T>
class Effect
{
public:
  Effect(T sampleRate) : mSampleRate(sampleRate) {}

  virtual T Process(T s) { return s; }

  virtual void ProcessStereo(T* s) {}

  virtual void SetSampleRate(T sampleRate)
  {
    mSampleRate = sampleRate;
  }

  virtual void SetMix(T wet) { mWetMix = wet; }

  virtual void SetParam1(T value) {}
  virtual void SetParam2(T value) {}
  virtual void SetParam3(T value) {}
  virtual void SetParam4(T value) { SetMix((T)value / 100.); }
  virtual void SetParam5(T value) {}
  virtual void SetParam6(T value) {}

protected:
  T mSampleRate;
  T mWetMix{ 0. };
};

#define DELAY_TEMPODIV_VALIST "1/64", "1/32", "1/16T", "1/16", "1/16D", "1/8T", "1/8", "1/8D", "1/4", "1/4D", "1/2", "1/1"

template<typename T>
class DelayEffect final : public Effect<T>
{
  enum EChannels
  {
    kLeft = 0,
    kRight,
    kMono
  };

public:
  DelayEffect(T sampleRate, T maxDelayMS = 5000.) :
    DelayEffect(sampleRate, maxDelayMS, nullptr)
  {
  }

  DelayEffect(double sampleRate, double maxDelayMS = 5000., ModMetronome* metronome=nullptr) :
    Effect<T>(sampleRate),
    mMaxDelayMS(maxDelayMS),
    mMaxDelay(static_cast<int>((mMaxDelayMS / 1000. + 1.)* mSampleRate)),
    mDelayL(mMaxDelay),
    mDelayR(mMaxDelay),
    mDelayLTime(mMaxDelay / 2),
    mDelayRTime(mMaxDelay / 2),
    mDelayLTimeMS(mMaxDelayMS / 2),
    mDelayRTimeMS(mMaxDelayMS / 2),
    mDelayLBeats(1.),
    mDelayRBeats(1.),
    mMetronome(metronome)
  {

  }
  // Helper functions
  virtual void SetParam1(T value) override
  {
    if (mTempoSync)
    {
      double qnScalar = LFO<T>::GetQNScalar(static_cast<LFO<T>::ETempoDivision>(Clip((int)value, 0, (int)LFO<T>::ETempoDivision::kNumDivisions)));
      double qnPerMeasure = 4. / mMetronome->mTSDenom * mMetronome->mTSNum;
      SetDelayTempo(1. / qnScalar / qnPerMeasure, 0);
    }
    else
      SetDelayMS(value, 0);
  }
  virtual void SetParam2(T value) override
  {
    if (mTempoSync)
    {
      double qnScalar = LFO<T>::GetQNScalar(static_cast<LFO<T>::ETempoDivision>(Clip((int)value, 0, (int)LFO<T>::ETempoDivision::kNumDivisions)));
      double qnPerMeasure = 4. / mMetronome->mTSDenom * mMetronome->mTSNum;
      SetDelayTempo(1. / qnScalar / qnPerMeasure, 1);
    }
    else
      SetDelayMS(value, 1);
  }
  virtual void SetParam3(T value) override { SetFeedback((T)value / 100.); }
  virtual void SetParam4(T value) override { SetGain((T)value / 100.); }
  virtual void SetParam5(T value) override { SetTempoSync(value > 0.5); }

  void SetDelayMS(T timeMS, int channel)
  {
    if (channel == kLeft || channel == kMono)
      mDelayLTimeMS = timeMS;
    if (channel == kRight || channel == kMono)
      mDelayRTimeMS = timeMS;
    CalculateDelaySamples();
  }

  void SetDelayTempo(T beatFraction, int channel, T tempo = 120.)
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
    if (mMetronome)
      mBPM = mMetronome->mTempo;
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

  // Returns only the wet signal
  inline T Process(T s) override
  {
    T left_out = mDelayL[mDelayLTime];
    T right_out = mDelayR[mDelayRTime];
    mDelayL.push(s + left_out * mFeedback);
    mDelayR.push(s + right_out * mFeedback);
    return left_out * mDelayLGain + right_out * mDelayRGain;
  }

  // Returns the mixed signal
  inline void ProcessStereo(T inputs[2])
  {
    const T left_out = mDelayL[mDelayLTime];
    const T right_out = mDelayR[mDelayRTime];
    mDelayL.push(inputs[0] + left_out * mFeedback);
    mDelayR.push(inputs[1] + right_out * mFeedback);
    inputs[0] += left_out * mDelayLGain;
    inputs[1] += right_out * mDelayRGain;
  }

private:
  const double mMaxDelayMS;
  int mMaxDelay;
  DelayLine mDelayL;
  DelayLine mDelayR;
  ModMetronome* mMetronome; // For tempo sync

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

template<typename T>
class SampleAndHold : public Effect<T>
{
public:
  SampleAndHold(T sampleRate) : Effect(sampleRate) {}

  void SetParam1(T value) override { SetRateMS(value); }
  void SetParam2(T value) override { SetDecay((T)value / 100.); }
  void SetParam3(T value) override { SetNoise((T)value / 100.); }

  T Process(T s) override
  {
    T out = mHold;
    if (mSampleCounter++ >= mRate)
    {
      out = mHold = s;
      mSampleCounter = 0;
    }
    return s + mWetMix * (out - s);
  }

  void ProcessStereo(T* s) override
  {
    s[0] = Process(s[0]);
    s[1] = s[0];
  }

  void SetRateMS(T rate)
  {
    mRate = static_cast<int>(rate / 1000. * mSampleRate);
  }

  void SetDecay(T decay)
  {
    mDecay = decay;
  }

  void SetNoise(T noise)
  {
    mNoise = noise;
  }

protected:
  int mRate{ 1 }; // samples
  int mSampleCounter{ 0 };
  T mHold{ 0. };
  T mDecay{ 0. };
  T mNoise{ 0. };
};