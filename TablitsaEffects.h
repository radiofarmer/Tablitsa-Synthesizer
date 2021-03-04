#pragma once

#include "SignalProcessing.h"
#include "Filter.h"
#include "Modulators.h"

#include <mutex>

template<typename T>
class Effect
{
public:
  Effect(T sampleRate) : mSampleRate(sampleRate) {}

  virtual T Process(T s);

  virtual void ProcessStereo(T* s);

  virtual void SetSampleRate(T sampleRate)
  {
    mSampleRate = sampleRate;
  }

  virtual void SetMix(T wet) { mMix = wet; }

  virtual void SetParam1(T value) {}
  virtual void SetParam2(T value) {}
  virtual void SetParam3(T value) {}
  virtual void SetParam4(T value) { SetMix(value / (T)100.); }
  virtual void SetParam5(T value) {}
  virtual void SetParam6(T value) {}

protected:
  T mSampleRate;
  T mMix{ 0. };
};

#define DELAY_TEMPODIV_VALIST "1/64", "1/32", "1/16T", "1/16", "1/16D", "1/8T", "1/8", "1/8D", "1/4", "1/4D", "1/2", "1/1"
#define TABLITSA_MAX_DELAY_MS 10000.
#define TABLITSA_MAX_DELAY_SAMP (int)480000

template<typename T, int MaxDelay=TABLITSA_MAX_DELAY_SAMP>
class DelayEffect final : public Effect<T>
{
  enum EChannels
  {
    kLeft = 0,
    kRight,
    kMono
  };

public:
  enum ETempoDivision
  {
    k64th = 0,   // 1 sixty fourth of a beat
    k32nd,       // 1 thirty second of a beat
    k16thT,      // 1 sixteenth note triplet
    k16th,       // 1 sixteenth note
    k16thD,      // 1 dotted sixteenth note
    k8thT,       // 1 eighth note triplet
    k8th,        // 1 eighth note
    k8thD,       // 1 dotted eighth note
    k4th,        // 1 quater note a.k.a 1 beat @ 4/4
    k4thD,       // 1 dotted beat @ 4/4
    k2th,        // 2 beats @ 4/4
    k1,          // 1 bar @ 4/4
    kNumDivisions
  };

  DelayEffect(T sampleRate) :
    DelayEffect(sampleRate, nullptr)
  {
  }

  DelayEffect(double sampleRate, ModMetronome* metronome=nullptr) :
    Effect<T>(sampleRate),
    mMaxDelayMS(TABLITSA_MAX_DELAY_MS),
    mMaxDelay(MaxDelay),
    mDelayLTime(mMaxDelay / 2),
    mDelayRTime(mMaxDelay / 2),
    mDelayLTimeMS(mMaxDelayMS / 2),
    mDelayRTimeMS(mMaxDelayMS / 2),
    mDelayLBeats(1.),
    mDelayRBeats(1.),
    mMetronome(metronome)
  {
    Reset();
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
  virtual void SetParam3(T value) override { SetFeedback(value / (T)100.); }
  virtual void SetParam4(T value) override { SetGain(value / (T)100.); }
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
    mDelayLTime = std::min(mDelayLTime, TABLITSA_MAX_DELAY_SAMP);
    mDelayRTime = std::min(mDelayRTime, TABLITSA_MAX_DELAY_SAMP);
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

  void Reset()
  {
    mDelayL.reset();
    mDelayR.reset();
  }

private:
  const T mMaxDelayMS;
  int mMaxDelay;

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

  DelayLine<T, TABLITSA_MAX_DELAY_SAMP>  mDelayL;
  DelayLine<T, TABLITSA_MAX_DELAY_SAMP> mDelayR;
  ModMetronome* mMetronome; // For tempo sync
};

template<typename T>
class SampleAndHold : public Effect<T>
{
  struct xor128
  {
    uint32_t w, x, y, z;
  };

public:
  SampleAndHold(T sampleRate) : Effect<T>(sampleRate) {}

  void SetParam1(T value) override { SetRateMS(value); }
  void SetParam2(T value) override { SetDecay(value / (T)100.); }
  void SetParam3(T value) override { SetJitter(value / (T)100.); }

  T Process(T s) override
  {
    T out = mHold;
    mHold *= mDecay;
    if (mSampleCounter++ >= mRateAdj)
    {
      out = mHold = s;
      mSampleCounter = 0;
      mRateAdj = mRate + (XOR_Shift(mRandGen) >> mJitter);
    }
    out = mMix * (out - s);
    return s + out;
  }

  void ProcessStereo(T* s) override
  {
    T out = mHold;
    mHold *= mDecay;
    if (mSampleCounter++ >= mRateAdj)
    {
      out = mHold = (s[0] + s[1]) / (T)2;
      mSampleCounter = 0;
      mRateAdj = mRate + (XOR_Shift(mRandGen) >> mJitter);
    }
    s[0] += mMix * (out - s[0]);
    s[1] += mMix * (out - s[1]);
  }

  void SetRateMS(T rate)
  {
    mRate = static_cast<int>(rate / 1000. * mSampleRate);
    mRateAdj = mRate + (XOR_Shift(mRandGen) >> mJitter);
  }

  void SetDecay(T decay)
  {
    mDecay = 1. - decay;
  }

  void SetJitter(T noise)
  {
    mJitter = static_cast<int>((T)31 * ((T)1 - noise * 0.25));
  }

  static inline uint32_t XOR_Shift(xor128& state)
  {
    const uint32_t first = state.w;
    uint32_t last = state.z;
    state.z = state.y;
    state.y = state.x;
    state.x = first;
    last ^= first << 11;
    last ^= first >> 8;
    return state.w = last ^ first ^ (first >> 19);
  }

protected:
  T mHold{ 0. };
  T mDecay{ 0. };
  unsigned int mRate{ 1 }; // samples
  unsigned int mSampleCounter{ 0 };
  unsigned int mJitter{ 0 };
  unsigned int mRateAdj{ 1 };
  xor128 mRandGen{ 0xAF31, 0x1234, 0xFF2E, 0xDCBA };
};

#define WAVESHAPE_TYPES "Sine", "Parabolic", "Hyp. Tan.", "Soft Clip", "Hard Clip"

enum EWaveshaperMode
{
  kWaveshapeSine,
  kWaveshapeParabola,
  kWaveshapeTanh,
  kWaveshapeSoft,
  kWaveshapeHard,
  kNumWaveshaperModes
};

template<typename T>
class Waveshaper : public Effect<T>
{
  static constexpr T piOver2{ (T)1.57079632679 };
  static constexpr T pi{ (T)3.14159265359 };
  static constexpr T MaxGain{ (T)3. };

  static inline T SineShaper(T x, T gain)
  {
    x *= gain;
    T x2 = std::copysign(x, piOver2 - x);
    return std::copysign(x2 - x2 * x2 * x2 / (T)6 + x2 * x2 * x2 * x2 * x2 / (T)120, x) / gain;
  }

  static inline T TanhShaper(T x, T gain)
  {
    x *= gain;
    return std::tanh(x);
  }

  static inline T ParabolicShaper(T x, T gain)
  {
    x *= gain;
    return copysign(x * x, x);
  }

  static inline T SoftClipShaper(T x, T gain)
  {
    return SoftClip<T>(x, gain) / gain;
  }

  static inline T HardClipShaper(T x, T gain)
  {
    return std::copysign(std::min(std::abs(x) * gain, 1. / (gain * 2)), x);
  }

public:
  Waveshaper(T sampleRate, EWaveshaperMode mode = kWaveshapeSine) : Effect<T>(sampleRate), mShaperMode(mode) {}

  virtual void SetParam1(T value) override { SetMode(static_cast<EWaveshaperMode>(value)); }
  virtual void SetParam2(T value) override { SetGain(value / (T)100); }

  T Process(T s) override
  {
    return mShaperFunc(s, mGain);
  }

  void ProcessStereo(T* s) override
  {
    std::lock_guard<std::mutex> lg(mFuncMutex);
    s[0] += mMix * (mShaperFunc(s[0], mGain) - s[0]);
    s[1] += mMix * (mShaperFunc(s[1], mGain) - s[1]);
  }

  void SetMode(EWaveshaperMode mode)
  {
    std::lock_guard<std::mutex> lg(mFuncMutex); // To prevent calling an empty `std::function`
    switch (mode)
    {
    case kWaveshapeSine:
      mShaperFunc = &Waveshaper::SineShaper;
      break;
    case kWaveshapeParabola:
      mShaperFunc = &Waveshaper::ParabolicShaper;
      break;
    case kWaveshapeTanh:
      mShaperFunc = &Waveshaper::TanhShaper;
      break;
    case kWaveshapeSoft:
      mShaperFunc = &Waveshaper::SoftClipShaper;
      break;
    case kWaveshapeHard:
      mShaperFunc = &Waveshaper::HardClipShaper;
      break;
    default:
      mShaperFunc = &Waveshaper::SineShaper;
      break;
    }
  }

  void SetGain(T gain)
  {
    mGain = 1. + gain * MaxGain;
  }

protected:
  EWaveshaperMode mShaperMode;
  T mGain{ (T)1 };
  std::function<T(T, T)> mShaperFunc{ &Waveshaper::SineShaper };
  std::mutex mFuncMutex;
};