#pragma once

#include "SignalProcessing.h"
#include "Filter.h"
#include "Modulators.h"

#include "vectormath_trig.h"
#include "vectormath_hyp.h"

#include <mutex>

template<class T>
struct StereoSample
{
  T l;
  T r;
};

template<typename T, class V=Vec4d>
class Effect
{
public:
  Effect(T sampleRate) : mSampleRate(sampleRate), mVectorSize(V().size())
  {
  }

  virtual T Process(T s) { return s; }
  virtual void ProcessStereo(T* s) {}

  virtual V __vectorcall Process(V& s) { return s; }
  virtual void ProcessStereo_Vector(StereoSample<V>& s) {}

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
  const int mVectorSize;
};

#define DELAY_TEMPODIV_VALIST "1/64", "1/32", "1/16T", "1/16", "1/16D", "1/8T", "1/8", "1/8D", "1/4", "1/4D", "1/2", "1/1"
#define TABLITSA_MAX_DELAY_SAMP (int)524288
#define TABLITSA_MAX_DELAY_MS 524288. / 48000. * 1000.

template<typename T, int MaxDelay=TABLITSA_MAX_DELAY_SAMP, class V = Vec4d>
class DelayEffect final : public Effect<T, V>
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
  T Process(T s) override
  {
    T left_out = mDelayL[mDelayLTime];
    T right_out = mDelayR[mDelayRTime];
    mDelayL.push(s + left_out * mFeedback);
    mDelayR.push(s + right_out * mFeedback);
    return left_out * mDelayLGain + right_out * mDelayRGain;
  }

  // Returns the mixed signal
  void ProcessStereo(T* s)
  {
    const T left_out = mDelayL[mDelayLTime];
    const T right_out = mDelayR[mDelayRTime];
    mDelayL.push(s[0] + left_out * mFeedback);
    mDelayR.push(s[1] + right_out * mFeedback);
    s[0] += left_out * mDelayLGain;
    s[1] += right_out * mDelayRGain;
  }

  void ProcessStereo_Vector(StereoSample<V>& s) override
  {
    V l_out = mDelayL.v_at(mDelayLTime);
    V r_out = mDelayR.v_at(mDelayRTime);
    mDelayL.push<V>(s.l);
    mDelayR.push<V>(s.r);
    s.l += l_out;
    s.r += r_out;
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

  T mFeedback{ 0. };
  T mBPM{ 120. };
  
  // Delay time mode
  bool mTempoSync{ false };

  DelayLine<T, TABLITSA_MAX_DELAY_SAMP>  mDelayL;
  DelayLine<T, TABLITSA_MAX_DELAY_SAMP> mDelayR;
  ModMetronome* mMetronome; // For tempo sync
};

template<typename T, class V = Vec4d>
class SampleAndHold : public Effect<T, V>
{
  typedef typename deduce_vector_from<V>::int_vec Vi;

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
    mHold = mHold + mDecay * mRateAdj / mSampleCounter * (out - mHold);
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
    mHold += mDecay * (s[0] - mHold);
    if (mSampleCounter++ >= mRateAdj)
    {
      out = mHold = (s[0] + s[1]) / (T)2;
      mSampleCounter = 0;
      mRateAdj = mRate + (XOR_Shift(mRandGen) >> mJitter);
    }
    s[0] += mMix * (out - s[0]);
    s[1] += mMix * (out - s[1]);
  }

  void ProcessStereo_Vector(StereoSample<V>& s) override
  {
    StereoSample<V> out{ mHold_v.l, mHold_v.r };
    mHold_v.l += mDecay * (s.l - mHold_v.l);
    mHold_v.r += mDecay * (s.r - mHold_v.r);
    if (mSampleCounter >= mRateAdj)
    {
      out.l = mHold_v.l = s.l;
      out.r = mHold_v.r = s.r;
      mSampleCounter = 0;
      mRateAdj = mRate + (XOR_Shift(mRandGen) >> mJitter);
    }
    else
      mSampleCounter += deduce_vector_from<V>::v_size;
    s.l += mMix * (out.l - s.l);
    s.r += mMix * (out.r - s.r);
  }

  void SetRateMS(T rate)
  {
    mRate = static_cast<int>(rate / 1000. * mSampleRate);
    mRateAdj = mRate + (XOR_Shift(mRandGen) >> mJitter);
  }

  void SetDecay(T decay)
  {
    mDecay = decay / mRate;
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

  V mSampleCounter_v = V(0., 1., 2., 3);
  StereoSample<V> mHold_v{ V(0.), V(0.) };
};

#define WAVESHAPE_TYPES "Sine", "Parabolic", "Cubic", "Hyp. Tan.", "Soft Clip", "Hard Clip"

enum EWaveshaperMode
{
  kWaveshapeSine,
  kWaveshapeParabola,
  kWaveshapeCubic,
  kWaveshapeTanh,
  kWaveshapeSoft,
  kWaveshapeHard,
  kNumWaveshaperModes
};

template<typename T, class V = Vec4d>
class Waveshaper : public Effect<T, V>
{
  static constexpr T piOver2{ (T)1.57079632679 };
  static constexpr T pi{ (T)3.14159265359 };

  /* Singleton Versions */

  static inline T SineShaper(T x, const T gain, const T gainCeil=1.)
  {
    x *= gain;
    T x2 = std::copysign(x, piOver2 - x);
    return std::copysign(x2 - x2 * x2 * x2 / (T)6 + x2 * x2 * x2 * x2 * x2 / (T)120, x) * 0.9;
  }

  static inline T ParabolicShaper(T x, const T gain, const T gainCeil = 1.)
  {
    x *= gain;
    return SoftClipShaper(copysign(x * x, x), 1., gainCeil);
  }

  static inline T TanhShaper(T x, const T gain, const T gainCeil = 1.)
  {
    return std::tanh(x * gain) * gainCeil;
  }

  static inline T CubicShaper(T x, const T gain, const T gainCeil = 1.)
  {
    x *= gain;
    return SoftClipShaper(1.5 * x + 0.5 * x * x * x, 1., gainCeil);
  }

  static inline T SoftClipShaper(T x, const T gain, const T gainCeil = 1.)
  {
    return SoftClip<T>(x, gain) * gainCeil;
  }

  static inline T HardClipShaper(T x, const T gain, const T gainCeil = 1.)
  {
    return std::copysign(std::min(std::abs(x) * gain, gainCeil), x);
  }

  /* Vector versions */

  static inline V __vectorcall VSineShaper(const V& x, const T gain, const T gainCeil = 1.)
  {
    return sin(x * gain) * gainCeil;
  }

  static inline V __vectorcall VParabolicShaper(const V& x, const T gain, const T gainCeil = 1.)
  {
    const V x2 = pow(x * gain, 2) * gain * sign(x);
    return x2 * gainCeil;
  }

  static inline V __vectorcall VCubicShaper(const V& x, const T gain, const T gainCeil = 1.)
  {
    const V x3 = pow(x * gain, 3) * gain;
    return x3 * gainCeil;
  }

  static inline V __vectorcall VTanhShaper(const V& x, const T gain, const T gainCeil = 1.)
  {
    return tanh(x * gain) * gainCeil;
  }

  static inline V VSoftClipShaper(const V& x, const T gain, const T gainCeil = 1.)
  {
    return SoftClip<V, T>(x, gain) * gainCeil;
  }

  static inline V __vectorcall VHardClipShaper(const V& x, const T gain, const T gainCeil = 1.)
  {
    return sign(x) * min(abs(x * gain), V(gainCeil));
  }

public:
  Waveshaper(T sampleRate, T maxGain=2., EWaveshaperMode mode = kWaveshapeSine) :
    Effect<T>(sampleRate), mMaxGain(maxGain), mShaperMode(mode) {}

  virtual void SetParam1(T value) override
  {
    SetMode(static_cast<EWaveshaperMode>(value + 0.01));
  }
  virtual void SetParam2(T value) override
  {
    SetGain(value / (T)100);
    SetThreshold((T) - value / 100.)
  }

  T Process(T s) override
  {
    return mShaperFunc(s, mGain, mThresh);
  }

  inline T DoProcess(T s)
  {
    return mMix * (mShaperFunc(s, mGain, mThresh) - s);
  }

  inline V __vectorcall DoProcess_Vector(V& s)
  {
    return mMix * (mShaperFunc_Vector(s, mGain, mThresh) - s);
  }

  void ProcessStereo(T* s) override
  {
    std::lock_guard<std::mutex> lg(mFuncMutex);
    s[0] += DoProcess(s[0]);
    s[1] += DoProcess(s[1]);
  }

  void ProcessStereo_Vector(StereoSample<V>& s)
  {
    std::lock_guard<std::mutex> lg(mFuncMutex);
    s.l += DoProcess_Vector(s.l);
    s.r += DoProcess_Vector(s.r);
  }

  inline void SetMode(EWaveshaperMode mode)
  {
    std::lock_guard<std::mutex> lg(mFuncMutex); // To prevent calling an empty `std::function`
    switch (mode)
    {
    case kWaveshapeParabola:
    {
      mShaperFunc = &Waveshaper::ParabolicShaper;
      mShaperFunc_Vector = &Waveshaper::VParabolicShaper;
      break;
    }
    case kWaveshapeCubic:
    {
      mShaperFunc = &Waveshaper::CubicShaper;
      mShaperFunc_Vector = &Waveshaper::VCubicShaper;
      break;
    }
    case kWaveshapeTanh:
    {
      mShaperFunc = &Waveshaper::TanhShaper;
      mShaperFunc_Vector = &Waveshaper::VTanhShaper;
      break;
    }
    case kWaveshapeSoft:
    {
      mShaperFunc = &Waveshaper::SoftClipShaper;
      mShaperFunc_Vector = &Waveshaper::VSoftClipShaper;
      break;
    }
    case kWaveshapeHard:
    {
      mShaperFunc = &Waveshaper::HardClipShaper;
      mShaperFunc_Vector = &Waveshaper::VHardClipShaper;
      break;
    }
    case kWaveshapeSine:
    default:
    {
      mShaperFunc = &Waveshaper::SineShaper;
      mShaperFunc_Vector = &Waveshaper::VSineShaper;
      break;
    }
    }
  }

  inline void SetGain(const T gain)
  {
    mGain = (T)1. + gain * mMaxGain;
  }

  inline void SetThreshold(T thresh)
  {
    mThresh = thresh;
  }

protected:
  const T mMaxGain; // Can be used to normalize inputs, since a waveshaper's behavior is amplitude-dependent
  T mGain{ (T)1 };
  T mThresh{ 0.5 };
  EWaveshaperMode mShaperMode;
  std::function<T(T, const T, const T)> mShaperFunc{ &Waveshaper::SineShaper };
  std::function<V(const V&, const T, const T)> mShaperFunc_Vector{ &Waveshaper::VSineShaper };
  std::mutex mFuncMutex;
};

/* 3-Band EQ */
template<typename T, class V = Vec4d>
class EQ3Effect : public Effect<T, V>
{
  static constexpr T pi{ 3.14159265359 };
  static constexpr T twoPi{ 6.28318530718 };
  static constexpr T denorm_fix = (T)(1.0 / 4294967295.0);

  struct EQState
  {
    // Low-shelf frequency and poles
    T lf;
    T lfp0;
    T lfp1;
    T lfp2;
    T lfp3;

    // High-shelf frequency and poles
    T hf;
    T hfp0;
    T hfp1;
    T hfp2;
    T hfp3;

    // Gain
    T lg;
    T mg;
    T hg;
  };

public:
  EQ3Effect(T sampleRate) : Effect<T>(sampleRate)
  {
    memset(&mStateL, static_cast<int>((T)0), sizeof(mStateL));
    memset(&mStateR, static_cast<int>((T)0), sizeof(mStateR));
    SetMidFreq(0.25);
  }

  void SetSampleRate(T sampleRate) override
  {
    SetMidFreq(std::min(mMidFreq * mSampleRate / sampleRate, 0.99));
    Effect<T>::SetSampleRate(sampleRate);
  }

  void SetParam1(T value) override { SetLowGain(value); }
  void SetParam2(T value) override { SetMidGain(value); }
  void SetParam3(T value) override { SetMidFreq(value * (T)0.00249); } // Accepts values between 0 and 100
  void SetParam4(T value) override { SetHighGain(value); }

  inline void SetLowGain(T lg) { mStateL.lg = lg; mStateR.lg = lg; }
  inline void SetMidGain(T mg) { mStateL.mg = mg; mStateR.mg = mg; }
  inline void SetHighGain(T hg) { mStateL.hg = hg; mStateR.hg = hg; }

  inline void SetMidFreq(T freqNorm)
  {
    T lowFreq = std::max(freqNorm - mHalfMidBand, 0.01);
    T highFreq = std::min(freqNorm + mHalfMidBand, 0.99);
    mStateL.lf = (T)2 * std::sin(pi * lowFreq);
    mStateL.hf = (T)2 * std::sin(pi * highFreq);
    mStateR.lf = (T)2 * std::sin(pi * lowFreq);
    mStateR.hf = (T)2 * std::sin(pi * highFreq);
    mMidFreq = freqNorm;
  }

  inline T DoProcess(EQState& state, DelayLine<T, 4>& z, T s)
  {
    T l, m, h;
    // Lowpass Filter
    state.lfp0 += state.lf * (s - state.lfp0) + denorm_fix;
    state.lfp1 += state.lf * (state.lfp0 - state.lfp1);
    state.lfp2 += state.lf * (state.lfp1 - state.lfp2);
    state.lfp3 += state.lf * (state.lfp2 - state.lfp3);
    l = state.lfp3;

    // Highpass Filter
    state.hfp0 += state.hf * (s - state.hfp0) + denorm_fix;
    state.hfp1 += state.hf * (state.hfp0 - state.hfp1);
    state.hfp2 += state.hf * (state.hfp1 - state.hfp2);
    state.hfp3 += state.hf * (state.hfp2 - state.hfp3);
    h = z[0] - state.hfp3;

    // Bandpass Filter
    m = z[0] - (l + h);

    // Apply gain
    l *= state.lg;
    m *= state.mg;
    h *= state.hg;

    // Update delay line
    z.push(s);

    return l + m + h;
  }

  void ProcessStereo(T* s) override
  {
    s[0] = DoProcess(mStateL, mZ0, s[0]);
    s[1] = DoProcess(mStateR, mZ1, s[1]);
  }

protected:
  EQState mStateL;
  EQState mStateR;
  T mMidFreq{ 1. };
  T mHalfMidBand{ 0.05 }; // Half the proportion of the normalized frequency range occupied by the mid band

  DelayLine<T, 4> mZ0;
  DelayLine<T, 4> mZ1;
};