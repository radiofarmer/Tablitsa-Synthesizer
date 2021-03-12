#pragma once

#include "radiofarmerDSP.h"
#include "FastOversampler.h"

#include "Filter.h"
#include "Modulators.h"

#include "vectormath_trig.h"
#include "vectormath_hyp.h"

#include <mutex>

using namespace radiofarmer;

template<typename T, class V=Vec4d>
class Effect
{
public:
  Effect(T sampleRate) : mSampleRate(sampleRate), mVectorSize(V().size())
  {
  }

  virtual T Process(T s) { return s; }
  virtual V __vectorcall Process(V s) { return s; }

  // Scalar and single StereoSample processing
  virtual void ProcessStereo(T* s) {}
  virtual void ProcessStereo(StereoSample<T>& s) {}

  // Vector Processing
  virtual V __vectorcall Process(V& s) { return s; }
  virtual void ProcessStereo_Vector(StereoSample<V>& s) {}


  // Block Processing
  virtual void ProcessBlock(T* inputs, T* outputs, const int nFrames)
  {
    for (int i{ 0 }; i < nFrames; ++i)
    {
      outputs[i] = Process(inputs[i]);
    }
  }

  // Parameters
  virtual void SetSampleRate(T sampleRate)
  {
    mSampleRate = sampleRate;
  }

  virtual void SetMix(T wet) { mMix = wet; }

  virtual void SetParam1(T value) {}
  virtual void SetParam2(T value) {}
  virtual void SetParam3(T value) {}
  virtual void SetParam4(T value) { SetMix(value); }
  virtual void SetParam5(T value) {}
  virtual void SetParam6(T value) {}

  // TODO make the individual SetParam...() functions non-virtual to improve performance.
  // Or just move the entire contents of those functions into this one (only for voice effects)
  virtual void SetContinuousParams(const T p1, const T p2, const T p3, const T p4)
  {
    SetParam1(p1);
    SetParam2(p2);
    SetParam3(p3);
    SetParam4(p4);
  }

  void ResizeBuffers(const int blockSize)
  {
    mOversampler.ResizeBuffers(blockSize);
  }

protected:
  T mSampleRate;
  T mMix{ 0. };
  const int mVectorSize;

  FastOversampler<T> mOversampler;
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
    Effect<T, V>(sampleRate),
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
  virtual void SetParam3(T value) override { SetFeedback(value); }
  virtual void SetParam4(T value) override { SetGain(value); }
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
    Effect<T, V>::SetSampleRate(sampleRate);
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
  void ProcessStereo(T* s) override
  {
    const T left_out = mDelayL[mDelayLTime];
    const T right_out = mDelayR[mDelayRTime];
    mDelayL.push(s[0] + left_out * mFeedback);
    mDelayR.push(s[1] + right_out * mFeedback);
    s[0] += left_out * mDelayLGain;
    s[1] += right_out * mDelayRGain;
  }

  void ProcessStereo(StereoSample<T>& s) override
  {
    const T left_out = mDelayL[mDelayLTime];
    const T right_out = mDelayR[mDelayRTime];
    mDelayL.push(s.l + left_out * mFeedback);
    mDelayR.push(s.r + right_out * mFeedback);
    s.l += left_out * mDelayLGain;
    s.r += right_out * mDelayRGain;
  }

  void ProcessStereo_Vector(StereoSample<V>& s) override
  {
    V l_out = mDelayL.v_at(mDelayLTime);
    V r_out = mDelayR.v_at(mDelayRTime);
    mDelayL.push(s.l);
    mDelayR.push(s.r);
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

  DelayLine<TABLITSA_MAX_DELAY_SAMP>  mDelayL;
  DelayLine<TABLITSA_MAX_DELAY_SAMP> mDelayR;
  ModMetronome* mMetronome; // For tempo sync
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
  EQ3Effect(T sampleRate) : Effect<T, V>(sampleRate)
  {
    memset(&mStateL, static_cast<int>((T)0), sizeof(mStateL));
    memset(&mStateR, static_cast<int>((T)0), sizeof(mStateR));
    SetMidFreq(0.25);
  }

  void SetSampleRate(T sampleRate) override
  {
    SetMidFreq(std::min(mMidFreq * mSampleRate / sampleRate, 0.99));
    Effect<T, V>::SetSampleRate(sampleRate);
  }

  void SetParam1(T value) override { SetLowGain(value); }
  void SetParam2(T value) override { SetMidGain(value); }
  void SetParam3(T value) override { SetMidFreq(value * (T)0.249); }
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

  inline T DoProcess(EQState& state, DelayLine<4>& z, T s)
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

  void ProcessStereo(StereoSample<T>& s) override
  {
    s.l = DoProcess(mStateL, mZ0, s.l);
    s.r = DoProcess(mStateR, mZ1, s.r);
  }

protected:
  EQState mStateL;
  EQState mStateR;
  T mMidFreq{ 1. };
  T mHalfMidBand{ 0.05 }; // Half the proportion of the normalized frequency range occupied by the mid band

  DelayLine<4> mZ0;
  DelayLine<4> mZ1;
};

template<typename T, class V = Vec4d>
class SampleAndHold : public Effect<T, V>
{
  typedef typename deduce_vector_from<V>::int_vec Vi;

  struct xor128
  {
    uint32_t w, x, y, z;
  };

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

public:
  SampleAndHold(T sampleRate) : Effect<T, V>(sampleRate)
  {
    // These must be dynamically allocated for proper alignment, to avoid heap corruption
    mSampleCounter_v = new V(0., 1., 2., 3);
    mHold_v = new StereoSample<V>(V(0.), V(0.));
  }

  void SetParam1(T value) override { SetRateMS(value); }
  void SetParam2(T value) override { SetDecay(value); }
  void SetParam3(T value) override { SetJitter(value); }

  void SetContinousParams(T param1, T param2, T param3, T param4)
  {
    SetRateMS(param1);
    SetDecay(param2);
    SetJitter(param3);

    mMix = param4;
  }

  T Process(T s) override
  {
    T out = mHold;
    mHold += mDecay * (s - mHold);
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

  void ProcessStereo(StereoSample<T>& s) override
  {
    StereoSample<T> out = mStereoHold;
    mStereoHold += (s - mStereoHold) * mDecay;
    if (mSampleCounter++ >= mRateAdj)
    {
      out = mStereoHold = s;
      mSampleCounter = 0;
      mRateAdj = mRate + (XOR_Shift(mRandGen) >> mJitter);
    }
    s += (out - s) * mMix;
  }

  void ProcessStereo_Vector(StereoSample<V>& s) override
  {
    StereoSample<V> out{ mHold_v->l, mHold_v->r };
    mHold_v->l += mDecay * (s.l - mHold_v->l);
    mHold_v->r += mDecay * (s.r - mHold_v->r);
    if (mSampleCounter >= mRateAdj)
    {
      out.l = mHold_v->l = s.l;
      out.r = mHold_v->r = s.r;
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

  ~SampleAndHold()
  {
    delete mSampleCounter;
    delete mHold_v;
  }

protected:
  T mHold{ 0. };
  T mDecay{ 0. };
  StereoSample<T> mStereoHold{ 0., 0. };
  unsigned int mRate{ 1 }; // samples
  unsigned int mSampleCounter{ 0 };
  unsigned int mJitter{ 0 };
  unsigned int mRateAdj{ 1 };
  xor128 mRandGen{ 0xAF31, 0x1234, 0xFF2E, 0xDCBA };

  V* mSampleCounter_v;
  StereoSample<V>* mHold_v;
};

template<typename T, class V=Vec4d>
class ReverbEffect : public Effect<T, V>
{
  static constexpr T MaxDelayMS = 100.;
  static constexpr T MinDelayMS = 5.;
  static constexpr T MinFeedback = 0.65;
  static constexpr T FeedbackRange = 0.99 - MinFeedback;

public:
  ReverbEffect(T sampleRate, T maxDelayMS = 50., T minDelayMS = 10., T maxFeedback = 0.75, T minFeedback = 0.65, T gain = 0.5) :
    Effect<T, V>(sampleRate),
    mReverb(sampleRate, maxDelayMS, minDelayMS, maxFeedback, minFeedback, gain)
  {}

  void SetSampleRate(T sampleRate) override
  {
    Effect<T, V>::SetSampleRate(sampleRate);
    mReverb.SetSampleRate(mSampleRate);
  }

  void SetParam1(T value) override
  {
    mReverb.SetDelay(value * (MaxDelayMS - MinDelayMS), MinDelayMS, true);
  }
  void SetParam2(T value) override
  {
    mReverb.SetFeedback(MinFeedback + value * FeedbackRange, MinFeedback, true);
  }
  void SetParam4(T value) override
  {
    mReverb.SetGain(value);
  }

  T Process(T s) override
  {
    return mReverb.Process(s);
  }

  void ProcessStereo(T* s) override
  {
    StereoSample<T> s2{ s[0], s[1] };
    mReverb.ProcessStereo(s2);
    s[0] = s2.l; s[1] = s2.r;
  }

  void ProcessStereo(StereoSample<T>& s) override
  {
    mReverb.ProcessStereo(s);
  }

protected:
  CascadeReverb<6, 2> mReverb;
};

template<typename T, class V=Vec4d>
class Reverb2Effect : public Effect<T, V>
{
public:
  Reverb2Effect(T sampleRate) : Effect<T, V>(sampleRate), mReverb(sampleRate)
  {
  }

  void SetParam1(T value) override { mReverb.SetDiffusion(value * 0.75); mReverb.SetEarlyReflectionsLevel(0.5 - value * 0.25); }
  void SetParam2(T value) override { mReverb.SetDamping(value * 0.2);  }
  void SetParam3(T value) override { mReverb.SetColor(value * 0.7); }
  void SetParam4(T value) override { mReverb.SetMixLevel(value); }

  void ProcessStereo(StereoSample<T>& s)
  {
    mReverb.ProcessStereo(s);
  }

protected:
  UFDNReverb mReverb;
};

template<typename T, class V = Vec4d>
class Texturizer : public Effect<T, V>
{

public:
  Texturizer(T sampleRate) : Effect<T, V>(sampleRate)
  {
    mPk[0].SetPeakGain(4.);
    mPk[1].SetPeakGain(4.);
  }

  void SetSampleRate(T sampleRate)
  {
    Effect<T, V>::SetSampleRate(sampleRate);
    for (int i{ 0 }; i < 2; ++i)
      mAP[i].SetSampleRate(sampleRate);
  }

  void SetContinuousParams(T param1, T param2, T param3, T param4) override
  {
    mFc = param1;
    mQ = param2;
    for (int i{ 0 }; i < 2; ++i)
    {
      mPk[i].SetCutoff(param3);
    }
    for (int i{ 0 }; i < 2; ++i)
    {
      mPk[i].SetPeakGain(param4);
    }
    SetFilterCoefs();

    mMix = param4;
  }

  void SetFilterCoefs() {
    mAP[0].SetCoefs(mFc, mG);
    mAP[1].SetCoefs(mFc, mG);
  }

  T Process(T s) override
  {
    return mAP[0].Process(s) + 0.1 * std::sin(mPk[0].Process(s) * 63. * mQ);
  }

  V __vectorcall Process(V s) override
  {
    return mAP[0].Process4(s) + 0.1 * sin(mPk[0].Process4(s) * 6.28 * mQ);
  }

  void ProcessStereo(StereoSample<T>& s)
  {
    s.l = mAP[0].Process(s.l) + 0.2 * std::sin(mPk[0].Process(s.l) * 63. * mQ);
    s.r = mAP[1].Process(s.r) + 0.2 * std::sin(mPk[1].Process(s.r) * 63. * mQ);
  }

  void ProcessBlock(T* inputs, T* outputs, const int nFrames) override
  {
    UpsampleBlock<4>(mOversampler, inputs, outputs, nFrames);
    for (int i{ 0 }; i < nFrames * 4; i += 4)
    {
      V s4 = Process(V().load(mOversampler.mOutputSource->Get() + i));
      s4.store(mOversampler.mOutputSource->Get() + i);
    }
    DownsampleBlock<4>(mOversampler, outputs, nFrames);
  }

protected:
  DelayLine<4> mZ;
  T mFc{ 1. };
  T mQ{ 0. };
  T mG{ 1. };
  AllpassLadder<4> mAP[2];
  PeakOnePole mPk[2];
};

#define WAVESHAPE_TYPES "Sine", "Parabolic", "Hyp. Tan.", "Soft Clip"

enum EWaveshaperMode
{
  kWaveshapeSine,
  kWaveshapeParabola,
  kWaveshapeTanh,
  kWaveshapeSoft,
  kNumWaveshaperModes
};

template<typename T, class V = Vec4d>
class Waveshaper : public Effect<T, V>
{
  static constexpr T piOver2{ (T)1.57079632679 };
  static constexpr T pi{ (T)3.14159265359 };

  /* Singleton Versions */

  static inline T SineShaper(T x, const T gain)
  {
    x *= gain;
    const T x3 = x * x * x;
    return x - 0.16666666666 * x3 + 0.008333333333333 * x3 * x * x;
  }

  static inline T ParabolicShaper(T x, const T gain)
  {
    x = SoftClipShaper(x, gain);
    return std::copysign(x * x, x);
  }

  static inline T TanhShaper(T x, const T gain)
  {
    return std::tanh(x * gain);
  }

  static inline T SoftClipShaper(T x, const T gain)
  {
    return SoftClip<T, 5>(x * gain);
  }

  static inline T HardClipShaper(T x, const T gain)
  {
    return std::copysign(std::min(std::abs(x) * gain), x);
  }

  /* Vector versions */

  static inline V __vectorcall VSineShaper(const V& x, const T gain)
  {
    return sin(x * gain);
  }

  static inline V __vectorcall VParabolicShaper(const V& x, const T gain)
  {
    V x_clip = VSoftClipShaper(x, gain);
    return pow(x_clip, 2) * sign(x);
  }

  static inline V __vectorcall VTanhShaper(const V& x, const T gain)
  {
    return tanh(x * gain);
  }

  static inline V __vectorcall VSoftClipShaper(const V& x, const T gain)
  {
    return SoftClip<V, 5>(x * gain);
  }

  static inline V __vectorcall VHardClipShaper(const V& x, const T gain)
  {
    return sign(x) * min(abs(x * gain));
  }

public:
  Waveshaper(T sampleRate, T maxGain = 2., EWaveshaperMode mode = kWaveshapeSine) :
    Effect<T, V>(sampleRate),
    mMaxGain(maxGain), mMaxGainCeil(1. / mMaxGain),
    mShaperMode(mode) {}

  virtual void SetParam1(T value) override
  {
    SetMode(static_cast<EWaveshaperMode>(value + 0.01));
  }
  virtual void SetParam2(T value) override
  {
    SetGain(value);
    SetThreshold(1. - (mMaxGain - 1.) * (value * mMaxGainCeil));
  }

  void SetContinuousParams(T param1, T param2, T param3, T param4) override
  {
    SetMode(static_cast<EWaveshaperMode>(param1 + 0.01));

    SetGain(param2);
    SetThreshold(1. - (mMaxGain - 1.) * (param2 * mMaxGainCeil));

    mMix = param4;
  }

  T Process(T s) override
  {
    return s + mMix * mThresh * (DoProcess(s) - s);
  }

  inline T DoProcess(T s)
  {
    return mShaperFunc(s, mGain);
  }

  inline V __vectorcall DoProcess_Vector(V& s)
  {
    return (mShaperFunc_Vector(s, mGain) - s);
  }

  void ProcessStereo(T* s) override
  {
    //std::lock_guard<std::mutex> lg(mFuncMutex);
    s[0] += mMix * mThresh * (DoProcess(s[0]) - s[0]);
    s[1] += mMix * mThresh * (DoProcess(s[1]) - s[1]);
  }

  void ProcessStereo(StereoSample<T>& s) override
  {
    //std::lock_guard<std::mutex> lg(mFuncMutex);
    s.l += mMix * mThresh * (DoProcess(s.l) - s.l);
    s.r += mMix * mThresh * (DoProcess(s.r) - s.r);
  }

  void ProcessStereo_Vector(StereoSample<V>& s)
  {
    //std::lock_guard<std::mutex> lg(mFuncMutex);
    s.l += mMix * mThresh * (DoProcess_Vector(s.l) - s.l);
    s.r += mMix * mThresh * (DoProcess_Vector(s.r) - s.r);
  }

  inline void SetMode(EWaveshaperMode mode)
  {
    //std::lock_guard<std::mutex> lg(mFuncMutex); // To prevent calling an empty `std::function`
    switch (mode)
    {
    case kWaveshapeParabola:
    {
      mShaperFunc = &Waveshaper::ParabolicShaper;
      mShaperFunc_Vector = &Waveshaper::VParabolicShaper;
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
  const T mMaxGainCeil;
  T mGain{ (T)1 };
  T mThresh{ 0.5 };
  EWaveshaperMode mShaperMode;
  std::function<T(T, const T)> mShaperFunc{ &Waveshaper::SineShaper };
  std::function<V(const V&, const T)> mShaperFunc_Vector{ &Waveshaper::VSineShaper };
  std::mutex mFuncMutex;
};  