#pragma once

#include "radiofarmerDSP.h"
#include "FastOversampler.h"

#include "Filter.h"
#include "Modulators.h"
#include "TablitsaOscillators.h"

#include "vectormath_trig.h"
#include "vectormath_hyp.h"

#include <mutex>

using namespace radiofarmer;

template<typename T, class V=Vec4d>
class Effect
{
public:
  Effect(T sampleRate, int oversampling=1) :
    mSampleRate(sampleRate),
    mOversampling(oversampling),
    mOsFreqScalar(1. / static_cast<T>(oversampling)),
    mVectorSize(V().size())
  {
  }

  virtual void SetSampleRate(T sampleRate, int oversampling=1)
  {
    mSampleRate = sampleRate;
    mOversampling = oversampling;
    mOsFreqScalar = 1. / static_cast<T>(oversampling);
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

  virtual void SetMix(T wet) { mMix = wet; }

  virtual void SetParam1(T value) {}
  virtual void SetParam2(T value) {}
  virtual void SetParam3(T value) {}
  virtual void SetParam4(T value) { SetMix(value); }
  virtual void SetParam5(T value) {}
  virtual void SetParam6(T value) {}

  // TODO make the individual SetParam...() functions non-virtual to improve performance.
  // Or just move the entire contents of those functions into this one (only for voice effects)
  virtual void SetContinuousParams(const T p1, const T p2, const T p3, const T p4, const T pitch)
  {
    SetParam1(p1);
    SetParam2(p2);
    SetParam3(p3);
    SetParam4(p4);
  }

protected:
  T mSampleRate;
  int mOversampling;
  T mOsFreqScalar;
  T mMix{ 0. };
  const int mVectorSize;
};

/* Coefficient Modulation ("Super Ring") */
template<typename T, class V=Vec4d>
class CMEffect : public Effect<T, V>
{
public:
  CMEffect(T sampleRate) : Effect<T, V>(sampleRate) {}

  void SetContinuousParams(const T p1, const T p2, const T p3, const T p4, const T pitch)
  {
    mModDepth = p1;
    mModRate = 440. * std::pow(2., pitch + p2); // p2 is scaled in octaves
    mOsc.SetFreqCPS(mModRate);
    mOsc.SetPhaseOffset(p3);
    mMix = p4;
  }

  T Process(T s) override
  {
    return mCMFilters.Process(s, mModDepth);
  }

  void ProcessBlock(T* inputs, T* outputs, const int nFrames)
  {
    T oscVals[4];
    for (int i{ 0 }; i < nFrames; i += 4)
    {
      mOsc.Process_Vector().store(oscVals);
#pragma clang loop unroll(full)
      for (int ii{ 0 }; ii < 4; ++ii)
      {
        outputs[i + ii] = inputs[i + ii] + mMix * (mCMFilters.Process(inputs[i + ii], mModDepth * oscVals[ii]) - inputs[i + ii]);
      }
    }
  }

private:
  T mModDepth{ 0. };
  T mModRate;
  VectorOscillator<T> mOsc;
  ModulatedAllpass<5> mCMFilters;
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

  void SetSampleRate(T sampleRate, int oversampling=1) override
  {
    Effect<T, V>::SetSampleRate(sampleRate, oversampling);
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

/* Distortion */
#define DISTORTION_TYPES "Fuzz", "Asymm.", "Soft Clip 1", "Soft Clip 2", "Hard Clip"
enum EDistortion
{
  kFuzz,
  kAsymmetric,
  kSoftClip1,
  kSoftClip2,
  kHardClip,
  kNumDistortionModes
};

template<typename T, class V = Vec4d>
class DistortionEffect final : public Effect<T, V>
{
private:
  static constexpr int BufferLength{ 1024 };
  static constexpr T BufferScale{ 1. / BufferLength };
  static constexpr T FilterModCenterFreqHz{ 600. };
  static constexpr T MaxFilterModHz{ 100. };

  inline T FuzzDistortion(T x) const
  {
    const T x_abs = std::abs(x);
    return x * (x_abs + mGain) / (x_abs * (x_abs + mGain - 1.) + 1.);
  }

  inline T ExpDistortion(T x) const
  {
    return std::copysign(1. - std::exp(-1. * (std::abs(x) * mGain)), x);
  }

  inline T AsymmDistortion(T x) const
  {
    x *= mGain;
    if (x >= 0.)
    {
      T tanhx = std::tanh(x);
      return tanhx * tanhx;
    }
    else
      return std::atan(x);
  }

  inline T SoftClipDistortion(T x) const
  {
    return SoftClip<T, 5>(x * mGain);
  }

  // Bram de Jong soft saturation algorithm (from MusicDSP via https://github.com/johannesmenzel/SRPlugins/wiki/DSP-ALGORITHMS-Saturation)
  inline T SoftClipBdJDistortion(T x) const
  {
    x *= mGain;
    if (std::abs(x) <= 1.)
      return (x / (x * x + 1.)) * 2.;
    else
      return std::copysign(1., x);
  }

  inline T HardClipDistortion(T x) const
  {
    return std::clamp(x * mGain, -1., 1.);
  }

public:
  DistortionEffect(T sampleRate) : Effect<T, V>(sampleRate) {}

  void SetContinuousParams(const T p1, const T p2, const T p3, const T p4, const T pitch) override
  {
    mGain = std::pow(10., p1 * 2.);
    mFilterMod = p2 * mMaxFilterMod;
    mColorFilter.SetCutoff(p3 * 0.1 * mOsFreqScalar);
    mMix = p4;
  }

  void SetSampleRate(T sampleRate, int oversampling=1) override
  {
    Effect<T, V>::SetSampleRate(sampleRate, oversampling); // Sets `mOversampling` and `mOsFreqScalar`
    mNoiseFilter.SetSampleRate(sampleRate * oversampling);
    mColorFilter.SetSampleRate(sampleRate * oversampling);
    mFilterOsc.SetSampleRate(sampleRate * oversampling);
    mFilterModCenterFreqNorm = FilterModCenterFreqHz / sampleRate;
    mMaxFilterMod = std::min(mFilterModCenterFreqNorm * 0.95, MaxFilterModHz / sampleRate);
  }

  inline T DoProcess(T s, const int type)
  {
    switch (type)
    {
    case kFuzz:
      return FuzzDistortion(s);
    case kAsymmetric:
      return AsymmDistortion(s);
    case kSoftClip1:
      return SoftClipDistortion(s);
    case kSoftClip2:
      return SoftClipBdJDistortion(s);
    case kHardClip:
      return HardClipDistortion(s);
    default:
      return s;
    }
  }

  void ProcessBlock(T* inputs, T* outputs, const int nFrames) override
  {
    mNoiseFilter.SetCutoff(mFilterModCenterFreqNorm + mFilterMod * mFilterOsc.x * ((T)(std::rand() / RAND_MAX) * 0.1 + 1.));
    for (int i{ 0 }; i < nFrames; ++i)
    {
#ifdef DISTORTION_GAIN_NORM
      const T ipt_filtered = mColorFilter.ProcessBP(inputs[i]);

      // Update average input volume
      mAvgIn += BufferScale * (std::abs(inputs[i]) - mInHist.last());
      mInHist.push(std::abs(inputs[i]));

      // Calculate output (and update average output volume)
      const T out = FuzzDistortion(mNoiseFilter.ProcessAP(inputs[i] + ipt_filtered));
      mAvgOut += BufferScale * (std::abs(out) - mOutHist.last());
      mOutHist.push(std::abs(out));

      mFilterOsc.Step();
      outputs[i] = inputs[i] + mMix * (out - inputs[i]);
#else
      const T ipt_filtered = mColorFilter.ProcessBP(inputs[i]);
      // Calculate output
      const T out = FuzzDistortion(mNoiseFilter.ProcessAP(inputs[i] + ipt_filtered));
      mFilterOsc.Step();
      outputs[i] = inputs[i] + mMix * (out - inputs[i]);
#endif
    }
  }

private:
  T mGain;
  T mAvgIn{ 0.25 };
  T mAvgOut{ 1. };
  int mType{ EDistortion::kFuzz };
  DelayLine<BufferLength> mInHist;
  DelayLine<BufferLength> mOutHist;

  TwoPoleTPTFilter mNoiseFilter{ DEFAULT_SAMPLE_RATE, 0.25, 0.5 };
  TwoPoleTPTFilter mColorFilter{ DEFAULT_SAMPLE_RATE, 0.25 };

  // Noise/fuzz filter controls
  T mFilterMod{ 0. };
  sine_osc_nd mFilterOsc{ 100., DEFAULT_SAMPLE_RATE };
  T mFilterModCenterFreqNorm{ 600. / DEFAULT_SAMPLE_RATE };
  T mMaxFilterMod{ std::min(mFilterModCenterFreqNorm * 0.9, MaxFilterModHz / DEFAULT_SAMPLE_RATE) };
};


/* 3-Band EQ */
template<typename T, class V = Vec4d>
class EQ3Effect final : public Effect<T, V>
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

  void SetSampleRate(T sampleRate, int oversampling) override
  {
    SetMidFreq(std::min(mMidFreq * mSampleRate / sampleRate, 0.99));
    Effect<T, V>::SetSampleRate(sampleRate * oversampling);
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

private:
  EQState mStateL;
  EQState mStateR;
  T mMidFreq{ 1. };
  T mHalfMidBand{ 0.05 }; // Half the proportion of the normalized frequency range occupied by the mid band

  DelayLine<4> mZ0;
  DelayLine<4> mZ1;
};

template<typename T, class V=Vec4d>
class ReverbEffect final : public Effect<T, V>
{
  static constexpr T MaxDelayMS = 100.;
  static constexpr T MinDelayMS = 5.;
  static constexpr T MinFeedback = 0.25;
  static constexpr T FeedbackRange = 0.99 - MinFeedback;

public:
  ReverbEffect(T sampleRate, T maxDelayMS = 50., T minDelayMS = 10., T maxFeedback = 0.75, T minFeedback = 0.65, T gain = 0.5) :
    Effect<T, V>(sampleRate),
    mReverb(sampleRate, maxDelayMS, minDelayMS, maxFeedback, minFeedback, gain)
  {}

  void SetSampleRate(T sampleRate, int oversampling) override
  {
    Effect<T, V>::SetSampleRate(sampleRate, oversampling);
    mReverb.SetSampleRate(mSampleRate * oversampling);
  }

  void SetParam1(T value) override
  {
    mReverb.SetDelay((MinDelayMS + 20.) + value * (MaxDelayMS - MinDelayMS), MinDelayMS, true);
  }
  void SetParam2(T value) override
  {
    mReverb.SetFeedback(MinFeedback + (1. - value) * FeedbackRange, MinFeedback, true);
  }
  void SetParam3(T value) override
  {
    mReverb.SetGain((1. - value) * 0.95);
  }
  void SetParam4(T value) override
  {
    mReverb.SetMix(value);
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

private:
  CascadeReverb<6, 2> mReverb;
};

template<typename T, class V=Vec4d>
class Reverb2Effect final : public Effect<T, V>
{
public:
  Reverb2Effect(T sampleRate) : Effect<T, V>(sampleRate), mReverb(sampleRate)
  {
  }

  void SetParam1(T value) override { mReverb.SetDiffusion(value * 0.6); mReverb.SetEarlyReflectionsLevel(0.05 + value * 0.1); }
  void SetParam2(T value) override { mReverb.SetDamping(value * value * 0.2);  }
  void SetParam3(T value) override { mReverb.SetColor(value * 0.7); }
  void SetParam4(T value) override { mReverb.SetMixLevel(value); }

  void ProcessStereo(StereoSample<T>& s)
  {
    mReverb.ProcessStereo(s);
  }

private:
  UFDNReverb mReverb;
};


template<typename T, class V = Vec4d>
class SampleAndHold final : public Effect<T, V>
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

  void SetContinousParams(const T p1, const T p2, const T p3, const T p4, const T pitch)
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

private:
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

template<typename T, class V = Vec4d>
class Texturizer final : public Effect<T, V>
{

public:
  Texturizer(T sampleRate) : Effect<T, V>(sampleRate)
  {
    mPk[0].SetPeakGain(4.);
    mPk[1].SetPeakGain(4.);
  }

  void SetSampleRate(T sampleRate, int oversampling=1) override
  {
    Effect<T, V>::SetSampleRate(sampleRate, oversampling);
    for (int i{ 0 }; i < 2; ++i)
      mAP[i].SetSampleRate(sampleRate);
  }

  void SetContinuousParams(const T p1, const T p2, const T p3, const T p4, const T pitch) override
  {
    mFc = p1 * mOsFreqScalar;
    mQ = p2;
    for (int i{ 0 }; i < 2; ++i)
    {
      mPk[i].SetCutoff(p3 * mOsFreqScalar);
    }
    for (int i{ 0 }; i < 2; ++i)
    {
      mPk[i].SetPeakGain(p4);
    }
    SetFilterCoefs();

    mMix = p4;
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
    return mAP[0].Process4(s) + 0.1 * sin(mPk[0].Process4(s) * 63. * mQ);
  }

  void ProcessStereo(StereoSample<T>& s)
  {
    s.l = mAP[0].Process(s.l) + 0.2 * std::sin(mPk[0].Process(s.l) * 63. * mQ);
    s.r = mAP[1].Process(s.r) + 0.2 * std::sin(mPk[1].Process(s.r) * 63. * mQ);
  }

  void ProcessBlock(T* inputs, T* outputs, const int nFrames) override
  {
    for (int i{ 0 }; i < nFrames; i += 4)
    {
      V s4 = Process(V().load(inputs + i));
      s4.store(outputs + i);
    }
  }

private:
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
class Waveshaper final : public Effect<T, V>
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

  void SetContinuousParams(const T p1, const T p2, const T p3, const T p4, const T pitch) override
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

private:
  const T mMaxGain; // Can be used to normalize inputs, since a waveshaper's behavior is amplitude-dependent
  const T mMaxGainCeil;
  T mGain{ (T)1 };
  T mThresh{ 0.5 };
  EWaveshaperMode mShaperMode;
  std::function<T(T, const T)> mShaperFunc{ &Waveshaper::SineShaper };
  std::function<V(const V&, const T)> mShaperFunc_Vector{ &Waveshaper::VSineShaper };
  std::mutex mFuncMutex;
};  