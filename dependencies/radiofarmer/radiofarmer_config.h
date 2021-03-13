#pragma once

#include <vectorclass.h>
#include <algorithm>

#define BEGIN_DSP_NAMESPACE namespace radiofarmer {
#define END_DSP_NAMESPACE }

BEGIN_DSP_NAMESPACE

#ifdef EFFECT_SAMPLE_FLOAT
typedef float sample_t;
typedef Vec4f sample_v;
typedef Vec4i sample_vi;
#else
typedef double sample_t;
typedef Vec4d sample_v;
typedef Vec4q sample_vi;
#endif

constexpr sample_t DEFAULT_SRATE = (sample_t)44100;

template<class T=sample_t>
struct StereoSample
{
  T l;
  T r;

  StereoSample(T l=0., T r=0.) : l(l), r(r) {}

  // Copy contstructors
  StereoSample& operator=(const StereoSample& s)
  {
    l = s.l;
    r = s.r;
    return *this;
  }

  StereoSample& operator=(const T s)
  {
    l = r = s;
    return *this;
  }

  StereoSample& operator=(const int i)
  {
    l = r = (T)i;
    return *this;
  }

  StereoSample operator+(const StereoSample& s)
  {
    return StereoSample(l + s.l, r + s.r);
  }

  StereoSample operator-(const StereoSample& s)
  {
    return StereoSample(l - s.l, r - s.r);
  }

  StereoSample& operator+=(const StereoSample& s)
  {
    l += s.l;
    r += s.r;
    return *this;
  }

  StereoSample& operator-=(const StereoSample& s)
  {
    l -= s.l;
    r -= s.r;
    return *this;
  }

  StereoSample& operator*(const T s)
  {
    l *= s;
    r *= s;
    return *this;
  }

  StereoSample& operator/(const T s)
  {
    l /= s;
    r /= s;
    return *this;
  }
};

template<typename T=sample_t>
struct StereoBuffer
{
  T* const L;
  T* const R;
  const unsigned int N;

  /** Create a buffer using explicit pointers to each channel.
  @param l Pointer to the left channel buffer.
  @param r Pointer to the right channel buffer.
  @param nFrames Length of the buffer (in samples). If this is set to -1, the length will be 2^32 - 1, and arithmetic operations will not be available. */
  StereoBuffer(T* l, T* r, const int nFrames=-1) : L(l), R(r), N(nFrames) {}

  /* Create a buffer based on the address of the first sample, and the buffer length */
  StereoBuffer(T* buf, const int nFrames) : L(buf), R(buf + nFrames), N(nFrames) {}

  StereoBuffer& operator+(T* buf)
  {
    assert(static_cast<signed int>(nFrames) >= 0 && "Cannot add to buffer with an undefined length");
    for (int i{ 0 }; i < nFrames; ++i)
    {
      L[i] += buf[0][i];
      R[i] += buf[1][i];
    }
    return *this;
  }

  StereoBuffer& operator-(T* buf)
  {
    assert(static_cast<signed int>(nFrames) >= 0 && "Cannot subtract from buffer with an undefined length");
    for (int i{ 0 }; i < nFrames; ++i)
    {
      L[i] -= buf[0][i];
      R[i] -= buf[1][i];
    }
    return *this;
  }

  StereoBuffer& operator*(T* buf)
  {
    assert(static_cast<signed int>(nFrames) >= 0 && "Cannot multiply buffer with an undefined length");
    for (int i{ 0 }; i < nFrames; ++i)
    {
      L[i] *= buf[0][i];
      R[i] *= buf[1][i];
    }
    return *this;
  }

  StereoBuffer& operator/(T* buf)
  {
    assert(static_cast<signed int>(nFrames) >= 0 && "Cannot divide buffer with an undefined length");
    for (int i{ 0 }; i < nFrames; ++i)
    {
      L[i] /= buf[0][i];
      R[i] /= buf[1][i];
    }
    return *this;
  }

  template<class V=sample_v>
  StereoSample<V> GetVector(const int offset = 0)
  {
    return StereoSample<V>(V().load(&L[offset]), V().load(&R[offset]));
  }

  StereoSample<T> operator[](const int idx)
  {
    return StereoSample<T>(L[idx], R[idx]);
  }
};

END_DSP_NAMESPACE