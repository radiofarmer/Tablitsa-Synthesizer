#pragma once
// Required for Oversampler.h:
#include "IPlugMidi.h"
#include <assert.h>
#include "Oversampler.h"

using namespace iplug;
using namespace hiir;

static constexpr double coeffs2x[12] = { 0.036681502163648017, 0.13654762463195794, 0.27463175937945444, 0.42313861743656711, 0.56109869787919531, 0.67754004997416184, 0.76974183386322703, 0.83988962484963892, 0.89226081800387902, 0.9315419599631839, 0.96209454837808417, 0.98781637073289585 };
static constexpr double coeffs4x[4] = { 0.041893991997656171, 0.16890348243995201, 0.39056077292116603, 0.74389574826847926 };
static constexpr double coeffs8x[3] = { 0.055748680811302048, 0.24305119574153072, 0.64669913119268196 };
static constexpr double coeffs16x[2] = { 0.10717745346023573, 0.53091435354504557 };

template<typename T>
struct FastOversampler
{
  /* Upsamplers */
  Upsampler2xFPU<12, T> mUpsampler2x;
  Upsampler2xFPU<4, T> mUpsampler4x;
  Upsampler2xFPU<3, T> mUpsampler8x;
  Upsampler2xFPU<2, T> mUpsampler16x;

  /* DownSamplers */
  Downsampler2xFPU<12, T> mDownsampler2x;
  Downsampler2xFPU<4, T> mDownsampler4x;
  Downsampler2xFPU<3, T> mDownsampler8x;
  Downsampler2xFPU<2, T> mDownsampler16x;

  WDL_TypedBuf<T> mUp2Buf;
  WDL_TypedBuf<T> mDown2Buf;
  WDL_TypedBuf<T> mUp4Buf;
  WDL_TypedBuf<T> mDown4Buf;
  WDL_TypedBuf<T> mUp8Buf;
  WDL_TypedBuf<T> mDown8Buf;
  WDL_TypedBuf<T> mUp16Buf;
  WDL_TypedBuf<T> mDown16Buf;

  WDL_TypedBuf<T>* mInputSource;
  WDL_TypedBuf<T>* mOutputSource;

  FastOversampler()
  {
    mUpsampler2x.set_coefs(coeffs2x);
    mDownsampler2x.set_coefs(coeffs2x);

    mUpsampler4x.set_coefs(coeffs4x);
    mDownsampler4x.set_coefs(coeffs4x);

    mUpsampler8x.set_coefs(coeffs8x);
    mDownsampler8x.set_coefs(coeffs8x);

    mUpsampler16x.set_coefs(coeffs16x);
    mDownsampler16x.set_coefs(coeffs16x);

    mOutputSource = &mUp2Buf;
  }

  void ResizeBuffers(const int blockSize)
  {
    mUp2Buf.Resize(blockSize * 2);
    mDown2Buf.Resize(blockSize * 2);
    mUp4Buf.Resize(blockSize * 4);
    mDown4Buf.Resize(blockSize * 4);
    mUp8Buf.Resize(blockSize * 8);
    mDown8Buf.Resize(blockSize * 8);
    mUp16Buf.Resize(blockSize * 16);
    mDown16Buf.Resize(blockSize * 16);
  }

  void SetOutputSource(WDL_TypedBuf<T>* src)
  {
    mOutputSource = src;
  }
};

typedef double sample;

template<int Factor>
void UpsampleBlock(FastOversampler<sample>& oversampler, const sample* inputs, sample* outputs, const int nFrames);

template<int Factor>
void DownsampleBlock(FastOversampler<sample>& oversampler, sample* outputs, const int nFrames);