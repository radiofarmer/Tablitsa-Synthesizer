#include "FastOversampler.h"

/* 1x (no) Oversampling */
template<>
void UpsampleBlock<1>(FastOversampler<sample>& oversampler, const sample* inputs, const int nFrames)
{
  // Set output source to inputs
  oversampler.mOutputSource->Set(inputs, nFrames);
}

template<>
void DownsampleBlock<1>(FastOversampler<sample>& oversampler, sample* outputs, const int nFrames)
{
  memcpy(outputs, oversampler.mOutputSource->Get(), nFrames * sizeof(outputs[0]));
}

/* 2x Oversampling */

template<>
void UpsampleBlock<2>(FastOversampler<sample>& oversampler, const sample* inputs, const int nFrames)
{
  oversampler.mUpsampler2x.process_block(oversampler.mUp2Buf.Get(), inputs, nFrames);
  oversampler.SetOutputSource(&oversampler.mUp2Buf);
}

template<>
void DownsampleBlock<2>(FastOversampler<sample>& oversampler, sample* outputs, const int nFrames)
{
  oversampler.mDownsampler2x.process_block(outputs, oversampler.mOutputSource->Get(), nFrames);
}

/* 4x Oversampling */

template<>
void UpsampleBlock<4>(FastOversampler<sample>& oversampler, const sample* inputs, const int nFrames)
{
  UpsampleBlock<2>(oversampler, inputs, nFrames);
  oversampler.mUpsampler4x.process_block(oversampler.mUp4Buf.Get(), oversampler.mOutputSource->Get(), nFrames * 2);
  oversampler.SetOutputSource(&oversampler.mUp4Buf);
}

template<>
void DownsampleBlock<4>(FastOversampler<sample>& oversampler, sample* outputs, const int nFrames)
{
  oversampler.mDownsampler4x.process_block(oversampler.mDown4Buf.Get(), oversampler.mOutputSource->Get(), nFrames * 2);
  oversampler.mDownsampler2x.process_block(outputs, oversampler.mDown4Buf.Get(), nFrames);
}