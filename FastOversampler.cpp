#include "FastOversampler.h"

template<>
void UpsampleBlock<2>(FastOversampler<sample>& oversampler, const sample* inputs, sample* outputs, const int nFrames)
{
  oversampler.mUpsampler2x.process_block(oversampler.mUp2Buf.Get(), inputs, nFrames);
  oversampler.SetOutputSource(&oversampler.mUp2Buf);
}

template<>
void DownsampleBlock<2>(FastOversampler<sample>& oversampler, sample* outputs, const int nFrames)
{
  oversampler.mDownsampler2x.process_block(outputs, oversampler.mOutputSource->Get(), nFrames);
}


template<>
void UpsampleBlock<4>(FastOversampler<sample>& oversampler, const sample* inputs, sample* outputs, const int nFrames)
{
  UpsampleBlock<2>(oversampler, inputs, outputs, nFrames);
  oversampler.mUpsampler4x.process_block(oversampler.mUp4Buf.Get(), oversampler.mOutputSource->Get(), nFrames * 2);
  oversampler.SetOutputSource(&oversampler.mUp4Buf);
}

template<>
void DownsampleBlock<4>(FastOversampler<sample>& oversampler, sample* outputs, const int nFrames)
{
  oversampler.mDownsampler4x.process_block(oversampler.mDown4Buf.Get(), oversampler.mOutputSource->Get(), nFrames * 2);
  oversampler.mDownsampler2x.process_block(outputs, oversampler.mDown4Buf.Get(), nFrames);
}