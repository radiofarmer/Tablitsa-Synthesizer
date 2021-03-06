#pragma once

#include <vectorclass.h>

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
};

END_DSP_NAMESPACE