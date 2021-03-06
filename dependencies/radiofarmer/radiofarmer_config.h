#pragma once

#include <vectorclass.h>

#define BEGIN_DSP_NAMESPACE namespace radiofarmer {
#define END_DSP_NAMESPACE }

BEGIN_DSP_NAMESPACE

#ifdef EFFECT_SAMPLE_FLOAT
typedef float sample_t;
typedef Vec4f sample_v;
#else
typedef double sample_t;
typedef Vec4d sample_v;
#endif

END_DSP_NAMESPACE