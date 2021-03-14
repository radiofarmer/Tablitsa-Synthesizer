#pragma once

#include "radiofarmer_config.h"
#include "fastmath_tables.h"

BEGIN_DSP_NAMESPACE

#define FAST_MATH_MACROS
#define UNITBIT32 1572864.  /* 3*2^19; bit 32 has place value 1 */
#define HIOFFSET 1
#define LOWOFFSET 0

union Union64 {
  double d;
  int i[2];
} __declspec(align(16));

inline int to_int_in_range(const double x, const int range)
{
  const int rangeM1 = range - 1;

  union Union64 u;
  u.d = UNITBIT32;
  const int normhipart = u.i[HIOFFSET];

  u.d = x * rangeM1 + (double)UNITBIT32;
  return u.i[HIOFFSET] & (range - 1); 
}

// Convert a double between -range and range to an integer between 0 and range, with the fractional part in `frac`
inline int to_int_in_range(const double x, const int range, double& frac)
{
  const int rangeM1 = range - 1;

  union Union64 u;
  u.d = UNITBIT32;
  const int normhipart = u.i[HIOFFSET];

  u.d = x * rangeM1 + (double)UNITBIT32;
  const int result = u.i[HIOFFSET] & (range - 1);
  u.i[HIOFFSET] = normhipart;
  frac = u.d - UNITBIT32;

  return result;
}

/* Returns sin(x * pi/2) for values between -1. and 1. */
inline sample_t fast_sin_norm(double x)
{
  const int idx{ to_int_in_range(0.5 * (x + 1.), 512) };
  return SineTable512[idx];
}

inline sample_t fast_cos_norm(double x)
{
  const int idx{ to_int_in_range(0.5 * (x + 1.), 512) };
  return CosineTable512[idx];
}

inline sample_t fast_tanh(sample_t x)
{
  return x / (1. + std::abs(2 * x));
}

END_DSP_NAMESPACE