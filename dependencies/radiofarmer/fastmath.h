#pragma once

#include "radiofarmer_config.h"
#include "fastmath_tables.h"

#include <cmath>

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

/* Non-deterministic sine oscillator */
struct sine_osc_nd
{
  double k; // omega = sqrt(k/m), m=1
  double x;
  double dxdt;
  double d2xdt2;
  double dt;
  double omega;

  sine_osc_nd(double freqHz, double sampleRate=DEFAULT_SRATE) :
    x(1.), dxdt(std::cos(x)), d2xdt2(-std::sin(dxdt)),
    k(std::pow(6.28318530718 * freqHz, 2.)),
    dt(1. / sampleRate),
    omega(6.28318530718 * freqHz)
  {
  }

  inline void SetFreq(double freqHz)
  {
    omega = 6.28318530718 * freqHz;
    k = omega * omega;
  }

  inline void SetSampleRate(double sampleRate)
  {
    dt = 1. / sampleRate;
  }

  inline double Step()
  {
    d2xdt2 = -k * (x + dxdt * dt);
    double dxdt_new = dxdt + d2xdt2 * dt;
    x += dxdt * dt;
    dxdt = dxdt_new;
    return x;
  }
};

END_DSP_NAMESPACE