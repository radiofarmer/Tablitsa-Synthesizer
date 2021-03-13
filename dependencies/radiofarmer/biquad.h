#pragma once
/**
* Adaptation of C code by Tom St. Denis, following the Biquad cookbook formula by Robert Bristow-Johson.
*
* See radiofarmer/third-party/biquad.c for original source code.
*/

#include "radiofarmer_config.h"

#include <math.h>
#include <stdlib.h>

BEGIN_DSP_NAMESPACE

#ifndef M_LN2
#define M_LN2	   0.69314718055994530942
#endif

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

/* filter types */
enum EBiquadTypes {
  LPF, /* low pass filter */
  HPF, /* High pass filter */
  BPF, /* band pass filter */
  NOTCH, /* Notch Filter */
  PEQ, /* Peaking band EQ filter */
  LSH, /* Low shelf filter */
  HSH /* High shelf filter */
};

/* this holds the data required to update samples thru a filter */
struct biquad
{
  sample_t a0, a1, a2, a3, a4;
  sample_t x1, x2, y1, y2;
};

class Biquad
{
public:
  Biquad(int type=LPF, sample_t dbGain=0., sample_t freq=0.5, sample_t bandwidth=2.);

  void CalculateCoefficients(sample_t dbGain, sample_t freq, sample_t bandwidth);

  sample_t Process(sample_t s);

protected:
  const int mType;
  biquad mBq;
};



END_DSP_NAMESPACE