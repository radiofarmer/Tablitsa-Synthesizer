#include "biquad.h"

using namespace radiofarmer;

/* sets up a BiQuad Filter */
Biquad::Biquad(int type, sample_t dbGain, sample_t freq, sample_t bandwidth) : mType(type)
{
  CalculateCoefficients(dbGain, freq, bandwidth);
}

void Biquad::CalculateCoefficients(sample_t dbGain, sample_t freq, sample_t bandwidth)
{
  sample_t A, omega, sn, cs, alpha, beta;
  sample_t a0, a1, a2, b0, b1, b2;

  /* setup variables */
  A = pow(10, dbGain / 40);
  omega = 2 * M_PI * freq;
  sn = sin(omega);
  cs = cos(omega);
  alpha = sn * sinh(M_LN2 / 2 * bandwidth * omega / sn);
  beta = sqrt(A + A);

  switch (mType) {
  case LPF:
  {
    b0 = (1 - cs) / 2;
    b1 = 1 - cs;
    b2 = (1 - cs) / 2;
    a0 = 1 + alpha;
    a1 = -2 * cs;
    a2 = 1 - alpha;
    break;
  }
  case HPF:
  {
    b0 = (1 + cs) / 2;
    b1 = -(1 + cs);
    b2 = (1 + cs) / 2;
    a0 = 1 + alpha;
    a1 = -2 * cs;
    a2 = 1 - alpha;
    break;
  }
  case BPF:
  {
    b0 = alpha;
    b1 = 0;
    b2 = -alpha;
    a0 = 1 + alpha;
    a1 = -2 * cs;
    a2 = 1 - alpha;
    break;
  }
  case NOTCH:
  {
    b0 = 1;
    b1 = -2 * cs;
    b2 = 1;
    a0 = 1 + alpha;
    a1 = -2 * cs;
    a2 = 1 - alpha;
    break;
  }
  case PEQ:
  {
    b0 = 1 + (alpha * A);
    b1 = -2 * cs;
    b2 = 1 - (alpha * A);
    a0 = 1 + (alpha / A);
    a1 = -2 * cs;
    a2 = 1 - (alpha / A);
    break;
  }
  case LSH:
  {
    b0 = A * ((A + 1) - (A - 1) * cs + beta * sn);
    b1 = 2 * A * ((A - 1) - (A + 1) * cs);
    b2 = A * ((A + 1) - (A - 1) * cs - beta * sn);
    a0 = (A + 1) + (A - 1) * cs + beta * sn;
    a1 = -2 * ((A - 1) + (A + 1) * cs);
    a2 = (A + 1) + (A - 1) * cs - beta * sn;
    break;
  }
  case HSH:
  {
    b0 = A * ((A + 1) + (A - 1) * cs + beta * sn);
    b1 = -2 * A * ((A - 1) + (A + 1) * cs);
    b2 = A * ((A + 1) + (A - 1) * cs - beta * sn);
    a0 = (A + 1) - (A - 1) * cs + beta * sn;
    a1 = 2 * ((A - 1) - (A + 1) * cs);
    a2 = (A + 1) - (A - 1) * cs - beta * sn;
    break;
  }
  default:
  {
    b0 = b1 = b2 = a1 = a2 = 0.;
    a0 = 1.;
  }
  }

  /* precompute the coefficients */
  mBq.a0 = b0 / a0;
  mBq.a1 = b1 / a0;
  mBq.a2 = b2 / a0;
  mBq.a3 = a1 / a0;
  mBq.a4 = a2 / a0;

  /* zero initial samples */
  mBq.x1 = mBq.x2 = 0;
  mBq.y1 = mBq.y2 = 0;
}

/* Computes a BiQuad filter on a sample */
sample_t Biquad::Process(sample_t s)
{
  sample_t result;

  /* compute result */
  result = mBq.a0 * s + mBq.a1 * mBq.x1 + mBq.a2 * mBq.x2 -
    mBq.a3 * mBq.y1 - mBq.a4 * mBq.y2;

  /* shift x1 to x2, sample to x1 */
  mBq.x2 = mBq.x1;
  mBq.x1 = s;

  /* shift y1 to y2, result to y1 */
  mBq.y2 = mBq.y1;
  mBq.y1 = result;

  return result;
}