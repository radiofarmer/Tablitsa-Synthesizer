#include "waveshaper.h"

using namespace radiofarmer;

template<>
sample_t SoftClip<sample_t, 5>(sample_t x, sample_t threshold)
{
  sample_t x3 = x * x * x;
  sample_t x5 = x3 * x * x;
  sample_t pn = 0.375 * x5 - 1.25 * x3 + 1.875 * x;
  return std::clamp(pn, -threshold, threshold);
}

template<>
sample_v SoftClip<sample_v, 5>(sample_v x, sample_v threshold)
{
  sample_v x3 = pow(x, 3);
  sample_v x5 = pow(x3, 2);
  sample_v pn = 0.375 * x5 - 1.25 * x3 + 1.875 * x;
  return max(min(pn, threshold), -threshold);
}