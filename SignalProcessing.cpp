#include "SignalProcessing.h"
#include "vectormath_hyp.h"

template<>
extern inline double SoftClip<double>(double s, double gain)
{
  double s_abs = std::abs(s * gain);
  double nonlin = 2. - 3. * s_abs * s_abs;
  return std::copysign((4. - nonlin * nonlin) / 5., s);
}

template<>
extern inline float SoftClip<float>(float s, float gain)
{
  float s_abs = std::abs(s * gain);
  float nonlin = 2.f - 3.f * s_abs * s_abs;
  return std::copysignf((4.f - nonlin * nonlin) / 5.f, s);
}

template<>
extern inline Vec4d __vectorcall SoftClip<Vec4d, double>(const Vec4d& s, double gain)
{
  Vec4d s_abs = abs(s);
  Vec4d nonlin = 2. - 3. * pow(s_abs, 2);
  return sign(s) * (4. - pow(nonlin, 2) * 0.2);
}


