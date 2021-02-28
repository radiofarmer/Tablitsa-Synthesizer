#include "SignalProcessing.h"

template<>
double SoftClip(double s, double gain)
{
  double s_abs = std::tanh(std::abs(s * gain));
  double nonlin = (double)2. - (double)3. * s_abs * s_abs;
  return std::copysign((3. - nonlin * nonlin) / 4., s);
}
