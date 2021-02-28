#include "SignalProcessing.h"

template<>
double SoftClip(double s, double gain)
{
  double s_abs = std::tanh(std::abs(s * gain));
  double nonlin = 2. - 3. * s_abs * s_abs;
  return std::copysign((4. - nonlin * nonlin) / 5., s);
}
