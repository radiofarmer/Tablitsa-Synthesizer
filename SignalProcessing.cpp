#include "SignalProcessing.h"

double SoftClip(double s, double gain)
{
  /*double s_abs = std::abs(s);
  double s2 = (2 - 3 * std::tanh(s_abs * gain));
  double sum = (2 * s_abs * std::tanh(s_abs * gain) + (3 - s2 * s2) / 3 + std::tanh(s_abs * gain * gain)) / 3.;*/
  double s_abs = std::tanh(std::abs(s * gain));
  double nonlin = 2. - 3 * s_abs * s_abs;
  return std::copysign((3. - nonlin * nonlin) / 4., s);
}
