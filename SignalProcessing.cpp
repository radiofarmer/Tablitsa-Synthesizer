#include "SignalProcessing.h"

template<typename T>
T SoftClip<T>(T s, T gain)
{
  /*T s_abs = std::abs(s);
  T s2 = (2 - 3 * std::tanh(s_abs * gain));
  T sum = (2 * s_abs * std::tanh(s_abs * gain) + (3 - s2 * s2) / 3 + std::tanh(s_abs * gain * gain)) / 3.;*/
  T s_abs = std::tanh(std::abs(s * gain));
  T nonlin = (T)2. - (T)3. * s_abs * s_abs;
  return std::copysign((3. - nonlin * nonlin) / 4., s);
}

template double SoftClip<double>(double s, double gain);

template float SoftClip<float>(float s, float gain);
