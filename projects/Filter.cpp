#include "Filter.h"

template<typename T>
inline void SVF2<T>::SetCutoff(double cutoffNorm)
{
  mFc = std::max(std::min(cutoffNorm, 0.99), 0.001);
}