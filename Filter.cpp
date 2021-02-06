#include "Filter.h"
#include "IPlug_include_in_plug_src.h"

/* Comb Filter */

template class CombFilter<double>;

CombFilter<double>::CombFilter(double sampleRate, double feedforward, double feedback, bool cutoffIsNormalized, double delayLength) :
  Filter<T>(sampleRate, feedforward, feedback, cutoffIsNormalized), mDelayIn(delayLength), mDelayOut(delayLength)
{

}

inline void CombFilter<double>::SetDrive(double delayLength)
{
  int dl_samples = static_cast<int>(delayLength) * mMaxDelay;
  mDelayIn.SetDelay(dl_samples);
  mDelayOut.SetDelay(dl_samples);
}

inline double CombFilter<double>::Process(double s)
{
  double out = s + mFF * mDelayIn[0] - mFB * mDelayOut[0];
  mDelayIn.push(s);
  mDelayOut.push(out);
  return out;
}