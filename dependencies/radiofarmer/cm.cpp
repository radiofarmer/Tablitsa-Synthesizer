#include "cm.h"

BEGIN_DSP_NAMESPACE

AllpassCascade::AllpassCascade(int numStages) : mNumStages(numStages)
{
  for (int i{ 0 }; i < numStages; ++i)
  {
    mStages[i].x1 = 0.;
    mStages[i].y1 = 0.;
  }
}

sample_t AllpassCascade::Process(sample_t s, sample_t m)
{
  const int n = mNumStages;
  ap_stage* stage = &mStages[0];

  sample_t y{};
  for (int i{ 0 }; i < n; ++i)
  {
    y = stage->x1 - m * (s - stage->y1);
    stage->x1 = s;
    stage->y1 = y;
    s = stage->y1;
    stage++;
  }
  return y;
}

END_DSP_NAMESPACE