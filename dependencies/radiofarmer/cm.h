#pragma once

#include "radiofarmer_config.h"

BEGIN_DSP_NAMESPACE

#define CM_MAX_STAGES 100

struct ap_stage
{
  sample_t x1;
  sample_t y1;
};

template<int NumStages>
class ModulatedAllpass
{
public:
  ModulatedAllpass()
  {
    for (int i{ 0 }; i < NumStages; ++i)
    {
      mStages[i].x1 = 0.;
      mStages[i].y1 = 0.;
    }
  }

  inline sample_t Process(sample_t s, sample_t m)
  {
    ap_stage* stage = &mStages[0];

    sample_t y{};
#pragma clang loop unroll(full)
    for (int i{ 0 }; i < NumStages; ++i)
    {
      y = stage->x1 - m * (s - stage->y1);
      stage->x1 = s;
      stage->y1 = y;
      s = stage->y1;
      stage++;
    }
    return y;
  }

private:
  ap_stage mStages[CM_MAX_STAGES]{};
};

END_DSP_NAMESPACE