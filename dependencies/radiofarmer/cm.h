#pragma once

#include "radiofarmer_config.h"

BEGIN_DSP_NAMESPACE

#define CM_MAX_STAGES 100

struct ap_stage
{
  sample_t x1;
  sample_t y1;
};

class AllpassCascade
{
public:
  AllpassCascade(int numStages);

  sample_t Process(sample_t s, sample_t m);

private:
  const int mNumStages;
  ap_stage mStages[CM_MAX_STAGES]{};
};

END_DSP_NAMESPACE