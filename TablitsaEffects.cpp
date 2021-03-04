#include "TablitsaEffects.h"
#include "VectorFunctions.h"

using namespace iplug;

template<>
Vec4d Effect<Vec4d>::Process(Vec4d s)
{
  return s;
}

template<>
sample Effect<sample>::Process(sample s)
{
  return s;
}

template<>
void Effect<sample>::ProcessStereo(sample* s) {}