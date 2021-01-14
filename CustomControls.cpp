#include "CustomControls.h"

template <int MAXNC>
SequencerControl<MAXNC>::SequencerControl(const IRECT& bounds, const char* label, const IVStyle& style, int nSteps, EDirection dir) :
  IVMultiSliderControl(bounds, label, style, nSteps, dir)
{

}