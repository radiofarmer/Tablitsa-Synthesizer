#include "TablitsaControls.h"

IVModKnobControl::IVModKnobControl(const IRECT& bounds, int paramIdx, const char* label, const IVStyle& style, bool valueIsEditable, double gearing)
  : IVModKnobControl(bounds, paramIdx, paramIdx + 1, label, style, valueIsEditable, gearing)
{
}

IVModKnobControl::IVModKnobControl(const IRECT& bounds, int paramIdx, int modStartIdx, const char* label, const IVStyle& style, bool valueIsEditable, double gearing)
  : IVKnobControl(bounds, paramIdx, label, style, valueIsEditable), mDefaultColor{ GetColor(kFG) }, mGearing(gearing), mModParamIdx(modStartIdx)
{
  GetMouseDblAsSingleClick();
}

void IVModKnobControl::OnMouseDown(float x, float y, const IMouseMod& mod)
{

  if (!mod.L)
  {
    LoadModParams();
    mMouseDown = !mMouseDown;
    /* By default, center-clicking causes the control to be captured such that it still responds to the mouse wheel when
    the mouse is not actually over it. ReleaseMouseCapture() empties the captured-control queue. */
    GetUI()->ReleaseMouseCapture();
  }
  else
    IVKnobControl::OnMouseDown(x, y, mod);
}