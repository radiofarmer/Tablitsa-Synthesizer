#include "IGraphics.h"
#include "Tablitsa.h"
#include "TablitsaControls.h"

using namespace igraphics;
using namespace iplug;


void TablitsaSliderControl::DrawTrack(IGraphics& g, const IRECT& filledArea)
{
  const float extra = mHandleInsideTrack ? mHandleSize : 0.f;
  const IRECT adjustedTrackBounds = mDirection == EDirection::Vertical ? mTrackBounds.GetVPadded(extra) : mTrackBounds.GetHPadded(extra);
  // Padd the filled area less, to account for the asymmetric rectangular handle
  const IRECT adjustedFillBounds = mDirection == EDirection::Vertical ? filledArea.GetVPadded(extra / 2.f) : filledArea.GetHPadded(extra / 2.f);
  const float cr = GetRoundedCornerRadius(mTrackBounds);

  g.FillRoundRect(GetColor(kSH), adjustedTrackBounds, cr, &mBlend);
  g.FillRoundRect(GetColor(kX1), adjustedFillBounds, cr, &mBlend);

  if (mStyle.drawFrame)
    g.DrawRoundRect(GetColor(kFR), adjustedTrackBounds, cr, &mBlend, mStyle.frameThickness);
}

ModSliderControl::ModSliderControl(const IRECT& bounds, int paramIdx, const char* label, const IVStyle& style, bool valueIsEditable, EDirection dir, double gearing, float handleSize, float trackSize, bool handleInsideTrack) :
  TablitsaSliderControl(bounds, paramIdx, label, style, valueIsEditable, dir, gearing, handleSize, trackSize, handleInsideTrack)
{
  mShape = EVShape::Rectangle;
  Toggle();
  SetActionFunction([this](IControl* pControl) {
    // Update control to whose parameter this slider is linked
    IControl* pTarget = GetUI()->GetControl(mTarget);
    if (pTarget)
      pTarget->SetDirty();
    });
}


/* Preset Selection Control */

PresetSelector::PresetSelector(const IRECT& bounds, IPopupMenuControl* menu, std::initializer_list<char*> defaultPresets) :
  ICaptionControl(bounds, kNoParameter, TABLITSA_TEXT, TABLITSA_STYLE.colorSpec.GetColor(EVColor::kBG)), mMenu(menu)
{
  for (auto p : defaultPresets)
  {
    mDefaultPresets.push_back(std::string(p));
  }
  // LoadUserPresets(...);
  mAllPresets.insert(mAllPresets.begin(), mDefaultPresets.begin(), mDefaultPresets.end());
  mAllPresets.insert(mAllPresets.end(), mUserPresets.begin(), mUserPresets.end());
}

void PresetSelector::LoadUserPresets(std::initializer_list<char*> userPresets)
{
  for (auto p : userPresets)
  {
    mUserPresets.push_back(std::string(p));
  }
}

void PresetSelector::Draw(IGraphics& g)
{
  const IParam* pParam = GetParam();

  mStr.Set(mDefaultPresets[0].c_str());

  ITextControl::Draw(g);

  if (mTri.W() > 0.f)
  {
    g.FillTriangle(mMouseIsOver ? mTriangleMouseOverColor : mTriangleColor, mTri.L, mTri.T, mTri.R, mTri.T, mTri.MW(), mTri.B, GetMouseIsOver() ? 0 : &BLEND_50);
  }
}