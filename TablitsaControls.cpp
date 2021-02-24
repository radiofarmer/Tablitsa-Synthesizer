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