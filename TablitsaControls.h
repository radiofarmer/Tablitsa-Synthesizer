#pragma once
#include "PeriodicTable.h"
#include "Modulation.h"

#include "IGraphics.h"
#include "IControls.h"

BEGIN_IPLUG_NAMESPACE
BEGIN_IGRAPHICS_NAMESPACE

const IText TABLITSA_TEXT = IText().WithFGColor(COLOR_WHITE);
const IVStyle TABLITSA_STYLE = IVStyle(DEFAULT_SHOW_LABEL,
  DEFAULT_SHOW_VALUE,
  /* Background       Foreground      Pressed                   Frame            Highlight    Shadow           Extra 1          Extra 2          Extra 3        */
  { DEFAULT_BGCOLOR, COLOR_DARK_GRAY, IColor(255, 225, 0, 190), DEFAULT_FRCOLOR, COLOR_WHITE, DEFAULT_SHCOLOR, DEFAULT_X1COLOR, DEFAULT_X2COLOR, DEFAULT_X3COLOR },
  DEFAULT_LABEL_TEXT.WithFGColor(COLOR_WHITE),
  DEFAULT_VALUE_TEXT.WithFGColor(COLOR_WHITE),
  DEFAULT_HIDE_CURSOR,
  DEFAULT_DRAW_FRAME,
  false, // Draw shadows
  DEFAULT_EMBOSS,
  DEFAULT_ROUNDNESS,
  DEFAULT_FRAME_THICKNESS,
  DEFAULT_SHADOW_OFFSET,
  DEFAULT_WIDGET_FRAC,
  DEFAULT_WIDGET_ANGLE);


const IVColorSpec knobColorSpec = IVColorSpec{
  COLOR_TRANSPARENT,
  COLOR_BLACK,
  COLOR_WHITE.WithOpacity(0.75),
  COLOR_WHITE.WithOpacity(0.5),
  DEFAULT_HLCOLOR,
  DEFAULT_SHCOLOR,
  DEFAULT_X1COLOR,
  DEFAULT_X2COLOR,
  DEFAULT_X3COLOR
};

const IVStyle modKnobStyle{ TABLITSA_STYLE.WithColors(knobColorSpec).WithLabelText(TABLITSA_STYLE.labelText.WithSize(17.f)) };

const IVStyle toggleStyle{ TABLITSA_STYLE.WithDrawFrame(true).WithColor(EVColor::kFG, COLOR_TRANSPARENT).WithShowLabel(false) };

const IText dropdownText{ DEFAULT_TEXT.WithFGColor(COLOR_WHITE) };

const IVStyle pushButtonStyle{ TABLITSA_STYLE.WithLabelText(TABLITSA_STYLE.labelText.WithVAlign(EVAlign::Middle)).WithColor(EVColor::kHL, TABLITSA_STYLE.colorSpec.GetColor(EVColor::kPR).WithContrast(0.5)) };

class TablitsaSliderControl : public IVSliderControl
{
public:
  TablitsaSliderControl(const IRECT& bounds, int paramIdx = kNoParameter, const char* label = "", const IVStyle& style = TABLITSA_STYLE.WithShowValue(true), bool valueIsEditable = false, EDirection dir = EDirection::Vertical, double gearing = DEFAULT_GEARING, float handleSize = 8.f, float trackSize = 2.f, bool handleInsideTrack = true) :
    IVSliderControl(bounds, paramIdx, label, style, valueIsEditable, dir, gearing, handleSize, trackSize, handleInsideTrack)
  {
    mShape = EVShape::Rectangle;
  }

  void OnAttached() override
  {
    SetColor(EVColor::kX1, COLOR_WHITE);
  }

  void DrawTrack(IGraphics& g, const IRECT& filledArea) override;

  void DrawHandle(IGraphics& g, const IRECT& bounds) override
  {
    const IRECT boundsAdj = mDirection == EDirection::Vertical ? bounds.GetMidVPadded(bounds.H() / 4.f) : bounds.GetMidHPadded(bounds.W() / 4.f);
    DrawPressableShape(g, mShape, boundsAdj, mMouseDown, mMouseIsOver, IsDisabled());
  }
};

class ModSliderControl : public TablitsaSliderControl
{
public:
  ModSliderControl(const IRECT& bounds, int paramIdx = kNoParameter, const char* label = "", const IVStyle& style = TABLITSA_STYLE.WithShowValue(true), bool valueIsEditable = false, EDirection dir = EDirection::Vertical, double gearing = DEFAULT_GEARING, float handleSize = 8.f, float trackSize = 2.f, bool handleInsideTrack = true);

  void Toggle(int targetParam=-1)
  {
    if (GetParamIdx() == kNoParameter)
    {
      SetValueStr("N/A");
      SetDisabled(true);
    }
    else
      SetDisabled(false);
    mTarget = targetParam;
  }

  // Overridden drawing function to draw the filled area starting in the middle of the track
  void DrawWidget(IGraphics& g) override
  {
    float value = (float)GetValue(); // NB: Value is normalized to between 0. and 1.
    const IRECT handleBounds = (GetParamIdx() == kNoParameter) ? mTrackBounds.FracRect(mDirection, 0.f) : mTrackBounds.FracRect(mDirection, value);
    const IRECT filledTrack = mDirection == EDirection::Vertical ? (GetParamIdx() == kNoParameter) ?
        handleBounds : ((value >= 0.5f) ?
        mTrackBounds.GetGridCell(0, 0, 2, 1).FracRect(mDirection, 2.f * (value - 0.5f)) :
        mTrackBounds.GetGridCell(1, 0, 2, 1).FracRect(mDirection, 2.f * (0.5 - value), true)) : // <- Vertical
      (GetParamIdx() == kNoParameter) ?
        handleBounds : ((value >= 0.5f) ?
        mTrackBounds.GetGridCell(0, 1, 1, 2).FracRect(mDirection, 2.f * (value - 0.5f)) :
        mTrackBounds.GetGridCell(0, 0, 1, 2).FracRect(mDirection, 2.f * (0.5 - value), true)); // <- Horizontal

    if (mTrackSize > 0.f)
      DrawTrack(g, filledTrack);

    float cx, cy;

    const float offset = (mStyle.drawShadows && mShape != EVShape::Ellipse /* TODO? */) ? mStyle.shadowOffset * 0.5f : 0.f;

    if (mDirection == EDirection::Vertical)
    {
      cx = handleBounds.MW() + offset;
      cy = handleBounds.T;
    }
    else
    {
      cx = handleBounds.R;
      cy = handleBounds.MH() + offset;
    }

    if (mHandleSize > 0.f)
    {
      DrawHandle(g, { cx - mHandleSize, cy - mHandleSize, cx + mHandleSize, cy + mHandleSize });
    }
  }

  void DrawTrack(IGraphics& g, const IRECT& filledArea) override
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

private:
  int mTarget{ -1 }; // The parameter currently linked to the modulation controls
};

class TablitsaIVKnobControl : public IVKnobControl
{
public:
  /* Create a knob control with a modulatable value */
  TablitsaIVKnobControl(const IRECT& bounds, int paramIdx, const char* label = "", const IVStyle& style = modKnobStyle, bool valueIsEditable = false, double gearing = DEFAULT_GEARING)
    : IVKnobControl(bounds, paramIdx, label, style, valueIsEditable, false, -135.f, 135.f, -135.f, EDirection::Vertical, gearing, 2.f)
  {
  }

  // Spin-Up Electron Pointer
  void DrawPointer(IGraphics& g, float angle, float cx, float cy, float radius)
  {
    const IColor pointerColor = GetColor(kFR).WithOpacity(1.f);
    g.DrawRadialLine(pointerColor, cx, cy, angle, mInnerPointerFrac * radius, mOuterPointerFrac * radius, &mBlend, mPointerThickness);
    float data1[2][2];
    float data2[2][2];
    iplug::igraphics::RadialPoints(angle, cx, cy, mInnerPointerFrac * radius, mOuterPointerFrac * radius, 2, data1);
    iplug::igraphics::RadialPoints(angle - 20., cx, cy, mInnerPointerFrac * radius, mOuterPointerFrac * radius * 0.65f, 2, data2);
    g.DrawLine(pointerColor, data1[1][0], data1[1][1], data2[1][0], data2[1][1], &mBlend, mPointerThickness);
  }
};

/* A clone of the normal knob control, but with the ability to receive modulation parameters. */
class TablitsaIVModKnobControl : public TablitsaIVKnobControl
{
public:
  /* Create a knob control with a modulatable value */

  TablitsaIVModKnobControl(const IRECT& bounds, int paramIdx, const char* label = "", const IVStyle& style = modKnobStyle, bool valueIsEditable = false, double gearing = DEFAULT_GEARING)
    : TablitsaIVModKnobControl(bounds, paramIdx, paramIdx + 1, label, style, valueIsEditable, gearing)
  {
  }

  TablitsaIVModKnobControl(const IRECT& bounds, int paramIdx, int modStartIdx, const char* label = "", const IVStyle& style = TABLITSA_STYLE, bool valueIsEditable = false, double gearing = DEFAULT_GEARING)
    : TablitsaIVKnobControl(bounds, paramIdx, label, style, valueIsEditable, gearing),
    mDefaultColor{ GetColor(kFG) },
    mModParamIdx(modStartIdx),
    mBaseTrack(style.colorSpec.GetColor(EVColor::kFR))
  {
    GetMouseDblAsSingleClick();
  }

  /* Get modulator values from a different parameter. (For mutually-exclusive parameters) */
  void GetModulationFrom(int paramIdx)
  {
    mModParamIdx = paramIdx + 1;
  }

  void ColorSwap()
  {
    if (mActive && GetActiveIdx() == GetUI()->GetControlIdx(this))
      SetColor(kFG, GetColor(kPR));
    else
      SetColor(kFG, mDefaultColor);
  }

  void OnMouseDown(float x, float y, const IMouseMod& mod) override
  {
    if (!mod.L || (mod.L && mod.A))
    {
      LoadModParams();
      mMouseDown = !mMouseDown;
      /* By default, center-clicking causes the control to be captured such that it still responds to the mouse wheel when
      the mouse is not actually over it. ReleaseMouseCapture() empties the captured-control queue. */
      GetUI()->ReleaseMouseCapture();
      GetUI()->SetAllControlsDirty();
    }
    else
      IVKnobControl::OnMouseDown(x, y, mod);
  }
  void OnMouseOut() override
  {
    IVKnobControl::OnMouseOut();
  }

  void OnMouseWheel(float x, float y, const IMouseMod& mod, float d) override
  {
    if (mMouseIsOver)
      IVKnobControl::OnMouseWheel(x, y, mod, d);
  }

  void LoadModParams()
  {
    // TODO make a list of Control Tags that can be looped through
    if (mActive && GetActiveIdx() == GetUI()->GetControlIdx(this))
    {
      for (int i{ kCtrlTagEnv1Depth }; i <= kCtrlTagRndDepth; ++i)
      {
        GetUI()->GetControlWithTag(i)->SetParamIdx(kNoParameter);
        dynamic_cast<ModSliderControl*>(GetUI()->GetControlWithTag(i))->Toggle(-1);
      }
      SetActiveIdx(false);
      mActive = false;
    }
    else
    {
      // Set all modulator sliders to the values of the currently-selected modulated parameter
      for (int i{ 0 }; i <= (kCtrlTagRndDepth - kCtrlTagEnv1Depth); ++i)
      {
        GetUI()->GetControlWithTag(kCtrlTagEnv1Depth + i)->SetParamIdx(mModParamIdx + i);
        dynamic_cast<ModSliderControl*>(GetUI()->GetControlWithTag(kCtrlTagEnv1Depth + i))->Toggle(GetUI()->GetControlIdx(this));
      }
      SetActiveIdx(true);
      mActive = true;
    }
    // Send values and change this control's active state
    GetDelegate()->SendCurrentParamValuesFromDelegate();
    //mActive = !mActive;
  }

  void DrawWidget(IGraphics& g) override
  {
    float widgetRadius; // The radius out to the indicator track arc

    if (mWidgetBounds.W() > mWidgetBounds.H())
      widgetRadius = (mWidgetBounds.H() / 2.f);
    else
      widgetRadius = (mWidgetBounds.W() / 2.f);

    const float cx = mWidgetBounds.MW(), cy = mWidgetBounds.MH();

    widgetRadius -= (mTrackSize / 2.f);

    IRECT knobHandleBounds = mWidgetBounds.GetCentredInside((widgetRadius - mTrackToHandleDistance) * 2.f);
    const float angle = mAngle1 + (static_cast<float>(GetValue()) * (mAngle2 - mAngle1));
    DrawPressableShape(g, EVShape::Ellipse, knobHandleBounds, mActive && (GetActiveIdx() == GetUI()->GetControlIdx(this)), mMouseIsOver, IsDisabled());
    DrawIndicatorTrack(g, angle, cx, cy, widgetRadius);
    DrawPointer(g, angle, cx, cy, knobHandleBounds.W() / 2.f);
  }

  void Draw(IGraphics& g) override
  {
    DrawBackground(g, mRECT);
    DrawLabel(g);
    DrawWidget(g);
    DrawValue(g, mValueMouseOver);
  }

  void DrawIndicatorTrack(IGraphics& g, float angle, float cx, float cy, float radius) override
  {
    // Set the origin of the track arch to the center of the range for dials with negative minimum values
    if (GetParam()->GetMin() < 0.)
      mAnchorAngle = 0.;
    if (mTrackSize > 0.f)
    {
      g.DrawArc(mBaseTrack, cx, cy, radius, angle >= mAnchorAngle ? mAnchorAngle : mAnchorAngle - (mAnchorAngle - angle), angle >= mAnchorAngle ? angle : mAnchorAngle, &mBlend, mTrackSize);

      // Envelope 1
      float env1Val = static_cast<float>(GetDelegate()->GetParam(mModParamIdx)->Value());
      if (env1Val)
      {
        float modAngle = std::max(std::min(env1Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[0][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 3);
      }

      // Envelope 2
      radius -= mTrackToHandleDistance / 2.f;
      float env2Val = static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 1)->Value());
      if (env2Val)
      {
        float modAngle = std::max(std::min(env2Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[1][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 2);
      }

      // AmpEnv
      radius -= mTrackToHandleDistance / 4.f;
      float ampEnvVal = static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 2)->Value());
      if (ampEnvVal)
      {
        float modAngle = std::max(std::min(ampEnvVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[2][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
      }

      // LFO 1
      radius -= 3.f * mTrackToHandleDistance / 4.f;
      float lfo1Val = static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 3)->Value());
      if (lfo1Val)
      {
        float modAngle = std::max(std::min(lfo1Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[3][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 3);
        modAngle = std::max(std::min(-1.f * lfo1Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[3][1], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 2);
      }

      // LFO 2
      radius += mTrackToHandleDistance / 4.f;
      float lfo2Val = static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 4)->Value());
      if (lfo2Val)
      {
        float modAngle = std::max(std::min(lfo2Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[4][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 2);
        modAngle = std::max(std::min(-1.f * lfo2Val * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[4][1], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
      }

      // Sequencer
      radius += mTrackToHandleDistance / 4.f;
      float seqVal = static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 5)->Value());
      if (seqVal)
      {
        float modAngle = std::max(std::min(seqVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[5][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
      }

      // Velocity
      radius -= mTrackToHandleDistance / 4.f;
      float velVal = static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 6)->Value());
      if (velVal)
      {
        float modAngle = std::max(std::min(velVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[6][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
      }

      // Keytrack
      float ktrVal = static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 7)->Value());
      if (ktrVal)
      {
        float modAngle = std::max(std::min(ktrVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[7][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
      }

      // Trigger Random
      float rndVal = static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 8)->Value());
      if (rndVal)
      {
        float modAngle = std::max(std::min(rndVal * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
        g.DrawArc(mModArcColor[8][0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
      }
    }
  }

  inline int GetActiveIdx()
  {
    return dynamic_cast<Tablitsa*>(GetDelegate())->GetActiveModIdx();
  }

  inline void SetActiveIdx(int idx)
  {
    return dynamic_cast<Tablitsa*>(GetDelegate())->SetActiveModIdx(idx);
  }

  inline void SetActiveIdx(bool isActive)
  {
    return dynamic_cast<Tablitsa*>(GetDelegate())->SetActiveModIdx(isActive ? GetUI()->GetControlIdx(this) : -1);
  }

protected:
  std::vector<ISVG*> mSVG;
  const IColor mDefaultColor;
  const IColor mModArcColor[9][2]{
    {{200, 200, 100, 100}, {}},
    {{200, 225, 150, 100}, {}},
    {{200, 250, 200, 100}, {}},
    {{200, 160, 0, 225}, {200, 200, 0, 235}},
    {{200, 225, 0, 190}, {200, 235, 0, 220}},
    {{200, 0, 250, 100}, {}},
    {{150, 0, 100, 255}, {}},
    {{150, 200, 0, 200}, {}},
    {{150, 0, 255, 255}, {}}
  };
  const IColor mWhite{ 200, 255, 255, 255 };
  const IColor mBaseTrack;
  bool mActive{ false };
  int mModParamIdx;
};


class PeriodicTable : public IControl
{
  int ElementCoords[118][2]{ ELEMENT_COORDS };
  int LaAcCoords[30][2]{ LA_AC_COORDS };
  char* ElementNames[118]{ ELEMENT_NAMES_HYPHENATED };
  char* ElementSymbols[118]{ ELEMENT_SYMBOLS };
  static constexpr float TableWidth{ 753.f };
  static constexpr float TableHeight{ 290.f };
  static constexpr float LaTLHC[2]{120.33f, 351.62f };
  static constexpr float LaAcWidth{ 626.26f };
  static constexpr float LaAcHeight{ 79.86f };

  IColor ElementIconColor[2]{ { 255, 200, 200, 0 },
    {255, 0, 100, 225} };
  IColor LabelColor{ 255, 200, 0, 150 };
  IText ElementIconText{ 16, EAlign::Center, ElementIconColor[0] };
  IText LabelText{16, EAlign::Center, LabelColor};

public:
  PeriodicTable(const IRECT& bounds, const ISVG& svg, int paramIdx) : PeriodicTable(bounds, svg)
  {
    SetParamIdx(paramIdx);
  }

  PeriodicTable(const IRECT& bounds, const ISVG& svg, const std::initializer_list<int> params) : PeriodicTable(bounds, svg)
  {
    SetNVals(static_cast<int>(params.size()));
    int valIdx = 0;
    for (auto param : params)
    {
      SetParamIdx(param, valIdx++);
    }
  }

  PeriodicTable(const IRECT& bounds, const ISVG& svg) : IControl(bounds), mSVG(svg)
  {
#define TABLE_ADJ 1.5;
    // The amount by which the SVG is scaled in order to fit in the control
    mScaleFact = std::min((mRECT.W() - 80.f) / mSVG.W(), mRECT.H() / mSVG.H());
    // Table height in plugin (without lanthanides and actinides)
    mTableWidth = std::ceil(TableWidth * mScaleFact) + TABLE_ADJ;
    mTableHeight = std::ceil(TableHeight * mScaleFact) + TABLE_ADJ;
    mTableTLHC[0] = mRECT.L + 37.f * mScaleFact;
    mTableTLHC[1] = mRECT.T + 37.f * mScaleFact;
    // Lanthanides and Actinides
    mLaAcWidth = std::ceil(LaAcWidth * mScaleFact) + TABLE_ADJ;
    mLaAcHeight = std::ceil(LaAcHeight * mScaleFact) + TABLE_ADJ;
    mLaTLHC[0] = mRECT.L + LaTLHC[0] * mScaleFact;
    mLaTLHC[1] = mRECT.T + LaTLHC[1] * mScaleFact;
    mTableBounds = IRECT(mTableTLHC[0], mTableTLHC[1], mTableTLHC[0] + mTableWidth, mTableTLHC[1] + mTableHeight);
    mLaAcBounds = IRECT(mLaTLHC[0], mLaTLHC[1], mLaTLHC[0] + mLaAcWidth, mLaTLHC[1] + mLaAcHeight);
  }

  void LoadValues()
  {
    // TODO: Make the number of elements a template parameter to make it easier to add an oscillator in the future
    for (auto i{0}; i < 2; ++i)
      mSelectedElements[i] = GetParam(i)->Value();
    SetDirty(false);
  }

  void SetSelectedElement(int atomicNumber, int elemIdx)
  {
    assert(atomicNumber >= 1 && atomicNumber <= 118);
    mSelectedElements[elemIdx] = atomicNumber;
  }
  
  void Draw(IGraphics& g)
  {
    //g.FillRoundRect(IColor(255, 20, 0, 45), mRECT);
    // Highlight box occupied by mouse
    if (mCurrentElementCoords)
    {
      if (!((mCurrentElement >= 57 && mCurrentElement <= 71) || (mCurrentElement >= 89 && mCurrentElement <= 103)))
        g.FillRect(IColor(150, 200, 200, 0), mTableBounds.GetReducedFromLeft(std::floor(*(mCurrentElementCoords + 1) / 18. * mTableWidth)) \
          .GetReducedFromTop(std::floor(*mCurrentElementCoords / 7.f * mTableHeight)) \
          .GetFromLeft(E_CELL_WIDTH) \
          .GetFromTop(E_CELL_HEIGHT));
      else
        g.FillRect(IColor(150, 200, 200, 0), mLaAcBounds.GetReducedFromLeft(std::floor(*(mCurrentElementCoords + 1) / 15. * mLaAcWidth)) \
          .GetReducedFromTop(std::floor(*mCurrentElementCoords / 2.f * mLaAcHeight)) \
          .GetFromLeft(E_CELL_WIDTH) \
          .GetFromTop(E_CELL_HEIGHT));
    };
    for (int i{ 0 }; i < 2; ++i)
    {
      auto e = mSelectedElements[i];
      if (!((e >= 57 && e <= 71) || (e >= 89 && e <= 103)))
      {
        int* eCoords = ElementCoords[e - 1];
        g.FillRect(ElementIconColor[i].WithOpacity(0.8), mTableBounds.GetReducedFromLeft(std::floor(*(eCoords + 1) / 18. * mTableWidth)) \
          .GetReducedFromTop(std::floor(*eCoords / 7.f * mTableHeight)) \
          .GetFromLeft(E_CELL_WIDTH) \
          .GetFromTop(E_CELL_HEIGHT));
      }
      else
      {
        // Convert atomic number to index in Lanthanide-Actinide table coordinates. Note that Number{Ac} - 74 corresponds to index 15, the first column in the second row
        int* eCoords = (e >= 89) ? LaAcCoords[e - 74] : LaAcCoords[e - 57];
        g.FillRect(ElementIconColor[i].WithOpacity(0.8), mLaAcBounds.GetReducedFromLeft(std::floor(*(eCoords + 1) / 15. * mLaAcWidth)) \
          .GetReducedFromTop(std::floor(*eCoords / 2.f * mLaAcHeight)) \
          .GetFromLeft(E_CELL_WIDTH) \
          .GetFromTop(E_CELL_HEIGHT));
      }
    };

    // Draw current wavetable icons
    const IRECT activeElem1 = mRECT.GetFromRight(80.f).GetGridCell(0, 0, 2, 1).GetCentredInside(70.f, 80.f);
    const IRECT activeElem2 = mRECT.GetFromRight(80.f).GetGridCell(1, 0, 2, 1).GetCentredInside(70.f, 80.f);
    DrawElement(g, activeElem1, mSelectedElements[0], 0);
    DrawElement(g, activeElem2, mSelectedElements[1], 1);

    g.DrawSVG(mSVG, mRECT);
  }

  void DrawElement(IGraphics& g, const IRECT& bounds, int atomicNumber, int idx)
  {
    assert(atomicNumber > 0);
    IColor col{ mTableLoading[idx] ? ElementIconColor[idx].WithOpacity(0.5f) : ElementIconColor[idx] };

    g.DrawRect(col, bounds);
    g.DrawText(LabelText, (atomicNumber) ? "Wavetable 2" : "Wavetable 1", bounds.GetVShifted(-20.f).GetFromTop(20.f));
    g.DrawText(ElementIconText.WithSize(42).WithFGColor(col), ElementSymbols[atomicNumber - 1], bounds.GetReducedFromTop(25.f).GetReducedFromBottom(30.f));

    // Check for multiline name and draw name
    std::string eName = ElementNames[atomicNumber - 1];
    const size_t breakpoint{ eName.find("\n") };
    if (breakpoint != std::string::npos)
    {
      g.DrawText(ElementIconText.WithSize(16).WithFGColor(col), (" " + eName.substr(0, breakpoint)).c_str(), bounds.GetFromBottom(25.f).GetFromTop(12.f));
      g.DrawText(ElementIconText.WithSize(16).WithFGColor(col), eName.substr(breakpoint).c_str(), bounds.GetFromBottom(12.f).GetFromTop(12.f));
    }
    else
      g.DrawText(ElementIconText.WithSize(16).WithFGColor(col), eName.c_str(), bounds.GetFromBottom(25.f).GetFromTop(12.f));
    // Draw atomic number
    g.DrawText(ElementIconText.WithSize(24).WithFGColor(col), std::to_string(atomicNumber).c_str(), bounds.GetFromTop(24.f));
  }
  
  void OnMouseOver(float x, float y, const IMouseMod& mod)
  {
    if (x < mTableTLHC[0] || y < mTableTLHC[1])
    {
      mCurrentElementCoords = nullptr;
      mCurrentElement = -1;
      goto EndFunction;
    }
    float x_adj{ (x - mTableTLHC[0]) / mTableWidth };
    float y_adj{ (y - mTableTLHC[1]) / mTableHeight };
    int atomicNumber{ 1 };
    for (auto e : ElementCoords)
    {
      // Check main table (excluding lanthanides and actinides)
      if (e[0] < 7.)
      {
        if (static_cast<int>(y_adj * 7.f) == e[0] && static_cast<int>(x_adj * 18.f) == e[1])
        {
          mCurrentElementCoords = e;
          mCurrentElement = atomicNumber;
          goto EndFunction;
        }
      }
      // Increment the current atomic number and reset the Current Element index and pointer
      atomicNumber++;
      mCurrentElementCoords = nullptr;
      mCurrentElement = -1;
    }
    x_adj = (x - mLaTLHC[0]) / mLaAcWidth;
    y_adj = (y - mLaTLHC[1]) / mLaAcHeight;
    atomicNumber = 57;
    if (x < mLaTLHC[0] || y < mLaTLHC[1])
    {
      mCurrentElementCoords = nullptr;
      mCurrentElement = -1;
      goto EndFunction;
    }
    // Check the lanthanide and actinide table
    for (auto e : LaAcCoords)
    {
      if (static_cast<int>(y_adj * 2.f) == e[0] && static_cast<int>(x_adj * 15.f) == e[1])
      {
        mCurrentElementCoords = e;
        mCurrentElement = atomicNumber;
        break;
      }
      if (++atomicNumber > 71 && atomicNumber < 89)
        atomicNumber = 89;
      mCurrentElementCoords = nullptr;
      mCurrentElement = -1;
    }
  EndFunction:
    if (mCurrentElement > -1 && mIsDragging > -1)
    {
      mSelectedElements[mIsDragging] = mCurrentElement;
    }
    SetDirty(false);
  }
  
  void OnMouseOut()
  {
  }

  void OnMouseDrag(float x, float y, float dx, float dy, const IMouseMod& mod)
  {
    OnMouseOver(x, y, mod);
  }

  void OnMouseDown(float x, float y, const IMouseMod& mod)
  {
    for (int i{ 0 }; i < 2; ++i)
    {
      if (mCurrentElement == mSelectedElements[i])
      {
        mIsDragging = i;
      }
    }
  }

  void OnMouseUp(float x, float y, const IMouseMod& mod)
  {
    if (mIsDragging > -1)
    {
      // Send new wavetable value (param method)
      SetValue(static_cast<double>(mSelectedElements[mIsDragging]) / 118., mIsDragging);

      // Send new wavetable value (hidden param method)
      double newTableIdx = static_cast<double>(mSelectedElements[mIsDragging]) / 118.;
      GetUI()->GetDelegate()->SendArbitraryMsgFromUI(mIsDragging == 0 ? kMsgWavetable1Changed : kMsgWavetable2Changed, kCtrlTagPeriodicTable, sizeof(double), &newTableIdx);
//      GetDelegate()->SendParameterValueFromUI(GetParamIdx(mIsDragging), static_cast<double>(mSelectedElements[mIsDragging] - 1) / 118.);
//      GetDelegate()->OnParamChange(GetParamIdx(mIsDragging));
    };
    mIsDragging = -1;
    SetDirty(false);
  }
  
  void SetTableLoading(const bool isLoading, const int tableIdx)
  {
    mTableLoading[tableIdx] = isLoading;
  }

protected:
  ISVG mSVG;
  int* mCurrentElementCoords{ nullptr };
  int mCurrentElement{ -1 };
  int mSelectedElements[2]{ 1, 2 };
  float mScaleFact{ 1. };
  float mTableTLHC[2]{}; // Periodic table top left-hand corner
  float mLaTLHC[2]{}; // Lanthanum TLHC
  float mTableWidth{};
  float mTableHeight{};
  float mLaAcWidth{};
  float mLaAcHeight{};
  IRECT mTableBounds;
  IRECT mLaAcBounds;
  int mIsDragging{ -1 };
  bool mTableLoading[2]{ false, false };
};

/* Custom Toggle Button */

class TablitsaVToggleButton : public IVToggleControl
{
public:
  TablitsaVToggleButton(const IRECT& bounds, int paramIdx, const char* text = "", const IVStyle& style = TABLITSA_STYLE.WithDrawShadows(false)) :
    TablitsaVToggleButton(bounds, paramIdx, "", style, text, text) {}

  TablitsaVToggleButton(const IRECT& bounds, int paramIdx = kNoParameter, const char* label = "", const IVStyle& style = TABLITSA_STYLE.WithDrawShadows(false), const char* offText = "OFF", const char* onText = "ON") :
    IVToggleControl(bounds, paramIdx, label, style, offText, onText)
  {

  }

protected:
  IBitmap mBitmap;
};

#define GROUP_BOX_ROUNDING 0.1f

const IVStyle TABLITSA_GROUPBOX_STYLE = TABLITSA_STYLE.WithRoundness(GROUP_BOX_ROUNDING).WithColor(EVColor::kFR, COLOR_WHITE).WithLabelText(DEFAULT_LABEL_TEXT.WithFGColor(IColor(255, 220, 200, 0))).WithDrawShadows(false);

class TablitsaVGroupControl : public IVGroupControl
{
public:
  TablitsaVGroupControl(const char* label="", const char* groupName="", float padL=0.f, float padT=0.f, float padR=0.f, float padB=0.f, const IVStyle& style = TABLITSA_GROUPBOX_STYLE) :
    IVGroupControl(label, groupName, padL, padT, padR, padB, style)
  {
    mLabelOffset = 0.f;
  }

  TablitsaVGroupControl(const IRECT& bounds, const char* label="", float labelOffset=0.f, const IVStyle& style=TABLITSA_GROUPBOX_STYLE) :
    IVGroupControl(bounds, label, labelOffset, style) {}
};

class TablitsaVTabBox : public TablitsaVGroupControl
{
public:
  TablitsaVTabBox(const IRECT& bounds, const char* label = "", const char* groupName="", bool active=false, float labelOffset = 0.f, const IVStyle& style = TABLITSA_GROUPBOX_STYLE) :
    TablitsaVGroupControl(bounds, label, labelOffset, style), mActive(active)
  {
    mIgnoreMouse = false;
    mGroupName.Set(groupName, MAX_PARAM_GROUP_LEN);
  }

  void OnInit() override {}

  void OnResize() override
  {
    TablitsaVGroupControl::OnResize();
    IRECT textBounds;
    GetUI()->MeasureText(mText, mLabelStr.Get(), textBounds);
    mWidgetBounds.ReduceFromTop(textBounds.H() / 1.5);
  }

  void SetActive(bool active = true)
  {
    mActive = active;
    // Show/hide controls in this tab view
    GetUI()->ForControlInGroup(mGroupName.Get(), [active](IControl& control) { control.Hide(!active); });
  }

  const IRECT& GetLabelBounds() const
  {
    return mLabelBounds;
  }

  const bool IsActive() const
  {
    return mActive;
  }

  const bool MouseIsOverLabel() const
  {
    return mMouseIsOverLabel;
  }

  void OnMouseOverLabel(float x, float y, const IMouseMod& mod)
  {
    bool prev = mMouseIsOverLabel;
    mMouseIsOverLabel = true;
    if (prev == false)
      SetDirty(false);
  }

  void OnMouseOutFromLabel()
  {
    bool prev = mMouseIsOverLabel;
    mMouseIsOverLabel = false;
    if (prev == true)
      SetDirty(false);
  }

  void Draw(IGraphics& g) override
  {
    DrawWidget(g);
    DrawLabel(g);
  }

  void DrawWidget(IGraphics& g) override
  {
    const float cr = GetRoundedCornerRadius(mWidgetBounds);
    const float ft = mStyle.frameThickness;
    const float hft = ft / 2.f;

    int nPaths = /*mStyle.drawShadows ? 2 :*/ 1;

    auto b = mWidgetBounds.GetPadded(/*mStyle.drawShadows ? -mStyle.shadowOffset :*/ 0.f);

    auto labelT = mLabelBounds.Empty() ? mRECT.MH() : mLabelBounds.T;
    auto labelB = mLabelBounds.Empty() ? mRECT.MH() : mLabelBounds.B;
    auto labelR = mLabelBounds.Empty() ? mRECT.MW() : mLabelBounds.R;
    auto labelL = mLabelBounds.Empty() ? mRECT.MW() : mLabelBounds.L;

    for (int i = 0; i < nPaths; i++)
    {
      const float offset = i == 0 ? 0.f : mStyle.shadowOffset;
      g.PathClear();
      // Path around label
      g.PathMoveTo(labelL, b.T + hft - offset);
      g.PathArc(labelL + cr, labelT + cr + hft - offset, cr, 270.f, 0.f);
      g.PathArc(labelR - cr, labelT + cr + hft - offset, cr, 0.f, 90.f);
      g.PathLineTo(labelR, b.T + hft - offset);
      if (mActive)
      {
        // Path around control
        g.PathArc(b.R - cr - hft - offset, b.T + cr + hft - offset, cr, 0.f, 90.f);
        g.PathArc(b.R - cr - hft - offset, b.B - cr - hft - offset, cr, 90.f, 180.f);
        g.PathArc(b.L + cr + hft - offset, b.B - cr - hft - offset, cr, 180.f, 270.f);
        g.PathArc(b.L + cr + hft - offset, b.T + cr + hft - offset, cr, 270.f, 360.f);
        g.PathLineTo(labelL, b.T + hft - offset);
        g.PathFill(IPattern(GetColor(EVColor::kFG).WithOpacity(0.5)), IFillOptions(true), &mBlend);
        g.PathStroke(mStyle.drawShadows ? GetColor(i == 0 ? kSH : kFR) : GetColor(kFR), ft);
      }
      else
      {
        g.PathFill(IPattern(GetColor(mMouseIsOverLabel ? EVColor::kHL : EVColor::kBG).WithOpacity(0.85)), IFillOptions(true), &mBlend);
        g.PathStroke(mStyle.drawShadows ? GetColor(i == 0 ? kSH : kFR) : GetColor(kFR), ft);
      }
    }
  }

  void DrawLabel(IGraphics& g)
  {
    IText labelText = mStyle.labelText.WithFGColor(mActive ? COLOR_WHITE : mMouseIsOverLabel ? GetColor(EVColor::kFG) : COLOR_WHITE.WithOpacity(0.8));
    if (mLabelBounds.H() && mStyle.showLabel)
    {
      IBlend blend = mControl->GetBlend();
      g.DrawText(labelText, mLabelStr.Get(), mLabelBounds, &blend);
    }
  }

protected:
  bool mActive;
  bool mMouseIsOverLabel{ false };
};

template <int MAXNC = 1>
class SequencerControl final : public IVMultiSliderControl<MAXNC>
{
public:
  SequencerControl(const IRECT& bounds, const char* label, const IVStyle& style = TABLITSA_STYLE.WithDrawFrame(false), int nSteps = 0, EDirection dir = EDirection::Vertical) :
    IVMultiSliderControl<MAXNC>(bounds, label, style, nSteps, dir)
  {
    SetColor(kX1, GetColor(kPR)); // Set active step color to the pressed step color
    SetColor(kHL, GetColor(kPR).WithOpacity(0.2)); // Background highlight
  }

  void DrawWidget(IGraphics& g) override
  {
    const int nVals = NVals();
    int nSteps = mNSteps + mZeroValueStepHasBounds;

    // Step labels
    IText labelText(mStyle.labelText);
    IRECT textBounds;
    g.MeasureText(labelText, "10", textBounds);
    const float stepHeight = IVTrackControlBase::mTrackBounds.Get()[0].H() / nSteps;
    // Make sure label fits inside the step borders
    while (textBounds.H() >= stepHeight)
    {
      labelText = labelText.WithSize(labelText.mSize - 0.5f);
      g.MeasureText(labelText, "10", textBounds);
    }

    for (int ch = 0; ch < nVals; ch++)
    {
      if (GetStepped())
      {
        for (int step{ 0 }; step < nSteps; ++step)
        {
          const IRECT trackBounds = IVTrackControlBase::mTrackBounds.Get()[ch];
          if (mDirection == EDirection::Vertical)
          {
            float stepSpan = trackBounds.H() / nSteps;
            const IRECT stepBounds = trackBounds.SubRect(EDirection::Vertical, nSteps, step);
            g.PathClear();
            g.PathMoveTo(trackBounds.L + 1, trackBounds.B - step * stepSpan);
            g.PathLineTo(trackBounds.R - 1, trackBounds.B - step * stepSpan);
            g.PathStroke(mStepMarkerColor, 1.f);
            if (step % 2 == 0)
              g.DrawText(labelText, std::to_string(nSteps - step - 1).c_str(), stepBounds);
          }
        }
      }
      DrawTrack(g, IVTrackControlBase::mTrackBounds.Get()[ch], ch);
    }
  }

private:
  const IColor mStepMarkerColor{ 50, 0, 0, 0 };
  const IColor mLabelColor = mStepMarkerColor.WithOpacity(0.75);
};

#define MAX_TAB_LABEL_LENGTH 32

class TablitsaEffectBankControl : public IControl
{
  struct TabText
  {
    char mText[MAX_TAB_LABEL_LENGTH];
  };

public:
  TablitsaEffectBankControl(const IRECT& bounds, std::initializer_list<char*> labels, std::initializer_list<char*> groupNames = { "" }, const IVStyle& style = TABLITSA_GROUPBOX_STYLE, const int maxTabs = 10) :
    IControl(bounds, kNoParameter), mMaxTabs(maxTabs), mStyle(style)
  {
    for (auto l : labels)
    {
      int nLabels = mLabels.GetSize();
      if (nLabels >= mMaxTabs)
        break;

      mLabels.Resize(nLabels + 1);
      TabText* pLab = mLabels.Get() + nLabels;
      strcpy(pLab->mText, l);
    }

    for (auto g : groupNames)
    {
      int nGroups = mGroups.GetSize();
      mGroups.Resize(nGroups + 1);
      TabText* pName = mGroups.Get() + nGroups;
      strcpy(pName->mText, g);
    }

    // If more labels were supplied than groups, append empty strings to the group  names
    if (labels.size() > groupNames.size())
    {
      for (auto i{ groupNames.size() }; i < labels.size(); ++i)
      {
        int nGroups = mGroups.GetSize();
        mGroups.Resize(nGroups + 1);
        TabText* pName = mGroups.Get() + nGroups;
        strcpy(pName->mText, "");
      }
    }
  }

  void OnMouseDown(float x, float y, const IMouseMod& mod) override
  {
    for (int i{ 0 }; i < mMaxTabs; ++i)
    {
      if (mTabs[i])
      {
        TablitsaVTabBox* tab = dynamic_cast<TablitsaVTabBox*>(mTabs[i]);
        if (!tab->IsActive() && tab->GetLabelBounds().Contains(x, y))
        {
          TabChanged(i);
        }
      }
    }
  }

  void OnMouseOver(float x, float y, const IMouseMod& mod) override
  {
    for (int i{ 0 }; i < mMaxTabs; ++i)
    {
      if (mTabs[i])
      {
        TablitsaVTabBox* tab = dynamic_cast<TablitsaVTabBox*>(mTabs[i]);
        if (!tab->IsActive() && tab->GetLabelBounds().Contains(x, y))
        {
          tab->OnMouseOverLabel(x, y, mod);
        }
        else
          tab->OnMouseOutFromLabel();
      }
    }
  }

  void OnMouseOut() override
  {
    for (int i{ 0 }; i < mMaxTabs; ++i)
    {
      if (mTabs[i])
      {
        TablitsaVTabBox* tab = dynamic_cast<TablitsaVTabBox*>(mTabs[i]);
        if (!tab->IsActive())
        {
          tab->OnMouseOutFromLabel();
        }
      }
    }
  }

  void Draw(IGraphics& g) override {}

  void TabChanged(int newIdx)
  {
    for (int i{ 0 }; i < mMaxTabs; ++i)
    {
      if (mTabs[i])
      {
        TablitsaVTabBox* tab = dynamic_cast<TablitsaVTabBox*>(mTabs[i]);
        tab->SetActive(i == newIdx);
      }
    }
  }

  void OnAttached() override
  {
    const int nLabels = mLabels.GetSize();
    float labelOffset = 0.f;
    for (int i{ 0 }; i < nLabels; ++i)
    {
      mTabs[i] = GetUI()->AttachControl(new TablitsaVTabBox(IRECT(mRECT), mLabels.Get()[i].mText, mGroups.Get()[i].mText, i==0, labelOffset, mStyle));
      labelOffset += dynamic_cast<TablitsaVTabBox*>(mTabs[i])->GetLabelBounds().W();
      // TODO: allow tabs to overlap if they don't all fit in the available space
    }
  }

private:
  const int mMaxTabs{ 10 };
  int MaxLabelLength{ 20 };
  WDL_TypedBuf<TabText> mLabels;
  WDL_TypedBuf<TabText> mGroups;

  IVStyle mStyle;
  IControl* mTabs[10]{ nullptr };
};

class DropdownListControl : public ICaptionControl
{
public:
  DropdownListControl(const IRECT& bounds, std::initializer_list<char*> options, const IText& text = TABLITSA_TEXT, const IColor& bgColor = DEFAULT_BGCOLOR, bool showLabel = false);

  void Draw(IGraphics& g) override;
  void OnResize() override;
  void OnMouseDown(float x, float y, const IMouseMod& mod) override;

  void AttachPopupMenu()
  {
    mMenu = new IPopupMenuControl();
    GetUI()->AttachControl(mMenu);
  }

protected:
  std::vector<std::string> mOptions;
  IPopupMenu mPopupMenu;
  IPopupMenuControl* mMenu;
  int mCurrentIdx{ 0 };
};

class PresetSelector : public ICaptionControl
{
public:
  PresetSelector(const IRECT& bounds, IPopupMenuControl* menu, std::initializer_list<char*> defaultPresets = { "" });

  void LoadUserPresets(std::initializer_list<char*> userPresets);
  void Draw(IGraphics& g) override;

protected:
  std::vector<std::string> mDefaultPresets;
  std::vector<std::string> mUserPresets;
  std::vector<std::string> mAllPresets;
  IPopupMenuControl* mMenu;
};

END_IPLUG_NAMESPACE
END_IGRAPHICS_NAMESPACE