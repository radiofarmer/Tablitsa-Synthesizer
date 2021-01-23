#pragma once
#include "Tablitsa.h"
#include "PeriodicTable.h"

/* A clone of the normal knob control, but with the ability to receive modulation parameters. */
class IVModKnobControl : public IVKnobControl
{
  static int mActiveIdx; // The parameter index of the IVModKnobControl currently linked to the modulator controls

public:
  /* Create a knob control with a modulatable value */
  IVModKnobControl(const IRECT& bounds, int paramIdx, const char* label = "", const IVStyle& style = DEFAULT_STYLE, bool valueIsEditable = false, double gearing = DEFAULT_GEARING)
    : IVKnobControl(bounds, paramIdx, label, style, valueIsEditable), mDefaultColor{ GetColor(kFG) }, mGearing(gearing), mModParamIdx(paramIdx)
  {
    GetMouseDblAsSingleClick();
  }

  /* Get modulator values from a different parameter. (For mutually-exclusive parameters) */
  void GetModulationFrom(int paramIdx)
  {
    mModParamIdx = paramIdx;
  }

  void OnMouseDown(float x, float y, const IMouseMod& mod) override
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
    if (mActive)
    {
      GetUI()->GetControlWithTag(kCtrlTagEnv1Depth)->SetParamIdx(kNoParameter);
      GetUI()->GetControlWithTag(kCtrlTagEnv2Depth)->SetParamIdx(kNoParameter);
      GetUI()->GetControlWithTag(kCtrlTagAmpEnvDepth)->SetParamIdx(kNoParameter);
      GetUI()->GetControlWithTag(kCtrlTagLFO1Depth)->SetParamIdx(kNoParameter);
      GetUI()->GetControlWithTag(kCtrlTagLFO2Depth)->SetParamIdx(kNoParameter);
    }
    else
    {
      // Set all modulator sliders to the values of the currently-selected modulated parameter
      GetUI()->GetControlWithTag(kCtrlTagEnv1Depth)->SetParamIdx(mModParamIdx + 1);
      GetUI()->GetControlWithTag(kCtrlTagEnv2Depth)->SetParamIdx(mModParamIdx + 2);
      GetUI()->GetControlWithTag(kCtrlTagAmpEnvDepth)->SetParamIdx(mModParamIdx + 3);
      GetUI()->GetControlWithTag(kCtrlTagLFO1Depth)->SetParamIdx(mModParamIdx + 4);
      GetUI()->GetControlWithTag(kCtrlTagLFO2Depth)->SetParamIdx(mModParamIdx + 5);
      mActiveIdx = GetParamIdx();
    }
    // Send values and change this control's active state
    GetDelegate()->SendCurrentParamValuesFromDelegate();
    mActive = !mActive;
  }

  void Draw(IGraphics& g) override
  {
    if (mActive && mActiveIdx == GetParamIdx())
      SetColor(kFG, GetColor(kPR));
    else
      SetColor(kFG, mDefaultColor);
    DrawBackground(g, mRECT);
    DrawLabel(g);
    DrawWidget(g);
    DrawValue(g, mValueMouseOver);
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
    DrawPressableShape(g, /*mShape*/ EVShape::Ellipse, knobHandleBounds, mMouseDown, mMouseIsOver, IsDisabled());
    DrawIndicatorTrack(g, angle, cx, cy, widgetRadius);
    DrawPointer(g, angle, cx, cy, knobHandleBounds.W() / 2.f);
  }

  void DrawIndicatorTrack(IGraphics& g, float angle, float cx, float cy, float radius) override
  {
    // Set the origin of the track arch to the center of the range for dials with negative minimum values
    if (GetParam()->GetMin() < 0.)
      mAnchorAngle = 0.;
    if (mTrackSize > 0.f)
    {
      g.DrawArc(IColor(100, 0, 0, 0), cx, cy, radius, angle >= mAnchorAngle ? mAnchorAngle : mAnchorAngle - (mAnchorAngle - angle), angle >= mAnchorAngle ? angle : mAnchorAngle, &mBlend, mTrackSize);

      // Envelope 1
      float modAngle = std::max(std::min(static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 1)->Value()) * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[0], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 3);
      // Envelope 2
      radius -= mTrackToHandleDistance / 2.f;
      modAngle = std::max(std::min(static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 2)->Value()) * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[1], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 2);
      // AmpEnv
      radius -= mTrackToHandleDistance / 4.f;
      modAngle = std::max(std::min(static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 3)->Value()) * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[2], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
      // LFO 1
      radius -= 3.f * mTrackToHandleDistance / 4.f;
      modAngle = std::max(std::min(static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 4)->Value()) * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[3], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 2);
      modAngle = std::max(std::min(static_cast<float>(-1. * GetDelegate()->GetParam(mModParamIdx + 4)->Value()) * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(IColor::LinearInterpolateBetween(mModArcColor[3], IColor(200, 255, 255, 255), 0.4), cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize * 2);
      // LFO 2
      radius += mTrackToHandleDistance / 4.f;
      modAngle = std::max(std::min(static_cast<float>(GetDelegate()->GetParam(mModParamIdx + 5)->Value()) * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(mModArcColor[4], cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
      modAngle = std::max(std::min(static_cast<float>(-1. * GetDelegate()->GetParam(mModParamIdx + 5)->Value()) * (mAngle2 - mAngle1) + angle, mAngle2), mAngle1);
      g.DrawArc(IColor::LinearInterpolateBetween(mModArcColor[4], IColor(200, 255, 255, 255), 0.4), cx, cy, radius, modAngle >= angle ? angle : modAngle, modAngle >= angle ? modAngle : angle, &mBlend, mTrackSize);
    }
  }

protected:
  const IColor mDefaultColor;
  static inline const IColor mModArcColor[5]{
    {200, 200, 100, 100},
    {200, 225, 150, 100},
    {200, 250, 200, 100},
    {200, 160, 0, 225},
    {200, 225, 0, 190},
  };
  bool mActive{ false };
  double mGearing;
  int mModParamIdx;
};

class ModSliderControl : public IVSliderControl
{
public:
  ModSliderControl(const IRECT& bounds, int paramIdx = kNoParameter, const char* label = "", const IVStyle& style = DEFAULT_STYLE, bool valueIsEditable = false, EDirection dir = EDirection::Vertical, double gearing = DEFAULT_GEARING, float handleSize = 8.f, float trackSize = 2.f, bool handleInsideTrack = true) :
    IVSliderControl(bounds, paramIdx, label, style, valueIsEditable, dir, gearing, handleSize, trackSize, handleInsideTrack)
  {
  }

  void Draw(IGraphics& g) override
  {
    if (GetParamIdx() == kNoParameter)
    {
      SetValueStr("N/A");
      SetShowValue(true);
      SetDisabled(true);
    }
    else
      SetDisabled(false);
    IVSliderControl::Draw(g);
  }

  void OnMouseDown(float x, float y, const IMouseMod& mod) override
  {
    if (!mod.L)
    {
      SetValue(0.5);
      ISliderControlBase::mMouseDown = true;
      SetDirty(true);
      IControl::OnMouseDown(x, y, mod);
    }
    else
      IVSliderControl::OnMouseDown(x, y, mod);
  }

  // Overridden drawing function to draw the filled area starting in the middle of the track
  void DrawWidget(IGraphics& g) override
  {
    float value = (float)GetValue(); // NB: Value is normalized to between 0. and 1.
    const IRECT handleBounds = (GetParamIdx() == kNoParameter) ? mTrackBounds.FracRect(mDirection, 0.f) : mTrackBounds.FracRect(mDirection, value);
    const IRECT filledTrack = (GetParamIdx() == kNoParameter) ?  handleBounds : (value >= 0.5f) ? mTrackBounds.GetGridCell(0, 0, 2, 1).FracRect(mDirection, 2.f * (value - 0.5f)) : mTrackBounds.GetGridCell(1, 0, 2, 1).FracRect(mDirection, 2.f * (0.5 - value), true);
    
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
    const IRECT adjustedFillBounds = filledArea;
    //const IRECT adjustedFillBounds = mDirection == EDirection::Vertical ? filledArea.GetVPadded(extra) : filledArea.GetHPadded(extra);
    const float cr = GetRoundedCornerRadius(mTrackBounds);

    g.FillRoundRect(GetColor(kSH), adjustedTrackBounds, cr, &mBlend);
    g.FillRoundRect(GetColor(kX1), adjustedFillBounds, cr, &mBlend);

    if (mStyle.drawFrame)
      g.DrawRoundRect(GetColor(kFR), adjustedTrackBounds, cr, &mBlend, mStyle.frameThickness);
  }

private:
};

int IVModKnobControl::mActiveIdx = -1;

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

  void Draw(IGraphics& g)
  {
    g.FillRoundRect(IColor(255, 20, 0, 45), mRECT);
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
    DrawElement(g, mRECT.GetFromRight(80.f).GetGridCell(0, 0, 2, 1).GetCentredInside(70.f, 80.f), mSelectedElements[0], 0);
    DrawElement(g, mRECT.GetFromRight(80.f).GetGridCell(1, 0, 2, 1).GetCentredInside(70.f, 80.f), mSelectedElements[1], 1);

    g.DrawSVG(mSVG, mRECT);
  }

  void DrawElement(IGraphics& g, IRECT& bounds, int atomicNumber, int idx)
  {
    assert(atomicNumber > 0);
    IColor col{ TablitsaDSP<double>::tableLoading[idx] ? ElementIconColor[idx].WithOpacity(0.5f) : ElementIconColor[idx] };

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
      // Send new wavetable value
      SetValue(static_cast<double>(mSelectedElements[mIsDragging]) / 118., mIsDragging);
//      GetDelegate()->SendParameterValueFromUI(GetParamIdx(mIsDragging), static_cast<double>(mSelectedElements[mIsDragging] - 1) / 118.);
//      GetDelegate()->OnParamChange(GetParamIdx(mIsDragging));
    };
    mIsDragging = -1;
    SetDirty(true);
  }

private:
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
};

template <int MAXNC = 1>
class SequencerControl : public IVMultiSliderControl<MAXNC>
{
public:
  SequencerControl(const IRECT& bounds, const char* label, const IVStyle& style = DEFAULT_STYLE, int nSteps = 0, EDirection dir = EDirection::Vertical) :
    IVMultiSliderControl<MAXNC>(bounds, label, style, nSteps, dir)
  {

  }
};