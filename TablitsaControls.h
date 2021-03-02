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

const IVStyle modKnobStyle{ TABLITSA_STYLE.WithShowLabel(true).WithColors(knobColorSpec).WithLabelText(TABLITSA_STYLE.labelText.WithSize(17.f)) };

const IVStyle toggleStyle{ TABLITSA_STYLE.WithDrawFrame(true).WithColor(EVColor::kFG, COLOR_TRANSPARENT).WithColor(EVColor::kHL, IColor(100, 220, 0, 180)).WithShowLabel(false) };

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
  void DrawWidget(IGraphics& g) override;
  void DrawTrack(IGraphics& g, const IRECT& filledArea) override;

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
  void DrawPointer(IGraphics& g, float angle, float cx, float cy, float radius);
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

  void LoadModParams();

  void DrawWidget(IGraphics& g) override;

  void Draw(IGraphics& g) override
  {
    DrawBackground(g, mRECT);
    DrawLabel(g);
    DrawWidget(g);
    DrawValue(g, mValueMouseOver);
  }

  void DrawIndicatorTrack(IGraphics& g, float angle, float cx, float cy, float radius) override;

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
    {{255, 255, 50, 50}, {}},
    {{255, 225, 150, 100}, {}},
    {{255, 255, 255, 100}, {}},
    {{255, 180, 0, 225}, {200, 200, 0, 235}},
    {{200, 50, 0, 220}, {200, 235, 0, 220}},
    {{200, 0, 250, 100}, {}},
    {{255, 0, 100, 255}, {}},
    {{150, 200, 0, 200}, {}},
    {{100, 255, 0, 0}, {}}
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
    mTableWidth = std::ceilf(TableWidth * mScaleFact) + TABLE_ADJ;
    mTableHeight = std::ceilf(TableHeight * mScaleFact) + TABLE_ADJ;
    mTableTLHC[0] = mRECT.L + 37.f * mScaleFact;
    mTableTLHC[1] = mRECT.T + 37.f * mScaleFact;
    // Lanthanides and Actinides
    mLaAcWidth = std::ceilf(LaAcWidth * mScaleFact) + TABLE_ADJ;
    mLaAcHeight = std::ceilf(LaAcHeight * mScaleFact) + TABLE_ADJ;
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
  
  void Draw(IGraphics& g);
  void DrawElement(IGraphics& g, const IRECT& bounds, int atomicNumber, int idx);
  void OnMouseOver(float x, float y, const IMouseMod& mod);
  void OnMouseOut() {}
  void OnMouseDrag(float x, float y, float dx, float dy, const IMouseMod& mod) { OnMouseOver(x, y, mod); }
  void OnMouseDown(float x, float y, const IMouseMod& mod);
  void OnMouseUp(float x, float y, const IMouseMod& mod);
  void SetTableLoading(const bool isLoading, const int tableIdx) { mTableLoading[tableIdx] = isLoading; }

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

template <int MAXNC = 1>
class SequencerControl final : public IVMultiSliderControl<MAXNC>
{
public:
  SequencerControl(const IRECT& bounds, const char* label, const IVStyle& style = TABLITSA_STYLE.WithDrawFrame(false), int nSteps = 0, EDirection dir = EDirection::Vertical) :
    IVMultiSliderControl<MAXNC>(bounds, label, style, nSteps, dir)
  {
    SetColor(kX1, GetColor(kPR)); // Set active step color to the pressed step color
    SetColor(kHL, GetColor(kPR).WithOpacity(0.2f)); // Background highlight
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

#define GROUP_BOX_ROUNDING 0.1f

const IVStyle TABLITSA_GROUPBOX_STYLE = TABLITSA_STYLE.WithRoundness(GROUP_BOX_ROUNDING).WithColor(EVColor::kFR, COLOR_WHITE).WithLabelText(DEFAULT_LABEL_TEXT.WithFGColor(IColor(255, 220, 200, 0))).WithDrawShadows(false);

class TablitsaVGroupControl : public IVGroupControl
{
public:
  TablitsaVGroupControl(const char* label = "", const char* groupName = "", float padL = 0.f, float padT = 0.f, float padR = 0.f, float padB = 0.f, const IVStyle& style = TABLITSA_GROUPBOX_STYLE) :
    IVGroupControl(label, groupName, padL, padT, padR, padB, style)
  {
    mLabelOffset = 0.f;
  }

  TablitsaVGroupControl(const IRECT& bounds, const char* label = "", float labelOffset = 0.f, const IVStyle& style = TABLITSA_GROUPBOX_STYLE) :
    IVGroupControl(bounds, label, labelOffset, style) {}
};

#define MAX_TAB_LABEL_LENGTH 32

class TablitsaVTabBox : public TablitsaVGroupControl
{
public:
  TablitsaVTabBox(const IRECT& bounds, const char* label = "", const char* groupName = "", bool active = false, float labelOffset = 0.f, const IVStyle& style = TABLITSA_GROUPBOX_STYLE) :
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

  void SetActive(bool active = true);

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

  void DrawWidget(IGraphics& g) override;

  void DrawLabel(IGraphics& g)
  {
    IText labelText = mStyle.labelText.WithFGColor(mActive ? COLOR_WHITE : mMouseIsOverLabel ? GetColor(EVColor::kFG) : COLOR_WHITE.WithOpacity(0.8));
    if (mLabelBounds.H() && mStyle.showLabel)
    {
      IBlend blend = mControl->GetBlend();
      g.DrawText(labelText, mLabelStr.Get(), mLabelBounds, &blend);
    }
  }

  void SetGroupName(const char* newGroupName);

protected:
  bool mActive;
  bool mMouseIsOverLabel{ false };
};

class TablitsaEffectBankControl : public IControl
{
  struct TabText
  {
    char mText[MAX_TAB_LABEL_LENGTH];
  };

public:
  TablitsaEffectBankControl(const IRECT& bounds, std::initializer_list<char*> labels, std::initializer_list<char*> groupNames = { "" }, const IVStyle& style = TABLITSA_GROUPBOX_STYLE, const int maxTabs = 10);

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

  void TabChanged(int newIdx);

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

  /* TablitsaEffectsBankControl */
  void SetTabGroup(int tabIdx, const char* newGroupName)
  {
    dynamic_cast<TablitsaVTabBox*>(mTabs[tabIdx])->SetGroupName(newGroupName);
  }

private:
  const int mMaxTabs{ 10 };
  int MaxLabelLength{ 20 };
  WDL_TypedBuf<TabText> mLabels;
  WDL_TypedBuf<TabText> mGroups;

  IVStyle mStyle;
  IControl* mTabs[10]{ nullptr }; // Needs to be an IControl pointer in order to attach (or maybe you can down-cast to IControl*?)
};

class DropdownListControl : public ICaptionControl
{
public:
  DropdownListControl(const IRECT& bounds, std::initializer_list<char*> options, const IText& text = TABLITSA_TEXT, const IColor& bgColor = DEFAULT_BGCOLOR, bool showLabel = false);

  void Draw(IGraphics& g) override;
  void OnResize() override;
  void OnMouseDown(float x, float y, const IMouseMod& mod) override;

  /* DropdownList Control */
  int GetCurrentIndex() { return mCurrentIdx; }
  const char* GetSelectedString() { return mOptions[mCurrentIdx].c_str(); }

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

class ModPlotControl : public virtual IVPlotControl
{
public:
  ModPlotControl(const IRECT& bounds, double* table, const int tableSize, const int numPoints, const IColor& color = TABLITSA_STYLE.colorSpec.mColors[EVColor::kHL], const IVStyle& style = TABLITSA_STYLE, float gearing=4.f);

  ModPlotControl(const IRECT& bounds, int paramIdx, double* table, const int tableSize, const int numPoints, const IColor& color = TABLITSA_STYLE.colorSpec.mColors[EVColor::kHL], const IVStyle& style = TABLITSA_STYLE, float gearing = 4.f);

  void SetPlotTable(const double* pTable);
  void OnMouseDrag(float x, float y, float dX, float dY, const IMouseMod& mod) override;

  ~ModPlotControl() { if (mEmptyTable) delete[] mEmptyTable; }

protected:
  const double* mTable;
  double* mEmptyTable{ nullptr };
  float mGearing;
  unsigned int mTablePhase{ 0 };
  unsigned int mTableSize;
};

END_IPLUG_NAMESPACE
END_IGRAPHICS_NAMESPACE