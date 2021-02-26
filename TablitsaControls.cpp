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

/* Tab Controls */

void TablitsaVTabBox::SetActive(bool active)
{
  mActive = active;
  // Show/hide controls in this tab view
  GetUI()->ForControlInGroup(mGroupName.Get(), [active](IControl& control) {
      control.Hide(!active);
    });
}

void TablitsaVTabBox::SetGroupName(const char* newGroupName)
{
  // Hide old group
  GetUI()->ForControlInGroup(mGroupName.Get(), [](IControl& control) {
    if (!control.IsHidden())
      control.Hide(true);
    });
  mGroupName.Set(newGroupName);
  // show new group
  GetUI()->ForControlInGroup(mGroupName.Get(), [](IControl& control) {
    if (control.IsHidden())
      control.Hide(false);
    control.SetDirty(true); // Trigger control show/hide functions
    });
}

/* Effect Bank */

TablitsaEffectBankControl::TablitsaEffectBankControl(const IRECT& bounds, std::initializer_list<char*> labels, std::initializer_list<char*> groupNames, const IVStyle& style, const int maxTabs) :
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
      strcpy(pName->mText, mGroups.Get()[i].mText);
    }
  }
}

void TablitsaEffectBankControl::TabChanged(int newIdx)
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

/* Generic Dropdown List */

DropdownListControl::DropdownListControl(const IRECT& bounds, std::initializer_list<char*> options, const IText& text, const IColor& bgColor, bool showLabel) :
  ICaptionControl(bounds, kNoParameter, text, bgColor, showLabel)
{
  for (auto s : options)
  {
    mOptions.push_back(std::string(s));
  }

  mPopupMenu.SetFunction([this](IPopupMenu* pControl) {
    if (pControl->GetChosenItemIdx() >= 0)
      mCurrentIdx = pControl->GetChosenItemIdx();
    this->SetDirty(true); // Tablitsa: Trigger action function to update the tab control
    });
}

void DropdownListControl::Draw(IGraphics& g)
{
  mStr.Set(mOptions[mCurrentIdx].c_str());

  ITextControl::Draw(g);

  if (mTri.W() > 0.f)
  {
    g.FillTriangle(mMouseIsOver ? mTriangleMouseOverColor : mTriangleColor, mTri.L, mTri.T, mTri.R, mTri.T, mTri.MW(), mTri.B, GetMouseIsOver() ? 0 : &BLEND_50);
  }
}

void DropdownListControl::OnMouseDown(float x, float y, const IMouseMod& mod)
{
  mPopupMenu.Clear();
  int nDisplayTexts = mOptions.size();
  // Fill the menu
  for (int i = 0; i < nDisplayTexts; ++i)
  {
    const char* str = &mOptions[i].c_str()[0];
    // TODO: what if two parameters have the same text?
    if (!strcmp(str, mOptions[mCurrentIdx].c_str())) // strings are equal
      mPopupMenu.AddItem(new IPopupMenu::Item(str, IPopupMenu::Item::kChecked), -1);
    else // not equal
      mPopupMenu.AddItem(new IPopupMenu::Item(str), -1);

    mPopupMenu.SetRootTitle(mOptions[mCurrentIdx].c_str());
  }
  mMenu->CreatePopupMenu(mPopupMenu, mRECT);
}

void DropdownListControl::OnResize()
{
  if (mOptions.size() != 0)
  {
    mTri = mRECT.FracRectHorizontal(0.2f, true).GetCentredInside(IRECT(0, 0, 8, 5)); //TODO: This seems rubbish
  }
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


DropdownListControl::DropdownListControl(const IRECT& bounds, std::initializer_list<char*> options, const IText& text, const IColor& bgColor, bool showLabel) :
  ICaptionControl(bounds, kNoParameter, text, bgColor, showLabel)
{
  for (auto s : options)
  {
    mOptions.push_back(std::string(s));
  }

  mPopupMenu.SetFunction([this](IPopupMenu* pControl) {
    mCurrentIdx = pControl->GetChosenItemIdx();
    });
}

void DropdownListControl::Draw(IGraphics& g)
{
  mStr.Set(mOptions[0].c_str());

  ITextControl::Draw(g);

  if (mTri.W() > 0.f)
  {
    g.FillTriangle(mMouseIsOver ? mTriangleMouseOverColor : mTriangleColor, mTri.L, mTri.T, mTri.R, mTri.T, mTri.MW(), mTri.B, GetMouseIsOver() ? 0 : &BLEND_50);
  }
}

void DropdownListControl::OnMouseDown(float x, float y, const IMouseMod& mod)
{
  mPopupMenu.Clear();
  int nDisplayTexts = mOptions.size();
  // Fill the menu
  for (int i = 0; i < nDisplayTexts; ++i)
  {
    const char* str = &mOptions[i].c_str()[0];
    // TODO: what if two parameters have the same text?
    if (!strcmp(str, mOptions[mCurrentIdx].c_str())) // strings are equal
      mPopupMenu.AddItem(new IPopupMenu::Item(str, IPopupMenu::Item::kChecked), -1);
    else // not equal
      mPopupMenu.AddItem(new IPopupMenu::Item(str), -1);

    mPopupMenu.SetRootTitle(mOptions[mCurrentIdx].c_str());
  }
  mMenu->CreatePopupMenu(mPopupMenu, mRECT);
}

void DropdownListControl::OnResize()
{
  if (mOptions.size() != 0)
  {
    mTri = mRECT.FracRectHorizontal(0.2f, true).GetCentredInside(IRECT(0, 0, 8, 5)); //TODO: This seems rubbish
  }
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