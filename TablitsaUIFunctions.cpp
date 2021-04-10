#include "TablitsaUIFunctions.h"


void PercentDisplayFunc(double value, WDL_String& str)
{
  int precision = (value < 0.1) + (value < 0.01);
  double displayValue = value * 100.;
  str.SetFormatted(MAX_PARAM_DISPLAY_LEN, "%.*f", precision, displayValue);
  str.Append(" %");
}

void DelayDisplayFunc(double value, WDL_String& str)
{
  int precision = 2;
  double displayValue = value;
  str.SetFormatted(MAX_PARAM_DISPLAY_LEN, "%.*f", precision, displayValue);
  str.Append(" ms");
}

void SetAllEffectControlsDirty(IGraphics* pGraphics, int idx)
{
  std::vector<IControl*> controls{
    pGraphics->GetControlWithTag(idx),
    pGraphics->GetControlWithTag(idx + 1),
    pGraphics->GetControlWithTag(idx + 2),
    pGraphics->GetControlWithTag(idx + 3),
    pGraphics->GetControlWithTag(idx + 4),
    pGraphics->GetControlWithTag(idx + 5)
  };
  for (auto* ctrl : controls)
    ctrl->SetDirty(true);
}

/* COEFICIENT MODULATOR (VOICE) */

void InitCoefModUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset)
{
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  // Save current values
  double paramVals[4]{};
  if (!reset)
  {
    for (int i{0}; i < 4; ++i)
      paramVals[i] = pPlugin->GetParam(params[i])->Value();
  }
  pPlugin->GetParam(params[0])->InitDouble(paramNames[0], 0., 0., 1., 0.01);
  pPlugin->GetParam(params[1])->InitDouble(paramNames[1], 0., -2., 2., 0.01, "8va.");
  pPlugin->GetParam(params[2])->InitDouble(paramNames[2], 0., 0., 1., 0.01, "Cycles");
  pPlugin->GetParam(params[3])->InitDouble(paramNames[3], 0., 0., 1., 0.01);
  // Reload original values if not resetting
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      pPlugin->GetParam(params[i])->Set(paramVals[i]);
  }

  pPlugin->GetParam(params[0])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[1])->SetDisplayFunc(nullptr);
  pPlugin->GetParam(params[2])->SetDisplayFunc(nullptr);
  pPlugin->GetParam(params[3])->SetDisplayFunc(PercentDisplayFunc);

  // Reset param 1 display texts
  auto* p1 = pPlugin->GetParam(params[0]);
  p1->SetDisplayPrecision(2);

  for (int i{ 0 }; i < TABLITSA_EFFECT_PARAMS; ++i)
    controls[i]->SetValue(pPlugin->GetParam(params[i])->GetNormalized());

  controls[0]->SetDisabled(false);
  controls[1]->SetDisabled(false);
  controls[2]->SetDisabled(false);
  controls[3]->SetDisabled(false);
  pGraphics->HideControl(params[4], true);
  pGraphics->HideControl(params[5], true);
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Mod Depth");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Pitch");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Phase");
  dynamic_cast<IVKnobControl*>(controls[3])->SetLabelStr("Mix");
  // Modulation ON for knob 1
  dynamic_cast<TablitsaIVModKnobControl*>(controls[0])->EnableModulation(true);
  // Toggle action functions
  controls[4]->SetActionFunction(nullptr);
  controls[5]->SetActionFunction(nullptr);
  controls[4]->Hide(true);
  controls[5]->Hide(true);
}

/* DELAY (MASTER) */

void DelayTempoSyncToggle(Plugin* pPlugin, IGraphics* pGraphics, const std::vector<int>& params, bool isTempoSync, const bool reset)
{
  dynamic_cast<Tablitsa*>(pPlugin)->SetDelayTempoSync(isTempoSync);
  if (isTempoSync)
  {
    const IParam* oldParams[]{ pPlugin->GetParam(params[0]), pPlugin->GetParam(params[1]) };
    pPlugin->GetParam(params[0])->InitEnum(pPlugin->GetParam(params[0])->GetName(), reset ? DelayEffect<sample>::k8th : oldParams[0]->Value(), { DELAY_TEMPODIV_VALIST });
    pPlugin->GetParam(params[1])->InitEnum(pPlugin->GetParam(params[1])->GetName(), reset ? DelayEffect<sample>::k8th : oldParams[1]->Value(), { DELAY_TEMPODIV_VALIST });
  }
  else
  {
    const IParam* oldParams[]{ pPlugin->GetParam(params[0]), pPlugin->GetParam(params[1]) };
    pPlugin->GetParam(params[0])->InitDouble(oldParams[0]->GetName(), reset ? 100. : oldParams[0]->Value(), 1., TABLITSA_MAX_EFFECT_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
    pPlugin->GetParam(params[1])->InitDouble(oldParams[1]->GetName(), reset ? 100. : oldParams[1]->Value(), 1., TABLITSA_MAX_EFFECT_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
  }
  pGraphics->SetAllControlsDirty();
}

void DelayStereoToggle(IGraphics* pGraphics, const std::vector<IControl*>& controls, const bool mono)
{
  auto leftChannelLock = IActionFunction([pGraphics, controls](IControl* pControl) {
    // Set right channel control's value to this control's value
    controls[1]->SetValue(pControl->GetValue());
    pGraphics->GetDelegate()->SendParameterValueFromUI(controls[1]->GetParamIdx(), pControl->GetValue());
    controls[1]->SetDirty(false);
    });
  auto rightChannelLock = IActionFunction([pGraphics, controls](IControl* pControl) {
    // Set right channel control's value to this control's value
    controls[0]->SetValue(pControl->GetValue());
    pGraphics->GetDelegate()->SendParameterValueFromUI(controls[0]->GetParamIdx(), pControl->GetValue());
    controls[0]->SetDirty(false);
    });

  if (mono)
  {
    controls[0]->SetActionFunction(leftChannelLock);
    controls[1]->SetActionFunction(rightChannelLock);
  }
  else
  {
    controls[0]->SetActionFunction([pGraphics](IControl* pControl) {});
    controls[1]->SetActionFunction([pGraphics](IControl* pControl) {});
  }
  pGraphics->SetAllControlsDirty();
}

void InitDelayUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset)
{
  // Save current values
  double paramVals[4]{};
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      paramVals[i] = pPlugin->GetParam(params[i])->Value();
  }
  // Adjust param ranges
  pPlugin->GetParam(params[4])->InitBool(paramNames[2], reset ? false : pPlugin->GetParam(params[4])->Value() > 0.5);
  pPlugin->GetParam(params[5])->InitBool(paramNames[3], reset ? false : pPlugin->GetParam(params[5])->Value() > 0.5);
  bool tempoSync = pPlugin->GetParam(params[4])->Value();
  if (tempoSync)
  {
    pPlugin->GetParam(params[0])->InitEnum(paramNames[0], DelayEffect<sample>::k8th, { DELAY_TEMPODIV_VALIST });
    pPlugin->GetParam(params[1])->InitEnum(paramNames[1], DelayEffect<sample>::k8th, { DELAY_TEMPODIV_VALIST });
  }
  else
  {
    pPlugin->GetParam(params[0])->InitDouble(paramNames[0], 100., 1., TABLITSA_MAX_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
    pPlugin->GetParam(params[1])->InitDouble(paramNames[1], 100., 1., TABLITSA_MAX_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
  }
  pPlugin->GetParam(params[2])->InitDouble(paramNames[2], 0., 0., 1., 0.01);
  pPlugin->GetParam(params[3])->InitDouble(paramNames[3], 0., 0., 1., 0.01);
  // Reload original values if not resetting
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      pPlugin->GetParam(params[i])->Set(paramVals[i]);
  }

  pPlugin->GetParam(params[2])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[3])->SetDisplayFunc(PercentDisplayFunc);

  for (int i{ 0 }; i < TABLITSA_EFFECT_PARAMS; ++i)
    controls[i]->SetValue(pPlugin->GetParam(params[i])->GetNormalized());

  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Left");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Right");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Feedback");
  dynamic_cast<IVKnobControl*>(controls[3])->SetLabelStr("Mix");

  // Toggle Control Action Functions
  controls[4]->SetActionFunction(
    [pGraphics, pPlugin, params, controls](IControl* pControl) {
      bool delayIsSynced = pControl->GetValue() > 0.5;
      std::vector<int> syncParams{ params[0], params[1] };
      DelayTempoSyncToggle(pPlugin, pGraphics, syncParams, delayIsSynced);
    });
  controls[5]->SetActionFunction(
    [pGraphics, controls](IControl* pControl) {
      bool mono = pControl->GetValue() > 0.5;
      std::vector<IControl*> channelControls{ controls[0], controls[1] };
      DelayStereoToggle(pGraphics, channelControls, mono);
    });
  // Show toggle buttons
  controls[4]->Hide(false);
  controls[5]->Hide(false);
  controls[4]->SetDisabled(false);
  controls[5]->SetDisabled(false);
  // Enable all knobs
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  for (auto* knob : allKnobs)
    knob->SetDisabled(false);
}

/* DISTORTION (VOICE) */

void InitDistortionUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset)
{
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  for (auto* knob : allKnobs)
    knob->SetDisabled(false);
  // Save current values
  double paramVals[4]{};
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      paramVals[i] = pPlugin->GetParam(params[i])->Value();
  }
  pPlugin->GetParam(params[0])->InitDouble(paramNames[0], 0., 0., 1., 0.01);
  pPlugin->GetParam(params[1])->InitDouble(paramNames[1], 0., 0., 1., 0.01);
  pPlugin->GetParam(params[2])->InitDouble(paramNames[2], 0., TABLITSA_DISTORTION_FREQ_LOW, TABLITSA_DISTORTION_FREQ_HIGH, 0.01, "", 0, "", IParam::ShapeExp());
  pPlugin->GetParam(params[3])->InitDouble(paramNames[3], 0., 0., 1., 0.01);
  // Reload original values if not resetting
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      pPlugin->GetParam(params[i])->Set(paramVals[i]);
  }

  pGraphics->HideControl(params[4], true);
  pGraphics->HideControl(params[5], true);

  pPlugin->GetParam(params[0])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[1])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[2])->SetDisplayFunc([pPlugin](double value, WDL_String& str) {
    str.SetFormatted(MAX_PARAM_DISPLAY_LEN, "%.0f", pPlugin->GetSampleRate() * value);
    str.Append(" Hz");
    });
  pPlugin->GetParam(params[3])->SetDisplayFunc(PercentDisplayFunc);

  for (int i{ 0 }; i < TABLITSA_EFFECT_PARAMS; ++i)
    controls[i]->SetValue(pPlugin->GetParam(params[i])->GetNormalized());

  controls[1]->SetDisabled(false);
  controls[2]->SetDisabled(false);
  // Labels
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Gain");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Fuzz");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Color");
  dynamic_cast<IVKnobControl*>(controls[3])->SetLabelStr("Mix");
  // Modulation off for knob 1
  dynamic_cast<TablitsaIVModKnobControl*>(controls[0])->EnableModulation(true);
  // Toggle action functions
  controls[4]->SetActionFunction(nullptr);
  controls[5]->SetActionFunction(nullptr);
  controls[4]->Hide(true);
  controls[5]->Hide(true);
}

/* EQ (MASTER) */

void InitEqualizerUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset)
{
  // Save current values
  double paramVals[4]{};
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      paramVals[i] = pPlugin->GetParam(params[i])->Value();
  }
  pPlugin->GetParam(params[0])->InitDouble(paramNames[0], 1., 0., 2., 0.01);
  pPlugin->GetParam(params[1])->InitDouble(paramNames[1], 1., 0., 2., 0.01);
  pPlugin->GetParam(params[2])->InitDouble(paramNames[2], 0.5, 0., 1., 0.01);
  pPlugin->GetParam(params[3])->InitDouble(paramNames[3], 1., 0., 2., 0.01);
  // Reload original values if not resetting
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      pPlugin->GetParam(params[i])->Set(paramVals[i]);
  }

  pPlugin->GetParam(params[2])->SetDisplayFunc(PercentDisplayFunc);

  for (int i{ 0 }; i < TABLITSA_EFFECT_PARAMS; ++i)
    controls[i]->SetValue(pPlugin->GetParam(params[i])->GetNormalized());

  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Low Level");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Mid Level");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Mid Pos");
  dynamic_cast<IVKnobControl*>(controls[3])->SetLabelStr("High Level");

  controls[2]->Hide(false);
  controls[3]->Hide(false);
  controls[0]->SetDisabled(false);
  controls[1]->SetDisabled(false);
  controls[2]->SetDisabled(false);
  controls[3]->SetDisabled(false);
  // Toggles
  controls[4]->Hide(true);
  controls[5]->Hide(true);
  controls[4]->SetActionFunction(nullptr);
  controls[5]->SetActionFunction(nullptr);
}

/* REVERB 1 & 2 (MASTER) */

void InitReverbUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset)
{
  // Save current values
  double paramVals[4]{};
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      paramVals[i] = pPlugin->GetParam(params[i])->Value();
  }
  pPlugin->GetParam(params[0])->InitDouble(paramNames[0], 0.5, 0., 1., 0.01);
  pPlugin->GetParam(params[1])->InitDouble(paramNames[1], 0.5, 0., 1., 0.01);
  pPlugin->GetParam(params[2])->InitDouble(paramNames[2], 0.5, 0., 1., 0.01);
  pPlugin->GetParam(params[3])->InitDouble(paramNames[3], 0., 0., 1., 0.01);
  // Reload original values if not resetting
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      pPlugin->GetParam(params[i])->Set(paramVals[i]);
  }

  pPlugin->GetParam(params[0])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[1])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[2])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[3])->SetDisplayFunc(PercentDisplayFunc);

  for (int i{ 0 }; i < TABLITSA_EFFECT_PARAMS; ++i)
    controls[i]->SetValue(pPlugin->GetParam(params[i])->GetNormalized());

  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Echo Time");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Room Size");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Damping");
  dynamic_cast<IVKnobControl*>(controls[3])->SetLabelStr("Mix");

  controls[2]->Hide(false);
  controls[3]->Hide(false);
  controls[0]->SetDisabled(false);
  controls[1]->SetDisabled(false);
  controls[2]->SetDisabled(false);
  controls[3]->SetDisabled(false);
  // Toggles
  controls[4]->Hide(true);
  controls[5]->Hide(true);
  controls[4]->SetActionFunction(nullptr);
  controls[5]->SetActionFunction(nullptr);
}

void InitReverb2UI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset)
{
  double paramVals[4]{};
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      paramVals[i] = pPlugin->GetParam(params[i])->Value();
  }
  pPlugin->GetParam(params[0])->InitDouble(paramNames[0], 0.5, 0., 1., 0.01);
  pPlugin->GetParam(params[1])->InitDouble(paramNames[1], 0.5, 0., 1., 0.01);
  pPlugin->GetParam(params[2])->InitDouble(paramNames[2], 1., 0., 1., 0.01);
  pPlugin->GetParam(params[3])->InitDouble(paramNames[3], 0., 0., 1., 0.01);
  // Reload original values if not resetting
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      pPlugin->GetParam(params[i])->Set(paramVals[i]);
  }

  // Display funcs
  pPlugin->GetParam(params[0])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[2])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[3])->SetDisplayFunc(PercentDisplayFunc);

  for (int i{ 0 }; i < TABLITSA_EFFECT_PARAMS; ++i)
    controls[i]->SetValue(pPlugin->GetParam(params[i])->GetNormalized());

  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Echo");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Damping");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("High-End");
  dynamic_cast<IVKnobControl*>(controls[3])->SetLabelStr("Mix");

  controls[2]->Hide(false);
  controls[3]->Hide(false);

  controls[0]->SetDisabled(false);
  controls[1]->SetDisabled(false);
  controls[2]->SetDisabled(false);
  controls[3]->SetDisabled(false);
  // Toggles
  controls[4]->Hide(true);
  controls[5]->Hide(true);
  controls[4]->SetActionFunction(nullptr);
  controls[5]->SetActionFunction(nullptr);
}

/* SAMPLE AND HOLD (VOICE) */

void InitSampleAndHoldUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset)
{
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  // Save current values
  double paramVals[4]{};
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      paramVals[i] = pPlugin->GetParam(params[i])->Value();
  }
  pPlugin->GetParam(params[0])->InitDouble(paramNames[0], 10., TABLITSA_SAH_MIN_MS, TABLITSA_SAH_MAX_MS, 0.01, "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
  pPlugin->GetParam(params[1])->InitDouble(paramNames[1], 0., 0., 1., 0.01, "", 0, "", IParam::ShapePowCurve(2.));
  pPlugin->GetParam(params[2])->InitDouble(paramNames[2], 0., 0., 1., 0.01);
  pPlugin->GetParam(params[3])->InitDouble(paramNames[3], 0., 0., 1., 0.01);
  // Reload original values if not resetting
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      pPlugin->GetParam(params[i])->Set(paramVals[i]);
  }

  pPlugin->GetParam(params[0])->SetDisplayFunc(nullptr);
  pPlugin->GetParam(params[1])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[2])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[3])->SetDisplayFunc(PercentDisplayFunc);

  // Reset param 1 display texts
  auto* p1 = pPlugin->GetParam(params[0]);
  p1->SetDisplayPrecision(2);

  for (int i{ 0 }; i < TABLITSA_EFFECT_PARAMS; ++i)
    controls[i]->SetValue(pPlugin->GetParam(params[i])->GetNormalized());

  pGraphics->HideControl(params[4], true);
  pGraphics->HideControl(params[5], true);
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Rate");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Decimate");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Noise");
  dynamic_cast<IVKnobControl*>(controls[3])->SetLabelStr("Mix");
  // Modulation ON for knob 1
  dynamic_cast<TablitsaIVModKnobControl*>(controls[0])->EnableModulation(true);
  // Toggle action functions
  controls[4]->SetActionFunction(nullptr);
  controls[5]->SetActionFunction(nullptr);
  controls[4]->Hide(true);
  controls[5]->Hide(true);
  for (auto* knob : allKnobs)
    knob->SetDisabled(false);
}

/* TEXTURIZER (VOICE) */

void InitTexturizerUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset)
{
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  for (auto* knob : allKnobs)
    knob->SetDisabled(false);

  const double minResFreq = 100. / pPlugin->GetSampleRate();
  const double maxResFreq = 10000. / pPlugin->GetSampleRate();

  // Save current values
  double paramVals[4]{};
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      paramVals[i] = pPlugin->GetParam(params[i])->Value();
  }
  pPlugin->GetParam(params[0])->InitDouble(paramNames[0], 0., 0., 1., 0.01);
  pPlugin->GetParam(params[1])->InitDouble(paramNames[1], 0., 0., 1., 0.01);
  pPlugin->GetParam(params[2])->InitDouble(paramNames[2], 0., minResFreq, maxResFreq, 0.01, "", 0, "", IParam::ShapeExp());
  pPlugin->GetParam(params[3])->InitDouble(paramNames[3], 0., 0., 1., 0.01);
  // Reload original values if not resetting
  if (!reset)
  {
    for (int i{ 0 }; i < 4; ++i)
      pPlugin->GetParam(params[i])->Set(paramVals[i]);
  }

  pGraphics->HideControl(params[4], true);
  pGraphics->HideControl(params[5], true);

  pPlugin->GetParam(params[0])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[1])->SetDisplayFunc(PercentDisplayFunc);
  pPlugin->GetParam(params[2])->SetDisplayFunc([pPlugin, minResFreq, maxResFreq](double value, WDL_String& str) {
    str.SetFormatted(MAX_PARAM_DISPLAY_LEN, "%.0f", pPlugin->GetSampleRate() * value);
    str.Append(" Hz");
    });
  pPlugin->GetParam(params[3])->SetDisplayFunc(PercentDisplayFunc);

  for (int i{ 0 }; i < TABLITSA_EFFECT_PARAMS; ++i)
    controls[i]->SetValue(pPlugin->GetParam(params[i])->GetNormalized());

  // All knobs used
  controls[0]->SetDisabled(false);
  controls[1]->SetDisabled(false);
  controls[2]->SetDisabled(false);
  controls[3]->SetDisabled(false);
  // Labels
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Cutoff");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Drive");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Res Freq");
  dynamic_cast<IVKnobControl*>(controls[3])->SetLabelStr("Res Amt");
  // Modulation on for knob 1
  dynamic_cast<TablitsaIVModKnobControl*>(controls[0])->EnableModulation(true);
  // Toggle action functions
  controls[4]->SetActionFunction(nullptr);
  controls[5]->SetActionFunction(nullptr);
  controls[4]->Hide(true);
  controls[5]->Hide(true);
}

/* DEFAULT (BOTH) */

void InitDefaultUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames)
{
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr(" ");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr(" ");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr(" ");
  for (auto ctrl : controls)
    ctrl->SetDisabled(true);
  // Toggle action functions
  for (auto* ctrl : controls)
    ctrl->SetActionFunction(nullptr);
}

void SwapMasterEffectsUI(int effectSlot, IControl* pEffectsList, IGraphics* pGraphics, Plugin* pPlugin, const bool reset)
{
  constexpr int numEffectParams = 6;
  int msg{ kMsgMasterEffect1Changed + effectSlot };
  std::vector<int> params{
    kParamMasterEffect1Param1 + effectSlot * numEffectParams,
    kParamMasterEffect1Param2 + effectSlot * numEffectParams,
    kParamMasterEffect1Param3 + effectSlot * numEffectParams,
    kParamMasterEffect1Param4 + effectSlot * numEffectParams,
    kParamMasterEffect1Param5 + effectSlot * numEffectParams,
    kParamMasterEffect1Param6 + effectSlot * numEffectParams
  };

  // Pointers to the knobs that control the selected effect slot
  IControl* knob1 = pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob1);
  IControl* knob2 = pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob2);
  IControl* knob3 = pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob3);
  IControl* knob4 = pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob4);
  IControl* toggle1 = pGraphics->GetControlWithTag(kCtrlTagMasterEffectsToggle1);
  IControl* toggle2 = pGraphics->GetControlWithTag(kCtrlTagMasterEffectsToggle2);
  std::vector<IControl*> allKnobs{ knob1, knob2, knob3, knob4 };
  std::vector<IControl*> controls{ knob1, knob2, knob3, knob4, toggle1, toggle2 };

  // Swap out effect parameters (and rearrange controls if necessary)
  int effectIdx = dynamic_cast<DropdownListControl*>(pEffectsList)->GetCurrentIndex(); // ID number of the effect
  switch (effectIdx)
  {
  case kDelayEffect:
  {
    std::vector<char*> paramNames{ "Delay L", "Delay R" , "Delay Feedback", "Delay Wet Mix" };
    InitDelayUI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  case kEQEffect:
  {
    std::vector<char*> paramNames{ "EQ Low Gain", "EQ Mid Gain" , "EQ Mid Position", "EQ High Gain" };
    InitEqualizerUI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  case kReverbEffect:
  {
    std::vector<char*> paramNames{ "Reverb Size", "Reverb Gain" , "Master Effect Parameter 3", "Reverb Mix" };
    InitReverbUI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  case kReverb2Effect:
  {
    std::vector<char*> paramNames{ "Size", "Damping" , "Color", "Level" };
    InitReverb2UI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  default:
  {
    std::vector<char*> paramNames{ "MasterEffect 1 Parameter 1", "MasterEffect 1 Parameter 2" , "MasterEffect 1 Parameter 3", "MasterEffect 1 Parameter 4" };
    InitDefaultUI(pPlugin, pGraphics, controls, params, paramNames);
    break;
  }
  }
  if (reset)
    pPlugin->SendArbitraryMsgFromUI(msg, kNoTag, sizeof(effectIdx), reinterpret_cast<void*>(&effectIdx)); // Effects must be swapped before OnParamChange is called
  for (auto* knob : allKnobs)
    knob->SetDirty(true);
  dynamic_cast<Tablitsa*>(pPlugin)->RefreshEffectBankControl();
}

void SwapMasterEffectsUI(IControl* pEffectsList, IGraphics* pGraphics, Plugin* pPlugin, const bool reset)
{
  int effectSlot = static_cast<int>(pGraphics->GetControlWithTag(kCtrlTagMasterEffectsSwitch)->GetValue() * (TABLITSA_MAX_MASTER_EFFECTS - 1)); // ID number of the slot into which the effect will be inserted
  // Save the ID number of newly inserted effect in the pPluginin class effects bank list:
  int effectIdx = dynamic_cast<DropdownListControl*>(pEffectsList)->GetCurrentIndex(); // ID number of the effect
  dynamic_cast<Tablitsa*>(pPlugin)->SetMasterFXSlot(effectSlot, static_cast<EMasterEffectTypes>(effectIdx));
  SwapMasterEffectsUI(effectSlot, pEffectsList, pGraphics, pPlugin, reset);
}

void SwapVoiceEffectsUI(int effectSlot, IControl* pEffectsList, IGraphics* pGraphics, Plugin* pPlugin, const bool reset)
{
  constexpr int numEffectParams = kParamVoiceEffect2Param1 - kParamVoiceEffect1Param1;

  int msg{ kMsgVoiceEffect1Changed + effectSlot };
  std::vector<int> params{
    kParamVoiceEffect1Param1 + effectSlot * numEffectParams,
    kParamVoiceEffect1Param2 + effectSlot * numEffectParams,
    kParamVoiceEffect1Param3 + effectSlot * numEffectParams,
    kParamVoiceEffect1Param4 + effectSlot * numEffectParams,
    kParamVoiceEffect1Param5 + effectSlot * numEffectParams,
    kParamVoiceEffect1Param6 + effectSlot * numEffectParams
  };

  // Pointers to the knobs that control the selected effect slot
  IControl* knob1 = pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob1);
  IControl* knob2 = pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob2);
  IControl* knob3 = pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob3);
  IControl* knob4 = pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob4);
  IControl* toggle1 = pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsToggle1);
  IControl* toggle2 = pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsToggle2);
  std::vector<IControl*> allKnobs{ knob1, knob2, knob3, knob4 };
  std::vector<IControl*> controls{ knob1, knob2, knob3, knob4, toggle1, toggle2 };

  // Swap out effect parameters (and rearrange controls if necessary)
  int effectIdx = dynamic_cast<DropdownListControl*>(pEffectsList)->GetCurrentIndex(); // ID number of the effect
  switch (effectIdx)
  {
  case kDistortionEffect:
  {
    std::vector<char*> paramNames{ "Distortion Type", "Distortion Gain" , "Distortion Noise", "Distortion Wet Mix" };
    InitDistortionUI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  case kSampleAndHoldEffect:
  {
    std::vector<char*> paramNames{ "Sample & Hold Rate", "Sample & Hold Bit Crush" , "Sample & Hold Noise", "Sample & Hold Wet Mix" };
    InitSampleAndHoldUI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  case kTexturizerEffect:
  {
    std::vector<char*> paramNames{ "Texturizer Cutoff", "Texturizer Drive" , "Texturizer Res. Freq.", "Texturizer Res. Amt." };
    InitTexturizerUI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  case kCMEffect:
  {
    std::vector<char*> paramNames{ "SuperRing Depth", "Voice Effect Param 2" , "Voice Effect Param 3", "SuperRing Mix" };
    InitCoefModUI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  default:
  {
    std::vector<char*> paramNames{ "Voice Effect 1 Parameter 1", "Voice Effect 1 Parameter 2" , "Voice Effect 1 Parameter 3", "Voice Effect 1 Parameter 4" };
    InitDefaultUI(pPlugin, pGraphics, controls, params, paramNames);
    break;
  }
  }
  if (reset)
    pPlugin->SendArbitraryMsgFromUI(msg, kNoTag, sizeof(effectIdx), reinterpret_cast<void*>(&effectIdx)); // Effects must be swapped before OnParamChange is called
  for (auto* knob : allKnobs)
    knob->SetDirty(true);
  dynamic_cast<Tablitsa*>(pPlugin)->RefreshEffectBankControl();
}

void SwapVoiceEffectsUI(IControl* pEffectsList, IGraphics* pGraphics, Plugin* pPlugin, const bool reset)
{
  constexpr int numEffectParams = kParamVoiceEffect2Param1 - kParamVoiceEffect1Param1;

  int effectIdx = dynamic_cast<DropdownListControl*>(pEffectsList)->GetCurrentIndex(); // ID number of the effect
  // Save the ID number of newly inserted effect in the pPluginin class effects bank list:
  int effectSlot = static_cast<int>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsSwitch)->GetValue() * (TABLITSA_MAX_VOICE_EFFECTS - 1)); // ID number of the slot into which the effect will be inserted
  dynamic_cast<Tablitsa*>(pPlugin)->SetVoiceFXSlot(effectSlot, static_cast<EVoiceEffectTypes>(effectIdx));
  SwapVoiceEffectsUI(effectSlot, pEffectsList, pGraphics, pPlugin, reset);
}

/* Other functions */

void ResetFilterControls(IGraphics* pGraphics, const int fltFlags)
{
  if (fltFlags & 1)
  {
    // Filter 1
    pGraphics->GetControlWithTag(kCtrlTagFilter1Mode)->SetDisabled(true);
    pGraphics->HideControl(kParamFilter1Cutoff, false);
    pGraphics->HideControl(kParamFilter1FF, true);
    pGraphics->HideControl(kParamFilter1Resonance, false);
    pGraphics->HideControl(kParamFilter1FB, true);
    // Drive/delay switch
    pGraphics->HideControl(kParamFilter1Drive, false);
    pGraphics->HideControl(kParamFilter1Delay, true);
  }
  if (fltFlags & 2)
  {
    // Filter 2
    pGraphics->GetControlWithTag(kCtrlTagFilter2Mode)->SetDisabled(true);
    pGraphics->HideControl(kParamFilter2Cutoff, false);
    pGraphics->HideControl(kParamFilter2FF, true);
    pGraphics->HideControl(kParamFilter2Resonance, false);
    pGraphics->HideControl(kParamFilter2FB, true);
    // Drive/delay switch
    pGraphics->HideControl(kParamFilter2Drive, false);
    pGraphics->HideControl(kParamFilter2Delay, true);
  }
}