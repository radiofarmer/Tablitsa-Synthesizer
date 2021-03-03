#define TABLITSA_MAX_EFFECT_DELAY_MS 5000.

void PercentDisplayFunc(double value, WDL_String& str)
{
  int precision = 1;
  double displayValue = value;
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

void DelayTempoSyncToggle(Plugin* plug, IGraphics* pGraphics, const std::vector<int>& params, bool isTempoSync, const bool reset = false)
{
  if (isTempoSync)
  {
    const IParam* oldParams[]{ plug->GetParam(params[0]), plug->GetParam(params[1]) };
    plug->GetParam(params[0])->InitEnum(plug->GetParam(params[0])->GetName(), reset ? LFO<>::k8th : oldParams[0]->Value(), { DELAY_TEMPODIV_VALIST });
    plug->GetParam(params[1])->InitEnum(plug->GetParam(params[1])->GetName(), reset ? LFO<>::k8th : oldParams[1]->Value(), { DELAY_TEMPODIV_VALIST });
  }
  else
  {
    const IParam* oldParams[]{ plug->GetParam(params[0]), plug->GetParam(params[1]) };
    plug->GetParam(params[0])->InitDouble(oldParams[0]->GetName(), reset ? 100. : oldParams[0]->Value(), 1., TABLITSA_MAX_EFFECT_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
    plug->GetParam(params[1])->InitDouble(oldParams[1]->GetName(), reset ? 100. : oldParams[1]->Value(), 1., TABLITSA_MAX_EFFECT_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
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

void InitDelayUI(Plugin* plug, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset=true)
{
  // Adjust param ranges
  plug->GetParam(params[0])->InitDouble(paramNames[0], reset ? 100. : plug->GetParam(params[0])->Value(), 1., TABLITSA_MAX_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
  plug->GetParam(params[1])->InitDouble(paramNames[1], reset ? 100. : plug->GetParam(params[1])->Value(), 1., TABLITSA_MAX_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
  plug->GetParam(params[2])->InitPercentage(paramNames[2], reset ? 0. : plug->GetParam(params[2])->Value());
  plug->GetParam(params[3])->InitPercentage(paramNames[3], reset ? 0. : plug->GetParam(params[3])->Value());
  plug->GetParam(params[4])->InitBool(paramNames[2], reset ? false : plug->GetParam(params[4])->Value());
  plug->GetParam(params[5])->InitBool(paramNames[3], reset ? false : plug->GetParam(params[5])->Value());

  pGraphics->HideControl(params[4], false);
  pGraphics->HideControl(params[5], true);
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Left");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Right");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Feedback");

  // Toggle Control Action Functions
  controls[4]->SetActionFunction(
    [pGraphics, plug, params, controls](IControl* pControl) {
      bool delayIsSynced = pControl->GetValue() > 0.5;
      std::vector<int> syncParams{ params[0], params[1] };
      DelayTempoSyncToggle(plug, pGraphics, syncParams, delayIsSynced);
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
  controls[4]->SetDirty(true);
  controls[5]->SetDirty(true);
  // Enable all knobs
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  for (auto* knob : allKnobs)
    knob->SetDisabled(false);
}
void InitWaveshaperUI(Plugin* plug, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset=true)
{
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  for (auto* knob : allKnobs)
    knob->SetDisabled(false);
  plug->GetParam(params[0])->InitEnum(paramNames[0], reset ? EWaveshaperMode::kWaveshapeSine : plug->GetParam(params[0])->Value(), { WAVESHAPE_TYPES });
  plug->GetParam(params[1])->InitPercentage(paramNames[1], reset ? 0. : plug->GetParam(params[1])->Value());
  controls[2]->SetDisabled(true);
  plug->GetParam(params[3])->InitPercentage(paramNames[3], reset ? 0. : plug->GetParam(params[3])->Value());
  pGraphics->HideControl(params[4], true);
  pGraphics->HideControl(params[5], true);
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Type");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Gain");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("");
  // Toggle action functions
  controls[4]->SetActionFunction(nullptr);
  controls[5]->SetActionFunction(nullptr);
  controls[4]->Hide(true);
  controls[5]->Hide(true);
}

void InitSampleAndHoldUI(Plugin* plug, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset=true)
{
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  plug->GetParam(params[0])->InitDouble(paramNames[0], reset ? 10. : plug->GetParam(params[0])->Value(), 1., 20., 0.1, "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
  plug->GetParam(params[1])->InitPercentage(paramNames[1], reset ? 0. : plug->GetParam(params[1])->Value());
  plug->GetParam(params[2])->InitPercentage(paramNames[2], reset ? 10. : plug->GetParam(params[2])->Value());
  plug->GetParam(params[3])->InitPercentage(paramNames[3], reset ? 10. : plug->GetParam(params[3])->Value());
  pGraphics->HideControl(params[4], true);
  pGraphics->HideControl(params[5], true);
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Rate");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Decay");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Noise");
  // Toggle action functions
  controls[4]->SetActionFunction(nullptr);
  controls[5]->SetActionFunction(nullptr);
  controls[4]->Hide(true);
  controls[5]->Hide(true);
  for (auto* knob : allKnobs)
    knob->SetDisabled(false);
}

void InitDefaultUI(Plugin* plug, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames)
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

void SwapMasterEffectsUI(IControl* pEffectsList, IGraphics* pGraphics, Plugin* pPlugin, const bool reset)
{
  constexpr int numEffectParams = 6;

  int effectIdx = dynamic_cast<DropdownListControl*>(pEffectsList)->GetCurrentIndex(); // ID number of the effect
  int effectSlot = static_cast<int>(pGraphics->GetControlWithTag(kCtrlTagMasterEffectsSwitch)->GetValue() * (TABLITSA_MAX_MASTER_EFFECTS - 1)); // ID number of the slot into which the effect will be inserted
  // Save the ID number of newly inserted effect in the plugin class effects bank list:
  dynamic_cast<Tablitsa*>(pPlugin)->SetMasterFXSlot(effectSlot, static_cast<EMasterEffectTypes>(effectIdx));

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
  switch (effectIdx)
  {
  case kDelayEffect:
  {
    std::vector<char*> paramNames{ "Delay L", "Delay R" , "Delay Feedback", "Delay Wet Mix" };
    InitDelayUI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  default:
  {
    std::vector<char*> paramNames{ "MasterEffect 1 Parameter 1", "MasterEffect 1 Parameter 2" , "MasterEffect 1 Parameter 3", "MasterEffect 1 Parameter 4" };
    InitDefaultUI(pPlugin, pGraphics, controls, params, paramNames);
    break;
  }
  }
  pPlugin->SendArbitraryMsgFromUI(msg, kNoTag, sizeof(effectIdx), reinterpret_cast<void*>(&effectIdx)); // Effects must be swapped before OnParamChange is called
  for (auto* knob : allKnobs)
    knob->SetDirty(true);
}

void SwapVoiceEffectsUI(IControl* pEffectsList, IGraphics* pGraphics, Plugin* pPlugin, const bool reset)
{
  constexpr int numEffectParams = kParamVoiceEffect2Param1 - kParamVoiceEffect1Param1;

  int effectIdx = dynamic_cast<DropdownListControl*>(pEffectsList)->GetCurrentIndex(); // ID number of the effect
  int effectSlot = static_cast<int>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsSwitch)->GetValue() * (TABLITSA_MAX_VOICE_EFFECTS - 1)); // ID number of the slot into which the effect will be inserted
  // Save the ID number of newly inserted effect in the plugin class effects bank list:
  dynamic_cast<Tablitsa*>(pPlugin)->SetVoiceFXSlot(effectSlot, static_cast<EVoiceEffectTypes>(effectIdx));

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
  switch (effectIdx)
  {
  case kWaveshaperEffect:
  {
    std::vector<char*> paramNames{ "Waveshaper Mode", "Waveshaper Gain" , "", "Waveshaper Wet Mix" };
    InitWaveshaperUI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  case kSampleAndHoldEffect:
  {
    std::vector<char*> paramNames{ "Sample & Hold Rate", "Sample & Hold Decay" , "Sample & Hold Noise", "Sample & Hold Wet Mix" };
    InitSampleAndHoldUI(pPlugin, pGraphics, controls, params, paramNames, reset);
    break;
  }
  default:
  {
    std::vector<char*> paramNames{ "Voice Effect 1 Parameter 1", "Voice Effect 1 Parameter 2" , "Voice Effect 1 Parameter 3", "Voice Effect 1 Parameter 4" };
    InitDefaultUI(pPlugin, pGraphics, controls, params, paramNames);
    break;
  }
  }
  pPlugin->SendArbitraryMsgFromUI(msg, kNoTag, sizeof(effectIdx), reinterpret_cast<void*>(&effectIdx)); // Effects must be swapped before OnParamChange is called
  for (auto* knob : allKnobs)
    knob->SetDirty(true);
}