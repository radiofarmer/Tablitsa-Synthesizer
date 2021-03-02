void SetAllEffectControlsDirty(IGraphics* pGraphics, int idx)
{
  std::vector<IControl*> controls{
    pGraphics->GetControlWithTag(idx),
    pGraphics->GetControlWithTag(idx + 1),
    pGraphics->GetControlWithTag(idx + 2),
    pGraphics->GetControlWithTag(idx + 3)/*,
    pGraphics->GetControlWithTag(idx + 4),
    pGraphics->GetControlWithTag(idx + 5)*/
  };
  for (auto* ctrl : controls)
    ctrl->SetDirty(true);
}

void InitDelayUI(Plugin* plug, IGraphics* g, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames)
{
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  plug->GetParam(params[0])->InitDouble(paramNames[0], 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
  plug->GetParam(params[1])->InitDouble(paramNames[1], 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
  plug->GetParam(params[2])->InitPercentage(paramNames[2], 0.);
  plug->GetParam(params[3])->InitPercentage(paramNames[3], 0.);
  g->HideControl(params[4], false);
  g->HideControl(params[5], true);
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Left");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Right");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Feedback");
  controls[4]->Hide(false);
  controls[5]->Hide(false);
  controls[4]->SetDirty(true);
  controls[5]->SetDirty(true);
  for (auto* knob : allKnobs)
    knob->SetDisabled(false);
}
void InitWaveshaperUI(Plugin* plug, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames)
{
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  for (auto* knob : allKnobs)
    knob->SetDisabled(false);
  plug->GetParam(params[0])->InitEnum(paramNames[0], EWaveshaperMode::kWaveshapeSine, { WAVESHAPE_TYPES });
  plug->GetParam(params[1])->InitPercentage(paramNames[1], 0.);
  controls[2]->SetDisabled(true);
  plug->GetParam(kParamEffect1Param4)->InitPercentage(paramNames[3], 0.);
  pGraphics->HideControl(params[4], true);
  pGraphics->HideControl(params[5], true);
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Type");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Gain");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("");
  controls[4]->Hide(true);
  controls[5]->Hide(true);
}

void InitSampleAndHoldUI(Plugin* plug, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames)
{
  std::vector<IControl*> allKnobs;
  allKnobs.insert(allKnobs.begin(), controls.begin(), controls.begin() + 4);
  plug->GetParam(params[0])->InitDouble(paramNames[0], 10., 1., 20., 0.1, "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
  plug->GetParam(params[1])->InitDouble(paramNames[1], 0., 1., 100., 1., "%", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
  plug->GetParam(params[2])->InitPercentage(paramNames[2], 0.);
  plug->GetParam(params[3])->InitPercentage(paramNames[3], 0.);
  pGraphics->HideControl(params[4], true);
  pGraphics->HideControl(params[5], true);
  dynamic_cast<IVKnobControl*>(controls[0])->SetLabelStr("Rate");
  dynamic_cast<IVKnobControl*>(controls[1])->SetLabelStr("Decay");
  dynamic_cast<IVKnobControl*>(controls[2])->SetLabelStr("Noise");
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
  controls[0]->SetDisabled(true);
  controls[1]->SetDisabled(true);
  controls[2]->SetDisabled(true);
  controls[3]->SetDisabled(true);
}

void DelayTempoSyncToggle(Plugin* plug, IGraphics* pGraphics, const std::vector<int>& params, bool isTempoSync, int effectListCtrlTag)
{
  int currentEffect = dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(effectListCtrlTag))->GetCurrentIndex();
  if (currentEffect == kDelayEffect)
  {
    if (isTempoSync)
    {
      plug->GetParam(params[0])->InitEnum(plug->GetParam(params[0])->GetName(), LFO<>::k8th, { DELAY_TEMPODIV_VALIST });
      plug->GetParam(params[1])->InitEnum(plug->GetParam(params[1])->GetName(), LFO<>::k8th, { DELAY_TEMPODIV_VALIST });
    }
    else
    {
      plug->GetParam(params[0])->InitDouble(plug->GetParam(params[0])->GetName(), 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
      plug->GetParam(params[1])->InitDouble(plug->GetParam(params[1])->GetName(), 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
    }
    pGraphics->SetAllControlsDirty();
  }
}

void DelayStereoToggle(IGraphics* pGraphics, const std::vector<int>& ctrlTags, const bool mono, IActionFunction channelLockFunction, int effectListCtrlTag)
{
  int currentEffect = dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(effectListCtrlTag))->GetCurrentIndex();
  if (currentEffect == kDelayEffect)
  {
    if (mono)
    {
      pGraphics->GetControlWithTag(ctrlTags[0])->SetActionFunction(channelLockFunction);
      pGraphics->GetControlWithTag(ctrlTags[1])->SetActionFunction(channelLockFunction);
    }
    else
    {
      pGraphics->GetControlWithTag(ctrlTags[0])->SetActionFunction([pGraphics](IControl* pControl) {});
      pGraphics->GetControlWithTag(ctrlTags[1])->SetActionFunction([pGraphics](IControl* pControl) {});
    }
    pGraphics->SetAllControlsDirty();
  }
}