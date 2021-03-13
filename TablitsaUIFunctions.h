#include "Tablitsa.h"

#define TABLITSA_MAX_EFFECT_DELAY_MS 5000.
#define TABLITSA_EFFECT_PARAMS 6

void PercentDisplayFunc(double value, WDL_String& str);

void DelayDisplayFunc(double value, WDL_String& str);

void SetAllEffectControlsDirty(IGraphics* pGraphics, int idx);


/* DELAY (MASTER) */

void DelayTempoSyncToggle(Plugin* pPlugin, IGraphics* pGraphics, const std::vector<int>& params, bool isTempoSync, const bool reset = false);

void DelayStereoToggle(IGraphics* pGraphics, const std::vector<IControl*>& controls, const bool mono);

void InitDelayUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset = true);

/* DISTORTION (VOICE) */

void InitDistortionUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset = false);

/* EQ (MASTER) */

void InitEqualizerUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset = true);

/* REVERB 1 & 2 (MASTER) */

void InitReverbUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset = true);

void InitReverb2UI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset = true);

/* SAMPLE AND HOLD (VOICE) */

void InitSampleAndHoldUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset = true);

/* TEXTURIZER (VOICE) */

void InitTexturizerUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset = true);

/* WAVESHAERP (Voice) */

void InitWaveshaperUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames, const bool reset = true);

/* DEFAULT (BOTH) */

void InitDefaultUI(Plugin* pPlugin, IGraphics* pGraphics, std::vector<IControl*> controls, const std::vector<int>& params, const std::vector<char*>& paramNames);

void SwapMasterEffectsUI(int effectSlot, IControl* pEffectsList, IGraphics* pGraphics, Plugin* pPlugin, const bool reset);

void SwapMasterEffectsUI(IControl* pEffectsList, IGraphics* pGraphics, Plugin* pPlugin, const bool reset);

void SwapVoiceEffectsUI(int effectSlot, IControl* pEffectsList, IGraphics* pGraphics, Plugin* pPlugin, const bool reset);

void SwapVoiceEffectsUI(IControl* pEffectsList, IGraphics* pGraphics, Plugin* pPlugin, const bool reset);