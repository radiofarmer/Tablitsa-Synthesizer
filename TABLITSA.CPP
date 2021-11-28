
#include "Tablitsa.h"
#include "IPlug_include_in_plug_src.h"
#include "TablitsaUIFunctions.h"
#include "Presets.h"

#include <atlstr.h>

#define UI_NROWS 4
#define UI_NCOLS 3
#define KEYBOARD
#define FILTER_PARTITION 0.8
#define ENVELOPE_PARTITION 0.8
#define LFO_PARTITION 1.4

// Load/save presets with file dialog or internal controls
#define LOAD_WITH_DIALOG 1

Tablitsa::Tablitsa(const InstanceInfo& info)
: Plugin(info, MakeConfig(kNumParams, kNumPresets))
{
  /*Master Controls*/
#if VST3_API
  GetParam(kParamGain)->InitGain("Gain", 0., -24., 12.);
#else
  GetParam(kParamGain)->InitGain("Gain", -24., -72., 0.);
#endif
  GetParam(kParamNoteGlideTime)->InitDouble("Note Glide Time", 0., 0.0, 3000., 5., "ms", 0, "Portamento", IParam::ShapePowCurve(3.));
  GetParam(kParamNoteGlideRate)->InitDouble("Note Glide Rate", 12., 0.5, 128., 0.1, "St/s", 0, "Portamento", IParam::ShapePowCurve(3.));
  GetParam(kParamPortamentoMode)->InitBool("Equal Glide Time", true);
  GetParam(kParamMonophonic)->InitBool("Monophonic", false);
  GetParam(kParamUnisonVoices)->InitInt("Unison Voices", 1, 1, 8);
  GetParam(kParamUnisonDetune)->InitDouble("Unison Detune", 0., 0., 12., 0.01, "St.");
  GetParam(kParamUnisonChord)->InitEnum("Chord", 0, { UNISON_CHORD_LIST });
  GetParam(kParamStereoSpread)->InitDouble("Stereo Spread", 0., -180., 180., 0.5, "°");
  GetParam(kParamVibratoSpeed)->InitDouble("Vibrato Speed", 5., 0., 10., 0.01, "Hz", 0, "Vibrato");
  GetParam(kParamVibratoDepth)->InitDouble("Vibrato Depth", 0., 0., 50., 0.1, "Cts", 0, "Vibrato", IParam::ShapePowCurve(2.));
  // Pan (modulatable)
  GetParam(kParamPan)->InitDouble("Pan", 0., -90., 90., 0.5, "°");
  GetParam(kParamPanEnv1)->InitDouble("Pan Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPanEnv2)->InitDouble("Pan Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPanAmpEnv)->InitDouble("Pan AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamPanLFO1)->InitDouble("Pan LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPanLFO2)->InitDouble("Pan LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPanSeq)->InitDouble("Pan Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamPanVel)->InitDouble("Pan Velocity Mod", 0., -1., 1., 0.01);
  GetParam(kParamPanKTk)->InitDouble("Pan Ketrack Mod", 0., -1., 1., 0.01);
  GetParam(kParamPanRnd)->InitDouble("Pan Random Mod", 0., -1., 1., 0.01);
  /*Modulators*/
  // Envelope 1
  GetParam(kParamEnv1Attack)->InitDouble("Env1 Attack", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv1Decay)->InitDouble("Env1 Decay", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv1Sustain)->InitDouble("Env1 Sustain", 50., 0., 100., 1, "%", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv1Release)->InitDouble("Env1 Release", 10., 2., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv1Velocity)->InitDouble("Env1 Velocity Mod", 0., 0., 1., 0.01);
  GetParam(kParamEnv1DecayCurve)->InitPercentage("Env1 Decay Curve", 100.);
  GetParam(kParamEnv1ReleaseCurve)->InitPercentage("Env1 Release Curve", 50.);
  // Env1 Meta-Modulation
  GetParam(kParamEnv1SustainEnv1)->InitDouble("Env1 Sustain Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv1SustainEnv2)->InitDouble("Env1 Sustain Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv1SustainAmpEnv)->InitDouble("Env1 Sustain AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv1SustainLFO1)->InitDouble("Env1 Sustain LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv1SustainLFO2)->InitDouble("Env1 Sustain LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv1SustainSeq)->InitDouble("Env1 Sustain Sequencer Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv1SustainVel)->InitDouble("Env1 Sustain Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv1SustainKTk)->InitDouble("Env1 Sustain Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv1SustainRnd)->InitDouble("Env1 Sustain Depth", 0., -1., 1., 0.01);
  // Envelope 2
  GetParam(kParamEnv2Attack)->InitDouble("Env2 Attack", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv2Decay)->InitDouble("Env2 Decay", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv2Sustain)->InitDouble("Env2 Sustain", 50., 0., 100., 1, "%", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv2Release)->InitDouble("Env2 Release", 10., 2., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv2Velocity)->InitDouble("Env2 Velocity Mod", 0., 0., 1., 0.01);
  GetParam(kParamEnv2DecayCurve)->InitPercentage("Env1 Decay Curve", 100.);
  GetParam(kParamEnv2ReleaseCurve)->InitPercentage("Env1 Release Curve", 50.);
  // Env2 Meta-Modulation
  GetParam(kParamEnv2SustainEnv1)->InitDouble("Env2 Sustain Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv2SustainEnv2)->InitDouble("Env2 Sustain Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv2SustainAmpEnv)->InitDouble("Env2 Sustain AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv2SustainLFO1)->InitDouble("Env2 Sustain LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv2SustainLFO2)->InitDouble("Env2 Sustain LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv2SustainSeq)->InitDouble("Env2 Sustain Sequencer Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv2SustainVel)->InitDouble("Env2 Sustain Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv2SustainKTk)->InitDouble("Env2 Sustain Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamEnv2SustainRnd)->InitDouble("Env2 Sustain Depth", 0., -1., 1., 0.01);
  // Amplitude Envelopoe
  GetParam(kParamAmpEnvAttack)->InitDouble("AmpEnv Attack", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamAmpEnvDecay)->InitDouble("AmpEnv Decay", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamAmpEnvSustain)->InitDouble("AmpEnv Sustain", 50., 0., 100., 1, "%", IParam::kFlagsNone, "ADSR");
  GetParam(kParamAmpEnvRelease)->InitDouble("AmpEnv Release", 10., 2., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR");
  GetParam(kParamAmpEnvVelocity)->InitDouble("AmpEnv Velocity Mod", 1., 0., 1., 0.01);
  GetParam(kParamAmpEnvDecayCurve)->InitPercentage("Env1 Decay Curve", 100.);
  GetParam(kParamAmpEnvReleaseCurve)->InitPercentage("Env1 Release Curve", 50.);
  // Amp Envelope Meta-Modulation
  GetParam(kParamAmpEnvSustainEnv1)->InitDouble("AmpEnv Sustain Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamAmpEnvSustainEnv2)->InitDouble("AmpEnv Sustain Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamAmpEnvSustainAmpEnv)->InitDouble("AmpEnv Sustain AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamAmpEnvSustainLFO1)->InitDouble("AmpEnv Sustain LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamAmpEnvSustainLFO2)->InitDouble("AmpEnv Sustain LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamAmpEnvSustainSeq)->InitDouble("AmpEnv Sustain Sequencer Depth", 0., -1., 1., 0.01);
  GetParam(kParamAmpEnvSustainVel)->InitDouble("AmpEnv Sustain Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamAmpEnvSustainKTk)->InitDouble("AmpEnv Sustain Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamAmpEnvSustainRnd)->InitDouble("AmpEnv Sustain Depth", 0., -1., 1., 0.01);
  // LFO1 Meta-Modulation
  GetParam(kParamLFO1RateHzEnv1)->InitDouble("LFO1 Amp Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzEnv2)->InitDouble("LFO1 Amp Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzAmpEnv)->InitDouble("LFO1 Amp AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzLFO1)->InitDouble("LFO1 Amp LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzLFO2)->InitDouble("LFO1 Amp LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzVel)->InitDouble("LFO1 Amp Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzKTk)->InitDouble("LFO1 Amp Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzRnd)->InitDouble("LFO1 Amp Depth", 0., -1., 1., 0.01);
  // LFO 1
  GetParam(kParamLFO1Shape)->InitEnum("LFO1 Shape", LFO<>::kTriangle, {TABLITSA_LFO_SHAPE_VALIST});
  GetParam(kParamLFO1Phase)->InitDouble("LFO1 Phase", 0., 0., 1., 0.01);
  GetParam(kParamLFO1RateHz)->InitFrequency("LFO1 Rate", 1., 0.01, 40.);
  GetParam(kParamLFO1RateTempo)->InitEnum("LFO1 Rate", LFO<>::k1, {LFO_TEMPODIV_VALIST});
  GetParam(kParamLFO1RateMode)->InitBool("LFO1 Sync", false);
  GetParam(kParamLFO1Restart)->InitBool("LFO1 Retrigger", false);
  GetParam(kParamLFO1Amp)->InitDouble("LFO Amp", 1., -1., 1., 0.01, "", IParam::kFlagsNone, "LFO");
  // LFO1 Meta-Modulation
  GetParam(kParamLFO1RateHzEnv1)->InitDouble("LFO1 Rate (Hz) Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzEnv2)->InitDouble("LFO1 Rate (Hz) Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzAmpEnv)->InitDouble("LFO1 Rate (Hz) AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzLFO1)->InitDouble("LFO1 Rate (Hz) LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzLFO2)->InitDouble("LFO1 Rate (Hz) LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzVel)->InitDouble("LFO1 Rate (Hz) Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzKTk)->InitDouble("LFO1 Rate (Hz) Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1RateHzRnd)->InitDouble("LFO1 Rate (Hz) Random Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1AmpEnv1)->InitDouble("LFO1 Amp Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1AmpEnv2)->InitDouble("LFO1 Amp Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1AmpAmpEnv)->InitDouble("LFO1 Amp AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1AmpLFO1)->InitDouble("LFO1 Amp LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1AmpLFO2)->InitDouble("LFO1 Amp LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1AmpVel)->InitDouble("LFO1 Amp Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1AmpKTk)->InitDouble("LFO1 Amp Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO1AmpRnd)->InitDouble("LFO1 Amp Random Depth", 0., -1., 1., 0.01);
  // LFO 2
  GetParam(kParamLFO2Shape)->InitEnum("LFO2 Shape", LFO<>::kTriangle, { TABLITSA_LFO_SHAPE_VALIST });
  GetParam(kParamLFO2Phase)->InitDouble("LFO2 Phase", 0., 0., 1., 0.01);
  GetParam(kParamLFO2RateHz)->InitFrequency("LFO2 Rate", 1., 0.01, 40.);
  GetParam(kParamLFO2RateTempo)->InitEnum("LFO2 Rate", LFO<>::k1, { LFO_TEMPODIV_VALIST });
  GetParam(kParamLFO2RateMode)->InitBool("LFO2 Sync", false);
  GetParam(kParamLFO2Restart)->InitBool("LFO2 Retrigger", false);
  GetParam(kParamLFO2Amp)->InitDouble("LFO Amp", 1., -1., 1., 0.01, "", IParam::kFlagsNone, "LFO");
  // LFO2 Meta-Modulation
  GetParam(kParamLFO2RateHzEnv1)->InitDouble("LFO2 Amp Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2RateHzEnv2)->InitDouble("LFO2 Amp Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2RateHzAmpEnv)->InitDouble("LFO2 Amp AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2RateHzLFO1)->InitDouble("LFO2 Amp LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2RateHzLFO2)->InitDouble("LFO2 Amp LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2RateHzVel)->InitDouble("LFO2 Amp Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2RateHzKTk)->InitDouble("LFO2 Amp Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2RateHzRnd)->InitDouble("LFO2 Amp Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2AmpEnv1)->InitDouble("LFO2 Amp Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2AmpEnv2)->InitDouble("LFO2 Amp Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2AmpAmpEnv)->InitDouble("LFO2 Amp AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2AmpLFO1)->InitDouble("LFO2 Amp LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2AmpLFO2)->InitDouble("LFO2 Amp LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2AmpVel)->InitDouble("LFO2 Amp Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2AmpKTk)->InitDouble("LFO2 Amp Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamLFO2AmpRnd)->InitDouble("LFO2 Amp Random Depth", 0., -1., 1., 0.01);
  // Sequencer
  GetParam(kParamSequencerRateHz)->InitFrequency("Sequencer Rate", 1., 0.01, 40.);
  GetParam(kParamSequencerStepMode)->InitEnum("Sequencer Step Mode", Sequencer<iplug::sample>::EStepMode::kAuto, { "Auto", "Gate", "Trigger" });
  GetParam(kParamSequencerRateTempo)->InitEnum("Sequencer Rate", LFO<>::k4th, { LFO_TEMPODIV_VALIST });
  GetParam(kParamSequencerRateMode)->InitBool("Sequencer Sync", true);
  GetParam(kParamSequencerRestart)->InitBool("Sequencer Retrigger", false);
  GetParam(kParamSequencerSteps)->InitInt("Sequencer Length", 16, 1, 16);
  GetParam(kParamSequencerGlide)->InitPercentage("Sequencer Glide");
//  GetParam(kParamSequencerCurve)->InitEnum()
  GetParam(kParamSequencerAmp)->InitDouble("Seqeuencer Amp", 1., 0., 1., 0.01, "", IParam::kFlagsNone, "Sequencer");
  // Sequencer Meta-Modulation
  GetParam(kParamSequencerRateHzEnv1)->InitDouble("Sequencer Rate (Hz) Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerRateHzEnv2)->InitDouble("Sequencer Rate (Hz) Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerRateHzAmpEnv)->InitDouble("Sequencer Rate (Hz) AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerRateHzLFO1)->InitDouble("Sequencer Rate (Hz) LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerRateHzLFO2)->InitDouble("Sequencer Rate (Hz) LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerRateHzSeq)->InitDouble("Sequencer Rate (Hz) Sequencer Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerRateHzVel)->InitDouble("Sequencer Rate (Hz) Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerRateHzKTk)->InitDouble("Sequencer Rate (Hz) Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerRateHzRnd)->InitDouble("Sequencer Rate (Hz) Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerAmpEnv1)->InitDouble("Sequencer Amp Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerAmpEnv2)->InitDouble("Sequencer Amp Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerAmpAmpEnv)->InitDouble("Sequencer Amp AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerAmpLFO1)->InitDouble("Sequencer Amp LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerAmpLFO2)->InitDouble("Sequencer Amp LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerAmpSeq)->InitDouble("Sequencer Amp Sequencer Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerAmpVel)->InitDouble("Sequencer Amp Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerAmpKTk)->InitDouble("Sequencer Amp Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamSequencerAmpRnd)->InitDouble("Sequencer Amp Random Depth", 0., -1., 1., 0.01);

  /*
  Wavetables
  */
  // Wavetable 1 Pitch
  GetParam(kParamWavetable1Pitch)->InitDouble("WT1 Pitch", 0., -24., 24., 0.1);
  GetParam(kParamWavetable1PitchEnv1)->InitDouble("WT1 Pitch Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PitchEnv2)->InitDouble("WT1 Pitch Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PitchAmpEnv)->InitDouble("WT1 Pitch AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PitchLFO1)->InitDouble("WT1 Pitch LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PitchLFO2)->InitDouble("WT1 Pitch LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PitchSeq)->InitDouble("WT1 Pitch Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PitchVel)->InitDouble("WT1 Pitch Velocity Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PitchKTk)->InitDouble("WT1 Pitch Ketrack Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PitchRnd)->InitDouble("WT1 Pitch Random Mod", 0., -1., 1., 0.01);
  // Wavetable 1 Position
  GetParam(kParamWavetable1Pos)->InitDouble("WT1 Pos", 0., 0., 1., 0.01);
  GetParam(kParamWavetable1PosEnv1)->InitDouble("WT1 Pos Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosEnv2)->InitDouble("WT1 Pos Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosAmpEnv)->InitDouble("WT1 Pos AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosLFO1)->InitDouble("WT1 Pos LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosLFO2)->InitDouble("WT1 Pos L Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosSeq)->InitDouble("WT1 Pos Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosVel)->InitDouble("WT1 Pos Velocity Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosKTk)->InitDouble("WT1 Pos Keytrack Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosRnd)->InitDouble("WT1 Pos Random Mod", 0., -1., 1., 0.01);
  // Wavetable 1 Bend
  GetParam(kParamWavetable1Bend)->InitDouble("WT1 Bend", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1BendEnv1)->InitDouble("WT1 Bend Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1BendEnv2)->InitDouble("WT1 Bend Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1BendAmpEnv)->InitDouble("WT1 Bend AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1BendLFO1)->InitDouble("WT1 Bend LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1BendLFO2)->InitDouble("WT1 Bend LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1BendSeq)->InitDouble("WT1 Bend Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1BendVel)->InitDouble("WT1 Bend Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1BendKTk)->InitDouble("WT1 Bend KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1BendRnd)->InitDouble("WT1 Bend Rnd Mod", 0., -1., 1., 0.01);
  // Wavetable 1 Formant
  GetParam(kParamWavetable1Formant)->InitPercentage("WT1 Formant", 0.);
  GetParam(kParamWavetable1FormantEnv1)->InitDouble("WT1 Formant Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1FormantEnv2)->InitDouble("WT1 Formant Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1FormantAmpEnv)->InitDouble("WT1 Formant FormantEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1FormantLFO1)->InitDouble("WT1 Formant LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1FormantLFO2)->InitDouble("WT1 Formant LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1FormantSeq)->InitDouble("WT1 Formant Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1FormantVel)->InitDouble("WT1 Formant Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1FormantKTk)->InitDouble("WT1 Formant KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1FormantRnd)->InitDouble("WT1 Formant Rnd Mod", 0., -1., 1., 0.01);
  // Wavetable 1 Amp
  GetParam(kParamWavetable1Amp)->InitDouble("WT1 Amp", 0.5, 0., 1., 0.01);
  GetParam(kParamWavetable1AmpEnv1)->InitDouble("WT1 Amp Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1AmpEnv2)->InitDouble("WT1 Amp Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1AmpAmpEnv)->InitDouble("WT1 Amp AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1AmpLFO1)->InitDouble("WT1 Amp LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1AmpLFO2)->InitDouble("WT1 Amp LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1AmpSeq)->InitDouble("WT1 Amp Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1AmpVel)->InitDouble("WT1 Amp Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1AmpKTk)->InitDouble("WT1 Amp KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1AmpRnd)->InitDouble("WT1 Amp Rnd Mod", 0., -1., 1., 0.01);

  // Wavetable 2 Pitch
  GetParam(kParamWavetable2Pitch)->InitDouble("WT2 Pitch", 0., -24., 24., 0.1);
  GetParam(kParamWavetable2PitchEnv1)->InitDouble("WT2 Pitch Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PitchEnv2)->InitDouble("WT2 Pitch Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PitchAmpEnv)->InitDouble("WT2 Pitch AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PitchLFO1)->InitDouble("WT2 Pitch LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PitchLFO2)->InitDouble("WT2 Pitch LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PitchSeq)->InitDouble("WT2 Pitch Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PitchVel)->InitDouble("WT2 Pitch Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PitchKTk)->InitDouble("WT2 Pitch LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PitchRnd)->InitDouble("WT2 Pitch Seq Mod", 0., -1., 1., 0.01);
  // Wavetable 2 Position
  GetParam(kParamWavetable2Pos)->InitDouble("WT2 Pos", 0., 0., 1., 0.01);
  GetParam(kParamWavetable2PosEnv1)->InitDouble("WT2 Pos Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PosEnv2)->InitDouble("WT2 Pos Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PosAmpEnv)->InitDouble("WT2 Pos AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PosLFO1)->InitDouble("WT2 Pos LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PosLFO2)->InitDouble("WT2 Pos L Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PosSeq)->InitDouble("WT2 Pos Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PosVel)->InitDouble("WT2 Pos Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PosKTk)->InitDouble("WT2 Pos KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2PosRnd)->InitDouble("WT2 Pos Rnd Mod", 0., -1., 1., 0.01);
  // Wavetable 2 Bend
  GetParam(kParamWavetable2Bend)->InitDouble("WT2 Bend", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2BendEnv1)->InitDouble("WT2 Bend Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2BendEnv2)->InitDouble("WT2 Bend Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2BendAmpEnv)->InitDouble("WT2 Bend AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2BendLFO1)->InitDouble("WT2 Bend LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2BendLFO2)->InitDouble("WT2 Bend LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2BendSeq)->InitDouble("WT2 Bend Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2BendVel)->InitDouble("WT2 Bend Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2BendKTk)->InitDouble("WT2 Bend KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2BendRnd)->InitDouble("WT2 Bend Rnd Mod", 0., -1., 1., 0.01);
  // Wavetable 2 Formant
  GetParam(kParamWavetable2Formant)->InitPercentage("WT2 Formant", 0.);
  GetParam(kParamWavetable2FormantEnv1)->InitDouble("WT2 Formant Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2FormantEnv2)->InitDouble("WT2 Formant Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2FormantAmpEnv)->InitDouble("WT2 Formant FormantEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2FormantLFO1)->InitDouble("WT2 Formant LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2FormantLFO2)->InitDouble("WT2 Formant LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2FormantSeq)->InitDouble("WT2 Formant Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2FormantVel)->InitDouble("WT2 Formant Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2FormantKTk)->InitDouble("WT2 Formant KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2FormantRnd)->InitDouble("WT2 Formant Rnd Mod", 0., -1., 1., 0.01);
  // Wavetable 2 Amp
  GetParam(kParamWavetable2Amp)->InitDouble("WT2 Amp", 0.5, 0., 1., 0.01);
  GetParam(kParamWavetable2AmpEnv1)->InitDouble("WT2 Amp Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2AmpEnv2)->InitDouble("WT2 Amp Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2AmpAmpEnv)->InitDouble("WT2 Amp AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2AmpLFO1)->InitDouble("WT2 Amp LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2AmpLFO2)->InitDouble("WT2 Amp LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2AmpSeq)->InitDouble("WT2 Amp Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2AmpVel)->InitDouble("WT2 Amp Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2AmpKTk)->InitDouble("WT2 Amp KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2AmpRnd)->InitDouble("WT2 Amp Rnd Mod", 0., -1., 1., 0.01);

  // Modulation Oscillator
  GetParam(kParamOscModulator)->InitEnum("Oscillator Modulation Mode", 0, { "Phase", "Ring" });
  GetParam(kParamOsc1PM)->InitBool("Osc1 Phase Modulation", false);
  GetParam(kParamOsc1RM)->InitBool("Osc1 Ring Modulation", false);
  GetParam(kParamOsc2PM)->InitBool("Osc2 Phase Modulation", false);
  GetParam(kParamOsc2RM)->InitBool("Osc2 Ring Modulation", false);
  // Phase Mod Freq
  GetParam(kParamPhaseModFreq)->InitDouble("Phase Modulator Freq", 0., -24., 24., 0.01);
  GetParam(kParamPhaseModFreqEnv1)->InitDouble("Phase Mod Frequency Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModFreqEnv2)->InitDouble("Phase Mod Frequency Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModFreqAmpEnv)->InitDouble("Phase Mod Frequency AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModFreqLFO1)->InitDouble("Phase Mod Frequency LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModFreqLFO2)->InitDouble("Phase Mod Frequency LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModFreqSeq)->InitDouble("Phase Mod Frequency Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModFreqVel)->InitDouble("Phase Mod Frequency Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModFreqKTk)->InitDouble("Phase Mod Frequency KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModFreqRnd)->InitDouble("Phase Mod Frequency Rnd Mod", 0., -1., 1., 0.01);
  // Phase Mod Amount
  GetParam(kParamPhaseModAmount)->InitPercentage("Phase Modulator Amount", 25., 0., 100.);
  GetParam(kParamPhaseModAmountEnv1)->InitDouble("Phase Mod Depth Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModAmountEnv2)->InitDouble("Phase Mod Depth Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModAmountAmpEnv)->InitDouble("Phase Mod Depth AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModAmountLFO1)->InitDouble("Phase Mod Depth LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModAmountLFO2)->InitDouble("Phase Mod Depth LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModAmountSeq)->InitDouble("Phase Mod Depth Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModAmountVel)->InitDouble("Phase Mod Depth Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModAmountKTk)->InitDouble("Phase Mod Depth KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamPhaseModAmountRnd)->InitDouble("Phase Mod Depth Rnd Mod", 0., -1., 1., 0.01);
  // Ring Mod Freq
  GetParam(kParamRingModFreq)->InitDouble("Ring Modulator Freq", 0., -24., 24., 0.01);
  GetParam(kParamRingModFreqEnv1)->InitDouble("Ring Mod Frequency Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModFreqEnv2)->InitDouble("Ring Mod Frequency Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModFreqAmpEnv)->InitDouble("Ring Mod Frequency AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModFreqLFO1)->InitDouble("Ring Mod Frequency LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModFreqLFO2)->InitDouble("Ring Mod Frequency LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModFreqSeq)->InitDouble("Ring Mod Frequency Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModFreqVel)->InitDouble("Ring Mod Frequency Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModFreqKTk)->InitDouble("Ring Mod Frequency KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModFreqRnd)->InitDouble("Ring Mod Frequency Rnd Mod", 0., -1., 1., 0.01);
  // Ring Mod Amount
  GetParam(kParamRingModAmount)->InitPercentage("Ring Modulator Amount", 25., 0., 100.);
  GetParam(kParamRingModAmount)->InitPercentage("Ring Modulator Amount", 25., 0., 100.);
  GetParam(kParamRingModAmountEnv1)->InitDouble("Ring Mod Depth Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModAmountEnv2)->InitDouble("Ring Mod Depth Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModAmountAmpEnv)->InitDouble("Ring Mod Depth AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModAmountLFO1)->InitDouble("Ring Mod Depth LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModAmountLFO2)->InitDouble("Ring Mod Depth LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModAmountSeq)->InitDouble("Ring Mod Depth Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModAmountVel)->InitDouble("Ring Mod Depth Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModAmountKTk)->InitDouble("Ring Mod Depth KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamRingModAmountRnd)->InitDouble("Ring Mod Depth Rnd Mod", 0., -1., 1., 0.01);

  // Filter 1
  GetParam(kParamFilter1Type)->InitEnum("Filter 1 Type", kNoFilter, { FILTER_TYPE_LIST });
  GetParam(kParamFilter1ModeVSF)->InitEnum("Filter 1 Type", kNoFilter, { FILTER_MODE_LIST_VSF });
  GetParam(kParamFilter1ModeMoog)->InitEnum("Filter 1 Type", kNoFilter, { FILTER_MODE_LIST_MOOG });
  GetParam(kParamFilter1ModeComb)->InitEnum("Filter 1 Type", kNoFilter, { FILTER_MODE_LIST_COMB });
  GetParam(kParamFilter1Cutoff)->InitFrequency("Filter 1 Cutoff", 20000., 20., 20000., 1.);
  GetParam(kParamFilter1CutoffEnv1)->InitDouble("Filter 1 Cutoff Env1 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter1CutoffEnv2)->InitDouble("Filter 1 Cutoff Env2 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter1CutoffAmpEnv)->InitDouble("Filter 1 Cutoff AmpEnv Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter1CutoffLFO1)->InitDouble("Filter 1 Cutoff LFO1 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter1CutoffLFO2)->InitDouble("Filter 1 Cutoff LFO2 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter1CutoffSeq)->InitDouble("Filter 1 Cutoff Seq Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter1CutoffVel)->InitDouble("Filter 1 Cutoff LFO1 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter1CutoffKTk)->InitDouble("Filter 1 Cutoff LFO2 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter1CutoffRnd)->InitDouble("Filter 1 Cutoff Seq Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter1Resonance)->InitPercentage("Filter 1 Resonance", 0.1, 0., 99., 1);
  GetParam(kParamFilter1ResonanceEnv1)->InitDouble("Filter 1 Resonance Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1ResonanceEnv2)->InitDouble("Filter 1 Resonance Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1ResonanceAmpEnv)->InitDouble("Filter 1 Resonance AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1ResonanceLFO1)->InitDouble("Filter 1 Resonance LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1ResonanceLFO2)->InitDouble("Filter 1 Resonance LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1ResonanceSeq)->InitDouble("Filter 1 Resonance Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1ResonanceVel)->InitDouble("Filter 1 Resonance LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1ResonanceKTk)->InitDouble("Filter 1 Resonance LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1ResonanceRnd)->InitDouble("Filter 1 Resonance Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1Drive)->InitPercentage("Filter 1 Drive", 0.1, 0., 100., 1);
  GetParam(kParamFilter1DriveEnv1)->InitDouble("Filter 1 Drive Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1DriveEnv2)->InitDouble("Filter 1 Drive Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1DriveAmpEnv)->InitDouble("Filter 1 Drive AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1DriveLFO1)->InitDouble("Filter 1 Drive LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1DriveLFO2)->InitDouble("Filter 1 Drive LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1DriveSeq)->InitDouble("Filter 1 Drive Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1DriveVel)->InitDouble("Filter 1 Drive LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1DriveKTk)->InitDouble("Filter 1 Drive LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter1DriveRnd)->InitDouble("Filter 1 Drive Seq Depth", 0., -1., 1., 0.01);
  // Comb-specific parameters with different ranges (the modulation parameters are common between all filter types, as these always between -1 and 1)
  GetParam(kParamFilter1FF)->InitPercentage("Filter 1 FF");
  GetParam(kParamFilter1FB)->InitPercentage("Filter 1 FB");
  GetParam(kParamFilter1Delay)->InitDouble("Filter 1 Delay", 0.1, 0.1, 20., 0.01, "ms");
  // Oscillator sends
  GetParam(kParamFilter1Osc1Send)->InitGain("Filter 1 Osc 1 Send", 0., SEND_DB_FLOOR, 0.);
  GetParam(kParamFilter1Osc2Send)->InitGain("Filter 1 Osc 2 Send", SEND_DB_FLOOR, SEND_DB_FLOOR, 0.);

  // Filter 2
  GetParam(kParamFilter2Type)->InitEnum("Filter 2 Type", kNoFilter, { FILTER_TYPE_LIST });
  GetParam(kParamFilter2ModeVSF)->InitEnum("Filter 2 Type", kNoFilter, { FILTER_MODE_LIST_VSF });
  GetParam(kParamFilter2ModeMoog)->InitEnum("Filter 2 Type", kNoFilter, { FILTER_MODE_LIST_MOOG });
  GetParam(kParamFilter2ModeComb)->InitEnum("Filter 2 Type", kNoFilter, { FILTER_MODE_LIST_COMB });
  GetParam(kParamFilter2Cutoff)->InitFrequency("Filter 2 Cutoff", 20000., 20., 20000., 1.);
  GetParam(kParamFilter2CutoffEnv1)->InitDouble("Filter 2 Cutoff Env1 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter2CutoffEnv2)->InitDouble("Filter 2 Cutoff Env2 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter2CutoffAmpEnv)->InitDouble("Filter 2 Cutoff AmpEnv Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter2CutoffLFO1)->InitDouble("Filter 2 Cutoff LFO1 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter2CutoffLFO2)->InitDouble("Filter 2 Cutoff LFO2 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter2CutoffSeq)->InitDouble("Filter 2 Cutoff Seq Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter2CutoffVel)->InitDouble("Filter 2 Cutoff LFO1 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter2CutoffKTk)->InitDouble("Filter 2 Cutoff LFO2 Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter2CutoffRnd)->InitDouble("Filter 2 Cutoff Seq Depth", 0., -1., 1., 0.01, "", IParam::kFlagsNone, "Mod Depth", ModShapePowCurve(3.));
  GetParam(kParamFilter2Resonance)->InitPercentage("Filter 2 Resonance", 0.1, 0., 99., 1);
  GetParam(kParamFilter2ResonanceEnv1)->InitDouble("Filter 2 Resonance Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2ResonanceEnv2)->InitDouble("Filter 2 Resonance Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2ResonanceAmpEnv)->InitDouble("Filter 2 Resonance AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2ResonanceLFO1)->InitDouble("Filter 2 Resonance LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2ResonanceLFO2)->InitDouble("Filter 2 Resonance LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2ResonanceSeq)->InitDouble("Filter 2 Resonance Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2ResonanceVel)->InitDouble("Filter 2 Resonance LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2ResonanceKTk)->InitDouble("Filter 2 Resonance LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2ResonanceRnd)->InitDouble("Filter 2 Resonance Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2Drive)->InitPercentage("Filter 2 Drive", 0.1, 0., 100., 1);
  GetParam(kParamFilter2DriveEnv1)->InitDouble("Filter 2 Drive Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2DriveEnv2)->InitDouble("Filter 2 Drive Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2DriveAmpEnv)->InitDouble("Filter 2 Drive AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2DriveLFO1)->InitDouble("Filter 2 Drive LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2DriveLFO2)->InitDouble("Filter 2 Drive LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2DriveSeq)->InitDouble("Filter 2 Drive Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2DriveVel)->InitDouble("Filter 2 Drive Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2DriveKTk)->InitDouble("Filter 2 Drive Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2DriveRnd)->InitDouble("Filter 2 Drive Random Depth", 0., -1., 1., 0.01);
  // Comb-specific parameters with different ranges (the modulation parameters are common between all filter types, as these always between -1 and 1)
  GetParam(kParamFilter2FF)->InitPercentage("Filter 2 FF");
  GetParam(kParamFilter2FB)->InitPercentage("Filter 2 FB");
  GetParam(kParamFilter2Delay)->InitDouble("Filter 2 Delay", 0.1, 0.1, 20., 0.01, "ms");
  // Oscillator sends
  GetParam(kParamFilter2Osc1Send)->InitGain("Filter 2 Osc 1 Send", SEND_DB_FLOOR, SEND_DB_FLOOR, 0.);
  GetParam(kParamFilter2Osc2Send)->InitGain("Filter 2 Osc 2 Send", 0., SEND_DB_FLOOR, 0.);

  // Effect Parameters
  // Voice Effects
  GetParam(kParamVoiceEffect1Param1)->InitDouble("Voice Effect 1 Parameter 1", 0., 0., 1., 0.01);
  GetParam(kParamVoiceEffect1Param1Env1)->InitDouble("Voice Effect 1 Param 1 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param1Env2)->InitDouble("Voice Effect 1 Param 1 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param1AmpEnv)->InitDouble("Voice Effect 1 Param 1 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param1LFO1)->InitDouble("Voice Effect 1 Param 1 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param1LFO2)->InitDouble("Voice Effect 1 Param 1 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param1Seq)->InitDouble("Voice Effect 1 Param 1 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param1Vel)->InitDouble("Voice Effect 1 Param 1 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param1KTk)->InitDouble("Voice Effect 1 Param 1 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param1Rnd)->InitDouble("Voice Effect 1 Param 1 Random Depth", 0., -1., 1., 0.01);
  // P2
  GetParam(kParamVoiceEffect1Param2)->InitDouble("VoiceEffect 1 Parameter 2", 0., 0., 1., 0.01);
  GetParam(kParamVoiceEffect1Param2Env1)->InitDouble("Voice Effect 1 Param 2 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param2Env2)->InitDouble("Voice Effect 1 Param 2 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param2AmpEnv)->InitDouble("Voice Effect 1 Param 2 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param2LFO1)->InitDouble("Voice Effect 1 Param 2 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param2LFO2)->InitDouble("Voice Effect 1 Param 2 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param2Seq)->InitDouble("Voice Effect 1 Param 2 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param2Vel)->InitDouble("Voice Effect 1 Param 2 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param2KTk)->InitDouble("Voice Effect 1 Param 2 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param2Rnd)->InitDouble("Voice Effect 1 Param 2 Random Depth", 0., -1., 1., 0.01);
  // P3
  GetParam(kParamVoiceEffect1Param3)->InitPercentage("VoiceEffect 1 Parameter 3", 0.);
  GetParam(kParamVoiceEffect1Param3Env1)->InitDouble("Voice Effect 1 Param 3 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param3Env2)->InitDouble("Voice Effect 1 Param 3 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param3AmpEnv)->InitDouble("Voice Effect 1 Param 3 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param3LFO1)->InitDouble("Voice Effect 1 Param 3 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param3LFO2)->InitDouble("Voice Effect 1 Param 3 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param3Seq)->InitDouble("Voice Effect 1 Param 3 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param3Vel)->InitDouble("Voice Effect 1 Param 3 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param3KTk)->InitDouble("Voice Effect 1 Param 3 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param3Rnd)->InitDouble("Voice Effect 1 Param 3 Random Depth", 0., -1., 1., 0.01);
  // P4
  GetParam(kParamVoiceEffect1Param4)->InitPercentage("VoiceEffect 1 Parameter 4", 0.);
  GetParam(kParamVoiceEffect1Param4Env1)->InitDouble("Voice Effect 1 Param 4 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param4Env2)->InitDouble("Voice Effect 1 Param 4 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param4AmpEnv)->InitDouble("Voice Effect 1 Param 4 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param4LFO1)->InitDouble("Voice Effect 1 Param 4 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param4LFO2)->InitDouble("Voice Effect 1 Param 4 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param4Seq)->InitDouble("Voice Effect 1 Param 4 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param4Vel)->InitDouble("Voice Effect 1 Param 4 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param4KTk)->InitDouble("Voice Effect 1 Param 4 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param4Rnd)->InitDouble("Voice Effect 1 Param 4 Random Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect1Param5)->InitBool("VoiceEffect 1 Parameter 5", false);
  GetParam(kParamVoiceEffect1Param6)->InitBool("VoiceEffect 1 Parameter 6", false);

  // P1
  GetParam(kParamVoiceEffect2Param1)->InitDouble("Voice Effect 2 Parameter 1", 0., 0., 1., 0.01);
  GetParam(kParamVoiceEffect2Param1Env1)->InitDouble("Voice Effect 2 Param 1 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param1Env2)->InitDouble("Voice Effect 2 Param 1 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param1AmpEnv)->InitDouble("Voice Effect 2 Param 1 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param1LFO1)->InitDouble("Voice Effect 2 Param 1 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param1LFO2)->InitDouble("Voice Effect 2 Param 1 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param1Seq)->InitDouble("Voice Effect 2 Param 1 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param1Vel)->InitDouble("Voice Effect 2 Param 1 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param1KTk)->InitDouble("Voice Effect 2 Param 1 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param1Rnd)->InitDouble("Voice Effect 2 Param 1 Random Depth", 0., -1., 1., 0.01);
  // P2
  GetParam(kParamVoiceEffect2Param2)->InitDouble("Voice Effect 2 Parameter 2", 0., 0., 1., 0.01);
  GetParam(kParamVoiceEffect2Param2Env1)->InitDouble("Voice Effect 2 Param 2 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param2Env2)->InitDouble("Voice Effect 2 Param 2 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param2AmpEnv)->InitDouble("Voice Effect 2 Param 2 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param2LFO1)->InitDouble("Voice Effect 2 Param 2 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param2LFO2)->InitDouble("Voice Effect 2 Param 2 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param2Seq)->InitDouble("Voice Effect 2 Param 2 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param2Vel)->InitDouble("Voice Effect 2 Param 2 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param2KTk)->InitDouble("Voice Effect 2 Param 2 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param2Rnd)->InitDouble("Voice Effect 2 Param 2 Random Depth", 0., -1., 1., 0.01);
  // P3
  GetParam(kParamVoiceEffect2Param3)->InitPercentage("VoiceEffect 2 Parameter 3", 0.);
  GetParam(kParamVoiceEffect2Param3Env1)->InitDouble("Voice Effect 2 Param 3 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param3Env2)->InitDouble("Voice Effect 2 Param 3 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param3AmpEnv)->InitDouble("Voice Effect 2 Param 3 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param3LFO1)->InitDouble("Voice Effect 2 Param 3 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param3LFO2)->InitDouble("Voice Effect 2 Param 3 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param3Seq)->InitDouble("Voice Effect 2 Param 3 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param3Vel)->InitDouble("Voice Effect 2 Param 3 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param3KTk)->InitDouble("Voice Effect 2 Param 3 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param3Rnd)->InitDouble("Voice Effect 2 Param 3 Random Depth", 0., -1., 1., 0.01);
  // P4
  GetParam(kParamVoiceEffect2Param4)->InitPercentage("Voice Effect 2 Parameter 3", 0.);
  GetParam(kParamVoiceEffect2Param4Env1)->InitDouble("Voice Effect 2 Param 4 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param4Env2)->InitDouble("Voice Effect 2 Param 4 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param4AmpEnv)->InitDouble("Voice Effect 2 Param 4 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param4LFO1)->InitDouble("Voice Effect 2 Param 4 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param4LFO2)->InitDouble("Voice Effect 2 Param 4 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param4Seq)->InitDouble("Voice Effect 2 Param 4 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param4Vel)->InitDouble("Voice Effect 2 Param 4 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param4KTk)->InitDouble("Voice Effect 2 Param 4 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect2Param4Rnd)->InitDouble("Voice Effect 2 Param 4 Random Depth", 0., -1., 1., 0.01);
  // P5&6
  GetParam(kParamVoiceEffect2Param5)->InitBool("Voice Effect 2 Parameter 5", false);
  GetParam(kParamVoiceEffect2Param6)->InitBool("Voice Effect 2 Parameter 6", false);

  GetParam(kParamVoiceEffect3Param1)->InitDouble("Voice Effect 3 Parameter 1", 0., 0., 1., 0.01);
  GetParam(kParamVoiceEffect3Param1Env1)->InitDouble("Voice Effect 3 Param 1 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param1Env2)->InitDouble("Voice Effect 3 Param 1 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param1AmpEnv)->InitDouble("Voice Effect 3 Param 1 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param1LFO1)->InitDouble("Voice Effect 3 Param 1 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param1LFO2)->InitDouble("Voice Effect 3 Param 1 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param1Seq)->InitDouble("Voice Effect 3 Param 1 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param1Vel)->InitDouble("Voice Effect 3 Param 1 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param1KTk)->InitDouble("Voice Effect 3 Param 1 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param1Rnd)->InitDouble("Voice Effect 3 Param 1 Random Depth", 0., -1., 1., 0.01);
  // P2
  GetParam(kParamVoiceEffect3Param2)->InitDouble("VoiceEffect 3 Parameter 2", 0., 0., 1., 0.01);
  GetParam(kParamVoiceEffect3Param2Env1)->InitDouble("Voice Effect 3 Param 2 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param2Env2)->InitDouble("Voice Effect 3 Param 2 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param2AmpEnv)->InitDouble("Voice Effect 3 Param 2 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param2LFO1)->InitDouble("Voice Effect 3 Param 2 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param2LFO2)->InitDouble("Voice Effect 3 Param 2 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param2Seq)->InitDouble("Voice Effect 3 Param 2 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param2Vel)->InitDouble("Voice Effect 3 Param 2 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param2KTk)->InitDouble("Voice Effect 3 Param 2 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param2Rnd)->InitDouble("Voice Effect 3 Param 2 Random Depth", 0., -1., 1., 0.01);
  // P3
  GetParam(kParamVoiceEffect3Param3)->InitPercentage("Voice Effect 3 Parameter 3", 0.);
  GetParam(kParamVoiceEffect3Param3Env1)->InitDouble("Voice Effect 3 Param 3 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param3Env2)->InitDouble("Voice Effect 3 Param 3 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param3AmpEnv)->InitDouble("Voice Effect 3 Param 3 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param3LFO1)->InitDouble("Voice Effect 3 Param 3 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param3LFO2)->InitDouble("Voice Effect 3 Param 3 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param3Seq)->InitDouble("Voice Effect 3 Param 3 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param3Vel)->InitDouble("Voice Effect 3 Param 3 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param3KTk)->InitDouble("Voice Effect 3 Param 3 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param3Rnd)->InitDouble("Voice Effect 3 Param 3 Random Depth", 0., -1., 1., 0.01);
  // P4
  GetParam(kParamVoiceEffect3Param4)->InitPercentage("Voice Effect 3 Parameter 4", 0.);
  GetParam(kParamVoiceEffect3Param4Env1)->InitDouble("Voice Effect 3 Param 4 Env1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param4Env2)->InitDouble("Voice Effect 3 Param 4 Env2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param4AmpEnv)->InitDouble("Voice Effect 3 Param 4 AmpEnv Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param4LFO1)->InitDouble("Voice Effect 3 Param 4 LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param4LFO2)->InitDouble("Voice Effect 3 Param 4 LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param4Seq)->InitDouble("Voice Effect 3 Param 4 Seq Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param4Vel)->InitDouble("Voice Effect 3 Param 4 Velocity Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param4KTk)->InitDouble("Voice Effect 3 Param 4 Keytrack Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param4Rnd)->InitDouble("Voice Effect 3 Param 4 Random Depth", 0., -1., 1., 0.01);
  GetParam(kParamVoiceEffect3Param5)->InitBool("VoiceEffect 3 Parameter 5", false);
  GetParam(kParamVoiceEffect3Param6)->InitBool("VoiceEffect 3 Parameter 6", false);

  // Master Effects
  GetParam(kParamMasterEffect1Param1)->InitDouble("Master Effect 1 Parameter 1", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect1Param2)->InitDouble("Master Effect 1 Parameter 2", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect1Param3)->InitPercentage("Master Effect 1 Parameter 3", 0.);
  GetParam(kParamMasterEffect1Param4)->InitPercentage("Master Effect 1 Parameter 4", 0.);
  GetParam(kParamMasterEffect1Param5)->InitBool("Master Effect 1 Parameter 5", false);
  GetParam(kParamMasterEffect1Param6)->InitBool("Master Effect 1 Parameter 6", false);

  GetParam(kParamMasterEffect2Param1)->InitDouble("Master Effect 2 Parameter 1", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect2Param2)->InitDouble("Master Effect 2 Parameter 2", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect2Param3)->InitPercentage("Master Effect 2 Parameter 3", 0.);
  GetParam(kParamMasterEffect2Param4)->InitPercentage("Master Effect 2 Parameter 4", 0.);
  GetParam(kParamMasterEffect2Param5)->InitBool("Master Effect 2 Parameter 5", false);
  GetParam(kParamMasterEffect2Param6)->InitBool("Master Effect 2 Parameter 6", false);

  GetParam(kParamMasterEffect3Param1)->InitDouble("Master Effect 3 Parameter 1", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect3Param2)->InitDouble("Master Effect 3 Parameter 2", 0., 0., 100., 0.1);
  GetParam(kParamMasterEffect3Param3)->InitPercentage("Master Effect 3 Parameter 3", 0.);
  GetParam(kParamMasterEffect3Param4)->InitPercentage("Master Effect 3 Parameter 4", 0.);
  GetParam(kParamMasterEffect3Param5)->InitBool("Master Effect 3 Parameter 5", false);
  GetParam(kParamMasterEffect3Param6)->InitBool("Master Effect 3 Parameter 6", false);

  GetParam(kParamOsc1FilterBypass)->InitDouble("Osc 1 Filter Bypass", 0., 0., 1., 0.01, "", 0, "Bypass", IParam::ShapePowCurve(0.5));
  GetParam(kParamOsc2FilterBypass)->InitDouble("Osc 2 Filter Bypass", 0., 0., 1., 0.01, "", 0, "Bypass", IParam::ShapePowCurve(0.5));
  GetParam(kParamOsc1EffectBypass)->InitDouble("Osc 1 Effect Bypass", 0., 0., 1., 0.01, "", 0, "Bypass", IParam::ShapePowCurve(0.5));
  GetParam(kParamOsc2EffectBypass)->InitDouble("Osc 2 Effect Bypass", 0., 0., 1., 0.01, "", 0, "Bypass", IParam::ShapePowCurve(0.5));

#if VST3_API
#define KEYBOARD_HEIGHT 0.f
#else
#define KEYBOARD_HEIGHT 150.f
#endif
#if IPLUG_EDITOR // http://bit.ly/2S64BDd
  mMakeGraphicsFunc = [&]() {
    return MakeGraphics(*this, PLUG_WIDTH, (int)(PLUG_HEIGHT + KEYBOARD_HEIGHT), PLUG_FPS, GetScaleForScreen((int)(PLUG_WIDTH), (int)(PLUG_HEIGHT + KEYBOARD_HEIGHT)));
  };
  mLayoutFunc = [&](IGraphics* pGraphics) {

    pGraphics->AttachCornerResizer(EUIResizerMode::Scale, false);
    pGraphics->AttachBackground(BG_FN);
    pGraphics->EnableMouseOver(true);
    pGraphics->EnableMultiTouch(true);

#ifdef OS_WEB
    pGraphics->AttachPopupMenuControl();
#endif

    //    pGraphics->EnableLiveEdit(true);
    pGraphics->LoadFont("Roboto-Regular", ROBOTO_FN);
    constexpr float plugBorders = 10.f;
    const IRECT b = pGraphics->GetBounds().GetPadded(-plugBorders);
    constexpr float dialSize{ 80.f };
    IRECT globalControls = b.GetOffset(0, 0, 0, -KEYBOARD_HEIGHT).GetHPadded(-5.f).GetReducedFromTop(15.f);
    const IRECT appControls = b.GetFromTop(30.f);
#ifndef VST3_API
    IRECT keyboardBounds = b.GetFromBottom((float)KEYBOARD_HEIGHT - plugBorders).GetPadded(-5.f, -25.f, -5.f, -2.f);
    IRECT wheelsBounds = keyboardBounds.ReduceFromLeft(100.f).GetPadded(-10.f);

    // Keyboard
    pGraphics->AttachControl(new IVKeyboardControl(keyboardBounds, 19, 79), kCtrlTagKeyboard);
    // Mod wheels
    pGraphics->AttachControl(new IWheelControl(wheelsBounds.FracRectHorizontal(0.5)), kCtrlTagBender);
    pGraphics->AttachControl(new IWheelControl(wheelsBounds.FracRectHorizontal(0.5, true), IMidiMsg::EControlChangeMsg::kModWheel));
    //    pGraphics->AttachControl(new IVMultiSliderControl<4>(b.GetGridCell(0, 2, 2).GetPadded(-30), "", TABLITSA_STYLE, kParamAttack, EDirection::Vertical, 0.f, 1.f));

#endif
    const ISVG tableImg = pGraphics->LoadSVG(PERIODIC_TABLE_FN);
    const ISVG routingImg = pGraphics->LoadSVG(ROUTING_FN);
    const IBitmap logoImg = pGraphics->LoadBitmap(LOGO_FN);

    // Periodic table
    pGraphics->AttachControl(new PeriodicTable(globalControls.GetReducedFromBottom(globalControls.H() * 0.5f).GetFromLeft(720.f).GetPadded(-30.f).GetVShifted(20.f), tableImg, { kNoParameter }),
      kCtrlTagPeriodicTable);
    int default_elements[]{ 1, 2 };
    dynamic_cast<PeriodicTable*>(pGraphics->GetControlWithTag(kCtrlTagPeriodicTable))->LoadValues(default_elements);
    pGraphics->AttachControl(new IBitmapControl(globalControls.GetReducedFromTop(90.f).GetReducedFromLeft(110.f).GetFromTLHC(280.f, 80.), logoImg));

    // Preset controls
    const IRECT presetPanel = globalControls.GetFromTop(50.f).GetFromLeft(720.f).GetHPadded(-20.f);
    IActionFunction savePresetFunction = IActionFunction([this, pGraphics](IControl* pControl) {
      bool saveBackup = false;
      SendArbitraryMsgFromUI(kMsgSavePreset, kNoTag, sizeof(bool), reinterpret_cast<const void*>(&saveBackup));
      dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagDefaultPresetList))->SetCustomStr(mPresetName);
      });
    IActionFunction loadPresetFunction = IActionFunction([this](IControl* pControl) {
      bool loadBackup = false;
      SendArbitraryMsgFromUI(kMsgLoadPreset, kNoTag, sizeof(bool), reinterpret_cast<const void*>(&loadBackup));
      });
    IActionFunction resetToDefaultFunction = IActionFunction([this, pGraphics](IControl* pControl) {
      IByteChunk currentState;
      SerializeState(currentState);
      SavePreset(currentState, "", true);
      LoadDefaultState();
      strcpy(mPresetName, "Default");
      dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagDefaultPresetList))->SetCustomStr(mPresetName);
      ResetFilterControls(pGraphics, 3);
      });
    IActionFunction loadBackupFunction = IActionFunction([this](IControl* pControl) {
      IByteChunk presetData = LoadPreset("", true);
      if (presetData.Size() == 0)
        return;
      UnserializeState(presetData, 0);
      UpdateUIControls();
      });
    const IVStyle ioButtonStyle = pushButtonStyle.WithColor(EVColor::kFG, COLOR_WHITE.WithOpacity(0.2));
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(0, 1, 2, 8).GetPadded(-2.f).GetHShifted(25.f), SplashClickActionFunc, "Save", ioButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(savePresetFunction);
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(0, 2, 2, 8).GetPadded(-2.f).GetHShifted(25.f), SplashClickActionFunc, "Load", ioButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(loadPresetFunction);
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(1, 1, 2, 8).GetPadded(-2.f).GetHShifted(25.f), SplashClickActionFunc, "Reset", ioButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(resetToDefaultFunction);
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(1, 2, 2, 8).GetPadded(-2.f).GetHShifted(25.f), SplashClickActionFunc, "Backup", ioButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(loadBackupFunction);
    // Built-in file selection
    pGraphics->AttachControl(new ITextControl(presetPanel.SubRectVertical(2, 0).GetFromLeft(110.f).GetPadded(-2.f), "Default Presets:", TABLITSA_TEXT), kNoTag, "Presets");
    pGraphics->AttachControl(new DropdownListControl(presetPanel.SubRectVertical(2, 1).GetFromLeft(110.f).GetPadded(-2.f), PRESET_NAMES), kCtrlTagDefaultPresetList, "Presets")->SetActionFunction(
      [this, pGraphics](IControl* pControl) {
        mPresetID = dynamic_cast<DropdownListControl*>(pControl)->GetCurrentIndex();
        dynamic_cast<DropdownListControl*>(pControl)->SetCustomStr("");
        // Load Resource
        HMODULE hInstance = static_cast<HMODULE>(pGraphics->GetWinModuleHandle());
#ifdef VST3_API
        std::wstring presetName = PRESET_ID_LIST_W[mPresetID];
        HRSRC hRes{ FindResource(hInstance, presetName.c_str(), (LPWSTR)RT_RCDATA) };
        HGLOBAL pData{ LoadResource(hInstance, hRes) };
#else
        std::string presetName = PRESET_ID_LIST[mPresetID];
        HRSRC hRes{ FindResource(hInstance, presetName.c_str(), RT_RCDATA) };
        HGLOBAL pData{ LoadResource(hInstance, hRes) };
#endif
        memcpy(mPresetName, PRESET_ID_LIST[mPresetID], PRESET_NAME_CHAR_LENGTH);
        if (!pData)
          return;
        char* resLock{ static_cast<char*>(LockResource(pData)) };
        const int byteLen = SizeofResource(hInstance, hRes);
        // Set State
        IByteChunk preset;
        preset.PutBytes(pData, byteLen);
        UnserializeState(preset, 0);
        UpdateUIControls();
        FreeResource(pData);
      }
    );
    dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagDefaultPresetList))->SetCurrentIdx(std::find(PRESET_ID_LIST, PRESET_ID_LIST + N_PRESETS, "Default") - PRESET_ID_LIST, false);
#if !LOAD_WITH_DIALOG
    const IRECT presetMenuBounds = presetPanel.GetGridCell(1, 1, 2, 4).GetVShifted(presetPanel.GetGridCell(1, 1, 2, 5).H()).GetPadded(0.f, 0.f, 100.f, 100.f);
    IPopupMenuControl* presetSelectionMenu = new IPopupMenuControl(kNoParameter, TABLITSA_TEXT);
    pGraphics->AttachControl(new PresetSelector(presetPanel.GetGridCell(1, 1, 2, 4).GetHPadded(-10.f), presetSelectionMenu, {"Preset 1", "Preset 2"}));
#endif
    // Wavetable controls
    const IRECT osc1Panel = globalControls.GetReducedFromLeft(720.f).GetReducedFromBottom(globalControls.H() / 2.f).GetGridCell(0, 0, 2, 3).GetPadded(-10.f).GetHPadded(-12.f);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc1Panel.GetGridCell(0, 0, 2, 2).GetCentredInside(dialSize), kParamWavetable1Pitch, "Pitch", modKnobStyle, true), kNoTag, "Osc1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc1Panel.GetGridCell(0, 1, 2, 2).GetCentredInside(dialSize), kParamWavetable1Amp, "Amp"), kNoTag, "Osc1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc1Panel.GetGridCell(1, 0, 2, 3).GetCentredInside(dialSize), kParamWavetable1Pos, "Valency"), kNoTag, "Osc1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc1Panel.GetGridCell(1, 1, 2, 3).GetCentredInside(dialSize), kParamWavetable1Bend, "Isotope"), kNoTag, "Osc1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc1Panel.GetGridCell(1, 2, 2, 3).GetCentredInside(dialSize), kParamWavetable1Formant, "Excitation"), kNoTag, "Osc1");
    pGraphics->AttachControl(new TablitsaVGroupControl("Wavetable 1", "Osc1", -5.f, 15.f, -5.f, 2.f));

    const IRECT osc2Panel = globalControls.GetReducedFromLeft(720.f).GetReducedFromBottom(globalControls.H() / 2.f).GetGridCell(1, 0, 2, 3).GetPadded(-10.f).GetHPadded(-12.f);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc2Panel.GetGridCell(0, 0, 2, 2).GetCentredInside(dialSize), kParamWavetable2Pitch, "Pitch", modKnobStyle, true), kNoTag, "Osc2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc2Panel.GetGridCell(0, 1, 2, 2).GetCentredInside(dialSize), kParamWavetable2Amp, "Amp"), kNoTag, "Osc2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc2Panel.GetGridCell(1, 0, 2, 3).GetCentredInside(dialSize), kParamWavetable2Pos, "Valency"), kNoTag, "Osc2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc2Panel.GetGridCell(1, 1, 2, 3).GetCentredInside(dialSize), kParamWavetable2Bend, "Isotope"), kNoTag, "Osc2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc2Panel.GetGridCell(1, 2, 2, 3).GetCentredInside(dialSize), kParamWavetable2Formant, "Excitation"), kNoTag, "Osc2");
    pGraphics->AttachControl(new TablitsaVGroupControl("Wavetable 2", "Osc2", -5.f, 15.f, -5.f, 2.f));

    /* Master Controls */
    const IRECT masterPanel = globalControls.GetReducedFromLeft(720.f + osc1Panel.W() + 30.f).GetReducedFromBottom(globalControls.H() / 2.f).GetPadded(-10.f);
    // Gain, Pan, and Volume Meter
    pGraphics->AttachControl(new TablitsaIVKnobControl(masterPanel.GetGridCell(0, 0, 4, 6).GetVShifted(-15.f).GetCentredInside(dialSize), kParamGain, "Gain"), kNoTag, "Master");
    // pGraphics->AttachControl(new IBKnobControl(masterPanel.GetGridCell(0, 0, 4, 6).GetCentredInside(dialSize), dialImg, kParamGain), kNoTag, "Master");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(masterPanel.GetGridCell(0, 1, 4, 6).GetVShifted(-15.f).GetCentredInside(dialSize), kParamPan, "Pan"), kNoTag, "Master");
    // Static modulators
    const IRECT staticModPanel = masterPanel.GetGridCell(1, 0, 2, 3).GetHAlignedTo(masterPanel.GetGridCell(0, 0, 1, 6), EAlign::Near).GetReducedFromRight(25.f).GetPadded(-12.f);
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetGridCell(0, 0, 3, 1), kNoParameter, "Velocity", TABLITSA_STYLE, false, EDirection::Horizontal), kCtrlTagVelDepth, "StaticMods");
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetGridCell(1, 0, 3, 1), kNoParameter, "Keytrack", TABLITSA_STYLE, false, EDirection::Horizontal), kCtrlTagKTkDepth, "StaticMods");
    pGraphics->AttachControl(new ModSliderControl(staticModPanel.GetGridCell(2, 0, 3, 1), kNoParameter, "Trg. Rnd.", TABLITSA_STYLE, false, EDirection::Horizontal), kCtrlTagRndDepth, "StaticMods");
    pGraphics->AttachControl(new TablitsaVGroupControl("Midi Mods", "StaticMods", 5.f, 20.f, 5.f, 15.f));
    // Polyphony and Portamento
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(3, 1, 12, 6).GetVShifted(-15.f).GetVPadded(-2.f), kParamPortamentoMode, "", toggleStyle, "Rate", "Time"),
      kCtrlTagGlideMode, "Master")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
        bool constantTime = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamNoteGlideTime, !constantTime);
        pGraphics->HideControl(kParamNoteGlideRate, constantTime);
        pGraphics->SetAllControlsDirty();
        });
    // Mono/Poly and Portamento Rate/Time switches
    pGraphics->AttachControl(new TablitsaIVKnobControl(masterPanel.GetGridCell(1, 0, 4, 6).GetVShifted(-15.f).GetPadded(-5.f).GetCentredInside(dialSize), kParamNoteGlideRate, "Glide Rate"), kNoTag, "Master")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVKnobControl(masterPanel.GetGridCell(1, 0, 4, 6).GetVShifted(-15.f).GetPadded(-5.f).GetCentredInside(dialSize), kParamNoteGlideTime, "Glide Time"), kNoTag, "Master")->DisablePrompt(false);
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(4, 1, 12, 6).GetVShifted(-15.f).GetVPadded(-2.f), kParamMonophonic, "", toggleStyle, "Poly", "Mono"));
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(5, 1, 12, 6).GetVShifted(-15.f).GetVPadded(-2.f), kParamLegato, "", toggleStyle, "Retrigger", "Legato"));
    // Unison and Stereo Width
    const IRECT unisonPanel = masterPanel.GetReducedFromLeft(masterPanel.W() / 3.f).GetReducedFromRight(masterPanel.W() / 6.f).GetGridCell(0, 0, 2, 1).GetVShifted(-12.f).GetHPadded(-10.f).GetHShifted(2.f);
    pGraphics->AttachControl(new TablitsaIVKnobControl(unisonPanel.GetGridCell(0, 0, 2, 3).GetHPadded(-10.f).GetHShifted(5.f), kParamUnisonVoices, "Voices", modKnobStyle), kNoTag, "Unison");
    pGraphics->AttachControl(new TablitsaIVKnobControl(unisonPanel.GetGridCell(1, 0, 2, 2).GetCentredInside(dialSize), kParamUnisonDetune, "Detune", modKnobStyle), kNoTag, "Unison");
    pGraphics->AttachControl(new TablitsaIVKnobControl(unisonPanel.FracRectVertical(0.5f, true).GetGridCell(0, 1, 1, 3).GetHPadded(-10.f).GetHShifted(-10.f), kParamUnisonChord, "Chord", modKnobStyle), kNoTag, "Unison");
    pGraphics->AttachControl(new TablitsaSliderControl(unisonPanel.FracRectVertical(0.5f, true).FracRectHorizontal(2.f / 3, true).GetGridCell(0, 1, 1, 3), kParamVibratoSpeed, "Rate"), kNoTag, "Unison");
    pGraphics->AttachControl(new TablitsaSliderControl(unisonPanel.FracRectVertical(0.5f, true).FracRectHorizontal(2.f / 3, true).GetGridCell(0, 2, 1, 3), kParamVibratoDepth, "Depth "), kNoTag, "Unison");
    pGraphics->AttachControl(new TablitsaIVKnobControl(unisonPanel.GetGridCell(1, 1, 2, 2).GetPadded(-5.f).GetCentredInside(dialSize), kParamStereoSpread, "Spread"), kNoTag, "Unison");
    pGraphics->AttachControl(new TablitsaVGroupControl("Unison", "Unison", 0.f, 13.f, 0.f, 0.f));
    pGraphics->AttachControl(new IVLEDMeterControl<2>(masterPanel.GetGridCell(0, 5, 2, 6).GetVShifted(-15.f).GetHShifted(-4.f).GetVPadded(-8.f)), kCtrlTagMeter, "Master");

    /* Effects Panel */
    const IRECT effectsPanel = IRECT(masterPanel.GetGridCell(1, 1, 2, 7).R, masterPanel.GetGridCell(1, 0, 2, 7).T, masterPanel.GetGridCell(1, 6, 2, 7).R, masterPanel.GetGridCell(1, 6, 2, 7).B).GetHPadded(-5.f).GetHShifted(2.f);
    const IRECT effectsControlsPanel = effectsPanel.GetReducedFromTop(40.f);

    IActionFunction voiceFxSlotChanged = [this, pGraphics](IControl* pControl) {
      int currentSlot = static_cast<int>(pControl->GetValue() * (TABLITSA_MAX_VOICE_EFFECTS - 1));
      char str[2];
      sprintf(str, "%d", currentSlot + 1);
      dynamic_cast<IVSlideSwitchControl*>(pControl)->SetValueStr(str);

      // Set the controls' parameters to those of the current effect slot
      constexpr int numEffectParams = kParamVoiceEffect2Param1 - kParamVoiceEffect1Param1;
      dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob1))->SetParamIdx(kParamVoiceEffect1Param1 + currentSlot * numEffectParams);
      dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob2))->SetParamIdx(kParamVoiceEffect1Param2 + currentSlot * numEffectParams);
      dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob3))->SetParamIdx(kParamVoiceEffect1Param3 + currentSlot * numEffectParams);
      dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsKnob4))->SetParamIdx(kParamVoiceEffect1Param4 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsToggle1)->SetParamIdx(kParamVoiceEffect1Param5 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsToggle2)->SetParamIdx(kParamVoiceEffect1Param6 + currentSlot * numEffectParams);

      for (int tag{ kCtrlTagVoiceEffectsKnob1 }; tag <= kCtrlTagVoiceEffectsKnob4; ++tag)
      {
        TablitsaIVModKnobControl* modKnob = dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(tag));
        if (modKnob->IsActive())
         modKnob->LoadModParams();
      }

      // Set the dropdown list's index to that of the effect which is currently active in this slot
      dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsList))->SetCurrentIdx(mVoiceEffectSlots[currentSlot], false);
      SwapVoiceEffectsUI(dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsList)), pGraphics, this, false); // Set parameters, keeping previous values
    };
    IActionFunction masterFxSlotChanged = [this, pGraphics](IControl* pControl) {
      int currentSlot = static_cast<int>(pControl->GetValue() * (TABLITSA_MAX_MASTER_EFFECTS - 1));
      char str[2];
      sprintf(str, "%d", currentSlot + 1);
      dynamic_cast<IVSlideSwitchControl*>(pControl)->SetValueStr(str);

      // Set the controls' parameters to those of the current effect slot
      constexpr int numEffectParams = 6;
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob1)->SetParamIdx(kParamMasterEffect1Param1 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob2)->SetParamIdx(kParamMasterEffect1Param2 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob3)->SetParamIdx(kParamMasterEffect1Param3 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsKnob4)->SetParamIdx(kParamMasterEffect1Param4 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsToggle1)->SetParamIdx(kParamMasterEffect1Param5 + currentSlot * numEffectParams);
      pGraphics->GetControlWithTag(kCtrlTagMasterEffectsToggle2)->SetParamIdx(kParamMasterEffect1Param6 + currentSlot * numEffectParams);

      // Set the dropdown list's index to that of the effect which is currently active in this slot
      dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagMasterEffectsList))->SetCurrentIdx(mMasterEffectSlots[currentSlot], false);
      SwapMasterEffectsUI(dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagMasterEffectsList)), pGraphics, this, false); // Set parameters, keeping previous values
    };

    /* Effects Tab Control */
    pGraphics->AttachControl(new TablitsaEffectBankControl(effectsPanel, { "Voice FX", "Master FX", "Routing" }, { "E1", "E2", "E3" }), kCtrlTagEffectBank, "Effects")->SetActionFunction(
      [this, pGraphics](IControl* pControl) {
        dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagVoiceEffectsList))->Collapse();
        dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagMasterEffectsList))->Collapse();
        mCurrentEffectsTab = dynamic_cast<TablitsaEffectBankControl*>(pControl)->GetActiveTabIdx();
      });
    
    /* Voice Effects */
    pGraphics->AttachControl(new ITextControl(effectsPanel.GetGridCell(0, 2, 1, 3).GetFromTop(60.f).GetHShifted(-55.f), "FX Slot:", TABLITSA_TEXT), kNoTag, "E1");
    pGraphics->AttachControl(new IVSlideSwitchControl(effectsPanel.GetGridCell(0, 2, 1, 3).GetFromTop(40.f).GetReducedFromTop(22.f).GetHPadded(-20.f), voiceFxSlotChanged, "", toggleStyle, true, EDirection::Horizontal, TABLITSA_MAX_VOICE_EFFECTS, 0), kCtrlTagVoiceEffectsSwitch, "E1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 0, 1, 4).GetCentredInside(dialSize), kParamVoiceEffect1Param1, "Left", modKnobStyle), kCtrlTagVoiceEffectsKnob1, "E1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 1, 1, 4).GetCentredInside(dialSize), kParamVoiceEffect1Param2, "Right", modKnobStyle), kCtrlTagVoiceEffectsKnob2, "E1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 2, 1, 4).GetCentredInside(dialSize), kParamVoiceEffect1Param3, "Feedback", modKnobStyle), kCtrlTagVoiceEffectsKnob3, "E1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 3, 1, 4).GetCentredInside(dialSize), kParamVoiceEffect1Param4, "Mix", modKnobStyle), kCtrlTagVoiceEffectsKnob4, "E1");
    pGraphics->AttachControl(new TablitsaVToggleButton(effectsControlsPanel.FracRectVertical(1.f / 3).GetGridCell(0, 1, 1, 4).GetPadded(-10.f), kParamVoiceEffect1Param5, "", toggleStyle), kCtrlTagVoiceEffectsToggle1, "E1");
    pGraphics->AttachControl(new TablitsaVToggleButton(effectsControlsPanel.FracRectVertical(1.f / 3).GetGridCell(0, 2, 1, 4).GetPadded(-10.f), kParamVoiceEffect1Param6, "", toggleStyle), kCtrlTagVoiceEffectsToggle2, "E1");
    // Selector list
    pGraphics->AttachControl(new DropdownListControl(effectsPanel.FracRectHorizontal(0.4f, false).GetFromTop(50.f).GetReducedFromTop(20.f).GetPadded(-2.f), TABLITSA_VOICE_EFFECTS_LIST), kCtrlTagVoiceEffectsList, "E1")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        SwapVoiceEffectsUI(pControl, pGraphics, this, true);
      });

    /* MASTER EFFECTS */
    pGraphics->AttachControl(new ITextControl(effectsPanel.GetGridCell(0, 2, 1, 3).GetFromTop(60.f).GetHShifted(-55.f), "FX Slot:", TABLITSA_TEXT), kNoTag, "E2");
    pGraphics->AttachControl(new IVSlideSwitchControl(effectsPanel.GetGridCell(0, 2, 1, 3).GetFromTop(40.f).GetReducedFromTop(22.f).GetHPadded(-20.f), masterFxSlotChanged, "", toggleStyle, true, EDirection::Horizontal, 3, 0), kCtrlTagMasterEffectsSwitch, "E2");
    pGraphics->AttachControl(new TablitsaIVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 0, 1, 4).GetCentredInside(dialSize), kParamMasterEffect1Param1, "Left", modKnobStyle), kCtrlTagMasterEffectsKnob1, "E2");
    pGraphics->AttachControl(new TablitsaIVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 1, 1, 4).GetCentredInside(dialSize), kParamMasterEffect1Param2, "Right", modKnobStyle), kCtrlTagMasterEffectsKnob2, "E2");
    pGraphics->AttachControl(new TablitsaIVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 2, 1, 4).GetCentredInside(dialSize), kParamMasterEffect1Param3, "Feedback", modKnobStyle), kCtrlTagMasterEffectsKnob3, "E2");
    pGraphics->AttachControl(new TablitsaIVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 3, 1, 4).GetCentredInside(dialSize), kParamMasterEffect1Param4, "Mix", modKnobStyle), kCtrlTagMasterEffectsKnob4, "E2");
    pGraphics->AttachControl(new TablitsaVToggleButton(effectsControlsPanel.FracRectVertical(1.f / 3).GetGridCell(0, 1, 1, 4).GetPadded(-10.f), kParamMasterEffect1Param5, "", toggleStyle), kCtrlTagMasterEffectsToggle1, "E2");
    pGraphics->AttachControl(new TablitsaVToggleButton(effectsControlsPanel.FracRectVertical(1.f / 3).GetGridCell(0, 2, 1, 4).GetPadded(-10.f), kNoParameter, "", toggleStyle), kCtrlTagMasterEffectsToggle2, "E2");
    // Selector List
    pGraphics->AttachControl(new DropdownListControl(effectsPanel.FracRectHorizontal(0.4f, false).GetFromTop(50.f).GetReducedFromTop(20.f).GetPadded(-2.f), TABLITSA_MASTER_EFFECTS_LIST), kCtrlTagMasterEffectsList, "E2")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        SwapMasterEffectsUI(pControl, pGraphics, this, true);
      });

    /* Routing */
    const IRECT routingPanel = effectsControlsPanel.GetPadded(-5.f, 10.f, -5.f, 0.f);
    pGraphics->AttachControl(new ISVGControl(routingPanel, routingImg), kNoTag, "E3");
    pGraphics->AttachControl(new TablitsaSliderControl(routingPanel.GetGridCell(0, 1, 5, 4).GetHPadded(-6.f).GetVShifted(-10.f), kParamOsc1FilterBypass, "Osc1 Filter BP", TABLITSA_STYLE.WithShowValue(false).WithShowLabel(false), false, EDirection::Horizontal), kNoTag, "E3");
    pGraphics->AttachControl(new TablitsaSliderControl(routingPanel.GetGridCell(1, 1, 5, 4).GetHPadded(-6.f).GetVShifted(-5.f), kParamOsc2FilterBypass, "Osc2 Filter BP", TABLITSA_STYLE.WithShowValue(false).WithShowLabel(false), false, EDirection::Horizontal), kNoTag, "E3");
    pGraphics->AttachControl(new TablitsaSliderControl(routingPanel.GetGridCell(3, 1, 5, 4).GetHShifted(-10.f), kParamOsc2EffectBypass, "Osc2 Effect BP", TABLITSA_STYLE.WithShowValue(false).WithShowLabel(false), false, EDirection::Horizontal), kNoTag, "E3");
    pGraphics->AttachControl(new TablitsaSliderControl(routingPanel.GetGridCell(4, 1, 5, 4).GetHShifted(-10.f).GetVShifted(2.f), kParamOsc1EffectBypass, "Osc1 Effect BP", TABLITSA_STYLE.WithShowValue(false).WithShowLabel(false), false, EDirection::Horizontal), kNoTag, "E3");

    /* LOWER-HALF CONTROLS */
    constexpr float modulePad[4]{-10.f, -5.f, -10.f, -15.f};

    /* Waveform Modulation Panel */
    const IRECT lowerPanel = globalControls.GetReducedFromTop(globalControls.H() * 0.5f).GetPadded(-10.f).GetReducedFromRight(250.f);
    const IRECT oscModPanel = lowerPanel.GetGridCell(0, 0, 3, 3);
    pGraphics->AttachControl(new ICaptionControl(oscModPanel.GetVPadded(-10.f).GetGridCell(0, 0, 3, 3).GetVShifted(10.f), kParamOscModulator, dropdownText), kNoTag, "OscMod")->SetActionFunction(
      [pGraphics](IControl* pControl) {
        bool phaseMod = pControl->GetValue() < 0.5;
        // Update parameter index
        pGraphics->GetControlWithTag(kCtrlTagOscModFreq)->SetParamIdx(phaseMod ? kParamPhaseModFreq : kParamRingModFreq);
        pGraphics->GetControlWithTag(kCtrlTagOscModAmt)->SetParamIdx(phaseMod ? kParamPhaseModAmount : kParamRingModAmount);
        // Update modulation index
        dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagOscModFreq))->GetModulationFrom(phaseMod ? kParamPhaseModFreq : kParamRingModFreq);
        dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagOscModAmt))->GetModulationFrom(phaseMod ? kParamPhaseModAmount : kParamRingModAmount);
        IControl* pOsc1Switch = pGraphics->GetControlWithTag(kCtrlTagOsc1ModSwitch);
        pOsc1Switch->SetParamIdx(phaseMod ? kParamOsc1PM : kParamOsc1RM);
        IControl* pOsc2Switch = pGraphics->GetControlWithTag(kCtrlTagOsc2ModSwitch);
        pOsc2Switch->SetParamIdx(phaseMod ? kParamOsc2PM : kParamOsc2RM);
        // Update controls
        pControl->GetDelegate()->SendCurrentParamValuesFromDelegate();
        pOsc2Switch->SetDirty(false);
        pOsc1Switch->SetDirty(false);
      });
    pGraphics->AttachControl(new IVToggleControl(oscModPanel.GetVPadded(-10.f).GetVPadded(-10.f).GetGridCell(1, 0, 3, 3).GetMidHPadded(35.f), kParamOsc1PM, "Osc1", toggleStyle, "Osc 1 OFF", "Osc 1 ON"), kCtrlTagOsc1ModSwitch, "OscMod");
    pGraphics->AttachControl(new IVToggleControl(oscModPanel.GetVPadded(-10.f).GetVPadded(-10.f).GetGridCell(2, 0, 3, 3).GetMidHPadded(35.f), kParamOsc2PM, "Osc1", toggleStyle, "Osc 2 OFF", "Osc 2 ON"), kCtrlTagOsc2ModSwitch, "OscMod");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(oscModPanel.GetVPadded(-10.f).GetGridCell(0, 1, 1, 3).GetCentredInside(dialSize), kParamPhaseModFreq, "Pitch"), kCtrlTagOscModFreq, "OscMod");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(oscModPanel.GetGridCell(0, 2, 1, 3).GetCentredInside(dialSize), kParamPhaseModAmount, "Amount"), kCtrlTagOscModAmt, "OscMod");
    pGraphics->AttachControl(new TablitsaVGroupControl(oscModPanel.GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "Oscillator Modulation"));
    
    /* Filter Panel */
    /* Filter 1 */
    const IRECT filter1Panel = lowerPanel.GetGridCell(1, 0, 3, 3).GetReducedFromRight(10.f);
    pGraphics->AttachControl(new ICaptionControl(filter1Panel.GetGridCell(0, 0, 3, 4).GetVShifted(20.f), kParamFilter1Type, dropdownText), kCtrlTagFilter1Type, "Flt1")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        IControl* filterModeControl{ pGraphics->GetControlWithTag(kCtrlTagFilter1Mode) }; // Control containing list of filter modes for a given type
        int filterModeParam{ kParamFilter1ModeVSF + std::max(static_cast<int>(pControl->GetValue() * (kNumFilters - 1.) - 1), 0) };
        // Get range of possible modes
        filterModeControl->SetParamIdx(filterModeParam);
        filterModeControl->SetValue(GetParam(filterModeParam)->GetNormalized());
        // Update Controls
        if (pControl->GetValue() == kNoFilter)
          ResetFilterControls(pGraphics, 1);
        else if (pControl->GetValue() * ((int)kNumFilters - 1) == kComb)
        {
          filterModeControl->SetDisabled(true);
          // Swap Cutoff/Res/Drive for FF/FB/Delay
          pGraphics->HideControl(kParamFilter1Cutoff, true);
          pGraphics->HideControl(kParamFilter1FF, false);
          pGraphics->HideControl(kParamFilter1Resonance, true);
          pGraphics->HideControl(kParamFilter1FB, false);
          // Drive/delay switch
          pGraphics->HideControl(kParamFilter1Drive, true);
          pGraphics->HideControl(kParamFilter1Delay, false);
        }
        else
        {
          filterModeControl->SetDisabled(false);
          pGraphics->HideControl(kParamFilter1Cutoff, false);
          pGraphics->HideControl(kParamFilter1FF, true);
          pGraphics->HideControl(kParamFilter1Resonance, false);
          pGraphics->HideControl(kParamFilter1FB, true);
          // Drive/delay switch
          pGraphics->HideControl(kParamFilter1Drive, false);
          pGraphics->HideControl(kParamFilter1Delay, true);
        }
        filterModeControl->SetDirty(true);
      }
    );
    pGraphics->AttachControl(new ICaptionControl(filter1Panel.GetGridCell(1, 0, 3, 4).GetFromTop(20.f).GetVShifted(20.f), kParamFilter1ModeVSF, dropdownText), kCtrlTagFilter1Mode, "Flt1")->
      SetDisabled(static_cast<int>(pGraphics->GetControlWithTag(kCtrlTagFilter1Type)->GetValue()) == kNoFilter);
    // Cutoff/Resonance/Drive
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 1, 1, 5).GetCentredInside(dialSize), kParamFilter1Cutoff, "Cutoff"), kCtrlTagFilter1Cutoff, "Flt1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 2, 1, 5).GetCentredInside(dialSize), kParamFilter1Resonance, "Res."), kCtrlTagFilter1Resonance, "Flt1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 3, 1, 5).GetCentredInside(dialSize), kParamFilter1Drive, "Drive"), kCtrlTagFilter1Drive, "Flt1");
    pGraphics->AttachControl(new TablitsaSliderControl(filter1Panel.GetGridCell(0, 8, 1, 10).GetVPadded(-20.f), kParamFilter1Osc1Send, "O1", modKnobStyle.WithShowValue(false)), kCtrlTagFilter1Osc1, "Flt1");
    pGraphics->AttachControl(new TablitsaSliderControl(filter1Panel.GetGridCell(0, 9, 1, 10).GetVPadded(-20.f), kParamFilter1Osc2Send, "O2", modKnobStyle.WithShowValue(false)), kCtrlTagFilter1Osc2, "Flt1");
    // FF/FB/Delay Time (comb filter)
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 1, 1, 5).GetCentredInside(dialSize), kParamFilter1FF, "F/F"), kCtrlTagFilter1FF, "Flt1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 2, 1, 5).GetCentredInside(dialSize), kParamFilter1FB, "F/B"), kCtrlTagFilter1FB, "Flt1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 3, 1, 5).GetCentredInside(dialSize), kParamFilter1Delay, "Delay"), kCtrlTagFilter1Delay, "Flt1")->Hide(true);

    pGraphics->AttachControl(new TablitsaVGroupControl(filter1Panel.GetPadded(modulePad[0], modulePad[1], modulePad[2] + 10.f, modulePad[3]), "Filter 1"));

    // Set comb filter controls' modulation parameters to the normal cutoff/res/drive controls' ones
    dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1FF))->GetModulationFrom(kParamFilter1Cutoff);
    dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1FB))->GetModulationFrom(kParamFilter1Resonance);
    dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1Delay))->GetModulationFrom(kParamFilter1Drive);

    /* Filter 2 */
    const IRECT filter2Panel = lowerPanel.GetGridCell(2, 0, 3, 3).GetReducedFromRight(10.f);
    pGraphics->AttachControl(new ICaptionControl(filter2Panel.GetGridCell(0, 0, 3, 4).GetVShifted(20.f), kParamFilter2Type, dropdownText), kCtrlTagFilter2Type, "Flt2")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        IControl* filterModeControl{ pGraphics->GetControlWithTag(kCtrlTagFilter2Mode) };
        int filterModeParam{ kParamFilter2ModeVSF + std::max(static_cast<int>(pControl->GetValue() * (kNumFilters - 1.) - 1), 0) };
        filterModeControl->SetParamIdx(filterModeParam);
        filterModeControl->SetValue(GetParam(filterModeParam)->GetNormalized());
        if (pControl->GetValue() == kNoFilter)
          ResetFilterControls(pGraphics, 2);
        else if (pControl->GetValue() * ((int)kNumFilters - 1) == kComb)
        {
          filterModeControl->SetDisabled(true);
          // Swap Cutoff/Res/Drive for FF/FB/Delay
          pGraphics->HideControl(kParamFilter2Cutoff, true);
          pGraphics->HideControl(kParamFilter2FF, false);
          pGraphics->HideControl(kParamFilter2Resonance, true);
          pGraphics->HideControl(kParamFilter2FB, false);
          pGraphics->HideControl(kParamFilter2Drive, true);
          pGraphics->HideControl(kParamFilter2Delay, false);
        }
        else
        {
          filterModeControl->SetDisabled(false);
          pGraphics->HideControl(kParamFilter2Cutoff, false);
          pGraphics->HideControl(kParamFilter2FF, true);
          pGraphics->HideControl(kParamFilter2Resonance, false);
          pGraphics->HideControl(kParamFilter2FB, true);
          pGraphics->HideControl(kParamFilter2Drive, false);
          pGraphics->HideControl(kParamFilter2Delay, true);
        }
        filterModeControl->SetDirty(true);
      }
    );
    pGraphics->AttachControl(new ICaptionControl(filter2Panel.GetGridCell(1, 0, 3, 4).GetFromTop(20.f).GetVShifted(20.f), kParamFilter2ModeVSF, dropdownText), kCtrlTagFilter2Mode, "Flt2")->
      SetDisabled(static_cast<int>(pGraphics->GetControlWithTag(kCtrlTagFilter2Type)->GetValue()) == kNoFilter);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 1, 1, 5).GetCentredInside(dialSize), kParamFilter2Cutoff, "Cutoff"), kNoTag, "Flt2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 2, 1, 5).GetCentredInside(dialSize), kParamFilter2Resonance, "Res."), kNoTag, "Flt2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 3, 1, 5).GetCentredInside(dialSize), kParamFilter2Drive, "Drive"), kNoTag, "Flt2");
    pGraphics->AttachControl(new TablitsaSliderControl(filter2Panel.GetGridCell(0, 8, 1, 10).GetVPadded(-20.f), kParamFilter2Osc1Send, "O1", modKnobStyle.WithShowValue(false)), kCtrlTagFilter2Osc1, "Flt2");
    pGraphics->AttachControl(new TablitsaSliderControl(filter2Panel.GetGridCell(0, 9, 1, 10).GetVPadded(-20.f), kParamFilter2Osc2Send, "O2", modKnobStyle.WithShowValue(false)), kCtrlTagFilter2Osc2, "Flt2");
    // FF/FB/Delay Time (comb filter)
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 1, 1, 5).GetCentredInside(dialSize), kParamFilter2FF, "F/F"), kCtrlTagFilter2FF, "Flt1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 2, 1, 5).GetCentredInside(dialSize), kParamFilter2FB, "F/B"), kCtrlTagFilter2FB, "Flt1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 3, 1, 5).GetCentredInside(dialSize), kParamFilter2Delay, "Delay"), kCtrlTagFilter2Delay, "Flt1")->Hide(true);

    pGraphics->AttachControl(new TablitsaVGroupControl(filter2Panel.GetPadded(modulePad[0], modulePad[1], modulePad[2] + 10.f, modulePad[3]), "Filter 2"));

    // Set comb filter controls' modulation parameters to the normal cutoff/res/drive controls' ones
    dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter2FF))->GetModulationFrom(kParamFilter2Cutoff);
    dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter2FB))->GetModulationFrom(kParamFilter2Resonance);
    dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter2Delay))->GetModulationFrom(kParamFilter2Drive);

#pragma mark Modulation globalControls -

    // Envelope Sliders
    const IRECT envPanel = lowerPanel.GetGridCell(0, 1, 1, 3);
    constexpr float envDialShift{ 5.f };
    // Envelope 1
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 0, 3, 6).GetCentredInside(dialSize), kParamEnv1Attack, "Att", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 1, 3, 6).GetCentredInside(dialSize), kParamEnv1Decay, "Dec", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 1, 3, 6).GetCentredInside(dialSize), kParamEnv1DecayCurve, "Exp Dec", modKnobStyle), kNoTag, "Env1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 2, 3, 6).GetCentredInside(dialSize), kParamEnv1Sustain, kParamEnv1SustainEnv1, "Sus", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 3, 3, 6).GetCentredInside(dialSize), kParamEnv1Release, "Rel", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 3, 3, 6).GetCentredInside(dialSize), kParamEnv1ReleaseCurve, "Exp Rel", modKnobStyle), kNoTag, "Env1")->Hide(true);
    pGraphics->AttachControl(new TablitsaSliderControl(envPanel.GetGridCell(0, 4, 3, 6).GetVPadded(-15.f), kParamEnv1Velocity, "Velocity", TABLITSA_STYLE, false, EDirection::Vertical, 4.f, 8.f, 2.f, true), kNoTag, "Env1");
    pGraphics->AttachControl(new ModSliderControl(envPanel.GetGridCell(0, 5, 3, 6).GetVPadded(-15.f), kNoParameter, "Depth"), kCtrlTagEnv1Depth, "Env1");
    pGraphics->AttachControl(new TablitsaVGroupControl(envPanel.GetGridCell(0, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "Env 1"));
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetHShifted(envDialShift).GetGridCell(0, 1, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamEnv1Decay, editCurve);
        pGraphics->HideControl(kParamEnv1DecayCurve, !editCurve);
        });
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetGridCell(0, 3, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamEnv1Release, editCurve);
        pGraphics->HideControl(kParamEnv1ReleaseCurve, !editCurve);
        });
    // Envelope 2
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 0, 3, 6).GetCentredInside(dialSize), kParamEnv2Attack, "Att", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 1, 3, 6).GetCentredInside(dialSize), kParamEnv2Decay, "Dec", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 1, 3, 6).GetCentredInside(dialSize), kParamEnv2DecayCurve, "Exp Dec", modKnobStyle), kNoTag, "Env2")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 2, 3, 6).GetCentredInside(dialSize), kParamEnv2Sustain, kParamEnv2SustainEnv1, "Sus", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 3, 3, 6).GetCentredInside(dialSize), kParamEnv2Release, "Rel", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 3, 3, 6).GetCentredInside(dialSize), kParamEnv2ReleaseCurve, "Exp Rel", modKnobStyle), kNoTag, "Env2")->Hide(true);
    pGraphics->AttachControl(new TablitsaSliderControl(envPanel.GetGridCell(1, 4, 3, 6).GetVPadded(-15.f), kParamEnv2Velocity, "Velocity", TABLITSA_STYLE, false, EDirection::Vertical, 4.f, 8.f, 2.f, true), kNoTag, "Env2");
    pGraphics->AttachControl(new ModSliderControl(envPanel.GetGridCell(1, 5, 3, 6).GetVPadded(-15.f), kNoParameter, "Depth"), kCtrlTagEnv2Depth, "Env2");
    pGraphics->AttachControl(new TablitsaVGroupControl(envPanel.GetGridCell(1, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "Env 2"));
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetHShifted(envDialShift).GetGridCell(1, 1, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamEnv2Decay, editCurve);
        pGraphics->HideControl(kParamEnv2DecayCurve, !editCurve);
        });
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetGridCell(1, 3, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamEnv2Release, editCurve);
        pGraphics->HideControl(kParamEnv2ReleaseCurve, !editCurve);
        });
    // Amplitude Envelope
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 0, 3, 6).GetCentredInside(dialSize), kParamAmpEnvAttack, "Att", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 1, 3, 6).GetCentredInside(dialSize), kParamAmpEnvDecay, "Dec", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 1, 3, 6).GetCentredInside(dialSize), kParamAmpEnvDecayCurve, "Dec", modKnobStyle), kNoTag, "AmpEnv")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 2, 3, 6).GetCentredInside(dialSize), kParamAmpEnvSustain, kParamAmpEnvSustainEnv1, "Sus", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 3, 3, 6).GetCentredInside(dialSize), kParamAmpEnvRelease, "Rel", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new TablitsaIVKnobControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 3, 3, 6).GetCentredInside(dialSize), kParamAmpEnvReleaseCurve, "Exp Rel", modKnobStyle), kNoTag, "AmpEnv")->Hide(true);
    pGraphics->AttachControl(new TablitsaSliderControl(envPanel.GetGridCell(2, 4, 3, 6).GetVPadded(-15.f), kParamAmpEnvVelocity, "Velocity", TABLITSA_STYLE, false, EDirection::Vertical, 4.f, 8.f, 2.f, true), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new ModSliderControl(envPanel.GetGridCell(2, 5, 3, 6).GetVPadded(-15.f), kNoParameter, "Depth"), kCtrlTagAmpEnvDepth, "AmpEnv");
    pGraphics->AttachControl(new TablitsaVGroupControl(envPanel.GetGridCell(2, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "Amp Env"));
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetHShifted(envDialShift).GetGridCell(2, 1, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamAmpEnvDecay, editCurve);
        pGraphics->HideControl(kParamAmpEnvDecayCurve, !editCurve);
        });
    pGraphics->AttachControl(new IVToggleControl(envPanel.GetGridCell(2, 3, 3, 6).GetVShifted(10.f).GetFromBottom(21.f), kNoParameter, "", toggleStyle, "Time", "Curve"),
      kNoTag)->SetActionFunction([pGraphics](IControl* pControl) {
        bool editCurve = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamAmpEnvRelease, editCurve);
        pGraphics->HideControl(kParamAmpEnvReleaseCurve, !editCurve);
        });

    /* LFO Panels */
    const IRECT lfoPanel = lowerPanel.GetGridCell(0, 2, 1, 3);
    // LFO1
    pGraphics->AttachControl(new TablitsaIVModKnobControl(lfoPanel.GetGridCell(0, 0, 3, 5).GetCentredInside(dialSize), kParamLFO1RateHz, "Rate"), kNoTag, "LFO1");
    pGraphics->AttachControl(new TablitsaIVKnobControl(lfoPanel.GetGridCell(0, 0, 3, 5).GetCentredInside(dialSize), kParamLFO1RateTempo, "Rate"), kNoTag, "LFO1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(lfoPanel.GetGridCell(0, 2, 3, 5).GetCentredInside(dialSize), kParamLFO1Amp, "Amp"), kNoTag, "LFO1");
    pGraphics->AttachControl(new ICaptionControl(lfoPanel.GetGridCell(0, 1, 6, 5).GetVPadded(-20.f), kParamLFO1Shape, dropdownText), kNoTag, "LFO1")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        SendArbitraryMsgFromUI(kMsgUpdateLFO1Plot);
      }
    );
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(1, 1, 9, 5).GetVPadded(-6.f), kParamLFO1RateMode, "", toggleStyle, "Sync Off", "Sync On"),
      kCtrlTagLFO1RateMode, "LFO1")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
      bool sync = pControl->GetValue() > 0.5;
      pGraphics->HideControl(kParamLFO1RateHz, sync);
      pGraphics->HideControl(kParamLFO1RateTempo, !sync);
      //pGraphics->SetAllControlsDirty();
    });
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(2, 1, 9, 5).GetVPadded(-6.f).GetVShifted(-15.f), kParamLFO1Restart, "Restart", toggleStyle, "Free Run", "Retrigger"), kNoTag, "LFO1");
    pGraphics->AttachControl(new ModSliderControl(lfoPanel.GetGridCell(0, 4, 3, 5).GetVPadded(-18.f), kNoParameter, "Depth"), kCtrlTagLFO1Depth, "LFO1");
    pGraphics->AttachControl(new ModPlotControl(lfoPanel.GetGridCell(0, 3, 3, 5).GetVPadded(-20.f), kParamLFO1Phase, nullptr, 1024, 25), kCtrlTagLFO1Plot, "LFO1");
    pGraphics->AttachControl(new TablitsaVGroupControl(lfoPanel.GetGridCell(0, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "LFO 1"));

    // LFO 2
    pGraphics->AttachControl(new TablitsaIVModKnobControl(lfoPanel.GetGridCell(1, 0, 3, 5).GetCentredInside(dialSize), kParamLFO2RateHz, "Rate"), kNoTag, "LFO2");
    pGraphics->AttachControl(new TablitsaIVKnobControl(lfoPanel.GetGridCell(1, 0, 3, 5).GetCentredInside(dialSize), kParamLFO2RateTempo, "Rate"), kNoTag, "LFO2")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(lfoPanel.GetGridCell(1, 2, 3, 5).GetCentredInside(dialSize), kParamLFO2Amp, "Amp"), kNoTag, "LFO2");
    pGraphics->AttachControl(new ICaptionControl(lfoPanel.GetGridCell(2, 1, 6, 5).GetVPadded(-20.f), kParamLFO2Shape, dropdownText), kNoTag, "LFO2")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        SendArbitraryMsgFromUI(kMsgUpdateLFO2Plot);
      }
    );
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(4, 1, 9, 5).GetVPadded(-6.f), kParamLFO2RateMode, "", toggleStyle, "Sync Off", "Sync On"),
      kCtrlTagLFO2RateMode, "LFO2")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
        bool sync = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamLFO2RateHz, sync);
        pGraphics->HideControl(kParamLFO2RateTempo, !sync);
        //pGraphics->SetAllControlsDirty();
        });
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(5, 1, 9, 5).GetVPadded(-6.f).GetVShifted(-15.f), kParamLFO2Restart, "Restart", toggleStyle, "Free Run", "Retrigger"), kNoTag, "LFO2");
    pGraphics->AttachControl(new ModSliderControl(lfoPanel.GetGridCell(1, 4, 3, 5).GetVPadded(-18.f), kNoParameter, "Depth"), kCtrlTagLFO2Depth, "LFO2");
    pGraphics->AttachControl(new ModPlotControl(lfoPanel.GetGridCell(1, 3, 3, 5).GetVPadded(-20.f), nullptr, 1024, 25), kCtrlTagLFO2Plot, "LFO1");
    pGraphics->AttachControl(new TablitsaVGroupControl(lfoPanel.GetGridCell(1, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "LFO 2"));

    /* Sequencer */
    const IRECT sequencerDisplayPanel = lowerPanel.GetGridCell(2, 2, 3, 3).GetFromLeft(PLUG_WIDTH - lfoPanel.L).GetVShifted(9.f).GetPadded(10.f, 12.f, 0.f, 0.f);
    const IRECT sequencerParamsPanel = IRECT(lfoPanel.R, lfoPanel.T, globalControls.R, lowerPanel.B - sequencerDisplayPanel.H()).GetPadded(-5.f, -5.f, -5.f, -10.f);

    // Sequencer Action Functions
    auto stepRandomizeAction = [pGraphics](IControl* pCaller) {
      SplashClickActionFunc(pCaller);
      pGraphics->GetDelegate()->SendArbitraryMsgFromUI(kMsgRandomizeSequencer);
    };

    pGraphics->AttachControl(new TablitsaIVKnobControl(sequencerParamsPanel.GetGridCell(0, 0, 2, 3).GetCentredInside(dialSize), kParamSequencerSteps, "Steps"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(sequencerParamsPanel.GetGridCell(0, 1, 2, 3).GetCentredInside(dialSize), kParamSequencerRateHz, "Rate"), kNoTag, "SequencerParams")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVKnobControl(sequencerParamsPanel.GetGridCell(0, 1, 2, 3).GetCentredInside(dialSize), kParamSequencerRateTempo, "Rate"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(sequencerParamsPanel.GetGridCell(0, 2, 2, 3).GetCentredInside(dialSize), kParamSequencerAmp, "Amp"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new TablitsaIVKnobControl(sequencerParamsPanel.GetGridCell(1, 0, 2, 3).GetCentredInside(dialSize), kParamSequencerGlide, "Glide"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(4, 2, 7, 3).GetPadded(-15.f, -2.f, -15.f, -2.f), kParamSequencerRateMode, "", toggleStyle, "Sync Off", "Sync On"),
      kCtrlTagSequencerRateMode, "Sequencer")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
        bool sync = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamSequencerRateTempo, !sync);
        pGraphics->HideControl(kParamSequencerRateHz, sync);
        pGraphics->SetAllControlsDirty();
        });
    pGraphics->AttachControl(new SequencerControl<kNumSeqSteps>(sequencerDisplayPanel.GetReducedFromRight(dialSize).GetPadded(-20.f, -10.f, -20.f, -16.f), ""),
      kCtrlTagSequencer, "Sequencer")->SetActionFunction([pGraphics](IControl* pCaller) {
      // Get all values
      /* NB: In the SnapToMouse() function, values for sliders which fall between MouseMove events are interpolated using the preceeding and following slider.
      If this function were to save only the last slider hit (as in the iPlugChunks example at the time of this comment), sliders with interpolated values will not
      have their values stored and serialized. */
      for(int i(0); i < pCaller->NVals(); ++i)
      {
        // Bundle the slider index with the slider value
        ParamTuple data = { i, pCaller->GetValue(i) };
        // Send the new value to be evaluated in the OnMessage function
        pGraphics->GetDelegate()->SendArbitraryMsgFromUI(kMsgSeqSliderChanged, kCtrlTagSequencer, sizeof(ParamTuple), &data);
      }
      });
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(5, 2, 7, 3).GetPadded(-15.f, -2.f, -15.f, -2.f), kParamSequencerRestart, "", toggleStyle, "Free Run", "Restart"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new ModSliderControl(sequencerDisplayPanel.GetGridCell(0, 5, 1, 6).GetVPadded(-18.f).GetVShifted(-5.f).GetMidHPadded(20.f), kNoParameter, "Depth"), kCtrlTagSequencerDepth, "Sequencer");
    pGraphics->AttachControl(new IVRadioButtonControl(sequencerParamsPanel.SubRectVertical(2, 1).SubRectHorizontal(3, 1).GetReducedFromBottom(sequencerParamsPanel.H() / 7.f).GetHPadded(-15.f), kParamSequencerStepMode, { "Auto", "Gate", "Trig." }, "Step Mode", TABLITSA_STYLE.WithLabelText(TABLITSA_TEXT.WithSize(13.f)).WithColor(EVColor::kX1, IColor(150, 155, 155, 155))), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(6, 1, 7, 3).GetPadded(-15.f, -2.f, -15.f, -2.f), kNoParameter, "", toggleStyle, "Cont.", "Quant."), kCtrlTagSequencerQuant, "SequencerParams")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        mSequencerIsQuantized = pControl->GetValue();
        int nSteps = mSequencerIsQuantized > 0.5 ? 12 : 0;
        dynamic_cast<IVTrackControlBase*>(pGraphics->GetControlWithTag(kCtrlTagSequencer))->SetNSteps(nSteps);
      });
    pGraphics->AttachControl(new IVButtonControl(sequencerParamsPanel.GetGridCell(6, 2, 7, 3).GetPadded(-15.f, -2.f, -15.f, -2.f), stepRandomizeAction, "Randomize", pushButtonStyle.WithLabelText(pushButtonStyle.labelText.WithSize(12.f))), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new TablitsaVGroupControl("Sequencer", "SequencerParams", 2.f, 15.f, 2.f, 5.f));

    
#ifndef VST3_API
    pGraphics->AttachControl(new IVButtonControl(b.GetFromBottom(KEYBOARD_HEIGHT + 15.f).GetFromRight(170.f).GetFromTop(30.).GetVPadded(-3.f), SplashClickActionFunc,
      "Show/Hide Keyboard", TABLITSA_STYLE.WithColor(kFG, COLOR_WHITE).WithLabelText({15.f, EVAlign::Middle})))->SetAnimationEndActionFunction(
      [pGraphics, globalControls](IControl* pCaller) {
        static bool hide = false;
        pGraphics->GetControlWithTag(kCtrlTagKeyboard)->Hide(hide = !hide);
        pGraphics->Resize((int)PLUG_WIDTH, hide ? (int)PLUG_HEIGHT : (float)PLUG_HEIGHT + (float)KEYBOARD_HEIGHT, pGraphics->GetDrawScale());
    });
    pGraphics->SetQwertyMidiKeyHandlerFunc([pGraphics](const IMidiMsg& msg) {
                                              dynamic_cast<IVKeyboardControl*>(pGraphics->GetControlWithTag(kCtrlTagKeyboard))->SetNoteFromMidi(msg.NoteNumber(), msg.StatusMsg() == IMidiMsg::kNoteOn);
                                           });
#endif
#ifdef OS_IOS
    if(!IsAuv3AppExtension())
    {
      pGraphics->AttachControl(new IVButtonControl(b.GetFromTRHC(100, 100), [pGraphics](IControl* pCaller) {
                               dynamic_cast<IGraphicsIOS*>(pGraphics)->LaunchBluetoothMidiDialog(pCaller->GetRECT().L, pCaller->GetRECT().MH());
                               SplashClickActionFunc(pCaller);
                             }, "BTMIDI"));
    }
#endif

    pGraphics->AttachPopupMenuControl();

#ifndef VST3_API
    pGraphics->Resize(PLUG_WIDTH, PLUG_HEIGHT + (int)KEYBOARD_HEIGHT, 0.8, true);
#else
    pGraphics->Resize(PLUG_WIDTH, PLUG_HEIGHT, 0.8, true);
#endif
  };
#endif
}

bool Tablitsa::SerializeState(IByteChunk& chunk) const
{
  bool savedOK = true;
  int version{ GetPluginVersion(false) };
  savedOK &= chunk.Put(&version) > 0;
  savedOK &= chunk.Put(&mPresetID) > 0;
  for (int i{ 0 }; i < PRESET_NAME_CHAR_LENGTH; ++i)
    savedOK &= chunk.Put(mPresetName + i) > 0;
  
  for (int i = 0; i < kNumSeqSteps; i++)
  {
    savedOK &= chunk.Put(&mDSP.mSeqSteps[i]) > 0;
  }
  savedOK &= chunk.Put(&mDSP.mLoadedWavetables[0]) > 0;
  savedOK &= chunk.Put(&mDSP.mLoadedWavetables[1]) > 0;
  savedOK &= chunk.Put(&mSequencerIsQuantized) > 0;

  // Save effect panel states
  for (int i{ 0 }; i < TABLITSA_MAX_VOICE_EFFECTS; ++i)
  {
    savedOK &= chunk.Put(&mVoiceEffectSlots[i]) > 0;
  }
  for (int i{ 0 }; i < TABLITSA_MAX_MASTER_EFFECTS; ++i)
  {
    savedOK &= chunk.Put(&mMasterEffectSlots[i]) > 0;
    savedOK &= chunk.Put(&mDelayTempoSync[i]) > 0;
  }
  savedOK &= SerializeParams(chunk);
  return savedOK; // must remember to call SerializeParams at the end
}

// this over-ridden method is called when the host is trying to load the plug-in state and you need to unpack the data into your algorithm
int Tablitsa::UnserializeState(const IByteChunk& chunk, int startPos)
{
  startPos = CheckVersion(chunk);
  startPos = chunk.Get(&mPresetID, startPos);
  for (int i{ 0 }; i < PRESET_NAME_CHAR_LENGTH; ++i)
    startPos = chunk.Get(mPresetName + i, startPos);

  double v = 0.;
  // unserialize the steps state before unserializing the regular params
  for (int i = 0; i < kNumSeqSteps; i++)
  {
    startPos = chunk.Get(&v, startPos);
    mDSP.mSeqSteps[i] = v;
  }
  // Wavetable 1
  startPos = chunk.Get(&v, startPos);
  mDSP.mLoadedWavetables[0] = v;
  // Wavetable 2
  startPos = chunk.Get(&v, startPos);
  mDSP.mLoadedWavetables[1] = v;

  // Notify user if improper wavetable values are loaded
  if (mDSP.mLoadedWavetables[0] - static_cast<int>(mDSP.mLoadedWavetables[0]) > 0.01 || mDSP.mLoadedWavetables[1] - static_cast<int>(mDSP.mLoadedWavetables[1]) > 0.01)
  {
    mDSP.mLoadedWavetables[0] = 1.;
    mDSP.mLoadedWavetables[1] = 2.;
    bool abort{ ShowLoadErrorMessageBox() };
    if (abort)
    {
      if (mStateBackup.Size() == 0)
      {
        LoadDefaultState();
        return -1;
      }
      else
        return UnserializeState(mStateBackup, 0);
    }
      
  }
  if (!LoadWavetables())
    return -1;

  // Sequencer quantization
  startPos = chunk.Get(&v, startPos);
  mSequencerIsQuantized = v;
  
  // Load effect panel states
  for (int i{ 0 }; i < TABLITSA_MAX_VOICE_EFFECTS; ++i)
  {
    startPos = chunk.Get(&mVoiceEffectSlots[i], startPos);

    switch (mVoiceEffectSlots[i])
    {
    case kDistortionEffect:
    {
      GetParam(kParamVoiceEffect1Param1 + kNumVoiceEffectParams * i)->InitDouble("Distortion Param 1", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param2 + kNumVoiceEffectParams * i)->InitDouble("Distortion Param 2", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param3 + kNumVoiceEffectParams * i)->InitDouble("Distortion Param 3", 0., TABLITSA_DISTORTION_FREQ_LOW, TABLITSA_DISTORTION_FREQ_HIGH, 0.01);
      GetParam(kParamVoiceEffect1Param4 + kNumVoiceEffectParams * i)->InitDouble("Distortion Param 4", 0., 0., 1., 0.01);
      break;
    }
    case kSampleAndHoldEffect:
    {
      GetParam(kParamVoiceEffect1Param1 + kNumVoiceEffectParams * i)->InitDouble("Sample & Hold Rate", 10., TABLITSA_SAH_MIN_MS, TABLITSA_SAH_MAX_MS, 0.01, "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
      GetParam(kParamVoiceEffect1Param2 + kNumVoiceEffectParams * i)->InitDouble("Sample & Hold Decay", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param3 + kNumVoiceEffectParams * i)->InitDouble("Sample & Hold Noise", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param4 + kNumVoiceEffectParams * i)->InitDouble("Sample & Hold Mix", 0., 0., 1., 0.01);
      break;
    }
    case kTexturizerEffect:
    {
      GetParam(kParamVoiceEffect1Param1 + kNumVoiceEffectParams * i)->InitDouble("Texturizer Cutoff", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param2 + kNumVoiceEffectParams * i)->InitDouble("Texturizer Drive", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param3 + kNumVoiceEffectParams * i)->InitDouble("Texturizer Res. Freq.", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param4 + kNumVoiceEffectParams * i)->InitDouble("Texturizer Res. Amt.", 0., 0., 1., 0.01);
      break;
    }
    case kCMEffect:
    {
      GetParam(kParamVoiceEffect1Param1 + kNumVoiceEffectParams * i)->InitDouble("SuperRing Mod Amount", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param2 + kNumVoiceEffectParams * i)->InitDouble("SuperRing Pitch", 0., -2., 2., 0.01, "8v.");
      GetParam(kParamVoiceEffect1Param3 + kNumVoiceEffectParams * i)->InitDouble("SuperRing Cutoff", 0., 0., 1., 0.01);
      GetParam(kParamVoiceEffect1Param4 + kNumVoiceEffectParams * i)->InitDouble("SuperRing Mix", 0., 0., 1., 0.01);
      break;
    }
    }
    // Send message to change effect BEFORE parameters are loaded, so that the effect receives the proper parameter values
    SendArbitraryMsgFromUI(kMsgVoiceEffect1Changed + i, kNoTag, sizeof(mVoiceEffectSlots[i]), reinterpret_cast<void*>(&mVoiceEffectSlots[i]));
  }
  for (int i{ 0 }; i < TABLITSA_MAX_MASTER_EFFECTS; ++i)
  {
    // Get the effect ID for each slot
    startPos = chunk.Get(&mMasterEffectSlots[i], startPos);
    startPos = chunk.Get(&mDelayTempoSync[i], startPos);
    // Give the params the proper ranges
    switch (mMasterEffectSlots[i])
    {
    case kDelayEffect:
    {
      if (mDelayTempoSync[i])
      {
        GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitEnum("Delay L (beats)", DelayEffect<sample>::k8th, { DELAY_TEMPODIV_VALIST });
        GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitEnum("Delay R (beats)", DelayEffect<sample>::k8th, { DELAY_TEMPODIV_VALIST });
      }
      else
      {
        GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitDouble("Delay L (ms)", 100., 1., TABLITSA_MAX_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
        GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitDouble("Delay R (ms)", 100., 1., TABLITSA_MAX_DELAY_MS, 1., "ms", IParam::kFlagsNone, "Effect", IParam::ShapePowCurve(3.));
      }
      GetParam(kParamMasterEffect1Param3 + kNumMasterEffectParams * i)->InitDouble("Delay Feedback", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param4 + kNumMasterEffectParams * i)->InitDouble("Delay Mix", 0., 0., 1., 0.01);
      break;
    }
    case kEQEffect:
    {
      GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitDouble("EQ Low Gain ", 1., 0., 2., 0.01);
      GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitDouble("EQ Mid Gain", 1., 0., 2., 0.01);
      GetParam(kParamMasterEffect1Param3 + kNumMasterEffectParams * i)->InitDouble("EQ Mid Freq", 0.5, 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param4 + kNumMasterEffectParams * i)->InitDouble("EQ High Gain", 1., 0., 2., 0.01);
      break;
    }
    case kReverbEffect:
    {
      GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 1", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 2", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param3 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 3", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param4 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 4", 0., 0., 1., 0.01);
      break;
    }
    case kReverb2Effect:
    {
      GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 1", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 2", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param3 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 3", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param4 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 4", 0., 0., 1., 0.01);
      break;
    }
    default:
    {
      GetParam(kParamMasterEffect1Param1 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 1", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param2 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 2", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param3 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 3", 0., 0., 1., 0.01);
      GetParam(kParamMasterEffect1Param4 + kNumMasterEffectParams * i)->InitDouble("Master Effect Parameter 4", 0., 0., 1., 0.01);
      break;
    }
    }
    SendArbitraryMsgFromUI(kMsgMasterEffect1Changed + i, kNoTag, sizeof(mMasterEffectSlots[i]), reinterpret_cast<void*>(&mMasterEffectSlots[i]));

    // Tempo Sync needs to be set before time values are sent to the DSP thread, but after the effect is loaded
    ENTER_PARAMS_MUTEX;
    GetParam(kParamMasterEffect1Param5 + kNumMasterEffectParams * i)->Set(mDelayTempoSync[i] ? 1. : 0.);
    OnParamChange(kParamMasterEffect1Param5 + kNumMasterEffectParams * i);
    LEAVE_PARAMS_MUTEX;
  }

  // must remember to call UnserializeParams at the end
  return UnserializeParams(chunk, startPos);
}

void Tablitsa::LoadDefaultState()
{
  // unserialize the steps state before unserializing the regular params
  for (int i = 0; i < kNumSeqSteps; i++)
  {
    mDSP.mSeqSteps[i] = 0.;
  }
  // Wavetable 1
  mDSP.mLoadedWavetables[0] = 1.;
  // Wavetable 2
  mDSP.mLoadedWavetables[1] = 2.;
  LoadWavetables();
  // Sequencer quantization
  mSequencerIsQuantized = 0.;

  for (int i{ 0 }; i < TABLITSA_MAX_VOICE_EFFECTS; ++i)
  {
    SetVoiceFXSlot(i, kNoVoiceEffect);
    SendArbitraryMsgFromUI(kMsgVoiceEffect1Changed + i, kNoTag, sizeof(mVoiceEffectSlots[i]), reinterpret_cast<void*>(&mVoiceEffectSlots[i]));
  }
  for (int i{ 0 }; i < TABLITSA_MAX_MASTER_EFFECTS; ++i)
  {
    SetMasterFXSlot(i, kNoMasterEffect);
    SendArbitraryMsgFromUI(kMsgMasterEffect1Changed + i, kNoTag, sizeof(mMasterEffectSlots[i]), reinterpret_cast<void*>(&mMasterEffectSlots[i]));
  }

  DefaultParamValues();
  OnParamReset(kDelegate);

  UpdateUIControls();
  dynamic_cast<DropdownListControl*>(GetUI()->GetControlWithTag(kCtrlTagDefaultPresetList))->SetCurrentIdx(std::find(PRESET_ID_LIST, PRESET_ID_LIST + N_PRESETS, "Default") - PRESET_ID_LIST, false);
}

int Tablitsa::CheckVersion(const IByteChunk& presetData)
{
  int presetVersion{};
  int bytePos = presetData.Get(&presetVersion, 0);
  if (presetVersion == GetPluginVersion(false))
    return bytePos;
  else
    return 0;
}

void Tablitsa::OnUIOpen()
{
  UpdateUIControls();
}

void Tablitsa::UpdateUIControls()
{
  //SendCurrentParamValuesFromDelegate();
  dynamic_cast<DropdownListControl*>(GetUI()->GetControlWithTag(kCtrlTagDefaultPresetList))->SetCustomStr(mPresetName);
  // Sequencer
  auto* pSequencer = GetUI()->GetControlWithTag(kCtrlTagSequencer);

  for (int i = 0; i < kNumSeqSteps; i++)
  {
    pSequencer->SetValue(mDSP.mSeqSteps[i], i);
  }

  // Wavetables
  auto* pPeriodicTable = dynamic_cast<PeriodicTable*>(GetUI()->GetControlWithTag(kCtrlTagPeriodicTable));
  pPeriodicTable->SetSelectedElement(static_cast<int>(mDSP.mLoadedWavetables[0]), 0);
  pPeriodicTable->SetSelectedElement(static_cast<int>(mDSP.mLoadedWavetables[1]), 1);
  GetUI()->GetControlWithTag(kCtrlTagSequencerQuant)->SetValue(mSequencerIsQuantized);

  // LFOPlots
  dynamic_cast<ModPlotControl*>(GetUI()->GetControlWithTag(kCtrlTagLFO1Plot))->SetPlotTable(static_cast<const double*>(mDSP.mGlobalLFO1.GetCurrentTable()));
  dynamic_cast<ModPlotControl*>(GetUI()->GetControlWithTag(kCtrlTagLFO2Plot))->SetPlotTable(static_cast<const double*>(mDSP.mGlobalLFO2.GetCurrentTable()));

  SendCurrentParamValuesFromDelegate();
  GetUI()->SetAllControlsDirty();

  for (auto i{ 0 }; i < std::size(kStartupTriggerControls); ++i)
  {
    IControl* pControl = GetUI()->GetControlWithTag(kStartupTriggerControls[i]);
    if (!pControl->IsHidden())
      pControl->SetDirty(true);
  }

  // Set up master effect controls
  DropdownListControl* effectList = dynamic_cast<DropdownListControl*>(GetUI()->GetControlWithTag(kCtrlTagMasterEffectsList));
  // Set up effect slots (without reseting the parameters!)
  for (int i{ TABLITSA_MAX_MASTER_EFFECTS - 1}; i >= 0; --i)
  {
    int effectIdx = mMasterEffectSlots[i];
    effectList->SetCurrentIdx(effectIdx);
  }
  effectList->SetCurrentIdx(static_cast<int>(GetUI()->GetControlWithTag(kCtrlTagMasterEffectsSwitch)->GetValue() * (TABLITSA_MAX_MASTER_EFFECTS - 1)));
  IControl* slotSwitch = GetUI()->GetControlWithTag(kCtrlTagMasterEffectsSwitch);
  slotSwitch->SetValueToDefault();
  slotSwitch->SetDirty(true);

  // Set up voice effect controls
  effectList = dynamic_cast<DropdownListControl*>(GetUI()->GetControlWithTag(kCtrlTagVoiceEffectsList));
  // Set up effect slots (without reseting the parameters!)
  for (int i{ TABLITSA_MAX_VOICE_EFFECTS - 1 }; i >= 0; --i)
  {
    int effectIdx = mVoiceEffectSlots[i];
    effectList->SetCurrentIdx(effectIdx);
  }
  effectList->SetCurrentIdx(static_cast<int>(GetUI()->GetControlWithTag(kCtrlTagVoiceEffectsSwitch)->GetValue() * (TABLITSA_MAX_VOICE_EFFECTS - 1)));
  slotSwitch = GetUI()->GetControlWithTag(kCtrlTagVoiceEffectsSwitch);
  slotSwitch->SetValueToDefault();
  slotSwitch->SetDirty(true);
}

/* implement this and return true to trigger your custom about box, when someone clicks about in the menu of a standalone app or VST3 plugin */
bool Tablitsa::OnHostRequestingAboutBox() { return false; } // See IPlugAPP_dialog.cpp

/* implement this and return true to trigger your custom help info, when someone clicks help in the menu of a standalone app or VST3 plugin */
bool Tablitsa::OnHostRequestingProductHelp() { return false; }

int Tablitsa::GetActiveModIdx() const
{
  return mActiveModIdx;
}

void Tablitsa::RefreshEffectBankControl()
{
  dynamic_cast<TablitsaEffectBankControl*>(GetUI()->GetControlWithTag(kCtrlTagEffectBank))->TabChanged(mCurrentEffectsTab, false); // Prevent controls in hidden panels from being shown
}

#if IPLUG_DSP
void Tablitsa::ProcessBlock(sample** inputs, sample** outputs, int nFrames)
{
  mDSP.ProcessBlock(nullptr, outputs, 2, nFrames, mTimeInfo.mPPQPos, mTimeInfo.mTransportIsRunning, mTimeInfo.mTempo);
  mMeterSender.ProcessBlock(outputs, nFrames, kCtrlTagMeter);
//  mLFO1VisSender.PushData({ kCtrlTagLFO1Vis, {float(mDSP.mGlobalLFO1.GetLastOutput())}});
//  mLFO2VisSender.PushData({ kCtrlTagLFO2Vis, {float(TablitsaDSP<sample>::mGlobalLFO2.GetLastOutput())} });
}

void Tablitsa::OnIdle()
{
  mMeterSender.TransmitData(*this);
//  mLFO1VisSender.TransmitData(*this);
//  mLFO2VisSender.TransmitData(*this);

  // Send time signature to DSP class
  GetTimeSig(mDSP.mTSNum, mDSP.mTSDenom);
  mDSP.mTempo = GetTempo();

  // Only allow the sequencer to retrigger gates in mono mode
  if (mDSP.mActiveSequencer && mDSP.mMono && mDSP.mGlobalSequencer.GetMode() == Sequencer<sample>::kTriggerGate)
  {
    bool retrigger = mDSP.mActiveSequencer->GetTrigger();
    if (retrigger)
    {
      mDSP.ForEachVoice([](TablitsaDSP<sample>::Voice& voice) {
        if (voice.GetBusy() && !voice.JustTriggered())
          voice.Retrigger();
        });
    }
  }

  if (GetUI())
  {
    // Update Sequencer Position
    if (mDSP.GetSequencerStep() != mDSP.mSeqPos)
    {
      mDSP.mSeqPos = mDSP.GetSequencerStep();
      SendControlMsgFromDelegate(kCtrlTagSequencer, IVMultiSliderControl<>::kMsgTagSetHighlight, sizeof(int), &mDSP.mSeqPos);
    }
  }
}

void Tablitsa::OnReset()
{
  mDSP.Reset(GetSampleRate(), GetBlockSize());
}

void Tablitsa::ProcessMidiMsg(const IMidiMsg& msg)
{
  TRACE;
  
  int status = msg.StatusMsg();
  
  switch (status)
  {
    case IMidiMsg::kNoteOn:
      if (mDSP.mTableLoading) // Don't accept new notes when swapping wavetables
        return;
    case IMidiMsg::kNoteOff:
    case IMidiMsg::kPolyAftertouch:
    case IMidiMsg::kControlChange:
    case IMidiMsg::kProgramChange:
    case IMidiMsg::kChannelAftertouch:
    case IMidiMsg::kPitchWheel:
    {
      goto handle;
    }
    default:
      return;
  }
  
handle:
  mDSP.ProcessMidiMsg(msg);
  SendMidiMsg(msg);

  // Send additional messages for unison voices
  mDSP.ResetDetune();
  for (int i{ 1 }; i < mDSP.mUnisonVoices; ++i)
  {
    IMidiMsg copy_msg = IMidiMsg(msg);
    if (status == kNoteOnAction)
      copy_msg.MakeNoteOnMsg(msg.NoteNumber(), msg.Velocity(), msg.mOffset + 1, msg.Channel());
    else if (status == kNoteOffAction)
      copy_msg.MakeNoteOffMsg(msg.NoteNumber(), msg.mOffset + 1, msg.Channel());
    mDSP.ProcessMidiMsg(msg);
  }
}

void Tablitsa::OnParamChange(int paramIdx)
{
  ENTER_PARAMS_MUTEX
  mDSP.SetParam(paramIdx, GetParam(paramIdx)->Value());
  LEAVE_PARAMS_MUTEX
}

bool Tablitsa::OnMessage(int msgTag, int ctrlTag, int dataSize, const void* pData)
{
  // TODO: switch statement
  if(ctrlTag == kCtrlTagBender && msgTag == IWheelControl::kMessageTagSetPitchBendRange)
  {
    const int bendRange = *static_cast<const int*>(pData);
    mDSP.mSynth.SetPitchBendRange(bendRange);
  }
  else
  {
    switch (msgTag)
    {
    case kMsgSavePreset:
    {
      auto* backup = reinterpret_cast<const bool*>(pData);
      std::string filePath;

      bool useSystemDialog{ (bool)LOAD_WITH_DIALOG };
      if (useSystemDialog)
      {
#pragma mark - Windows
        WDL_String fileStr("NewPreset");
        WDL_String pathStr(GetDataPath("\\Tablitsa\\presets\\").c_str());
        GetUI()->PromptForFile(fileStr, pathStr, EFileAction::Save, "TPST");
        // Get and store preset name
        std::string path = fileStr.Get();
        size_t nameStart{ path.find_last_of("\\") + 1 };
        size_t nameEnd{ path.find_last_of(".") };
        strcpy(mPresetName, path.substr(nameStart, nameEnd - nameStart).c_str());
        // Abort if no filename was submitted
        filePath = fileStr.Get();
        if (filePath == "\0")
          break;
      }
      else
      {
        filePath = GetDataPath("\\Tablitsa\\presets\\") + "TEST.TPST";
      }
      IByteChunk newPreset;
      SerializeState(newPreset);
      SavePreset(newPreset, filePath.c_str(), *backup);

      break;
    }
    case kMsgLoadPreset:
    {
      auto* backup = reinterpret_cast<const bool*>(pData);
      std::string filePath;
      if (*backup)
      {
        filePath = GetDataPath("\\Tablitsa\\presets\\") + "_tablitsa.bkp";
      }
      else if ((bool)LOAD_WITH_DIALOG)
      {
        WDL_String fileStr("UserPreset");
        WDL_String pathStr(GetDataPath("\\Tablitsa\\presets\\").c_str());
        GetUI()->PromptForFile(fileStr, pathStr, EFileAction::Open, "TPST");
        filePath = fileStr.Get();
        if (filePath == "\0")
          break;
      }
      else
        filePath = GetDataPath("\\Tablitsa\\presets\\") + "TEST.TPST";

      IByteChunk presetData = LoadPreset(filePath.c_str(), *backup);
      if (presetData.Size() == 0)
        break;
      UnserializeState(presetData, 0);
      UpdateUIControls();
      break;
    }
    case kMsgWavetable1Changed:
    {
      PeriodicTable* pTableControl = dynamic_cast<PeriodicTable*>(GetUI()->GetControlWithTag(kCtrlTagPeriodicTable));
      pTableControl->SetTableLoading(true, 0);
      const double wtIdx = *static_cast<const double*>(pData); //Semi-Normalized value (between 1/118 and 1)
      mDSP.mLoadedWavetables[0] = wtIdx * 118.; // Corresponds to the atomic number (not the 0-indexed vector item)
      try
      {
        mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[0]) - 1, 0);
      }
      catch (...)
      {
        GetUI()->ShowMessageBox("Wavetable file not found. Check AppData directory.", "Error", EMsgBoxType::kMB_OK);
      }
      pTableControl->SetTableLoading(false, 0);
      break;
    }
    case kMsgWavetable2Changed:
    {
      PeriodicTable* pTableControl = dynamic_cast<PeriodicTable*>(GetUI()->GetControlWithTag(kCtrlTagPeriodicTable));
      pTableControl->SetTableLoading(true, 1);
      const double wtIdx = *static_cast<const double*>(pData);
      mDSP.mLoadedWavetables[1] = wtIdx * 118.;
      mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[1]) - 1, 1);
      pTableControl->SetTableLoading(false, 1);
      break;
    }
    case kMsgUpdateLFO1Plot:
    {
      dynamic_cast<ModPlotControl*>(GetUI()->GetControlWithTag(kCtrlTagLFO1Plot))->SetPlotTable(static_cast<const double*>(mDSP.mGlobalLFO1.GetCurrentTable()));
      break;
    }
    case kMsgUpdateLFO2Plot:
    {
      dynamic_cast<ModPlotControl*>(GetUI()->GetControlWithTag(kCtrlTagLFO2Plot))->SetPlotTable(static_cast<const double*>(mDSP.mGlobalLFO2.GetCurrentTable()));
      break;
    }
    case kMsgSeqSliderChanged:
    {
      // Read incoming data as a tuple
      auto* pConverted = reinterpret_cast<const ParamTuple*>(pData);
      mDSP.mSeqSteps[pConverted->idx] = pConverted->value;
      break;
    }
    case kMsgRandomizeSequencer:
    {
      auto* pSequencer = dynamic_cast<SequencerControl<kNumSeqSteps>*>(GetUI()->GetControlWithTag(kCtrlTagSequencer));
      for (auto i{ 0 }; i < kNumSeqSteps; ++i)
      {
        mDSP.mSeqSteps[i] = static_cast<double>(std::rand() % 1000) / 1000.;
        pSequencer->SetValue(mDSP.mSeqSteps[i], i);
      }
      pSequencer->SetDirty(true);
      break;
    }
    case kMsgVoiceEffect1Changed:
    case kMsgVoiceEffect2Changed:
    case kMsgVoiceEffect3Changed:
    {
      //std::lock_guard<std::mutex> lg(mDSP.mProcMutex);
      int effectIdx = *reinterpret_cast<const int*>(pData);
      switch (effectIdx)
      {
      case kDistortionEffect:
        mDSP.ForEachVoice([msgTag, this](TablitsaDSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kDistortionEffect);
          });
        break;
      case kSampleAndHoldEffect:
        mDSP.ForEachVoice([msgTag, this](TablitsaDSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kSampleAndHoldEffect);
          });
        break;
      case kTexturizerEffect:
        mDSP.ForEachVoice([msgTag, this](TablitsaDSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kTexturizerEffect);
          });
        break;
      case kCMEffect:
        mDSP.ForEachVoice([msgTag, this](TablitsaDSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kCMEffect);
          });
        break;
      case kWaveFolderEffect:
        mDSP.ForEachVoice([msgTag, this](TablitsaDSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kWaveFolderEffect);
          });
        break;
      default:
        mDSP.ForEachVoice([msgTag, this](TablitsaDSP<sample>::Voice& voice) {
          voice.SetEffect(msgTag - kMsgVoiceEffect1Changed, kNoVoiceEffect);
          });
        break;
      }
      SendCurrentParamValuesFromDelegate();
      break;
    }
    case kMsgMasterEffect1Changed:
    case kMsgMasterEffect2Changed:
    case kMsgMasterEffect3Changed:
    {
      std::lock_guard<std::mutex> lg(mDSP.mEffectMutex);
      int effectIdx = *reinterpret_cast<const int*>(pData);
      delete mDSP.mEffects[msgTag - kMsgMasterEffect1Changed];
      switch (effectIdx)
      {
      case kChorusEffect:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new ChorusEffect<iplug::sample>(GetSampleRate());
        break;
      case kDelayEffect:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new DelayEffect<iplug::sample>(GetSampleRate(), &mDSP.mGlobalMetronome);
        break;
      case kEQEffect:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new EQ3Effect<iplug::sample>(GetSampleRate());
        break;
      case kReverbEffect:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new ReverbEffect<iplug::sample>(GetSampleRate());
        break;
      case kReverb2Effect:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new Reverb2Effect<iplug::sample>(GetSampleRate());
        break;
      default:
        mDSP.mEffects[msgTag - kMsgMasterEffect1Changed] = new Effect<iplug::sample>(GetSampleRate());
        break;
      }
      SendCurrentParamValuesFromDelegate();
      break;
    }
    default:
      break;
    }
  }
  return false;
}

void Tablitsa::SetActiveModIdx(int idx)
{
  int prev = mActiveModIdx;
  mActiveModIdx = idx;
  if (idx > -1)
    dynamic_cast<TablitsaIVModKnobControl*>(GetUI()->GetControl(idx))->ColorSwap();
  if (prev > -1)
    dynamic_cast<TablitsaIVModKnobControl*>(GetUI()->GetControl(prev))->ColorSwap();
}

bool Tablitsa::LoadWavetables()
{
  for (int i{ 0 }; i < 2; ++i)
  {
    try
    {
      mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[i]) - 1, i);
    }
    catch (const std::out_of_range& e)
    {
      mDSP.mLoadedWavetables[i] = static_cast<double>(i + 1);
      dynamic_cast<PeriodicTable*>(GetUI()->GetControlWithTag(kCtrlTagPeriodicTable))->SetSelectedElement(static_cast<int>(mDSP.mLoadedWavetables[i]), i);
      mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[i]) - 1, i);

      bool abort{ ShowLoadErrorMessageBox() };
      if (abort)
      {
        LoadDefaultState();
        return false;
      }
    }
  }
  return true;
}

void Tablitsa::SavePreset(IByteChunk& byteData, const char* filename, bool isBackup)
{
  assert(byteData.Size() != 0);
  if (isBackup)
  {
    mStateBackup = byteData;
    return;
  }

  std::ofstream presetFile(filename, std::ios::binary | std::ios::trunc);

  if (!presetFile.is_open())
  {
    GetUI()->ShowMessageBox("Cannot access selected file path.", "Error", EMsgBoxType::kMB_OK);
    presetFile.close();
    return;
  }

  //char* allBytes = new char[byteData.Size()];
  std::vector<char> allBytes(byteData.Size());
  byteData.GetBytes(&allBytes[0], byteData.Size(), 0);
  presetFile.write(&allBytes[0], byteData.Size());
  presetFile.close();
  //delete[] allBytes;
}

IByteChunk Tablitsa::LoadPreset(const char* filename, bool isBackup)
{
  if (isBackup)
    return mStateBackup;
  else
  {
    mStateBackup.Clear();
    SerializeState(mStateBackup);
  }

  IByteChunk presetData;

  std::vector<char> fileData = ReadAllBytes(filename);

  if (fileData.size() == 0)
  {
    GetUI()->ShowMessageBox("Preset file not found or invalid format", "Error", EMsgBoxType::kMB_OK);
    return presetData;
  }
  
  presetData.PutBytes(&fileData[0], fileData.size());

  return presetData;
}

bool Tablitsa::ShowLoadErrorMessageBox()
{
  bool abort;
  GetUI()->ShowMessageBox("Invalid wavetable codes. Preset file may be corrupted or incompatible with the plugin version. Continue loading preset?", "Warning", EMsgBoxType::kMB_YESNO,
    [&abort](EMsgBoxResult result) {
      abort = (result == EMsgBoxResult::kNO);
    });
  return abort;
}
#endif

std::string GetDataPath(char* appendPath)
{
  USES_CONVERSION;
  TCHAR szPath[MAX_PATH];
  if (SUCCEEDED(SHGetFolderPath(NULL, CSIDL_APPDATA, NULL, 0, szPath)))
  {
    PathAppend(szPath, A2T(appendPath));
  }
#if VST3_API
  // In the standalone app builds, `TCHAR` is type `char`. In VST3 builds, it's `wchar_t`, requiring a different conversion method.
  std::wstring wpath(szPath);
  std::string path = std::wstring_convert<std::codecvt_utf8<wchar_t>>().to_bytes(wpath);
#else
  std::string path = szPath;
#endif
  return path;
}

std::vector<char> ReadAllBytes(const char* fname)
{
  std::ifstream ifs(fname, std::ios::binary | std::ios::ate);
  std::ifstream::pos_type pos = ifs.tellg();

  if (!ifs.is_open())
  {
    return std::vector<char>(0);
  }

  std::vector<char> allBytes(pos);
  if (pos == 0)
    return allBytes;
  ifs.seekg(0, std::ios::beg);
  ifs.read(&allBytes[0], pos);
  ifs.close();
  return allBytes;
}
