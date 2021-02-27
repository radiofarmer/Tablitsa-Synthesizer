
#include "Tablitsa.h"
#include "IPlug_include_in_plug_src.h"

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
  GetParam(kParamGain)->InitGain("Gain", -24., -72., 0.);
  GetParam(kParamPan)->InitDouble("Pan", 0., -90., 90., 0.5, "°");
  GetParam(kParamNoteGlideTime)->InitMilliseconds("Note Glide Time", 0., 0.0, 3000.);
  GetParam(kParamNoteGlideRate)->InitDouble("Note Glide Rate", 12., 0.5, 128., 0.1, "St/s");
  GetParam(kParamPortamentoMode)->InitBool("Equal Glide Time", true);
  GetParam(kParamMonophonic)->InitBool("Monophonic", false);
  GetParam(kParamUnisonVoices)->InitInt("Unison Voices", 1, 1, 8);
  GetParam(kParamUnisonDetune)->InitDouble("Unison Detune", 0., 0., 12., 0.01, "St.");
  GetParam(kParamUnisonChord)->InitEnum("Chord", 0, { UNISON_CHORD_LIST });
  GetParam(kParamStereoSpread)->InitDouble("Stereo Spread", 0., -180., 180., 0.5, "°");
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
  GetParam(kParamLFO1Shape)->InitEnum("LFO1 Shape", LFO<>::kTriangle, {LFO_SHAPE_VALIST});
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
  GetParam(kParamLFO2Shape)->InitEnum("LFO2 Shape", LFO<>::kTriangle, { LFO_SHAPE_VALIST });
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
  // Wavetable selection
  GetParam(kParamWavetable1)->InitEnum("Wavetable 1", 1, { ELEMENT_NAMES });
  GetParam(kParamWavetable2)->InitEnum("Wavetable 2", 2, { ELEMENT_NAMES });
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
  GetParam(kParamWavetable1Formant)->InitDouble("WT1 Formant", 1., 0., 1., 0.01);
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
  GetParam(kParamWavetable2Formant)->InitDouble("WT2 Formant", 1., 0., 1., 0.01);
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
  GetParam(kParamFilter1Resonance)->InitPercentage("Filter 1 Resonance", 0.1, 0., 100., 1);
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
  GetParam(kParamFilter1Delay)->InitInt("Filter 1 Delay", 8, 1, 512);

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
  GetParam(kParamFilter2Resonance)->InitPercentage("Filter 2 Resonance", 0.1, 0., 100., 1);
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
  GetParam(kParamFilter2DriveVel)->InitDouble("Filter 2 Drive LFO1 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2DriveKTk)->InitDouble("Filter 2 Drive LFO2 Depth", 0., -1., 1., 0.01);
  GetParam(kParamFilter2DriveRnd)->InitDouble("Filter 2 Drive Seq Depth", 0., -1., 1., 0.01);
  // Comb-specific parameters with different ranges (the modulation parameters are common between all filter types, as these always between -1 and 1)
  GetParam(kParamFilter2FF)->InitPercentage("Filter 2 FF");
  GetParam(kParamFilter2FB)->InitPercentage("Filter 2 FB");
  GetParam(kParamFilter2Delay)->InitInt("Filter 2 Delay", 8, 1, 512);

  // Effect Parameters
  GetParam(kParamEffect1Param1)->InitDouble("Effect 1 Parameter 1", 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Effect1", IParam::ShapePowCurve(3.));
  GetParam(kParamEffect1Param2)->InitDouble("Effect 1 Parameter 2", 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Effect1", IParam::ShapePowCurve(3.));
  GetParam(kParamEffect1Param3)->InitPercentage("Effect 1 Parameter 3", 0.);
  GetParam(kParamEffect1Param4)->InitPercentage("Effect 1 Parameter 4", 0.);
  GetParam(kParamEffect1Param5)->InitBool("Effect 1 Parameter 5", false);
  //GetParam(kParamDelayTimeLBeats)->InitEnum("Stereo Delay Time L (tempo)", LFO<>::k8th, { DELAY_TEMPODIV_VALIST });
  //GetParam(kParamDelayTimeRBeats)->InitEnum("Stereo Delay Time R (tempo)", LFO<>::k8th, { DELAY_TEMPODIV_VALIST });

  GetParam(kParamEffect2Param1)->InitDouble("Effect 1 Parameter 1", 10., 1., 20., 0.1, "ms", IParam::kFlagsNone, "Effect1", IParam::ShapePowCurve(3.));
  GetParam(kParamEffect2Param2)->InitPercentage("Effect 1 Parameter 2", 0.);
  GetParam(kParamEffect2Param3)->InitPercentage("Sample-and-Hold Noise");
  GetParam(kParamEffect2Param4)->InitPercentage("Sample-and-Hold Mix", 0.);

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
    IRECT keyboardBounds = b.GetFromBottom((float)KEYBOARD_HEIGHT - plugBorders).GetPadded(-5.f, -25.f, -5.f, -10.f);
    IRECT wheelsBounds = keyboardBounds.ReduceFromLeft(100.f).GetPadded(-10.f);

    // Keyboard
    pGraphics->AttachControl(new IVKeyboardControl(keyboardBounds, 19, 79), kCtrlTagKeyboard);
    // Mod wheels
    pGraphics->AttachControl(new IWheelControl(wheelsBounds.FracRectHorizontal(0.5)), kCtrlTagBender);
    pGraphics->AttachControl(new IWheelControl(wheelsBounds.FracRectHorizontal(0.5, true), IMidiMsg::EControlChangeMsg::kModWheel));
    //    pGraphics->AttachControl(new IVMultiSliderControl<4>(b.GetGridCell(0, 2, 2).GetPadded(-30), "", TABLITSA_STYLE, kParamAttack, EDirection::Vertical, 0.f, 1.f));

#endif
    const ISVG tableImg = pGraphics->LoadSVG(PERIODIC_TABLE_FN);
    const ISVG incrButtonOn = pGraphics->LoadSVG(INCR_BUTTON_ON);
    const ISVG incrButtonOff = pGraphics->LoadSVG(INCR_BUTTON_OFF);
    const ISVG decrButtonOn = pGraphics->LoadSVG(DECR_BUTTON_ON_FN);
    const ISVG decrButtonOff = pGraphics->LoadSVG(DECR_BUTTON_OFF_FN);
    const IBitmap logoImg = pGraphics->LoadBitmap(LOGO_FN);
    const IBitmap dialImg = pGraphics->LoadBitmap(DIAL_FN, 64);

    // Periodic table
    pGraphics->AttachControl(new PeriodicTable(globalControls.GetReducedFromBottom(globalControls.H() * 0.5f).GetFromLeft(720.f).GetPadded(-30.f).GetVShifted(20.f), tableImg, { kParamWavetable1, kParamWavetable2 }),
      kCtrlTagPeriodicTable);
    dynamic_cast<PeriodicTable*>(pGraphics->GetControlWithTag(kCtrlTagPeriodicTable))->LoadValues();
    pGraphics->AttachControl(new IBitmapControl(globalControls.GetReducedFromTop(90.f).GetReducedFromLeft(115.f).GetFromTLHC(280.f, 80.), logoImg));

    // Preset controls
    const IRECT presetPanel = globalControls.GetFromTop(50.f).GetFromLeft(720.f).GetHPadded(-60.f);
    IActionFunction savePresetFunction = IActionFunction([this](IControl* pControl) {
      bool saveBackup = false;
      SendArbitraryMsgFromUI(kMsgSavePreset, kNoTag, sizeof(bool), reinterpret_cast<const void*>(&saveBackup));
      });
    IActionFunction loadPresetFunction = IActionFunction([this](IControl* pControl) {
      bool loadBackup = false;
      SendArbitraryMsgFromUI(kMsgLoadPreset, kNoTag, sizeof(bool), reinterpret_cast<const void*>(&loadBackup));
      });
    IActionFunction resetToDefaultFunction = IActionFunction([this](IControl* pControl) {
      // Save current state
      std::string backupPath = GetDataPath("\\Tablitsa\\presets\\") + "_backup.TPST";
      IByteChunk currentState;
      SerializeState(currentState);
      SavePreset(currentState, backupPath.c_str(), true);
      // Load default state
      std::string defaultPath = GetDataPath("\\Tablitsa\\presets\\") + "_Tablitsa_Default.TPST";
      IByteChunk presetData = LoadPreset(defaultPath.c_str(), false);
      if (presetData.Size() == 0)
        return;
      UnserializeState(presetData, 0);
      UpdateUIControls();
      });
    IActionFunction loadBackupFunction = IActionFunction([this](IControl* pControl) {
      std::string defaultPath = GetDataPath("\\Tablitsa\\presets\\") + "_backup.TPST";
      IByteChunk presetData = LoadPreset(defaultPath.c_str(), false);
      if (presetData.Size() == 0)
        return;
      UnserializeState(presetData, 0);
      UpdateUIControls();
      });
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(0, 0, 2, 4).GetHPadded(-20.f), SplashClickActionFunc, "Save", pushButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(savePresetFunction);
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(0, 1, 2, 4).GetHPadded(-20.f), SplashClickActionFunc, "Load", pushButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(loadPresetFunction);
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(1, 0, 2, 4).GetHPadded(-20.f), SplashClickActionFunc, "Reset", pushButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(resetToDefaultFunction);
    pGraphics->AttachControl(new IVButtonControl(presetPanel.GetGridCell(1, 1, 2, 4).GetHPadded(-20.f), SplashClickActionFunc, "Load Backup", pushButtonStyle), kNoTag, "Presets")->SetAnimationEndActionFunction(loadBackupFunction);
    // Built-in file selection
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
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc1Panel.GetGridCell(1, 1, 2, 3).GetCentredInside(dialSize), kParamWavetable1Bend, "Energy"), kNoTag, "Osc1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc1Panel.GetGridCell(1, 2, 2, 3).GetCentredInside(dialSize), kParamWavetable1Formant, "Mass"), kNoTag, "Osc1");
    pGraphics->AttachControl(new TablitsaVGroupControl("Wavetable 1", "Osc1", -5.f, 15.f, -5.f, 2.f));

    const IRECT osc2Panel = globalControls.GetReducedFromLeft(720.f).GetReducedFromBottom(globalControls.H() / 2.f).GetGridCell(1, 0, 2, 3).GetPadded(-10.f).GetHPadded(-12.f);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc2Panel.GetGridCell(0, 0, 2, 2).GetCentredInside(dialSize), kParamWavetable2Pitch, "Pitch", modKnobStyle, true), kNoTag, "Osc2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc2Panel.GetGridCell(0, 1, 2, 2).GetCentredInside(dialSize), kParamWavetable2Amp, "Amp"), kNoTag, "Osc2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc2Panel.GetGridCell(1, 0, 2, 3).GetCentredInside(dialSize), kParamWavetable2Pos, "Valency"), kNoTag, "Osc2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc2Panel.GetGridCell(1, 1, 2, 3).GetCentredInside(dialSize), kParamWavetable2Bend, "Energy"), kNoTag, "Osc2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(osc2Panel.GetGridCell(1, 2, 2, 3).GetCentredInside(dialSize), kParamWavetable2Formant, "Mass"), kNoTag, "Osc2");
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
    pGraphics->AttachControl(new TablitsaIVKnobControl(unisonPanel.GetGridCell(0, 0, 2, 2).GetCentredInside(dialSize), kParamUnisonVoices, "№ Voices", modKnobStyle), kNoTag, "Unison");
    pGraphics->AttachControl(new TablitsaIVKnobControl(unisonPanel.GetGridCell(1, 0, 2, 2).GetCentredInside(dialSize), kParamUnisonDetune, "Detune", modKnobStyle), kNoTag, "Unison");
    pGraphics->AttachControl(new TablitsaIVKnobControl(unisonPanel.GetGridCell(0, 1, 2, 2).GetCentredInside(dialSize), kParamUnisonChord, "Chord", modKnobStyle), kNoTag, "Unison");
    pGraphics->AttachControl(new TablitsaIVKnobControl(unisonPanel.GetGridCell(1, 1, 2, 2).GetPadded(-5.f).GetCentredInside(dialSize), kParamStereoSpread, "Spread"), kNoTag, "Unison");
    pGraphics->AttachControl(new TablitsaVGroupControl("Unison", "Unison", 0.f, 13.f, 0.f, 0.f));
    pGraphics->AttachControl(new IVLEDMeterControl<2>(masterPanel.GetGridCell(0, 5, 2, 6).GetVShifted(-15.f).GetHShifted(-4.f).GetVPadded(-8.f)), kCtrlTagMeter, "Master");

    /* Effects Panel */
    const IRECT effectsPanel = IRECT(masterPanel.GetGridCell(1, 1, 2, 7).R, masterPanel.GetGridCell(1, 0, 2, 7).T, masterPanel.GetGridCell(1, 6, 2, 7).R, masterPanel.GetGridCell(1, 6, 2, 7).B).GetHPadded(-5.f).GetHShifted(2.f);

    // Effects tab control
    pGraphics->AttachControl(new TablitsaEffectBankControl(effectsPanel, { "Effect 1", "Effect1", "Effect3" }, { "E1", "E2", "E3" }), kCtrlTagEffectBank, "Effects");
    // Effect 1
    auto channelLockFunction = IActionFunction([pGraphics, this](IControl* pControl) {
      // Tag of the other channel's control
      int ctrlTag = pControl->GetTag() == kCtrlTagEffect1Knob1 || pControl->GetTag() == kCtrlTagEffect1Knob1 ? pControl->GetTag() + 1 : pControl->GetTag() - 1;
      pGraphics->GetControlWithTag(ctrlTag)->SetValue(pControl->GetValue());
      GetUI()->GetDelegate()->SendParameterValueFromUI(pGraphics->GetControlWithTag(ctrlTag)->GetParamIdx(), pControl->GetValue());
      pGraphics->GetControlWithTag(ctrlTag)->SetDirty(false);
      });
    // Selector list
    pGraphics->AttachControl(new DropdownListControl(effectsPanel.FracRectHorizontal(0.4f, false).GetFromTop(40.f).GetReducedFromTop(20.f).GetPadded(-2.f), TABLITSA_EFFECTS_LIST), kCtrlTagEffect1Tab, "E1")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        // Swap out effect parameters (and rearrange controls if necessary)
        int effectIdx = dynamic_cast<DropdownListControl*>(pControl)->GetCurrentIndex();
        IControl* knob1 = pGraphics->GetControlWithTag(kCtrlTagEffect1Knob1);
        IControl* knob2 = pGraphics->GetControlWithTag(kCtrlTagEffect1Knob2);
        IControl* knob3 = pGraphics->GetControlWithTag(kCtrlTagEffect1Knob3);
        IControl* knob4 = pGraphics->GetControlWithTag(kCtrlTagEffect1Knob4);
        IControl* toggle1 = pGraphics->GetControlWithTag(kCtrlTagEffect1Toggle1);
        IControl* toggle2 = pGraphics->GetControlWithTag(kCtrlTagEffect1Toggle2);
        switch (effectIdx)
        {
        case kDelayEffect:
        {
          GetParam(kParamEffect1Param1)->InitDouble("Effect 1 Parameter 1", 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Effect1", IParam::ShapePowCurve(3.));
          GetParam(kParamEffect1Param2)->InitDouble("Effect 1 Parameter 2", 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Effect1", IParam::ShapePowCurve(3.));
          pGraphics->HideControl(kParamEffect1Param5, false);
          pGraphics->HideControl(kParamEffect1Param6, true);
          dynamic_cast<IVKnobControl*>(knob3)->SetLabelStr("L");
          dynamic_cast<IVKnobControl*>(knob3)->SetLabelStr("R");
          dynamic_cast<IVKnobControl*>(knob3)->SetLabelStr("Feedback");
          toggle1->Hide(false);
          toggle2->Hide(false);
          knob1->SetDisabled(false);
          knob2->SetDisabled(false);
          knob3->SetDisabled(false);
          knob4->SetDisabled(false);
          toggle1->SetDirty(true);
          toggle2->SetDirty(true);
          break;
        }
        case kSampleAndHoldEffect:
        {
          GetParam(kParamEffect1Param1)->InitDouble("Effect 1 Parameter 1", 10., 1., 20., 0.1, "ms", IParam::kFlagsNone, "Effect1", IParam::ShapePowCurve(3.));
          GetParam(kParamEffect1Param2)->InitDouble("Effect 1 Parameter 2", 0., 1., 100., 1., "%", IParam::kFlagsNone, "Effect1", IParam::ShapePowCurve(3.));
          knob1->SetParamIdx(kParamEffect1Param1);
          pGraphics->HideControl(kParamEffect1Param5, true);
          pGraphics->HideControl(kParamEffect1Param6, true);
          dynamic_cast<IVKnobControl*>(knob1)->SetLabelStr("Rate");
          dynamic_cast<IVKnobControl*>(knob2)->SetLabelStr("Decay");
          dynamic_cast<IVKnobControl*>(knob3)->SetLabelStr("Noise");
          toggle1->Hide(true);
          toggle2->Hide(true);
          knob1->SetDisabled(false);
          knob2->SetDisabled(false);
          knob3->SetDisabled(false);
          knob4->SetDisabled(false);
          break;
        }
        default:
        {
          dynamic_cast<IVKnobControl*>(knob1)->SetLabelStr("");
          dynamic_cast<IVKnobControl*>(knob2)->SetLabelStr("");
          dynamic_cast<IVKnobControl*>(knob3)->SetLabelStr("");
          knob1->SetDisabled(true);
          knob2->SetDisabled(true);
          knob3->SetDisabled(true);
          knob4->SetDisabled(true);
        }
        }
        SendArbitraryMsgFromUI(kMsgEffect1Changed, kNoTag, sizeof(effectIdx), reinterpret_cast<void*>(&effectIdx)); // Effects must be swapped before OnParamChange is called
        knob1->SetDirty(true);
        knob2->SetDirty(true);
        knob3->SetDirty(true);
        knob4->SetDirty(true);
      }
    );
    const IRECT effectsControlsPanel = effectsPanel.GetReducedFromTop(40.f);
    pGraphics->AttachControl(new TablitsaIVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 0, 1, 4).GetCentredInside(dialSize), kParamEffect1Param1, "Left", modKnobStyle), kCtrlTagEffect1Knob1, "E1");
    pGraphics->AttachControl(new TablitsaIVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 1, 1, 4).GetCentredInside(dialSize), kParamEffect1Param2, "Right", modKnobStyle), kCtrlTagEffect1Knob2, "E1");
    pGraphics->AttachControl(new TablitsaIVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 2, 1, 4).GetCentredInside(dialSize), kParamEffect1Param3, "Feedback", modKnobStyle), kCtrlTagEffect1Knob3, "E1");
    pGraphics->AttachControl(new TablitsaIVKnobControl(effectsControlsPanel.FracRectVertical(2.f / 3, true).GetGridCell(0, 3, 1, 4).GetCentredInside(dialSize), kParamEffect1Param4, "Mix", modKnobStyle), kCtrlTagEffect1Knob4, "E1");
    pGraphics->AttachControl(new IVToggleControl(effectsControlsPanel.FracRectVertical(1.f / 3).GetGridCell(0,  1, 1, 4).GetPadded(-10.f), kParamEffect1Param5, "", toggleStyle, "Time", "Tempo"),
      kCtrlTagEffect1Toggle1, "E1")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        mDelayIsSynced = pControl->GetValue();
        bool isTempoSync = mDelayIsSynced > 0.5;
        if (isTempoSync)
        {
          GetParam(kParamEffect1Param1)->InitEnum("Stereo Delay Time L (tempo)", LFO<>::k8th, { DELAY_TEMPODIV_VALIST });
          GetParam(kParamEffect1Param2)->InitEnum("Stereo Delay Time R (tempo)", LFO<>::k8th, { DELAY_TEMPODIV_VALIST });
        }
        else
        {
          GetParam(kParamEffect1Param1)->InitDouble("Stereo Delay Time L (ms)", 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Delay", IParam::ShapePowCurve(3.));
          GetParam(kParamEffect1Param2)->InitDouble("Stereo Delay Time R (ms)", 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Delay", IParam::ShapePowCurve(3.));
        }
        pGraphics->SetAllControlsDirty();
      }
    );
    pGraphics->AttachControl(new IVToggleControl(effectsControlsPanel.FracRectVertical(1.f / 3).GetGridCell(0, 2, 1, 4).GetPadded(-10.f), kNoParameter, "", toggleStyle, "Stereo", "Mono"), kCtrlTagEffect1Toggle2, "E1")->SetActionFunction(
      [pGraphics, channelLockFunction](IControl* pControl) {
        bool mono = pControl->GetValue() > 0.5;
        if (mono)
        {
          pGraphics->GetControlWithTag(kCtrlTagEffect1Knob1)->SetActionFunction(channelLockFunction);
          pGraphics->GetControlWithTag(kCtrlTagEffect1Knob2)->SetActionFunction(channelLockFunction);
        }
        else
        {
          pGraphics->GetControlWithTag(kCtrlTagEffect1Knob1)->SetActionFunction([pGraphics](IControl* pControl) {});
          pGraphics->GetControlWithTag(kCtrlTagEffect1Knob2)->SetActionFunction([pGraphics](IControl* pControl) {});
        }
        pGraphics->SetAllControlsDirty();
      }
    );
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
    auto setFilterOscillator = IActionFunction([this](IControl* pControl) {
      if (pControl->GetTag() == kCtrlTagFilter1Osc)
      {
        GetUI()->GetControlWithTag(kCtrlTagFilter2Osc)->SetValue(1. - pControl->GetValue());
        GetUI()->GetControlWithTag(kCtrlTagFilter2Osc)->SetDirty(false);
      }
      else
      {
        GetUI()->GetControlWithTag(kCtrlTagFilter1Osc)->SetValue(1. - pControl->GetValue());
        GetUI()->GetControlWithTag(kCtrlTagFilter1Osc)->SetDirty(false);
      }
      SendArbitraryMsgFromUI(kMsgFilterOscChanged, pControl->GetTag());
    });
    /* Filter 1 */
    const IRECT filter1Panel = lowerPanel.GetGridCell(1, 0, 3, 3);
    pGraphics->AttachControl(new ICaptionControl(filter1Panel.GetGridCell(0, 0, 3, 4).GetVShifted(20.f), kParamFilter1Type, dropdownText), kCtrlTagFilter1Type, "Flt1")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        IControl* filterModeControl{ pGraphics->GetControlWithTag(kCtrlTagFilter1Mode) }; // Control containing list of filter modes for a given type
        int filterModeParam{ kParamFilter1ModeVSF + std::max(static_cast<int>(pControl->GetValue() * (kNumFilters - 1.) - 1), 0) };
        // Get range of possible modes
        filterModeControl->SetParamIdx(filterModeParam);
        filterModeControl->SetValue(std::min(filterModeControl->GetValue(), GetParam(filterModeParam)->GetMax()));
        // Update Controls
        if (pControl->GetValue() == kNoFilter)
          filterModeControl->SetDisabled(true);
        else if (pControl->GetValue() * (kNumFilters - 1) == kComb)
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
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 1, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter1Cutoff, "Cutoff"), kCtrlTagFilter1Cutoff, "Flt1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 2, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter1Resonance, "Res."), kCtrlTagFilter1Resonance, "Flt1");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 3, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter1Drive, "Drive"), kCtrlTagFilter1Drive, "Flt1");
    pGraphics->AttachControl(new IVTabSwitchControl(filter1Panel.GetGridCell(0, 4, 1, 5).GetPadded(-10.f).GetVPadded(-20.f).GetHShifted(-6.f), setFilterOscillator, { "Osc 1", "Osc 2" }, "", toggleStyle, EVShape::Rectangle, EDirection::Vertical), kCtrlTagFilter1Osc, "Flt1");
    // FF/FB/Delay Time (comb filter)
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 1, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter1FF, "F/F"), kCtrlTagFilter1FF, "Flt1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 2, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter1FB, "F/B"), kCtrlTagFilter1FB, "Flt1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter1Panel.GetGridCell(0, 3, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter1Delay, "Delay"), kCtrlTagFilter1Delay, "Flt1")->Hide(true);

    pGraphics->AttachControl(new TablitsaVGroupControl(filter1Panel.GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "Filter 1"));

    // Set comb filter controls' modulation parameters to the normal cutoff/res/drive controls' ones
    dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1FF))->GetModulationFrom(kParamFilter1Cutoff);
    dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1FB))->GetModulationFrom(kParamFilter1Resonance);
    dynamic_cast<TablitsaIVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1Delay))->GetModulationFrom(kParamFilter1Drive);

    /* Filter 2 */
    const IRECT filter2Panel = lowerPanel.GetGridCell(2, 0, 3, 3);
    pGraphics->AttachControl(new ICaptionControl(filter2Panel.GetGridCell(0, 0, 3, 4).GetVShifted(20.f), kParamFilter2Type, dropdownText), kCtrlTagFilter2Type, "Flt2")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        IControl* filterModeControl{ pGraphics->GetControlWithTag(kCtrlTagFilter2Mode) };
        int filterModeParam{ kParamFilter2ModeVSF + std::max(static_cast<int>(pControl->GetValue() * (kNumFilters - 1.) - 1), 0) };
        filterModeControl->SetParamIdx(filterModeParam);
        filterModeControl->SetValue(std::min(filterModeControl->GetValue(), GetParam(filterModeParam)->GetMax()));
        if (pControl->GetValue() == kNoFilter || pControl->GetValue() == kComb)
          filterModeControl->SetDisabled(true);
        else if (pControl->GetValue() * (kNumFilters - 1) == kComb)
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
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 1, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter2Cutoff, "Cutoff"), kNoTag, "Flt2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 2, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter2Resonance, "Res."), kNoTag, "Flt2");
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 3, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter2Drive, "Drive"), kNoTag, "Flt2");
    pGraphics->AttachControl(new IVTabSwitchControl(filter2Panel.GetGridCell(0, 4, 1, 5).GetPadded(-10.f).GetVPadded(-20.f).GetHShifted(-6.f), setFilterOscillator, { "Osc 1", "Osc 2" }, "", toggleStyle, EVShape::Rectangle, EDirection::Vertical), kCtrlTagFilter2Osc, "Flt2")->SetValue(1.);
    // FF/FB/Delay Time (comb filter)
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 1, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter2FF, "F/F"), kCtrlTagFilter2FF, "Flt1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 2, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter2FB, "F/B"), kCtrlTagFilter2FB, "Flt1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(filter2Panel.GetGridCell(0, 3, 1, 5).GetHShifted(16.f).GetCentredInside(dialSize), kParamFilter2Delay, "Delay"), kCtrlTagFilter2Delay, "Flt1")->Hide(true);

    pGraphics->AttachControl(new TablitsaVGroupControl(filter2Panel.GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "Filter 2"));

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

    /*
    LFO Panels
    */
    const IRECT lfoPanel = lowerPanel.GetGridCell(0, 2, 1, 3);
    // LFO1
    pGraphics->AttachControl(new TablitsaIVModKnobControl(lfoPanel.GetGridCell(0, 0, 3, 5).GetCentredInside(dialSize), kParamLFO1RateHz, "Rate"), kNoTag, "LFO1");
    pGraphics->AttachControl(new TablitsaIVKnobControl(lfoPanel.GetGridCell(0, 0, 3, 5).GetCentredInside(dialSize), kParamLFO1RateTempo, "Rate"), kNoTag, "LFO1")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(lfoPanel.GetGridCell(0, 2, 3, 5).GetCentredInside(dialSize), kParamLFO1Amp, "Amp"), kNoTag, "LFO1");
    pGraphics->AttachControl(new ICaptionControl(lfoPanel.GetGridCell(0, 1, 6, 5).GetVPadded(-20.f), kParamLFO1Shape, dropdownText), kNoTag, "LFO1")->DisablePrompt(false);
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(1, 1, 9, 5).GetVPadded(-6.f), kParamLFO1RateMode, "", toggleStyle, "Sync Off", "Sync On"),
      kCtrlTagLFO1RateMode, "LFO1")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
      bool sync = pControl->GetValue() > 0.5;
      pGraphics->HideControl(kParamLFO1RateHz, sync);
      pGraphics->HideControl(kParamLFO1RateTempo, !sync);
      //pGraphics->SetAllControlsDirty();
    });
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(2, 1, 9, 5).GetVPadded(-6.f).GetVShifted(-15.f), kParamLFO1Restart, "Restart", toggleStyle, "Free Run", "Retrigger"), kNoTag, "LFO1");
    pGraphics->AttachControl(new ModSliderControl(lfoPanel.GetGridCell(0, 4, 3, 5).GetVPadded(-18.f), kNoParameter, "Depth"), kCtrlTagLFO1Depth, "LFO1");
    pGraphics->AttachControl(new TablitsaVGroupControl(lfoPanel.GetGridCell(0, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "LFO 1"));

    // LFO 2
    pGraphics->AttachControl(new TablitsaIVModKnobControl(lfoPanel.GetGridCell(1, 0, 3, 5).GetCentredInside(dialSize), kParamLFO2RateHz, "Rate"), kNoTag, "LFO2");
    pGraphics->AttachControl(new TablitsaIVKnobControl(lfoPanel.GetGridCell(1, 0, 3, 5).GetCentredInside(dialSize), kParamLFO2RateTempo, "Rate"), kNoTag, "LFO2")->Hide(true);
    pGraphics->AttachControl(new TablitsaIVModKnobControl(lfoPanel.GetGridCell(1, 2, 3, 5).GetCentredInside(dialSize), kParamLFO2Amp, "Amp"), kNoTag, "LFO2");
    pGraphics->AttachControl(new ICaptionControl(lfoPanel.GetGridCell(2, 1, 6, 5).GetVPadded(-20.f), kParamLFO2Shape, dropdownText), kNoTag, "LFO2")->DisablePrompt(false);
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(4, 1, 9, 5).GetVPadded(-6.f), kParamLFO2RateMode, "", toggleStyle, "Sync Off", "Sync On"),
      kCtrlTagLFO2RateMode, "LFO2")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
        bool sync = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamLFO2RateHz, sync);
        pGraphics->HideControl(kParamLFO2RateTempo, !sync);
        //pGraphics->SetAllControlsDirty();
        });
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(5, 1, 9, 5).GetVPadded(-6.f).GetVShifted(-15.f), kParamLFO2Restart, "Restart", toggleStyle, "Free Run", "Retrigger"), kNoTag, "LFO2");
//    pGraphics->AttachControl(new IVPlotControl(lfoPanel.GetGridCell(1, 3, 3, 5).GetVPadded(-20.f), ...), kCtrlTagLFO2Vis, "LFO2");
    pGraphics->AttachControl(new ModSliderControl(lfoPanel.GetGridCell(1, 4, 3, 5).GetVPadded(-18.f), kNoParameter, "Depth"), kCtrlTagLFO2Depth, "LFO2");
    pGraphics->AttachControl(new TablitsaVGroupControl(lfoPanel.GetGridCell(1, 0, 3, 1).GetPadded(modulePad[0], modulePad[1], modulePad[2], modulePad[3]), "LFO 2"));

    /* Sequencer */
    const IRECT sequencerDisplayPanel = lowerPanel.GetGridCell(2, 2, 3, 3).GetFromLeft(PLUG_WIDTH - lfoPanel.L).GetVShifted(9.f).GetPadded(10.f, 12.f, 0.f, 0.f);
    const IRECT sequencerParamsPanel = IRECT(lowerPanel.R, lowerPanel.T, globalControls.R, lowerPanel.B - lowerPanel.H() / 3.f).GetHShifted(-2.f);

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
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(4, 1, 7, 3).GetPadded(-10.f, -2.f, -10.f, -2.f), kParamSequencerRateMode, "", toggleStyle, "Sync Off", "Sync On"),
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
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(5, 1, 7, 3).GetPadded(-10.f, -2.f, -10.f, -2.f), kParamSequencerRestart, "", toggleStyle, "Free Run", "Restart"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new ModSliderControl(sequencerDisplayPanel.GetGridCell(0, 5, 1, 6).GetVPadded(-18.f).GetVShifted(-5.f).GetMidHPadded(20.f), kNoParameter, "Depth"), kCtrlTagSequencerDepth, "Sequencer");
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(4, 2, 7, 3).GetPadded(-15.f, -2.f, -15.f, -2.f), kNoParameter, "", toggleStyle, "Cont.", "Quant."), kCtrlTagSequencerQuant, "SequencerParams")->SetActionFunction(
      [pGraphics, this](IControl* pControl) {
        mSequencerIsQuantized = pControl->GetValue();
        int nSteps = mSequencerIsQuantized > 0.5 ? 12 : 0;
        dynamic_cast<IVTrackControlBase*>(pGraphics->GetControlWithTag(kCtrlTagSequencer))->SetNSteps(nSteps);
      });
    pGraphics->AttachControl(new IVButtonControl(sequencerParamsPanel.GetGridCell(5, 2, 7, 3).GetPadded(-15.f, -2.f, -15.f, -2.f), stepRandomizeAction, "Randomize", pushButtonStyle.WithLabelText(pushButtonStyle.labelText.WithSize(12.f))), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new TablitsaVGroupControl(sequencerParamsPanel.GetPadded(modulePad[0], modulePad[1] + 5.f, modulePad[2], modulePad[3]), "Sequencer"));

    
#ifndef VST3_API
    pGraphics->AttachControl(new IVButtonControl(b.GetFromBottom(KEYBOARD_HEIGHT + 15.f).GetFromRight(170.f).GetFromTop(30.), SplashClickActionFunc,
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

    // Popup Menus
    dynamic_cast<DropdownListControl*>(pGraphics->GetControlWithTag(kCtrlTagEffect1Tab))->AttachPopupMenu();

  };
#endif
}

bool Tablitsa::SerializeState(IByteChunk& chunk) const
{
  // serialize the multislider state before serializing the regular params
  for (int i = 0; i < kNumSeqSteps; i++)
  {
    chunk.Put(&mDSP.mSeqSteps[i]);
  }
  chunk.Put(&mDSP.mLoadedWavetables[0]);
  chunk.Put(&mDSP.mLoadedWavetables[1]);
  chunk.Put(&mSequencerIsQuantized);
  chunk.Put(&mDelayIsSynced);
  
  return SerializeParams(chunk); // must remember to call SerializeParams at the end
}

// this over-ridden method is called when the host is trying to load the plug-in state and you need to unpack the data into your algorithm
int Tablitsa::UnserializeState(const IByteChunk& chunk, int startPos)
{
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
  // Sequencer quantization
  startPos = chunk.Get(&v, startPos);
  mSequencerIsQuantized = v;
  // Delay tempo sync
  startPos = chunk.Get(&v, startPos);
  mDelayIsSynced = v;

  assert(mDSP.mLoadedWavetables[0] == static_cast<int>(mDSP.mLoadedWavetables[0]));
  assert(mDSP.mLoadedWavetables[1] == static_cast<int>(mDSP.mLoadedWavetables[1]));
  
  // If UI exists
  if (GetUI())
    UpdateUIControls(); 

  // must remember to call UnserializeParams at the end
  return UnserializeParams(chunk, startPos);
}

void Tablitsa::OnUIOpen()
{
  UpdateUIControls();
}

void Tablitsa::UpdateUIControls()
{
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
  mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[0]) - 1, 0);
  mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[1]) - 1, 1);
  GetUI()->GetControlWithTag(kCtrlTagSequencerQuant)->SetValue(mSequencerIsQuantized);

  SendCurrentParamValuesFromDelegate();
  GetUI()->SetAllControlsDirty();

  // Trigger the ActionFunctions of controls that change the user interface
  /*GetUI()->GetControlWithTag(kCtrlTagLFO1RateMode)->SetDirty(true);
  GetUI()->GetControlWithTag(kCtrlTagLFO2RateMode)->SetDirty(true);
  GetUI()->GetControlWithTag(kCtrlTagSequencerRateMode)->SetDirty(true);
  GetUI()->GetControlWithTag(kCtrlTagSequencerQuant)->SetDirty(true);
  GetUI()->GetControlWithTag(kCtrlTagFilter1Type)->SetDirty(true);
  GetUI()->GetControlWithTag(kCtrlTagFilter2Type)->SetDirty(true);
  GetUI()->GetControlWithTag(kCtrlTagGlideMode)->SetDirty(true);*/
  for (auto i{ 0 }; i < std::size(kStartupTriggerControls); ++i)
  {
    IControl* pControl = GetUI()->GetControlWithTag(kStartupTriggerControls[i]);
    if (!pControl->IsHidden())
      pControl->SetDirty(true);
  }
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

  if (mDSP.mStepPos != mDSP.mPrevPos)
  {
    mDSP.mPrevPos = mDSP.mStepPos;
    SendControlMsgFromDelegate(kCtrlTagSequencer, IVMultiSliderControl<>::kMsgTagSetHighlight, sizeof(int), &mDSP.mPrevPos);
  }

  // Update Sequencer Position
  if (mDSP.GetSequencerStep() != mDSP.mSeqPos)
  {
    mDSP.mSeqPos = mDSP.GetSequencerStep();
    SendControlMsgFromDelegate(kCtrlTagSequencer, IVMultiSliderControl<>::kMsgTagSetHighlight, sizeof(int), &mDSP.mSeqPos);
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
      IByteChunk newPreset;
      SerializeState(newPreset);

      bool useSystemDialog{ (bool)LOAD_WITH_DIALOG };
      if (useSystemDialog)
      {
        WDL_String fileStr("NewPreset");
        WDL_String pathStr(GetDataPath("\\Tablitsa\\presets\\").c_str());
        GetUI()->PromptForFile(fileStr, pathStr, EFileAction::Save, "TPST");
        filePath = fileStr.Get();
        if (filePath == "\0")
          break;
      }
      else
      {
        filePath = GetDataPath("\\Tablitsa\\presets\\") + "TEST.TPST";
      }

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
      mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[0]) - 1, 0);
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
    case kMsgFilterOscChanged:
    {
      IControl* pControl = GetUI()->GetControlWithTag(ctrlTag);
      if (ctrlTag == kCtrlTagFilter1Osc)
      {
        mDSP.UpdateFilterSource(0, static_cast<int>(pControl->GetValue()));
        mDSP.UpdateFilterSource(1, 1 - static_cast<int>(pControl->GetValue()));
      }
      else
      {
        mDSP.UpdateFilterSource(1, static_cast<int>(pControl->GetValue()));
        mDSP.UpdateFilterSource(0, 1 - static_cast<int>(pControl->GetValue()));
      }
      break;
    }
    case kMsgEffect1Changed:
    case kMsgEffect2Changed:
    case kMsgEffect3Changed:
    {
      int effectIdx = *reinterpret_cast<const int*>(pData);
      switch (effectIdx)
      {
      case kDelayEffect:
        mDSP.mEffects[msgTag - kMsgEffect1Changed] = new DelayEffect<iplug::sample, TABLITSA_MAX_DELAY_SAMP>(GetSampleRate(), &mDSP.mGlobalMetronome);
        break;
      case kSampleAndHoldEffect:
        mDSP.mEffects[msgTag - kMsgEffect1Changed] = new SampleAndHold<iplug::sample>(GetSampleRate());
        break;
      default:
        mDSP.mEffects[msgTag - kMsgEffect1Changed] = new Effect<iplug::sample>(GetSampleRate());
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

/* implement this and return true to trigger your custom about box, when someone clicks about in the menu of a standalone app or VST3 plugin */
bool Tablitsa::OnHostRequestingAboutBox() { return false; } // See IPlugAPP_dialog.cpp

/* implement this and return true to trigger your custom help info, when someone clicks help in the menu of a standalone app or VST3 plugin */
bool Tablitsa::OnHostRequestingProductHelp() { return false; }

int Tablitsa::GetActiveModIdx() const
{
  return mActiveModIdx;
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

void Tablitsa::SavePreset(IByteChunk& byteData, const char* filename, bool isBackup)
{
  assert(byteData.Size() != 0);

  if (!strcmp(filename, "_Tablitsa_Default.TPST"))
  {
    GetUI()->ShowMessageBox("Cannot overwrite default preset", "Error", EMsgBoxType::kMB_OK);
  }

  std::ofstream presetFile(filename);

  if (!presetFile.is_open())
  {
    if (isBackup)
      GetUI()->ShowMessageBox("Could not open backup file. Check directory and try again.", "Error", EMsgBoxType::kMB_OK);
    else
      GetUI()->ShowMessageBox("Cannot access selected file path.", "Error", EMsgBoxType::kMB_OK);
    presetFile.close();
    return;
  }

  unsigned char* allBytes = static_cast<unsigned char*>(byteData.GetData());
  presetFile.write(reinterpret_cast<char*>(allBytes), byteData.Size());
  presetFile.close();
}

IByteChunk Tablitsa::LoadPreset(const char* filename, bool isBackup)
{
  IByteChunk presetData;

  std::vector<char> fileData = ReadAllBytes(filename);

  if (fileData.size() == 0)
  {
    if (isBackup)
      GetUI()->ShowMessageBox("Backup file not found", "Error", EMsgBoxType::kMB_OK);
    else
      GetUI()->ShowMessageBox("Preset file not found", "Error", EMsgBoxType::kMB_OK);
    return presetData;
  }

  for (auto b : fileData)
    presetData.Put(&b);

  return presetData;
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
  ifs.seekg(0, std::ios::beg);
  // TODO: Exception handling here in case the vector is empty
  ifs.read(&allBytes[0], pos);
  ifs.close();
  return allBytes;
}
