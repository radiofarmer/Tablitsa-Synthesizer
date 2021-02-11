
#include "Tablitsa.h"
#include "IPlug_include_in_plug_src.h"
#include "LFO.h"
#include "CustomControls.h"

#define UI_NROWS 4
#define UI_NCOLS 3
#define KEYBOARD
#define FILTER_PARTITION 0.8
#define ENVELOPE_PARTITION 0.8
#define LFO_PARTITION 1.4


const IVStyle modKnobStyle{ DEFAULT_STYLE.WithLabelText(17.) };

Tablitsa::Tablitsa(const InstanceInfo& info)
: Plugin(info, MakeConfig(kNumParams, kNumPresets))
{
  /*Master Controls*/
#ifdef VST3_API
  GetParam(kParamGain)->InitDouble("Gain", 25., 0., 100.0, 0.01, "%", IParam::kFlagsNone, "Master");
#else
  GetParam(kParamGain)->InitGain("Gain", -24., -72., 0.);
#endif
  GetParam(kParamPan)->InitAngleDegrees("Pan", 0., -90., 90.);
  GetParam(kParamNoteGlideTime)->InitMilliseconds("Note Glide Time", 0., 0.0, 3000.);
  GetParam(kParamNoteGlideRate)->InitDouble("Note Glide Rate", 12., 0.5, 128., 0.1, "St/s");
  GetParam(kParamPortamentoMode)->InitBool("Equal Glide Time", true);
  GetParam(kParamMonophonic)->InitBool("Monophonic", false);
  GetParam(kParamUnisonVoices)->InitInt("Unison Voices", 1, 1, 8);
  GetParam(kParamUnisonDetune)->InitDouble("Unison Detune", 0., 0., 12., 0.01, "St.");
  GetParam(kParamUnisonChord)->InitEnum("Chord", 0, { UNISON_CHORD_LIST });
  /*Modulators*/
  GetParam(kParamEnv1Attack)->InitDouble("Env1 Attack", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv1Decay)->InitDouble("Env1 Decay", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv1Sustain)->InitDouble("Env1 Sustain", 50., 0., 100., 1, "%", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv1Release)->InitDouble("Env1 Release", 10., 2., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv1Velocity)->InitDouble("Env1 Velocity Mod", 0., 0., 1., 0.01);
  GetParam(kParamEnv2Attack)->InitDouble("Env2 Attack", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv2Decay)->InitDouble("Env2 Decay", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamEnv2Sustain)->InitDouble("Env2 Sustain", 50., 0., 100., 1, "%", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv2Release)->InitDouble("Env2 Release", 10., 2., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR");
  GetParam(kParamEnv2Velocity)->InitDouble("Env2 Velocity Mod", 0., 0., 1., 0.01);
  GetParam(kParamAmpEnvAttack)->InitDouble("AmpEnv Attack", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamAmpEnvDecay)->InitDouble("AmpEnv Decay", 10., 1., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR", IParam::ShapePowCurve(3.));
  GetParam(kParamAmpEnvSustain)->InitDouble("AmpEnv Sustain", 50., 0., 100., 1, "%", IParam::kFlagsNone, "ADSR");
  GetParam(kParamAmpEnvRelease)->InitDouble("AmpEnv Release", 10., 2., 3000., 0.1, "ms", IParam::kFlagsNone, "ADSR");
  GetParam(kParamAmpEnvVelocity)->InitDouble("AmpEnv Velocity Mod", 1., 0., 1., 0.01);
  GetParam(kParamLFO1Shape)->InitEnum("LFO1 Shape", LFO<>::kTriangle, {LFO_SHAPE_VALIST});
  GetParam(kParamLFO1RateHz)->InitFrequency("LFO1 Rate", 1., 0.01, 40.);
  GetParam(kParamLFO1RateTempo)->InitEnum("LFO1 Rate", LFO<>::k1, {LFO_TEMPODIV_VALIST});
  GetParam(kParamLFO1RateMode)->InitBool("LFO1 Sync", false);
  GetParam(kParamLFO1Restart)->InitBool("LFO1 Retrigger", false);
  GetParam(kParamLFO1Amp)->InitDouble("LFO Amp", 1., -1., 1., 0.01, "", IParam::kFlagsNone, "LFO");
  GetParam(kParamLFO2Shape)->InitEnum("LFO2 Shape", LFO<>::kTriangle, { LFO_SHAPE_VALIST });
  GetParam(kParamLFO2RateHz)->InitFrequency("LFO2 Rate", 1., 0.01, 40.);
  GetParam(kParamLFO2RateTempo)->InitEnum("LFO2 Rate", LFO<>::k1, { LFO_TEMPODIV_VALIST });
  GetParam(kParamLFO2RateMode)->InitBool("LFO2 Sync", false);
  GetParam(kParamLFO2Restart)->InitBool("LFO2 Retrigger", false);
  GetParam(kParamLFO2Amp)->InitDouble("LFO Amp", 1., -1., 1., 0.01, "", IParam::kFlagsNone, "LFO");
  GetParam(kParamSequencerRateHz)->InitFrequency("Sequencer Rate", 1., 0.01, 40.);
  GetParam(kParamSequencerRateTempo)->InitEnum("Sequencer Rate", LFO<>::k4th, { LFO_TEMPODIV_VALIST });
  GetParam(kParamSequencerRateMode)->InitBool("Sequencer Sync", true);
  GetParam(kParamSequencerRestart)->InitBool("Sequencer Retrigger", false);
  GetParam(kParamSequencerSteps)->InitInt("Sequencer Length", 16, 1, 16);
  GetParam(kParamSequencerGlide)->InitPercentage("Sequencer Glide");
//  GetParam(kParamSequencerCurve)->InitEnum()
  GetParam(kParamSequencerAmp)->InitDouble("Seqeuencer Amp", 1., 0., 1., 0.01, "", IParam::kFlagsNone, "Sequencer");

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
  GetParam(kParamWavetable1PitchVel)->InitDouble("WT1 Pitch Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PitchKTk)->InitDouble("WT1 Pitch LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PitchRnd)->InitDouble("WT1 Pitch Seq Mod", 0., -1., 1., 0.01);
  // Wavetable 1 Position
  GetParam(kParamWavetable1Pos)->InitDouble("WT1 Pos", 0., 0., 1., 0.01);
  GetParam(kParamWavetable1PosEnv1)->InitDouble("WT1 Pos Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosEnv2)->InitDouble("WT1 Pos Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosAmpEnv)->InitDouble("WT1 Pos AmpEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosLFO1)->InitDouble("WT1 Pos LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosLFO2)->InitDouble("WT1 Pos L Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosSeq)->InitDouble("WT1 Pos Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosVel)->InitDouble("WT1 Pos Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosKTk)->InitDouble("WT1 Pos KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1PosRnd)->InitDouble("WT1 Pos Rnd Mod", 0., -1., 1., 0.01);
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
  // Wavetable 1 Sub
  GetParam(kParamWavetable1Sub)->InitDouble("WT1 Sub", 0., 0., 1., 0.01);
  GetParam(kParamWavetable1SubEnv1)->InitDouble("WT1 Sub Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1SubEnv2)->InitDouble("WT1 Sub Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1SubAmpEnv)->InitDouble("WT1 Sub SubEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1SubLFO1)->InitDouble("WT1 Sub LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1SubLFO2)->InitDouble("WT1 Sub LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1SubSeq)->InitDouble("WT1 Sub Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1SubVel)->InitDouble("WT1 Sub Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1SubKTk)->InitDouble("WT1 Sub KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable1SubRnd)->InitDouble("WT1 Sub Rnd Mod", 0., -1., 1., 0.01);
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
  // Wavetable 2 Sub
  GetParam(kParamWavetable2Sub)->InitDouble("WT2 Sub", 0., 0., 1., 0.01);
  GetParam(kParamWavetable2SubEnv1)->InitDouble("WT2 Sub Env1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2SubEnv2)->InitDouble("WT2 Sub Env2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2SubAmpEnv)->InitDouble("WT2 Sub SubEnv Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2SubLFO1)->InitDouble("WT2 Sub LFO1 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2SubLFO2)->InitDouble("WT2 Sub LFO2 Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2SubSeq)->InitDouble("WT2 Sub Seq Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2SubVel)->InitDouble("WT2 Sub Vel Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2SubKTk)->InitDouble("WT2 Sub KTk Mod", 0., -1., 1., 0.01);
  GetParam(kParamWavetable2SubRnd)->InitDouble("WT2 Sub Rnd Mod", 0., -1., 1., 0.01);
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
  GetParam(kParamPhaseModAmount)->InitPercentage("Phase Modulator Amount", 25., 0., 100., 1.);
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
  GetParam(kParamRingModAmount)->InitPercentage("Ring Modulator Amount", 25., 0., 100., 1.);
  GetParam(kParamRingModAmount)->InitPercentage("Ring Modulator Amount", 25., 0., 100., 1.);
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
  GetParam(kParamDelayTimeLMilliseconds)->InitDouble("Stereo Delay Time L (ms)", 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Delay", IParam::ShapePowCurve(3.));
  GetParam(kParamDelayTimeRMilliseconds)->InitDouble("Stereo Delay Time R (ms)", 100., 1., 5000., 1., "ms", IParam::kFlagsNone, "Delay", IParam::ShapePowCurve(3.));
  GetParam(kParamDelayTimeLBeats)->InitEnum("Stereo Delay Time L (tempo)", LFO<>::k8th, { DELAY_TEMPODIV_VALIST });
  GetParam(kParamDelayTimeRBeats)->InitEnum("Stereo Delay Time R (tempo)", LFO<>::k8th, { DELAY_TEMPODIV_VALIST });
  GetParam(kParamDelayTimeMode)->InitBool("Delay Tempo Sync", false);
  GetParam(kParamDelayFeedback)->InitPercentage("Delay Feedback");
  GetParam(kParamDelayMix)->InitPercentage("Delay Mix", 0.);

#ifndef VST3_API
#define KEYBOARD_HEIGHT 150.f
#else
#define KEYBOARD_HEIGHT 0.f
#endif
#if IPLUG_EDITOR // http://bit.ly/2S64BDd
  mMakeGraphicsFunc = [&]() {
    return MakeGraphics(*this, PLUG_WIDTH, PLUG_HEIGHT + KEYBOARD_HEIGHT, PLUG_FPS, GetScaleForScreen(PLUG_WIDTH, PLUG_HEIGHT));
  };
  
  mLayoutFunc = [&](IGraphics* pGraphics) {
    pGraphics->AttachCornerResizer(EUIResizerMode::Scale, false);
    pGraphics->AttachPanelBackground(COLOR_GRAY);
    pGraphics->EnableMouseOver(true);
    pGraphics->EnableMultiTouch(true);

#ifdef OS_WEB
    pGraphics->AttachPopupMenuControl();
#endif

    //    pGraphics->EnableLiveEdit(true);
    pGraphics->LoadFont("Roboto-Regular", ROBOTO_FN);
    const IRECT b = pGraphics->GetBounds().GetPadded(-10.f);
    constexpr float dialSize{ 80.f };
    IRECT controls = b.GetOffset(0, 0, 0, -KEYBOARD_HEIGHT).GetHPadded(-5.f);
#ifndef VST3_API
    IRECT keyboardBounds = b.GetFromBottom((float)KEYBOARD_HEIGHT).GetPadded(-10.f);
    IRECT wheelsBounds = keyboardBounds.ReduceFromLeft(100.f).GetPadded(-10.f);

    // Keyboard
    pGraphics->AttachControl(new IVKeyboardControl(keyboardBounds, 36), kCtrlTagKeyboard);
    // Mod wheels
    pGraphics->AttachControl(new IWheelControl(wheelsBounds.FracRectHorizontal(0.5)), kCtrlTagBender);
    pGraphics->AttachControl(new IWheelControl(wheelsBounds.FracRectHorizontal(0.5, true), IMidiMsg::EControlChangeMsg::kModWheel));
    //    pGraphics->AttachControl(new IVMultiSliderControl<4>(b.GetGridCell(0, 2, 2).GetPadded(-30), "", DEFAULT_STYLE, kParamAttack, EDirection::Vertical, 0.f, 1.f));

#endif
    const ISVG tableImg = pGraphics->LoadSVG(PERIODIC_TABLE_FN);
    pGraphics->AttachControl(new PeriodicTable(controls.GetReducedFromBottom(controls.H() * 0.5f).GetFromLeft(720.f).GetPadded(-30.f).GetVShifted(20.f), tableImg, { kParamWavetable1, kParamWavetable2 }),
      kCtrlTagPeriodicTable); // TODO: Set action function to store non-modulatable parameter value
    dynamic_cast<PeriodicTable*>(pGraphics->GetControlWithTag(kCtrlTagPeriodicTable))->LoadValues();

    // Wavetable controls
    const IRECT osc1Panel = controls.GetReducedFromLeft(720.f).GetReducedFromBottom(controls.H() / 2.f).GetGridCell(0, 0, 2, 3).GetPadded(-10.f).GetHPadded(-12.f);
    pGraphics->AttachControl(new IVModKnobControl(osc1Panel.GetGridCell(0, 0, 2, 2).GetCentredInside(dialSize), kParamWavetable1Pitch, "Pitch", DEFAULT_STYLE, true, 0.05), kNoTag, "Osc1");
    pGraphics->AttachControl(new IVModKnobControl(osc1Panel.GetGridCell(0, 1, 2, 2).GetCentredInside(dialSize), kParamWavetable1Amp, "Amp"), kNoTag, "Osc1");
    pGraphics->AttachControl(new IVModKnobControl(osc1Panel.GetGridCell(1, 0, 2, 3).GetCentredInside(dialSize), kParamWavetable1Pos, "Valency"), kNoTag, "Osc1");
    pGraphics->AttachControl(new IVModKnobControl(osc1Panel.GetGridCell(1, 1, 2, 3).GetCentredInside(dialSize), kParamWavetable1Bend, "Energy"), kNoTag, "Osc1");
    pGraphics->AttachControl(new IVModKnobControl(osc1Panel.GetGridCell(1, 2, 2, 3).GetCentredInside(dialSize), kParamWavetable1Sub, "Mass"), kNoTag, "Osc1");
    pGraphics->AttachControl(new IVGroupControl("Wavetable 1", "Osc1", 0.f, 15.f, 0.f, 2.f));

    const IRECT osc2Panel = controls.GetReducedFromLeft(720.f).GetReducedFromBottom(controls.H() / 2.f).GetGridCell(1, 0, 2, 3).GetPadded(-10.f).GetHPadded(-12.f);
    pGraphics->AttachControl(new IVModKnobControl(osc2Panel.GetGridCell(0, 0, 2, 2).GetCentredInside(dialSize), kParamWavetable2Pitch, "Pitch", DEFAULT_STYLE, true, 0.05), kNoTag, "Osc2");
    pGraphics->AttachControl(new IVModKnobControl(osc2Panel.GetGridCell(0, 1, 2, 2).GetCentredInside(dialSize), kParamWavetable2Amp, "Amp"), kNoTag, "Osc2");
    pGraphics->AttachControl(new IVModKnobControl(osc2Panel.GetGridCell(1, 0, 2, 3).GetCentredInside(dialSize), kParamWavetable2Pos, "Valency"), kNoTag, "Osc2");
    pGraphics->AttachControl(new IVModKnobControl(osc2Panel.GetGridCell(1, 1, 2, 3).GetCentredInside(dialSize), kParamWavetable2Bend, "Energy"), kNoTag, "Osc2");
    pGraphics->AttachControl(new IVModKnobControl(osc2Panel.GetGridCell(1, 2, 2, 3).GetCentredInside(dialSize), kParamWavetable2Sub, "Mass"), kNoTag, "Osc2");
    pGraphics->AttachControl(new IVGroupControl("Wavetable 2", "Osc2", 0.f, 15.f, 0.f, 2.f));

    /* Master Controls */
    const IRECT masterPanel = controls.GetReducedFromLeft(720.f + osc1Panel.W() + 30.f).GetReducedFromBottom(controls.H() / 2.f).GetPadded(-10.f);
    // Gain, Pan, and Volume Meter
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(0, 0, 4, 6).GetCentredInside(dialSize), kParamGain, "Gain"), kNoTag, "Master");
    pGraphics->AttachControl(new IVModKnobControl(masterPanel.GetGridCell(0, 1, 4, 6).GetCentredInside(dialSize), kParamPan, "Pan"), kNoTag, "Master");
    // Static modulators
    pGraphics->AttachControl(new ModSliderControl(masterPanel.GetGridCell(1, 0, 2, 1).GetVPadded(-16.f).GetGridCell(0, 0, 3, 3).GetHPadded(-24.f).GetHShifted(-15.f), kNoParameter, "Velocity", DEFAULT_STYLE, false, EDirection::Horizontal), kCtrlTagVelDepth, "StaticMods");
    pGraphics->AttachControl(new ModSliderControl(masterPanel.GetGridCell(1, 0, 2, 1).GetVPadded(-16.f).GetGridCell(1, 0, 3, 3).GetHPadded(-24.f).GetHShifted(-15.f), kNoParameter, "Keytrack", DEFAULT_STYLE, false, EDirection::Horizontal), kCtrlTagKTkDepth, "StaticMods");
    pGraphics->AttachControl(new ModSliderControl(masterPanel.GetGridCell(1, 0, 2, 1).GetVPadded(-16.f).GetGridCell(2, 0, 3, 3).GetHPadded(-24.f).GetHShifted(-15.f), kNoParameter, "Trg. Rnd.", DEFAULT_STYLE, false, EDirection::Horizontal), kCtrlTagRndDepth, "StaticMods");
    pGraphics->AttachControl(new IVGroupControl("Midi Mods", "StaticMods", 5.f, 15.f, 5.f, 12.f));
    // Polyphony and Portamento
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(3, 1, 12, 6).GetVPadded(-2.f), kParamPortamentoMode, "", DEFAULT_STYLE.WithShowLabel(false), "Rate", "Time"),
      kNoTag, "Master")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
        bool constantTime = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamNoteGlideTime, !constantTime);
        pGraphics->HideControl(kParamNoteGlideRate, constantTime);
        pGraphics->SetAllControlsDirty();
        });
    // Mono/Poly and Portamento Rate/Time switches
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(1, 0, 4, 6).GetCentredInside(dialSize), kParamNoteGlideRate, "Glide Rate"), kNoTag, "Master")->Hide(true);
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(1, 0, 4, 6).GetCentredInside(dialSize), kParamNoteGlideTime, "Glide Time"), kNoTag, "Master")->DisablePrompt(false);
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(4, 1, 12, 6).GetVPadded(-2.f), kParamMonophonic, "", DEFAULT_STYLE.WithShowLabel(false), "Poly", "Mono"));
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(5, 1, 12, 6).GetVPadded(-2.f), kParamLegato, "", DEFAULT_STYLE.WithShowLabel(false), "Retrigger", "Legato"));
    pGraphics->AttachControl(new IVLEDMeterControl<2>(masterPanel.GetGridCell(0, 5, 2, 6).GetVPadded(-20.f)), kCtrlTagMeter, "Master");
    // Unison and Stereo Width
    pGraphics->AttachControl(new IVNumberBoxControl(masterPanel.GetGridCell(0, 2, 4, 6).GetHPadded(-10.f).GetVPadded(-18.f), kParamUnisonVoices, nullptr, "Unison", DEFAULT_STYLE.WithValueText(IText(18.f))), kNoTag, "Stereo");
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(0, 3, 4, 6).GetCentredInside(dialSize), kParamUnisonDetune, "Detune"), kNoTag, "Stereo");
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(0, 4, 4, 6).GetCentredInside(dialSize), kParamUnisonChord, "Chord"), kNoTag, "Stereo");

    /* Effects Panel */
    auto channelLockFunction = IActionFunction([pGraphics, this](IControl* pControl) {
      int ctrlTag = pControl->GetTag() == kCtrlTagDelayLMilliseconds || pControl->GetTag() == kCtrlTagDelayLBeats ? pControl->GetTag() + 1 : pControl->GetTag() - 1;
      pGraphics->GetControlWithTag(ctrlTag)->SetValue(pControl->GetValue());
      GetUI()->GetDelegate()->SendParameterValueFromUI(pGraphics->GetControlWithTag(ctrlTag)->GetParamIdx(), pControl->GetValue());
      pGraphics->GetControlWithTag(ctrlTag)->SetDirty(false);
      });
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(2, 2, 4, 7).GetVShifted(5.f).GetVPadded(-12.f).GetCentredInside(dialSize), kParamDelayTimeLMilliseconds, "Left"), kCtrlTagDelayLMilliseconds, "Delay");
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(2, 3, 4, 7).GetVShifted(5.f).GetVPadded(-12.f).GetCentredInside(dialSize), kParamDelayTimeRMilliseconds, "Right"), kCtrlTagDelayRMilliseconds, "Delay");
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(2, 2, 4, 7).GetVShifted(5.f).GetVPadded(-12.f).GetCentredInside(dialSize), kParamDelayTimeLBeats, "Left"), kCtrlTagDelayLBeats, "Delay")->Hide(true);
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(2, 3, 4, 7).GetVShifted(5.f).GetVPadded(-12.f).GetCentredInside(dialSize), kParamDelayTimeRBeats, "Right"), kCtrlTagDelayRBeats, "Delay")->Hide(true);
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(2, 4, 4, 7).GetVShifted(5.f).GetVPadded(-12.f).GetCentredInside(dialSize), kParamDelayFeedback, "Feedback"), kNoTag, "Delay");
    pGraphics->AttachControl(new IVKnobControl(masterPanel.GetGridCell(2, 5, 4, 7).GetVShifted(5.f).GetVPadded(-12.f).GetCentredInside(dialSize), kParamDelayMix, "Mix"), kNoTag, "Delay");
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(2, 6, 4, 7).GetVPadded(-8.f).SubRect(EDirection::Vertical, 2, 0).GetVPadded(-2.f), kParamDelayTimeMode, "", DEFAULT_STYLE.WithShowLabel(false), "Time", "Tempo"), kNoTag, "Delay")->SetActionFunction(
      [pGraphics](IControl* pControl) {
        bool isTempoSync = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamDelayTimeLMilliseconds, isTempoSync);
        pGraphics->HideControl(kParamDelayTimeRMilliseconds, isTempoSync);
        pGraphics->HideControl(kParamDelayTimeLBeats, !isTempoSync);
        pGraphics->HideControl(kParamDelayTimeRBeats, !isTempoSync);
        pGraphics->SetAllControlsDirty();
      }
    );
    pGraphics->AttachControl(new IVToggleControl(masterPanel.GetGridCell(2, 6, 4, 7) .GetVPadded(-8.f).SubRect(EDirection::Vertical, 2, 1).GetVPadded(-2.f), kNoParameter, "", DEFAULT_STYLE, "Stereo", "Mono"), kNoTag, "Delay")->SetActionFunction(
      [pGraphics, channelLockFunction](IControl* pControl) {
        bool mono = pControl->GetValue() > 0.5;
        if (mono)
        {
          pGraphics->GetControlWithTag(kCtrlTagDelayLMilliseconds)->SetActionFunction(channelLockFunction);
          pGraphics->GetControlWithTag(kCtrlTagDelayRMilliseconds)->SetActionFunction(channelLockFunction);
          pGraphics->GetControlWithTag(kCtrlTagDelayLBeats)->SetActionFunction(channelLockFunction);
          pGraphics->GetControlWithTag(kCtrlTagDelayRBeats)->SetActionFunction(channelLockFunction);
        }
        else
        {
          pGraphics->GetControlWithTag(kCtrlTagDelayLMilliseconds)->SetActionFunction([pGraphics](IControl* pControl) {});
          pGraphics->GetControlWithTag(kCtrlTagDelayRMilliseconds)->SetActionFunction([pGraphics](IControl* pControl) {});
          pGraphics->GetControlWithTag(kCtrlTagDelayLBeats)->SetActionFunction([pGraphics](IControl* pControl) {});
          pGraphics->GetControlWithTag(kCtrlTagDelayRBeats)->SetActionFunction([pGraphics](IControl* pControl) {});
        }
        pGraphics->SetAllControlsDirty();
      }
    );;
    pGraphics->AttachControl(new IVGroupControl("Delay", "Delay", -10.f, 15.f, 0.f, 0.f));

    /* Waveform Modulation Panel */
    const IRECT lowerPanel = controls.GetReducedFromTop(controls.H() * 0.5f).GetPadded(-10.f).GetReducedFromRight(250.f);
    const IRECT oscModPanel = lowerPanel.GetGridCell(0, 0, 3, 3).GetPadded(-2.f);
    pGraphics->AttachControl(new ICaptionControl(oscModPanel.GetVPadded(-10.f).GetGridCell(0, 0, 3, 3).GetVShifted(10.f), kParamOscModulator), kNoTag, "OscMod")->SetActionFunction(
      [pGraphics](IControl* pControl) {
        bool phaseMod = pControl->GetValue() < 0.5;
        pGraphics->GetControlWithTag(kCtrlTagOscModFreq)->SetParamIdx(phaseMod ? kParamPhaseModFreq : kParamRingModFreq);
        pGraphics->GetControlWithTag(kCtrlTagOscModAmt)->SetParamIdx(phaseMod ? kParamPhaseModAmount : kParamRingModAmount);
        IControl* pOsc1Switch = pGraphics->GetControlWithTag(kCtrlTagOsc1ModSwitch);
        pOsc1Switch->SetParamIdx(phaseMod ? kParamOsc1PM : kParamOsc1RM);
        IControl* pOsc2Switch = pGraphics->GetControlWithTag(kCtrlTagOsc2ModSwitch);
        pOsc2Switch->SetParamIdx(phaseMod ? kParamOsc2PM : kParamOsc2RM);
        // Update controls
        pControl->GetDelegate()->SendCurrentParamValuesFromDelegate();
        pOsc2Switch->SetDirty(true);
        pOsc1Switch->SetDirty(true);
      }
      );
    pGraphics->AttachControl(new IVToggleControl(oscModPanel.GetVPadded(-10.f).GetVPadded(-10.f).GetGridCell(1, 0, 3, 3).GetMidHPadded(35.f), kParamOsc1PM, "Osc1", DEFAULT_STYLE.WithShowLabel(false), "Osc 1 OFF", "Osc 1 ON"), kCtrlTagOsc1ModSwitch, "OscMod");
    pGraphics->AttachControl(new IVToggleControl(oscModPanel.GetVPadded(-10.f).GetVPadded(-10.f).GetGridCell(2, 0, 3, 3).GetMidHPadded(35.f), kParamOsc2PM, "Osc1", DEFAULT_STYLE.WithShowLabel(false), "Osc 2 OFF", "Osc 2 ON"), kCtrlTagOsc2ModSwitch, "OscMod");
    pGraphics->AttachControl(new IVModKnobControl(oscModPanel.GetVPadded(-10.f).GetGridCell(0, 1, 1, 3).GetCentredInside(dialSize).GetCentredInside(dialSize), kParamPhaseModFreq, "Pitch"), kCtrlTagOscModFreq, "OscMod");
    pGraphics->AttachControl(new IVModKnobControl(oscModPanel.GetGridCell(0, 2, 1, 3).GetCentredInside(dialSize).GetCentredInside(dialSize), kParamPhaseModAmount, "Amount"), kCtrlTagOscModAmt, "OscMod");
    pGraphics->AttachControl(new IVGroupControl("Oscillator Modulation", "OscMod", 0.f, 10.f, 0.f, 0.f));
    
    /*
    Filter Panel
    */
    const IRECT filter1Panel = lowerPanel.GetGridCell(1, 0, 3, 3).GetPadded(-2.f);
    /* Filter 1 */
    pGraphics->AttachControl(new ICaptionControl(filter1Panel.GetGridCell(0, 0, 3, 4).GetVShifted(20.f), kParamFilter1Type), kCtrlTagFilter1Type, "Flt1")->SetActionFunction(
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
          pGraphics->HideControl(kParamFilter1Drive, false);
          pGraphics->HideControl(kParamFilter1Delay, true);
        }
        filterModeControl->SetDirty(true);
      }
    );
    pGraphics->AttachControl(new ICaptionControl(filter1Panel.GetGridCell(1, 0, 3, 4).GetFromTop(20.f).GetVShifted(20.f), kParamFilter1ModeVSF), kCtrlTagFilter1Mode, "Flt1")->
      SetDisabled(static_cast<int>(pGraphics->GetControlWithTag(kCtrlTagFilter1Type)->GetValue()) == kNoFilter);
    // Cutoff/Resonance/Drive
    pGraphics->AttachControl(new IVModKnobControl(filter1Panel.GetGridCell(0, 1, 1, 4).GetCentredInside(dialSize), kParamFilter1Cutoff, "Cutoff"), kCtrlTagFilter1Cutoff, "Flt1");
    pGraphics->AttachControl(new IVModKnobControl(filter1Panel.GetGridCell(0, 2, 1, 4).GetCentredInside(dialSize), kParamFilter1Resonance, "Res."), kCtrlTagFilter1Resonance, "Flt1");
    pGraphics->AttachControl(new IVModKnobControl(filter1Panel.GetGridCell(0, 3, 1, 4).GetCentredInside(dialSize), kParamFilter1Drive, "Drive"), kCtrlTagFilter1Drive, "Flt1");
    // FF/FB/Delay Time (comb filter)
    pGraphics->AttachControl(new IVModKnobControl(filter1Panel.GetGridCell(0, 1, 1, 4).GetCentredInside(dialSize), kParamFilter1FF, "F/F"), kCtrlTagFilter1FF, "Flt1")->Hide(true);
    pGraphics->AttachControl(new IVModKnobControl(filter1Panel.GetGridCell(0, 2, 1, 4).GetCentredInside(dialSize), kParamFilter1FB, "F/B"), kCtrlTagFilter1FB, "Flt1")->Hide(true);
    pGraphics->AttachControl(new IVModKnobControl(filter1Panel.GetGridCell(0, 3, 1, 4).GetCentredInside(dialSize), kParamFilter1Delay, "Delay"), kCtrlTagFilter1Delay, "Flt1")->Hide(true);

    pGraphics->AttachControl(new IVGroupControl("Filter 1", "Flt1", 0.f, 15.f, 0.f, 0.f));

    // Set comb filter controls' modulation parameters to the normal cutoff/res/drive controls' ones
    dynamic_cast<IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1FF))->GetModulationFrom(kParamFilter1Cutoff);
    dynamic_cast<IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1FB))->GetModulationFrom(kParamFilter1Resonance);
    dynamic_cast<IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter1Delay))->GetModulationFrom(kParamFilter1Drive);

    /* Filter 2 */
    const IRECT filter2Panel = lowerPanel.GetGridCell(2, 0, 3, 3).GetPadded(-2.f);
    pGraphics->AttachControl(new ICaptionControl(filter2Panel.GetGridCell(0, 0, 3, 4).GetVShifted(20.f), kParamFilter2Type), kCtrlTagFilter2Type, "Flt2")->SetActionFunction(
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
    pGraphics->AttachControl(new ICaptionControl(filter2Panel.GetGridCell(1, 0, 3, 4).GetFromTop(20.f).GetVShifted(20.f), kParamFilter2ModeVSF), kCtrlTagFilter2Mode, "Flt2")->
      SetDisabled(static_cast<int>(pGraphics->GetControlWithTag(kCtrlTagFilter2Type)->GetValue()) == kNoFilter);
    pGraphics->AttachControl(new IVModKnobControl(filter2Panel.GetGridCell(0, 1, 1, 4).GetCentredInside(dialSize), kParamFilter2Cutoff, "Cutoff"), kNoTag, "Flt2");
    pGraphics->AttachControl(new IVModKnobControl(filter2Panel.GetGridCell(0, 2, 1, 4).GetCentredInside(dialSize), kParamFilter2Resonance, "Res."), kNoTag, "Flt2");
    pGraphics->AttachControl(new IVModKnobControl(filter2Panel.GetGridCell(0, 3, 1, 4).GetCentredInside(dialSize), kParamFilter2Drive, "Drive"), kNoTag, "Flt2");
    // FF/FB/Delay Time (comb filter)
    pGraphics->AttachControl(new IVModKnobControl(filter2Panel.GetGridCell(0, 1, 1, 4).GetCentredInside(dialSize), kParamFilter2FF, "F/F"), kCtrlTagFilter2FF, "Flt1")->Hide(true);
    pGraphics->AttachControl(new IVModKnobControl(filter2Panel.GetGridCell(0, 2, 1, 4).GetCentredInside(dialSize), kParamFilter2FB, "F/B"), kCtrlTagFilter2FB, "Flt1")->Hide(true);
    pGraphics->AttachControl(new IVModKnobControl(filter2Panel.GetGridCell(0, 3, 1, 4).GetCentredInside(dialSize), kParamFilter2Delay, "Delay"), kCtrlTagFilter2Delay, "Flt1")->Hide(true);

    pGraphics->AttachControl(new IVGroupControl("Filter 2", "Flt2", 0.f, 15.f, 0.f, 0.f));

    // Set comb filter controls' modulation parameters to the normal cutoff/res/drive controls' ones
    dynamic_cast<IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter2FF))->GetModulationFrom(kParamFilter2Cutoff);
    dynamic_cast<IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter2FB))->GetModulationFrom(kParamFilter2Resonance);
    dynamic_cast<IVModKnobControl*>(pGraphics->GetControlWithTag(kCtrlTagFilter2Delay))->GetModulationFrom(kParamFilter2Drive);

#pragma mark Modulation controls -

    // Envelope Sliders
    const IRECT envPanel = lowerPanel.GetGridCell(0, 1, 1, 3).GetPadded(-2.f);
    // Envelope 1
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(0, 0, 3, 6).GetCentredInside(dialSize), kParamEnv1Attack, "Att", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(0, 1, 3, 6).GetCentredInside(dialSize), kParamEnv1Decay, "Dec", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(0, 2, 3, 6).GetCentredInside(dialSize), kParamEnv1Sustain, "Sus", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(0, 3, 3, 6).GetCentredInside(dialSize), kParamEnv1Release, "Rel", modKnobStyle), kNoTag, "Env1");
    pGraphics->AttachControl(new IVSliderControl(envPanel.GetGridCell(0, 4, 3, 6).GetVPadded(-15.f), kParamEnv1Velocity, "Velocity", modKnobStyle, false, EDirection::Vertical, 4.f, 8.f, 2.f, true), kNoTag, "Env1");
    pGraphics->AttachControl(new ModSliderControl(envPanel.GetGridCell(0, 5, 3, 6).GetVPadded(-15.f), kNoParameter, "Depth", modKnobStyle), kCtrlTagEnv1Depth, "Env1");
    pGraphics->AttachControl(new IVGroupControl("Env 1", "Env1", -2.f, 10.f, 0.f, 0.f));
    // Envelope 2
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(1, 0, 3, 6).GetCentredInside(dialSize), kParamEnv2Attack, "Att", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(1, 1, 3, 6).GetCentredInside(dialSize), kParamEnv2Decay, "Dec", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(1, 2, 3, 6).GetCentredInside(dialSize), kParamEnv2Sustain, "Sus", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(1, 3, 3, 6).GetCentredInside(dialSize), kParamEnv2Release, "Rel", modKnobStyle), kNoTag, "Env2");
    pGraphics->AttachControl(new IVSliderControl(envPanel.GetGridCell(1, 4, 3, 6).GetVPadded(-15.f), kParamEnv2Velocity, "Velocity", modKnobStyle, false, EDirection::Vertical, 4.f, 8.f, 2.f, true), kNoTag, "Env2");
    pGraphics->AttachControl(new ModSliderControl(envPanel.GetGridCell(1, 5, 3, 6).GetVPadded(-15.f), kNoParameter, "Depth", modKnobStyle), kCtrlTagEnv2Depth, "Env2");
    pGraphics->AttachControl(new IVGroupControl("Env 2", "Env2", -2.f, 10.f, 0.f, 0.f));
    // Amplitude Envelope
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(2, 0, 3, 6).GetCentredInside(dialSize), kParamAmpEnvAttack, "Att", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(2, 1, 3, 6).GetCentredInside(dialSize), kParamAmpEnvDecay, "Dec", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(2, 2, 3, 6).GetCentredInside(dialSize), kParamAmpEnvSustain, "Sus", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new IVKnobControl(envPanel.GetGridCell(2, 3, 3, 6).GetCentredInside(dialSize), kParamAmpEnvRelease, "Rel", modKnobStyle), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new IVSliderControl(envPanel.GetGridCell(2, 4, 3, 6).GetVPadded(-15.f), kParamAmpEnvVelocity, "Velocity", modKnobStyle, false, EDirection::Vertical, 4.f, 8.f, 2.f, true), kNoTag, "AmpEnv");
    pGraphics->AttachControl(new ModSliderControl(envPanel.GetGridCell(2, 5, 3, 6).GetVPadded(-15.f), kNoParameter, "Depth", modKnobStyle), kCtrlTagAmpEnvDepth, "AmpEnv");
    pGraphics->AttachControl(new IVGroupControl("Amp Env", "AmpEnv", -2.f, 10.f, 0.f, 0.f));

    /*
    LFO Panels
    */
    const IRECT lfoPanel = lowerPanel.GetGridCell(0, 2, 1, 3).GetPadded(-4.f);
    // LFO1
    pGraphics->AttachControl(new IVKnobControl(lfoPanel.GetGridCell(0, 0, 3, 5).GetCentredInside(dialSize), kParamLFO1RateHz, "Rate"), kNoTag, "LFO1");
    pGraphics->AttachControl(new IVKnobControl(lfoPanel.GetGridCell(0, 0, 3, 5).GetCentredInside(dialSize), kParamLFO1RateTempo, "Rate"), kNoTag, "LFO1")->Hide(true);
    pGraphics->AttachControl(new IVKnobControl(lfoPanel.GetGridCell(0, 2, 3, 5).GetCentredInside(dialSize), kParamLFO1Amp, "Amp"), kNoTag, "LFO1");
    pGraphics->AttachControl(new ICaptionControl(lfoPanel.GetGridCell(0, 1, 6, 5).GetVPadded(-20.f), kParamLFO1Shape), kNoTag, "LFO1")->DisablePrompt(false);
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(1, 1, 9, 5).GetVPadded(-6.f), kParamLFO1RateMode, "", DEFAULT_STYLE.WithShowLabel(false), "Sync Off", "Sync On"),
      kNoTag, "LFO1")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
      bool sync = pControl->GetValue() > 0.5;
      pGraphics->HideControl(kParamLFO1RateHz, sync);
      pGraphics->HideControl(kParamLFO1RateTempo, !sync);
      pGraphics->SetAllControlsDirty();
    });
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(2, 1, 9, 5).GetVPadded(-6.f).GetVShifted(-15.f), kParamLFO1Restart, "Restart", DEFAULT_STYLE.WithShowLabel(false), "Free Run", "Retrigger"), kNoTag, "LFO1");
    pGraphics->AttachControl(new IVDisplayControl(lfoPanel.GetGridCell(0, 3, 3, 5).GetVPadded(-20.f)));
    pGraphics->AttachControl(new ModSliderControl(lfoPanel.GetGridCell(0, 4, 3, 5).GetVPadded(-15.f), kNoParameter, "Depth"), kCtrlTagLFO1Depth, "LFO1");
    //pGraphics->AttachControl(new IVDisplayControl(lfoPanel.GetGridCell(1, 1, 2, 3).Union(lfoPanel.GetGridCell(1, 2, 2, 3)), "", DEFAULT_STYLE, EDirection::Horizontal, 0.f, 1.f, 0.f, 1024), kCtrlTagLFOVis, "LFO");
    pGraphics->AttachControl(new IVGroupControl("LFO 1", "LFO1", 0.f, 10.f, 0.f, 0.f));

    // LFO 2
    pGraphics->AttachControl(new IVKnobControl(lfoPanel.GetGridCell(1, 0, 3, 5).GetCentredInside(dialSize), kParamLFO2RateHz, "Rate"), kNoTag, "LFO2");
    pGraphics->AttachControl(new IVKnobControl(lfoPanel.GetGridCell(1, 0, 3, 5).GetCentredInside(dialSize), kParamLFO2RateTempo, "Rate"), kNoTag, "LFO2")->Hide(true);
    pGraphics->AttachControl(new IVKnobControl(lfoPanel.GetGridCell(1, 2, 3, 5).GetCentredInside(dialSize), kParamLFO2Amp, "Amp"), kNoTag, "LFO2");
    pGraphics->AttachControl(new ICaptionControl(lfoPanel.GetGridCell(2, 1, 6, 5).GetVPadded(-20.f), kParamLFO2Shape), kNoTag, "LFO2")->DisablePrompt(false);
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(4, 1, 9, 5).GetVPadded(-6.f), kParamLFO2RateMode, "", DEFAULT_STYLE.WithShowLabel(false), "Sync Off", "Sync On"),
      kNoTag, "LFO2")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
        bool sync = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamLFO2RateHz, sync);
        pGraphics->HideControl(kParamLFO2RateTempo, !sync);
        pGraphics->SetAllControlsDirty();
        });
    pGraphics->AttachControl(new IVToggleControl(lfoPanel.GetGridCell(5, 1, 9, 5).GetVPadded(-6.f).GetVShifted(-15.f), kParamLFO2Restart, "Restart", DEFAULT_STYLE.WithShowLabel(false), "Free Run", "Retrigger"), kNoTag, "LFO2");
    pGraphics->AttachControl(new IVDisplayControl(lfoPanel.GetGridCell(1, 3, 3, 5).GetVPadded(-20.f)));
    pGraphics->AttachControl(new ModSliderControl(lfoPanel.GetGridCell(1, 4, 3, 5).GetVPadded(-15.f), kNoParameter, "Depth"), kCtrlTagLFO2Depth, "LFO2");
    pGraphics->AttachControl(new IVGroupControl("LFO 2", "LFO2", 0.f, 10.f, 0.f, 0.f));

    /* Sequencer */
    const IRECT sequencerDisplayPanel = lowerPanel.GetGridCell(2, 2, 3, 3).GetFromLeft(PLUG_WIDTH - lfoPanel.L).GetPadded(0.f, 10.f, 0.f, 0.f);
    const IRECT sequencerParamsPanel = controls.GetFromBottom(PLUG_HEIGHT - lfoPanel.T).GetFromRight(PLUG_WIDTH - lfoPanel.R).GetReducedFromBottom(lowerPanel.H() / 3).GetPadded(-20.f, -15.f, -20.f, -15.f);

    // Sequencer Action Functions
    auto stepRandomizeAction = [pGraphics](IControl* pCaller) {
      SplashClickActionFunc(pCaller);
      pGraphics->GetDelegate()->SendArbitraryMsgFromUI(kMsgRandomizeSequencer);
    };

    pGraphics->AttachControl(new IVKnobControl(sequencerParamsPanel.GetGridCell(0, 0, 2, 3).GetCentredInside(dialSize), kParamSequencerSteps, "Steps"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new IVKnobControl(sequencerParamsPanel.GetGridCell(0, 1, 2, 3).GetCentredInside(dialSize), kParamSequencerRateHz, "Rate"), kNoTag, "SequencerParams")->Hide(true);
    pGraphics->AttachControl(new IVKnobControl(sequencerParamsPanel.GetGridCell(0, 1, 2, 3).GetCentredInside(dialSize), kParamSequencerRateTempo, "Rate"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new IVKnobControl(sequencerParamsPanel.GetGridCell(0, 2, 2, 3).GetCentredInside(dialSize), kParamSequencerAmp, "Amp"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new IVKnobControl(sequencerParamsPanel.GetGridCell(1, 0, 2, 3).GetCentredInside(dialSize), kParamSequencerGlide, "Glide"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(4, 1, 7, 3).GetHPadded(-5.f), kParamSequencerRateMode, "", DEFAULT_STYLE.WithShowLabel(false), "Sync Off", "Sync On"),
      kNoTag, "Sequencer")->SetAnimationEndActionFunction([pGraphics](IControl* pControl) {
        bool sync = pControl->GetValue() > 0.5;
        pGraphics->HideControl(kParamSequencerRateTempo, !sync);
        pGraphics->HideControl(kParamSequencerRateHz, sync);
        pGraphics->SetAllControlsDirty();
        });
    pGraphics->AttachControl(new SequencerControl<kNumSeqSteps>(sequencerDisplayPanel.GetReducedFromRight(dialSize).GetVPadded(-10.f), ""),
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
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(5, 1, 7, 3).GetHPadded(-5.f), kParamSequencerRestart, "", DEFAULT_STYLE.WithShowLabel(false), "Free Run", "Restart"), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new ModSliderControl(sequencerDisplayPanel.GetGridCell(0, 5, 1, 6).GetVPadded(-2.f).GetVShifted(-5.f).GetMidHPadded(20.f), kNoParameter, "Depth"), kCtrlTagSequencerDepth, "Sequencer");
    pGraphics->AttachControl(new IVToggleControl(sequencerParamsPanel.GetGridCell(4, 2, 7, 3).GetHPadded(-5.f), kNoParameter, "", DEFAULT_STYLE.WithShowLabel(false), "Quant.", "Cont."), kNoTag, "SequencerParams")->SetActionFunction(
      [pGraphics](IControl* pControl) {
        int nSteps = pControl->GetValue() > 0.5 ? 12 : 0;
        dynamic_cast<IVTrackControlBase*>(pGraphics->GetControlWithTag(kCtrlTagSequencer))->SetNSteps(nSteps);
      });
    pGraphics->AttachControl(new IVButtonControl(sequencerParamsPanel.GetGridCell(5, 2, 7, 3).GetHPadded(-5.f), stepRandomizeAction, "Randomize", DEFAULT_STYLE.WithLabelText(IText(12.f))), kNoTag, "SequencerParams");
    pGraphics->AttachControl(new IVGroupControl("Sequencer", "SequencerParams", -5.f, 15.f, -5.f, 5.f));

    
#ifndef VST3_API
    pGraphics->AttachControl(new IVButtonControl(keyboardBounds.GetFromTRHC(200, 30).GetTranslated(0, -30.f), SplashClickActionFunc,
      "Show/Hide Keyboard", DEFAULT_STYLE.WithColor(kFG, COLOR_WHITE).WithLabelText({15.f, EVAlign::Middle})))->SetAnimationEndActionFunction(
      [pGraphics, controls](IControl* pCaller) {
        static bool hide = false;
        pGraphics->GetControlWithTag(kCtrlTagKeyboard)->Hide(hide = !hide);
        pGraphics->Resize((float)PLUG_WIDTH, hide ? (float)PLUG_HEIGHT : (float)PLUG_HEIGHT + (float)KEYBOARD_HEIGHT, pGraphics->GetDrawScale());
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
  startPos = chunk.Get(&v, startPos);
  mDSP.mLoadedWavetables[0] = v;
  startPos = chunk.Get(&v, startPos);
  mDSP.mLoadedWavetables[1] = v;

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

  SendCurrentParamValuesFromDelegate();
  GetUI()->SetAllControlsDirty();
}

#if IPLUG_DSP
void Tablitsa::ProcessBlock(sample** inputs, sample** outputs, int nFrames)
{
  mDSP.ProcessBlock(nullptr, outputs, 2, nFrames, mTimeInfo.mPPQPos, mTimeInfo.mTransportIsRunning);
  mMeterSender.ProcessBlock(outputs, nFrames, kCtrlTagMeter);
  // mLFOVisSender.PushData({kCtrlTagLFOVis, {float(mDSP.mLFO.GetLastOutput())}});
  // TODO: If no voices are active, increment the phases of free-run modulators
}

void Tablitsa::OnIdle()
{
  mMeterSender.TransmitData(*this);

  // Send time signature to DSP class
  GetTimeSig(mDSP.mTSNum, mDSP.mTSDenom);
  mDSP.mTempo = GetTempo();

  if (mDSP.mStepPos != mDSP.mPrevPos)
  {
    mDSP.mPrevPos = mDSP.mStepPos;
    SendControlMsgFromDelegate(kCtrlTagSequencer, IVMultiSliderControl<>::kMsgTagSetHighlight, sizeof(int), &mDSP.mPrevPos);
  }

  //mLFOVisSender.TransmitData(*this);
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
  mDSP.SetParam(paramIdx, GetParam(paramIdx)->Value());
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
      case kMsgWavetable1Changed:
      {
        const double wtIdx = *static_cast<const double*>(pData); //Semi-Normalized value (between 1/118 and 1)
        mDSP.mLoadedWavetables[0] = wtIdx * 118.; // Corresponds to the atomic number (not the 0-indexed vector item)
        mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[0]) - 1, 0);
        break;
      }
      case kMsgWavetable2Changed:
      {
        const double wtIdx = *static_cast<const double*>(pData);
        mDSP.mLoadedWavetables[1] = wtIdx * 118.;
        mDSP.UpdateOscillatorWavetable(static_cast<int>(mDSP.mLoadedWavetables[1]) - 1, 1);
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
      default:
        break;
    }
  }
  return false;
}
#endif
