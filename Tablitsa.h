#pragma once

#define DEBUG_VECTOR 0

#if !_DEBUG || DEBUG_VECTOR
#define VECTOR
//#define VECTOR_VOICE_EFFECTS_TEST
//#define VECTOR_MASTER_EFFECTS_TEST
#endif

#include "IPlug_include_in_plug_hdr.h"
#include "IControls.h"

#define PRESET_NAME_CHAR_LENGTH 32
#define TABLITSA_MAX_VOICE_EFFECTS 3
#define TABLITSA_MAX_MASTER_EFFECTS 3
#define TABLITSA_VOICE_EFFECTS_LIST {"None", "Sample & Hold", "Texturizer", "Distortion", "Super Ring"}
#define TABLITSA_MASTER_EFFECTS_LIST {"None", "Delay", "EQ", "Reverb 1", "Reverb 2"}

const int kNumPresets = 1;
constexpr int kNumVoices = 16;
constexpr int kMaxUnisonVoices = 8;
constexpr int kNumSeqSteps = 16;

enum EParams
{
  kParamGain = 0,
  kParamPan,
  kParamPanEnv1,
  kParamPanEnv2,
  kParamPanAmpEnv,
  kParamPanLFO1,
  kParamPanLFO2,
  kParamPanSeq,
  kParamPanVel,
  kParamPanKTk,
  kParamPanRnd,
  kParamNoteGlideTime,
  kParamNoteGlideRate,
  kParamPortamentoMode,
  kParamMonophonic,
  kParamLegato,
  kParamUnisonVoices,
  kParamUnisonDetune,
  kParamUnisonChord,
  kParamStereoSpread,
  kParamAmpEnvAttack,
  kParamAmpEnvDecay,
  kParamAmpEnvSustain,
  kParamAmpEnvRelease,
  kParamAmpEnvVelocity,
  kParamAmpEnvSustainEnv1,
  kParamAmpEnvSustainEnv2,
  kParamAmpEnvSustainAmpEnv,
  kParamAmpEnvSustainLFO1,
  kParamAmpEnvSustainLFO2,
  kParamAmpEnvSustainSeq,
  kParamAmpEnvSustainVel,
  kParamAmpEnvSustainKTk,
  kParamAmpEnvSustainRnd,
  kParamAmpEnvDecayCurve,
  kParamAmpEnvReleaseCurve,
  kParamEnv1Attack,
  kParamEnv1Decay,
  kParamEnv1Sustain,
  kParamEnv1Release,
  kParamEnv1Velocity,
  kParamEnv1SustainEnv1,
  kParamEnv1SustainEnv2,
  kParamEnv1SustainAmpEnv,
  kParamEnv1SustainLFO1,
  kParamEnv1SustainLFO2,
  kParamEnv1SustainSeq,
  kParamEnv1SustainVel,
  kParamEnv1SustainKTk,
  kParamEnv1SustainRnd,
  kParamEnv1DecayCurve,
  kParamEnv1ReleaseCurve,
  kParamEnv1Level,
  kParamEnv2Attack,
  kParamEnv2Decay,
  kParamEnv2Sustain,
  kParamEnv2Release,
  kParamEnv2Velocity,
  kParamEnv2SustainEnv1,
  kParamEnv2SustainEnv2,
  kParamEnv2SustainAmpEnv,
  kParamEnv2SustainLFO1,
  kParamEnv2SustainLFO2,
  kParamEnv2SustainSeq,
  kParamEnv2SustainVel,
  kParamEnv2SustainKTk,
  kParamEnv2SustainRnd,
  kParamEnv2DecayCurve,
  kParamEnv2ReleaseCurve,
  kParamLFO1Shape,
  kParamLFO1Phase,
  kParamLFO1RateHz,
  kParamLFO1RateHzEnv1,
  kParamLFO1RateHzEnv2,
  kParamLFO1RateHzAmpEnv,
  kParamLFO1RateHzLFO1,
  kParamLFO1RateHzLFO2,
  kParamLFO1RateHzSeq,
  kParamLFO1RateHzVel,
  kParamLFO1RateHzKTk,
  kParamLFO1RateHzRnd,
  kParamLFO1RateTempo,
  kParamLFO1Amp,
  kParamLFO1AmpEnv1,
  kParamLFO1AmpEnv2,
  kParamLFO1AmpAmpEnv,
  kParamLFO1AmpLFO1,
  kParamLFO1AmpLFO2,
  kParamLFO1AmpSeq,
  kParamLFO1AmpVel,
  kParamLFO1AmpKTk,
  kParamLFO1AmpRnd,
  kParamLFO1RateMode,
  kParamLFO1Restart,
  kParamLFO2Shape,
  kParamLFO2Phase,
  kParamLFO2RateHz,
  kParamLFO2RateHzEnv1,
  kParamLFO2RateHzEnv2,
  kParamLFO2RateHzAmpEnv,
  kParamLFO2RateHzLFO1,
  kParamLFO2RateHzLFO2, 
  kParamLFO2RateHzSeq,
  kParamLFO2RateHzVel,
  kParamLFO2RateHzKTk,
  kParamLFO2RateHzRnd,
  kParamLFO2RateTempo,
  kParamLFO2Amp,
  kParamLFO2AmpEnv1,
  kParamLFO2AmpEnv2,
  kParamLFO2AmpAmpEnv,
  kParamLFO2AmpLFO1,
  kParamLFO2AmpLFO2,
  kParamLFO2AmpSeq,
  kParamLFO2AmpVel,
  kParamLFO2AmpKTk,
  kParamLFO2AmpRnd,
  kParamLFO2RateMode,
  kParamLFO2Restart,
  kParamSequencerSteps,
  kParamSequencerStepMode,
  kParamSequencerGlide,
  kParamSequencerCurve,
  kParamSequencerRateHz,
  kParamSequencerRateHzEnv1,
  kParamSequencerRateHzEnv2,
  kParamSequencerRateHzAmpEnv,
  kParamSequencerRateHzLFO1,
  kParamSequencerRateHzLFO2,
  kParamSequencerRateHzSeq, // Not used; just here for indexing purposes
  kParamSequencerRateHzVel,
  kParamSequencerRateHzKTk,
  kParamSequencerRateHzRnd,
  kParamSequencerRateTempo,
  kParamSequencerRateMode,
  kParamSequencerRestart,
  kParamSequencerAmp,
  kParamSequencerAmpEnv1,
  kParamSequencerAmpEnv2,
  kParamSequencerAmpAmpEnv,
  kParamSequencerAmpLFO1,
  kParamSequencerAmpLFO2,
  kParamSequencerAmpSeq, // Not used; just here for indexing purposes
  kParamSequencerAmpVel,
  kParamSequencerAmpKTk,
  kParamSequencerAmpRnd,
  kParamWavetable1Pitch, // Voice-specific (polyphonic) parameters
  kParamWavetable1PitchEnv1,
  kParamWavetable1PitchEnv2,
  kParamWavetable1PitchAmpEnv,
  kParamWavetable1PitchLFO1,
  kParamWavetable1PitchLFO2,
  kParamWavetable1PitchSeq,
  kParamWavetable1PitchVel,
  kParamWavetable1PitchKTk,
  kParamWavetable1PitchRnd,
  kParamWavetable1Pos,
  kParamWavetable1PosEnv1,
  kParamWavetable1PosEnv2,
  kParamWavetable1PosAmpEnv,
  kParamWavetable1PosLFO1,
  kParamWavetable1PosLFO2,
  kParamWavetable1PosSeq,
  kParamWavetable1PosVel,
  kParamWavetable1PosKTk,
  kParamWavetable1PosRnd,
  kParamWavetable1Bend,
  kParamWavetable1BendEnv1,
  kParamWavetable1BendEnv2,
  kParamWavetable1BendAmpEnv,
  kParamWavetable1BendLFO1,
  kParamWavetable1BendLFO2,
  kParamWavetable1BendSeq,
  kParamWavetable1BendVel,
  kParamWavetable1BendKTk,
  kParamWavetable1BendRnd,
  kParamWavetable1Amp,
  kParamWavetable1AmpEnv1,
  kParamWavetable1AmpEnv2,
  kParamWavetable1AmpAmpEnv,
  kParamWavetable1AmpLFO1,
  kParamWavetable1AmpLFO2,
  kParamWavetable1AmpSeq,
  kParamWavetable1AmpVel,
  kParamWavetable1AmpKTk,
  kParamWavetable1AmpRnd,
  kParamWavetable1Formant,
  kParamWavetable1FormantEnv1,
  kParamWavetable1FormantEnv2,
  kParamWavetable1FormantAmpEnv,
  kParamWavetable1FormantLFO1,
  kParamWavetable1FormantLFO2,
  kParamWavetable1FormantSeq,
  kParamWavetable1FormantVel,
  kParamWavetable1FormantKTk,
  kParamWavetable1FormantRnd,
  kParamWavetable2Pitch,
  kParamWavetable2PitchEnv1,
  kParamWavetable2PitchEnv2,
  kParamWavetable2PitchAmpEnv,
  kParamWavetable2PitchLFO1,
  kParamWavetable2PitchLFO2,
  kParamWavetable2PitchSeq,
  kParamWavetable2PitchVel,
  kParamWavetable2PitchKTk,
  kParamWavetable2PitchRnd,
  kParamWavetable2Pos,
  kParamWavetable2PosEnv1,
  kParamWavetable2PosEnv2,
  kParamWavetable2PosAmpEnv,
  kParamWavetable2PosLFO1,
  kParamWavetable2PosLFO2,
  kParamWavetable2PosSeq,
  kParamWavetable2PosVel,
  kParamWavetable2PosKTk,
  kParamWavetable2PosRnd,
  kParamWavetable2Bend,
  kParamWavetable2BendEnv1,
  kParamWavetable2BendEnv2,
  kParamWavetable2BendAmpEnv,
  kParamWavetable2BendLFO1,
  kParamWavetable2BendLFO2,
  kParamWavetable2BendSeq,
  kParamWavetable2BendVel,
  kParamWavetable2BendKTk,
  kParamWavetable2BendRnd,
  kParamWavetable2Amp,
  kParamWavetable2AmpEnv1,
  kParamWavetable2AmpEnv2,
  kParamWavetable2AmpAmpEnv,
  kParamWavetable2AmpLFO1,
  kParamWavetable2AmpLFO2,
  kParamWavetable2AmpSeq,
  kParamWavetable2AmpVel,
  kParamWavetable2AmpKTk,
  kParamWavetable2AmpRnd,
  kParamWavetable2Formant,
  kParamWavetable2FormantEnv1,
  kParamWavetable2FormantEnv2,
  kParamWavetable2FormantAmpEnv,
  kParamWavetable2FormantLFO1,
  kParamWavetable2FormantLFO2,
  kParamWavetable2FormantSeq,
  kParamWavetable2FormantVel,
  kParamWavetable2FormantKTk,
  kParamWavetable2FormantRnd,
  kParamFilter1Cutoff, // Filter 1
  kParamFilter1CutoffEnv1,
  kParamFilter1CutoffEnv2,
  kParamFilter1CutoffAmpEnv,
  kParamFilter1CutoffLFO1,
  kParamFilter1CutoffLFO2,
  kParamFilter1CutoffSeq,
  kParamFilter1CutoffVel,
  kParamFilter1CutoffKTk,
  kParamFilter1CutoffRnd,
  kParamFilter1Resonance,
  kParamFilter1ResonanceEnv1,
  kParamFilter1ResonanceEnv2,
  kParamFilter1ResonanceAmpEnv,
  kParamFilter1ResonanceLFO1,
  kParamFilter1ResonanceLFO2,
  kParamFilter1ResonanceSeq,
  kParamFilter1ResonanceVel,
  kParamFilter1ResonanceKTk,
  kParamFilter1ResonanceRnd,
  kParamFilter1Drive,
  kParamFilter1DriveEnv1,
  kParamFilter1DriveEnv2,
  kParamFilter1DriveAmpEnv,
  kParamFilter1DriveLFO1,
  kParamFilter1DriveLFO2,
  kParamFilter1DriveSeq,
  kParamFilter1DriveVel,
  kParamFilter1DriveKTk,
  kParamFilter1DriveRnd,
  kParamFilter1FF,
  kParamFilter1FB,
  kParamFilter1Delay,
  kParamFilter1DelayEnv1,
  kParamFilter1DelayEnv2,
  kParamFilter1DelayAmpEnv,
  kParamFilter1DelayLFO1,
  kParamFilter1DelayLFO2,
  kParamFilter1DelaySeq,
  kParamFilter1DelayVel,
  kParamFilter1DelayKTk,
  kParamFilter1DelayRnd,
  kParamFilter1Type,
  kParamFilter1ModeVSF,
  kParamFilter1ModeMoog,
  kParamFilter1ModeComb,
  kParamFilter1Osc1Send,
  kParamFilter1Osc2Send, // !Filter 1
  kParamFilter2Cutoff, // Filter 2
  kParamFilter2CutoffEnv1,
  kParamFilter2CutoffEnv2,
  kParamFilter2CutoffAmpEnv,
  kParamFilter2CutoffLFO1,
  kParamFilter2CutoffLFO2,
  kParamFilter2CutoffSeq,
  kParamFilter2CutoffVel,
  kParamFilter2CutoffKTk,
  kParamFilter2CutoffRnd,
  kParamFilter2Resonance,
  kParamFilter2ResonanceEnv1,
  kParamFilter2ResonanceEnv2,
  kParamFilter2ResonanceAmpEnv,
  kParamFilter2ResonanceLFO1,
  kParamFilter2ResonanceLFO2,
  kParamFilter2ResonanceSeq,
  kParamFilter2ResonanceVel,
  kParamFilter2ResonanceKTk,
  kParamFilter2ResonanceRnd,
  kParamFilter2Drive,
  kParamFilter2DriveEnv1,
  kParamFilter2DriveEnv2,
  kParamFilter2DriveAmpEnv,
  kParamFilter2DriveLFO1,
  kParamFilter2DriveLFO2,
  kParamFilter2DriveSeq,
  kParamFilter2DriveVel,
  kParamFilter2DriveKTk,
  kParamFilter2DriveRnd,
  kParamFilter2FF,
  kParamFilter2FB,
  kParamFilter2Delay,
  kParamFilter2DelayEnv1,
  kParamFilter2DelayEnv2,
  kParamFilter2DelayAmpEnv,
  kParamFilter2DelayLFO1,
  kParamFilter2DelayLFO2,
  kParamFilter2DelaySeq,
  kParamFilter2DelayVel,
  kParamFilter2DelayKTk,
  kParamFilter2DelayRnd,
  kParamFilter2Type,
  kParamFilter2ModeVSF,
  kParamFilter2ModeMoog,
  kParamFilter2ModeComb,
  kParamFilter2Osc1Send,
  kParamFilter2Osc2Send, // !Filter 2
  kParamOscModulator,
  kParamOsc1PM,
  kParamOsc1RM,
  kParamOsc2PM,
  kParamOsc2RM,
  kParamPhaseModFreq, // !Phase Mod
  kParamPhaseModFreqEnv1,
  kParamPhaseModFreqEnv2,
  kParamPhaseModFreqAmpEnv,
  kParamPhaseModFreqLFO1,
  kParamPhaseModFreqLFO2,
  kParamPhaseModFreqSeq,
  kParamPhaseModFreqVel,
  kParamPhaseModFreqKTk,
  kParamPhaseModFreqRnd,
  kParamPhaseModAmount,
  kParamPhaseModAmountEnv1,
  kParamPhaseModAmountEnv2,
  kParamPhaseModAmountAmpEnv,
  kParamPhaseModAmountLFO1,
  kParamPhaseModAmountLFO2,
  kParamPhaseModAmountSeq,
  kParamPhaseModAmountVel,
  kParamPhaseModAmountKTk,
  kParamPhaseModAmountRnd, // !Phase Mod
  kParamRingModFreq, // Ring Mod
  kParamRingModFreqEnv1,
  kParamRingModFreqEnv2,
  kParamRingModFreqAmpEnv,
  kParamRingModFreqLFO1,
  kParamRingModFreqLFO2,
  kParamRingModFreqSeq,
  kParamRingModFreqVel,
  kParamRingModFreqKTk,
  kParamRingModFreqRnd,
  kParamRingModAmount,
  kParamRingModAmountEnv1,
  kParamRingModAmountEnv2,
  kParamRingModAmountAmpEnv,
  kParamRingModAmountLFO1,
  kParamRingModAmountLFO2,
  kParamRingModAmountSeq,
  kParamRingModAmountVel,
  kParamRingModAmountKTk,
  kParamRingModAmountRnd, // !Ring Mod
  kParamVoiceEffect1Param1, // Voice Effects
  kParamVoiceEffect1Param1Env1,
  kParamVoiceEffect1Param1Env2,
  kParamVoiceEffect1Param1AmpEnv,
  kParamVoiceEffect1Param1LFO1,
  kParamVoiceEffect1Param1LFO2,
  kParamVoiceEffect1Param1Seq,
  kParamVoiceEffect1Param1Vel,
  kParamVoiceEffect1Param1KTk,
  kParamVoiceEffect1Param1Rnd,
  kParamVoiceEffect1Param2,
  kParamVoiceEffect1Param2Env1,
  kParamVoiceEffect1Param2Env2,
  kParamVoiceEffect1Param2AmpEnv,
  kParamVoiceEffect1Param2LFO1,
  kParamVoiceEffect1Param2LFO2,
  kParamVoiceEffect1Param2Seq,
  kParamVoiceEffect1Param2Vel,
  kParamVoiceEffect1Param2KTk,
  kParamVoiceEffect1Param2Rnd,
  kParamVoiceEffect1Param3,
  kParamVoiceEffect1Param3Env1,
  kParamVoiceEffect1Param3Env2,
  kParamVoiceEffect1Param3AmpEnv,
  kParamVoiceEffect1Param3LFO1,
  kParamVoiceEffect1Param3LFO2,
  kParamVoiceEffect1Param3Seq,
  kParamVoiceEffect1Param3Vel,
  kParamVoiceEffect1Param3KTk,
  kParamVoiceEffect1Param3Rnd,
  kParamVoiceEffect1Param4,
  kParamVoiceEffect1Param4Env1,
  kParamVoiceEffect1Param4Env2,
  kParamVoiceEffect1Param4AmpEnv,
  kParamVoiceEffect1Param4LFO1,
  kParamVoiceEffect1Param4LFO2,
  kParamVoiceEffect1Param4Seq,
  kParamVoiceEffect1Param4Vel,
  kParamVoiceEffect1Param4KTk,
  kParamVoiceEffect1Param4Rnd, // !Voice Effect 1
  kParamVoiceEffect1Param5,
  kParamVoiceEffect1Param6,
  kParamVoiceEffect2Param1, // Voice Effect 2
  kParamVoiceEffect2Param1Env1, 
  kParamVoiceEffect2Param1Env2,
  kParamVoiceEffect2Param1AmpEnv,
  kParamVoiceEffect2Param1LFO1,
  kParamVoiceEffect2Param1LFO2,
  kParamVoiceEffect2Param1Seq,
  kParamVoiceEffect2Param1Vel,
  kParamVoiceEffect2Param1KTk,
  kParamVoiceEffect2Param1Rnd,
  kParamVoiceEffect2Param2,
  kParamVoiceEffect2Param2Env1,
  kParamVoiceEffect2Param2Env2,
  kParamVoiceEffect2Param2AmpEnv,
  kParamVoiceEffect2Param2LFO1,
  kParamVoiceEffect2Param2LFO2,
  kParamVoiceEffect2Param2Seq,
  kParamVoiceEffect2Param2Vel,
  kParamVoiceEffect2Param2KTk,
  kParamVoiceEffect2Param2Rnd,
  kParamVoiceEffect2Param3,
  kParamVoiceEffect2Param3Env1,
  kParamVoiceEffect2Param3Env2,
  kParamVoiceEffect2Param3AmpEnv,
  kParamVoiceEffect2Param3LFO1,
  kParamVoiceEffect2Param3LFO2,
  kParamVoiceEffect2Param3Seq,
  kParamVoiceEffect2Param3Vel,
  kParamVoiceEffect2Param3KTk,
  kParamVoiceEffect2Param3Rnd,
  kParamVoiceEffect2Param4,
  kParamVoiceEffect2Param4Env1,
  kParamVoiceEffect2Param4Env2,
  kParamVoiceEffect2Param4AmpEnv,
  kParamVoiceEffect2Param4LFO1,
  kParamVoiceEffect2Param4LFO2,
  kParamVoiceEffect2Param4Seq,
  kParamVoiceEffect2Param4Vel,
  kParamVoiceEffect2Param4KTk,
  kParamVoiceEffect2Param4Rnd,
  kParamVoiceEffect2Param5,
  kParamVoiceEffect2Param6, // !Voice Effect 2
  kParamVoiceEffect3Param1, // Voice Effect 3
  kParamVoiceEffect3Param1Env1,
  kParamVoiceEffect3Param1Env2,
  kParamVoiceEffect3Param1AmpEnv,
  kParamVoiceEffect3Param1LFO1,
  kParamVoiceEffect3Param1LFO2,
  kParamVoiceEffect3Param1Seq,
  kParamVoiceEffect3Param1Vel,
  kParamVoiceEffect3Param1KTk,
  kParamVoiceEffect3Param1Rnd,
  kParamVoiceEffect3Param2,
  kParamVoiceEffect3Param2Env1,
  kParamVoiceEffect3Param2Env2,
  kParamVoiceEffect3Param2AmpEnv,
  kParamVoiceEffect3Param2LFO1,
  kParamVoiceEffect3Param2LFO2,
  kParamVoiceEffect3Param2Seq,
  kParamVoiceEffect3Param2Vel,
  kParamVoiceEffect3Param2KTk,
  kParamVoiceEffect3Param2Rnd,
  kParamVoiceEffect3Param3,
  kParamVoiceEffect3Param3Env1,
  kParamVoiceEffect3Param3Env2,
  kParamVoiceEffect3Param3AmpEnv,
  kParamVoiceEffect3Param3LFO1,
  kParamVoiceEffect3Param3LFO2,
  kParamVoiceEffect3Param3Seq,
  kParamVoiceEffect3Param3Vel,
  kParamVoiceEffect3Param3KTk,
  kParamVoiceEffect3Param3Rnd,
  kParamVoiceEffect3Param4,
  kParamVoiceEffect3Param4Env1,
  kParamVoiceEffect3Param4Env2,
  kParamVoiceEffect3Param4AmpEnv,
  kParamVoiceEffect3Param4LFO1,
  kParamVoiceEffect3Param4LFO2,
  kParamVoiceEffect3Param4Seq,
  kParamVoiceEffect3Param4Vel,
  kParamVoiceEffect3Param4KTk,
  kParamVoiceEffect3Param4Rnd,
  kParamVoiceEffect3Param5,
  kParamVoiceEffect3Param6, // !Voice Effect 3
  kParamMasterEffect1Param1, // Master Effects
  kParamMasterEffect1Param2,
  kParamMasterEffect1Param3,
  kParamMasterEffect1Param4,
  kParamMasterEffect1Param5,
  kParamMasterEffect1Param6,
  kParamMasterEffect2Param1,
  kParamMasterEffect2Param2,
  kParamMasterEffect2Param3,
  kParamMasterEffect2Param4,
  kParamMasterEffect2Param5,
  kParamMasterEffect2Param6,
  kParamMasterEffect3Param1,
  kParamMasterEffect3Param2,
  kParamMasterEffect3Param3,
  kParamMasterEffect3Param4,
  kParamMasterEffect3Param5,
  kParamMasterEffect3Param6,
  kParamOsc1FilterBypass,
  kParamOsc2FilterBypass,
  kParamOsc1EffectBypass,
  kParamOsc2EffectBypass,
  kNumParams
};

constexpr int kNumVoiceEffectParams = kParamVoiceEffect2Param1 - kParamVoiceEffect1Param1;
constexpr int kNumMasterEffectParams = kParamMasterEffect2Param1 - kParamMasterEffect1Param1;

#if IPLUG_DSP
// will use EParams in Tablitsa_DSP.h
#include "Tablitsa_DSP.h"
#endif

enum EControlTags
{
  kCtrlTagPeriodicTable = 0,
  kCtrlTagMeter,
  kCtrlTagLFO1Vis,
  kCtrlTagLFO2Vis,
  kCtrlTagScope,
  kCtrlTagRTText,
  kCtrlTagKeyboard,
  kCtrlTagBender,
  kCtrlTagGlideMode,
  kCtrlTagEnv1Depth, // Modulator depths
  kCtrlTagEnv2Depth,
  kCtrlTagAmpEnvDepth,
  kCtrlTagLFO1Depth,
  kCtrlTagLFO2Depth,
  kCtrlTagSequencerDepth,
  kCtrlTagVelDepth,
  kCtrlTagKTkDepth,
  kCtrlTagRndDepth, // !Modulator depths
  kCtrlTagLFO1RateMode,
  kCtrlTagLFO2RateMode,
  kCtrlTagSequencerRateMode,
  kCtrlTagSequencerQuant,
  kCtrlTagFilter1Cutoff,
  kCtrlTagFilter1Resonance,
  kCtrlTagFilter1Drive,
  kCtrlTagFilter1FF,
  kCtrlTagFilter1FB,
  kCtrlTagFilter1Delay,
  kCtrlTagFilter1Mode,
  kCtrlTagFilter1Type,
  kCtrlTagFilter1Osc1,
  kCtrlTagFilter1Osc2,
  kCtrlTagFilter2Cutoff,
  kCtrlTagFilter2Resonance,
  kCtrlTagFilter2Drive,
  kCtrlTagFilter2FF,
  kCtrlTagFilter2FB,
  kCtrlTagFilter2Delay,
  kCtrlTagFilter2Mode,
  kCtrlTagFilter2Type,
  kCtrlTagFilter2Osc1,
  kCtrlTagFilter2Osc2,
  kCtrlTagOscModFreq,
  kCtrlTagOscModAmt,
  kCtrlTagOsc1ModSwitch,
  kCtrlTagOsc2ModSwitch,
  kCtrlTagLFO1Plot,
  kCtrlTagLFO2Plot,
  kCtrlTagSequencer,
  kCtrlTagEffectBank,
  kCtrlTagVoiceEffectsList, // Voice effect controls
  kCtrlTagVoiceEffectsSwitch,
  kCtrlTagVoiceEffectsKnob1,
  kCtrlTagVoiceEffectsKnob2,
  kCtrlTagVoiceEffectsKnob3,
  kCtrlTagVoiceEffectsKnob4,
  kCtrlTagVoiceEffectsToggle1,
  kCtrlTagVoiceEffectsToggle2,
  kCtrlTagMasterEffectsList, // Master effect controls
  kCtrlTagMasterEffectsSwitch,
  kCtrlTagMasterEffectsKnob1,
  kCtrlTagMasterEffectsKnob2,
  kCtrlTagMasterEffectsKnob3,
  kCtrlTagMasterEffectsKnob4,
  kCtrlTagMasterEffectsToggle1,
  kCtrlTagMasterEffectsToggle2,
  kCtrlTagDefaultPresetList,
  kNumCtrlTags
};

enum EMsgTags
{
  kMsgWavetable1Changed = 0,
  kMsgWavetable2Changed,
  kMsgSeqSliderChanged,
  kMsgUpdateLFO1Plot,
  kMsgUpdateLFO2Plot,
  kMsgRandomizeSequencer,
  kMsgSavePreset,
  kMsgLoadPreset,
  kMsgLoadDefaultPreset,
  kMsgVoiceEffect1Changed,
  kMsgVoiceEffect2Changed,
  kMsgVoiceEffect3Changed,
  kMsgMasterEffect1Changed,
  kMsgMasterEffect2Changed,
  kMsgMasterEffect3Changed
};

enum EModulators
{
  kEnv1 = 0,
  kEnv2,
  kAmpEnv,
  kLFO1,
  kLFO2,
  kSequencer,
  kVelocity,
  kKeytrack,
  kTriggerRandom,
  kNumMods
};

// This must be in the same order as the effect labels. (If you change the effect order, this is all you need to modify)
enum EVoiceEffectTypes
{
  kNoVoiceEffect=0,
  kSampleAndHoldEffect,
  kTexturizerEffect,
  kDistortionEffect,
  kCMEffect,
  kLimiterEffect,
  kNumVoiceEffectTypes
};

enum EMasterEffectTypes
{
  kNoMasterEffect = 0,
  kDelayEffect,
  kEQEffect,
  kReverbEffect,
  kReverb2Effect,
  kNumMasterEffectTypes
};

constexpr EControlTags kStartupTriggerControls[]{
  kCtrlTagLFO1RateMode,
  kCtrlTagLFO2RateMode,
  kCtrlTagSequencerRateMode,
  kCtrlTagSequencerQuant,
//  kCtrlTagFilter1Mode,
//  kCtrlTagFilter2Mode,
  kCtrlTagFilter1Type,
  kCtrlTagFilter2Type,
  kCtrlTagGlideMode
};

using namespace iplug;
using namespace igraphics;

class Tablitsa final : public Plugin
{
public:
  Tablitsa(const InstanceInfo& info);

#if IPLUG_DSP // http://bit.ly/2S64BDd
public:
  void ProcessBlock(sample** inputs, sample** outputs, int nFrames) override;
  void ProcessMidiMsg(const IMidiMsg& msg) override;
  void OnReset() override;
  void OnParamChange(int paramIdx) override;
  void OnIdle() override;
  bool OnMessage(int msgTag, int ctrlTag, int dataSize, const void* pData) override;
  
  void OnUIOpen() override;
  bool SerializeState(IByteChunk& chunk) const override;
  int UnserializeState(const IByteChunk& chunk, int startPos) override;
  void UpdateUIControls();
  bool LoadWavetables();
  /* implement this and return true to trigger your custom about box, when someone clicks about in the menu of a standalone app or VST3 plugin */
  bool OnHostRequestingAboutBox() override; // See IPlugAPP_dialog.cpp
  /* implement this and return true to trigger your custom help info, when someone clicks help in the menu of a standalone app or VST3 plugin */
  bool OnHostRequestingProductHelp() override;

  IByteChunk LoadPreset(const char* filename="UserPreset", bool isBackup=false);
  void SavePreset(IByteChunk& byteData, const char* filename = "UserPreset", bool isBackup=false);
  void LoadDefaultState();
  int CheckVersion(const IByteChunk& presetData);
  bool ShowLoadErrorMessageBox();

  void SetMasterFXSlot(int slotIdx, int masterEffectIdx) { mCurrentMasterFXSlot = slotIdx; mMasterEffectSlots[slotIdx] = masterEffectIdx; }
  void SetVoiceFXSlot(int slotIdx, int voiceEffectIdx) { mCurrentVoiceFXSlot = slotIdx; mVoiceEffectSlots[slotIdx] = voiceEffectIdx; }

  // Store the status of the tempo sync control, required for setting the correct parameter range during preset loading
  void SetDelayTempoSync(int slotIdx, bool tempoSync) { mDelayTempoSync[slotIdx] = tempoSync; }
  void SetDelayTempoSync(bool tempoSync) { SetDelayTempoSync(mCurrentMasterFXSlot, tempoSync); }

  int GetActiveModIdx() const;
  void SetActiveModIdx(int idx);
  int GetFirstModCtrlTag() const { return kCtrlTagEnv1Depth; }
  int GetLastModCtrlTag() const { return kCtrlTagRndDepth; }

  void RefreshEffectBankControl();

private:
  TablitsaDSP<sample> mDSP {kNumVoices}; // sample is an alias for double
  IPeakSender<2> mMeterSender;
  IControl* mActiveControl{};
  IByteChunk mStateBackup; // Stores the last user state before reseting, so that it may be restored

  // UI Status variables not stored as parameters
  char mPresetName[PRESET_NAME_CHAR_LENGTH]{};
  int mPresetID{};
  int mVoiceEffectSlots[TABLITSA_MAX_VOICE_EFFECTS]{}; // Holds the ID number of the effect in each slot for the voice effects
  int mMasterEffectSlots[TABLITSA_MAX_MASTER_EFFECTS]{}; // Holds the ID number of the effect in each slot for the master effects
  int mCurrentVoiceFXSlot{ 0 }; // The effect slot current open for editing. Controled by the Slide-Switch controls
  int mCurrentMasterFXSlot{ 0 };
  int mCurrentEffectsTab{ 1 };
  bool mDelayTempoSync[TABLITSA_MAX_MASTER_EFFECTS]{};
  bool mMonoDelay[TABLITSA_MAX_MASTER_EFFECTS]{};
  double mSequencerIsQuantized{ 0. };

  int mActiveModIdx{ -1 };
#endif
};

std::string GetDataPath(char* appendPath="\\");

std::vector<char> ReadAllBytes(const char* fname); // For reading preset files of arbitrary length

#include "TablitsaControls.h"