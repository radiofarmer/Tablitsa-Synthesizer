#pragma once

#include "IPlug_include_in_plug_hdr.h"
#include "IControls.h"
#include <vectorclass.h>


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
  kParamAmpEnvAttack,
  kParamAmpEnvDecay,
  kParamAmpEnvSustain,
  kParamAmpEnvRelease,
  kParamAmpEnvVelocity,
  kParamEnv1Attack,
  kParamEnv1Decay,
  kParamEnv1Sustain,
  kParamEnv1Release,
  kParamEnv1Velocity,
  kParamEnv1Level,
  kParamEnv2Attack,
  kParamEnv2Decay,
  kParamEnv2Sustain,
  kParamEnv2Release,
  kParamEnv2Velocity,
  kParamLFO1Shape,
  kParamLFO1RateHz,
  kParamLFO1RateTempo,
  kParamLFO1Amp,
  kParamLFO1RateMode,
  kParamLFO1Restart,
  kParamLFO2Shape,
  kParamLFO2RateHz,
  kParamLFO2RateTempo,
  kParamLFO2Amp,
  kParamLFO2RateMode,
  kParamLFO2Restart,
  kParamSequencerSteps,
  kParamSequencerGlide,
  kParamSequencerCurve,
  kParamSequencerRateHz,
  kParamSequencerRateTempo,
  kParamSequencerRateMode,
  kParamSequencerRestart,
  kParamSequencerAmp,
  kParamWavetable1,
  kParamWavetable2,
  kParamWavetable1Pitch,
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
  kParamWavetable1Sub,
  kParamWavetable1SubEnv1,
  kParamWavetable1SubEnv2,
  kParamWavetable1SubAmpEnv,
  kParamWavetable1SubLFO1,
  kParamWavetable1SubLFO2,
  kParamWavetable1SubSeq,
  kParamWavetable1SubVel,
  kParamWavetable1SubKTk,
  kParamWavetable1SubRnd,
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
  kParamWavetable2Sub,
  kParamWavetable2SubEnv1,
  kParamWavetable2SubEnv2,
  kParamWavetable2SubAmpEnv,
  kParamWavetable2SubLFO1,
  kParamWavetable2SubLFO2,
  kParamWavetable2SubSeq,
  kParamWavetable2SubVel,
  kParamWavetable2SubKTk,
  kParamWavetable2SubRnd,
  kParamFilter1Cutoff,
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
  kParamFilter1Type,
  kParamFilter1ModeVSF,
  kParamFilter1ModeMoog,
  kParamFilter1ModeComb,
  kParamFilter2Cutoff,
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
  kParamFilter2Type,
  kParamFilter2ModeVSF,
  kParamFilter2ModeMoog,
  kParamFilter2ModeComb,
  kParamOscModulator,
  kParamOsc1PM,
  kParamOsc1RM,
  kParamOsc2PM,
  kParamOsc2RM,
  kParamPhaseModFreq,
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
  kParamPhaseModAmountRnd,
  kParamRingModFreq,
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
  kParamRingModAmountRnd,
  kParamDelayTimeLMilliseconds,
  kParamDelayTimeRMilliseconds,
  kParamDelayTimeLBeats,
  kParamDelayTimeRBeats,
  kParamDelayTimeMode,
  kParamDelayFeedback,
  kParamDelayMix,
  kNumParams
};

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
  kCtrlTagEnv1Depth, // Modulator depths
  kCtrlTagEnv2Depth,
  kCtrlTagAmpEnvDepth,
  kCtrlTagLFO1Depth,
  kCtrlTagLFO2Depth,
  kCtrlTagSequencerDepth,
  kCtrlTagVelDepth,
  kCtrlTagKTkDepth,
  kCtrlTagRndDepth, // !Modulator depths
  kCtrlTagFilter1Cutoff,
  kCtrlTagFilter1Resonance,
  kCtrlTagFilter1Drive,
  kCtrlTagFilter1FF,
  kCtrlTagFilter1FB,
  kCtrlTagFilter1Delay,
  kCtrlTagFilter1Mode,
  kCtrlTagFilter1Type,
  kCtrlTagFilter2Cutoff,
  kCtrlTagFilter2Resonance,
  kCtrlTagFilter2Drive,
  kCtrlTagFilter2FF,
  kCtrlTagFilter2FB,
  kCtrlTagFilter2Delay,
  kCtrlTagFilter2Mode,
  kCtrlTagFilter2Type,
  kCtrlTagOscModFreq,
  kCtrlTagOscModAmt,
  kCtrlTagOsc1ModSwitch,
  kCtrlTagOsc2ModSwitch,
  kCtrlTagSequencer,
  kCtrlTagDelayLMilliseconds,
  kCtrlTagDelayRMilliseconds,
  kCtrlTagDelayLBeats,
  kCtrlTagDelayRBeats,
  kNumCtrlTags
};

enum EMsgTags
{
  kMsgWavetable1Changed = 0,
  kMsgWavetable2Changed,
  kMsgSeqSliderChanged,
  kMsgRandomizeSequencer,
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
  
private:
  TablitsaDSP<sample> mDSP {kNumVoices}; // sample is an alias for double
  IPeakSender<2> mMeterSender;
  ISender<1> mLFO1VisSender;
  ISender<1> mLFO2VisSender;
  IControl* mActiveControl{};
#endif
};
