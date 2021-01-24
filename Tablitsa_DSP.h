#pragma once

#include "Oscillator.h"
#include "MidiSynth.h"
#include "ADSREnvelope.h"
#include "Smoothers.h"
#include "Wavetable.h"
#include "LFO.h"
#include "Filter.h"
#include "PeriodicTable.h"

#ifndef _DEBUG
  #define VECTOR
#endif
#ifdef VECTOR
  #ifdef OVERSAMPLING
    #define FRAME_INTERVAL 2
  #else
    #define FRAME_INTERVAL 4
  #endif
#else
  #define FRAME_INTERVAL 1
#endif

using namespace iplug;

/*
Global Modulations (smoothers): These values are computed once per sample and sent to all voices 
*/
enum EModulations
{
  kModGainSmoother = 0,
  kModEnv1SustainSmoother,
  kModEnv2SustainSmoother,
  kModAmpEnvSustainSmoother,
  kModWavetable1PitchSmoother,
  kModWavetable1PosSmoother,
  kModWavetable1BendSmoother,
  kModWavetable1SubSmoother,
  kModWavetable1AmpSmoother,
  kModWavetable2PitchSmoother,
  kModWavetable2PosSmoother,
  kModWavetable2BendSmoother,
  kModWavetable2SubSmoother,
  kModWavetable2AmpSmoother,
  kModFilter1CutoffSmoother,
  kModFilter1ResonanceSmoother,
  kModFilter1DriveSmoother,
  kModFilter2CutoffSmoother,
  kModFilter2ResonanceSmoother,
  kModFilter2DriveSmoother,
  kModPhaseModFreqSmoother,
  kModPhaseModAmtSmoother,
  kModRingModFreqSmoother,
  kModRingModAmtSmoother,
  kNumModulations,
};

/*
Per-Voice Modulations: These values are calculated based on the current value of the modulator object which each voice owns.
They are calculated once per sample per voice.

NOTE: These need to be in the same order as the global modulations (the smoothers) after a certain point (e.g. after the envelope sustain smoothers).
This needs to be fixed eventually so that both sets of modulations operate off the same list.
*/
enum EVoiceModParams
{
  kVWavetable1PitchOffset = 0,
  kVWavetable1Position,
  kVWavetable1Bend,
  kVWavetable1Sub,
  kVWavetable1Amp,
  kVWavetable2PitchOffset,
  kVWavetable2Position,
  kVWavetable2Bend,
  kVWavetable2Sub,
  kVWavetable2Amp,
  kVFilter1Cutoff,
  kVFilter1Resonance,
  kVFilter1Drive,
  kVFilter2Cutoff,
  kVFilter2Resonance,
  kVFilter2Drive,
  kVPhaseModFreq,
  kVPhaseModAmt,
  kVRingModFreq,
  kVRingModAmt,
  kNumVoiceModulations
};

// See Modulators.h for the enumeration of the number of modulators

template<typename T>
class TablitsaDSP
{
public:
#pragma mark - Voice
  class Voice : public SynthVoice
  {
  public:
    Voice(std::vector<std::string> programList, int id=0) : Voice()
    {
    }

    Voice(int id=0)
      : mID(id),
      mAmpEnv("gain", [&]() { mOsc1.Reset(); }),
      mEnv1("env1", [&]() { mOsc1.Reset(); }) // capture ok on RT thread?
    {
//      DBGMSG("new Voice: %i control inputs.\n", static_cast<int>(mInputs.size()));
      mModulators.AddModulator(&mEnv1);
      mModulators.AddModulator(&mEnv2);
      mModulators.AddModulator(&mAmpEnv);
      mModulators.AddModulator(&mLFO1);
      mModulators.AddModulator(&mLFO2);

//      mVParameterModulators.SetList(mVoiceModParams);
    }

    bool GetBusy() const override
    {
      return mAmpEnv.GetBusy();
    }

    void Trigger(double level, bool isRetrigger) override
    {
      mOsc1.Reset();
      mOsc2.Reset();
      for (auto f : mFilters)
      {
        f->Reset();
      }

      // Reset LFOs
      if (mLFO1Restart)
        mLFO1.Reset();

      if (mLFO2Restart)
        mLFO2.Reset();

      if (isRetrigger)
      {
        mAmpEnv.Retrigger(level);
        mEnv1.Retrigger(level);
        mEnv2.Retrigger(level);
      }
      else
      {
        mAmpEnv.Start(level);
        mEnv1.Start(level);
        mEnv2.Start(level);
      }
    }
    
    void Release() override
    {
      mAmpEnv.Release();
      mEnv1.Release();
      mEnv2.Release();
    }

    void ProcessSamplesAccumulating(T** inputs, T** outputs, int nInputs, int nOutputs, int startIdx, int nFrames) override
    {
      // inputs to the synthesizer can just fetch a value every block, like this:
//      double gate = mInputs[kVoiceControlGate].endValue;
      double pitch = mInputs[kVoiceControlPitch].endValue;
//      double pitchBend = mInputs[kVoiceControlPitchBend].endValue;
      // or write the entire control ramp to a buffer, like this, to get sample-accurate ramps:
//      mInputs[kVoiceControlTimbre].Write(mTimbreBuffer.Get(), startIdx, nFrames);

      mModulators.ProcessBlock(&(inputs[kModEnv1SustainSmoother]), nFrames);
      mVoiceModParams.ProcessBlock(&inputs[kModWavetable1PitchSmoother], mModulators.GetList(), mVModulations.GetList(), nFrames);

      const double phaseModFreqFact = pow(2., inputs[kModPhaseModFreqSmoother][0] / 12.);
      const double ringModFreqFact = pow(2., inputs[kModRingModFreqSmoother][0] / 12.);

      // make sound output for each output channel
      for(auto i = startIdx; i < startIdx + nFrames; i += FRAME_INTERVAL)
      {
        int bufferIdx = i - startIdx;
//        float noise = mTimbreBuffer.Get()[i] * Rand();
        double ampEnvVal{ mModulators.GetList()[2][bufferIdx] }; // Calculated for easy access
        double modVals[]{ mModulators.GetList()[0][bufferIdx],
          mModulators.GetList()[1][bufferIdx],
          ampEnvVal,
          mModulators.GetList()[3][bufferIdx],
          mModulators.GetList()[4][bufferIdx] };
        ParameterModulator::SetModValues(modVals);

        // Oscillator Parameters
        double osc1Freq = 440. * pow(2., pitch + mVModulations.GetList()[kVWavetable1PitchOffset][bufferIdx] / 12.);
        mOsc1.SetWtPosition(1 - mVModulations.GetList()[kVWavetable1Position][bufferIdx]); // Wavetable 1 Position
        mOsc1.SetWtBend(mVModulations.GetList()[kVWavetable1Bend][bufferIdx]); // Wavetable 1 Bend
        mOsc1Sub.SetLevel(mVModulations.GetList()[kVWavetable1Sub][bufferIdx]);
        mOsc1.SetPhaseModulation(mVModulations.GetList()[kVPhaseModAmt][bufferIdx], phaseModFreqFact);
        mOsc1.SetRingModulation(mVModulations.GetList()[kVRingModFreq][bufferIdx], ringModFreqFact);

        double osc2Freq = 440. * pow(2., pitch + mVModulations.GetList()[kVWavetable2PitchOffset][bufferIdx] / 12.);
        mOsc2.SetWtPosition(1 - mVModulations.GetList()[kVWavetable2Position][bufferIdx]); // Wavetable 2 Position
        mOsc2.SetWtBend(mVModulations.GetList()[kVWavetable2Bend][bufferIdx]); // Wavetable 2 Bend
        mOsc2Sub.SetLevel(mVModulations.GetList()[kVWavetable2Sub][bufferIdx]);
        mOsc2.SetPhaseModulation(mVModulations.GetList()[kVPhaseModAmt][bufferIdx], phaseModFreqFact);
        mOsc2.SetRingModulation(mVModulations.GetList()[kVRingModFreq][bufferIdx], ringModFreqFact);
        
        mFilters.at(0)->SetCutoff(mVModulations.GetList()[kVFilter1Cutoff][bufferIdx]); // Filter 1 Cutoff
        mFilters.at(0)->SetQ(mVModulations.GetList()[kVFilter1Resonance][bufferIdx]); // Filter 1 Resonance
        mFilters.at(0)->SetDrive(mVModulations.GetList()[kVFilter1Drive][bufferIdx]); // Filter 1 Drive

        mFilters.at(1)->SetCutoff(mVModulations.GetList()[kVFilter2Cutoff][bufferIdx]); // Filter 2 Cutoff
        mFilters.at(1)->SetQ(mVModulations.GetList()[kVFilter2Resonance][bufferIdx]); // Filter 2 Resonance
        mFilters.at(1)->SetDrive(mVModulations.GetList()[kVFilter2Drive][bufferIdx]); // Filter 2 Drive
        
        // Signal Processing
        std::array<T, OUTPUT_SIZE> osc1Output{ mOsc1.ProcessMultiple(osc1Freq) };
        std::array<T, OUTPUT_SIZE> osc2Output{ mOsc2.ProcessMultiple(osc2Freq) };
        
       for (auto j = 0; j < FRAME_INTERVAL; ++j)
       {
         osc1Output[j] *= mVModulations.GetList()[kVWavetable1Amp][bufferIdx];
         osc2Output[j] *= mVModulations.GetList()[kVWavetable2Amp][bufferIdx];
         osc1Output[j] = mOsc1Sub.Process(osc1Output[j]);
         double filter1Output = mFilters.at(0)->Process(osc1Output[j]);
         double filter2Output = mFilters.at(1)->Process(osc2Output[j]);
         double output_summed = filter1Output + filter2Output;
         outputs[0][i + j] += output_summed * ampEnvVal * mGain;
         outputs[1][i + j] = outputs[0][i + j];
       }
      }
    }

    void SetSampleRateAndBlockSize(double sampleRate, int blockSize) override
    {
      mOsc1.SetSampleRate(sampleRate);
      mOsc2.SetSampleRate(sampleRate);
      mAmpEnv.SetSampleRate(sampleRate);
      mEnv1.SetSampleRate(sampleRate);
      mEnv2.SetSampleRate(sampleRate);
      mLFO1.SetSampleRate(sampleRate);
      mLFO2.SetSampleRate(sampleRate);
      
      mVModulationsData.Resize(blockSize * kNumModulations);
      mVModulations.Empty();

      mModulators.EmptyAndResize(blockSize, kNumMods);

      for (auto i = 0; i < kNumVoiceModulations; i++)
      {
        mVModulations.Add(mVModulationsData.Get() + static_cast<size_t>(blockSize) * i);
      }
      
    }

    void SetProgramNumber(int pgm) override
    {
      //TODO:
    }

    // this is called by the VoiceAllocator to set generic control values.
    void SetControl(int controlNumber, float value) override
    {
      //TODO:
    }

    void SetLFORateCPS(int lfo, double value)
    {
      switch (lfo) {
      case 0:
        mLFO1.SetFreqCPS(value);
        break;
      case 1:
        mLFO2.SetFreqCPS(value);
        break;
      default:
        break;
      };
    }

    void SetLFORateTempo(int lfo, double value)
    {
      switch(lfo) {
      case 0:
        mLFO1.SetQNScalarFromDivision(static_cast<int>(value));
        break;
      case 1:
        mLFO2.SetQNScalarFromDivision(static_cast<int>(value));
        break;
      default:
        break;
      };
    }

    static inline void SetTempoAndBeat(double qnPos = 0., bool transportIsRunning = false, double tempo = 120.)
    {
      mQNPos = qnPos;
      mTransportIsRunning = transportIsRunning;
      mTempo = tempo;
      FastLFO<T>::SetTempoAndBeat(mQNPos, mTransportIsRunning, mTempo);
    }

    void SetFilterType(int filter, int filterType)
    {
      // Indices of cutoff/res/drive for the given filter
      switch (filterType) {
      case kVSF:
      {
        TablitsaDSP<T>::mCombOn = false;
        mVoiceModParams[filter ? kVFilter2Cutoff : kVFilter1Cutoff].SetMinMax(0.001, 0.49);
        mVoiceModParams[filter ? kVFilter2Resonance : kVFilter1Resonance].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Drive : kVFilter1Drive].SetMinMax(0., 1);
        mFilters.at(filter) = new SVF2<T>(mOsc1.GetSampleRate());
        break;
      }
      case kMoog:
      {
        TablitsaDSP<T>::mCombOn = false;
        mVoiceModParams[filter ? kVFilter2Cutoff : kVFilter1Cutoff].SetMinMax(0.001, 0.49);
        mVoiceModParams[filter ? kVFilter2Resonance : kVFilter1Resonance].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Drive : kVFilter1Drive].SetMinMax(0., 1);
        mFilters.at(filter) = new MoogLadder<T>(mOsc1.GetSampleRate());
        break;
      }
      case kComb:
      {
        TablitsaDSP<T>::mCombOn = true;
        mVoiceModParams[filter ? kVFilter2Cutoff : kVFilter1Cutoff].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Resonance : kVFilter1Resonance].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Drive : kVFilter1Drive].SetMinMax(0., (double)COMB_MAX_DELAY);
        mFilters.at(filter) = new CombFilter<T>(mOsc1.GetSampleRate());
        break;
      }
      default:
        mFilters.at(filter) = new NullFilter<T>(mOsc1.GetSampleRate());
        break;

      }
    }

    void UpdateVoiceParam(int voiceParam, int modIdx, double value)
    {
      mVoiceModParams[voiceParam].SetValue(modIdx - 1, value); // Eventually adjust the value of modIdx in the SetParam switch statement rather than here
    }

  public:
    WavetableOscillator<T> mOsc1{ 0 };
    WavetableOscillator<T> mOsc2{ 1 };
    BassBoost<T> mOsc1Sub;
    BassBoost<T> mOsc2Sub;

    ADSREnvelope<T> mEnv1;
    ADSREnvelope<T> mEnv2;
    ADSREnvelope<T> mAmpEnv;
    FastLFO<T> mLFO1;
    FastLFO<T> mLFO2;
    ModulatorList<T, ADSREnvelope, FastLFO> mModulators;

    bool mLFO1Restart{ false };
    bool mLFO2Restart{ false };
    int mFilterUpdateFreq{ 2 };

    std::vector<Filter<T>*> mFilters{ new NullFilter<T>(), new NullFilter<T>() };

    WDL_PtrList<T> mVModulations; // Pointers to modulator buffers
    WDL_TypedBuf<T> mVModulationsData; // Modulator buffer sample data

  private:
//    WDL_TypedBuf<float> mTimbreBuffer;
    ModulatedParameterList<T, kNumVoiceModulations> mVoiceModParams{
      new ParameterModulator(-24., 24.), /* Wavetable 1 Pitch Offset */
      new ParameterModulator(0., 1.), /* Wavetable 1 Position */
      new ParameterModulator(-1., 1.), /* Wavetable 1 Bend */
      new ParameterModulator(0., 1.), /* Wavetable 1 Sub */
      new ParameterModulator(0., 1.), /* Wavetable 1 Amp */
      new ParameterModulator(-24., 24.), /* Wavetable 2 Pitch Offset */
      new ParameterModulator(0., 1.), /* Wavetable 2 Position */
      new ParameterModulator(-1., 1.), /* Wavetable 2 Bend */
      new ParameterModulator(0., 1.), /* Wavetable 2 Sub */
      new ParameterModulator(0., 1.), /* Wavetable 2 Amp */
      new ParameterModulator(0.001, 0.5, true), /* Filter 1 Cutoff */
      new ParameterModulator(0., 1.), /* Filter 1 Resonance */
      new ParameterModulator(0., 1.), /*Filter 1 Drive */
      new ParameterModulator(0.001, 0.5, true), /* Filter 2 Cutoff */
      new ParameterModulator(0., 1.), /* Filter 2 Resonance */
      new ParameterModulator(0., 1.), /*Filter 2 Drive */
      new ParameterModulator(-24., 24.), /* Phase Mod Frequency*/
      new ParameterModulator(0., 1.), /* Phase Mod Depth */
      new ParameterModulator(-24., 24.), /* Ring Mod Frequency*/
      new ParameterModulator(0., 1.) /* Ring Mod Depth */ };

    static inline double mTempo{ 120. };
    static inline bool mTransportIsRunning{ false };
    static inline double mQNPos{ 0. };

    LogParamSmooth<T> mFilter1Smoother{ 10. };
    LogParamSmooth<T> mFilter2Smoother{ 10. };

    int mID;

    // noise generator for test
    uint32_t mRandSeed = 0;

    // return single-precision floating point number on [-1, 1]
    float Rand()
    {
      mRandSeed = mRandSeed * 0x0019660D + 0x3C6EF35F;
      uint32_t temp = ((mRandSeed >> 9) & 0x007FFFFF) | 0x3F800000;
      return (*reinterpret_cast<float*>(&temp))*2.f - 3.f;
    }

  };

/* end Voice class */

public:
#pragma mark -
  TablitsaDSP(int nVoices)
  {
    for (auto i = 0; i < nVoices; i++)
    {
      // add a voice to Zone 0.
      mSynthVoices.push_back(new Voice(i));
      mSynth.AddVoice(mSynthVoices.at(i), 0);
    }
    // some MidiSynth API examples:
    // mSynth.SetKeyToPitchFn([](int k){return (k - 69.)/24.;}); // quarter-tone scale
    // mSynth.SetNoteGlideTime(0.5); // portamento
  }

  void ProcessBlock(T** inputs, T** outputs, int nOutputs, int nFrames, double qnPos = 0., bool transportIsRunning = false, double tempo = 120.)
  {
    // clear outputs
    for(auto i = 0; i < nOutputs; i++)
    {
      memset(outputs[i], 0, nFrames * sizeof(T));
    }

    mParamSmoother.ProcessBlock(mParamsToSmooth, mModulations.GetList(), nFrames); // Populate modulations list (to be sent to mSynth as inputs)
    Voice::SetTempoAndBeat(qnPos, transportIsRunning, tempo);
    mSynth.ProcessBlock(mModulations.GetList(), outputs, 0, nOutputs, nFrames);

    for(int s=0; s < nFrames;s++)
    {
      T smoothedGain = mModulations.GetList()[kModGainSmoother][s];
      outputs[0][s] *= smoothedGain;
      outputs[1][s] *= smoothedGain;
    }
  }

  void Reset(double sampleRate, int blockSize)
  {
    mSampleRate = sampleRate;
    mSynth.SetSampleRateAndBlockSize(sampleRate, blockSize);
    mSynth.Reset();
    mSynth.ForEachVoice([sampleRate](SynthVoice& voice) {
      for (auto f : dynamic_cast<TablitsaDSP::Voice&>(voice).mFilters)
        f->UpdateSampleRate(sampleRate);
      });
    mModulationsData.Resize(blockSize * kNumModulations);
    mModulations.Empty();
    
    for(auto i = 0; i < kNumModulations; i++)
    {
      mModulations.Add(mModulationsData.Get() + static_cast<size_t>(blockSize) * i);
    }
  }

  void ProcessMidiMsg(const IMidiMsg& msg)
  {
    mSynth.AddMidiMsgToQueue(msg);
  }

  void UpdateOscillatorWavetable(int wtIdx, int oscIdx)
  {
    mSynth.Reset();
    tableLoading[oscIdx] = true; // NB: this variable lets the PeriodicTable control know whether to display the selected element in the loading (faded) state or not
    WtFile wtFile{ mWavetables.at(wtIdx) };
    WavetableOscillator<T>::LoadNewTable(wtFile, oscIdx);
    if (oscIdx == 0)
    {
      SendParam([this, oscIdx, &wtFile](Voice* voice) {
        voice->mOsc1.SetWavetable(WavetableOscillator<T>::LoadedTables[oscIdx]);
        voice->mOsc1.ReloadLUT();
        });
      WavetableOscillator<T>::NotifyLoaded(oscIdx);
    }
    else
    {
      SendParam([this, oscIdx, &wtFile](Voice* voice) {
        voice->mOsc2.SetWavetable(WavetableOscillator<T>::LoadedTables[oscIdx]);
        voice->mOsc2.ReloadLUT();
        });
      WavetableOscillator<T>::NotifyLoaded(oscIdx);
    }
    tableLoading[oscIdx] = false;
  }

  inline void SendParam(std::function<void(Voice* voice)> func)
  {
    for (Voice* v : mSynthVoices)
    {
      func(v);
    }
  }

  inline void SendData(std::function<void(Voice* voice)> func)
  {
    for (Voice* v : mSynthVoices)
    {
      func(v);
    }
  }

  void SetParam(int paramIdx, double value)
  {
    using EEnvStage = ADSREnvelope<sample>::EStage;
    
    switch (paramIdx) {
      case kParamNoteGlideTime:
        mSynth.SetNoteGlideTime(value / 1000.);
        break;
      case kParamGain:
        mParamsToSmooth[kModGainSmoother] = (T) value / 100.;
        break;
      case kParamEnv1Sustain:
        mParamsToSmooth[kModEnv1SustainSmoother] = (T)value / 100.;
        break;
      case kParamEnv1Attack:
      case kParamEnv1Decay:
      case kParamEnv1Release:
      {
        EEnvStage stage = static_cast<EEnvStage>(EEnvStage::kAttack + (paramIdx - kParamEnv1Attack));
        mSynth.ForEachVoice([stage, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mEnv1.SetStageTime(stage, value);
          });
        break;
      }
      case kParamEnv2Sustain:
        mParamsToSmooth[kModEnv2SustainSmoother] = (T)value / 100.;
        break;
      case kParamEnv2Attack:
      case kParamEnv2Decay:
      case kParamEnv2Release:
      {
        EEnvStage stage = static_cast<EEnvStage>(EEnvStage::kAttack + (paramIdx - kParamEnv2Attack));
        mSynth.ForEachVoice([stage, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mEnv2.SetStageTime(stage, value);
          });
        break;
      }
      case kParamAmpEnvSustain:
        mParamsToSmooth[kModAmpEnvSustainSmoother] = (T) value / 100.;
        break;
      case kParamAmpEnvAttack:
      case kParamAmpEnvDecay:
      case kParamAmpEnvRelease:
      {
        EEnvStage stage = static_cast<EEnvStage>(EEnvStage::kAttack + (paramIdx - kParamAmpEnvAttack));
        mSynth.ForEachVoice([stage, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mAmpEnv.SetStageTime(stage, value);
        });
        break;
      }
      case kParamLFO1Amp:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO1.SetScalar(value);
          });
        break;
      case kParamLFO1RateTempo:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).SetLFORateTempo(0, value);
          });
        break;
      case kParamLFO1RateHz:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).SetLFORateCPS(0, value);
          });
        break;
      case kParamLFO1RateMode:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO1.SetRateMode(value > 0.5);
          });
        break;
      case kParamLFO1Shape:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO1.SetShape(static_cast<int>(value));
          });
        break;
      case kParamLFO1Restart:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO1Restart = (value > 0.5);
          });
        break;
      case kParamLFO2Amp:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO2.SetScalar(value);
          });
        break;
      case kParamLFO2RateTempo:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).SetLFORateTempo(1, value);
          });
        break;
      case kParamLFO2RateHz:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).SetLFORateCPS(1, value);
          });
        break;
      case kParamLFO2RateMode:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO2.SetRateMode(value > 0.5);
          });
        break;
      case kParamLFO2Shape:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO2.SetShape(static_cast<int>(value));
          });
        break;
      case kParamLFO2Restart:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO2Restart = (value > 0.5);
          });
        break;
      case kParamWavetable1:
      {
        value = std::max(1., value);
        UpdateOscillatorWavetable(static_cast<int>(value) - 1, 0);
        break;
      }
      case kParamWavetable1Pitch:
        mParamsToSmooth[kModWavetable1PitchSmoother] = value;
        break;
      case kParamWavetable1PitchEnv1:
      case kParamWavetable1PitchEnv2:
      case kParamWavetable1PitchLFO1:
      case kParamWavetable1PitchLFO2:
      case kParamWavetable1PitchSeq:
      {
        const int modIdx = paramIdx - kParamWavetable1Pitch;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceParam(kVWavetable1PitchOffset, modIdx, value);
          });
        break;
      }
      case kParamWavetable1Amp:
        mParamsToSmooth[kModWavetable1AmpSmoother] = value;
        break;
      case kParamWavetable1AmpEnv1:
      case kParamWavetable1AmpEnv2:
      case kParamWavetable1AmpLFO1:
      case kParamWavetable1AmpLFO2:
      case kParamWavetable1AmpSeq:
      {
        const int modIdx = paramIdx - kParamWavetable1Amp;
        SendParam([value, modIdx](Voice* voice) {
          voice->UpdateVoiceParam(kVWavetable1Amp, modIdx, value);
          });
        break;
      }
      case kParamWavetable1Pos:
        mParamsToSmooth[kModWavetable1PosSmoother] = (T)value;
        break;
      case kParamWavetable1PosEnv1:
      case kParamWavetable1PosEnv2:
      case kParamWavetable1PosLFO1:
      case kParamWavetable1PosLFO2:
      case kParamWavetable1PosSeq:
      {
        const int modIdx = paramIdx - kParamWavetable1Pos;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceParam(kVWavetable1Position, modIdx, value);
          });
        break;
      }
      case kParamWavetable1Bend:
        mParamsToSmooth[kModWavetable1BendSmoother] = value;
        break;
      case kParamWavetable1BendEnv1:
      case kParamWavetable1BendEnv2:
      case kParamWavetable1BendLFO1:
      case kParamWavetable1BendLFO2:
      case kParamWavetable1BendSeq:
      {
        const int modIdx = paramIdx - kParamWavetable1Bend;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVWavetable1Bend, modIdx, value);
          });
        break;
      }
      case kParamWavetable1Sub:
        mParamsToSmooth[kModWavetable1SubSmoother] = value;
        break;
      case kParamWavetable1SubEnv1:
      case kParamWavetable1SubEnv2:
      case kParamWavetable1SubLFO1:
      case kParamWavetable1SubLFO2:
      case kParamWavetable1SubSeq:
      {
        const int modIdx = paramIdx - kParamWavetable1Sub;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVWavetable1Sub, modIdx, value);
          });
        break;
      }
      case kParamWavetable2:
      {
        value = std::max(1., value);
        UpdateOscillatorWavetable(static_cast<int>(value) - 1, 1);
        break;
      }
      case kParamWavetable2Pitch:
        mParamsToSmooth[kModWavetable2PitchSmoother] = value;
        break;
      case kParamWavetable2PitchEnv1:
      case kParamWavetable2PitchEnv2:
      case kParamWavetable2PitchLFO1:
      case kParamWavetable2PitchLFO2:
      case kParamWavetable2PitchSeq:
      {
        const int modIdx = paramIdx - kParamWavetable2Pitch;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVWavetable2PitchOffset, modIdx, value);
          });
        break;
      }
      case kParamWavetable2Amp:
        mParamsToSmooth[kModWavetable2AmpSmoother] = value;
        break;
      case kParamWavetable2AmpEnv1:
      case kParamWavetable2AmpEnv2:
      case kParamWavetable2AmpLFO1:
      case kParamWavetable2AmpLFO2:
      case kParamWavetable2AmpSeq:
      {
        const int modIdx = paramIdx - kParamWavetable2Amp;
        SendParam([value, modIdx](Voice* voice) {
          voice->UpdateVoiceParam(kVWavetable2Amp, modIdx, value);
          });
          break;
      }
      case kParamWavetable2Pos:
        mParamsToSmooth[kModWavetable2PosSmoother] = value;
        break;
      case kParamWavetable2PosEnv1:
      case kParamWavetable2PosEnv2:
      case kParamWavetable2PosLFO1:
      case kParamWavetable2PosLFO2:
      case kParamWavetable2PosSeq:
      {
        const int modIdx = paramIdx - kParamWavetable2Pos;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVWavetable2Position, modIdx, value);
          });
        break;
      }
      case kParamWavetable2Bend:
        mParamsToSmooth[kModWavetable2BendSmoother] = value;
        break;
      case kParamWavetable2BendEnv1:
      case kParamWavetable2BendEnv2:
      case kParamWavetable2BendLFO1:
      case kParamWavetable2BendLFO2:
      case kParamWavetable2BendSeq:
      {
        const int modIdx = paramIdx - kParamWavetable2Bend;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVWavetable2Bend, modIdx, value);
          });
        break;
      }
      case kParamFilter1Type:
        SendParam([value](Voice* voice) {
          voice->SetFilterType(0, static_cast<int>(value));
          });
        break;
      case kParamFilter1ModeVSF:
      case kParamFilter1ModeMoog:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mFilters.at(0)->SetMode(static_cast<int>(value));
          });
        break;
      case kParamFilter1Cutoff:
        mParamsToSmooth[kModFilter1CutoffSmoother] = value / mSampleRate;
        break;
      case kParamFilter1FF:
      {
        if (mCombOn)
          mParamsToSmooth[kModFilter1CutoffSmoother] = value / 100.;
        break;
      }
      case kParamFilter1CutoffEnv1:
      case kParamFilter1CutoffEnv2:
      case kParamFilter1CutoffAmpEnv:
      case kParamFilter1CutoffLFO1:
      case kParamFilter1CutoffLFO2:
      {
        const int modIdx = paramIdx - kParamFilter1Cutoff;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVFilter1Cutoff, modIdx, value);
          });
        break;
      }
      case kParamFilter1Resonance:
      case kParamFilter1FB:
        mParamsToSmooth[kModFilter1ResonanceSmoother] = value / 100.;
        break;
      case kParamFilter1ResonanceEnv1:
      case kParamFilter1ResonanceEnv2:
      case kParamFilter1ResonanceAmpEnv:
      case kParamFilter1ResonanceLFO1:
      case kParamFilter1ResonanceLFO2:
      {
        const int modIdx = paramIdx - kParamFilter1Resonance;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceParam(kVFilter1Resonance, modIdx, value);
          });
        break;
      }
      case kParamFilter1Drive:
        value /= 100.;
      case kParamFilter1Delay:
        mParamsToSmooth[kModFilter1DriveSmoother] = value;
        break;
      case kParamFilter1DriveEnv1:
      case kParamFilter1DriveEnv2:
      case kParamFilter1DriveAmpEnv:
      case kParamFilter1DriveLFO1:
      case kParamFilter1DriveLFO2:
      {
        const int modIdx = paramIdx - kParamFilter1Drive;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVFilter1Drive, modIdx, value);
          });
        break;
      }
      case kParamFilter2Type:
        SendParam([value](Voice* voice) {
          voice->SetFilterType(1, static_cast<int>(value));
          });
        break;
      case kParamFilter2ModeVSF:
      case kParamFilter2ModeMoog:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mFilters.at(1)->SetMode(static_cast<int>(value));
          });
        break;
      case kParamFilter2Cutoff:
        mParamsToSmooth[kModFilter2CutoffSmoother] = value / mSampleRate;
        break;
      case kParamFilter2FF:
      {
        if (mCombOn)
          mParamsToSmooth[kModFilter2CutoffSmoother] = value / 100.;
        break;
      }
      case kParamFilter2CutoffEnv1:
      case kParamFilter2CutoffEnv2:
      case kParamFilter2CutoffAmpEnv:
      case kParamFilter2CutoffLFO1:
      case kParamFilter2CutoffLFO2:
      {
        const int modIdx = paramIdx - kParamFilter2Cutoff;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVFilter2Cutoff, modIdx, value);
          });
        break;
      }
      case kParamFilter2Resonance:
      case kParamFilter2FB:
        mParamsToSmooth[kModFilter2ResonanceSmoother] = value / 100.;
        break;
      case kParamFilter2ResonanceEnv1:
      case kParamFilter2ResonanceEnv2:
      case kParamFilter2ResonanceAmpEnv:
      case kParamFilter2ResonanceLFO1:
      case kParamFilter2ResonanceLFO2:
      {
        const int modIdx = paramIdx - kParamFilter2Resonance;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVFilter2Resonance, modIdx, value / 100.);
          });
        break;
      }
      case kParamFilter2Drive:
        value /= 100.;
      case kParamFilter2Delay:
        mParamsToSmooth[kModFilter2DriveSmoother] = value;
        break;
      case kParamFilter2DriveEnv1:
      case kParamFilter2DriveEnv2:
      case kParamFilter2DriveAmpEnv:
      case kParamFilter2DriveLFO1:
      case kParamFilter2DriveLFO2:
      {
        const int modIdx = paramIdx - kParamFilter2Drive;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVFilter2Drive, modIdx, value / 100.);
          });
        break;
      }
      case kParamOsc1PM:
      {
        const bool phaseModOn = value > 0.5;
        SendParam([phaseModOn](Voice* voice) {
          voice->mOsc1.SetPhaseModulation(phaseModOn);
          });
        break;
      }
      case kParamOsc2PM:
      {
        const bool phaseModOn = value > 0.5;
        SendParam([phaseModOn](Voice* voice) {
          voice->mOsc2.SetPhaseModulation(phaseModOn);
          });
        break;
      }
      case kParamPhaseModFreq:
        mParamsToSmooth[kModPhaseModFreqSmoother] = value;
        break;
      case kParamPhaseModFreqEnv1:
      case kParamPhaseModFreqEnv2:
      case kParamPhaseModFreqAmpEnv:
      case kParamPhaseModFreqLFO1:
      case kParamPhaseModFreqLFO2:
      case kParamPhaseModFreqSeq:
      {
        const int modIdx = paramIdx - kParamPhaseModFreq;
        SendParam([value, modIdx](Voice* voice) {
          voice->UpdateVoiceParam(kVPhaseModFreq, modIdx, value);
          });
        break;
      }
      case kParamPhaseModAmount:
        mParamsToSmooth[kModPhaseModAmtSmoother] = value / 100.;
        break;
      case kParamPhaseModAmountEnv1:
      case kParamPhaseModAmountEnv2:
      case kParamPhaseModAmountAmpEnv:
      case kParamPhaseModAmountLFO1:
      case kParamPhaseModAmountLFO2:
      case kParamPhaseModAmountSeq:
      {
        const int modIdx = paramIdx - kParamPhaseModAmount;
        SendParam([value, modIdx](Voice* voice) {
          voice->UpdateVoiceParam(kVPhaseModAmt, modIdx, value);
          });
        break;
      }
      case kParamOsc1RM:
      {
        const bool ringModOn = value > 0.5;
        SendParam([ringModOn](Voice* voice) {
          voice->mOsc1.SetRingModulation(ringModOn);
          });
        break;
      }
      case kParamOsc2RM:
      {
        const bool ringModOn = value > 0.5;
        SendParam([ringModOn](Voice* voice) {
          voice->mOsc2.SetRingModulation(ringModOn);
          });
        break;
      }
      case kParamRingModFreq:
        mParamsToSmooth[kModRingModFreqSmoother] = value;
        break;
      case kParamRingModFreqEnv1:
      case kParamRingModFreqEnv2:
      case kParamRingModFreqAmpEnv:
      case kParamRingModFreqLFO1:
      case kParamRingModFreqLFO2:
      case kParamRingModFreqSeq:
      {
        const int modIdx = paramIdx - kParamRingModFreq;
        SendParam([value, modIdx](Voice* voice) {
          voice->UpdateVoiceParam(kVRingModFreq, modIdx, value);
          });
        break;
      }
      case kParamRingModAmount:
        mParamsToSmooth[kModRingModAmtSmoother] = value / 100.;
        break;
      case kParamRingModAmountEnv1:
      case kParamRingModAmountEnv2:
      case kParamRingModAmountAmpEnv:
      case kParamRingModAmountLFO1:
      case kParamRingModAmountLFO2:
      case kParamRingModAmountSeq:
      {
        const int modIdx = paramIdx - kParamRingModAmount;
        SendParam([value, modIdx](Voice* voice) {
          voice->UpdateVoiceParam(kVRingModAmt, modIdx, value);
          });
        break;
      }
      default:
        break;
    }
  }
  
public:
  MidiSynth mSynth { VoiceAllocator::kPolyModePoly, MidiSynth::kDefaultBlockSize };
  WDL_TypedBuf<T> mModulationsData; // Sample data for global modulations (e.g. smoothed sustain)
  WDL_PtrList<T> mModulations; // Ptrlist for global modulations
  LogParamSmooth<T, kNumModulations> mParamSmoother;
  sample mParamsToSmooth[kNumModulations];
  LFO<T> mLFO;
  std::vector<std::string> mWavetables{ ELEMENT_NAMES };
  std::vector<Voice*> mSynthVoices;

  double mSampleRate{ 44100. };
  static inline bool mCombOn{ false };

  static inline bool tableLoading[2]{ true, true };

  // Non-modulatable parameters
  double mLoadedWavetables[2]{ 1., 2. }; // Integer indices of current wavetables
  double mSeqSteps[kNumSeqSteps]{}; // Value of each step in the sequencer
};