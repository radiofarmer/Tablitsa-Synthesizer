#pragma once

#include "PeriodicTable.h"
#include "Effects.h"
#include "TablitsaOscillators.h"
#include "Modulation.h"


#include "IPlugConstants.h"
#include "Oscillator.h"
#include "MidiSynth.h"
#include "ADSREnvelope.h"
#include "Smoothers.h"
#include "LFO.h"

#ifdef VECTOR
 #define FRAME_INTERVAL OUTPUT_SIZE
#else
  #define FRAME_INTERVAL 1
#endif

constexpr double kMaxEnvTimeScalar = 0.5;

using namespace iplug;

/*
Global Modulations (smoothers): These values are computed once per sample and sent to all voices 
*/
enum EModulations
{
  kModGainSmoother = 0,
  kModPanSmoother,
  kModEnv1SustainSmoother,
  kModEnv2SustainSmoother,
  kModAmpEnvSustainSmoother,
  kModWavetable1PitchSmoother,
  kModWavetable1PosSmoother,
  kModWavetable1BendSmoother,
  kModWavetable1FormantSmoother,
  kModWavetable1AmpSmoother,
  kModWavetable2PitchSmoother,
  kModWavetable2PosSmoother,
  kModWavetable2BendSmoother,
  kModWavetable2FormantSmoother,
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
  kVWavetable1Formant,
  kVWavetable1Amp,
  kVWavetable2PitchOffset,
  kVWavetable2Position,
  kVWavetable2Bend,
  kVWavetable2Formant,
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

enum EVoiceMetaModParams
{
  kVEnv1Sustain,
  kVEnv2Sustain,
  kVAmpEnvSustain,
  kVLFO1RateHz,
  kVLFO1Amp,
  kVLFO2RateHz,
  kVLFO2Amp,
  kVSequencerRateHz,
  kVSequencerAmp,
  kNumVoiceMetaModulations
};

// See Modulators.h for the enumeration of the number of modulators

#define UNISON_CHORD_LIST "None", "8va", "M7", "D7", "D7 6/5", "D7 4/3", "m7", "ø7", "dim7"

enum EUnisonChords
{
  kNoChord = -1,
  kOctaves = 0,
  kMaj7,
  kDom7,
  kDom7FirstInv,
  kDom7SecondInv,
  kMin7,
  kHalfDim7,
  kDim7
};

struct ChannelParam
{
  double L{ 0. };
  double R{ 0. };
  double LR[2]{ 0. };

  void Set(double l, double r)
  {
    LR[0] = L = l;
    LR[1] = R = r;
  }

  double& operator [](int ch)
  {
    return LR[ch];
  }
};

struct UnisonVoiceManager
{
  int mMaxVoices;

  // Detuning
  double mMaxDetune; // octaves
  int mNVoices{ 1 };
  int mDetuneVoiceIdx{ 0 };
  int mChord{ EUnisonChords::kNoChord };
  double* mDetuneBuf;

  // Panning
  int mPanVoiceIdx{ 0 };
  double mMaxPan; // [-1., 1.]
  ChannelParam* mPanBuf;

  const double mUnisonInvervals[8][5]{
    {0., 1., 0., 2., 0.}, // Octaves
    {0., 7. / 12, 4. / 12, 11. / 12, 1.}, // Major 7
    {0., 7. / 12, 4. / 12, 10. / 12, 1.}, // Dominant 7
    {0., 7. / 12, 4. / 12, -2. / 12, 1.}, // Dominant 7 1st Inversion
    {0., -5. / 12, -2. / 12, 4. / 12, 1.}, // Dominant 7 2nd Inversion
    {0., 7. / 12, 3. / 12, 10. / 12, 1.}, // Minor 7
    {0., 3. / 12., 6. / 12, 10. / 12, 1.}, // Half-Diminished 7
    {0., 3. / 12., 6. / 12, 9. / 12, 1.} // Full-Diminished 7
  };

  UnisonVoiceManager(int maxVoices, double maxDetuneSemitones=1.) : mMaxVoices(maxVoices), mMaxDetune(maxDetuneSemitones)
  {
    mDetuneBuf = new double[mMaxVoices] {0.};
    mPanBuf = new ChannelParam[mMaxVoices]{};
    mPanBuf[0].Set(1., 1.);
  }

  void SetNVoices(int nVoices)
  {
    mNVoices = std::min(nVoices, mMaxVoices);
    RefillBuffers();
  }

  void SetMaxDetune(double maxDetuneSemitones, bool refillBuffer=true)
  {
    mMaxDetune = maxDetuneSemitones / 12;
    if (refillBuffer)
      RefillDetuneBuffer();
  }

  void SetMaxPan(double maxPanDegrees, bool refillBuffer=true)
  {
    mMaxPan = maxPanDegrees / 180.;
    if (refillBuffer)
      RefillPanBuffer();
  }

  void Reset()
  {
    mDetuneVoiceIdx = 0;
    mPanVoiceIdx = 0;
  }

  void SetChord(int chord)
  {
    mChord = chord;
  }

  void RefillBuffers()
  {
    RefillDetuneBuffer();
    RefillPanBuffer();
  }

  inline void RefillDetuneBuffer()
  {
    int nVoices = std::max(2, mNVoices);

    // Pitch detuning
    if (mChord == EUnisonChords::kNoChord)
    {
      for (auto i{ 0 }; i < mNVoices; ++i)
        mDetuneBuf[i] = mMaxDetune * static_cast<double>(-i % mNVoices) / (static_cast<double>(nVoices) - 1.) * -1;
    }
    else
    {
      for (auto i{ 0 }; i < mNVoices; ++i)
        mDetuneBuf[i] = mUnisonInvervals[mChord][i % 5] + mMaxDetune * (static_cast<double>(std::rand() % 100) / 50 - 1.);
    }
  }

  inline void RefillPanBuffer()
  {
    constexpr double sqrt2{ 1.45 };
    // Panning (in progress)
    double totalPan = std::abs(mMaxPan);
    if (mNVoices > 1)
    {
      double pan = totalPan;
      mPanBuf[0].Set(
        (1. + std::copysign(pan, mMaxPan)) / sqrt2,
        (1. - std::copysign(pan, mMaxPan)) / sqrt2
      );
    }
    else
      mPanBuf[0].Set(1., 1.);
    for (auto i{ 1 }; i < mNVoices; ++i)
    {
      double pan;
      if (i % 2)
        pan = totalPan / i; // Odd-numbered voices: Pan by the same magnitude as the last voice, but in the opposite direction
      else
        pan = totalPan / (i + 1); // Even-numbered voices: Pan within the maximum pan range, to ever-smaller extents
      mPanBuf[i].Set(
        (1. - std::copysign(pan, mMaxPan)) / sqrt2,
        (1. + std::copysign(pan, mMaxPan)) / sqrt2
      );
    }
  }

  double DetuneNext()
  {
    int next = mDetuneVoiceIdx++ % mNVoices;
    mDetuneVoiceIdx %= mNVoices;
    return mDetuneBuf[next];
  }

  void SetNextPan(double* channelPan)
  {
    int next = mPanVoiceIdx++ % mNVoices;
    mPanVoiceIdx %= mNVoices;
    channelPan[0] = mPanBuf[next][0];
    channelPan[1] = mPanBuf[next][1];
  }

  ~UnisonVoiceManager()
  {
    delete[] mDetuneBuf;
  }
};

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

    Voice(TablitsaDSP<T>* master, int id = 0)
      : mMaster(master), mID(id),
      mAmpEnv("gain", [&]() { mOsc1.Reset(); }),
      mEnv1("env1", [&]() { mOsc1.Reset(); }),
      mLFO1(&GetMaster()->mGlobalMetronome),
      mLFO2(&GetMaster()->mGlobalMetronome),
      mSequencer(&GetMaster()->mGlobalMetronome, GetMaster()->mSeqSteps) // capture ok on RT thread?
    {
//      DBGMSG("new Voice: %i control inputs.\n", static_cast<int>(mInputs.size()));

      // Meta-Modulation Parameters
      std::vector<ParameterModulator<6, 3>*> env1Mods{ &mVoiceMetaModParams[kVEnv1Sustain] };
      std::vector<ParameterModulator<6, 3>*> env2Mods{ &mVoiceMetaModParams[kVEnv2Sustain] };
      std::vector<ParameterModulator<6, 3>*> ampEnvMods{ &mVoiceMetaModParams[kVAmpEnvSustain] };
      std::vector<ParameterModulator<6, 3>*> lfo1Mods{ &mVoiceMetaModParams[kVLFO1RateHz], &mVoiceMetaModParams[kVLFO1Amp] };
      std::vector<ParameterModulator<6, 3>*> lfo2Mods{ &mVoiceMetaModParams[kVLFO2RateHz], &mVoiceMetaModParams[kVLFO2Amp] };
      std::vector<ParameterModulator<6, 3>*> seqMods{ &mVoiceMetaModParams[kVSequencerRateHz], &mVoiceMetaModParams[kVSequencerAmp] };

      mModulators.AddModulator(&mEnv1, env1Mods);
      mModulators.AddModulator(&mEnv2, env2Mods);
      mModulators.AddModulator(&mAmpEnv, ampEnvMods);
      mModulators.AddModulator(&mLFO1, lfo1Mods);
      mModulators.AddModulator(&mLFO2, lfo2Mods);
      mModulators.AddModulator(&mSequencer, seqMods);
      mAmpEnv.Kill(true); // Force amplitude envelopes to start in the "Idle" stage

      // Fill the envelope queues for legato mode with null pointers
      GetMaster()->AmpEnvQueue.push_back(nullptr);
      GetMaster()->Env1Queue.push_back(nullptr);
      GetMaster()->Env2Queue.push_back(nullptr);
    }

    bool GetBusy() const override
    {
      return mAmpEnv.GetBusy();
    }

    UnisonVoiceManager& GetDetuner()
    {
      return GetMaster()->mDetuner;
    }

    void ResetUnisonParams()
    {
      mDetune = GetDetuner().DetuneNext();
      GetDetuner().SetNextPan(mPan);
    }

    void Trigger(double level, bool isRetrigger) override
    {
      mOsc1.Reset();
      mOsc2.Reset();

      mVelocity = level; // TODO: Handling of different velocity settings (i.e. which envelopes are affected by velocity)
      mTriggerRand = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);

      ResetUnisonParams();

      for (auto f : mFilters)
      {
        f->Reset();
      }

      // Reset LFOs and Sequencer
      if (mLFO1Restart)
      {
        mModulators.ReplaceModulator(&mLFO1, 0);
        mLFO1.Reset();
      }
      else
      {
        mModulators.ReplaceModulator(dynamic_cast<FastLFO<T>*>(&(GetMaster()->mGlobalLFO1)), 0);
      }

      if (mLFO2Restart)
      {
        mModulators.ReplaceModulator(&mLFO2, 1);
        mLFO2.Reset();
      }
      else
      {
        mModulators.ReplaceModulator(dynamic_cast<FastLFO<T>*>(&GetMaster()->mGlobalLFO2), 1);
      }

      if (mSequencerRestart)
      {
        mModulators.ReplaceModulator(dynamic_cast<Sequencer<T>*>(&mSequencer), 0);
        mSequencer.Reset();
        // Update sequencer display with this voice's phase
        GetMaster()->mActiveSequencer = &mSequencer;
      }
      else
      {
        mModulators.ReplaceModulator(&GetMaster()->mGlobalSequencer, 0);
        // Update sequencer display with this voice's phase
        GetMaster()->mActiveSequencer = &GetMaster()->mGlobalSequencer;
      }


      if (isRetrigger)
      {
        mAmpEnv.Retrigger(level);
        mEnv1.Retrigger(1., 1. - mEnv1VelocityMod);
        mEnv2.Retrigger(1.);
      }
      else if (!mLegato)
      {
        double velSubtr = 1. - level;
        mAmpEnv.Start(1. - velSubtr * mAmpEnvVelocityMod, 1. - mAmpEnvVelocityMod * kMaxEnvTimeScalar * level);
        mEnv1.Start(1. - velSubtr * mEnv1VelocityMod, 1. - mEnv1VelocityMod * kMaxEnvTimeScalar * level);
        mEnv2.Start(1. - velSubtr * mEnv2VelocityMod, 1. - mEnv2VelocityMod * kMaxEnvTimeScalar * level);
      }
      else
      {
        // Check for currently active envelopes
        Envelope<T>* masterAmpEnv = nullptr;
        Envelope<T>* masterEnv1 = nullptr;
        Envelope<T>* masterEnv2 = nullptr;
        for (auto i{ 0 }; i < GetMaster()->AmpEnvQueue.size(); ++i)
        {
          if (GetMaster()->AmpEnvQueue[i])
          {
            masterAmpEnv = GetMaster()->AmpEnvQueue[i];
            masterEnv1 = GetMaster()->Env1Queue[i];
            masterEnv2 = GetMaster()->Env2Queue[i];
          }
        }
        // If active envelopes were found, sync this voice's envelopes to them
        if (masterAmpEnv)
        {
          double velSubtr = 1. - level;
          mAmpEnv.StartAt(1. - velSubtr * mAmpEnvVelocityMod,
            masterAmpEnv->GetValue(), masterAmpEnv->GetPrevResult(), masterAmpEnv->GetStage(),
            1. - mAmpEnvVelocityMod * kMaxEnvTimeScalar * level);
          mEnv1.StartAt(1. - velSubtr * mEnv1VelocityMod,
            masterEnv1->GetValue(), masterEnv1->GetPrevResult(), masterEnv1->GetStage(),
            1. - mEnv1VelocityMod * kMaxEnvTimeScalar * level);
          mEnv2.StartAt(1. - velSubtr * mEnv2VelocityMod,
            masterEnv2->GetValue(), masterEnv2->GetPrevResult(), masterEnv2->GetStage(),
            1. - mEnv2VelocityMod * kMaxEnvTimeScalar * level);

        }
        else
        {
          double velSubtr = 1. - level;
          mAmpEnv.Start(1. - velSubtr * mAmpEnvVelocityMod, 1. - mAmpEnvVelocityMod * kMaxEnvTimeScalar * level);
          mEnv1.Start(1. - velSubtr * mEnv1VelocityMod, 1. - mEnv1VelocityMod * kMaxEnvTimeScalar * level);
          mEnv2.Start(1. - velSubtr * mEnv2VelocityMod, 1. - mEnv2VelocityMod * kMaxEnvTimeScalar * level);
        }
        // Sync the master envelopes to this voice's envelopes
        GetMaster()->Env1Queue[mID] = &mEnv1;
        GetMaster()->Env2Queue[mID] = &mEnv2;
        GetMaster()->AmpEnvQueue[mID] = &mAmpEnv;
      }
    }
    
    void Release() override
    {
      mAmpEnv.Release();
      mEnv1.Release();
      mEnv2.Release();
      // Remove this voice's envelopes from the envelope queue
      GetMaster()->AmpEnvQueue[mID] = nullptr;
      GetMaster()->Env1Queue[mID] = nullptr;
      GetMaster()->Env2Queue[mID] = nullptr;
    }

    void ProcessSamplesAccumulating(T** inputs, T** outputs, int nInputs, int nOutputs, int startIdx, int nFrames) override
    {
      constexpr double ktShelfIncr{ 0.001 };
      // inputs to the synthesizer can just fetch a value every block, like this:
//      double gate = mInputs[kVoiceControlGate].endValue;
      double pitch = mInputs[kVoiceControlPitch].endValue + mDetune; // pitch = (MidiKey - 69) / 12
//      double pitchBend = mInputs[kVoiceControlPitchBend].endValue;
      // or write the entire control ramp to a buffer, like this, to get sample-accurate ramps:
//      mInputs[kVoiceControlTimbre].Write(mTimbreBuffer.Get(), startIdx, nFrames);

      // Set the static (note-constant) modulator values
      T keytrack = (pitch * 12. + 69.) / 128.;
      T staticMods[]{ mVelocity, keytrack, mTriggerRand };
      mVoiceModParams.SetStaticModulation(staticMods);
      // Write ramps for modulators
      mVoiceMetaModParams[kVEnv1Sustain].SetInitialValue(inputs[kModEnv1SustainSmoother][0]);
      mVoiceMetaModParams[kVEnv2Sustain].SetInitialValue(inputs[kModEnv2SustainSmoother][0]);
      mVoiceMetaModParams[kVAmpEnvSustain].SetInitialValue(inputs[kModAmpEnvSustainSmoother][0]);
      mModulators.MetaProcessBlock_Fast(&(inputs[kModEnv1SustainSmoother]), nFrames); // Modulate the modulators
      // mModulators.ProcessBlock(&(inputs[kModEnv1SustainSmoother]), nFrames);
      // Apply modulation ramps to all modulated parameters
#if 0
      mVoiceModParams.ProcessBlockVec4d(&inputs[kModWavetable1PitchSmoother], mModulators.GetList(), mVModulations.GetList(), nFrames);
#else
      mVoiceModParams.ProcessBlock(&inputs[kModWavetable1PitchSmoother], mModulators.GetList(), mVModulations.GetList(), nFrames);
#endif

      const double phaseModFreqFact = pow(2., mVModulations.GetList()[kVPhaseModFreq][0] / 12.);
      const double ringModFreqFact = pow(2., mVModulations.GetList()[kVRingModFreq][0] / 12.);

      const int osc1Filter = mFilterRouting[0];
      const int osc2Filter = mFilterRouting[1];

      // make sound output for each output channel
      for(auto i = startIdx; i < startIdx + nFrames; i += FRAME_INTERVAL)
      {
        int bufferIdx = i - startIdx;
//        float noise = mTimbreBuffer.Get()[i] * Rand();
        double ampEnvVal{ mModulators.GetList()[2][bufferIdx] }; // Calculated for easy access

        // Oscillator Parameters
        double osc1Freq = 440. * pow(2., pitch + mVModulations.GetList()[kVWavetable1PitchOffset][bufferIdx] / 12.);
        mOsc1.SetWtPosition(1 - mVModulations.GetList()[kVWavetable1Position][bufferIdx]); // Wavetable 1 Position
        mOsc1.SetWtBend(mVModulations.GetList()[kVWavetable1Bend][bufferIdx]); // Wavetable 1 Bend
        mOsc1Sat.SetLevel(mVModulations.GetList()[kVWavetable1Formant][bufferIdx], keytrack * ktShelfIncr);
        mOsc1.SetPhaseModulation(mVModulations.GetList()[kVPhaseModAmt][bufferIdx], osc1Freq * phaseModFreqFact);
        mOsc1.SetRingModulation(mVModulations.GetList()[kVRingModAmt][bufferIdx], osc1Freq * ringModFreqFact);

        double osc2Freq = 440. * pow(2., pitch + mVModulations.GetList()[kVWavetable2PitchOffset][bufferIdx] / 12.);
        mOsc2.SetWtPosition(1 - mVModulations.GetList()[kVWavetable2Position][bufferIdx]); // Wavetable 2 Position
        mOsc2.SetWtBend(mVModulations.GetList()[kVWavetable2Bend][bufferIdx]); // Wavetable 2 Bend
        mOsc2Sat.SetLevel(mVModulations.GetList()[kVWavetable2Formant][bufferIdx], keytrack * ktShelfIncr);
        mOsc2.SetPhaseModulation(mVModulations.GetList()[kVPhaseModAmt][bufferIdx], osc2Freq * phaseModFreqFact);
        mOsc2.SetRingModulation(mVModulations.GetList()[kVRingModAmt][bufferIdx], osc2Freq * ringModFreqFact);
        
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
         // Saturation
         osc1Output[j] = mOsc1Sat.Process(osc1Output[j]); 
         osc2Output[j] = mOsc2Sat.Process(osc2Output[j]);
         // Filters
         double osc1FilterOutput = mFilters.at(osc1Filter)->Process(osc1Output[j]);
         double osc2FilterOutput = mFilters.at(osc2Filter)->Process(osc2Output[j]);
         double output_summed = osc1FilterOutput * mVModulations.GetList()[kVWavetable1Amp][bufferIdx] + osc2FilterOutput * mVModulations.GetList()[kVWavetable2Amp][bufferIdx];
         double output_scaled = output_summed * ampEnvVal * mGain;
         outputs[0][i + j] += output_scaled * mPan[0];
         outputs[1][i + j] += output_scaled * mPan[1];
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

      mFilters[0]->SetSampleRate(sampleRate);
      mFilters[1]->SetSampleRate(sampleRate);

      mOsc1Sat.SetSampleRate(sampleRate);
      mOsc2Sat.SetSampleRate(sampleRate);
      
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

    void SetFilterType(int filter, int filterType)
    {
      // Indices of cutoff/res/drive for the given filter
      switch (filterType) {
      case kVSF:
      {
        mVoiceModParams[filter ? kVFilter2Cutoff : kVFilter1Cutoff].SetMinMax(0.001, 0.49);
        mVoiceModParams[filter ? kVFilter2Resonance : kVFilter1Resonance].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Drive : kVFilter1Drive].SetMinMax(0., 1.);
        mFilters.at(filter) = new SVF2<T>(mOsc1.GetSampleRate());
        break;
      }
      case kMoog:
      {
        mVoiceModParams[filter ? kVFilter2Cutoff : kVFilter1Cutoff].SetMinMax(0.001, 0.49);
        mVoiceModParams[filter ? kVFilter2Resonance : kVFilter1Resonance].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Drive : kVFilter1Drive].SetMinMax(0., 1.);
        mFilters.at(filter) = new MoogLadder<T>(mOsc1.GetSampleRate());
        break;
      }
      case kComb:
      {
        mVoiceModParams[filter ? kVFilter2Cutoff : kVFilter1Cutoff].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Resonance : kVFilter1Resonance].SetMinMax(0., 1.);
        mVoiceModParams[filter ? kVFilter2Drive : kVFilter1Drive].SetMinMax(0., 1.);
        mFilters.at(filter) = new CombFilter<T>(mOsc1.GetSampleRate());
        break;
      }
      default:
        mFilters.at(filter) = new NullFilter<T>(mOsc1.GetSampleRate());
        break;

      }
    }

    /* Route the given oscillator through the given filter */
    void UpdateFilterSource(int filterIdx, int oscIdx)
    {
      mFilterRouting[oscIdx] = filterIdx;
    }

    /* Update polyphonic modulation depths */
    void UpdateVoiceParam(int voiceParam, int modIdx, double value)
    {
      mVoiceModParams[voiceParam].SetValue(modIdx - 1, value); // Eventually adjust the value of modIdx in the SetParam switch statement rather than here
    }

    /* Update meta-modulation depths */
    void UpdateVoiceModulatorParam(int paramIdx, int modIdx, double value)
    {
      if (modIdx == 0)
        mVoiceMetaModParams[paramIdx].SetInitialValue(value);
      else
        mVoiceMetaModParams[paramIdx].SetValue(modIdx - 1, value);
    }

    TablitsaDSP<T>* GetMaster() const
    {
      return mMaster;
    }

  public:
    TablitsaDSP<T>* mMaster;

    WavetableOscillator<T> mOsc1{ 0, "Hydrogen" };
    WavetableOscillator<T> mOsc2{ 1, "Helium" };
    SaturationEQ<T> mOsc1Sat;
    SaturationEQ<T> mOsc2Sat;

    // Static Modulators
    T mKey{ 69. };
    T mVelocity{ 1. };
    T mTriggerRand{ 0.5 };

    // Dynamic Modulators
    Envelope<T> mEnv1;
    Envelope<T> mEnv2;
    Envelope<T> mAmpEnv;
    FastLFO<T> mLFO1;
    FastLFO<T> mLFO2;
    Sequencer<T, kNumSeqSteps> mSequencer;
    ModulatorList<T, Envelope<T>, FastLFO<T>, Sequencer<T>, 6, 3> mModulators;

    // Status parameters
    bool mLFO1Restart{ false };
    bool mLFO2Restart{ false };
    bool mSequencerRestart{ false };
    bool mLegato{ false };

    // Static Modulators
    double mEnv1VelocityMod{ 0. };
    double mEnv2VelocityMod{ 0. };
    double mAmpEnvVelocityMod{ 1. };

    // Filters
    std::vector<Filter<T>*> mFilters{ new NullFilter<T>(), new NullFilter<T>() };
    int mFilterRouting[2]{ 0, 1 };

    WDL_PtrList<T> mVModulations; // Pointers to modulator buffers
    WDL_TypedBuf<T> mVModulationsData; // Modulator buffer sample data

  private:
//    WDL_TypedBuf<float> mTimbreBuffer;
    ModulatedParameterList<T, kNumVoiceModulations> mVoiceModParams{
      new ParameterModulator<>(-24., 24., "Wt1 Pitch Offset"),
      new ParameterModulator<>(0., 1., "Wt1 Position"),
      new ParameterModulator<>(-1., 1., "Wt1 Bend"),
      new ParameterModulator<>(0., 1., "Wt1 Formant"),
      new ParameterModulator<>(0., 1., "Wt1 Amp"),
      new ParameterModulator<>(-24., 24., "Wt1 Pitch Offset"),
      new ParameterModulator<>(0., 1., "Wt2 Position"),
      new ParameterModulator<>(-1., 1., "Wt2 Bend"),
      new ParameterModulator<>(0., 1., "Wt2 Formant"),
      new ParameterModulator<>(0., 1., "Wt2 Amp"),
      new ParameterModulator<>(0.001, 0.5, "Flt1 Cutoff", true),
      new ParameterModulator<>(0., 1., "Flt1 Resonance"),
      new ParameterModulator<>(0., 1., "Flt1 Drive"),
      new ParameterModulator<>(0.001, 0.5, "Flt2 Cutoff", true),
      new ParameterModulator<>(0., 1., "Flt2 Resonance"),
      new ParameterModulator<>(0., 1., "Flt2 Drive"),
      new ParameterModulator<>(-24., 24., "Phase Mod Freq"),
      new ParameterModulator<>(0., 1., "Phase Mod Depth"),
      new ParameterModulator<>(-24., 24., "Ring Mod Freq"),
      new ParameterModulator<>(0., 1., "Ring Mod Depth") };

    // Modulator parameters that can themselves be modulated
    ModulatedParameterList<T, kNumVoiceMetaModulations> mVoiceMetaModParams{
      new ParameterModulator<>(0., 1., "Env1 Sustain"),
      new ParameterModulator<>(0., 1., "Env2 Sustain"),
      new ParameterModulator<>(0., 1., "AmpEnv Sustain"),
      new ParameterModulator<>(0.01, 40., "LFO1 Rate Hz", true),
      new ParameterModulator<>(-1., 1., "LFO1 Amp"),
      new ParameterModulator<>(0.01, 40., "LFO2 Rate Hz", true),
      new ParameterModulator<>(-1., 1., "LFO2 Amp"),
      new ParameterModulator<>(0.01, 40., "Sequencer Rate Hz", true),
      new ParameterModulator<>(0., 1., "Sequencer Amp"),

    };

    // Unison parameters
    double mDetune{ 0. };
    double mPan[2]{ 1., 1. };

    // Sample and Beat data
    double mTempo{ 120. };
    bool mTransportIsRunning{ false };
    double mQNPos{ 0. };

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
    assert((int)kNumVoiceModulations == (int)(kNumModulations - kModWavetable1PitchSmoother));
    for (auto i = 0; i < nVoices; i++)
    {
      // add a voice to Zone 0.
      mSynthVoices.push_back(new Voice(this, i));
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
    // Process global modulators
    mGlobalLFO1.FillBuffer(nFrames);
    mGlobalLFO2.FillBuffer(nFrames);
    mGlobalSequencer.FillBuffer(nFrames);

    // Process voices
    mParamSmoother.ProcessBlock(mParamsToSmooth, mModulations.GetList(), nFrames); // Populate modulations list (to be sent to mSynth as inputs)
    SetTempoAndBeat(qnPos, transportIsRunning, tempo);
    mSynth.ProcessBlock(mModulations.GetList(), outputs, 0, nOutputs, nFrames);

    for(int s=0; s < nFrames;s++)
    {
      const T smoothedGain = mModulations.GetList()[kModGainSmoother][s];

      // Master effects processing
      T delay[2]{ outputs[0][s], outputs[1][s] };
      mDelayEffect.ProcessStereo(delay);
      //const T delay[2]{ mDelayEffect.Process(outputs[0][s]), mDelayEffect.Process(outputs[0][s]) };

      outputs[0][s] += delay[0];
      outputs[1][s] += delay[1];

      outputs[0][s] *= smoothedGain * (2 - mModulations.GetList()[kModPanSmoother][s]);
      outputs[1][s] *= smoothedGain * mModulations.GetList()[kModPanSmoother][s];
    }
  }

  void Reset(double sampleRate, int blockSize)
  {
    mSampleRate = sampleRate;
    mSynth.SetSampleRateAndBlockSize(sampleRate, blockSize);
    ResetAllVoices();
    mSynth.ForEachVoice([sampleRate](SynthVoice& voice) {
      for (auto f : dynamic_cast<TablitsaDSP::Voice&>(voice).mFilters)
        f->SetSampleRate(sampleRate);
      });

    // Global modulators
    mGlobalLFO1.SetSampleRate(sampleRate);
    mGlobalLFO2.SetSampleRate(sampleRate);
    mGlobalSequencer.SetSampleRate(sampleRate);
    mGlobalLFO1.Resize(blockSize);
    mGlobalLFO2.Resize(blockSize);
    mGlobalSequencer.Resize(blockSize);

    // Param Smoother list
    mModulationsData.Resize(blockSize * kNumModulations);
    mModulations.Empty();

    // Effects
    mDelayEffect.SetSampleRate(sampleRate);
    
    for(auto i = 0; i < kNumModulations; i++)
    {
      mModulations.Add(mModulationsData.Get() + static_cast<size_t>(blockSize) * i);
    }
  }

  void ResetAllVoices()
  {
    mSynth.Reset(); // Note: this sets all envelopes to the release stage, meaning all voices are technically active
    mSynth.ForEachVoice([](SynthVoice& voice) {
      dynamic_cast<TablitsaDSP::Voice&>(voice).mAmpEnv.Kill(true); // Set all amp envelopes to idle
      });
  }

  void ProcessMidiMsg(const IMidiMsg& msg, float detune=0.f)
  {
    if (!mConstantGlideTime && msg.StatusMsg() == IMidiMsg::EStatusMsg::kNoteOn)
    {
      mSynth.SetNoteGlideTime(std::abs(static_cast<double>(msg.NoteNumber()) - mLastNoteOn) / mGlideRateScalar);
      mLastNoteOn = static_cast<double>(msg.NoteNumber());
    }
    mSynth.AddMidiMsgToQueue(msg);
  }

  inline void SetTempoAndBeat(double qnPos, bool transportIsRunning, double tempo)
  {
    mGlobalMetronome.Set(qnPos, tempo, transportIsRunning);
  }

  void UpdateOscillatorWavetable(int wtIdx, int oscIdx)
  {
    ResetAllVoices();
    WtFile wtFile{ mWavetables.at(wtIdx) };

    if (oscIdx == 0)
    {
      mSynth.ForEachVoice([&wtFile, oscIdx](SynthVoice& voice) {
        dynamic_cast<Voice&>(voice).mOsc1.LoadNewTable(wtFile, oscIdx);
        });
      SendParam([this, oscIdx, &wtFile](Voice* voice) {
        voice->mOsc1.SetWavetable(oscIdx);
        voice->mOsc1.ReloadLUT();
        voice->mOsc1.NotifyLoaded();
        });
    }
    else
    {
      mSynth.ForEachVoice([&wtFile, oscIdx](SynthVoice& voice) {
        dynamic_cast<Voice&>(voice).mOsc2.LoadNewTable(wtFile, oscIdx);
        });
      SendParam([this, oscIdx, &wtFile](Voice* voice) {
        voice->mOsc2.SetWavetable(oscIdx);
        voice->mOsc2.ReloadLUT();
        voice->mOsc2.NotifyLoaded();
        });
    }
  }

  void ResetDetune()
  {
    mDetuner.RefillBuffers();
  }

  int GetSequencerStep()
  {
    if (mActiveSequencer)
      return mActiveSequencer->GetCurrentStep();
    else
      return 0;
  }

  void UpdateFilterSource(int filterIdx, int oscIdx)
  {
    mSynth.ForEachVoice([filterIdx, oscIdx](SynthVoice& voice) {
      dynamic_cast<TablitsaDSP<T>::Voice&>(voice).UpdateFilterSource(filterIdx, oscIdx);
      });
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
      case kParamNoteGlideRate:
        mGlideRateScalar = value;
        break;
      case kParamPortamentoMode:
        mConstantGlideTime = value > 0.5;
        break;
      case kParamMonophonic:
      {
        mMono = value > 0.5;
        ResetAllVoices();
        mSynth.SetPolyMode(value > 0.5 ? VoiceAllocator::EPolyMode::kPolyModeMono : VoiceAllocator::EPolyMode::kPolyModePoly);
        break;
      }
      case kParamGain:
      {
#ifdef VST3_API
        mParamsToSmooth[kModGainSmoother] = (T)value / 100.;
#else
        mParamsToSmooth[kModGainSmoother] = std::pow(10., value / 20.);
#endif
        break;
      }
      case kParamPan:
        mParamsToSmooth[kModPanSmoother] = (T)value / 90. + 1.;
        break;
      case kParamPanEnv1:
      case kParamPanEnv2:
      case kParamPanLFO1:
      case kParamPanLFO2:
      case kParamPanSeq:
      case kParamPanVel:
      case kParamPanKTk:
      case kParamPanRnd:
      {
        const int modIdx = paramIdx - kParamPan;
        // TODO : Make this voice-specific?
        break;
      }
      case kParamEnv1Sustain:
      {
        mParamsToSmooth[kModEnv1SustainSmoother] = (T)value / 100.;
        SendParam([paramIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVEnv1Sustain, 0, (T)value / 100.);
          });
        break;
      }
      case kParamEnv1SustainEnv1:
      case kParamEnv1SustainEnv2:
      case kParamEnv1SustainAmpEnv:
      case kParamEnv1SustainLFO1:
      case kParamEnv1SustainLFO2:
      case kParamEnv1SustainSeq:
      case kParamEnv1SustainVel:
      case kParamEnv1SustainKTk:
      case kParamEnv1SustainRnd:
      {
        const int modIdx = paramIdx - kParamEnv1Velocity;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVEnv1Sustain, modIdx, value);
          });
        break;
      }
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
      case kParamEnv1Velocity:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mEnv1VelocityMod = (T)value;
          });
        break;
      case kParamEnv1DecayCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mEnv1.SetStageCurve(EEnvStage::kDecay, (T)value / 100.);
          });
        break;
      case kParamEnv1ReleaseCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mEnv1.SetStageCurve(EEnvStage::kRelease, (T)value / 100.);
          });
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
      case kParamEnv2Sustain:
      {
        mParamsToSmooth[kModEnv2SustainSmoother] = (T)value / 100.;
        SendParam([paramIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVEnv2Sustain, 0, (T)value / 100.);
          });
        break;
      }
      case kParamEnv2SustainEnv1:
      case kParamEnv2SustainEnv2:
      case kParamEnv2SustainAmpEnv:
      case kParamEnv2SustainLFO1:
      case kParamEnv2SustainLFO2:
      case kParamEnv2SustainSeq:
      case kParamEnv2SustainVel:
      case kParamEnv2SustainKTk:
      case kParamEnv2SustainRnd:
      {
        const int modIdx = paramIdx - kParamEnv2Velocity;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVEnv2Sustain, modIdx, value);
          });
        break;
      }
      case kParamEnv2Velocity:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mEnv2VelocityMod = (T)value;
          });
        break;
      case kParamEnv2DecayCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mEnv2.SetStageCurve(EEnvStage::kDecay, (T)value / 100.);
          });
        break;
      case kParamEnv2ReleaseCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mEnv2.SetStageCurve(EEnvStage::kRelease, (T)value / 100.);
          });
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
      case kParamAmpEnvSustain:
      {
        mParamsToSmooth[kModAmpEnvSustainSmoother] = (T)value / 100.;
        SendParam([paramIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVAmpEnvSustain, 0, (T)value / 100.);
          });
        break;
      }
      case kParamAmpEnvSustainEnv1:
      case kParamAmpEnvSustainEnv2:
      case kParamAmpEnvSustainAmpEnv:
      case kParamAmpEnvSustainLFO1:
      case kParamAmpEnvSustainLFO2:
      case kParamAmpEnvSustainSeq:
      case kParamAmpEnvSustainVel:
      case kParamAmpEnvSustainKTk:
      case kParamAmpEnvSustainRnd:
      {
        const int modIdx = paramIdx - kParamAmpEnvVelocity;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVAmpEnvSustain, modIdx, value);
          });
        break;
      }
      case kParamAmpEnvVelocity:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mAmpEnvVelocityMod = (T)value;
          });
        break;
      case kParamAmpEnvDecayCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mAmpEnv.SetStageCurve(EEnvStage::kDecay, (T)value / 100.);
          });
        break;
      case kParamAmpEnvReleaseCurve:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mAmpEnv.SetStageCurve(EEnvStage::kRelease, (T)value / 100.);
          });
        break;
      case kParamLFO1Amp:
      {
        mGlobalLFO1.SetScalar(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO1.SetScalar(value);
          });
        [[fallthrough]];
      }
      case kParamLFO1AmpEnv1:
      case kParamLFO1AmpEnv2:
      case kParamLFO1AmpAmpEnv:
      case kParamLFO1AmpLFO2:
      case kParamLFO1AmpSeq:
      case kParamLFO1AmpVel:
      case kParamLFO1AmpKTk:
      case kParamLFO1AmpRnd:
      {
        const int modIdx = paramIdx - kParamLFO1Amp;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVLFO1Amp, modIdx, value);
          });
        break;
      }
      case kParamLFO1RateTempo:
      {
        mGlobalLFO1.FastLFO<T>::SetQNScalarFromDivision(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO1.SetQNScalarFromDivision(static_cast<int>(value));
          });
        break;
      }
      case kParamLFO1RateHz:
      {
        mGlobalLFO1.FastLFO<T>::SetFreqCPS(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO1.SetFreqCPS(value);
          });
      }
      case kParamLFO1RateHzEnv1:
      case kParamLFO1RateHzEnv2:
      case kParamLFO1RateHzAmpEnv:
      case kParamLFO1RateHzLFO1:
      case kParamLFO1RateHzLFO2:
      case kParamLFO1RateHzSeq:
      case kParamLFO1RateHzVel:
      case kParamLFO1RateHzKTk:
      case kParamLFO1RateHzRnd:
      {
        const int modIdx = paramIdx - kParamLFO1RateHz;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVLFO1RateHz, modIdx, value);
          });
        break;
      }
      case kParamLFO1RateMode:
      {
        mGlobalLFO1.SetRateMode(value > 0.5);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO1.SetRateMode(value > 0.5);
          });
        break;
      }
      case kParamLFO1Shape:
      {
        mGlobalLFO1.SetShape(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO1.SetShape(static_cast<int>(value));
          });
        break;
      }
      case kParamLFO1Restart:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO1Restart = (value > 0.5);
          });
        break;
      case kParamLFO2Amp:
      {
        mGlobalLFO2.SetScalar(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO2.SetScalar(value);
          });
        [[fallthrough]];
      }
      case kParamLFO2AmpEnv1:
      case kParamLFO2AmpEnv2:
      case kParamLFO2AmpAmpEnv:
      case kParamLFO2AmpLFO1:
      case kParamLFO2AmpLFO2:
      case kParamLFO2AmpSeq:
      case kParamLFO2AmpVel:
      case kParamLFO2AmpKTk:
      case kParamLFO2AmpRnd:
      {
        const int modIdx = paramIdx - kParamLFO2Amp;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVLFO2Amp, modIdx, value);
          });
        break;
      }
      case kParamLFO2RateTempo:
      {
        mGlobalLFO2.FastLFO<T>::SetQNScalarFromDivision(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO2.SetQNScalarFromDivision(static_cast<int>(value));
          });
        break;
      }
      case kParamLFO2RateHz:
      {
        mGlobalLFO2.SetFreqCPS(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO2.SetFreqCPS(value);
          });
        [[fallthrough]];
      }
      case kParamLFO2RateHzEnv1:
      case kParamLFO2RateHzEnv2:
      case kParamLFO2RateHzAmpEnv:
      case kParamLFO2RateHzLFO1:
      case kParamLFO2RateHzLFO2:
      case kParamLFO2RateHzSeq:
      case kParamLFO2RateHzVel:
      case kParamLFO2RateHzKTk:
      case kParamLFO2RateHzRnd:
      {
        const int modIdx = paramIdx - kParamLFO2RateHz;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVLFO2RateHz, modIdx, value);
          });
        break;
      }
      case kParamLFO2RateMode:
      {
        mGlobalLFO2.SetRateMode(value > 0.5);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO2.SetRateMode(value > 0.5);
          });
        break;
      }
      case kParamLFO2Shape:
      {
        mGlobalLFO2.SetShape(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO2.SetShape(static_cast<int>(value));
          });
        break;
      }
      case kParamLFO2Restart:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLFO2Restart = (value > 0.5);
          });
        break;
      case kParamSequencerAmp:
      {
        mGlobalSequencer.SetScalar(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mSequencer.SetScalar(value);
          });
        [[fallthrough]];
      }
      case kParamSequencerAmpEnv1:
      case kParamSequencerAmpEnv2:
      case kParamSequencerAmpAmpEnv:
      case kParamSequencerAmpLFO1:
      case kParamSequencerAmpLFO2:
      case kParamSequencerAmpSeq:
      case kParamSequencerAmpVel:
      case kParamSequencerAmpKTk:
      case kParamSequencerAmpRnd:
      {
        const int modIdx = paramIdx - kParamSequencerAmp;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVSequencerAmp, modIdx, value);
          });
        break;
      }
      case kParamSequencerSteps:
      {
        mGlobalSequencer.SetLength(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mSequencer.SetLength(static_cast<int>(value));
          });
        break;
      }
      case kParamSequencerRateTempo:
      {
        mGlobalSequencer.Sequencer<T>::SetQNScalarFromDivision(static_cast<int>(value));
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mSequencer.SetQNScalarFromDivision(static_cast<int>(value));
          });
        break;
      }
      case kParamSequencerRateHz:
      {
        mGlobalSequencer.SetFreqCPS(value);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mSequencer.SetFreqCPS(value);
          });
        [[fallthrough]];
      }
      case kParamSequencerRateHzEnv1:
      case kParamSequencerRateHzEnv2:
      case kParamSequencerRateHzAmpEnv:
      case kParamSequencerRateHzLFO1:
      case kParamSequencerRateHzLFO2:
      case kParamSequencerRateHzSeq:
      case kParamSequencerRateHzVel:
      case kParamSequencerRateHzKTk:
      case kParamSequencerRateHzRnd:
      {
        const int modIdx = paramIdx - kParamSequencerRateHz;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceModulatorParam(kVSequencerRateHz, modIdx, value);
          });
        break;
      }
      case kParamSequencerRateMode:
      {
        mGlobalSequencer.SetRateMode(value > 0.5);
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mSequencer.SetRateMode(value > 0.5);
          });
        break;
      }
      case kParamSequencerRestart:
      {
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mSequencerRestart = (value > 0.5);
          });
        // Toggle between using a master/static phase to update the Sequencer display, and using the phase of the last-triggered voice
        break;
      }
      case kParamSequencerGlide:
      {
        T glideNorm = (T)value / 100.;
        mGlobalSequencer.SetGlide(glideNorm);
        mSynth.ForEachVoice([glideNorm](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mSequencer.SetGlide(glideNorm);
          });
      }
      case kParamLegato:
      {
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mLegato = value > 0.5;
          });
        if (!(value > 0.5))
        {
          std::fill(AmpEnvQueue.begin(), AmpEnvQueue.end(), nullptr);
          std::fill(Env1Queue.begin(), Env1Queue.end(), nullptr);
          std::fill(Env2Queue.begin(), Env2Queue.end(), nullptr);
        }
        break;
      }
      case kParamUnisonVoices:
      {
        int voices = static_cast<int>(value);
        mDetuner.SetNVoices(voices);
        mSynth.SetMonoUnison(voices);
        // If the unison voice number has increased, start new voices
        if (voices > mUnisonVoices)
        {
          for (int i{ mUnisonVoices }; i < voices; ++i)
          {
            // TODO
          }
        }
        else
        {
          // TODO: When the unison number is descreased, only stop any excess voices. (Will probably require keeping track of unison voices in the VoiceAllocator)
          mSynth.Reset();
        }
        mUnisonVoices = voices;
        break;
      }
      case kParamUnisonDetune:
      {
        mDetuner.SetMaxDetune(value);
        mDetuner.Reset();
        // Send new values to voices
        mSynth.ForEachVoice([](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP<T>::Voice&>(voice).ResetUnisonParams();
          });
        break;
      }
      case kParamUnisonChord:
        mDetuner.SetChord(static_cast<int>(value) + EUnisonChords::kNoChord);
        break;
      case kParamStereoSpread:
      {
        if (mDetuner.mNVoices == 1)
        {
          mStereoWidth = (T)value;
          mDetuner.SetMaxPan(mStereoWidth);
        }
        else
        {
          mStereoWidth = (T)value;
          mDetuner.SetMaxPan(mStereoWidth);
        }
        // Send new values to voices
        mDetuner.Reset();
        mSynth.ForEachVoice([](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP<T>::Voice&>(voice).ResetUnisonParams();
          });
        break;
      }
      case kParamWavetable1Pitch:
        mParamsToSmooth[kModWavetable1PitchSmoother] = value;
        break;
      case kParamWavetable1PitchEnv1:
      case kParamWavetable1PitchEnv2:
      case kParamWavetable1PitchAmpEnv:
      case kParamWavetable1PitchLFO1:
      case kParamWavetable1PitchLFO2:
      case kParamWavetable1PitchSeq:
      case kParamWavetable1PitchVel:
      case kParamWavetable1PitchKTk:
      case kParamWavetable1PitchRnd:
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
      case kParamWavetable1AmpAmpEnv:
      case kParamWavetable1AmpLFO1:
      case kParamWavetable1AmpLFO2:
      case kParamWavetable1AmpSeq:
      case kParamWavetable1AmpVel:
      case kParamWavetable1AmpKTk:
      case kParamWavetable1AmpRnd:
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
      case kParamWavetable1PosAmpEnv:
      case kParamWavetable1PosLFO1:
      case kParamWavetable1PosLFO2:
      case kParamWavetable1PosSeq:
      case kParamWavetable1PosVel:
      case kParamWavetable1PosKTk:
      case kParamWavetable1PosRnd:
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
      case kParamWavetable1BendAmpEnv:
      case kParamWavetable1BendLFO1:
      case kParamWavetable1BendLFO2:
      case kParamWavetable1BendSeq:
      case kParamWavetable1BendVel:
      case kParamWavetable1BendKTk:
      case kParamWavetable1BendRnd:
      {
        const int modIdx = paramIdx - kParamWavetable1Bend;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVWavetable1Bend, modIdx, value);
          });
        break;
      }
      case kParamWavetable1Formant:
        mParamsToSmooth[kModWavetable1FormantSmoother] = value;
        break;
      case kParamWavetable1FormantEnv1:
      case kParamWavetable1FormantEnv2:
      case kParamWavetable1FormantAmpEnv:
      case kParamWavetable1FormantLFO1:
      case kParamWavetable1FormantLFO2:
      case kParamWavetable1FormantSeq:
      case kParamWavetable1FormantVel:
      case kParamWavetable1FormantKTk:
      case kParamWavetable1FormantRnd:
      {
        const int modIdx = paramIdx - kParamWavetable1Formant;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVWavetable1Formant, modIdx, value);
          });
        break;
      }
      case kParamWavetable2:
      {
        break;
      }
      case kParamWavetable2Pitch:
        mParamsToSmooth[kModWavetable2PitchSmoother] = value;
        break;
      case kParamWavetable2PitchEnv1:
      case kParamWavetable2PitchEnv2:
      case kParamWavetable2PitchAmpEnv:
      case kParamWavetable2PitchLFO1:
      case kParamWavetable2PitchLFO2:
      case kParamWavetable2PitchSeq:
      case kParamWavetable2PitchVel:
      case kParamWavetable2PitchKTk:
      case kParamWavetable2PitchRnd:
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
      case kParamWavetable2AmpAmpEnv:
      case kParamWavetable2AmpLFO1:
      case kParamWavetable2AmpLFO2:
      case kParamWavetable2AmpSeq:
      case kParamWavetable2AmpVel:
      case kParamWavetable2AmpKTk:
      case kParamWavetable2AmpRnd:
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
      case kParamWavetable2PosAmpEnv:
      case kParamWavetable2PosLFO1:
      case kParamWavetable2PosLFO2:
      case kParamWavetable2PosSeq:
      case kParamWavetable2PosVel:
      case kParamWavetable2PosKTk:
      case kParamWavetable2PosRnd:
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
      case kParamWavetable2BendAmpEnv:
      case kParamWavetable2BendLFO1:
      case kParamWavetable2BendLFO2:
      case kParamWavetable2BendSeq:
      case kParamWavetable2BendVel:
      case kParamWavetable2BendKTk:
      case kParamWavetable2BendRnd:
      {
        const int modIdx = paramIdx - kParamWavetable2Bend;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVWavetable2Bend, modIdx, value);
          });
        break;
      }
      case kParamWavetable2Formant:
        mParamsToSmooth[kModWavetable2FormantSmoother] = value;
        break;
      case kParamWavetable2FormantEnv1:
      case kParamWavetable2FormantEnv2:
      case kParamWavetable2FormantAmpEnv:
      case kParamWavetable2FormantLFO1:
      case kParamWavetable2FormantLFO2:
      case kParamWavetable2FormantSeq:
      case kParamWavetable2FormantVel:
      case kParamWavetable2FormantKTk:
      case kParamWavetable2FormantRnd:
      {
        const int modIdx = paramIdx - kParamWavetable2Formant;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVWavetable2Formant, modIdx, value);
          });
        break;
      }
      case kParamFilter1Type:
      {
        mFilter1Comb = static_cast<int>(value) == kComb;
        SendParam([value](Voice* voice) {
          voice->SetFilterType(0, static_cast<int>(value));
          });
        break;
      }
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
        if (mFilter1Comb)
          mParamsToSmooth[kModFilter1CutoffSmoother] = value / 100.;
        break;
      }
      case kParamFilter1CutoffEnv1:
      case kParamFilter1CutoffEnv2:
      case kParamFilter1CutoffAmpEnv:
      case kParamFilter1CutoffLFO1:
      case kParamFilter1CutoffLFO2:
      case kParamFilter1CutoffSeq:
      case kParamFilter1CutoffVel:
      case kParamFilter1CutoffKTk:
      case kParamFilter1CutoffRnd:
      {
        const int modIdx = paramIdx - kParamFilter1Cutoff;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVFilter1Cutoff, modIdx, value);
          });
        break;
      }
      case kParamFilter1Resonance:
        mParamsToSmooth[kModFilter1ResonanceSmoother] = value / 100.;
        break;
      case kParamFilter1FB:
      {
        if (mFilter1Comb)
          mParamsToSmooth[kModFilter1ResonanceSmoother] = value / 100.;
        break;
      }
      case kParamFilter1ResonanceEnv1:
      case kParamFilter1ResonanceEnv2:
      case kParamFilter1ResonanceAmpEnv:
      case kParamFilter1ResonanceLFO1:
      case kParamFilter1ResonanceLFO2:
      case kParamFilter1ResonanceSeq:
      case kParamFilter1ResonanceVel:
      case kParamFilter1ResonanceKTk:
      case kParamFilter1ResonanceRnd:
      {
        const int modIdx = paramIdx - kParamFilter1Resonance;
        SendParam([paramIdx, modIdx, value](Voice* voice) {
          voice->UpdateVoiceParam(kVFilter1Resonance, modIdx, value);
          });
        break;
      }
      case kParamFilter1Drive:
      case kParamFilter1Delay:
      {
        mParamsToSmooth[kModFilter1DriveSmoother] = value / (mFilter1Comb ? (double)COMB_MAX_DELAY : 100.);
        break;
      }
      case kParamFilter1DriveEnv1:
      case kParamFilter1DriveEnv2:
      case kParamFilter1DriveAmpEnv:
      case kParamFilter1DriveLFO1:
      case kParamFilter1DriveLFO2:
      case kParamFilter1DriveSeq:
      case kParamFilter1DriveVel:
      case kParamFilter1DriveKTk:
      case kParamFilter1DriveRnd:
      {
        const int modIdx = paramIdx - kParamFilter1Drive;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVFilter1Drive, modIdx, value);
          });
        break;
      }
      case kParamFilter2Type:
      {
        mFilter2Comb = static_cast<int>(value) == kComb;
        SendParam([value](Voice* voice) {
          voice->SetFilterType(1, static_cast<int>(value));
          });
        break;
      }
      case kParamFilter2ModeVSF:
      case kParamFilter2ModeMoog:
        mSynth.ForEachVoice([value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).mFilters.at(1)->SetMode(static_cast<int>(value));
          });
        break;
      case kParamFilter2Cutoff:
        mParamsToSmooth[kModFilter2CutoffSmoother] = value / (mFilter2Comb ? 100. : mSampleRate);
        break;
      case kParamFilter2FF:
      {
        if (mFilter2Comb)
          mParamsToSmooth[kModFilter2CutoffSmoother] = value / (mFilter2Comb ? 100. : mSampleRate);
        break;
      }
      case kParamFilter2CutoffEnv1:
      case kParamFilter2CutoffEnv2:
      case kParamFilter2CutoffAmpEnv:
      case kParamFilter2CutoffLFO1:
      case kParamFilter2CutoffLFO2:
      case kParamFilter2CutoffSeq:
      case kParamFilter2CutoffVel:
      case kParamFilter2CutoffKTk:
      case kParamFilter2CutoffRnd:
      {
        const int modIdx = paramIdx - kParamFilter2Cutoff;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVFilter2Cutoff, modIdx, value);
          });
        break;
      }
      case kParamFilter2Resonance:
        mParamsToSmooth[kModFilter2ResonanceSmoother] = value / 100.;
        break;
      case kParamFilter2FB:
      {
        if (mFilter2Comb)
          mParamsToSmooth[kModFilter2ResonanceSmoother] = value / 100.;
        break;
      }
      case kParamFilter2ResonanceEnv1:
      case kParamFilter2ResonanceEnv2:
      case kParamFilter2ResonanceAmpEnv:
      case kParamFilter2ResonanceLFO1:
      case kParamFilter2ResonanceLFO2:
      case kParamFilter2ResonanceSeq:
      case kParamFilter2ResonanceVel:
      case kParamFilter2ResonanceKTk:
      case kParamFilter2ResonanceRnd:
      {
        const int modIdx = paramIdx - kParamFilter2Resonance;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVFilter2Resonance, modIdx, value / 100.);
          });
        break;
      }
      case kParamFilter2Drive:
      case kParamFilter2Delay:
        mParamsToSmooth[kModFilter2DriveSmoother] = value / (mFilter2Comb ? (double)COMB_MAX_DELAY : 100.);
        break;
      case kParamFilter2DriveEnv1:
      case kParamFilter2DriveEnv2:
      case kParamFilter2DriveAmpEnv:
      case kParamFilter2DriveLFO1:
      case kParamFilter2DriveLFO2:
      case kParamFilter2DriveSeq:
      case kParamFilter2DriveVel:
      case kParamFilter2DriveKTk:
      case kParamFilter2DriveRnd:
      {
        const int modIdx = paramIdx - kParamFilter2Drive;
        mSynth.ForEachVoice([paramIdx, modIdx, value](SynthVoice& voice) {
          dynamic_cast<TablitsaDSP::Voice&>(voice).UpdateVoiceParam(kVFilter2Drive, modIdx, value);
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
      case kParamPhaseModFreqVel:
      case kParamPhaseModFreqKTk:
      case kParamPhaseModFreqRnd:
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
      case kParamPhaseModAmountVel:
      case kParamPhaseModAmountKTk:
      case kParamPhaseModAmountRnd:
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
      case kParamRingModFreqVel:
      case kParamRingModFreqKTk:
      case kParamRingModFreqRnd:
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
      case kParamRingModAmountVel:
      case kParamRingModAmountKTk:
      case kParamRingModAmountRnd:
      {
        const int modIdx = paramIdx - kParamRingModAmount;
        SendParam([value, modIdx](Voice* voice) {
          voice->UpdateVoiceParam(kVRingModAmt, modIdx, value);
          });
        break;
      }
      case kParamDelayTimeMode:
        mDelayEffect.SetTempoSync(value > 0.5);
        break;
      case kParamDelayTimeLMilliseconds:
      case kParamDelayTimeRMilliseconds:
        mDelayEffect.SetDelayMS(value, paramIdx - kParamDelayTimeLMilliseconds);
        break;
      case kParamDelayTimeLBeats:
      case kParamDelayTimeRBeats:
      {
        double qnScalar = LFO<T>::GetQNScalar(static_cast<LFO<T>::ETempoDivision>(Clip((int)value, 0, (int)LFO<T>::ETempoDivision::kNumDivisions)));
        double qnPerMeasure = 4. / mTSDenom * mTSNum;
        mDelayEffect.SetDelayTempo(1. / qnScalar / qnPerMeasure, paramIdx - kParamDelayTimeLBeats, mTempo);
        break;
      }
      case kParamDelayFeedback:
        mDelayEffect.SetFeedback((T)value / 100.);
        break;
      case kParamDelayMix:
        mDelayEffect.SetGain((T)value / 100.);
        break;
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
  std::vector<std::string> mWavetables{ ELEMENT_NAMES };
  std::vector<Voice*> mSynthVoices;

  // Polyphonic/Monophonic
  bool mMono{ false };
  int mUnisonVoices{ 1 };
  float mUnisonDetune{ 0.f };
  UnisonVoiceManager mDetuner{ kMaxUnisonVoices };
  double mStereoWidth{ 0. };

  // Portamento Parameters
  bool mConstantGlideTime{ true };
  double mGlideRateScalar{ 1. }; // Semitones per second
  double mLastNoteOn{ 0 };

  // Audio processing and musical parameters
  double mSampleRate{ DEFAULT_SAMPLE_RATE };
  int mTSNum{ 4 };
  int mTSDenom{ 4 };
  double mTempo{ DEFAULT_TEMPO };

  // Status Variables
  int mSeqPos{ 0 };
  Sequencer<T, kNumSeqSteps>* mActiveSequencer{ nullptr };
  int mFilter1Osc{ 0 }; // Oscillator which provides the filter's input
  int mFilter2Osc{ 1 };
  bool mFilter1Comb{ false }; // Set to `true` when Filter 1 is a comb filter. Used for scaling delay/drive values by the proper amount.
  bool mFilter2Comb{ false }; // Set to `true` when Filter 2 is a comb filter. Used for scaling delay/drive values by the proper amount.

  // Non-modulatable parameters
  double mLoadedWavetables[2]{ 1., 2. }; // Integer indices (stored as double) of current wavetables
  double mSeqSteps[kNumSeqSteps]{}; // Value of each step in the sequencer
  int mStepPos{ 0 };
  int mPrevPos{ -1 };

  // Effects
  std::vector<Effect<T>*> mEffects;
  DelayEffect<T> mDelayEffect{ DEFAULT_SAMPLE_RATE, DEFAULT_SAMPLE_RATE * 12. };

  // Pointers to master modulators, for free-run and legato modes
  std::vector<Envelope<T>*> Env1Queue;
  std::vector<Envelope<T>*> Env2Queue;
  std::vector<Envelope<T>*> AmpEnvQueue;
  FastLFO<T>* mMasterLFO1{ nullptr }; // The last-triggered `mLFO1`, which "owns" the master phase
  FastLFO<T>* mMasterLFO2{ nullptr }; // The last-triggered `mLFO2`, which "owns" the master phase
  Sequencer<T>* mMasterSeq{ nullptr }; // The last-triggered `mSequencer`, which "owns" the master phase

  // Global Modulators
  ModMetronome mGlobalMetronome;
  GlobalModulator<T, FastLFO<T>> mGlobalLFO1{ &mGlobalMetronome };
  GlobalModulator<T, FastLFO<T>> mGlobalLFO2{ &mGlobalMetronome };
  GlobalModulator<T, Sequencer<T>> mGlobalSequencer{ &mGlobalMetronome, mSeqSteps };
};