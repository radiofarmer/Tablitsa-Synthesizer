#include "reverb.h"

#include <cmath>

using namespace radiofarmer;

UDFNReverb::UDFNReverb(sample_t sampleRate, sample_t lpFreq) : mSampleRate(sampleRate)
{
  // Filters
  SetHFDecay(lpFreq);
  SetAAPFilters();
  SetSampleRate(mSampleRate);
}

void UDFNReverb::SetAAPFilters()
{
  // Early Allpass
  for (int i{ 0 }; i < 2; ++i)
  {
    mEarlyAP[i].SetLPGain(m_a);
    mEarlyAP[i].SetFeedbackGain(m_g);
  }

  // Late Allpass
  for (int i{ 0 }; i < 2; ++i)
  {
    for (int j{ 0 }; j < NLateAPFilters; ++j)
    {
      mLateAP[i][j].SetLPGain(m_a);
      mLateAP[i][j].SetFeedbackGain(m_g);
    }
  }
}

void UDFNReverb::SetSampleRate(sample_t sampleRate)
{
  mSampleRate = sampleRate;

  mEarlyAP[0].SetSampleRate(mSampleRate);
  mEarlyAP[1].SetSampleRate(mSampleRate);

  for (int i{ 0 }; i < 2; ++i)
  {
    for (int j{ 0 }; j < NLateAPFilters; ++j)
    {
      mLateAP[i][j].SetSampleRate(mSampleRate);
    }
  }

  // Set sample delay times
  SetDelays();
}

void UDFNReverb::SetDelays()
{
  // Early Delays
  for (int i{ 0 }; i < 2; ++i)
  {
    for (int j{ 0 }; j < NEarlyDelays; ++j)
    {
      mEarlyDelays[i][j].SetDelay(static_cast<int>(mEarlyDelayTimes[i] * mSampleRate / 1000.));
    }
  }
  // Early Allpass
  for (int i{ 0 }; i < 2; ++i)
  {
    mEarlyAP[i].SetDelayMS(mEarlyDelayTimes[NEarlyDelays]);
    mEarlyAP[i].SetLPF(0., mLPCutoff, 2.);
  }

  // Late Allpass
  for (int i{ 0 }; i < 2; ++i)
  {
    for (int j{ 0 }; j < NLateAPFilters; ++j)
    {
      const int delayIdx = j + static_cast<int>(j >= LateDelayPos) + 1;
      mLateAP[i][j].SetDelayMS(mLateDelayTimes[delayIdx]);
      mLateAP[i][j].SetLPF(0., mLPCutoff, 2.);
    }
  }

  // Late delays
  for (int i{ 0 }; i < 2; ++i)
  {
    mLateDelays[i].SetDelay(static_cast<int>(mLateDelayTimes[LateDelayPos] * mSampleRate / 1000.));
  }
}

/* Reverberation Parameters */

// Diffusion
void UDFNReverb::SetDiffusion(const sample_t diff)
{
  m_g = diff;
  CalcNorm();
  SetAAPFilters();
}

// Damping (HF attenuation)
void UDFNReverb::SetDamping(const sample_t damp)
{
  m_a = 1. - std::max(damp * 0.99, 0.01);
  CalcNorm();
  SetAAPFilters();
}

void UDFNReverb::SetColor(const sample_t freqNorm)
{
  SetHFDecay(freqNorm);
  for (int i{ 0 }; i < NLateAPFilters + 1; ++i)
  {
    m_t[i] = 1. + std::pow(static_cast<sample_t>(i) / static_cast<sample_t>(NLateAPFilters + 2) - freqNorm, 2.);
  }
}

void UDFNReverb::SetHFDecay(const sample_t freqNorm)
{
  mLPCutoff = freqNorm;
  sample_t slope = 2.;
  for (int i{ 0 }; i < 2; ++i)
  {
    sample_t fc = std::clamp(freqNorm, 0.001, 0.95) * 0.5;
    mEarlyLPF[i].CalculateCoefficients(0., fc, slope);
    for (int j{ 0 }; j < NLateAPFilters; ++j)
    {
      mLateAP[i][j].SetLPF(0., fc, slope);
      fc *= 0.95;
      if (j == LateDelayPos)
      {
        mLateLPF[i].CalculateCoefficients(0., fc, slope);
        fc *= 0.95;
      }
    }
  }
}

void UDFNReverb::SetDecayTime(const sample_t tNorm)
{
  sample_t scale = 1. + tNorm;
  for (int i{ 0 }; i < NEarlyDelays + 1; ++i)
  {
    mEarlyDelayTimes[i] = DefaultEarlyDelayTimes[i] * scale * Primes[i] / 10.;
  }

  for (int i{ 0 }; i < NLateAPFilters + 1; ++i)
  {
    mLateDelayTimes[i] = DefaultLateDelayTimes[i] * scale * Primes[i] / 5.;
  }
  SetDelays();
}

void UDFNReverb::CalcNorm()
{
  // This would need to be an array if the different AAPF's have different lowpass gains
  sample_t c_i = m_g * m_g + (1 - m_g * m_g) * ((m_a * m_a) / (1 - m_a * m_a * m_g * m_g));

  sample_t A = 2 * m_g * m_g * std::pow(c_i, static_cast<sample_t>(NLateAPFilters + 1)); // assumes mono parameters
  // Add up energy gain values from filters 1-4
  sample_t B{};
  for (int t{ 0 }; t < NLateAPFilters - 1; ++t)
  {
    B += m_t[t] * m_t[t] * std::pow(c_i, static_cast<sample_t>(t));
  }
  // Account for the delay gain added between filters 4 and 5
  B += m_t[4] * m_t[4] * m_g * std::pow(c_i, static_cast<sample_t>(4));
  B += m_t[5] * m_t[5] * m_g * m_g * std::pow(c_i, static_cast<sample_t>(5));

  mNorm = (sample_t(1.) - A) / B;
}

/* Processing */

void UDFNReverb::ProcessStereo(StereoSample<sample_t>& s)
{
  /* Early Reflections */
  sample_t s_in[]{
    mEarlyLPF[0].Process(s.l),
    mEarlyLPF[1].Process(s.r)
  };
  // Early delays
  sample_t delays[2]{};
  for (int i{ 0 }; i < 2; ++i)
  {
    for (int j{ 0 }; j < NEarlyDelays; ++j)
    {
      delays[i] += mEarlyDelays[i][j].at(0);
      mEarlyDelays[i][j].push(s_in[i]);
    }
  }

  // Late Reverb
  // First four allpass filters
  sample_t late_in[]{ delays[0], delays[1] };
  sample_t late_out[2]{};
  for (int i{ 0 }; i < 2; ++i)
  {
    for (int j{ 0 }; j < NLateAPFilters - 1; ++j)
    {
      late_in[i] = mLateAP[i][j].Process(late_in[i]);
      late_out[i] += late_in[i] * m_t[i];
    }
  }

  // Delay and LPF
  for (int i{ 0 }; i < 2; ++i)
  {
    sample_t temp = late_in[i];
    late_in[i] = mLateLPF[i].Process(mLateDelays[i].at(0)) * mG[i];
    mLateDelays[i].push(temp);

    // Last APF (run twice)
    late_in[i] = mLateAP[i][NLateAPFilters - 1].Process(late_in[i]);
    late_out[i] += late_in[i] * m_t[NLateAPFilters - 1];
    late_in[i] = mLateAP[i][NLateAPFilters - 1].Process(late_in[i]);
    late_out[i] += late_in[i] * m_t[NLateAPFilters];
    late_out[i] *= mNorm;
  }

  StereoSample<sample_t> out{ late_out[0] + mEarlyAP[0].Process(delays[0]) * mERMix,
    late_out[1] + mEarlyAP[1].Process(delays[1]) * mERMix };

    s = s + (out - s) * mMix;
}