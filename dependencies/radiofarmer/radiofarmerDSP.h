
#include "vectorclass.h"

#include "radiofarmer_config.h"

#include "allpass.h"
#include "delayline.h"
#include "cm.h"
#include "reverb.h"
#include "waveguide.h"
#include "waveshaper.h"

/*
class VectorOscillator
{
  union tabfudge
  {
    double d;
    int i[2];
  }; // todo: aligned

public:
  VectorOscillator(double startPhase = 0.) : mPhase(startPhase) {}

  inline Vec4d __vectorcall ProcessBlock4_Vector()
  {
    double phase = mPhase + (double)UNITBIT32;
    const double phaseIncr = mPhaseIncr * tableSize;
    Vec4d vPhase = phase + phaseIncr * mIncrVector;

    // Read 4 doubles as 8 integers (signed ints are used, but the AND operations later make this irrelevant)
    Vec8i viPhase = reinterpret_i(vPhase);
    // Upper 32 bits of 3*2^19 in ______ indices, 0xFFFF in _____: i.e. 0xFFFF, 0x18, 0xFFFF, ...
    Vec8i normhipart = blend8<HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V, 8>(reinterpret_i(Vec4d((double)UNITBIT32)), Vec8i(0xFFFFFFFF));
    // Mask the 8-item vector of 32-bit ints with one less than the table size, pad the upper bits (lower indices) with zeros, and reinterpret as a 4-item vector of 64-bit ints
    Vec8i offsets32 = viPhase & tableSizeM1;
    Vec4q offsets = reinterpret_i(permute8<HIOFFSET_V, -1, HIOFFSET_V + 2, -1, HIOFFSET_V + 4, -1, HIOFFSET_V + 6, -1>(offsets32));
    // Force the double to wrap. (AND the upper bits with the upper 32 bits of UNITBIT32)
    viPhase &= normhipart;
    Vec4d frac = reinterpret_d(viPhase) - (double)UNITBIT32; // get the fractional part
    const Vec4d f1 = lookup<tableSize>(offsets, mLUT);
    const Vec4d f2 = lookup<tableSize>(offsets, mLUT + 1);
    Vec4d output = mul_add((f2 - f1), frac, f1);

    // Restore mPhase
    tabfudge tf;
    phase += phaseIncr * 4.;
    tf.d = UNITBIT32 * tableSize;
    const int normhipart2 = tf.i[HIOFFSET_V];
    tf.d = phase + (UNITBIT32 * tableSize - UNITBIT32); // Remove the offset we introduced at the start of UNITBIT32.
    tf.i[HIOFFSET_V] = normhipart2;
    mPhase = tf.d - UNITBIT32 * tableSize;
    return output;
  }

  inline Vec4d __vectorcall Process_Vector()
  {
    //    double output[4];
    //    FastSinOscillator<T>::ProcessBlock(output, 4);
    //    return Vec4d().load_a(output);
    return ProcessBlock4_Vector();
  }

private:
  double mPhase;
  double mPhaseIncr;

  const Vec4d mIncrVector = Vec4d(0., 1., 2., 3.);

  static constexpr sample mLUT[512];
  static constexpr int tableSize{ 512 };
  static constexpr int tableSizeM1{ 511 };
}; //todo: aligned*/