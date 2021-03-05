#include "VectorFunctions.h"

void transpose16f(float matrix[4][4])
{
  Vec16f r, c;
  r.load(&matrix[0][0]);
  c = permute16<0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15>(r);
  c.store(&matrix[0][0]);
}

void transpose8d(double matrix[4][2])
{
  Vec8d r, c;
  r.load(&matrix[0][0]);
  c = permute8<0, 4, 1, 5, 2, 6, 3, 7>(r);
  c.store(&matrix[0][0]);
}

extern inline Vec4q __vectorcall to_int64_in_range(Vec4d v, int scale)
{
  v = v * scale + (double)UNITBIT32;
  const int scaleM1 = scale - 1;

  // Read 4 doubles as 8 integers (signed ints are used, but the AND operations later make this irrelevant)
  Vec8i v_reinterp = reinterpret_i(v);
  Vec8i normhipart = blend8<8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V, 8, HIOFFSET_V>(reinterpret_i(Vec4d((double)UNITBIT32)), Vec8i(0));
  // Mask the 8-item vector of 32-bit ints with one less than the table size, pad the upper bits (lower indices) with zeros, and reinterpret as a 4-item vector of 64-bit ints
  Vec8i vInt32 = v_reinterp & scaleM1;
  Vec4q vInt64 = reinterpret_i(permute8<HIOFFSET_V, -1, HIOFFSET_V + 2, -1, HIOFFSET_V + 4, -1, HIOFFSET_V + 6, -1>(vInt32));
  return vInt64;
}

template<>
extern inline Vec4d __vectorcall sign<Vec4d>(const Vec4d& v)
{
  constexpr int32_t sign_bit = 0x80000000; // Bit mask for the sign bit in the upper 32 bits of a 64-bit float
  constexpr int32_t minus1_float64_upper32 = 0xbff00000; // Upper 32 bits of the representation of -1 as a 32-bit float
  constexpr int32_t msb_zero = 0x7fffffff; // Zero in most-significant bit, one in all other bits

  // Double-precision floating point representation of -1., stored as 32-bit integers
  const Vec8i minus1_vec = Vec8i(0, minus1_float64_upper32, 0, minus1_float64_upper32, 0, minus1_float64_upper32, 0, minus1_float64_upper32);

  Vec8i v_reinterp = reinterpret_i(v); // Get 4 64-bit floats as 8 32-bit integers
  Vec8i sign_bits = v_reinterp & Vec8i(0, sign_bit, 0, sign_bit, 0, sign_bit, 0, sign_bit); // AND-out everything but the most significant bit (the sign)
  Vec8i sign_scalar_int = (sign_bits | msb_zero) & minus1_vec; // OR the sign bits with a zero followed by 63 ones, then AND that with floating-point -1
  Vec4d sign_scalar = reinterpret_d(sign_scalar_int);
  return sign_scalar;
}