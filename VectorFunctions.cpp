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


template<>
extern inline void IIR_2pole_coefficients_4(double a0, double a1, double a2, Vec4d& c1, Vec4d& c2, Vec4d& c3, Vec4d& c4, Vec4d& c5, Vec4d& c6)
{
  double a_rec1[4]{
    1.,
    a1,
    a1 * a1 + a2,
    a1 * a2
  };
  double a_rec2[4]{
    a1,
    a1 * a1 + a2,
    a1 * a1 * a1 + 2 * a1 * a2,
    a1 * a1 * a2 + a2 * a2
  };
  double a_rec3[4]{
    a1 * a1 + a2,
    a1 * a1 * a1 + 2 * a1 * a2,
    a1 * a1 * a1 * a1 + 3 * a1 * a1 * a2 + a2 * a2,
    a1 * a1 * a1 * a2 + 2 * a1 * a2 * a2
  };
  // Column vectors of coefficients (ALG2)
  c1 = Vec4d(a_rec3[0], a_rec2[0], a_rec1[0], 0.) * a0; // multiplied by x[n + 1]
  c2 = Vec4d(a_rec3[1], a_rec2[1], a_rec1[1], 1.) * a0; // " " x[n]
  c3 = Vec4d(a_rec3[2], a_rec2[2], a_rec1[2], a1); // " " w[n - 1]
  c4 = Vec4d(a_rec3[3], a_rec2[3], a_rec1[3], a2); // " " w[n - 2]

  c5 = Vec4d(a1, 1., 0., 0.) * a0; // multiplied by x[n + 2]
  c6 = Vec4d(a0, 0., 0., 0.); // " " x[n + 3]
}

/** Evaluate 4 samples of an a 2-pole IIR filter
* @param x Four input samples
* @param z Delay line (only the first two items are used)
* @param a Feedback coefficients
* @param b Feedforward coefficients
*/
template<>
extern inline Vec4d __vectorcall IIR_2pole_4(const Vec4d& x, Vec4d& z, const Vec4d* c, const Vec4d& b)
{
  // All delay terms:
  Vec4d y_in = c[5] * x[3] + c[4] * x[2] + c[0] * x[1] + c[1] * x[0] + c[2] * z[0] + c[3] * z[1];

  const Vec4d y_out = Vec4d(
    horizontal_add(blend4<3, 4, 5, V_DC>(y_in, z) * b), // y[n] = b0 * w[n] + b1 * w[n-1] + b2 * w[n-2]
    horizontal_add(blend4<2, 3, 4, V_DC>(y_in, z) * b), // y[n+1] = b0 * w[n+1] + b1 * w[n] + b2 * w[n-1]
    horizontal_add(permute4<1, 2, 3, V_DC>(y_in) * b),  // y[n+2] = b0 * w[n+2] + b1 * w[n + 1] + b2 * w[n]
    horizontal_add(permute4<0, 1, 2, V_DC>(y_in) * b) // y[n+3] = b0 * w[n+3] + b1 * w[n+2] + b2 * w[n+1]
  ); // earliest -> latest

  z = permute4<0, 1, V_DC, V_DC>(y_in); // w[n+3] -> w[n-1]; w[n+2] -> w[n-2]; Order is latest -> earliest
  return y_out;
}