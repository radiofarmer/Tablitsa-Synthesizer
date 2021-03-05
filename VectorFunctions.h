#pragma once
#include <vectorclass.h>

#define UNITBIT32 1572864.
#define HIOFFSET_V 1

template<class V>
struct deduce_vector_from
{
};

template<>
struct deduce_vector_from<double>
{
  typedef Vec2d double_vec2;
  typedef Vec4d double_vec4;
  typedef Vec8d double_vec8;
};

template<>
struct deduce_vector_from<float>
{
  typedef Vec4f float_vec4;
  typedef Vec8f float_vec8;
  typedef Vec16f float_vec16;
};

template<>
struct deduce_vector_from<Vec8d>
{
  typedef Vec8q int_vec;
  typedef Vec16f float_vec;
  static constexpr int v_size = 8;
};

template<>
struct deduce_vector_from<Vec4d>
{
  typedef Vec4q int_vec;
  typedef Vec8f float_vec;
  static constexpr int v_size = 4;
};

template<>
struct deduce_vector_from<Vec8f>
{
  typedef Vec4i int_vec;
  static constexpr int v_size = 8;
};

void transpose16f(float matrix[4][4]);

void transpose4d(double matrix[2][2]);

/* Get the lower 32 bits of a double-precision floating-point vector with absolute values between 0 and 2^19 as a vector fo 64-bit integers betwee 0 and `scale`. */
extern inline Vec4q __vectorcall to_int64_in_range(Vec4d v, int scale);

template<typename Vd>
extern inline Vd __vectorcall sign(const Vd& v);


template<typename T, class V>
extern inline void IIR_2pole_coefficients_4(T a0, T a1, T a2, V& c1, V& c2, V& c3, V& c4, V& c5, V& c6);

template<typename V>
extern inline V __vectorcall IIR_2pole_4(const V& x, V& z, const V* c, const V& b);