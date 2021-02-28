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