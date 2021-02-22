#include <vectorclass.h>

#define UNITBIT32 1572864.
#define HIOFFSET_V 1

void transpose16f(float matrix[4][4]);

void transpose4d(double matrix[2][2]);

extern inline Vec4q __vectorcall to_int64_in_range(Vec4d v, int scale);