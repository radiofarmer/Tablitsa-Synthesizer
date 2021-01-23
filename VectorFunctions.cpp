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