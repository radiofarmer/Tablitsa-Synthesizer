#include "SignalProcessing.h"

#include <assert.h>

DelayLine::DelayLine(const int length) : mLength(length)
{
  assert(mLength > 1);
  mBuffer = new double[mLength] {};
}

void DelayLine::reset()
{
  delete[] mBuffer;
  mBuffer = new double[mLength] {};
}

void DelayLine::SetDelay(const int d)
{
  mLength = d;
}

inline double DelayLine::at(const int offset)
{
  int readPoint{ mRead - offset };
  if (readPoint < 0)
    readPoint += mLength;
  return mBuffer[readPoint];
}

double DelayLine::operator[](const int idx)
{
  return at(idx);
}
