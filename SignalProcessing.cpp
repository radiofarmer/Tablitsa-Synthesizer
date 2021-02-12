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
