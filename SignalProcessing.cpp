#include "SignalProcessing.h"

template<typename T>
DelayLine<T>::~DelayLine()
{
  delete mBuffer;
};

/*
Push a new value onto the delay line and adjust the read and write heads
*/
template<typename T>
void DelayLine<T>::push(T s)
{
  mWrite++ = s;
  mRead++;
  if (mWrite == mLength)
    mWrite = mBuffer;
  if (mRead == mLength)
    mRead = mBuffer;
}

template<typename T>
T DelayLine<T>::at(size_t offset)
{
  return (mRead - mBuffer)[offset];
}