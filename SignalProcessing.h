#pragma once

#include "VectorFunctions.h"

#include <cmath>
#include <assert.h>

template<typename T>
extern inline T SoftClip(T s, T gain = (T)1);

template<typename T1, typename T2>
extern inline T1 __vectorcall SoftClip(const T1& s, T2 gain = (T2)1);


