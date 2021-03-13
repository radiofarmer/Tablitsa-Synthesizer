#pragma once

#include "radiofarmer_config.h"

template<typename T, int order=5>
T SoftClip(T x, T threshold=1.);