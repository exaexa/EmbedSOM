#ifndef use_intrins_h

#ifdef __SSE4_2__
#define USE_INTRINS
#endif

#ifdef USE_INTRINS
#include <xmmintrin.h>
#endif

#endif // use_intrins_h
