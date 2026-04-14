#ifndef LEAP_COMPAT_H_
#define LEAP_COMPAT_H_

#if defined(__cplusplus)
#define LEAP_ALIGNAS(bytes) alignas(bytes)
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
#define LEAP_ALIGNAS(bytes) _Alignas(bytes)
#else
#define LEAP_ALIGNAS(bytes) __attribute__((aligned(bytes)))
#endif

#endif /* LEAP_COMPAT_H_ */
