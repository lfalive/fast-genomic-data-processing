/* The MIT License

   Copyright (c) 2008, by Attractive Chaos <attractivechaos@aol.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
  2008-11-16 (0.1.4):

    * Fixed a bug in introsort() that happens in rare cases.

  2008-11-05 (0.1.3):

    * Fixed a bug in introsort() for complex comparisons.

	* Fixed a bug in mergesort(). The previous version is not stable.

  2008-09-15 (0.1.2):

	* Accelerated introsort. On my Mac (not on another Linux machine),
	  my implementation is as fast as std::sort on random input.

	* Added combsort and in introsort, switch to combsort if the
	  recursion is too deep.

  2008-09-13 (0.1.1):

	* Added k-small algorithm

  2008-09-05 (0.1.0):

	* Initial version

*/

#ifndef AC_KSORT_H
#define AC_KSORT_H

#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <memory>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

template <typename T, typename Compare>
static inline void insertsort_new(T* s, T* t, Compare cmp) {
    for (T* i = s + 1; i < t; ++i) {
        for (T* j = i; j > s && cmp(*j, *(j - 1)); --j) {
            std::swap(*j, *(j - 1));
        }
    }
}

template <typename T, typename Compare>
void combsort_new(size_t n, T* a, Compare cmp) {
    const double shrink_factor = 1.2473309501039786540366528676643;
    bool do_swap;
    size_t gap = n;

    do {
        if (gap > 2) {
            gap = static_cast<size_t>(gap / shrink_factor);
            if (gap == 9 || gap == 10) gap = 11;
        }
        do_swap = false;
        for (T* i = a; i < a + n - gap; ++i) {
            T* j = i + gap;
            if (cmp(*j, *i)) {
                std::swap(*i, *j);
                do_swap = true;
            }
        }
    } while (do_swap || gap > 2);

    if (gap != 1) insertsort_new(a, a + n, cmp);
}

template <typename T, typename Compare>
void introsort_new(size_t n, T* a, Compare cmp) {
    if (n <= 1) return;

    int d = 2;
    T rp,*s, *t, *i, *j, *k;
    struct ks_isort_stack_t {
        T *left, *right;
        int depth;
    } *top, *stack;

    if (n == 2) {
        if (cmp(a[1], a[0])) std::swap(a[0], a[1]);
        return;
    }

    while (1ul << d < n) ++d;
    stack = new ks_isort_stack_t[(sizeof(size_t) * d) + 2];
    top = stack; s = a; t = a + (n - 1); d <<= 1;

    while (true) {
        if (s < t) {
            if (--d == 0) {
                combsort_new(t - s + 1, s, cmp);
                t = s;
                continue;
            }

            i = s; j = t; k = i + ((j - i) >> 1) + 1;

            if (cmp(*k, *i)) {
                if (cmp(*k, *j)) k = j;
            } else k = cmp(*j, *i) ? i : j;
            rp = *k;
            if (k != t) std::swap(*k, *t);

            for (;;) {
                do ++i; while (cmp(*i, rp));
                do --j; while (i <= j && cmp(rp, *j));
                if (j <= i) break;
                std::swap(*i, *j);
            }
            std::swap(*i, *t);

            if (i - s > t - i) {
                if (i - s > 16) {
                    top->left = s; top->right = i - 1; top->depth = d; ++top;
                }
                s = t - i > 16 ? i + 1 : t;
            } else {
                if (t - i > 16) {
                    top->left = i + 1; top->right = t; top->depth = d; ++top;
                }
                t = i - s > 16 ? i - 1 : s;
            }
        } else {
            if (top == stack) {
                delete[] stack;
                insertsort_new(a, a + n, cmp);
                return;
            } else {
                --top; s = top->left; t = top->right; d = top->depth;
            }
        }
    }
}

typedef struct {
    void *left, *right;
    int depth;
} ks_isort_stack_t;

#define KSORT_INIT(name, type_t, __sort_lt)								\
	static inline void __ks_insertsort_##name(type_t *s, type_t *t)		\
	{																	\
		type_t *i, *j;										\
		for (i = s + 1; i < t; ++i)										\
			for (j = i; j > s && __sort_lt(*j, *(j-1)); --j) {          \
                std::swap(*j, *(j-1));                                  \
			}															\
	}																	\
	void ks_combsort_##name(size_t n, type_t a[])						\
	{																	\
		const double shrink_factor = 1.2473309501039786540366528676643; \
		bool do_swap;													\
		size_t gap = n;													\
		type_t *i, *j;												\
		do {															\
			if (gap > 2) {												\
				gap = (size_t)(gap / shrink_factor);					\
				if (gap == 9 || gap == 10) gap = 11;					\
			}															\
			do_swap = false;												\
			for (i = a; i < a + n - gap; ++i) {							\
				j = i + gap;											\
				if (__sort_lt(*j, *i)) {                       \
                    std::swap(*i, *j);                          \
					do_swap = true;										\
				}														\
			}															\
		} while (do_swap || gap > 2);									\
		if (gap != 1) __ks_insertsort_##name(a, a + n);					\
	}                                                 \
    void ks_introsort_##name(size_t n, type_t a[])						\
	{																	\
		int d;															\
		ks_isort_stack_t *top, *stack;									\
		type_t rp;											\
		type_t *s, *t, *i, *j, *k;										\
																		\
		if (n <= 1) return;												\
		else if (n == 2) {												\
			if (__sort_lt(a[1], a[0])) { std::swap(a[0], a[1]); }       \
			return;														\
		}																\
		for (d = 2; 1ul<<d < n; ++d);									\
		stack = new ks_isort_stack_t[(sizeof(size_t) * d) + 2];         \
		top = stack; s = a; t = a + (n-1); d <<= 1;						\
		while (1) {														\
			if (s < t) {												\
				if (--d == 0) {											\
					ks_combsort_##name(t - s + 1, s);					\
					t = s;												\
					continue;											\
				}														\
				i = s; j = t; k = i + ((j-i)>>1) + 1;					\
				if (__sort_lt(*k, *i)) {								\
					if (__sort_lt(*k, *j)) k = j;						\
				} else k = __sort_lt(*j, *i)? i : j;					\
				rp = *k;												\
				if (k != t) { std::swap(*k, *t); }	                    \
				for (;;) {												\
					do ++i; while (__sort_lt(*i, rp));					\
					do --j; while (i <= j && __sort_lt(rp, *j));		\
					if (j <= i) break;                                  \
					std::swap(*i, *j);                                  \
				}														\
				std::swap(*i, *t);					                    \
				if (i-s > t-i) {										\
					if (i-s > 16) { top->left = s; top->right = i-1; top->depth = d; ++top; } \
					s = t-i > 16? i+1 : t;								\
				} else {												\
					if (t-i > 16) { top->left = i+1; top->right = t; top->depth = d; ++top; } \
					t = i-s > 16? i-1 : s;								\
				}														\
			} else {													\
				if (top == stack) {										\
					delete[] stack;								        \
					__ks_insertsort_##name(a, a+n);						\
					return;												\
				} else { --top; s = (type_t*)top->left; t = (type_t*)top->right; d = top->depth; } \
			}															\
		}																\
	}

#define ks_introsort(name, n, a) ks_introsort_##name(n, a)
#define ks_lt_generic(a, b) ((a) < (b))

#endif
