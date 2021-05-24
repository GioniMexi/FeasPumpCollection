/**
 * @file asserter.h
 * @brief assert utilities
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#ifndef ASSERTER_H
#define ASSERTER_H

#include <cassert>

#ifndef NDEBUG
#define DOMINIQS_ASSERT(condition) do { assert((condition)); } while(false)
#define DOMINIQS_ASSERT_EQUAL(first, second) do { assert((first)==(second)); } while(false)
#else
#define DOMINIQS_ASSERT(unused) do {} while(false)
#define DOMINIQS_ASSERT_EQUAL(unused1, unused2) do {} while(false)
#endif

#endif /* ASSERTER_H */
