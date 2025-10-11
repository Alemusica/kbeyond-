#pragma once

#include <cstdlib>
#include <iostream>

namespace kbeyond::dsp::detail {
inline void dsp_assert_fail(const char *expr,
                            const char *file,
                            int line,
                            const char *message) {
    std::cerr << "[dsp_assert] Assertion failed: " << expr << " at " << file << ':' << line;
    if (message && *message)
        std::cerr << " -> " << message;
    std::cerr << std::endl;
    std::abort();
}
} // namespace kbeyond::dsp::detail

#ifndef NDEBUG
#define dsp_assert(cond)                                                                           \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            kbeyond::dsp::detail::dsp_assert_fail(#cond, __FILE__, __LINE__, nullptr);             \
        }                                                                                          \
    } while (0)
#define dsp_assert_msg(cond, msg)                                                                  \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            kbeyond::dsp::detail::dsp_assert_fail(#cond, __FILE__, __LINE__, (msg));               \
        }                                                                                          \
    } while (0)
#else
#define dsp_assert(cond) ((void)0)
#define dsp_assert_msg(cond, msg) ((void)0)
#endif

