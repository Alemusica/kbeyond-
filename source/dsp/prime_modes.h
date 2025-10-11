#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace prime_modes {

enum class Pattern {
    Prime,
    Aureo,
    Plastica,
    PrimeAureo
};

Pattern pattern_from_string(const char *name, Pattern fallback = Pattern::Aureo) noexcept;

const char *to_string(Pattern pattern) noexcept;

std::vector<double> generate_pattern(Pattern pattern, std::size_t count, std::uint32_t seed);

} // namespace prime_modes

