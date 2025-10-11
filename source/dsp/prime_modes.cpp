#include "prime_modes.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace prime_modes {
namespace {

double fract(double x) noexcept {
    return x - std::floor(x);
}

unsigned next_prime(unsigned candidate) noexcept {
    if (candidate <= 2U)
        return 2U;
    if ((candidate & 1U) == 0U)
        ++candidate;
    while (true) {
        bool isPrime = true;
        for (unsigned d = 3U; d * d <= candidate; d += 2U) {
            if (candidate % d == 0U) {
                isPrime = false;
                break;
            }
        }
        if (isPrime)
            return candidate;
        candidate += 2U;
    }
}

void enforce_bounds(std::vector<double> &values) {
    if (values.empty())
        return;
    const double eps = 1.0 / (static_cast<double>(values.size()) * 128.0);
    double minBound = eps;
    double maxBound = 1.0 - eps * static_cast<double>(values.size());
    if (maxBound <= minBound)
        maxBound = 0.999999;
    values[0] = std::clamp(values[0], minBound, maxBound);
    for (std::size_t i = 1; i < values.size(); ++i) {
        double minVal = values[i - 1] + eps;
        double maxVal = 1.0 - eps * static_cast<double>(values.size() - i);
        if (maxVal <= minVal)
            maxVal = minVal + eps;
        values[i] = std::clamp(values[i], minVal, maxVal);
    }
}

std::uint32_t advance_seed(std::uint32_t state) noexcept {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state ? state : 0x9E3779B9u;
}

void apply_seed(std::vector<double> &values, std::uint32_t seed) {
    if (values.empty())
        return;
    std::uint32_t state = seed ? seed : 0xA511E9B5u;
    const double jitterBase = 1.0 / (static_cast<double>(values.size()) * 256.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        state = advance_seed(state + static_cast<std::uint32_t>(i * 97U));
        double r = static_cast<double>(state & 0xFFFFu) / 65535.0;
        double jitter = (r - 0.5) * jitterBase;
        values[i] += jitter;
    }
    std::sort(values.begin(), values.end());
    enforce_bounds(values);
}

std::vector<double> make_prime_sequence(std::size_t count) {
    std::vector<double> seq;
    seq.reserve(count);
    unsigned prime = 2U;
    for (std::size_t i = 0; i < count; ++i) {
        if (i == 0)
            prime = 2U;
        else
            prime = next_prime(prime + 1U);
        seq.push_back(static_cast<double>(prime));
    }
    double maxVal = seq.empty() ? 1.0 : seq.back();
    if (maxVal <= 0.0)
        maxVal = 1.0;
    for (double &v : seq)
        v = v / maxVal;
    return seq;
}

std::vector<double> make_aureo_sequence(std::size_t count) {
    std::vector<double> seq;
    seq.reserve(count);
    const double phi = 1.6180339887498948482;
    const double phi2 = phi * phi;
    for (std::size_t i = 0; i < count; ++i) {
        double n = static_cast<double>(i + 1);
        double a = fract(n * phi);
        double b = fract(n * phi2);
        seq.push_back(0.5 * (a + b));
    }
    std::sort(seq.begin(), seq.end());
    return seq;
}

std::vector<double> make_plastica_sequence(std::size_t count) {
    std::vector<double> seq;
    seq.reserve(count);
    const double rho = 1.3247179572447458;
    const double rho2 = rho * rho;
    const double rho3 = rho2 * rho;
    for (std::size_t i = 0; i < count; ++i) {
        double n = static_cast<double>(i + 1);
        double a = fract(n * rho);
        double b = fract(n * rho2);
        double c = fract(n * rho3);
        seq.push_back((a + b + c) / 3.0);
    }
    std::sort(seq.begin(), seq.end());
    return seq;
}

std::vector<double> make_prime_aureo_sequence(std::size_t count) {
    std::vector<double> seq;
    seq.reserve(count);
    const double phi = 1.6180339887498948482;
    unsigned prime = 2U;
    for (std::size_t i = 0; i < count; ++i) {
        if (i == 0)
            prime = 2U;
        else
            prime = next_prime(prime + 1U);
        double v = fract(static_cast<double>(prime) * phi);
        seq.push_back(v);
    }
    std::sort(seq.begin(), seq.end());
    return seq;
}

} // namespace

Pattern pattern_from_string(const char *name, Pattern fallback) noexcept {
    if (!name)
        return fallback;
    if (std::strcmp(name, "prime") == 0)
        return Pattern::Prime;
    if (std::strcmp(name, "aureo") == 0)
        return Pattern::Aureo;
    if (std::strcmp(name, "plastica") == 0)
        return Pattern::Plastica;
    if (std::strcmp(name, "prime_aureo") == 0)
        return Pattern::PrimeAureo;
    return fallback;
}

const char *to_string(Pattern pattern) noexcept {
    switch (pattern) {
    case Pattern::Prime:
        return "prime";
    case Pattern::Aureo:
        return "aureo";
    case Pattern::Plastica:
        return "plastica";
    case Pattern::PrimeAureo:
        return "prime_aureo";
    default:
        return "aureo";
    }
}

std::vector<double> generate_pattern(Pattern pattern, std::size_t count, std::uint32_t seed) {
    std::vector<double> values;
    if (count == 0)
        return values;
    switch (pattern) {
    case Pattern::Prime:
        values = make_prime_sequence(count);
        break;
    case Pattern::Aureo:
        values = make_aureo_sequence(count);
        break;
    case Pattern::Plastica:
        values = make_plastica_sequence(count);
        break;
    case Pattern::PrimeAureo:
        values = make_prime_aureo_sequence(count);
        break;
    default:
        values = make_aureo_sequence(count);
        break;
    }
    apply_seed(values, seed);
    return values;
}

} // namespace prime_modes

