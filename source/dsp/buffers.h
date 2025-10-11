#pragma once

#include <vector>
#include <cstddef>

namespace kbeyond::dsp {

class DelayLine {
public:
    DelayLine() = default;

    void setup(std::size_t len);
    double readFrac(double delaySamp) const;
    double readInt(long delaySamp) const;
    void write(double value);
    std::size_t size() const noexcept { return buffer_.size(); }

private:
    std::vector<double> buffer_{};
    long writeIndex_ = 0;
};

} // namespace kbeyond::dsp

