#include "buffers.h"

#include <algorithm>
#include <cmath>

namespace kbeyond::dsp {

void DelayLine::setup(std::size_t len) {
    buffer_.assign(len, 0.0);
    writeIndex_ = 0;
}

double DelayLine::readFrac(double delaySamp) const {
    const std::size_t len = buffer_.size();
    if (len == 0)
        return 0.0;
    double readPos = static_cast<double>(writeIndex_) - delaySamp;
    readPos -= std::floor(readPos / static_cast<double>(len)) * static_cast<double>(len);
    const long i0 = static_cast<long>(readPos);
    const long i1 = (i0 + 1) % static_cast<long>(len);
    const double frac = readPos - static_cast<double>(i0);
    const double s0 = buffer_[static_cast<std::size_t>(i0)];
    const double s1 = buffer_[static_cast<std::size_t>(i1)];
    return s0 + (s1 - s0) * frac;
}

double DelayLine::readInt(long delaySamp) const {
    const std::size_t len = buffer_.size();
    if (len == 0)
        return 0.0;
    long rp = writeIndex_ - delaySamp;
    while (rp < 0)
        rp += static_cast<long>(len);
    return buffer_[static_cast<std::size_t>(rp % static_cast<long>(len))];
}

void DelayLine::write(double value) {
    if (buffer_.empty())
        return;
    buffer_[static_cast<std::size_t>(writeIndex_)] = value;
    writeIndex_ = (writeIndex_ + 1) % static_cast<long>(buffer_.size());
}

} // namespace kbeyond::dsp

