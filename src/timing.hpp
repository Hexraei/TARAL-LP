#pragma once
#include <chrono>

// Wall-clock stopwatch on the steady clock. Construction starts it; calling it returns the seconds elapsed.
struct Stopwatch {
    using Clock = std::chrono::steady_clock;
    Clock::time_point start = Clock::now();
    double operator()() const { return std::chrono::duration<double>(Clock::now() - start).count(); }
};
