#ifndef TIMER_H_
#define TIMER_H_
#include <iostream>
#include <chrono>

class Timer {
    std::chrono::time_point<std::chrono::high_resolution_clock> st, ed;
    std::chrono::duration<double> elapsed;
    double t;

public:
    Timer();
    void start();
    void stop_and_log(char const *label);
};

#endif